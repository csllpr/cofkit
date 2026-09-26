from __future__ import annotations

from dataclasses import dataclass, field
from typing import Mapping

from .geometry import (
    Mat3,
    Vec3,
    add,
    cross,
    dot,
    matmul_vec,
    norm,
    orthogonal_component,
    rotation_from_frame_to_axes,
    safe_normalize,
    scale,
    sub,
)
from .linkage_geometry import effective_motif_origin
from .model import AssemblyState, MonomerSpec, Pose, ReactionTemplate
from .scoring import BridgeGeometryReport, CandidateScorer
from .search import AssignmentOutcome


@dataclass(frozen=True)
class OptimizerConfig:
    """Iteration caps and step sizes for the refinement pass.

    There is no c-axis move by design: optimization runs on the pre-stacking
    embedding state (vacuum-slab/implicit-repeat cell) — stacking expansion
    happens downstream — so only the lateral cell axes are meaningful here.
    """

    max_iterations: int = 8
    translation_step: float = 0.35
    cell_scale_step: float = 0.5
    min_lateral_scale: float = 0.75
    max_lateral_scale: float = 1.25


@dataclass(frozen=True)
class OptimizationResult:
    state: AssemblyState
    metrics: Mapping[str, object] = field(default_factory=dict)


class ContinuousOptimizer:
    """Applies a small dependency-free bridge-geometry refinement pass.

    Runs on the pre-stacking embedding state (vacuum-slab/implicit-repeat
    cell); stacking expansion happens downstream of every call site
    (``engine.py`` / ``batch.py`` feed it the embedder's state directly), so
    the lateral cell move intentionally leaves the c axis untouched.
    """

    def __init__(
        self,
        config: OptimizerConfig | None = None,
        scorer: CandidateScorer | None = None,
    ):
        self.config = config or OptimizerConfig()
        self.scorer = scorer or CandidateScorer()

    def optimize(
        self,
        outcome: AssignmentOutcome,
        state: AssemblyState,
        monomer_specs: Mapping[str, MonomerSpec],
        templates: Mapping[str, ReactionTemplate],
    ) -> OptimizationResult:
        initial_report = self.scorer.bridge_geometry_report(outcome, state, monomer_specs, templates)
        best_state = state
        best_report = initial_report
        accepted_iterations = 0

        if not initial_report.event_metrics:
            return OptimizationResult(
                state=state,
                metrics={
                    "enabled": False,
                    "iterations": 0,
                    "accepted_iterations": 0,
                    "initial_residual": initial_report.total_residual,
                    "final_residual": initial_report.total_residual,
                    "improved": False,
                },
            )

        working_state = state
        for _ in range(self.config.max_iterations):
            proposal = self._propose_state(working_state, outcome, monomer_specs, templates)
            proposal_report = self.scorer.bridge_geometry_report(outcome, proposal, monomer_specs, templates)
            if proposal_report.total_residual + 1e-9 < best_report.total_residual:
                best_state = proposal
                best_report = proposal_report
                working_state = proposal
                accepted_iterations += 1
            else:
                working_state = best_state

        return OptimizationResult(
            state=best_state,
            metrics={
                "enabled": True,
                "iterations": self.config.max_iterations,
                "accepted_iterations": accepted_iterations,
                "initial_residual": initial_report.total_residual,
                "final_residual": best_report.total_residual,
                "improved": best_report.total_residual + 1e-9 < initial_report.total_residual,
            },
        )

    def _propose_state(
        self,
        state: AssemblyState,
        outcome: AssignmentOutcome,
        monomer_specs: Mapping[str, MonomerSpec],
        templates: Mapping[str, ReactionTemplate],
    ) -> AssemblyState:
        working_state = state
        working_report = self.scorer.bridge_geometry_report(outcome, working_state, monomer_specs, templates)
        refiners = (
            self._scale_lateral_cell,
            self._refine_translations,
            self._refine_orientations,
        )
        for refiner in refiners:
            proposal = refiner(working_state, outcome, monomer_specs, templates)
            proposal_report = self.scorer.bridge_geometry_report(outcome, proposal, monomer_specs, templates)
            if proposal_report.total_residual + 1e-9 < working_report.total_residual:
                working_state = proposal
                working_report = proposal_report
        return working_state

    def _scale_lateral_cell(
        self,
        state: AssemblyState,
        outcome: AssignmentOutcome,
        monomer_specs: Mapping[str, MonomerSpec],
        templates: Mapping[str, ReactionTemplate],
    ) -> AssemblyState:
        distance_ratios: list[float] = []
        for event in outcome.events:
            template = templates[event.template_id]
            if template.topology_role != "bridge" or len(event.participants) != 2:
                continue
            first, second = event.participants
            pose1 = state.monomer_poses[first.monomer_instance_id]
            pose2 = state.monomer_poses[second.monomer_instance_id]
            motif1 = monomer_specs[first.monomer_id].motif_by_id(first.motif_id)
            motif2 = monomer_specs[second.monomer_id].motif_by_id(second.motif_id)
            origin1 = self._world_motif_origin(
                state.cell,
                pose1,
                effective_motif_origin(event.template_id, monomer_specs[first.monomer_id], motif1),
                first.periodic_image,
            )
            origin2 = self._world_motif_origin(
                state.cell,
                pose2,
                effective_motif_origin(event.template_id, monomer_specs[second.monomer_id], motif2),
                second.periodic_image,
            )
            separation = norm(sub(origin2, origin1))
            if separation < 1e-8:
                continue
            distance_ratios.append(self.scorer._target_distance(template) / separation)

        if not distance_ratios:
            return state

        avg_ratio = sum(distance_ratios) / len(distance_ratios)
        step_ratio = 1.0 + (avg_ratio - 1.0) * self.config.cell_scale_step
        step_ratio = max(self.config.min_lateral_scale, min(self.config.max_lateral_scale, step_ratio))

        cell = (
            scale(state.cell[0], step_ratio),
            scale(state.cell[1], step_ratio),
            state.cell[2],
        )
        monomer_poses = {
            instance_id: Pose(
                translation=(
                    pose.translation[0] * step_ratio,
                    pose.translation[1] * step_ratio,
                    pose.translation[2],
                ),
                rotation_matrix=pose.rotation_matrix,
            )
            for instance_id, pose in state.monomer_poses.items()
        }
        return AssemblyState(
            cell=cell,
            monomer_poses=monomer_poses,
            torsions=state.torsions,
            layer_offsets=state.layer_offsets,
            stacking_state=state.stacking_state,
        )

    def _refine_translations(
        self,
        state: AssemblyState,
        outcome: AssignmentOutcome,
        monomer_specs: Mapping[str, MonomerSpec],
        templates: Mapping[str, ReactionTemplate],
    ) -> AssemblyState:
        translation_updates = {
            instance_id: (0.0, 0.0, 0.0) for instance_id in state.monomer_poses
        }
        counts = {instance_id: 0 for instance_id in state.monomer_poses}

        for event in outcome.events:
            template = templates[event.template_id]
            if template.topology_role != "bridge" or len(event.participants) != 2:
                continue

            first, second = event.participants
            pose1 = state.monomer_poses[first.monomer_instance_id]
            pose2 = state.monomer_poses[second.monomer_instance_id]
            motif1 = monomer_specs[first.monomer_id].motif_by_id(first.motif_id)
            motif2 = monomer_specs[second.monomer_id].motif_by_id(second.motif_id)
            origin1 = self._world_motif_origin(
                state.cell,
                pose1,
                effective_motif_origin(event.template_id, monomer_specs[first.monomer_id], motif1),
                first.periodic_image,
            )
            origin2 = self._world_motif_origin(
                state.cell,
                pose2,
                effective_motif_origin(event.template_id, monomer_specs[second.monomer_id], motif2),
                second.periodic_image,
            )

            delta = sub(origin2, origin1)
            direction = self._safe_normalize(delta)
            target = self.scorer._target_distance(template)
            distance_error = norm(delta) - target

            normal1 = self._safe_normalize(matmul_vec(pose1.rotation_matrix, motif1.frame.normal))
            normal2 = self._safe_normalize(matmul_vec(pose2.rotation_matrix, motif2.frame.normal))
            plane_normal = self._safe_normalize(add(normal1, normal2))
            if norm(add(normal1, normal2)) < 1e-8:
                plane_normal = normal1
            planarity_error = dot(delta, plane_normal)

            correction = add(
                scale(direction, 0.5 * distance_error * self.config.translation_step),
                scale(plane_normal, 0.35 * planarity_error * self.config.translation_step),
            )
            translation_updates[first.monomer_instance_id] = add(
                translation_updates[first.monomer_instance_id],
                correction,
            )
            translation_updates[second.monomer_instance_id] = add(
                translation_updates[second.monomer_instance_id],
                scale(correction, -1.0),
            )
            counts[first.monomer_instance_id] += 1
            counts[second.monomer_instance_id] += 1

        monomer_poses: dict[str, Pose] = {}
        for instance_id, pose in state.monomer_poses.items():
            update = translation_updates[instance_id]
            count = max(1, counts[instance_id])
            monomer_poses[instance_id] = Pose(
                translation=add(pose.translation, scale(update, 1.0 / count)),
                rotation_matrix=pose.rotation_matrix,
            )

        return AssemblyState(
            cell=state.cell,
            monomer_poses=monomer_poses,
            torsions=state.torsions,
            layer_offsets=state.layer_offsets,
            stacking_state=state.stacking_state,
        )

    def _refine_orientations(
        self,
        state: AssemblyState,
        outcome: AssignmentOutcome,
        monomer_specs: Mapping[str, MonomerSpec],
        templates: Mapping[str, ReactionTemplate],
    ) -> AssemblyState:
        event_lookup = {metrics.event_id: metrics for metrics in self.scorer.bridge_geometry_report(outcome, state, monomer_specs, templates).event_metrics}
        connected_events: dict[str, list[tuple[float, object, object, ReactionTemplate]]] = {}
        for event in outcome.events:
            template = templates[event.template_id]
            if template.topology_role != "bridge" or len(event.participants) != 2:
                continue
            metrics = event_lookup.get(event.id)
            if metrics is None:
                continue
            first, second = event.participants
            connected_events.setdefault(first.monomer_instance_id, []).append((metrics.total_residual, event, first, template))
            connected_events.setdefault(second.monomer_instance_id, []).append((metrics.total_residual, event, second, template))

        monomer_poses = dict(state.monomer_poses)
        for instance_id, entries in connected_events.items():
            pose = monomer_poses[instance_id]
            proposals: list[tuple[float, Mat3]] = []
            pivot_residual = -1.0
            pivot_attachment: Vec3 | None = None
            pivot_local_origin: Vec3 | None = None
            pivot_image: tuple[int, int, int] = (0, 0, 0)
            for residual, event, participant, _template in entries:
                other = event.participants[0] if event.participants[1] == participant else event.participants[1]
                other_pose = monomer_poses[other.monomer_instance_id]
                motif = monomer_specs[participant.monomer_id].motif_by_id(participant.motif_id)
                other_motif = monomer_specs[other.monomer_id].motif_by_id(other.motif_id)
                local_origin = effective_motif_origin(event.template_id, monomer_specs[participant.monomer_id], motif)
                origin = self._world_motif_origin(
                    state.cell,
                    pose,
                    local_origin,
                    participant.periodic_image,
                )
                other_origin = self._world_motif_origin(
                    state.cell,
                    other_pose,
                    effective_motif_origin(event.template_id, monomer_specs[other.monomer_id], other_motif),
                    other.periodic_image,
                )
                target_primary = self._safe_normalize(sub(other_origin, origin))

                other_normal = self._safe_normalize(matmul_vec(other_pose.rotation_matrix, other_motif.frame.normal))
                target_normal = self._orthogonal_component(other_normal, target_primary)
                if norm(target_normal) < 1e-8:
                    current_normal = matmul_vec(pose.rotation_matrix, motif.frame.normal)
                    target_normal = self._orthogonal_component(current_normal, target_primary)
                if norm(target_normal) < 1e-8:
                    target_normal = (0.0, 0.0, 1.0)

                rotation = rotation_from_frame_to_axes(
                    motif.frame,
                    target_primary,
                    self._safe_normalize(target_normal),
                )
                proposals.append((max(residual, 0.0), rotation))
                if residual > pivot_residual:
                    pivot_residual = residual
                    pivot_attachment = origin
                    pivot_local_origin = local_origin
                    pivot_image = participant.periodic_image

            if pivot_attachment is None or pivot_local_origin is None:
                continue
            # Blend every incident event's desired rotation, weighted by the
            # event's residual, instead of steering by the single worst event.
            rotation = self._blend_rotations(proposals)
            if rotation is None:
                rotation = max(proposals, key=lambda item: item[0])[1]
            # Rotate about the pivot motif's bridge attachment point (the
            # worst-residual event's) so the rotation does not displace the
            # endpoint it is aligning: translation' keeps
            # translation' + R'·local_origin + image_offset == attachment.
            # Other motifs' attachment points move; the greedy acceptance test
            # arbitrates.
            image_offset = add(
                add(scale(state.cell[0], pivot_image[0]), scale(state.cell[1], pivot_image[1])),
                scale(state.cell[2], pivot_image[2]),
            )
            new_translation = sub(
                sub(pivot_attachment, image_offset),
                matmul_vec(rotation, pivot_local_origin),
            )
            monomer_poses[instance_id] = Pose(
                translation=new_translation,
                rotation_matrix=rotation,
            )

        return AssemblyState(
            cell=state.cell,
            monomer_poses=monomer_poses,
            torsions=state.torsions,
            layer_offsets=state.layer_offsets,
            stacking_state=state.stacking_state,
        )

    def _blend_rotations(self, proposals: list[tuple[float, Mat3]]) -> Mat3 | None:
        """Weighted average of desired rotations, re-orthonormalized.

        Weights are the per-event residuals. The element-wise mean of rotation
        matrices is not itself a rotation, so the averaged rows are
        orthonormalized by Gram–Schmidt (row three completed with the cross
        product). Returns ``None`` when the mean is degenerate — strongly
        opposing proposals cancel — in which case the caller falls back to the
        worst-event proposal.
        """
        total_weight = sum(weight for weight, _ in proposals)
        if total_weight < 1e-12:
            return None
        averaged = tuple(
            tuple(
                sum(weight * rotation[i][j] for weight, rotation in proposals) / total_weight
                for j in range(3)
            )
            for i in range(3)
        )
        row1 = averaged[0]
        if norm(row1) < 1e-8:
            return None
        u1 = scale(row1, 1.0 / norm(row1))
        row2_orthogonal = orthogonal_component(averaged[1], u1, axis_is_unit=True)
        if norm(row2_orthogonal) < 1e-8:
            return None
        u2 = scale(row2_orthogonal, 1.0 / norm(row2_orthogonal))
        return (u1, u2, cross(u1, u2))

    def _world_motif_origin(
        self,
        cell: tuple[Vec3, Vec3, Vec3],
        pose: Pose,
        local_origin: Vec3,
        periodic_image: tuple[int, int, int],
    ) -> Vec3:
        rotated = matmul_vec(pose.rotation_matrix, local_origin)
        return add(
            add(pose.translation, rotated),
            add(
                add(scale(cell[0], periodic_image[0]), scale(cell[1], periodic_image[1])),
                scale(cell[2], periodic_image[2]),
            ),
        )

    def _orthogonal_component(self, vector: Vec3, axis: Vec3) -> Vec3:
        # Call sites pass pre-normalized axes (unit vectors built above).
        return orthogonal_component(vector, axis, axis_is_unit=True)

    def _safe_normalize(self, vector: Vec3) -> Vec3:
        # Fallback (1, 0, 0) is a fabricated +x direction preserved from the
        # pre-consolidation implementation; it feeds rotation targets, so a
        # triggered fallback indicates a data problem upstream (audit A10).
        return safe_normalize(vector, fallback=(1.0, 0.0, 0.0))
