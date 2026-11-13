from __future__ import annotations

from dataclasses import dataclass, field, replace
from math import isfinite, sqrt
from sys import float_info
from typing import Mapping

from .geometry import (
    Mat3,
    Vec3,
    add,
    cross,
    covariance_eigenpairs,
    dot,
    frame_axes,
    matmul_vec,
    norm,
    normalize,
    orthogonal_component,
    rotation_from_frame_to_axes,
    scale,
    sub,
    transpose,
)
from .linkage_geometry import effective_motif_origin
from .model import AssemblyState, MonomerSpec, Pose, ReactionTemplate
from .scoring import BridgeGeometryReport, CandidateScorer
from .search import AssignmentOutcome

# Heuristic — pending calibration: weight of the out-of-plane (planarity)
# term relative to the along-bridge distance term in the translation
# correction of ContinuousOptimizer._refine_translations. Its value
# coincides numerically with OptimizerConfig.translation_step's default but
# is a separate quantity; the two must not be merged.
_PLANARITY_CORRECTION_WEIGHT = 0.35

# Derived: normal equations lose conditioning quadratically; sqrt(machine
# epsilon) is a relative pivot guard, independent of coordinate units.
_CELL_FIT_RELATIVE_TOLERANCE = sqrt(float_info.epsilon)
# Derived: normalizing the difference of nearly equal unit directions loses
# angular precision at sqrt(machine epsilon).
_CELL_FIT_DIRECTION_TOLERANCE = sqrt(float_info.epsilon)
# Derived: binary backtracking uses step_next = step / 2.
_CELL_FIT_BACKTRACK_FACTOR = 0.5
# Heuristic — pending calibration: numerical tie guard for the surrogate
# residual, far below geometric screening tolerances.
_OPTIMIZER_RESIDUAL_IMPROVEMENT_EPSILON = 1e-9
# Heuristic — pending calibration: numerical cancellation guards for the
# residual-weighted rotation blend (not geometric acceptance thresholds).
_ROTATION_BLEND_MIN_AXIS_NORM = 1e-8
_ROTATION_BLEND_MIN_TOTAL_WEIGHT = 1e-12


@dataclass(frozen=True)
class OptimizerConfig:
    """Iteration caps and step sizes for the refinement pass.

    A 2D seed varies a/b within their plane and preserves c (slab padding).
    A 3D seed varies all three lattice vectors. These are geometric closure
    proposals, not physical cell relaxation.
    """

    # Heuristic — pending calibration: small greedy budget; every iteration
    # re-scores all bridge events, so the cap keeps the pass cheap.
    max_iterations: int = 8
    # Heuristic — pending calibration: damping fraction of each per-event
    # error applied per step (sub-unit for stability of the greedy loop).
    translation_step: float = 0.35
    # Heuristic — pending calibration: damping fraction of the mean
    # target/actual distance ratio applied per uniform active-cell step.
    cell_scale_step: float = 0.5
    # Heuristic — pending calibration: uniform scale / principal-stretch
    # bounds per cell proposal. Legacy field names remain compatible.
    min_lateral_scale: float = 0.75
    max_lateral_scale: float = 1.25
    cell_shape_refinement: bool = True
    coupled_pose_cell_refinement: bool = True
    # Derived: a unit step applies the fitted lattice itself. Backtracking
    # enforces trust bounds and closure preservation, rather than assuming a
    # damped chord between bridge vectors preserves their lengths.
    cell_fit_step: float = 1.0
    # Heuristic — pending calibration: bounded halvings of a proposed affine
    # cell step to keep its principal stretches within the existing scale
    # limits. Rank-deficient fits retain the scalar/pose refinement path.
    cell_fit_backtracking_steps: int = 8


@dataclass(frozen=True)
class OptimizationResult:
    state: AssemblyState
    metrics: Mapping[str, object] = field(default_factory=dict)


class ContinuousOptimizer:
    """Applies a small dependency-free bridge-geometry refinement pass.

    Stacking expansion happens downstream. Only 2D cells have a protected
    padding axis; pre-stacking 3D nets still need all lattice directions to
    satisfy periodic linkage constraints.
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
            if proposal_report.total_residual + _OPTIMIZER_RESIDUAL_IMPROVEMENT_EPSILON < best_report.total_residual:
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
                "improved": best_report.total_residual + _OPTIMIZER_RESIDUAL_IMPROVEMENT_EPSILON < initial_report.total_residual,
                "cell_shape_refinement": self.config.cell_shape_refinement,
                "cell_fit_step": self.config.cell_fit_step,
                "coupled_pose_cell_refinement": self.config.coupled_pose_cell_refinement,
                "active_cell_axes": ("a", "b", "c")[:self._active_cell_dimensions(outcome)],
                "initial_cell": state.cell,
                "final_cell": best_state.cell,
                "physical_relaxation": False,
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
            (self._fit_cell_shape, True),
            (self._scale_lateral_cell, False),
            (self._refine_translations, False),
            (self._refine_orientations, False),
        )
        for refiner, protect_bonds in refiners:
            proposal = refiner(working_state, outcome, monomer_specs, templates)
            proposal_report = self.scorer.bridge_geometry_report(outcome, proposal, monomer_specs, templates)
            if self._improves_report(working_report, proposal_report, protect_bonds=protect_bonds):
                working_state = proposal
                working_report = proposal_report
        if self.config.coupled_pose_cell_refinement:
            # A rigid rotation can be locally worse until its changed site
            # offsets are accommodated by cell/center placement. Evaluate
            # that coupled proposal as a whole, using the unchanged scorer;
            # no intermediate state is promoted on its own.
            coupled = self._refine_orientations(state, outcome, monomer_specs, templates)
            coupled = self._fit_cell_shape(coupled, outcome, monomer_specs, templates)
            coupled = self._refine_translations(coupled, outcome, monomer_specs, templates)
            coupled_report = self.scorer.bridge_geometry_report(outcome, coupled, monomer_specs, templates)
            if self._improves_report(working_report, coupled_report, protect_bonds=True):
                working_state = coupled
        return working_state

    @staticmethod
    def _improves_report(before, after, *, protect_bonds=False):
        if after.total_residual + _OPTIMIZER_RESIDUAL_IMPROVEMENT_EPSILON >= before.total_residual:
            return False
        if protect_bonds:
            # Cell/compound proposals cannot purchase a better angular
            # surrogate by opening an already-closed periodic linkage.
            previous = {m.event_id: m.distance_residual for m in before.event_metrics}
            if any(m.distance_residual > previous[m.event_id] + _OPTIMIZER_RESIDUAL_IMPROVEMENT_EPSILON for m in after.event_metrics):
                return False
        return True

    @staticmethod
    def _active_cell_dimensions(outcome: AssignmentOutcome) -> int:
        topology = outcome.assignment_plan.net_plan.topology
        # Unspecified topology retains the established slab convention.
        return 3 if topology is not None and topology.dimensionality == "3D" else 2

    @staticmethod
    def _fractional(point: Vec3, cell: tuple[Vec3, Vec3, Vec3]) -> Vec3:
        a, b, c = cell
        volume = dot(a, cross(b, c))
        return (
            dot(point, cross(b, c)) / volume,
            dot(point, cross(c, a)) / volume,
            dot(point, cross(a, b)) / volume,
        )

    def _with_cell(self, state: AssemblyState, cell: tuple[Vec3, Vec3, Vec3]) -> AssemblyState:
        # Move molecular centers affinely, preserving their fractional
        # positions and every rigid monomer's internal Cartesian geometry.
        poses = {}
        for iid, pose in state.monomer_poses.items():
            fractional = self._fractional(pose.translation, state.cell)
            center = tuple(sum(f * vector[i] for f, vector in zip(fractional, cell)) for i in range(3))
            poses[iid] = replace(pose, translation=center)
        return replace(state, cell=cell, monomer_poses=poses)

    @staticmethod
    def _solve_cell_fit(matrix, rhs):
        size = len(matrix)
        scale_max = max((abs(x) for row in matrix for x in row), default=0.0)
        if not isfinite(scale_max) or scale_max == 0:
            return None
        augmented = [list(row) + list(target) for row, target in zip(matrix, rhs)]
        for column in range(size):
            pivot = max(range(column, size), key=lambda row: abs(augmented[row][column]))
            if abs(augmented[pivot][column]) <= _CELL_FIT_RELATIVE_TOLERANCE * scale_max:
                return None
            augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
            divisor = augmented[column][column]
            augmented[column] = [x / divisor for x in augmented[column]]
            for row in range(size):
                if row == column:
                    continue
                factor = augmented[row][column]
                augmented[row] = [x - factor * y for x, y in zip(augmented[row], augmented[column])]
        solution = tuple(tuple(row[size:]) for row in augmented)
        return solution if all(isfinite(x) for row in solution for x in row) else None

    def _cell_step_is_bounded(self, original, proposed):
        if dot(original[0], cross(original[1], original[2])) * dot(proposed[0], cross(proposed[1], proposed[2])) <= 0:
            return False
        columns = []
        for axis in range(3):
            unit = tuple(float(i == axis) for i in range(3))
            fractional = self._fractional(unit, original)
            columns.append(tuple(sum(f * vector[i] for f, vector in zip(fractional, proposed)) for i in range(3)))
        # Derived: the covariance of +/-F e_i is F F^T / 3. Its eigenvalues
        # therefore give the squared principal cell stretches divided by 3.
        pairs = covariance_eigenpairs(tuple(v for column in columns for v in (column, scale(column, -1.0))))
        if pairs is None:
            return False
        stretches = tuple(sqrt(max(0.0, 3.0 * value)) for value in pairs.eigenvalues)
        # The legacy field names remain compatible; the limits now bound
        # active-cell affine proposals, including the third axis of 3D nets.
        return min(stretches) >= self.config.min_lateral_scale and max(stretches) <= self.config.max_lateral_scale

    def _fit_cell_shape(self, state, outcome, monomer_specs, templates):
        """Fit independent cell vectors to the current linkage directions.

        With centers fixed in fractional coordinates, each bridge obeys
        sum_k q_k a_k + (R_j p_j - R_i p_i) = target * direction,
        where q = f_j + image_j - f_i - image_i. Solve the small linear
        least-squares system for a/b (2D) or a/b/c (3D), then damp/bound the
        proposal. The unchanged full scorer decides whether to accept it.
        """
        if not self.config.cell_shape_refinement:
            return state
        a, b, c = state.cell
        volume = dot(a, cross(b, c))
        if not isfinite(volume) or abs(volume) <= _CELL_FIT_RELATIVE_TOLERANCE * norm(a) * norm(b) * norm(c):
            return state
        active = self._active_cell_dimensions(outcome)
        matrix = [[0.0] * active for _ in range(active)]
        rhs = [[0.0] * 3 for _ in range(active)]
        fractions = {iid: self._fractional(p.translation, state.cell) for iid, p in state.monomer_poses.items()}
        for event in outcome.events:
            template = templates[event.template_id]
            if template.topology_role != "bridge" or len(event.participants) != 2:
                continue
            first, second = event.participants
            p1, p2 = (state.monomer_poses[p.monomer_instance_id] for p in (first, second))
            local = []
            primaries = []
            for ref, pose in ((first, p1), (second, p2)):
                spec = monomer_specs[ref.monomer_id]
                motif = spec.motif_by_id(ref.motif_id)
                local.append(matmul_vec(pose.rotation_matrix, effective_motif_origin(event.template_id, spec, motif)))
                primaries.append(self._safe_normalize(matmul_vec(pose.rotation_matrix, motif.frame.primary)))
            coefficients = tuple(f2 + i2 - f1 - i1 for f1, i1, f2, i2 in zip(fractions[first.monomer_instance_id], first.periodic_image, fractions[second.monomer_instance_id], second.periodic_image))
            offset = sub(local[1], local[0])
            delta = add(offset, tuple(sum(q * v[i] for q, v in zip(coefficients, state.cell)) for i in range(3)))
            length = norm(delta)
            if length == 0:
                continue
            # Derived: for unit bridge direction v the existing alignment
            # residual is 2 - v.(u_first - u_second). Its minimizing direction
            # is the normalized difference of the two outward primaries.
            # Coincident priors carry no directional information.
            direction = sub(primaries[0], primaries[1])
            direction = scale(direction, 1.0 / norm(direction)) if norm(direction) > _CELL_FIT_DIRECTION_TOLERANCE else scale(delta, 1.0 / length)
            target = sub(scale(direction, self.scorer._target_distance(template)), offset)
            for k in range(active, 3):
                target = sub(target, scale(state.cell[k], coefficients[k]))
            for row in range(active):
                for column in range(active):
                    matrix[row][column] += coefficients[row] * coefficients[column]
                for axis in range(3):
                    rhs[row][axis] += coefficients[row] * target[axis]
        fitted = self._solve_cell_fit(matrix, rhs)
        if fitted is None:
            return state
        if active == 2:
            normal = normalize(cross(a, b))
            fitted = tuple(sub(v, scale(normal, dot(v, normal))) for v in fitted)
        step = self.config.cell_fit_step
        before = self.scorer.bridge_geometry_report(outcome, state, monomer_specs, templates)
        for _ in range(max(0, self.config.cell_fit_backtracking_steps)):
            proposed = tuple(add(v, scale(sub(fitted[i], v), step)) if i < active else v for i, v in enumerate(state.cell))
            if self._cell_step_is_bounded(state.cell, proposed):
                trial = self._with_cell(state, proposed)
                report = self.scorer.bridge_geometry_report(outcome, trial, monomer_specs, templates)
                if self._improves_report(before, report, protect_bonds=True):
                    return trial
            step *= _CELL_FIT_BACKTRACK_FACTOR
        return state

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

        active = self._active_cell_dimensions(outcome)
        cell = tuple(scale(v, step_ratio) if i < active else v for i, v in enumerate(state.cell))
        return self._with_cell(state, cell)

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

            normal1 = matmul_vec(pose1.rotation_matrix, motif1.plane_normal)
            normal2 = matmul_vec(pose2.rotation_matrix, motif2.plane_normal)
            # Zero physical normals are the honest "no plane prior" marker
            # (non-planar/collinear/degenerate monomer conformer): steer the
            # planarity correction along whichever plane priors exist, and
            # skip it entirely when neither motif carries one rather than
            # fabricating a plane (impact-review claim T1-8).
            available_normals = tuple(normal for normal in (normal1, normal2) if norm(normal) >= 1e-8)
            planarity_correction = (0.0, 0.0, 0.0)
            if available_normals:
                unit_normals = tuple(normalize(normal) for normal in available_normals)
                normal_sum = unit_normals[0]
                for extra in unit_normals[1:]:
                    normal_sum = add(normal_sum, extra)
                if norm(normal_sum) < 1e-8:
                    # Antiparallel motif normals (a flat bridge) carry no plane
                    # direction; steer the planarity correction along normal1.
                    plane_normal = unit_normals[0]
                else:
                    plane_normal = self._safe_normalize(normal_sum)
                planarity_error = dot(delta, plane_normal)
                planarity_correction = scale(
                    plane_normal,
                    _PLANARITY_CORRECTION_WEIGHT * planarity_error * self.config.translation_step,
                )

            correction = add(
                scale(direction, 0.5 * distance_error * self.config.translation_step),
                planarity_correction,
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

                other_normal = matmul_vec(other_pose.rotation_matrix, other_motif.frame.normal)
                target_normal = (0.0, 0.0, 0.0)
                if norm(other_normal) >= 1e-8:
                    target_normal = self._orthogonal_component(self._safe_normalize(other_normal), target_primary)
                if norm(target_normal) < 1e-8 and norm(motif.frame.normal) >= 1e-8:
                    current_normal = matmul_vec(pose.rotation_matrix, motif.frame.normal)
                    target_normal = self._orthogonal_component(current_normal, target_primary)
                if norm(target_normal) < 1e-8:
                    # Neither motif supplies a numerical orientation normal:
                    # align the primary axis and preserve the current roll.
                    # frame_axes completes a numerical frame; scoring and
                    # translation terms use the separate physical priors.
                    current_roll_axis = matmul_vec(pose.rotation_matrix, frame_axes(motif.frame)[2])
                    target_normal = self._orthogonal_component(current_roll_axis, target_primary)
                if norm(target_normal) < 1e-8:
                    target_normal = (0.0, 0.0, 1.0) if abs(target_primary[2]) < 0.9 else (1.0, 0.0, 0.0)

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
        matrices is not itself a rotation, so the averaged columns (body axes
        in world coordinates) are orthonormalized by Gram–Schmidt. This is
        covariant under a rigid rotation of the whole periodic assembly;
        orthonormalizing rows would depend on the world's coordinate axes.
        Returns ``None`` when the mean is degenerate — strongly
        opposing proposals cancel — in which case the caller falls back to the
        worst-event proposal.
        """
        total_weight = sum(weight for weight, _ in proposals)
        if total_weight < _ROTATION_BLEND_MIN_TOTAL_WEIGHT:
            return None
        averaged = tuple(
            tuple(
                sum(weight * rotation[i][j] for weight, rotation in proposals) / total_weight
                for j in range(3)
            )
            for i in range(3)
        )
        columns = transpose(averaged)
        first = columns[0]
        if norm(first) < _ROTATION_BLEND_MIN_AXIS_NORM:
            return None
        u1 = scale(first, 1.0 / norm(first))
        second = orthogonal_component(columns[1], u1, axis_is_unit=True)
        if norm(second) < _ROTATION_BLEND_MIN_AXIS_NORM:
            return None
        u2 = scale(second, 1.0 / norm(second))
        return transpose((u1, u2, cross(u1, u2)))

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
        # Degenerate inputs here mean collapsed geometry (coincident motif
        # origins, a zero-length frame axis). Probe evidence (audit A10
        # follow-up): the pre-consolidation silent (1,0,0) fallback never
        # fired on any live path, so a degenerate vector now raises instead
        # of fabricating a direction; batch loops absorb the failure per
        # record via the standard error convention.
        if norm(vector) < 1e-8:
            raise ValueError(
                f"{type(self).__name__}: degenerate (near-zero-length) vector has no "
                f"direction to normalize (vector={vector!r}); refusing to fabricate a fallback"
            )
        return normalize(vector)
