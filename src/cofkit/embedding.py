from __future__ import annotations

from dataclasses import dataclass, field
from math import atan2, cos, pi, sin, sqrt
from typing import Mapping

from .geometry import (
    Mat3,
    Vec3,
    add,
    centroid,
    dot,
    mat3_identity,
    matmul_vec,
    norm,
    normalize,
    rotation_from_frame_to_axes,
    scale,
    sub,
)
from .model import AssemblyState, MonomerSpec, Pose, ReactionTemplate
from .model import ReactionEvent
from .planner import TopologyHint
from .reactions import bridge_target_distance
from .search import AssignmentOutcome
from .single_node_topologies import resolve_single_node_topology_layout
from .single_node_topologies_3d import resolve_three_d_single_node_topology_layout


@dataclass(frozen=True)
class EmbeddingResult:
    state: AssemblyState
    metadata: Mapping[str, object] = field(default_factory=dict)


@dataclass(frozen=True)
class EmbeddingConfig:
    default_layer_spacing: float = 8.0
    default_lateral_span: float = 30.0
    bridge_target_distance: float = 1.4


def _resolve_supported_single_node_layout(topology_id: str):
    try:
        return resolve_single_node_topology_layout(topology_id)
    except (KeyError, ValueError):
        return resolve_three_d_single_node_topology_layout(topology_id)


class PeriodicEmbedder:
    """Builds a simple initial periodic placement for monomer instances."""

    def __init__(self, config: EmbeddingConfig | None = None):
        self.config = config or EmbeddingConfig()

    def embed(
        self,
        outcome: AssignmentOutcome,
        monomer_specs: Mapping[str, MonomerSpec],
        templates: Mapping[str, ReactionTemplate],
    ) -> EmbeddingResult:
        topology = outcome.assignment_plan.net_plan.topology
        self._validate_single_node_builder_support(topology, outcome, monomer_specs)
        if self._is_single_node_bipartite_case(topology, outcome, monomer_specs):
            return self._embed_single_node_bipartite(outcome, monomer_specs, templates)
        template_distances = [
            self._template_target_distance(templates[event.template_id]) for event in outcome.events
        ]
        target_distance = sum(template_distances) / len(template_distances) if template_distances else self.config.bridge_target_distance
        cell = self._build_cell(topology, target_distance, len(outcome.monomer_instances))
        centers = self._initial_centers(topology, cell, len(outcome.monomer_instances))

        monomer_poses: dict[str, Pose] = {}
        pose_details: dict[str, object] = {}
        for index, instance in enumerate(outcome.monomer_instances):
            spec = monomer_specs[instance.monomer_id]
            center = centers.get(instance.id, self._fallback_center(index, target_distance))
            rotation = self._orientation_for_instance(instance.id, center, centers, outcome, monomer_specs)
            monomer_poses[instance.id] = Pose(translation=center, rotation_matrix=rotation)
            pose_details[instance.id] = {
                "slot_id": instance.metadata.get("slot_id"),
                "translation": center,
                "rotation_matrix": rotation,
                "motif_count": len(spec.motifs),
            }
        monomer_poses = self._refine_poses_from_events(monomer_poses, outcome, monomer_specs, templates)
        for instance_id, pose in monomer_poses.items():
            pose_details[instance_id] = dict(pose_details[instance_id])
            pose_details[instance_id]["translation"] = pose.translation

        state = AssemblyState(
            cell=cell,
            monomer_poses=monomer_poses,
            stacking_state="disabled",
        )
        metadata = {
            "mode": "topology-guided" if topology is not None else "motif-guided",
            "topology": topology.id if topology is not None else None,
            "target_distance": target_distance,
            "cell_kind": self._cell_kind(topology),
            "stacking_enabled": False,
            "poses": pose_details,
        }
        return EmbeddingResult(state=state, metadata=metadata)

    def _embed_single_node_bipartite(
        self,
        outcome: AssignmentOutcome,
        monomer_specs: Mapping[str, MonomerSpec],
        templates: Mapping[str, ReactionTemplate],
    ) -> EmbeddingResult:
        topology = outcome.assignment_plan.net_plan.topology
        template_distances = [
            self._template_target_distance(templates[event.template_id]) for event in outcome.events
        ]
        bridge_target = sum(template_distances) / len(template_distances) if template_distances else self.config.bridge_target_distance
        instances = outcome.monomer_instances
        topology_layout = self._single_node_topology_layout(topology, outcome, monomer_specs)
        canonical_a = topology_layout.directions if topology_layout is not None else self._trigonal_directions(pi / 6.0)
        canonical_b = tuple(scale(direction, -1.0) for direction in canonical_a)

        pose_details: dict[str, object] = {}
        projected_offsets: dict[str, tuple[float, ...]] = {}
        rotations: dict[str, Mat3] = {}

        for instance, target_directions in zip(instances, (canonical_a, canonical_b)):
            spec = monomer_specs[instance.monomer_id]
            rotation, offsets = self._rotation_for_planar_motifs(spec, target_directions)
            rotations[instance.id] = rotation
            projected_offsets[instance.id] = offsets

        edge_reactive_site_distances: tuple[float, ...] = ()
        if topology_layout is not None:
            single_node_layout = self._generalized_single_node_layout(
                outcome,
                monomer_specs,
                rotations,
                bridge_target,
            )
            if single_node_layout is not None:
                cell, centers, edge_reactive_site_distances = single_node_layout
            else:
                nearest_neighbor_distance = max(
                    bridge_target + max(0.0, first_offset) + max(0.0, second_offset)
                    for first_offset, second_offset in zip(
                        projected_offsets[instances[0].id],
                        projected_offsets[instances[1].id],
                    )
                )
                span = sqrt(3.0) * nearest_neighbor_distance
                cell = (
                    (span, 0.0, 0.0),
                    (0.5 * span, sqrt(3.0) * 0.5 * span, 0.0),
                    (0.0, 0.0, self.config.default_layer_spacing),
                )
                centers = self._single_node_bipartite_centers(topology, cell)
                edge_reactive_site_distances = (nearest_neighbor_distance,)
        else:
            nearest_neighbor_distance = max(
                bridge_target + max(0.0, first_offset) + max(0.0, second_offset)
                for first_offset, second_offset in zip(
                    projected_offsets[instances[0].id],
                    projected_offsets[instances[1].id],
                )
            )
            cell = self._build_cell(topology, nearest_neighbor_distance, len(instances))
            centers = self._single_node_bipartite_centers(topology, cell)
            edge_reactive_site_distances = (nearest_neighbor_distance,)

        monomer_poses: dict[str, Pose] = {}
        for index, instance in enumerate(instances):
            center = centers[instance.id]
            rotation = rotations[instance.id]
            monomer_poses[instance.id] = Pose(translation=center, rotation_matrix=rotation)
            pose_details[instance.id] = {
                "slot_id": instance.metadata.get("slot_id"),
                "translation": center,
                "rotation_matrix": rotation,
                "motif_count": len(monomer_specs[instance.monomer_id].motifs),
                "radial_offsets": tuple(round(offset, 6) for offset in projected_offsets[instance.id]),
                "sublattice": "A" if index == 0 else "B",
            }

        state = AssemblyState(cell=cell, monomer_poses=monomer_poses, stacking_state="disabled")
        metadata = {
            "mode": "topology-guided",
            "topology": topology.id if topology is not None else None,
            "target_distance": bridge_target,
            "reactive_site_distance": sum(edge_reactive_site_distances) / len(edge_reactive_site_distances),
            "edge_reactive_site_distances": tuple(round(value, 6) for value in edge_reactive_site_distances),
            "cell_kind": self._single_node_bipartite_cell_kind(topology, cell),
            "stacking_enabled": False,
            "placement_mode": "single-node-bipartite",
            "topology_family": "single-node-2d" if topology_layout is not None else None,
            "poses": pose_details,
        }
        return EmbeddingResult(state=state, metadata=metadata)

    def _build_cell(
        self,
        topology: TopologyHint | None,
        target_distance: float,
        n_instances: int,
    ) -> tuple[Vec3, Vec3, Vec3]:
        span = max(self.config.default_lateral_span * 0.4, target_distance * max(6, n_instances * 2))
        if topology is not None and topology.id == "hcb":
            return (
                (span, 0.0, 0.0),
                (0.5 * span, sqrt(3.0) * 0.5 * span, 0.0),
                (0.0, 0.0, self.config.default_layer_spacing),
            )
        return (
            (span, 0.0, 0.0),
            (0.0, span, 0.0),
            (0.0, 0.0, self.config.default_layer_spacing),
        )

    def _initial_centers(
        self,
        topology: TopologyHint | None,
        cell: tuple[Vec3, Vec3, Vec3],
        n_instances: int,
    ) -> dict[str, Vec3]:
        centers: list[Vec3] = []
        if topology is not None and topology.id == "hcb" and n_instances == 2:
            a, b, _ = cell
            centers = [
                scale(add(a, b), 1.0 / 3.0),
                scale(add(scale(a, 2.0), scale(b, 2.0)), 1.0 / 3.0),
            ]
        elif topology is not None and topology.id == "sql" and n_instances == 2:
            a, b, _ = cell
            centers = [
                scale(add(a, b), 0.25),
                scale(add(scale(a, 3.0), scale(b, 3.0)), 0.25),
            ]
        else:
            centers = [self._fallback_center(i, self.config.bridge_target_distance * 4.0) for i in range(n_instances)]
        return {f"m{i+1}": center for i, center in enumerate(centers)}

    def _refine_poses_from_events(
        self,
        poses: Mapping[str, Pose],
        outcome: AssignmentOutcome,
        monomer_specs: Mapping[str, MonomerSpec],
        templates: Mapping[str, ReactionTemplate],
    ) -> dict[str, Pose]:
        refined = dict(poses)
        if not outcome.events:
            return refined

        for event in outcome.events:
            if len(event.participants) != 2:
                continue
            template = templates[event.template_id]
            target_distance = self._template_target_distance(template)
            first, second = event.participants
            pose1 = refined[first.monomer_instance_id]
            pose2 = refined[second.monomer_instance_id]
            motif1 = monomer_specs[first.monomer_id].motif_by_id(first.motif_id)
            motif2 = monomer_specs[second.monomer_id].motif_by_id(second.motif_id)
            origin1 = self._world_position(pose1, motif1.frame.origin)
            origin2 = self._world_position(pose2, motif2.frame.origin)
            midpoint = centroid((origin1, origin2))
            direction = self._direction_between(origin1, origin2)
            if origin1 == origin2:
                direction = self._event_direction(event, monomer_specs)
            desired1 = add(midpoint, scale(direction, -0.5 * target_distance))
            desired2 = add(midpoint, scale(direction, 0.5 * target_distance))
            refined[first.monomer_instance_id] = Pose(
                translation=add(pose1.translation, sub(desired1, origin1)),
                rotation_matrix=pose1.rotation_matrix,
            )
            refined[second.monomer_instance_id] = Pose(
                translation=add(pose2.translation, sub(desired2, origin2)),
                rotation_matrix=pose2.rotation_matrix,
            )
        return refined

    def _orientation_for_instance(
        self,
        instance_id: str,
        center: Vec3,
        centers: Mapping[str, Vec3],
        outcome: AssignmentOutcome,
        monomer_specs: Mapping[str, MonomerSpec],
    ) -> Mat3:
        motif_frames = []
        target_primary = (1.0, 0.0, 0.0)
        for event in outcome.events:
            for participant in event.participants:
                if participant.monomer_instance_id != instance_id:
                    continue
                motif = monomer_specs[participant.monomer_id].motif_by_id(participant.motif_id)
                motif_frames.append(motif.frame)
                for other in event.participants:
                    if other.monomer_instance_id != instance_id:
                        partner_center = centers.get(other.monomer_instance_id)
                        if partner_center is not None:
                            target_primary = self._direction_between(center, partner_center)
                            break
        if not motif_frames:
            return mat3_identity()

        frame = motif_frames[0]
        target_normal = frame.normal
        return rotation_from_frame_to_axes(frame, target_primary, target_normal)

    def _direction_between(self, center: Vec3, partner_center: Vec3) -> Vec3:
        delta = (partner_center[0] - center[0], partner_center[1] - center[1], partner_center[2] - center[2])
        if delta == (0.0, 0.0, 0.0):
            return (1.0, 0.0, 0.0)
        length = sqrt(delta[0] ** 2 + delta[1] ** 2 + delta[2] ** 2)
        return (delta[0] / length, delta[1] / length, delta[2] / length)

    def _event_direction(
        self,
        event: ReactionEvent,
        monomer_specs: Mapping[str, MonomerSpec],
    ) -> Vec3:
        first, second = event.participants
        motif1 = monomer_specs[first.monomer_id].motif_by_id(first.motif_id)
        motif2 = monomer_specs[second.monomer_id].motif_by_id(second.motif_id)
        avg = (
            motif1.frame.primary[0] - motif2.frame.primary[0],
            motif1.frame.primary[1] - motif2.frame.primary[1],
            motif1.frame.primary[2] - motif2.frame.primary[2],
        )
        length = sqrt(avg[0] ** 2 + avg[1] ** 2 + avg[2] ** 2)
        if length < 1e-6:
            return (1.0, 0.0, 0.0)
        return (avg[0] / length, avg[1] / length, avg[2] / length)

    def _template_target_distance(self, template: ReactionTemplate) -> float:
        return bridge_target_distance(template, default_bridge_distance=self.config.bridge_target_distance)

    def _fallback_center(self, index: int, spacing: float) -> Vec3:
        return (index * spacing, 0.0, 0.0)

    def _cell_kind(self, topology: TopologyHint | None) -> str:
        if topology is not None and topology.id == "hcb":
            return "hexagonal"
        return "orthogonal"

    def _is_single_node_bipartite_case(
        self,
        topology: TopologyHint | None,
        outcome: AssignmentOutcome,
        monomer_specs: Mapping[str, MonomerSpec],
    ) -> bool:
        if topology is None or len(outcome.monomer_instances) != 2:
            return False
        if len({instance.monomer_id for instance in outcome.monomer_instances}) != 2:
            return False
        node_definitions = int(topology.metadata.get("n_node_definitions", len(topology.node_coordination) or 0))
        if node_definitions != 1:
            return False
        connectivities = {len(monomer_specs[instance.monomer_id].motifs) for instance in outcome.monomer_instances}
        return len(connectivities) == 1

    def _validate_single_node_builder_support(
        self,
        topology: TopologyHint | None,
        outcome: AssignmentOutcome,
        monomer_specs: Mapping[str, MonomerSpec],
    ) -> None:
        if topology is None or len(outcome.monomer_instances) != 2:
            return
        if int(topology.metadata.get("n_node_definitions", len(topology.node_coordination) or 0)) != 1:
            return
        layout = _resolve_supported_single_node_layout(topology.id)
        if layout.placement_model != "p1-two-node":
            raise ValueError(
                f"topology {topology.id!r} requires an expanded single-node builder; direct engine support is "
                "currently limited to two-node P1 single-node nets"
            )

    def _single_node_bipartite_centers(
        self,
        topology: TopologyHint | None,
        cell: tuple[Vec3, Vec3, Vec3],
    ) -> dict[str, Vec3]:
        if topology is not None and topology.id == "hcb":
            a, b, _ = cell
            return {
                "m1": scale(add(a, b), 1.0 / 3.0),
                "m2": scale(add(scale(a, 2.0), scale(b, 2.0)), 1.0 / 3.0),
            }
        a, _, _ = cell
        return {"m1": (0.0, 0.0, 0.0), "m2": (norm(a) / 2.0, 0.0, 0.0)}

    def _generalized_single_node_layout(
        self,
        outcome: AssignmentOutcome,
        monomer_specs: Mapping[str, MonomerSpec],
        rotations: Mapping[str, Mat3],
        bridge_target: float,
    ) -> tuple[tuple[Vec3, Vec3, Vec3], dict[str, Vec3], tuple[float, ...]] | None:
        instances = outcome.monomer_instances
        if len(instances) != 2:
            return None

        rows: list[tuple[tuple[int, int, int], Vec3]] = []
        first_instance_id = instances[0].id
        second_instance_id = instances[1].id

        for event in outcome.events:
            if len(event.participants) != 2:
                continue
            participant_by_instance = {participant.monomer_instance_id: participant for participant in event.participants}
            first = participant_by_instance.get(first_instance_id)
            second = participant_by_instance.get(second_instance_id)
            if first is None or second is None:
                return None

            first_spec = monomer_specs[first.monomer_id]
            second_spec = monomer_specs[second.monomer_id]
            first_motif = first_spec.motif_by_id(first.motif_id)
            second_motif = second_spec.motif_by_id(second.motif_id)
            direction = self._single_node_bridge_direction(
                first_motif=first_motif,
                first_rotation=rotations[first_instance_id],
                second_motif=second_motif,
                second_rotation=rotations[second_instance_id],
            )
            edge_vector = self._single_node_bridge_edge_vector(
                first_motif=first_motif,
                first_rotation=rotations[first_instance_id],
                second_motif=second_motif,
                second_rotation=rotations[second_instance_id],
                direction=direction,
                bridge_target=bridge_target,
            )
            image_delta = (
                second.periodic_image[0] - first.periodic_image[0],
                second.periodic_image[1] - first.periodic_image[1],
                second.periodic_image[2] - first.periodic_image[2],
            )
            rows.append((image_delta, edge_vector))

        if len(rows) < 3:
            return None

        reference_image, reference_vector = rows[0]
        cell_a = None
        cell_b = None
        for image_delta, edge_vector in rows[1:]:
            relative_image = (
                image_delta[0] - reference_image[0],
                image_delta[1] - reference_image[1],
                image_delta[2] - reference_image[2],
            )
            delta_vector = sub(edge_vector, reference_vector)
            if relative_image[0] != 0 and relative_image[1] == 0 and cell_a is None:
                cell_a = scale(delta_vector, 1.0 / relative_image[0])
            elif relative_image[1] != 0 and relative_image[0] == 0 and cell_b is None:
                cell_b = scale(delta_vector, 1.0 / relative_image[1])

        if cell_a is None or cell_b is None:
            return None

        cell = (
            cell_a,
            cell_b,
            (0.0, 0.0, self.config.default_layer_spacing),
        )
        reference_offset = add(scale(cell_a, reference_image[0]), scale(cell_b, reference_image[1]))
        center_delta = sub(reference_vector, reference_offset)
        center_a = scale(add(cell_a, cell_b), 1.0 / 3.0)
        centers = {
            first_instance_id: center_a,
            second_instance_id: add(center_a, center_delta),
        }
        return cell, centers, tuple(norm(edge_vector) for _, edge_vector in rows)

    def _single_node_bridge_direction(
        self,
        *,
        first_motif,
        first_rotation: Mat3,
        second_motif,
        second_rotation: Mat3,
    ) -> Vec3:
        first_direction = self._safe_normalize(matmul_vec(first_rotation, self._motif_direction_seed(first_motif)))
        second_direction = scale(
            self._safe_normalize(matmul_vec(second_rotation, self._motif_direction_seed(second_motif))),
            -1.0,
        )
        combined = add(first_direction, second_direction)
        if norm(combined) < 1e-8:
            combined = first_direction
        return self._safe_normalize(combined)

    def _single_node_bridge_edge_vector(
        self,
        *,
        first_motif,
        first_rotation: Mat3,
        second_motif,
        second_rotation: Mat3,
        direction: Vec3,
        bridge_target: float,
    ) -> Vec3:
        first_origin = matmul_vec(first_rotation, first_motif.frame.origin)
        second_origin = matmul_vec(second_rotation, second_motif.frame.origin)
        return add(scale(direction, bridge_target), sub(first_origin, second_origin))

    def _motif_direction_seed(self, motif) -> Vec3:
        if norm(motif.frame.origin) > 1e-6:
            return motif.frame.origin
        return motif.frame.primary

    def _single_node_bipartite_cell_kind(
        self,
        topology: TopologyHint | None,
        cell: tuple[Vec3, Vec3, Vec3],
    ) -> str:
        first, second, _ = cell
        first_norm = norm(first)
        second_norm = norm(second)
        if first_norm < 1e-8 or second_norm < 1e-8:
            return "oblique"
        cosine = dot(first, second) / (first_norm * second_norm)
        if abs(first_norm - second_norm) < 1e-3 and abs(cosine) < 1e-3:
            return "square"
        if abs(first_norm - second_norm) < 1e-3 and abs(cosine - 0.5) < 1e-3:
            return "hexagonal"
        return "oblique"

    def _single_node_topology_layout(
        self,
        topology: TopologyHint | None,
        outcome: AssignmentOutcome,
        monomer_specs: Mapping[str, MonomerSpec],
    ):
        if topology is None or len(outcome.monomer_instances) != 2:
            return None
        if int(topology.metadata.get("n_node_definitions", len(topology.node_coordination) or 0)) != 1:
            return None
        connectivities = {len(monomer_specs[instance.monomer_id].motifs) for instance in outcome.monomer_instances}
        if len(connectivities) != 1:
            return None
        try:
            layout = resolve_single_node_topology_layout(topology.id)
        except (KeyError, ValueError):
            return None
        if layout.connectivity != next(iter(connectivities)):
            return None
        return layout

    def _trigonal_directions(self, start_angle: float) -> tuple[Vec3, Vec3, Vec3]:
        return tuple(
            (cos(start_angle + offset), sin(start_angle + offset), 0.0)
            for offset in (0.0, 2.0 * pi / 3.0, -2.0 * pi / 3.0)
        )  # type: ignore[return-value]

    def _rotation_for_planar_motifs(
        self,
        spec: MonomerSpec,
        target_directions: tuple[Vec3, ...],
    ) -> tuple[Mat3, tuple[float, ...]]:
        motif_vectors = [
            motif.frame.origin if norm(motif.frame.origin) > 1e-6 else motif.frame.primary
            for motif in spec.motifs
        ]
        if len(motif_vectors) != len(target_directions):
            return mat3_identity(), tuple(0.0 for _ in target_directions)

        local_angles = [atan2(vector[1], vector[0]) for vector in motif_vectors]
        target_angles = [atan2(direction[1], direction[0]) for direction in target_directions]
        diffs = [
            self._wrap_angle(target_angle - local_angle)
            for local_angle, target_angle in zip(local_angles, target_angles)
        ]
        best_angle = sum(diffs) / len(diffs)

        rotation = (
            (cos(best_angle), -sin(best_angle), 0.0),
            (sin(best_angle), cos(best_angle), 0.0),
            (0.0, 0.0, 1.0),
        )
        projected = tuple(
            dot(matmul_vec(rotation, vector), direction)
            for vector, direction in zip(motif_vectors, target_directions)
        )
        return rotation, projected

    def _wrap_angle(self, value: float) -> float:
        while value <= -pi:
            value += 2.0 * pi
        while value > pi:
            value -= 2.0 * pi
        return value

    def _world_position(self, pose: Pose, local_position: Vec3) -> Vec3:
        return add(pose.translation, matmul_vec(pose.rotation_matrix, local_position))

    def _safe_normalize(self, vector: Vec3) -> Vec3:
        if norm(vector) < 1e-8:
            return (1.0, 0.0, 0.0)
        return normalize(vector)
