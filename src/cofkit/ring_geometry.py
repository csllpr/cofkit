from __future__ import annotations

from dataclasses import dataclass, field
from math import atan2, pi
from typing import Mapping

from .geometry import Vec3, add, cross, dot, matmul_vec, norm, normalize, scale, sub
from .model import AssemblyState, MonomerSpec, Pose, ReactionEvent


@dataclass(frozen=True)
class RingGeometryProfile:
    template_id: str
    ring_atom_bond_length: float
    radial_tolerance: float = 0.18
    planarity_tolerance: float = 0.12
    angular_tolerance_degrees: float = 8.0

    @property
    def participant_radius(self) -> float:
        # Both supported products are regular, alternating six-membered rings.
        return self.ring_atom_bond_length


@dataclass(frozen=True)
class RingEventGeometry:
    event_id: str
    radial_rms: float
    radial_max: float
    planarity_rms: float
    planarity_max: float
    angular_rms_degrees: float
    angular_max_degrees: float
    participant_instances: tuple[str, ...]

    @property
    def residual(self) -> float:
        return self.radial_rms + self.planarity_rms + self.angular_rms_degrees / 30.0


@dataclass(frozen=True)
class RingGeometryReport:
    event_metrics: tuple[RingEventGeometry, ...]
    total_residual: float

    def as_dict(self) -> dict[str, object]:
        return {
            "total_residual": self.total_residual,
            "events": [
                {
                    "event_id": metric.event_id,
                    "radial_rms": metric.radial_rms,
                    "radial_max": metric.radial_max,
                    "planarity_rms": metric.planarity_rms,
                    "planarity_max": metric.planarity_max,
                    "angular_rms_degrees": metric.angular_rms_degrees,
                    "angular_max_degrees": metric.angular_max_degrees,
                    "participant_instances": metric.participant_instances,
                }
                for metric in self.event_metrics
            ],
        }


@dataclass(frozen=True)
class RingValidationResult:
    classification: str
    reasons: tuple[str, ...]
    metrics: Mapping[str, object] = field(default_factory=dict)


@dataclass(frozen=True)
class RingOptimizationResult:
    state: AssemblyState
    metrics: Mapping[str, object] = field(default_factory=dict)


def ring_geometry_profile(template_id: str) -> RingGeometryProfile:
    try:
        return {
            "boroxine_trimerization": RingGeometryProfile(template_id, 1.38),
            "triazine_trimerization": RingGeometryProfile(template_id, 1.35),
        }[template_id]
    except KeyError as exc:
        raise KeyError(f"no ring geometry profile for template {template_id!r}") from exc


def fractional_to_cartesian(fractional: Vec3, cell: tuple[Vec3, Vec3, Vec3]) -> Vec3:
    return add(add(scale(cell[0], fractional[0]), scale(cell[1], fractional[1])), scale(cell[2], fractional[2]))


def motif_world_position(
    event_ref,
    state: AssemblyState,
    monomer_specs: Mapping[str, MonomerSpec],
) -> Vec3:
    motif = monomer_specs[event_ref.monomer_id].motif_by_id(event_ref.motif_id)
    pose = state.monomer_poses[event_ref.monomer_instance_id]
    image = fractional_to_cartesian(tuple(float(v) for v in event_ref.periodic_image), state.cell)
    return add(add(matmul_vec(pose.rotation_matrix, motif.frame.origin), pose.translation), image)


def ring_geometry_report(
    events: tuple[ReactionEvent, ...],
    state: AssemblyState,
    monomer_specs: Mapping[str, MonomerSpec],
) -> RingGeometryReport:
    metrics: list[RingEventGeometry] = []
    for event in events:
        if len(event.participants) != 3 or "ring_center_fractional" not in event.metadata:
            continue
        profile = ring_geometry_profile(event.template_id)
        center_fractional = tuple(float(v) for v in event.metadata["ring_center_fractional"])
        center = fractional_to_cartesian(center_fractional, state.cell)
        center = add(center, _event_ring_center_offset(event))
        normal = normalize(tuple(float(v) for v in event.metadata.get("ring_normal", (0.0, 0.0, 1.0))))
        points = tuple(motif_world_position(ref, state, monomer_specs) for ref in event.participants)
        radial_errors: list[float] = []
        plane_errors: list[float] = []
        projected: list[Vec3] = []
        for point in points:
            offset = sub(point, center)
            plane_error = dot(offset, normal)
            in_plane = sub(offset, scale(normal, plane_error))
            radial_errors.append(norm(in_plane) - profile.participant_radius)
            plane_errors.append(plane_error)
            projected.append(normalize(in_plane))

        basis_x = projected[0]
        basis_y = normalize(cross(normal, basis_x))
        angles = sorted((atan2(dot(v, basis_y), dot(v, basis_x)) % (2.0 * pi)) for v in projected)
        gaps = tuple(
            ((angles[(index + 1) % 3] - angles[index]) % (2.0 * pi))
            for index in range(3)
        )
        angular_errors = tuple((gap - 2.0 * pi / 3.0) * 180.0 / pi for gap in gaps)
        metric = RingEventGeometry(
            event_id=event.id,
            radial_rms=_rms(radial_errors),
            radial_max=max(abs(value) for value in radial_errors),
            planarity_rms=_rms(plane_errors),
            planarity_max=max(abs(value) for value in plane_errors),
            angular_rms_degrees=_rms(angular_errors),
            angular_max_degrees=max(abs(value) for value in angular_errors),
            participant_instances=tuple(
                f"{ref.monomer_instance_id}@{ref.periodic_image}" for ref in event.participants
            ),
        )
        metrics.append(metric)
    return RingGeometryReport(tuple(metrics), sum(metric.residual for metric in metrics))


class RingGeometryOptimizer:
    """Refines shared precursor translations against all incident ring constraints."""

    def __init__(self, max_iterations: int = 6, translation_step: float = 0.5):
        self.max_iterations = max_iterations
        self.translation_step = translation_step

    def optimize(
        self,
        events: tuple[ReactionEvent, ...],
        state: AssemblyState,
        monomer_specs: Mapping[str, MonomerSpec],
    ) -> RingOptimizationResult:
        initial = ring_geometry_report(events, state, monomer_specs)
        best_state = state
        best = initial
        accepted = 0
        for _ in range(self.max_iterations):
            proposal = self._translation_proposal(events, best_state, monomer_specs)
            report = ring_geometry_report(events, proposal, monomer_specs)
            if report.total_residual + 1e-10 < best.total_residual:
                best_state, best = proposal, report
                accepted += 1
            else:
                break
        return RingOptimizationResult(
            state=best_state,
            metrics={
                "enabled": True,
                "iterations": self.max_iterations,
                "accepted_iterations": accepted,
                "initial_residual": initial.total_residual,
                "final_residual": best.total_residual,
                "improved": best.total_residual + 1e-10 < initial.total_residual,
            },
        )

    def _translation_proposal(
        self,
        events: tuple[ReactionEvent, ...],
        state: AssemblyState,
        monomer_specs: Mapping[str, MonomerSpec],
    ) -> AssemblyState:
        updates: dict[str, Vec3] = {instance_id: (0.0, 0.0, 0.0) for instance_id in state.monomer_poses}
        counts: dict[str, int] = {instance_id: 0 for instance_id in state.monomer_poses}
        for event in events:
            if len(event.participants) != 3 or "ring_center_fractional" not in event.metadata:
                continue
            profile = ring_geometry_profile(event.template_id)
            center = add(
                fractional_to_cartesian(tuple(event.metadata["ring_center_fractional"]), state.cell),
                _event_ring_center_offset(event),
            )
            normal = normalize(tuple(event.metadata.get("ring_normal", (0.0, 0.0, 1.0))))
            for ref in event.participants:
                point = motif_world_position(ref, state, monomer_specs)
                offset = sub(point, center)
                in_plane = sub(offset, scale(normal, dot(offset, normal)))
                if norm(in_plane) < 1e-10:
                    continue
                desired = add(center, scale(normalize(in_plane), profile.participant_radius))
                correction = sub(desired, point)
                instance_id = ref.monomer_instance_id
                updates[instance_id] = add(updates[instance_id], correction)
                counts[instance_id] += 1
        poses: dict[str, Pose] = {}
        for instance_id, pose in state.monomer_poses.items():
            count = counts[instance_id]
            update = updates[instance_id] if count == 0 else scale(updates[instance_id], self.translation_step / count)
            poses[instance_id] = Pose(add(pose.translation, update), pose.rotation_matrix)
        return AssemblyState(
            cell=state.cell,
            monomer_poses=poses,
            torsions=state.torsions,
            layer_offsets=state.layer_offsets,
            stacking_state=state.stacking_state,
        )


def validate_ring_geometry(
    events: tuple[ReactionEvent, ...],
    state: AssemblyState,
    monomer_specs: Mapping[str, MonomerSpec],
) -> RingValidationResult:
    report = ring_geometry_report(events, state, monomer_specs)
    reasons: list[str] = []
    for event, metric in zip((event for event in events if len(event.participants) == 3), report.event_metrics):
        profile = ring_geometry_profile(event.template_id)
        if len(set(metric.participant_instances)) != 3:
            reasons.append(f"{event.id}: ring participants are not three distinct monomer instances")
        if metric.radial_max > profile.radial_tolerance:
            reasons.append(f"{event.id}: radial residual {metric.radial_max:.3f} A exceeds {profile.radial_tolerance:.3f} A")
        if metric.planarity_max > profile.planarity_tolerance:
            reasons.append(f"{event.id}: planarity residual {metric.planarity_max:.3f} A exceeds {profile.planarity_tolerance:.3f} A")
        if metric.angular_max_degrees > profile.angular_tolerance_degrees:
            reasons.append(
                f"{event.id}: angular residual {metric.angular_max_degrees:.2f} deg exceeds "
                f"{profile.angular_tolerance_degrees:.2f} deg"
            )
    if len(report.event_metrics) != len(events):
        reasons.append("one or more events lack ring geometry metadata")
    return RingValidationResult(
        classification="accepted" if not reasons else "rejected",
        reasons=tuple(reasons),
        metrics=report.as_dict(),
    )


def _rms(values) -> float:
    values = tuple(float(value) for value in values)
    return (sum(value * value for value in values) / len(values)) ** 0.5 if values else 0.0


def _event_ring_center_offset(event: ReactionEvent) -> Vec3:
    raw_offset = event.metadata.get("ring_center_offset_cartesian", (0.0, 0.0, 0.0))
    if not isinstance(raw_offset, (tuple, list)) or len(raw_offset) != 3:
        raise ValueError(f"event {event.id!r} has invalid ring_center_offset_cartesian metadata")
    return tuple(float(value) for value in raw_offset)  # type: ignore[return-value]
