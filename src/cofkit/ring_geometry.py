from __future__ import annotations

from dataclasses import dataclass, field
from math import acos, atan2, degrees, pi
from typing import Mapping

from .geometry import ANGULAR_RESIDUAL_DOWN_WEIGHT, Vec3, add, cross, dot, matmul_vec, norm, normalize, scale, sub
from .model import AssemblyState, MonomerSpec, Pose, ReactionEvent

# Cited from standard aromatic ring bond lengths: boroxine B-O ~1.38
# angstrom (boroxine B3O3 ring), triazine C-N ~1.35 angstrom
# (1,3,5-triazine ring). This module owns ring-template geometry, so these
# are the single owner for both the ring profiles below and the
# bridge_target_distance / default ring radius consumers in reactions.py
# and reaction_realization.py.
BOROXINE_BO_BOND_LENGTH = 1.38
TRIAZINE_CN_BOND_LENGTH = 1.35

# Numerical guard for "did the residual improve" comparisons.
_RESIDUAL_IMPROVEMENT_EPSILON = 1e-10

# Derived: the supported ring products place a trigonal sp2 participant atom
# (boroxine B, triazine C) at each ring vertex. Its three coplanar bond
# directions are spaced 360/3 = 120 degrees, so the exocyclic ring-monomer
# attachment bond should sit ~120 degrees from each of the two ring bonds.
RING_ATTACHMENT_IDEAL_ANGLE_DEGREES = 120.0


@dataclass(frozen=True)
class RingGeometryProfile:
    template_id: str
    ring_atom_bond_length: float
    # Heuristic — pending calibration: acceptance tolerances for the ring
    # geometry validation (validate_ring_geometry): max radial deviation of a
    # participant from the ring radius (A), max out-of-plane deviation (A),
    # and max deviation from the ideal 120-degree angular gap.
    radial_tolerance: float = 0.18
    planarity_tolerance: float = 0.12
    angular_tolerance_degrees: float = 8.0
    # Heuristic — pending calibration: exocyclic attachment acceptance
    # tolerances for validate_ring_geometry's attachment channel. The
    # deviation is |measured - RING_ATTACHMENT_IDEAL_ANGLE_DEGREES| for the
    # two angles between the exocyclic bond and the ring bonds at each ring
    # participant. Accepted aromatic sp2 attachments sit within ~1 degree of
    # ideal after assembly, so 15 degrees already indicates strained seed
    # geometry that a relaxation pass might repair -> warning tier. Beyond
    # 40 degrees the precursor's motif spacing is incompatible with a regular
    # ring (the reproduced ortho-diboronic-acid pathology deviates by ~55
    # degrees: 64.3/175.3 degrees instead of 120/120) and is classified as
    # rejection, not repairable strain.
    attachment_warning_deviation_degrees: float = 15.0
    attachment_rejection_deviation_degrees: float = 40.0
    # Heuristic — pending calibration: exocyclic bond-length tolerances as a
    # fraction of the precursor's own internal reactive-to-anchor bond length
    # (the target is derived per attachment from the monomer geometry, not an
    # external constant). The current realizers keep this distance almost
    # exactly (they move the participant atom by at most the arrangement
    # radial tolerance), so these guard future realizers/placements that move
    # participants further.
    attachment_warning_distance_fraction: float = 0.15
    attachment_rejection_distance_fraction: float = 0.35

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
        return self.radial_rms + self.planarity_rms + self.angular_rms_degrees / ANGULAR_RESIDUAL_DOWN_WEIGHT


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
class RingAttachmentMeasurement:
    """Exocyclic attachment geometry at one ring participant.

    The ring bonds and the exocyclic bond are evaluated at the participant's
    realized position (the ring vertex the realizer moves the reactive atom
    to), so the angles describe the geometry the CIF export writes. The two
    ``ring_angles_degrees`` are the angles between the exocyclic bond and each
    of the two ring bonds at the participant; both should sit near
    ``RING_ATTACHMENT_IDEAL_ANGLE_DEGREES`` for an unstrained sp2 attachment.
    """

    event_id: str
    participant: str
    exocyclic_distance: float
    precursor_distance: float
    ring_angles_degrees: tuple[float, float]


# Attachment report statuses, aligned with the validation coverage
# vocabulary (cofkit.validation): the measurement ran for every participant
# ("measured"), ran partially or not at all because atom metadata/positions
# were absent ("missing_data"), or the events carry no ring geometry at all
# ("not_applicable"). Missing data is reported explicitly and never counted
# as a pass.
ATTACHMENT_MEASURED = "measured"
ATTACHMENT_MISSING_DATA = "missing_data"
ATTACHMENT_NOT_APPLICABLE = "not_applicable"


@dataclass(frozen=True)
class RingAttachmentReport:
    status: str
    classification: str
    measurements: tuple[RingAttachmentMeasurement, ...] = ()
    reasons: tuple[str, ...] = ()
    n_unmeasured_participants: int = 0

    def as_dict(self) -> dict[str, object]:
        return {
            "status": self.status,
            "classification": self.classification,
            "reasons": self.reasons,
            "n_unmeasured_participants": self.n_unmeasured_participants,
            "measurements": [
                {
                    "event_id": measurement.event_id,
                    "participant": measurement.participant,
                    "exocyclic_distance": measurement.exocyclic_distance,
                    "precursor_distance": measurement.precursor_distance,
                    "ring_angles_degrees": measurement.ring_angles_degrees,
                }
                for measurement in self.measurements
            ],
        }


@dataclass(frozen=True)
class RingValidationResult:
    classification: str
    reasons: tuple[str, ...]
    metrics: Mapping[str, object] = field(default_factory=dict)
    # Ring-arrangement and attachment verdicts are kept separately
    # interpretable: the arrangement verdict covers the ring participants
    # themselves (radial/planarity/angular-gap placement), the attachment
    # verdict the exocyclic ring-monomer bonds/angles. `classification` is
    # the combined verdict: "rejected" if either channel rejects, "warning"
    # if the attachment channel warns, otherwise "accepted". Passing either
    # channel is a geometric screening result only — it does not certify
    # equilibrium geometry.
    arrangement_classification: str = "accepted"
    attachment_classification: str = "not_applicable"
    attachment: RingAttachmentReport | None = None


@dataclass(frozen=True)
class RingOptimizationResult:
    state: AssemblyState
    metrics: Mapping[str, object] = field(default_factory=dict)


def ring_geometry_profile(template_id: str) -> RingGeometryProfile:
    try:
        return {
            "boroxine_trimerization": RingGeometryProfile(template_id, BOROXINE_BO_BOND_LENGTH),
            "triazine_trimerization": RingGeometryProfile(template_id, TRIAZINE_CN_BOND_LENGTH),
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


def _atom_world_position(
    event_ref,
    state: AssemblyState,
    monomer_specs: Mapping[str, MonomerSpec],
    atom_id: int,
) -> Vec3:
    spec = monomer_specs[event_ref.monomer_id]
    pose = state.monomer_poses[event_ref.monomer_instance_id]
    image = fractional_to_cartesian(tuple(float(v) for v in event_ref.periodic_image), state.cell)
    return add(add(matmul_vec(pose.rotation_matrix, spec.atom_positions[atom_id]), pose.translation), image)


def _angle_degrees(first: Vec3, second: Vec3) -> float:
    cosine = dot(first, second) / (norm(first) * norm(second))
    return degrees(acos(max(-1.0, min(1.0, cosine))))


def ring_attachment_report(
    events: tuple[ReactionEvent, ...],
    state: AssemblyState,
    monomer_specs: Mapping[str, MonomerSpec],
) -> RingAttachmentReport:
    """Measure the exocyclic attachment geometry of each ring participant.

    For every three-participant ring event, the participant reactive atom
    (boroxine B, triazine C) is placed by the realizer at the ring vertex
    ``center + radius * radial_direction``, and the two ring neighbours at
    ``center + radius * normalize(radial_i + radial_{i+1})`` — the same
    construction ``reaction_realization`` uses — so the angles between the
    exocyclic anchor bond and each ring bond describe the exported geometry.
    Participants whose motifs lack ``reactive_atom_id``/``anchor_atom_id``
    metadata or whose monomer has no atom positions are counted as unmeasured
    (status ``missing_data``); unmeasured participants are never treated as
    passing.
    """
    measurements: list[RingAttachmentMeasurement] = []
    n_unmeasured = 0
    n_ring_events = 0
    for event in events:
        if len(event.participants) != 3 or "ring_center_fractional" not in event.metadata:
            continue
        n_ring_events += 1
        profile = ring_geometry_profile(event.template_id)
        center = fractional_to_cartesian(tuple(float(v) for v in event.metadata["ring_center_fractional"]), state.cell)
        center = add(center, _event_ring_center_offset(event))
        normal = normalize(tuple(float(v) for v in event.metadata.get("ring_normal", (0.0, 0.0, 1.0))))
        radius = float(event.metadata.get("ring_atom_bond_length", profile.ring_atom_bond_length))
        refs = event.participants
        radial_directions: list[Vec3] = []
        participant_atom_ids: list[tuple[int, int] | None] = []
        for ref in refs:
            motif = monomer_specs[ref.monomer_id].motif_by_id(ref.motif_id)
            reactive_atom_id = motif.metadata.get("reactive_atom_id")
            anchor_atom_id = motif.metadata.get("anchor_atom_id")
            spec = monomer_specs[ref.monomer_id]
            if (
                not isinstance(reactive_atom_id, int)
                or not isinstance(anchor_atom_id, int)
                or not spec.atom_positions
            ):
                participant_atom_ids.append(None)
            else:
                participant_atom_ids.append((reactive_atom_id, anchor_atom_id))
            atom_id = reactive_atom_id if isinstance(reactive_atom_id, int) else motif.atom_ids[0]
            if spec.atom_positions:
                current = _atom_world_position(ref, state, monomer_specs, atom_id)
            else:
                current = motif_world_position(ref, state, monomer_specs)
            offset = sub(current, center)
            in_plane = sub(offset, scale(normal, dot(offset, normal)))
            if norm(in_plane) < 1e-10:
                participant_atom_ids[-1] = None
                radial_directions.append((0.0, 0.0, 0.0))
                continue
            radial_directions.append(normalize(in_plane))
        targets = tuple(add(center, scale(direction, radius)) for direction in radial_directions)
        intermediates: list[Vec3] = []
        for index in range(3):
            direction_sum = add(radial_directions[index], radial_directions[(index + 1) % 3])
            if norm(direction_sum) < 1e-10:
                intermediates.append(targets[index])
                continue
            intermediates.append(add(center, scale(normalize(direction_sum), radius)))
        for index, ref in enumerate(refs):
            atom_ids = participant_atom_ids[index]
            participant_label = f"{ref.monomer_instance_id}@{ref.periodic_image}"
            if atom_ids is None:
                n_unmeasured += 1
                continue
            reactive_atom_id, anchor_atom_id = atom_ids
            spec = monomer_specs[ref.monomer_id]
            anchor_world = _atom_world_position(ref, state, monomer_specs, anchor_atom_id)
            exocyclic = sub(anchor_world, targets[index])
            ring_vectors = (
                sub(intermediates[index], targets[index]),
                sub(intermediates[(index - 1) % 3], targets[index]),
            )
            if norm(exocyclic) < 1e-10 or min(norm(vector) for vector in ring_vectors) < 1e-10:
                n_unmeasured += 1
                continue
            measurements.append(
                RingAttachmentMeasurement(
                    event_id=event.id,
                    participant=participant_label,
                    exocyclic_distance=norm(exocyclic),
                    precursor_distance=norm(
                        sub(spec.atom_positions[anchor_atom_id], spec.atom_positions[reactive_atom_id])
                    ),
                    ring_angles_degrees=tuple(
                        _angle_degrees(exocyclic, vector) for vector in ring_vectors
                    ),
                )
            )
    if n_ring_events == 0:
        return RingAttachmentReport(status=ATTACHMENT_NOT_APPLICABLE, classification="not_applicable")
    status = ATTACHMENT_MEASURED if n_unmeasured == 0 else ATTACHMENT_MISSING_DATA
    if not measurements:
        # Nothing could be measured: report honestly, never as a pass.
        return RingAttachmentReport(
            status=status,
            classification="not_applicable",
            n_unmeasured_participants=n_unmeasured,
        )
    template_by_event_id = {event.id: event.template_id for event in events}
    reasons: list[str] = []
    classification = "accepted"
    for measurement in measurements:
        profile = ring_geometry_profile(template_by_event_id[measurement.event_id])
        worst_deviation = max(
            abs(angle - RING_ATTACHMENT_IDEAL_ANGLE_DEGREES) for angle in measurement.ring_angles_degrees
        )
        distance_fraction = (
            abs(measurement.exocyclic_distance - measurement.precursor_distance) / measurement.precursor_distance
            if measurement.precursor_distance > 1e-10
            else 0.0
        )
        label = f"{measurement.event_id} {measurement.participant}"
        angles_text = (
            f"exocyclic attachment angles "
            f"{measurement.ring_angles_degrees[0]:.2f}/{measurement.ring_angles_degrees[1]:.2f} deg "
            f"deviate {worst_deviation:.2f} deg from the ideal "
            f"{RING_ATTACHMENT_IDEAL_ANGLE_DEGREES:.1f} deg"
        )
        distance_text = (
            f"exocyclic bond length {measurement.exocyclic_distance:.3f} A deviates "
            f"{100.0 * distance_fraction:.1f}% from the precursor {measurement.precursor_distance:.3f} A"
        )
        if worst_deviation > profile.attachment_rejection_deviation_degrees:
            reasons.append(
                f"{label}: {angles_text}, beyond the rejection tolerance "
                f"{profile.attachment_rejection_deviation_degrees:.2f} deg; the precursor attachment "
                "geometry is incompatible with a regular ring"
            )
            classification = "rejected"
        elif worst_deviation > profile.attachment_warning_deviation_degrees:
            reasons.append(
                f"{label}: {angles_text}, beyond the warning tolerance "
                f"{profile.attachment_warning_deviation_degrees:.2f} deg; strained seed attachment geometry"
            )
            if classification == "accepted":
                classification = "warning"
        if distance_fraction > profile.attachment_rejection_distance_fraction:
            reasons.append(
                f"{label}: {distance_text}, beyond the rejection fraction "
                f"{profile.attachment_rejection_distance_fraction:.2f}"
            )
            classification = "rejected"
        elif distance_fraction > profile.attachment_warning_distance_fraction:
            reasons.append(
                f"{label}: {distance_text}, beyond the warning fraction "
                f"{profile.attachment_warning_distance_fraction:.2f}"
            )
            if classification == "accepted":
                classification = "warning"
    return RingAttachmentReport(
        status=status,
        classification=classification,
        measurements=tuple(measurements),
        reasons=tuple(reasons),
        n_unmeasured_participants=n_unmeasured,
    )


class RingGeometryOptimizer:
    """Refines shared precursor translations against all incident ring constraints."""

    # Heuristic — pending calibration: small greedy budget and per-step
    # damping for the translation proposals; the loop breaks early on the
    # first non-improving proposal, so the cap is rarely reached.
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
            if report.total_residual + _RESIDUAL_IMPROVEMENT_EPSILON < best.total_residual:
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
                "improved": best.total_residual + _RESIDUAL_IMPROVEMENT_EPSILON < initial.total_residual,
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
    arrangement_classification = "accepted" if not reasons else "rejected"
    attachment = ring_attachment_report(events, state, monomer_specs)
    if attachment.classification == "rejected" or arrangement_classification == "rejected":
        classification = "rejected"
    elif attachment.classification == "warning":
        classification = "warning"
    else:
        classification = "accepted"
    return RingValidationResult(
        classification=classification,
        reasons=tuple(reasons) + attachment.reasons,
        metrics=report.as_dict(),
        arrangement_classification=arrangement_classification,
        attachment_classification=attachment.classification,
        attachment=attachment,
    )


def _rms(values) -> float:
    values = tuple(float(value) for value in values)
    return (sum(value * value for value in values) / len(values)) ** 0.5 if values else 0.0


def _event_ring_center_offset(event: ReactionEvent) -> Vec3:
    raw_offset = event.metadata.get("ring_center_offset_cartesian", (0.0, 0.0, 0.0))
    if not isinstance(raw_offset, (tuple, list)) or len(raw_offset) != 3:
        raise ValueError(f"event {event.id!r} has invalid ring_center_offset_cartesian metadata")
    return tuple(float(value) for value in raw_offset)  # type: ignore[return-value]
