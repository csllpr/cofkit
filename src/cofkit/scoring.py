from __future__ import annotations

from dataclasses import dataclass
from math import fabs
from typing import Mapping

from .geometry import Vec3, add, distance, dot, matmul_vec, normalize, safe_normalize, scale, sub
from .linkage_geometry import effective_motif_origin
from .model import AssemblyState, MonomerSpec, ReactionTemplate
from .reactions import bridge_target_distance
from .search import AssignmentOutcome


@dataclass(frozen=True)
class BridgeEventMetrics:
    event_id: str
    template_id: str
    target_distance: float
    actual_distance: float
    distance_residual: float
    planarity_residual: float
    alignment_residual: float
    normal_misalignment_residual: float
    normal_alignment: float
    total_residual: float


@dataclass(frozen=True)
class BridgeGeometryReport:
    """Residual aggregation contract: ``total_residual`` is the **sum** of the
    per-event ``total_residual`` values over all bridge events — never a mean.

    Ranking is the only place a mean is derived, and it does so by dividing
    the sum by the number of per-event metric entries
    (``cofkit.model.residual_ranking_key``); the validator likewise computes
    its own mean from the per-event metrics.  Consumers that change the event
    count (e.g. ``stacking.py`` duplicating a layer's events) must recompute
    the sum from the per-event metrics so the derived mean stays consistent.
    """

    event_metrics: tuple[BridgeEventMetrics, ...] = ()
    total_residual: float = 0.0


class CandidateScorer:
    """Computes per-bridge-event geometry residuals from assignment and geometry.

    The per-event :meth:`bridge_geometry_report` residuals are the live metric:
    they drive the continuous optimizer and the coarse structure validator, and
    ranking uses the mean per-bridge-event geometry residual (see
    ``cofkit.model.candidate_ranking_key``). :meth:`scoring_metadata` packages
    those residuals for candidate metadata; the exported
    ``bridge_geometry_residual`` is the sum over events per the
    :class:`BridgeGeometryReport` contract.
    """

    def scoring_metadata(
        self,
        outcome: AssignmentOutcome,
        state: AssemblyState,
        monomer_specs: Mapping[str, MonomerSpec],
        templates: Mapping[str, ReactionTemplate],
    ) -> dict[str, object]:
        topology = outcome.assignment_plan.net_plan.topology
        bridge_report = self.bridge_geometry_report(outcome, state, monomer_specs, templates)
        return {
            "n_unreacted_motifs": len(outcome.unreacted_motifs),
            "topology": topology.id if topology is not None else None,
            "bridge_geometry_residual": bridge_report.total_residual,
            "bridge_event_metrics": tuple(
                {
                    "event_id": metrics.event_id,
                    "template_id": metrics.template_id,
                    "target_distance": metrics.target_distance,
                    "actual_distance": metrics.actual_distance,
                    "distance_residual": metrics.distance_residual,
                    "planarity_residual": metrics.planarity_residual,
                    "alignment_residual": metrics.alignment_residual,
                    "normal_misalignment_residual": metrics.normal_misalignment_residual,
                    "normal_alignment": metrics.normal_alignment,
                    "total_residual": metrics.total_residual,
                }
                for metrics in bridge_report.event_metrics
            ),
        }

    def bridge_geometry_report(
        self,
        outcome: AssignmentOutcome,
        state: AssemblyState,
        monomer_specs: Mapping[str, MonomerSpec],
        templates: Mapping[str, ReactionTemplate],
    ) -> BridgeGeometryReport:
        event_metrics: list[BridgeEventMetrics] = []
        total_residual = 0.0
        for event in outcome.events:
            template = templates[event.template_id]
            if len(event.participants) == 2:
                first, second = event.participants
                pose1 = state.monomer_poses[first.monomer_instance_id]
                pose2 = state.monomer_poses[second.monomer_instance_id]
                motif1 = monomer_specs[first.monomer_id].motif_by_id(first.motif_id)
                motif2 = monomer_specs[second.monomer_id].motif_by_id(second.motif_id)
                origin1 = self._world_motif_origin(
                    state.cell,
                    pose1.translation,
                    pose1.rotation_matrix,
                    effective_motif_origin(event.template_id, monomer_specs[first.monomer_id], motif1),
                    first.periodic_image,
                )
                origin2 = self._world_motif_origin(
                    state.cell,
                    pose2.translation,
                    pose2.rotation_matrix,
                    effective_motif_origin(event.template_id, monomer_specs[second.monomer_id], motif2),
                    second.periodic_image,
                )
                separation = distance(origin1, origin2)
                target = self._target_distance(template)
                distance_residual = fabs(separation - target)

                normal1 = matmul_vec(pose1.rotation_matrix, motif1.frame.normal)
                normal2 = matmul_vec(pose2.rotation_matrix, motif2.frame.normal)
                normal_alignment = max(-1.0, min(1.0, dot(normalize(normal1), normalize(normal2))))

                bridge_vector = sub(origin2, origin1)
                unit_vector = self._safe_normalize(bridge_vector)
                primary1 = matmul_vec(pose1.rotation_matrix, motif1.frame.primary)
                primary2 = matmul_vec(pose2.rotation_matrix, motif2.frame.primary)
                alignment_residual = max(0.0, 1.0 - dot(self._safe_normalize(primary1), unit_vector))
                alignment_residual += max(0.0, 1.0 - dot(self._safe_normalize(primary2), self._invert(unit_vector)))

                planarity_residual = fabs(dot(unit_vector, self._safe_normalize(normal1)))
                planarity_residual += fabs(dot(unit_vector, self._safe_normalize(normal2)))
                # This is only a torsion-style surrogate from motif-frame normals. It rewards
                # consistent local bridge planes for planar/restricted templates without claiming
                # any atomistic dihedral relaxation or force-field meaning.
                normal_misalignment_residual = self._normal_misalignment_residual(template, normal_alignment)

                total_event_residual = (
                    distance_residual
                    + 0.5 * planarity_residual
                    + 0.5 * alignment_residual
                    + normal_misalignment_residual
                )
                total_residual += total_event_residual
                event_metrics.append(
                    BridgeEventMetrics(
                        event_id=event.id,
                        template_id=event.template_id,
                        target_distance=target,
                        actual_distance=separation,
                        distance_residual=distance_residual,
                        planarity_residual=planarity_residual,
                        alignment_residual=alignment_residual,
                        normal_misalignment_residual=normal_misalignment_residual,
                        normal_alignment=normal_alignment,
                        total_residual=total_event_residual,
                    )
                )
        return BridgeGeometryReport(
            event_metrics=tuple(event_metrics),
            total_residual=total_residual,
        )

    def _target_distance(self, template: ReactionTemplate) -> float:
        return bridge_target_distance(template)

    def _world_motif_origin(
        self,
        cell: tuple[Vec3, Vec3, Vec3],
        translation: tuple[float, float, float],
        rotation: tuple[tuple[float, float, float], ...],
        local_origin: tuple[float, float, float],
        periodic_image: tuple[int, int, int],
    ) -> tuple[float, float, float]:
        rotated = matmul_vec(rotation, local_origin)
        imaged = (
            translation[0] + rotated[0],
            translation[1] + rotated[1],
            translation[2] + rotated[2],
        )
        return add(
            imaged,
            add(
                add(scale(cell[0], periodic_image[0]), scale(cell[1], periodic_image[1])),
                scale(cell[2], periodic_image[2]),
            ),
        )

    def _safe_normalize(self, vector: Vec3) -> Vec3:
        return safe_normalize(vector, fallback=(1.0, 0.0, 0.0))

    def _invert(self, vector: Vec3) -> Vec3:
        return (-vector[0], -vector[1], -vector[2])

    def _normal_misalignment_residual(self, template: ReactionTemplate, normal_alignment: float) -> float:
        return self._normal_misalignment_weight(template) * 0.5 * (1.0 - normal_alignment)

    def _normal_misalignment_weight(self, template: ReactionTemplate) -> float:
        if template.planarity_prior == "planar" and template.torsion_prior in {"restricted", "locked"}:
            return 1.0
        if template.planarity_prior in {"planar", "semi-planar"} or template.torsion_prior in {"moderate", "restricted"}:
            return 0.35
        return 0.0
