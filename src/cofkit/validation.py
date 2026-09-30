from __future__ import annotations

import json
import os
import shutil
import sys
from collections import Counter, defaultdict
from concurrent.futures import FIRST_COMPLETED, Future, ThreadPoolExecutor, wait
from dataclasses import dataclass, field
from math import radians, sin
from pathlib import Path
from typing import Any, Mapping

try:  # pragma: no cover - exercised in integration environments
    import gemmi
except ImportError:  # pragma: no cover - import guard for incomplete environments
    gemmi = None

from .cif_checks import cif_value_str
from .periodic_geometry import images_within, p1_shift
from .reactions import linkage_profile
from .topologies import get_topology_hint
from .vdw import (
    DEFAULT_HARD_MIN_NONBONDED_PLAIN_DISTANCE,
    DEFAULT_NONBONDED_SEARCH_RADIUS,
    DEFAULT_SEVERE_OVERLAP_BONDED_DISTANCE,
    DEFAULT_VDW_CLASH_RATIO,
    assess_pair,
)


@dataclass(frozen=True)
class CoarseValidationThresholds:
    # Bridge residuals are |measured - target| inter-monomer linkage bond
    # distances in angstrom, measured from the final exported CIF coordinates
    # (assembly-time seed metrics are informational only); ratios are
    # measured / target (dimensionless).  All threshold values in this block
    # are heuristic — pending calibration.
    # Warning tier: worst single bridge residual before a warning is raised.
    warning_max_bridge_distance_residual: float = 0.75
    # Warning tier: mean bridge residual across all events.
    warning_mean_bridge_distance_residual: float = 0.35
    # Per-event residual above which a bridge counts as "bad" for the
    # bad-fraction metric gated by warning_max_bad_bridge_fraction.
    warning_bad_bridge_distance_residual: float = 0.50
    # Warning tier: fraction of bad bridge events before a warning is raised.
    warning_max_bad_bridge_fraction: float = 0.25
    # Hard-hard tier (blocks CIF export): largest actual bridge distance in
    # angstrom; far above genuine linkage bond lengths (~1.3-1.6 A for the
    # supported templates), so only disconnected or mis-assembled pairs trip it.
    hard_hard_max_bridge_distance: float = 2.5
    # Hard tier: worst single bridge residual.
    hard_max_bridge_distance_residual: float = 1.00
    # Hard tier: mean bridge residual across all events.
    hard_mean_bridge_distance_residual: float = 0.60
    # Hard tier: actual/target ratio below which a bridge is too short.
    hard_min_bridge_distance_ratio: float = 0.70
    # Hard tier: actual/target ratio above which a bridge is too long.
    hard_max_bridge_distance_ratio: float = 1.60
    # Coarse plain-distance backstop: hard clash floor for non-excluded
    # heavy-heavy pairs, and the severe-overlap floor for graph-excluded 1-3
    # pairs (a heavy 1-3 contact below bond-length scale means a broken graph).
    # Heuristic — pending calibration: 1.05 A sits just below the shortest
    # genuine heavy-heavy bond lengths (~1.16 A for C#N), so only fused or
    # collapsed geometry can trip the floor.
    min_nonbonded_heavy_distance: float = 1.05
    # W4.1: primary clash criterion; flag a non-excluded pair when
    # d / (r_vdw_i + r_vdw_j) < this value (Bondi radii from cofkit.vdw).
    min_nonbonded_heavy_vdw_ratio: float = DEFAULT_VDW_CLASH_RATIO
    # W4.1: search radius for the contact scan; covers the largest flaggable
    # vdW sum in the supported table (0.75 * 2 * 2.10 = 3.15 A for Si...Si)
    # plus margin.
    nonbonded_heavy_search_radius: float = DEFAULT_NONBONDED_SEARCH_RADIUS
    # Severe-overlap floor (plain distance) applied to non-excluded and to
    # graph-excluded heavy-heavy 1-4 pairs; flags independent of the radius
    # table.  Cis/gauche 1-4 contacts can legitimately sit below
    # 0.75 * sum(r_vdw) (~2.55 A for C...C) but never below this floor.
    hard_min_nonbonded_heavy_plain_distance: float = DEFAULT_HARD_MIN_NONBONDED_PLAIN_DISTANCE
    # Severe-overlap floor for directly bonded (1-2) pairs: below this the
    # "bond" is fused nuclei, i.e. a broken structure.  Stays below genuine
    # heavy-atom bonds (~1.16 A for C#N) and X-H bonds (~1.0 A).
    severe_overlap_bonded_distance: float = DEFAULT_SEVERE_OVERLAP_BONDED_DISTANCE
    # Degenerate-cell floors (heuristic — pending calibration): a 2D cell
    # under 10 A^2 or a 3D cell under 20 A^3 cannot contain even one monomer;
    # real candidate cells are orders of magnitude larger (the placeholder
    # lateral span alone is 30 A, i.e. 900 A^2).
    min_2d_cell_area: float = 10.0
    min_3d_cell_volume: float = 20.0
    # Per-template acceptable distance windows for the realized inter-monomer
    # linkage bonds measured from the exported CIF coordinates and periodic
    # bond images (measure_inter_instance_bond_distances), keyed by template
    # id. The boronate ester window brackets the ~1.47 angstrom B-O bonds of
    # the five-membered ring; templates whose profile bridge_target_distance
    # IS the realized bond length are validated by residual/ratio thresholds
    # against that target instead.
    realized_bridge_bond_distance_windows: Mapping[str, tuple[float, float]] = field(
        default_factory=lambda: {"boronate_ester_bridge": (1.25, 1.65)}
    )
    # When metadata-level hard reasons exist, skip only the metadata-dependent
    # CIF checks (template bond windows); the geometry checks (clash scan,
    # cell degeneracy, instance graph) always run.
    skip_cif_checks_when_metadata_invalid: bool = True


@dataclass(frozen=True)
class CoarseValidationReport:
    classification: str
    is_valid: bool | None
    passes_hard_validation: bool | None
    blocks_cif_export: bool = False
    reasons: tuple[str, ...] = ()
    hard_hard_invalid_reasons: tuple[str, ...] = ()
    warning_reasons: tuple[str, ...] = ()
    hard_invalid_reasons: tuple[str, ...] = ()
    needs_optimization_reasons: tuple[str, ...] = ()
    unmeasured_required_checks: tuple[str, ...] = ()
    coverage: Mapping[str, str] = field(default_factory=dict)
    metrics: Mapping[str, object] = field(default_factory=dict)


# Coverage statuses for the per-check map on CoarseValidationReport.coverage.
# They deliberately distinguish "the check ran and produced values" from "the
# check ran to completion and found nothing" from "the check could not run".
CHECK_MEASURED = "measured"
# A completed bounded search/measurement that found zero items (e.g. a contact
# scan with no non-excluded neighbor inside the search radius). Distinct from
# missing data: the measurement was performed.
CHECK_NO_CONTACTS = "no_contacts"
# A required check could not be evaluated because its inputs are absent
# (missing CIF, missing bond loop, no linkage target/window for the template).
CHECK_MISSING_DATA = "missing_data"
# The check does not apply to this record (e.g. bridge geometry when the
# assembly metadata reports no bridge events and no inter-monomer bonds).
CHECK_NOT_APPLICABLE = "not_applicable"
# The check was deliberately skipped because metadata-level hard reasons
# already determine the verdict (see
# CoarseValidationThresholds.skip_cif_checks_when_metadata_invalid).
CHECK_SKIPPED = "skipped"

# Checks whose absence blocks a "valid" verdict: without them the export is
# unvalidated, not valid. ring_geometry (ring-participant arrangement) and
# ring_attachment (exocyclic ring-monomer bonds/angles) are covered by their
# own metadata channel and are not part of the required set; a missing
# attachment channel on a ring record is reported as missing_data coverage
# but does not by itself make the record unvalidated.
REQUIRED_COVERAGE_CHECKS = ("cell_geometry", "instance_graph", "contact_scan", "bridge_geometry")


@dataclass(frozen=True)
class RealizedBridgeBondMeasurement:
    """One inter-monomer linkage bond measured from exported CIF coordinates.

    The distance is recomputed from fractional coordinates and the explicit
    periodic bond-image codes in the ``_geom_bond`` loop, never trusted from
    the ``_geom_bond_distance`` column, so repaired/relaxed coordinates are
    measured as exported.
    """

    label_1: str
    label_2: str
    image: tuple[int, int, int]
    distance: float


def _instance_id(atom_label: str) -> str:
    # Atom labels are written as "{instance_id}_{symbol}{n}" and instance
    # ids themselves may contain underscores (e.g. "bex_node"), so split
    # off only the trailing symbol tag. Matches decompose.py:_instance_id.
    return str(atom_label).rsplit("_", 1)[0]


def measure_inter_instance_bond_distances(block) -> tuple[RealizedBridgeBondMeasurement, ...]:
    """Measure inter-monomer (bridge) bond lengths from a CIF block.

    Bonds are taken from the ``_geom_bond`` loop; a bond is a bridge bond when
    its two atom labels belong to different monomer instances. Distances are
    recomputed from the fractional coordinates and the P1 image codes in
    ``_geom_bond_site_symmetry_1/2`` (via ``periodic_geometry.p1_shift``), so
    bonds crossing a cell boundary measure correctly. Rows referencing unknown
    atom labels are skipped.
    """
    small = gemmi.make_small_structure_from_block(block)
    positions = {
        str(site.label): (site.fract.x, site.fract.y, site.fract.z) for site in small.sites
    }
    labels1 = block.find_loop("_geom_bond_atom_site_label_1")
    labels2 = block.find_loop("_geom_bond_atom_site_label_2")
    sym1 = block.find_loop("_geom_bond_site_symmetry_1")
    sym2 = block.find_loop("_geom_bond_site_symmetry_2")
    measurements: list[RealizedBridgeBondMeasurement] = []
    for index in range(min(len(labels1), len(labels2))):
        label_a = cif_value_str(labels1[index])
        label_b = cif_value_str(labels2[index])
        if _instance_id(label_a) == _instance_id(label_b):
            continue
        fract_a = positions.get(label_a)
        fract_b = positions.get(label_b)
        if fract_a is None or fract_b is None:
            continue
        image_a = p1_shift(cif_value_str(sym1[index]) if len(sym1) else ".")
        image_b = p1_shift(cif_value_str(sym2[index]) if len(sym2) else ".")
        shift = tuple(b - a for a, b in zip(image_a, image_b))
        vector = gemmi.Fractional(
            fract_b[0] + shift[0] - fract_a[0],
            fract_b[1] + shift[1] - fract_a[1],
            fract_b[2] + shift[2] - fract_a[2],
        )
        measurements.append(
            RealizedBridgeBondMeasurement(
                label_1=label_a,
                label_2=label_b,
                image=shift,
                distance=float(small.cell.orthogonalize(vector).length()),
            )
        )
    return tuple(measurements)


@dataclass(frozen=True)
class BatchOutputClassificationSummary:
    source_dir: str
    output_dir: str
    total_structures: int
    valid_structures: int
    warning_structures: int
    needs_optimization_structures: int = 0
    hard_hard_invalid_structures: int = 0
    hard_invalid_structures: int = 0
    warning_reason_counts: Mapping[str, int] = field(default_factory=dict)
    needs_optimization_reason_counts: Mapping[str, int] = field(default_factory=dict)
    hard_hard_invalid_reason_counts: Mapping[str, int] = field(default_factory=dict)
    hard_invalid_reason_counts: Mapping[str, int] = field(default_factory=dict)
    unvalidated_structures: int = 0
    unvalidated_reason_counts: Mapping[str, int] = field(default_factory=dict)
    classification_manifest_path: str = ""
    valid_manifest_path: str = ""
    warning_manifest_path: str = ""
    needs_optimization_manifest_path: str = ""
    hard_hard_invalid_manifest_path: str = ""
    hard_invalid_manifest_path: str = ""
    unvalidated_manifest_path: str = ""

    @property
    def invalid_structures(self) -> int:
        return self.hard_hard_invalid_structures + self.hard_invalid_structures


class CoarseStructureValidator:
    _REPAIRABLE_GEOMETRY_REASONS = frozenset(
        {
            "bridge_distance_residual_max_hard",
            "bridge_distance_residual_mean_hard",
            "bridge_distance_too_short",
            "bridge_distance_too_long",
            "realized_bridge_bond_distance",
        }
    )

    def __init__(self, thresholds: CoarseValidationThresholds | None = None) -> None:
        self.thresholds = thresholds or CoarseValidationThresholds()

    def validate_manifest_record(
        self,
        record: Mapping[str, object],
        *,
        source_root: str | Path | None = None,
    ) -> CoarseValidationReport:
        warning_reasons: list[str] = []
        hard_hard_invalid_reasons: list[str] = []
        hard_invalid_reasons: list[str] = []
        metrics: dict[str, object] = {}
        coverage: dict[str, str] = {}

        metadata = self._mapping(record.get("metadata"))
        graph_summary = self._mapping(metadata.get("graph_summary"))
        score_metadata = self._mapping(metadata.get("score_metadata"))
        ring_validation = self._mapping(metadata.get("ring_validation"))
        stacking_metadata = self._mapping(metadata.get("stacking"))
        n_unreacted_motifs = self._n_unreacted_motifs(record, score_metadata)
        metrics["n_unreacted_motifs"] = n_unreacted_motifs
        if n_unreacted_motifs > 0:
            hard_invalid_reasons.append("unreacted_motifs")

        if ring_validation:
            # Arrangement and attachment verdicts are separate channels:
            # arrangement covers the ring participants themselves, attachment
            # the exocyclic ring-monomer bonds/angles. Legacy records carry
            # only the combined "classification"; treat it as the arrangement
            # verdict and report attachment coverage as missing data.
            ring_classification = self._string(ring_validation.get("classification"))
            arrangement_classification = self._string(ring_validation.get("arrangement_classification"))
            if arrangement_classification is None:
                arrangement_classification = ring_classification
            attachment_classification = self._string(ring_validation.get("attachment_classification"))
            attachment_status = self._string(ring_validation.get("attachment_status"))
            ring_reasons = tuple(str(reason) for reason in ring_validation.get("reasons", ()))
            metrics["ring_validation_classification"] = ring_classification
            metrics["ring_validation_reasons"] = ring_reasons
            metrics["ring_geometry"] = ring_validation.get("metrics", {})
            coverage["ring_geometry"] = CHECK_MEASURED
            if arrangement_classification not in {None, "accepted", "valid"}:
                hard_invalid_reasons.append("ring_geometry_invalid")
            attachment_block = self._mapping(ring_validation.get("attachment"))
            if attachment_classification is None:
                coverage["ring_attachment"] = CHECK_MISSING_DATA
                metrics["ring_attachment_coverage_detail"] = (
                    "ring_validation metadata carries no attachment channel (record predates "
                    "attachment measurement or the attachment report was dropped)"
                )
            elif attachment_status == CHECK_MISSING_DATA or attachment_block.get("status") == CHECK_MISSING_DATA:
                coverage["ring_attachment"] = CHECK_MISSING_DATA
                metrics["ring_attachment_classification"] = attachment_classification
                metrics["ring_attachment_coverage_detail"] = (
                    f"{attachment_block.get('n_unmeasured_participants', '?')} ring participant(s) "
                    "lacked atom metadata/positions for attachment measurement"
                )
                metrics["ring_attachment"] = attachment_block
            elif attachment_classification == "not_applicable":
                coverage["ring_attachment"] = CHECK_NOT_APPLICABLE
                metrics["ring_attachment_classification"] = attachment_classification
            else:
                coverage["ring_attachment"] = CHECK_MEASURED
                metrics["ring_attachment_classification"] = attachment_classification
                metrics["ring_attachment"] = attachment_block
                if attachment_classification == "rejected":
                    hard_invalid_reasons.append("ring_attachment_invalid")
                elif attachment_classification == "warning":
                    warning_reasons.append("ring_attachment_strained")
        else:
            coverage["ring_geometry"] = CHECK_NOT_APPLICABLE
            coverage["ring_attachment"] = CHECK_NOT_APPLICABLE

        monomer_geometry_warnings = tuple(
            str(warning) for warning in metadata.get("monomer_geometry_warnings", ()) or ()
        )
        if monomer_geometry_warnings:
            warning_reasons.append("monomer_geometry_degraded")
            metrics["monomer_geometry_degraded_details"] = monomer_geometry_warnings

        # Assembly-time ("seed") bridge metrics are informational only: they
        # describe the embedding/assembly coordinates, not the exported
        # structure. Validation verdicts on bridge geometry are derived from
        # the final CIF coordinates measured below, so a repair pass that
        # moves atoms cannot inherit a stale passing (or failing) verdict.
        bridge_metrics = tuple(self._mapping(item) for item in score_metadata.get("bridge_event_metrics", ()))
        raw_distance_residuals = tuple(self._distance_residual(item) for item in bridge_metrics)
        seed_distance_residuals = tuple(value for value in raw_distance_residuals if value is not None)
        n_missing_distance_data = len(raw_distance_residuals) - len(seed_distance_residuals)
        metrics["n_bridge_events"] = len(bridge_metrics)
        if n_missing_distance_data:
            metrics["n_bridge_events_missing_distance_data"] = n_missing_distance_data
        if seed_distance_residuals:
            metrics["seed_max_bridge_distance_residual"] = max(seed_distance_residuals)
            metrics["seed_mean_bridge_distance_residual"] = (
                sum(seed_distance_residuals) / len(seed_distance_residuals)
            )
        seed_targets = tuple(
            sorted(
                {
                    float(item["target_distance"])
                    for item in bridge_metrics
                    if item.get("target_distance") is not None
                }
            )
        )
        if seed_targets:
            metrics["seed_bridge_target_distances"] = seed_targets
        metrics["bridge_metrics_source"] = None

        template_id = self._record_template_id(metadata, graph_summary)

        cif_path = self._resolve_cif_path(record, source_root=source_root)
        metrics["cif_path"] = str(cif_path) if cif_path is not None else None
        if cif_path is None or not cif_path.is_file():
            hard_invalid_reasons.append("cif_missing")
            for check in REQUIRED_COVERAGE_CHECKS:
                coverage.setdefault(check, CHECK_MISSING_DATA)
        else:
            skip_metadata_checks = bool(
                self.thresholds.skip_cif_checks_when_metadata_invalid and hard_invalid_reasons
            )
            cif_metrics, cif_reasons, cif_warnings = self._validate_cif(
                cif_path,
                topology_id=self._string(record.get("topology_id")),
                stacking_metadata=stacking_metadata,
                template_id=template_id,
                skip_metadata_checks=skip_metadata_checks,
            )
            metrics.update(cif_metrics)
            hard_invalid_reasons.extend(cif_reasons)
            warning_reasons.extend(cif_warnings)
            cif_parsed = "cif_parse_error" not in cif_metrics
            if cif_parsed:
                coverage["cell_geometry"] = CHECK_MEASURED
                coverage["instance_graph"] = CHECK_MEASURED
                coverage["contact_scan"] = (
                    CHECK_MEASURED
                    if int(cif_metrics.get("n_nonbonded_pairs_scanned", 0)) > 0
                    else CHECK_NO_CONTACTS
                )
                coverage["bridge_geometry"] = self._apply_final_bridge_geometry_verdicts(
                    metrics,
                    template_id=template_id,
                    seed_targets=seed_targets,
                    n_bridge_events=len(bridge_metrics),
                    window_check_skipped=bool(cif_metrics.get("cif_metadata_checks_skipped")),
                    warning_reasons=warning_reasons,
                    hard_invalid_reasons=hard_invalid_reasons,
                    hard_hard_invalid_reasons=hard_hard_invalid_reasons,
                )
            else:
                for check in REQUIRED_COVERAGE_CHECKS:
                    coverage.setdefault(check, CHECK_MISSING_DATA)

        normalized_warning_reasons = tuple(dict.fromkeys(warning_reasons))
        normalized_hard_hard_invalid_reasons = tuple(dict.fromkeys(hard_hard_invalid_reasons))
        normalized_hard_invalid_reasons = tuple(dict.fromkeys(hard_invalid_reasons))
        needs_optimization_reasons: tuple[str, ...] = ()
        remaining_hard_invalid_reasons = normalized_hard_invalid_reasons
        if (
            not normalized_hard_hard_invalid_reasons
            and normalized_hard_invalid_reasons
            and all(reason in self._REPAIRABLE_GEOMETRY_REASONS for reason in normalized_hard_invalid_reasons)
        ):
            needs_optimization_reasons = normalized_hard_invalid_reasons
            remaining_hard_invalid_reasons = ()

        unmeasured_required_checks = tuple(
            check for check in REQUIRED_COVERAGE_CHECKS if coverage.get(check) == CHECK_MISSING_DATA
        )
        if unmeasured_required_checks:
            metrics["unmeasured_required_checks"] = unmeasured_required_checks

        normalized_reasons = normalized_hard_hard_invalid_reasons + remaining_hard_invalid_reasons + tuple(
            reason for reason in needs_optimization_reasons if reason not in remaining_hard_invalid_reasons
        ) + (("unmeasured_required_checks",) if unmeasured_required_checks else ()) + tuple(
            reason for reason in normalized_warning_reasons if reason not in normalized_hard_invalid_reasons
        )
        classification = "valid"
        if normalized_hard_hard_invalid_reasons:
            classification = "hard_hard_invalid"
        elif remaining_hard_invalid_reasons:
            classification = "hard_invalid"
        elif needs_optimization_reasons:
            classification = "needs_optimization"
        elif unmeasured_required_checks:
            # Required checks could not be evaluated: the export is honestly
            # unvalidated, never silently valid.
            classification = "unvalidated"
        elif normalized_warning_reasons:
            classification = "warning"
        return CoarseValidationReport(
            classification=classification,
            is_valid=True if classification == "valid" else (None if classification == "unvalidated" else False),
            passes_hard_validation=(
                classification in {"valid", "warning"} if classification != "unvalidated" else None
            ),
            blocks_cif_export=bool(normalized_hard_hard_invalid_reasons),
            reasons=normalized_reasons,
            hard_hard_invalid_reasons=normalized_hard_hard_invalid_reasons,
            warning_reasons=normalized_warning_reasons,
            hard_invalid_reasons=remaining_hard_invalid_reasons,
            needs_optimization_reasons=needs_optimization_reasons,
            unmeasured_required_checks=unmeasured_required_checks,
            coverage=coverage,
            metrics=metrics,
        )

    def _record_template_id(
        self,
        metadata: Mapping[str, object],
        graph_summary: Mapping[str, object],
    ) -> str | None:
        template_id = self._string(metadata.get("template_id"))
        if template_id is not None:
            return template_id
        reaction_templates = graph_summary.get("reaction_templates")
        if isinstance(reaction_templates, Mapping) and len(reaction_templates) == 1:
            return self._string(next(iter(reaction_templates)))
        return None

    def _apply_final_bridge_geometry_verdicts(
        self,
        metrics: dict[str, object],
        *,
        template_id: str | None,
        seed_targets: tuple[float, ...],
        n_bridge_events: int,
        window_check_skipped: bool,
        warning_reasons: list[str],
        hard_invalid_reasons: list[str],
        hard_hard_invalid_reasons: list[str],
    ) -> str:
        """Derive bridge-geometry verdicts from the final CIF coordinates.

        Returns the coverage status for the ``bridge_geometry`` check. The
        measurements in ``metrics["realized_bridge_bond_measurements"]`` are
        recomputed from the exported coordinates and periodic bond images (see
        ``measure_inter_instance_bond_distances``), so verdicts reflect the
        structure as written, including after any repair pass moved atoms.
        """
        measurements = tuple(
            self._mapping(item) for item in metrics.get("realized_bridge_bond_measurements", ())
        )
        distances = tuple(
            float(item["distance"]) for item in measurements if item.get("distance") is not None
        )
        thresholds = self.thresholds
        if distances:
            # Template-independent backstop: a realized inter-monomer "bond"
            # beyond the export limit is disconnected/mis-assembled geometry
            # no matter which linkage family produced it.
            max_actual_distance = max(distances)
            metrics["max_actual_bridge_distance"] = max_actual_distance
            if max_actual_distance >= thresholds.hard_hard_max_bridge_distance:
                hard_hard_invalid_reasons.append("bridge_distance_exceeds_cif_export_limit")
        bond_window = thresholds.realized_bridge_bond_distance_windows.get(template_id or "")
        if bond_window is not None and window_check_skipped:
            return CHECK_SKIPPED
        if bond_window is not None:
            # Window verdict (realized_bridge_bond_distance) is applied by
            # _validate_cif against the same measured distances.
            if distances:
                metrics["bridge_metrics_source"] = "final_cif"
                return CHECK_MEASURED
            if n_bridge_events:
                metrics["bridge_geometry_coverage_detail"] = "no inter-monomer bonds in the CIF bond loop"
                return CHECK_MISSING_DATA
            return CHECK_NOT_APPLICABLE

        profile = linkage_profile(template_id) if template_id is not None else None
        target_distance: float | None = None
        target_source: str | None = None
        if profile is not None:
            target_distance = float(profile.bridge_target_distance)
            target_source = "linkage_profile"
        elif len(seed_targets) == 1:
            # Fallback for records without resolvable template metadata: the
            # assembly-time target is a template property, not a coordinate
            # measurement, so it stays valid after repair.
            target_distance = float(seed_targets[0])
            target_source = "seed_metadata"
        if target_source is not None:
            metrics["bridge_target_source"] = target_source

        if not distances:
            if n_bridge_events:
                metrics["bridge_geometry_coverage_detail"] = "no inter-monomer bonds in the CIF bond loop"
                return CHECK_MISSING_DATA
            return CHECK_NOT_APPLICABLE
        if target_distance is None or target_distance <= 0.0:
            metrics["bridge_geometry_coverage_detail"] = (
                f"no realized bridge-bond target or distance window for template {template_id!r}"
            )
            return CHECK_MISSING_DATA

        residuals = tuple(abs(distance - target_distance) for distance in distances)
        ratios = tuple(distance / target_distance for distance in distances)
        max_residual = max(residuals)
        mean_residual = sum(residuals) / len(residuals)
        bad_fraction = (
            sum(value > thresholds.warning_bad_bridge_distance_residual for value in residuals)
            / len(residuals)
        )
        metrics["bridge_metrics_source"] = "final_cif"
        metrics["n_measured_bridge_bonds"] = len(distances)
        metrics["bridge_target_distance"] = target_distance
        metrics["max_bridge_distance_residual"] = max_residual
        metrics["mean_bridge_distance_residual"] = mean_residual
        metrics["bad_bridge_event_fraction"] = bad_fraction
        metrics["min_bridge_distance_ratio"] = min(ratios)
        metrics["max_bridge_distance_ratio"] = max(ratios)
        if max_residual > thresholds.warning_max_bridge_distance_residual:
            warning_reasons.append("bridge_distance_residual_max")
        if mean_residual > thresholds.warning_mean_bridge_distance_residual:
            warning_reasons.append("bridge_distance_residual_mean")
        if bad_fraction > thresholds.warning_max_bad_bridge_fraction:
            warning_reasons.append("bridge_distance_residual_fraction")
        if max_residual > thresholds.hard_max_bridge_distance_residual:
            hard_invalid_reasons.append("bridge_distance_residual_max_hard")
        if mean_residual > thresholds.hard_mean_bridge_distance_residual:
            hard_invalid_reasons.append("bridge_distance_residual_mean_hard")
        if min(ratios) < thresholds.hard_min_bridge_distance_ratio:
            hard_invalid_reasons.append("bridge_distance_too_short")
        if max(ratios) > thresholds.hard_max_bridge_distance_ratio:
            hard_invalid_reasons.append("bridge_distance_too_long")
        return CHECK_MEASURED

    def _validate_cif(
        self,
        cif_path: Path,
        *,
        topology_id: str | None,
        stacking_metadata: Mapping[str, object] | None = None,
        template_id: str | None = None,
        skip_metadata_checks: bool = False,
    ) -> tuple[dict[str, object], tuple[str, ...], tuple[str, ...]]:
        if gemmi is None:  # pragma: no cover - import guard for incomplete environments
            raise ModuleNotFoundError(
                "gemmi is required for CIF-backed coarse validation. Install gemmi into this environment."
            )

        metrics: dict[str, object] = {}
        reasons: list[str] = []
        warnings: list[str] = []
        try:
            block = gemmi.cif.read_file(str(cif_path)).sole_block()
            small = gemmi.make_small_structure_from_block(block)
        except Exception as exc:  # pragma: no cover - defensive
            metrics["cif_parse_error"] = f"{type(exc).__name__}: {exc}"
            return metrics, ("cif_parse_error",), ()

        cell = small.cell
        dimensionality = self._topology_dimensionality(topology_id)
        cell_area_2d = cell.a * cell.b * abs(sin(radians(cell.gamma)))
        cell_volume_3d = float(cell.volume)
        metrics["cell_area_2d"] = cell_area_2d
        metrics["cell_volume_3d"] = cell_volume_3d
        metrics["topology_dimensionality"] = dimensionality
        if dimensionality == "2D":
            if cell_area_2d < self.thresholds.min_2d_cell_area:
                reasons.append("degenerate_cell")
        elif cell_volume_3d < self.thresholds.min_3d_cell_volume:
            reasons.append("degenerate_cell")

        bonded_pairs = self._bonded_pairs(block)
        instance_graph = self._instance_graph(small, bonded_pairs)
        metrics.update(instance_graph["metrics"])
        if int(instance_graph["metrics"]["n_instance_components"]) > 1:
            if self._stacked_multilayer_components_are_allowed(
                stacking_metadata,
                instance_graph["metrics"],
            ):
                metrics["disconnected_instance_graph_allowed"] = True
            else:
                reasons.append("disconnected_instance_graph")

        contact_scan = self._nonbonded_contact_scan(small, *self._bond_graph_exclusions(block))
        metrics.update(contact_scan)
        if int(contact_scan["n_heavy_atom_clash_pairs"]) > 0:
            reasons.append("heavy_atom_clash")
        if int(contact_scan["n_excluded_pair_severe_overlaps"]) > 0:
            reasons.append("excluded_pair_severe_overlap")
        if int(contact_scan["n_hydrogen_atom_clash_pairs"]) > 0:
            warnings.append("hydrogen_atom_clash")

        # Final-geometry bridge measurement: inter-monomer bond distances
        # recomputed from the exported coordinates and periodic bond images.
        # Verdicts are derived by validate_manifest_record (residual vs the
        # linkage target) except the per-template distance-window check, which
        # stays here because it is metadata-dependent and can be skipped.
        measurements = measure_inter_instance_bond_distances(block)
        metrics["realized_bridge_bond_count"] = len(measurements)
        metrics["realized_bridge_bond_measurements"] = tuple(
            {
                "label_1": measurement.label_1,
                "label_2": measurement.label_2,
                "image": tuple(measurement.image),
                "distance": measurement.distance,
            }
            for measurement in measurements
        )
        if skip_metadata_checks:
            bond_window = None
            metrics["cif_metadata_checks_skipped"] = True
        else:
            bond_window = self.thresholds.realized_bridge_bond_distance_windows.get(template_id or "")
        if bond_window is not None and measurements:
            realized_distances = tuple(measurement.distance for measurement in measurements)
            min_realized = min(realized_distances)
            max_realized = max(realized_distances)
            metrics["realized_bridge_bond_distance_min"] = min_realized
            metrics["realized_bridge_bond_distance_max"] = max_realized
            if min_realized < bond_window[0] or max_realized > bond_window[1]:
                reasons.append("realized_bridge_bond_distance")

        return metrics, tuple(dict.fromkeys(reasons)), tuple(dict.fromkeys(warnings))

    def _bonded_pairs(self, block) -> set[frozenset[str]]:
        label_1 = block.find_loop("_geom_bond_atom_site_label_1")
        label_2 = block.find_loop("_geom_bond_atom_site_label_2")
        if len(label_1) == 0 or len(label_2) == 0:
            return set()
        return {
            frozenset((cif_value_str(label_1[index]), cif_value_str(label_2[index])))
            for index in range(min(len(label_1), len(label_2)))
        }

    def _instance_graph(
        self,
        small,
        bonded_pairs: set[frozenset[str]],
    ) -> dict[str, object]:
        instance_ids = {self._instance_id(site.label) for site in small.sites}
        adjacency: dict[str, set[str]] = defaultdict(set)
        inter_instance_edges = 0
        for pair in bonded_pairs:
            if len(pair) != 2:
                continue
            label_1, label_2 = tuple(pair)
            instance_1 = self._instance_id(label_1)
            instance_2 = self._instance_id(label_2)
            if instance_1 == instance_2:
                continue
            if instance_2 not in adjacency[instance_1]:
                inter_instance_edges += 1
            adjacency[instance_1].add(instance_2)
            adjacency[instance_2].add(instance_1)

        components = 0
        component_sizes: list[int] = []
        seen: set[str] = set()
        for instance_id in instance_ids:
            if instance_id in seen:
                continue
            components += 1
            stack = [instance_id]
            seen.add(instance_id)
            component_size = 0
            while stack:
                current = stack.pop()
                component_size += 1
                for other in adjacency.get(current, ()):
                    if other not in seen:
                        seen.add(other)
                        stack.append(other)
            component_sizes.append(component_size)
        return {
            "metrics": {
                "n_instance_nodes": len(instance_ids),
                "n_instance_components": components if instance_ids else 0,
                "n_inter_instance_edges": inter_instance_edges,
                "instance_component_sizes": tuple(sorted(component_sizes)),
            }
        }

    def _stacked_multilayer_components_are_allowed(
        self,
        stacking_metadata: Mapping[str, object] | None,
        instance_graph_metrics: Mapping[str, object],
    ) -> bool:
        metadata = self._mapping(stacking_metadata)
        try:
            layer_count = int(metadata.get("layer_count"))
        except (TypeError, ValueError):
            return False
        if layer_count <= 1:
            return False
        try:
            n_components = int(instance_graph_metrics.get("n_instance_components"))
        except (TypeError, ValueError):
            return False
        if n_components != layer_count:
            return False
        raw_sizes = instance_graph_metrics.get("instance_component_sizes", ())
        if not isinstance(raw_sizes, (list, tuple)):
            return False
        try:
            component_sizes = tuple(int(size) for size in raw_sizes)
        except (TypeError, ValueError):
            return False
        if len(component_sizes) != layer_count:
            return False
        return len(set(component_sizes)) == 1

    def _bond_graph_exclusions(
        self, block
    ) -> tuple[
        set[tuple[str, str, tuple[int, int, int]]],
        set[tuple[str, str, tuple[int, int, int]]],
        set[tuple[str, str, tuple[int, int, int]]],
    ]:
        """Periodic-image-aware 1-2 / 1-3 / 1-4 exclusion sets from the bond loop.

        Each set holds directed ``(label, other_label, relative_image)`` triples
        meaning "``other_label`` at ``relative_image`` is graph-related to
        ``label`` at the home image": 1-2 from the ``_geom_bond`` loop directly,
        1-3 from pairs of bonded neighbors around a common center (relative
        image ``tb - ta``), and 1-4 from bond paths ``a--i--j--b`` (relative
        image ``s_ij + tb - ta``).  The image bookkeeping mirrors
        ``soft_relax._parse_system`` and composes correctly for bonds that
        cross a cell boundary.
        """
        labels1 = block.find_loop("_geom_bond_atom_site_label_1")
        labels2 = block.find_loop("_geom_bond_atom_site_label_2")
        sym1 = block.find_loop("_geom_bond_site_symmetry_1")
        sym2 = block.find_loop("_geom_bond_site_symmetry_2")
        bonded: set[tuple[str, str, tuple[int, int, int]]] = set()
        # Directed adjacency: adjacency[x] holds (y, t) meaning "y at image t
        # is bonded to x at the home image".
        adjacency: dict[str, set[tuple[str, tuple[int, int, int]]]] = defaultdict(set)
        for i in range(min(len(labels1), len(labels2))):
            label_a = cif_value_str(labels1[i])
            label_b = cif_value_str(labels2[i])
            first = p1_shift(cif_value_str(sym1[i]) if len(sym1) else ".")
            second = p1_shift(cif_value_str(sym2[i]) if len(sym2) else ".")
            shift = tuple(b - a for a, b in zip(first, second))
            bonded.add((label_a, label_b, shift))
            bonded.add((label_b, label_a, tuple(-v for v in shift)))
            adjacency[label_a].add((label_b, shift))
            adjacency[label_b].add((label_a, tuple(-v for v in shift)))

        excluded13: set[tuple[str, str, tuple[int, int, int]]] = set()
        for neighbors in adjacency.values():
            ordered = sorted(neighbors)
            for x, (label_a, ta) in enumerate(ordered):
                for label_b, tb in ordered[x + 1 :]:
                    if label_a == label_b and ta == tb:
                        continue
                    # Relative image of b seen from a is tb - ta.
                    relative = tuple(vb - va for va, vb in zip(ta, tb))
                    excluded13.add((label_a, label_b, relative))
                    excluded13.add((label_b, label_a, tuple(-v for v in relative)))

        excluded14: set[tuple[str, str, tuple[int, int, int]]] = set()
        for label_i, neighbors_i in adjacency.items():
            for label_j, s_ij in neighbors_i:
                for label_a, ta in adjacency.get(label_i, ()):
                    if label_a == label_j:
                        continue
                    for label_b, tb in adjacency.get(label_j, ()):
                        if label_b == label_i:
                            continue
                        # Relative image of b seen from a is (s_ij + tb) - ta.
                        relative = tuple(s + vtb - vta for s, vtb, vta in zip(s_ij, tb, ta))
                        excluded14.add((label_a, label_b, relative))
                        excluded14.add((label_b, label_a, tuple(-v for v in relative)))
        return bonded, excluded13, excluded14

    def _nonbonded_contact_scan(
        self,
        small,
        bonded: set[tuple[str, str, tuple[int, int, int]]],
        excluded13: set[tuple[str, str, tuple[int, int, int]]],
        excluded14: set[tuple[str, str, tuple[int, int, int]]],
    ) -> dict[str, object]:
        """vdW-aware periodic contact scan with bond-graph exclusions (W4.1).

        Criterion: a pair that is not graph-excluded is a clash when
        ``d / (r_vdw_i + r_vdw_j) < min_nonbonded_heavy_vdw_ratio`` (Bondi
        radii from ``cofkit.vdw``) or — heavy-heavy pairs only — when the
        plain distance falls below the severe-overlap floor
        ``hard_min_nonbonded_heavy_plain_distance`` (radius-table-independent
        backstop).  ``nonbonded_heavy_search_radius`` is sized to cover the
        largest flaggable vdW sum for supported elements.

        Bond-graph exclusion policy (periodic-image-aware label/image triples,
        built by ``_bond_graph_exclusions``):

        * 1-2 (bonded) pairs are excluded from the ratio check; below
          ``severe_overlap_bonded_distance`` they are reported as
          ``excluded_pair_severe_overlap`` — a "bonded" pair at fused-nuclei
          distance means a broken structure, and the exclusion must not hide it.
        * 1-3 (angle) pairs are excluded from the ratio check: ordinary angle
          neighbors (benzene 1-3 C...C at 2.42 A, H-C-H at ~1.9 A) sit below
          the naive 0.75-ratio cutoff.  Below ``min_nonbonded_heavy_distance``
          a heavy 1-3 pair is a severe overlap.
        * 1-4 (torsion) pairs are excluded from the ratio check: cis/gauche
          torsion contacts legitimately sit below 0.75 * sum(r_vdw).  They
          remain subject to the severe-overlap floor —
          ``hard_min_nonbonded_heavy_plain_distance`` for heavy-heavy 1-4
          pairs, ``min_nonbonded_heavy_distance`` for hydrogen-involving ones
          (eclipsed X-H...H-C torsion contacts can sit near 2.3 A).

        Hydrogen-involving contacts are reported in a separate metric channel
        (``min_nonbonded_hydrogen_*``) and flagged as the warning reason
        ``hydrogen_atom_clash`` rather than ``heavy_atom_clash``.

        ``n_nonbonded_pairs_scanned`` counts the non-excluded pairs that
        reached the assessment criterion; zero means the bounded search
        completed without finding any assessable neighbor (coverage status
        ``no_contacts``), which is distinct from the scan not running at all.
        """
        thresholds = self.thresholds
        search_radius = thresholds.nonbonded_heavy_search_radius
        search = gemmi.NeighborSearch(small, search_radius).populate(include_h=True)
        pairs_scanned = 0
        min_heavy_distance: float | None = None
        min_heavy_ratio: float | None = None
        min_hydrogen_distance: float | None = None
        min_hydrogen_ratio: float | None = None
        heavy_clash_details: list[dict[str, object]] = []
        hydrogen_clash_details: list[dict[str, object]] = []
        severe_overlap_details: list[dict[str, object]] = []
        for index, site in enumerate(small.sites):
            candidates = {(index, 0)}
            candidates.update((int(mark.atom_idx), int(mark.image_idx))
                for mark in search.find_site_neighbors(site, min_dist=0, max_dist=search_radius))
            for other_index, image_index in candidates:
                other = small.sites[other_index]
                position = other.fract
                image_translation = (0, 0, 0)
                if image_index:
                    op = small.cell.images[image_index - 1]
                    position = op.apply(position)
                    origin = op.apply(gemmi.Fractional(0.0, 0.0, 0.0))
                    image_translation = (round(origin.x), round(origin.y), round(origin.z))
                for shift, distance in images_within(small.cell, site.fract, position, search_radius):
                    # Compose the NeighborSearch image with the images_within
                    # shift so the bond-graph lookup sees the true relative
                    # image even for bonds crossing a cell boundary.
                    relative_image = (
                        image_translation[0] + shift[0],
                        image_translation[1] + shift[1],
                        image_translation[2] + shift[2],
                    )
                    if index == other_index and relative_image == (0, 0, 0):
                        continue
                    involves_hydrogen = site.element.is_hydrogen or other.element.is_hydrogen
                    triple = (site.label, other.label, relative_image)
                    if triple in bonded:
                        if distance < thresholds.severe_overlap_bonded_distance:
                            severe_overlap_details.append(
                                self._contact_detail(site, other, relative_image, distance, None, "bonded_1-2")
                            )
                        continue
                    if triple in excluded13:
                        if distance < thresholds.min_nonbonded_heavy_distance:
                            severe_overlap_details.append(
                                self._contact_detail(site, other, relative_image, distance, None, "angle_1-3")
                            )
                        continue
                    if triple in excluded14:
                        floor = (
                            thresholds.min_nonbonded_heavy_distance
                            if involves_hydrogen
                            else thresholds.hard_min_nonbonded_heavy_plain_distance
                        )
                        if distance < floor:
                            severe_overlap_details.append(
                                self._contact_detail(site, other, relative_image, distance, None, "torsion_1-4")
                            )
                        continue
                    assessment = assess_pair(
                        distance,
                        str(site.element.name),
                        str(other.element.name),
                        ratio_threshold=thresholds.min_nonbonded_heavy_vdw_ratio,
                        plain_floor=thresholds.hard_min_nonbonded_heavy_plain_distance,
                    )
                    pairs_scanned += 1
                    if involves_hydrogen:
                        if min_hydrogen_distance is None or distance < min_hydrogen_distance:
                            min_hydrogen_distance = distance
                        if min_hydrogen_ratio is None or assessment.vdw_ratio < min_hydrogen_ratio:
                            min_hydrogen_ratio = assessment.vdw_ratio
                        if assessment.clash:
                            hydrogen_clash_details.append(
                                self._contact_detail(site, other, relative_image, distance, assessment.vdw_ratio, None)
                            )
                    else:
                        if min_heavy_distance is None or distance < min_heavy_distance:
                            min_heavy_distance = distance
                        if min_heavy_ratio is None or assessment.vdw_ratio < min_heavy_ratio:
                            min_heavy_ratio = assessment.vdw_ratio
                        if assessment.clash:
                            heavy_clash_details.append(
                                self._contact_detail(site, other, relative_image, distance, assessment.vdw_ratio, None)
                            )

        # Every contact is found once per direction; report unordered counts.
        return {
            "min_nonbonded_heavy_distance": min_heavy_distance,
            "min_nonbonded_heavy_vdw_ratio": min_heavy_ratio,
            "min_nonbonded_hydrogen_distance": min_hydrogen_distance,
            "min_nonbonded_hydrogen_vdw_ratio": min_hydrogen_ratio,
            "n_nonbonded_pairs_scanned": pairs_scanned // 2,
            "n_heavy_atom_clash_pairs": len(heavy_clash_details) // 2,
            "heavy_atom_clash_min": min(heavy_clash_details, key=lambda item: item["vdw_ratio"], default=None),
            "n_hydrogen_atom_clash_pairs": len(hydrogen_clash_details) // 2,
            "hydrogen_atom_clash_min": min(hydrogen_clash_details, key=lambda item: item["vdw_ratio"], default=None),
            "n_excluded_pair_severe_overlaps": len(severe_overlap_details) // 2,
            "excluded_pair_severe_overlap_min": min(
                severe_overlap_details, key=lambda item: item["distance"], default=None
            ),
        }

    @staticmethod
    def _contact_detail(site, other, relative_image, distance, vdw_ratio, exclusion) -> dict[str, object]:
        detail: dict[str, object] = {
            "labels": (str(site.label), str(other.label)),
            "image": tuple(int(v) for v in relative_image),
            "distance": float(distance),
        }
        if vdw_ratio is not None:
            detail["vdw_ratio"] = float(vdw_ratio)
        if exclusion is not None:
            detail["exclusion"] = exclusion
        return detail

    def _resolve_cif_path(
        self,
        record: Mapping[str, object],
        *,
        source_root: str | Path | None,
    ) -> Path | None:
        raw = self._string(record.get("cif_path"))
        if not raw:
            return None
        path = Path(raw)
        if path.is_absolute():
            return path
        if source_root is not None:
            root = Path(source_root)
            candidate = root / path
            if candidate.exists():
                return candidate
            candidate = root.parent / path
            if candidate.exists():
                return candidate
        return path

    def _topology_dimensionality(self, topology_id: str | None) -> str | None:
        if not topology_id:
            return None
        try:
            return get_topology_hint(topology_id).dimensionality
        except KeyError:
            return None
        except Exception as exc:
            print(
                f"warning: topology dimensionality lookup failed for {topology_id!r}; "
                f"falling back to the 3D cell-volume check: {type(exc).__name__}: {exc}",
                file=sys.stderr,
            )
            return None

    def _instance_id(self, atom_label: str) -> str:
        return _instance_id(atom_label)

    def _n_unreacted_motifs(
        self,
        record: Mapping[str, object],
        score_metadata: Mapping[str, object],
    ) -> int:
        raw = score_metadata.get("n_unreacted_motifs")
        if raw is not None:
            try:
                return int(raw)
            except (TypeError, ValueError):
                return 0
        for flag in record.get("flags", ()):
            text = str(flag)
            if text.startswith("unreacted_motifs:"):
                try:
                    return int(text.split(":", 1)[1])
                except ValueError:
                    return 0
        return 0

    def _distance_residual(self, item: Mapping[str, object]) -> float | None:
        raw = item.get("distance_residual")
        if raw is not None:
            return float(raw)
        actual = item.get("actual_distance")
        target = item.get("target_distance")
        if actual is None or target is None:
            # Missing data must not masquerade as a perfect bridge.
            return None
        return abs(float(actual) - float(target))

    def _actual_distance_bounds(self, item: Mapping[str, object]) -> tuple[float, float] | None:
        actual = item.get("actual_distance")
        target = item.get("target_distance")
        if actual is None or target is None:
            return None
        target_value = float(target)
        if target_value <= 0.0:
            return None
        actual_value = float(actual)
        ratio = actual_value / target_value
        return ratio, ratio

    def _mapping(self, value: object) -> Mapping[str, object]:
        return value if isinstance(value, Mapping) else {}

    def _string(self, value: object) -> str | None:
        if value is None:
            return None
        text = str(value)
        return text if text else None


def classify_batch_output(
    source_dir: str | Path,
    output_dir: str | Path,
    *,
    thresholds: CoarseValidationThresholds | None = None,
    link_mode: str = "symlink",
    max_workers: int | None = None,
    max_structures: int | None = None,
) -> BatchOutputClassificationSummary:
    validator = CoarseStructureValidator(thresholds=thresholds)
    source_root = Path(source_dir)
    manifest_path = source_root / "manifest.jsonl"
    if not manifest_path.is_file():
        raise FileNotFoundError(f"expected manifest at {manifest_path}")

    output_root = Path(output_dir)
    output_root.mkdir(parents=True, exist_ok=True)
    valid_root = output_root / "valid" / "cifs"
    warning_root = output_root / "warning"
    needs_optimization_root = output_root / "needs_optimization"
    hard_hard_invalid_root = output_root / "hard_hard_invalid"
    hard_invalid_root = output_root / "hard_invalid"
    unvalidated_root = output_root / "unvalidated"
    valid_root.mkdir(parents=True, exist_ok=True)
    warning_root.mkdir(parents=True, exist_ok=True)
    needs_optimization_root.mkdir(parents=True, exist_ok=True)
    hard_hard_invalid_root.mkdir(parents=True, exist_ok=True)
    hard_invalid_root.mkdir(parents=True, exist_ok=True)
    unvalidated_root.mkdir(parents=True, exist_ok=True)

    classification_manifest_path = output_root / "classification_manifest.jsonl"
    valid_manifest_path = output_root / "valid" / "manifest.jsonl"
    warning_manifest_path = output_root / "warning" / "manifest.jsonl"
    needs_optimization_manifest_path = output_root / "needs_optimization" / "manifest.jsonl"
    hard_hard_invalid_manifest_path = output_root / "hard_hard_invalid" / "manifest.jsonl"
    hard_invalid_manifest_path = output_root / "hard_invalid" / "manifest.jsonl"
    unvalidated_manifest_path = output_root / "unvalidated" / "manifest.jsonl"

    max_workers = max_workers or min(8, os.cpu_count() or 1)
    pending: dict[Future[tuple[dict[str, object], CoarseValidationReport]], dict[str, object]] = {}
    warning_reason_counts: Counter[str] = Counter()
    needs_optimization_reason_counts: Counter[str] = Counter()
    hard_hard_invalid_reason_counts: Counter[str] = Counter()
    hard_invalid_reason_counts: Counter[str] = Counter()
    unvalidated_reason_counts: Counter[str] = Counter()
    total_structures = 0
    valid_structures = 0
    warning_structures = 0
    needs_optimization_structures = 0
    hard_hard_invalid_structures = 0
    hard_invalid_structures = 0
    unvalidated_structures = 0

    with (
        manifest_path.open(encoding="utf-8") as manifest,
        classification_manifest_path.open("w", encoding="utf-8") as classification_manifest,
        valid_manifest_path.open("w", encoding="utf-8") as valid_manifest,
        warning_manifest_path.open("w", encoding="utf-8") as warning_manifest,
        needs_optimization_manifest_path.open("w", encoding="utf-8") as needs_optimization_manifest,
        hard_hard_invalid_manifest_path.open("w", encoding="utf-8") as hard_hard_invalid_manifest,
        hard_invalid_manifest_path.open("w", encoding="utf-8") as hard_invalid_manifest,
        unvalidated_manifest_path.open("w", encoding="utf-8") as unvalidated_manifest,
        ThreadPoolExecutor(max_workers=max_workers) as executor,
    ):
        for line in manifest:
            if max_structures is not None and total_structures >= max_structures:
                break
            if not line.strip():
                continue
            record = json.loads(line)
            if str(record.get("status")) != "ok":
                continue
            future = executor.submit(_validate_record_worker, record, str(source_root), validator.thresholds)
            pending[future] = record
            total_structures += 1
            if len(pending) >= max_workers * 4:
                (
                    valid_structures,
                    warning_structures,
                    needs_optimization_structures,
                    hard_hard_invalid_structures,
                    hard_invalid_structures,
                    unvalidated_structures,
                ) = _drain_completed(
                    pending,
                    classification_manifest,
                    valid_manifest,
                    warning_manifest,
                    needs_optimization_manifest,
                    hard_hard_invalid_manifest,
                    hard_invalid_manifest,
                    unvalidated_manifest,
                    valid_root,
                    warning_root,
                    needs_optimization_root,
                    hard_hard_invalid_root,
                    hard_invalid_root,
                    unvalidated_root,
                    warning_reason_counts,
                    needs_optimization_reason_counts,
                    hard_hard_invalid_reason_counts,
                    hard_invalid_reason_counts,
                    unvalidated_reason_counts,
                    link_mode,
                    valid_count=valid_structures,
                    warning_count=warning_structures,
                    needs_optimization_count=needs_optimization_structures,
                    hard_hard_invalid_count=hard_hard_invalid_structures,
                    hard_invalid_count=hard_invalid_structures,
                    unvalidated_count=unvalidated_structures,
                )

        while pending:
            (
                valid_structures,
                warning_structures,
                needs_optimization_structures,
                hard_hard_invalid_structures,
                hard_invalid_structures,
                unvalidated_structures,
            ) = _drain_completed(
                pending,
                classification_manifest,
                valid_manifest,
                warning_manifest,
                needs_optimization_manifest,
                hard_hard_invalid_manifest,
                hard_invalid_manifest,
                unvalidated_manifest,
                valid_root,
                warning_root,
                needs_optimization_root,
                hard_hard_invalid_root,
                hard_invalid_root,
                unvalidated_root,
                warning_reason_counts,
                needs_optimization_reason_counts,
                hard_hard_invalid_reason_counts,
                hard_invalid_reason_counts,
                unvalidated_reason_counts,
                link_mode,
                valid_count=valid_structures,
                warning_count=warning_structures,
                needs_optimization_count=needs_optimization_structures,
                hard_hard_invalid_count=hard_hard_invalid_structures,
                hard_invalid_count=hard_invalid_structures,
                unvalidated_count=unvalidated_structures,
            )

    summary = BatchOutputClassificationSummary(
        source_dir=str(source_root),
        output_dir=str(output_root),
        total_structures=total_structures,
        valid_structures=valid_structures,
        warning_structures=warning_structures,
        needs_optimization_structures=needs_optimization_structures,
        hard_hard_invalid_structures=hard_hard_invalid_structures,
        hard_invalid_structures=hard_invalid_structures,
        warning_reason_counts=dict(warning_reason_counts),
        needs_optimization_reason_counts=dict(needs_optimization_reason_counts),
        hard_hard_invalid_reason_counts=dict(hard_hard_invalid_reason_counts),
        hard_invalid_reason_counts=dict(hard_invalid_reason_counts),
        unvalidated_structures=unvalidated_structures,
        unvalidated_reason_counts=dict(unvalidated_reason_counts),
        classification_manifest_path=str(classification_manifest_path),
        valid_manifest_path=str(valid_manifest_path),
        warning_manifest_path=str(warning_manifest_path),
        needs_optimization_manifest_path=str(needs_optimization_manifest_path),
        hard_hard_invalid_manifest_path=str(hard_hard_invalid_manifest_path),
        hard_invalid_manifest_path=str(hard_invalid_manifest_path),
        unvalidated_manifest_path=str(unvalidated_manifest_path),
    )
    _write_classification_summary(output_root / "summary.md", summary, validator.thresholds)
    return summary


def _drain_completed(
    pending: dict[Future[tuple[dict[str, object], CoarseValidationReport]], dict[str, object]],
    classification_manifest,
    valid_manifest,
    warning_manifest,
    needs_optimization_manifest,
    hard_hard_invalid_manifest,
    hard_invalid_manifest,
    unvalidated_manifest,
    valid_root: Path,
    warning_root: Path,
    needs_optimization_root: Path,
    hard_hard_invalid_root: Path,
    hard_invalid_root: Path,
    unvalidated_root: Path,
    warning_reason_counts: Counter[str],
    needs_optimization_reason_counts: Counter[str],
    hard_hard_invalid_reason_counts: Counter[str],
    hard_invalid_reason_counts: Counter[str],
    unvalidated_reason_counts: Counter[str],
    link_mode: str,
    *,
    valid_count: int,
    warning_count: int,
    needs_optimization_count: int,
    hard_hard_invalid_count: int,
    hard_invalid_count: int,
    unvalidated_count: int,
) -> tuple[int, int, int, int, int, int]:
    done, _not_done = wait(set(pending), return_when=FIRST_COMPLETED)
    for future in done:
        record, report = future.result()
        del pending[future]
        enriched = dict(record)
        enriched["validation"] = {
            "classification": report.classification,
            "is_valid": report.is_valid,
            "passes_hard_validation": report.passes_hard_validation,
            "blocks_cif_export": report.blocks_cif_export,
            "reasons": list(report.reasons),
            "hard_hard_invalid_reasons": list(report.hard_hard_invalid_reasons),
            "warning_reasons": list(report.warning_reasons),
            "hard_invalid_reasons": list(report.hard_invalid_reasons),
            "needs_optimization_reasons": list(report.needs_optimization_reasons),
            "unmeasured_required_checks": list(report.unmeasured_required_checks),
            "coverage": dict(report.coverage),
            "metrics": _json_safe(report.metrics),
        }
        classification_manifest.write(json.dumps(_json_safe(enriched), sort_keys=True) + "\n")
        cif_path = report.metrics.get("cif_path")
        source = Path(str(cif_path)) if cif_path is not None else None
        if report.classification == "unvalidated":
            unvalidated_count += 1
            unvalidated_manifest.write(json.dumps(_json_safe(enriched), sort_keys=True) + "\n")
            for reason in report.reasons:
                unvalidated_reason_counts[reason] += 1
            if source is not None and source.is_file():
                _materialize_link(source, unvalidated_root / "cifs" / source.name, link_mode)
            continue
        if report.is_valid:
            valid_count += 1
            valid_manifest.write(json.dumps(_json_safe(enriched), sort_keys=True) + "\n")
            if source is not None and source.is_file():
                _materialize_link(source, valid_root / source.name, link_mode)
            continue
        if report.classification == "warning":
            warning_count += 1
            warning_manifest.write(json.dumps(_json_safe(enriched), sort_keys=True) + "\n")
            for reason in report.warning_reasons:
                warning_reason_counts[reason] += 1
            if source is not None and source.is_file():
                _materialize_link(source, warning_root / "cifs" / source.name, link_mode)
                for reason in report.warning_reasons:
                    _materialize_link(source, warning_root / "reasons" / reason / source.name, link_mode)
            continue
        if report.classification == "needs_optimization":
            needs_optimization_count += 1
            needs_optimization_manifest.write(json.dumps(_json_safe(enriched), sort_keys=True) + "\n")
            for reason in report.needs_optimization_reasons:
                needs_optimization_reason_counts[reason] += 1
            if source is not None and source.is_file():
                _materialize_link(source, needs_optimization_root / "cifs" / source.name, link_mode)
                for reason in report.needs_optimization_reasons:
                    _materialize_link(source, needs_optimization_root / "reasons" / reason / source.name, link_mode)
            continue
        if report.classification == "hard_hard_invalid":
            hard_hard_invalid_count += 1
            hard_hard_invalid_manifest.write(json.dumps(_json_safe(enriched), sort_keys=True) + "\n")
            for reason in report.hard_hard_invalid_reasons:
                hard_hard_invalid_reason_counts[reason] += 1
            if source is not None and source.is_file():
                _materialize_link(source, hard_hard_invalid_root / "cifs" / source.name, link_mode)
                for reason in report.hard_hard_invalid_reasons:
                    _materialize_link(source, hard_hard_invalid_root / "reasons" / reason / source.name, link_mode)
            continue
        hard_invalid_count += 1
        hard_invalid_manifest.write(json.dumps(_json_safe(enriched), sort_keys=True) + "\n")
        for reason in report.hard_invalid_reasons:
            hard_invalid_reason_counts[reason] += 1
        if source is not None and source.is_file():
            _materialize_link(source, hard_invalid_root / "cifs" / source.name, link_mode)
            for reason in report.hard_invalid_reasons:
                _materialize_link(source, hard_invalid_root / "reasons" / reason / source.name, link_mode)
    return valid_count, warning_count, needs_optimization_count, hard_hard_invalid_count, hard_invalid_count, unvalidated_count


def _validate_record_worker(
    record: dict[str, object],
    source_root: str,
    thresholds: CoarseValidationThresholds,
) -> tuple[dict[str, object], CoarseValidationReport]:
    validator = CoarseStructureValidator(thresholds=thresholds)
    return record, validator.validate_manifest_record(record, source_root=source_root)


def _materialize_link(source: Path, destination: Path, mode: str) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() or destination.is_symlink():
        return
    if mode == "symlink":
        destination.symlink_to(source.resolve())
        return
    if mode == "hardlink":
        os.link(source, destination)
        return
    if mode == "copy":
        shutil.copy2(source, destination)
        return
    raise ValueError(f"unsupported link mode {mode!r}")


def _write_classification_summary(
    path: Path,
    summary: BatchOutputClassificationSummary,
    thresholds: CoarseValidationThresholds,
) -> None:
    lines = [
        "# Batch output coarse validation summary",
        "",
        f"- Source directory: `{summary.source_dir}`",
        f"- Output directory: `{summary.output_dir}`",
        f"- Total structures classified: {summary.total_structures}",
        f"- Valid structures: {summary.valid_structures}",
        f"- Warning structures: {summary.warning_structures}",
        f"- Needs-optimization structures: {summary.needs_optimization_structures}",
        f"- Hard-hard-invalid structures: {summary.hard_hard_invalid_structures}",
        f"- Hard-invalid structures: {summary.hard_invalid_structures}",
        f"- Unvalidated structures (required checks unmeasured): {summary.unvalidated_structures}",
        f"- Classification manifest: `{summary.classification_manifest_path}`",
        f"- Valid manifest: `{summary.valid_manifest_path}`",
        f"- Warning manifest: `{summary.warning_manifest_path}`",
        f"- Needs-optimization manifest: `{summary.needs_optimization_manifest_path}`",
        f"- Hard-hard-invalid manifest: `{summary.hard_hard_invalid_manifest_path}`",
        f"- Hard-invalid manifest: `{summary.hard_invalid_manifest_path}`",
        f"- Unvalidated manifest: `{summary.unvalidated_manifest_path}`",
        "",
        "## Thresholds",
        "",
        f"- warning max bridge distance residual: {thresholds.warning_max_bridge_distance_residual:.3f} A",
        f"- warning mean bridge distance residual: {thresholds.warning_mean_bridge_distance_residual:.3f} A",
        f"- warning bad bridge event cutoff: {thresholds.warning_bad_bridge_distance_residual:.3f} A",
        f"- warning bad bridge event fraction cutoff: {thresholds.warning_max_bad_bridge_fraction:.3f}",
        f"- hard-hard max bridge distance: {thresholds.hard_hard_max_bridge_distance:.3f} A",
        f"- hard max bridge distance residual: {thresholds.hard_max_bridge_distance_residual:.3f} A",
        f"- hard mean bridge distance residual: {thresholds.hard_mean_bridge_distance_residual:.3f} A",
        f"- hard min bridge distance ratio: {thresholds.hard_min_bridge_distance_ratio:.3f}",
        f"- hard max bridge distance ratio: {thresholds.hard_max_bridge_distance_ratio:.3f}",
        f"- minimum nonbonded heavy-atom distance: {thresholds.min_nonbonded_heavy_distance:.3f} A",
        f"- nonbonded clash vdW-ratio threshold: {thresholds.min_nonbonded_heavy_vdw_ratio:.3f}",
        f"- nonbonded contact search radius: {thresholds.nonbonded_heavy_search_radius:.3f} A",
        f"- heavy severe-overlap plain-distance floor: {thresholds.hard_min_nonbonded_heavy_plain_distance:.3f} A",
        f"- bonded-pair severe-overlap floor: {thresholds.severe_overlap_bonded_distance:.3f} A",
        f"- minimum 2D cell area: {thresholds.min_2d_cell_area:.3f} A^2",
        f"- minimum 3D cell volume: {thresholds.min_3d_cell_volume:.3f} A^3",
        "",
        "## Warning Reason Counts",
        "",
    ]
    if not summary.warning_reason_counts:
        lines.append("- No warning structures were detected.")
    else:
        for reason, count in sorted(summary.warning_reason_counts.items(), key=lambda item: (-item[1], item[0])):
            lines.append(f"- `{reason}`: {count}")
    lines.extend(
        [
            "",
            "## Needs-Optimization Reason Counts",
            "",
        ]
    )
    if not summary.needs_optimization_reason_counts:
        lines.append("- No needs-optimization structures were detected.")
    else:
        for reason, count in sorted(summary.needs_optimization_reason_counts.items(), key=lambda item: (-item[1], item[0])):
            lines.append(f"- `{reason}`: {count}")
    lines.extend(
        [
            "",
            "## Hard-Hard-Invalid Reason Counts",
            "",
        ]
    )
    if not summary.hard_hard_invalid_reason_counts:
        lines.append("- No hard-hard-invalid structures were detected.")
    else:
        for reason, count in sorted(summary.hard_hard_invalid_reason_counts.items(), key=lambda item: (-item[1], item[0])):
            lines.append(f"- `{reason}`: {count}")
    lines.extend(
        [
            "",
            "## Hard-Invalid Reason Counts",
            "",
        ]
    )
    if not summary.hard_invalid_reason_counts:
        lines.append("- No hard-invalid structures were detected.")
    else:
        for reason, count in sorted(summary.hard_invalid_reason_counts.items(), key=lambda item: (-item[1], item[0])):
            lines.append(f"- `{reason}`: {count}")
    lines.extend(
        [
            "",
            "## Unvalidated Reason Counts",
            "",
        ]
    )
    if not summary.unvalidated_reason_counts:
        lines.append("- No unvalidated structures were detected.")
    else:
        for reason, count in sorted(summary.unvalidated_reason_counts.items(), key=lambda item: (-item[1], item[0])):
            lines.append(f"- `{reason}`: {count}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _json_safe(value: object) -> object:
    if isinstance(value, Mapping):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    return value


__all__ = [
    "BatchOutputClassificationSummary",
    "CHECK_MEASURED",
    "CHECK_MISSING_DATA",
    "CHECK_NOT_APPLICABLE",
    "CHECK_NO_CONTACTS",
    "CHECK_SKIPPED",
    "CoarseStructureValidator",
    "CoarseValidationReport",
    "CoarseValidationThresholds",
    "REQUIRED_COVERAGE_CHECKS",
    "RealizedBridgeBondMeasurement",
    "classify_batch_output",
    "measure_inter_instance_bond_distances",
]
