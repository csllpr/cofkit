from __future__ import annotations

import json
import os
import shutil
from collections import Counter, defaultdict
from concurrent.futures import FIRST_COMPLETED, Future, ThreadPoolExecutor, wait
from dataclasses import dataclass, field
from math import radians, sin
from pathlib import Path
from typing import Any, Mapping

try:  # pragma: no cover - exercised in integration environments
    import gemmi
except ImportError:  # pragma: no cover - optional dependency
    gemmi = None

from .topologies import get_topology_hint


@dataclass(frozen=True)
class CoarseValidationThresholds:
    warning_max_bridge_distance_residual: float = 0.75
    warning_mean_bridge_distance_residual: float = 0.35
    warning_bad_bridge_distance_residual: float = 0.50
    warning_max_bad_bridge_fraction: float = 0.25
    hard_hard_max_bridge_distance: float = 2.5
    hard_max_bridge_distance_residual: float = 1.00
    hard_mean_bridge_distance_residual: float = 0.60
    hard_min_bridge_distance_ratio: float = 0.70
    hard_max_bridge_distance_ratio: float = 1.60
    min_nonbonded_heavy_distance: float = 1.05
    min_2d_cell_area: float = 10.0
    min_3d_cell_volume: float = 20.0
    skip_cif_checks_when_metadata_invalid: bool = True


@dataclass(frozen=True)
class CoarseValidationReport:
    classification: str
    is_valid: bool
    passes_hard_validation: bool
    blocks_cif_export: bool = False
    reasons: tuple[str, ...] = ()
    hard_hard_invalid_reasons: tuple[str, ...] = ()
    warning_reasons: tuple[str, ...] = ()
    hard_invalid_reasons: tuple[str, ...] = ()
    metrics: Mapping[str, object] = field(default_factory=dict)


@dataclass(frozen=True)
class BatchOutputClassificationSummary:
    source_dir: str
    output_dir: str
    total_structures: int
    valid_structures: int
    warning_structures: int
    hard_hard_invalid_structures: int
    hard_invalid_structures: int
    warning_reason_counts: Mapping[str, int] = field(default_factory=dict)
    hard_hard_invalid_reason_counts: Mapping[str, int] = field(default_factory=dict)
    hard_invalid_reason_counts: Mapping[str, int] = field(default_factory=dict)
    classification_manifest_path: str = ""
    valid_manifest_path: str = ""
    warning_manifest_path: str = ""
    hard_hard_invalid_manifest_path: str = ""
    hard_invalid_manifest_path: str = ""

    @property
    def invalid_structures(self) -> int:
        return self.hard_hard_invalid_structures + self.hard_invalid_structures


class CoarseStructureValidator:
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

        metadata = self._mapping(record.get("metadata"))
        score_metadata = self._mapping(metadata.get("score_metadata"))
        n_unreacted_motifs = self._n_unreacted_motifs(record, score_metadata)
        metrics["n_unreacted_motifs"] = n_unreacted_motifs
        if n_unreacted_motifs > 0:
            hard_invalid_reasons.append("unreacted_motifs")

        bridge_metrics = tuple(self._mapping(item) for item in score_metadata.get("bridge_event_metrics", ()))
        distance_residuals = tuple(self._distance_residual(item) for item in bridge_metrics)
        actual_distance_bounds = tuple(self._actual_distance_bounds(item) for item in bridge_metrics)
        if distance_residuals:
            max_residual = max(distance_residuals)
            mean_residual = sum(distance_residuals) / len(distance_residuals)
            bad_fraction = sum(
                value > self.thresholds.warning_bad_bridge_distance_residual for value in distance_residuals
            ) / len(distance_residuals)
        else:
            max_residual = 0.0
            mean_residual = 0.0
            bad_fraction = 0.0
        metrics["n_bridge_events"] = len(distance_residuals)
        metrics["max_bridge_distance_residual"] = max_residual
        metrics["mean_bridge_distance_residual"] = mean_residual
        metrics["bad_bridge_event_fraction"] = bad_fraction
        if max_residual > self.thresholds.warning_max_bridge_distance_residual:
            warning_reasons.append("bridge_distance_residual_max")
        if mean_residual > self.thresholds.warning_mean_bridge_distance_residual:
            warning_reasons.append("bridge_distance_residual_mean")
        if bad_fraction > self.thresholds.warning_max_bad_bridge_fraction:
            warning_reasons.append("bridge_distance_residual_fraction")
        if max_residual > self.thresholds.hard_max_bridge_distance_residual:
            hard_invalid_reasons.append("bridge_distance_residual_max_hard")
        if mean_residual > self.thresholds.hard_mean_bridge_distance_residual:
            hard_invalid_reasons.append("bridge_distance_residual_mean_hard")
        min_ratio = None
        max_ratio = None
        max_actual_distance = None
        if actual_distance_bounds:
            finite_bounds = tuple(bounds for bounds in actual_distance_bounds if bounds is not None)
            if finite_bounds:
                min_ratio = min(bounds[0] for bounds in finite_bounds)
                max_ratio = max(bounds[1] for bounds in finite_bounds)
        actual_distances = tuple(
            float(item.get("actual_distance"))
            for item in bridge_metrics
            if item.get("actual_distance") is not None
        )
        if actual_distances:
            max_actual_distance = max(actual_distances)
        metrics["min_bridge_distance_ratio"] = min_ratio
        metrics["max_bridge_distance_ratio"] = max_ratio
        metrics["max_actual_bridge_distance"] = max_actual_distance
        if (
            max_actual_distance is not None
            and max_actual_distance >= self.thresholds.hard_hard_max_bridge_distance
        ):
            hard_hard_invalid_reasons.append("bridge_distance_exceeds_cif_export_limit")
        if min_ratio is not None and min_ratio < self.thresholds.hard_min_bridge_distance_ratio:
            hard_invalid_reasons.append("bridge_distance_too_short")
        if max_ratio is not None and max_ratio > self.thresholds.hard_max_bridge_distance_ratio:
            hard_invalid_reasons.append("bridge_distance_too_long")

        cif_path = self._resolve_cif_path(record, source_root=source_root)
        metrics["cif_path"] = str(cif_path) if cif_path is not None else None
        if cif_path is None or not cif_path.is_file():
            hard_invalid_reasons.append("cif_missing")
        elif not (self.thresholds.skip_cif_checks_when_metadata_invalid and hard_invalid_reasons):
            cif_metrics, cif_reasons = self._validate_cif(cif_path, topology_id=self._string(record.get("topology_id")))
            metrics.update(cif_metrics)
            hard_invalid_reasons.extend(cif_reasons)
        else:
            metrics["cif_checks_skipped"] = True

        normalized_warning_reasons = tuple(dict.fromkeys(warning_reasons))
        normalized_hard_hard_invalid_reasons = tuple(dict.fromkeys(hard_hard_invalid_reasons))
        normalized_hard_invalid_reasons = tuple(dict.fromkeys(hard_invalid_reasons))
        normalized_reasons = normalized_hard_hard_invalid_reasons + normalized_hard_invalid_reasons + tuple(
            reason for reason in normalized_warning_reasons if reason not in normalized_hard_invalid_reasons
        )
        classification = "valid"
        if normalized_hard_hard_invalid_reasons:
            classification = "hard_hard_invalid"
        elif normalized_hard_invalid_reasons:
            classification = "hard_invalid"
        elif normalized_warning_reasons:
            classification = "warning"
        return CoarseValidationReport(
            classification=classification,
            is_valid=classification == "valid",
            passes_hard_validation=classification not in {"hard_invalid", "hard_hard_invalid"},
            blocks_cif_export=bool(normalized_hard_hard_invalid_reasons),
            reasons=normalized_reasons,
            hard_hard_invalid_reasons=normalized_hard_hard_invalid_reasons,
            warning_reasons=normalized_warning_reasons,
            hard_invalid_reasons=normalized_hard_invalid_reasons,
            metrics=metrics,
        )

    def _validate_cif(
        self,
        cif_path: Path,
        *,
        topology_id: str | None,
    ) -> tuple[dict[str, object], tuple[str, ...]]:
        if gemmi is None:  # pragma: no cover - optional dependency
            raise ModuleNotFoundError(
                "gemmi is required for CIF-backed coarse validation. Install the 'crystal' extra."
            )

        metrics: dict[str, object] = {}
        reasons: list[str] = []
        try:
            block = gemmi.cif.read_file(str(cif_path)).sole_block()
            small = gemmi.make_small_structure_from_block(block)
        except Exception as exc:  # pragma: no cover - defensive
            metrics["cif_parse_error"] = f"{type(exc).__name__}: {exc}"
            return metrics, ("cif_parse_error",)

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
            reasons.append("disconnected_instance_graph")

        clash_distance = self._min_nonbonded_heavy_distance_below_cutoff(small, bonded_pairs)
        metrics["min_nonbonded_heavy_distance"] = clash_distance
        if clash_distance is not None and clash_distance < self.thresholds.min_nonbonded_heavy_distance:
            reasons.append("heavy_atom_clash")

        return metrics, tuple(dict.fromkeys(reasons))

    def _bonded_pairs(self, block) -> set[frozenset[str]]:
        label_1 = block.find_loop("_geom_bond_atom_site_label_1")
        label_2 = block.find_loop("_geom_bond_atom_site_label_2")
        if len(label_1) == 0 or len(label_2) == 0:
            return set()
        return {
            frozenset((str(label_1[index]), str(label_2[index])))
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
        seen: set[str] = set()
        for instance_id in instance_ids:
            if instance_id in seen:
                continue
            components += 1
            stack = [instance_id]
            seen.add(instance_id)
            while stack:
                current = stack.pop()
                for other in adjacency.get(current, ()):
                    if other not in seen:
                        seen.add(other)
                        stack.append(other)
        return {
            "metrics": {
                "n_instance_nodes": len(instance_ids),
                "n_instance_components": components if instance_ids else 0,
                "n_inter_instance_edges": inter_instance_edges,
            }
        }

    def _min_nonbonded_heavy_distance_below_cutoff(
        self,
        small,
        bonded_pairs: set[frozenset[str]],
    ) -> float | None:
        cutoff = self.thresholds.min_nonbonded_heavy_distance
        search = gemmi.NeighborSearch(small, cutoff).populate(include_h=False)
        minimum: float | None = None
        for index, site in enumerate(small.sites):
            if site.element.is_hydrogen:
                continue
            origin = site.orth(small.cell)
            for mark in search.find_site_neighbors(site, min_dist=0.001, max_dist=cutoff):
                other_index = int(mark.atom_idx)
                if other_index <= index:
                    continue
                other = small.sites[other_index]
                if other.element.is_hydrogen:
                    continue
                if frozenset((site.label, other.label)) in bonded_pairs:
                    continue
                distance = origin.dist(mark.pos)
                if minimum is None or distance < minimum:
                    minimum = distance
        return minimum

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
        if path.exists():
            return path
        return path

    def _topology_dimensionality(self, topology_id: str | None) -> str | None:
        if not topology_id:
            return None
        try:
            return get_topology_hint(topology_id).dimensionality
        except Exception:
            return None

    def _instance_id(self, atom_label: str) -> str:
        head, _separator, _tail = str(atom_label).partition("_")
        return head or str(atom_label)

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

    def _distance_residual(self, item: Mapping[str, object]) -> float:
        raw = item.get("distance_residual")
        if raw is not None:
            return float(raw)
        actual = item.get("actual_distance")
        target = item.get("target_distance")
        if actual is None or target is None:
            return 0.0
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
    hard_hard_invalid_root = output_root / "hard_hard_invalid"
    hard_invalid_root = output_root / "hard_invalid"
    valid_root.mkdir(parents=True, exist_ok=True)
    warning_root.mkdir(parents=True, exist_ok=True)
    hard_hard_invalid_root.mkdir(parents=True, exist_ok=True)
    hard_invalid_root.mkdir(parents=True, exist_ok=True)

    classification_manifest_path = output_root / "classification_manifest.jsonl"
    valid_manifest_path = output_root / "valid" / "manifest.jsonl"
    warning_manifest_path = output_root / "warning" / "manifest.jsonl"
    hard_hard_invalid_manifest_path = output_root / "hard_hard_invalid" / "manifest.jsonl"
    hard_invalid_manifest_path = output_root / "hard_invalid" / "manifest.jsonl"

    max_workers = max_workers or min(8, os.cpu_count() or 1)
    pending: dict[Future[tuple[dict[str, object], CoarseValidationReport]], dict[str, object]] = {}
    warning_reason_counts: Counter[str] = Counter()
    hard_hard_invalid_reason_counts: Counter[str] = Counter()
    hard_invalid_reason_counts: Counter[str] = Counter()
    total_structures = 0
    valid_structures = 0
    warning_structures = 0
    hard_hard_invalid_structures = 0
    hard_invalid_structures = 0

    with (
        manifest_path.open(encoding="utf-8") as manifest,
        classification_manifest_path.open("w", encoding="utf-8") as classification_manifest,
        valid_manifest_path.open("w", encoding="utf-8") as valid_manifest,
        warning_manifest_path.open("w", encoding="utf-8") as warning_manifest,
        hard_hard_invalid_manifest_path.open("w", encoding="utf-8") as hard_hard_invalid_manifest,
        hard_invalid_manifest_path.open("w", encoding="utf-8") as hard_invalid_manifest,
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
                valid_structures, warning_structures, hard_hard_invalid_structures, hard_invalid_structures = _drain_completed(
                    pending,
                    classification_manifest,
                    valid_manifest,
                    warning_manifest,
                    hard_hard_invalid_manifest,
                    hard_invalid_manifest,
                    valid_root,
                    warning_root,
                    hard_hard_invalid_root,
                    hard_invalid_root,
                    warning_reason_counts,
                    hard_hard_invalid_reason_counts,
                    hard_invalid_reason_counts,
                    link_mode,
                    valid_count=valid_structures,
                    warning_count=warning_structures,
                    hard_hard_invalid_count=hard_hard_invalid_structures,
                    hard_invalid_count=hard_invalid_structures,
                )

        while pending:
            valid_structures, warning_structures, hard_hard_invalid_structures, hard_invalid_structures = _drain_completed(
                pending,
                classification_manifest,
                valid_manifest,
                warning_manifest,
                hard_hard_invalid_manifest,
                hard_invalid_manifest,
                valid_root,
                warning_root,
                hard_hard_invalid_root,
                hard_invalid_root,
                warning_reason_counts,
                hard_hard_invalid_reason_counts,
                hard_invalid_reason_counts,
                link_mode,
                valid_count=valid_structures,
                warning_count=warning_structures,
                hard_hard_invalid_count=hard_hard_invalid_structures,
                hard_invalid_count=hard_invalid_structures,
            )

    summary = BatchOutputClassificationSummary(
        source_dir=str(source_root),
        output_dir=str(output_root),
        total_structures=total_structures,
        valid_structures=valid_structures,
        warning_structures=warning_structures,
        hard_hard_invalid_structures=hard_hard_invalid_structures,
        hard_invalid_structures=hard_invalid_structures,
        warning_reason_counts=dict(warning_reason_counts),
        hard_hard_invalid_reason_counts=dict(hard_hard_invalid_reason_counts),
        hard_invalid_reason_counts=dict(hard_invalid_reason_counts),
        classification_manifest_path=str(classification_manifest_path),
        valid_manifest_path=str(valid_manifest_path),
        warning_manifest_path=str(warning_manifest_path),
        hard_hard_invalid_manifest_path=str(hard_hard_invalid_manifest_path),
        hard_invalid_manifest_path=str(hard_invalid_manifest_path),
    )
    _write_classification_summary(output_root / "summary.md", summary, validator.thresholds)
    return summary


def _drain_completed(
    pending: dict[Future[tuple[dict[str, object], CoarseValidationReport]], dict[str, object]],
    classification_manifest,
    valid_manifest,
    warning_manifest,
    hard_hard_invalid_manifest,
    hard_invalid_manifest,
    valid_root: Path,
    warning_root: Path,
    hard_hard_invalid_root: Path,
    hard_invalid_root: Path,
    warning_reason_counts: Counter[str],
    hard_hard_invalid_reason_counts: Counter[str],
    hard_invalid_reason_counts: Counter[str],
    link_mode: str,
    *,
    valid_count: int,
    warning_count: int,
    hard_hard_invalid_count: int,
    hard_invalid_count: int,
) -> tuple[int, int, int, int]:
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
            "metrics": _json_safe(report.metrics),
        }
        classification_manifest.write(json.dumps(_json_safe(enriched), sort_keys=True) + "\n")
        cif_path = report.metrics.get("cif_path")
        source = Path(str(cif_path)) if cif_path is not None else None
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
    return valid_count, warning_count, hard_hard_invalid_count, hard_invalid_count


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
        f"- Hard-hard-invalid structures: {summary.hard_hard_invalid_structures}",
        f"- Hard-invalid structures: {summary.hard_invalid_structures}",
        f"- Classification manifest: `{summary.classification_manifest_path}`",
        f"- Valid manifest: `{summary.valid_manifest_path}`",
        f"- Warning manifest: `{summary.warning_manifest_path}`",
        f"- Hard-hard-invalid manifest: `{summary.hard_hard_invalid_manifest_path}`",
        f"- Hard-invalid manifest: `{summary.hard_invalid_manifest_path}`",
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
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _json_safe(value: object) -> object:
    if isinstance(value, Mapping):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    return value


__all__ = [
    "BatchOutputClassificationSummary",
    "CoarseStructureValidator",
    "CoarseValidationReport",
    "CoarseValidationThresholds",
    "classify_batch_output",
]
