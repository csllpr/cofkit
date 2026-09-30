from __future__ import annotations

from collections.abc import Mapping as ABCMapping, Sequence as ABCSequence
from dataclasses import dataclass, field
from typing import Mapping

from .model import MonomerSpec


@dataclass(frozen=True)
class BatchMonomerRecord:
    id: str
    name: str
    smiles: str
    motif_kind: str
    expected_connectivity: int
    source_path: str = ""
    source_line: int = 0
    metadata: Mapping[str, object] = field(default_factory=dict)


@dataclass(frozen=True)
class BuiltBatchMonomer:
    record: BatchMonomerRecord
    monomer: MonomerSpec | None = None
    error: str | None = None

    @property
    def ok(self) -> bool:
        return self.monomer is not None and self.error is None


@dataclass(frozen=True)
class BatchPairSummary:
    structure_id: str
    pair_id: str
    pair_mode: str
    status: str
    reactant_a_record_id: str
    reactant_b_record_id: str
    reactant_a_connectivity: int = 0
    reactant_b_connectivity: int = 0
    topology_id: str | None = None
    score: float | None = None
    flags: tuple[str, ...] = ()
    cif_path: str | None = None
    metadata: Mapping[str, object] = field(default_factory=dict)

    @property
    def reactant_roles(self) -> tuple[str, ...]:
        value = self.metadata.get("reactant_roles")
        if isinstance(value, ABCSequence) and not isinstance(value, (str, bytes)):
            roles = tuple(str(item) for item in value)
            if len(roles) == 2:
                return roles
        return ("reactant_a", "reactant_b")

    @property
    def reactant_record_ids(self) -> Mapping[str, str]:
        mapping = self.metadata.get("reactant_record_ids")
        if isinstance(mapping, ABCMapping):
            return {str(key): str(value) for key, value in mapping.items()}
        first_role, second_role = self.reactant_roles
        return {
            first_role: self.reactant_a_record_id,
            second_role: self.reactant_b_record_id,
        }

    @property
    def reactant_connectivities(self) -> Mapping[str, int]:
        mapping = self.metadata.get("reactant_connectivities")
        if isinstance(mapping, ABCMapping):
            return {str(key): int(value) for key, value in mapping.items()}
        first_role, second_role = self.reactant_roles
        return {
            first_role: self.reactant_a_connectivity,
            second_role: self.reactant_b_connectivity,
        }

    @property
    def reactant_node_shapes(self) -> Mapping[str, str]:
        """4-connecting reactant node-shape labels (square/rectangular/etc.)."""
        mapping = self.metadata.get("reactant_node_shapes")
        if isinstance(mapping, ABCMapping):
            return {str(key): str(value) for key, value in mapping.items()}
        return {}

    @property
    def amine_record_id(self) -> str:
        return str(self.reactant_record_ids.get("amine", self.reactant_a_record_id))

    @property
    def aldehyde_record_id(self) -> str:
        return str(self.reactant_record_ids.get("aldehyde", self.reactant_b_record_id))

    @property
    def amine_connectivity(self) -> int:
        return int(self.reactant_connectivities.get("amine", self.reactant_a_connectivity))

    @property
    def aldehyde_connectivity(self) -> int:
        return int(self.reactant_connectivities.get("aldehyde", self.reactant_b_connectivity))

    @property
    def validation_classification(self) -> str | None:
        """Coarse-validation classification of the exported structure.

        One of valid / warning / needs_optimization / hard_invalid /
        hard_hard_invalid / unvalidated; None when no validation record was
        attached (e.g. CIF export disabled or blocked before validation).
        """
        validation = self.metadata.get("validation")
        if isinstance(validation, ABCMapping):
            value = validation.get("classification")
            return str(value) if value is not None else None
        return None

    @property
    def validation_coverage(self) -> Mapping[str, str]:
        """Per-check coverage statuses from validation.

        Values are the ``cofkit.validation`` coverage statuses:
        ``measured``, ``no_contacts`` (a completed search found nothing),
        ``missing_data`` (a required check could not be evaluated),
        ``not_applicable``, or ``skipped``.
        """
        validation = self.metadata.get("validation")
        if isinstance(validation, ABCMapping):
            coverage = validation.get("coverage")
            if isinstance(coverage, ABCMapping):
                return {str(key): str(value) for key, value in coverage.items()}
        return {}

    @property
    def unmeasured_required_checks(self) -> tuple[str, ...]:
        """Required validation checks that could not be measured (missing data)."""
        validation = self.metadata.get("validation")
        if isinstance(validation, ABCMapping):
            checks = validation.get("unmeasured_required_checks")
            if isinstance(checks, ABCSequence) and not isinstance(checks, (str, bytes)):
                return tuple(str(check) for check in checks)
        return ()


@dataclass(frozen=True)
class BatchRunSummary:
    input_dir: str
    output_dir: str
    attempted_pairs: int
    successful_pairs: int
    attempted_structures: int
    successful_structures: int
    cifs_written: int
    built_monomers: int
    failed_monomers: int
    build_failures: Mapping[str, str] = field(default_factory=dict)
    # Every non-"ok" manifest record, structure_id -> "TypeName: message".
    # Covers pair/structure-level failures (unsupported connectivity,
    # generation, realization/CIF export, pair-task failures); monomer build
    # failures are additionally keyed by monomer record in build_failures.
    record_failures: Mapping[str, str] = field(default_factory=dict)
    mode_counts: Mapping[str, int] = field(default_factory=dict)
    topology_counts: Mapping[str, int] = field(default_factory=dict)
    # Per-classification counts of validation outcomes for structures with a
    # validation record (valid / warning / needs_optimization / hard_invalid /
    # hard_hard_invalid / unvalidated). Construction success (status == "ok")
    # is not a validated yield; these counts are the honest destination.
    validation_counts: Mapping[str, int] = field(default_factory=dict)
    geometry_repair_counts: Mapping[str, int] = field(default_factory=dict)
    geometry_repair_revalidation_counts: Mapping[str, int] = field(default_factory=dict)
    geometry_repair_failed_records_path: str | None = None
    manifest_path: str | None = None
    top_results: tuple[BatchPairSummary, ...] = ()
