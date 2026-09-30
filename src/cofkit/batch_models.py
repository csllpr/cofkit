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
class MonomerConformerProvenance:
    """Effective conformer-construction settings and outcome for one built
    monomer (impact-review claims T1-7 / T3.5 / T3.11, action A18).

    ``requested_num_conformers`` is the user/config-requested budget;
    ``effective_num_conformers`` is the budget actually handed to the builder
    after the shape-aware ensemble floor. The ``actual_*`` / outcome fields
    come from the built monomer's metadata and are None when the build failed
    or a custom builder records no conformer metadata.
    """

    requested_num_conformers: int
    effective_num_conformers: int
    random_seed: int
    shape_aware_requested: bool
    shape_aware_applied: bool
    builder: str
    cofkit_version: str
    actual_num_conformers: int | None = None
    selected_conformer_id: int | None = None
    conformer_selection: str | None = None
    conformer_selection_note: str | None = None
    embedding_method: str | None = None
    embedding_fallback: bool | None = None
    rdkit_version: str | None = None


@dataclass(frozen=True)
class ConformerConstructionSettings:
    """Run-level effective conformer-construction settings (A18).

    Recorded on ``BatchRunSummary.conformer_settings`` and rendered into
    ``summary.md`` / the CLI batch summary so the durable record alone
    explains the two-conformer autodetection clamp, shape-aware ensemble
    increases, and where per-monomer selected-conformer provenance lives.
    """

    requested_num_conformers: int
    autodetect_num_conformers: int
    shape_aware_conformer: bool
    shape_aware_ensemble_floor: int
    shape_aware_min_motifs: int
    shape_aware_validated_max_motifs: int
    random_seed: int
    cofkit_version: str

    def explanation_lines(self) -> tuple[str, ...]:
        if self.shape_aware_conformer:
            shape_line = (
                f"Shape-aware conformer selection: enabled for monomers with "
                f"{self.shape_aware_min_motifs}..{self.shape_aware_validated_max_motifs} motifs; "
                f"their embedding budget is raised to at least "
                f"{self.shape_aware_ensemble_floor} conformers and the conformer whose motif "
                f"origins best form a regular planar polygon is selected instead of the "
                f"lowest-energy one"
            )
        else:
            shape_line = (
                "Shape-aware conformer selection: disabled; the lowest-energy conformer "
                "is selected for every monomer"
            )
        return (
            f"Requested conformer budget per monomer: {self.requested_num_conformers}",
            (
                f"Motif-kind autodetection probes embed at most {self.autodetect_num_conformers} "
                f"conformer(s) (clamped from the requested budget); the actual monomer build "
                f"re-embeds at the full budget above"
            ),
            shape_line,
            f"RDKit embedding random seed: {self.random_seed}",
            f"Implementation: cofkit {self.cofkit_version}",
            (
                "Per-monomer requested/effective budgets, embedded conformer counts, selected "
                "conformer ids, selection mode, and embedding fallback status are recorded in "
                "monomers.jsonl (`conformer_provenance`) and per structure in manifest.jsonl "
                "(`metadata.reactant_conformer_provenance`)"
            ),
        )


@dataclass(frozen=True)
class BuiltBatchMonomer:
    record: BatchMonomerRecord
    monomer: MonomerSpec | None = None
    error: str | None = None
    # Effective conformer-construction settings and build outcome (A18);
    # None for monomers built before this provenance was attached.
    conformer_provenance: MonomerConformerProvenance | None = None

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

    @property
    def reactant_conformer_provenance(self) -> Mapping[str, object]:
        """Per-reactant selected-conformer provenance (A18).

        Keyed by monomer id; each value carries the requested/effective
        conformer budget, embedding seed, selection mode, selected conformer
        id, and embedding fallback status recorded at monomer build time.
        Empty for records that predate this metadata.
        """
        mapping = self.metadata.get("reactant_conformer_provenance")
        if isinstance(mapping, ABCMapping):
            return {str(key): value for key, value in mapping.items()}
        return {}


@dataclass(frozen=True)
class BatchRunSummary:
    input_dir: str
    output_dir: str
    attempted_pairs: int
    successful_pairs: int
    attempted_structures: int
    # Construction successes (status == "ok") only — NOT a validated yield.
    # Use the derived properties below (constructed/exported/screened/
    # unvalidated) for reporting; validation_counts is the quality breakdown.
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
    # Durable per-monomer record ledger (monomers.jsonl): every library record
    # with its source, detection metadata (including autodetect per-kind
    # failure causes and overlap warnings otherwise confined to stderr), and
    # its build outcome.
    monomer_records_path: str | None = None
    # Run-level effective conformer-construction settings (A18): requested vs
    # autodetect budgets, the shape-aware gate/ensemble floor, the embedding
    # seed, and the cofkit version. None for summaries built before this
    # provenance existed.
    conformer_settings: ConformerConstructionSettings | None = None
    top_results: tuple[BatchPairSummary, ...] = ()

    @property
    def constructed_structures(self) -> int:
        """Structures whose construction succeeded (status "ok").

        This is identical to the legacy `successful_structures` field; the
        new name states explicitly that it is a construction count, not a
        validated yield.
        """
        return self.successful_structures

    @property
    def exported_structures(self) -> int:
        """Structures with a written CIF (identical to `cifs_written`)."""
        return self.cifs_written

    @property
    def screened_structures(self) -> int:
        """Constructed structures that received a validation record."""
        return sum(int(count) for count in self.validation_counts.values())

    @property
    def unvalidated_structures(self) -> int:
        """Structures whose required validation checks could not be measured."""
        return int(self.validation_counts.get("unvalidated", 0))
