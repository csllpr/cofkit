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
    mode_counts: Mapping[str, int] = field(default_factory=dict)
    topology_counts: Mapping[str, int] = field(default_factory=dict)
    manifest_path: str | None = None
    top_results: tuple[BatchPairSummary, ...] = ()
