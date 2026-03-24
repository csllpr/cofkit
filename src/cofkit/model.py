from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from typing import Any, Mapping

from .geometry import Frame

Image3 = tuple[int, int, int]


@dataclass(frozen=True)
class ReactiveMotif:
    id: str
    kind: str
    atom_ids: tuple[int, ...]
    frame: Frame
    valence: int = 1
    symmetry_order: int = 1
    planarity_class: str = "free"
    allowed_reaction_templates: tuple[str, ...] = ()
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def accepts(self, template_id: str) -> bool:
        return not self.allowed_reaction_templates or template_id in self.allowed_reaction_templates


@dataclass(frozen=True)
class MonomerSpec:
    id: str
    name: str
    motifs: tuple[ReactiveMotif, ...]
    formal_charge: int = 0
    rigid_groups: tuple[tuple[int, ...], ...] = ()
    rotatable_bonds: tuple[tuple[int, int], ...] = ()
    conformer_ids: tuple[str, ...] = ()
    atom_symbols: tuple[str, ...] = ()
    atom_positions: tuple[tuple[float, float, float], ...] = ()
    bonds: tuple[tuple[int, int, float], ...] = ()
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.atom_symbols and len(self.atom_symbols) != len(self.atom_positions):
            raise ValueError("atom_symbols and atom_positions must have the same length")
        atom_count = len(self.atom_symbols)
        for atom_id_1, atom_id_2, _ in self.bonds:
            if atom_id_1 < 0 or atom_id_2 < 0 or atom_id_1 >= atom_count or atom_id_2 >= atom_count:
                raise ValueError("bond atom ids must reference atom_symbols/atom_positions entries")

    def motif_by_id(self, motif_id: str) -> ReactiveMotif:
        for motif in self.motifs:
            if motif.id == motif_id:
                return motif
        raise KeyError(f"monomer {self.id!r} has no motif {motif_id!r}")

    @property
    def motif_kinds(self) -> tuple[str, ...]:
        return tuple(m.kind for m in self.motifs)

    @property
    def motif_counts(self) -> Counter[str]:
        return Counter(self.motif_kinds)


@dataclass(frozen=True)
class ReactionTemplate:
    id: str
    arity: int
    reactant_motif_kinds: tuple[str, ...]
    product_name: str
    topology_role: str = "bridge"
    planarity_prior: str = "free"
    torsion_prior: str = "free"
    allowed_dimensionalities: tuple[str, ...] = ("2D", "3D")
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def matches(self, motif_kinds: tuple[str, ...]) -> bool:
        if len(motif_kinds) != self.arity:
            return False
        return Counter(motif_kinds) == Counter(self.reactant_motif_kinds)

    def allows_dimensionality(self, dimensionality: str) -> bool:
        return dimensionality in self.allowed_dimensionalities


@dataclass(frozen=True)
class MotifRef:
    monomer_instance_id: str
    monomer_id: str
    motif_id: str
    periodic_image: Image3 = (0, 0, 0)


@dataclass(frozen=True)
class ReactionEvent:
    id: str
    template_id: str
    participants: tuple[MotifRef, ...]
    product_state: str = "default"
    metadata: Mapping[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class MonomerInstance:
    id: str
    monomer_id: str
    periodic_image: Image3 = (0, 0, 0)
    conformer_id: str | None = None
    metadata: Mapping[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class Pose:
    translation: tuple[float, float, float] = (0.0, 0.0, 0.0)
    rotation_matrix: tuple[tuple[float, float, float], ...] = (
        (1.0, 0.0, 0.0),
        (0.0, 1.0, 0.0),
        (0.0, 0.0, 1.0),
    )


@dataclass(frozen=True)
class AssemblyState:
    cell: tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]] = (
        (30.0, 0.0, 0.0),
        (0.0, 30.0, 0.0),
        (0.0, 0.0, 8.0),
    )
    monomer_poses: Mapping[str, Pose] = field(default_factory=dict)
    torsions: Mapping[str, float] = field(default_factory=dict)
    layer_offsets: Mapping[str, tuple[float, float, float]] = field(default_factory=dict)
    stacking_state: str = "unassigned"


@dataclass(frozen=True)
class Candidate:
    id: str
    score: float
    state: AssemblyState
    events: tuple[ReactionEvent, ...]
    flags: tuple[str, ...] = ()
    metadata: Mapping[str, Any] = field(default_factory=dict)


@dataclass
class CandidateEnsemble:
    candidates: list[Candidate] = field(default_factory=list)

    def add(self, candidate: Candidate) -> None:
        self.candidates.append(candidate)
        self.candidates.sort(key=lambda c: c.score, reverse=True)

    def top(self, n: int = 5) -> list[Candidate]:
        return list(self.candidates[:n])
