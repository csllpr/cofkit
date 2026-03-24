from __future__ import annotations

from dataclasses import dataclass, field
from typing import Mapping


@dataclass(frozen=True)
class MotifKindDefinition:
    kind: str
    id_prefix: str
    cif_symbol: str = "C"
    allowed_reaction_templates: tuple[str, ...] = ()
    rdkit_smarts: str | None = None
    metadata: Mapping[str, object] = field(default_factory=dict)


class MotifKindRegistry:
    def __init__(self, definitions: tuple[MotifKindDefinition, ...] = ()) -> None:
        self._definitions: dict[str, MotifKindDefinition] = {}
        for definition in definitions:
            self.register(definition)

    def register(self, definition: MotifKindDefinition) -> None:
        self._definitions[definition.kind] = definition

    def get(self, kind: str) -> MotifKindDefinition:
        try:
            return self._definitions[kind]
        except KeyError as exc:
            raise KeyError(f"unknown motif kind {kind!r}") from exc

    def supported_kinds(self) -> tuple[str, ...]:
        return tuple(sorted(self._definitions))

    @classmethod
    def builtin(cls) -> "MotifKindRegistry":
        return cls(_builtin_motif_kind_definitions())


def default_motif_kind_registry() -> MotifKindRegistry:
    return MotifKindRegistry.builtin()


def motif_pseudo_atom_symbol(kind: str, *, registry: MotifKindRegistry | None = None) -> str:
    lookup = registry or default_motif_kind_registry()
    try:
        return lookup.get(kind).cif_symbol
    except KeyError:
        return "C"


def _builtin_motif_kind_definitions() -> tuple[MotifKindDefinition, ...]:
    return (
        MotifKindDefinition(
            kind="amine",
            id_prefix="ami",
            cif_symbol="N",
            allowed_reaction_templates=("imine_bridge", "keto_enamine_bridge"),
            rdkit_smarts="[NX3;H2;!$(N[C,S]=O)]-[#6]",
        ),
        MotifKindDefinition(
            kind="aldehyde",
            id_prefix="ald",
            cif_symbol="C",
            allowed_reaction_templates=("imine_bridge", "hydrazone_bridge", "vinylene_bridge"),
            rdkit_smarts="[CX3H1](=[OX1])-[#6]",
        ),
        MotifKindDefinition(
            kind="hydrazine",
            id_prefix="hyd",
            cif_symbol="N",
        ),
        MotifKindDefinition(
            kind="hydrazide",
            id_prefix="hdz",
            cif_symbol="N",
            allowed_reaction_templates=("hydrazone_bridge",),
            rdkit_smarts="[NX3H2][NX3H][CX3](=[OX1])[#6]",
        ),
        MotifKindDefinition(
            kind="boronic_acid",
            id_prefix="bor",
            cif_symbol="B",
            allowed_reaction_templates=("boronate_ester_bridge", "boroxine_trimerization"),
            rdkit_smarts="[#6]-[BX3]([OX2H])[OX2H]",
        ),
        MotifKindDefinition(
            kind="catechol",
            id_prefix="cat",
            cif_symbol="O",
            allowed_reaction_templates=("boronate_ester_bridge",),
            rdkit_smarts="[OX2H]-[c;r6]",
        ),
        MotifKindDefinition(
            kind="keto_aldehyde",
            id_prefix="kal",
            cif_symbol="C",
            allowed_reaction_templates=("keto_enamine_bridge",),
            rdkit_smarts="[CX3H1](=[OX1])[c]",
        ),
        MotifKindDefinition(
            kind="activated_methylene",
            id_prefix="act",
            cif_symbol="C",
            allowed_reaction_templates=("vinylene_bridge",),
            rdkit_smarts="[CH2,CH3]",
        ),
        MotifKindDefinition(
            kind="nitrile",
            id_prefix="nit",
            cif_symbol="N",
            allowed_reaction_templates=("triazine_trimerization",),
        ),
    )
