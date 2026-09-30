from __future__ import annotations

import math
from typing import Callable, Mapping, Sequence

from ..geometry import Frame, centroid, normalize, vec3
from ..model import MonomerSpec, ReactiveMotif
from .motif_registry import MotifKindDefinition, MotifKindRegistry, default_motif_kind_registry
from .molecule import Molecule


MotifDetectionHandler = Callable[[Molecule, str, MotifKindDefinition], MonomerSpec]


# Cited: Cordero et al., "Covalent radii revisited", Dalton Trans., 2008,
# 2832-2838 (DOI 10.1039/B801115J), single-bond covalent radii in angstrom
# (C is the sp3 value). Single shared owner for both bond-detection call sites
# (build-side motif detection here, decompose-side in cofkit.decompose) — one
# table, one value per element (magic-number plan W1.4/W3.1).
COVALENT_RADII_SINGLE_BOND_ANGSTROM: dict[int, float] = {
    1: 0.31,  # H
    5: 0.84,  # B
    6: 0.76,  # C
    7: 0.71,  # N
    8: 0.66,  # O
    9: 0.57,  # F
    15: 1.07,  # P
    16: 1.05,  # S
    17: 1.02,  # Cl
}

# Fixed atomic-number facts used to derive the symbol-keyed view of the shared
# table above; not tunables.
_SYMBOL_TO_ATOMIC_NUMBER = {
    "H": 1,
    "B": 5,
    "C": 6,
    "N": 7,
    "O": 8,
    "F": 9,
    "P": 15,
    "S": 16,
    "Cl": 17,
}
COVALENT_RADII_BY_SYMBOL: dict[str, float] = {
    symbol: COVALENT_RADII_SINGLE_BOND_ANGSTROM[z]
    for symbol, z in _SYMBOL_TO_ATOMIC_NUMBER.items()
}

# Heuristic — pending calibration: build-side fallback covalent radius (A) for
# elements missing from COVALENT_RADII_SINGLE_BOND_ANGSTROM. Tighter than the
# decompose-side fallback (cofkit.decompose.DECOMPOSE_COVALENT_RADIUS_FALLBACK_ANGSTROM)
# because the build side detects motifs on monomers with known chemistry.
# Per-context fallback values are an explicit owner decision (2026-09-28,
# MAGIC_NUMBER_FIX_PLAN.md W3.1): decompose may face unknown input structures,
# not only COFKit-built files, so the two sides intentionally differ.
BUILD_SIDE_COVALENT_RADIUS_FALLBACK_ANGSTROM = 0.7

# Heuristic — pending calibration: build-side bond-detection tolerance, a
# multiplicative slack on the covalent-radius sum. Build-side monomers have
# known chemistry, so one uniform factor suffices. The split from decompose's
# additive tolerance model is an explicit owner decision (2026-09-28,
# MAGIC_NUMBER_FIX_PLAN.md W3.1): decompose keeps its own, more generous
# additive model because it faces unknown input structures.
BUILD_SIDE_BOND_TOLERANCE_FACTOR = 1.3


class MotifDetector:
    """Discovers reactive motifs on a MonomerSpec based on geometry heuristics.

    This is a lightweight, dependency-free fallback detector, not the
    practical ingestion path: the CLI, batch, and library routes detect
    motifs through the RDKit/SMARTS handlers in ``cofkit.chem.rdkit``. Known
    limits (impact-review claim T4-12):

    - Built-in handlers cover only ``"amine"`` and ``"aldehyde"``; every
      other motif kind raises ``ValueError("unsupported motif kind")`` unless
      a custom handler is registered via ``register``.
    - Bonding is perceived purely from interatomic distances
      (``COVALENT_RADII_BY_SYMBOL`` sums scaled by
      ``BUILD_SIDE_BOND_TOLERANCE_FACTOR``): no bond orders, no aromaticity,
      no periodic images, and no hydrogen suppression beyond the explicit
      neighbor counting in the handlers.
    - The returned ``MonomerSpec`` carries ``molecule.bonds`` verbatim, and
      ``Molecule.from_xyz`` never populates bonds, so detector-built specs
      from XYZ input have an empty bond graph. Downstream atomistic
      realization and bond-aware consumers cannot use such specs directly.
    """

    def __init__(
        self,
        *,
        registry: MotifKindRegistry | None = None,
        handlers: Mapping[str, MotifDetectionHandler] | None = None,
    ) -> None:
        self.registry = registry or default_motif_kind_registry()
        self._handlers: dict[str, MotifDetectionHandler] = dict(handlers or {})

    def register(self, motif_kind: str, handler: MotifDetectionHandler) -> None:
        self._handlers[motif_kind] = handler

    def detect(self, molecule: Molecule, monomer_id: str, motif_kind: str) -> MonomerSpec:
        definition = self.registry.get(motif_kind)
        try:
            handler = self._handlers[motif_kind]
        except KeyError as exc:
            raise ValueError(f"unsupported motif kind {motif_kind!r}") from exc
        return handler(molecule, monomer_id, definition)

    @classmethod
    def builtin(cls, *, registry: MotifKindRegistry | None = None) -> "MotifDetector":
        detector = cls(registry=registry)
        detector.register("amine", _detect_amine)
        detector.register("aldehyde", _detect_aldehyde)
        return detector

    @staticmethod
    def detect_amine(molecule: Molecule, monomer_id: str) -> MonomerSpec:
        """Finds primary amines (-NH2)."""
        return MotifDetector.builtin().detect(molecule, monomer_id, "amine")

    @staticmethod
    def detect_aldehyde(molecule: Molecule, monomer_id: str) -> MonomerSpec:
        """Finds aldehydes (-CHO)."""
        return MotifDetector.builtin().detect(molecule, monomer_id, "aldehyde")

    @staticmethod
    def _find_neighbors(atom_idx: int, molecule: Molecule) -> list[int]:
        # Simple distance-based neighbor detection (very basic)
        neighbors = []

        p1 = molecule.positions[atom_idx]
        r1 = COVALENT_RADII_BY_SYMBOL.get(molecule.symbols[atom_idx], BUILD_SIDE_COVALENT_RADIUS_FALLBACK_ANGSTROM)

        for j, p2 in enumerate(molecule.positions):
            if j == atom_idx:
                continue
            r2 = COVALENT_RADII_BY_SYMBOL.get(molecule.symbols[j], BUILD_SIDE_COVALENT_RADIUS_FALLBACK_ANGSTROM)
            dist = math.sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2 + (p1[2]-p2[2])**2)

            if dist < (r1 + r2) * BUILD_SIDE_BOND_TOLERANCE_FACTOR:
                neighbors.append(j)

        return neighbors


def _detect_amine(molecule: Molecule, monomer_id: str, definition: MotifKindDefinition) -> MonomerSpec:
    motifs = []
    for i, sym in enumerate(molecule.symbols):
        if sym != "N":
            continue

        neighbors = MotifDetector._find_neighbors(i, molecule)
        h_neighbors = [n for n in neighbors if molecule.symbols[n] == "H"]
        heavy_neighbors = [n for n in neighbors if molecule.symbols[n] != "H"]

        if len(h_neighbors) >= 2 and len(heavy_neighbors) == 1:
            anchor = heavy_neighbors[0]
            origin = molecule.positions[i]
            anchor_pos = molecule.positions[anchor]
            h1_pos = molecule.positions[h_neighbors[0]]
            h2_pos = molecule.positions[h_neighbors[1]]

            v_anchor_to_n = vec3(
                origin[0] - anchor_pos[0],
                origin[1] - anchor_pos[1],
                origin[2] - anchor_pos[2],
            )
            primary = normalize(v_anchor_to_n)

            v1 = vec3(h1_pos[0] - origin[0], h1_pos[1] - origin[1], h1_pos[2] - origin[2])
            v2 = vec3(h2_pos[0] - origin[0], h2_pos[1] - origin[1], h2_pos[2] - origin[2])
            nx = v1[1] * v2[2] - v1[2] * v2[1]
            ny = v1[2] * v2[0] - v1[0] * v2[2]
            nz = v1[0] * v2[1] - v1[1] * v2[0]
            normal = normalize((nx, ny, nz))

            motifs.append(
                ReactiveMotif(
                    id=f"{definition.kind}_{i}",
                    kind=definition.kind,
                    atom_ids=(i, *h_neighbors),
                    frame=Frame(origin=origin, primary=primary, normal=normal),
                    allowed_reaction_templates=definition.allowed_reaction_templates,
                    metadata={"reactive_atom_id": i, "anchor_atom_id": anchor},
                )
            )

    return MonomerSpec(
        id=monomer_id,
        name=monomer_id,
        motifs=tuple(motifs),
        atom_symbols=molecule.symbols,
        atom_positions=molecule.positions,
        bonds=molecule.bonds,
    )


def _detect_aldehyde(molecule: Molecule, monomer_id: str, definition: MotifKindDefinition) -> MonomerSpec:
    motifs = []
    for i, sym in enumerate(molecule.symbols):
        if sym != "C":
            continue

        neighbors = MotifDetector._find_neighbors(i, molecule)
        o_neighbors = [n for n in neighbors if molecule.symbols[n] == "O"]
        h_neighbors = [n for n in neighbors if molecule.symbols[n] == "H"]
        r_neighbors = [n for n in neighbors if molecule.symbols[n] not in ("O", "H")]

        if len(o_neighbors) == 1 and len(h_neighbors) == 1 and len(r_neighbors) == 1:
            anchor = r_neighbors[0]
            origin = molecule.positions[i]
            anchor_pos = molecule.positions[anchor]
            o_pos = molecule.positions[o_neighbors[0]]

            v_anchor_to_c = vec3(
                origin[0] - anchor_pos[0],
                origin[1] - anchor_pos[1],
                origin[2] - anchor_pos[2],
            )
            primary = normalize(v_anchor_to_c)

            v_c_to_o = vec3(
                o_pos[0] - origin[0],
                o_pos[1] - origin[1],
                o_pos[2] - origin[2],
            )

            nx = primary[1] * v_c_to_o[2] - primary[2] * v_c_to_o[1]
            ny = primary[2] * v_c_to_o[0] - primary[0] * v_c_to_o[2]
            nz = primary[0] * v_c_to_o[1] - primary[1] * v_c_to_o[0]

            n_len = math.sqrt(nx * nx + ny * ny + nz * nz)
            if n_len > 1e-4:
                normal = (nx / n_len, ny / n_len, nz / n_len)
            else:
                normal = (0.0, 0.0, 1.0)

            motifs.append(
                ReactiveMotif(
                    id=f"{definition.kind}_{i}",
                    kind=definition.kind,
                    atom_ids=(i, o_neighbors[0], h_neighbors[0]),
                    frame=Frame(origin=origin, primary=primary, normal=normal),
                    allowed_reaction_templates=definition.allowed_reaction_templates,
                    metadata={"reactive_atom_id": i, "anchor_atom_id": anchor},
                )
            )

    return MonomerSpec(
        id=monomer_id,
        name=monomer_id,
        motifs=tuple(motifs),
        atom_symbols=molecule.symbols,
        atom_positions=molecule.positions,
        bonds=molecule.bonds,
    )
