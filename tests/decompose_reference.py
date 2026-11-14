"""Frozen exhaustive distance-candidate loop, before spatial acceleration.

Image selection, distance rules and pruning are shared unchanged policy.
"""

from cofkit.decompose import (
    BondCandidate, _fractional_from_cartesian_row_norms,
    _minimum_image_bond_geometry, _is_plausible_bond_distance,
    _prune_distance_inferred_bond_candidates,
)
from cofkit.decompose_cif import PeriodicCifAtoms


def exhaustive_bond_candidates(atoms: PeriodicCifAtoms) -> tuple[BondCandidate, ...]:
    candidates: list[BondCandidate] = []
    reciprocal_row_norms = _fractional_from_cartesian_row_norms(atoms.cell_basis)
    for idx_1, symbol_1 in enumerate(atoms.symbols):
        for idx_2 in range(idx_1 + 1, len(atoms.symbols)):
            symbol_2 = atoms.symbols[idx_2]
            image, distance = _minimum_image_bond_geometry(
                atoms,
                idx_1,
                idx_2,
                reciprocal_row_norms=reciprocal_row_norms,
            )
            if _is_plausible_bond_distance(symbol_1, symbol_2, distance):
                candidates.append(BondCandidate(idx_1, idx_2, distance, explicit_order=None, periodic_image=image))
    return _prune_distance_inferred_bond_candidates(atoms, tuple(candidates))

