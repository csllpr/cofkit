"""Distance bonding retains the exhaustive algorithm's exact output."""

from dataclasses import replace
from itertools import product
import math
import random
from unittest.mock import patch

import gemmi
import pytest

from cofkit import decompose
from cofkit.cofid import parse_cofid
from cofkit.decompose_cif import PeriodicCifAtoms, _basis_from_gemmi_cell, read_periodic_cif_atoms
from decompose_reference import exhaustive_bond_candidates
from soft_relax_cof_cases import COF_CASES, replicate_cif, write_cof_seed


def atoms_for(cell, positions, symbols=None):
    return PeriodicCifAtoms(
        symbols=tuple(symbols or ["C"] * len(positions)),
        fractional_positions=tuple(tuple(p) for p in positions),
        cartesian_positions=tuple(tuple(cell.orthogonalize(gemmi.Fractional(*p))) for p in positions),
        cell_basis=_basis_from_gemmi_cell(cell),
        info={"_atom_site_label": tuple(f"site{i}" for i in range(len(positions)))},
    )


@pytest.mark.parametrize("parameters", [
    (40, 50, 60, 90, 90, 90), (12, 14, 15, 70, 80, 45),
    (2, 3, 4, 90, 90, 90), (45, 6, 7, 90, 90, 110),
    (24, 24, 20, 90, 90, 5),
])
def test_random_cells_and_unwrapped_images_match_exhaustive(parameters):
    rng = random.Random(72134)
    cell = gemmi.UnitCell(*parameters)
    for iteration in range(20):
        positions = [
            [rng.random() + rng.choice((0, 0, -3, 1e6, -1e12)) for _ in range(3)]
            for _ in range(16)
        ]
        positions[1] = positions[0][:]  # coincident sites retain the short-distance floor
        symbols = [rng.choice(("C", "N", "O", "B", "H", "S", "Xe", "La")) for _ in positions]
        atoms = atoms_for(cell, positions, symbols)
        assert decompose._distance_inferred_bond_candidates(atoms) == exhaustive_bond_candidates(atoms), iteration


@pytest.mark.parametrize("translation", [0, 1e6, 1e9, 1e12])
def test_exact_bond_window_and_bin_boundaries(translation):
    cutoff = 2 * decompose._covalent_radius(6) + decompose.DECOMPOSE_BOND_TOLERANCE_HEAVY_ANGSTROM
    length = 10 * cutoff
    for _ in range(8):
        length = math.nextafter(length, math.inf)
    fractions = [0, 0.1, 0.2, 0.3, 0.9]
    positions = [(translation + .1, .3, .4)]
    for x in fractions:
        positions.extend((f, .3, .4) for f in (math.nextafter(x, -math.inf), x, math.nextafter(x, math.inf)))
    atoms = atoms_for(gemmi.UnitCell(length, length, length, 90, 90, 90), positions)
    assert decompose._distance_inferred_bond_candidates(atoms) == exhaustive_bond_candidates(atoms)


@pytest.mark.parametrize("n", [0, 1])
def test_empty_or_single_atom_has_no_inferred_self_bond(n):
    atoms = atoms_for(gemmi.UnitCell(1, 1, 1, 90, 90, 90), [(0, 0, 0)] * n)
    assert decompose._distance_inferred_bond_candidates(atoms) == exhaustive_bond_candidates(atoms) == ()


def test_sparse_system_avoids_exhaustive_exact_distance_calls():
    positions = [
        ((20 * x + d) / 80, 20 * y / 80, 20 * z / 80)
        for x, y, z in product(range(4), repeat=3) for d in (0, .8)
    ]
    atoms = atoms_for(gemmi.UnitCell(80, 80, 80, 90, 90, 90), positions)
    expected = exhaustive_bond_candidates(atoms)
    with patch.object(decompose, "_minimum_image_bond_geometry", wraps=decompose._minimum_image_bond_geometry) as distance:
        actual = decompose._distance_inferred_bond_candidates(atoms)
    assert actual == expected
    assert len(actual) == 64
    assert distance.call_count <= 4 * len(atoms)
    assert distance.call_count < len(atoms) * (len(atoms) - 1) // 20


@pytest.mark.parametrize("case", COF_CASES, ids=lambda case: case["id"])
def test_common_cof_full_distance_decomposition_matches_exhaustive(tmp_path, case):
    path = tmp_path / "seed.cif"
    write_cof_seed(case, path)
    atoms = read_periodic_cif_atoms(path)
    assert decompose._distance_inferred_bond_candidates(atoms) == exhaustive_bond_candidates(atoms)
    topology = parse_cofid(case["cofid"]).topology
    actual = decompose.decompose_cif_to_cofid(path, topology=topology, bond_mode="distance")
    with patch.object(decompose, "_distance_inferred_bond_candidates", exhaustive_bond_candidates):
        expected = decompose.decompose_cif_to_cofid(path, topology=topology, bond_mode="distance")
    # Includes verdict, precursors, graph gains, ledger and attachment recipe.
    assert actual.to_dict() == expected.to_dict()


def test_cof_supercell_pair_order_and_pruning_match_exhaustive(tmp_path):
    seed, supercell = tmp_path / "seed.cif", tmp_path / "supercell.cif"
    write_cof_seed(next(case for case in COF_CASES if case["id"] == "cof-5"), seed)
    replicate_cif(seed, supercell, (2, 2, 1))
    atoms = read_periodic_cif_atoms(supercell)
    assert decompose._distance_inferred_bond_candidates(atoms) == exhaustive_bond_candidates(atoms)
    order = list(range(len(atoms)))
    random.Random(378).shuffle(order)
    permuted = replace(atoms, symbols=tuple(atoms.symbols[i] for i in order),
                       fractional_positions=tuple(atoms.fractional_positions[i] for i in order))
    assert decompose._distance_inferred_bond_candidates(permuted) == exhaustive_bond_candidates(permuted)
