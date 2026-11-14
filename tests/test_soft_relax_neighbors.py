"""Spatial candidate completeness against the frozen exhaustive pair search."""

from itertools import product
import math
import random
from unittest.mock import patch

import gemmi
import pytest

from cofkit import soft_relax
from cofkit.soft_relax import SoftRelaxConfig, _System
from soft_relax_reference import exhaustive_pair_list


def _system(cell, fractional, symbols=None):
    return _System(
        labels=[f"a{i}" for i in range(len(fractional))],
        symbols=list(symbols or ["C"] * len(fractional)),
        frac=[list(f) for f in fractional],
        cell=cell,
        orth=cell.orth.mat.tolist(),
        bonds=[], angles13=[], excluded=set(), excluded13=set(),
        excluded14=set(), bonded_label_images=set(),
    )


def _check(system, config):
    expected = exhaustive_pair_list(system, config)
    actual = soft_relax._build_pair_list(system, config)
    assert actual == expected
    cart = [soft_relax._frac_to_cart(system.orth, f) for f in system.frac]
    # Identical ordering should yield bit-for-bit equality, not merely close
    # forces that could hide an omitted weak contact or a changed image tie.
    expected_energy, expected_forces = soft_relax._energy_and_forces(
        cart, [], expected, config.repulsion_coefficients[-1], config
    )
    energy, forces = soft_relax._energy_and_forces(
        cart, [], actual, config.repulsion_coefficients[-1], config
    )
    assert energy == expected_energy
    assert forces == expected_forces
    return actual, energy, forces


@pytest.mark.parametrize("cell_parameters", [
    (40, 50, 60, 90, 90, 90),
    (12, 14, 15, 70, 80, 45),
    (3.5, 4, 6, 90, 90, 90),
    (45, 8, 7, 90, 90, 110),
    (30, 30, 18, 90, 90, 20),
])
@pytest.mark.parametrize("config", [
    SoftRelaxConfig(),
    SoftRelaxConfig(hbond_guard=False),
    SoftRelaxConfig(repel_hydrogen=False),
    SoftRelaxConfig(neighbor_margin=0),
    SoftRelaxConfig(repulsion_scale=0.65, neighbor_margin=0.4, exclude_one_four=False),
])
def test_randomized_pair_energy_force_equivalence(cell_parameters, config):
    rng = random.Random(18937)
    cell = gemmi.UnitCell(*cell_parameters)
    for _ in range(3):
        fractional = [[rng.uniform(-2.5, 3.5) for _ in range(3)] for _ in range(18)]
        fractional[0] = [0.2, 0.2, 0.2]
        fractional[1] = [0.2 + 0.5 / cell.a, 0.2, 0.2]
        fractional[2] = fractional[0][:]
        symbols = ["C", "C", "N"] + [rng.choice(("C", "N", "O", "H", "B", "Xe")) for _ in range(15)]
        system = _system(cell, fractional, symbols)
        _, energy, forces = _check(system, config)
        assert energy > 0
        assert any(any(component for component in force) for force in forces)


@pytest.mark.parametrize("n", [0, 1])
@pytest.mark.parametrize("length", [0.7, 40])
def test_empty_and_single_atom_periodic_self_images(n, length):
    system = _system(gemmi.UnitCell(length, 20, 20, 90, 90, 90), [(0.1, 0.2, 0.3)] * n)
    pairs, _, _ = _check(system, SoftRelaxConfig())
    if n and length < 1:
        assert pairs and pairs[0][:2] == (0, 0)
        # The two equally near images keep the old lexicographic tie break.
        assert pairs[0][2] == [-length, 0, 0]


def test_cutoff_and_bin_boundaries_are_conservative():
    config = SoftRelaxConfig(neighbor_margin=0)
    cutoff = 2 * soft_relax._element_radius("C", config) * config.repulsion_scale
    length = 6 * cutoff
    cell = gemmi.UnitCell(length, length, length, 90, 90, 90)
    positions = [0, math.nextafter(0, -math.inf), 1]
    for boundary in (cutoff / length, 1 / 5, 2 / 5, 1 - cutoff / length):
        positions.extend((math.nextafter(boundary, -math.inf), boundary, math.nextafter(boundary, math.inf)))
    system = _system(cell, [(p, 0.3, 0.4) for p in positions])
    _check(system, config)


def test_individual_periodic_translations_and_atom_reordering():
    cell = gemmi.UnitCell(19, 23, 27, 67, 82, 41)
    fractional = [(0.01, 0.2, 0.4), (0.99, 0.2, 0.4), (0.5, 0.8, 0.7), (0.51, 0.8, 0.7)]
    for order in ((0, 1, 2, 3), (3, 1, 0, 2)):
        for translated in (False, True):
            positions = [
                [fractional[i][k] + ((i + 1) * (-1) ** k if translated else 0) for k in range(3)]
                for i in order
            ]
            _check(_system(cell, positions), SoftRelaxConfig())


@pytest.mark.parametrize("translation", [1e6, 1e9, 1e12])
def test_unwrapped_large_images_do_not_lose_cutoff_boundary_pairs(translation):
    config = SoftRelaxConfig(neighbor_margin=0)
    cutoff = 2 * soft_relax._element_radius("C", config) * config.repulsion_scale
    length = 10 * cutoff
    for _ in range(8):
        length = math.nextafter(length, math.inf)
    # Subtracting an unwrapped endpoint can round the reference separation
    # inward even when separately wrapped bin coordinates lie just farther
    # apart. Candidate padding must cover that cancellation, not only a ULP
    # of the Cartesian cutoff itself.
    fractional = [(translation + 0.1, 0.3, 0.4), (0.2, 0.3, 0.4)]
    fractional.extend((0.6, 0.6, i / 20) for i in range(10))
    _check(_system(gemmi.UnitCell(length, length, length, 90, 90, 90), fractional), config)


@pytest.mark.parametrize("cell_parameters", [
    (40, 50, 60, 90, 90, 90), (19, 23, 27, 67, 82, 41),
    (24, 24, 20, 90, 90, 5),
])
def test_randomized_large_images_and_sharp_skew_cells(cell_parameters):
    rng = random.Random(81267)
    cell = gemmi.UnitCell(*cell_parameters)
    for _ in range(60):
        fractional = [
            [rng.random() + rng.choice((0, 0, 0, -3, 1e6, -1e12)) for _ in range(3)]
            for _ in range(12)
        ]
        config = SoftRelaxConfig(neighbor_margin=rng.choice((0, 0.2, 2)))
        _check(_system(cell, fractional, ["C", "N", "O", "H"] * 3), config)


@pytest.mark.parametrize("exclusion", ["excluded", "excluded13", "excluded14"])
@pytest.mark.parametrize("distance", [0.6, 1.8])
def test_periodic_exclusions_keep_image_and_deep_overlap_rules(exclusion, distance):
    system = _system(
        gemmi.UnitCell(20, 20, 20, 90, 90, 90),
        [(0.99, 0.5, 0.5), (0.99 + distance / 20 - 1, 0.5, 0.5)],
    )
    # A different image of the same labels must not exclude this contact.
    getattr(system, exclusion).add((0, 1, (0, 0, 0)))
    pairs, _, _ = _check(system, SoftRelaxConfig())
    assert len(pairs) == 1
    getattr(system, exclusion).add((0, 1, (1, 0, 0)))
    pairs, _, _ = _check(system, SoftRelaxConfig())
    assert bool(pairs) == (exclusion != "excluded" and distance < 1.2)


@pytest.mark.parametrize("guard", [False, True])
@pytest.mark.parametrize("repel_hydrogen", [False, True])
def test_hydrogen_guard_and_radius_modes(guard, repel_hydrogen):
    system = _system(
        gemmi.UnitCell(6, 15, 15, 90, 90, 90),
        [(1.2 / 6, 0.5, 0.5), (3 / 6, 0.5, 0.5)], ["H", "N"],
    )
    config = SoftRelaxConfig(hbond_guard=guard, repel_hydrogen=repel_hydrogen)
    pairs, _, _ = _check(system, config)
    assert bool(pairs) == (not guard)


def test_sparse_pair_search_avoids_quadratic_exact_image_calls():
    positions = [
        ((20 * x + offset) / 80, 20 * y / 80, 20 * z / 80)
        for x, y, z in product(range(4), repeat=3) for offset in (0, 0.8)
    ]
    system = _system(gemmi.UnitCell(80, 80, 80, 90, 90, 90), positions)
    config = SoftRelaxConfig()
    expected = exhaustive_pair_list(system, config)
    with patch.object(soft_relax, "images_within", wraps=soft_relax.images_within) as images:
        pairs = soft_relax._build_pair_list(system, config)
    assert pairs == expected
    assert len(pairs) == 64
    assert images.call_count <= 4 * len(positions)
    assert images.call_count < len(positions) * (len(positions) + 1) // 20


def test_nonpositive_cutoff_preserves_exhaustive_behavior():
    system = _system(gemmi.UnitCell(20, 20, 20, 90, 90, 90), [(0, 0, 0), (0.2, 0, 0)])
    config = SoftRelaxConfig(repulsion_scale=0, neighbor_margin=0)
    _check(system, config)
