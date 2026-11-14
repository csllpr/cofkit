"""End-to-end spatial/exhaustive equivalence on common COF precursor models."""

from collections import Counter
from math import prod
from unittest.mock import patch

import gemmi
import pytest

from cofkit import soft_relax
from soft_relax_cof_cases import (
    COF_CASES, CROWDED_LAYER_CASE_IDS, normalized_report,
    replicate_cif, write_cof_seed, write_crowded_layers,
)
from soft_relax_reference import exhaustive_pair_list


def _compare_relaxations(source, tmp_path):
    config = soft_relax.SoftRelaxConfig()
    block = gemmi.cif.read_file(str(source)).sole_block()
    system, _ = soft_relax._parse_system(block, config)
    reference = exhaustive_pair_list(system, config)
    assert soft_relax._build_pair_list(system, config) == reference
    cart = [soft_relax._frac_to_cart(system.orth, f) for f in system.frac]
    assert soft_relax._energy_and_forces(cart, [], reference, config.repulsion_coefficients[-1], config) == soft_relax._energy_and_forces(
        cart, [], soft_relax._build_pair_list(system, config), config.repulsion_coefficients[-1], config
    )
    with patch.object(soft_relax, "_build_pair_list", exhaustive_pair_list):
        baseline = soft_relax.relax_cif_clashes(source, tmp_path / "exhaustive.cif", config)
    current = soft_relax.relax_cif_clashes(source, tmp_path / "spatial.cif", config)
    assert normalized_report(current) == normalized_report(baseline)
    assert (tmp_path / "spatial.cif").read_bytes() == (tmp_path / "exhaustive.cif").read_bytes()
    relaxed = gemmi.cif.read_file(current.output_path).sole_block()
    assert len(relaxed.find_values("_atom_site_label")) == len(block.find_values("_atom_site_label"))
    assert len(relaxed.find_values("_geom_bond_atom_site_label_1")) == len(block.find_values("_geom_bond_atom_site_label_1"))
    return current


@pytest.mark.parametrize("case", COF_CASES, ids=lambda case: case["id"])
def test_common_cof_matches_exhaustive_search_through_relaxation(case, tmp_path):
    source = tmp_path / "seed.cif"
    write_cof_seed(case, source)
    _compare_relaxations(source, tmp_path)


@pytest.mark.parametrize("case_id", CROWDED_LAYER_CASE_IDS)
def test_crowded_common_cof_layers_match_during_active_relaxation(case_id, tmp_path):
    case = next(case for case in COF_CASES if case["id"] == case_id)
    seed, source = tmp_path / "layer.cif", tmp_path / "crowded.cif"
    write_cof_seed(case, seed)
    write_crowded_layers(seed, source)
    report = _compare_relaxations(source, tmp_path)
    assert report.clashes_before > 0
    assert report.clashes_after < report.clashes_before


@pytest.mark.parametrize("case_id,repeats", [
    ("cof-5", (2, 2, 1)), ("cof-300-single-net", (2, 2, 2)),
])
def test_benchmark_supercells_preserve_bonds_and_precursor_instances(case_id, repeats, tmp_path):
    from cofkit.validation import measure_inter_instance_bond_distances

    case = next(case for case in COF_CASES if case["id"] == case_id)
    source, target = tmp_path / "seed.cif", tmp_path / "supercell.cif"
    write_cof_seed(case, source)
    replicate_cif(source, target, repeats)
    blocks = [gemmi.cif.read_file(str(p)).sole_block() for p in (source, target)]
    systems = [soft_relax._parse_system(block, soft_relax.SoftRelaxConfig())[0] for block in blocks]
    first, second = systems
    assert Counter(second.symbols) == Counter({s: n * prod(repeats) for s, n in Counter(first.symbols).items()})
    assert len(second.bonds) == len(first.bonds) * prod(repeats)
    assert len(measure_inter_instance_bond_distances(blocks[1])) == len(measure_inter_instance_bond_distances(blocks[0])) * prod(repeats)

    def distances(system):
        cart = [soft_relax._frac_to_cart(system.orth, f) for f in system.frac]
        bonds = [((i, j), soft_relax._frac_to_cart(system.orth, shift)) for i, j, shift in system.bonds]
        return soft_relax._bond_distances(cart, bonds)

    assert sorted(distances(second)) == pytest.approx(sorted(distances(first) * prod(repeats)), abs=1e-8)
