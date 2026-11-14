"""CIF evidence must not be silently changed before decomposition."""

from itertools import permutations
from dataclasses import replace
from math import sqrt

import pytest

from cofkit.decompose import _build_bonded_mol, _explicit_bond_candidates, decompose_cif_to_cofid
from cofkit.decompose_cif import CifInputError, PeriodicCifAtoms


def two_atom_input(**columns):
    return PeriodicCifAtoms(
        symbols=("C", "N"), fractional_positions=((0, 0, 0), (0.13, 0, 0)),
        cartesian_positions=((0, 0, 0), (1.3, 0, 0)),
        cell_basis=((10, 0, 0), (0, 10, 0), (0, 0, 10)),
        info={"_atom_site_label": ("C1", "N1"), **columns},
    )


@pytest.mark.parametrize("orders", [("S", "D"), ("D", "S")])
@pytest.mark.parametrize("reverse", [False, True])
def test_duplicate_edge_conflicting_orders_rejected(orders, reverse):
    atoms = two_atom_input(_ccdc_geom_bond_type=orders)
    a, b = (("C1", "N1"), ("N1", "C1")) if reverse else (("C1", "C1"), ("N1", "N1"))
    with pytest.raises(CifInputError, match="conflicting bond orders"):
        _explicit_bond_candidates(atoms, {"C1": 0, "N1": 1}, a, b)


def test_duplicate_missing_order_cannot_override_explicit_evidence():
    for rows in permutations((("?", "1.31"), ("D", "1.30"), ("D", "1.32"))):
        orders, distances = zip(*rows)
        atoms = two_atom_input(_ccdc_geom_bond_type=orders, _geom_bond_distance=distances)
        candidates, missing = _explicit_bond_candidates(
            atoms, {"C1": 0, "N1": 1}, ("C1",) * 3, ("N1",) * 3,
        )
        assert len(candidates) == 1
        assert candidates[0].explicit_order == 2.0
        assert candidates[0].distance == 1.30
        assert missing == 0


def test_duplicate_site_labels_rejected_before_graph_building():
    with pytest.raises(CifInputError, match="labels must be unique"):
        _build_bonded_mol(two_atom_input(_atom_site_label=("C1", "C1")))


def test_unknown_bond_endpoint_is_not_silently_discarded():
    with pytest.raises(CifInputError, match="unknown atom-site label"):
        _explicit_bond_candidates(two_atom_input(), {"C1": 0, "N1": 1}, ("missing",), ("N1",))


def test_parallel_images_are_kept_as_distinct_edges():
    atoms = two_atom_input(
        _ccdc_geom_bond_type=("S", "S"),
        _geom_bond_site_symmetry_2=("1_555", "1_655"),
    )
    candidates, _ = _explicit_bond_candidates(atoms, {"C1": 0, "N1": 1}, ("C1",) * 2, ("N1",) * 2)
    assert [c.periodic_image for c in candidates] == [(0, 0, 0), (1, 0, 0)]


@pytest.mark.parametrize("orders", [("S", "D"), ("D", "S")])
def test_public_result_distinguishes_invalid_input_from_internal_error(tmp_path, orders):
    path = tmp_path / "conflict.cif"
    path.write_text("""data_conflict
_cell_length_a 10
_cell_length_b 10
_cell_length_c 10
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
C1 C 0 0 0
N1 N .13 0 0
loop_
_geom_bond_atom_site_label_1
_geom_bond_atom_site_label_2
_ccdc_geom_bond_type
""" + "".join(f"C1 N1 {order}\n" for order in orders))
    result = decompose_cif_to_cofid(path)
    assert result.status == "error"
    assert result.cofid is None
    assert result.metadata["event_status"] == "FAILED_INPUT_VALIDATION"
    assert "conflicting bond orders" in result.reason


@pytest.mark.parametrize("reported", [None, "?", ".", "invalid", "nan", "inf", "0", "-1"])
@pytest.mark.parametrize("reverse", [False, True])
def test_missing_distance_uses_declared_periodic_edge(reported, reverse):
    columns = {
        "_geom_bond_site_symmetry_1": ("1_655" if reverse else "1_555",),
        "_geom_bond_site_symmetry_2": ("1_555" if reverse else "1_655",),
    }
    if reported is not None:
        columns["_geom_bond_distance"] = (reported,)
    atoms = two_atom_input(**columns)
    a, b = ("N1", "C1") if reverse else ("C1", "N1")
    (candidate,), _ = _explicit_bond_candidates(atoms, {"C1": 0, "N1": 1}, (a,), (b,))
    assert candidate.periodic_image == (1, 0, 0)
    assert candidate.distance == pytest.approx(11.3)


def test_missing_distance_accounts_for_both_symmetry_codes_and_triclinic_cell():
    atoms = replace(two_atom_input(
        _geom_bond_site_symmetry_1=("1_455",),
        _geom_bond_site_symmetry_2=("1_565",),
    ), cell_basis=((10, 0, 0), (3, 8, 0), (0, 0, 10)))
    (candidate,), _ = _explicit_bond_candidates(atoms, {"C1": 0, "N1": 1}, ("C1",), ("N1",))
    assert candidate.periodic_image == (1, 1, 0)
    assert candidate.distance == pytest.approx(sqrt(14.3 ** 2 + 8 ** 2))


def test_usable_reported_distance_retains_existing_evidence_policy():
    atoms = two_atom_input(_geom_bond_distance=("1.31(2)",), _geom_bond_site_symmetry_2=("1_655",))
    (candidate,), _ = _explicit_bond_candidates(atoms, {"C1": 0, "N1": 1}, ("C1",), ("N1",))
    assert candidate.distance == 1.31


def test_missing_distance_without_image_hint_keeps_nearest_image_selection():
    atoms = replace(two_atom_input(), fractional_positions=((0, 0, 0), (.87, 0, 0)))
    (candidate,), _ = _explicit_bond_candidates(atoms, {"C1": 0, "N1": 1}, ("C1",), ("N1",))
    assert candidate.periodic_image == (-1, 0, 0)
    assert candidate.distance == pytest.approx(1.3)


def test_bond_order_inference_uses_length_of_declared_edge():
    atoms = two_atom_input(
        _geom_bond_atom_site_label_1=("C1",),
        _geom_bond_atom_site_label_2=("N1",),
        _geom_bond_site_symmetry_2=("1_655",),
    )
    build = _build_bonded_mol(atoms)
    assert build.candidates[0].distance == pytest.approx(11.3)
    assert build.mol.GetBondWithIdx(0).GetBondTypeAsDouble() == 1.0
