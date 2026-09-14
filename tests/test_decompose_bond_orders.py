"""ReDD-COFFEE regressions with independently labeled linkage edges."""

import gzip
import json
from dataclasses import replace
from pathlib import Path
import random
from unittest.mock import patch

import pytest
from rdkit import Chem

from cofkit import decompose as legacy
from cofkit import decompose_events as events
from cofkit.decompose_bond_orders import normalize_imine_bond_orders


FIXTURES = Path(__file__).parent / "fixtures" / "redd_coffee"
SAMPLES = json.loads((FIXTURES / "manifest.json").read_text())


def _signature(mol):
    return (
        tuple((a.GetAtomicNum(), a.GetFormalCharge(), a.GetNumExplicitHs()) for a in mol.GetAtoms()),
        tuple((b.GetBeginAtomIdx(), b.GetEndAtomIdx(), b.GetBondTypeAsDouble()) for b in mol.GetBonds()),
    )


def _valences(mol):
    return tuple(sum(b.GetBondTypeAsDouble() for b in a.GetBonds()) for a in mol.GetAtoms())


def _raw_build(path):
    with patch("cofkit.decompose_bond_orders.normalize_imine_bond_orders", return_value={"changed_bonds": []}):
        return legacy._build_bonded_mol(legacy.read_periodic_cif_atoms(path))


@pytest.fixture(params=SAMPLES, ids=lambda sample: sample["label"])
def sample(request, tmp_path):
    record = request.param
    path = tmp_path / record["source"]
    path.write_bytes(gzip.decompress((FIXTURES / (record["label"] + ".cif.gz")).read_bytes()))
    return record, path


def test_recovers_labeled_linkages_and_binary_precursors(sample):
    record, path = sample
    raw = _raw_build(path)
    repaired = legacy._build_bonded_mol(legacy.read_periodic_cif_atoms(path))
    assert raw.mol.GetNumBonds() == repaired.mol.GetNumBonds() == record["bonds"]
    assert raw.mol.GetNumAtoms() == repaired.mol.GetNumAtoms() == record["atoms"]
    assert _signature(raw.mol)[0] == _signature(repaired.mol)[0]
    assert _valences(raw.mol) == _valences(repaired.mol)
    assert [(c.atom_idx_1, c.atom_idx_2, c.periodic_image) for c in raw.candidates] == [
        (c.atom_idx_1, c.atom_idx_2, c.periodic_image) for c in repaired.candidates
    ]
    report = repaired.metadata["imine_bond_order_normalization"]
    assert report["restored_imine_bonds"] == record["restored"]
    assert report["unresolved_components"] == 0

    detected = events.detect_linkage_events(repaired)
    imines = tuple(event for event in detected.events if event.family == "imine")
    # The oracle is the CHK force-field serial boundary, not the repair's
    # distance rules, matching result, or desired species count.
    assert {tuple(sorted(event.atoms)) for event in imines} == {
        tuple(pair) for pair in record["expected_cuts"]
    }
    result = events._cut_and_reconstruct(repaired, imines, legacy._IMINE_SPEC)
    assert result.status == "ok", result.errors
    assert {m.reactive_group for m in result.monomers} == {"amine", "aldehyde"}
    assert len(result.monomers) == 2
    assert {
        (m.reactive_group, m.connectivity, m.amount, m.canonical_smiles) for m in result.monomers
    } == {
        (m["role"], m["connectivity"], m["amount"], m["smiles"]) for m in record["monomers"]
    }
    for candidate in repaired.candidates:
        assert candidate.explicit_order == repaired.mol.GetBondBetweenAtoms(
            candidate.atom_idx_1, candidate.atom_idx_2,
        ).GetBondTypeAsDouble()

    before = _signature(repaired.mol)
    again = normalize_imine_bond_orders(repaired.mol, {
        frozenset((c.atom_idx_1, c.atom_idx_2)): c.distance for c in repaired.candidates
    })
    assert again["changed_bonds"] == []
    assert before == _signature(repaired.mol)


def test_public_api_retains_repaired_pair(sample):
    record, path = sample
    # Supplying the dataset's known net skips open-ended topology search; a
    # topology validation rejection must still retain the repaired monomers.
    result = legacy.decompose_cif_to_cofid(path, linkage="imine", topology=record["source"].split("_")[0])
    assert len(result.monomers) == 2, result.reason
    assert result.metadata["event_status"] in {events.EVENT_STATUS_COMPLETE, events.EVENT_STATUS_TOPOLOGY}


def test_atom_row_order_does_not_change_precursors(tmp_path):
    record = SAMPLES[0]
    path = tmp_path / "A.cif"
    path.write_bytes(gzip.decompress((FIXTURES / "A.cif.gz").read_bytes()))
    atoms = legacy.read_periodic_cif_atoms(path)
    order = list(range(len(atoms)))
    random.Random(531).shuffle(order)
    permuted = replace(
        atoms,
        symbols=tuple(atoms.symbols[i] for i in order),
        fractional_positions=tuple(atoms.fractional_positions[i] for i in order),
        cartesian_positions=tuple(atoms.cartesian_positions[i] for i in order),
        info={**atoms.info, "_atom_site_label": tuple(atoms.info["_atom_site_label"][i] for i in order)},
    )
    build = legacy._build_bonded_mol(permuted)
    detected = events.detect_linkage_events(build)
    imines = tuple(event for event in detected.events if event.family == "imine")
    result = events._cut_and_reconstruct(build, imines, legacy._IMINE_SPEC)
    assert result.status == "ok", result.errors
    assert {m.canonical_smiles for m in result.monomers} == {m["smiles"] for m in record["monomers"]}


@pytest.mark.parametrize("label", ["A", "B", "C"])
def test_atom_row_order_does_not_change_bond_orders(tmp_path, label):
    """The geometry tie-break makes the repaired assignment structural."""

    record = next(candidate for candidate in SAMPLES if candidate["label"] == label)
    path = tmp_path / record["source"]
    path.write_bytes(gzip.decompress((FIXTURES / f"{label}.cif.gz").read_bytes()))
    atoms = legacy.read_periodic_cif_atoms(path)
    reference = _bond_orders(legacy._build_bonded_mol(atoms).mol)
    for seed in (531, 20260914):
        order = list(range(len(atoms)))
        random.Random(seed).shuffle(order)
        permuted = replace(
            atoms,
            symbols=tuple(atoms.symbols[i] for i in order),
            fractional_positions=tuple(atoms.fractional_positions[i] for i in order),
            cartesian_positions=tuple(atoms.cartesian_positions[i] for i in order),
            info={**atoms.info, "_atom_site_label": tuple(atoms.info["_atom_site_label"][i] for i in order)},
        )
        mapped = {}
        for pair, order_value in _bond_orders(legacy._build_bonded_mol(permuted).mol).items():
            first, second = tuple(pair)
            mapped[frozenset((order[first], order[second]))] = order_value
        assert mapped == reference


def test_small_heteroaromatic_ring_is_not_an_imine_linkage():
    mol = Chem.AddHs(Chem.MolFromSmiles("c1ccncc1"))
    Chem.Kekulize(mol, clearAromaticFlags=True)
    distances = {}
    for bond in mol.GetBonds():
        # Deliberately misleading geometry must not turn a heteroaromatic
        # intramonomer bond into a candidate linkage.
        pair = frozenset((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))
        distances[pair] = 1.44 if bond.GetBondTypeAsDouble() == 2.0 else 1.28
    before = _signature(mol)
    assert normalize_imine_bond_orders(mol, distances)["changed_bonds"] == []
    assert _signature(mol) == before


@pytest.mark.parametrize("guard", ["no_distances", "ambiguous_lengths", "no_explicit_hydrogens", "charged_nitrogens"])
def test_uncertain_chemistry_is_unchanged(guard, tmp_path):
    path = tmp_path / "A.cif"
    path.write_bytes(gzip.decompress((FIXTURES / "A.cif.gz").read_bytes()))
    build = _raw_build(path)
    mol = Chem.Mol(build.mol)
    distances = {frozenset((c.atom_idx_1, c.atom_idx_2)): c.distance for c in build.candidates}
    if guard == "no_distances":
        distances = {}
    elif guard == "ambiguous_lengths":
        distances = {edge: 1.32 for edge in distances}
    elif guard == "no_explicit_hydrogens":
        # Keep indexing fixed but make the H evidence unavailable.
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 1:
                atom.SetAtomicNum(0)
    else:
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 7:
                atom.SetFormalCharge(1)
    before = _signature(mol)
    assert normalize_imine_bond_orders(mol, distances)["changed_bonds"] == []
    assert _signature(mol) == before


def test_failed_perfect_matching_does_not_partially_mutate(tmp_path):
    path = tmp_path / "A.cif"
    path.write_bytes(gzip.decompress((FIXTURES / "A.cif.gz").read_bytes()))
    build = _raw_build(path)
    # Forbid every alternative to the original matching. This simulates an
    # inconsistent/under-specified conjugated component at the solver seam.
    with patch("cofkit.decompose_bond_orders.nx.max_weight_matching", return_value=set()):
        before = _signature(build.mol)
        report = normalize_imine_bond_orders(build.mol, {
            frozenset((c.atom_idx_1, c.atom_idx_2)): c.distance for c in build.candidates
        })
    assert report["unresolved_components"] > 0
    assert report["changed_bonds"] == []
    assert _signature(build.mol) == before


def _bond_orders(mol):
    return {
        frozenset((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())): int(bond.GetBondTypeAsDouble())
        for bond in mol.GetBonds()
    }


def test_changed_bonds_diagnostics_are_complete(sample):
    record, path = sample
    raw = _raw_build(path)
    repaired = legacy._build_bonded_mol(legacy.read_periodic_cif_atoms(path))
    report = repaired.metadata["imine_bond_order_normalization"]
    before, after = _bond_orders(raw.mol), _bond_orders(repaired.mol)
    actual = {
        (tuple(sorted(pair)), before[pair], after[pair])
        for pair in before
        if before[pair] != after[pair]
    }
    reported = {
        (tuple(change["atoms"]), change["before"], change["after"])
        for change in report["changed_bonds"]
    }
    assert reported == actual
    # The matching re-solves the whole conjugated component, so C=C and C-N
    # double bonds far from the imine links move too.
    assert len(actual) > report["restored_imine_bonds"]


_QUINOID_PROBE_ATOMS = (
    "C1 C 0.00 0.00 0.00 1.0",
    "H2 H 0.10 0.00 0.00 1.0",
    "N3 N 0.25 0.00 0.00 1.0",
    "C4 C 0.35 0.00 0.00 1.0",
    "C5 C 0.45 0.10 0.00 1.0",
    "C6 C 0.55 0.00 0.00 1.0",
    "C7 C 0.65 0.10 0.00 1.0",
    "C8 C 0.75 0.00 0.00 1.0",
    "C9 C 0.85 0.10 0.00 1.0",
)


def _write_quinoid_probe_cif(path, *, parallel_image_distance):
    """A quinoid-misassigned CH=N anchor with a six-carbon alternating tail."""

    lines = [
        "data_quinoid_probe",
        "_cell_length_a 10",
        "_cell_length_b 11",
        "_cell_length_c 12",
        "_cell_angle_alpha 90",
        "_cell_angle_beta 90",
        "_cell_angle_gamma 90",
        "_space_group_name_H-M_alt 'P 1'",
        "loop_",
        "_space_group_symop_operation_xyz",
        "'x,y,z'",
        "loop_",
        "_atom_site_label",
        "_atom_site_type_symbol",
        "_atom_site_fract_x",
        "_atom_site_fract_y",
        "_atom_site_fract_z",
        "_atom_site_occupancy",
        *_QUINOID_PROBE_ATOMS,
        "loop_",
        "_geom_bond_atom_site_label_1",
        "_geom_bond_atom_site_label_2",
        "_geom_bond_site_symmetry_1",
        "_geom_bond_site_symmetry_2",
        "_geom_bond_distance",
        "_ccdc_geom_bond_type",
        # The short CH=N contact is single while the longer N-C anchors carry
        # the double bonds, and the C1-N3 row is repeated for a second image.
        "C1 H2 . . 1.09 S",
        "C1 N3 . . 1.28 S",
        f"C1 N3 . 1_655 {parallel_image_distance:.4f} S",
        "C1 C4 . . 1.35 D",
        "N3 C9 . . 1.41 D",
        "C4 C5 . . 1.46 S",
        "C5 C6 . . 1.35 D",
        "C6 C7 . . 1.46 S",
        "C7 C8 . . 1.35 D",
        "C8 C9 . . 1.46 S",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _raw_build_and_distances(path):
    build = _raw_build(path)
    distances = {
        frozenset((c.atom_idx_1, c.atom_idx_2)): c.distance for c in build.candidates
    }
    return build, distances


@pytest.mark.parametrize(
    ("parallel_image_distance", "repaired"),
    [(1.2805, True), (1.60, False)],
)
def test_parallel_image_length_conflicts_exclude_distance_evidence(
    tmp_path, parallel_image_distance, repaired
):
    path = tmp_path / "quinoid_probe.cif"
    _write_quinoid_probe_cif(path, parallel_image_distance=parallel_image_distance)
    build = legacy._build_bonded_mol(legacy.read_periodic_cif_atoms(path))
    report = build.metadata["imine_bond_order_normalization"]
    # Both image rows stay in the candidate set; only the evidence differs.
    assert sum(1 for c in build.candidates if {c.atom_idx_1, c.atom_idx_2} == {0, 2}) == 2
    assert report["candidate_imine_bonds"] == (1 if repaired else 0)
    assert report["restored_imine_bonds"] == (1 if repaired else 0)
    assert report["unresolved_components"] == 0
    assert build.mol.GetBondBetweenAtoms(0, 2).GetBondTypeAsDouble() == (2.0 if repaired else 1.0)


def test_same_instance_contacts_are_not_treated_as_linkages(tmp_path):
    # cofkit-exported CIFs label atoms `m1_C1`, so atoms of one monomer share a
    # non-empty instance id; such contacts must not become repair targets.
    path = tmp_path / "quinoid_probe.cif"
    _write_quinoid_probe_cif(path, parallel_image_distance=1.2805)
    build, distances = _raw_build_and_distances(path)
    assert build.mol.GetBondBetweenAtoms(0, 2).GetBondTypeAsDouble() == 1.0
    for atom_idx in (0, 2):
        build.mol.GetAtomWithIdx(atom_idx).SetProp("instance_id", "m1")
    report = normalize_imine_bond_orders(build.mol, distances)
    assert report["candidate_imine_bonds"] == 0
    assert report["changed_bonds"] == []
    assert build.mol.GetBondBetweenAtoms(0, 2).GetBondTypeAsDouble() == 1.0
