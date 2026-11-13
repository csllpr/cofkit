"""Decomposition must confine bond-order changes to the inverse reaction."""

from dataclasses import replace

import pytest
from rdkit import Chem

from cofkit import decompose as legacy
from cofkit import decompose_events as events
from test_decompose import _write_molecular_probe_cif


@pytest.mark.parametrize(
    "smiles",
    [
        "CN=CC1=CC=CC(=O)C1=O",  # exocyclic imine next to a ring carbonyl
        "CN=CC1C(=O)C=CC=C1",  # explicit-H completion must also be preserved
        "O=C1C=CN=CC1",  # endocyclic C=N is part of the precursor scaffold
        "CNC=C1C(=O)C=CC=C1",  # the keto-enamine representation
        "CC(=O)Nc1ccc(C=NC)cc1",  # unrelated amide control
    ],
)
def test_graph_loading_preserves_supplied_carbonyl_and_imine_orders(tmp_path, smiles):
    source = Chem.MolFromSmiles(smiles)
    path = _write_molecular_probe_cif(smiles, tmp_path / "probe.cif")
    built = legacy._build_bonded_mol(legacy.read_periodic_cif_atoms(path))
    expected = {
        frozenset((b.GetBeginAtomIdx(), b.GetEndAtomIdx())): b.GetBondTypeAsDouble()
        for b in source.GetBonds()
    }
    actual = {
        frozenset((b.GetBeginAtomIdx(), b.GetEndAtomIdx())): b.GetBondTypeAsDouble()
        for b in built.mol.GetBonds()
    }
    assert actual == expected
    # Sanitization must not encounter valence errors introduced by loading.
    Chem.SanitizeMol(built.mol)
    assert Chem.MolToSmiles(built.mol) == Chem.MolToSmiles(source)
    assert built.metadata["n_beta_ketoenamine_bond_pairs_normalized"] == 0
    assert built.metadata["beta_ketoenamine_bond_order_policy"] == "detect_without_pre_cut_rewrite"
    assert [c.explicit_order for c in built.candidates] == [
        expected[frozenset((c.atom_idx_1, c.atom_idx_2))] for c in built.candidates
    ]


@pytest.mark.parametrize("smiles", ["CN=CC1C(=O)C=CC=C1", "CNC=C1C(=O)C=CC=C1"])
def test_both_beta_ketoenamine_representations_are_detected_without_rewriting(tmp_path, smiles):
    path = _write_molecular_probe_cif(smiles, tmp_path / "probe.cif")
    built = legacy._build_bonded_mol(legacy.read_periodic_cif_atoms(path))
    before = Chem.MolToSmiles(built.mol)
    event_detection = events.detect_linkage_events(built)
    bken = [e for e in event_detection.events if e.family == "bken"]
    assert len(bken) == 1
    legacy_detection = legacy._classify_nitrogen_linkage_bonds(built.mol)
    assert len(legacy_detection.bken_bond_indices) == 1
    assert Chem.MolToSmiles(built.mol) == before


@pytest.mark.parametrize(
    "scaffold",
    [
        "CC1CCC(=O)CC1",  # remote cyclic ketone
        "CN1CCCC1=O",  # remote lactam
        "CC1=CC(=O)C=CC1=O",  # remote quinone
        "c1ccc(C=O)cc1",  # pendant aldehyde
        "CC(=O)C",  # pendant acyclic ketone
    ],
)
def test_keto_aldehyde_repair_preserves_remote_scaffold_carbonyls(scaffold):
    product_fragment = Chem.MolFromSmiles("C=C1C(=O)C=CC=C1" + scaffold)
    assert product_fragment is not None
    # The cut linkage carbon is the initial exocyclic carbon, not any of the
    # pendant scaffold's carbonyls. Only its adjacent carbocycle is restored.
    product_fragment.GetAtomWithIdx(0).SetProp("cofkit_decompose_role", "keto_aldehyde")
    repaired = legacy._restore_keto_aldehyde_fragment(product_fragment)
    monomer = legacy._finalize_repaired_fragment(repaired, "keto_aldehyde")
    expected = Chem.MolFromSmiles("O=Cc1c(O)cccc1" + scaffold)
    assert monomer.canonical_smiles == Chem.MolToSmiles(expected)
    assert monomer.connectivity == 1


def test_invalid_explicit_bond_orders_are_not_silently_guessed(tmp_path):
    # An erroneous C=O in an ether gives O valence three. Removing the broad
    # rewrite is not permission to invent a different input structure.
    path = _write_molecular_probe_cif("COC", tmp_path / "invalid.cif")
    text = path.read_text()
    rows = text.splitlines()
    first_bond = next(i for i, row in enumerate(rows) if row.startswith("C1 O2 "))
    rows[first_bond] = rows[first_bond].rsplit(" ", 1)[0] + " D"
    path.write_text("\n".join(rows) + "\n")
    built = legacy._build_bonded_mol(legacy.read_periodic_cif_atoms(path))
    assert built.mol.GetBondBetweenAtoms(0, 1).GetBondTypeAsDouble() == 2.0
    assert any(p.GetType() == "AtomValenceException" for p in Chem.DetectChemistryProblems(built.mol))


def _competing_bond_order_graphs(tmp_path):
    # Two individually valid representations share the same periodic heavy
    # edges, but their explicit reaction bonds encode different families.
    paths = [
        _write_molecular_probe_cif(smiles, tmp_path / f"probe_{i}.cif")
        for i, smiles in enumerate(("CNC=C1C(=O)C=CC=C1", "CN=Cc1c(O)cccc1"))
    ]
    primary, fallback = (
        legacy._build_bonded_mol(legacy.read_periodic_cif_atoms(path)) for path in paths
    )
    primary_events = events.detect_linkage_events(primary).events
    assert [e.family for e in primary_events] == ["bken"]
    assert [e.family for e in events.detect_linkage_events(fallback).events] == ["imine"]
    candidate = legacy.CifDecompositionResult(
        status="success", input_cif=str(paths[0]), topology="hcb", linkage="imine",
        cofid="candidate retained for review",
    )
    return paths[0], primary, fallback, primary_events, candidate


def test_distance_fallback_retains_conflicting_explicit_linkage_evidence(tmp_path):
    _, primary, fallback, primary_events, candidate = _competing_bond_order_graphs(tmp_path)
    conflicts = events._distance_fallback_bond_order_conflicts(primary, fallback, primary_events, candidate)
    assert len(conflicts) == 1
    assert conflicts[0]["explicit_family"] == "bken"
    assert conflicts[0]["fallback_family"] == "imine"
    change = conflicts[0]["changed_cut_bonds"]
    assert len(change) == 1
    assert change[0]["explicit_order"] == 1.0
    assert change[0]["inferred_order"] == 2.0


@pytest.mark.parametrize("condition", ["missing_orders", "distance_primary", "low_confidence", "same_family", "failed_candidate", "different_periodic_graph"])
def test_distance_fallback_conflict_guard_requires_all_its_evidence(tmp_path, condition):
    _, primary, fallback, primary_events, candidate = _competing_bond_order_graphs(tmp_path)
    if condition == "missing_orders":
        primary = replace(primary, metadata={**primary.metadata, "n_missing_cif_bond_orders": 1})
    elif condition == "distance_primary":
        primary = replace(primary, metadata={**primary.metadata, "bond_source": "distance_inferred"})
    elif condition == "low_confidence":
        primary_events = tuple(replace(e, confidence="low") for e in primary_events)
    elif condition == "same_family":
        candidate = replace(candidate, linkage="bken")
    elif condition == "failed_candidate":
        candidate = replace(candidate, status="skipped", cofid=None)
    else:
        edge = fallback.candidates[0]
        fallback = replace(fallback, candidates=(*fallback.candidates, replace(edge, periodic_image=(1, 0, 0))))
    assert events._distance_fallback_bond_order_conflicts(primary, fallback, primary_events, candidate) == []


def test_public_fallback_reports_ambiguity_and_retains_candidate(tmp_path, monkeypatch):
    path, primary, fallback, _, candidate = _competing_bond_order_graphs(tmp_path)
    failed = legacy.CifDecompositionResult(
        status="skipped", input_cif=str(path), topology="hcb", linkage="bken",
        reason="incomplete endpoint accounting",
        metadata={"event_status": events.EVENT_STATUS_ENDPOINT, "successful_hypothesis_count": 0},
    )
    builds = iter((primary, fallback))
    results = iter((failed, candidate))
    monkeypatch.setattr(legacy, "_build_bonded_mol", lambda *a, **k: next(builds))
    # Isolate the selection boundary: both graph interpretations have already
    # been evaluated; a sole buildable fallback is insufficient to disambiguate.
    monkeypatch.setattr(events, "_generate_and_validate_hypotheses", lambda *a, **k: ((), {}))
    monkeypatch.setattr(events, "_select_event_result", lambda *a, **k: next(results))
    result = events.decompose_cif_to_cofid_event(path, topology="hcb")
    assert result.status == "ambiguous"
    assert result.cofid is None
    assert not result.ok
    assert result.metadata["event_status"] == events.EVENT_STATUS_BOND_ORDER_AMBIGUOUS
    assert result.metadata["bond_order_evidence_conflicts"]
    assert result.metadata["bond_graph_fallback"]["applied"] is False
    assert result.metadata["bond_graph_fallback"]["fallback_cofid"] == candidate.cofid
