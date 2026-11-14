"""Charge presence survives cancellation within or between ledger entries."""

import pytest
from rdkit import Chem

from cofkit.decompose_ledger import (
    LedgerGuest, LedgerUnaccountedFragment, build_atom_ledger, element_counts_for_atoms,
    net_formal_charge_for_atoms,
)


@pytest.mark.parametrize("smiles", ["[Na+].[Cl-]", "[NH3+]CC(=O)[O-]"])
@pytest.mark.parametrize("bucket", ["guests", "residue", "unaccounted"])
def test_neutral_buckets_do_not_hide_charged_atoms(smiles, bucket):
    mol = Chem.MolFromSmiles(smiles)
    fragments = Chem.GetMolFrags(mol)
    common = [dict(
        fragment_id=i, atom_indices=indices, atom_count=len(indices),
        elements=element_counts_for_atoms(mol, indices),
        net_formal_charge=net_formal_charge_for_atoms(mol, indices),
    ) for i, indices in enumerate(fragments)]
    guests = tuple(LedgerGuest(**record, molecular_formula="fixture") for record in common) if bucket == "guests" else ()
    unaccounted = tuple(LedgerUnaccountedFragment(**record, reason="fixture") for record in common) if bucket == "unaccounted" else ()
    residue = tuple(range(mol.GetNumAtoms())) if bucket == "residue" else ()
    ledger = build_atom_ledger(
        family="imine", mol=mol, recovered_fragments=(), guests=guests,
        residue_atom_indices=residue, unaccounted=unaccounted,
        identity_evidence={"formal_charge_column_present": True},
    )
    report = ledger.charge_scope
    assert report.input_net_formal_charge == 0
    assert report.input_charged_atom_count == 2
    assert report.net_charge_delta == 0
    assert report.status == ("guest_localized" if bucket == "guests" else "limited")
    assert report.notes
    assert ledger.balance.balanced
