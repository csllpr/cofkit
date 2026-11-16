"""Equivalent CIF aromatic notation must not change UFF/DREIDING types."""

from dataclasses import replace

import pytest
from rdkit import Chem
from rdkit.Chem import rdDepictor

from cofkit import lammps as lmp
from cofkit.bond_types import bond_order_to_cif_type


def _cif_from_molecule(tmp_path, molecule, name):
    path = tmp_path / f"{name}.cif"
    conformer = molecule.GetConformer()
    lines = [
        "data_typing", "_cell_length_a 40", "_cell_length_b 40", "_cell_length_c 40",
        "_cell_angle_alpha 90", "_cell_angle_beta 90", "_cell_angle_gamma 90",
        "_space_group_name_H-M_alt 'P 1'", "_space_group_IT_number 1",
        "loop_", "_atom_site_label", "_atom_site_type_symbol",
        "_atom_site_fract_x", "_atom_site_fract_y", "_atom_site_fract_z",
    ]
    for atom in molecule.GetAtoms():
        i = atom.GetIdx()
        point = conformer.GetAtomPosition(i)
        lines.append(f"a{i} {atom.GetSymbol()} {(point.x + 20)/40} {(point.y + 20)/40} 0.5")
    lines.extend([
        "loop_", "_geom_bond_atom_site_label_1", "_geom_bond_atom_site_label_2",
        "_geom_bond_site_symmetry_1", "_geom_bond_site_symmetry_2",
        "_geom_bond_distance", "_ccdc_geom_bond_type",
    ])
    for bond in molecule.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        distance = (conformer.GetAtomPosition(i) - conformer.GetAtomPosition(j)).Length()
        lines.append(f"a{i} a{j} . . {distance} {bond_order_to_cif_type(bond.GetBondTypeAsDouble())}")
    path.write_text("\n".join(lines) + "\n")
    return path, lmp._parse_explicit_bond_cif(path)


def _types(parsed):
    images, _ = lmp._compute_unwrapped_atom_images(parsed)
    molecule = lmp._build_openbabel_molecule(parsed, lmp._unwrapped_cartesian_positions(parsed, images))
    return lmp._assign_openbabel_uff_atom_types(molecule), lmp._assign_dreiding_atom_types(molecule)


@pytest.mark.parametrize("smiles", [
    "c1ccc(-c2ncncn2)cc1",  # Mixed benzene/triazine notation hid triazine aromaticity.
    "B1(OC2=CC=CC=C2O1)c1ccccc1",  # Boronate O_R versus O_2.
    "c1ccc(B2OB(c3ccccc3)OB(c3ccccc3)O2)cc1",  # Boroxine O_R versus O_3.
    "c1ccc(-c2cc[nH]c2)cc1",  # Keep the explicit pyrrolic N-H.
    "c1ccccc1C#CC=C",  # Nonaromatic triple/double bonds must survive.
])
def test_typing_is_invariant_to_aromatic_kekule_and_mixed_notation(tmp_path, smiles):
    aromatic = Chem.AddHs(Chem.MolFromSmiles(smiles))
    rdDepictor.Compute2DCoords(aromatic)
    kekule = Chem.Mol(aromatic)
    Chem.Kekulize(kekule, clearAromaticFlags=True)
    mixed = Chem.Mol(kekule)
    benzene = next(ring for ring in aromatic.GetRingInfo().BondRings()
                   if len(ring) == 6 and all(
                       aromatic.GetBondWithIdx(i).GetBeginAtom().GetSymbol() == "C"
                       and aromatic.GetBondWithIdx(i).GetEndAtom().GetSymbol() == "C"
                       for i in ring))
    for i in benzene:
        bond = mixed.GetBondWithIdx(i)
        bond.SetBondType(Chem.BondType.AROMATIC)
        bond.SetIsAromatic(True)
        bond.GetBeginAtom().SetIsAromatic(True)
        bond.GetEndAtom().SetIsAromatic(True)
    types = []
    normalized = []
    for name, mol in [("aromatic", aromatic), ("kekule", kekule), ("mixed", mixed)]:
        path, parsed = _cif_from_molecule(tmp_path, mol, name)
        before = path.read_bytes()
        original_bonds = parsed.bonds
        types.append(_types(parsed))
        normalized.append(lmp._forcefield_typing_bond_orders(parsed))
        assert parsed.bonds == original_bonds
        assert path.read_bytes() == before
        for bond, order in zip(parsed.bonds, normalized[-1]):
            if order != 1.5:
                assert order == bond.bond_order
    assert normalized[0] == normalized[1] == normalized[2]
    assert types[0] == types[1] == types[2]
    if "ncncn" in smiles:
        assert {types[0][1][a.GetIdx() + 1] for a in aromatic.GetAtoms() if a.GetSymbol() == "N"} == {"N_R"}
    order = list(reversed(range(mixed.GetNumAtoms())))
    _, reordered = _cif_from_molecule(tmp_path, Chem.RenumberAtoms(mixed, order), "reordered")
    for expected, actual in zip(types[0], _types(reordered)):
        assert {order[i - 1] + 1: value for i, value in actual.items()} == expected


def test_nonaromatic_ring_typing_is_independent_of_aromatic_spectator(tmp_path):
    core = Chem.AddHs(Chem.MolFromSmiles("OB1OB(O)OB(O)O1"))
    rdDepictor.Compute2DCoords(core)
    _, parsed = _cif_from_molecule(tmp_path, core, "core")
    spectator = Chem.AddHs(Chem.MolFromSmiles("c1ccccc1"))
    rdDepictor.Compute2DCoords(spectator)
    _, combined = _cif_from_molecule(tmp_path, Chem.CombineMols(core, spectator), "spectator")
    for expected, actual in zip(_types(parsed), _types(combined)):
        assert {i: value for i, value in actual.items() if i <= core.GetNumAtoms()} == expected


@pytest.mark.parametrize("order, message", [
    (None, "requires explicit bond orders"),
    (2.5, "Unsupported explicit bond order"),
])
def test_typing_requires_supported_explicit_bond_orders(tmp_path, order, message):
    molecule = Chem.AddHs(Chem.MolFromSmiles("CC"))
    rdDepictor.Compute2DCoords(molecule)
    _, parsed = _cif_from_molecule(tmp_path, molecule, "orders")
    invalid = replace(parsed, bonds=(replace(parsed.bonds[0], bond_order=order), *parsed.bonds[1:]))
    with pytest.raises(lmp.LammpsInputError, match=message):
        _types(invalid)


def test_typing_normalization_reports_unresolvable_aromatic_notation(tmp_path):
    molecule = Chem.AddHs(Chem.MolFromSmiles("C1CC1"))
    rdDepictor.Compute2DCoords(molecule)
    _, parsed = _cif_from_molecule(tmp_path, molecule, "invalid")
    # An aromatic bond outside any ring has no equivalent Kekule notation.
    invalid = replace(parsed, bonds=tuple(
        replace(bond, bond_order=1.5) if bond.atom_id_2 > 3 else bond for bond in parsed.bonds
    ))
    with pytest.raises(lmp.LammpsInputError, match="Cannot normalize aromatic bond notation"):
        _types(invalid)
