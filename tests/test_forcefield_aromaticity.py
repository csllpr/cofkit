"""Equivalent CIF aromatic notation must not change UFF/DREIDING types."""

from dataclasses import replace

import pytest
from rdkit import Chem
from rdkit.Chem import rdDepictor

from cofkit import lammps as lmp
from cofkit.bond_types import bond_order_to_cif_type
from cofkit.periodic_aromaticity import normalize_periodic_bond_orders


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


def test_winding_polyene_is_not_a_molecular_benzene_ring(tmp_path):
    molecule = Chem.AddHs(Chem.MolFromSmiles("c1ccccc1"))
    rdDepictor.Compute2DCoords(molecule)
    Chem.Kekulize(molecule, clearAromaticFlags=True)
    _, parsed = _cif_from_molecule(tmp_path, molecule, "winding")
    # The last C-C edge connects to the NEXT cell: the infinite graph is a
    # polyene chain, despite the six-cycle in its finite quotient.
    bonds = tuple(replace(b, shift_2=(1, 0, 0)) if {b.atom_id_1, b.atom_id_2} == {1, 6} else b
                  for b in parsed.bonds)
    chain = replace(parsed, bonds=bonds)
    assert lmp._forcefield_typing_bond_orders(chain) == tuple(b.bond_order for b in bonds)
    assert {t for t in _types(chain)[1].values() if t.startswith("C")} == {"C_2"}


def test_real_aromatic_ring_survives_boundary_crossing_and_basis_change(tmp_path):
    molecule = Chem.AddHs(Chem.MolFromSmiles("c1ccc2[nH]ccc2c1"))
    rdDepictor.Compute2DCoords(molecule)
    Chem.Kekulize(molecule, clearAromaticFlags=True)
    _, parsed = _cif_from_molecule(tmp_path, molecule, "boundary")
    # A change of atom representatives can move edges across arbitrarily many
    # cells, without changing the infinite graph or any molecular ring.
    gauges = {a.atom_id: ((a.atom_id % 3) * 7, -(a.atom_id % 2) * 5, 0) for a in parsed.atoms}
    moved = replace(
        parsed,
        atoms=tuple(replace(a, fractional_position=lmp._add_shift(a.fractional_position, gauges[a.atom_id]))
                    for a in parsed.atoms),
        bonds=tuple(replace(b, shift_1=lmp._sub_shift(b.shift_1, gauges[b.atom_id_1]),
                            shift_2=lmp._sub_shift(b.shift_2, gauges[b.atom_id_2])) for b in parsed.bonds),
    )
    assert _types(moved) == _types(parsed)
    edges = [(b.atom_id_1 - 1, b.atom_id_2 - 1, lmp._sub_shift(b.shift_2, b.shift_1), b.bond_order)
             for b in moved.bonds]
    rebased = [(u, v, (g[0] + g[1], g[1], g[2]), order) for u, v, g, order in edges]
    symbols = [a.symbol for a in moved.atoms]
    assert normalize_periodic_bond_orders(symbols, edges) == normalize_periodic_bond_orders(symbols, rebased)


def _explicit_supercell(parsed, factors):
    import itertools
    cells = list(itertools.product(*(range(x) for x in factors)))
    lookup = {cell: i for i, cell in enumerate(cells)}
    n = len(parsed.atoms)
    atoms, bonds = [], []
    for cell in cells:
        atoms.extend(replace(a, atom_id=lookup[cell] * n + a.atom_id,
                             fractional_position=tuple((a.fractional_position[i] + cell[i]) / factors[i]
                                                       for i in range(3))) for a in parsed.atoms)
        for b in parsed.bonds:
            first = tuple(cell[i] + b.shift_1[i] for i in range(3))
            second = tuple(cell[i] + b.shift_2[i] for i in range(3))
            bonds.append(replace(
                b, bond_id=len(bonds) + 1,
                atom_id_1=lookup[tuple(first[i] % factors[i] for i in range(3))] * n + b.atom_id_1,
                atom_id_2=lookup[tuple(second[i] % factors[i] for i in range(3))] * n + b.atom_id_2,
                shift_1=tuple(first[i] // factors[i] for i in range(3)),
                shift_2=tuple(second[i] // factors[i] for i in range(3)),
            ))
    return replace(parsed, atoms=tuple(atoms), bonds=tuple(bonds),
                   lammps_basis=tuple(tuple(x * factors[i] for x in row)
                                      for i, row in enumerate(parsed.lammps_basis)))


@pytest.mark.parametrize("factors", [(2, 1, 1), (2, 3, 1)])
def test_typing_and_bond_comparison_are_invariant_to_polyene_supercells(tmp_path, factors):
    molecule = Chem.AddHs(Chem.MolFromSmiles("c1ccccc1"))
    rdDepictor.Compute2DCoords(molecule)
    Chem.Kekulize(molecule, clearAromaticFlags=True)
    _, parsed = _cif_from_molecule(tmp_path, molecule, "supercell")
    parsed = replace(parsed, bonds=tuple(
        replace(b, shift_2=(1, 0, 0)) if {b.atom_id_1, b.atom_id_2} == {1, 6} else b for b in parsed.bonds
    ))
    big = _explicit_supercell(parsed, factors)
    expected = _types(parsed)
    actual = _types(big)
    for small_types, big_types in zip(expected, actual):
        assert all(value == small_types[(i - 1) % len(parsed.atoms) + 1] for i, value in big_types.items())
    assert lmp._forcefield_typing_bond_orders(big) == lmp._forcefield_typing_bond_orders(parsed) * (len(big.bonds) // len(parsed.bonds))
