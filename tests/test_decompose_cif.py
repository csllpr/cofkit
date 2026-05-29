import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from cofkit.decompose_cif import (
    atoms_have_explicit_cif_bonds,
    prepare_periodic_cif_atoms,
    read_periodic_cif_atoms,
)


EXAMPLE_CIF = (
    Path(__file__).resolve().parents[1]
    / "reference_repositories"
    / "decofpose-refactored"
    / "examples"
    / "dia_samples"
    / "dia+ALD_4_0.xyz+AMI_2_149.xyz.cif"
)


def _write_explicit_bond_cif(path: Path) -> None:
    path.write_text(
        "\n".join(
            [
                "data_demo",
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
                "C1 C 0.0000 0.0000 0.0000 1.0",
                "N1 N 0.1500 0.0000 0.0000 1.0",
                "loop_",
                "_geom_bond_atom_site_label_1",
                "_geom_bond_atom_site_label_2",
                "_geom_bond_site_symmetry_1",
                "_geom_bond_site_symmetry_2",
                "_geom_bond_distance",
                "_ccdc_geom_bond_type",
                "C1 N1 . . 1.500 S",
                "",
            ]
        ),
        encoding="utf-8",
    )


class DecomposeCifExtractionTests(unittest.TestCase):
    def test_read_periodic_cif_atoms_extracts_sites_without_ase(self):
        atoms = read_periodic_cif_atoms(EXAMPLE_CIF)

        self.assertEqual(len(atoms), 552)
        self.assertEqual(atoms.data_name, "I")
        self.assertEqual(atoms.symbols[:3], ("C", "C", "C"))
        self.assertEqual(atoms.info["_atom_site_label"][:3], ("C1", "C2", "C3"))
        self.assertAlmostEqual(atoms.cell_parameters[0], 38.6606, places=4)
        self.assertAlmostEqual(atoms.cell_parameters[5], 90.3064, places=4)
        self.assertFalse(atoms_have_explicit_cif_bonds(atoms))

    def test_read_periodic_cif_atoms_preserves_explicit_bond_tags(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            cif_path = Path(temp_dir) / "explicit.cif"
            _write_explicit_bond_cif(cif_path)

            atoms = read_periodic_cif_atoms(cif_path)

        self.assertEqual(len(atoms), 2)
        self.assertEqual(atoms.get_chemical_symbols(), ["C", "N"])
        self.assertEqual(atoms.get_atomic_numbers(), [6, 7])
        self.assertEqual(atoms.info["_geom_bond_atom_site_label_1"], ("C1",))
        self.assertEqual(atoms.info["_geom_bond_atom_site_label_2"], ("N1",))
        self.assertEqual(atoms.info["_ccdc_geom_bond_type"], ("S",))
        self.assertTrue(atoms_have_explicit_cif_bonds(atoms))
        self.assertAlmostEqual(atoms.get_positions()[1][0], 1.5, places=6)

    def test_prepare_periodic_cif_atoms_keeps_explicit_bond_cif_in_original_cell(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            cif_path = temp_path / "explicit.cif"
            primitive_dir = temp_path / "primitive"
            _write_explicit_bond_cif(cif_path)

            atoms, used_explicit_bonds = prepare_periodic_cif_atoms(cif_path, primitive_dir)

        self.assertTrue(used_explicit_bonds)
        self.assertEqual(len(atoms), 2)
        self.assertTrue(atoms_have_explicit_cif_bonds(atoms))
        self.assertFalse((primitive_dir / "explicit_primitive.cif").exists())

    def test_prepare_periodic_cif_atoms_writes_primitive_fallback_without_ase(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            primitive_dir = Path(temp_dir) / "primitive"

            atoms, used_explicit_bonds = prepare_periodic_cif_atoms(EXAMPLE_CIF, primitive_dir)

            self.assertFalse(used_explicit_bonds)
            self.assertEqual(len(atoms), 552)
            self.assertFalse(atoms_have_explicit_cif_bonds(atoms))
            self.assertTrue((primitive_dir / f"{EXAMPLE_CIF.stem}_primitive.cif").is_file())

    def test_periodic_cif_atoms_repeat_uses_original_atom_order_per_image(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            cif_path = Path(temp_dir) / "explicit.cif"
            _write_explicit_bond_cif(cif_path)
            atoms = read_periodic_cif_atoms(cif_path)

        repeated = atoms.repeat((3, 1, 1))

        self.assertEqual(len(repeated), 6)
        self.assertEqual(repeated.symbols[:4], ("C", "N", "C", "N"))
        self.assertAlmostEqual(repeated.cell_parameters[0], 30.0, places=6)
        self.assertAlmostEqual(repeated.fractional_positions[0][0], 0.0, places=6)
        self.assertAlmostEqual(repeated.fractional_positions[1][0], 0.05, places=6)
        self.assertAlmostEqual(repeated.fractional_positions[2][0], 1.0 / 3.0, places=6)


if __name__ == "__main__":
    unittest.main()
