import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from cofkit.guest_restart import (
    GuestRestartError,
    load_lammps_guest_force_field_assets,
    parse_lammps_guest_restart_snapshot,
)


class GuestRestartTests(unittest.TestCase):
    def test_parse_binary_guest_lammps_snapshot_from_packaged_force_fields(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            snapshot_path = temp_path / "result_5.data"
            snapshot_path.write_text(
                "\n".join(
                    [
                        "gRASPA movie snapshot",
                        "",
                        "2 atoms",
                        "2 atom types",
                        "",
                        "Masses",
                        "",
                        "1 131.293 # Xe",
                        "2 83.798 # Kr",
                        "",
                        "Atoms # full",
                        "",
                        "1 1 1 0.0 1.0 2.0 3.0 0 0 0",
                        "2 2 2 0.0 4.0 5.0 6.0 0 0 0",
                        "",
                    ]
                ),
                encoding="utf-8",
            )

            templates, sites = load_lammps_guest_force_field_assets(("Xe", "Kr"))
            state = parse_lammps_guest_restart_snapshot(snapshot_path, templates=templates, sites=sites)

        self.assertEqual(state.n_atoms, 2)
        self.assertEqual(state.components, ("Xe", "Kr"))
        self.assertEqual([atom.site_label for atom in state.atoms], ["Xe", "Kr"])
        self.assertEqual(state.atoms[0].x, 1.0)
        self.assertEqual(state.atoms[1].z, 6.0)
        self.assertGreater(state.site_by_label()["Xe"].epsilon_kcal_per_mol, 0.0)

    def test_parse_molecular_atoms_snapshot_coordinates_without_charge_column(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            snapshot_path = temp_path / "result_1.data"
            snapshot_path.write_text(
                "\n".join(
                    [
                        "gRASPA movie snapshot",
                        "",
                        "1 atoms",
                        "1 atom types",
                        "",
                        "Masses",
                        "",
                        "1 131.293 # Xe",
                        "",
                        "Atoms # molecular",
                        "",
                        "1 1 1 1.25 2.50 3.75 0 0 0",
                        "",
                    ]
                ),
                encoding="utf-8",
            )

            templates, sites = load_lammps_guest_force_field_assets(("Xe",))
            state = parse_lammps_guest_restart_snapshot(snapshot_path, templates=templates, sites=sites)

        self.assertEqual(state.n_atoms, 1)
        self.assertEqual((state.atoms[0].x, state.atoms[0].y, state.atoms[0].z), (1.25, 2.50, 3.75))

    def test_zero_mass_pseudo_sites_are_rejected_for_lammps_restart(self):
        with self.assertRaises(GuestRestartError) as raised:
            load_lammps_guest_force_field_assets(("TIP4P",))

        self.assertIn("zero or negative mass", str(raised.exception))


if __name__ == "__main__":
    unittest.main()
