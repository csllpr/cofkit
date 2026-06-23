import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from cofkit.guest_restart import (
    GuestRestartError,
    build_lammps_guest_restart_state_from_lammps_md_result,
    load_lammps_guest_force_field_assets,
    parse_lammps_guest_restart_snapshot,
    write_graspa_restart_file,
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

    def test_parse_guest_snapshot_keeps_lammps_supercell_box(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            snapshot_path = temp_path / "result_5.data"
            snapshot_path.write_text(
                "\n".join(
                    [
                        "gRASPA movie snapshot",
                        "",
                        "1 atoms",
                        "1 atom types",
                        "0 20 xlo xhi",
                        "0 30 ylo yhi",
                        "0 40 zlo zhi",
                        "5 0.25 -0.5 xy xz yz",
                        "",
                        "Masses",
                        "",
                        "1 131.293 # Xe",
                        "",
                        "Atoms # full",
                        "",
                        "1 1 1 0.0 12.0 15.0 20.0 0 0 0",
                        "",
                    ]
                ),
                encoding="utf-8",
            )

            templates, sites = load_lammps_guest_force_field_assets(("Xe",))
            state = parse_lammps_guest_restart_snapshot(snapshot_path, templates=templates, sites=sites)

        self.assertIsNotNone(state.snapshot_cell)
        assert state.snapshot_cell is not None
        self.assertEqual(state.snapshot_cell.basis[0], (20.0, 0.0, 0.0))
        self.assertEqual(state.snapshot_cell.basis[1], (5.0, 30.0, 0.0))
        self.assertEqual(state.snapshot_cell.basis[2], (0.25, -0.5, 40.0))
        self.assertEqual(state.to_dict()["snapshot_cell"]["unit_cells"], [1, 1, 1])

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

    def test_parse_molecular_atoms_snapshot_with_trailing_label(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            snapshot_path = temp_path / "result_1.data"
            snapshot_path.write_text(
                "\n".join(
                    [
                        "LAMMPS molecular snapshot",
                        "",
                        "1 atoms",
                        "1 atom types",
                        "",
                        "Masses",
                        "",
                        "1 131.293",
                        "",
                        "Atoms # molecular",
                        "",
                        "1 1 1 1.25 2.50 3.75 Xe",
                        "",
                    ]
                ),
                encoding="utf-8",
            )

            templates, sites = load_lammps_guest_force_field_assets(("Xe",))
            state = parse_lammps_guest_restart_snapshot(snapshot_path, templates=templates, sites=sites)

        self.assertEqual(state.n_atoms, 1)
        self.assertEqual((state.atoms[0].x, state.atoms[0].y, state.atoms[0].z), (1.25, 2.50, 3.75))

    def test_parse_bare_atoms_header_with_charge_column_from_graspa_movie(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            snapshot_path = temp_path / "result_0.data"
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
                        "Atoms",
                        "",
                        "1 1 1 0.0 26.94846895 10.74633654 12.64790086 # Xe Xe",
                        "",
                    ]
                ),
                encoding="utf-8",
            )

            templates, sites = load_lammps_guest_force_field_assets(("Xe",))
            state = parse_lammps_guest_restart_snapshot(snapshot_path, templates=templates, sites=sites)

        self.assertEqual(state.n_atoms, 1)
        self.assertEqual(
            (state.atoms[0].x, state.atoms[0].y, state.atoms[0].z),
            (26.94846895, 10.74633654, 12.64790086),
        )

    def test_zero_mass_pseudo_sites_are_rejected_for_lammps_restart(self):
        with self.assertRaises(GuestRestartError) as raised:
            load_lammps_guest_force_field_assets(("TIP4P",))

        self.assertIn("zero or negative mass", str(raised.exception))

    def test_lammps_md_guest_coordinates_write_graspa_restart_file(self):
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
            previous_state = parse_lammps_guest_restart_snapshot(snapshot_path, templates=templates, sites=sites)

            data_path = temp_path / "lammps_md_input.data"
            data_path.write_text(
                "\n".join(
                    [
                        "LAMMPS data",
                        "",
                        "5 atoms",
                        "4 atom types",
                        "",
                        "Masses",
                        "",
                        "1 12.011 # C_3",
                        "2 15.999 # O_3",
                        "3 131.293 # Xe",
                        "4 83.798 # Kr",
                        "",
                        "Atoms # full",
                        "",
                        "1 1 1 0.0 0.0 0.0 0.0 0 0 0 # C_3",
                        "2 1 2 0.0 1.0 0.0 0.0 0 0 0 # O_3",
                        "3 1 2 0.0 0.0 1.0 0.0 0 0 0 # O_3",
                        "4 2 3 0.0 1.25 2.50 3.75 0 0 0 # Xe Xe",
                        "5 3 4 0.0 4.25 5.50 6.75 0 0 0 # Kr Kr",
                        "",
                    ]
                ),
                encoding="utf-8",
            )
            dump_path = temp_path / "lammps_md_trajectory.lammpstrj"
            dump_path.write_text(
                "\n".join(
                    [
                        "ITEM: TIMESTEP",
                        "10",
                        "ITEM: NUMBER OF ATOMS",
                        "5",
                        "ITEM: BOX BOUNDS pp pp pp",
                        "0 20",
                        "0 21",
                        "0 22",
                        "ITEM: ATOMS id type x y z",
                        "1 1 0.0 0.0 0.0",
                        "2 2 1.0 0.0 0.0",
                        "3 2 0.0 1.0 0.0",
                        "4 3 1.25 2.50 3.75",
                        "5 4 4.25 5.50 6.75",
                        "",
                    ]
                ),
                encoding="utf-8",
            )

            md_state, cell = build_lammps_guest_restart_state_from_lammps_md_result(
                SimpleNamespace(lammps_data_path=str(data_path), lammps_dump_path=str(dump_path)),
                previous_guest_restart_state=previous_state,
            )
            restart_result = write_graspa_restart_file(
                md_state,
                temp_path / "restartfile",
                cell=cell,
                component_order=("Xe", "Kr"),
            )
            restart_text = Path(restart_result.restart_file_path).read_text(encoding="utf-8")

        self.assertEqual(md_state.n_atoms, 2)
        self.assertEqual(md_state.components, ("Xe", "Kr"))
        self.assertEqual((md_state.atoms[0].x, md_state.atoms[0].y, md_state.atoms[0].z), (1.25, 2.50, 3.75))
        self.assertEqual(cell.lengths, (20.0, 21.0, 22.0))
        self.assertEqual(restart_result.n_adsorbate_molecules, 2)
        self.assertEqual(restart_result.components, ("Xe", "Kr"))
        self.assertIn("Components: 2 (Adsorbates 2, Cations 0)", restart_text)
        self.assertIn("Components 0 (Xe)", restart_text)
        self.assertIn("Components 1 (Kr)", restart_text)
        self.assertIn("Adsorbate-atom-position: 0 0 1.250000 2.500000 3.750000", restart_text)
        self.assertIn("Adsorbate-atom-position: 0 0 4.250000 5.500000 6.750000", restart_text)
        self.assertIn("Adsorbate-atom-charge: 0 0 0.000000", restart_text)
        self.assertIn("Adsorbate-atom-scaling: 0 0 1", restart_text)
        self.assertIn("Adsorbate-atom-fixed: 0 0 0  0  0", restart_text)


if __name__ == "__main__":
    unittest.main()
