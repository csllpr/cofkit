import json
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace


from cofkit.guest_restart import (
    GuestRestartError,
    build_lammps_guest_restart_state_from_lammps_md_result,
    guest_site_model_diagnostics,
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

            templates, sites = load_lammps_guest_force_field_assets(("Xe_GENERICMOFS", "Kr_GENERICMOFS"))
            state = parse_lammps_guest_restart_snapshot(snapshot_path, templates=templates, sites=sites)

        self.assertEqual(state.n_atoms, 2)
        self.assertEqual(state.components, ("Xe_GENERICMOFS", "Kr_GENERICMOFS"))
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

            templates, sites = load_lammps_guest_force_field_assets(("Xe_GENERICMOFS",))
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

            templates, sites = load_lammps_guest_force_field_assets(("Xe_GENERICMOFS",))
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

            templates, sites = load_lammps_guest_force_field_assets(("Xe_GENERICMOFS",))
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

            templates, sites = load_lammps_guest_force_field_assets(("Xe_GENERICMOFS",))
            state = parse_lammps_guest_restart_snapshot(snapshot_path, templates=templates, sites=sites)

        self.assertEqual(state.n_atoms, 1)
        self.assertEqual(
            (state.atoms[0].x, state.atoms[0].y, state.atoms[0].z),
            (26.94846895, 10.74633654, 12.64790086),
        )

    def test_parse_graspa_movie_component_then_site_comment(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            snapshot_path = temp_path / "result_0.data"
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
                        "Atoms",
                        "",
                        "1 1 1 0.0 1.0 2.0 3.0 # Xe_GENERICMOFS Xe",
                        "2 2 2 0.0 4.0 5.0 6.0 # Kr_GENERICMOFS Kr",
                        "",
                    ]
                ),
                encoding="utf-8",
            )

            templates, sites = load_lammps_guest_force_field_assets(("Xe_GENERICMOFS", "Kr_GENERICMOFS"))
            state = parse_lammps_guest_restart_snapshot(snapshot_path, templates=templates, sites=sites)

        self.assertEqual(state.n_atoms, 2)
        self.assertEqual(state.components, ("Xe_GENERICMOFS", "Kr_GENERICMOFS"))
        self.assertEqual([atom.site_label for atom in state.atoms], ["Xe", "Kr"])

    def test_parse_empty_guest_population_as_valid_restart_state(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            snapshot_path = temp_path / "result_0.data"
            snapshot_path.write_text(
                "\n".join(
                    [
                        "gRASPA movie snapshot",
                        "",
                        "2 atoms",
                        "4 atom types",
                        "",
                        "Masses",
                        "",
                        "1 12.011 # C",
                        "2 1.008 # H",
                        "3 131.293 # Xe",
                        "4 83.798 # Kr",
                        "",
                        "Atoms",
                        "",
                        "1 1 1 0.0 1.0 2.0 3.0 # C",
                        "2 1 2 0.0 4.0 5.0 6.0 # H",
                        "",
                    ]
                ),
                encoding="utf-8",
            )

            templates, sites = load_lammps_guest_force_field_assets(("Xe_GENERICMOFS", "Kr_GENERICMOFS"))
            state = parse_lammps_guest_restart_snapshot(snapshot_path, templates=templates, sites=sites)

        self.assertEqual(state.n_atoms, 0)
        self.assertEqual(state.components, ())
        self.assertTrue(
            any("adsorbate population is being treated as empty" in warning for warning in state.warnings)
        )
        self.assertEqual(state.skipped_unknown_site_counts, (("C", 1), ("H", 1)))
        self.assertTrue(
            any("counts by label: C: 1, H: 1" in warning for warning in state.warnings)
        )

    def test_parse_empty_guest_population_rejects_incomplete_atom_table(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            snapshot_path = temp_path / "result_0.data"
            snapshot_path.write_text(
                "\n".join(
                    [
                        "gRASPA movie snapshot",
                        "",
                        "2 atoms",
                        "4 atom types",
                        "",
                        "Masses",
                        "",
                        "1 12.011 # C",
                        "2 1.008 # H",
                        "3 131.293 # Xe",
                        "4 83.798 # Kr",
                        "",
                        "Atoms",
                        "",
                        "1 1 1 0.0 1.0 2.0 3.0 # C",
                        "2 1 2 0.0 invalid 5.0 6.0 # H",
                        "",
                    ]
                ),
                encoding="utf-8",
            )

            templates, sites = load_lammps_guest_force_field_assets(("Xe_GENERICMOFS", "Kr_GENERICMOFS"))
            with self.assertRaises(GuestRestartError) as raised:
                parse_lammps_guest_restart_snapshot(snapshot_path, templates=templates, sites=sites)

        self.assertIn("declares 2 atoms, but 1 atom rows were parsed", str(raised.exception))

    def test_parse_empty_guest_population_requires_requested_guest_types(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            snapshot_path = temp_path / "result_0.data"
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
                        "1 12.011 # C",
                        "2 1.008 # H",
                        "",
                        "Atoms",
                        "",
                        "1 1 1 0.0 1.0 2.0 3.0 # C",
                        "2 1 2 0.0 4.0 5.0 6.0 # H",
                        "",
                    ]
                ),
                encoding="utf-8",
            )

            templates, sites = load_lammps_guest_force_field_assets(("Xe_GENERICMOFS", "Kr_GENERICMOFS"))
            with self.assertRaises(GuestRestartError) as raised:
                parse_lammps_guest_restart_snapshot(snapshot_path, templates=templates, sites=sites)

        self.assertIn("does not declare the requested guest site type(s) 'Kr', 'Xe'", str(raised.exception))

    def test_zero_mass_pseudo_sites_are_rejected_for_lammps_restart(self):
        with self.assertRaises(GuestRestartError) as raised:
            load_lammps_guest_force_field_assets(("TIP4P_DREIDING",))

        self.assertIn("zero or negative mass", str(raised.exception))

    def test_parse_complete_co2_multisite_snapshot_from_packaged_force_fields(self):
        # Regression for repeated site labels: CO2 templates enumerate
        # O_co2 twice, which must not make the component ambiguous.
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            snapshot_path = temp_path / "result_0.data"
            snapshot_path.write_text(
                "\n".join(
                    [
                        "gRASPA movie snapshot",
                        "",
                        "6 atoms",
                        "2 atom types",
                        "",
                        "Masses",
                        "",
                        "1 15.9994 # O_co2",
                        "2 12.0107 # C_co2",
                        "",
                        "Atoms # full",
                        "",
                        "1 1 1 -0.35 1.0 2.0 3.16 0 0 0 # CO2_DREIDING O_co2",
                        "2 1 2 0.70 1.0 2.0 2.0 0 0 0 # CO2_DREIDING C_co2",
                        "3 1 1 -0.35 1.0 2.0 0.84 0 0 0 # CO2_DREIDING O_co2",
                        "4 2 1 -0.35 5.0 6.0 7.16 0 0 0 # CO2_DREIDING O_co2",
                        "5 2 2 0.70 5.0 6.0 6.0 0 0 0 # CO2_DREIDING C_co2",
                        "6 2 1 -0.35 5.0 6.0 4.84 0 0 0 # CO2_DREIDING O_co2",
                        "",
                    ]
                ),
                encoding="utf-8",
            )

            templates, sites = load_lammps_guest_force_field_assets(("CO2_DREIDING",))
            state = parse_lammps_guest_restart_snapshot(snapshot_path, templates=templates, sites=sites)

        self.assertEqual(state.n_atoms, 6)
        self.assertEqual(state.components, ("CO2_DREIDING",))
        self.assertEqual([atom.site_label for atom in state.atoms[:3]], ["O_co2", "C_co2", "O_co2"])
        self.assertEqual([atom.molecule_key for atom in state.atoms], ["CO2_DREIDING:1"] * 3 + ["CO2_DREIDING:2"] * 3)
        self.assertEqual(state.warnings, ())
        self.assertEqual(state.skipped_molecules, ())
        self.assertEqual(state.skipped_unknown_site_counts, ())
        self.assertEqual(state.n_skipped_ambiguous_atoms, 0)
        self.assertGreater(state.site_by_label()["O_co2"].mass, 0.0)

    def test_parse_complete_so2_massive_multisite_snapshot(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            snapshot_path = temp_path / "result_0.data"
            snapshot_path.write_text(
                "\n".join(
                    [
                        "gRASPA movie snapshot",
                        "",
                        "3 atoms",
                        "2 atom types",
                        "",
                        "Masses",
                        "",
                        "1 15.9994 # O_so2",
                        "2 32.065 # S_so2",
                        "",
                        "Atoms # full",
                        "",
                        "1 1 1 -0.295 0.0 0.507 0.0 0 0 0 # SO2_DREIDING O_so2",
                        "2 1 2 0.59 0.8619 0.0 0.0 0 0 0 # SO2_DREIDING S_so2",
                        "3 1 1 -0.295 1.7239 0.507 0.0 0 0 0 # SO2_DREIDING O_so2",
                        "",
                    ]
                ),
                encoding="utf-8",
            )

            templates, sites = load_lammps_guest_force_field_assets(("SO2_DREIDING",))
            state = parse_lammps_guest_restart_snapshot(snapshot_path, templates=templates, sites=sites)

        self.assertEqual(state.n_atoms, 3)
        self.assertEqual(state.components, ("SO2_DREIDING",))
        self.assertEqual([atom.site_label for atom in state.atoms], ["O_so2", "S_so2", "O_so2"])
        self.assertEqual(state.warnings, ())

    def test_mislabeled_guest_rows_are_explicitly_accounted(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            snapshot_path = temp_path / "result_0.data"
            snapshot_path.write_text(
                "\n".join(
                    [
                        "gRASPA movie snapshot",
                        "",
                        "5 atoms",
                        "4 atom types",
                        "",
                        "Masses",
                        "",
                        "1 15.9994 # O_co2",
                        "2 12.0107 # C_co2",
                        "3 12.011 # C_3",
                        "4 15.9994 # O_co2_typo",
                        "",
                        "Atoms # full",
                        "",
                        "1 1 1 -0.35 1.0 2.0 3.16 0 0 0 # CO2_DREIDING O_co2",
                        "2 1 2 0.70 1.0 2.0 2.0 0 0 0 # CO2_DREIDING C_co2",
                        "3 1 1 -0.35 1.0 2.0 0.84 0 0 0 # CO2_DREIDING O_co2",
                        "4 2 3 0.0 0.0 0.0 0.0 0 0 0 # C_3",
                        "5 2 4 -0.35 4.0 4.0 4.0 0 0 0 # CO2_DREIDING O_co2_typo",
                        "",
                    ]
                ),
                encoding="utf-8",
            )

            templates, sites = load_lammps_guest_force_field_assets(("CO2_DREIDING",))
            state = parse_lammps_guest_restart_snapshot(snapshot_path, templates=templates, sites=sites)

        # The complete CO2 molecule survives; the framework row and the
        # mislabeled O_co2_typo row are excluded but explicitly accounted.
        self.assertEqual(state.n_atoms, 3)
        self.assertEqual(state.skipped_unknown_site_counts, (("C_3", 1), ("O_co2_typo", 1)))
        self.assertTrue(
            any("counts by label: C_3: 1, O_co2_typo: 1" in warning for warning in state.warnings)
        )

    def test_incomplete_guest_molecule_is_explicitly_accounted(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            snapshot_path = temp_path / "result_0.data"
            snapshot_path.write_text(
                "\n".join(
                    [
                        "gRASPA movie snapshot",
                        "",
                        "5 atoms",
                        "2 atom types",
                        "",
                        "Masses",
                        "",
                        "1 15.9994 # O_co2",
                        "2 12.0107 # C_co2",
                        "",
                        "Atoms # full",
                        "",
                        "1 1 1 -0.35 1.0 2.0 3.16 0 0 0 # CO2_DREIDING O_co2",
                        "2 1 2 0.70 1.0 2.0 2.0 0 0 0 # CO2_DREIDING C_co2",
                        "3 1 1 -0.35 1.0 2.0 0.84 0 0 0 # CO2_DREIDING O_co2",
                        "4 2 1 -0.35 5.0 6.0 7.16 0 0 0 # CO2_DREIDING O_co2",
                        "5 2 2 0.70 5.0 6.0 6.0 0 0 0 # CO2_DREIDING C_co2",
                        "",
                    ]
                ),
                encoding="utf-8",
            )

            templates, sites = load_lammps_guest_force_field_assets(("CO2_DREIDING",))
            state = parse_lammps_guest_restart_snapshot(snapshot_path, templates=templates, sites=sites)

        self.assertEqual(state.n_atoms, 3)
        self.assertEqual(len(state.skipped_molecules), 1)
        skipped = state.skipped_molecules[0]
        self.assertEqual(skipped.component, "CO2_DREIDING")
        self.assertEqual(skipped.molecule_key, "CO2_DREIDING:2")
        self.assertEqual(skipped.n_sites_found, 2)
        self.assertEqual(skipped.n_sites_expected, 3)
        self.assertTrue(
            any("Skipped incomplete CO2_DREIDING molecule" in warning for warning in state.warnings)
        )
        state_dict = state.to_dict()
        self.assertEqual(state_dict["skipped_molecules"][0]["n_sites_found"], 2)

    def test_lammps_md_handoff_conserves_co2_population(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            snapshot_path = temp_path / "result_0.data"
            snapshot_path.write_text(
                "\n".join(
                    [
                        "gRASPA movie snapshot",
                        "",
                        "3 atoms",
                        "2 atom types",
                        "",
                        "Masses",
                        "",
                        "1 15.9994 # O_co2",
                        "2 12.0107 # C_co2",
                        "",
                        "Atoms # full",
                        "",
                        "1 1 1 -0.35 1.0 2.0 3.16 0 0 0 # CO2_DREIDING O_co2",
                        "2 1 2 0.70 1.0 2.0 2.0 0 0 0 # CO2_DREIDING C_co2",
                        "3 1 1 -0.35 1.0 2.0 0.84 0 0 0 # CO2_DREIDING O_co2",
                        "",
                    ]
                ),
                encoding="utf-8",
            )
            templates, sites = load_lammps_guest_force_field_assets(("CO2_DREIDING",))
            previous_state = parse_lammps_guest_restart_snapshot(snapshot_path, templates=templates, sites=sites)

            data_path = temp_path / "lammps_md_input.data"
            data_path.write_text(
                "\n".join(
                    [
                        "LAMMPS data",
                        "",
                        "6 atoms",
                        "4 atom types",
                        "",
                        "Masses",
                        "",
                        "1 12.011 # C_3",
                        "2 15.999 # O_3",
                        "3 15.9994 # O_co2",
                        "4 12.0107 # C_co2",
                        "",
                        "Atoms # full",
                        "",
                        "1 1 1 0.0 0.0 0.0 0.0 0 0 0 # C_3",
                        "2 1 2 0.0 1.0 0.0 0.0 0 0 0 # O_3",
                        "3 1 2 0.0 0.0 1.0 0.0 0 0 0 # O_3",
                        "4 2 3 -0.35 1.1 2.1 3.2 0 0 0 # O_co2 CO2_DREIDING",
                        "5 2 4 0.70 1.1 2.1 2.1 0 0 0 # C_co2 CO2_DREIDING",
                        "6 2 3 -0.35 1.1 2.1 0.9 0 0 0 # O_co2 CO2_DREIDING",
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
                        "6",
                        "ITEM: BOX BOUNDS pp pp pp",
                        "0 20",
                        "0 21",
                        "0 22",
                        "ITEM: ATOMS id type x y z",
                        "1 1 0.0 0.0 0.0",
                        "2 2 1.0 0.0 0.0",
                        "3 2 0.0 1.0 0.0",
                        "4 3 1.1 2.1 3.2",
                        "5 4 1.1 2.1 2.1",
                        "6 3 1.1 2.1 0.9",
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
                component_order=("CO2_DREIDING",),
            )
            restart_text = Path(restart_result.restart_file_path).read_text(encoding="utf-8")

        self.assertEqual(md_state.n_atoms, 3)
        self.assertEqual(md_state.components, ("CO2_DREIDING",))
        self.assertEqual([atom.site_label for atom in md_state.atoms], ["O_co2", "C_co2", "O_co2"])
        self.assertEqual(md_state.skipped_molecules, ())
        self.assertFalse(any("handoff recovered" in warning for warning in md_state.warnings))
        self.assertEqual(restart_result.n_adsorbate_atoms, 3)
        self.assertEqual(restart_result.n_adsorbate_molecules, 1)
        self.assertIn("Components: 1 (Adsorbates 1, Cations 0)", restart_text)
        self.assertIn("Component: 0   Adsorbate 1 molecules of CO2_DREIDING", restart_text)
        self.assertEqual(restart_text.count("Adsorbate-atom-position:"), 3)
        self.assertEqual(restart_result.warnings, ())

    def test_lammps_md_handoff_accounts_dropped_guest_atoms(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            snapshot_path = temp_path / "result_0.data"
            snapshot_path.write_text(
                "\n".join(
                    [
                        "gRASPA movie snapshot",
                        "",
                        "6 atoms",
                        "2 atom types",
                        "",
                        "Masses",
                        "",
                        "1 15.9994 # O_co2",
                        "2 12.0107 # C_co2",
                        "",
                        "Atoms # full",
                        "",
                        "1 1 1 -0.35 1.0 2.0 3.16 0 0 0 # CO2_DREIDING O_co2",
                        "2 1 2 0.70 1.0 2.0 2.0 0 0 0 # CO2_DREIDING C_co2",
                        "3 1 1 -0.35 1.0 2.0 0.84 0 0 0 # CO2_DREIDING O_co2",
                        "4 2 1 -0.35 5.0 6.0 7.16 0 0 0 # CO2_DREIDING O_co2",
                        "5 2 2 0.70 5.0 6.0 6.0 0 0 0 # CO2_DREIDING C_co2",
                        "6 2 1 -0.35 5.0 6.0 4.84 0 0 0 # CO2_DREIDING O_co2",
                        "",
                    ]
                ),
                encoding="utf-8",
            )
            templates, sites = load_lammps_guest_force_field_assets(("CO2_DREIDING",))
            previous_state = parse_lammps_guest_restart_snapshot(snapshot_path, templates=templates, sites=sites)

            # The second CO2 molecule lost a site in the MD data file: it must
            # be accounted, not silently dropped.
            data_path = temp_path / "lammps_md_input.data"
            data_path.write_text(
                "\n".join(
                    [
                        "LAMMPS data",
                        "",
                        "8 atoms",
                        "4 atom types",
                        "",
                        "Masses",
                        "",
                        "1 12.011 # C_3",
                        "2 15.999 # O_3",
                        "3 15.9994 # O_co2",
                        "4 12.0107 # C_co2",
                        "",
                        "Atoms # full",
                        "",
                        "1 1 1 0.0 0.0 0.0 0.0 0 0 0 # C_3",
                        "2 1 2 0.0 1.0 0.0 0.0 0 0 0 # O_3",
                        "3 1 2 0.0 0.0 1.0 0.0 0 0 0 # O_3",
                        "4 2 3 -0.35 1.1 2.1 3.2 0 0 0 # O_co2 CO2_DREIDING",
                        "5 2 4 0.70 1.1 2.1 2.1 0 0 0 # C_co2 CO2_DREIDING",
                        "6 2 3 -0.35 1.1 2.1 0.9 0 0 0 # O_co2 CO2_DREIDING",
                        "7 3 3 -0.35 5.1 6.1 7.2 0 0 0 # O_co2 CO2_DREIDING",
                        "8 3 4 0.70 5.1 6.1 6.1 0 0 0 # C_co2 CO2_DREIDING",
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
                        "8",
                        "ITEM: BOX BOUNDS pp pp pp",
                        "0 20",
                        "0 21",
                        "0 22",
                        "ITEM: ATOMS id type x y z",
                        "1 1 0.0 0.0 0.0",
                        "2 2 1.0 0.0 0.0",
                        "3 2 0.0 1.0 0.0",
                        "4 3 1.1 2.1 3.2",
                        "5 4 1.1 2.1 2.1",
                        "6 3 1.1 2.1 0.9",
                        "7 3 5.1 6.1 7.2",
                        "8 4 5.1 6.1 6.1",
                        "",
                    ]
                ),
                encoding="utf-8",
            )

            md_state, _cell = build_lammps_guest_restart_state_from_lammps_md_result(
                SimpleNamespace(lammps_data_path=str(data_path), lammps_dump_path=str(dump_path)),
                previous_guest_restart_state=previous_state,
            )

        self.assertEqual(md_state.n_atoms, 3)
        self.assertEqual(len(md_state.skipped_molecules), 1)
        skipped = md_state.skipped_molecules[0]
        self.assertEqual(skipped.molecule_key, "CO2_DREIDING:md:3")
        self.assertEqual(skipped.n_sites_found, 2)
        self.assertTrue(
            any("incomplete CO2_DREIDING molecule" in warning for warning in md_state.warnings)
        )
        self.assertTrue(
            any("handoff recovered 3 of 6 guest atoms" in warning for warning in md_state.warnings)
        )

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
            templates, sites = load_lammps_guest_force_field_assets(("Xe_GENERICMOFS", "Kr_GENERICMOFS"))
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
                component_order=("Xe_GENERICMOFS", "Kr_GENERICMOFS"),
            )
            restart_text = Path(restart_result.restart_file_path).read_text(encoding="utf-8")

        self.assertEqual(md_state.n_atoms, 2)
        self.assertEqual(md_state.components, ("Xe_GENERICMOFS", "Kr_GENERICMOFS"))
        self.assertEqual((md_state.atoms[0].x, md_state.atoms[0].y, md_state.atoms[0].z), (1.25, 2.50, 3.75))
        self.assertEqual(cell.lengths, (20.0, 21.0, 22.0))
        self.assertEqual(restart_result.n_adsorbate_molecules, 2)
        self.assertEqual(restart_result.components, ("Xe_GENERICMOFS", "Kr_GENERICMOFS"))
        self.assertIn("Components: 2 (Adsorbates 2, Cations 0)", restart_text)
        self.assertIn("Components 0 (Xe_GENERICMOFS)", restart_text)
        self.assertIn("Components 1 (Kr_GENERICMOFS)", restart_text)
        self.assertIn("Adsorbate-atom-position: 0 0 1.250000 2.500000 3.750000", restart_text)
        self.assertIn("Adsorbate-atom-position: 0 0 4.250000 5.500000 6.750000", restart_text)
        self.assertIn("Adsorbate-atom-charge: 0 0 0.000000", restart_text)
        self.assertIn("Adsorbate-atom-scaling: 0 0 1", restart_text)
        self.assertIn("Adsorbate-atom-fixed: 0 0 0  0  0", restart_text)

    def test_feynman_hibbs_guest_row_is_diagnosed_not_converted_silently(self):
        # Impact-review claim T4-20 / action A11: a massive external guest
        # parameterized as feynman-hibbs-lennard-jones is staged into LAMMPS
        # as plain classical LJ numbers; that conversion must be an explicit
        # diagnostic, never silent.
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            bundle_path = _write_guest_bundle(
                temp_path / "h2fh_bundle.json",
                name="H2FH_CUSTOM",
                site_label="H2_fh",
                element="H",
                mass=2.016,
                charge=0.0,
                mixing_rule_row="H2_fh feynman-hibbs-lennard-jones 36.7000 2.95800",
            )
            templates, sites = load_lammps_guest_force_field_assets(
                ("H2FH_CUSTOM",), guest_bundles=[str(bundle_path)]
            )
            snapshot_path = temp_path / "result_5.data"
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
                        "1 2.016 # H2_fh",
                        "",
                        "Atoms # full",
                        "",
                        "1 1 1 0.0 1.0 2.0 3.0 0 0 0",
                        "",
                    ]
                ),
                encoding="utf-8",
            )
            state = parse_lammps_guest_restart_snapshot(snapshot_path, templates=templates, sites=sites)

        self.assertEqual(sites[0].interaction, "feynman-hibbs-lennard-jones")
        self.assertEqual(sites[0].epsilon_k, 36.7)
        diagnostics = guest_site_model_diagnostics(sites)
        self.assertEqual(len(diagnostics), 1)
        self.assertIn("Feynman-Hibbs", diagnostics[0])
        self.assertIn("classical Lennard-Jones", diagnostics[0])
        self.assertIn("'H2_fh'", diagnostics[0])
        self.assertTrue(
            any("Feynman-Hibbs" in warning for warning in state.warnings),
            msg=f"state warnings did not carry the conversion diagnostic: {state.warnings}",
        )
        self.assertEqual(state.to_dict()["sites"][0]["interaction"], "feynman-hibbs-lennard-jones")

    def test_conflicting_guest_overrides_are_diagnosed(self):
        # Impact-review claim T4-24 / action A11: bundle lammps-section
        # overrides that disagree with the raspa-section rows staged for the
        # MC engine make the two hybrid legs run different guest models;
        # that conflict must be an explicit diagnostic.
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            bundle_path = _write_guest_bundle(
                temp_path / "xe_conflict_bundle.json",
                name="XE_CONFLICT",
                site_label="Xe_c",
                element="Xe",
                mass=131.293,
                charge=0.0,
                mixing_rule_row="Xe_c lennard-jones 221.0 4.01000",
                lammps_overrides={
                    "charges": {"Xe_c": -0.05},
                    "pair_coeff_rows": ["Xe_c lennard-jones 230.0 4.01000"],
                },
            )
            _templates, sites = load_lammps_guest_force_field_assets(
                ("XE_CONFLICT",), guest_bundles=[str(bundle_path)]
            )

        site = sites[0]
        self.assertEqual(site.charge, -0.05)
        self.assertEqual(site.epsilon_k, 230.0)
        self.assertEqual(len(site.override_conflicts), 2)
        self.assertTrue(any("charge override -0.05 differs from the RASPA-side" in c for c in site.override_conflicts))
        self.assertTrue(any("epsilon_k=230.0" in c and "epsilon_k=221.0" in c for c in site.override_conflicts))
        diagnostics = guest_site_model_diagnostics(sites)
        self.assertEqual(len(diagnostics), 2)
        self.assertTrue(all("different guest model values" in d for d in diagnostics))

    def test_matching_guest_overrides_produce_no_conflict_diagnostics(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            bundle_path = _write_guest_bundle(
                temp_path / "xe_ok_bundle.json",
                name="XE_OK",
                site_label="Xe_o",
                element="Xe",
                mass=131.293,
                charge=0.0,
                mixing_rule_row="Xe_o lennard-jones 221.0 4.01000",
                lammps_overrides={
                    "masses": {"Xe_o": 131.293},
                    "charges": {"Xe_o": 0.0},
                    "pair_coeff_rows": ["Xe_o lennard-jones 221.0 4.01000"],
                },
            )
            _templates, sites = load_lammps_guest_force_field_assets(
                ("XE_OK",), guest_bundles=[str(bundle_path)]
            )

        self.assertEqual(sites[0].override_conflicts, ())
        self.assertEqual(guest_site_model_diagnostics(sites), ())


def _write_guest_bundle(
    path: Path,
    *,
    name: str,
    site_label: str,
    element: str,
    mass: float,
    charge: float,
    mixing_rule_row: str,
    lammps_overrides: dict | None = None,
) -> Path:
    path.write_text(
        json.dumps(
            {
                "version": 2,
                "name": name,
                "parameter_family": "dreiding",
                "parameter_source": "test guest model",
                "raspa": {
                    "molecule_definition": (
                        "# critical constants: Temperature [T], Pressure [Pa], and Acentric factor [-]\n"
                        "100.0\n"
                        "1000000.0\n"
                        "0.0\n"
                        "#Number Of Atoms\n"
                        " 1\n"
                        "# Number of groups\n"
                        "1\n"
                        "# test-group\n"
                        "rigid\n"
                        "# number of atoms\n"
                        "1\n"
                        "# atomic positions\n"
                        f"0 {site_label}     0.0 0.0 0.0\n"
                        "# Chiral centers Bond  BondDipoles Bend  UrayBradley InvBend  Torsion Imp. Torsion Bond/Bond Stretch/Bend Bend/Bend Stretch/Torsion Bend/Torsion IntraVDW IntraCoulomb\n"
                        "               0    0            0    0            0       0        0            0         0            0         0               0            0        0            0\n"
                        "# Number of config moves\n"
                        "0\n"
                    ),
                    "pseudo_atom_rows": [
                        f"{site_label}      yes     {element}     {element}     0          {mass}    {charge}      0.0          1.0      0.720  0            0           relative           0",
                    ],
                    "mixing_rule_rows": [mixing_rule_row],
                },
                "lammps": {"units": "real", "atom_style": "full", **(lammps_overrides or {})},
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    return path


if __name__ == "__main__":
    unittest.main()
