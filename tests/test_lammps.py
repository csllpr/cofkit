import contextlib
import io
import json
import os
import stat
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from cofkit.cli import main as cli_main
import cofkit.lammps as lammps_module
from cofkit.lammps import (
    COFKIT_LMP_ENV_VAR,
    LammpsConfigurationError,
    LammpsInputError,
    LammpsOptimizationSettings,
    optimize_cif_with_lammps,
    resolve_lammps_binary,
)


class LammpsTests(unittest.TestCase):
    def test_resolve_lammps_binary_requires_configuration_without_env_or_default(self):
        with patch.object(lammps_module, "DEFAULT_LAMMPS_BINARY", Path("/definitely/missing/lmp_mpi")):
            with patch.dict(os.environ, {}, clear=True):
                with self.assertRaises(LammpsConfigurationError) as raised:
                    resolve_lammps_binary()

        self.assertIn(COFKIT_LMP_ENV_VAR, str(raised.exception))

    def test_optimize_cif_with_lammps_writes_report_and_optimized_cif(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            fake_binary = self._write_fake_lammps_binary(temp_path / "lmp_fake")
            cif_path = temp_path / "example.cif"
            cif_path.write_text(self._example_cif_text(), encoding="utf-8")
            output_dir = temp_path / "lammps_out"

            with patch.dict(os.environ, {COFKIT_LMP_ENV_VAR: str(fake_binary)}):
                result = optimize_cif_with_lammps(
                    cif_path,
                    output_dir=output_dir,
                    settings=LammpsOptimizationSettings(forcefield="uff"),
                )

            self.assertTrue(Path(result.optimized_cif).is_file())
            self.assertTrue(Path(result.report_path).is_file())
            self.assertEqual(result.n_atoms, 3)
            self.assertEqual(result.n_bonds, 2)
            self.assertEqual(result.n_angles, 1)
            self.assertEqual(result.n_dihedrals, 0)
            self.assertEqual(result.n_impropers, 0)
            self.assertEqual(result.forcefield_backend, "uff_openbabel_explicit_graph_pymatgen")
            self.assertIn("WARNING: fake warning from test binary", result.warnings)

            optimized_text = Path(result.optimized_cif).read_text(encoding="utf-8")
            self.assertIn("_ccdc_geom_bond_type", optimized_text)
            self.assertIn("a2 C 0.200000 0.125000 0.100000 1.00", optimized_text)
            self.assertIn("a1 a2 . . 1.030776 D", optimized_text)

            report = json.loads(Path(result.report_path).read_text(encoding="utf-8"))
            self.assertEqual(report["n_angles"], 1)
            self.assertEqual(report["n_dihedrals"], 0)
            self.assertEqual(report["settings"]["forcefield"], "uff")
            self.assertEqual(report["n_atom_types"], 2)
            self.assertTrue(report["atom_type_symbols"]["1"])

    def test_optimize_cif_with_lammps_defaults_to_uff(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            fake_binary = self._write_fake_lammps_binary(temp_path / "lmp_fake")
            cif_path = temp_path / "example.cif"
            cif_path.write_text(self._example_cif_text(), encoding="utf-8")

            with patch.dict(os.environ, {COFKIT_LMP_ENV_VAR: str(fake_binary)}):
                result = optimize_cif_with_lammps(cif_path, output_dir=temp_path / "default_uff_out")
            script_text = Path(result.lammps_input_script_path).read_text(encoding="utf-8")

        self.assertEqual(result.settings.forcefield, "uff")
        self.assertEqual(result.forcefield_backend, "uff_openbabel_explicit_graph_pymatgen")
        self.assertEqual(result.settings.pre_minimization_steps, 10000)
        self.assertTrue(result.settings.two_stage_protocol)
        self.assertTrue(result.settings.relax_cell)
        self.assertEqual(result.settings.max_iterations, 200000)
        self.assertEqual(result.settings.max_evaluations, 2000000)
        self.assertIn("# pre_minimization", script_text)
        self.assertIn("run 10000", script_text)
        self.assertIn("# stage2", script_text)
        self.assertIn("# box_relax", script_text)
        self.assertEqual(script_text.count("fix cofkit_hold all spring/self"), 1)
        self.assertEqual(script_text.count("minimize "), 3)
        self.assertIn("fix cofkit_boxrelax all box/relax aniso 0 vmax 0.001", script_text)

    def test_lammps_script_includes_fix_modify_for_restraint_energy(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            fake_binary = self._write_fake_lammps_binary(temp_path / "lmp_fake")
            cif_path = temp_path / "example.cif"
            cif_path.write_text(self._example_cif_text(), encoding="utf-8")

            with patch.dict(os.environ, {COFKIT_LMP_ENV_VAR: str(fake_binary)}):
                result = optimize_cif_with_lammps(
                    cif_path,
                    output_dir=temp_path / "fix_modify_out",
                    settings=LammpsOptimizationSettings(forcefield="uff"),
                )

            script_text = Path(result.lammps_input_script_path).read_text(encoding="utf-8")
            self.assertIn("fix cofkit_hold all spring/self 0.200000", script_text)
            self.assertIn("fix_modify cofkit_hold energy yes", script_text)

    def test_lammps_script_supports_two_stage_protocol_and_min_modify(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            fake_binary = self._write_fake_lammps_binary(temp_path / "lmp_fake")
            cif_path = temp_path / "example.cif"
            cif_path.write_text(self._example_cif_text(), encoding="utf-8")
            settings = LammpsOptimizationSettings(
                forcefield="uff",
                two_stage_protocol=True,
                stage2_position_restraint_force_constant=0.05,
                stage2_min_style="cg",
                relax_cell=False,
                timestep=0.5,
                min_modify_dmax=0.15,
                min_modify_line="quadratic",
                min_modify_norm="max",
                min_modify_fire_integrator="verlet",
                min_modify_fire_tmax=4.0,
                min_modify_fire_abcfire=True,
            )

            with patch.dict(os.environ, {COFKIT_LMP_ENV_VAR: str(fake_binary)}):
                result = optimize_cif_with_lammps(
                    cif_path,
                    output_dir=temp_path / "two_stage_out",
                    settings=settings,
                )

            script_text = Path(result.lammps_input_script_path).read_text(encoding="utf-8")
            self.assertIn("timestep 0.5", script_text)
            self.assertIn("# stage1", script_text)
            self.assertIn("# stage2", script_text)
            self.assertIn("min_style fire", script_text)
            self.assertIn("min_style cg", script_text)
            self.assertIn("fix cofkit_hold all spring/self 0.200000", script_text)
            self.assertIn("fix cofkit_hold all spring/self 0.050000", script_text)
            self.assertIn("min_modify dmax 0.15 line quadratic norm max integrator verlet tmax 4 abcfire yes", script_text)
            self.assertEqual(script_text.count("minimize "), 2)

    def test_lammps_script_supports_optional_pre_minimization_prerun(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            fake_binary = self._write_fake_lammps_binary(temp_path / "lmp_fake")
            cif_path = temp_path / "example.cif"
            cif_path.write_text(self._example_cif_text(), encoding="utf-8")
            settings = LammpsOptimizationSettings(
                forcefield="uff",
                pre_minimization_steps=25,
                pre_minimization_temperature=350.0,
                pre_minimization_damping=50.0,
                pre_minimization_seed=13579,
                pre_minimization_displacement_limit=0.05,
            )

            with patch.dict(os.environ, {COFKIT_LMP_ENV_VAR: str(fake_binary)}):
                result = optimize_cif_with_lammps(
                    cif_path,
                    output_dir=temp_path / "prerun_out",
                    settings=settings,
                )

            script_text = Path(result.lammps_input_script_path).read_text(encoding="utf-8")
            self.assertIn("# pre_minimization", script_text)
            self.assertIn("velocity all create 350 13579 mom yes rot yes dist gaussian", script_text)
            self.assertIn("fix cofkit_prerun_nve all nve/limit 0.05", script_text)
            self.assertIn("fix cofkit_prerun_langevin all langevin 350 350 50 118308", script_text)
            self.assertIn("run 25", script_text)
            self.assertIn("unfix cofkit_prerun_langevin", script_text)
            self.assertIn("unfix cofkit_prerun_nve", script_text)

    def test_lammps_script_supports_final_box_relax_stage(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            fake_binary = self._write_fake_lammps_binary(temp_path / "lmp_fake")
            cif_path = temp_path / "example.cif"
            cif_path.write_text(self._example_cif_text(), encoding="utf-8")
            settings = LammpsOptimizationSettings(
                forcefield="uff",
                two_stage_protocol=True,
                stage2_position_restraint_force_constant=0.0,
                relax_cell=True,
                box_relax_mode="auto",
                box_relax_vmax=0.002,
                box_relax_nreset=25,
                box_relax_min_style="cg",
            )

            with patch.dict(os.environ, {COFKIT_LMP_ENV_VAR: str(fake_binary)}):
                result = optimize_cif_with_lammps(
                    cif_path,
                    output_dir=temp_path / "box_relax_out",
                    settings=settings,
                )

            script_text = Path(result.lammps_input_script_path).read_text(encoding="utf-8")
            self.assertIn("# box_relax", script_text)
            self.assertIn("fix cofkit_boxrelax all box/relax aniso 0 vmax 0.002 nreset 25", script_text)
            self.assertIn("min_style cg", script_text)
            self.assertIn("unfix cofkit_boxrelax", script_text)
            self.assertEqual(script_text.count("minimize "), 3)

    def test_relax_cell_rejects_fire_box_relax_minimizer(self):
        with self.assertRaises(ValueError):
            optimize_cif_with_lammps(
                Path("/definitely/missing.cif"),
                settings=LammpsOptimizationSettings(
                    forcefield="uff",
                    relax_cell=True,
                    box_relax_min_style="fire",
                ),
            )

    def test_resolve_lammps_binary_accepts_executable_name_on_path(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            fake_binary = self._write_fake_lammps_binary(temp_path / "lmp_fake")
            with patch.object(lammps_module, "DEFAULT_LAMMPS_BINARY", Path("/definitely/missing/lmp_mpi")):
                with patch.dict(
                    os.environ,
                    {
                        COFKIT_LMP_ENV_VAR: "lmp_fake",
                        "PATH": f"{temp_path}{os.pathsep}{os.environ.get('PATH', '')}",
                    },
                    clear=False,
                ):
                    resolved = resolve_lammps_binary()

            self.assertEqual(resolved, fake_binary.resolve())

    def test_lammps_defaults_omp_num_threads_to_half_core_count(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            fake_binary = self._write_fake_lammps_binary(temp_path / "lmp_fake")
            cif_path = temp_path / "example.cif"
            cif_path.write_text(self._example_cif_text(), encoding="utf-8")
            output_dir = temp_path / "omp_default_out"

            with patch.object(lammps_module.os, "cpu_count", return_value=10):
                with patch.dict(os.environ, {COFKIT_LMP_ENV_VAR: str(fake_binary)}, clear=False):
                    result = optimize_cif_with_lammps(
                        cif_path,
                        output_dir=output_dir,
                        settings=LammpsOptimizationSettings(forcefield="uff"),
                    )

            stdout_text = Path(result.stdout_log_path).read_text(encoding="utf-8")
            self.assertIn("OMP_NUM_THREADS=5", stdout_text)

    def test_lammps_preserves_explicit_omp_num_threads(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            fake_binary = self._write_fake_lammps_binary(temp_path / "lmp_fake")
            cif_path = temp_path / "example.cif"
            cif_path.write_text(self._example_cif_text(), encoding="utf-8")
            output_dir = temp_path / "omp_override_out"

            with patch.object(lammps_module.os, "cpu_count", return_value=10):
                with patch.dict(
                    os.environ,
                    {COFKIT_LMP_ENV_VAR: str(fake_binary), "OMP_NUM_THREADS": "7"},
                    clear=False,
                ):
                    result = optimize_cif_with_lammps(
                        cif_path,
                        output_dir=output_dir,
                        settings=LammpsOptimizationSettings(forcefield="uff"),
                    )

            stdout_text = Path(result.stdout_log_path).read_text(encoding="utf-8")
            self.assertIn("OMP_NUM_THREADS=7", stdout_text)

    def test_calculate_lammps_optimize_cli_prints_json_report(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            fake_binary = self._write_fake_lammps_binary(temp_path / "lmp_fake")
            cif_path = temp_path / "example.cif"
            cif_path.write_text(self._example_cif_text(), encoding="utf-8")
            output_dir = temp_path / "cli_out"
            buffer = io.StringIO()

            with patch.dict(os.environ, {COFKIT_LMP_ENV_VAR: str(fake_binary)}):
                with contextlib.redirect_stdout(buffer):
                    cli_main(
                        [
                            "calculate",
                            "lammps-optimize",
                            str(cif_path),
                            "--output-dir",
                            str(output_dir),
                            "--forcefield",
                            "uff",
                            "--json",
                        ]
                    )

            report = json.loads(buffer.getvalue())
            self.assertEqual(report["output_dir"], str(output_dir.resolve()))
            self.assertEqual(report["n_bonds"], 2)
            self.assertEqual(report["n_angles"], 1)
            self.assertEqual(report["settings"]["forcefield"], "uff")
            self.assertTrue(report["settings"]["two_stage_protocol"])
            self.assertTrue(report["settings"]["relax_cell"])

    def test_output_cif_writes_inferred_periodic_bond_symmetry_when_needed(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            fake_binary = self._write_fake_lammps_binary(temp_path / "lmp_fake")
            cif_path = temp_path / "periodic_example.cif"
            cif_path.write_text(self._periodic_example_cif_text(), encoding="utf-8")

            with patch.dict(os.environ, {COFKIT_LMP_ENV_VAR: str(fake_binary)}):
                result = optimize_cif_with_lammps(
                    cif_path,
                    output_dir=temp_path / "periodic_out",
                    settings=LammpsOptimizationSettings(forcefield="uff"),
                )

            optimized_text = Path(result.optimized_cif).read_text(encoding="utf-8")
            self.assertIn("a1 a2 . 1_655 1.030776 S", optimized_text)

    def test_output_cif_recomputes_explicit_periodic_bond_symmetry_from_final_geometry(self):
        cif_path = Path(__file__).resolve().parent / "fixtures" / "tabf3__tctb__hcb_periodic_conflict.cif"
        parsed = lammps_module._parse_explicit_bond_cif(cif_path)
        final_positions = {atom.atom_id: atom.fractional_position for atom in parsed.atoms}
        atom_id_by_label = {atom.label: atom.atom_id for atom in parsed.atoms}
        final_positions[atom_id_by_label["m2_C14"]] = (0.520042, 0.014438, 0.415150)
        final_positions[atom_id_by_label["m1_N6"]] = (0.408258, 0.104422, 0.485174)

        rendered = lammps_module._render_optimized_cif(parsed, final_positions)
        self.assertIn("m2_C14 m1_N6 . .", rendered)
        self.assertNotIn("m2_C14 m1_N6 . 1_565", rendered)

        with tempfile.TemporaryDirectory() as temp_dir:
            rendered_path = Path(temp_dir) / "rendered.cif"
            rendered_path.write_text(rendered, encoding="utf-8")
            reparsed = lammps_module._parse_explicit_bond_cif(rendered_path)

        for bond in reparsed.bonds:
            if {bond.label_1, bond.label_2} != {"m2_C14", "m1_N6"}:
                continue
            distance = lammps_module._bond_distance(
                reparsed.atoms[bond.atom_id_1 - 1].fractional_position,
                reparsed.atoms[bond.atom_id_2 - 1].fractional_position,
                bond.shift_1,
                bond.shift_2,
                reparsed.lammps_basis,
            )
            self.assertLess(distance, 2.0)
            break
        else:
            self.fail("Rendered CIF did not preserve the m2_C14-m1_N6 bond.")

    def test_output_cif_uses_relaxed_dump_cell_and_origin(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            cif_path = temp_path / "example.cif"
            cif_path.write_text(self._example_cif_text(), encoding="utf-8")
            parsed = lammps_module._parse_explicit_bond_cif(cif_path)
            basis = lammps_module._lammps_basis_from_cell(8.0, 9.0, 10.0, 90.0, 90.0, 60.0)
            origin = (0.5, 1.0, 2.0)
            a_vec, b_vec, c_vec = basis
            xy = b_vec[0]
            xz = c_vec[0]
            yz = c_vec[1]
            xlo = origin[0]
            xhi = origin[0] + a_vec[0]
            ylo = origin[1]
            yhi = origin[1] + a_vec[1] + b_vec[1]
            zlo = origin[2]
            zhi = origin[2] + c_vec[2]
            xlo_bound = xlo + min(0.0, xy, xz, xy + xz)
            xhi_bound = xhi + max(0.0, xy, xz, xy + xz)
            ylo_bound = ylo + min(0.0, yz)
            yhi_bound = yhi + max(0.0, yz)

            atom_lines: list[str] = []
            for atom in parsed.atoms:
                x, y, z = lammps_module._fractional_to_lammps(atom.fractional_position, basis)
                atom_lines.append(
                    f"{atom.atom_id} {x + origin[0]:.6f} {y + origin[1]:.6f} {z + origin[2]:.6f}"
                )

            dump_path = temp_path / "relaxed.lammpstrj"
            dump_path.write_text(
                "\n".join(
                    [
                        "ITEM: TIMESTEP",
                        "0",
                        "ITEM: NUMBER OF ATOMS",
                        str(len(parsed.atoms)),
                        "ITEM: BOX BOUNDS xy xz yz pp pp pp",
                        f"{xlo_bound:.6f} {xhi_bound:.6f} {xy:.6f}",
                        f"{ylo_bound:.6f} {yhi_bound:.6f} {xz:.6f}",
                        f"{zlo:.6f} {zhi:.6f} {yz:.6f}",
                        "ITEM: ATOMS id x y z",
                        *atom_lines,
                        "",
                    ]
                ),
                encoding="utf-8",
            )

            frame = lammps_module._parse_lammps_dump_last_frame(dump_path, expected_atoms=len(parsed.atoms))
            final_positions = lammps_module._cartesian_positions_to_fractional(
                frame.lammps_origin,
                frame.lammps_basis,
                frame.cartesian_positions,
            )
            rendered = lammps_module._render_optimized_cif(
                parsed,
                final_positions,
                cell_parameters=frame.cell_parameters,
                basis=frame.lammps_basis,
            )

            self.assertIn("_cell_length_a 8.000000", rendered)
            self.assertIn("_cell_length_b 9.000000", rendered)
            self.assertIn("_cell_length_c 10.000000", rendered)
            self.assertIn("a2 C 0.200000 0.100000 0.100000 1.00", rendered)
            self.assertIn("a1 a2 . . 0.800000 D", rendered)

            rendered_path = temp_path / "rendered.cif"
            rendered_path.write_text(rendered, encoding="utf-8")
            reparsed = lammps_module._parse_explicit_bond_cif(rendered_path)
            self.assertAlmostEqual(reparsed.cell_parameters[0], 8.0, places=6)
            self.assertAlmostEqual(reparsed.cell_parameters[1], 9.0, places=6)
            self.assertAlmostEqual(reparsed.cell_parameters[2], 10.0, places=6)
            self.assertAlmostEqual(reparsed.cell_parameters[5], 60.0, places=5)

    def test_lammps_keeps_periodic_primitive_model_for_conflicted_primitive_bond_graph(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            fake_binary = self._write_fake_lammps_binary(temp_path / "lmp_fake")
            cif_path = temp_path / "conflicted_ring.cif"
            cif_path.write_text(self._conflicted_ring_cif_text(), encoding="utf-8")

            with patch.dict(os.environ, {COFKIT_LMP_ENV_VAR: str(fake_binary)}):
                result = optimize_cif_with_lammps(
                    cif_path,
                    output_dir=temp_path / "conflicted_out",
                    settings=LammpsOptimizationSettings(forcefield="uff"),
                )

            self.assertEqual(result.n_atoms, 4)
            self.assertNotIn("replicated finite-cover cluster", " ".join(result.warnings))

            script_text = Path(result.lammps_input_script_path).read_text(encoding="utf-8")
            self.assertIn("boundary p p p", script_text)

            optimized_text = Path(result.optimized_cif).read_text(encoding="utf-8")
            self.assertIn("a1 C", optimized_text)
            self.assertIn("a4 C", optimized_text)

    def test_conflicted_periodic_primitive_graph_preserves_requested_prerun(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            fake_binary = self._write_fake_lammps_binary(temp_path / "lmp_fake")
            cif_path = temp_path / "conflicted_ring.cif"
            cif_path.write_text(self._conflicted_ring_cif_text(), encoding="utf-8")

            with patch.dict(os.environ, {COFKIT_LMP_ENV_VAR: str(fake_binary)}):
                result = optimize_cif_with_lammps(
                    cif_path,
                    output_dir=temp_path / "conflicted_prerun_out",
                    settings=LammpsOptimizationSettings(
                        forcefield="uff",
                        pre_minimization_steps=5,
                        timestep=0.5,
                    ),
                )

            script_text = Path(result.lammps_input_script_path).read_text(encoding="utf-8")
            self.assertIn("boundary p p p", script_text)
            self.assertIn("# pre_minimization", script_text)
            self.assertIn("run 5", script_text)
            self.assertEqual(result.settings.pre_minimization_steps, 5)
            self.assertNotIn("Skipped pre-minimization MD", " ".join(result.warnings))

    def test_conflicted_periodic_fixture_stays_primitive_periodic_model(self):
        cif_path = Path(__file__).resolve().parent / "fixtures" / "tabf3__tctb__hcb_periodic_conflict.cif"
        parsed = lammps_module._parse_explicit_bond_cif(cif_path)
        _atom_images, image_conflicts = lammps_module._compute_unwrapped_atom_images(parsed)
        self.assertGreater(image_conflicts, 0)
        model = lammps_module._build_optimization_model(
            parsed,
            settings=LammpsOptimizationSettings(forcefield="uff"),
        )

        self.assertEqual(model.boundary, "p p p")
        self.assertIsNone(model.primitive_parsed)
        self.assertEqual(model.parsed, parsed)
        self.assertEqual(len(model.output_atom_id_by_model_atom_id), len(parsed.atoms))

    def test_lammps_rejects_legacy_cif_without_explicit_bond_type(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            fake_binary = self._write_fake_lammps_binary(temp_path / "lmp_fake")
            cif_path = temp_path / "legacy_example.cif"
            cif_path.write_text(self._example_cif_text(include_bond_type=False), encoding="utf-8")

            with patch.dict(os.environ, {COFKIT_LMP_ENV_VAR: str(fake_binary)}):
                with self.assertRaises(LammpsInputError) as raised:
                    optimize_cif_with_lammps(
                        cif_path,
                        output_dir=temp_path / "legacy_out",
                        settings=LammpsOptimizationSettings(forcefield="uff"),
                    )

            self.assertIn("_ccdc_geom_bond_type", str(raised.exception))

    def test_lammps_data_can_include_dihedrals_and_impropers_from_explicit_cif(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            fake_binary = self._write_fake_lammps_binary(temp_path / "lmp_fake")
            cif_path = temp_path / "multi_term_example.cif"
            cif_path.write_text(self._multi_term_cif_text(), encoding="utf-8")

            with patch.dict(os.environ, {COFKIT_LMP_ENV_VAR: str(fake_binary)}):
                result = optimize_cif_with_lammps(
                    cif_path,
                    output_dir=temp_path / "multi_term_out",
                    settings=LammpsOptimizationSettings(forcefield="uff"),
                )

            script_text = Path(result.lammps_input_script_path).read_text(encoding="utf-8")
            data_text = Path(result.lammps_data_path).read_text(encoding="utf-8")
            self.assertIn("angle_style hybrid cosine/periodic", script_text)
            self.assertIn("dihedral_style harmonic", script_text)
            self.assertIn("improper_style fourier", script_text)
            self.assertIn("Dihedral Coeffs", data_text)
            self.assertIn("Improper Coeffs", data_text)
            self.assertIn("Dihedrals", data_text)
            self.assertIn("Impropers", data_text)
            self.assertGreaterEqual(result.n_dihedrals, 1)
            self.assertGreaterEqual(result.n_impropers, 1)

    def test_calculate_help_lists_lammps_optimize(self):
        buffer = io.StringIO()
        with self.assertRaises(SystemExit), contextlib.redirect_stdout(buffer):
            cli_main(["calculate", "--help"])

        self.assertIn("lammps-optimize", buffer.getvalue())

    def _write_fake_lammps_binary(self, path: Path) -> Path:
        path.write_text(
            "#!/usr/bin/env python\n"
            "from __future__ import annotations\n"
            "import sys\n"
            "import os\n"
            "from pathlib import Path\n"
            "\n"
            "args = sys.argv[1:]\n"
            "input_script = None\n"
            "log_path = None\n"
            "index = 0\n"
            "while index < len(args):\n"
            "    if args[index] == '-in':\n"
            "        input_script = Path(args[index + 1])\n"
            "        index += 2\n"
            "    elif args[index] == '-log':\n"
            "        log_path = Path(args[index + 1])\n"
            "        index += 2\n"
            "    else:\n"
            "        index += 1\n"
            "\n"
            "if input_script is None:\n"
            "    sys.stderr.write('missing -in\\n')\n"
            "    sys.exit(2)\n"
            "\n"
            "script_lines = input_script.read_text(encoding='utf-8').splitlines()\n"
            "data_path = None\n"
            "dump_path = None\n"
            "for line in script_lines:\n"
            "    parts = line.split()\n"
            "    if not parts:\n"
            "        continue\n"
            "    if parts[0] == 'read_data':\n"
            "        data_path = Path(parts[1])\n"
            "    elif parts[0] == 'dump':\n"
            "        dump_path = Path(parts[5])\n"
            "\n"
            "if data_path is None or dump_path is None:\n"
            "    sys.stderr.write('missing data or dump path\\n')\n"
            "    sys.exit(3)\n"
            "\n"
            "atoms = []\n"
            "in_atoms = False\n"
            "for raw in data_path.read_text(encoding='utf-8').splitlines():\n"
            "    line = raw.strip()\n"
            "    if line.startswith('Atoms'):\n"
            "        in_atoms = True\n"
            "        continue\n"
            "    if not in_atoms:\n"
            "        continue\n"
            "    if not line:\n"
            "        continue\n"
            "    if line.startswith('Bonds'):\n"
            "        break\n"
            "    parts = line.split()\n"
            "    if len(parts) >= 6:\n"
            "        atom_id = int(parts[0])\n"
            "        x = float(parts[3])\n"
            "        y = float(parts[4])\n"
            "        z = float(parts[5])\n"
            "        if atom_id == 2:\n"
            "            y += 0.25\n"
            "        atoms.append((atom_id, x, y, z))\n"
            "\n"
            "dump_path.parent.mkdir(parents=True, exist_ok=True)\n"
            "dump_path.write_text(\n"
            "    'ITEM: TIMESTEP\\n0\\n'\n"
            "    f'ITEM: NUMBER OF ATOMS\\n{len(atoms)}\\n'\n"
            "    'ITEM: BOX BOUNDS xy xz yz pp pp pp\\n'\n"
            "    '0 10 0\\n0 10 0\\n0 10 0\\n'\n"
            "    'ITEM: ATOMS id x y z\\n'\n"
            "    + ''.join(f'{atom_id} {x:.6f} {y:.6f} {z:.6f}\\n' for atom_id, x, y, z in atoms),\n"
            "    encoding='utf-8',\n"
            ")\n"
            "if log_path is not None:\n"
            "    log_path.write_text('LAMMPS fake log\\n', encoding='utf-8')\n"
            "sys.stdout.write(f\"OMP_NUM_THREADS={os.environ.get('OMP_NUM_THREADS', '')}\\n\")\n"
            "sys.stderr.write('WARNING: fake warning from test binary\\n')\n"
            "sys.exit(0)\n",
            encoding="utf-8",
        )
        path.chmod(path.stat().st_mode | stat.S_IXUSR)
        return path

    def _example_cif_text(self, *, include_bond_type: bool = True) -> str:
        bond_type_header = "_ccdc_geom_bond_type\n" if include_bond_type else ""
        return (
            "data_example\n"
            "_audit_creation_method 'cofkit test'\n"
            "_space_group_name_H-M_alt 'P 1'\n"
            "_space_group_IT_number 1\n"
            "_cell_length_a 10.000000\n"
            "_cell_length_b 10.000000\n"
            "_cell_length_c 10.000000\n"
            "_cell_angle_alpha 90.000000\n"
            "_cell_angle_beta 90.000000\n"
            "_cell_angle_gamma 90.000000\n"
            "\n"
            "loop_\n"
            "_space_group_symop_operation_xyz\n"
            "'x,y,z'\n"
            "\n"
            "loop_\n"
            "_atom_site_label\n"
            "_atom_site_type_symbol\n"
            "_atom_site_fract_x\n"
            "_atom_site_fract_y\n"
            "_atom_site_fract_z\n"
            "_atom_site_occupancy\n"
            "a1 C 0.100000 0.100000 0.100000 1.00\n"
            "a2 C 0.200000 0.100000 0.100000 1.00\n"
            "a3 O 0.300000 0.100000 0.100000 1.00\n"
            "\n"
            "loop_\n"
            "_geom_bond_atom_site_label_1\n"
            "_geom_bond_atom_site_label_2\n"
            "_geom_bond_site_symmetry_1\n"
            "_geom_bond_site_symmetry_2\n"
            "_geom_bond_distance\n"
            f"{bond_type_header}"
            f"a1 a2 . . 1.000000{' D' if include_bond_type else ''}\n"
            f"a2 a3 . . 1.000000{' S' if include_bond_type else ''}\n"
        )

    def _periodic_example_cif_text(self) -> str:
        return (
            "data_periodic_example\n"
            "_audit_creation_method 'cofkit test'\n"
            "_space_group_name_H-M_alt 'P 1'\n"
            "_space_group_IT_number 1\n"
            "_cell_length_a 10.000000\n"
            "_cell_length_b 10.000000\n"
            "_cell_length_c 10.000000\n"
            "_cell_angle_alpha 90.000000\n"
            "_cell_angle_beta 90.000000\n"
            "_cell_angle_gamma 90.000000\n"
            "\n"
            "loop_\n"
            "_space_group_symop_operation_xyz\n"
            "'x,y,z'\n"
            "\n"
            "loop_\n"
            "_atom_site_label\n"
            "_atom_site_type_symbol\n"
            "_atom_site_fract_x\n"
            "_atom_site_fract_y\n"
            "_atom_site_fract_z\n"
            "_atom_site_occupancy\n"
            "a1 C 0.950000 0.100000 0.100000 1.00\n"
            "a2 C 0.050000 0.100000 0.100000 1.00\n"
            "\n"
            "loop_\n"
            "_geom_bond_atom_site_label_1\n"
            "_geom_bond_atom_site_label_2\n"
            "_geom_bond_site_symmetry_1\n"
            "_geom_bond_site_symmetry_2\n"
            "_geom_bond_distance\n"
            "_ccdc_geom_bond_type\n"
            "a1 a2 . . 1.000000 S\n"
        )

    def _conflicted_ring_cif_text(self) -> str:
        return (
            "data_conflicted_ring\n"
            "_audit_creation_method 'cofkit test'\n"
            "_space_group_name_H-M_alt 'P 1'\n"
            "_space_group_IT_number 1\n"
            "_cell_length_a 10.000000\n"
            "_cell_length_b 10.000000\n"
            "_cell_length_c 10.000000\n"
            "_cell_angle_alpha 90.000000\n"
            "_cell_angle_beta 90.000000\n"
            "_cell_angle_gamma 90.000000\n"
            "\n"
            "loop_\n"
            "_space_group_symop_operation_xyz\n"
            "'x,y,z'\n"
            "\n"
            "loop_\n"
            "_atom_site_label\n"
            "_atom_site_type_symbol\n"
            "_atom_site_fract_x\n"
            "_atom_site_fract_y\n"
            "_atom_site_fract_z\n"
            "_atom_site_occupancy\n"
            "a1 C 0.900000 0.100000 0.100000 1.00\n"
            "a2 C 0.050000 0.100000 0.100000 1.00\n"
            "a3 C 0.050000 0.250000 0.100000 1.00\n"
            "a4 C 0.900000 0.250000 0.100000 1.00\n"
            "\n"
            "loop_\n"
            "_geom_bond_atom_site_label_1\n"
            "_geom_bond_atom_site_label_2\n"
            "_geom_bond_site_symmetry_1\n"
            "_geom_bond_site_symmetry_2\n"
            "_geom_bond_distance\n"
            "_ccdc_geom_bond_type\n"
            "a1 a2 . 1_655 1.500000 S\n"
            "a2 a3 . . 1.500000 S\n"
            "a3 a4 . 1_655 1.500000 S\n"
            "a4 a1 . . 1.500000 S\n"
        )

    def _multi_term_cif_text(self) -> str:
        return (
            "data_multi_term_example\n"
            "_audit_creation_method 'cofkit test'\n"
            "_space_group_name_H-M_alt 'P 1'\n"
            "_space_group_IT_number 1\n"
            "_cell_length_a 20.000000\n"
            "_cell_length_b 20.000000\n"
            "_cell_length_c 20.000000\n"
            "_cell_angle_alpha 90.000000\n"
            "_cell_angle_beta 90.000000\n"
            "_cell_angle_gamma 90.000000\n"
            "\n"
            "loop_\n"
            "_space_group_symop_operation_xyz\n"
            "'x,y,z'\n"
            "\n"
            "loop_\n"
            "_atom_site_label\n"
            "_atom_site_type_symbol\n"
            "_atom_site_fract_x\n"
            "_atom_site_fract_y\n"
            "_atom_site_fract_z\n"
            "_atom_site_occupancy\n"
            "a1 N 0.100000 0.100000 0.100000 1.00\n"
            "a2 C 0.160000 0.100000 0.100000 1.00\n"
            "a3 O 0.220000 0.100000 0.100000 1.00\n"
            "a4 C 0.160000 0.160000 0.100000 1.00\n"
            "a5 C 0.220000 0.160000 0.100000 1.00\n"
            "\n"
            "loop_\n"
            "_geom_bond_atom_site_label_1\n"
            "_geom_bond_atom_site_label_2\n"
            "_geom_bond_site_symmetry_1\n"
            "_geom_bond_site_symmetry_2\n"
            "_geom_bond_distance\n"
            "_ccdc_geom_bond_type\n"
            "a1 a2 . . 1.000000 S\n"
            "a2 a3 . . 1.000000 D\n"
            "a2 a4 . . 1.000000 S\n"
            "a4 a5 . . 1.000000 S\n"
        )
