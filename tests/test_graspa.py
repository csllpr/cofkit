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
from cofkit.graspa import (
    COFKIT_EQEQ_ENV_VAR,
    COFKIT_GRASPA_ENV_VAR,
    EqeqChargeSettings,
    GraspaConfigurationError,
    GraspaIsothermSettings,
    GraspaWidomSettings,
    resolve_eqeq_binary,
    resolve_graspa_binary,
    run_graspa_isotherm_workflow,
    run_graspa_widom_workflow,
)


class GraspaWidomTests(unittest.TestCase):
    def test_resolve_eqeq_binary_requires_configuration(self):
        with patch.dict(os.environ, {}, clear=True):
            with self.assertRaises(GraspaConfigurationError) as raised:
                resolve_eqeq_binary()

        self.assertIn(COFKIT_EQEQ_ENV_VAR, str(raised.exception))

    def test_resolve_graspa_binary_requires_configuration(self):
        with patch.dict(os.environ, {}, clear=True):
            with self.assertRaises(GraspaConfigurationError) as raised:
                resolve_graspa_binary()

        self.assertIn(COFKIT_GRASPA_ENV_VAR, str(raised.exception))

    def test_run_graspa_widom_workflow_prepares_inputs_and_parses_results(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            fake_eqeq = self._write_fake_eqeq_binary(temp_path / "eqeq_fake")
            fake_graspa = self._write_fake_graspa_binary(temp_path / "graspa_fake")
            cif_path = temp_path / "example_framework.cif"
            cif_path.write_text(
                "data_example\n"
                "_cell_length_a 26.0\n"
                "_cell_length_b 13.0\n"
                "_cell_length_c 9.0\n",
                encoding="utf-8",
            )

            with patch.dict(
                os.environ,
                {
                    COFKIT_EQEQ_ENV_VAR: str(fake_eqeq),
                    COFKIT_GRASPA_ENV_VAR: str(fake_graspa),
                },
                clear=False,
            ):
                result = run_graspa_widom_workflow(
                    cif_path,
                    output_dir=temp_path / "widom_out",
                    eqeq_settings=EqeqChargeSettings(),
                    widom_settings=GraspaWidomSettings(),
                    graspa_timeout_seconds=30.0,
                )

            self.assertEqual(result.unit_cells, (1, 2, 3))
            self.assertEqual(len(result.component_results), 7)
            self.assertTrue(Path(result.eqeq_charged_cif).is_file())
            self.assertTrue(Path(result.results_csv_path).is_file())
            self.assertTrue(Path(result.report_path).is_file())
            self.assertFalse((Path(result.eqeq_run_dir) / "data_C.cif").exists())
            self.assertTrue((Path(result.widom_run_dir) / "Xe.def").is_file())
            self.assertTrue((Path(result.widom_run_dir) / "Kr.def").is_file())

            simulation_input = Path(result.simulation_input_path).read_text(encoding="utf-8")
            self.assertIn("UseGPUReduction yes", simulation_input)
            self.assertIn("UseFastHostRNG yes", simulation_input)
            self.assertIn("FrameworkName framework", simulation_input)
            self.assertIn("UnitCells 0 1 2 3", simulation_input)
            self.assertIn("NumberOfBlocks 5", simulation_input)
            self.assertIn("NumberOfProductionCycles     2000000", simulation_input)
            self.assertIn("Component 4 MoleculeName              SO2", simulation_input)
            self.assertIn("Component 5 MoleculeName              Xe", simulation_input)
            self.assertIn("Component 6 MoleculeName              Kr", simulation_input)

            eqeq_audit = json.loads((Path(result.eqeq_run_dir) / "eqeq_invocation.json").read_text(encoding="utf-8"))
            self.assertEqual(eqeq_audit["argv"][0], "example_framework.cif")
            self.assertEqual(eqeq_audit["argv"][4], "ewald")

            graspa_audit = json.loads((Path(result.widom_run_dir) / "graspa_invocation.json").read_text(encoding="utf-8"))
            self.assertEqual(graspa_audit["framework_name"], "framework")
            self.assertEqual(graspa_audit["unit_cells"], [1, 2, 3])
            self.assertEqual(graspa_audit["components"], ["TIP4P", "CO2", "H2", "N2", "SO2", "Xe", "Kr"])

            csv_text = Path(result.results_csv_path).read_text(encoding="utf-8")
            self.assertIn("component,widom_energy,widom_energy_errorbar,henry,henry_errorbar", csv_text)
            self.assertIn("CO2", csv_text)
            self.assertIn("Xe", csv_text)
            self.assertIn("Kr", csv_text)

            report = json.loads(Path(result.report_path).read_text(encoding="utf-8"))
            self.assertEqual(report["unit_cells"], [1, 2, 3])
            self.assertTrue(report["widom_settings"]["use_gpu_reduction"])
            self.assertTrue(report["widom_settings"]["use_fast_host_rng"])
            self.assertEqual(report["widom_settings"]["number_of_blocks"], 5)
            self.assertEqual(report["widom_settings"]["production_cycles"], 2000000)
            self.assertIsNone(report["widom_settings"]["widom_moves_per_component"])
            self.assertEqual(report["component_results"][1]["component"], "CO2")
            self.assertEqual(report["component_results"][1]["henry"], 2e-05)
            self.assertEqual(report["component_results"][5]["component"], "Xe")
            self.assertEqual(report["component_results"][5]["henry"], 6e-05)
            self.assertEqual(report["component_results"][6]["component"], "Kr")
            self.assertEqual(report["component_results"][6]["henry"], 7e-05)

    def test_calculate_graspa_widom_cli_prints_json_report(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            fake_eqeq = self._write_fake_eqeq_binary(temp_path / "eqeq_fake")
            fake_graspa = self._write_fake_graspa_binary(temp_path / "graspa_fake")
            cif_path = temp_path / "cli_example.cif"
            cif_path.write_text(
                "data_example\n"
                "_cell_length_a 26.0\n"
                "_cell_length_b 13.0\n"
                "_cell_length_c 9.0\n",
                encoding="utf-8",
            )
            output_dir = temp_path / "cli_widom"
            buffer = io.StringIO()

            with patch.dict(
                os.environ,
                {
                    COFKIT_EQEQ_ENV_VAR: str(fake_eqeq),
                    COFKIT_GRASPA_ENV_VAR: str(fake_graspa),
                },
                clear=False,
            ):
                with contextlib.redirect_stdout(buffer):
                    cli_main(
                        [
                            "calculate",
                            "graspa-widom",
                            str(cif_path),
                            "--output-dir",
                            str(output_dir),
                            "--component",
                            "CO2",
                            "--component",
                            "Xe",
                            "--widom-moves-per-component",
                            "2500",
                            "--json",
                        ]
                    )

            report = json.loads(buffer.getvalue())
            self.assertEqual(report["output_dir"], str(output_dir.resolve()))
            self.assertEqual(report["unit_cells"], [1, 2, 3])
            self.assertEqual(report["widom_settings"]["widom_moves_per_component"], 2500)
            self.assertEqual(report["widom_settings"]["production_cycles"], 5000)
            self.assertEqual(len(report["component_results"]), 2)
            self.assertEqual(report["component_results"][0]["component"], "CO2")
            self.assertEqual(report["component_results"][1]["component"], "Xe")

    def test_run_graspa_isotherm_workflow_prepares_inputs_and_parses_results(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            fake_eqeq = self._write_fake_eqeq_binary(temp_path / "eqeq_fake")
            fake_graspa = self._write_fake_graspa_binary(temp_path / "graspa_fake")
            cif_path = temp_path / "adsorption_framework.cif"
            cif_path.write_text(
                "data_example\n"
                "_cell_length_a 26.0\n"
                "_cell_length_b 13.0\n"
                "_cell_length_c 9.0\n",
                encoding="utf-8",
            )

            with patch.dict(
                os.environ,
                {
                    COFKIT_EQEQ_ENV_VAR: str(fake_eqeq),
                    COFKIT_GRASPA_ENV_VAR: str(fake_graspa),
                },
                clear=False,
            ):
                result = run_graspa_isotherm_workflow(
                    cif_path,
                    output_dir=temp_path / "isotherm_out",
                    eqeq_settings=EqeqChargeSettings(),
                    isotherm_settings=GraspaIsothermSettings(component="CO2", pressures=(10000.0, 500000.0)),
                    graspa_timeout_seconds=30.0,
                )

            self.assertEqual(result.unit_cells, (1, 2, 3))
            self.assertEqual(len(result.point_results), 2)
            self.assertTrue(Path(result.eqeq_charged_cif).is_file())
            self.assertTrue(Path(result.results_csv_path).is_file())
            self.assertTrue(Path(result.report_path).is_file())
            self.assertTrue(Path(result.isotherm_framework_cif).is_file())

            first_point = result.point_results[0]
            second_point = result.point_results[1]
            self.assertEqual(first_point.component, "CO2")
            self.assertEqual(first_point.pressure, 10000.0)
            self.assertAlmostEqual(first_point.loading_mol_per_kg, 0.025)
            self.assertAlmostEqual(first_point.loading_g_per_l, 0.25)
            self.assertAlmostEqual(first_point.heat_of_adsorption_kj_per_mol, -20.1)
            self.assertEqual(second_point.pressure, 500000.0)
            self.assertAlmostEqual(second_point.loading_mol_per_kg, 1.25)
            self.assertAlmostEqual(second_point.loading_g_per_l, 12.5)
            self.assertAlmostEqual(second_point.heat_of_adsorption_kj_per_mol, -25.0)

            simulation_input = Path(first_point.simulation_input_path).read_text(encoding="utf-8")
            self.assertIn("NumberOfBlocks 5", simulation_input)
            self.assertIn("Component 0 MoleculeName             CO2", simulation_input)
            self.assertIn("Pressure     10000", simulation_input)
            self.assertIn("TranslationProbability   1", simulation_input)
            self.assertIn("RotationProbability      1", simulation_input)
            self.assertIn("ReinsertionProbability   1", simulation_input)
            self.assertIn("SwapProbability          1", simulation_input)
            self.assertNotIn("WidomProbability", simulation_input)

            csv_text = Path(result.results_csv_path).read_text(encoding="utf-8")
            self.assertIn("pressure_pa", csv_text)
            self.assertIn("loading_mol_per_kg", csv_text)
            self.assertIn("heat_of_adsorption_kj_per_mol", csv_text)

            report = json.loads(Path(result.report_path).read_text(encoding="utf-8"))
            self.assertEqual(report["unit_cells"], [1, 2, 3])
            self.assertEqual(report["isotherm_settings"]["component"], "CO2")
            self.assertEqual(report["isotherm_settings"]["pressures"], [10000.0, 500000.0])
            self.assertEqual(report["point_results"][0]["loading_mol_per_kg"], 0.025)
            self.assertEqual(report["point_results"][0]["loading_mmol_per_g"], 0.025)
            self.assertEqual(report["point_results"][1]["loading_mol_per_kg"], 1.25)

    def test_calculate_graspa_isotherm_cli_prints_json_report(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            fake_eqeq = self._write_fake_eqeq_binary(temp_path / "eqeq_fake")
            fake_graspa = self._write_fake_graspa_binary(temp_path / "graspa_fake")
            cif_path = temp_path / "cli_adsorption_example.cif"
            cif_path.write_text(
                "data_example\n"
                "_cell_length_a 26.0\n"
                "_cell_length_b 13.0\n"
                "_cell_length_c 9.0\n",
                encoding="utf-8",
            )
            output_dir = temp_path / "cli_isotherm"
            buffer = io.StringIO()

            with patch.dict(
                os.environ,
                {
                    COFKIT_EQEQ_ENV_VAR: str(fake_eqeq),
                    COFKIT_GRASPA_ENV_VAR: str(fake_graspa),
                },
                clear=False,
            ):
                with contextlib.redirect_stdout(buffer):
                    cli_main(
                        [
                            "calculate",
                            "graspa-isotherm",
                            str(cif_path),
                            "--output-dir",
                            str(output_dir),
                            "--component",
                            "Xe",
                            "--pressure",
                            "10000",
                            "--pressure",
                            "100000",
                            "--fugacity-coefficient",
                            "PR-EOS",
                            "--json",
                        ]
                    )

            report = json.loads(buffer.getvalue())
            self.assertEqual(report["output_dir"], str(output_dir.resolve()))
            self.assertEqual(report["unit_cells"], [1, 2, 3])
            self.assertEqual(report["isotherm_settings"]["component"], "Xe")
            self.assertEqual(report["isotherm_settings"]["fugacity_coefficient"], "PR-EOS")
            self.assertEqual(len(report["point_results"]), 2)
            self.assertEqual(report["point_results"][0]["pressure"], 10000.0)
            self.assertEqual(report["point_results"][1]["pressure"], 100000.0)

            simulation_input = Path(report["point_results"][0]["simulation_input_path"]).read_text(encoding="utf-8")
            self.assertIn("FugacityCoefficient      PR-EOS", simulation_input)
            self.assertNotIn("RotationProbability", simulation_input)

    def test_calculate_help_lists_graspa_widom(self):
        buffer = io.StringIO()
        with self.assertRaises(SystemExit), contextlib.redirect_stdout(buffer):
            cli_main(["calculate", "--help"])

        self.assertIn("graspa-widom", buffer.getvalue())
        self.assertIn("graspa-isotherm", buffer.getvalue())

    def _write_fake_eqeq_binary(self, path: Path) -> Path:
        path.write_text(
            f"#!{sys.executable}\n"
            "from __future__ import annotations\n"
            "import json\n"
            "import sys\n"
            "from pathlib import Path\n"
            "\n"
            "args = sys.argv[1:]\n"
            "if len(args) < 8:\n"
            "    sys.stderr.write('insufficient arguments\\n')\n"
            "    sys.exit(2)\n"
            "\n"
            "input_name = args[0]\n"
            "lambda_value = float(args[1])\n"
            "h_i0 = float(args[2])\n"
            "method = args[4]\n"
            "input_path = Path(input_name)\n"
            "if not input_path.is_file():\n"
            "    sys.stderr.write('missing input\\n')\n"
            "    sys.exit(3)\n"
            "\n"
            "stem = f'{input_name}_EQeq_{method}_{lambda_value:.2f}_{h_i0:.2f}'\n"
            "input_text = input_path.read_text(encoding='utf-8')\n"
            "Path(stem + '.cif').write_text(input_text + '\\n# fake charged output\\n', encoding='utf-8')\n"
            "Path(stem + '.json').write_text(json.dumps({'argv': args}, indent=2), encoding='utf-8')\n"
            "Path('eqeq_invocation.json').write_text(json.dumps({'argv': args}, indent=2), encoding='utf-8')\n"
            "sys.stdout.write('fake eqeq stdout\\n')\n"
            "sys.stderr.write('fake eqeq stderr\\n')\n",
            encoding="utf-8",
        )
        path.chmod(path.stat().st_mode | stat.S_IEXEC)
        return path

    def _write_fake_graspa_binary(self, path: Path) -> Path:
        path.write_text(
            f"#!{sys.executable}\n"
            "from __future__ import annotations\n"
            "import json\n"
            "import re\n"
            "import sys\n"
            "from pathlib import Path\n"
            "\n"
            "simulation_input = Path('simulation.input')\n"
            "if not simulation_input.is_file():\n"
            "    sys.stderr.write('missing simulation.input\\n')\n"
            "    sys.exit(2)\n"
            "\n"
            "text = simulation_input.read_text(encoding='utf-8')\n"
            "framework_match = re.search(r'^FrameworkName\\s+(\\S+)', text, re.MULTILINE)\n"
            "unit_cells_match = re.search(r'^UnitCells\\s+0\\s+(\\d+)\\s+(\\d+)\\s+(\\d+)', text, re.MULTILINE)\n"
            "pressure_match = re.search(r'^Pressure\\s+(\\S+)', text, re.MULTILINE)\n"
            "components = re.findall(r'^Component\\s+\\d+\\s+MoleculeName\\s+(\\S+)', text, re.MULTILINE)\n"
            "if framework_match is None or unit_cells_match is None or pressure_match is None:\n"
            "    sys.stderr.write('missing framework/unitcells config\\n')\n"
            "    sys.exit(3)\n"
            "\n"
            "framework_name = framework_match.group(1)\n"
            "pressure = float(pressure_match.group(1))\n"
            "mode = 'widom' if 'WidomProbability' in text else 'isotherm'\n"
            "framework_cif = Path(framework_name + '.cif')\n"
            "if not framework_cif.is_file():\n"
            "    sys.stderr.write('missing framework cif\\n')\n"
            "    sys.exit(4)\n"
            "\n"
            "unit_cells = [int(unit_cells_match.group(1)), int(unit_cells_match.group(2)), int(unit_cells_match.group(3))]\n"
            "Path('graspa_invocation.json').write_text(\n"
            "    json.dumps(\n"
            "        {\n"
            "            'mode': mode,\n"
            "            'framework_name': framework_name,\n"
            "            'unit_cells': unit_cells,\n"
            "            'pressure': pressure,\n"
            "            'components': components,\n"
            "        },\n"
            "        indent=2,\n"
            "    ),\n"
            "    encoding='utf-8',\n"
            ")\n"
            "\n"
            "output_dir = Path('Output')\n"
            "output_dir.mkdir(exist_ok=True)\n"
            "data_path = output_dir / 'System_0_fake.data'\n"
            "if mode == 'widom':\n"
            "    sections = []\n"
            "    for index, component in enumerate(components):\n"
            "        energy = -10.0 - index\n"
            "        energy_err = 0.1 + index * 0.01\n"
            "        henry = (index + 1) * 1.0e-5\n"
            "        henry_err = (index + 1) * 1.0e-6\n"
            "        sections.append(\n"
            "            '================Rosenbluth Summary For Component '\n"
            "            f'[{index}] ({component})================\\n'\n"
            "            f'Averaged Excess Chemical Potential: {energy:.2f} +/- {energy_err:.2f}\\n'\n"
            "            f'Averaged Henry Coefficient [mol/kg/Pa]: {henry:.8e} +/- {henry_err:.8e}\\n'\n"
            "        )\n"
            "    data_path.write_text('\\n'.join(sections), encoding='utf-8')\n"
            "else:\n"
            "    component = components[0]\n"
            "    loading_mol = pressure / 100000.0 * 0.25\n"
            "    loading_err = loading_mol * 0.1\n"
            "    loading_gpl = loading_mol * 10.0\n"
            "    loading_gpl_err = loading_err * 10.0\n"
            "    heat = -20.0 - pressure / 100000.0\n"
            "    heat_err = 0.5\n"
            "    data_path.write_text(\n"
            "        '==============================================================\\n'\n"
            "        '=====================BLOCK AVERAGES (LOADING: mol/kg)=============\\n'\n"
            "        f'COMPONENT [1] ({component})\\n'\n"
            "        f'Overall: Average: {loading_mol:.5f}, ErrorBar: {loading_err:.5f}\\n'\n"
            "        '----------------------------------------------------------\\n'\n"
            "        '==============================================================\\n'\n"
            "        '=====================BLOCK AVERAGES (LOADING: g/L)=============\\n'\n"
            "        f'COMPONENT [1] ({component})\\n'\n"
            "        f'Overall: Average: {loading_gpl:.5f}, ErrorBar: {loading_gpl_err:.5f}\\n'\n"
            "        '----------------------------------------------------------\\n'\n"
            "        '==============================================================\\n'\n"
            "        '============= BLOCK AVERAGES (HEAT OF ADSORPTION: kJ/mol) =========\\n'\n"
            "        f'COMPONENT [1] ({component})\\n'\n"
            "        f'Overall: Average: {heat:.5f}, ErrorBar: {heat_err:.5f}\\n'\n"
            "        '-----------------------------\\n',\n"
            "        encoding='utf-8',\n"
            "    )\n"
            "sys.stdout.write('fake graspa stdout\\n')\n",
            encoding="utf-8",
        )
        path.chmod(path.stat().st_mode | stat.S_IEXEC)
        return path


if __name__ == "__main__":
    unittest.main()
