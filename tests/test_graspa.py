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
    GraspaWidomSettings,
    resolve_eqeq_binary,
    resolve_graspa_binary,
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
            self.assertEqual(len(result.component_results), 5)
            self.assertTrue(Path(result.eqeq_charged_cif).is_file())
            self.assertTrue(Path(result.results_csv_path).is_file())
            self.assertTrue(Path(result.report_path).is_file())
            self.assertFalse((Path(result.eqeq_run_dir) / "data_C.cif").exists())

            simulation_input = Path(result.simulation_input_path).read_text(encoding="utf-8")
            self.assertIn("UseGPUReduction yes", simulation_input)
            self.assertIn("UseFastHostRNG yes", simulation_input)
            self.assertIn("FrameworkName framework", simulation_input)
            self.assertIn("UnitCells 0 1 2 3", simulation_input)
            self.assertIn("Component 4 MoleculeName              SO2", simulation_input)

            eqeq_audit = json.loads((Path(result.eqeq_run_dir) / "eqeq_invocation.json").read_text(encoding="utf-8"))
            self.assertEqual(eqeq_audit["argv"][0], "example_framework.cif")
            self.assertEqual(eqeq_audit["argv"][4], "ewald")

            graspa_audit = json.loads((Path(result.widom_run_dir) / "graspa_invocation.json").read_text(encoding="utf-8"))
            self.assertEqual(graspa_audit["framework_name"], "framework")
            self.assertEqual(graspa_audit["unit_cells"], [1, 2, 3])
            self.assertEqual(graspa_audit["components"], ["TIP4P", "CO2", "H2", "N2", "SO2"])

            csv_text = Path(result.results_csv_path).read_text(encoding="utf-8")
            self.assertIn("component,widom_energy,widom_energy_errorbar,henry,henry_errorbar", csv_text)
            self.assertIn("CO2", csv_text)

            report = json.loads(Path(result.report_path).read_text(encoding="utf-8"))
            self.assertEqual(report["unit_cells"], [1, 2, 3])
            self.assertTrue(report["widom_settings"]["use_gpu_reduction"])
            self.assertTrue(report["widom_settings"]["use_fast_host_rng"])
            self.assertEqual(report["component_results"][1]["component"], "CO2")
            self.assertEqual(report["component_results"][1]["henry"], 2e-05)

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
                            "--production-cycles",
                            "5000",
                            "--json",
                        ]
                    )

            report = json.loads(buffer.getvalue())
            self.assertEqual(report["output_dir"], str(output_dir.resolve()))
            self.assertEqual(report["unit_cells"], [1, 2, 3])
            self.assertEqual(report["widom_settings"]["production_cycles"], 5000)
            self.assertEqual(len(report["component_results"]), 5)
            self.assertEqual(report["component_results"][0]["component"], "TIP4P")

    def test_calculate_help_lists_graspa_widom(self):
        buffer = io.StringIO()
        with self.assertRaises(SystemExit), contextlib.redirect_stdout(buffer):
            cli_main(["calculate", "--help"])

        self.assertIn("graspa-widom", buffer.getvalue())

    def _write_fake_eqeq_binary(self, path: Path) -> Path:
        path.write_text(
            "#!/usr/bin/env python\n"
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
            "#!/usr/bin/env python\n"
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
            "components = re.findall(r'^Component\\s+\\d+\\s+MoleculeName\\s+(\\S+)', text, re.MULTILINE)\n"
            "if framework_match is None or unit_cells_match is None:\n"
            "    sys.stderr.write('missing framework/unitcells config\\n')\n"
            "    sys.exit(3)\n"
            "\n"
            "framework_name = framework_match.group(1)\n"
            "framework_cif = Path(framework_name + '.cif')\n"
            "if not framework_cif.is_file():\n"
            "    sys.stderr.write('missing framework cif\\n')\n"
            "    sys.exit(4)\n"
            "\n"
            "unit_cells = [int(unit_cells_match.group(1)), int(unit_cells_match.group(2)), int(unit_cells_match.group(3))]\n"
            "Path('graspa_invocation.json').write_text(\n"
            "    json.dumps(\n"
            "        {\n"
            "            'framework_name': framework_name,\n"
            "            'unit_cells': unit_cells,\n"
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
            "sections = []\n"
            "for index, component in enumerate(components):\n"
            "    energy = -10.0 - index\n"
            "    energy_err = 0.1 + index * 0.01\n"
            "    henry = (index + 1) * 1.0e-5\n"
            "    henry_err = (index + 1) * 1.0e-6\n"
            "    sections.append(\n"
            "        '================Rosenbluth Summary For Component '\n"
            "        f'[{index}] ({component})================\\n'\n"
            "        f'Averaged Excess Chemical Potential: {energy:.2f} +/- {energy_err:.2f}\\n'\n"
            "        f'Averaged Henry Coefficient [mol/kg/Pa]: {henry:.8e} +/- {henry_err:.8e}\\n'\n"
            "    )\n"
            "data_path.write_text('\\n'.join(sections), encoding='utf-8')\n"
            "sys.stdout.write('fake graspa stdout\\n')\n",
            encoding="utf-8",
        )
        path.chmod(path.stat().st_mode | stat.S_IEXEC)
        return path


if __name__ == "__main__":
    unittest.main()
