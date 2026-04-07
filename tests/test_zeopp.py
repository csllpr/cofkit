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
from cofkit.zeopp import (
    COFKIT_ZEOPP_ENV_VAR,
    ZeoppConfigurationError,
    analyze_zeopp_pore_properties,
    resolve_zeopp_binary,
)


class ZeoppTests(unittest.TestCase):
    def test_resolve_zeopp_binary_requires_configuration(self):
        with patch.dict(os.environ, {}, clear=True):
            with self.assertRaises(ZeoppConfigurationError) as raised:
                resolve_zeopp_binary()

        self.assertIn(COFKIT_ZEOPP_ENV_VAR, str(raised.exception))

    def test_analyze_zeopp_pore_properties_parses_baseline_and_probe_scans(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            fake_binary = self._write_fake_zeopp_binary(temp_path / "network")
            cif_path = temp_path / "example.cif"
            cif_path.write_text("data_example\n", encoding="utf-8")
            output_dir = temp_path / "zeopp_out"

            with patch.dict(os.environ, {COFKIT_ZEOPP_ENV_VAR: str(fake_binary)}):
                result = analyze_zeopp_pore_properties(
                    cif_path,
                    output_dir=output_dir,
                    probe_radii=(1.2, 1.86),
                    surface_samples_per_atom=300,
                    volume_samples_total=7000,
                )

            basic = result.baseline.basic_pore_properties
            self.assertEqual(basic.largest_included_sphere, 12.34)
            self.assertEqual(basic.largest_free_sphere, 10.11)
            self.assertEqual(basic.axis_aligned_free_sphere["a"], 7.1)
            self.assertEqual(basic.axis_aligned_included_sphere_along_free_path["c"], 8.3)

            point_probe_channels = result.baseline.point_probe_channels
            self.assertEqual(point_probe_channels.n_channels, 1)
            self.assertEqual(point_probe_channels.n_pockets, 0)
            self.assertEqual(point_probe_channels.probe_radius, 0.0)

            point_probe_surface = result.baseline.point_probe_surface_area
            self.assertEqual(point_probe_surface.accessible_surface_area_a2, 1000.0)
            self.assertEqual(point_probe_surface.n_channels, 1)

            point_probe_volume = result.baseline.point_probe_volume
            self.assertEqual(point_probe_volume.accessible_volume_a3, 5000.0)
            self.assertEqual(point_probe_volume.n_channels, 1)

            self.assertEqual(len(result.probe_scans), 2)
            first_scan = result.probe_scans[0]
            self.assertEqual(first_scan.status, "ok")
            self.assertEqual(first_scan.settings.probe_radius, 1.2)
            self.assertEqual(first_scan.channel_summary.n_channels, 2)
            self.assertEqual(first_scan.channel_summary.n_pockets, 1)
            self.assertEqual(first_scan.accessibility.n_voronoi_nodes, 10)
            self.assertEqual(first_scan.accessibility.n_accessible_nodes, 6)
            self.assertAlmostEqual(first_scan.surface_area.accessible_surface_area_a2, 880.0)
            self.assertAlmostEqual(first_scan.volume.accessible_volume_a3, 4400.0)

            second_scan = result.probe_scans[1]
            self.assertEqual(second_scan.status, "ok")
            self.assertAlmostEqual(second_scan.surface_area.accessible_surface_area_a2, 814.0)
            self.assertAlmostEqual(second_scan.volume.accessible_volume_a3, 4070.0)

            self.assertTrue(Path(result.report_path).is_file())
            report = json.loads(Path(result.report_path).read_text(encoding="utf-8"))
            self.assertEqual(report["baseline"]["basic_pore_properties"]["largest_free_sphere"], 10.11)
            self.assertEqual(report["probe_scans"][0]["accessibility"]["n_accessible_nodes"], 6)

    def test_probe_scan_failures_are_recorded_without_aborting_default_run(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            fake_binary = self._write_fake_zeopp_binary(temp_path / "network")
            cif_path = temp_path / "example.cif"
            cif_path.write_text("data_example\n", encoding="utf-8")

            with patch.dict(os.environ, {COFKIT_ZEOPP_ENV_VAR: str(fake_binary)}):
                result = analyze_zeopp_pore_properties(
                    cif_path,
                    output_dir=temp_path / "zeopp_out",
                    probe_radii=(2.5,),
                )

            self.assertEqual(result.probe_scans[0].status, "error")
            self.assertIn("sa:", result.probe_scans[0].error)
            self.assertIn("vol:", result.probe_scans[0].error)
            self.assertIsNotNone(result.probe_scans[0].channel_summary)
            self.assertIsNotNone(result.probe_scans[0].accessibility)
            self.assertIsNone(result.probe_scans[0].surface_area)
            self.assertIsNone(result.probe_scans[0].volume)

    def test_analyze_zeopp_cli_prints_json_report(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            fake_binary = self._write_fake_zeopp_binary(temp_path / "network")
            cif_path = temp_path / "example.cif"
            cif_path.write_text("data_example\n", encoding="utf-8")
            output_dir = temp_path / "cli_zeopp"
            buffer = io.StringIO()

            with patch.dict(os.environ, {COFKIT_ZEOPP_ENV_VAR: str(fake_binary)}):
                with contextlib.redirect_stdout(buffer):
                    cli_main(
                        [
                            "analyze",
                            "zeopp",
                            str(cif_path),
                            "--output-dir",
                            str(output_dir),
                            "--probe-radius",
                            "1.2",
                            "--probe-radius",
                            "1.86",
                            "--json",
                        ]
                    )

            report = json.loads(buffer.getvalue())
            self.assertEqual(report["baseline"]["point_probe_surface_area"]["accessible_surface_area_a2"], 1000.0)
            self.assertEqual(len(report["probe_scans"]), 2)
            self.assertEqual(report["probe_scans"][1]["settings"]["probe_radius"], 1.86)
            self.assertEqual(report["probe_scans"][1]["channel_summary"]["n_channels"], 2)
            self.assertEqual(report["output_dir"], str(output_dir.resolve()))

    def test_analyze_help_lists_zeopp(self):
        buffer = io.StringIO()
        with self.assertRaises(SystemExit), contextlib.redirect_stdout(buffer):
            cli_main(["analyze", "--help"])

        self.assertIn("zeopp", buffer.getvalue())

    def _write_fake_zeopp_binary(self, path: Path) -> Path:
        path.write_text(
            f"#!{sys.executable}\n"
            "from __future__ import annotations\n"
            "import sys\n"
            "from pathlib import Path\n"
            "\n"
            "args = sys.argv[1:]\n"
            "if not args:\n"
            "    sys.exit(2)\n"
            "\n"
            "command = None\n"
            "payload = []\n"
            "index = 0\n"
            "while index < len(args):\n"
            "    value = args[index]\n"
            "    if value.startswith('-') and command is None:\n"
            "        command = value\n"
            "        if value in {'-res', '-resex'}:\n"
            "            payload = args[index + 1 : index + 2]\n"
            "            index += 2\n"
            "        elif value == '-chan':\n"
            "            payload = args[index + 1 : index + 3]\n"
            "            index += 3\n"
            "        elif value in {'-sa', '-vol'}:\n"
            "            payload = args[index + 1 : index + 5]\n"
            "            index += 5\n"
            "        elif value == '-axs':\n"
            "            payload = args[index + 1 : index + 3]\n"
            "            index += 3\n"
            "        else:\n"
            "            sys.stderr.write('unsupported arguments\\n')\n"
            "            sys.exit(2)\n"
            "    else:\n"
            "        break\n"
            "\n"
            "input_path = Path(args[-1])\n"
            "if not input_path.is_file():\n"
            "    sys.stderr.write('missing input\\n')\n"
            "    sys.exit(4)\n"
            "\n"
            "def write_text(path_str: str, text: str) -> None:\n"
            "    path = Path(path_str)\n"
            "    path.parent.mkdir(parents=True, exist_ok=True)\n"
            "    path.write_text(text, encoding='utf-8')\n"
            "\n"
            "if command == '-res':\n"
            "    output_path = payload[0]\n"
            "    write_text(output_path, f'{output_path}    12.34 10.11 9.87\\n')\n"
            "elif command == '-resex':\n"
            "    output_path = payload[0]\n"
            "    write_text(output_path, f'{output_path}    12.34 10.11 9.87  7.10  7.20  7.30  8.10  8.20  8.30\\n')\n"
            "elif command == '-chan':\n"
            "    probe = float(payload[0])\n"
            "    output_path = payload[1]\n"
            "    if probe == 0:\n"
            "        sys.stdout.write('Identified 1 channels and 0 pockets.\\n')\n"
            "        write_text(\n"
            "            output_path,\n"
            "            f'{output_path}   1 channels identified of dimensionality 3\\n'\n"
            "            'Channel  0  12.34  10.11  9.87\\n'\n"
            "            f'{output_path} summary(Max_of_columns_above)   12.34 10.11  9.87  probe_rad: 0  probe_diam: 0\\n',\n"
            "        )\n"
            "    else:\n"
            "        sys.stdout.write('Identified 2 channels and 1 pockets.\\n')\n"
            "        lis = round(9.0 + probe, 2)\n"
            "        lfs = round(8.0 + probe, 2)\n"
            "        lif = round(7.0 + probe, 2)\n"
            "        write_text(\n"
            "            output_path,\n"
            "            f'{output_path}   2 channels identified of dimensionality 2 3\\n'\n"
            "            f'Channel  0  {lis:.2f}  {lfs:.2f}  {lif:.2f}\\n'\n"
            "            f'Channel  1  {lis - 1:.2f}  {lfs - 1:.2f}  {lif - 1:.2f}\\n'\n"
            "            f'{output_path} summary(Max_of_columns_above)   {lis:.2f} {lfs:.2f}  {lif:.2f}  probe_rad: {probe:g}  probe_diam: {2 * probe:g}\\n',\n"
            "        )\n"
            "elif command == '-sa':\n"
            "    channel_radius = float(payload[0])\n"
            "    probe = float(payload[1])\n"
            "    output_path = payload[3]\n"
            "    if probe > 2.0:\n"
            "        sys.stderr.write('surface area failure\\n')\n"
            "        sys.exit(134)\n"
            "    accessible = round(1000.0 - probe * 100.0, 2)\n"
            "    inaccessible = round(probe * 25.0, 2)\n"
            "    n_channels = 1 if probe == 0 else 2\n"
            "    n_pockets = 0 if probe == 0 else 1\n"
            "    channel_values = ' '.join(str(round(accessible / max(n_channels, 1), 2)) for _ in range(n_channels))\n"
            "    pocket_values = '' if n_pockets == 0 else str(round(inaccessible, 2))\n"
            "    write_text(\n"
            "        output_path,\n"
            "        f'@ {output_path} Unitcell_volume: 321.0   Density: 0.12345   ASA_A^2: {accessible:.2f} ASA_m^2/cm^3: {round(accessible / 2.0, 2):.2f} ASA_m^2/g: {round(accessible / 3.0, 2):.2f} NASA_A^2: {inaccessible:.2f} NASA_m^2/cm^3: {round(inaccessible / 2.0, 2):.2f} NASA_m^2/g: {round(inaccessible / 3.0, 2):.2f}\\n'\n"
            "        f'Number_of_channels: {n_channels} Channel_surface_area_A^2: {channel_values}\\n'\n"
            "        f'Number_of_pockets: {n_pockets} Pocket_surface_area_A^2: {pocket_values}\\n',\n"
            "    )\n"
            "elif command == '-vol':\n"
            "    channel_radius = float(payload[0])\n"
            "    probe = float(payload[1])\n"
            "    output_path = payload[3]\n"
            "    if probe > 2.0:\n"
            "        sys.stderr.write('volume failure\\n')\n"
            "        sys.exit(134)\n"
            "    accessible = round(5000.0 - probe * 500.0, 2)\n"
            "    inaccessible = round(probe * 50.0, 2)\n"
            "    n_channels = 1 if probe == 0 else 2\n"
            "    n_pockets = 0 if probe == 0 else 1\n"
            "    channel_values = ' '.join(str(round(accessible / max(n_channels, 1), 2)) for _ in range(n_channels))\n"
            "    pocket_values = '' if n_pockets == 0 else str(round(inaccessible, 2))\n"
            "    write_text(\n"
            "        output_path,\n"
            "        f'@ {output_path} Unitcell_volume: 321.0   Density: 0.12345   AV_A^3: {accessible:.2f} AV_Volume_fraction: {round(accessible / 10000.0, 5):.5f} AV_cm^3/g: {round(accessible / 345.0, 5):.5f} NAV_A^3: {inaccessible:.2f} NAV_Volume_fraction: {round(inaccessible / 10000.0, 5):.5f} NAV_cm^3/g: {round(inaccessible / 345.0, 5):.5f}\\n'\n"
            "        f'Number_of_channels: {n_channels} Channel_volume_A^3: {channel_values}\\n'\n"
            "        f'Number_of_pockets: {n_pockets} Pocket_volume_A^3: {pocket_values}\\n',\n"
            "    )\n"
            "elif command == '-axs':\n"
            "    probe = float(payload[0])\n"
            "    output_path = payload[1]\n"
            "    n_true = 8 if probe == 0 else 6\n"
            "    n_false = 2 if probe == 0 else 4\n"
            "    write_text(output_path, ''.join(['true\\n'] * n_true + ['false\\n'] * n_false))\n"
            "else:\n"
            "    sys.stderr.write('unsupported arguments\\n')\n"
            "    sys.exit(2)\n",
            encoding="utf-8",
        )
        path.chmod(path.stat().st_mode | stat.S_IXUSR)
        return path


if __name__ == "__main__":
    unittest.main()
