import contextlib
import io
import json
import math
import os
import re
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
    COFKIT_RASPA2_ENV_VAR,
    EqeqChargeSettings,
    GraspaConfigurationError,
    GraspaMixtureComponentSettings,
    GraspaMixtureSettings,
    GraspaIsothermSettings,
    GraspaWidomSettings,
    PACKAGED_GUEST_FORCEFIELD_METADATA,
    resolve_eqeq_binary,
    resolve_graspa_binary,
    resolve_raspa2_binary,
    run_graspa_mixture_workflow,
    run_graspa_isotherm_workflow,
    run_graspa_widom_workflow,
)
from cofkit.guest_bundles import GuestBundleError, load_guest_bundle


class GraspaWidomTests(unittest.TestCase):
    def test_load_guest_bundle_rejects_legacy_schema_without_forcefield_metadata(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            bundle_path = Path(temp_dir) / "legacy_guest.json"
            bundle_path.write_text(json.dumps({"version": 1, "name": "CH4"}), encoding="utf-8")

            with self.assertRaises(GuestBundleError) as raised:
                load_guest_bundle(bundle_path)

        self.assertIn("expected 2", str(raised.exception))
        self.assertIn("force-field provenance", str(raised.exception))

    def test_load_guest_bundle_requires_forcefield_metadata(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            bundle_path = Path(temp_dir) / "missing_forcefield_metadata.json"
            bundle_path.write_text(json.dumps({"version": 2, "name": "CH4"}), encoding="utf-8")

            with self.assertRaises(GuestBundleError) as raised:
                load_guest_bundle(bundle_path)

        self.assertIn("parameter_family", str(raised.exception))

    def test_load_guest_bundle_requires_synchronized_lammps_section(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            bundle_path = temp_path / "bad_guest.json"
            bundle_path.write_text(
                json.dumps(
                    {
                        "version": 2,
                        "name": "CH4",
                        "parameter_family": "dreiding",
                        "parameter_source": "test methane model",
                        "compatible_framework_forcefields": ["dreiding"],
                        "raspa": {
                            "molecule_definition": "# methane\n",
                            "pseudo_atom_rows": [
                                "C_ch4 yes C C 0 12.01070 0.0 0.0 1.0 0.720 0 0 relative 0",
                            ],
                            "mixing_rule_rows": [
                                "C_ch4 lennard-jones 148.0 3.73",
                            ],
                        },
                    }
                ),
                encoding="utf-8",
            )

            with self.assertRaises(GuestBundleError) as raised:
                load_guest_bundle(bundle_path)

            self.assertIn("lammps", str(raised.exception))

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

    def test_resolve_raspa2_binary_requires_configuration(self):
        with patch.dict(os.environ, {}, clear=True):
            with self.assertRaises(GraspaConfigurationError) as raised:
                resolve_raspa2_binary()

        self.assertIn(COFKIT_RASPA2_ENV_VAR, str(raised.exception))

    def test_run_graspa_widom_workflow_prepares_inputs_and_parses_results(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            fake_eqeq = self._write_fake_eqeq_binary(temp_path / "eqeq_fake", strip_leading_cofid_comment=True)
            fake_graspa = self._write_fake_graspa_binary(temp_path / "graspa_fake")
            cif_path = temp_path / "example_framework.cif"
            cofid = "3:amine:Nc1ccc(-c2cc(-c3ccc(N)cc3)cc(-c3ccc(N)cc3)c2)cc1.2:aldehyde:O=Cc1ccc(C=O)cc1&&hcb&&imine"
            cif_path.write_text(
                f"# COFid: {cofid}\n"
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
            self.assertEqual(Path(result.eqeq_charged_cif).read_text(encoding="utf-8").splitlines()[0], f"# COFid: {cofid}")
            self.assertEqual(
                Path(result.widom_framework_cif).read_text(encoding="utf-8").splitlines()[0],
                f"# COFid: {cofid}",
            )
            mixing_rules_text = (Path(result.widom_run_dir) / "force_field_mixing_rules.def").read_text(encoding="utf-8")
            self.assertIn("// standard DREIDING Tables I-II (Mayo et al., J. Phys. Chem. 1990)", mixing_rules_text)
            self.assertRegex(
                mixing_rules_text,
                re.compile(r"^Br\s+lennard-jones\s+186\.1924\s+3\.51905\b", re.MULTILINE),
            )
            self.assertNotIn("126.2900", mixing_rules_text)
            guest_metadata = json.loads(
                (Path(result.widom_run_dir) / "guest_forcefield_metadata.json").read_text(encoding="utf-8")
            )
            self.assertEqual(guest_metadata["framework_forcefield"], "dreiding")
            self.assertEqual(
                guest_metadata["framework_forcefield_id"],
                "dreiding-standard-1990-cofkit-1.0",
            )
            framework_metadata_path = Path(result.widom_run_dir) / "framework_forcefield_metadata.json"
            framework_metadata = json.loads(framework_metadata_path.read_text(encoding="utf-8"))
            self.assertEqual(framework_metadata["schema_version"], 1)
            self.assertEqual(framework_metadata["id"], "dreiding-standard-1990-cofkit-1.0")
            self.assertEqual(result.framework_forcefield_metadata, framework_metadata)
            self.assertEqual(
                [(guest["name"], guest["parameter_family"]) for guest in guest_metadata["guests"]],
                [
                    ("TIP4P", "dreiding"),
                    ("CO2", "dreiding"),
                    ("H2", "dreiding"),
                    ("N2", "dreiding"),
                    ("SO2", "dreiding"),
                    ("Xe", "uff"),
                    ("Kr", "uff"),
                ],
            )

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
            self.assertEqual(
                report["guest_forcefield_metadata_paths"],
                [str(Path(result.widom_run_dir) / "guest_forcefield_metadata.json")],
            )
            self.assertEqual(
                report["framework_forcefield_metadata_paths"],
                [str(framework_metadata_path)],
            )
            self.assertEqual(report["framework_forcefield_metadata"], framework_metadata)
            self.assertEqual(report["component_results"][1]["component"], "CO2")
            self.assertEqual(report["component_results"][1]["henry"], 2e-05)
            self.assertEqual(report["component_results"][5]["component"], "Xe")
            self.assertEqual(report["component_results"][5]["henry"], 6e-05)
            self.assertEqual(report["component_results"][6]["component"], "Kr")
            self.assertEqual(report["component_results"][6]["henry"], 7e-05)

    def test_run_graspa_isotherm_workflow_accepts_external_guest_bundle(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            fake_eqeq = self._write_fake_eqeq_binary(temp_path / "eqeq_fake")
            fake_graspa = self._write_fake_graspa_binary(temp_path / "graspa_fake")
            bundle_path = self._write_guest_bundle(temp_path / "ch4_bundle.json")
            cif_path = temp_path / "external_guest_framework.cif"
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
                    output_dir=temp_path / "external_guest_out",
                    eqeq_settings=EqeqChargeSettings(),
                    isotherm_settings=GraspaIsothermSettings(
                        component="methane",
                        guest_bundles=(str(bundle_path),),
                        pressures=(100000.0,),
                    ),
                    graspa_timeout_seconds=30.0,
                )

            self.assertEqual(result.isotherm_settings.component, "CH4")
            point = result.point_results[0]
            run_dir = Path(point.pressure_run_dir)
            self.assertTrue((run_dir / "CH4.def").is_file())
            simulation_input = Path(point.simulation_input_path).read_text(encoding="utf-8")
            self.assertIn("Component 0 MoleculeName             CH4", simulation_input)
            self.assertNotIn("RotationProbability", simulation_input)
            self.assertAlmostEqual(point.loading_mol_per_kg, 0.25)

            pseudo_atoms_text = (run_dir / "pseudo_atoms.def").read_text(encoding="utf-8")
            self.assertIn("#number of pseudo atoms\n24", pseudo_atoms_text)
            self.assertIn("C_ch4", pseudo_atoms_text)
            mixing_rules_text = (run_dir / "force_field_mixing_rules.def").read_text(encoding="utf-8")
            self.assertIn("C_ch4          lennard-jones   148.0000   3.73000", mixing_rules_text)
            metadata = json.loads((run_dir / "guest_forcefield_metadata.json").read_text(encoding="utf-8"))
            self.assertEqual(metadata["guests"][0]["source"], "guest_bundle")
            self.assertEqual(metadata["guests"][0]["parameter_family"], "dreiding")
            self.assertEqual(metadata["guests"][0]["parameter_source"], "test methane model")

            report = json.loads(Path(result.report_path).read_text(encoding="utf-8"))
            self.assertEqual(report["isotherm_settings"]["component"], "CH4")
            self.assertEqual(report["isotherm_settings"]["guest_bundles"], [str(bundle_path)])
            self.assertEqual(
                report["guest_forcefield_metadata_paths"],
                [str(run_dir / "guest_forcefield_metadata.json")],
            )

    def test_run_graspa_widom_workflow_preserves_stacking_suffix_in_eqeq_outputs(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            fake_eqeq = self._write_fake_eqeq_binary(temp_path / "eqeq_fake", strip_leading_cofid_comment=True)
            fake_graspa = self._write_fake_graspa_binary(temp_path / "graspa_fake")
            cif_path = temp_path / "stacked_framework.cif"
            cofid = "3:amine:Nc1ccc(-c2cc(-c3ccc(N)cc3)cc(-c3ccc(N)cc3)c2)cc1.2:aldehyde:O=Cc1ccc(C=O)cc1&&hcb&&imine"
            cif_path.write_text(
                f"# COFid: {cofid} stacking=AA\n"
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
                    output_dir=temp_path / "stacked_widom_out",
                    eqeq_settings=EqeqChargeSettings(),
                    widom_settings=GraspaWidomSettings(components=("Xe",)),
                    graspa_timeout_seconds=30.0,
                )

            self.assertEqual(
                Path(result.eqeq_charged_cif).read_text(encoding="utf-8").splitlines()[0],
                f"# COFid: {cofid} stacking=AA",
            )
            self.assertEqual(
                Path(result.widom_framework_cif).read_text(encoding="utf-8").splitlines()[0],
                f"# COFid: {cofid} stacking=AA",
            )

    def test_packaged_guest_forcefield_metadata_declares_native_and_compatible_families(self):
        metadata = {entry.name: entry for entry in PACKAGED_GUEST_FORCEFIELD_METADATA}

        self.assertEqual(set(metadata), {"TIP4P", "CO2", "H2", "N2", "SO2", "Xe", "Kr"})
        for component in ("TIP4P", "CO2", "H2", "N2", "SO2"):
            self.assertEqual(metadata[component].parameter_family, "dreiding")
            self.assertEqual(metadata[component].compatible_framework_forcefields, ("dreiding",))
        for component in ("Xe", "Kr"):
            self.assertEqual(metadata[component].parameter_family, "uff")
            self.assertEqual(metadata[component].parameter_source, "RASPA GenericMOFs force field.")
            self.assertEqual(metadata[component].compatible_framework_forcefields, ("uff", "dreiding"))
            self.assertIn("distinct from", metadata[component].provenance_note)

    def test_run_graspa_widom_workflow_rejects_dreiding_guest_with_uff_framework(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            cif_path = Path(temp_dir) / "uff_incompatible_framework.cif"
            cif_path.write_text("data_example\n", encoding="utf-8")

            with self.assertRaises(ValueError) as raised:
                run_graspa_widom_workflow(
                    cif_path,
                    widom_settings=GraspaWidomSettings(components=("CO2",), forcefield="uff"),
                )

        self.assertIn("CO2", str(raised.exception))
        self.assertIn("not registered as compatible", str(raised.exception))

    def test_run_graspa_isotherm_rejects_incompatible_external_guest_bundle(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            cif_path = temp_path / "uff_external_framework.cif"
            cif_path.write_text("data_example\n", encoding="utf-8")
            bundle_path = self._write_guest_bundle(temp_path / "ch4_bundle.json")

            with self.assertRaises(ValueError) as raised:
                run_graspa_isotherm_workflow(
                    cif_path,
                    isotherm_settings=GraspaIsothermSettings(
                        component="CH4",
                        guest_bundles=(str(bundle_path),),
                        forcefield="uff",
                    ),
                )

        self.assertIn("CH4", str(raised.exception))
        self.assertIn("not registered as compatible", str(raised.exception))

    def test_run_graspa_mixture_rejects_incompatible_packaged_guest(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            cif_path = Path(temp_dir) / "uff_mixture_framework.cif"
            cif_path.write_text("data_example\n", encoding="utf-8")

            with self.assertRaises(ValueError) as raised:
                run_graspa_mixture_workflow(
                    cif_path,
                    mixture_settings=GraspaMixtureSettings(
                        components=(
                            GraspaMixtureComponentSettings(component="CO2", mol_fraction=0.5),
                            GraspaMixtureComponentSettings(component="Xe", mol_fraction=0.5),
                        ),
                        forcefield="uff",
                    ),
                )

        self.assertIn("CO2", str(raised.exception))
        self.assertIn("not registered as compatible", str(raised.exception))

    def test_run_graspa_widom_workflow_supports_uff_forcefield_assets_for_xe(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            fake_eqeq = self._write_fake_eqeq_binary(temp_path / "eqeq_fake", strip_leading_cofid_comment=True)
            fake_graspa = self._write_fake_graspa_binary(temp_path / "graspa_fake")
            cif_path = temp_path / "uff_framework.cif"
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
                    output_dir=temp_path / "uff_widom_out",
                    eqeq_settings=EqeqChargeSettings(),
                    widom_settings=GraspaWidomSettings(components=("Xe",), forcefield="uff"),
                    graspa_timeout_seconds=30.0,
                )

            mixing_rules_text = (Path(result.widom_run_dir) / "force_field_mixing_rules.def").read_text(encoding="utf-8")
            self.assertEqual(result.widom_settings.forcefield, "uff")
            self.assertIn("// UFF from bundled Open Babel UFF.prm", mixing_rules_text)
            self.assertIn("C              lennard-jones    52.8384", mixing_rules_text)
            self.assertIn("Br             lennard-jones   126.3089", mixing_rules_text)
            metadata = json.loads(
                (Path(result.widom_run_dir) / "guest_forcefield_metadata.json").read_text(encoding="utf-8")
            )
            self.assertEqual(metadata["guests"][0]["parameter_family"], "uff")
            self.assertEqual(metadata["guests"][0]["parameter_source"], "RASPA GenericMOFs force field.")
            self.assertEqual(metadata["guests"][0]["compatible_framework_forcefields"], ["uff", "dreiding"])
            self.assertIn("distinct from", metadata["guests"][0]["provenance_note"])

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
                            "--forcefield",
                            "uff-openbabel-3.1.0-cofkit-1.0",
                            "--component",
                            "Xe",
                            "--component",
                            "Kr",
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
            self.assertEqual(
                report["widom_settings"]["forcefield"],
                "uff-openbabel-3.1.0-cofkit-1.0",
            )
            self.assertEqual(report["framework_forcefield_metadata"]["family"], "uff")
            self.assertEqual(len(report["component_results"]), 2)
            self.assertEqual(report["component_results"][0]["component"], "Xe")
            self.assertEqual(report["component_results"][1]["component"], "Kr")

    def test_run_raspa2_widom_workflow_prepares_inputs_and_parses_results(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            fake_eqeq = self._write_fake_eqeq_binary(temp_path / "eqeq_fake")
            fake_raspa2 = self._write_fake_raspa2_binary(temp_path / "simulate")
            cif_path = temp_path / "raspa2_widom_framework.cif"
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
                    COFKIT_RASPA2_ENV_VAR: str(fake_raspa2),
                },
                clear=False,
            ):
                result = run_graspa_widom_workflow(
                    cif_path,
                    output_dir=temp_path / "raspa2_widom_out",
                    eqeq_settings=EqeqChargeSettings(),
                    widom_settings=GraspaWidomSettings(components=("CO2", "Xe"), backend="raspa2"),
                    graspa_timeout_seconds=30.0,
                )

            self.assertEqual(result.raspa_backend, "raspa2")
            self.assertEqual(result.unit_cells, (1, 2, 3))
            self.assertEqual([row.component for row in result.component_results], ["CO2", "Xe"])
            self.assertAlmostEqual(result.component_results[0].henry, 1.0e-5)
            self.assertTrue(Path(result.data_file_paths[0]).parent.name == "System_0")

            simulation_input = Path(result.simulation_input_path).read_text(encoding="utf-8")
            self.assertIn("SimulationType                MonteCarlo", simulation_input)
            self.assertIn("NumberOfCycles                2000000", simulation_input)
            self.assertIn("Forcefield                    COFKit", simulation_input)
            self.assertIn("Framework 0", simulation_input)
            self.assertIn("UnitCells 1 2 3", simulation_input)
            self.assertIn("ExternalTemperature 300", simulation_input)
            self.assertIn("MoleculeDefinition        COFKit", simulation_input)
            self.assertNotIn("UseGPUReduction", simulation_input)
            self.assertNotIn("UnitCells 0", simulation_input)

            report = json.loads(Path(result.report_path).read_text(encoding="utf-8"))
            self.assertEqual(report["raspa_backend"], "raspa2")
            self.assertEqual(report["widom_settings"]["backend"], "raspa2")
            self.assertEqual(report["raspa_binary"], str(fake_raspa2.resolve()))

    def test_run_graspa_isotherm_workflow_prepares_inputs_and_parses_results(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            fake_eqeq = self._write_fake_eqeq_binary(temp_path / "eqeq_fake", strip_leading_cofid_comment=True)
            fake_graspa = self._write_fake_graspa_binary(temp_path / "graspa_fake")
            cif_path = temp_path / "adsorption_framework.cif"
            cofid = "3:aldehyde:O=Cc1cc(C=O)cc(C=O)c1.2:hydrazide:CCOc1cc(C(=O)NN)cc(C(=O)NN)c1OCC&&hcb&&hydrazone"
            cif_path.write_text(
                f"# COFid: {cofid}\n"
                "data_example\n"
                "_cell_length_a 26.0\n"
                "_cell_length_b 13.0\n"
                "_cell_length_c 9.0\n",
                encoding="utf-8",
            )
            initial_restart_file = temp_path / "restartfile"
            initial_restart_file.write_text("initial restart sentinel\n", encoding="utf-8")

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
                    initial_restart_file=initial_restart_file,
                    graspa_timeout_seconds=30.0,
                )

            self.assertEqual(result.unit_cells, (1, 2, 3))
            self.assertEqual(len(result.point_results), 2)
            self.assertTrue(Path(result.eqeq_charged_cif).is_file())
            self.assertTrue(Path(result.results_csv_path).is_file())
            self.assertTrue(Path(result.report_path).is_file())
            self.assertTrue(Path(result.isotherm_framework_cif).is_file())
            self.assertEqual(Path(result.eqeq_charged_cif).read_text(encoding="utf-8").splitlines()[0], f"# COFid: {cofid}")
            self.assertEqual(
                Path(result.isotherm_framework_cif).read_text(encoding="utf-8").splitlines()[0],
                f"# COFid: {cofid}",
            )

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
            self.assertIsNotNone(first_point.initial_restart_file_path)
            self.assertEqual(
                Path(first_point.initial_restart_file_path).read_text(encoding="utf-8"),
                "initial restart sentinel\n",
            )
            self.assertEqual(Path(first_point.initial_restart_file_path).name, "restartfile")
            self.assertEqual(Path(first_point.initial_restart_file_path).parent.parts[-2:], ("RestartInitial", "System_0"))

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
            self.assertTrue(Path(report["point_results"][0]["initial_restart_file_path"]).is_file())

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

    def test_calculate_graspa_isotherm_cli_supports_raspa2_backend(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            fake_eqeq = self._write_fake_eqeq_binary(temp_path / "eqeq_fake")
            fake_raspa2 = self._write_fake_raspa2_binary(temp_path / "simulate")
            cif_path = temp_path / "cli_raspa2_adsorption_example.cif"
            cif_path.write_text(
                "data_example\n"
                "_cell_length_a 26.0\n"
                "_cell_length_b 13.0\n"
                "_cell_length_c 9.0\n",
                encoding="utf-8",
            )
            output_dir = temp_path / "cli_raspa2_isotherm"
            buffer = io.StringIO()

            with patch.dict(
                os.environ,
                {COFKIT_EQEQ_ENV_VAR: str(fake_eqeq)},
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
                            "--backend",
                            "raspa2",
                            "--raspa2-path",
                            str(fake_raspa2),
                            "--component",
                            "CO2",
                            "--pressure",
                            "100000",
                            "--json",
                        ]
                    )

            report = json.loads(buffer.getvalue())
            self.assertEqual(report["raspa_backend"], "raspa2")
            self.assertEqual(report["isotherm_settings"]["backend"], "raspa2")
            self.assertEqual(report["point_results"][0]["loading_mol_per_kg"], 0.25)
            self.assertIsNone(report["point_results"][0]["loading_g_per_l"])
            self.assertEqual(report["point_results"][0]["heat_of_adsorption_kj_per_mol"], -21.0)

            simulation_input = Path(report["point_results"][0]["simulation_input_path"]).read_text(encoding="utf-8")
            self.assertIn("ExternalPressure 100000", simulation_input)
            self.assertIn("NumberOfCycles                200000", simulation_input)
            self.assertIn("MoleculeDefinition       COFKit", simulation_input)
            self.assertNotIn("Pressure     100000", simulation_input)

    def test_run_graspa_mixture_workflow_prepares_inputs_and_parses_selectivity(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            fake_eqeq = self._write_fake_eqeq_binary(temp_path / "eqeq_fake", strip_leading_cofid_comment=True)
            fake_graspa = self._write_fake_graspa_binary(temp_path / "graspa_fake")
            cif_path = temp_path / "mixture_framework.cif"
            cofid = "3:amine:Nc1ccc(-c2ccc(N)cc2)cc1.2:aldehyde:O=Cc1ccc(C=O)cc1&&sql&&imine"
            cif_path.write_text(
                f"# COFid: {cofid}\n"
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
                result = run_graspa_mixture_workflow(
                    cif_path,
                    output_dir=temp_path / "mixture_out",
                    eqeq_settings=EqeqChargeSettings(),
                    mixture_settings=GraspaMixtureSettings(
                        components=(
                            GraspaMixtureComponentSettings("Kr", 0.1, fugacity_coefficient="PR-EOS"),
                            GraspaMixtureComponentSettings("Xe", 0.9, fugacity_coefficient="PR-EOS"),
                        ),
                        pressures=(10000.0, 100000.0),
                    ),
                    graspa_timeout_seconds=30.0,
                )

            self.assertEqual(result.unit_cells, (1, 2, 3))
            self.assertEqual(len(result.point_results), 2)
            self.assertTrue(Path(result.eqeq_charged_cif).is_file())
            self.assertTrue(Path(result.component_results_csv_path).is_file())
            self.assertTrue(Path(result.selectivity_results_csv_path).is_file())
            self.assertTrue(Path(result.report_path).is_file())
            self.assertEqual(Path(result.eqeq_charged_cif).read_text(encoding="utf-8").splitlines()[0], f"# COFid: {cofid}")
            self.assertEqual(
                Path(result.mixture_framework_cif).read_text(encoding="utf-8").splitlines()[0],
                f"# COFid: {cofid}",
            )

            first_point = result.point_results[0]
            self.assertEqual(first_point.pressure, 10000.0)
            self.assertEqual(first_point.component_results[0].component, "Kr")
            self.assertAlmostEqual(first_point.component_results[0].loading_mol_per_kg, 0.01)
            self.assertAlmostEqual(first_point.component_results[1].loading_mol_per_kg, 0.18)
            self.assertAlmostEqual(first_point.component_results[0].adsorbed_mol_fraction, 0.01 / 0.19)
            self.assertAlmostEqual(first_point.component_results[1].adsorbed_mol_fraction, 0.18 / 0.19)
            self.assertAlmostEqual(first_point.selectivity_results[0].selectivity, 0.5)
            self.assertAlmostEqual(first_point.selectivity_results[1].selectivity, 2.0)
            self.assertAlmostEqual(first_point.selectivity_results[1].selectivity_errorbar, 2.0 * (2.0 ** 0.5) * 0.1)

            simulation_input = Path(first_point.simulation_input_path).read_text(encoding="utf-8")
            self.assertIn("MolFraction              0.1", simulation_input)
            self.assertIn("MolFraction              0.9", simulation_input)
            self.assertIn("IdentityChangeProbability 1", simulation_input)
            self.assertIn("FugacityCoefficient      PR-EOS", simulation_input)
            self.assertNotIn("RotationProbability", simulation_input)

            component_csv_text = Path(result.component_results_csv_path).read_text(encoding="utf-8")
            self.assertIn("adsorbed_mol_fraction", component_csv_text)
            self.assertIn("Kr", component_csv_text)
            self.assertIn("Xe", component_csv_text)
            selectivity_csv_text = Path(result.selectivity_results_csv_path).read_text(encoding="utf-8")
            self.assertIn("selectivity", selectivity_csv_text)
            self.assertIn("Kr", selectivity_csv_text)
            self.assertIn("Xe", selectivity_csv_text)

            report = json.loads(Path(result.report_path).read_text(encoding="utf-8"))
            self.assertEqual(report["unit_cells"], [1, 2, 3])
            self.assertEqual(len(report["point_results"]), 2)
            self.assertEqual(report["point_results"][0]["component_results"][0]["feed_mol_fraction"], 0.1)
            self.assertEqual(report["point_results"][0]["selectivity_results"][1]["selectivity"], 2.0)

    def test_calculate_graspa_mixture_cli_prints_json_report(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            fake_eqeq = self._write_fake_eqeq_binary(temp_path / "eqeq_fake")
            fake_graspa = self._write_fake_graspa_binary(temp_path / "graspa_fake")
            cif_path = temp_path / "cli_mixture_example.cif"
            cif_path.write_text(
                "data_example\n"
                "_cell_length_a 26.0\n"
                "_cell_length_b 13.0\n"
                "_cell_length_c 9.0\n",
                encoding="utf-8",
            )
            output_dir = temp_path / "cli_mixture"
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
                            "graspa-mixture",
                            str(cif_path),
                            "--output-dir",
                            str(output_dir),
                            "--component",
                            "CO2:0.25",
                            "--component",
                            "N2:0.75",
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
            self.assertEqual(len(report["mixture_settings"]["components"]), 2)
            self.assertEqual(report["mixture_settings"]["components"][0]["component"], "CO2")
            self.assertEqual(report["mixture_settings"]["components"][1]["component"], "N2")
            self.assertEqual(len(report["point_results"]), 2)
            self.assertEqual(len(report["point_results"][0]["component_results"]), 2)
            self.assertEqual(len(report["point_results"][0]["selectivity_results"]), 2)
            self.assertAlmostEqual(report["point_results"][0]["selectivity_results"][1]["selectivity"], 2.0)

            simulation_input = Path(report["point_results"][0]["simulation_input_path"]).read_text(encoding="utf-8")
            self.assertIn("MolFraction              0.25", simulation_input)
            self.assertIn("MolFraction              0.75", simulation_input)
            self.assertIn("IdentityChangeProbability 1", simulation_input)
            self.assertIn("RotationProbability      1", simulation_input)
            self.assertIn("FugacityCoefficient      PR-EOS", simulation_input)

    def test_run_raspa2_mixture_workflow_prepares_inputs_and_parses_selectivity(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            fake_eqeq = self._write_fake_eqeq_binary(temp_path / "eqeq_fake")
            fake_raspa2 = self._write_fake_raspa2_binary(temp_path / "simulate")
            cif_path = temp_path / "raspa2_mixture_framework.cif"
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
                    COFKIT_RASPA2_ENV_VAR: str(fake_raspa2),
                },
                clear=False,
            ):
                result = run_graspa_mixture_workflow(
                    cif_path,
                    output_dir=temp_path / "raspa2_mixture_out",
                    eqeq_settings=EqeqChargeSettings(),
                    mixture_settings=GraspaMixtureSettings(
                        components=(
                            GraspaMixtureComponentSettings("Kr", 0.1),
                            GraspaMixtureComponentSettings("Xe", 0.9),
                        ),
                        pressures=(100000.0,),
                        backend="raspa2",
                    ),
                    graspa_timeout_seconds=30.0,
                )

            self.assertEqual(result.raspa_backend, "raspa2")
            first_point = result.point_results[0]
            self.assertAlmostEqual(first_point.component_results[0].loading_mol_per_kg, 0.1)
            self.assertAlmostEqual(first_point.component_results[1].loading_mol_per_kg, 1.8)
            self.assertAlmostEqual(first_point.selectivity_results[0].selectivity, 0.5)
            self.assertTrue(math.isnan(first_point.component_results[0].loading_g_per_l))
            self.assertTrue(math.isnan(first_point.component_results[0].heat_of_adsorption_kj_per_mol))

            simulation_input = Path(first_point.simulation_input_path).read_text(encoding="utf-8")
            self.assertIn("ExternalPressure 100000", simulation_input)
            self.assertIn("MolFraction              0.1", simulation_input)
            self.assertIn("MolFraction              0.9", simulation_input)
            self.assertIn("IdentityChangeProbability 1", simulation_input)
            self.assertNotIn("UseGPUReduction", simulation_input)

    def test_calculate_help_lists_graspa_widom(self):
        buffer = io.StringIO()
        with self.assertRaises(SystemExit), contextlib.redirect_stdout(buffer):
            cli_main(["calculate", "--help"])

        self.assertIn("graspa-widom", buffer.getvalue())
        self.assertIn("graspa-isotherm", buffer.getvalue())
        self.assertIn("graspa-mixture", buffer.getvalue())
        self.assertIn("hybrid-mdmc", buffer.getvalue())

    def _write_guest_bundle(self, path: Path) -> Path:
        path.write_text(
            json.dumps(
                {
                    "version": 2,
                    "name": "CH4",
                    "aliases": ["methane"],
                    "rotatable": False,
                    "parameter_family": "dreiding",
                    "parameter_source": "test methane model",
                    "compatible_framework_forcefields": ["dreiding"],
                    "raspa": {
                        "molecule_definition": (
                            "# critical constants: Temperature [T], Pressure [Pa], and Acentric factor [-]\n"
                            "190.564\n"
                            "4599000.0\n"
                            "0.011\n"
                            "#Number Of Atoms\n"
                            " 1\n"
                            "# Number of groups\n"
                            "1\n"
                            "# methane-group\n"
                            "rigid\n"
                            "# number of atoms\n"
                            "1\n"
                            "# atomic positions\n"
                            "0 C_ch4     0.0 0.0 0.0\n"
                            "# Chiral centers Bond  BondDipoles Bend  UrayBradley InvBend  Torsion Imp. Torsion Bond/Bond Stretch/Bend Bend/Bend Stretch/Torsion Bend/Torsion IntraVDW IntraCoulomb\n"
                            "               0    0            0    0            0       0        0            0         0            0         0               0            0        0            0\n"
                            "# Number of config moves\n"
                            "0\n"
                        ),
                        "pseudo_atom_rows": [
                            "C_ch4      yes     C     C     0          16.04300    0.0      0.0          1.0      0.720  0            0           relative           0",
                        ],
                        "mixing_rule_rows": [
                            "C_ch4          lennard-jones   148.0000   3.73000      // test methane guest bundle",
                        ],
                    },
                    "lammps": {
                        "units": "real",
                        "atom_style": "full",
                        "site_order": ["C_ch4"],
                        "pair_style": "lj/cut/coul/long",
                        "pair_coeff_rows": [
                            "C_ch4 lennard-jones 148.0000 3.73000",
                        ],
                    },
                },
                indent=2,
            ),
            encoding="utf-8",
        )
        return path

    def _write_fake_eqeq_binary(self, path: Path, *, strip_leading_cofid_comment: bool = False) -> Path:
        path.write_text(
            f"#!{sys.executable}\n"
            "from __future__ import annotations\n"
            "import json\n"
            "import sys\n"
            "from pathlib import Path\n"
            "\n"
            f"strip_leading_cofid_comment = {strip_leading_cofid_comment!r}\n"
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
            "if strip_leading_cofid_comment and input_text.startswith('# COFid: '):\n"
            "    _, _, input_text = input_text.partition('\\n')\n"
            "Path(stem + '.cif').write_text(input_text + '\\n# fake charged output\\n', encoding='utf-8')\n"
            "Path(stem + '.json').write_text(json.dumps({'argv': args}, indent=2), encoding='utf-8')\n"
            "Path('eqeq_invocation.json').write_text(json.dumps({'argv': args}, indent=2), encoding='utf-8')\n"
            "sys.stdout.write('fake eqeq stdout\\n')\n"
            "sys.stderr.write('fake eqeq stderr\\n')\n",
            encoding="utf-8",
        )
        path.chmod(path.stat().st_mode | stat.S_IEXEC)
        return path

    def _write_fake_raspa2_binary(self, path: Path) -> Path:
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
            "if 'SimulationType                MonteCarlo' not in text:\n"
            "    sys.stderr.write('missing RASPA2 simulation type\\n')\n"
            "    sys.exit(3)\n"
            "framework_match = re.search(r'^FrameworkName\\s+(\\S+)', text, re.MULTILINE)\n"
            "unit_cells_match = re.search(r'^UnitCells\\s+(\\d+)\\s+(\\d+)\\s+(\\d+)', text, re.MULTILINE)\n"
            "pressure_match = re.search(r'^ExternalPressure\\s+(\\S+)', text, re.MULTILINE)\n"
            "components = re.findall(r'^Component\\s+\\d+\\s+MoleculeName\\s+(\\S+)', text, re.MULTILINE)\n"
            "mol_fractions = [float(value) for value in re.findall(r'^\\s*MolFraction\\s+(\\S+)', text, re.MULTILINE)]\n"
            "if framework_match is None or unit_cells_match is None or not components:\n"
            "    sys.stderr.write('missing framework/unitcells/components config\\n')\n"
            "    sys.exit(4)\n"
            "framework_name = framework_match.group(1)\n"
            "pressure = float(pressure_match.group(1)) if pressure_match is not None else 0.0\n"
            "if 'WidomProbability' in text:\n"
            "    mode = 'widom'\n"
            "elif mol_fractions and len(components) > 1:\n"
            "    mode = 'mixture'\n"
            "else:\n"
            "    mode = 'isotherm'\n"
            "framework_cif = Path(framework_name + '.cif')\n"
            "if not framework_cif.is_file():\n"
            "    sys.stderr.write('missing framework cif\\n')\n"
            "    sys.exit(5)\n"
            "unit_cells = [int(unit_cells_match.group(1)), int(unit_cells_match.group(2)), int(unit_cells_match.group(3))]\n"
            "Path('raspa2_invocation.json').write_text(\n"
            "    json.dumps(\n"
            "        {\n"
            "            'mode': mode,\n"
            "            'framework_name': framework_name,\n"
            "            'unit_cells': unit_cells,\n"
            "            'pressure': pressure,\n"
            "            'components': components,\n"
            "            'mol_fractions': mol_fractions,\n"
            "        },\n"
            "        indent=2,\n"
            "    ),\n"
            "    encoding='utf-8',\n"
            ")\n"
            "output_dir = Path('Output') / 'System_0'\n"
            "output_dir.mkdir(parents=True, exist_ok=True)\n"
            "data_path = output_dir / f'output_{framework_name}_{unit_cells[0]}.{unit_cells[1]}.{unit_cells[2]}_fake.data'\n"
            "if mode == 'widom':\n"
            "    sections = ['Average Widom excess contribution:', '==================================']\n"
            "    for index, component in enumerate(components):\n"
            "        energy = -10.0 - index\n"
            "        energy_err = 0.1 + index * 0.01\n"
            "        sections.append(f'[{component}] Average Widom excess chemical potential:   {energy:.2f} +/- {energy_err:.6f} [-]')\n"
            "    sections.extend(['', 'Average Henry coefficient:', '=========================='])\n"
            "    for index, component in enumerate(components):\n"
            "        henry = (index + 1) * 1.0e-5\n"
            "        henry_err = (index + 1) * 1.0e-6\n"
            "        sections.append(f'[{component}] Average Henry coefficient:  {henry:.8e} +/- {henry_err:.8e} [mol/kg/Pa]')\n"
            "    data_path.write_text('\\n'.join(sections) + '\\n', encoding='utf-8')\n"
            "elif mode == 'mixture':\n"
            "    sections = ['Number of molecules:', '====================', '']\n"
            "    for index, (component, mol_fraction) in enumerate(zip(components, mol_fractions), start=1):\n"
            "        loading_mol = pressure / 100000.0 * mol_fraction * index\n"
            "        loading_err = loading_mol * 0.1\n"
            "        sections.extend([\n"
            "            f'Component {index - 1} [{component}]',\n"
            "            '-------------------------------------------------------------',\n"
            "            f'Average loading absolute [mol/kg framework]       {loading_mol:.10f} +/-       {loading_err:.10f} [-]',\n"
            "            '',\n"
            "        ])\n"
            "    sections.append('Average Widom Rosenbluth factor:')\n"
            "    data_path.write_text('\\n'.join(sections) + '\\n', encoding='utf-8')\n"
            "else:\n"
            "    loading_mol = pressure / 100000.0 * 0.25\n"
            "    loading_err = loading_mol * 0.1\n"
            "    heat = -20.0 - pressure / 100000.0\n"
            "    data_path.write_text(\n"
            "        'Number of molecules:\\n'\n"
            "        '====================\\n\\n'\n"
            "        f'Average loading absolute [mol/kg framework]       {loading_mol:.10f} +/-       {loading_err:.10f} [-]\\n\\n'\n"
            "        'Enthalpy of adsorption:\\n'\n"
            "        '=======================\\n\\n'\n"
            "        '\\tTotal enthalpy of adsorption\\n'\n"
            "        '\\t----------------------------\\n'\n"
            "        f'\\tAverage          {-3200.0:.5f} +/-          {52.0:.6f} [K]\\n'\n"
            "        f'\\t                   {heat:.5f} +/-           {0.5:.6f} [KJ/MOL]\\n',\n"
            "        encoding='utf-8',\n"
            "    )\n"
            "sys.stdout.write('fake raspa2 stdout\\n')\n",
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
            "mol_fractions = [float(value) for value in re.findall(r'^\\s*MolFraction\\s+(\\S+)', text, re.MULTILINE)]\n"
            "if framework_match is None or unit_cells_match is None or pressure_match is None:\n"
            "    sys.stderr.write('missing framework/unitcells config\\n')\n"
            "    sys.exit(3)\n"
            "\n"
            "framework_name = framework_match.group(1)\n"
            "pressure = float(pressure_match.group(1))\n"
            "if 'WidomProbability' in text:\n"
            "    mode = 'widom'\n"
            "elif mol_fractions and len(components) > 1:\n"
            "    mode = 'mixture'\n"
            "else:\n"
            "    mode = 'isotherm'\n"
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
            "            'mol_fractions': mol_fractions,\n"
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
            "elif mode == 'mixture':\n"
            "    sections = [\n"
            "        '================= MOL FRACTIONS =================',\n"
            "    ]\n"
            "    for index, (component, mol_fraction) in enumerate(zip(components, mol_fractions), start=1):\n"
            "        sections.append(f'Component [{index}] ({component}), Mol Fraction: {mol_fraction:.5f}')\n"
            "    sections.append('=================================================')\n"
            "    sections.append('==============================================================')\n"
            "    sections.append('============= BLOCK AVERAGES (HEAT OF ADSORPTION: kJ/mol) =========')\n"
            "    for index, (component, mol_fraction) in enumerate(zip(components, mol_fractions), start=1):\n"
            "        heat = -15.0 - index - pressure / 100000.0\n"
            "        heat_err = 0.5\n"
            "        sections.extend([\n"
            "            f'COMPONENT [{index}] ({component})',\n"
            "            f'Overall: Average: {heat:.5f}, ErrorBar: {heat_err:.5f}',\n"
            "            '-----------------------------',\n"
            "        ])\n"
            "    sections.append('==============================================================')\n"
            "    sections.append('=====================BLOCK AVERAGES (LOADING: mol/kg)=============')\n"
            "    for index, (component, mol_fraction) in enumerate(zip(components, mol_fractions), start=1):\n"
            "        loading_mol = pressure / 100000.0 * mol_fraction * index\n"
            "        loading_err = loading_mol * 0.1\n"
            "        sections.extend([\n"
            "            f'COMPONENT [{index}] ({component})',\n"
            "            f'Overall: Average: {loading_mol:.5f}, ErrorBar: {loading_err:.5f}',\n"
            "            '----------------------------------------------------------',\n"
            "        ])\n"
            "    sections.append('==============================================================')\n"
            "    sections.append('=====================BLOCK AVERAGES (LOADING: g/L)=============')\n"
            "    for index, (component, mol_fraction) in enumerate(zip(components, mol_fractions), start=1):\n"
            "        loading_mol = pressure / 100000.0 * mol_fraction * index\n"
            "        loading_err = loading_mol * 0.1\n"
            "        loading_gpl = loading_mol * 10.0\n"
            "        loading_gpl_err = loading_err * 10.0\n"
            "        sections.extend([\n"
            "            f'COMPONENT [{index}] ({component})',\n"
            "            f'Overall: Average: {loading_gpl:.5f}, ErrorBar: {loading_gpl_err:.5f}',\n"
            "            '----------------------------------------------------------',\n"
            "        ])\n"
            "    data_path.write_text('\\n'.join(sections) + '\\n', encoding='utf-8')\n"
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
