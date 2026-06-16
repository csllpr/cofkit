import json
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from cofkit.graspa import (
    EqeqChargeSettings,
    GraspaIsothermPointResult,
    GraspaIsothermResult,
    GraspaIsothermSettings,
    GraspaMixtureComponentSettings,
)
from cofkit.hybrid_mdmc import HybridMdMcSettings, run_hybrid_mdmc_workflow
from cofkit.lammps import LammpsMdResult, LammpsMdSettings


class HybridMdMcTests(unittest.TestCase):
    def test_hybrid_mdmc_cycles_framework_snapshots_between_md_and_gcmc(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            input_cif = temp_path / "framework.cif"
            input_cif.write_text("data_framework\n", encoding="utf-8")
            lammps_inputs: list[Path] = []
            gcmc_inputs: list[Path] = []

            def fake_lammps(cif_path, *, output_dir, settings, **kwargs):
                input_path = Path(cif_path)
                lammps_inputs.append(input_path)
                run_dir = Path(output_dir)
                run_dir.mkdir(parents=True, exist_ok=True)
                output_cif = run_dir / f"{input_path.stem}_md.cif"
                output_cif.write_text(f"data_cycle_{len(lammps_inputs)}\n", encoding="utf-8")
                return _fake_lammps_result(input_path, run_dir, output_cif, settings)

            def fake_isotherm(cif_path, *, output_dir, isotherm_settings, **kwargs):
                input_path = Path(cif_path)
                gcmc_inputs.append(input_path)
                run_dir = Path(output_dir)
                run_dir.mkdir(parents=True, exist_ok=True)
                return _fake_isotherm_result(input_path, run_dir, isotherm_settings.component)

            with patch("cofkit.hybrid_mdmc.run_lammps_md_on_cif", side_effect=fake_lammps):
                with patch("cofkit.hybrid_mdmc.run_graspa_isotherm_workflow", side_effect=fake_isotherm):
                    result = run_hybrid_mdmc_workflow(
                        input_cif,
                        output_dir=temp_path / "hybrid_out",
                        settings=HybridMdMcSettings(
                            cycles=2,
                            pressure=100000.0,
                            components=(GraspaMixtureComponentSettings(component="CO2", mol_fraction=1.0),),
                            initialization_cycles=1,
                            equilibration_cycles=1,
                            production_cycles=2,
                        ),
                        lammps_md_settings=LammpsMdSettings(
                            forcefield="uff",
                            charge_model="none",
                            steps=5,
                        ),
                        lammps_eqeq_settings=EqeqChargeSettings(),
                        raspa_eqeq_settings=EqeqChargeSettings(),
                    )

            self.assertEqual(lammps_inputs[0], input_cif.resolve())
            self.assertEqual(lammps_inputs[1], Path(result.cycle_results[0].output_framework_cif))
            self.assertEqual(gcmc_inputs[0], Path(result.cycle_results[0].output_framework_cif))
            self.assertEqual(gcmc_inputs[1], Path(result.cycle_results[1].output_framework_cif))
            self.assertEqual(result.final_framework_cif, result.cycle_results[1].output_framework_cif)
            self.assertEqual(result.cycle_results[0].gcmc_result_type, "isotherm")
            self.assertIn("exchange_mode='framework'", result.warnings[0])

            report = json.loads(Path(result.report_path).read_text(encoding="utf-8"))
            self.assertEqual(report["settings"]["exchange_mode"], "framework")
            self.assertEqual(len(report["cycle_results"]), 2)
            self.assertEqual(report["cycle_results"][0]["gcmc_result_type"], "isotherm")
            self.assertEqual(report["cycle_results"][1]["lammps_md_result"]["settings"]["steps"], 5)

    def test_hybrid_mdmc_guest_restart_carries_binary_gcmc_snapshot_to_next_md(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            input_cif = temp_path / "framework.cif"
            input_cif.write_text("data_framework\n", encoding="utf-8")
            lammps_guest_states = []
            gcmc_initial_restart_files = []
            gcmc_restart_flags = []
            restart_text = ""

            def fake_lammps(cif_path, *, output_dir, settings, guest_restart_state=None, **kwargs):
                input_path = Path(cif_path)
                lammps_guest_states.append(guest_restart_state)
                run_dir = Path(output_dir)
                run_dir.mkdir(parents=True, exist_ok=True)
                output_cif = run_dir / f"{input_path.stem}_md.cif"
                output_cif.write_text(f"data_cycle_{len(lammps_guest_states)}\n", encoding="utf-8")
                return _fake_lammps_result(
                    input_path,
                    run_dir,
                    output_cif,
                    settings,
                    guest_restart_state=guest_restart_state,
                )

            def fake_mixture(cif_path, *, output_dir, mixture_settings, initial_restart_file=None, **kwargs):
                gcmc_initial_restart_files.append(initial_restart_file)
                gcmc_restart_flags.append(mixture_settings.restart_file)
                run_dir = Path(output_dir)
                pressure_run_dir = run_dir / "mixture" / "pressure_100000"
                movie_dir = pressure_run_dir / "Movies" / "System_0"
                movie_dir.mkdir(parents=True, exist_ok=True)
                (movie_dir / "result_2.data").write_text(
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
                return _FakeMixtureResult(run_dir, pressure_run_dir, mixture_settings)

            with patch("cofkit.hybrid_mdmc.run_lammps_md_on_cif", side_effect=fake_lammps):
                with patch("cofkit.hybrid_mdmc.run_graspa_mixture_workflow", side_effect=fake_mixture):
                    result = run_hybrid_mdmc_workflow(
                        input_cif,
                        output_dir=temp_path / "hybrid_guest_out",
                        settings=HybridMdMcSettings(
                            cycles=2,
                            exchange_mode="guest_restart",
                            pressure=100000.0,
                            components=(
                                GraspaMixtureComponentSettings(component="Xe", mol_fraction=0.5),
                                GraspaMixtureComponentSettings(component="Kr", mol_fraction=0.5),
                            ),
                            initialization_cycles=1,
                            equilibration_cycles=1,
                            production_cycles=2,
                        ),
                        lammps_md_settings=LammpsMdSettings(
                            forcefield="uff",
                            charge_model="none",
                            steps=5,
                        ),
                        lammps_eqeq_settings=EqeqChargeSettings(),
                        raspa_eqeq_settings=EqeqChargeSettings(),
                    )
                    report = json.loads(Path(result.report_path).read_text(encoding="utf-8"))
                    restart_text = Path(gcmc_initial_restart_files[1]).read_text(encoding="utf-8")

        self.assertIsNone(lammps_guest_states[0])
        self.assertEqual(lammps_guest_states[1].n_atoms, 2)
        self.assertEqual(lammps_guest_states[1].components, ("Xe", "Kr"))
        self.assertEqual(result.cycle_results[0].n_output_guest_atoms, 2)
        self.assertEqual(result.cycle_results[1].n_input_guest_atoms, 2)
        self.assertEqual(result.cycle_results[1].lammps_md_result.n_guest_atoms, 2)
        self.assertEqual(result.cycle_results[1].input_guest_components, ("Xe", "Kr"))
        self.assertIsNone(gcmc_initial_restart_files[0])
        self.assertIsNotNone(gcmc_initial_restart_files[1])
        self.assertEqual(gcmc_restart_flags, [False, True])
        self.assertEqual(result.cycle_results[0].gcmc_initial_restart_file_path, None)
        self.assertEqual(result.cycle_results[1].n_md_output_guest_atoms, 2)
        self.assertEqual(result.cycle_results[1].md_output_guest_components, ("Xe", "Kr"))
        self.assertEqual(result.cycle_results[1].gcmc_initial_restart_file_path, gcmc_initial_restart_files[1])
        self.assertIn("Components: 2 (Adsorbates 2, Cations 0)", restart_text)
        self.assertIn("Adsorbate-atom-position: 0 0 1.250000 2.500000 3.750000", restart_text)
        self.assertIn("Adsorbate-atom-position: 0 0 4.250000 5.500000 6.750000", restart_text)
        self.assertIn("exchange_mode='guest_restart'", result.warnings[0])
        self.assertEqual(report["settings"]["exchange_mode"], "guest_restart")
        self.assertEqual(report["cycle_results"][0]["n_output_guest_atoms"], 2)
        self.assertEqual(report["cycle_results"][1]["n_input_guest_atoms"], 2)
        self.assertEqual(report["cycle_results"][1]["n_md_output_guest_atoms"], 2)
        self.assertEqual(report["cycle_results"][1]["gcmc_initial_restart_file_path"], gcmc_initial_restart_files[1])

    def test_hybrid_mdmc_guest_restart_rejects_raspa2_backend(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            input_cif = temp_path / "framework.cif"
            input_cif.write_text("data_framework\n", encoding="utf-8")

            with self.assertRaises(ValueError) as raised:
                run_hybrid_mdmc_workflow(
                    input_cif,
                    output_dir=temp_path / "hybrid_guest_out",
                    settings=HybridMdMcSettings(
                        exchange_mode="guest_restart",
                        raspa_backend="raspa2",
                        components=(GraspaMixtureComponentSettings(component="Xe", mol_fraction=1.0),),
                    ),
                    lammps_md_settings=LammpsMdSettings(forcefield="uff", charge_model="none", steps=5),
                    lammps_eqeq_settings=EqeqChargeSettings(),
                    raspa_eqeq_settings=EqeqChargeSettings(),
                )

        self.assertIn("RASPA2 restart staging is not yet supported", str(raised.exception))


def _fake_lammps_result(
    input_path: Path,
    run_dir: Path,
    output_cif: Path,
    settings: LammpsMdSettings,
    *,
    guest_restart_state=None,
) -> LammpsMdResult:
    _write_fake_lammps_md_files(run_dir, guest_restart_state=guest_restart_state)
    report_path = run_dir / "lammps_md_report.json"
    result = LammpsMdResult(
        input_cif=str(input_path),
        output_cif=str(output_cif),
        output_dir=str(run_dir),
        lammps_binary="/fake/lmp",
        lammps_data_path=str(run_dir / "lammps_md_input.data"),
        lammps_input_script_path=str(run_dir / "lammps_md.in"),
        lammps_dump_path=str(run_dir / "lammps_md_trajectory.lammpstrj"),
        lammps_log_path=str(run_dir / "lammps_md.log"),
        stdout_log_path=str(run_dir / "lammps_md.stdout.log"),
        stderr_log_path=str(run_dir / "lammps_md.stderr.log"),
        report_path=str(report_path),
        n_atoms=3,
        n_bonds=2,
        n_angles=1,
        n_dihedrals=0,
        n_impropers=0,
        n_atom_types=2,
        n_bond_types=2,
        n_angle_types=1,
        n_dihedral_types=0,
        n_improper_types=0,
        atom_type_symbols={1: "C_3", 2: "O_3"},
        settings=settings,
        forcefield_backend="fake",
        parameter_sources={},
        charge_model=settings.charge_model,
        n_charged_atoms=0,
        net_charge=None,
        n_guest_atoms=guest_restart_state.n_atoms if guest_restart_state is not None else 0,
        n_total_atoms=3 + (guest_restart_state.n_atoms if guest_restart_state is not None else 0),
        guest_components=guest_restart_state.components if guest_restart_state is not None else (),
        guest_restart_source_path=(
            guest_restart_state.source_snapshot_path if guest_restart_state is not None else None
        ),
    )
    report_path.write_text(json.dumps(result.to_dict(), indent=2), encoding="utf-8")
    return result


def _write_fake_lammps_md_files(run_dir: Path, *, guest_restart_state=None) -> None:
    data_path = run_dir / "lammps_md_input.data"
    dump_path = run_dir / "lammps_md_trajectory.lammpstrj"
    n_guest_atoms = guest_restart_state.n_atoms if guest_restart_state is not None else 0
    site_by_label = guest_restart_state.site_by_label() if guest_restart_state is not None else {}
    guest_sites = tuple(dict.fromkeys(atom.site_label for atom in guest_restart_state.atoms)) if guest_restart_state else ()
    guest_type_by_site = {site_label: index for index, site_label in enumerate(guest_sites, start=3)}
    mass_rows = [
        "1 12.011 # C_3",
        "2 15.999 # O_3",
        *(
            f"{type_id} {site_by_label[site_label].mass:.8g} # {site_label}"
            for site_label, type_id in guest_type_by_site.items()
        ),
    ]
    atom_rows = [
        "1 1 1 0.0 0.0 0.0 0.0 0 0 0 # C_3",
        "2 1 2 0.0 1.0 0.0 0.0 0 0 0 # O_3",
        "3 1 2 0.0 0.0 1.0 0.0 0 0 0 # O_3",
    ]
    dump_rows = [
        "1 1 0.0 0.0 0.0",
        "2 2 1.0 0.0 0.0",
        "3 2 0.0 1.0 0.0",
    ]
    if guest_restart_state is not None:
        molecule_atoms_by_key = {}
        for atom in guest_restart_state.atoms:
            molecule_atoms_by_key.setdefault(atom.molecule_key, []).append(atom)
        next_atom_id = 4
        next_molecule_id = 2
        for molecule_key, molecule_atoms in molecule_atoms_by_key.items():
            del molecule_key
            molecule_id = next_molecule_id
            next_molecule_id += 1
            for atom in sorted(molecule_atoms, key=lambda item: item.site_index):
                site = site_by_label[atom.site_label]
                atom_type = guest_type_by_site[atom.site_label]
                x = atom.x + 0.25
                y = atom.y + 0.50
                z = atom.z + 0.75
                atom_rows.append(
                    (
                        f"{next_atom_id} {molecule_id} {atom_type} {site.charge:.8g} "
                        f"{x:.10g} {y:.10g} {z:.10g} 0 0 0 # {atom.site_label} {atom.component}"
                    )
                )
                dump_rows.append(f"{next_atom_id} {atom_type} {x:.10g} {y:.10g} {z:.10g}")
                next_atom_id += 1
    data_path.write_text(
        "\n".join(
            [
                "LAMMPS data",
                "",
                f"{3 + n_guest_atoms} atoms",
                f"{2 + len(guest_sites)} atom types",
                "",
                "Masses",
                "",
                *mass_rows,
                "",
                "Atoms # full",
                "",
                *atom_rows,
                "",
            ]
        ),
        encoding="utf-8",
    )
    dump_path.write_text(
        "\n".join(
            [
                "ITEM: TIMESTEP",
                "5",
                "ITEM: NUMBER OF ATOMS",
                str(3 + n_guest_atoms),
                "ITEM: BOX BOUNDS pp pp pp",
                "0 20",
                "0 21",
                "0 22",
                "ITEM: ATOMS id type x y z",
                *dump_rows,
                "",
            ]
        ),
        encoding="utf-8",
    )


class _FakeMixtureResult:
    def __init__(self, run_dir: Path, pressure_run_dir: Path, mixture_settings):
        self.output_dir = str(run_dir)
        self.mixture_root_dir = str(run_dir / "mixture")
        self.point_results = (SimpleNamespace(pressure_run_dir=str(pressure_run_dir)),)
        self.mixture_settings = mixture_settings

    def to_dict(self):
        return {
            "output_dir": self.output_dir,
            "mixture_root_dir": self.mixture_root_dir,
            "point_results": [{"pressure_run_dir": self.point_results[0].pressure_run_dir}],
            "mixture_settings": self.mixture_settings.to_dict(),
        }


def _fake_isotherm_result(input_path: Path, run_dir: Path, component: str) -> GraspaIsothermResult:
    point = GraspaIsothermPointResult(
        component=component,
        pressure=100000.0,
        pressure_run_dir=str(run_dir / "isotherm" / "pressure_100000"),
        simulation_input_path=str(run_dir / "isotherm" / "pressure_100000" / "simulation.input"),
        graspa_stdout_log_path=str(run_dir / "graspa.stdout.log"),
        graspa_stderr_log_path=str(run_dir / "graspa.stderr.log"),
        data_file_paths=(str(run_dir / "Output" / "System_0_fake.data"),),
        source_data_file=str(run_dir / "Output" / "System_0_fake.data"),
        loading_mol_per_kg=0.25,
        loading_mol_per_kg_errorbar=0.01,
        loading_g_per_l=2.5,
        loading_g_per_l_errorbar=0.1,
        heat_of_adsorption_kj_per_mol=-20.0,
        heat_of_adsorption_kj_per_mol_errorbar=0.5,
    )
    return GraspaIsothermResult(
        input_cif=str(input_path),
        output_dir=str(run_dir),
        eqeq_binary="/fake/eqeq",
        graspa_binary="/fake/graspa",
        raspa_backend="graspa",
        eqeq_run_dir=str(run_dir / "eqeq"),
        isotherm_root_dir=str(run_dir / "isotherm"),
        eqeq_input_cif=str(run_dir / "eqeq" / input_path.name),
        eqeq_charged_cif=str(run_dir / "eqeq" / f"{input_path.name}_charged.cif"),
        isotherm_framework_cif=str(run_dir / "isotherm" / "framework.cif"),
        eqeq_json_output_path=None,
        eqeq_stdout_log_path=str(run_dir / "eqeq.stdout.log"),
        eqeq_stderr_log_path=str(run_dir / "eqeq.stderr.log"),
        results_csv_path=str(run_dir / "isotherm" / "results.csv"),
        report_path=str(run_dir / "graspa_isotherm_report.json"),
        unit_cells=(1, 1, 1),
        eqeq_settings=EqeqChargeSettings(),
        isotherm_settings=GraspaIsothermSettings(component=component, pressures=(100000.0,)),
        point_results=(point,),
    )


if __name__ == "__main__":
    unittest.main()
