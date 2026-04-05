from __future__ import annotations

import argparse
import json

from .graspa import (
    EqeqChargeSettings,
    EqeqExecutionError,
    GraspaConfigurationError,
    GraspaError,
    GraspaExecutionError,
    GraspaParseError,
    GraspaWidomSettings,
    run_graspa_widom_workflow,
)
from .lammps import (
    LammpsConfigurationError,
    LammpsError,
    LammpsExecutionError,
    LammpsInputError,
    LammpsOptimizationSettings,
    LammpsParseError,
    optimize_cif_with_lammps,
)


def add_calculate_group(subparsers) -> None:
    parser = subparsers.add_parser(
        "calculate",
        help="Run external calculation-tool workflows.",
        description="External calculation-tool integrations and optimization workflows.",
    )
    _set_help_default(parser)

    calculate_subparsers = parser.add_subparsers(dest="calculate_command")
    _add_lammps_optimize_parser(calculate_subparsers)
    _add_graspa_widom_parser(calculate_subparsers)


def _set_help_default(parser: argparse.ArgumentParser) -> None:
    def _show_help(_: argparse.Namespace) -> None:
        parser.print_help()

    parser.set_defaults(func=_show_help)


def _add_lammps_optimize_parser(subparsers) -> None:
    parser = subparsers.add_parser(
        "lammps-optimize",
        help="Run a UFF-backed LAMMPS optimization for an explicit-bond CIF.",
        description=(
            "Convert one explicit-bond P1 CIF into a bonded LAMMPS model, populate the LAMMPS data file with "
            "the current UFF backend, run a local minimization workflow, and write an updated CIF with the "
            "same bond topology."
        ),
    )
    parser.add_argument("cif_path", help="Input CIF file to optimize.")
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Directory for the generated LAMMPS inputs, logs, dump, JSON report, and optimized CIF.",
    )
    parser.add_argument(
        "--lmp-path",
        default=None,
        help="Optional explicit path to the LAMMPS executable or alias. Defaults to COFKIT_LMP_PATH or the dev fallback.",
    )
    parser.add_argument(
        "--forcefield",
        default="uff",
        choices=("uff",),
        help=(
            "Forcefield backend to use for the LAMMPS data file. Default: uff. "
            "UFF is the only implemented backend in the current public workflow."
        ),
    )
    parser.add_argument(
        "--pair-cutoff",
        type=float,
        default=12.0,
        help="Global LJ cutoff in angstrom for the selected forcefield backend. Default: 12.0.",
    )
    parser.add_argument(
        "--position-restraint-force-constant",
        type=float,
        default=0.20,
        help="Stage-1 spring/self restraint in kcal/mol/A^2. Default: 0.20.",
    )
    parser.add_argument(
        "--pre-minimization-steps",
        type=int,
        default=10000,
        help=(
            "Restrained pre-minimization MD run length in timesteps. "
            "Default: 10000. Set to 0 to disable the prerun."
        ),
    )
    parser.add_argument(
        "--pre-minimization-temperature",
        type=float,
        default=300.0,
        help="Target temperature in K for the optional prerun Langevin stage. Default: 300.0.",
    )
    parser.add_argument(
        "--pre-minimization-damping",
        type=float,
        default=100.0,
        help="Langevin damping parameter in fs for the optional prerun stage. Default: 100.0.",
    )
    parser.add_argument(
        "--pre-minimization-seed",
        type=int,
        default=246813,
        help="Random seed for the optional prerun velocity/Langevin setup. Default: 246813.",
    )
    parser.add_argument(
        "--pre-minimization-displacement-limit",
        type=float,
        default=0.10,
        help="Per-step displacement cap in angstrom for the optional prerun nve/limit stage. Default: 0.10.",
    )
    parser.add_argument(
        "--two-stage",
        action="store_true",
        help="Run a second minimization stage with a weaker or user-specified restraint.",
    )
    parser.add_argument(
        "--stage2-position-restraint-force-constant",
        type=float,
        default=None,
        help=(
            "Optional stage-2 spring/self restraint in kcal/mol/A^2. "
            "Default when --two-stage is enabled: min(stage1*0.25, 0.05)."
        ),
    )
    parser.add_argument(
        "--energy-tolerance",
        type=float,
        default=1.0e-6,
        help="Stage-1 LAMMPS minimization energy tolerance. Default: 1e-6.",
    )
    parser.add_argument(
        "--force-tolerance",
        type=float,
        default=1.0e-6,
        help="Stage-1 LAMMPS minimization force tolerance. Default: 1e-6.",
    )
    parser.add_argument(
        "--max-iterations",
        type=int,
        default=200000,
        help="Stage-1 maximum LAMMPS minimization iterations. Default: 200000.",
    )
    parser.add_argument(
        "--max-evaluations",
        type=int,
        default=2000000,
        help="Stage-1 maximum LAMMPS minimization force evaluations. Default: 2000000.",
    )
    parser.add_argument(
        "--min-style",
        default="fire",
        help="Stage-1 LAMMPS minimization style. Default: fire.",
    )
    parser.add_argument(
        "--stage2-energy-tolerance",
        type=float,
        default=None,
        help="Optional stage-2 energy tolerance. Defaults to the stage-1 value.",
    )
    parser.add_argument(
        "--stage2-force-tolerance",
        type=float,
        default=None,
        help="Optional stage-2 force tolerance. Defaults to the stage-1 value.",
    )
    parser.add_argument(
        "--stage2-max-iterations",
        type=int,
        default=None,
        help="Optional stage-2 iteration cap. Defaults to the stage-1 value.",
    )
    parser.add_argument(
        "--stage2-max-evaluations",
        type=int,
        default=None,
        help="Optional stage-2 force-evaluation cap. Defaults to the stage-1 value.",
    )
    parser.add_argument(
        "--stage2-min-style",
        default=None,
        help="Optional stage-2 minimization style. Defaults to the stage-1 style.",
    )
    parser.add_argument(
        "--timestep",
        type=float,
        default=None,
        help="Optional global LAMMPS timestep in fs for FIRE or other timestep-aware minimizers.",
    )
    parser.add_argument(
        "--min-modify-dmax",
        type=float,
        default=None,
        help="Optional min_modify dmax value in angstrom.",
    )
    parser.add_argument(
        "--min-modify-line",
        choices=("backtrack", "quadratic", "forcezero", "spin_cubic", "spin_none"),
        default=None,
        help="Optional min_modify line setting.",
    )
    parser.add_argument(
        "--min-modify-norm",
        choices=("two", "inf", "max"),
        default=None,
        help="Optional min_modify norm setting.",
    )
    parser.add_argument(
        "--min-modify-fire-integrator",
        choices=("eulerimplicit", "verlet", "leapfrog", "eulerexplicit"),
        default=None,
        help="Optional FIRE integrator selection for min_modify.",
    )
    parser.add_argument(
        "--min-modify-fire-tmax",
        type=float,
        default=None,
        help="Optional FIRE tmax multiplier for min_modify.",
    )
    parser.add_argument(
        "--min-modify-fire-abcfire",
        action="store_true",
        help="Enable the ABC-FIRE min_modify variant for FIRE stages.",
    )
    parser.add_argument(
        "--relax-cell",
        action="store_true",
        help="Append a final box/relax minimization stage using a box-relax-compatible minimizer.",
    )
    parser.add_argument(
        "--box-relax-mode",
        choices=("auto", "iso", "aniso", "tri"),
        default="auto",
        help="Cell-relax mode for fix box/relax. Default: auto.",
    )
    parser.add_argument(
        "--box-relax-target-pressure",
        type=float,
        default=0.0,
        help="Target pressure for fix box/relax. Default: 0.0.",
    )
    parser.add_argument(
        "--box-relax-vmax",
        type=float,
        default=0.001,
        help="fix box/relax vmax setting. Default: 0.001.",
    )
    parser.add_argument(
        "--box-relax-nreset",
        type=int,
        default=None,
        help="Optional fix box/relax nreset value.",
    )
    parser.add_argument(
        "--box-relax-min-style",
        default="cg",
        help="Minimization style for the final box/relax stage. Default: cg.",
    )
    parser.add_argument(
        "--box-relax-energy-tolerance",
        type=float,
        default=None,
        help="Optional final box-relax energy tolerance. Defaults to the previous stage value.",
    )
    parser.add_argument(
        "--box-relax-force-tolerance",
        type=float,
        default=None,
        help="Optional final box-relax force tolerance. Defaults to the previous stage value.",
    )
    parser.add_argument(
        "--box-relax-max-iterations",
        type=int,
        default=None,
        help="Optional final box-relax iteration cap. Defaults to the previous stage value.",
    )
    parser.add_argument(
        "--box-relax-max-evaluations",
        type=int,
        default=None,
        help="Optional final box-relax force-evaluation cap. Defaults to the previous stage value.",
    )
    parser.add_argument(
        "--timeout-seconds",
        type=float,
        default=300.0,
        help="Per-run timeout for the LAMMPS subprocess. Default: 300.",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Print the full result report as JSON instead of a short human-readable summary.",
    )
    parser.set_defaults(func=_run_lammps_optimize)


def _run_lammps_optimize(args: argparse.Namespace) -> None:
    settings = LammpsOptimizationSettings(
        forcefield=args.forcefield,
        pair_cutoff=args.pair_cutoff,
        position_restraint_force_constant=args.position_restraint_force_constant,
        pre_minimization_steps=args.pre_minimization_steps,
        pre_minimization_temperature=args.pre_minimization_temperature,
        pre_minimization_damping=args.pre_minimization_damping,
        pre_minimization_seed=args.pre_minimization_seed,
        pre_minimization_displacement_limit=args.pre_minimization_displacement_limit,
        two_stage_protocol=args.two_stage,
        stage2_position_restraint_force_constant=args.stage2_position_restraint_force_constant,
        energy_tolerance=args.energy_tolerance,
        force_tolerance=args.force_tolerance,
        max_iterations=args.max_iterations,
        max_evaluations=args.max_evaluations,
        min_style=args.min_style,
        stage2_energy_tolerance=args.stage2_energy_tolerance,
        stage2_force_tolerance=args.stage2_force_tolerance,
        stage2_max_iterations=args.stage2_max_iterations,
        stage2_max_evaluations=args.stage2_max_evaluations,
        stage2_min_style=args.stage2_min_style,
        timestep=args.timestep,
        min_modify_dmax=args.min_modify_dmax,
        min_modify_line=args.min_modify_line,
        min_modify_norm=args.min_modify_norm,
        min_modify_fire_integrator=args.min_modify_fire_integrator,
        min_modify_fire_tmax=args.min_modify_fire_tmax,
        min_modify_fire_abcfire=True if args.min_modify_fire_abcfire else None,
        relax_cell=args.relax_cell,
        box_relax_mode=args.box_relax_mode,
        box_relax_target_pressure=args.box_relax_target_pressure,
        box_relax_vmax=args.box_relax_vmax,
        box_relax_nreset=args.box_relax_nreset,
        box_relax_min_style=args.box_relax_min_style,
        box_relax_energy_tolerance=args.box_relax_energy_tolerance,
        box_relax_force_tolerance=args.box_relax_force_tolerance,
        box_relax_max_iterations=args.box_relax_max_iterations,
        box_relax_max_evaluations=args.box_relax_max_evaluations,
    )
    try:
        result = optimize_cif_with_lammps(
            args.cif_path,
            output_dir=args.output_dir,
            lmp_path=args.lmp_path,
            settings=settings,
            timeout_seconds=args.timeout_seconds,
        )
    except (
        FileNotFoundError,
        ValueError,
        LammpsConfigurationError,
        LammpsInputError,
        LammpsExecutionError,
        LammpsParseError,
        LammpsError,
    ) as exc:
        raise SystemExit(f"error: {exc}") from exc

    if args.json:
        print(json.dumps(result.to_dict(), indent=2))
        return

    print("input_cif:", result.input_cif)
    print("lammps_binary:", result.lammps_binary)
    print("output_dir:", result.output_dir)
    print("optimized_cif:", result.optimized_cif)
    print("n_atoms:", result.n_atoms)
    print("n_bonds:", result.n_bonds)
    print("n_angles:", result.n_angles)
    print("n_atom_types:", result.n_atom_types)
    print("n_bond_types:", result.n_bond_types)
    print("n_angle_types:", result.n_angle_types)
    print("forcefield_backend:", result.forcefield_backend)
    print("warnings:", list(result.warnings))
    print("report_path:", result.report_path)


def _add_graspa_widom_parser(subparsers) -> None:
    parser = subparsers.add_parser(
        "graspa-widom",
        help="Run an EQeq-to-gRASPA Widom insertion workflow for one CIF file.",
        description=(
            "Assign framework charges with EQeq, prepare a gRASPA Widom insertion "
            "run directory, execute gRASPA, and parse Henry coefficients plus "
            "Widom excess chemical potentials into a JSON report."
        ),
    )
    parser.add_argument("cif_path", help="Input CIF file to charge and screen with gRASPA Widom insertion.")
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Directory for EQeq artifacts, the prepared gRASPA run tree, parsed results, and the JSON report.",
    )
    parser.add_argument(
        "--eqeq-path",
        default=None,
        help="Optional explicit path to the EQeq executable. Defaults to COFKIT_EQEQ_PATH.",
    )
    parser.add_argument(
        "--graspa-path",
        default=None,
        help="Optional explicit path to the gRASPA executable. Defaults to COFKIT_GRASPA_PATH.",
    )
    parser.add_argument(
        "--eqeq-lambda",
        type=float,
        default=1.2,
        help="EQeq dielectric screening parameter. Default: 1.2.",
    )
    parser.add_argument(
        "--eqeq-h-i0",
        type=float,
        default=-2.0,
        help="EQeq hydrogen electron affinity parameter. Default: -2.0.",
    )
    parser.add_argument(
        "--eqeq-charge-precision",
        type=int,
        default=3,
        help="Number of digits for EQeq point charges. Default: 3.",
    )
    parser.add_argument(
        "--eqeq-method",
        choices=("ewald", "nonperiodic"),
        default="ewald",
        help="EQeq Coulomb treatment. Default: ewald.",
    )
    parser.add_argument(
        "--eqeq-real-space-cells",
        type=int,
        default=2,
        help="EQeq real-space expansion count. Default: 2.",
    )
    parser.add_argument(
        "--eqeq-reciprocal-space-cells",
        type=int,
        default=2,
        help="EQeq reciprocal-space expansion count. Default: 2.",
    )
    parser.add_argument(
        "--eqeq-eta",
        type=float,
        default=50.0,
        help="EQeq Ewald splitting parameter. Default: 50.0.",
    )
    parser.add_argument(
        "--temperature",
        type=float,
        default=300.0,
        help="gRASPA Widom temperature in K. Default: 300.0.",
    )
    parser.add_argument(
        "--pressure",
        type=float,
        default=100000.0,
        help="gRASPA Widom pressure in Pa. Default: 100000.0.",
    )
    parser.add_argument(
        "--initialization-cycles",
        type=int,
        default=0,
        help="gRASPA NumberOfInitializationCycles. Default: 0.",
    )
    parser.add_argument(
        "--equilibration-cycles",
        type=int,
        default=0,
        help="gRASPA NumberOfEquilibrationCycles. Default: 0.",
    )
    parser.add_argument(
        "--production-cycles",
        type=int,
        default=2000000,
        help="gRASPA NumberOfProductionCycles. Default: 2000000.",
    )
    parser.add_argument(
        "--trial-positions",
        type=int,
        default=10,
        help="gRASPA NumberOfTrialPositions. Default: 10.",
    )
    parser.add_argument(
        "--trial-orientations",
        type=int,
        default=10,
        help="gRASPA NumberOfTrialOrientations. Default: 10.",
    )
    parser.add_argument(
        "--cutoff-vdw",
        type=float,
        default=12.8,
        help="gRASPA CutOffVDW in angstrom. Default: 12.8.",
    )
    parser.add_argument(
        "--cutoff-coulomb",
        type=float,
        default=12.8,
        help="gRASPA CutOffCoulomb in angstrom. Default: 12.8.",
    )
    parser.add_argument(
        "--ewald-precision",
        type=float,
        default=1.0e-6,
        help="gRASPA EwaldPrecision. Default: 1e-6.",
    )
    parser.add_argument(
        "--eqeq-timeout-seconds",
        type=float,
        default=300.0,
        help="Optional timeout for the EQeq subprocess. Default: 300.",
    )
    parser.add_argument(
        "--graspa-timeout-seconds",
        type=float,
        default=None,
        help="Optional timeout for the gRASPA subprocess. Default: no timeout.",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Print the full Widom workflow report as JSON instead of a short human-readable summary.",
    )
    parser.set_defaults(func=_run_graspa_widom)


def _run_graspa_widom(args: argparse.Namespace) -> None:
    eqeq_settings = EqeqChargeSettings(
        lambda_value=args.eqeq_lambda,
        hydrogen_electron_affinity=args.eqeq_h_i0,
        charge_precision=args.eqeq_charge_precision,
        method=args.eqeq_method,
        real_space_cells=args.eqeq_real_space_cells,
        reciprocal_space_cells=args.eqeq_reciprocal_space_cells,
        eta=args.eqeq_eta,
    )
    widom_settings = GraspaWidomSettings(
        temperature=args.temperature,
        pressure=args.pressure,
        initialization_cycles=args.initialization_cycles,
        equilibration_cycles=args.equilibration_cycles,
        production_cycles=args.production_cycles,
        number_of_trial_positions=args.trial_positions,
        number_of_trial_orientations=args.trial_orientations,
        cutoff_vdw=args.cutoff_vdw,
        cutoff_coulomb=args.cutoff_coulomb,
        ewald_precision=args.ewald_precision,
    )
    try:
        result = run_graspa_widom_workflow(
            args.cif_path,
            output_dir=args.output_dir,
            eqeq_path=args.eqeq_path,
            graspa_path=args.graspa_path,
            eqeq_settings=eqeq_settings,
            widom_settings=widom_settings,
            eqeq_timeout_seconds=args.eqeq_timeout_seconds,
            graspa_timeout_seconds=args.graspa_timeout_seconds,
        )
    except (
        FileNotFoundError,
        ValueError,
        EqeqExecutionError,
        GraspaConfigurationError,
        GraspaExecutionError,
        GraspaParseError,
        GraspaError,
    ) as exc:
        raise SystemExit(f"error: {exc}") from exc

    if args.json:
        print(json.dumps(result.to_dict(), indent=2))
        return

    print("input_cif:", result.input_cif)
    print("eqeq_binary:", result.eqeq_binary)
    print("graspa_binary:", result.graspa_binary)
    print("output_dir:", result.output_dir)
    print("eqeq_charged_cif:", result.eqeq_charged_cif)
    print("widom_framework_cif:", result.widom_framework_cif)
    print("unit_cells:", list(result.unit_cells))
    print("components:", [row.component for row in result.component_results])
    print("results_csv_path:", result.results_csv_path)
    print("warnings:", list(result.warnings))
    print("report_path:", result.report_path)
