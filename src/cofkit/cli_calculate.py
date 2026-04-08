from __future__ import annotations

import argparse
import json

from .graspa import (
    AVAILABLE_WIDOM_COMPONENTS,
    DEFAULT_WIDOM_MOVES_PER_COMPONENT,
    EqeqChargeSettings,
    EqeqExecutionError,
    GraspaConfigurationError,
    GraspaError,
    GraspaExecutionError,
    GraspaMixtureComponentSettings,
    GraspaMixtureSettings,
    GraspaIsothermSettings,
    GraspaParseError,
    GraspaWidomSettings,
    run_graspa_mixture_workflow,
    run_graspa_isotherm_workflow,
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


_GRASPA_COMPONENT_NAME_MAP = {name.casefold(): name for name in AVAILABLE_WIDOM_COMPONENTS}


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
    _add_graspa_isotherm_parser(calculate_subparsers)
    _add_graspa_mixture_parser(calculate_subparsers)


def _set_help_default(parser: argparse.ArgumentParser) -> None:
    def _show_help(_: argparse.Namespace) -> None:
        parser.print_help()

    parser.set_defaults(func=_show_help)


def _parse_graspa_component_name(raw_value: str) -> str:
    component_name = raw_value.strip()
    if component_name == "":
        raise argparse.ArgumentTypeError("gRASPA component names must not be blank.")
    normalized = _GRASPA_COMPONENT_NAME_MAP.get(component_name.casefold())
    if normalized is None:
        supported = ", ".join(AVAILABLE_WIDOM_COMPONENTS)
        raise argparse.ArgumentTypeError(
            f"Unsupported gRASPA component {raw_value!r}. Supported components: {supported}."
        )
    return normalized


def _parse_graspa_fugacity_coefficient(raw_value: str) -> float | str:
    value = raw_value.strip()
    if value == "":
        raise argparse.ArgumentTypeError("Fugacity coefficient values must not be blank.")
    if value.casefold() == "pr-eos":
        return "PR-EOS"
    try:
        parsed = float(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            "Fugacity coefficient must be a float or the literal 'PR-EOS'."
        ) from exc
    if parsed <= 0.0:
        raise argparse.ArgumentTypeError("Fugacity coefficient must be positive.")
    return parsed


def _parse_graspa_mixture_component_spec(raw_value: str) -> tuple[str, float]:
    value = raw_value.strip()
    if value == "":
        raise argparse.ArgumentTypeError("Mixture component specifications must not be blank.")
    component_name, separator, raw_fraction = value.partition(":")
    if separator == "":
        raise argparse.ArgumentTypeError(
            "Mixture components must use NAME:FRACTION syntax, for example CO2:0.15."
        )
    component = _parse_graspa_component_name(component_name)
    try:
        mol_fraction = float(raw_fraction)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("Mixture component fractions must be numeric.") from exc
    if mol_fraction <= 0.0:
        raise argparse.ArgumentTypeError("Mixture component fractions must be positive.")
    return component, mol_fraction


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
        help="Optional explicit path to the LAMMPS executable or alias. Defaults to COFKIT_LMP_PATH.",
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
        "--charge-model",
        default="eqeq",
        choices=("none", "eqeq"),
        help=(
            "Charge assignment path for the LAMMPS export. "
            "Default: eqeq. Use none only when you explicitly want an uncharged export."
        ),
    )
    parser.add_argument(
        "--eqeq-path",
        default=None,
        help="Optional explicit path to the EQeq executable. Defaults to COFKIT_EQEQ_PATH when charges are enabled.",
    )
    parser.add_argument(
        "--eqeq-lambda",
        type=float,
        default=1.2,
        help="EQeq lambda parameter for the default LAMMPS charge-assignment stage. Default: 1.2.",
    )
    parser.add_argument(
        "--eqeq-h-i0",
        type=float,
        default=-2.0,
        help="EQeq hydrogen electron affinity for the default LAMMPS charge-assignment stage. Default: -2.0.",
    )
    parser.add_argument(
        "--eqeq-charge-precision",
        type=int,
        default=3,
        help="Number of digits for EQeq point charges in the default LAMMPS charge stage. Default: 3.",
    )
    parser.add_argument(
        "--eqeq-method",
        choices=("ewald", "nonperiodic"),
        default="ewald",
        help="EQeq method for the default LAMMPS charge stage. Default: ewald.",
    )
    parser.add_argument(
        "--eqeq-real-space-cells",
        type=int,
        default=2,
        help="EQeq real-space image count for the default LAMMPS charge stage. Default: 2.",
    )
    parser.add_argument(
        "--eqeq-reciprocal-space-cells",
        type=int,
        default=2,
        help="EQeq reciprocal-space image count for the default LAMMPS charge stage. Default: 2.",
    )
    parser.add_argument(
        "--eqeq-eta",
        type=float,
        default=50.0,
        help="EQeq eta parameter for the default LAMMPS charge stage. Default: 50.0.",
    )
    parser.add_argument(
        "--pair-cutoff",
        type=float,
        default=12.0,
        help="Global LJ cutoff in angstrom for the selected forcefield backend. Default: 12.0.",
    )
    parser.add_argument(
        "--coulomb-cutoff",
        type=float,
        default=12.0,
        help="Coulomb cutoff in angstrom when the LAMMPS export includes charges. Default: 12.0.",
    )
    parser.add_argument(
        "--ewald-precision",
        type=float,
        default=1.0e-6,
        help="LAMMPS Ewald precision when the export includes periodic charges. Default: 1e-6.",
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
        help="Target temperature in K for the default prerun Langevin stage. Default: 300.0.",
    )
    parser.add_argument(
        "--pre-minimization-damping",
        type=float,
        default=100.0,
        help="Langevin damping parameter in fs for the default prerun stage. Default: 100.0.",
    )
    parser.add_argument(
        "--pre-minimization-seed",
        type=int,
        default=246813,
        help="Random seed for the default prerun velocity/Langevin setup. Default: 246813.",
    )
    parser.add_argument(
        "--pre-minimization-displacement-limit",
        type=float,
        default=0.10,
        help="Per-step displacement cap in angstrom for the default prerun nve/limit stage. Default: 0.10.",
    )
    parser.add_argument(
        "--two-stage",
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "Enable or disable a second minimization stage. Default: enabled. "
            "When enabled and no explicit stage-2 restraint is provided, stage 2 is unrestrained."
        ),
    )
    parser.add_argument(
        "--stage2-position-restraint-force-constant",
        type=float,
        default=None,
        help=(
            "Optional stage-2 spring/self restraint in kcal/mol/A^2. "
            "Default when stage 2 is enabled: 0.0."
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
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "Enable or disable a final box/relax minimization stage using a box-relax-compatible minimizer. "
            "Default: enabled."
        ),
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
        "--eqeq-timeout-seconds",
        type=float,
        default=300.0,
        help="Per-run timeout for the default EQeq charge-assignment subprocess. Default: 300.",
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
        charge_model=args.charge_model,
        pair_cutoff=args.pair_cutoff,
        coulomb_cutoff=args.coulomb_cutoff,
        ewald_precision=args.ewald_precision,
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
    eqeq_settings = EqeqChargeSettings(
        lambda_value=args.eqeq_lambda,
        hydrogen_electron_affinity=args.eqeq_h_i0,
        charge_precision=args.eqeq_charge_precision,
        method=args.eqeq_method,
        real_space_cells=args.eqeq_real_space_cells,
        reciprocal_space_cells=args.eqeq_reciprocal_space_cells,
        eta=args.eqeq_eta,
    )
    try:
        result = optimize_cif_with_lammps(
            args.cif_path,
            output_dir=args.output_dir,
            lmp_path=args.lmp_path,
            eqeq_path=args.eqeq_path,
            settings=settings,
            timeout_seconds=args.timeout_seconds,
            eqeq_settings=eqeq_settings,
            eqeq_timeout_seconds=args.eqeq_timeout_seconds,
        )
    except (
        FileNotFoundError,
        ValueError,
        GraspaConfigurationError,
        EqeqExecutionError,
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
    print("charge_model:", result.charge_model)
    print("n_charged_atoms:", result.n_charged_atoms)
    print("net_charge:", result.net_charge)
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
        "--component",
        dest="components",
        action="append",
        type=_parse_graspa_component_name,
        default=None,
        metavar="NAME",
        help=(
            "Activate one packaged Widom probe component. Repeat for multiple components. "
            f"Supported: {', '.join(AVAILABLE_WIDOM_COMPONENTS)}."
        ),
    )
    parser.add_argument(
        "--all-components",
        action="store_true",
        help="Activate all packaged Widom probe components.",
    )
    parser.add_argument(
        "--widom-moves-per-component",
        type=int,
        default=DEFAULT_WIDOM_MOVES_PER_COMPONENT,
        help=(
            "Target Widom moves per active component. cofkit multiplies this by the number of "
            f"selected components to derive NumberOfProductionCycles unless --production-cycles is set. "
            f"Default: {DEFAULT_WIDOM_MOVES_PER_COMPONENT}."
        ),
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
        default=None,
        help=(
            "Optional explicit gRASPA NumberOfProductionCycles override. "
            "When omitted, cofkit derives the total from --widom-moves-per-component "
            "times the number of active components."
        ),
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
    if args.production_cycles is None and args.widom_moves_per_component <= 0:
        raise SystemExit("error: --widom-moves-per-component must be positive.")

    selected_components: list[str] = []
    if args.all_components:
        selected_components.extend(AVAILABLE_WIDOM_COMPONENTS)
    if args.components:
        selected_components.extend(args.components)
    deduped_components = tuple(dict.fromkeys(selected_components))
    if not deduped_components:
        raise SystemExit("error: select at least one Widom component with --component or --all-components.")

    if args.production_cycles is None:
        production_cycles = args.widom_moves_per_component * len(deduped_components)
        widom_moves_per_component = args.widom_moves_per_component
    else:
        production_cycles = args.production_cycles
        widom_moves_per_component = None

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
        components=deduped_components,
        temperature=args.temperature,
        pressure=args.pressure,
        initialization_cycles=args.initialization_cycles,
        equilibration_cycles=args.equilibration_cycles,
        production_cycles=production_cycles,
        widom_moves_per_component=widom_moves_per_component,
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


def _add_graspa_isotherm_parser(subparsers) -> None:
    parser = subparsers.add_parser(
        "graspa-isotherm",
        help="Run a single-component EQeq-to-gRASPA GCMC adsorption isotherm workflow for one CIF file.",
        description=(
            "Assign framework charges with EQeq, prepare one gRASPA GCMC adsorption run per pressure point, "
            "execute those runs, and aggregate absolute loading plus heat-of-adsorption summaries into a JSON report."
        ),
    )
    parser.add_argument("cif_path", help="Input CIF file to charge and screen with gRASPA GCMC adsorption.")
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Directory for EQeq artifacts, per-pressure gRASPA runs, parsed results, and the JSON report.",
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
        "--component",
        required=True,
        type=_parse_graspa_component_name,
        metavar="NAME",
        help=f"Select one packaged adsorption component. Supported: {', '.join(AVAILABLE_WIDOM_COMPONENTS)}.",
    )
    parser.add_argument(
        "--pressure",
        dest="pressures",
        action="append",
        type=float,
        default=None,
        metavar="PA",
        help="Add one pressure point in Pa. Repeat for an isotherm sweep.",
    )
    parser.add_argument(
        "--temperature",
        type=float,
        default=298.0,
        help="gRASPA adsorption temperature in K. Default: 298.0.",
    )
    parser.add_argument(
        "--fugacity-coefficient",
        type=_parse_graspa_fugacity_coefficient,
        default=1.0,
        metavar="VALUE",
        help="Adsorbate fugacity coefficient as a positive float or PR-EOS. Default: 1.0.",
    )
    parser.add_argument(
        "--initialization-cycles",
        type=int,
        default=50000,
        help="gRASPA NumberOfInitializationCycles. Default: 50000.",
    )
    parser.add_argument(
        "--equilibration-cycles",
        type=int,
        default=50000,
        help="gRASPA NumberOfEquilibrationCycles. Default: 50000.",
    )
    parser.add_argument(
        "--production-cycles",
        type=int,
        default=200000,
        help="gRASPA NumberOfProductionCycles per pressure point. Default: 200000.",
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
        help="Optional timeout for each gRASPA subprocess. Default: no timeout.",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Print the full adsorption-isotherm workflow report as JSON instead of a short human-readable summary.",
    )
    parser.set_defaults(func=_run_graspa_isotherm)


def _run_graspa_isotherm(args: argparse.Namespace) -> None:
    if not args.pressures:
        raise SystemExit("error: supply at least one pressure point with --pressure.")

    eqeq_settings = EqeqChargeSettings(
        lambda_value=args.eqeq_lambda,
        hydrogen_electron_affinity=args.eqeq_h_i0,
        charge_precision=args.eqeq_charge_precision,
        method=args.eqeq_method,
        real_space_cells=args.eqeq_real_space_cells,
        reciprocal_space_cells=args.eqeq_reciprocal_space_cells,
        eta=args.eqeq_eta,
    )
    isotherm_settings = GraspaIsothermSettings(
        component=args.component,
        pressures=tuple(args.pressures),
        fugacity_coefficient=args.fugacity_coefficient,
        temperature=args.temperature,
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
        result = run_graspa_isotherm_workflow(
            args.cif_path,
            output_dir=args.output_dir,
            eqeq_path=args.eqeq_path,
            graspa_path=args.graspa_path,
            eqeq_settings=eqeq_settings,
            isotherm_settings=isotherm_settings,
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
    print("isotherm_framework_cif:", result.isotherm_framework_cif)
    print("unit_cells:", list(result.unit_cells))
    print("component:", result.isotherm_settings.component)
    print("pressures_pa:", [point.pressure for point in result.point_results])
    print("results_csv_path:", result.results_csv_path)
    print("warnings:", list(result.warnings))
    print("report_path:", result.report_path)


def _add_graspa_mixture_parser(subparsers) -> None:
    parser = subparsers.add_parser(
        "graspa-mixture",
        help="Run an EQeq-to-gRASPA multi-component adsorption/selectivity workflow for one CIF file.",
        description=(
            "Assign framework charges with EQeq, prepare one gRASPA mixture GCMC adsorption run per pressure point, "
            "execute those runs, and aggregate component-resolved loadings plus pairwise selectivities into a JSON report."
        ),
    )
    parser.add_argument("cif_path", help="Input CIF file to screen with gRASPA mixture adsorption.")
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Directory for EQeq artifacts, per-pressure gRASPA runs, parsed results, and the JSON report.",
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
        "--component",
        dest="components",
        action="append",
        type=_parse_graspa_mixture_component_spec,
        default=None,
        metavar="NAME:FRACTION",
        help=(
            "Add one mixture component and its feed mol fraction. Repeat for each component. "
            f"Supported names: {', '.join(AVAILABLE_WIDOM_COMPONENTS)}."
        ),
    )
    parser.add_argument(
        "--pressure",
        dest="pressures",
        action="append",
        type=float,
        default=None,
        metavar="PA",
        help="Add one pressure point in Pa. Repeat for a mixture isotherm sweep.",
    )
    parser.add_argument(
        "--temperature",
        type=float,
        default=298.0,
        help="gRASPA adsorption temperature in K. Default: 298.0.",
    )
    parser.add_argument(
        "--fugacity-coefficient",
        type=_parse_graspa_fugacity_coefficient,
        default=1.0,
        metavar="VALUE",
        help="Component fugacity coefficient as a positive float or PR-EOS. Applied to every component. Default: 1.0.",
    )
    parser.add_argument(
        "--initialization-cycles",
        type=int,
        default=50000,
        help="gRASPA NumberOfInitializationCycles. Default: 50000.",
    )
    parser.add_argument(
        "--equilibration-cycles",
        type=int,
        default=50000,
        help="gRASPA NumberOfEquilibrationCycles. Default: 50000.",
    )
    parser.add_argument(
        "--production-cycles",
        type=int,
        default=200000,
        help="gRASPA NumberOfProductionCycles per pressure point. Default: 200000.",
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
        "--translation-probability",
        type=float,
        default=1.0,
        help="Per-component TranslationProbability. Default: 1.0.",
    )
    parser.add_argument(
        "--rotation-probability",
        type=float,
        default=1.0,
        help="Per-component RotationProbability for rotatable components. Default: 1.0.",
    )
    parser.add_argument(
        "--reinsertion-probability",
        type=float,
        default=1.0,
        help="Per-component ReinsertionProbability. Default: 1.0.",
    )
    parser.add_argument(
        "--identity-change-probability",
        type=float,
        default=1.0,
        help="Per-component IdentityChangeProbability. Default: 1.0.",
    )
    parser.add_argument(
        "--swap-probability",
        type=float,
        default=1.0,
        help="Per-component SwapProbability. Default: 1.0.",
    )
    parser.add_argument(
        "--create-number-of-molecules",
        type=int,
        default=0,
        help="Per-component CreateNumberOfMolecules. Default: 0.",
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
        help="Optional timeout for each gRASPA subprocess. Default: no timeout.",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Print the full mixture adsorption/selectivity report as JSON instead of a short human-readable summary.",
    )
    parser.set_defaults(func=_run_graspa_mixture)


def _run_graspa_mixture(args: argparse.Namespace) -> None:
    if not args.components or len(args.components) < 2:
        raise SystemExit("error: supply at least two mixture components with --component NAME:FRACTION.")
    if not args.pressures:
        raise SystemExit("error: supply at least one pressure point with --pressure.")

    eqeq_settings = EqeqChargeSettings(
        lambda_value=args.eqeq_lambda,
        hydrogen_electron_affinity=args.eqeq_h_i0,
        charge_precision=args.eqeq_charge_precision,
        method=args.eqeq_method,
        real_space_cells=args.eqeq_real_space_cells,
        reciprocal_space_cells=args.eqeq_reciprocal_space_cells,
        eta=args.eqeq_eta,
    )
    mixture_components = tuple(
        GraspaMixtureComponentSettings(
            component=component_name,
            mol_fraction=mol_fraction,
            fugacity_coefficient=args.fugacity_coefficient,
            translation_probability=args.translation_probability,
            rotation_probability=args.rotation_probability,
            reinsertion_probability=args.reinsertion_probability,
            identity_change_probability=args.identity_change_probability,
            swap_probability=args.swap_probability,
            create_number_of_molecules=args.create_number_of_molecules,
        )
        for component_name, mol_fraction in args.components
    )
    mixture_settings = GraspaMixtureSettings(
        components=mixture_components,
        pressures=tuple(args.pressures),
        temperature=args.temperature,
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
        result = run_graspa_mixture_workflow(
            args.cif_path,
            output_dir=args.output_dir,
            eqeq_path=args.eqeq_path,
            graspa_path=args.graspa_path,
            eqeq_settings=eqeq_settings,
            mixture_settings=mixture_settings,
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
    print("mixture_framework_cif:", result.mixture_framework_cif)
    print("unit_cells:", list(result.unit_cells))
    print(
        "components:",
        [
            f"{component.component}:{component.mol_fraction:g}"
            for component in result.mixture_settings.components
        ],
    )
    print("pressures_pa:", [point.pressure for point in result.point_results])
    print("component_results_csv_path:", result.component_results_csv_path)
    print("selectivity_results_csv_path:", result.selectivity_results_csv_path)
    print("warnings:", list(result.warnings))
    print("report_path:", result.report_path)
