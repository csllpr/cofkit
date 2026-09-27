"""Wiring assertions for CLI argparse defaults.

Guardrail W5.2 of agent-docs/MAGIC_NUMBER_FIX_PLAN.md: every rewired argparse
default in the analyze/calculate CLI must equal the owning config dataclass
field default. These tests fail if someone retypes a literal in the future.
"""

from cofkit.cli import build_parser
from cofkit.graspa import (
    DEFAULT_CUTOFF_ANGSTROM,
    DEFAULT_EWALD_PRECISION,
    DEFAULT_OVERLAP_CRITERIA,
    DEFAULT_RASPA_BACKEND,
    EqeqChargeSettings,
    GraspaIsothermSettings,
    GraspaMixtureComponentSettings,
    GraspaMixtureSettings,
    GraspaWidomSettings,
)
from cofkit.hybrid_mdmc import HybridMdMcSettings
from cofkit.lammps import LammpsMdSettings, LammpsOptimizationSettings
from cofkit.validation import CoarseValidationThresholds


def _parse(*argv: str):
    return build_parser().parse_args(list(argv))


def _assert_eqeq_defaults(args) -> None:
    defaults = EqeqChargeSettings()
    assert args.eqeq_lambda == defaults.lambda_value
    assert args.eqeq_h_i0 == defaults.hydrogen_electron_affinity
    assert args.eqeq_charge_precision == defaults.charge_precision
    assert args.eqeq_target_charge == defaults.target_charge
    assert args.eqeq_net_charge_tolerance == defaults.net_charge_tolerance
    assert args.eqeq_method == defaults.method
    assert args.eqeq_real_space_cells == defaults.real_space_cells
    assert args.eqeq_reciprocal_space_cells == defaults.reciprocal_space_cells
    assert args.eqeq_eta == defaults.eta


def test_analyze_classify_output_defaults_match_thresholds():
    args = _parse("analyze", "classify-output", "batch_out")
    thresholds = CoarseValidationThresholds()
    assert args.warning_max_bridge_distance_residual == thresholds.warning_max_bridge_distance_residual
    assert args.warning_mean_bridge_distance_residual == thresholds.warning_mean_bridge_distance_residual
    assert args.warning_bad_bridge_distance_residual == thresholds.warning_bad_bridge_distance_residual
    assert args.warning_max_bad_bridge_fraction == thresholds.warning_max_bad_bridge_fraction
    assert args.hard_hard_max_bridge_distance == thresholds.hard_hard_max_bridge_distance
    assert args.hard_max_bridge_distance_residual == thresholds.hard_max_bridge_distance_residual
    assert args.hard_mean_bridge_distance_residual == thresholds.hard_mean_bridge_distance_residual
    assert args.hard_min_bridge_distance_ratio == thresholds.hard_min_bridge_distance_ratio
    assert args.hard_max_bridge_distance_ratio == thresholds.hard_max_bridge_distance_ratio
    assert args.min_nonbonded_heavy_distance == thresholds.min_nonbonded_heavy_distance
    assert args.min_2d_cell_area == thresholds.min_2d_cell_area
    assert args.min_3d_cell_volume == thresholds.min_3d_cell_volume


def test_lammps_optimize_defaults_match_settings():
    args = _parse("calculate", "lammps-optimize", "in.cif")
    settings = LammpsOptimizationSettings()
    assert args.forcefield == settings.forcefield
    assert args.charge_model == settings.charge_model
    assert args.dreiding_hbond == settings.dreiding_hbond
    assert args.pair_cutoff == settings.pair_cutoff
    assert args.coulomb_cutoff == settings.coulomb_cutoff
    assert args.ewald_precision == settings.ewald_precision
    assert args.position_restraint_force_constant == settings.position_restraint_force_constant
    assert args.pre_minimization_mode == settings.pre_minimization_mode
    assert args.pre_minimization_steps == settings.pre_minimization_steps
    assert args.pre_minimization_temperature == settings.pre_minimization_temperature
    assert args.pre_minimization_damping == settings.pre_minimization_damping
    assert args.pre_minimization_seed == settings.pre_minimization_seed
    assert args.pre_minimization_displacement_limit == settings.pre_minimization_displacement_limit
    assert args.soft_pre_minimization_cutoff == settings.soft_pre_minimization_cutoff
    assert args.soft_pre_minimization_coefficients == settings.soft_pre_minimization_coefficients
    assert args.soft_pre_minimization_min_style == settings.soft_pre_minimization_min_style
    assert args.soft_pre_minimization_max_iterations == settings.soft_pre_minimization_max_iterations
    assert args.soft_pre_minimization_max_evaluations == settings.soft_pre_minimization_max_evaluations
    assert args.two_stage == settings.two_stage_protocol
    assert args.dump_interval == settings.dump_interval
    assert args.pressure_tolerance == settings.pressure_tolerance
    assert args.energy_tolerance == settings.energy_tolerance
    assert args.force_tolerance == settings.force_tolerance
    assert args.max_iterations == settings.max_iterations
    assert args.max_evaluations == settings.max_evaluations
    assert args.min_style == settings.min_style
    assert args.relax_cell == settings.relax_cell
    assert args.box_relax_mode == settings.box_relax_mode
    assert args.box_relax_target_pressure == settings.box_relax_target_pressure
    assert args.box_relax_vmax == settings.box_relax_vmax
    assert args.box_relax_min_style == settings.box_relax_min_style
    assert args.box_relax_force_tolerance == settings.box_relax_force_tolerance
    _assert_eqeq_defaults(args)


def test_graspa_widom_defaults_match_settings():
    args = _parse("calculate", "graspa-widom", "in.cif")
    settings = GraspaWidomSettings()
    assert args.backend == DEFAULT_RASPA_BACKEND
    assert args.forcefield == settings.forcefield
    assert args.temperature == settings.temperature
    assert args.pressure == settings.pressure
    assert args.initialization_cycles == settings.initialization_cycles
    assert args.equilibration_cycles == settings.equilibration_cycles
    assert args.trial_positions == settings.number_of_trial_positions
    assert args.trial_orientations == settings.number_of_trial_orientations
    assert args.cutoff_vdw == settings.cutoff_vdw
    assert args.cutoff_coulomb == settings.cutoff_coulomb
    assert args.ewald_precision == settings.ewald_precision
    _assert_eqeq_defaults(args)


def test_graspa_isotherm_defaults_match_settings():
    args = _parse("calculate", "graspa-isotherm", "in.cif", "--component", "CO2_DREIDING")
    settings = GraspaIsothermSettings()
    assert args.backend == DEFAULT_RASPA_BACKEND
    assert args.forcefield == settings.forcefield
    assert args.temperature == settings.temperature
    assert args.fugacity_coefficient == settings.fugacity_coefficient
    assert args.initialization_cycles == settings.initialization_cycles
    assert args.equilibration_cycles == settings.equilibration_cycles
    assert args.production_cycles == settings.production_cycles
    assert args.trial_positions == settings.number_of_trial_positions
    assert args.trial_orientations == settings.number_of_trial_orientations
    assert args.cutoff_vdw == settings.cutoff_vdw
    assert args.cutoff_coulomb == settings.cutoff_coulomb
    assert args.ewald_precision == settings.ewald_precision
    _assert_eqeq_defaults(args)


def test_graspa_mixture_defaults_match_settings():
    args = _parse("calculate", "graspa-mixture", "in.cif")
    mixture_fields = GraspaMixtureSettings.__dataclass_fields__
    component_fields = GraspaMixtureComponentSettings.__dataclass_fields__
    assert args.backend == DEFAULT_RASPA_BACKEND
    assert args.forcefield == mixture_fields["forcefield"].default
    assert args.temperature == mixture_fields["temperature"].default
    assert args.initialization_cycles == mixture_fields["initialization_cycles"].default
    assert args.equilibration_cycles == mixture_fields["equilibration_cycles"].default
    assert args.production_cycles == mixture_fields["production_cycles"].default
    assert args.trial_positions == mixture_fields["number_of_trial_positions"].default
    assert args.trial_orientations == mixture_fields["number_of_trial_orientations"].default
    assert args.translation_probability == component_fields["translation_probability"].default
    assert args.rotation_probability == component_fields["rotation_probability"].default
    assert args.reinsertion_probability == component_fields["reinsertion_probability"].default
    assert args.identity_change_probability == component_fields["identity_change_probability"].default
    assert args.swap_probability == component_fields["swap_probability"].default
    assert args.create_number_of_molecules == component_fields["create_number_of_molecules"].default
    assert args.cutoff_vdw == mixture_fields["cutoff_vdw"].default
    assert args.cutoff_coulomb == mixture_fields["cutoff_coulomb"].default
    assert args.ewald_precision == mixture_fields["ewald_precision"].default
    _assert_eqeq_defaults(args)


def test_hybrid_mdmc_defaults_match_settings():
    args = _parse("calculate", "hybrid-mdmc", "in.cif")
    hybrid = HybridMdMcSettings()
    md = LammpsMdSettings()
    assert args.backend == DEFAULT_RASPA_BACKEND
    assert args.cycles == hybrid.cycles
    assert args.exchange_mode == hybrid.exchange_mode
    assert args.lammps_forcefield == md.forcefield
    assert args.charge_model == md.charge_model
    assert args.dreiding_hbond == md.dreiding_hbond
    assert args.raspa_forcefield == hybrid.raspa_forcefield
    assert args.pressure == hybrid.pressure
    assert args.temperature == hybrid.temperature
    assert args.md_steps == md.steps
    assert args.md_timestep == md.timestep
    assert args.md_ensemble == md.ensemble
    assert args.md_thermostat_damping == md.thermostat_damping
    assert args.md_seed == md.velocity_seed
    assert args.md_dump_interval == md.dump_interval
    assert args.md_position_restraint_force_constant == md.position_restraint_force_constant
    assert args.pair_cutoff == md.pair_cutoff
    assert args.coulomb_cutoff == md.coulomb_cutoff
    assert args.ewald_precision == hybrid.ewald_precision
    assert args.gcmc_initialization_cycles == hybrid.initialization_cycles
    assert args.gcmc_equilibration_cycles == hybrid.equilibration_cycles
    assert args.gcmc_production_cycles == hybrid.production_cycles
    assert args.trial_positions == hybrid.number_of_trial_positions
    assert args.trial_orientations == hybrid.number_of_trial_orientations
    assert args.cutoff_vdw == hybrid.cutoff_vdw
    assert args.cutoff_coulomb == hybrid.cutoff_coulomb
    _assert_eqeq_defaults(args)


def test_widom_and_adsorption_temperature_defaults_unify_at_298():
    # Owner decision 2026-09-28 (MAGIC_NUMBER_FIX_PLAN.md W3.3): the former
    # intentional 300 K Widom / 298 K adsorption split is retired; all
    # calculate workflows default to 298 K, each via its own settings owner.
    widom_args = _parse("calculate", "graspa-widom", "in.cif")
    isotherm_args = _parse("calculate", "graspa-isotherm", "in.cif", "--component", "CO2_DREIDING")
    mixture_args = _parse("calculate", "graspa-mixture", "in.cif")
    hybrid_args = _parse("calculate", "hybrid-mdmc", "in.cif")
    assert widom_args.temperature == GraspaWidomSettings().temperature == 298.0
    assert isotherm_args.temperature == GraspaIsothermSettings().temperature == 298.0
    assert mixture_args.temperature == GraspaMixtureSettings.__dataclass_fields__["temperature"].default == 298.0
    assert hybrid_args.temperature == HybridMdMcSettings().temperature == 298.0


def test_graspa_settings_classes_reference_shared_constants():
    for settings_class in (GraspaWidomSettings, GraspaIsothermSettings, GraspaMixtureSettings):
        fields = settings_class.__dataclass_fields__
        assert fields["cutoff_vdw"].default is DEFAULT_CUTOFF_ANGSTROM
        assert fields["cutoff_coulomb"].default is DEFAULT_CUTOFF_ANGSTROM
        assert fields["overlap_criteria"].default is DEFAULT_OVERLAP_CRITERIA
        assert fields["ewald_precision"].default is DEFAULT_EWALD_PRECISION
