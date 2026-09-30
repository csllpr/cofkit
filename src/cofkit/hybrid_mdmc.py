from __future__ import annotations

import json
import sys
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Literal

from .calculation_io import atomic_write_text, begin_attempt, derive_seed, finish_attempt, validate_numbers
from .graspa import (
    _canonical_component_name,
    _guest_nonbonded_convention,
    _load_graspa_guest_bundles,
    _canonicalize_isotherm_settings,
    _canonicalize_mixture_settings,
    _validate_isotherm_settings,
    _validate_mixture_settings,
    DEFAULT_RASPA_BACKEND,
    EqeqChargeSettings,
    GraspaIsothermResult,
    GraspaIsothermSettings,
    GraspaMixtureComponentSettings,
    GraspaMixtureResult,
    GraspaMixtureSettings,
    run_graspa_isotherm_workflow,
    run_graspa_mixture_workflow,
)
from .guest_restart import (
    GuestRestartError,
    LammpsGuestRestartState,
    build_lammps_guest_restart_state_from_gcmc_result,
    build_lammps_guest_restart_state_from_lammps_md_result,
    guest_site_model_diagnostics,
    load_lammps_guest_force_field_assets,
    write_graspa_restart_file,
)
from .guest_forcefields import packaged_guest_forcefield_catalog
from .lammps import LammpsMdResult, LammpsMdSettings, run_lammps_md_on_cif


HybridExchangeMode = Literal["framework", "guest_restart"]


@dataclass(frozen=True)
class HybridMdMcSettings:
    # Heuristic — pending calibration: number of alternating LAMMPS-MD /
    # RASPA-family-GCMC exchange cycles.
    cycles: int = 3
    exchange_mode: HybridExchangeMode = "framework"
    # derived: 100 kPa = 1 bar in the Pascal units RASPA-family inputs use;
    # ambient-pressure adsorption conditions.
    pressure: float = 100_000.0
    components: tuple[GraspaMixtureComponentSettings, ...] = (
        GraspaMixtureComponentSettings(component="CO2_DREIDING", mol_fraction=1.0),
    )
    guest_bundles: tuple[str, ...] = ()
    raspa_backend: str = DEFAULT_RASPA_BACKEND
    raspa_forcefield: str = "dreiding"
    # Mirrors the unified 298 K owner decision documented on
    # GraspaWidomSettings.temperature (graspa.py); per-command owner.
    temperature: float = 298.0
    # Heuristic — pending calibration: GCMC cycle budgets per exchange cycle,
    # following the RASPA2/gRASPA example-input convention (same pattern as
    # the run-control defaults in graspa.py).
    initialization_cycles: int = 50_000
    equilibration_cycles: int = 50_000
    production_cycles: int = 200_000
    # Heuristic — pending calibration: trial positions/orientations per MC
    # move; same convention as the graspa.py settings dataclasses.
    number_of_trial_positions: int = 10
    number_of_trial_orientations: int = 10
    # Deliberately retyped as this dataclass's own defaults (per-command
    # owner): they intentionally mirror graspa.py DEFAULT_CUTOFF_ANGSTROM /
    # DEFAULT_EWALD_PRECISION and must be changed in both places together.
    cutoff_vdw: float = 12.8
    cutoff_coulomb: float = 12.8
    ewald_precision: float = 1.0e-6

    def to_dict(self) -> dict[str, object]:
        return {
            "cycles": self.cycles,
            "exchange_mode": self.exchange_mode,
            "pressure": self.pressure,
            "components": [component.to_dict() for component in self.components],
            "guest_bundles": list(self.guest_bundles),
            "raspa_backend": self.raspa_backend,
            "raspa_forcefield": self.raspa_forcefield,
            "temperature": self.temperature,
            "initialization_cycles": self.initialization_cycles,
            "equilibration_cycles": self.equilibration_cycles,
            "production_cycles": self.production_cycles,
            "number_of_trial_positions": self.number_of_trial_positions,
            "number_of_trial_orientations": self.number_of_trial_orientations,
            "cutoff_vdw": self.cutoff_vdw,
            "cutoff_coulomb": self.cutoff_coulomb,
            "ewald_precision": self.ewald_precision,
        }


@dataclass(frozen=True)
class HybridMdMcCycleResult:
    cycle: int
    input_framework_cif: str
    lammps_md_result: LammpsMdResult
    gcmc_result_type: str
    gcmc_result: GraspaIsothermResult | GraspaMixtureResult
    output_framework_cif: str
    input_guest_restart_source_path: str | None = None
    md_output_guest_restart_source_path: str | None = None
    gcmc_initial_restart_file_path: str | None = None
    output_guest_restart_source_path: str | None = None
    n_input_guest_atoms: int = 0
    n_md_output_guest_atoms: int = 0
    n_output_guest_atoms: int = 0
    input_guest_components: tuple[str, ...] = ()
    md_output_guest_components: tuple[str, ...] = ()
    output_guest_components: tuple[str, ...] = ()
    warnings: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, object]:
        return {
            "cycle": self.cycle,
            "input_framework_cif": self.input_framework_cif,
            "lammps_md_result": self.lammps_md_result.to_dict(),
            "gcmc_result_type": self.gcmc_result_type,
            "gcmc_result": self.gcmc_result.to_dict(),
            "output_framework_cif": self.output_framework_cif,
            "input_guest_restart_source_path": self.input_guest_restart_source_path,
            "md_output_guest_restart_source_path": self.md_output_guest_restart_source_path,
            "gcmc_initial_restart_file_path": self.gcmc_initial_restart_file_path,
            "output_guest_restart_source_path": self.output_guest_restart_source_path,
            "n_input_guest_atoms": self.n_input_guest_atoms,
            "n_md_output_guest_atoms": self.n_md_output_guest_atoms,
            "n_output_guest_atoms": self.n_output_guest_atoms,
            "input_guest_components": list(self.input_guest_components),
            "md_output_guest_components": list(self.md_output_guest_components),
            "output_guest_components": list(self.output_guest_components),
            "warnings": list(self.warnings),
        }


# The hybrid workflow never samples one shared Hamiltonian: fixed-length
# LAMMPS MD segments alternate with grand-canonical gRASPA/RASPA2 segments and
# exchange snapshots. This label is the machine-readable statement of that
# contract in every hybrid report.
HYBRID_MODEL_CONTRACT_KIND = "approximate_alternating_md_mc"

# Run-status vocabulary of the "status" field in hybrid_mdmc_report.json:
# "in_progress" for the partial report (re)written after each completed
# cycle, "completed" when every requested cycle finished, "failed" when a
# cycle raised — in which case `failure` records the cycle, stage, and cause
# and no dependent cycle was run afterwards.
HYBRID_RUN_STATUS_IN_PROGRESS = "in_progress"
HYBRID_RUN_STATUS_COMPLETED = "completed"
HYBRID_RUN_STATUS_FAILED = "failed"


@dataclass(frozen=True)
class HybridMdMcFailure:
    """Where and why a hybrid run stopped.

    `stage` is one of "lammps_md", "md_to_gcmc_restart", "gcmc",
    "gcmc_to_lammps_restart", or "cycle_record"; `error` is the cause as
    ``TypeName: message``.
    """

    cycle: int
    stage: str
    error: str

    def to_dict(self) -> dict[str, object]:
        return {
            "cycle": self.cycle,
            "stage": self.stage,
            "error": self.error,
        }


@dataclass(frozen=True)
class HybridModelContractDifference:
    aspect: str
    md_value: object
    mc_value: object
    note: str

    def to_dict(self) -> dict[str, object]:
        return {
            "aspect": self.aspect,
            "md_value": self.md_value,
            "mc_value": self.mc_value,
            "note": self.note,
        }


@dataclass(frozen=True)
class HybridInteractionSide:
    engine: str
    role: str
    forcefield: str
    vdw_cutoff_angstrom: float
    coulomb_cutoff_angstrom: float
    ewald_precision: float
    temperature_k: float
    vdw_treatment: str
    tail_corrections: bool
    mixing_rule: str
    framework_lj_typing: str
    framework_rigidity: str
    guest_presence: str
    charge_treatment: str

    def to_dict(self) -> dict[str, object]:
        return {
            "engine": self.engine,
            "role": self.role,
            "forcefield": self.forcefield,
            "vdw_cutoff_angstrom": self.vdw_cutoff_angstrom,
            "coulomb_cutoff_angstrom": self.coulomb_cutoff_angstrom,
            "ewald_precision": self.ewald_precision,
            "temperature_k": self.temperature_k,
            "vdw_treatment": self.vdw_treatment,
            "tail_corrections": self.tail_corrections,
            "mixing_rule": self.mixing_rule,
            "framework_lj_typing": self.framework_lj_typing,
            "framework_rigidity": self.framework_rigidity,
            "guest_presence": self.guest_presence,
            "charge_treatment": self.charge_treatment,
        }


@dataclass(frozen=True)
class HybridModelContract:
    """Honest statement of what the alternating MD/MC workflow simulates.

    `differences` lists every interaction/constraint aspect where the LAMMPS
    MD leg and the gRASPA/RASPA2 MC leg do not run equivalent models;
    `guest_model_diagnostics` carries detected guest-model conversions
    (Feynman-Hibbs staged as classical LJ) and conflicting bundle overrides.
    """

    workflow_kind: str
    description: str
    md_side: HybridInteractionSide
    mc_side: HybridInteractionSide
    differences: tuple[HybridModelContractDifference, ...]
    guest_model_diagnostics: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, object]:
        return {
            "workflow_kind": self.workflow_kind,
            "description": self.description,
            "md_side": self.md_side.to_dict(),
            "mc_side": self.mc_side.to_dict(),
            "differences": [difference.to_dict() for difference in self.differences],
            "guest_model_diagnostics": list(self.guest_model_diagnostics),
        }


@dataclass(frozen=True)
class HybridMdMcResult:
    input_cif: str
    output_dir: str
    report_path: str
    final_framework_cif: str
    settings: HybridMdMcSettings
    lammps_md_settings: LammpsMdSettings
    lammps_eqeq_settings: EqeqChargeSettings
    raspa_eqeq_settings: EqeqChargeSettings
    cycle_results: tuple[HybridMdMcCycleResult, ...]
    model_contract: HybridModelContract
    warnings: tuple[str, ...] = ()
    # "in_progress" (partial report between cycles) / "completed" / "failed";
    # see HYBRID_RUN_STATUS_*.
    status: str = HYBRID_RUN_STATUS_COMPLETED
    failure: HybridMdMcFailure | None = None

    def to_dict(self) -> dict[str, object]:
        return {
            "input_cif": self.input_cif,
            "output_dir": self.output_dir,
            "report_path": self.report_path,
            "final_framework_cif": self.final_framework_cif,
            "settings": self.settings.to_dict(),
            "lammps_md_settings": self.lammps_md_settings.to_dict(),
            "lammps_eqeq_settings": self.lammps_eqeq_settings.to_dict(),
            "raspa_eqeq_settings": self.raspa_eqeq_settings.to_dict(),
            "cycle_results": [cycle.to_dict() for cycle in self.cycle_results],
            "model_contract": self.model_contract.to_dict(),
            "warnings": list(self.warnings),
            "status": self.status,
            "failure": self.failure.to_dict() if self.failure is not None else None,
        }


def run_hybrid_mdmc_workflow(
    cif_path: str | Path,
    *,
    output_dir: str | Path | None = None,
    lmp_path: str | Path | None = None,
    eqeq_path: str | Path | None = None,
    graspa_path: str | Path | None = None,
    raspa_path: str | Path | None = None,
    raspa2_path: str | Path | None = None,
    settings: HybridMdMcSettings | None = None,
    lammps_md_settings: LammpsMdSettings | None = None,
    lammps_eqeq_settings: EqeqChargeSettings | None = None,
    raspa_eqeq_settings: EqeqChargeSettings | None = None,
    lammps_timeout_seconds: float = 300.0,
    eqeq_timeout_seconds: float | None = 300.0,
    graspa_timeout_seconds: float | None = None,
) -> HybridMdMcResult:
    settings = settings or HybridMdMcSettings()
    lammps_md_settings = lammps_md_settings or LammpsMdSettings()
    lammps_eqeq_settings = lammps_eqeq_settings or EqeqChargeSettings()
    raspa_eqeq_settings = raspa_eqeq_settings or EqeqChargeSettings()
    _validate_hybrid_settings(settings)

    input_path = Path(cif_path).expanduser().resolve()
    if not input_path.is_file():
        raise FileNotFoundError(f"CIF file does not exist: {input_path}")
    run_dir = (
        input_path.parent / f"{input_path.stem}_hybrid_mdmc"
        if output_dir is None
        else Path(output_dir).expanduser().resolve()
    )
    run_dir = begin_attempt(run_dir, input_path=input_path, settings=settings, binary=None)

    warnings = _hybrid_exchange_warnings(settings)
    model_contract = _build_hybrid_model_contract(
        settings,
        lammps_md_settings,
        lammps_eqeq_settings=lammps_eqeq_settings,
        raspa_eqeq_settings=raspa_eqeq_settings,
    )
    warnings = tuple(
        dict.fromkeys(
            (
                *warnings,
                _model_contract_summary_warning(model_contract),
                *model_contract.guest_model_diagnostics,
            )
        )
    )

    current_cif = input_path
    current_guest_restart_state: LammpsGuestRestartState | None = None
    cycle_results: list[HybridMdMcCycleResult] = []
    report_path = run_dir / "hybrid_mdmc_report.json"

    def _write_report(status: str, failure: HybridMdMcFailure | None) -> HybridMdMcResult:
        report_result = HybridMdMcResult(
            input_cif=str(input_path),
            output_dir=str(run_dir),
            report_path=str(report_path),
            final_framework_cif=str(current_cif),
            settings=settings,
            lammps_md_settings=lammps_md_settings,
            lammps_eqeq_settings=lammps_eqeq_settings,
            raspa_eqeq_settings=raspa_eqeq_settings,
            cycle_results=tuple(cycle_results),
            model_contract=model_contract,
            warnings=warnings,
            status=status,
            failure=failure,
        )
        atomic_write_text(report_path, json.dumps(report_result.to_dict(), indent=2, allow_nan=False))
        return report_result

    # Stage tracking: on a failure the partial report records the failed stage
    # so the cause is attributable (MD leg, MD->MC restart staging, GCMC leg,
    # or MC->MD restart parsing), not just the cycle number.
    cycle = 0
    stage = "setup"
    try:
        for cycle in range(1, settings.cycles + 1):
            cycle_dir = run_dir / f"cycle_{cycle:03d}"
            input_guest_restart_state = current_guest_restart_state
            stage = "lammps_md"
            lammps_result = run_lammps_md_on_cif(
                current_cif,
                output_dir=cycle_dir / "1.lammps_md",
                lmp_path=lmp_path,
                eqeq_path=eqeq_path,
                settings=replace(lammps_md_settings, velocity_seed=derive_seed(lammps_md_settings.velocity_seed, "md", cycle)),
                timeout_seconds=lammps_timeout_seconds,
                eqeq_settings=lammps_eqeq_settings,
                eqeq_timeout_seconds=eqeq_timeout_seconds,
                guest_restart_state=input_guest_restart_state,
            )
            md_framework_cif = Path(lammps_result.output_cif)
            md_output_guest_restart_state: LammpsGuestRestartState | None = None
            gcmc_initial_restart_file_path: str | None = None
            cycle_warnings: list[str] = []
            if settings.exchange_mode == "guest_restart" and input_guest_restart_state is not None:
                stage = "md_to_gcmc_restart"
                try:
                    md_output_guest_restart_state, restart_cell = build_lammps_guest_restart_state_from_lammps_md_result(
                        lammps_result,
                        previous_guest_restart_state=input_guest_restart_state,
                    )
                    restart_file_result = write_graspa_restart_file(
                        md_output_guest_restart_state,
                        cycle_dir / "md_to_gcmc_restartfile",
                        cell=restart_cell,
                        component_order=tuple(component.component for component in settings.components),
                    )
                except GuestRestartError as exc:
                    raise GuestRestartError(
                        f"Failed to convert LAMMPS MD guest coordinates to a gRASPA initial restart file "
                        f"for hybrid cycle {cycle}: {exc}"
                    ) from exc
                gcmc_initial_restart_file_path = restart_file_result.restart_file_path
                cycle_warnings.extend(md_output_guest_restart_state.warnings)
                cycle_warnings.extend(restart_file_result.warnings)
            gcmc_result_type: str
            stage = "gcmc"
            if len(settings.components) == 1:
                component = settings.components[0]
                gcmc_result_type = "isotherm"
                isotherm_settings = replace(_isotherm_settings_from_hybrid(settings, component),
                    random_seed=derive_seed(lammps_md_settings.velocity_seed, "mc", cycle))
                if settings.exchange_mode == "guest_restart":
                    isotherm_settings = replace(
                        isotherm_settings,
                        restart_file=gcmc_initial_restart_file_path is not None,
                    )
                gcmc_result = run_graspa_isotherm_workflow(
                    md_framework_cif,
                    output_dir=cycle_dir / "2.gcmc",
                    eqeq_path=eqeq_path,
                    graspa_path=graspa_path,
                    raspa_path=raspa_path,
                    raspa2_path=raspa2_path,
                    eqeq_settings=raspa_eqeq_settings,
                    isotherm_settings=isotherm_settings,
                    initial_restart_file=gcmc_initial_restart_file_path,
                    eqeq_timeout_seconds=eqeq_timeout_seconds,
                    graspa_timeout_seconds=graspa_timeout_seconds,
                )
            else:
                gcmc_result_type = "mixture"
                mixture_settings = replace(_mixture_settings_from_hybrid(settings),
                    random_seed=derive_seed(lammps_md_settings.velocity_seed, "mc", cycle))
                if settings.exchange_mode == "guest_restart":
                    mixture_settings = replace(
                        mixture_settings,
                        restart_file=gcmc_initial_restart_file_path is not None,
                    )
                gcmc_result = run_graspa_mixture_workflow(
                    md_framework_cif,
                    output_dir=cycle_dir / "2.gcmc",
                    eqeq_path=eqeq_path,
                    graspa_path=graspa_path,
                    raspa_path=raspa_path,
                    raspa2_path=raspa2_path,
                    eqeq_settings=raspa_eqeq_settings,
                    mixture_settings=mixture_settings,
                    initial_restart_file=gcmc_initial_restart_file_path,
                    eqeq_timeout_seconds=eqeq_timeout_seconds,
                    graspa_timeout_seconds=graspa_timeout_seconds,
                )
            output_guest_restart_state: LammpsGuestRestartState | None = None
            if settings.exchange_mode == "guest_restart":
                stage = "gcmc_to_lammps_restart"
                output_guest_restart_state = build_lammps_guest_restart_state_from_gcmc_result(
                    gcmc_result,
                    components=tuple(component.component for component in settings.components),
                    guest_bundles=settings.guest_bundles,
                )
                cycle_warnings.extend(output_guest_restart_state.warnings)
            stage = "cycle_record"
            cycle_results.append(
                HybridMdMcCycleResult(
                    cycle=cycle,
                    input_framework_cif=str(current_cif),
                    lammps_md_result=lammps_result,
                    gcmc_result_type=gcmc_result_type,
                    gcmc_result=gcmc_result,
                    output_framework_cif=str(md_framework_cif),
                    input_guest_restart_source_path=(
                        input_guest_restart_state.source_snapshot_path
                        if input_guest_restart_state is not None
                        else None
                    ),
                    md_output_guest_restart_source_path=(
                        md_output_guest_restart_state.source_snapshot_path
                        if md_output_guest_restart_state is not None
                        else None
                    ),
                    gcmc_initial_restart_file_path=gcmc_initial_restart_file_path,
                    output_guest_restart_source_path=(
                        output_guest_restart_state.source_snapshot_path
                        if output_guest_restart_state is not None
                        else None
                    ),
                    n_input_guest_atoms=input_guest_restart_state.n_atoms if input_guest_restart_state is not None else 0,
                    n_md_output_guest_atoms=(
                        md_output_guest_restart_state.n_atoms if md_output_guest_restart_state is not None else 0
                    ),
                    n_output_guest_atoms=output_guest_restart_state.n_atoms if output_guest_restart_state is not None else 0,
                    input_guest_components=(
                        input_guest_restart_state.components if input_guest_restart_state is not None else ()
                    ),
                    md_output_guest_components=(
                        md_output_guest_restart_state.components if md_output_guest_restart_state is not None else ()
                    ),
                    output_guest_components=(
                        output_guest_restart_state.components if output_guest_restart_state is not None else ()
                    ),
                    warnings=tuple(dict.fromkeys(cycle_warnings)),
                )
            )
            current_cif = md_framework_cif
            current_guest_restart_state = (
                output_guest_restart_state
                if output_guest_restart_state is not None and output_guest_restart_state.n_atoms > 0
                else None
            )
            # Partial report after every completed cycle: a failure in a later
            # cycle must not lose the completed cycles or their warnings.
            _write_report(HYBRID_RUN_STATUS_IN_PROGRESS, None)
    except Exception as exc:
        # Cycles are dependent simulation states: stop the run instead of
        # continuing from an undefined state, and preserve what completed.
        failure = HybridMdMcFailure(
            cycle=cycle,
            stage=stage,
            error=f"{type(exc).__name__}: {exc}",
        )
        try:
            _write_report(HYBRID_RUN_STATUS_FAILED, failure)
        except OSError as report_exc:
            print(
                f"warning: could not write the partial hybrid MD/MC report to {report_path} "
                f"({type(report_exc).__name__}: {report_exc})",
                file=sys.stderr,
            )
        else:
            print(
                f"warning: hybrid MD/MC run failed at cycle {failure.cycle} stage {failure.stage!r} "
                f"({failure.error}); partial report preserving {len(cycle_results)} completed "
                f"cycle(s) written to {report_path}. Dependent cycles were not run.",
                file=sys.stderr,
            )
        raise

    result = _write_report(HYBRID_RUN_STATUS_COMPLETED, None)
    finish_attempt(run_dir)
    return result


def _validate_hybrid_settings(settings: HybridMdMcSettings) -> None:
    validate_numbers(settings)
    if settings.cycles <= 0:
        raise ValueError("cycles must be positive.")
    if settings.exchange_mode not in {"framework", "guest_restart"}:
        raise ValueError("hybrid exchange_mode must be one of: framework, guest_restart.")
    if settings.exchange_mode == "guest_restart" and _normalize_hybrid_raspa_backend(settings.raspa_backend) != "graspa":
        raise ValueError(
            "hybrid exchange_mode='guest_restart' requires raspa_backend='graspa' because post-MD guest "
            "coordinates are staged through gRASPA RestartInitial files; RASPA2 restart staging is not yet supported."
        )
    if settings.exchange_mode == "guest_restart":
        packaged = {name.casefold(): metadata for name, metadata in packaged_guest_forcefield_catalog().items()}
        unsupported = [
            component.component
            for component in settings.components
            if (
                (metadata := packaged.get(component.component.strip().casefold())) is not None
                and not metadata.supports_hybrid_guest_restart
            )
        ]
        if unsupported:
            raise ValueError(
                "hybrid exchange_mode='guest_restart' is not supported for packaged guest model(s): "
                f"{', '.join(unsupported)}. Use exchange_mode='framework' for RASPA-only GCMC coupling."
            )
    if settings.pressure <= 0.0:
        raise ValueError("pressure must be positive.")
    if not settings.components:
        raise ValueError("At least one hybrid MD/MC component must be configured.")
    seen: set[str] = set()
    for component in settings.components:
        if not component.component.strip():
            raise ValueError("component names must not be blank.")
        if component.component in seen:
            raise ValueError(f"Duplicate hybrid MD/MC component {component.component!r} is not allowed.")
        seen.add(component.component)
        if component.mol_fraction <= 0.0:
            raise ValueError("component mol_fraction values must be positive.")
    if settings.temperature <= 0.0:
        raise ValueError("temperature must be positive.")
    if settings.initialization_cycles < 0:
        raise ValueError("initialization_cycles must be non-negative.")
    if settings.equilibration_cycles < 0:
        raise ValueError("equilibration_cycles must be non-negative.")
    if settings.production_cycles <= 0:
        raise ValueError("production_cycles must be positive.")
    if settings.number_of_trial_positions <= 0:
        raise ValueError("number_of_trial_positions must be positive.")
    if settings.number_of_trial_orientations <= 0:
        raise ValueError("number_of_trial_orientations must be positive.")
    if settings.cutoff_vdw <= 0.0:
        raise ValueError("cutoff_vdw must be positive.")
    if settings.cutoff_coulomb <= 0.0:
        raise ValueError("cutoff_coulomb must be positive.")
    if settings.ewald_precision <= 0.0:
        raise ValueError("ewald_precision must be positive.")
    bundles = _load_graspa_guest_bundles(settings.guest_bundles)
    if len(settings.components) == 1:
        adsorption = _canonicalize_isotherm_settings(_isotherm_settings_from_hybrid(settings, settings.components[0]), bundles)
        _validate_isotherm_settings(adsorption, guest_bundles=bundles)
    else:
        adsorption = _canonicalize_mixture_settings(_mixture_settings_from_hybrid(settings), bundles)
        _validate_mixture_settings(adsorption, guest_bundles=bundles)


def _hybrid_exchange_warnings(settings: HybridMdMcSettings) -> tuple[str, ...]:
    if settings.exchange_mode == "framework":
        return (
            "Hybrid exchange_mode='framework' alternates LAMMPS framework MD with gRASPA/RASPA2 GCMC on the "
            "updated framework CIF. Guest molecule coordinates from GCMC are not reinjected into the next "
            "LAMMPS segment in this mode.",
        )
    return (
        "Hybrid exchange_mode='guest_restart' feeds the final GCMC guest restart/movie snapshot into the following "
        "LAMMPS MD segment, then writes post-MD guest coordinates to gRASPA RestartInitial/System_0/restartfile "
        "for the next MC segment. Cycle 1 starts framework-only unless a future API supplies an initial guest restart.",
        "RASPA2 MD-to-MC guest restart staging is not yet supported; guest_restart currently requires the gRASPA backend.",
    )


def _normalize_hybrid_raspa_backend(backend: str) -> str:
    return backend.strip().lower().replace("-", "").replace("_", "")


def _lammps_charge_treatment_description(lammps_md_settings: LammpsMdSettings) -> str:
    if lammps_md_settings.charge_model == "none":
        return "uncharged (charge_model='none')"
    return "EQeq point charges with Ewald kspace (charge_model='eqeq')"


def _build_hybrid_model_contract(
    settings: HybridMdMcSettings,
    lammps_md_settings: LammpsMdSettings,
    *,
    lammps_eqeq_settings: EqeqChargeSettings,
    raspa_eqeq_settings: EqeqChargeSettings,
) -> HybridModelContract:
    bundles = _load_graspa_guest_bundles(settings.guest_bundles)
    component_names = tuple(
        _canonical_component_name(component.component, bundles) for component in settings.components
    )
    mc_vdw_treatment, mc_tail_corrections, mc_mixing_rule = _guest_nonbonded_convention(component_names, bundles)
    guest_restart = settings.exchange_mode == "guest_restart"
    md_charge = _lammps_charge_treatment_description(lammps_md_settings)
    mc_charge = (
        "framework: EQeq point charges from the staged charged CIF (UseChargesFromCIFFile yes); "
        "guests: pseudo-atom model charges from the selected guest models (charge_method Ewald)"
    )
    md_side = HybridInteractionSide(
        engine="lammps",
        role="md",
        forcefield=lammps_md_settings.forcefield,
        vdw_cutoff_angstrom=lammps_md_settings.pair_cutoff,
        coulomb_cutoff_angstrom=lammps_md_settings.coulomb_cutoff,
        ewald_precision=lammps_md_settings.ewald_precision,
        temperature_k=lammps_md_settings.temperature,
        # The LAMMPS staging always emits lj/cut-family pair styles:
        # truncated, unshifted, without tail corrections.
        vdw_treatment="truncated",
        tail_corrections=False,
        # Guest PairIJ rows mix epsilon geometrically and sigma
        # arithmetically, matching the RASPA-side lorentz_berthelot rule.
        mixing_rule="lorentz_berthelot",
        framework_lj_typing=f"full {lammps_md_settings.forcefield.strip().lower()} per-atom-type assignment",
        framework_rigidity="flexible (bonded force field, molecular dynamics)",
        guest_presence=(
            "guests injected as massive sites held by stiff harmonic springs (approximately rigid)"
            if guest_restart
            else "none: guests never enter the MD leg in exchange_mode='framework'"
        ),
        charge_treatment=md_charge,
    )
    mc_side = HybridInteractionSide(
        engine=_normalize_hybrid_raspa_backend(settings.raspa_backend),
        role="mc",
        forcefield=settings.raspa_forcefield,
        vdw_cutoff_angstrom=settings.cutoff_vdw,
        coulomb_cutoff_angstrom=settings.cutoff_coulomb,
        ewald_precision=settings.ewald_precision,
        temperature_k=settings.temperature,
        vdw_treatment=mc_vdw_treatment,
        tail_corrections=mc_tail_corrections,
        mixing_rule=mc_mixing_rule,
        framework_lj_typing=(
            "one representative Lennard-Jones row per element "
            "(element-keyed; not synchronized with the LAMMPS per-type assignment)"
        ),
        framework_rigidity="rigid (framework fixed at the staged CIF)",
        guest_presence="rigid guest molecules sampled grand-canonically",
        charge_treatment=mc_charge,
    )

    differences: list[HybridModelContractDifference] = [
        HybridModelContractDifference(
            aspect="ensemble_and_sampling",
            md_value=f"{lammps_md_settings.ensemble} molecular dynamics (fixed atom count)",
            mc_value="grand-canonical Monte Carlo (fluctuating guest count)",
            note=(
                "The workflow alternates the two legs and exchanges snapshots; it is not a "
                "consistently sampled shared ensemble, and cycle-to-cycle observables are not "
                "equilibration evidence for one joint Hamiltonian."
            ),
        ),
        HybridModelContractDifference(
            aspect="framework_rigidity",
            md_value=md_side.framework_rigidity,
            mc_value=mc_side.framework_rigidity,
            note="The MC leg evaluates guest insertion on a frozen framework while the MD leg moves it.",
        ),
        HybridModelContractDifference(
            aspect="framework_lj_typing",
            md_value=md_side.framework_lj_typing,
            mc_value=mc_side.framework_lj_typing,
            note=(
                "Force-field types of the same element with different Lennard-Jones parameters "
                "collapse to one representative on the MC leg."
            ),
        ),
    ]
    if guest_restart:
        differences.append(
            HybridModelContractDifference(
                aspect="guest_rigidity",
                md_value=md_side.guest_presence,
                mc_value=mc_side.guest_presence,
                note="MD-side guests are held by finite harmonic springs; MC-side guests are exactly rigid.",
            )
        )
    else:
        differences.append(
            HybridModelContractDifference(
                aspect="guests_in_md_leg",
                md_value=md_side.guest_presence,
                mc_value=mc_side.guest_presence,
                note="Guest-framework and guest-guest interactions only act on the MC leg.",
            )
        )
    if lammps_md_settings.forcefield.strip().lower() != settings.raspa_forcefield.strip().lower():
        differences.append(
            HybridModelContractDifference(
                aspect="forcefield_family",
                md_value=lammps_md_settings.forcefield,
                mc_value=settings.raspa_forcefield,
                note="The two legs parameterize the framework with different force-field families.",
            )
        )
    numeric_aspects = (
        ("cutoff_vdw", lammps_md_settings.pair_cutoff, settings.cutoff_vdw, "angstrom"),
        ("cutoff_coulomb", lammps_md_settings.coulomb_cutoff, settings.cutoff_coulomb, "angstrom"),
        ("ewald_precision", lammps_md_settings.ewald_precision, settings.ewald_precision, None),
        ("temperature_k", lammps_md_settings.temperature, settings.temperature, "K"),
    )
    for aspect, md_value, mc_value, unit in numeric_aspects:
        if md_value != mc_value:
            unit_note = f" {unit}" if unit else ""
            differences.append(
                HybridModelContractDifference(
                    aspect=aspect,
                    md_value=md_value,
                    mc_value=mc_value,
                    note=f"The two legs truncate/temper interactions at different{unit_note} settings.",
                )
            )
    if mc_vdw_treatment != md_side.vdw_treatment:
        differences.append(
            HybridModelContractDifference(
                aspect="vdw_treatment",
                md_value=md_side.vdw_treatment,
                mc_value=mc_vdw_treatment,
                note="The MC leg applies a shifted potential while the LAMMPS leg runs truncated unshifted lj/cut.",
            )
        )
    if mc_tail_corrections != md_side.tail_corrections:
        differences.append(
            HybridModelContractDifference(
                aspect="tail_corrections",
                md_value=md_side.tail_corrections,
                mc_value=mc_tail_corrections,
                note="The MC leg applies analytic tail corrections that the LAMMPS leg does not include.",
            )
        )
    if lammps_md_settings.charge_model == "none":
        differences.append(
            HybridModelContractDifference(
                aspect="charge_treatment",
                md_value=md_charge,
                mc_value=mc_charge,
                note="The MC leg is charged while the MD leg runs without electrostatics.",
            )
        )
    else:
        eqeq_differences = sorted(
            key
            for key, value in lammps_eqeq_settings.to_dict().items()
            if raspa_eqeq_settings.to_dict().get(key) != value
        )
        if eqeq_differences:
            differences.append(
                HybridModelContractDifference(
                    aspect="eqeq_settings",
                    md_value={key: lammps_eqeq_settings.to_dict()[key] for key in eqeq_differences},
                    mc_value={key: raspa_eqeq_settings.to_dict()[key] for key in eqeq_differences},
                    note="The two legs equilibrate framework charges with different EQeq settings.",
                )
            )

    guest_model_diagnostics: tuple[str, ...] = ()
    if guest_restart:
        try:
            _templates, guest_sites = load_lammps_guest_force_field_assets(
                component_names, guest_bundles=settings.guest_bundles
            )
        except GuestRestartError:
            # Staging failures are raised with full context at the actual
            # handoff; the contract only reports successful asset loads.
            guest_sites = ()
        guest_model_diagnostics = guest_site_model_diagnostics(guest_sites)

    return HybridModelContract(
        workflow_kind=HYBRID_MODEL_CONTRACT_KIND,
        description=(
            "The hybrid workflow alternates fixed-length LAMMPS MD segments with gRASPA/RASPA2 GCMC "
            "segments and exchanges framework (and optionally guest) snapshots between them. This is an "
            "approximate alternating workflow between two different interaction/constraint models, not a "
            "consistently sampled ensemble of one shared Hamiltonian."
        ),
        md_side=md_side,
        mc_side=mc_side,
        differences=tuple(differences),
        guest_model_diagnostics=guest_model_diagnostics,
    )


def _model_contract_summary_warning(contract: HybridModelContract) -> str:
    aspects = ", ".join(difference.aspect for difference in contract.differences)
    return (
        f"Hybrid model contract: approximate alternating LAMMPS-MD / {contract.mc_side.engine} GCMC "
        "workflow, not a consistently sampled ensemble; "
        f"{len(contract.differences)} interaction/constraint aspect(s) differ between the MD and MC legs "
        f"({aspects}). See model_contract in hybrid_mdmc_report.json."
    )


def _hybrid_guest_snapshot_movies_every(settings: HybridMdMcSettings) -> int | None:
    if settings.exchange_mode != "guest_restart":
        return None
    # derived: gRASPA writes Movies/System_0/result_<cycle>.data during the
    # production phase whenever cycle % MoviesEvery == 0 with 0-based cycle
    # indexing (gRASPA axpy.cu GatherStatisticsDuringSimulation), so
    # MoviesEvery = production_cycles - 1 guarantees a guest snapshot at the
    # final production cycle for the MC->MD handoff, alongside the cycle-0
    # snapshot; clamped to >= 1 for single-cycle production runs.
    return max(1, settings.production_cycles - 1)


def _isotherm_settings_from_hybrid(
    settings: HybridMdMcSettings,
    component: GraspaMixtureComponentSettings,
) -> GraspaIsothermSettings:
    return GraspaIsothermSettings(
        component=component.component,
        guest_bundles=settings.guest_bundles,
        pressures=(settings.pressure,),
        fugacity_coefficient=component.fugacity_coefficient,
        backend=settings.raspa_backend,
        forcefield=settings.raspa_forcefield,
        temperature=settings.temperature,
        initialization_cycles=settings.initialization_cycles,
        equilibration_cycles=settings.equilibration_cycles,
        production_cycles=settings.production_cycles,
        restart_file=settings.exchange_mode == "guest_restart",
        movies_every=_hybrid_guest_snapshot_movies_every(settings),
        number_of_trial_positions=settings.number_of_trial_positions,
        number_of_trial_orientations=settings.number_of_trial_orientations,
        cutoff_vdw=settings.cutoff_vdw,
        cutoff_coulomb=settings.cutoff_coulomb,
        ewald_precision=settings.ewald_precision,
        translation_probability=component.translation_probability,
        rotation_probability=component.rotation_probability,
        reinsertion_probability=component.reinsertion_probability,
        swap_probability=component.swap_probability,
        create_number_of_molecules=component.create_number_of_molecules,
    )


def _mixture_settings_from_hybrid(settings: HybridMdMcSettings) -> GraspaMixtureSettings:
    return GraspaMixtureSettings(
        components=settings.components,
        guest_bundles=settings.guest_bundles,
        pressures=(settings.pressure,),
        backend=settings.raspa_backend,
        forcefield=settings.raspa_forcefield,
        temperature=settings.temperature,
        initialization_cycles=settings.initialization_cycles,
        equilibration_cycles=settings.equilibration_cycles,
        production_cycles=settings.production_cycles,
        restart_file=settings.exchange_mode == "guest_restart",
        movies_every=_hybrid_guest_snapshot_movies_every(settings),
        number_of_trial_positions=settings.number_of_trial_positions,
        number_of_trial_orientations=settings.number_of_trial_orientations,
        cutoff_vdw=settings.cutoff_vdw,
        cutoff_coulomb=settings.cutoff_coulomb,
        ewald_precision=settings.ewald_precision,
    )


__all__ = [
    "HYBRID_MODEL_CONTRACT_KIND",
    "HYBRID_RUN_STATUS_COMPLETED",
    "HYBRID_RUN_STATUS_FAILED",
    "HYBRID_RUN_STATUS_IN_PROGRESS",
    "HybridInteractionSide",
    "HybridMdMcCycleResult",
    "HybridMdMcFailure",
    "HybridMdMcResult",
    "HybridMdMcSettings",
    "HybridModelContract",
    "HybridModelContractDifference",
    "run_hybrid_mdmc_workflow",
]
