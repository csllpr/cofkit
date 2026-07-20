from __future__ import annotations

import json
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Literal

from .graspa import (
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
    write_graspa_restart_file,
)
from .guest_forcefields import packaged_guest_forcefield_catalog
from .lammps import LammpsMdResult, LammpsMdSettings, run_lammps_md_on_cif


HybridExchangeMode = Literal["framework", "guest_restart"]


@dataclass(frozen=True)
class HybridMdMcSettings:
    cycles: int = 3
    exchange_mode: HybridExchangeMode = "framework"
    pressure: float = 100_000.0
    components: tuple[GraspaMixtureComponentSettings, ...] = (
        GraspaMixtureComponentSettings(component="CO2_DREIDING", mol_fraction=1.0),
    )
    guest_bundles: tuple[str, ...] = ()
    raspa_backend: str = DEFAULT_RASPA_BACKEND
    raspa_forcefield: str = "dreiding"
    temperature: float = 298.0
    initialization_cycles: int = 50_000
    equilibration_cycles: int = 50_000
    production_cycles: int = 200_000
    number_of_trial_positions: int = 10
    number_of_trial_orientations: int = 10
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
    warnings: tuple[str, ...] = ()

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
            "warnings": list(self.warnings),
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
    run_dir.mkdir(parents=True, exist_ok=True)

    warnings = _hybrid_exchange_warnings(settings)

    current_cif = input_path
    current_guest_restart_state: LammpsGuestRestartState | None = None
    cycle_results: list[HybridMdMcCycleResult] = []
    for cycle in range(1, settings.cycles + 1):
        cycle_dir = run_dir / f"cycle_{cycle:03d}"
        input_guest_restart_state = current_guest_restart_state
        lammps_result = run_lammps_md_on_cif(
            current_cif,
            output_dir=cycle_dir / "1.lammps_md",
            lmp_path=lmp_path,
            eqeq_path=eqeq_path,
            settings=lammps_md_settings,
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
        if len(settings.components) == 1:
            component = settings.components[0]
            gcmc_result_type = "isotherm"
            isotherm_settings = _isotherm_settings_from_hybrid(settings, component)
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
            mixture_settings = _mixture_settings_from_hybrid(settings)
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
            output_guest_restart_state = build_lammps_guest_restart_state_from_gcmc_result(
                gcmc_result,
                components=tuple(component.component for component in settings.components),
                guest_bundles=settings.guest_bundles,
            )
            cycle_warnings.extend(output_guest_restart_state.warnings)
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
        current_guest_restart_state = output_guest_restart_state

    report_path = run_dir / "hybrid_mdmc_report.json"
    result = HybridMdMcResult(
        input_cif=str(input_path),
        output_dir=str(run_dir),
        report_path=str(report_path),
        final_framework_cif=str(current_cif),
        settings=settings,
        lammps_md_settings=lammps_md_settings,
        lammps_eqeq_settings=lammps_eqeq_settings,
        raspa_eqeq_settings=raspa_eqeq_settings,
        cycle_results=tuple(cycle_results),
        warnings=warnings,
    )
    report_path.write_text(json.dumps(result.to_dict(), indent=2), encoding="utf-8")
    return result


def _validate_hybrid_settings(settings: HybridMdMcSettings) -> None:
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
        number_of_trial_positions=settings.number_of_trial_positions,
        number_of_trial_orientations=settings.number_of_trial_orientations,
        cutoff_vdw=settings.cutoff_vdw,
        cutoff_coulomb=settings.cutoff_coulomb,
        ewald_precision=settings.ewald_precision,
    )


__all__ = [
    "HybridMdMcCycleResult",
    "HybridMdMcResult",
    "HybridMdMcSettings",
    "run_hybrid_mdmc_workflow",
]
