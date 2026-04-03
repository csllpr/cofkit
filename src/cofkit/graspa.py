from __future__ import annotations

import csv
import json
import math
import os
import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence


COFKIT_EQEQ_ENV_VAR = "COFKIT_EQEQ_PATH"
COFKIT_GRASPA_ENV_VAR = "COFKIT_GRASPA_PATH"
DEFAULT_EQEQ_BINARY: Path | None = None
DEFAULT_GRASPA_BINARY: Path | None = None
_SUPPORTED_EQEQ_METHODS = {"ewald", "nonperiodic"}
_DEFAULT_WIDOM_COMPONENTS = ("TIP4P", "CO2", "H2", "N2", "SO2")
_NUMBER_RE = re.compile(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?")
_FLOAT_TOKEN_PATTERN = r"(?:[-+]?(?:nan|inf)|[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)"
_WIDOM_RESULT_RE = re.compile(
    r"=+Rosenbluth Summary For Component \[\d+\] \((.+?)\)=+.*?"
    r"Averaged Excess Chemical Potential:\s*"
    rf"({_FLOAT_TOKEN_PATTERN})\s*\+/-\s*"
    rf"({_FLOAT_TOKEN_PATTERN})"
    r".*?Averaged Henry Coefficient \[mol/kg/Pa\]:\s*"
    rf"({_FLOAT_TOKEN_PATTERN})\s*\+/-\s*"
    rf"({_FLOAT_TOKEN_PATTERN})",
    re.DOTALL | re.IGNORECASE,
)


class GraspaError(RuntimeError):
    """Base error for the gRASPA Widom workflow wrapper."""


class GraspaConfigurationError(GraspaError):
    """Raised when external executables or packaged assets cannot be resolved."""


class EqeqExecutionError(GraspaError):
    """Raised when the EQeq subprocess fails or does not produce the charged CIF."""


class GraspaExecutionError(GraspaError):
    """Raised when the gRASPA subprocess exits unsuccessfully."""


class GraspaParseError(GraspaError):
    """Raised when cofkit cannot extract Widom results from gRASPA outputs."""


@dataclass(frozen=True)
class EqeqChargeSettings:
    lambda_value: float = 1.2
    hydrogen_electron_affinity: float = -2.0
    charge_precision: int = 3
    method: str = "ewald"
    real_space_cells: int = 2
    reciprocal_space_cells: int = 2
    eta: float = 50.0

    def to_dict(self) -> dict[str, object]:
        return {
            "lambda_value": self.lambda_value,
            "hydrogen_electron_affinity": self.hydrogen_electron_affinity,
            "charge_precision": self.charge_precision,
            "method": self.method,
            "real_space_cells": self.real_space_cells,
            "reciprocal_space_cells": self.reciprocal_space_cells,
            "eta": self.eta,
        }


@dataclass(frozen=True)
class GraspaWidomSettings:
    components: tuple[str, ...] = _DEFAULT_WIDOM_COMPONENTS
    use_gpu_reduction: bool = False
    use_flag: bool = True
    initialization_cycles: int = 0
    equilibration_cycles: int = 0
    production_cycles: int = 10_000_000
    use_max_step: bool = True
    max_step_per_cycle: int = 1
    use_charges_from_cif_file: bool = True
    restart_file: bool = False
    random_seed: int = 0
    number_of_trial_positions: int = 10
    number_of_trial_orientations: int = 10
    number_of_blocks: int = 1
    adsorbate_allocate_space: int = 10_240
    number_of_simulations: int = 1
    single_simulation: bool = True
    different_frameworks: bool = True
    input_file_type: str = "cif"
    framework_name: str = "framework"
    charge_method: str = "Ewald"
    temperature: float = 300.0
    pressure: float = 100_000.0
    overlap_criteria: float = 1.0e5
    cutoff_vdw: float = 12.8
    cutoff_coulomb: float = 12.8
    ewald_precision: float = 1.0e-6
    use_dnn_for_host_guest: bool = False
    save_output_to_file: bool = True

    def to_dict(self) -> dict[str, object]:
        return {
            "components": list(self.components),
            "use_gpu_reduction": self.use_gpu_reduction,
            "use_flag": self.use_flag,
            "initialization_cycles": self.initialization_cycles,
            "equilibration_cycles": self.equilibration_cycles,
            "production_cycles": self.production_cycles,
            "use_max_step": self.use_max_step,
            "max_step_per_cycle": self.max_step_per_cycle,
            "use_charges_from_cif_file": self.use_charges_from_cif_file,
            "restart_file": self.restart_file,
            "random_seed": self.random_seed,
            "number_of_trial_positions": self.number_of_trial_positions,
            "number_of_trial_orientations": self.number_of_trial_orientations,
            "number_of_blocks": self.number_of_blocks,
            "adsorbate_allocate_space": self.adsorbate_allocate_space,
            "number_of_simulations": self.number_of_simulations,
            "single_simulation": self.single_simulation,
            "different_frameworks": self.different_frameworks,
            "input_file_type": self.input_file_type,
            "framework_name": self.framework_name,
            "charge_method": self.charge_method,
            "temperature": self.temperature,
            "pressure": self.pressure,
            "overlap_criteria": self.overlap_criteria,
            "cutoff_vdw": self.cutoff_vdw,
            "cutoff_coulomb": self.cutoff_coulomb,
            "ewald_precision": self.ewald_precision,
            "use_dnn_for_host_guest": self.use_dnn_for_host_guest,
            "save_output_to_file": self.save_output_to_file,
        }


@dataclass(frozen=True)
class GraspaWidomComponentResult:
    component: str
    widom_energy: float
    widom_energy_errorbar: float
    henry: float
    henry_errorbar: float
    source_data_file: str

    def to_dict(self) -> dict[str, object]:
        return {
            "component": self.component,
            "widom_energy": _json_safe_float(self.widom_energy),
            "widom_energy_errorbar": _json_safe_float(self.widom_energy_errorbar),
            "henry": _json_safe_float(self.henry),
            "henry_errorbar": _json_safe_float(self.henry_errorbar),
            "source_data_file": self.source_data_file,
        }


@dataclass(frozen=True)
class GraspaWidomResult:
    input_cif: str
    output_dir: str
    eqeq_binary: str
    graspa_binary: str
    eqeq_run_dir: str
    widom_run_dir: str
    eqeq_input_cif: str
    eqeq_charged_cif: str
    widom_framework_cif: str
    eqeq_json_output_path: str | None
    eqeq_stdout_log_path: str
    eqeq_stderr_log_path: str
    simulation_input_path: str
    graspa_stdout_log_path: str
    graspa_stderr_log_path: str
    widom_output_dir: str
    data_file_paths: tuple[str, ...]
    results_csv_path: str
    report_path: str
    unit_cells: tuple[int, int, int]
    eqeq_settings: EqeqChargeSettings
    widom_settings: GraspaWidomSettings
    component_results: tuple[GraspaWidomComponentResult, ...]
    warnings: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, object]:
        return {
            "input_cif": self.input_cif,
            "output_dir": self.output_dir,
            "eqeq_binary": self.eqeq_binary,
            "graspa_binary": self.graspa_binary,
            "eqeq_run_dir": self.eqeq_run_dir,
            "widom_run_dir": self.widom_run_dir,
            "eqeq_input_cif": self.eqeq_input_cif,
            "eqeq_charged_cif": self.eqeq_charged_cif,
            "widom_framework_cif": self.widom_framework_cif,
            "eqeq_json_output_path": self.eqeq_json_output_path,
            "eqeq_stdout_log_path": self.eqeq_stdout_log_path,
            "eqeq_stderr_log_path": self.eqeq_stderr_log_path,
            "simulation_input_path": self.simulation_input_path,
            "graspa_stdout_log_path": self.graspa_stdout_log_path,
            "graspa_stderr_log_path": self.graspa_stderr_log_path,
            "widom_output_dir": self.widom_output_dir,
            "data_file_paths": list(self.data_file_paths),
            "results_csv_path": self.results_csv_path,
            "report_path": self.report_path,
            "unit_cells": list(self.unit_cells),
            "eqeq_settings": self.eqeq_settings.to_dict(),
            "widom_settings": self.widom_settings.to_dict(),
            "component_results": [component.to_dict() for component in self.component_results],
            "warnings": list(self.warnings),
        }


def resolve_eqeq_binary(eqeq_path: str | Path | None = None) -> Path:
    return _resolve_binary(
        explicit_path=eqeq_path,
        env_var=COFKIT_EQEQ_ENV_VAR,
        default_path=DEFAULT_EQEQ_BINARY,
        display_name="EQeq",
    )


def resolve_graspa_binary(graspa_path: str | Path | None = None) -> Path:
    return _resolve_binary(
        explicit_path=graspa_path,
        env_var=COFKIT_GRASPA_ENV_VAR,
        default_path=DEFAULT_GRASPA_BINARY,
        display_name="gRASPA",
    )


def run_graspa_widom_workflow(
    cif_path: str | Path,
    *,
    output_dir: str | Path | None = None,
    eqeq_path: str | Path | None = None,
    graspa_path: str | Path | None = None,
    eqeq_settings: EqeqChargeSettings | None = None,
    widom_settings: GraspaWidomSettings | None = None,
    eqeq_timeout_seconds: float | None = 300.0,
    graspa_timeout_seconds: float | None = None,
) -> GraspaWidomResult:
    input_path = Path(cif_path).expanduser().resolve()
    if not input_path.is_file():
        raise FileNotFoundError(f"CIF file does not exist: {input_path}")

    eqeq_settings = eqeq_settings or EqeqChargeSettings()
    widom_settings = widom_settings or GraspaWidomSettings()
    _validate_eqeq_settings(eqeq_settings)
    _validate_widom_settings(widom_settings)
    eqeq_timeout_seconds = _normalize_timeout(eqeq_timeout_seconds, "eqeq_timeout_seconds")
    graspa_timeout_seconds = _normalize_timeout(graspa_timeout_seconds, "graspa_timeout_seconds")

    eqeq_binary = resolve_eqeq_binary(eqeq_path)
    graspa_binary = resolve_graspa_binary(graspa_path)

    run_dir = _resolve_output_dir(input_path, output_dir)
    eqeq_run_dir = run_dir / "eqeq"
    widom_run_dir = run_dir / "widom"
    eqeq_run_dir.mkdir(parents=True, exist_ok=True)
    widom_run_dir.mkdir(parents=True, exist_ok=True)

    eqeq_input_path = eqeq_run_dir / input_path.name
    shutil.copy2(input_path, eqeq_input_path)
    eqeq_stdout_log_path = eqeq_run_dir / "eqeq.stdout.log"
    eqeq_stderr_log_path = eqeq_run_dir / "eqeq.stderr.log"

    eqeq_command = [
        str(eqeq_binary),
        eqeq_input_path.name,
        _format_cli_number(eqeq_settings.lambda_value),
        _format_cli_number(eqeq_settings.hydrogen_electron_affinity),
        str(eqeq_settings.charge_precision),
        eqeq_settings.method,
        str(eqeq_settings.real_space_cells),
        str(eqeq_settings.reciprocal_space_cells),
        _format_cli_number(eqeq_settings.eta),
    ]
    try:
        with eqeq_stdout_log_path.open("w", encoding="utf-8") as stdout_handle:
            with eqeq_stderr_log_path.open("w", encoding="utf-8") as stderr_handle:
                eqeq_completed = subprocess.run(
                    eqeq_command,
                    cwd=eqeq_run_dir,
                    check=False,
                    stdout=stdout_handle,
                    stderr=stderr_handle,
                    text=True,
                    timeout=eqeq_timeout_seconds,
                )
    except subprocess.TimeoutExpired as exc:
        raise EqeqExecutionError(
            f"EQeq timed out after {eqeq_timeout_seconds} seconds. "
            f"See {eqeq_stdout_log_path} and {eqeq_stderr_log_path}."
        ) from exc

    eqeq_output_stem = _expected_eqeq_output_stem(eqeq_input_path.name, eqeq_settings)
    eqeq_charged_cif_path = eqeq_run_dir / f"{eqeq_output_stem}.cif"
    eqeq_json_output_path = eqeq_run_dir / f"{eqeq_output_stem}.json"

    if eqeq_completed.returncode != 0:
        raise EqeqExecutionError(
            f"EQeq failed with exit code {eqeq_completed.returncode}. "
            f"See {eqeq_stdout_log_path} and {eqeq_stderr_log_path}."
        )
    if not eqeq_charged_cif_path.is_file():
        raise EqeqExecutionError(
            f"EQeq completed without writing the expected charged CIF: {eqeq_charged_cif_path}"
        )

    widom_framework_cif_path = widom_run_dir / f"{widom_settings.framework_name}.cif"
    shutil.copy2(eqeq_charged_cif_path, widom_framework_cif_path)
    _copy_widom_template_assets(widom_run_dir, widom_settings.components)

    unit_cells = _compute_unit_cells_from_cif(
        widom_framework_cif_path,
        cutoff=max(widom_settings.cutoff_vdw, widom_settings.cutoff_coulomb),
    )
    simulation_input_path = widom_run_dir / "simulation.input"
    simulation_input_path.write_text(
        _render_widom_simulation_input(widom_settings, unit_cells=unit_cells),
        encoding="utf-8",
    )

    graspa_stdout_log_path = widom_run_dir / "graspa.stdout.log"
    graspa_stderr_log_path = widom_run_dir / "graspa.stderr.log"
    try:
        with graspa_stdout_log_path.open("w", encoding="utf-8") as stdout_handle:
            with graspa_stderr_log_path.open("w", encoding="utf-8") as stderr_handle:
                graspa_completed = subprocess.run(
                    [str(graspa_binary)],
                    cwd=widom_run_dir,
                    check=False,
                    stdout=stdout_handle,
                    stderr=stderr_handle,
                    text=True,
                    timeout=graspa_timeout_seconds,
                )
    except subprocess.TimeoutExpired as exc:
        raise GraspaExecutionError(
            f"gRASPA timed out after {graspa_timeout_seconds} seconds. "
            f"See {graspa_stdout_log_path} and {graspa_stderr_log_path}."
        ) from exc

    if graspa_completed.returncode != 0:
        raise GraspaExecutionError(
            f"gRASPA failed with exit code {graspa_completed.returncode}. "
            f"See {graspa_stdout_log_path} and {graspa_stderr_log_path}."
        )

    widom_output_dir = widom_run_dir / "Output"
    data_file_paths = tuple(str(path) for path in sorted(widom_output_dir.glob("*.data")))
    if not data_file_paths:
        raise GraspaParseError(
            f"gRASPA completed without producing any Output/*.data files under {widom_output_dir}"
        )

    component_results = _parse_widom_result_files(tuple(Path(path) for path in data_file_paths))
    if not component_results:
        raise GraspaParseError(
            f"No Widom component summaries could be parsed from {len(data_file_paths)} data file(s)."
        )

    results_csv_path = widom_output_dir / "results.csv"
    _write_results_csv(component_results, results_csv_path)

    warnings: list[str] = []
    if not eqeq_json_output_path.is_file():
        warnings.append("EQeq did not write the companion JSON output file.")
    if len(data_file_paths) > 1:
        warnings.append(f"Parsed Widom results from {len(data_file_paths)} data files.")
    if graspa_stderr_log_path.read_text(encoding="utf-8", errors="replace").strip():
        warnings.append("gRASPA wrote content to stderr; inspect the stderr log if needed.")
    if any(
        not math.isfinite(value)
        for component in component_results
        for value in (
            component.widom_energy,
            component.widom_energy_errorbar,
            component.henry,
            component.henry_errorbar,
        )
    ):
        warnings.append("One or more Widom summary values are non-finite; inspect the raw gRASPA data output.")

    report_path = run_dir / "graspa_widom_report.json"
    result = GraspaWidomResult(
        input_cif=str(input_path),
        output_dir=str(run_dir),
        eqeq_binary=str(eqeq_binary),
        graspa_binary=str(graspa_binary),
        eqeq_run_dir=str(eqeq_run_dir),
        widom_run_dir=str(widom_run_dir),
        eqeq_input_cif=str(eqeq_input_path),
        eqeq_charged_cif=str(eqeq_charged_cif_path),
        widom_framework_cif=str(widom_framework_cif_path),
        eqeq_json_output_path=str(eqeq_json_output_path) if eqeq_json_output_path.is_file() else None,
        eqeq_stdout_log_path=str(eqeq_stdout_log_path),
        eqeq_stderr_log_path=str(eqeq_stderr_log_path),
        simulation_input_path=str(simulation_input_path),
        graspa_stdout_log_path=str(graspa_stdout_log_path),
        graspa_stderr_log_path=str(graspa_stderr_log_path),
        widom_output_dir=str(widom_output_dir),
        data_file_paths=data_file_paths,
        results_csv_path=str(results_csv_path),
        report_path=str(report_path),
        unit_cells=unit_cells,
        eqeq_settings=eqeq_settings,
        widom_settings=widom_settings,
        component_results=component_results,
        warnings=tuple(warnings),
    )
    report_path.write_text(json.dumps(result.to_dict(), indent=2, allow_nan=False), encoding="utf-8")
    return result


def _resolve_binary(
    *,
    explicit_path: str | Path | None,
    env_var: str,
    default_path: Path | None,
    display_name: str,
) -> Path:
    raw_value = str(explicit_path) if explicit_path is not None else os.environ.get(env_var)
    if raw_value:
        candidate = Path(raw_value).expanduser()
        if not candidate.is_file():
            resolved = shutil.which(raw_value)
            if resolved:
                candidate = Path(resolved)
    elif default_path is not None and default_path.is_file():
        candidate = default_path
    else:
        raise GraspaConfigurationError(
            f"{display_name} binary is not configured. Set {env_var} or provide an explicit path."
        )

    if not candidate.is_file():
        raise GraspaConfigurationError(f"{display_name} binary does not exist: {candidate}")
    if not os.access(candidate, os.X_OK):
        raise GraspaConfigurationError(f"{display_name} binary is not executable: {candidate}")
    return candidate.resolve()


def _normalize_timeout(value: float | None, field_name: str) -> float | None:
    if value is None:
        return None
    if value <= 0.0:
        raise ValueError(f"{field_name} must be positive when provided.")
    return float(value)


def _validate_eqeq_settings(settings: EqeqChargeSettings) -> None:
    if settings.method not in _SUPPORTED_EQEQ_METHODS:
        raise ValueError(
            f"Unsupported EQeq method {settings.method!r}. "
            f"Supported methods: {sorted(_SUPPORTED_EQEQ_METHODS)}."
        )
    if settings.charge_precision < 0:
        raise ValueError("charge_precision must be non-negative.")
    if settings.real_space_cells < 0 or settings.reciprocal_space_cells < 0:
        raise ValueError("real_space_cells and reciprocal_space_cells must be non-negative.")
    if settings.eta <= 0.0:
        raise ValueError("eta must be positive.")


def _validate_widom_settings(settings: GraspaWidomSettings) -> None:
    if not settings.components:
        raise ValueError("At least one Widom component must be configured.")
    if settings.framework_name.strip() == "":
        raise ValueError("framework_name must not be blank.")
    if settings.input_file_type.lower() != "cif":
        raise ValueError("Only CIF gRASPA framework inputs are supported in the current pipeline.")
    if settings.initialization_cycles < 0:
        raise ValueError("initialization_cycles must be non-negative.")
    if settings.equilibration_cycles < 0:
        raise ValueError("equilibration_cycles must be non-negative.")
    if settings.production_cycles <= 0:
        raise ValueError("production_cycles must be positive.")
    if settings.number_of_trial_positions <= 0 or settings.number_of_trial_orientations <= 0:
        raise ValueError("number_of_trial_positions and number_of_trial_orientations must be positive.")
    if settings.number_of_blocks <= 0:
        raise ValueError("number_of_blocks must be positive.")
    if settings.adsorbate_allocate_space <= 0:
        raise ValueError("adsorbate_allocate_space must be positive.")
    if settings.number_of_simulations <= 0:
        raise ValueError("number_of_simulations must be positive.")
    if settings.cutoff_vdw <= 0.0 or settings.cutoff_coulomb <= 0.0:
        raise ValueError("cutoff_vdw and cutoff_coulomb must be positive.")
    if settings.ewald_precision <= 0.0:
        raise ValueError("ewald_precision must be positive.")
    if settings.temperature <= 0.0:
        raise ValueError("temperature must be positive.")
    if settings.pressure < 0.0:
        raise ValueError("pressure must be non-negative.")
    if not settings.save_output_to_file:
        raise ValueError("save_output_to_file must remain enabled for result parsing in the current workflow.")


def _resolve_output_dir(input_path: Path, output_dir: str | Path | None) -> Path:
    if output_dir is not None:
        return Path(output_dir).expanduser().resolve()
    return input_path.parent / f"{input_path.stem}_graspa_widom"


def _expected_eqeq_output_stem(input_filename: str, settings: EqeqChargeSettings) -> str:
    return (
        f"{input_filename}_EQeq_{settings.method}_"
        f"{settings.lambda_value:.2f}_{settings.hydrogen_electron_affinity:.2f}"
    )


def _format_cli_number(value: float) -> str:
    return f"{value:g}"


def _widom_template_dir() -> Path:
    template_dir = Path(__file__).resolve().parent / "data" / "graspa" / "widom_template"
    if not template_dir.is_dir():
        raise GraspaConfigurationError(f"Packaged Widom template directory is missing: {template_dir}")
    return template_dir


def _copy_widom_template_assets(destination_dir: Path, components: Sequence[str]) -> None:
    template_dir = _widom_template_dir()
    common_files = ("force_field.def", "force_field_mixing_rules.def", "pseudo_atoms.def")
    required_names = [*common_files, *(f"{name}.def" for name in components)]
    for filename in required_names:
        source_path = template_dir / filename
        if not source_path.is_file():
            raise GraspaConfigurationError(f"Packaged Widom asset is missing: {source_path}")
        shutil.copy2(source_path, destination_dir / filename)


def _compute_unit_cells_from_cif(cif_path: Path, *, cutoff: float) -> tuple[int, int, int]:
    lengths = _extract_cell_lengths(cif_path)
    missing = [key for key in ("_cell_length_a", "_cell_length_b", "_cell_length_c") if key not in lengths]
    if missing:
        raise GraspaParseError(
            f"Missing cell lengths {', '.join(missing)} in CIF file required for Widom UnitCells calculation: {cif_path}"
        )
    return (
        _compute_supercell(lengths["_cell_length_a"], cutoff),
        _compute_supercell(lengths["_cell_length_b"], cutoff),
        _compute_supercell(lengths["_cell_length_c"], cutoff),
    )


def _extract_cell_lengths(cif_path: Path) -> dict[str, float]:
    lengths: dict[str, float] = {}
    with cif_path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            stripped = raw_line.strip()
            if not stripped:
                continue
            lower = stripped.lower()
            for tag in ("_cell_length_a", "_cell_length_b", "_cell_length_c"):
                if lower.startswith(tag):
                    parts = stripped.split()
                    if len(parts) >= 2:
                        value = _parse_numeric_value(parts[1])
                        if value is not None:
                            lengths[tag] = value
            if len(lengths) == 3:
                break
    return lengths


def _parse_numeric_value(raw_value: str) -> float | None:
    cleaned = re.sub(r"\([^)]*\)", "", raw_value)
    match = _NUMBER_RE.search(cleaned)
    if not match:
        return None
    try:
        return float(match.group(0))
    except ValueError:
        return None


def _compute_supercell(cell_length: float, cutoff: float) -> int:
    if cell_length <= 0.0:
        return 1
    return max(int((2.0 * cutoff + cell_length - 1e-9) / cell_length), 1)


def _render_widom_simulation_input(
    settings: GraspaWidomSettings,
    *,
    unit_cells: tuple[int, int, int],
) -> str:
    lines = [
        f"UseGPUReduction {_bool_token(settings.use_gpu_reduction)}",
        f"Useflag {_bool_token(settings.use_flag)}",
        "",
        f"NumberOfInitializationCycles {settings.initialization_cycles}",
        f"NumberOfEquilibrationCycles  {settings.equilibration_cycles}",
        f"NumberOfProductionCycles     {settings.production_cycles}",
        "",
        f"UseMaxStep  {_bool_token(settings.use_max_step)}",
        f"MaxStepPerCycle {settings.max_step_per_cycle}",
        "",
        f"UseChargesFromCIFFile {_bool_token(settings.use_charges_from_cif_file)}",
        "",
        f"RestartFile {_bool_token(settings.restart_file)}",
        f"RandomSeed  {settings.random_seed}",
        "",
        f"NumberOfTrialPositions {settings.number_of_trial_positions}",
        f"NumberOfTrialOrientations {settings.number_of_trial_orientations}",
        "",
        f"NumberOfBlocks {settings.number_of_blocks}",
        f"AdsorbateAllocateSpace {settings.adsorbate_allocate_space}",
        f"NumberOfSimulations {settings.number_of_simulations}",
        f"SingleSimulation {_bool_token(settings.single_simulation)}",
        "",
        f"DifferentFrameworks {_bool_token(settings.different_frameworks)}",
        f"InputFileType {settings.input_file_type}",
        f"FrameworkName {settings.framework_name}",
        f"UnitCells 0 {unit_cells[0]} {unit_cells[1]} {unit_cells[2]}",
        "",
        f"ChargeMethod {settings.charge_method}",
        f"Temperature  {_format_simulation_number(settings.temperature)}",
        f"Pressure     {_format_simulation_number(settings.pressure)}",
        "",
        f"OverlapCriteria {_format_simulation_number(settings.overlap_criteria)}",
        f"CutOffVDW {_format_simulation_number(settings.cutoff_vdw)}",
        f"CutOffCoulomb {_format_simulation_number(settings.cutoff_coulomb)}",
        f"EwaldPrecision {_format_simulation_number(settings.ewald_precision)}",
        "",
        f"UseDNNforHostGuest {_bool_token(settings.use_dnn_for_host_guest)}",
        "",
        f"SaveOutputToFile {_bool_token(settings.save_output_to_file)}",
        "",
    ]
    for index, component in enumerate(settings.components):
        lines.extend(
            [
                f"Component {index} MoleculeName              {component}",
                "            IdealGasRosenbluthWeight 1.0",
                "            FugacityCoefficient      1.0",
                "            WidomProbability         1.0",
                "            CreateNumberOfMolecules  0",
                "",
            ]
        )
    return "\n".join(lines).rstrip() + "\n"


def _bool_token(value: bool) -> str:
    return "yes" if value else "no"


def _format_simulation_number(value: float) -> str:
    return f"{value:g}"


def _json_safe_float(value: float) -> float | None:
    if math.isfinite(value):
        return value
    return None


def _parse_widom_result_files(data_paths: Sequence[Path]) -> tuple[GraspaWidomComponentResult, ...]:
    results: list[GraspaWidomComponentResult] = []
    for path in data_paths:
        content = path.read_text(encoding="utf-8", errors="replace")
        for match in _WIDOM_RESULT_RE.finditer(content):
            results.append(
                GraspaWidomComponentResult(
                    component=match.group(1),
                    widom_energy=float(match.group(2)),
                    widom_energy_errorbar=float(match.group(3)),
                    henry=float(match.group(4)),
                    henry_errorbar=float(match.group(5)),
                    source_data_file=str(path),
                )
            )
    return tuple(results)


def _write_results_csv(
    component_results: Sequence[GraspaWidomComponentResult],
    output_path: Path,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "component",
                "widom_energy",
                "widom_energy_errorbar",
                "henry",
                "henry_errorbar",
            ],
        )
        writer.writeheader()
        for result in component_results:
            writer.writerow(
                {
                    "component": result.component,
                    "widom_energy": result.widom_energy,
                    "widom_energy_errorbar": result.widom_energy_errorbar,
                    "henry": result.henry,
                    "henry_errorbar": result.henry_errorbar,
                }
            )
