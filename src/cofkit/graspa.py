from __future__ import annotations

from .calculation_io import atomic_write_text, begin_attempt, derive_seed, finish_attempt, record_execution, validate_numbers

import csv
import json
import math
import os
import re
import shutil
import subprocess
from dataclasses import dataclass, field, replace
from functools import lru_cache
from pathlib import Path
from typing import Mapping, Sequence

from ._dreiding_reference import (
    DREIDING_FRAMEWORK_TYPE_BY_ELEMENT,
    DREIDING_PARAMETERS,
    DREIDING_REFERENCE_SOURCE,
)
from .periodic_geometry import face_widths
from .cif_checks import read_ordered_structure, validate_charge_assignment
from .cofid import ensure_cif_has_cofid_comment, read_cofid_comment_from_cif
from .guest_bundles import GuestBundle, GuestBundleError, load_guest_bundles
from .guest_forcefields import (
    PACKAGED_GUEST_FORCEFIELD_SCHEMA_VERSION,
    load_packaged_guest_forcefield_metadata,
    packaged_guest_forcefield_catalog,
)
from .forcefields import (
    ForceFieldMetadataError,
    packaged_forcefield_artifact_path,
    resolve_forcefield_metadata,
    supported_forcefield_families,
)


COFKIT_EQEQ_ENV_VAR = "COFKIT_EQEQ_PATH"
COFKIT_GRASPA_ENV_VAR = "COFKIT_GRASPA_PATH"
COFKIT_RASPA2_ENV_VAR = "COFKIT_RASPA2_PATH"
DEFAULT_EQEQ_BINARY: Path | None = None
DEFAULT_GRASPA_BINARY: Path | None = None
DEFAULT_RASPA2_BINARY: Path | None = None
_SUPPORTED_EQEQ_METHODS = {"ewald", "nonperiodic"}
_SUPPORTED_GRASPA_FORCEFIELDS = supported_forcefield_families()
SUPPORTED_RASPA_BACKENDS = ("graspa", "raspa2")
DEFAULT_RASPA_BACKEND = "graspa"
PACKAGED_GUEST_FORCEFIELD_METADATA = load_packaged_guest_forcefield_metadata()
DEFAULT_WIDOM_COMPONENTS = tuple(
    metadata.name for metadata in PACKAGED_GUEST_FORCEFIELD_METADATA if metadata.default_widom and metadata.name != "H2_DREIDING"
)
AVAILABLE_WIDOM_COMPONENTS = tuple(metadata.name for metadata in PACKAGED_GUEST_FORCEFIELD_METADATA)
# Heuristic — pending calibration: per-component Widom insertion budget, on
# the order of the 2_000_000-cycle Widom production budget divided across the
# default component set.
DEFAULT_WIDOM_MOVES_PER_COMPONENT = 285_715
# RASPA-family run-control defaults shared by the Widom/isotherm/mixture
# settings dataclasses below (heuristic — pending calibration; values follow
# the RASPA2/gRASPA example-input convention).
DEFAULT_CUTOFF_ANGSTROM = 12.8
DEFAULT_OVERLAP_CRITERIA = 1.0e5
DEFAULT_EWALD_PRECISION = 1.0e-6
# derived: kcal/mol -> K via 1/R with R = 1.987191e-3 kcal/(mol*K); converts
# DREIDING/UFF LJ well depths into the Kelvin units gRASPA mixing rules expect.
_KCAL_PER_MOL_TO_KELVIN = 503.222681
_RASPA2_LOCAL_FORCEFIELD_NAME = "COFKit"
_RASPA2_LOCAL_MOLECULE_DEFINITION = "COFKit"
_NUMBER_RE = re.compile(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?")
_FLOAT_TOKEN_PATTERN = r"(?:[-+]?(?:nan|inf)|[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)"
# derived: Lennard-Jones equilibrium distance r_min = 2^(1/6) * sigma; converts
# UFF-tabulated r_min values into the sigma RASPA mixing-rule rows expect.
_UFF_RMIN_TO_SIGMA_FACTOR = 2.0 ** (1.0 / 6.0)
_UFF_FRAMEWORK_TYPE_BY_ELEMENT = {
    "H": "H_",
    "B": "B_2",
    "C": "C_3",
    "N": "N_3",
    "O": "O_3",
    "F": "F_",
    "P": "P_3+3",
    "S": "S_3+2",
    "Cl": "Cl",
    "Br": "Br",
}
_GRASPA_FRAMEWORK_ELEMENT_ORDER = ("H", "B", "C", "N", "O", "F", "P", "S", "Cl", "Br")
_GRASPA_SHARED_MIXING_RULE_LINES = (
    "Ow_T4          lennard-jones    77.99902   3.1536      // Jorgensen et al. J. Chem. Phys. (1983).",
    "Hw_T4          lennard-jones    0.0 0.0                                    // Jorgensen et al. J. Chem. Phys. (1983).",
    "M_T4           lennard-jones    0.0 0.0                                    // Jorgensen et al. J. Chem. Phys. (1983).",
    "O_co2          lennard-jones    79.0000   3.05000      // Potoff and Siepmann, AIChE J. 47, 1676-1682 (2001).",
    "C_co2          lennard-jones    27.0000   2.80000      // Potoff and Siepmann, AIChE J. 47, 1676-1682 (2001).",
    "N_n2           lennard-jones    38.2980   3.30600      // Martin-Calvo et al., PCCP 2011, 13, 11165-11174.",
    "N_com          lennard-jones    0.0 0.0                                   // Martin-Calvo et al., PCCP 2011, 13, 11165-11174.",
    "H_h2           lennard-jones    0.0 0.0                                    // Garcia-Holley et al., ACS Energy Lett. 2018, 3, 748-754.",
    "H2_com         feynman-hibbs-lennard-jones    36.7000   2.95800    1.008  // Garcia-Holley et al., ACS Energy Lett. 2018, 3, 748-754.",
    "S_so2          lennard-jones    73.8       3.39        // J. Phys. Chem. B 2011, 115, 4949-4954.",
    "O_so2          lennard-jones    79.0       3.05        // J. Phys. Chem. B 2011, 115, 4949-4954.",
    "Xe             lennard-jones   221.0000   4.10000      // RASPA GenericMOFs force field.",
    "Kr             lennard-jones   166.4000   3.63600      // RASPA GenericMOFs force field.",
)
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
_RASPA2_WIDOM_EXCESS_RE = re.compile(
    r"\[(?P<component>[^\]]+)\]\s+Average Widom excess chemical potential:\s*"
    rf"(?P<value>{_FLOAT_TOKEN_PATTERN})\s*\+/-\s*"
    rf"(?P<error>{_FLOAT_TOKEN_PATTERN})",
    re.IGNORECASE,
)
_RASPA2_HENRY_RE = re.compile(
    r"\[(?P<component>[^\]]+)\]\s+Average Henry coefficient:\s*"
    rf"(?P<value>{_FLOAT_TOKEN_PATTERN})\s*\+/-\s*"
    rf"(?P<error>{_FLOAT_TOKEN_PATTERN})\s*\[mol/kg/Pa\]",
    re.IGNORECASE,
)
_RASPA2_COMPONENT_HEADER_RE = re.compile(r"Component\s+\d+\s+\[([^\]]+)\]", re.IGNORECASE)
_RASPA2_AVERAGE_LOADING_MOL_KG_RE = re.compile(
    r"Average loading absolute \[mol/kg framework\]\s+"
    rf"({_FLOAT_TOKEN_PATTERN})\s*\+/-\s*({_FLOAT_TOKEN_PATTERN})",
    re.IGNORECASE,
)
_RASPA2_ENTHALPY_KJMOL_RE = re.compile(
    r"Enthalpy of adsorption:.*?Total enthalpy of adsorption.*?"
    rf"\n\s*({_FLOAT_TOKEN_PATTERN})\s*\+/-\s*({_FLOAT_TOKEN_PATTERN})\s*\[KJ/MOL\]",
    re.DOTALL | re.IGNORECASE,
)
# Content markers identifying which backend output schema a *.data file uses.
# Detection is deliberately content-based (never filename-based): gRASPA writes
# "Rosenbluth Summary For Component" Widom blocks and "BLOCK AVERAGES (...)"
# sections, while RASPA2 writes "[name] Average ..." rows and
# "Number of molecules:" sections. A file containing markers from both backends
# is a mixed-schema artifact (e.g. concatenated outputs) and must be parsed
# completely or diagnosed explicitly, never silently truncated.
_GRASPA_SCHEMA_MARKERS = (
    "rosenbluth summary for component",
    "block averages (",
)
_RASPA2_SCHEMA_MARKERS = (
    "average widom excess chemical potential:",
    "average henry coefficient:",
    "average loading absolute [mol/kg framework]",
)


def _detect_output_backend_schemas(content: str) -> tuple[str, ...]:
    """Return the backend output schemas whose markers appear in content."""
    lowered = content.lower()
    schemas: list[str] = []
    if any(marker in lowered for marker in _GRASPA_SCHEMA_MARKERS):
        schemas.append("graspa")
    if any(marker in lowered for marker in _RASPA2_SCHEMA_MARKERS):
        schemas.append("raspa2")
    return tuple(schemas)


def _value_tuples_conflict(
    first: tuple[float, float],
    second: tuple[float, float],
) -> bool:
    """True when two (value, errorbar) pairs carry different information."""
    for first_value, second_value in zip(first, second, strict=True):
        if math.isnan(first_value) and math.isnan(second_value):
            continue
        if not math.isclose(first_value, second_value, rel_tol=1e-9, abs_tol=1e-12):
            return True
    return False


class GraspaError(RuntimeError):
    """Base error for the gRASPA workflow wrappers."""


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
    # Cited: standard parameter values of the published extended charge
    # equilibration (EQeq) method and its reference implementation (Wilmer,
    # Kim & Snurr, J. Phys. Chem. Lett. 2012, 3, 2506-2511): dielectric
    # screening lambda = 1.2 (corresponds to an effective dielectric constant
    # of ~1.67) and hydrogen electron affinity -2.0 eV.
    lambda_value: float = 1.2
    hydrogen_electron_affinity: float = -2.0
    # Heuristic — pending calibration: decimal digits written for point
    # charges; more than the EQeq reference default for downstream reuse.
    charge_precision: int = 6
    # Total framework charge the equilibration targets; 0.0 = neutral.
    target_charge: float = 0.0
    # Heuristic — pending calibration: cofkit-side acceptance tolerance on
    # the net charge of the EQeq output.
    net_charge_tolerance: float = 1e-3
    method: str = "ewald"
    # Cited: EQeq reference-implementation CLI defaults — real-space and
    # reciprocal-space expansion cells for the Ewald summation.
    real_space_cells: int = 2
    reciprocal_space_cells: int = 2
    # Cited: EQeq reference-implementation CLI default for the atomic
    # hardness scaling parameter (eta).
    eta: float = 50.0

    def to_dict(self) -> dict[str, object]:
        return {
            "lambda_value": self.lambda_value,
            "hydrogen_electron_affinity": self.hydrogen_electron_affinity,
            "charge_precision": self.charge_precision,
            "target_charge": self.target_charge,
            "net_charge_tolerance": self.net_charge_tolerance,
            "method": self.method,
            "real_space_cells": self.real_space_cells,
            "reciprocal_space_cells": self.reciprocal_space_cells,
            "eta": self.eta,
        }


@dataclass(frozen=True)
class EqeqChargeResult:
    input_cif: str
    output_dir: str
    eqeq_binary: str
    eqeq_input_cif: str
    eqeq_charged_cif: str
    eqeq_json_output_path: str | None
    eqeq_stdout_log_path: str
    eqeq_stderr_log_path: str
    settings: EqeqChargeSettings

    def to_dict(self) -> dict[str, object]:
        return {
            "input_cif": self.input_cif,
            "output_dir": self.output_dir,
            "eqeq_binary": self.eqeq_binary,
            "eqeq_input_cif": self.eqeq_input_cif,
            "eqeq_charged_cif": self.eqeq_charged_cif,
            "eqeq_json_output_path": self.eqeq_json_output_path,
            "eqeq_stdout_log_path": self.eqeq_stdout_log_path,
            "eqeq_stderr_log_path": self.eqeq_stderr_log_path,
            "settings": self.settings.to_dict(),
        }


@dataclass(frozen=True)
class GraspaWidomSettings:
    components: tuple[str, ...] = DEFAULT_WIDOM_COMPONENTS
    guest_bundles: tuple[str, ...] = ()
    backend: str = DEFAULT_RASPA_BACKEND
    forcefield: str = "dreiding"
    use_gpu_reduction: bool = True
    use_fast_host_rng: bool = True
    use_flag: bool = True
    initialization_cycles: int = 0
    equilibration_cycles: int = 0
    production_cycles: int = 2_000_000
    widom_moves_per_component: int | None = None
    use_max_step: bool = True
    max_step_per_cycle: int = 1
    use_charges_from_cif_file: bool = True
    restart_file: bool = False
    random_seed: int = 0
    number_of_trial_positions: int = 10
    number_of_trial_orientations: int = 10
    number_of_blocks: int = 5
    adsorbate_allocate_space: int = 10_240
    number_of_simulations: int = 1
    single_simulation: bool = True
    different_frameworks: bool = True
    input_file_type: str = "cif"
    framework_name: str = "framework"
    charge_method: str = "Ewald"
    # Owner decision 2026-09-28 (MAGIC_NUMBER_FIX_PLAN.md W3.3): default
    # temperature unified at 298 K across all calculate workflows; the former
    # 300 K Widom / 298 K adsorption split is retired.
    temperature: float = 298.0
    pressure: float = 100_000.0
    overlap_criteria: float = DEFAULT_OVERLAP_CRITERIA
    cutoff_vdw: float = DEFAULT_CUTOFF_ANGSTROM
    cutoff_coulomb: float = DEFAULT_CUTOFF_ANGSTROM
    ewald_precision: float = DEFAULT_EWALD_PRECISION
    use_dnn_for_host_guest: bool = False
    save_output_to_file: bool = True

    def to_dict(self) -> dict[str, object]:
        return {
            "components": list(self.components),
            "guest_bundles": list(self.guest_bundles),
            "backend": self.backend,
            "forcefield": self.forcefield,
            "use_gpu_reduction": self.use_gpu_reduction,
            "use_fast_host_rng": self.use_fast_host_rng,
            "use_flag": self.use_flag,
            "initialization_cycles": self.initialization_cycles,
            "equilibration_cycles": self.equilibration_cycles,
            "production_cycles": self.production_cycles,
            "widom_moves_per_component": self.widom_moves_per_component,
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
    # Output schema ("graspa" or "raspa2") the component summary was parsed
    # from; recorded so mixed-schema inputs keep per-component provenance.
    source_schema: str | None = None

    def to_dict(self) -> dict[str, object]:
        return {
            "component": self.component,
            "widom_energy": _json_safe_float(self.widom_energy),
            "widom_energy_errorbar": _json_safe_float(self.widom_energy_errorbar),
            "henry": _json_safe_float(self.henry),
            "henry_errorbar": _json_safe_float(self.henry_errorbar),
            "source_data_file": self.source_data_file,
            "source_schema": self.source_schema,
        }


@dataclass(frozen=True)
class GraspaWidomResult:
    input_cif: str
    output_dir: str
    eqeq_binary: str
    graspa_binary: str
    raspa_backend: str
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
    framework_forcefield_metadata: dict[str, object] = field(default_factory=dict)

    def to_dict(self) -> dict[str, object]:
        return {
            "input_cif": self.input_cif,
            "output_dir": self.output_dir,
            "eqeq_binary": self.eqeq_binary,
            "graspa_binary": self.graspa_binary,
            "raspa_binary": self.graspa_binary,
            "raspa_backend": self.raspa_backend,
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
            "framework_forcefield_metadata": self.framework_forcefield_metadata,
            "framework_forcefield_metadata_paths": [
                str(Path(self.widom_run_dir) / "framework_forcefield_metadata.json")
            ],
            "guest_forcefield_metadata_paths": [
                str(Path(self.widom_run_dir) / "guest_forcefield_metadata.json")
            ],
            "component_results": [component.to_dict() for component in self.component_results],
            "warnings": list(self.warnings),
        }


@dataclass(frozen=True)
class GraspaIsothermSettings:
    component: str = "CO2_DREIDING"
    guest_bundles: tuple[str, ...] = ()
    pressures: tuple[float, ...] = (10_000.0, 100_000.0, 1_000_000.0)
    fugacity_coefficient: float | str = "PR-EOS"
    backend: str = DEFAULT_RASPA_BACKEND
    forcefield: str = "dreiding"
    use_gpu_reduction: bool = False
    use_fast_host_rng: bool = True
    use_flag: bool = True
    initialization_cycles: int = 50_000
    equilibration_cycles: int = 50_000
    production_cycles: int = 200_000
    use_max_step: bool = False
    max_step_per_cycle: int = 1
    use_charges_from_cif_file: bool = True
    restart_file: bool = False
    random_seed: int = 0
    number_of_trial_positions: int = 10
    number_of_trial_orientations: int = 10
    number_of_blocks: int = 5
    adsorbate_allocate_space: int = 10_240
    number_of_simulations: int = 1
    single_simulation: bool = True
    different_frameworks: bool = True
    input_file_type: str = "cif"
    framework_name: str = "framework"
    charge_method: str = "Ewald"
    temperature: float = 298.0
    overlap_criteria: float = DEFAULT_OVERLAP_CRITERIA
    cutoff_vdw: float = DEFAULT_CUTOFF_ANGSTROM
    cutoff_coulomb: float = DEFAULT_CUTOFF_ANGSTROM
    ewald_precision: float = DEFAULT_EWALD_PRECISION
    translation_probability: float = 1.0
    rotation_probability: float = 1.0
    reinsertion_probability: float = 1.0
    swap_probability: float = 1.0
    create_number_of_molecules: int = 0
    # None keeps the backend default snapshot cadence (cited: gRASPA writes
    # Movies/System_0/result_<cycle>.data during production whenever
    # cycle % MoviesEvery == 0, with 0-based cycles and a default
    # MoviesEvery of 5000 — gRASPA data_struct.h). Setting this renders an
    # explicit `MoviesEvery` line; only the gRASPA backend supports it.
    movies_every: int | None = None
    use_dnn_for_host_guest: bool = False
    save_output_to_file: bool = True

    def to_dict(self) -> dict[str, object]:
        return {
            "component": self.component,
            "guest_bundles": list(self.guest_bundles),
            "pressures": list(self.pressures),
            "fugacity_coefficient": self.fugacity_coefficient,
            "backend": self.backend,
            "forcefield": self.forcefield,
            "use_gpu_reduction": self.use_gpu_reduction,
            "use_fast_host_rng": self.use_fast_host_rng,
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
            "overlap_criteria": self.overlap_criteria,
            "cutoff_vdw": self.cutoff_vdw,
            "cutoff_coulomb": self.cutoff_coulomb,
            "ewald_precision": self.ewald_precision,
            "translation_probability": self.translation_probability,
            "rotation_probability": self.rotation_probability,
            "reinsertion_probability": self.reinsertion_probability,
            "swap_probability": self.swap_probability,
            "create_number_of_molecules": self.create_number_of_molecules,
            "movies_every": self.movies_every,
            "use_dnn_for_host_guest": self.use_dnn_for_host_guest,
            "save_output_to_file": self.save_output_to_file,
        }


@dataclass(frozen=True)
class GraspaIsothermPointResult:
    component: str
    pressure: float
    pressure_run_dir: str
    simulation_input_path: str
    graspa_stdout_log_path: str
    graspa_stderr_log_path: str
    data_file_paths: tuple[str, ...]
    source_data_file: str
    loading_mol_per_kg: float
    loading_mol_per_kg_errorbar: float
    loading_g_per_l: float
    loading_g_per_l_errorbar: float
    heat_of_adsorption_kj_per_mol: float
    heat_of_adsorption_kj_per_mol_errorbar: float
    initial_restart_file_path: str | None = None
    # Output schema ("graspa" or "raspa2") the summary was parsed from.
    source_schema: str | None = None
    warnings: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, object]:
        return {
            "component": self.component,
            "pressure": self.pressure,
            "pressure_run_dir": self.pressure_run_dir,
            "simulation_input_path": self.simulation_input_path,
            "graspa_stdout_log_path": self.graspa_stdout_log_path,
            "graspa_stderr_log_path": self.graspa_stderr_log_path,
            "data_file_paths": list(self.data_file_paths),
            "source_data_file": self.source_data_file,
            "source_schema": self.source_schema,
            "initial_restart_file_path": self.initial_restart_file_path,
            "loading_mol_per_kg": _json_safe_float(self.loading_mol_per_kg),
            "loading_mol_per_kg_errorbar": _json_safe_float(self.loading_mol_per_kg_errorbar),
            "loading_mmol_per_g": _json_safe_float(self.loading_mol_per_kg),
            "loading_mmol_per_g_errorbar": _json_safe_float(self.loading_mol_per_kg_errorbar),
            "loading_g_per_l": _json_safe_float(self.loading_g_per_l),
            "loading_g_per_l_errorbar": _json_safe_float(self.loading_g_per_l_errorbar),
            "heat_of_adsorption_kj_per_mol": _json_safe_float(self.heat_of_adsorption_kj_per_mol),
            "heat_of_adsorption_kj_per_mol_errorbar": _json_safe_float(self.heat_of_adsorption_kj_per_mol_errorbar),
            "warnings": list(self.warnings),
        }


@dataclass(frozen=True)
class GraspaIsothermResult:
    input_cif: str
    output_dir: str
    eqeq_binary: str
    graspa_binary: str
    raspa_backend: str
    eqeq_run_dir: str
    isotherm_root_dir: str
    eqeq_input_cif: str
    eqeq_charged_cif: str
    isotherm_framework_cif: str
    eqeq_json_output_path: str | None
    eqeq_stdout_log_path: str
    eqeq_stderr_log_path: str
    results_csv_path: str
    report_path: str
    unit_cells: tuple[int, int, int]
    eqeq_settings: EqeqChargeSettings
    isotherm_settings: GraspaIsothermSettings
    point_results: tuple[GraspaIsothermPointResult, ...]
    warnings: tuple[str, ...] = ()
    framework_forcefield_metadata: dict[str, object] = field(default_factory=dict)

    def to_dict(self) -> dict[str, object]:
        return {
            "input_cif": self.input_cif,
            "output_dir": self.output_dir,
            "eqeq_binary": self.eqeq_binary,
            "graspa_binary": self.graspa_binary,
            "raspa_binary": self.graspa_binary,
            "raspa_backend": self.raspa_backend,
            "eqeq_run_dir": self.eqeq_run_dir,
            "isotherm_root_dir": self.isotherm_root_dir,
            "eqeq_input_cif": self.eqeq_input_cif,
            "eqeq_charged_cif": self.eqeq_charged_cif,
            "isotherm_framework_cif": self.isotherm_framework_cif,
            "eqeq_json_output_path": self.eqeq_json_output_path,
            "eqeq_stdout_log_path": self.eqeq_stdout_log_path,
            "eqeq_stderr_log_path": self.eqeq_stderr_log_path,
            "results_csv_path": self.results_csv_path,
            "report_path": self.report_path,
            "unit_cells": list(self.unit_cells),
            "eqeq_settings": self.eqeq_settings.to_dict(),
            "isotherm_settings": self.isotherm_settings.to_dict(),
            "framework_forcefield_metadata": self.framework_forcefield_metadata,
            "framework_forcefield_metadata_paths": [
                str(Path(point.pressure_run_dir) / "framework_forcefield_metadata.json")
                for point in self.point_results
            ],
            "guest_forcefield_metadata_paths": [
                str(Path(point.pressure_run_dir) / "guest_forcefield_metadata.json") for point in self.point_results
            ],
            "point_results": [point.to_dict() for point in self.point_results],
            "warnings": list(self.warnings),
        }


@dataclass(frozen=True)
class GraspaMixtureComponentSettings:
    component: str
    mol_fraction: float
    fugacity_coefficient: float | str = "PR-EOS"
    translation_probability: float = 1.0
    rotation_probability: float = 1.0
    reinsertion_probability: float = 1.0
    identity_change_probability: float = 1.0
    swap_probability: float = 1.0
    create_number_of_molecules: int = 0

    def to_dict(self) -> dict[str, object]:
        return {
            "component": self.component,
            "mol_fraction": self.mol_fraction,
            "fugacity_coefficient": self.fugacity_coefficient,
            "translation_probability": self.translation_probability,
            "rotation_probability": self.rotation_probability,
            "reinsertion_probability": self.reinsertion_probability,
            "identity_change_probability": self.identity_change_probability,
            "swap_probability": self.swap_probability,
            "create_number_of_molecules": self.create_number_of_molecules,
        }


@dataclass(frozen=True)
class GraspaMixtureSettings:
    components: tuple[GraspaMixtureComponentSettings, ...]
    guest_bundles: tuple[str, ...] = ()
    pressures: tuple[float, ...] = (100_000.0,)
    backend: str = DEFAULT_RASPA_BACKEND
    forcefield: str = "dreiding"
    use_gpu_reduction: bool = False
    use_fast_host_rng: bool = True
    use_flag: bool = True
    initialization_cycles: int = 50_000
    equilibration_cycles: int = 50_000
    production_cycles: int = 200_000
    use_max_step: bool = False
    max_step_per_cycle: int = 1
    use_charges_from_cif_file: bool = True
    restart_file: bool = False
    random_seed: int = 0
    number_of_trial_positions: int = 10
    number_of_trial_orientations: int = 10
    number_of_blocks: int = 5
    adsorbate_allocate_space: int = 10_240
    number_of_simulations: int = 1
    single_simulation: bool = True
    different_frameworks: bool = True
    input_file_type: str = "cif"
    framework_name: str = "framework"
    charge_method: str = "Ewald"
    temperature: float = 298.0
    overlap_criteria: float = DEFAULT_OVERLAP_CRITERIA
    cutoff_vdw: float = DEFAULT_CUTOFF_ANGSTROM
    cutoff_coulomb: float = DEFAULT_CUTOFF_ANGSTROM
    ewald_precision: float = DEFAULT_EWALD_PRECISION
    # None keeps the backend default snapshot cadence (cited: gRASPA writes
    # Movies/System_0/result_<cycle>.data during production whenever
    # cycle % MoviesEvery == 0, with 0-based cycles and a default
    # MoviesEvery of 5000 — gRASPA data_struct.h). Setting this renders an
    # explicit `MoviesEvery` line; only the gRASPA backend supports it.
    movies_every: int | None = None
    use_dnn_for_host_guest: bool = False
    save_output_to_file: bool = True

    def to_dict(self) -> dict[str, object]:
        return {
            "components": [component.to_dict() for component in self.components],
            "guest_bundles": list(self.guest_bundles),
            "pressures": list(self.pressures),
            "backend": self.backend,
            "forcefield": self.forcefield,
            "use_gpu_reduction": self.use_gpu_reduction,
            "use_fast_host_rng": self.use_fast_host_rng,
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
            "overlap_criteria": self.overlap_criteria,
            "cutoff_vdw": self.cutoff_vdw,
            "cutoff_coulomb": self.cutoff_coulomb,
            "ewald_precision": self.ewald_precision,
            "movies_every": self.movies_every,
            "use_dnn_for_host_guest": self.use_dnn_for_host_guest,
            "save_output_to_file": self.save_output_to_file,
        }


@dataclass(frozen=True)
class GraspaMixturePointComponentResult:
    component: str
    feed_mol_fraction: float
    adsorbed_mol_fraction: float
    loading_mol_per_kg: float
    loading_mol_per_kg_errorbar: float
    loading_g_per_l: float
    loading_g_per_l_errorbar: float
    heat_of_adsorption_kj_per_mol: float
    heat_of_adsorption_kj_per_mol_errorbar: float

    def to_dict(self) -> dict[str, object]:
        return {
            "component": self.component,
            "feed_mol_fraction": _json_safe_float(self.feed_mol_fraction),
            "adsorbed_mol_fraction": _json_safe_float(self.adsorbed_mol_fraction),
            "loading_mol_per_kg": _json_safe_float(self.loading_mol_per_kg),
            "loading_mol_per_kg_errorbar": _json_safe_float(self.loading_mol_per_kg_errorbar),
            "loading_mmol_per_g": _json_safe_float(self.loading_mol_per_kg),
            "loading_mmol_per_g_errorbar": _json_safe_float(self.loading_mol_per_kg_errorbar),
            "loading_g_per_l": _json_safe_float(self.loading_g_per_l),
            "loading_g_per_l_errorbar": _json_safe_float(self.loading_g_per_l_errorbar),
            "heat_of_adsorption_kj_per_mol": _json_safe_float(self.heat_of_adsorption_kj_per_mol),
            "heat_of_adsorption_kj_per_mol_errorbar": _json_safe_float(
                self.heat_of_adsorption_kj_per_mol_errorbar
            ),
        }


@dataclass(frozen=True)
class GraspaMixtureSelectivityResult:
    numerator_component: str
    denominator_component: str
    feed_mol_fraction_ratio: float
    adsorbed_mol_fraction_ratio: float
    selectivity: float
    selectivity_errorbar: float

    def to_dict(self) -> dict[str, object]:
        return {
            "uncertainty_status": "unavailable_without_paired_loading_statistics",
            "numerator_component": self.numerator_component,
            "denominator_component": self.denominator_component,
            "feed_mol_fraction_ratio": _json_safe_float(self.feed_mol_fraction_ratio),
            "adsorbed_mol_fraction_ratio": _json_safe_float(self.adsorbed_mol_fraction_ratio),
            "selectivity": _json_safe_float(self.selectivity),
            "selectivity_errorbar": _json_safe_float(self.selectivity_errorbar),
        }


@dataclass(frozen=True)
class GraspaMixturePointResult:
    pressure: float
    pressure_run_dir: str
    simulation_input_path: str
    graspa_stdout_log_path: str
    graspa_stderr_log_path: str
    data_file_paths: tuple[str, ...]
    source_data_file: str
    component_results: tuple[GraspaMixturePointComponentResult, ...]
    selectivity_results: tuple[GraspaMixtureSelectivityResult, ...]
    initial_restart_file_path: str | None = None
    warnings: tuple[str, ...] = ()
    # Per-component output schema ("graspa" or "raspa2") provenance; a mixed
    # value set means the source file combined both backend layouts.
    component_source_schemas: Mapping[str, str] = field(default_factory=dict)

    def to_dict(self) -> dict[str, object]:
        return {
            "pressure": self.pressure,
            "pressure_run_dir": self.pressure_run_dir,
            "simulation_input_path": self.simulation_input_path,
            "graspa_stdout_log_path": self.graspa_stdout_log_path,
            "graspa_stderr_log_path": self.graspa_stderr_log_path,
            "data_file_paths": list(self.data_file_paths),
            "source_data_file": self.source_data_file,
            "initial_restart_file_path": self.initial_restart_file_path,
            "component_results": [component.to_dict() for component in self.component_results],
            "selectivity_results": [result.to_dict() for result in self.selectivity_results],
            "warnings": list(self.warnings),
            "component_source_schemas": dict(self.component_source_schemas),
        }


@dataclass(frozen=True)
class GraspaMixtureResult:
    input_cif: str
    output_dir: str
    eqeq_binary: str
    graspa_binary: str
    raspa_backend: str
    eqeq_run_dir: str
    mixture_root_dir: str
    eqeq_input_cif: str
    eqeq_charged_cif: str
    mixture_framework_cif: str
    eqeq_json_output_path: str | None
    eqeq_stdout_log_path: str
    eqeq_stderr_log_path: str
    component_results_csv_path: str
    selectivity_results_csv_path: str
    report_path: str
    unit_cells: tuple[int, int, int]
    eqeq_settings: EqeqChargeSettings
    mixture_settings: GraspaMixtureSettings
    point_results: tuple[GraspaMixturePointResult, ...]
    warnings: tuple[str, ...] = ()
    framework_forcefield_metadata: dict[str, object] = field(default_factory=dict)

    def to_dict(self) -> dict[str, object]:
        return {
            "input_cif": self.input_cif,
            "output_dir": self.output_dir,
            "eqeq_binary": self.eqeq_binary,
            "graspa_binary": self.graspa_binary,
            "raspa_binary": self.graspa_binary,
            "raspa_backend": self.raspa_backend,
            "eqeq_run_dir": self.eqeq_run_dir,
            "mixture_root_dir": self.mixture_root_dir,
            "eqeq_input_cif": self.eqeq_input_cif,
            "eqeq_charged_cif": self.eqeq_charged_cif,
            "mixture_framework_cif": self.mixture_framework_cif,
            "eqeq_json_output_path": self.eqeq_json_output_path,
            "eqeq_stdout_log_path": self.eqeq_stdout_log_path,
            "eqeq_stderr_log_path": self.eqeq_stderr_log_path,
            "component_results_csv_path": self.component_results_csv_path,
            "selectivity_results_csv_path": self.selectivity_results_csv_path,
            "report_path": self.report_path,
            "unit_cells": list(self.unit_cells),
            "eqeq_settings": self.eqeq_settings.to_dict(),
            "mixture_settings": self.mixture_settings.to_dict(),
            "framework_forcefield_metadata": self.framework_forcefield_metadata,
            "framework_forcefield_metadata_paths": [
                str(Path(point.pressure_run_dir) / "framework_forcefield_metadata.json")
                for point in self.point_results
            ],
            "guest_forcefield_metadata_paths": [
                str(Path(point.pressure_run_dir) / "guest_forcefield_metadata.json") for point in self.point_results
            ],
            "point_results": [point.to_dict() for point in self.point_results],
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


def resolve_raspa2_binary(raspa2_path: str | Path | None = None) -> Path:
    return _resolve_binary(
        explicit_path=raspa2_path,
        env_var=COFKIT_RASPA2_ENV_VAR,
        default_path=DEFAULT_RASPA2_BINARY,
        display_name="RASPA2",
    )


def resolve_raspa_backend_binary(
    backend: str,
    *,
    raspa_path: str | Path | None = None,
    graspa_path: str | Path | None = None,
    raspa2_path: str | Path | None = None,
) -> Path:
    normalized = _normalize_raspa_backend_name(backend)
    if normalized == "graspa":
        return resolve_graspa_binary(raspa_path if raspa_path is not None else graspa_path)
    if normalized == "raspa2":
        return resolve_raspa2_binary(raspa_path if raspa_path is not None else raspa2_path)
    raise ValueError(
        f"Unsupported RASPA backend {backend!r}. Expected one of: {', '.join(SUPPORTED_RASPA_BACKENDS)}."
    )


def assign_eqeq_charges_to_cif(
    cif_path: str | Path,
    *,
    output_dir: str | Path | None = None,
    eqeq_path: str | Path | None = None,
    settings: EqeqChargeSettings | None = None,
    timeout_seconds: float | None = 300.0,
) -> EqeqChargeResult:
    input_path = Path(cif_path).expanduser().resolve()
    if not input_path.is_file():
        raise FileNotFoundError(f"CIF file does not exist: {input_path}")
    input_cofid_comment = read_cofid_comment_from_cif(input_path)

    settings = settings or EqeqChargeSettings()
    _validate_eqeq_settings(settings)
    timeout_seconds = _normalize_timeout(timeout_seconds, "timeout_seconds")
    read_ordered_structure(input_path)
    eqeq_binary = resolve_eqeq_binary(eqeq_path)

    run_dir = (
        Path(output_dir).expanduser().resolve()
        if output_dir is not None
        else input_path.parent / f"{input_path.stem}_eqeq"
    )
    run_dir = begin_attempt(run_dir, input_path=input_path, settings=settings, binary=eqeq_binary)

    eqeq_input_path = run_dir / input_path.name
    shutil.copy2(input_path, eqeq_input_path)
    eqeq_stdout_log_path = run_dir / "eqeq.stdout.log"
    eqeq_stderr_log_path = run_dir / "eqeq.stderr.log"

    eqeq_command = [
        str(eqeq_binary),
        eqeq_input_path.name,
        _format_cli_number(settings.lambda_value),
        _format_cli_number(settings.hydrogen_electron_affinity),
        str(settings.charge_precision),
        settings.method,
        str(settings.real_space_cells),
        str(settings.reciprocal_space_cells),
        _format_cli_number(settings.eta),
    ]
    record_execution(run_dir, eqeq_command)
    try:
        with eqeq_stdout_log_path.open("w", encoding="utf-8") as stdout_handle:
            with eqeq_stderr_log_path.open("w", encoding="utf-8") as stderr_handle:
                eqeq_completed = subprocess.run(
                    eqeq_command,
                    cwd=run_dir,
                    check=False,
                    stdout=stdout_handle,
                    stderr=stderr_handle,
                    text=True,
                    timeout=timeout_seconds,
                )
    except subprocess.TimeoutExpired as exc:
        raise EqeqExecutionError(
            f"EQeq timed out after {timeout_seconds} seconds. "
            f"See {eqeq_stdout_log_path} and {eqeq_stderr_log_path}."
        ) from exc

    eqeq_output_stem = _expected_eqeq_output_stem(eqeq_input_path.name, settings)
    eqeq_charged_cif_path = run_dir / f"{eqeq_output_stem}.cif"
    eqeq_json_output_path = run_dir / f"{eqeq_output_stem}.json"

    if eqeq_completed.returncode != 0:
        raise EqeqExecutionError(
            f"EQeq failed with exit code {eqeq_completed.returncode}. "
            f"See {eqeq_stdout_log_path} and {eqeq_stderr_log_path}."
        )
    if not eqeq_charged_cif_path.is_file():
        raise EqeqExecutionError(
            f"EQeq completed without writing the expected charged CIF: {eqeq_charged_cif_path}"
        )
    try:
        charge_checks = validate_charge_assignment(input_path, eqeq_charged_cif_path,
            target_charge=settings.target_charge, tolerance=settings.net_charge_tolerance)
    except ValueError as exc:
        raise EqeqExecutionError(str(exc)) from exc
    atomic_write_text(run_dir / "charge_validation.json", json.dumps(charge_checks, indent=2))
    ensure_cif_has_cofid_comment(
        eqeq_charged_cif_path,
        None if input_cofid_comment is None else input_cofid_comment.cofid,
        suffix=None if input_cofid_comment is None else input_cofid_comment.suffix,
    )

    finish_attempt(run_dir)
    return EqeqChargeResult(
        input_cif=str(input_path),
        output_dir=str(run_dir),
        eqeq_binary=str(eqeq_binary),
        eqeq_input_cif=str(eqeq_input_path),
        eqeq_charged_cif=str(eqeq_charged_cif_path),
        eqeq_json_output_path=str(eqeq_json_output_path) if eqeq_json_output_path.is_file() else None,
        eqeq_stdout_log_path=str(eqeq_stdout_log_path),
        eqeq_stderr_log_path=str(eqeq_stderr_log_path),
        settings=settings,
    )


def run_graspa_widom_workflow(
    cif_path: str | Path,
    *,
    output_dir: str | Path | None = None,
    eqeq_path: str | Path | None = None,
    graspa_path: str | Path | None = None,
    raspa_path: str | Path | None = None,
    raspa2_path: str | Path | None = None,
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
    guest_bundles = _load_graspa_guest_bundles(widom_settings.guest_bundles)
    widom_settings = _canonicalize_widom_settings(widom_settings, guest_bundles)
    _validate_widom_settings(widom_settings, guest_bundles=guest_bundles)
    graspa_timeout_seconds = _normalize_timeout(graspa_timeout_seconds, "graspa_timeout_seconds")

    raspa_backend = _normalize_raspa_backend_name(widom_settings.backend)
    raspa_display_name = _raspa_backend_display_name(raspa_backend)
    graspa_binary = resolve_raspa_backend_binary(
        raspa_backend,
        raspa_path=raspa_path,
        graspa_path=graspa_path,
        raspa2_path=raspa2_path,
    )

    run_dir = _resolve_output_dir(input_path, output_dir, suffix=f"_{raspa_backend}_widom")
    run_dir = begin_attempt(run_dir, input_path=input_path, settings=widom_settings, binary=graspa_binary)
    eqeq_run_dir = run_dir / "eqeq"
    widom_run_dir = run_dir / "widom"
    widom_run_dir.mkdir(parents=True, exist_ok=True)
    eqeq_result = assign_eqeq_charges_to_cif(
        input_path,
        output_dir=eqeq_run_dir,
        eqeq_path=eqeq_path,
        settings=eqeq_settings,
        timeout_seconds=eqeq_timeout_seconds,
    )

    widom_framework_cif_path = widom_run_dir / f"{widom_settings.framework_name}.cif"
    _write_graspa_framework_cif(eqeq_result.eqeq_charged_cif, widom_framework_cif_path)
    _copy_widom_template_assets(
        widom_run_dir,
        widom_settings.components,
        forcefield=widom_settings.forcefield,
        guest_bundles=guest_bundles,
    )

    unit_cells = _compute_unit_cells_from_cif(
        widom_framework_cif_path,
        cutoff=max(widom_settings.cutoff_vdw, widom_settings.cutoff_coulomb),
    )
    simulation_input_path = widom_run_dir / "simulation.input"
    simulation_input_path.write_text(
        _render_widom_simulation_input(widom_settings, unit_cells=unit_cells, guest_bundles=guest_bundles),
        encoding="utf-8",
    )

    graspa_stdout_log_path = widom_run_dir / f"{raspa_backend}.stdout.log"
    graspa_stderr_log_path = widom_run_dir / f"{raspa_backend}.stderr.log"
    record_execution(widom_run_dir, [str(graspa_binary)])
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
            f"{raspa_display_name} timed out after {graspa_timeout_seconds} seconds. "
            f"See {graspa_stdout_log_path} and {graspa_stderr_log_path}."
        ) from exc

    if graspa_completed.returncode != 0:
        raise GraspaExecutionError(
            f"{raspa_display_name} failed with exit code {graspa_completed.returncode}. "
            f"See {graspa_stdout_log_path} and {graspa_stderr_log_path}."
        )

    widom_output_dir = widom_run_dir / "Output"
    data_file_paths = _collect_raspa_data_file_paths(widom_output_dir)
    if not data_file_paths:
        raise GraspaParseError(
            f"{raspa_display_name} completed without producing any Output/**/*.data files under {widom_output_dir}"
        )

    parse_diagnostics: list[str] = []
    component_results = _parse_widom_result_files(
        tuple(Path(path) for path in data_file_paths),
        diagnostics=parse_diagnostics,
    )
    if not component_results:
        raise GraspaParseError(
            f"No Widom component summaries could be parsed from {len(data_file_paths)} data file(s)."
        )

    requested_components = tuple(widom_settings.components)
    parsed_components = tuple(result.component for result in component_results)
    if set(parsed_components) != set(requested_components) or len(parsed_components) != len(requested_components):
        missing_components = sorted(set(requested_components) - set(parsed_components))
        unexpected_components = sorted(set(parsed_components) - set(requested_components))
        duplicate_components = sorted(
            {component for component in parsed_components if parsed_components.count(component) > 1}
        )
        raise GraspaParseError(
            "Widom summaries do not match the requested components exactly: "
            f"requested {list(requested_components)!r}; missing {missing_components!r}; "
            f"unexpected {unexpected_components!r}; duplicated {duplicate_components!r}."
        )
    results_csv_path = widom_output_dir / "results.csv"
    _write_results_csv(component_results, results_csv_path)

    warnings: list[str] = []
    warnings.extend(parse_diagnostics)
    if widom_settings.components == DEFAULT_WIDOM_COMPONENTS:
        warnings.append("The default probe set excludes H2_DREIDING: its Feynman-Hibbs potential is unsupported by upstream gRASPA.")
    if eqeq_result.eqeq_json_output_path is None:
        warnings.append("EQeq did not write the companion JSON output file.")
    if len(data_file_paths) > 1:
        warnings.append(f"Parsed Widom results from {len(data_file_paths)} data files.")
    if graspa_stderr_log_path.read_text(encoding="utf-8", errors="replace").strip():
        warnings.append(f"{raspa_display_name} wrote content to stderr; inspect the stderr log if needed.")
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
        warnings.append(f"One or more Widom summary values are non-finite; inspect the raw {raspa_display_name} data output.")

    report_path = run_dir / "graspa_widom_report.json"
    result = GraspaWidomResult(
        input_cif=str(input_path),
        output_dir=str(run_dir),
        eqeq_binary=eqeq_result.eqeq_binary,
        graspa_binary=str(graspa_binary),
        raspa_backend=raspa_backend,
        eqeq_run_dir=eqeq_result.output_dir,
        widom_run_dir=str(widom_run_dir),
        eqeq_input_cif=eqeq_result.eqeq_input_cif,
        eqeq_charged_cif=eqeq_result.eqeq_charged_cif,
        widom_framework_cif=str(widom_framework_cif_path),
        eqeq_json_output_path=eqeq_result.eqeq_json_output_path,
        eqeq_stdout_log_path=eqeq_result.eqeq_stdout_log_path,
        eqeq_stderr_log_path=eqeq_result.eqeq_stderr_log_path,
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
        framework_forcefield_metadata=resolve_forcefield_metadata(widom_settings.forcefield).to_dict(),
    )
    atomic_write_text(report_path, json.dumps(result.to_dict(), indent=2, allow_nan=False))
    finish_attempt(run_dir)
    return result


def run_graspa_isotherm_workflow(
    cif_path: str | Path,
    *,
    output_dir: str | Path | None = None,
    eqeq_path: str | Path | None = None,
    graspa_path: str | Path | None = None,
    raspa_path: str | Path | None = None,
    raspa2_path: str | Path | None = None,
    eqeq_settings: EqeqChargeSettings | None = None,
    isotherm_settings: GraspaIsothermSettings | None = None,
    initial_restart_file: str | Path | None = None,
    eqeq_timeout_seconds: float | None = 300.0,
    graspa_timeout_seconds: float | None = None,
) -> GraspaIsothermResult:
    input_path = Path(cif_path).expanduser().resolve()
    if not input_path.is_file():
        raise FileNotFoundError(f"CIF file does not exist: {input_path}")

    eqeq_settings = eqeq_settings or EqeqChargeSettings()
    isotherm_settings = isotherm_settings or GraspaIsothermSettings()
    guest_bundles = _load_graspa_guest_bundles(isotherm_settings.guest_bundles)
    isotherm_settings = _canonicalize_isotherm_settings(isotherm_settings, guest_bundles)
    _validate_isotherm_settings(isotherm_settings, guest_bundles=guest_bundles)
    graspa_timeout_seconds = _normalize_timeout(graspa_timeout_seconds, "graspa_timeout_seconds")

    raspa_backend = _normalize_raspa_backend_name(isotherm_settings.backend)
    raspa_display_name = _raspa_backend_display_name(raspa_backend)
    if initial_restart_file is not None and raspa_backend != "graspa":
        raise GraspaConfigurationError(
            "Initial guest restart files are currently supported only for the gRASPA backend; "
            "RASPA2 restart-file staging has not been validated."
        )
    graspa_binary = resolve_raspa_backend_binary(
        raspa_backend,
        raspa_path=raspa_path,
        graspa_path=graspa_path,
        raspa2_path=raspa2_path,
    )

    run_dir = _resolve_output_dir(input_path, output_dir, suffix=f"_{raspa_backend}_isotherm")
    run_dir = begin_attempt(run_dir, input_path=input_path, settings=isotherm_settings, binary=graspa_binary)
    eqeq_run_dir = run_dir / "eqeq"
    isotherm_root_dir = run_dir / "isotherm"
    isotherm_root_dir.mkdir(parents=True, exist_ok=True)
    eqeq_result = assign_eqeq_charges_to_cif(
        input_path,
        output_dir=eqeq_run_dir,
        eqeq_path=eqeq_path,
        settings=eqeq_settings,
        timeout_seconds=eqeq_timeout_seconds,
    )

    isotherm_framework_cif_path = isotherm_root_dir / f"{isotherm_settings.framework_name}.cif"
    _write_graspa_framework_cif(eqeq_result.eqeq_charged_cif, isotherm_framework_cif_path)

    unit_cells = _compute_unit_cells_from_cif(
        isotherm_framework_cif_path,
        cutoff=max(isotherm_settings.cutoff_vdw, isotherm_settings.cutoff_coulomb),
    )

    point_results: list[GraspaIsothermPointResult] = []
    warnings: list[str] = []

    for index, pressure in enumerate(isotherm_settings.pressures):
        point_settings = replace(isotherm_settings, random_seed=derive_seed(isotherm_settings.random_seed, "pressure", index, pressure))
        pressure_run_dir = isotherm_root_dir / _format_pressure_run_dir_name(index, pressure)
        pressure_run_dir.mkdir(parents=True, exist_ok=True)
        pressure_framework_cif_path = pressure_run_dir / f"{isotherm_settings.framework_name}.cif"
        shutil.copy2(isotherm_framework_cif_path, pressure_framework_cif_path)
        _copy_widom_template_assets(
            pressure_run_dir,
            (isotherm_settings.component,),
            forcefield=isotherm_settings.forcefield,
            guest_bundles=guest_bundles,
        )
        staged_initial_restart_file = (
            _stage_graspa_initial_restart_file(pressure_run_dir, initial_restart_file)
            if initial_restart_file is not None
            else None
        )

        simulation_input_path = pressure_run_dir / "simulation.input"
        simulation_input_path.write_text(
            _render_isotherm_simulation_input(
                point_settings,
                unit_cells=unit_cells,
                pressure=pressure,
                guest_bundles=guest_bundles,
            ),
            encoding="utf-8",
        )

        graspa_stdout_log_path = pressure_run_dir / f"{raspa_backend}.stdout.log"
        graspa_stderr_log_path = pressure_run_dir / f"{raspa_backend}.stderr.log"
        record_execution(pressure_run_dir, [str(graspa_binary)])
        try:
            with graspa_stdout_log_path.open("w", encoding="utf-8") as stdout_handle:
                with graspa_stderr_log_path.open("w", encoding="utf-8") as stderr_handle:
                    graspa_completed = subprocess.run(
                        [str(graspa_binary)],
                        cwd=pressure_run_dir,
                        check=False,
                        stdout=stdout_handle,
                        stderr=stderr_handle,
                        text=True,
                        timeout=graspa_timeout_seconds,
                    )
        except subprocess.TimeoutExpired as exc:
            raise GraspaExecutionError(
                f"{raspa_display_name} timed out after {graspa_timeout_seconds} seconds while running "
                f"pressure {pressure:g} Pa. See {graspa_stdout_log_path} and {graspa_stderr_log_path}."
            ) from exc

        if graspa_completed.returncode != 0:
            raise GraspaExecutionError(
                f"{raspa_display_name} failed with exit code {graspa_completed.returncode} at pressure {pressure:g} Pa. "
                f"See {graspa_stdout_log_path} and {graspa_stderr_log_path}."
            )

        data_file_paths = _collect_raspa_data_file_paths(pressure_run_dir / "Output")
        if not data_file_paths:
            raise GraspaParseError(
                f"{raspa_display_name} completed without producing any Output/**/*.data files under {pressure_run_dir / 'Output'} "
                f"for pressure {pressure:g} Pa."
            )

        point_result = _parse_isotherm_result_files(
            tuple(Path(path) for path in data_file_paths),
            component=isotherm_settings.component,
            pressure=pressure,
            pressure_run_dir=pressure_run_dir,
            simulation_input_path=simulation_input_path,
            graspa_stdout_log_path=graspa_stdout_log_path,
            graspa_stderr_log_path=graspa_stderr_log_path,
        )
        if staged_initial_restart_file is not None:
            point_result = replace(
                point_result,
                initial_restart_file_path=str(staged_initial_restart_file),
            )
        point_results.append(point_result)
        warnings.extend(point_result.warnings)

        if len(data_file_paths) > 1:
            warnings.append(
                f"Parsed isotherm point {pressure:g} Pa from {len(data_file_paths)} data files."
            )
        if graspa_stderr_log_path.read_text(encoding="utf-8", errors="replace").strip():
            warnings.append(
                f"{raspa_display_name} wrote content to stderr for pressure {pressure:g} Pa; inspect the stderr log if needed."
            )
        if any(
            not math.isfinite(value)
            for value in (
                point_result.loading_mol_per_kg,
                point_result.loading_mol_per_kg_errorbar,
                point_result.loading_g_per_l,
                point_result.loading_g_per_l_errorbar,
                point_result.heat_of_adsorption_kj_per_mol,
                point_result.heat_of_adsorption_kj_per_mol_errorbar,
            )
        ):
            warnings.append(
                f"One or more parsed isotherm values are non-finite at pressure {pressure:g} Pa; inspect the raw {raspa_display_name} data output."
            )

    if eqeq_result.eqeq_json_output_path is None:
        warnings.append("EQeq did not write the companion JSON output file.")

    results_csv_path = isotherm_root_dir / "results.csv"
    _write_isotherm_results_csv(point_results, results_csv_path)

    report_path = run_dir / "graspa_isotherm_report.json"
    result = GraspaIsothermResult(
        input_cif=str(input_path),
        output_dir=str(run_dir),
        eqeq_binary=eqeq_result.eqeq_binary,
        graspa_binary=str(graspa_binary),
        raspa_backend=raspa_backend,
        eqeq_run_dir=eqeq_result.output_dir,
        isotherm_root_dir=str(isotherm_root_dir),
        eqeq_input_cif=eqeq_result.eqeq_input_cif,
        eqeq_charged_cif=eqeq_result.eqeq_charged_cif,
        isotherm_framework_cif=str(isotherm_framework_cif_path),
        eqeq_json_output_path=eqeq_result.eqeq_json_output_path,
        eqeq_stdout_log_path=eqeq_result.eqeq_stdout_log_path,
        eqeq_stderr_log_path=eqeq_result.eqeq_stderr_log_path,
        results_csv_path=str(results_csv_path),
        report_path=str(report_path),
        unit_cells=unit_cells,
        eqeq_settings=eqeq_settings,
        isotherm_settings=isotherm_settings,
        point_results=tuple(point_results),
        warnings=tuple(warnings),
        framework_forcefield_metadata=resolve_forcefield_metadata(isotherm_settings.forcefield).to_dict(),
    )
    atomic_write_text(report_path, json.dumps(result.to_dict(), indent=2, allow_nan=False))
    finish_attempt(run_dir)
    return result


def run_graspa_mixture_workflow(
    cif_path: str | Path,
    *,
    output_dir: str | Path | None = None,
    eqeq_path: str | Path | None = None,
    graspa_path: str | Path | None = None,
    raspa_path: str | Path | None = None,
    raspa2_path: str | Path | None = None,
    eqeq_settings: EqeqChargeSettings | None = None,
    mixture_settings: GraspaMixtureSettings | None = None,
    initial_restart_file: str | Path | None = None,
    eqeq_timeout_seconds: float | None = 300.0,
    graspa_timeout_seconds: float | None = None,
) -> GraspaMixtureResult:
    input_path = Path(cif_path).expanduser().resolve()
    if not input_path.is_file():
        raise FileNotFoundError(f"CIF file does not exist: {input_path}")

    eqeq_settings = eqeq_settings or EqeqChargeSettings()
    if mixture_settings is None:
        raise ValueError("mixture_settings must be provided for the gRASPA mixture workflow.")
    guest_bundles = _load_graspa_guest_bundles(mixture_settings.guest_bundles)
    mixture_settings = _canonicalize_mixture_settings(mixture_settings, guest_bundles)
    _validate_mixture_settings(mixture_settings, guest_bundles=guest_bundles)
    graspa_timeout_seconds = _normalize_timeout(graspa_timeout_seconds, "graspa_timeout_seconds")

    raspa_backend = _normalize_raspa_backend_name(mixture_settings.backend)
    raspa_display_name = _raspa_backend_display_name(raspa_backend)
    if initial_restart_file is not None and raspa_backend != "graspa":
        raise GraspaConfigurationError(
            "Initial guest restart files are currently supported only for the gRASPA backend; "
            "RASPA2 restart-file staging has not been validated."
        )
    graspa_binary = resolve_raspa_backend_binary(
        raspa_backend,
        raspa_path=raspa_path,
        graspa_path=graspa_path,
        raspa2_path=raspa2_path,
    )

    run_dir = _resolve_output_dir(input_path, output_dir, suffix=f"_{raspa_backend}_mixture")
    run_dir = begin_attempt(run_dir, input_path=input_path, settings=mixture_settings, binary=graspa_binary)
    eqeq_run_dir = run_dir / "eqeq"
    mixture_root_dir = run_dir / "mixture"
    mixture_root_dir.mkdir(parents=True, exist_ok=True)
    eqeq_result = assign_eqeq_charges_to_cif(
        input_path,
        output_dir=eqeq_run_dir,
        eqeq_path=eqeq_path,
        settings=eqeq_settings,
        timeout_seconds=eqeq_timeout_seconds,
    )

    mixture_framework_cif_path = mixture_root_dir / f"{mixture_settings.framework_name}.cif"
    _write_graspa_framework_cif(eqeq_result.eqeq_charged_cif, mixture_framework_cif_path)

    unit_cells = _compute_unit_cells_from_cif(
        mixture_framework_cif_path,
        cutoff=max(mixture_settings.cutoff_vdw, mixture_settings.cutoff_coulomb),
    )

    point_results: list[GraspaMixturePointResult] = []
    warnings: list[str] = []
    expected_components = tuple(component.component for component in mixture_settings.components)
    raw_mol_fraction_total = sum(component.mol_fraction for component in mixture_settings.components)
    if not math.isclose(raw_mol_fraction_total, 1.0, rel_tol=1e-9, abs_tol=1e-9):
        warnings.append(
            f"Mixture component mol fractions sum to {raw_mol_fraction_total:g}; gRASPA normalizes them internally and cofkit reports normalized feed fractions."
        )

    for index, pressure in enumerate(mixture_settings.pressures):
        point_settings = replace(mixture_settings, random_seed=derive_seed(mixture_settings.random_seed, "pressure", index, pressure))
        pressure_run_dir = mixture_root_dir / _format_pressure_run_dir_name(index, pressure)
        pressure_run_dir.mkdir(parents=True, exist_ok=True)
        pressure_framework_cif_path = pressure_run_dir / f"{mixture_settings.framework_name}.cif"
        shutil.copy2(mixture_framework_cif_path, pressure_framework_cif_path)
        _copy_widom_template_assets(
            pressure_run_dir,
            expected_components,
            forcefield=mixture_settings.forcefield,
            guest_bundles=guest_bundles,
        )
        staged_initial_restart_file = (
            _stage_graspa_initial_restart_file(pressure_run_dir, initial_restart_file)
            if initial_restart_file is not None
            else None
        )

        simulation_input_path = pressure_run_dir / "simulation.input"
        simulation_input_path.write_text(
            _render_mixture_simulation_input(
                point_settings,
                unit_cells=unit_cells,
                pressure=pressure,
                guest_bundles=guest_bundles,
            ),
            encoding="utf-8",
        )

        graspa_stdout_log_path = pressure_run_dir / f"{raspa_backend}.stdout.log"
        graspa_stderr_log_path = pressure_run_dir / f"{raspa_backend}.stderr.log"
        record_execution(pressure_run_dir, [str(graspa_binary)])
        try:
            with graspa_stdout_log_path.open("w", encoding="utf-8") as stdout_handle:
                with graspa_stderr_log_path.open("w", encoding="utf-8") as stderr_handle:
                    graspa_completed = subprocess.run(
                        [str(graspa_binary)],
                        cwd=pressure_run_dir,
                        check=False,
                        stdout=stdout_handle,
                        stderr=stderr_handle,
                        text=True,
                        timeout=graspa_timeout_seconds,
                    )
        except subprocess.TimeoutExpired as exc:
            raise GraspaExecutionError(
                f"{raspa_display_name} timed out after {graspa_timeout_seconds} seconds while running "
                f"mixture pressure {pressure:g} Pa. See {graspa_stdout_log_path} and {graspa_stderr_log_path}."
            ) from exc

        if graspa_completed.returncode != 0:
            raise GraspaExecutionError(
                f"{raspa_display_name} failed with exit code {graspa_completed.returncode} at mixture pressure {pressure:g} Pa. "
                f"See {graspa_stdout_log_path} and {graspa_stderr_log_path}."
            )

        data_file_paths = _collect_raspa_data_file_paths(pressure_run_dir / "Output")
        if not data_file_paths:
            raise GraspaParseError(
                f"{raspa_display_name} completed without producing any Output/**/*.data files under {pressure_run_dir / 'Output'} "
                f"for mixture pressure {pressure:g} Pa."
            )

        point_result = _parse_mixture_result_files(
            tuple(Path(path) for path in data_file_paths),
            mixture_settings=mixture_settings,
            pressure=pressure,
            pressure_run_dir=pressure_run_dir,
            simulation_input_path=simulation_input_path,
            graspa_stdout_log_path=graspa_stdout_log_path,
            graspa_stderr_log_path=graspa_stderr_log_path,
        )
        if staged_initial_restart_file is not None:
            point_result = replace(
                point_result,
                initial_restart_file_path=str(staged_initial_restart_file),
            )
        point_results.append(point_result)
        warnings.extend(point_result.warnings)

        if len(data_file_paths) > 1:
            warnings.append(
                f"Parsed mixture point {pressure:g} Pa from {len(data_file_paths)} data files."
            )
        if graspa_stderr_log_path.read_text(encoding="utf-8", errors="replace").strip():
            warnings.append(
                f"{raspa_display_name} wrote content to stderr for mixture pressure {pressure:g} Pa; inspect the stderr log if needed."
            )
        if any(
            not math.isfinite(value)
            for component_result in point_result.component_results
            for value in (
                component_result.feed_mol_fraction,
                component_result.adsorbed_mol_fraction,
                component_result.loading_mol_per_kg,
                component_result.loading_mol_per_kg_errorbar,
                component_result.loading_g_per_l,
                component_result.loading_g_per_l_errorbar,
                component_result.heat_of_adsorption_kj_per_mol,
                component_result.heat_of_adsorption_kj_per_mol_errorbar,
            )
        ):
            warnings.append(
                f"One or more parsed mixture component values are non-finite at pressure {pressure:g} Pa; inspect the raw {raspa_display_name} data output."
            )
        if any(
            not math.isfinite(value)
            for selectivity_result in point_result.selectivity_results
            for value in (
                selectivity_result.feed_mol_fraction_ratio,
                selectivity_result.adsorbed_mol_fraction_ratio,
                selectivity_result.selectivity,
            )
        ):
            warnings.append(
                f"One or more parsed mixture selectivity values are non-finite at pressure {pressure:g} Pa; inspect the raw {raspa_display_name} data output."
            )

    if eqeq_result.eqeq_json_output_path is None:
        warnings.append("EQeq did not write the companion JSON output file.")

    component_results_csv_path = mixture_root_dir / "component_results.csv"
    _write_mixture_component_results_csv(point_results, component_results_csv_path)
    warnings.append("Selectivity uncertainty is unavailable: backend marginal loading errors do not provide paired covariance/statistics.")
    selectivity_results_csv_path = mixture_root_dir / "selectivity_results.csv"
    _write_mixture_selectivity_results_csv(point_results, selectivity_results_csv_path)

    report_path = run_dir / "graspa_mixture_report.json"
    result = GraspaMixtureResult(
        input_cif=str(input_path),
        output_dir=str(run_dir),
        eqeq_binary=eqeq_result.eqeq_binary,
        graspa_binary=str(graspa_binary),
        raspa_backend=raspa_backend,
        eqeq_run_dir=eqeq_result.output_dir,
        mixture_root_dir=str(mixture_root_dir),
        eqeq_input_cif=eqeq_result.eqeq_input_cif,
        eqeq_charged_cif=eqeq_result.eqeq_charged_cif,
        mixture_framework_cif=str(mixture_framework_cif_path),
        eqeq_json_output_path=eqeq_result.eqeq_json_output_path,
        eqeq_stdout_log_path=eqeq_result.eqeq_stdout_log_path,
        eqeq_stderr_log_path=eqeq_result.eqeq_stderr_log_path,
        component_results_csv_path=str(component_results_csv_path),
        selectivity_results_csv_path=str(selectivity_results_csv_path),
        report_path=str(report_path),
        unit_cells=unit_cells,
        eqeq_settings=eqeq_settings,
        mixture_settings=mixture_settings,
        point_results=tuple(point_results),
        warnings=tuple(warnings),
        framework_forcefield_metadata=resolve_forcefield_metadata(mixture_settings.forcefield).to_dict(),
    )
    atomic_write_text(report_path, json.dumps(result.to_dict(), indent=2, allow_nan=False))
    finish_attempt(run_dir)
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
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(f"{field_name} must be positive when provided.")
    return float(value)


def _validate_eqeq_settings(settings: EqeqChargeSettings) -> None:
    validate_numbers(settings)
    if settings.method not in _SUPPORTED_EQEQ_METHODS:
        raise ValueError(
            f"Unsupported EQeq method {settings.method!r}. "
            f"Supported methods: {sorted(_SUPPORTED_EQEQ_METHODS)}."
        )
    if settings.charge_precision < 0:
        raise ValueError("charge_precision must be non-negative.")
    if settings.real_space_cells < 0 or settings.reciprocal_space_cells < 0:
        raise ValueError("real_space_cells and reciprocal_space_cells must be non-negative.")
    if settings.net_charge_tolerance <= 0:
        raise ValueError("net_charge_tolerance must be positive.")
    if settings.lambda_value <= 0:
        raise ValueError("lambda_value must be positive.")
    if settings.eta <= 0.0:
        raise ValueError("eta must be positive.")


def _normalize_graspa_forcefield_name(forcefield: str) -> str:
    try:
        return resolve_forcefield_metadata(forcefield).family
    except ForceFieldMetadataError:
        return forcefield.strip().lower().replace("-", "_")


def _normalize_raspa_backend_name(backend: str) -> str:
    return backend.strip().lower().replace("-", "").replace("_", "")


def _raspa_backend_display_name(backend: str) -> str:
    normalized = _normalize_raspa_backend_name(backend)
    if normalized == "graspa":
        return "gRASPA"
    if normalized == "raspa2":
        return "RASPA2"
    return backend


def _validate_raspa_backend(backend: str) -> str:
    normalized = _normalize_raspa_backend_name(backend)
    if normalized not in SUPPORTED_RASPA_BACKENDS:
        raise ValueError(
            f"Unsupported RASPA backend {backend!r}. Expected one of: {', '.join(SUPPORTED_RASPA_BACKENDS)}."
        )
    return normalized


def _load_graspa_guest_bundles(paths: Sequence[str]) -> tuple[GuestBundle, ...]:
    try:
        bundles = load_guest_bundles(paths)
    except GuestBundleError as exc:
        raise GraspaConfigurationError(str(exc)) from exc
    reserved = {component.casefold(): component for component in AVAILABLE_WIDOM_COMPONENTS}
    for bundle in bundles:
        for name in bundle.names():
            reserved_name = reserved.get(name.casefold())
            if reserved_name is not None:
                raise GraspaConfigurationError(
                    f"Guest bundle component name or alias {name!r} conflicts with packaged component {reserved_name!r}."
                )
    return bundles


def _canonicalize_widom_settings(
    settings: GraspaWidomSettings,
    guest_bundles: Sequence[GuestBundle],
) -> GraspaWidomSettings:
    return replace(
        settings,
        components=tuple(_canonical_component_name(component, guest_bundles) for component in settings.components),
    )


def _canonicalize_isotherm_settings(
    settings: GraspaIsothermSettings,
    guest_bundles: Sequence[GuestBundle],
) -> GraspaIsothermSettings:
    return replace(settings, component=_canonical_component_name(settings.component, guest_bundles))


def _canonicalize_mixture_settings(
    settings: GraspaMixtureSettings,
    guest_bundles: Sequence[GuestBundle],
) -> GraspaMixtureSettings:
    return replace(
        settings,
        components=tuple(
            replace(component, component=_canonical_component_name(component.component, guest_bundles))
            for component in settings.components
        ),
    )


def _canonical_component_name(component: str, guest_bundles: Sequence[GuestBundle]) -> str:
    candidate = component.strip()
    packaged = {name.casefold(): name for name in AVAILABLE_WIDOM_COMPONENTS}.get(candidate.casefold())
    if packaged is not None:
        return packaged
    bundle = _guest_bundle_catalog(guest_bundles).get(candidate.casefold())
    if bundle is not None:
        return bundle.name
    return candidate


def _guest_bundle_catalog(guest_bundles: Sequence[GuestBundle]) -> dict[str, GuestBundle]:
    catalog: dict[str, GuestBundle] = {}
    for bundle in guest_bundles:
        for name in bundle.names():
            catalog[name.casefold()] = bundle
    return catalog


def _supported_graspa_component_names(guest_bundles: Sequence[GuestBundle]) -> tuple[str, ...]:
    return (*AVAILABLE_WIDOM_COMPONENTS, *(bundle.name for bundle in guest_bundles))


def _is_supported_graspa_component(component: str, guest_bundles: Sequence[GuestBundle]) -> bool:
    return component in AVAILABLE_WIDOM_COMPONENTS or any(component == bundle.name for bundle in guest_bundles)


def _guest_forcefield_metadata_dict(
    component: str,
    guest_bundles: Sequence[GuestBundle],
) -> dict[str, object]:
    packaged = packaged_guest_forcefield_catalog().get(component)
    if packaged is not None:
        return {"source": "packaged", **packaged.to_dict()}
    for bundle in guest_bundles:
        if component == bundle.name:
            return {
                "source": "guest_bundle",
                "source_path": bundle.source_path,
                "name": bundle.name,
                "parameter_family": bundle.parameter_family,
                "parameter_source": bundle.parameter_source,
                "rotatable": bundle.rotatable,
                "vdw_treatment": bundle.vdw_treatment,
                "tail_corrections": bundle.tail_corrections,
                "mixing_rule": bundle.mixing_rule,
            }
    raise GraspaConfigurationError(f"No guest force-field metadata is registered for component {component!r}.")


def _guest_nonbonded_convention(
    components: Sequence[str],
    guest_bundles: Sequence[GuestBundle],
) -> tuple[str, bool, str]:
    by_convention: dict[tuple[str, bool, str], list[str]] = {}
    for component in components:
        metadata = _guest_forcefield_metadata_dict(component, guest_bundles)
        convention = (
            str(metadata.get("vdw_treatment", "truncated")),
            bool(metadata.get("tail_corrections", False)),
            str(metadata.get("mixing_rule", "lorentz_berthelot")),
        )
        by_convention.setdefault(convention, []).append(component)
    if len(by_convention) != 1:
        details = "; ".join(
            f"{', '.join(names)}: vdw={vdw}, tail_corrections={tail}, mixing_rule={mixing}"
            for (vdw, tail, mixing), names in by_convention.items()
        )
        raise ValueError(
            "Selected guest components require incompatible global RASPA nonbonded conventions. "
            f"Choose models with one shared convention; {details}."
        )
    return next(iter(by_convention))


def _component_is_rotatable(component: str, guest_bundles: Sequence[GuestBundle]) -> bool:
    for bundle in guest_bundles:
        if component == bundle.name:
            return bundle.rotatable
    packaged = packaged_guest_forcefield_catalog().get(component)
    if packaged is not None:
        return packaged.rotatable
    return True


def _validate_eos_component(component: str, guest_bundles: Sequence[GuestBundle]) -> None:
    bundle = next((b for b in guest_bundles if b.name == component), None)
    definition = (bundle.raspa.molecule_definition_text if bundle else
                  (_widom_template_dir() / f"{component}.def").read_text())
    rows = [line.strip() for line in definition.splitlines() if line.strip() and not line.lstrip().startswith("#")]
    try:
        tc, pc, omega = (float(row) for row in rows[:3])
    except (ValueError, TypeError) as exc:
        raise ValueError(f"PR-EOS requires critical constants for {component}; supply an explicit fugacity coefficient.") from exc
    if not all(map(math.isfinite, (tc, pc, omega))) or tc <= 0 or pc <= 0:
        raise ValueError(f"Invalid critical constants for {component}; PR-EOS cannot be used.")


def _validate_calculation_contract(settings, guest_bundles: Sequence[GuestBundle]) -> None:
    backend = _validate_raspa_backend(settings.backend)
    if not re.fullmatch(r"[A-Za-z0-9_-]+", settings.framework_name):
        raise ValueError("framework_name must contain only letters, digits, underscores, or hyphens.")
    if settings.number_of_simulations != 1 or not settings.single_simulation:
        raise ValueError("Result parsing supports exactly one simulation per calculation.")
    if settings.max_step_per_cycle <= 0:
        raise ValueError("max_step_per_cycle must be positive.")
    if settings.random_seed < 0:
        raise ValueError("random_seed must be non-negative.")
    if settings.number_of_blocks > settings.production_cycles:
        raise ValueError("number_of_blocks cannot exceed production_cycles.")
    movies_every = getattr(settings, "movies_every", None)
    if movies_every is not None:
        if movies_every <= 0:
            raise ValueError("movies_every must be positive when provided.")
        if backend != "graspa":
            raise ValueError("movies_every is currently supported only for the gRASPA backend.")
    if isinstance(settings, GraspaIsothermSettings):
        components = (settings.component,)
        moves = (settings,)
    elif isinstance(settings, GraspaMixtureSettings):
        components = tuple(item.component for item in settings.components)
        moves = settings.components
    else:
        components = settings.components
        moves = ()
    if len(set(components)) != len(components):
        raise ValueError("Duplicate components are not supported.")
    for item in moves:
        if sum((item.translation_probability, item.rotation_probability,
                item.reinsertion_probability, item.swap_probability)) <= 0:
            raise ValueError("At least one sampling move weight must be positive.")
        if item.swap_probability <= 0:
            raise ValueError("Grand-canonical adsorption requires positive swap_probability.")
    for item in moves:
        if isinstance(item.fugacity_coefficient, str):
            _validate_eos_component(item.component, guest_bundles)
    if backend == "graspa":
        if "H2_DREIDING" in components:
            raise GraspaConfigurationError("H2_DREIDING requires Feynman-Hibbs interactions, unsupported by the verified upstream gRASPA parser. Select a verified supporting backend explicitly.")
        for bundle in guest_bundles:
            if bundle.name in components and any(row.split()[1] != "lennard-jones" for row in bundle.raspa.mixing_rule_rows):
                raise GraspaConfigurationError("gRASPA guest bundles currently require verified lennard-jones interactions.")


def _validate_widom_settings(
    settings: GraspaWidomSettings,
    *,
    guest_bundles: Sequence[GuestBundle] = (),
) -> None:
    validate_numbers(settings)
    _validate_calculation_contract(settings, guest_bundles)
    forcefield = _normalize_graspa_forcefield_name(settings.forcefield)
    if forcefield not in _SUPPORTED_GRASPA_FORCEFIELDS:
        raise ValueError(
            f"Unsupported gRASPA forcefield {settings.forcefield!r}. "
            f"Expected one of: {', '.join(_SUPPORTED_GRASPA_FORCEFIELDS)}."
        )
    if not settings.components:
        raise ValueError("At least one Widom component must be configured.")
    for component in settings.components:
        if not _is_supported_graspa_component(component, guest_bundles):
            supported = ", ".join(_supported_graspa_component_names(guest_bundles))
            raise ValueError(f"Unsupported gRASPA component {component!r}. Supported components: {supported}.")
    _guest_nonbonded_convention(settings.components, guest_bundles)
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
    if settings.widom_moves_per_component is not None and settings.widom_moves_per_component <= 0:
        raise ValueError("widom_moves_per_component must be positive when provided.")
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


def _validate_isotherm_settings(
    settings: GraspaIsothermSettings,
    *,
    guest_bundles: Sequence[GuestBundle] = (),
) -> None:
    validate_numbers(settings)
    _validate_calculation_contract(settings, guest_bundles)
    forcefield = _normalize_graspa_forcefield_name(settings.forcefield)
    if forcefield not in _SUPPORTED_GRASPA_FORCEFIELDS:
        raise ValueError(
            f"Unsupported gRASPA forcefield {settings.forcefield!r}. "
            f"Expected one of: {', '.join(_SUPPORTED_GRASPA_FORCEFIELDS)}."
        )
    if not _is_supported_graspa_component(settings.component, guest_bundles):
        supported = ", ".join(_supported_graspa_component_names(guest_bundles))
        raise ValueError(f"Unsupported gRASPA component {settings.component!r}. Supported components: {supported}.")
    _guest_nonbonded_convention((settings.component,), guest_bundles)
    if not settings.pressures:
        raise ValueError("At least one pressure point must be configured.")
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
    if settings.temperature <= 0.0:
        raise ValueError("temperature must be positive.")
    if any(pressure < 0.0 for pressure in settings.pressures):
        raise ValueError("pressure points must be non-negative.")
    if settings.cutoff_vdw <= 0.0 or settings.cutoff_coulomb <= 0.0:
        raise ValueError("cutoff_vdw and cutoff_coulomb must be positive.")
    if settings.ewald_precision <= 0.0:
        raise ValueError("ewald_precision must be positive.")
    for name, value in (
        ("translation_probability", settings.translation_probability),
        ("rotation_probability", settings.rotation_probability),
        ("reinsertion_probability", settings.reinsertion_probability),
        ("swap_probability", settings.swap_probability),
    ):
        if value < 0.0:
            raise ValueError(f"{name} must be non-negative.")
    if settings.create_number_of_molecules < 0:
        raise ValueError("create_number_of_molecules must be non-negative.")
    if not settings.save_output_to_file:
        raise ValueError("save_output_to_file must remain enabled for result parsing in the current workflow.")
    _validate_fugacity_coefficient(settings.fugacity_coefficient)


def _validate_mixture_settings(
    settings: GraspaMixtureSettings,
    *,
    guest_bundles: Sequence[GuestBundle] = (),
) -> None:
    validate_numbers(settings)
    _validate_calculation_contract(settings, guest_bundles)
    forcefield = _normalize_graspa_forcefield_name(settings.forcefield)
    if forcefield not in _SUPPORTED_GRASPA_FORCEFIELDS:
        raise ValueError(
            f"Unsupported gRASPA forcefield {settings.forcefield!r}. "
            f"Expected one of: {', '.join(_SUPPORTED_GRASPA_FORCEFIELDS)}."
        )
    if len(settings.components) < 2:
        raise ValueError("At least two adsorbate components must be configured for mixture adsorption.")
    seen_components: set[str] = set()
    mol_fraction_total = 0.0
    for component in settings.components:
        if not _is_supported_graspa_component(component.component, guest_bundles):
            supported = ", ".join(_supported_graspa_component_names(guest_bundles))
            raise ValueError(
                f"Unsupported gRASPA component {component.component!r}. Supported components: {supported}."
            )
        if component.component in seen_components:
            raise ValueError(f"Duplicate gRASPA mixture component {component.component!r} is not allowed.")
        seen_components.add(component.component)
        if component.mol_fraction <= 0.0:
            raise ValueError("mixture component mol_fraction values must be positive.")
        mol_fraction_total += component.mol_fraction
        _validate_fugacity_coefficient(component.fugacity_coefficient)
        for name, value in (
            ("translation_probability", component.translation_probability),
            ("rotation_probability", component.rotation_probability),
            ("reinsertion_probability", component.reinsertion_probability),
            ("identity_change_probability", component.identity_change_probability),
            ("swap_probability", component.swap_probability),
        ):
            if value < 0.0:
                raise ValueError(f"{name} must be non-negative.")
        if component.create_number_of_molecules < 0:
            raise ValueError("create_number_of_molecules must be non-negative.")

    if mol_fraction_total <= 0.0:
        raise ValueError("The total mixture mol_fraction must be positive.")
    _guest_nonbonded_convention(
        tuple(component.component for component in settings.components),
        guest_bundles,
    )
    if not settings.pressures:
        raise ValueError("At least one pressure point must be configured.")
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
    if settings.temperature <= 0.0:
        raise ValueError("temperature must be positive.")
    if any(pressure < 0.0 for pressure in settings.pressures):
        raise ValueError("pressure points must be non-negative.")
    if settings.cutoff_vdw <= 0.0 or settings.cutoff_coulomb <= 0.0:
        raise ValueError("cutoff_vdw and cutoff_coulomb must be positive.")
    if settings.ewald_precision <= 0.0:
        raise ValueError("ewald_precision must be positive.")
    if not settings.save_output_to_file:
        raise ValueError("save_output_to_file must remain enabled for result parsing in the current workflow.")


def _resolve_output_dir(input_path: Path, output_dir: str | Path | None, *, suffix: str) -> Path:
    if output_dir is not None:
        return Path(output_dir).expanduser().resolve()
    return input_path.parent / f"{input_path.stem}{suffix}"


def _stage_graspa_initial_restart_file(
    pressure_run_dir: Path,
    initial_restart_file: str | Path,
) -> Path:
    source_path = Path(initial_restart_file).expanduser().resolve()
    if not source_path.is_file():
        raise GraspaConfigurationError(f"Initial gRASPA restart file does not exist: {source_path}")
    staged_path = pressure_run_dir / "RestartInitial" / "System_0" / "restartfile"
    staged_path.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source_path, staged_path)
    return staged_path


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


def _copy_widom_template_assets(
    destination_dir: Path,
    components: Sequence[str],
    *,
    forcefield: str,
    guest_bundles: Sequence[GuestBundle] = (),
) -> None:
    template_dir = _widom_template_dir()
    force_field_source_path = template_dir / "force_field.def"
    if not force_field_source_path.is_file():
        raise GraspaConfigurationError(f"Packaged Widom asset is missing: {force_field_source_path}")
    shutil.copy2(force_field_source_path, destination_dir / "force_field.def")

    selected_bundles = _guest_bundles_for_components(components, guest_bundles)
    selected_bundle_names = {bundle.name for bundle in selected_bundles}
    packaged_catalog = packaged_guest_forcefield_catalog()
    selected_packaged = tuple(
        packaged_catalog[component]
        for component in components
        if component not in selected_bundle_names
    )
    for component in components:
        if component in selected_bundle_names:
            continue
        source_path = template_dir / f"{component}.def"
        if not source_path.is_file():
            raise GraspaConfigurationError(f"Packaged Widom asset is missing: {source_path}")
        shutil.copy2(source_path, destination_dir / f"{component}.def")
    for bundle in selected_bundles:
        (destination_dir / f"{bundle.name}.def").write_text(
            bundle.raspa.molecule_definition_text.rstrip() + "\n",
            encoding="utf-8",
        )

    pseudo_atoms_source_path = template_dir / "pseudo_atoms.def"
    if not pseudo_atoms_source_path.is_file():
        raise GraspaConfigurationError(f"Packaged Widom asset is missing: {pseudo_atoms_source_path}")
    pseudo_atom_rows = (
        tuple(row for metadata in selected_packaged for row in metadata.pseudo_atom_rows)
        + tuple(row for bundle in selected_bundles for row in bundle.raspa.pseudo_atom_rows)
    )
    (destination_dir / "pseudo_atoms.def").write_text(
        _append_pseudo_atom_rows(
            pseudo_atoms_source_path.read_text(encoding="utf-8"),
            pseudo_atom_rows,
        ),
        encoding="utf-8",
    )
    guest_mixing_rule_rows = (
        tuple(row for metadata in selected_packaged for row in metadata.mixing_rule_rows)
        + tuple(row for bundle in selected_bundles for row in bundle.raspa.mixing_rule_rows)
    )
    vdw_treatment, tail_corrections, mixing_rule = _guest_nonbonded_convention(
        components,
        guest_bundles,
    )
    (destination_dir / "force_field_mixing_rules.def").write_text(
        _render_graspa_force_field_mixing_rules(
            forcefield,
            guest_mixing_rule_rows=guest_mixing_rule_rows,
            vdw_treatment=vdw_treatment,
            tail_corrections=tail_corrections,
            mixing_rule=mixing_rule,
        ),
        encoding="utf-8",
    )
    framework_forcefield_metadata = resolve_forcefield_metadata(forcefield)
    (destination_dir / "framework_forcefield_metadata.json").write_text(
        json.dumps(framework_forcefield_metadata.to_dict(), indent=2) + "\n",
        encoding="utf-8",
    )
    (destination_dir / "guest_forcefield_metadata.json").write_text(
        json.dumps(
            {
                "schema_version": PACKAGED_GUEST_FORCEFIELD_SCHEMA_VERSION,
                "framework_forcefield": framework_forcefield_metadata.family,
                "framework_forcefield_id": framework_forcefield_metadata.id,
                "guests": [
                    _guest_forcefield_metadata_dict(component, guest_bundles) for component in components
                ],
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )


def _guest_bundles_for_components(
    components: Sequence[str],
    guest_bundles: Sequence[GuestBundle],
) -> tuple[GuestBundle, ...]:
    selected: list[GuestBundle] = []
    seen: set[str] = set()
    for component in components:
        for bundle in guest_bundles:
            if component == bundle.name and bundle.name not in seen:
                selected.append(bundle)
                seen.add(bundle.name)
                break
    return tuple(selected)


def _append_pseudo_atom_rows(base_text: str, rows: Sequence[str]) -> str:
    if not rows:
        return base_text.rstrip() + "\n"
    lines = base_text.splitlines()
    count_index = _first_integer_line_index(lines)
    if count_index is None:
        raise GraspaConfigurationError("Packaged pseudo_atoms.def is missing the pseudo-atom count line.")
    try:
        base_count = int(lines[count_index].strip())
    except ValueError as exc:
        raise GraspaConfigurationError("Packaged pseudo_atoms.def has an invalid pseudo-atom count line.") from exc

    existing_types = {
        row_type
        for line in lines
        if (row_type := _row_type_token(line)) is not None
    }
    appended_types: set[str] = set()
    for row in rows:
        row_type = _row_type_token(row)
        if row_type is None:
            raise GraspaConfigurationError(f"Invalid guest pseudo-atom row: {row!r}")
        if row_type in existing_types or row_type in appended_types:
            raise GraspaConfigurationError(f"Duplicate gRASPA pseudo-atom type in guest bundles: {row_type}")
        appended_types.add(row_type)

    lines[count_index] = str(base_count + len(rows))
    # NOTE: gRASPA parses pseudo_atoms.def with a strict line counter (one atom
    # per line, count must match the header); an empty separator line here
    # shifts every appended row and causes an out-of-bounds crash in gRASPA.
    lines.extend(row.strip() for row in rows)
    return "\n".join(lines).rstrip() + "\n"


def _first_integer_line_index(lines: Sequence[str]) -> int | None:
    for index, line in enumerate(lines):
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        try:
            int(stripped)
        except ValueError:
            continue
        return index
    return None


def _row_type_token(row: str) -> str | None:
    stripped = row.strip()
    if not stripped or stripped.startswith("#"):
        return None
    if len(stripped.split()) == 1:
        try:
            int(stripped)
        except ValueError:
            pass
        else:
            return None
    return stripped.split()[0]


def _bundled_uff_parameter_file() -> Path:
    try:
        return packaged_forcefield_artifact_path("uff", "parameter_table")
    except ForceFieldMetadataError as exc:
        raise GraspaConfigurationError(str(exc)) from exc


@dataclass(frozen=True)
class _UffFrameworkParameters:
    x1: float
    d1: float


@lru_cache(maxsize=1)
def _load_uff_framework_parameters() -> dict[str, _UffFrameworkParameters]:
    parameters: dict[str, _UffFrameworkParameters] = {}
    for raw_line in _bundled_uff_parameter_file().read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw_line.strip()
        if not line.startswith("param "):
            continue
        parts = line.split()
        if len(parts) < 13:
            continue
        atom_type = parts[1]
        if atom_type not in _UFF_FRAMEWORK_TYPE_BY_ELEMENT.values():
            continue
        parameters[atom_type] = _UffFrameworkParameters(x1=float(parts[4]), d1=float(parts[5]))
    missing = [atom_type for atom_type in _UFF_FRAMEWORK_TYPE_BY_ELEMENT.values() if atom_type not in parameters]
    if missing:
        raise GraspaConfigurationError(
            "Packaged UFF.prm is missing one or more framework representative atom types needed for gRASPA: "
            + ", ".join(sorted(missing))
        )
    return parameters


def _graspa_framework_rows_for_forcefield(forcefield: str) -> list[str]:
    normalized = _normalize_graspa_forcefield_name(forcefield)
    rows: list[str] = []
    if normalized == "dreiding":
        for element in _GRASPA_FRAMEWORK_ELEMENT_ORDER:
            atom_type = DREIDING_FRAMEWORK_TYPE_BY_ELEMENT[element]
            parameters = DREIDING_PARAMETERS[atom_type]
            epsilon_kelvin = parameters.d0 * _KCAL_PER_MOL_TO_KELVIN
            sigma = parameters.r0 / _UFF_RMIN_TO_SIGMA_FACTOR
            rows.append(
                f"{element:<14} lennard-jones {epsilon_kelvin:10.4f} {sigma:10.5f}      "
                f"// {DREIDING_REFERENCE_SOURCE}"
            )
        return rows
    if normalized == "uff":
        uff_parameters = _load_uff_framework_parameters()
        for element in _GRASPA_FRAMEWORK_ELEMENT_ORDER:
            atom_type = _UFF_FRAMEWORK_TYPE_BY_ELEMENT[element]
            parameters = uff_parameters[atom_type]
            epsilon_kelvin = parameters.d1 * _KCAL_PER_MOL_TO_KELVIN
            sigma = parameters.x1 / _UFF_RMIN_TO_SIGMA_FACTOR
            rows.append(
                f"{element:<14} lennard-jones {epsilon_kelvin:10.4f} {sigma:10.5f}      "
                "// UFF from bundled Open Babel UFF.prm"
            )
        return rows
    raise GraspaConfigurationError(
        f"Unsupported gRASPA forcefield backend requested: {forcefield!r}"
    )


def _render_graspa_force_field_mixing_rules(
    forcefield: str,
    *,
    guest_mixing_rule_rows: Sequence[str] = (),
    vdw_treatment: str = "truncated",
    tail_corrections: bool = False,
    mixing_rule: str = "lorentz_berthelot",
) -> str:
    framework_rows = _graspa_framework_rows_for_forcefield(forcefield)
    base_rows = [*framework_rows, *_GRASPA_SHARED_MIXING_RULE_LINES]
    _validate_guest_mixing_rule_rows(base_rows, guest_mixing_rule_rows)
    interaction_rows = [*base_rows, *(row.strip() for row in guest_mixing_rule_rows)]
    return "\n".join(
        [
            "# general rule for shifted vs truncated",
            vdw_treatment,
            "# general rule tailcorrections",
            _bool_token(tail_corrections),
            "# number of defined interactions",
            str(len(interaction_rows)),
            "# type,        interaction,  epsilon [K],  sigma [A] IMPORTANT: define shortest matches first, so that more specific ones overwrites these",
            *interaction_rows,
            "",
            "# general mixing rule for Lennard-Jones",
            mixing_rule.replace("_", "-").title(),
            "",
        ]
    )


def _validate_guest_mixing_rule_rows(base_rows: Sequence[str], guest_rows: Sequence[str]) -> None:
    existing_types = {
        row_type
        for row in base_rows
        if (row_type := _row_type_token(row)) is not None
    }
    appended_types: set[str] = set()
    for row in guest_rows:
        row_type = _row_type_token(row)
        if row_type is None:
            raise GraspaConfigurationError(f"Invalid guest force-field mixing-rule row: {row!r}")
        if row_type in existing_types or row_type in appended_types:
            raise GraspaConfigurationError(f"Duplicate gRASPA mixing-rule type in guest bundles: {row_type}")
        appended_types.add(row_type)


def _write_graspa_framework_cif(source: Path, destination: Path) -> None:
    """Copy the charged CIF for gRASPA with pseudo-atom-compatible labels.

    Charge validation restores the caller's original atom labels in the
    user-facing charged CIF, but gRASPA matches `_atom_site_label` against
    pseudo_atoms.def entries, and the framework force-field assets name
    pseudo atoms by element symbol. Rewrite the labels to their type symbols
    in the simulation input copy only.
    """
    import gemmi
    source = Path(source)
    destination = Path(destination)
    header_comments: list[str] = []
    with source.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                header_comments.append(line.rstrip("\n"))
            elif line.strip():
                break
    document = gemmi.cif.read_file(str(source))
    block = document.sole_block()
    labels = block.find_loop("_atom_site_label")
    type_symbols = block.find_loop("_atom_site_type_symbol")
    if not labels or not type_symbols or len(labels) != len(type_symbols):
        raise GraspaConfigurationError(
            "Charged CIF must provide matching _atom_site_label and "
            "_atom_site_type_symbol loops."
        )
    for index in range(len(labels)):
        labels[index] = gemmi.cif.quote(gemmi.cif.as_string(type_symbols[index]))
    document.write_file(str(destination))
    if header_comments:
        body = destination.read_text(encoding="utf-8")
        destination.write_text(
            "\n".join(header_comments) + "\n" + body, encoding="utf-8"
        )


def _compute_unit_cells_from_cif(cif_path: Path, *, cutoff: float) -> tuple[int, int, int]:
    import gemmi
    block = gemmi.cif.read_file(str(cif_path)).sole_block()
    tags = ("length_a", "length_b", "length_c", "angle_alpha", "angle_beta", "angle_gamma")
    parameters = [gemmi.cif.as_number(block.find_value("_cell_" + tag) or "?") for tag in tags]
    if not all(math.isfinite(v) and v > 0 for v in parameters):
        raise GraspaParseError("Supercell calculation requires all six finite, positive cell parameters.")
    if not all(0 < angle < 180 for angle in parameters[3:]):
        raise GraspaParseError("Cell angles must be strictly between 0 and 180 degrees.")
    try:
        widths = face_widths(gemmi.UnitCell(*parameters))
    except ValueError as exc:
        raise GraspaParseError(str(exc)) from exc
    return tuple(_compute_supercell(width, cutoff) for width in widths)


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
    if not all(math.isfinite(x) and x > 0 for x in (cell_length, cutoff)):
        raise ValueError("Cell width and cutoff must be finite and positive.")
    # Conservative ceiling: never round a slightly undersized box down.
    return max(math.ceil(2.0 * cutoff / cell_length), 1)


def _render_widom_simulation_input(
    settings: GraspaWidomSettings,
    *,
    unit_cells: tuple[int, int, int],
    guest_bundles: Sequence[GuestBundle] = (),
) -> str:
    if _normalize_raspa_backend_name(settings.backend) == "raspa2":
        return _render_raspa2_widom_simulation_input(
            settings,
            unit_cells=unit_cells,
            guest_bundles=guest_bundles,
        )

    lines = [
        f"UseGPUReduction {_bool_token(settings.use_gpu_reduction)}",
        f"UseFastHostRNG {_bool_token(settings.use_fast_host_rng)}",
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


def _render_isotherm_simulation_input(
    settings: GraspaIsothermSettings,
    *,
    unit_cells: tuple[int, int, int],
    pressure: float,
    guest_bundles: Sequence[GuestBundle] = (),
) -> str:
    if pressure == 0:
        settings = replace(settings, fugacity_coefficient=1.0)
    if _normalize_raspa_backend_name(settings.backend) == "raspa2":
        return _render_raspa2_isotherm_simulation_input(
            settings,
            unit_cells=unit_cells,
            pressure=pressure,
            guest_bundles=guest_bundles,
        )

    lines = [
        f"UseGPUReduction {_bool_token(settings.use_gpu_reduction)}",
        f"UseFastHostRNG {_bool_token(settings.use_fast_host_rng)}",
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
        f"Pressure     {_format_simulation_number(pressure)}",
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
        f"Component 0 MoleculeName             {settings.component}",
        "            IdealGasRosenbluthWeight 1.0",
        f"            FugacityCoefficient      {_render_fugacity_coefficient(settings.fugacity_coefficient)}",
        f"            TranslationProbability   {_format_simulation_number(settings.translation_probability)}",
    ]
    if settings.movies_every is not None:
        # gRASPA scans simulation.input line-by-line for MoviesEvery
        # (read_movies_stats_print), so keyword placement is order-independent.
        lines.extend([f"MoviesEvery {settings.movies_every}", ""])
    if _component_is_rotatable(settings.component, guest_bundles) and settings.rotation_probability > 0.0:
        lines.append(
            f"            RotationProbability      {_format_simulation_number(settings.rotation_probability)}"
        )
    lines.extend(
        [
            f"            ReinsertionProbability   {_format_simulation_number(settings.reinsertion_probability)}",
            f"            SwapProbability          {_format_simulation_number(settings.swap_probability)}",
            f"            CreateNumberOfMolecules  {settings.create_number_of_molecules}",
        ]
    )
    return "\n".join(lines).rstrip() + "\n"


def _render_mixture_simulation_input(
    settings: GraspaMixtureSettings,
    *,
    unit_cells: tuple[int, int, int],
    pressure: float,
    guest_bundles: Sequence[GuestBundle] = (),
) -> str:
    if pressure == 0:
        settings = replace(settings, components=tuple(replace(c, fugacity_coefficient=1.0) for c in settings.components))
    if _normalize_raspa_backend_name(settings.backend) == "raspa2":
        return _render_raspa2_mixture_simulation_input(
            settings,
            unit_cells=unit_cells,
            pressure=pressure,
            guest_bundles=guest_bundles,
        )

    lines = [
        f"UseGPUReduction {_bool_token(settings.use_gpu_reduction)}",
        f"UseFastHostRNG {_bool_token(settings.use_fast_host_rng)}",
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
        f"Pressure     {_format_simulation_number(pressure)}",
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
    if settings.movies_every is not None:
        # gRASPA scans simulation.input line-by-line for MoviesEvery
        # (read_movies_stats_print), so keyword placement is order-independent.
        lines.extend([f"MoviesEvery {settings.movies_every}", ""])
    for index, component in enumerate(settings.components):
        lines.extend(
            [
                f"Component {index} MoleculeName             {component.component}",
                "            IdealGasRosenbluthWeight 1.0",
                f"            FugacityCoefficient      {_render_fugacity_coefficient(component.fugacity_coefficient)}",
                f"            MolFraction              {_format_simulation_number(component.mol_fraction)}",
                f"            TranslationProbability   {_format_simulation_number(component.translation_probability)}",
            ]
        )
        if _component_is_rotatable(component.component, guest_bundles) and component.rotation_probability > 0.0:
            lines.append(
                f"            RotationProbability      {_format_simulation_number(component.rotation_probability)}"
            )
        lines.extend(
            [
                f"            ReinsertionProbability   {_format_simulation_number(component.reinsertion_probability)}",
                f"            IdentityChangeProbability {_format_simulation_number(component.identity_change_probability)}",
                f"            SwapProbability          {_format_simulation_number(component.swap_probability)}",
                f"            CreateNumberOfMolecules  {component.create_number_of_molecules}",
                "",
            ]
        )
    return "\n".join(lines).rstrip() + "\n"


def _render_raspa2_widom_simulation_input(
    settings: GraspaWidomSettings,
    *,
    unit_cells: tuple[int, int, int],
    guest_bundles: Sequence[GuestBundle] = (),
) -> str:
    lines = _render_raspa2_common_lines(
        random_seed=settings.random_seed,
        number_of_cycles=settings.production_cycles,
        initialization_cycles=settings.initialization_cycles,
        equilibration_cycles=settings.equilibration_cycles,
        number_of_blocks=settings.number_of_blocks,
        restart_file=settings.restart_file,
        use_charges_from_cif_file=settings.use_charges_from_cif_file,
        framework_name=settings.framework_name,
        input_file_type=settings.input_file_type,
        unit_cells=unit_cells,
        charge_method=settings.charge_method,
        temperature=settings.temperature,
        pressure=None,
        cutoff_vdw=settings.cutoff_vdw,
        cutoff_coulomb=settings.cutoff_coulomb,
        ewald_precision=settings.ewald_precision,
        overlap_criteria=settings.overlap_criteria,
        number_of_trial_positions=settings.number_of_trial_positions,
        number_of_trial_orientations=settings.number_of_trial_orientations,
    )
    lines.append("")
    for index, component in enumerate(settings.components):
        lines.extend(
            [
                f"Component {index} MoleculeName              {component}",
                f"            MoleculeDefinition        {_RASPA2_LOCAL_MOLECULE_DEFINITION}",
                "            IdealGasRosenbluthWeight  1.0",
                "            WidomProbability          1.0",
                "            CreateNumberOfMolecules   0",
                "",
            ]
        )
    return "\n".join(lines).rstrip() + "\n"


def _render_raspa2_isotherm_simulation_input(
    settings: GraspaIsothermSettings,
    *,
    unit_cells: tuple[int, int, int],
    pressure: float,
    guest_bundles: Sequence[GuestBundle] = (),
) -> str:
    lines = _render_raspa2_common_lines(
        random_seed=settings.random_seed,
        number_of_cycles=settings.production_cycles,
        initialization_cycles=settings.initialization_cycles,
        equilibration_cycles=settings.equilibration_cycles,
        number_of_blocks=settings.number_of_blocks,
        restart_file=settings.restart_file,
        use_charges_from_cif_file=settings.use_charges_from_cif_file,
        framework_name=settings.framework_name,
        input_file_type=settings.input_file_type,
        unit_cells=unit_cells,
        charge_method=settings.charge_method,
        temperature=settings.temperature,
        pressure=pressure,
        cutoff_vdw=settings.cutoff_vdw,
        cutoff_coulomb=settings.cutoff_coulomb,
        ewald_precision=settings.ewald_precision,
        overlap_criteria=settings.overlap_criteria,
        number_of_trial_positions=settings.number_of_trial_positions,
        number_of_trial_orientations=settings.number_of_trial_orientations,
    )
    lines.append("")
    lines.extend(
        [
            f"Component 0 MoleculeName             {settings.component}",
            f"            MoleculeDefinition       {_RASPA2_LOCAL_MOLECULE_DEFINITION}",
        ]
    )
    rendered_fugacity = _render_raspa2_fugacity_coefficient(settings.fugacity_coefficient)
    if rendered_fugacity is not None:
        lines.append(f"            FugacityCoefficient      {rendered_fugacity}")
    lines.append(f"            TranslationProbability   {_format_simulation_number(settings.translation_probability)}")
    if _component_is_rotatable(settings.component, guest_bundles) and settings.rotation_probability > 0.0:
        lines.append(
            f"            RotationProbability      {_format_simulation_number(settings.rotation_probability)}"
        )
    lines.extend(
        [
            f"            ReinsertionProbability   {_format_simulation_number(settings.reinsertion_probability)}",
            f"            SwapProbability          {_format_simulation_number(settings.swap_probability)}",
            f"            CreateNumberOfMolecules  {settings.create_number_of_molecules}",
        ]
    )
    return "\n".join(lines).rstrip() + "\n"


def _render_raspa2_mixture_simulation_input(
    settings: GraspaMixtureSettings,
    *,
    unit_cells: tuple[int, int, int],
    pressure: float,
    guest_bundles: Sequence[GuestBundle] = (),
) -> str:
    lines = _render_raspa2_common_lines(
        random_seed=settings.random_seed,
        number_of_cycles=settings.production_cycles,
        initialization_cycles=settings.initialization_cycles,
        equilibration_cycles=settings.equilibration_cycles,
        number_of_blocks=settings.number_of_blocks,
        restart_file=settings.restart_file,
        use_charges_from_cif_file=settings.use_charges_from_cif_file,
        framework_name=settings.framework_name,
        input_file_type=settings.input_file_type,
        unit_cells=unit_cells,
        charge_method=settings.charge_method,
        temperature=settings.temperature,
        pressure=pressure,
        cutoff_vdw=settings.cutoff_vdw,
        cutoff_coulomb=settings.cutoff_coulomb,
        ewald_precision=settings.ewald_precision,
        overlap_criteria=settings.overlap_criteria,
        number_of_trial_positions=settings.number_of_trial_positions,
        number_of_trial_orientations=settings.number_of_trial_orientations,
    )
    lines.append("")
    for index, component in enumerate(settings.components):
        lines.extend(
            [
                f"Component {index} MoleculeName             {component.component}",
                f"            MoleculeDefinition       {_RASPA2_LOCAL_MOLECULE_DEFINITION}",
            ]
        )
        rendered_fugacity = _render_raspa2_fugacity_coefficient(component.fugacity_coefficient)
        if rendered_fugacity is not None:
            lines.append(f"            FugacityCoefficient      {rendered_fugacity}")
        lines.extend(
            [
                f"            MolFraction              {_format_simulation_number(component.mol_fraction)}",
                f"            TranslationProbability   {_format_simulation_number(component.translation_probability)}",
            ]
        )
        if _component_is_rotatable(component.component, guest_bundles) and component.rotation_probability > 0.0:
            lines.append(
                f"            RotationProbability      {_format_simulation_number(component.rotation_probability)}"
            )
        lines.extend(
            [
                f"            ReinsertionProbability   {_format_simulation_number(component.reinsertion_probability)}",
                f"            IdentityChangeProbability {_format_simulation_number(component.identity_change_probability)}",
                f"            SwapProbability          {_format_simulation_number(component.swap_probability)}",
                f"            CreateNumberOfMolecules  {component.create_number_of_molecules}",
                "",
            ]
        )
    return "\n".join(lines).rstrip() + "\n"


def _render_raspa2_common_lines(
    *,
    random_seed: int,
    number_of_cycles: int,
    initialization_cycles: int,
    equilibration_cycles: int,
    number_of_blocks: int,
    restart_file: bool,
    use_charges_from_cif_file: bool,
    framework_name: str,
    input_file_type: str,
    unit_cells: tuple[int, int, int],
    charge_method: str,
    temperature: float,
    pressure: float | None,
    cutoff_vdw: float,
    cutoff_coulomb: float,
    ewald_precision: float,
    overlap_criteria: float,
    number_of_trial_positions: int,
    number_of_trial_orientations: int,
) -> list[str]:
    lines = [
        "SimulationType                MonteCarlo",
        f"RandomSeed                    {random_seed}",
        f"NumberOfCycles                {number_of_cycles}",
        f"NumberOfInitializationCycles  {initialization_cycles}",
    ]
    if equilibration_cycles > 0:
        lines.append(f"NumberOfEquilibrationCycles   {equilibration_cycles}")
    lines.extend(
        [
            f"PrintEvery                    {_raspa2_print_every(number_of_cycles, number_of_blocks)}",
            f"RestartFile                   {_bool_token(restart_file)}",
            "",
            f"Forcefield                    {_RASPA2_LOCAL_FORCEFIELD_NAME}",
            f"UseChargesFromCIFFile         {_bool_token(use_charges_from_cif_file)}",
            "",
            "Framework 0",
            f"FrameworkName {framework_name}",
            f"InputFileType {input_file_type}",
            f"UnitCells {unit_cells[0]} {unit_cells[1]} {unit_cells[2]}",
            f"ExternalTemperature {_format_simulation_number(temperature)}",
        ]
    )
    if pressure is not None:
        lines.append(f"ExternalPressure {_format_simulation_number(pressure)}")
    lines.extend(
        [
            "",
            f"ChargeMethod {charge_method}",
            f"CutOffVDW {_format_simulation_number(cutoff_vdw)}",
            f"CutOffCoulomb {_format_simulation_number(cutoff_coulomb)}",
            f"EwaldPrecision {_format_simulation_number(ewald_precision)}",
            f"EnergyOverlapCriteria {_format_simulation_number(overlap_criteria)}",
            "",
            f"NumberOfTrialPositions {number_of_trial_positions}",
            f"NumberOfTrialPositionsWidom {number_of_trial_positions}",
            f"NumberOfTrialPositionsSwap {number_of_trial_positions}",
            f"NumberOfTrialPositionsReinsertion {number_of_trial_positions}",
            f"NumberOfTrialPositionsForTheFirstBead {number_of_trial_orientations}",
        ]
    )
    return lines


def _raspa2_print_every(number_of_cycles: int, number_of_blocks: int) -> int:
    if number_of_blocks <= 0:
        return max(number_of_cycles, 1)
    return max(number_of_cycles // number_of_blocks, 1)


def _bool_token(value: bool) -> str:
    return "yes" if value else "no"


def _format_simulation_number(value: float) -> str:
    return f"{value:g}"


def _validate_fugacity_coefficient(value: float | str) -> None:
    if isinstance(value, str):
        if value.strip().casefold() != "pr-eos":
            raise ValueError("fugacity_coefficient must be a float or the literal 'PR-EOS'.")
        return
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError("fugacity_coefficient must be positive when provided as a float.")


def _render_fugacity_coefficient(value: float | str) -> str:
    _validate_fugacity_coefficient(value)
    if isinstance(value, str):
        return "PR-EOS"
    return _format_simulation_number(value)


def _render_raspa2_fugacity_coefficient(value: float | str) -> str | None:
    _validate_fugacity_coefficient(value)
    if isinstance(value, str):
        return None
    return _format_simulation_number(value)


def _format_pressure_run_dir_name(index: int, pressure: float) -> str:
    raw = _format_simulation_number(pressure)
    sanitized = re.sub(r"[^0-9A-Za-z]+", "_", raw).strip("_")
    if sanitized == "":
        sanitized = "0"
    return f"pressure_{index + 1:04d}_{sanitized}Pa"


def _collect_raspa_data_file_paths(output_dir: Path) -> tuple[str, ...]:
    paths = tuple(str(path) for path in sorted(output_dir.rglob("*.data")))
    if len(paths) > 1:
        raise GraspaParseError("Multiple simulation output files are ambiguous; one simulation per attempt is required.")
    return paths


def _json_safe_float(value: float) -> float | None:
    if math.isfinite(value):
        return value
    return None


def _parse_widom_result_files(
    data_paths: Sequence[Path],
    *,
    diagnostics: list[str] | None = None,
) -> tuple[GraspaWidomComponentResult, ...]:
    results: list[GraspaWidomComponentResult] = []
    for path in data_paths:
        content = path.read_text(encoding="utf-8", errors="replace")
        schemas = _detect_output_backend_schemas(content)
        graspa_results: list[GraspaWidomComponentResult] = []
        for match in _WIDOM_RESULT_RE.finditer(content):
            graspa_results.append(
                GraspaWidomComponentResult(
                    component=match.group(1),
                    widom_energy=float(match.group(2)),
                    widom_energy_errorbar=float(match.group(3)),
                    henry=float(match.group(4)),
                    henry_errorbar=float(match.group(5)),
                    source_data_file=str(path),
                    source_schema="graspa",
                )
            )
        raspa2_results = _parse_raspa2_widom_result_file(path, content)
        results.extend(
            _merge_component_results_across_schemas(
                path,
                schemas=schemas,
                graspa_results=graspa_results,
                raspa2_results=raspa2_results,
                diagnostics=diagnostics,
            )
        )
    return tuple(results)


def _merge_component_results_across_schemas(
    path: Path,
    *,
    schemas: Sequence[str],
    graspa_results: Sequence[GraspaWidomComponentResult],
    raspa2_results: Sequence[GraspaWidomComponentResult],
    diagnostics: list[str] | None,
) -> list[GraspaWidomComponentResult]:
    """Merge per-schema Widom parses of one file without losing components.

    A file with markers from both backend layouts is parsed under both
    schemas: disjoint component coverage merges completely (with a recorded
    diagnostic), identical duplicate blocks are deduplicated, and conflicting
    duplicate blocks are an explicit ambiguity error rather than a silent
    choice.
    """
    if not graspa_results:
        if raspa2_results:
            return list(raspa2_results)
        return []
    if not raspa2_results:
        if "raspa2" in schemas and diagnostics is not None:
            diagnostics.append(
                f"Mixed-schema markers in {path}: RASPA2-style markers are present alongside "
                f"{len(graspa_results)} parsed gRASPA component block(s), but no RASPA2 component "
                "summaries were parseable; if the file combines backends, some RASPA2 components "
                "may be unrecoverable."
            )
        return list(graspa_results)

    graspa_components = [result.component for result in graspa_results]
    raspa2_components = [result.component for result in raspa2_results]
    if len(graspa_components) != len(set(graspa_components)) or len(raspa2_components) != len(
        set(raspa2_components)
    ):
        # Within-schema duplicate component blocks: keep every row verbatim so
        # the workflow's strict exact-component check sees the duplication
        # instead of the merge hiding it.
        return [*graspa_results, *raspa2_results]

    merged: dict[str, GraspaWidomComponentResult] = {}
    for result in (*graspa_results, *raspa2_results):
        existing = merged.get(result.component)
        if existing is None:
            merged[result.component] = result
            continue
        conflicting = _value_tuples_conflict(
            (existing.widom_energy, existing.widom_energy_errorbar),
            (result.widom_energy, result.widom_energy_errorbar),
        ) or _value_tuples_conflict(
            (existing.henry, existing.henry_errorbar),
            (result.henry, result.henry_errorbar),
        )
        if conflicting:
            raise GraspaParseError(
                f"Ambiguous mixed-schema output file {path}: component {result.component!r} appears in "
                "both gRASPA and RASPA2 summary blocks with conflicting values; the authoritative "
                "backend schema cannot be determined from content, so no component is dropped or "
                "picked silently. Provide single-backend output per file."
            )
        if diagnostics is not None:
            diagnostics.append(
                f"Mixed-schema file {path}: component {result.component!r} appears in both gRASPA "
                "and RASPA2 summary blocks with identical values; keeping one copy."
            )
    if set(graspa_result.component for graspa_result in graspa_results).isdisjoint(
        raspa2_result.component for raspa2_result in raspa2_results
    ) and diagnostics is not None:
        diagnostics.append(
            f"Mixed-schema file {path}: parsed {len(graspa_results)} component(s) from gRASPA "
            f"blocks and {len(raspa2_results)} component(s) from RASPA2 blocks; all blocks were "
            "recovered and merged."
        )
    return list(merged.values())


def _parse_raspa2_widom_result_file(path: Path, content: str) -> list[GraspaWidomComponentResult]:
    widom_by_component = {
        match.group("component"): (
            _parse_float_token(match.group("value")),
            _parse_float_token(match.group("error")),
        )
        for match in _RASPA2_WIDOM_EXCESS_RE.finditer(content)
    }
    results: list[GraspaWidomComponentResult] = []
    for match in _RASPA2_HENRY_RE.finditer(content):
        component = match.group("component")
        widom_energy, widom_energy_errorbar = widom_by_component.get(component, (math.nan, math.nan))
        results.append(
            GraspaWidomComponentResult(
                component=component,
                widom_energy=widom_energy,
                widom_energy_errorbar=widom_energy_errorbar,
                henry=_parse_float_token(match.group("value")),
                henry_errorbar=_parse_float_token(match.group("error")),
                source_data_file=str(path),
                source_schema="raspa2",
            )
        )
    return results


def _parse_isotherm_result_files(
    data_paths: Sequence[Path],
    *,
    component: str,
    pressure: float,
    pressure_run_dir: Path,
    simulation_input_path: Path,
    graspa_stdout_log_path: Path,
    graspa_stderr_log_path: Path,
) -> GraspaIsothermPointResult:
    for path in data_paths:
        parsed = _parse_isotherm_result_file(path, component=component)
        if parsed is None:
            continue
        parsed_values, source_schema, point_warnings = parsed
        return GraspaIsothermPointResult(
            component=component,
            pressure=pressure,
            pressure_run_dir=str(pressure_run_dir),
            simulation_input_path=str(simulation_input_path),
            graspa_stdout_log_path=str(graspa_stdout_log_path),
            graspa_stderr_log_path=str(graspa_stderr_log_path),
            data_file_paths=tuple(str(item) for item in data_paths),
            source_data_file=str(path),
            source_schema=source_schema,
            warnings=tuple(point_warnings),
            loading_mol_per_kg=parsed_values["loading_mol_per_kg"][0],
            loading_mol_per_kg_errorbar=parsed_values["loading_mol_per_kg"][1],
            loading_g_per_l=parsed_values["loading_g_per_l"][0],
            loading_g_per_l_errorbar=parsed_values["loading_g_per_l"][1],
            heat_of_adsorption_kj_per_mol=parsed_values["heat_of_adsorption_kj_per_mol"][0],
            heat_of_adsorption_kj_per_mol_errorbar=parsed_values["heat_of_adsorption_kj_per_mol"][1],
        )
    raise GraspaParseError(
        f"No adsorption isotherm summaries for component {component!r} could be parsed from "
        f"{len(data_paths)} data file(s) at pressure {pressure:g} Pa."
    )


def _parse_mixture_result_files(
    data_paths: Sequence[Path],
    *,
    mixture_settings: GraspaMixtureSettings,
    pressure: float,
    pressure_run_dir: Path,
    simulation_input_path: Path,
    graspa_stdout_log_path: Path,
    graspa_stderr_log_path: Path,
) -> GraspaMixturePointResult:
    expected_components = tuple(component.component for component in mixture_settings.components)
    normalized_feed_mol_fractions = _normalized_mixture_mol_fractions(mixture_settings.components)
    section_titles = {
        "loading_mol_per_kg": "BLOCK AVERAGES (LOADING: mol/kg)",
        "loading_g_per_l": "BLOCK AVERAGES (LOADING: g/L)",
        "heat_of_adsorption_kj_per_mol": "BLOCK AVERAGES (HEAT OF ADSORPTION: kJ/mol)",
    }

    unresolved_file_details: list[str] = []
    for path in data_paths:
        content = path.read_text(encoding="utf-8", errors="replace")
        schemas = _detect_output_backend_schemas(content)
        parsed_sections = _parse_component_average_sections(content, section_titles=section_titles)
        graspa_values_by_component = {
            component_name: values
            for component_name in expected_components
            if (
                values := _resolve_component_summary_values(parsed_sections, component=component_name)
            )
            is not None
        }
        raspa2_loading_by_component = _parse_raspa2_component_loading_summaries(content)
        raspa2_values_by_component = {
            component_name: _raspa2_mixture_component_values_from_loading(
                raspa2_loading_by_component[component_name]
            )
            for component_name in expected_components
            if component_name in raspa2_loading_by_component
        }

        conflicting_components = sorted(
            component_name
            for component_name in set(graspa_values_by_component) & set(raspa2_values_by_component)
            if _value_tuples_conflict(
                graspa_values_by_component[component_name]["loading_mol_per_kg"],
                raspa2_values_by_component[component_name]["loading_mol_per_kg"],
            )
        )
        if conflicting_components:
            raise GraspaParseError(
                f"Ambiguous mixed-schema output file {path}: component(s) {conflicting_components!r} "
                "have conflicting gRASPA and RASPA2 loading summaries; the authoritative backend "
                "schema cannot be determined from content, so no component is dropped or picked "
                "silently. Provide single-backend output per file."
            )

        component_source_schemas: dict[str, str] = {}
        resolved_components: list[tuple[str, dict[str, tuple[float, float]]]] = []
        for component_name in expected_components:
            if component_name in graspa_values_by_component:
                resolved_components.append((component_name, graspa_values_by_component[component_name]))
                component_source_schemas[component_name] = "graspa"
            elif component_name in raspa2_values_by_component:
                resolved_components.append((component_name, raspa2_values_by_component[component_name]))
                component_source_schemas[component_name] = "raspa2"
        missing_components = sorted(set(expected_components) - set(component_source_schemas))
        if missing_components:
            unresolved_file_details.append(
                f"{path}: gRASPA block sections provide {sorted(graspa_values_by_component)!r}, "
                f"RASPA2 loading rows provide {sorted(raspa2_values_by_component)!r}; "
                f"missing {missing_components!r}."
            )
            continue

        point_warnings: list[str] = []
        if len(set(component_source_schemas.values())) > 1:
            point_warnings.append(
                f"Mixed-schema file {path}: component coverage spans both backend layouts "
                f"({dict(sorted(component_source_schemas.items()))}); every requested component "
                "was recovered, but g/L and heat-of-adsorption values are unavailable for the "
                "RASPA2-layout components. This usually means concatenated or renamed backend "
                "output; verify the run provenance."
            )
        elif raspa2_values_by_component and "graspa" in schemas:
            point_warnings.append(
                f"Mixed-schema markers in {path}: components were parsed from the RASPA2 layout, "
                "but gRASPA block-average markers are also present without complete parseable "
                "sections; if the file combines backends, some gRASPA data may be unrecoverable."
            )

        loading_by_component = {
            component_name: values["loading_mol_per_kg"][0]
            for component_name, values in resolved_components
        }
        nonfinite_loading_components = tuple(
            component_name
            for component_name, loading in loading_by_component.items()
            if not math.isfinite(loading)
        )
        if nonfinite_loading_components:
            # A complete adsorbed composition requires every required component
            # loading; normalizing over only the finite subset would fabricate
            # full fractions (e.g. 1.0) for the measured components, so the
            # composition stays unavailable and per-component loadings remain
            # marginal measurements.
            total_loading_mol_per_kg = math.nan
            point_warnings.append(
                f"Adsorbed mixture composition is unavailable at pressure {pressure:g} Pa: "
                f"non-finite loading for component(s) {list(nonfinite_loading_components)!r} means the "
                "normalization denominator is incomplete; per-component loadings are retained as "
                "marginal measurements."
            )
        else:
            total_loading_mol_per_kg = sum(loading_by_component.values())
            if total_loading_mol_per_kg <= 0.0:
                point_warnings.append(
                    f"Adsorbed mixture composition is unavailable at pressure {pressure:g} Pa: "
                    "the total adsorbed loading is zero (measured zero loadings, not missing data), "
                    "so mole fractions are undefined."
                )

        component_results = tuple(
            GraspaMixturePointComponentResult(
                component=component_name,
                feed_mol_fraction=normalized_feed_mol_fractions[component_name],
                adsorbed_mol_fraction=_compute_adsorbed_mol_fraction(
                    loading_by_component[component_name],
                    total_loading_mol_per_kg,
                ),
                loading_mol_per_kg=values["loading_mol_per_kg"][0],
                loading_mol_per_kg_errorbar=values["loading_mol_per_kg"][1],
                loading_g_per_l=values["loading_g_per_l"][0],
                loading_g_per_l_errorbar=values["loading_g_per_l"][1],
                heat_of_adsorption_kj_per_mol=values["heat_of_adsorption_kj_per_mol"][0],
                heat_of_adsorption_kj_per_mol_errorbar=values["heat_of_adsorption_kj_per_mol"][1],
            )
            for component_name, values in resolved_components
        )
        selectivity_results = tuple(
            GraspaMixtureSelectivityResult(
                numerator_component=numerator_component,
                denominator_component=denominator_component,
                feed_mol_fraction_ratio=_safe_ratio(
                    normalized_feed_mol_fractions[numerator_component],
                    normalized_feed_mol_fractions[denominator_component],
                ),
                adsorbed_mol_fraction_ratio=_safe_ratio(
                    loading_by_component[numerator_component],
                    loading_by_component[denominator_component],
                ),
                selectivity=_compute_selectivity(
                    numerator_loading=loading_by_component[numerator_component],
                    denominator_loading=loading_by_component[denominator_component],
                    numerator_feed_mol_fraction=normalized_feed_mol_fractions[numerator_component],
                    denominator_feed_mol_fraction=normalized_feed_mol_fractions[denominator_component],
                ),
                selectivity_errorbar=math.nan,
            )
            for numerator_component in expected_components
            for denominator_component in expected_components
            if numerator_component != denominator_component
        )

        return GraspaMixturePointResult(
            pressure=pressure,
            pressure_run_dir=str(pressure_run_dir),
            simulation_input_path=str(simulation_input_path),
            graspa_stdout_log_path=str(graspa_stdout_log_path),
            graspa_stderr_log_path=str(graspa_stderr_log_path),
            data_file_paths=tuple(str(item) for item in data_paths),
            source_data_file=str(path),
            component_results=component_results,
            selectivity_results=selectivity_results,
            warnings=tuple(point_warnings),
            component_source_schemas=component_source_schemas,
        )

    detail_suffix = ""
    if unresolved_file_details:
        detail_suffix = " Per-file parse coverage: " + " ".join(unresolved_file_details)
    raise GraspaParseError(
        f"No mixture adsorption summaries for components {list(expected_components)!r} could be parsed from "
        f"{len(data_paths)} data file(s) at pressure {pressure:g} Pa.{detail_suffix}"
    )


def _parse_isotherm_result_file(
    path: Path,
    *,
    component: str,
) -> tuple[dict[str, tuple[float, float]], str, list[str]] | None:
    """Parse one isotherm data file under both backend output schemas.

    Returns (values, source_schema, diagnostics) or None when neither schema
    yields the requested component. Mixed-schema files are resolved
    explicitly: agreeing blocks merge with a diagnostic, conflicting blocks
    raise instead of silently picking one backend.
    """
    content = path.read_text(encoding="utf-8", errors="replace")
    section_titles = {
        "loading_mol_per_kg": "BLOCK AVERAGES (LOADING: mol/kg)",
        "loading_g_per_l": "BLOCK AVERAGES (LOADING: g/L)",
        "heat_of_adsorption_kj_per_mol": "BLOCK AVERAGES (HEAT OF ADSORPTION: kJ/mol)",
    }
    schemas = _detect_output_backend_schemas(content)
    parsed_sections = _parse_component_average_sections(content, section_titles=section_titles)
    graspa_values = _resolve_component_summary_values(parsed_sections, component=component)
    raspa2_values = _parse_raspa2_single_component_adsorption_result(content, component=component)

    if graspa_values is None:
        if raspa2_values is None:
            return None
        warnings = _unparsed_schema_marker_warnings(path, schemas=schemas, parsed_schema="raspa2")
        return raspa2_values, "raspa2", warnings
    if raspa2_values is None:
        warnings = _unparsed_schema_marker_warnings(path, schemas=schemas, parsed_schema="graspa")
        return graspa_values, "graspa", warnings

    if _value_tuples_conflict(
        graspa_values["loading_mol_per_kg"],
        raspa2_values["loading_mol_per_kg"],
    ):
        raise GraspaParseError(
            f"Ambiguous mixed-schema output file {path}: component {component!r} has conflicting "
            "gRASPA and RASPA2 loading summaries; the authoritative backend schema cannot be "
            "determined from content, so no value is picked silently. Provide single-backend "
            "output per file."
        )
    warnings = [
        f"Mixed-schema file {path}: component {component!r} appears in both gRASPA and RASPA2 "
        "summary blocks with the same loading; using the gRASPA block-average sections, which "
        "also carry g/L and heat-of-adsorption values."
    ]
    return graspa_values, "graspa", warnings


def _unparsed_schema_marker_warnings(
    path: Path,
    *,
    schemas: Sequence[str],
    parsed_schema: str,
) -> list[str]:
    """Warn when a file carries markers of a schema that produced no parse."""
    other_schemas = [schema for schema in schemas if schema != parsed_schema]
    if not other_schemas:
        return []
    return [
        f"Mixed-schema markers in {path}: parsed component summaries from the {parsed_schema} "
        f"layout, but {', '.join(other_schemas)} markers are also present without parseable "
        "component summaries; if the file combines backends, some components may be unrecoverable."
    ]


def _parse_raspa2_single_component_adsorption_result(
    content: str,
    *,
    component: str,
) -> dict[str, tuple[float, float]] | None:
    loading_by_component = _parse_raspa2_component_loading_summaries(content, default_component=component)
    loading = loading_by_component.get(component)
    if loading is None:
        return None
    enthalpy_match = _RASPA2_ENTHALPY_KJMOL_RE.search(content)
    if enthalpy_match is None:
        enthalpy = (math.nan, math.nan)
    else:
        enthalpy = (
            _parse_float_token(enthalpy_match.group(1)),
            _parse_float_token(enthalpy_match.group(2)),
        )
    return {
        "loading_mol_per_kg": loading,
        "loading_g_per_l": (math.nan, math.nan),
        "heat_of_adsorption_kj_per_mol": enthalpy,
    }


def _raspa2_mixture_component_values_from_loading(
    loading: tuple[float, float],
) -> dict[str, tuple[float, float]]:
    return {
        "loading_mol_per_kg": loading,
        "loading_g_per_l": (math.nan, math.nan),
        "heat_of_adsorption_kj_per_mol": (math.nan, math.nan),
    }


def _parse_raspa2_component_loading_summaries(
    content: str,
    *,
    default_component: str | None = None,
) -> dict[str, tuple[float, float]]:
    loading_by_component: dict[str, tuple[float, float]] = {}
    active_component = default_component
    in_number_of_molecules = False
    for raw_line in content.splitlines():
        line = raw_line.strip()
        if line == "Number of molecules:":
            in_number_of_molecules = True
            active_component = default_component
            continue
        if not in_number_of_molecules:
            continue
        if line.startswith("Average Widom") or line.startswith("Enthalpy of adsorption:"):
            break
        component_match = _RASPA2_COMPONENT_HEADER_RE.match(line)
        if component_match is not None:
            active_component = component_match.group(1)
            continue
        loading_match = _RASPA2_AVERAGE_LOADING_MOL_KG_RE.search(line)
        if loading_match is None or active_component is None:
            continue
        loading_by_component[active_component] = (
            _parse_float_token(loading_match.group(1)),
            _parse_float_token(loading_match.group(2)),
        )
    return loading_by_component


def _parse_component_average_sections(
    content: str,
    *,
    section_titles: dict[str, str],
) -> dict[str, dict[str, tuple[float, float]]]:
    component_values: dict[str, dict[str, tuple[float, float]]] = {
        key: {} for key in section_titles
    }

    active_section: str | None = None
    active_component: str | None = None
    for raw_line in content.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        matched_section = None
        for section_key, title_fragment in section_titles.items():
            if title_fragment in line:
                matched_section = section_key
                break
        if matched_section is not None:
            active_section = matched_section
            active_component = None
            continue
        if "BLOCK AVERAGES (" in line and matched_section is None:
            active_section = None
            active_component = None
            continue
        if active_section is None:
            continue
        if line.startswith("COMPONENT ["):
            match = re.match(r"COMPONENT \[\d+\] \((.+)\)", line)
            active_component = match.group(1) if match is not None else None
            continue
        if active_component is None:
            continue
        if line.startswith("Overall: Average:"):
            match = re.match(
                rf"Overall: Average:\s*({_FLOAT_TOKEN_PATTERN}),\s*ErrorBar:\s*({_FLOAT_TOKEN_PATTERN})",
                line,
                re.IGNORECASE,
            )
            if match is None:
                continue
            component_values[active_section][active_component] = (
                _parse_float_token(match.group(1)),
                _parse_float_token(match.group(2)),
            )
            active_component = None
    return component_values


def _resolve_component_summary_values(
    parsed_sections: dict[str, dict[str, tuple[float, float]]],
    *,
    component: str,
) -> dict[str, tuple[float, float]] | None:
    parsed = {
        section: values.get(component)
        for section, values in parsed_sections.items()
    }
    if parsed.get("loading_mol_per_kg") is None:
        return None

    resolved: dict[str, tuple[float, float]] = {}
    for section, value in parsed.items():
        if value is None:
            resolved[section] = (math.nan, math.nan)
        else:
            resolved[section] = value
    return resolved


def _normalized_mixture_mol_fractions(
    components: Sequence[GraspaMixtureComponentSettings],
) -> dict[str, float]:
    total = sum(component.mol_fraction for component in components)
    if total <= 0.0:
        raise ValueError("The total mixture mol_fraction must be positive.")
    return {
        component.component: component.mol_fraction / total
        for component in components
    }


def _compute_adsorbed_mol_fraction(loading: float, total_loading: float) -> float:
    return _safe_ratio(loading, total_loading)


def _compute_selectivity(
    *,
    numerator_loading: float,
    denominator_loading: float,
    numerator_feed_mol_fraction: float,
    denominator_feed_mol_fraction: float,
) -> float:
    adsorbed_ratio = _safe_ratio(numerator_loading, denominator_loading)
    feed_ratio = _safe_ratio(numerator_feed_mol_fraction, denominator_feed_mol_fraction)
    return _safe_ratio(adsorbed_ratio, feed_ratio)


def _safe_ratio(numerator: float, denominator: float) -> float:
    if not math.isfinite(numerator) or not math.isfinite(denominator):
        return math.nan
    if denominator <= 0.0:
        return math.nan
    if numerator < 0.0:
        return math.nan
    return numerator / denominator


def _parse_float_token(raw_value: str) -> float:
    lowered = raw_value.strip().lower()
    if lowered in {"nan", "+nan", "-nan"}:
        return math.nan
    if lowered in {"inf", "+inf"}:
        return math.inf
    if lowered == "-inf":
        return -math.inf
    return float(raw_value)


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


def _write_isotherm_results_csv(
    point_results: Sequence[GraspaIsothermPointResult],
    output_path: Path,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "component",
                "pressure_pa",
                "loading_mol_per_kg",
                "loading_mol_per_kg_errorbar",
                "loading_mmol_per_g",
                "loading_mmol_per_g_errorbar",
                "loading_g_per_l",
                "loading_g_per_l_errorbar",
                "heat_of_adsorption_kj_per_mol",
                "heat_of_adsorption_kj_per_mol_errorbar",
                "source_data_file",
            ],
        )
        writer.writeheader()
        for result in point_results:
            writer.writerow(
                {
                    "component": result.component,
                    "pressure_pa": result.pressure,
                    "loading_mol_per_kg": result.loading_mol_per_kg,
                    "loading_mol_per_kg_errorbar": result.loading_mol_per_kg_errorbar,
                    "loading_mmol_per_g": result.loading_mol_per_kg,
                    "loading_mmol_per_g_errorbar": result.loading_mol_per_kg_errorbar,
                    "loading_g_per_l": result.loading_g_per_l,
                    "loading_g_per_l_errorbar": result.loading_g_per_l_errorbar,
                    "heat_of_adsorption_kj_per_mol": result.heat_of_adsorption_kj_per_mol,
                    "heat_of_adsorption_kj_per_mol_errorbar": result.heat_of_adsorption_kj_per_mol_errorbar,
                    "source_data_file": result.source_data_file,
                }
            )


def _write_mixture_component_results_csv(
    point_results: Sequence[GraspaMixturePointResult],
    output_path: Path,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "pressure_pa",
                "component",
                "feed_mol_fraction",
                "adsorbed_mol_fraction",
                "loading_mol_per_kg",
                "loading_mol_per_kg_errorbar",
                "loading_mmol_per_g",
                "loading_mmol_per_g_errorbar",
                "loading_g_per_l",
                "loading_g_per_l_errorbar",
                "heat_of_adsorption_kj_per_mol",
                "heat_of_adsorption_kj_per_mol_errorbar",
                "source_data_file",
            ],
        )
        writer.writeheader()
        for point_result in point_results:
            for component_result in point_result.component_results:
                writer.writerow(
                    {
                        "pressure_pa": point_result.pressure,
                        "component": component_result.component,
                        "feed_mol_fraction": component_result.feed_mol_fraction,
                        "adsorbed_mol_fraction": component_result.adsorbed_mol_fraction,
                        "loading_mol_per_kg": component_result.loading_mol_per_kg,
                        "loading_mol_per_kg_errorbar": component_result.loading_mol_per_kg_errorbar,
                        "loading_mmol_per_g": component_result.loading_mol_per_kg,
                        "loading_mmol_per_g_errorbar": component_result.loading_mol_per_kg_errorbar,
                        "loading_g_per_l": component_result.loading_g_per_l,
                        "loading_g_per_l_errorbar": component_result.loading_g_per_l_errorbar,
                        "heat_of_adsorption_kj_per_mol": component_result.heat_of_adsorption_kj_per_mol,
                        "heat_of_adsorption_kj_per_mol_errorbar": (
                            component_result.heat_of_adsorption_kj_per_mol_errorbar
                        ),
                        "source_data_file": point_result.source_data_file,
                    }
                )


def _write_mixture_selectivity_results_csv(
    point_results: Sequence[GraspaMixturePointResult],
    output_path: Path,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "pressure_pa",
                "numerator_component",
                "denominator_component",
                "feed_mol_fraction_ratio",
                "adsorbed_mol_fraction_ratio",
                "selectivity",
                "selectivity_errorbar",
                "source_data_file",
            ],
        )
        writer.writeheader()
        for point_result in point_results:
            for selectivity_result in point_result.selectivity_results:
                writer.writerow(
                    {
                        "pressure_pa": point_result.pressure,
                        "numerator_component": selectivity_result.numerator_component,
                        "denominator_component": selectivity_result.denominator_component,
                        "feed_mol_fraction_ratio": selectivity_result.feed_mol_fraction_ratio,
                        "adsorbed_mol_fraction_ratio": selectivity_result.adsorbed_mol_fraction_ratio,
                        "selectivity": selectivity_result.selectivity,
                        "selectivity_errorbar": selectivity_result.selectivity_errorbar,
                        "source_data_file": point_result.source_data_file,
                    }
                )
