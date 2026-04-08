from __future__ import annotations

import hashlib
import json
import math
import os
import shutil
import subprocess
from dataclasses import asdict, dataclass, replace
from functools import lru_cache
from pathlib import Path
from typing import Sequence

from .cofid import cofid_comment_line, read_cofid_comment_from_cif
from .graspa import (
    EqeqChargeResult,
    EqeqChargeSettings,
    assign_eqeq_charges_to_cif,
)

try:
    import pandas as pd
except ImportError:  # pragma: no cover - exercised in environments without pandas
    pd = None

from .bond_types import (
    bond_order_to_cif_type,
    cif_type_to_bond_order,
    is_aromatic_bond_order,
    normalize_bond_order,
)

try:
    import gemmi
except ImportError:  # pragma: no cover - exercised in environments without gemmi
    gemmi = None

try:
    from openbabel import openbabel as ob
except ImportError:  # pragma: no cover - exercised in environments without Open Babel
    ob = None

try:
    from pymatgen.io.lammps.data import LammpsBox, LammpsData
except ImportError:  # pragma: no cover - exercised in environments without pymatgen
    LammpsBox = None
    LammpsData = None


COFKIT_LMP_ENV_VAR = "COFKIT_LMP_PATH"
DEFAULT_LAMMPS_BINARY: Path | None = None
_SUPPORTED_FORCEFIELDS = ("uff",)
_SUPPORTED_CHARGE_MODELS = ("none", "eqeq")
_UFF_RMIN_TO_SIGMA_FACTOR = 2.0 ** (1.0 / 6.0)
_MIN_MODIFY_LINE_OPTIONS = {"backtrack", "quadratic", "forcezero", "spin_cubic", "spin_none"}
_MIN_MODIFY_NORM_OPTIONS = {"two", "inf", "max"}
_MIN_MODIFY_FIRE_INTEGRATOR_OPTIONS = {"eulerimplicit", "verlet", "leapfrog", "eulerexplicit"}
_BOX_RELAX_MODE_OPTIONS = {"auto", "iso", "aniso", "tri"}
_UNSUPPORTED_BOX_RELAX_MIN_STYLES = {"quickmin", "fire", "hftn", "cg/kk"}
_PINNED_UFF_PARAMETER_SHA256 = "934cb0e2ee1ef2102b2ae8dba74b1e9853b299bd7e26f695fb1fe6e55727bdc5"


class LammpsError(RuntimeError):
    """Base error for LAMMPS wrapper failures."""


class LammpsConfigurationError(LammpsError):
    """Raised when the LAMMPS binary or required Python support cannot be resolved."""


class LammpsInputError(LammpsError):
    """Raised when the input CIF cannot be converted into the supported LAMMPS workflow."""


class LammpsExecutionError(LammpsError):
    """Raised when the LAMMPS subprocess exits unsuccessfully."""


class LammpsParseError(LammpsError):
    """Raised when cofkit cannot parse the generated LAMMPS outputs."""


@dataclass(frozen=True)
class LammpsOptimizationSettings:
    forcefield: str = "uff"
    charge_model: str = "eqeq"
    pair_cutoff: float = 12.0
    coulomb_cutoff: float = 12.0
    ewald_precision: float = 1.0e-6
    position_restraint_force_constant: float = 0.20
    pre_minimization_steps: int = 10000
    pre_minimization_temperature: float = 300.0
    pre_minimization_damping: float = 100.0
    pre_minimization_seed: int = 246813
    pre_minimization_displacement_limit: float = 0.10
    two_stage_protocol: bool = True
    stage2_position_restraint_force_constant: float | None = None
    energy_tolerance: float = 1.0e-6
    force_tolerance: float = 1.0e-6
    max_iterations: int = 200000
    max_evaluations: int = 2000000
    min_style: str = "fire"
    stage2_energy_tolerance: float | None = None
    stage2_force_tolerance: float | None = None
    stage2_max_iterations: int | None = None
    stage2_max_evaluations: int | None = None
    stage2_min_style: str | None = None
    timestep: float | None = None
    min_modify_dmax: float | None = None
    min_modify_line: str | None = None
    min_modify_norm: str | None = None
    min_modify_fire_integrator: str | None = None
    min_modify_fire_tmax: float | None = None
    min_modify_fire_abcfire: bool | None = None
    relax_cell: bool = True
    box_relax_mode: str = "auto"
    box_relax_target_pressure: float = 0.0
    box_relax_vmax: float = 0.001
    box_relax_nreset: int | None = None
    box_relax_min_style: str = "cg"
    box_relax_energy_tolerance: float | None = None
    box_relax_force_tolerance: float | None = None
    box_relax_max_iterations: int | None = None
    box_relax_max_evaluations: int | None = None

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class LammpsOptimizationResult:
    input_cif: str
    optimized_cif: str
    output_dir: str
    lammps_binary: str
    lammps_data_path: str
    lammps_input_script_path: str
    lammps_dump_path: str
    lammps_log_path: str
    stdout_log_path: str
    stderr_log_path: str
    report_path: str
    n_atoms: int
    n_bonds: int
    n_angles: int
    n_dihedrals: int
    n_impropers: int
    n_atom_types: int
    n_bond_types: int
    n_angle_types: int
    n_dihedral_types: int
    n_improper_types: int
    atom_type_symbols: dict[int, str]
    settings: LammpsOptimizationSettings
    forcefield_backend: str
    parameter_sources: dict[str, str]
    charge_model: str
    n_charged_atoms: int
    net_charge: float | None
    eqeq_binary: str | None = None
    eqeq_run_dir: str | None = None
    eqeq_charged_cif: str | None = None
    eqeq_stdout_log_path: str | None = None
    eqeq_stderr_log_path: str | None = None
    eqeq_json_output_path: str | None = None
    warnings: tuple[str, ...] = ()


    def to_dict(self) -> dict[str, object]:
        data = asdict(self)
        data["warnings"] = list(self.warnings)
        return data


@dataclass(frozen=True)
class _AtomChargeExtractionResult:
    charges_by_label: dict[str, float]
    used_atom_order_fallback: bool = False


@dataclass(frozen=True)
class _CifAtomRecord:
    atom_id: int
    label: str
    symbol: str
    occupancy: float
    fractional_position: tuple[float, float, float]
    charge: float | None = None


@dataclass(frozen=True)
class _CifBondRecord:
    bond_id: int
    atom_id_1: int
    atom_id_2: int
    label_1: str
    label_2: str
    symmetry_1: str
    symmetry_2: str
    shift_1: tuple[int, int, int]
    shift_2: tuple[int, int, int]
    equilibrium_distance: float
    bond_order: float | None


@dataclass(frozen=True)
class _CifAngleRecord:
    angle_id: int
    atom_id_1: int
    atom_id_2: int
    atom_id_3: int
    equilibrium_degrees: float


@dataclass(frozen=True)
class _CifDihedralRecord:
    dihedral_id: int
    atom_id_1: int
    atom_id_2: int
    atom_id_3: int
    atom_id_4: int


@dataclass(frozen=True)
class _CifImproperRecord:
    improper_id: int
    atom_id_1: int
    atom_id_2: int
    atom_id_3: int
    atom_id_4: int


@dataclass(frozen=True)
class _ParsedExplicitBondCif:
    source_path: Path
    data_name: str
    cell_parameters: tuple[float, float, float, float, float, float]
    lammps_basis: tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]]
    atoms: tuple[_CifAtomRecord, ...]
    bonds: tuple[_CifBondRecord, ...]
    angles: tuple[_CifAngleRecord, ...]
    dihedrals: tuple[_CifDihedralRecord, ...]
    impropers: tuple[_CifImproperRecord, ...]


@dataclass(frozen=True)
class _PreparedLammpsSystem:
    data_text: str
    has_image_flags: bool
    has_charges: bool
    atom_type_symbols: dict[int, str]
    n_bond_types: int
    n_angle_types: int
    n_dihedral_types: int
    n_improper_types: int
    n_dihedrals: int
    n_impropers: int
    angle_styles: tuple[str, ...]
    dihedral_styles: tuple[str, ...]
    improper_styles: tuple[str, ...]
    forcefield_backend: str
    parameter_sources: dict[str, str]
    n_charged_atoms: int
    net_charge: float | None
    warnings: tuple[str, ...]


@dataclass(frozen=True)
class _OptimizationModel:
    parsed: _ParsedExplicitBondCif
    boundary: str
    warnings: tuple[str, ...] = ()


@dataclass(frozen=True)
class _UffAtomParameters:
    r1: float
    theta0: float
    x1: float
    d1: float
    zeta: float
    z1: float
    vi: float
    uj: float
    xi: float
    hard: float
    radius: float


@dataclass(frozen=True)
class _LammpsMinimizationStage:
    label: str
    min_style: str
    energy_tolerance: float
    force_tolerance: float
    max_iterations: int
    max_evaluations: int
    position_restraint_force_constant: float
    relax_cell: bool = False


@dataclass(frozen=True)
class _LammpsDumpFrame:
    cartesian_positions: dict[int, tuple[float, float, float]]
    lammps_origin: tuple[float, float, float]
    lammps_basis: tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]]
    cell_parameters: tuple[float, float, float, float, float, float]


def resolve_lammps_binary(lmp_path: str | Path | None = None) -> Path:
    raw_value = str(lmp_path) if lmp_path is not None else os.environ.get(COFKIT_LMP_ENV_VAR)
    if raw_value:
        path = Path(raw_value).expanduser()
        if not path.is_file():
            resolved = shutil.which(raw_value)
            if resolved:
                path = Path(resolved)
    elif DEFAULT_LAMMPS_BINARY is not None and DEFAULT_LAMMPS_BINARY.is_file():
        path = DEFAULT_LAMMPS_BINARY
    else:
        raise LammpsConfigurationError(
            f"LAMMPS binary is not configured. Set {COFKIT_LMP_ENV_VAR} or provide --lmp-path."
        )

    if not path.is_file():
        raise LammpsConfigurationError(f"LAMMPS binary does not exist: {path}")
    if not os.access(path, os.X_OK):
        raise LammpsConfigurationError(f"LAMMPS binary is not executable: {path}")
    return path.resolve()


def optimize_cif_with_lammps(
    cif_path: str | Path,
    *,
    output_dir: str | Path | None = None,
    lmp_path: str | Path | None = None,
    eqeq_path: str | Path | None = None,
    settings: LammpsOptimizationSettings | None = None,
    timeout_seconds: float = 300.0,
    eqeq_settings: EqeqChargeSettings | None = None,
    eqeq_timeout_seconds: float | None = 300.0,
) -> LammpsOptimizationResult:
    if timeout_seconds <= 0.0:
        raise ValueError("timeout_seconds must be positive.")
    if eqeq_timeout_seconds is not None and eqeq_timeout_seconds <= 0.0:
        raise ValueError("eqeq_timeout_seconds must be positive when provided.")
    settings = settings or LammpsOptimizationSettings()
    _validate_settings(settings)

    input_path = Path(cif_path).expanduser().resolve()
    if not input_path.is_file():
        raise FileNotFoundError(f"CIF file does not exist: {input_path}")
    input_cofid_comment = read_cofid_comment_from_cif(input_path)

    binary = resolve_lammps_binary(lmp_path)
    run_dir = _resolve_output_dir(input_path, output_dir)
    run_dir.mkdir(parents=True, exist_ok=True)

    parsed = _parse_explicit_bond_cif(input_path)
    eqeq_result: EqeqChargeResult | None = None
    warnings: list[str] = []
    if _normalize_charge_model_name(settings.charge_model) == "eqeq":
        eqeq_result = assign_eqeq_charges_to_cif(
            input_path,
            output_dir=run_dir / "eqeq",
            eqeq_path=eqeq_path,
            settings=eqeq_settings,
            timeout_seconds=eqeq_timeout_seconds,
        )
        extracted_charges = _extract_atom_site_charges_from_cif(
            Path(eqeq_result.eqeq_charged_cif),
            expected_labels=[atom.label for atom in parsed.atoms],
        )
        if extracted_charges.used_atom_order_fallback:
            warnings.append(
                "EQeq output atom labels did not match the input CIF labels, so cofkit mapped charges by atom order. "
                "Verify the EQeq output preserved atom ordering before trusting the assigned charges."
            )
        parsed = _with_atom_charges(
            parsed,
            extracted_charges.charges_by_label,
            source_path=Path(eqeq_result.eqeq_charged_cif),
        )
        warnings.append(
            "EQeq charges were assigned on the input CIF and carried through the LAMMPS optimization. "
            "If you need post-relaxation charges, rerun EQeq on the optimized CIF."
        )
    optimization_model = _build_optimization_model(parsed, settings=settings)
    effective_settings = settings
    prepared = _prepare_lammps_system(optimization_model.parsed, settings=effective_settings)

    data_path = run_dir / "lammps_input.data"
    input_script_path = run_dir / "lammps_minimize.in"
    dump_path = run_dir / "lammps_trajectory.lammpstrj"
    log_path = run_dir / "lammps.log"
    stdout_log_path = run_dir / "lammps.stdout.log"
    stderr_log_path = run_dir / "lammps.stderr.log"
    optimized_cif_path = run_dir / f"{input_path.stem}_lammps_optimized.cif"
    report_path = run_dir / "lammps_report.json"

    data_path.write_text(prepared.data_text, encoding="utf-8")
    input_script_path.write_text(
        _render_lammps_input_script(
            data_path=data_path,
            dump_path=dump_path,
            settings=effective_settings,
            parsed=optimization_model.parsed,
            prepared=prepared,
            boundary=optimization_model.boundary,
        ),
        encoding="utf-8",
    )

    execution_warnings = _run_lammps(
        binary=binary,
        input_script_path=input_script_path,
        log_path=log_path,
        stdout_log_path=stdout_log_path,
        stderr_log_path=stderr_log_path,
        timeout_seconds=timeout_seconds,
    )
    final_frame = _parse_lammps_dump_last_frame(
        dump_path,
        expected_atoms=len(optimization_model.parsed.atoms),
    )
    model_fractional_positions = _cartesian_positions_to_fractional(
        final_frame.lammps_origin,
        final_frame.lammps_basis,
        final_frame.cartesian_positions,
    )
    final_fractional_positions = model_fractional_positions
    optimized_cif_path.write_text(
        _render_optimized_cif(
            parsed,
            final_fractional_positions,
            cell_parameters=final_frame.cell_parameters,
            basis=final_frame.lammps_basis,
            cofid=None if input_cofid_comment is None else input_cofid_comment.cofid,
            cofid_comment_suffix=None if input_cofid_comment is None else input_cofid_comment.suffix,
        ),
        encoding="utf-8",
    )

    parameter_sources = dict(prepared.parameter_sources)
    if eqeq_result is not None:
        parameter_sources["charge_assignment"] = "EQeq charges assigned from the input CIF before LAMMPS export"
    all_warnings = tuple(warnings) + tuple(optimization_model.warnings) + tuple(prepared.warnings) + tuple(execution_warnings)
    resolved_charge_model = (
        "eqeq" if eqeq_result is not None else ("input_cif" if prepared.has_charges else "none")
    )
    result = LammpsOptimizationResult(
        input_cif=str(input_path),
        optimized_cif=str(optimized_cif_path),
        output_dir=str(run_dir),
        lammps_binary=str(binary),
        lammps_data_path=str(data_path),
        lammps_input_script_path=str(input_script_path),
        lammps_dump_path=str(dump_path),
        lammps_log_path=str(log_path),
        stdout_log_path=str(stdout_log_path),
        stderr_log_path=str(stderr_log_path),
        report_path=str(report_path),
        n_atoms=len(optimization_model.parsed.atoms),
        n_bonds=len(optimization_model.parsed.bonds),
        n_angles=len(optimization_model.parsed.angles),
        n_dihedrals=prepared.n_dihedrals,
        n_impropers=prepared.n_impropers,
        n_atom_types=len(prepared.atom_type_symbols),
        n_bond_types=prepared.n_bond_types,
        n_angle_types=prepared.n_angle_types,
        n_dihedral_types=prepared.n_dihedral_types,
        n_improper_types=prepared.n_improper_types,
        atom_type_symbols=prepared.atom_type_symbols,
        settings=effective_settings,
        forcefield_backend=prepared.forcefield_backend,
        parameter_sources=parameter_sources,
        charge_model=resolved_charge_model,
        n_charged_atoms=prepared.n_charged_atoms,
        net_charge=prepared.net_charge,
        eqeq_binary=eqeq_result.eqeq_binary if eqeq_result is not None else None,
        eqeq_run_dir=eqeq_result.output_dir if eqeq_result is not None else None,
        eqeq_charged_cif=eqeq_result.eqeq_charged_cif if eqeq_result is not None else None,
        eqeq_stdout_log_path=eqeq_result.eqeq_stdout_log_path if eqeq_result is not None else None,
        eqeq_stderr_log_path=eqeq_result.eqeq_stderr_log_path if eqeq_result is not None else None,
        eqeq_json_output_path=eqeq_result.eqeq_json_output_path if eqeq_result is not None else None,
        warnings=all_warnings,
    )
    report_path.write_text(json.dumps(result.to_dict(), indent=2), encoding="utf-8")
    return result


def _validate_settings(settings: LammpsOptimizationSettings) -> None:
    forcefield = _normalize_forcefield_name(settings.forcefield)
    if forcefield not in _SUPPORTED_FORCEFIELDS:
        raise ValueError(
            f"Unsupported forcefield {settings.forcefield!r}. Expected one of: {', '.join(_SUPPORTED_FORCEFIELDS)}."
        )
    charge_model = _normalize_charge_model_name(settings.charge_model)
    if charge_model not in _SUPPORTED_CHARGE_MODELS:
        raise ValueError(
            f"Unsupported charge_model {settings.charge_model!r}. Expected one of: {', '.join(_SUPPORTED_CHARGE_MODELS)}."
        )
    if settings.pair_cutoff <= 0.0:
        raise ValueError("pair_cutoff must be positive.")
    if settings.coulomb_cutoff <= 0.0:
        raise ValueError("coulomb_cutoff must be positive.")
    if settings.ewald_precision <= 0.0:
        raise ValueError("ewald_precision must be positive.")
    if settings.position_restraint_force_constant < 0.0:
        raise ValueError("position_restraint_force_constant must be non-negative.")
    if settings.pre_minimization_steps < 0:
        raise ValueError("pre_minimization_steps must be non-negative.")
    if settings.pre_minimization_steps > 0 and settings.pre_minimization_temperature <= 0.0:
        raise ValueError("pre_minimization_temperature must be positive when pre_minimization_steps is enabled.")
    if settings.pre_minimization_steps > 0 and settings.pre_minimization_damping <= 0.0:
        raise ValueError("pre_minimization_damping must be positive when pre_minimization_steps is enabled.")
    if settings.pre_minimization_steps > 0 and settings.pre_minimization_seed <= 0:
        raise ValueError("pre_minimization_seed must be a positive integer when pre_minimization_steps is enabled.")
    if settings.pre_minimization_steps > 0 and settings.pre_minimization_displacement_limit <= 0.0:
        raise ValueError(
            "pre_minimization_displacement_limit must be positive when pre_minimization_steps is enabled."
        )
    if settings.stage2_position_restraint_force_constant is not None and settings.stage2_position_restraint_force_constant < 0.0:
        raise ValueError("stage2_position_restraint_force_constant must be non-negative when provided.")
    if settings.energy_tolerance <= 0.0:
        raise ValueError("energy_tolerance must be positive.")
    if settings.force_tolerance <= 0.0:
        raise ValueError("force_tolerance must be positive.")
    if settings.max_iterations <= 0:
        raise ValueError("max_iterations must be positive.")
    if settings.max_evaluations <= 0:
        raise ValueError("max_evaluations must be positive.")
    if not settings.min_style:
        raise ValueError("min_style must be non-empty.")
    if settings.stage2_energy_tolerance is not None and settings.stage2_energy_tolerance <= 0.0:
        raise ValueError("stage2_energy_tolerance must be positive when provided.")
    if settings.stage2_force_tolerance is not None and settings.stage2_force_tolerance <= 0.0:
        raise ValueError("stage2_force_tolerance must be positive when provided.")
    if settings.stage2_max_iterations is not None and settings.stage2_max_iterations <= 0:
        raise ValueError("stage2_max_iterations must be positive when provided.")
    if settings.stage2_max_evaluations is not None and settings.stage2_max_evaluations <= 0:
        raise ValueError("stage2_max_evaluations must be positive when provided.")
    if settings.stage2_min_style is not None and not settings.stage2_min_style:
        raise ValueError("stage2_min_style must be non-empty when provided.")
    if settings.timestep is not None and settings.timestep <= 0.0:
        raise ValueError("timestep must be positive when provided.")
    if settings.min_modify_dmax is not None and settings.min_modify_dmax <= 0.0:
        raise ValueError("min_modify_dmax must be positive when provided.")
    if settings.min_modify_line is not None and settings.min_modify_line not in _MIN_MODIFY_LINE_OPTIONS:
        raise ValueError(
            "min_modify_line must be one of: " + ", ".join(sorted(_MIN_MODIFY_LINE_OPTIONS))
        )
    if settings.min_modify_norm is not None and settings.min_modify_norm not in _MIN_MODIFY_NORM_OPTIONS:
        raise ValueError(
            "min_modify_norm must be one of: " + ", ".join(sorted(_MIN_MODIFY_NORM_OPTIONS))
        )
    if settings.min_modify_fire_integrator is not None and settings.min_modify_fire_integrator not in _MIN_MODIFY_FIRE_INTEGRATOR_OPTIONS:
        raise ValueError(
            "min_modify_fire_integrator must be one of: "
            + ", ".join(sorted(_MIN_MODIFY_FIRE_INTEGRATOR_OPTIONS))
        )
    if settings.min_modify_fire_tmax is not None and settings.min_modify_fire_tmax <= 0.0:
        raise ValueError("min_modify_fire_tmax must be positive when provided.")
    if settings.box_relax_mode not in _BOX_RELAX_MODE_OPTIONS:
        raise ValueError("box_relax_mode must be one of: " + ", ".join(sorted(_BOX_RELAX_MODE_OPTIONS)))
    if settings.box_relax_vmax <= 0.0:
        raise ValueError("box_relax_vmax must be positive.")
    if settings.box_relax_nreset is not None and settings.box_relax_nreset <= 0:
        raise ValueError("box_relax_nreset must be positive when provided.")
    if not settings.box_relax_min_style:
        raise ValueError("box_relax_min_style must be non-empty.")
    if settings.box_relax_energy_tolerance is not None and settings.box_relax_energy_tolerance <= 0.0:
        raise ValueError("box_relax_energy_tolerance must be positive when provided.")
    if settings.box_relax_force_tolerance is not None and settings.box_relax_force_tolerance <= 0.0:
        raise ValueError("box_relax_force_tolerance must be positive when provided.")
    if settings.box_relax_max_iterations is not None and settings.box_relax_max_iterations <= 0:
        raise ValueError("box_relax_max_iterations must be positive when provided.")
    if settings.box_relax_max_evaluations is not None and settings.box_relax_max_evaluations <= 0:
        raise ValueError("box_relax_max_evaluations must be positive when provided.")
    if settings.relax_cell and _is_box_relax_incompatible_min_style(settings.box_relax_min_style):
        raise ValueError(
            "box_relax_min_style is incompatible with LAMMPS fix box/relax. "
            "Choose a minimizer such as cg or sd instead."
        )

def _build_optimization_model(
    parsed: _ParsedExplicitBondCif,
    *,
    settings: LammpsOptimizationSettings,
) -> _OptimizationModel:
    _ = settings
    return _OptimizationModel(parsed=parsed, boundary="p p p")


def _normalize_forcefield_name(forcefield: str) -> str:
    return forcefield.strip().lower().replace("-", "_")


def _normalize_charge_model_name(charge_model: str) -> str:
    return charge_model.strip().lower().replace("-", "_")


def _normalize_min_style_name(min_style: str) -> str:
    return min_style.strip().lower()


def _is_fire_min_style(min_style: str) -> bool:
    return _normalize_min_style_name(min_style).startswith("fire")


def _is_box_relax_incompatible_min_style(min_style: str) -> bool:
    normalized = _normalize_min_style_name(min_style)
    return normalized in _UNSUPPORTED_BOX_RELAX_MIN_STYLES or normalized.startswith("fire")


def _prepare_lammps_system(
    parsed: _ParsedExplicitBondCif,
    *,
    settings: LammpsOptimizationSettings,
) -> _PreparedLammpsSystem:
    forcefield = _normalize_forcefield_name(settings.forcefield)
    if forcefield == "uff":
        return _prepare_uff_lammps_system(parsed)
    raise LammpsConfigurationError(f"Unsupported forcefield backend requested: {settings.forcefield!r}")


def _prepare_uff_lammps_system(parsed: _ParsedExplicitBondCif) -> _PreparedLammpsSystem:
    _require_uff_support()
    if any(bond.bond_order is None for bond in parsed.bonds):
        raise LammpsInputError(
            "UFF-backed LAMMPS optimization now requires explicit bond orders for every bond. "
            "Populate `_ccdc_geom_bond_type` in the CIF bond loop before running `cofkit calculate lammps-optimize`."
        )
    has_any_charges = any(atom.charge is not None for atom in parsed.atoms)
    has_all_charges = all(atom.charge is not None for atom in parsed.atoms)
    if has_any_charges and not has_all_charges:
        raise LammpsInputError(
            "LAMMPS export requires either a complete set of atom charges or no atom charges at all."
        )

    atom_images, image_conflicts = _compute_unwrapped_atom_images(parsed)
    typing_cartesian_positions = _unwrapped_cartesian_positions(parsed, atom_images)
    wrapped_cartesian_positions = {
        atom.atom_id: _fractional_to_lammps(atom.fractional_position, parsed.lammps_basis) for atom in parsed.atoms
    }
    ob_molecule = _build_openbabel_molecule(parsed, typing_cartesian_positions)
    atom_type_by_atom_id = _assign_openbabel_uff_atom_types(ob_molecule)
    uff_parameters = _load_uff_parameters()
    (
        ordered_type_labels,
        atom_type_ids,
        representative_symbols,
    ) = _ordered_unique_uff_atom_types(parsed.atoms, atom_type_by_atom_id)
    pairij_coeff_rows = _build_uff_pairij_rows(
        ordered_type_labels=ordered_type_labels,
        atom_type_ids=atom_type_ids,
        parameters=uff_parameters,
    )
    bond_topology_rows, bond_coeff_rows, bond_reference = _build_uff_bond_terms(
        parsed=parsed,
        atom_type_by_atom_id=atom_type_by_atom_id,
        parameters=uff_parameters,
    )
    angle_topology_rows, angle_coeff_rows = _build_uff_angle_terms(
        parsed=parsed,
        atom_type_by_atom_id=atom_type_by_atom_id,
        parameters=uff_parameters,
        bond_reference=bond_reference,
    )
    dihedral_topology_rows, dihedral_coeff_rows = _build_uff_dihedral_terms(
        parsed=parsed,
        atom_type_by_atom_id=atom_type_by_atom_id,
        parameters=uff_parameters,
    )
    improper_topology_rows, improper_coeff_rows = _build_uff_improper_terms(
        parsed=parsed,
        atom_type_by_atom_id=atom_type_by_atom_id,
        parameters=uff_parameters,
    )

    a_vec, b_vec, c_vec = parsed.lammps_basis
    box = LammpsBox(
        bounds=[
            [0.0, a_vec[0]],
            [0.0, b_vec[1]],
            [0.0, c_vec[2]],
        ],
        tilt=[b_vec[0], c_vec[0], c_vec[1]],
    )
    masses = pd.DataFrame(
        {"mass": [gemmi.Element(representative_symbols[label]).weight for label in ordered_type_labels]},
        index=range(1, len(ordered_type_labels) + 1),
    )
    atoms = pd.DataFrame(
        [
            (
                {
                    "molecule-ID": 1,
                    "type": atom_type_ids[atom_type_by_atom_id[atom.atom_id]],
                    "q": float(atom.charge),
                    "x": wrapped_cartesian_positions[atom.atom_id][0],
                    "y": wrapped_cartesian_positions[atom.atom_id][1],
                    "z": wrapped_cartesian_positions[atom.atom_id][2],
                }
                if has_all_charges
                else {
                    "molecule-ID": 1,
                    "type": atom_type_ids[atom_type_by_atom_id[atom.atom_id]],
                    "x": wrapped_cartesian_positions[atom.atom_id][0],
                    "y": wrapped_cartesian_positions[atom.atom_id][1],
                    "z": wrapped_cartesian_positions[atom.atom_id][2],
                }
            )
            for atom in parsed.atoms
        ],
        index=range(1, len(parsed.atoms) + 1),
    )
    force_field: dict[str, pd.DataFrame] = {
        "PairIJ Coeffs": pd.DataFrame(
            pairij_coeff_rows,
            columns=["id1", "id2", "coeff1", "coeff2"],
        ),
        "Bond Coeffs": _coeff_rows_to_dataframe(bond_coeff_rows),
    }
    if angle_coeff_rows:
        force_field["Angle Coeffs"] = _coeff_rows_to_dataframe(angle_coeff_rows)
    if dihedral_coeff_rows:
        force_field["Dihedral Coeffs"] = _coeff_rows_to_dataframe(dihedral_coeff_rows)
    if improper_coeff_rows:
        force_field["Improper Coeffs"] = _coeff_rows_to_dataframe(improper_coeff_rows)
    topology: dict[str, pd.DataFrame] = {
        "Bonds": _topology_rows_to_dataframe(bond_topology_rows, ("type", "atom1", "atom2")),
    }
    if angle_topology_rows:
        topology["Angles"] = _topology_rows_to_dataframe(angle_topology_rows, ("type", "atom1", "atom2", "atom3"))
    if dihedral_topology_rows:
        topology["Dihedrals"] = _topology_rows_to_dataframe(
            dihedral_topology_rows,
            ("type", "atom1", "atom2", "atom3", "atom4"),
        )
    if improper_topology_rows:
        topology["Impropers"] = _topology_rows_to_dataframe(
            improper_topology_rows,
            ("type", "atom1", "atom2", "atom3", "atom4"),
        )
    lammps_data = LammpsData(
        box=box,
        masses=masses,
        atoms=atoms,
        force_field=force_field,
        topology=topology,
        atom_style="full" if has_all_charges else "molecular",
    )
    data_text = lammps_data.get_str(distance=10, charge=4, hybrid=True)
    data_text = data_text.replace(
        "Generated by pymatgen.io.lammps.data.LammpsData",
        "LAMMPS data file written by cofkit via pymatgen/Open Babel/UFF",
        1,
    )
    data_text = data_text.replace(
        "\nAtoms\n\n",
        "\nAtoms # full\n\n" if has_all_charges else "\nAtoms # molecular\n\n",
        1,
    )
    has_image_flags = image_conflicts == 0
    if has_image_flags:
        data_text = _inject_atom_image_flags(data_text, atom_images)

    warnings: list[str] = []
    if image_conflicts > 0:
        warnings.append(
            "Periodic spanning-tree unwrapping for UFF atom typing encountered one or more cell-cycle conflicts. "
            "cofkit kept the wrapped CIF coordinates for the LAMMPS data file and used the spanning tree only for "
            "local UFF atom typing."
        )
    parameter_sources = {
        "forcefield": "UFF",
        "atom_typing": "Open Babel UFF atom types from the explicit CIF bond graph",
        "bonded_parameters": "cofkit UFF formulas adapted from reference_repositories/lammps_interface",
        "nonbond_parameters": "UFF.prm plus Lorentz-Berthelot mixing",
    }
    uff_parameter_file = _bundled_uff_parameter_file()
    parameter_sources["reference_parameter_file"] = str(uff_parameter_file)
    parameter_sources["reference_parameter_sha256"] = _bundled_uff_parameter_sha256()
    parameter_sources["reference_logic"] = "reference_repositories/lammps_interface"
    if has_all_charges:
        parameter_sources["electrostatics"] = "Explicit atom charges carried into the LAMMPS data file"
        parameter_sources["charge_assignment"] = "Atom charges read from the CIF atom-site charge loop"

    return _PreparedLammpsSystem(
        data_text=data_text,
        has_image_flags=has_image_flags,
        has_charges=has_all_charges,
        atom_type_symbols={index: label for index, label in enumerate(ordered_type_labels, start=1)},
        n_bond_types=len(bond_coeff_rows),
        n_angle_types=len(angle_coeff_rows),
        n_dihedral_types=len(dihedral_coeff_rows),
        n_improper_types=len(improper_coeff_rows),
        n_dihedrals=len(dihedral_topology_rows),
        n_impropers=len(improper_topology_rows),
        angle_styles=_ordered_coeff_styles(angle_coeff_rows),
        dihedral_styles=("harmonic",) if dihedral_coeff_rows else (),
        improper_styles=("fourier",) if improper_coeff_rows else (),
        forcefield_backend="uff_openbabel_explicit_graph_pymatgen",
        parameter_sources=parameter_sources,
        n_charged_atoms=len(parsed.atoms) if has_all_charges else 0,
        net_charge=(sum(float(atom.charge) for atom in parsed.atoms) if has_all_charges else None),
        warnings=tuple(warnings),
    )


def _inject_atom_image_flags(
    data_text: str,
    atom_images: Mapping[int, tuple[int, int, int]],
) -> str:
    lines = data_text.splitlines()
    atoms_header_index = None
    for candidate in ("Atoms # molecular", "Atoms # full"):
        try:
            atoms_header_index = lines.index(candidate)
            break
        except ValueError:
            continue
    if atoms_header_index is None:
        return data_text
    atom_line_index = atoms_header_index + 2
    while atom_line_index < len(lines):
        line = lines[atom_line_index]
        if not line.strip():
            break
        fields = line.split()
        atom_id = int(fields[0])
        image = atom_images.get(atom_id, (0, 0, 0))
        lines[atom_line_index] = f"{line} {image[0]} {image[1]} {image[2]}"
        atom_line_index += 1
    return "\n".join(lines) + ("\n" if data_text.endswith("\n") else "")


def _bundled_uff_parameter_file() -> Path:
    path = Path(__file__).resolve().parent / "data" / "forcefields" / "openbabel" / "3.1.0" / "UFF.prm"
    if not path.is_file():
        raise LammpsConfigurationError(f"Bundled UFF.prm is missing from the package data: {path}")
    return path


@lru_cache(maxsize=1)
def _bundled_uff_parameter_sha256() -> str:
    digest = hashlib.sha256(_bundled_uff_parameter_file().read_bytes()).hexdigest()
    if digest != _PINNED_UFF_PARAMETER_SHA256:
        raise LammpsConfigurationError(
            "Bundled UFF.prm does not match the pinned Open Babel reference hash. "
            f"Expected {_PINNED_UFF_PARAMETER_SHA256}, observed {digest}."
        )
    return digest


def _require_uff_support() -> None:
    missing: list[str] = []
    if gemmi is None:
        missing.append("gemmi")
    if ob is None:
        missing.append("openbabel")
    if any(item is None for item in (LammpsBox, LammpsData)):
        missing.append("pymatgen")
    if pd is None:
        missing.append("pandas")
    try:
        _bundled_uff_parameter_sha256()
    except LammpsConfigurationError:
        missing.append("bundled UFF.prm")
    if missing:
        raise LammpsConfigurationError(
            "UFF-backed LAMMPS optimization requires these Python packages in the current environment: "
            + ", ".join(missing)
        )


def _ordered_unique_uff_atom_types(
    atoms: tuple[_CifAtomRecord, ...],
    atom_type_by_atom_id: dict[int, str],
) -> tuple[list[str], dict[str, int], dict[str, str]]:
    ordered_types: list[str] = []
    atom_type_ids: dict[str, int] = {}
    representative_symbols: dict[str, str] = {}
    for atom in atoms:
        atom_type = atom_type_by_atom_id[atom.atom_id]
        if atom_type in atom_type_ids:
            continue
        ordered_types.append(atom_type)
        atom_type_ids[atom_type] = len(ordered_types)
        representative_symbols[atom_type] = atom.symbol
    return ordered_types, atom_type_ids, representative_symbols


@lru_cache(maxsize=1)
def _load_uff_parameters() -> dict[str, _UffAtomParameters]:
    parameter_file = _bundled_uff_parameter_file()
    _bundled_uff_parameter_sha256()
    parameters: dict[str, _UffAtomParameters] = {}
    for raw_line in parameter_file.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw_line.strip()
        if not line.startswith("param "):
            continue
        parts = line.split()
        if len(parts) < 13:
            continue
        parameters[parts[1]] = _UffAtomParameters(*(float(value) for value in parts[2:13]))
    if not parameters:
        raise LammpsConfigurationError(f"No UFF parameter rows were parsed from {parameter_file}.")
    return parameters


def _build_uff_pairij_rows(
    *,
    ordered_type_labels: list[str],
    atom_type_ids: dict[str, int],
    parameters: dict[str, _UffAtomParameters],
) -> list[list[float]]:
    rows: list[list[float]] = []
    for left_index, left_label in enumerate(ordered_type_labels):
        left_params = _uff_parameters_for_type(left_label, parameters)
        for right_label in ordered_type_labels[left_index:]:
            right_params = _uff_parameters_for_type(right_label, parameters)
            epsilon = math.sqrt(left_params.d1 * right_params.d1)
            sigma = ((left_params.x1 / _UFF_RMIN_TO_SIGMA_FACTOR) + (right_params.x1 / _UFF_RMIN_TO_SIGMA_FACTOR)) / 2.0
            rows.append([atom_type_ids[left_label], atom_type_ids[right_label], float(epsilon), float(sigma)])
    return rows


def _build_uff_bond_terms(
    *,
    parsed: _ParsedExplicitBondCif,
    atom_type_by_atom_id: dict[int, str],
    parameters: dict[str, _UffAtomParameters],
) -> tuple[list[list[int]], list[list[float]], dict[tuple[int, int], tuple[float, float]]]:
    type_ids: dict[tuple[str, str, float], int] = {}
    coeff_rows: list[list[float]] = []
    topology_rows: list[list[int]] = []
    bond_reference: dict[tuple[int, int], tuple[float, float]] = {}
    for bond in parsed.bonds:
        assert bond.bond_order is not None
        left_type = atom_type_by_atom_id[bond.atom_id_1]
        right_type = atom_type_by_atom_id[bond.atom_id_2]
        bond_order = normalize_bond_order(bond.bond_order)
        coefficient_key = _canonical_bond_term_key(left_type, right_type, bond_order)
        type_id = type_ids.get(coefficient_key)
        if type_id is None:
            force_constant, equilibrium_distance = _compute_uff_bond_coefficients(
                left_type=left_type,
                right_type=right_type,
                bond_order=bond_order,
                parameters=parameters,
            )
            coeff_rows.append([force_constant, equilibrium_distance])
            type_id = len(coeff_rows)
            type_ids[coefficient_key] = type_id
        else:
            equilibrium_distance = coeff_rows[type_id - 1][1]
        topology_rows.append([type_id, bond.atom_id_1, bond.atom_id_2])
        bond_reference[_canonical_atom_pair(bond.atom_id_1, bond.atom_id_2)] = (bond_order, equilibrium_distance)
    return topology_rows, coeff_rows, bond_reference


def _build_uff_angle_terms(
    *,
    parsed: _ParsedExplicitBondCif,
    atom_type_by_atom_id: dict[int, str],
    parameters: dict[str, _UffAtomParameters],
    bond_reference: dict[tuple[int, int], tuple[float, float]],
) -> tuple[list[list[int]], list[list[object]]]:
    type_ids: dict[tuple[object, ...], int] = {}
    coeff_rows: list[list[object]] = []
    topology_rows: list[list[int]] = []
    for angle in parsed.angles:
        coeffs = _compute_uff_angle_coefficients(
            angle=angle,
            atom_type_by_atom_id=atom_type_by_atom_id,
            parameters=parameters,
            bond_reference=bond_reference,
        )
        coefficient_key = _freeze_coeff_values(coeffs)
        type_id = type_ids.get(coefficient_key)
        if type_id is None:
            coeff_rows.append(list(coeffs))
            type_id = len(coeff_rows)
            type_ids[coefficient_key] = type_id
        topology_rows.append([type_id, angle.atom_id_1, angle.atom_id_2, angle.atom_id_3])
    return topology_rows, coeff_rows


def _build_uff_dihedral_terms(
    *,
    parsed: _ParsedExplicitBondCif,
    atom_type_by_atom_id: dict[int, str],
    parameters: dict[str, _UffAtomParameters],
) -> tuple[list[list[int]], list[list[float]]]:
    type_ids: dict[tuple[object, ...], int] = {}
    coeff_rows: list[list[float]] = []
    topology_rows: list[list[int]] = []
    for dihedral in parsed.dihedrals:
        coeffs = _compute_uff_dihedral_coefficients(
            dihedral=dihedral,
            parsed=parsed,
            atom_type_by_atom_id=atom_type_by_atom_id,
            parameters=parameters,
        )
        if coeffs is None:
            continue
        coefficient_key = _freeze_coeff_values(coeffs)
        type_id = type_ids.get(coefficient_key)
        if type_id is None:
            coeff_rows.append(list(coeffs))
            type_id = len(coeff_rows)
            type_ids[coefficient_key] = type_id
        topology_rows.append([type_id, dihedral.atom_id_1, dihedral.atom_id_2, dihedral.atom_id_3, dihedral.atom_id_4])
    return topology_rows, coeff_rows


def _build_uff_improper_terms(
    *,
    parsed: _ParsedExplicitBondCif,
    atom_type_by_atom_id: dict[int, str],
    parameters: dict[str, _UffAtomParameters],
) -> tuple[list[list[int]], list[list[float]]]:
    type_ids: dict[tuple[object, ...], int] = {}
    coeff_rows: list[list[float]] = []
    topology_rows: list[list[int]] = []
    degrees = _degree_by_atom(parsed)
    for improper in parsed.impropers:
        coeffs = _compute_uff_improper_coefficients(
            improper=improper,
            parsed=parsed,
            atom_type_by_atom_id=atom_type_by_atom_id,
            degrees=degrees,
        )
        if coeffs is None:
            continue
        coefficient_key = _freeze_coeff_values(coeffs)
        type_id = type_ids.get(coefficient_key)
        if type_id is None:
            coeff_rows.append(list(coeffs))
            type_id = len(coeff_rows)
            type_ids[coefficient_key] = type_id
        topology_rows.append([type_id, improper.atom_id_1, improper.atom_id_2, improper.atom_id_3, improper.atom_id_4])
    return topology_rows, coeff_rows


def _coeff_rows_to_dataframe(rows: Sequence[Sequence[object]]) -> pd.DataFrame:
    width = max(len(row) for row in rows)
    dataframe = pd.DataFrame(list(rows), columns=[f"coeff{i}" for i in range(1, width + 1)])
    dataframe.index = range(1, len(dataframe) + 1)
    return dataframe


def _topology_rows_to_dataframe(rows: Sequence[Sequence[int]], columns: Sequence[str]) -> pd.DataFrame:
    dataframe = pd.DataFrame(list(rows), columns=list(columns))
    dataframe.index = range(1, len(dataframe) + 1)
    return dataframe


def _compute_uff_bond_coefficients(
    *,
    left_type: str,
    right_type: str,
    bond_order: float,
    parameters: dict[str, _UffAtomParameters],
) -> tuple[float, float]:
    left_parameters = _uff_parameters_for_type(left_type, parameters)
    right_parameters = _uff_parameters_for_type(right_type, parameters)
    rbo = -0.1332 * (left_parameters.r1 + right_parameters.r1) * math.log(bond_order)
    ren = (
        left_parameters.r1
        * right_parameters.r1
        * ((math.sqrt(left_parameters.xi) - math.sqrt(right_parameters.xi)) ** 2)
        / (left_parameters.xi * left_parameters.r1 + right_parameters.xi * right_parameters.r1)
    )
    equilibrium_distance = left_parameters.r1 + right_parameters.r1 + rbo - ren
    force_constant = 664.12 * (left_parameters.z1 * right_parameters.z1) / (equilibrium_distance**3) / 2.0
    return float(force_constant), float(equilibrium_distance)


def _compute_uff_angle_coefficients(
    *,
    angle: _CifAngleRecord,
    atom_type_by_atom_id: dict[int, str],
    parameters: dict[str, _UffAtomParameters],
    bond_reference: dict[tuple[int, int], tuple[float, float]],
) -> tuple[object, ...]:
    left_type = atom_type_by_atom_id[angle.atom_id_1]
    center_type = atom_type_by_atom_id[angle.atom_id_2]
    right_type = atom_type_by_atom_id[angle.atom_id_3]
    center_parameters = _uff_parameters_for_type(center_type, parameters)
    theta0 = center_parameters.theta0
    angle_family = _uff_angle_family(center_type)
    left_bond = bond_reference[_canonical_atom_pair(angle.atom_id_1, angle.atom_id_2)]
    right_bond = bond_reference[_canonical_atom_pair(angle.atom_id_2, angle.atom_id_3)]
    r_ab = left_bond[1]
    r_bc = right_bond[1]
    cos_theta0 = math.cos(math.radians(theta0))
    za = _uff_parameters_for_type(left_type, parameters).z1
    zc = _uff_parameters_for_type(right_type, parameters).z1
    r_ac = math.sqrt(r_ab * r_ab + r_bc * r_bc - 2.0 * r_ab * r_bc * cos_theta0)
    beta = 664.12 / r_ab / r_bc
    ka = beta * (za * zc / (r_ac**5.0)) * r_ab * r_bc
    ka *= 3.0 * r_ab * r_bc * (1.0 - cos_theta0 * cos_theta0) - r_ac * r_ac * cos_theta0

    if angle_family in {"linear", "trigonal-planar", "square-planar", "octahedral"} or (
        angle_family == "tetrahedral" and int(theta0) == 90
    ):
        if angle_family == "linear":
            return ("cosine/periodic", float(ka), 1, 1)
        if angle_family == "tetrahedral":
            return ("cosine/periodic", float(ka), -1, 2)
        if angle_family == "trigonal-planar":
            return ("cosine/periodic", float(ka), -1, 3)
        return ("cosine/periodic", float(ka), 1, 4)

    sin_theta0 = math.sin(math.radians(theta0))
    if abs(sin_theta0) <= 1.0e-12:
        raise LammpsInputError(
            f"UFF angle coefficient generation became singular for angle {angle.atom_id_1}-{angle.atom_id_2}-{angle.atom_id_3}."
        )
    c2 = 1.0 / (4.0 * sin_theta0 * sin_theta0)
    c1 = -4.0 * c2 * cos_theta0
    c0 = c2 * (2.0 * cos_theta0 * cos_theta0 + 1.0)
    return ("fourier", float(ka), float(c0), float(c1), float(c2))


def _compute_uff_dihedral_coefficients(
    *,
    dihedral: _CifDihedralRecord,
    parsed: _ParsedExplicitBondCif,
    atom_type_by_atom_id: dict[int, str],
    parameters: dict[str, _UffAtomParameters],
) -> tuple[float, int, int] | None:
    left_center_type = atom_type_by_atom_id[dihedral.atom_id_2]
    right_center_type = atom_type_by_atom_id[dihedral.atom_id_3]
    left_center_symbol = parsed.atoms[dihedral.atom_id_2 - 1].symbol
    right_center_symbol = parsed.atoms[dihedral.atom_id_3 - 1].symbol
    if gemmi.Element(left_center_symbol).is_metal or gemmi.Element(right_center_symbol).is_metal:
        return None

    left_family = _uff_hybridization_family(left_center_type)
    right_family = _uff_hybridization_family(right_center_type)
    torsion_order = _bond_order_for_pair(parsed, dihedral.atom_id_2, dihedral.atom_id_3)
    degree_by_atom = _degree_by_atom(parsed)
    multiplicity = float((degree_by_atom[dihedral.atom_id_2] - 1) * (degree_by_atom[dihedral.atom_id_3] - 1))
    if multiplicity <= 0.0:
        return None

    left_atomic_number = ob.GetAtomicNum(left_center_symbol)
    right_atomic_number = ob.GetAtomicNum(right_center_symbol)
    group6_atomic_numbers = {8, 16, 34, 52, 84}
    left_sp2_like = left_family in {"sp2", "aromatic"}
    right_sp2_like = right_family in {"sp2", "aromatic"}
    all_sp3 = left_family == "sp3" and right_family == "sp3"
    all_sp2 = left_sp2_like and right_sp2_like
    mixed_case = (left_sp2_like and right_family == "sp3") or (left_family == "sp3" and right_sp2_like)

    phase = 0.0
    periodicity = 0
    amplitude = 0.0
    if all_sp3:
        phase = 60.0
        periodicity = 3
        vi = _uff_parameters_for_type(left_center_type, parameters).vi
        vj = _uff_parameters_for_type(right_center_type, parameters).vi
        if left_atomic_number == 8:
            vi = 2.0
            periodicity = 2
            phase = 90.0
        elif left_atomic_number in group6_atomic_numbers - {8}:
            vi = 6.8
            periodicity = 2
            phase = 90.0
        if right_atomic_number == 8:
            vj = 2.0
            periodicity = 2
            phase = 90.0
        elif right_atomic_number in group6_atomic_numbers - {8}:
            vj = 6.8
            periodicity = 2
            phase = 90.0
        amplitude = math.sqrt(vi * vj)
    elif all_sp2:
        ui = _uff_parameters_for_type(left_center_type, parameters).uj
        uj = _uff_parameters_for_type(right_center_type, parameters).uj
        phase = 180.0
        periodicity = 2
        amplitude = 5.0 * math.sqrt(ui * uj) * (1.0 + 4.18 * math.log(torsion_order))
    elif mixed_case:
        phase = 180.0
        periodicity = 3
        amplitude = 2.0
        if right_family == "sp3" and right_atomic_number in group6_atomic_numbers:
            periodicity = 2
            phase = 90.0
        elif left_family == "sp3" and left_atomic_number in group6_atomic_numbers:
            periodicity = 2
            phase = 90.0
        if periodicity == 2:
            ui = _uff_parameters_for_type(left_center_type, parameters).uj
            uj = _uff_parameters_for_type(right_center_type, parameters).uj
            amplitude = 5.0 * math.sqrt(ui * uj) * (1.0 + 4.18 * math.log(torsion_order))
    else:
        return None

    amplitude /= multiplicity
    if abs(amplitude) <= 1.0e-12:
        return None
    d_value = int(round(-math.cos(math.radians(periodicity * phase))))
    return float(0.5 * amplitude), d_value, int(periodicity)


def _compute_uff_improper_coefficients(
    *,
    improper: _CifImproperRecord,
    parsed: _ParsedExplicitBondCif,
    atom_type_by_atom_id: dict[int, str],
    degrees: dict[int, int],
) -> tuple[float, float, float, float, int] | None:
    center_type = atom_type_by_atom_id[improper.atom_id_2]
    center_atomic_number = ob.GetAtomicNum(parsed.atoms[improper.atom_id_2 - 1].symbol)
    allowed_atomic_numbers = {6, 7, 8, 15, 33, 51, 83}
    if center_atomic_number not in allowed_atomic_numbers:
        return None

    if center_type in {"N_3", "N_2", "N_R", "O_2", "O_R"}:
        return (2.0, 1.0, -1.0, 0.0, 0)
    if center_type in {"P_3+3", "As3+3", "Sb3+3", "Bi3+3"}:
        if center_type == "P_3+3":
            phi = math.radians(84.4339)
        elif center_type == "As3+3":
            phi = math.radians(86.9735)
        elif center_type == "Sb3+3":
            phi = math.radians(87.7047)
        else:
            phi = math.radians(90.0)
        c1 = -4.0 * math.cos(phi)
        c2 = 1.0
        c0 = -c1 * math.cos(phi) + c2 * math.cos(2.0 * phi)
        return (22.0 / 3.0, float(c0), float(c1), float(c2), 0)
    if center_type in {"C_2", "C_R"}:
        force_constant = 2.0
        for atom_id in (improper.atom_id_1, improper.atom_id_3, improper.atom_id_4):
            if atom_type_by_atom_id[atom_id] == "O_2" and degrees[atom_id] == 1:
                force_constant = 50.0 / 3.0
                break
        return (force_constant, 1.0, -1.0, 0.0, 0)
    return None


def _uff_parameters_for_type(
    atom_type: str,
    parameters: dict[str, _UffAtomParameters],
) -> _UffAtomParameters:
    try:
        return parameters[atom_type]
    except KeyError as exc:
        raise LammpsInputError(f"UFF parameters are unavailable for atom type {atom_type!r}.") from exc


def _canonical_atom_pair(left_atom_id: int, right_atom_id: int) -> tuple[int, int]:
    return (left_atom_id, right_atom_id) if left_atom_id <= right_atom_id else (right_atom_id, left_atom_id)


def _canonical_bond_term_key(left_type: str, right_type: str, bond_order: float) -> tuple[str, str, float]:
    forward = (left_type, right_type, bond_order)
    reverse = (right_type, left_type, bond_order)
    return min(forward, reverse)


def _freeze_coeff_values(values: Sequence[object]) -> tuple[object, ...]:
    frozen: list[object] = []
    for value in values:
        if isinstance(value, float):
            frozen.append(round(value, 8))
        else:
            frozen.append(value)
    return tuple(frozen)


def _ordered_coeff_styles(rows: Sequence[Sequence[object]]) -> tuple[str, ...]:
    styles: list[str] = []
    for row in rows:
        if not row or not isinstance(row[0], str):
            continue
        style = row[0]
        if style not in styles:
            styles.append(style)
    return tuple(styles)


def _uff_angle_family(atom_type: str) -> str:
    if atom_type.endswith("_R") or atom_type.endswith("R"):
        return "trigonal-planar"
    if len(atom_type) < 3:
        return "linear"
    coordination = atom_type[2]
    if coordination == "1":
        return "linear"
    if coordination in {"2", "R"}:
        return "trigonal-planar"
    if coordination == "3":
        return "tetrahedral"
    if coordination == "4":
        return "square-planar"
    if coordination == "5":
        return "trigonal-bipyramidal"
    if coordination == "6":
        return "octahedral"
    if coordination == "8":
        return "cubic-antiprism"
    return "linear"


def _uff_hybridization_family(atom_type: str) -> str:
    if atom_type.endswith("_R") or atom_type.endswith("R"):
        return "aromatic"
    if len(atom_type) < 3:
        return "sp"
    coordination = atom_type[2]
    if coordination == "1":
        return "sp"
    if coordination == "2":
        return "sp2"
    if coordination in {"3", "4", "5", "6", "8"}:
        return "sp3"
    return "sp3"


def _degree_by_atom(parsed: _ParsedExplicitBondCif) -> dict[int, int]:
    degrees = {atom.atom_id: 0 for atom in parsed.atoms}
    for bond in parsed.bonds:
        degrees[bond.atom_id_1] += 1
        degrees[bond.atom_id_2] += 1
    return degrees


def _bond_order_for_pair(parsed: _ParsedExplicitBondCif, left_atom_id: int, right_atom_id: int) -> float:
    pair = _canonical_atom_pair(left_atom_id, right_atom_id)
    for bond in parsed.bonds:
        if _canonical_atom_pair(bond.atom_id_1, bond.atom_id_2) != pair:
            continue
        if bond.bond_order is None:
            raise LammpsInputError(
                f"Missing explicit bond order for bond {bond.label_1}-{bond.label_2} while assigning UFF torsions."
            )
        return normalize_bond_order(bond.bond_order)
    raise LammpsInputError(f"Internal error: could not resolve bond order for atom pair {pair!r}.")


def _compute_unwrapped_atom_images(parsed: _ParsedExplicitBondCif) -> tuple[dict[int, tuple[int, int, int]], int]:
    neighbors: dict[int, list[tuple[int, tuple[int, int, int]]]] = {atom.atom_id: [] for atom in parsed.atoms}
    for bond in parsed.bonds:
        delta = _sub_shift(bond.shift_2, bond.shift_1)
        neighbors[bond.atom_id_1].append((bond.atom_id_2, delta))
        neighbors[bond.atom_id_2].append((bond.atom_id_1, _neg_shift(delta)))

    images: dict[int, tuple[int, int, int]] = {}
    conflict_count = 0
    for atom in parsed.atoms:
        if atom.atom_id in images:
            continue
        images[atom.atom_id] = (0, 0, 0)
        stack = [atom.atom_id]
        while stack:
            current = stack.pop()
            current_image = images[current]
            for neighbor_id, delta in neighbors[current]:
                proposed = _add_shift(current_image, delta)
                if neighbor_id not in images:
                    images[neighbor_id] = proposed
                    stack.append(neighbor_id)
                    continue
                if images[neighbor_id] != proposed:
                    conflict_count += 1
    return images, conflict_count


def _unwrapped_cartesian_positions(
    parsed: _ParsedExplicitBondCif,
    atom_images: dict[int, tuple[int, int, int]],
) -> dict[int, tuple[float, float, float]]:
    positions: dict[int, tuple[float, float, float]] = {}
    for atom in parsed.atoms:
        image = atom_images[atom.atom_id]
        positions[atom.atom_id] = _fractional_to_lammps(
            (
                atom.fractional_position[0] + image[0],
                atom.fractional_position[1] + image[1],
                atom.fractional_position[2] + image[2],
            ),
            parsed.lammps_basis,
        )
    return positions


def _build_openbabel_molecule(
    parsed: _ParsedExplicitBondCif,
    cartesian_positions: dict[int, tuple[float, float, float]],
):
    molecule = ob.OBMol()
    for atom in parsed.atoms:
        atomic_number = ob.GetAtomicNum(atom.symbol)
        if atomic_number <= 0:
            raise LammpsInputError(f"Open Babel could not resolve an atomic number for element symbol {atom.symbol!r}.")
        ob_atom = molecule.NewAtom()
        ob_atom.SetAtomicNum(atomic_number)
        x, y, z = cartesian_positions[atom.atom_id]
        ob_atom.SetVector(x, y, z)
    aromatic_atom_ids: set[int] = set()
    for bond in parsed.bonds:
        assert bond.bond_order is not None
        bond_order = 1
        if not is_aromatic_bond_order(bond.bond_order):
            bond_order = max(1, int(round(normalize_bond_order(bond.bond_order))))
        if not molecule.AddBond(bond.atom_id_1, bond.atom_id_2, bond_order):
            raise LammpsInputError(
                f"Open Babel rejected explicit bond {bond.label_1}-{bond.label_2} while preparing the UFF model."
            )
        if is_aromatic_bond_order(bond.bond_order):
            ob_bond = molecule.GetBond(bond.atom_id_1, bond.atom_id_2)
            if ob_bond is None:
                raise LammpsInputError(
                    f"Open Babel could not resolve explicit aromatic bond {bond.label_1}-{bond.label_2}."
                )
            ob_bond.SetAromatic(True)
            aromatic_atom_ids.add(bond.atom_id_1)
            aromatic_atom_ids.add(bond.atom_id_2)
    for atom_id in aromatic_atom_ids:
        ob_atom = molecule.GetAtom(atom_id)
        if ob_atom is not None:
            ob_atom.SetAromatic(True)
    if aromatic_atom_ids:
        molecule.SetAromaticPerceived(True)
        molecule.FindRingAtomsAndBonds()
        molecule.SetRingAtomsAndBondsPerceived(True)
    return molecule


def _assign_openbabel_uff_atom_types(ob_molecule) -> dict[int, str]:
    forcefield = ob.OBForceField.FindForceField("UFF")
    if not forcefield:
        raise LammpsConfigurationError(
            "Open Babel is installed, but the UFF forcefield backend is unavailable in this build."
        )
    if not forcefield.Setup(ob_molecule):
        raise LammpsInputError(
            "Open Babel could not set up a UFF-typed model from this CIF-derived explicit-bond graph."
        )
    forcefield.GetAtomTypes(ob_molecule)

    atom_types: dict[int, str] = {}
    for atom in ob.OBMolAtomIter(ob_molecule):
        data = atom.GetData("FFAtomType")
        atom_type = data.GetValue() if data is not None else ""
        if not atom_type:
            raise LammpsInputError(
                f"Open Babel did not assign a UFF atom type to atom index {atom.GetIdx()}."
            )
        atom_types[atom.GetIdx()] = atom_type
    return atom_types


def _parse_explicit_bond_cif(cif_path: Path) -> _ParsedExplicitBondCif:
    if gemmi is None:
        raise LammpsConfigurationError(
            "gemmi is required for CIF-backed LAMMPS optimization. Install gemmi into this environment."
        )

    block = gemmi.cif.read_file(str(cif_path)).sole_block()
    _require_p1_block(block, cif_path)
    small = gemmi.make_small_structure_from_block(block)
    if len(small.sites) == 0:
        raise LammpsInputError(f"CIF file does not contain any atom sites: {cif_path}")

    a, b, c = float(small.cell.a), float(small.cell.b), float(small.cell.c)
    alpha, beta, gamma = float(small.cell.alpha), float(small.cell.beta), float(small.cell.gamma)
    basis = _lammps_basis_from_cell(a, b, c, alpha, beta, gamma)

    atoms: list[_CifAtomRecord] = []
    label_to_atom_id: dict[str, int] = {}
    for atom_id, site in enumerate(small.sites, start=1):
        label = str(site.label or f"a{atom_id}")
        if label in label_to_atom_id:
            raise LammpsInputError(f"CIF atom labels must be unique for LAMMPS optimization: {label}")
        symbol = str(site.element.name)
        if not symbol:
            raise LammpsInputError(f"CIF atom {label} is missing a valid element symbol.")
        label_to_atom_id[label] = atom_id
        atoms.append(
            _CifAtomRecord(
                atom_id=atom_id,
                label=label,
                symbol=symbol,
                occupancy=float(site.occ),
                fractional_position=(float(site.fract.x), float(site.fract.y), float(site.fract.z)),
            )
        )
    atom_charge_result = _extract_atom_site_charges_from_block(block, label_to_atom_id.keys())
    if atom_charge_result.charges_by_label:
        atoms = [
            replace(atom, charge=atom_charge_result.charges_by_label[atom.label])
            for atom in atoms
        ]
    positions = {atom.atom_id: atom.fractional_position for atom in atoms}

    bond_table = block.find(
        [
            "_geom_bond_atom_site_label_1",
            "_geom_bond_atom_site_label_2",
            "_geom_bond_site_symmetry_1",
            "_geom_bond_site_symmetry_2",
            "_geom_bond_distance",
        ]
    )
    bond_type_loop = block.find_loop("_ccdc_geom_bond_type")
    if len(bond_type_loop) == 0:
        bond_type_loop = block.find_loop("_geom_bond_type")
    if len(bond_table) == 0:
        raise LammpsInputError(
            f"CIF file does not contain an explicit _geom_bond_* loop and cannot drive the LAMMPS bonded cleanup: {cif_path}"
        )

    bonds: list[_CifBondRecord] = []
    for bond_id in range(len(bond_table)):
        row = bond_table[bond_id]
        label_1 = str(row[0]).strip()
        label_2 = str(row[1]).strip()
        symmetry_1 = str(row[2]).strip() or "."
        symmetry_2 = str(row[3]).strip() or "."
        try:
            equilibrium_distance = float(row[4])
        except ValueError as exc:
            raise LammpsInputError(
                f"Bond distance for {label_1}-{label_2} in {cif_path} is not numeric: {row[4]!r}"
            ) from exc

        if label_1 not in label_to_atom_id or label_2 not in label_to_atom_id:
            raise LammpsInputError(
                f"Bond loop references unknown atom labels {label_1!r} and/or {label_2!r} in {cif_path}"
            )

        effective_shift_1, effective_shift_2, equilibrium_distance = _resolve_effective_bond_geometry(
            left_fractional=positions[label_to_atom_id[label_1]],
            right_fractional=positions[label_to_atom_id[label_2]],
            raw_shift_1=_parse_p1_symmetry_shift(symmetry_1),
            raw_shift_2=_parse_p1_symmetry_shift(symmetry_2),
            raw_symmetry_1=symmetry_1,
            raw_symmetry_2=symmetry_2,
            basis=basis,
        )
        bond_order = cif_type_to_bond_order(str(bond_type_loop[bond_id]).strip()) if len(bond_type_loop) > bond_id else None
        bonds.append(
            _CifBondRecord(
                bond_id=bond_id + 1,
                atom_id_1=label_to_atom_id[label_1],
                atom_id_2=label_to_atom_id[label_2],
                label_1=label_1,
                label_2=label_2,
                symmetry_1=symmetry_1,
                symmetry_2=symmetry_2,
                shift_1=effective_shift_1,
                shift_2=effective_shift_2,
                equilibrium_distance=equilibrium_distance,
                bond_order=bond_order,
            )
        )

    angles = _derive_angles(tuple(atoms), tuple(bonds), basis)
    dihedrals = _derive_dihedrals(tuple(atoms), tuple(bonds))
    impropers = _derive_impropers(tuple(atoms), tuple(bonds))
    return _ParsedExplicitBondCif(
        source_path=cif_path,
        data_name=str(block.name or cif_path.stem),
        cell_parameters=(a, b, c, alpha, beta, gamma),
        lammps_basis=basis,
        atoms=tuple(atoms),
        bonds=tuple(bonds),
        angles=angles,
        dihedrals=dihedrals,
        impropers=impropers,
    )


def _require_p1_block(block, cif_path: Path) -> None:
    if gemmi is None:  # pragma: no cover - guarded earlier
        return

    number_value = (block.find_value("_space_group_IT_number") or "").strip()
    if number_value and number_value != "1":
        raise LammpsInputError(
            f"Initial LAMMPS optimization support is limited to P1 CIFs with explicit bond loops. "
            f"{cif_path} declares _space_group_IT_number {number_value!r}."
        )

    name_value = (
        block.find_value("_space_group_name_H-M_alt")
        or block.find_value("_symmetry_space_group_name_H-M")
        or ""
    ).strip().replace(" ", "").upper()
    if name_value and name_value not in {"P1", "P1'", "'P1'", '"P1"'}:
        raise LammpsInputError(
            f"Initial LAMMPS optimization support is limited to P1 CIFs with explicit bond loops. "
            f"{cif_path} declares space group {name_value!r}."
        )


def _parse_p1_symmetry_shift(value: str) -> tuple[int, int, int]:
    raw = value.strip()
    if not raw or raw == ".":
        return (0, 0, 0)
    if "_" not in raw:
        if raw == "1":
            return (0, 0, 0)
        raise LammpsInputError(
            f"Unsupported CIF bond symmetry token {value!r}. Initial LAMMPS support expects P1-style identity operators."
        )
    operation, translation = raw.split("_", 1)
    if operation != "1":
        raise LammpsInputError(
            f"Unsupported CIF bond symmetry token {value!r}. Initial LAMMPS support expects only identity operator 1."
        )
    if len(translation) != 3 or not translation.isdigit():
        raise LammpsInputError(f"Unsupported CIF bond symmetry translation {value!r}.")
    return tuple(int(ch) - 5 for ch in translation)  # type: ignore[return-value]


def _extract_atom_site_charges_from_cif(
    cif_path: Path,
    expected_labels: Sequence[str] | None = None,
) -> _AtomChargeExtractionResult:
    if gemmi is None:
        raise LammpsConfigurationError(
            "gemmi is required for CIF-backed charge extraction in the LAMMPS workflow."
        )
    block = gemmi.cif.read_file(str(cif_path)).sole_block()
    return _extract_atom_site_charges_from_block(block, expected_labels)


def _extract_atom_site_charges_from_block(
    block,
    expected_labels: Sequence[str] | None = None,
) -> _AtomChargeExtractionResult:
    charge_column_names = ("_atom_site_charge", "_atom_site_partial_charge")
    expected_label_sequence = tuple(str(label) for label in (expected_labels or ()))
    expected_label_set = set(expected_label_sequence)
    for charge_column_name in charge_column_names:
        table = block.find(["_atom_site_label", charge_column_name])
        if len(table) == 0:
            continue
        charges_by_label: dict[str, float] = {}
        charges_in_row_order: list[float] = []
        duplicate_labels = False
        for row_index in range(len(table)):
            row = table[row_index]
            label = str(row[0]).strip()
            if not label:
                raise LammpsInputError(f"CIF charge loop {charge_column_name} contains an empty atom label.")
            try:
                charge = float(row[1])
            except ValueError as exc:
                raise LammpsInputError(
                    f"CIF charge loop {charge_column_name} contains a non-numeric value for atom {label!r}: {row[1]!r}"
                ) from exc
            if label in charges_by_label:
                duplicate_labels = True
            charges_by_label[label] = charge
            charges_in_row_order.append(charge)
        if expected_label_set:
            missing = sorted(expected_label_set - charges_by_label.keys())
            if missing:
                matching_labels = expected_label_set.intersection(charges_by_label)
                if len(charges_in_row_order) == len(expected_label_sequence) and (
                    duplicate_labels or not matching_labels
                ):
                    return _AtomChargeExtractionResult(
                        charges_by_label={
                            label: charges_in_row_order[index]
                            for index, label in enumerate(expected_label_sequence)
                        },
                        used_atom_order_fallback=True,
                    )
                raise LammpsInputError(
                    f"CIF charge loop {charge_column_name} is missing charges for atom labels: {', '.join(missing)}"
                )
        return _AtomChargeExtractionResult(charges_by_label=charges_by_label)
    return _AtomChargeExtractionResult(charges_by_label={})


def _with_atom_charges(
    parsed: _ParsedExplicitBondCif,
    charge_by_label: dict[str, float],
    *,
    source_path: Path,
) -> _ParsedExplicitBondCif:
    missing = [atom.label for atom in parsed.atoms if atom.label not in charge_by_label]
    if missing:
        raise LammpsInputError(
            f"Charge source {source_path} is missing atom charges for labels: {', '.join(sorted(missing))}"
        )
    atoms = tuple(
        replace(atom, charge=float(charge_by_label[atom.label]))
        for atom in parsed.atoms
    )
    return replace(parsed, atoms=atoms)


def _derive_angles(
    atoms: tuple[_CifAtomRecord, ...],
    bonds: tuple[_CifBondRecord, ...],
    basis: tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]],
) -> tuple[_CifAngleRecord, ...]:
    positions = {atom.atom_id: atom.fractional_position for atom in atoms}
    neighbors: dict[int, list[tuple[int, tuple[int, int, int]]]] = {atom.atom_id: [] for atom in atoms}
    for bond in bonds:
        offset_1_to_2 = _sub_shift(bond.shift_2, bond.shift_1)
        neighbors[bond.atom_id_1].append((bond.atom_id_2, offset_1_to_2))
        neighbors[bond.atom_id_2].append((bond.atom_id_1, _neg_shift(offset_1_to_2)))

    angles: list[_CifAngleRecord] = []
    seen: set[tuple[int, int, int]] = set()
    for center_id, entries in neighbors.items():
        ordered = sorted(entries, key=lambda item: item[0])
        for left_index in range(len(ordered)):
            for right_index in range(left_index + 1, len(ordered)):
                left_atom_id, left_shift = ordered[left_index]
                right_atom_id, right_shift = ordered[right_index]
                key = (left_atom_id, center_id, right_atom_id)
                if key in seen:
                    continue
                seen.add(key)
                center = positions[center_id]
                left_vector = _fractional_vector(
                    positions[left_atom_id],
                    center,
                    left_shift,
                )
                right_vector = _fractional_vector(
                    positions[right_atom_id],
                    center,
                    right_shift,
                )
                left_cart = _fractional_vector_to_cartesian(left_vector, basis)
                right_cart = _fractional_vector_to_cartesian(right_vector, basis)
                angle = _angle_degrees(left_cart, right_cart)
                angles.append(
                    _CifAngleRecord(
                        angle_id=len(angles) + 1,
                        atom_id_1=left_atom_id,
                        atom_id_2=center_id,
                        atom_id_3=right_atom_id,
                        equilibrium_degrees=angle,
                    )
                )
    return tuple(angles)


def _derive_dihedrals(
    atoms: tuple[_CifAtomRecord, ...],
    bonds: tuple[_CifBondRecord, ...],
) -> tuple[_CifDihedralRecord, ...]:
    neighbors: dict[int, list[int]] = {atom.atom_id: [] for atom in atoms}
    for bond in bonds:
        neighbors[bond.atom_id_1].append(bond.atom_id_2)
        neighbors[bond.atom_id_2].append(bond.atom_id_1)

    dihedrals: list[_CifDihedralRecord] = []
    seen: set[tuple[int, int, int, int]] = set()
    for bond in bonds:
        left_center = bond.atom_id_1
        right_center = bond.atom_id_2
        left_neighbors = [atom_id for atom_id in neighbors[left_center] if atom_id != right_center]
        right_neighbors = [atom_id for atom_id in neighbors[right_center] if atom_id != left_center]
        for left_atom_id in sorted(left_neighbors):
            for right_atom_id in sorted(right_neighbors):
                if left_atom_id == right_atom_id:
                    continue
                forward = (left_atom_id, left_center, right_center, right_atom_id)
                reverse = tuple(reversed(forward))
                key = min(forward, reverse)
                if key in seen:
                    continue
                seen.add(key)
                dihedrals.append(
                    _CifDihedralRecord(
                        dihedral_id=len(dihedrals) + 1,
                        atom_id_1=forward[0],
                        atom_id_2=forward[1],
                        atom_id_3=forward[2],
                        atom_id_4=forward[3],
                    )
                )
    return tuple(dihedrals)


def _derive_impropers(
    atoms: tuple[_CifAtomRecord, ...],
    bonds: tuple[_CifBondRecord, ...],
) -> tuple[_CifImproperRecord, ...]:
    neighbors: dict[int, list[int]] = {atom.atom_id: [] for atom in atoms}
    for bond in bonds:
        neighbors[bond.atom_id_1].append(bond.atom_id_2)
        neighbors[bond.atom_id_2].append(bond.atom_id_1)

    impropers: list[_CifImproperRecord] = []
    for center_id, entries in sorted(neighbors.items()):
        ordered = sorted(entries)
        if len(ordered) != 3:
            continue
        permutations = [
            (ordered[0], ordered[1], ordered[2]),
            (ordered[1], ordered[0], ordered[2]),
            (ordered[2], ordered[0], ordered[1]),
        ]
        for atom_id_1, atom_id_3, atom_id_4 in permutations:
            impropers.append(
                _CifImproperRecord(
                    improper_id=len(impropers) + 1,
                    atom_id_1=atom_id_1,
                    atom_id_2=center_id,
                    atom_id_3=atom_id_3,
                    atom_id_4=atom_id_4,
                )
            )
    return tuple(impropers)


def _resolve_effective_bond_geometry(
    *,
    left_fractional: tuple[float, float, float],
    right_fractional: tuple[float, float, float],
    raw_shift_1: tuple[int, int, int],
    raw_shift_2: tuple[int, int, int],
    raw_symmetry_1: str,
    raw_symmetry_2: str,
    basis: tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]],
) -> tuple[tuple[int, int, int], tuple[int, int, int], float]:
    if raw_symmetry_1 != "." or raw_symmetry_2 != ".":
        equilibrium_distance = _bond_distance(
            left_fractional,
            right_fractional,
            raw_shift_1,
            raw_shift_2,
            basis,
        )
        return raw_shift_1, raw_shift_2, equilibrium_distance

    _shift_1, best_shift_2 = _closest_periodic_shifts(
        left_fractional=left_fractional,
        right_fractional=right_fractional,
        basis=basis,
    )
    best_distance = _bond_distance(
        left_fractional,
        right_fractional,
        (0, 0, 0),
        best_shift_2,
        basis,
    )
    return (0, 0, 0), best_shift_2, best_distance


def _closest_periodic_shifts(
    *,
    left_fractional: tuple[float, float, float],
    right_fractional: tuple[float, float, float],
    basis: tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]],
) -> tuple[tuple[int, int, int], tuple[int, int, int]]:
    best_shift_2 = (0, 0, 0)
    best_distance = float("inf")
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for dz in (-1, 0, 1):
                trial_shift_2 = (dx, dy, dz)
                distance = _bond_distance(
                    left_fractional,
                    right_fractional,
                    (0, 0, 0),
                    trial_shift_2,
                    basis,
                )
                if distance + 1e-9 < best_distance:
                    best_distance = distance
                    best_shift_2 = trial_shift_2
                elif abs(distance - best_distance) <= 1e-9 and trial_shift_2 == (0, 0, 0):
                    best_shift_2 = trial_shift_2
    return (0, 0, 0), best_shift_2


def _default_stage2_restraint_force_constant(settings: LammpsOptimizationSettings) -> float:
    return 0.0


def _build_minimization_stages(
    settings: LammpsOptimizationSettings,
) -> list[_LammpsMinimizationStage]:
    stages = [
        _LammpsMinimizationStage(
            label="stage1",
            min_style=settings.min_style,
            energy_tolerance=settings.energy_tolerance,
            force_tolerance=settings.force_tolerance,
            max_iterations=settings.max_iterations,
            max_evaluations=settings.max_evaluations,
            position_restraint_force_constant=settings.position_restraint_force_constant,
        )
    ]
    if settings.two_stage_protocol:
        stages.append(
            _LammpsMinimizationStage(
                label="stage2",
                min_style=settings.stage2_min_style or settings.min_style,
                energy_tolerance=settings.stage2_energy_tolerance or settings.energy_tolerance,
                force_tolerance=settings.stage2_force_tolerance or settings.force_tolerance,
                max_iterations=settings.stage2_max_iterations or settings.max_iterations,
                max_evaluations=settings.stage2_max_evaluations or settings.max_evaluations,
                position_restraint_force_constant=(
                    settings.stage2_position_restraint_force_constant
                    if settings.stage2_position_restraint_force_constant is not None
                    else _default_stage2_restraint_force_constant(settings)
                ),
            )
        )
    if settings.relax_cell:
        last_stage = stages[-1]
        stages.append(
            _LammpsMinimizationStage(
                label="box_relax",
                min_style=settings.box_relax_min_style,
                energy_tolerance=settings.box_relax_energy_tolerance or last_stage.energy_tolerance,
                force_tolerance=settings.box_relax_force_tolerance or last_stage.force_tolerance,
                max_iterations=settings.box_relax_max_iterations or last_stage.max_iterations,
                max_evaluations=settings.box_relax_max_evaluations or last_stage.max_evaluations,
                position_restraint_force_constant=last_stage.position_restraint_force_constant,
                relax_cell=True,
            )
        )
    return stages


def _render_min_modify_line(settings: LammpsOptimizationSettings, *, min_style: str) -> str | None:
    tokens: list[str] = []
    if settings.min_modify_dmax is not None:
        tokens.extend(["dmax", f"{settings.min_modify_dmax:.8g}"])
    if settings.min_modify_line is not None:
        tokens.extend(["line", settings.min_modify_line])
    if settings.min_modify_norm is not None:
        tokens.extend(["norm", settings.min_modify_norm])
    if _is_fire_min_style(min_style):
        if settings.min_modify_fire_integrator is not None:
            tokens.extend(["integrator", settings.min_modify_fire_integrator])
        if settings.min_modify_fire_tmax is not None:
            tokens.extend(["tmax", f"{settings.min_modify_fire_tmax:.8g}"])
        if settings.min_modify_fire_abcfire is not None:
            tokens.extend(["abcfire", "yes" if settings.min_modify_fire_abcfire else "no"])
    if not tokens:
        return None
    return "min_modify " + " ".join(tokens)


def _box_relax_mode_for_structure(
    settings: LammpsOptimizationSettings,
    parsed: _ParsedExplicitBondCif,
) -> str:
    if settings.box_relax_mode != "auto":
        return settings.box_relax_mode
    alpha, beta, gamma = parsed.cell_parameters[3:]
    right_angle_tolerance = 1.0e-5
    if (
        abs(alpha - 90.0) <= right_angle_tolerance
        and abs(beta - 90.0) <= right_angle_tolerance
        and abs(gamma - 90.0) <= right_angle_tolerance
    ):
        return "aniso"
    return "tri"


def _render_box_relax_fix_line(
    settings: LammpsOptimizationSettings,
    parsed: _ParsedExplicitBondCif,
) -> str:
    mode = _box_relax_mode_for_structure(settings, parsed)
    tokens = [
        "fix",
        "cofkit_boxrelax",
        "all",
        "box/relax",
        mode,
        f"{settings.box_relax_target_pressure:.8g}",
        "vmax",
        f"{settings.box_relax_vmax:.8g}",
    ]
    if settings.box_relax_nreset is not None:
        tokens.extend(["nreset", str(settings.box_relax_nreset)])
    return " ".join(tokens)


def _has_nonperiodic_boundary(boundary: str) -> bool:
    return any(token.lower() == "f" for token in boundary.split())


def _render_lammps_input_script(
    *,
    data_path: Path,
    dump_path: Path,
    settings: LammpsOptimizationSettings,
    parsed: _ParsedExplicitBondCif,
    prepared: _PreparedLammpsSystem,
    boundary: str,
) -> str:
    thermo_terms = ["step", "pe", "ebond"]
    minimization_stages = _build_minimization_stages(settings)
    periodic_electrostatics = prepared.has_charges and not _has_nonperiodic_boundary(boundary)
    pair_style_line = (
        f"pair_style lj/cut/coul/long {settings.pair_cutoff:.6f} {settings.coulomb_cutoff:.6f}"
        if periodic_electrostatics
        else (
            f"pair_style lj/cut/coul/cut {settings.pair_cutoff:.6f} {settings.coulomb_cutoff:.6f}"
            if prepared.has_charges
            else f"pair_style lj/cut {settings.pair_cutoff:.6f}"
        )
    )
    special_bonds_line = (
        "special_bonds lj/coul 0.0 0.0 1.0" if prepared.has_charges else "special_bonds lj 0.0 0.0 1.0"
    )
    lines = [
        "units real",
        "atom_style full" if prepared.has_charges else "atom_style molecular",
        f"boundary {boundary}",
        pair_style_line,
        "pair_modify shift yes",
        special_bonds_line,
        "bond_style harmonic",
    ]
    if prepared.angle_styles:
        lines.append(_render_style_line("angle_style", prepared.angle_styles, force_hybrid=True))
        thermo_terms.append("eangle")
    if prepared.dihedral_styles:
        lines.append(_render_style_line("dihedral_style", prepared.dihedral_styles))
        thermo_terms.append("edihed")
    if prepared.improper_styles:
        lines.append(_render_style_line("improper_style", prepared.improper_styles))
        thermo_terms.append("eimp")
    thermo_terms.append("evdwl")
    if prepared.has_charges:
        thermo_terms.append("ecoul")
        if periodic_electrostatics:
            thermo_terms.append("elong")
    thermo_terms.append("press")
    lines.extend(
        [
            f"read_data {data_path}",
            *(["kspace_style ewald " + f"{settings.ewald_precision:.8g}"] if periodic_electrostatics else []),
            "neighbor 2.0 bin",
            "neigh_modify every 1 delay 0 check yes",
        ]
    )
    if settings.timestep is not None:
        lines.append(f"timestep {settings.timestep:.8g}")
    lines.extend(
        [
            "thermo 50",
            "thermo_style custom " + " ".join(thermo_terms),
            f"dump cofkit_dump all custom 1 {dump_path} id x y z",
            "dump_modify cofkit_dump sort id",
        ]
    )
    active_restraint_force_constant: float | None = None
    if settings.pre_minimization_steps > 0:
        initial_restraint_force_constant = minimization_stages[0].position_restraint_force_constant
        lines.append("# pre_minimization")
        if initial_restraint_force_constant > 0.0:
            lines.append(f"fix cofkit_hold all spring/self {initial_restraint_force_constant:.6f}")
            lines.append("fix_modify cofkit_hold energy yes")
            active_restraint_force_constant = initial_restraint_force_constant
        lines.extend(
            [
                (
                    f"velocity all create {settings.pre_minimization_temperature:.8g} "
                    f"{settings.pre_minimization_seed} mom yes rot yes dist gaussian"
                ),
                (
                    f"fix cofkit_prerun_nve all nve/limit "
                    f"{settings.pre_minimization_displacement_limit:.8g}"
                ),
                (
                    f"fix cofkit_prerun_langevin all langevin "
                    f"{settings.pre_minimization_temperature:.8g} "
                    f"{settings.pre_minimization_temperature:.8g} "
                    f"{settings.pre_minimization_damping:.8g} "
                    f"{settings.pre_minimization_seed + 104729}"
                ),
                f"run {settings.pre_minimization_steps}",
                "unfix cofkit_prerun_langevin",
                "unfix cofkit_prerun_nve",
            ]
        )
    for stage in minimization_stages:
        lines.append(f"# {stage.label}")
        if active_restraint_force_constant is not None and (
            abs(active_restraint_force_constant - stage.position_restraint_force_constant) > 1.0e-12
        ):
            lines.append("unfix cofkit_hold")
            active_restraint_force_constant = None
        if stage.position_restraint_force_constant > 0.0 and active_restraint_force_constant is None:
            lines.append(f"fix cofkit_hold all spring/self {stage.position_restraint_force_constant:.6f}")
            lines.append("fix_modify cofkit_hold energy yes")
            active_restraint_force_constant = stage.position_restraint_force_constant
        if stage.relax_cell:
            lines.append(_render_box_relax_fix_line(settings, parsed))
        lines.append(f"min_style {stage.min_style}")
        min_modify_line = _render_min_modify_line(settings, min_style=stage.min_style)
        if min_modify_line is not None:
            lines.append(min_modify_line)
        lines.append(
            f"minimize {stage.energy_tolerance:.8g} {stage.force_tolerance:.8g} "
            f"{stage.max_iterations} {stage.max_evaluations}"
        )
        if stage.relax_cell:
            lines.append("unfix cofkit_boxrelax")
    lines.append("undump cofkit_dump")
    if active_restraint_force_constant is not None:
        lines.append("unfix cofkit_hold")
    lines.append("")
    return "\n".join(lines)


def _run_lammps(
    *,
    binary: Path,
    input_script_path: Path,
    log_path: Path,
    stdout_log_path: Path,
    stderr_log_path: Path,
    timeout_seconds: float,
) -> tuple[str, ...]:
    command = [
        str(binary),
        "-in",
        str(input_script_path),
        "-log",
        str(log_path),
        "-screen",
        "none",
    ]
    environment = os.environ.copy()
    if not str(environment.get("OMP_NUM_THREADS", "")).strip():
        environment["OMP_NUM_THREADS"] = str(_default_lammps_omp_num_threads())
    try:
        completed = subprocess.run(
            command,
            check=False,
            capture_output=True,
            text=True,
            timeout=timeout_seconds,
            env=environment,
        )
    except subprocess.TimeoutExpired as exc:
        stdout = exc.stdout if isinstance(exc.stdout, str) else (exc.stdout or b"").decode("utf-8", errors="replace")
        stderr = exc.stderr if isinstance(exc.stderr, str) else (exc.stderr or b"").decode("utf-8", errors="replace")
        stdout_log_path.write_text(stdout, encoding="utf-8")
        stderr_log_path.write_text(stderr, encoding="utf-8")
        raise LammpsExecutionError(
            f"LAMMPS command timed out after {timeout_seconds:g} seconds: {' '.join(command)}"
        ) from exc

    stdout_log_path.write_text(completed.stdout or "", encoding="utf-8")
    stderr_log_path.write_text(completed.stderr or "", encoding="utf-8")
    if completed.returncode != 0:
        raise LammpsExecutionError(
            f"LAMMPS command failed with exit code {completed.returncode}: {' '.join(command)}"
        )

    warnings: list[str] = []
    merged = "\n".join(filter(None, (completed.stdout, completed.stderr)))
    for line in merged.splitlines():
        stripped = line.strip()
        if stripped.startswith("WARNING:"):
            warnings.append(stripped)
    return tuple(warnings)


def _default_lammps_omp_num_threads() -> int:
    cpu_count = os.cpu_count() or 1
    return max(1, cpu_count // 2)


def _render_style_line(keyword: str, styles: Sequence[str], *, force_hybrid: bool = False) -> str:
    unique_styles = tuple(dict.fromkeys(styles))
    if len(unique_styles) == 1 and not force_hybrid:
        return f"{keyword} {unique_styles[0]}"
    return f"{keyword} hybrid {' '.join(unique_styles)}"


def _parse_lammps_dump_last_frame(
    dump_path: Path,
    *,
    expected_atoms: int,
) -> _LammpsDumpFrame:
    if not dump_path.is_file():
        raise LammpsParseError(f"LAMMPS dump file was not created: {dump_path}")

    lines = dump_path.read_text(encoding="utf-8").splitlines()
    index = 0
    last_frame: dict[int, tuple[float, float, float]] | None = None
    last_origin: tuple[float, float, float] | None = None
    last_basis: tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]] | None = None
    last_cell_parameters: tuple[float, float, float, float, float, float] | None = None
    while index < len(lines):
        if lines[index].strip() != "ITEM: TIMESTEP":
            raise LammpsParseError(f"Unexpected LAMMPS dump format in {dump_path}: missing ITEM: TIMESTEP")
        index += 1
        if index >= len(lines):
            break
        index += 1
        if index >= len(lines) or lines[index].strip() != "ITEM: NUMBER OF ATOMS":
            raise LammpsParseError(f"Unexpected LAMMPS dump format in {dump_path}: missing atom count header")
        index += 1
        if index >= len(lines):
            raise LammpsParseError(f"Unexpected end of LAMMPS dump in {dump_path}")
        try:
            n_atoms = int(lines[index].strip())
        except ValueError as exc:
            raise LammpsParseError(f"Invalid atom count in LAMMPS dump {dump_path}: {lines[index]!r}") from exc
        index += 1
        if index >= len(lines) or not lines[index].startswith("ITEM: BOX BOUNDS"):
            raise LammpsParseError(f"Unexpected LAMMPS dump format in {dump_path}: missing box-bounds header")
        origin, basis, cell_parameters = _parse_lammps_dump_box(lines[index : index + 4], dump_path=dump_path)
        index += 4
        if index > len(lines) or not lines[index - 1]:
            pass
        if index >= len(lines) or not lines[index].startswith("ITEM: ATOMS"):
            raise LammpsParseError(f"Unexpected LAMMPS dump format in {dump_path}: missing atom header")
        header = lines[index].split()[2:]
        index += 1
        expected_header = ["id", "x", "y", "z"]
        if header != expected_header:
            raise LammpsParseError(
                f"Unexpected LAMMPS dump atom columns in {dump_path}: expected {expected_header}, got {header}"
            )
        frame: dict[int, tuple[float, float, float]] = {}
        for _ in range(n_atoms):
            if index >= len(lines):
                raise LammpsParseError(f"Unexpected end of LAMMPS dump in {dump_path}")
            parts = lines[index].split()
            index += 1
            if len(parts) != 4:
                raise LammpsParseError(f"Unexpected LAMMPS dump atom row in {dump_path}: {lines[index - 1]!r}")
            atom_id = int(parts[0])
            frame[atom_id] = (float(parts[1]), float(parts[2]), float(parts[3]))
        last_frame = frame
        last_origin = origin
        last_basis = basis
        last_cell_parameters = cell_parameters

    if last_frame is None or last_origin is None or last_basis is None or last_cell_parameters is None:
        raise LammpsParseError(f"LAMMPS dump file did not contain any frames: {dump_path}")
    if len(last_frame) != expected_atoms:
        raise LammpsParseError(
            f"LAMMPS dump file {dump_path} contained {len(last_frame)} atoms in the final frame, expected {expected_atoms}"
        )
    return _LammpsDumpFrame(
        cartesian_positions=last_frame,
        lammps_origin=last_origin,
        lammps_basis=last_basis,
        cell_parameters=last_cell_parameters,
    )


def _cartesian_positions_to_fractional(
    origin: tuple[float, float, float],
    basis: tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]],
    final_cartesian_positions: dict[int, tuple[float, float, float]],
) -> dict[int, tuple[float, float, float]]:
    a_vec, b_vec, c_vec = basis
    ax = a_vec[0]
    by = b_vec[1]
    xy = b_vec[0]
    xz = c_vec[0]
    yz = c_vec[1]
    cz = c_vec[2]
    origin_x, origin_y, origin_z = origin

    positions: dict[int, tuple[float, float, float]] = {}
    for atom_id, (x, y, z) in final_cartesian_positions.items():
        if abs(cz) < 1e-12 or abs(by) < 1e-12 or abs(ax) < 1e-12:
            raise LammpsParseError("LAMMPS cell basis is singular and cannot be converted back to fractional coordinates.")
        relative_x = x - origin_x
        relative_y = y - origin_y
        relative_z = z - origin_z
        fz = relative_z / cz
        fy = (relative_y - fz * yz) / by
        fx = (relative_x - fy * xy - fz * xz) / ax
        positions[atom_id] = (_wrap_fraction(fx), _wrap_fraction(fy), _wrap_fraction(fz))
    return positions


def _render_optimized_cif(
    parsed: _ParsedExplicitBondCif,
    final_fractional_positions: dict[int, tuple[float, float, float]],
    *,
    cell_parameters: tuple[float, float, float, float, float, float] | None = None,
    basis: tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]] | None = None,
    cofid: str | None = None,
    cofid_comment_suffix: str | None = None,
) -> str:
    active_cell_parameters = cell_parameters or parsed.cell_parameters
    active_basis = basis or parsed.lammps_basis
    a, b, c, alpha, beta, gamma = active_cell_parameters
    has_charges = any(atom.charge is not None for atom in parsed.atoms)
    lines: list[str] = []
    if cofid:
        lines.append(cofid_comment_line(cofid, suffix=cofid_comment_suffix))
    lines.extend(
        [
            f"data_{_sanitize_data_name(parsed.source_path.stem + '_lammps_optimized')}",
            "# CIF generated by cofkit LAMMPS optimization wrapper",
            "_audit_creation_method 'cofkit LAMMPS optimizer'",
            "_space_group_name_H-M_alt 'P 1'",
            "_space_group_IT_number 1",
            f"_cell_length_a {a:.6f}",
            f"_cell_length_b {b:.6f}",
            f"_cell_length_c {c:.6f}",
            f"_cell_angle_alpha {alpha:.6f}",
            f"_cell_angle_beta {beta:.6f}",
            f"_cell_angle_gamma {gamma:.6f}",
            "",
            "loop_",
            "_space_group_symop_operation_xyz",
            "'x,y,z'",
            "",
            "loop_",
            "_atom_site_label",
            "_atom_site_type_symbol",
            "_atom_site_fract_x",
            "_atom_site_fract_y",
            "_atom_site_fract_z",
            "_atom_site_occupancy",
        ]
    )
    if has_charges:
        lines.append("_atom_site_charge")
    for atom in parsed.atoms:
        fractional = final_fractional_positions[atom.atom_id]
        line = (
            f"{atom.label} {atom.symbol} {fractional[0]:.6f} {fractional[1]:.6f} "
            f"{fractional[2]:.6f} {atom.occupancy:.2f}"
        )
        if has_charges:
            line += f" {float(atom.charge):.6f}"
        lines.append(line)

    lines.extend(
        [
            "",
            "loop_",
            "_geom_bond_atom_site_label_1",
            "_geom_bond_atom_site_label_2",
            "_geom_bond_site_symmetry_1",
            "_geom_bond_site_symmetry_2",
            "_geom_bond_distance",
            "_ccdc_geom_bond_type",
        ]
    )
    for bond in parsed.bonds:
        fractional_1 = final_fractional_positions[bond.atom_id_1]
        fractional_2 = final_fractional_positions[bond.atom_id_2]
        output_shift_1, output_shift_2 = _closest_periodic_shifts(
            left_fractional=fractional_1,
            right_fractional=fractional_2,
            basis=active_basis,
        )
        symmetry_1 = _format_p1_symmetry_shift(output_shift_1)
        symmetry_2 = _format_p1_symmetry_shift(output_shift_2)
        distance = _bond_distance(
            fractional_1,
            fractional_2,
            output_shift_1,
            output_shift_2,
            active_basis,
        )
        lines.append(
            f"{bond.label_1} {bond.label_2} {symmetry_1} {symmetry_2} {distance:.6f} "
            f"{bond_order_to_cif_type(bond.bond_order if bond.bond_order is not None else 1.0)}"
        )
    lines.append("")
    return "\n".join(lines)


def _resolve_output_dir(input_path: Path, output_dir: str | Path | None) -> Path:
    if output_dir is None:
        return input_path.parent / f"{input_path.stem}_lammps"
    return Path(output_dir).expanduser().resolve()


def _lammps_basis_from_cell(
    a: float,
    b: float,
    c: float,
    alpha_degrees: float,
    beta_degrees: float,
    gamma_degrees: float,
) -> tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]]:
    alpha = math.radians(alpha_degrees)
    beta = math.radians(beta_degrees)
    gamma = math.radians(gamma_degrees)
    xy = b * math.cos(gamma)
    xz = c * math.cos(beta)
    by = b * math.sin(gamma)
    if abs(by) < 1e-12:
        raise LammpsInputError("CIF cell gamma angle is degenerate and cannot be converted into a LAMMPS triclinic box.")
    yz = (b * c * math.cos(alpha) - xy * xz) / by
    cz_sq = max(c * c - xz * xz - yz * yz, 0.0)
    cz = math.sqrt(cz_sq)
    return (
        (a, 0.0, 0.0),
        (xy, by, 0.0),
        (xz, yz, cz),
    )


def _cell_parameters_from_lammps_basis(
    basis: tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]],
) -> tuple[float, float, float, float, float, float]:
    a_vec, b_vec, c_vec = basis
    a = _norm(a_vec)
    b = _norm(b_vec)
    c = _norm(c_vec)
    if a < 1e-12 or b < 1e-12 or c < 1e-12:
        raise LammpsParseError("LAMMPS dump box basis is singular and cannot be converted into CIF cell parameters.")
    alpha = math.degrees(math.acos(max(-1.0, min(1.0, _dot(b_vec, c_vec) / (b * c)))))
    beta = math.degrees(math.acos(max(-1.0, min(1.0, _dot(a_vec, c_vec) / (a * c)))))
    gamma = math.degrees(math.acos(max(-1.0, min(1.0, _dot(a_vec, b_vec) / (a * b)))))
    return (a, b, c, alpha, beta, gamma)


def _parse_lammps_dump_box(
    lines: Sequence[str],
    *,
    dump_path: Path,
) -> tuple[
    tuple[float, float, float],
    tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]],
    tuple[float, float, float, float, float, float],
]:
    if len(lines) != 4:
        raise LammpsParseError(f"Unexpected LAMMPS dump box block in {dump_path}: expected 4 lines, got {len(lines)}")
    header = lines[0].split()
    if header[:3] != ["ITEM:", "BOX", "BOUNDS"]:
        raise LammpsParseError(f"Unexpected LAMMPS dump box header in {dump_path}: {lines[0]!r}")
    triclinic = "xy" in header and "xz" in header and "yz" in header
    try:
        bounds = [tuple(float(value) for value in line.split()) for line in lines[1:4]]
    except ValueError as exc:
        raise LammpsParseError(f"Invalid numeric box-bounds line in {dump_path}") from exc
    if triclinic:
        if any(len(row) != 3 for row in bounds):
            raise LammpsParseError(f"Unexpected triclinic box-bounds format in {dump_path}")
        (xlo_bound, xhi_bound, xy), (ylo_bound, yhi_bound, xz), (zlo_bound, zhi_bound, yz) = bounds
        xlo = xlo_bound - min(0.0, xy, xz, xy + xz)
        xhi = xhi_bound - max(0.0, xy, xz, xy + xz)
        ylo = ylo_bound - min(0.0, yz)
        yhi = yhi_bound - max(0.0, yz)
        zlo = zlo_bound
        zhi = zhi_bound
        basis = (
            (xhi - xlo, 0.0, 0.0),
            (xy, yhi - ylo, 0.0),
            (xz, yz, zhi - zlo),
        )
        origin = (xlo, ylo, zlo)
    else:
        if any(len(row) != 2 for row in bounds):
            raise LammpsParseError(f"Unexpected orthogonal box-bounds format in {dump_path}")
        (xlo, xhi), (ylo, yhi), (zlo, zhi) = bounds
        basis = (
            (xhi - xlo, 0.0, 0.0),
            (0.0, yhi - ylo, 0.0),
            (0.0, 0.0, zhi - zlo),
        )
        origin = (xlo, ylo, zlo)
    return origin, basis, _cell_parameters_from_lammps_basis(basis)


def _fractional_to_lammps(
    fractional_position: tuple[float, float, float],
    basis: tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]],
) -> tuple[float, float, float]:
    fx, fy, fz = fractional_position
    a_vec, b_vec, c_vec = basis
    return (
        fx * a_vec[0] + fy * b_vec[0] + fz * c_vec[0],
        fx * a_vec[1] + fy * b_vec[1] + fz * c_vec[1],
        fx * a_vec[2] + fy * b_vec[2] + fz * c_vec[2],
    )


def _fractional_vector(
    left: tuple[float, float, float],
    center: tuple[float, float, float],
    shift: tuple[int, int, int],
) -> tuple[float, float, float]:
    return (
        left[0] + shift[0] - center[0],
        left[1] + shift[1] - center[1],
        left[2] + shift[2] - center[2],
    )


def _fractional_vector_to_cartesian(
    fractional_vector: tuple[float, float, float],
    basis: tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]],
) -> tuple[float, float, float]:
    return _fractional_to_lammps(fractional_vector, basis)


def _angle_degrees(left: tuple[float, float, float], right: tuple[float, float, float]) -> float:
    left_norm = _norm(left)
    right_norm = _norm(right)
    if left_norm < 1e-12 or right_norm < 1e-12:
        return 180.0
    cosine = _dot(left, right) / (left_norm * right_norm)
    cosine = max(-1.0, min(1.0, cosine))
    return math.degrees(math.acos(cosine))


def _bond_distance(
    fractional_1: tuple[float, float, float],
    fractional_2: tuple[float, float, float],
    shift_1: tuple[int, int, int],
    shift_2: tuple[int, int, int],
    basis: tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]],
) -> float:
    vector = (
        (fractional_2[0] + shift_2[0]) - (fractional_1[0] + shift_1[0]),
        (fractional_2[1] + shift_2[1]) - (fractional_1[1] + shift_1[1]),
        (fractional_2[2] + shift_2[2]) - (fractional_1[2] + shift_1[2]),
    )
    return _norm(_fractional_vector_to_cartesian(vector, basis))

def _add_shift(left: tuple[int, int, int], right: tuple[int, int, int]) -> tuple[int, int, int]:
    return (left[0] + right[0], left[1] + right[1], left[2] + right[2])


def _sub_shift(left: tuple[int, int, int], right: tuple[int, int, int]) -> tuple[int, int, int]:
    return (left[0] - right[0], left[1] - right[1], left[2] - right[2])


def _neg_shift(shift: tuple[int, int, int]) -> tuple[int, int, int]:
    return (-shift[0], -shift[1], -shift[2])


def _dot(left: tuple[float, float, float], right: tuple[float, float, float]) -> float:
    return left[0] * right[0] + left[1] * right[1] + left[2] * right[2]


def _norm(vector: tuple[float, float, float]) -> float:
    return math.sqrt(_dot(vector, vector))


def _wrap_fraction(value: float) -> float:
    wrapped = value % 1.0
    if wrapped < 0.0:
        wrapped += 1.0
    if wrapped >= 1.0 - 1e-12:
        return 0.0
    return wrapped


def _format_p1_symmetry_shift(shift: tuple[int, int, int]) -> str:
    if shift == (0, 0, 0):
        return "."
    return "1_" + "".join(str(component + 5) for component in shift)


def _sanitize_data_name(value: str) -> str:
    sanitized = []
    for char in value:
        if char.isalnum() or char in {"_", "-"}:
            sanitized.append(char)
        else:
            sanitized.append("_")
    collapsed = "".join(sanitized).strip("_")
    return collapsed or "cofkit_lammps"
