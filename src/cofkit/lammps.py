from __future__ import annotations

import itertools
import json
import math
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

try:
    import gemmi
except ImportError:  # pragma: no cover - exercised in environments without gemmi
    gemmi = None

try:
    from openbabel import openbabel as ob
except ImportError:  # pragma: no cover - exercised in environments without Open Babel
    ob = None

try:
    from pymatgen.core import Lattice, Structure
    from pymatgen.io.lammps.data import ForceField, LammpsBox, LammpsData, Topology
except ImportError:  # pragma: no cover - exercised in environments without pymatgen
    ForceField = None
    LammpsBox = None
    LammpsData = None
    Lattice = None
    Structure = None
    Topology = None

try:
    from rdkit import Chem
    from rdkit.Chem import rdForceFieldHelpers
except ImportError:  # pragma: no cover - exercised in environments without RDKit
    Chem = None
    rdForceFieldHelpers = None


COFKIT_LMP_ENV_VAR = "COFKIT_LMP_PATH"
DEFAULT_LAMMPS_BINARY = Path("/opt/homebrew/bin/lmp_mpi")
_SUPPORTED_FORCEFIELDS = ("uff",)
_UFF_RMIN_TO_SIGMA_FACTOR = 2.0 ** (1.0 / 6.0)
_MIN_MODIFY_LINE_OPTIONS = {"backtrack", "quadratic", "forcezero", "spin_cubic", "spin_none"}
_MIN_MODIFY_NORM_OPTIONS = {"two", "inf", "max"}
_MIN_MODIFY_FIRE_INTEGRATOR_OPTIONS = {"eulerimplicit", "verlet", "leapfrog", "eulerexplicit"}
_BOX_RELAX_MODE_OPTIONS = {"auto", "iso", "aniso", "tri"}
_UNSUPPORTED_BOX_RELAX_MIN_STYLES = {"quickmin", "fire", "hftn", "cg/kk"}


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
    pair_cutoff: float = 12.0
    position_restraint_force_constant: float = 0.20
    two_stage_protocol: bool = False
    stage2_position_restraint_force_constant: float | None = None
    energy_tolerance: float = 1.0e-6
    force_tolerance: float = 1.0e-6
    max_iterations: int = 500
    max_evaluations: int = 5000
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
    relax_cell: bool = False
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
        return {
            "forcefield": self.forcefield,
            "pair_cutoff": self.pair_cutoff,
            "position_restraint_force_constant": self.position_restraint_force_constant,
            "two_stage_protocol": self.two_stage_protocol,
            "stage2_position_restraint_force_constant": self.stage2_position_restraint_force_constant,
            "energy_tolerance": self.energy_tolerance,
            "force_tolerance": self.force_tolerance,
            "max_iterations": self.max_iterations,
            "max_evaluations": self.max_evaluations,
            "min_style": self.min_style,
            "stage2_energy_tolerance": self.stage2_energy_tolerance,
            "stage2_force_tolerance": self.stage2_force_tolerance,
            "stage2_max_iterations": self.stage2_max_iterations,
            "stage2_max_evaluations": self.stage2_max_evaluations,
            "stage2_min_style": self.stage2_min_style,
            "timestep": self.timestep,
            "min_modify_dmax": self.min_modify_dmax,
            "min_modify_line": self.min_modify_line,
            "min_modify_norm": self.min_modify_norm,
            "min_modify_fire_integrator": self.min_modify_fire_integrator,
            "min_modify_fire_tmax": self.min_modify_fire_tmax,
            "min_modify_fire_abcfire": self.min_modify_fire_abcfire,
            "relax_cell": self.relax_cell,
            "box_relax_mode": self.box_relax_mode,
            "box_relax_target_pressure": self.box_relax_target_pressure,
            "box_relax_vmax": self.box_relax_vmax,
            "box_relax_nreset": self.box_relax_nreset,
            "box_relax_min_style": self.box_relax_min_style,
            "box_relax_energy_tolerance": self.box_relax_energy_tolerance,
            "box_relax_force_tolerance": self.box_relax_force_tolerance,
            "box_relax_max_iterations": self.box_relax_max_iterations,
            "box_relax_max_evaluations": self.box_relax_max_evaluations,
        }


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
    n_atom_types: int
    n_bond_types: int
    n_angle_types: int
    atom_type_symbols: dict[int, str]
    settings: LammpsOptimizationSettings
    forcefield_backend: str
    parameter_sources: dict[str, str]
    warnings: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, object]:
        return {
            "input_cif": self.input_cif,
            "optimized_cif": self.optimized_cif,
            "output_dir": self.output_dir,
            "lammps_binary": self.lammps_binary,
            "lammps_data_path": self.lammps_data_path,
            "lammps_input_script_path": self.lammps_input_script_path,
            "lammps_dump_path": self.lammps_dump_path,
            "lammps_log_path": self.lammps_log_path,
            "stdout_log_path": self.stdout_log_path,
            "stderr_log_path": self.stderr_log_path,
            "report_path": self.report_path,
            "n_atoms": self.n_atoms,
            "n_bonds": self.n_bonds,
            "n_angles": self.n_angles,
            "n_atom_types": self.n_atom_types,
            "n_bond_types": self.n_bond_types,
            "n_angle_types": self.n_angle_types,
            "atom_type_symbols": dict(self.atom_type_symbols),
            "settings": self.settings.to_dict(),
            "forcefield_backend": self.forcefield_backend,
            "parameter_sources": dict(self.parameter_sources),
            "warnings": list(self.warnings),
        }


@dataclass(frozen=True)
class _CifAtomRecord:
    atom_id: int
    label: str
    symbol: str
    occupancy: float
    fractional_position: tuple[float, float, float]


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


@dataclass(frozen=True)
class _CifAngleRecord:
    angle_id: int
    atom_id_1: int
    atom_id_2: int
    atom_id_3: int
    equilibrium_degrees: float


@dataclass(frozen=True)
class _ParsedExplicitBondCif:
    source_path: Path
    data_name: str
    cell_parameters: tuple[float, float, float, float, float, float]
    lammps_basis: tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]]
    atoms: tuple[_CifAtomRecord, ...]
    bonds: tuple[_CifBondRecord, ...]
    angles: tuple[_CifAngleRecord, ...]


@dataclass(frozen=True)
class _PreparedLammpsSystem:
    data_text: str
    atom_type_symbols: dict[int, str]
    n_bond_types: int
    n_angle_types: int
    forcefield_backend: str
    parameter_sources: dict[str, str]
    warnings: tuple[str, ...]


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


def resolve_lammps_binary(lmp_path: str | Path | None = None) -> Path:
    raw_value = str(lmp_path) if lmp_path is not None else os.environ.get(COFKIT_LMP_ENV_VAR)
    if raw_value:
        path = Path(raw_value).expanduser()
        if not path.is_file():
            resolved = shutil.which(raw_value)
            if resolved:
                path = Path(resolved)
    elif DEFAULT_LAMMPS_BINARY.is_file():
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
    settings: LammpsOptimizationSettings | None = None,
    timeout_seconds: float = 300.0,
) -> LammpsOptimizationResult:
    if timeout_seconds <= 0.0:
        raise ValueError("timeout_seconds must be positive.")
    settings = settings or LammpsOptimizationSettings()
    _validate_settings(settings)

    input_path = Path(cif_path).expanduser().resolve()
    if not input_path.is_file():
        raise FileNotFoundError(f"CIF file does not exist: {input_path}")

    binary = resolve_lammps_binary(lmp_path)
    parsed = _parse_explicit_bond_cif(input_path)
    prepared = _prepare_lammps_system(parsed, settings=settings)
    run_dir = _resolve_output_dir(input_path, output_dir)
    run_dir.mkdir(parents=True, exist_ok=True)

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
            settings=settings,
            parsed=parsed,
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
    final_cartesian_positions = _parse_lammps_dump_last_frame(dump_path, expected_atoms=len(parsed.atoms))
    final_fractional_positions = _cartesian_positions_to_fractional(parsed, final_cartesian_positions)
    optimized_cif_path.write_text(
        _render_optimized_cif(parsed, final_fractional_positions),
        encoding="utf-8",
    )

    all_warnings = tuple(prepared.warnings) + tuple(execution_warnings)
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
        n_atoms=len(parsed.atoms),
        n_bonds=len(parsed.bonds),
        n_angles=len(parsed.angles),
        n_atom_types=len(prepared.atom_type_symbols),
        n_bond_types=prepared.n_bond_types,
        n_angle_types=prepared.n_angle_types,
        atom_type_symbols=prepared.atom_type_symbols,
        settings=settings,
        forcefield_backend=prepared.forcefield_backend,
        parameter_sources=prepared.parameter_sources,
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
    if settings.pair_cutoff <= 0.0:
        raise ValueError("pair_cutoff must be positive.")
    if settings.position_restraint_force_constant < 0.0:
        raise ValueError("position_restraint_force_constant must be non-negative.")
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


def _normalize_forcefield_name(forcefield: str) -> str:
    return forcefield.strip().lower().replace("-", "_")


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

    atom_images, image_conflicts = _compute_unwrapped_atom_images(parsed)
    typing_cartesian_positions = _unwrapped_cartesian_positions(parsed, atom_images)
    wrapped_cartesian_positions = {
        atom.atom_id: _fractional_to_lammps(atom.fractional_position, parsed.lammps_basis) for atom in parsed.atoms
    }
    ob_molecule = _build_openbabel_molecule(parsed, typing_cartesian_positions)
    atom_type_by_atom_id = _assign_openbabel_uff_atom_types(ob_molecule)
    rdkit_molecule = _build_rdkit_molecule_from_openbabel(ob_molecule)
    if not rdForceFieldHelpers.UFFHasAllMoleculeParams(rdkit_molecule):
        raise LammpsInputError(
            f"RDKit UFF parameterization is incomplete for {parsed.source_path}. "
            "The CIF-derived explicit-bond graph could not be fully typed by UFF."
        )

    ordered_type_labels, representative_atom_ids, representative_symbols = _ordered_representative_atom_types(
        parsed.atoms,
        atom_type_by_atom_id,
    )
    pairij_coeffs = _build_uff_pairij_coeffs(rdkit_molecule, ordered_type_labels, representative_atom_ids)
    bond_coeffs, angle_coeffs = _build_uff_topology_coeffs(parsed, rdkit_molecule, atom_type_by_atom_id)

    lattice = Lattice.from_parameters(*parsed.cell_parameters)
    structure = Structure(
        lattice,
        [atom.symbol for atom in parsed.atoms],
        [wrapped_cartesian_positions[atom.atom_id] for atom in parsed.atoms],
        coords_are_cartesian=True,
        to_unit_cell=False,
        site_properties={"ff_map": [atom_type_by_atom_id[atom.atom_id] for atom in parsed.atoms]},
        labels=[atom.label for atom in parsed.atoms],
    )
    topology_sections: dict[str, list[list[int]]] = {
        "Bonds": [[bond.atom_id_1 - 1, bond.atom_id_2 - 1] for bond in parsed.bonds],
    }
    if parsed.angles:
        topology_sections["Angles"] = [
            [angle.atom_id_1 - 1, angle.atom_id_2 - 1, angle.atom_id_3 - 1] for angle in parsed.angles
        ]

    a_vec, b_vec, c_vec = parsed.lammps_basis
    box = LammpsBox(
        bounds=[
            [0.0, a_vec[0]],
            [0.0, b_vec[1]],
            [0.0, c_vec[2]],
        ],
        tilt=[b_vec[0], c_vec[0], c_vec[1]],
    )
    topo_coeffs: dict[str, list[dict[str, object]]] = {}
    if bond_coeffs:
        topo_coeffs["Bond Coeffs"] = [
            {"coeffs": [coeffs[0], coeffs[1]], "types": [bond_type]} for bond_type, coeffs in bond_coeffs
        ]
    if angle_coeffs:
        topo_coeffs["Angle Coeffs"] = [
            {"coeffs": [coeffs[0], coeffs[1]], "types": [angle_type]} for angle_type, coeffs in angle_coeffs
        ]
    forcefield = ForceField(
        mass_info=[(label, representative_symbols[label]) for label in ordered_type_labels],
        nonbond_coeffs=pairij_coeffs,
        topo_coeffs=topo_coeffs or None,
    )
    topology = Topology(structure, ff_label="ff_map", topologies=topology_sections)
    lammps_data = LammpsData.from_ff_and_topologies(box, forcefield, [topology], atom_style="molecular")
    data_text = lammps_data.get_str(distance=10, charge=4, hybrid=False)
    data_text = data_text.replace(
        "Generated by pymatgen.io.lammps.data.LammpsData",
        "LAMMPS data file written by cofkit via pymatgen/Open Babel/RDKit",
        1,
    )
    data_text = data_text.replace("\nAtoms\n\n", "\nAtoms # molecular\n\n", 1)

    warnings: list[str] = [
        "UFF-backed LAMMPS export currently includes bonded stretch/bend and van der Waals terms. "
        "Dihedral and improper terms are not emitted yet.",
    ]
    if image_conflicts > 0:
        warnings.append(
            "Periodic spanning-tree unwrapping for UFF atom typing encountered one or more cell-cycle conflicts. "
            "cofkit kept the wrapped CIF coordinates for the LAMMPS data file and used the spanning tree only for "
            "local bond-order perception."
        )
    parameter_sources = {
        "forcefield": "UFF",
        "atom_typing": "Open Babel UFF",
        "bonded_parameters": "RDKit UFF parameter helpers",
        "nonbond_parameters": "RDKit UFF parameter helpers",
    }
    uff_parameter_file = _discover_uff_parameter_file()
    if uff_parameter_file is not None:
        parameter_sources["reference_parameter_file"] = str(uff_parameter_file)

    return _PreparedLammpsSystem(
        data_text=data_text,
        atom_type_symbols={index: label for index, label in enumerate(ordered_type_labels, start=1)},
        n_bond_types=len(bond_coeffs),
        n_angle_types=len(angle_coeffs),
        forcefield_backend="uff_openbabel_rdkit_pymatgen",
        parameter_sources=parameter_sources,
        warnings=tuple(warnings),
    )


def _discover_uff_parameter_file() -> Path | None:
    candidate_roots = [
        Path(sys.prefix) / "share",
        Path(sys.base_prefix) / "share",
        Path("/opt/homebrew/share"),
        Path("/usr/local/share"),
    ]
    seen: set[Path] = set()
    for root in candidate_roots:
        if not root.exists():
            continue
        for path in root.glob("openbabel*/**/UFF.prm"):
            resolved = path.resolve()
            if resolved in seen:
                continue
            seen.add(resolved)
            return resolved
    return None


def _require_uff_support() -> None:
    missing: list[str] = []
    if gemmi is None:
        missing.append("gemmi")
    if ob is None:
        missing.append("openbabel")
    if Chem is None or rdForceFieldHelpers is None:
        missing.append("rdkit")
    if any(item is None for item in (ForceField, LammpsBox, LammpsData, Lattice, Structure, Topology)):
        missing.append("pymatgen")
    if missing:
        raise LammpsConfigurationError(
            "UFF-backed LAMMPS optimization requires these Python packages in the current environment: "
            + ", ".join(missing)
        )


def _ordered_representative_atom_types(
    atoms: tuple[_CifAtomRecord, ...],
    atom_type_by_atom_id: dict[int, str],
) -> tuple[list[str], dict[str, int], dict[str, str]]:
    ordered_types: list[str] = []
    representative_atom_ids: dict[str, int] = {}
    representative_symbols: dict[str, str] = {}
    for atom in atoms:
        atom_type = atom_type_by_atom_id[atom.atom_id]
        if atom_type in representative_atom_ids:
            continue
        ordered_types.append(atom_type)
        representative_atom_ids[atom_type] = atom.atom_id
        representative_symbols[atom_type] = atom.symbol
    return ordered_types, representative_atom_ids, representative_symbols


def _build_uff_pairij_coeffs(
    rdkit_molecule,
    ordered_type_labels: list[str],
    representative_atom_ids: dict[str, int],
) -> list[list[float]]:
    coeffs: list[list[float]] = []
    for left_index, left_label in enumerate(ordered_type_labels):
        for right_label in ordered_type_labels[left_index:]:
            params = rdForceFieldHelpers.GetUFFVdWParams(
                rdkit_molecule,
                representative_atom_ids[left_label] - 1,
                representative_atom_ids[right_label] - 1,
            )
            if params is None:
                raise LammpsInputError(
                    f"UFF van der Waals parameters are unavailable for atom types {left_label!r} and {right_label!r}."
                )
            x_ij, d_ij = params
            coeffs.append([float(d_ij), float(x_ij) / _UFF_RMIN_TO_SIGMA_FACTOR])
    return coeffs


def _build_uff_topology_coeffs(
    parsed: _ParsedExplicitBondCif,
    rdkit_molecule,
    atom_type_by_atom_id: dict[int, str],
) -> tuple[list[tuple[tuple[str, str], tuple[float, float]]], list[tuple[tuple[str, str, str], tuple[float, float]]]]:
    bond_coeffs_by_type: dict[tuple[str, str], tuple[float, float]] = {}
    angle_coeffs_by_type: dict[tuple[str, str, str], tuple[float, float]] = {}

    for bond in parsed.bonds:
        left_type = atom_type_by_atom_id[bond.atom_id_1]
        right_type = atom_type_by_atom_id[bond.atom_id_2]
        bond_type = _canonical_bond_type(left_type, right_type)
        params = rdForceFieldHelpers.GetUFFBondStretchParams(
            rdkit_molecule,
            bond.atom_id_1 - 1,
            bond.atom_id_2 - 1,
        )
        if params is None:
            raise LammpsInputError(
                f"UFF bond parameters are unavailable for bond {bond.label_1}-{bond.label_2} "
                f"with atom types {left_type!r}/{right_type!r}."
            )
        kb, r0 = params
        _store_unique_coefficients(
            bond_coeffs_by_type,
            bond_type,
            (0.5 * float(kb), float(r0)),
            kind="bond",
        )

    for angle in parsed.angles:
        left_type = atom_type_by_atom_id[angle.atom_id_1]
        center_type = atom_type_by_atom_id[angle.atom_id_2]
        right_type = atom_type_by_atom_id[angle.atom_id_3]
        angle_type = _canonical_angle_type(left_type, center_type, right_type)
        params = rdForceFieldHelpers.GetUFFAngleBendParams(
            rdkit_molecule,
            angle.atom_id_1 - 1,
            angle.atom_id_2 - 1,
            angle.atom_id_3 - 1,
        )
        if params is None:
            raise LammpsInputError(
                "UFF angle parameters are unavailable for angle "
                f"{angle.atom_id_1}-{angle.atom_id_2}-{angle.atom_id_3} with atom types "
                f"{left_type!r}/{center_type!r}/{right_type!r}."
            )
        ka, theta0 = params
        _store_unique_coefficients(
            angle_coeffs_by_type,
            angle_type,
            (0.5 * float(ka), float(theta0)),
            kind="angle",
        )

    return list(bond_coeffs_by_type.items()), list(angle_coeffs_by_type.items())


def _store_unique_coefficients(
    table: dict[tuple[str, ...], tuple[float, float]],
    key: tuple[str, ...],
    values: tuple[float, float],
    *,
    kind: str,
) -> None:
    existing = table.get(key)
    if existing is None:
        table[key] = values
        return
    if _coefficients_close(existing, values):
        return
    raise LammpsInputError(
        f"Encountered conflicting {kind} coefficients for forcefield type {key!r}: "
        f"{existing!r} vs {values!r}."
    )


def _coefficients_close(left: tuple[float, float], right: tuple[float, float]) -> bool:
    for a, b in zip(left, right):
        absolute_delta = abs(a - b)
        relative_scale = max(abs(a), abs(b), 1.0)
        if absolute_delta > 5.0e-2 and absolute_delta / relative_scale > 2.0e-2:
            return False
    return True


def _canonical_bond_type(left_type: str, right_type: str) -> tuple[str, str]:
    return min((left_type, right_type), (right_type, left_type))


def _canonical_angle_type(left_type: str, center_type: str, right_type: str) -> tuple[str, str, str]:
    forward = (left_type, center_type, right_type)
    reverse = (right_type, center_type, left_type)
    return min(forward, reverse)


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
    for bond in parsed.bonds:
        if not molecule.AddBond(bond.atom_id_1, bond.atom_id_2, 1):
            raise LammpsInputError(
                f"Open Babel rejected explicit bond {bond.label_1}-{bond.label_2} while preparing the UFF model."
            )
    molecule.PerceiveBondOrders()
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


def _build_rdkit_molecule_from_openbabel(ob_molecule):
    converter = ob.OBConversion()
    if not converter.SetOutFormat("mol"):
        raise LammpsConfigurationError("Open Babel does not support MOL export in this environment.")
    mol_block = converter.WriteString(ob_molecule)
    rdkit_molecule = Chem.MolFromMolBlock(mol_block, sanitize=True, removeHs=False)
    if rdkit_molecule is None:
        raise LammpsInputError(
            "RDKit could not reconstruct a sanitized molecule from the Open Babel bond-order perception step."
        )
    return rdkit_molecule


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
            reported_distance=equilibrium_distance,
        )
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
            )
        )

    angles = _derive_angles(tuple(atoms), tuple(bonds), basis)
    return _ParsedExplicitBondCif(
        source_path=cif_path,
        data_name=str(block.name or cif_path.stem),
        cell_parameters=(a, b, c, alpha, beta, gamma),
        lammps_basis=basis,
        atoms=tuple(atoms),
        bonds=tuple(bonds),
        angles=angles,
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


def _resolve_effective_bond_geometry(
    *,
    left_fractional: tuple[float, float, float],
    right_fractional: tuple[float, float, float],
    raw_shift_1: tuple[int, int, int],
    raw_shift_2: tuple[int, int, int],
    raw_symmetry_1: str,
    raw_symmetry_2: str,
    basis: tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]],
    reported_distance: float,
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
    if settings.position_restraint_force_constant <= 0.0:
        return 0.0
    return min(settings.position_restraint_force_constant * 0.25, 0.05)


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


def _render_lammps_input_script(
    *,
    data_path: Path,
    dump_path: Path,
    settings: LammpsOptimizationSettings,
    parsed: _ParsedExplicitBondCif,
) -> str:
    forcefield = _normalize_forcefield_name(settings.forcefield)
    thermo_terms = ["step", "pe", "ebond"]
    minimization_stages = _build_minimization_stages(settings)
    lines = [
        "units real",
        "atom_style molecular",
        "boundary p p p",
        f"pair_style lj/cut {settings.pair_cutoff:.6f}",
        "pair_modify shift yes",
        "special_bonds lj 0.0 0.0 1.0",
        "bond_style harmonic",
    ]
    if parsed.angles:
        lines.append("angle_style harmonic")
        thermo_terms.append("eangle")
    thermo_terms.extend(["evdwl", "press"])
    lines.extend(
        [
            f"read_data {data_path}",
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
    if forcefield == "uff":
        return "\n".join(lines)
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


def _parse_lammps_dump_last_frame(
    dump_path: Path,
    *,
    expected_atoms: int,
) -> dict[int, tuple[float, float, float]]:
    if not dump_path.is_file():
        raise LammpsParseError(f"LAMMPS dump file was not created: {dump_path}")

    lines = dump_path.read_text(encoding="utf-8").splitlines()
    index = 0
    last_frame: dict[int, tuple[float, float, float]] | None = None
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

    if last_frame is None:
        raise LammpsParseError(f"LAMMPS dump file did not contain any frames: {dump_path}")
    if len(last_frame) != expected_atoms:
        raise LammpsParseError(
            f"LAMMPS dump file {dump_path} contained {len(last_frame)} atoms in the final frame, expected {expected_atoms}"
        )
    return last_frame


def _cartesian_positions_to_fractional(
    parsed: _ParsedExplicitBondCif,
    final_cartesian_positions: dict[int, tuple[float, float, float]],
) -> dict[int, tuple[float, float, float]]:
    a_vec, b_vec, c_vec = parsed.lammps_basis
    ax = a_vec[0]
    by = b_vec[1]
    xy = b_vec[0]
    xz = c_vec[0]
    yz = c_vec[1]
    cz = c_vec[2]

    positions: dict[int, tuple[float, float, float]] = {}
    for atom_id, (x, y, z) in final_cartesian_positions.items():
        if abs(cz) < 1e-12 or abs(by) < 1e-12 or abs(ax) < 1e-12:
            raise LammpsParseError("LAMMPS cell basis is singular and cannot be converted back to fractional coordinates.")
        fz = z / cz
        fy = (y - fz * yz) / by
        fx = (x - fy * xy - fz * xz) / ax
        positions[atom_id] = (_wrap_fraction(fx), _wrap_fraction(fy), _wrap_fraction(fz))
    return positions


def _render_optimized_cif(
    parsed: _ParsedExplicitBondCif,
    final_fractional_positions: dict[int, tuple[float, float, float]],
) -> str:
    a, b, c, alpha, beta, gamma = parsed.cell_parameters
    lines = [
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
    for atom in parsed.atoms:
        fractional = final_fractional_positions[atom.atom_id]
        lines.append(
            f"{atom.label} {atom.symbol} {fractional[0]:.6f} {fractional[1]:.6f} {fractional[2]:.6f} {atom.occupancy:.2f}"
        )

    lines.extend(
        [
            "",
            "loop_",
            "_geom_bond_atom_site_label_1",
            "_geom_bond_atom_site_label_2",
            "_geom_bond_site_symmetry_1",
            "_geom_bond_site_symmetry_2",
            "_geom_bond_distance",
        ]
    )
    for bond in parsed.bonds:
        fractional_1 = final_fractional_positions[bond.atom_id_1]
        fractional_2 = final_fractional_positions[bond.atom_id_2]
        if bond.symmetry_1 == "." and bond.symmetry_2 == ".":
            output_shift_1, output_shift_2 = _closest_periodic_shifts(
                left_fractional=fractional_1,
                right_fractional=fractional_2,
                basis=parsed.lammps_basis,
            )
            symmetry_1 = _format_p1_symmetry_shift(output_shift_1)
            symmetry_2 = _format_p1_symmetry_shift(output_shift_2)
        else:
            output_shift_1 = bond.shift_1
            output_shift_2 = bond.shift_2
            symmetry_1 = bond.symmetry_1 if bond.symmetry_1 != "." else _format_p1_symmetry_shift(output_shift_1)
            symmetry_2 = bond.symmetry_2 if bond.symmetry_2 != "." else _format_p1_symmetry_shift(output_shift_2)
        distance = _bond_distance(
            fractional_1,
            fractional_2,
            output_shift_1,
            output_shift_2,
            parsed.lammps_basis,
        )
        lines.append(
            f"{bond.label_1} {bond.label_2} {symmetry_1} {symmetry_2} {distance:.6f}"
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


def _atomic_mass(symbol: str) -> float:
    if gemmi is None:  # pragma: no cover - guarded earlier
        return 12.011
    mass = float(gemmi.Element(symbol).weight)
    return mass if mass > 0.0 else 12.011


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
