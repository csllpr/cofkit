from __future__ import annotations

import re
from dataclasses import dataclass, field
from math import acos, degrees, sqrt
from pathlib import Path
from typing import Mapping

try:
    import gemmi
except ImportError:  # pragma: no cover - exercised in incomplete environments
    gemmi = None


Vec3 = tuple[float, float, float]
CellBasis = tuple[Vec3, Vec3, Vec3]


_ELEMENT_RE = re.compile(r"^[^A-Za-z]*([A-Za-z]{1,2})")
_UNCERTAINTY_RE = re.compile(r"^([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?)(?:\(\d+\))?$")


@dataclass(frozen=True)
class PeriodicCifAtoms:
    """Minimal no-ASE periodic atom container for CIF decomposition workflows."""

    symbols: tuple[str, ...]
    fractional_positions: tuple[Vec3, ...]
    cartesian_positions: tuple[Vec3, ...]
    cell_basis: CellBasis
    info: Mapping[str, tuple[str, ...]] = field(default_factory=dict)
    data_name: str = ""
    source_path: str = ""

    def __post_init__(self) -> None:
        if len(self.symbols) != len(self.fractional_positions):
            raise ValueError("symbols and fractional_positions must have the same length")
        if len(self.symbols) != len(self.cartesian_positions):
            raise ValueError("symbols and cartesian_positions must have the same length")

    def __len__(self) -> int:
        return len(self.symbols)

    def get_chemical_symbols(self) -> list[str]:
        return list(self.symbols)

    def get_scaled_positions(self, wrap: bool = True) -> list[Vec3]:
        if not wrap:
            return list(self.fractional_positions)
        return [tuple(coord % 1.0 for coord in position) for position in self.fractional_positions]  # type: ignore[list-item]

    def get_positions(self) -> list[Vec3]:
        return list(self.cartesian_positions)

    def get_atomic_numbers(self) -> list[int]:
        return [_atomic_number(symbol) for symbol in self.symbols]

    def repeat(self, repeat: tuple[int, int, int]) -> "PeriodicCifAtoms":
        repeat = _validate_repeat(repeat)
        repeated_symbols: list[str] = []
        repeated_fractional: list[Vec3] = []
        repeated_cartesian: list[Vec3] = []
        nx, ny, nz = repeat
        repeated_basis = (
            _scale_vec(self.cell_basis[0], nx),
            _scale_vec(self.cell_basis[1], ny),
            _scale_vec(self.cell_basis[2], nz),
        )

        for ix in range(nx):
            for iy in range(ny):
                for iz in range(nz):
                    for symbol, fractional in zip(self.symbols, self.fractional_positions):
                        new_fractional = (
                            (fractional[0] + ix) / nx,
                            (fractional[1] + iy) / ny,
                            (fractional[2] + iz) / nz,
                        )
                        repeated_symbols.append(symbol)
                        repeated_fractional.append(new_fractional)
                        repeated_cartesian.append(_fractional_to_cartesian(repeated_basis, new_fractional))

        return PeriodicCifAtoms(
            symbols=tuple(repeated_symbols),
            fractional_positions=tuple(repeated_fractional),
            cartesian_positions=tuple(repeated_cartesian),
            cell_basis=repeated_basis,
            info={},
            data_name=self.data_name,
            source_path=self.source_path,
        )

    @property
    def cell_parameters(self) -> tuple[float, float, float, float, float, float]:
        return _cell_parameters_from_basis(self.cell_basis)


def read_periodic_cif_atoms(cif_path: str | Path) -> PeriodicCifAtoms:
    if gemmi is None:
        raise ModuleNotFoundError("gemmi is required for CIF extraction.")

    input_path = Path(cif_path)
    document = gemmi.cif.read_file(str(input_path))
    for block in document:
        if _block_has_required_atom_sites(block):
            return _atoms_from_gemmi_block(block, source_path=input_path)
    raise ValueError(f"No structural CIF block with fractional atom sites found in {input_path}")


def atoms_have_explicit_cif_bonds(atoms: PeriodicCifAtoms) -> bool:
    return bool(atoms.info.get("_geom_bond_atom_site_label_1")) and bool(
        atoms.info.get("_geom_bond_atom_site_label_2")
    )


def get_primitive_periodic_cif_atoms(
    cif_path: str | Path,
    *,
    output_path: str | Path | None = None,
    symprec: float = 0.1,
) -> PeriodicCifAtoms:
    try:
        from pymatgen.core import Structure
        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    except ImportError as exc:  # pragma: no cover - project dependency guard
        raise ModuleNotFoundError("pymatgen is required for primitive CIF extraction.") from exc

    structure = Structure.from_file(str(cif_path))
    primitive = SpacegroupAnalyzer(structure, symprec=symprec).get_primitive_standard_structure()
    if output_path is not None:
        out_path = Path(output_path)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        primitive.to(filename=str(out_path))
    return periodic_cif_atoms_from_pymatgen_structure(
        primitive,
        data_name=Path(cif_path).stem,
        source_path=str(cif_path),
    )


def prepare_periodic_cif_atoms(
    cif_path: str | Path,
    primitive_dir: str | Path | None = None,
    symprec: float = 0.1,
) -> tuple[PeriodicCifAtoms, bool]:
    input_path = Path(cif_path)
    try:
        source_atoms = read_periodic_cif_atoms(input_path)
    except Exception:
        source_atoms = None

    if source_atoms is not None and atoms_have_explicit_cif_bonds(source_atoms):
        return source_atoms, True

    output_path = None
    if primitive_dir is not None:
        output_path = Path(primitive_dir) / f"{input_path.stem}_primitive.cif"
    return get_primitive_periodic_cif_atoms(input_path, output_path=output_path, symprec=symprec), False


def periodic_cif_atoms_from_pymatgen_structure(
    structure,
    *,
    data_name: str = "",
    source_path: str = "",
) -> PeriodicCifAtoms:
    symbols: list[str] = []
    fractional_positions: list[Vec3] = []
    for site in structure:
        symbols.append(_symbol_from_pymatgen_site(site))
        fractional_positions.append(tuple(float(value) for value in site.frac_coords))  # type: ignore[arg-type]

    basis = tuple(
        tuple(float(value) for value in vector)
        for vector in structure.lattice.matrix
    )
    labels = tuple(f"{symbol}{index}" for index, symbol in enumerate(symbols, start=1))
    return PeriodicCifAtoms(
        symbols=tuple(symbols),
        fractional_positions=tuple(fractional_positions),
        cartesian_positions=tuple(_fractional_to_cartesian(basis, position) for position in fractional_positions),
        cell_basis=basis,  # type: ignore[arg-type]
        info={"_atom_site_label": labels},
        data_name=data_name,
        source_path=source_path,
    )


def _atoms_from_gemmi_block(block, *, source_path: Path) -> PeriodicCifAtoms:
    labels = _required_loop(block, "_atom_site_label")
    type_symbols = _optional_loop(block, "_atom_site_type_symbol")
    fract_x = _required_loop(block, "_atom_site_fract_x")
    fract_y = _required_loop(block, "_atom_site_fract_y")
    fract_z = _required_loop(block, "_atom_site_fract_z")
    row_count = len(labels)
    if not (len(fract_x) == len(fract_y) == len(fract_z) == row_count):
        raise ValueError("CIF atom-site label and fractional-coordinate columns have different lengths.")
    if type_symbols and len(type_symbols) != row_count:
        raise ValueError("CIF atom-site type-symbol column length does not match atom labels.")

    cell = _gemmi_cell_from_block(block)
    basis = _basis_from_gemmi_cell(cell)
    symbols: list[str] = []
    fractional_positions: list[Vec3] = []
    for index in range(row_count):
        raw_symbol = type_symbols[index] if type_symbols else labels[index]
        symbols.append(_normalize_element_symbol(raw_symbol))
        fractional_positions.append(
            (
                _parse_cif_float(fract_x[index]),
                _parse_cif_float(fract_y[index]),
                _parse_cif_float(fract_z[index]),
            )
        )

    info = _collect_cif_info(block)
    if "_atom_site_label" not in info:
        info["_atom_site_label"] = tuple(_normalize_cif_value(label) for label in labels)

    return PeriodicCifAtoms(
        symbols=tuple(symbols),
        fractional_positions=tuple(fractional_positions),
        cartesian_positions=tuple(_fractional_to_cartesian(basis, position) for position in fractional_positions),
        cell_basis=basis,
        info=info,
        data_name=str(block.name or source_path.stem),
        source_path=str(source_path),
    )


def _collect_cif_info(block) -> dict[str, tuple[str, ...]]:
    tags = (
        "_atom_site_label",
        "_space_group_symop_operation_xyz",
        "_symmetry_equiv_pos_as_xyz",
        "_geom_bond_atom_site_label_1",
        "_geom_bond_atom_site_label_2",
        "_geom_bond_site_symmetry_1",
        "_geom_bond_site_symmetry_2",
        "_geom_bond_distance",
        "_ccdc_geom_bond_type",
        "_geom_bond_type",
    )
    info: dict[str, tuple[str, ...]] = {}
    for tag in tags:
        values = _optional_loop(block, tag)
        if values:
            info[tag] = tuple(_normalize_cif_value(value) for value in values)
    return info


def _block_has_required_atom_sites(block) -> bool:
    return (
        bool(_optional_loop(block, "_atom_site_label"))
        and bool(_optional_loop(block, "_atom_site_fract_x"))
        and bool(_optional_loop(block, "_atom_site_fract_y"))
        and bool(_optional_loop(block, "_atom_site_fract_z"))
        and bool(block.find_value("_cell_length_a"))
    )


def _required_loop(block, tag: str) -> tuple[str, ...]:
    values = _optional_loop(block, tag)
    if not values:
        raise ValueError(f"CIF block is missing required loop value {tag}.")
    return values


def _optional_loop(block, tag: str) -> tuple[str, ...]:
    values = block.find_loop(tag)
    if len(values) == 0:
        return ()
    return tuple(str(value) for value in values)


def _gemmi_cell_from_block(block):
    assert gemmi is not None
    return gemmi.UnitCell(
        _parse_cif_float(block.find_value("_cell_length_a")),
        _parse_cif_float(block.find_value("_cell_length_b")),
        _parse_cif_float(block.find_value("_cell_length_c")),
        _parse_cif_float(block.find_value("_cell_angle_alpha")),
        _parse_cif_float(block.find_value("_cell_angle_beta")),
        _parse_cif_float(block.find_value("_cell_angle_gamma")),
    )


def _basis_from_gemmi_cell(cell) -> CellBasis:
    origin = cell.orthogonalize(gemmi.Fractional(0.0, 0.0, 0.0))
    basis = []
    for fractional in (
        gemmi.Fractional(1.0, 0.0, 0.0),
        gemmi.Fractional(0.0, 1.0, 0.0),
        gemmi.Fractional(0.0, 0.0, 1.0),
    ):
        point = cell.orthogonalize(fractional)
        basis.append((float(point.x - origin.x), float(point.y - origin.y), float(point.z - origin.z)))
    return tuple(basis)  # type: ignore[return-value]


def _fractional_to_cartesian(basis: CellBasis, fractional: Vec3) -> Vec3:
    return (
        fractional[0] * basis[0][0] + fractional[1] * basis[1][0] + fractional[2] * basis[2][0],
        fractional[0] * basis[0][1] + fractional[1] * basis[1][1] + fractional[2] * basis[2][1],
        fractional[0] * basis[0][2] + fractional[1] * basis[1][2] + fractional[2] * basis[2][2],
    )


def _cell_parameters_from_basis(basis: CellBasis) -> tuple[float, float, float, float, float, float]:
    a_vec, b_vec, c_vec = basis
    a = _norm(a_vec)
    b = _norm(b_vec)
    c = _norm(c_vec)
    return (
        a,
        b,
        c,
        _angle_degrees(b_vec, c_vec),
        _angle_degrees(a_vec, c_vec),
        _angle_degrees(a_vec, b_vec),
    )


def _angle_degrees(first: Vec3, second: Vec3) -> float:
    denominator = _norm(first) * _norm(second)
    if denominator <= 0.0:
        raise ValueError("cell basis vector length must be positive")
    cosine = max(-1.0, min(1.0, _dot(first, second) / denominator))
    return degrees(acos(cosine))


def _dot(first: Vec3, second: Vec3) -> float:
    return first[0] * second[0] + first[1] * second[1] + first[2] * second[2]


def _norm(vector: Vec3) -> float:
    return sqrt(_dot(vector, vector))


def _scale_vec(vector: Vec3, factor: int) -> Vec3:
    return (vector[0] * factor, vector[1] * factor, vector[2] * factor)


def _validate_repeat(repeat: tuple[int, int, int]) -> tuple[int, int, int]:
    if len(repeat) != 3:
        raise ValueError("repeat must contain exactly three integers")
    normalized = tuple(int(value) for value in repeat)
    if any(value <= 0 for value in normalized):
        raise ValueError("repeat values must be positive integers")
    return normalized  # type: ignore[return-value]


def _parse_cif_float(value: object) -> float:
    text = _normalize_cif_value(value)
    if text in {"", ".", "?"}:
        raise ValueError(f"CIF numeric value is missing: {value!r}")
    if "/" in text and not any(marker in text for marker in ("(", ")", "e", "E")):
        numerator, denominator = text.split("/", 1)
        return float(numerator) / float(denominator)
    match = _UNCERTAINTY_RE.match(text)
    if match:
        text = match.group(1)
    return float(text)


def _normalize_cif_value(value: object) -> str:
    return str(value).strip().strip("'\"")


def _normalize_element_symbol(value: object) -> str:
    text = _normalize_cif_value(value)
    if not text:
        raise ValueError("CIF atom site is missing an element symbol.")
    candidates = [text]
    match = _ELEMENT_RE.match(text)
    if match:
        raw = match.group(1)
        candidates.append(raw[:1].upper() + raw[1:].lower())
        candidates.append(raw[:1].upper())
    for candidate in candidates:
        try:
            element = gemmi.Element(candidate)
        except Exception:
            continue
        if int(element.atomic_number) > 0:
            return str(element.name)
    raise ValueError(f"Could not infer element symbol from CIF value {value!r}.")


def _atomic_number(symbol: str) -> int:
    if gemmi is None:
        raise ModuleNotFoundError("gemmi is required for atomic-number lookup.")
    element = gemmi.Element(symbol)
    atomic_number = int(element.atomic_number)
    if atomic_number <= 0:
        raise ValueError(f"unknown element symbol {symbol!r}")
    return atomic_number


def _symbol_from_pymatgen_site(site) -> str:
    if getattr(site, "is_ordered", True):
        return _normalize_element_symbol(site.specie.symbol)
    elements = list(site.species.elements)
    if not elements:
        raise ValueError("pymatgen site does not contain any element species")
    return _normalize_element_symbol(elements[0].symbol)


__all__ = [
    "PeriodicCifAtoms",
    "atoms_have_explicit_cif_bonds",
    "get_primitive_periodic_cif_atoms",
    "periodic_cif_atoms_from_pymatgen_structure",
    "prepare_periodic_cif_atoms",
    "read_periodic_cif_atoms",
]
