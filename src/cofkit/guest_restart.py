from __future__ import annotations

import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Sequence

from .guest_bundles import GuestBundle, GuestBundleError, load_guest_bundles


KCAL_PER_MOL_PER_K = 0.00198720425864083


class GuestRestartError(RuntimeError):
    """Raised when a GCMC guest snapshot cannot be converted for LAMMPS restart."""


@dataclass(frozen=True)
class LammpsGuestSite:
    label: str
    element: str
    mass: float
    charge: float
    epsilon_k: float
    sigma: float

    @property
    def epsilon_kcal_per_mol(self) -> float:
        return self.epsilon_k * KCAL_PER_MOL_PER_K

    def to_dict(self) -> dict[str, object]:
        return {
            "label": self.label,
            "element": self.element,
            "mass": self.mass,
            "charge": self.charge,
            "epsilon_k": self.epsilon_k,
            "epsilon_kcal_per_mol": self.epsilon_kcal_per_mol,
            "sigma": self.sigma,
        }


@dataclass(frozen=True)
class LammpsGuestTemplate:
    component: str
    site_labels: tuple[str, ...]
    relative_positions: tuple[tuple[float, float, float], ...]
    bonds: tuple[tuple[int, int], ...]
    bond_force_constant: float = 1000.0
    angle_force_constant: float = 100.0

    def to_dict(self) -> dict[str, object]:
        return {
            "component": self.component,
            "site_labels": list(self.site_labels),
            "relative_positions": [list(position) for position in self.relative_positions],
            "bonds": [list(bond) for bond in self.bonds],
            "bond_force_constant": self.bond_force_constant,
            "angle_force_constant": self.angle_force_constant,
        }


@dataclass(frozen=True)
class LammpsGuestSnapshotAtom:
    component: str
    molecule_key: str
    site_index: int
    site_label: str
    x: float
    y: float
    z: float

    def to_dict(self) -> dict[str, object]:
        return {
            "component": self.component,
            "molecule_key": self.molecule_key,
            "site_index": self.site_index,
            "site_label": self.site_label,
            "x": self.x,
            "y": self.y,
            "z": self.z,
        }


@dataclass(frozen=True)
class LammpsGuestRestartState:
    source_snapshot_path: str
    atoms: tuple[LammpsGuestSnapshotAtom, ...]
    sites: tuple[LammpsGuestSite, ...]
    templates: tuple[LammpsGuestTemplate, ...]
    warnings: tuple[str, ...] = ()

    @property
    def n_atoms(self) -> int:
        return len(self.atoms)

    @property
    def components(self) -> tuple[str, ...]:
        return tuple(dict.fromkeys(atom.component for atom in self.atoms))

    def site_by_label(self) -> dict[str, LammpsGuestSite]:
        return {site.label: site for site in self.sites}

    def template_by_component(self) -> dict[str, LammpsGuestTemplate]:
        return {template.component: template for template in self.templates}

    def to_dict(self) -> dict[str, object]:
        return {
            "source_snapshot_path": self.source_snapshot_path,
            "n_atoms": self.n_atoms,
            "components": list(self.components),
            "atoms": [atom.to_dict() for atom in self.atoms],
            "sites": [site.to_dict() for site in self.sites],
            "templates": [template.to_dict() for template in self.templates],
            "warnings": list(self.warnings),
        }


@dataclass(frozen=True)
class LammpsGuestRestartCell:
    unit_cells: tuple[int, int, int]
    basis: tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]]
    lengths: tuple[float, float, float]
    angles: tuple[float, float, float]

    def to_dict(self) -> dict[str, object]:
        return {
            "unit_cells": list(self.unit_cells),
            "basis": [list(vector) for vector in self.basis],
            "lengths": list(self.lengths),
            "angles": list(self.angles),
        }


@dataclass(frozen=True)
class GraspaRestartFileResult:
    restart_file_path: str
    n_adsorbate_atoms: int
    n_adsorbate_molecules: int
    components: tuple[str, ...]
    warnings: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, object]:
        return {
            "restart_file_path": self.restart_file_path,
            "n_adsorbate_atoms": self.n_adsorbate_atoms,
            "n_adsorbate_molecules": self.n_adsorbate_molecules,
            "components": list(self.components),
            "warnings": list(self.warnings),
        }


def build_lammps_guest_restart_state_from_gcmc_result(
    gcmc_result: object,
    *,
    components: Sequence[str],
    guest_bundles: Sequence[str] = (),
) -> LammpsGuestRestartState:
    snapshot_path = find_latest_gcmc_movie_snapshot(gcmc_result)
    if snapshot_path is None:
        raise GuestRestartError(
            "No gRASPA/RASPA2 movie snapshot was found under the GCMC run directory. "
            "Guest-containing hybrid exchange requires a LAMMPS-style Movies/System_0/result_*.data file."
        )
    templates, sites = load_lammps_guest_force_field_assets(components, guest_bundles=guest_bundles)
    return parse_lammps_guest_restart_snapshot(snapshot_path, templates=templates, sites=sites)


def build_lammps_guest_restart_state_from_lammps_md_result(
    lammps_md_result: object,
    *,
    previous_guest_restart_state: LammpsGuestRestartState,
) -> tuple[LammpsGuestRestartState, LammpsGuestRestartCell]:
    data_path = Path(getattr(lammps_md_result, "lammps_data_path")).expanduser().resolve()
    dump_path = Path(getattr(lammps_md_result, "lammps_dump_path")).expanduser().resolve()
    if not data_path.is_file():
        raise GuestRestartError(f"LAMMPS MD data file does not exist: {data_path}")
    positions_by_atom_id, cell = parse_lammps_md_dump_guest_positions(dump_path)
    identities = _parse_lammps_guest_atom_identities_from_data(
        data_path.read_text(encoding="utf-8", errors="replace"),
        previous_guest_restart_state=previous_guest_restart_state,
    )
    atoms: list[LammpsGuestSnapshotAtom] = []
    warnings: list[str] = []
    for identity in identities:
        position = positions_by_atom_id.get(identity["atom_id"])
        if position is None:
            warnings.append(
                f"LAMMPS dump did not contain final coordinates for guest atom id {identity['atom_id']}; skipped."
            )
            continue
        x, y, z = position
        atoms.append(
            LammpsGuestSnapshotAtom(
                component=identity["component"],
                molecule_key=identity["molecule_key"],
                site_index=identity["site_index"],
                site_label=identity["site_label"],
                x=x,
                y=y,
                z=z,
            )
        )
    if not atoms:
        raise GuestRestartError(
            f"No guest atoms from the previous restart state were recovered from LAMMPS MD dump {dump_path}."
        )
    state = LammpsGuestRestartState(
        source_snapshot_path=str(dump_path),
        atoms=tuple(atoms),
        sites=previous_guest_restart_state.sites,
        templates=previous_guest_restart_state.templates,
        warnings=tuple(dict.fromkeys((*previous_guest_restart_state.warnings, *warnings))),
    )
    return state, cell


def parse_lammps_md_dump_guest_positions(
    dump_path: str | Path,
) -> tuple[dict[int, tuple[float, float, float]], LammpsGuestRestartCell]:
    path = Path(dump_path).expanduser().resolve()
    if not path.is_file():
        raise GuestRestartError(f"LAMMPS MD dump file does not exist: {path}")
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    index = 0
    last_positions: dict[int, tuple[float, float, float]] | None = None
    last_cell: LammpsGuestRestartCell | None = None
    while index < len(lines):
        if lines[index].strip() != "ITEM: TIMESTEP":
            raise GuestRestartError(f"Unexpected LAMMPS dump format in {path}: missing ITEM: TIMESTEP")
        index += 2
        if index >= len(lines) or lines[index].strip() != "ITEM: NUMBER OF ATOMS":
            raise GuestRestartError(f"Unexpected LAMMPS dump format in {path}: missing atom count header")
        index += 1
        try:
            n_atoms = int(lines[index].strip())
        except (IndexError, ValueError) as exc:
            raise GuestRestartError(f"Invalid atom count in LAMMPS dump {path}.") from exc
        index += 1
        if index >= len(lines) or not lines[index].startswith("ITEM: BOX BOUNDS"):
            raise GuestRestartError(f"Unexpected LAMMPS dump format in {path}: missing box-bounds header")
        last_cell = _parse_lammps_dump_cell(lines[index : index + 4], dump_path=path)
        index += 4
        if index >= len(lines) or not lines[index].startswith("ITEM: ATOMS"):
            raise GuestRestartError(f"Unexpected LAMMPS dump format in {path}: missing atom header")
        header = lines[index].split()[2:]
        index += 1
        column_by_name = {name: column_index for column_index, name in enumerate(header)}
        required = ("id", "x", "y", "z")
        if any(name not in column_by_name for name in required):
            raise GuestRestartError(
                f"Unexpected LAMMPS dump atom columns in {path}: expected at least {required}, got {header}."
            )
        positions: dict[int, tuple[float, float, float]] = {}
        for _ in range(n_atoms):
            if index >= len(lines):
                raise GuestRestartError(f"Unexpected end of LAMMPS dump in {path}")
            parts = lines[index].split()
            index += 1
            try:
                atom_id = int(parts[column_by_name["id"]])
                positions[atom_id] = (
                    float(parts[column_by_name["x"]]),
                    float(parts[column_by_name["y"]]),
                    float(parts[column_by_name["z"]]),
                )
            except (IndexError, ValueError) as exc:
                raise GuestRestartError(f"Unexpected LAMMPS dump atom row in {path}: {lines[index - 1]!r}") from exc
        last_positions = positions
    if last_positions is None or last_cell is None:
        raise GuestRestartError(f"LAMMPS dump file did not contain any frames: {path}")
    return last_positions, last_cell


def write_graspa_restart_file(
    restart_state: LammpsGuestRestartState,
    output_path: str | Path,
    *,
    cell: LammpsGuestRestartCell,
    component_order: Sequence[str],
) -> GraspaRestartFileResult:
    path = Path(output_path).expanduser().resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    site_by_label = restart_state.site_by_label()
    template_by_component = restart_state.template_by_component()
    grouped = _group_guest_restart_atoms_by_component_and_molecule(restart_state)
    ordered_components = tuple(
        dict.fromkeys(
            [
                *(component for component in component_order if component in template_by_component),
                *(component for component in grouped if component not in component_order),
            ]
        )
    )
    n_molecules = sum(len(grouped.get(component, ())) for component in ordered_components)
    warnings: list[str] = []
    if n_molecules == 0:
        raise GuestRestartError("Cannot write a gRASPA restart file without adsorbate molecules.")

    lines: list[str] = []
    a_vec, b_vec, c_vec = cell.basis
    a_len, b_len, c_len = cell.lengths
    alpha, beta, gamma = cell.angles
    lines.extend(
        [
            "Cell info:",
            "========================================================================",
            f"number-of-unit-cells: {cell.unit_cells[0]} {cell.unit_cells[1]} {cell.unit_cells[2]}",
            "unit-cell-vector-a: " + " ".join(_format_restart_float(value) for value in a_vec),
            "unit-cell-vector-b: " + " ".join(_format_restart_float(value) for value in b_vec),
            "unit-cell-vector-c: " + " ".join(_format_restart_float(value) for value in c_vec),
            "",
            "cell-vector-a: " + " ".join(_format_restart_float(value) for value in a_vec),
            "cell-vector-b: " + " ".join(_format_restart_float(value) for value in b_vec),
            "cell-vector-c: " + " ".join(_format_restart_float(value) for value in c_vec),
            "",
            "cell-lengths: "
            + " ".join(_format_restart_float(value) for value in (a_len, b_len, c_len)),
            "cell-angles: "
            + " ".join(_format_restart_float(value) for value in (alpha, beta, gamma)),
            "",
            "",
            "Maximum changes for MC-moves:",
            "========================================================================",
            "Maximum-volume-change: 0.006250",
            "Maximum-Gibbs-volume-change: 0.025000",
            (
                "Maximum-box-shape-change: 0.100000 0.100000 0.100000, "
                "0.100000 0.100000 0.100000, 0.100000 0.100000 0.100000"
            ),
            "",
            "Acceptance targets for MC-moves:",
            "========================================================================",
            "Target-volume-change: 0.500000",
            "Target-box-shape-change: 0.500000",
            "Target-Gibbs-volume-change: 0.500000",
            "",
            f"Components: {len(ordered_components)} (Adsorbates {n_molecules}, Cations 0)",
            "========================================================================",
        ]
    )
    for component_index, component in enumerate(ordered_components):
        lines.extend(
            [
                f"Components {component_index} ({component})",
                "",
                (
                    f"Maximum-translation-change component {component_index}: "
                    f"{a_len / 10.0:.6f} {b_len / 10.0:.6f} {c_len / 10.0:.6f}"
                ),
                f"Maximum-translation-in-plane-change component {component_index}: 0.000000,0.000000,0.000000",
                f"Maximum-rotation-change component {component_index}: 0.523583 0.523583 0.523583",
                "",
            ]
        )
    lines.extend(["Reactions: 0", ""])
    for component_index, component in enumerate(ordered_components):
        molecules = grouped.get(component, ())
        template = template_by_component[component]
        lines.extend(
            [
                f"Component: {component_index}   Adsorbate {len(molecules)} molecules of {component}",
                "------------------------------------------------------------------------",
            ]
        )
        for molecule_index, molecule_atoms in enumerate(molecules):
            if len(molecule_atoms) != len(template.site_labels):
                warnings.append(
                    f"Component {component} molecule {molecule_index} has {len(molecule_atoms)} sites, "
                    f"expected {len(template.site_labels)}."
                )
            for atom in molecule_atoms:
                lines.append(
                    (
                        f"Adsorbate-atom-position: {molecule_index} {atom.site_index} "
                        f"{_format_restart_float(atom.x)} {_format_restart_float(atom.y)} {_format_restart_float(atom.z)}"
                    )
                )
        for molecule_index, molecule_atoms in enumerate(molecules):
            for atom in molecule_atoms:
                lines.append(f"Adsorbate-atom-velocity: {molecule_index} {atom.site_index} 0  0  0")
        for molecule_index, molecule_atoms in enumerate(molecules):
            for atom in molecule_atoms:
                lines.append(f"Adsorbate-atom-force: {molecule_index} {atom.site_index} 0  0  0")
        for molecule_index, molecule_atoms in enumerate(molecules):
            for atom in molecule_atoms:
                charge = site_by_label[atom.site_label].charge
                lines.append(
                    f"Adsorbate-atom-charge: {molecule_index} {atom.site_index} {_format_restart_float(charge)}"
                )
        for molecule_index, molecule_atoms in enumerate(molecules):
            for atom in molecule_atoms:
                lines.append(f"Adsorbate-atom-scaling: {molecule_index} {atom.site_index} 1")
        for molecule_index, molecule_atoms in enumerate(molecules):
            for atom in molecule_atoms:
                lines.append(f"Adsorbate-atom-fixed: {molecule_index} {atom.site_index} 0  0  0")
        lines.append("")
    path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")
    return GraspaRestartFileResult(
        restart_file_path=str(path),
        n_adsorbate_atoms=restart_state.n_atoms,
        n_adsorbate_molecules=n_molecules,
        components=ordered_components,
        warnings=tuple(dict.fromkeys(warnings)),
    )


def find_latest_gcmc_movie_snapshot(gcmc_result: object) -> Path | None:
    run_dirs: list[Path] = []
    for point in getattr(gcmc_result, "point_results", ()):
        pressure_run_dir = getattr(point, "pressure_run_dir", None)
        if pressure_run_dir is not None:
            run_dirs.append(Path(pressure_run_dir))
    for attr_name in ("widom_run_dir", "isotherm_root_dir", "mixture_root_dir", "output_dir"):
        attr_value = getattr(gcmc_result, attr_name, None)
        if attr_value is not None:
            run_dirs.append(Path(attr_value))

    candidates: list[Path] = []
    for run_dir in dict.fromkeys(path for path in run_dirs if path is not None):
        for root in (
            run_dir / "Movies" / "System_0",
            run_dir / "Movies",
            run_dir / "Movie" / "System_0",
            run_dir / "Movie",
        ):
            if root.is_dir():
                candidates.extend(path for path in root.rglob("*.data") if path.is_file())
    if not candidates:
        return None
    return max(candidates, key=_movie_sort_key)


def load_lammps_guest_force_field_assets(
    components: Sequence[str],
    *,
    guest_bundles: Sequence[str] = (),
) -> tuple[tuple[LammpsGuestTemplate, ...], tuple[LammpsGuestSite, ...]]:
    requested = tuple(dict.fromkeys(component.strip() for component in components if component.strip()))
    if not requested:
        raise GuestRestartError("At least one guest component is required for guest restart conversion.")

    bundles = _load_guest_bundle_catalog(guest_bundles)
    packaged_site_rows, packaged_mixing_rows = _load_packaged_pseudo_and_mixing_rows()
    templates: list[LammpsGuestTemplate] = []
    sites_by_label: dict[str, LammpsGuestSite] = {}
    for component in requested:
        bundle = bundles.get(component.casefold())
        if bundle is None:
            molecule_definition = _load_packaged_molecule_definition(component)
            pseudo_rows = packaged_site_rows
            mixing_rows = packaged_mixing_rows
            lammps_settings: Mapping[str, object] = {}
        else:
            molecule_definition = bundle.raspa.molecule_definition_text
            pseudo_rows = (*packaged_site_rows, *bundle.raspa.pseudo_atom_rows)
            mixing_rows = (*packaged_mixing_rows, *bundle.raspa.mixing_rule_rows)
            lammps_settings = bundle.lammps
        template = _parse_guest_template(component, molecule_definition, lammps_settings=lammps_settings)
        templates.append(template)
        row_sites = _parse_guest_sites(
            template.site_labels,
            pseudo_rows=pseudo_rows,
            mixing_rows=mixing_rows,
            lammps_settings=lammps_settings,
        )
        for site in row_sites:
            previous = sites_by_label.get(site.label)
            if previous is not None and previous != site:
                raise GuestRestartError(f"Conflicting guest force-field rows for pseudo atom {site.label!r}.")
            sites_by_label[site.label] = site

    for site in sites_by_label.values():
        if site.mass <= 0.0:
            raise GuestRestartError(
                f"Guest pseudo atom {site.label!r} has zero or negative mass. "
                "LAMMPS guest restart currently requires massive interaction sites."
            )
        if site.epsilon_k < 0.0 or site.sigma < 0.0:
            raise GuestRestartError(f"Guest pseudo atom {site.label!r} has invalid Lennard-Jones parameters.")
    return tuple(templates), tuple(sites_by_label.values())


def parse_lammps_guest_restart_snapshot(
    snapshot_path: str | Path,
    *,
    templates: Sequence[LammpsGuestTemplate],
    sites: Sequence[LammpsGuestSite],
) -> LammpsGuestRestartState:
    path = Path(snapshot_path).expanduser().resolve()
    if not path.is_file():
        raise GuestRestartError(f"Guest restart snapshot does not exist: {path}")
    text = path.read_text(encoding="utf-8", errors="replace")
    known_sites = {site.label for site in sites}
    if not known_sites:
        raise GuestRestartError("No guest force-field sites were provided for snapshot parsing.")
    type_labels = _parse_lammps_mass_type_labels(text)
    raw_atoms = _parse_lammps_data_atom_rows(text, type_labels=type_labels)

    site_to_components: dict[str, list[str]] = {}
    template_by_component = {template.component: template for template in templates}
    for template in templates:
        for site_label in template.site_labels:
            site_to_components.setdefault(site_label, []).append(template.component)

    grouped: dict[tuple[str, str], list[tuple[str, float, float, float]]] = {}
    warnings: list[str] = []
    for row in raw_atoms:
        label = row["label"]
        if label not in known_sites:
            continue
        component_candidates = site_to_components.get(label, [])
        if len(component_candidates) != 1:
            warnings.append(
                f"Snapshot atom with pseudo atom {label!r} could not be assigned uniquely to a component; skipped."
            )
            continue
        component = component_candidates[0]
        molecule_key = f"{component}:{row['molecule_id']}"
        grouped.setdefault((component, molecule_key), []).append((label, row["x"], row["y"], row["z"]))

    atoms: list[LammpsGuestSnapshotAtom] = []
    for (component, molecule_key), rows in grouped.items():
        template = template_by_component[component]
        if len(rows) != len(template.site_labels):
            warnings.append(
                f"Skipped incomplete {component} molecule {molecule_key!r}: "
                f"found {len(rows)} sites, expected {len(template.site_labels)}."
            )
            continue
        used = [False] * len(rows)
        ordered_rows: list[tuple[str, float, float, float]] = []
        for site_label in template.site_labels:
            match_index = next(
                (index for index, row in enumerate(rows) if not used[index] and row[0] == site_label),
                None,
            )
            if match_index is None:
                warnings.append(
                    f"Skipped {component} molecule {molecule_key!r}: missing expected site {site_label!r}."
                )
                ordered_rows = []
                break
            used[match_index] = True
            ordered_rows.append(rows[match_index])
        for site_index, (site_label, x, y, z) in enumerate(ordered_rows):
            atoms.append(
                LammpsGuestSnapshotAtom(
                    component=component,
                    molecule_key=molecule_key,
                    site_index=site_index,
                    site_label=site_label,
                    x=x,
                    y=y,
                    z=z,
                )
            )

    if not atoms:
        raise GuestRestartError(
            f"No guest atoms from the requested components were parsed from snapshot {path}."
        )
    return LammpsGuestRestartState(
        source_snapshot_path=str(path),
        atoms=tuple(atoms),
        sites=tuple(sites),
        templates=tuple(templates),
        warnings=tuple(dict.fromkeys(warnings)),
    )


def _load_guest_bundle_catalog(paths: Sequence[str]) -> dict[str, GuestBundle]:
    try:
        bundles = load_guest_bundles(paths)
    except GuestBundleError as exc:
        raise GuestRestartError(str(exc)) from exc
    catalog: dict[str, GuestBundle] = {}
    for bundle in bundles:
        for name in bundle.names():
            catalog[name.casefold()] = bundle
    return catalog


def _load_packaged_molecule_definition(component: str) -> str:
    path = _packaged_guest_template_dir() / f"{component}.def"
    if not path.is_file():
        raise GuestRestartError(
            f"Packaged guest component {component!r} is not available for LAMMPS guest restart. "
            "Provide an explicit --guest-bundle with RASPA and LAMMPS guest data."
        )
    return path.read_text(encoding="utf-8")


def _load_packaged_pseudo_and_mixing_rows() -> tuple[tuple[str, ...], tuple[str, ...]]:
    template_dir = _packaged_guest_template_dir()
    pseudo_path = template_dir / "pseudo_atoms.def"
    mixing_path = template_dir / "force_field_mixing_rules.def"
    if not pseudo_path.is_file() or not mixing_path.is_file():
        raise GuestRestartError(f"Packaged guest force-field assets are missing from {template_dir}.")
    return _extract_data_rows(pseudo_path.read_text(encoding="utf-8")), _extract_data_rows(
        mixing_path.read_text(encoding="utf-8")
    )


def _packaged_guest_template_dir() -> Path:
    return Path(__file__).resolve().parent / "data" / "graspa" / "widom_template"


def _parse_guest_template(
    component: str,
    molecule_definition: str,
    *,
    lammps_settings: Mapping[str, object],
) -> LammpsGuestTemplate:
    positions: list[tuple[str, tuple[float, float, float]]] = []
    bonds: list[tuple[int, int]] = []
    in_positions = False
    in_bonds = False
    for raw_line in molecule_definition.splitlines():
        line = _strip_inline_comment(raw_line).strip()
        lowered = line.casefold()
        if not line:
            continue
        if "atomic positions" in lowered:
            in_positions = True
            in_bonds = False
            continue
        if "bond stretch" in lowered:
            in_positions = False
            in_bonds = True
            continue
        if line.startswith("#"):
            if in_positions and positions:
                in_positions = False
            elif in_bonds and bonds:
                in_bonds = False
            continue
        parts = line.split()
        if in_positions and len(parts) >= 5 and _is_int_token(parts[0]):
            positions.append((parts[1], (float(parts[2]), float(parts[3]), float(parts[4]))))
            continue
        if in_bonds and len(parts) >= 2 and _is_int_token(parts[0]) and _is_int_token(parts[1]):
            bonds.append((int(parts[0]), int(parts[1])))

    if not positions:
        raise GuestRestartError(f"No atomic positions were parsed from molecule definition for {component!r}.")
    site_labels = tuple(label for label, _position in positions)
    site_order = lammps_settings.get("site_order")
    if site_order is not None:
        if not isinstance(site_order, Sequence) or isinstance(site_order, (str, bytes)):
            raise GuestRestartError(f"Guest bundle lammps.site_order for {component!r} must be a list of strings.")
        normalized_order = tuple(str(label) for label in site_order)
        if sorted(normalized_order) != sorted(site_labels):
            raise GuestRestartError(
                f"Guest bundle lammps.site_order for {component!r} must match the RASPA molecule sites."
            )
        position_by_label = {label: position for label, position in positions}
        site_labels = normalized_order
        positions_tuple = tuple(position_by_label[label] for label in site_labels)
    else:
        positions_tuple = tuple(position for _label, position in positions)

    return LammpsGuestTemplate(
        component=component,
        site_labels=site_labels,
        relative_positions=positions_tuple,
        bonds=tuple(bonds),
        bond_force_constant=float(lammps_settings.get("bond_force_constant", 1000.0)),
        angle_force_constant=float(lammps_settings.get("angle_force_constant", 100.0)),
    )


def _parse_guest_sites(
    site_labels: Sequence[str],
    *,
    pseudo_rows: Sequence[str],
    mixing_rows: Sequence[str],
    lammps_settings: Mapping[str, object],
) -> tuple[LammpsGuestSite, ...]:
    pseudo_by_label = _parse_pseudo_atom_rows(pseudo_rows)
    mixing_by_label = _parse_mixing_rule_rows(mixing_rows)
    lammps_masses = _mapping_of_floats(lammps_settings.get("masses"), "masses")
    lammps_charges = _mapping_of_floats(lammps_settings.get("charges"), "charges")
    lammps_pair_rows = _parse_mixing_rule_rows(tuple(str(row) for row in lammps_settings.get("pair_coeff_rows", ())))
    sites: list[LammpsGuestSite] = []
    for label in site_labels:
        pseudo = pseudo_by_label.get(label)
        if pseudo is None:
            raise GuestRestartError(f"Missing pseudo-atom row for guest site {label!r}.")
        mixing = lammps_pair_rows.get(label) or mixing_by_label.get(label)
        if mixing is None:
            raise GuestRestartError(f"Missing Lennard-Jones row for guest site {label!r}.")
        sites.append(
            LammpsGuestSite(
                label=label,
                element=pseudo["element"],
                mass=lammps_masses.get(label, pseudo["mass"]),
                charge=lammps_charges.get(label, pseudo["charge"]),
                epsilon_k=mixing["epsilon_k"],
                sigma=mixing["sigma"],
            )
        )
    return tuple(sites)


def _parse_pseudo_atom_rows(rows: Sequence[str]) -> dict[str, dict[str, float | str]]:
    parsed: dict[str, dict[str, float | str]] = {}
    for row in rows:
        parts = _strip_inline_comment(row).split()
        if len(parts) < 7 or parts[0].startswith("#"):
            continue
        try:
            parsed[parts[0]] = {
                "element": parts[2],
                "mass": float(parts[5]),
                "charge": float(parts[6]),
            }
        except ValueError:
            continue
    return parsed


def _parse_mixing_rule_rows(rows: Sequence[str]) -> dict[str, dict[str, float]]:
    parsed: dict[str, dict[str, float]] = {}
    for row in rows:
        parts = _strip_inline_comment(row).split()
        if len(parts) < 4 or parts[0].startswith("#"):
            continue
        interaction = parts[1].casefold()
        if interaction not in {"lennard-jones", "feynman-hibbs-lennard-jones"}:
            continue
        try:
            parsed[parts[0]] = {"epsilon_k": float(parts[2]), "sigma": float(parts[3])}
        except ValueError:
            continue
    return parsed


def _parse_lammps_mass_type_labels(text: str) -> dict[int, str]:
    masses = _section_rows(text, "Masses")
    labels: dict[int, str] = {}
    for row in masses:
        before, _hash, after = row.partition("#")
        parts = before.split()
        if len(parts) < 2 or not _is_int_token(parts[0]):
            continue
        comment_parts = after.split()
        if comment_parts:
            labels[int(parts[0])] = comment_parts[0]
    return labels


def _parse_lammps_data_atom_rows(text: str, *, type_labels: Mapping[int, str]) -> list[dict[str, object]]:
    section_header, rows = _atom_section_rows(text)
    atom_style = "full" if "full" in section_header.casefold() else "molecular"
    parsed: list[dict[str, object]] = []
    for row in rows:
        before, _hash, after = row.partition("#")
        parts = before.split()
        if len(parts) < 6 or not _is_int_token(parts[0]) or not _is_int_token(parts[1]):
            continue
        type_id = int(parts[2]) if _is_int_token(parts[2]) else None
        if atom_style == "full":
            if len(parts) < 7:
                continue
            coord_offset = 4
        else:
            coord_offset = 3
        try:
            x = float(parts[coord_offset])
            y = float(parts[coord_offset + 1])
            z = float(parts[coord_offset + 2])
        except ValueError:
            continue
        label = None
        component = None
        comment_parts = after.split()
        if comment_parts:
            label = comment_parts[0]
            if len(comment_parts) > 1:
                component = comment_parts[1]
        if label is None and type_id is not None:
            label = type_labels.get(type_id)
        if label is None:
            tail = parts[-1]
            if not _is_float_token(tail):
                label = tail
        if label is None:
            continue
        parsed.append(
            {
                "atom_id": int(parts[0]),
                "molecule_id": parts[1],
                "type_id": type_id,
                "label": label,
                "component": component,
                "x": x,
                "y": y,
                "z": z,
            }
        )
    return parsed


def _parse_lammps_guest_atom_identities_from_data(
    text: str,
    *,
    previous_guest_restart_state: LammpsGuestRestartState,
) -> tuple[dict[str, object], ...]:
    known_sites = {site.label for site in previous_guest_restart_state.sites}
    type_labels = _parse_lammps_mass_type_labels(text)
    raw_atoms = _parse_lammps_data_atom_rows(text, type_labels=type_labels)
    template_by_component = previous_guest_restart_state.template_by_component()
    site_to_components: dict[str, list[str]] = {}
    for template in previous_guest_restart_state.templates:
        for site_label in template.site_labels:
            site_to_components.setdefault(site_label, []).append(template.component)

    grouped: dict[tuple[str, str], list[dict[str, object]]] = {}
    for row in raw_atoms:
        site_label = str(row["label"])
        if site_label not in known_sites:
            continue
        component = row.get("component")
        if isinstance(component, str) and component in template_by_component:
            resolved_component = component
        else:
            candidates = site_to_components.get(site_label, [])
            if len(candidates) != 1:
                continue
            resolved_component = candidates[0]
        molecule_id = str(row["molecule_id"])
        grouped.setdefault((resolved_component, molecule_id), []).append(row)

    identities: list[dict[str, object]] = []
    for (component, molecule_id), rows in grouped.items():
        template = template_by_component[component]
        ordered_rows = sorted(rows, key=lambda row: int(row["atom_id"]))
        if len(ordered_rows) != len(template.site_labels):
            continue
        for site_index, row in enumerate(ordered_rows):
            expected_label = template.site_labels[site_index]
            if str(row["label"]) != expected_label:
                raise GuestRestartError(
                    f"LAMMPS MD data guest atom order for component {component!r}, molecule {molecule_id!r} "
                    f"does not match the molecule definition: expected {expected_label!r}, got {row['label']!r}."
                )
            identities.append(
                {
                    "atom_id": int(row["atom_id"]),
                    "component": component,
                    "molecule_key": f"{component}:md:{molecule_id}",
                    "site_index": site_index,
                    "site_label": expected_label,
                }
            )
    if not identities:
        raise GuestRestartError("No guest atom identities were recovered from the LAMMPS MD data file.")
    return tuple(identities)


def _atom_section_rows(text: str) -> tuple[str, list[str]]:
    lines = text.splitlines()
    for index, line in enumerate(lines):
        if line.strip().startswith("Atoms"):
            rows: list[str] = []
            cursor = index + 1
            while cursor < len(lines) and not lines[cursor].strip():
                cursor += 1
            while cursor < len(lines):
                candidate = lines[cursor]
                if _is_lammps_section_header(candidate):
                    break
                if candidate.strip():
                    rows.append(candidate)
                cursor += 1
            return line.strip(), rows
    raise GuestRestartError("LAMMPS-style guest snapshot does not contain an Atoms section.")


def _section_rows(text: str, section_name: str) -> list[str]:
    lines = text.splitlines()
    for index, line in enumerate(lines):
        if line.strip() == section_name:
            rows: list[str] = []
            cursor = index + 1
            while cursor < len(lines) and not lines[cursor].strip():
                cursor += 1
            while cursor < len(lines):
                candidate = lines[cursor]
                if _is_lammps_section_header(candidate):
                    break
                if candidate.strip():
                    rows.append(candidate)
                cursor += 1
            return rows
    return []


def _parse_lammps_dump_cell(lines: Sequence[str], *, dump_path: Path) -> LammpsGuestRestartCell:
    if len(lines) != 4:
        raise GuestRestartError(f"Unexpected LAMMPS dump box block in {dump_path}: expected 4 lines, got {len(lines)}")
    header = lines[0].split()
    if header[:3] != ["ITEM:", "BOX", "BOUNDS"]:
        raise GuestRestartError(f"Unexpected LAMMPS dump box header in {dump_path}: {lines[0]!r}")
    try:
        bounds = [tuple(float(value) for value in line.split()) for line in lines[1:4]]
    except ValueError as exc:
        raise GuestRestartError(f"Invalid numeric box-bounds line in {dump_path}") from exc
    triclinic = "xy" in header and "xz" in header and "yz" in header
    if triclinic:
        if any(len(row) != 3 for row in bounds):
            raise GuestRestartError(f"Unexpected triclinic box-bounds format in {dump_path}")
        (xlo_bound, xhi_bound, xy), (ylo_bound, yhi_bound, xz), (zlo_bound, zhi_bound, yz) = bounds
        xlo = xlo_bound - min(0.0, xy, xz, xy + xz)
        xhi = xhi_bound - max(0.0, xy, xz, xy + xz)
        ylo = ylo_bound - min(0.0, yz)
        yhi = yhi_bound - max(0.0, yz)
        zlo = zlo_bound
        zhi = zhi_bound
        _ = (xlo, ylo, zlo)
        basis = (
            (xhi - xlo, 0.0, 0.0),
            (xy, yhi - ylo, 0.0),
            (xz, yz, zhi - zlo),
        )
    else:
        if any(len(row) != 2 for row in bounds):
            raise GuestRestartError(f"Unexpected orthogonal box-bounds format in {dump_path}")
        (xlo, xhi), (ylo, yhi), (zlo, zhi) = bounds
        _ = (xlo, ylo, zlo)
        basis = (
            (xhi - xlo, 0.0, 0.0),
            (0.0, yhi - ylo, 0.0),
            (0.0, 0.0, zhi - zlo),
        )
    lengths = tuple(_vector_norm(vector) for vector in basis)
    if any(length <= 1.0e-12 for length in lengths):
        raise GuestRestartError(f"LAMMPS dump box in {dump_path} is singular.")
    angles = (
        _angle_degrees(basis[1], basis[2]),
        _angle_degrees(basis[0], basis[2]),
        _angle_degrees(basis[0], basis[1]),
    )
    return LammpsGuestRestartCell(
        unit_cells=(1, 1, 1),
        basis=basis,
        lengths=(float(lengths[0]), float(lengths[1]), float(lengths[2])),
        angles=angles,
    )


def _group_guest_restart_atoms_by_component_and_molecule(
    restart_state: LammpsGuestRestartState,
) -> dict[str, tuple[tuple[LammpsGuestSnapshotAtom, ...], ...]]:
    grouped: dict[str, dict[str, list[LammpsGuestSnapshotAtom]]] = {}
    for atom in restart_state.atoms:
        grouped.setdefault(atom.component, {}).setdefault(atom.molecule_key, []).append(atom)
    result: dict[str, tuple[tuple[LammpsGuestSnapshotAtom, ...], ...]] = {}
    for component, molecule_map in grouped.items():
        result[component] = tuple(
            tuple(sorted(atoms, key=lambda atom: atom.site_index))
            for _molecule_key, atoms in sorted(molecule_map.items())
        )
    return result


_LAMMPS_SECTION_NAMES = {
    "Masses",
    "Pair Coeffs",
    "PairIJ Coeffs",
    "Bond Coeffs",
    "Angle Coeffs",
    "Dihedral Coeffs",
    "Improper Coeffs",
    "Atoms",
    "Bonds",
    "Angles",
    "Dihedrals",
    "Impropers",
    "Velocities",
}


def _is_lammps_section_header(line: str) -> bool:
    stripped = line.strip()
    if not stripped:
        return False
    return any(stripped == name or stripped.startswith(name + " #") for name in _LAMMPS_SECTION_NAMES)


def _extract_data_rows(text: str) -> tuple[str, ...]:
    rows: list[str] = []
    for raw_line in text.splitlines():
        line = _strip_inline_comment(raw_line).strip()
        if not line or line.startswith("#") or line.startswith("//"):
            continue
        rows.append(line)
    return tuple(rows)


def _strip_inline_comment(line: str) -> str:
    return line.split("//", 1)[0].strip()


def _mapping_of_floats(raw: object, field_name: str) -> dict[str, float]:
    if raw is None:
        return {}
    if not isinstance(raw, Mapping):
        raise GuestRestartError(f"Guest bundle lammps.{field_name} must be an object mapping site labels to numbers.")
    result: dict[str, float] = {}
    for key, value in raw.items():
        try:
            result[str(key)] = float(value)
        except (TypeError, ValueError) as exc:
            raise GuestRestartError(f"Guest bundle lammps.{field_name}.{key} must be numeric.") from exc
    return result


def _movie_sort_key(path: Path) -> tuple[int, str]:
    match = re.search(r"(\d+)(?=\.data$)", path.name)
    return (int(match.group(1)) if match else -1, str(path))


def _format_restart_float(value: float) -> str:
    return f"{float(value):.6f}"


def _vector_norm(vector: tuple[float, float, float]) -> float:
    return math.sqrt(sum(value * value for value in vector))


def _dot(left: tuple[float, float, float], right: tuple[float, float, float]) -> float:
    return sum(left[index] * right[index] for index in range(3))


def _angle_degrees(left: tuple[float, float, float], right: tuple[float, float, float]) -> float:
    left_norm = _vector_norm(left)
    right_norm = _vector_norm(right)
    if left_norm <= 1.0e-12 or right_norm <= 1.0e-12:
        return 0.0
    cosine = max(-1.0, min(1.0, _dot(left, right) / (left_norm * right_norm)))
    return math.degrees(math.acos(cosine))


def _is_int_token(value: str) -> bool:
    try:
        int(value)
    except ValueError:
        return False
    return True


def _is_float_token(value: str) -> bool:
    try:
        number = float(value)
    except ValueError:
        return False
    return math.isfinite(number)


__all__ = [
    "GuestRestartError",
    "GraspaRestartFileResult",
    "KCAL_PER_MOL_PER_K",
    "LammpsGuestRestartCell",
    "LammpsGuestRestartState",
    "LammpsGuestSite",
    "LammpsGuestSnapshotAtom",
    "LammpsGuestTemplate",
    "build_lammps_guest_restart_state_from_lammps_md_result",
    "build_lammps_guest_restart_state_from_gcmc_result",
    "find_latest_gcmc_movie_snapshot",
    "load_lammps_guest_force_field_assets",
    "parse_lammps_md_dump_guest_positions",
    "parse_lammps_guest_restart_snapshot",
    "write_graspa_restart_file",
]
