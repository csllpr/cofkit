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
        comment_parts = after.split()
        if comment_parts:
            label = comment_parts[0]
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
                "x": x,
                "y": y,
                "z": z,
            }
        )
    return parsed


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
    "KCAL_PER_MOL_PER_K",
    "LammpsGuestRestartState",
    "LammpsGuestSite",
    "LammpsGuestSnapshotAtom",
    "LammpsGuestTemplate",
    "build_lammps_guest_restart_state_from_gcmc_result",
    "find_latest_gcmc_movie_snapshot",
    "load_lammps_guest_force_field_assets",
    "parse_lammps_guest_restart_snapshot",
]
