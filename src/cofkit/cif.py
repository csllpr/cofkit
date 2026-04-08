from __future__ import annotations

from dataclasses import dataclass, field
from math import acos, degrees, floor
from dataclasses import replace
from pathlib import Path
from typing import Iterable, Mapping

from .bond_types import bond_order_to_cif_type, normalize_bond_order
from .chem.motif_registry import motif_pseudo_atom_symbol
from .cofid import cofid_comment_line
from .geometry import Vec3, add, cross, dot, matmul_vec, norm, scale
from .model import Candidate, MonomerSpec, Pose
from .reaction_realization import RealizedBond, ReactionRealizationResult, ReactionRealizer


@dataclass(frozen=True)
class AtomSite:
    label: str
    type_symbol: str
    fract_x: float
    fract_y: float
    fract_z: float
    occupancy: float = 1.0


@dataclass(frozen=True)
class CIFExportResult:
    text: str
    mode: str
    n_sites: int
    data_name: str
    metadata: Mapping[str, object] = field(default_factory=dict)


class CIFWriter:
    """Writes legal P1 CIF files from a cofkit candidate.

    Export modes:
    - atomistic: all monomers provide atom_symbols + atom_positions
    - mixed: some monomers provide atomistic coordinates, others fall back to coarse pseudo-sites
    - coarse: no monomer provides atomistic coordinates, so monomer centers / motif origins are exported
    """

    def __init__(self, reaction_realizer: ReactionRealizer | None = None) -> None:
        self.reaction_realizer = reaction_realizer or ReactionRealizer()

    def export_candidate(
        self,
        candidate: Candidate,
        monomer_specs: Mapping[str, MonomerSpec] | Iterable[MonomerSpec],
        data_name: str | None = None,
        cofid: str | None = None,
        cofid_comment_suffix: str | None = None,
    ) -> CIFExportResult:
        spec_map = self._coerce_spec_map(monomer_specs)
        instance_to_monomer = self._instance_to_monomer(candidate)
        cell = candidate.state.cell
        sites, bonds, mode, export_metadata = self._build_sites(candidate, spec_map, instance_to_monomer)
        rendered = self._render_cif(
            data_name=self._sanitize_data_name(data_name or candidate.id),
            cell=cell,
            sites=sites,
            bonds=bonds,
            mode=mode,
            metadata=export_metadata,
            cofid=cofid,
            cofid_comment_suffix=cofid_comment_suffix,
        )
        return CIFExportResult(
            text=rendered,
            mode=mode,
            n_sites=len(sites),
            data_name=self._sanitize_data_name(data_name or candidate.id),
            metadata=export_metadata,
        )

    def write_candidate(
        self,
        path: str | Path,
        candidate: Candidate,
        monomer_specs: Mapping[str, MonomerSpec] | Iterable[MonomerSpec],
        data_name: str | None = None,
        cofid: str | None = None,
        cofid_comment_suffix: str | None = None,
    ) -> CIFExportResult:
        result = self.export_candidate(
            candidate,
            monomer_specs,
            data_name=data_name,
            cofid=cofid,
            cofid_comment_suffix=cofid_comment_suffix,
        )
        out_path = Path(path)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(result.text)
        return result

    def _build_sites(
        self,
        candidate: Candidate,
        monomer_specs: Mapping[str, MonomerSpec],
        instance_to_monomer: Mapping[str, str],
    ) -> tuple[list[AtomSite], list[RealizedBond], str, dict[str, object]]:
        cell = candidate.state.cell
        sites: list[AtomSite] = []
        bonds: list[RealizedBond] = []
        atomistic_instances = 0
        coarse_instances = 0
        realization = self.reaction_realizer.realize(candidate, monomer_specs, instance_to_monomer)
        atomistic_shift = self._atomistic_c_axis_shift(candidate, monomer_specs, instance_to_monomer, realization)

        for instance_id, pose in candidate.state.monomer_poses.items():
            monomer_id = instance_to_monomer[instance_id]
            monomer = monomer_specs[monomer_id]
            if monomer.atom_symbols and monomer.atom_positions:
                atomistic_instances += 1
                sites.extend(self._atomistic_sites(instance_id, pose, monomer, cell, realization, atomistic_shift))
            else:
                coarse_instances += 1
                sites.extend(self._coarse_sites(instance_id, pose, monomer, cell))

        if atomistic_instances and coarse_instances:
            mode = "mixed"
        elif atomistic_instances:
            mode = "atomistic"
        else:
            mode = "coarse"

        metadata = {
            "atomistic_instances": atomistic_instances,
            "coarse_instances": coarse_instances,
            "n_monomer_instances": len(candidate.state.monomer_poses),
        }
        bonds.extend(self._atomistic_bonds(candidate, monomer_specs, instance_to_monomer, realization, cell, atomistic_shift))
        if realization is not None:
            metadata["reaction_realization"] = dict(realization.metadata)
        return sites, bonds, mode, metadata

    def _atomistic_sites(
        self,
        instance_id: str,
        pose: Pose,
        monomer: MonomerSpec,
        cell: tuple[Vec3, Vec3, Vec3],
        realization: ReactionRealizationResult | None,
        atomistic_shift: Vec3,
    ) -> list[AtomSite]:
        sites: list[AtomSite] = []
        realized_atoms = None if realization is None else realization.atoms_by_instance.get(instance_id)
        atom_rows = (
            ((atom.label, atom.symbol, atom.local_position) for atom in realized_atoms)
            if realized_atoms is not None
            else (
                (
                    f"{instance_id}_{symbol}{idx}",
                    symbol,
                    local_pos,
                )
                for idx, (symbol, local_pos) in enumerate(zip(monomer.atom_symbols, monomer.atom_positions), start=1)
            )
        )
        for label, symbol, local_pos in atom_rows:
            world = add(self._world_position(pose, local_pos), atomistic_shift)
            frac = self._wrap_fractional(self._cartesian_to_fractional(cell, world))
            sites.append(
                AtomSite(
                    label=label,
                    type_symbol=symbol,
                    fract_x=frac[0],
                    fract_y=frac[1],
                    fract_z=frac[2],
                )
            )
        return sites

    def _coarse_sites(
        self,
        instance_id: str,
        pose: Pose,
        monomer: MonomerSpec,
        cell: tuple[Vec3, Vec3, Vec3],
    ) -> list[AtomSite]:
        sites: list[AtomSite] = []
        center_frac = self._wrap_fractional(self._cartesian_to_fractional(cell, pose.translation))
        sites.append(
            AtomSite(
                label=f"{instance_id}_CTR",
                type_symbol="C",
                fract_x=center_frac[0],
                fract_y=center_frac[1],
                fract_z=center_frac[2],
            )
        )
        for idx, motif in enumerate(monomer.motifs, start=1):
            world = self._world_position(pose, motif.frame.origin)
            frac = self._wrap_fractional(self._cartesian_to_fractional(cell, world))
            sites.append(
                AtomSite(
                    label=f"{instance_id}_M{idx}",
                    type_symbol=self._motif_symbol(motif.kind),
                    fract_x=frac[0],
                    fract_y=frac[1],
                    fract_z=frac[2],
                )
            )
        return sites

    def _render_cif(
        self,
        data_name: str,
        cell: tuple[Vec3, Vec3, Vec3],
        sites: list[AtomSite],
        bonds: list[RealizedBond],
        mode: str,
        metadata: Mapping[str, object],
        cofid: str | None = None,
        cofid_comment_suffix: str | None = None,
    ) -> str:
        a, b, c, alpha, beta, gamma = self._cell_parameters(cell)
        lines = [
            *([cofid_comment_line(cofid, suffix=cofid_comment_suffix)] if cofid else []),
            f"data_{data_name}",
            "# CIF generated by cofkit",
            f"# export_mode: {mode}",
            f"# atomistic_instances: {metadata.get('atomistic_instances', 0)}",
            f"# coarse_instances: {metadata.get('coarse_instances', 0)}",
            "_audit_creation_method 'cofkit CIFWriter'",
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
        for site in sites:
            lines.append(
                f"{site.label} {site.type_symbol} {site.fract_x:.6f} {site.fract_y:.6f} {site.fract_z:.6f} {site.occupancy:.2f}"
            )
        reaction_realization = metadata.get("reaction_realization")
        if isinstance(reaction_realization, Mapping):
            notes = reaction_realization.get("notes", ())
            if notes:
                lines.append("")
                for note in notes:
                    lines.append(f"# reaction_realization_note: {note}")
        if bonds:
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
            for bond in bonds:
                lines.append(
                    f"{bond.label_1} {bond.label_2} {bond.symmetry_1} {bond.symmetry_2} {bond.distance:.6f} "
                    f"{bond_order_to_cif_type(bond.bond_order)}"
                )
        lines.append("")
        return "\n".join(lines)

    def _atomistic_bonds(
        self,
        candidate: Candidate,
        monomer_specs: Mapping[str, MonomerSpec],
        instance_to_monomer: Mapping[str, str],
        realization: ReactionRealizationResult | None,
        cell: tuple[Vec3, Vec3, Vec3],
        atomistic_shift: Vec3,
    ) -> list[RealizedBond]:
        bonds: list[RealizedBond] = []
        exported_labels: set[str] = set()
        fractional_by_label: dict[str, Vec3] = {}
        removed_atom_ids_by_instance = {} if realization is None else realization.removed_atom_ids_by_instance

        for instance_id, monomer_id in instance_to_monomer.items():
            monomer = monomer_specs[monomer_id]
            if not (monomer.atom_symbols and monomer.atom_positions):
                continue
            realized_atoms = None if realization is None else realization.atoms_by_instance.get(instance_id)
            if realized_atoms is not None:
                atom_labels = {
                    atom.atom_id: atom.label
                    for atom in realized_atoms
                }
                atom_positions = {
                    atom.atom_id: atom.local_position
                    for atom in realized_atoms
                }
            else:
                removed_atom_ids = set(removed_atom_ids_by_instance.get(instance_id, ()))
                atom_labels = {
                    atom_id: ReactionRealizer.atom_label(instance_id, symbol, atom_id)
                    for atom_id, symbol in enumerate(monomer.atom_symbols)
                    if atom_id not in removed_atom_ids
                }
                atom_positions = {
                    atom_id: monomer.atom_positions[atom_id]
                    for atom_id in atom_labels
                }
            exported_labels.update(atom_labels.values())
            pose = candidate.state.monomer_poses[instance_id]
            for atom_id, label in atom_labels.items():
                world = add(self._world_position(pose, atom_positions[atom_id]), atomistic_shift)
                fractional_by_label[label] = self._wrap_fractional(self._cartesian_to_fractional(cell, world))
            for atom_id_1, atom_id_2, bond_order in monomer.bonds:
                if atom_id_1 not in atom_labels or atom_id_2 not in atom_labels:
                    continue
                label_1 = atom_labels[atom_id_1]
                label_2 = atom_labels[atom_id_2]
                relative_shift, distance = self._minimum_image_bond_geometry(
                    cell,
                    fractional_by_label[label_1],
                    fractional_by_label[label_2],
                )
                bonds.append(
                    RealizedBond(
                        label_1=label_1,
                        label_2=label_2,
                        distance=distance,
                        symmetry_1=".",
                        symmetry_2=self._format_p1_symmetry_shift(relative_shift),
                        bond_order=normalize_bond_order(bond_order),
                    )
                )

        if realization is not None:
            for bond in realization.bonds:
                if bond.label_1 in exported_labels and bond.label_2 in exported_labels:
                    fractional_1 = fractional_by_label.get(bond.label_1)
                    fractional_2 = fractional_by_label.get(bond.label_2)
                    if fractional_1 is None or fractional_2 is None:
                        bonds.append(bond)
                        continue
                    relative_shift, distance = self._minimum_image_bond_geometry(
                        cell,
                        fractional_1,
                        fractional_2,
                    )
                    bonds.append(
                        replace(
                            bond,
                            distance=distance,
                            symmetry_1=".",
                            symmetry_2=self._format_p1_symmetry_shift(relative_shift),
                        )
                    )
        return self._dedupe_bonds(bonds)

    def _atomistic_c_axis_shift(
        self,
        candidate: Candidate,
        monomer_specs: Mapping[str, MonomerSpec],
        instance_to_monomer: Mapping[str, str],
        realization: ReactionRealizationResult | None,
    ) -> Vec3:
        fractional_z_values: list[float] = []
        cell = candidate.state.cell
        for instance_id, pose in candidate.state.monomer_poses.items():
            monomer = monomer_specs[instance_to_monomer[instance_id]]
            if not (monomer.atom_symbols and monomer.atom_positions):
                continue
            realized_atoms = None if realization is None else realization.atoms_by_instance.get(instance_id)
            atom_positions = (
                (atom.local_position for atom in realized_atoms)
                if realized_atoms is not None
                else monomer.atom_positions
            )
            for local_position in atom_positions:
                world = self._world_position(pose, local_position)
                fractional_z_values.append(self._cartesian_to_fractional(cell, world)[2])
        if not fractional_z_values:
            return (0.0, 0.0, 0.0)
        min_z = min(fractional_z_values)
        max_z = max(fractional_z_values)
        if max_z - min_z >= 0.5:
            return (0.0, 0.0, 0.0)
        shift_fractional_z = 0.5 - 0.5 * (min_z + max_z)
        return (
            cell[2][0] * shift_fractional_z,
            cell[2][1] * shift_fractional_z,
            cell[2][2] * shift_fractional_z,
        )

    def _dedupe_bonds(self, bonds: list[RealizedBond]) -> list[RealizedBond]:
        unique: list[RealizedBond] = []
        index_by_key: dict[tuple[str, str, str, str], int] = {}
        for bond in bonds:
            key = (bond.label_1, bond.label_2, bond.symmetry_1, bond.symmetry_2)
            reverse_key = (bond.label_2, bond.label_1, bond.symmetry_2, bond.symmetry_1)
            existing_index = index_by_key.get(key)
            if existing_index is None:
                existing_index = index_by_key.get(reverse_key)
            if existing_index is None:
                index_by_key[key] = len(unique)
                unique.append(bond)
                continue
            unique[existing_index] = bond
            index_by_key[key] = existing_index
        return unique

    def _cell_parameters(self, cell: tuple[Vec3, Vec3, Vec3]) -> tuple[float, float, float, float, float, float]:
        a_vec, b_vec, c_vec = cell
        a = norm(a_vec)
        b = norm(b_vec)
        c = norm(c_vec)
        alpha = self._angle_degrees(b_vec, c_vec)
        beta = self._angle_degrees(a_vec, c_vec)
        gamma = self._angle_degrees(a_vec, b_vec)
        return a, b, c, alpha, beta, gamma

    def _angle_degrees(self, u: Vec3, v: Vec3) -> float:
        cosine = dot(u, v) / max(norm(u) * norm(v), 1e-12)
        cosine = max(-1.0, min(1.0, cosine))
        return degrees(acos(cosine))

    def _cartesian_to_fractional(self, cell: tuple[Vec3, Vec3, Vec3], point: Vec3) -> Vec3:
        a_vec, b_vec, c_vec = cell
        volume = dot(a_vec, cross(b_vec, c_vec))
        if abs(volume) < 1e-12:
            raise ValueError("cell vectors are singular; cannot convert to fractional coordinates")
        return (
            dot(point, cross(b_vec, c_vec)) / volume,
            dot(point, cross(c_vec, a_vec)) / volume,
            dot(point, cross(a_vec, b_vec)) / volume,
        )

    def _wrap_fractional(self, frac: Vec3) -> Vec3:
        return tuple(value - floor(value) for value in frac)  # type: ignore[return-value]

    def _format_p1_symmetry_shift(self, shift: tuple[int, int, int]) -> str:
        if shift == (0, 0, 0):
            return "."
        if any(component < -4 or component > 4 for component in shift):
            raise ValueError(f"periodic image shift {shift!r} exceeds CIF P1 formatter range")
        return "1_" + "".join(str(component + 5) for component in shift)

    def _minimum_image_bond_geometry(
        self,
        cell: tuple[Vec3, Vec3, Vec3],
        left_fractional: Vec3,
        right_fractional: Vec3,
    ) -> tuple[tuple[int, int, int], float]:
        best_shift = (0, 0, 0)
        best_distance = float("inf")
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    shift = (dx, dy, dz)
                    vector = (
                        (right_fractional[0] + dx) - left_fractional[0],
                        (right_fractional[1] + dy) - left_fractional[1],
                        (right_fractional[2] + dz) - left_fractional[2],
                    )
                    distance = norm(self._fractional_vector_to_cartesian(cell, vector))
                    if distance + 1e-9 < best_distance:
                        best_distance = distance
                        best_shift = shift
                    elif abs(distance - best_distance) <= 1e-9 and shift == (0, 0, 0):
                        best_shift = shift
        return best_shift, best_distance

    def _fractional_vector_to_cartesian(self, cell: tuple[Vec3, Vec3, Vec3], vector: Vec3) -> Vec3:
        a_vec, b_vec, c_vec = cell
        return add(add(scale(a_vec, vector[0]), scale(b_vec, vector[1])), scale(c_vec, vector[2]))

    def _world_position(self, pose: Pose, local_position: Vec3) -> Vec3:
        return add(pose.translation, matmul_vec(pose.rotation_matrix, local_position))

    def _distance(self, left: Vec3, right: Vec3) -> float:
        dx = left[0] - right[0]
        dy = left[1] - right[1]
        dz = left[2] - right[2]
        return (dx * dx + dy * dy + dz * dz) ** 0.5

    def _motif_symbol(self, kind: str) -> str:
        return motif_pseudo_atom_symbol(kind)

    def _instance_to_monomer(self, candidate: Candidate) -> dict[str, str]:
        raw = candidate.metadata.get("instance_to_monomer")
        if isinstance(raw, dict) and raw:
            return {str(k): str(v) for k, v in raw.items()}

        assignment = candidate.metadata.get("assignment")
        if not isinstance(assignment, dict) or not assignment:
            raise KeyError("candidate metadata does not contain instance_to_monomer or assignment mapping")

        monomer_ids = list(assignment.values())
        instance_ids = sorted(candidate.state.monomer_poses.keys(), key=self._instance_sort_key)
        if len(instance_ids) != len(monomer_ids):
            raise ValueError("candidate metadata assignment does not align with monomer poses")
        return {instance_id: str(monomer_id) for instance_id, monomer_id in zip(instance_ids, monomer_ids)}

    def _instance_sort_key(self, instance_id: str) -> tuple[int, str]:
        digits = "".join(ch for ch in instance_id if ch.isdigit())
        return (int(digits) if digits else 0, instance_id)

    def _coerce_spec_map(self, monomer_specs: Mapping[str, MonomerSpec] | Iterable[MonomerSpec]) -> dict[str, MonomerSpec]:
        if isinstance(monomer_specs, Mapping):
            return {str(key): value for key, value in monomer_specs.items()}
        return {spec.id: spec for spec in monomer_specs}

    def _sanitize_data_name(self, data_name: str) -> str:
        cleaned = "".join(ch if ch.isalnum() or ch in {"_", "-"} else "_" for ch in data_name.strip())
        return cleaned or "cofkit_candidate"


def candidate_to_cif(
    candidate: Candidate,
    monomer_specs: Mapping[str, MonomerSpec] | Iterable[MonomerSpec],
    data_name: str | None = None,
    cofid: str | None = None,
    cofid_comment_suffix: str | None = None,
) -> str:
    return CIFWriter().export_candidate(
        candidate,
        monomer_specs,
        data_name=data_name,
        cofid=cofid,
        cofid_comment_suffix=cofid_comment_suffix,
    ).text


def write_candidate_cif(
    path: str | Path,
    candidate: Candidate,
    monomer_specs: Mapping[str, MonomerSpec] | Iterable[MonomerSpec],
    data_name: str | None = None,
    cofid: str | None = None,
    cofid_comment_suffix: str | None = None,
) -> CIFExportResult:
    return CIFWriter().write_candidate(
        path,
        candidate,
        monomer_specs,
        data_name=data_name,
        cofid=cofid,
        cofid_comment_suffix=cofid_comment_suffix,
    )
