from __future__ import annotations

from dataclasses import dataclass, field
from math import atan2
from typing import Callable, Mapping

from ..geometry import Frame, centroid, cross, dot, normalize, sub
from ..model import MonomerSpec, ReactiveMotif
from .motif_registry import MotifKindDefinition, MotifKindRegistry, default_motif_kind_registry

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except ImportError:  # pragma: no cover - handled at call sites
    Chem = None
    AllChem = None


@dataclass(frozen=True)
class _DetectedMotif:
    reactive_atom_id: int
    anchor_atom_id: int
    atom_ids: tuple[int, ...]
    origin: tuple[float, float, float]
    anchor: tuple[float, float, float]
    metadata: Mapping[str, object] = field(default_factory=dict)


RDKitMatchHandler = Callable[[object, object, tuple[int, ...], MotifKindDefinition], _DetectedMotif | None]
RDKitPostprocessHandler = Callable[[object, tuple[_DetectedMotif, ...], MotifKindDefinition], tuple[_DetectedMotif, ...]]


class RDKitMotifBuilder:
    def __init__(
        self,
        *,
        motif_registry: MotifKindRegistry | None = None,
        match_handlers: Mapping[str, RDKitMatchHandler] | None = None,
        postprocess_handlers: Mapping[str, RDKitPostprocessHandler] | None = None,
    ) -> None:
        self.motif_registry = motif_registry or default_motif_kind_registry()
        self._match_handlers: dict[str, RDKitMatchHandler] = dict(match_handlers or {})
        self._postprocess_handlers: dict[str, RDKitPostprocessHandler] = dict(postprocess_handlers or {})

    def register_match_handler(self, motif_kind: str, handler: RDKitMatchHandler) -> None:
        self._match_handlers[motif_kind] = handler

    def register_postprocess_handler(self, motif_kind: str, handler: RDKitPostprocessHandler) -> None:
        self._postprocess_handlers[motif_kind] = handler

    def supported_motif_kinds(self) -> tuple[str, ...]:
        return tuple(sorted(self._match_handlers))

    def build_monomer(
        self,
        monomer_id: str,
        name: str,
        smiles: str,
        motif_kind: str,
        *,
        num_conformers: int = 8,
        random_seed: int = 0xC0F,
    ) -> MonomerSpec:
        if Chem is None or AllChem is None:
            raise RuntimeError("RDKit is required for build_rdkit_monomer()")
        definition = self.motif_registry.get(motif_kind)
        if definition.rdkit_smarts is None:
            raise ValueError(f"motif kind {motif_kind!r} has no RDKit SMARTS configuration")

        base = Chem.MolFromSmiles(smiles)
        if base is None:
            raise ValueError(f"RDKit could not parse SMILES for {monomer_id!r}")
        molecule = Chem.AddHs(base)
        conformer_ids = _embed_conformers(molecule, num_conformers=num_conformers, random_seed=random_seed)
        best_conf_id, best_energy, forcefield = _optimize_conformers(molecule, conformer_ids)
        conformer = molecule.GetConformer(best_conf_id)

        detected = self._detect_motifs(molecule, conformer, definition)
        if not detected:
            raise ValueError(f"no {motif_kind!r} motifs detected in {monomer_id!r}")

        atom_positions, motifs, plane_normal = _build_geometry(detected, molecule, conformer, definition)
        atom_symbols = tuple(atom.GetSymbol() for atom in molecule.GetAtoms())
        bonds = tuple(
            (
                int(bond.GetBeginAtomIdx()),
                int(bond.GetEndAtomIdx()),
                float(bond.GetBondTypeAsDouble()),
            )
            for bond in molecule.GetBonds()
        )
        return MonomerSpec(
            id=monomer_id,
            name=name,
            motifs=motifs,
            conformer_ids=(f"rdkit-conf-{best_conf_id}",),
            atom_symbols=atom_symbols,
            atom_positions=atom_positions,
            bonds=bonds,
            metadata={
                "source_smiles": smiles,
                "geometry_mode": "rdkit-etkdg",
                "motif_detection": f"rdkit-smarts:{definition.rdkit_smarts}",
                "motif_kind": motif_kind,
                "n_atoms": molecule.GetNumAtoms(),
                "n_heavy_atoms": base.GetNumAtoms(),
                "n_conformers": len(conformer_ids),
                "selected_conformer_id": best_conf_id,
                "forcefield": forcefield,
                "selected_conformer_energy": best_energy,
                "plane_normal": plane_normal,
            },
        )

    def _detect_motifs(self, molecule, conformer, definition: MotifKindDefinition) -> tuple[_DetectedMotif, ...]:
        assert definition.rdkit_smarts is not None
        pattern = Chem.MolFromSmarts(definition.rdkit_smarts)
        if pattern is None:
            raise ValueError(f"invalid SMARTS for motif kind {definition.kind!r}")
        try:
            handler = self._match_handlers[definition.kind]
        except KeyError as exc:
            raise ValueError(f"motif kind {definition.kind!r} has no RDKit match handler") from exc
        matches = molecule.GetSubstructMatches(pattern)
        detected_by_atoms: dict[tuple[int, ...], _DetectedMotif] = {}
        for match in matches:
            detected = handler(molecule, conformer, tuple(int(atom_id) for atom_id in match), definition)
            if detected is None:
                continue
            detected_by_atoms.setdefault(tuple(sorted(detected.atom_ids)), detected)
        return tuple(
            self._postprocess_detected_motifs(
                molecule,
                tuple(
                    sorted(
                        detected_by_atoms.values(),
                        key=lambda item: (item.reactive_atom_id, item.anchor_atom_id, item.atom_ids),
                    )
                ),
                definition,
            )
        )

    def _postprocess_detected_motifs(
        self,
        molecule,
        detected: tuple[_DetectedMotif, ...],
        definition: MotifKindDefinition,
    ) -> tuple[_DetectedMotif, ...]:
        handler = self._postprocess_handlers.get(definition.kind)
        if handler is None:
            return detected
        return handler(molecule, detected, definition)

    @classmethod
    def builtin(cls, *, motif_registry: MotifKindRegistry | None = None) -> "RDKitMotifBuilder":
        builder = cls(motif_registry=motif_registry)
        builder.register_match_handler("amine", _interpret_primary_amine_match)
        builder.register_match_handler("aldehyde", _interpret_aldehyde_match)
        builder.register_match_handler("hydrazide", _interpret_hydrazide_match)
        builder.register_match_handler("boronic_acid", _interpret_boronic_acid_match)
        builder.register_match_handler("catechol", _interpret_catechol_match)
        builder.register_match_handler("keto_aldehyde", _interpret_keto_aldehyde_match)
        builder.register_match_handler("activated_methylene", _interpret_activated_methylene_match)
        builder.register_postprocess_handler("keto_aldehyde", _postprocess_keto_aldehyde_matches)
        return builder


def build_rdkit_monomer(
    monomer_id: str,
    name: str,
    smiles: str,
    motif_kind: str,
    *,
    num_conformers: int = 8,
    random_seed: int = 0xC0F,
    motif_registry: MotifKindRegistry | None = None,
    builder: RDKitMotifBuilder | None = None,
) -> MonomerSpec:
    effective_builder = builder or RDKitMotifBuilder.builtin(motif_registry=motif_registry)
    return effective_builder.build_monomer(
        monomer_id,
        name,
        smiles,
        motif_kind,
        num_conformers=num_conformers,
        random_seed=random_seed,
    )


def _embed_conformers(molecule, *, num_conformers: int, random_seed: int) -> tuple[int, ...]:
    params = AllChem.ETKDGv3()
    params.randomSeed = random_seed
    params.pruneRmsThresh = 0.2
    params.useSmallRingTorsions = True
    conf_ids = AllChem.EmbedMultipleConfs(molecule, numConfs=max(1, num_conformers), params=params)
    if not conf_ids:
        raise ValueError("RDKit conformer embedding failed")
    return tuple(int(conf_id) for conf_id in conf_ids)


def _optimize_conformers(molecule, conformer_ids: tuple[int, ...]) -> tuple[int, float, str]:
    if AllChem.MMFFHasAllMoleculeParams(molecule):
        props = AllChem.MMFFGetMoleculeProperties(molecule)
        if props is not None:
            best = None
            for conf_id in conformer_ids:
                field = AllChem.MMFFGetMoleculeForceField(molecule, props, confId=conf_id)
                if field is None:
                    continue
                field.Minimize(maxIts=500)
                energy = float(field.CalcEnergy())
                if best is None or energy < best[1]:
                    best = (conf_id, energy, "MMFF")
            if best is not None:
                return best

    best = None
    for conf_id in conformer_ids:
        field = AllChem.UFFGetMoleculeForceField(molecule, confId=conf_id)
        if field is None:
            continue
        field.Minimize(maxIts=500)
        energy = float(field.CalcEnergy())
        if best is None or energy < best[1]:
            best = (conf_id, energy, "UFF")
    if best is None:
        raise ValueError("RDKit force-field optimization failed")
    return best


def _build_geometry(
    detected: tuple[_DetectedMotif, ...],
    molecule,
    conformer,
    definition: MotifKindDefinition,
) -> tuple[tuple[tuple[float, float, float], ...], tuple[ReactiveMotif, ...], tuple[float, float, float]]:
    points = tuple(_conformer_point(conformer, atom.GetIdx()) for atom in molecule.GetAtoms())
    center = centroid(points)
    plane_normal = _plane_normal(detected)
    provisional_primary = _project_onto_plane(sub(detected[0].origin, detected[0].anchor), plane_normal)
    if _is_near_zero(provisional_primary):
        provisional_primary = (1.0, 0.0, 0.0)
    x_axis = normalize(provisional_primary)
    y_axis = normalize(cross(plane_normal, x_axis))

    def transform(point: tuple[float, float, float]) -> tuple[float, float, float]:
        shifted = sub(point, center)
        return (
            dot(shifted, x_axis),
            dot(shifted, y_axis),
            dot(shifted, plane_normal),
        )

    atom_positions = tuple(transform(point) for point in points)
    motif_rows: list[tuple[float, ReactiveMotif]] = []
    for detected_motif in detected:
        origin = transform(detected_motif.origin)
        anchor = transform(detected_motif.anchor)
        primary = _project_onto_plane(sub(origin, anchor), (0.0, 0.0, 1.0))
        if _is_near_zero(primary):
            primary = (1.0, 0.0, 0.0)
        frame = Frame(origin=origin, primary=normalize(primary), normal=(0.0, 0.0, 1.0))
        angle = atan2(origin[1], origin[0])
        motif_rows.append(
            (
                angle,
                ReactiveMotif(
                    id=f"{definition.id_prefix}{len(motif_rows) + 1}",
                    kind=definition.kind,
                    atom_ids=detected_motif.atom_ids,
                    frame=frame,
                    allowed_reaction_templates=definition.allowed_reaction_templates,
                    metadata={
                        **detected_motif.metadata,
                        "reactive_atom_id": detected_motif.reactive_atom_id,
                        "anchor_atom_id": detected_motif.anchor_atom_id,
                    },
                ),
            )
        )

    motifs = tuple(
        ReactiveMotif(
            id=f"{definition.id_prefix}{index}",
            kind=motif.kind,
            atom_ids=motif.atom_ids,
            frame=motif.frame,
            valence=motif.valence,
            symmetry_order=motif.symmetry_order,
            planarity_class=motif.planarity_class,
            allowed_reaction_templates=motif.allowed_reaction_templates,
            metadata=motif.metadata,
        )
        for index, (_, motif) in enumerate(sorted(motif_rows, key=lambda item: item[0]), start=1)
    )
    return atom_positions, motifs, plane_normal


def _interpret_primary_amine_match(molecule, conformer, match: tuple[int, ...], definition: MotifKindDefinition) -> _DetectedMotif:
    del definition
    reactive_atom_id, anchor_atom_id = int(match[0]), int(match[1])
    atom_ids = [reactive_atom_id, anchor_atom_id]
    atom_ids.extend(
        neighbor.GetIdx()
        for neighbor in molecule.GetAtomWithIdx(reactive_atom_id).GetNeighbors()
        if neighbor.GetAtomicNum() == 1
    )
    atom_ids = sorted(set(atom_ids))
    return _DetectedMotif(
        reactive_atom_id=reactive_atom_id,
        anchor_atom_id=anchor_atom_id,
        atom_ids=tuple(atom_ids),
        origin=_conformer_point(conformer, reactive_atom_id),
        anchor=_conformer_point(conformer, anchor_atom_id),
    )


def _interpret_aldehyde_match(molecule, conformer, match: tuple[int, ...], definition: MotifKindDefinition) -> _DetectedMotif:
    del definition
    reactive_atom_id, oxygen_atom_id, anchor_atom_id = (int(match[0]), int(match[1]), int(match[2]))
    atom_ids = [reactive_atom_id, oxygen_atom_id, anchor_atom_id]
    atom_ids.extend(
        neighbor.GetIdx()
        for neighbor in molecule.GetAtomWithIdx(reactive_atom_id).GetNeighbors()
        if neighbor.GetAtomicNum() == 1
    )
    atom_ids = sorted(set(atom_ids))
    return _DetectedMotif(
        reactive_atom_id=reactive_atom_id,
        anchor_atom_id=anchor_atom_id,
        atom_ids=tuple(atom_ids),
        origin=_conformer_point(conformer, reactive_atom_id),
        anchor=_conformer_point(conformer, anchor_atom_id),
    )


def _interpret_hydrazide_match(molecule, conformer, match: tuple[int, ...], definition: MotifKindDefinition) -> _DetectedMotif:
    del definition
    terminal_nitrogen_atom_id, internal_nitrogen_atom_id, carbonyl_carbon_atom_id, carbonyl_oxygen_atom_id, anchor_atom_id = (
        int(match[0]),
        int(match[1]),
        int(match[2]),
        int(match[3]),
        int(match[4]),
    )
    hydrogen_atom_ids = _attached_hydrogen_ids(molecule, terminal_nitrogen_atom_id)
    atom_ids = sorted(
        {
            terminal_nitrogen_atom_id,
            internal_nitrogen_atom_id,
            carbonyl_carbon_atom_id,
            carbonyl_oxygen_atom_id,
            anchor_atom_id,
            *hydrogen_atom_ids,
        }
    )
    return _DetectedMotif(
        reactive_atom_id=terminal_nitrogen_atom_id,
        anchor_atom_id=internal_nitrogen_atom_id,
        atom_ids=tuple(atom_ids),
        origin=_conformer_point(conformer, terminal_nitrogen_atom_id),
        anchor=_conformer_point(conformer, internal_nitrogen_atom_id),
        metadata={
            "internal_nitrogen_atom_id": internal_nitrogen_atom_id,
            "carbonyl_carbon_atom_id": carbonyl_carbon_atom_id,
            "carbonyl_oxygen_atom_id": carbonyl_oxygen_atom_id,
            "hydrogen_atom_ids": hydrogen_atom_ids,
        },
    )


def _interpret_boronic_acid_match(molecule, conformer, match: tuple[int, ...], definition: MotifKindDefinition) -> _DetectedMotif:
    del definition
    anchor_atom_id, boron_atom_id, oxygen_atom_id_1, oxygen_atom_id_2 = (
        int(match[0]),
        int(match[1]),
        int(match[2]),
        int(match[3]),
    )
    hydrogen_atom_ids = tuple(
        hydrogen_atom_id
        for oxygen_atom_id in (oxygen_atom_id_1, oxygen_atom_id_2)
        for hydrogen_atom_id in _attached_hydrogen_ids(molecule, oxygen_atom_id)
    )
    atom_ids = sorted(
        {
            anchor_atom_id,
            boron_atom_id,
            oxygen_atom_id_1,
            oxygen_atom_id_2,
            *hydrogen_atom_ids,
        }
    )
    return _DetectedMotif(
        reactive_atom_id=boron_atom_id,
        anchor_atom_id=anchor_atom_id,
        atom_ids=tuple(atom_ids),
        origin=_conformer_point(conformer, boron_atom_id),
        anchor=_conformer_point(conformer, anchor_atom_id),
        metadata={
            "oxygen_atom_ids": (oxygen_atom_id_1, oxygen_atom_id_2),
            "hydrogen_atom_ids": hydrogen_atom_ids,
        },
    )


def _interpret_catechol_match(molecule, conformer, match: tuple[int, ...], definition: MotifKindDefinition) -> _DetectedMotif | None:
    del definition
    oxygen_atom_id_1, carbon_atom_id_1 = (int(match[0]), int(match[1]))
    neighbor_pairs = []
    for aromatic_neighbor in molecule.GetAtomWithIdx(carbon_atom_id_1).GetNeighbors():
        carbon_atom_id_2 = int(aromatic_neighbor.GetIdx())
        if carbon_atom_id_2 == oxygen_atom_id_1:
            continue
        if aromatic_neighbor.GetAtomicNum() != 6 or not aromatic_neighbor.GetIsAromatic():
            continue
        oxygen_atom_id_2 = _aromatic_hydroxyl_substituent(molecule, carbon_atom_id_2, excluded_atom_id=carbon_atom_id_1)
        if oxygen_atom_id_2 is None:
            continue
        neighbor_pairs.append((carbon_atom_id_2, oxygen_atom_id_2))
    if not neighbor_pairs:
        return None
    carbon_atom_id_2, oxygen_atom_id_2 = min(neighbor_pairs)
    if carbon_atom_id_1 > carbon_atom_id_2:
        return None
    hydrogen_atom_ids_1 = _attached_hydrogen_ids(molecule, oxygen_atom_id_1)
    hydrogen_atom_ids_2 = _attached_hydrogen_ids(molecule, oxygen_atom_id_2)
    if len(hydrogen_atom_ids_1) != 1 or len(hydrogen_atom_ids_2) != 1:
        return None
    hydrogen_atom_ids = hydrogen_atom_ids_1 + hydrogen_atom_ids_2
    origin = centroid(
        (
            _conformer_point(conformer, oxygen_atom_id_1),
            _conformer_point(conformer, oxygen_atom_id_2),
        )
    )
    anchor = centroid(
        (
            _conformer_point(conformer, carbon_atom_id_1),
            _conformer_point(conformer, carbon_atom_id_2),
        )
    )
    atom_ids = sorted(
        {
            oxygen_atom_id_1,
            oxygen_atom_id_2,
            carbon_atom_id_1,
            carbon_atom_id_2,
            *hydrogen_atom_ids,
        }
    )
    return _DetectedMotif(
        reactive_atom_id=min(oxygen_atom_id_1, oxygen_atom_id_2),
        anchor_atom_id=min(carbon_atom_id_1, carbon_atom_id_2),
        atom_ids=tuple(atom_ids),
        origin=origin,
        anchor=anchor,
        metadata={
            "reactive_atom_ids": (oxygen_atom_id_1, oxygen_atom_id_2),
            "anchor_atom_ids": (carbon_atom_id_1, carbon_atom_id_2),
            "hydrogen_atom_ids": hydrogen_atom_ids,
        },
    )


def _interpret_keto_aldehyde_match(molecule, conformer, match: tuple[int, ...], definition: MotifKindDefinition) -> _DetectedMotif | None:
    del definition
    reactive_atom_id, oxygen_atom_id, anchor_atom_id = (int(match[0]), int(match[1]), int(match[2]))
    ortho_hydroxyl = _find_ortho_hydroxyl(molecule, anchor_atom_id)
    if ortho_hydroxyl is None:
        return None
    hydroxyl_oxygen_atom_id, hydroxyl_hydrogen_atom_id, hydroxyl_anchor_atom_id = ortho_hydroxyl
    atom_ids = [reactive_atom_id, oxygen_atom_id, anchor_atom_id, hydroxyl_oxygen_atom_id, hydroxyl_hydrogen_atom_id]
    atom_ids.extend(_attached_hydrogen_ids(molecule, reactive_atom_id))
    atom_ids = sorted(set(atom_ids))
    return _DetectedMotif(
        reactive_atom_id=reactive_atom_id,
        anchor_atom_id=anchor_atom_id,
        atom_ids=tuple(atom_ids),
        origin=_conformer_point(conformer, reactive_atom_id),
        anchor=_conformer_point(conformer, anchor_atom_id),
        metadata={
            "aldehyde_oxygen_atom_id": oxygen_atom_id,
            "ortho_hydroxyl_oxygen_atom_id": hydroxyl_oxygen_atom_id,
            "ortho_hydroxyl_hydrogen_atom_id": hydroxyl_hydrogen_atom_id,
            "ortho_hydroxyl_anchor_atom_id": hydroxyl_anchor_atom_id,
        },
    )


def _postprocess_keto_aldehyde_matches(
    molecule,
    detected: tuple[_DetectedMotif, ...],
    definition: MotifKindDefinition,
) -> tuple[_DetectedMotif, ...]:
    del definition
    return _assign_unique_keto_aldehyde_hydroxyls(molecule, detected)


def _interpret_activated_methylene_match(molecule, conformer, match: tuple[int, ...], definition: MotifKindDefinition) -> _DetectedMotif | None:
    del definition
    reactive_atom_id = int(match[0])
    hydrogen_atom_ids = _attached_hydrogen_ids(molecule, reactive_atom_id)
    if len(hydrogen_atom_ids) < 2:
        return None
    activating_neighbor_ids = tuple(
        neighbor.GetIdx()
        for neighbor in molecule.GetAtomWithIdx(reactive_atom_id).GetNeighbors()
        if neighbor.GetAtomicNum() != 1 and _is_electron_withdrawing_anchor(neighbor, reactive_atom_id)
    )
    if not activating_neighbor_ids:
        return None
    anchor_atom_id = min(activating_neighbor_ids)
    anchor = centroid(tuple(_conformer_point(conformer, atom_id) for atom_id in activating_neighbor_ids))
    atom_ids = sorted({reactive_atom_id, *hydrogen_atom_ids, *activating_neighbor_ids})
    return _DetectedMotif(
        reactive_atom_id=reactive_atom_id,
        anchor_atom_id=anchor_atom_id,
        atom_ids=tuple(atom_ids),
        origin=_conformer_point(conformer, reactive_atom_id),
        anchor=anchor,
        metadata={
            "hydrogen_atom_ids": hydrogen_atom_ids,
            "activator_atom_ids": activating_neighbor_ids,
        },
    )


def _plane_normal(detected: tuple[_DetectedMotif, ...]) -> tuple[float, float, float]:
    anchor_points = [item.anchor for item in detected]
    if len(anchor_points) >= 3:
        first, second, third = anchor_points[:3]
        normal = cross(sub(second, first), sub(third, first))
        if not _is_near_zero(normal):
            normalized = normalize(normal)
            if normalized[2] < 0.0:
                return (-normalized[0], -normalized[1], -normalized[2])
            return normalized
    return (0.0, 0.0, 1.0)


def _project_onto_plane(
    vector: tuple[float, float, float],
    plane_normal: tuple[float, float, float],
) -> tuple[float, float, float]:
    scale = dot(vector, plane_normal)
    return (
        vector[0] - scale * plane_normal[0],
        vector[1] - scale * plane_normal[1],
        vector[2] - scale * plane_normal[2],
    )


def _is_near_zero(vector: tuple[float, float, float]) -> bool:
    return abs(vector[0]) + abs(vector[1]) + abs(vector[2]) < 1e-8


def _conformer_point(conformer, atom_id: int) -> tuple[float, float, float]:
    position = conformer.GetAtomPosition(atom_id)
    return (float(position.x), float(position.y), float(position.z))


def _attached_hydrogen_ids(molecule, atom_id: int) -> tuple[int, ...]:
    return tuple(
        int(neighbor.GetIdx())
        for neighbor in molecule.GetAtomWithIdx(atom_id).GetNeighbors()
        if neighbor.GetAtomicNum() == 1
    )


def _find_ortho_hydroxyl(molecule, anchor_atom_id: int) -> tuple[int, int, int] | None:
    options = _find_ortho_hydroxyl_options(molecule, anchor_atom_id)
    return options[0] if options else None


def _find_ortho_hydroxyl_options(molecule, anchor_atom_id: int) -> tuple[tuple[int, int, int], ...]:
    anchor_atom = molecule.GetAtomWithIdx(anchor_atom_id)
    options: list[tuple[int, int, int]] = []
    for aromatic_neighbor in anchor_atom.GetNeighbors():
        if aromatic_neighbor.GetAtomicNum() != 6 or not aromatic_neighbor.GetIsAromatic():
            continue
        oxygen_atom_id = _aromatic_hydroxyl_substituent(molecule, aromatic_neighbor.GetIdx(), excluded_atom_id=anchor_atom_id)
        if oxygen_atom_id is None:
            continue
        hydrogen_atom_ids = _attached_hydrogen_ids(molecule, oxygen_atom_id)
        if len(hydrogen_atom_ids) == 1:
            options.append((oxygen_atom_id, hydrogen_atom_ids[0], aromatic_neighbor.GetIdx()))
    return tuple(sorted(set(options)))


def _assign_unique_keto_aldehyde_hydroxyls(molecule, detected: tuple[_DetectedMotif, ...]) -> tuple[_DetectedMotif, ...]:
    if len(detected) <= 1:
        return detected

    options_by_index = [
        _find_ortho_hydroxyl_options(molecule, motif.anchor_atom_id)
        for motif in detected
    ]
    if any(not options for options in options_by_index):
        return detected

    ordered_indices = tuple(sorted(range(len(detected)), key=lambda idx: (len(options_by_index[idx]), detected[idx].reactive_atom_id)))
    best_assignment: dict[int, tuple[int, int, int]] | None = None
    best_score: tuple[int, tuple[int, ...]] | None = None

    def _search(position: int, used_oxygen_ids: set[int], assignment: dict[int, tuple[int, int, int]]) -> None:
        nonlocal best_assignment, best_score
        if position == len(ordered_indices):
            chosen_oxygen_ids = tuple(assignment[idx][0] for idx in range(len(detected)))
            score = (len(set(chosen_oxygen_ids)), tuple(chosen_oxygen_ids))
            if best_score is None or score > best_score:
                best_score = score
                best_assignment = dict(assignment)
            return

        detected_index = ordered_indices[position]
        options = sorted(options_by_index[detected_index], key=lambda option: (option[0] in used_oxygen_ids, option[0], option[2]))
        for option in options:
            assignment[detected_index] = option
            next_used = set(used_oxygen_ids)
            next_used.add(option[0])
            _search(position + 1, next_used, assignment)
            del assignment[detected_index]

    _search(0, set(), {})
    if best_assignment is None:
        return detected

    resolved: list[_DetectedMotif] = []
    for index, motif in enumerate(detected):
        oxygen_atom_id, hydrogen_atom_id, hydroxyl_anchor_atom_id = best_assignment[index]
        original_oxygen_atom_id = motif.metadata.get("ortho_hydroxyl_oxygen_atom_id")
        original_hydrogen_atom_id = motif.metadata.get("ortho_hydroxyl_hydrogen_atom_id")
        atom_ids = set(motif.atom_ids)
        if isinstance(original_oxygen_atom_id, int):
            atom_ids.discard(original_oxygen_atom_id)
        if isinstance(original_hydrogen_atom_id, int):
            atom_ids.discard(original_hydrogen_atom_id)
        atom_ids.update((oxygen_atom_id, hydrogen_atom_id, hydroxyl_anchor_atom_id))
        metadata = dict(motif.metadata)
        metadata["ortho_hydroxyl_oxygen_atom_id"] = oxygen_atom_id
        metadata["ortho_hydroxyl_hydrogen_atom_id"] = hydrogen_atom_id
        metadata["ortho_hydroxyl_anchor_atom_id"] = hydroxyl_anchor_atom_id
        resolved.append(
            _DetectedMotif(
                reactive_atom_id=motif.reactive_atom_id,
                anchor_atom_id=motif.anchor_atom_id,
                atom_ids=tuple(sorted(atom_ids)),
                origin=motif.origin,
                anchor=motif.anchor,
                metadata=metadata,
            )
        )
    return tuple(resolved)


def _is_electron_withdrawing_anchor(atom, excluded_neighbor_id: int) -> bool:
    if atom.GetAtomicNum() != 6:
        return False
    if atom.GetIsAromatic():
        aromatic_nitrogen_neighbors = sum(
            1
            for neighbor in atom.GetNeighbors()
            if neighbor.GetIdx() != excluded_neighbor_id and neighbor.GetAtomicNum() == 7 and neighbor.GetIsAromatic()
        )
        if aromatic_nitrogen_neighbors >= 2:
            return True
    for bond in atom.GetBonds():
        other = bond.GetOtherAtom(atom)
        if other.GetIdx() == excluded_neighbor_id:
            continue
        bond_order = float(bond.GetBondTypeAsDouble())
        if bond_order >= 2.0 and other.GetAtomicNum() in (7, 8, 16):
            return True
        if bond_order >= 3.0 and other.GetAtomicNum() == 7:
            return True
    return False


def _aromatic_hydroxyl_substituent(molecule, carbon_atom_id: int, *, excluded_atom_id: int) -> int | None:
    for substituent in molecule.GetAtomWithIdx(carbon_atom_id).GetNeighbors():
        if substituent.GetIdx() == excluded_atom_id or substituent.GetAtomicNum() != 8:
            continue
        if len(_attached_hydrogen_ids(molecule, substituent.GetIdx())) == 1:
            return int(substituent.GetIdx())
    return None
