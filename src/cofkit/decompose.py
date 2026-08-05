"""CIF-to-COFid decomposition helpers.

The decomposition approach here is adapted from the
deCOFpose project: https://github.com/r-fedorov/deCOFpose
"""

from __future__ import annotations

import re
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from math import sqrt
from pathlib import Path
from typing import Callable, Mapping

from .bond_types import cif_type_to_bond_order, is_aromatic_bond_order
from .cofid import COFidMonomer, canonicalize_smiles, parse_cofid, read_cofid_from_cif, serialize_cofid
from .decompose_cif import PeriodicCifAtoms, read_periodic_cif_atoms
from .topologies import default_topology_repository
from .topology_symmetry import expand_topology_definition

try:
    from rdkit import Chem
except ImportError:  # pragma: no cover - project dependency guard
    Chem = None


@dataclass(frozen=True)
class DecomposedMonomer:
    connectivity: int
    reactive_group: str
    canonical_smiles: str
    amount: int = 1

    def to_cofid_monomer(self) -> COFidMonomer:
        return COFidMonomer(
            connectivity=self.connectivity,
            reactive_group=self.reactive_group,
            canonical_smiles=self.canonical_smiles,
        )


@dataclass(frozen=True)
class CifDecompositionResult:
    status: str
    input_cif: str
    topology: str | None
    linkage: str
    cofid: str | None = None
    monomers: tuple[DecomposedMonomer, ...] = ()
    reason: str | None = None
    metadata: Mapping[str, object] = field(default_factory=dict)

    @property
    def ok(self) -> bool:
        return self.status == "success" and self.cofid is not None

    def to_dict(self) -> dict[str, object]:
        return {
            "status": self.status,
            "input_cif": self.input_cif,
            "topology": self.topology,
            "linkage": self.linkage,
            "cofid": self.cofid,
            "monomers": [
                {
                    "connectivity": monomer.connectivity,
                    "reactive_group": monomer.reactive_group,
                    "canonical_smiles": monomer.canonical_smiles,
                    "amount": monomer.amount,
                }
                for monomer in self.monomers
            ],
            "reason": self.reason,
            "metadata": dict(self.metadata),
        }


FragmentRepairer = Callable[[object], object]
LinkageMarker = Callable[[object], tuple[int, ...]]
ConnectivityCounter = Callable[[object, str], int]
Vec3 = tuple[float, float, float]
_CIF_NUMBER_RE = re.compile(r"^([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?)(?:\(\d+\))?$")


@dataclass(frozen=True)
class LinkageDecompositionSpec:
    linkage_code: str
    template_id: str
    roles: tuple[str, ...]
    marker: LinkageMarker
    repairers: Mapping[str, FragmentRepairer] = field(default_factory=dict)
    connectivity_counter: ConnectivityCounter | None = None

    @property
    def metadata_key(self) -> str:
        return f"n_{self.linkage_code}_linkage_bonds"


@dataclass(frozen=True)
class RingDecompositionSpec:
    linkage_code: str
    template_id: str
    reactive_group: str
    anchor_atomic_num: int
    hetero_atomic_num: int
    repairer: FragmentRepairer

    @property
    def metadata_key(self) -> str:
        return f"n_{self.linkage_code}_rings"


DecompositionSpec = LinkageDecompositionSpec | RingDecompositionSpec


@dataclass(frozen=True)
class BondCandidate:
    atom_idx_1: int
    atom_idx_2: int
    distance: float
    explicit_order: float | None = None
    periodic_image: tuple[int, int, int] = (0, 0, 0)


@dataclass(frozen=True)
class BondedMolBuildResult:
    mol: object
    metadata: Mapping[str, object] = field(default_factory=dict)
    candidates: tuple[BondCandidate, ...] = ()


@dataclass(frozen=True)
class TopologyDetectionCandidate:
    topology: str
    score: float
    confidence: str
    reason: str
    metadata: Mapping[str, object] = field(default_factory=dict)

    def to_dict(self) -> dict[str, object]:
        return {
            "topology": self.topology,
            "score": self.score,
            "confidence": self.confidence,
            "reason": self.reason,
            "metadata": dict(self.metadata),
        }


@dataclass(frozen=True)
class TopologyDetectionResult:
    status: str
    selected_topology: str | None = None
    confidence: str = "failed"
    reason: str | None = None
    candidates: tuple[TopologyDetectionCandidate, ...] = ()
    metadata: Mapping[str, object] = field(default_factory=dict)

    @property
    def ok(self) -> bool:
        return self.status == "success" and self.selected_topology is not None

    def to_dict(self) -> dict[str, object]:
        return {
            "status": self.status,
            "selected_topology": self.selected_topology,
            "confidence": self.confidence,
            "reason": self.reason,
            "candidates": [candidate.to_dict() for candidate in self.candidates],
            "metadata": dict(self.metadata),
        }


@dataclass(frozen=True)
class LinkageTopologyGraph:
    node_connectivities: tuple[int, ...]
    gain_edges: tuple[tuple[int, int, tuple[int, int, int]], ...]
    dimensionality_hint: str | None

    @property
    def node_count(self) -> int:
        return len(self.node_connectivities)

    @property
    def edge_count(self) -> int:
        return len(self.gain_edges)

    def to_metadata(self) -> dict[str, object]:
        return {
            "node_connectivities": list(self.node_connectivities),
            "n_nodes": self.node_count,
            "n_edges": self.edge_count,
            "dimensionality_hint": self.dimensionality_hint,
        }


@dataclass(frozen=True)
class LinkageDecompositionDetails:
    monomers: tuple[DecomposedMonomer, ...]
    metadata: dict[str, object]
    topology_graph: LinkageTopologyGraph | None = None


@dataclass(frozen=True)
class NitrogenLinkageBondClassification:
    raw_cn_bond_indices: tuple[int, ...]
    azine_bond_indices: tuple[int, ...]
    hydrazone_bond_indices: tuple[int, ...]
    bken_bond_indices: tuple[int, ...]
    imine_bond_indices: tuple[int, ...]
    cross_branch_conflict_bond_indices: tuple[int, ...] = ()

    def bond_indices_for(self, linkage: str) -> tuple[int, ...]:
        return {
            "azine": self.azine_bond_indices,
            "hydrazone": self.hydrazone_bond_indices,
            "bken": self.bken_bond_indices,
            "imine": self.imine_bond_indices,
        }.get(str(linkage), ())

    def to_metadata(self) -> dict[str, object]:
        assigned_counts = {
            "azine": len(self.azine_bond_indices),
            "hydrazone": len(self.hydrazone_bond_indices),
            "bken": len(self.bken_bond_indices),
            "imine": len(self.imine_bond_indices),
        }
        detected_families = [
            linkage
            for linkage, count in assigned_counts.items()
            if count > 0
        ]
        return {
            "branches": {
                "nitrogen_nitrogen": ["azine", "hydrazone", "imine"],
                "keto_enamine": ["bken", "imine"],
            },
            "raw_cn_bond_count": len(self.raw_cn_bond_indices),
            "assigned_bond_counts": assigned_counts,
            "detected_families": detected_families,
            "cross_branch_conflict_bond_count": len(self.cross_branch_conflict_bond_indices),
            "resolution_status": (
                "cross_branch_ambiguous"
                if self.cross_branch_conflict_bond_indices
                else "resolved"
            ),
        }


@dataclass(frozen=True)
class VinyleneLinkageBondClassification:
    formal_cc_double_bond_indices: tuple[int, ...]
    same_instance_rejected_bond_indices: tuple[int, ...]
    small_ring_rejected_bond_indices: tuple[int, ...]
    endpoint_rejected_bond_indices: tuple[int, ...]
    accepted_assignments: tuple[tuple[int, int, int], ...]

    @property
    def accepted_bond_indices(self) -> tuple[int, ...]:
        return tuple(assignment[0] for assignment in self.accepted_assignments)

    def to_metadata(self) -> dict[str, object]:
        return {
            "formal_cc_double_bond_count": len(self.formal_cc_double_bond_indices),
            "same_instance_rejected_bond_count": len(self.same_instance_rejected_bond_indices),
            "small_ring_rejected_bond_count": len(self.small_ring_rejected_bond_indices),
            "endpoint_rejected_bond_count": len(self.endpoint_rejected_bond_indices),
            "accepted_bond_count": len(self.accepted_assignments),
            "excluded_ring_sizes": [5, 6],
            "endpoint_requirement": (
                "acyclic/exocyclic C=C with carbon anchors on both endpoints and exactly one "
                "higher-scoring activated-methylene endpoint"
            ),
        }


@dataclass(frozen=True)
class RingDecompositionEvent:
    anchor_atom_indices: tuple[int, int, int]
    anchor_images: tuple[tuple[int, int, int], tuple[int, int, int], tuple[int, int, int]]


def decompose_cif_to_cofid(
    cif_path: str | Path,
    *,
    topology: str | None = None,
    linkage: str = "imine",
    bond_mode: str = "auto",
) -> CifDecompositionResult:
    input_path = Path(cif_path)
    spec = _resolve_linkage_spec(linkage)
    if spec is None:
        return CifDecompositionResult(
            status="skipped",
            input_cif=str(input_path),
            topology=topology,
            linkage=linkage,
            reason=(
                f"unsupported linkage {linkage!r}; current CIF decomposition supports "
                f"{tuple(sorted(_DECOMPOSITION_LINKAGE_ALIASES))!r}"
            ),
        )

    try:
        if bond_mode not in {"auto", "distance"}:
            raise ValueError("bond_mode must be 'auto' or 'distance'")
        atoms = read_periodic_cif_atoms(input_path)
        details = _decompose_atoms(atoms, spec, bond_mode=bond_mode)
        monomers = details.monomers
        metadata = dict(details.metadata)
        if not monomers:
            reason = f"no {spec.linkage_code} monomers were recovered"
            nitrogen_detection = metadata.get("nitrogen_linkage_detection")
            if (
                isinstance(nitrogen_detection, Mapping)
                and int(nitrogen_detection.get("cross_branch_conflict_bond_count", 0)) > 0
            ):
                reason = (
                    "nitrogen linkage classification is ambiguous across the "
                    "N-N and beta-ketoenamine priority branches"
                )
            return CifDecompositionResult(
                status="skipped",
                input_cif=str(input_path),
                topology=topology,
                linkage=spec.linkage_code,
                reason=reason,
                metadata=metadata,
            )
        selected_topology = topology
        if selected_topology is None:
            detection = detect_cif_topology(
                input_path,
                linkage=spec.linkage_code,
                bond_mode=bond_mode,
                atoms=atoms,
                details=details,
            )
            metadata["topology_detection"] = detection.to_dict()
            if not detection.ok:
                return CifDecompositionResult(
                    status="skipped",
                    input_cif=str(input_path),
                    topology=None,
                    linkage=spec.linkage_code,
                    reason=detection.reason or "topology auto-detection failed",
                    monomers=monomers,
                    metadata=metadata,
                )
            selected_topology = str(detection.selected_topology)
        cofid_monomers = tuple(
            sorted(
                (monomer.to_cofid_monomer() for monomer in monomers),
                key=lambda monomer: (-monomer.connectivity, monomer.canonical_smiles),
            )
        )
        cofid = serialize_cofid(monomers=cofid_monomers, topology=selected_topology, linkage=spec.linkage_code)
        return CifDecompositionResult(
            status="success",
            input_cif=str(input_path),
            topology=selected_topology,
            linkage=spec.linkage_code,
            cofid=cofid,
            monomers=monomers,
            metadata=metadata,
        )
    except Exception as exc:
        return CifDecompositionResult(
            status="skipped",
            input_cif=str(input_path),
            topology=topology,
            linkage=spec.linkage_code,
            reason=f"{type(exc).__name__}: {exc}",
        )


def _decompose_linkage_atoms(
    atoms: PeriodicCifAtoms,
    spec: LinkageDecompositionSpec,
    *,
    bond_mode: str,
) -> LinkageDecompositionDetails:
    if Chem is None:
        raise RuntimeError("RDKit is required for CIF decomposition.")
    build_result = _build_bonded_mol(atoms, bond_mode=bond_mode)
    mol = build_result.mol
    nitrogen_detection_metadata: dict[str, object] = {}
    if spec.linkage_code in {"azine", "hydrazone", "bken", "imine"}:
        nitrogen_detection_metadata = {
            "nitrogen_linkage_detection": _classify_nitrogen_linkage_bonds(mol).to_metadata()
        }
    vinylene_detection_metadata: dict[str, object] = {}
    if spec.linkage_code == "vinylene":
        vinylene_detection_metadata = {
            "vinylene_linkage_detection": _classify_vinylene_linkage_bonds(mol).to_metadata()
        }
    linkage_bond_indices = spec.marker(mol)
    if not linkage_bond_indices:
        return LinkageDecompositionDetails((), {
            **dict(build_result.metadata),
            **nitrogen_detection_metadata,
            **vinylene_detection_metadata,
            "n_atoms": mol.GetNumAtoms(),
            "n_bonds": mol.GetNumBonds(),
            spec.metadata_key: 0,
        })

    editable = Chem.RWMol(mol)
    linkage_bond_atom_pairs: list[tuple[int, int]] = []
    for bond_idx in sorted(linkage_bond_indices, reverse=True):
        bond = editable.GetBondWithIdx(bond_idx)
        linkage_bond_atom_pairs.append((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))
        editable.RemoveBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())

    cut_mol = editable.GetMol()
    fragment_members = Chem.GetMolFrags(cut_mol, asMols=False, sanitizeFrags=False)
    atom_to_fragment = {
        int(atom_idx): fragment_idx
        for fragment_idx, atom_indices in enumerate(fragment_members)
        for atom_idx in atom_indices
    }
    fragments = Chem.GetMolFrags(cut_mol, asMols=True, sanitizeFrags=False)
    monomers_by_key: dict[tuple[str, str, int], tuple[DecomposedMonomer, int]] = {}
    skipped_fragments = 0
    for fragment in fragments:
        monomer = _repair_linkage_fragment_to_monomer(fragment, spec)
        if monomer is None:
            skipped_fragments += 1
            continue
        key = (monomer.reactive_group, monomer.canonical_smiles, monomer.connectivity)
        previous = monomers_by_key.get(key)
        if previous is None:
            monomers_by_key[key] = (monomer, 1)
        else:
            monomers_by_key[key] = (previous[0], previous[1] + 1)

    monomers = tuple(
        sorted(
            (
                DecomposedMonomer(
                    connectivity=monomer.connectivity,
                    reactive_group=monomer.reactive_group,
                    canonical_smiles=monomer.canonical_smiles,
                    amount=amount,
                )
                for monomer, amount in monomers_by_key.values()
            ),
            key=lambda monomer: (-monomer.connectivity, monomer.canonical_smiles),
        )
    )
    topology_graph = _linkage_topology_graph(
        mol,
        linkage_bond_atom_pairs,
        atom_to_fragment,
        build_result.candidates,
    )
    metadata = {
        **dict(build_result.metadata),
        **nitrogen_detection_metadata,
        **vinylene_detection_metadata,
        "n_atoms": mol.GetNumAtoms(),
        "n_bonds": mol.GetNumBonds(),
        spec.metadata_key: len(linkage_bond_indices),
        "n_fragments_after_cut": len(fragments),
        "n_unique_monomers": len(monomers),
        "n_skipped_fragments": skipped_fragments,
    }
    if topology_graph is not None:
        metadata["topology_graph"] = topology_graph.to_metadata()
    return LinkageDecompositionDetails(monomers, metadata, topology_graph)


def _decompose_atoms(
    atoms: PeriodicCifAtoms,
    spec: DecompositionSpec,
    *,
    bond_mode: str,
) -> LinkageDecompositionDetails:
    if isinstance(spec, RingDecompositionSpec):
        return _decompose_ring_atoms(atoms, spec, bond_mode=bond_mode)
    return _decompose_linkage_atoms(atoms, spec, bond_mode=bond_mode)


def _decompose_ring_atoms(
    atoms: PeriodicCifAtoms,
    spec: RingDecompositionSpec,
    *,
    bond_mode: str,
) -> LinkageDecompositionDetails:
    if Chem is None:
        raise RuntimeError("RDKit is required for CIF decomposition.")
    build_result = _build_bonded_mol(atoms, bond_mode=bond_mode)
    mol = build_result.mol
    ring_events, ring_bond_indices = _mark_ring_decomposition_events(mol, build_result.candidates, spec)
    if not ring_events:
        return LinkageDecompositionDetails((), {
            **dict(build_result.metadata),
            "n_atoms": mol.GetNumAtoms(),
            "n_bonds": mol.GetNumBonds(),
            spec.metadata_key: 0,
            f"n_{spec.linkage_code}_ring_bonds": 0,
        })

    editable = Chem.RWMol(mol)
    for bond_idx in sorted(ring_bond_indices, reverse=True):
        bond = editable.GetBondWithIdx(bond_idx)
        editable.RemoveBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())

    cut_mol = editable.GetMol()
    fragment_members = Chem.GetMolFrags(cut_mol, asMols=False, sanitizeFrags=False)
    atom_to_fragment = {
        int(atom_idx): fragment_idx
        for fragment_idx, atom_indices in enumerate(fragment_members)
        for atom_idx in atom_indices
    }
    fragments = Chem.GetMolFrags(cut_mol, asMols=True, sanitizeFrags=False)
    monomers_by_key: dict[tuple[str, str, int], tuple[DecomposedMonomer, int]] = {}
    skipped_fragments = 0
    for fragment in fragments:
        monomer = _repair_ring_fragment_to_monomer(fragment, spec)
        if monomer is None:
            skipped_fragments += 1
            continue
        key = (monomer.reactive_group, monomer.canonical_smiles, monomer.connectivity)
        previous = monomers_by_key.get(key)
        if previous is None:
            monomers_by_key[key] = (monomer, 1)
        else:
            monomers_by_key[key] = (previous[0], previous[1] + 1)

    monomers = tuple(
        sorted(
            (
                DecomposedMonomer(
                    connectivity=monomer.connectivity,
                    reactive_group=monomer.reactive_group,
                    canonical_smiles=monomer.canonical_smiles,
                    amount=amount,
                )
                for monomer, amount in monomers_by_key.values()
            ),
            key=lambda monomer: (-monomer.connectivity, monomer.canonical_smiles),
        )
    )
    topology_graph, topology_component_count = _ring_topology_graph(ring_events, atom_to_fragment)
    metadata = {
        **dict(build_result.metadata),
        "decomposition_strategy": "ring",
        "n_atoms": mol.GetNumAtoms(),
        "n_bonds": mol.GetNumBonds(),
        spec.metadata_key: len(ring_events),
        f"n_{spec.linkage_code}_ring_bonds": len(ring_bond_indices),
        "n_fragments_after_cut": len(fragments),
        "n_unique_monomers": len(monomers),
        "n_skipped_fragments": skipped_fragments,
        "n_topology_components": topology_component_count,
    }
    if topology_graph is not None:
        metadata["topology_graph"] = topology_graph.to_metadata()
    return LinkageDecompositionDetails(monomers, metadata, topology_graph)


def _mark_ring_decomposition_events(
    mol,
    candidates: tuple[BondCandidate, ...],
    spec: RingDecompositionSpec,
) -> tuple[tuple[RingDecompositionEvent, ...], tuple[int, ...]]:
    Chem.GetSymmSSSR(mol)
    candidate_by_pair = {
        frozenset((candidate.atom_idx_1, candidate.atom_idx_2)): candidate
        for candidate in candidates
    }
    events: list[RingDecompositionEvent] = []
    marked_bonds: set[int] = set()
    used_atoms: set[int] = set()
    seen_rings: set[frozenset[int]] = set()
    for raw_ring in mol.GetRingInfo().AtomRings():
        ring = tuple(int(atom_idx) for atom_idx in raw_ring)
        ring_key = frozenset(ring)
        if ring_key in seen_rings or not _matches_ring_decomposition_pattern(mol, ring, spec):
            continue
        seen_rings.add(ring_key)
        if used_atoms.intersection(ring):
            continue

        ring_bonds: list[int] = []
        for index, atom_idx in enumerate(ring):
            next_idx = ring[(index + 1) % len(ring)]
            bond = mol.GetBondBetweenAtoms(atom_idx, next_idx)
            if bond is None:
                ring_bonds = []
                break
            ring_bonds.append(int(bond.GetIdx()))
        if len(ring_bonds) != 6:
            continue

        atom_images = _unwrap_ring_atom_images(ring, candidate_by_pair)
        if atom_images is None:
            continue
        anchors = tuple(atom_idx for atom_idx in ring if mol.GetAtomWithIdx(atom_idx).GetAtomicNum() == spec.anchor_atomic_num)
        if len(anchors) != 3:
            continue
        for atom_idx in anchors:
            mol.GetAtomWithIdx(atom_idx).SetProp("cofkit_decompose_role", spec.reactive_group)
        events.append(
            RingDecompositionEvent(
                anchor_atom_indices=anchors,  # type: ignore[arg-type]
                anchor_images=tuple(atom_images[atom_idx] for atom_idx in anchors),  # type: ignore[arg-type]
            )
        )
        marked_bonds.update(ring_bonds)
        used_atoms.update(ring)
    return tuple(events), tuple(sorted(marked_bonds))


def _matches_ring_decomposition_pattern(
    mol,
    ring: tuple[int, ...],
    spec: RingDecompositionSpec,
) -> bool:
    if len(ring) != 6 or len(set(ring)) != 6:
        return False
    ring_set = set(ring)
    expected = (spec.anchor_atomic_num, spec.hetero_atomic_num) * 3
    atomic_nums = tuple(mol.GetAtomWithIdx(atom_idx).GetAtomicNum() for atom_idx in ring)
    if atomic_nums not in {expected, tuple(reversed(expected))}:
        return False
    for atom_idx in ring:
        atom = mol.GetAtomWithIdx(atom_idx)
        ring_neighbors = sum(1 for neighbor in atom.GetNeighbors() if neighbor.GetIdx() in ring_set)
        external_neighbors = atom.GetDegree() - ring_neighbors
        if ring_neighbors != 2:
            return False
        if atom.GetAtomicNum() == spec.anchor_atomic_num and external_neighbors != 1:
            return False
        if atom.GetAtomicNum() == spec.hetero_atomic_num and external_neighbors != 0:
            return False
    return True


def _unwrap_ring_atom_images(
    ring: tuple[int, ...],
    candidate_by_pair: Mapping[frozenset[int], BondCandidate],
) -> dict[int, tuple[int, int, int]] | None:
    images = {ring[0]: (0, 0, 0)}
    for index in range(len(ring) - 1):
        current = ring[index]
        next_idx = ring[index + 1]
        candidate = candidate_by_pair.get(frozenset((current, next_idx)))
        step = (0, 0, 0) if candidate is None else _candidate_image_from_to(candidate, current, next_idx)
        current_image = images[current]
        images[next_idx] = tuple(current_image[axis] + step[axis] for axis in range(3))  # type: ignore[assignment]
    closing_candidate = candidate_by_pair.get(frozenset((ring[-1], ring[0])))
    closing_step = (
        (0, 0, 0)
        if closing_candidate is None
        else _candidate_image_from_to(closing_candidate, ring[-1], ring[0])
    )
    closing_image = tuple(images[ring[-1]][axis] + closing_step[axis] for axis in range(3))
    if closing_image != (0, 0, 0):
        return None
    return images


def _candidate_image_from_to(
    candidate: BondCandidate,
    start_atom_idx: int,
    end_atom_idx: int,
) -> tuple[int, int, int]:
    if (candidate.atom_idx_1, candidate.atom_idx_2) == (start_atom_idx, end_atom_idx):
        return candidate.periodic_image
    return _negate_image(candidate.periodic_image)


def _repair_ring_fragment_to_monomer(
    fragment,
    spec: RingDecompositionSpec,
) -> DecomposedMonomer | None:
    endpoint_roles = {
        atom.GetProp("cofkit_decompose_role")
        for atom in fragment.GetAtoms()
        if atom.HasProp("cofkit_decompose_role")
    }
    if endpoint_roles != {spec.reactive_group}:
        return None
    return _finalize_repaired_fragment(
        spec.repairer(fragment),
        spec.reactive_group,
        connectivity_counter=_connectivity_count,
    )


def _ring_topology_graph(
    events: tuple[RingDecompositionEvent, ...],
    atom_to_fragment: Mapping[int, int],
) -> tuple[LinkageTopologyGraph | None, int]:
    incidences_by_fragment: dict[int, list[tuple[int, tuple[int, int, int]]]] = defaultdict(list)
    for event_index, event in enumerate(events):
        for atom_idx, image in zip(event.anchor_atom_indices, event.anchor_images):
            fragment_idx = atom_to_fragment.get(atom_idx)
            if fragment_idx is not None:
                incidences_by_fragment[fragment_idx].append((event_index, image))
    if not incidences_by_fragment:
        return None, 0

    retained_fragments = tuple(
        sorted(fragment_idx for fragment_idx, incidences in incidences_by_fragment.items() if len(incidences) != 2)
    )
    fragment_node_by_id = {
        fragment_idx: len(events) + offset
        for offset, fragment_idx in enumerate(retained_fragments)
    }
    node_connectivities = [3] * len(events)
    node_connectivities.extend(len(incidences_by_fragment[fragment_idx]) for fragment_idx in retained_fragments)
    edges: list[tuple[int, int, tuple[int, int, int]]] = []

    for fragment_idx, incidences in sorted(incidences_by_fragment.items()):
        if len(incidences) == 2:
            (first_event, first_image), (second_event, second_image) = incidences
            image = tuple(first_image[axis] - second_image[axis] for axis in range(3))
            start, end = first_event, second_event
            if end < start:
                start, end = end, start
                image = _negate_image(image)  # type: ignore[arg-type]
            edges.append((start, end, image))  # type: ignore[arg-type]
            continue
        fragment_node = fragment_node_by_id[fragment_idx]
        for event_index, image in incidences:
            edges.append((event_index, fragment_node, image))

    if not edges:
        return None, 0
    graph = LinkageTopologyGraph(
        node_connectivities=tuple(int(value) for value in node_connectivities),
        gain_edges=tuple(sorted(edges)),
        dimensionality_hint="2D" if all(image[2] == 0 for _, _, image in edges) else "3D",
    )
    return _collapse_equivalent_topology_components(graph)


def _collapse_equivalent_topology_components(
    graph: LinkageTopologyGraph,
) -> tuple[LinkageTopologyGraph, int]:
    adjacency: dict[int, set[int]] = defaultdict(set)
    for start, end, _image in graph.gain_edges:
        adjacency[start].add(end)
        adjacency[end].add(start)
    remaining = set(range(graph.node_count))
    components: list[set[int]] = []
    while remaining:
        seed = min(remaining)
        component: set[int] = set()
        stack = [seed]
        while stack:
            node = stack.pop()
            if node in component:
                continue
            component.add(node)
            stack.extend(adjacency.get(node, ()))
        remaining.difference_update(component)
        components.append(component)
    if len(components) <= 1:
        return graph, len(components)

    component_graphs = tuple(_topology_component_graph(graph, component) for component in components)
    reference = component_graphs[0]
    if all(_topology_components_are_equivalent(reference, candidate) for candidate in component_graphs[1:]):
        return reference, len(component_graphs)
    return graph, len(component_graphs)


def _topology_components_are_equivalent(
    reference: LinkageTopologyGraph,
    candidate: LinkageTopologyGraph,
) -> bool:
    if (
        reference.node_connectivities == candidate.node_connectivities
        and _edge_counter(reference.gain_edges, compare_gains=False)
        == _edge_counter(candidate.gain_edges, compare_gains=False)
    ):
        return True
    return _periodic_graphs_match(reference, candidate, compare_gains=False)


def _topology_component_graph(
    graph: LinkageTopologyGraph,
    component: set[int],
) -> LinkageTopologyGraph:
    node_map = {old: new for new, old in enumerate(sorted(component))}
    edges: list[tuple[int, int, tuple[int, int, int]]] = []
    for start, end, image in graph.gain_edges:
        if start not in component or end not in component:
            continue
        new_start = node_map[start]
        new_end = node_map[end]
        if new_end < new_start:
            new_start, new_end = new_end, new_start
            image = _negate_image(image)
        edges.append((new_start, new_end, image))
    return LinkageTopologyGraph(
        node_connectivities=tuple(graph.node_connectivities[old] for old in sorted(component)),
        gain_edges=tuple(sorted(edges)),
        dimensionality_hint=graph.dimensionality_hint,
    )


def _build_bonded_mol(atoms: PeriodicCifAtoms, *, bond_mode: str = "auto") -> BondedMolBuildResult:
    labels = tuple(atoms.info.get("_atom_site_label", ()))
    if len(labels) != len(atoms):
        raise ValueError("decomposition requires CIF atom labels aligned with atom sites")

    rw_mol = Chem.RWMol()
    for label, symbol in zip(labels, atoms.symbols):
        atom = Chem.Atom(symbol)
        atom.SetProp("cif_label", str(label))
        atom.SetProp("instance_id", _instance_id(str(label)))
        rw_mol.AddAtom(atom)

    candidates, metadata = _bond_candidates_from_cif_or_geometry(atoms, labels, bond_mode=bond_mode)
    if not candidates:
        raise ValueError("decomposition could not infer any covalent bonds from CIF atom positions")

    added_pairs: set[frozenset[int]] = set()
    for candidate in candidates:
        key = frozenset((candidate.atom_idx_1, candidate.atom_idx_2))
        if key in added_pairs:
            continue
        rw_mol.AddBond(candidate.atom_idx_1, candidate.atom_idx_2, Chem.BondType.SINGLE)
        added_pairs.add(key)

    mol = rw_mol.GetMol()
    ring_bonds = _ring_bond_keys(mol)
    inferred_order_count = 0
    for candidate in candidates:
        bond = mol.GetBondBetweenAtoms(candidate.atom_idx_1, candidate.atom_idx_2)
        if bond is None:
            continue
        order = candidate.explicit_order
        if order is None:
            order = _infer_bond_order(mol, candidate, ring_bonds)
            inferred_order_count += 1
        _apply_bond_order(mol, bond, order)

    mol.UpdatePropertyCache(strict=False)
    return BondedMolBuildResult(
        mol=mol,
        metadata={
            **metadata,
            "n_bond_orders_inferred": inferred_order_count,
        },
        candidates=tuple(candidates),
    )


def detect_cif_topology(
    cif_path: str | Path,
    *,
    linkage: str = "imine",
    bond_mode: str = "auto",
    atoms: PeriodicCifAtoms | None = None,
    details: LinkageDecompositionDetails | None = None,
) -> TopologyDetectionResult:
    input_path = Path(cif_path)
    comment_cofid = read_cofid_from_cif(input_path)
    if comment_cofid:
        try:
            parsed = parse_cofid(comment_cofid)
        except ValueError:
            pass
        else:
            return TopologyDetectionResult(
                status="success",
                selected_topology=parsed.topology,
                confidence="exact",
                reason="topology recovered from embedded COFid comment",
                metadata={"source": "cofid_comment", "cofid": comment_cofid},
            )

    spec = _resolve_linkage_spec(linkage)
    if spec is None:
        return TopologyDetectionResult(
            status="failed",
            reason=f"unsupported linkage {linkage!r}",
            metadata={"source": "linkage_decomposition"},
        )

    try:
        if details is None:
            atoms = atoms or read_periodic_cif_atoms(input_path)
            details = _decompose_atoms(atoms, spec, bond_mode=bond_mode)
    except Exception as exc:
        return TopologyDetectionResult(
            status="failed",
            reason=f"{type(exc).__name__}: {exc}",
            metadata={"source": "linkage_decomposition"},
        )

    graph = details.topology_graph
    if graph is None or graph.node_count == 0 or graph.edge_count == 0:
        return TopologyDetectionResult(
            status="failed",
            reason="could not extract a periodic linkage graph from the CIF",
            metadata={
                "source": "linkage_graph",
                "decomposition": dict(details.metadata),
            },
        )

    candidates = _rank_topology_candidates(graph)
    if not candidates:
        return TopologyDetectionResult(
            status="failed",
            reason="no topology candidates matched the recovered linkage graph connectivities",
            metadata={"source": "topology_repository", "topology_graph": graph.to_metadata()},
        )

    best = candidates[0]
    runner_up = candidates[1] if len(candidates) > 1 else None
    decisive = best.confidence in {"exact", "high"} and (
        runner_up is None or best.score - runner_up.score >= 0.05
    )
    if decisive:
        return TopologyDetectionResult(
            status="success",
            selected_topology=best.topology,
            confidence=best.confidence,
            reason=best.reason,
            candidates=candidates[:10],
            metadata={"source": "periodic_linkage_graph", "topology_graph": graph.to_metadata()},
        )

    ambiguous = ", ".join(candidate.topology for candidate in candidates[:5])
    return TopologyDetectionResult(
        status="ambiguous",
        confidence=best.confidence,
        reason=f"topology auto-detection is ambiguous among: {ambiguous}",
        candidates=candidates[:10],
        metadata={"source": "periodic_linkage_graph", "topology_graph": graph.to_metadata()},
    )


def _linkage_topology_graph(
    mol,
    linkage_bond_atom_pairs: list[tuple[int, int]],
    atom_to_fragment: Mapping[int, int],
    candidates: tuple[BondCandidate, ...],
) -> LinkageTopologyGraph | None:
    candidate_by_pair = {
        frozenset((candidate.atom_idx_1, candidate.atom_idx_2)): candidate
        for candidate in candidates
    }
    participating_fragments: set[int] = set()
    raw_edges: list[tuple[int, int, tuple[int, int, int]]] = []
    for atom_idx_1, atom_idx_2 in linkage_bond_atom_pairs:
        fragment_1 = atom_to_fragment.get(atom_idx_1)
        fragment_2 = atom_to_fragment.get(atom_idx_2)
        if fragment_1 is None or fragment_2 is None or fragment_1 == fragment_2:
            continue
        candidate = candidate_by_pair.get(frozenset((atom_idx_1, atom_idx_2)))
        image = (0, 0, 0) if candidate is None else candidate.periodic_image
        if candidate is not None and (candidate.atom_idx_1, candidate.atom_idx_2) != (atom_idx_1, atom_idx_2):
            image = _negate_image(image)
        participating_fragments.update((fragment_1, fragment_2))
        raw_edges.append((fragment_1, fragment_2, image))

    if not raw_edges:
        return None

    node_map = {fragment_id: index for index, fragment_id in enumerate(sorted(participating_fragments))}
    degrees: Counter[int] = Counter()
    edges: list[tuple[int, int, tuple[int, int, int]]] = []
    for fragment_1, fragment_2, image in raw_edges:
        start = node_map[fragment_1]
        end = node_map[fragment_2]
        degrees[start] += 1
        degrees[end] += 1
        if end < start:
            start, end = end, start
            image = _negate_image(image)
        edges.append((start, end, image))

    dimensionality_hint = "2D" if all(edge[2][2] == 0 for edge in edges) else "3D"
    return LinkageTopologyGraph(
        node_connectivities=tuple(int(degrees[index]) for index in range(len(node_map))),
        gain_edges=tuple(sorted(edges)),
        dimensionality_hint=dimensionality_hint,
    )


def _rank_topology_candidates(graph: LinkageTopologyGraph) -> tuple[TopologyDetectionCandidate, ...]:
    repository = default_topology_repository()
    entries = repository.list_index(dimensionality=graph.dimensionality_hint)
    if not entries:
        entries = repository.list_index()

    observed_connectivities = tuple(sorted(set(graph.node_connectivities)))
    ranked: list[TopologyDetectionCandidate] = []
    for entry in entries:
        entry_connectivities = tuple(sorted(set(int(value) for value in entry.node_connectivities)))
        if entry_connectivities != observed_connectivities:
            continue
        score = 0.35
        reasons = ["node connectivities match"]
        metadata: dict[str, object] = {
            "node_connectivities": list(entry.node_connectivities),
            "dimensionality": entry.dimensionality,
        }
        if graph.dimensionality_hint == entry.dimensionality:
            score += 0.10
            reasons.append("dimensionality hint matches")
        try:
            topology_graph = _topology_definition_graph(repository.load(entry.id))
        except Exception as exc:
            metadata["topology_load_error"] = f"{type(exc).__name__}: {exc}"
            topology_graph = None

        confidence = "low"
        if topology_graph is not None:
            metadata["definition_graph"] = topology_graph.to_metadata()
            if graph.node_count == topology_graph.node_count and graph.edge_count == topology_graph.edge_count:
                if _periodic_graphs_match(graph, topology_graph, compare_gains=True):
                    score += 0.55
                    confidence = "exact"
                    reasons.append("periodic quotient graph matches")
                elif _periodic_graphs_match(graph, topology_graph, compare_gains=False):
                    score += 0.45
                    confidence = "high"
                    reasons.append("quotient graph matches without periodic-gain comparison")
            else:
                if topology_graph.node_count > 0 and graph.node_count % topology_graph.node_count == 0:
                    score += 0.10
                    reasons.append("observed node count is a topology supercell multiple")
                if topology_graph.edge_count > 0 and graph.edge_count % topology_graph.edge_count == 0:
                    score += 0.10
                    reasons.append("observed edge count is a topology supercell multiple")
                if _degree_histogram_matches_supercell(graph, topology_graph):
                    score += 0.20
                    reasons.append("degree histogram is consistent with a topology supercell")
        if confidence == "low":
            if score >= 0.85:
                confidence = "high"
            elif score >= 0.65:
                confidence = "medium"
        ranked.append(
            TopologyDetectionCandidate(
                topology=entry.id,
                score=round(score, 6),
                confidence=confidence,
                reason="; ".join(reasons),
                metadata=metadata,
            )
        )
    return tuple(sorted(ranked, key=lambda candidate: (-candidate.score, candidate.topology)))


def _topology_definition_graph(definition) -> LinkageTopologyGraph:
    expanded = expand_topology_definition(definition)
    node_index_by_id = {node.id: index for index, node in enumerate(expanded.node_sites)}
    edges: list[tuple[int, int, tuple[int, int, int]]] = []
    for edge in expanded.edge_sites:
        start = node_index_by_id[edge.start_node_id]
        end = node_index_by_id[edge.end_node_id]
        image = _pad_image(edge.end_image)
        if end < start:
            start, end = end, start
            image = _negate_image(image)
        edges.append((start, end, image))
    return LinkageTopologyGraph(
        node_connectivities=tuple(int(node.connectivity) for node in expanded.node_sites),
        gain_edges=tuple(sorted(edges)),
        dimensionality_hint=definition.dimensionality,
    )


def _periodic_graphs_match(
    observed: LinkageTopologyGraph,
    expected: LinkageTopologyGraph,
    *,
    compare_gains: bool,
) -> bool:
    if observed.node_count != expected.node_count or observed.edge_count != expected.edge_count:
        return False
    if sorted(observed.node_connectivities) != sorted(expected.node_connectivities):
        return False
    if observed.node_count > 10:
        return False

    expected_pair_counts = _pair_edge_counter(expected.gain_edges)
    expected_edge_counter = _edge_counter(expected.gain_edges, compare_gains=compare_gains)
    observed_nodes = sorted(
        range(observed.node_count),
        key=lambda index: (-observed.node_connectivities[index], index),
    )
    expected_by_connectivity: dict[int, list[int]] = defaultdict(list)
    for index, connectivity in enumerate(expected.node_connectivities):
        expected_by_connectivity[int(connectivity)].append(index)

    mapping: dict[int, int] = {}
    used_expected: set[int] = set()

    def backtrack(position: int) -> bool:
        if position == len(observed_nodes):
            mapped_edges = []
            for start, end, image in observed.gain_edges:
                mapped_start = mapping[start]
                mapped_end = mapping[end]
                mapped_image = image
                if mapped_end < mapped_start:
                    mapped_start, mapped_end = mapped_end, mapped_start
                    mapped_image = _negate_image(mapped_image)
                mapped_edges.append((mapped_start, mapped_end, mapped_image))
            return _edge_counter(tuple(mapped_edges), compare_gains=compare_gains) == expected_edge_counter

        observed_node = observed_nodes[position]
        connectivity = int(observed.node_connectivities[observed_node])
        for expected_node in expected_by_connectivity.get(connectivity, ()):
            if expected_node in used_expected:
                continue
            if not _partial_mapping_is_consistent(
                observed,
                expected_pair_counts,
                mapping,
                observed_node,
                expected_node,
            ):
                continue
            mapping[observed_node] = expected_node
            used_expected.add(expected_node)
            if backtrack(position + 1):
                return True
            used_expected.remove(expected_node)
            del mapping[observed_node]
        return False

    return backtrack(0)


def _partial_mapping_is_consistent(
    observed: LinkageTopologyGraph,
    expected_pair_counts: Counter[tuple[int, int]],
    mapping: Mapping[int, int],
    observed_node: int,
    expected_node: int,
) -> bool:
    for start, end, _image in observed.gain_edges:
        other = None
        if start == observed_node and end in mapping:
            other = end
        elif end == observed_node and start in mapping:
            other = start
        if other is None:
            continue
        mapped_pair = tuple(sorted((expected_node, mapping[other])))
        observed_count = sum(
            1
            for edge_start, edge_end, _ in observed.gain_edges
            if {edge_start, edge_end} == {observed_node, other}
        )
        if expected_pair_counts[mapped_pair] < observed_count:
            return False
    return True


def _degree_histogram_matches_supercell(
    observed: LinkageTopologyGraph,
    expected: LinkageTopologyGraph,
) -> bool:
    observed_counts = Counter(observed.node_connectivities)
    expected_counts = Counter(expected.node_connectivities)
    if not expected_counts:
        return False
    multiplier: int | None = None
    for connectivity, expected_count in expected_counts.items():
        observed_count = observed_counts.get(connectivity, 0)
        if observed_count == 0 or observed_count % expected_count != 0:
            return False
        current = observed_count // expected_count
        if multiplier is None:
            multiplier = current
        elif multiplier != current:
            return False
    return sum(observed_counts.values()) == sum(expected_counts.values()) * (multiplier or 0)


def _pair_edge_counter(edges: tuple[tuple[int, int, tuple[int, int, int]], ...]) -> Counter[tuple[int, int]]:
    return Counter(tuple(sorted((start, end))) for start, end, _image in edges)


def _edge_counter(
    edges: tuple[tuple[int, int, tuple[int, int, int]], ...],
    *,
    compare_gains: bool,
) -> Counter[tuple[object, ...]]:
    if compare_gains:
        return Counter((start, end, image) for start, end, image in edges)
    return Counter((start, end) for start, end, _image in edges)


def _pad_image(image: tuple[int, ...]) -> tuple[int, int, int]:
    if len(image) >= 3:
        return (int(image[0]), int(image[1]), int(image[2]))
    if len(image) == 2:
        return (int(image[0]), int(image[1]), 0)
    if len(image) == 1:
        return (int(image[0]), 0, 0)
    return (0, 0, 0)


def _negate_image(image: tuple[int, int, int]) -> tuple[int, int, int]:
    return (-image[0], -image[1], -image[2])


def _bond_candidates_from_cif_or_geometry(
    atoms: PeriodicCifAtoms,
    labels: tuple[str, ...],
    *,
    bond_mode: str,
) -> tuple[tuple[BondCandidate, ...], dict[str, object]]:
    label_to_idx = {label: index for index, label in enumerate(labels)}
    bond_labels_1 = tuple(atoms.info.get("_geom_bond_atom_site_label_1", ()))
    bond_labels_2 = tuple(atoms.info.get("_geom_bond_atom_site_label_2", ()))
    if bond_mode == "auto" and bond_labels_1 and bond_labels_2:
        candidates, n_missing_orders = _explicit_bond_candidates(atoms, label_to_idx, bond_labels_1, bond_labels_2)
        return candidates, {
            "bond_mode": bond_mode,
            "bond_source": "explicit_cif",
            "n_explicit_cif_bonds": len(candidates),
            "n_missing_cif_bond_orders": n_missing_orders,
            "n_distance_inferred_bonds": 0,
        }

    candidates = _distance_inferred_bond_candidates(atoms)
    return candidates, {
        "bond_mode": bond_mode,
        "bond_source": "distance_inferred",
        "n_explicit_cif_bonds": 0,
        "n_missing_cif_bond_orders": 0,
        "n_distance_inferred_bonds": len(candidates),
    }


def _explicit_bond_candidates(
    atoms: PeriodicCifAtoms,
    label_to_idx: Mapping[str, int],
    bond_labels_1: tuple[str, ...],
    bond_labels_2: tuple[str, ...],
) -> tuple[tuple[BondCandidate, ...], int]:
    bond_types = tuple(atoms.info.get("_ccdc_geom_bond_type") or atoms.info.get("_geom_bond_type") or ())
    bond_distances = tuple(atoms.info.get("_geom_bond_distance") or ())
    symmetries_1 = tuple(atoms.info.get("_geom_bond_site_symmetry_1") or ())
    symmetries_2 = tuple(atoms.info.get("_geom_bond_site_symmetry_2") or ())
    candidates: list[BondCandidate] = []
    added_pairs: set[frozenset[int]] = set()
    for row_idx, (label_1, label_2) in enumerate(zip(bond_labels_1, bond_labels_2)):
        if label_1 not in label_to_idx or label_2 not in label_to_idx:
            continue
        idx_1 = label_to_idx[label_1]
        idx_2 = label_to_idx[label_2]
        if idx_1 == idx_2:
            continue
        key = frozenset((idx_1, idx_2))
        if key in added_pairs:
            continue
        order = cif_type_to_bond_order(bond_types[row_idx] if row_idx < len(bond_types) else None)
        distance = _explicit_bond_distance(atoms, idx_1, idx_2, bond_distances, row_idx)
        image = _explicit_bond_periodic_image(symmetries_1, symmetries_2, row_idx)
        candidates.append(BondCandidate(idx_1, idx_2, distance, explicit_order=order, periodic_image=image))
        added_pairs.add(key)
    return tuple(candidates), sum(1 for candidate in candidates if candidate.explicit_order is None)


def _explicit_bond_periodic_image(
    symmetries_1: tuple[str, ...],
    symmetries_2: tuple[str, ...],
    row_idx: int,
) -> tuple[int, int, int]:
    first = _parse_p1_symmetry_shift(symmetries_1[row_idx] if row_idx < len(symmetries_1) else None)
    second = _parse_p1_symmetry_shift(symmetries_2[row_idx] if row_idx < len(symmetries_2) else None)
    return (second[0] - first[0], second[1] - first[1], second[2] - first[2])


def _parse_p1_symmetry_shift(value: object | None) -> tuple[int, int, int]:
    text = "" if value is None else str(value).strip().strip("'\"")
    if not text or text in {".", "?", "1"}:
        return (0, 0, 0)
    match = re.match(r"^1_([0-9])([0-9])([0-9])$", text)
    if match is None:
        return (0, 0, 0)
    return tuple(int(component) - 5 for component in match.groups())  # type: ignore[return-value]


def _explicit_bond_distance(
    atoms: PeriodicCifAtoms,
    idx_1: int,
    idx_2: int,
    bond_distances: tuple[str, ...],
    row_idx: int,
) -> float:
    if row_idx < len(bond_distances):
        try:
            return _parse_cif_float(bond_distances[row_idx])
        except ValueError:
            pass
    return _minimum_image_distance(atoms, idx_1, idx_2)


def _distance_inferred_bond_candidates(atoms: PeriodicCifAtoms) -> tuple[BondCandidate, ...]:
    candidates: list[BondCandidate] = []
    for idx_1, symbol_1 in enumerate(atoms.symbols):
        for idx_2 in range(idx_1 + 1, len(atoms.symbols)):
            symbol_2 = atoms.symbols[idx_2]
            image, distance = _minimum_image_bond_geometry(atoms, idx_1, idx_2)
            if _is_plausible_bond_distance(symbol_1, symbol_2, distance):
                candidates.append(BondCandidate(idx_1, idx_2, distance, explicit_order=None, periodic_image=image))
    return _prune_distance_inferred_bond_candidates(atoms, tuple(candidates))


def _prune_distance_inferred_bond_candidates(
    atoms: PeriodicCifAtoms,
    candidates: tuple[BondCandidate, ...],
) -> tuple[BondCandidate, ...]:
    if not candidates:
        return ()
    active: set[int] = set(range(len(candidates)))
    incident: dict[int, set[int]] = defaultdict(set)
    for candidate_index, candidate in enumerate(candidates):
        incident[candidate.atom_idx_1].add(candidate_index)
        incident[candidate.atom_idx_2].add(candidate_index)

    changed = True
    while changed:
        changed = False
        for atom_idx, symbol in enumerate(atoms.symbols):
            max_degree = _max_distance_inferred_degree(symbol)
            active_incident = [candidate_index for candidate_index in incident.get(atom_idx, ()) if candidate_index in active]
            while len(active_incident) > max_degree:
                remove_index = max(
                    active_incident,
                    key=lambda candidate_index: (
                        candidates[candidate_index].distance,
                        _distance_prune_priority(atoms, candidates[candidate_index], atom_idx),
                    ),
                )
                active.remove(remove_index)
                active_incident.remove(remove_index)
                changed = True
    return tuple(candidate for candidate_index, candidate in enumerate(candidates) if candidate_index in active)


def _max_distance_inferred_degree(symbol: str) -> int:
    atomic_number = _atomic_number(symbol)
    return {
        1: 1,
        5: 3,
        6: 4,
        7: 4,
        8: 2,
        9: 1,
        15: 5,
        16: 6,
        17: 1,
        35: 1,
        53: 1,
    }.get(atomic_number, 6)


def _distance_prune_priority(
    atoms: PeriodicCifAtoms,
    candidate: BondCandidate,
    atom_idx: int,
) -> int:
    other_idx = candidate.atom_idx_2 if candidate.atom_idx_1 == atom_idx else candidate.atom_idx_1
    pair = frozenset((_atomic_number(atoms.symbols[atom_idx]), _atomic_number(atoms.symbols[other_idx])))
    if pair in {frozenset((6, 7)), frozenset((6, 8)), frozenset((6, 6))}:
        return 0
    return 1


def _ring_bond_keys(mol) -> set[frozenset[int]]:
    ring_bonds: set[frozenset[int]] = set()
    Chem.GetSymmSSSR(mol)
    for bond_ring in mol.GetRingInfo().BondRings():
        if len(bond_ring) not in {5, 6}:
            continue
        for bond_idx in bond_ring:
            bond = mol.GetBondWithIdx(int(bond_idx))
            ring_bonds.add(frozenset((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())))
    return ring_bonds


def _infer_bond_order(mol, candidate: BondCandidate, ring_bonds: set[frozenset[int]]) -> float:
    atom_1 = mol.GetAtomWithIdx(candidate.atom_idx_1)
    atom_2 = mol.GetAtomWithIdx(candidate.atom_idx_2)
    z1 = atom_1.GetAtomicNum()
    z2 = atom_2.GetAtomicNum()
    pair = frozenset((z1, z2))
    distance = candidate.distance
    ring_key = frozenset((candidate.atom_idx_1, candidate.atom_idx_2))

    instance_1 = atom_1.GetProp("instance_id")
    instance_2 = atom_2.GetProp("instance_id")
    same_instance = bool(instance_1 and instance_2 and instance_1 == instance_2)
    known_different_instance = bool(instance_1 and instance_2 and instance_1 != instance_2)
    if pair == frozenset((6, 7)):
        if distance <= 1.22:
            return 3.0
        if known_different_instance and distance <= 1.48:
            return 2.0
        if not same_instance and distance <= 1.34:
            return 2.0
    if ring_key in ring_bonds and pair in {frozenset((6, 6)), frozenset((6, 7)), frozenset((6, 16))}:
        if 1.20 <= distance <= 1.50:
            return 1.5

    if 1 in pair:
        return 1.0
    if 5 in pair:
        return 1.0
    if pair == frozenset((6, 8)):
        return 2.0 if distance <= 1.30 else 1.0
    if pair == frozenset((6, 7)):
        return 2.0 if distance <= 1.28 else 1.0
    if pair == frozenset((6, 6)):
        if distance <= 1.24:
            return 3.0
        if known_different_instance and distance <= 1.46:
            return 2.0
        return 2.0 if distance <= 1.37 else 1.0
    return 1.0


def _apply_bond_order(mol, bond, order: float) -> None:
    bond.SetBondType(_rdkit_bond_type(order))
    if is_aromatic_bond_order(order):
        bond.SetIsAromatic(True)
        bond.GetBeginAtom().SetIsAromatic(True)
        bond.GetEndAtom().SetIsAromatic(True)


def _eligible_carbon_nitrogen_double_bonds(mol) -> dict[int, tuple[int, int]]:
    eligible: dict[int, tuple[int, int]] = {}
    for bond in mol.GetBonds():
        if abs(float(bond.GetBondTypeAsDouble()) - 2.0) > 1.0e-6:
            continue
        first = bond.GetBeginAtom()
        second = bond.GetEndAtom()
        if {first.GetAtomicNum(), second.GetAtomicNum()} != {6, 7}:
            continue
        first_instance = first.GetProp("instance_id")
        second_instance = second.GetProp("instance_id")
        if first_instance and second_instance and first_instance == second_instance:
            continue
        carbon = first if first.GetAtomicNum() == 6 else second
        nitrogen = first if first.GetAtomicNum() == 7 else second
        eligible[bond.GetIdx()] = (carbon.GetIdx(), nitrogen.GetIdx())
    return eligible


def _classify_nitrogen_linkage_bonds(mol) -> NitrogenLinkageBondClassification:
    eligible = _eligible_carbon_nitrogen_double_bonds(mol)
    eligible_by_nitrogen: dict[int, list[int]] = defaultdict(list)
    for bond_idx, (_carbon_idx, nitrogen_idx) in eligible.items():
        eligible_by_nitrogen[nitrogen_idx].append(bond_idx)

    azine_candidates: set[int] = set()
    hydrazone_candidates: set[int] = set()
    bken_candidates: set[int] = set()
    for bond_idx, (carbon_idx, nitrogen_idx) in eligible.items():
        nitrogen = mol.GetAtomWithIdx(nitrogen_idx)
        for neighbor in nitrogen.GetNeighbors():
            if neighbor.GetAtomicNum() != 7:
                continue
            nn_bond = mol.GetBondBetweenAtoms(nitrogen_idx, neighbor.GetIdx())
            if nn_bond is None or abs(float(nn_bond.GetBondTypeAsDouble()) - 1.0) > 1.0e-6:
                continue
            partner_cn_bonds = eligible_by_nitrogen.get(neighbor.GetIdx(), ())
            if partner_cn_bonds:
                azine_candidates.add(bond_idx)
                azine_candidates.update(partner_cn_bonds)
                continue
            if _nitrogen_is_bonded_to_carbonyl(mol, neighbor.GetIdx(), exclude_atom_idx=nitrogen_idx):
                hydrazone_candidates.add(bond_idx)

        if _carbon_has_beta_ketoenamine_environment(mol, carbon_idx, exclude_atom_idx=nitrogen_idx):
            bken_candidates.add(bond_idx)

    assigned_azine: set[int] = set()
    assigned_hydrazone: set[int] = set()
    assigned_bken: set[int] = set()
    assigned_imine: set[int] = set()
    cross_branch_conflicts: set[int] = set()
    for bond_idx in eligible:
        nn_assignment = (
            "azine"
            if bond_idx in azine_candidates
            else "hydrazone"
            if bond_idx in hydrazone_candidates
            else None
        )
        keto_assignment = "bken" if bond_idx in bken_candidates else None
        if nn_assignment is not None and keto_assignment is not None:
            cross_branch_conflicts.add(bond_idx)
        elif nn_assignment == "azine":
            assigned_azine.add(bond_idx)
        elif nn_assignment == "hydrazone":
            assigned_hydrazone.add(bond_idx)
        elif keto_assignment == "bken":
            assigned_bken.add(bond_idx)
        else:
            assigned_imine.add(bond_idx)

    return NitrogenLinkageBondClassification(
        raw_cn_bond_indices=tuple(sorted(eligible)),
        azine_bond_indices=tuple(sorted(assigned_azine)),
        hydrazone_bond_indices=tuple(sorted(assigned_hydrazone)),
        bken_bond_indices=tuple(sorted(assigned_bken)),
        imine_bond_indices=tuple(sorted(assigned_imine)),
        cross_branch_conflict_bond_indices=tuple(sorted(cross_branch_conflicts)),
    )


def _nitrogen_is_bonded_to_carbonyl(mol, nitrogen_idx: int, *, exclude_atom_idx: int) -> bool:
    nitrogen = mol.GetAtomWithIdx(nitrogen_idx)
    for carbon in nitrogen.GetNeighbors():
        if carbon.GetIdx() == exclude_atom_idx or carbon.GetAtomicNum() != 6:
            continue
        nc_bond = mol.GetBondBetweenAtoms(nitrogen_idx, carbon.GetIdx())
        if nc_bond is None or abs(float(nc_bond.GetBondTypeAsDouble()) - 1.0) > 1.0e-6:
            continue
        if _atom_has_double_bonded_oxygen(mol, carbon.GetIdx(), exclude_atom_idx=nitrogen_idx):
            return True
    return False


def _carbon_has_beta_ketoenamine_environment(mol, carbon_idx: int, *, exclude_atom_idx: int) -> bool:
    carbon = mol.GetAtomWithIdx(carbon_idx)
    for ring_anchor in carbon.GetNeighbors():
        if ring_anchor.GetIdx() == exclude_atom_idx or ring_anchor.GetAtomicNum() != 6:
            continue
        anchor_bond = mol.GetBondBetweenAtoms(carbon_idx, ring_anchor.GetIdx())
        if anchor_bond is None or abs(float(anchor_bond.GetBondTypeAsDouble()) - 1.0) > 1.0e-6:
            continue
        if not ring_anchor.IsInRing():
            continue
        for carbonyl in ring_anchor.GetNeighbors():
            if carbonyl.GetIdx() == carbon_idx or carbonyl.GetAtomicNum() != 6:
                continue
            ring_bond = mol.GetBondBetweenAtoms(ring_anchor.GetIdx(), carbonyl.GetIdx())
            if ring_bond is None or not ring_bond.IsInRing() or not carbonyl.IsInRing():
                continue
            if _atom_has_double_bonded_oxygen(mol, carbonyl.GetIdx(), exclude_atom_idx=ring_anchor.GetIdx()):
                return True
    return False


def _atom_has_double_bonded_oxygen(mol, atom_idx: int, *, exclude_atom_idx: int) -> bool:
    atom = mol.GetAtomWithIdx(atom_idx)
    for oxygen in atom.GetNeighbors():
        if oxygen.GetIdx() == exclude_atom_idx or oxygen.GetAtomicNum() != 8:
            continue
        bond = mol.GetBondBetweenAtoms(atom_idx, oxygen.GetIdx())
        if bond is not None and abs(float(bond.GetBondTypeAsDouble()) - 2.0) <= 1.0e-6:
            return True
    return False


def _mark_classified_nitrogen_linkage_bonds(
    mol,
    *,
    linkage: str,
    carbon_role: str,
    nitrogen_role: str,
) -> tuple[int, ...]:
    classification = _classify_nitrogen_linkage_bonds(mol)
    if classification.cross_branch_conflict_bond_indices:
        return ()
    bond_indices = classification.bond_indices_for(linkage)
    for bond_idx in bond_indices:
        bond = mol.GetBondWithIdx(bond_idx)
        first = bond.GetBeginAtom()
        second = bond.GetEndAtom()
        carbon = first if first.GetAtomicNum() == 6 else second
        hetero = first if first.GetAtomicNum() == 7 else second
        carbon.SetProp("cofkit_decompose_role", carbon_role)
        hetero.SetProp("cofkit_decompose_role", nitrogen_role)
    return bond_indices


def _mark_imine_linkage_bonds(mol) -> tuple[int, ...]:
    return _mark_classified_nitrogen_linkage_bonds(
        mol,
        linkage="imine",
        carbon_role="aldehyde",
        nitrogen_role="amine",
    )


def _mark_hydrazone_linkage_bonds(mol) -> tuple[int, ...]:
    return _mark_classified_nitrogen_linkage_bonds(
        mol,
        linkage="hydrazone",
        carbon_role="aldehyde",
        nitrogen_role="hydrazide",
    )


def _mark_azine_linkage_bonds(mol) -> tuple[int, ...]:
    return _mark_classified_nitrogen_linkage_bonds(
        mol,
        linkage="azine",
        carbon_role="aldehyde",
        nitrogen_role="hydrazine",
    )


def _mark_keto_enamine_linkage_bonds(mol) -> tuple[int, ...]:
    return _mark_classified_nitrogen_linkage_bonds(
        mol,
        linkage="bken",
        carbon_role="keto_aldehyde",
        nitrogen_role="amine",
    )


def _mark_vinylene_linkage_bonds(mol) -> tuple[int, ...]:
    classification = _classify_vinylene_linkage_bonds(mol)
    for _bond_idx, aldehyde_idx, activated_idx in classification.accepted_assignments:
        mol.GetAtomWithIdx(aldehyde_idx).SetProp("cofkit_decompose_role", "aldehyde")
        mol.GetAtomWithIdx(activated_idx).SetProp("cofkit_decompose_role", "activated_methylene")
    return classification.accepted_bond_indices


def _classify_vinylene_linkage_bonds(mol) -> VinyleneLinkageBondClassification:
    formal_cc_double_bonds: list[int] = []
    same_instance_rejected: list[int] = []
    small_ring_rejected: list[int] = []
    endpoint_rejected: list[int] = []
    accepted: list[tuple[int, int, int]] = []
    small_ring_bond_indices = {
        int(bond_idx)
        for ring in mol.GetRingInfo().BondRings()
        if len(ring) in {5, 6}
        for bond_idx in ring
    }
    for bond in mol.GetBonds():
        if abs(float(bond.GetBondTypeAsDouble()) - 2.0) > 1.0e-6:
            continue
        first = bond.GetBeginAtom()
        second = bond.GetEndAtom()
        if first.GetAtomicNum() != 6 or second.GetAtomicNum() != 6:
            continue
        bond_idx = bond.GetIdx()
        formal_cc_double_bonds.append(bond_idx)
        first_instance = first.GetProp("instance_id")
        second_instance = second.GetProp("instance_id")
        if first_instance and second_instance and first_instance == second_instance:
            same_instance_rejected.append(bond_idx)
            continue
        if bond_idx in small_ring_bond_indices:
            small_ring_rejected.append(bond_idx)
            continue
        first_idx = first.GetIdx()
        second_idx = second.GetIdx()
        first_anchors = _vinylene_carbon_anchor_indices(mol, first_idx, exclude_atom_idx=second_idx)
        second_anchors = _vinylene_carbon_anchor_indices(mol, second_idx, exclude_atom_idx=first_idx)
        if not first_anchors or not second_anchors:
            endpoint_rejected.append(bond_idx)
            continue
        first_activation = _vinylene_activation_score(mol, first_anchors, exclude_atom_idx=first_idx)
        second_activation = _vinylene_activation_score(mol, second_anchors, exclude_atom_idx=second_idx)
        if max(first_activation, second_activation) <= 0 or first_activation == second_activation:
            endpoint_rejected.append(bond_idx)
            continue
        if first_activation > second_activation:
            activated_idx, aldehyde_idx = first_idx, second_idx
        else:
            activated_idx, aldehyde_idx = second_idx, first_idx
        accepted.append((bond_idx, aldehyde_idx, activated_idx))
    return VinyleneLinkageBondClassification(
        formal_cc_double_bond_indices=tuple(sorted(formal_cc_double_bonds)),
        same_instance_rejected_bond_indices=tuple(sorted(same_instance_rejected)),
        small_ring_rejected_bond_indices=tuple(sorted(small_ring_rejected)),
        endpoint_rejected_bond_indices=tuple(sorted(endpoint_rejected)),
        accepted_assignments=tuple(sorted(accepted)),
    )


def _vinylene_carbon_anchor_indices(mol, atom_idx: int, *, exclude_atom_idx: int) -> tuple[int, ...]:
    return tuple(
        sorted(
            neighbor.GetIdx()
            for neighbor in mol.GetAtomWithIdx(atom_idx).GetNeighbors()
            if neighbor.GetIdx() != exclude_atom_idx and neighbor.GetAtomicNum() == 6
        )
    )


def _vinylene_activation_score(mol, anchor_indices: tuple[int, ...], *, exclude_atom_idx: int) -> int:
    best = 0
    for anchor_idx in anchor_indices:
        anchor = mol.GetAtomWithIdx(anchor_idx)
        nitrogen_neighbors = sum(
            neighbor.GetAtomicNum() == 7
            for neighbor in anchor.GetNeighbors()
            if neighbor.GetIdx() != exclude_atom_idx
        )
        if nitrogen_neighbors >= 2:
            best = max(best, 4)
        for bond in anchor.GetBonds():
            neighbor = bond.GetOtherAtom(anchor)
            if neighbor.GetIdx() == exclude_atom_idx:
                continue
            bond_order = float(bond.GetBondTypeAsDouble())
            if bond_order >= 3.0 and neighbor.GetAtomicNum() == 7:
                best = max(best, 4)
            elif bond_order >= 2.0 and neighbor.GetAtomicNum() in {7, 8, 16}:
                best = max(best, 3)
    return best


def _mark_boronate_ester_linkage_bonds(mol) -> tuple[int, ...]:
    bond_indices: list[int] = []
    for bond in mol.GetBonds():
        if abs(float(bond.GetBondTypeAsDouble()) - 1.0) > 1.0e-6:
            continue
        first = bond.GetBeginAtom()
        second = bond.GetEndAtom()
        atomic_nums = {first.GetAtomicNum(), second.GetAtomicNum()}
        if atomic_nums != {5, 8}:
            continue
        first_instance = first.GetProp("instance_id")
        second_instance = second.GetProp("instance_id")
        if first_instance and second_instance and first_instance == second_instance:
            continue
        boron = first if first.GetAtomicNum() == 5 else second
        oxygen = first if first.GetAtomicNum() == 8 else second
        boron.SetProp("cofkit_decompose_role", "boronic_acid")
        oxygen.SetProp("cofkit_decompose_role", "catechol")
        bond_indices.append(bond.GetIdx())
    return tuple(bond_indices)


def _repair_linkage_fragment_to_monomer(
    fragment,
    spec: LinkageDecompositionSpec,
) -> DecomposedMonomer | None:
    endpoint_roles = {
        atom.GetProp("cofkit_decompose_role")
        for atom in fragment.GetAtoms()
        if atom.HasProp("cofkit_decompose_role")
    }
    if len(endpoint_roles) != 1:
        return None
    reactive_group = next(iter(endpoint_roles))
    if reactive_group not in spec.roles:
        return None
    repaired = spec.repairers.get(reactive_group, _identity_fragment_repairer)(fragment)
    return _finalize_repaired_fragment(
        repaired,
        reactive_group,
        connectivity_counter=spec.connectivity_counter,
    )


def _identity_fragment_repairer(fragment):
    return fragment


def _restore_aldehyde_oxygens(fragment):
    return _restore_double_oxygens(fragment, "aldehyde")


def _restore_keto_aldehyde_fragment(fragment):
    editable = Chem.RWMol(fragment)
    for bond in tuple(editable.GetBonds()):
        if bond.GetBondType() != Chem.BondType.DOUBLE:
            continue
        first = bond.GetBeginAtom()
        second = bond.GetEndAtom()
        if {first.GetAtomicNum(), second.GetAtomicNum()} != {6, 8}:
            continue
        carbon = first if first.GetAtomicNum() == 6 else second
        if carbon.GetIsAromatic():
            bond.SetBondType(Chem.BondType.SINGLE)
    return _restore_double_oxygens(editable.GetMol(), "keto_aldehyde")


def _restore_boronic_acid_oxygens(fragment):
    editable = Chem.RWMol(fragment)
    endpoint_indices = [
        atom.GetIdx()
        for atom in editable.GetAtoms()
        if atom.HasProp("cofkit_decompose_role") and atom.GetProp("cofkit_decompose_role") == "boronic_acid"
    ]
    for boron_idx in endpoint_indices:
        boron = editable.GetAtomWithIdx(boron_idx)
        oxygen_neighbor_count = sum(1 for neighbor in boron.GetNeighbors() if neighbor.GetAtomicNum() == 8)
        for _ in range(max(0, 2 - oxygen_neighbor_count)):
            oxygen_idx = editable.AddAtom(Chem.Atom("O"))
            editable.AddBond(boron_idx, oxygen_idx, Chem.BondType.SINGLE)
    return editable.GetMol()


def _restore_nitrile_nitrogens(fragment):
    editable = Chem.RWMol(fragment)
    endpoint_indices = [
        atom.GetIdx()
        for atom in editable.GetAtoms()
        if atom.HasProp("cofkit_decompose_role") and atom.GetProp("cofkit_decompose_role") == "nitrile"
    ]
    for carbon_idx in endpoint_indices:
        carbon = editable.GetAtomWithIdx(carbon_idx)
        carbon.SetIsAromatic(False)
        nitrogen_idx = editable.AddAtom(Chem.Atom("N"))
        editable.AddBond(carbon_idx, nitrogen_idx, Chem.BondType.TRIPLE)
    return editable.GetMol()


def _restore_double_oxygens(fragment, role: str):
    editable = Chem.RWMol(fragment)
    endpoint_indices = [
        atom.GetIdx()
        for atom in editable.GetAtoms()
        if atom.HasProp("cofkit_decompose_role") and atom.GetProp("cofkit_decompose_role") == role
    ]
    for carbon_idx in endpoint_indices:
        carbon = editable.GetAtomWithIdx(carbon_idx)
        if _has_double_bonded_oxygen(carbon):
            continue
        oxygen_idx = editable.AddAtom(Chem.Atom("O"))
        editable.AddBond(carbon_idx, oxygen_idx, Chem.BondType.DOUBLE)
    return editable.GetMol()


def _finalize_repaired_fragment(
    fragment,
    reactive_group: str,
    *,
    connectivity_counter: ConnectivityCounter | None = None,
) -> DecomposedMonomer | None:
    endpoint_count = (
        connectivity_counter(fragment, reactive_group)
        if connectivity_counter is not None
        else _default_connectivity_count(fragment, reactive_group)
    )
    if endpoint_count <= 0:
        return None

    mol = Chem.Mol(fragment)
    for atom in mol.GetAtoms():
        atom.SetFormalCharge(0)
        atom.SetNumRadicalElectrons(0)
        if atom.HasProp("cofkit_decompose_role"):
            atom.ClearProp("cofkit_decompose_role")
        if atom.HasProp("cif_label"):
            atom.ClearProp("cif_label")
        if atom.HasProp("instance_id"):
            atom.ClearProp("instance_id")
    Chem.SanitizeMol(mol)
    mol = Chem.RemoveHs(mol)
    Chem.SanitizeMol(mol)
    return DecomposedMonomer(
        connectivity=endpoint_count,
        reactive_group=reactive_group,
        canonical_smiles=canonicalize_smiles(Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)),
    )


def _default_connectivity_count(fragment, reactive_group: str) -> int:
    return sum(
        1
        for atom in fragment.GetAtoms()
        if atom.HasProp("cofkit_decompose_role") and atom.GetProp("cofkit_decompose_role") == reactive_group
    )


def _connectivity_count(fragment, reactive_group: str) -> int:
    count = _default_connectivity_count(fragment, reactive_group)
    if reactive_group == "catechol":
        if count % 2 != 0:
            raise ValueError("catechol decomposition found an odd number of marked catechol oxygens")
        return count // 2
    return count


def _has_double_bonded_oxygen(atom) -> bool:
    return any(
        bond.GetBondType() == Chem.BondType.DOUBLE and bond.GetOtherAtom(atom).GetAtomicNum() == 8
        for bond in atom.GetBonds()
    )


def _minimum_image_distance(atoms: PeriodicCifAtoms, idx_1: int, idx_2: int) -> float:
    return _minimum_image_bond_geometry(atoms, idx_1, idx_2)[1]


def _minimum_image_bond_geometry(
    atoms: PeriodicCifAtoms,
    idx_1: int,
    idx_2: int,
) -> tuple[tuple[int, int, int], float]:
    frac_1 = atoms.fractional_positions[idx_1]
    frac_2 = atoms.fractional_positions[idx_2]
    delta = tuple(frac_2[axis] - frac_1[axis] for axis in range(3))
    image = tuple(-int(round(value)) for value in delta)
    minimum_delta = tuple(delta[axis] + image[axis] for axis in range(3))
    cartesian = _fractional_delta_to_cartesian(atoms.cell_basis, minimum_delta)  # type: ignore[arg-type]
    return image, _norm(cartesian)


def _fractional_delta_to_cartesian(cell_basis: tuple[Vec3, Vec3, Vec3], delta: Vec3) -> Vec3:
    return (
        cell_basis[0][0] * delta[0] + cell_basis[1][0] * delta[1] + cell_basis[2][0] * delta[2],
        cell_basis[0][1] * delta[0] + cell_basis[1][1] * delta[1] + cell_basis[2][1] * delta[2],
        cell_basis[0][2] * delta[0] + cell_basis[1][2] * delta[1] + cell_basis[2][2] * delta[2],
    )


def _norm(vector: Vec3) -> float:
    return sqrt(vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2])


def _is_plausible_bond_distance(symbol_1: str, symbol_2: str, distance: float) -> bool:
    if distance < 0.35:
        return False
    z1 = _atomic_number(symbol_1)
    z2 = _atomic_number(symbol_2)
    if z1 == 1 and z2 == 1:
        return False
    if {z1, z2} == {5, 8}:
        return distance <= 2.05
    radius_sum = _covalent_radius(z1) + _covalent_radius(z2)
    tolerance = 0.30 if 1 in {z1, z2} else 0.45
    return distance <= radius_sum + tolerance


def _atomic_number(symbol: str) -> int:
    periodic_table = Chem.GetPeriodicTable()
    atomic_number = int(periodic_table.GetAtomicNumber(str(symbol)))
    if atomic_number <= 0:
        raise ValueError(f"cannot infer atomic number for element symbol {symbol!r}")
    return atomic_number


def _covalent_radius(atomic_number: int) -> float:
    radius = float(Chem.GetPeriodicTable().GetRcovalent(atomic_number))
    if radius > 0:
        return radius
    fallback = {
        1: 0.31,
        5: 0.84,
        6: 0.76,
        7: 0.71,
        8: 0.66,
        9: 0.57,
        15: 1.07,
        16: 1.05,
        17: 1.02,
    }
    return fallback.get(atomic_number, 0.75)


def _parse_cif_float(value: object) -> float:
    text = str(value).strip().strip("'\"")
    match = _CIF_NUMBER_RE.match(text)
    if match is None:
        raise ValueError(f"invalid CIF numeric value {value!r}")
    return float(match.group(1))


def _rdkit_bond_type(order: float):
    if is_aromatic_bond_order(order):
        return Chem.BondType.AROMATIC
    if abs(order - 2.0) <= 1.0e-6:
        return Chem.BondType.DOUBLE
    if abs(order - 3.0) <= 1.0e-6:
        return Chem.BondType.TRIPLE
    return Chem.BondType.SINGLE


def _instance_id(label: str) -> str:
    if "_" not in label:
        return ""
    return label.rsplit("_", 1)[0]


def _resolve_linkage_spec(linkage: str) -> DecompositionSpec | None:
    return _DECOMPOSITION_LINKAGE_ALIASES.get(str(linkage).strip().lower())


_IMINE_SPEC = LinkageDecompositionSpec(
    linkage_code="imine",
    template_id="imine_bridge",
    roles=("amine", "aldehyde"),
    marker=_mark_imine_linkage_bonds,
    repairers={"aldehyde": _restore_aldehyde_oxygens},
    connectivity_counter=_connectivity_count,
)

_HYDRAZONE_SPEC = LinkageDecompositionSpec(
    linkage_code="hydrazone",
    template_id="hydrazone_bridge",
    roles=("hydrazide", "aldehyde"),
    marker=_mark_hydrazone_linkage_bonds,
    repairers={"aldehyde": _restore_aldehyde_oxygens},
    connectivity_counter=_connectivity_count,
)

_AZINE_SPEC = LinkageDecompositionSpec(
    linkage_code="azine",
    template_id="azine_bridge",
    roles=("hydrazine", "aldehyde"),
    marker=_mark_azine_linkage_bonds,
    repairers={"aldehyde": _restore_aldehyde_oxygens},
    connectivity_counter=_connectivity_count,
)

_BORONATE_ESTER_SPEC = LinkageDecompositionSpec(
    linkage_code="boest",
    template_id="boronate_ester_bridge",
    roles=("boronic_acid", "catechol"),
    marker=_mark_boronate_ester_linkage_bonds,
    repairers={"boronic_acid": _restore_boronic_acid_oxygens},
    connectivity_counter=_connectivity_count,
)

_KETO_ENAMINE_SPEC = LinkageDecompositionSpec(
    linkage_code="bken",
    template_id="keto_enamine_bridge",
    roles=("amine", "keto_aldehyde"),
    marker=_mark_keto_enamine_linkage_bonds,
    repairers={"keto_aldehyde": _restore_keto_aldehyde_fragment},
    connectivity_counter=_connectivity_count,
)

_VINYLENE_SPEC = LinkageDecompositionSpec(
    linkage_code="vinylene",
    template_id="vinylene_bridge",
    roles=("activated_methylene", "aldehyde"),
    marker=_mark_vinylene_linkage_bonds,
    repairers={"aldehyde": _restore_aldehyde_oxygens},
    connectivity_counter=_connectivity_count,
)

_BOROXINE_SPEC = RingDecompositionSpec(
    linkage_code="boroxine",
    template_id="boroxine_trimerization",
    reactive_group="boronic_acid",
    anchor_atomic_num=5,
    hetero_atomic_num=8,
    repairer=_restore_boronic_acid_oxygens,
)

_TRIAZINE_SPEC = RingDecompositionSpec(
    linkage_code="triazine",
    template_id="triazine_trimerization",
    reactive_group="nitrile",
    anchor_atomic_num=6,
    hetero_atomic_num=7,
    repairer=_restore_nitrile_nitrogens,
)

_DECOMPOSITION_LINKAGE_ALIASES: Mapping[str, DecompositionSpec] = {
    "imine": _IMINE_SPEC,
    "imine_bridge": _IMINE_SPEC,
    "hydrazone": _HYDRAZONE_SPEC,
    "hydrazone_bridge": _HYDRAZONE_SPEC,
    "azine": _AZINE_SPEC,
    "azine_bridge": _AZINE_SPEC,
    "boest": _BORONATE_ESTER_SPEC,
    "boronate_ester": _BORONATE_ESTER_SPEC,
    "boronate_ester_bridge": _BORONATE_ESTER_SPEC,
    "bken": _KETO_ENAMINE_SPEC,
    "beta-ketoenamine": _KETO_ENAMINE_SPEC,
    "beta_ketoenamine": _KETO_ENAMINE_SPEC,
    "keto_enamine_bridge": _KETO_ENAMINE_SPEC,
    "vinylene": _VINYLENE_SPEC,
    "vinylene_bridge": _VINYLENE_SPEC,
    "boroxine": _BOROXINE_SPEC,
    "boroxine_trimerization": _BOROXINE_SPEC,
    "triazine": _TRIAZINE_SPEC,
    "triazine_trimerization": _TRIAZINE_SPEC,
}


__all__ = [
    "CifDecompositionResult",
    "DecomposedMonomer",
    "TopologyDetectionCandidate",
    "TopologyDetectionResult",
    "detect_cif_topology",
    "decompose_cif_to_cofid",
]
