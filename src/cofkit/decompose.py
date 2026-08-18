"""CIF-to-COFid decomposition helpers.

The decomposition approach here is adapted from the
deCOFpose project: https://github.com/r-fedorov/deCOFpose
"""

from __future__ import annotations

import re
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from fractions import Fraction
from itertools import combinations, permutations, product
from math import ceil, floor, gcd, sqrt
from pathlib import Path
from typing import Callable, Mapping

from .bond_types import cif_type_to_bond_order, is_aromatic_bond_order
from .chem.rdkit import detect_rdkit_motif_count
from .cofid import (
    COFidMonomer,
    canonicalize_smiles,
    cofid_to_build_request,
    parse_cofid,
    read_cofid_from_cif,
    serialize_cofid,
)
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
            "periodic_rank": _periodic_gain_rank(self.node_count, self.gain_edges),
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
    beta_ketoenamine_single_bond_indices: tuple[int, ...]
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
            "beta_ketoenamine_single_bond_count": len(self.beta_ketoenamine_single_bond_indices),
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
    boron_linkage_rejected_bond_indices: tuple[int, ...]
    recognized_boron_linkages: tuple[str, ...]
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
            "boron_linkage_rejected_bond_count": len(self.boron_linkage_rejected_bond_indices),
            "recognized_boron_linkages": list(self.recognized_boron_linkages),
            "boron_linkage_override_applied": bool(self.boron_linkage_rejected_bond_indices),
            "accepted_bond_count": len(self.accepted_assignments),
            "excluded_ring_sizes": [5, 6],
            "linkage_priority": "candidate-local; coexisting boron motifs do not suppress vinylene",
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
    linkage: str = "auto",
    bond_mode: str = "auto",
    decomposition_mode: str = "event",
) -> CifDecompositionResult:
    input_path = Path(cif_path)
    normalized_decomposition_mode = str(decomposition_mode).strip().lower()
    if normalized_decomposition_mode == "event":
        # Imported lazily so the event pipeline can reuse the stable graph,
        # precursor, and topology primitives in this module without creating an
        # import cycle.  Event mode is the default after its CoRE-COFs labelled
        # benchmark; callers can still request the compatibility path explicitly.
        from .decompose_events import decompose_cif_to_cofid_event

        return decompose_cif_to_cofid_event(
            input_path,
            topology=topology,
            linkage=linkage,
            bond_mode=bond_mode,
        )
    if normalized_decomposition_mode != "legacy":
        return CifDecompositionResult(
            status="error",
            input_cif=str(input_path),
            topology=None if topology is None else str(topology),
            linkage=str(linkage),
            reason=(
                "decomposition_mode must be 'legacy' or 'event'; "
                f"received {decomposition_mode!r}"
            ),
            metadata={"decomposition_mode": normalized_decomposition_mode},
        )
    normalized_linkage = str(linkage).strip().lower()
    if normalized_linkage == "auto":
        return _decompose_cif_with_auto_linkage(
            input_path,
            topology=topology,
            bond_mode=bond_mode,
        )
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

    metadata: dict[str, object] = {}
    normalized_topology: str | None = None
    try:
        if bond_mode not in {"auto", "distance"}:
            raise ValueError("bond_mode must be 'auto' or 'distance'")
        if topology is not None:
            normalized_topology = str(topology).strip().lower()
            if not normalized_topology:
                raise ValueError("topology must not be blank")
            default_topology_repository().get_index_entry(normalized_topology)
        atoms = read_periodic_cif_atoms(input_path)
        _validate_supported_cif_periodicity(atoms)
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
            triazine_resolution = metadata.get("triazine_linkage_resolution")
            if (
                isinstance(triazine_resolution, Mapping)
                and bool(triazine_resolution.get("ambiguous", False))
            ):
                competing = ", ".join(
                    str(linkage)
                    for linkage in triazine_resolution.get("coexisting_linkages", ())
                )
                reason = (
                    "triazine linkage assignment is ambiguous because a complete competing "
                    f"periodic decomposition was recovered ({competing})"
                )
            return CifDecompositionResult(
                status="skipped",
                input_cif=str(input_path),
                topology=normalized_topology,
                linkage=spec.linkage_code,
                reason=reason,
                metadata=metadata,
            )
        recovery_error, recovery_metadata = _validate_recovered_precursors(monomers, spec)
        metadata["precursor_recovery_validation"] = recovery_metadata
        if recovery_error is not None:
            return CifDecompositionResult(
                status="skipped",
                input_cif=str(input_path),
                topology=normalized_topology,
                linkage=spec.linkage_code,
                reason=recovery_error,
                monomers=monomers,
                metadata=metadata,
            )
        selected_topology = normalized_topology
        if selected_topology is not None:
            topology_error, topology_validation = _validate_selected_topology(
                details.topology_graph,
                selected_topology,
                monomers,
                spec,
            )
            metadata["topology_validation"] = topology_validation
            if topology_error is not None:
                return CifDecompositionResult(
                    status="skipped",
                    input_cif=str(input_path),
                    topology=selected_topology,
                    linkage=spec.linkage_code,
                    reason=topology_error,
                    monomers=monomers,
                    metadata=metadata,
                )
        else:
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
        parse_cofid(cofid)
        cofid_to_build_request(cofid)
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
            status="error",
            input_cif=str(input_path),
            topology=normalized_topology if normalized_topology is not None else (
                None if topology is None else str(topology)
            ),
            linkage=spec.linkage_code,
            reason=f"{type(exc).__name__}: {exc}",
            metadata=metadata,
        )


def _decompose_cif_with_auto_linkage(
    input_path: Path,
    *,
    topology: str | None,
    bond_mode: str,
) -> CifDecompositionResult:
    canonical_linkages = (
        "azine",
        "hydrazone",
        "bken",
        "imine",
        "boest",
        "vinylene",
        "boroxine",
        "triazine",
    )
    evaluated = tuple(
        decompose_cif_to_cofid(
            input_path,
            topology=topology,
            linkage=candidate_linkage,
            bond_mode=bond_mode,
            decomposition_mode="legacy",
        )
        for candidate_linkage in canonical_linkages
    )
    successful = tuple(result for result in evaluated if result.ok)
    detection_metadata = {
        "requested_linkage": "auto",
        "successful_linkages": [result.linkage for result in successful],
        "candidates": [
            {
                "linkage": result.linkage,
                "status": result.status,
                "reason": result.reason,
                "cofid": result.cofid,
            }
            for result in evaluated
        ],
    }
    if len(successful) == 1:
        selected = successful[0]
        return CifDecompositionResult(
            status=selected.status,
            input_cif=selected.input_cif,
            topology=selected.topology,
            linkage=selected.linkage,
            cofid=selected.cofid,
            monomers=selected.monomers,
            reason=selected.reason,
            metadata={**dict(selected.metadata), "linkage_detection": detection_metadata},
        )
    normalized_topology = None if topology is None else str(topology).strip().lower() or None
    if len(successful) > 1:
        names = ", ".join(result.linkage for result in successful)
        return CifDecompositionResult(
            status="ambiguous",
            input_cif=str(input_path),
            topology=normalized_topology,
            linkage="auto",
            reason=f"linkage auto-detection is ambiguous among: {names}",
            metadata={"linkage_detection": detection_metadata},
        )
    all_errors = all(result.status == "error" for result in evaluated)
    common_reasons = {result.reason for result in evaluated}
    reason = "no supported linkage produced a complete periodic precursor decomposition"
    if all_errors and len(common_reasons) == 1:
        reason = next(iter(common_reasons)) or reason
    return CifDecompositionResult(
        status="error" if all_errors else "skipped",
        input_cif=str(input_path),
        topology=normalized_topology,
        linkage="auto",
        reason=reason,
        metadata={"linkage_detection": detection_metadata},
    )


def _validate_recovered_precursors(
    monomers: tuple[DecomposedMonomer, ...],
    spec: DecompositionSpec,
) -> tuple[str | None, dict[str, object]]:
    expected_roles = (
        tuple(spec.roles)
        if isinstance(spec, LinkageDecompositionSpec)
        else (spec.reactive_group,)
    )
    actual_roles = tuple(monomer.reactive_group for monomer in monomers)
    expected_counts = Counter(expected_roles)
    actual_counts = Counter(actual_roles)
    metadata: dict[str, object] = {
        "status": "valid",
        "expected_monomer_count": len(expected_roles),
        "actual_monomer_count": len(monomers),
        "expected_role_counts": dict(sorted(expected_counts.items())),
        "actual_role_counts": dict(sorted(actual_counts.items())),
    }
    if len(monomers) != len(expected_roles) or actual_counts != expected_counts:
        metadata["status"] = "incomplete"
        return (
            f"incomplete {spec.linkage_code} precursor recovery: expected "
            f"{len(expected_roles)} monomer block(s) with roles {dict(expected_counts)!r}, "
            f"recovered {len(monomers)} with roles {dict(actual_counts)!r}",
            metadata,
        )
    motif_validation: list[dict[str, object]] = []
    for monomer in monomers:
        validation: dict[str, object] = {
            "reactive_group": monomer.reactive_group,
            "canonical_smiles": monomer.canonical_smiles,
            "declared_connectivity": monomer.connectivity,
        }
        try:
            detected_connectivity = detect_rdkit_motif_count(
                monomer.canonical_smiles,
                monomer.reactive_group,
            )
        except Exception as exc:
            validation.update({
                "status": "invalid",
                "reason": f"{type(exc).__name__}: {exc}",
            })
            motif_validation.append(validation)
            metadata["status"] = "chemically_invalid"
            metadata["motif_validation"] = motif_validation
            return (
                f"recovered {monomer.reactive_group} precursor does not expose a buildable "
                f"{monomer.reactive_group} motif: {type(exc).__name__}: {exc}",
                metadata,
            )
        validation["detected_connectivity"] = detected_connectivity
        validation["status"] = (
            "valid" if detected_connectivity == monomer.connectivity else "connectivity_mismatch"
        )
        motif_validation.append(validation)
        if detected_connectivity != monomer.connectivity:
            metadata["status"] = "chemically_invalid"
            metadata["motif_validation"] = motif_validation
            return (
                f"recovered {monomer.reactive_group} precursor connectivity mismatch: "
                f"decomposition declared {monomer.connectivity}, motif detection found "
                f"{detected_connectivity}",
                metadata,
            )
    metadata["motif_validation"] = motif_validation
    return None, metadata


def _validate_selected_topology(
    graph: LinkageTopologyGraph | None,
    topology: str,
    monomers: tuple[DecomposedMonomer, ...],
    spec: DecompositionSpec,
) -> tuple[str | None, dict[str, object]]:
    repository = default_topology_repository()
    entry = repository.get_index_entry(topology)
    expected_rank = 3 if entry.dimensionality == "3D" else 2
    periodic_rank = 0 if graph is None else _periodic_gain_rank(graph.node_count, graph.gain_edges)
    metadata: dict[str, object] = {
        "status": "valid",
        "topology": topology,
        "topology_dimensionality": entry.dimensionality,
        "expected_periodic_rank": expected_rank,
        "observed_periodic_rank": periodic_rank,
    }
    if graph is None or periodic_rank != expected_rank:
        metadata["status"] = "incompatible"
        return (
            f"recovered linkage graph has periodic rank {periodic_rank}, but topology "
            f"{topology!r} requires rank {expected_rank}",
            metadata,
        )

    connectivities = tuple(sorted((monomer.connectivity for monomer in monomers), reverse=True))
    metadata["recovered_monomer_connectivities"] = list(connectivities)

    compatible = False
    compatibility_source = ""
    if isinstance(spec, RingDecompositionSpec):
        ranked_candidate = next(
            iter(_rank_topology_candidates(graph, topology_ids=frozenset((topology,)))),
            None,
        )
        if ranked_candidate is not None:
            metadata["graph_candidate"] = ranked_candidate.to_dict()
        compatible = ranked_candidate is not None and ranked_candidate.confidence in {"exact", "high"}
        compatibility_source = "periodic_graph"
    elif len(connectivities) == 2:
        high, low = connectivities
        if low == 2 and high >= 3:
            mode = f"{high}+{low}"
            allowed_modes = tuple(str(value) for value in entry.metadata.get("two_monomer_node_linker_modes", ()))
            compatibility_source = "two_monomer_node_linker_modes"
        else:
            # Node-node roles are interchangeable.  Topology metadata stores
            # these modes in ascending order (for example ``3+4``), unlike
            # node-linker modes where the node intentionally comes first.
            mode = f"{low}+{high}"
            allowed_modes = tuple(str(value) for value in entry.metadata.get("two_monomer_node_node_modes", ()))
            compatibility_source = "two_monomer_node_node_modes"
        metadata["connectivity_mode"] = mode
        metadata["allowed_connectivity_modes"] = list(allowed_modes)
        compatible = mode in allowed_modes
        if topology == "bex" and connectivities == (4, 4):
            compatible = True
            compatibility_source = "decorated_bex_4+4"
    metadata["compatibility_source"] = compatibility_source
    if not compatible:
        metadata["status"] = "incompatible"
        return (
            f"recovered precursor connectivities {connectivities!r} and periodic linkage graph "
            f"are incompatible with topology {topology!r}",
            metadata,
        )
    return None, metadata


def _decompose_linkage_atoms(
    atoms: PeriodicCifAtoms,
    spec: LinkageDecompositionSpec,
    *,
    bond_mode: str,
) -> LinkageDecompositionDetails:
    if Chem is None:
        raise RuntimeError("RDKit is required for CIF decomposition.")
    build_result = _build_bonded_mol(atoms, bond_mode=bond_mode)
    return _decompose_linkage_build_result(build_result, spec)


def _decompose_linkage_build_result(
    build_result: BondedMolBuildResult,
    spec: LinkageDecompositionSpec,
) -> LinkageDecompositionDetails:
    mol = build_result.mol
    nitrogen_detection_metadata: dict[str, object] = {}
    nitrogen_classification: NitrogenLinkageBondClassification | None = None
    if spec.linkage_code in {"azine", "hydrazone", "bken", "imine"}:
        nitrogen_classification = _classify_nitrogen_linkage_bonds(mol)
        nitrogen_detection_metadata = {
            "nitrogen_linkage_detection": nitrogen_classification.to_metadata()
        }
    vinylene_detection_metadata: dict[str, object] = {}
    vinylene_classification: VinyleneLinkageBondClassification | None = None
    if spec.linkage_code == "vinylene":
        vinylene_classification = _classify_vinylene_linkage_bonds(mol)
        vinylene_detection_metadata = {
            "vinylene_linkage_detection": vinylene_classification.to_metadata()
        }
    if nitrogen_classification is not None:
        carbon_role, nitrogen_role = {
            "imine": ("aldehyde", "amine"),
            "hydrazone": ("aldehyde", "hydrazide"),
            "azine": ("aldehyde", "hydrazine"),
            "bken": ("keto_aldehyde", "amine"),
        }[spec.linkage_code]
        linkage_bond_indices = _apply_nitrogen_linkage_classification(
            mol,
            classification=nitrogen_classification,
            linkage=spec.linkage_code,
            carbon_role=carbon_role,
            nitrogen_role=nitrogen_role,
        )
    elif vinylene_classification is not None:
        linkage_bond_indices = _apply_vinylene_linkage_classification(mol, vinylene_classification)
    else:
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

    triazine_resolution_metadata: dict[str, object] = {}
    if spec.linkage_code == "triazine":
        recognized_linkages, resolution = _classify_triazine_linkage_priority(
            build_result,
        )
        triazine_resolution_metadata = {"triazine_linkage_resolution": resolution}
        if recognized_linkages:
            return LinkageDecompositionDetails((), {
                **dict(build_result.metadata),
                **triazine_resolution_metadata,
                "decomposition_strategy": "ring",
                "n_atoms": mol.GetNumAtoms(),
                "n_bonds": mol.GetNumBonds(),
                spec.metadata_key: len(ring_events),
                f"n_{spec.linkage_code}_ring_bonds": len(ring_bond_indices),
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
    atom_image_potentials = _fragment_atom_image_potentials(atom_to_fragment, build_result.candidates)
    topology_graph, topology_component_count = _ring_topology_graph(
        ring_events,
        atom_to_fragment,
        atom_image_potentials,
    )
    metadata = {
        **dict(build_result.metadata),
        **triazine_resolution_metadata,
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


def _classify_triazine_linkage_priority(
    build_result: BondedMolBuildResult,
) -> tuple[tuple[str, ...], dict[str, object]]:
    recognized: list[str] = []
    evaluations: dict[str, object] = {}
    for linkage, candidate_spec in (("imine", _IMINE_SPEC), ("vinylene", _VINYLENE_SPEC)):
        try:
            candidate_mol = Chem.Mol(build_result.mol)
            for atom in candidate_mol.GetAtoms():
                if atom.HasProp("cofkit_decompose_role"):
                    atom.ClearProp("cofkit_decompose_role")
            details = _decompose_linkage_build_result(
                BondedMolBuildResult(
                    mol=candidate_mol,
                    metadata=build_result.metadata,
                    candidates=build_result.candidates,
                ),
                candidate_spec,
            )
        except Exception as exc:
            evaluations[linkage] = {
                "status": "error",
                "reason": f"{type(exc).__name__}: {exc}",
            }
            continue
        recovery_error, recovery_metadata = _validate_recovered_precursors(
            details.monomers,
            candidate_spec,
        )
        graph = details.topology_graph
        periodic_rank = (
            0
            if graph is None
            else _periodic_gain_rank(graph.node_count, graph.gain_edges)
        )
        is_recognized = (
            bool(details.monomers)
            and recovery_error is None
            and periodic_rank in {2, 3}
        )
        evaluation: dict[str, object] = {
            "status": "recognized" if is_recognized else "not_recognized",
            "n_linkage_bonds": int(details.metadata.get(candidate_spec.metadata_key, 0)),
            "n_unique_monomers": len(details.monomers),
            "periodic_rank": periodic_rank,
            "precursor_recovery_validation": recovery_metadata,
        }
        if recovery_error is not None:
            evaluation["reason"] = recovery_error
        if linkage == "vinylene":
            vinylene_detection = details.metadata.get("vinylene_linkage_detection")
            if isinstance(vinylene_detection, Mapping):
                evaluation["boron_linkage_override_applied"] = bool(
                    vinylene_detection.get("boron_linkage_override_applied", False)
                )
        evaluations[linkage] = evaluation
        if is_recognized:
            recognized.append(linkage)
    recognized_linkages = tuple(recognized)
    return recognized_linkages, {
        "resolution_strategy": "independent periodic-backbone validation",
        "recognized_higher_priority_linkages": list(recognized_linkages),
        "coexisting_linkages": list(recognized_linkages),
        "ambiguous": bool(recognized_linkages),
        "override_applied": False,
        "evaluations": evaluations,
    }


def _mark_ring_decomposition_events(
    mol,
    candidates: tuple[BondCandidate, ...],
    spec: RingDecompositionSpec,
) -> tuple[tuple[RingDecompositionEvent, ...], tuple[int, ...]]:
    Chem.GetSymmSSSR(mol)
    candidates_by_pair: dict[frozenset[int], list[BondCandidate]] = defaultdict(list)
    for candidate in candidates:
        candidates_by_pair[frozenset((candidate.atom_idx_1, candidate.atom_idx_2))].append(candidate)
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

        atom_images = _unwrap_ring_atom_images(ring, candidates_by_pair)
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
    candidates_by_pair: Mapping[frozenset[int], list[BondCandidate]],
) -> dict[int, tuple[int, int, int]] | None:
    edge_steps: list[tuple[tuple[int, int, int], ...]] = []
    for index, current in enumerate(ring):
        next_idx = ring[(index + 1) % len(ring)]
        candidates = candidates_by_pair.get(frozenset((current, next_idx)), ())
        steps = tuple(
            sorted(
                {
                    _candidate_image_from_to(candidate, current, next_idx)
                    for candidate in candidates
                }
            )
        )
        edge_steps.append(steps or ((0, 0, 0),))

    images: dict[int, tuple[int, int, int]] = {ring[0]: (0, 0, 0)}

    def assign_edge(edge_index: int, current_image: tuple[int, int, int]) -> bool:
        if edge_index == len(ring):
            return current_image == (0, 0, 0)
        next_idx = ring[(edge_index + 1) % len(ring)]
        for step in edge_steps[edge_index]:
            next_image = tuple(current_image[axis] + step[axis] for axis in range(3))
            if edge_index == len(ring) - 1:
                if next_image == (0, 0, 0):
                    return True
                continue
            images[next_idx] = next_image  # type: ignore[assignment]
            if assign_edge(edge_index + 1, next_image):
                return True
            images.pop(next_idx, None)
        return False

    if not assign_edge(0, (0, 0, 0)):
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
    atom_image_potentials: Mapping[int, tuple[int, int, int]] | None,
) -> tuple[LinkageTopologyGraph | None, int]:
    if atom_image_potentials is None:
        return None, 0
    incidences_by_fragment: dict[int, list[tuple[int, tuple[int, int, int]]]] = defaultdict(list)
    for event_index, event in enumerate(events):
        for atom_idx, image in zip(event.anchor_atom_indices, event.anchor_images):
            fragment_idx = atom_to_fragment.get(atom_idx)
            if fragment_idx is not None:
                atom_potential = atom_image_potentials.get(atom_idx, (0, 0, 0))
                fragment_root_image = tuple(
                    image[axis] - atom_potential[axis]
                    for axis in range(3)
                )
                incidences_by_fragment[fragment_idx].append((event_index, fragment_root_image))
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
        dimensionality_hint=_periodic_dimension_hint(len(node_connectivities), tuple(edges)),
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


def _validate_supported_cif_periodicity(atoms: PeriodicCifAtoms) -> None:
    operation_values = (
        tuple(atoms.info.get("_space_group_symop_operation_xyz", ()))
        or tuple(atoms.info.get("_symmetry_equiv_pos_as_xyz", ()))
    )
    normalized_operations = {
        str(value).strip().strip("'\"").lower().replace(" ", "").replace("+", "")
        for value in operation_values
        if str(value).strip().strip("'\"") not in {"", ".", "?"}
    }
    if normalized_operations and normalized_operations != {"x,y,z"}:
        raise ValueError(
            "decomposition currently requires a P1-expanded CIF; non-identity "
            "space-group operations were found"
        )

    number_values = (
        tuple(atoms.info.get("_space_group_IT_number", ()))
        or tuple(atoms.info.get("_symmetry_Int_Tables_number", ()))
    )
    if number_values:
        try:
            space_group_number = int(float(str(number_values[0]).strip().strip("'\"")))
        except ValueError as exc:
            raise ValueError(f"invalid CIF space-group number {number_values[0]!r}") from exc
        if space_group_number != 1:
            raise ValueError(
                "decomposition currently requires a P1-expanded CIF; "
                f"space-group number {space_group_number} was declared"
            )
    name_values = (
        tuple(atoms.info.get("_space_group_name_H-M_alt", ()))
        or tuple(atoms.info.get("_symmetry_space_group_name_H-M", ()))
    )
    if name_values:
        normalized_name = re.sub(
            r"[\s_'\"]+",
            "",
            str(name_values[0]).strip().lower(),
        )
        if normalized_name not in {"p1"}:
            raise ValueError(
                "decomposition currently requires a P1-expanded CIF; "
                f"space group {name_values[0]!r} was declared"
            )


def _build_bonded_mol(atoms: PeriodicCifAtoms, *, bond_mode: str = "auto") -> BondedMolBuildResult:
    labels = tuple(atoms.info.get("_atom_site_label", ()))
    if len(labels) != len(atoms):
        raise ValueError("decomposition requires CIF atom labels aligned with atom sites")

    rw_mol = Chem.RWMol()
    formal_charges = tuple(atoms.info.get("_atom_site_pdbx_formal_charge", ()))
    if formal_charges and len(formal_charges) != len(atoms):
        raise ValueError("CIF formal-charge column length does not match atom labels")
    for atom_index, (label, symbol) in enumerate(zip(labels, atoms.symbols)):
        atom = Chem.Atom(symbol)
        if formal_charges:
            atom.SetFormalCharge(_parse_cif_formal_charge(formal_charges[atom_index]))
        atom.SetProp("cif_label", str(label))
        atom.SetProp("instance_id", _instance_id(str(label)))
        rw_mol.AddAtom(atom)

    candidates, metadata = _bond_candidates_from_cif_or_geometry(atoms, labels, bond_mode=bond_mode)
    if not candidates:
        raise ValueError("decomposition could not infer any covalent bonds from CIF atom positions")

    candidates_by_pair: dict[frozenset[int], list[BondCandidate]] = defaultdict(list)
    for candidate in candidates:
        key = frozenset((candidate.atom_idx_1, candidate.atom_idx_2))
        candidates_by_pair[key].append(candidate)
        if len(candidates_by_pair[key]) == 1:
            rw_mol.AddBond(candidate.atom_idx_1, candidate.atom_idx_2, Chem.BondType.SINGLE)

    mol = rw_mol.GetMol()
    ring_bonds = _ring_bond_keys(mol)
    inferred_order_count = 0
    for pair_candidates in candidates_by_pair.values():
        representative = min(pair_candidates, key=lambda candidate: candidate.distance)
        bond = mol.GetBondBetweenAtoms(representative.atom_idx_1, representative.atom_idx_2)
        if bond is None:
            continue
        explicit_orders = {
            float(candidate.explicit_order)
            for candidate in pair_candidates
            if candidate.explicit_order is not None
        }
        if len(explicit_orders) > 1:
            raise ValueError(
                "periodic CIF bond rows assign conflicting bond orders to the same atom pair"
            )
        order = next(iter(explicit_orders), None)
        if order is None:
            order = _infer_bond_order(mol, representative, ring_bonds)
            inferred_order_count += 1
        _apply_bond_order(mol, bond, order)

    normalized_bken_bond_count = _normalize_beta_ketoenamine_bond_orders(mol)
    mol.UpdatePropertyCache(strict=False)
    return BondedMolBuildResult(
        mol=mol,
        metadata={
            **metadata,
            "n_bond_orders_inferred": inferred_order_count,
            "n_beta_ketoenamine_bond_pairs_normalized": normalized_bken_bond_count,
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
    if bond_mode not in {"auto", "distance"}:
        return TopologyDetectionResult(
            status="failed",
            reason="bond_mode must be 'auto' or 'distance'",
            metadata={"source": "input_validation"},
        )
    comment_cofid = read_cofid_from_cif(input_path)
    parsed_comment = None
    if comment_cofid:
        try:
            parsed_comment = parse_cofid(comment_cofid)
        except ValueError:
            pass

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
            _validate_supported_cif_periodicity(atoms)
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

    periodic_rank = _periodic_gain_rank(graph.node_count, graph.gain_edges)
    if periodic_rank not in {2, 3}:
        return TopologyDetectionResult(
            status="failed",
            reason=(
                f"recovered linkage graph has periodic rank {periodic_rank}; "
                "a COF topology requires rank 2 or 3"
            ),
            metadata={"source": "linkage_graph", "topology_graph": graph.to_metadata()},
        )

    comment_diagnostic: dict[str, object] | None = None
    if parsed_comment is not None:
        comment_spec = _resolve_linkage_spec(parsed_comment.linkage)
        if comment_spec is not None and comment_spec.linkage_code == spec.linkage_code:
            topology_error, topology_validation = _validate_selected_topology(
                graph,
                parsed_comment.topology,
                details.monomers,
                spec,
            )
            if topology_error is None:
                return TopologyDetectionResult(
                    status="success",
                    selected_topology=parsed_comment.topology,
                    confidence="exact",
                    reason="topology recovered from a graph-compatible embedded COFid comment",
                    metadata={
                        "source": "cofid_comment",
                        "cofid": comment_cofid,
                        "topology_validation": topology_validation,
                    },
                )
            comment_diagnostic = {
                "status": "rejected",
                "cofid": comment_cofid,
                "reason": topology_error,
                "topology_validation": topology_validation,
            }
        else:
            comment_diagnostic = {
                "status": "ignored",
                "cofid": comment_cofid,
                "reason": "embedded COFid linkage does not match the requested linkage family",
            }

    candidates = _rank_topology_candidates(graph)
    if not candidates:
        return TopologyDetectionResult(
            status="failed",
            reason="no topology candidates matched the recovered linkage graph connectivities",
            metadata={
                "source": "topology_repository",
                "topology_graph": graph.to_metadata(),
                **({"embedded_cofid_comment": comment_diagnostic} if comment_diagnostic else {}),
            },
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
            metadata={
                "source": "periodic_linkage_graph",
                "topology_graph": graph.to_metadata(),
                **({"embedded_cofid_comment": comment_diagnostic} if comment_diagnostic else {}),
            },
        )

    ambiguous = ", ".join(candidate.topology for candidate in candidates[:5])
    return TopologyDetectionResult(
        status="ambiguous",
        confidence=best.confidence,
        reason=f"topology auto-detection is ambiguous among: {ambiguous}",
        candidates=candidates[:10],
        metadata={
            "source": "periodic_linkage_graph",
            "topology_graph": graph.to_metadata(),
            **({"embedded_cofid_comment": comment_diagnostic} if comment_diagnostic else {}),
        },
    )


def _linkage_topology_graph(
    mol,
    linkage_bond_atom_pairs: list[tuple[int, int]],
    atom_to_fragment: Mapping[int, int],
    candidates: tuple[BondCandidate, ...],
) -> LinkageTopologyGraph | None:
    candidates_by_pair: dict[frozenset[int], list[BondCandidate]] = defaultdict(list)
    for candidate in candidates:
        candidates_by_pair[frozenset((candidate.atom_idx_1, candidate.atom_idx_2))].append(candidate)
    atom_image_potentials = _fragment_atom_image_potentials(atom_to_fragment, candidates)
    if atom_image_potentials is None:
        return None
    participating_fragments: set[int] = set()
    raw_edges: list[tuple[int, int, tuple[int, int, int]]] = []
    for atom_idx_1, atom_idx_2 in linkage_bond_atom_pairs:
        fragment_1 = atom_to_fragment.get(atom_idx_1)
        fragment_2 = atom_to_fragment.get(atom_idx_2)
        if fragment_1 is None or fragment_2 is None or fragment_1 == fragment_2:
            continue
        participating_fragments.update((fragment_1, fragment_2))
        pair_candidates = candidates_by_pair.get(frozenset((atom_idx_1, atom_idx_2)), ())
        if not pair_candidates:
            image = tuple(
                atom_image_potentials.get(atom_idx_1, (0, 0, 0))[axis]
                - atom_image_potentials.get(atom_idx_2, (0, 0, 0))[axis]
                for axis in range(3)
            )
            raw_edges.append((fragment_1, fragment_2, image))
            continue
        for candidate in pair_candidates:
            image = candidate.periodic_image
            if (candidate.atom_idx_1, candidate.atom_idx_2) != (atom_idx_1, atom_idx_2):
                image = _negate_image(image)
            image = tuple(
                atom_image_potentials.get(atom_idx_1, (0, 0, 0))[axis]
                + image[axis]
                - atom_image_potentials.get(atom_idx_2, (0, 0, 0))[axis]
                for axis in range(3)
            )
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

    dimensionality_hint = _periodic_dimension_hint(len(node_map), tuple(edges))
    return LinkageTopologyGraph(
        node_connectivities=tuple(int(degrees[index]) for index in range(len(node_map))),
        gain_edges=tuple(sorted(edges)),
        dimensionality_hint=dimensionality_hint,
    )


def _fragment_atom_image_potentials(
    atom_to_fragment: Mapping[int, int],
    candidates: tuple[BondCandidate, ...],
) -> dict[int, tuple[int, int, int]] | None:
    """Unwrap each cut precursor fragment relative to a stable fragment root.

    CIF atom coordinates may wrap a finite precursor across cell boundaries.  A
    linkage edge gain must therefore include the internal image offsets between
    each endpoint and its fragment root; using only the linkage bond's symmetry
    code incorrectly turns many valid periodic frameworks into rank-zero graphs.
    A conflicting internal cycle means the alleged precursor fragment is itself
    periodic, so it cannot be represented as a finite COFid monomer.
    """

    adjacency: dict[int, list[tuple[int, tuple[int, int, int]]]] = defaultdict(list)
    for candidate in candidates:
        fragment_1 = atom_to_fragment.get(candidate.atom_idx_1)
        fragment_2 = atom_to_fragment.get(candidate.atom_idx_2)
        if fragment_1 is None or fragment_1 != fragment_2:
            continue
        adjacency[candidate.atom_idx_1].append((candidate.atom_idx_2, candidate.periodic_image))
        adjacency[candidate.atom_idx_2].append((candidate.atom_idx_1, _negate_image(candidate.periodic_image)))

    potentials: dict[int, tuple[int, int, int]] = {}
    atoms_by_fragment: dict[int, list[int]] = defaultdict(list)
    for atom_idx, fragment_idx in atom_to_fragment.items():
        atoms_by_fragment[int(fragment_idx)].append(int(atom_idx))
    for fragment_idx in sorted(atoms_by_fragment):
        atoms = sorted(atoms_by_fragment[fragment_idx])
        if not atoms:
            continue
        root = atoms[0]
        potentials[root] = (0, 0, 0)
        stack = [root]
        while stack:
            atom_idx = stack.pop()
            atom_potential = potentials[atom_idx]
            for neighbor_idx, step in adjacency.get(atom_idx, ()):
                expected = tuple(
                    atom_potential[axis] + step[axis]
                    for axis in range(3)
                )
                previous = potentials.get(neighbor_idx)
                if previous is None:
                    potentials[neighbor_idx] = expected  # type: ignore[assignment]
                    stack.append(neighbor_idx)
                elif previous != expected:
                    return None
        for atom_idx in atoms:
            potentials.setdefault(atom_idx, (0, 0, 0))
    return potentials


def _periodic_dimension_hint(
    node_count: int,
    edges: tuple[tuple[int, int, tuple[int, int, int]], ...],
) -> str | None:
    rank = _periodic_gain_rank(node_count, edges)
    if rank == 2:
        return "2D"
    if rank == 3:
        return "3D"
    return None


def _periodic_gain_rank(
    node_count: int,
    edges: tuple[tuple[int, int, tuple[int, int, int]], ...],
) -> int:
    if node_count <= 0 or not edges:
        return 0
    adjacency: dict[int, list[tuple[int, tuple[int, int, int]]]] = defaultdict(list)
    for start, end, image in edges:
        adjacency[start].append((end, image))
        adjacency[end].append((start, _negate_image(image)))

    potentials: dict[int, tuple[int, int, int]] = {}
    for seed in range(node_count):
        if seed in potentials:
            continue
        potentials[seed] = (0, 0, 0)
        stack = [seed]
        while stack:
            node = stack.pop()
            node_potential = potentials[node]
            for neighbor, step in adjacency.get(node, ()):
                if neighbor in potentials:
                    continue
                potentials[neighbor] = tuple(
                    node_potential[axis] + step[axis]
                    for axis in range(3)
                )  # type: ignore[assignment]
                stack.append(neighbor)

    cycle_gains = []
    for start, end, image in edges:
        start_potential = potentials[start]
        end_potential = potentials[end]
        cycle_gains.append(
            tuple(
                start_potential[axis] + image[axis] - end_potential[axis]
                for axis in range(3)
            )
        )
    return _integer_vector_rank(tuple(cycle_gains))


def _integer_vector_rank(vectors: tuple[tuple[int, int, int], ...]) -> int:
    nonzero = [vector for vector in vectors if vector != (0, 0, 0)]
    if not nonzero:
        return 0
    first = nonzero[0]
    second = next((vector for vector in nonzero[1:] if _cross(first, vector) != (0, 0, 0)), None)
    if second is None:
        return 1
    normal = _cross(first, second)
    if any(_dot_int(normal, vector) != 0 for vector in nonzero):
        return 3
    return 2


def _cross(
    first: tuple[int, int, int],
    second: tuple[int, int, int],
) -> tuple[int, int, int]:
    return (
        first[1] * second[2] - first[2] * second[1],
        first[2] * second[0] - first[0] * second[2],
        first[0] * second[1] - first[1] * second[0],
    )


def _dot_int(first: tuple[int, int, int], second: tuple[int, int, int]) -> int:
    return first[0] * second[0] + first[1] * second[1] + first[2] * second[2]


def _rank_topology_candidates(
    graph: LinkageTopologyGraph,
    *,
    topology_ids: frozenset[str] | None = None,
) -> tuple[TopologyDetectionCandidate, ...]:
    repository = default_topology_repository()
    entries = repository.list_index(dimensionality=graph.dimensionality_hint)
    if not entries:
        entries = repository.list_index()
    if topology_ids is not None:
        entries = tuple(entry for entry in entries if entry.id in topology_ids)

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
            if not compare_gains:
                return _edge_counter(tuple(mapped_edges), compare_gains=False) == expected_edge_counter
            return _periodic_edge_gains_match(
                tuple(mapped_edges),
                expected.gain_edges,
            )

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


def _periodic_edge_gains_match(
    observed_edges: tuple[tuple[int, int, tuple[int, int, int]], ...],
    expected_edges: tuple[tuple[int, int, tuple[int, int, int]], ...],
) -> bool:
    if _pair_edge_counter(observed_edges) != _pair_edge_counter(expected_edges):
        return False
    tested_transforms: set[tuple[tuple[int, int, int], ...]] = set()
    for axis_order in permutations((0, 1, 2)):
        for signs in product((-1, 1), repeat=3):
            transform = tuple(
                tuple(
                    signs[row] if axis_order[row] == column else 0
                    for column in range(3)
                )
                for row in range(3)
            )
            tested_transforms.add(transform)  # type: ignore[arg-type]
            transformed = _transform_periodic_edges(observed_edges, transform)  # type: ignore[arg-type]
            if _edge_gains_match_up_to_vertex_switching(transformed, expected_edges):
                return True
    for transform in _candidate_unimodular_gain_transforms(observed_edges, expected_edges):
        if transform in tested_transforms:
            continue
        tested_transforms.add(transform)
        transformed = _transform_periodic_edges(observed_edges, transform)
        if _edge_gains_match_up_to_vertex_switching(transformed, expected_edges):
            return True
    return False


def _transform_periodic_edges(
    edges: tuple[tuple[int, int, tuple[int, int, int]], ...],
    transform: tuple[tuple[int, int, int], ...],
) -> tuple[tuple[int, int, tuple[int, int, int]], ...]:
    return tuple(
        (
            start,
            end,
            tuple(
                sum(transform[row][column] * image[column] for column in range(3))
                for row in range(3)
            ),
        )
        for start, end, image in edges
    )


def _candidate_unimodular_gain_transforms(
    observed_edges: tuple[tuple[int, int, tuple[int, int, int]], ...],
    expected_edges: tuple[tuple[int, int, tuple[int, int, int]], ...],
) -> tuple[tuple[tuple[int, int, int], ...], ...]:
    observed_vectors = _gain_transform_invariant_vectors(observed_edges)
    expected_vectors = _gain_transform_invariant_vectors(expected_edges)
    observed_rank = _integer_vector_rank(observed_vectors)
    expected_rank = _integer_vector_rank(expected_vectors)
    if observed_rank != expected_rank or observed_rank < 2:
        return ()

    rank = observed_rank
    observed_bases = _independent_gain_bases(observed_vectors, rank, limit=32)
    expected_bases = _independent_gain_bases(expected_vectors, rank, limit=32)
    transforms: set[tuple[tuple[int, int, int], ...]] = set()
    for observed_basis in observed_bases:
        for expected_basis in expected_bases:
            for ordered_expected in permutations(expected_basis):
                transform = _unimodular_transform_between_bases(observed_basis, ordered_expected)
                if transform is not None:
                    transforms.add(transform)
    return tuple(sorted(transforms))


def _gain_transform_invariant_vectors(
    edges: tuple[tuple[int, int, tuple[int, int, int]], ...],
) -> tuple[tuple[int, int, int], ...]:
    parallel_differences = _parallel_edge_gain_differences(edges)
    if _integer_vector_rank(parallel_differences) >= 2:
        return parallel_differences
    pair_counts = _pair_edge_counter(edges)
    if any(count > 1 for count in pair_counts.values()):
        return parallel_differences
    cycle_gains = _simple_graph_cycle_gain_vectors(edges)
    return tuple(sorted(set(parallel_differences).union(cycle_gains)))


def _parallel_edge_gain_differences(
    edges: tuple[tuple[int, int, tuple[int, int, int]], ...],
) -> tuple[tuple[int, int, int], ...]:
    gains_by_pair: dict[tuple[int, int], set[tuple[int, int, int]]] = defaultdict(set)
    for start, end, image in edges:
        gains_by_pair[(start, end)].add(image)
    differences: set[tuple[int, int, int]] = set()
    for gains in gains_by_pair.values():
        for first, second in combinations(sorted(gains), 2):
            difference = tuple(first[axis] - second[axis] for axis in range(3))
            if difference != (0, 0, 0):
                differences.add(difference)  # type: ignore[arg-type]
                differences.add(_negate_image(difference))  # type: ignore[arg-type]
    return tuple(sorted(differences))


def _simple_graph_cycle_gain_vectors(
    edges: tuple[tuple[int, int, tuple[int, int, int]], ...],
) -> tuple[tuple[int, int, int], ...]:
    if not edges:
        return ()
    adjacency: dict[int, list[tuple[int, int, tuple[int, int, int]]]] = defaultdict(list)
    for edge_index, (start, end, image) in enumerate(edges):
        adjacency[start].append((end, edge_index, image))
        adjacency[end].append((start, edge_index, _negate_image(image)))
    potentials: dict[int, tuple[int, int, int]] = {}
    tree_edges: set[int] = set()
    for seed in sorted(adjacency):
        if seed in potentials:
            continue
        potentials[seed] = (0, 0, 0)
        stack = [seed]
        while stack:
            node = stack.pop()
            for neighbor, edge_index, step in sorted(adjacency[node]):
                if neighbor in potentials:
                    continue
                node_potential = potentials[node]
                potentials[neighbor] = tuple(
                    node_potential[axis] + step[axis]
                    for axis in range(3)
                )  # type: ignore[assignment]
                tree_edges.add(edge_index)
                stack.append(neighbor)

    cycle_gains: set[tuple[int, int, int]] = set()
    for edge_index, (start, end, image) in enumerate(edges):
        if edge_index in tree_edges:
            continue
        gain = tuple(
            potentials[start][axis] + image[axis] - potentials[end][axis]
            for axis in range(3)
        )
        if gain != (0, 0, 0):
            cycle_gains.add(gain)  # type: ignore[arg-type]
            cycle_gains.add(_negate_image(gain))  # type: ignore[arg-type]
    return tuple(sorted(cycle_gains))


def _independent_gain_bases(
    vectors: tuple[tuple[int, int, int], ...],
    rank: int,
    *,
    limit: int,
) -> tuple[tuple[tuple[int, int, int], ...], ...]:
    bases: list[tuple[tuple[int, int, int], ...]] = []
    for basis in combinations(vectors, rank):
        if _integer_vector_rank(tuple(basis)) != rank:
            continue
        bases.append(tuple(basis))
        if len(bases) >= limit:
            break
    return tuple(bases)


def _unimodular_transform_between_bases(
    observed_basis: tuple[tuple[int, int, int], ...],
    expected_basis: tuple[tuple[int, int, int], ...],
) -> tuple[tuple[int, int, int], ...] | None:
    if len(observed_basis) != len(expected_basis) or len(observed_basis) not in {2, 3}:
        return None
    observed_columns = observed_basis
    expected_columns = expected_basis
    if len(observed_basis) == 2:
        observed_completion = _complete_rank_two_integer_basis(observed_basis[0], observed_basis[1])
        expected_completion = _complete_rank_two_integer_basis(expected_basis[0], expected_basis[1])
        if observed_completion is None or expected_completion is None:
            return None
        observed_columns = (*observed_basis, observed_completion)
        expected_columns = (*expected_basis, expected_completion)

    observed_matrix = _matrix_from_columns(observed_columns)  # type: ignore[arg-type]
    expected_matrix = _matrix_from_columns(expected_columns)  # type: ignore[arg-type]
    inverse_observed = _inverse_integer_matrix_as_fractions(observed_matrix)
    if inverse_observed is None:
        return None
    product_matrix = tuple(
        tuple(
            sum(
                Fraction(expected_matrix[row][index]) * inverse_observed[index][column]
                for index in range(3)
            )
            for column in range(3)
        )
        for row in range(3)
    )
    if any(value.denominator != 1 for row in product_matrix for value in row):
        return None
    transform = tuple(
        tuple(int(value) for value in row)
        for row in product_matrix
    )
    if abs(_determinant_3x3_int(transform)) != 1:
        return None
    return transform  # type: ignore[return-value]


def _complete_rank_two_integer_basis(
    first: tuple[int, int, int],
    second: tuple[int, int, int],
) -> tuple[int, int, int] | None:
    normal = _cross(first, second)
    divisor = gcd(gcd(abs(normal[0]), abs(normal[1])), abs(normal[2]))
    if divisor == 0:
        return None
    coefficients = _bezout_coefficients_three(*normal)
    if _dot_int(normal, coefficients) != divisor:
        return None
    return coefficients


def _bezout_coefficients_three(first: int, second: int, third: int) -> tuple[int, int, int]:
    gcd_first_second, first_coefficient, second_coefficient = _extended_gcd(first, second)
    common_gcd, pair_coefficient, third_coefficient = _extended_gcd(gcd_first_second, third)
    if common_gcd == 0:
        return (0, 0, 0)
    return (
        first_coefficient * pair_coefficient,
        second_coefficient * pair_coefficient,
        third_coefficient,
    )


def _extended_gcd(first: int, second: int) -> tuple[int, int, int]:
    if second == 0:
        sign = 1 if first >= 0 else -1
        return abs(first), sign, 0
    divisor, remainder = divmod(first, second)
    common_gcd, coefficient_second, coefficient_remainder = _extended_gcd(second, remainder)
    return (
        common_gcd,
        coefficient_remainder,
        coefficient_second - divisor * coefficient_remainder,
    )


def _matrix_from_columns(
    columns: tuple[tuple[int, int, int], tuple[int, int, int], tuple[int, int, int]],
) -> tuple[tuple[int, int, int], ...]:
    return tuple(
        tuple(columns[column][row] for column in range(3))
        for row in range(3)
    )


def _inverse_integer_matrix_as_fractions(
    matrix: tuple[tuple[int, int, int], ...],
) -> tuple[tuple[Fraction, Fraction, Fraction], ...] | None:
    (a, b, c), (d, e, f), (g, h, i) = matrix
    determinant = _determinant_3x3_int(matrix)
    if determinant == 0:
        return None
    return (
        (
            Fraction(e * i - f * h, determinant),
            Fraction(c * h - b * i, determinant),
            Fraction(b * f - c * e, determinant),
        ),
        (
            Fraction(f * g - d * i, determinant),
            Fraction(a * i - c * g, determinant),
            Fraction(c * d - a * f, determinant),
        ),
        (
            Fraction(d * h - e * g, determinant),
            Fraction(b * g - a * h, determinant),
            Fraction(a * e - b * d, determinant),
        ),
    )


def _determinant_3x3_int(matrix: tuple[tuple[int, int, int], ...]) -> int:
    (a, b, c), (d, e, f), (g, h, i) = matrix
    return a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)


def _edge_gains_match_up_to_vertex_switching(
    observed_edges: tuple[tuple[int, int, tuple[int, int, int]], ...],
    expected_edges: tuple[tuple[int, int, tuple[int, int, int]], ...],
) -> bool:
    observed_by_pair: dict[tuple[int, int], Counter[tuple[int, int, int]]] = defaultdict(Counter)
    expected_by_pair: dict[tuple[int, int], Counter[tuple[int, int, int]]] = defaultdict(Counter)
    for start, end, image in observed_edges:
        observed_by_pair[(start, end)][image] += 1
    for start, end, image in expected_edges:
        expected_by_pair[(start, end)][image] += 1
    if set(observed_by_pair) != set(expected_by_pair):
        return False

    constraints: list[tuple[int, int, tuple[tuple[int, int, int], ...]]] = []
    for (start, end), observed_counter in observed_by_pair.items():
        expected_counter = expected_by_pair[(start, end)]
        if sum(observed_counter.values()) != sum(expected_counter.values()):
            return False
        if start == end:
            if observed_counter != expected_counter:
                return False
            continue
        observed_seed = next(iter(observed_counter))
        differences: set[tuple[int, int, int]] = set()
        for expected_image in expected_counter:
            difference = tuple(expected_image[axis] - observed_seed[axis] for axis in range(3))
            shifted = Counter(
                {
                    tuple(image[axis] + difference[axis] for axis in range(3)): count
                    for image, count in observed_counter.items()
                }
            )
            if shifted == expected_counter:
                differences.add(difference)  # type: ignore[arg-type]
        if not differences:
            return False
        constraints.append((start, end, tuple(sorted(differences))))

    constraints.sort(key=lambda item: (len(item[2]), item[0], item[1]))
    offsets: dict[int, tuple[int, int, int]] = {}

    def solve(position: int) -> bool:
        if position == len(constraints):
            return True
        start, end, differences = constraints[position]
        start_offset = offsets.get(start)
        end_offset = offsets.get(end)
        for difference in differences:
            assigned: list[int] = []
            if start_offset is not None and end_offset is not None:
                actual = tuple(start_offset[axis] - end_offset[axis] for axis in range(3))
                if actual != difference:
                    continue
            elif start_offset is not None:
                offsets[end] = tuple(start_offset[axis] - difference[axis] for axis in range(3))  # type: ignore[assignment]
                assigned.append(end)
            elif end_offset is not None:
                offsets[start] = tuple(end_offset[axis] + difference[axis] for axis in range(3))  # type: ignore[assignment]
                assigned.append(start)
            else:
                offsets[start] = (0, 0, 0)
                offsets[end] = _negate_image(difference)
                assigned.extend((start, end))
            if solve(position + 1):
                return True
            for node in assigned:
                offsets.pop(node, None)
        return False

    return solve(0)


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
    if len(bond_labels_1) != len(bond_labels_2):
        raise ValueError("CIF explicit-bond endpoint columns have different lengths")
    bond_types = tuple(atoms.info.get("_ccdc_geom_bond_type") or atoms.info.get("_geom_bond_type") or ())
    bond_distances = tuple(atoms.info.get("_geom_bond_distance") or ())
    symmetries_1 = tuple(atoms.info.get("_geom_bond_site_symmetry_1") or ())
    symmetries_2 = tuple(atoms.info.get("_geom_bond_site_symmetry_2") or ())
    row_count = len(bond_labels_1)
    for column_name, values in (
        ("bond type", bond_types),
        ("bond distance", bond_distances),
        ("first bond symmetry", symmetries_1),
        ("second bond symmetry", symmetries_2),
    ):
        if values and len(values) != row_count:
            raise ValueError(
                f"CIF explicit {column_name} column has {len(values)} rows; expected {row_count}"
            )
    candidates: list[BondCandidate] = []
    added_edges: set[tuple[int, int, tuple[int, int, int]]] = set()
    for row_idx, (label_1, label_2) in enumerate(zip(bond_labels_1, bond_labels_2)):
        if label_1 not in label_to_idx or label_2 not in label_to_idx:
            continue
        idx_1 = label_to_idx[label_1]
        idx_2 = label_to_idx[label_2]
        symmetry_1 = symmetries_1[row_idx] if row_idx < len(symmetries_1) else None
        symmetry_2 = symmetries_2[row_idx] if row_idx < len(symmetries_2) else None
        if _has_cif_symmetry_hint(symmetry_1) or _has_cif_symmetry_hint(symmetry_2):
            image = _explicit_bond_periodic_image(symmetries_1, symmetries_2, row_idx)
        else:
            image, _distance = _minimum_image_bond_geometry(atoms, idx_1, idx_2)
        if idx_1 == idx_2:
            if image != (0, 0, 0):
                raise ValueError(
                    "self-image CIF bonds require a P1 supercell expansion before decomposition"
                )
            continue
        if idx_2 < idx_1:
            idx_1, idx_2 = idx_2, idx_1
            image = _negate_image(image)
        edge_key = (idx_1, idx_2, image)
        if edge_key in added_edges:
            continue
        order = cif_type_to_bond_order(bond_types[row_idx] if row_idx < len(bond_types) else None)
        distance = _explicit_bond_distance(atoms, idx_1, idx_2, bond_distances, row_idx)
        candidates.append(BondCandidate(idx_1, idx_2, distance, explicit_order=order, periodic_image=image))
        added_edges.add(edge_key)
    return tuple(candidates), sum(1 for candidate in candidates if candidate.explicit_order is None)


def _explicit_bond_periodic_image(
    symmetries_1: tuple[str, ...],
    symmetries_2: tuple[str, ...],
    row_idx: int,
) -> tuple[int, int, int]:
    first = _parse_p1_symmetry_shift(symmetries_1[row_idx] if row_idx < len(symmetries_1) else None)
    second = _parse_p1_symmetry_shift(symmetries_2[row_idx] if row_idx < len(symmetries_2) else None)
    return (second[0] - first[0], second[1] - first[1], second[2] - first[2])


def _has_cif_symmetry_hint(value: object | None) -> bool:
    text = "" if value is None else str(value).strip().strip("'\"")
    return text not in {"", ".", "?", "1"}


def _parse_p1_symmetry_shift(value: object | None) -> tuple[int, int, int]:
    text = "" if value is None else str(value).strip().strip("'\"")
    if not text or text in {".", "?", "1"}:
        return (0, 0, 0)
    match = re.match(r"^1_([0-9])([0-9])([0-9])$", text)
    if match is None:
        raise ValueError(
            f"unsupported CIF bond symmetry code {value!r}; a P1-expanded CIF is required"
        )
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
    reciprocal_row_norms = _fractional_from_cartesian_row_norms(atoms.cell_basis)
    for idx_1, symbol_1 in enumerate(atoms.symbols):
        for idx_2 in range(idx_1 + 1, len(atoms.symbols)):
            symbol_2 = atoms.symbols[idx_2]
            image, distance = _minimum_image_bond_geometry(
                atoms,
                idx_1,
                idx_2,
                reciprocal_row_norms=reciprocal_row_norms,
            )
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


def _normalize_beta_ketoenamine_bond_orders(mol) -> int:
    """Canonicalize enol-imine-like input to the keto-enamine C-N/C=C form."""

    normalized = 0
    for cn_bond in tuple(mol.GetBonds()):
        if float(cn_bond.GetBondTypeAsDouble()) not in {1.0, 2.0}:
            continue
        first = cn_bond.GetBeginAtom()
        second = cn_bond.GetEndAtom()
        if {first.GetAtomicNum(), second.GetAtomicNum()} != {6, 7}:
            continue
        first_instance = first.GetProp("instance_id")
        second_instance = second.GetProp("instance_id")
        if first_instance and second_instance and first_instance == second_instance:
            continue
        carbon = first if first.GetAtomicNum() == 6 else second
        nitrogen = first if first.GetAtomicNum() == 7 else second
        for ring_anchor in carbon.GetNeighbors():
            if ring_anchor.GetIdx() == nitrogen.GetIdx() or ring_anchor.GetAtomicNum() != 6:
                continue
            anchor_bond = mol.GetBondBetweenAtoms(carbon.GetIdx(), ring_anchor.GetIdx())
            if anchor_bond is None or float(anchor_bond.GetBondTypeAsDouble()) not in {1.0, 2.0}:
                continue
            if not ring_anchor.IsInRing():
                continue
            if not any(
                carbonyl.GetIdx() != carbon.GetIdx()
                and carbonyl.GetAtomicNum() == 6
                and mol.GetBondBetweenAtoms(ring_anchor.GetIdx(), carbonyl.GetIdx()) is not None
                and mol.GetBondBetweenAtoms(ring_anchor.GetIdx(), carbonyl.GetIdx()).IsInRing()
                and carbonyl.IsInRing()
                and _atom_has_double_bonded_oxygen(
                    mol,
                    carbonyl.GetIdx(),
                    exclude_atom_idx=ring_anchor.GetIdx(),
                )
                for carbonyl in ring_anchor.GetNeighbors()
            ):
                continue
            changed = False
            if abs(float(cn_bond.GetBondTypeAsDouble()) - 1.0) > 1.0e-6:
                cn_bond.SetIsAromatic(False)
                cn_bond.SetBondType(Chem.BondType.SINGLE)
                changed = True
            if abs(float(anchor_bond.GetBondTypeAsDouble()) - 2.0) > 1.0e-6:
                anchor_bond.SetIsAromatic(False)
                anchor_bond.SetBondType(Chem.BondType.DOUBLE)
                changed = True
            if changed:
                normalized += 1
            break
    return normalized


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


def _eligible_beta_ketoenamine_single_bonds(mol) -> dict[int, tuple[int, int]]:
    eligible: dict[int, tuple[int, int]] = {}
    for bond in mol.GetBonds():
        if abs(float(bond.GetBondTypeAsDouble()) - 1.0) > 1.0e-6:
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
        if _carbon_has_beta_ketoenamine_environment(
            mol,
            carbon.GetIdx(),
            exclude_atom_idx=nitrogen.GetIdx(),
            anchor_bond_order=2.0,
        ):
            eligible[bond.GetIdx()] = (carbon.GetIdx(), nitrogen.GetIdx())
    return eligible


def _classify_nitrogen_linkage_bonds(mol) -> NitrogenLinkageBondClassification:
    eligible = _eligible_carbon_nitrogen_double_bonds(mol)
    eligible_bken_single = _eligible_beta_ketoenamine_single_bonds(mol)
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

        if _carbon_has_beta_ketoenamine_environment(
            mol,
            carbon_idx,
            exclude_atom_idx=nitrogen_idx,
            anchor_bond_order=1.0,
        ):
            bken_candidates.add(bond_idx)
    bken_candidates.update(eligible_bken_single)

    assigned_azine: set[int] = set()
    assigned_hydrazone: set[int] = set()
    assigned_bken: set[int] = set()
    assigned_imine: set[int] = set()
    cross_branch_conflicts: set[int] = set()
    for bond_idx in set(eligible) | set(eligible_bken_single):
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
        raw_cn_bond_indices=tuple(sorted(set(eligible) | set(eligible_bken_single))),
        beta_ketoenamine_single_bond_indices=tuple(sorted(eligible_bken_single)),
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


def _carbon_has_beta_ketoenamine_environment(
    mol,
    carbon_idx: int,
    *,
    exclude_atom_idx: int,
    anchor_bond_order: float,
) -> bool:
    carbon = mol.GetAtomWithIdx(carbon_idx)
    for ring_anchor in carbon.GetNeighbors():
        if ring_anchor.GetIdx() == exclude_atom_idx or ring_anchor.GetAtomicNum() != 6:
            continue
        anchor_bond = mol.GetBondBetweenAtoms(carbon_idx, ring_anchor.GetIdx())
        if anchor_bond is None or abs(float(anchor_bond.GetBondTypeAsDouble()) - anchor_bond_order) > 1.0e-6:
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
    return _apply_nitrogen_linkage_classification(
        mol,
        classification=classification,
        linkage=linkage,
        carbon_role=carbon_role,
        nitrogen_role=nitrogen_role,
    )


def _apply_nitrogen_linkage_classification(
    mol,
    *,
    classification: NitrogenLinkageBondClassification,
    linkage: str,
    carbon_role: str,
    nitrogen_role: str,
) -> tuple[int, ...]:
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
    return _apply_vinylene_linkage_classification(mol, classification)


def _apply_vinylene_linkage_classification(
    mol,
    classification: VinyleneLinkageBondClassification,
) -> tuple[int, ...]:
    for _bond_idx, aldehyde_idx, activated_idx in classification.accepted_assignments:
        mol.GetAtomWithIdx(aldehyde_idx).SetProp("cofkit_decompose_role", "aldehyde")
        mol.GetAtomWithIdx(activated_idx).SetProp("cofkit_decompose_role", "activated_methylene")
    return classification.accepted_bond_indices


def _classify_vinylene_linkage_bonds(mol) -> VinyleneLinkageBondClassification:
    formal_cc_double_bonds: list[int] = []
    same_instance_rejected: list[int] = []
    small_ring_rejected: list[int] = []
    endpoint_rejected: list[int] = []
    structurally_accepted: list[tuple[int, int, int]] = []
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
        structurally_accepted.append((bond_idx, aldehyde_idx, activated_idx))

    recognized_boron_linkages = _recognized_boron_linkages(mol)
    boron_linkage_rejected: tuple[int, ...] = ()
    accepted = tuple(sorted(structurally_accepted))
    return VinyleneLinkageBondClassification(
        formal_cc_double_bond_indices=tuple(sorted(formal_cc_double_bonds)),
        same_instance_rejected_bond_indices=tuple(sorted(same_instance_rejected)),
        small_ring_rejected_bond_indices=tuple(sorted(small_ring_rejected)),
        endpoint_rejected_bond_indices=tuple(sorted(endpoint_rejected)),
        boron_linkage_rejected_bond_indices=boron_linkage_rejected,
        recognized_boron_linkages=recognized_boron_linkages,
        accepted_assignments=accepted,
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


def _boroxine_ring_bond_indices(mol) -> tuple[int, ...]:
    Chem.GetSymmSSSR(mol)
    bond_indices: set[int] = set()
    for raw_ring in mol.GetRingInfo().AtomRings():
        ring = tuple(int(atom_idx) for atom_idx in raw_ring)
        if not _matches_ring_decomposition_pattern(mol, ring, _BOROXINE_SPEC):
            continue
        for index, atom_idx in enumerate(ring):
            next_idx = ring[(index + 1) % len(ring)]
            bond = mol.GetBondBetweenAtoms(atom_idx, next_idx)
            if bond is not None:
                bond_indices.add(int(bond.GetIdx()))
    return tuple(sorted(bond_indices))


def _boronate_ester_linkage_bond_indices(mol) -> tuple[int, ...]:
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
        bond_indices.append(bond.GetIdx())
    return tuple(bond_indices)


def _recognized_boron_linkages(mol) -> tuple[str, ...]:
    boroxine_bonds = set(_boroxine_ring_bond_indices(mol))
    boronate_ester_bonds = set(_boronate_ester_linkage_bond_indices(mol)) - boroxine_bonds
    recognized: list[str] = []
    if boronate_ester_bonds:
        recognized.append("boest")
    if boroxine_bonds:
        recognized.append("boroxine")
    return tuple(recognized)


def _mark_boronate_ester_linkage_bonds(mol) -> tuple[int, ...]:
    bond_indices = _boronate_ester_linkage_bond_indices(mol)
    for bond_idx in bond_indices:
        bond = mol.GetBondWithIdx(bond_idx)
        first = bond.GetBeginAtom()
        second = bond.GetEndAtom()
        boron = first if first.GetAtomicNum() == 5 else second
        oxygen = first if first.GetAtomicNum() == 8 else second
        boron.SetProp("cofkit_decompose_role", "boronic_acid")
        oxygen.SetProp("cofkit_decompose_role", "catechol")
    return bond_indices


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


def _restore_keto_aldehyde_fragment(fragment, *, endpoint_role: str = "keto_aldehyde"):
    editable = Chem.RWMol(fragment)
    endpoint_indices = [
        atom.GetIdx()
        for atom in editable.GetAtoms()
        if atom.HasProp("cofkit_decompose_role")
        and atom.GetProp("cofkit_decompose_role") == endpoint_role
    ]
    rings_to_rearomatize: set[frozenset[int]] = set()
    Chem.GetSymmSSSR(editable)
    atom_rings = tuple(tuple(int(atom_idx) for atom_idx in ring) for ring in editable.GetRingInfo().AtomRings())
    for carbon_idx in endpoint_indices:
        carbon = editable.GetAtomWithIdx(carbon_idx)
        for neighbor in carbon.GetNeighbors():
            if neighbor.GetAtomicNum() != 6:
                continue
            for ring in atom_rings:
                if (
                    neighbor.GetIdx() in ring
                    and len(ring) == 6
                    and all(editable.GetAtomWithIdx(atom_idx).GetAtomicNum() == 6 for atom_idx in ring)
                ):
                    rings_to_rearomatize.add(frozenset(ring))
            bond = editable.GetBondBetweenAtoms(carbon_idx, neighbor.GetIdx())
            if bond is not None and abs(float(bond.GetBondTypeAsDouble()) - 2.0) <= 1.0e-6:
                bond.SetIsAromatic(False)
                bond.SetBondType(Chem.BondType.SINGLE)
    for bond in tuple(editable.GetBonds()):
        if bond.GetBondType() != Chem.BondType.DOUBLE:
            continue
        first = bond.GetBeginAtom()
        second = bond.GetEndAtom()
        if {first.GetAtomicNum(), second.GetAtomicNum()} != {6, 8}:
            continue
        carbon = first if first.GetAtomicNum() == 6 else second
        if carbon.GetIsAromatic() or carbon.IsInRing():
            bond.SetBondType(Chem.BondType.SINGLE)
    for ring in rings_to_rearomatize:
        ring_atoms = tuple(sorted(ring))
        for atom_idx in ring_atoms:
            editable.GetAtomWithIdx(atom_idx).SetIsAromatic(True)
        for bond in editable.GetBonds():
            if bond.GetBeginAtomIdx() in ring and bond.GetEndAtomIdx() in ring:
                bond.SetBondType(Chem.BondType.AROMATIC)
                bond.SetIsAromatic(True)
    return _restore_double_oxygens(editable.GetMol(), endpoint_role)


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
        if atom.HasProp("cofkit_decompose_role"):
            atom.ClearProp("cofkit_decompose_role")
        if atom.HasProp("cif_label"):
            atom.ClearProp("cif_label")
        if atom.HasProp("instance_id"):
            atom.ClearProp("instance_id")
    # CIFs commonly carry explicit hydrogens on linkage atoms.  Bond cutting
    # and precursor-group restoration change the heavy-atom valence first, so
    # sanitizing while those product-state hydrogens are still attached can
    # raise a spurious valence error (for example N(4), O(3), or C(5)).  Remove
    # them without sanitization and let RDKit assign precursor-state implicit
    # hydrogens during the single final sanitization.
    # Remove malformed bridging hydrogens explicitly.  RDKit intentionally
    # preserves degree-two H atoms, but in a distance-bonded CIF graph these
    # are almost always duplicated/disordered proton positions rather than
    # covalent framework bridges.
    editable = Chem.RWMol(mol)
    for atom_idx in sorted(
        (
            atom.GetIdx()
            for atom in editable.GetAtoms()
            if atom.GetAtomicNum() == 1 and atom.GetDegree() != 1
        ),
        reverse=True,
    ):
        editable.RemoveAtom(atom_idx)
    mol = Chem.RemoveHs(editable.GetMol(), sanitize=False)
    # CIF formal-charge fields are frequently absent even when the normalized
    # bond graph contains a conventional tetravalent iminium/pyridinium N.
    # Restoring the charge is deterministic from that local valence and avoids
    # rejecting an otherwise valid precursor for N(4).
    for atom in mol.GetAtoms():
        explicit_valence = sum(float(bond.GetBondTypeAsDouble()) for bond in atom.GetBonds())
        if (
            atom.GetAtomicNum() == 7
            and atom.GetFormalCharge() == 0
            and 3.9 <= explicit_valence <= 4.1
        ):
            atom.SetFormalCharge(1)
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
    *,
    reciprocal_row_norms: tuple[float, float, float] | None = None,
) -> tuple[tuple[int, int, int], float]:
    frac_1 = atoms.fractional_positions[idx_1]
    frac_2 = atoms.fractional_positions[idx_2]
    delta = tuple(frac_2[axis] - frac_1[axis] for axis in range(3))
    initial_image = tuple(-int(round(value)) for value in delta)
    initial_delta = tuple(delta[axis] + initial_image[axis] for axis in range(3))
    initial_cartesian = _fractional_delta_to_cartesian(atoms.cell_basis, initial_delta)  # type: ignore[arg-type]
    best_image = initial_image
    best_distance_squared = sum(component * component for component in initial_cartesian)

    reciprocal_row_norms = reciprocal_row_norms or _fractional_from_cartesian_row_norms(atoms.cell_basis)
    best_distance = sqrt(best_distance_squared)
    image_ranges: list[range] = []
    for axis in range(3):
        radius = reciprocal_row_norms[axis] * best_distance + 1.0e-10
        lower = ceil(-delta[axis] - radius)
        upper = floor(-delta[axis] + radius)
        lower = min(lower, initial_image[axis])
        upper = max(upper, initial_image[axis])
        image_ranges.append(range(lower, upper + 1))

    for image in product(*image_ranges):
        candidate_delta = tuple(delta[axis] + image[axis] for axis in range(3))
        cartesian = _fractional_delta_to_cartesian(atoms.cell_basis, candidate_delta)  # type: ignore[arg-type]
        distance_squared = sum(component * component for component in cartesian)
        if distance_squared + 1.0e-12 < best_distance_squared:
            best_image = image  # type: ignore[assignment]
            best_distance_squared = distance_squared
    return best_image, sqrt(best_distance_squared)


def _fractional_from_cartesian_row_norms(
    cell_basis: tuple[Vec3, Vec3, Vec3],
) -> tuple[float, float, float]:
    cartesian_from_fractional = (
        (cell_basis[0][0], cell_basis[1][0], cell_basis[2][0]),
        (cell_basis[0][1], cell_basis[1][1], cell_basis[2][1]),
        (cell_basis[0][2], cell_basis[1][2], cell_basis[2][2]),
    )
    inverse = _inverse_3x3(cartesian_from_fractional)
    return tuple(sqrt(sum(component * component for component in row)) for row in inverse)  # type: ignore[return-value]


def _inverse_3x3(matrix: tuple[Vec3, Vec3, Vec3]) -> tuple[Vec3, Vec3, Vec3]:
    (a, b, c), (d, e, f), (g, h, i) = matrix
    determinant = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)
    if abs(determinant) <= 1.0e-12:
        raise ValueError("CIF cell basis is singular")
    inverse_determinant = 1.0 / determinant
    return (
        (
            (e * i - f * h) * inverse_determinant,
            (c * h - b * i) * inverse_determinant,
            (b * f - c * e) * inverse_determinant,
        ),
        (
            (f * g - d * i) * inverse_determinant,
            (a * i - c * g) * inverse_determinant,
            (c * d - a * f) * inverse_determinant,
        ),
        (
            (d * h - e * g) * inverse_determinant,
            (b * g - a * h) * inverse_determinant,
            (a * e - b * d) * inverse_determinant,
        ),
    )


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


def _parse_cif_formal_charge(value: object) -> int:
    text = str(value).strip().strip("'\"")
    if text in {"", ".", "?"}:
        return 0
    trailing_sign = re.fullmatch(r"(\d+)([+-])", text)
    if trailing_sign is not None:
        magnitude = int(trailing_sign.group(1))
        return magnitude if trailing_sign.group(2) == "+" else -magnitude
    try:
        charge = float(text)
    except ValueError as exc:
        raise ValueError(f"invalid CIF formal charge {value!r}") from exc
    rounded = round(charge)
    if abs(charge - rounded) > 1.0e-8:
        raise ValueError(f"CIF formal charge must be integral, got {value!r}")
    return int(rounded)


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
