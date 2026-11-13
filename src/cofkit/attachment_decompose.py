"""Recover attachment recipes while event fragments still have atom identities."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import replace
import hashlib
import json
from pathlib import Path

import networkx as nx
from rdkit import Chem, rdBase

from .attachment import (
    AttachmentAtom,
    AttachmentBond,
    AttachmentEvent,
    AttachmentInstance,
    AttachmentPosition,
    AttachmentRecipe,
    AttachmentSite,
    canonical_attachment_bond,
)
from .cofid import cofid_to_build_request
from ._version import __version__
from .geometry import cross, norm, normalize, sub


def _remap(source_smiles, target_smiles):
    """A chosen full molecular correspondence, never independent port ranks."""
    source, target = (
        Chem.MolFromSmiles(source_smiles),
        Chem.MolFromSmiles(target_smiles),
    )
    if source_smiles == target_smiles:
        return dict(enumerate(range(source.GetNumAtoms())))

    def graph(mol):
        g = nx.Graph()
        for atom in mol.GetAtoms():
            g.add_node(atom.GetIdx(), element=atom.GetSymbol())
        g.add_edges_from(
            (b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds()
        )
        return g

    matcher = nx.algorithms.isomorphism.GraphMatcher(
        graph(source),
        graph(target),
        node_match=lambda a, b: a["element"] == b["element"],
    )
    mapping = next(matcher.isomorphisms_iter(), None)
    if mapping is None:
        raise ValueError(
            "recovered precursor has no element-graph correspondence to the selected COFid monomer"
        )
    return mapping


def extract_attachment_recipe(
    input_path, atoms, build_result, hypothesis
) -> AttachmentRecipe:
    from . import decompose as legacy
    from .cif import read_c_axis_semantics_from_cif
    from .decompose_events import _precursor_identity_report

    request = cofid_to_build_request(hypothesis.cofid)
    roles = {r.fragment_id: r for r in hypothesis.roles}
    instance_ids = {f: f"f{f}" for f in roles}
    prototype_indices = {}
    source_map = {}
    for fragment, role in roles.items():
        choices = [
            (i, m)
            for i, m in enumerate(request.monomers)
            if m.reactive_group == role.role
        ]
        if len(choices) != 1:
            raise ValueError(
                "attachment precursor role is not uniquely represented by the COFid"
            )
        index, monomer = choices[0]
        if not role.source_atom_map:
            raise ValueError(
                role.source_atom_map_error
                or "precursor repair did not retain a complete canonical atom correspondence"
            )
        remapping = _remap(role.canonical_smiles, monomer.canonical_smiles)
        prototype_indices[fragment] = index
        for source_atom, canonical in role.source_atom_map.items():
            source_map[source_atom] = (instance_ids[fragment], remapping[canonical])
    cut_pairs = {
        frozenset(
            (
                build_result.mol.GetBondWithIdx(i).GetBeginAtomIdx(),
                build_result.mol.GetBondWithIdx(i).GetEndAtomIdx(),
            )
        )
        for event in hypothesis.events
        for i in event.cut_bonds
    }
    # In a higher-connected ring precursor, cut product bonds can join
    # periodic copies of the same fragment. They are not internal precursor
    # edges and must not be reintroduced during unwrapping.
    internal_candidates = tuple(
        b
        for b in build_result.candidates
        if frozenset((b.atom_idx_1, b.atom_idx_2)) not in cut_pairs
    )
    potentials = legacy._fragment_atom_image_potentials(
        hypothesis.fragment_by_atom, internal_candidates
    )
    if potentials is None:
        raise ValueError("precursor fragments cannot be unwrapped consistently")
    potentials = dict(potentials)
    candidates_by_pair = defaultdict(list)
    for bond in build_result.candidates:
        candidates_by_pair[frozenset((bond.atom_idx_1, bond.atom_idx_2))].append(bond)

    def image(a, b):
        candidates = candidates_by_pair[frozenset((a, b))]
        if len(candidates) != 1:
            raise ValueError(
                "attachment endpoints need an explicit supercell to separate parallel-image bonds"
            )
        candidate = candidates[0]
        return (
            candidate.periodic_image
            if candidate.atom_idx_1 == a
            else tuple(-x for x in candidate.periodic_image)
        )

    def cart(frac):
        return tuple(
            sum(frac[j] * atoms.cell_basis[j][i] for j in range(3)) for i in range(3)
        )

    def site(source_atoms, periodic_image=(0, 0, 0)):
        refs = [source_map[a] for a in source_atoms]
        if len({a[0] for a in refs}) != 1:
            raise ValueError("one reactive site spans several precursor fragments")
        return AttachmentSite(
            refs[0][0], tuple(sorted({a[1] for a in refs})), periodic_image
        )

    attachment_events = []
    for event in hypothesis.events:
        if event.family in {"boroxine", "triazine"}:
            ring = tuple(int(a) for a in event.metadata["ring_atoms"])
            ring_images = legacy._unwrap_ring_atom_images(ring, candidates_by_pair)
            if ring_images is None:
                raise ValueError("ring event has inconsistent periodic images")
            anchors = tuple(int(a) for a in event.metadata["anchor_atom_indices"])
            if event.family == "triazine":
                # Follow C -> retained N -> next C; the two alternating
                # nitrile interpretations require opposite ring directions.
                retained = {
                    frozenset(
                        (
                            build_result.mol.GetBondWithIdx(i).GetBeginAtomIdx(),
                            build_result.mol.GetBondWithIdx(i).GetEndAtomIdx(),
                        )
                    )
                    for i in event.metadata["retained_bonds"]
                }
                ordered = [anchors[0]]
                for _ in range(len(anchors) - 1):
                    carbon = ordered[-1]
                    nitrogen = next(
                        a
                        for pair in retained
                        if carbon in pair
                        for a in pair
                        if a != carbon
                    )
                    ordered.append(
                        next(
                            a.GetIdx()
                            for a in build_result.mol.GetAtomWithIdx(
                                nitrogen
                            ).GetNeighbors()
                            if a.GetIdx() in anchors and a.GetIdx() != carbon
                        )
                    )
                anchors = tuple(ordered)
            metadata = {}
            if event.family == "boroxine":
                retained_oxygens = []
                for i, boron in enumerate(anchors):
                    next_boron = anchors[(i + 1) % len(anchors)]
                    oxygen = next(
                        a
                        for a in ring
                        if build_result.mol.GetAtomWithIdx(a).GetSymbol() == "O"
                        and build_result.mol.GetBondBetweenAtoms(a, boron) is not None
                        and build_result.mol.GetBondBetweenAtoms(a, next_boron)
                        is not None
                    )
                    iid, canonical_boron = source_map[boron]
                    fragment = hypothesis.fragment_by_atom[boron]
                    precursor = Chem.MolFromSmiles(
                        request.monomers[prototype_indices[fragment]].canonical_smiles
                    )
                    canonical_oxygen = min(
                        a.GetIdx()
                        for a in precursor.GetAtomWithIdx(
                            canonical_boron
                        ).GetNeighbors()
                        if a.GetSymbol() == "O"
                    )
                    source_map[oxygen] = (iid, canonical_oxygen)
                    potentials[oxygen] = tuple(
                        x + y for x, y in zip(potentials[boron], image(boron, oxygen))
                    )
                    retained_oxygens.append(canonical_oxygen)
                metadata["attachment_retained_oxygen_ids"] = retained_oxygens
            points = [
                cart(
                    tuple(
                        f + g
                        for f, g in zip(atoms.fractional_positions[a], ring_images[a])
                    )
                )
                for a in anchors
            ]
            center = tuple(
                sum(atoms.fractional_positions[a][i] + ring_images[a][i] for a in ring)
                / len(ring)
                for i in range(3)
            )
            normal = cross(sub(points[1], points[0]), sub(points[2], points[0]))
            if norm(normal) == 0:
                raise ValueError("ring event has a degenerate source plane")
            metadata.update(
                ring_center_fractional=center, ring_normal=normalize(normal)
            )
            participants = tuple(
                site((a,), tuple(g - p for g, p in zip(ring_images[a], potentials[a])))
                for a in anchors
            )
            attachment_events.append(
                AttachmentEvent(
                    event.event_id,
                    request.template_id,
                    participants,
                    event.event_id,
                    metadata,
                )
            )
            continue

        cut_pairs = [
            (
                build_result.mol.GetBondWithIdx(i).GetBeginAtomIdx(),
                build_result.mol.GetBondWithIdx(i).GetEndAtomIdx(),
            )
            for i in event.cut_bonds
        ]
        groups = (
            [[pair] for pair in cut_pairs] if event.family == "azine" else [cut_pairs]
        )
        for group_index, pairs in enumerate(groups):
            endpoints = defaultdict(set)
            for a, b in pairs:
                endpoints[hypothesis.fragment_by_atom[a]].add(a)
                endpoints[hypothesis.fragment_by_atom[b]].add(b)
            if len(endpoints) != 2:
                raise ValueError(
                    "binary attachment event does not join two precursor copies"
                )
            left, right = sorted(endpoints)
            gains = set()
            for a, b in pairs:
                if hypothesis.fragment_by_atom[a] != left:
                    a, b = b, a
                gains.add(
                    tuple(
                        x + g - y
                        for x, g, y in zip(potentials[a], image(a, b), potentials[b])
                    )
                )
            if len(gains) != 1:
                raise ValueError(
                    "one binary event joins inconsistent periodic precursor copies"
                )
            participants = (
                site(endpoints[left]),
                site(endpoints[right], next(iter(gains))),
            )
            attachment_events.append(
                AttachmentEvent(
                    f"{event.event_id}:{group_index}",
                    request.template_id,
                    participants,
                    event.event_id,
                )
            )

    # Every event-bearing product heavy atom must be mapped. Disconnected
    # guests remain in the existing atom ledger, outside the framework recipe.
    framework_atoms = {
        a
        for a, f in hypothesis.fragment_by_atom.items()
        if f in roles and build_result.mol.GetAtomWithIdx(a).GetAtomicNum() != 1
    }
    if not framework_atoms <= source_map.keys():
        raise ValueError("precursor repair lost source framework atom identities")
    if len(set(source_map.values())) != len(source_map):
        raise ValueError("source heavy atoms map to a repeated precursor atom")
    labels = atoms.info["_atom_site_label"]
    positions = defaultdict(list)
    reference_atoms = []
    for source, (iid, canonical) in sorted(source_map.items()):
        if build_result.mol.GetAtomWithIdx(source).GetAtomicNum() == 1:
            continue
        pos = cart(
            tuple(
                f + g
                for f, g in zip(atoms.fractional_positions[source], potentials[source])
            )
        )
        positions[iid].append(AttachmentPosition(canonical, pos, labels[source]))
        reference_atoms.append(AttachmentAtom(iid, canonical, atoms.symbols[source]))
    bonds = set()
    for b in build_result.candidates:
        a, z = b.atom_idx_1, b.atom_idx_2
        if a in source_map and z in source_map:
            gain = tuple(
                x + g - y
                for x, g, y in zip(potentials[a], b.periodic_image, potentials[z])
            )
            bonds.add(canonical_attachment_bond(source_map[a], source_map[z], gain))
    recipe = AttachmentRecipe(
        cofid=request.cofid,
        cell=atoms.cell_basis,
        instances=tuple(
            AttachmentInstance(
                instance_ids[f],
                prototype_indices[f],
                tuple(sorted(positions[instance_ids[f]], key=lambda p: p.atom_index)),
            )
            for f in sorted(roles)
        ),
        events=tuple(attachment_events),
        reference_atoms=tuple(reference_atoms),
        reference_bonds=tuple(AttachmentBond(*b) for b in sorted(bonds)),
        provenance={
            "source": "decomposition",
            "cofkit_version": __version__,
            "rdkit_version": rdBase.rdkitVersion,
            "canonical_atom_index_basis": "atom order in parsed COFid canonical precursor SMILES",
            "precursor_identity": _precursor_identity_report(hypothesis.metadata),
            "input_name": Path(input_path).name,
            "source_sha256": hashlib.sha256(Path(input_path).read_bytes()).hexdigest(),
            "bond_source": build_result.metadata.get("bond_source"),
            "source_graph_sha256": hashlib.sha256(
                json.dumps(
                    {
                        "atoms": sorted((a, atoms.symbols[a]) for a in source_map),
                        "bonds": sorted(
                            (b.atom_idx_1, b.atom_idx_2, b.periodic_image)
                            for b in build_result.candidates
                            if b.atom_idx_1 in source_map and b.atom_idx_2 in source_map
                        ),
                    },
                    sort_keys=True,
                ).encode()
            ).hexdigest(),
            "experimental_status": "not_asserted",
            "identity_scope": "periodic_heavy_graph",
            "component_policy": "all_recovered_framework_instances; guests excluded",
            "geometry": "source_fragment_positions_are_placement_hints",
        },
    )
    c_axis_semantics = read_c_axis_semantics_from_cif(input_path)
    provenance = dict(recipe.provenance)
    if c_axis_semantics is not None:
        provenance["c_axis_semantics"] = c_axis_semantics
    if hypothesis.atom_ledger is not None:
        provenance.update(
            excluded_guest_atom_count=hypothesis.atom_ledger.guest_atom_count,
            input_charge_scope=hypothesis.atom_ledger.charge_scope.status,
            input_stereo_scope=hypothesis.atom_ledger.stereo_scope.status,
        )
    recipe = replace(recipe, provenance=provenance)
    recipe.validate()
    return recipe
