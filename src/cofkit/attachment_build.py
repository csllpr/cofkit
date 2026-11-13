"""Build a fixed periodic attachment assignment using the existing realizers."""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Mapping

import gemmi
import networkx as nx
from rdkit import Chem

from .attachment import (
    AttachmentInput,
    AttachmentBond,
    compare_attachment_graphs,
    load_attachment_recipe,
)
from .cofid import canonicalize_smiles, cofid_to_build_request
from .geometry import fit_rigid_transform
from .model import (
    AssemblyState,
    Candidate,
    MonomerInstance,
    MonomerSpec,
    MotifRef,
    Pose,
    ReactionEvent,
)
from .periodic_geometry import p1_shift
from .planner import AssignmentPlan, NetPlan
from .product_graph import PeriodicProductGraph
from .reactions import ReactionLibrary
from .scoring import CandidateScorer
from .search import AssignmentOutcome
from .topologies import get_topology_hint


@dataclass(frozen=True)
class AttachmentBuildResult:
    candidate: Candidate
    monomers: tuple[MonomerSpec, ...]
    graph: PeriodicProductGraph
    outcome: AssignmentOutcome


def _atom_mapping(canonical_smiles, spec):
    canonical = Chem.MolFromSmiles(canonical_smiles)
    source = nx.Graph()
    for atom in canonical.GetAtoms():
        source.add_node(atom.GetIdx(), element=atom.GetSymbol())
    source.add_edges_from(
        (b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in canonical.GetBonds()
    )
    target = nx.Graph()
    for i, symbol in enumerate(spec.atom_symbols):
        if symbol != "H":
            target.add_node(i, element=symbol)
    target.add_edges_from(
        (a, b) for a, b, _ in spec.bonds if a in target and b in target
    )
    matcher = nx.algorithms.isomorphism.GraphMatcher(
        source, target, node_match=lambda a, b: a["element"] == b["element"]
    )
    mapping = next(matcher.isomorphisms_iter(), None)
    if mapping is None:
        raise ValueError(
            f"monomer {spec.id!r} does not have the attachment precursor's heavy graph"
        )
    return mapping


def _bind_monomers(recipe, monomers):
    request = cofid_to_build_request(recipe.cofid)
    if monomers is None:
        from .chem.rdkit import build_rdkit_monomer

        monomers = tuple(
            build_rdkit_monomer(
                f"attachment_m{i}", f"attachment_m{i}", m.canonical_smiles, m.motif_kind
            )
            for i, m in enumerate(request.monomers)
        )
    specs = (
        tuple(monomers.values()) if isinstance(monomers, Mapping) else tuple(monomers)
    )
    bound = []
    maps = []
    for monomer in request.monomers:
        choices = []
        for spec in specs:
            smiles = spec.metadata.get("canonical_smiles") or spec.metadata.get(
                "source_smiles"
            )
            if (
                isinstance(smiles, str)
                and canonicalize_smiles(smiles) == monomer.canonical_smiles
                and set(spec.motif_kinds) == {monomer.motif_kind}
            ):
                choices.append(spec)
        if len(choices) != 1:
            raise ValueError(
                "attachment recipe requires exactly one matching atomistic monomer for each COFid block"
            )
        spec = choices[0]
        if len(spec.motifs) != monomer.connectivity:
            raise ValueError("attachment monomer motif count disagrees with the COFid")
        bound.append(spec)
        maps.append(_atom_mapping(monomer.canonical_smiles, spec))
    if len({s.id for s in bound}) != len(bound):
        raise ValueError("attachment precursor ids must be distinct")
    return request, tuple(bound), maps


def _resolve_site(site, spec, mapping):
    atoms = {mapping[a] for a in site.atom_indices}
    choices = []
    for motif in spec.motifs:
        primary = motif.metadata.get("reactive_atom_id")
        reactive = motif.metadata.get("reactive_atom_ids")
        if len(atoms) == 1:
            matches = primary in atoms
        else:
            matches = bool(reactive) and atoms == set(reactive)
        if matches:
            choices.append(motif)
    if len(choices) != 1:
        raise ValueError(
            f"attachment site {site.instance_id}:{site.atom_indices} does not identify exactly one reactive motif"
        )
    return choices[0]


def build_from_attachment(
    attachment: AttachmentInput,
    monomers=None,
    *,
    candidate_id: str = "attachment-candidate-1",
    initial_state: AssemblyState | None = None,
) -> AttachmentBuildResult:
    """Honor a recipe before fitting poses; no site assignment is optimized away.

    Source fragment coordinates are optional rigid-placement hints, not a
    relaxation result. Exact connectivity verification is independent of pose
    quality and ignores bond orders, hydrogens, charges and stereochemistry.
    """
    recipe = load_attachment_recipe(attachment)
    request, specs, maps = _bind_monomers(recipe, monomers)
    by_instance = {i.id: i for i in recipe.instances}
    spec_map = {s.id: s for s in specs}
    instance_specs = {i.id: specs[i.monomer_index] for i in recipe.instances}
    poses = {}
    fit_errors = {}
    for instance in recipe.instances:
        spec = instance_specs[instance.id]
        if not instance.seed_positions:
            if initial_state is not None and instance.id in initial_state.monomer_poses:
                if initial_state.cell != recipe.cell:
                    raise ValueError(
                        "initial attachment state must use the recipe cell"
                    )
                poses[instance.id] = initial_state.monomer_poses[instance.id]
                continue
            raise ValueError(
                f"attachment instance {instance.id!r} has no placement hints; supply seed_positions or an initial_state"
            )
        mapping = maps[instance.monomer_index]
        reference = [p.position for p in instance.seed_positions]
        probe = [
            spec.atom_positions[mapping[p.atom_index]] for p in instance.seed_positions
        ]
        fit = fit_rigid_transform(tuple(reference), tuple(probe))
        poses[instance.id] = Pose(fit.translation, fit.rotation)
        fit_errors[instance.id] = fit.rms

    graph = PeriodicProductGraph()
    instances = tuple(
        MonomerInstance(i.id, instance_specs[i.id].id) for i in recipe.instances
    )
    for instance in instances:
        graph.add_monomer(instance)
    template = ReactionLibrary.builtin().get(request.template_id)
    events = []
    for raw in recipe.events:
        refs = []
        metadata = dict(raw.metadata)
        for site in raw.participants:
            instance = by_instance[site.instance_id]
            spec = instance_specs[site.instance_id]
            motif = _resolve_site(site, spec, maps[instance.monomer_index])
            refs.append(
                MotifRef(site.instance_id, spec.id, motif.id, site.periodic_image)
            )
        if "attachment_retained_oxygen_ids" in metadata:
            ids = metadata["attachment_retained_oxygen_ids"]
            if len(ids) != len(refs):
                raise ValueError(
                    "ring attachment oxygen choices do not match the participants"
                )
            metadata["attachment_retained_oxygen_ids"] = tuple(
                maps[by_instance[s.instance_id].monomer_index][a]
                for s, a in zip(raw.participants, ids)
            )
        event = ReactionEvent(
            raw.id,
            raw.template_id,
            tuple(refs),
            metadata={**metadata, "attachment_source_event_id": raw.source_event_id},
        )
        graph.add_reaction_event(event, spec_map, {template.id: template})
        events.append(event)
    unreacted = graph.unreacted_motifs(spec_map)
    if unreacted:
        raise ValueError("attachment recipe leaves precursor motifs unreacted")
    plan = AssignmentPlan(
        NetPlan(
            get_topology_hint(request.topology),
            tuple(s.id for s in specs),
            (template.id,),
        ),
        {i.id: i.monomer_id for i in instances},
    )
    outcome = AssignmentOutcome(
        plan, instances, tuple(events), (), sum(len(e.participants) for e in events)
    )
    state = AssemblyState(
        cell=recipe.cell, monomer_poses=poses, stacking_state="attachment_cell"
    )
    flags = [
        "attachment_constrained",
        "source_placement_hints",
        "stacking_defined_by_attachment",
    ]
    metadata = {
        "instance_to_monomer": {i.id: i.monomer_id for i in instances},
        "graph_summary": graph.summary(),
        "net_plan": {
            "topology": request.topology,
            "metadata": {"planning_mode": "attachment_recipe"},
        },
        "assignment": dict(plan.slot_to_monomer),
        "embedding": {
            "method": "attachment_source_rigid_fit",
            "physical_relaxation": False,
            "fit_rms_angstrom": fit_errors,
        },
        "optimization": {
            "enabled": False,
            "reason": "attachment connectivity is fixed; source fragment fits are seed placement",
        },
        "attachment": {
            "status": "applied",
            "recipe": recipe.to_dict(),
            "recipe_fingerprint": recipe.fingerprint,
            "provenance": dict(recipe.provenance),
        },
    }
    c_axis_semantics = recipe.provenance.get("c_axis_semantics")
    if c_axis_semantics is not None:
        metadata["embedding"]["c_axis_semantics"] = str(c_axis_semantics)
    candidate = Candidate(
        candidate_id, None, state, tuple(events), tuple(flags), metadata
    )
    if template.topology_role == "ring":
        from .ring_geometry import ring_geometry_report, validate_ring_geometry
        from .build_workflows.ring_forming import (
            _ring_validation_flags,
            _ring_validation_metadata,
        )

        geometry = ring_geometry_report(tuple(events), state, spec_map)
        validation = validate_ring_geometry(tuple(events), state, spec_map)
        flags.extend(_ring_validation_flags(validation))
        metadata.update(
            ring_validation=_ring_validation_metadata(validation),
            score_metadata={
                "ring_geometry": geometry.as_dict(),
                "n_unreacted_motifs": 0,
            },
        )
        candidate = replace(candidate, flags=tuple(flags), metadata=metadata)
    else:
        metadata["score_metadata"] = CandidateScorer().scoring_metadata(
            outcome, state, spec_map, {template.id: template}
        )
    verification = verify_attachment_candidate(recipe, candidate, specs, maps=maps)
    if verification["status"] == "mismatch":
        raise ValueError(
            f"attachment realization does not preserve its reference graph: {verification['reason']}"
        )
    metadata["attachment"]["verification"] = verification
    return AttachmentBuildResult(
        replace(candidate, metadata=metadata), specs, graph, outcome
    )


def verify_attachment_candidate(
    recipe, candidate, monomers, *, maps=None, cif_text=None
):
    """Audit the actual CIF writer's output, including retained product atoms."""
    from .cif import CIFWriter
    from .reaction_realization import ReactionRealizer

    recipe = load_attachment_recipe(recipe)
    request, specs, bound_maps = _bind_monomers(recipe, monomers)
    maps = bound_maps if maps is None else maps
    if not recipe.reference_atoms:
        return {
            "status": "not_checked",
            "scope": "periodic_heavy_graph",
            "reason": "recipe supplied no source reference graph",
        }
    label_map = {}
    for instance in recipe.instances:
        spec = specs[instance.monomer_index]
        for canonical, atom in maps[instance.monomer_index].items():
            label_map[
                ReactionRealizer.atom_label(instance.id, spec.atom_symbols[atom], atom)
            ] = (instance.id, canonical)
    text = (
        CIFWriter().export_candidate(candidate, specs, cofid=request.cofid).text
        if cif_text is None
        else cif_text
    )
    block = gemmi.cif.read_string(text).sole_block()
    actual_atoms = {}
    for label, symbol in block.find(["_atom_site_label", "_atom_site_type_symbol"]):
        if str(symbol) == "H":
            continue
        if str(label) not in label_map:
            return {
                "status": "mismatch",
                "reason": "realization produced an unmapped heavy atom",
                "scope": "periodic_heavy_graph",
            }
        actual_atoms[label_map[str(label)]] = str(symbol)
    actual_bonds = []
    for a, b, s1, s2 in block.find(
        [
            "_geom_bond_atom_site_label_1",
            "_geom_bond_atom_site_label_2",
            "_geom_bond_site_symmetry_1",
            "_geom_bond_site_symmetry_2",
        ]
    ):
        if str(a) not in label_map or str(b) not in label_map:
            continue
        actual_bonds.append(
            AttachmentBond(
                label_map[str(a)],
                label_map[str(b)],
                tuple(y - x for x, y in zip(p1_shift(str(s1)), p1_shift(str(s2)))),
            )
        )
    return compare_attachment_graphs(
        {a.key: a.element for a in recipe.reference_atoms},
        recipe.reference_bonds,
        actual_atoms,
        actual_bonds,
    )


def verify_attachment_cif(recipe, cif_path, monomers=None):
    """Check a written CIF against a recipe's known atom correspondence."""
    from pathlib import Path

    return verify_attachment_candidate(
        recipe, None, monomers, cif_text=Path(cif_path).read_text()
    )
