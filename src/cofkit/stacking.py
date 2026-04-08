from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, replace
from typing import Mapping

from .geometry import Vec3, add, norm, normalize, scale
from .model import AssemblyState, Candidate, Pose, ReactionEvent
from .topologies import get_topology_hint


_STACKING_LAYER_SUFFIXES: tuple[str, str] = ("L0", "L1")


@dataclass(frozen=True)
class LayerRegistry:
    id: str
    lateral_shift: tuple[float, float] = (0.0, 0.0)
    interlayer_distance: float = 3.4


class StackingExplorer:
    """Enumerates built-in bilayer registries for exported 2D COF candidates.

    `lateral_shift` is expressed as a fractional shift in the in-plane `a`/`b`
    basis and is applied to the second layer in an `A/B` repeating bilayer cell.
    """

    def __init__(self, *, cell_kind: str | None = None) -> None:
        self.cell_kind = str(cell_kind or "").strip().lower() or None

    def enumerate_registries(self) -> tuple[LayerRegistry, ...]:
        if self.cell_kind == "hexagonal":
            return (
                LayerRegistry(id="AA", lateral_shift=(0.0, 0.0), interlayer_distance=3.4),
                LayerRegistry(id="AB", lateral_shift=(1.0 / 3.0, 1.0 / 3.0), interlayer_distance=3.5),
                LayerRegistry(id="slipped", lateral_shift=(0.5, 0.0), interlayer_distance=3.6),
            )
        return (
            LayerRegistry(id="AA", lateral_shift=(0.0, 0.0), interlayer_distance=3.4),
            LayerRegistry(id="AB", lateral_shift=(0.5, 0.5), interlayer_distance=3.5),
            LayerRegistry(id="slipped", lateral_shift=(0.5, 0.0), interlayer_distance=3.6),
        )

    def resolve_registries(self, registry_ids: tuple[str, ...] = ()) -> tuple[LayerRegistry, ...]:
        available = {registry.id.casefold(): registry for registry in self.enumerate_registries()}
        if not registry_ids:
            return tuple(available.values())

        resolved: list[LayerRegistry] = []
        missing: list[str] = []
        for raw_registry_id in registry_ids:
            registry = available.get(str(raw_registry_id).casefold())
            if registry is None:
                missing.append(str(raw_registry_id))
                continue
            if registry.id not in {item.id for item in resolved}:
                resolved.append(registry)
        if missing:
            supported = ", ".join(sorted(registry.id for registry in available.values()))
            raise ValueError(
                f"Unknown stacking registry ids {tuple(missing)!r}. Supported stacking registries: {supported}."
            )
        return tuple(resolved)


def enumerate_candidate_stackings(
    candidate: Candidate,
    *,
    registry_ids: tuple[str, ...] = (),
) -> tuple[Candidate, ...]:
    if not registry_ids or not _is_eligible_2d_candidate(candidate):
        return (candidate,)
    explorer = StackingExplorer(cell_kind=_candidate_cell_kind(candidate))
    registries = explorer.resolve_registries(registry_ids)
    return tuple(_apply_layer_registry(candidate, registry) for registry in registries)


def stacking_comment_suffix(registry: LayerRegistry) -> str:
    return f"stacking={registry.id}"


def _apply_layer_registry(candidate: Candidate, registry: LayerRegistry) -> Candidate:
    state = candidate.state
    if len(_STACKING_LAYER_SUFFIXES) != 2:
        raise ValueError("the current stacking enumerator supports exactly two layers")
    if registry.interlayer_distance <= 0.0:
        raise ValueError("interlayer_distance must be positive")

    base_cell = state.cell
    c_direction = _safe_normalize(base_cell[2])
    new_c = scale(c_direction, 2.0 * registry.interlayer_distance)
    layer_offsets = (
        scale(new_c, 0.25),
        add(scale(new_c, 0.75), _fractional_shift_in_plane(base_cell, registry.lateral_shift)),
    )

    instance_to_monomer = _string_mapping(candidate.metadata.get("instance_to_monomer"))
    instance_to_slot = _string_mapping(candidate.metadata.get("instance_to_slot"))
    monomer_poses: dict[str, Pose] = {}
    stacked_layer_offsets: dict[str, tuple[float, float, float]] = {}
    stacked_instance_to_monomer: dict[str, str] = {}
    stacked_instance_to_slot: dict[str, str] = {}
    base_pose_details = _mapping(candidate.metadata.get("embedding")).get("poses")
    stacked_pose_details: dict[str, object] = {}

    for instance_id, pose in state.monomer_poses.items():
        pose_detail = _mapping(base_pose_details).get(instance_id) if isinstance(base_pose_details, Mapping) else None
        for layer_index, layer_suffix in enumerate(_STACKING_LAYER_SUFFIXES):
            stacked_instance_id = f"{instance_id}{layer_suffix}"
            offset = layer_offsets[layer_index]
            translation = add(pose.translation, offset)
            monomer_poses[stacked_instance_id] = Pose(
                translation=translation,
                rotation_matrix=pose.rotation_matrix,
            )
            stacked_layer_offsets[stacked_instance_id] = offset
            if instance_id in instance_to_monomer:
                stacked_instance_to_monomer[stacked_instance_id] = instance_to_monomer[instance_id]
            if instance_id in instance_to_slot:
                stacked_instance_to_slot[stacked_instance_id] = f"{instance_to_slot[instance_id]}@{layer_suffix}"
            if isinstance(pose_detail, Mapping):
                stacked_pose_detail = dict(pose_detail)
                stacked_pose_detail["translation"] = translation
                stacked_pose_detail["layer_index"] = layer_index
                stacked_pose_detail["stacking_registry"] = registry.id
                stacked_pose_detail["layer_offset"] = offset
                stacked_pose_details[stacked_instance_id] = stacked_pose_detail

    stacked_events: list[ReactionEvent] = []
    for event in candidate.events:
        for layer_index, layer_suffix in enumerate(_STACKING_LAYER_SUFFIXES):
            stacked_events.append(
                ReactionEvent(
                    id=f"{event.id}{layer_suffix}",
                    template_id=event.template_id,
                    participants=tuple(
                        replace(
                            participant,
                            monomer_instance_id=f"{participant.monomer_instance_id}{layer_suffix}",
                        )
                        for participant in event.participants
                    ),
                    product_state=event.product_state,
                    metadata={
                        **dict(event.metadata),
                        "layer_index": layer_index,
                        "stacking_registry": registry.id,
                    },
                )
            )

    reaction_templates = Counter(event.template_id for event in stacked_events)
    graph_summary = {
        "n_monomer_instances": len(monomer_poses),
        "n_reaction_events": len(stacked_events),
        "reaction_templates": dict(reaction_templates),
    }

    embedding = dict(_mapping(candidate.metadata.get("embedding")))
    if stacked_pose_details:
        embedding["poses"] = stacked_pose_details
    embedding["stacking_enabled"] = True
    embedding["stacking"] = _stacking_metadata(registry)

    score_metadata = dict(_mapping(candidate.metadata.get("score_metadata")))
    bridge_event_metrics = tuple(
        item
        for item in score_metadata.get("bridge_event_metrics", ())
        if isinstance(item, Mapping)
    )
    if bridge_event_metrics:
        duplicated_metrics = []
        for layer_index, layer_suffix in enumerate(_STACKING_LAYER_SUFFIXES):
            for metric in bridge_event_metrics:
                duplicated_metric = dict(metric)
                event_id = duplicated_metric.get("event_id")
                if event_id is not None:
                    duplicated_metric["event_id"] = f"{event_id}{layer_suffix}"
                duplicated_metric["layer_index"] = layer_index
                duplicated_metric["stacking_registry"] = registry.id
                duplicated_metrics.append(duplicated_metric)
        score_metadata["bridge_event_metrics"] = tuple(duplicated_metrics)
    score_metadata["stacking_considered"] = False

    flags = tuple(
        dict.fromkeys(
            tuple(flag for flag in candidate.flags if str(flag) != "stacking_disabled")
            + ("stacked_2d", f"stacking:{registry.id}")
        )
    )

    metadata = {
        **dict(candidate.metadata),
        "graph_summary": graph_summary,
        "instance_to_monomer": stacked_instance_to_monomer,
        "instance_to_slot": stacked_instance_to_slot,
        "embedding": embedding,
        "score_metadata": score_metadata,
        "stacking_mode": "enumerated",
        "stacking": {
            **_stacking_metadata(registry),
            "comment_suffix": stacking_comment_suffix(registry),
            "source_candidate_id": candidate.id,
        },
    }
    return replace(
        candidate,
        id=f"{candidate.id}__{registry.id}",
        state=AssemblyState(
            cell=(base_cell[0], base_cell[1], new_c),
            monomer_poses=monomer_poses,
            torsions=state.torsions,
            layer_offsets=stacked_layer_offsets,
            stacking_state=registry.id,
        ),
        events=tuple(stacked_events),
        flags=flags,
        metadata=metadata,
    )


def _candidate_cell_kind(candidate: Candidate) -> str | None:
    embedding = _mapping(candidate.metadata.get("embedding"))
    cell_kind = embedding.get("cell_kind")
    if cell_kind is None:
        return None
    text = str(cell_kind).strip()
    return text or None


def _is_eligible_2d_candidate(candidate: Candidate) -> bool:
    topology_id = _candidate_topology_id(candidate)
    if topology_id is None:
        return False
    try:
        return get_topology_hint(topology_id).dimensionality == "2D"
    except Exception:
        return False


def _candidate_topology_id(candidate: Candidate) -> str | None:
    net_plan = _mapping(candidate.metadata.get("net_plan"))
    topology_id = net_plan.get("topology")
    if topology_id is None:
        return None
    text = str(topology_id).strip()
    return text or None


def _fractional_shift_in_plane(
    cell: tuple[Vec3, Vec3, Vec3],
    shift: tuple[float, float],
) -> Vec3:
    return add(scale(cell[0], float(shift[0])), scale(cell[1], float(shift[1])))


def _stacking_metadata(registry: LayerRegistry) -> dict[str, object]:
    return {
        "id": registry.id,
        "interlayer_distance": float(registry.interlayer_distance),
        "lateral_shift_fractional": (
            float(registry.lateral_shift[0]),
            float(registry.lateral_shift[1]),
        ),
        "layer_count": 2,
    }


def _mapping(value: object) -> Mapping[str, object]:
    return value if isinstance(value, Mapping) else {}


def _string_mapping(value: object) -> dict[str, str]:
    if not isinstance(value, Mapping):
        return {}
    return {str(key): str(item) for key, item in value.items()}


def _safe_normalize(vector: Vec3) -> Vec3:
    if norm(vector) < 1.0e-8:
        return (0.0, 0.0, 1.0)
    return normalize(vector)


__all__ = [
    "LayerRegistry",
    "StackingExplorer",
    "enumerate_candidate_stackings",
    "stacking_comment_suffix",
]
