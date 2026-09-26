from __future__ import annotations

import sys
from collections import Counter
from collections.abc import Callable
from dataclasses import dataclass, replace
from math import acos, degrees
from typing import Mapping

from .geometry import (
    LayerSpanReport,
    Vec3,
    add,
    classify_2d_cell,
    dot,
    layer_normal_axis,
    matmul_vec,
    measure_layer_z_span,
    norm,
    scale,
)
from .model import AssemblyState, Candidate, Pose, ReactionEvent
from .topologies import get_topology_hint
from .vdw import DEFAULT_NONBONDED_SEARCH_RADIUS, min_periodic_pair_contact


_STACKING_LAYER_SUFFIXES: tuple[str, str] = ("L0", "L1")

# Base-cell c tilts at or below this angle (degrees) are numerical noise: the
# stacked cell is still reported as built along the base c direction.
_C_TILT_PROVENANCE_TOLERANCE_DEG = 1e-4

_UNSET = object()


@dataclass(frozen=True)
class LayerRegistry:
    """Describes one named stacking registry for a bilayer COF export.

    Semantics
    ---------
    ``interlayer_distance`` here is a **nuclear-plane clearance**: the gap
    between the extreme nuclear planes of adjacent layers.  It is *not* the
    mean-plane repeat published in the literature.  The centre-to-centre
    distance between adjacent layers is::

        c2c = interlayer_clearance + layer_z_span

    and the full c axis of the bilayer cell is::

        c = 2 * c2c

    built along the layer normal ``n = normalize(a × b)`` of the base cell, so
    the normal repeat of the shipped cell is exactly ``c2c`` even when the
    base cell's c was tilted (the tilt is recorded as provenance in the
    stacking metadata).  The quantity to compare to experimental interlayer
    spacings is ``c / 2``.

    Default clearances (3.4 / 3.5 / 3.6 Å) are heuristic nuclear-plane
    clearances, registry-dependent; they are not fitted to specific chemistries.
    """

    id: str
    lateral_shift: tuple[float, float] = (0.0, 0.0)
    interlayer_distance: float = 3.4  # nuclear-plane clearance (see docstring)


class StackingExplorer:
    """Enumerates built-in bilayer registries for exported 2D COF candidates.

    ``lateral_shift`` is expressed as a fractional shift in the in-plane a/b
    basis and is applied to the second layer.

    For hexagonal cells the AB shift depends on the cell setting:
    * 60° setting (γ ≈ 60°): AB shift = (1/3, 1/3)  — Bernal, vertex→pore
    * 120° setting (γ ≈ 120°): AB shift = (1/3, 2/3) — same geometry, different basis

    ``slipped`` and ``AB_square`` are setting-independent.
    """

    def __init__(
        self,
        *,
        cell_kind: str | None = None,
        cell_setting: str | None = None,
    ) -> None:
        self.cell_kind = str(cell_kind or "").strip().lower() or None
        self.cell_setting = str(cell_setting or "").strip().lower() or None

    def enumerate_registries(self) -> tuple[LayerRegistry, ...]:
        if self.cell_kind == "hexagonal":
            # Select AB shift by cell setting (W2.2)
            if self.cell_setting == "120deg":
                ab_shift = (1.0 / 3.0, 2.0 / 3.0)
            else:
                # 60deg setting (default) or unknown
                ab_shift = (1.0 / 3.0, 1.0 / 3.0)
            return (
                LayerRegistry(id="AA", lateral_shift=(0.0, 0.0), interlayer_distance=3.4),
                LayerRegistry(id="AB", lateral_shift=ab_shift, interlayer_distance=3.5),
                LayerRegistry(id="slipped", lateral_shift=(0.5, 0.0), interlayer_distance=3.6),
            )
        return (
            LayerRegistry(id="AA", lateral_shift=(0.0, 0.0), interlayer_distance=3.4),
            LayerRegistry(id="AB", lateral_shift=(0.5, 0.5), interlayer_distance=3.5),  # hollow-site registry
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
    monomer_specs: Mapping[str, object] | None = None,
) -> tuple[Candidate, ...]:
    """Expand *candidate* into one variant per requested stacking registry.

    Parameters
    ----------
    registry_ids:
        Registry names to apply.  Empty → return the candidate unchanged.
    monomer_specs:
        Monomer id → MonomerSpec mapping, used by the span measurer.  When
        ``None`` and registries are requested the call proceeds with a
        ``layer_z_span`` of 0.0 and records ``stacking_skipped:span-unavailable``
        on each variant (W4.3).
    """
    if not registry_ids:
        return (candidate,)
    eligible, skip_reason = _eligibility_check(candidate)
    if not eligible:
        # Attach skip flag and warn (W4.3)
        reason_tag = f"stacking_skipped:{skip_reason}"
        flags = tuple(dict.fromkeys(candidate.flags + (reason_tag,)))
        print(
            f"warning: stacking requested but skipped for candidate {candidate.id!r}: {skip_reason}",
            file=sys.stderr,
        )
        return (replace(candidate, flags=flags),)

    cell_kind, cell_setting = _candidate_cell_classification(candidate)
    explorer = StackingExplorer(cell_kind=cell_kind, cell_setting=cell_setting)
    registries = explorer.resolve_registries(registry_ids)
    return tuple(
        _apply_layer_registry(candidate, registry, monomer_specs=monomer_specs)
        for registry in registries
    )


def stacking_comment_suffix(registry: LayerRegistry) -> str:
    return f"stacking={registry.id}"


def _apply_layer_registry(
    candidate: Candidate,
    registry: LayerRegistry,
    *,
    monomer_specs: Mapping[str, object] | None = None,
) -> Candidate:
    state = candidate.state
    if len(_STACKING_LAYER_SUFFIXES) != 2:
        raise ValueError("the current stacking enumerator supports exactly two layers")
    if registry.interlayer_distance <= 0.0:
        raise ValueError("interlayer_distance must be positive")

    base_cell = state.cell
    # Measure the span along the layer normal n = normalize(a × b): a c_hat
    # span is not invariant under in-plane periodic images when c is tilted
    # (stacking review point 2).  The same axis also builds new_c below, so
    # the measured span and the shipped cell share one reference frame; for
    # degenerate in-plane vectors layer_normal_axis warns and falls back to
    # the c direction, preserving the prior cell-construction behaviour.
    span_axis, span_axis_label = layer_normal_axis(base_cell[0], base_cell[1], fallback_axis=base_cell[2])

    # W1.2/W1.3: measure span fresh; use metadata only as cross-check
    cell_kind, cell_setting = _candidate_cell_classification(candidate)
    stacking_warnings: list[str] = []
    # Realize the product atoms once; the same realization feeds both the span
    # measurement and the W4.2 atomistic interlayer-contact self-check.
    base_instance_to_monomer = _string_mapping(candidate.metadata.get("instance_to_monomer"))
    realization = _realize_candidate(candidate, monomer_specs, base_instance_to_monomer, warnings=stacking_warnings)
    span_report = _measure_span(
        candidate, monomer_specs, span_axis, axis_label=span_axis_label,
        warnings=stacking_warnings, realization=realization,
    )
    layer_z_span = span_report.span
    metadata_span = _metadata_layer_z_span(candidate)

    if span_report.mode == "unavailable":
        stacking_warnings.append("layer_z_span unavailable; using 0.0")
        print(
            f"warning: stacking for {candidate.id!r}: layer_z_span could not be measured; using 0.0",
            file=sys.stderr,
        )

    if metadata_span is not None and layer_z_span > 0.0:
        if abs(metadata_span - layer_z_span) > 0.5:
            stacking_warnings.append(
                f"metadata layer_z_span {metadata_span:.3f} Å differs from measured {layer_z_span:.3f} Å"
            )
            print(
                f"warning: stacking for {candidate.id!r}: metadata layer_z_span "
                f"{metadata_span:.3f} Å differs from freshly measured {layer_z_span:.3f} Å",
                file=sys.stderr,
            )

    center_to_center_distance = registry.interlayer_distance + layer_z_span

    # W4.2 self-check: c2c must exceed layer_z_span
    if center_to_center_distance <= layer_z_span:
        stacking_warnings.append(
            f"c2c ({center_to_center_distance:.3f}) ≤ layer_z_span ({layer_z_span:.3f}); layers may interpenetrate"
        )
        print(
            f"warning: stacking for {candidate.id!r}: c2c distance {center_to_center_distance:.3f} Å "
            f"≤ layer_z_span {layer_z_span:.3f} Å; layers may interpenetrate",
            file=sys.stderr,
        )

    # W4.4/W1.2 (route A): build the stacked c along the layer normal, so the
    # advertised derivation (c2c = interlayer_clearance + layer_z_span) is the
    # true normal repeat of the shipped cell even when the base cell's c was
    # tilted.  A tilted base c is not a defect in the output — the geometry is
    # correct by construction — but the basis change is recorded as provenance
    # (stacking metadata: c_axis_basis / c_axis_orthogonalized /
    # base_c_tilt_degrees) and noted on stderr.
    base_c_tilt_degrees: float | None = None
    if span_axis_label == "layer_normal":
        base_c_length = norm(base_cell[2])
        if base_c_length >= 1e-8:
            cosine = min(1.0, abs(dot(base_cell[2], span_axis)) / base_c_length)
            base_c_tilt_degrees = degrees(acos(cosine))
            if base_c_tilt_degrees > _C_TILT_PROVENANCE_TOLERANCE_DEG:
                note = (
                    f"base cell c axis tilted {base_c_tilt_degrees:.3f}° from the layer normal; "
                    "stacked c axis rebuilt along the layer normal"
                )
                stacking_warnings.append(note)
                print(
                    f"warning: stacking for {candidate.id!r}: {note}",
                    file=sys.stderr,
                )

    new_c = scale(span_axis, 2.0 * center_to_center_distance)
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
            event_metadata = {
                **dict(event.metadata),
                "layer_index": layer_index,
                "stacking_registry": registry.id,
            }
            if "ring_center_fractional" in event.metadata:
                # T1: assert the w component is zero so fractional reuse is exact
                rcf = event.metadata["ring_center_fractional"]
                if isinstance(rcf, (list, tuple)) and len(rcf) >= 3:
                    w = float(rcf[2])
                    if abs(w) > 1e-9:
                        stacking_warnings.append(
                            f"ring_center_fractional w={w:.3e} for event {event.id!r}; "
                            "fractional reuse may be inexact"
                        )
                        print(
                            f"warning: stacking for {candidate.id!r}: "
                            f"ring_center_fractional w={w:.3e} for event {event.id!r}; "
                            "fractional reuse into new cell may be inexact",
                            file=sys.stderr,
                        )
                event_metadata["ring_center_offset_cartesian"] = layer_offsets[layer_index]
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
                    metadata=event_metadata,
                )
            )

    # T2: torsion keys are passthrough-only (no per-layer suffixing yet; document as such)
    # When a consumer needs per-layer torsions, suffix keys here.

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
    # The stacked c axis is a real periodic repeat (2·c2c along the layer
    # normal), not the base candidate's vacuum-slab padding (W1.4).
    embedding["c_axis_semantics"] = "periodic_bilayer"
    embedding["stacking"] = _stacking_metadata(
        registry, layer_z_span, span_report, center_to_center_distance, cell_kind, cell_setting,
        c_axis_basis=span_axis_label, base_c_tilt_degrees=base_c_tilt_degrees,
    )

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

    # Residual contract (scoring.BridgeGeometryReport): aggregates are sums
    # over per-event metrics.  The events were just duplicated, so recompute
    # the sum from the duplicated per-event metrics in this one place; the
    # per-event mean derived downstream (model.residual_ranking_key) is then
    # invariant under stacking expansion.
    if bridge_event_metrics:
        recomputed = _recompute_total_residual(
            duplicated_metrics,
            lambda metric: float(metric.get("total_residual")),
            context=f"bridge_geometry_residual for candidate {candidate.id!r}",
        )
        if recomputed is not None:
            score_metadata["bridge_geometry_residual"] = recomputed

    ring_geometry = _mapping(score_metadata.get("ring_geometry"))
    if ring_geometry:
        score_metadata["ring_geometry"] = _duplicate_ring_geometry_metrics(ring_geometry, registry.id)

    ring_validation = dict(_mapping(candidate.metadata.get("ring_validation")))
    ring_validation_metrics = _mapping(ring_validation.get("metrics"))
    if ring_validation_metrics:
        ring_validation["metrics"] = _duplicate_ring_geometry_metrics(ring_validation_metrics, registry.id)
    if ring_validation:
        base_reasons = tuple(str(reason) for reason in ring_validation.get("reasons", ()))
        if base_reasons:
            ring_validation["reasons"] = tuple(
                f"{layer_suffix}: {reason}"
                for layer_suffix in _STACKING_LAYER_SUFFIXES
                for reason in base_reasons
            )
        ring_validation["stacking_registry"] = registry.id
        ring_validation["layer_count"] = 2

    # W4.2 self-check: atomistic interlayer contact between the two stacked
    # layers, measured on the realized product atoms (the same realization as
    # the span measurement / CIF export) across BOTH periodic c galleries
    # (images (0,0,±1)); the displayed layer pair alone can miss the minimum.
    # Per-candidate isolation: a measurement failure warns and degrades to
    # "no contact metadata" rather than aborting the stacking expansion.
    new_cell = (base_cell[0], base_cell[1], new_c)
    contact_info: dict[str, object] | None = None
    try:
        contact_info = _measure_interlayer_contact(
            monomer_poses, stacked_instance_to_monomer, monomer_specs, realization, new_cell
        )
    except Exception as exc:  # noqa: BLE001 - per-candidate failure isolation
        reason = f"{type(exc).__name__}: {exc}"
        stacking_warnings.append(f"atomistic interlayer contact measurement failed ({reason})")
        print(
            f"warning: stacking for {candidate.id!r} registry {registry.id!r}: "
            f"atomistic interlayer contact measurement failed ({reason}); "
            "no contact self-check available for this variant",
            file=sys.stderr,
        )

    stacking_clash = False
    if contact_info is not None and contact_info.get("clash"):
        stacking_clash = True
        contact_atoms = contact_info["atoms"]
        clash_detail = (
            f"stacking clash: min_interlayer_contact={contact_info['distance']:.3f} Å "
            f"(vdW ratio {contact_info['vdw_ratio']:.2f}, atoms {contact_atoms[0]}..{contact_atoms[1]}, "
            f"image {contact_info['image']})"
        )
        stacking_warnings.append(clash_detail)
        print(
            f"warning: stacking for {candidate.id!r} registry {registry.id!r}: "
            f"{clash_detail} (interpenetrating layers)",
            file=sys.stderr,
        )

    flags = tuple(
        dict.fromkeys(
            tuple(flag for flag in candidate.flags if str(flag) != "stacking_disabled")
            + ("stacked_2d", f"stacking:{registry.id}")
            + (("stacking_clash",) if stacking_clash else ())
        )
    )

    stacking_meta = {
        **_stacking_metadata(
            registry, layer_z_span, span_report, center_to_center_distance, cell_kind, cell_setting,
            c_axis_basis=span_axis_label, base_c_tilt_degrees=base_c_tilt_degrees,
        ),
        "layer_z_span": layer_z_span,
        "center_to_center_distance": center_to_center_distance,
        "comment_suffix": stacking_comment_suffix(registry),
        "source_candidate_id": candidate.id,
        "warnings": list(stacking_warnings),
    }
    if contact_info is not None:
        # Atomistic measurement; None distance means "no contact below the
        # cutoff", i.e. a lower bound, not a missing value.
        stacking_meta["min_interlayer_contact"] = contact_info["distance"]
        stacking_meta["min_interlayer_contact_mode"] = contact_info["mode"]
        stacking_meta["min_interlayer_contact_cutoff"] = contact_info["cutoff"]
        if contact_info["distance"] is not None:
            stacking_meta["min_interlayer_contact_atoms"] = contact_info["atoms"]
            stacking_meta["min_interlayer_contact_image"] = contact_info["image"]
            stacking_meta["min_interlayer_contact_involves_hydrogen"] = contact_info["involves_hydrogen"]
            stacking_meta["min_interlayer_contact_vdw_ratio"] = contact_info["vdw_ratio"]

    metadata = {
        **dict(candidate.metadata),
        "graph_summary": graph_summary,
        "instance_to_monomer": stacked_instance_to_monomer,
        "instance_to_slot": stacked_instance_to_slot,
        "embedding": embedding,
        "score_metadata": score_metadata,
        **({"ring_validation": ring_validation} if ring_validation else {}),
        "stacking_mode": "enumerated",
        "stacking": stacking_meta,
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


def _candidate_cell_classification(candidate: Candidate) -> tuple[str | None, str | None]:
    """Return (cell_kind, cell_setting) using the shared classify_2d_cell logic.

    Falls back to the embedding metadata ``cell_kind`` when the cell vectors
    are not available.
    """
    cell = getattr(getattr(candidate, "state", None), "cell", None)
    if cell is not None:
        try:
            kind, setting = classify_2d_cell(cell)
            if kind != "oblique":
                return kind, setting or None
        except (TypeError, ValueError) as exc:
            print(
                f"warning: stacking: cell classification failed for candidate "
                f"{getattr(candidate, 'id', '?')!r} ({type(exc).__name__}: {exc}); "
                "falling back to embedding metadata",
                file=sys.stderr,
            )

    embedding = _mapping(candidate.metadata.get("embedding"))
    raw_kind = embedding.get("cell_kind")
    if raw_kind is None:
        return None, None
    text = str(raw_kind).strip()
    return (text or None), None


def _is_eligible_2d_candidate(candidate: Candidate) -> bool:
    topology_id = _candidate_topology_id(candidate)
    if topology_id is None:
        return False
    try:
        return get_topology_hint(topology_id).dimensionality == "2D"
    except Exception:
        return False


def _eligibility_check(candidate: Candidate) -> tuple[bool, str]:
    """Return (eligible, reason_tag) for stacking eligibility."""
    topology_id = _candidate_topology_id(candidate)
    if topology_id is None:
        return False, "no-net-plan"
    try:
        hint = get_topology_hint(topology_id)
    except Exception:
        return False, "unknown-topology"
    if hint.dimensionality != "2D":
        return False, "not-2d"
    return True, ""


def _candidate_topology_id(candidate: Candidate) -> str | None:
    net_plan = _mapping(candidate.metadata.get("net_plan"))
    topology_id = net_plan.get("topology")
    if topology_id is None:
        return None
    text = str(topology_id).strip()
    return text or None


def _metadata_layer_z_span(candidate: Candidate) -> float | None:
    """Read layer_z_span from embedding metadata; return None if absent."""
    embedding = _mapping(candidate.metadata.get("embedding"))
    value = embedding.get("layer_z_span")
    if value is None:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _realize_candidate(
    candidate: Candidate,
    monomer_specs: Mapping[str, object] | None,
    instance_to_monomer: Mapping[str, str],
    *,
    warnings: list[str] | None = None,
) -> object | None:
    """Realize the product atoms once; warn + None on failure.

    The result feeds both the span measurement and the W4.2 atomistic
    interlayer-contact self-check, so both measure the same atoms that CIF
    export realizes.  ``None`` (no events applied or realization failure)
    downgrades callers to precursor coordinates; failures warn on stderr.
    """
    if monomer_specs is None:
        return None
    try:
        from .reaction_realization import ReactionRealizer

        return ReactionRealizer().realize(candidate, monomer_specs, instance_to_monomer)
    except (AttributeError, IndexError, KeyError, TypeError, ValueError) as exc:
        reason = f"{type(exc).__name__}: {exc}"
        if warnings is not None:
            warnings.append(
                f"atomistic product realization failed ({reason}); using precursor coordinates"
            )
        print(
            f"warning: stacking for {getattr(candidate, 'id', '?')!r}: atomistic product realization failed "
            f"({reason}); falling back to precursor coordinates",
            file=sys.stderr,
        )
        return None


def _measure_span(
    candidate: Candidate,
    monomer_specs: Mapping[str, object] | None,
    axis: Vec3,
    *,
    axis_label: str,
    warnings: list[str] | None = None,
    realization: object = _UNSET,
) -> object:
    """Measure layer z-span using the shared geometry.measure_layer_z_span.

    ``realization`` may carry a pre-computed ``ReactionRealizationResult``
    (or None) so the span and the interlayer-contact self-check share one
    realization; the default ``_UNSET`` realizes here.
    """
    if monomer_specs is None:
        # Fall back to the value recorded in embedding metadata; provenance
        # (mode/axis) is propagated from that metadata when present, never
        # stamped as a fresh measurement.
        meta_span = _metadata_layer_z_span(candidate)
        if meta_span is not None and meta_span > 0.0:
            embedding = _mapping(candidate.metadata.get("embedding"))
            meta_mode = str(embedding.get("layer_z_span_mode") or "embedding_metadata")
            meta_axis = str(embedding.get("layer_z_span_axis") or "unknown")
            return LayerSpanReport(
                span=meta_span,
                mode=meta_mode,
                axis=meta_axis,
                n_atoms=0,
                z_translation_included=True,
            )
        return LayerSpanReport(
            span=0.0,
            mode="unavailable",
            axis=axis_label,
            n_atoms=0,
            z_translation_included=False,
        )

    instance_to_monomer = _string_mapping(candidate.metadata.get("instance_to_monomer"))
    poses = candidate.state.monomer_poses

    if realization is _UNSET:
        realization = _realize_candidate(candidate, monomer_specs, instance_to_monomer, warnings=warnings)

    return measure_layer_z_span(
        poses=poses,
        monomer_specs=monomer_specs,
        instance_to_monomer=instance_to_monomer,
        axis=axis,
        axis_label=axis_label,
        realization=realization,
    )


def _measure_interlayer_contact(
    monomer_poses: dict[str, Pose],
    stacked_instance_to_monomer: dict[str, str],
    monomer_specs: Mapping[str, object] | None,
    realization: object | None,
    new_cell: tuple[Vec3, Vec3, Vec3],
) -> dict[str, object] | None:
    """Atomistic minimum interlayer contact between the two stacked layers.

    Measures on the realized product atoms (the same realization the span
    measurement and CIF export use), falling back per instance to precursor
    coordinates; the reported ``mode`` honestly distinguishes
    ``atomistic_product`` / ``precursor_coordinates`` / ``mixed``.  Every
    lattice translation within the cutoff competes — including both periodic
    galleries across the c boundary (images ``(0, 0, ±1)``), whose minimum the
    displayed layer pair alone can miss.

    Cross-layer pairs are never bonded in stacked 2D COFs, so no bond-graph
    exclusions are applied; the clash decision uses the shared vdW-ratio
    criterion plus the heavy-pair plain-distance backstop (``cofkit.vdw``),
    the same criterion as validation.

    Returns None when there are no atoms to measure (e.g. no monomer_specs).
    A ``distance`` of None in the result means "no contact below the cutoff"
    — a lower-bound statement, not a missing measurement.
    """
    if monomer_specs is None:
        return None

    layer_atoms: dict[str, list[tuple[str, str, Vec3]]] = {
        suffix: [] for suffix in _STACKING_LAYER_SUFFIXES
    }
    n_product = 0
    n_precursor = 0
    for instance_id, pose in monomer_poses.items():
        suffix = instance_id[-2:]
        if suffix not in layer_atoms:
            continue
        base_instance_id = instance_id[: -len(suffix)]
        atoms: list[tuple[str, str, Vec3]] | None = None
        if realization is not None:
            realized_atoms = getattr(realization, "atoms_by_instance", {}).get(base_instance_id)
            if realized_atoms:
                prefix = f"{base_instance_id}_"
                atoms = [
                    (
                        f"{instance_id}_{atom.label[len(prefix):]}"
                        if atom.label.startswith(prefix)
                        else f"{instance_id}_{atom.label}",
                        atom.symbol,
                        atom.local_position,
                    )
                    for atom in realized_atoms
                ]
        if atoms is None:
            monomer_id = stacked_instance_to_monomer.get(instance_id)
            monomer = monomer_specs.get(monomer_id) if monomer_id else None
            if monomer is None or not getattr(monomer, "atom_symbols", ()):
                continue
            atoms = [
                (f"{instance_id}_{symbol}{atom_index + 1}", symbol, position)
                for atom_index, (symbol, position) in enumerate(
                    zip(monomer.atom_symbols, monomer.atom_positions)
                )
            ]
            n_precursor += len(atoms)
        else:
            n_product += len(atoms)
        for label, symbol, local_position in atoms:
            world = add(matmul_vec(pose.rotation_matrix, local_position), pose.translation)
            layer_atoms[suffix].append((label, symbol, world))

    atoms_l0 = layer_atoms[_STACKING_LAYER_SUFFIXES[0]]
    atoms_l1 = layer_atoms[_STACKING_LAYER_SUFFIXES[1]]
    if not atoms_l0 or not atoms_l1:
        return None
    if n_product and n_precursor:
        mode = "mixed"
    elif n_product:
        mode = "atomistic_product"
    else:
        mode = "precursor_coordinates"

    cutoff = DEFAULT_NONBONDED_SEARCH_RADIUS
    record = min_periodic_pair_contact(new_cell, atoms_l0, atoms_l1, cutoff=cutoff)
    info: dict[str, object] = {"mode": mode, "cutoff": cutoff}
    if record is None:
        info["distance"] = None
        info["clash"] = False
        return info
    info.update(
        {
            "distance": record.distance,
            "vdw_ratio": record.vdw_ratio,
            "atoms": (record.label_i, record.label_j),
            "image": record.image,
            "involves_hydrogen": record.involves_hydrogen,
            "clash": record.clash,
        }
    )
    return info


def _fractional_shift_in_plane(
    cell: tuple[Vec3, Vec3, Vec3],
    shift: tuple[float, float],
) -> Vec3:
    return add(scale(cell[0], float(shift[0])), scale(cell[1], float(shift[1])))


def _stacking_metadata(
    registry: LayerRegistry,
    layer_z_span: float,
    span_report: object,
    center_to_center_distance: float,
    cell_kind: str | None,
    cell_setting: str | None,
    *,
    c_axis_basis: str,
    base_c_tilt_degrees: float | None,
) -> dict[str, object]:
    """Build the canonical stacking metadata dict (W3.1)."""
    span_mode = getattr(span_report, "mode", "")
    span_axis = getattr(span_report, "axis", "unknown") or "unknown"
    cell_classification: dict[str, object] = {}
    if cell_kind:
        cell_classification["kind"] = cell_kind
    if cell_setting:
        cell_classification["setting"] = cell_setting

    c_axis_orthogonalized = (
        base_c_tilt_degrees is not None
        and base_c_tilt_degrees > _C_TILT_PROVENANCE_TOLERANCE_DEG
    )
    axis_phrase = (
        "layer normal"
        if c_axis_basis == "layer_normal"
        else "base-cell c direction (degenerate in-plane fallback)"
    )

    metadata: dict[str, object] = {
        "id": registry.id,
        # Canonical keys (W3.1)
        "interlayer_clearance": float(registry.interlayer_distance),
        "registry_shift_fractional": (
            float(registry.lateral_shift[0]),
            float(registry.lateral_shift[1]),
        ),
        "cell_classification": cell_classification,
        "layer_z_span": layer_z_span,
        "layer_z_span_mode": span_mode,
        "layer_z_span_axis": span_axis,
        "center_to_center_distance": center_to_center_distance,
        "derivation": (
            "c2c = interlayer_clearance + layer_z_span; "
            f"c = layer_count * c2c along the {axis_phrase}"
        ),
        "layer_count": 2,
        # c-axis provenance: which axis the stacked c was built along, and
        # whether that rebuilt (orthogonalized) a tilted base-cell c.
        "c_axis_basis": c_axis_basis,
        "c_axis_orthogonalized": c_axis_orthogonalized,
    }
    if base_c_tilt_degrees is not None:
        metadata["base_c_tilt_degrees"] = base_c_tilt_degrees
    return metadata


def _duplicate_ring_geometry_metrics(metrics: Mapping[str, object], registry_id: str) -> dict[str, object]:
    event_metrics = tuple(item for item in metrics.get("events", ()) if isinstance(item, Mapping))
    duplicated_events: list[dict[str, object]] = []
    for layer_index, layer_suffix in enumerate(_STACKING_LAYER_SUFFIXES):
        for metric in event_metrics:
            duplicated = dict(metric)
            event_id = duplicated.get("event_id")
            if event_id is not None:
                duplicated["event_id"] = f"{event_id}{layer_suffix}"
            duplicated["layer_index"] = layer_index
            duplicated["stacking_registry"] = registry_id
            duplicated_events.append(duplicated)
    result = dict(metrics)
    result["events"] = tuple(duplicated_events)
    if duplicated_events:
        # Same sum-over-events contract as the bridge residuals; the ring
        # per-event residual mirrors ring_geometry.RingEventGeometry.residual.
        recomputed = _recompute_total_residual(
            tuple(duplicated_events),
            _ring_event_residual,
            context="ring_geometry total_residual",
        )
        if recomputed is not None:
            result["total_residual"] = recomputed
    return result


def _ring_event_residual(metric: Mapping[str, object]) -> float:
    return (
        float(metric.get("radial_rms"))
        + float(metric.get("planarity_rms"))
        + float(metric.get("angular_rms_degrees")) / 30.0
    )


def _recompute_total_residual(
    event_metrics: tuple[Mapping[str, object], ...],
    extract: Callable[[Mapping[str, object]], float],
    *,
    context: str,
) -> float | None:
    """Sum per-event residuals; warn and return None on non-numeric data."""
    total = 0.0
    for metric in event_metrics:
        try:
            total += extract(metric)
        except (TypeError, ValueError) as exc:
            print(
                f"warning: stacking: cannot recompute {context}: non-numeric "
                f"per-event residual data ({exc}); keeping the pre-stacking "
                "aggregate unchanged",
                file=sys.stderr,
            )
            return None
    return total


def _mapping(value: object) -> Mapping[str, object]:
    return value if isinstance(value, Mapping) else {}


def _string_mapping(value: object) -> dict[str, str]:
    if not isinstance(value, Mapping):
        return {}
    return {str(key): str(item) for key, item in value.items()}


__all__ = [
    "LayerRegistry",
    "StackingExplorer",
    "enumerate_candidate_stackings",
    "stacking_comment_suffix",
]
