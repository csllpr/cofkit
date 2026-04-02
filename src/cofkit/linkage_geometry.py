from __future__ import annotations

from .geometry import Vec3, norm, scale, sub
from .model import MonomerSpec


IMINE_EFFECTIVE_ORIGIN_RETRACTION_FRACTION = 0.11
AZINE_EFFECTIVE_ORIGIN_RETRACTION_FRACTION = 0.08


def effective_motif_origin(
    template_id: str | None,
    monomer: MonomerSpec,
    motif,
) -> Vec3:
    origin = motif.frame.origin
    if template_id not in {"imine_bridge", "azine_bridge"}:
        return origin
    if not monomer.atom_positions:
        return origin
    if template_id == "imine_bridge" and motif.kind not in {"amine", "aldehyde"}:
        return origin
    if template_id == "azine_bridge" and motif.kind not in {"hydrazine", "aldehyde"}:
        return origin

    reactive_atom_id = motif.metadata.get("reactive_atom_id")
    anchor_atom_id = motif.metadata.get("anchor_atom_id")
    if not isinstance(reactive_atom_id, int) or not isinstance(anchor_atom_id, int):
        return origin
    if reactive_atom_id >= len(monomer.atom_positions) or anchor_atom_id >= len(monomer.atom_positions):
        return origin

    reactive_position = monomer.atom_positions[reactive_atom_id]
    anchor_position = monomer.atom_positions[anchor_atom_id]
    anchor_to_reactive = sub(reactive_position, anchor_position)
    if norm(anchor_to_reactive) < 1e-8:
        return origin
    retraction_fraction = (
        IMINE_EFFECTIVE_ORIGIN_RETRACTION_FRACTION
        if template_id == "imine_bridge"
        else AZINE_EFFECTIVE_ORIGIN_RETRACTION_FRACTION
    )
    return sub(
        reactive_position,
        scale(anchor_to_reactive, retraction_fraction),
    )
