"""Structural and charge checks shared by external calculation adapters."""

from __future__ import annotations

import math
from pathlib import Path

import gemmi


def read_ordered_structure(path: Path):
    block = gemmi.cif.read_file(str(path)).sole_block()
    for tag in (
        "length_a",
        "length_b",
        "length_c",
        "angle_alpha",
        "angle_beta",
        "angle_gamma",
    ):
        value = gemmi.cif.as_number(block.find_value("_cell_" + tag) or "?")
        if not math.isfinite(value) or value <= 0:
            raise ValueError(
                f"Calculation requires an explicit, finite positive _cell_{tag}."
            )
    small = gemmi.make_small_structure_from_block(block)
    cell = small.cell
    parameters = (cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma)
    if (
        not all(math.isfinite(v) and v > 0 for v in parameters)
        or not math.isfinite(cell.volume)
        or cell.volume <= 0
    ):
        raise ValueError("Calculation requires a finite, nondegenerate unit cell.")
    if not all(0 < angle < 180 for angle in (cell.alpha, cell.beta, cell.gamma)):
        raise ValueError("Cell angles must be strictly between 0 and 180 degrees.")
    if not small.sites:
        raise ValueError("Calculation requires explicit atom sites.")
    for site in small.sites:
        if not math.isfinite(site.occ) or site.occ != 1.0:
            raise ValueError(
                "Fractional occupancy/disorder is unsupported; supply a resolved configuration."
            )
        if not all(math.isfinite(v) for v in site.fract):
            raise ValueError(f"Nonfinite coordinates for atom {site.label!r}.")
        if site.element.atomic_number == 0:
            raise ValueError(f"Unknown element for atom {site.label!r}.")
    return block, small


def validate_charge_assignment(
    source: Path, charged: Path, *, target_charge: float, tolerance: float
) -> dict[str, object]:
    """Verify a bijection before restoring input labels in the charged CIF.

    Geometry tolerance (1e-4 Å) accommodates decimal CIF serialization, not
    structural relaxation. EQeq is required to leave the geometry unchanged.
    """
    _, original = read_ordered_structure(source)
    block, result = read_ordered_structure(charged)
    if len(original.sites) != len(result.sites):
        raise ValueError("EQeq changed the number of atom sites.")
    for a, b in zip(original.cell.parameters, result.cell.parameters):
        if not math.isclose(a, b, rel_tol=1e-8, abs_tol=1e-8):
            raise ValueError("EQeq changed the unit cell.")
    charges = block.find_loop("_atom_site_charge")
    if not charges:
        charges = block.find_loop("_atom_site_partial_charge")
    if len(charges) != len(result.sites):
        raise ValueError("EQeq must supply one charge for every atom.")
    values = [gemmi.cif.as_number(raw) for raw in charges]
    if not all(math.isfinite(q) for q in values):
        raise ValueError("EQeq charges must all be finite.")
    total = math.fsum(values)
    if abs(total - target_charge) > tolerance:
        raise ValueError(
            f"EQeq net charge {total:g} e differs from target {target_charge:g} e by more than {tolerance:g} e; increase charge precision or review the model."
        )
    original_labels = [site.label for site in original.sites]
    if len(set(original_labels)) != len(original_labels):
        raise ValueError("Input atom labels must be unique for charge assignment.")
    output_labels = block.find_loop("_atom_site_label")
    by_label = {site.label: j for j, site in enumerate(original.sites)}
    search = None
    used: set[int] = set()
    for i, atom in enumerate(result.sites):
        if atom.label in by_label:
            candidates = [by_label[atom.label]]
        else:
            if search is None:
                # 1 Å bins control search cost, independently of the 1e-4 Å
                # serialization tolerance. Avoid a quadratic atom-pair scan.
                search = gemmi.NeighborSearch(original, 1.0).populate(include_h=True)
            candidates = sorted(
                {
                    int(mark.atom_idx)
                    for mark in search.find_site_neighbors(
                        atom, min_dist=0, max_dist=1e-4
                    )
                }
            )
        candidates = [
            j
            for j in candidates
            if atom.element == original.sites[j].element
            and original.cell.find_nearest_pbc_image(
                atom.fract, original.sites[j].fract, 0
            ).dist()
            <= 1e-4
        ]
        if len(candidates) != 1 or candidates[0] in used:
            raise ValueError(
                "EQeq atom mapping is ambiguous or geometry/elements changed."
            )
        j = candidates[0]
        used.add(j)
        output_labels[i] = gemmi.cif.quote(original.sites[j].label)
    document = gemmi.cif.Document()
    document.add_copied_block(block)
    document.write_file(str(charged))
    return {
        "net_charge": total,
        "target_charge": target_charge,
        "tolerance_e": tolerance,
        "mapping": "verified_elements_and_periodic_coordinates",
    }
