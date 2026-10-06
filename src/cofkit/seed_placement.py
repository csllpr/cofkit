"""Bounded atomistic screening of rigid seed placements before stacking.

This is construction, not physical relaxation: only linker roll and lateral
placement/cell span change. Reaction realization keeps owning bridge geometry;
the validator keeps owning periodic bond exclusions and contact thresholds.
"""

from __future__ import annotations

import sys
from dataclasses import asdict, dataclass, replace
from math import isfinite, pi
from typing import Mapping

import gemmi

from .batch_models import SeedPlacementDiagnostics, SeedPlacementMove
from .cif import CIFWriter
from .geometry import (
    add,
    cross,
    dot,
    matmul,
    matmul_vec,
    norm,
    rotation_about_axis,
    scale,
    sub,
)
from .model import Candidate, MonomerSpec
from .reactions import bridge_target_distance
from .topologies import get_topology_hint
from .validation import CoarseStructureValidator, measure_inter_instance_bond_distances


@dataclass(frozen=True)
class SeedPlacementConfig:
    enabled: bool = True
    # Derived: twelve evenly spaced rolls cover 360 degrees at 30-degree
    # spacing. These are construction alternatives, not a torsion potential.
    roll_samples: int = 12
    # Heuristic — pending calibration: two coordinate-descent sweeps allow
    # interacting linkers to respond to one another without unbounded search.
    roll_sweeps: int = 2
    # Heuristic — pending calibration: cap full realization/contact scans
    # per seed, including large nets with many independent linker copies.
    max_evaluations: int = 48
    # Heuristic — pending calibration: search up to ten percent expansion,
    # with the smallest successful one-percent lateral step preferred.
    lateral_scale_step: float = 0.01
    lateral_scale_steps: int = 10
    # Heuristic — pending calibration: tolerate CIF coordinate rounding only;
    # candidates cannot trade clashes for stretched realized bridge bonds.
    bond_residual_tolerance_angstrom: float = 0.002

    def __post_init__(self):
        if (
            self.roll_samples < 1
            or self.roll_sweeps < 0
            or self.max_evaluations < 1
            or self.lateral_scale_steps < 0
        ):
            raise ValueError(
                "seed placement search budgets must be nonnegative, with at least one sample/evaluation"
            )
        if (
            not isfinite(self.lateral_scale_step)
            or not isfinite(self.bond_residual_tolerance_angstrom)
            or self.lateral_scale_step < 0
            or self.bond_residual_tolerance_angstrom < 0
        ):
            raise ValueError(
                "seed placement scale step and bond tolerance must be nonnegative"
            )


_SCREENED_TEMPLATES = frozenset(
    {"imine_bridge", "azine_bridge", "hydrazone_bridge", "vinylene_bridge"}
)


def _roll(candidate, monomer, instance_ids, angle):
    first, second = (motif.frame.origin for motif in monomer.motifs)
    rotation = rotation_about_axis(sub(second, first), angle)
    poses = dict(candidate.state.monomer_poses)
    for instance_id in instance_ids:
        pose = poses[instance_id]
        # Rotate about the reactive-origin line: both connection origins
        # stay fixed, as do all internal distances of the rigid monomer.
        local_shift = sub(first, matmul_vec(rotation, first))
        poses[instance_id] = replace(
            pose,
            rotation_matrix=matmul(pose.rotation_matrix, rotation),
            translation=add(
                pose.translation, matmul_vec(pose.rotation_matrix, local_shift)
            ),
        )
    return replace(candidate, state=replace(candidate.state, monomer_poses=poses))


def _expand(candidate, factor):
    a, b, c = candidate.state.cell
    determinant = dot(a, cross(b, c))
    poses = {}
    for instance_id, pose in candidate.state.monomer_poses.items():
        fa = dot(pose.translation, cross(b, c)) / determinant
        fb = dot(pose.translation, cross(c, a)) / determinant
        lateral = add(scale(a, fa), scale(b, fb))
        poses[instance_id] = replace(
            pose, translation=add(pose.translation, scale(lateral, factor - 1.0))
        )
    return replace(
        candidate,
        state=replace(
            candidate.state,
            cell=(scale(a, factor), scale(b, factor), c),
            monomer_poses=poses,
        ),
    )


def refine_seed_placement(
    candidate: Candidate,
    monomer_specs: Mapping[str, MonomerSpec],
    *,
    config: SeedPlacementConfig | None = None,
) -> Candidate:
    """Keep a contact-improving placement, or preserve the original on failure.

    The screen runs on the unstacked node/linker layer. Node/node assemblies,
    untested reaction families, coarse monomers, and 3D nets retain their
    existing construction path: irregular multi-port nodes need a joint
    attachment solve rather than this rigid-linker search.
    Trials use final realized coordinates and periodic graph exclusions, with
    no filesystem artifacts or external simulation executables.
    """
    config = config or SeedPlacementConfig()
    if not config.enabled:
        return candidate
    template_ids = {event.template_id for event in candidate.events}
    topology_id = candidate.metadata.get("net_plan", {}).get("topology")
    if (
        len(template_ids) != 1
        or not template_ids <= _SCREENED_TEMPLATES
        or not topology_id
    ):
        return candidate
    try:
        dimensionality = get_topology_hint(topology_id).dimensionality
    except (KeyError, ValueError):
        return candidate
    if dimensionality != "2D":
        return candidate
    if not all(spec.atom_positions for spec in monomer_specs.values()):
        return candidate
    if not any(len(spec.motifs) == 2 for spec in monomer_specs.values()):
        return candidate
    validator = CoarseStructureValidator()
    writer = CIFWriter()
    target = bridge_target_distance(next(iter(template_ids)))
    evaluations = 0

    def evaluate(trial):
        nonlocal evaluations
        export = writer.export_candidate(trial, monomer_specs)
        if export.mode != "atomistic":
            raise ValueError("seed contact screening requires atomistic realization")
        block = gemmi.cif.read_string(export.text).sole_block()
        small = gemmi.make_small_structure_from_block(block)
        scan = validator._nonbonded_contact_scan(
            small, *validator._bond_graph_exclusions(block)
        )
        bonds = measure_inter_instance_bond_distances(block)
        if not bonds:
            raise ValueError("seed contact screening found no realized bridge bonds")
        # Protect each realized edge, not just the worst residual. Otherwise
        # an already stretched edge could mask damage to a different bridge.
        residuals = {
            (bond.label_1, bond.label_2, tuple(bond.image)): abs(bond.distance - target)
            for bond in bonds
        }
        count = int(scan["n_heavy_atom_clash_pairs"]) + int(
            scan["n_excluded_pair_severe_overlaps"]
        )
        clearance = 1.0
        if scan["heavy_atom_clash_min"]:
            clearance = min(
                clearance,
                scan["heavy_atom_clash_min"]["vdw_ratio"]
                / validator.thresholds.min_nonbonded_heavy_vdw_ratio,
            )
        detail = scan["excluded_pair_severe_overlap_min"]
        if detail:
            hydrogen = any(
                label.rsplit("_", 1)[-1].startswith(("H", "D"))
                for label in detail["labels"]
            )
            floor = validator.thresholds.min_nonbonded_heavy_distance
            if detail["exclusion"] == "bonded_1-2":
                floor = validator.thresholds.severe_overlap_bonded_distance
            elif detail["exclusion"] == "torsion_1-4" and not hydrogen:
                floor = validator.thresholds.hard_min_nonbonded_heavy_plain_distance
            clearance = min(clearance, detail["distance"] / floor)
        evaluations += 1
        return (count, -clearance), residuals, scan

    try:
        initial_key, initial_residuals, initial_scan = evaluate(candidate)
        if initial_key[0] == 0:
            return candidate
        best, best_key, best_scan = candidate, initial_key, initial_scan
        moves = []

        def consider(trial, move, previous_moves):
            nonlocal best, best_key, best_scan, moves
            if evaluations >= config.max_evaluations:
                return
            key, residuals, scan = evaluate(trial)
            protected_bonds = residuals.keys() == initial_residuals.keys() and all(
                residual
                <= initial_residuals[edge] + config.bond_residual_tolerance_angstrom
                for edge, residual in residuals.items()
            )
            if protected_bonds and key < best_key:
                best, best_key, best_scan = trial, key, scan
                # Alternatives in a sweep all start from the same base.
                # Record only the lineage of the selected placement, rather
                # than implying every accepted trial was applied in sequence.
                moves = [*previous_moves, move]

        def expand_span():
            base = best
            base_moves = tuple(moves)
            base_factor = norm(base.state.cell[0]) / norm(candidate.state.cell[0])
            for step in range(1, config.lateral_scale_steps + 1):
                if best_key[0] == 0 or evaluations >= config.max_evaluations:
                    break
                factor = 1.0 + config.lateral_scale_step * step
                if factor <= base_factor:
                    continue
                consider(
                    _expand(base, factor / base_factor),
                    SeedPlacementMove(kind="lateral_expansion", factor=factor),
                    base_moves,
                )

        # A small span correction is the cheapest way to remove clashes
        # caused by anticipated ideal bridge bends at crowded attachments.
        expand_span()
        instance_map = candidate.metadata["instance_to_monomer"]
        for _ in range(config.roll_sweeps):
            if best_key[0] == 0 or evaluations >= config.max_evaluations:
                break
            for monomer_id, spec in monomer_specs.items():
                if (
                    len(spec.motifs) != 2
                    or norm(
                        sub(spec.motifs[1].frame.origin, spec.motifs[0].frame.origin)
                    )
                    == 0
                ):
                    continue
                ids = tuple(i for i, m in instance_map.items() if m == monomer_id)
                groups = [ids] + [(i,) for i in ids] if len(ids) > 1 else [ids]
                for group in groups:
                    base = best
                    base_moves = tuple(moves)
                    for sample in range(1, config.roll_samples):
                        angle = 2.0 * pi * sample / config.roll_samples
                        consider(
                            _roll(base, spec, group, angle),
                            SeedPlacementMove(
                                kind="linker_roll",
                                instances=group,
                                angle_degrees=360.0 * sample / config.roll_samples,
                            ),
                            base_moves,
                        )
                        if best_key[0] == 0 or evaluations >= config.max_evaluations:
                            break
                    if best_key[0] == 0 or evaluations >= config.max_evaluations:
                        break
                if best_key[0] == 0 or evaluations >= config.max_evaluations:
                    break

        expand_span()
        diagnostics = SeedPlacementDiagnostics(
            evaluations=evaluations,
            moves=tuple(moves),
            initial_heavy_clashes=initial_scan["n_heavy_atom_clash_pairs"],
            initial_excluded_overlaps=initial_scan["n_excluded_pair_severe_overlaps"],
            final_heavy_clashes=best_scan["n_heavy_atom_clash_pairs"],
            final_excluded_overlaps=best_scan["n_excluded_pair_severe_overlaps"],
            resolved=best_key[0] == 0,
            warnings=()
            if best_key[0] == 0
            else (
                "bounded seed placement could not clear all reported contacts; final validation remains authoritative",
            ),
        )
        return replace(
            best, metadata={**best.metadata, "seed_placement": asdict(diagnostics)}
        )
    except Exception as exc:
        warning = f"{type(exc).__name__}: {exc}"
        print(
            f"warning: seed placement screening failed ({warning}); retaining the original seed",
            file=sys.stderr,
        )
        return replace(
            candidate,
            metadata={
                **candidate.metadata,
                "seed_placement": asdict(
                    SeedPlacementDiagnostics(
                        evaluations=evaluations, warnings=(warning,)
                    )
                ),
            },
        )
