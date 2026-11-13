"""Ring attachment repair requires rigid rotations as well as translations."""

from dataclasses import replace
from math import cos, pi, sin

import pytest

from cofkit.geometry import (
    Frame,
    add,
    matmul,
    matmul_vec,
    norm,
    rotation_about_axis,
    scale,
    sub,
)
from cofkit.model import (
    AssemblyState,
    MonomerSpec,
    MotifRef,
    Pose,
    ReactiveMotif,
    ReactionEvent,
)
from cofkit.ring_geometry import (
    ATTACHMENT_MISSING_DATA,
    RingGeometryOptimizer,
    ring_attachment_report,
    ring_geometry_profile,
    ring_geometry_report,
    validate_ring_geometry,
)


def ring_case(turn_degrees=30.0, radial_offset=0.0):
    template = "triazine_trimerization"
    profile = ring_geometry_profile(template)
    # Minimal nitrile/exocyclic geometry. The off-axis fourth site tests that
    # fitting a collinear reactive/anchor pair preserves unobserved body roll.
    positions = (
        (0.0, 0.0, 0.0),
        (1.4, 0.0, 0.0),
        (-1.16, 0.0, 0.0),
        (2.1, 1.2124355653, 0.0),
    )
    motif = ReactiveMotif(
        "n",
        "nitrile",
        (0, 1, 2),
        Frame((0.0, 0.0, 0.0), (-1.0, 0.0, 0.0), (0.0, 0.0, 1.0)),
        metadata={"reactive_atom_id": 0, "anchor_atom_id": 1},
    )
    spec = MonomerSpec(
        "nitrile",
        "local nitrile attachment",
        (motif,),
        atom_symbols=("C", "C", "N", "C"),
        atom_positions=positions,
        bonds=((0, 1, 1.0), (0, 2, 3.0), (1, 3, 1.0)),
    )
    cell = ((20.0, 0.0, 0.0), (0.0, 20.0, 0.0), (0.0, 0.0, 8.0))
    center = (10.0, 10.0, 4.0)
    poses = {}
    for i in range(3):
        angle = 2 * pi * i / 3
        outward = (cos(angle), sin(angle), 0.0)
        rotate = matmul(
            rotation_about_axis((0.0, 0.0, 1.0), angle + turn_degrees * pi / 180),
            rotation_about_axis((1.0, 0.0, 0.0), 0.4),
        )
        poses[f"m{i}"] = Pose(
            add(center, scale(outward, profile.participant_radius + radial_offset)),
            rotate,
        )
    event = ReactionEvent(
        "r",
        template,
        tuple(MotifRef(f"m{i}", spec.id, motif.id) for i in range(3)),
        metadata={
            "ring_center_fractional": (0.5, 0.5, 0.5),
            "ring_normal": (0.0, 0.0, 1.0),
        },
    )
    return (event,), AssemblyState(cell=cell, monomer_poses=poses), {spec.id: spec}


def test_rigid_pose_repairs_attachment_when_ring_arrangement_is_already_exact():
    events, state, specs = ring_case()
    initial = ring_geometry_report(events, state, specs)
    assert initial.total_residual < 1e-10
    assert (
        validate_ring_geometry(events, state, specs).attachment_classification
        == "warning"
    )
    translation_only = RingGeometryOptimizer(rigid_pose_refinement=False).optimize(
        events, state, specs
    )
    assert translation_only.metrics["accepted_iterations"] == 0
    refined = RingGeometryOptimizer().optimize(events, state, specs)
    assert refined.metrics["accepted_rigid_moves"] > 0
    assert refined.metrics["improved"]
    assert refined.metrics["final_attachment_max_deviation_degrees"] < 1e-8
    assert refined.metrics["physical_relaxation"] is False
    assert ring_geometry_report(events, refined.state, specs).total_residual < 1e-10
    assert (
        validate_ring_geometry(events, refined.state, specs).classification
        == "accepted"
    )


def test_collinear_attachment_fit_preserves_roll_and_internal_distances():
    events, state, specs = ring_case()
    result = RingGeometryOptimizer().optimize(events, state, specs)
    spec = specs["nitrile"]
    for i in range(3):
        pose = result.state.monomer_poses[f"m{i}"]
        expected_rotation = matmul(
            rotation_about_axis((0.0, 0.0, 1.0), 2 * pi * i / 3),
            rotation_about_axis((1.0, 0.0, 0.0), 0.4),
        )
        assert matmul_vec(
            pose.rotation_matrix, spec.atom_positions[3]
        ) == pytest.approx(
            matmul_vec(expected_rotation, spec.atom_positions[3]), abs=1e-8
        )
        before = [
            norm(sub(a, b)) for a in spec.atom_positions for b in spec.atom_positions
        ]
        world = [
            add(pose.translation, matmul_vec(pose.rotation_matrix, p))
            for p in spec.atom_positions
        ]
        after = [norm(sub(a, b)) for a in world for b in world]
        assert after == pytest.approx(before)


def test_arrangement_improvement_cannot_hide_worse_exocyclic_geometry():
    events, before, specs = ring_case(turn_degrees=30.0, radial_offset=0.1)
    _, after, _ = ring_case(turn_degrees=50.0)
    old = ring_geometry_report(events, before, specs)
    new = ring_geometry_report(events, after, specs)
    assert new.total_residual < old.total_residual
    optimizer = RingGeometryOptimizer()
    old_attachment = ring_attachment_report(events, before, specs)
    new_attachment = ring_attachment_report(events, after, specs)
    assert optimizer._max_attachment_deviation(
        new_attachment
    ) > optimizer._max_attachment_deviation(old_attachment)
    assert not optimizer._preserves_both_channels(
        old, old_attachment, new, new_attachment
    )


def test_ring_pose_refinement_rotates_covariantly():
    events, state, specs = ring_case()
    rotation = rotation_about_axis((1.0, 2.0, 3.0), 0.7)
    rotated_events = tuple(
        replace(
            e,
            metadata={
                **e.metadata,
                "ring_normal": matmul_vec(rotation, tuple(e.metadata["ring_normal"])),
            },
        )
        for e in events
    )
    rotated_state = replace(
        state,
        cell=tuple(matmul_vec(rotation, v) for v in state.cell),
        monomer_poses={
            iid: replace(
                p,
                translation=matmul_vec(rotation, p.translation),
                rotation_matrix=matmul(rotation, p.rotation_matrix),
            )
            for iid, p in state.monomer_poses.items()
        },
    )
    optimizer = RingGeometryOptimizer()
    reference = optimizer.optimize(events, state, specs).state
    transformed = optimizer.optimize(rotated_events, rotated_state, specs).state
    for iid, p in reference.monomer_poses.items():
        q = transformed.monomer_poses[iid]
        for atom in specs["nitrile"].atom_positions:
            expected = matmul_vec(
                rotation, add(p.translation, matmul_vec(p.rotation_matrix, atom))
            )
            assert add(
                q.translation, matmul_vec(q.rotation_matrix, atom)
            ) == pytest.approx(expected, abs=1e-8)


def test_missing_attachment_data_does_not_become_a_measured_pass():
    events, state, specs = ring_case()
    spec = specs["nitrile"]
    specs = {
        spec.id: replace(
            spec, motifs=tuple(replace(m, metadata={}) for m in spec.motifs)
        )
    }
    result = RingGeometryOptimizer().optimize(events, state, specs)
    report = ring_attachment_report(events, result.state, specs)
    assert report.status == ATTACHMENT_MISSING_DATA
    assert report.n_unmeasured_participants == 3
    assert result.metrics["final_attachment_max_deviation_degrees"] is None


def test_optional_rigid_fit_failure_keeps_the_existing_seed(monkeypatch):
    events, state, specs = ring_case()

    def fail(*args, **kwargs):
        raise RuntimeError("synthetic alignment failure")

    monkeypatch.setattr("cofkit.ring_geometry.fit_rigid_transform", fail)
    result = RingGeometryOptimizer().optimize(events, state, specs)
    assert result.state == state
    assert result.metrics["warnings"] == ("RuntimeError: synthetic alignment failure",)
