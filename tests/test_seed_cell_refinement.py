"""Analytically solvable rigid periodic closure, independent of benchmark CIFs."""

from dataclasses import replace

import pytest

from cofkit.geometry import (
    dot,
    Frame,
    add,
    matmul,
    matmul_vec,
    norm,
    rotation_about_axis,
    scale,
)
from cofkit.model import (
    AssemblyState,
    MonomerInstance,
    MonomerSpec,
    MotifRef,
    Pose,
    ReactiveMotif,
    ReactionEvent,
)
from cofkit.optimizer import ContinuousOptimizer, OptimizerConfig
from cofkit.planner import AssignmentPlan, NetPlan, TopologyHint
from cofkit.reactions import ReactionLibrary, bridge_target_distance
from cofkit.scoring import CandidateScorer
from cofkit.search import AssignmentOutcome


def grid_case(dimensions=2, shear=0.0):
    # Each periodic bridge joins two opposite sites of the same rigid node.
    # Exact closure requires |axis_k| = 2 * radius_k + linkage target.
    radii = (3.0, 1.0) if dimensions == 2 else (2.0, 3.0, 4.0)
    template = ReactionLibrary.builtin().get("imine_bridge")
    target = bridge_target_distance(template)
    desired = tuple(2 * r + target for r in radii)
    seed = (
        (desired[0] + 1, desired[1] - 0.4)
        if dimensions == 2
        else (desired[0] + 0.4, desired[1] - 0.4, desired[2] + 1)
    )
    cell = (
        (seed[0], 0.0, 0.0),
        (shear, seed[1], 0.0),
        (0.0, 0.0, 8.0 if dimensions == 2 else seed[2]),
    )
    motifs, events = [], []
    for axis, radius in enumerate(radii):
        vector = tuple(float(i == axis) for i in range(3))
        normal = (0.0, 1.0, 0.0) if axis == 2 else (0.0, 0.0, 1.0)
        for sign, kind in ((1, "amine"), (-1, "aldehyde")):
            motifs.append(
                ReactiveMotif(
                    f"p{axis}_{sign}",
                    kind,
                    (),
                    Frame(scale(vector, sign * radius), scale(vector, sign), normal),
                    metadata={
                        "plane_normal": (0.0, 0.0, 1.0) if dimensions == 2 else None
                    },
                )
            )
        events.append(
            ReactionEvent(
                f"e{axis}",
                template.id,
                (
                    MotifRef("m", "node", f"p{axis}_1"),
                    MotifRef(
                        "m",
                        "node",
                        f"p{axis}_-1",
                        tuple(int(i == axis) for i in range(3)),
                    ),
                ),
            )
        )
    specs = {"node": MonomerSpec("node", "rigid node", tuple(motifs))}
    templates = {template.id: template}
    plan = AssignmentPlan(
        NetPlan(
            TopologyHint("sql" if dimensions == 2 else "pcu", f"{dimensions}D"),
            ("node",),
            (template.id,),
        ),
        {"m": "node"},
    )
    outcome = AssignmentOutcome(
        plan, (MonomerInstance("m", "node"),), tuple(events), (), 2 * dimensions
    )
    center = tuple(0.3 * sum(v[i] for v in cell) for i in range(3))
    state = AssemblyState(cell=cell, monomer_poses={"m": Pose(center)})
    return specs, templates, outcome, state, desired


@pytest.mark.parametrize("dimensions", [2, 3])
@pytest.mark.parametrize("shear", [0.0, 0.8])
def test_cell_fit_closes_independent_axes_and_shear(dimensions, shear):
    specs, templates, outcome, state, desired = grid_case(dimensions, shear)
    scorer = CandidateScorer()
    fit = ContinuousOptimizer(scorer=scorer).optimize(outcome, state, specs, templates)
    uniform = ContinuousOptimizer(
        OptimizerConfig(
            cell_shape_refinement=False, coupled_pose_cell_refinement=False
        ),
        scorer,
    ).optimize(outcome, state, specs, templates)
    report = scorer.bridge_geometry_report(outcome, fit.state, specs, templates)
    assert max(m.distance_residual for m in report.event_metrics) < 0.01
    assert report.total_residual < 0.02
    assert uniform.metrics["final_residual"] > 0.05
    # A global rotation is a free physical gauge in 3D: compare lengths and
    # angles, rather than requiring one arbitrary world-coordinate matrix.
    assert norm(fit.state.cell[0]) == pytest.approx(desired[0], abs=0.01)
    assert norm(fit.state.cell[1]) == pytest.approx(desired[1], abs=0.01)
    assert (
        abs(
            dot(fit.state.cell[0], fit.state.cell[1])
            / norm(fit.state.cell[0])
            / norm(fit.state.cell[1])
        )
        < 0.002
    )
    if dimensions == 2:
        assert fit.state.cell[2] == state.cell[2]
        assert fit.metrics["active_cell_axes"] == ("a", "b")
    else:
        assert norm(fit.state.cell[2]) == pytest.approx(desired[2], abs=0.01)
        assert fit.metrics["active_cell_axes"] == ("a", "b", "c")
    assert fit.metrics["physical_relaxation"] is False


def _rotate_state(state, rotation):
    return replace(
        state,
        cell=tuple(matmul_vec(rotation, v) for v in state.cell),
        monomer_poses={
            iid: replace(
                pose,
                translation=matmul_vec(rotation, pose.translation),
                rotation_matrix=matmul(rotation, pose.rotation_matrix),
            )
            for iid, pose in state.monomer_poses.items()
        },
    )


@pytest.mark.parametrize("dimensions", [2, 3])
def test_cell_refinement_rotates_covariantly(dimensions):
    specs, templates, outcome, state, _ = grid_case(dimensions, 0.8)
    rotation = rotation_about_axis((1.0, 2.0, 3.0), 0.7)
    optimizer = ContinuousOptimizer()
    reference = optimizer.optimize(outcome, state, specs, templates).state
    rotated = optimizer.optimize(
        outcome, _rotate_state(state, rotation), specs, templates
    ).state
    for before, after in zip(reference.cell, rotated.cell):
        assert after == pytest.approx(matmul_vec(rotation, before), abs=1e-8)
    for iid in reference.monomer_poses:
        assert rotated.monomer_poses[iid].translation == pytest.approx(
            matmul_vec(rotation, reference.monomer_poses[iid].translation), abs=1e-8
        )
    if dimensions == 2:
        assert rotated.cell[2] == _rotate_state(state, rotation).cell[2]


def test_uniform_fallback_moves_centers_in_the_actual_cell_plane():
    specs, templates, outcome, state, _ = grid_case()
    rotation = rotation_about_axis((1.0, 2.0, 3.0), 0.7)
    optimizer = ContinuousOptimizer(OptimizerConfig(cell_shape_refinement=False))
    reference = optimizer._scale_lateral_cell(state, outcome, specs, templates)
    rotated = optimizer._scale_lateral_cell(
        _rotate_state(state, rotation), outcome, specs, templates
    )
    assert rotated.monomer_poses["m"].translation == pytest.approx(
        matmul_vec(rotation, reference.monomer_poses["m"].translation), abs=1e-10
    )


def test_rank_deficient_fit_does_not_invent_a_cell_axis():
    specs, templates, outcome, state, _ = grid_case(3)
    incomplete = replace(outcome, events=outcome.events[:2])
    assert (
        ContinuousOptimizer()._fit_cell_shape(state, incomplete, specs, templates)
        is state
    )


def test_affine_cell_trust_bounds_prevent_collapse_and_inversion():
    optimizer = ContinuousOptimizer()
    cell = ((4.0, 0.0, 0.0), (0.0, 5.0, 0.0), (0.0, 0.0, 6.0))
    assert optimizer._cell_step_is_bounded(cell, cell)
    assert not optimizer._cell_step_is_bounded(
        cell, (cell[0], cell[1], scale(cell[2], -1.0))
    )
    assert not optimizer._cell_step_is_bounded(
        cell, (cell[0], (3.99, 0.01, 0.0), cell[2])
    )
    assert not optimizer._cell_step_is_bounded(cell, tuple(scale(v, 2.0) for v in cell))


def test_source_rigid_geometry_is_not_deformed():
    specs, templates, outcome, state, _ = grid_case(3)
    positions = ((1.0, 2.0, 3.0), (-2.0, 1.0, 0.0), (0.0, -1.0, 1.0))
    specs = {
        "node": replace(
            specs["node"], atom_symbols=("C", "C", "C"), atom_positions=positions
        )
    }
    before = tuple(
        norm(tuple(a - b for a, b in zip(x, y))) for x in positions for y in positions
    )
    optimized = ContinuousOptimizer().optimize(outcome, state, specs, templates)
    pose = optimized.state.monomer_poses["m"]
    world = tuple(
        add(matmul_vec(pose.rotation_matrix, p), pose.translation) for p in positions
    )
    after = tuple(
        norm(tuple(a - b for a, b in zip(x, y))) for x in world for y in world
    )
    assert after == pytest.approx(before)
    assert specs["node"].atom_positions == positions


@pytest.mark.parametrize("angle", [0.3, 0.7, 1.1])
@pytest.mark.parametrize("displacement", [(0.0, 0.0, 0.0), (0.5, 0.5, 0.0)])
def test_coupled_pose_proposal_improves_shared_linker_closure(angle, displacement):
    from test_embedding_scoring import build_two_event_linker_case

    specs, templates, outcome, state = build_two_event_linker_case()
    poses = dict(state.monomer_poses)
    linker = poses["mL"]
    poses["mL"] = replace(
        linker,
        rotation_matrix=matmul(
            rotation_about_axis((0.0, 0.0, 1.0), angle), linker.rotation_matrix
        ),
        translation=add(linker.translation, displacement),
    )
    state = replace(state, monomer_poses=poses)
    greedy = ContinuousOptimizer(
        OptimizerConfig(max_iterations=1, coupled_pose_cell_refinement=False)
    ).optimize(outcome, state, specs, templates)
    coupled = ContinuousOptimizer(OptimizerConfig(max_iterations=1)).optimize(
        outcome, state, specs, templates
    )
    assert coupled.metrics["final_residual"] <= greedy.metrics["final_residual"]
    if any(displacement):
        assert coupled.metrics["final_residual"] < greedy.metrics["final_residual"]
    else:
        # Pure orientation repair must not open the initially exact bonds.
        report = CandidateScorer().bridge_geometry_report(outcome, coupled.state, specs, templates)
        assert max(m.distance_residual for m in report.event_metrics) < 1e-9
    assert coupled.metrics["physical_relaxation"] is False
