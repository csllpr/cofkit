"""Construction invariants and bounded/failure-isolated seed screening."""

from dataclasses import replace
from math import pi
from unittest.mock import patch

import pytest

from cofkit import (
    AssemblyState,
    BatchGenerationConfig,
    BatchMonomerRecord,
    BatchStructureGenerator,
    Candidate,
    Frame,
    MonomerSpec,
    Pose,
    ReactiveMotif,
    SeedPlacementConfig,
)
from cofkit.cofid import cofid_to_build_request
from cofkit.geometry import add, distance, matmul_vec
from cofkit.seed_placement import _expand, _roll, refine_seed_placement


def _candidate():
    return Candidate(
        id="seed",
        score=None,
        state=AssemblyState(
            cell=((20.0, 0.0, 0.0), (5.0, 20.0, 0.0), (0.0, 0.0, 8.0)),
            monomer_poses={"m1": Pose(translation=(7.0, 6.0, 2.0))},
        ),
        events=(),
    )


def test_roll_preserves_connection_origins_and_rigid_internal_geometry():
    spec = MonomerSpec(
        id="linker",
        name="linker",
        motifs=(
            ReactiveMotif("a", "amine", (0,), Frame.xy((1.0, 2.0, 3.0))),
            ReactiveMotif("b", "amine", (1,), Frame.xy((5.0, 4.0, 3.0))),
        ),
        atom_symbols=("N", "N", "C"),
        atom_positions=((1.0, 2.0, 3.0), (5.0, 4.0, 3.0), (3.0, 1.0, 5.0)),
    )
    original = _candidate()
    rolled = _roll(original, spec, ("m1",), pi / 3)

    def world(candidate, point):
        pose = candidate.state.monomer_poses["m1"]
        return add(pose.translation, matmul_vec(pose.rotation_matrix, point))

    for motif in spec.motifs:
        assert world(rolled, motif.frame.origin) == pytest.approx(
            world(original, motif.frame.origin)
        )
    assert distance(
        world(rolled, spec.atom_positions[0]), world(rolled, spec.atom_positions[2])
    ) == pytest.approx(distance(spec.atom_positions[0], spec.atom_positions[2]))
    assert rolled.state.cell == original.state.cell
    assert world(rolled, spec.atom_positions[2]) != pytest.approx(
        world(original, spec.atom_positions[2])
    )


def test_lateral_expansion_preserves_c_and_fractional_placement_in_oblique_cell():
    original = _candidate()
    expanded = _expand(original, 1.05)
    assert expanded.state.cell[2] == original.state.cell[2]
    assert expanded.state.cell[0] == pytest.approx((21.0, 0.0, 0.0))
    assert expanded.state.cell[1] == pytest.approx((5.25, 21.0, 0.0))
    assert expanded.state.monomer_poses["m1"].translation == pytest.approx(
        (7.35, 6.3, 2.0)
    )
    assert (
        expanded.state.monomer_poses["m1"].rotation_matrix
        == original.state.monomer_poses["m1"].rotation_matrix
    )


def _azine_generator(config):
    request = cofid_to_build_request(
        "3:aldehyde:O=Cc1c(O)c(C=O)c(O)c(C=O)c1O.2:hydrazine:NN&&hcb&&azine"
    )
    records = tuple(
        BatchMonomerRecord(
            id=f"reactant_{index}",
            name=f"reactant_{index}",
            smiles=monomer.canonical_smiles,
            motif_kind=monomer.motif_kind,
            expected_connectivity=monomer.connectivity,
        )
        for index, monomer in enumerate(request.monomers)
    )
    generator = BatchStructureGenerator(
        config=BatchGenerationConfig(
            allowed_reactions=(request.template_id,),
            topology_ids=(request.topology,),
            enumerate_all_topologies=False,
            rdkit_num_conformers=2,
            seed_placement_config=config,
        )
    )
    return generator, records


def test_evaluation_budget_keeps_unresolved_seed_and_records_diagnostics():
    generator, records = _azine_generator(SeedPlacementConfig(max_evaluations=1))
    summary, candidate = generator.generate_pair_candidate(*records)
    assert summary.status == "ok"
    assert candidate is not None
    diagnostics = summary.metadata["seed_placement"]
    assert diagnostics["evaluations"] == 1
    assert not diagnostics["resolved"]
    assert (
        diagnostics["initial_excluded_overlaps"]
        == diagnostics["final_excluded_overlaps"]
        == 6
    )
    assert diagnostics["moves"] == ()
    assert diagnostics["warnings"]
    assert not diagnostics["physical_relaxation"]


def test_screen_failure_is_nonfatal_and_preserves_original_state(capsys):
    generator, records = _azine_generator(SeedPlacementConfig(enabled=False))
    _, candidate = generator.generate_pair_candidate(*records)
    specs = {record.id: generator.build_monomer(record).monomer for record in records}
    with patch(
        "cofkit.seed_placement.CIFWriter.export_candidate",
        side_effect=RuntimeError("failed trial"),
    ):
        screened = refine_seed_placement(candidate, specs)
    assert screened.state is candidate.state
    assert screened.events == candidate.events
    assert screened.metadata["seed_placement"]["warnings"] == (
        "RuntimeError: failed trial",
    )
    assert "warning: seed placement screening failed" in capsys.readouterr().err


def test_custom_motifs_keep_legacy_prior_and_can_explicitly_disable_it():
    motif = ReactiveMotif("a", "amine", (), Frame.xy())
    assert motif.plane_normal == motif.frame.normal
    uncertain = replace(motif, metadata={"plane_normal": None})
    assert uncertain.plane_normal == (0.0, 0.0, 0.0)
    assert uncertain.frame.normal == (0.0, 0.0, 1.0)
