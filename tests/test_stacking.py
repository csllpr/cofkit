"""Live-machinery stacking tests: residual aggregation contract (A9/W5.1).

The residual contract (scoring.BridgeGeometryReport) is that
``bridge_geometry_residual`` / ``ring_geometry.total_residual`` are *sums*
over per-event metrics; ranking derives a per-event mean by dividing by the
event count (model.residual_ranking_key).  Stacking duplicates a layer's
events, so these tests pin the ranking parity the doubling hack used to
approximate: the normalized residual and the unreacted-motif tie-breaker
must be invariant under stacking expansion (the candidate-id tie-breaker
intentionally changes).
"""

import contextlib
import io
from math import cos, radians, sin

import pytest
from unittest import mock

from cofkit.geometry import cross, dot, norm
from cofkit.model import (
    AssemblyState,
    Candidate,
    MonomerSpec,
    MotifRef,
    Pose,
    ReactionEvent,
    candidate_ranking_key,
    order_candidates,
)
from cofkit.stacking import LayerRegistry, _apply_layer_registry, enumerate_candidate_stackings
from cofkit.vdw import assess_pair, min_periodic_pair_contact, vdw_radius

_CELL = ((15.0, 0.0, 0.0), (0.0, 15.0, 0.0), (0.0, 0.0, 8.0))


def _bridge_candidate(candidate_id, event_residuals, n_unreacted=1):
    events = tuple(
        ReactionEvent(
            id=f"e{k}",
            template_id="imine_bridge",
            participants=(
                MotifRef(f"i{2 * k}", "m1", "mo1"),
                MotifRef(f"i{2 * k + 1}", "m2", "mo1"),
            ),
        )
        for k in range(len(event_residuals))
    )
    score_metadata = {
        "n_unreacted_motifs": n_unreacted,
        "topology": "hcb",
        # Sum over events, per the BridgeGeometryReport contract.
        "bridge_geometry_residual": sum(event_residuals),
        "bridge_event_metrics": tuple(
            {
                "event_id": f"e{k}",
                "template_id": "imine_bridge",
                "total_residual": residual,
                "distance_residual": residual,
            }
            for k, residual in enumerate(event_residuals)
        ),
    }
    return Candidate(
        id=candidate_id,
        score=None,
        state=AssemblyState(cell=_CELL),
        events=events,
        metadata={
            "net_plan": {"topology": "hcb"},
            "embedding": {"layer_z_span": 3.5, "layer_z_span_mode": "embedding_metadata"},
            "score_metadata": score_metadata,
        },
    )


def _expand(candidate, registry="AA"):
    expanded = enumerate_candidate_stackings(candidate, registry_ids=(registry,))
    assert len(expanded) == 1
    return expanded[0]


def test_stacking_expansion_keeps_ranking_residual_and_tiebreaker_invariant():
    base = _bridge_candidate("cand", [0.1, 0.2])
    stacked = _expand(base)

    base_key = candidate_ranking_key(base)
    stacked_key = candidate_ranking_key(stacked)

    # Normalized (mean per-event) residual and unreacted-motif tie-breaker
    # are invariant; only the candidate-id tie-breaker changes.
    assert stacked_key[0] == pytest.approx(base_key[0])
    assert stacked_key[1] == base_key[1]
    assert stacked_key[2] != base_key[2]

    # The exported aggregate is the honest sum over the duplicated events.
    stacked_metrics = stacked.metadata["score_metadata"]["bridge_event_metrics"]
    assert len(stacked_metrics) == 2 * len(base.events)
    assert stacked.metadata["score_metadata"]["bridge_geometry_residual"] == pytest.approx(
        sum(metric["total_residual"] for metric in stacked_metrics)
    )


def test_stacking_expansion_ring_geometry_residual_parity():
    events = (
        ReactionEvent(
            id="r1",
            template_id="boroxine_ring",
            participants=(MotifRef("i0", "m1", "mo1"),),
        ),
    )
    ring_geometry = {
        "total_residual": 0.1 + 0.1 + 9.0 / 30.0,
        "events": [
            {
                "event_id": "r1",
                "radial_rms": 0.1,
                "planarity_rms": 0.1,
                "angular_rms_degrees": 9.0,
            }
        ],
    }
    base = Candidate(
        id="ring",
        score=None,
        state=AssemblyState(cell=_CELL),
        events=events,
        metadata={
            "net_plan": {"topology": "hcb"},
            "embedding": {"layer_z_span": 3.5},
            "score_metadata": {"n_unreacted_motifs": 0, "ring_geometry": ring_geometry},
        },
    )
    stacked = _expand(base)

    base_key = candidate_ranking_key(base)
    stacked_key = candidate_ranking_key(stacked)
    assert stacked_key[0] == pytest.approx(base_key[0])
    assert stacked_key[1] == base_key[1]

    stacked_ring = stacked.metadata["score_metadata"]["ring_geometry"]
    assert len(stacked_ring["events"]) == 2
    assert stacked_ring["total_residual"] == pytest.approx(2.0 * ring_geometry["total_residual"])


def test_expanded_candidate_orders_against_unexpanded_by_id_tiebreaker():
    base = _bridge_candidate("cand", [0.1, 0.2])
    stacked = _expand(base)

    # Equal normalized residuals: ordering is decided by the id tie-breaker.
    assert [candidate.id for candidate in order_candidates([stacked, base])] == ["cand", "cand__AA"]

    # A genuinely worse unexpanded candidate still ranks after the stacked one.
    worse = _bridge_candidate("worse", [0.4, 0.5])
    ordered = [candidate.id for candidate in order_candidates([worse, stacked, base])]
    assert ordered == ["cand", "cand__AA", "worse"]


def test_stacking_residual_recompute_warns_and_keeps_aggregate_on_non_numeric():
    base = _bridge_candidate("cand", [0.1, 0.2])
    metrics = list(base.metadata["score_metadata"]["bridge_event_metrics"])
    metrics[0] = {**metrics[0], "total_residual": "n/a"}
    score_metadata = {**base.metadata["score_metadata"], "bridge_event_metrics": tuple(metrics)}
    base = Candidate(
        id=base.id,
        score=base.score,
        state=base.state,
        events=base.events,
        flags=base.flags,
        metadata={**base.metadata, "score_metadata": score_metadata},
    )

    stderr = io.StringIO()
    with contextlib.redirect_stderr(stderr):
        stacked = _expand(base)

    assert "warning: stacking: cannot recompute bridge_geometry_residual" in stderr.getvalue()
    # Aggregate left unchanged rather than silently wrong.
    assert stacked.metadata["score_metadata"]["bridge_geometry_residual"] == pytest.approx(0.3)


# ---------------------------------------------------------------------------
# W4.2 atomistic interlayer-contact self-check (audit #12)
# ---------------------------------------------------------------------------


def _stackable_candidate(monomer, translation=(0.0, 0.0, 0.0), cell=_CELL, events=()):
    return Candidate(
        id="stackme",
        score=None,
        state=AssemblyState(
            cell=cell,
            monomer_poses={"i0": Pose(translation=translation)},
        ),
        events=events,
        metadata={
            "net_plan": {"topology": "hcb"},
            "instance_to_monomer": {"i0": monomer.id},
            "embedding": {},
        },
    )


def _rod_monomer(monomer_id="m1", symbols=("C", "C"), positions=((0.0, 0.0, 0.0), (0.0, 0.0, 3.0))):
    return MonomerSpec(
        id=monomer_id,
        name=monomer_id,
        motifs=(),
        atom_symbols=symbols,
        atom_positions=positions,
    )


def test_stacking_self_check_flags_interpenetrated_bilayer():
    """Plan W7 test 7: a synthetically interpenetrated bilayer (sub-vdW
    cross-layer contact) must raise the stacking_clash flag."""
    monomer = _rod_monomer()
    candidate = _stackable_candidate(monomer)

    stderr = io.StringIO()
    with contextlib.redirect_stderr(stderr):
        stacked = _apply_layer_registry(
            candidate,
            LayerRegistry(id="AA", lateral_shift=(0.0, 0.0), interlayer_distance=0.8),
            monomer_specs={"m1": monomer},
        )

    meta = stacked.metadata["stacking"]
    assert "stacking_clash" in stacked.flags
    assert meta["min_interlayer_contact"] == pytest.approx(0.8)
    assert meta["min_interlayer_contact_mode"] == "precursor_coordinates"
    assert meta["min_interlayer_contact_involves_hydrogen"] is False
    assert meta["min_interlayer_contact_vdw_ratio"] == pytest.approx(0.8 / 3.4)
    assert tuple(meta["min_interlayer_contact_image"]) in {(0, 0, 0), (0, 0, -1)}
    assert meta["min_interlayer_contact_atoms"][0].startswith("i0L")
    assert "stacking clash" in stderr.getvalue()


def test_stacking_self_check_clean_bilayer_records_atomistic_contact():
    """A well-separated AA bilayer records the atomistic contact metadata
    without a clash flag (public enumeration path)."""
    monomer = _rod_monomer()
    candidate = _stackable_candidate(monomer)

    (stacked,) = enumerate_candidate_stackings(candidate, registry_ids=("AA",), monomer_specs={"m1": monomer})

    meta = stacked.metadata["stacking"]
    assert "stacking_clash" not in stacked.flags
    assert meta["min_interlayer_contact"] == pytest.approx(3.4)
    assert meta["min_interlayer_contact_mode"] == "precursor_coordinates"
    assert meta["min_interlayer_contact_involves_hydrogen"] is False
    assert meta["min_interlayer_contact_cutoff"] == pytest.approx(3.5)


def test_stacking_contact_records_hydrogen_involvement_without_clash():
    monomer = _rod_monomer(
        monomer_id="mh", symbols=("C", "H"), positions=((0.0, 0.0, 0.0), (0.0, 0.0, 2.7))
    )
    candidate = _stackable_candidate(monomer)

    stacked = _apply_layer_registry(
        candidate,
        LayerRegistry(id="AA", lateral_shift=(0.0, 0.0), interlayer_distance=3.4),
        monomer_specs={"mh": monomer},
    )

    meta = stacked.metadata["stacking"]
    assert meta["min_interlayer_contact"] == pytest.approx(3.4)
    assert meta["min_interlayer_contact_involves_hydrogen"] is True
    assert "stacking_clash" not in stacked.flags


def test_stacking_contact_hydrogen_involving_clash_uses_ratio_criterion():
    # H...C contact at 2.0 A: below the 2.2 A heavy floor (which must NOT
    # apply to H-involving pairs) but ratio 2.0 / 2.9 = 0.69 < 0.75 -> clash.
    monomer = _rod_monomer(
        monomer_id="mh", symbols=("C", "H"), positions=((0.0, 0.0, 0.0), (0.0, 0.0, 2.7))
    )
    candidate = _stackable_candidate(monomer)

    stacked = _apply_layer_registry(
        candidate,
        LayerRegistry(id="AA", lateral_shift=(0.0, 0.0), interlayer_distance=2.0),
        monomer_specs={"mh": monomer},
    )

    meta = stacked.metadata["stacking"]
    assert meta["min_interlayer_contact"] == pytest.approx(2.0)
    assert meta["min_interlayer_contact_involves_hydrogen"] is True
    assert "stacking_clash" in stacked.flags


def test_min_periodic_pair_contact_finds_cross_boundary_gallery_minimum():
    """Regression: the true minimum can live only across the c boundary.

    The displayed pair is 8.0 A apart (beyond the cutoff); the contact exists
    solely via image (0, 0, -1).  Checking only the displayed layer pair —
    the old behaviour — would miss it.
    """
    cell = ((10.0, 0.0, 0.0), (0.0, 10.0, 0.0), (0.0, 0.0, 10.0))
    atoms_i = [("a_C1", "C", (0.0, 0.0, 1.0))]
    atoms_j = [("b_C1", "C", (0.0, 0.0, 9.0))]

    record = min_periodic_pair_contact(cell, atoms_i, atoms_j, cutoff=3.5)

    assert record is not None
    assert record.distance == pytest.approx(2.0)
    assert record.image == (0, 0, -1)
    assert record.clash
    assert not record.involves_hydrogen


def test_min_periodic_pair_contact_returns_none_as_cutoff_bound():
    cell = ((10.0, 0.0, 0.0), (0.0, 10.0, 0.0), (0.0, 0.0, 10.0))
    atoms_i = [("a_C1", "C", (0.0, 0.0, 1.0))]
    atoms_j = [("b_C1", "C", (0.0, 0.0, 6.0))]

    assert min_periodic_pair_contact(cell, atoms_i, atoms_j, cutoff=3.5) is None


def test_assess_pair_shared_criterion():
    assert assess_pair(2.0, "C", "C").clash  # ratio 0.588 < 0.75
    assert not assess_pair(2.6, "C", "C").clash  # ratio 0.765
    assert assess_pair(1.7, "H", "H").clash  # ratio 0.708 < 0.75
    assert not assess_pair(2.1, "H", "H").clash  # ratio 0.875; the 2.2 A floor is heavy-only
    assert assess_pair(2.1, "H", "H").involves_hydrogen
    assert not assess_pair(2.6, "C", "C").involves_hydrogen


def test_vdw_radius_fallback_warns_once_per_symbol():
    assert vdw_radius("C") == pytest.approx(1.70)
    stderr = io.StringIO()
    with contextlib.redirect_stderr(stderr):
        radius = vdw_radius("Xe")
    assert radius == pytest.approx(1.70)  # documented fallback
    assert "warning: vdw: no Bondi vdW radius for element 'Xe'" in stderr.getvalue()


def test_stacking_contact_failure_warns_and_degrades_without_raising():
    """Per-candidate isolation: a measurement failure degrades the metadata
    and flags instead of aborting the stacking expansion."""
    monomer = _rod_monomer()
    candidate = _stackable_candidate(monomer)

    stderr = io.StringIO()
    with mock.patch("cofkit.stacking.min_periodic_pair_contact", side_effect=ValueError("boom")):
        with contextlib.redirect_stderr(stderr):
            stacked = _apply_layer_registry(
                candidate,
                LayerRegistry(id="AA", lateral_shift=(0.0, 0.0), interlayer_distance=3.4),
                monomer_specs={"m1": monomer},
            )

    meta = stacked.metadata["stacking"]
    assert "stacking_clash" not in stacked.flags
    assert "min_interlayer_contact" not in meta
    assert any("interlayer contact measurement failed" in warning for warning in meta["warnings"])
    assert "interlayer contact measurement failed (ValueError: boom)" in stderr.getvalue()


def test_stacking_geometry_comment_labels_atomistic_contact():
    from cofkit.cif import CIFWriter

    metadata = {
        "stacking": {
            "id": "AA",
            "interlayer_clearance": 3.4,
            "layer_z_span": 3.0,
            "center_to_center_distance": 6.4,
            "layer_count": 2,
            "min_interlayer_contact": 3.4,
            "min_interlayer_contact_mode": "atomistic_product",
            "min_interlayer_contact_atoms": ("i0L0_C1", "i0L1_C4"),
            "min_interlayer_contact_image": (0, 0, -1),
            "min_interlayer_contact_involves_hydrogen": False,
        }
    }

    comment = CIFWriter()._stacking_geometry_comment(metadata)

    assert "min_interlayer_contact=3.4" in comment
    assert "atomistic contact, mode=atomistic_product" in comment
    assert "atoms=i0L0_C1..i0L1_C4" in comment
    assert "image=(0, 0, -1)" in comment
    assert "involves_hydrogen=False" in comment


def test_stacking_geometry_comment_reports_cutoff_bound_when_no_contact():
    from cofkit.cif import CIFWriter

    metadata = {
        "stacking": {
            "id": "AA",
            "min_interlayer_contact": None,
            "min_interlayer_contact_mode": "atomistic_product",
            "min_interlayer_contact_cutoff": 3.5,
        }
    }

    comment = CIFWriter()._stacking_geometry_comment(metadata)

    assert "min_interlayer_contact>3.5" in comment
    assert "none below cutoff" in comment


# ---------------------------------------------------------------------------
# Tilted-c base cells: stacked c is built along the layer normal (route A,
# audit #13 / stacking review point 2)
# ---------------------------------------------------------------------------

# c tilted 1° away from the ab-plane normal (real fitted cells tilt ~0.5°).
_TILTED_CELL = (
    (15.0, 0.0, 0.0),
    (0.0, 15.0, 0.0),
    (8.0 * sin(radians(1.0)), 0.0, 8.0 * cos(radians(1.0))),
)


def _expand_tilted(monomer, cell=_TILTED_CELL, registry_distance=3.4, events=()):
    candidate = _stackable_candidate(monomer, cell=cell, events=events)
    stderr = io.StringIO()
    with contextlib.redirect_stderr(stderr):
        stacked = _apply_layer_registry(
            candidate,
            LayerRegistry(id="AA", lateral_shift=(0.0, 0.0), interlayer_distance=registry_distance),
            monomer_specs={monomer.id: monomer},
        )
    return stacked, stderr.getvalue()


def test_stacking_tilted_base_cell_builds_c_along_layer_normal():
    monomer = _rod_monomer()
    stacked, stderr = _expand_tilted(monomer)

    # a × b of the tilted cell is exactly +z; the stacked c must lie along it
    # with length 2 * c2c (span of the rod along z is 3.0 Å).
    expected_c2c = 3.4 + 3.0
    assert stacked.state.cell[2] == pytest.approx((0.0, 0.0, 2.0 * expected_c2c))
    assert stacked.state.cell[0] == _TILTED_CELL[0]
    assert stacked.state.cell[1] == _TILTED_CELL[1]

    meta = stacked.metadata["stacking"]
    assert meta["c_axis_basis"] == "layer_normal"
    assert meta["c_axis_orthogonalized"] is True
    assert meta["base_c_tilt_degrees"] == pytest.approx(1.0)
    assert "along the layer normal" in meta["derivation"]
    assert any("rebuilt along the layer normal" in warning for warning in meta["warnings"])
    assert "rebuilt along the layer normal" in stderr

    # Layer offsets ride the orthogonalized c; lateral positions are untouched.
    assert stacked.state.monomer_poses["i0L0"].translation == pytest.approx((0.0, 0.0, 0.25 * 2.0 * expected_c2c))
    assert stacked.state.monomer_poses["i0L1"].translation == pytest.approx((0.0, 0.0, 0.75 * 2.0 * expected_c2c))


def test_stacking_tilted_base_cell_realizes_advertised_normal_clearance():
    """The atomistic interlayer contact equals the advertised clearance —
    the derivation c2c = clearance + span is true in the shipped cell."""
    monomer = _rod_monomer()
    stacked, _ = _expand_tilted(monomer)

    meta = stacked.metadata["stacking"]
    assert meta["layer_z_span"] == pytest.approx(3.0)
    assert meta["center_to_center_distance"] == pytest.approx(
        meta["interlayer_clearance"] + meta["layer_z_span"]
    )
    assert meta["min_interlayer_contact"] == pytest.approx(meta["interlayer_clearance"])
    assert "stacking_clash" not in stacked.flags


def test_stacking_orthogonal_base_cell_is_bit_identical_to_prior_construction():
    monomer = _rod_monomer()
    candidate = _stackable_candidate(monomer)  # orthogonal _CELL

    stderr = io.StringIO()
    with contextlib.redirect_stderr(stderr):
        stacked = _apply_layer_registry(
            candidate,
            LayerRegistry(id="AA", lateral_shift=(0.0, 0.0), interlayer_distance=3.4),
            monomer_specs={"m1": monomer},
        )

    # Identical floats to the legacy scale(c_hat, 2*c2c) construction.
    assert stacked.state.cell[2] == (0.0, 0.0, 2.0 * (3.4 + 3.0))
    assert stacked.state.monomer_poses["i0L0"].translation == (0.0, 0.0, 3.2)

    meta = stacked.metadata["stacking"]
    assert meta["c_axis_basis"] == "layer_normal"
    assert meta["c_axis_orthogonalized"] is False
    assert meta["base_c_tilt_degrees"] == pytest.approx(0.0, abs=1e-9)
    assert not any("rebuilt along the layer normal" in warning for warning in meta["warnings"])
    assert "rebuilt along the layer normal" not in stderr.getvalue()


def test_stacking_degenerate_ab_vectors_fall_back_to_c_direction_with_warning():
    monomer = _rod_monomer()
    degenerate_cell = ((15.0, 0.0, 0.0), (30.0, 0.0, 0.0), (1.0, 0.0, 8.0))  # collinear a, b
    stacked, stderr = _expand_tilted(monomer, cell=degenerate_cell)

    # Prior behaviour preserved: new c is the base c direction scaled to 2*c2c.
    new_c = stacked.state.cell[2]
    base_c = degenerate_cell[2]
    assert norm(cross(new_c, base_c)) == pytest.approx(0.0, abs=1e-9)
    assert dot(new_c, base_c) > 0.0
    assert norm(new_c) == pytest.approx(2.0 * stacked.metadata["stacking"]["center_to_center_distance"])

    meta = stacked.metadata["stacking"]
    assert meta["c_axis_basis"] == "c_hat"
    assert meta["c_axis_orthogonalized"] is False
    assert "base_c_tilt_degrees" not in meta
    assert "degenerate in-plane fallback" in meta["derivation"]
    assert "degenerate in-plane cell vectors" in stderr


def test_stacking_ring_center_nonzero_w_still_warns_when_basis_is_orthogonalized():
    monomer = _rod_monomer()
    events = (
        ReactionEvent(
            id="r1",
            template_id="boroxine_ring",
            participants=(MotifRef("i0", "m1", "mo1"),),
            metadata={"ring_center_fractional": (0.5, 0.5, 0.1)},
        ),
    )
    stacked, stderr = _expand_tilted(monomer, events=events)

    assert "ring_center_fractional w=1.000e-01" in stderr
    assert any(
        "ring_center_fractional w=1.000e-01" in warning
        for warning in stacked.metadata["stacking"]["warnings"]
    )
    # The cartesian offset rides the orthogonalized c for both layers.
    new_c = stacked.state.cell[2]
    offsets = {
        event.id: event.metadata["ring_center_offset_cartesian"]
        for event in stacked.events
    }
    assert offsets["r1L0"] == pytest.approx(tuple(0.25 * component for component in new_c))
    assert offsets["r1L1"] == pytest.approx(tuple(0.75 * component for component in new_c))


def test_stacking_geometry_comment_reports_c_axis_orthogonalization():
    from cofkit.cif import CIFWriter

    metadata = {
        "stacking": {
            "id": "AA",
            "c_axis_basis": "layer_normal",
            "c_axis_orthogonalized": True,
            "base_c_tilt_degrees": 0.525,
        }
    }

    comment = CIFWriter()._stacking_geometry_comment(metadata)

    assert "c_rebuilt_along_layer_normal(base_c_tilt=0.525deg)" in comment


# ---------------------------------------------------------------------------
# W3.2 wiring: the stacking metadata block must reach the exported CIF
# ---------------------------------------------------------------------------


def test_stacked_export_renders_stacking_geometry_comment_from_real_metadata():
    """A stacked export carries the stacking-geometry derivation comment with
    the real registry/shift/clearance/contact values, alongside the
    periodic_bilayer c-axis-semantics label.  Previously the export metadata
    never carried the stacking block and the comment was dead in production."""
    from cofkit.cif import CIFWriter

    monomer = _rod_monomer()
    candidate = _stackable_candidate(monomer)
    (stacked,) = enumerate_candidate_stackings(
        candidate, registry_ids=("AA",), monomer_specs={"m1": monomer}
    )

    export = CIFWriter().export_candidate(stacked, {"m1": monomer})

    comment_lines = [line for line in export.text.splitlines() if line.startswith("# stacking-geometry:")]
    assert len(comment_lines) == 1
    comment = comment_lines[0]
    meta = stacked.metadata["stacking"]
    assert "registry=AA" in comment
    assert "shift_frac=0,0" in comment
    assert "interlayer_clearance=3.4" in comment
    assert f"layer_z_span={meta['layer_z_span']:.4g}" in comment
    assert f"c2c={meta['center_to_center_distance']:.4g}" in comment
    assert "c=2*c2c" in comment
    assert "min_interlayer_contact=3.4" in comment
    assert "atomistic contact, mode=precursor_coordinates" in comment
    assert "atoms=" in comment
    assert "image=" in comment
    assert "involves_hydrogen=" in comment

    # Coexists with the c-axis-semantics label (W1.4).
    assert "# c-axis-semantics: periodic_bilayer" in export.text
    assert export.metadata["stacking"]["id"] == "AA"


def test_monolayer_export_omits_stacking_geometry_comment():
    from cofkit.cif import CIFWriter

    monomer = _rod_monomer()
    candidate = _stackable_candidate(monomer)

    export = CIFWriter().export_candidate(candidate, {"m1": monomer})

    assert "# stacking-geometry:" not in export.text
    assert "stacking" not in export.metadata


def test_stacked_export_comment_reports_c_axis_orthogonalization_from_real_metadata():
    from cofkit.cif import CIFWriter

    monomer = _rod_monomer()
    stacked, _ = _expand_tilted(monomer)

    export = CIFWriter().export_candidate(stacked, {"m1": monomer})

    comment = next(
        line for line in export.text.splitlines() if line.startswith("# stacking-geometry:")
    )
    assert "c_rebuilt_along_layer_normal(base_c_tilt=1deg)" in comment
