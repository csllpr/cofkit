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

import pytest

from cofkit.model import (
    AssemblyState,
    Candidate,
    MotifRef,
    ReactionEvent,
    candidate_ranking_key,
    order_candidates,
)
from cofkit.stacking import enumerate_candidate_stackings

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
