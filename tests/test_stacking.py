"""Tests for stacking.py: interlayer-offset computation, stacking-pattern resolution."""

import pytest
from cofkit.stacking import (
    resolve_stacking_pattern,
    compute_interlayer_offset,
)


def test_resolve_stacking_pattern_explicit():
    assert resolve_stacking_pattern("AA", None, None) == "AA"
    assert resolve_stacking_pattern("AB", None, None) == "AB"
    assert resolve_stacking_pattern("ABC", None, None) == "ABC"


def test_resolve_stacking_pattern_default_for_topology():
    """M2.1: sql → AA, kgm → AB, others → AA."""
    assert resolve_stacking_pattern(None, "sql", None) == "AA"
    assert resolve_stacking_pattern(None, "kgm", None) == "AB"
    assert resolve_stacking_pattern(None, "hcb", None) == "AA"
    assert resolve_stacking_pattern(None, "fxt", None) == "AA"


def test_resolve_stacking_pattern_from_assembly_mode():
    """M2.2: assembly_mode='serrated' → AB, others → AA."""
    assert resolve_stacking_pattern(None, "hcb", "serrated") == "AB"
    assert resolve_stacking_pattern(None, "sql", "serrated") == "AB"
    assert resolve_stacking_pattern(None, "hcb", "default") == "AA"
    assert resolve_stacking_pattern(None, None, "serrated") == "AB"


def test_resolve_stacking_pattern_explicit_wins():
    """M2.3: explicit pattern overrides topology and assembly_mode."""
    assert resolve_stacking_pattern("ABC", "kgm", "serrated") == "ABC"
    assert resolve_stacking_pattern("AA", "kgm", None) == "AA"


def test_compute_interlayer_offset_aa_zero():
    """M3.1: AA → zero offset."""
    assert compute_interlayer_offset("AA", 10.0, 15.0) == 0.0


def test_compute_interlayer_offset_ab_half_diagonal():
    """M3.2: AB → half the in-plane diagonal."""
    a, b = 10.0, 15.0
    expected = 0.5 * (a**2 + b**2) ** 0.5
    assert compute_interlayer_offset("AB", a, b) == pytest.approx(expected, rel=1e-9)


def test_compute_interlayer_offset_abc_third_diagonal():
    """M3.3: ABC → one-third the diagonal."""
    a, b = 12.0, 16.0
    expected = (1.0 / 3.0) * (a**2 + b**2) ** 0.5
    assert compute_interlayer_offset("ABC", a, b) == pytest.approx(expected, rel=1e-9)


def test_compute_interlayer_offset_unknown_fallback():
    """M3.4: unknown pattern → zero (defensive fallback)."""
    assert compute_interlayer_offset("XYZ", 10.0, 15.0) == 0.0
    assert compute_interlayer_offset("", 10.0, 15.0) == 0.0
