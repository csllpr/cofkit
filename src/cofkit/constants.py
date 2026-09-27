"""Shared, single-owner constants for cell defaults.

This module must import nothing from cofkit so any module may depend on it
without creating an import cycle.
"""

# Vacuum-slab padding along c for single-layer (monolayer) cells, chosen so a
# single layer is effectively non-periodic along c under current calculation
# practices. Heuristic — pending calibration; value unified 2026-09 per W1.4
# of the stacking plan. NOT an interlayer distance; stacking registries carry
# their own clearances.
DEFAULT_MONOLAYER_C_ANGSTROM = 8.0

# Default lateral span of the square placeholder cell built before the cell is
# fitted to the assembled structure. Heuristic — pending calibration.
DEFAULT_LATERAL_SPAN_ANGSTROM = 30.0
