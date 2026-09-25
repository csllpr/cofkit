"""Shared van der Waals contact measurement for validation and stacking.

This module is the single source of truth for the sub-vdW clash criterion
used by both the coarse structure validator (``validation.py``) and the
stacking self-check (``stacking.py``):

* Bondi vdW radii table with an explicit supported-element set and a
  documented fallback policy (``vdw_radius``).
* The shared pair criterion (``assess_pair``): a pair is a clash when
  ``d / (r_vdw_i + r_vdw_j) < DEFAULT_VDW_CLASH_RATIO`` or — for
  heavy-heavy pairs only — when the plain distance falls below
  ``DEFAULT_HARD_MIN_NONBONDED_PLAIN_DISTANCE`` (a severe-overlap floor that
  does not depend on the radius table).  Hydrogen-involving pairs use the
  ratio criterion alone; the heavy plain-distance floor would false-positive
  ordinary pore-facing and stacked C-H / H...H contacts, whose vdW sums are
  smaller.
* ``min_periodic_pair_contact``: the minimum-distance pair between two atom
  sets under all periodic images within a cutoff (via
  ``periodic_geometry.images_within``, which enumerates every lattice
  translation inside the cutoff sphere — including both ``(0, 0, ±1)``
  galleries across a c boundary, not just the displayed pair).

The radii here are Bondi vdW radii; they are deliberately NOT the DREIDING
force-field radii from ``_dreiding_reference.py`` — force-field LJ radii are
not interchangeable with vdW contact radii (soft-relax keeps using its own
DREIDING-derived repulsion radii).

gemmi is imported lazily inside ``min_periodic_pair_contact`` so that the
radii/criterion half of the module stays importable in gemmi-less
environments.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass
from math import acos, degrees
from typing import Iterable, Mapping

from .geometry import Vec3, dot, norm

# Bondi vdW radii (A): Bondi, J. Phys. Chem. 1964, 68, 441, extended by his
# 1966 compilation for B, Br, I.
BONDI_VDW_RADII: Mapping[str, float] = {
    "H": 1.20,
    "B": 1.92,
    "C": 1.70,
    "N": 1.55,
    "O": 1.52,
    "F": 1.47,
    "Si": 2.10,
    "P": 1.80,
    "S": 1.80,
    "Cl": 1.75,
    "Br": 1.85,
    "I": 1.98,
}

SUPPORTED_VDW_ELEMENTS: frozenset[str] = frozenset(BONDI_VDW_RADII)

# Fallback policy for unsupported elements: warn on stderr (once per symbol)
# and use a carbon-like midpoint radius.  1.70 A keeps the ratio criterion
# meaningful for organic-adjacent elements without silently dropping the
# check; callers that need exact radii for exotic elements should extend the
# table above instead.
FALLBACK_VDW_RADIUS: float = 1.70

# Shared clash criterion defaults.  CoarseValidationThresholds defaults its
# corresponding fields to these values, and stacking uses them directly, so
# both consumers flag the same contacts.
DEFAULT_VDW_CLASH_RATIO: float = 0.75
DEFAULT_HARD_MIN_NONBONDED_PLAIN_DISTANCE: float = 2.2
# Search radius: the largest flaggable pair distance for supported elements
# is 0.75 * 2 * 2.10 = 3.15 A (Si...Si); 3.5 A covers it with margin.
DEFAULT_NONBONDED_SEARCH_RADIUS: float = 3.5
# Floor for directly bonded (1-2) pairs: below this the "bond" is fused
# nuclei — a broken structure.  Must stay below genuine heavy-atom bonds
# (~1.16 A for C#N) and X-H bonds (~1.0 A).
DEFAULT_SEVERE_OVERLAP_BONDED_DISTANCE: float = 0.7

_HYDROGEN_SYMBOLS = frozenset({"H", "D"})

_warned_fallback_symbols: set[str] = set()


def vdw_radius(symbol: str) -> float:
    """Bondi vdW radius (A) for *symbol*; warn + fallback for unsupported ones."""
    radius = BONDI_VDW_RADII.get(symbol)
    if radius is not None:
        return radius
    if symbol not in _warned_fallback_symbols:
        _warned_fallback_symbols.add(symbol)
        print(
            f"warning: vdw: no Bondi vdW radius for element {symbol!r}; "
            f"using the fallback radius {FALLBACK_VDW_RADIUS:.2f} A",
            file=sys.stderr,
        )
    return FALLBACK_VDW_RADIUS


def is_hydrogen_symbol(symbol: str) -> bool:
    return symbol in _HYDROGEN_SYMBOLS


@dataclass(frozen=True)
class PairAssessment:
    """Result of applying the shared sub-vdW clash criterion to one pair."""

    distance: float
    vdw_ratio: float
    involves_hydrogen: bool
    clash: bool


def assess_pair(
    distance: float,
    symbol_i: str,
    symbol_j: str,
    *,
    ratio_threshold: float = DEFAULT_VDW_CLASH_RATIO,
    plain_floor: float = DEFAULT_HARD_MIN_NONBONDED_PLAIN_DISTANCE,
) -> PairAssessment:
    """Apply the shared clash criterion to one pair at *distance*.

    Clash when ``distance / (r_i + r_j) < ratio_threshold`` or — heavy-heavy
    pairs only — ``distance < plain_floor``.  Hydrogen-involving pairs skip
    the plain floor (see module docstring).
    """
    involves_hydrogen = is_hydrogen_symbol(symbol_i) or is_hydrogen_symbol(symbol_j)
    ratio = distance / (vdw_radius(symbol_i) + vdw_radius(symbol_j))
    clash = ratio < ratio_threshold or (not involves_hydrogen and distance < plain_floor)
    return PairAssessment(
        distance=distance,
        vdw_ratio=ratio,
        involves_hydrogen=involves_hydrogen,
        clash=clash,
    )


@dataclass(frozen=True)
class ContactRecord:
    """One minimum-distance periodic contact between two atom sets.

    ``image`` is the integer lattice translation applied to the second atom's
    home-cell position, e.g. ``(0, 0, -1)`` for a contact across the -c cell
    boundary.  ``clash`` applies the shared criterion (``assess_pair``) with
    the thresholds passed to ``min_periodic_pair_contact``.
    """

    distance: float
    vdw_ratio: float
    label_i: str
    label_j: str
    symbol_i: str
    symbol_j: str
    image: tuple[int, int, int]
    involves_hydrogen: bool
    clash: bool


def _unit_cell_from_vectors(cell_vectors: tuple[Vec3, Vec3, Vec3]):
    """Build a gemmi.UnitCell from Cartesian cell vectors."""
    import gemmi

    a_vec, b_vec, c_vec = cell_vectors
    a_len, b_len, c_len = norm(a_vec), norm(b_vec), norm(c_vec)
    if min(a_len, b_len, c_len) < 1e-8:
        raise ValueError("degenerate cell vector in contact measurement")
    alpha = degrees(acos(max(-1.0, min(1.0, dot(b_vec, c_vec) / (b_len * c_len)))))
    beta = degrees(acos(max(-1.0, min(1.0, dot(a_vec, c_vec) / (a_len * c_len)))))
    gamma = degrees(acos(max(-1.0, min(1.0, dot(a_vec, b_vec) / (a_len * b_len)))))
    return gemmi.UnitCell(a_len, b_len, c_len, alpha, beta, gamma)


def min_periodic_pair_contact(
    cell_vectors: tuple[Vec3, Vec3, Vec3],
    atoms_i: Iterable[tuple[str, str, Vec3]],
    atoms_j: Iterable[tuple[str, str, Vec3]],
    *,
    cutoff: float = DEFAULT_NONBONDED_SEARCH_RADIUS,
    ratio_threshold: float = DEFAULT_VDW_CLASH_RATIO,
    plain_floor: float = DEFAULT_HARD_MIN_NONBONDED_PLAIN_DISTANCE,
) -> ContactRecord | None:
    """Minimum-distance pair between two Cartesian atom sets, periodic-aware.

    ``atoms_*`` are ``(label, element_symbol, (x, y, z))`` triples in the
    Cartesian frame of *cell_vectors*.  Every lattice translation inside the
    cutoff sphere is enumerated (``periodic_geometry.images_within``), so
    contacts across cell boundaries — both ``(0, 0, ±1)`` galleries for a
    layered cell — compete on equal footing with the displayed pair.

    Returns ``None`` when no pair lies within *cutoff*; since *cutoff* is
    sized above the largest flaggable contact distance, ``None`` is a
    ``> cutoff`` bound statement, not a missing measurement.
    """
    import gemmi

    from .periodic_geometry import images_within

    cell = _unit_cell_from_vectors(cell_vectors)
    atoms_j_list = [
        (label, symbol, cell.fractionalize(gemmi.Position(*position)))
        for label, symbol, position in atoms_j
    ]
    best: tuple[float, tuple[int, int, int], str, str, str, str] | None = None
    for label_i, symbol_i, position_i in atoms_i:
        fract_i = cell.fractionalize(gemmi.Position(*position_i))
        reference = (fract_i.x, fract_i.y, fract_i.z)
        for label_j, symbol_j, fract_j in atoms_j_list:
            other = (fract_j.x, fract_j.y, fract_j.z)
            for shift, distance in images_within(cell, reference, other, cutoff):
                if best is None or distance < best[0]:
                    best = (distance, shift, label_i, label_j, symbol_i, symbol_j)
    if best is None:
        return None
    distance, image, label_i, label_j, symbol_i, symbol_j = best
    assessment = assess_pair(
        distance,
        symbol_i,
        symbol_j,
        ratio_threshold=ratio_threshold,
        plain_floor=plain_floor,
    )
    return ContactRecord(
        distance=distance,
        vdw_ratio=assessment.vdw_ratio,
        label_i=label_i,
        label_j=label_j,
        symbol_i=symbol_i,
        symbol_j=symbol_j,
        image=image,
        involves_hydrogen=assessment.involves_hydrogen,
        clash=assessment.clash,
    )


__all__ = [
    "BONDI_VDW_RADII",
    "ContactRecord",
    "DEFAULT_HARD_MIN_NONBONDED_PLAIN_DISTANCE",
    "DEFAULT_NONBONDED_SEARCH_RADIUS",
    "DEFAULT_SEVERE_OVERLAP_BONDED_DISTANCE",
    "DEFAULT_VDW_CLASH_RATIO",
    "FALLBACK_VDW_RADIUS",
    "PairAssessment",
    "SUPPORTED_VDW_ELEMENTS",
    "assess_pair",
    "is_hydrogen_symbol",
    "min_periodic_pair_contact",
    "vdw_radius",
]
