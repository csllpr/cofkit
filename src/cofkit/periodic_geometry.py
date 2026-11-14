"""Cell widths and image-aware distances for periodic calculations."""

from collections.abc import Iterator, Sequence
from itertools import product
import math

import gemmi


def face_widths(cell) -> tuple[float, float, float]:
    """Perpendicular widths in Å, from reciprocal vectors (without 2π)."""
    if not math.isfinite(cell.volume) or cell.volume <= 0:
        raise ValueError("Unit cell must have positive finite volume.")
    widths = tuple(
        1 / math.sqrt(sum(x * x for x in row)) for row in cell.frac.mat.tolist()
    )
    if not all(math.isfinite(v) and v > 0 for v in widths):
        raise ValueError("Unit cell has invalid perpendicular widths.")
    return widths


def images_within(cell, reference, other, cutoff):
    """Enumerate all translations inside a sphere, including in small cells.

    Reciprocal widths bound the fractional components of any Cartesian vector
    inside the sphere. This avoids a fixed ±1 image-search assumption.
    """
    delta = tuple(b - a for a, b in zip(reference, other))
    bounds = [
        range(math.ceil(-d - cutoff / w), math.floor(-d + cutoff / w) + 1)
        for d, w in zip(delta, face_widths(cell))
    ]
    for shift in product(*bounds):
        vector = gemmi.Fractional(*(d + n for d, n in zip(delta, shift)))
        distance = cell.orthogonalize(vector).length()
        if distance <= cutoff:
            yield shift, distance


def p1_shift(token: str) -> tuple[int, int, int]:
    if token in ("", ".", "?", "1", "1_555"):
        return (0, 0, 0)
    if not token.startswith("1_") or len(token) != 5 or not token[2:].isdigit():
        raise ValueError(f"Expected P1 bond-image translation, got {token!r}.")
    return tuple(int(v) - 5 for v in token[2:])


def format_p1_shift(shift: tuple[int, int, int]) -> str:
    """Encode a translation without silently emitting an invalid CIF code."""
    if shift == (0, 0, 0):
        return "."
    # Derived: each CIF translation digit encodes shift + 5 in one digit.
    digits = tuple(component + 5 for component in shift)
    if any(not isinstance(digit, int) or not 0 <= digit <= 9 for digit in digits):
        raise ValueError(f"periodic image shift {shift!r} exceeds CIF P1 formatter range")
    return "1_" + "".join(str(digit) for digit in digits)


def rebase_bond_shift(shift, first_fractional, second_fractional) -> tuple[int, int, int]:
    """Preserve an edge when its unwrapped endpoint coordinates are wrapped.

    Derived: f'_i = f_i - floor(f_i), so the vector f_j + g - f_i
    becomes f'_j + [g + floor(f_j) - floor(f_i)] - f'_i.
    This is a change of atom representatives, never a nearest-image search.
    """
    return tuple(
        gain + math.floor(second) - math.floor(first)
        for gain, first, second in zip(shift, first_fractional, second_fractional)
    )


def nearby_atom_pairs(
    fractional_positions: Sequence[Sequence[float]], widths: Sequence[float],
    cutoff: float, *, include_self: bool = False,
) -> Iterator[tuple[int, int]]:
    """Conservative periodic cell-list candidates in the original i,j order.

    Derived: a Cartesian separation <= cutoff has fractional component
    |delta_k| <= cutoff / face_width_k. Choose fractional bins at least that
    wide, so only the same and adjacent bins (with periodic wrapping) can
    contain neighbors. Reciprocal face widths make this true for skew cells
    too. Exact distances, image ties and exclusions remain downstream.
    """
    n = len(fractional_positions)
    if not n:
        return
    if cutoff <= 0 or not math.isfinite(cutoff):
        # Preserve the old exact-search behavior for unusual caller settings;
        # no spatial subdivision is defined for a nonpositive/infinite radius.
        for i in range(n):
            for j in range(i if include_self else i + 1, n):
                yield i, j
        return
    # Derived: the exact search subtracts two original fractional endpoints.
    # Allow one ULP per endpoint before choosing bins, since large unwrapped
    # image coordinates can lose more precision than a wrapped coordinate or
    # the cutoff alone. This only broadens candidates, never the exact cutoff.
    # More than n bins along an axis is unnecessary for this sparse index.
    counts = []
    for axis, width in enumerate(widths):
        rounding = 2 * max(math.ulp(f[axis]) for f in fractional_positions)
        fractional_cutoff = math.nextafter(cutoff / width + rounding, math.inf)
        counts.append(max(1, math.floor(min(n, 1 / fractional_cutoff))))
    # Derived: the three adjacent indices cover every bin on an axis with
    # <=3 bins. If this holds on all axes, indexing cannot prune any pair;
    # keep the cheaper ordered loop for small/dense cells.
    if all(count <= 3 for count in counts):
        for i in range(n):
            for j in range(i if include_self else i + 1, n):
                yield i, j
        return
    bins: dict[tuple[int, int, int], list[int]] = {}
    keys = []
    for i, fractional in enumerate(fractional_positions):
        key = tuple(
            math.floor((coordinate % 1.0) * count) % count
            for coordinate, count in zip(fractional, counts)
        )
        keys.append(key)
        bins.setdefault(key, []).append(i)
    for i, key in enumerate(keys):
        # Sets collapse repeated bins when an axis has only one or two bins.
        adjacent = (
            {(coordinate + offset) % count for offset in (-1, 0, 1)}
            for coordinate, count in zip(key, counts)
        )
        neighbors = {
            j for neighbor_key in product(*adjacent)
            for j in bins.get(neighbor_key, ()) if j > i or (include_self and j == i)
        }
        # Self pairs are optional; exact image selection belongs to the caller.
        for j in sorted(neighbors):
            yield i, j
