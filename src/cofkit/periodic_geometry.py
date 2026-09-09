"""Cell widths and image-aware distances for periodic calculations."""

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
