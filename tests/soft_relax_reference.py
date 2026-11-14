"""Frozen exhaustive soft-relax pair search for independent regression checks.

Copied from the pre-spatial-search implementation. Keep the all-pairs loop:
this is an oracle for candidate completeness, ordering, and image selection,
not a production backend. Physical validity is a separate question.
"""

from cofkit.periodic_geometry import images_within
from cofkit.soft_relax import SoftRelaxConfig, _System, _element_radius, _frac_to_cart


def exhaustive_pair_list(
    system: _System, config: SoftRelaxConfig
) -> list[tuple[int, int, list[float], float]]:
    radii = [_element_radius(symbol, config) for symbol in system.symbols]
    if not config.repel_hydrogen:
        radii = [0.0 if s == "H" else r for s, r in zip(system.symbols, radii)]
    else:
        radii = [
            config.hydrogen_radius if s == "H" else r
            for s, r in zip(system.symbols, radii)
        ]
    max_r0 = 2.0 * max(radii, default=config.fallback_radius) * config.repulsion_scale
    cutoff = max_r0 + config.neighbor_margin
    hetero_elements = set(config.hbond_elements)
    hb_lo, hb_hi = config.hbond_window

    def hbond_guarded(i: int, j: int, distance: float) -> bool:
        if not config.hbond_guard or not (hb_lo <= distance <= hb_hi):
            return False
        if system.symbols[i] == "H":
            _h, other = i, j
        elif system.symbols[j] == "H":
            _h, other = j, i
        else:
            return False
        return system.symbols[other] in hetero_elements

    pairs: list[tuple[int, int, list[float], float]] = []
    n = len(system.frac)
    for i in range(n):
        for j in range(i, n):
            r0 = (radii[i] + radii[j]) * config.repulsion_scale
            if r0 <= 1e-9:
                continue
            # images_within yields in lexicographic shift order, not distance
            # order: select the minimum-distance image first so the clash-floor
            # and hbond-guard decisions below are evaluated at the nearest
            # image, not at an arbitrary in-cutoff one.
            nearest: tuple[float, tuple[int, int, int]] | None = None
            for shift, distance in images_within(system.cell, system.frac[i], system.frac[j], cutoff):
                if i == j and shift == (0, 0, 0):
                    continue
                if nearest is None or distance < nearest[0]:
                    nearest = (distance, shift)
            if nearest is None:
                continue
            distance, shift = nearest
            if (i, j, shift) in system.excluded:
                continue
            if (i, j, shift) in system.excluded13 and distance >= config.one_four_clash_floor:
                continue
            if (i, j, shift) in system.excluded14 and distance >= config.one_four_clash_floor:
                continue
            if hbond_guarded(i, j, distance):
                continue
            shift_cart = _frac_to_cart(system.orth, shift)
            # nearest image only; the margin covers drift during the run
            pairs.append((i, j, shift_cart, r0))
    return pairs
