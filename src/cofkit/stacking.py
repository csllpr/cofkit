from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class LayerRegistry:
    id: str
    lateral_shift: tuple[float, float] = (0.0, 0.0)
    interlayer_distance: float = 3.4


class StackingExplorer:
    """Placeholder enumerator for 2D COF layer registries.

    The real implementation will score registry choices against steric clashes,
    pore alignment, and approximate pi-stacking plausibility.
    """

    def enumerate_registries(self) -> tuple[LayerRegistry, ...]:
        return (
            LayerRegistry(id="AA", lateral_shift=(0.0, 0.0), interlayer_distance=3.4),
            LayerRegistry(id="AB", lateral_shift=(0.5, 0.5), interlayer_distance=3.5),
            LayerRegistry(id="slipped", lateral_shift=(1.2, 0.6), interlayer_distance=3.6),
        )
