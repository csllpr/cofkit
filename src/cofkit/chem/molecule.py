from __future__ import annotations

from dataclasses import dataclass
from typing import Any


@dataclass(frozen=True)
class Molecule:
    """A minimal dependency-free representation of a molecule."""

    symbols: tuple[str, ...]
    positions: tuple[tuple[float, float, float], ...]
    bonds: tuple[tuple[int, int, float], ...] = ()

    @classmethod
    def from_xyz(cls, text: str) -> "Molecule":
        lines = [line.strip() for line in text.split("\n") if line.strip()]
        if not lines:
            raise ValueError("empty XYZ text")
        
        n_atoms = int(lines[0])
        if len(lines) < n_atoms + 2:
            raise ValueError("XYZ text truncated")
            
        symbols = []
        positions = []
        for line in lines[2:2+n_atoms]:
            parts = line.split()
            if len(parts) >= 4:
                symbols.append(parts[0])
                positions.append((float(parts[1]), float(parts[2]), float(parts[3])))
                
        return cls(symbols=tuple(symbols), positions=tuple(positions))

    def __len__(self) -> int:
        return len(self.symbols)
