import math
from typing import Mapping

from ..geometry import Frame, normalize, vec3
from ..model import MonomerSpec, ReactionTemplate


class LinkageGeometryModel:
    """Computes initial guesses for linkage geometry."""

    @staticmethod
    def initialize_imine_bridge(
        amine_frame: Frame, aldehyde_frame: Frame
    ) -> dict[str, float]:
        """Provides target parameters for an imine bridge."""
        return {
            "target_bond_length": 1.28,  # C=N double bond
            "target_c_n_c_angle_deg": 120.0,
            "planarity_weight": 1.0,
        }

    @staticmethod
    def initialize_boroxine_ring(
        b1_frame: Frame, b2_frame: Frame, b3_frame: Frame
    ) -> dict[str, float]:
        """Provides target parameters for a boroxine ring."""
        return {
            "target_b_o_bond_length": 1.38,
            "target_b_o_b_angle_deg": 120.0,
            "target_o_b_o_angle_deg": 120.0,
            "planarity_weight": 2.0,
        }


def orient_frames_for_bridge(f1: Frame, f2: Frame, target_distance: float) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
    """Given two frames, computes the translation to align them head-to-head.
    
    This is a simplistic geometric initialization that assumes f1 stays at origin
    and f2 is moved/rotated to face it exactly at `target_distance`.
    """
    # For a perfect head-to-head bridge:
    # f1.primary should point exactly opposite to f2.primary
    # f1.normal should align with f2.normal (if we want planar conjugation)
    
    # Place f2 origin at f1.origin + f1.primary * target_distance
    t2 = (
        f1.origin[0] + f1.primary[0] * target_distance,
        f1.origin[1] + f1.primary[1] * target_distance,
        f1.origin[2] + f1.primary[2] * target_distance,
    )
    
    # In a full optimizer, we'd return rigid body transforms.
    # Here we just return target translations to prove the abstraction.
    t1 = (0.0, 0.0, 0.0)
    
    return t1, t2
