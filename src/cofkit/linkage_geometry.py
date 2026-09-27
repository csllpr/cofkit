from __future__ import annotations

from math import atan2, cos, radians, sin, sqrt

from .geometry import Vec3, norm, scale, sub
from .model import MonomerSpec
from .reactions import bridge_geometry_priors, bridge_target_distance

# Single-owner imine/azine bridge geometry (2026-09-27 rework).
#
# Two compensators used to fight imine collinearity independently: a
# hand-tuned embedding-time motif-origin retraction (0.11 for imine, 0.08
# for azine) and a realization-time 721-step angular scan whose objective
# was displacement-dominated (exported angles landed near 131 degrees
# regardless of the angle target). Both are now derived from the template
# profile's geometry priors (BridgeGeometryPriors in reactions.py):
#
# * required_bridge_span() below solves the anchor-to-anchor distance at
#   which the bridge chain A-C-N-B closes exactly with both interior
#   angles at their priors, in the E (zig-zag) configuration with C and N
#   on opposite sides of the anchor axis.
# * The realization-time constructor (reaction_realization.py) solves the
#   same quadrilateral: exactly when the placed anchor distance matches
#   the required span, best-effort (least-squares on the two angle
#   residuals, lengths held exact) with an honest diagnostic note
#   otherwise.
# * The embedding-time retraction below anticipates that constructor:
#   under ideal aligned placement (effective origins at
#   bridge_target_distance separation, anchor->reactive bonds collinear
#   with the bridge axis) the placed anchor distance is
#
#       d(f) = bridge_target + (1 - f) * (l_C + l_N)
#
#   for a retraction fraction f applied on both sides. Setting
#   d(f) = required_bridge_span gives
#
#       f = (bridge_target + l_C + l_N - d*) / (l_C + l_N).
#
#   Each motif computes f independently with the symmetric-event
#   approximation l_C = l_N = l_own (its own measured |anchor->reactive|);
#   for typical aryl side lengths (1.40-1.45 angstrom on both sides) the
#   resulting total retraction matches the joint optimum to ~1 percent.
#
#   Azine is the asymmetric case: its realization fit pins both nitrogens
#   on the N-N axis, so the constructed endpoint is the triangle
#   A-C-N with |A-C| = l_C, |C-N| = bridge_target and the third side set
#   by placement. Consistency means |A-N| equals the triangle-closure
#   distance r* = sqrt(l_C^2 + l_CN^2 - 2 l_C l_CN cos theta_C); the C-N-N
#   angle is then a derived consequence, not independently fittable. The
#   retraction solves the same symmetric-split equation against r*:
#
#       f = (bridge_target + l_own - r*) / (2 * l_own).
#
# Derived values for representative aryl side lengths, versus the retired
# hand-tuned constants: imine l=1.45 -> f = 0.164 (was 0.11); azine
# l=1.45 -> f = 0.127 (was 0.08). The hand-tuned values systematically
# under-retract: they were tuned against the old displacement-dominated
# scan, which never actually reached the angle priors, so full prior
# anticipation was never wanted. The derived values are larger because
# they place the anchors where the constructor can actually close at the
# priors.

# Anchor distances within this tolerance of required_bridge_span are
# closed exactly at the angle priors; larger deviations take the
# best-effort path with a deviation diagnostic.
BRIDGE_SPAN_EXACT_TOLERANCE_ANGSTROM = 0.02


def required_bridge_span(
    outer_length_a: float,
    bridge_length: float,
    outer_length_b: float,
    angle_a_deg: float,
    angle_b_deg: float,
) -> float | None:
    """Anchor distance that makes the bridge quadrilateral close exactly.

    Chain A-X-Y-B with |A-X| = outer_length_a, |X-Y| = bridge_length,
    |Y-B| = outer_length_b and interior angles angle_a_deg at X, angle_b_deg
    at Y, in the E (zig-zag) configuration: A=(0,0), B=(d,0), X above the
    axis, Y below. Writing the exterior angles tA = 180 - angle_a_deg and
    tB = 180 - angle_b_deg, the segment directions are a, a - tA and
    a - tA + tB, and lateral closure is linear in sin a / cos a:

        S sin a - D cos a = 0
        S = outer_a + bridge cos tA + outer_b cos(tA - tB)
        D = bridge sin tA + outer_b sin(tA - tB)

    so a = atan2(D, S) and
    d = outer_a cos a + bridge cos(a - tA) + outer_b cos(a - tA + tB).

    Returns None when no zig-zag solution exists (the solved chain fails
    to put X and Y on opposite sides of the axis, or d <= 0).
    """
    lengths = (outer_length_a, bridge_length, outer_length_b)
    if any(length < 1e-8 for length in lengths):
        return None
    exterior_a = radians(180.0 - angle_a_deg)
    exterior_b = radians(180.0 - angle_b_deg)
    sin_a = sin(exterior_a)
    if sin_a <= 1e-12:
        return None
    s_coeff = outer_length_a + bridge_length * cos(exterior_a) + outer_length_b * cos(exterior_a - exterior_b)
    d_coeff = bridge_length * sin_a + outer_length_b * sin(exterior_a - exterior_b)
    if d_coeff <= 0.0:
        return None
    alpha = atan2(d_coeff, s_coeff)
    if not (1e-9 < alpha < exterior_a - 1e-9):
        # X must sit above the axis and the X->Y segment must head downward.
        return None
    x_height = outer_length_a * sin(alpha)
    y_height = x_height + bridge_length * sin(alpha - exterior_a)
    if y_height >= 0.0:
        return None
    span = (
        outer_length_a * cos(alpha)
        + bridge_length * cos(alpha - exterior_a)
        + outer_length_b * cos(alpha - exterior_a + exterior_b)
    )
    if span <= 0.0:
        return None
    return span


def derived_origin_retraction_fraction(template_id: str, side_length: float) -> float:
    """Embedding-time motif-origin retraction fraction derived from the
    template profile's geometry priors (see the module docstring for the
    derivation). Returns 0.0 for templates without priors or when the
    prior-consistent span needs no retraction."""
    if side_length < 1e-8 or not template_id:
        return 0.0
    priors = bridge_geometry_priors(template_id)
    if priors is None or priors.carbon_angle_deg is None:
        return 0.0
    target = bridge_target_distance(template_id)
    if template_id == "azine_bridge":
        # Asymmetric case: the fit pins the nitrogens on the N-N axis, so
        # consistency is the endpoint-triangle closure at the carbon-angle
        # prior (module docstring).
        cosine = max(-1.0, min(1.0, cos(radians(priors.carbon_angle_deg))))
        closure = sqrt(side_length * side_length + target * target - 2.0 * side_length * target * cosine)
        naive = target + side_length
        fraction = (naive - closure) / (2.0 * side_length)
    else:
        if priors.nitrogen_angle_deg is None:
            return 0.0
        span = required_bridge_span(
            side_length,
            target,
            side_length,
            priors.carbon_angle_deg,
            priors.nitrogen_angle_deg,
        )
        if span is None:
            return 0.0
        naive = target + 2.0 * side_length
        fraction = (naive - span) / (2.0 * side_length)
    return min(max(fraction, 0.0), 0.9)


# Realized B-O bond length target for the five-membered boronate ester ring
# (B-O-C-C-O): 1.44 angstrom matches both the measured CoRE-COF baseline ring
# mean (1.43 angstrom) and the DREIDING B_2-O_3 equilibrium (0.79+0.66-0.01).
# Distinct from the boronate profile's bridge_target_distance, which is a
# boron-to-oxygen-centroid placement distance, not a bond length.
BORONATE_ESTER_BOND_TARGET_DISTANCE = 1.44

# Measured O-B-O angle inside the baseline boronate ester rings; the closure
# fit pulls the rigid catechol oxygens inward toward this value instead of
# leaving the ring at the free-catechol opening (~138 degrees).
BORONATE_ESTER_OBO_TARGET_ANGLE_DEG = 112.4

# Exact-closure tolerance for the realization-time boronate ester angle solve
# (derived/numerical): the golden-section solve in reaction_realization.py
# minimizes the squared O-B-O angle residual over the feasible boron-rotation
# interval and converges the residual to float64 noise (~1e-12 degrees)
# whenever exact closure at the prior is geometrically feasible. Residuals at
# or below this threshold are reported as exact closures; the value sits far
# above numerical noise and far below any chemically meaningful deviation.
BORONATE_OBO_EXACT_TOLERANCE_DEG = 0.01


def effective_motif_origin(
    template_id: str | None,
    monomer: MonomerSpec,
    motif,
) -> Vec3:
    origin = motif.frame.origin
    if template_id not in {"imine_bridge", "azine_bridge"}:
        return origin
    if not monomer.atom_positions:
        return origin
    if template_id == "imine_bridge" and motif.kind not in {"amine", "aldehyde"}:
        return origin
    if template_id == "azine_bridge" and motif.kind not in {"hydrazine", "aldehyde"}:
        return origin

    reactive_atom_id = motif.metadata.get("reactive_atom_id")
    anchor_atom_id = motif.metadata.get("anchor_atom_id")
    if not isinstance(reactive_atom_id, int) or not isinstance(anchor_atom_id, int):
        return origin
    if reactive_atom_id >= len(monomer.atom_positions) or anchor_atom_id >= len(monomer.atom_positions):
        return origin

    reactive_position = monomer.atom_positions[reactive_atom_id]
    anchor_position = monomer.atom_positions[anchor_atom_id]
    anchor_to_reactive = sub(reactive_position, anchor_position)
    side_length = norm(anchor_to_reactive)
    if side_length < 1e-8:
        return origin
    retraction_fraction = derived_origin_retraction_fraction(template_id, side_length)
    return sub(
        reactive_position,
        scale(anchor_to_reactive, retraction_fraction),
    )
