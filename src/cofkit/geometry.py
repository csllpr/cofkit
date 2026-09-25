from __future__ import annotations

import sys
from dataclasses import dataclass
from math import acos, cos, degrees, radians, sin, sqrt
from typing import Iterable, Mapping, Sequence

Vec3 = tuple[float, float, float]
Mat3 = tuple[Vec3, Vec3, Vec3]


def vec3(x: float, y: float, z: float) -> Vec3:
    return (float(x), float(y), float(z))


def add(a: Vec3, b: Vec3) -> Vec3:
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def sub(a: Vec3, b: Vec3) -> Vec3:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def scale(v: Vec3, factor: float) -> Vec3:
    return (v[0] * factor, v[1] * factor, v[2] * factor)


def dot(a: Vec3, b: Vec3) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def cross(a: Vec3, b: Vec3) -> Vec3:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def norm(v: Vec3) -> float:
    return sqrt(dot(v, v))


def normalize(v: Vec3) -> Vec3:
    length = norm(v)
    if length == 0.0:
        raise ValueError("cannot normalize a zero-length vector")
    return scale(v, 1.0 / length)


def centroid(points: Iterable[Vec3]) -> Vec3:
    pts = list(points)
    if not pts:
        raise ValueError("cannot compute centroid of an empty point set")
    inv = 1.0 / len(pts)
    return (
        sum(p[0] for p in pts) * inv,
        sum(p[1] for p in pts) * inv,
        sum(p[2] for p in pts) * inv,
    )


def distance(a: Vec3, b: Vec3) -> float:
    return norm(sub(a, b))


def safe_normalize(
    vector: Vec3,
    *,
    fallback: Vec3 = (0.0, 0.0, 1.0),
    warn_context: str | None = None,
) -> Vec3:
    """Normalize *vector*; return *fallback* verbatim when it is degenerate.

    Fallback policy (explicit per call site): a (near-)zero vector has no
    direction, so the caller-supplied *fallback* is returned unchanged.  It
    is a fabricated direction, not something derived from the input — callers
    whose fallback feeds a physically meaningful axis should treat a triggered
    fallback as a data problem.  When *warn_context* is given, a ``warning:``
    line naming the context is printed to stderr; otherwise the fallback is
    silent by the caller's choice (e.g. hot loops that handle degeneracy
    upstream).
    """
    length = norm(vector)
    if length < 1e-8:
        if warn_context is not None:
            print(
                f"warning: {warn_context}; using fallback direction {fallback}",
                file=sys.stderr,
            )
        return fallback
    return scale(vector, 1.0 / length)


def orthogonal_component(vector: Vec3, axis: Vec3, *, axis_is_unit: bool = False) -> Vec3:
    """Return the component of *vector* orthogonal to *axis*.

    Degenerate-axis policy: a (near-)zero axis defines no projection
    direction, so *vector* is returned unchanged.  Pass
    ``axis_is_unit=True`` when the caller guarantees ``norm(axis) == 1`` to
    skip re-normalization (the projection formula is only correct for a unit
    axis).
    """
    if not axis_is_unit:
        length = norm(axis)
        if length < 1e-8:
            return vector
        axis = scale(axis, 1.0 / length)
    return sub(vector, scale(axis, dot(vector, axis)))


def angle_degrees(
    point_a: Sequence[float],
    point_b: Sequence[float],
    point_c: Sequence[float],
    *,
    on_degenerate: float | None,
) -> float | None:
    """Angle ∠abc (vertex at *point_b*) in degrees; works for 2D or 3D points.

    Degenerate policy: when either arm has (near-)zero length the angle is
    undefined, and the caller-chosen *on_degenerate* value is returned
    (``180.0`` where a straight angle is the safe neutral default, ``None``
    where the caller must handle the failure explicitly).
    """
    ba = tuple(x - y for x, y in zip(point_a, point_b))
    bc = tuple(x - y for x, y in zip(point_c, point_b))
    norm_ba = sqrt(sum(x * x for x in ba))
    norm_bc = sqrt(sum(x * x for x in bc))
    if norm_ba < 1e-8 or norm_bc < 1e-8:
        return on_degenerate
    cosine = sum(x * y for x, y in zip(ba, bc)) / (norm_ba * norm_bc)
    return degrees(acos(max(-1.0, min(1.0, cosine))))


def mat3_identity() -> Mat3:
    return (
        (1.0, 0.0, 0.0),
        (0.0, 1.0, 0.0),
        (0.0, 0.0, 1.0),
    )


def transpose(m: Mat3) -> Mat3:
    return (
        (m[0][0], m[1][0], m[2][0]),
        (m[0][1], m[1][1], m[2][1]),
        (m[0][2], m[1][2], m[2][2]),
    )


def matmul(a: Mat3, b: Mat3) -> Mat3:
    b_t = transpose(b)
    return tuple(
        tuple(dot(row, col) for col in b_t)
        for row in a
    )  # type: ignore[return-value]


def matmul_vec(m: Mat3, v: Vec3) -> Vec3:
    return (
        dot(m[0], v),
        dot(m[1], v),
        dot(m[2], v),
    )


def frame_axes(frame: "Frame") -> Mat3:
    primary = normalize(frame.primary)
    normal_seed = sub(frame.normal, scale(primary, dot(frame.normal, primary)))
    if norm(normal_seed) < 1e-8:
        fallback = (0.0, 0.0, 1.0) if abs(primary[2]) < 0.9 else (1.0, 0.0, 0.0)
        normal_seed = sub(fallback, scale(primary, dot(fallback, primary)))
    normal = normalize(normal_seed)
    secondary = normalize(cross(normal, primary))
    normal = normalize(cross(primary, secondary))
    return (primary, secondary, normal)


def rotation_from_frame_to_axes(frame: "Frame", target_primary: Vec3, target_normal: Vec3) -> Mat3:
    source = frame_axes(frame)
    target = frame_axes(Frame(origin=(0.0, 0.0, 0.0), primary=target_primary, normal=target_normal))
    return matmul(transpose(target), source)


@dataclass(frozen=True)
class Frame:
    origin: Vec3
    primary: Vec3
    normal: Vec3

    def normalized(self) -> "Frame":
        return Frame(
            origin=self.origin,
            primary=normalize(self.primary),
            normal=normalize(self.normal),
        )

    @staticmethod
    def xy(origin: Vec3 = (0.0, 0.0, 0.0)) -> "Frame":
        return Frame(origin=origin, primary=(1.0, 0.0, 0.0), normal=(0.0, 0.0, 1.0))

    @staticmethod
    def yz(origin: Vec3 = (0.0, 0.0, 0.0)) -> "Frame":
        return Frame(origin=origin, primary=(0.0, 1.0, 0.0), normal=(1.0, 0.0, 0.0))

    @staticmethod
    def zx(origin: Vec3 = (0.0, 0.0, 0.0)) -> "Frame":
        return Frame(origin=origin, primary=(0.0, 0.0, 1.0), normal=(0.0, 1.0, 0.0))


# ---------------------------------------------------------------------------
# Stacking support: layer-span measurement and cell classification
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class LayerSpanReport:
    """Result of measuring the z-extent of a layer's atoms.

    ``axis`` is a *label* honestly naming the measurement axis actually used
    (e.g. ``"layer_normal"``, ``"c_hat"``, or ``"unknown"`` for values taken
    from upstream metadata), never a hard-coded claim.
    """

    span: float
    mode: str           # "atomistic_product" | "precursor_coordinates" | "embedding_metadata" | "unavailable"
    axis: str           # label of the axis actually used: "layer_normal" | "c_hat" | "unknown"
    n_atoms: int
    z_translation_included: bool


def layer_normal_axis(
    a_vec: Vec3,
    b_vec: Vec3,
    *,
    fallback_axis: Vec3 = (0.0, 0.0, 1.0),
) -> tuple[Vec3, str]:
    """Return ``(unit_axis, label)`` for layer z-span measurement.

    The preferred axis is the layer normal ``n = normalize(a × b)`` (label
    ``"layer_normal"``): a span measured along ``c_hat`` is not invariant
    under in-plane periodic images when c is tilted away from the ab normal
    (real fitted cells tilt ~0.5°), whereas the layer normal is.

    Degenerate policy: when ``a × b`` is (near-)zero (collinear or missing
    in-plane vectors) a ``warning:`` is printed and the normalized
    *fallback_axis* is used with the honest label ``"c_hat"``; if that is
    degenerate too, ``(0, 0, 1)`` is returned with the same label.
    """
    normal = cross(a_vec, b_vec)
    if norm(normal) >= 1e-8:
        return normalize(normal), "layer_normal"
    print(
        "warning: geometry: degenerate in-plane cell vectors (a × b ≈ 0); "
        "falling back to the fallback axis for layer-span measurement",
        file=sys.stderr,
    )
    if norm(fallback_axis) >= 1e-8:
        return normalize(fallback_axis), "c_hat"
    return (0.0, 0.0, 1.0), "c_hat"


def measure_layer_z_span(
    poses: "Mapping[str, object]",
    monomer_specs: "Mapping[str, object]",
    instance_to_monomer: "Mapping[str, str]",
    *,
    axis: Vec3 = (0.0, 0.0, 1.0),
    axis_label: str = "c_hat",
    realization: object = None,
) -> LayerSpanReport:
    """Measure the layer z-span (nuclear-plane extent along *axis*).

    The axis choice is explicit: pass ``axis``/``axis_label`` from
    :func:`layer_normal_axis` when a cell is available (preferred — invariant
    under in-plane periodic images in tilted cells), or an explicit c-axis
    unit vector with ``axis_label="c_hat"`` when only a c axis exists.

    Works in one of two modes depending on what is available:

    * ``atomistic_product``: realized atom positions from ``realization``
      (preferred — matches the ring-forming path).
    * ``precursor_coordinates``: monomer ``atom_positions`` rotated and
      translated via the pose.

    The pose translation **is** included: different instances of a layer can
    sit at different heights along the axis, and omitting the translation
    collapses their contributions (the documented root cause of shipped
    interpenetrated stacks).

    Parameters
    ----------
    poses:
        Mapping ``instance_id → Pose``-like object with ``.translation`` and
        ``.rotation_matrix``.
    monomer_specs:
        Mapping ``monomer_id → MonomerSpec``-like with ``.atom_positions``.
    instance_to_monomer:
        Maps instance ids to monomer ids.
    axis:
        Unit vector along the measurement axis (already normalised by
        caller; see :func:`layer_normal_axis`).
    axis_label:
        Honest label for *axis*, reported verbatim in the report.
    realization:
        Optional ``ReactionRealizationResult``-like with
        ``.atoms_by_instance`` mapping ``instance_id → sequence of atoms``
        where each atom has ``.local_position``.
    """
    z_values: list[float] = []
    used_atomistic = False

    for instance_id, pose in poses.items():
        monomer_id = instance_to_monomer.get(instance_id)
        if monomer_id is None:
            continue
        monomer = monomer_specs.get(monomer_id)
        if monomer is None:
            continue

        # Prefer realized product atoms
        realized_atoms = None
        if realization is not None:
            atoms_by_instance = getattr(realization, "atoms_by_instance", {})
            realized_atoms = atoms_by_instance.get(instance_id)

        translation = getattr(pose, "translation", (0.0, 0.0, 0.0))
        rotation = getattr(pose, "rotation_matrix", ((1, 0, 0), (0, 1, 0), (0, 0, 1)))

        if realized_atoms is not None:
            positions = [atom.local_position for atom in realized_atoms]
            used_atomistic = True
        else:
            positions = list(getattr(monomer, "atom_positions", ()))

        for local_pos in positions:
            world = add(matmul_vec(rotation, local_pos), translation)
            z_values.append(dot(world, axis))

    if not z_values:
        return LayerSpanReport(
            span=0.0,
            mode="unavailable",
            axis=axis_label,
            n_atoms=0,
            z_translation_included=True,
        )

    span = max(z_values) - min(z_values)
    mode = "atomistic_product" if used_atomistic else "precursor_coordinates"
    return LayerSpanReport(
        span=span,
        mode=mode,
        axis=axis_label,
        n_atoms=len(z_values),
        z_translation_included=True,
    )


def classify_2d_cell(
    cell: "tuple[Vec3, Vec3, Vec3]",
    *,
    angle_tolerance_deg: float = 0.5,
    length_rtol: float = 0.01,
) -> "tuple[str, str]":
    """Classify a 2D cell into ``(kind, setting)``.

    Returns
    -------
    kind:
        ``"hexagonal"``, ``"square"``, ``"orthogonal"``, or ``"oblique"``.
    setting:
        For hexagonal cells: ``"60deg"`` (γ ≈ 60°) or ``"120deg"`` (γ ≈ 120°).
        Empty string for all other kinds.

    The looser ±0.5° angle tolerance (vs. the old 1e-3 cosine tolerance) means
    a fitted cell with γ = 60.09° is correctly classified as hexagonal rather
    than falling through to ``"oblique"``.
    """
    a_vec, b_vec, _ = cell
    a_len = norm(a_vec)
    b_len = norm(b_vec)
    if a_len < 1e-8 or b_len < 1e-8:
        return "oblique", ""

    cosine = dot(a_vec, b_vec) / (a_len * b_len)
    cosine = max(-1.0, min(1.0, cosine))
    gamma_deg = degrees(acos(cosine))
    lengths_equal = abs(a_len - b_len) / max(a_len, b_len) < length_rtol

    if lengths_equal and abs(gamma_deg - 60.0) <= angle_tolerance_deg:
        return "hexagonal", "60deg"
    if lengths_equal and abs(gamma_deg - 120.0) <= angle_tolerance_deg:
        return "hexagonal", "120deg"
    if lengths_equal and abs(gamma_deg - 90.0) <= angle_tolerance_deg:
        return "square", ""
    if abs(gamma_deg - 90.0) <= angle_tolerance_deg:
        return "orthogonal", ""
    return "oblique", ""


def classify_2d_cell_parameters(
    cell_parameters: "tuple[float, ...]",
    *,
    angle_tolerance_deg: float = 0.5,
    length_rtol: float = 0.01,
) -> "tuple[str, str]":
    """Classify a 2D cell from RCSR-style cell *parameters*.

    Accepts ``(a, b, c, alpha, beta, gamma)`` or ``(a, b, gamma)``; only the
    in-plane metric (a, b, γ) is used.  This is an adapter over
    :func:`classify_2d_cell` for topology *definition* cells, which are stored
    as parameters in the γ = 120° convention — the same tolerances and the
    same ``(kind, setting)`` contract apply, so a definition cell and a built
    cell classify identically.  Returns ``("unknown", "")`` when the
    parameter count is not recognized.
    """
    if len(cell_parameters) >= 6:
        a, b, _c, _alpha, _beta, gamma = cell_parameters[:6]
    elif len(cell_parameters) == 3:
        a, b, gamma = cell_parameters
    else:
        return "unknown", ""
    gamma_radians = radians(gamma)
    cell = (
        (a, 0.0, 0.0),
        (b * cos(gamma_radians), b * sin(gamma_radians), 0.0),
        (0.0, 0.0, 1.0),
    )
    return classify_2d_cell(
        cell,
        angle_tolerance_deg=angle_tolerance_deg,
        length_rtol=length_rtol,
    )
