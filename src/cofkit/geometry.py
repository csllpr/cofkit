from __future__ import annotations

from dataclasses import dataclass
from math import sqrt
from typing import Iterable

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
