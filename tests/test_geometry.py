import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from cofkit.geometry import Frame, frame_axes, matmul_vec, norm, rotation_from_frame_to_axes, sub


class GeometryTests(unittest.TestCase):
    def test_rotation_from_frame_to_axes_maps_orthonormalized_axes(self):
        source_frame = Frame(origin=(0.0, 0.0, 0.0), primary=(1.0, 2.0, 0.1), normal=(0.2, -0.1, 1.0))
        target_frame = Frame(origin=(0.0, 0.0, 0.0), primary=(0.0, 1.0, 0.2), normal=(0.0, 0.0, 1.0))

        rotation = rotation_from_frame_to_axes(
            source_frame,
            target_frame.primary,
            target_frame.normal,
        )
        source_axes = frame_axes(source_frame)
        target_axes = frame_axes(target_frame)

        self.assertLess(norm(sub(matmul_vec(rotation, source_axes[0]), target_axes[0])), 1e-8)
        self.assertLess(norm(sub(matmul_vec(rotation, source_axes[1]), target_axes[1])), 1e-8)
        self.assertLess(norm(sub(matmul_vec(rotation, source_axes[2]), target_axes[2])), 1e-8)


if __name__ == "__main__":
    unittest.main()
