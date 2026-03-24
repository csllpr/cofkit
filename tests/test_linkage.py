import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from cofkit.chem import orient_frames_for_bridge
from cofkit.geometry import Frame


class LinkageTests(unittest.TestCase):
    def test_orient_frames(self):
        f1 = Frame.xy()
        f2 = Frame.xy()
        _, t2 = orient_frames_for_bridge(f1, f2, target_distance=1.28)
        self.assertEqual(t2, (1.28, 0.0, 0.0))


if __name__ == "__main__":
    unittest.main()
