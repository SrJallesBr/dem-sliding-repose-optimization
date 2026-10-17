import math
import sys
import unittest
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT))
from optimizer.measurements import measure_v_sliding_angle


class ToeAnchoredSlidingMeasurementTests(unittest.TestCase):
    def test_excludes_flat_top_and_anchors_each_face_at_its_lowest_surface_point(self):
        expected = {"left": 30.0, "right": 45.0}
        points = []
        toe_z = {"left": -0.1, "right": 0.1}
        for side, angle in expected.items():
            sign = -1.0 if side == "left" else 1.0
            for step in range(1, 81):
                distance = step / 80.0 * 1.2
                height = math.tan(math.radians(angle)) * min(distance, 0.65)
                z = toe_z[side] + sign * distance
                for x in (-0.01, 0.0, 0.01):
                    points.append((x, height, z))
            for x in (-0.01, 0.0, 0.01):
                points.append((x, 0.0, toe_z[side]))

        # Extra upper particles represent a cap/lobe that does not belong to the
        # straight pile face and should not enter the fitted segment.
        for side, angle in expected.items():
            sign = -1.0 if side == "left" else 1.0
            for step in range(48, 81):
                distance = step / 80.0 * 1.2
                z = toe_z[side] + sign * distance
                points.append((0.0, math.tan(math.radians(angle)) * 0.65 + 0.25, z))

        measured = measure_v_sliding_angle(
            points,
            horizontal_axis="z",
            section_half_width=0.02,
            mesh_bounds={"x": (-0.5, 0.5), "z": (-1.4, 1.4)},
        )
        self.assertAlmostEqual(measured["left_angle_deg"], expected["left"], delta=1.0)
        self.assertAlmostEqual(measured["right_angle_deg"], expected["right"], delta=1.0)
        self.assertAlmostEqual(measured["toe_left"][0], toe_z["left"], delta=0.06)
        self.assertAlmostEqual(measured["toe_right"][0], toe_z["right"], delta=0.06)
        self.assertLess(measured["left_fit_end_distance"], 0.9)
        self.assertLess(measured["right_fit_end_distance"], 0.9)
        self.assertGreaterEqual(len(measured["left_fit_points"]), 12)
        self.assertGreaterEqual(len(measured["right_fit_points"]), 12)
        self.assertGreater(len(measured["left_excluded_points"]), 0)
        self.assertGreater(len(measured["right_excluded_points"]), 0)
        self.assertLess(max(abs(x - measured["toe_left"][0]) for x, _ in measured["left_fit_points"]), 0.9)
        self.assertLess(max(abs(x - measured["toe_right"][0]) for x, _ in measured["right_fit_points"]), 0.9)
        self.assertAlmostEqual(measured["mean_angle_deg"], 37.5, delta=1.0)


if __name__ == "__main__":
    unittest.main()
