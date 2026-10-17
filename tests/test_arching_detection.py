import unittest

from optimizer.measurements import detect_vertical_profile_arching


class VerticalArchingDetectionTests(unittest.TestCase):
    @staticmethod
    def profile(distances, particle_radius=0.1):
        points = []
        radii = []
        band_height = 2.0 * particle_radius
        for index, distance in enumerate(distances):
            y = (index + 0.5) * band_height
            for x in (-distance, distance):
                points.append((x, y, 0.0))
                radii.append(particle_radius)
        return points, radii

    def test_particle_radii_contribute_to_arching_surface_profile(self):
        particle_radii = [0.05, 0.07, 0.09, 0.11, 0.13, 0.15, 0.17]
        band_height = 2.0 * 0.11
        points = []
        radii = []
        for index, radius in enumerate(particle_radii):
            y = (index + 0.5) * band_height
            for x in (-0.5, 0.5):
                points.append((x, y, 0.0))
                radii.append(radius)
        result = detect_vertical_profile_arching(
            points, radii, horizontal_axis="x", section_half_width=0.05,
            min_consecutive_segments=2,
        )
        self.assertTrue(result["detected"])
        self.assertAlmostEqual(result["sides"]["right"]["profile"][0]["outward_distance"], 0.55)
        self.assertAlmostEqual(result["sides"]["left"]["profile"][0]["outward_distance"], 0.55)

    def test_sustained_outward_turn_above_ninety_degrees_is_arching(self):
        # Both side profiles first narrow normally, then turn outward across
        # multiple successive diameter-height bands.
        points, radii = self.profile([1.0, 0.8, 0.6, 0.4, 0.6, 0.8])
        result = detect_vertical_profile_arching(
            points, radii, horizontal_axis="x", section_half_width=0.05,
            min_consecutive_segments=2,
        )
        self.assertTrue(result["detected"])
        for side in ("left", "right"):
            diagnosis = result["sides"][side]
            self.assertGreater(diagnosis["max_directed_angle_deg"], 90.0)
            self.assertGreaterEqual(diagnosis["consecutive_overhang_segments"], 2)

    def test_normal_inward_sloping_faces_and_single_noisy_outward_step_are_not_arching(self):
        points, radii = self.profile([1.0, 0.8, 0.6, 0.4, 0.6, 0.4])
        result = detect_vertical_profile_arching(
            points, radii, horizontal_axis="x", section_half_width=0.05,
            min_consecutive_segments=2,
        )
        self.assertFalse(result["detected"])
        self.assertTrue(all(side["consecutive_overhang_segments"] < 2 for side in result["sides"].values()))


if __name__ == "__main__":
    unittest.main()
