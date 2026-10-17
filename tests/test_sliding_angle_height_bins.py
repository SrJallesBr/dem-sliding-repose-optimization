import math
import sys
import unittest
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT))
from optimizer.measurements import measure_v_sliding_angle_by_height, parse_vtk_scalar_field


class HeightBandCenterlineTests(unittest.TestCase):
    def test_reads_per_particle_radii_from_vtk_field(self):
        vtk = "POINT_DATA 3\nFIELD FieldData 1\nradius 1 3 double\n0.125 0.1375 0.2625\n\n"
        self.assertEqual(parse_vtk_scalar_field(vtk, "radius"), [0.125, 0.1375, 0.2625])

    def test_reads_exact_declared_values_when_followed_by_another_field_without_blank_line(self):
        vtk = "POINT_DATA 3\nFIELD FieldData 2\nradius 1 3 double\n0.125 0.1375 0.2625\nomega 3 3 double\n0 0 0 1 1 1 2 2 2\n"
        self.assertEqual(parse_vtk_scalar_field(vtk, "radius"), [0.125, 0.1375, 0.2625])

    def test_individual_radii_shift_profile_to_particle_outer_surface(self):
        points = []
        radii = []
        radii_by_band = [0.05, 0.15] * 6
        center_x = 0.0
        for side_sign in (-1.0, 1.0):
            for i, radius in enumerate(radii_by_band):
                y = 0.1 + i * 0.2
                center_distance = 0.1 + i * 0.1
                points.append((center_x, y, side_sign * center_distance))
                radii.append(radius)
        measured = measure_v_sliding_angle_by_height(
            points, particle_radii=radii, horizontal_axis="z",
            mesh_bounds={"x": (-0.5, 0.5), "z": (-2.0, 2.0)},
            min_fit_points=5, min_r2=0.0,
        )
        # Fit coordinates must describe the particle's exposed outer edge,
        # not only the centers; alternating radii make the effect observable.
        self.assertAlmostEqual(measured["left_fit_points"][0][0], -0.15)
        self.assertAlmostEqual(measured["right_fit_points"][0][0], 0.15)
        self.assertAlmostEqual(measured["height_bin_size"], 0.2)
        self.assertAlmostEqual(measured["median_radius"], 0.1)

    def test_tracks_nearest_particle_per_band_and_keeps_top_jump_downweighted(self):
        points = []
        radii = []
        expected_angle = 45.0
        radius = 0.1
        base_y = 0.05
        base_z = {"left": -0.05, "right": 0.05}

        for side in ("left", "right"):
            sign = -1.0 if side == "left" else 1.0
            for i in range(31):
                y = base_y + i * radius + (0.01 if i else 0.0)
                z = base_z[side] + sign * i * radius
                for x in (-0.01, 0.0, 0.01):
                    points.append((x, y, z))
                    radii.append(radius)
                if i > 0:
                    # A farther particle at the same height must not replace
                    # the surface particle nearest the centerline Z=0.
                    points.append((0.0, y, sign * (abs(z) + 0.45)))
                    radii.append(radius)

            # One upper-band candidate jumps far beyond the face: reject it,
            # then continue tracking the next valid height band.
            y = base_y + 32 * radius + 0.03
            points.append((0.0, y, sign * 6.0))
            radii.append(radius)

        measured = measure_v_sliding_angle_by_height(
            points,
            particle_radii=radii,
            horizontal_axis="z",
            mesh_bounds={"x": (-0.5, 0.5), "z": (-7.0, 7.0)},
        )
        self.assertAlmostEqual(measured["left_angle_deg"], expected_angle, delta=1.0)
        self.assertAlmostEqual(measured["right_angle_deg"], expected_angle, delta=1.0)
        self.assertAlmostEqual(measured["mean_angle_deg"], expected_angle, delta=1.0)
        self.assertGreaterEqual(len(measured["left_fit_points"]), 15)
        self.assertGreaterEqual(len(measured["right_fit_points"]), 15)
        self.assertAlmostEqual(measured["base_left"][0], base_z["left"] - radius, delta=1e-6)
        self.assertAlmostEqual(measured["base_right"][0], base_z["right"] + radius, delta=1e-6)
        self.assertEqual(len(measured["left_discarded_points"]), 0)
        self.assertEqual(len(measured["right_discarded_points"]), 0)
        self.assertEqual(len(measured["left_downweighted_points"]), 1)
        self.assertEqual(len(measured["right_downweighted_points"]), 1)
        self.assertEqual(measured["height_bin_size"], 2 * radius)
        self.assertEqual(measured["jump_limit"], 5 * radius)

    def test_robust_fit_preserves_a_continuous_shallow_face(self):
        points = []
        radii = []
        radius = 0.05
        angle = 15.0
        base_y = 0.01
        base_outward = 0.15
        for index in range(12):
            y = base_y + 2 * radius * index
            outward = base_outward + (y - base_y) / math.tan(math.radians(angle))
            for sign in (-1.0, 1.0):
                points.append((0.0, y, sign * outward))
                radii.append(radius)

        measured = measure_v_sliding_angle_by_height(
            points, particle_radii=radii, horizontal_axis="z",
            mesh_bounds={"x": (-0.5, 0.5), "z": (-6.0, 6.0)},
        )
        self.assertAlmostEqual(measured["left_angle_deg"], angle, delta=1.0)
        self.assertAlmostEqual(measured["right_angle_deg"], angle, delta=1.0)

    def test_full_face_fit_includes_supported_upper_profile_points(self):
        points = []
        radii = []
        radius = 0.05
        expected_angle = 35.0
        base_y = 0.05
        base_distance = 0.12
        slope = math.tan(math.radians(expected_angle))
        last_distance = None
        for index in range(12):
            y = base_y + 2 * radius * index
            distance = base_distance + (y - base_y) / slope
            last_distance = distance
            for sign in (-1.0, 1.0):
                points.append((0.0, y, sign * distance))
                radii.append(radius)
        for index in range(1, 6):
            y = base_y + 2 * radius * (11 + index)
            for sign in (-1.0, 1.0):
                points.append((0.0, y, sign * last_distance))
                radii.append(radius)

        measured = measure_v_sliding_angle_by_height(
            points, particle_radii=radii, horizontal_axis="z",
            mesh_bounds={"x": (-0.5, 0.5), "z": (-3.0, 3.0)},
        )

        self.assertGreater(measured["left_angle_deg"], 0.0)
        self.assertGreater(measured["right_angle_deg"], 0.0)
        self.assertEqual(len(measured["left_fit_points"]), 17)
        self.assertEqual(len(measured["right_fit_points"]), 17)
        self.assertAlmostEqual(max(y for _, y in measured["left_fit_points"]), base_y + 16 * 2 * radius)
        self.assertAlmostEqual(max(y for _, y in measured["right_fit_points"]), base_y + 16 * 2 * radius)

    def test_moderate_rogue_remains_in_total_fit_with_reduced_robust_weight(self):
        points = []
        radii = []
        radius = 0.125
        rogue = None
        for i in range(9):
            y = radius + i * 2 * radius
            left_z = -(5.0 - i * 2 * radius)
            right_z = 5.0 - i * 2 * radius
            points.extend([(0.0, y, left_z), (0.0, y, right_z)])
            radii.extend([radius, radius])
            if i == 4:
                # This single particle is selected because it is nearer the
                # centerline, although it is off the otherwise straight face.
                rogue = (right_z - 0.6, y)
                points.append((0.0, y, rogue[0]))
                radii.append(radius)

        measured = measure_v_sliding_angle_by_height(
            points, particle_radii=radii, horizontal_axis="z",
            mesh_bounds={"x": (-0.5, 0.5), "z": (-6.0, 6.0)},
        )
        self.assertAlmostEqual(measured["left_angle_deg"], 45.0, delta=0.5)
        self.assertAlmostEqual(measured["right_angle_deg"], 45.0, delta=2.0)
        self.assertIn((rogue[0] + radius, rogue[1]), measured["right_fit_points"])
        self.assertIn((rogue[0] + radius, rogue[1]), measured["right_downweighted_points"])


    def test_low_isolated_particle_does_not_control_face_base(self):
        points = []
        radii = []
        radius = 0.125
        base_y = 0.125
        for i in range(12):
            y = base_y + i * 2 * radius
            outward = 0.5 + (y - base_y)
            points.extend([(0.0, y, -outward), (0.0, y, outward)])
            radii.extend([radius, radius])
        # One low, near-center particle is not part of either exposed face.
        points.append((0.0, -0.1, 0.001))
        radii.append(radius)
        measured = measure_v_sliding_angle_by_height(
            points, particle_radii=radii, horizontal_axis="z",
            mesh_bounds={"x": (-0.5, 0.5), "z": (-6.0, 6.0)},
        )
        self.assertAlmostEqual(measured["base_left"][0], -0.625, delta=1e-6)
        self.assertAlmostEqual(measured["base_left"][1], base_y, delta=1e-6)
        self.assertAlmostEqual(measured["base_right"][0], 0.625, delta=1e-6)
        self.assertAlmostEqual(measured["base_right"][1], base_y, delta=1e-6)
        self.assertAlmostEqual(measured["left_angle_deg"], 45.0, delta=0.5)
        self.assertAlmostEqual(measured["right_angle_deg"], 45.0, delta=0.5)
        self.assertAlmostEqual(measured["mean_angle_deg"], 45.0, delta=0.5)
    def test_single_near_centerline_particle_with_moderate_profile_error_is_kept(self):
        points = []
        radii = []
        radius = 0.1
        base_y = 0.1
        expected_angle = 45.0
        rogue = None
        for index in range(12):
            y = base_y + index * 2 * radius
            surface_distance = 0.5 + index * 2 * radius
            for sign in (-1.0, 1.0):
                center_distance = surface_distance - radius
                points.append((0.0, y, sign * center_distance))
                radii.append(radius)
            if index == 5:
                # A moderate inward-shifted sphere is nearer the centerline,
                # so the legacy per-band selector chooses it over the true face.
                rogue = (surface_distance - radius - 0.25, y)
                points.append((0.0, rogue[1], rogue[0]))
                radii.append(radius)

        measured = measure_v_sliding_angle_by_height(
            points, particle_radii=radii, horizontal_axis="z",
            mesh_bounds={"x": (-0.5, 0.5), "z": (-6.0, 6.0)},
        )
        self.assertAlmostEqual(measured["right_angle_deg"], expected_angle, delta=5.0)
        self.assertEqual(len(measured["right_fit_points"]), 12)
        self.assertIn((rogue[0] + radius, rogue[1]), measured["right_fit_points"])


    def test_full_profile_fit_retains_every_supported_band_from_base_to_top(self):
        points = []
        radii = []
        radius = 0.1
        angle = 45.0
        base_y = 0.1
        for index in range(12):
            y = base_y + 2 * radius * index
            outward = 0.5 + 2 * radius * index
            left_outward = outward
            right_outward = outward
            if index == 6:
                # A modest profile deviation should not make an otherwise
                # supported face lose that band from the total base-to-top fit.
                right_outward -= 0.35
            points.extend([(0.0, y, -(left_outward - radius)),
                           (0.0, y, right_outward - radius)])
            radii.extend([radius, radius])
        measured = measure_v_sliding_angle_by_height(
            points, particle_radii=radii, horizontal_axis="z",
            mesh_bounds={"x": (-0.5, 0.5), "z": (-6.0, 6.0)},
        )
        self.assertEqual(len(measured["left_fit_points"]), 12)
        self.assertEqual(len(measured["right_fit_points"]), 12)
        self.assertAlmostEqual(measured["left_fit_points"][0][1], base_y)
        self.assertAlmostEqual(measured["right_fit_points"][0][1], base_y)
        self.assertAlmostEqual(max(y for _, y in measured["left_fit_points"]), base_y + 2 * radius * 11)
        self.assertAlmostEqual(max(y for _, y in measured["right_fit_points"]), base_y + 2 * radius * 11)


if __name__ == "__main__":
    unittest.main()
