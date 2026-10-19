import math
import sys
import tempfile
import unittest
from unittest.mock import patch
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT))

from optimizer.calibration_controller import CalibrationController, GeometryDiagnostic, measure_pile_angle
from optimizer.measurements import measure_repose_eight_ray_least_squares, measure_v_sliding_angle, parse_ascii_stl_bounds


class ArticleReposeMeasurementTests(unittest.TestCase):
    def test_sliding_measurement_recovers_both_pile_angles_from_horizontal_x_axis(self):
        expected_left, expected_right = 30.0, 45.0
        points = []
        for step in range(1, 61):
            for x_sign, angle in ((-1.0, expected_left), (1.0, expected_right)):
                x = x_sign * step / 60.0
                y = math.tan(math.radians(angle)) * abs(x)
                for z in (-0.01, 0.0, 0.01):
                    points.append((x, y, z))
        result = measure_v_sliding_angle(points, horizontal_axis="x")
        self.assertAlmostEqual(result["left_angle_deg"], expected_left, places=6)
        self.assertAlmostEqual(result["right_angle_deg"], expected_right, places=6)
        self.assertAlmostEqual(result["mean_angle_deg"], 37.5, places=6)
        self.assertAlmostEqual(result["asymmetry_deg"], 15.0, places=6)

    def test_sliding_measurement_recovers_both_pile_angles_when_z_is_horizontal(self):
        expected_left, expected_right = 30.0, 45.0
        points = []
        for step in range(1, 61):
            for sign, angle in ((-1.0, expected_left), (1.0, expected_right)):
                z = sign * step / 60.0
                y = math.tan(math.radians(angle)) * abs(z)
                for x in (-0.01, 0.0, 0.01):
                    points.append((x, y, z))
        result = measure_v_sliding_angle(points, horizontal_axis="z", section_half_width=0.01)
        self.assertAlmostEqual(result["left_angle_deg"], expected_left, places=6)
        self.assertAlmostEqual(result["right_angle_deg"], expected_right, places=6)

    def test_sliding_measurement_uses_mesh_center_for_zy_section_and_face_split(self):
        expected_left, expected_right = 25.0, 40.0
        points = []
        for step in range(1, 61):
            for sign, angle in ((-1.0, expected_left), (1.0, expected_right)):
                z = 7.0 + sign * step / 60.0
                y = 1.0 + math.tan(math.radians(angle)) * abs(z - 7.0)
                for x in (4.99, 5.0, 5.01):
                    points.append((x, y, z))
        result = measure_v_sliding_angle(
            points, horizontal_axis="z", section_half_width=0.02,
            mesh_bounds={"x": (0.0, 10.0), "z": (0.75, 13.25)},
        )
        self.assertEqual(result["section_center"], 5.0)
        self.assertEqual(result["horizontal_center"], 7.0)
        self.assertAlmostEqual(result["left_angle_deg"], expected_left, places=6)
        self.assertAlmostEqual(result["right_angle_deg"], expected_right, places=6)

    def test_ascii_mesh_bounds_include_nbr_scaled_geometry(self):
        stl = """solid ASCII
 facet normal 0 0 1
 outer loop
  vertex -6.25 1.0 -6.25
  vertex 6.25 1.0 6.25
  vertex 5.0 4.0 -6.25
 endloop
 endfacet
endsolid
"""
        self.assertEqual(parse_ascii_stl_bounds(stl), {"x": (-6.25, 6.25), "y": (1.0, 4.0), "z": (-6.25, 6.25)})

    def test_controller_uses_latest_scaled_box_mesh_for_measurement_bounds(self):
        with tempfile.TemporaryDirectory() as temporary:
            post = Path(temporary)
            for step, extent in ((90000, 5.0), (100000, 6.25)):
                stl = f"solid ASCII\n vertex {-extent} 1 {-extent}\n vertex {extent} 2 {extent}\nendsolid\n"
                (post / f"box_{step}.stl").write_text(stl, encoding="utf-8")
            self.assertEqual(
                CalibrationController._mesh_bounds_for_snapshot(post / "deslizamento_final.vtk"),
                {"x": (-6.25, 6.25), "z": (-6.25, 6.25)},
            )

    def test_actual_ten_second_bridge_shape_is_diagnostic_not_a_score(self):
        snapshot = PROJECT / "smoke_angle_horizontal_20261013_173200" / "trial_0001" / "Ang Deslizamento" / "post" / "deslizamento_final.vtk"
        if not snapshot.is_file():
            self.skipTest("real 10 s smoke snapshot is not present")
        with tempfile.TemporaryDirectory() as temporary:
            controller = CalibrationController(PROJECT, run_root=Path(temporary), solver_command=["true"])
            with self.assertRaisesRegex(GeometryDiagnostic, "bridge/arching suspected"):
                controller._measure_snapshot(snapshot, "discharge")

    def test_sequential_opening_snapshots_preserve_user_phase_order(self):
        from optimizer.opening_protocol import inspect_case,CASE_DIRECTORY,INPUT_NAME
        deck=(PROJECT/CASE_DIRECTORY/INPUT_NAME).read_text(encoding='utf-8')
        repose_dump=deck.index('post/repouso_final.vtk')
        self.assertLess(repose_dump,deck.index('fix move_base_esq'))
        self.assertGreater(deck.index('post/descarga_final.vtk'),deck.index('run ${steps_pos_descarga}'))
        self.assertEqual(inspect_case(PROJECT/CASE_DIRECTORY)['snapshot_times_s'],{'repose':13.0,'discharge':23.0})

    def test_repose_surface_profiles_use_each_particle_radius(self):
        angle_deg = 30.0
        slope = math.tan(math.radians(angle_deg))
        points = []
        particle_radii = []
        for radial_index in range(1, 15):
            radial = radial_index / 15.0
            height = 1.0 - slope * radial
            for direction_index in range(72):
                theta = 2.0 * math.pi * direction_index / 72.0
                particle_radius = 0.01 + 0.002 * (radial_index % 5)
                points.append((radial * math.cos(theta), height, radial * math.sin(theta)))
                particle_radii.append(particle_radius)
        result = measure_repose_eight_ray_least_squares(
            points, particle_radii=particle_radii, radial_bins=12, min_points_per_ray=8,
        )
        self.assertGreater(result["median_particle_radius"], 0)
        self.assertAlmostEqual(result["median_particle_radius"], 0.014)
        for profile in result["direction_profiles"]:
            for evidence in profile["sampling_evidence"]:
                h, y = evidence["surface"]
                center_h, center_y = evidence["center"]
                radius = evidence["section_radius"]
                self.assertAlmostEqual((h-center_h)**2+(y-center_y)**2, radius**2, places=10)
                self.assertGreaterEqual(center_y, profile["profile_start_height"])
                self.assertAlmostEqual(evidence["distance_to_axis"], min(item["distance_to_axis"] for item in evidence["eligible_candidates"]))

    def test_repose_snapshot_passes_vtk_particle_radii_into_measurement(self):
        points = []
        particle_radii = []
        for radial_index in range(1, 15):
            radial = radial_index / 15.0
            height = 1.0 - math.tan(math.radians(30.0)) * radial
            for direction_index in range(72):
                theta = 2.0 * math.pi * direction_index / 72.0
                points.append((radial * math.cos(theta), height, radial * math.sin(theta)))
                particle_radii.append(0.01 + 0.002 * (radial_index % 5))
        with tempfile.TemporaryDirectory() as temporary:
            snapshot = Path(temporary) / "trial_0001" / "Ang Repouso" / "post" / "repouso_10s.vtk"
            snapshot.parent.mkdir(parents=True)
            snapshot.write_text("synthetic VTK fixture", encoding="utf-8")
            controller = CalibrationController(PROJECT, run_root=Path(temporary) / "runs", solver_command=["true"])
            with patch("optimizer.calibration_controller.parse_vtk_points", return_value=points), patch("optimizer.calibration_controller.parse_vtk_scalar_field", return_value=particle_radii) as parse_radii:
                result = controller._measure_snapshot(snapshot, "repose")
        parse_radii.assert_called_once_with("synthetic VTK fixture", "radius")
        self.assertAlmostEqual(result["median_particle_radius"], 0.014)

    def test_eight_radial_profiles_recover_cone_angle_with_about_1000_particles(self):
        angle_deg = 30.0
        slope = math.tan(math.radians(angle_deg))
        points = []
        radii = []
        for radial_index in range(1, 15):
            radius = radial_index / 15.0
            height = 1.0 - slope * radius
            for direction_index in range(72):
                theta = 2.0 * math.pi * direction_index / 72.0
                for band in range(1, 2):
                    points.append((radius * math.cos(theta), height - 0.002 * band, radius * math.sin(theta)))
                    radii.append(0.014)
        result = measure_repose_eight_ray_least_squares(points, radial_bins=12, min_points_per_ray=8, particle_radii=radii)
        self.assertEqual(len(result["direction_angles_deg"]), 8)
        self.assertAlmostEqual(result["mean_angle_deg"], angle_deg, delta=0.2)
        self.assertLess(result["stddev_deg"], 1e-9)
        controller_result = measure_pile_angle(points, particle_radii=radii)
        self.assertEqual(len(controller_result["direction_angles_deg"]), 8)
        self.assertAlmostEqual(controller_result["mean_angle_deg"], angle_deg, delta=0.2)

    def test_isolated_particle_outside_heap_does_not_expand_radial_fit_domain(self):
        angle_deg = 30.0
        slope = math.tan(math.radians(angle_deg))
        points = []
        radii = []
        for radial_index in range(1, 15):
            radius = radial_index / 15.0
            height = 1.0 - slope * radius
            for direction_index in range(72):
                theta = 2.0 * math.pi * direction_index / 72.0
                points.append((radius * math.cos(theta), height, radius * math.sin(theta)))
                radii.append(0.014)
        points.append((4.0, 0.0, 0.0))
        radii.append(0.014)
        result = measure_repose_eight_ray_least_squares(points, radial_bins=12, min_points_per_ray=8, particle_radii=radii)
        self.assertAlmostEqual(result["mean_angle_deg"], angle_deg, delta=0.2)
        self.assertEqual(result["discarded_outliers"], 1)


if __name__ == "__main__":
    unittest.main()
