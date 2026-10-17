import math
import re
import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path
from unittest.mock import patch

from optimizer.calibration_controller import CalibrationController, MeasurementRejected, discharge_cut_result_fields
from optimizer.measurements import UnmeasurableAngle


class MeasurementCaptureAndDumpCleanupTests(unittest.TestCase):
    def _snapshot(self, root: Path, test_name: str, filename: str) -> Path:
        trial = root / "trial_0001"
        test_dir = trial / test_name
        post = test_dir / "post"
        post.mkdir(parents=True)
        ratio = 0.125 / 0.00325
        points = []
        for offset in (0.0, 0.1, -0.1):
            center = offset * ratio
            for index in range(3):
                y = 0.2 + 0.2 * index
                outer = 0.2 + 0.1 * index
                points.extend(((center, y, -outer), (center, y, outer)))
        coordinates = " ".join(str(value) for point in points for value in point)
        radii = " ".join("0.01" for _ in points)
        path = post / filename
        vtk_lines = [
            "# vtk DataFile Version 2.0", "fixture", "ASCII", "DATASET POLYDATA",
            f"POINTS {len(points)} float", coordinates,
            f"POINT_DATA {len(points)}", "FIELD FieldData 1",
            f"radius 1 {len(points)} float", radii,
        ]
        path.write_text(chr(10).join(vtk_lines) + chr(10), encoding="utf-8")
        deck_name = "in.repouso" if test_name == "Ang Repouso" else "in.deslizamento"
        deck_lines = [
            "variable radius1 equal 0.250/2",
            "variable radiusnbr equal 0.0065/2",
            "variable nbrratio equal v_radius1/v_radiusnbr",
        ]
        (test_dir / deck_name).write_text(chr(10).join(deck_lines) + chr(10), encoding="utf-8")
        return path

    def test_sliding_deck_has_only_the_measurement_particle_snapshot(self):
        deck = (Path(__file__).parents[1] / "Ang Deslizamento" / "in.deslizamento").read_text(encoding="utf-8")
        self.assertEqual(re.findall(r"(?im)^\s*dump\s+", deck), [])
        self.assertEqual(len(re.findall(r"(?im)^\s*write_dump\s+all\s+custom/vtk\s+post/deslizamento_final\.vtk\b", deck)), 1)
        self.assertNotIn("box_*.stl", deck)
        self.assertNotIn("walls_*.stl", deck)

    def test_mesh_bounds_are_derived_from_source_stl_and_trial_scale(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            snapshot = self._snapshot(root, "Ang Deslizamento", "deslizamento_final.vtk")
            test_dir = snapshot.parents[1]
            mesh_dir = test_dir / "meshes"
            mesh_dir.mkdir()
            (mesh_dir / "box.stl").write_text(
                "solid box\n facet normal 0 0 0\n outer loop\n"
                "vertex -0.125 1 -0.25\nvertex 0.125 1.5 0.25\nvertex 0 1.2 0\n"
                "endloop\nendfacet\nendsolid box\n",
                encoding="utf-8",
            )
            (test_dir / "in.deslizamento").write_text(
                "variable radiusnbr equal 0.0065/2\n"
                "variable radius1 equal 0.250/2\n"
                "variable nbrratio equal v_radius1/v_radiusnbr\n"
                "fix box all mesh/surface file meshes/box.stl type 2 scale $(v_nbrratio)\n",
                encoding="utf-8",
            )
            bounds = CalibrationController._mesh_bounds_for_snapshot(snapshot)
        self.assertAlmostEqual(bounds["x"][0], -4.8076923077, places=7)
        self.assertAlmostEqual(bounds["x"][1], 4.8076923077, places=7)
        self.assertAlmostEqual(bounds["z"][0], -9.6153846154, places=7)
        self.assertAlmostEqual(bounds["z"][1], 9.6153846154, places=7)

    def test_valid_discharge_and_eight_ray_repose_write_trial_svgs(self):
        ratio = 0.125 / 0.00325
        discharge_result = {
            "horizontal_axis": "z", "centerline": 0.0,
            "left_angle_deg": 30.0, "right_angle_deg": 45.0, "mean_angle_deg": 37.5,
            "height_bin_size": 0.02, "median_radius": 0.01,
            "left_r2": 1.0, "right_r2": 1.0,
            "left_candidate_count": 2, "right_candidate_count": 2,
            "left_downweighted_points": [], "right_downweighted_points": [],
            "left_fit_points": [(-0.2, 0.2), (-0.1, 0.3)],
            "right_fit_points": [(0.1, 0.2), (0.2, 0.4)],
            "left_excluded_points": [], "right_excluded_points": [],
            "toe_left": (-0.2, 0.2), "toe_right": (0.1, 0.2),
            "left_slope": 1.0, "right_slope": 2.0,
        }
        discharge_result["sections"] = []
        for physical_center in (0.0, 0.1, -0.1):
            section = dict(discharge_result)
            section.update({
                "section_center_physical": physical_center,
                "section_center": physical_center * ratio,
                "section_half_width": 0.05 * ratio,
                "base_left": discharge_result["toe_left"],
                "base_right": discharge_result["toe_right"],
            })
            discharge_result["sections"].append(section)
        profiles = [
            {"envelope": [(0.1, 0.9), (0.2, 0.8)], "center_envelope": [(0.1, 0.91)],
             "angle_deg": 30.0 + i, "intercept": 1.0, "slope": -0.5,
             "selected_particle_points": [], "outlier_points": []}
            for i in range(8)
        ]
        pair_profiles = []
        for i in range(4):
            pos, neg = profiles[i], profiles[i + 4]
            pair_profiles.append({
                "axis_angle_deg": 45.0 * i,
                "positive_profile": pos, "negative_profile": neg,
                "positive_angle_deg": pos["angle_deg"], "negative_angle_deg": neg["angle_deg"],
                "mean_angle_deg": (pos["angle_deg"] + neg["angle_deg"]) / 2,
                "positive_fit_points": pos["envelope"], "negative_fit_points": neg["envelope"],
            })
        repose_result = {
            "direction_profiles": profiles, "direction_angles_deg": [p["angle_deg"] for p in profiles],
            "pair_profiles": pair_profiles, "center_x": 0.0, "center_z": 0.0,
            "mean_angle_deg": 33.5, "stddev_deg": 2.291,
            "median_particle_radius": 0.01, "discarded_outliers": 0,
        }
        with tempfile.TemporaryDirectory(dir="/mnt/c/DEM") as temp:
            discharge = self._snapshot(Path(temp), "Ang Deslizamento", "deslizamento_final.vtk")
            repose = self._snapshot(Path(temp), "Ang Repouso", "repouso_10s.vtk")
            controller = object.__new__(CalibrationController)
            with patch.object(CalibrationController, "_mesh_bounds_for_snapshot", return_value={"x": (-1, 1), "z": (-1, 1)}), \
                 patch("optimizer.calibration_controller.measure_v_sliding_angle_three_cuts", return_value=discharge_result), \
                 patch("optimizer.calibration_controller.measure_pile_angle", return_value=repose_result):
                controller._measure_snapshot(discharge, "discharge")
                controller._measure_snapshot(repose, "repose")
            trial_dir = discharge.parents[2]
            discharge_svg = trial_dir / "medicao_descarga.svg"
            repose_svg = trial_dir / "medicao_repouso_oito_direcoes.svg"
            discharge_png = trial_dir / "medicao_descarga.png"
            repose_png = trial_dir / "medicao_repouso_oito_direcoes.png"
            ET.parse(discharge_svg)
            ET.parse(repose_svg)
            self.assertTrue(discharge_png.is_file())
            self.assertTrue(repose_png.is_file())
            repose_text = repose_svg.read_text(encoding="utf-8")
            self.assertEqual(repose_text.count("CORTE DIAMETRAL"), 4)
            self.assertEqual(repose_text.count("LADO +"), 4)
            self.assertEqual(repose_text.count("LADO -"), 4)
            self.assertIn("EIXO CENTRAL (0,0)", repose_text)
            self.assertIn("REFERENCIA HORIZONTAL Y=0", repose_text)
            self.assertIn("REFERENCIA HORIZONTAL", repose_text)
            self.assertIn("TRIAL", repose_text)

    def test_discharge_result_fields_follow_physical_cut_coordinates(self):
        measurement = {
            "sections": [
                {
                    "section_center_physical": center,
                    "mean_angle_deg": angle,
                    "left_angle_deg": angle - 1.0,
                    "right_angle_deg": angle + 1.0,
                    "band_clipped": center != 0.0,
                    "effective_band_physical": (max(center - 0.075, -0.12), min(center + 0.075, 0.12)),
                    "left_segments": [],
                    "right_segments": [],
                }
                for center, angle in ((-0.1, 10.0), (0.0, 20.0), (0.1, 30.0))
            ]
        }
        fields = discharge_cut_result_fields(measurement)
        self.assertEqual(fields["discharge_x0_deg"], 20.0)
        self.assertEqual(fields["discharge_xplus01_deg"], 30.0)
        self.assertEqual(fields["discharge_xminus01_deg"], 10.0)
        self.assertEqual(fields["discharge_x0_center_physical_m"], 0.0)
        self.assertEqual(fields["discharge_xplus01_center_physical_m"], 0.1)
        self.assertEqual(fields["discharge_xminus01_center_physical_m"], -0.1)
        self.assertEqual(fields["discharge_xminus01_band_min_m"], -0.12)
        self.assertTrue(fields["discharge_xplus01_band_clipped"])

    def test_unmeasurable_discharge_uses_official_three_cut_figure_without_angle(self):
        with tempfile.TemporaryDirectory() as temp:
            snapshot = self._snapshot(Path(temp), "Ang Deslizamento", "deslizamento_final.vtk")
            controller = object.__new__(CalibrationController)
            with (
                patch.object(CalibrationController, "_mesh_bounds_for_snapshot", return_value={"x": (-5, 5), "z": (-2, 2)}),
                patch("optimizer.calibration_controller.measure_v_sliding_angle_three_cuts", side_effect=UnmeasurableAngle("insufficient segment support on left pile face")),
            ):
                with self.assertRaises(MeasurementRejected) as raised:
                    controller._measure_snapshot(snapshot, "discharge")
            self.assertNotIn("arching", str(raised.exception).lower())
            trial_dir = snapshot.parents[2]
            official_svg = trial_dir / "medicao_descarga.svg"
            ET.parse(official_svg)
            text = official_svg.read_text(encoding="utf-8")
            self.assertEqual(text.count("CORTE X FISICO"), 3)
            self.assertIn("NO VALID ANGLE", text)
            self.assertIn("insufficient segment support", text.lower())
            self.assertNotIn("MEDIA DO CORTE:", text)
            self.assertNotIn("data-kind=\"fit-", text)
            self.assertTrue((trial_dir / "medicao_descarga.png").is_file())
            self.assertFalse((trial_dir / "medicao_descarga_diagnostico.png").exists())


if __name__ == "__main__":
    unittest.main()
