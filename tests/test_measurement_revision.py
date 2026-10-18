import math
import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path

from optimizer.measurements import (
    measure_repose_eight_ray_least_squares,
    measure_v_sliding_angle_by_height,
    measure_v_sliding_angle_three_cuts,
)
from optimizer.measurement_visuals import write_discharge_svg, write_repose_svg, write_measurement_diagnostic_svg, write_repose_cutoff_diagnostic_svg, write_repose_failure_diagnostic_svg
from optimizer.calibration_controller import CalibrationController, discharge_cut_result_fields, RESULT_FIELDS, repose_pair_result_fields


class MeasurementRevisionTests(unittest.TestCase):
    @staticmethod
    def _repose_pile():
        points = []
        radii = []
        slope = math.tan(math.radians(30.0))
        max_radius = 0.48
        for radial_index in range(1, 25):
            radius = radial_index * max_radius / 24
            height = slope * (max_radius - radius)
            for direction_index in range(72):
                theta = 2.0 * math.pi * direction_index / 72.0
                points.append((radius * math.cos(theta), height, radius * math.sin(theta)))
                radii.append(0.01)
        return points, radii

    @staticmethod
    def _three_cut_cloud():
        points = []
        radii = []
        scale = 20.0
        for physical_center in (-0.1, 0.0, 0.1):
            section_center = physical_center * scale
            for sign, angle in ((-1.0, 30.0), (1.0, 45.0)):
                for index in range(12):
                    radius = 0.01
                    y = 0.05 + 2 * radius * index
                    outer = 0.2 + (y - 0.05) / math.tan(math.radians(angle))
                    points.append((section_center, y, sign * (outer - radius)))
                    radii.append(radius)
        return points, radii

    def test_discharge_uses_three_physical_centers_and_scaled_half_width(self):
        points, radii = self._three_cut_cloud()
        result = measure_v_sliding_angle_three_cuts(
            points,
            radii,
            horizontal_axis="z",
            mesh_bounds={"x": (-5.0, 5.0), "z": (-5.0, 5.0)},
            nbrratio=20.0,
        )
        self.assertEqual(len(result["sections"]), 3)
        self.assertEqual([section["section_center_physical"] for section in result["sections"]], [-0.1, 0.0, 0.1])
        self.assertEqual([section["section_center"] for section in result["sections"]], [-2.0, 0.0, 2.0])
        self.assertEqual([section["section_half_width"] for section in result["sections"]], [1.5, 1.5, 1.5])
        self.assertEqual([section["section_half_width_physical"] for section in result["sections"]], [0.075, 0.075, 0.075])
        for section in result["sections"]:
            self.assertAlmostEqual(section["left_angle_deg"], 30.0, delta=1.0)
            self.assertAlmostEqual(section["right_angle_deg"], 45.0, delta=1.0)
            self.assertEqual(len(section["left_selected_particle_centers"]), len(section["left_fit_points"]))
            self.assertEqual(len(section["right_selected_particle_centers"]), len(section["right_fit_points"]))
        self.assertAlmostEqual(result["mean_angle_deg"], 37.5, delta=1.0)

    def test_discharge_applies_repose_minimum_height_and_preserves_two_piece_fits(self):
        points, radii = self._three_cut_cloud_with_reversal()
        result = measure_v_sliding_angle_three_cuts(
            points, radii, horizontal_axis="z",
            mesh_bounds={"x": (-5.0, 5.0), "z": (-5.0, 5.0)}, nbrratio=20.0,
        )
        self.assertEqual(result["minimum_height_particle_diameters"], 2.0)
        for section in result["sections"]:
            for side in ("left", "right"):
                self.assertAlmostEqual(section[f"profile_start_height_{side}"], 0.09)
                self.assertTrue(all(y >= section[f"profile_start_height_{side}"] for _, y in section[f"{side}_selected_particle_centers"]))
                self.assertEqual(section[f"{side}_segment_count"], 2)
                self.assertEqual(len(section[f"{side}_segments"]), 2)

    def test_discharge_sustained_reversal_reports_two_segment_angles_and_scalar(self):
        points, radii = self._three_cut_cloud_with_reversal()
        result = measure_v_sliding_angle_three_cuts(
            points, radii, horizontal_axis="z",
            mesh_bounds={"x": (-5.0, 5.0), "z": (-5.0, 5.0)}, nbrratio=20.0,
        )
        for section in result["sections"]:
            self.assertEqual(section["left_segment_count"], 2)
            self.assertEqual(section["right_segment_count"], 2)
            for side, expected in (("left", (30.0, 45.0)), ("right", (30.0, 45.0))):
                segments = section[f"{side}_segments"]
                self.assertEqual(len(segments), 2)
                self.assertEqual([item["fit_support"] for item in segments], [6, 5])
                for segment, angle in zip(segments, expected):
                    self.assertAlmostEqual(segment["angle_deg"], angle, delta=1.5)
        self.assertAlmostEqual(result["sections"][0]["mean_angle_deg"], 37.5, delta=1.5)
        csv_fields = discharge_cut_result_fields(result)
        self.assertIn("discharge_x0_left_segment1_deg", RESULT_FIELDS)
        self.assertIn("discharge_x0_left_segment1_count", RESULT_FIELDS)
        self.assertEqual(csv_fields["discharge_x0_left_segment1_count"], 2)
        self.assertIn("discharge_x0_right_segment2_support", RESULT_FIELDS)
        self.assertAlmostEqual(csv_fields["discharge_x0_left_segment1_deg"], result["sections"][0]["left_segments"][0]["angle_deg"])
        self.assertEqual(csv_fields["discharge_x0_left_segment2_support"], 5)
        self.assertEqual(csv_fields["discharge_x0_minimum_height_particle_diameters"], 2.0)
        self.assertEqual(csv_fields["discharge_x0_left_profile_start_height_m"], 0.09)
        self.assertIn("discharge_x0_left_profile_start_height_m", RESULT_FIELDS)
        with tempfile.TemporaryDirectory() as temp_dir:
            output = Path(temp_dir) / "two-segment.svg"
            write_discharge_svg(output, points, result, "trial", "nbrratio=20")
            svg = output.read_text(encoding="utf-8")
            root = ET.parse(output).getroot()
        fit_lines = [el for el in root.iter() if el.attrib.get("data-kind", "").startswith("fit-") and "segment" in el.attrib.get("data-kind", "")]
        labels = [el for el in root.iter() if el.attrib.get("data-kind") == "segment-label"]
        self.assertEqual(len(fit_lines), 12)
        self.assertEqual(len(labels), 6)
        self.assertTrue(all("SEG1" in el.text and "SEG2" in el.text for el in labels))
        self.assertIn("MEDIA DO CORTE", svg)

    def test_discharge_reversal_without_two_supported_segments_is_unmeasurable(self):
        points, radii = self._three_cut_cloud_with_reversal(bands=8)
        with self.assertRaisesRegex(ValueError, "segment support"):
            measure_v_sliding_angle_three_cuts(
                points, radii, horizontal_axis="z",
                mesh_bounds={"x": (-5.0, 5.0), "z": (-5.0, 5.0)}, nbrratio=20.0,
            )

    @staticmethod
    def _three_cut_cloud_with_reversal(bands=12):
        points, radii = [], []
        for physical_center in (-0.1, 0.0, 0.1):
            center_x = physical_center * 20.0
            for sign in (-1.0, 1.0):
                for index in range(bands):
                    y = 0.05 + 0.02 * index
                    if index <= bands // 2:
                        outward = 0.25 + (y - 0.05) / math.tan(math.radians(30.0))
                    else:
                        crest_y = 0.05 + 0.02 * (bands // 2)
                        crest = 0.25 + (crest_y - 0.05) / math.tan(math.radians(30.0))
                        outward = crest - (y - crest_y) / math.tan(math.radians(45.0))
                    points.append((center_x, y, sign * (outward - 0.01)))
                    radii.append(0.01)
        return points, radii

    def test_discharge_reports_outer_bands_clipped_by_mesh_bounds(self):
        points, radii = self._three_cut_cloud()
        result = measure_v_sliding_angle_three_cuts(
            points,
            radii,
            horizontal_axis="z",
            mesh_bounds={"x": (-2.4, 2.4), "z": (-5.0, 5.0)},
            nbrratio=20.0,
        )
        clipped_by_center = {section["section_center_physical"]: section["band_clipped"] for section in result["sections"]}
        self.assertTrue(clipped_by_center[-0.1])
        self.assertFalse(clipped_by_center[0.0])
        self.assertTrue(clipped_by_center[0.1])
        bands_by_center = {section["section_center_physical"]: section["effective_band_physical"] for section in result["sections"]}
        self.assertEqual(bands_by_center[-0.1], (-0.12, -0.025))
        self.assertEqual(bands_by_center[0.0], (-0.075, 0.075))
        self.assertEqual(bands_by_center[0.1], (0.025, 0.12))
        self.assertLess(bands_by_center[0.0][0], bands_by_center[0.1][0])
        self.assertGreater(bands_by_center[0.0][1], bands_by_center[0.1][0])
        with tempfile.TemporaryDirectory() as temp_dir:
            output = Path(temp_dir) / "clipped.svg"
            write_discharge_svg(output, points, result, "trial", "nbrratio=20")
            svg = output.read_text(encoding="utf-8")
        self.assertIn("BANDA EFETIVA", svg)
        self.assertIn("CLIPADA PELA MALHA", svg)
        self.assertEqual(svg.count('data-kind="minimum-height-cutoff"'), 6)
        self.assertIn("FIT ESQ Y=", svg)

    def test_clipped_discharge_band_excludes_particles_beyond_scaled_mesh_bounds(self):
        points, radii = self._three_cut_cloud()
        expected = measure_v_sliding_angle_three_cuts(
            points, radii, horizontal_axis="z",
            mesh_bounds={"x": (-2.4, 2.4), "z": (-5.0, 5.0)}, nbrratio=20.0,
        )
        contaminated_points = list(points)
        contaminated_radii = list(radii)
        for index in range(12):
            y = 0.05 + 2 * 0.01 * index
            for outside_x in (-2.6, 2.6):
                for sign in (-1.0, 1.0):
                    contaminated_points.append((outside_x, y, sign * 0.01))
                    contaminated_radii.append(0.01)
        actual = measure_v_sliding_angle_three_cuts(
            contaminated_points, contaminated_radii, horizontal_axis="z",
            mesh_bounds={"x": (-2.4, 2.4), "z": (-5.0, 5.0)}, nbrratio=20.0,
        )
        for expected_section, actual_section in zip(expected["sections"], actual["sections"]):
            self.assertEqual(actual_section["section_center_physical"], expected_section["section_center_physical"])
            self.assertAlmostEqual(actual_section["left_angle_deg"], expected_section["left_angle_deg"], places=8)
            self.assertAlmostEqual(actual_section["right_angle_deg"], expected_section["right_angle_deg"], places=8)
            self.assertEqual(actual_section["section_normal_bounds_sim"], actual_section["effective_band_sim"])
        with tempfile.TemporaryDirectory() as temp_dir:
            output = Path(temp_dir) / "clipped-effective-bands.svg"
            write_discharge_svg(output, contaminated_points, actual, "trial", "nbrratio=20")
            svg = output.read_text(encoding="utf-8")
        self.assertEqual(svg.count('data-kind="particle"'), 60)

    def test_repose_selection_uses_outer_face_not_nearest_centerline_particle(self):
        points = []
        radii = []
        radius = 0.01
        slope = math.tan(math.radians(30.0))
        for index in range(13):
            y = 0.01 + 2.0 * radius * index
            outer = 0.8 - (y - 0.01) / slope
            for side in (-1.0, 1.0):
                points.append((0.0, y, side * outer))
                points.append((0.0, y, side * 0.1))
                radii.extend((radius, radius))
        measured = measure_v_sliding_angle_by_height(
            points,
            radii,
            horizontal_axis="z",
            mesh_bounds={"x": (-radius, radius), "z": (-1.0, 1.0)},
            section_half_width=radius,
            section_center=0.0,
            profile_selection="outer_envelope",
            min_fit_points=8,
        )
        self.assertAlmostEqual(measured["left_angle_deg"], 30.0, delta=1.0)
        self.assertAlmostEqual(measured["right_angle_deg"], 30.0, delta=1.0)
        self.assertGreater(min(abs(z) for z, _ in measured["left_selected_particle_centers"]), 0.3)
        self.assertGreater(min(abs(z) for z, _ in measured["right_selected_particle_centers"]), 0.3)

    def test_base_edge_nearest_repose_uses_huber_for_a_single_reversal_outlier(self):
        points = []
        radii = []
        radius = 0.01
        slope = math.tan(math.radians(30.0))
        for band in range(36):
            y = 0.01 + band * 1.5 * radius
            z = 0.50 - (y - 0.01) / slope
            points.append((0.0, y, z))
            radii.append(radius)
            if 10 <= band <= 12:
                # One selected sample lies near the basal-edge axis and creates
                # a local reversal; Huber should suppress it, not reject the face.
                points.append((0.0, y + 0.002, 0.49))
                radii.append(radius)
        result = measure_v_sliding_angle_by_height(
            points,
            radii,
            horizontal_axis="z",
            mesh_bounds={"x": (-radius, radius), "z": (-0.5, 0.5)},
            section_half_width=radius,
            section_center=0.0,
            profile_selection="base_edge_nearest",
            height_band_factor=1.5,
            minimum_height_particle_diameters=2.0,
            min_fit_points=8,
        )
        self.assertAlmostEqual(result["right_angle_deg"], 30.0, delta=2.0)
        self.assertGreater(len(result["right_downweighted_points"]), 0)
        self.assertLess(result["base_right"][1], result["profile_start_height_right"])
        self.assertAlmostEqual(result["base_right"][0], 0.51, places=3)
        self.assertTrue(all(y >= result["profile_start_height_right"] for _, y in result["right_selected_particle_centers"]))

    def test_repose_uses_discharge_base_anchored_robust_fit_on_four_diameters(self):
        points = []
        radii = []
        slope = math.tan(math.radians(30.0))
        max_radius = 0.48
        for radial_index in range(1, 25):
            radius = radial_index * max_radius / 24
            height = slope * (max_radius - radius)
            for direction_index in range(72):
                theta = 2.0 * math.pi * direction_index / 72.0
                points.append((radius * math.cos(theta), height, radius * math.sin(theta)))
                radii.append(0.01)
        result = measure_repose_eight_ray_least_squares(points, particle_radii=radii, height_band_factor=1.5)
        baseline = measure_repose_eight_ray_least_squares(points, particle_radii=radii, height_band_factor=2.0)
        self.assertEqual(result["profile_method"], "repose_horizontal_surface_consensus")
        self.assertEqual(result["repose_axis_reference"], "section_external_vertical")
        self.assertEqual(result["axis_angles_deg"], [0.0, 45.0, 90.0, 135.0])
        self.assertEqual(len(result["pair_profiles"]), 4)
        self.assertEqual(len(result["direction_angles_deg"]), 8)
        for face_index, profile in enumerate(result["direction_profiles"]):
            self.assertEqual(profile["profile_selection"], "base_edge_nearest")
            self.assertAlmostEqual(profile["height_band_size"], 1.5 * result["median_particle_radius"])
            self.assertGreater(profile["fit_support"], baseline["direction_profiles"][face_index]["fit_support"])
            first_x, first_y, first_z = profile["selected_particle_points"][0]
            direction = math.radians(profile["direction_angle_deg"])
            first_projection = first_x * math.cos(direction) + first_z * math.sin(direction)
            self.assertGreater(abs(profile["base_axis_coordinate"]), abs(first_projection))
            self.assertAlmostEqual(profile["angle_deg"], 30.0, delta=1.0)
            self.assertGreaterEqual(profile["fit_support"], 5)
            self.assertLess(profile["base"][1], profile["profile_start_height"])
            self.assertGreaterEqual(
                min(y for _, y in profile["selected_particle_centers"]),
                profile["profile_start_height"],
            )
            fit_start_h, fit_start_y = profile["fit_start"]
            self.assertGreaterEqual(fit_start_y, profile["profile_start_height"])
            self.assertAlmostEqual(
                profile["intercept"], fit_start_y - profile["slope"] * abs(fit_start_h), places=9
            )
            self.assertAlmostEqual(profile["envelope"][0][1], fit_start_y, places=9)
        self.assertEqual(result["horizontal_reference_y"], 0.0)
        self.assertAlmostEqual(result["mean_angle_deg"], 30.0, delta=1.0)

    def test_repose_pairs_opposite_faces_with_the_discharge_profile_fit(self):
        points, radii = self._repose_pile()
        result = measure_repose_eight_ray_least_squares(
            points, radial_bins=12, min_points_per_ray=8, particle_radii=radii
        )
        self.assertEqual((result["center_x"], result["center_z"]), (0.0, 0.0))
        self.assertEqual(result["profile_method"], "repose_horizontal_surface_consensus")
        self.assertEqual(len(result["pair_profiles"]), 4)
        self.assertEqual(len(result["direction_profiles"]), 8)
        for pair in result["pair_profiles"]:
            self.assertAlmostEqual(pair["positive_angle_deg"], 30.0, delta=1.0)
            self.assertAlmostEqual(pair["negative_angle_deg"], 30.0, delta=1.0)
            self.assertAlmostEqual(pair["mean_angle_deg"], 30.0, delta=1.0)
            self.assertGreaterEqual(pair["positive_profile"]["fit_support"], 8)
            self.assertGreaterEqual(pair["negative_profile"]["fit_support"], 8)

    def test_repose_svg_groups_four_diameters_with_both_sides_and_zero_references(self):
        points, radii = self._repose_pile()
        measurement = measure_repose_eight_ray_least_squares(
            points, radial_bins=12, min_points_per_ray=8, particle_radii=radii
        )
        measurement["direction_profiles"][0]["downweighted_particle_points"] = [measurement["direction_profiles"][0]["selected_particle_points"][0]]
        measurement["direction_profiles"][0]["downweighted_count"] = 1
        measurement["discarded_particle_points"] = [(0.0, -0.01, 0.0), (0.0, 0.2, 0.0)]
        for pair in measurement["pair_profiles"]:
            pair["section_discarded_particle_points"] = measurement["discarded_particle_points"]
        with tempfile.TemporaryDirectory() as temp_dir:
            output = Path(temp_dir) / "repose.svg"
            write_repose_svg(output, measurement, "trial_0045", "nbrratio=25", points=points)
            svg = output.read_text(encoding="utf-8")
            root = ET.parse(output).getroot()
        self.assertEqual(svg.count("CORTE DIAMETRAL"), 4)
        self.assertEqual(svg.count("LADO +"), 4)
        self.assertEqual(svg.count("LADO -"), 4)
        self.assertIn("EIXO CENTRAL (0,0)", svg)
        self.assertIn("EIXOS EXTERNOS DA SECAO", svg)
        self.assertIn("PASSO EM Y: 0.25xRAIO 3D MEDIANO", svg)
        self.assertIn("CORTE HORIZONTAL: 2.0 DIAMETROS; AJUSTE COMECA ACIMA DO CORTE", svg)
        self.assertIn("REFERENCIA HORIZONTAL Y=0", svg)
        self.assertIn("CINZA: ESFERAS DA SECAO", svg)
        self.assertIn("BASE TRANSLUCIDA: FORA DO FIT", svg)
        self.assertIn("X LARANJA: EXCLUIDOS DO FIT", svg)
        toe_axes = [element for element in root.iter() if element.tag.endswith("line") and element.attrib.get("data-kind") == "toe-axis"]
        self.assertEqual(len(toe_axes), 8)
        selected_markers = [element for element in root.iter() if element.attrib.get("data-kind") == "selected"]
        self.assertEqual(len(selected_markers), sum(profile["candidate_count"] for profile in measurement["direction_profiles"]))
        outlier_markers = [element for element in root.iter() if element.attrib.get("data-kind") == "outlier"]
        self.assertGreaterEqual(len(outlier_markers), 1)
        cutoff_lines = [element for element in root.iter() if element.attrib.get("data-kind") == "minimum-height-cutoff"]
        self.assertEqual(len(cutoff_lines), 8)
        paired_fit_lines = [el for el in root.iter() if el.attrib.get("data-kind", "").startswith("fit-")]
        self.assertTrue(paired_fit_lines)
        self.assertTrue(all(el.attrib.get("clip-path") for el in paired_fit_lines))
        self.assertEqual({el.attrib["data-cutoff-side"] for el in cutoff_lines}, {"left", "right"})
        for fit_line in paired_fit_lines:
            cutoff_y = float(fit_line.attrib["data-cutoff-screen-y"])
            self.assertLessEqual(float(fit_line.attrib["y1"]), cutoff_y + 1e-6)
            self.assertLessEqual(float(fit_line.attrib["y2"]), cutoff_y + 1e-6)
        particle_markers = [element for element in root.iter() if element.attrib.get("data-kind") == "particle"]
        self.assertTrue(particle_markers)
        global_outliers = [element for element in root.iter() if element.attrib.get("data-kind") == "global-outlier"]
        self.assertTrue(global_outliers)
        cutoff_markers = particle_markers + selected_markers + outlier_markers + global_outliers
        self.assertTrue(all(float(el.attrib["data-particle-y"]) >= float(el.attrib["data-cutoff-y"]) for el in cutoff_markers if "data-particle-y" in el.attrib))
        self.assertEqual(len(global_outliers), 4)  # only the above-cutoff diagnostic, once per panel
        ordered_nodes = list(root.iter())
        selected_index = next(index for index, el in enumerate(ordered_nodes) if el.attrib.get("data-kind") == "selected")
        outlier_index = next(index for index, el in enumerate(ordered_nodes) if el.attrib.get("data-kind") == "outlier")
        self.assertGreater(outlier_index, selected_index)
        circles = [el for el in root.iter() if el.tag.endswith("circle")]
        self.assertGreater(len(circles), sum(pair["section_particle_count"] for pair in measurement["pair_profiles"]))

    def test_discharge_svg_combines_all_three_cuts_in_one_figure(self):
        points, radii = self._three_cut_cloud()
        measurement = measure_v_sliding_angle_three_cuts(
            points,
            radii,
            horizontal_axis="z",
            mesh_bounds={"x": (-5.0, 5.0), "z": (-5.0, 5.0)},
            nbrratio=20.0,
        )
        measurement["sections"][0]["left_downweighted_points"] = [measurement["sections"][0]["left_selected_particle_centers"][0]]
        with tempfile.TemporaryDirectory() as temp_dir:
            output = Path(temp_dir) / "discharge.svg"
            write_discharge_svg(output, points, measurement, "trial_0045", "nbrratio=20")
            svg = output.read_text(encoding="utf-8")
            root = ET.parse(output).getroot()
        self.assertEqual(svg.count("CORTE X FISICO"), 3)
        self.assertIn("3 CORTES FISICOS: X=-0.1 M, X=+0.0 M, X=+0.1 M", svg)
        self.assertIn("MEIA-LARGURA=0.075 M", svg)
        self.assertIn("X=-0.100 M", svg)
        self.assertIn("X=+0.000 M", svg)
        self.assertIn("X=+0.100 M", svg)
        fit_colors = {"#d93535", "#276bd2"}
        fit_lines = [el for el in root.iter() if el.tag.endswith("line") and el.attrib.get("data-kind", "").startswith("fit-") and el.attrib.get("stroke") in fit_colors]
        self.assertEqual(len(fit_lines), 6)
        cutoff_lines = [el for el in root.iter() if el.attrib.get("data-kind") == "minimum-height-cutoff"]
        self.assertEqual(len(cutoff_lines), 6)
        self.assertIn("FIT ESQ Y=", svg)
        self.assertIn("BASE TRANSLUCIDA", svg)
        self.assertIn("BASE TRANSLUCIDA: REFERENCIA GEOMETRICA; FORA DOS FITS", svg)
        ordered_nodes = list(root.iter())
        selected_index = next(index for index, el in enumerate(ordered_nodes) if el.attrib.get("data-kind") == "selected")
        outlier_index = next(index for index, el in enumerate(ordered_nodes) if el.attrib.get("data-kind") == "outlier")
        self.assertGreater(outlier_index, selected_index)

    def test_repose_pair_fields_preserve_four_means_and_eight_face_angles(self):
        points, radii = self._repose_pile()
        measurement = measure_repose_eight_ray_least_squares(points, particle_radii=radii)
        fields = repose_pair_result_fields(measurement)
        self.assertEqual(len(fields), 12)
        for pair in measurement["pair_profiles"]:
            axis = int(pair["axis_angle_deg"])
            self.assertAlmostEqual(fields[f"repose_axis{axis}_deg"], pair["mean_angle_deg"])
            self.assertAlmostEqual(fields[f"repose_axis{axis}_plus_deg"], pair["positive_angle_deg"])
            self.assertAlmostEqual(fields[f"repose_axis{axis}_minus_deg"], pair["negative_angle_deg"])

    def test_repose_failure_diagnostic_hides_all_particles_below_cutoff(self):
        points = [(0.0, 0.01, 0.0), (0.1, 0.05, 0.0), (0.2, 0.12, 0.0)]
        with tempfile.TemporaryDirectory() as temp_dir:
            output = Path(temp_dir) / "diagnostic.svg"
            write_measurement_diagnostic_svg(
                output, points, "repose", "trial", "no valid face", minimum_height=0.05
            )
            root = ET.parse(output).getroot()
        particles = [el for el in root.iter() if el.attrib.get("data-kind") == "diagnostic-particle"]
        self.assertEqual([float(el.attrib["data-particle-y"]) for el in particles], [0.05, 0.12])
        cutoff = [el for el in root.iter() if el.attrib.get("data-kind") == "minimum-height-cutoff"]
        self.assertEqual(len(cutoff), 1)

    def test_standard_four_cut_diagnostic_hides_particles_below_cutoff(self):
        points = [(0.0, 0.01, 0.0), (0.1, 0.05, 0.0), (0.2, 0.12, 0.0)]
        with tempfile.TemporaryDirectory() as temp_dir:
            output = Path(temp_dir) / "repose-four-cuts.svg"
            write_repose_cutoff_diagnostic_svg(output, points, "trial_0045", 0.05, "reversal")
            root = ET.parse(output).getroot()
        particles = [el for el in root.iter() if el.attrib.get("data-kind") == "diagnostic-particle"]
        self.assertEqual([float(el.attrib["data-particle-y"]) for el in particles], [0.05, 0.12] * 4)
        self.assertEqual(len([el for el in root.iter() if el.attrib.get("data-kind") == "minimum-height-cutoff"]), 4)
        self.assertIn("SEM ANGULO VALIDO", "".join(root.itertext()))

    def test_repose_failure_uses_standard_four_cut_layout_without_angle(self):
        points = [(0.0, 0.01, 0.0), (0.1, 0.05, 0.0), (-0.1, 0.08, 0.0)]
        with tempfile.TemporaryDirectory() as temp_dir:
            output = Path(temp_dir) / "repose-failure.svg"
            write_repose_failure_diagnostic_svg(output, points, "trial_0022", "insufficient height bands on left face")
            root = ET.parse(output).getroot()
            svg = output.read_text(encoding="utf-8")
        self.assertEqual(svg.count("CORTE DIAMETRAL"), 4)
        self.assertEqual(len([el for el in root.iter() if el.attrib.get("data-kind") == "particle"]), 12)
        self.assertIn("SEM ANGULO VALIDO", "".join(root.itertext()))
        self.assertIn("insufficient height bands on left face", "".join(root.itertext()))
        self.assertNotIn("ALPHA", svg)
        self.assertNotIn("data-kind=\"fit", svg)

    def test_deck_nbrratio_expression_is_evaluated(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            deck = Path(temp_dir) / "in.deslizamento"
            deck.write_text(chr(10).join([
                "variable radius1 equal 0.125",
                "variable radiusnbr equal 0.0065/2",
                "variable nbrratio equal v_radius1/v_radiusnbr",
            ]) + chr(10), encoding="utf-8")
            self.assertAlmostEqual(CalibrationController._nbr_ratio_from_deck(deck), 0.125 / 0.00325)


if __name__ == "__main__":
    unittest.main()
