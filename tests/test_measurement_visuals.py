import math
import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path

from optimizer.measurements import measure_repose_eight_ray_least_squares, measure_v_sliding_angle_by_height
from optimizer.measurement_visuals import write_discharge_svg, write_repose_svg


class MeasurementVisualTests(unittest.TestCase):
    def test_repose_svg_contains_four_paired_axes_and_eight_individual_faces(self):
        points = []
        slope = math.tan(math.radians(30.0))
        for radial_index in range(1, 15):
            radius = radial_index / 15.0
            height = 1.0 - slope * radius
            for direction_index in range(72):
                theta = 2.0 * math.pi * direction_index / 72.0
                points.append((radius * math.cos(theta), height, radius * math.sin(theta)))
        measurement = measure_repose_eight_ray_least_squares(points, radial_bins=12, min_points_per_ray=8, particle_radii=[0.01] * len(points))
        with tempfile.TemporaryDirectory() as temp_dir:
            output = Path(temp_dir) / "repose.svg"
            write_repose_svg(output, measurement, "trial_0001", "radiusnbr=0.0065/2; nbrratio=v_radius1/v_radiusnbr")
            text = output.read_text(encoding="utf-8")
            svg_root = ET.parse(output).getroot()
        fit_line = next(element for element in svg_root.iter() if element.tag.endswith("line") and element.attrib.get("stroke") == "#d93535")
        screen_angle = abs(math.degrees(math.atan2(float(fit_line.attrib["y2"]) - float(fit_line.attrib["y1"]), float(fit_line.attrib["x2"]) - float(fit_line.attrib["x1"]))))
        self.assertAlmostEqual(screen_angle, measurement["direction_profiles"][0]["angle_deg"], delta=0.1)
        self.assertEqual(text.count("CORTE DIAMETRAL"), 4)
        self.assertEqual(text.count("LADO +"), 4)
        self.assertEqual(text.count("LADO -"), 4)
        self.assertEqual(len(measurement["direction_profiles"]), 8)
        self.assertIn('width="1800" height="1160"', text)
        self.assertIn("TRIAL 0001 - REPOUSO", text)
        self.assertIn("ANGULO DE CADA FACE RELATIVO A Y=0", text)
        self.assertIn("REFERENCIA HORIZONTAL Y=0", text)
        self.assertIn("MEDIA GLOBAL", text)
        self.assertIn("NBR:", text)
        self.assertIn("TODAS AS PARTICULAS: CINZA", text)
        self.assertIn("BASE ABAIXO DO CORTE", text)
        self.assertIn('data-kind="basal-reference"', text)
        basal = [element for element in svg_root.iter() if element.attrib.get("data-kind") == "basal-reference"]
        self.assertTrue(basal)
        self.assertTrue(all(float(element.attrib["fill-opacity"]) < 0.3 for element in basal))
        for angle in measurement["direction_angles_deg"]:
            self.assertIn(f"{angle:.2f} GRAUS", text)

    def test_discharge_svg_matches_v5_style_and_shows_both_base_anchored_faces(self):
        points = []
        radii = []
        radius = 0.01
        base_y = 0.05
        for sign, angle in ((-1.0, 30.0), (1.0, 45.0)):
            for index in range(12):
                y = base_y + 2 * radius * index
                outer = 0.2 + (y - base_y) / math.tan(math.radians(angle))
                points.append((0.0, y, sign * (outer - radius)))
                radii.append(radius)
        measurement = measure_v_sliding_angle_by_height(
            points, radii, horizontal_axis="z",
            mesh_bounds={"x": (-0.5, 0.5), "z": (-6.0, 6.0)},
        )
        with tempfile.TemporaryDirectory() as temp_dir:
            output = Path(temp_dir) / "discharge.svg"
            write_discharge_svg(output, points, measurement, "trial_0001", "radiusnbr=0.0065/2; nbrratio=v_radius1/v_radiusnbr")
            text = output.read_text(encoding="utf-8")
            svg_root = ET.parse(output).getroot()
        fit_lines = [element for element in svg_root.iter() if element.tag.endswith("line") and element.attrib.get("stroke") in {"#d93535", "#276bd2"}]
        self.assertEqual(len(fit_lines), 2)
        expected_angles = sorted((measurement["left_angle_deg"], measurement["right_angle_deg"]))
        screen_angles = sorted(math.degrees(math.atan2(abs(float(line.attrib["y2"]) - float(line.attrib["y1"])), abs(float(line.attrib["x2"]) - float(line.attrib["x1"])))) for line in fit_lines)
        for screen_angle, expected_angle in zip(screen_angles, expected_angles):
            self.assertAlmostEqual(screen_angle, expected_angle, delta=0.1)
        self.assertIn('width="1440" height="1000"', text)
        self.assertIn("TRIAL 0001 - DESLIZAMENTO", text)
        self.assertIn("CORTE CENTRAL", text)
        self.assertIn("FAIXAS DY=DIAMETRO MEDIO", text)
        self.assertIn("ALPHA", text)
        self.assertIn("AJUSTE ROBUSTO DE TODOS OS PONTOS DE PERFIL DA BASE AO TOPO", text)
        self.assertIn("INCLUIDOS COM PESO ROBUSTO MENOR", text)
        self.assertIn("ESQ ", text)
        self.assertIn("DIR ", text)
        self.assertIn("MEDIA ", text)
        self.assertIn("TRACEJADO: EIXO Z=0", text)
        self.assertIn("NBR", text)
        self.assertIn("#d93535", text)
        self.assertIn("#276bd2", text)

    def test_svg_to_png_renders_requested_dimensions_with_chrome(self):
        from optimizer.measurement_visuals import write_svg_png
        with tempfile.TemporaryDirectory(dir="/mnt/c/DEM") as temp_dir:
            source = Path(temp_dir) / "source.svg"
            target = Path(temp_dir) / "source.png"
            source.write_text('<svg xmlns="http://www.w3.org/2000/svg" width="1440" height="1000"><rect width="100%" height="100%" fill="#fff"/></svg>', encoding="utf-8")
            write_svg_png(source, target)
            payload = target.read_bytes()
        self.assertEqual(payload[:8], b"\x89PNG\r\n\x1a\n")
        import struct
        self.assertEqual(struct.unpack(">II", payload[16:24]), (1440, 1000))
        self.assertGreater(len(payload), 1000)


if __name__ == "__main__":
    unittest.main()
