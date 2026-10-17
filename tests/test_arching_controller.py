import tempfile
import unittest
from pathlib import Path

from optimizer.calibration_controller import CalibrationController, MeasurementRejected


class ArchingControllerDiagnosticTests(unittest.TestCase):
    def test_successful_linear_fit_is_not_rejected_by_arching_precheck(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            project = Path(temp_dir) / "project"
            trial = Path(temp_dir) / "run" / "trial_0001"
            post = trial / "Ang Deslizamento" / "post"
            post.mkdir(parents=True)
            (trial / "Ang Deslizamento" / "in.deslizamento").write_text(
                "variable radius1 equal 0.250/2\nvariable radiusnbr equal 0.0065/2\nvariable nbrratio equal v_radius1/v_radiusnbr\n",
                encoding="utf-8",
            )
            points = []
            ratio = 0.125 / 0.00325
            for offset in (0.0, 0.1, -0.1):
                center = offset * ratio
                for index in range(10):
                    y = (index + 0.5) * 0.2
                    distance = 1.0 + index * 0.2
                    points.extend(((center, y, -distance), (center, y, distance)))
            coords = " ".join(str(value) for point in points for value in point)
            radii = " ".join("0.1" for _ in points)
            vtk = f"# vtk DataFile Version 3.0\nfixture\nASCII\nDATASET POLYDATA\nPOINTS {len(points)} double\n{coords}\nPOINT_DATA {len(points)}\nFIELD FieldData 1\nradius 1 {len(points)} double\n{radii}\n"
            snapshot = post / "deslizamento_final.vtk"
            snapshot.write_text(vtk, encoding="utf-8")
            (post / "box_100.stl").write_text(
                "solid box\nvertex -4.8076923077 0 -9.6153846154\nvertex 4.8076923077 0 9.6153846154\nendsolid\n",
                encoding="utf-8",
            )
            (project / "Ang Deslizamento").mkdir(parents=True)
            (project / "Ang Deslizamento" / "in.deslizamento").write_text("", encoding="utf-8")
            (project / "configuracao_otimizacao.txt").write_text('{"parameters":{"sf11":{"min":0.1,"max":1.0}}}', encoding="utf-8")
            controller = CalibrationController(project, run_root=Path(temp_dir) / "run", solver_command=["true"])
            measured = controller._measure_snapshot(snapshot, "discharge")
            self.assertAlmostEqual(measured["mean_angle_deg"], 45.0, delta=0.5)
            self.assertTrue((trial / "medicao_descarga.svg").is_file())

    def test_diverging_inner_faces_without_fit_support_are_not_labeled_arching(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            project = Path(temp_dir) / "project"
            trial = Path(temp_dir) / "run" / "trial_0001"
            post = trial / "Ang Deslizamento" / "post"
            post.mkdir(parents=True)
            (trial / "Ang Deslizamento" / "in.deslizamento").write_text(
                "variable radius1 equal 0.250/2\nvariable radiusnbr equal 0.0065/2\nvariable nbrratio equal v_radius1/v_radiusnbr\n",
                encoding="utf-8",
            )
            # A narrow central opening and two continuous inner faces that
            # diverge outward as height rises; no particles bridge the opening.
            distances = [1.0, 0.8, 0.6, 0.4, 0.6, 0.8]
            points = []
            ratio = 0.125 / 0.00325
            for offset in (0.0, 0.1, -0.1):
                center = offset * ratio
                for index, distance in enumerate(distances):
                    y = (index + 0.5) * 0.2
                    points.extend(((center, y, -distance), (center, y, distance)))
            coords = " ".join(str(value) for point in points for value in point)
            radii = " ".join("0.1" for _ in points)
            vtk = f"# vtk DataFile Version 3.0\nfixture\nASCII\nDATASET POLYDATA\nPOINTS {len(points)} double\n{coords}\nPOINT_DATA {len(points)}\nFIELD FieldData 1\nradius 1 {len(points)} double\n{radii}\n"
            snapshot = post / "deslizamento_final.vtk"
            snapshot.write_text(vtk, encoding="utf-8")
            (post / "box_100.stl").write_text(
                "solid box\nvertex -4.8076923077 0 -9.6153846154\nvertex 4.8076923077 0 9.6153846154\nendsolid\n",
                encoding="utf-8",
            )
            (project / "Ang Deslizamento").mkdir(parents=True)
            (project / "Ang Deslizamento" / "in.deslizamento").write_text("", encoding="utf-8")
            (project / "configuracao_otimizacao.txt").write_text('{"parameters":{"sf11":{"min":0.1,"max":1.0}}}', encoding="utf-8")
            controller = CalibrationController(project, run_root=Path(temp_dir) / "run", solver_command=["true"])
            with self.assertRaises(MeasurementRejected) as raised:
                controller._measure_snapshot(snapshot, "discharge")
            self.assertIn("angle unavailable", str(raised.exception).lower())
            self.assertNotIn("arching", str(raised.exception).lower())


if __name__ == "__main__":
    unittest.main()
