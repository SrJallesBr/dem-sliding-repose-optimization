import copy
import math
import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path

from optimizer.measurements import measure_v_sliding_angle_by_height, measure_v_sliding_angle_three_cuts
from optimizer.measurement_visuals import write_discharge_svg, write_discharge_failure_svg


def cloud(centers):
    points = []
    for center in centers:
        for sign, angle in ((-1, 30), (1, 45)):
            for i in range(16):
                y = 0.05 + 0.02 * i
                outer = 0.2 + (y - 0.05) / math.tan(math.radians(angle))
                points.append((center, y, sign * (outer - 0.01)))
    return points, [0.01] * len(points)


class DischargeParticleRadiiTests(unittest.TestCase):
    def test_three_cuts_render_individual_projected_radii_without_changing_measurement(self):
        points, radii = cloud((-2, 0, 2))
        result = measure_v_sliding_angle_three_cuts(points, radii, horizontal_axis='z',
                    mesh_bounds={'x': (-5, 5), 'z': (-5, 5)}, nbrratio=20)
        original = copy.deepcopy(result)
        display_radii = [0.01 if i % 2 else 0.02 for i in range(len(points))]
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / 'discharge.svg'
            write_discharge_svg(path, points, result, 'trial', particle_radii=display_radii)
            root = ET.parse(path).getroot()
        spheres = [e for e in root.iter() if e.get('data-kind') == 'sphere-projection']
        self.assertEqual(len(spheres), len(points))
        self.assertEqual({float(e.get('data-radius-sim')) for e in spheres}, {0.01, 0.02})
        for sphere in spheres:
            index = int(sphere.get('data-particle-index'))
            self.assertEqual(float(sphere.get('data-radius-sim')), display_radii[index])
            self.assertAlmostEqual(float(sphere.get('r')), display_radii[index] * float(sphere.get('data-px-per-unit')), delta=0.011)
            if sphere.get('data-basal-reference') == 'true':
                self.assertLess(float(sphere.get('stroke-opacity')), 0.3)
        self.assertTrue(any(e.get('data-basal-reference') == 'true' for e in spheres))
        self.assertEqual(result, original)


    def test_single_cut_also_uses_particle_radii(self):
        points, radii = cloud((0,))
        result = measure_v_sliding_angle_by_height(points, radii, horizontal_axis='z',
                    mesh_bounds={'x': (-5, 5), 'z': (-5, 5)})
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / 'single.svg'
            write_discharge_svg(path, points, result, 'trial', particle_radii=radii)
            spheres = [e for e in ET.parse(path).getroot().iter() if e.get('data-kind') == 'sphere-projection']
        self.assertEqual(len(spheres), len(points))


    def test_single_cut_renders_only_effective_slab_with_matching_particle_radii(self):
        base_points, _ = cloud((0,))
        for axis in ('z', 'x'):
            with self.subTest(horizontal_axis=axis):
                points, radii = [], []
                for i, (_, y, horizontal) in enumerate(base_points):
                    # Interleave nominal-only centers so eligible indices are
                    # not a prefix; include both effective slab boundaries.
                    for normal, radius in ((0.08, 0.04), ((-0.025, 0.0, 0.025)[i % 3], 0.01 + 0.001 * (i % 2))):
                        points.append((normal, y, horizontal) if axis == 'z' else (horizontal, y, normal))
                        radii.append(radius)
                result = measure_v_sliding_angle_by_height(
                    points, radii, horizontal_axis=axis,
                    mesh_bounds={'x': (-5, 5), 'z': (-5, 5)},
                    section_half_width=0.1, section_center=0.0,
                    section_normal_bounds_sim=(-0.025, 0.025))
                original = copy.deepcopy(result)
                expected_indices = list(range(1, len(points), 2))
                self.assertEqual(result['section_particle_count'], len(expected_indices))
                with tempfile.TemporaryDirectory() as temp:
                    path = Path(temp) / 'effective-slab.svg'
                    write_discharge_svg(path, points, result, 'trial', particle_radii=radii)
                    root = ET.parse(path).getroot()
                spheres = [e for e in root.iter() if e.get('data-kind') == 'sphere-projection']
                self.assertEqual([int(e.get('data-particle-index')) for e in spheres], expected_indices)
                center_marks = [e for e in root.iter() if e.tag.endswith('circle') and e.get('fill-opacity') == '.76']
                self.assertEqual(len(center_marks), len(expected_indices))
                for sphere, center_mark, index in zip(spheres, center_marks, expected_indices):
                    self.assertEqual(float(sphere.get('data-radius-sim')), radii[index])
                    self.assertAlmostEqual(float(sphere.get('r')), radii[index] * float(sphere.get('data-px-per-unit')), delta=0.011)
                    self.assertEqual((sphere.get('cx'), sphere.get('cy')), (center_mark.get('cx'), center_mark.get('cy')))
                self.assertEqual(result, original)


    def test_rejected_three_cut_diagnostic_renders_actual_radii(self):
        points, radii = cloud((-2, 0, 2))
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / 'rejected.svg'
            write_discharge_failure_svg(path, points, 'trial', 'unsupported face',
                horizontal_axis='z', nbrratio=20,
                mesh_bounds={'x': (-5, 5), 'z': (-5, 5)}, particle_radii=radii)
            spheres = [e for e in ET.parse(path).getroot().iter() if e.get('data-kind') == 'sphere-projection']
        self.assertTrue(spheres)
        for sphere in spheres:
            self.assertAlmostEqual(float(sphere.get('r')), radii[int(sphere.get('data-particle-index'))] * float(sphere.get('data-px-per-unit')), delta=0.011)


    def test_controller_passes_snapshot_radii_to_success_renderer(self):
        from unittest.mock import patch
        from optimizer.calibration_controller import CalibrationController
        points, radii = cloud((-2, 0, 2))
        payload = ('# vtk DataFile Version 2.0\nTEST\nASCII\nDATASET POLYDATA\n'
                   f'POINTS {len(points)} double\n' + '\n'.join(' '.join(map(str, p)) for p in points)
                   + f'\nPOINT_DATA {len(points)}\nFIELD FieldData 1\nradius 1 {len(points)} double\n'
                   + '\n'.join(map(str, radii)) + '\n')
        with tempfile.TemporaryDirectory() as temp:
            snapshot = Path(temp) / 'trial_0001/Ang Deslizamento/post/deslizamento_final.vtk'
            snapshot.parent.mkdir(parents=True)
            snapshot.write_text(payload)
            controller = object.__new__(CalibrationController)
            with patch.object(controller, '_mesh_bounds_for_snapshot', return_value={'x': (-5, 5), 'z': (-5, 5)}), \
                 patch.object(controller, '_nbr_ratio_from_deck', return_value=20), \
                 patch.object(controller, '_nbr_note_for_snapshot', return_value='test'), \
                 patch('optimizer.calibration_controller.write_svg_png'):
                controller._measure_snapshot(snapshot, 'discharge')
            svg = snapshot.parents[2] / 'medicao_descarga.svg'
            spheres = [e for e in ET.parse(svg).getroot().iter() if e.get('data-kind') == 'sphere-projection']
        self.assertEqual(len(spheres), len(points))


if __name__ == '__main__':
    unittest.main()
