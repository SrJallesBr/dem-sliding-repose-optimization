import math
import tempfile
import unittest
from pathlib import Path
from xml.etree import ElementTree as ET

from optimizer.measurements import measure_repose_eight_ray_least_squares
from optimizer.measurement_visuals import write_repose_svg


class ReposeSurfaceRegressionTests(unittest.TestCase):
    @staticmethod
    def cone_with_cutoff_gap():
        points, radii = [], []
        slope = math.tan(math.radians(30.0))
        for band in range(18):
            y = 0.009 + 0.015 * band
            radial = 0.50 - (y - 0.009) / slope
            for azimuth in range(72):
                theta = math.radians(azimuth * 5.0)
                points.append((radial * math.cos(theta), y, radial * math.sin(theta)))
                radii.append(0.01)
        # The true exterior particle below the center cutoff still intersects
        # the first sampled height. An internal particle center above cutoff
        # must not become a mandatory anchor that makes the real face outliers.
        points.append((0.0, 0.050, 0.15))
        radii.append(0.01)
        return points, radii

    def test_internal_first_partial_band_cannot_anchor_the_whole_face(self):
        points, radii = self.cone_with_cutoff_gap()
        result = measure_repose_eight_ray_least_squares(points, particle_radii=radii)
        pair = result['pair_profiles'][2]
        face = pair['positive_profile']
        self.assertAlmostEqual(face['angle_deg'], 30.0, delta=2.0)
        self.assertLessEqual(face['downweighted_count'], 2)
        self.assertGreater(face['r2'], 0.95)

    def test_true_surface_is_not_downweighted_to_accommodate_two_internal_start_samples(self):
        from optimizer.repose_profile import fit_repose_profile
        slope = math.tan(math.radians(30.0))
        samples = [(4.0-y/slope, y) for y in (0.6+i*0.125 for i in range(13))]
        samples[0] = (samples[0][0]-1.0, samples[0][1])
        samples[1] = (samples[1][0]-1.0, samples[1][1])
        fitted = fit_repose_profile(samples, 0.125)
        self.assertAlmostEqual(fitted['angle_deg'], 30.0, delta=0.5)
        self.assertEqual(fitted['downweighted_points'], samples[:2])

    def test_column_cannot_be_returned_as_a_valid_repose_angle(self):
        from optimizer.measurements import UnmeasurableAngle
        points, radii = [], []
        for band in range(28):
            for azimuth in range(72):
                theta = math.radians(azimuth*5)
                points.append((0.5*math.cos(theta), 0.01+band*0.015, 0.5*math.sin(theta)))
                radii.append(0.01)
        with self.assertRaises(UnmeasurableAngle) as caught:
            measure_repose_eight_ray_least_squares(points,particle_radii=radii)
        self.assertIsNone(caught.exception.diagnostic_result['mean_angle_deg'])
        self.assertEqual(len(caught.exception.diagnostic_result['pair_profiles']),4)

    def test_sustained_middle_departure_is_not_disguised_as_isolated_outliers(self):
        from optimizer.repose_profile import fit_repose_profile
        samples = [(4.0-y/math.tan(math.radians(30)),y) for y in (0.6+i*0.125 for i in range(16))]
        for i in (6,7,8):
            samples[i] = (samples[i][0]+1.0,samples[i][1])
        result = fit_repose_profile(samples,0.125)
        self.assertFalse(result['valid'])
        self.assertIn('sustained', result['rejection_reason'])

    def test_invalid_faces_keep_the_standard_layout_but_no_angles_or_fit_lines(self):
        from optimizer.measurements import UnmeasurableAngle
        points,radii=[],[]
        for band in range(28):
            for azimuth in range(72):
                theta=math.radians(azimuth*5)
                points.append((0.5*math.cos(theta),0.01+band*0.015,0.5*math.sin(theta)))
                radii.append(0.01)
        try:
            measure_repose_eight_ray_least_squares(points,particle_radii=radii)
        except UnmeasurableAngle as exc:
            result=exc.diagnostic_result
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'diagnostic.svg'
            write_repose_svg(path,result,'column','',points=points)
            root=ET.parse(path).getroot()
        self.assertIn('SEM ANGULO VALIDO',' '.join(root.itertext()))
        self.assertFalse(any(el.attrib.get('data-kind','').startswith('fit-') for el in root.iter()))
        self.assertEqual(sum(el.attrib.get('data-kind')=='toe-axis' for el in root.iter()),8)
        self.assertTrue(any(el.attrib.get('data-kind')=='selected-profile-trace' for el in root.iter()))

    def test_controller_uses_the_same_diagnostic_metadata_after_rejection(self):
        from unittest.mock import patch
        from optimizer.calibration_controller import CalibrationController, GeometryDiagnostic
        from optimizer.measurements import UnmeasurableAngle
        points,radii=self.cone_with_cutoff_gap()
        diagnostic=measure_repose_eight_ray_least_squares(points,particle_radii=radii)
        diagnostic['mean_angle_deg']=None
        diagnostic['stddev_deg']=None
        diagnostic['direction_angles_deg']=[None]*8
        for pair in diagnostic['pair_profiles']:
            pair['positive_angle_deg']=pair['negative_angle_deg']=pair['mean_angle_deg']=None
            for key in ('positive_profile','negative_profile'):
                pair[key]['valid']=False
                pair[key]['angle_deg']=None
                pair[key]['rejection_reason']='fixture invalid face'
        error=UnmeasurableAngle('fixture invalid face')
        error.diagnostic_result=diagnostic
        with tempfile.TemporaryDirectory() as directory:
            trial=Path(directory)/'trial_0001'
            snapshot=trial/'Ang Repouso/post/repouso_10s.vtk'
            snapshot.parent.mkdir(parents=True)
            snapshot.write_text('POINTS 1 float\n0 1 0\nPOINT_DATA 1\nFIELD FieldData 1\nradius 1 1 float\n0.01\n')
            controller=CalibrationController.__new__(CalibrationController)
            with patch('optimizer.calibration_controller.measure_pile_angle',side_effect=error), patch('optimizer.calibration_controller.write_svg_png'):
                with self.assertRaises(GeometryDiagnostic):
                    controller._measure_snapshot(snapshot,'repose')
            root=ET.parse(trial/'medicao_repouso_diagnostico.svg').getroot()
        self.assertEqual(sum(el.attrib.get('data-kind')=='toe-axis' for el in root.iter()),8)
        self.assertTrue(any(el.attrib.get('data-kind')=='selected-profile-trace' for el in root.iter()))

    def test_one_particle_cannot_supply_independent_face_support(self):
        from optimizer.repose_profile import fit_repose_profile
        samples=[(4-y/math.tan(math.radians(30)),y) for y in (0.6+i*0.1 for i in range(12))]
        result=fit_repose_profile(samples,0.125,particle_indices=[42]*len(samples))
        self.assertFalse(result['valid'])
        self.assertIn('unique particle', result['rejection_reason'])

    def test_low_quality_positive_r2_is_not_a_valid_face(self):
        from optimizer.repose_profile import fit_repose_profile
        samples=[(3-0.8*(i*0.125)+0.2*(-1)**i, i*0.125) for i in range(16)]
        result=fit_repose_profile(samples,0.125)
        self.assertGreater(result['r2'],0.0)
        self.assertLess(result['r2'],0.90)
        self.assertFalse(result['valid'])
        self.assertIn('orthogonal residual',result['rejection_reason'])

    def test_missing_height_range_cannot_be_bridged_as_a_supported_face(self):
        from optimizer.repose_profile import fit_repose_profile
        heights=[0.1*i for i in range(6)]+[2+0.1*i for i in range(6)]
        result=fit_repose_profile([(4-y/math.tan(math.radians(30)),y) for y in heights],0.125)
        self.assertFalse(result['valid'])
        self.assertIn('height gap',result['rejection_reason'])

    def test_sustained_endpoint_departure_is_not_trimmed_into_a_valid_face(self):
        from optimizer.repose_profile import fit_repose_profile
        samples=[(4-y/math.tan(math.radians(30)),y) for y in (0.6+i*0.125 for i in range(16))]
        for i in (0,1,2):
            samples[i]=(samples[i][0]-1.0,samples[i][1])
        result=fit_repose_profile(samples,0.125)
        self.assertFalse(result['valid'])
        self.assertIn('sustained',result['rejection_reason'])

    def test_repeated_samples_of_one_particle_do_not_change_weighted_r2(self):
        from optimizer.repose_profile import fit_repose_profile
        samples=[(4-1.7*(i*0.1)+(0.2 if i==5 else 0),i*0.1) for i in range(14)]
        base=fit_repose_profile(samples,0.125,particle_indices=list(range(14)))
        repeated=sorted(list(zip(samples,range(14)))+[(samples[5],5)]*30,key=lambda item:item[0][1])
        changed=fit_repose_profile([p for p,index in repeated],0.125,particle_indices=[index for p,index in repeated])
        self.assertAlmostEqual(base['angle_deg'],changed['angle_deg'],places=10)
        self.assertAlmostEqual(base['r2'],changed['r2'],places=10)
        self.assertEqual(base['valid'],changed['valid'])

    def test_sphere_crossing_center_contributes_to_both_faces_and_its_top_cap(self):
        from optimizer.repose_profile import sample_repose_section
        points=[(0.0,0.1,0.0)]
        radii=[0.06]
        for sign in (-1,1):
            for k in range(7):
                points.append((0.0,k*0.02,sign*0.03))
                radii.append(0.01)
        result=sample_repose_section(points,radii,list(range(len(points))),0.01,1.0,0.0,5)
        for side in ('left','right'):
            at_height=[e for e in result[side+'_sampling_evidence'] if abs(e['height']-0.1)<1e-9]
            self.assertEqual(at_height[0]['particle_index'],0)
            self.assertGreater(max(y for h,y in result[side+'_fit_points']),0.12)

    def test_refining_the_same_real_surface_does_not_reweight_particle_votes(self):
        import json
        from optimizer.repose_profile import fit_repose_profile
        fixture=json.loads((Path(__file__).parent/'fixtures/trial0007_surface_resolution.json').read_text())
        angles=[]
        for grid in fixture['grids']:
            result=fit_repose_profile(grid['points'],fixture['median_radius'],grid['particle_indices'])
            self.assertTrue(result['valid'])
            angles.append(result['angle_deg'])
        self.assertLess(max(angles)-min(angles),1.0)

    def test_reduced_failure_diagnostic_still_uses_the_actual_section_planes(self):
        from optimizer.measurement_visuals import write_repose_failure_diagnostic_svg
        points=[(0.0,1.0,0.0),(0.0,1.0,1.0)]
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'reduced.svg'
            write_repose_failure_diagnostic_svg(path,points,'reduced','missing support',particle_radii=[0.01,0.01])
            root=ET.parse(path).getroot()
        self.assertEqual(sum(el.attrib.get('data-kind')=='particle' for el in root.iter()),5)

    def test_production_and_explicit_historical_preview_have_distinct_provenance_only(self):
        import copy
        from unittest.mock import patch
        from optimizer.calibration_controller import CalibrationController

        points, radii = self.cone_with_cutoff_gap()
        expected = measure_repose_eight_ray_least_squares(points, particle_radii=radii)
        with tempfile.TemporaryDirectory() as directory:
            trial = Path(directory) / 'trial_0001'
            snapshot = trial / 'Ang Repouso/post/repouso_10s.vtk'
            snapshot.parent.mkdir(parents=True)
            snapshot.write_text(
                f'POINTS {len(points)} double\n'
                + '\n'.join(' '.join(map(str, point)) for point in points)
                + f'\nPOINT_DATA {len(points)}\nFIELD FieldData 1\nradius 1 {len(points)} double\n'
                + '\n'.join(map(str, radii)) + '\n')
            controller = CalibrationController.__new__(CalibrationController)
            with patch('optimizer.calibration_controller.write_svg_png'):
                result = controller._measure_snapshot(snapshot, 'repose')
            production = ET.parse(trial / 'medicao_repouso_oito_direcoes.svg').getroot()
            production_text = ' '.join(production.itertext())
            self.assertIn('MEDICAO DO TRIAL', production_text)
            self.assertNotIn('PREVIA DE DUMP HISTORICO', production_text)
            self.assertNotIn('REPROCESSAMENTO SOMENTE LEITURA', production_text)
            self.assertEqual(result, expected)
            original = copy.deepcopy(result)
            preview_path = trial / 'preview.svg'
            write_repose_svg(preview_path, result, trial.name,
                             controller._nbr_note_for_snapshot(snapshot),
                             points=points, historical_preview=True)
            preview = ET.parse(preview_path).getroot()
            self.assertIn('PREVIA DE DUMP HISTORICO; REPROCESSAMENTO SOMENTE LEITURA', ' '.join(preview.itertext()))
            self.assertNotIn('MEDICAO DO TRIAL', ' '.join(preview.itertext()))
            self.assertEqual(result, original)
            self.assertEqual(production.attrib, preview.attrib)
            self.assertEqual((production.get('width'), production.get('height')), ('1800', '1160'))
            # The context changes only the footer text, not numbers or layout.
            production_footer = next(e for e in production if e.get('y') == '1122')
            preview_footer = next(e for e in preview if e.get('y') == '1122')
            self.assertEqual(production_footer.attrib, preview_footer.attrib)
            production_footer.text = preview_footer.text = ''
            self.assertEqual(ET.tostring(production), ET.tostring(preview))

    def test_gray_cloud_contains_only_particles_intersecting_each_plane(self):
        points, radii = self.cone_with_cutoff_gap()
        result = measure_repose_eight_ray_least_squares(points, particle_radii=radii)
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / 'repose.svg'
            write_repose_svg(output, result, 'regression', '', points=points)
            root = ET.parse(output).getroot()
        actual = sum(el.attrib.get('data-kind') in {'particle', 'basal-reference'} for el in root.iter())
        expected = sum(pair['section_particle_count'] for pair in result['pair_profiles'])
        self.assertEqual(actual, expected)


if __name__ == '__main__':
    unittest.main()
