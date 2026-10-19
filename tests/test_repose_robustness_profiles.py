"""Regression evidence for candidate REPOSE policy; no solver or source-run writes."""
import json
import math
import unittest
from pathlib import Path
from optimizer.repose_profile import fit_repose_profile

FIXTURE = Path(__file__).parent/'fixtures/test_repose_robustness_saved_profiles.json'

class ReposeRobustnessTests(unittest.TestCase):
    def test_three_particles_over_one_radius_are_not_sustained_geometry(self):
        trial=json.loads(FIXTURE.read_text())['trials']['trial_0026']
        face=trial['faces'][0]
        fitted=fit_repose_profile(face['points'],trial['median_radius'],face['particle_indices'])
        self.assertTrue(fitted['valid'], fitted['rejection_reason'])
        self.assertAlmostEqual(fitted['angle_deg'],36.18901934883905,places=9)
        run=next(run for run in fitted['departure_runs'] if run['distinct_particles']==3)
        self.assertAlmostEqual(run['height_span_radii'],1.0)
        self.assertFalse(run['sustained'])

    def test_real_trial30_low_vertical_r2_has_good_orthogonal_quality(self):
        trial=json.loads(FIXTURE.read_text())['trials']['trial_0030']
        face=next(f for f in trial['faces'] if f['direction_deg']==45)
        fitted=fit_repose_profile(face['points'],trial['median_radius'],face['particle_indices'])
        self.assertTrue(fitted['valid'], fitted['rejection_reason'])
        self.assertLess(fitted['r2'],0.90)  # preserve the old axis-dependent diagnostic
        self.assertLess(fitted['orthogonal_rms_radii'],1.0)
        self.assertGreater(fitted['orthogonal_quality'],0.90)
        self.assertGreater(fitted['fitted_horizontal_span_radii'],2.0)
        self.assertEqual(fitted['validity_policy']['r2_y_on_h_role'],'diagnostic_only')

    def test_short_large_coherent_step_remains_invalid(self):
        # Three independently supported samples, two radii in height, but an
        # eight-radius horizontal step: not grain-scale roughness (old guard).
        points=[(4-y/math.tan(math.radians(30)),y) for y in (0.6+i*0.125 for i in range(16))]
        for i in (6,7,8): points[i]=(points[i][0]+1.0,points[i][1])
        fitted=fit_repose_profile(points,0.125)
        self.assertFalse(fitted['valid'])
        self.assertIn('sustained',fitted['rejection_reason'])
        self.assertTrue(any(r['large_coherent_step'] for r in fitted['departure_runs']))

    def test_adjacent_outlier_cannot_mask_large_coherent_step(self):
        # A moderate shoulder or opposite-sign neighbor must not hide the
        # independently supported, four-radius departure at indices 6-8.
        for neighbor in (5,9):
            for offset in (0.5,-1.0):
                with self.subTest(neighbor=neighbor,offset=offset):
                    points=[(4-y/math.tan(math.radians(30)),y) for y in (0.6+i*0.125 for i in range(16))]
                    for i in (6,7,8): points[i]=(points[i][0]+1.0,points[i][1])
                    points[neighbor]=(points[neighbor][0]+offset,points[neighbor][1])
                    fitted=fit_repose_profile(points,0.125)
                    self.assertFalse(fitted['valid'])
                    self.assertIn('sustained',fitted['rejection_reason'])
                    self.assertTrue(any(r['large_coherent_step'] for r in fitted['departure_runs']))

    def test_saved_piles_recover_without_accepting_whole_column_trials(self):
        trials=json.loads(FIXTURE.read_text())['trials']
        for name,trial in trials.items():
            with self.subTest(trial=name):
                valid=[]
                for face in trial['faces']:
                    fitted=fit_repose_profile(face['points'],trial['median_radius'],face['particle_indices'])
                    sign=1 if face['direction_deg'] in (0,45,90,135) else -1
                    valid.append(fitted['valid'] and sign*fitted['horizontal_gradient']<0)
                self.assertEqual(all(valid),trial['expected_all_faces_measurable'])

    def test_vertical_and_unresolved_horizontal_drift_remain_invalid(self):
        for gradient in (0.0,0.001):
            with self.subTest(gradient=gradient):
                fitted=fit_repose_profile([(2-gradient*y,y) for y in (i*.1 for i in range(25))],.125)
                self.assertFalse(fitted['valid'])
                self.assertIn('near-vertical',fitted['rejection_reason'])

    def test_resolved_near_vertical_inclined_face_is_not_arbitrarily_capped(self):
        slope=math.tan(math.radians(80))
        fitted=fit_repose_profile([(2-y/slope,y) for y in (i*.025 for i in range(100))],.01)
        self.assertTrue(fitted['valid'],fitted['rejection_reason'])
        self.assertAlmostEqual(fitted['angle_deg'],80.0,places=9)

    def test_long_reversal_and_large_gap_remain_invalid(self):
        points=[(5-y,y) for y in (i*.1 for i in range(40))]
        for i in range(12,24):points[i]=(points[i][0]+1,points[i][1])
        fitted=fit_repose_profile(points,.125)
        self.assertFalse(fitted['valid'])
        self.assertIn('sustained',fitted['rejection_reason'])
        heights=[i*.1 for i in range(6)]+[2+i*.1 for i in range(6)]
        fitted=fit_repose_profile([(5-y,y) for y in heights],.125)
        self.assertFalse(fitted['valid'])
        self.assertIn('height gap',fitted['rejection_reason'])

    def test_physical_support_cannot_be_created_by_dense_height_sampling(self):
        fitted=fit_repose_profile([(3-i*.01,i*.01) for i in range(20)],.125)
        self.assertFalse(fitted['valid'])
        self.assertIn('physical height',fitted['rejection_reason'])

    def test_legacy_r2_failure_still_rejects_bad_geometry_for_geometric_reason(self):
        points=[(3-.8*(i*.125)+.2*(-1)**i,i*.125) for i in range(16)]
        fitted=fit_repose_profile(points,.125)
        self.assertFalse(fitted['valid'])
        self.assertIn('orthogonal residual',fitted['rejection_reason'])
        self.assertGreater(fitted['orthogonal_rms_radii'],1)

    def test_repeated_heights_preserve_quality_and_zero_weight_semantics(self):
        trial=json.loads(FIXTURE.read_text())['trials']['trial_0026'];face=trial['faces'][0]
        a=fit_repose_profile(face['points'],trial['median_radius'],face['particle_indices'])
        pairs=list(zip(face['points'],face['particle_indices']))
        pairs=sorted(pairs+[pairs[30]]*20,key=lambda p:p[0][1])
        b=fit_repose_profile([p for p,i in pairs],trial['median_radius'],[i for p,i in pairs])
        self.assertEqual(a['valid'],b['valid'])
        self.assertAlmostEqual(a['angle_deg'],b['angle_deg'],places=9)
        self.assertAlmostEqual(a['orthogonal_quality'],b['orthogonal_quality'],places=9)
        self.assertAlmostEqual(a['orthogonal_rms_radii'],b['orthogonal_rms_radii'],places=9)
        self.assertTrue(all(w==0 for w,k in zip(b['fit_sample_weights'],b['inlier_mask']) if not k))
