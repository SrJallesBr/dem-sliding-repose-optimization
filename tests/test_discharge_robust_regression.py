"""Discharge regressions distilled from read-only saved-dump profile evidence."""
import json
import unittest
from pathlib import Path
from optimizer.measurements import measure_v_sliding_angle_by_height

FIXTURES=json.loads((Path(__file__).parent/'fixtures/discharge_robust_profiles.json').read_text())

def replay_selected_profile(tid, scale=1.):
    f=FIXTURES[str(tid)]; r=f['median_radius']; points=[]; radii=[]
    for sign in [-1,1]:
        points.append((0.,f['base_height']*scale,sign*9.*scale));radii.append(r*scale)
        for _,z,y,radius in f['surface_samples']:
            points.append((0.,y*scale,sign*(abs(z)-radius)*scale));radii.append(radius*scale)
    return measure_v_sliding_angle_by_height(points,radii,horizontal_axis='z',minimum_height_particle_diameters=2.)

class DischargeRobustRegressions(unittest.TestCase):
    def test_subradius_alternating_deltas_are_not_a_sustained_reversal(self):
        # Trial13 has two tiny inward edges separated by a much larger outward
        # edge. Raw 2-of-3 signs falsely choose a 3/13 split and reject it.
        result=replay_selected_profile(13)
        self.assertEqual(result['left_segment_count'],1)
        self.assertEqual(result['right_segment_count'],1)

    def test_shallow_but_sustained_particle_scale_return_keeps_two_segments(self):
        # Each 3-band net displacement is below one diameter, but the complete
        # supported lobe spans 4.8 radii. A local threshold must not erase it.
        points=[]; radii=[]
        for sign in [-1,1]:
            for i in range(17):
                outward=.5+.009375*min(i,16-i)
                # Binary-exact height spacing isolates detector behavior from
                # the legacy selector's floating-point bin-boundary issue.
                points.append((0.,.0625+.03125*i,sign*(outward-.015625)))
                radii.append(.015625)
        result=measure_v_sliding_angle_by_height(points,radii,horizontal_axis='z')
        self.assertEqual(result['left_segment_count'],2)
        self.assertEqual(result['right_segment_count'],2)
        self.assertTrue(all(segment['fit_support']>=5 for segment in result['left_segments']))

    def test_true_return_on_a_short_plateau_is_not_hidden_by_zero_deltas(self):
        # A flat turning region is still a genuine free-face change. The lower
        # piece owns the plateau; both flanks have ample independent support.
        points=[]; radii=[]
        outward=[.5+.05*i for i in range(7)]+[.8,.8]+[.8-.05*i for i in range(1,8)]
        for sign in [-1,1]:
            for i,value in enumerate(outward):
                points.append((0.,.0625+.03125*i,sign*(value-.015625)))
                radii.append(.015625)
        result=measure_v_sliding_angle_by_height(points,radii,horizontal_axis='z')
        self.assertEqual(result['left_segment_count'],2)
        self.assertEqual(result['right_segment_count'],2)

    def test_saved_small_excursions_recover_without_relaxing_support(self):
        for tid in [2,11,13]:
            with self.subTest(trial=tid):
                result=replay_selected_profile(tid)
                self.assertEqual(result['left_segment_count'],1)
                self.assertEqual(result['right_segment_count'],1)
                self.assertGreaterEqual(result['left_candidate_count'],5)

    def test_saved_true_unsupported_return_still_rejects(self):
        with self.assertRaisesRegex(ValueError,'insufficient segment support'):
            replay_selected_profile(10)

    def test_saved_trial7_face_retains_exact_baseline_fit(self):
        result=replay_selected_profile(7)
        self.assertAlmostEqual(result['left_angle_deg'],43.606752535808745,places=11)
        self.assertEqual(result['left_candidate_count'],19)
        self.assertEqual(result['left_segment_count'],1)
        self.assertTrue(all(y>=result['profile_start_height_left'] for _,y in result['left_fit_points']))

    def test_physical_noise_prominence_is_scale_invariant_with_individual_radii(self):
        reference=replay_selected_profile(2)
        for scale in [.01,10,100]:
            with self.subTest(scale=scale):
                result=replay_selected_profile(2,scale)
                self.assertEqual(result['left_segment_count'],reference['left_segment_count'])
                self.assertAlmostEqual(result['left_angle_deg'],reference['left_angle_deg'],places=9)

    def test_trial29_supported_turn_vertex_can_anchor_the_upper_segment(self):
        # Same physical extremum, two disjoint assignments of its vertex:
        # lower owns it -> 20/4 (invalid); upper owns it -> 19/5 (supported).
        # Do not shift the turn to an unrelated weaker lobe or fall back.
        result=replay_selected_profile(29)
        self.assertEqual(result['left_segment_count'],2)
        self.assertEqual([s['fit_support'] for s in result['left_segments']],[19,5])

if __name__=='__main__': unittest.main()
