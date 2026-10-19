import json,math,unittest
from pathlib import Path
from optimizer import opening_measurements

class TrueSurfaceTests(unittest.TestCase):
    def test_real_baseline_inner_contour_and_five_unique_grains_not_coarse_bins(self):
        d=json.loads((Path(__file__).parent/'fixtures/jenike_opening_baseline.json').read_text())
        r=opening_measurements.measure_opening_section(d['points'],d['radii'],d['protocol'],minimum_height_particle_diameters=2.0)
        self.assertEqual(r['section_count'],1)
        self.assertEqual(r['slice_method'],'central_slab_projection')
        self.assertAlmostEqual(r['section_half_width_physical'],.075)
        self.assertAlmostEqual(r['left_angle_deg'],d['expected_left'],places=8)
        self.assertAlmostEqual(r['right_angle_deg'],d['expected_right'],places=8)
        for side in ('left','right'):
            self.assertGreaterEqual(r[side+'_unique_particle_support'],5)
            for e in r[side+'_sampling_evidence']:
                h,y=e['surface'];ch,cy=e['center']
                self.assertAlmostEqual((h-ch)**2+(y-cy)**2,e['radius']**2,places=9)
                self.assertGreaterEqual(cy,r['profile_start_height_'+side])

    def test_isolated_axis_crossings_do_not_become_a_sustained_bridge(self):
        self.assertFalse(opening_measurements.has_sustained_bridge([0.,1.],.03125,.125))
        self.assertTrue(opening_measurements.has_sustained_bridge([i*.03125 for i in range(10)],.03125,.125))

    def test_real_faces_survive_two_disconnected_axis_obstructions(self):
        d=json.loads((Path(__file__).parent/'fixtures/jenike_opening_baseline.json').read_text())
        radius=.125
        cutoff=min(p[1] for p in d['points'] if p[2]>0 and p[1]>=0)+4*radius
        for k in (4,16):
            d['points'].append([0.,cutoff+k*.25*radius,.001]);d['radii'].append(.009)
        r=opening_measurements.measure_opening_section(d['points'],d['radii'],d['protocol'],minimum_height_particle_diameters=2.0)
        self.assertTrue(r['valid'])
        self.assertGreaterEqual(r['right_unique_particle_support'],5)
        self.assertEqual(len(r['right_blocked_heights']),2)

    def test_oversized_sampling_grid_is_rejected_before_consensus_cost(self):
        from unittest.mock import patch
        from optimizer.measurements import UnmeasurableAngle
        d=json.loads((Path(__file__).parent/'fixtures/jenike_opening_baseline.json').read_text())
        with patch('optimizer.repose_profile.fit_repose_profile',side_effect=AssertionError('consensus must not run')):
            with self.assertRaisesRegex(UnmeasurableAngle,'512'):
                opening_measurements.measure_opening_section(d['points'],d['radii'],d['protocol'],height_band_factor=.01)

if __name__=='__main__':unittest.main()
