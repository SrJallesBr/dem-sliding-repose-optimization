import json, unittest
from pathlib import Path
from optimizer.opening_measurements import measure_opening_section

class FullHeightOpeningTests(unittest.TestCase):
    def fixture(self):
        return json.loads((Path(__file__).parent/'fixtures/jenike_opening_baseline.json').read_text())

    def test_default_uses_entire_height_with_no_two_diameter_lower_cut(self):
        d=self.fixture(); r=measure_opening_section(d['points'],d['radii'],d['protocol'])
        self.assertEqual(r['minimum_height_particle_diameters'],0.)
        for side,sign in [('left',-1),('right',1)]:
            cloud=[p for p in r['section_points'] if sign*p[2]>0]
            bottom=min(p[1] for p in cloud)
            self.assertAlmostEqual(r['profile_start_height_'+side],bottom)
            self.assertAlmostEqual(min(y for h,y in r[side+'_profile_points']),bottom)
            self.assertLess(min(y for h,y in r[side+'_fit_points']),bottom+2*r['median_radius'])
            self.assertGreaterEqual(r[side+'_unique_particle_support'],5)

    def test_real_crest_caps_are_visible_but_not_used_as_slope_support(self):
        d=self.fixture();r=measure_opening_section(d['points'],d['radii'],d['protocol'])
        for side,crest in [('left',{895,769,643}),('right',{1879})]:
            ids=r[side+'_sample_particle_indices'];mask=r[side+'_fit']['inlier_mask']
            self.assertTrue(crest.issubset(set(ids)))
            self.assertTrue(crest.isdisjoint({pid for pid,keep in zip(ids,mask) if keep}))
            self.assertGreater(len(r[side+'_downweighted_points']),0)
            self.assertGreaterEqual(r[side+'_unique_particle_support'],5)
        # Supported high particles still define the face: not a blanket top cut.
        self.assertIn(1015,{i for i,k in zip(r['left_sample_particle_indices'],r['left_fit']['inlier_mask']) if k})
        self.assertIn(1750,{i for i,k in zip(r['right_sample_particle_indices'],r['right_fit']['inlier_mask']) if k})

    def test_active_configuration_has_no_discharge_lower_cut_and_gate_stays_closed(self):
        from tempfile import TemporaryDirectory
        from optimizer.opening_controller import OpeningController
        root=Path(__file__).resolve().parents[1]
        with TemporaryDirectory() as directory:
            c=OpeningController(root,Path(directory))
            self.assertEqual(c.config['measurement']['minimum_height_particle_diameters'],0.)
            self.assertIs(c.config['measurement_approved'],False)

if __name__=='__main__': unittest.main()
