import math,unittest
from optimizer import opening_measurements

class SectionTests(unittest.TestCase):
    def test_inner_face_radius_correction_and_single_plane(self):
        pts=[];radii=[];m=math.tan(math.radians(30))
        for side in (-1,1):
            for i in range(24):
                y=.2+i*.26;r=.125 if i%2 else .1625;h=side*(.5+y/m+r)
                pts.append((0,y,h));radii.append(r)
        protocol={'nbrratio':31.25,'plane_x_sim':0.0,'final_base_bounds_sim':{'x':[-100,100],'z':[-100,100]},'snapshot_times_s':{'repose':13.,'discharge':23.}}
        result=opening_measurements.measure_opening_section_legacy_height(pts,radii,protocol)
        self.assertEqual(result['slice_method'],'sphere_plane_intersection')
        self.assertEqual(result['section_count'],1)
        self.assertAlmostEqual(result['mean_angle_deg'],30.,places=8)
        for side in ('left','right'):
            for h,y in result[side+'_fit_points']:
                self.assertAlmostEqual(abs(h),.5+y/m,places=8)
        self.assertEqual(result['response_definition'],'post_opening_inner_free_surface')
    def test_sampling_scale_uses_sphere_radius_not_small_intersection_caps(self):
        pts=[];radii=[]
        for sign in (-1,1):
            for i in range(25):
                y=.2+i*.26;pts.append((.075,y,sign*(1+y*2+.1)));radii.append(.125)
        protocol={'nbrratio':31.25,'plane_x_sim':0.,'final_base_bounds_sim':{'x':[-100,100],'z':[-100,100]},'snapshot_times_s':{'repose':13.,'discharge':23.}}
        r=opening_measurements.measure_opening_section_legacy_height(pts,radii,protocol)
        self.assertAlmostEqual(r['height_bin_size'],.25)
        self.assertAlmostEqual(r['section_radii'][0],.1)
        self.assertEqual(r['section_original_particle_indices'],list(range(len(pts))))

    def test_particles_below_support_surface_are_not_a_false_basal_anchor(self):
        pts=[];radii=[]
        for sign in (-1,1):
            for i in range(24):
                y=.2+i*.26;pts.append((0,y,sign*(1+y*2+.125)));radii.append(.125)
            for i in range(5):pts.append((0,-.5+i*.08,sign*.03));radii.append(.125)
        protocol={'nbrratio':31.25,'plane_x_sim':0.,'final_base_bounds_sim':{'x':[-100,100],'z':[-100,100]},'snapshot_times_s':{'repose':13.,'discharge':23.}}
        r=opening_measurements.measure_opening_section_legacy_height(pts,radii,protocol)
        self.assertTrue(all(p[1]>=0 for p in r['section_points']))
        self.assertEqual(r['excluded_below_support_centers'],10)

if __name__=='__main__':unittest.main()
