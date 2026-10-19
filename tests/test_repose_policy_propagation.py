import unittest
from optimizer.measurements import measure_repose_eight_ray_least_squares
import test_repose_surface_regression as fixtures

class ReposePolicyPropagationTests(unittest.TestCase):
    def test_public_measurement_records_actual_geometric_acceptance_policy(self):
        points,radii=fixtures.ReposeSurfaceRegressionTests.cone_with_cutoff_gap()
        result=measure_repose_eight_ray_least_squares(points,particle_radii=radii)
        policy=result['validity_policy']
        self.assertNotIn('minimum_inlier_r2',policy)
        self.assertEqual(policy['version'],'repose-geometric-quality-v1')
        self.assertEqual(policy['minimum_fit_samples'],8)
        for face in result['direction_profiles']:
            self.assertEqual(face['validity_policy']['version'],policy['version'])
            self.assertEqual(face['r2_role'],'diagnostic_y_on_horizontal_not_acceptance_gate')
            self.assertGreater(face['orthogonal_quality'],0.9)
            self.assertLess(face['orthogonal_rms_radii'],1.)
            self.assertEqual(len(face['fit_sample_weights']),len(face['profile_points']))
            self.assertFalse(face['fit_quality_warning'])

if __name__=='__main__':unittest.main()
