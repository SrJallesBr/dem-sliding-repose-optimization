import sys
import tempfile
import unittest
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT))

from optimizer.cohesion_guard import CohesionLimitTracker


class CohesionLimitTrackerTests(unittest.TestCase):
    def test_problematic_value_brackets_and_recommends_below_confirmed_failure(self):
        tracker = CohesionLimitTracker(minimum=0.0, maximum=10_000_000.0, max_probes=4)
        tracker.observe_valid(1000.0, trial_id=1, parameters={"sf11": 0.5, "rf11": 0.8})
        tracker.observe_problematic(5_770_000.0, trial_id=2, parameters={"sf11": 0.5, "rf11": 0.8}, reason="cohesive bridge")

        self.assertEqual(tracker.recommended_ceiling, 1000.0)
        self.assertEqual(tracker.next_probe(), 2_885_500.0)

    def test_binary_probes_move_valid_and_problematic_bracket_and_export_report(self):
        tracker = CohesionLimitTracker(minimum=0.0, maximum=10_000_000.0, max_probes=3)
        params = {"sf11": 0.5, "rf11": 0.8}
        tracker.observe_valid(1000.0, trial_id=1, parameters=params)
        tracker.observe_problematic(5_770_000.0, trial_id=2, parameters=params, reason="cohesive bridge")
        candidate = tracker.next_probe()
        tracker.observe_valid(candidate, trial_id=3, parameters=params)

        report = tracker.report()
        self.assertAlmostEqual(report["highest_measurable_CED11"], candidate)
        self.assertEqual(report["lowest_problematic_CED11"], 5_770_000.0)
        self.assertEqual(report["next_probe_CED11"], 4_327_750.0)
        self.assertEqual(report["recommended_ceiling_CED11"], candidate)

    def test_tracker_ignores_valid_values_above_a_known_problem_point_for_recommendation(self):
        tracker = CohesionLimitTracker(minimum=0.0, maximum=10.0, max_probes=4)
        tracker.observe_valid(9.0, trial_id=1, parameters={})
        tracker.observe_problematic(8.0, trial_id=2, parameters={}, reason="bad profile")
        self.assertEqual(tracker.recommended_ceiling, 0.0)
        self.assertIsNone(tracker.next_probe())


if __name__ == "__main__":
    unittest.main()
