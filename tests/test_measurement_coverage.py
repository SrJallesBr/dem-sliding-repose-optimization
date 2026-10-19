import unittest
from optimizer.measurement_coverage import summarize_coverage

class CoverageTests(unittest.TestCase):
    def test_missing_second_snapshot_is_not_silently_removed_from_trial_denominator(self):
        report=summarize_coverage([
            {'trial':'1','discharge':{'available':True,'valid':True},'repose':{'available':True,'valid':True}},
            {'trial':'2','discharge':{'available':True,'valid':True},'repose':{'available':False,'valid':False}},
            {'trial':'3','discharge':{'available':True,'valid':False},'repose':{'available':True,'valid':False}},
        ])
        self.assertEqual(report['completed_discharge_trials'],3)
        self.assertEqual(report['paired_valid_trials'],1)
        self.assertEqual(report['missing_repose_trials'],['2'])
        self.assertAlmostEqual(report['paired_coverage_all_completed_discharge_pct'],100/3)
        self.assertEqual(report['paired_snapshots_available'],2)
        self.assertEqual(report['paired_coverage_available_pairs_pct'],50)
        self.assertFalse(report['target_reached'])

    def test_target_requires_valid_pairs_not_individual_success_or_empty_input(self):
        self.assertFalse(summarize_coverage([])['target_reached'])
        rows=[{'trial':str(i),'discharge':{'available':True,'valid':True},'repose':{'available':True,'valid':i!=0}} for i in range(10)]
        report=summarize_coverage(rows)
        self.assertEqual(report['paired_coverage_all_completed_discharge_pct'],90)
        self.assertTrue(report['target_reached'])

if __name__=='__main__': unittest.main()
