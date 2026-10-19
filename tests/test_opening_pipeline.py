"""Software scheduler test only: analytic responses, NOT DEM results."""
import contextlib,io,tempfile,unittest
from collections import Counter
from pathlib import Path
from optimizer.opening_controller import OpeningController,OPENING_FIELDS
ROOT=Path(__file__).resolve().parents[1]

class SyntheticPipelineTests(unittest.TestCase):
    def test_doe_refinement_validation_and_one_seed_pair_per_trial(self):
        class AnalyticController(OpeningController):
            def execute_trial(self,trial_id,parameters,phase='doe',seeds=None):
                seeds=seeds or self._trial_seeds(trial_id)
                d=70+parameters['sf11']*2+parameters['CED11']/1e7
                r=38+parameters['rf11']+parameters['CED11']/2e7
                row={**parameters,'trial_id':trial_id,'phase':phase,'status':'complete','discharge_deg':d,'repose_deg':r,'discharge_error':abs(d-self.config['target_discharge_deg']),'repose_error':abs(r-self.config['target_repose_deg']),'distribution_seed':seeds['distribution'],'insertion_seed':seeds['insertion']}
                self._append_csv('results.csv',OPENING_FIELDS,row)
                return row
        with tempfile.TemporaryDirectory() as directory:
            c=AnalyticController(ROOT,Path(directory))
            with contextlib.redirect_stdout(io.StringIO()):rows=c.run_campaign()
            phases=Counter(r['phase'] for r in rows)
            for phase,field in [('doe','doe_trials'),('refinement','refinement_trials'),('validation','validation_trials')]:self.assertEqual(phases[phase],c.config['campaign'][field])
            self.assertEqual(len(rows),sum(c.config['campaign'][f] for f in ['doe_trials','refinement_trials','validation_trials']))
            self.assertEqual(len({r['trial_id'] for r in rows}),len(rows))
            self.assertTrue(all(r['score']!='' for r in rows))
            self.assertEqual(len({r['insertion_seed'] for r in rows}),len(rows))
if __name__=='__main__':unittest.main()
