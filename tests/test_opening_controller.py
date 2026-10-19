import csv,json,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
from optimizer.calibration_controller import GeometryDiagnostic
from optimizer import opening_controller
ROOT=Path(__file__).resolve().parents[1]

class CombinedControllerTests(unittest.TestCase):
    def test_one_solver_produces_two_ordered_snapshots_and_no_three_cut_csv(self):
        with tempfile.TemporaryDirectory() as directory:
            c=opening_controller.OpeningController(ROOT,Path(directory))
            self.assertEqual(c.preflight()['protocol'],'repose_then_base_opening')
            events=[]
            repose={'mean_angle_deg':40.,'pair_profiles':[{'axis_angle_deg':a,'mean_angle_deg':40.,'positive_angle_deg':40.,'negative_angle_deg':40.} for a in [0,45,90,135]]}
            discharge={'mean_angle_deg':70.,'left_angle_deg':69.,'right_angle_deg':71.,'asymmetry_deg':2.,'left_segments':[],'right_segments':[],'section_count':1,'snapshot_time_s':23.}
            def run(*args):events.append(('solver',args[1:]))
            def measure(path,kind):
                events.append((kind,path.name));return repose if kind=='repose' else discharge
            with patch.object(c,'_run_input',side_effect=run) as solver,patch.object(c,'_measure_snapshot',side_effect=measure):
                r=c.execute_trial(1,{'sf11':.5,'rf11':.8,'CED11':1000.})
            self.assertEqual(solver.call_count,1)
            self.assertEqual([e[0] for e in events],['solver','repose','discharge'])
            self.assertEqual(events[1][1],'repouso_final.vtk');self.assertEqual(events[2][1],'descarga_final.vtk')
            self.assertEqual(r['status'],'complete');self.assertEqual(r['discharge_sections'],1)
            with (Path(directory)/'results.csv').open() as f:columns=next(csv.reader(f))
            self.assertFalse(any('xplus' in n or 'xminus' in n for n in columns))
            case=Path(directory)/'trial_0001'/'Ang Repouso e Desgarga'
            self.assertTrue((case/'in.repousodesgarga').is_file())
            self.assertFalse((case/'post').exists())
            self.assertFalse(list(case.glob('restart*')))
    def test_geometry_failure_preserves_other_response_and_does_not_rerun(self):
        with tempfile.TemporaryDirectory() as directory:
            c=opening_controller.OpeningController(ROOT,Path(directory))
            def measure(path,kind):
                if kind=='repose':raise GeometryDiagnostic('synthetic unsupported repose')
                return {'mean_angle_deg':70.}
            with patch.object(c,'_run_input') as solver,patch.object(c,'_measure_snapshot',side_effect=measure) as reader:
                with self.assertRaises(GeometryDiagnostic):c.execute_trial(1,{'sf11':.5,'rf11':.8,'CED11':1000.})
            self.assertEqual(solver.call_count,1);self.assertEqual(reader.call_count,2)
            data=json.loads((Path(directory)/'trial_0001/partial_responses.json').read_text())
            self.assertEqual(data['angles_deg'],{'discharge':70.})
            self.assertFalse((Path(directory)/'results.csv').exists())

    def test_nonzero_process_exit_is_solver_failure_not_geometry(self):
        import sys
        with tempfile.TemporaryDirectory() as directory:
            c=opening_controller.OpeningController(ROOT,Path(directory),solver_command=[sys.executable,'-c','import sys;sys.exit(7)'])
            with patch.object(c,'_measure_snapshot') as reader:
                with self.assertRaisesRegex(RuntimeError,'exit code 7'):c.execute_trial(1,{'sf11':.5,'rf11':.8,'CED11':1000.})
            reader.assert_not_called()
            with (Path(directory)/'status.csv').open() as f:rows=list(csv.DictReader(f))
            self.assertEqual(rows[-1]['status'],'failed')
            self.assertFalse((Path(directory)/'results.csv').exists())

    def test_artifact_failures_attempt_both_responses_and_preserve_failed_status(self):
        for failed_kind in ('repose','discharge'):
            with self.subTest(failed_kind=failed_kind),tempfile.TemporaryDirectory() as directory:
                c=opening_controller.OpeningController(ROOT,Path(directory))
                def measure(path,kind):
                    if kind==failed_kind:raise RuntimeError('synthetic corrupt '+kind)
                    return {'mean_angle_deg':40. if kind=='repose' else 73.4}
                with patch.object(c,'_run_input') as solver,patch.object(c,'_measure_snapshot',side_effect=measure) as reader:
                    with self.assertRaisesRegex(RuntimeError,'synthetic corrupt') as caught:
                        c.execute_trial(1,{'sf11':.5,'rf11':.8,'CED11':1000.})
                self.assertNotIsInstance(caught.exception,GeometryDiagnostic)
                self.assertEqual(solver.call_count,1)
                self.assertEqual([call.args[1] for call in reader.call_args_list],['repose','discharge'])
                data=json.loads((Path(directory)/'trial_0001/partial_responses.json').read_text())
                other='discharge' if failed_kind=='repose' else 'repose'
                self.assertEqual(data['angles_deg'],{other:73.4 if other=='discharge' else 40.})
                self.assertEqual(data['failures'][failed_kind]['type'],'RuntimeError')
                self.assertEqual(data['failures'][failed_kind]['category'],'artifact_or_postprocess')
                self.assertIn('synthetic corrupt',data['failures'][failed_kind]['message'])
                with (Path(directory)/'status.csv').open() as f:self.assertEqual(list(csv.DictReader(f))[-1]['status'],'failed')
                self.assertFalse((Path(directory)/'results.csv').exists())

    def test_mixed_failures_are_artifact_failure_not_geometry_only(self):
        for errors in [(GeometryDiagnostic('geometry'),RuntimeError('artifact')),(RuntimeError('artifact'),GeometryDiagnostic('geometry'))]:
            with self.subTest(errors=errors),tempfile.TemporaryDirectory() as directory:
                c=opening_controller.OpeningController(ROOT,Path(directory))
                with patch.object(c,'_run_input') as solver,patch.object(c,'_measure_snapshot',side_effect=errors) as reader:
                    with self.assertRaises(RuntimeError) as caught:c.execute_trial(1,{'sf11':.5,'rf11':.8,'CED11':1000.})
                self.assertNotIsInstance(caught.exception,GeometryDiagnostic)
                self.assertEqual(solver.call_count,1);self.assertEqual(reader.call_count,2)
                data=json.loads((Path(directory)/'trial_0001/partial_responses.json').read_text())
                self.assertEqual(data['angles_deg'],{})
                self.assertEqual({failure['category'] for failure in data['failures'].values()},{'geometry','artifact_or_postprocess'})
                with (Path(directory)/'status.csv').open() as f:self.assertEqual(list(csv.DictReader(f))[-1]['status'],'failed')

    def test_missing_and_invalid_real_artifacts_do_not_skip_other_snapshot(self):
        for raw in (None,b'',b'corrupt VTK'):
            with self.subTest(raw=raw),tempfile.TemporaryDirectory() as directory:
                c=opening_controller.OpeningController(ROOT,Path(directory));measure_snapshot=c._measure_snapshot
                def run(workspace,*args):
                    path=workspace/opening_controller.CASE_DIRECTORY/opening_controller.SNAPSHOTS['repose']
                    path.parent.mkdir(parents=True)
                    if raw is not None:path.write_bytes(raw)
                def measure(path,kind):
                    return measure_snapshot(path,kind) if kind=='repose' else {'mean_angle_deg':73.4}
                with patch.object(c,'_run_input',side_effect=run) as solver,patch.object(c,'_measure_snapshot',side_effect=measure) as reader:
                    with self.assertRaises(RuntimeError):c.execute_trial(1,{'sf11':.5,'rf11':.8,'CED11':1000.})
                self.assertEqual(solver.call_count,1);self.assertEqual(reader.call_count,2)
                data=json.loads((Path(directory)/'trial_0001/partial_responses.json').read_text())
                self.assertEqual(data['angles_deg'],{'discharge':73.4})
                self.assertEqual(data['failures']['repose']['type'],'FileNotFoundError' if raw is None else 'RuntimeError')

    def test_render_failure_preserves_computed_angle_and_attempts_other_snapshot(self):
        with tempfile.TemporaryDirectory() as directory:
            c=opening_controller.OpeningController(ROOT,Path(directory));measure_snapshot=c._measure_snapshot
            def run(workspace,*args):
                path=workspace/opening_controller.CASE_DIRECTORY/opening_controller.SNAPSHOTS['repose']
                path.parent.mkdir(parents=True)
                path.write_text('# vtk DataFile Version 3.0\nsynthetic\nASCII\nDATASET POLYDATA\nPOINTS 1 float\n0 1 0\nPOINT_DATA 1\nFIELD FieldData 1\nradius 1 1 float\n0.1\n')
            def measure(path,kind):
                return measure_snapshot(path,kind) if kind=='repose' else {'mean_angle_deg':73.4}
            with patch.object(c,'_run_input',side_effect=run) as solver,patch.object(c,'_measure_snapshot',side_effect=measure) as reader,patch.object(opening_controller,'measure_repose_eight_ray_least_squares',return_value={'mean_angle_deg':40.}),patch.object(opening_controller,'write_repose_svg',side_effect=RuntimeError('render failed')):
                with self.assertRaisesRegex(RuntimeError,'render failed'):c.execute_trial(1,{'sf11':.5,'rf11':.8,'CED11':1000.})
            self.assertEqual(solver.call_count,1);self.assertEqual(reader.call_count,2)
            data=json.loads((Path(directory)/'trial_0001/partial_responses.json').read_text())
            self.assertEqual(data['angles_deg'],{'repose':40.,'discharge':73.4})
            self.assertEqual(data['failures']['repose']['type'],'RuntimeError')
            self.assertFalse((Path(directory)/'results.csv').exists())

    def test_preflight_rejects_unsupported_face_and_invalid_measurement_scales(self):
        invalid=[('face','outer'),('section_half_width_physical_m',0.),('section_half_width_physical_m',-1.),
            ('section_half_width_physical_m',float('nan')),('section_half_width_physical_m',float('inf')),
            ('height_band_factor',0.),('height_band_factor',-1.),('height_band_factor',float('nan')),('height_band_factor',float('inf')),
            ('minimum_height_particle_diameters',-1.),('minimum_height_particle_diameters',float('nan')),('minimum_height_particle_diameters',float('inf'))]
        for key,value in invalid:
            with self.subTest(key=key,value=value),tempfile.TemporaryDirectory() as directory:
                c=opening_controller.OpeningController(ROOT,Path(directory));c.config['measurement'][key]=value
                with self.assertRaisesRegex(RuntimeError,'preflight failed'):c.preflight()

    def test_preflight_accepts_zero_height_cutoff(self):
        with tempfile.TemporaryDirectory() as directory:
            c=opening_controller.OpeningController(ROOT,Path(directory));c.config['measurement']['minimum_height_particle_diameters']=0.
            self.assertEqual(c.preflight()['status'],'ok')

    def test_live_phase_uses_same_combined_log(self):
        with tempfile.TemporaryDirectory() as directory:
            c=opening_controller.OpeningController(ROOT,Path(directory));trial=c.prepare_trial(1,{'sf11':.5,'rf11':.8,'CED11':1000.})
            (trial/'sequential.log').write_text('Step Time Atoms KinEng\n150000 15 5 1\n')
            state={1:{'status':'preparing'}};c.refresh_live_state(Path(directory),state)
            self.assertEqual(state[1]['status'],'opening');self.assertEqual(state[1]['simulation_time_s'],15.)

    def test_unapproved_measurement_preview_blocks_campaign_before_solver(self):
        with tempfile.TemporaryDirectory() as directory, patch.object(opening_controller.OpeningController,'run_campaign') as campaign:
            with self.assertRaisesRegex(RuntimeError,'Approve the single-cut Jenike measurement preview'):
                opening_controller.main(['--project-root',str(ROOT),'--run-root',str(Path(directory)/'campaign'),'--run'])
            campaign.assert_not_called()

if __name__=='__main__':unittest.main()
