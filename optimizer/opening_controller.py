"""Optimization for one sequential repose -> split-base opening LIGGGHTS case."""
from __future__ import annotations
import argparse,csv,hashlib,json,math,os,re,shutil,sys
from datetime import datetime
from pathlib import Path
from optimizer.calibration_controller import (CalibrationController,GeometryDiagnostic,PARAMETER_LINE,SEED_PATTERNS,STATUS_FIELDS,_utc_now,rewrite_input_parameters,rewrite_input_seeds,repose_pair_result_fields,simulation_time_seconds)
from optimizer.article_optimizer import weighted_normalized_score
from optimizer.opening_protocol import CASE_DIRECTORY,INPUT_NAME,SNAPSHOTS,inspect_case
from optimizer.opening_measurements import measure_opening_section
from optimizer.measurements import UnmeasurableAngle,parse_vtk_points,parse_vtk_scalar_field,measure_repose_eight_ray_least_squares
from optimizer.measurement_visuals import write_repose_svg,write_repose_failure_diagnostic_svg,write_svg_png

OPENING_FIELDS=('trial_id','phase','status','sf11','rf11','CED11','distribution_seed','insertion_seed','repose_deg','discharge_deg','repose_time_s','discharge_time_s','discharge_left_deg','discharge_right_deg','discharge_asymmetry_deg','discharge_sections','discharge_plane','discharge_plane_x_m','discharge_slice_method','discharge_band_half_width_m','discharge_profile_method',*(f'repose_axis{a}{suffix}_deg' for a in [0,45,90,135] for suffix in ['','_plus','_minus']),*(f'discharge_{side}_segment{n}_{f}' for side in ['left','right'] for n in [1,2] for f in ['deg','support']),'discharge_error','repose_error','score','workspace')

class OpeningController(CalibrationController):
    def preflight(self):
        case=self.project_root/CASE_DIRECTORY
        protocol=inspect_case(case);text=(case/INPUT_NAME).read_text(encoding='utf-8');problems=[]
        if self.config.get('protocol')!='repose_then_base_opening':problems.append('wrong optimization protocol')
        measurement=self.config.get('measurement',{})
        if measurement.get('discharge_plane')!='ZY' or measurement.get('slice_method')!='central_slab_projection' or float(measurement.get('plane_x_sim',0))!=protocol['plane_x_sim']:problems.append('measurement must use one central ZY observation band')
        if measurement.get('face')!='inner_toward_opening':problems.append('measurement face must be inner_toward_opening')
        for name,allow_zero in [('section_half_width_physical_m',False),('height_band_factor',False),('minimum_height_particle_diameters',True)]:
            try:number=float(measurement.get(name))
            except (TypeError,ValueError):number=float('nan')
            if not math.isfinite(number) or number<0 or (not allow_zero and number==0):problems.append('invalid measurement '+name)
        for name,spec in self.config['parameters'].items():
            if len(re.findall(PARAMETER_LINE.format(name=re.escape(name)),text,re.MULTILINE))!=1:problems.append('missing/ambiguous tunable '+name)
            if not (math.isfinite(float(spec['min'])) and math.isfinite(float(spec['max'])) and float(spec['min'])<float(spec['max']) and float(spec['min'])<=float(spec['initial'])<=float(spec['max'])):problems.append('invalid parameter bounds '+name)
            if not math.isclose(float(spec['initial']),protocol['parameters'][name],abs_tol=1e-12):problems.append('initial differs from source deck '+name)
        for name,pattern in SEED_PATTERNS.items():
            if len(pattern.findall(text))!=1:problems.append('missing/ambiguous seed '+name)
        ranks=int(self.config.get('mpi_processes',1));grid=re.search(r'(?m)^\s*processors\s+(\d+)\s+(\d+)\s+(\d+)',text)
        if not grid or math.prod(map(int,grid.groups()))!=ranks or ranks!=int(self.config.get('mpi_processes_per_trial',ranks)):problems.append('MPI processors mismatch')
        available=len(os.sched_getaffinity(0)) if hasattr(os,'sched_getaffinity') else (os.cpu_count() or 1)
        if not 1<=int(self.config['parallel_trials'])<=available:problems.append('worker count exceeds available CPUs')
        if os.name!='nt' and not Path(self.config.get('solver','/usr/local/bin/liggghts')).is_file() and self.solver_command is None:problems.append('solver not found')
        if problems:raise RuntimeError('preflight failed: '+'; '.join(problems))
        return {'status':'ok','protocol':'repose_then_base_opening','case':CASE_DIRECTORY,'input':INPUT_NAME,'solvers_per_trial':1,'snapshots':SNAPSHOTS,'snapshot_times_s':protocol['snapshot_times_s'],'discharge_plane':protocol['measurement_plane'],'nbrratio':protocol['nbrratio'],'opening_gap_physical_m':protocol['opening_gap_physical_m'],'parallel_trials':self.config['parallel_trials'],'ranks_per_trial':ranks,'available_cpus':available,'campaign':self.config['campaign'],'targets_protocol_confirmed':bool(self.config.get('targets_protocol_confirmed',False)),'measurement_approved':bool(self.config.get('measurement_approved',False)),'discharge_observation':measurement}

    def prepare_trial(self,trial_id,parameters,seeds=None):
        workspace=self.run_root/f'trial_{trial_id:04d}'
        workspace.mkdir(parents=True,exist_ok=False);case=workspace/CASE_DIRECTORY;case.mkdir()
        source=self.project_root/CASE_DIRECTORY
        shutil.copy2(source/INPUT_NAME,case/INPUT_NAME)
        shutil.copytree(source/'meshes',case/'meshes')
        seeds=seeds or self._trial_seeds(trial_id)
        rewrite_input_parameters(case/INPUT_NAME,parameters);rewrite_input_seeds(case/INPUT_NAME,seeds)
        protocol=inspect_case(case)
        (workspace/'parameters.json').write_text(json.dumps({'parameters':parameters,'seeds':seeds,'protocol':'repose_then_base_opening','source_input':str(source/INPUT_NAME),'source_input_sha256':hashlib.sha256((source/INPUT_NAME).read_bytes()).hexdigest(),'trial_input_sha256':hashlib.sha256((case/INPUT_NAME).read_bytes()).hexdigest(),'snapshot_times_s':protocol['snapshot_times_s']},indent=2))
        return workspace

    def refresh_live_state(self,run_root,active):
        for trial_id,state in active.items():
            trial=Path(run_root)/f'trial_{trial_id:04d}';log=trial/'sequential.log'
            if not log.is_file():continue
            t=simulation_time_seconds(trial/CASE_DIRECTORY/INPUT_NAME,log)
            state['simulation_time_s']=t
            state['status']='repose' if t<13. else ('opening' if t<16. else 'post_opening')

    @staticmethod
    def _nbr_note_for_snapshot(path):
        case=Path(path).parents[1]
        if not (case/INPUT_NAME).is_file():return 'Sequential source deck unavailable'
        p=inspect_case(case)
        return f'nbrratio={p["nbrratio"]:.8g}; repose13s -> opening13-16s -> discharge23s; ZY X=0'

    def _measure_snapshot(self,path,kind):
        path=Path(path);raw=path.read_bytes()
        if not raw or b'\x00' in raw:raise RuntimeError(f'corrupt or incomplete VTK (empty/NUL bytes): {path}')
        try:
            text=raw.decode('utf-8');points=parse_vtk_points(text);radii=parse_vtk_scalar_field(text,'radius')
        except (UnicodeError,ValueError) as exc:raise RuntimeError(f'invalid snapshot artifact {path}: {exc}') from exc
        if len(points)!=len(radii) or not points:raise RuntimeError('snapshot point/radius counts mismatch or empty')
        workspace=path.parents[2];label=workspace.name;output=workspace/('medicao_repouso_oito_direcoes.svg' if kind=='repose' else 'medicao_descarga.svg')
        protocol=inspect_case(path.parents[1]);result=None
        try:
            if kind=='repose':result=measure_repose_eight_ray_least_squares(points,particle_radii=radii)
            elif kind=='discharge':
                settings=self.config.get('measurement',{})
                result=measure_opening_section(points,radii,protocol,settings.get('minimum_height_particle_diameters',0.0),section_half_width_physical=settings.get('section_half_width_physical_m',.075),height_band_factor=settings.get('height_band_factor',.25))
            else:raise ValueError('unknown response '+kind)
        except UnmeasurableAngle as exc:
            result=getattr(exc,'diagnostic_result',None)
            if kind=='repose':
                if result is not None:write_repose_svg(output,result,label,str(exc),points=points,historical_preview=False)
                else:write_repose_failure_diagnostic_svg(output,points,label,str(exc),particle_radii=radii)
            else:
                from optimizer.opening_visuals import write_opening_diagnostic
                write_opening_diagnostic(output,points,radii,protocol,label,str(exc),diagnostic_result=result,settings=self.config.get('measurement',{}))
            write_svg_png(output,output.with_suffix('.png'))
            (workspace/f'{kind}_measurement.json').write_text(json.dumps({'valid':False,'reason':str(exc),'source_sha256':hashlib.sha256(raw).hexdigest(),'snapshot_time_s':protocol['snapshot_times_s'][kind],'diagnostic_result':result},indent=2))
            raise GeometryDiagnostic(f'{kind}: surface geometry measurement unavailable: {exc}') from exc
        # Preserve the numerical result even when a subsequent renderer fails.
        try:
            (workspace/f'{kind}_measurement.json').write_text(json.dumps({'valid':True,'source_sha256':hashlib.sha256(raw).hexdigest(),'snapshot_time_s':protocol['snapshot_times_s'][kind],'result':result},indent=2))
            if kind=='repose':write_repose_svg(output,result,label,self._nbr_note_for_snapshot(path),points=points,historical_preview=False)
            else:
                from optimizer.opening_visuals import write_opening_svg
                write_opening_svg(output,result,protocol,label)
            write_svg_png(output,output.with_suffix('.png'))
        except Exception as exc:
            exc.measurement_result=result
            raise
        return result

    def execute_trial(self,trial_id,parameters,phase='doe',seeds=None):
        seeds=seeds or self._trial_seeds(trial_id);started=_utc_now();workspace=self.prepare_trial(trial_id,parameters,seeds)
        base={'trial_id':trial_id,'phase':phase,'workspace':str(workspace),'started_at':started}
        self._record_status(**base,status='running',stage='sequential',message='one solver; repose then opening',finished_at='')
        try:
            self._run_input(workspace,CASE_DIRECTORY,INPUT_NAME,'sequential')
            responses={};failures={};artifact_errors=[]
            for kind in ('repose','discharge'):
                try:responses[kind]=self._measure_snapshot(workspace/CASE_DIRECTORY/SNAPSHOTS[kind],kind)
                except Exception as exc:
                    geometry=isinstance(exc,GeometryDiagnostic)
                    failures[kind]={'type':type(exc).__name__,'category':'geometry' if geometry else 'artifact_or_postprocess','message':str(exc)}
                    available=getattr(exc,'measurement_result',None)
                    if available is not None:responses[kind]=available
                    if not geometry:artifact_errors.append(exc)
            if failures:
                (workspace/'partial_responses.json').write_text(json.dumps({'valid_pair':False,'angles_deg':{k:v['mean_angle_deg'] for k,v in responses.items()},'failures':failures},indent=2))
                if artifact_errors:raise artifact_errors[0]
                raise GeometryDiagnostic('; '.join(f['message'] for f in failures.values()))
            repose,discharge=responses['repose'],responses['discharge']
            row={**parameters,'trial_id':trial_id,'phase':phase,'status':'complete','workspace':str(workspace),'distribution_seed':seeds['distribution'],'insertion_seed':seeds['insertion'],'repose_deg':repose['mean_angle_deg'],'discharge_deg':discharge['mean_angle_deg'],'repose_time_s':13.,'discharge_time_s':23.,'discharge_left_deg':discharge['left_angle_deg'],'discharge_right_deg':discharge['right_angle_deg'],'discharge_asymmetry_deg':discharge['asymmetry_deg'],'discharge_sections':1,'discharge_plane':'ZY','discharge_plane_x_m':0.,'discharge_slice_method':discharge.get('slice_method','central_slab_projection'),'discharge_band_half_width_m':discharge.get('section_half_width_physical',.075),'discharge_profile_method':discharge.get('profile_method','inner_circle_surface_consensus'),**repose_pair_result_fields(repose)}
            for side in ('left','right'):
                for n,seg in enumerate(discharge.get(side+'_segments',[]),start=1):
                    row[f'discharge_{side}_segment{n}_deg']=seg['angle_deg'];row[f'discharge_{side}_segment{n}_support']=seg['fit_support']
            for kind in ('repose','discharge'):row[kind+'_error']=abs(row[kind+'_deg']-self.config['target_'+kind+'_deg'])
            self._append_csv('results.csv',OPENING_FIELDS,row);self._record_status(**base,status='complete',stage='measure',message='two snapshots from same simulation measured',finished_at=_utc_now())
            return row
        except GeometryDiagnostic as exc:
            self._record_status(**base,status='diagnostic',stage='measurement',message=str(exc),finished_at=_utc_now());raise
        except Exception as exc:
            self._record_status(**base,status='failed',stage='artifact_or_solver',message=str(exc),finished_at=_utc_now());raise RuntimeError(f'trial {trial_id}: {exc}') from exc

    def score_completed_results(self):
        path=self.run_root/'results.csv'
        if not path.is_file():return []
        with path.open(newline='') as f:rows=list(csv.DictReader(f))
        if not rows:return []
        maxima={k:max(float(r[k+'_error']) for r in rows) for k in ('repose','discharge')}
        for r in rows:r['score']=weighted_normalized_score({k:float(r[k+'_error']) for k in maxima},maxima,self.config['weights'])
        with path.open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=OPENING_FIELDS);w.writeheader();w.writerows(rows)
        return rows

def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__,allow_abbrev=False)
    parser.add_argument('--project-root',type=Path,default=Path(__file__).resolve().parents[1]);parser.add_argument('--run-root',type=Path)
    modes=parser.add_mutually_exclusive_group(required=True)
    for flag in ('preflight','prepare-only','run'):modes.add_argument('--'+flag,action='store_true')
    parser.add_argument('--parallel-trials',type=int);args=parser.parse_args(argv)
    run=args.run_root or args.project_root/(('calibration_runs_' if args.run else 'verification_')+datetime.now().strftime('%Y%m%d_%H%M%S'))
    if args.run and run.exists() and any(run.iterdir()):raise RuntimeError('fresh run-root required; existing results will not be reused')
    c=OpeningController(args.project_root,run,parallel_trials_override=args.parallel_trials)
    report=c.preflight();print(json.dumps(report,indent=2))
    if args.preflight:return 0
    if args.prepare_only:
        workspace=c.prepare_trial(1,{n:s['initial'] for n,s in c.config['parameters'].items()});print(workspace);return 0
    if not c.config.get('targets_protocol_confirmed',False):raise RuntimeError('Confirm experimental targets for the new Jenike opening protocol before --run; preserved targets are provisional.')
    if not c.config.get('measurement_approved',False):raise RuntimeError('Approve the single-cut Jenike measurement preview before --run; no campaign started.')
    provenance={str(p.relative_to(c.project_root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [c.project_root/CASE_DIRECTORY/INPUT_NAME,*sorted((c.project_root/CASE_DIRECTORY/'meshes').glob('*.stl')),*sorted((c.project_root/'optimizer').glob('*.py')),c.project_root/'configuracao_otimizacao.txt']}
    (run/'source_provenance.json').write_text(json.dumps({'protocol':report,'source_root':str(c.project_root),'sha256':provenance},indent=2))
    rows=c.run_campaign();print(json.dumps({'status':'complete','rows':len(rows),'run_root':str(run)}));return 0
if __name__=='__main__':sys.exit(main())
