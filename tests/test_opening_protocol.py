from pathlib import Path
from unittest.mock import patch
import math,unittest
from optimizer import opening_protocol
ROOT=Path(__file__).resolve().parents[1]
CASE=ROOT/'Ang Repouso e Desgarga'

class ProtocolTests(unittest.TestCase):
    def test_original_mesh_motion_and_snapshot_times(self):
        p=opening_protocol.inspect_case(CASE)
        self.assertEqual(p['snapshot_times_s'],{'repose':13.0,'discharge':23.0})
        self.assertAlmostEqual(p['nbrratio'],31.25)
        self.assertAlmostEqual(p['opening_displacement_sim'],2.34375)
        self.assertAlmostEqual(p['opening_gap_physical_m'],.15)
        self.assertAlmostEqual(p['meshes']['base_repouso_esq']['final_bounds']['z'][1],-2.34375)
        self.assertAlmostEqual(p['meshes']['base_repouso_dir']['final_bounds']['z'][0],2.34375)
        self.assertEqual(p['measurement_plane'],'ZY at X=0')

    def inspect_mutation(self,old,new):
        read_text=Path.read_text
        text=(CASE/opening_protocol.INPUT_NAME).read_text(encoding='utf-8')
        self.assertIn(old,text)
        mutated=text.replace(old,new)
        def read(path,*args,**kwargs):
            return mutated if path==CASE/opening_protocol.INPUT_NAME else read_text(path,*args,**kwargs)
        with patch.object(Path,'read_text',new=read):
            return opening_protocol.inspect_case(CASE)

    def test_reject_base_mesh_translation(self):
        with self.assertRaises(ValueError):
            self.inspect_mutation('scale $(0.4*v_nbrratio)','scale $(0.4*v_nbrratio) move 0 1 0')

    def test_reject_left_base_motion_in_x(self):
        with self.assertRaises(ValueError):
            self.inspect_mutation('linear 0.0 0.0 -${v_base_z}','linear -${v_base_z} 0.0 0.0')

    def test_reject_shortened_executed_opening_run(self):
        with self.assertRaises(ValueError):
            self.inspect_mutation('run ${steps_mov_base}','run 1')

    def test_reject_gravity_flip(self):
        with self.assertRaises(ValueError):
            self.inspect_mutation('gravity 9.81 vector 0.0 -1.0 0.0','gravity 9.81 vector 0.0 1.0 0.0')

    def test_reject_other_mesh_transforms(self):
        for mesh,scale in [('tubo','1.5'),('inlet','1.5'),('base_repouso_esq','0.4')]:
            with self.subTest(mesh=mesh),self.assertRaises(ValueError):
                self.inspect_mutation(f'file meshes/{mesh}.stl type '+('4' if mesh=='tubo' else '1' if mesh=='inlet' else '2')+f' scale $({scale}*v_nbrratio)',
                    f'file meshes/{mesh}.stl type '+('4' if mesh=='tubo' else '1' if mesh=='inlet' else '2')+f' scale $({scale}*v_nbrratio) rotate axis 0 1 0 angle 90')

    def test_reject_motion_target_identity_change(self):
        with self.assertRaises(ValueError):
            self.inspect_mutation('move_base_esq all move/mesh mesh base_repouso_esq','move_base_esq all move/mesh mesh base_repouso_dir')

    def test_reject_missing_motion_unfix(self):
        with self.assertRaises(ValueError):
            self.inspect_mutation('unfix move_base_esq','# removed unfix move_base_esq')

    def test_reject_late_motion_unfix(self):
        with self.assertRaises(ValueError):
            self.inspect_mutation('unfix move_tubo','run 1\nunfix move_tubo')

    def test_reject_unsupported_execution_control(self):
        for command in ['jump SELF phase','include other.deck','reset_timestep 0','run 1 upto']:
            with self.subTest(command=command),self.assertRaises(ValueError):
                self.inspect_mutation('run ${steps_mov_base}',command)

    def test_reject_changed_timestep_or_mesh_file_identity(self):
        for old,new in [('timestep ${dt}','timestep 0.001'),('file meshes/base_repouso_esq.stl','file meshes/base_repouso_dir.stl')]:
            with self.subTest(mutation=new),self.assertRaises(ValueError):self.inspect_mutation(old,new)

    def test_reject_noncanonical_phase_times_even_if_declarations_agree(self):
        text=(CASE/opening_protocol.INPUT_NAME).read_text(encoding='utf-8')
        mutated=text.replace('variable tempo_mov_base equal 3.0','variable tempo_mov_base equal 4.0').replace('variable tempo_repouso_pos_descarga equal 7.0','variable tempo_repouso_pos_descarga equal 6.0')
        with self.assertRaisesRegex(ValueError,'executed phase times'):
            self.inspect_mutation(text,mutated)

    def test_evaluates_executed_run_expression(self):
        p=self.inspect_mutation('run ${steps_mov_base}','run $(v_steps_mov_base)')
        self.assertEqual(p['snapshot_times_s'],{'repose':13.,'discharge':23.})
        self.assertEqual(p['opening_end_time_s'],16.)

if __name__=='__main__':unittest.main()
