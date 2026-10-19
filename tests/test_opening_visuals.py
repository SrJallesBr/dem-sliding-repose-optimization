import math,tempfile,unittest,xml.etree.ElementTree as ET
from pathlib import Path
from optimizer import opening_visuals

class OpeningVisualTests(unittest.TestCase):
    def test_each_segment_keeps_its_own_visible_angle_not_averaged_slope(self):
        result={'section_points':[(0,.2,-2),(0,2.,-1),(0,.2,2),(0,2,1)],'section_radii':[.1]*4,'median_radius':.1,'profile_start_height_left':.4,'profile_start_height_right':.4,'mean_angle_deg':37.5,'left_angle_deg':37.5,'right_angle_deg':37.5}
        for side,sign in [('left',-1),('right',1)]:
            result[side+'_fit_points']=[(sign,.5),(sign*2,1.)]
            result[side+'_segments']=[{'base':(sign,.5),'slope':sign*math.tan(math.radians(a)),'angle_deg':a,'fit_points':[(sign,.5),(sign*1.5,.8)],'fit_support':5,'r2':.95} for a in [30.,45.]]
        p={'nbrratio':31.25,'snapshot_times_s':{'repose':13.,'discharge':23.},'opening_gap_physical_m':.15}
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'opening.svg';opening_visuals.write_opening_svg(path,result,p,'fixture')
            root=ET.parse(path).getroot()
        fits=[el for el in root.iter() if el.attrib.get('data-kind')=='opening-fit']
        self.assertEqual(len(fits),4)
        for el in fits:
            a=el.attrib
            actual=math.degrees(math.atan2(abs(float(a['y2'])-float(a['y1'])),abs(float(a['x2'])-float(a['x1']))))
            self.assertAlmostEqual(actual,float(a['data-angle']),places=5)
        self.assertIn('UM CORTE',' '.join(root.itertext()))
        self.assertNotIn('3 CORTES',' '.join(root.itertext()))
    def _fixture(self):
        result={'section_points':[(0,.2,-2),(0,2.,-1),(0,.2,2),(0,2.,1)],
                'section_radii':[.1]*4,'profile_start_height_left':.4,
                'profile_start_height_right':.4,'minimum_height_particle_diameters':0.,
                'slice_method':'central_slab_projection','section_half_width_physical':.075,
                'mean_angle_deg':37.5,'left_angle_deg':30.,'right_angle_deg':45.}
        for side,sign,angle in [('left',-1,30.),('right',1,45.)]:
            result[side+'_fit_points']=[(sign,.5),(sign*1.5,.8)]
            result[side+'_segments']=[{'base':(sign,.5),'slope':sign*math.tan(math.radians(angle)),
                'angle_deg':angle,'fit_points':result[side+'_fit_points'],'fit_support':2,'r2':.95}]
        protocol={'nbrratio':31.25,'plane_x_sim':0.,'snapshot_times_s':{'repose':13.,'discharge':23.},
                  'opening_gap_physical_m':.15,'meshes':{n:{'initial_bounds':{'y':[0.,0.],
                  'z':[-3.,-1.] if n.endswith('esq') else [1.,3.]}}
                  for n in ['base_repouso_esq','base_repouso_dir']}}
        return result,protocol

    def _render(self,result,protocol,reason=None):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'opening.svg'
            if reason is None:opening_visuals.write_opening_svg(path,result,protocol,'trial_0001')
            else:opening_visuals.write_opening_diagnostic(path,[],[],protocol,'trial_0001',reason,diagnostic_result=result)
            return ET.parse(path).getroot()

    def test_zero_cutoff_uses_one_airy_monospace_repose_style_panel(self):
        result,protocol=self._fixture();root=self._render(result,protocol)
        texts=[el for el in root.iter() if el.tag.endswith('text')]
        self.assertTrue(all(el.attrib.get('font-family')=='monospace' for el in texts))
        label=' '.join(root.itertext())
        for expected in ['UM CORTE CENTRAL: ZY','DUAS FACES INTERNAS','SEM CORTE INFERIOR',
                         'MEIA-LARGURA=0.075 m','23.0 s','Z FISICO (m)','Y FISICO (m)']:
            self.assertIn(expected,label)
        self.assertNotIn('Cutoff',label);self.assertNotIn('FORA DO FIT',label)
        kinds=[el.attrib.get('data-kind') for el in root.iter()]
        self.assertEqual(kinds.count('measurement-panel'),1)
        self.assertNotIn('cutoff',kinds);self.assertNotIn('basal-reference',kinds)
        particles=[el for el in root.iter() if el.attrib.get('data-kind')=='section-particle']
        self.assertEqual(len(particles),4)
        self.assertTrue(all(float(el.attrib['fill-opacity'])<=.07 for el in particles))
        self.assertEqual(len({el.attrib['fill-opacity'] for el in particles}),1)
        self.assertTrue(all(el.attrib['stroke']=='#7d8da1' for el in particles))
        self.assertIn('center-axis',kinds);self.assertIn('support-reference',kinds)
        self.assertEqual(kinds.count('basal-axis'),2)
        self.assertIn('REFERENCIA GEOMETRICA',label)

    def test_unique_particle_support_is_not_the_sample_count(self):
        result,protocol=self._fixture()
        for side in ('left','right'):
            result[side+'_unique_particle_support']=5
            result[side+'_fit_points'] *= 6
            result[side+'_profile_points']=[(1.,.5)]*28
            result[side+'_sample_particle_indices']=list(range(7))*4
            result[side+'_segments'][0].update(fit_support=12,unique_particle_support=5)
        root=self._render(result,protocol)
        summaries=[el for el in root.iter() if el.attrib.get('data-kind')=='face-summary']
        self.assertEqual(len(summaries),2)
        for el in summaries:
            self.assertIn('SUPORTE UNICO 5/7',el.text)
            self.assertIn('AMOSTRAS 12/28',el.text)
            self.assertIn('R2 0.950',el.text)
        segments=[el for el in root.iter() if el.attrib.get('data-kind')=='segment-summary']
        self.assertEqual(len(segments),2)
        self.assertTrue(all('UNICAS 5 | AMOSTRAS 12' in el.text for el in segments))
        self.assertIn('DESVIO ENTRE 2 FACES: 7.50 GRAUS',' '.join(root.itertext()))
        for side in ('left','right'):
            result.pop(side+'_unique_particle_support')
            result.pop(side+'_sample_particle_indices')
            result[side+'_segments'][0].pop('unique_particle_support')
        fallback=self._render(result,protocol)
        self.assertIn('SUPORTE UNICO N/D',' '.join(fallback.itertext()))

    def test_default_diagnostic_also_has_no_lower_cut(self):
        result,protocol=self._fixture()
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'failure.svg'
            opening_visuals.write_opening_diagnostic(path,result['section_points'],result['section_radii'],protocol,'trial_0001','fixture')
            root=ET.parse(path).getroot()
        self.assertIn('SEM CORTE INFERIOR',' '.join(root.itertext()))
        self.assertFalse(any(e.attrib.get('data-kind') in ('cutoff','basal-reference','opening-fit') for e in root.iter()))

    def test_failure_keeps_excluded_markers_and_face_reasons_without_fits(self):
        result,protocol=self._fixture()
        for side in ('left','right'):
            result[side+'_downweighted_points']=[(1.,1.)]
            result[side+'_fit']={'rejection_reason':'unsupported fixture face'}
            result[side+'_segments']=[]
        root=self._render(result,protocol,'unsupported fixture face')
        kinds=[e.attrib.get('data-kind') for e in root.iter()]
        self.assertEqual(kinds.count('excluded-sample'),2)
        self.assertNotIn('opening-fit',kinds)
        self.assertEqual(kinds.count('face-summary'),2)
        self.assertNotIn('30.00 GRAUS',' '.join(root.itertext()))
        self.assertIn('unsupported fixture face',' '.join(root.itertext()))
        reasons=[float(e.attrib['y']) for e in root.iter() if e.attrib.get('data-kind')=='face-reason']
        footer=next(e for e in root.iter() if (e.text or '').startswith('DIAGNOSTICO:'))
        self.assertGreater(float(footer.attrib['y']),max(reasons)+25)

    def test_reduced_diagnostic_honors_support_height_and_configured_cutoff(self):
        protocol={'nbrratio':31.25,'plane_x_sim':0.,'snapshot_times_s':{'repose':13.,'discharge':23.},'opening_gap_physical_m':.15,'meshes':{n:{'initial_bounds':{'y':[1.,1.]}} for n in ['base_repouso_esq','base_repouso_dir']}}
        points=[(0,1.125,-2),(0,1.125,2),(0,.1,-.1),(0,.1,.1)]
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'failure.svg'
            opening_visuals.write_opening_diagnostic(path,points,[.125]*4,protocol,'diagnostic','test',settings={'section_half_width_physical_m':.075,'minimum_height_particle_diameters':1.})
            root=ET.parse(path).getroot()
        cuts=[float(el.attrib['data-height-sim']) for el in root.iter() if el.attrib.get('data-kind')=='cutoff']
        self.assertEqual(cuts,[1.375,1.375])
        self.assertFalse(any(el.attrib.get('data-kind')=='opening-fit' for el in root.iter()))

if __name__=='__main__':unittest.main()
