"""Read the supplied sequential deck and moved STL geometry without executing it."""
import ast,math,re
from pathlib import Path
from optimizer.measurements import parse_ascii_stl_bounds
CASE_DIRECTORY='Ang Repouso e Desgarga'
INPUT_NAME='in.repousodesgarga'
SNAPSHOTS={'repose':'post/repouso_final.vtk','discharge':'post/descarga_final.vtk'}

def deck_values(text):
    definitions={}
    for line in text.splitlines():
        match=re.match(r'\s*variable\s+(\w+)\s+(?:index|equal)\s+([^#]+)',line)
        if match:
            if match[1] in definitions:raise ValueError('duplicate variable '+match[1])
            definitions[match[1]]=match[2].strip()
    def evaluate(expression,stack=()):
        expression=re.sub(r'\$\{(\w+)\}',lambda m:'v_'+m[1],expression.strip())
        if expression.startswith('$(') and expression.endswith(')'):expression=expression[2:-1]
        def visit(n):
            if isinstance(n,ast.Expression):return visit(n.body)
            if isinstance(n,ast.Constant) and type(n.value) in (int,float):return float(n.value)
            if isinstance(n,ast.Name) and n.id.startswith('v_'):
                name=n.id[2:]
                if name in stack or name not in definitions:raise ValueError('unknown/cyclic variable '+name)
                return evaluate(definitions[name],stack+(name,))
            if isinstance(n,ast.UnaryOp) and isinstance(n.op,(ast.USub,ast.UAdd)):
                return (-1 if isinstance(n.op,ast.USub) else 1)*visit(n.operand)
            if isinstance(n,ast.BinOp) and isinstance(n.op,(ast.Add,ast.Sub,ast.Mult,ast.Div,ast.Pow)):
                a,b=visit(n.left),visit(n.right)
                if isinstance(n.op,ast.Add):return a+b
                if isinstance(n.op,ast.Sub):return a-b
                if isinstance(n.op,ast.Mult):return a*b
                if isinstance(n.op,ast.Div):return a/b
                if abs(b)>16:raise ValueError('excessive exponent')
                return a**b
            if isinstance(n,ast.Call) and isinstance(n.func,ast.Name) and n.func.id=='round' and len(n.args)==1 and not n.keywords:return float(round(visit(n.args[0])))
            raise ValueError('unsupported expression '+expression)
        result=visit(ast.parse(expression,mode='eval'))
        if not math.isfinite(result):raise ValueError('nonfinite expression '+expression)
        return result
    return definitions,evaluate

def inspect_case(case):
    case=Path(case);text=(case/INPUT_NAME).read_text(encoding='utf-8');definitions,evaluate=deck_values(text)
    value=lambda name:evaluate('v_'+name)
    ratio=value('nbrratio');dt=None;elapsed=0.;times={};opening_end=None
    if ratio<=0:raise ValueError('nonpositive NBR ratio')
    mesh_specs={'tubo':('mesh/surface','meshes/tubo.stl','4'),
        'base_repouso_esq':('mesh/surface','meshes/base_repouso_esq.stl','2'),
        'base_repouso_dir':('mesh/surface','meshes/base_repouso_dir.stl','2'),
        'ins_mesh':('mesh/surface/planar','meshes/inlet.stl','1')}
    motions={'move_tubo':('tubo',(0.,value('v_tubo_y'),0.)),
        'move_base_esq':('base_repouso_esq',(0.,0.,-value('v_base_z'))),
        'move_base_dir':('base_repouso_dir',(0.,0.,value('v_base_z')))}
    expected=[('run','fillsteps'),('unfix','ins'),('run','steps_ate_subida_tubo'),
        ('fix','move_tubo'),('run','steps_tubo_subida'),('unfix','move_tubo'),
        ('snapshot','repose'),('fix','move_base_esq'),('fix','move_base_dir'),
        ('run','steps_mov_base'),('unfix','move_base_esq'),('unfix','move_base_dir'),
        ('run','steps_pos_descarga'),('snapshot','discharge')]
    # This is a deliberately narrow static interpreter, not a general LIGGGHTS parser.
    # Anything that could bypass the straight-line execution model is unsupported.
    setup_commands={'atom_style','atom_modify','boundary','units','communicate','newton',
        'processors','region','create_box','neighbor','neigh_modify','pair_style','pair_coeff',
        'thermo','thermo_style','thermo_modify'}
    setup_fixes={'nve/sphere','check/timestep/gran','property/global','wall/gran',
        'particletemplate/sphere','particledistribution/discrete','insert/stream'}
    meshes={};active={};events=[];fix_ids=set();gravity=False
    def event(action,name):
        item=(action,name)
        if len(events)>=len(expected) or item!=expected[len(events)]:
            raise ValueError('unsupported phase/move/unfix sequence: '+str(item))
        events.append(item)
    for raw in text.splitlines():
        line=raw.split('#',1)[0].strip()
        if not line:continue
        tokens=re.findall(r'\$\([^)]*\)|\S+',line);command=tokens[0]
        if command=='variable':
            if not re.fullmatch(r'variable\s+\w+\s+(?:index|equal)\s+[^\n]+',line):raise ValueError('unsupported variable declaration')
            continue
        if command in setup_commands:
            if events:raise ValueError('setup command after execution started: '+command)
            continue
        if command=='shell':
            if line!='shell "mkdir -p post"' or events:raise ValueError('unsupported shell command')
            continue
        if command=='timestep':
            if len(tokens)!=2 or dt is not None or events:raise ValueError('unsupported timestep command')
            dt=evaluate(tokens[1])
            if dt<=0 or not math.isclose(dt,value('dt'),rel_tol=0.,abs_tol=1e-12):raise ValueError('executed timestep differs from dt')
            continue
        if command=='fix':
            if len(tokens)<4 or tokens[2]!='all' or tokens[1] in fix_ids:raise ValueError('unsupported/duplicate fix')
            name,style=tokens[1],tokens[3];fix_ids.add(name)
            if style.startswith('mesh/surface'):
                if events or name not in mesh_specs:raise ValueError('unsupported mesh identity '+name)
                expected_style,rel,material=mesh_specs[name]
                if len(tokens)!=10 or tokens[3:9]!=[expected_style,'file',rel,'type',material,'scale']:
                    raise ValueError('unsupported mesh declaration/transform '+name)
                path=case/rel
                if not path.resolve().is_relative_to(case.resolve()):raise ValueError('mesh outside case')
                scale=evaluate(tokens[9])
                if scale<=0:raise ValueError('nonpositive mesh scale')
                bounds={axis:[lo*scale,hi*scale] for axis,(lo,hi) in parse_ascii_stl_bounds(path.read_text()).items()}
                meshes[name]={'file':rel,'scale':scale,'initial_bounds':bounds,'offset_sim':dict.fromkeys(('x','y','z'),0.)}
            elif style=='gravity':
                if events or gravity or len(tokens)!=9 or tokens[5]!='vector':raise ValueError('unsupported gravity command')
                magnitude=evaluate(tokens[4]);vector=tuple(evaluate(v) for v in tokens[6:])
                if not math.isclose(magnitude,9.81,rel_tol=0.,abs_tol=1e-12) or vector!=(0.,-1.,0.):raise ValueError('gravity must be 9.81 in negative Y')
                gravity=True
            elif style=='move/mesh':
                if name not in motions or len(tokens)!=10 or tokens[4]!='mesh' or tokens[6]!='linear':raise ValueError('unsupported mesh motion')
                target,vector=motions[name];actual=tuple(evaluate(v) for v in tokens[7:])
                if tokens[5]!=target or target not in meshes or any(not math.isclose(a,b,rel_tol=0.,abs_tol=1e-12) for a,b in zip(actual,vector)):
                    raise ValueError('unexpected move/mesh identity/vector '+name)
                if (name=='move_tubo' and actual[1]<=0) or (name!='move_tubo' and value('v_base_z')<=0):raise ValueError('invalid motion direction')
                event('fix',name);active[name]=(target,actual)
            elif style not in setup_fixes or events:raise ValueError('unsupported fix style/phase '+style)
            continue
        if command=='unfix':
            if len(tokens)!=2:raise ValueError('unsupported unfix')
            name=tokens[1];event('unfix',name)
            if name not in fix_ids:raise ValueError('unfix of unknown fix '+name)
            if name in motions:active.pop(name)
            fix_ids.remove(name)
            continue
        if command=='run':
            if len(tokens)!=2 or dt is None or not gravity or set(meshes)!=set(mesh_specs):raise ValueError('unsupported/incomplete executed run')
            if len(events)>=len(expected) or expected[len(events)][0]!='run':raise ValueError('unexpected executed run')
            steps=evaluate(tokens[1]);phase=expected[len(events)][1]
            if steps<=0 or not steps.is_integer() or not math.isclose(steps,value(phase),rel_tol=0.,abs_tol=1e-7):raise ValueError('executed run duration differs from '+phase)
            event('run',phase);duration=steps*dt;elapsed+=duration
            for target,vector in active.values():
                for axis,velocity in zip(('x','y','z'),vector):meshes[target]['offset_sim'][axis]+=velocity*duration
            if phase=='steps_mov_base':opening_end=elapsed
            continue
        if command=='write_dump':
            if len(tokens)<5 or tokens[1:3]!=['all','custom/vtk'] or tokens[3] not in SNAPSHOTS.values():raise ValueError('unsupported snapshot')
            kind=next(k for k,v in SNAPSHOTS.items() if v==tokens[3]);event('snapshot',kind);times[kind]=elapsed
            continue
        raise ValueError('unsupported protocol command '+command)
    if events!=expected or active:raise ValueError('incomplete executed protocol')
    if any(not math.isclose(actual,wanted,rel_tol=0.,abs_tol=1e-9) for actual,wanted in [(times['repose'],13.),(opening_end,16.),(times['discharge'],23.),(elapsed,value('simulationtime')),(times['repose'],value('tempo_fim_subida_tubo'))]):
        raise ValueError('executed phase times must be repose=13, opening end=16, discharge=23 seconds')
    if not math.isclose(meshes['tubo']['offset_sim']['y'],value('desloc_tubo_y'),rel_tol=0.,abs_tol=1e-9):raise ValueError('executed tube displacement differs from declaration')
    for mesh in meshes.values():
        mesh['final_bounds']={a:[v+mesh['offset_sim'][a] for v in vs] for a,vs in mesh['initial_bounds'].items()}
    move=meshes['base_repouso_dir']['offset_sim']['z']
    bases=[meshes[n] for n in ['base_repouso_esq','base_repouso_dir']]
    bounds={a:[min(b['final_bounds'][a][0] for b in bases),max(b['final_bounds'][a][1] for b in bases)] for a in ('x','z')}
    center=(bounds['x'][0]+bounds['x'][1])/2
    if not math.isclose(center,0.,abs_tol=1e-8):raise ValueError('single-plane center requires review for non-centered geometry')
    return {'nbrratio':ratio,'snapshot_times_s':times,'opening_end_time_s':opening_end,'opening_displacement_sim':move,'opening_gap_physical_m':(meshes['base_repouso_dir']['final_bounds']['z'][0]-meshes['base_repouso_esq']['final_bounds']['z'][1])/ratio,'meshes':meshes,'final_base_bounds_sim':bounds,'measurement_plane':'ZY at X=0','plane_x_sim':center,'parameters':{n:value(n) for n in ('sf11','rf11','CED11')},'dt':value('dt')}
