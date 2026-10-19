"""Single central ZY opening evidence in the paired repose visual language."""
import html
import math
import statistics
from pathlib import Path


def _section(points, radii, center=0.):
    pairs = [((center, float(p[1]), float(p[2])), math.sqrt(r*r-(p[0]-center)**2))
             for p, r in zip(points, radii) if abs(p[0]-center) < r]
    return [p for p, r in pairs], [r for p, r in pairs]


def write_opening_svg(path, result, protocol, label):
    _render(path, result, protocol, label, '')


def write_opening_diagnostic(path, points, radii, protocol, label, reason,
                             diagnostic_result=None, settings=None):
    if diagnostic_result is not None:
        return _render(path, diagnostic_result, protocol, label, reason)
    settings = settings or {}
    center = protocol['plane_x_sim']
    half = settings.get('section_half_width_physical_m', .075)*protocol['nbrratio']
    meshes = protocol.get('meshes', {})
    support = max((meshes[n]['initial_bounds']['y'][1]
                   for n in ('base_repouso_esq', 'base_repouso_dir') if n in meshes), default=0.)
    minimum = settings.get('minimum_height_particle_diameters', 0.)
    selected = [(p, r) for p, r in zip(points, radii)
                if abs(p[0]-center) <= half and p[1] >= support-1e-9]
    section = [(center, p[1], p[2]) for p, r in selected]
    rr = [r for p, r in selected]
    median = statistics.median(rr) if rr else 0.
    result = {'section_points': section, 'section_radii': rr, 'median_radius': median,
              'slice_method': 'central_slab_projection',
              'section_half_width_physical': settings.get('section_half_width_physical_m', .075),
              'minimum_height_particle_diameters': minimum, 'support_surface_y_sim': support}
    for side, sign in [('left', -1), ('right', 1)]:
        ys = [p[1] for p in section if sign*p[2] > 0]
        result['profile_start_height_'+side] = min(ys)+minimum*2*median if ys else support
    _render(path, result, protocol, label, reason)


def _render(path, result, protocol, label, reason):
    points, radii = result['section_points'], result['section_radii']
    ratio = protocol['nbrratio']
    minimum = result.get('minimum_height_particle_diameters', 0.)
    navy, green, orange = '#263449', '#14966c', '#ef8c24'
    sides = [('left', '#d93535', 'ESQUERDA'), ('right', '#276bd2', 'DIREITA')]
    meshes = protocol.get('meshes', {})
    support = result.get('support_surface_y_sim', max(
        (meshes[n]['initial_bounds']['y'][1] for n in ('base_repouso_esq', 'base_repouso_dir')
         if n in meshes), default=0.))
    axes = []
    for side, name in [('left', 'base_repouso_esq'), ('right', 'base_repouso_dir')]:
        mesh = meshes.get(name, {})
        bounds = mesh.get('final_bounds', mesh.get('initial_bounds', {}))
        if 'z' in bounds:
            axes.append((side, bounds['z'][1 if side == 'left' else 0]))
    if points:
        xlo = min([p[2]-r for p, r in zip(points, radii)]+[0.]+[h for _, h in axes])
        xhi = max([p[2]+r for p, r in zip(points, radii)]+[0.]+[h for _, h in axes])
        ylo = min([p[1]-r for p, r in zip(points, radii)]+[support])
        yhi = max([p[1]+r for p, r in zip(points, radii)]+[support])
        dx, dy = max(xhi-xlo, 1e-9), max(yhi-ylo, 1e-9)
        xlo -= dx*.04; xhi += dx*.04
        ylo -= dy*.06; yhi += dy*.06
        scale = min(1260/(xhi-xlo), 620/(yhi-ylo))
        plot_w, plot_h = (xhi-xlo)*scale, (yhi-ylo)*scale
    else:
        xlo, xhi, ylo, yhi, scale, plot_w, plot_h = -1., 1., 0., 1., 1., 1260., 260.
    left, top = (1440-plot_w)/2, 220.
    bottom = top+plot_h
    footer_y = bottom+102+sum(34+(22 if reason else 22*len(result.get(side+'_segments', []))) for side, _, _ in sides)
    height = math.ceil(footer_y+135)
    X = lambda h: left+(h-xlo)*scale
    Y = lambda y: top+(yhi-y)*scale
    out = [f'<svg xmlns="http://www.w3.org/2000/svg" width="1440" height="{height}" viewBox="0 0 1440 {height}">',
           '<rect width="100%" height="100%" fill="#f7f9fc"/>']

    def text(x, y, value, size=16, color=navy, kind=None, weight='normal'):
        data = f' data-kind="{kind}"' if kind else ''
        out.append(f'<text{data} x="{x:.2f}" y="{y:.2f}" font-size="{size}" '
                   f'font-family="monospace" font-weight="{weight}" fill="{color}">{html.escape(str(value))}</text>')

    text(48, 46, label.upper().replace('TRIAL_', 'TRIAL ')+' - DESCARGA / DESLIZAMENTO', 30, weight='bold')
    text(48, 80, 'UM CORTE CENTRAL: ZY | DUAS FACES INTERNAS | ESCALA Z/Y IGUAL', 19)
    text(48, 109, f'FAIXA X PROJETADA | MEIA-LARGURA={result.get("section_half_width_physical", .075):.3f} m | '
         f'SNAPSHOT SALVO: {result.get("snapshot_time_s", protocol["snapshot_times_s"]["discharge"]):.1f} s', 17)
    mode = 'SEM CORTE INFERIOR' if minimum == 0 else f'CORTE INFERIOR: {minimum:g} DIAMETROS ACIMA DO CENTRO BASAL'
    text(48, 137, mode+' | CONTORNO MATERIAL-AR DO VAO CENTRAL', 16)
    text(48, 164, 'ANGULO RELATIVO A HORIZONTAL | INTERCEPTO LIVRE; REFERENCIAS NAO ANCORAM O FIT', 15)
    if reason:
        text(48, 193, 'SEM ANGULO VALIDO: '+reason[:128], 15, '#b42318')
    else:
        text(48, 193, f'ESQ {result["left_angle_deg"]:.2f} GRAUS | DIR {result["right_angle_deg"]:.2f} GRAUS', 18)
    out.append(f'<rect data-kind="measurement-panel" x="{left-20}" y="{top-12}" width="{plot_w+40}" height="{plot_h+60}" rx="6" fill="#ffffff" stroke="#cbd5e1" stroke-width="2"/>')
    out.append(f'<rect data-kind="plot-area" x="{left}" y="{top}" width="{plot_w}" height="{plot_h}" fill="#f8fafc" stroke="{navy}" stroke-width="1.5"/>')
    if points:
        for i in range(1, 5):
            gy = top+plot_h*i/5
            out.append(f'<line x1="{left}" y1="{gy}" x2="{left+plot_w}" y2="{gy}" stroke="#d8e0e9" stroke-width="1"/>')
        for p, r in zip(points, radii):
            side = 'left' if p[2] < 0 else 'right'
            basal = minimum > 0 and p[1] < result.get('profile_start_height_'+side, support)
            out.append(f'<circle data-kind="{"basal-reference" if basal else "section-particle"}" cx="{X(p[2])}" cy="{Y(p[1])}" r="{r*scale}" fill="#7d8da1" fill-opacity="{.04 if basal else .07}" stroke="#7d8da1" stroke-opacity="{.12 if basal else .4}" stroke-width="1"/>')
        out.append(f'<line data-kind="center-axis" x1="{X(0)}" y1="{top}" x2="{X(0)}" y2="{bottom}" stroke="#14816b" stroke-width="2" stroke-dasharray="10 7"/>')
        out.append(f'<line data-kind="support-reference" data-height-sim="{support}" x1="{left}" y1="{Y(support)}" x2="{left+plot_w}" y2="{Y(support)}" stroke="{green}" stroke-width="2" stroke-dasharray="10 7"/>')
        for side, h in axes:
            out.append(f'<line data-kind="basal-axis" data-side="{side}" data-coordinate-sim="{h}" x1="{X(h)}" y1="{Y(support)}" x2="{X(h)}" y2="{top}" stroke="{orange}" stroke-width="1.5" stroke-dasharray="6 5"/>')
            out.append(f'<circle data-kind="geometric-reference" cx="{X(h)}" cy="{Y(support)}" r="5" fill="#ffd84d" stroke="{navy}" stroke-width="1.5"/>')
        summary_y = bottom+80
        for si, (side, color, title) in enumerate(sides):
            cutoff = result.get('profile_start_height_'+side, support)
            if minimum > 0:
                lo, hi = (xlo, 0) if side == 'left' else (0, xhi)
                out.append(f'<line data-kind="cutoff" data-height-sim="{cutoff}" x1="{X(lo)}" y1="{Y(cutoff)}" x2="{X(hi)}" y2="{Y(cutoff)}" stroke="#b42318" stroke-dasharray="8 5"/>')
            trace = result.get(side+'_profile_points', [])
            if trace:
                poly = ' '.join(f'{X(h)},{Y(y)}' for h, y in trace)
                out.append(f'<polyline data-kind="selected-profile-trace" points="{poly}" stroke="{green}" stroke-width="1" stroke-dasharray="3 4" fill="none"/>')
            for h, y in result.get(side+'_downweighted_points', []):
                xx, yy = X(h), Y(y)
                out.append(f'<path data-kind="excluded-sample" d="M {xx-5} {yy-5} l 10 10 M {xx-5} {yy+5} l 10 -10" stroke="{orange}" stroke-width="2.5"/>')
            for h, y in result.get(side+'_fit_points', []):
                out.append(f'<circle data-kind="selected-sample" cx="{X(h)}" cy="{Y(y)}" r="4" fill="{green}" stroke="#ffffff" stroke-width="1"/>')
            info = []
            for i, seg in enumerate(result.get(side+'_segments', []) if not reason else []):
                bx, by = seg['base']; m = seg['slope']; b = seg.get('intercept', by-m*bx)
                y1 = max(cutoff, by, min(p[1] for p in seg['fit_points']))
                y2 = max(p[1] for p in seg['fit_points'])
                if abs(m) > 1e-12 and y2 > y1:
                    h1, h2 = (y1-b)/m, (y2-b)/m
                    out.append(f'<line data-kind="opening-fit" data-side="{side}" data-segment="{i+1}" data-angle="{seg["angle_deg"]}" x1="{X(h1)}" y1="{Y(y1)}" x2="{X(h2)}" y2="{Y(y2)}" stroke="{color}" stroke-width="4"/>')
                info.append(f'S{i+1}: {seg["angle_deg"]:.2f} GRAUS | R2 {seg.get("r2", 0.):.3f} | '
                            f'UNICAS {seg.get("unique_particle_support", "N/D")} | AMOSTRAS {seg.get("fit_support", len(seg["fit_points"]))}')
            segments = result.get(side+'_segments', [])
            unique = result.get(side+'_unique_particle_support', 'N/D')
            ids = result.get(side+'_sample_particle_indices')
            candidates = str(len(set(ids))) if ids is not None else 'N/D'
            samples = len(result.get(side+'_fit_points', []))
            total = len(result.get(side+'_profile_points', result.get(side+'_fit_points', [])))
            r2 = result.get(side+'_r2')
            if r2 is None and len(segments) == 1:
                r2 = segments[0].get('r2')
            quality = f'{r2:.3f}' if r2 is not None else 'N/D'
            angle = result.get(side+'_angle_deg') if not reason else None
            angle_label = f'{angle:.2f} GRAUS' if angle is not None else 'SEM ANGULO VALIDO'
            text(48, summary_y, f'{title}: {angle_label} | R2 {quality} | SUPORTE UNICO {unique}/{candidates} | AMOSTRAS {samples}/{total}',
                 15, color, 'face-summary')
            for index, segment_text in enumerate(info):
                text(68, summary_y+22*(index+1), segment_text, 14, color, 'segment-summary')
            if reason:
                face_reason = (result.get(side+'_fit') or {}).get('rejection_reason') or reason
                text(68, summary_y+22, str(face_reason)[:140], 12, '#b42318', 'face-reason')
                summary_y += 22
            summary_y += 34+22*len(info)
    else:
        text(90, top+130, 'SEM PARTICULAS NA FAIXA CENTRAL; SEM SUPERFICIE OBSERVADA', 20, '#b42318')
    text(left, bottom+24, 'Z FISICO (m) | Y FISICO (m) | COORDENADAS DEM / NBRRATIO', 13)
    if not reason:
        spread = result.get('stddev_deg')
        if spread is None:
            spread = statistics.pstdev([result['left_angle_deg'], result['right_angle_deg']])
        text(48, footer_y, f'MEDIA GLOBAL: {result["mean_angle_deg"]:.2f} GRAUS | DESVIO ENTRE 2 FACES: {spread:.2f} GRAUS', 19, weight='bold')
    else:
        text(48, footer_y, 'DIAGNOSTICO: SEM MEDIA GLOBAL / SEM AJUSTE / NAO USAR NO OBJETIVO', 17, '#b42318')
    text(48, footer_y+31, 'CINZA: ESFERAS PROJETADAS | VERDE: PERFIL | X LARANJA: EXCLUIDOS | VERMELHO/AZUL: FITS', 14)
    text(48, footer_y+57, 'REFERENCIA GEOMETRICA: TRACEJADO VERDE = EIXO Z=0 E HORIZONTAL DO SUPORTE; LARANJA = BORDAS DO VAO', 13)
    text(48, footer_y+81, f'ESTADO POS-ABERTURA; NAO E INICIO DE DESLIZAMENTO | NBRRATIO={ratio:.8g} | GEOMETRIA NAO E PONTO DE AJUSTE', 13)
    out.append('</svg>')
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text('\n'.join(out), encoding='utf-8')
