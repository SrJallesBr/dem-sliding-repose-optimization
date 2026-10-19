"""Single central ZY observation of the supplied opening test; repose is unchanged."""
import math,statistics
from optimizer.measurements import UnmeasurableAngle,measure_v_sliding_angle_by_height

def measure_opening_section_legacy_height(points,radii,protocol,minimum_height_particle_diameters=2.0):
    if len(points)!=len(radii):raise ValueError('one radius per particle required')
    center=float(protocol['plane_x_sim']);section=[];section_radii=[];ids=[]
    meshes=protocol.get('meshes',{})
    support_y=max((meshes[n]['initial_bounds']['y'][1] for n in ('base_repouso_esq','base_repouso_dir') if n in meshes),default=0.)
    below=0
    for i,(point,radius) in enumerate(zip(points,radii)):
        if len(point)!=3 or any(not math.isfinite(float(v)) for v in point) or not math.isfinite(radius) or radius<=0:raise ValueError('finite XYZ and positive radii required')
        normal=point[0]-center
        if abs(normal)>=radius:continue
        if point[1]<support_y-1e-9:
            below+=1;continue
        section.append((center,float(point[1]),float(point[2])))
        section_radii.append(math.sqrt(radius*radius-normal*normal));ids.append(i)
    if len(section)<10:raise UnmeasurableAngle('single ZY plane has fewer than ten sphere intersections')
    result=measure_v_sliding_angle_by_height(section,section_radii,horizontal_axis='z',mesh_bounds=protocol['final_base_bounds_sim'],section_center=center,section_half_width=1e-9,minimum_height_particle_diameters=minimum_height_particle_diameters,exposed_surface='inner',reference_radius=statistics.median(radii[i] for i in ids))
    result.update(response_definition='post_opening_inner_free_surface',slice_method='sphere_plane_intersection',section_count=1,section_center_physical=center/protocol['nbrratio'],section_half_width_physical=0.,nbrratio=protocol['nbrratio'],section_original_particle_indices=ids,section_points=section,section_radii=section_radii,snapshot_time_s=protocol['snapshot_times_s']['discharge'],exposed_surface='inner',support_surface_y_sim=support_y,excluded_below_support_centers=below)
    return result


def has_sustained_bridge(heights,step,radius):
    start=previous=None
    for h in sorted(heights):
        if previous is None or h-previous>step*(1+1e-8):start=h
        if h-start>=2*radius-1e-10*radius:return True
        previous=h
    return False


def measure_opening_section(points,radii,protocol,minimum_height_particle_diameters=0.0,section_half_width_physical=.075,height_band_factor=.25):
    """One central slab, projected circles, full-height inner contour by default.

    Retains the previous discharge observation half-width, but only X=0.
    Uses physical surface arcs and the unchanged geometric-consensus fit engine.
    This is not a zero-thickness sphere-plane observation or an OLS TLS fit.
    """
    from optimizer.repose_profile import fit_repose_profile
    if len(points)!=len(radii) or not points:raise ValueError('one radius per particle and nonempty cloud required')
    if not math.isfinite(section_half_width_physical) or section_half_width_physical<=0 or not math.isfinite(height_band_factor) or height_band_factor<=0:raise ValueError('positive finite observation/sampling scales required')
    ratio=float(protocol['nbrratio']);center=float(protocol['plane_x_sim']);half=section_half_width_physical*ratio
    selected=[]
    meshes=protocol.get('meshes',{})
    support=max((meshes[n]['initial_bounds']['y'][1] for n in ('base_repouso_esq','base_repouso_dir') if n in meshes),default=0.)
    for i,(p,r) in enumerate(zip(points,radii)):
        if len(p)!=3 or any(not math.isfinite(float(v)) for v in p) or not math.isfinite(float(r)) or r<=0:raise ValueError('finite XYZ and positive radii required')
        if abs(p[0]-center)<=half and p[1]>=support-1e-9:selected.append((i,p,float(r)))
    if len(selected)<10:raise UnmeasurableAngle('insufficient material in the single central observation band')
    median=statistics.median(r for _,p,r in selected);dy=height_band_factor*median
    result={'response_definition':'post_opening_inner_free_surface','profile_method':'inner_circle_surface_consensus','slice_method':'central_slab_projection','section_count':1,'horizontal_axis':'z','section_center':center,'section_center_physical':center/ratio,'section_half_width_physical':section_half_width_physical,'section_half_width':half,'section_normal_bounds_sim':[center-half,center+half],'section_original_particle_indices':[i for i,p,r in selected],'section_points':[(center,p[1],p[2]) for i,p,r in selected],'section_radii':[r for i,p,r in selected],'median_radius':median,'height_bin_size':dy,'height_band_factor':height_band_factor,'minimum_height_particle_diameters':minimum_height_particle_diameters,'snapshot_time_s':protocol['snapshot_times_s']['discharge'],'exposed_surface':'inner','support_surface_y_sim':support,'nbrratio':ratio,'centerline':0.0}
    reasons=[]
    for side,sign in [('left',-1),('right',1)]:
        cloud=[(i,p,r) for i,p,r in selected if sign*p[2]>0]
        if not cloud:raise UnmeasurableAngle('no material on '+side+' side')
        cutoff=min(p[1] for i,p,r in cloud)+minimum_height_particle_diameters*2*median
        eligible=[(i,p,r) for i,p,r in cloud if p[1]>=cutoff]
        result['profile_start_height_'+side]=cutoff
        samples=[];sample_ids=[];evidence=[];blocked=[]
        if eligible:
            top=max(p[1]+r for i,p,r in eligible)
            levels=(top-cutoff)/dy
            if not math.isfinite(levels) or levels>512:
                raise UnmeasurableAngle('height grid exceeds 512 levels; review physical sampling resolution')
            for k in range(1+int(math.floor(levels))):
                y=cutoff+k*dy;options=[];crossing=False
                for i,p,r in eligible:
                    if abs(y-p[1])>r:continue
                    h=p[2]-sign*math.sqrt(max(0.,r*r-(y-p[1])**2))
                    if sign*h<=0:crossing=True;continue
                    options.append((sign*h,i,h,p,r))
                if crossing:
                    blocked.append(y);continue
                if not options:continue
                _,i,h,p,r=min(options,key=lambda v:(v[0],v[1]))
                samples.append((h,y));sample_ids.append(i);evidence.append({'particle_index':i,'surface':[h,y],'center':[p[2],p[1]],'radius':r,'distance_to_axis':abs(h)})
        if len(samples)<5:
            fit=None;reason='insufficient above-cutoff inner contour support'
        else:
            fit=fit_repose_profile(samples,median,sample_ids)
            reason=fit['rejection_reason']
            if not reason and sign*fit['horizontal_gradient']<=0:reason='inner face does not retreat from gap with increasing height'
            if has_sustained_bridge(blocked,dy,median):reason='material bridges the centerline over a sustained height range'
        result[side+'_blocked_heights']=blocked
        result[side+'_profile_points']=samples;result[side+'_sampling_evidence']=evidence;result[side+'_sample_particle_indices']=sample_ids
        result[side+'_unique_particle_support']=len({pid for pid,keep in zip(sample_ids,fit['inlier_mask']) if keep}) if fit else 0
        result[side+'_segments']=[];result[side+'_fit_points']=fit['inlier_points'] if fit else []
        result[side+'_downweighted_points']=fit['downweighted_points'] if fit else []
        result[side+'_angle_deg']=fit['angle_deg'] if fit and not reason else None
        result[side+'_r2']=fit['r2'] if fit else None
        result[side+'_fit']=fit
        if reason:reasons.append(side+': '+reason)
        elif fit:
            start=min(y for h,y in fit['inlier_points']);base=[(start-fit['intercept'])/fit['slope'],start]
            result[side+'_segments']=[{'base':base,'slope':fit['slope'],'intercept':fit['intercept'],'angle_deg':fit['angle_deg'],'fit_points':fit['inlier_points'],'fit_support':len(fit['inlier_points']),'unique_particle_support':result[side+'_unique_particle_support'],'r2':fit['r2'],'fit_method':'free_intercept_surface_consensus','validity_policy':fit['validity_policy']}]
    result['valid']=not reasons
    result['mean_angle_deg']=(result['left_angle_deg']+result['right_angle_deg'])/2 if not reasons else None
    result['asymmetry_deg']=abs(result['left_angle_deg']-result['right_angle_deg']) if not reasons else None
    if reasons:
        exc=UnmeasurableAngle('; '.join(reasons));exc.diagnostic_result=result;raise exc
    return result
