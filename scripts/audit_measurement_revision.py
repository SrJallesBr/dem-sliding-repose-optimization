from __future__ import annotations

import importlib.util
import json
import statistics
import sys
from pathlib import Path

CANDIDATE = Path(r"/mnt/c/DEM/ANGLETOE_V2_v5_equal_scale_15s_12workers_sf11_2_CED11_1e7_measurement_v2_candidate_20261016")
OLD_RUN = Path(r"/mnt/c/DEM/ANGLETOE_V2_v5_equal_scale_15s_12workers_20261015/calibration_runs_20261016_081048")
sys.path.insert(0, str(CANDIDATE))

from optimizer.measurements import (  # noqa: E402
    UnmeasurableAngle,
    measure_repose_eight_ray_least_squares,
    measure_v_sliding_angle_three_cuts,
    parse_vtk_points,
    parse_vtk_scalar_field,
)
from optimizer.measurement_visuals import (  # noqa: E402
    write_discharge_svg,
    write_repose_cutoff_diagnostic_svg,
    write_repose_svg,
    write_svg_png,
)
from optimizer.calibration_controller import CalibrationController  # noqa: E402

baseline_path = OLD_RUN.parents[0] / "optimizer/measurements.py"
baseline_spec = importlib.util.spec_from_file_location("baseline_measurements", baseline_path)
if baseline_spec is None or baseline_spec.loader is None:
    raise RuntimeError(f"cannot load baseline estimator from {baseline_path}")
baseline_measurements = importlib.util.module_from_spec(baseline_spec)
baseline_spec.loader.exec_module(baseline_measurements)

rows = []
preview_files = []
for trial_id in (1, 23, 45):
    trial_dir = OLD_RUN / f"trial_{trial_id:04d}"
    discharge_path = trial_dir / "Ang Deslizamento/post/deslizamento_final.vtk"
    repose_path = trial_dir / "Ang Repouso/post/repouso_10s.vtk"
    discharge_deck = trial_dir / "Ang Deslizamento/in.deslizamento"
    if not discharge_path.is_file() or not repose_path.is_file():
        raise FileNotFoundError(f"missing stable snapshot for trial {trial_id}: {discharge_path} or {repose_path}")

    d_text = discharge_path.read_text(encoding="utf-8", errors="replace")
    d_points = parse_vtk_points(d_text)
    d_radii = parse_vtk_scalar_field(d_text, "radius")
    ratio = CalibrationController._nbr_ratio_from_deck(discharge_deck)
    bounds = CalibrationController._mesh_bounds_for_snapshot(discharge_path)
    d_result = measure_v_sliding_angle_three_cuts(
        d_points,
        d_radii,
        horizontal_axis="z",
        mesh_bounds=bounds,
        nbrratio=ratio,
    )
    row = {
        "trial": trial_id,
        "nbrratio_discharge": ratio,
        "mesh_x_bounds": list(bounds["x"]),
        "cut_centers_sim": [section["section_center"] for section in d_result["sections"]],
        "cut_centers_physical_m": [section["section_center_physical"] for section in d_result["sections"]],
        "cut_half_width_sim": [section["section_half_width"] for section in d_result["sections"]],
        "effective_band_physical_m": [section["effective_band_physical"] for section in d_result["sections"]],
        "band_clipped_by_mesh": [section["band_clipped"] for section in d_result["sections"]],
        "discharge_mean_deg": d_result["mean_angle_deg"],
        "discharge_sections": [
            {
                "center_m": section["section_center_physical"],
                "left_deg": section["left_angle_deg"],
                "right_deg": section["right_angle_deg"],
                "mean_deg": section["mean_angle_deg"],
                "left_fit_points": len(section["left_fit_points"]),
                "right_fit_points": len(section["right_fit_points"]),
                "left_r2": section["left_r2"],
                "right_r2": section["right_r2"],
            }
            for section in d_result["sections"]
        ],
    }

    r_text = repose_path.read_text(encoding="utf-8", errors="replace")
    r_points = parse_vtk_points(r_text)
    r_radii = parse_vtk_scalar_field(r_text, "radius")
    baseline_repose = baseline_measurements.measure_repose_eight_ray_least_squares(r_points, particle_radii=r_radii)
    row["repose_baseline_mean_deg"] = baseline_repose["mean_angle_deg"]
    row["repose_baseline_measurable_faces"] = len(baseline_repose["direction_profiles"])
    try:
        r_result = measure_repose_eight_ray_least_squares(r_points, particle_radii=r_radii)
    except (UnmeasurableAngle, ValueError) as exc:
        r_result = None
        row["repose_candidate_status"] = "unmeasurable"
        row["repose_candidate_error"] = str(exc)
        row["repose_candidate_minus_baseline_deg"] = None
        row["repose_measurable_faces"] = 0
        row["repose_outliers_global"] = None
        row["repose_downweighted_by_face"] = []
        row["repose_support_by_face"] = []
    else:
        row["repose_candidate_status"] = "measurable"
        row["repose_candidate_minus_baseline_deg"] = r_result["mean_angle_deg"] - baseline_repose["mean_angle_deg"]
        row["repose_mean_deg"] = r_result["mean_angle_deg"]
        row["repose_face_angles_deg"] = r_result["direction_angles_deg"]
        row["repose_measurable_faces"] = len(r_result["direction_profiles"])
        row["repose_outliers_global"] = r_result["discarded_outliers"]
        row["repose_downweighted_by_face"] = [p["downweighted_count"] for p in r_result["direction_profiles"]]
        row["repose_support_by_face"] = [p["fit_support"] for p in r_result["direction_profiles"]]
    rows.append(row)

    if trial_id == 45:
        preview = CANDIDATE / "preview_trial_0045"
        preview.mkdir(parents=True, exist_ok=True)
        d_svg = preview / "medicao_descarga_tres_cortes.svg"
        d_png = preview / "medicao_descarga_tres_cortes.png"
        write_discharge_svg(d_svg, d_points, d_result, "DOE trial 45 — dump final", "NBR e coordenadas do deck trial-local")
        write_svg_png(d_svg, d_png)
        if r_result is not None:
            r_svg = preview / "medicao_repouso_quatro_cortes.svg"
            r_png = preview / "medicao_repouso_quatro_cortes.png"
            write_repose_svg(r_svg, r_result, "DOE trial 45 — repouso 10 s", "NBR e escala do dump original", points=r_points, historical_preview=True)
        else:
            r_svg = preview / "medicao_repouso_quatro_cortes.svg"
            r_png = preview / "medicao_repouso_quatro_cortes.png"
            diagnostic_cutoff = min(point[1] for point in r_points) + 4.0 * statistics.median(r_radii)
            write_repose_cutoff_diagnostic_svg(
                r_svg, r_points, "DOE trial 45 — repouso 10 s",
                diagnostic_cutoff, row["repose_candidate_error"], particle_radii=r_radii,
            )
        write_svg_png(r_svg, r_png)
        preview_files = [str(d_png), str(r_png), str(d_svg), str(r_svg)]

print(json.dumps({"rows": rows, "preview_files": preview_files}, indent=2, ensure_ascii=False))
