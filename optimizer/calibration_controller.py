"""Controlador portável para a campanha de calibração DEM de três respostas."""

from __future__ import annotations

import argparse
import ast
import csv
import json
import math
import os
import random
import re
import shutil
import subprocess
import sys
import time
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from datetime import datetime, timezone
from pathlib import Path
from threading import Lock

from optimizer.article_optimizer import latin_hypercube, recommend_campaign_size, weighted_normalized_score
from optimizer.cohesion_guard import CohesionLimitTracker
from optimizer.measurements import UnmeasurableAngle, measure_repose_eight_ray_least_squares, measure_v_sliding_angle, measure_v_sliding_angle_by_height, measure_v_sliding_angle_three_cuts, parse_ascii_stl_bounds, parse_vtk_points, parse_vtk_scalar_field
from optimizer.measurement_visuals import write_discharge_failure_svg, write_discharge_svg, write_measurement_diagnostic_svg, write_repose_failure_diagnostic_svg, write_repose_svg, write_svg_png

TESTS = (
    ("discharge", "Ang Deslizamento", "in.deslizamento", "post/deslizamento_final.vtk"),
    ("repose", "Ang Repouso", "in.repouso", "post/repouso_10s.vtk"),
)
PARAMETER_LINE = r"^(\s*variable\s+{name}\s+index\s+)(\S+)(.*)$"
SEED_PATTERNS = {
    "distribution": re.compile(r"(particledistribution/discrete\s+)(\d+)", re.IGNORECASE),
    "insertion": re.compile(r"(insert/stream\s+seed\s+)(\d+)", re.IGNORECASE),
}
STATUS_FIELDS = ("trial_id", "phase", "status", "stage", "message", "workspace", "started_at", "finished_at")
RESULT_FIELDS = ("trial_id", "phase", "status", "sf11", "rf11", "rf12", "CED11", "distribution_seed", "insertion_seed", "discharge_deg", "discharge_left_deg", "discharge_right_deg", "discharge_asymmetry_deg", "discharge_x0_deg", "discharge_x0_left_deg", "discharge_x0_right_deg", "discharge_xplus01_deg", "discharge_xplus01_left_deg", "discharge_xplus01_right_deg", "discharge_xminus01_deg", "discharge_xminus01_left_deg", "discharge_xminus01_right_deg", "discharge_x0_center_physical_m", "discharge_x0_band_clipped", "discharge_x0_band_min_m", "discharge_x0_band_max_m", "discharge_xplus01_center_physical_m", "discharge_xplus01_band_clipped", "discharge_xplus01_band_min_m", "discharge_xplus01_band_max_m", "discharge_xminus01_center_physical_m", "discharge_xminus01_band_clipped", "discharge_xminus01_band_min_m", "discharge_xminus01_band_max_m", *(f"discharge_{cut}_minimum_height_particle_diameters" for cut in ("x0", "xplus01", "xminus01")), *(f"discharge_{cut}_{side}_profile_start_height_m" for cut in ("x0", "xplus01", "xminus01") for side in ("left", "right")), "repose_deg", "repose_axis0_deg", "repose_axis0_plus_deg", "repose_axis0_minus_deg", "repose_axis45_deg", "repose_axis45_plus_deg", "repose_axis45_minus_deg", "repose_axis90_deg", "repose_axis90_plus_deg", "repose_axis90_minus_deg", "repose_axis135_deg", "repose_axis135_plus_deg", "repose_axis135_minus_deg", "discharge_error", "repose_error", "score", "workspace", *(f"discharge_{cut}_{side}_segment{segment}_{field}" for cut in ("x0", "xplus01", "xminus01") for side in ("left", "right") for segment in (1, 2) for field in ("count", "deg", "support")))
DIAGNOSTIC_FIELDS = ("trial_id", "phase", "status", "reason", "parameters_json", "workspace")
LIVE_STATUS_FIELDS = ("trial", "phase", "status", "simulation_time_s", "started_at", "elapsed_seconds", "updated_at")


def simulation_time_seconds(input_file: Path, log_file: Path) -> float:
    """Calcula o último tempo físico a partir de Step × dt no log do LIGGGHTS."""
    try:
        text = input_file.read_text(encoding="utf-8")
        match = re.search(r"(?m)^\s*variable\s+dt\s+equal\s+([^\s#]+)", text)
        if not match:
            return 0.0
        timestep = float(match.group(1))
        latest_step = None
        thermo = False
        for line in log_file.read_text(encoding="utf-8", errors="replace").splitlines():
            if "Step" in line and "Atoms" in line:
                thermo = True
                continue
            if thermo:
                step = re.match(r"^\s*(\d+)\s+", line)
                if step:
                    latest_step = int(step.group(1))
        return round((latest_step or 0) * timestep, 1)
    except (OSError, ValueError):
        return 0.0


def refresh_live_state(run_root: Path, active: dict[int, dict]) -> None:
    """Atualiza estágio e tempo físico a partir dos logs produzidos por cada trial ativo."""
    for trial_id, state in active.items():
        workspace = run_root / f"trial_{trial_id:04d}"
        v_input = workspace / "Ang Deslizamento" / "in.deslizamento"
        repose_input = workspace / "Ang Repouso" / "in.repouso"
        v_log = workspace / "discharge.log"
        repose_log = workspace / "repose.log"
        if repose_log.is_file():
            state["status"] = "repose"
            state["simulation_time_s"] = simulation_time_seconds(repose_input, repose_log)
        elif v_log.is_file():
            state["status"] = "discharge"
            state["simulation_time_s"] = simulation_time_seconds(v_input, v_log)


def write_live_status_csv(root: Path, active: dict[int, dict], now: float | None = None) -> None:
    """Grava status_execucao.csv atomicamente, com tentativas para pastas sincronizadas."""
    now = time.time() if now is None else now
    target = root / "status_execucao.csv"
    temporary = root / ".status_execucao.csv.tmp"
    for attempt in range(3):
        try:
            with temporary.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=LIVE_STATUS_FIELDS)
                writer.writeheader()
                for trial_id, state in sorted(active.items()):
                    started = state["started_at"]
                    writer.writerow({"trial": f"trial_{trial_id:04d}", "phase": state["phase"], "status": state["status"], "simulation_time_s": f"{state.get('simulation_time_s', 0.0):.1f}", "started_at": time.strftime("%Y-%m-%dT%H:%M:%S", time.localtime(started)), "elapsed_seconds": str(int(now - started)), "updated_at": time.strftime("%Y-%m-%dT%H:%M:%S", time.localtime(now))})
            temporary.replace(target)
            return
        except PermissionError:
            if attempt == 2:
                return
            time.sleep(0.25 * (attempt + 1))


def render_dashboard(active: dict[int, dict], completed: int, maximum: int, now: float | None = None, interval: int = 10) -> None:
    """Exibe o painel ao vivo no terminal, no mesmo formato do otimizador anterior."""
    now = time.time() if now is None else now
    print("=== PAINEL DA OTIMIZACAO DEM ===")
    print(f"Concluidos: {completed}/{maximum} | Atualiza a cada {interval}s | {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(now))}")
    print(f"{'TRIAL':<18} {'FASE':<22} {'STATUS':<18} {'SIM_TIME':<10} TEMPO REAL")
    print("-" * 82)
    for trial_id, state in sorted(active.items()):
        sim_time = f"{state.get('simulation_time_s', 0.0):.1f}".replace(".", ",")
        elapsed = int(now - state["started_at"])
        hours, remainder = divmod(elapsed, 3600)
        minutes, seconds = divmod(remainder, 60)
        print(f"{('trial_' + f'{trial_id:04d}'):<18} {state['phase']:<22} {state['status']:<18} {sim_time + 's':<10} {hours:02d}:{minutes:02d}:{seconds:02d}")
    print("\nCSV ao vivo: status_execucao.csv | Resultados finalizados: results.csv", flush=True)


def is_valid_liggghts_seed(value: int) -> bool:
    """LIGGGHTS exige sementes primas maiores que 10000; manter abaixo de 900 milhões é compatível com o gerador interno."""
    value = int(value)
    if value <= 10_000 or value >= 900_000_000:
        return False
    if value % 2 == 0:
        return False
    divisor = 3
    while divisor * divisor <= value:
        if value % divisor == 0:
            return False
        divisor += 2
    return True


def next_liggghts_prime(candidate: int) -> int:
    """Retorna a primeira semente prima válida no intervalo aceito pelo LIGGGHTS."""
    candidate = max(10_001, int(candidate) | 1)
    while candidate < 900_000_000:
        if is_valid_liggghts_seed(candidate):
            return candidate
        candidate += 2
    raise RuntimeError("could not allocate a valid LIGGGHTS random seed")


def load_configuration(path: Path) -> dict:
    """Lê a configuração raiz e deriva o tamanho da campanha pelos parâmetros."""
    config = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(config.get("parameters"), dict) or not config["parameters"]:
        raise ValueError("configuracao_otimizacao.txt requires non-empty parameters")
    defaults = recommend_campaign_size(config["parameters"])
    explicit_campaign = config.get("campaign", {})
    if explicit_campaign and not isinstance(explicit_campaign, dict):
        raise ValueError("campaign must be an object when provided")
    legacy_campaign = {name: config[name] for name in ("doe_trials", "refinement_trials", "validation_trials", "max_trials") if name in config}
    config["campaign"] = {**defaults, **legacy_campaign, **explicit_campaign}
    required_total = sum(int(config["campaign"][name]) for name in ("doe_trials", "refinement_trials", "validation_trials"))
    configured_sizes = {name for name in ("doe_trials", "refinement_trials", "validation_trials", "max_trials") if name in legacy_campaign or name in explicit_campaign}
    if {"doe_trials", "refinement_trials", "validation_trials", "max_trials"}.issubset(configured_sizes) and int(config["campaign"]["max_trials"]) != required_total:
        raise ValueError("campaign.max_trials must equal doe_trials + refinement_trials + validation_trials")
    return config


def rewrite_input_parameters(path: Path, parameters: dict[str, float]) -> None:
    """Substitui somente variáveis index explicitamente solicitadas no input copiado."""
    text = Path(path).read_text(encoding="utf-8")
    for name, value in parameters.items():
        pattern = re.compile(PARAMETER_LINE.format(name=re.escape(name)), re.MULTILINE)
        text, count = pattern.subn(lambda match: f"{match.group(1)}{value:.12g}{match.group(3)}", text, count=1)
        if count != 1:
            raise ValueError(f"variable {name!r} was not found exactly once in {path}")
    Path(path).write_text(text, encoding="utf-8")


def rewrite_input_seeds(path: Path, seeds: dict[str, int]) -> None:
    """Replaces the two stochastic seed sites in a copied input only."""
    text = Path(path).read_text(encoding="utf-8")
    for name, pattern in SEED_PATTERNS.items():
        if name not in seeds:
            raise ValueError(f"missing {name} seed")
        text, count = pattern.subn(lambda match: f"{match.group(1)}{int(seeds[name])}", text, count=1)
        if count != 1:
            raise ValueError(f"{name} seed site was not found exactly once in {path}")
    Path(path).write_text(text, encoding="utf-8")


def measure_pile_angle(points, bins=24, x_half_width=0.02, particle_radii=None) -> dict:
    """Measure four diametral repose cuts with eight base-anchored robust faces."""
    return measure_repose_eight_ray_least_squares(points, particle_radii=particle_radii)


def discharge_cut_result_fields(measurement: dict) -> dict:
    """Persist cut measurements by their physical coordinates, not section list order."""
    cut_names = {0.0: "x0", 0.1: "xplus01", -0.1: "xminus01"}
    sections_by_cut = {}
    for section in measurement.get("sections", []):
        if "section_center_physical" not in section:
            raise ValueError("discharge section is missing section_center_physical")
        center = float(section["section_center_physical"])
        matching = [value for value in cut_names if math.isclose(center, value, rel_tol=0.0, abs_tol=1e-9)]
        if len(matching) != 1:
            raise ValueError(f"unexpected discharge section center: {center}")
        cut_name = cut_names[matching[0]]
        if cut_name in sections_by_cut:
            raise ValueError(f"duplicate discharge section center: {center}")
        sections_by_cut[cut_name] = section
    if set(sections_by_cut) != set(cut_names.values()):
        raise ValueError("discharge results must contain the three physical cuts at -0.1, 0.0, and +0.1 m")
    fields = {}
    for cut_name, section in sections_by_cut.items():
        fields[f"discharge_{cut_name}_minimum_height_particle_diameters"] = float(section["minimum_height_particle_diameters"])
        fields[f"discharge_{cut_name}_left_profile_start_height_m"] = float(section["profile_start_height_left"])
        fields[f"discharge_{cut_name}_right_profile_start_height_m"] = float(section["profile_start_height_right"])
        fields[f"discharge_{cut_name}_deg"] = section["mean_angle_deg"]
        fields[f"discharge_{cut_name}_left_deg"] = section["left_angle_deg"]
        fields[f"discharge_{cut_name}_right_deg"] = section["right_angle_deg"]
        fields[f"discharge_{cut_name}_center_physical_m"] = float(section["section_center_physical"])
        fields[f"discharge_{cut_name}_band_clipped"] = section["band_clipped"]
        fields[f"discharge_{cut_name}_band_min_m"] = section["effective_band_physical"][0]
        fields[f"discharge_{cut_name}_band_max_m"] = section["effective_band_physical"][1]
        for side in ("left", "right"):
            segments = section[f"{side}_segments"]
            fields[f"discharge_{cut_name}_{side}_segment1_count"] = len(segments)
            for index in (1, 2):
                segment = segments[index - 1] if index <= len(segments) else None
                fields[f"discharge_{cut_name}_{side}_segment{index}_deg"] = segment["angle_deg"] if segment else ""
                fields[f"discharge_{cut_name}_{side}_segment{index}_support"] = segment["fit_support"] if segment else ""
    return fields


def repose_pair_result_fields(measurement: dict) -> dict:
    """Return machine-readable pair means and individual face angles for four repose axes."""
    expected_axes = {0, 45, 90, 135}
    pairs = measurement.get("pair_profiles", [])
    fields = {}
    for pair in pairs:
        axis = int(round(float(pair["axis_angle_deg"])))
        if axis not in expected_axes:
            raise ValueError(f"unexpected repose diametral axis: {axis}")
        fields[f"repose_axis{axis}_deg"] = float(pair["mean_angle_deg"])
        fields[f"repose_axis{axis}_plus_deg"] = float(pair["positive_angle_deg"])
        fields[f"repose_axis{axis}_minus_deg"] = float(pair["negative_angle_deg"])
    if {int(key.split("axis", 1)[1].split("_", 1)[0]) for key in fields if key.endswith("_deg")} != expected_axes:
        raise ValueError("repose result must contain exactly four diametral pairs")
    return fields


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


class GeometryDiagnostic(RuntimeError):
    """The solver completed but the resulting geometry does not define the requested angle."""


class MeasurementRejected(RuntimeError):
    """Identifica uma geometria de saída inadequada para leitura, apesar de o solver ter concluído."""


class CalibrationController:
    def __init__(self, project_root: Path, run_root: Path | None = None, solver_command=None, parallel_trials_override: int | None = None):
        self.project_root = Path(project_root).resolve()
        self.config = load_configuration(self.project_root / "configuracao_otimizacao.txt")
        if parallel_trials_override is not None:
            configured_workers = int(self.config.get("parallel_trials", 1))
            parallel_trials_override = int(parallel_trials_override)
            if parallel_trials_override < 1:
                raise ValueError("parallel_trials override must be positive")
            if parallel_trials_override > configured_workers:
                raise ValueError(f"parallel_trials override={parallel_trials_override} exceeds configured ceiling ({configured_workers})")
            self.config["parallel_trials"] = parallel_trials_override
        self.run_root = Path(run_root or self.project_root / "calibration_runs").resolve()
        self.run_root.mkdir(parents=True, exist_ok=True)
        self.solver_command = list(solver_command) if solver_command else None
        self._csv_lock = Lock()
        self._console_lock = Lock()
        self._cohesion_lock = Lock()
        guard = self.config.get("cohesion_guard", {})
        cohesion_bounds = self.config.get("parameters", {}).get("CED11")
        self._cohesion_tracker = CohesionLimitTracker(
            float(cohesion_bounds["min"]),
            float(cohesion_bounds["max"]),
            max_probes=int(guard.get("max_boundary_probes", 8)),
        ) if guard.get("enabled", True) and cohesion_bounds else None
        self._boundary_search_report = None
        self._next_trial_id = 1

    @property
    def parameter_names(self):
        return tuple(self.config["parameters"])

    def _write_cohesion_report(self, stage: str = "running") -> None:
        if self._cohesion_tracker is None:
            return
        with self._cohesion_lock:
            report = self._cohesion_tracker.report()
            report.update({
                "status": stage,
                "interpretation": "A geometry diagnostic is evidence of a non-measurable response, not proof that cohesion alone caused it.",
                "boundary_search": self._boundary_search_report,
            })
            if self._boundary_search_report and self._boundary_search_report.get("status") == "complete":
                report["recommended_ceiling_CED11"] = self._boundary_search_report["recommended_ceiling_CED11"]
                report["lowest_problematic_CED11"] = self._boundary_search_report["lowest_problematic_CED11"]
                report["highest_problematic_CED11"] = self._boundary_search_report["highest_problematic_CED11"]
                report["highest_measurable_CED11"] = self._boundary_search_report["highest_measurable_CED11"]
            target = self.run_root / "cohesion_limits.json"
            temporary = target.with_suffix(".json.tmp")
            temporary.write_text(json.dumps(report, indent=2, sort_keys=True), encoding="utf-8")
            temporary.replace(target)

    def _observe_cohesion(self, value: float, trial_id: int, parameters: dict, measurable: bool, reason: str = "") -> None:
        if self._cohesion_tracker is None:
            return
        with self._cohesion_lock:
            if measurable:
                self._cohesion_tracker.observe_valid(value, trial_id, parameters)
            else:
                self._cohesion_tracker.observe_problematic(value, trial_id, parameters, reason)
        self._write_cohesion_report()

    def _solver_invocation(self, cwd: Path, input_name: str) -> list[str]:
        if self.solver_command is not None:
            return [*self.solver_command, "-in", input_name]
        solver = self.config.get("solver", "/usr/local/bin/liggghts")
        mpi_processes = int(self.config.get("mpi_processes", 1))
        unbuffered_solver = ["stdbuf", "-oL", "-eL", solver]
        if os.name == "nt":
            distro = self.config.get("wsl_distribution", "Ubuntu-24.04")
            linux_cwd = subprocess.check_output(["wsl.exe", "-d", distro, "--", "wslpath", "-a", str(cwd)], text=True).strip()
            if mpi_processes > 1:
                return ["wsl.exe", "-d", distro, "--cd", linux_cwd, "--", "mpirun", "--allow-run-as-root", "--oversubscribe", "-np", str(mpi_processes), *unbuffered_solver, "-in", input_name]
            return ["wsl.exe", "-d", distro, "--cd", linux_cwd, "--", *unbuffered_solver, "-in", input_name]
        if mpi_processes > 1:
            return ["mpirun", "-np", str(mpi_processes), *unbuffered_solver, "-in", input_name]
        return [*unbuffered_solver, "-in", input_name]

    def preflight(self) -> dict:
        problems = []
        for _, directory, input_name, snapshot in TESTS:
            source = self.project_root / directory
            input_path = source / input_name
            if not input_path.is_file():
                problems.append(f"missing input: {input_path}")
                continue
            text = input_path.read_text(encoding="utf-8")
            for name in self.parameter_names:
                if not re.search(PARAMETER_LINE.format(name=re.escape(name)), text, re.MULTILINE):
                    problems.append(f"missing variable {name} in {input_path}")
            for name, expected in self.config.get("fixed_parameters", {}).items():
                fixed_pattern = re.compile(rf"(?m)^\s*variable\s+{re.escape(name)}\s+equal\s+(\S+)")
                fixed_values = fixed_pattern.findall(text)
                if len(fixed_values) != 1:
                    problems.append(f"missing or ambiguous fixed variable {name} in {input_path}")
                else:
                    try:
                        if not math.isclose(float(fixed_values[0]), float(expected), rel_tol=0.0, abs_tol=1e-12):
                            problems.append(f"fixed variable {name} is {fixed_values[0]}, expected {expected}, in {input_path}")
                    except ValueError:
                        problems.append(f"fixed variable {name} is not numeric in {input_path}")
            for seed_name, pattern in SEED_PATTERNS.items():
                if len(pattern.findall(text)) != 1:
                    problems.append(f"missing or ambiguous {seed_name} seed site in {input_path}")
            if snapshot not in text:
                problems.append(f"missing write_dump for {snapshot} in {input_path}")
        mpi_processes = int(self.config.get("mpi_processes", self.config.get("mpi_processes_per_trial", 1)))
        ranks_per_trial = int(self.config.get("mpi_processes_per_trial", mpi_processes))
        if mpi_processes < 1 or ranks_per_trial != mpi_processes:
            problems.append(f"MPI ranks mismatch: mpi_processes={mpi_processes}, mpi_processes_per_trial={ranks_per_trial}")
        for _, directory, input_name, _ in TESTS:
            deck = (self.project_root / directory / input_name).read_text(encoding="utf-8")
            grid = re.search(r"(?m)^\s*processors\s+(\d+)\s+(\d+)\s+(\d+)\b", deck)
            if not grid or math.prod(map(int, grid.groups())) != mpi_processes:
                problems.append(f"processors grid in {directory}/{input_name} must multiply to {mpi_processes}")
        available_cpus = os.cpu_count()
        requested_workers = int(self.config.get("parallel_trials", 1))
        if available_cpus is not None and requested_workers > available_cpus:
            problems.append(f"parallel_trials={requested_workers} exceeds CPUs visible to Python/WSL ({available_cpus})")
        if problems:
            raise RuntimeError("preflight failed: " + "; ".join(problems))
        return {"status": "ok", "campaign": self.config["campaign"], "parameters": list(self.parameter_names), "parallel_trials": requested_workers, "available_cpus": available_cpus}

    def prepare_trial(self, trial_id: int, parameters: dict[str, float], seeds: dict[str, int] | None = None) -> Path:
        workspace = self.run_root / f"trial_{trial_id:04d}"
        if workspace.exists():
            raise FileExistsError(f"trial workspace already exists: {workspace}")
        workspace.mkdir(parents=True)
        seeds = seeds or self._trial_seeds(trial_id)
        for _, directory, input_name, _ in TESTS:
            shutil.copytree(self.project_root / directory, workspace / directory, ignore=shutil.ignore_patterns("post", "log.*", "screen.*", "*.pyc", "__pycache__"))
            rewrite_input_parameters(workspace / directory / input_name, parameters)
            rewrite_input_seeds(workspace / directory / input_name, seeds)
        (workspace / "parameters.json").write_text(json.dumps({"parameters": parameters, "fixed_parameters": self.config.get("fixed_parameters", {}), "seeds": seeds}, indent=2), encoding="utf-8")
        return workspace

    def _append_csv(self, filename: str, fields, row: dict) -> None:
        path = self.run_root / filename
        with self._csv_lock:
            write_header = not path.exists()
            with path.open("a", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
                if write_header:
                    writer.writeheader()
                writer.writerow({key: row.get(key, "") for key in fields})

    def _record_status(self, **row) -> None:
        self._append_csv("status.csv", STATUS_FIELDS, row)

    def _run_input(self, workspace: Path, directory: str, input_name: str, stage: str) -> None:
        trial_dir = workspace / directory
        command = self._solver_invocation(trial_dir, input_name)
        log_path = workspace / f"{stage}.log"
        with log_path.open("w", encoding="utf-8") as log:
            process = subprocess.Popen(command, cwd=trial_dir, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1)
            assert process.stdout is not None
            with process.stdout:
                for line in process.stdout:
                    log.write(line)
                    log.flush()
            returncode = process.wait()
        if returncode != 0:
            raise RuntimeError(f"{stage} solver failed with exit code {returncode}; see {log_path}")

    @staticmethod
    def _mesh_bounds_for_snapshot(snapshot_path: Path) -> dict[str, tuple[float, float]]:
        """Get X/Z bounds from one legacy emitted STL or derive from the trial's source mesh and scale."""
        post_dir = Path(snapshot_path).parent
        candidates = []
        for path in post_dir.glob("box_*.stl"):
            match = re.fullmatch(r"box_(\d+)\.stl", path.name)
            if match:
                candidates.append((int(match.group(1)), path))
        if candidates:
            _, latest = max(candidates, key=lambda item: item[0])
            bounds = parse_ascii_stl_bounds(latest.read_text(encoding="utf-8", errors="replace"))
            return {axis: bounds[axis] for axis in ("x", "z")}

        test_dir = Path(snapshot_path).parents[1]
        deck_path = test_dir / "in.deslizamento"
        if not deck_path.is_file():
            raise RuntimeError(f"discharge deck missing beside measurement snapshot: {deck_path}")
        deck = deck_path.read_text(encoding="utf-8", errors="replace")
        box_fix = re.search(
            r"(?im)^\s*fix\s+box\s+\S+\s+mesh/surface\s+file\s+(\S+)\s+.*?\bscale\s+(\S+)",
            deck,
        )
        if not box_fix:
            raise RuntimeError(f"could not find scaled box mesh declaration in {deck_path}")

        definitions = {}
        for line in deck.splitlines():
            line = line.split("#", 1)[0]
            match = re.match(r"\s*variable\s+(\w+)\s+(?:equal|index)\s+(.+?)\s*$", line, re.IGNORECASE)
            if match:
                definitions[match.group(1)] = match.group(2).strip()

        def evaluate(expression: str, stack: tuple[str, ...] = ()) -> float:
            expression = re.sub(r"\$\((.*?)\)", r"\1", expression.strip())
            tree = ast.parse(expression, mode="eval")

            def visit(node) -> float:
                if isinstance(node, ast.Expression):
                    return visit(node.body)
                if isinstance(node, ast.Constant) and isinstance(node.value, (int, float)):
                    return float(node.value)
                if isinstance(node, ast.Name) and node.id.startswith("v_"):
                    name = node.id[2:]
                    if name in stack or name not in definitions:
                        raise ValueError(f"unknown or cyclic LIGGGHTS variable {node.id}")
                    return evaluate(definitions[name], stack + (name,))
                if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.UAdd, ast.USub)):
                    value = visit(node.operand)
                    return value if isinstance(node.op, ast.UAdd) else -value
                if isinstance(node, ast.BinOp) and isinstance(node.op, (ast.Add, ast.Sub, ast.Mult, ast.Div, ast.Pow)):
                    left, right = visit(node.left), visit(node.right)
                    if isinstance(node.op, ast.Add): return left + right
                    if isinstance(node.op, ast.Sub): return left - right
                    if isinstance(node.op, ast.Mult): return left * right
                    if isinstance(node.op, ast.Div): return left / right
                    return left ** right
                raise ValueError(f"unsupported LIGGGHTS scale expression: {expression}")

            return visit(tree)

        source_mesh = (test_dir / box_fix.group(1)).resolve()
        scale = evaluate(box_fix.group(2))
        if not source_mesh.is_file() or not math.isfinite(scale) or scale <= 0:
            raise RuntimeError(f"invalid source box mesh or scale: {source_mesh}, scale={scale}")
        bounds = parse_ascii_stl_bounds(source_mesh.read_text(encoding="utf-8", errors="replace"))
        bounds = {axis: (low * scale, high * scale) for axis, (low, high) in bounds.items()}
        return {axis: bounds[axis] for axis in ("x", "z")}

    @staticmethod
    def _nbr_ratio_from_deck(deck_path: Path) -> float:
        """Evaluate the trial-local LIGGGHTS nbrratio expression without executing deck code."""
        if not deck_path.is_file():
            raise RuntimeError(f"NBR deck missing: {deck_path}")
        definitions = {}
        for line in deck_path.read_text(encoding="utf-8", errors="replace").splitlines():
            line = line.split("#", 1)[0]
            match = re.match(r"\s*variable\s+(\w+)\s+(?:equal|index)\s+(.+?)\s*$", line, re.IGNORECASE)
            if match:
                definitions[match.group(1)] = match.group(2).strip()

        def evaluate(expression: str, stack: tuple[str, ...] = ()) -> float:
            expression = re.sub(r"\$\((.*?)\)", r"\1", expression.strip())
            tree = ast.parse(expression, mode="eval")

            def visit(node) -> float:
                if isinstance(node, ast.Expression):
                    return visit(node.body)
                if isinstance(node, ast.Constant) and isinstance(node.value, (int, float)):
                    return float(node.value)
                if isinstance(node, ast.Name) and node.id.startswith("v_"):
                    name = node.id[2:]
                    if name in stack or name not in definitions:
                        raise ValueError(f"unknown or cyclic LIGGGHTS variable {node.id}")
                    return evaluate(definitions[name], stack + (name,))
                if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.UAdd, ast.USub)):
                    value = visit(node.operand)
                    return value if isinstance(node.op, ast.UAdd) else -value
                if isinstance(node, ast.BinOp) and isinstance(node.op, (ast.Add, ast.Sub, ast.Mult, ast.Div, ast.Pow)):
                    left, right = visit(node.left), visit(node.right)
                    if isinstance(node.op, ast.Add): return left + right
                    if isinstance(node.op, ast.Sub): return left - right
                    if isinstance(node.op, ast.Mult): return left * right
                    if isinstance(node.op, ast.Div): return left / right
                    return left ** right
                raise ValueError(f"unsupported LIGGGHTS NBR expression: {expression}")

            return visit(tree)

        ratio = evaluate("v_nbrratio")
        if not math.isfinite(ratio) or ratio <= 0:
            raise ValueError(f"invalid nbrratio in {deck_path}: {ratio}")
        return ratio

    @staticmethod
    def _nbr_note_for_snapshot(snapshot_path: Path) -> str:
        """Record the NBR radius and ratio expressions from this trial's own input deck."""
        case_dir = snapshot_path.parents[1]
        deck_name = "in.repouso" if case_dir.name.lower().startswith("ang repouso") else "in.deslizamento"
        deck = case_dir / deck_name
        if not deck.is_file():
            return "NBR scaling: deck unavailable"
        found = {}
        for line in deck.read_text(encoding="utf-8", errors="replace").splitlines():
            match = re.match(r"\s*variable\s+(radiusnbr|nbrratio)\s+(?:equal|index)\s+([^#]+)", line, re.IGNORECASE)
            if match:
                found[match.group(1).lower()] = match.group(2).strip()
        return "NBR deck values: " + "; ".join(f"{name}={found[name]}" for name in ("radiusnbr", "nbrratio") if name in found) if found else "NBR scaling values not found in deck"

    def _measure_snapshot(self, path: Path, kind: str) -> dict:
        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError(f"measurement artifact missing: {path}")
        snapshot_text = path.read_text(encoding="utf-8")
        points = parse_vtk_points(snapshot_text)
        radii = None
        trial_dir = path.parents[2]
        trial_label = trial_dir.name

        def save_diagnostic(reason: str, measurement_diagnostic=None) -> None:
            axis = "z" if kind == "discharge" else "radial"
            output = trial_dir / ("medicao_descarga_diagnostico.svg" if kind == "discharge" else "medicao_repouso_diagnostico.svg")
            try:
                if kind == "repose":
                    output = trial_dir / "medicao_repouso_diagnostico.svg"
                    if measurement_diagnostic is not None:
                        write_repose_svg(output, measurement_diagnostic, trial_label, "medicao rejeitada; perfil diagnostico", points=points, historical_preview=False)
                    else:
                        write_repose_failure_diagnostic_svg(output, points, trial_label, reason, particle_radii=radii)
                else:
                    output = trial_dir / "medicao_descarga.svg"
                    ratio = self._nbr_ratio_from_deck(path.parents[1] / "in.deslizamento")
                    bounds = self._mesh_bounds_for_snapshot(path)
                    write_discharge_failure_svg(
                        output, points, trial_label, reason,
                        horizontal_axis="z", nbrratio=ratio, mesh_bounds=bounds,
                        particle_radii=radii, minimum_height_particle_diameters=2.0,
                    )
                png_output = output.with_suffix(".png")
                write_svg_png(output, png_output)
            except Exception as render_exc:
                note = output.with_suffix(output.suffix + ".txt")
                note.write_text(f"Measurement diagnostic: {reason}\nVisual generation failed: {type(render_exc).__name__}: {render_exc}\n", encoding="utf-8")

        try:
            if kind == "discharge":
                radii = parse_vtk_scalar_field(snapshot_text, "radius")
                mesh_bounds = self._mesh_bounds_for_snapshot(path)
                nbr_ratio = self._nbr_ratio_from_deck(path.parents[1] / "in.deslizamento")
                result = measure_v_sliding_angle_three_cuts(
                    points,
                    radii,
                    horizontal_axis="z",
                    mesh_bounds=mesh_bounds,
                    nbrratio=nbr_ratio,
                )

            else:
                radii = parse_vtk_scalar_field(snapshot_text, "radius")
                result = measure_pile_angle(points, particle_radii=radii)
        except GeometryDiagnostic as exc:
            save_diagnostic(str(exc))
            raise
        except UnmeasurableAngle as exc:
            if kind == "discharge":
                diagnostic = MeasurementRejected(
                    f"{kind} angle unavailable: no reliable free-face fit: {exc}"
                )
            else:
                diagnostic = GeometryDiagnostic(
                    f"{kind} angle unavailable: profile is not a reliable repose face; cohesive bridge/arching suspected, not confirmed: {exc}"
                )
            save_diagnostic(str(diagnostic), getattr(exc, "diagnostic_result", None))
            raise diagnostic from exc
        except ValueError as exc:
            diagnostic = MeasurementRejected(f"{kind} measurement rejected: {exc}")
            save_diagnostic(str(diagnostic))
            raise diagnostic from exc
        nbr_note = self._nbr_note_for_snapshot(path)
        if kind == "discharge":
            svg_path = trial_dir / "medicao_descarga.svg"
            png_path = trial_dir / "medicao_descarga.png"
            write_discharge_svg(svg_path, points, result, trial_label, nbr_note, particle_radii=radii)
            write_svg_png(svg_path, png_path)
        else:
            svg_path = trial_dir / "medicao_repouso_oito_direcoes.svg"
            png_path = trial_dir / "medicao_repouso_oito_direcoes.png"
            write_repose_svg(svg_path, result, trial_label, nbr_note, points=points, historical_preview=False)
            write_svg_png(svg_path, png_path)
        return result

    def execute_trial(self, trial_id: int, parameters: dict[str, float], phase: str = "doe", seeds: dict[str, int] | None = None) -> dict:
        started = _utc_now()
        seeds = seeds or self._trial_seeds(trial_id)
        workspace = self.prepare_trial(trial_id, parameters, seeds)
        base = {"trial_id": trial_id, "phase": phase, "workspace": str(workspace), "started_at": started}
        self._record_status(**base, status="running", stage="prepare", message="workspace created", finished_at="")
        try:
            self._run_input(workspace, "Ang Deslizamento", "in.deslizamento", "discharge")
            # Independent cases: preserve both snapshots before post-processing.
            self._run_input(workspace, "Ang Repouso", "in.repouso", "repose")
            discharge_measurement = self._measure_snapshot(workspace / "Ang Deslizamento/post/deslizamento_final.vtk", "discharge")
            v_angle = discharge_measurement["mean_angle_deg"]
            repose_measurement = self._measure_snapshot(workspace / "Ang Repouso/post/repouso_10s.vtk", "repose")
            repose_angle = repose_measurement["mean_angle_deg"]
            repose_cut_results = repose_pair_result_fields(repose_measurement)
            targets = {"discharge": self.config["target_discharge_deg"], "repose": self.config["target_repose_deg"]}
            values = {"discharge": v_angle, "repose": repose_angle}
            errors = {name: abs(values[name] - targets[name]) for name in targets}
            discharge_cut_results = discharge_cut_result_fields(discharge_measurement)
            result = {**self.config.get("fixed_parameters", {}), **parameters, "trial_id": trial_id, "phase": phase, "status": "complete", "workspace": str(workspace), "distribution_seed": seeds["distribution"], "insertion_seed": seeds["insertion"], "discharge_left_deg": discharge_measurement["left_angle_deg"], "discharge_right_deg": discharge_measurement["right_angle_deg"], "discharge_asymmetry_deg": discharge_measurement["asymmetry_deg"], **discharge_cut_results, **repose_cut_results, **{f"{name}_deg": values[name] for name in values}, **{f"{name}_error": errors[name] for name in errors}}
            self._append_csv("results.csv", RESULT_FIELDS, result)
            self._record_status(**base, status="complete", stage="measure", message="all artifacts measured", finished_at=_utc_now())
            return result
        except GeometryDiagnostic as exc:
            self._record_status(**base, status="diagnostic", stage="measurement", message=str(exc), finished_at=_utc_now())
            raise
        except MeasurementRejected as exc:
            self._record_status(**base, status="rejected", stage="measurement", message=str(exc), finished_at=_utc_now())
            raise
        except Exception as exc:
            self._record_status(**base, status="failed", stage="failed", message=str(exc), finished_at=_utc_now())
            raise RuntimeError(f"trial {trial_id} stopped after failed measurement or solver: {exc}") from exc

    def score_completed_results(self) -> list[dict]:
        path = self.run_root / "results.csv"
        if not path.is_file():
            return []
        with path.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        if not rows:
            return []
        metrics = ("discharge", "repose")
        maxima = {metric: max(float(row[f"{metric}_error"]) for row in rows) for metric in metrics}
        for row in rows:
            errors = {metric: float(row[f"{metric}_error"]) for metric in metrics}
            row["score"] = weighted_normalized_score(errors, maxima, self.config["weights"])
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=RESULT_FIELDS)
            writer.writeheader()
            writer.writerows(rows)
        return rows

    def _trial_seeds(self, trial_id: int) -> dict[str, int]:
        rng = random.Random(int(self.config.get("random_seed", 1)) + trial_id * 1009)
        distribution = next_liggghts_prime(rng.randrange(10_001, 899_000_000))
        insertion = next_liggghts_prime(rng.randrange(10_001, 899_000_000))
        if insertion == distribution:
            insertion = next_liggghts_prime(insertion + 2)
        return {"distribution": distribution, "insertion": insertion}

    def _point_key(self, parameters: dict[str, float]) -> tuple[tuple[str, float], ...]:
        """Cria uma chave estável para retomar pontos DOE sem repetir parâmetros já encerrados."""
        return tuple((name, round(float(parameters[name]), 12)) for name in self.parameter_names)

    def _is_within_active_bounds(self, parameters: dict[str, float]) -> bool:
        """Aceita somente pontos dentro dos limites ativos, inclusive nas extremidades."""
        return all(float(spec["min"]) <= float(parameters[name]) <= float(spec["max"]) for name, spec in self.config["parameters"].items())

    def _calibrate_cohesion_ceiling(self) -> None:
        """Bracket a non-measurable cohesion value at the same sf11/rf11 pair before refinement."""
        if self._cohesion_tracker is None:
            return
        if self._boundary_search_report and self._boundary_search_report.get("status") == "complete":
            return
        prior = self.run_root / "cohesion_boundary_search.json"
        if prior.is_file():
            saved = json.loads(prior.read_text(encoding="utf-8"))
            if saved.get("status") == "complete":
                self._boundary_search_report = saved
                ceiling = saved.get("recommended_ceiling_CED11")
                if ceiling is not None:
                    self.config["parameters"]["CED11"]["max"] = float(ceiling)
                self._write_cohesion_report("boundary_search_complete")
                return
        bad = next((row for row in sorted(self._cohesion_tracker.problematic, key=lambda row: row["CED11"]) if "sf11" in row["parameters"] and "rf11" in row["parameters"]), None)
        if bad is None:
            self._boundary_search_report = {"status": "not_needed", "detail": "No cohesion-associated non-measurable geometry was observed."}
            self._write_cohesion_report("no_problem_observed")
            return
        p = bad["parameters"]
        lower = float(self.config["parameters"]["CED11"]["min"])
        upper = float(bad["CED11"])
        pair = {name: p[name] for name in ("sf11", "rf11")}
        existing = [row for row in self._cohesion_tracker.valid if all(abs(row["parameters"].get(name, float("inf")) - value) <= 1e-10 for name, value in pair.items()) and lower <= row["CED11"] < upper]
        if existing:
            low_row = max(existing, key=lambda row: row["CED11"])
            low_value = float(low_row["CED11"])
        else:
            baseline = {**p, "CED11": lower}
            completed = self._run_points([baseline], "cohesion_boundary", self._next_trial_id)
            if not completed:
                self._boundary_search_report = {"status": "no_valid_lower_bracket", "problematic_CED11": upper, "baseline_CED11": lower, "parameters": pair}
                prior.write_text(json.dumps(self._boundary_search_report, indent=2), encoding="utf-8")
                self._write_cohesion_report("boundary_search_blocked")
                return
            low_value = lower
        search = CohesionLimitTracker(lower, upper, max_probes=int(self.config.get("cohesion_guard", {}).get("max_boundary_probes", 8)))
        search.observe_valid(low_value, 0, pair)
        search.observe_problematic(upper, int(bad["trial_id"]), pair, bad.get("reason", "geometry not measurable"))
        probes = []
        while (candidate := search.next_probe()) is not None:
            point = {**p, "CED11": candidate}
            before_bad_count = len(self._cohesion_tracker.problematic)
            completed = self._run_points([point], "cohesion_boundary", self._next_trial_id)
            if completed:
                search.observe_valid(candidate, int(completed[0]["trial_id"]), point)
                probes.append({"CED11": candidate, "status": "measurable", "trial_id": completed[0]["trial_id"]})
            elif len(self._cohesion_tracker.problematic) > before_bad_count:
                latest = self._cohesion_tracker.problematic[-1]
                search.observe_problematic(candidate, int(latest["trial_id"]), point, latest["reason"])
                probes.append({"CED11": candidate, "status": "not_measurable", "trial_id": latest["trial_id"]})
            else:
                probes.append({"CED11": candidate, "status": "rejected_or_solver_failure"})
                break
        ceiling = search.recommended_ceiling
        self.config["parameters"]["CED11"]["max"] = ceiling
        self._boundary_search_report = {
            "status": "complete",
            "parameters": pair,
            "initial_problematic_CED11": upper,
            "highest_measurable_CED11": search.highest_valid_below_problem,
            "lowest_problematic_CED11": search.lowest_problematic,
            "highest_problematic_CED11": search.highest_problematic,
            "recommended_ceiling_CED11": ceiling,
            "remaining_interval_CED11": search.lowest_problematic - ceiling if search.lowest_problematic is not None else None,
            "probes": probes,
            "interpretation": "The ceiling is the highest tested measurable value below the lowest tested problematic value for this sf11/rf11 pair; it is an empirical guard, not a universal material law.",
        }
        prior.write_text(json.dumps(self._boundary_search_report, indent=2, sort_keys=True), encoding="utf-8")
        self._write_cohesion_report("boundary_search_complete")

    def _doe_backfill_batch_size(self, completed_count: int) -> int:
        """Mantém um lote cheio de contingência até atingir o mínimo de DOE completo."""
        minimum = int(self.config.get("minimum_doe_complete", self.config["campaign"]["doe_trials"]))
        if completed_count >= minimum:
            return 0
        return max(minimum - completed_count, int(self.config.get("doe_contingency_batch_size", self.config.get("parallel_trials", 1))))

    def _resume_doe_state(self) -> tuple[set[tuple[tuple[str, float], ...]], dict[tuple[tuple[str, float], ...], int], int]:
        """Lê resultados e rejeições persistidos para retomar apenas pontos DOE pendentes."""
        completed = set()
        rejected_counts = {}
        maximum_trial_id = 0
        results_path = self.run_root / "results.csv"
        if results_path.is_file():
            with results_path.open(newline="", encoding="utf-8") as handle:
                for row in csv.DictReader(handle):
                    if row.get("phase") == "doe" and row.get("status") == "complete":
                        completed.add(self._point_key(row))
                        maximum_trial_id = max(maximum_trial_id, int(row["trial_id"]))
        status_path = self.run_root / "status.csv"
        if status_path.is_file():
            with status_path.open(newline="", encoding="utf-8") as handle:
                for row in csv.DictReader(handle):
                    maximum_trial_id = max(maximum_trial_id, int(row["trial_id"]))
                    if row.get("phase") != "doe" or row.get("status") != "rejected":
                        continue
                    parameters_path = Path(row["workspace"]) / "parameters.json"
                    if not parameters_path.is_file():
                        continue
                    parameters = json.loads(parameters_path.read_text(encoding="utf-8")).get("parameters", {})
                    key = self._point_key(parameters)
                    rejected_counts[key] = rejected_counts.get(key, 0) + 1
        diagnostics_path = self.run_root / "diagnostics.csv"
        if diagnostics_path.is_file():
            with diagnostics_path.open(newline="", encoding="utf-8") as handle:
                for row in csv.DictReader(handle):
                    if row.get("phase") != "doe" or not row.get("parameters_json"):
                        continue
                    parameters = json.loads(row["parameters_json"])
                    completed.add(self._point_key(parameters))
                    maximum_trial_id = max(maximum_trial_id, int(row["trial_id"]))
        return completed, rejected_counts, maximum_trial_id + 1

    def refresh_live_state(self, run_root: Path, active: dict[int, dict]) -> None:
        refresh_live_state(run_root, active)

    def _run_points(self, points: list[dict], phase: str, first_trial_id: int, initial_retries: dict | None = None) -> list[dict]:
        workers = min(len(points), int(self.config.get("parallel_trials", 1)))
        if workers < 1:
            raise ValueError("parallel_trials must be positive")
        results = []
        point_iterator = iter(points)
        next_trial_id = first_trial_id
        retry_limit = int(self.config.get("measurement_retry_limit", 2))
        retries = dict(initial_retries or {})
        active = {}
        dashboard_interval = 10
        maximum = self.config["campaign"]["max_trials"]
        with ThreadPoolExecutor(max_workers=workers) as executor:
            pending = {}
            def submit(point):
                nonlocal next_trial_id
                trial_id = next_trial_id; next_trial_id += 1
                active[trial_id] = {"phase": phase, "status": "preparing", "started_at": time.time(), "simulation_time_s": 0.0}
                pending[executor.submit(self.execute_trial, trial_id, point, phase, self._trial_seeds(trial_id))] = (trial_id, point)

            for _ in range(workers):
                try:
                    submit(next(point_iterator))
                except StopIteration:
                    break
            while pending:
                self.refresh_live_state(self.run_root, active)
                write_live_status_csv(self.run_root, active)
                render_dashboard(active, first_trial_id - 1 + len(results), maximum, interval=dashboard_interval)
                completed, _ = wait(pending, timeout=dashboard_interval, return_when=FIRST_COMPLETED)
                for future in completed:
                    trial_id, point = pending.pop(future)
                    try:
                        result = future.result()
                        results.append(result)
                        if "CED11" in point:
                            self._observe_cohesion(float(point["CED11"]), trial_id, point, measurable=True)
                        active.pop(trial_id, None)
                    except GeometryDiagnostic as exc:
                        active.pop(trial_id, None)
                        if "CED11" in point:
                            self._observe_cohesion(float(point["CED11"]), trial_id, point, measurable=False, reason=str(exc))
                        self._append_csv("diagnostics.csv", DIAGNOSTIC_FIELDS, {"trial_id": trial_id, "phase": phase, "status": "geometry_unmeasurable", "reason": str(exc), "parameters_json": json.dumps(point, sort_keys=True), "workspace": str(self.run_root / f"trial_{trial_id:04d}")})
                        try:
                            submit(next(point_iterator))
                        except StopIteration:
                            pass
                        continue
                    except MeasurementRejected:
                        active.pop(trial_id, None)
                        key = tuple(sorted((name, float(value)) for name, value in point.items()))
                        retries[key] = retries.get(key, 0) + 1
                        if retries[key] > retry_limit:
                            write_live_status_csv(self.run_root, active)
                            try:
                                submit(next(point_iterator))
                            except StopIteration:
                                pass
                            continue
                        submit(point)
                        continue
                    except Exception as exc:
                        active.pop(trial_id, None)
                        self._record_status(trial_id=trial_id, phase=phase, status="failed", stage="failed", message=str(exc), workspace=str(self.run_root / f"trial_{trial_id:04d}"), started_at="", finished_at=_utc_now())
                        try:
                            submit(next(point_iterator))
                        except StopIteration:
                            pass
                        continue
                    try:
                        submit(next(point_iterator))
                    except StopIteration:
                        continue
            self.refresh_live_state(self.run_root, active)
            write_live_status_csv(self.run_root, active)
            render_dashboard(active, first_trial_id - 1 + len(results), maximum, interval=dashboard_interval)
        self._next_trial_id = next_trial_id
        return results

    def _quadratic_features(self, point: dict[str, float]) -> list[float]:
        values = [(float(point[name]) - self.config["parameters"][name]["min"]) / (self.config["parameters"][name]["max"] - self.config["parameters"][name]["min"]) for name in self.parameter_names]
        return [1.0, *values, *(value * value for value in values), *(values[i] * values[j] for i in range(len(values)) for j in range(i + 1, len(values)))]

    @staticmethod
    def _solve_linear_system(matrix: list[list[float]], vector: list[float]) -> list[float]:
        augmented = [row[:] + [value] for row, value in zip(matrix, vector)]
        for column in range(len(vector)):
            pivot = max(range(column, len(vector)), key=lambda row: abs(augmented[row][column]))
            if abs(augmented[pivot][column]) < 1e-12:
                raise ValueError("quadratic surrogate is rank deficient")
            augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
            divisor = augmented[column][column]
            augmented[column] = [value / divisor for value in augmented[column]]
            for row in range(len(vector)):
                if row != column:
                    factor = augmented[row][column]
                    augmented[row] = [value - factor * base for value, base in zip(augmented[row], augmented[column])]
        return [row[-1] for row in augmented]

    def nominate_refinement_candidates(self, completed_rows: list[dict]) -> list[dict]:
        rows = [row for row in completed_rows if row.get("status") == "complete" and row.get("score", "") != "" and row.get("phase") != "cohesion_boundary" and self._is_within_active_bounds(row)]
        initial = {name: self.config["parameters"][name]["initial"] for name in self.parameter_names}
        terms = len(self._quadratic_features(initial))
        if len(rows) < terms:
            raise ValueError(f"quadratic surrogate needs at least {terms} scored trials")
        features = [self._quadratic_features(row) for row in rows]
        scores = [float(row["score"]) for row in rows]
        normal = [[sum(row[i] * row[j] for row in features) + (1e-10 if i == j else 0.0) for j in range(terms)] for i in range(terms)]
        coefficients = self._solve_linear_system(normal, [sum(row[i] * score for row, score in zip(features, scores)) for i in range(terms)])
        bounds = {name: (spec["min"], spec["max"]) for name, spec in self.config["parameters"].items()}
        pool = latin_hypercube(int(self.config.get("model_candidate_pool", 1500)), bounds, int(self.config.get("random_seed", 1)) + 7919)
        observed = {tuple(round(float(row[name]), 12) for name in self.parameter_names) for row in rows}
        selected = []
        for _, point in sorted(((sum(value * coefficient for value, coefficient in zip(self._quadratic_features(point), coefficients)), point) for point in pool), key=lambda item: item[0]):
            key = tuple(round(float(point[name]), 12) for name in self.parameter_names)
            if key not in observed:
                selected.append(point)
                observed.add(key)
            if len(selected) == self.config["campaign"]["refinement_trials"]:
                return selected
        raise RuntimeError("could not nominate enough distinct refinement candidates")

    def run_doe(self) -> list[dict]:
        scored = self.score_completed_results()
        minimum = int(self.config.get("minimum_doe_complete", self.config["campaign"]["doe_trials"]))
        if sum(row.get("phase") == "doe" and row.get("status") == "complete" for row in scored) >= minimum:
            _, _, self._next_trial_id = self._resume_doe_state()
            self._calibrate_cohesion_ceiling()
            return scored
        bounds = {name: (spec["min"], spec["max"]) for name, spec in self.config["parameters"].items()}
        manifest_path = self.run_root / "doe_manifest_original.json"
        if manifest_path.is_file():
            points = json.loads(manifest_path.read_text(encoding="utf-8"))["points"]
        else:
            points = latin_hypercube(self.config["campaign"]["doe_trials"], bounds, self.config.get("random_seed", 1))
        completed, rejected_counts, next_trial_id = self._resume_doe_state()
        exhaustion_count = int(self.config.get("measurement_retry_limit", 2)) + 1
        pending = [point for point in points if self._is_within_active_bounds(point) and self._point_key(point) not in completed and rejected_counts.get(self._point_key(point), 0) < exhaustion_count]
        needed = sum(1 for point in points if not self._is_within_active_bounds(point) and self._point_key(point) not in completed)
        observed = {self._point_key(point) for point in points} | completed
        replacements = []
        if needed > 0:
            for point in latin_hypercube(max(50, needed * 10), bounds, int(self.config.get("random_seed", 1)) + 104729):
                if self._point_key(point) not in observed:
                    replacements.append(point)
                    observed.add(self._point_key(point))
                if len(replacements) == needed:
                    break
        if len(replacements) != needed:
            raise RuntimeError("could not generate enough bounded DOE replacement points")
        pending.extend(replacements)
        initial_retries = {self._point_key(point): rejected_counts.get(self._point_key(point), 0) - 1 for point in pending if rejected_counts.get(self._point_key(point), 0)}
        self._run_points(pending, "doe", next_trial_id, initial_retries)
        scored = self.score_completed_results()
        while True:
            complete_doe = [row for row in scored if row.get("phase") == "doe" and row.get("status") == "complete"]
            batch_size = self._doe_backfill_batch_size(len(complete_doe))
            if batch_size == 0:
                self._calibrate_cohesion_ceiling()
                return scored
            completed, rejected_counts, next_trial_id = self._resume_doe_state()
            observed = {self._point_key(point) for point in points} | completed | set(rejected_counts)
            backfill = []
            seed = int(self.config.get("random_seed", 1)) + 209759 + next_trial_id
            for point in latin_hypercube(max(50, batch_size * 10), bounds, seed):
                if self._point_key(point) not in observed:
                    backfill.append(point)
                    observed.add(self._point_key(point))
                if len(backfill) == batch_size:
                    break
            if len(backfill) != batch_size:
                raise RuntimeError("could not generate enough DOE contingency points")
            self._run_points(backfill, "doe", next_trial_id)
            scored = self.score_completed_results()

    def run_campaign(self) -> list[dict]:
        doe = self.run_doe()
        self._run_points(self.nominate_refinement_candidates(doe), "refinement", self._next_trial_id)
        scored = self.score_completed_results()
        best = min(scored, key=lambda row: (float(row["score"]), int(row["trial_id"])))
        validation = [{name: float(best[name]) for name in self.parameter_names} for _ in range(self.config["campaign"]["validation_trials"])]
        self._run_points(validation, "validation", self._next_trial_id)
        final_results = self.score_completed_results()
        self._write_cohesion_report("campaign_complete")
        return final_results

    def smoke(self, trials: int = 3) -> list[dict]:
        if trials < 1:
            raise ValueError("smoke trials must be positive")
        bounds = {name: (spec["min"], spec["max"]) for name, spec in self.config["parameters"].items()}
        points = latin_hypercube(trials, bounds, self.config.get("random_seed", 1))
        workers = min(trials, int(self.config.get("parallel_trials", trials)))
        results = []
        with ThreadPoolExecutor(max_workers=workers) as executor:
            pending = {}
            point_iterator = iter(enumerate(points, start=1))
            for _ in range(workers):
                try:
                    trial_id, point = next(point_iterator)
                except StopIteration:
                    break
                pending[executor.submit(self.execute_trial, trial_id, point)] = trial_id
            while pending:
                completed, _ = wait(pending, return_when=FIRST_COMPLETED)
                for future in completed:
                    pending.pop(future)
                    try:
                        results.append(future.result())
                    except Exception:
                        for queued in pending:
                            queued.cancel()
                        raise
                    try:
                        trial_id, point = next(point_iterator)
                    except StopIteration:
                        continue
                    pending[executor.submit(self.execute_trial, trial_id, point)] = trial_id
        return self.score_completed_results()


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, allow_abbrev=False)
    parser.add_argument("--project-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--run-root", type=Path)
    parser.add_argument("--preflight", action="store_true")
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--run", action="store_true")
    parser.add_argument("--smoke-trials", type=int, default=3)
    parser.add_argument("--parallel-trials", type=int, help="override concurrency without modifying the saved config")
    parser.add_argument("--auto-workers", action="store_true", help="use the smaller of configured max workers and CPUs visible to Python/WSL")
    args = parser.parse_args(argv)
    if args.auto_workers and args.parallel_trials is not None:
        parser.error("--auto-workers cannot be combined with --parallel-trials")
    workers_override = args.parallel_trials
    if args.auto_workers:
        configured_workers = int(load_configuration(args.project_root / "configuracao_otimizacao.txt").get("parallel_trials", 1))
        workers_override = min(configured_workers, max(1, os.cpu_count() or 1))
    controller = CalibrationController(args.project_root, args.run_root, parallel_trials_override=workers_override)
    preflight = controller.preflight()
    if args.preflight:
        print(json.dumps(preflight, indent=2))
        return 0
    if args.smoke:
        print(json.dumps(controller.smoke(args.smoke_trials), indent=2))
        return 0
    if args.run:
        results = controller.run_campaign()
        print(json.dumps({
            "status": "complete",
            "doe_trials": controller.config["campaign"]["doe_trials"],
            "refinement_trials": controller.config["campaign"]["refinement_trials"],
            "validation_trials": controller.config["campaign"]["validation_trials"],
            "completed_results": len(results),
            "cohesion_limits_report": str(controller.run_root / "cohesion_limits.json"),
            "recommended_ceiling_CED11": (controller._boundary_search_report or {}).get("recommended_ceiling_CED11"),
            "lowest_problematic_CED11": (controller._boundary_search_report or {}).get("lowest_problematic_CED11"),
            "highest_problematic_CED11": (controller._boundary_search_report or {}).get("highest_problematic_CED11"),
        }, indent=2))
        return 0
    parser.error("choose --preflight, --smoke, or --run")


if __name__ == "__main__":
    # This project uses the supplied sequential opening protocol. The class
    # above remains a tested compatibility/scheduler base, not the CLI entry.
    from optimizer.opening_controller import main as opening_main
    sys.exit(opening_main())
