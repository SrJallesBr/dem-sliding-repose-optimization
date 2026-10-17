import csv
import io
import json
import shutil
import sys
import tempfile
import threading
import time
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from unittest import mock

PROJECT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT))

from optimizer.calibration_controller import CalibrationController, GeometryDiagnostic, MeasurementRejected, is_valid_liggghts_seed, main, render_dashboard, write_live_status_csv


class ParallelCampaignTests(unittest.TestCase):
    def test_preflight_rejects_more_workers_than_visible_cpus(self):
        root = Path(tempfile.mkdtemp())
        config_path = root / "configuracao_otimizacao.txt"
        shutil.copy2(PROJECT / "configuracao_otimizacao.txt", config_path)
        config = json.loads(config_path.read_text(encoding="utf-8"))
        config["parallel_trials"] = 20
        config_path.write_text(json.dumps(config), encoding="utf-8")
        shutil.copytree(PROJECT / "Ang Deslizamento", root / "Ang Deslizamento")
        shutil.copytree(PROJECT / "Ang Repouso", root / "Ang Repouso")
        controller = CalibrationController(root, run_root=root / "runs", solver_command=["true"])
        with mock.patch("optimizer.calibration_controller.os.cpu_count", return_value=12):
            with self.assertRaisesRegex(RuntimeError, "parallel_trials=20 exceeds CPUs visible"):
                controller.preflight()

    def test_auto_workers_caps_configured_ceiling_to_visible_cpus_without_editing_config(self):
        root = Path(tempfile.mkdtemp())
        config_path = root / "configuracao_otimizacao.txt"
        shutil.copy2(PROJECT / "configuracao_otimizacao.txt", config_path)
        config = json.loads(config_path.read_text(encoding="utf-8"))
        config["parallel_trials"] = 20
        config_path.write_text(json.dumps(config), encoding="utf-8")
        shutil.copytree(PROJECT / "Ang Deslizamento", root / "Ang Deslizamento")
        shutil.copytree(PROJECT / "Ang Repouso", root / "Ang Repouso")
        output = io.StringIO()
        with mock.patch("optimizer.calibration_controller.os.cpu_count", return_value=12), redirect_stdout(output):
            exit_code = main(["--project-root", str(root), "--run-root", str(root / "runs"), "--preflight", "--auto-workers"])
        result = json.loads(output.getvalue())
        self.assertEqual(exit_code, 0)
        self.assertEqual(result["parallel_trials"], 12)
        self.assertEqual(result["available_cpus"], 12)
        self.assertEqual(json.loads(config_path.read_text(encoding="utf-8"))["parallel_trials"], 20)

    def test_preflight_rejects_processor_grid_mismatch_for_single_rank(self):
        root = Path(tempfile.mkdtemp())
        shutil.copy2(PROJECT / "configuracao_otimizacao.txt", root / "configuracao_otimizacao.txt")
        shutil.copytree(PROJECT / "Ang Deslizamento", root / "Ang Deslizamento")
        shutil.copytree(PROJECT / "Ang Repouso", root / "Ang Repouso")
        input_path = root / "Ang Deslizamento" / "in.deslizamento"
        input_path.write_text(input_path.read_text(encoding="utf-8").replace("processors    1 1 1", "processors    2 1 1", 1), encoding="utf-8")
        controller = CalibrationController(root, run_root=root / "runs", solver_command=["true"])
        with self.assertRaisesRegex(RuntimeError, "processors grid in Ang Deslizamento/in.deslizamento must multiply to 1"):
            controller.preflight()

    def test_trial_seeds_are_valid_liggghts_primes(self):
        root = Path(tempfile.mkdtemp())
        config = {"parameters": {"sf11": {"min": 0.3, "max": 1.1, "initial": 0.95}}, "mpi_processes_per_trial": 1, "parallel_trials": 1, "random_seed": 86028157}
        (root / "configuracao_otimizacao.txt").write_text(json.dumps(config), encoding="utf-8")
        controller = CalibrationController(root, run_root=root / "runs", solver_command=["true"])
        for trial_id in range(1, 66):
            seeds = controller._trial_seeds(trial_id)
            self.assertNotEqual(seeds["distribution"], seeds["insertion"])
            self.assertTrue(is_valid_liggghts_seed(seeds["distribution"]))
            self.assertTrue(is_valid_liggghts_seed(seeds["insertion"]))

    def test_live_dashboard_renders_table_and_writes_status_csv(self):
        root = Path(tempfile.mkdtemp())
        active = {
            7: {"phase": "doe", "status": "repose_drawdown", "started_at": 100.0, "simulation_time_s": 13.2},
        }
        write_live_status_csv(root, active, now=165.0)
        with (root / "status_execucao.csv").open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(rows[0]["trial"], "trial_0007")
        self.assertEqual(rows[0]["simulation_time_s"], "13.2")
        output = io.StringIO()
        with redirect_stdout(output):
            render_dashboard(active, completed=6, maximum=65, now=165.0, interval=10)
        table = output.getvalue()
        self.assertIn("PAINEL DA OTIMIZACAO DEM", table)
        self.assertIn("trial_0007", table)
        self.assertIn("13,2s", table)
        self.assertIn("00:01:05", table)

    def test_run_input_keeps_solver_output_in_log_without_printing_it(self):
        root = Path(tempfile.mkdtemp())
        config = {"parameters": {"sf11": {"min": 0.2, "max": 2.0}}}
        (root / "configuracao_otimizacao.txt").write_text(json.dumps(config), encoding="utf-8")
        controller = CalibrationController(root, run_root=root / "runs", solver_command=[sys.executable, "-c", "print('solver chatter')"])
        workspace = root / "trial_0001"
        (workspace / "case").mkdir(parents=True)
        output = io.StringIO()
        with redirect_stdout(output):
            controller._run_input(workspace, "case", "in.deck", "discharge")
        self.assertEqual(output.getvalue(), "")
        self.assertIn("solver chatter", (workspace / "discharge.log").read_text(encoding="utf-8"))

    def test_run_points_limits_concurrency_to_configured_parallel_trials(self):
        root = Path(tempfile.mkdtemp())
        config = {
            "parameters": {"sf11": {"min": 0.3, "max": 1.1, "initial": 0.95}},
            "mpi_processes_per_trial": 1,
            "parallel_trials": 12,
            "random_seed": 7,
            "campaign": {"max_trials": 25},
        }
        (root / "configuracao_otimizacao.txt").write_text(json.dumps(config), encoding="utf-8")
        controller = CalibrationController(root, run_root=root / "runs", solver_command=["true"])
        active = 0
        maximum_active = 0
        lock = threading.Lock()

        def execute_trial(trial_id, parameters, phase="doe", seeds=None):
            nonlocal active, maximum_active
            with lock:
                active += 1
                maximum_active = max(maximum_active, active)
            time.sleep(0.12)
            with lock:
                active -= 1
            return {"trial_id": trial_id, "phase": phase, "status": "complete", **parameters}

        controller.execute_trial = execute_trial
        points = [{"sf11": 0.3 + index * 0.01} for index in range(25)]
        controller._run_points(points, "doe", 1)
        self.assertEqual(maximum_active, 12)

    def test_run_points_refills_a_freed_worker_without_waiting_for_the_cohort(self):
        root = Path(tempfile.mkdtemp())
        config = {
            "parameters": {"sf11": {"min": 0.2, "max": 2.0}},
            "mpi_processes_per_trial": 1,
            "parallel_trials": 2,
            "random_seed": 7,
            "campaign": {"max_trials": 3},
        }
        (root / "configuracao_otimizacao.txt").write_text(json.dumps(config), encoding="utf-8")
        controller = CalibrationController(root, run_root=root / "runs", solver_command=["true"])
        events = {}
        lock = threading.Lock()

        def execute_trial(trial_id, parameters, phase="doe", seeds=None):
            started = time.monotonic()
            with lock:
                events[trial_id] = {"started": started}
            time.sleep(0.04 if trial_id == 1 else 0.30)
            with lock:
                events[trial_id]["finished"] = time.monotonic()
            return {"trial_id": trial_id, "phase": phase, "status": "complete", **parameters}

        controller.execute_trial = execute_trial
        controller._run_points([{"sf11": 0.4}, {"sf11": 0.8}, {"sf11": 1.2}], "doe", 1)
        self.assertLess(events[3]["started"], events[2]["finished"])

    def test_live_status_removes_completed_trials(self):
        root = Path(tempfile.mkdtemp())
        config = {
            "parameters": {"sf11": {"min": 0.3, "max": 1.1, "initial": 0.95}},
            "mpi_processes_per_trial": 1,
            "parallel_trials": 1,
            "random_seed": 7,
            "campaign": {"max_trials": 2},
        }
        (root / "configuracao_otimizacao.txt").write_text(json.dumps(config), encoding="utf-8")
        controller = CalibrationController(root, run_root=root / "runs", solver_command=["true"])
        controller.execute_trial = lambda trial_id, parameters, phase="doe", seeds=None: {"trial_id": trial_id, "phase": phase, "status": "complete", **parameters}
        controller._run_points([{"sf11": 0.3}, {"sf11": 0.4}], "doe", 1)
        with (root / "runs" / "status_execucao.csv").open(newline="", encoding="utf-8") as handle:
            self.assertEqual(list(csv.DictReader(handle)), [])

    def test_run_doe_resumes_without_repeating_completed_or_exhausted_points(self):
        root = Path(tempfile.mkdtemp())
        config = {
            "parameters": {"sf11": {"min": 0.3, "max": 1.1, "initial": 0.95}},
            "mpi_processes_per_trial": 1,
            "parallel_trials": 1,
            "random_seed": 7,
            "measurement_retry_limit": 2,
            "minimum_doe_complete": 3,
            "campaign": {"max_trials": 3, "doe_trials": 3},
        }
        (root / "configuracao_otimizacao.txt").write_text(json.dumps(config), encoding="utf-8")
        controller = CalibrationController(root, run_root=root / "runs", solver_command=["true"])
        doe_trials = controller.config["campaign"]["doe_trials"]
        points = __import__("optimizer.article_optimizer", fromlist=["latin_hypercube"]).latin_hypercube(doe_trials, {"sf11": (0.3, 1.1)}, 7)
        controller._append_csv("results.csv", __import__("optimizer.calibration_controller", fromlist=["RESULT_FIELDS"]).RESULT_FIELDS, {"trial_id": 1, "phase": "doe", "status": "complete", "sf11": points[0]["sf11"]})
        for trial_id in (2, 3, 4):
            workspace = controller.run_root / f"trial_{trial_id:04d}"
            workspace.mkdir()
            (workspace / "parameters.json").write_text(json.dumps({"parameters": points[1], "seeds": {}}), encoding="utf-8")
            controller._record_status(trial_id=trial_id, phase="doe", status="rejected", stage="measurement", message="drawdown measurement rejected", workspace=str(workspace), started_at="", finished_at="")
        scheduled = []
        controller._run_points = lambda run_points, phase, first_trial_id, initial_retries=None: scheduled.extend(run_points)
        score_calls = [0]
        def scored_results():
            score_calls[0] += 1
            return [{"status": "complete", "phase": "doe"}] if score_calls[0] == 1 else [{"status": "complete", "phase": "doe"} for _ in range(3)]
        controller.score_completed_results = scored_results

        controller.run_doe()

        self.assertEqual(scheduled, points[2:])

    def test_run_doe_skips_pending_points_when_minimum_is_already_complete(self):
        root = Path(tempfile.mkdtemp())
        config = {
            "parameters": {"sf11": {"min": 0.3, "max": 1.1, "initial": 0.95}},
            "minimum_doe_complete": 2,
            "campaign": {"doe_trials": 3, "max_trials": 3},
        }
        (root / "configuracao_otimizacao.txt").write_text(json.dumps(config), encoding="utf-8")
        controller = CalibrationController(root, run_root=root / "runs", solver_command=["true"])
        completed = [{"phase": "doe", "status": "complete"}, {"phase": "doe", "status": "complete"}]
        controller.score_completed_results = lambda: completed
        controller._run_points = lambda *args, **kwargs: self.fail("DOE should not restart pending points")
        controller._resume_doe_state = lambda: (set(), {}, 47)
        self.assertIs(controller.run_doe(), completed)
        self.assertEqual(controller._next_trial_id, 47)

    def test_doe_backfill_keeps_six_workers_when_two_results_are_missing(self):
        root = Path(tempfile.mkdtemp())
        config = {
            "parameters": {"sf11": {"min": 0.3, "max": 1.1, "initial": 0.95}},
            "minimum_doe_complete": 45,
            "doe_contingency_batch_size": 6,
        }
        (root / "configuracao_otimizacao.txt").write_text(json.dumps(config), encoding="utf-8")
        controller = CalibrationController(root, run_root=root / "runs", solver_command=["true"])
        self.assertEqual(controller._doe_backfill_batch_size(43), 6)
        self.assertEqual(controller._doe_backfill_batch_size(45), 0)

    def test_active_bounds_exclude_low_rf12_without_discarding_existing_results(self):
        root = Path(tempfile.mkdtemp())
        config = {
            "parameters": {"rf12": {"min": 0.15, "max": 0.30, "initial": 0.15}},
            "mpi_processes_per_trial": 1,
            "parallel_trials": 1,
            "random_seed": 7,
        }
        (root / "configuracao_otimizacao.txt").write_text(json.dumps(config), encoding="utf-8")
        controller = CalibrationController(root, run_root=root / "runs", solver_command=["true"])
        self.assertFalse(controller._is_within_active_bounds({"rf12": 0.149999}))
        self.assertTrue(controller._is_within_active_bounds({"rf12": 0.15}))
        self.assertTrue(controller._is_within_active_bounds({"rf12": 0.30}))

    def test_run_points_continues_after_exhausted_measurement_retries(self):
        root = Path(tempfile.mkdtemp())
        config = {"parameters": {"sf11": {"min": 0.3, "max": 1.1, "initial": 0.95}}, "mpi_processes_per_trial": 1, "parallel_trials": 1, "random_seed": 7, "measurement_retry_limit": 0, "campaign": {"max_trials": 2}}
        (root / "configuracao_otimizacao.txt").write_text(json.dumps(config), encoding="utf-8")
        controller = CalibrationController(root, run_root=root / "runs", solver_command=["true"])
        calls = []

        def execute_trial(trial_id, parameters, phase="doe", seeds=None):
            calls.append(parameters["sf11"])
            if parameters["sf11"] == 0.30:
                raise MeasurementRejected("drawdown measurement rejected")
            return {"trial_id": trial_id, "phase": phase, "status": "complete", **parameters}

        controller.execute_trial = execute_trial
        results = controller._run_points([{"sf11": 0.30}, {"sf11": 0.31}], "doe", 1)

        self.assertEqual(calls, [0.30, 0.31])
        self.assertEqual([result["sf11"] for result in results], [0.31])

    def test_run_points_retries_rejected_point_and_keeps_workers_occupied(self):
        root = Path(tempfile.mkdtemp())
        config = {
            "parameters": {"sf11": {"min": 0.3, "max": 1.1, "initial": 0.95}},
            "mpi_processes_per_trial": 1,
            "parallel_trials": 2,
            "random_seed": 7,
            "campaign": {"max_trials": 3},
            "measurement_retry_limit": 2,
        }
        (root / "configuracao_otimizacao.txt").write_text(json.dumps(config), encoding="utf-8")
        controller = CalibrationController(root, run_root=root / "runs", solver_command=["true"])
        calls = []
        rejected_once = set()

        def execute_trial(trial_id, parameters, phase="doe", seeds=None):
            calls.append((trial_id, parameters["sf11"]))
            time.sleep(0.01)
            if parameters["sf11"] == 0.31 and parameters["sf11"] not in rejected_once:
                rejected_once.add(parameters["sf11"])
                raise MeasurementRejected("insufficient radial support")
            return {"trial_id": trial_id, "phase": phase, "status": "complete", **parameters}

        controller.execute_trial = execute_trial
        results = controller._run_points([{"sf11": 0.30}, {"sf11": 0.31}, {"sf11": 0.32}], "doe", 1)
        self.assertEqual(len(results), 3)
        self.assertEqual(sorted(result["sf11"] for result in results), [0.30, 0.31, 0.32])
        self.assertEqual(sum(value == 0.31 for _, value in calls), 2)
        self.assertEqual(len({trial_id for trial_id, _ in calls}), 4)
    def test_geometry_diagnostic_is_saved_and_never_retried(self):
        root = Path(tempfile.mkdtemp())
        config = {"parameters": {"CED11": {"min": 0.0, "max": 10.0, "initial": 1.0}}, "mpi_processes_per_trial": 1, "parallel_trials": 1, "random_seed": 7, "campaign": {"max_trials": 2}}
        (root / "configuracao_otimizacao.txt").write_text(json.dumps(config), encoding="utf-8")
        controller = CalibrationController(root, run_root=root / "runs", solver_command=["true"])
        calls = []

        def execute_trial(trial_id, parameters, phase="doe", seeds=None):
            calls.append(parameters["CED11"])
            if parameters["CED11"] == 9.0:
                raise GeometryDiagnostic("bridge-like profile; angle unavailable")
            return {"trial_id": trial_id, "phase": phase, "status": "complete", **parameters}

        controller.execute_trial = execute_trial
        results = controller._run_points([{"CED11": 9.0}, {"CED11": 1.0}], "doe", 1)

        self.assertEqual(calls, [9.0, 1.0])
        self.assertEqual([row["CED11"] for row in results], [1.0])
        diagnostics = (root / "runs" / "diagnostics.csv").read_text(encoding="utf-8")
        self.assertIn("bridge-like profile; angle unavailable", diagnostics)
        completed, _, _ = controller._resume_doe_state()
        self.assertIn(controller._point_key({"CED11": 9.0}), completed)
        report = json.loads((root / "runs" / "cohesion_limits.json").read_text(encoding="utf-8"))
        self.assertEqual(report["problematic_observations"][0]["CED11"], 9.0)
        self.assertEqual(report["lowest_problematic_CED11"], 9.0)

    def test_successful_trials_are_recorded_as_measurable_cohesion_values(self):
        root = Path(tempfile.mkdtemp())
        config = {"parameters": {"CED11": {"min": 0.0, "max": 10.0, "initial": 1.0}}, "mpi_processes_per_trial": 1, "parallel_trials": 1, "campaign": {"max_trials": 1}}
        (root / "configuracao_otimizacao.txt").write_text(json.dumps(config), encoding="utf-8")
        controller = CalibrationController(root, run_root=root / "runs", solver_command=["true"])
        controller.execute_trial = lambda trial_id, parameters, phase="doe", seeds=None: {"trial_id": trial_id, "phase": phase, "status": "complete", **parameters}
        controller._run_points([{"CED11": 2.5}], "doe", 1)
        report = json.loads((root / "runs" / "cohesion_limits.json").read_text(encoding="utf-8"))
        self.assertEqual(report["highest_measurable_CED11"], 2.5)

    def test_cohesion_boundary_search_sets_ceiling_to_highest_measurable_probe(self):
        root = Path(tempfile.mkdtemp())
        config = {
            "parameters": {
                "sf11": {"min": 0.2, "max": 2.0, "initial": 0.5},
                "rf11": {"min": 0.2, "max": 2.0, "initial": 0.8},
                "CED11": {"min": 0.0, "max": 10.0, "initial": 1.0},
            },
            "cohesion_guard": {"enabled": True, "max_boundary_probes": 4},
            "mpi_processes_per_trial": 1,
            "parallel_trials": 1,
            "campaign": {"max_trials": 10},
        }
        (root / "configuracao_otimizacao.txt").write_text(json.dumps(config), encoding="utf-8")
        controller = CalibrationController(root, run_root=root / "runs", solver_command=["true"])

        def execute_trial(trial_id, parameters, phase="doe", seeds=None):
            if parameters["CED11"] > 5.0:
                raise GeometryDiagnostic("cohesive arch; angle not measurable")
            return {"trial_id": trial_id, "phase": phase, "status": "complete", **parameters}

        controller.execute_trial = execute_trial
        with redirect_stdout(io.StringIO()):
            controller._run_points([{"sf11": 0.5, "rf11": 0.8, "CED11": 10.0}], "doe", 1)
            controller._calibrate_cohesion_ceiling()

        report = json.loads((root / "runs" / "cohesion_boundary_search.json").read_text(encoding="utf-8"))
        self.assertEqual(report["status"], "complete")
        self.assertEqual(report["highest_measurable_CED11"], 5.0)
        self.assertEqual(report["lowest_problematic_CED11"], 5.625)
        self.assertEqual(controller.config["parameters"]["CED11"]["max"], 5.0)

    def test_explicit_campaign_sizes_override_defaults_and_require_consistent_total(self):
        root = Path(tempfile.mkdtemp())
        config = {
            "parameters": {"sf11": {"min": 0.3, "max": 1.1, "initial": 0.95}},
            "doe_trials": 45,
            "refinement_trials": 48,
            "validation_trials": 15,
            "max_trials": 108,
        }
        (root / "configuracao_otimizacao.txt").write_text(json.dumps(config), encoding="utf-8")
        controller = CalibrationController(root, run_root=root / "runs", solver_command=["true"])
        self.assertEqual(controller.config["campaign"], {"quadratic_terms": 3, "doe_trials": 45, "refinement_trials": 48, "validation_trials": 15, "max_trials": 108})
        config["max_trials"] = 107
        (root / "configuracao_otimizacao.txt").write_text(json.dumps(config), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "max_trials"):
            CalibrationController(root, run_root=root / "bad-runs", solver_command=["true"])


if __name__ == "__main__":
    unittest.main()
