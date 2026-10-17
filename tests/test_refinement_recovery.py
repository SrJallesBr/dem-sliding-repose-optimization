import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

PROJECT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT))

from optimizer.refinement_recovery import build_recovery_plan, pending_refinement_jobs


class RefinementRecoveryTests(unittest.TestCase):
    def make_run(self, run_root):
        fields = ["trial_id", "phase", "status", "sf11", "rf11", "rf12", "CED11", "score"]
        results = []
        statuses = []
        candidates = [
            {"sf11": 0.3 + i * 0.01, "rf11": 0.1 + i * 0.001, "rf12": 0.21 + i * 0.001, "CED11": 100 + i}
            for i in range(48)
        ]
        for trial_id in range(1, 46):
            results.append({"trial_id": str(trial_id), "phase": "doe", "status": "complete", "sf11": "0.5", "rf11": "0.2", "rf12": "0.24", "CED11": "500", "score": "0.1"})
            statuses.append({"trial_id": str(trial_id), "phase": "doe", "status": "complete", "workspace": "", "stage": "measure"})
        for trial_id in range(46, 70):
            point = candidates[trial_id - 46]
            results.append({"trial_id": str(trial_id), "phase": "refinement", "status": "complete", **{k: str(v) for k, v in point.items()}, "score": "0.2"})
            workspace = run_root / f"trial_{trial_id:04d}"
            workspace.mkdir()
            (workspace / "parameters.json").write_text(json.dumps({"parameters": point, "seeds": {}}), encoding="utf-8")
            statuses.append({"trial_id": str(trial_id), "phase": "refinement", "status": "complete", "workspace": str(workspace), "stage": "measure"})
        for trial_id in range(70, 82):
            point = candidates[trial_id - 46]
            workspace = run_root / f"trial_{trial_id:04d}"
            workspace.mkdir()
            (workspace / "parameters.json").write_text(json.dumps({"parameters": point, "seeds": {}}), encoding="utf-8")
            statuses.append({"trial_id": str(trial_id), "phase": "refinement", "status": "running", "workspace": str(workspace), "stage": "prepare"})
        for name, rows, names in (("results.csv", results, fields), ("status.csv", statuses, ["trial_id", "phase", "status", "workspace", "stage"])):
            with (run_root / name).open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=names)
                writer.writeheader()
                writer.writerows(rows)
        return candidates

    def test_recovery_plan_preserves_old_trials_and_maps_24_remaining_refinements_to_fresh_ids(self):
        run_root = Path(tempfile.mkdtemp())
        candidates = self.make_run(run_root)
        controller = SimpleNamespace(
            parameter_names=("sf11", "rf11", "rf12", "CED11"),
            nominate_refinement_candidates=lambda rows: candidates,
            config={"campaign": {"validation_trials": 15}},
            _trial_seeds=lambda trial_id: {"distribution": 10001 + trial_id, "insertion": 20011 + trial_id},
        )

        plan = build_recovery_plan(controller, run_root)

        self.assertEqual(len(plan["refinement_jobs"]), 24)
        self.assertEqual([job["trial_id"] for job in plan["refinement_jobs"]], list(range(94, 118)))
        self.assertEqual(plan["refinement_jobs"][0]["source_trial_id"], 70)
        self.assertEqual(plan["refinement_jobs"][11]["source_trial_id"], 81)
        self.assertEqual(plan["refinement_jobs"][12]["source_trial_id"], 82)
        self.assertEqual(plan["refinement_jobs"][23]["source_trial_id"], 93)
        self.assertEqual(plan["refinement_jobs"][0]["parameters"], candidates[24])
        self.assertEqual(plan["validation_trial_ids"], list(range(118, 133)))

    def test_resume_plan_skips_a_candidate_already_recorded_complete(self):
        jobs = [
            {"trial_id": 94, "parameters": {"sf11": 0.6, "rf11": 0.2}},
            {"trial_id": 95, "parameters": {"sf11": 0.7, "rf11": 0.3}},
        ]
        results = [{"phase": "refinement", "status": "complete", "sf11": "0.6", "rf11": "0.2"}]

        pending = pending_refinement_jobs(results, jobs, ("sf11", "rf11"))

        self.assertEqual([job["trial_id"] for job in pending], [95])

    def test_recovery_plan_refuses_to_overwrite_existing_new_trial_directory(self):
        run_root = Path(tempfile.mkdtemp())
        candidates = self.make_run(run_root)
        (run_root / "trial_0094").mkdir()
        controller = SimpleNamespace(
            parameter_names=("sf11", "rf11", "rf12", "CED11"),
            nominate_refinement_candidates=lambda rows: candidates,
            config={"campaign": {"validation_trials": 15}},
            _trial_seeds=lambda trial_id: {"distribution": 10001 + trial_id, "insertion": 20011 + trial_id},
        )

        with self.assertRaises(FileExistsError):
            build_recovery_plan(controller, run_root)


if __name__ == "__main__":
    unittest.main()
