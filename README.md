# DEM optimization: sliding and angle of repose

This repository contains the LIGGGHTS optimization project for two paired responses:

1. **Sliding / discharge measurement** — the requested sliding test, measured from the V-test deck. Its standard result image uses three physical X-cuts (−0.10, 0.00 and +0.10 m), reports both faces per cut, and records cut/overall means.
2. **Angle of repose** — measured independently from the repose deck and recorded as the second response for each trial.

The attached reference image is the sliding/discharge measurement figure; the repose measurement is a separate artifact. Targets are objectives, not guaranteed outcomes. The measurement protocol, physical cut mapping, clipping behavior, and validity limits are documented in [`docs/measurement-protocol.md`](docs/measurement-protocol.md).

## Campaign settings

- 12 concurrent trials × 1 LIGGGHTS rank per trial = 12 solver processes maximum, subject to the runnable queue.
- 45 DOE points, 48 refinement points and 15 independent validation points.
- The controller refills a worker as trials complete. A run creates a new timestamped `calibration_runs_*` directory.
- The Windows launcher is `run_visible.cmd`; it starts the controller in WSL `Ubuntu-24.04` and leaves the live dashboard in the terminal.
- Current parameter bounds and targets are in `configuracao_otimizacao.txt`.

Do not commit campaign outputs, dumps, restart files, logs or local audit CSVs. `.gitignore` excludes these by default; review staged files before every push.

## Validation

In WSL, from the repository root:

```bash
python3 -m unittest discover -s tests -v
python3 -u -m optimizer.calibration_controller --preflight --parallel-trials 12 --run-root /tmp/dem_preflight_12
```

Preflight validates setup; it does not start a campaign. `run_visible.cmd` starts a fresh optimization campaign.

## Requirements and limits

Windows 11 with WSL2/Ubuntu-24.04, Python 3 and the compatible LIGGGHTS executable at `/usr/local/bin/liggghts`. A 15 s discharge duration is a fixed maximum, not proof the material settled. Inspect each run's logs, status and measurement artifacts before interpreting its results.

## Repository contents

- `optimizer/` — controller, measurement, scoring, recovery and rendering code.
- `Ang Deslizamento/` and `Ang Repouso/` — input decks and STL assets.
- `tests/` — regression and controller tests.
- `scripts/` — audit helper scripts.
- `configuracao_otimizacao.txt` — campaign parameters and trial counts.
