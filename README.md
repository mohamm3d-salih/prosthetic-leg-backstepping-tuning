# Backstepping gain tuning for a belt-driven underactuated prosthetic leg

MATLAB code and results for the paper

> M. Salih, S. K. Kadhim, A. A. Oglah, "Does the Metaheuristic Matter? Tuning the Gains of a
> Backstepping Controller for a Belt-Driven Underactuated Prosthetic Leg on a Noisy Digital Loop"
> (under review).

The leg has a thigh and a shank. A single motor on the pelvis drives the knee through a belt, so the
hip is passive and the system has two degrees of freedom but only one input. A backstepping
controller tracks the shank angle; its two gains `c3` and `c4` are tuned by six metaheuristics
(DE, PSO, BA, GWO, CMA-ES and PSO with pattern search) on a 1 kHz digital loop with sensor noise,
a measurement delay and a heavier shank, and compared with trial-and-error gains, an analytical
rule, a systematic manual tuning and the Nelder-Mead method.

Everything is simulation. No toolbox is required; the Parallel Computing Toolbox speeds things up
if you have it. The results were produced with MATLAB R2023b on Windows.

## Quick start

```matlab
cd path/to/this/folder
run_paper1('check')      % about 5 minutes: checks the simulator and estimates run times
run_paper1('stats')      % statistics of the 90 tuning runs stored in results/
```

Other options:

| Command | What it does | Time (8 workers) |
|---|---|---|
| `run_paper1('check')` | quick checks and run-time estimate | ~5 min |
| `run_paper1('landscape')` | cost map over a 20 x 20 grid of (c3, c4) | ~12 min |
| `run_paper1('tune')` | six optimisers x 15 runs (finished runs are skipped) | ~2-3 h |
| `run_paper1('stats')` | statistics of the optimisers | ~1 min |
| `run_paper1('evaluate')` | tuned gains vs trial and error: noise, transients, robustness, gait, Monte Carlo | ~15 min |
| `run_paper1('extra')` | alternatives to metaheuristics, cost weights, actuator lag, belt compliance, velocity noise, knee flexion | ~25 min |
| `run_paper1` | all of the above | several hours |

`results/` already contains every run used in the paper, so the statistics and the evaluation can
be repeated without tuning again. To repeat the tuning from scratch, move the `p1_run_*.mat` files
out of `results/` first. Each step writes a log to `logs/`.

## Folders

```
run_paper1.m         entry point
src/                 the steps of the study (p1_step0 ... p1_step5) and helpers
lib/model/           leg model, reference trajectories, cost and performance measures
lib/control/         control laws (only backstepping, type 0, is used here)
lib/simulation/      digital-loop simulator (ZOH at 1 kHz, noise, delay, RK4 plant)
lib/optimizers/      DE, PSO, BA, GWO, CMA-ES, PSO + pattern search
data/mc_samples.csv  the 300 random legs and conditions of the Monte Carlo test
results/             results of the paper (JSON, CSV and one .mat file per tuning run)
```

The prefixes `b2_`, `n2_` and `r2_` in `lib/` come from earlier versions of the simulator and are
kept so that the code matches the logs of the paper.

## Main settings

All settings are in `src/p1_settings.m` and `lib/model/b2_settings.m`:

- sampling period 1 ms, RK4 plant with 5 sub-steps
- angle noise 1e-3 rad, velocity noise 1e-2 rad/s (gyroscope), fixed noise seed 2026 during tuning
- cost: mean over three scenarios of `int t|e| dt + 1e-4 int u^2 dt`, T = 20 s
  (test leg; test leg with 3 ms delay; shank 50 % heavier with twice the noise)
- search box c3 in [1, 500], c4 in [1, 1500]; population 8; 12 generations (104 evaluations); 15 runs

## Result files

| File | Content |
|---|---|
| `p1_step0_check.json` | simulator check and timing |
| `p1_landscape.csv/.json` | cost map |
| `p1_run_<alg>_<seed>.mat` | one tuning run: best gains, cost, convergence history |
| `p1_stats.json`, `p1_convergence.csv` | statistics of the optimisers |
| `p1_evaluate.json`, `p1_robustness.csv`, `p1_mc.csv`, `p1_timehist_*.csv` | evaluation of the designs |
| `p1_revision.json`, `p1_scan.csv`, `p1_baseline_hist.csv` | alternatives and practical limits |

## Citation

If you use this code, please cite the paper above (the full reference will be added after
publication). See also `CITATION.cff`.

## License

MIT, see `LICENSE`.
