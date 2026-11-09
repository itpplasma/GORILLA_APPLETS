#!/usr/bin/env python3
"""
B08 convergence campaign for the thermal-force applet (i_option = 16).

Scans (each a subcommand, results appended to <work_dir>/b08_<scan>.csv):

  dt     B07 screening drift versus collision step: W1+, W3+, W10+ in H+,
         grad T = 100 eV/m e_x, for the guiding centre with the background at
         the guiding centre ('gc_gc') and on the gyro orbit ('gc_orbit'), and
         the full orbit as a reference ('full'). Drift = least-squares slope of
         the paired difference <X>(t) over Omega_a t = 300..700 (after the
         initial transient); one row per seed.
  models The three models on one footing: same scan as 'dt' over
         dt = 0.001..0.1 tau_s, but the drift is fitted over t = 3..7 tau_s,
         after the initial transient for every charge state (the transient
         lasts about 2 tau_s, i.e. Omega_a t ~ 670 for W1+).
  relax  B04 relaxation of the drifting bi-Maxwellian W3+ population versus
         collision step: p_par, p_perp, u_par at t = 1, 3, 10 tau_s.
  n      Statistics versus ion number: W3+ gc_orbit at dt = 0.02 tau_s, average
         drift <X>(t_end)/t_end of the gradient run alone and of the paired
         difference, with the reported standard error, to compare with the
         seed-to-seed scatter.

Usage: python3 run_b08.py <dt|models|relax|n|all> [work_dir]   (default ./run)
Then:  python3 plot_b08.py [work_dir]
Binary: ../../BUILD/gorilla_applets_main.x or $GORILLA_APPLETS_BIN;
threads: $OMP_NUM_THREADS (default 8, physical cores).
"""
import csv
import os
import shutil
import subprocess
import sys
from pathlib import Path

import f90nml
import numpy as np

HERE = Path(__file__).resolve().parent
TEMPLATE = HERE.parent / "thermal_force"
BINARY = Path(os.environ.get("GORILLA_APPLETS_BIN", HERE.parent.parent / "BUILD" / "gorilla_applets_main.x"))
ENV = dict(os.environ, OMP_NUM_THREADS=os.environ.get("OMP_NUM_THREADS", "8"))

DT_GRID = [0.001, 0.005, 0.01, 0.025, 0.05, 0.10, 0.15, 0.25, 0.50]
MODELS = {"gc_gc": {"orbit_model": "guiding_centre", "collision_background_at": "guiding_centre"},
          "gc_orbit": {"orbit_model": "guiding_centre", "collision_background_at": "gyro_orbit"},
          "full": {"orbit_model": "full_orbit"}}
FIT_WINDOW = [300.0 + 50.0 * k for k in range(9)]      # Omega_a t = 300, 350, ..., 700

MAXWELLIAN = {"test_distribution": "maxwellian", "test_temperature_ev": 50.0,
              "test_drift_cm_s": [0.0, 0.0, 0.0]}
SCREENING = {**MAXWELLIAN, "driver": "slab", "grad_temperature_ev_cm": [1.0, 0.0, 0.0],
             "output_time_unit": "gyro", "output_times": [0.0] + FIT_WINDOW,
             "boole_control_run": True, "output_file": "out.csv"}


def run_case(case_dir, overrides):
    if case_dir.exists():
        shutil.rmtree(case_dir)
    case_dir.mkdir(parents=True)
    for f in ("gorilla_applets.inp", "gorilla.inp", "tetra_grid.inp"):
        shutil.copy(TEMPLATE / f, case_dir / f)
    nml = f90nml.read(TEMPLATE / "thermal_force.inp")
    nml["thermal_force_nml"].update(overrides)
    nml.write(case_dir / "thermal_force.inp")
    run = subprocess.run([str(BINARY)], cwd=case_dir, capture_output=True, text=True, env=ENV)
    if run.returncode != 0:
        print(run.stdout, run.stderr)
        sys.exit(f"applet failed in {case_dir}")
    with open(case_dir / "out.csv") as f:
        rows = [{k: (v.strip() if k == "run" else float(v)) for k, v in r.items()}
                for r in csv.DictReader(line for line in f if not line.startswith("#"))]
    shutil.rmtree(case_dir)
    return rows


def slope(rows, run, key, column="omega_t", start=FIT_WINDOW[0]):
    """Least-squares slope of <key>(t) over the fit window (column >= start), cm -> m per s."""
    pts = [(r["t_s"], r[key] / 100.0) for r in rows if r["run"] == run and r[column] >= start - 1e-9]
    t, y = np.array(pts).T
    return np.polyfit(t, y, 1)[0]


def append(path, header, row):
    new = not path.exists()
    with open(path, "a", newline="") as f:
        w = csv.writer(f)
        if new:
            w.writerow(header)
        w.writerow(row)


def scan_dt(work, seeds=4, n=10000):
    out = work / "b08_dt.csv"
    header = ["model", "Z", "dt_over_tau_s", "omega_dt", "seed", "n", "dXdt_m_s", "dYdt_m_s", "omega_tau_s"]
    for model, opts in MODELS.items():
        for z in (1, 3, 10):
            for dt in DT_GRID:
                for seed in range(1, seeds + 1):
                    rows = run_case(work / "tmp", {**SCREENING, **opts, "test_charge_number": float(z),
                                                   "dt_coll_over_tau_s": dt, "n_particles": n, "seed": seed})
                    omega_tau = header_value(work, rows)
                    append(out, header, [model, z, dt, omega_tau * dt, seed, n,
                                         slope(rows, "difference", "mean_x"), slope(rows, "difference", "mean_y"),
                                         omega_tau])
                print(f"dt scan: {model} Z={z} dt={dt} done", flush=True)


def scan_models(work, seeds=4, n=10000):
    out = work / "b08_models.csv"
    header = ["model", "Z", "dt_over_tau_s", "omega_dt", "seed", "n", "dXdt_m_s", "dYdt_m_s", "omega_tau_s"]
    window = [3.0 + 0.5 * k for k in range(9)]           # t / tau_s = 3, 3.5, ..., 7
    for model, opts in MODELS.items():
        for z in (1, 3, 10):
            for dt in (0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1):
                for seed in range(1, seeds + 1):
                    rows = run_case(work / "tmp", {**SCREENING, **opts, "test_charge_number": float(z),
                                                   "output_time_unit": "tau_s", "output_times": [0.0] + window,
                                                   "dt_coll_over_tau_s": dt, "n_particles": n, "seed": seed})
                    omega_tau = header_value(work, rows)
                    append(out, header, [model, z, dt, omega_tau * dt, seed, n,
                                         slope(rows, "difference", "mean_x", "t_over_tau_s", window[0]),
                                         slope(rows, "difference", "mean_y", "t_over_tau_s", window[0]),
                                         omega_tau])
                print(f"models scan: {model} Z={z} dt={dt} done", flush=True)


def header_value(work, rows):
    """Omega_a tau_s from the output times: omega_t/(t/tau_s)."""
    r = next(r for r in rows if r["t_over_tau_s"] > 0)
    return r["omega_t"] / r["t_over_tau_s"]


def scan_relax(work, seeds=4, n=20000):
    out = work / "b08_relax.csv"
    header = ["dt_over_tau_s", "seed", "n", "t_over_tau_s", "u_par_cm_s", "se_u_par", "p_par_over_T", "se_p_par",
              "p_perp_over_T", "se_p_perp"]
    t_b = 50.0 * 1.6022e-12
    for dt in DT_GRID:
        for seed in range(1, seeds + 1):
            rows = run_case(work / "tmp", {"driver": "velocity", "test_distribution": "drifting_bimaxwellian",
                                           "output_time_unit": "tau_s", "output_times": [1.0, 3.0, 10.0],
                                           "dt_coll_over_tau_s": dt, "n_particles": n, "seed": seed,
                                           "boole_control_run": False, "output_file": "out.csv"})
            for r in rows:
                append(out, header, [dt, seed, n, r["t_over_tau_s"], r["u_par"], r["se_u_par"],
                                     r["p_par"] / t_b, r["se_p_par"] / t_b, r["p_perp"] / t_b, r["se_p_perp"] / t_b])
        print(f"relax scan: dt={dt} done", flush=True)


def scan_n(work):
    out = work / "b08_n.csv"
    header = ["n", "seed", "dXdt_gradient_m_s", "se_gradient_m_s", "dXdt_difference_m_s", "se_difference_m_s"]
    for n, seeds in ((1000, 32), (10000, 16), (100000, 8), (1000000, 4)):
        for seed in range(1, seeds + 1):
            rows = run_case(work / "tmp", {**SCREENING, **MODELS["gc_orbit"], "dt_coll_over_tau_s": 0.02,
                                           "n_particles": n, "seed": seed})
            # average drift over 0..t_end and its reported standard error
            end = {r["run"]: r for r in rows if abs(r["omega_t"] - FIT_WINDOW[-1]) < 1e-6}
            v = {k: (r["mean_x"] / 100.0 / r["t_s"], r["se_x"] / 100.0 / r["t_s"]) for k, r in end.items()}
            append(out, header, [n, seed, *v["gradient"], *v["difference"]])
        print(f"n scan: n={n} done", flush=True)


def main():
    if len(sys.argv) < 2 or sys.argv[1] not in ("dt", "models", "relax", "n", "all"):
        sys.exit(__doc__)
    work = Path(sys.argv[2]) if len(sys.argv) > 2 else HERE / "run"
    work.mkdir(parents=True, exist_ok=True)
    scans = {"dt": scan_dt, "models": scan_models, "relax": scan_relax, "n": scan_n}
    for name in (scans if sys.argv[1] == "all" else [sys.argv[1]]):
        f = work / f"b08_{name}.csv"
        if f.exists():
            f.unlink()
        scans[name](work)


if __name__ == "__main__":
    main()
