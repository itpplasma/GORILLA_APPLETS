#!/usr/bin/env python3
"""
Temperature-screening test of Homma & Hatayama (2014), Sec. 3.2, Figs. 1-4,
with the thermal-force applet (i_option = 16, slab driver).

W3+ test ions in H+ (n = 1e20 m^-3, T = 50 eV at the origin, lnL = 15),
B = 1 T e_z, grad T = 300 eV/m e_x, 1e4 ions starting at the origin from the
50 eV Maxwellian, collision steps 0.05 and 0.15 tau_s, up to Omega_a t = 700.
Each step size runs with both orbit models of the slab driver:
  full_orbit      Homma's own model (exact helix between collision kicks)
  guiding_centre  the GORILLA model (random gyrophase, kick, guiding-centre shift)
The first two ions of the full-orbit runs are traced for Figs. 1 and 3.

Usage: python3 run_homma2014.py [work_dir]   (default: ./run)
Then:  python3 plot_homma2014.py [work_dir]
The applet binary is ../../BUILD/gorilla_applets_main.x or $GORILLA_APPLETS_BIN.
"""
import os
import shutil
import subprocess
import sys
from pathlib import Path

import f90nml

HERE = Path(__file__).resolve().parent
TEMPLATE = HERE.parent / "thermal_force"
BINARY = Path(os.environ.get("GORILLA_APPLETS_BIN", HERE.parent.parent / "BUILD" / "gorilla_applets_main.x"))

STEPS = [0.05, 0.15]                      # dt_coll / tau_s, Homma Figs. 1-2 and 3-4
MODELS = ["full_orbit", "guiding_centre"]

CASE = {
    "driver": "slab",
    "test_distribution": "maxwellian",
    "test_temperature_ev": 50.0,
    "test_drift_cm_s": [0.0, 0.0, 0.0],
    "grad_temperature_ev_cm": [3.0, 0.0, 0.0],      # 300 eV/m
    "collision_background_at": "gyro_orbit",
    "output_time_unit": "gyro",
    "output_times": [5.0 * k for k in range(141)],  # Omega_a t = 0, 5, ..., 700
    "n_particles": 10000,
    "boole_control_run": True,
    "seed": 2014,
    "trajectory_dt_gyro": 0.2,
}


def main():
    work = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE / "run"
    for dt in STEPS:
        for model in MODELS:
            case_dir = work / f"{model}_dt{dt}"
            if case_dir.exists():
                shutil.rmtree(case_dir)
            case_dir.mkdir(parents=True)
            for f in ("gorilla_applets.inp", "gorilla.inp", "tetra_grid.inp"):
                shutil.copy(TEMPLATE / f, case_dir / f)
            nml = f90nml.read(TEMPLATE / "thermal_force.inp")
            nml["thermal_force_nml"].update(CASE)
            nml["thermal_force_nml"].update({
                "orbit_model": model,
                "dt_coll_over_tau_s": dt,
                "n_trajectories": 2 if model == "full_orbit" else 0,
                "output_file": "screening.csv",
            })
            nml.write(case_dir / "thermal_force.inp")
            run = subprocess.run([str(BINARY)], cwd=case_dir, capture_output=True, text=True)
            if run.returncode != 0:
                print(run.stdout, run.stderr)
                sys.exit(f"applet failed in {case_dir}")
            print(f"done: {case_dir}")


if __name__ == "__main__":
    main()
