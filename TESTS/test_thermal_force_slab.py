#!/usr/bin/env python3
"""
CTest for the slab driver of the thermal-force applet (i_option = 16).

W3+ in H+ (n = 1e14 cm^-3, T = 50 eV, B = 1e4 G e_z), grad T = 100 eV/m e_x,
ions start at the origin from the 50 eV Maxwellian; collision step
0.01 tau_s, mean drift between Omega_a t = 300 and 700 from the paired
difference to the grad T = 0 control. The window starts after the initial
transient: all ions start at the same particle position, which correlates
gyrophase and guiding-centre position in the full orbit until collisions
decorrelate the velocities (tau_s = 110/Omega_a); the full orbit then runs up
to 10 um ahead of the guiding centre before both drift alike. Against our
first-order B07 targets
(TESTS/REFERENCE/homma/B07_targets.csv):

  guiding_centre, background at the guiding centre  ->  gc_dXdt
  guiding_centre, background on the gyro orbit       ->  gyro_orbit_dXdt
  full_orbit                                         ->  gyro_orbit_dXdt

dY/dt must vanish for the guiding-centre model. The full orbit carries a small
physical transverse drift of order (nu/Omega_a) dX/dt, which is bounded only.
Tolerance: 4 standard errors plus 1 % (second-order effects in grad T, the
explicit collision step and, for the full orbit, the aliasing between
collision step and gyration (Omega_a dt = 1.1) are below that here).
"""
import csv
import os
import shutil
import subprocess
import sys
from pathlib import Path

try:
    import f90nml
except ImportError:
    sys.exit("f90nml is required to run this test (pip install f90nml)")


def target(root, column):
    with open(root / "TESTS" / "REFERENCE" / "homma" / "B07_targets.csv") as f:
        for r in csv.DictReader(line for line in f if not line.startswith("#")):
            if int(r["Z"]) == 3 and float(r["gradT_x_eV_m"]) == 100.0:
                return float(r[column])
    raise RuntimeError("no W3+ 100 eV/m row in B07_targets.csv")


def run_case(root, binary, work, name, overrides):
    case_dir = work / name
    if case_dir.exists():
        shutil.rmtree(case_dir)
    example = root / "EXAMPLES" / "thermal_force"
    case_dir.mkdir(parents=True)
    for f in ("gorilla_applets.inp", "gorilla.inp", "tetra_grid.inp"):
        shutil.copy(example / f, case_dir / f)
    nml = f90nml.read(example / "thermal_force.inp")
    nml["thermal_force_nml"].update({
        "driver": "slab", "test_distribution": "maxwellian", "test_temperature_ev": 50.0,
        "test_drift_cm_s": [0.0, 0.0, 0.0], "grad_temperature_ev_cm": [1.0, 0.0, 0.0],
        "output_time_unit": "gyro", "output_times": [0.0, 300.0, 700.0], "dt_coll_over_tau_s": 0.01,
        "n_particles": 6000, "boole_control_run": True, "output_file": "slab.csv"})
    nml["thermal_force_nml"].update(overrides)
    nml.write(case_dir / "thermal_force.inp")
    run = subprocess.run([str(binary)], cwd=case_dir, capture_output=True, text=True)
    if run.returncode != 0:
        print(run.stdout, run.stderr)
        raise RuntimeError(f"{name}: applet failed")
    with open(case_dir / "slab.csv") as f:
        rows = [r for r in csv.DictReader(line for line in f if not line.startswith("#"))
                if r["run"].strip() == "difference"]
    a, b = rows[1], rows[2]
    dt = float(b["t_s"]) - float(a["t_s"])
    # slope between the two times, cm -> m; the standard error is bounded by
    # the sum of the two (the positions are correlated)
    return {"mean_x": (float(b["mean_x"]) - float(a["mean_x"])) / 100.0 / dt,
            "se_x": (float(b["se_x"]) + float(a["se_x"])) / 100.0 / dt,
            "mean_y": (float(b["mean_y"]) - float(a["mean_y"])) / 100.0 / dt,
            "se_y": (float(b["se_y"]) + float(a["se_y"])) / 100.0 / dt}


def main():
    root = Path(os.environ["APPLETS_ROOT"])
    binary = Path(os.environ["GORILLA_APPLETS_BIN"])
    work = Path(os.environ["WORK_DIR"])
    failures = []

    cases = [
        ("gc_at_guiding_centre", {"orbit_model": "guiding_centre", "collision_background_at": "guiding_centre"},
         "gc_dXdt_m_s"),
        ("gc_on_gyro_orbit", {"orbit_model": "guiding_centre", "collision_background_at": "gyro_orbit"},
         "gyro_orbit_dXdt_m_s"),
        ("full_orbit", {"orbit_model": "full_orbit"}, "gyro_orbit_dXdt_m_s"),
    ]
    for name, overrides, column in cases:
        d = run_case(root, binary, work, name, overrides)
        expect = target(root, column)
        tol = 4 * d["se_x"] + 0.01 * abs(expect)
        ok = abs(d["mean_x"] - expect) <= tol
        print(f"{'ok  ' if ok else 'FAIL'} {name}: dX/dt = {d['mean_x']:.4f} +- {d['se_x']:.4f} m/s, "
              f"target {expect:.4f} ({column})")
        if not ok:
            failures.append(name)
        if name.startswith("gc"):
            ok = abs(d["mean_y"]) <= 4 * d["se_y"] + 1e-3 * abs(expect)
        else:
            ok = abs(d["mean_y"]) <= 0.05 * abs(expect)
        print(f"{'ok  ' if ok else 'FAIL'} {name}: dY/dt = {d['mean_y']:.4f} +- {d['se_y']:.4f} m/s")
        if not ok:
            failures.append(name + " dY/dt")

    if failures:
        print(f"FAIL: {len(failures)} checks failed")
        return 1
    print("PASS: thermal-force slab driver")
    return 0


if __name__ == "__main__":
    sys.exit(main())
