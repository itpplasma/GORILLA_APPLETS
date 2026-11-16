#!/usr/bin/env python3
"""
CTest for the tokamak driver of the thermal-force applet (i_option = 16).

Runs EXAMPLES/thermal_force_tokamak: W3+ in H+ (n = 1e14 cm^-3, T = 50 eV at
the start), grad T = 100 eV/m along R, GORILLA's analytic circular tokamak with
R0 = 1 km, a = 1 m, B0 = 1 T, q = 10 (one toroidal field period of 2 pi/1000).
Guiding-centre orbits by GORILLA's pusher, collision step 0.01 tau_s, drift of
the paired difference to the grad T = 0 control between Omega_a t = 300 and 700.

  collision shift on   d<R>/dt = gyro_orbit_dXdt of TESTS/REFERENCE/homma/B07_targets.csv
                       (slab mapping e_x <-> e_R), d<Z>/dt ~ 0
  collision shift off  no screening drift: d<R>/dt = 0 (velocity-only collisions)
  field periods        collisionless (n = 1e6 cm^-3, dt = 1e-6 tau_s) to Omega_a t = 7000:
                       the toroidal spread std(R_s dphi) = sqrt(<v_par^2>) t over about
                       3.6 field periods, i.e. the periodic phi crossings are counted

Tolerance: 4 standard errors plus 2 % (second-order effects in grad T, the
explicit collision step, and the curvature and grad-B drift, which the control
run subtracts to first order).
"""
import csv
import math
import os
import shutil
import subprocess
import sys
from pathlib import Path

try:
    import f90nml
except ImportError:
    sys.exit("f90nml is required to run this test (pip install f90nml)")

EXAMPLE_FILES = ("gorilla_applets.inp", "gorilla.inp", "tetra_grid.inp", "field_divB0.inp", "preload_for_SYNCH.inp")


def target(root):
    with open(root / "TESTS" / "REFERENCE" / "homma" / "B07_targets.csv") as f:
        for r in csv.DictReader(line for line in f if not line.startswith("#")):
            if int(r["Z"]) == 3 and float(r["gradT_x_eV_m"]) == 100.0:
                return float(r["gyro_orbit_dXdt_m_s"])
    raise RuntimeError("no W3+ 100 eV/m row in B07_targets.csv")


def run_rows(root, binary, work, name, overrides):
    case_dir = work / name
    if case_dir.exists():
        shutil.rmtree(case_dir)
    example = root / "EXAMPLES" / "thermal_force_tokamak"
    case_dir.mkdir(parents=True)
    for f in EXAMPLE_FILES:
        if (example / f).exists():
            shutil.copy(example / f, case_dir / f)
    nml = f90nml.read(example / "thermal_force.inp")
    nml["thermal_force_nml"].update({"output_times": [0.0, 300.0, 700.0], "n_particles": 1500,
                                     "output_file": "out.csv"})
    nml["thermal_force_nml"].update(overrides)
    nml.write(case_dir / "thermal_force.inp")
    run = subprocess.run([str(binary)], cwd=case_dir, capture_output=True, text=True)
    if run.returncode != 0:
        print(run.stdout[-3000:], run.stderr[-3000:])
        raise RuntimeError(f"{name}: applet failed")
    with open(case_dir / "out.csv") as f:
        lines = f.readlines()
    lost = [line for line in lines if line.startswith("# lost_ions_gradient_control:")]
    rows = list(csv.DictReader(line for line in lines if not line.startswith("#")))
    return rows, (lost[0].split(":")[1].strip() if lost else "?")


def run_case(root, binary, work, name, overrides):
    rows, lost = run_rows(root, binary, work, name, overrides)
    rows = [r for r in rows if r["run"].strip() == "difference"]
    a, b = rows[1], rows[2]
    dt = float(b["t_s"]) - float(a["t_s"])
    # columns: x = R - R_s, z = Z - Z_s (cm); slope between the two times in m/s
    return {"dR": (float(b["mean_x"]) - float(a["mean_x"])) / 100.0 / dt,
            "se_R": (float(b["se_x"]) + float(a["se_x"])) / 100.0 / dt,
            "dZ": (float(b["mean_z"]) - float(a["mean_z"])) / 100.0 / dt,
            "se_Z": (float(b["se_z"]) + float(a["se_z"])) / 100.0 / dt,
            "lost": lost}


def main():
    root = Path(os.environ["APPLETS_ROOT"])
    binary = Path(os.environ["GORILLA_APPLETS_BIN"])
    work = Path(os.environ["WORK_DIR"])
    failures = []
    expect = target(root)

    d = run_case(root, binary, work, "shift_on", {"boole_collision_shift": True})
    ok = abs(d["dR"] - expect) <= 4 * d["se_R"] + 0.02 * abs(expect)
    print(f"{'ok  ' if ok else 'FAIL'} shift on: dR/dt = {d['dR']:.4f} +- {d['se_R']:.4f} m/s, target {expect:.4f}, "
          f"lost ions {d['lost']}")
    if not ok:
        failures.append("shift on dR/dt")
    ok = abs(d["dZ"]) <= 4 * d["se_Z"] + 0.02 * abs(expect)
    print(f"{'ok  ' if ok else 'FAIL'} shift on: dZ/dt = {d['dZ']:.4f} +- {d['se_Z']:.4f} m/s")
    if not ok:
        failures.append("shift on dZ/dt")

    d = run_case(root, binary, work, "shift_off", {"boole_collision_shift": False})
    ok = abs(d["dR"]) <= 4 * d["se_R"] + 0.02 * abs(expect)
    print(f"{'ok  ' if ok else 'FAIL'} shift off: dR/dt = {d['dR']:.4f} +- {d['se_R']:.4f} m/s, expected 0")
    if not ok:
        failures.append("shift off dR/dt")

    n = 4000
    rows, _ = run_rows(root, binary, work, "field_periods", {
        "density_cm3": 1.0e6, "dt_coll_over_tau_s": 1.0e-6, "grad_temperature_ev_cm": [0.0, 0.0, 0.0],
        "boole_control_run": False, "n_particles": n, "output_times": [0.0, 7000.0]})
    r = [x for x in rows if x["run"].strip() == "gradient"][-1]
    t = float(r["t_s"])
    spread = float(r["se_y"]) * math.sqrt(n)                               # cm
    v_par_rms = math.sqrt(2.0 * float(r["mean_energy"]) / (183.84 * 1.66053906660e-24) / 3.0)
    ratio = spread / t / v_par_rms
    ok = abs(ratio - 1.0) <= 0.05
    print(f"{'ok  ' if ok else 'FAIL'} field periods: std(R_s dphi) = {spread / 100:.2f} m over 6.29 m periods, "
          f"/ (sqrt(<v_par^2>) t) = {ratio:.4f}")
    if not ok:
        failures.append("field periods")

    if failures:
        print(f"FAIL: {len(failures)} checks failed")
        return 1
    print("PASS: thermal-force tokamak driver")
    return 0


if __name__ == "__main__":
    sys.exit(main())
