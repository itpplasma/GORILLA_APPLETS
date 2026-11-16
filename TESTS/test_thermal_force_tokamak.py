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
  flux profiles        EXAMPLES/thermal_force_tokamak_flux (GORILLA example 8, R0 = 1.7 m,
                       B0 = 2 T, T_i(r) = 50 eV + 100 eV/m (r - 25 cm)): screening drift
                       dr/dt = (d<s>/dt)/(ds/dr) over 0..2 tau_s, ds/dr from GORILLA's
                       flux_functions.dat, against the target times (1 T/B)^2 at the start
                       (tolerance 4 standard errors plus 3 %)
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


def run_rows(root, binary, work, name, overrides, example_name="thermal_force_tokamak"):
    case_dir = work / name
    if case_dir.exists():
        shutil.rmtree(case_dir)
    example = root / "EXAMPLES" / example_name
    case_dir.mkdir(parents=True)
    for f in EXAMPLE_FILES:
        if (example / f).exists():
            shutil.copy(example / f, case_dir / f)
    if (example / "profiles").is_dir():
        shutil.copytree(example / "profiles", case_dir / "profiles")
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

    # flux-surface profiles in a realistic circular tokamak
    rows, lost = run_rows(root, binary, work, "flux_profiles", {
        "output_time_unit": "tau_s", "output_times": [0.0, 2.0], "n_particles": 2000},
        example_name="thermal_force_tokamak_flux")
    d = [x for x in rows if x["run"].strip() == "difference"][-1]
    t = float(d["t_s"])
    r_grid, psi_tor = [], []
    with open(work / "flux_profiles" / "flux_functions.dat") as f:
        for line in f:
            if line.strip().startswith("#") or not line.strip():
                continue
            cols = [float(c) for c in line.split()]
            r_grid.append(cols[1]); psi_tor.append(cols[4])
    r_s = 25.0
    i = min(range(len(r_grid)), key=lambda k: abs(r_grid[k] - r_s))
    ds_dr = (psi_tor[i + 1] - psi_tor[i - 1]) / (r_grid[i + 1] - r_grid[i - 1]) / psi_tor[-1]   # 1/cm
    b_start = 2.0e4 * 170.0 / 195.0                                                       # G
    expect_flux = expect * (1.0e4 / b_start) ** 2
    v_r = float(d["mean_s"]) / ds_dr / 100.0 / t                                          # m/s
    se_r = float(d["se_s"]) / ds_dr / 100.0 / t
    ok = abs(v_r - expect_flux) <= 4 * se_r + 0.03 * abs(expect_flux)
    print(f"{'ok  ' if ok else 'FAIL'} flux profiles: dr/dt = {v_r:.4f} +- {se_r:.4f} m/s, target {expect_flux:.4f} "
          f"(B = {b_start/1e4:.3f} T), lost ions {lost}")
    if not ok:
        failures.append("flux profiles")

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
