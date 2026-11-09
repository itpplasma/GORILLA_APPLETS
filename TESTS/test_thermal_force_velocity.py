#!/usr/bin/env python3
"""
CTest for the velocity driver of the thermal-force applet (i_option = 16).

Runs EXAMPLES/thermal_force with overrides (W3+ in H+, n = 1e14 cm^-3,
T = 50 eV, B = 1e4 G e_z) and checks:

  equilibrium  a Maxwellian at the background temperature stays one over 5
               tau_s: p_par, p_perp per ion equal T_b, <v^4> = 15 (T_b/m_a)^2,
               mean velocity zero (fluctuation-dissipation balance of A and D)
  B05          parallel gradient 5 eV/m: the paired difference to the grad T = 0
               control gives <v_z>(t) = a_z t with a = <F_T>/m_a, the Maxwellian
               thermal-force moment of TESTS/REFERENCE/homma/B03_oracle.csv
  B06          perpendicular gradient 100 eV/m: with gyration the mean
               perpendicular velocity is (sin(W t) a + (1 - cos(W t)) a x b)/W

The results do not depend on the number of OpenMP threads (one random stream
per test ion). Paths come from TESTS/CMakeLists.txt via environment variables.
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

EV2ERG = 1.6022e-12                 # GORILLA constants_mod
AMU_G = 1.66053906660e-24
M_A = 183.84 * AMU_G
T_B = 50.0 * EV2ERG


def run_case(root, binary, work, name, overrides):
    case_dir = work / name
    if case_dir.exists():
        shutil.rmtree(case_dir)
    example = root / "EXAMPLES" / "thermal_force"
    case_dir.mkdir(parents=True)
    for f in ("gorilla_applets.inp", "gorilla.inp", "tetra_grid.inp"):
        shutil.copy(example / f, case_dir / f)
    nml = f90nml.read(example / "thermal_force.inp")
    nml["thermal_force_nml"].update(overrides)
    nml.write(case_dir / "thermal_force.inp")
    run = subprocess.run([str(binary)], cwd=case_dir, capture_output=True, text=True)
    if run.returncode != 0:
        print(run.stdout, run.stderr)
        raise RuntimeError(f"{name}: applet failed")
    rows = {}
    with open(case_dir / "thermal_force_velocity.csv") as f:
        reader = csv.DictReader(line for line in f if not line.startswith("#"))
        for r in reader:
            rows.setdefault(r["run"], []).append({k: float(v) for k, v in r.items() if k != "run"})
    return rows


def thermal_acceleration(root, grad):
    """<F_T>/m_a in cm/s^2 for the Maxwellian W3+ row of the B03 oracle (SI)."""
    with open(root / "TESTS" / "REFERENCE" / "homma" / "B03_oracle.csv") as f:
        reader = csv.DictReader(line for line in f if not line.startswith("#"))
        for r in reader:
            if (r["distribution"] == "maxwellian" and int(r["Z"]) == 3
                    and abs(float(r["ma_amu"]) - 183.84) < 1e-9 and float(r["n_m3"]) == 1e20
                    and float(r["T_eV"]) == 50.0
                    and [float(r[f"gradT_{c}_eV_m"]) for c in "xyz"] == grad):
                m_a_si = 183.84 * 1.66053906660e-27
                return [float(r[f"FT_{c}"]) / m_a_si * 100.0 for c in "xyz"]
    raise RuntimeError(f"no B03 oracle row for grad {grad}")


def check(label, value, expected, tol, failures):
    ok = abs(value - expected) <= tol
    print(f"{'ok  ' if ok else 'FAIL'} {label}: {value:.6e} expected {expected:.6e} (tolerance {tol:.2e})")
    if not ok:
        failures.append(label)


def main():
    root = Path(os.environ["APPLETS_ROOT"])
    binary = Path(os.environ["GORILLA_APPLETS_BIN"])
    work = Path(os.environ["WORK_DIR"])
    failures = []

    # equilibrium: Maxwellian at T_b, 5 tau_s
    rows = run_case(root, binary, work, "equilibrium", {
        "test_distribution": "maxwellian", "test_temperature_ev": 50.0,
        "test_drift_cm_s": [0.0, 0.0, 0.0], "grad_temperature_ev_cm": [0.0, 0.0, 0.0],
        "output_times": [0.0, 5.0], "dt_coll_over_tau_s": 0.01,
        "n_particles": 20000, "boole_control_run": False})
    end = rows["gradient"][-1]
    v4_eq = 15.0 * (T_B / M_A) ** 2
    check("equilibrium p_par/T_b", end["p_par"] / T_B, 1.0, 4 * end["se_p_par"] / T_B, failures)
    check("equilibrium p_perp/T_b", end["p_perp"] / T_B, 1.0, 4 * end["se_p_perp"] / T_B, failures)
    check("equilibrium v4/v4_eq", end["v4"] / v4_eq, 1.0, 4 * end["se_v4"] / v4_eq, failures)
    check("equilibrium u_par [cm/s]", end["u_par"], 0.0, 4 * end["se_u_par"], failures)

    # B05 and B06: short window, paired difference to the control run
    window = {"test_distribution": "maxwellian", "test_temperature_ev": 50.0,
              "test_drift_cm_s": [0.0, 0.0, 0.0], "output_time_unit": "gyro",
              "output_times": [0.0, 0.0025, 0.005, 0.01], "dt_coll_over_tau_s": 0.001,
              "n_particles": 20000, "boole_control_run": True}
    omega = 3.0 * 4.8032e-10 * 1.0e4 / (M_A * 2.9979e10)

    rows = run_case(root, binary, work, "B05", {**window, "grad_temperature_ev_cm": [0.0, 0.0, 0.05]})
    a = thermal_acceleration(root, [0.0, 0.0, 5.0])
    for r in rows["difference"][1:]:
        t = r["t_s"]
        check(f"B05 <dv_z> at Omega t = {r['omega_t']:.4f} [cm/s]", r["mean_vz"], a[2] * t,
              max(4 * r["se_vz"], 1e-3 * abs(a[2] * t)), failures)

    rows = run_case(root, binary, work, "B06", {**window, "grad_temperature_ev_cm": [1.0, 0.0, 0.0]})
    a = thermal_acceleration(root, [100.0, 0.0, 0.0])
    for r in rows["difference"][1:]:
        t = r["t_s"]
        s, c = math.sin(omega * t) / omega, (1.0 - math.cos(omega * t)) / omega
        expect = [s * a[0] + c * a[1], s * a[1] - c * a[0]]     # a x e_z = (a_y, -a_x)
        scale = math.hypot(*expect)
        check(f"B06 <dv_x> at Omega t = {r['omega_t']:.4f} [cm/s]", r["mean_vx"], expect[0],
              max(4 * r["se_vx"], 1e-3 * scale), failures)
        check(f"B06 <dv_y> at Omega t = {r['omega_t']:.4f} [cm/s]", r["mean_vy"], expect[1],
              max(4 * r["se_vy"], 1e-3 * scale), failures)

    if failures:
        print(f"FAIL: {len(failures)} checks failed")
        return 1
    print("PASS: thermal-force velocity driver")
    return 0


if __name__ == "__main__":
    sys.exit(main())
