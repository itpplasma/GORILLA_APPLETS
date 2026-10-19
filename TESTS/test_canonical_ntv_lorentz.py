#!/usr/bin/env python3
"""Gates of the monoenergetic Lorentz mode of the canonical NTV estimator.

Analytic circular tokamak (grid_kind 5, template rmp_fm_stationarity), one
harmonic dB1 = sigma eps B(x) cos(m theta - n phi) relative to the local |B|,
monoenergetic loading, flat temperature and a density exponential in the
poloidal flux, so d ln f0/dp_phi = kappa is a constant. Lorentz pitch-angle
kicks at the constant rate ntv_lorentz_nu, MH corrected to f0.

Gates:
  (a) delta-f / Green-Kubo identity: w = -kappa X with X = int -dH1/dphi dt,
      so the trapezoidal torque integral is tq_deltaf = -kappa X^2/2;
  (b) checkpoints are filled with the end values after the end of the
      (physical-time) trace;
  (c) weight control: reset period beyond the trace time leaves the settled
      torque integral equal to tq_deltaf; a reset period of a quarter of the
      trace bounds the weights (rms w_end smaller than without resets);
  (d) a harmonic that is not periodic on the mesh field period
      (n_tor = 3 with 2 field periods) is rejected at initialisation.
"""

import math
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

DENSITY = 1.2e14
ENERGY_EV = 1000.0
A_MINOR = 69.991288
B0 = 17991.7336           # |B0| of the template [G]
Q0, Q1 = 1.0, 2.0         # template q = q0 + q1 (r/a)^2
TRACE = 2.0e-5
N_CHECK = 4


def env_path(name):
    value = os.environ.get(name)
    if not value:
        sys.exit(f"missing environment variable {name}")
    return Path(value)


APPLETS_ROOT = env_path("APPLETS_ROOT")
BIN = env_path("GORILLA_APPLETS_BIN")
WORK_DIR = env_path("WORK_DIR")
TEMPLATE = APPLETS_ROOT / "TESTS" / "rmp_fm_stationarity"


def set_nml(path, values):
    text = path.read_text()
    for key, value in values.items():
        text, count = re.subn(
            rf"(^\s*{key}\s*=\s*)[^,\n]*", rf"\g<1>{value}", text, flags=re.M
        )
        if count != 1:
            sys.exit(f"{path.name}: expected one {key} entry, found {count}")
    path.write_text(text)


def psi_abs(r):
    """|psi_pol| of the template field, B0 a^2/(2 q1) ln(q/q0)."""
    return B0 * A_MINOR**2 / (2.0 * Q1) * math.log(1.0 + Q1 / Q0 * (r / A_MINOR) ** 2)


def prepare(case_dir, n_particles, reset_dt=0.0, n_tor=2, n_periods=2):
    if case_dir.exists():
        shutil.rmtree(case_dir)
    shutil.copytree(TEMPLATE, case_dir)
    profiles = case_dir / "profiles"
    profiles.mkdir()
    radii = [0.1 + 0.25 * i for i in range(300)]
    r_ref = 0.75 * A_MINOR
    psi_l = 0.3 * psi_abs(A_MINOR)
    dens = [DENSITY * math.exp(-(psi_abs(r) - psi_abs(r_ref)) / psi_l) for r in radii]
    (profiles / "n.dat").write_text(
        "\n".join(f"{r:.8e} {d:.12e}" for r, d in zip(radii, dens)) + "\n")
    for name in ("Te.dat", "Ti.dat"):
        (profiles / name).write_text(
            "\n".join(f"{r:.8e} {ENERGY_EV:.8e}" for r in radii) + "\n")
    (profiles / "Er.dat").write_text("\n".join(f"{r:.6e} 0.0" for r in radii) + "\n")
    (case_dir / "seed.inp").write_text(
        "33\n" + " ".join(str(2000 + 11 * k) for k in range(33)) + "\n")
    set_nml(case_dir / "rmp_response_currents.inp",
            {"n_particles": n_particles, "time_step": f"{TRACE:.6e}",
             "energy_eV": f"{ENERGY_EV:.6f}", "energy_dist_kind": "'mono'",
             "boole_antithetic_variate": ".false.",
             "boole_collisions": ".true.", "boole_nonlinear_weight": ".false."})
    set_nml(case_dir / "gorilla.inp", {"boole_time_hamiltonian": ".true."})
    set_nml(case_dir / "tetra_grid.inp", {"n_field_periods_manual": n_periods})
    with open(case_dir / "rmp_response_currents.inp", "a") as f:
        f.write(f"""
&canonical_ntv_nml
  boole_canonical_ntv = .true. ,
  ntv_m_pol = 3 ,
  ntv_n_tor = {n_tor} ,
  ntv_amplitude = 1.0 ,
  ntv_bmn = 1.0e-3 ,
  ntv_bmn_relative = .true. ,
  ntv_lorentz_nu = 2.0e4 ,
  ntv_lorentz_dt = 5.0e-7 ,
  ntv_n_check = {N_CHECK} ,
  ntv_check_dt = {TRACE / N_CHECK!r} ,
  ntv_w_reset_dt = {reset_dt!r} ,
  ntv_w_reset_burn = 0.0 ,
/
""")


def run(case_dir, threads, expect_fail=False):
    env = dict(os.environ, OMP_NUM_THREADS=str(threads))
    with open(case_dir / "run.log", "w") as log:
        result = subprocess.run([str(BIN)], cwd=case_dir, env=env,
                                stdout=log, stderr=subprocess.STDOUT)
    if expect_fail:
        return result.returncode
    if result.returncode != 0:
        sys.exit(f"applet failed in {case_dir} (rc={result.returncode})")
    return None


def read_table(path, names=None):
    rows = []
    for line in path.read_text().splitlines():
        if line.startswith("# n ") and names is None:
            names = line[2:].split()
        elif line.strip() and not line.startswith("#"):
            rows.append([float(v) for v in line.split()])
    return names, rows


def markers(case_dir):
    names, rows = read_table(case_dir / "canonical_ntv_markers.dat")
    return {k: [r[i] for r in rows] for i, k in enumerate(names)}


def main():
    WORK_DIR.mkdir(parents=True, exist_ok=True)
    threads = int(os.environ.get("OMP_NUM_THREADS", "4"))
    n_particles = int(os.environ.get("NTV_N_PARTICLES", "128"))
    failures = []

    base = WORK_DIR / "lorentz_base"
    prepare(base, n_particles)
    run(base, threads)
    m = markers(base)
    kappa = m["dlnf0_dpphi_birth"]
    kbar = sum(kappa) / len(kappa)
    spread = max(abs(k - kbar) for k in kappa) / abs(kbar)
    tq_scale = max(abs(t) for t in m["tq_deltaf"])
    resid = max(abs(t + 0.5 * kbar * x * x) for t, x in zip(m["tq_deltaf"], m["dh_int"]))
    print(f"(a) kappa spread {spread:.2e}; max |tq + kappa X^2/2| / max|tq| = "
          f"{resid / tq_scale:.2e}")
    if not (tq_scale > 0.0 and math.isfinite(tq_scale)):
        failures.append("delta-f torque zero or not finite")
    if resid > max(1.0e-6, 20.0 * spread) * tq_scale:
        failures.append("tq_deltaf differs from -kappa X^2/2")
    log = (base / "run.log").read_text()
    rej = re.search(r"Canonical NTV MH collisions: proposals =\s*(\S+), rejected "
                    r"fraction =\s*(\S+)", log)
    if not rej or float(rej.group(1)) <= 0.0 or float(rej.group(2)) > 0.1:
        failures.append("Lorentz MH statistics missing or rejection above 10%")
    else:
        print(f"    Lorentz MH proposals {rej.group(1)}, rejected {rej.group(2)}")

    # (b) checkpoints
    # Checkpoints are in physical (Hamiltonian) time; a marker whose physical
    # trace time ends before the last checkpoint gets the end values there.
    _, rows = read_table(base / "canonical_ntv_checkpoints.dat", names=[])
    worst, n_end = 0.0, 0
    k = N_CHECK
    for row, x, t in zip(rows, m["dh_int"], m["tq_deltaf"]):
        if row[1] >= TRACE:
            continue
        n_end += 1
        worst = max(worst, abs(row[1 + k] - x), abs(row[1 + 2 * k] - t))
    print(f"(b) {n_end} markers end before the last checkpoint; max |diff| {worst:.2e}")
    if n_end == 0 or worst != 0.0:
        failures.append("last checkpoint differs from the end values")

    # (c) weight control
    long = WORK_DIR / "lorentz_reset_long"
    prepare(long, n_particles, reset_dt=10.0 * TRACE)
    run(long, threads)
    _, rrows = read_table(long / "canonical_ntv_reset.dat", names=[])
    ml = markers(long)
    diff = max(abs(r[1] - t) for r, t in zip(rrows, ml["tq_deltaf"]))
    print(f"(c) reset beyond trace: max |tq_settled - tq_deltaf| = {diff:.2e}")
    if diff > 1.0e-12 * tq_scale:
        failures.append("settled torque differs from tq_deltaf without resets")
    short = WORK_DIR / "lorentz_reset_quarter"
    prepare(short, n_particles, reset_dt=0.25 * TRACE)
    run(short, threads)
    ms = markers(short)
    rms0 = math.sqrt(sum(w * w for w in m["w_end"]) / len(m["w_end"]))
    rms1 = math.sqrt(sum(w * w for w in ms["w_end"]) / len(ms["w_end"]))
    print(f"    rms w_end without / with resets: {rms0:.3e} / {rms1:.3e}")
    if not rms1 < rms0:
        failures.append("weight resets do not bound the weights")

    # (d) harmonic not periodic on the mesh period
    bad = WORK_DIR / "lorentz_bad_period"
    prepare(bad, 8, n_tor=3, n_periods=2)
    rc = run(bad, threads, expect_fail=True)
    print(f"(d) n_tor = 3 on a 2-period mesh: rc = {rc}")
    if rc == 0 or "multiple of n_field_periods" not in (bad / "run.log").read_text():
        failures.append("non-periodic harmonic was not rejected")

    if failures:
        print("FAIL: " + "; ".join(failures))
        return 1
    print("PASS: canonical NTV Lorentz-mode gates")
    return 0


if __name__ == "__main__":
    sys.exit(main())
