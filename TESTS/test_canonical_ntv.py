#!/usr/bin/env python3
"""Gates of the canonical delta-f / full-f p_phi torque estimator.

Runs the RMP applet (i_option 15) in the analytic circular tokamak
(grid_kind 5) with the canonical f0(p_phi, H0) loader and the canonical NTV
estimator (canonical_ntv_nml), density and temperature gradients, and one
harmonic dB1 = sigma b_mn cos(m theta - n phi) evaluated analytically along
the unperturbed mesh-field orbits.

Gates:
  (a) sigma = 0 without collisions: the delta-f weights stay exactly zero;
      the p_phi drift along the orbits on the tetrahedral mesh is reported
      as an equivalent radial shift of the canonical flux label s* and must
      stay below DRIFT_SSTAR_MAX (piecewise-linear field, no exact
      toroidal symmetry inside a cell);
  (c) with MH collisions, common random numbers and amplitudes +sigma,
      -sigma, 2 sigma: the delta-f torque sum P int w (-dH1/dphi) dt is
      nonzero, even in sigma (odd part vanishes) and quadratic
      (T(2 sigma) = 4 T(sigma)); int -dH1/dphi dt is odd; the orbit and
      collisional p_phi tallies are independent of sigma.

Paths come from CMake through APPLETS_ROOT, GORILLA_APPLETS_BIN and WORK_DIR.
NTV_N_PARTICLES and NTV_TIME_STEP override the case size.
"""

import math
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

DENSITY = 1.2e14
TEMPERATURE_EV = 1000.0
A_MINOR = 69.991288
DRIFT_SSTAR_MAX = 2.0e-3
REL_TOL = 1.0e-9


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


def prepare(case_dir, n_particles, time_step, sigma, collisions):
    if case_dir.exists():
        shutil.rmtree(case_dir)
    shutil.copytree(TEMPLATE, case_dir)
    profiles = case_dir / "profiles"
    profiles.mkdir()
    radii = [0.1 + 0.25 * i for i in range(300)]
    for name, value in (("n.dat", DENSITY), ("Te.dat", TEMPERATURE_EV),
                        ("Ti.dat", TEMPERATURE_EV)):
        # Gradients in n and T so that d ln f0/dp_phi is nonzero.
        lines = [f"{r:.6e} {value * (1.0 - 0.6 * (r / A_MINOR) ** 2):.6e}"
                 for r in radii]
        (profiles / name).write_text("\n".join(lines) + "\n")
    (profiles / "Er.dat").write_text(
        "\n".join(f"{r:.6e} 0.0" for r in radii) + "\n")
    (case_dir / "seed.inp").write_text(
        "33\n" + " ".join(str(1000 + 7 * k) for k in range(33)) + "\n")
    set_nml(case_dir / "rmp_response_currents.inp",
            {"n_particles": n_particles, "time_step": time_step,
             "boole_collisions": ".true." if collisions else ".false.",
             "boole_nonlinear_weight": ".false."})
    set_nml(case_dir / "gorilla.inp", {"boole_time_hamiltonian": ".true."})
    with open(case_dir / "rmp_response_currents.inp", "a") as f:
        f.write(f"""
&canonical_ntv_nml
  boole_canonical_ntv = .true. ,
  ntv_m_pol = 3 ,
  ntv_n_tor = 2 ,
  ntv_amplitude = {sigma!r} ,
  ntv_bmn = 10.0 ,
  ntv_phimn = 0.0 ,
  ntv_theta_kind = 1 ,
/
""")


def run(case_dir, threads):
    env = dict(os.environ, OMP_NUM_THREADS=str(threads))
    with open(case_dir / "run.log", "w") as log:
        result = subprocess.run([str(BIN)], cwd=case_dir, env=env,
                                stdout=log, stderr=subprocess.STDOUT)
    if result.returncode != 0:
        sys.exit(f"applet failed in {case_dir} (rc={result.returncode})")
    return read_markers(case_dir / "canonical_ntv_markers.dat")


def read_markers(path):
    names = None
    rows = []
    for line in path.read_text().splitlines():
        if line.startswith("# n "):
            names = line[2:].split()
        elif line.strip() and not line.startswith("#"):
            rows.append([float(v) for v in line.split()])
    if names is None or not rows:
        sys.exit(f"{path}: no marker table")
    return {k: [r[i] for r in rows] for i, k in enumerate(names)}


def psum(m, key):
    return sum(p * v for p, v in zip(m["P"], m[key]))


def close(a, b, scale):
    return abs(a - b) <= REL_TOL * scale


def main():
    WORK_DIR.mkdir(parents=True, exist_ok=True)
    threads = int(os.environ.get("OMP_NUM_THREADS", "4"))
    n_particles = int(os.environ.get("NTV_N_PARTICLES", "256"))
    time_step = os.environ.get("NTV_TIME_STEP", "2.0d-6")
    failures = []

    # (a) unperturbed, collisionless
    m0 = run_case("sigma0_nocoll", n_particles, time_step, 0.0, False,
                  threads)
    max_w = max(abs(w) for w in m0["w_end"])
    dsstar = [abs(b - a) for a, b, lost in
              zip(m0["sstar_birth"], m0["sstar_end"], m0["lost"]) if lost == 0]
    rms = math.sqrt(sum(d * d for d in dsstar) / len(dsstar))
    print(f"(a) max |w| = {max_w:.3e}; p_phi drift as |Delta s*|: "
          f"max {max(dsstar):.3e}, rms {rms:.3e} over {len(dsstar)} markers")
    if max_w != 0.0:
        failures.append(f"unperturbed weights nonzero ({max_w:.3e})")
    if max(dsstar) > DRIFT_SSTAR_MAX:
        failures.append(f"p_phi drift on the mesh too large ({max(dsstar):.3e})")
    if any(abs(c) > 0.0 for c in m0["dpphi_coll"]):
        failures.append("collisional p_phi tally nonzero without collisions")

    # (c) amplitude scaling with common random numbers and collisions
    sig = 1.0
    runs = {s: run_case(f"sigma_{s:+.1f}", n_particles, time_step, s, True,
                        threads) for s in (sig, -sig, 2.0 * sig)}
    t1, tm, t2 = (psum(runs[s], "tq_deltaf") for s in (sig, -sig, 2 * sig))
    h1, hm = (psum(runs[s], "dh_int") for s in (sig, -sig))
    scale = abs(t1)
    print(f"(c) T(+s) = {t1:.10e}, T(-s) = {tm:.10e}, T(2s) = {t2:.10e}")
    print(f"    odd part (T(+s)-T(-s))/2 / T = {(t1 - tm) / 2 / scale:.3e}; "
          f"T(2s)/T(s) = {t2 / t1:.12f}")
    print(f"    sum P int -dH1/dphi dt: +s {h1:.6e}, -s {hm:.6e}")
    if not (t1 != 0.0 and math.isfinite(t1)):
        failures.append("delta-f torque zero or not finite")
    if not close(t1, tm, scale):
        failures.append("delta-f torque has an odd part")
    if not close(t2, 4.0 * t1, 4.0 * scale):
        failures.append("delta-f torque not quadratic in sigma")
    if not close(h1, -hm, abs(h1)):
        failures.append("int -dH1/dphi dt not odd in sigma")
    for key in ("dpphi_orbit", "dpphi_coll", "pphi_end"):
        if runs[sig][key] != runs[-sig][key] or runs[sig][key] != runs[2 * sig][key]:
            failures.append(f"{key} depends on sigma (orbits must be unperturbed)")
    ncoll = sum(1 for c in runs[sig]["dpphi_coll"] if c != 0.0)
    print(f"    markers with collisional p_phi change: {ncoll}")
    if ncoll == 0:
        failures.append("no collisional p_phi change recorded")

    if failures:
        print("FAIL: " + "; ".join(failures))
        return 1
    print("PASS: canonical NTV estimator gates")
    return 0


def run_case(tag, n_particles, time_step, sigma, collisions, threads):
    case = WORK_DIR / tag
    prepare(case, n_particles, time_step, sigma, collisions)
    return run(case, threads)


if __name__ == "__main__":
    sys.exit(main())
