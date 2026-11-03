#!/usr/bin/env python3
"""Stationarity of an unperturbed Maxwellian under orbits and collisions.

Runs the RMP applet (i_option 15) in the axisymmetric analytic circular
tokamak (grid_kind 5) without perturbation and with flat density and
temperature, so that the local Maxwellian is an exact equilibrium of the
guiding-centre motion and of the mode-5 OU collision operator with
Metropolis-Hastings correction. Markers are loaded with the nonlinear
loader (uniform physical volume, untruncated Maxwellian, antithetic pairs).

Checks, with independent analytic oracles and one marker per antithetic pair:
  * the loaded and the final markers follow the analytic window density in
    rho^2 and geometric poloidal angle, E/T ~ Gamma(3/2) and uniform pitch
    (chi-square tests inside a band inset from the window edges), after
    more than eight collision times;
  * no marker is lost and the nonlinear weights stay at round-off;
  * a fixed external seed.inp reproduces the run bit for bit with several
    OpenMP threads, and a different seed.inp changes it.

Paths come from CMake through APPLETS_ROOT, GORILLA_APPLETS_BIN and WORK_DIR.
FM_N_PARTICLES and FM_TIME_STEP override the main case size (cluster runs).
"""

import math
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

ME = 9.1094e-28  # g, GORILLA constants_mod
EV2ERG = 1.6022e-12
TEMPERATURE_EV = 1000.0
DENSITY = 1.2e14
P_MIN = 1.0e-4
MIN_COLLISION_TIMES = 8.0
RHO_MARGIN = 1.5  # cm


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


def seed_values(seed):
    state = seed
    values = []
    for _ in range(33):
        state = (6364136223846793005 * state + 1442695040888963407) % 2**64
        values.append(1 + (state >> 33) % 2147483646)
    return values


def prepare(case_dir, seed, n_particles, time_step):
    if case_dir.exists():
        shutil.rmtree(case_dir)
    shutil.copytree(TEMPLATE, case_dir)
    profiles = case_dir / "profiles"
    profiles.mkdir()
    radii = [0.1 + 0.25 * i for i in range(300)]
    for name, value in (("n.dat", DENSITY), ("Te.dat", TEMPERATURE_EV),
                        ("Ti.dat", TEMPERATURE_EV), ("Er.dat", 0.0)):
        lines = [f"{r:.6e} {value:.6e}" for r in radii]
        (profiles / name).write_text("\n".join(lines) + "\n")
    values = seed_values(seed)
    (case_dir / "seed.inp").write_text(
        "33\n" + " ".join(str(v) for v in values) + "\n")
    set_nml(case_dir / "rmp_response_currents.inp",
            {"n_particles": n_particles, "time_step": time_step})


def run(case_dir, threads):
    env = dict(os.environ, OMP_NUM_THREADS=str(threads))
    with open(case_dir / "run.log", "w") as log:
        result = subprocess.run([str(BIN)], cwd=case_dir, env=env,
                                stdout=log, stderr=subprocess.STDOUT)
    if result.returncode != 0:
        sys.exit(f"applet failed in {case_dir} (rc={result.returncode})")
    return (case_dir / "run.log").read_text()


def log_value(log, pattern):
    match = re.search(pattern + r"\s*([-+0-9.Ee]+)", log)
    if not match:
        sys.exit(f"log entry not found: {pattern}")
    return float(match.group(1))


def gammq(a, x):
    """Regularized upper incomplete gamma function Q(a, x)."""
    if x <= 0.0:
        return 1.0
    gln = math.lgamma(a)
    if x < a + 1.0:
        term = total = 1.0 / a
        ap = a
        for _ in range(10000):
            ap += 1.0
            term *= x / ap
            total += term
            if abs(term) < abs(total) * 1e-15:
                break
        return 1.0 - total * math.exp(-x + a * math.log(x) - gln)
    b = x + 1.0 - a
    c = 1.0 / 1e-300
    d = 1.0 / b
    h = d
    for i in range(1, 10000):
        an = -i * (i - a)
        b += 2.0
        d = an * d + b
        d = 1e-300 if abs(d) < 1e-300 else d
        c = b + an / c
        c = 1e-300 if abs(c) < 1e-300 else c
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < 1e-15:
            break
    return math.exp(-x + a * math.log(x) - gln) * h


def chi2_pvalue(samples, cdf, edges):
    """Chi-square goodness of fit of samples against cdf on bin edges."""
    counts = [0] * (len(edges) - 1)
    outside = 0
    for value in samples:
        for k in range(len(edges) - 1):
            if edges[k] <= value < edges[k + 1]:
                counts[k] += 1
                break
        else:
            outside += 1
    probs = [cdf(edges[k + 1]) - cdf(edges[k]) for k in range(len(edges) - 1)]
    n = len(samples)
    expected = [n * p for p in probs] + [n * (1.0 - sum(probs))]
    observed = counts + [outside]
    chi2 = 0.0
    dof = -1
    for obs, exp in zip(observed, expected):
        if exp > 1e-12:
            chi2 += (obs - exp) ** 2 / exp
            dof += 1
        elif obs > 0:
            return chi2, 0.0
    return chi2, gammq(0.5 * dof, 0.5 * chi2)


def maxwell_energy_cdf(x):
    """CDF of E/T for a 3D Maxwellian, Gamma(3/2, 1)."""
    if x <= 0.0:
        return 0.0
    if math.isinf(x):
        return 1.0
    return math.erf(math.sqrt(x)) - 2.0 * math.sqrt(x / math.pi) * math.exp(-x)


def window_geometry(case_dir):
    grid = (case_dir / "tetra_grid.inp").read_text()
    rmp = (case_dir / "rmp_response_currents.inp").read_text()

    def number(text, key):
        match = re.search(rf"^\s*{key}\s*=\s*([-+0-9.eEdD]+)", text, flags=re.M)
        return float(match.group(1).lower().replace("d", "e"))

    r0 = number(grid, "R0_analytic_circ")
    a = number(grid, "a_analytic_circ")
    s_lo = number(rmp, "s_inner_sample")
    s_hi = number(rmp, "s_outer_sample")
    edge = a * a / (r0 + math.sqrt(r0 * r0 - a * a))
    rho2_lo = s_lo * edge * (2.0 * r0 - s_lo * edge)
    rho2_hi = s_hi * edge * (2.0 * r0 - s_hi * edge)
    return r0, rho2_lo, rho2_hi


def read_loaded(case_dir):
    markers = {}
    for line in (case_dir / "loaded_maxwellian.dat").read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        f = line.split()
        markers[int(f[0])] = (float(f[1]), float(f[3]), float(f[5]), float(f[6]))
    return markers


def read_final(case_dir):
    markers = {}
    for line in (case_dir / "exit_data.dat").read_text().splitlines():
        f = line.split()
        if not f:
            continue
        vpar, vperp = float(f[6]), float(f[7])
        v2 = vpar * vpar + vperp * vperp
        energy = 0.5 * ME * v2 / (TEMPERATURE_EV * EV2ERG)
        pitch = vpar / math.sqrt(v2) if v2 > 0.0 else 0.0
        markers[int(f[0])] = (float(f[3]), float(f[5]), energy, pitch,
                              int(f[1]))
    return markers


def check_distribution(label, markers, geometry, failures):
    r0, rho2_lo, rho2_hi = geometry
    # The window is a geometric annulus, but the invariant surfaces of the
    # piecewise-linear mesh field are polygons (sagitta ~0.5 cm here), so the
    # density next to the window edges relaxes. Test inside a band inset by
    # RHO_MARGIN, which exceeds the largest observed radial excursion.
    band_lo = (math.sqrt(rho2_lo) + RHO_MARGIN) ** 2
    band_hi = (math.sqrt(rho2_hi) - RHO_MARGIN) ** 2
    odd = [markers[n] for n in sorted(markers) if n % 2 == 1]
    rho2_all = [(m[0] - r0) ** 2 + m[1] ** 2 for m in odd]
    band = [m for m, r2 in zip(odd, rho2_all) if band_lo <= r2 < band_hi]
    rho2 = [(m[0] - r0) ** 2 + m[1] ** 2 for m in band]
    theta = [math.atan2(m[1], m[0] - r0) for m in band]
    energy = [m[2] for m in band]
    pitch = [m[3] for m in band]

    nb = 10
    rho2_edges = [band_lo + (band_hi - band_lo) * k / nb for k in range(nb + 1)]

    def rho2_cdf(x):
        # Uniform physical volume makes rho^2 uniform over the whole window;
        # markers outside the band form the remaining chi-square cell.
        return min(max((x - rho2_lo) / (rho2_hi - rho2_lo), 0.0), 1.0)

    # Density in (rho^2, theta) proportional to R = r0 + rho cos(theta).
    rho3 = band_hi ** 1.5 - band_lo ** 1.5
    c_amp = (2.0 / 3.0) * rho3 / (r0 * (band_hi - band_lo))

    def theta_cdf(x):
        return ((x + math.pi) + c_amp * math.sin(x)) / (2.0 * math.pi)

    theta_edges = [-math.pi + 2.0 * math.pi * k / 12 for k in range(13)]
    energy_edges = [0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 4.0,
                    5.0, 7.0, math.inf]
    pitch_edges = [-1.0 + 0.2 * k for k in range(11)]
    pitch_edges[-1] = 1.0 + 1e-12

    def pitch_cdf(x):
        return min(max(0.5 * (x + 1.0), 0.0), 1.0)

    tests = (("rho^2", rho2_all, rho2_cdf, rho2_edges),
             ("theta", theta, theta_cdf, theta_edges),
             ("E/T", energy, maxwell_energy_cdf, energy_edges),
             ("pitch", pitch, pitch_cdf, pitch_edges))
    for name, samples, cdf, edges in tests:
        chi2, p = chi2_pvalue(samples, cdf, edges)
        status = "ok" if p > P_MIN else "FAIL"
        print(f"{label:8s} {name:6s} N={len(samples):6d} chi2={chi2:8.2f} "
              f"p={p:.3e} {status}")
        if p <= P_MIN:
            failures.append(f"{label} {name} p={p:.3e}")
    mean_energy = sum(energy) / len(energy)
    print(f"{label:8s} <E/T>={mean_energy:.4f} (expected 1.5, "
          f"standard error {math.sqrt(1.5 / len(energy)):.4f})")


def main():
    WORK_DIR.mkdir(parents=True, exist_ok=True)
    threads = int(os.environ.get("OMP_NUM_THREADS", "4"))
    n_particles = int(os.environ.get("FM_N_PARTICLES", "2048"))
    time_step = os.environ.get("FM_TIME_STEP", "2.0d-5")
    failures = []

    case = WORK_DIR / "stationarity"
    prepare(case, 20261019, n_particles, time_step)
    log = run(case, threads)
    collision_steps = log_value(log, r"average number of integration steps =")
    ou_dtau = 0.1
    print(f"collision times per marker ~ {collision_steps * ou_dtau:.2f}")
    if collision_steps * ou_dtau < MIN_COLLISION_TIMES:
        failures.append("fewer than eight collision times")
    max_w = log_value(log, r"max \|w\| =")
    print(f"max |w| = {max_w:.3e}")
    if max_w > 1e-10:
        failures.append(f"nonlinear weights drift, max |w| = {max_w:.3e}")
    match = re.search(r"Lost markers: before burn-in = (\d+), after = (\d+)", log)
    if not match or match.group(1) != "0" or match.group(2) != "0":
        failures.append("markers lost in the unperturbed field")

    geometry = window_geometry(case)
    loaded = read_loaded(case)
    final = read_final(case)
    if len(loaded) != n_particles or len(final) != n_particles:
        failures.append("marker files incomplete")
    if any(m[4] != 0 for m in final.values()):
        failures.append("exit data reports lost markers")
    check_distribution("loaded", loaded, geometry, failures)
    check_distribution("final", final, geometry, failures)

    # External seed: bit reproducibility with several threads, and sensitivity.
    outputs = []
    for tag, seed in (("seed_a1", 11), ("seed_a2", 11), ("seed_b", 12)):
        small = WORK_DIR / tag
        prepare(small, seed, 64, "2.0d-6")
        run(small, 2)
        outputs.append(((small / "loaded_maxwellian.dat").read_bytes(),
                        (small / "exit_data.dat").read_bytes()))
    if outputs[0] != outputs[1]:
        failures.append("same seed.inp does not reproduce the run")
    if outputs[0][0] == outputs[2][0] or outputs[0][1] == outputs[2][1]:
        failures.append("different seed.inp does not change the run")
    print("seed reproducibility:", "ok" if outputs[0] == outputs[1] else "FAIL")
    print("seed sensitivity:", "ok" if outputs[0][1] != outputs[2][1] else "FAIL")

    if failures:
        print("FAIL: " + "; ".join(failures))
        return 1
    print("PASS: unperturbed Maxwellian stationary with collisions; seed honoured")
    return 0


if __name__ == "__main__":
    sys.exit(main())
