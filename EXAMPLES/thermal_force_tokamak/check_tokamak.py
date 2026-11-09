#!/usr/bin/env python3
"""
Checks of the tokamak driver of the thermal-force applet against analytic
results and against the slab driver (uniform field), W in H+, 50 eV.

  orbits     collisionless (n = 1e6 cm^-3 and dt = 1e-6 tau_s, so the one
             collision at t = 0 is negligible): vertical guiding-centre velocity
             (<Z(t)> - <Z(0)>)/t against grad-B plus curvature drift
             c m (<v_par^2> + <v_perp^2>/2)/(q B0 R0) plus the vertical part of
             the parallel streaming, <v_par> B_Z/B_phi = <v_par> r/(q R0) at the
             outboard midplane; the moments are those of the sampled ions
             (energy, and the spread of R_s dphi = v_par t). The drift depends
             on m/q, i.e. on the species override.
  screening  drift d<R>/dt (paired difference, Omega_a t = 300..700) for
             W1+, W3+, W10+ at 100 eV/m against our first-order targets and
             the slab driver, with the collision shift on and off
  steps      the same drift versus collision step, tokamak and slab
  diffusion  spread of R and Z (control run, grad T = 0) versus time against the
             slab's x and y: classical diffusion from the guiding-centre shift

Usage: python3 check_tokamak.py [work_dir]   (default ./run); writes
<work_dir>/tokamak_checks.png and prints the numbers.
"""
import csv
import math
import os
import shutil
import subprocess
import sys
from pathlib import Path

import f90nml
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
BINARY = Path(os.environ.get("GORILLA_APPLETS_BIN", ROOT / "BUILD" / "gorilla_applets_main.x"))
ENV = dict(os.environ, OMP_NUM_THREADS=os.environ.get("OMP_NUM_THREADS", "8"))
TOKAMAK_FILES = ("gorilla_applets.inp", "gorilla.inp", "tetra_grid.inp", "field_divB0.inp", "preload_for_SYNCH.inp")
SLAB_FILES = ("gorilla_applets.inp", "gorilla.inp", "tetra_grid.inp")

C, ECHARGE, EV2ERG, AMU = 2.9979e10, 4.8032e-10, 1.6022e-12, 1.66053906660e-24   # GORILLA constants, CODATA amu
R_S, B0, R0 = 100050.0, 1.0e4, 1.0e5

BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, MUTED, GRID, SURFACE = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"


def run(case_dir, template, files, overrides):
    if case_dir.exists():
        shutil.rmtree(case_dir)
    case_dir.mkdir(parents=True)
    for f in files:
        if (template / f).exists():
            shutil.copy(template / f, case_dir / f)
    nml = f90nml.read(template / "thermal_force.inp")
    nml["thermal_force_nml"].update({"output_file": "out.csv", **overrides})
    nml.write(case_dir / "thermal_force.inp")
    res = subprocess.run([str(BINARY)], cwd=case_dir, capture_output=True, text=True, env=ENV)
    if res.returncode != 0:
        print(res.stdout[-2000:], res.stderr[-2000:])
        sys.exit(f"applet failed in {case_dir}")
    with open(case_dir / "out.csv") as f:
        rows = [{k: (v.strip() if k == "run" else float(v)) for k, v in r.items()}
                for r in csv.DictReader(line for line in f if not line.startswith("#"))]
    shutil.rmtree(case_dir)
    return rows


def tokamak(work, **kw):
    return run(work / "tmp", HERE, TOKAMAK_FILES, kw)


def slab(work, **kw):
    base = {"driver": "slab", "orbit_model": "guiding_centre", "test_distribution": "maxwellian",
            "test_temperature_ev": 50.0, "test_drift_cm_s": [0.0, 0.0, 0.0], "start_x_cm": [0.0, 0.0, 0.0],
            "b_gauss": [0.0, 0.0, B0]}
    return run(work / "tmp", ROOT / "EXAMPLES" / "thermal_force", SLAB_FILES, {**base, **kw})


def drift(rows, key="mean_x", run_name="difference"):
    """Slope between the last two output times (cm -> m per s) and its bound on the standard error."""
    r = [x for x in rows if x["run"] == run_name]
    a, b = r[-2], r[-1]
    dt = b["t_s"] - a["t_s"]
    se_key = "se_" + key.split("_")[1]
    return (b[key] - a[key]) / 100.0 / dt, (b[se_key] + a[se_key]) / 100.0 / dt


def targets():
    t = {}
    with open(ROOT / "TESTS" / "REFERENCE" / "homma" / "B07_targets.csv") as f:
        for r in csv.DictReader(line for line in f if not line.startswith("#")):
            if float(r["gradT_x_eV_m"]) == 100.0:
                t[int(r["Z"])] = float(r["gyro_orbit_dXdt_m_s"])
    return t


def style(ax, title):
    ax.set_facecolor(SURFACE)
    ax.grid(color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=8.5)
    ax.set_title(title, fontsize=10, color=INK, loc="left")


def main():
    work = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE / "run"
    work.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 9), facecolor=SURFACE)
    times = [0.0] + [50.0 * k for k in range(1, 15)]           # Omega_a t = 0, 50, ..., 700

    # 1. collisionless orbits: species override and pusher
    ax = axes[0, 0]
    n = 20000
    t_energy = 50.0 * EV2ERG
    r_minor, q_safety = R_S - R0, 10.0
    for z, color in ((1, BLUE), (3, ORANGE), (10, AQUA)):
        rows = tokamak(work, test_charge_number=float(z), density_cm3=1.0e6, dt_coll_over_tau_s=1.0e-6,
                       grad_temperature_ev_cm=[0.0, 0.0, 0.0], boole_control_run=False, n_particles=n,
                       output_times=times)
        rows = [r for r in rows if r["run"] == "gradient"]
        r0, rows = rows[0], rows[1:]
        m, q = 183.84 * AMU, z * ECHARGE
        ratio, stream = [], []
        for r in rows:
            t = r["t_s"]
            vz = (r["mean_z"] - r0["mean_z"]) / t                                   # cm/s
            vpar = (r["mean_y"] - r0["mean_y"]) / t
            vpar2 = (r["se_y"] * math.sqrt(n) / t) ** 2 + vpar ** 2               # <v_par^2>
            v2 = 2.0 * r["mean_energy"] / m                                       # <v^2>
            expect = C * m * (vpar2 + 0.5 * (v2 - vpar2)) / (q * B0 * R0) + vpar * r_minor / (q_safety * R0)
            ratio.append(vz / expect)
            stream.append(math.sqrt(vpar2) / math.sqrt(t_energy / m))
        wt = [r["omega_t"] for r in rows]
        ax.plot(wt, ratio, "-o", color=color, markersize=3, linewidth=1.4, label=f"W{z}+: vertical drift")
        ax.plot(wt, stream, "--", color=color, linewidth=1.2, label=f"W{z}+: sqrt(<v_par^2>) ÷ sqrt(T/m)")
        print(f"orbits W{z}+: vertical velocity / (drift + streaming) = {np.mean(ratio[-5:]):.4f}; "
              f"parallel streaming sqrt(<v_par^2>)/sqrt(T/m) = {stream[-1]:.4f}")
    ax.axhline(1.0, color=INK, linewidth=0.8)
    ax.set_ylim(0.95, 1.03)
    ax.set_xlabel("Ω_a t")
    ax.set_ylabel("measured / expected")
    ax.legend(frameon=False, fontsize=7, labelcolor=INK, ncol=2, loc="lower right")
    style(ax, "1. Collisionless orbits: ∇B+curvature drift, parallel streaming")

    # 2. screening drift for three charge states, shift on and off
    ax = axes[0, 1]
    tg = targets()
    common = {"output_time_unit": "gyro", "output_times": [0.0, 300.0, 700.0], "dt_coll_over_tau_s": 0.01,
              "grad_temperature_ev_cm": [1.0, 0.0, 0.0], "boole_control_run": True}
    for i, z in enumerate((1, 3, 10)):
        n = 6000
        tok_on = drift(tokamak(work, test_charge_number=float(z), n_particles=n, **common))
        tok_off = drift(tokamak(work, test_charge_number=float(z), n_particles=n, boole_collision_shift=False,
                                **common))
        sl = drift(slab(work, test_charge_number=float(z), n_particles=n, **common))
        for k, (label, (v, s), color, marker) in enumerate((("tokamak, shift on", tok_on, BLUE, "s"),
                                                             ("slab, shift on", sl, AQUA, "^"),
                                                             ("tokamak, shift off", tok_off, ORANGE, "o"))):
            ax.errorbar(i + 0.15 * (k - 1), v / tg[z], yerr=s / abs(tg[z]), color=color, marker=marker,
                        markersize=6, capsize=2, linestyle="none", label=label if i == 0 else None)
        print(f"screening W{z}+: target {tg[z]:.4f}, tokamak {tok_on[0]:.4f} +- {tok_on[1]:.4f}, "
              f"slab {sl[0]:.4f} +- {sl[1]:.4f}, tokamak without shift {tok_off[0]:.4f} +- {tok_off[1]:.4f} m/s")
    ax.axhline(1.0, color=INK, linestyle="--", linewidth=1.0, label="first-order target (collisions at X+ρ)")
    ax.axhline(0.0, color=MUTED, linewidth=0.8)
    ax.set_xticks([0, 1, 2], ["W1+", "W3+", "W10+"])
    ax.set_ylabel("d⟨R⟩/dt ÷ target")
    ax.legend(frameon=False, fontsize=7.5, labelcolor=INK, loc="center right")
    style(ax, "2. Screening drift at 100 eV/m, Ω_a t = 300–700")

    # 3. step-size dependence, tokamak versus slab, W3+
    ax = axes[1, 0]
    for name, fn, color, marker in (("tokamak", tokamak, BLUE, "s"), ("slab", slab, AQUA, "^")):
        dts, vs, ss = [], [], []
        for dt in (0.005, 0.01, 0.02, 0.05, 0.1, 0.15):
            v, s = drift(fn(work, **{**common, "dt_coll_over_tau_s": dt}, n_particles=6000))
            dts.append(dt); vs.append(v); ss.append(s)
        ax.errorbar(dts, (np.array(vs) / tg[3] - 1) * 100, yerr=np.array(ss) / abs(tg[3]) * 100, color=color,
                    marker=marker, markersize=5, linewidth=1.4, capsize=2, label=name)
        print(f"steps {name}: " + ", ".join(f"{d}: {(v/tg[3]-1)*100:+.2f}%" for d, v in zip(dts, vs)))
    ax.axhline(0.0, color=INK, linestyle="--", linewidth=1.0)
    ax.set_xscale("log")
    ax.set_xlabel("collision step Δt / τ_s")
    ax.set_ylabel("d⟨R⟩/dt relative to target [%]")
    ax.legend(frameon=False, fontsize=8, labelcolor=INK)
    style(ax, "3. Step-size dependence, W3+, 100 eV/m")

    # 4. classical diffusion from the guiding-centre shift (control run, grad T = 0)
    ax = axes[1, 1]
    n = 20000
    common_d = {"output_time_unit": "gyro", "output_times": times, "dt_coll_over_tau_s": 0.01,
                "grad_temperature_ev_cm": [0.0, 0.0, 0.0], "boole_control_run": False, "n_particles": n}
    tok = [r for r in tokamak(work, **common_d) if r["run"] == "gradient"]
    sl = [r for r in slab(work, **common_d) if r["run"] == "gradient"]
    wt = [r["omega_t"] for r in tok]
    var = lambda rows, key: np.array([(r[key] * math.sqrt(n)) ** 2 / 1e4 * 1e6 for r in rows])   # mm^2
    ax.plot(wt, var(tok, "se_x"), "-", color=BLUE, linewidth=1.6, label="tokamak: var R")
    ax.plot(wt, var(tok, "se_z"), "-", color=ORANGE, linewidth=1.6, label="tokamak: var Z")
    ax.plot(wt, var(sl, "se_x"), "o", color=BLUE, markersize=4, label="slab: var x")
    ax.plot(wt, var(sl, "se_y"), "o", color=ORANGE, markersize=4, label="slab: var y")
    print(f"diffusion at Omega t = 700: var R {var(tok, 'se_x')[-1]:.4f}, var Z {var(tok, 'se_z')[-1]:.4f} mm^2; "
          f"slab var x {var(sl, 'se_x')[-1]:.4f}, var y {var(sl, 'se_y')[-1]:.4f} mm^2")
    ax.set_xlabel("Ω_a t")
    ax.set_ylabel("variance of the guiding-centre position [mm²]")
    ax.legend(frameon=False, fontsize=8, labelcolor=INK)
    style(ax, "4. Classical diffusion from the guiding-centre shift, ∇T = 0")

    fig.suptitle("Tokamak driver (GORILLA analytic circular tokamak, R₀ = 1 km) versus analytic results "
                 "and the slab driver", fontsize=11, color=INK)
    fig.tight_layout()
    out = work / "tokamak_checks.png"
    fig.savefig(out, dpi=140)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
