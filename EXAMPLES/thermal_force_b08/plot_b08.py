#!/usr/bin/env python3
"""
Plots and convergence fits for the B08 campaign of run_b08.py.

  b08_dt.png       B07 drift versus collision step, relative to our first-order
                   targets (TESTS/REFERENCE/homma/B07_targets.csv), and the
                   transverse drift; error bars from the seed-to-seed scatter
  b08_relax.png    B04 moments at 1, 3, 10 tau_s versus collision step
  b08_n.png        standard error of the drift versus ion number: reported by
                   the applet and from the seed-to-seed scatter
  b08_summary.txt  first-order fits d(dt) = d0 + c dt of the drift over
                   dt <= 0.05 tau_s (full orbit: also Omega_a dt <= 1)

Usage: python3 plot_b08.py [work_dir]   (default ./run)
"""
import csv
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
TARGETS = HERE.parent.parent / "TESTS" / "REFERENCE" / "homma" / "B07_targets.csv"

# reference palette (dataviz skill), light mode
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, MUTED, GRID, SURFACE = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
MODELS = {"gc_gc": ("guiding centre, background at X", ORANGE, "o", "gc_dXdt_m_s"),
          "gc_orbit": ("guiding centre, background at X+ρ", BLUE, "s", "gyro_orbit_dXdt_m_s"),
          "full": ("full orbit (reference)", AQUA, "^", "gyro_orbit_dXdt_m_s")}


def read(path):
    with open(path) as f:
        return [{k: (v if k == "model" else float(v)) for k, v in r.items()} for r in csv.DictReader(f)]


def targets():
    t = {}
    with open(TARGETS) as f:
        for r in csv.DictReader(line for line in f if not line.startswith("#")):
            if float(r["gradT_x_eV_m"]) == 100.0:
                t[int(r["Z"])] = {k: float(v) for k, v in r.items()}
    return t


def style(ax):
    ax.set_facecolor(SURFACE)
    ax.grid(color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=8.5)


def seed_stats(values):
    v = np.array(values)
    return v.mean(), v.std(ddof=1) / np.sqrt(len(v))




def plot_dt(work, summary):
    rows = read(work / "b08_dt.csv")
    tg = targets()
    groups = defaultdict(list)
    for r in rows:
        groups[(r["model"], int(r["Z"]), r["dt_over_tau_s"])].append(r)
    fig, axes = plt.subplots(2, 3, figsize=(13, 7.5), facecolor=SURFACE, sharex=True)
    for col, z in enumerate((1, 3, 10)):
        ax_x, ax_y = axes[0, col], axes[1, col]
        omega_tau = None
        for model, (label, color, marker, tcol) in MODELS.items():
            target = tg[z][tcol]
            dts = sorted({k[2] for k in groups if k[0] == model and k[1] == z})
            if not dts:
                continue
            mx, sx, my, sy = [], [], [], []
            for dt in dts:
                g = groups[(model, z, dt)]
                omega_tau = g[0]["omega_tau_s"]
                m, s = seed_stats([r["dXdt_m_s"] for r in g])
                mx.append(m); sx.append(s)
                m, s = seed_stats([r["dYdt_m_s"] for r in g])
                my.append(m); sy.append(s)
            dts, mx, sx, my, sy = map(np.array, (dts, mx, sx, my, sy))
            scale = 100.0 / abs(target)
            ax_x.errorbar(dts, (mx / target - 1.0) * 100.0, yerr=sx * scale, color=color, marker=marker,
                          markersize=5, linewidth=1.4, capsize=2, label=label)
            ax_y.errorbar(dts, my * scale, yerr=sy * scale, color=color, marker=marker, markersize=5,
                          linewidth=1.4, capsize=2, label=label)
            # first-order fit d0 + c dt (Euler-Maruyama); full orbit only where Omega dt <= 1
            sel = dts <= 0.05
            if model == "full":
                sel &= dts * omega_tau <= 1.0
            if sel.sum() >= 3:
                w = 1.0 / np.maximum(sx[sel], 1e-12) ** 2
                a = np.vstack([np.ones(sel.sum()), dts[sel]]).T
                cov = np.linalg.inv(a.T @ (w[:, None] * a))
                d0, c = cov @ (a.T @ (w * mx[sel]))
                chi2 = np.sum(w * (mx[sel] - d0 - c * dts[sel]) ** 2) / max(sel.sum() - 2, 1)
                summary.append(f"Z={z:2d} {model:8s} d0 = {d0:.5f} +- {np.sqrt(cov[0, 0]):.5f} m/s "
                               f"(d0/target - 1 = {d0/target-1:+.4f} +- {np.sqrt(cov[0, 0])/abs(target):.4f}), "
                               f"bias per 0.01 tau_s = {0.01*c/abs(target)*100:+.3f} +- "
                               f"{0.01*np.sqrt(cov[1, 1])/abs(target)*100:.3f} %, chi2/dof = {chi2:.1f}, "
                               f"points {sel.sum()}")
        for ax in (ax_x, ax_y):
            ax.axvline(1.0 / omega_tau, color=MUTED, linestyle=":", linewidth=1.0)
            ax.set_xscale("log")
            style(ax)
        ax_x.axhline(0.0, color=MUTED, linewidth=0.8)
        ax_y.axhline(0.0, color=MUTED, linewidth=0.8)
        ax_x.text(1.0 / omega_tau, ax_x.get_ylim()[1], " Ω_a Δt = 1", color=MUTED, fontsize=8, va="top")
        ax_x.set_title(f"W{z}+ (Ω_a τ_s = {omega_tau:.0f})", fontsize=10, color=INK, loc="left")
        ax_y.set_xlabel("collision step Δt / τ_s")
    axes[0, 0].set_ylabel("d⟨X⟩/dt relative to first-order target [%]")
    axes[1, 0].set_ylabel("d⟨Y⟩/dt relative to |target| [%]")
    axes[0, 0].legend(frameon=False, fontsize=8, labelcolor=INK)
    fig.suptitle("B08: screening drift versus collision step, ∇T = 100 eV/m, paired difference, "
                 "slope over Ω_a t = 300–700", fontsize=11, color=INK)
    fig.tight_layout()
    fig.savefig(work / "b08_dt.png", dpi=140)


def plot_relax(work):
    rows = read(work / "b08_relax.csv")
    groups = defaultdict(list)
    for r in rows:
        groups[(r["t_over_tau_s"], r["dt_over_tau_s"])].append(r)
    fig, axes = plt.subplots(1, 3, figsize=(13, 4), facecolor=SURFACE)
    colors = {1.0: BLUE, 3.0: ORANGE, 10.0: AQUA}
    for key, ylabel, ax in (("p_par_over_T", "p_∥ / (n T_b)", axes[0]), ("p_perp_over_T", "p_⊥ / (n T_b)", axes[1]),
                            ("u_par_cm_s", "u_∥ [m/s]", axes[2])):
        for t, color in colors.items():
            dts = sorted({k[1] for k in groups if k[0] == t})
            m, s = zip(*[seed_stats([r[key] for r in groups[(t, dt)]]) for dt in dts])
            m, s = np.array(m), np.array(s)
            if key == "u_par_cm_s":
                m, s = m / 100.0, s / 100.0
            ax.errorbar(dts, m, yerr=s, color=color, marker="o", markersize=4, linewidth=1.4, capsize=2,
                        label=f"t = {t:g} τ_s")
        ax.set_xscale("log")
        ax.set_xlabel("collision step Δt / τ_s")
        ax.set_ylabel(ylabel)
        style(ax)
    axes[0].axhline(1.0, color=MUTED, linewidth=0.8)
    axes[1].axhline(1.0, color=MUTED, linewidth=0.8)
    axes[0].legend(frameon=False, fontsize=8, labelcolor=INK)
    fig.suptitle("B08: B04 relaxation of the drifting bi-Maxwellian W3+ versus collision step "
                 "(velocity driver)", fontsize=11, color=INK)
    fig.tight_layout()
    fig.savefig(work / "b08_relax.png", dpi=140)


def plot_n(work, summary):
    rows = read(work / "b08_n.csv")
    ns = sorted({r["n"] for r in rows})
    fig, ax = plt.subplots(figsize=(6.5, 4.5), facecolor=SURFACE)
    for run, color, label in (("gradient", BLUE, "gradient run alone"), ("difference", ORANGE, "paired difference")):
        reported, scatter = [], []
        for n in ns:
            g = [r for r in rows if r["n"] == n]
            reported.append(np.mean([r[f"se_{run}_m_s"] for r in g]))
            scatter.append(np.std([r[f"dXdt_{run}_m_s"] for r in g], ddof=1))
            summary.append(f"n={int(n):8d} {run:10s} reported se = {reported[-1]:.5f} m/s, seed scatter = "
                           f"{scatter[-1]:.5f} m/s ({len(g)} seeds), mean = "
                           f"{np.mean([r[f'dXdt_{run}_m_s'] for r in g]):.5f} m/s")
        ax.plot(ns, reported, "-", color=color, linewidth=1.6, label=f"{label}: reported")
        ax.plot(ns, scatter, "o", color=color, markersize=6, label=f"{label}: seed scatter")
        ref = reported[0] * np.sqrt(ns[0] / np.array(ns))
        ax.plot(ns, ref, ":", color=MUTED, linewidth=1.0)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("ion number N")
    ax.set_ylabel("standard error of ⟨X⟩/t [m/s]")
    ax.set_title("B08: statistics, W3+, 100 eV/m, Δt = 0.02 τ_s (dotted: N^-1/2)", fontsize=10, color=INK,
                 loc="left")
    ax.legend(frameon=False, fontsize=8, labelcolor=INK)
    style(ax)
    fig.tight_layout()
    fig.savefig(work / "b08_n.png", dpi=140)


def main():
    work = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE / "run"
    summary = []
    if (work / "b08_dt.csv").exists():
        plot_dt(work, summary)
    if (work / "b08_relax.csv").exists():
        plot_relax(work)
    if (work / "b08_n.csv").exists():
        plot_n(work, summary)
    (work / "b08_summary.txt").write_text("\n".join(summary) + "\n")
    print("\n".join(summary))


if __name__ == "__main__":
    main()
