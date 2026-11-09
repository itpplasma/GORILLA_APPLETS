#!/usr/bin/env python3
"""
Plot the counterparts of Homma & Hatayama (2014) Figs. 1-4 from the runs of
run_homma2014.py.

  Fig. 1 / 3  XY trajectories of two test ions (full orbit), dt = 0.05 / 0.15 tau_s
  Fig. 2 / 4  mean positions <X>, <Y> versus Omega_a t for the full orbit
              (Homma's model) and the guiding centre (GORILLA model), with
              Eq. (8) (Trubnikov collision time, as in Homma) X = v_scr t, Y = 0

Means are those of the gradient run, as in the paper (no control subtraction).
Usage: python3 plot_homma2014.py [work_dir]   (default: ./run)
Writes <work_dir>/homma2014_figs_1-4.png.
"""
import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
TARGETS = HERE.parent.parent / "TESTS" / "REFERENCE" / "homma" / "B07_targets.csv"

# reference palette (dataviz skill), light mode
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"
MODELS = {"full_orbit": ("full orbit (Homma model)", BLUE),
          "guiding_centre": ("guiding centre (GORILLA model)", ORANGE)}


def read_csv(path):
    with open(path) as f:
        return list(csv.DictReader(line for line in f if not line.startswith("#")))


def eq8_trubnikov():
    for r in read_csv(TARGETS):
        if int(r["Z"]) == 3 and float(r["gradT_x_eV_m"]) == 300.0:
            return float(r["eq8_trubnikov_m_s"]), float(r["gyro_orbit_dXdt_m_s"])
    raise RuntimeError("no W3+ 300 eV/m row in B07_targets.csv")


def style(ax):
    ax.grid(color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=9)
    ax.xaxis.label.set_color(INK)
    ax.yaxis.label.set_color(INK)


def plot_trajectories(ax, case_dir, dt, fig_no):
    rows = read_csv(case_dir / "thermal_force_trajectories.csv")
    for pid, color, ls, label in ((1, BLUE, "-", "ion 1"), (2, AQUA, "--", "ion 2")):
        pts = sorted((float(r["t_s"]), float(r["x"]) * 10, float(r["y"]) * 10)
                     for r in rows if int(r["particle"]) == pid)
        ax.plot([p[1] for p in pts], [p[2] for p in pts], ls, color=color, linewidth=1.0, label=label)
    ax.plot([0], [0], "o", color=INK, markersize=4)
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_xlabel("X [mm]")
    ax.set_ylabel("Y [mm]")
    ax.set_title(f"Fig. {fig_no}: two test-ion trajectories, Δt = {dt} τ_s", fontsize=10, color=INK, loc="left")
    ax.legend(frameon=False, fontsize=8, labelcolor=INK)
    style(ax)


def plot_means(ax, work, dt, fig_no, v_scr):
    t_end = 0.0
    for model, (label, color) in MODELS.items():
        rows = [r for r in read_csv(work / f"{model}_dt{dt}" / "screening.csv") if r["run"].strip() == "gradient"]
        wt = [float(r["omega_t"]) for r in rows]
        t_end = max(t_end, float(rows[-1]["t_s"]))
        ax.plot(wt, [float(r["mean_x"]) * 10 for r in rows], "-", color=color, linewidth=1.6,
                label=f"⟨X⟩ {label}")
        ax.plot(wt, [float(r["mean_y"]) * 10 for r in rows], "--", color=color, linewidth=1.2,
                label=f"⟨Y⟩ {label}")
    wt_end = 700.0
    ax.plot([0, wt_end], [0, v_scr * t_end * 1e3], ":", color=INK, linewidth=1.4,
            label=f"Eq. (8): X = v_scr t, v_scr = {v_scr:.2f} m/s")
    ax.axhline(0.0, color=MUTED, linewidth=0.8)
    ax.set_xlim(0, wt_end)
    ax.set_xlabel("Ω_a t")
    ax.set_ylabel("mean position [mm]")
    ax.set_title(f"Fig. {fig_no}: ⟨X⟩ and ⟨Y⟩ of 10⁴ ions, Δt = {dt} τ_s", fontsize=10, color=INK, loc="left")
    ax.legend(frameon=False, fontsize=7.5, labelcolor=INK, loc="lower left")
    style(ax)


def main():
    work = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE / "run"
    v_scr, _ = eq8_trubnikov()
    fig, axes = plt.subplots(2, 2, figsize=(11, 9), facecolor="#fcfcfb")
    for row, (dt, figs) in enumerate(((0.05, (1, 2)), (0.15, (3, 4)))):
        plot_trajectories(axes[row, 0], work / f"full_orbit_dt{dt}", dt, figs[0])
        plot_means(axes[row, 1], work, dt, figs[1], v_scr)
    for ax in axes.flat:
        ax.set_facecolor("#fcfcfb")
    fig.suptitle("Homma & Hatayama (2014) temperature-screening test: W³⁺ in H⁺, ∇T = 300 eV/m e_X, B = 1 T e_Z",
                 fontsize=11, color=INK)
    fig.tight_layout()
    out = work / "homma2014_figs_1-4.png"
    fig.savefig(out, dpi=150)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
