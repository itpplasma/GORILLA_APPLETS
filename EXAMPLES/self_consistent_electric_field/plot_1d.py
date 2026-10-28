#!/usr/bin/env python3
"""1D plots of a self-consistent electric field run (Python port of the
"%% 1d plot" section of self_consistent_ef_plots.m).

Figures:
  1. particle source across source-update iterations (particle_source_<k>.dat)
  2. 1D electric potential of one source step (phi_elec_1d_<source_step>_<i>.dat)
  3. 1D ion/electron densities of one source step (one_d_densities<source_step>_<i>.dat)
     together with the desired density profile

By default the last source step and the last electric-potential update found in
the run folder are plotted, so the script works for both source-tuning runs and
fixed-source runs (n_source_updates = 0). Use the options to pick other ones.

Usage:
  python3 plot_1d.py [run_folder] [--source-step K] [--n-low I] [--n-high I]
                     [--n-source-low K] [--n-source-high K] [--save]
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def available_indices(folder: Path, pattern: str) -> list[tuple[int, ...]]:
    """Integer groups of all files in folder whose names fully match pattern."""
    regex = re.compile(pattern)
    hits = []
    for path in folder.iterdir():
        m = regex.fullmatch(path.name)
        if m:
            hits.append(tuple(int(g) for g in m.groups()))
    return sorted(hits)


def desired_density(s: np.ndarray) -> np.ndarray:
    """Target density profile in cm^-3 (same formula as calc_desired_density_profile)."""
    dens_e20 = 0.37 * (0.74 + 0.26 * (1.0 - s**2.5)**1.5 - 0.06 * (1.0 - np.exp(-4.0 * s**2)))
    return dens_e20 * 1.0e14


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("folder", nargs="?", type=Path, default=Path.cwd(),
                        help="run folder (default: current directory)")
    parser.add_argument("--source-step", type=int,
                        help="source-update iteration to plot (default: last available)")
    parser.add_argument("--n-low", type=int,
                        help="first electric-potential update to overlay (default: --n-high)")
    parser.add_argument("--n-high", type=int,
                        help="last electric-potential update to overlay (default: last available)")
    parser.add_argument("--n-source-low", type=int,
                        help="first particle_source_<k>.dat to overlay (default: last 11)")
    parser.add_argument("--n-source-high", type=int,
                        help="last particle_source_<k>.dat to overlay (default: last available)")
    parser.add_argument("--save", action="store_true",
                        help="write PNGs into the run folder instead of showing the figures")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    folder = args.folder.resolve()
    if not folder.is_dir():
        sys.exit(f"run folder not found: {folder}")

    phi_indices = available_indices(folder, r"phi_elec_1d_(\d+)_(\d+)\.dat")
    source_indices = [k for (k,) in available_indices(folder, r"particle_source_(\d+)\.dat")]
    if not phi_indices or not source_indices:
        sys.exit(f"no phi_elec_1d_*.dat / particle_source_*.dat files in {folder}")

    source_step = args.source_step or max(k for k, _ in phi_indices)
    n_high = args.n_high or max(i for k, i in phi_indices if k == source_step)
    n_low = args.n_low or n_high
    n_source_high = args.n_source_high or max(source_indices)
    n_source_low = args.n_source_low or max(min(source_indices), n_source_high - 10)

    # Particle source across source-update iterations
    fig_source, ax = plt.subplots()
    colors = plt.cm.jet(np.linspace(0.0, 1.0, n_source_high - n_source_low + 1))
    for color, k in zip(colors, range(n_source_low, n_source_high + 1)):
        data = np.loadtxt(folder / f"particle_source_{k}.dat")
        ax.plot(data[:, 0], data[:, 1], color=color, label=f"source_step = {k}")
    ax.grid(True)
    ax.set_xlabel("s")
    ax.set_ylabel("particle source")
    ax.set_title("particle source across source-update iterations")
    ax.legend(loc="best")

    s_vertices = data[:, 0]
    s_layers = 0.5 * (s_vertices[:-1] + s_vertices[1:])

    # Electric potential: first update in red, updates n_low..n_high in blue
    fig_potential, ax = plt.subplots()
    start_potential = np.loadtxt(folder / f"phi_elec_1d_{source_step}_1.dat")
    ax.plot(s_vertices, start_potential, "r", label="update 1")
    for i in range(n_low, n_high + 1):
        potential = np.loadtxt(folder / f"phi_elec_1d_{source_step}_{i}.dat")
        ax.plot(s_vertices, potential, "b", label=f"update {n_low}-{n_high}" if i == n_low else None)
    ax.grid(True)
    ax.set_xlabel("s")
    ax.set_ylabel("potential")
    ax.set_title(f"electric potential, source_step = {source_step}")
    ax.legend(loc="best")

    # Densities: first update (ions red, electrons blue), updates n_low..n_high
    # (ions magenta, electrons black), desired profile green
    fig_density, ax = plt.subplots()
    one_d = np.loadtxt(folder / f"one_d_densities{source_step}_1.dat")
    ax.plot(s_layers, one_d[:, 1], "b", label="electrons, update 1")
    ax.plot(s_layers, one_d[:, 0], "r", label="ions, update 1")
    ax.plot(s_layers, desired_density(s_layers), "g", label="desired density")
    for i in range(n_low, n_high + 1):
        one_d = np.loadtxt(folder / f"one_d_densities{source_step}_{i}.dat")
        first = i == n_low
        ax.plot(s_layers, one_d[:, 0], "m", label=f"ions, update {n_low}-{n_high}" if first else None)
        ax.plot(s_layers, one_d[:, 1], "k", label=f"electrons, update {n_low}-{n_high}" if first else None)
    ax.grid(True)
    ax.set_xlabel("s")
    ax.set_ylabel(r"density [cm$^{-3}$]")
    ax.set_title(f"1D densities, source_step = {source_step}")
    ax.legend(loc="best")

    if args.save:
        for fig, name in [(fig_source, "particle_source"), (fig_potential, "phi_elec_1d"),
                          (fig_density, "one_d_densities")]:
            out = folder / f"{name}_plot.png"
            fig.savefig(out, dpi=150, bbox_inches="tight")
            print(f"wrote {out}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
