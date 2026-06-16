#!/usr/bin/env python3
"""Plot the mono-energetic transport outputs produced by run.sh.

Always emits the diffusion-coefficient figure when nustar_diffcoef_std.dat is
present. When idiffcoef_output >= 2 the raw (Delta s)^2 trends are also plotted
(one figure per nu*).
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

try:
    import f90nml
except ImportError:
    sys.exit("f90nml is required (pip install f90nml)")


HERE = Path(__file__).resolve().parent
RUN = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else HERE / "run_output"
if not RUN.is_dir():
    sys.exit(f"run output directory not found: {RUN}")

PLOTS = RUN / "plots"
PLOTS.mkdir(exist_ok=True)

nml = f90nml.read(str(RUN / "mono_energetic_transp_coef.inp"))["transpcoefnml"]
v_E = float(nml.get("v_e", 0.0))
nu_star_start = float(nml.get("nu_star_start", 1.0))
nu_exp_basis = float(nml.get("nu_exp_basis", 0.5))


def plot_diffusion_coefficient() -> None:
    data = np.loadtxt(RUN / "nustar_diffcoef_std.dat")
    nu_star, d11, sigma = data[:, 0], data[:, 1], data[:, 2]
    fig, ax = plt.subplots()
    ax.errorbar(nu_star, d11, yerr=1.96 * sigma, marker="o", linestyle="-",
                label=f"v_E = {v_E:g}")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$\nu^*$")
    ax.set_ylabel(r"$D_{11}$ [1/s]")
    ax.set_title("Mono-energetic radial diffusion coefficient")
    ax.legend()
    ax.grid(True, which="both", alpha=0.3)
    out = PLOTS / "diffusion_coefficients.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


def split_blocks(path: Path) -> list[np.ndarray]:
    """Split psi2.dat / std_psi2.dat into per-nu* blocks.

    The Fortran writer separates blocks with a blank line. np.loadtxt would
    silently merge them, so we parse manually.
    """
    blocks: list[list[list[float]]] = [[]]
    for line in path.read_text().splitlines():
        if not line.strip():
            if blocks[-1]:
                blocks.append([])
            continue
        blocks[-1].append([float(x) for x in line.split()])
    return [np.asarray(b) for b in blocks if b]


def plot_delta_s_squared_trends() -> None:
    psi2_path = RUN / "psi2.dat"
    std_path = RUN / "std_psi2.dat"
    if not psi2_path.exists() or not std_path.exists():
        return
    psi2_blocks = split_blocks(psi2_path)
    std_blocks = split_blocks(std_path)
    if len(psi2_blocks) != len(std_blocks):
        sys.exit(f"block count mismatch: {len(psi2_blocks)} vs {len(std_blocks)}")
    for i, (psi2, std) in enumerate(zip(psi2_blocks, std_blocks)):
        nu_star = nu_star_start * nu_exp_basis ** i
        fig, ax = plt.subplots()
        ax.errorbar(psi2[:, 0], psi2[:, 1], yerr=1.96 * std[:, 1],
                    marker="x", linestyle="none",
                    label=f"v_E = {v_E:g}\n$\\nu^*$ = {nu_star:g}")
        ax.set_xlabel("t [s]")
        ax.set_ylabel(r"$\langle (\Delta s)^2 \rangle$")
        ax.set_title(r"Trend of $\langle (\Delta s)^2 \rangle$")
        ax.legend()
        ax.grid(True, alpha=0.3)
        out = PLOTS / f"psi2_{i + 1}.png"
        fig.savefig(out, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"wrote {out}")


if (RUN / "nustar_diffcoef_std.dat").exists():
    plot_diffusion_coefficient()
plot_delta_s_squared_trends()
