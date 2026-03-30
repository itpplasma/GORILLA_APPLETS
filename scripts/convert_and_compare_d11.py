#!/usr/bin/env python3
"""Compare GORILLA KM D11 against NEO-2-QL and NEOART for AUG 39461."""

import argparse
import sys

import h5py
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
from scipy.io import netcdf_file


def load_gorilla_csv(path):
    data = np.genfromtxt(path, delimiter=",", names=True)
    return data


def load_neo2(path):
    with h5py.File(path, "r") as f:
        boozer_s = f["boozer_s"][:]
        num_spec = int(f["num_spec"][0])
        d11_ax = f["D11_AX"][:].reshape(boozer_s.size, num_spec, num_spec)
        av_nabla_stor = f["av_nabla_stor"][:]
        z_spec = f["z_spec"][0, :]
    return boozer_s, d11_ax, av_nabla_stor, z_spec


def load_neoart(path, charge_stages=(10,)):
    neoart = netcdf_file(path, "r")
    rho_pol = neoart.variables["rho_poloidal_grid"][:].copy()
    d_neo = neoart.variables["D_NEO"][:][0].copy()
    result = []
    for stage in charge_stages:
        result.append(d_neo[0][stage - 1])
    neoart.close()
    return rho_pol, np.array(result)


def convert_gorilla_d11(d11_s, gorilla_s, av_nabla_stor, neo2_s):
    av_nabla_interp = interp1d(neo2_s, av_nabla_stor, kind="cubic", fill_value="extrapolate")
    av_nabla_at_gorilla = av_nabla_interp(gorilla_s)
    d11_r = d11_s / av_nabla_at_gorilla**2
    return d11_r


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gorilla_csv", help="GORILLA km_d11_profile.csv")
    parser.add_argument("neo2_h5", help="NEO-2 neon_discharge_out.h5")
    parser.add_argument("--neoart-nc", help="NEOART neon_discharge_neoart.nc")
    parser.add_argument("--species-index", type=int, default=0,
                        help="NEO-2 species index for D11 diagonal (0=e, 1=D, 2=Ne10)")
    parser.add_argument("--output", default="/tmp/gorilla_d11_comparison.png")
    parser.add_argument("--show", action="store_true")
    args = parser.parse_args()

    gorilla = load_gorilla_csv(args.gorilla_csv)
    neo2_s, d11_ax_neo2, av_nabla_stor, z_spec = load_neo2(args.neo2_h5)
    gorilla_s = gorilla["s"]
    gorilla_d11_r = convert_gorilla_d11(gorilla["d11_s"], gorilla_s, av_nabla_stor, neo2_s)
    gorilla_convolved_r = None
    if "d11_convolved_s" in gorilla.dtype.names:
        gorilla_convolved_r = convert_gorilla_d11(
            gorilla["d11_convolved_s"], gorilla_s, av_nabla_stor, neo2_s
        )

    sp = args.species_index
    species_labels = {0: "e", 1: "D", 2: "Ne$^{10+}$", 3: "Ne$^{9+}$", 4: "Ne$^{8+}$"}
    species_label = species_labels.get(sp, f"species {sp}")

    fontsize = 18
    markersize = 8
    linewidth = 2
    figsize = (14, 10)

    fig, axes = plt.subplots(3, 1, figsize=(14, 13), height_ratios=[3, 1.2, 1])

    # D11 comparison
    ax = axes[0]
    ax.plot(neo2_s, d11_ax_neo2[:, sp, sp] * 1e-4, "o-",
            label=f"NEO-2-QL {species_label}", markersize=markersize, linewidth=linewidth)
    ax.plot(gorilla_s, gorilla_d11_r * 1e-4, "s-",
            label=f"GORILLA KM {species_label}", markersize=markersize, linewidth=linewidth)
    if gorilla_convolved_r is not None:
        ax.plot(gorilla_s, gorilla_convolved_r * 1e-4, "d--",
                label=f"GORILLA convolved {species_label}",
                markersize=markersize - 1, linewidth=linewidth)

    if args.neoart_nc:
        rho_pol_neoart, d11_neoart = load_neoart(args.neoart_nc, charge_stages=[10])
        # NEOART uses rho_pol, but we plot vs boozer_s, so skip NEOART for now
        # unless we have a mapping. Instead, just note it in the legend.
        ax.set_xlabel(r"$s_\mathrm{tor}$ (boozer_s)", fontsize=fontsize)
    else:
        ax.set_xlabel(r"$s_\mathrm{tor}$ (boozer_s)", fontsize=fontsize)

    ax.set_ylabel(r"$D_{11}$ [m$^2$/s]", fontsize=fontsize)
    ax.set_yscale("log")
    ax.legend(fontsize=fontsize - 2)
    ax.grid(True, alpha=0.3)
    ax.tick_params(axis="both", labelsize=fontsize - 2)
    ax.set_title("GORILLA vs NEO-2-QL: mono-energetic $D_{11}$", fontsize=fontsize)

    # Ratio to NEO-2
    ax2 = axes[1]
    neo2_diag_interp = interp1d(
        neo2_s, d11_ax_neo2[:, sp, sp] * 1e-4, kind="linear", fill_value="extrapolate"
    )
    ratio_mono = (gorilla_d11_r * 1e-4) / neo2_diag_interp(gorilla_s)
    ax2.plot(gorilla_s, ratio_mono, "s-", markersize=markersize, linewidth=linewidth,
             label="GORILLA / NEO-2")
    if gorilla_convolved_r is not None:
        ratio_conv = (gorilla_convolved_r * 1e-4) / neo2_diag_interp(gorilla_s)
        ax2.plot(gorilla_s, ratio_conv, "d--", markersize=markersize - 1, linewidth=linewidth,
                 label="GORILLA convolved / NEO-2")
    ax2.axhline(1.0, color="k", linewidth=0.5)
    ax2.set_xlabel(r"$s_\mathrm{tor}$ (boozer_s)", fontsize=fontsize)
    ax2.set_ylabel("ratio", fontsize=fontsize)
    ax2.legend(fontsize=fontsize - 2)
    ax2.grid(True, alpha=0.3)
    ax2.tick_params(axis="both", labelsize=fontsize - 2)

    # Drift coefficient A(s)
    ax3 = axes[2]
    ax3.plot(gorilla_s, gorilla["convection_A"], "s-", markersize=markersize, linewidth=linewidth,
             label="GORILLA A(s)")
    ax3.axhline(0, color="k", linewidth=0.5)
    ax3.set_xlabel(r"$s_\mathrm{tor}$ (boozer_s)", fontsize=fontsize)
    ax3.set_ylabel("A(s) [1/s]", fontsize=fontsize)
    ax3.legend(fontsize=fontsize - 2)
    ax3.grid(True, alpha=0.3)
    ax3.tick_params(axis="both", labelsize=fontsize - 2)

    plt.tight_layout()
    plt.savefig(args.output, dpi=150, bbox_inches="tight")
    print(f"Plot saved to {args.output}")

    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
