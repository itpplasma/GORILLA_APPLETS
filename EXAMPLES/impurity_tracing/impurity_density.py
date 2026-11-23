#!/usr/bin/env python3
"""Absolute impurity density from an impurity tracing run (GORILLA_APPLETS i_option = 17).

usage: impurity_density.py RUNDIR [--prefix impurity_tracing] [--source 1e19] [--phi 0]
                           [--cut R1 Z1 R2 Z2] [--npts 400] [--out FILE.png]

Steady state for a constant source of S particles/s (--source), each of the N markers standing for
S/N particles/s:  n_k(cell) = (S/N) * T_k(cell) / V(cell),  T_k = residence time of charge state k
summed over the markers (<prefix>_residence.dat), V from <prefix>_cells.dat. Only <prefix>_* files
are read, no mesh files.

residence_output = '2d': cells are 2D triangles, V = ring volume of one field period (toroidal
average). residence_output = '3d': cells are prisms of one toroidal slice; the map is the plane
phi = --phi [deg]: inside a slice that slice, on a slice boundary (e.g. phi = 0) the two adjacent
slices combined (time and volume summed).

Plots: total density (log) and one map per charge state; with --cut, density per charge state
along the straight line (R1, Z1) -> (R2, Z2) [cm] vs distance from (R1, Z1) (also written as text,
<out>_cut.txt). Markers that end at t_max_s are counted only up to t_max_s: their fraction is
printed and should be small for a steady state.
"""
import argparse
import os

import numpy as np


def read_residence(path):
    with open(path) as f:
        f.readline(); f.readline()
        head = f.readline().split()
        ntri, n_part, t_max = int(head[0]), int(head[1]), float(head[2])
        ncol = int(head[3]) if len(head) > 3 else 1
        nsl = int(head[4]) if len(head) > 4 else 1
        res = np.array(f.read().split(), dtype=float).reshape(nsl, ntri, ncol)
    return res, n_part, t_max


def read_cells(path):
    with open(path) as f:
        f.readline(); f.readline()
        head = f.readline().split()
        ntri, nphi, nfp, dphi = int(head[0]), int(head[1]), int(head[2]), float(head[3])
    c = np.loadtxt(path, skiprows=3, ndmin=2)
    assert len(c) == ntri
    return dict(region=c[:, 1].astype(int) == 1, tri=c[:, 2:8].reshape(-1, 3, 2), area=c[:, 8], rc=c[:, 9],
                v_prism=c[:, 10], v_ring=c[:, 11], nphi=nphi, nfp=nfp, dphi=dphi)


def plane_density(res, cells, n_part, source, phi_deg):
    """density [m^-3] per triangle and charge state in the plane phi (3d) or toroidal average (2d)"""
    nsl = res.shape[0]
    if nsl == 1:
        t, v, label = res[0], cells["v_ring"], "toroidal average"
    else:
        dphi = cells["dphi"]
        period = nsl * dphi
        x = (np.deg2rad(phi_deg) % period) / dphi
        j = int(np.floor(x + 1e-9)) % nsl
        if abs(x - round(x)) < 1e-6:          # on a slice boundary: the two adjacent slices
            jm = (j - 1) % nsl
            t, v = res[j] + res[jm], 2 * cells["v_prism"]
            label = f"phi = {phi_deg:g} deg (slices {jm + 1}, {j + 1})"
        else:
            t, v = res[j], cells["v_prism"]
            label = f"phi = {phi_deg:g} deg (slice {j + 1})"
    n = source / n_part * t / v[:, None] * 1e6        # cm^-3 -> m^-3
    return n, label


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("rundir")
    ap.add_argument("--prefix", default="impurity_tracing")
    ap.add_argument("--source", type=float, default=1.0, help="source rate [particles/s] (default 1)")
    ap.add_argument("--phi", type=float, default=0.0, help="plane phi [deg] for 3d output (default 0)")
    ap.add_argument("--cut", type=float, nargs=4, metavar=("R1", "Z1", "R2", "Z2"), default=None)
    ap.add_argument("--npts", type=int, default=400)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    pre = os.path.join(a.rundir, a.prefix)
    out = a.out or pre + "_density.png"

    res, n_part, t_max = read_residence(pre + "_residence.dat")
    cells = read_cells(pre + "_cells.dat")
    ex = np.loadtxt(pre + "_exits.dat", comments="#", ndmin=2)
    n_tmax = int(np.sum(ex[:, 1] == 0))
    print(f"{n_part} markers; {n_tmax} ({100 * n_tmax / n_part:.2f}%) reached t_max = {t_max:g} s"
          + ("  <- not a steady state" if n_tmax > 0.01 * n_part else ""))
    n, label = plane_density(res, cells, n_part, a.source, a.phi)
    ncol = n.shape[1]
    states = [k for k in range(ncol) if res[..., k].sum() > 0]
    names = {k: (f"charge {k}" if ncol > 1 else "density") for k in range(ncol)}
    ntot = n.sum(axis=1)
    print(f"{label}: max total density {ntot.max():.3e} m^-3 (source {a.source:g} /s)")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    from matplotlib.tri import Triangulation

    # shared vertices for tripcolor / TriFinder
    pts = cells["tri"].reshape(-1, 2)
    uniq, inv = np.unique(np.round(pts, 9), axis=0, return_inverse=True)
    tris = inv.reshape(-1, 3)
    reg = cells["region"]
    tr = Triangulation(uniq[:, 0], uniq[:, 1], tris)

    panels = [("total", ntot)] + ([(names[k], n[:, k]) for k in states] if ncol > 1 else [])
    fig, axs = plt.subplots(1, len(panels), figsize=(3.8 * len(panels), 6.5), squeeze=False)
    for ax, (title, val) in zip(axs[0], panels):
        pos = val[reg & (val > 0)]
        if len(pos) == 0:
            ax.set_title(title + " (empty)"); continue
        norm = LogNorm(vmin=max(pos.min(), pos.max() * 1e-4), vmax=pos.max())
        m = reg & (val > 0)
        tc = ax.tripcolor(uniq[:, 0], uniq[:, 1], tris[m], facecolors=val[m], norm=norm, cmap="viridis",
                          edgecolors="none")
        ax.triplot(uniq[:, 0], uniq[:, 1], tris[reg], lw=0.05, color="0.8")
        fig.colorbar(tc, ax=ax, shrink=0.6, label="n [m$^{-3}$]")
        if a.cut:
            ax.plot([a.cut[0], a.cut[2]], [a.cut[1], a.cut[3]], "r-", lw=1)
        ax.set_aspect("equal"); ax.set_xlabel("R [cm]"); ax.set_title(title, fontsize=10)
    axs[0][0].set_ylabel("Z [cm]")
    fig.suptitle(f"{os.path.basename(os.path.abspath(a.rundir))}: {label}, source {a.source:g} /s, "
                 f"{n_part} markers", fontsize=10)
    fig.tight_layout()
    fig.savefig(out, dpi=110)
    print("plot:", out)

    if a.cut:
        R1, Z1, R2, Z2 = a.cut
        s = np.linspace(0.0, 1.0, a.npts)
        Rp, Zp = R1 + s * (R2 - R1), Z1 + s * (Z2 - Z1)
        it = tr.get_trifinder()(Rp, Zp)
        ok = it >= 0
        dist = s * np.hypot(R2 - R1, Z2 - Z1)
        prof = np.full((a.npts, ncol), np.nan)
        prof[ok] = n[it[ok]]
        txt = out.replace(".png", "_cut.txt")
        np.savetxt(txt, np.c_[dist, Rp, Zp, it, prof],
                   header=f"{label}; distance[cm] R[cm] Z[cm] triangle(0-based, -1 outside) "
                          f"n[m^-3] per charge state 0..{ncol - 1}")
        fig2, ax2 = plt.subplots(figsize=(7, 4.5))
        for k in states:
            ax2.plot(dist, prof[:, k], label=names[k])
        ax2.set_xlabel(f"distance from ({R1:g}, {Z1:g}) cm [cm]")
        ax2.set_ylabel("n [m$^{-3}$]")
        ax2.set_title(f"cut ({R1:g},{Z1:g}) -> ({R2:g},{Z2:g}) cm, {label}", fontsize=9)
        ax2.legend()
        fig2.tight_layout()
        out2 = out.replace(".png", "_cut.png")
        fig2.savefig(out2, dpi=110)
        print("plot:", out2, " data:", txt)


if __name__ == "__main__":
    main()
