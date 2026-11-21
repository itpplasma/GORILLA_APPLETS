#!/usr/bin/env python3
"""Plots for the impurity tracing (GORILLA_APPLETS i_option = 17).

usage: region_plots.py RUNDIR [--mesh MESHDIR] [--prefix impurity_tracing] [--bin 2.0] [--out FILE.png]

RUNDIR holds <prefix>_exits.dat and <prefix>_residence.dat; MESHDIR (default RUNDIR/mesh) the
extended mesh with s_edge_region. Four panels:
  (a) density per unit source rate, n/S = (time per triangle)/(N * V_triangle) [cm^-3 per marker/s],
      V_triangle = 2 pi R_centroid * area (axisymmetric ring);
  (b) the region boundaries as 1D contours coloured by hits per unit length (per marker);
  (c), (d) hits per unit length vs arc length along the inner (s = s_edge) and outer (wall) contour,
      arc length counter-clockwise in (R, Z) from the outer midplane (Z of the magnetic axis).
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from region_check import boundary_edges, load, loops, nml_s_edge, region_flags  # noqa: E402


def ordered_loop(edges, P):
    lps, _ = loops(edges)
    if len(lps) != 1:
        raise SystemExit(f"expected one closed boundary loop, found {len(lps)}")
    v = np.array(lps[0])
    xy = P[v]
    area = 0.5 * np.sum(xy[:, 0] * np.roll(xy[:, 1], -1) - np.roll(xy[:, 0], -1) * xy[:, 1])
    if area < 0:
        v = v[::-1]
    return v


def start_at_omp(v, P, z_axis):
    """Rotate the loop to start at the segment that crosses Z = z_axis with the largest R; return the
    closed polyline (first point on Z = z_axis) and its cumulative arc length."""
    xy = P[v]
    nxt = np.roll(xy, -1, axis=0)
    best, k_best, p_best = -np.inf, None, None
    for k in range(len(xy)):
        z0, z1 = xy[k, 1] - z_axis, nxt[k, 1] - z_axis
        if z0 <= 0 < z1:  # upward crossing (counter-clockwise on the low-field side)
            t = -z0 / (z1 - z0)
            p = xy[k] + t * (nxt[k] - xy[k])
            if p[0] > best:
                best, k_best, p_best = p[0], k, p
    poly = np.vstack([p_best, np.roll(xy, -(k_best + 1), axis=0), p_best])
    s = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(poly, axis=0), axis=1))])
    return poly, s


def project(points, poly, s):
    """Arc length of the nearest point on the closed polyline for each point."""
    a, b = poly[:-1], poly[1:]
    ab = b - a
    L2 = np.maximum(np.sum(ab * ab, axis=1), 1e-30)
    out = np.empty(len(points))
    dist = np.empty(len(points))
    for i, p in enumerate(points):
        t = np.clip(np.sum((p - a) * ab, axis=1) / L2, 0, 1)
        q = a + t[:, None] * ab
        d = np.linalg.norm(q - p, axis=1)
        k = np.argmin(d)
        out[i] = s[k] + t[k] * np.sqrt(L2[k])
        dist[i] = d[k]
    return out, dist


def landmarks(poly, s):
    """Arc length of top, high-field-side midplane and bottom of a contour."""
    zc = poly[0, 1]
    ktop, kbot = np.argmax(poly[:, 1]), np.argmin(poly[:, 1])
    # high-field side: the crossing of Z = Z_axis with smallest R
    k_hfs, rmin = None, np.inf
    for k in range(len(poly) - 1):
        if (poly[k, 1] - zc) * (poly[k + 1, 1] - zc) <= 0 and min(poly[k, 0], poly[k + 1, 0]) < rmin:
            rmin, k_hfs = min(poly[k, 0], poly[k + 1, 0]), k
    return {"top": s[ktop], "HFS": s[k_hfs], "bottom": s[kbot]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("rundir")
    ap.add_argument("--mesh", default=None)
    ap.add_argument("--prefix", default="impurity_tracing")
    ap.add_argument("--bin", type=float, default=2.0, help="arc-length bin [cm] for (c), (d)")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    meshdir = a.mesh or os.path.join(a.rundir, "mesh")
    out = a.out or os.path.join(a.rundir, f"{a.prefix}_plots.png")

    P, T, s_v, reg = load(meshdir)
    s_edge = nml_s_edge(meshdir)
    flag = region_flags(T, s_v, reg, s_edge)
    allb = set(boundary_edges(T, np.ones(len(T), bool)))
    rb = boundary_edges(T, flag)
    outer_e = [e for e in rb if e in allb]
    inner_e = [e for e in rb if e not in allb]
    core_v = (s_v > 0)
    z_axis = P[core_v][np.argmin(s_v[core_v])].mean() if False else P[s_v == s_v[core_v].min()][:, 1].mean()

    ex = np.loadtxt(os.path.join(a.rundir, f"{a.prefix}_exits.dat"))
    status = ex[:, 1].astype(int)
    n_mark = len(ex)
    with open(os.path.join(a.rundir, f"{a.prefix}_residence.dat")) as f:
        f.readline(); f.readline()
        head = f.readline().split()
        ntri, n_part, t_max = int(head[0]), int(head[1]), float(head[2])
        ncol = int(head[3]) if len(head) > 3 else 1
        res_states = np.array(f.read().split(), dtype=float).reshape(ntri, ncol)
        res = res_states.sum(axis=1)
    assert ntri == len(T) and n_part == n_mark

    # (a) density per unit source rate
    tri = P[T]
    area = 0.5 * np.abs((tri[:, 1, 0] - tri[:, 0, 0]) * (tri[:, 2, 1] - tri[:, 0, 1])
                        - (tri[:, 2, 0] - tri[:, 0, 0]) * (tri[:, 1, 1] - tri[:, 0, 1]))
    vol = 2 * np.pi * tri[:, :, 0].mean(axis=1) * area
    dens = res / (n_mark * vol)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    from matplotlib.colors import LogNorm

    fig = plt.figure(figsize=(14, 10))
    gs = fig.add_gridspec(2, 3, width_ratios=[1, 1, 1.4])
    axa = fig.add_subplot(gs[:, 0])
    axb = fig.add_subplot(gs[:, 1])
    axc = fig.add_subplot(gs[0, 2])
    axd = fig.add_subplot(gs[1, 2])

    pos = dens[dens > 0]
    norm = LogNorm(vmin=max(pos.min(), pos.max() * 1e-5), vmax=pos.max()) if len(pos) else None
    m = np.where(flag & (dens > 0))[0]
    tc = axa.tripcolor(P[:, 0], P[:, 1], T[m], facecolors=dens[m], norm=norm, cmap="viridis", edgecolors="none")
    fig.colorbar(tc, ax=axa, shrink=0.6, label="n / S  [cm$^{-3}$ per marker/s]")
    start = ex[:, 8:10]
    axa.plot(start[:, 0], start[:, 1], "r+", ms=8, label="start")
    for ax in (axa, axb):
        ax.add_collection(LineCollection(P[np.array(outer_e)], colors="0.3", lw=0.6))
        ax.add_collection(LineCollection(P[np.array(inner_e)], colors="0.3", lw=0.6))
        ax.set_aspect("equal")
        ax.set_xlabel("R [cm]")
        ax.autoscale_view()
    axa.set_ylabel("Z [cm]")
    axa.set_title(f"(a) density per unit source, s_edge = {s_edge:g}")
    axa.legend(loc="upper right", fontsize=8)

    # boundary hits
    hit_sets = []
    for e_list, st, name, ax1d in ((inner_e, 2, "inner boundary (s = s_edge)", axc),
                                    (outer_e, 1, "outer boundary (wall)", axd)):
        v = ordered_loop(e_list, P)
        poly, s = start_at_omp(v, P, z_axis)
        pts = ex[status == st][:, [3, 5]]
        sh, dist = project(pts, poly, s) if len(pts) else (np.array([]), np.array([]))
        if len(dist) and dist.max() > 0.5:
            print(f"warning: {name}: exit point up to {dist.max():.2f} cm from the contour")
        # per segment hits per unit length (per marker), for the coloured contour
        seg_len = np.diff(s)
        kseg = np.clip(np.searchsorted(s, sh, side="right") - 1, 0, len(seg_len) - 1)
        per_seg = np.bincount(kseg, minlength=len(seg_len)) / np.maximum(seg_len, 1e-12) / n_mark
        hit_sets.append((poly, per_seg))
        # 1D histogram
        nb = max(1, int(round(s[-1] / a.bin)))
        edges = np.linspace(0, s[-1], nb + 1)
        h, _ = np.histogram(sh, bins=edges)
        ax1d.stairs(h / np.diff(edges) / n_mark, edges, fill=True, alpha=0.5, color="C3" if st == 2 else "C0")
        ax1d.stairs(h / np.diff(edges) / n_mark, edges, color="k", lw=0.6)
        for lab, sl in {"OMP": 0.0, **landmarks(poly, s)}.items():
            ax1d.axvline(sl, color="0.5", ls="--", lw=0.7)
            ax1d.text(sl, 1.0, lab, transform=ax1d.get_xaxis_transform(), ha="center", va="bottom", fontsize=8)
        ax1d.set_xlim(0, s[-1])
        ax1d.set_xlabel("arc length from outer midplane, counter-clockwise [cm]")
        ax1d.set_ylabel("hits per cm per marker")
        ax1d.set_title(f"({'c' if st == 2 else 'd'}) {name}: {len(pts)}/{n_mark} markers "
                       f"({100 * len(pts) / n_mark:.1f}%), length {s[-1]:.0f} cm", fontsize=9, pad=14)

    vmax = max(ps.max() for _, ps in hit_sets)
    vpos = [ps[ps > 0].min() for _, ps in hit_sets if np.any(ps > 0)]
    cnorm = LogNorm(vmin=min(vpos), vmax=vmax) if vpos else None
    for poly, ps in hit_sets:
        segs = np.stack([poly[:-1], poly[1:]], axis=1)
        lc = LineCollection(segs, array=np.where(ps > 0, ps, np.nan), cmap="plasma", norm=cnorm, lw=4)
        lc.cmap.set_bad("0.85")
        axb.add_collection(lc)
    axb.triplot(P[:, 0], P[:, 1], T[flag], lw=0.08, color="0.75")
    fig.colorbar(lc, ax=axb, shrink=0.6, label="hits per cm per marker (grey: none)")
    axb.set_title("(b) boundary hits")

    nst = {k: int(np.sum(status == k)) for k in (1, 2, 0, 3, -1)}
    fig.suptitle(f"{n_mark} markers, t_max = {t_max:g} s: wall {nst[1]}, inner boundary {nst[2]}, "
                 f"time limit {nst[0]}, pusher removed {nst[3]}, started outside {nst[-1]}", fontsize=11)
    fig.tight_layout()
    fig.savefig(out, dpi=120)
    print("plot:", out)

    # density per charge state (charge-state mode: one residence column per state 0 ... Z_nuc)
    if ncol > 1:
        states = [z for z in range(ncol) if res_states[:, z].sum() > 0]
        fig2, axs = plt.subplots(1, len(states), figsize=(3.6 * len(states), 6.5), squeeze=False)
        dens_z = res_states / (n_mark * vol[:, None])
        pos = dens_z[dens_z > 0]
        nz = LogNorm(vmin=max(pos.min(), pos.max() * 1e-5), vmax=pos.max())
        for ax, z in zip(axs[0], states):
            mz = np.where(flag & (dens_z[:, z] > 0))[0]
            tcz = ax.tripcolor(P[:, 0], P[:, 1], T[mz], facecolors=dens_z[mz, z], norm=nz, cmap="viridis",
                               edgecolors="none")
            ax.add_collection(LineCollection(P[np.array(outer_e)], colors="0.3", lw=0.6))
            ax.add_collection(LineCollection(P[np.array(inner_e)], colors="0.3", lw=0.6))
            ax.set_aspect("equal"); ax.autoscale_view(); ax.set_xlabel("R [cm]")
            ax.set_title(f"charge {z}: {100 * res_states[:, z].sum() / res_states.sum():.1f}% of time", fontsize=9)
        axs[0][0].set_ylabel("Z [cm]")
        fig2.colorbar(tcz, ax=axs[0].tolist(), shrink=0.7, label="n / S  [cm$^{-3}$ per marker/s]")
        out2 = out.replace(".png", "_charge_states.png")
        fig2.savefig(out2, dpi=110)
        print("plot:", out2)
    t_exit = ex[status > 0, 2]
    if len(t_exit):
        print(f"exit times: median {np.median(t_exit):.3e} s, p90 {np.percentile(t_exit, 90):.3e} s")


if __name__ == "__main__":
    main()
