#!/usr/bin/env python3
"""Check and plot the edge region of an extended AUG mesh.

usage: region_check.py MESHDIR [S_EDGE] [OUT.png]
S_EDGE defaults to s_edge_region from MESHDIR/edge_mesh.nml.
Region rule (same as GORILLA tetra_grid_mod, tetra_in_edge_region):
  s_edge = 0      -> every triangle
  core triangle   -> mean vertex s > s_edge
  edge (region 1) -> yes, axis fill (region 2) -> no
Checks: the region's inner boundary consists only of mesh edges whose two vertices lie on s = s_edge,
and it forms one closed loop; the outer boundary is the mesh boundary (wall).
"""
import re
import sys
from collections import defaultdict

import numpy as np


def load(meshdir):
    P = np.loadtxt(f"{meshdir}/knots_aug_edge.dat", skiprows=1)
    T = np.loadtxt(f"{meshdir}/triangles_aug_edge.dat", skiprows=1, dtype=int) - 1
    with open(f"{meshdir}/edge_mesh_core.dat") as f:
        nv, nt = map(int, f.readline().split())
        rest = np.array(f.read().split(), dtype=float)
    st = rest[: 2 * nv].reshape(nv, 2)
    reg = rest[2 * nv:].astype(int)
    assert len(reg) == nt == len(T) and nv == len(P)
    return P, T, st[:, 0], reg


def nml_s_edge(meshdir):
    txt = open(f"{meshdir}/edge_mesh.nml").read()
    m = re.search(r"s_edge_region\s*=\s*([-+0-9.eEdD]+)", txt)
    return float(m.group(1).replace("d", "e").replace("D", "e")) if m else 0.0


def region_flags(T, s, reg, s_edge):
    if s_edge <= 0:
        return np.ones(len(T), bool)
    core = reg == 0
    flag = reg == 1
    flag[core] = s[T[core]].sum(axis=1) > 3 * s_edge  # triangle: 3 vertices (tetra rule uses 4, same ring set)
    return flag


def boundary_edges(T, mask):
    cnt = defaultdict(int)
    for tri in T[mask]:
        for a, b in ((0, 1), (1, 2), (2, 0)):
            cnt[tuple(sorted((tri[a], tri[b])))] += 1
    return [e for e, c in cnt.items() if c == 1]


def loops(edges):
    adj = defaultdict(list)
    for a, b in edges:
        adj[a].append(b)
        adj[b].append(a)
    seen, out = set(), []
    for start in adj:
        if start in seen:
            continue
        loop, prev, cur = [start], None, start
        seen.add(start)
        while True:
            nxt = [v for v in adj[cur] if v != prev]
            if not nxt or nxt[0] == start:
                break
            prev, cur = cur, nxt[0]
            if cur in seen:
                break
            seen.add(cur)
            loop.append(cur)
        out.append(loop)
    return out, adj


def main():
    meshdir = sys.argv[1]
    s_edge = float(sys.argv[2]) if len(sys.argv) > 2 and sys.argv[2] != "-" else nml_s_edge(meshdir)
    out = sys.argv[3] if len(sys.argv) > 3 else f"{meshdir}/region.png"
    P, T, s, reg = load(meshdir)
    flag = region_flags(T, s, reg, s_edge)
    allb = set(boundary_edges(T, np.ones(len(T), bool)))
    rb = boundary_edges(T, flag)
    outer = [e for e in rb if e in allb]
    inner = [e for e in rb if e not in allb]
    print(f"{meshdir}: s_edge = {s_edge}")
    print(f"  triangles in region: {flag.sum()} of {len(T)} (core {np.sum(flag & (reg == 0))}, "
          f"edge {np.sum(flag & (reg == 1))}, axis fill {np.sum(flag & (reg == 2))})")
    print(f"  region boundary: {len(outer)} edges on the wall, {len(inner)} inner edges")
    ok = True
    if s_edge > 0 and inner:
        sv = s[np.array(inner)]
        dev = np.abs(sv - s_edge).max()
        print(f"  inner boundary vertices: max |s - s_edge| = {dev:.2e}")
        ok &= dev < 1e-12
        lps, adj = loops(inner)
        deg = {len(v) for v in adj.values()}
        print(f"  inner boundary: {len(lps)} loop(s), vertex degrees {sorted(deg)}, "
              f"{len(adj)} vertices")
        ok &= len(lps) == 1 and deg == {2}
        # ring gaps around s_edge (radial spacing of the core rings)
        ring_s = np.unique(np.round(s[s > 0], 14))
        k = np.searchsorted(ring_s, s_edge)
        gaps = np.diff(ring_s[max(k - 3, 0): k + 4])
        print(f"  ring gaps around s_edge: {np.array2string(gaps, precision=5)}")
        # smallest triangle angle in the core rings next to the boundary
    elif s_edge > 0:
        ok = False
    print("  CHECK", "OK" if ok else "FAILED")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    fig, axs = plt.subplots(1, 2, figsize=(11, 7), gridspec_kw={"width_ratios": [1, 1]})
    for ax, zoom in zip(axs, (False, True)):
        ax.triplot(P[:, 0], P[:, 1], T, lw=0.15, color="0.6")
        ax.tripcolor(P[:, 0], P[:, 1], T, facecolors=flag.astype(float), cmap="Blues", vmin=-0.3, vmax=1.6,
                     alpha=0.6, edgecolors="none")
        if inner:
            ax.add_collection(LineCollection(P[np.array(inner)], colors="C3", lw=1.2, label="inner boundary"))
        ax.add_collection(LineCollection(P[np.array(outer)], colors="k", lw=1.0, label="wall"))
        ax.set_aspect("equal")
        ax.set_xlabel("R [cm]")
        ax.set_ylabel("Z [cm]")
        if zoom and inner:
            c = P[np.array(inner)].reshape(-1, 2)
            i = np.argmax(c[:, 0])
            ax.set_xlim(c[i, 0] - 8, c[i, 0] + 4)
            ax.set_ylim(c[i, 1] - 6, c[i, 1] + 6)
            ax.set_title("outer midplane, zoom")
        else:
            ax.set_title(f"edge region (shaded), s_edge = {s_edge:g}")
            ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    fig.savefig(out, dpi=130)
    print("  plot:", out)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
