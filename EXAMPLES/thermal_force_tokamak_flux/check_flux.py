#!/usr/bin/env python3
"""
Checks of the tokamak driver with flux-surface profiles in GORILLA's example-8
circular tokamak (R0 = 1.7 m, a = 0.5 m, B0 = 2 T, q = 1.1 ... 3.1), W3+ in H+.

  mesh       screening drift dr/dt over 0..2 tau_s (paired difference) and the
             radial spread of the control run, for three mesh resolutions

Usage: python3 check_flux.py [work_dir]   (default ./run)
"""
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "thermal_force_tokamak"))
import check_tokamak as ck  # noqa: E402  (run helper)

FILES = ("gorilla_applets.inp", "gorilla.inp", "tetra_grid.inp", "field_divB0.inp", "preload_for_SYNCH.inp")
R0, R_S, R_MINOR, B0 = 170.0, 195.0, 25.0, 2.0e4


def ds_dr(r, R0=170.0, a=50.0):
    """d s/d r [1/cm] for GORILLA's analytic circular tokamak, where the flux label is
    s = (R0 - sqrt(R0^2 - r^2))/(R0 - sqrt(R0^2 - a^2)) (field_analytic_circ_mod)."""
    return r / math.sqrt(R0**2 - r**2) / (R0 - math.sqrt(R0**2 - a**2))


def template_with_mesh(work, n1, n2, n3):
    t = work / f"template_{n1}_{n2}_{n3}"
    t.mkdir(parents=True, exist_ok=True)
    for f in FILES + ("thermal_force.inp",):
        if (HERE / f).exists():
            (t / f).write_text((HERE / f).read_text())
    grid = (t / "tetra_grid.inp").read_text()
    for key, val in (("n1", n1), ("n2", n2), ("n3", n3)):
        lines = [f" {key:<24}= {val} ," if l.strip().startswith(key + " ") else l for l in grid.splitlines()]
        grid = "\n".join(lines) + "\n"
    (t / "tetra_grid.inp").write_text(grid)
    if not (t / "profiles").exists():
        (t / "profiles").symlink_to(HERE / "profiles")
    return t


def run_mesh(work, mesh, n=4000):
    t = template_with_mesh(work, *mesh)
    case = work / "case"
    rows = ck.run(case, t, FILES, {"output_time_unit": "tau_s", "output_times": [0.0, 2.0], "n_particles": n,
                                   "output_file": "out.csv"})
    return rows


def main():
    work = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE / "run"
    work.mkdir(parents=True, exist_ok=True)
    print("mesh (n1, n2, n3)   screening dr/dt [m/s]     control: radial spread growth var(r)/t [cm^2/s]")
    for mesh in ((20, 40, 20), (40, 80, 40), (80, 160, 80)):
        rows = run_mesh(work, mesh)
        d = [r for r in rows if r["run"] == "difference"][-1]
        c0, c1 = [r for r in rows if r["run"] == "control"][0], [r for r in rows if r["run"] == "control"][-1]
        t = d["t_s"]
        k = ds_dr(R_MINOR)
        v = d["mean_s"] / k / 100.0 / t
        se = d["se_s"] / k / 100.0 / t
        var0 = (c0["se_s"] * math.sqrt(4000) / k) ** 2
        var1 = (c1["se_s"] * math.sqrt(4000) / k) ** 2
        print(f"{str(mesh):18s}  {v:8.4f} +- {se:.4f}        {(var1 - var0) / t:10.1f}")


if __name__ == "__main__":
    main()


def flat_template(work):
    """Example mesh with flat profiles (T = 50 eV, n = 1e14 cm^-3): pure diffusion."""
    t = template_with_mesh(work / "flat", 40, 80, 40)
    prof = t / "profiles"
    if prof.is_symlink():
        prof.unlink()
    prof.mkdir(exist_ok=True)
    for name, val in (("n.dat", 1.0e14), ("Te.dat", 50.0), ("Ti.dat", 50.0), ("Er.dat", 0.0)):
        with open(prof / name, "w") as f:
            f.write("# r_eff [cm], value\n")
            for r in range(0, 51, 5):
                f.write(f"{float(r):.6e} {val:.6e}\n")
    return t


def diffusion(work, shift, t_end=300.0, n=2000, dt=0.02):
    """Radial diffusion coefficient D_r [m^2/s] from var(r)(t) over the second half of the run."""
    t = flat_template(work)
    times = [0.0] + [t_end * k / 12 for k in range(1, 13)]
    rows = ck.run(work / "case", t, FILES, {"output_time_unit": "tau_s", "output_times": times, "n_particles": n,
                                            "boole_control_run": False, "boole_collision_shift": shift,
                                            "dt_coll_over_tau_s": dt, "output_file": "out.csv"})
    rows = [r for r in rows if r["run"] == "gradient"]
    k = ds_dr(R_MINOR)
    ts = [r["t_s"] for r in rows]
    var = [(r["se_s"] * math.sqrt(n) / k / 100.0) ** 2 for r in rows]          # m^2
    mean = [r["mean_s"] / k / 100.0 for r in rows]
    half = len(ts) // 2
    import numpy as np
    slope_var = np.polyfit(ts[half:], var[half:], 1)[0]
    slope_mean = np.polyfit(ts[half:], mean[half:], 1)[0]
    return slope_var / 2.0, slope_mean, ts, var
