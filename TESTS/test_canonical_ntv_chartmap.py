#!/usr/bin/env python3
"""Gates of the Boozer chartmap route (grid_kind 6) of the canonical NTV estimator.

The chartmap is the exact Boozer chart of GORILLA's analytic circular field
(write_circular_boozer_chartmap.x), so the analytic route (grid_kind 5,
cylindrical coordinates, harmonic in the straight-field-line angle) and the
chartmap route describe the same physics. Harmonic dB1 = b cos(m theta* - n
phi), m = 3, n = 2; profiles with n and T gradients; 1 keV-scale electrons.

Gates:
  (a) chartmap, sigma = 0, no collisions: delta-f weights exactly zero; the
      mesh p_phi drift as |Delta s*| (reported; G's grid_kind 5 numbers for
      comparison) below DRIFT_SSTAR_MAX.
  (b) chartmap: harmonics from the perturbed chartmap file (Fourier
      projection) give the same delta-f torque as the analytic harmonic of
      the same amplitude (same markers): relative difference < 1e-6.
  (c) analytic route vs chartmap route, same harmonic, MH collisions: the
      delta-f torque estimates V/N sum P tq_deltaf agree within 4 combined
      standard errors (different marker sets).
  (d) full-f on the perturbed chartmap meshes (+sigma, -sigma, 0; common
      starts and per-marker collision streams) vs delta-f on the
      axisymmetric mesh with the harmonic of the +sigma file:
      linear order, per marker: the odd part of the full-f p_phi change and of
      int -dH/dphi|mesh dt correlates with int -dH1/dphi dt of delta-f
      (corr > 0.98, slope 0.85..1.05; reported for the total p_phi change,
      whose collisional part decorrelates; the mesh resolves the m = 3
      harmonic with n3 = 24 points, which lowers its derivative by ~5%);
      quadratic order, statistically: [dp(+) + dp(-) - 2 dp(0)]/2 agrees with
      the delta-f torque within 4 combined standard errors. At this size the
      full-f error is ~100x the delta-f one (phase-mixed second-order orbit
      response), so the quadratic gate is weak; NTV_N_PARTICLES and
      NTV_TIME_STEP scale it up.
  (e) optional wall contour (loss_wall_file, circle r = 60 cm inside the
      spawn window): markers born outside are lost there (loss_kind 2),
      markers well inside are not, and every loss lies at the wall.

Paths come from CMake (APPLETS_ROOT, GORILLA_APPLETS_BIN, CHARTMAP_WRITER_BIN,
WORK_DIR). NTV_N_PARTICLES, NTV_TIME_STEP, NTV_N2 and NTV_GATES ('abcde')
override the case.
"""

import math
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

DENSITY = 1.2e14
TEMPERATURE_EV = 1000.0
A_MINOR = 69.991288
DRIFT_SSTAR_MAX = 2.0e-3
M_POL, N_TOR, B_MN = 3, 2, 10.0
FULLF_BMN = 40.0
N_SIGMA = 4.0


def env_path(name):
    value = os.environ.get(name)
    if not value:
        sys.exit(f"missing environment variable {name}")
    return Path(value)


APPLETS_ROOT = env_path("APPLETS_ROOT")
BIN = env_path("GORILLA_APPLETS_BIN")
WRITER = env_path("CHARTMAP_WRITER_BIN")
WORK_DIR = env_path("WORK_DIR")
TEMPLATE = APPLETS_ROOT / "TESTS" / "rmp_fm_stationarity"


def set_nml(path, values):
    text = path.read_text()
    for key, value in values.items():
        text, count = re.subn(
            rf"(^\s*{key}\s*=\s*)[^,\n]*", rf"\g<1>{value}", text, flags=re.M
        )
        if count == 0:
            text = re.sub(r"^\s*/\s*$", f" {key} = {value} ,\n/\n", text,
                          count=1, flags=re.M)
        elif count != 1:
            sys.exit(f"{path.name}: expected one {key} entry, found {count}")
    path.write_text(text)


def write_chartmap(path, sigma_b):
    subprocess.run([str(WRITER), str(path), "1.0", repr(sigma_b), str(M_POL),
                    str(N_TOR), "33", "64", "48"], check=True,
                   stdout=subprocess.DEVNULL)


def prepare(case_dir, chart, n_particles, time_step, collisions, ntv, mesh=None,
            n2=24):
    if case_dir.exists():
        shutil.rmtree(case_dir)
    shutil.copytree(TEMPLATE, case_dir)
    profiles = case_dir / "profiles"
    profiles.mkdir()
    radii = [0.1 + 0.25 * i for i in range(300)]
    for name, value in (("n.dat", DENSITY), ("Te.dat", TEMPERATURE_EV),
                        ("Ti.dat", TEMPERATURE_EV)):
        lines = [f"{r:.6e} {value * (1.0 - 0.6 * (r / A_MINOR) ** 2):.6e}"
                 for r in radii]
        (profiles / name).write_text("\n".join(lines) + "\n")
    (profiles / "Er.dat").write_text(
        "\n".join(f"{r:.6e} 0.0" for r in radii) + "\n")
    (case_dir / "seed.inp").write_text(
        "33\n" + " ".join(str(1000 + 7 * k) for k in range(33)) + "\n")
    # The refined sqrt(g) prism volumes are cylindrical-only (moments only).
    set_nml(case_dir / "rmp_response_currents.inp",
            {"n_particles": n_particles, "time_step": time_step,
             "boole_collisions": ".true." if collisions else ".false.",
             "boole_nonlinear_weight": ".false.",
             "boole_refined_sqrt_g": ".true." if chart == "analytic" else ".false."})
    gorilla = {"boole_time_hamiltonian": ".true."}
    grid = {"n2": n2}
    if chart == "chartmap":
        gorilla["coord_system"] = 2
        # E_r = 0: the prescribed-E_r loader of GORILLA serves EFIT/analytic
        # grids only (it evaluates the EFIT field).
        gorilla["er_profile_file"] = "''"
        # netcdf_filename holds 100 characters: link the mesh into the case.
        (case_dir / "mesh.nc").symlink_to(Path(mesh).resolve())
        grid.update({"grid_kind": 6, "netcdf_filename": "'mesh.nc'"})
    set_nml(case_dir / "gorilla.inp", gorilla)
    set_nml(case_dir / "tetra_grid.inp", grid)
    entries = {"boole_canonical_ntv": ".true.", "ntv_m_pol": M_POL,
               "ntv_n_tor": N_TOR, "ntv_theta_kind": 1}
    entries.update(ntv)
    body = "".join(f"  {k} = {v} ,\n" for k, v in entries.items())
    with open(case_dir / "rmp_response_currents.inp", "a") as f:
        f.write(f"\n&canonical_ntv_nml\n{body}/\n")


def run(case_dir, threads):
    env = dict(os.environ, OMP_NUM_THREADS=str(threads))
    with open(case_dir / "run.log", "w") as log:
        result = subprocess.run([str(BIN)], cwd=case_dir, env=env,
                                stdout=log, stderr=subprocess.STDOUT)
    if result.returncode != 0:
        sys.exit(f"applet failed in {case_dir} (rc={result.returncode})")
    markers = read_markers(case_dir / "canonical_ntv_markers.dat")
    markers["volume"] = read_volume(case_dir / "canonical_ntv_markers.dat")
    return markers


def read_volume(path):
    for line in path.read_text().splitlines():
        if line.startswith("# spawn_volume_cm3"):
            return float(line.split()[2])
    sys.exit(f"{path}: no spawn volume")


def read_markers(path):
    names = None
    rows = []
    for line in path.read_text().splitlines():
        if line.startswith("# n "):
            names = line[2:].split()
        elif line.strip() and not line.startswith("#"):
            rows.append([float(v) for v in line.split()])
    if names is None or not rows:
        sys.exit(f"{path}: no marker table")
    return {k: [r[i] for r in rows] for i, k in enumerate(names)}


def mean_err(values):
    n = len(values)
    mean = sum(values) / n
    var = sum((v - mean) ** 2 for v in values) / (n - 1)
    return mean, math.sqrt(var / n)


def corr_slope(y, x):
    mx, my = sum(x) / len(x), sum(y) / len(y)
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    sxx = sum((a - mx) ** 2 for a in x)
    syy = sum((b - my) ** 2 for b in y)
    return sxy / math.sqrt(sxx * syy), sum(a * b for a, b in zip(x, y)) / sum(
        a * a for a in x)


def torque_samples(m):
    """Per-marker samples of V sum_i P_i tq_i / N (mean = estimate)."""
    return [m["volume"] * p * t for p, t in zip(m["P"], m["tq_deltaf"])]


def main():
    WORK_DIR.mkdir(parents=True, exist_ok=True)
    threads = int(os.environ.get("OMP_NUM_THREADS", "4"))
    npart = int(os.environ.get("NTV_N_PARTICLES", "512"))
    dt = os.environ.get("NTV_TIME_STEP", "2.0d-6")
    n2 = int(os.environ.get("NTV_N2", "24"))
    gates = os.environ.get("NTV_GATES", "abcde")
    failures = []

    zero = WORK_DIR / "circ_zero.nc"
    pert = WORK_DIR / "circ_b10.nc"
    write_chartmap(zero, 0.0)
    write_chartmap(pert, B_MN)

    def case(tag, chart, coll, ntv, mesh=zero, n=npart):
        d = WORK_DIR / tag
        prepare(d, chart, n, dt, coll, ntv, mesh=mesh, n2=n2)
        return run(d, threads)

    analytic = {"ntv_amplitude": 1.0, "ntv_bmn": B_MN}
    from_file = {"ntv_amplitude": 1.0, "ntv_chartmap_file": f"'{pert}'"}

    if "a" in gates:
        m0 = case("a_sigma0", "chartmap", False, {"ntv_amplitude": 0.0})
        max_w = max(abs(w) for w in m0["w_end"])
        ds = [abs(b - a) for a, b, lost in
              zip(m0["sstar_birth"], m0["sstar_end"], m0["lost"]) if lost == 0]
        rms = math.sqrt(sum(d * d for d in ds) / len(ds))
        print(f"(a) n2={n2}, t={dt}: max |w| = {max_w:.3e}; p_phi drift |Delta s*|: "
              f"max {max(ds):.3e}, rms {rms:.3e} over {len(ds)} markers "
              f"(G, grid_kind 5 at 2e-4 s: rms 5.7e-5/1.3e-5/4.2e-6 for n2=24/48/96)")
        if max_w != 0.0:
            failures.append(f"unperturbed weights nonzero ({max_w:.3e})")
        if max(ds) > DRIFT_SSTAR_MAX:
            failures.append(f"p_phi drift on the chartmap mesh too large ({max(ds):.3e})")

    if "b" in gates or "c" in gates:
        mc = case("b_chartmap_analytic", "chartmap", True, analytic)
    if "b" in gates:
        mf = case("b_chartmap_file", "chartmap", True, from_file)
        t_a = sum(p * t for p, t in zip(mc["P"], mc["tq_deltaf"]))
        t_f = sum(p * t for p, t in zip(mf["P"], mf["tq_deltaf"]))
        rel = abs(t_f - t_a) / abs(t_a)
        print(f"(b) sum P tq: analytic harmonic {t_a:.10e}, chartmap file {t_f:.10e}, "
              f"rel diff {rel:.2e}")
        if not rel < 1.0e-6:
            failures.append(f"chartmap-file harmonic differs from analytic ({rel:.2e})")

    if "c" in gates:
        m5 = case("c_analytic_route", "analytic", True, analytic)
        t5, e5 = mean_err(torque_samples(m5))
        t6, e6 = mean_err(torque_samples(mc))
        z = abs(t5 - t6) / math.hypot(e5, e6)
        print(f"(c) V/N sum P tq: analytic route {t5:.4e} +- {e5:.1e}, chartmap "
              f"{t6:.4e} +- {e6:.1e}; ratio {t6 / t5:.3f}, |diff|/err {z:.2f}; "
              f"volumes {m5['volume']:.6e} {mc['volume']:.6e}")
        if not z < N_SIGMA:
            failures.append(f"analytic and chartmap routes disagree ({z:.2f} sigma)")
        if abs(mc["volume"] / m5["volume"] - 1.0) > 3.0e-3:
            failures.append("chartmap J0 volume differs from the analytic volume")

    if "d" in gates:
        meshes = {}
        for label, b in (("p", FULLF_BMN), ("m", -FULLF_BMN)):
            meshes[label] = WORK_DIR / f"circ_full_{label}.nc"
            write_chartmap(meshes[label], b)
        meshes["0"] = zero
        coll = os.environ.get("NTV_FULLF_COLLISIONS", "1") == "1"
        full = {k: case(f"d_fullf_{k}", "chartmap", coll, {"ntv_amplitude": 0.0},
                        mesh=v) for k, v in meshes.items()}
        df = case("d_deltaf", "chartmap", coll,
                  {"ntv_amplitude": 1.0, "ntv_chartmap_file": f"'{meshes['p']}'"})
        same_start = all(full[k]["s_birth"] == full["0"]["s_birth"] for k in "pm")
        if not same_start:
            failures.append("full-f runs do not share marker starts")
        vol = full["0"]["volume"]

        def weighted(m, key):
            if key == "dp":
                return [p * (b - a) for p, a, b in
                        zip(m["P"], m["pphi_birth"], m["pphi_end"])]
            return [p * v for p, v in zip(m["P"], m[key])]

        # Quadratic part of the orbit tally int -dH/dphi|mesh dt (the NTV part;
        # the collisional p_phi change is a separate tally) and, for reference,
        # of the total p_phi change.
        est = {}
        for key in ("dh_mesh", "dp"):
            a, b, c = (weighted(full[k], key) for k in "pm0")
            est[key] = mean_err([vol * 0.5 * (x + y - 2.0 * w)
                                 for x, y, w in zip(a, b, c)])
        tq, eq = est["dh_mesh"]
        td, ed = mean_err(torque_samples(df))
        z = abs(tq - td) / math.hypot(eq, ed)
        lost = [sum(m["lost"]) for m in full.values()]
        print(f"(d) full-f quadratic part of V P dh_mesh {tq:.4e} +- {eq:.1e} (of V P dp "
              f"{est['dp'][0]:.4e} +- {est['dp'][1]:.1e}), delta-f {td:.4e} +- {ed:.1e}; "
              f"|diff|/err {z:.2f}; lost (+,-,0) {lost}")
        if not z < N_SIGMA:
            failures.append(f"full-f and delta-f disagree ({z:.2f} sigma)")
        # Linear order per marker: the odd part of the full-f p_phi change and of
        # int -dH/dphi|mesh dt against int -dH1/dphi dt of delta-f (same starts).
        lin = df["dh_int"]
        # The total p_phi change includes the collisional tally, whose second
        # runs lose alignment as the orbits separate: reported only.
        for key, cmin in (("dh_mesh", 0.98), ("dp", None)):
            if key == "dp":
                odd_i = [0.5 * ((a1 - a0) - (b1 - b0)) for a0, a1, b0, b1 in
                         zip(full["p"]["pphi_birth"], full["p"]["pphi_end"],
                             full["m"]["pphi_birth"], full["m"]["pphi_end"])]
            else:
                odd_i = [0.5 * (a - b) for a, b in zip(full["p"][key], full["m"][key])]
            c, slope = corr_slope(odd_i, lin)
            print(f"    linear order per marker, {key}: corr {c:.5f}, slope {slope:.4f}")
            if cmin is not None and not (c > cmin and 0.85 < slope < 1.05):
                failures.append(f"full-f linear response ({key}) differs from delta-f")

    if "e" in gates:
        # Wall: circle r_w = 60 cm around the axis inside the spawn window.
        r0, r_w = 165.0, 60.0
        s_edge = r0 - math.sqrt(r0 ** 2 - A_MINOR ** 2)
        s_w = (r0 - math.sqrt(r0 ** 2 - r_w ** 2)) / s_edge
        wall = WORK_DIR / "wall_circle.dat"
        wall.write_text("".join(
            f"{r0 + r_w * math.cos(2 * math.pi * k / 256):.8f} "
            f"{r_w * math.sin(2 * math.pi * k / 256):.8f}\n" for k in range(256)))
        mw = case("e_wall", "chartmap", True,
                  {"ntv_amplitude": 0.0, "loss_wall_file": f"'{wall}'"})
        rows = list(zip(mw["s_birth"], mw["s_end"], mw["lost"], mw["loss_kind"]))
        n_out = sum(1 for sb, _, _, _ in rows if sb > s_w + 0.01)
        bad = [r for r in rows if (r[0] > s_w + 0.01 and r[3] != 2)
               or (r[0] < s_w - 0.02 and r[2] != 0)
               or (r[2] != 0 and (r[3] != 2 or r[1] < s_w - 0.01))]
        print(f"(e) wall r = {r_w} cm (s = {s_w:.4f}): {int(sum(mw['lost']))} lost, "
              f"{n_out} born outside; inconsistent markers: {len(bad)}")
        if n_out == 0 or bad:
            failures.append("wall contour loss boundary inconsistent")

    if failures:
        print("FAIL: " + "; ".join(failures))
        return 1
    print("PASS: canonical NTV chartmap route gates")
    return 0


if __name__ == "__main__":
    sys.exit(main())
