#!/usr/bin/env python3
"""
CTest driver: radial potential Phi0(s) on the vertices of a Boozer chartmap
mesh (grid_kind = 6, flux_mesh_potential_mod).

The chartmap is GORILLA's analytic circular tokamak written in Boozer
coordinates (write_circular_boozer_chartmap.x). With the vertex potential
Phi0(s) = p0 (s - 0.8 s^2) the ExB drift changes the field-line label
alpha = theta - iota*phi at the rate
    dalpha/dt = c dPhi0/ds / torflux = -iota Omega_E,   Omega_E = c dPhi0/dA_phi,
on every surface (A_theta = torflux s, dA_phi/ds = -iota torflux). Omega_E
changes sign at s = 0.625. The potential is linear in each tetrahedron, so
the mesh rate is the secant slope of the radial cell of the orbit.

Checks, for deeply trapped (no streaming) 1 eV deuterons on several surfaces,
with the rate of the same orbit without potential subtracted:
  1. the measured rate equals the cell secant c dPhi0/ds / torflux (2e-3; the
     residual is the angular mesh error, at s0 = 0.4: 1.6e-3 for n2 = n3 = 30,
     2e-4 for 60; at s0 = 0.2 1.2e-3 for 60);
  2. Omega_E = -rate/iota equals the analytic c Phi0'(s0)/(dA_phi/ds) within
     the linear-element error |Phi0''| ds_cell / 2 (+1e-3);
  3. Phi0 = 0 via set_flux_mesh_potential leaves the orbit bitwise unchanged.
"""
import os
import shutil
import subprocess
import sys
from pathlib import Path

try:
    import f90nml
except ImportError:
    sys.exit("f90nml is required to run this test (pip install f90nml)")

P0 = 2.0  # statV
ENERGY_EV = 1.0
SURFACES = (0.2, 0.4, 0.8)
T_STEP, N_STEPS = 2.0e-6, 400
REL_TOL = 1.0e-3
CELL_TOL = 2.0e-3


def env_path(name):
    value = os.environ.get(name)
    if not value:
        sys.exit(f"missing required env var: {name}")
    return Path(value)


GORILLA_ROOT = env_path("GORILLA_ROOT")
DRIVER = env_path("DRIVER_BIN")
WRITER = env_path("CHARTMAP_WRITER_BIN")
WORK_DIR = env_path("WORK_DIR")

if WORK_DIR.exists():
    shutil.rmtree(WORK_DIR)
WORK_DIR.mkdir(parents=True)
subprocess.run([str(WRITER), "chartmap.nc", "0.0", "0.0", "3", "2", "65", "64", "48"],
               cwd=WORK_DIR, check=True, stdout=subprocess.DEVNULL)

gorilla = f90nml.read(str(GORILLA_ROOT / "INPUT" / "gorilla.inp"))
tetra_grid = f90nml.read(str(GORILLA_ROOT / "INPUT" / "tetra_grid.inp"))
for nml in (gorilla, tetra_grid):
    nml.end_comma = True
g = gorilla["gorillanml"]
g["eps_phi"] = 0.0
g["coord_system"] = 2
g["ispecies"] = 2
g["boole_periodic_relocation"] = True
g["ipusher"] = 2
g["poly_order"] = 2
g["boole_grid_for_find_tetra"] = False
g["boole_adaptive_time_steps"] = False
g["er_profile_file"] = ""
t = tetra_grid["tetra_grid_nml"]
t["grid_kind"] = 6
t["n1"] = 100
t["n2"] = 60
t["n3"] = 60
t["boole_n_field_periods"] = True
t["sfc_s_min"] = 1.0e-1
t["netcdf_filename"] = "chartmap.nc"
gorilla.write(str(WORK_DIR / "gorilla.inp"), force=True)
tetra_grid.write(str(WORK_DIR / "tetra_grid.inp"), force=True)
shutil.copy(GORILLA_ROOT / "INPUT" / "field_divB0.inp", WORK_DIR / "field_divB0.inp")


def run(mode, s0):
    out = subprocess.run(
        [str(DRIVER), mode, repr(s0), repr(ENERGY_EV), repr(T_STEP), str(N_STEPS), repr(P0)],
        cwd=WORK_DIR, check=True, capture_output=True, text=True).stdout
    result = {}
    for line in out.splitlines():
        if line.startswith("RESULT"):
            _, key, value = line.split()
            result[key] = float(value)
        elif line.startswith("FINAL"):
            result["final"] = line.split()[1:]
    return result


def phi0(s):
    return P0 * (s - 0.8 * s**2)


failures = []
for s0 in SURFACES:
    ref = run("none", s0)
    zero = run("zero", s0)
    prof = run("profile", s0)
    if zero["final"] != ref["final"]:
        failures.append(f"s0={s0}: Phi0 = 0 changed the orbit {zero['final']} vs {ref['final']}")
    c, torflux, iota = prof["clight"], prof["torflux"], prof["iota"]
    s_lo, s_hi = prof["s_lo"], prof["s_hi"]
    if not (s_lo <= prof["s_min"] and prof["s_max"] <= s_hi):
        failures.append(f"s0={s0}: orbit left its radial cell [{s_lo}, {s_hi}]")
    rate = prof["alpha_rate"] - ref["alpha_rate"]
    secant = c * (phi0(s_hi) - phi0(s_lo)) / (s_hi - s_lo) / torflux
    rel_cell = abs(rate - secant) / abs(secant)
    omega_mesh = -rate / iota
    dphi0 = P0 * (1.0 - 1.6 * s0)
    omega_exact = c * dphi0 / (-iota * torflux)
    bound = 0.8 * P0 * (s_hi - s_lo) / abs(dphi0) + REL_TOL
    rel_exact = abs(omega_mesh - omega_exact) / abs(omega_exact)
    print(f"s0={s0:.3f} cell=[{s_lo:.4f},{s_hi:.4f}] ds_orbit={prof['s_max'] - prof['s_min']:.1e} "
          f"rate(Phi=0)={ref['alpha_rate']:.3e} ExB rate={rate:.6e} secant={secant:.6e} "
          f"rel={rel_cell:.1e}; Omega_E={omega_mesh:.6e} exact={omega_exact:.6e} "
          f"rel={rel_exact:.1e} (bound {bound:.1e})")
    if not rel_cell < CELL_TOL:
        failures.append(f"s0={s0}: ExB rate {rate} differs from cell secant {secant}")
    if not rel_exact < bound:
        failures.append(f"s0={s0}: Omega_E {omega_mesh} differs from analytic {omega_exact}")

if failures:
    sys.exit("FAIL: " + "; ".join(failures))
print("PASS: chartmap mesh potential Phi0(s) gives Omega_E(s); Phi0 = 0 is bitwise neutral")
