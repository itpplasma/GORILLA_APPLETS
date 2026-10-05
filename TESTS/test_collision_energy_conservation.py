#!/usr/bin/env python3
"""
CTest driver: energy and momentum conservation of the conserving background
collision operator (boole_preserve_energy_and_momentum_during_collisions).

Runs gorilla_applets_main.x in the Boltzmann test (i_option = 9) with alpha
test particles (4 amp, Z = 2) colliding with a deuterium background (2 amp;
electrons are not updated by the conserving operator). The opt-in diagnostic
boole_write_collision_energy_diagnostics writes collision_energy_diagnostics.dat.

Conservation law as implemented, per background ion species j and tetrahedron
k, with r = c%weight_factor * w_n * 1e-4 (an internal normalisation, not a
physical particle count):

    E_res = sum_k (3/2) T_jk + (1/2) m_j u_jk^2
    every collision:  dE_res = -r * d_eps_marker,  dP_res = -r * m_t * dvpar

The check compares the reservoir change with the exchange accumulated at
collision time with the r actually used (the "collision ledger"):

    rel_err_E = |dE_res + sum r d_eps| / sum r |d_eps|
    rel_err_P = |dP_res + sum r m_t dvpar| / sum r m_t |dvpar|

The denominators are the exchanged amounts, not the total reservoir energy:
the reservoir holds sum (3/2) T over all cells, which is orders of magnitude
larger than the exchange and would hide any error.

Tolerance: the ledger does not involve the orbit pusher, so the only error is
floating-point round-off. Reservoir changes are summed per cell, and the eV <->
erg round trip of each collision perturbs T by ~1e-16 relative, i.e. ~1e-11 of
a single exchange. TOL = 1e-8 leaves three orders of magnitude of headroom
while any real bookkeeping error (missing mass ratio, erg added to eV) gives
O(1e-2..1).

A second check uses the marker states instead of the ledger, with the final
per-marker ratio r_n:

    rel_err_state = |d(sum_n r_n eps_n) + dE_res| / sum r |d_eps|

It also fails if the first collision of a marker uses a different r than the
later ones (weights not yet initialised before the first collision gave ~3e-4).
Unlike the ledger it includes the pusher's energy drift: a control run with
boole_collisions = .false. (env BOOLE_COLLISIONS=0) gives a marker energy drift
of ~2e-8 of the exchanged energy, so TOL_STATE = 1e-6.

Paths are passed in by TESTS/CMakeLists.txt via environment variables.
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

TOL = 1.0e-8
TOL_STATE = 1.0e-6


def env_path(name: str) -> Path:
    value = os.environ.get(name)
    if not value:
        sys.exit(f"missing required env var: {name}")
    return Path(value)


APPLETS_ROOT = env_path("APPLETS_ROOT")
GORILLA_ROOT = env_path("GORILLA_ROOT")
BINARY = env_path("GORILLA_APPLETS_BIN")
WORK_DIR = env_path("WORK_DIR")
BOOLE_COLLISIONS = os.environ.get("BOOLE_COLLISIONS", "1") == "1"

# Fresh work dir for every test run
if WORK_DIR.exists():
    shutil.rmtree(WORK_DIR)
WORK_DIR.mkdir(parents=True)

# Load blueprints from their conceptual owners
gorilla = f90nml.read(str(GORILLA_ROOT / "INPUT" / "gorilla.inp"))
tetra_grid = f90nml.read(str(GORILLA_ROOT / "INPUT" / "tetra_grid.inp"))
gorilla_applets = f90nml.read(str(APPLETS_ROOT / "INPUT" / "gorilla_applets.inp"))
boltzmann = f90nml.read(str(APPLETS_ROOT / "INPUT" / "boltzmann.inp"))
for nml in (gorilla, tetra_grid, gorilla_applets, boltzmann):
    nml.end_comma = True

# Core integrator: no electrostatic potential, so kinetic energy is the whole
# particle energy
gorilla["gorillanml"]["eps_phi"] = 0.0
gorilla["gorillanml"]["coord_system"] = 1            # (R, phi, Z)
gorilla["gorillanml"]["ispecies"] = 3                # alpha particles
gorilla["gorillanml"]["boole_periodic_relocation"] = True
gorilla["gorillanml"]["ipusher"] = 2                 # polynomial pusher
gorilla["gorillanml"]["poly_order"] = 2
gorilla["gorillanml"]["boole_guess"] = True
gorilla["gorillanml"]["boole_grid_for_find_tetra"] = False
gorilla["gorillanml"]["boole_adaptive_time_steps"] = False

# Tetrahedral grid: rectangular grid on the shared EFIT test equilibrium
tetra_grid["tetra_grid_nml"]["grid_kind"] = 1
tetra_grid["tetra_grid_nml"]["n1"] = 10
tetra_grid["tetra_grid_nml"]["n2"] = 10
tetra_grid["tetra_grid_nml"]["n3"] = 10
tetra_grid["tetra_grid_nml"]["boole_n_field_periods"] = True
tetra_grid["tetra_grid_nml"]["g_file_filename"] = "MHD_EQUILIBRIA/g_file_for_test"
tetra_grid["tetra_grid_nml"]["convex_wall_filename"] = "MHD_EQUILIBRIA/convex_wall_for_test.dat"

# Applets dispatcher: Boltzmann test
gorilla_applets["gorilla_applets_nml"]["i_option"] = 9

# Boltzmann test: few mono-energetic markers, short trace, uniform background
b = boltzmann["boltzmann_nml"]
b["time_step"] = 1.0e-4
b["energy_ev"] = 3.5e3
b["n_particles"] = 32
b["density"] = 1.0e14
b["boole_squared_moments"] = False
b["boole_point_source"] = False
b["boole_collisions"] = BOOLE_COLLISIONS
b["boole_precalc_collisions"] = True
b["boole_preserve_energy_and_momentum_during_collisions"] = True
b["n_background_density_updates"] = 0
b["boole_refined_sqrt_g"] = True
b["boole_monoenergetic"] = True
b["boole_linear_density_simulation"] = False
b["boole_antithetic_variate"] = False
b["boole_linear_temperature_simulation"] = False
b["i_integrator_type"] = 1
b["seed_option"] = 2
b["boole_write_vertex_indices"] = False
b["boole_write_vertex_coordinates"] = False
b["boole_write_prism_volumes"] = False
b["boole_write_refined_prism_volumes"] = False
b["boole_write_boltzmann_density"] = False
b["boole_write_electric_potential"] = False
b["boole_write_moments"] = False
b["boole_write_fourier_moments"] = False
b["boole_write_exit_data"] = True
b["boole_write_grid_data"] = False
b["boole_write_collision_energy_diagnostics"] = True

# Deterministic seed for reproducible output
shutil.copy(APPLETS_ROOT / "INPUT" / "seed.inp", WORK_DIR / "seed.inp")

# Axisymmetric equilibrium only (ipert = 0 in the upstream blueprint)
shutil.copy(GORILLA_ROOT / "INPUT" / "field_divB0.inp", WORK_DIR / "field_divB0.inp")

# Symlink the shared test equilibrium into the work dir
(WORK_DIR / "MHD_EQUILIBRIA").symlink_to(
    GORILLA_ROOT / "MHD_EQUILIBRIA", target_is_directory=True)

# Write namelists
gorilla.write(str(WORK_DIR / "gorilla.inp"), force=True)
tetra_grid.write(str(WORK_DIR / "tetra_grid.inp"), force=True)
gorilla_applets.write(str(WORK_DIR / "gorilla_applets.inp"), force=True)
boltzmann.write(str(WORK_DIR / "boltzmann.inp"), force=True)

# Single thread: the eV <-> erg conversion of the background temperature in the
# collision operator races with unsynchronised reads by other threads
env = dict(os.environ, OMP_NUM_THREADS="1")

print("=== collision energy conservation: i_option=9 ===", flush=True)
result = subprocess.run(
    [str(BINARY)],
    cwd=WORK_DIR,
    env=env,
    capture_output=True,
    text=True,
    check=False,
)
print(result.stdout[-2000:])
if result.stderr:
    print(result.stderr, file=sys.stderr)
if result.returncode != 0:
    sys.exit(f"FAIL: gorilla_applets_main.x exited with code {result.returncode}")

diag_file = WORK_DIR / "collision_energy_diagnostics.dat"
if not diag_file.exists():
    sys.exit("FAIL: collision_energy_diagnostics.dat was not written")
diag = {}
for line in diag_file.read_text().splitlines():
    key, value = line.split("=")
    diag[key.strip()] = float(value)
for key, value in diag.items():
    print(f"{key:34s} {value: .16e}")

marker_drift = abs(diag["marker_energy_final"] - diag["marker_energy_initial"]) \
    / diag["marker_energy_initial"]
print(f"relative marker energy change      {marker_drift: .3e}")

if not BOOLE_COLLISIONS:
    # Control run: only reports the pusher's own energy drift
    print("CONTROL: collisions off, nothing to check")
    sys.exit(0)

if diag["n_collisions"] == 0 or diag["collision_energy_exchange_abs"] == 0.0:
    sys.exit("FAIL: no collisional exchange happened, test setup is broken")
if diag["n_lost_markers"] > 0:
    print(f"note: {int(diag['n_lost_markers'])} marker(s) lost; "
          "they stop colliding and do not affect the ledger")

rel_err_energy = abs(diag["delta_reservoir_energy"] + diag["collision_energy_exchange"]) \
    / diag["collision_energy_exchange_abs"]
rel_err_momentum = abs(diag["delta_reservoir_momentum"] + diag["collision_momentum_exchange"]) \
    / diag["collision_momentum_exchange_abs"]
print(f"rel_err_energy                     {rel_err_energy: .3e}  (tol {TOL:.0e})")
print(f"rel_err_momentum                   {rel_err_momentum: .3e}  (tol {TOL:.0e})")
rel_err_state = abs(diag["marker_energy_final"] - diag["marker_energy_initial"]
                    + diag["delta_reservoir_energy"]) / diag["collision_energy_exchange_abs"]
print(f"rel_err_state                      {rel_err_state: .3e}  (tol {TOL_STATE:.0e})")

failures = []
if not rel_err_energy < TOL:
    failures.append(f"energy not conserved: rel_err_energy = {rel_err_energy:.3e}")
if not rel_err_momentum < TOL:
    failures.append(f"momentum not conserved: rel_err_momentum = {rel_err_momentum:.3e}")
if not rel_err_state < TOL_STATE:
    failures.append(f"marker + reservoir energy not conserved: rel_err_state = {rel_err_state:.3e}")
if failures:
    sys.exit("FAIL: " + "; ".join(failures))

print("PASS: collision reservoir conserves energy and parallel momentum")
