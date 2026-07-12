#!/usr/bin/env python3
"""
CTest driver for the helical-core CI smoke test.

Reads GORILLA-core blueprints from ../GORILLA/INPUT/ and applet blueprints
from GORILLA_APPLETS/INPUT/, sets every key the test depends on explicitly,
then invokes the binary with i_option = 14 (helical core particle tracing).

The test uses the shared field-aligned EFIT test equilibrium
(grid_kind = 2 with MHD_EQUILIBRIA/g_file_for_test and the SYNCH preload),
the configuration the shipped example runs; the rectangular EFIT grid
aborts inside the option-14 volume integrals. It covers the option-14
execution path, not helical-core physics: the perturbation stays off,
matching the shipped example.

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


def env_path(name: str) -> Path:
    value = os.environ.get(name)
    if not value:
        sys.exit(f"missing required env var: {name}")
    return Path(value)


APPLETS_ROOT = env_path("APPLETS_ROOT")
GORILLA_ROOT = env_path("GORILLA_ROOT")
BINARY = env_path("GORILLA_APPLETS_BIN")
WORK_DIR = env_path("WORK_DIR")

# Fresh work dir for every test run
if WORK_DIR.exists():
    shutil.rmtree(WORK_DIR)
WORK_DIR.mkdir(parents=True)

# Load blueprints from their conceptual owners
gorilla = f90nml.read(str(GORILLA_ROOT / "INPUT" / "gorilla.inp"))
tetra_grid = f90nml.read(str(GORILLA_ROOT / "INPUT" / "tetra_grid.inp"))
gorilla_applets = f90nml.read(str(APPLETS_ROOT / "INPUT" / "gorilla_applets.inp"))
helical_core = f90nml.read(str(APPLETS_ROOT / "INPUT" / "helical_core.inp"))
for nml in (gorilla, tetra_grid, gorilla_applets, helical_core):
    nml.end_comma = True

# Core integrator (cylindrical coords, polynomial pusher, cheap order for CI)
gorilla["gorillanml"]["eps_phi"] = 0.0
gorilla["gorillanml"]["coord_system"] = 1            # (R, phi, Z)
gorilla["gorillanml"]["ispecies"] = 2                # deuterium
gorilla["gorillanml"]["boole_periodic_relocation"] = True
gorilla["gorillanml"]["ipusher"] = 2                 # polynomial pusher
gorilla["gorillanml"]["poly_order"] = 2
gorilla["gorillanml"]["boole_guess"] = True
gorilla["gorillanml"]["boole_grid_for_find_tetra"] = False
gorilla["gorillanml"]["boole_adaptive_time_steps"] = False

# Tetrahedral grid: field-aligned EFIT as in the shipped example, small
tetra_grid["tetra_grid_nml"]["grid_kind"] = 2        # field-aligned EFIT
tetra_grid["tetra_grid_nml"]["n1"] = 20
tetra_grid["tetra_grid_nml"]["n2"] = 30
tetra_grid["tetra_grid_nml"]["n3"] = 50
tetra_grid["tetra_grid_nml"]["boole_n_field_periods"] = True
tetra_grid["tetra_grid_nml"]["g_file_filename"] = "MHD_EQUILIBRIA/g_file_for_test"
tetra_grid["tetra_grid_nml"]["convex_wall_filename"] = "MHD_EQUILIBRIA/convex_wall_for_test.dat"

# Applets dispatcher
gorilla_applets["gorilla_applets_nml"]["i_option"] = 14

# Helical core: few antithetic particle pairs, short trace, delta-f weights
helical_core["helical_core_nml"]["time_step"] = 1.0e-4
helical_core["helical_core_nml"]["energy_ev"] = 1.0e3
helical_core["helical_core_nml"]["n_species"] = 1
helical_core["helical_core_nml"]["n_particles"] = 4
helical_core["helical_core_nml"]["boole_squared_moments"] = False
helical_core["helical_core_nml"]["boole_point_source"] = False
helical_core["helical_core_nml"]["boole_collisions"] = True
helical_core["helical_core_nml"]["boole_precalc_collisions"] = False
helical_core["helical_core_nml"]["boole_refined_sqrt_g"] = True
helical_core["helical_core_nml"]["boole_monoenergetic"] = True
helical_core["helical_core_nml"]["boole_linear_density_simulation"] = False
helical_core["helical_core_nml"]["boole_linear_temperature_simulation"] = False
helical_core["helical_core_nml"]["boole_eliminate_particles_outside_flux"] = False
helical_core["helical_core_nml"]["boole_delta_f"] = True
helical_core["helical_core_nml"]["boole_antithetic_variate"] = True
helical_core["helical_core_nml"]["i_integrator_type"] = 1
helical_core["helical_core_nml"]["seed_option"] = 2
helical_core["helical_core_nml"]["boole_write_moments"] = True
helical_core["helical_core_nml"]["boole_write_fourier_moments"] = True
helical_core["helical_core_nml"]["boole_write_exit_data"] = True
helical_core["helical_core_nml"]["boole_write_grid_data"] = True

# Deterministic seed for reproducible CI output
(WORK_DIR / "seed.inp").write_text("8\n  1 2 3 4 5 6 7 8\n")

# Axisymmetric equilibrium only (ipert = 0 in the upstream blueprint);
# the field-aligned grid additionally needs the SYNCH preload namelist.
shutil.copy(GORILLA_ROOT / "INPUT" / "field_divB0.inp", WORK_DIR / "field_divB0.inp")
shutil.copy(GORILLA_ROOT / "INPUT" / "preload_for_SYNCH.inp",
            WORK_DIR / "preload_for_SYNCH.inp")

# Symlink the shared test equilibrium into the work dir
(WORK_DIR / "MHD_EQUILIBRIA").symlink_to(
    GORILLA_ROOT / "MHD_EQUILIBRIA", target_is_directory=True)

# Write namelists
gorilla.write(str(WORK_DIR / "gorilla.inp"), force=True)
tetra_grid.write(str(WORK_DIR / "tetra_grid.inp"), force=True)
gorilla_applets.write(str(WORK_DIR / "gorilla_applets.inp"), force=True)
helical_core.write(str(WORK_DIR / "helical_core.inp"), force=True)

print("=== helical core: i_option=14 ===", flush=True)
result = subprocess.run(
    [str(BINARY)],
    cwd=WORK_DIR,
    capture_output=True,
    text=True,
    check=False,
)
print(result.stdout)
if result.stderr:
    print(result.stderr, file=sys.stderr)
if result.returncode != 0:
    sys.exit(f"FAIL: gorilla_applets_main.x exited with code {result.returncode}")

# Smoke check: the final summary stanza must be present
expected_marker = "Number of lost ions"
if expected_marker not in result.stdout:
    sys.exit(f"FAIL: expected '{expected_marker}' in stdout but it was missing")

print("PASS: helical core tracing completed and summary marker found")
