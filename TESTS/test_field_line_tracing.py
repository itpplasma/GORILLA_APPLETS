#!/usr/bin/env python3
"""
CTest driver for the field-line tracing CI test.

Reads GORILLA-core blueprints from ../GORILLA/INPUT/ and applet blueprints
from GORILLA_APPLETS/INPUT/, sets every key the test depends on explicitly,
then invokes the binary with i_option = 10 (field line tracing).

The test uses the bundled ASDEX g-file (grid_kind = 1, axisymmetric EFIT)
so it does not depend on the VMEC equilibrium files from the GORILLA repo.

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

ASDEX_DIR = APPLETS_ROOT / "INPUT" / "data_files" / "DATA" / "ASDEX"

# Fresh work dir for every test run
if WORK_DIR.exists():
    shutil.rmtree(WORK_DIR)
WORK_DIR.mkdir(parents=True)

# Load blueprints from their conceptual owners
gorilla = f90nml.read(str(GORILLA_ROOT / "INPUT" / "gorilla.inp"))
tetra_grid = f90nml.read(str(GORILLA_ROOT / "INPUT" / "tetra_grid.inp"))
gorilla_applets = f90nml.read(str(APPLETS_ROOT / "INPUT" / "gorilla_applets.inp"))
field_line_tracing = f90nml.read(str(APPLETS_ROOT / "INPUT" / "field_line_tracing.inp"))
for nml in (gorilla, tetra_grid, gorilla_applets, field_line_tracing):
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

# Tetrahedral grid: rectangular EFIT, small for CI
tetra_grid["tetra_grid_nml"]["grid_kind"] = 1        # rectangular EFIT
tetra_grid["tetra_grid_nml"]["n1"] = 20              # nR
tetra_grid["tetra_grid_nml"]["n2"] = 4               # nphi
tetra_grid["tetra_grid_nml"]["n3"] = 20              # nZ
tetra_grid["tetra_grid_nml"]["boole_n_field_periods"] = True
tetra_grid["tetra_grid_nml"]["g_file_filename"] = "DATA/ASDEX/g26884.4300"
tetra_grid["tetra_grid_nml"]["convex_wall_filename"] = "DATA/ASDEX/convexwall.dat"

# Applets dispatcher
gorilla_applets["gorilla_applets_nml"]["i_option"] = 10

# Field line tracing: cheap, no poincare/divertor/collisions
field_line_tracing["field_line_tracing_nml"]["time_step"] = 1.0e-4
field_line_tracing["field_line_tracing_nml"]["energy_ev"] = 3.5e3
field_line_tracing["field_line_tracing_nml"]["n_particles"] = 3
field_line_tracing["field_line_tracing_nml"]["boole_poincare_plot"] = False
field_line_tracing["field_line_tracing_nml"]["boole_divertor_intersection"] = False
field_line_tracing["field_line_tracing_nml"]["boole_point_source"] = False
field_line_tracing["field_line_tracing_nml"]["boole_refined_sqrt_g"] = True
field_line_tracing["field_line_tracing_nml"]["boole_monoenergetic"] = True
field_line_tracing["field_line_tracing_nml"]["boole_linear_density_simulation"] = False
field_line_tracing["field_line_tracing_nml"]["boole_linear_temperature_simulation"] = False
field_line_tracing["field_line_tracing_nml"]["seed_option"] = 2

# Deterministic seed for reproducible CI output
(WORK_DIR / "seed.inp").write_text("8\n  1 2 3 4 5 6 7 8\n")

# Axisymmetric equilibrium only (ipert = 0 in the upstream blueprint); gfile,
# convexfile and iaxieq are overridden from tetra_grid_settings_mod, so their
# values are read into a dummy and ignored.
shutil.copy(GORILLA_ROOT / "INPUT" / "field_divB0.inp", WORK_DIR / "field_divB0.inp")

# Symlink the bundled ASDEX g-file + convex wall into the work dir
(WORK_DIR / "DATA").mkdir()
(WORK_DIR / "DATA" / "ASDEX").symlink_to(ASDEX_DIR, target_is_directory=True)

# Write namelists
gorilla.write(str(WORK_DIR / "gorilla.inp"), force=True)
tetra_grid.write(str(WORK_DIR / "tetra_grid.inp"), force=True)
gorilla_applets.write(str(WORK_DIR / "gorilla_applets.inp"), force=True)
field_line_tracing.write(str(WORK_DIR / "field_line_tracing.inp"), force=True)

print("=== field line tracing: i_option=10 ===", flush=True)
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
expected_marker = "Number of lost particles"
if expected_marker not in result.stdout:
    sys.exit(f"FAIL: expected '{expected_marker}' in stdout but it was missing")

print(f"PASS: field line tracing completed and summary marker found")
