#!/usr/bin/env python3
"""
CTest driver for the alpha-particle lifetime CI test.

Reads GORILLA-core blueprints from ../GORILLA/INPUT/ and applet blueprints
from GORILLA_APPLETS/INPUT/, sets every key the test depends on explicitly,
then invokes the binary in two stages:

  1) i_option = 1  -> precompute fluxtubevolume.dat for the test grid
  2) i_option = 5  -> trace alpha particles and write
                      alpha_lifetime_gorilla.dat

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

MHD_DIR = APPLETS_ROOT / "EXAMPLES" / "ALPHA_LIFETIME" / "MHD_EQUILIBRIA"

# Fresh work dir for every test run
if WORK_DIR.exists():
    shutil.rmtree(WORK_DIR)
WORK_DIR.mkdir(parents=True)

# Load blueprints from their conceptual owners
gorilla = f90nml.read(str(GORILLA_ROOT / "INPUT" / "gorilla.inp"))
tetra_grid = f90nml.read(str(GORILLA_ROOT / "INPUT" / "tetra_grid.inp"))
gorilla_applets = f90nml.read(str(APPLETS_ROOT / "INPUT" / "gorilla_applets.inp"))
alpha_lifetime = f90nml.read(str(APPLETS_ROOT / "INPUT" / "alpha_lifetime.inp"))
for nml in (gorilla, tetra_grid, gorilla_applets, alpha_lifetime):
    nml.end_comma = True

# Core integrator
gorilla["gorillanml"]["eps_phi"] = 0.0
gorilla["gorillanml"]["coord_system"] = 2          # (s, theta, phi)
gorilla["gorillanml"]["ispecies"] = 3              # alpha particle
gorilla["gorillanml"]["boole_periodic_relocation"] = True
gorilla["gorillanml"]["ipusher"] = 2               # polynomial pusher
gorilla["gorillanml"]["poly_order"] = 2            # cheap order for CI
gorilla["gorillanml"]["boole_grid_for_find_tetra"] = False
gorilla["gorillanml"]["boole_adaptive_time_steps"] = False

# Tetrahedral grid (small for CI)
tetra_grid["tetra_grid_nml"]["grid_kind"] = 3      # VMEC field-aligned
tetra_grid["tetra_grid_nml"]["n1"] = 10            # ns
tetra_grid["tetra_grid_nml"]["n2"] = 10            # nphi
tetra_grid["tetra_grid_nml"]["n3"] = 10            # ntheta
tetra_grid["tetra_grid_nml"]["boole_n_field_periods"] = True
tetra_grid["tetra_grid_nml"]["sfc_s_min"] = 1.0e-1
tetra_grid["tetra_grid_nml"]["i_radial_spacing"] = 1
tetra_grid["tetra_grid_nml"]["theta_geom_flux"] = 1
tetra_grid["tetra_grid_nml"]["theta0_at_xpoint"] = True
tetra_grid["tetra_grid_nml"]["g_file_filename"] = "MHD_EQUILIBRIA/g_file_for_test"
tetra_grid["tetra_grid_nml"]["convex_wall_filename"] = "MHD_EQUILIBRIA/convex_wall_for_test.dat"
tetra_grid["tetra_grid_nml"]["netcdf_filename"] = "MHD_EQUILIBRIA/netcdf_file_for_test.nc"

# Applets dispatcher / fluxtube precomputation
gorilla_applets["gorilla_applets_nml"]["filename_fluxtv_precomp"] = "fluxtubevolume.dat"
gorilla_applets["gorilla_applets_nml"]["filename_fluxtv_load"] = "fluxtubevolume.dat"
gorilla_applets["gorilla_applets_nml"]["start_pos_x1"] = 0.6
gorilla_applets["gorilla_applets_nml"]["start_pos_x2"] = 0.00013
gorilla_applets["gorilla_applets_nml"]["start_pos_x3"] = 0.00013
gorilla_applets["gorilla_applets_nml"]["t_step_fluxtv"] = 1.0e-3
gorilla_applets["gorilla_applets_nml"]["nt_steps_fluxtv"] = 200
gorilla_applets["gorilla_applets_nml"]["energy_ev_fluxtv"] = 3.0e3

# Alpha lifetime tracing
alpha_lifetime["alpha_lifetimenml"]["time_step"] = 1.0e-5
alpha_lifetime["alpha_lifetimenml"]["energy_ev"] = 3.5e6
alpha_lifetime["alpha_lifetimenml"]["n_particles"] = 10
alpha_lifetime["alpha_lifetimenml"]["i_integrator_type"] = 1
alpha_lifetime["alpha_lifetimenml"]["seed_option"] = 2
alpha_lifetime["alpha_lifetimenml"]["boole_random_precalc"] = True
alpha_lifetime["alpha_lifetimenml"]["filename_alpha_lifetime"] = "alpha_lifetime_gorilla.dat"

# Deterministic seed for reproducible CI output
(WORK_DIR / "seed.inp").write_text("8\n  1 2 3 4 5 6 7 8\n")

# Symlink MHD equilibrium data into the work dir
(WORK_DIR / "MHD_EQUILIBRIA").symlink_to(MHD_DIR, target_is_directory=True)

# Write namelists (gorilla_applets.inp is rewritten between stages)
gorilla.write(str(WORK_DIR / "gorilla.inp"), force=True)
tetra_grid.write(str(WORK_DIR / "tetra_grid.inp"), force=True)
alpha_lifetime.write(str(WORK_DIR / "alpha_lifetime.inp"), force=True)


def run_stage(i_option: int, label: str) -> None:
    gorilla_applets["gorilla_applets_nml"]["i_option"] = i_option
    gorilla_applets.write(str(WORK_DIR / "gorilla_applets.inp"), force=True)
    print(f"=== {label}: i_option={i_option} ===", flush=True)
    subprocess.run([str(BINARY)], cwd=WORK_DIR, check=True)


run_stage(1, "fluxtube precomputation")
flux_file = WORK_DIR / "fluxtubevolume.dat"
if not flux_file.exists() or flux_file.stat().st_size == 0:
    sys.exit(f"FAIL: {flux_file} not produced by stage 1")

run_stage(5, "alpha lifetime tracing")
out_file = WORK_DIR / "alpha_lifetime_gorilla.dat"
if not out_file.exists() or out_file.stat().st_size == 0:
    sys.exit(f"FAIL: {out_file} not produced by stage 2")

with out_file.open() as fh:
    n_lines = sum(1 for line in fh if line.strip())
expected = int(alpha_lifetime["alpha_lifetimenml"]["n_particles"])
if n_lines != expected:
    sys.exit(f"FAIL: expected {expected} confinement times, got {n_lines}")

print(f"PASS: {n_lines} confinement times written to {out_file.name}")
