#!/usr/bin/env python3
"""
CTest driver for the energy and momentum conserving collisions of helical core.

A cold deuterium beam (pitch +1, density 0.05 of the background) is pushed in four
time slices through a field-aligned grid of the shared EFIT test equilibrium, with
the background reservoir taking up the momentum and energy the markers lose in
collisions. The background exchange is summed in thread order, so its rounding
depends on the number of threads; the test therefore checks physics properties
of collision_relaxation.dat instead of reference files:

- every slice writes a row, no background update is rejected (T <= 0) and no marker
  is reflected across the axis hole;
- total energy of markers, background and lost markers is conserved;
- the beam slows down and the background flow grows in the beam direction.

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
N_SLICES = 4

if WORK_DIR.exists():
    shutil.rmtree(WORK_DIR)
WORK_DIR.mkdir(parents=True)

gorilla = f90nml.read(str(GORILLA_ROOT / "INPUT" / "gorilla.inp"))
tetra_grid = f90nml.read(str(GORILLA_ROOT / "INPUT" / "tetra_grid.inp"))
gorilla_applets = f90nml.read(str(APPLETS_ROOT / "INPUT" / "gorilla_applets.inp"))
helical_core = f90nml.read(str(APPLETS_ROOT / "INPUT" / "helical_core.inp"))
for nml in (gorilla, tetra_grid, gorilla_applets, helical_core):
    nml.end_comma = True

g = gorilla["gorillanml"]
g["eps_phi"] = 0.0
g["coord_system"] = 1            # (R, phi, Z)
g["ispecies"] = 2                # deuterium
g["boole_periodic_relocation"] = True
g["ipusher"] = 2                 # polynomial pusher
g["poly_order"] = 2
g["boole_guess"] = True
g["boole_grid_for_find_tetra"] = True
g["boole_adaptive_time_steps"] = False

# Coarse field-aligned grid with a tiny axis hole (no reflections across the hole,
# which would change the marker energy); the small tetrahedra near the axis are
# merged into reservoir cells (reservoir_min_volume_fraction, default 0.5)
tg = tetra_grid["tetra_grid_nml"]
tg["grid_kind"] = 2
tg["n1"] = 10
tg["n2"] = 10
tg["n3"] = 20
tg["i_radial_spacing"] = 2
tg["sfc_s_min"] = 1.0e-7
tg["boole_n_field_periods"] = True
tg["g_file_filename"] = "MHD_EQUILIBRIA/g_file_for_test"
tg["theta_geom_flux"] = 2  # cylindrical coordinates: theta_geom_flux = 1 folds near the separatrix
tg["convex_wall_filename"] = "MHD_EQUILIBRIA/convex_wall_for_test.dat"

gorilla_applets["gorilla_applets_nml"]["i_option"] = 14

h = helical_core["helical_core_nml"]
h["time_step"] = 2.0e-5
h["energy_ev"] = 1.0e3
h["n_species"] = 1
h["n_particles"] = 2000
h["density"] = 3.0e14
h["boole_squared_moments"] = False
h["boole_point_source"] = False
h["boole_collisions"] = True
h["boole_preserve_energy_and_momentum_during_collisions"] = True
h["boole_refined_sqrt_g"] = True
h["boole_monoenergetic"] = True
h["boole_linear_density_simulation"] = False
h["boole_linear_temperature_simulation"] = False
h["boole_eliminate_particles_outside_flux"] = True
h["flux_threshold_for_elimination"] = 0.7
h["boole_delta_f"] = False
h["boole_antithetic_variate"] = False
h["boole_unidirectional_beam"] = True
h["beam_density_fraction"] = 0.05
h["n_time_slices"] = N_SLICES
h["i_integrator_type"] = 1
h["seed_option"] = 2
for key in ("boole_write_vertex_indices", "boole_write_vertex_coordinates", "boole_write_prism_volumes",
            "boole_write_refined_prism_volumes", "boole_write_moments", "boole_write_fourier_moments",
            "boole_write_exit_data", "boole_write_grid_data"):
    h[key] = False

(WORK_DIR / "seed.inp").write_text("8\n  1 2 3 4 5 6 7 8\n")
shutil.copy(GORILLA_ROOT / "INPUT" / "field_divB0.inp", WORK_DIR / "field_divB0.inp")
shutil.copy(GORILLA_ROOT / "INPUT" / "preload_for_SYNCH.inp", WORK_DIR / "preload_for_SYNCH.inp")
(WORK_DIR / "MHD_EQUILIBRIA").symlink_to(GORILLA_ROOT / "MHD_EQUILIBRIA", target_is_directory=True)

gorilla.write(str(WORK_DIR / "gorilla.inp"), force=True)
tetra_grid.write(str(WORK_DIR / "tetra_grid.inp"), force=True)
gorilla_applets.write(str(WORK_DIR / "gorilla_applets.inp"), force=True)
helical_core.write(str(WORK_DIR / "helical_core.inp"), force=True)

print("=== helical core, conserving collisions, beam relaxation ===", flush=True)
result = subprocess.run([str(BINARY)], cwd=WORK_DIR, capture_output=True, text=True, check=False, timeout=300)
print(result.stdout[-3000:])
if result.stderr:
    print(result.stderr, file=sys.stderr)
if result.returncode != 0:
    sys.exit(f"FAIL: gorilla_applets_main.x exited with code {result.returncode}")
if "background updates rejected" in result.stdout:
    sys.exit("FAIL: background updates were rejected (T <= 0), energy and momentum not conserved")
if "reservoir cells" not in result.stdout:
    sys.exit("FAIL: no reservoir cells were built near the axis")
hole_marker = "pushed across the inside hole ="
if hole_marker in result.stdout and int(result.stdout.split(hole_marker)[1].split()[0]) != 0:
    sys.exit("FAIL: markers were reflected across the axis hole, which does not conserve energy")

# Columns of collision_relaxation.dat (see write_collision_relaxation_diagnostics)
rows = []
for line in (WORK_DIR / "collision_relaxation.dat").read_text().splitlines():
    if line.strip() and not line.startswith("#"):
        rows.append([float(v) for v in line.split()])
if len(rows) != N_SLICES + 1:
    sys.exit(f"FAIL: expected {N_SLICES + 1} rows in collision_relaxation.dat, found {len(rows)}")
vpar_m = [r[1] for r in rows]
u_b = [r[2] for r in rows]
e_tot = [r[9] for r in rows]

rel_energy_error = abs(e_tot[-1] / e_tot[0] - 1.0)
print(f"<v_par>_m/v0: {vpar_m}")
print(f"<u_par>_b/v0: {u_b}")
print(f"relative energy error: {rel_energy_error:.2e}")
# rejected updates give about 1e-5, clean runs at most a few 1e-8
if rel_energy_error > 1.0e-6:
    sys.exit(f"FAIL: total energy not conserved, relative error {rel_energy_error:.2e}")
if not all(later < earlier for earlier, later in zip(vpar_m, vpar_m[1:])):
    sys.exit("FAIL: the beam does not slow down monotonically")
if not all(later > earlier for earlier, later in zip(u_b, u_b[1:])):
    sys.exit("FAIL: the background flow does not grow monotonically in the beam direction")

print("PASS: beam slows down, background spins up, energy conserved")
