#!/usr/bin/env python3
"""
CTest driver for the self-consistent electric field (SCEF) CI smoke test.

Uses the input files of EXAMPLES/self_consistent_electric_field/ as blueprints,
so the shipped example is checked to parse and run, and overrides only what is
needed for a short CI run: the shared VMEC test equilibrium instead of the
W7-X file the example points to, a small grid, few particles, a fresh
diffusion calibration for both species (no cached A_and_B_species<k>.dat) and
one source update after two electric potential updates. The fixed cost is
splining the VMEC equilibrium (about 6 s).

The outputs are compared against TESTS/REFERENCE/self_consistent_electric_field/.
The run is pinned to one OpenMP thread, because the random numbers a particle
draws depend on the thread that traces it. After an intended change of the
results, rerun the test and copy the files listed in REFERENCE_FILES from the
work dir into the reference dir.

Paths are passed in by TESTS/CMakeLists.txt via environment variables.
"""
import math
import os
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _ci_utils import StageTimer, compare_numeric_file

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
EXAMPLE_DIR = APPLETS_ROOT / "EXAMPLES" / "self_consistent_electric_field"
REFERENCE_DIR = APPLETS_ROOT / "TESTS" / "REFERENCE" / "self_consistent_electric_field"
REFERENCE_FILES = [
    "A_and_B_species1.dat",     # diffusion calibration
    "A_and_B_species2.dat",
    "phi_elec_1d_1_1.dat",      # electric potential updates
    "phi_elec_1d_1_2.dat",
    "one_d_densities1_1.dat",   # flux-surface densities
    "one_d_densities1_2.dat",
    "particle_source_1.dat",    # source before and after the source update
    "particle_source_2.dat",
]

# Fresh work dir for every test run
if WORK_DIR.exists():
    shutil.rmtree(WORK_DIR)
WORK_DIR.mkdir(parents=True)

# Load blueprints from the shipped example
gorilla = f90nml.read(str(EXAMPLE_DIR / "gorilla.inp"))
tetra_grid = f90nml.read(str(EXAMPLE_DIR / "tetra_grid.inp"))
gorilla_applets = f90nml.read(str(EXAMPLE_DIR / "gorilla_applets.inp"))
scef = f90nml.read(str(EXAMPLE_DIR / "self_consistent_ef.inp"))
for nml in (gorilla, tetra_grid, gorilla_applets, scef):
    nml.end_comma = True

if gorilla_applets["gorilla_applets_nml"]["i_option"] != 12:
    sys.exit("FAIL: the SCEF example must select i_option = 12")

# Tetrahedral grid: shared VMEC test equilibrium, small enough for CI
tetra_grid["tetra_grid_nml"]["netcdf_filename"] = "MHD_EQUILIBRIA/netcdf_file_for_test.nc"
tetra_grid["tetra_grid_nml"]["n1"] = 12
tetra_grid["tetra_grid_nml"]["n2"] = 10
tetra_grid["tetra_grid_nml"]["n3"] = 10

# SCEF: random walk for both species with a fresh, small diffusion calibration
nml = scef["self_consistent_ef_nml"]
nml["n_particles"] = 100
nml["n_particles_diffusion"] = 50
nml["diffusion_tau"] = 5.0e-6   # default 3.4e-2 s (ions) dominates the runtime
nml["boole_honest_tracing"] = [False, False]
nml["boole_recompute_D"] = [True, True]
nml["n_source_updates"] = 1
nml["n_electric_potential_updates"] = 2
nml["update_dimension"] = 1
nml["seed_option"] = 2

# Deterministic seed for reproducible CI output
(WORK_DIR / "seed.inp").write_text("8\n  1 2 3 4 5 6 7 8\n")

shutil.copy(GORILLA_ROOT / "INPUT" / "field_divB0.inp", WORK_DIR / "field_divB0.inp")
shutil.copy(GORILLA_ROOT / "INPUT" / "preload_for_SYNCH.inp",
            WORK_DIR / "preload_for_SYNCH.inp")
(WORK_DIR / "MHD_EQUILIBRIA").symlink_to(
    GORILLA_ROOT / "MHD_EQUILIBRIA", target_is_directory=True)

gorilla.write(str(WORK_DIR / "gorilla.inp"), force=True)
tetra_grid.write(str(WORK_DIR / "tetra_grid.inp"), force=True)
gorilla_applets.write(str(WORK_DIR / "gorilla_applets.inp"), force=True)
scef.write(str(WORK_DIR / "self_consistent_ef.inp"), force=True)

timer = StageTimer("self_consistent_electric_field", output_path=WORK_DIR / "timings.json")

print("=== self-consistent electric field: i_option=12 ===", flush=True)
with timer.time("scef"):
    result = subprocess.run(
        [str(BINARY)],
        cwd=WORK_DIR,
        env={**os.environ, "OMP_NUM_THREADS": "1"},
        capture_output=True,
        text=True,
        check=False,
        timeout=300,
    )
print(result.stdout)
if result.stderr:
    print(result.stderr, file=sys.stderr)
if result.returncode != 0:
    sys.exit(f"FAIL: gorilla_applets_main.x exited with code {result.returncode}")
# Input errors end in a plain Fortran stop, which exits with code 0
if "Program terminated" in result.stdout:
    sys.exit("FAIL: input check terminated the program")


def read_rows(name: str):
    path = WORK_DIR / name
    if not path.exists():
        sys.exit(f"FAIL: output file missing: {name}")
    rows = [[float(t) for t in ln.split()] for ln in path.read_text().splitlines() if ln.strip()]
    if not rows:
        sys.exit(f"FAIL: output file empty: {name}")
    if not all(math.isfinite(v) for row in rows for v in row):
        sys.exit(f"FAIL: non-finite value in {name}")
    return rows


for name in REFERENCE_FILES:
    read_rows(name)   # compare_numeric_file would let NaN through
    compare_numeric_file(WORK_DIR / name, REFERENCE_DIR / name, rtol=1.0e-10, atol=1.0e-14)

print("PASS: SCEF calibration, potential and source updates match the reference")
