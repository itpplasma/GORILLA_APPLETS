#!/usr/bin/env python3
"""
CTest driver for reference tests of the applets without other end-to-end coverage.

Usage: test_applet_reference.py <mode>, with <mode> one of the keys of MODES.

Uses the mode's input file from EXAMPLES/ as blueprint, so the shipped example
is checked to parse and run, and overrides only what is needed for a short CI
run: the shared EFIT test equilibrium without perturbation field (ipert = 0),
a small grid, few particles and a short tracing time. The outputs are
compared against TESTS/REFERENCE/<mode>/. The run is pinned to one OpenMP
thread. Each marker has its own random stream, so the marker histories do not
depend on the thread count, but sums over markers can differ by rounding. After an intended change of the results, run the test with
UPDATE_REFERENCE=1 in the environment to copy the files listed in the mode's
"reference_files" from the work dir into the reference dir, and commit them.

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


# Per mode: example dir, i_option, mode input file and namelist, grid and
# mode-namelist overrides, and the output files compared against the reference.
MODES = {
    "boltzmann": {
        "example": "boltzmann_test",
        "i_option": 9,
        "inp": "boltzmann.inp",
        "nml": "boltzmann_nml",
        "grid": {"grid_kind": 2, "n1": 20, "n2": 30, "n3": 50},
        "overrides": {
            "time_step": 1.0e-4,
            "n_particles": 20,
        },
        "reference_files": ["exit_data.dat", "fourier_moments_summary.dat"],
    },
    "divertor_heat_loads": {
        "example": "divertor_heat_loads_lambda_0.7",
        "i_option": 11,
        "inp": "divertor_heat_loads.inp",
        "nml": "divertor_heat_loads_nml",
        "grid": {"grid_kind": 1, "n1": 50, "n2": 30, "n3": 50},
        "overrides": {
            "time_step": 1.0e-4,
            "n_particles": 10,
            "n_poincare_mappings": 20,
            "n_mappings_ignored": 2,
        },
        "reference_files": ["exit_data.dat", "poincare_maps.dat", "divertor_intersections.dat",
                            "fourier_moments_summary.dat"],
    },
    "anomalous_transport": {
        "example": "anomalous_transport",
        "i_option": 13,
        "inp": "anomalous_transport.inp",
        "nml": "anomalous_transport_nml",
        "grid": {"grid_kind": 2, "n1": 20, "n2": 30, "n3": 50},
        "overrides": {
            "time_step": 5.0e-5,   # the anomalous kick dominates the runtime
            "n_particles": 10,
            "i_scan_option": 0,                        # the example's scan takes minutes
            "boole_calc_diffusion_coefficient": False,
            "anomalous_diffusion_coefficient": 1.0e2,  # exercise the kick (1e4 takes minutes)
        },
        "reference_files": ["exit_data.dat", "fourier_moments_summary.dat"],
    },
}


def env_path(name: str) -> Path:
    value = os.environ.get(name)
    if not value:
        sys.exit(f"missing required env var: {name}")
    return Path(value)


if len(sys.argv) != 2 or sys.argv[1] not in MODES:
    sys.exit(f"usage: {Path(__file__).name} <{'|'.join(MODES)}>")
MODE_NAME = sys.argv[1]
MODE = MODES[MODE_NAME]

APPLETS_ROOT = env_path("APPLETS_ROOT")
GORILLA_ROOT = env_path("GORILLA_ROOT")
BINARY = env_path("GORILLA_APPLETS_BIN")
WORK_DIR = env_path("WORK_DIR")
EXAMPLE_DIR = APPLETS_ROOT / "EXAMPLES" / MODE["example"]
REFERENCE_DIR = APPLETS_ROOT / "TESTS" / "REFERENCE" / MODE_NAME

# Fresh work dir for every test run
if WORK_DIR.exists():
    shutil.rmtree(WORK_DIR)
WORK_DIR.mkdir(parents=True)

# Load blueprints: core inputs from GORILLA, the dispatcher and the mode input
# from the shipped example
gorilla = f90nml.read(str(GORILLA_ROOT / "INPUT" / "gorilla.inp"))
tetra_grid = f90nml.read(str(GORILLA_ROOT / "INPUT" / "tetra_grid.inp"))
gorilla_applets = f90nml.read(str(EXAMPLE_DIR / "gorilla_applets.inp"))
mode_inp = f90nml.read(str(EXAMPLE_DIR / MODE["inp"]))
for nml in (gorilla, tetra_grid, gorilla_applets, mode_inp):
    nml.end_comma = True

if gorilla_applets["gorilla_applets_nml"]["i_option"] != MODE["i_option"]:
    sys.exit(f"FAIL: the {MODE['example']} example must select i_option = {MODE['i_option']}")

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
gorilla["gorillanml"]["boole_time_hamiltonian"] = True  # deposit the time-integral moment

# Tetrahedral grid: shared EFIT, small enough for CI
for key, value in MODE["grid"].items():
    tetra_grid["tetra_grid_nml"][key] = value
tetra_grid["tetra_grid_nml"]["boole_n_field_periods"] = True
tetra_grid["tetra_grid_nml"]["g_file_filename"] = "MHD_EQUILIBRIA/g_file_for_test"
tetra_grid["tetra_grid_nml"]["theta_geom_flux"] = 2  # cylindrical coordinates: theta_geom_flux = 1 folds near the separatrix
tetra_grid["tetra_grid_nml"]["convex_wall_filename"] = "MHD_EQUILIBRIA/convex_wall_for_test.dat"

# Fourier moments instead of the much larger per-prism moments
mode_inp[MODE["nml"]]["boole_write_moments"] = False
mode_inp[MODE["nml"]]["boole_write_fourier_moments"] = True
for key, value in MODE["overrides"].items():
    mode_inp[MODE["nml"]][key] = value

# Deterministic seed for reproducible CI output
(WORK_DIR / "seed.inp").write_text("8\n  1 2 3 4 5 6 7 8\n")

# Axisymmetric equilibrium only (ipert = 0 in the upstream blueprint);
# the field-aligned grid additionally needs the SYNCH preload namelist.
shutil.copy(GORILLA_ROOT / "INPUT" / "field_divB0.inp", WORK_DIR / "field_divB0.inp")
shutil.copy(GORILLA_ROOT / "INPUT" / "preload_for_SYNCH.inp",
            WORK_DIR / "preload_for_SYNCH.inp")
(WORK_DIR / "MHD_EQUILIBRIA").symlink_to(
    GORILLA_ROOT / "MHD_EQUILIBRIA", target_is_directory=True)

gorilla.write(str(WORK_DIR / "gorilla.inp"), force=True)
tetra_grid.write(str(WORK_DIR / "tetra_grid.inp"), force=True)
gorilla_applets.write(str(WORK_DIR / "gorilla_applets.inp"), force=True)
mode_inp.write(str(WORK_DIR / MODE["inp"]), force=True)

timer = StageTimer(MODE_NAME, output_path=WORK_DIR / "timings.json")

print(f"=== {MODE_NAME}: i_option={MODE['i_option']} ===", flush=True)
with timer.time(MODE_NAME):
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


def check_rows(name: str) -> None:
    path = WORK_DIR / name
    if not path.exists():
        sys.exit(f"FAIL: output file missing: {name}")
    rows = [ln.split() for ln in path.read_text().splitlines() if ln.strip()]
    if not rows:
        sys.exit(f"FAIL: output file empty: {name}")
    values = [float(t) for row in rows for t in row]
    if not all(math.isfinite(v) for v in values):
        sys.exit(f"FAIL: non-finite value in {name}")


def summarise_moments(name: str, summary_name: str) -> None:
    """Reduce a moments file (one row per triangle) to per-column sums, sums of
    absolute values, sums of squares and row-index weighted sums, which keeps
    the committed reference small but still sensitive to every deposited value."""
    rows = [[float(t) for t in ln.split()] for ln in (WORK_DIR / name).read_text().splitlines() if ln.strip()]
    columns = list(zip(*rows))
    stats = [
        [sum(c) for c in columns],
        [sum(abs(v) for v in c) for c in columns],
        [sum(v * v for v in c) for c in columns],
        [sum(i * v for i, v in enumerate(c, start=1)) for c in columns],
    ]
    (WORK_DIR / summary_name).write_text(
        "\n".join(" ".join(f"{v:.16e}" for v in row) for row in stats) + "\n")


if "fourier_moments_summary.dat" in MODE["reference_files"]:
    check_rows("fourier_moments.dat")
    summarise_moments("fourier_moments.dat", "fourier_moments_summary.dat")

if os.environ.get("UPDATE_REFERENCE") == "1":
    REFERENCE_DIR.mkdir(parents=True, exist_ok=True)
    for name in MODE["reference_files"]:
        check_rows(name)
        shutil.copy(WORK_DIR / name, REFERENCE_DIR / name)
    sys.exit(f"Reference files of {MODE_NAME} updated in {REFERENCE_DIR}; rerun without UPDATE_REFERENCE")

def compare_moments_summary(name: str, n_rows: int, rtol: float) -> None:
    """Compare a summary written by summarise_moments. The sums of absolute values and of
    squares are compared relative to themselves. The signed sums (row 1) and row-index
    weighted sums (row 4) can cancel to much less than their terms, so they are compared
    relative to the size of their terms, sum|v| and n_rows*sum|v| of the reference."""
    obs = [[float(t) for t in ln.split()] for ln in (WORK_DIR / name).read_text().splitlines() if ln.strip()]
    ref = [[float(t) for t in ln.split()] for ln in (REFERENCE_DIR / name).read_text().splitlines() if ln.strip()]
    if len(obs) != 4 or len(ref) != 4 or any(len(o) != len(r) for o, r in zip(obs, ref)):
        sys.exit(f"FAIL: shape mismatch in {name}")
    scale = [ref[1], ref[1], ref[2], [n_rows * v for v in ref[1]]]
    for i in range(4):
        for j, (a, b) in enumerate(zip(obs[i], ref[i])):
            tol = rtol * abs(scale[i][j]) + 1.0e-14
            if abs(a - b) > tol:
                sys.exit(f"FAIL: value mismatch in {name} line {i + 1} token {j + 1}: observed {a!r} "
                         f"vs reference {b!r} (|diff|={abs(a - b):.3e}, tol={tol:.3e})")
    print(f"PASS: {name} matches reference within tolerance.")


for name in MODE["reference_files"]:
    check_rows(name)   # compare_numeric_file would let NaN through
    if name == "fourier_moments_summary.dat":
        n_rows = sum(1 for ln in (WORK_DIR / "fourier_moments.dat").read_text().splitlines() if ln.strip())
        compare_moments_summary(name, n_rows, rtol=1.0e-10)
    else:
        compare_numeric_file(WORK_DIR / name, REFERENCE_DIR / name, rtol=1.0e-10, atol=1.0e-14)

print(f"PASS: {MODE_NAME} outputs match the reference")
