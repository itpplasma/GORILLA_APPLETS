#!/usr/bin/env python3
"""
CTest driver for the mono-energetic transport coefficient CI test.

Reads GORILLA-core blueprints from ../GORILLA/INPUT/ and applet blueprints
from GORILLA_APPLETS/INPUT/, sets every key the test depends on explicitly,
then invokes the binary in two stages:

  1) i_option = 1  -> precompute fluxtubevolume.dat for the test grid
  2) i_option = 2  -> single nu*: compute D11 and write
                      nustar_diffcoef_std.dat

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
COMPARE_BACKEND = os.environ.get("TRANSPORT_COMPARE_BACKEND", "cpu-batch")
if COMPARE_BACKEND not in ("cpu-batch", "cuda"):
    sys.exit("TRANSPORT_COMPARE_BACKEND must be cpu-batch or cuda")

MHD_DIR = GORILLA_ROOT / "MHD_EQUILIBRIA"
REFERENCE_DIR = APPLETS_ROOT / "TESTS" / "REFERENCE" / "mono_energetic_transport"

# Fresh work dir for every test run
if WORK_DIR.exists():
    shutil.rmtree(WORK_DIR)
WORK_DIR.mkdir(parents=True)

# Load blueprints from their conceptual owners
gorilla = f90nml.read(str(GORILLA_ROOT / "INPUT" / "gorilla.inp"))
tetra_grid = f90nml.read(str(GORILLA_ROOT / "INPUT" / "tetra_grid.inp"))
gorilla_applets = f90nml.read(str(APPLETS_ROOT / "INPUT" / "gorilla_applets.inp"))
mono = f90nml.read(str(APPLETS_ROOT / "INPUT" / "mono_energetic_transp_coef.inp"))
for nml in (gorilla, tetra_grid, gorilla_applets, mono):
    nml.end_comma = True

# Core integrator (flux coordinates, polynomial pusher, cheap order for CI)
gorilla["gorillanml"]["eps_phi"] = 0.0
gorilla["gorillanml"]["coord_system"] = 2          # (s, theta, phi)
gorilla["gorillanml"]["ispecies"] = 1              # electron
gorilla["gorillanml"]["boole_periodic_relocation"] = True
gorilla["gorillanml"]["ipusher"] = 2               # polynomial pusher
gorilla["gorillanml"]["poly_order"] = 2
gorilla["gorillanml"]["i_precomp"] = 0
gorilla["gorillanml"]["i_time_tracing_option"] = 1
gorilla["gorillanml"]["handover_processing_kind"] = 1
gorilla["gorillanml"]["boole_strong_electric_field"] = False
gorilla["gorillanml"]["boole_grid_for_find_tetra"] = False
gorilla["gorillanml"]["boole_adaptive_time_steps"] = False

# Tetrahedral grid: VMEC field-aligned, small for CI
tetra_grid["tetra_grid_nml"]["grid_kind"] = 3
tetra_grid["tetra_grid_nml"]["n1"] = 10            # ns
tetra_grid["tetra_grid_nml"]["n2"] = 10            # nphi
tetra_grid["tetra_grid_nml"]["n3"] = 10            # ntheta
tetra_grid["tetra_grid_nml"]["boole_n_field_periods"] = True
tetra_grid["tetra_grid_nml"]["sfc_s_min"] = 1.0e-1
tetra_grid["tetra_grid_nml"]["theta_geom_flux"] = 1
tetra_grid["tetra_grid_nml"]["theta0_at_xpoint"] = True
tetra_grid["tetra_grid_nml"]["g_file_filename"] = "MHD_EQUILIBRIA/g_file_for_test"
tetra_grid["tetra_grid_nml"]["convex_wall_filename"] = "MHD_EQUILIBRIA/convex_wall_for_test.dat"
tetra_grid["tetra_grid_nml"]["netcdf_filename"] = "MHD_EQUILIBRIA/netcdf_file_for_test.nc"

# Applets dispatcher / fluxtube precomputation
gorilla_applets["gorilla_applets_nml"]["filename_fluxtv_precomp"] = "fluxtubevolume.dat"
gorilla_applets["gorilla_applets_nml"]["filename_fluxtv_load"] = "fluxtubevolume.dat"
gorilla_applets["gorilla_applets_nml"]["start_pos_x1"] = 0.5
gorilla_applets["gorilla_applets_nml"]["start_pos_x2"] = 0.00013
gorilla_applets["gorilla_applets_nml"]["start_pos_x3"] = 0.00013
gorilla_applets["gorilla_applets_nml"]["t_step_fluxtv"] = 1.0e-3
gorilla_applets["gorilla_applets_nml"]["nt_steps_fluxtv"] = 200
gorilla_applets["gorilla_applets_nml"]["energy_ev_fluxtv"] = 3.0e3

# Mono-energetic transport: single nu* (i_option=2), no scan
mono["transpcoefnml"]["i_integrator_type"] = 1
mono["transpcoefnml"]["transport_backend"] = "cpu"
mono["transpcoefnml"]["n_particles"] = 3
mono["transpcoefnml"]["boole_collisions"] = True
mono["transpcoefnml"]["energy_ev"] = 3.0e3
mono["transpcoefnml"]["boole_random_precalc"] = True
mono["transpcoefnml"]["seed_option"] = 2
mono["transpcoefnml"]["random_seed_filename"] = "seed.inp"
mono["transpcoefnml"]["idiffcoef_output"] = 1
mono["transpcoefnml"]["filename_transp_diff_coef"] = "nustar_diffcoef_std.dat"
mono["transpcoefnml"]["nu_star"] = 0.1
mono["transpcoefnml"]["v_e"] = 0.0
# One transport timescale is enough for this deterministic output/interface test.
mono["transpcoefnml"]["flight_time_multiplier"] = 1.0
mono["transpcoefnml"]["boole_psi_mat"] = True
mono["transpcoefnml"]["boole_write_particle_histories"] = True
mono["transpcoefnml"]["filename_particle_histories"] = "particle_histories.dat"
mono["transpcoefnml"]["filename_transport_metadata"] = "transport_metadata.dat"

# Deterministic seed for reproducible CI output
(WORK_DIR / "seed.inp").write_text("8\n  1 2 3 4 5 6 7 8\n")

# Reuse the field_divB0.inp blueprint
shutil.copy(GORILLA_ROOT / "INPUT" / "field_divB0.inp", WORK_DIR / "field_divB0.inp")

# Symlink MHD equilibrium data into the work dir
(WORK_DIR / "MHD_EQUILIBRIA").symlink_to(MHD_DIR, target_is_directory=True)

# Write namelists (gorilla_applets.inp rewritten between stages)
gorilla.write(str(WORK_DIR / "gorilla.inp"), force=True)
tetra_grid.write(str(WORK_DIR / "tetra_grid.inp"), force=True)
mono.write(str(WORK_DIR / "mono_energetic_transp_coef.inp"), force=True)


timer = StageTimer("mono_energetic_transport", output_path=WORK_DIR / "timings.json")


def run_stage(i_option: int, label: str) -> None:
    gorilla_applets["gorilla_applets_nml"]["i_option"] = i_option
    gorilla_applets.write(str(WORK_DIR / "gorilla_applets.inp"), force=True)
    mono.write(str(WORK_DIR / "mono_energetic_transp_coef.inp"), force=True)
    print(f"=== {label}: i_option={i_option} ===", flush=True)
    with timer.time(label):
        subprocess.run([str(BINARY)], cwd=WORK_DIR, check=True)


run_stage(1, "fluxtube_precomputation")
flux_file = WORK_DIR / "fluxtubevolume.dat"
if not flux_file.exists() or flux_file.stat().st_size == 0:
    sys.exit(f"FAIL: {flux_file} not produced by stage 1")

run_stage(2, "transport_coefficient")
out_file = WORK_DIR / "nustar_diffcoef_std.dat"
if not out_file.exists() or out_file.stat().st_size == 0:
    sys.exit(f"FAIL: {out_file} not produced by stage 2")

with out_file.open() as fh:
    rows = [line.split() for line in fh if line.strip()]
if len(rows) != 1 or len(rows[0]) != 3:
    sys.exit(f"FAIL: expected one row of three numbers, got {rows}")
nu_star, d11, sigma = (float(x) for x in rows[0])
print(f"OBSERVED: nu*={nu_star:g}  D11={d11:g}  sigma={sigma:g}")

# Quantitative comparison against committed reference.
# Fluxtube geometry is deterministic: tight tolerance.
compare_numeric_file(
    flux_file,
    REFERENCE_DIR / "fluxtubevolume.dat",
    rtol=1.0e-12,
    atol=1.0e-14,
)
# D11 and sigma are reduction sums under OpenMP; floating-point order can
# differ by a few ULP. Use a relaxed relative tolerance.
compare_numeric_file(
    out_file,
    REFERENCE_DIR / "nustar_diffcoef_std.dat",
    rtol=1.0e-8,
    atol=1.0e-14,
)
compare_numeric_file(
    WORK_DIR / "n_lost_particles.dat",
    REFERENCE_DIR / "n_lost_particles.dat",
    rtol=0.0,
    atol=0.0,
)

with (WORK_DIR / "n_lost_particles.dat").open() as fh:
    n_lost = int(fh.read())
with (WORK_DIR / "loss_summary.dat").open() as fh:
    loss_rows = [line.split() for line in fh if line.strip()]
if len(loss_rows) != 1 or len(loss_rows[0]) != 3:
    sys.exit(f"FAIL: unexpected loss summary: {loss_rows}")
loss_nu_star, summary_lost, summary_total = loss_rows[0]
if float(loss_nu_star) != nu_star or int(summary_lost) != n_lost or int(summary_total) != 3:
    sys.exit(f"FAIL: inconsistent loss summary: {loss_rows[0]}")
with (WORK_DIR / "lost_particle_events.dat").open() as fh:
    event_rows = [line.split() for line in fh if line.strip()]
if len(event_rows) != n_lost:
    sys.exit(f"FAIL: expected {n_lost} loss events, got {len(event_rows)}")

with (WORK_DIR / "transport_metadata.dat").open() as fh:
    metadata_rows = [line.split() for line in fh if line.strip() and not line.startswith("#")]
if len(metadata_rows) != 1 or len(metadata_rows[0]) != 11:
    sys.exit(f"FAIL: unexpected transport metadata: {metadata_rows}")
meta = metadata_rows[0]
meta_nu_star, radius, aiota, q_saf, speed, collision_frequency, time_step = (
    float(value) for value in meta[:7]
)
if abs(meta_nu_star - nu_star) > 1.0e-15:
    sys.exit("FAIL: metadata nu* differs from transport output")
reconstructed_nu_star = radius * collision_frequency / (abs(aiota) * speed)
if abs(reconstructed_nu_star - nu_star) > 1.0e-14:
    sys.exit("FAIL: metadata does not obey standard nu* normalization")
if abs(q_saf * aiota - 1.0) > 1.0e-12:
    sys.exit("FAIL: metadata q and iota are not reciprocal")
if int(meta[10]) != 5:
    sys.exit("FAIL: metadata does not retain the configured field periods")

with (WORK_DIR / "particle_histories.dat").open() as fh:
    history_rows = [line.split() for line in fh if line.strip()]
if not history_rows or any(len(row) != 5 for row in history_rows):
    sys.exit("FAIL: malformed particle history output")
for row in history_rows:
    row_nu_star, marker, step, physical_time, _ = row
    if float(row_nu_star) != nu_star:
        sys.exit("FAIL: history row carries the wrong nu*")
    if abs(float(physical_time) - int(step) * time_step) > 1.0e-12:
        sys.exit(f"FAIL: inconsistent time for marker {marker}, step {step}")

print(f"PASS: nu*={nu_star:g}  D11={d11:g}  sigma={sigma:g} (matches reference)")


def snapshot_outputs(label: str, filenames: tuple[str, ...]) -> Path:
    snapshot = WORK_DIR / label
    snapshot.mkdir()
    for filename in filenames:
        shutil.copyfile(WORK_DIR / filename, snapshot / filename)
    return snapshot


def compare_outputs(reference: Path, filenames: tuple[str, ...]) -> None:
    for filename in filenames:
        observed = WORK_DIR / filename
        for line in observed.read_text().splitlines():
            if line.strip() and not line.lstrip().startswith("#") and not all(
                math.isfinite(float(token)) for token in line.split()
            ):
                sys.exit(f"FAIL: nonfinite value in {filename}")
        exact = filename in ("n_lost_particles.dat", "loss_summary.dat", "lost_particle_events.dat")
        compare_numeric_file(
            observed, reference / filename,
            rtol=0.0 if exact else 1.0e-8,
            atol=0.0 if exact else 1.0e-12,
        )


# The same seeded marker histories check the actual transport consumer, not
# merely the local arithmetic used by both backends.
transport_outputs = (
    "nustar_diffcoef_std.dat", "particle_histories.dat", "n_lost_particles.dat",
    "loss_summary.dat", "lost_particle_events.dat", "transport_metadata.dat",
)
scalar_reference = snapshot_outputs("scalar_lorentz", transport_outputs)
mono["transpcoefnml"]["transport_backend"] = COMPARE_BACKEND
run_stage(2, f"{COMPARE_BACKEND}_lorentz")
compare_outputs(scalar_reference, transport_outputs)

mono["transpcoefnml"]["boole_collisions"] = False
mono["transpcoefnml"]["transport_backend"] = "cpu"
run_stage(2, "scalar_collisionless")
scalar_reference = snapshot_outputs("scalar_collisionless", transport_outputs)
mono["transpcoefnml"]["transport_backend"] = COMPARE_BACKEND
run_stage(2, f"{COMPARE_BACKEND}_collisionless")
compare_outputs(scalar_reference, transport_outputs)

# Option 4 shares the tracer but not the monoenergetic dispatcher. Exercise
# its streamed moment path without storing a particle-by-time matrix.
mono["transpcoefnml"]["nt_steps_numerical_diff"] = 32
mono["transpcoefnml"]["total_mc_time"] = 32 * time_step
mono["transpcoefnml"]["idiffcoef_output"] = 3
mono["transpcoefnml"]["filename_numerical_diff_coef"] = "numerical_diffcoef.dat"
mono["transpcoefnml"]["filename_delta_s_squared"] = "mean_square_displacement.dat"
mono["transpcoefnml"]["filename_std_dvt_delta_s_squared"] = "mean_square_uncertainty.dat"
mono["transpcoefnml"]["boole_psi_mat"] = False
mono["transpcoefnml"]["boole_write_particle_histories"] = False
mono["transpcoefnml"]["transport_backend"] = "cpu"
run_stage(4, "scalar_numerical_diffusion")
numerical_outputs = (
    "numerical_diffcoef.dat", "mean_square_displacement.dat",
    "mean_square_uncertainty.dat", "n_lost_particles.dat",
)
scalar_reference = snapshot_outputs("scalar_numerical_diffusion", numerical_outputs)
mono["transpcoefnml"]["transport_backend"] = COMPARE_BACKEND
run_stage(4, f"{COMPARE_BACKEND}_numerical_diffusion")
compare_outputs(scalar_reference, numerical_outputs)
mono["transpcoefnml"]["boole_psi_mat"] = True
run_stage(4, f"{COMPARE_BACKEND}_numerical_diffusion_histories")
compare_outputs(scalar_reference, numerical_outputs)

# Move the published equilibrium's angular samples close to the outer boundary.
# Several markers then leave after retaining a nonempty history.
saved_mono = dict(mono["transpcoefnml"])
saved_fluxtube = flux_file.read_text()
flux_rows = saved_fluxtube.splitlines()
for row_index in range(2, len(flux_rows)):
    values = flux_rows[row_index].split()
    values[0] = "0.99999"
    flux_rows[row_index] = " ".join(values)
flux_file.write_text("\n".join(flux_rows) + "\n")
mono["transpcoefnml"].update(
    n_particles=16, boole_collisions=False, lambda_start=1.0,
    boole_write_particle_histories=True, idiffcoef_output=2,
    transport_backend="cpu",
)
run_stage(2, "scalar_boundary_loss")
lost = int((WORK_DIR / "n_lost_particles.dat").read_text())
if not 0 < lost < 16:
    sys.exit("FAIL: boundary-loss fixture needs both lost and surviving markers")
loss_outputs = (
    "particle_histories.dat", "n_lost_particles.dat", "loss_summary.dat",
    "transport_metadata.dat", "mean_square_displacement.dat",
    "mean_square_uncertainty.dat",
)
scalar_reference = snapshot_outputs(
    "scalar_boundary_loss", loss_outputs + ("lost_particle_events.dat",)
)
mono["transpcoefnml"]["transport_backend"] = COMPARE_BACKEND
run_stage(2, f"{COMPARE_BACKEND}_boundary_loss")
compare_outputs(scalar_reference, loss_outputs)
reference_events = [line.split() for line in (scalar_reference / "lost_particle_events.dat").read_text().splitlines()]
observed_events = [line.split() for line in (WORK_DIR / "lost_particle_events.dat").read_text().splitlines()]
if len(observed_events) != lost:
    sys.exit("FAIL: boundary-loss event count differs")
for observed, reference in zip(observed_events, reference_events):
    if observed[1:4] != reference[1:4]:
        sys.exit("FAIL: boundary-loss marker, reason or measurement step differs")
compare_numeric_file(
    WORK_DIR / "lost_particle_events.dat", scalar_reference / "lost_particle_events.dat",
    rtol=1.0e-8, atol=1.0e-12,
)

# Reconstruct the estimator from each marker's retained trajectory prefix.
initial_psi = {}
samples = {}
for line in (WORK_DIR / "particle_histories.dat").read_text().splitlines():
    _, marker, step, sample_time, psi = (float(value) for value in line.split())
    marker, step = int(marker), int(step)
    if step == 0:
        initial_psi[marker] = psi
    else:
        samples.setdefault(step, []).append((sample_time, (psi-initial_psi[marker])**2))
mean_rows = [line.split() for line in (WORK_DIR / "mean_square_displacement.dat").read_text().splitlines()][:-1]
error_rows = [line.split() for line in (WORK_DIR / "mean_square_uncertainty.dat").read_text().splitlines()][:-1]
if len(mean_rows) != len(samples) or len(error_rows) != len(samples):
    sys.exit("FAIL: moment output does not cover retained measurement steps")
for step, mean_row, error_row in zip(sorted(samples), mean_rows, error_rows):
    values = samples[step]
    second = math.fsum(value[1] for value in values)/len(values)
    fourth = math.fsum(value[1]**2 for value in values)/len(values)
    uncertainty = math.sqrt(max(0.0, fourth-second**2)/len(values))
    if not math.isclose(float(mean_row[0]), values[0][0], rel_tol=1.0e-12, abs_tol=1.0e-18):
        sys.exit("FAIL: retained-history measurement time differs")
    if not math.isclose(float(mean_row[1]), second, rel_tol=1.0e-8, abs_tol=1.0e-20):
        sys.exit("FAIL: second moment omits a retained trajectory prefix")
    if not math.isclose(float(error_row[1]), uncertainty, rel_tol=1.0e-8, abs_tol=1.0e-20):
        sys.exit("FAIL: uncertainty omits a retained trajectory prefix")
flux_file.write_text(saved_fluxtube)
mono["transpcoefnml"] = saved_mono

for key, invalid, message in (
    ("transport_backend", "unknown", "transport_backend must be"),
    ("i_integrator_type", 2, "batch transport requires i_integrator_type=1"),
    ("boole_random_precalc", False, "batch transport requires precomputed random numbers"),
):
    previous = mono["transpcoefnml"][key]
    mono["transpcoefnml"][key] = invalid
    mono.write(str(WORK_DIR / "mono_energetic_transp_coef.inp"), force=True)
    rejected = subprocess.run([str(BINARY)], cwd=WORK_DIR, capture_output=True, text=True)
    if rejected.returncode == 0 or message not in rejected.stdout + rejected.stderr:
        sys.exit(f"FAIL: unsupported {key}={invalid!r} was not rejected clearly")
    mono["transpcoefnml"][key] = previous
mono.write(str(WORK_DIR / "mono_energetic_transp_coef.inp"), force=True)

gorilla["gorillanml"]["boole_helical_pert"] = True
gorilla.write(str(WORK_DIR / "gorilla.inp"), force=True)
rejected = subprocess.run([str(BINARY)], cwd=WORK_DIR, capture_output=True, text=True)
if rejected.returncode == 0 or "batch transport does not support helical perturbations" not in (
    rejected.stdout + rejected.stderr
):
    sys.exit("FAIL: helical perturbations were not rejected clearly")
gorilla["gorillanml"]["boole_helical_pert"] = False
gorilla.write(str(WORK_DIR / "gorilla.inp"), force=True)

print(f"PASS: scalar/{COMPARE_BACKEND} seeded Lorentz and collisionless transport agree")
