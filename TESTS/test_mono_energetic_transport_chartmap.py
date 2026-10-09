#!/usr/bin/env python3
"""
CTest driver: mono-energetic transport on a Boozer chartmap (grid_kind = 6).

Writes a manufactured chartmap with constant iota = 0.4 and two field
periods, runs i_option = 1 (fluxtube volume) and i_option = 2 (single nu*),
and checks that transport_metadata.dat carries the chartmap iota and the
standard nu* normalization. D11 itself is not compared.

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

FIXTURE_IOTA = 0.4
FIXTURE_NFP = 2


def env_path(name: str) -> Path:
    value = os.environ.get(name)
    if not value:
        sys.exit(f"missing required env var: {name}")
    return Path(value)


APPLETS_ROOT = env_path("APPLETS_ROOT")
GORILLA_ROOT = env_path("GORILLA_ROOT")
BINARY = env_path("GORILLA_APPLETS_BIN")
FIXTURE_WRITER = env_path("CHARTMAP_FIXTURE_BIN")
WORK_DIR = env_path("WORK_DIR")

if WORK_DIR.exists():
    shutil.rmtree(WORK_DIR)
WORK_DIR.mkdir(parents=True)

subprocess.run([str(FIXTURE_WRITER), "chartmap.nc"], cwd=WORK_DIR, check=True)

gorilla = f90nml.read(str(GORILLA_ROOT / "INPUT" / "gorilla.inp"))
tetra_grid = f90nml.read(str(GORILLA_ROOT / "INPUT" / "tetra_grid.inp"))
gorilla_applets = f90nml.read(str(APPLETS_ROOT / "INPUT" / "gorilla_applets.inp"))
mono = f90nml.read(str(APPLETS_ROOT / "INPUT" / "mono_energetic_transp_coef.inp"))
for nml in (gorilla, tetra_grid, gorilla_applets, mono):
    nml.end_comma = True

gorilla["gorillanml"]["eps_phi"] = 0.0
gorilla["gorillanml"]["coord_system"] = 2          # (s, theta, phi)
gorilla["gorillanml"]["ispecies"] = 1              # electron
gorilla["gorillanml"]["boole_periodic_relocation"] = True
gorilla["gorillanml"]["ipusher"] = 2               # polynomial pusher
gorilla["gorillanml"]["poly_order"] = 2
gorilla["gorillanml"]["boole_grid_for_find_tetra"] = False
gorilla["gorillanml"]["boole_adaptive_time_steps"] = False

tetra_grid["tetra_grid_nml"]["grid_kind"] = 6
tetra_grid["tetra_grid_nml"]["n1"] = 10            # ns
tetra_grid["tetra_grid_nml"]["n2"] = 10            # nphi
tetra_grid["tetra_grid_nml"]["n3"] = 10            # ntheta
tetra_grid["tetra_grid_nml"]["boole_n_field_periods"] = True
tetra_grid["tetra_grid_nml"]["sfc_s_min"] = 1.0e-1
tetra_grid["tetra_grid_nml"]["netcdf_filename"] = "chartmap.nc"

gorilla_applets["gorilla_applets_nml"]["filename_fluxtv_precomp"] = "fluxtubevolume.dat"
gorilla_applets["gorilla_applets_nml"]["filename_fluxtv_load"] = "fluxtubevolume.dat"
gorilla_applets["gorilla_applets_nml"]["start_pos_x1"] = 0.5
gorilla_applets["gorilla_applets_nml"]["start_pos_x2"] = 0.00013
gorilla_applets["gorilla_applets_nml"]["start_pos_x3"] = 0.00013
gorilla_applets["gorilla_applets_nml"]["t_step_fluxtv"] = 1.0e-3
gorilla_applets["gorilla_applets_nml"]["nt_steps_fluxtv"] = 200
gorilla_applets["gorilla_applets_nml"]["energy_ev_fluxtv"] = 3.0e3

mono["transpcoefnml"]["i_integrator_type"] = 1
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
mono["transpcoefnml"]["flight_time_multiplier"] = 1.0
mono["transpcoefnml"]["boole_psi_mat"] = False
mono["transpcoefnml"]["boole_write_particle_histories"] = False
mono["transpcoefnml"]["filename_transport_metadata"] = "transport_metadata.dat"

(WORK_DIR / "seed.inp").write_text("8\n  1 2 3 4 5 6 7 8\n")
shutil.copy(GORILLA_ROOT / "INPUT" / "field_divB0.inp", WORK_DIR / "field_divB0.inp")

gorilla.write(str(WORK_DIR / "gorilla.inp"), force=True)
tetra_grid.write(str(WORK_DIR / "tetra_grid.inp"), force=True)
mono.write(str(WORK_DIR / "mono_energetic_transp_coef.inp"), force=True)


def run_stage(i_option: int, label: str) -> None:
    gorilla_applets["gorilla_applets_nml"]["i_option"] = i_option
    gorilla_applets.write(str(WORK_DIR / "gorilla_applets.inp"), force=True)
    print(f"=== {label}: i_option={i_option} ===", flush=True)
    subprocess.run([str(BINARY)], cwd=WORK_DIR, check=True)


run_stage(1, "fluxtube_precomputation")
run_stage(2, "transport_coefficient")

with (WORK_DIR / "transport_metadata.dat").open() as fh:
    rows = [line.split() for line in fh if line.strip() and not line.startswith("#")]
if len(rows) != 1 or len(rows[0]) != 11:
    sys.exit(f"FAIL: unexpected transport metadata: {rows}")
meta = rows[0]
nu_star, radius, aiota, q_saf, speed, collision_frequency = (float(v) for v in meta[:6])

if abs(aiota - FIXTURE_IOTA) > 1.0e-12:
    sys.exit(f"FAIL: metadata iota {aiota!r} differs from chartmap iota {FIXTURE_IOTA}")
if abs(q_saf * aiota - 1.0) > 1.0e-12:
    sys.exit("FAIL: metadata q and iota are not reciprocal")
reconstructed_nu_star = radius * collision_frequency / (abs(aiota) * speed)
if abs(reconstructed_nu_star - nu_star) > 1.0e-12 * nu_star:
    sys.exit("FAIL: metadata does not obey standard nu* normalization")
if int(meta[10]) != FIXTURE_NFP:
    sys.exit(f"FAIL: metadata field periods {meta[10]} differ from chartmap {FIXTURE_NFP}")

print(f"PASS: chartmap iota={aiota:g} q={q_saf:g} nu*={nu_star:g}")
