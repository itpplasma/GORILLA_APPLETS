#!/usr/bin/env python3
"""
CTest driver: ExB rotation from v_E in the mono-energetic transport applet on
a Boozer chartmap (grid_kind = 6).

The manufactured chartmap (write_chartmap_fixture.f90) has |B| = B_00 =
1000 G and psi_tor,a = 2e8 G cm^2. The applet defines v_E = c*E_r/(B_00*v)
with r = sqrt(s*psi_tor,a/(pi*B_00)). For Phi = Phi(s), the guiding-centre
ExB drift changes the Boozer field-line label alpha = theta - iota*phi at
the rate c*dPhi/dpsi_tor on every surface, while parallel motion leaves
alpha unchanged. The expected rate is therefore

    dalpha/dt = v_E*v*B_00/(psi_tor,a*ds/dr),  ds/dr = 2*sqrt(s*pi*B_00/psi_tor,a).

The test follows a pitch-0 orbit (deeply trapped limit; the fixture has
constant |B|) and co- and counter-passing orbits (pitch +-0.9) with v_E = 0
and v_E /= 0, and checks the difference against this rate. It also checks that v_E = 0
keeps eps_Phi from gorilla.inp, as before this option existed.
"""
import math
import os
import shutil
import subprocess
import sys
from pathlib import Path

try:
    import f90nml
except ImportError:
    sys.exit("f90nml is required to run this test (pip install f90nml)")

B00 = 1000.0
TORFLUX = 2.0e8
V_E = 1.0e-2
REL_TOL = 1.0e-3
# (pitch, time step [s], steps). The fixture's covariant B_theta, B_phi do
# not match its geometry, so its Jacobian is small and parallel streaming
# advances the Boozer angles at ~6e9 rad/s; passing orbits need tiny steps
# for unambiguous angle unwrapping. Pitch 0 never streams: |B| is constant,
# so there is no mirror force and the orbit only drifts (deeply trapped limit).
ORBITS = ((0.0, 1.0e-7, 200), (0.9, 5.0e-12, 4000), (-0.9, 5.0e-12, 4000))


def env_path(name: str) -> Path:
    value = os.environ.get(name)
    if not value:
        sys.exit(f"missing required env var: {name}")
    return Path(value)


APPLETS_ROOT = env_path("APPLETS_ROOT")
GORILLA_ROOT = env_path("GORILLA_ROOT")
BINARY = env_path("GORILLA_APPLETS_BIN")
DRIVER = env_path("EXB_DRIVER_BIN")
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
gorilla["gorillanml"]["coord_system"] = 2
gorilla["gorillanml"]["ispecies"] = 1              # electron
gorilla["gorillanml"]["boole_periodic_relocation"] = True
gorilla["gorillanml"]["ipusher"] = 2
gorilla["gorillanml"]["poly_order"] = 2
gorilla["gorillanml"]["boole_grid_for_find_tetra"] = False
gorilla["gorillanml"]["boole_adaptive_time_steps"] = False

tetra_grid["tetra_grid_nml"]["grid_kind"] = 6
tetra_grid["tetra_grid_nml"]["n1"] = 10
tetra_grid["tetra_grid_nml"]["n2"] = 10
tetra_grid["tetra_grid_nml"]["n3"] = 10
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

mono["transpcoefnml"]["energy_ev"] = 3.0e3
mono["transpcoefnml"]["boole_collisions"] = False

shutil.copy(GORILLA_ROOT / "INPUT" / "field_divB0.inp", WORK_DIR / "field_divB0.inp")
gorilla.write(str(WORK_DIR / "gorilla.inp"), force=True)
tetra_grid.write(str(WORK_DIR / "tetra_grid.inp"), force=True)

gorilla_applets["gorilla_applets_nml"]["i_option"] = 1
gorilla_applets.write(str(WORK_DIR / "gorilla_applets.inp"), force=True)
subprocess.run([str(BINARY)], cwd=WORK_DIR, check=True, stdout=subprocess.DEVNULL)


def run_driver(v_e: float, pitch: float, t_step: float, n_steps: int) -> dict:
    mono["transpcoefnml"]["v_e"] = v_e
    mono.write(str(WORK_DIR / "mono_energetic_transp_coef.inp"), force=True)
    out = subprocess.run(
        [str(DRIVER), repr(pitch), repr(t_step), str(n_steps)],
        cwd=WORK_DIR, check=True, capture_output=True, text=True,
    ).stdout
    result = {}
    for line in out.splitlines():
        if line.startswith("RESULT"):
            _, key, value = line.split()
            result[key] = float(value)
    return result


failures = []
for pitch, t_step, n_steps in ORBITS:
    ref = run_driver(0.0, pitch, t_step, n_steps)
    exb = run_driver(V_E, pitch, t_step, n_steps)
    if ref["eps_Phi"] != 0.0:
        failures.append(f"v_E = 0 changed eps_Phi from gorilla.inp: {ref['eps_Phi']}")
    s0 = exb["s0"]
    ds_dr = 2.0 * math.sqrt(s0 * math.pi * B00 / TORFLUX)
    expected = V_E * exb["vmod"] * B00 / (TORFLUX * ds_dr)
    measured = exb["alpha_rate"] - ref["alpha_rate"]
    rel = abs(measured - expected) / abs(expected)
    print(f"pitch={pitch:+.2f} s0={s0:.4f} turns={exb['phi_turns']:.2f} "
          f"ds={exb['s_excursion']:.2e} alpha_rate(v_E=0)={ref['alpha_rate']:.6e} "
          f"ExB rate={measured:.6e} expected={expected:.6e} rel_err={rel:.2e}")
    if not rel < REL_TOL:
        failures.append(f"pitch {pitch}: ExB rate {measured} differs from {expected}")

if failures:
    sys.exit("FAIL: " + "; ".join(failures))
print("PASS: chartmap ExB rotation matches v_E")
