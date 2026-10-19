#!/usr/bin/env python3
"""Switch gates of the collisional background drive (ntv_fax_file).

Same analytic circular case as test_canonical_ntv_lorentz.py, with an m = 0
harmonic (the drive needs b.grad = b^phi d/dphi). One thread, so that runs
with the same seed are reproducible marker by marker.

Gates:
  (a) a zero table changes nothing: the marker table is byte-identical to the
      run without ntv_fax_file (default off; the off path is the 565e693 one);
  (b) a non-zero table changes the weights but not the orbits: dh_int,
      dpphi_orbit and the end states are identical, w_end and tq_deltaf differ;
  (c) ntv_fax_file with m != 0 is rejected at initialisation.
"""

import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import test_canonical_ntv_lorentz as base  # noqa: E402

N_THETA, N_ETA = 7, 5


def write_table(path, value):
    import math
    vts = [-math.pi + 2.0 * math.pi * j / (N_THETA - 1) for j in range(N_THETA)]
    etas = [0.0, 0.5, 1.0, 1.2, 2.0]
    rows = [" ".join(f"{value * (1.0 + 0.1 * j):.16e}" for _ in etas)
            for j in range(N_THETA)]
    neg = [" ".join(f"{-value * (1.0 + 0.1 * j):.16e}" for _ in etas)
           for j in range(N_THETA)]
    path.write_text("# test table\n" + f"{N_THETA} {N_ETA}\n"
                    + " ".join(f"{v:.16e}" for v in vts) + "\n"
                    + " ".join(f"{e:.16e}" for e in etas) + "\n"
                    + "\n".join(rows + neg) + "\n")


def case(name, n_particles, m_pol=0, table=None):
    case_dir = base.WORK_DIR / name
    base.prepare(case_dir, n_particles)
    base.set_nml(case_dir / "rmp_response_currents.inp", {"ntv_m_pol": m_pol})
    if table is not None:
        write_table(case_dir / "fax_table.dat", table)
        inp = case_dir / "rmp_response_currents.inp"
        text = inp.read_text()
        key = "  ntv_w_reset_burn = 0.0 ,\n"
        if text.count(key) != 1:
            sys.exit("cannot place ntv_fax_file")
        inp.write_text(text.replace(key, key + "  ntv_fax_file = 'fax_table.dat' ,\n"))
    return case_dir


def main():
    base.WORK_DIR.mkdir(parents=True, exist_ok=True)
    n_particles = int(os.environ.get("NTV_N_PARTICLES", "64"))
    failures = []

    off = case("fax_off", n_particles)
    base.run(off, 1)
    zero = case("fax_zero", n_particles, table=0.0)
    base.run(zero, 1)
    same = (off / "canonical_ntv_markers.dat").read_bytes() == \
        (zero / "canonical_ntv_markers.dat").read_bytes()
    print(f"(a) zero table: marker table byte-identical = {same}")
    if not same:
        failures.append("zero table changed the marker table")
    if "collisional background" not in (zero / "run.log").read_text():
        failures.append("table not reported as loaded")

    on = case("fax_on", n_particles, table=0.3)
    base.run(on, 1)
    m0, m1 = base.markers(off), base.markers(on)
    orbit_keys = ("dh_int", "dpphi_orbit", "pphi_end", "s_end", "t_end", "n_push")
    orbits_same = all(m0[k] == m1[k] for k in orbit_keys)
    dw = max(abs(a - b) for a, b in zip(m0["w_end"], m1["w_end"]))
    wscale = max(abs(a) for a in m0["w_end"])
    dtq = sum(m1["tq_deltaf"]) - sum(m0["tq_deltaf"])
    print(f"(b) orbits identical = {orbits_same}; max |dw_end| / max|w_end| = "
          f"{dw / wscale:.3e}; sum tq change = {dtq:.3e}")
    if not orbits_same:
        failures.append("the drive changed the orbits")
    if not dw > 1.0e-6 * wscale or dtq == 0.0:
        failures.append("a non-zero table left the weights unchanged")

    bad = case("fax_m3", 8, m_pol=3, table=0.3)
    rc = base.run(bad, 1, expect_fail=True)
    print(f"(c) ntv_fax_file with m = 3: rc = {rc}")
    if rc == 0 or "needs ntv_m_pol = 0" not in (bad / "run.log").read_text():
        failures.append("m != 0 was not rejected")

    if failures:
        print("FAIL: " + "; ".join(failures))
        return 1
    print("PASS: canonical NTV collisional background switch gates")
    return 0


if __name__ == "__main__":
    sys.exit(main())
