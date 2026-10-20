#!/usr/bin/env python3
"""Eight-marker chartmap switch and spectral-derivative regression, m=3 n=2.

The table is artificial: these checks do not qualify a device f_AX background.
"""
import sys

import test_canonical_ntv_chartmap as base
from test_canonical_ntv_fax import write_table


def main():
    base.WORK_DIR.mkdir(parents=True, exist_ok=True)
    mesh = base.WORK_DIR / "axi.nc"
    pert = base.WORK_DIR / "pert.nc"
    base.write_chartmap(mesh, 0.0)
    base.write_chartmap(pert, base.B_MN)
    results = {}
    for name, table, spectrum in (("off", None, False), ("zero", 0.0, False),
                                  ("on", 0.3, False), ("spectrum", 0.3, True)):
        case = base.WORK_DIR / name
        options = dict(ntv_amplitude=1.0, ntv_bmn=base.B_MN)
        if table is not None:
            options["ntv_fax_file"] = "'fax.dat'"
        if spectrum:
            options.update(ntv_m_pol=0, ntv_m_max=3, ntv_chartmap_file="'pert.nc'")
        base.prepare(case, "chartmap", 8, "2.0d-7", False, options, mesh=mesh)
        if table is not None:
            write_table(case / "fax.dat", table)
        if spectrum:
            (case / "pert.nc").symlink_to(pert)
        results[name] = base.run(case, 1)
    off = (base.WORK_DIR / "off/canonical_ntv_markers.dat").read_bytes()
    zero = (base.WORK_DIR / "zero/canonical_ntv_markers.dat").read_bytes()
    if off != zero:
        sys.exit("zero f_AX table changed marker output")
    orbit_keys = ("dh_int", "dpphi_orbit", "pphi_end", "s_end", "t_end", "n_push")
    for key in orbit_keys:
        if results["off"][key] != results["on"][key]:
            sys.exit(f"f_AX changed orbit column {key}")
    delta = max(abs(a - b) for a, b in zip(results["off"]["w_end"],
                                         results["on"]["w_end"]))
    if delta == 0.0:
        sys.exit("nonzero f_AX table did not change weights")
    err = max(abs(a - b) for a, b in zip(results["on"]["w_end"],
                                       results["spectrum"]["w_end"]))
    scale = max(abs(v) for v in results["on"]["w_end"])
    print(f"m=3 f_AX max weight change = {delta:.12e}")
    print(f"spectrum / explicit harmonic max weight error = {err / scale:.12e}")
    if err > 1.0e-6 * scale:
        sys.exit("chartmap spectrum disagrees with explicit m=3 f_AX source")
    print("PASS: chartmap f_AX zero switch, fixed orbits, nonzero-m spectrum")
    return 0


if __name__ == "__main__":
    sys.exit(main())
