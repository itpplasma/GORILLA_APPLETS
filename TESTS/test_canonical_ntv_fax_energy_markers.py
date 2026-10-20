"""The real energy-table route changes weights and preserves the marker orbits."""
import os
from pathlib import Path
import sys
import test_canonical_ntv_lorentz as base

root = base.WORK_DIR
root.mkdir(parents=True, exist_ok=True)
table = Path(sys.argv[1]).resolve()
for name in ('off','energy'):
    case = root/name
    base.prepare(case,8)
    base.set_nml(case/"rmp_response_currents.inp", {"ntv_m_pol": 0})
    if name == 'energy':
        inp = case/'rmp_response_currents.inp'
        s = inp.read_text()
        pos = s.rfind('/')
        if pos < 0:
            raise SystemExit('missing namelist terminator')
        inp.write_text(s[:pos]+f"  ntv_fax_file = '{table}',\n"+s[pos:])
    base.run(case,1)
a,b = base.markers(root/'off'),base.markers(root/'energy')
for key in ('dh_int','dpphi_orbit','pphi_end','s_end','t_end','n_push'):
    if a[key] != b[key]:
        raise SystemExit(f'FAIL: energy background changed orbits: {key}')
import math
if not all(math.isfinite(x) for x in b['w_end']):
    raise SystemExit('FAIL: energy background produced nonfinite weights')
dw = max(abs(x-y) for x,y in zip(a['w_end'],b['w_end']))
if dw == 0:
    raise SystemExit('FAIL: nonzero native energy background did not change weights')
print(f'PASS: native energy table; orbits identical; max weight change {dw:.17e}')
