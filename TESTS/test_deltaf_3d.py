#!/usr/bin/env python3
"""Native 3D orbit, AX control-variate, f_AX, loss and periodicity gates."""
import json
import math
import os
from pathlib import Path
import re
import subprocess
import test_canonical_ntv_chartmap as base

ROOT = base.WORK_DIR
ROOT.mkdir(parents=True, exist_ok=True)


def chart(name, delta=0., n=18, nfp=18, rmp=0.):
    path = ROOT / (name + '.nc')
    subprocess.run([str(base.WRITER), str(path), '1', str(delta), '0', str(n),
                    '33', '96', '96' if nfp == 1 else '16', '3', '2' if rmp else '0',
                    str(rmp), str(nfp), '1'], check=True, stdout=subprocess.DEVNULL)
    return path


def table(path, amplitude):
    """Artificial total f_AX/f_M=a*lambda*cos(vartheta); one species/energy basis."""
    vt = [-math.pi + j * math.pi / 16 for j in range(33)]
    lam = [-1 + j / 16 for j in range(33)]
    data = [9.1094e-28, -1., 1.6022e-9]  # core electron and T=1keV, CGS
    data += [amplitude, 0., 0.] + [0., 20.] + [0.] + [1.] + vt + lam
    data += [l * math.cos(t) if k == 0 else 0. for k in range(3) for t in vt for l in lam]
    data += [math.cos(t) if k == 0 else 0. for k in range(3) for t in vt for l in lam]
    path.write_text('FAX_ENERGY_V1\n33 33 1 1 0 1 T\n' + '\n'.join(map(repr, data)) + '\n')
    return path


def prepare(name, mesh, extra=None, n=64, trace='2e-5'):
    case = ROOT / name
    ntv = {'ntv_3d_background': '.true.', 'ntv_mesh_phi0': '.true.',
           'ntv_amplitude': 0, 'ntv_m_pol': 0, 'ntv_n_tor': 0}
    ntv.update(extra or {})
    base.prepare(case, 'chartmap', n, trace, False, ntv, mesh, n2=48)
    base.set_nml(case / 'tetra_grid.inp', {'boole_n_field_periods': '.true.', 'n1': 16})
    return case


def run(case, threads=1):
    m = base.run(case, threads)
    assert all(math.isfinite(v) for k, values in m.items() if k != 'volume' for v in values)
    assert all(t >= 0 for t in m['t_end'])
    text = (case / 'run.log').read_text()
    m['seconds_per_marker'] = float(re.search(r'NTV orbit wall seconds per marker:\s*(\S+)',text)[1])
    return m


def same_orbits(a,b):
    for key in ('dpphi_orbit','dpphi_coll','pphi_end','s_end','t_end','n_push','lost'):
        assert a[key] == b[key], f'orbits changed: {key}'


def reject(case, text):
    result = subprocess.run([str(base.BIN)],cwd=case,capture_output=True,text=True,
                            env=dict(os.environ,OMP_NUM_THREADS='1'))
    assert result.returncode != 0 and text in result.stdout+result.stderr


def main():
    axis = chart('axis')
    ripple = chart('ripple',.01)
    coupled = chart('coupled',.01,n=4,nfp=1,rmp=10.)
    on = run(prepare('axis_3d',axis))
    off = run(prepare('axis_legacy',axis,{'ntv_3d_background':'.false.'}))
    same_orbits(on,off)
    assert max(map(abs,on['w_end'])) == 0
    # Legacy harmonic labels must not restrict a total-field spectrum.
    a = run(prepare('ripple_3d',ripple,{'ntv_n_tor':2,'ntv_m_pol':3}))
    assert max(map(abs,a['w_end'])) > 1e-10
    zero = table(ROOT/'fax_zero.dat',0.)
    finite = table(ROOT/'fax_finite.dat',.01)
    z = run(prepare('ripple_zero',ripple,{'ntv_fax_file':repr(str(zero))}))
    f = run(prepare('ripple_fax',ripple,{'ntv_fax_file':repr(str(finite))}))
    same_orbits(a,z)
    same_orbits(a,f)
    # A zero native table means f_AX=0; the no-table AX reference retains G.
    assert max(abs(x-y) for x,y in zip(a['w_end'],z['w_end'])) > 1e-10
    dw = max(abs(x-y) for x,y in zip(z['w_end'],f['w_end']))
    assert dw > 1e-9
    faxis = run(prepare('axis_fax',axis,{'ntv_fax_file':repr(str(finite))}))
    assert max(map(abs,faxis['w_end'])) == 0
    coupled_case = prepare('coupled_3d',coupled,n=64)
    base.set_nml(coupled_case/'tetra_grid.inp',{'n3':96})
    c = run(coupled_case)
    assert max(map(abs,c['w_end'])) > 1e-10
    # A domain exit must freeze all signed integrals after loss.
    loss_case = prepare('loss',coupled,{'ntv_n_check':4,'ntv_check_dt':5e-6})
    base.set_nml(loss_case/'tetra_grid.inp',{'n3':96})
    base.set_nml(loss_case/'rmp_response_currents.inp',{'s_outer_cut':.6})
    loss = run(loss_case)
    lost = [i for i,v in enumerate(loss['lost']) if v]
    assert lost and all(loss['loss_kind'][i] == 3 for i in lost)
    assert all(loss['t_end'][i] < 2e-5 for i in lost)
    # Thread independence covers both weight arrays and native-table lookup.
    f4 = run(prepare('ripple_fax_4',ripple,{'ntv_fax_file':repr(str(finite))}),4)
    for key in ('w_end','tq_deltaf','s_end','P'):
        assert f[key] == f4[key], key
    invalid = prepare('bad_period',coupled)
    base.set_nml(invalid/'tetra_grid.inp',{'boole_n_field_periods':'.false.',
                                         'n_field_periods_manual':18})
    reject(invalid,'mesh field periods must match chartmap_nfp')
    invalid = prepare('bad_drive',ripple,{'ntv_amplitude':1})
    reject(invalid,'perturbations belong in the mesh')
    result = {'axis_max_w':max(map(abs,on['w_end'])),
              'ripple_max_w':max(map(abs,a['w_end'])),
              'fax_max_weight_change':dw,'outer_losses':len(lost),
              'coupled_max_w':max(map(abs,c['w_end'])),
              'seconds_per_marker':{k:m['seconds_per_marker'] for k,m in
                                    [('axis',on),('ripple',a),('fax',f),('coupled',c)]}}
    (ROOT/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
    print('PASS: 3D total-field orbits, AX/f_AX reduction, losses and period rejection')


if __name__ == '__main__':
    main()
