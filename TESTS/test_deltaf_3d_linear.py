#!/usr/bin/env python3
"""Small RMP: compare the 3D local DKE with #82's canonical linear estimator."""
import json
import math
from pathlib import Path
import test_canonical_ntv_chartmap as b


def main():
    b.WORK_DIR.mkdir(parents=True,exist_ok=True)
    axis = b.WORK_DIR/'axis.nc'
    b.write_chartmap(axis,0.)
    results = {}
    for amplitude, mode in [(.01,'linear'),(.01,'dke'),(.02,'dke')]:
        mesh = b.WORK_DIR/f'mesh{amplitude}.nc'
        b.write_chartmap(mesh,amplitude)
        case = b.WORK_DIR/f'{mode}{amplitude}'
        ntv = {'ntv_3d_background':'.true.' if mode=='dke' else '.false.',
               'ntv_mesh_phi0':'.true.','ntv_amplitude':0 if mode=='dke' else 1,
               'ntv_bmn':amplitude}
        b.prepare(case,'chartmap',512,'2e-5',False,ntv,mesh if mode=='dke' else axis,n2=48)
        b.set_nml(case/'tetra_grid.inp',{'boole_n_field_periods':'.true.','n1':24,'n3':48})
        m = b.run(case,1)
        samples = [p*t for p,t in zip(m['P'],m['tq_deltaf'])]
        mean,se = b.mean_err(samples)
        results[f'{mode}_{amplitude}'] = {'mean':mean,'SE':se}
    relative = abs(results['dke_0.01']['mean']/results['linear_0.01']['mean']-1)
    quadratic = results['dke_0.02']['mean']/results['dke_0.01']['mean']
    # The two routes use raw-chart and tetra/harmonic derivatives respectively,
    # and canonical G is retained to first drift-kinetic order in the new route.
    assert relative < .05, relative
    assert abs(quadratic/4-1) < 1e-3, quadratic
    results['relative_linear_difference'] = relative
    results['amplitude_doubling_ratio'] = quadratic
    (b.WORK_DIR/'verification.json').write_text(json.dumps(results,indent=2)+'\n')
    print(json.dumps(results,indent=2))
    print('PASS: small-RMP linear limit (5% mesh/ordering tolerance), quadratic response')


if __name__ == '__main__':
    main()
