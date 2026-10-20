#!/usr/bin/env python3
"""A native ripple-trapped orbit and an AX orbit in a large-aspect-ratio chart."""
import math
import os
from pathlib import Path
import subprocess
import test_canonical_ntv_chartmap as base

ROOT = base.WORK_DIR
ROOT.mkdir(parents=True,exist_ok=True)
PROBE = Path(os.environ['ORBIT_PROBE_BIN'])


def prepare(name, delta, nr=33, nt=96, nz=32, r0=7000):
    mesh = ROOT/(name+'.nc')
    subprocess.run([str(base.WRITER),str(mesh),'1',str(delta),'0','18',str(nr),str(nt),
                    str(nz),'0','0','0','18','1',str(r0)],check=True,stdout=subprocess.DEVNULL)
    case = ROOT/name
    base.prepare(case,'chartmap',1,'1e-5',False,{'ntv_amplitude':0},mesh,n2=nt)
    base.set_nml(case/'tetra_grid.inp',{'boole_n_field_periods':'.true.','n1':32,'n3':nz})
    return case


def probe(case, n=1, steps=10000, dt=1e-6, nu=0., energy=1., theta=1.57, pitch=.01, seed=71):
    command = [str(PROBE),str(n),str(steps),str(dt),str(nu),str(energy),'0.6',str(theta),
               str(math.pi/18),str(pitch),str(seed)]
    with (case/'probe.log').open('w') as log:
        subprocess.run(command,cwd=case,stdout=log,stderr=subprocess.STDOUT,check=True,
                       env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'))
    return [[float(v) for v in line.split()] for line in (case/'probe.dat').read_text().splitlines()
            if line and not line.startswith('#')]


def main():
    well = probe(prepare('well',.004))[0]
    banana = probe(prepare('axisymmetric',0.))[0]
    assert well[1] == banana[1] == 0
    assert well[6] > 5 and well[8] < 2*math.pi/18, well
    assert banana[8] > 2*math.pi/18, banana
    print(f'Native ripple orbit: {int(well[6])} mirrors, phi span={well[8]:.8g} rad')
    print(f'Native AX orbit: {int(banana[6])} mirrors, phi span={banana[8]:.8g} rad')
    print('PASS: total-field native ripple trapping')


if __name__ == '__main__':
    main()
