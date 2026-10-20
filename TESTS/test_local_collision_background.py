"""Independent physical deflection-rate oracle for profile-fed D/T/e backgrounds.

Native CGS constants and pair Coulomb-log convention are retained. The exact
collision-frequency coefficient is checked without the native velocity kernel.
Deflection rates additionally use the declared native kernel approximation.
This is a rate/input oracle, not full-FP or transport qualification.
"""
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

DRIVER, WRITER, WORK = map(Path, sys.argv[1:])
TEMPLATE = Path(__file__).parent / 'rmp_fm_stationarity'
NE = float(os.environ.get('DENSITY_ORACLE_NE', '1.2e14'))
TI = float(os.environ.get('DENSITY_ORACLE_TI', '3500'))
TE = float(os.environ.get('DENSITY_ORACLE_TE', '4700'))
MASSES = [3.345228e-24, 5.017842e-24, 9.1091e-28]


def prepare(name, explicit=True, fractions=(.5, .5), electron=1.):
    root = WORK / name
    if root.exists():
        shutil.rmtree(root)
    shutil.copytree(TEMPLATE, root)
    (root / 'profiles').mkdir()
    # Linear profiles make the spline/density oracle exact between nodes;
    # composition varies radially and includes unequal ion fractions.
    nodes = [i / 8 for i in range(9)]
    (root / 'equil.dat').write_text(''.join(f'100 {s} 2 {s} {s}\n' for s in nodes))
    profiles = {'n.dat': (NE, 0.2*NE), 'Ti.dat': (TI, .1*TI),
                'Te.dat': (TE, -.1*TE), 'Er.dat': (0., 0.),
                'nD.dat': (fractions[0]*NE, .15*NE),
                'nT.dat': (fractions[1]*NE, .05*NE),
                'ne.dat': (electron*NE, electron*.2*NE)}
    for file, (value, slope) in profiles.items():
        (root / 'profiles' / file).write_text(''.join(
            f'{s:.17e} {value+slope*s:.17e}\n' for s in nodes))
    if explicit:
        (root / 'profiles/collision_species.dat').write_text(
            '3\n' + ''.join(f'{mass:.17e} {z} {file}\n' for mass, z, file in
            zip(MASSES, [1, 1, -1], ['nD.dat', 'nT.dat', 'ne.dat'])))
    text = (root / 'gorilla.inp').read_text()
    text = re.sub(r'(coord_system\s*=)\s*1', r'\1 2', text)
    text = re.sub(r"(er_profile_file\s*=)\s*'[^']*'", r"\1 ''", text)
    (root / 'gorilla.inp').write_text(text)
    text = (root / 'tetra_grid.inp').read_text()
    for key, value in dict(grid_kind=6, n1=5, n2=8, n3=12).items():
        text = re.sub(rf'({key}\s*=)\s*\d+', rf'\g<1> {value}', text)
    text = text.replace('&TETRA_GRID_NML', "&TETRA_GRID_NML\n netcdf_filename='chartmap.nc',")
    (root / 'tetra_grid.inp').write_text(text)
    with (root / 'writer.log').open('w') as log:
        subprocess.run([str(WRITER), 'chartmap.nc', '0', '10', '3', '2',
                        '513', '32', '2'], cwd=root, stdout=log,
                       stderr=subprocess.STDOUT, check=True)
    return root


def run(root, n, expected_error=None):
    with (root / 'run.log').open('w') as log:
        p = subprocess.run([str(DRIVER), str(n)], cwd=root, stdout=log,
                           stderr=subprocess.STDOUT)
    if expected_error:
        assert p.returncode != 0, root
        assert expected_error in (root / 'run.log').read_text(), root
        return
    assert p.returncode == 0, root / 'run.log'
    return [list(map(float, line.split())) for line in (root / 'rates.dat').read_text().splitlines()]


def rate(mt, mb, za, zb, nb, tb, exact_shape=False):
    # Direct Rutherford prefactor times the Maxwellian deflection shape,
    # with the native pair-log choice (T_test=T_background in collis_init).
    v = math.sqrt(2*3500*1.6022e-12/mt)
    x = v / math.sqrt(2*tb*1.6022e-12/mb)
    if za > 0 and zb > 0:
        logl = 23-math.log(za*zb/tb*math.sqrt(nb*(za*za+zb*zb)/tb))
    else:
        zi, mi = (za, mt) if za > 0 else (zb, mb)
        if tb < tb*9.1094e-28/mi:
            logl = 16-math.log(math.sqrt(nb)*tb**-1.5*zi**2*mi/1.6726e-24)
        elif tb < 10*zi**2:
            logl = 23-math.log(math.sqrt(nb)*zi*tb**-1.5)
        else:
            logl = 24-math.log(math.sqrt(nb)/tb)
    shape = math.erf(x)-(math.erf(x)-2*x*math.exp(-x*x)/math.sqrt(math.pi))/(2*x*x)
    if not exact_shape:
        if x < .01:
            shape = .75225278*x*(1-.2*x*x)
        elif x > 6:
            shape = 1-.5/(x*x)
        else:
            shape = math.erf(x)*(1-.5/(x*x))+math.exp(-x*x)/(1.7724538*x)
    prefactor = 4*math.pi*nb*(za*zb)**2*(4.8032e-10)**4*logl/(mt*mt*v**3)
    return prefactor*shape, prefactor*x/(2*v)


def main():
    pure = prepare('pure_D', explicit=False)
    run(pure, 2)  # Fortran checks bitwise against the frozen pure-plasma calculation.
    print('PASS: pure D/e collision coefficients and velocity/energy ratios bitwise unchanged')
    root = prepare('ITER_DT')
    rows = run(root, 3)
    err, coefficient_err, exact_shape_err = 0., 0., 0.
    for a, b, cell, s, mt, mb, za, zb, nb, tb, nu, coefficient in rows:
        expected_n = [(0.5+.15*s)*NE, (0.5+.05*s)*NE, (1+.2*s)*NE][int(b)-1]
        assert abs(nb/expected_n-1) < 2e-14
        assert abs(mb/MASSES[int(b)-1]-1) < 2e-15
        ref, coefficient_ref = rate(mt, mb, za, zb, nb, tb)
        err = max(err, abs(nu/ref-1))
        coefficient_err = max(coefficient_err, abs(coefficient/coefficient_ref-1))
        exact, _ = rate(mt, mb, za, zb, nb, tb, exact_shape=True)
        exact_shape_err = max(exact_shape_err, abs(nu/exact-1))
    assert coefficient_err < 2e-14, coefficient_err
    assert err < 1e-7, err  # Native onseff uses sqrt(pi)=1.7724538.
    print(f'PASS: independent collision-frequency coefficients max relative error {coefficient_err:.9e}, tolerance 2e-14')
    print(f'LIMIT: existing native kernel / exact Maxwellian shape error {exact_shape_err:.9e}')
    print(f'PASS: D/T/e pair deflection rates max relative error {err:.9e}, tolerance 1e-7')
    for a in (1, 2, 3):
        first = [r for r in rows if r[0] == a and r[2] == 1]
        print(f'TEST_SPECIES {a} cell1 s={first[0][3]:.17e} pair_rates_s^-1=' +
              ','.join(f'{r[-2]:.17e}' for r in first))
    # Rejection checks exercise the production loader, not just a separate validator.
    run(prepare('missing_mixture', explicit=False), 3, 'Mixture requires collision_species.dat')
    run(prepare('nonneutral', fractions=(.6, .5)), 3, 'violate quasineutrality')
    run(prepare('wrong_electron', fractions=(.6, .6), electron=1.2), 3, 'differs from n.dat')
    missing = prepare('missing_density')
    (missing / 'profiles/nT.dat').unlink()
    run(missing, 3, 'Cannot open required collision density profile')
    print('PASS: absent mixture/density, nonneutral and mismatched electron profiles rejected')


if __name__ == "__main__":
    main()
