"""A half-density D component halves canonical P and torque, leaving orbits fixed."""
import math

from test_canonical_ntv import WORK_DIR, prepare, read_markers, run, set_nml


def main():
    markers = []
    for name, mixture in [('pure', False), ('mixture', True)]:
        case = WORK_DIR / name
        prepare(case, 16, '2.0d-6', 1.0, False)
        set_nml(case / 'gorilla.inp', {'ispecies': 2})
        if mixture:
            profiles = case / 'profiles'
            source = (profiles / 'n.dat').read_text().splitlines()
            half = ''.join(f'{float(r):.17e} {float(n)*.5:.17e}\n'
                           for r, n in (row.split() for row in source))
            (profiles / 'nD.dat').write_text(half)
            (profiles / 'nT.dat').write_text(half)
            (profiles / 'ne.dat').write_text('\n'.join(source)+'\n')
            (profiles / 'collision_species.dat').write_text(
                '3\n3.3452e-24 1 nD.dat\n5.0178e-24 1 nT.dat\n'
                '9.1094e-28 -1 ne.dat\n')
            inp = case / 'rmp_response_currents.inp'
            inp.write_text(inp.read_text().replace('&rmp_response_currents_nml',
                '&rmp_response_currents_nml\n n_background_species=3,'))
        run(case, 1)
        markers.append(read_markers(case / 'canonical_ntv_markers.dat'))
    pure, mixed = markers
    for old, new in zip(pure['P'], mixed['P']):
        assert math.isclose(new, .5*old, rel_tol=2e-13)
    for key in ['w_end', 'tq_deltaf', 'dlnf0_dpphi_birth', 'pphi_end', 's_end']:
        for old, new in zip(pure[key], mixed[key]):
            assert math.isclose(new, old, rel_tol=2e-11, abs_tol=1e-25), key
    old_t = sum(p*t for p, t in zip(pure['P'], pure['tq_deltaf']))
    new_t = sum(p*t for p, t in zip(mixed['P'], mixed['tq_deltaf']))
    assert old_t != 0
    assert math.isclose(new_t, .5*old_t, rel_tol=2e-11)
    print('PASS: half-density canonical markers give half torque and identical orbits')


if __name__ == '__main__':
    main()
