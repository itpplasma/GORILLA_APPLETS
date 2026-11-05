program test_stost_cartesian

! B01 for GORILLA's existing Maxwellian operator (collis_init, coleff).
! stost advances p = |v|/v0 and lambda = v_par/|v| in tau = v0 t by
!   dp      = fpeff dtau + sqrt(2 dpp dtau) xi
!   dlambda = -2 lambda dhh dtau + sqrt(2 dhh (1 - lambda^2) dtau) xi
! Ito calculus for |v| and the pitch angle gives the Cartesian coefficients
!   D_par = v0^3 dpp,  D_perp = v0^3 p^2 dhh,
!   D_ij  = D_par vh_i vh_j + D_perp (delta_ij - vh_i vh_j),
!   A     = (v0^2 fpeff - 2 D_perp/|v|) vh,
! which are compared with the oracle of DERIVATIONS/homma_hatayama.wl.
! collis_init computes its own Coulomb logarithm; efcolf is linear in it and
! is rescaled to the benchmark value 15. GORILLA's constants_mod has 5-digit
! e, eV and m_p and onseff an 8-digit sqrt(pi), hence the loose tolerance.
! Argument: directory with the oracle files (TESTS/REFERENCE/homma).

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use collis_ions, only: collis_init, coleff, lambda_alpha_beta

    implicit none

    real(dp), parameter :: amu_g = 1.66053906660e-24_dp
    real(dp), parameter :: m_p_g = 1.67262192369e-24_dp
    real(dp), parameter :: tol = 1.0e-4_dp

    character(len=1024) :: ref_dir
    character(len=4096) :: line
    real(dp) :: m_test, z_test, mass(1), charge(1), dens(1), temp(1), e0, v0
    real(dp) :: efcolf(1), velrat(1), enrat(1), dpp(1), dhh(1), fpeff(1), lambda_coul
    real(dp) :: z, grad(3), x, dir(3), ref(18), p, v_abs, d_par, d_perp, a_v
    real(dp) :: a(3), d(3,3), err_a, err_d, max_err_a, max_err_d
    integer :: unit, ios, i, nrows
    logical :: all_passed

    if (command_argument_count() < 1) error stop 'usage: test_stost_cartesian.x <oracle directory>'
    call get_command_argument(1, ref_dir)

    ! W3+ (183.84 amu) in H+ at n = 1e20 m^-3, T = 50 eV; CGS inside GORILLA
    m_test = 183.84_dp*amu_g
    z_test = 3.0_dp
    mass = m_p_g
    charge = 1.0_dp
    dens = 1.0e14_dp
    temp = 50.0_dp
    e0 = 50.0_dp
    call collis_init(m_test, z_test, mass, charge, dens, temp, e0, v0, efcolf, velrat, enrat)
    call lambda_alpha_beta(z_test, charge(1), m_test, mass(1), temp(1), temp(1), dens(1), dens(1), lambda_coul)
    efcolf = efcolf*15.0_dp/lambda_coul

    open(newunit=unit, file=trim(ref_dir)//'/B01_oracle.csv', status='old', action='read')
    do
        read(unit, '(a)') line
        if (line(1:1) /= '#') exit
    end do

    all_passed = .true.
    max_err_a = 0.0_dp
    max_err_d = 0.0_dp
    nrows = 0
    do
        read(unit, *, iostat=ios) z, grad, x, dir, ref
        if (ios /= 0) exit
        if (x == 0.0_dp) cycle          ! (p, lambda) is singular at v = 0
        nrows = nrows + 1

        p = x/velrat(1)                  ! x = p velrat
        call coleff(efcolf, velrat, enrat, p, dpp, dhh, fpeff)
        v_abs = p*v0
        d_par = v0**3*dpp(1)
        d_perp = v0**3*p**2*dhh(1)
        a_v = v0**2*fpeff(1) - 2.0_dp*d_perp/v_abs

        a = a_v*dir*1.0e-2_dp            ! cm/s^2 -> m/s^2
        do i = 1, 3
            d(:,i) = (d_par - d_perp)*dir*dir(i)
            d(i,i) = d(i,i) + d_perp
        end do
        d = d*1.0e-4_dp                  ! cm^2/s^3 -> m^2/s^3

        err_a = maxval(abs(a - ref(1:3)))/maxval(abs(ref(1:3)))
        err_d = maxval(abs([d(1,1), d(1,2), d(1,3), d(2,2), d(2,3), d(3,3)] - ref(7:12))) &
              /maxval(abs(ref(7:12)))
        max_err_a = max(max_err_a, err_a)
        max_err_d = max(max_err_d, err_d)
        if (max(err_a, err_d) > tol) then
            print *, 'FAIL: x = ', x, ' dir = ', dir, ' rel. errors A, D = ', err_a, err_d
            all_passed = .false.
        end if
    end do
    close(unit)

    print '(a,i0,a,2es10.2)', 'stost in Cartesian form, B01: ', nrows, ' rows, max rel. error A, D = ', &
        max_err_a, max_err_d
    if (nrows == 0) all_passed = .false.
    if (.not. all_passed) then
        print *, 'FAIL: stost Cartesian test failed'
        error stop 1
    end if
    print *, 'PASS: stost Cartesian test passed'

end program test_stost_cartesian
