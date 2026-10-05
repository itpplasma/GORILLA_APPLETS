!Write a manufactured Boozer chartmap with constant rotational transform.
!Circular cross-section, R0 = 200 cm, minor radius 10 cm, Bmod = 1000 G,
!two field periods. A_phi is linear in s, so iota = -dA_phi/ds/torflux is
!exactly fixture_iota on every surface. Usage: write_chartmap_fixture.x <file>
program write_chartmap_fixture

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use netcdf

    implicit none

    real(dp), parameter :: pi = 4.0_dp*atan(1.0_dp)
    real(dp), parameter :: fixture_iota = 0.4_dp
    real(dp), parameter :: torflux = 2.0e8_dp
    real(dp), parameter :: major_radius = 200.0_dp, minor_radius = 10.0_dp
    integer, parameter :: nrho = 9, ntheta = 16, nzeta = 4, nfp = 2

    character(len=1024) :: filename
    integer :: ncid, dims(4), vars(12), ir, it, iz
    real(dp) :: rho(nrho), s(nrho), theta(ntheta), zeta(nzeta), radius
    real(dp), dimension(nrho,ntheta,nzeta) :: x, y, z, bmod
    real(dp) :: aphi(nrho), btheta(nrho), bphi(nrho)

    if (command_argument_count() /= 1) error stop 'usage: write_chartmap_fixture.x <file>'
    call get_command_argument(1, filename)

    do ir = 1, nrho
        rho(ir) = real(ir - 1, dp)/real(nrho - 1, dp)
    end do
    s = rho**2
    do it = 1, ntheta
        theta(it) = 2.0_dp*pi*real(it - 1, dp)/real(ntheta, dp)
    end do
    do iz = 1, nzeta
        zeta(iz) = 2.0_dp*pi/real(nfp, dp)*real(iz - 1, dp)/real(nzeta, dp)
    end do

    bmod = 1000.0_dp
    aphi = -fixture_iota*torflux*s
    btheta = 100.0_dp + 10.0_dp*rho
    bphi = 500.0_dp - 20.0_dp*rho
    do iz = 1, nzeta
        do it = 1, ntheta
            do ir = 1, nrho
                radius = major_radius + minor_radius*rho(ir)*cos(theta(it))
                x(ir, it, iz) = radius*cos(zeta(iz))
                y(ir, it, iz) = radius*sin(zeta(iz))
                z(ir, it, iz) = minor_radius*rho(ir)*sin(theta(it))
            end do
        end do
    end do

    call ok(nf90_create(trim(filename), nf90_clobber, ncid))
    call ok(nf90_def_dim(ncid, 'rho', nrho, dims(1)))
    call ok(nf90_def_dim(ncid, 's', nrho, dims(2)))
    call ok(nf90_def_dim(ncid, 'theta', ntheta, dims(3)))
    call ok(nf90_def_dim(ncid, 'zeta', nzeta, dims(4)))
    call ok(nf90_def_var(ncid, 'rho', nf90_double, [dims(1)], vars(1)))
    call ok(nf90_def_var(ncid, 's', nf90_double, [dims(2)], vars(2)))
    call ok(nf90_def_var(ncid, 'theta', nf90_double, [dims(3)], vars(3)))
    call ok(nf90_def_var(ncid, 'zeta', nf90_double, [dims(4)], vars(4)))
    call ok(nf90_def_var(ncid, 'A_phi', nf90_double, [dims(2)], vars(5)))
    call ok(nf90_def_var(ncid, 'B_theta', nf90_double, [dims(1)], vars(6)))
    call ok(nf90_def_var(ncid, 'B_phi', nf90_double, [dims(1)], vars(7)))
    call ok(nf90_def_var(ncid, 'x', nf90_double, dims([1, 3, 4]), vars(8)))
    call ok(nf90_def_var(ncid, 'y', nf90_double, dims([1, 3, 4]), vars(9)))
    call ok(nf90_def_var(ncid, 'z', nf90_double, dims([1, 3, 4]), vars(10)))
    call ok(nf90_def_var(ncid, 'Bmod', nf90_double, dims([1, 3, 4]), vars(11)))
    call ok(nf90_def_var(ncid, 'num_field_periods', nf90_int, vars(12)))
    call ok(nf90_put_att(ncid, nf90_global, 'torflux', torflux))
    call ok(nf90_put_att(ncid, nf90_global, 'boozer_field', 1))
    call ok(nf90_enddef(ncid))
    call ok(nf90_put_var(ncid, vars(1), rho))
    call ok(nf90_put_var(ncid, vars(2), s))
    call ok(nf90_put_var(ncid, vars(3), theta))
    call ok(nf90_put_var(ncid, vars(4), zeta))
    call ok(nf90_put_var(ncid, vars(5), aphi))
    call ok(nf90_put_var(ncid, vars(6), btheta))
    call ok(nf90_put_var(ncid, vars(7), bphi))
    call ok(nf90_put_var(ncid, vars(8), x))
    call ok(nf90_put_var(ncid, vars(9), y))
    call ok(nf90_put_var(ncid, vars(10), z))
    call ok(nf90_put_var(ncid, vars(11), bmod))
    call ok(nf90_put_var(ncid, vars(12), nfp))
    call ok(nf90_close(ncid))

contains

    subroutine ok(status)
        integer, intent(in) :: status
        if (status /= nf90_noerr) error stop nf90_strerror(status)
    end subroutine ok

end program write_chartmap_fixture
