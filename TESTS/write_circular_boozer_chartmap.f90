!> Writes the Boozer chartmap of GORILLA's analytic circular tokamak
!> (grid_kind 5: B_phi = B0 R0/R, B_pol = B0 r/(R q(r)), q = q0 + q1 (r/a)^2)
!> with an optional |B| perturbation, for the two-route tests of the
!> canonical NTV estimator.
!>
!> For this field (r, theta*, phi) are Boozer coordinates: with the straight-
!> field-line angle tan(theta*/2) = sqrt((1-eps)/(1+eps)) tan(theta/2),
!> eps = r/R0, and the geometric phi, the covariant components
!>   B_phi = B0 R0,  B_theta* = B0 r^2/(q R0 sqrt(1-eps^2)),  B_r = 0
!> are flux functions. With s = (R0 - sqrt(R0^2-r^2))/(R0 - sqrt(R0^2-a^2))
!> (normalised toroidal flux), A_theta = torflux s with
!> torflux = -B0 R0 (R0 - sqrt(R0^2-a^2)) (the (s, theta*, phi) chart is
!> left-handed) and A_phi = psi_pol_analytic_circ(r). The perturbation
!>   |B| += sigma [b cos(m theta* - n phi) + b2 sin(m2 theta* - n2 phi)]
!> leaves A and B_k axisymmetric (chartmap model).
!>
!> Usage: write_circular_boozer_chartmap.x file sigma b m n
!>            [nrho ntheta nzeta [m2 n2 b2 [nfp relative R0]]]]
!> Geometry: R0 = 165 cm, a = 69.991288 cm, B0 = -17991.7336 G, q0 = 1, q1 = 2
!> (the rmp_fm_stationarity template).
program write_circular_boozer_chartmap

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use netcdf
    use field_analytic_circ_mod, only: set_field_analytic_circ, &
        psi_pol_analytic_circ, q_analytic_circ

    implicit none

    real(dp), parameter :: a = 69.991288_dp, b0 = -17991.7336_dp
    real(dp) :: r0 = 165.0_dp
    real(dp), parameter :: q0 = 1.0_dp, q1 = 2.0_dp
    real(dp), parameter :: pi = acos(-1.0_dp)
    character(len=512) :: filename, arg
    real(dp) :: sigma, bmn, b2
    integer :: m, n, m2, n2, nr, nt, nz
    integer :: ir, it, iz, ncid, dims(4), vars(12), nfp
    logical :: relative
    real(dp) :: s_edge, torflux, r, eps, th, q
    real(dp), allocatable :: rho(:), s(:), theta(:), zeta(:), aphi(:), btheta(:), &
        bphi(:), x(:, :, :), y(:, :, :), z(:, :, :), bmod(:, :, :)

    if (command_argument_count() < 5) error stop &
        'usage: write_circular_boozer_chartmap.x file sigma b m n [nr nt nz [m2 n2 b2]]'
    call get_command_argument(1, filename)
    sigma = real_arg(2)
    bmn = real_arg(3)
    m = int_arg(4)
    n = int_arg(5)
    nr = 33
    nt = 64
    nz = 48
    m2 = 0
    n2 = 0
    b2 = 0.0_dp
    if (command_argument_count() >= 8) then
        nr = int_arg(6)
        nt = int_arg(7)
        nz = int_arg(8)
    end if
    if (command_argument_count() >= 11) then
        m2 = int_arg(9)
        n2 = int_arg(10)
        b2 = real_arg(11)
    end if

    nfp = 1
    relative = .false.
    if (command_argument_count() >= 12) nfp = int_arg(12)
    if (command_argument_count() >= 13) relative = int_arg(13) /= 0
    if (command_argument_count() >= 14) r0 = real_arg(14)
    if (nfp < 1 .or. r0 <= a) error stop 'invalid analytic chart geometry'
    if (mod(n,nfp) /= 0 .or. mod(n2,nfp) /= 0) &
        error stop 'harmonic is not periodic on the requested chart wedge'
    call set_field_analytic_circ(r0, a, b0, q0, q1, 2)
    s_edge = r0 - sqrt(r0**2 - a**2)
    torflux = -b0*r0*s_edge

    allocate(rho(nr), s(nr), theta(nt), zeta(nz), aphi(nr), btheta(nr), bphi(nr))
    allocate(x(nr, nt, nz), y(nr, nt, nz), z(nr, nt, nz), bmod(nr, nt, nz))
    do ir = 1, nr
        rho(ir) = real(ir - 1, dp)/real(nr - 1, dp)
    end do
    s = rho**2
    do it = 1, nt
        theta(it) = 2.0_dp*pi*real(it - 1, dp)/real(nt, dp)
    end do
    do iz = 1, nz
        zeta(iz) = 2.0_dp*pi*real(iz - 1, dp)/real(nz*nfp, dp)
    end do
    do ir = 1, nr
        r = sqrt(max(0.0_dp, r0**2 - (r0 - s(ir)*s_edge)**2))
        eps = r/r0
        q = q_analytic_circ(r)
        aphi(ir) = psi_pol_analytic_circ(r)
        btheta(ir) = b0*r**2/(q*r0*sqrt(1.0_dp - eps**2))
        bphi(ir) = b0*r0
        do it = 1, nt
            th = 2.0_dp*atan2(sqrt(1.0_dp + eps)*sin(0.5_dp*theta(it)), &
                              sqrt(1.0_dp - eps)*cos(0.5_dp*theta(it)))
            do iz = 1, nz
                x(ir, it, iz) = (r0 + r*cos(th))*cos(zeta(iz))
                y(ir, it, iz) = (r0 + r*cos(th))*sin(zeta(iz))
                z(ir, it, iz) = r*sin(th)
                bmod(ir, it, iz) = abs(b0)/(r0 + r*cos(th))*sqrt(r0**2 + (r/q)**2)
                if (relative) then
                    bmod(ir,it,iz) = bmod(ir,it,iz) &
                        * (1.0_dp+sigma*bmn*cos(real(m,dp)*theta(it)-real(n,dp)*zeta(iz))) &
                        + sigma*b2*sin(real(m2,dp)*theta(it)-real(n2,dp)*zeta(iz))
                else
                    bmod(ir,it,iz) = bmod(ir,it,iz) &
                        + sigma*(bmn*cos(real(m,dp)*theta(it)-real(n,dp)*zeta(iz)) &
                        + b2*sin(real(m2,dp)*theta(it)-real(n2,dp)*zeta(iz)))
                end if
            end do
        end do
    end do

    call nc(nf90_create(trim(filename), nf90_clobber, ncid))
    call nc(nf90_def_dim(ncid, "rho", nr, dims(1)))
    call nc(nf90_def_dim(ncid, "s", nr, dims(2)))
    call nc(nf90_def_dim(ncid, "theta", nt, dims(3)))
    call nc(nf90_def_dim(ncid, "zeta", nz, dims(4)))
    call nc(nf90_def_var(ncid, "rho", nf90_double, dims(1:1), vars(1)))
    call nc(nf90_def_var(ncid, "s", nf90_double, dims(2:2), vars(2)))
    call nc(nf90_def_var(ncid, "theta", nf90_double, dims(3:3), vars(3)))
    call nc(nf90_def_var(ncid, "zeta", nf90_double, dims(4:4), vars(4)))
    call nc(nf90_def_var(ncid, "A_phi", nf90_double, dims(2:2), vars(5)))
    call nc(nf90_def_var(ncid, "B_theta", nf90_double, dims(1:1), vars(6)))
    call nc(nf90_def_var(ncid, "B_phi", nf90_double, dims(1:1), vars(7)))
    call nc(nf90_def_var(ncid, "x", nf90_double, dims([1, 3, 4]), vars(8)))
    call nc(nf90_def_var(ncid, "y", nf90_double, dims([1, 3, 4]), vars(9)))
    call nc(nf90_def_var(ncid, "z", nf90_double, dims([1, 3, 4]), vars(10)))
    call nc(nf90_def_var(ncid, "Bmod", nf90_double, dims([1, 3, 4]), vars(11)))
    call nc(nf90_def_var(ncid, "num_field_periods", nf90_int, vars(12)))
    call nc(nf90_put_att(ncid, nf90_global, "torflux", torflux))
    call nc(nf90_put_att(ncid, nf90_global, "boozer_field", 1))
    call nc(nf90_enddef(ncid))
    call nc(nf90_put_var(ncid, vars(1), rho))
    call nc(nf90_put_var(ncid, vars(2), s))
    call nc(nf90_put_var(ncid, vars(3), theta))
    call nc(nf90_put_var(ncid, vars(4), zeta))
    call nc(nf90_put_var(ncid, vars(5), aphi))
    call nc(nf90_put_var(ncid, vars(6), btheta))
    call nc(nf90_put_var(ncid, vars(7), bphi))
    call nc(nf90_put_var(ncid, vars(8), x))
    call nc(nf90_put_var(ncid, vars(9), y))
    call nc(nf90_put_var(ncid, vars(10), z))
    call nc(nf90_put_var(ncid, vars(11), bmod))
    call nc(nf90_put_var(ncid, vars(12), nfp))
    call nc(nf90_close(ncid))

contains

    real(dp) function real_arg(i)
        integer, intent(in) :: i
        call get_command_argument(i, arg)
        read(arg, *) real_arg
    end function real_arg

    integer function int_arg(i)
        integer, intent(in) :: i
        call get_command_argument(i, arg)
        read(arg, *) int_arg
    end function int_arg

    subroutine nc(status)
        integer, intent(in) :: status
        if (status /= nf90_noerr) error stop nf90_strerror(status)
    end subroutine nc

end program write_circular_boozer_chartmap
