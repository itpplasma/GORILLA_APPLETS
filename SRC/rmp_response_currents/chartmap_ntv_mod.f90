!> Boozer chartmap (grid_kind 6) support for the canonical NTV estimator.
!>
!> The chartmap holds an axisymmetric vector potential A_theta = torflux s,
!> A_phi(s), axisymmetric covariant B_theta(rho), B_phi(rho) and a 3D |B| on
!> (rho, theta_B, zeta_B) nodes, rho = sqrt(s). This module reads the raw node
!> data (independently of the copy GORILLA interpolates) and provides
!>   * the flux table psi_pol = A_phi(s), psi_tor = torflux s, q, r_eff on the
!>     s nodes, written in the flux_functions.dat format, from which
!>     profile_data_mod builds the psi* -> s inverse;
!>   * Fourier harmonics b_mn(rho) of |B| for one toroidal mode number n,
!>       |B| = sum_m [bc_mn cos(m theta - n zeta) + bs_mn sin(m theta - n zeta)],
!>     by an exact discrete Fourier projection on the uniform angle nodes;
!>   * the zeta average B0(rho, theta) of |B| on the zeta nodes (the
!>     axisymmetric part for a perturbation without n = 0 components), the
!>     axisymmetric Boozer Jacobian J0 = (B_phi torflux - B_theta A_phi')/B0^2
!>     and the physical volume of an s window for marker loading. These use
!>     the chartmap GORILLA loaded for the mesh (same interpolation).
module chartmap_ntv_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use netcdf

    implicit none

    private
    public :: chartmap_raw_t, read_chartmap_raw, chartmap_harmonic, &
        write_chartmap_flux_functions, chartmap_q_on_s_nodes, chartmap_r_eff, &
        chartmap_mesh_setup, chartmap_point, chartmap_window_volume, &
        chartmap_mesh_loaded

    real(dp), parameter :: twopi = 8.0_dp*atan(1.0_dp)

    !> Raw node data of one chartmap file (unique periodic nodes only).
    type :: chartmap_raw_t
        integer :: nrho = 0, ns = 0, nth = 0, nze = 0, nfp = 1
        real(dp) :: torflux = 0.0_dp
        real(dp), allocatable :: rho(:), s(:), aphi(:), btheta(:), bphi(:)
        real(dp), allocatable :: theta(:), zeta(:)
        real(dp), allocatable :: bmod(:, :, :), r_cyl(:, :, :), z_cyl(:, :, :)
    end type chartmap_raw_t

    !> Raw data of the mesh chartmap (zeta nodes for B0, geometry for r_eff).
    type(chartmap_raw_t), save :: mesh_cm
    logical, save :: mesh_loaded = .false.

contains

!> Reads the mesh chartmap (GORILLA's netcdf_filename) and writes its flux
!> table to flux_file (the equil_mapping_file of the profiles).
subroutine chartmap_mesh_setup(mesh_file, flux_file)
    character(len=*), intent(in) :: mesh_file, flux_file

    call read_chartmap_raw(mesh_file, mesh_cm, with_geometry=.true.)
    mesh_loaded = .true.
    call write_chartmap_flux_functions(mesh_cm, flux_file)
    print *, 'Chartmap flux table written to ', trim(flux_file)
end subroutine chartmap_mesh_setup

logical function chartmap_mesh_loaded()
    chartmap_mesh_loaded = mesh_loaded
end function chartmap_mesh_loaded

!> Mesh-chartmap quantities at (s, theta, zeta): interpolated |B| (b), its
!> zeta average over the zeta nodes (b0), the axisymmetric Jacobian j0 (with
!> b0), the covariant B_phi and A_phi(s).
subroutine chartmap_point(s, theta, zeta, b, b0, j0, bphi_cov, aphi)
    use boozer_chartmap_mod, only: evaluate_boozer_chartmap
    real(dp), intent(in)  :: s, theta, zeta
    real(dp), intent(out) :: b, b0, j0, bphi_cov, aphi
    real(dp) :: a_s, a_t, h_s, h_t, h_p, sqrtg, bk, dummy(7)
    integer :: k

    if (.not. mesh_loaded) error stop 'chartmap_ntv: mesh chartmap not set up'
    call evaluate_boozer_chartmap(s, theta, zeta, a_s, a_t, aphi, h_s, h_t, h_p, b, sqrtg)
    b0 = 0.0_dp
    do k = 1, mesh_cm%nze
        call evaluate_boozer_chartmap(s, theta, mesh_cm%zeta(k), dummy(1), dummy(2), &
            dummy(3), dummy(4), dummy(5), dummy(6), bk, dummy(7))
        b0 = b0 + bk
    end do
    b0 = b0/real(mesh_cm%nze, dp)
    bphi_cov = h_p*b
    j0 = sqrtg*(b/b0)**2
end subroutine chartmap_point

!> Physical volume of the window s1 <= s <= s2 over one of nper toroidal
!> periods with the axisymmetric Jacobian |J0| (midpoint rule, ns_int x nth
!> cells), and the maximum of |J0| on the window (rejection bound).
subroutine chartmap_window_volume(s1, s2, nper, volume, j0_max)
    real(dp), intent(in)  :: s1, s2
    integer, intent(in)   :: nper
    real(dp), intent(out) :: volume, j0_max
    integer, parameter :: ns_int = 200
    integer :: i, k, nth
    real(dp) :: s, th, b, b0, j0, bphi, aphi, ds, dth

    if (.not. mesh_loaded) error stop 'chartmap_ntv: mesh chartmap not set up'
    nth = 2*mesh_cm%nth
    ds = (s2 - s1)/real(ns_int, dp)
    dth = twopi/real(nth, dp)
    volume = 0.0_dp
    j0_max = 0.0_dp
    do i = 1, ns_int
        s = s1 + (real(i, dp) - 0.5_dp)*ds
        do k = 1, nth
            th = (real(k, dp) - 0.5_dp)*dth
            call chartmap_point(s, th, 0.0_dp, b, b0, j0, bphi, aphi)
            volume = volume + abs(j0)
            j0_max = max(j0_max, abs(j0))
        end do
        call chartmap_point(s1 + (i - 1)*ds, 0.0_dp, 0.0_dp, b, b0, j0, bphi, aphi)
        j0_max = max(j0_max, abs(j0))
    end do
    volume = volume*ds*dth*twopi/real(nper, dp)
    j0_max = 1.02_dp*j0_max
end subroutine chartmap_window_volume

!> Reads a Boozer chartmap. A duplicated periodic endpoint is dropped.
subroutine read_chartmap_raw(filename, cm, with_geometry)
    character(len=*), intent(in) :: filename
    type(chartmap_raw_t), intent(out) :: cm
    logical, intent(in), optional :: with_geometry

    integer :: ncid, nth_file, nze_file, booz, status
    real(dp), allocatable :: theta(:), zeta(:), b3(:, :, :), x3(:, :, :), y3(:, :, :)
    logical :: geom

    geom = .true.
    if (present(with_geometry)) geom = with_geometry
    call nc(nf90_open(trim(filename), nf90_nowrite, ncid), 'open '//trim(filename))
    cm%nrho = dim_len(ncid, 'rho')
    cm%ns = dim_len(ncid, 's')
    nth_file = dim_len(ncid, 'theta')
    nze_file = dim_len(ncid, 'zeta')
    allocate(cm%rho(cm%nrho), cm%s(cm%ns), cm%aphi(cm%ns), cm%btheta(cm%nrho), &
             cm%bphi(cm%nrho), theta(nth_file), zeta(nze_file))
    allocate(b3(cm%nrho, nth_file, nze_file))
    call get1(ncid, 'rho', cm%rho)
    call get1(ncid, 's', cm%s)
    call get1(ncid, 'A_phi', cm%aphi)
    call get1(ncid, 'B_theta', cm%btheta)
    call get1(ncid, 'B_phi', cm%bphi)
    call get1(ncid, 'theta', theta)
    call get1(ncid, 'zeta', zeta)
    call get3(ncid, 'Bmod', b3)
    call nc(nf90_get_var(ncid, var_id(ncid, 'num_field_periods'), cm%nfp), &
            'read num_field_periods')
    call nc(nf90_get_att(ncid, nf90_global, 'torflux', cm%torflux), 'read torflux')
    booz = 0
    status = nf90_get_att(ncid, nf90_global, 'boozer_field', booz)
    if (status /= nf90_noerr .or. booz /= 1) &
        error stop 'chartmap_ntv: file is not an extended Boozer chartmap'
    if (cm%nfp < 1) error stop 'chartmap_ntv: invalid field-period count'
    cm%nth = unique_count(theta, twopi)
    cm%nze = unique_count(zeta, twopi/real(cm%nfp, dp))
    allocate(cm%theta(cm%nth), cm%zeta(cm%nze))
    cm%theta = theta(:cm%nth)
    cm%zeta = zeta(:cm%nze)
    allocate(cm%bmod(cm%nrho, cm%nth, cm%nze))
    cm%bmod = b3(:, :cm%nth, :cm%nze)
    if (geom) then
        allocate(x3(cm%nrho, nth_file, nze_file), y3(cm%nrho, nth_file, nze_file))
        call get3(ncid, 'x', x3)
        call get3(ncid, 'y', y3)
        call get3(ncid, 'z', b3)
        allocate(cm%r_cyl(cm%nrho, cm%nth, cm%nze), cm%z_cyl(cm%nrho, cm%nth, cm%nze))
        cm%r_cyl = sqrt(x3(:, :cm%nth, :cm%nze)**2 + y3(:, :cm%nth, :cm%nze)**2)
        cm%z_cyl = b3(:, :cm%nth, :cm%nze)
    end if
    call nc(nf90_close(ncid), 'close '//trim(filename))
    if (any(cm%s(2:) <= cm%s(:cm%ns - 1))) error stop 'chartmap_ntv: s must increase'
    if (any(cm%rho(2:) <= cm%rho(:cm%nrho - 1))) &
        error stop 'chartmap_ntv: rho must increase'
end subroutine read_chartmap_raw

!> Cosine and sine amplitudes of the (m, n) harmonic of |B| on the rho nodes.
!> n refers to the full torus (zeta period 2 pi/nfp); a mode with n not a
!> multiple of nfp is absent from the file and returns zero.
subroutine chartmap_harmonic(cm, m, n, bc, bs)
    type(chartmap_raw_t), intent(in) :: cm
    integer, intent(in) :: m, n
    real(dp), intent(out) :: bc(:), bs(:)

    integer :: it, iz
    real(dp) :: phase, norm

    bc = 0.0_dp
    bs = 0.0_dp
    if (mod(n, cm%nfp) /= 0) return
    if (2*abs(m) >= cm%nth .or. 2*abs(n/cm%nfp) >= cm%nze) &
        error stop 'chartmap_ntv: harmonic above the Nyquist limit of the chartmap'
    do iz = 1, cm%nze
        do it = 1, cm%nth
            phase = real(m, dp)*cm%theta(it) - real(n, dp)*cm%zeta(iz)
            bc = bc + cm%bmod(:, it, iz)*cos(phase)
            bs = bs + cm%bmod(:, it, iz)*sin(phase)
        end do
    end do
    norm = 2.0_dp/real(cm%nth*cm%nze, dp)
    if (m == 0 .and. n == 0) norm = 0.5_dp*norm
    bc = norm*bc
    bs = norm*bs
end subroutine chartmap_harmonic

!> Safety factor q = 1/iota = -torflux/(dA_phi/ds) on the s nodes
!> (centred differences of the piecewise-linear A_phi, one-sided at the ends).
subroutine chartmap_q_on_s_nodes(cm, q)
    type(chartmap_raw_t), intent(in) :: cm
    real(dp), intent(out) :: q(:)
    integer :: i, lo, hi

    do i = 1, cm%ns
        lo = max(1, i - 1)
        hi = min(cm%ns, i + 1)
        q(i) = -cm%torflux*(cm%s(hi) - cm%s(lo))/(cm%aphi(hi) - cm%aphi(lo))
    end do
end subroutine chartmap_q_on_s_nodes

!> Mean distance of the zeta = zeta(1) cross-section of each rho node from the
!> axis (the rho = rho(1) ring centre), and the outboard major radius R_beg.
subroutine chartmap_r_eff(cm, r_eff, r_beg)
    type(chartmap_raw_t), intent(in) :: cm
    real(dp), intent(out) :: r_eff(:), r_beg(:)
    real(dp) :: r_axis, z_axis
    integer :: i, it0

    if (.not. allocated(cm%r_cyl)) error stop 'chartmap_ntv: geometry not loaded'
    r_axis = sum(cm%r_cyl(1, :, 1))/real(cm%nth, dp)
    z_axis = sum(cm%z_cyl(1, :, 1))/real(cm%nth, dp)
    it0 = minloc(abs(modulo(cm%theta + 0.5_dp*twopi, twopi) - 0.5_dp*twopi), 1)
    do i = 1, cm%nrho
        r_eff(i) = sum(sqrt((cm%r_cyl(i, :, 1) - r_axis)**2 &
                            + (cm%z_cyl(i, :, 1) - z_axis)**2))/real(cm%nth, dp)
        r_beg(i) = cm%r_cyl(i, it0, 1)
    end do
end subroutine chartmap_r_eff

!> Writes the flux table (R_beg [cm], r_eff [cm], q, psi_pol, psi_tor) on the
!> s nodes in the flux_functions.dat format read by profile_data_mod.
!> psi_pol is the chartmap A_phi(s) (the A_phi of the mesh and of f0),
!> psi_tor = torflux s; geometric columns are interpolated in rho.
subroutine write_chartmap_flux_functions(cm, filename)
    type(chartmap_raw_t), intent(in) :: cm
    character(len=*), intent(in) :: filename
    real(dp) :: q(cm%ns), r_eff(cm%nrho), r_beg(cm%nrho)
    integer :: u, i

    call chartmap_q_on_s_nodes(cm, q)
    call chartmap_r_eff(cm, r_eff, r_beg)
    open(newunit=u, file=trim(filename), status='replace', action='write')
    write(u, '(a)') ' # R_beg [cm]  r_eff [cm]  q  psi_pol = A_phi  psi_tor = torflux s' &
        //'  (Boozer chartmap)'
    do i = 1, cm%ns
        write(u, '(5es24.16)') lin(cm%rho, r_beg, sqrt(cm%s(i))), &
            lin(cm%rho, r_eff, sqrt(cm%s(i))), q(i), cm%aphi(i), cm%torflux*cm%s(i)
    end do
    close(u)
end subroutine write_chartmap_flux_functions

pure function lin(x, y, xv) result(v)
    real(dp), intent(in) :: x(:), y(:), xv
    real(dp) :: v, u, xc
    integer :: i

    xc = min(x(size(x)), max(x(1), xv))
    i = 1
    do while (i < size(x) - 1 .and. x(i + 1) < xc)
        i = i + 1
    end do
    u = (xc - x(i))/(x(i + 1) - x(i))
    v = (1.0_dp - u)*y(i) + u*y(i + 1)
end function lin

integer function unique_count(axis, period) result(n)
    real(dp), intent(in) :: axis(:), period
    n = size(axis)
    if (n > 1) then
        if (abs(axis(n) - axis(1) - period) <= 1.0e-10_dp*period) n = n - 1
    end if
end function unique_count

integer function dim_len(ncid, name) result(n)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    integer :: id
    call nc(nf90_inq_dimid(ncid, name, id), 'dimension '//name)
    call nc(nf90_inquire_dimension(ncid, id, len=n), 'dimension '//name)
end function dim_len

integer function var_id(ncid, name) result(id)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    call nc(nf90_inq_varid(ncid, name, id), 'variable '//name)
end function var_id

subroutine get1(ncid, name, v)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    real(dp), intent(out) :: v(:)
    call nc(nf90_get_var(ncid, var_id(ncid, name), v), 'read '//name)
end subroutine get1

subroutine get3(ncid, name, v)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: name
    real(dp), intent(out) :: v(:, :, :)
    call nc(nf90_get_var(ncid, var_id(ncid, name), v), 'read '//name)
end subroutine get3

subroutine nc(status, what)
    integer, intent(in) :: status
    character(len=*), intent(in) :: what
    if (status /= nf90_noerr) then
        print *, 'chartmap_ntv: ', trim(what), ': ', trim(nf90_strerror(status))
        error stop 'chartmap_ntv: NetCDF error'
    end if
end subroutine nc

end module chartmap_ntv_mod
