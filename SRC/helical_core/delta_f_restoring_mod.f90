module delta_f_restoring_mod
!
! Momentum-, energy- and particle-restoring term for time-dependent delta-f with test-particle
! collisions (field-particle part of the linearised like-particle operator), C(f0, f1) ~ P f_M.
!
! The test-particle kicks (stost) remove momentum and energy from delta-f. Per conservation bin
! (n_psi rings in psi_N x n_theta poloidal sectors, toroidally summed) the losses of one time slice
!     dp = sum_kicks W m dv_par,   dE = sum_kicks W d(m v^2/2),   W = base * w (delta-f weight)
! are returned during the next slice through the weights of all markers in the bin:
!     dw = (1 - w1) P(v) dt,   P = c_1 s_1 + c_2 s_2 + c_3 s_3,
! with the shapes of the test-particle operator (Brunner et al. 1999, Eqs. 9-12; Lin et al. 1995):
!     s_1 = H(x) v_par/v_th   (momentum: drag nu_s(v) v_par),
!     s_2 = K(x) + H(x) - x^2 H(x)   (energy: rate of change of v^2 of a test particle),
!     s_3 = 1   (density),   x = v/v_th, v_th^2 = T/m.
! The coefficients are fixed by exact discrete conservation with the markers present in the bin at the
! end of the slice ("matrix method"): sum_k M_lk c_k = R_l,
!     M_lk = sum_j base_j (1 - w1_j) s_k(v_j) mu_l(v_j),  mu = (1, m v_par, m v^2/2),
!     R = (0, -dp, -dE) / dt_slice.
! Fallback when the bin has fewer than min_markers markers or the scaled 3x3 matrix is ill-conditioned:
! the analytic coefficients of Wang et al. 1999 / Lin et al. (continuum Maxwellian moments of the shapes,
! momentum and energy only): c_1 = -dp/dt / (N_b <s_1 m v_par>_M), c_2 = -dE/dt / (N_b <s_2 m v^2/2>_M),
! c_3 = 0, N_b = sum_j base_j (1 - w1_j). Empty bins are skipped (their losses are not returned).
!
! TODO (electrons): the kicks contain all background species; only the like-species part must be
! restored (the electron-ion part goes to the ions). For ions the ion-electron part is O(sqrt(me/mi)).
!
    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none
    private

    logical, public :: boole_restoring = .false.
    integer, public :: n_rest_psi = 30, n_rest_theta = 16, rest_min_markers = 50
    real(dp), public :: rest_max_condition = 1.0d8

    integer, public :: n_rest_bins = 0
    integer, dimension(:), allocatable, public :: tetr_bin          !conservation bin of every tetrahedron
    integer, dimension(:), allocatable, public :: marker_tetr       !tetrahedron of every marker at the end of a slice
    real(dp), dimension(:,:), allocatable, public :: rest_loss      !(2, bins): momentum and energy lost in the slice
    real(dp), dimension(:,:), allocatable, public :: rest_coef      !(3, bins): coefficients c_k for the next slice
    integer, public :: n_rest_exact = 0, n_rest_fallback = 0, n_rest_empty = 0
    !momentum and energy actually returned to delta-f by the restoring term in the current slice
    !(per thread, summed into rest_given at the end of the slice): check of the conservation
    real(dp), dimension(2), public :: given_thr = 0.0_dp, rest_given = 0.0_dp
    !$omp threadprivate(given_thr)

    real(dp) :: maxwell_hx2 = 0.0_dp, maxwell_s2x2 = 0.0_dp      !<H x^2>_M and <s_2 x^2>_M

    public :: setup_restoring, restoring_shapes, restoring_rate, compute_restoring_coefficients, &
              add_kick_to_losses, write_restoring_diagnostics

contains

! ====================================================================
subroutine setup_restoring(n_markers)
!
! Conservation bins of all tetrahedra (psi_N of the centroid, poloidal angle about the magnetic axis),
! per-marker tetrahedron array, and the Maxwellian averages for the fallback.
!
    use tetra_grid_mod, only: tetra_grid, verts_rphiz, ntetr
    use tetra_physics_mod, only: tetra_physics
    use utils_data_pre_and_post_processing_mod, only: normalised_poloidal_flux
    use constants, only: pi

    integer, intent(in) :: n_markers
    integer :: k, ip, it, i_axis, i
    real(dp), dimension(:), allocatable :: psin, rc, zc
    real(dp), dimension(3) :: centroid
    real(dp) :: theta, x, dx, w, h, kk

    n_rest_bins = n_rest_psi*n_rest_theta
    allocate(tetr_bin(ntetr), psin(ntetr), rc(ntetr), zc(ntetr))
    do k = 1, ntetr
        centroid = 0.25_dp*sum(verts_rphiz(:, tetra_grid(k)%ind_knot(1:4)), dim=2)
        psin(k) = normalised_poloidal_flux(tetra_physics(k)%Aphi1 + &
                                           sum(tetra_physics(k)%gAphi*(centroid - tetra_physics(k)%x1)))
        rc(k) = centroid(1)
        zc(k) = centroid(3)
    enddo
    i_axis = minloc(psin, dim=1)
    do k = 1, ntetr
        ip = min(max(int(psin(k)*n_rest_psi) + 1, 1), n_rest_psi)
        theta = atan2(zc(k) - zc(i_axis), rc(k) - rc(i_axis))
        it = min(max(int((theta + pi)/(2.0_dp*pi)*n_rest_theta) + 1, 1), n_rest_theta)
        tetr_bin(k) = (ip - 1)*n_rest_theta + it
    enddo
    deallocate(psin, rc, zc)

    allocate(marker_tetr(n_markers), rest_loss(2, n_rest_bins), rest_coef(3, n_rest_bins))
    marker_tetr = -1
    rest_loss = 0.0_dp
    rest_coef = 0.0_dp

    !Maxwellian averages over p(x) = sqrt(2/pi) x^2 exp(-x^2/2) (midpoint rule up to x = 12)
    maxwell_hx2 = 0.0_dp
    maxwell_s2x2 = 0.0_dp
    dx = 12.0_dp/20000.0_dp
    do i = 1, 20000
        x = (i - 0.5_dp)*dx
        w = sqrt(2.0_dp/pi)*x**2*exp(-0.5_dp*x**2)*dx
        call rosenbluth_h_k(x, h, kk)
        maxwell_hx2 = maxwell_hx2 + w*h*x**2
        maxwell_s2x2 = maxwell_s2x2 + w*(kk + h - x**2*h)*x**2
    enddo

    print*, 'delta-f restoring: ', n_rest_psi, ' x ', n_rest_theta, ' conservation bins; <H x^2>_M = ', maxwell_hx2, &
            ', <s_2 x^2>_M = ', maxwell_s2x2

end subroutine setup_restoring

! ====================================================================
pure subroutine rosenbluth_h_k(x, h, k)
!
! H(x), K(x) of Brunner et al. 1999, Eqs. (11)-(12), x = v/v_th, v_th^2 = T/m, with series for small x.
!
    real(dp), intent(in) :: x
    real(dp), intent(out) :: h, k
    real(dp), parameter :: c = 0.7978845608028654_dp   !sqrt(2/pi)
    real(dp) :: e, g

    if (x.lt.1.0d-2) then
        h = c*(1.0_dp/3.0_dp - 0.1_dp*x**2)
        k = c*(2.0_dp/3.0_dp - x**2/15.0_dp)
    else
        e = erf(x/sqrt(2.0_dp))
        g = c*x*exp(-0.5_dp*x**2)
        h = (e - g)/x**3
        k = ((x**2 - 1.0_dp)*e + g)/x**3
    endif

end subroutine rosenbluth_h_k

! ====================================================================
pure function restoring_shapes(vpar, vperp, temp_eV, mass) result(s)
!
! Shapes s_1..s_3 at velocity (vpar, vperp) for a Maxwellian of temperature temp_eV, particle mass [g].
!
    use constants, only: ev2erg

    real(dp), intent(in) :: vpar, vperp, temp_eV, mass
    real(dp), dimension(3) :: s
    real(dp) :: vth, x, h, k

    vth = sqrt(temp_eV*ev2erg/mass)
    x = sqrt(vpar**2 + vperp**2)/vth
    call rosenbluth_h_k(x, h, k)
    s(1) = h*vpar/vth
    s(2) = k + h - x**2*h
    s(3) = 1.0_dp

end function restoring_shapes

! ====================================================================
function restoring_rate(ind_tetr, x_rphiz, vpar, vperp, species) result(p_rate)
!
! P(v) [1/s] of the bin of tetrahedron ind_tetr at the guiding-center position x_rphiz (global R, phi, Z).
!
    use gorilla_applets_types_mod, only: start
    use vertex_data_mod, only: vertex_data_profiles

    integer, intent(in) :: ind_tetr, species
    real(dp), dimension(3), intent(in) :: x_rphiz
    real(dp), intent(in) :: vpar, vperp
    real(dp) :: p_rate
    real(dp) :: dens, temp_i, temp_e, temp
    integer :: b

    p_rate = 0.0_dp
    if (.not.boole_restoring) return
    if (ind_tetr.le.0) return
    b = tetr_bin(ind_tetr)
    if (all(rest_coef(:, b).eq.0.0_dp)) return
    call vertex_data_profiles(ind_tetr, x_rphiz(1), x_rphiz(3), dens, temp_i, temp_e)
    temp = merge(temp_i, temp_e, start%particle_charge(species) > 0.0_dp)
    p_rate = sum(rest_coef(:, b)*restoring_shapes(vpar, vperp, temp, start%particle_mass(species)))

end function restoring_rate

! ====================================================================
pure subroutine add_kick_to_losses(loss, ind_tetr, weight, mass, vpar0, vperp0, vpar1, vperp1)
!
! Adds the momentum and energy change of one collision kick of a marker with delta-f weight `weight`
! (physical particles) to the thread-local loss array of its bin.
!
    real(dp), dimension(:,:), intent(inout) :: loss
    integer, intent(in) :: ind_tetr
    real(dp), intent(in) :: weight, mass, vpar0, vperp0, vpar1, vperp1
    integer :: b

    if (ind_tetr.le.0) return
    b = tetr_bin(ind_tetr)
    loss(1, b) = loss(1, b) + weight*mass*(vpar1 - vpar0)
    loss(2, b) = loss(2, b) + weight*0.5_dp*mass*((vpar1**2 + vperp1**2) - (vpar0**2 + vperp0**2))

end subroutine add_kick_to_losses

! ====================================================================
subroutine compute_restoring_coefficients(dt_slice, species, base, w1)
!
! Coefficients c_k of every bin from the losses of the last slice (rest_loss, then reset).
! base(n): physical particles per marker, w1(n): second weight (0 without it).
!
    use gorilla_applets_types_mod, only: start
    use vertex_data_mod, only: vertex_data_profiles
    use constants, only: ev2erg

    real(dp), intent(in) :: dt_slice
    integer, intent(in) :: species
    real(dp), dimension(:), intent(in) :: base, w1

    real(dp), dimension(:,:,:), allocatable :: m_mat
    real(dp), dimension(:), allocatable :: n_b, t_b
    integer, dimension(:), allocatable :: cnt
    real(dp), dimension(3) :: s, mu, rhs, c
    real(dp) :: mass, vpar, vperp, v, dens, temp_i, temp_e, temp, a, vth
    integer :: n, b, k, l, ind
    logical :: ok

    mass = start%particle_mass(species)
    allocate(m_mat(3, 3, n_rest_bins), n_b(n_rest_bins), t_b(n_rest_bins), cnt(n_rest_bins))
    m_mat = 0.0_dp
    n_b = 0.0_dp
    t_b = 0.0_dp
    cnt = 0

    do n = 1, size(marker_tetr)
        ind = marker_tetr(n)
        if (ind.le.0) cycle
        if (start%lost(n, species)) cycle
        b = tetr_bin(ind)
        v = sqrt(2.0_dp*start%energy(n, species)*ev2erg/mass)
        vpar = start%pitch(n, species)*v
        vperp = sqrt(max(v**2 - vpar**2, 0.0_dp))
        call vertex_data_profiles(ind, start%x(1, n, species), start%x(3, n, species), dens, temp_i, temp_e)
        temp = merge(temp_i, temp_e, start%particle_charge(species) > 0.0_dp)
        a = base(n)*(1.0_dp - w1(n))
        s = restoring_shapes(vpar, vperp, temp, mass)
        mu = [1.0_dp, mass*vpar, 0.5_dp*mass*v**2]
        do l = 1, 3
            do k = 1, 3
                m_mat(l, k, b) = m_mat(l, k, b) + a*s(k)*mu(l)
            enddo
        enddo
        n_b(b) = n_b(b) + a
        t_b(b) = t_b(b) + a*temp
        cnt(b) = cnt(b) + 1
    enddo

    n_rest_exact = 0
    n_rest_fallback = 0
    n_rest_empty = 0
    do b = 1, n_rest_bins
        rhs = [0.0_dp, -rest_loss(1, b), -rest_loss(2, b)]/dt_slice
        if ((cnt(b).eq.0).or.(n_b(b).le.0.0_dp)) then
            rest_coef(:, b) = 0.0_dp
            if (any(rest_loss(:, b).ne.0.0_dp)) n_rest_empty = n_rest_empty + 1
            cycle
        endif
        ok = .false.
        if (cnt(b).ge.rest_min_markers) call solve_scaled_3x3(m_mat(:, :, b), rhs, rest_max_condition, c, ok)
        if (ok) then
            n_rest_exact = n_rest_exact + 1
        else
            !Wang/Lin: analytic Maxwellian moments of the shapes, momentum and energy only
            temp = t_b(b)/n_b(b)
            vth = sqrt(temp*ev2erg/mass)
            c(1) = rhs(2)/(n_b(b)*mass*vth*maxwell_hx2/3.0_dp)
            c(2) = rhs(3)/(n_b(b)*0.5_dp*mass*vth**2*maxwell_s2x2)
            c(3) = 0.0_dp
            n_rest_fallback = n_rest_fallback + 1
        endif
        rest_coef(:, b) = c
    enddo
    rest_loss = 0.0_dp

end subroutine compute_restoring_coefficients

! ====================================================================
pure subroutine solve_scaled_3x3(a, rhs, max_cond, c, ok)
!
! Solves a c = rhs after row and column equilibration (the rows have different units: number,
! momentum, energy); ok = .false. if a row or column vanishes or the condition number (infinity norm)
! of the equilibrated matrix exceeds max_cond.
!
    real(dp), dimension(3,3), intent(in) :: a
    real(dp), dimension(3), intent(in) :: rhs
    real(dp), intent(in) :: max_cond
    real(dp), dimension(3), intent(out) :: c
    logical, intent(out) :: ok

    real(dp), dimension(3,3) :: as, inv
    real(dp), dimension(3) :: dr, dc, y
    real(dp) :: det, norm_a, norm_inv
    integer :: i, j

    ok = .false.
    c = 0.0_dp
    do i = 1, 3
        dr(i) = maxval(abs(a(i, :)))
        if (.not.(dr(i).gt.0.0_dp)) return
    enddo
    do j = 1, 3
        dc(j) = maxval(abs(a(:, j))/dr)
        if (.not.(dc(j).gt.0.0_dp)) return
    enddo
    do i = 1, 3
        do j = 1, 3
            as(i, j) = a(i, j)/(dr(i)*dc(j))
        enddo
    enddo
    det = as(1,1)*(as(2,2)*as(3,3) - as(2,3)*as(3,2)) - as(1,2)*(as(2,1)*as(3,3) - as(2,3)*as(3,1)) &
        + as(1,3)*(as(2,1)*as(3,2) - as(2,2)*as(3,1))
    if (det.eq.0.0_dp) return
    inv(1,1) = as(2,2)*as(3,3) - as(2,3)*as(3,2)
    inv(1,2) = as(1,3)*as(3,2) - as(1,2)*as(3,3)
    inv(1,3) = as(1,2)*as(2,3) - as(1,3)*as(2,2)
    inv(2,1) = as(2,3)*as(3,1) - as(2,1)*as(3,3)
    inv(2,2) = as(1,1)*as(3,3) - as(1,3)*as(3,1)
    inv(2,3) = as(1,3)*as(2,1) - as(1,1)*as(2,3)
    inv(3,1) = as(2,1)*as(3,2) - as(2,2)*as(3,1)
    inv(3,2) = as(1,2)*as(3,1) - as(1,1)*as(3,2)
    inv(3,3) = as(1,1)*as(2,2) - as(1,2)*as(2,1)
    inv = inv/det
    norm_a = maxval(sum(abs(as), dim=2))
    norm_inv = maxval(sum(abs(inv), dim=2))
    if (norm_a*norm_inv.gt.max_cond) return
    !a = Dr as Dc  =>  as (Dc c) = Dr^-1 rhs
    y = matmul(inv, rhs/dr)
    c = y/dc
    ok = .true.

end subroutine solve_scaled_3x3

! ====================================================================
subroutine write_restoring_diagnostics(unit, i_slice, time, dp_tot, de_tot, dp_abs)
!
! One line per slice: total momentum and energy lost by delta-f in the slice (sum over bins), the sum of
! |dp| over bins, momentum and energy returned by the restoring term in the slice (from the coefficients
! of the previous slice), and the bin counts of the coefficient solve that followed.
!
    integer, intent(in) :: unit, i_slice
    real(dp), intent(in) :: time, dp_tot, de_tot, dp_abs

    write(unit, '(i6, 6es16.7, 3i8)') i_slice, time, dp_tot, de_tot, dp_abs, rest_given(1), rest_given(2), &
                                      n_rest_exact, n_rest_fallback, n_rest_empty

end subroutine write_restoring_diagnostics

end module delta_f_restoring_mod
