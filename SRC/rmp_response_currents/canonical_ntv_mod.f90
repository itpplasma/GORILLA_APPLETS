!> Canonical delta-f and full-f toroidal-momentum (NTV) estimators.
!>
!> Background: the canonical Maxwellian f0(p_phi, H0) of the axisymmetric field
!> (ln_f0_and_H in utils_rmp_response_currents_mod), loaded through the birth
!> factor P = f0/g0 of the nonlinear loader. Perturbation: one harmonic
!>   dB1 = b_mn(s) cos(m theta - n phi),  dPhi1 = Phi_mn(s) cos(m theta - n phi),
!> with axisymmetric covariant B_k and A (no delta A_par). Per push of a marker
!> along the orbit of the field on the mesh:
!>   * linear delta-f weight  dw/dt = (d ln f0/dp_phi)|_H0 dH1/dphi,
!>     canonical dH1/dphi = (mu + m v_par^2/B) dB1/dphi + e dPhi1/dphi
!>     (generated kernels deltaf_weight_kernel and deltaf_harmonic_kernel);
!>   * delta-f torque  int w (-dH1/dphi) dt;
!>   * full-f  int (-dH1/dphi) dt and the p_phi change along the orbit from the
!>     mesh A_phi and h_phi;
!>   * collisional p_phi change m h_phi Delta v_par of each kick (tallied
!>     separately; the kicks leave w unchanged, so Lorentz scattering is a
!>     momentum sink of the test species without a restoring term);
!>   * birth and end (loss) p_phi, s and time.
!> The per-marker table (one row per marker, numeric columns, header line) is
!> written for the external NTV / delta-B-loss / axisymmetric-loss
!> decomposition of the paired +-sigma/0 runs.
module canonical_ntv_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use deltaf_harmonic_kernel, only: evaluate_deltaf_harmonic
    use deltaf_weight_kernel, only: evaluate_deltaf_weight

    implicit none

    private
    public :: boole_canonical_ntv, ntv_m_pol, ntv_n_tor, ntv_amplitude, ntv_bmn, &
        ntv_phimn, ntv_profile_file, ntv_theta_kind, ntv_output_file
    public :: canonical_ntv_nml_read, canonical_ntv_init, canonical_ntv_allocate, &
        canonical_ntv_birth, canonical_ntv_push, canonical_ntv_collision, &
        canonical_ntv_end, canonical_ntv_write, ntv_harmonic_dphi, &
        ntv_theta_of_x, species_temperature_eV, ntv_point_t, ntv_eval_point
    public :: n_ntv_cols, ntv_col_names

    logical  :: boole_canonical_ntv = .false.
    integer  :: ntv_m_pol = 0
    integer  :: ntv_n_tor = 0
    !> Amplitude sigma multiplying b_mn and Phi_mn (CRN pairs: +sigma, -sigma, 0).
    real(dp) :: ntv_amplitude = 0.0_dp
    !> Constant b_mn [G] and Phi_mn [statV] used when ntv_profile_file is empty.
    real(dp) :: ntv_bmn = 0.0_dp
    real(dp) :: ntv_phimn = 0.0_dp
    !> Optional table: s, b_mn [G], Phi_mn [statV] (s = psi_tor/psi_tor_edge).
    character(len=512) :: ntv_profile_file = ''
    !> Poloidal angle of the harmonic (grid_kind 5): 1 straight-field-line angle
    !> of B_phi = B0 R0/R (as GORILLA eval_helical_pert_circ), 0 geometric angle.
    integer  :: ntv_theta_kind = 1
    character(len=512) :: ntv_output_file = 'canonical_ntv_markers.dat'

    namelist /canonical_ntv_nml/ boole_canonical_ntv, ntv_m_pol, ntv_n_tor, &
        ntv_amplitude, ntv_bmn, ntv_phimn, ntv_profile_file, ntv_theta_kind, &
        ntv_output_file

    integer, parameter :: n_ntv_cols = 17
    character(len=*), parameter :: ntv_col_names = &
        'n lost P pphi_birth s_birth sstar_birth w_end tq_deltaf dh_int ' // &
        'dpphi_orbit dpphi_coll pphi_end s_end sstar_end t_end n_push ' // &
        'dlnf0_dpphi_birth'
    integer, parameter :: c_n = 1, c_lost = 2, c_p = 3, c_pphi0 = 4, c_s0 = 5, &
        c_sstar0 = 6, c_w = 7, c_tq = 8, c_dh = 9, c_dporb = 10, c_dpcoll = 11, &
        c_pphi1 = 12, c_s1 = 13, c_sstar1 = 14, c_t1 = 15, c_npush = 16, &
        c_dlnf0 = 17

    !> Per-marker records (n_ntv_cols, n_particles); one marker per thread.
    real(dp), allocatable :: rec(:,:)

    integer :: n_tab = 0
    real(dp), allocatable :: tab_s(:), tab_b(:), tab_b2(:), tab_p(:), tab_p2(:)

    !> Local state of a marker at one point of the mesh field.
    type :: ntv_point_t
        real(dp) :: bmod, h_phi, aphi, pphi, s_star, mu
        real(dp) :: dlnf0_dpphi, dh1_dphi
    end type ntv_point_t

contains

!> Reads canonical_ntv_nml from unit file fname; a missing group keeps defaults.
subroutine canonical_ntv_nml_read(fname)
    character(len=*), intent(in) :: fname
    integer :: u, ios

    open(newunit=u, file=fname, status='old', action='read', iostat=ios)
    if (ios /= 0) return
    read(u, nml=canonical_ntv_nml, iostat=ios)
    close(u)
    if (ios > 0) error stop 'canonical_ntv_nml: invalid namelist entry'
    if (.not. boole_canonical_ntv) return
    print '(a, 2i5, a, es12.4, a, 2es12.4)', ' Canonical NTV estimator: (m,n) =', &
        ntv_m_pol, ntv_n_tor, ', sigma =', ntv_amplitude, ', b_mn, Phi_mn =', &
        ntv_bmn, ntv_phimn
end subroutine canonical_ntv_nml_read

!> Checks the configuration and loads the optional amplitude table.
subroutine canonical_ntv_init()
    use tetra_grid_settings_mod, only: grid_kind
    integer :: u, ios, i
    real(dp) :: row(3)
    character(len=512) :: line

    if (grid_kind /= 5) error stop &
        'canonical NTV: only the analytic circular chart (grid_kind 5) is wired'
    if (ntv_theta_kind /= 0 .and. ntv_theta_kind /= 1) &
        error stop 'canonical NTV: ntv_theta_kind must be 0 or 1'
    n_tab = 0
    if (len_trim(ntv_profile_file) == 0) return
    open(newunit=u, file=trim(ntv_profile_file), status='old', action='read', &
         iostat=ios)
    if (ios /= 0) error stop 'canonical NTV: cannot open ntv_profile_file'
    do
        read(u, '(a)', iostat=ios) line
        if (ios /= 0) exit
        line = adjustl(line)
        if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
        n_tab = n_tab + 1
    end do
    if (n_tab < 2) error stop 'canonical NTV: ntv_profile_file needs >= 2 rows'
    if (allocated(tab_s)) deallocate(tab_s, tab_b, tab_b2, tab_p, tab_p2)
    allocate(tab_s(n_tab), tab_b(n_tab), tab_b2(n_tab), tab_p(n_tab), tab_p2(n_tab))
    rewind(u)
    i = 0
    do
        read(u, '(a)', iostat=ios) line
        if (ios /= 0) exit
        line = adjustl(line)
        if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
        i = i + 1
        read(line, *) row
        tab_s(i) = row(1)
        tab_b(i) = row(2)
        tab_p(i) = row(3)
    end do
    close(u)
    if (any(tab_s(2:) <= tab_s(:n_tab - 1))) &
        error stop 'canonical NTV: s in ntv_profile_file must increase'
    call spline_second_derivatives(tab_s, tab_b, tab_b2)
    call spline_second_derivatives(tab_s, tab_p, tab_p2)
end subroutine canonical_ntv_init

subroutine canonical_ntv_allocate(n_particles)
    integer, intent(in) :: n_particles
    integer :: n

    if (allocated(rec)) deallocate(rec)
    allocate(rec(n_ntv_cols, n_particles))
    rec = 0.0_dp
    do n = 1, n_particles
        rec(c_n, n) = real(n, dp)
    end do
end subroutine canonical_ntv_allocate

!> Temperature of the test species' background: T_i for ions, T_e otherwise.
pure function species_temperature_eV(Te, Ti, charge) result(T)
    real(dp), intent(in) :: Te, Ti, charge
    real(dp) :: T
    if (charge > 0.0_dp) then
        T = Ti
    else
        T = Te
    end if
end function species_temperature_eV

!> Poloidal angle of the harmonic at the cylindrical position x = (R, phi, Z).
pure function ntv_theta_of_x(x) result(theta)
    use tetra_grid_settings_mod, only: R0_analytic_circ
    real(dp), intent(in) :: x(3)
    real(dp) :: theta, rho, th, eps

    rho = sqrt((x(1) - R0_analytic_circ)**2 + x(3)**2)
    th = atan2(x(3), x(1) - R0_analytic_circ)
    if (ntv_theta_kind == 0) then
        theta = th
    else
        eps = rho / R0_analytic_circ
        theta = 2.0_dp * atan2(sqrt(1.0_dp - eps) * sin(0.5_dp * th), &
                               sqrt(1.0_dp + eps) * cos(0.5_dp * th))
    end if
end function ntv_theta_of_x

!> dB1/dphi and dPhi1/dphi of the scaled harmonic at (s, theta, phi).
subroutine ntv_harmonic_dphi(s, theta, phi, db1_dphi, dphi1_dphi)
    real(dp), intent(in)  :: s, theta, phi
    real(dp), intent(out) :: db1_dphi, dphi1_dphi
    real(dp) :: bmn, phimn, db1, db1_ds, db1_dth, dp1, dp1_ds, dp1_dth

    if (n_tab > 0) then
        bmn = spline_value(tab_s, tab_b, tab_b2, s)
        phimn = spline_value(tab_s, tab_p, tab_p2, s)
    else
        bmn = ntv_bmn
        phimn = ntv_phimn
    end if
    bmn = ntv_amplitude * bmn
    phimn = ntv_amplitude * phimn
    ! Radial derivatives are not needed for dH1/dphi.
    call evaluate_deltaf_harmonic(real(ntv_m_pol, dp), real(ntv_n_tor, dp), bmn, &
        0.0_dp, phimn, 0.0_dp, theta, phi, db1, db1_ds, db1_dth, db1_dphi, dp1, &
        dp1_ds, dp1_dth, dphi1_dphi)
end subroutine ntv_harmonic_dphi

!> Local canonical quantities of a marker (vpar, perpinv) at x in cell ind_tetr.
!> psi0: unperturbed poloidal flux (the A_phi of f0), s0 = s(psi0), Phi0 the
!> profile potential at s0; s_mesh: flux label for the amplitude profile.
subroutine ntv_eval_point(ind_tetr, x, vpar, perpinv, species, psi0, s_mesh, pt)
    use tetra_physics_mod, only: tetra_physics
    use gorilla_applets_types_mod, only: start
    use profile_data_mod, only: eval_profiles, profile_values_t, &
        eval_s_and_ds_dpsi_from_psi_pol
    use constants, only: ev2erg, clight

    integer,  intent(in)  :: ind_tetr, species
    real(dp), intent(in)  :: x(3), vpar, perpinv, psi0, s_mesh
    type(ntv_point_t), intent(out) :: pt

    real(dp) :: z(3), mass, charge, K, s0, ds0, dsdpsi, psi_star, T_eV, dT_ds
    real(dp) :: db1_dphi, dphi1_dphi, ln_f0, dw_dt, torque, h0, dens, ddens, temp
    real(dp) :: dtemp, phi0, dphi0
    type(profile_values_t) :: pv0, pvs

    mass = start%particle_mass(species)
    charge = start%particle_charge(species)
    z = x - tetra_physics(ind_tetr)%x1
    pt%bmod = tetra_physics(ind_tetr)%bmod1 + sum(tetra_physics(ind_tetr)%gB * z)
    pt%h_phi = tetra_physics(ind_tetr)%h2_1 + sum(tetra_physics(ind_tetr)%gh2 * z)
    pt%aphi = tetra_physics(ind_tetr)%Aphi1 + sum(tetra_physics(ind_tetr)%gAphi * z)
    pt%pphi = charge * pt%aphi / clight + mass * vpar * pt%h_phi
    pt%mu = -mass * perpinv
    K = mass * (0.5_dp * vpar**2 - perpinv * pt%bmod)

    ! f0 on the canonical label of the unperturbed field (as ln_f0_and_H):
    ! psi* = psi0 + (c m/e) vpar h_phi, H0 = K + e Phi0(s0).
    psi_star = psi0 + start%cm_over_e(species) * vpar * pt%h_phi
    call eval_s_and_ds_dpsi_from_psi_pol(psi0, s0, ds0)
    call eval_s_and_ds_dpsi_from_psi_pol(psi_star, pt%s_star, dsdpsi)
    call eval_profiles(s0, pv0)
    call eval_profiles(pt%s_star, pvs)
    h0 = K + charge * pv0%Phi0
    T_eV = species_temperature_eV(pvs%Te, pvs%Ti, charge)
    dT_ds = species_temperature_eV(pvs%dTe_ds, pvs%dTi_ds, charge)
    dens = pvs%n_e
    ddens = pvs%dn_ds * dsdpsi
    temp = T_eV * ev2erg
    dtemp = dT_ds * ev2erg * dsdpsi
    phi0 = pvs%Phi0
    dphi0 = pvs%dPhi0_ds * dsdpsi

    call ntv_harmonic_dphi(s_mesh, ntv_theta_of_x(x), x(2), db1_dphi, dphi1_dphi)
    call evaluate_deltaf_weight(mass, charge, clight, pt%mu, vpar, pt%bmod, h0, &
        dens, ddens, temp, dtemp, phi0, dphi0, db1_dphi, dphi1_dphi, ln_f0, &
        pt%dlnf0_dpphi, pt%dh1_dphi, dw_dt, torque)
end subroutine ntv_eval_point

!> Records the birth state of marker n (P: birth factor f0/g0).
subroutine canonical_ntv_birth(n, P, pt, s0)
    integer,  intent(in) :: n
    real(dp), intent(in) :: P, s0
    type(ntv_point_t), intent(in) :: pt

    rec(c_p, n) = P
    rec(c_pphi0, n) = pt%pphi
    rec(c_s0, n) = s0
    rec(c_sstar0, n) = pt%s_star
    rec(c_dlnf0, n) = pt%dlnf0_dpphi
    rec(c_pphi1, n) = pt%pphi
    rec(c_s1, n) = s0
    rec(c_sstar1, n) = pt%s_star
end subroutine canonical_ntv_birth

!> One push: pa, pb at the start and end point (same cell), pm at the
!> midpoint; dt is the physical time of the push. Advances w (in/out).
subroutine canonical_ntv_push(n, pa, pm, pb, s_end, dt, w)
    integer,  intent(in) :: n
    type(ntv_point_t), intent(in) :: pa, pm, pb
    real(dp), intent(in) :: s_end, dt
    real(dp), intent(inout) :: w
    real(dp) :: w_old

    w_old = w
    w = w + pm%dlnf0_dpphi * pm%dh1_dphi * dt
    rec(c_tq, n) = rec(c_tq, n) - 0.5_dp * (w_old + w) * pm%dh1_dphi * dt
    rec(c_dh, n) = rec(c_dh, n) - pm%dh1_dphi * dt
    rec(c_dporb, n) = rec(c_dporb, n) + (pb%pphi - pa%pphi)
    rec(c_w, n) = w
    rec(c_pphi1, n) = pb%pphi
    rec(c_s1, n) = s_end
    rec(c_sstar1, n) = pb%s_star
    rec(c_npush, n) = rec(c_npush, n) + 1.0_dp
end subroutine canonical_ntv_push

!> Collisional p_phi change m h_phi (vpar_new - vpar_old) at fixed position.
subroutine canonical_ntv_collision(n, species, h_phi, vpar_old, vpar_new)
    use gorilla_applets_types_mod, only: start
    integer,  intent(in) :: n, species
    real(dp), intent(in) :: h_phi, vpar_old, vpar_new

    rec(c_dpcoll, n) = rec(c_dpcoll, n) &
        + start%particle_mass(species) * h_phi * (vpar_new - vpar_old)
    rec(c_pphi1, n) = rec(c_pphi1, n) &
        + start%particle_mass(species) * h_phi * (vpar_new - vpar_old)
end subroutine canonical_ntv_collision

subroutine canonical_ntv_end(n, lost, t_end)
    integer,  intent(in) :: n
    logical,  intent(in) :: lost
    real(dp), intent(in) :: t_end

    rec(c_lost, n) = merge(1.0_dp, 0.0_dp, lost)
    rec(c_t1, n) = t_end
end subroutine canonical_ntv_end

!> Writes the per-marker table and prints population sums (P-weighted).
subroutine canonical_ntv_write(volume, trace_time)
    real(dp), intent(in) :: volume, trace_time
    integer :: u, n, np
    real(dp) :: sum_tq, sum_dh, sum_dp, sum_dc, max_w, max_dp

    np = size(rec, 2)
    open(newunit=u, file=trim(ntv_output_file), status='replace', action='write')
    write(u, '(a, i0, a, i0, a, es24.16)') '# canonical NTV markers: m = ', &
        ntv_m_pol, ', n = ', ntv_n_tor, ', sigma = ', ntv_amplitude
    write(u, '(a, es24.16, a, es24.16)') '# spawn_volume_cm3 ', volume, &
        ' trace_time_s ', trace_time
    write(u, '(a)') '# torque density estimate = volume/N sum P tq_deltaf / ' // &
        'trace_time (dyn cm / cm^3 per marker population)'
    write(u, '(2a)') '# ', ntv_col_names
    do n = 1, np
        write(u, '(i10, 1x, i2, 15(1x, es24.16))') nint(rec(c_n, n)), &
            nint(rec(c_lost, n)), rec(c_p:c_t1, n), rec(c_npush:c_dlnf0, n)
    end do
    close(u)
    sum_tq = sum(rec(c_p, :) * rec(c_tq, :))
    sum_dh = sum(rec(c_p, :) * rec(c_dh, :))
    sum_dp = sum(rec(c_p, :) * rec(c_dporb, :))
    sum_dc = sum(rec(c_p, :) * rec(c_dpcoll, :))
    max_w = maxval(abs(rec(c_w, :)))
    max_dp = maxval(abs(rec(c_dporb, :)))
    print '(a, es24.16)', ' Canonical NTV: sum P tq_deltaf = ', sum_tq
    print '(a, es24.16)', ' Canonical NTV: sum P dh_int = ', sum_dh
    print '(a, es24.16)', ' Canonical NTV: sum P dpphi_orbit = ', sum_dp
    print '(a, es24.16)', ' Canonical NTV: sum P dpphi_coll = ', sum_dc
    print '(a, es24.16)', ' Canonical NTV: max |w| = ', max_w
    print '(a, es24.16)', ' Canonical NTV: max |dpphi_orbit| = ', max_dp
    print '(a, i0)', ' Canonical NTV: lost markers = ', nint(sum(rec(c_lost, :)))
end subroutine canonical_ntv_write

subroutine spline_second_derivatives(x, y, y2)
    real(dp), intent(in)  :: x(:), y(:)
    real(dp), intent(out) :: y2(:)
    real(dp) :: u(size(x)), sig, p
    integer :: i, n

    n = size(x)
    y2(1) = 0.0_dp
    u(1) = 0.0_dp
    do i = 2, n - 1
        sig = (x(i) - x(i - 1)) / (x(i + 1) - x(i - 1))
        p = sig * y2(i - 1) + 2.0_dp
        y2(i) = (sig - 1.0_dp) / p
        u(i) = (6.0_dp * ((y(i + 1) - y(i)) / (x(i + 1) - x(i)) - (y(i) - y(i - 1)) &
               / (x(i) - x(i - 1))) / (x(i + 1) - x(i - 1)) - sig * u(i - 1)) / p
    end do
    y2(n) = 0.0_dp
    do i = n - 1, 1, -1
        y2(i) = y2(i) * y2(i + 1) + u(i)
    end do
end subroutine spline_second_derivatives

!> Natural cubic spline value; zero outside the tabulated range.
pure function spline_value(x, y, y2, xv) result(val)
    real(dp), intent(in) :: x(:), y(:), y2(:), xv
    real(dp) :: val, h, a, b
    integer :: lo, hi, k

    val = 0.0_dp
    if (xv < x(1) .or. xv > x(size(x))) return
    lo = 1
    hi = size(x)
    do while (hi - lo > 1)
        k = (hi + lo) / 2
        if (x(k) > xv) then
            hi = k
        else
            lo = k
        end if
    end do
    h = x(hi) - x(lo)
    a = (x(hi) - xv) / h
    b = (xv - x(lo)) / h
    val = a * y(lo) + b * y(hi) + ((a**3 - a) * y2(lo) + (b**3 - b) * y2(hi)) &
          * h**2 / 6.0_dp
end function spline_value

end module canonical_ntv_mod
