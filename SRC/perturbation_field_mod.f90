!
! perturbation_field_mod.f90
!
! Loads a prescribed magnetic perturbation field delta_B_s for a single
! (m,n) Fourier harmonic from a text file and evaluates it along particle
! orbits.
!
! The perturbation is stored as a 1D radial profile (complex amplitude vs.
! the GORILLA flux coordinate s).  Along an orbit at position (s, theta, phi),
! the real perturbation is reconstructed as:
!
!   delta_B_s(s, theta, phi) = Re[ delta_B_hat(s) * exp(i*(m*theta + n*phi)) ]
!
! Input file format (three columns):
!   r_eff [cm]    Re(delta_B_s) [Gauss]    Im(delta_B_s) [Gauss]
!
! The r_eff -> s mapping uses the same equil_r_q_psi.dat as profile_data_mod.
!
!
! Coordinate note: the r_eff coordinate used throughout this module is the
! toroidal-flux-based effective minor radius, r_eff = sqrt(2*psi_tor/B_ref).
! It is NOT the geometric minor radius r_geom (arithmetic mean of flux-surface
! half-widths in the horizontal plane, ~50 cm at the AUG edge vs ~62 cm for
! r_eff).  All input files and namelist parameters that accept r in [cm] use
! r_eff unless explicitly stated otherwise.  The confusion between these two
! led to a factor-of-2 amplitude bug (fixed 2026-05-19) -- do not swap them.
!
module perturbation_field_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

    private
    public :: load_perturbation_field, init_constant_perturbation, &
              init_step_perturbation, load_eperp_field, &
              eval_delta_B_s, eval_delta_E_s, cleanup_perturbation_field
    public :: pert_m_mode, pert_n_mode

    ! Mode numbers
    integer :: pert_m_mode = 0
    integer :: pert_n_mode = 0

    ! Number of grid points
    integer :: ns_pert

    ! Spline data on s grid
    real(dp), allocatable :: s_pert_grid(:)
    real(dp), allocatable :: dB_re_spl(:), dB_re_dd(:)
    real(dp), allocatable :: dB_im_spl(:), dB_im_dd(:)

    logical :: perturbation_loaded = .false.

    ! Constant-amplitude mode: dB_s = dB_const * cos(m theta + n phi).
    ! Skips the radial spline and just returns dB_const at every s.
    logical :: use_constant_amplitude = .false.
    real(dp) :: dB_constant = 0.0_dp

    ! Step-function mode: Br = dB_constant inside [s_step_min, s_step_max], 0 outside.
    ! Activated by init_step_perturbation; s bounds are computed from r_eff inputs.
    logical  :: use_step_function = .false.
    real(dp) :: s_step_min = 0.0_dp
    real(dp) :: s_step_max = 1.0_dp

    ! E_perp perturbation (KIM electrostatic field, statV/cm).
    ! Loaded by load_eperp_field; evaluated by eval_delta_E_s.
    ! Mode numbers are shared with the B-field perturbation (pert_m_mode/pert_n_mode).
    logical  :: boole_eperp_loaded = .false.
    integer  :: ns_eperp = 0
    real(dp), allocatable :: s_eperp_grid(:)
    real(dp), allocatable :: eperp_re_spl(:), eperp_re_dd(:)
    real(dp), allocatable :: eperp_im_spl(:), eperp_im_dd(:)

    ! Diagnostic switch: if .true., skip the cos(m theta + n phi) factor
    ! and return dB_constant as a truly uniform (axisymmetric) value. Used
    ! to verify that the resonant structure in n=2 originates from the
    ! helical phase, not from any other coupling.
    logical, public :: boole_skip_phase = .false.

contains

! ============================================================
! Load perturbation field and build splines in s
! ============================================================
subroutine load_perturbation_field(delta_B_file, equil_mapping_file, m_mode, n_mode)

    use profile_data_mod, only: kim_b_ref, boole_kim_reff_coords

    character(len=*), intent(in) :: delta_B_file
    character(len=*), intent(in) :: equil_mapping_file
    integer, intent(in) :: m_mode, n_mode

    integer :: n_equil, n_raw, iunit, ios, i, j
    real(dp), allocatable :: r_equil(:), psi_tor_equil(:), r_map(:)
    real(dp), allocatable :: r_raw(:), dB_re_raw(:), dB_im_raw(:)
    character(len=512) :: line
    real(dp) :: psi_tor_edge, frac, r_val, re_val, im_val

    pert_m_mode = m_mode
    pert_n_mode = n_mode

    print *, ''
    print *, 'Loading perturbation field...'
    print *, '  File: ', trim(delta_B_file)
    print *, '  Mode: m = ', m_mode, ', n = ', n_mode

    ! --- Load equil mapping for r_eff -> s conversion ---
    call read_equil_mapping_pert(equil_mapping_file, n_equil, r_equil, psi_tor_equil)

    psi_tor_edge = psi_tor_equil(n_equil)
    ns_pert = n_equil
    allocate(s_pert_grid(ns_pert))
    if (abs(psi_tor_edge) > 0.0_dp) then
        s_pert_grid = psi_tor_equil / psi_tor_edge
    else
        do i = 1, ns_pert
            s_pert_grid(i) = real(i - 1, dp) / real(ns_pert - 1, dp)
        end do
    end if
    s_pert_grid(1) = 0.0_dp
    s_pert_grid(ns_pert) = 1.0_dp

    ! --- Read perturbation file ---
    iunit = 44
    open(unit=iunit, file=delta_B_file, status='old', action='read', iostat=ios)
    if (ios /= 0) then
        print *, 'ERROR: Cannot open ', trim(delta_B_file)
        stop
    end if

    ! Count data lines
    n_raw = 0
    do
        read(iunit, '(A)', iostat=ios) line
        if (ios /= 0) exit
        if (line(1:1) == '#' .or. len_trim(line) == 0) cycle
        n_raw = n_raw + 1
    end do
    close(iunit)

    allocate(r_raw(n_raw), dB_re_raw(n_raw), dB_im_raw(n_raw))

    open(unit=iunit, file=delta_B_file, status='old', action='read')
    i = 0
    do
        read(iunit, '(A)', iostat=ios) line
        if (ios /= 0) exit
        if (line(1:1) == '#' .or. len_trim(line) == 0) cycle
        i = i + 1
        read(line, *) r_raw(i), dB_re_raw(i), dB_im_raw(i)
    end do
    close(iunit)

    print *, '  Loaded ', n_raw, ' data points'

    ! --- Interpolate onto the GORILLA grid ---
    ! r_map(i) is the radius at grid point i that the FILE's r-column is matched
    ! against.  Legacy: r_geom (flux_functions col1).  KIM r_eff alignment:
    ! r_eff = sqrt(2*psi_tor/B_ref), so a KIM-coordinate file lands on the
    ! matching flux surface (see profile_data_mod note).
    allocate(r_map(n_equil))
    if (boole_kim_reff_coords .and. kim_b_ref > 0.0_dp) then
        r_map = sqrt(2.0_dp * abs(psi_tor_equil) / kim_b_ref)
    else
        r_map = r_equil
    end if

    allocate(dB_re_spl(ns_pert), dB_re_dd(ns_pert))
    allocate(dB_im_spl(ns_pert), dB_im_dd(ns_pert))

    do i = 1, ns_pert
        if (r_map(i) <= r_raw(1)) then
            dB_re_spl(i) = dB_re_raw(1)
            dB_im_spl(i) = dB_im_raw(1)
        else if (r_map(i) >= r_raw(n_raw)) then
            dB_re_spl(i) = dB_re_raw(n_raw)
            dB_im_spl(i) = dB_im_raw(n_raw)
        else
            do j = 1, n_raw - 1
                if (r_raw(j+1) >= r_map(i)) exit
            end do
            frac = (r_map(i) - r_raw(j)) / (r_raw(j+1) - r_raw(j))
            dB_re_spl(i) = dB_re_raw(j) + frac * (dB_re_raw(j+1) - dB_re_raw(j))
            dB_im_spl(i) = dB_im_raw(j) + frac * (dB_im_raw(j+1) - dB_im_raw(j))
        end if
    end do
    deallocate(r_map)

    ! --- Build cubic splines in s ---
    call spline_natural_pert(ns_pert, s_pert_grid, dB_re_spl, dB_re_dd)
    call spline_natural_pert(ns_pert, s_pert_grid, dB_im_spl, dB_im_dd)

    perturbation_loaded = .true.

    print *, '  Perturbation field splined on ', ns_pert, ' s-grid points'
    print *, '  |delta_B|(s=0) = ', sqrt(dB_re_spl(1)**2 + dB_im_spl(1)**2), ' Gauss'
    print *, ''

    deallocate(r_equil, psi_tor_equil, r_raw, dB_re_raw, dB_im_raw)

end subroutine load_perturbation_field

! ============================================================
! Constant-amplitude initialisation
!
! Sets up a radially constant perturbation envelope:
!   delta_B_s = dB_const * cos(m*theta + n*phi)
! ============================================================
subroutine init_constant_perturbation(dB_const, m_mode, n_mode)

    real(dp), intent(in) :: dB_const
    integer, intent(in) :: m_mode, n_mode

    use_constant_amplitude = .true.
    dB_constant = dB_const
    pert_m_mode = m_mode
    pert_n_mode = n_mode
    perturbation_loaded = .true.

    print *, ''
    print *, 'Perturbation field: constant amplitude mode'
    print *, '  dB_s_const = ', dB_const, ' Gauss'
    print *, '  Mode: m = ', m_mode, ', n = ', n_mode
    print *, ''

end subroutine init_constant_perturbation

! ============================================================
! Step-function-amplitude initialisation
!
! Sets up a step-function radial envelope:
!   delta_B_r(s) = dB_const  if s_step_min <= s <= s_step_max
!                = 0          otherwise
!
! center_reff and halfwidth_reff are in r_eff [cm].
! r_eff is the TOROIDAL-FLUX-BASED effective radius r_eff = sqrt(2*psi_tor/B_ref).
! It is NOT the geometric minor radius r_geom (~50 cm at AUG edge vs ~62 cm
! for r_eff).  Use the value from flux_functions.dat col1 (0-indexed).
!
! The r_eff -> s conversion reads equil_mapping_file (flux_functions.dat).
! ============================================================
subroutine init_step_perturbation(center_reff, halfwidth_reff, dB_const, &
                                  m_mode, n_mode, equil_mapping_file)

    real(dp),         intent(in) :: center_reff      ! r_eff [cm] — NOT r_geom
    real(dp),         intent(in) :: halfwidth_reff   ! r_eff [cm] — NOT r_geom
    real(dp),         intent(in) :: dB_const
    integer,          intent(in) :: m_mode, n_mode
    character(len=*), intent(in) :: equil_mapping_file

    integer  :: n_equil, i
    real(dp), allocatable :: r_equil(:), psi_tor_equil(:), s_equil(:)
    real(dp) :: psi_tor_edge, frac, r_lo, r_hi

    if (halfwidth_reff <= 0.0_dp) then
        print *, 'ERROR: step_halfwidth_reff must be > 0, got ', halfwidth_reff
        stop
    end if

    call read_equil_mapping_pert(equil_mapping_file, n_equil, r_equil, psi_tor_equil)

    psi_tor_edge = psi_tor_equil(n_equil)
    allocate(s_equil(n_equil))
    if (abs(psi_tor_edge) > 0.0_dp) then
        s_equil = psi_tor_equil / psi_tor_edge
    else
        do i = 1, n_equil
            s_equil(i) = real(i-1, dp) / real(n_equil-1, dp)
        end do
    end if
    s_equil(1) = 0.0_dp
    s_equil(n_equil) = 1.0_dp

    ! Convert center_reff - halfwidth_reff -> s_step_min
    r_lo = center_reff - halfwidth_reff
    if (r_lo <= r_equil(1)) then
        s_step_min = 0.0_dp
    else if (r_lo >= r_equil(n_equil)) then
        s_step_min = 1.0_dp
    else
        do i = 1, n_equil - 1
            if (r_equil(i+1) >= r_lo) exit
        end do
        frac = (r_lo - r_equil(i)) / (r_equil(i+1) - r_equil(i))
        s_step_min = s_equil(i) + frac * (s_equil(i+1) - s_equil(i))
    end if

    ! Convert center_reff + halfwidth_reff -> s_step_max
    r_hi = center_reff + halfwidth_reff
    if (r_hi <= r_equil(1)) then
        s_step_max = 0.0_dp
    else if (r_hi >= r_equil(n_equil)) then
        s_step_max = 1.0_dp
    else
        do i = 1, n_equil - 1
            if (r_equil(i+1) >= r_hi) exit
        end do
        frac = (r_hi - r_equil(i)) / (r_equil(i+1) - r_equil(i))
        s_step_max = s_equil(i) + frac * (s_equil(i+1) - s_equil(i))
    end if

    use_constant_amplitude = .true.
    use_step_function      = .true.
    dB_constant   = dB_const
    pert_m_mode   = m_mode
    pert_n_mode   = n_mode
    perturbation_loaded = .true.

    print *, ''
    print *, 'Perturbation field: step-function mode'
    print *, '  dB_const        = ', dB_const, ' Gauss'
    print *, '  centre  r_eff   = ', center_reff,   ' cm  (NOT r_geom)'
    print *, '  halfwidth r_eff = ', halfwidth_reff, ' cm  (NOT r_geom)'
    print *, '  s_step_min      = ', s_step_min
    print *, '  s_step_max      = ', s_step_max
    print *, '  Mode: m = ', m_mode, ', n = ', n_mode
    print *, ''

    deallocate(r_equil, psi_tor_equil, s_equil)

end subroutine init_step_perturbation

! ============================================================
! Load E_perp perturbation field from KIM output file and
! build cubic splines over s.
!
! File format (4 columns, KIM E_perp.dat convention):
!   r_eff [cm]   Re(E_perp) [statV/cm]   Im(E_perp) [statV/cm]   |E_perp| (ignored)
!
! r_eff is the toroidal-flux-based effective radius (NOT r_geom).
! Mode numbers (m, n) are reused from module state set by a prior
! B-field init call (init_constant_perturbation, init_step_perturbation,
! or load_perturbation_field).  Call this AFTER one of those.
! ============================================================
subroutine load_eperp_field(eperp_file, equil_mapping_file)

    use profile_data_mod, only: kim_b_ref, boole_kim_reff_coords

    character(len=*), intent(in) :: eperp_file
    character(len=*), intent(in) :: equil_mapping_file

    integer :: n_equil, n_raw, iunit, ios, i, j
    real(dp), allocatable :: r_equil(:), psi_tor_equil(:), r_map(:)
    real(dp), allocatable :: r_raw(:), eperp_re_raw(:), eperp_im_raw(:)
    character(len=512) :: line
    real(dp) :: psi_tor_edge, frac

    print *, ''
    print *, 'Loading E_perp perturbation field...'
    print *, '  File: ', trim(eperp_file)
    print *, '  Mode (shared with B-field): m = ', pert_m_mode, ', n = ', pert_n_mode

    ! --- Load equil mapping for r_eff -> s conversion ---
    call read_equil_mapping_pert(equil_mapping_file, n_equil, r_equil, psi_tor_equil)

    psi_tor_edge = psi_tor_equil(n_equil)
    ns_eperp = n_equil
    allocate(s_eperp_grid(ns_eperp))
    if (abs(psi_tor_edge) > 0.0_dp) then
        s_eperp_grid = psi_tor_equil / psi_tor_edge
    else
        do i = 1, ns_eperp
            s_eperp_grid(i) = real(i - 1, dp) / real(ns_eperp - 1, dp)
        end do
    end if
    s_eperp_grid(1) = 0.0_dp
    s_eperp_grid(ns_eperp) = 1.0_dp

    ! --- Read E_perp file (4-column KIM format) ---
    iunit = 46
    open(unit=iunit, file=eperp_file, status='old', action='read', iostat=ios)
    if (ios /= 0) then
        print *, 'ERROR: Cannot open ', trim(eperp_file)
        stop
    end if

    n_raw = 0
    do
        read(iunit, '(A)', iostat=ios) line
        if (ios /= 0) exit
        if (line(1:1) == '#' .or. len_trim(line) == 0) cycle
        n_raw = n_raw + 1
    end do
    close(iunit)

    allocate(r_raw(n_raw), eperp_re_raw(n_raw), eperp_im_raw(n_raw))

    open(unit=iunit, file=eperp_file, status='old', action='read')
    i = 0
    do
        read(iunit, '(A)', iostat=ios) line
        if (ios /= 0) exit
        if (line(1:1) == '#' .or. len_trim(line) == 0) cycle
        i = i + 1
        read(line, *) r_raw(i), eperp_re_raw(i), eperp_im_raw(i)
    end do
    close(iunit)

    print *, '  Loaded ', n_raw, ' data points'

    ! --- Interpolate onto the GORILLA grid ---
    ! KIM's E_perp.dat r-column is the toroidal-flux radius r_eff; match it
    ! against r_eff(s)=sqrt(2*psi_tor/B_ref) when alignment is on, else r_geom.
    allocate(r_map(n_equil))
    if (boole_kim_reff_coords .and. kim_b_ref > 0.0_dp) then
        r_map = sqrt(2.0_dp * abs(psi_tor_equil) / kim_b_ref)
    else
        r_map = r_equil
    end if

    allocate(eperp_re_spl(ns_eperp), eperp_re_dd(ns_eperp))
    allocate(eperp_im_spl(ns_eperp), eperp_im_dd(ns_eperp))

    do i = 1, ns_eperp
        if (r_map(i) <= r_raw(1)) then
            eperp_re_spl(i) = eperp_re_raw(1)
            eperp_im_spl(i) = eperp_im_raw(1)
        else if (r_map(i) >= r_raw(n_raw)) then
            eperp_re_spl(i) = eperp_re_raw(n_raw)
            eperp_im_spl(i) = eperp_im_raw(n_raw)
        else
            do j = 1, n_raw - 1
                if (r_raw(j+1) >= r_map(i)) exit
            end do
            frac = (r_map(i) - r_raw(j)) / (r_raw(j+1) - r_raw(j))
            eperp_re_spl(i) = eperp_re_raw(j) + frac * (eperp_re_raw(j+1) - eperp_re_raw(j))
            eperp_im_spl(i) = eperp_im_raw(j) + frac * (eperp_im_raw(j+1) - eperp_im_raw(j))
        end if
    end do
    deallocate(r_map)

    ! --- Build cubic splines in s ---
    call spline_natural_pert(ns_eperp, s_eperp_grid, eperp_re_spl, eperp_re_dd)
    call spline_natural_pert(ns_eperp, s_eperp_grid, eperp_im_spl, eperp_im_dd)

    boole_eperp_loaded = .true.

    print *, '  E_perp splined on ', ns_eperp, ' s-grid points'
    print *, '  |E_perp|(s=0) = ', &
             sqrt(eperp_re_spl(1)**2 + eperp_im_spl(1)**2), ' statV/cm'
    print *, ''

    deallocate(r_equil, psi_tor_equil, r_raw, eperp_re_raw, eperp_im_raw)

end subroutine load_eperp_field

! ============================================================
! Evaluate delta_E_perp at a particle position (s, theta, phi).
!
! Returns the complex E_perp amplitude including the helical phase:
!   dE_s = (ds/dr_eff)(s) * E_perp_hat(s) * exp(i*(m*theta + n*phi))
!
! Returns zero if load_eperp_field has not been called.
! ============================================================
subroutine eval_delta_E_s(s_val, theta, phi, dE_s)

    use profile_data_mod, only: eval_ds_dreff

    real(dp),    intent(in)  :: s_val, theta, phi
    complex(dp), intent(out) :: dE_s

    real(dp) :: s_clamped, eperp_re, eperp_im, phase, dummy, ds_dreff
    complex(dp), parameter :: i_imag = (0.0_dp, 1.0_dp)

    if (.not. boole_eperp_loaded) then
        dE_s = cmplx(0.0_dp, 0.0_dp, kind=dp)
        return
    end if

    s_clamped = max(0.0_dp, min(1.0_dp, s_val))
    phase = real(pert_m_mode, dp) * theta + real(pert_n_mode, dp) * phi
    ds_dreff = eval_ds_dreff(s_val)

    call splint_pert(ns_eperp, s_eperp_grid, eperp_re_spl, eperp_re_dd, s_clamped, eperp_re, dummy)
    call splint_pert(ns_eperp, s_eperp_grid, eperp_im_spl, eperp_im_dd, s_clamped, eperp_im, dummy)

    dE_s = ds_dreff * cmplx(eperp_re, eperp_im, kind=dp) * exp(i_imag * phase)

end subroutine eval_delta_E_s

! ============================================================
! Evaluate delta_B^s at a particle position (s, theta, phi).
!
! The user-supplied input (dB_constant or the file-loaded amplitude)
! is the radial perturbation delta_B_r in Gauss.  We convert to the
! contravariant component delta_B^s = (ds/dr_eff) * delta_B_r [Gauss/cm]
! using eval_ds_dreff(s) from profile_data_mod.
!
! Constant mode: dB_s = (ds/dr_eff)(s) * dB_r_const * exp(i(m theta + n phi))
! File mode:     dB_s = (ds/dr_eff)(s) * Re[ dB_hat_r(s) * exp(i(...)) ]
! ============================================================
subroutine eval_delta_B_s(s_val, theta, phi, dB_s)

    use profile_data_mod, only: eval_ds_dreff

    real(dp),    intent(in)  :: s_val, theta, phi
    complex(dp), intent(out) :: dB_s

    real(dp) :: s_clamped, dB_re, dB_im, phase, dummy, ds_dreff
    complex(dp), parameter :: i_imag = (0.0_dp, 1.0_dp)

    phase = real(pert_m_mode, dp) * theta + real(pert_n_mode, dp) * phi
    ds_dreff = eval_ds_dreff(s_val)

    if (use_constant_amplitude) then
        if (use_step_function) then
            if (s_val < s_step_min .or. s_val > s_step_max) then
                dB_s = cmplx(0.0_dp, 0.0_dp, kind=dp)
                return
            end if
        end if
        if (boole_skip_phase) then
            dB_s = cmplx(ds_dreff * dB_constant, 0.0_dp, kind=dp)
        else
            dB_s = ds_dreff * dB_constant * exp(i_imag * phase)
        end if
        return
    end if

    s_clamped = max(0.0_dp, min(1.0_dp, s_val))

    call splint_pert(ns_pert, s_pert_grid, dB_re_spl, dB_re_dd, s_clamped, dB_re, dummy)
    call splint_pert(ns_pert, s_pert_grid, dB_im_spl, dB_im_dd, s_clamped, dB_im, dummy)

    dB_s = ds_dreff * cmplx(dB_re, dB_im, kind=dp) * exp(i_imag * phase)

end subroutine eval_delta_B_s

! ============================================================
! Cleanup
! ============================================================
subroutine cleanup_perturbation_field()

    if (allocated(s_pert_grid)) deallocate(s_pert_grid)
    if (allocated(dB_re_spl))  deallocate(dB_re_spl)
    if (allocated(dB_re_dd))   deallocate(dB_re_dd)
    if (allocated(dB_im_spl))  deallocate(dB_im_spl)
    if (allocated(dB_im_dd))   deallocate(dB_im_dd)
    perturbation_loaded = .false.
    use_constant_amplitude = .false.
    dB_constant = 0.0_dp
    use_step_function = .false.
    s_step_min = 0.0_dp
    s_step_max = 1.0_dp
    if (allocated(s_eperp_grid))  deallocate(s_eperp_grid)
    if (allocated(eperp_re_spl))  deallocate(eperp_re_spl)
    if (allocated(eperp_re_dd))   deallocate(eperp_re_dd)
    if (allocated(eperp_im_spl))  deallocate(eperp_im_spl)
    if (allocated(eperp_im_dd))   deallocate(eperp_im_dd)
    boole_eperp_loaded = .false.
    ns_eperp = 0

end subroutine cleanup_perturbation_field

! ============================================================
! PRIVATE: Read equil mapping (simplified, only need r and psi_tor)
! ============================================================
subroutine read_equil_mapping_pert(filename, n, r_eff, psi_tor)

    character(len=*), intent(in) :: filename
    integer, intent(out) :: n
    real(dp), allocatable, intent(out) :: r_eff(:), psi_tor(:)

    ! Reads flux_functions.dat (5 columns, 0-indexed):
    !   col0: R_beg [cm]  (skipped)
    !   col1: r_eff [cm]  -- toroidal-flux-based effective radius, NOT r_geom
    !   col2: q_FL        (skipped)
    !   col3: psi_pol     (skipped)
    !   col4: psi_tor     [G cm^2]
    integer :: iunit, ios, i
    character(len=512) :: line
    real(dp) :: c0, c1, c2, c3, c4

    iunit = 45
    open(unit=iunit, file=filename, status='old', action='read', iostat=ios)
    if (ios /= 0) then
        print *, 'ERROR: Cannot open ', trim(filename)
        stop
    end if

    n = 0
    do
        read(iunit, '(A)', iostat=ios) line
        if (ios /= 0) exit
        if (line(1:1) == '#' .or. line(1:2) == ' #' .or. len_trim(line) == 0) cycle
        n = n + 1
    end do
    close(iunit)

    if (n == 0) then
        print *, 'ERROR: no data rows found in equil_mapping_file: ', trim(filename)
        stop
    end if

    allocate(r_eff(n), psi_tor(n))

    open(unit=iunit, file=filename, status='old', action='read')
    i = 0
    do
        read(iunit, '(A)', iostat=ios) line
        if (ios /= 0) exit
        if (line(1:1) == '#' .or. line(1:2) == ' #' .or. len_trim(line) == 0) cycle
        i = i + 1
        read(line, *) c0, c1, c2, c3, c4
        r_eff(i)   = c1   ! r_eff [cm] — toroidal-flux radius, NOT r_geom
        psi_tor(i) = c4   ! psi_tor [G cm^2]
    end do
    close(iunit)

end subroutine read_equil_mapping_pert

! ============================================================
! PRIVATE: Natural cubic spline (same algorithm as profile_data_mod)
! ============================================================
subroutine spline_natural_pert(n, x, y, y2)

    integer, intent(in) :: n
    real(dp), intent(in) :: x(n), y(n)
    real(dp), intent(out) :: y2(n)

    real(dp), allocatable :: u(:)
    real(dp) :: sig, p
    integer :: i

    allocate(u(n))

    y2(1) = 0.0_dp
    u(1) = 0.0_dp

    do i = 2, n - 1
        sig = (x(i) - x(i-1)) / (x(i+1) - x(i-1))
        p = sig * y2(i-1) + 2.0_dp
        y2(i) = (sig - 1.0_dp) / p
        u(i) = (6.0_dp * ((y(i+1) - y(i)) / (x(i+1) - x(i)) &
              - (y(i) - y(i-1)) / (x(i) - x(i-1))) &
              / (x(i+1) - x(i-1)) - sig * u(i-1)) / p
    end do

    y2(n) = 0.0_dp
    do i = n - 1, 1, -1
        y2(i) = y2(i) * y2(i+1) + u(i)
    end do

    deallocate(u)

end subroutine spline_natural_pert

! ============================================================
! PRIVATE: Cubic spline evaluation with derivative
! ============================================================
subroutine splint_pert(n, xa, ya, y2a, x, y, dydx)

    integer, intent(in) :: n
    real(dp), intent(in) :: xa(n), ya(n), y2a(n)
    real(dp), intent(in) :: x
    real(dp), intent(out) :: y, dydx

    integer :: klo, khi, k
    real(dp) :: h, a, b

    klo = 1
    khi = n
    do while (khi - klo > 1)
        k = (khi + klo) / 2
        if (xa(k) > x) then
            khi = k
        else
            klo = k
        end if
    end do

    h = xa(khi) - xa(klo)
    a = (xa(khi) - x) / h
    b = (x - xa(klo)) / h

    y = a * ya(klo) + b * ya(khi) + &
        ((a**3 - a) * y2a(klo) + (b**3 - b) * y2a(khi)) * h**2 / 6.0_dp

    dydx = (ya(khi) - ya(klo)) / h + &
           (-(3.0_dp * a**2 - 1.0_dp) * y2a(klo) + &
             (3.0_dp * b**2 - 1.0_dp) * y2a(khi)) * h / 6.0_dp

end subroutine splint_pert

end module perturbation_field_mod
