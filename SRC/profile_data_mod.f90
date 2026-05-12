!
! profile_data_mod.f90
!
! Loads plasma profiles (n, Te, Ti, Er) from two-column text files in KAMEL/CGS
! format (r_eff [cm] vs. value) and the equil_r_q_psi.dat coordinate mapping.
! Builds cubic splines in the GORILLA flux coordinate s (= psi_tor / psi_tor_edge)
! so that profiles and their s-derivatives can be evaluated along particle orbits.
!
! Units throughout:
!   - r_eff: cm
!   - n: 1/cm^3
!   - Te, Ti: eV  (internally converted to erg via ev2erg where needed)
!   - Er: statV/cm
!   - Phi_0: statV (electrostatic potential, integrated from Er)
!   - s: dimensionless (normalised toroidal flux, GORILLA's x(1))
!   - psi_pol, psi_tor: Gauss*cm^2
!
module profile_data_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

    private
    public :: load_profiles, eval_profiles, cleanup_profiles
    public :: eval_q, get_psi_tor_edge
    public :: eval_s_from_psi_pol, eval_ds_dreff
    public :: profile_values_t

    ! Result type returned by eval_profiles
    type :: profile_values_t
        real(dp) :: n_e        ! electron density [1/cm^3]
        real(dp) :: Te         ! electron temperature [eV]
        real(dp) :: Ti         ! ion temperature [eV]
        real(dp) :: Phi0       ! electrostatic potential [statV]
        real(dp) :: dn_ds      ! d(n_e)/ds
        real(dp) :: dTe_ds     ! d(Te)/ds
        real(dp) :: dTi_ds     ! d(Ti)/ds
        real(dp) :: dPhi0_ds   ! d(Phi_0)/ds
        real(dp) :: dlnn_ds    ! d(ln n_e)/ds
        real(dp) :: dlnTe_ds   ! d(ln Te)/ds
        real(dp) :: dlnTi_ds   ! d(ln Ti)/ds
    end type profile_values_t

    ! Number of grid points after mapping to s
    integer :: ns
    ! s grid and spline data (second derivatives for cubic spline)
    real(dp), allocatable :: s_grid(:)
    real(dp), allocatable :: ne_spl(:),  ne_dd(:)
    real(dp), allocatable :: Te_spl(:),  Te_dd(:)
    real(dp), allocatable :: Ti_spl(:),  Ti_dd(:)
    real(dp), allocatable :: Phi_spl(:), Phi_dd(:)

    ! q(s) spline built from flux_functions.dat. Needed to convert
    ! d/ds (where s = psi_tor/psi_tor_edge) to d/dpsi_pol via
    !     d/dpsi_pol = (q(s)/psi_tor_edge) d/ds.
    real(dp), allocatable :: q_spl(:),   q_dd(:)
    real(dp) :: psi_tor_edge_stored = 0.0_dp

    ! r_eff(s) spline built from flux_functions.dat (column "r" in cm).
    ! Used by eval_ds_dreff to convert delta_B_r (input as Gauss) to the
    ! contravariant component delta_B^s = (ds/dr_eff) * delta_B_r [Gauss/cm].
    real(dp), allocatable :: r_eff_spl(:), r_eff_dd(:)

    ! Map psi_pol -> s_true = psi_tor/psi_tor_edge so that an arbitrary
    ! poloidal-flux value (e.g. tetra_physics(...)%Aphi1 + gAphi . z_cell
    ! at a particle's position) can be converted to true s for profile
    ! lookups and s-window cuts.  Built in load_profiles, sorted in
    ! ascending psi_pol so binary search works.
    real(dp), allocatable :: psi_pol_sorted(:)
    real(dp), allocatable :: s_of_psi_pol_spl(:), s_of_psi_pol_dd(:)
    integer :: ns_psipol = 0

    logical :: profiles_loaded = .false.

contains

! ============================================================
! Load all profiles and build splines in s-coordinate
! ============================================================
subroutine load_profiles(profile_dir, equil_mapping_file)

    character(len=*), intent(in) :: profile_dir
    character(len=*), intent(in) :: equil_mapping_file

    ! Local arrays for equil mapping
    integer :: n_equil
    real(dp), allocatable :: r_equil(:), q_equil(:), psi_pol_equil(:), psi_tor_equil(:)
    ! Local arrays for raw profiles
    integer :: n_raw
    real(dp), allocatable :: r_raw(:), val_raw(:)
    ! Mapped profiles on s grid
    real(dp), allocatable :: ne_mapped(:), Te_mapped(:), Ti_mapped(:), Er_mapped(:), Phi_mapped(:)
    real(dp) :: psi_tor_edge, ds
    integer :: i

    print *, ''
    print *, 'Loading plasma profiles...'
    print *, '  Profile directory: ', trim(profile_dir)
    print *, '  Equil mapping:    ', trim(equil_mapping_file)

    ! --- Load equil_r_q_psi.dat ---
    call read_equil_mapping(equil_mapping_file, n_equil, r_equil, q_equil, &
                            psi_pol_equil, psi_tor_equil)

    ! Compute s = psi_tor / psi_tor_edge (normalised toroidal flux)
    psi_tor_edge = psi_tor_equil(n_equil)
    psi_tor_edge_stored = psi_tor_edge
    ns = n_equil
    allocate(s_grid(ns))
    if (abs(psi_tor_edge) > 0.0_dp) then
        s_grid = psi_tor_equil / psi_tor_edge
    else
        ! Fallback: uniform grid
        do i = 1, ns
            s_grid(i) = real(i - 1, dp) / real(ns - 1, dp)
        end do
    end if
    ! Ensure s_grid(1) = 0 and s_grid(ns) = 1
    s_grid(1) = 0.0_dp
    s_grid(ns) = 1.0_dp

    ! Store q(s) for the d/ds -> d/dpsi_pol chain rule used by the
    ! delta-f weight evolution.
    allocate(q_spl(ns), q_dd(ns))
    q_spl = q_equil
    call spline_natural(ns, s_grid, q_spl, q_dd)

    ! Store r_eff(s) [cm] so eval_ds_dreff can return the geometric factor
    ! that converts delta_B_r (Gauss) to delta_B^s (Gauss/cm).
    allocate(r_eff_spl(ns), r_eff_dd(ns))
    r_eff_spl = r_equil
    call spline_natural(ns, s_grid, r_eff_spl, r_eff_dd)

    ! Map psi_pol -> true s so the applet can convert local poloidal
    ! flux (from Aphi1 + gAphi . z_cell) to physical s for profile
    ! lookups and window cuts.  flux_functions.dat may have psi_pol
    ! decreasing with r; reverse if needed so the array is monotonically
    ! increasing for binary search.
    ns_psipol = ns
    if (allocated(psi_pol_sorted))   deallocate(psi_pol_sorted)
    if (allocated(s_of_psi_pol_spl)) deallocate(s_of_psi_pol_spl)
    if (allocated(s_of_psi_pol_dd))  deallocate(s_of_psi_pol_dd)
    allocate(psi_pol_sorted(ns_psipol))
    allocate(s_of_psi_pol_spl(ns_psipol))
    allocate(s_of_psi_pol_dd(ns_psipol))
    if (psi_pol_equil(ns_psipol) >= psi_pol_equil(1)) then
        psi_pol_sorted   = psi_pol_equil
        s_of_psi_pol_spl = s_grid
    else
        do i = 1, ns_psipol
            psi_pol_sorted(i)   = psi_pol_equil(ns_psipol - i + 1)
            s_of_psi_pol_spl(i) = s_grid(ns_psipol - i + 1)
        end do
    end if
    call spline_natural(ns_psipol, psi_pol_sorted, s_of_psi_pol_spl, s_of_psi_pol_dd)

    ! --- Load and map each profile onto s_grid via r_eff ---
    allocate(ne_mapped(ns), Te_mapped(ns), Ti_mapped(ns), Er_mapped(ns), Phi_mapped(ns))

    call read_and_map_profile(trim(profile_dir)//'/n.dat', r_equil, ns, ne_mapped)
    call read_and_map_profile(trim(profile_dir)//'/Te.dat', r_equil, ns, Te_mapped)
    call read_and_map_profile(trim(profile_dir)//'/Ti.dat', r_equil, ns, Ti_mapped)
    call read_and_map_profile(trim(profile_dir)//'/Er.dat', r_equil, ns, Er_mapped)

    ! --- Compute Phi_0(s) by integrating Er: Phi_0(s) = -integral Er dr ---
    ! Er = -dPhi/dr, so Phi(r) = Phi(0) - integral_0^r Er(r') dr'
    ! We set Phi(0) = 0
    Phi_mapped(1) = 0.0_dp
    do i = 2, ns
        ds = r_equil(i) - r_equil(i-1)
        Phi_mapped(i) = Phi_mapped(i-1) - Er_mapped(i) * ds
    end do

    ! --- Build cubic splines in s ---
    allocate(ne_spl(ns), ne_dd(ns))
    allocate(Te_spl(ns), Te_dd(ns))
    allocate(Ti_spl(ns), Ti_dd(ns))
    allocate(Phi_spl(ns), Phi_dd(ns))

    ne_spl  = ne_mapped
    Te_spl  = Te_mapped
    Ti_spl  = Ti_mapped
    Phi_spl = Phi_mapped

    call spline_natural(ns, s_grid, ne_spl, ne_dd)
    call spline_natural(ns, s_grid, Te_spl, Te_dd)
    call spline_natural(ns, s_grid, Ti_spl, Ti_dd)
    call spline_natural(ns, s_grid, Phi_spl, Phi_dd)

    profiles_loaded = .true.

    print *, '  Profiles loaded and splined on ', ns, ' s-grid points'
    print *, '  n_e(0) = ', ne_spl(1), ' 1/cm^3'
    print *, '  Te(0)  = ', Te_spl(1), ' eV'
    print *, '  Ti(0)  = ', Ti_spl(1), ' eV'
    print *, ''

    deallocate(r_equil, q_equil, psi_pol_equil, psi_tor_equil)
    deallocate(ne_mapped, Te_mapped, Ti_mapped, Er_mapped, Phi_mapped)

end subroutine load_profiles

! ============================================================
! Evaluate all profiles and their s-derivatives at a given s
! ============================================================
subroutine eval_profiles(s_val, pv)

    real(dp), intent(in)              :: s_val
    type(profile_values_t), intent(out) :: pv

    real(dp) :: s_clamped, ne_val, Te_val, Ti_val, Phi_val
    real(dp) :: dne, dTe, dTi, dPhi

    s_clamped = max(0.0_dp, min(1.0_dp, s_val))

    call splint_with_deriv(ns, s_grid, ne_spl,  ne_dd,  s_clamped, ne_val,  dne)
    call splint_with_deriv(ns, s_grid, Te_spl,  Te_dd,  s_clamped, Te_val,  dTe)
    call splint_with_deriv(ns, s_grid, Ti_spl,  Ti_dd,  s_clamped, Ti_val,  dTi)
    call splint_with_deriv(ns, s_grid, Phi_spl, Phi_dd, s_clamped, Phi_val, dPhi)

    ! Clamp to positive values for safety
    ne_val = max(ne_val, 1.0d0)
    Te_val = max(Te_val, 1.0d0)
    Ti_val = max(Ti_val, 1.0d0)

    pv%n_e     = ne_val
    pv%Te      = Te_val
    pv%Ti      = Ti_val
    pv%Phi0    = Phi_val
    pv%dn_ds   = dne
    pv%dTe_ds  = dTe
    pv%dTi_ds  = dTi
    pv%dPhi0_ds = dPhi

    ! Logarithmic derivatives
    pv%dlnn_ds  = dne / ne_val
    pv%dlnTe_ds = dTe / Te_val
    pv%dlnTi_ds = dTi / Ti_val

end subroutine eval_profiles

! ============================================================
! Cleanup allocated memory
! ============================================================
subroutine cleanup_profiles()

    if (allocated(s_grid))  deallocate(s_grid)
    if (allocated(ne_spl))  deallocate(ne_spl)
    if (allocated(ne_dd))   deallocate(ne_dd)
    if (allocated(Te_spl))  deallocate(Te_spl)
    if (allocated(Te_dd))   deallocate(Te_dd)
    if (allocated(Ti_spl))  deallocate(Ti_spl)
    if (allocated(Ti_dd))   deallocate(Ti_dd)
    if (allocated(Phi_spl)) deallocate(Phi_spl)
    if (allocated(Phi_dd))  deallocate(Phi_dd)
    if (allocated(q_spl))   deallocate(q_spl)
    if (allocated(q_dd))    deallocate(q_dd)
    if (allocated(r_eff_spl)) deallocate(r_eff_spl)
    if (allocated(r_eff_dd))  deallocate(r_eff_dd)
    if (allocated(psi_pol_sorted))   deallocate(psi_pol_sorted)
    if (allocated(s_of_psi_pol_spl)) deallocate(s_of_psi_pol_spl)
    if (allocated(s_of_psi_pol_dd))  deallocate(s_of_psi_pol_dd)
    psi_tor_edge_stored = 0.0_dp
    ns_psipol = 0
    profiles_loaded = .false.

end subroutine cleanup_profiles

! ============================================================
! Accessors: safety factor q(s) and psi_tor_edge, used by callers
! that need to convert between d/ds and d/dpsi_pol.
! ============================================================
real(dp) function eval_q(s_val)
    real(dp), intent(in) :: s_val
    real(dp) :: s_clamped, q_val, dq
    s_clamped = max(0.0_dp, min(1.0_dp, s_val))
    call splint_with_deriv(ns, s_grid, q_spl, q_dd, s_clamped, q_val, dq)
    eval_q = q_val
end function eval_q

real(dp) function get_psi_tor_edge()
    get_psi_tor_edge = psi_tor_edge_stored
end function get_psi_tor_edge

! ============================================================
! ds / dr_eff at a given s, in 1/cm.  Builds on the r_eff(s) spline
! from flux_functions.dat: derivative dr_eff/ds is taken from the spline
! and inverted.  Used by eval_delta_B_s to convert delta_B_r [Gauss] to
! delta_B^s [Gauss/cm] = (ds/dr_eff) * delta_B_r.
! ============================================================
real(dp) function eval_ds_dreff(s_val) result(ds_dreff)
    real(dp), intent(in) :: s_val
    real(dp) :: s_clamped, r_val, dr_ds

    if (ns < 2 .or. .not. profiles_loaded) then
        ds_dreff = 0.0_dp
        return
    end if
    s_clamped = max(0.0_dp, min(1.0_dp, s_val))
    call splint_with_deriv(ns, s_grid, r_eff_spl, r_eff_dd, s_clamped, r_val, dr_ds)
    if (dr_ds > 0.0_dp) then
        ds_dreff = 1.0_dp / dr_ds
    else
        ds_dreff = 0.0_dp
    end if
end function eval_ds_dreff

! ============================================================
! Convert a local poloidal flux value (e.g. as evaluated inside a
! GORILLA tetrahedron) to true normalised toroidal flux s. Clamped to
! [0, 1] so values outside the equilibrium support are pinned to the
! relevant boundary.
! ============================================================
real(dp) function eval_s_from_psi_pol(psi_pol_val) result(s_val)
    real(dp), intent(in) :: psi_pol_val
    real(dp) :: psi_lo, psi_hi, dummy

    if (ns_psipol < 2 .or. .not. profiles_loaded) then
        s_val = 0.0_dp
        return
    end if
    psi_lo = psi_pol_sorted(1)
    psi_hi = psi_pol_sorted(ns_psipol)

    if (psi_pol_val <= psi_lo) then
        s_val = s_of_psi_pol_spl(1)
    else if (psi_pol_val >= psi_hi) then
        s_val = s_of_psi_pol_spl(ns_psipol)
    else
        call splint_with_deriv(ns_psipol, psi_pol_sorted, s_of_psi_pol_spl, &
                               s_of_psi_pol_dd, psi_pol_val, s_val, dummy)
    end if
    s_val = max(0.0_dp, min(1.0_dp, s_val))
end function eval_s_from_psi_pol

! ============================================================
! PRIVATE: Read equil_r_q_psi.dat mapping file
! ============================================================
subroutine read_equil_mapping(filename, n, r_eff, q, psi_pol, psi_tor)

    ! Reads GORILLA's flux_functions.dat, which has a header line starting
    ! with " #" followed by data rows of five columns:
    !
    !   R_beg [cm]   r_eff [cm]   q   psi_pol [Gauss*cm^2]   psi_tor [Gauss*cm^2]
    !
    ! Only r_eff, q, psi_pol, psi_tor are returned.

    character(len=*), intent(in) :: filename
    integer, intent(out) :: n
    real(dp), allocatable, intent(out) :: r_eff(:), q(:), psi_pol(:), psi_tor(:)

    integer :: iunit, ios, i
    character(len=512) :: line
    real(dp) :: R_beg, r_val, q_val, psi_p, psi_t

    ! First pass: count data lines
    iunit = 42
    open(unit=iunit, file=filename, status='old', action='read', iostat=ios)
    if (ios /= 0) then
        print *, 'ERROR: Cannot open ', trim(filename)
        stop
    end if

    n = 0
    do
        read(iunit, '(A)', iostat=ios) line
        if (ios /= 0) exit
        if (line(1:2) == ' #' .or. line(1:1) == '#' .or. len_trim(line) == 0) cycle
        n = n + 1
    end do
    close(iunit)

    allocate(r_eff(n), q(n), psi_pol(n), psi_tor(n))

    ! Second pass: read data
    open(unit=iunit, file=filename, status='old', action='read')
    i = 0
    do
        read(iunit, '(A)', iostat=ios) line
        if (ios /= 0) exit
        if (line(1:2) == ' #' .or. line(1:1) == '#' .or. len_trim(line) == 0) cycle
        i = i + 1
        read(line, *) R_beg, r_val, q_val, psi_p, psi_t
        r_eff(i) = r_val
        q(i) = q_val
        psi_pol(i) = psi_p
        psi_tor(i) = psi_t
    end do
    close(iunit)

    print *, '  Loaded equil mapping: ', n, ' points, r_max = ', r_eff(n), ' cm'

end subroutine read_equil_mapping

! ============================================================
! PRIVATE: Read a two-column profile file and interpolate onto
!          the r_eff grid from the equil mapping
! ============================================================
subroutine read_and_map_profile(filename, r_target, n_target, val_mapped)

    character(len=*), intent(in) :: filename
    real(dp), intent(in) :: r_target(:)
    integer, intent(in) :: n_target
    real(dp), intent(out) :: val_mapped(:)

    integer :: iunit, ios, n_raw, i, j
    character(len=512) :: line
    real(dp), allocatable :: r_raw(:), v_raw(:)
    real(dp) :: frac

    ! First pass: count data lines
    iunit = 43
    open(unit=iunit, file=filename, status='old', action='read', iostat=ios)
    if (ios /= 0) then
        print *, '  WARNING: Cannot open ', trim(filename), ' -- using zeros'
        val_mapped = 0.0_dp
        return
    end if

    n_raw = 0
    do
        read(iunit, '(A)', iostat=ios) line
        if (ios /= 0) exit
        if (line(1:1) == '#' .or. len_trim(line) == 0) cycle
        n_raw = n_raw + 1
    end do
    close(iunit)

    allocate(r_raw(n_raw), v_raw(n_raw))

    ! Second pass: read data
    open(unit=iunit, file=filename, status='old', action='read')
    i = 0
    do
        read(iunit, '(A)', iostat=ios) line
        if (ios /= 0) exit
        if (line(1:1) == '#' .or. len_trim(line) == 0) cycle
        i = i + 1
        read(line, *) r_raw(i), v_raw(i)
    end do
    close(iunit)

    ! Linear interpolation onto r_target grid
    do i = 1, n_target
        if (r_target(i) <= r_raw(1)) then
            val_mapped(i) = v_raw(1)
        else if (r_target(i) >= r_raw(n_raw)) then
            val_mapped(i) = v_raw(n_raw)
        else
            ! Find bracket
            do j = 1, n_raw - 1
                if (r_raw(j+1) >= r_target(i)) exit
            end do
            frac = (r_target(i) - r_raw(j)) / (r_raw(j+1) - r_raw(j))
            val_mapped(i) = v_raw(j) + frac * (v_raw(j+1) - v_raw(j))
        end if
    end do

    print *, '  Loaded profile: ', trim(filename), ' (', n_raw, ' points)'

    deallocate(r_raw, v_raw)

end subroutine read_and_map_profile

! ============================================================
! PRIVATE: Natural cubic spline setup (second derivatives)
! Adapted from Numerical Recipes spline routine.
! ============================================================
subroutine spline_natural(n, x, y, y2)

    integer, intent(in) :: n
    real(dp), intent(in) :: x(n), y(n)
    real(dp), intent(out) :: y2(n)

    real(dp), allocatable :: u(:)
    real(dp) :: sig, p
    integer :: i

    allocate(u(n))

    ! Natural spline: y2(1) = y2(n) = 0
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

end subroutine spline_natural

! ============================================================
! PRIVATE: Cubic spline evaluation with derivative
! ============================================================
subroutine splint_with_deriv(n, xa, ya, y2a, x, y, dydx)

    integer, intent(in) :: n
    real(dp), intent(in) :: xa(n), ya(n), y2a(n)
    real(dp), intent(in) :: x
    real(dp), intent(out) :: y, dydx

    integer :: klo, khi, k
    real(dp) :: h, a, b

    ! Bisection to find interval
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

    ! Spline value
    y = a * ya(klo) + b * ya(khi) + &
        ((a**3 - a) * y2a(klo) + (b**3 - b) * y2a(khi)) * h**2 / 6.0_dp

    ! Spline derivative
    dydx = (ya(khi) - ya(klo)) / h + &
           (-(3.0_dp * a**2 - 1.0_dp) * y2a(klo) + &
             (3.0_dp * b**2 - 1.0_dp) * y2a(khi)) * h / 6.0_dp

end subroutine splint_with_deriv

end module profile_data_mod
