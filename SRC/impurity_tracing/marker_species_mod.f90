module marker_species_mod

! Charge states of one impurity element for markers whose species changes along the orbit
! (ionisation, recombination), e.g. Be0 ... Be4+ as in Harting et al., NME 42 (2025) 101887, Sec. 5.
!
! A marker carries its charge state iz (0 ... z_nuc). set_marker_charge sets GORILLA's per-thread
! species (tetra_physics_mod, threadprivate) for the pusher. change_charge_state applies a
! transition at a point of the orbit:
!   - the particle keeps its velocity (v_par, v_perp): the electron changes the ion momentum only
!     by O(m_e/m_ion), and kinetic energy and magnetic moment stay the same;
!   - the guiding centre jumps (boole_gc_jump): the particle position x = X + b x v/Omega is kept
!     and Omega = Z e B/(m c) changes, so X_new = X_old + (b x v_perp)(1/Omega_old - 1/Omega_new).
!     The gyrophase is not known in the guiding-centre model and is drawn at random.
! Transition rates [1/s] per charge state (ionisation iz -> iz+1, recombination iz -> iz-1):
!   - uniform: one value per charge state, given directly or from ADAS ADF11 files at the constant
!     background (adas_rates, adf11_mod);
!   - per cell: one value per charge state and cell, precomputed from cell values of n_e, T_e
!     (set_cell_rates, adas_cell_rates). A cell is any piece of the mesh on which the background is
!     constant (impurity tracing: a 2D triangle). The rate is then piecewise constant along the
!     orbit and the transition is sampled exactly: draw an optical depth tau = -ln(1-u)
!     (draw_optical_depth), add nu_cell * dt per piece of orbit, transition where the sum reaches
!     tau (time_to_event gives the time left in the current cell).
!     NOTE: for a background that varies inside a cell (e.g. linear per tetrahedron, as used by the
!     thermal-force collisions), the rates would instead have to be evaluated along the orbit at
!     the local n_e, T_e, with the sum of nu dt approximated per step (time-step error); not done.
! Neutrals (iz = 0) fly on straight lines (neutral_flight_mod); ion_to_neutral and neutral_to_ion convert
! between the guiding centre of the ion and the particle position and velocity of the neutral.
! Collision coefficients (collis_ions collis_init) are precomputed per charge state for a background
! of constant density and temperature.
! CGS units, temperatures in eV, masses in amu.

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

    private

    real(dp), parameter :: amu_g = 1.66053906660e-24_dp
    real(dp), parameter :: pi = 3.14159265358979323846264338327950288_dp

    type, public :: species_set_t
        character(len=8) :: element = ''
        integer :: z_nuc = 0
        real(dp) :: mass_amu = 0.0_dp, mass_g = 0.0_dp
        ! collision coefficients per background species (nb) and charge state (1:z_nuc)
        integer :: nb = 0
        real(dp) :: v0 = 0.0_dp
        real(dp), allocatable :: efcolf(:,:), velrat(:,:), enrat(:,:)
        ! transition rates per charge state (0:z_nuc) [1/s]
        real(dp), allocatable :: rate_ion(:), rate_rec(:)
        ! optional per-cell rates (0:z_nuc, ncell) [1/s]; used when allocated (cell index >= 1)
        real(dp), allocatable :: cell_rate_ion(:,:), cell_rate_rec(:,:)
        logical :: boole_gc_jump = .true.
    end type species_set_t

    public :: init_species_set, set_marker_charge, total_rate, next_event_time, choose_transition, &
              change_charge_state, field_direction, element_data, adas_rates, adas_cell_rates, &
              set_cell_rates, cell_rates, draw_optical_depth, time_to_event, neutral_to_ion, ion_to_neutral

contains

subroutine element_data(element, z_nuc, mass_amu)

    ! Nuclear charge and standard atomic mass of the supported elements

    character(len=*), intent(in) :: element
    integer, intent(out) :: z_nuc
    real(dp), intent(out) :: mass_amu

    select case (trim(adjustl(element)))
    case ('Be'); z_nuc = 4;  mass_amu = 9.0121831_dp
    case ('C');  z_nuc = 6;  mass_amu = 12.011_dp
    case ('N');  z_nuc = 7;  mass_amu = 14.007_dp
    case ('Ne'); z_nuc = 10; mass_amu = 20.1797_dp
    case ('Ar'); z_nuc = 18; mass_amu = 39.948_dp
    case ('W');  z_nuc = 74; mass_amu = 183.84_dp
    case default
        print *, 'marker species: unknown element ', trim(element), ' (known: Be, C, N, Ne, Ar, W)'
        error stop
    end select

end subroutine element_data

subroutine adas_rates(element, scd_file, acd_file, ne_cm3, te_ev, rate_ion, rate_rec)

    ! Uniform rates [1/s] per charge state 0 ... z_nuc from ADF11 files at one (n_e, T_e); see
    ! adas_cell_rates. An empty file name leaves the corresponding rates unchanged.

    character(len=*), intent(in) :: element, scd_file, acd_file
    real(dp), intent(in) :: ne_cm3, te_ev
    real(dp), intent(inout) :: rate_ion(0:), rate_rec(0:)
    real(dp), allocatable :: ri(:,:), rr(:,:)
    integer :: z_nuc
    real(dp) :: mass_amu

    call element_data(element, z_nuc, mass_amu)
    if (ubound(rate_ion, 1) < z_nuc .or. ubound(rate_rec, 1) < z_nuc) &
        error stop 'marker species: rate arrays shorter than the nuclear charge'
    allocate(ri(0:z_nuc, 1), rr(0:z_nuc, 1))
    ri(:, 1) = rate_ion(0:z_nuc)
    rr(:, 1) = rate_rec(0:z_nuc)
    call adas_cell_rates(element, scd_file, acd_file, [ne_cm3], [te_ev], ri, rr)
    rate_ion(0:z_nuc) = ri(:, 1)
    rate_rec(0:z_nuc) = rr(:, 1)

end subroutine adas_rates

subroutine adas_cell_rates(element, scd_file, acd_file, ne_cm3, te_ev, rate_ion, rate_rec)

    ! Rates [1/s] per charge state 0 ... z_nuc and cell from ADF11 files (each file read once):
    ! ionisation iz -> iz+1 from scd (block Z1 = iz+1), recombination iz -> iz-1 from acd
    ! (block Z1 = iz), rate = n_e <sigma v>(n_e, T_e) of the cell. An empty file name leaves the
    ! corresponding rates unchanged. Recombination to neutrals (iz = 1 -> 0) from block Z1 = 1.

    use adf11_mod, only: adf11_t, read_adf11, adf11_coef

    character(len=*), intent(in) :: element, scd_file, acd_file
    real(dp), intent(in) :: ne_cm3(:), te_ev(:)                  ! per cell
    real(dp), intent(inout) :: rate_ion(0:,:), rate_rec(0:,:)    ! (0:z_nuc, ncell)
    type(adf11_t) :: tab
    integer :: z_nuc, iz, ic
    real(dp) :: mass_amu

    call element_data(element, z_nuc, mass_amu)
    if (size(te_ev) /= size(ne_cm3) .or. size(rate_ion, 2) /= size(ne_cm3) .or. size(rate_rec, 2) /= size(ne_cm3) &
        .or. ubound(rate_ion, 1) < z_nuc .or. ubound(rate_rec, 1) < z_nuc) &
        error stop 'marker species: adas_cell_rates: inconsistent array sizes'
    if (len_trim(scd_file) > 0) then
        call read_adf11(trim(scd_file), tab)
        call check_table(tab)
        rate_ion(0:z_nuc, :) = 0.0_dp
        do ic = 1, size(ne_cm3)
            do iz = 0, z_nuc - 1
                rate_ion(iz, ic) = ne_cm3(ic)*adf11_coef(tab, iz + 1, ne_cm3(ic), te_ev(ic))
            end do
        end do
    end if
    if (len_trim(acd_file) > 0) then
        call read_adf11(trim(acd_file), tab)
        call check_table(tab)
        rate_rec(0:z_nuc, :) = 0.0_dp
        do ic = 1, size(ne_cm3)
            do iz = 1, z_nuc
                rate_rec(iz, ic) = ne_cm3(ic)*adf11_coef(tab, iz, ne_cm3(ic), te_ev(ic))
            end do
        end do
    end if

contains

    subroutine check_table(t)
        type(adf11_t), intent(in) :: t
        if (t%z_nuc /= z_nuc .or. t%z1_min /= 1 .or. t%z1_max /= z_nuc) then
            print *, 'marker species: ', trim(t%file), ' has Z = ', t%z_nuc, ', blocks Z1 = ', t%z1_min, ' ... ', &
                t%z1_max, '; element ', trim(element), ' needs Z = 1 ... ', z_nuc
            error stop
        end if
    end subroutine check_table

end subroutine adas_cell_rates

subroutine set_cell_rates(ss, rate_ion, rate_rec)

    ! Per-cell rates (0:z_nuc, ncell) [1/s]; same checks as the uniform rates of init_species_set

    type(species_set_t), intent(inout) :: ss
    real(dp), intent(in) :: rate_ion(0:,:), rate_rec(0:,:)

    if (ubound(rate_ion, 1) < ss%z_nuc .or. ubound(rate_rec, 1) < ss%z_nuc .or. size(rate_ion, 2) /= size(rate_rec, 2)) &
        error stop 'marker species: set_cell_rates: inconsistent array sizes'
    if (any(rate_ion(0:ss%z_nuc, :) < 0.0_dp) .or. any(rate_rec(0:ss%z_nuc, :) < 0.0_dp)) &
        error stop 'marker species: negative cell rate'
    if (any(rate_ion(ss%z_nuc, :) > 0.0_dp)) error stop 'marker species: ionisation of the bare nucleus'
    if (allocated(ss%cell_rate_ion)) deallocate(ss%cell_rate_ion, ss%cell_rate_rec)
    allocate(ss%cell_rate_ion(0:ss%z_nuc, size(rate_ion, 2)), ss%cell_rate_rec(0:ss%z_nuc, size(rate_rec, 2)))
    ss%cell_rate_ion(:,:) = rate_ion(0:ss%z_nuc, :)          ! explicit bounds: charge state 0 ... z_nuc
    ss%cell_rate_rec(:,:) = rate_rec(0:ss%z_nuc, :)

end subroutine set_cell_rates

pure subroutine cell_rates(ss, iz, icell, nu_ion, nu_rec)

    ! Rates out of charge state iz in cell icell (uniform rates if no cell rates are set or icell < 1)

    type(species_set_t), intent(in) :: ss
    integer, intent(in) :: iz, icell
    real(dp), intent(out) :: nu_ion, nu_rec

    if (allocated(ss%cell_rate_ion) .and. icell >= 1) then
        nu_ion = ss%cell_rate_ion(iz, icell)
        nu_rec = ss%cell_rate_rec(iz, icell)
    else
        nu_ion = ss%rate_ion(iz)
        nu_rec = ss%rate_rec(iz)
    end if

end subroutine cell_rates

function draw_optical_depth(rng) result(tau)

    ! Optical depth to the next transition, exponentially distributed with mean 1

    use thermal_force_rng_mod, only: rng_t, rng_uniform

    type(rng_t), intent(inout) :: rng
    real(dp) :: tau

    tau = -log(1.0_dp - rng_uniform(rng))

end function draw_optical_depth

pure function time_to_event(nu, tau_acc, tau) result(dt)

    ! Time until the accumulated optical depth tau_acc reaches tau at constant rate nu;
    ! huge() for nu = 0

    real(dp), intent(in) :: nu, tau_acc, tau
    real(dp) :: dt

    dt = huge(1.0_dp)
    if (nu > 0.0_dp) dt = max(0.0_dp, tau - tau_acc)/nu

end function time_to_event

subroutine init_species_set(ss, element, mb, zb, db, tb, temperature_ev, rate_ion, rate_rec, boole_gc_jump)

    ! mb, zb, db, tb: background species (mass [g], charge number, density [cm^-3], temperature [eV]);
    ! temperature_ev: velocity normalisation v0 = sqrt(2 T/m); rate_ion, rate_rec: rates [1/s]
    ! indexed by charge state 0 ... (at least z_nuc)

    use collis_ions, only: collis_init

    type(species_set_t), intent(out) :: ss
    character(len=*), intent(in) :: element
    real(dp), intent(in) :: mb(:), zb(:), db(:), tb(:), temperature_ev
    real(dp), intent(in) :: rate_ion(0:), rate_rec(0:)
    logical, intent(in) :: boole_gc_jump
    integer :: iz
    real(dp) :: v0

    ss%element = element
    call element_data(element, ss%z_nuc, ss%mass_amu)
    ss%mass_g = ss%mass_amu*amu_g
    ss%nb = size(mb)
    ss%boole_gc_jump = boole_gc_jump
    if (ubound(rate_ion, 1) < ss%z_nuc .or. ubound(rate_rec, 1) < ss%z_nuc) &
        error stop 'marker species: rate arrays shorter than the nuclear charge'
    allocate(ss%rate_ion(0:ss%z_nuc), ss%rate_rec(0:ss%z_nuc))
    ss%rate_ion = rate_ion(0:ss%z_nuc)
    ss%rate_rec = rate_rec(0:ss%z_nuc)
    if (any(ss%rate_ion < 0.0_dp) .or. any(ss%rate_rec < 0.0_dp)) error stop 'marker species: negative rate'
    if (ss%rate_ion(ss%z_nuc) > 0.0_dp) error stop 'marker species: ionisation of the bare nucleus'

    allocate(ss%efcolf(ss%nb, ss%z_nuc), ss%velrat(ss%nb, ss%z_nuc), ss%enrat(ss%nb, ss%z_nuc))
    do iz = 1, ss%z_nuc
        call collis_init(ss%mass_g, real(iz, dp), mb, zb, db, tb, temperature_ev, v0, &
                         ss%efcolf(:, iz), ss%velrat(:, iz), ss%enrat(:, iz))
    end do
    ss%v0 = v0

end subroutine init_species_set

subroutine set_marker_charge(ss, iz)

    ! GORILLA species of the calling thread for a marker in charge state iz

    use constants, only: echarge
    use tetra_physics_mod, only: set_particle_species

    type(species_set_t), intent(in) :: ss
    integer, intent(in) :: iz

    if (iz < 1 .or. iz > ss%z_nuc) error stop 'marker species: charge state outside 1 ... z_nuc (neutrals: neutral_flight_mod)'
    call set_particle_species(ss%mass_g, real(iz, dp)*echarge)

end subroutine set_marker_charge

pure function total_rate(ss, iz) result(nu)

    type(species_set_t), intent(in) :: ss
    integer, intent(in) :: iz
    real(dp) :: nu

    nu = ss%rate_ion(iz) + ss%rate_rec(iz)

end function total_rate

function next_event_time(ss, iz, t, rng) result(t_event)

    ! Absolute time of the next transition out of iz for a marker at time t: exponential waiting
    ! time with the total rate of iz. huge() if iz has no transitions (no random number drawn then).

    use thermal_force_rng_mod, only: rng_t, rng_uniform

    type(species_set_t), intent(in) :: ss
    integer, intent(in) :: iz
    real(dp), intent(in) :: t
    type(rng_t), intent(inout) :: rng
    real(dp) :: t_event, nu

    t_event = huge(1.0_dp)
    nu = total_rate(ss, iz)
    if (nu > 0.0_dp) t_event = t - log(1.0_dp - rng_uniform(rng))/nu

end function next_event_time

pure function choose_transition(ss, iz, u, icell) result(iz_new)

    ! Next charge state for a transition out of iz; u uniform in [0, 1); icell: cell of the
    ! transition (per-cell rates), absent or < 1 for the uniform rates

    type(species_set_t), intent(in) :: ss
    integer, intent(in) :: iz
    real(dp), intent(in) :: u
    integer, intent(in), optional :: icell
    integer :: iz_new
    real(dp) :: nu_ion, nu_rec

    if (present(icell)) then
        call cell_rates(ss, iz, icell, nu_ion, nu_rec)
    else
        call cell_rates(ss, iz, 0, nu_ion, nu_rec)
    end if
    if (u*(nu_ion + nu_rec) < nu_ion) then
        iz_new = iz + 1
    else
        iz_new = iz - 1
    end if

end function choose_transition

subroutine field_direction(x, ind_tetr, b_hat, bmod)

    ! Unit vector of B at x (R, phi, Z) in physical components (R, phi, Z), and |B| [G], from the
    ! linear interpolation of the covariant h_k in tetrahedron ind_tetr (coord_system = 1)

    use tetra_physics_mod, only: tetra_physics
    use supporting_functions_mod, only: bmod_func

    real(dp), intent(in) :: x(3)
    integer, intent(in) :: ind_tetr
    real(dp), intent(out) :: b_hat(3), bmod
    real(dp) :: z(3), h(3)

    z = x - tetra_physics(ind_tetr)%x1
    h(1) = tetra_physics(ind_tetr)%h1_1 + sum(tetra_physics(ind_tetr)%gh1*z)
    h(2) = (tetra_physics(ind_tetr)%h2_1 + sum(tetra_physics(ind_tetr)%gh2*z))/x(1)
    h(3) = tetra_physics(ind_tetr)%h3_1 + sum(tetra_physics(ind_tetr)%gh3*z)
    b_hat = h/norm2(h)
    bmod = bmod_func(z, ind_tetr)

end subroutine field_direction

subroutine change_charge_state(ss, iz_old, iz_new, phase_u, x, ind_tetr, iface, vpar, vperp, jump)

    ! Transition iz_old -> iz_new at the guiding-centre position x. Velocity kept; guiding-centre
    ! jump (if ss%boole_gc_jump) with random gyrophase 2 pi phase_u; GORILLA species updated.
    ! ind_tetr = -1 on return if the jump left the mesh. jump: physical displacement (R, phi, Z) [cm].

    use constants, only: clight, echarge
    use anomalous_transport_displacement_mod, only: displace_by_straight_line

    type(species_set_t), intent(in) :: ss
    integer, intent(in) :: iz_old, iz_new
    real(dp), intent(in) :: phase_u, vpar, vperp
    real(dp), intent(inout) :: x(3)
    integer, intent(inout) :: ind_tetr, iface
    real(dp), intent(out) :: jump(3)
    real(dp) :: b_hat(3), bmod, e1(3), e2(3), alpha, len

    jump = 0.0_dp
    if (ss%boole_gc_jump .and. iz_old > 0 .and. iz_new > 0) then
        call field_direction(x, ind_tetr, b_hat, bmod)
        call perp_basis(b_hat, e1, e2)
        alpha = 2.0_dp*pi*phase_u
        ! |b x v_perp| (1/Omega_old - 1/Omega_new), Omega = Z e B/(m c)
        len = vperp*ss%mass_g*clight/(echarge*bmod)*(1.0_dp/real(iz_old, dp) - 1.0_dp/real(iz_new, dp))
        jump = len*(cos(alpha)*e1 + sin(alpha)*e2)
        call displace_by_straight_line(x, ind_tetr, iface, [jump(1), jump(2)/x(1), jump(3)], vpar, vperp)
    end if
    if (iz_new >= 1) call set_marker_charge(ss, iz_new)

end subroutine change_charge_state

subroutine neutral_to_ion(ss, x, ind_tetr, iface, vcart, vpar, vperp, jump)

    ! Neutral -> charge 1 at the particle position x (R, phi, Z) with Cartesian velocity vcart: the ion
    ! keeps the velocity, v_par = v.b, v_perp = |v - v_par b|; guiding centre X = x - (b x v)/Omega,
    ! Omega = e B/(m c) (gyrophase known). x, ind_tetr, iface: guiding centre on return
    ! (ind_tetr = -1 if outside the mesh). jump: X - x, physical components (R, phi, Z) [cm].

    use constants, only: clight, echarge
    use anomalous_transport_displacement_mod, only: displace_by_straight_line

    type(species_set_t), intent(in) :: ss
    real(dp), intent(inout) :: x(3)
    integer, intent(inout) :: ind_tetr, iface
    real(dp), intent(in) :: vcart(3)
    real(dp), intent(out) :: vpar, vperp, jump(3)
    real(dp) :: b_hat(3), bmod, vloc(3), c, s, bxv(3)

    c = cos(x(2)); s = sin(x(2))
    vloc = [c*vcart(1) + s*vcart(2), -s*vcart(1) + c*vcart(2), vcart(3)]    ! along e_R, e_phi, e_Z at x
    call field_direction(x, ind_tetr, b_hat, bmod)
    vpar = dot_product(vloc, b_hat)
    vperp = norm2(vloc - vpar*b_hat)
    bxv = cross(b_hat, vloc)
    jump = -bxv*ss%mass_g*clight/(echarge*bmod)
    call displace_by_straight_line(x, ind_tetr, iface, [jump(1), jump(2)/x(1), jump(3)], vpar, vperp)
    call set_marker_charge(ss, 1)

end subroutine neutral_to_ion

subroutine ion_to_neutral(ss, iz_old, phase_u, x, ind_tetr, iface, vpar, vperp, vcart, jump)

    ! Charge iz_old -> neutral at the guiding centre x (R, phi, Z) of the ion: random gyrophase
    ! alpha = 2 pi phase_u, v_perp vector = v_perp (cos alpha e1 + sin alpha e2) (e1, e2 perpendicular to
    ! b), particle position x + (b x v)/Omega, Omega = iz_old e B/(m c). The neutral keeps the
    ! velocity: vcart (Cartesian). x, ind_tetr, iface: particle position on return (ind_tetr = -1 if
    ! outside the mesh). jump: particle - guiding centre, physical components (R, phi, Z) [cm].

    use constants, only: clight, echarge
    use anomalous_transport_displacement_mod, only: displace_by_straight_line

    type(species_set_t), intent(in) :: ss
    integer, intent(in) :: iz_old
    real(dp), intent(in) :: phase_u, vpar, vperp
    real(dp), intent(inout) :: x(3)
    integer, intent(inout) :: ind_tetr, iface
    real(dp), intent(out) :: vcart(3), jump(3)
    real(dp) :: b_hat(3), bmod, e1(3), e2(3), alpha, vloc(3), c, s, vperp_vec(3)

    call field_direction(x, ind_tetr, b_hat, bmod)
    call perp_basis(b_hat, e1, e2)
    alpha = 2.0_dp*pi*phase_u
    vperp_vec = vperp*(cos(alpha)*e1 + sin(alpha)*e2)
    vloc = vpar*b_hat + vperp_vec
    jump = cross(b_hat, vperp_vec)*ss%mass_g*clight/(real(iz_old, dp)*echarge*bmod)
    ! Cartesian velocity: vloc is given in the local basis at the guiding centre
    c = cos(x(2)); s = sin(x(2))
    vcart = [c*vloc(1) - s*vloc(2), s*vloc(1) + c*vloc(2), vloc(3)]
    call displace_by_straight_line(x, ind_tetr, iface, [jump(1), jump(2)/x(1), jump(3)], vpar, vperp)

end subroutine ion_to_neutral

pure function cross(a, b) result(c)

    real(dp), intent(in) :: a(3), b(3)
    real(dp) :: c(3)

    c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]

end function cross

pure subroutine perp_basis(b, e1, e2)

    ! Orthonormal e1, e2 perpendicular to the unit vector b, (e1, e2, b) right-handed

    real(dp), intent(in) :: b(3)
    real(dp), intent(out) :: e1(3), e2(3)
    real(dp) :: a(3)

    a = [1.0_dp, 0.0_dp, 0.0_dp]
    if (abs(b(1)) > 0.9_dp) a = [0.0_dp, 0.0_dp, 1.0_dp]
    e1 = a - dot_product(a, b)*b
    e1 = e1/norm2(e1)
    e2 = [b(2)*e1(3) - b(3)*e1(2), b(3)*e1(1) - b(1)*e1(3), b(1)*e1(2) - b(2)*e1(1)]

end subroutine perp_basis

end module marker_species_mod
