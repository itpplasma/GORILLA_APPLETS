module impurity_tracing_mod

! Impurity tracing (i_option = 17), input impurity_tracing.inp.
!
! Test ions (fixed species, e.g. W5+, or an element with charge-state changes) are followed as GORILLA guiding
! centres in the edge region of an extended EFIT mesh (tetra_grid.inp:
! boole_extend_to_wall, s_edge_region; GORILLA tetra_in_edge_region) until they
!   - hit the wall (leave the mesh)                       -> status 1
!   - cross the inner region boundary s = s_edge_region   -> status 2
!   - reach the tracing limit t_max_s                     -> status 0
!   - are removed by the pusher                           -> status 3
!   - start outside the region                            -> status -1
! Collisions: GORILLA's Maxwellian Coulomb operator (collis_ions stost, full
! operator: pitch-angle and energy scattering and drag) with a background of
! constant density and temperature (ions + electrons, at rest). Velocity-only
! collisions: no guiding-centre shift, no thermal force, no anomalous transport.
! Cross-field motion comes from the guiding-centre drifts only.
!
! Output:
!   impurity_tracing_exits.dat      one line per marker: status, time and position
!                                   (R, phi, Z) at the end, v_par, v_perp
!   impurity_tracing_residence.dat  time spent per 2D mesh triangle, summed over all
!                                   markers and toroidal slices [s]. For a source of
!                                   S markers/s, n(triangle) = S*time/(N*V_triangle).
! Charge-state changes (species_element /= ''): markers of one element (e.g. Be as in Harting et
! al. 2025, Sec. 5) change their charge state with constant ionisation/recombination rates per
! charge state (marker_species_mod): velocity kept, guiding-centre jump. Residence times are then
! written per charge state.
! Random numbers per marker (thermal_force_rng_mod): results do not depend on
! the number of threads. CGS units, temperatures in eV, masses in amu.

    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use thermal_force_rng_mod, only: rng_t, rng_seed, rng_normal3, rng_uniform
    use marker_species_mod, only: adas_rates, species_set_t, init_species_set, set_marker_charge, total_rate, &
                                  next_event_time, choose_transition, change_charge_state

    implicit none

    private

    real(dp), parameter :: amu_g = 1.66053906660e-24_dp
    real(dp), parameter :: me_g = 9.1093837015e-28_dp
    real(dp), parameter :: pi = 3.14159265358979323846264338327950288_dp

    integer, parameter :: max_z = 74
    integer, parameter :: st_time = 0, st_wall = 1, st_inner = 2, st_pusher = 3, st_outside = -1

    type :: impurity_tracing_input_t
        real(dp) :: test_mass_amu = 183.84_dp
        real(dp) :: test_charge_number = 5.0_dp
        real(dp) :: background_mass_amu = 2.01410177812_dp
        real(dp) :: background_charge_number = 1.0_dp
        real(dp) :: electron_density_cm3 = 1.0e13_dp   ! ion density = n_e / Z_b
        real(dp) :: temperature_ev = 20.0_dp           ! T_i = T_e = T_test at start
        logical :: boole_collisions = .true.
        logical :: boole_electron_collisions = .true.
        real(dp) :: dt_free_s = 1.0e-6_dp               ! step without collisions
        ! start: 'point' (start_R_cm, start_Z_cm, phi random) or 'region'
        ! (uniform in the volume of the region)
        character(len=32) :: start_mode = 'point'
        real(dp) :: start_R_cm = 0.0_dp
        real(dp) :: start_Z_cm = 0.0_dp
        real(dp) :: t_max_s = 1.0e-3_dp
        integer :: n_particles = 1000
        integer :: seed = 1
        character(len=256) :: output_prefix = 'impurity_tracing'
        ! charge-state changes: element ('' = off, fixed species test_mass_amu, test_charge_number)
        character(len=8) :: species_element = ''
        integer :: start_charge = 1
        real(dp) :: ionisation_rate_s(0:max_z) = 0.0_dp     ! per charge state [1/s]
        real(dp) :: recombination_rate_s(0:max_z) = 0.0_dp  ! per charge state [1/s]
        logical :: boole_gc_jump = .true.
        ! ADAS ADF11 files (e.g. scd96_be.dat, acd96_be.dat); if given, they replace the
        ! ionisation / recombination rates by n_e <sigma v>(n_e, T_e) of the background
        character(len=256) :: adas_scd_file = ''
        character(len=256) :: adas_acd_file = ''
    end type impurity_tracing_input_t

    public :: calc_impurity_tracing

contains

subroutine load_impurity_tracing_inp(inp)

    type(impurity_tracing_input_t), intent(out) :: inp
    real(dp) :: test_mass_amu, test_charge_number, background_mass_amu, background_charge_number, &
                electron_density_cm3, temperature_ev, dt_free_s, start_R_cm, start_Z_cm, t_max_s
    logical :: boole_collisions, boole_electron_collisions
    character(len=32) :: start_mode
    integer :: n_particles, seed, iu, ios
    character(len=256) :: output_prefix
    character(len=8) :: species_element
    integer :: start_charge
    real(dp) :: ionisation_rate_s(0:max_z), recombination_rate_s(0:max_z)
    logical :: boole_gc_jump
    character(len=256) :: adas_scd_file, adas_acd_file
    namelist /impurity_tracing_nml/ test_mass_amu, test_charge_number, background_mass_amu, &
        background_charge_number, electron_density_cm3, temperature_ev, boole_collisions, &
        boole_electron_collisions, dt_free_s, start_mode, start_R_cm, start_Z_cm, t_max_s, &
        n_particles, seed, output_prefix, species_element, start_charge, ionisation_rate_s, &
        recombination_rate_s, boole_gc_jump, adas_scd_file, adas_acd_file

    test_mass_amu = inp%test_mass_amu; test_charge_number = inp%test_charge_number
    background_mass_amu = inp%background_mass_amu; background_charge_number = inp%background_charge_number
    electron_density_cm3 = inp%electron_density_cm3; temperature_ev = inp%temperature_ev
    boole_collisions = inp%boole_collisions; boole_electron_collisions = inp%boole_electron_collisions
    dt_free_s = inp%dt_free_s; start_mode = inp%start_mode
    start_R_cm = inp%start_R_cm; start_Z_cm = inp%start_Z_cm; t_max_s = inp%t_max_s
    n_particles = inp%n_particles; seed = inp%seed; output_prefix = inp%output_prefix
    species_element = inp%species_element; start_charge = inp%start_charge
    ionisation_rate_s = inp%ionisation_rate_s; recombination_rate_s = inp%recombination_rate_s
    boole_gc_jump = inp%boole_gc_jump
    adas_scd_file = inp%adas_scd_file; adas_acd_file = inp%adas_acd_file

    open(newunit=iu, file='impurity_tracing.inp', status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'impurity tracing: cannot open impurity_tracing.inp'
    read(iu, nml=impurity_tracing_nml, iostat=ios)
    if (ios /= 0) error stop 'impurity tracing: error reading namelist impurity_tracing_nml in impurity_tracing.inp'
    close(iu)

    inp%test_mass_amu = test_mass_amu; inp%test_charge_number = test_charge_number
    inp%background_mass_amu = background_mass_amu; inp%background_charge_number = background_charge_number
    inp%electron_density_cm3 = electron_density_cm3; inp%temperature_ev = temperature_ev
    inp%boole_collisions = boole_collisions; inp%boole_electron_collisions = boole_electron_collisions
    inp%dt_free_s = dt_free_s; inp%start_mode = start_mode
    inp%start_R_cm = start_R_cm; inp%start_Z_cm = start_Z_cm; inp%t_max_s = t_max_s
    inp%n_particles = n_particles; inp%seed = seed; inp%output_prefix = output_prefix
    inp%species_element = species_element; inp%start_charge = start_charge
    inp%ionisation_rate_s = ionisation_rate_s; inp%recombination_rate_s = recombination_rate_s
    inp%boole_gc_jump = boole_gc_jump
    inp%adas_scd_file = adas_scd_file; inp%adas_acd_file = adas_acd_file

    if (trim(inp%start_mode) /= 'point' .and. trim(inp%start_mode) /= 'region') &
        error stop "impurity tracing: start_mode must be 'point' or 'region'"
    if (inp%t_max_s <= 0.0_dp .or. inp%n_particles < 1) error stop 'impurity tracing: need t_max_s > 0, n_particles >= 1'

end subroutine load_impurity_tracing_inp

subroutine calc_impurity_tracing()
    use tetra_physics_mod, only: cm_over_e, particle_charge, particle_mass

    use constants, only: clight, echarge, ev2erg
    use orbit_timestep_gorilla_mod, only: initialize_gorilla
    use tetra_physics_mod, only: particle_mass, particle_charge, cm_over_e
    use tetra_physics_poly_precomp_mod, only: make_precomp_poly
    use gorilla_settings_mod, only: coord_system
    use tetra_grid_mod, only: ntetr, tetra_in_edge_region, tetra_grid, verts_rphiz
    use tetra_grid_settings_mod, only: grid_size, n_field_periods, s_edge_region
    use collis_ions, only: collis_init
    use tetra_physics_mod, only: check_species_switch_allowed

    type(impurity_tracing_input_t) :: inp
    integer :: nphi, ntri, ip, itri, k, nb, iu
    integer, allocatable :: status(:)
    real(dp), allocatable :: final(:,:), residence(:,:), cum_weight(:), tri_rz(:,:,:)
    type(species_set_t) :: ss
    logical :: boole_species
    integer :: nstate, iz
    integer(int64) :: n_jumps
    real(dp) :: sum_jump
    real(dp) :: m_a, v0, mb(2), zb(2), db(2), tb(2), efcolf(2), velrat(2), enrat(2)
    integer :: n_status(-1:3)

    call load_impurity_tracing_inp(inp)
    if (coord_system /= 1) error stop 'impurity tracing: needs coord_system = 1 (R, phi, Z) in gorilla.inp'
    call initialize_gorilla()
    boole_species = (len_trim(inp%species_element) > 0)
    m_a = inp%test_mass_amu*amu_g
    particle_mass = m_a
    particle_charge = inp%test_charge_number*echarge
    cm_over_e = clight*particle_mass/particle_charge
    call make_precomp_poly()
    if (.not. allocated(tetra_in_edge_region)) &
        error stop 'impurity tracing: needs the extended mesh (tetra_grid.inp: grid_kind = 2, boole_extend_to_wall)'

    nphi = grid_size(2)
    if (mod(ntetr, 3*nphi) /= 0) error stop 'impurity tracing: ntetr is not 3 * triangles * slices'
    ntri = ntetr/(3*nphi)

    ! background: ions and (optionally) electrons, constant n and T
    nb = merge(2, 1, inp%boole_electron_collisions)
    mb = [inp%background_mass_amu*amu_g, me_g]
    zb = [inp%background_charge_number, -1.0_dp]
    db = [inp%electron_density_cm3/inp%background_charge_number, inp%electron_density_cm3]
    tb = inp%temperature_ev
    ! velocity normalisation v0 = sqrt(2 T/m_a)
    if (boole_species) then
        call check_species_switch_allowed()
        if (len_trim(inp%adas_scd_file) > 0 .or. len_trim(inp%adas_acd_file) > 0) then
            print '(a)', ' impurity tracing: rates from ADAS: '//trim(inp%adas_scd_file)//' '//trim(inp%adas_acd_file)
            call adas_rates(inp%species_element, inp%adas_scd_file, inp%adas_acd_file, inp%electron_density_cm3, &
                            inp%temperature_ev, inp%ionisation_rate_s, inp%recombination_rate_s)
        end if
        call init_species_set(ss, inp%species_element, mb(1:nb), zb(1:nb), db(1:nb), tb(1:nb), &
                              inp%temperature_ev, inp%ionisation_rate_s, inp%recombination_rate_s, inp%boole_gc_jump)
        m_a = ss%mass_g
        v0 = ss%v0
        call set_marker_charge(ss, inp%start_charge)
        nstate = ss%z_nuc + 1                 ! columns for charge states 0 ... z_nuc
        print '(a,a,a,i0,a,f8.4,a,i0,a,l1)', ' impurity tracing: element ', trim(ss%element), ', Z_nuc = ', ss%z_nuc, &
            ', m = ', ss%mass_amu, ' amu, start charge ', inp%start_charge, ', guiding-centre jump ', ss%boole_gc_jump
        do iz = 0, ss%z_nuc
            if (total_rate(ss, iz) > 0.0_dp) print '(a,i3,a,es11.4,a,es11.4)', '   charge ', iz, &
                ': ionisation rate [1/s] ', ss%rate_ion(iz), ', recombination rate [1/s] ', ss%rate_rec(iz)
        end do
    else
        call collis_init(m_a, inp%test_charge_number, mb(1:nb), zb(1:nb), db(1:nb), tb(1:nb), &
                         inp%temperature_ev, v0, efcolf(1:nb), velrat(1:nb), enrat(1:nb))
        nstate = 1
    end if

    ! 2D triangles of the region (first slice) for start_mode = 'region': weight = area * R_centroid
    allocate(cum_weight(ntri), tri_rz(2, 3, ntri))
    do itri = 1, ntri
        k = 3*(itri - 1) + 1          ! first tetrahedron of the prism of triangle itri in slice 1
        call prism_triangle(k, tri_rz(:,:,itri))
        cum_weight(itri) = 0.0_dp
        if (tetra_in_edge_region(k)) cum_weight(itri) = triangle_area(tri_rz(:,:,itri))*sum(tri_rz(1,:,itri))/3.0_dp
    end do
    do itri = 2, ntri
        cum_weight(itri) = cum_weight(itri) + cum_weight(itri - 1)
    end do
    if (cum_weight(ntri) <= 0.0_dp) error stop 'impurity tracing: the region is empty'

    print '(a,es12.4,a,i0,a,i0,a)', ' impurity tracing: s_edge_region = ', s_edge_region, ', ', &
        count(tetra_in_edge_region), ' of ', ntetr, ' tetrahedra'
    if (.not. boole_species) print '(a,f8.2,a,f5.1,a,es10.3,a,f7.2,a)', ' impurity tracing: test ion m = ', inp%test_mass_amu, ' amu, Z = ', &
        inp%test_charge_number, ', background n_e = ', inp%electron_density_cm3, ' cm^-3, T = ', inp%temperature_ev, ' eV'

    allocate(status(inp%n_particles), final(10, inp%n_particles), residence(ntri, nstate))
    residence = 0.0_dp
    n_jumps = 0
    sum_jump = 0.0_dp

    !$omp parallel do schedule(dynamic, 4) default(shared) private(ip) reduction(+:residence, n_jumps, sum_jump) &
    !$omp& COPYIN(cm_over_e,particle_charge,particle_mass)
    do ip = 1, inp%n_particles
        call trace_marker(ip, status(ip), final(:, ip), residence, n_jumps, sum_jump)
    end do
    !$omp end parallel do

    n_status = 0
    do ip = 1, inp%n_particles
        n_status(status(ip)) = n_status(status(ip)) + 1
    end do
    print '(a,5(a,i0))', ' impurity tracing: markers', ' wall ', n_status(st_wall), ', inner boundary ', &
        n_status(st_inner), ', time limit ', n_status(st_time), ', pusher removed ', n_status(st_pusher), &
        ', started outside ', n_status(st_outside)
    if (boole_species) then
        print '(a,15i8)', ' impurity tracing: final charge state counts (0 ... Z_nuc): ', &
            [(count(nint(final(9, :)) == iz), iz = 0, ss%z_nuc)]
        if (n_jumps > 0) print '(a,i0,a,es12.4,a)', ' impurity tracing: ', n_jumps, &
            ' charge-state transitions, mean guiding-centre jump ', sum_jump/real(n_jumps, dp), ' cm'
    end if

    open(newunit=iu, file=trim(inp%output_prefix)//'_exits.dat', status='replace', action='write')
    write(iu, '(a)') '# status (1 wall, 2 inner boundary, 0 time limit, 3 pusher removed, -1 started outside)'
    write(iu, '(a)') '# marker status t[s] R[cm] phi[rad] Z[cm] v_par[cm/s] v_perp[cm/s] R_start[cm] Z_start[cm]' // &
        ' charge n_transitions'
    write(iu, '(a,es24.16,a,es24.16,a,i0,a,es24.16)') '# s_edge_region = ', s_edge_region, ' t_max_s = ', &
        inp%t_max_s, ' n_particles = ', inp%n_particles, ' v0[cm/s] = ', v0
    do ip = 1, inp%n_particles
        write(iu, '(i10,i4,8es24.16,2i6)') ip, status(ip), final(1:8, ip), nint(final(9:10, ip))
    end do
    close(iu)

    open(newunit=iu, file=trim(inp%output_prefix)//'_residence.dat', status='replace', action='write')
    write(iu, '(a)') '# time spent per 2D triangle [s], summed over markers and toroidal slices; one column per'
    write(iu, '(a)') '# charge state 0 ... Z_nuc (charge-state mode) or one column (fixed species). ntri n_particles t_max_s ncol'
    write(iu, '(i10,i10,es24.16,i6)') ntri, inp%n_particles, inp%t_max_s, nstate
    do itri = 1, ntri
        write(iu, '(*(es24.16))') residence(itri, :)
    end do
    close(iu)

contains

    subroutine prism_triangle(ind_tetr, rz)
        ! (R, Z) of the 2D triangle of the prism that contains tetrahedron ind_tetr:
        ! the three distinct (R, Z) of its knots (axisymmetric extrusion, no shear)
        integer, intent(in) :: ind_tetr
        real(dp), intent(out) :: rz(2, 3)
        integer :: j, n, l
        real(dp) :: p(2)
        n = 0
        do j = 1, 4
            p = verts_rphiz([1, 3], tetra_grid(ind_tetr)%ind_knot(j))
            if (n > 0) then
                if (any([(all(abs(p - rz(:, l)) < 1.0e-9_dp), l = 1, n)])) cycle
            end if
            if (n == 3) error stop 'impurity tracing: prism tetrahedron with more than 3 distinct (R,Z) knots'
            n = n + 1
            rz(:, n) = p
        end do
        if (n /= 3) error stop 'impurity tracing: prism tetrahedron with less than 3 distinct (R,Z) knots'
    end subroutine prism_triangle

    pure function triangle_area(rz) result(a)
        real(dp), intent(in) :: rz(2, 3)
        real(dp) :: a
        a = 0.5_dp*abs((rz(1,2) - rz(1,1))*(rz(2,3) - rz(2,1)) - (rz(1,3) - rz(1,1))*(rz(2,2) - rz(2,1)))
    end function triangle_area

    subroutine trace_marker(ip, stat, fin, res, n_jumps, sum_jump)

        use find_tetra_mod, only: find_tetra
        use collis_ions, only: stost
        use tetra_physics_mod, only: tetra_physics
        use gorilla_settings_mod, only: optional_quantities_type
        use supporting_functions_mod, only: vperp_func
        use utils_orbit_timestep_mod, only: push_within_tetrahedron, initialize_constants_of_motion

        integer, intent(in) :: ip
        integer, intent(out) :: stat
        real(dp), intent(out) :: fin(10)
        real(dp), intent(inout) :: res(:,:)
        integer(int64), intent(inout) :: n_jumps
        real(dp), intent(inout) :: sum_jump

        type(rng_t) :: rng
        type(optional_quantities_type) :: optional_quantities
        real(dp) :: x(3), xi(3), vpar, vperp, v, t, h, dtau, zet(5), randnum(3), z_save(3), perpinv, &
                    t_remain, t_pass, x_start(3), u, w(3), t_event, t_stop, jump(3)
        integer :: ind_tetr, iface, ind_tetr_save, iper_phi, ierr, j, iz, icol, n_trans
        logical :: boole_t_finished, boole_pusher_failed

        call rng_seed(rng, int(inp%seed, int64), int(ip, int64))

        ! start position
        if (trim(inp%start_mode) == 'point') then
            x = [inp%start_R_cm, 0.0_dp, inp%start_Z_cm]
        else
            u = rng_uniform(rng)*cum_weight(ntri)
            j = min(ntri, max(1, count(cum_weight < u) + 1))
            w(1) = rng_uniform(rng); w(2) = rng_uniform(rng)
            if (w(1) + w(2) > 1.0_dp) w(1:2) = 1.0_dp - w(1:2)
            w(3) = 1.0_dp - w(1) - w(2)
            x = [sum(w*tri_rz(1,:,j)), 0.0_dp, sum(w*tri_rz(2,:,j))]
        end if
        x(2) = rng_uniform(rng)*2.0_dp*pi/n_field_periods
        x_start = x

        ! Maxwellian at the background temperature
        call rng_normal3(rng, xi)
        xi = xi*sqrt(inp%temperature_ev*ev2erg/m_a)
        vpar = xi(1)
        vperp = sqrt(xi(2)**2 + xi(3)**2)

        t = 0.0_dp
        stat = st_time
        n_trans = 0
        iz = 1
        icol = 1
        t_event = huge(1.0_dp)
        if (boole_species) then
            iz = inp%start_charge
            icol = iz + 1
            call set_marker_charge(ss, iz)     ! this thread's previous marker may have ended in another state
            t_event = next_event_time(ss, iz, t, rng)
        end if
        call find_tetra(x, vpar, vperp, ind_tetr, iface)
        if (ind_tetr == -1) then
            stat = st_outside
        else if (.not. tetra_in_edge_region(ind_tetr)) then
            stat = st_outside
        end if

        do while (stat == st_time .and. t < inp%t_max_s)
            t_stop = min(inp%t_max_s, t_event)
            ! collision, adaptive step dtau (length at v0) from stost, at most up to the next event
            if (inp%boole_collisions) then
                v = sqrt(vpar**2 + vperp**2)
                zet(1:3) = x
                zet(4) = v/v0
                zet(5) = vpar/v
                u = rng_uniform(rng); randnum(1) = sqrt(12.0_dp)*(u - 0.5_dp)
                randnum(2) = rng_uniform(rng)
                u = rng_uniform(rng); randnum(3) = sqrt(12.0_dp)*(u - 0.5_dp)
                if (boole_species) then
                    call stost(ss%efcolf(:, iz), ss%velrat(:, iz), ss%enrat(:, iz), zet, dtau, 1, ierr, &
                               (t_stop - t)*v0, randnum)
                else
                    call stost(efcolf(1:nb), velrat(1:nb), enrat(1:nb), zet, dtau, 1, ierr, &
                               (t_stop - t)*v0, randnum)
                end if
                vpar = zet(5)*zet(4)*v0
                vperp = sqrt(max(0.0_dp, 1.0_dp - zet(5)**2))*zet(4)*v0
                h = dtau/v0
            else
                h = min(inp%dt_free_s, t_stop - t)
            end if

            ! guiding-centre orbit over h, stop at the region boundary
            z_save = x - tetra_physics(ind_tetr)%x1
            call initialize_constants_of_motion(vperp, z_save, ind_tetr, perpinv)
            t_remain = h
            do
                ind_tetr_save = ind_tetr
                call push_within_tetrahedron(ind_tetr, iface, x, vpar, z_save, t_remain, t_pass, &
                                             boole_t_finished, iper_phi, optional_quantities, boole_pusher_failed)
                res(mod(ind_tetr_save - 1, 3*ntri)/3 + 1, icol) = res(mod(ind_tetr_save - 1, 3*ntri)/3 + 1, icol) + t_pass
                t = t + t_pass
                t_remain = t_remain - t_pass
                if (boole_pusher_failed) then
                    stat = st_pusher
                    exit
                end if
                vperp = vperp_func(z_save, perpinv, ind_tetr_save)
                if (ind_tetr == -1) then
                    stat = st_wall
                    exit
                end if
                if (.not. tetra_in_edge_region(ind_tetr)) then
                    stat = st_inner
                    exit
                end if
                if (boole_t_finished) exit
            end do

            ! charge-state transition when the event time is reached
            if (stat == st_time .and. boole_species .and. t >= t_event*(1.0_dp - 1.0e-12_dp)) then
                j = choose_transition(ss, iz, rng_uniform(rng))
                call change_charge_state(ss, iz, j, rng_uniform(rng), x, ind_tetr, iface, vpar, vperp, jump)
                n_jumps = n_jumps + 1
                sum_jump = sum_jump + norm2(jump)
                iz = j
                icol = iz + 1
                n_trans = n_trans + 1
                if (ind_tetr == -1) then
                    stat = st_wall
                else if (.not. tetra_in_edge_region(ind_tetr)) then
                    stat = st_inner
                end if
                t_event = next_event_time(ss, iz, t, rng)
            end if
        end do

        fin = [t, x(1), x(2), x(3), vpar, vperp, x_start(1), x_start(3), real(iz, dp), real(n_trans, dp)]


    end subroutine trace_marker

end subroutine calc_impurity_tracing

end module impurity_tracing_mod
