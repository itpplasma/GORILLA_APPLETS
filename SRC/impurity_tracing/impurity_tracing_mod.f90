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
! Neutrals (start_charge = 0): straight-line flight (neutral_flight_mod) with energy neutral_energy_ev
! in the direction neutral_direction (components along e_R, e_phi, e_Z at the start point) until
! ionisation with ionisation_rate_s(0) (rate per tetrahedron, track-length estimator), wall or t_max_s.
! At ionisation the ion keeps the neutral velocity; its guiding centre is X = x - (b x v)/Omega (the
! gyrophase is known here). Face test neutral_face_test = 'auto', 'cartesian' or 'hybrid' (see
! neutral_flight_mod; 'auto' uses cartesian if the face sag / inradius <= face_sag_tol).
! Random numbers per marker (thermal_force_rng_mod): results do not depend on
! the number of threads. CGS units, temperatures in eV, masses in amu.

    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use thermal_force_rng_mod, only: rng_t, rng_seed, rng_normal3, rng_uniform
    use marker_species_mod, only: adas_rates, adas_cell_rates, set_cell_rates, cell_rates, draw_optical_depth, &
                                  time_to_event, species_set_t, init_species_set, set_marker_charge, total_rate, &
                                  choose_transition, change_charge_state, field_direction
    use neutral_flight_mod, only: flight_result_t, choose_face_mode, fly_neutral, face_cartesian, &
                                  nf_ionised, nf_wall, nf_time, nf_error, cart_of_rphiz

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
        ! per-cell background for the rates (needs the ADAS files): one line 'n_e[cm^-3] T_e[eV]'
        ! per 2D triangle in mesh order (as in <prefix>_residence.dat), '#' lines are comments.
        ! Collisions still use electron_density_cm3, temperature_ev.
        character(len=256) :: background_cells_file = ''
        ! neutrals (start_charge = 0)
        real(dp) :: neutral_energy_ev = 25.0_dp
        real(dp) :: neutral_direction(3) = [-1.0_dp, 0.0_dp, 0.0_dp]   ! along e_R, e_phi, e_Z
        character(len=16) :: neutral_face_test = 'auto'
        real(dp) :: face_sag_tol = 0.1_dp
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
    character(len=256) :: adas_scd_file, adas_acd_file, background_cells_file
    real(dp) :: neutral_energy_ev, neutral_direction(3), face_sag_tol
    character(len=16) :: neutral_face_test
    namelist /impurity_tracing_nml/ test_mass_amu, test_charge_number, background_mass_amu, &
        background_charge_number, electron_density_cm3, temperature_ev, boole_collisions, &
        boole_electron_collisions, dt_free_s, start_mode, start_R_cm, start_Z_cm, t_max_s, &
        n_particles, seed, output_prefix, species_element, start_charge, ionisation_rate_s, &
        recombination_rate_s, boole_gc_jump, adas_scd_file, adas_acd_file, background_cells_file, &
        neutral_energy_ev, neutral_direction, neutral_face_test, face_sag_tol

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
    background_cells_file = inp%background_cells_file
    neutral_energy_ev = inp%neutral_energy_ev; neutral_direction = inp%neutral_direction
    neutral_face_test = inp%neutral_face_test; face_sag_tol = inp%face_sag_tol

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
    inp%background_cells_file = background_cells_file
    inp%neutral_energy_ev = neutral_energy_ev; inp%neutral_direction = neutral_direction
    inp%neutral_face_test = neutral_face_test; inp%face_sag_tol = face_sag_tol

    if (trim(inp%start_mode) /= 'point' .and. trim(inp%start_mode) /= 'region') &
        error stop "impurity tracing: start_mode must be 'point' or 'region'"
    if (inp%t_max_s <= 0.0_dp .or. inp%n_particles < 1) error stop 'impurity tracing: need t_max_s > 0, n_particles >= 1'
    if (inp%start_charge == 0) then
        if (len_trim(inp%species_element) == 0) error stop 'impurity tracing: neutrals (start_charge = 0) need species_element'
        if (inp%neutral_energy_ev <= 0.0_dp .or. norm2(inp%neutral_direction) <= 0.0_dp) &
            error stop 'impurity tracing: neutrals need neutral_energy_ev > 0 and a nonzero neutral_direction'
    end if

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
    ! neutrals
    logical :: boole_neutral
    integer :: face_mode, n_degenerate, n_neut(4), n_handover(0:4)
    integer, allocatable :: tri_of_tetr(:)
    real(dp) :: sag_ratio, sag_max, v_neutral
    real(dp), allocatable :: nu_tetr(:), ion_point(:,:)

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
        if (inp%start_charge >= 1) call set_marker_charge(ss, inp%start_charge)
        nstate = ss%z_nuc + 1                 ! columns for charge states 0 ... z_nuc
        print '(a,a,a,i0,a,f8.4,a,i0,a,l1)', ' impurity tracing: element ', trim(ss%element), ', Z_nuc = ', ss%z_nuc, &
            ', m = ', ss%mass_amu, ' amu, start charge ', inp%start_charge, ', guiding-centre jump ', ss%boole_gc_jump
        do iz = 0, ss%z_nuc
            if (total_rate(ss, iz) > 0.0_dp) print '(a,i3,a,es11.4,a,es11.4)', '   charge ', iz, &
                ': ionisation rate [1/s] ', ss%rate_ion(iz), ', recombination rate [1/s] ', ss%rate_rec(iz)
        end do
        if (len_trim(inp%background_cells_file) > 0) call load_cell_rates()
    else if (len_trim(inp%background_cells_file) > 0) then
        error stop 'impurity tracing: background_cells_file needs charge-state mode (species_element)'
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

    ! neutrals: face test, rate per tetrahedron (constant for now), triangle of each tetrahedron
    boole_neutral = boole_species .and. inp%start_charge == 0
    face_mode = 0
    v_neutral = 0.0_dp
    if (boole_neutral) then
        call choose_face_mode(inp%neutral_face_test, inp%face_sag_tol, face_mode, sag_ratio, sag_max, n_degenerate)
        print '(a,a,a,es10.3,a,es10.3,a,i0,a)', ' impurity tracing: neutral face test ', &
            trim(merge('cartesian', 'hybrid   ', face_mode == face_cartesian)), ' (max face sag ', sag_max, &
            ' cm, sag/inradius ', sag_ratio, ', ', n_degenerate, ' degenerate tetrahedra)'
        allocate(nu_tetr(ntetr), tri_of_tetr(ntetr))
        ! neutral ionisation rate per tetrahedron: rate of its cell (ADAS per triangle) or uniform
        do k = 1, ntetr
            nu_tetr(k) = ss%rate_ion(0)
            if (allocated(ss%cell_rate_ion)) nu_tetr(k) = ss%cell_rate_ion(0, mod(k - 1, 3*ntri)/3 + 1)
        end do
        do k = 1, ntetr
            tri_of_tetr(k) = mod(k - 1, 3*ntri)/3 + 1
        end do
        v_neutral = sqrt(2.0_dp*inp%neutral_energy_ev*ev2erg/ss%mass_g)
        print '(a,f8.2,a,es11.4,a,3f7.3,a,es11.4,a)', ' impurity tracing: neutrals ', inp%neutral_energy_ev, &
            ' eV (', v_neutral, ' cm/s), direction (e_R, e_phi, e_Z) ', inp%neutral_direction, &
            ', ionisation rate ', ss%rate_ion(0), ' 1/s'
    end if

    allocate(status(inp%n_particles), final(10, inp%n_particles), residence(ntri, nstate))
    allocate(ion_point(6, inp%n_particles))
    ion_point = 0.0_dp
    residence = 0.0_dp
    n_jumps = 0
    sum_jump = 0.0_dp
    n_neut = 0
    n_handover = 0

    !$omp parallel do schedule(dynamic, 4) default(shared) private(ip) &
    !$omp& reduction(+:residence, n_jumps, sum_jump, n_neut, n_handover) &
    !$omp& COPYIN(cm_over_e,particle_charge,particle_mass)
    do ip = 1, inp%n_particles
        call trace_marker(ip, status(ip), final(:, ip), residence, n_jumps, sum_jump, n_neut, n_handover, &
                          ion_point(:, ip))
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
        if (boole_neutral) then
            print '(a,4(a,i0))', ' impurity tracing: neutrals', ' ionised ', n_neut(nf_ionised), ', wall ', &
                n_neut(nf_wall), ', time limit ', n_neut(nf_time), ', flight error ', n_neut(nf_error)
            if (face_mode == face_cartesian) print '(a,5i9)', ' impurity tracing: cartesian handover ' // &
                '(direct, same, neighbour, find_tetra, failed): ', n_handover
        end if
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

    if (boole_neutral) then
        open(newunit=iu, file=trim(inp%output_prefix)//'_neutrals.dat', status='replace', action='write')
        write(iu, '(a)') '# end of the neutral flight per marker: flight status (1 ionised, 2 wall, 3 time limit, ' // &
            '4 error), t[s], particle position R[cm] phi[rad] Z[cm], tetrahedron'
        write(iu, '(a)') '# marker status t R phi Z ind_tetr'
        do ip = 1, inp%n_particles
            write(iu, '(i10,i4,4es24.16,i10)') ip, nint(ion_point(1, ip)), ion_point(2:5, ip), nint(ion_point(6, ip))
        end do
        close(iu)
    end if

    open(newunit=iu, file=trim(inp%output_prefix)//'_residence.dat', status='replace', action='write')
    write(iu, '(a)') '# time spent per 2D triangle [s], summed over markers and toroidal slices; one column per'
    write(iu, '(a)') '# charge state 0 ... Z_nuc (charge-state mode) or one column (fixed species). ntri n_particles t_max_s ncol'
    write(iu, '(i10,i10,es24.16,i6)') ntri, inp%n_particles, inp%t_max_s, nstate
    do itri = 1, ntri
        write(iu, '(*(es24.16))') residence(itri, :)
    end do
    close(iu)

contains

    subroutine load_cell_rates()
        ! n_e, T_e per 2D triangle from background_cells_file -> ADAS rates per triangle and
        ! charge state (uniform input rates where no ADAS file is given)
        real(dp), allocatable :: ne(:), te(:), ri(:,:), rr(:,:)
        character(len=512) :: line
        integer :: iu_c, ios, n
        if (len_trim(inp%adas_scd_file) == 0 .and. len_trim(inp%adas_acd_file) == 0) &
            error stop 'impurity tracing: background_cells_file needs adas_scd_file and/or adas_acd_file'
        allocate(ne(ntri), te(ntri))
        open(newunit=iu_c, file=trim(inp%background_cells_file), status='old', action='read', iostat=ios)
        if (ios /= 0) error stop 'impurity tracing: cannot open background_cells_file'
        n = 0
        do
            read(iu_c, '(a)', iostat=ios) line
            if (ios /= 0) exit
            line = adjustl(line)
            if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
            n = n + 1
            if (n > ntri) error stop 'impurity tracing: background_cells_file has more lines than mesh triangles'
            read(line, *, iostat=ios) ne(n), te(n)
            if (ios /= 0) error stop 'impurity tracing: background_cells_file: cannot read n_e T_e'
        end do
        close(iu_c)
        if (n /= ntri) then
            print '(a,i0,a,i0)', ' impurity tracing: background_cells_file has ', n, ' cells, mesh triangles: ', ntri
            error stop
        end if
        if (any(ne <= 0.0_dp) .or. any(te <= 0.0_dp)) error stop 'impurity tracing: background cells need n_e, T_e > 0'
        allocate(ri(0:ss%z_nuc, ntri), rr(0:ss%z_nuc, ntri))
        ri = spread(ss%rate_ion, 2, ntri)
        rr = spread(ss%rate_rec, 2, ntri)
        call adas_cell_rates(inp%species_element, inp%adas_scd_file, inp%adas_acd_file, ne, te, ri, rr)
        call set_cell_rates(ss, ri, rr)
        print '(a,a,a,i0,a)', ' impurity tracing: rates per triangle from ', trim(inp%background_cells_file), ' (', ntri, &
            ' cells; collisions keep the constant n_e, T)'
        print '(a,2es11.3,a,2es11.3)', '   n_e [cm^-3] min/max ', minval(ne), maxval(ne), ', T_e [eV] min/max ', &
            minval(te), maxval(te)
        do iz = 0, ss%z_nuc
            if (maxval(ri(iz, :) + rr(iz, :)) > 0.0_dp) print '(a,i3,a,2es11.3,a,2es11.3)', '   charge ', iz, &
                ': ionisation rate [1/s] min/max ', minval(ri(iz, :)), maxval(ri(iz, :)), &
                ', recombination min/max ', minval(rr(iz, :)), maxval(rr(iz, :))
        end do
    end subroutine load_cell_rates

    pure function cell_of_tetra(ind_tetr) result(icell)
        ! cell (2D triangle) of tetrahedron ind_tetr; 0 = uniform rates
        integer, intent(in) :: ind_tetr
        integer :: icell
        icell = 0
        if (allocated(ss%cell_rate_ion)) icell = mod(ind_tetr - 1, 3*ntri)/3 + 1
    end function cell_of_tetra

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

    subroutine trace_marker(ip, stat, fin, res, n_jumps, sum_jump, n_neut, n_handover, ion_pt)

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
        integer, intent(inout) :: n_neut(4), n_handover(0:4)
        real(dp), intent(out) :: ion_pt(6)

        type(rng_t) :: rng
        type(optional_quantities_type) :: optional_quantities
        real(dp) :: x(3), xi(3), vpar, vperp, v, t, h, dtau, zet(5), randnum(3), z_save(3), perpinv, &
                    t_remain, t_pass, x_start(3), u, w(3), t_stop, jump(3), tau, tau_acc, nu, nu_ion, nu_rec, &
                    t_push
        integer :: ind_tetr, iface, ind_tetr_save, iper_phi, ierr, j, iz, icol, n_trans
        type(flight_result_t) :: fl
        real(dp) :: dir(3), vcart(3), cphi, sphi
        logical :: boole_t_finished, boole_pusher_failed, boole_event

        call rng_seed(rng, int(inp%seed, int64), int(ip, int64))
        ion_pt = 0.0_dp

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
        ! transition when the accumulated optical depth tau_acc (sum of nu dt) reaches tau
        tau = huge(1.0_dp)
        tau_acc = 0.0_dp
        if (boole_species) then
            iz = inp%start_charge
            icol = iz + 1
            ! this thread's previous marker may have ended in another state
            if (iz >= 1) then
                call set_marker_charge(ss, iz)
                if (total_rate(ss, iz) > 0.0_dp .or. allocated(ss%cell_rate_ion)) tau = draw_optical_depth(rng)
            end if
        end if
        if (iz == 0) then
            ! neutral: velocity along neutral_direction (local e_R, e_phi, e_Z at the start point)
            dir = inp%neutral_direction/norm2(inp%neutral_direction)
            cphi = cos(x(2)); sphi = sin(x(2))
            vcart = v_neutral*[dir(1)*cphi - dir(2)*sphi, dir(1)*sphi + dir(2)*cphi, dir(3)]
            vpar = v_neutral; vperp = 0.0_dp
        end if
        call find_tetra(x, vpar, vperp, ind_tetr, iface)
        if (ind_tetr == -1) then
            stat = st_outside
        else if (.not. tetra_in_edge_region(ind_tetr)) then
            stat = st_outside
        end if
        if (stat == st_time .and. iz == 0) then
            call fly_neutral(face_mode, cart_of_rphiz(x), vcart, t, ind_tetr, inp%t_max_s, &
                             draw_optical_depth(rng), nu_tetr, fl, tri_of_tetr, res(:, 1))
            n_neut(fl%status) = n_neut(fl%status) + 1
            ion_pt = [real(fl%status, dp), fl%t, fl%x(1), fl%x(2), fl%x(3), real(fl%ind_tetr, dp)]
            if (face_mode == face_cartesian) n_handover(fl%handover) = n_handover(fl%handover) + 1
            t = fl%t
            x = fl%x
            ind_tetr = fl%ind_tetr
            iface = 0
            select case (fl%status)
            case (nf_wall)
                stat = st_wall
            case (nf_time)
                stat = st_time
            case (nf_error)
                stat = st_pusher
            case (nf_ionised)
                if (ind_tetr < 1) then
                    stat = st_pusher
                else
                    call ionise_neutral(x, ind_tetr, iface, vcart, vpar, vperp, jump)
                    n_jumps = n_jumps + 1
                    sum_jump = sum_jump + norm2(jump)
                    n_trans = n_trans + 1
                    iz = 1
                    icol = 2
                    call set_marker_charge(ss, iz)
                    tau_acc = 0.0_dp
                    tau = huge(1.0_dp)
                    if (total_rate(ss, iz) > 0.0_dp .or. allocated(ss%cell_rate_ion)) tau = draw_optical_depth(rng)
                    if (ind_tetr == -1) then
                        stat = st_wall
                    else if (.not. tetra_in_edge_region(ind_tetr)) then
                        stat = st_inner
                    end if
                end if
            end select
        end if

        do while (stat == st_time .and. t < inp%t_max_s)
            ! step end: t_max or the transition, if the rate of the current cell stayed constant
            nu = 0.0_dp
            if (boole_species) then
                call cell_rates(ss, iz, cell_of_tetra(ind_tetr), nu_ion, nu_rec)
                nu = nu_ion + nu_rec
            end if
            t_stop = inp%t_max_s
            if (nu > 0.0_dp) t_stop = min(t_stop, t + time_to_event(nu, tau_acc, tau))
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
            boole_event = .false.
            do
                ind_tetr_save = ind_tetr
                ! rate in this tetrahedron's cell; the push stops where the optical depth is reached
                nu = 0.0_dp
                if (boole_species) then
                    call cell_rates(ss, iz, cell_of_tetra(ind_tetr), nu_ion, nu_rec)
                    nu = nu_ion + nu_rec
                end if
                t_push = min(t_remain, time_to_event(nu, tau_acc, tau))
                call push_within_tetrahedron(ind_tetr, iface, x, vpar, z_save, t_push, t_pass, &
                                             boole_t_finished, iper_phi, optional_quantities, boole_pusher_failed)
                res(mod(ind_tetr_save - 1, 3*ntri)/3 + 1, icol) = res(mod(ind_tetr_save - 1, 3*ntri)/3 + 1, icol) + t_pass
                t = t + t_pass
                t_remain = t_remain - t_pass
                tau_acc = tau_acc + nu*t_pass
                ! transition reached inside this tetrahedron (orbit step ends early if t_push < t_remain)
                boole_event = boole_t_finished .and. nu > 0.0_dp .and. tau_acc >= tau*(1.0_dp - 1.0e-12_dp)
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

            ! charge-state transition when the optical depth is reached (branching with the cell rates)
            if (stat == st_time .and. boole_event) then
                j = choose_transition(ss, iz, rng_uniform(rng), cell_of_tetra(ind_tetr_save))
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
                tau_acc = 0.0_dp
                tau = huge(1.0_dp)
                if (total_rate(ss, iz) > 0.0_dp .or. allocated(ss%cell_rate_ion)) tau = draw_optical_depth(rng)
            end if
        end do

        fin = [t, x(1), x(2), x(3), vpar, vperp, x_start(1), x_start(3), real(iz, dp), real(n_trans, dp)]

    end subroutine trace_marker

    subroutine ionise_neutral(x, ind_tetr, iface, vcart, vpar, vperp, jump)

        ! Neutral -> charge 1 at the particle position x (R, phi, Z) with Cartesian velocity vcart:
        ! v_par = v.b, v_perp = |v - v_par b|, guiding centre X = x - (b x v)/Omega, Omega = e B/(m c).
        ! ind_tetr = -1 on return if the guiding centre lies outside the mesh.

        use anomalous_transport_displacement_mod, only: displace_by_straight_line

        real(dp), intent(inout) :: x(3)
        integer, intent(inout) :: ind_tetr, iface
        real(dp), intent(in) :: vcart(3)
        real(dp), intent(out) :: vpar, vperp, jump(3)
        real(dp) :: b_hat(3), bmod, vloc(3), c, s, bxv(3)

        ! velocity in the local physical basis (e_R, e_phi, e_Z) at x
        c = cos(x(2)); s = sin(x(2))
        vloc = [c*vcart(1) + s*vcart(2), -s*vcart(1) + c*vcart(2), vcart(3)]
        call field_direction(x, ind_tetr, b_hat, bmod)
        vpar = dot_product(vloc, b_hat)
        vperp = norm2(vloc - vpar*b_hat)
        bxv = [b_hat(2)*vloc(3) - b_hat(3)*vloc(2), b_hat(3)*vloc(1) - b_hat(1)*vloc(3), &
               b_hat(1)*vloc(2) - b_hat(2)*vloc(1)]
        jump = -bxv*ss%mass_g*clight/(echarge*bmod)
        call displace_by_straight_line(x, ind_tetr, iface, [jump(1), jump(2)/x(1), jump(3)], vpar, vperp)

    end subroutine ionise_neutral

end subroutine calc_impurity_tracing

end module impurity_tracing_mod
