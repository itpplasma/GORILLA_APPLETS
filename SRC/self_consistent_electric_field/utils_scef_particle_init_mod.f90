module utils_scef_particle_init_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

contains

subroutine calc_desired_density_profile

    !Fills desired_density%value on g%s_vertices and desired_density%grad on
    !the flux shells connecting them, from target_density (scef_target_profiles_mod).

    use gorilla_applets_types_mod, only: g, desired_density
    use tetra_grid_settings_mod, only: grid_size
    use scef_target_profiles_mod, only: target_density, rad_from_s

    integer :: ns

    if (.not.allocated(g%s_vertices)) then
        print*, 'Error, calc_desired_density_profile called before g%s_vertices was populated.'
        stop
    endif

    if (.not.allocated(desired_density%value)) allocate(desired_density%value(grid_size(1)+1))
    if (.not.allocated(desired_density%grad))  allocate(desired_density%grad(grid_size(1)))

    do ns = 1, grid_size(1)+1
        desired_density%value(ns) = target_density(rad_from_s(g%s_vertices(ns)))
    enddo

    do ns = 1, grid_size(1)
        desired_density%grad(ns) = (desired_density%value(ns+1) - desired_density%value(ns)) &
                                   /(g%s_vertices(ns+1) - g%s_vertices(ns))
    enddo

end subroutine calc_desired_density_profile

subroutine calc_starting_conditions

    use gorilla_applets_types_mod, only: in, start, particle_source
    use tetra_grid_settings_mod, only: sfc_s_min

    real(dp), dimension(:,:,:), allocatable                :: rand_matrix

    allocate(rand_matrix(5,in%num_particles,in%n_species))
    call RANDOM_NUMBER(rand_matrix)

    call allocate_start_type
    call allocate_weights
    call set_particle_type_specifications
    call init_marker_distributions(in%energy_eV)
    call set_starting_positions(rand_matrix)!,s0=sfc_s_min*1.1_dp)
    call set_rest_of_individual_particle_specifications(rand_matrix)
    if (allocated(particle_source%value)) call apply_particle_source_weights

end subroutine calc_starting_conditions

subroutine init_marker_distributions(temperature_eV)

    !Initialise the marker position/pitch/energy distributions so the central weight routine (and
    !the shared sampling helpers) can evaluate pdf_position/pdf_lambda/pdf_energy. Position: uniform
    !in flux coordinates on [sfc_s_min, 1] x [0, 2*pi] x [0, 2*pi/n_field_periods]; pitch: uniform
    !on [-1, 1]; energy: Maxwell sqrt(E)*exp(-E/T) on [0, epsilon_max*T] (eV).

    use gorilla_applets_types_mod, only: in, start
    use tetra_grid_settings_mod, only: sfc_s_min, n_field_periods
    use constants, only: pi
    use marker_distribution_mod, only: pdf_flat, pdf_maxwell, pdf_maxwell_integral, maxwell_temperature_eV, &
                                       init_distribution_3d, init_distribution_1d, cleanup_distribution_1d

    real(dp), intent(in) :: temperature_eV

    if (.not. start%dist_position%initialized) then
        call init_distribution_3d(start%dist_position, pdf_flat, &
                                  [sfc_s_min, 0.0_dp, 0.0_dp], &
                                  [1.0_dp, 2.0_dp*pi, 2.0_dp*pi/n_field_periods])
    endif

    if (.not. start%dist_lambda%initialized) then
        call init_distribution_1d(start%dist_lambda, pdf_flat, -1.0_dp, 1.0_dp)
    endif

    if (in%boole_monoenergetic) then
        if (.not. start%dist_energy%initialized) then
            call init_distribution_1d(start%dist_energy, pdf_flat, in%energy_eV, in%energy_eV)
        endif
    else
        !Re-initialise if temperature changed (e.g. between the main run and the diffusion calibration)
        if (start%dist_energy%initialized .and. maxwell_temperature_eV /= temperature_eV) then
            call cleanup_distribution_1d(start%dist_energy)
        endif
        if (.not. start%dist_energy%initialized) then
            maxwell_temperature_eV = temperature_eV
            call init_distribution_1d(start%dist_energy, pdf_maxwell, 0.0_dp, start%epsilon_max*temperature_eV, &
                                      normalization = pdf_maxwell_integral(start%epsilon_max*temperature_eV))
        endif
    endif

end subroutine init_marker_distributions

subroutine initialize_particle_source

    !Allocate particle_source and initialise it. With source tuning (n_source_updates > 0) the
    !desired density profile is the initial guess, so this must be called after
    !calc_desired_density_profile. Without tuning the source is fixed to the constant in%density,
    !which gives every marker the uniform physical density in%density.

    use gorilla_applets_types_mod, only: in, g, desired_density, particle_source
    use tetra_grid_settings_mod, only: grid_size

    integer :: k

    if (in%n_source_updates.gt.0 .and. .not.allocated(desired_density%value)) then
        print*, 'Error: initialize_particle_source called before calc_desired_density_profile.'
        stop
    endif

    if (.not.allocated(particle_source%value)) allocate(particle_source%value(grid_size(1)+1))
    if (.not.allocated(particle_source%grad))  allocate(particle_source%grad(grid_size(1)))

    if (in%n_source_updates.gt.0) then
        particle_source%value(:) = desired_density%value(:)
    else
        particle_source%value(:) = in%density
    endif

    do k = 1, grid_size(1)
        particle_source%grad(k) = (particle_source%value(k+1) - particle_source%value(k)) &
                                 /(g%s_vertices(k+1) - g%s_vertices(k))
    enddo

end subroutine initialize_particle_source

subroutine apply_particle_source_weights

    !Set the physical density represented by each marker to the particle-source profile evaluated at
    !the marker's s value. The central weight routine picks this up as `start%density(n, species)` and
    !multiplies by V_box * sqrt(g) (plus the velocity-space factors) to produce the marker weight.

    use gorilla_applets_types_mod, only: in, start, g, particle_source
    use tetra_grid_settings_mod, only: grid_size
    use binsrc_mod, only: binsrc

    integer :: n, species, k
    real(dp) :: s_marker

    do species = 1, in%n_species
        do n = 1, in%num_particles
            s_marker = start%x(1, n, species)
            call binsrc(g%s_vertices, 1, grid_size(1)+1, s_marker, k)
            k = k - 1
            k = max(1, min(k, grid_size(1)))
            start%density(n, species) = particle_source%value(k) + particle_source%grad(k)*(s_marker - g%s_vertices(k))
        enddo
    enddo

end subroutine apply_particle_source_weights

subroutine write_particle_source(source_step)

    !Write particle_source at all s-vertices to particle_source_<source_step>.dat.
    !Two columns: s_vertex, source_value.

    use gorilla_applets_types_mod, only: g, particle_source
    use tetra_grid_settings_mod, only: grid_size

    integer, intent(in) :: source_step
    integer :: k, file_id
    character(len=100) :: filename
    character(len=32) :: source_step_str

    write(source_step_str, '(I0)') source_step
    filename = 'particle_source_' // trim(source_step_str) // '.dat'
    open(newunit = file_id, file = filename)
    do k = 1, grid_size(1)+1
        write(file_id,*) g%s_vertices(k), particle_source%value(k)
    enddo
    close(file_id)

end subroutine write_particle_source


subroutine read_restart_particle_source

    !Start from the particle source of an earlier run (restart_source_file, as written by
    !write_particle_source: s_vertex, value on grid_size(1)+1 lines, same s-grid) instead of the target.

    use gorilla_applets_types_mod, only: in, g, particle_source
    use tetra_grid_settings_mod, only: grid_size

    integer :: k, file_id, ios
    real(dp) :: s_file

    open(newunit = file_id, file = trim(in%restart_source_file), status = 'old', action = 'read', iostat = ios)
    if (ios.ne.0) then
        print*, 'Error: cannot open restart_source_file ', trim(in%restart_source_file)
        stop
    endif
    do k = 1, grid_size(1)+1
        read(file_id, *, iostat = ios) s_file, particle_source%value(k)
        if ((ios.ne.0) .or. (abs(s_file - g%s_vertices(k)).gt.1.0d-6)) then
            print*, 'Error: restart_source_file ', trim(in%restart_source_file), ' does not match the s-grid of this run'
            stop
        endif
    enddo
    close(file_id)

    do k = 1, grid_size(1)
        particle_source%grad(k) = (particle_source%value(k+1) - particle_source%value(k)) &
                                 /(g%s_vertices(k+1) - g%s_vertices(k))
    enddo
    print*, 'restart: particle source read from ', trim(in%restart_source_file)

end subroutine read_restart_particle_source
subroutine update_particle_source

    !Richardson step S := S - alpha * (n - n_target) using the last inner-iteration
    !density in one_d%densities. Values live on s-vertices; the shell-averaged
    !one_d%densities is projected to vertices by midpoint averaging (nearest-shell
    !at the two boundaries). The outer 10% of the domain (or the outer buffer) keeps its source, because
    !particles only leave there and none enter, so n stays below n_target whatever
    !the source (source_buffer_mode = 1: the buffer is updated too, with |S| capped at
    !source_buffer_bound_factor times max |S| of the tuned region). Logs the |n - n_target| / |n_target|
    !residual over the tuned shells (never the buffer)
    !(the electron/ion mismatch is logged after every electric potential update).

    use gorilla_applets_types_mod, only: in, g, one_d, desired_density, particle_source
    use tetra_grid_settings_mod, only: grid_size, sfc_s_min
    use constants, only: pi

    integer :: k, ns_g
    real(dp) :: alpha, target_mismatch, denom, s_cutoff_high, s_bound
    real(dp), dimension(:), allocatable :: n_shell, n_vert, target_shell
    logical,  dimension(:), allocatable :: tuned_shell

    ns_g = grid_size(1)
    alpha = in%source_relaxation_factor
    s_cutoff_high = 1.0_dp - 0.1_dp*(1.0_dp - sfc_s_min)   !outer 10% of the domain is not tuned
    if (in%edge_buffer_width.gt.0.0_dp) s_cutoff_high = (1.0_dp - in%edge_buffer_width)**2   !nor the outer buffer

    allocate(n_shell(ns_g), n_vert(ns_g+1), target_shell(ns_g), tuned_shell(ns_g))

    n_shell(:) = 0.5_dp * (one_d%densities(:,1) + one_d%densities(:,2))

    n_vert(1)        = n_shell(1)
    do k = 2, ns_g
        n_vert(k) = 0.5_dp * (n_shell(k-1) + n_shell(k))
    enddo
    n_vert(ns_g+1) = n_shell(ns_g)

    if (in%source_buffer_mode.eq.1) then
        !Buffer updated like the rest, but |S| there is capped at source_buffer_bound_factor times the
        !largest |S| of the tuned region (n stays below the target in the buffer, so the plain update
        !would let the source grow without bound)
        particle_source%value = particle_source%value - alpha*(n_vert - desired_density%value)
        s_bound = in%source_buffer_bound_factor &
                  *maxval(abs(particle_source%value), mask=(g%s_vertices < s_cutoff_high))
        where (g%s_vertices >= s_cutoff_high)
            particle_source%value = sign(min(abs(particle_source%value), s_bound), particle_source%value)
        end where
    else
        where (g%s_vertices < s_cutoff_high)
            particle_source%value = particle_source%value - alpha*(n_vert - desired_density%value)
        end where
    endif


    do k = 1, ns_g
        particle_source%grad(k) = (particle_source%value(k+1) - particle_source%value(k)) &
                                 /(g%s_vertices(k+1) - g%s_vertices(k))
    enddo

    target_shell(:) = 0.5_dp*(desired_density%value(1:ns_g) + desired_density%value(2:ns_g+1))
    tuned_shell(:) = g%s_vertices(2:ns_g+1) < s_cutoff_high   !shells whose both vertices are tuned
    denom = sum(target_shell * g%s_shell_volumes, mask=tuned_shell)
    if (denom > 0.0_dp) then
        target_mismatch = sum(abs(n_shell - target_shell) * g%s_shell_volumes, mask=tuned_shell) / denom
    else
        target_mismatch = 0.0_dp
    endif
    one_d%source_residual = target_mismatch
    !its expected value from noise alone, with sigma(n_shell) = sqrt(sigma_i^2 + sigma_e^2)/2
    one_d%source_residual_noise = 0.0_dp
    if (allocated(one_d%sigma) .and. (denom > 0.0_dp)) one_d%source_residual_noise = sqrt(2.0_dp/pi) &
        *sum(0.5_dp*sqrt(one_d%sigma(:,1)**2 + one_d%sigma(:,2)**2) * g%s_shell_volumes, mask=tuned_shell)/denom

    print*, 'source update: relative |n - n_target| / |n_target| over the tuned shells = ', target_mismatch, &
            ', noise level = ', one_d%source_residual_noise

end subroutine update_particle_source

subroutine allocate_start_type(n_particles_in)

    use gorilla_applets_types_mod, only: start, in

    integer, intent(in), optional :: n_particles_in
    integer :: n_particles

    if (present(n_particles_in)) then
        n_particles = n_particles_in
    else
        n_particles = in%num_particles
    endif

    !before allocating, deallocate if necessary
    call deallocate_start_type

    allocate(start%x(3,n_particles,in%n_species))
    allocate(start%pitch(n_particles,in%n_species))
    allocate(start%energy(n_particles,in%n_species))
    allocate(start%jperp(n_particles,in%n_species))
    allocate(start%density(n_particles,in%n_species))
    allocate(start%lost(n_particles,in%n_species))
    allocate(start%particle_charge(in%n_species))
    allocate(start%particle_mass(in%n_species))
    allocate(start%cm_over_e(in%n_species))
    allocate(start%t(in%n_species))
    allocate(start%v0(in%n_species))

    start%density = in%density   !overwritten per marker in apply_particle_source_weights

end subroutine allocate_start_type

subroutine allocate_weights(n_particles_in)

    use gorilla_applets_types_mod, only: in, weights

    integer, intent(in), optional :: n_particles_in
    integer :: n_particles

    if (present(n_particles_in)) then
        n_particles = n_particles_in
    else
        n_particles = in%num_particles
    endif

    call deallocate_weights
    allocate(weights%w(n_particles,in%n_species))

end subroutine allocate_weights

subroutine deallocate_start_type

    use gorilla_applets_types_mod, only: start

    if (allocated(start%x))               deallocate(start%x)
    if (allocated(start%pitch))           deallocate(start%pitch)
    if (allocated(start%energy))          deallocate(start%energy)
    if (allocated(start%jperp))           deallocate(start%jperp)
    if (allocated(start%density))         deallocate(start%density)
    if (allocated(start%lost))            deallocate(start%lost)
    if (allocated(start%particle_charge)) deallocate(start%particle_charge)
    if (allocated(start%particle_mass))   deallocate(start%particle_mass)
    if (allocated(start%cm_over_e))       deallocate(start%cm_over_e)
    if (allocated(start%t))               deallocate(start%t)
    if (allocated(start%v0))              deallocate(start%v0)

end subroutine deallocate_start_type

subroutine deallocate_weights

    use gorilla_applets_types_mod, only: weights

    if (allocated(weights%w)) deallocate(weights%w)

end subroutine deallocate_weights

subroutine set_starting_positions(rand_matrix,species_in,s0)

    use gorilla_applets_types_mod, only: in, start, g
    use tetra_physics_mod, only: coord_system
    use tetra_grid_settings_mod, only: grid_kind
    use constants, only: pi
    use tetra_grid_settings_mod, only: sfc_s_min, n_field_periods

    real(dp), dimension(:,:,:), intent(in) :: rand_matrix
    integer, dimension(:), intent(in), optional :: species_in
    integer, dimension(:), allocatable :: species
    integer :: i
    real(dp), intent(in), optional :: s0

    if (present(species_in)) then
        allocate(species(size(species_in)))
        species = species_in
    else
        allocate(species(in%n_species))
        species = [(i,i=1,in%n_species)]
    endif

    start%x(1,:,species) = sfc_s_min + rand_matrix(1,:,:)*(1-sfc_s_min)
    if (present(s0).and.(.not.in%boole_static_ne))  start%x(1,:,species) = s0
    start%x(2,:,species) = 2*pi*rand_matrix(2,:,:) !theta
    start%x(3,:,species) = 2*pi/n_field_periods*rand_matrix(3,:,:) !phi

    !unless a single species is initiated, make elctrons and ions start at identical positions in real space
    if (size(species).gt.1) then
        start%x(:,:,in%n_species) = start%x(:,:,1)
    endif

end subroutine set_starting_positions

subroutine set_rest_of_individual_particle_specifications(rand_matrix,boole_diffusion_coefficient_in,species_in,n_particles_in)

    use gorilla_applets_types_mod, only: in, start, s
    use marker_distribution_mod, only: sample_array_1d, maxwell_temperature_eV
    use scef_target_profiles_mod, only: target_temperature, rad_from_s

    real(dp), dimension(:,:,:), intent(in) :: rand_matrix
    integer, dimension(:), intent(in), optional :: species_in
    logical, intent(in), optional :: boole_diffusion_coefficient_in
    logical :: boole_diffusion_coefficient
    integer, dimension(:), allocatable :: species
    integer :: i, isp, n
    integer, intent(in), optional :: n_particles_in
    integer :: n_particles

    if (present(species_in)) then
        allocate(species(size(species_in)))
        species = species_in
    else
        allocate(species(in%n_species))
        species = [(i,i=1,in%n_species)]
    endif

    if (present(n_particles_in)) then
        n_particles = n_particles_in
    else
        n_particles = in%num_particles
    endif

    boole_diffusion_coefficient = .false.   !not initialised in the declaration: that would imply SAVE
    if (present(boole_diffusion_coefficient_in)) boole_diffusion_coefficient = boole_diffusion_coefficient_in

    start%pitch(:,species) = 2*rand_matrix(4,:,:)-1 !pitch parameter
    start%energy(:,species) = in%energy_eV
    if ((.not. in%boole_monoenergetic).or.boole_diffusion_coefficient) then
        !Energies drawn through start%dist_energy (Maxwell sqrt(E)*exp(-E/T) in eV); the diffusion
        !calibration must re-initialise it at s%temperature before this call.
        do i = 1, size(species)
            isp = species(i)
            call sample_array_1d(start%dist_energy, n_particles, start%energy(:,isp))
        enddo
    endif

    if (in%boole_antithetic_variate) then
        start%x(:,1:n_particles:2,species) = start%x(:,2:n_particles:2,species)
        start%pitch(1:n_particles:2,species) = -start%pitch(2:n_particles:2,species)
        start%energy(1:n_particles:2,species) = start%energy(2:n_particles:2,species)
    endif

    !Markers carry the local target temperature: energies drawn at maxwell_temperature_eV are rescaled
    !by T(s)/maxwell_temperature_eV (calc_particle_weights_and_jperp rescales the marker pdf alike)
    if (.not. in%boole_monoenergetic) then
        do i = 1, size(species)
            isp = species(i)
            do n = 1, n_particles
                start%energy(n,isp) = start%energy(n,isp) &
                                      *target_temperature(rad_from_s(start%x(1,n,isp)), isp)/maxwell_temperature_eV
            enddo
        enddo
    endif

end subroutine set_rest_of_individual_particle_specifications

subroutine set_particle_type_specifications

    use gorilla_applets_types_mod, only: in, start
    use constants, only: echarge,ame,clight, ev2erg, amp
    use gorilla_settings_mod, only: ispecies

    real(dp) :: charge, mass, cm_over_e

    select case(ispecies)
        case(1) !electron
            charge = -echarge
            mass = ame
            cm_over_e = -clight*ame/echarge
        case(2) !deuterium ion
            charge = echarge
            mass = 2.d0*amp
            cm_over_e=2.d0*clight*amp/echarge
        case(3) !alpha particle
            charge = 2.d0*echarge
            mass = 4.d0*amp
            cm_over_e=2.d0*clight*amp/echarge
        case(4) !ionised tungsten
            charge = 74.d0*echarge
            mass = 184.d0*amp
            cm_over_e= 184.d0*clight*amp/(74.d0*echarge)
        case(5) !hydrogen ion
            charge = echarge
            mass = amp
            cm_over_e = clight*amp/echarge
    end select

    start%particle_charge = (/charge, -echarge/)
    start%particle_mass = (/mass, ame/)
    start%cm_over_e = (/cm_over_e, -clight*ame/echarge/)
    start%t = (/in%time_step, in%time_step/) !/42.0_dp
    if (in%ion_tracing_time_limit.gt.0.0_dp) start%t(1) = in%ion_tracing_time_limit
    if (in%boole_static_ne) start%t(2) = 0.0_dp

    start%v0 = sqrt(2.0_dp*in%energy_eV*ev2erg/start%particle_mass)
    start%epsilon_max = 16.0_dp

end subroutine set_particle_type_specifications

subroutine generate_distribution_x4_exp_neg_x(b, output_array)

     use gorilla_applets_types_mod, only: s, in

    real(dp), dimension(:,:), intent(out) :: output_array
    real(dp), intent(in) :: b

    real(dp), dimension(:), allocatable :: uniform_random
    real(dp) :: x, y, max_value
    integer :: j, accept_count, i, num_cols, row, num_rows

    num_cols = size(output_array, 1)
    num_rows = size(output_array, 2)
    allocate(uniform_random(num_cols))

    max_value = 4.0_dp**4.0_dp * (s%temperature/in%energy_eV)**4.0_dp*exp(-4.0_dp)

    do row = 1, num_rows
        accept_count = 0

        do while (accept_count < num_cols)
            call random_number(uniform_random)

            do j = 1, num_cols
                if (accept_count >= num_cols) exit

                x = b * uniform_random(j)
                y = x**4.0_dp * exp(-x/(s%temperature/in%energy_eV))

                call random_number(uniform_random(j))
                if (uniform_random(j) * max_value <= y) then
                    accept_count = accept_count + 1
                    output_array(accept_count, row) = x
                endif
            enddo
        enddo
    enddo

end subroutine generate_distribution_x4_exp_neg_x

subroutine generate_distribution_sqrt_x_exp_neg_x(xmax, output_array)

    real(dp), dimension(:,:), intent(out) :: output_array
    real(dp), intent(in) :: xmax

    real(dp), dimension(:), allocatable :: uniform_random
    real(dp) :: x, y, max_value
    integer :: j, accept_count, i, num_cols, row, num_rows

    num_cols = size(output_array, 1)
    num_rows = size(output_array, 2)
    allocate(uniform_random(num_cols))

    max_value = sqrt(0.5_dp) * exp(-0.5_dp)

    do row = 1, num_rows
        accept_count = 0

        do while (accept_count < num_cols)
            call random_number(uniform_random)

            do j = 1, num_cols
                if (accept_count >= num_cols) exit

                x = xmax * uniform_random(j)
                y = sqrt(x) * exp(-x)

                call random_number(uniform_random(j))
                if (uniform_random(j) * max_value <= y) then
                    accept_count = accept_count + 1
                    output_array(accept_count, row) = x
                endif
            enddo
        enddo
    enddo

end subroutine generate_distribution_sqrt_x_exp_neg_x

subroutine generate_marker_distribution(xmax, output_array)

    use binsrc_mod, only: binsrc

!for a velocity distribution according to (1+v^7)*v^2*exp(-v^2),
!we have an energy distribution according to (1+E^(7/2))*sqrt(E)*exp(-E)

    real(dp), dimension(:,:), intent(out) :: output_array
    real(dp), intent(in) :: xmax
    integer, parameter :: nsize=10000
    logical, save :: firstentry = .true.
    real(dp), dimension(0:nsize), save :: pdf
    integer :: i,j, k, num_cols, num_rows
    real(dp) :: x,xi
    real(dp), save :: hx

    if(firstentry) then
        hx=xmax/dble(nsize)
        firstentry = .false.
        pdf(0)=0.d0
        do i=1,nsize
            x=(dble(i)-0.5d0)*hx
            pdf(i)=pdf(i-1)+(1.d0+x**7)*x**2*exp(-x**2)
        enddo
        pdf=pdf/pdf(nsize)
    endif

    num_cols = size(output_array, 1)
    num_rows = size(output_array, 2)

    do i = 1,num_cols
        do j = 1, num_rows
            call random_number(xi)
            call binsrc(pdf,0,nsize,xi,k)
            output_array(i,j) = (dble(k)-0.5d0)*hx
        enddo
    enddo

end subroutine generate_marker_distribution

subroutine generate_distribution_one_minus_x2(output_array)

    use binsrc_mod, only: binsrc

    real(dp), dimension(:,:), intent(out) :: output_array

    integer, parameter :: nsize = 10000
    logical, save :: firstentry = .true.
    real(dp), dimension(0:nsize), save :: cdf
    real(dp), save :: hx

    integer :: i, j, k, num_cols, num_rows
    real(dp) :: xi, x, fx

    if (firstentry) then
        hx = 1.0_dp / dble(nsize)
        firstentry = .false.
        cdf(0) = 0.0_dp
        do i = 1, nsize
            x  = (dble(i) - 0.5_dp) * hx
            fx = 1.0_dp - x*x
            cdf(i) = cdf(i-1) + fx
        enddo
        cdf = cdf / cdf(nsize)
    endif

    num_cols = size(output_array, 1)
    num_rows = size(output_array, 2)

    do i = 1, num_cols
        do j = 1, num_rows
            call random_number(xi)
            call binsrc(cdf, 0, nsize, xi, k)
            output_array(i, j) = (dble(k) - 0.5_dp) * hx
        enddo
    enddo

end subroutine generate_distribution_one_minus_x2

end module utils_scef_particle_init_mod
