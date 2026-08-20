module utils_scef_particle_init_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

contains

subroutine calc_starting_conditions

    use gorilla_applets_types_mod, only: in, start
    use tetra_grid_settings_mod, only: sfc_s_min

    real(dp), dimension(:,:,:), allocatable                :: rand_matrix

    allocate(rand_matrix(5,in%num_particles,in%n_species))
    call RANDOM_NUMBER(rand_matrix)

    call allocate_start_type
    call allocate_weights
    call set_particle_type_specifications
    call set_starting_positions(rand_matrix)!,s0=sfc_s_min*1.1_dp)
    call set_rest_of_individual_particle_specifications(rand_matrix)

end subroutine calc_starting_conditions

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
    allocate(start%lost(n_particles,in%n_species))
    allocate(start%particle_charge(in%n_species))
    allocate(start%particle_mass(in%n_species))
    allocate(start%cm_over_e(in%n_species))
    allocate(start%t(in%n_species))
    allocate(start%v0(in%n_species))

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

    real(dp), dimension(:,:,:), intent(in) :: rand_matrix
    integer, dimension(:), intent(in), optional :: species_in
    logical, intent(in), optional :: boole_diffusion_coefficient_in
    logical :: boole_diffusion_coefficient=.false.
    integer, dimension(:), allocatable :: species
    integer :: i
    integer, intent(in), optional :: n_particles_in
    integer :: n_particles
    real(dp), dimension(:,:), allocatable :: radial_transport_energies
    real(dp) :: temperature

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

    if (present(boole_diffusion_coefficient_in)) boole_diffusion_coefficient = boole_diffusion_coefficient_in

    start%pitch(:,species) = 2*rand_matrix(4,:,:)-1 !pitch parameter
    start%energy(:,species) = in%energy_eV
    if ((.not. in%boole_monoenergetic).or.boole_diffusion_coefficient) then
        !start%energy(:,species) = start%epsilon_max*in%energy_eV*rand_matrix(5,:,:) !boltzmann energy distribution
        temperature = in%energy_eV
        if (boole_diffusion_coefficient) temperature = s%temperature
        allocate(radial_transport_energies(n_particles,size(species)))
        call generate_distribution_sqrt_x_exp_neg_x(start%epsilon_max,radial_transport_energies)
        !call generate_marker_distribution(start%epsilon_max,radial_transport_energies)
        start%energy(:,species) = temperature*radial_transport_energies !boltzmann energy distribution
    endif

    if (in%boole_antithetic_variate) then
        start%x(:,1:n_particles:2,species) = start%x(:,2:n_particles:2,species)
        start%pitch(1:n_particles:2,species) = -start%pitch(2:n_particles:2,species)
        start%energy(1:n_particles:2,species) = start%energy(2:n_particles:2,species)
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
    end select

    start%particle_charge = (/charge, -echarge/)
    start%particle_mass = (/mass, ame/)
    start%cm_over_e = (/cm_over_e, -clight*ame/echarge/)
    start%t = (/in%time_step, in%time_step/) !/42.0_dp
    if (in%boole_static_ne) start%t(2) = 0.0_dp

    start%v0 = sqrt(2.0_dp*in%energy_eV*ev2erg/start%particle_mass)
    start%epsilon_max = 16.0_dp

    call set_weight

end subroutine set_particle_type_specifications

subroutine set_weight

    use gorilla_applets_types_mod, only: in, weights
    use tetra_grid_settings_mod, only: sfc_s_min, n_field_periods
    use constants, only: pi

    weights%w = in%density*(1-sfc_s_min)*4*pi**2/n_field_periods

end subroutine set_weight

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
