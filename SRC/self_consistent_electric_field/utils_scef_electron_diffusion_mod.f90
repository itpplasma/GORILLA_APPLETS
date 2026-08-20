module utils_scef_electron_diffusion_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

contains

subroutine calc_electron_diffusion_coefficients !call this before the first ion pushing

    use gorilla_applets_types_mod, only: in, dc, start, s, g, exit_data
    use tetra_grid_settings_mod, only: grid_size
    use tetra_grid_mod, only: verts_sthetaphi
    use utils_data_pre_and_post_processing_mod, only: prepare_next_round_of_parallelised_particle_pushing, &
    calc_collision_coefficients_for_all_tetrahedra, normalise_prism_moments_and_prism_moments_squared, initialize_exit_data
    use llsq_mod, only: llsq
    use tetra_physics_mod, only: particle_mass,tetra_physics
    use russian_roulette_mod, only: prepare_russian_roulette
    use tetra_grid_settings_mod, only: sfc_s_min
    use utils_scef_particle_init_mod, only: allocate_start_type, set_particle_type_specifications, set_starting_positions, &
                                            set_rest_of_individual_particle_specifications, set_weight
    use utils_scef_particle_pushing_mod, only: parallelised_particle_pushing

    integer :: ns, i, n_particles, file_id, n_ignored, ns_start
    real(dp) :: extrapolation_factor, A, B, offset, tau_c_ei, v0, v_max, weights_before_redistribution
    real(dp) :: standard_deviation, dummy
    real(dp), dimension(:,:,:), allocatable :: rand_matrix
    character(len=100) :: filename, ns_str
    real(dp), dimension(:), allocatable :: data_for_diffusion_coefficient, data_for_convection_coefficient
    logical :: ignore_condition

    s%n_particles = 40000
    call random_number(dummy)

    if (.not.allocated(rand_matrix)) allocate(rand_matrix(5,s%n_particles,1))
    if (.not.allocated(dc%s_vertices)) allocate(dc%s_vertices(grid_size(1)+1))
    if (.not.allocated(dc%A)) allocate(dc%A(grid_size(1)+1)) !diffusion coefficients will be computed on every flux layer (including inner and outer boundaries)
    if (.not.allocated(dc%B)) allocate(dc%B(grid_size(1)+1)) !diffusion coefficients will be computed on every flux layer (including inner and outer boundaries)
    if (.not.allocated(dc%grad_A)) allocate(dc%grad_A(grid_size(1)))
    if (.not.allocated(dc%grad_B)) allocate(dc%grad_B(grid_size(1)))
    dc%s_vertices = verts_sthetaphi(1, [(grid_size(3)*(ns-1)+1,ns=1,grid_size(1)+1)])
    s%k = 1000
    s%j = 100
    if (.not.allocated(s%delta_s)) allocate(s%delta_s(s%k))
    if (.not.allocated(s%delta_s_squared)) allocate(s%delta_s_squared(s%k))
    if (.not.allocated(s%time)) allocate(s%time(s%k))
    if (.not.allocated(s%f_v)) allocate(s%f_v(s%k,s%j))
    if (.not.allocated(s%check)) allocate(s%check(s%k))
    if (.not.allocated(data_for_diffusion_coefficient)) allocate(data_for_diffusion_coefficient(s%n_particles))
    if (.not.allocated(data_for_convection_coefficient)) allocate(data_for_convection_coefficient(s%n_particles))
    if (.not.allocated(s%boole_large_distance)) allocate(s%boole_large_distance(s%n_particles))
    s%boole_large_distance = .false.

    call allocate_start_type(s%n_particles)
    call set_particle_type_specifications
    start%v0(2) = start%v0(2)*sqrt(s%temperature/in%energy_eV)
    call initialize_exit_data(s%n_particles)

    tau_c_ei = 1.7_dp*1.0d-4 !rough estimate from nrl formula booklet
    start%t(2) = 2.0d0*tau_c_ei !check afterwards if this was too little time

    s%time = [(start%t(2)/s%k*i,i = 1,s%k)]

    ns_start = 15

    do ns = ns_start,grid_size(1)

        print*, 'ns = ', ns, '/', grid_size(1)


        !Initiate electrons at the different flux surfaces leaving out the boundaries
        s%s0 = dc%s_vertices(ns)
        s%s0 = 0.5d0

        call RANDOM_NUMBER(rand_matrix)
        call set_starting_positions(rand_matrix,(/2/), s%s0)
        call set_rest_of_individual_particle_specifications(rand_matrix,boole_diffusion_coefficient_in = .true., &
                                                            species_in=(/2/),n_particles_in=s%n_particles)
        call set_weight


        call prepare_next_round_of_parallelised_particle_pushing(2)

        if (ns.eq.ns_start) then
            call calc_collision_coefficients_for_all_tetrahedra(2)

            v0 = start%v0(2)
            v_max = v0*sqrt(start%epsilon_max)
            weights_before_redistribution = g%total_volume*in%density
            call prepare_russian_roulette(v0,v_max,weights_before_redistribution,10)
        endif

        ! open(23,file = 'energy_before_pushing.dat')
        !     do i = 1,s%n_particles
        !         write(23,*) start%energy(i,2)
        !     enddo
        ! close(23)
        call parallelised_particle_pushing(species = 2,j = 1,boole_diffusion_coefficient = .true.,n_particles_in=s%n_particles)



        ! open(23,file = 'energy_after_pushing.dat')
        !     do i = 1,s%n_particles
        !         write(23,*) start%energy(i,2)
        !     enddo
        ! close(23)

        !Starting index (Throw away first 20 percent of values)
        i = ceiling(dble(s%k)*0.2_dp)


        call llsq ( int(s%k - i + 1, kind=8), s%time(i:s%k), s%delta_s(i:s%k), A, offset)

        data_for_diffusion_coefficient = s%delta_s_squared-2*A*s%time*s%delta_s+A**2*s%time**2
        call llsq ( int(s%k - i + 1, kind=8), s%time(i:s%k), data_for_diffusion_coefficient, B, offset)
        B = B/2

        ! n_ignored = 0
        ! data_for_convection_coefficient = exit_data%x(1,:,2)-s%s0
        ! do i = 1,s%n_particles
        !     ignore_condition = abs(data_for_convection_coefficient(i)).gt.0.99_dp*min(s%s0-sfc_s_min,1.0_dp-s%s0)
        !     if (ignore_condition) then
        !         data_for_convection_coefficient(i) = 0.0_dp
        !         n_ignored = n_ignored + 1
        !     endif
        ! enddo

        ! A = sum(exit_data%x(1,:,2)-s%s0)/(s%n_particles*start%t(2))
        ! data_for_diffusion_coefficient = abs(exit_data%x(1,:,2)-(s%s0+A*start%t(2)))
        ! call quicksort(data_for_diffusion_coefficient, 1, s%n_particles)
        ! standard_deviation = data_for_diffusion_coefficient(int(s%n_particles*0.682689_dp))
        ! B = standard_deviation**2/(2*start%t(2))


        dc%A(ns) = A
        dc%B(ns) = B

        print*, 'A = ', A
        print*, 'B = ', B
    write(ns_str, '(I0)') ns
    filename = 'exit_s_values' // trim(ns_str) // '.dat'
    open(newunit=file_id,file = filename)
        do i = 1,s%n_particles
            write(file_id,*) exit_data%x(1,i,2)
        enddo
    close(file_id)


    ! filename = 'data' // trim(ns_str) // '.dat'
    ! open(23,file = filename)
    !     do i = 1,s%k
    !         write(23,*) s%time(i), s%delta_s(i), s%delta_s_squared(i)
    !     enddo
    ! close(23)

    ! open(23,file = 'f_v.dat')
    !     do i = 1,s%j
    !         write(23,*) s%f_v(:,i)
    !     enddo
    ! close(23)

    !print*, weights%w(:,2)

    enddo

    !Extrapolate values to the boundaries
    extrapolation_factor = (dc%s_vertices(2)-dc%s_vertices(1))/(dc%s_vertices(3)-dc%s_vertices(2))
    dc%A(1) =     dc%A(2) + (dc%A(2)-dc%A(3))*extrapolation_factor
    dc%B(1) = max(dc%B(2) + (dc%B(2)-dc%B(3))*extrapolation_factor, 0.0_dp)
    extrapolation_factor = (dc%s_vertices(grid_size(1)+1)-dc%s_vertices(grid_size(1)))/ &
                           (dc%s_vertices(grid_size(1))-dc%s_vertices(grid_size(1)-1))
    dc%A(grid_size(1)+1) =     dc%A(grid_size(1)) + (dc%A(grid_size(1))-dc%A(grid_size(1)-1))*extrapolation_factor
    dc%B(grid_size(1)+1) = max(dc%B(grid_size(1)) + (dc%B(grid_size(1))-dc%B(grid_size(1)-1))*extrapolation_factor, 0.0_dp)

    !Calculate gradients
    do ns = 1,grid_size(1)
        dc%grad_A(ns) = (dc%A(ns+1)-dc%A(ns))/(dc%s_vertices(ns+1)-dc%s_vertices(ns))
        dc%grad_B(ns) = (dc%B(ns+1)-dc%B(ns))/(dc%s_vertices(ns+1)-dc%s_vertices(ns))
    enddo

    open(newunit=file_id,file = 'A_and_B.dat')
    do i = 1,grid_size(1)+1
        write(file_id,*) dc%A(i), dc%B(i)
    enddo
    close(file_id)

end subroutine calc_electron_diffusion_coefficients

subroutine calc_electron_density_via_random_walk(iteration_step) !call this after every ion pushing sequence

    use gorilla_applets_types_mod, only: in, time_t, dc, ep, g, output, start, exit_data, weights, filenames
    use tetra_grid_settings_mod, only: grid_size, sfc_s_min
    use binsrc_mod, only: binsrc
    use tetra_grid_mod, only: ntetr, verts_sthetaphi
    use utils_scef_electric_potential_mod, only: fill_vector_parts_with_value

    real(dp) :: delta_x, delta_t, xi, A, B, cell_size, B_fit, A_fit, p
    integer :: i, ns, k, num_steps_min, count_lost_particles, num_particles, particle_multiplication, iteration_step
    integer :: density_unit
    real(dp), dimension(:), allocatable :: electron_density, electron_prism_densities
    real(dp), dimension(:), allocatable :: position, weight, exit_time
    type(time_t) :: t
    logical :: boole_lost, boole_use_fit_function = .true.

    particle_multiplication = max(int(1.0d4/in%num_particles),1)
    num_particles = in%num_particles*particle_multiplication

    allocate(electron_density(grid_size(1)), electron_prism_densities(ntetr/3))
    electron_density = 0.0_dp

    !>This should later be an input
    allocate(position(num_particles))
    allocate(weight(num_particles))
    allocate(exit_time(num_particles))

    do i = 1,particle_multiplication
        position(i:num_particles:particle_multiplication) = start%x(1,:,2)
        weight(i:num_particles:particle_multiplication) = weights%w(:,2)
    enddo

    ! call RANDOM_NUMBER(position)
    ! call generate_distribution_one_minus_x2(position)
    ! allocate(dc%A(grid_size(1)+1)) !diffusion coefficients will be computed on every flux layer (including inner and outer boundaries)
    ! allocate(dc%B(grid_size(1)+1)) !diffusion coefficients will be computed on every flux layer (including inner and outer boundaries)
    ! allocate(dc%grad_A(grid_size(1)))
    ! allocate(dc%grad_B(grid_size(1)))
    ! allocate(dc%s_vertices(grid_size(1)+1))
    ! dc%s_vertices = verts_sthetaphi(1, [(grid_size(3)*(ns-1)+1,ns=1,grid_size(1)+1)])
    ! dc%A = 0.0_dp
    ! dc%B = 0.1_dp
    ! dc%grad_A = 0.0_dp
    ! dc%grad_B = 0.0_dp


    if (iteration_step.eq.1) then
        call get_diffusion_coefficient_data(boole_use_fit_function)
    else
        call calc_convection_coefficient_from_electric_field(boole_use_fit_function)
    endif

    t%step = in%time_step
    num_steps_min = 1000

count_lost_particles = 0
    do i = 1,num_particles
        t%confined = 0.0_dp
        boole_lost = .false.
        do while ((t%confined.lt.t%step).and.(.not.boole_lost))

            !binsrc finds k such that dc%s_vertices(k-1) < position(i) < dc%s_vertices(k)
            call binsrc(dc%s_vertices,1,grid_size(1)+1,position(i),k)
            k = k-1
            cell_size = dc%s_vertices(k+1) - dc%s_vertices(k)
            A = dc%A(k) + dc%grad_A(k)*(position(i)-dc%s_vertices(k))
            B = dc%B(k) + dc%grad_B(k)*(position(i)-dc%s_vertices(k))

            if (boole_use_fit_function) then
                B_fit = sum(dc%polynomial_coefficients_for_B*(/1.0_dp,position(i),position(i)**2,0.0_dp/))
                A_fit = sum(dc%polynomial_coefficients_for_A*(/1.0_dp,position(i),0.0_dp/))
                if (position(i).lt.dc%s_vertices(4)) then
                    p = (position(i) - sfc_s_min)/(dc%s_vertices(4) - sfc_s_min)
                    B = p*B_fit + (1-p)*B
                    A = p*A_fit + (1-p)*B/position(i) + A
                else
                    B = B_fit
                    A = A + A_fit
                endif
            endif


            delta_t = min(t%step/num_steps_min, cell_size**2/((abs(A)+abs(B))*4), t%step - t%confined) !control maximum possible jump

            electron_density(k) = electron_density(k) + delta_t*weight(i)
            t%confined = t%confined + delta_t

            call random_number(xi)
            xi = sqrt(12.0_dp)*(xi-0.5_dp)
            delta_x = sqrt(2*delta_t*B)*xi + A*delta_t
            !print*, 'diffusive part over convective part = ', A,B,abs(sqrt(2*delta_t*B)/(A*delta_t)), &
            !sqrt(2*delta_t*B)*xi/(A*delta_t)

            position(i) = position(i) + delta_x

            if (position(i).lt.sfc_s_min) then
                position(i) = 2.0_dp*sfc_s_min - position(i)
            endif
            if (position(i).gt.1.0_dp)    then
                boole_lost = .true.
                count_lost_particles = count_lost_particles+1
            endif
        enddo
        exit_time(i) = t%confined
    enddo

    if (count_lost_particles.lt.num_particles) print*, 'Warning: the tracing time (', t%step,'s) was so short that only ',&
                                                        count_lost_particles,'out of', num_particles, &
                                                        'electrons left the computation domain'

    open(newunit=density_unit,file=filenames%electron_density)
    do ns = 1, grid_size(1)
        write(density_unit,*) sfc_s_min + (1.0_dp - sfc_s_min) * (ns - 0.5_dp) / (grid_size(1)+1), electron_density(ns)
    enddo
    close(density_unit)

    do i = 1,grid_size(1)
        electron_density(i) = electron_density(i)/(ep%s_shell_volumes(i)*t%step*num_particles)
        call fill_vector_parts_with_value(electron_prism_densities, g%prisms_per_flux_tube(i,:), electron_density(i))
    enddo


    !This is a bit cumbersome having to use electron_prism_densities and a loop over all prisms instead of simply putting
    !output%prism_moments into fill_vector_parts_with_value, but that does not work because the latter does not accept
    !a complex argument. Think about a solution to this problem
    do i = 1,ntetr/3
        output%prism_moments(1,i,2) = complex(electron_prism_densities(i), 0.0_dp)
    enddo

    print*, 'Total tracing time of all electrons divided by number of particles is: ', &
    sum(exit_time)/num_particles, 's'
    !print*, 'exit times are : ', exit_data%t_confined(:,2)

end subroutine calc_electron_density_via_random_walk

subroutine get_diffusion_coefficient_data(boole_use_fit_function)

    use gorilla_applets_types_mod, only: dc
    use tetra_grid_settings_mod, only: grid_size
    use tetra_grid_mod, only: verts_sthetaphi
    use utils_polyfit_mod, only: quadratic_fit

    character(len=100) :: filename
    integer :: id, ns, lower_bound, upper_bound
    real(dp), dimension(:), allocatable :: x,y
    logical :: success, boole_use_fit_function
    real(dp) :: a0, a1, a2


    if (.not.allocated(dc%A)) allocate(dc%A(grid_size(1)+1))
    if (.not.allocated(dc%A_from_first_run)) allocate(dc%A_from_first_run(grid_size(1)+1))
    if (.not.allocated(dc%B)) allocate(dc%B(grid_size(1)+1))
    if (.not.allocated(dc%grad_A)) allocate(dc%grad_A(grid_size(1)))
    if (.not.allocated(dc%grad_B)) allocate(dc%grad_B(grid_size(1)))
    if (.not.allocated(dc%s_vertices)) allocate(dc%s_vertices(grid_size(1)+1))

    filename = 'A_and_B.dat'
    open(newunit=id,file = filename,action='read')
    do ns = 1,grid_size(1)+1
        read(id,*) dc%A(ns), dc%B(ns)
    enddo

    dc%s_vertices = verts_sthetaphi(1, [(grid_size(3)*(ns-1)+1,ns=1,grid_size(1)+1)])

    dc%A_from_first_run = dc%A
    if (boole_use_fit_function) dc%A = 0.0_dp

    do ns = 1,grid_size(1)
        dc%grad_A(ns) = (dc%A(ns+1)-dc%A(ns))/(dc%s_vertices(ns+1)-dc%s_vertices(ns))
        dc%grad_B(ns) = (dc%B(ns+1)-dc%B(ns))/(dc%s_vertices(ns+1)-dc%s_vertices(ns))
    enddo



    lower_bound = 2
    upper_bound = grid_size(1) - 4

    if (.not.allocated(x)) allocate(x(upper_bound-lower_bound+1))
    if (.not.allocated(y)) allocate(y(upper_bound-lower_bound+1))

    x = dc%s_vertices(lower_bound:upper_bound)
    y = dc%B(lower_bound:upper_bound)

    call quadratic_fit(size(x), x, y, a0, a1, a2, success) ! Least-squares fit of y ≈ a2*x^2 + a1*x + a0

    dc%polynomial_coefficients_for_B = (/a0,a1,a2,0.0_dp/)
    dc%polynomial_coefficients_for_A = (/a1,2*a2,0.0_dp/)

end subroutine get_diffusion_coefficient_data

subroutine calc_convection_coefficient_from_electric_field(boole_use_fit_function)

    use tetra_physics_mod, only: tetra_physics
    use constants, only: echarge, ev2erg
    use gorilla_applets_types_mod, only: dc, in, s
    use tetra_grid_settings_mod, only: grid_size

    logical :: boole_use_fit_function
    integer :: i, ns
    real(dp) :: electric_field

    do i = 1,grid_size(1)+1
        if (i.eq.1) then
            electric_field = -tetra_physics(1)%gPhi(1)
        elseif (i.eq.grid_size(1)+1) then
            electric_field = -tetra_physics((grid_size(1)-1)*6*grid_size(2)+1)%gPhi(1)
        else
            electric_field = -0.5_dp*(tetra_physics((i-2)*6*grid_size(2)+1)%gPhi(1) + &
                                      tetra_physics((i-1)*6*grid_size(2)+1)%gPhi(1))
        endif
        dc%A(i) = -dc%B(i)*echarge*electric_field/(s%temperature*ev2erg)
        if (.not.boole_use_fit_function) dc%A(i) = dc%A(i) + dc%A_from_first_run(i)
    enddo

    do ns = 1,grid_size(1)
        dc%grad_A(ns) = (dc%A(ns+1)-dc%A(ns))/(dc%s_vertices(ns+1)-dc%s_vertices(ns))
        dc%grad_B(ns) = (dc%B(ns+1)-dc%B(ns))/(dc%s_vertices(ns+1)-dc%s_vertices(ns))
    enddo

end subroutine calc_convection_coefficient_from_electric_field

subroutine sort_array(array)
  real(dp), dimension(:) :: array
  integer, parameter :: n = 6
  real :: a(n) = [3.5, 1.2, -4.0, 7.8, 0.0, 2.1]
  integer :: i, j
  real(dp) :: temp

  ! simple bubble sort (ascending)
  do i = 1, size(array)
     do j = i+1, size(array)
        if (array(j) < array(i)) then
           temp = array(i)
           array(i) = array(j)
           array(j) = temp
        end if
     end do
  end do

end subroutine sort_array

recursive subroutine quicksort(arr, left, right)
real(dp), intent(inout) :: arr(:)
integer, intent(in) :: left, right
integer :: i, j
real(dp) :: pivot, temp

i = left
j = right
pivot = arr((left+right)/2)

do
    do while (arr(i) < pivot)
        i = i + 1
    end do
    do while (arr(j) > pivot)
        j = j - 1
    end do
    if (i <= j) then
        temp = arr(i)
        arr(i) = arr(j)
        arr(j) = temp
        i = i + 1
        j = j - 1
    end if
    if (i > j) exit
end do

if (left < j) call quicksort(arr, left, j)
if (i < right) call quicksort(arr, i, right)
end subroutine quicksort

end module utils_scef_electron_diffusion_mod
