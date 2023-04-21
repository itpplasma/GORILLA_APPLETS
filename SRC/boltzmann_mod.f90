module boltzmann_mod
    implicit none
!
    private
!
    double precision :: time_step,energy_eV, max_poloidal_flux
    integer, dimension(:,:), allocatable :: tetra_indices_per_prism
    double precision, dimension(:), allocatable :: prism_volumes, refined_prism_volumes, elec_pot_vec, n_b
    double precision, dimension(:,:), allocatable :: verts, sqrt_g, r_integrand_constants
    complex, dimension(:,:,:), allocatable :: tetr_moments, prism_moments
    complex, dimension(:,:,:,:), allocatable :: moments_in_frequency_space, prism_moments_statistics
    integer :: i_integrator_type, seed_option, n_moments, n_species, n_prisms, num_particles, &
               & ind_a, ind_b, ind_c, n_pushings, counter_phi_0_mappings, n_simulations
    double precision :: n_particles, density, constant_part_of_weights
    integer, dimension(4) :: moments_selector
    complex, dimension(:,:), allocatable :: weights
    double precision, dimension(:), allocatable :: J_perp, poloidal_flux
    logical :: boole_random_precalc, boole_refined_sqrt_g, boole_boltzmann_energies
    character(1024) :: filename_dwell_times, filename_starting_conditions, filename_vertex_coordinates, &
    & filename_vertex_indices
    integer :: n_fourier_modes, n_triangles
    logical :: boole_linear_density_simulation = .true., boole_antithetic_variate = .true.
!
    !Namelist for boltzmann input
    NAMELIST /boltzmann_nml/ time_step,energy_eV,n_particles,n_simulations,density,boole_refined_sqrt_g,boole_boltzmann_energies, &
    & i_integrator_type,seed_option,boole_random_precalc,filename_dwell_times,filename_starting_conditions, &
    & filename_vertex_coordinates, filename_vertex_indices
!
    public :: calc_boltzmann
!    
contains
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine load_boltzmann_inp()
!
    open(unit=71, file='boltzmann.inp', status='unknown')
    read(71,nml=boltzmann_nml)
    close(71)
!    
    print *,'GORILLA: Loaded input data from boltzmann.inp'
end subroutine load_boltzmann_inp
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine calc_starting_conditions(vmod,start_pos_pitch_mat)
!
    use constants, only: pi, ev2erg
    use tetra_grid_mod, only: verts_rphiz, verts_sthetaphi, ntetr
    use find_tetra_mod, only: find_tetra
    use tetra_physics_mod, only: coord_system
!
    implicit none
    double precision, intent(in)                                   :: vmod
    double precision, dimension(:,:), allocatable, intent(out)     :: start_pos_pitch_mat
    double precision                                               :: rand_scalar, vpar, vperp
    double precision                                               :: amin, amax, cmin, cmax
    double precision, dimension(:), allocatable                    :: rand_vector
    double precision, dimension(:,:), allocatable                  :: rand_matrix
    integer                                                        :: i
    double precision, dimension(3)                                 :: x
    integer                                                        :: ind_tetr_out,iface
!
!!!!comment out the following section to make starting conditions really random!!!
!
    integer,dimension(:), allocatable                              :: seed
    integer                                                        :: n,j
!
    open(unit = 85, file='seed.inp', status='old',action = 'read')
    read(85,*) n
    allocate(seed(n))
    read(85,*) seed
    close(85)
    CALL RANDOM_SEED (PUT=seed)
    deallocate(seed)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    allocate(start_pos_pitch_mat(5,num_particles*n_simulations))
    allocate(rand_vector(num_particles*n_simulations))
    allocate(rand_matrix(5,num_particles*n_simulations))
!
    start_pos_pitch_mat = 0
!
    ind_a = 1 !(R in cylindrical coordinates, s in flux coordinates)
    ind_b = 2 !(phi in cylindrical and flux coordinates)
    ind_c = 3 !(z in cylindrical coordinates, theta in flux coordinates)
!
    if (coord_system.eq.2) then
        ind_b = 3
        ind_c = 2
    endif
!
    if (coord_system.eq.1) allocate(verts(size(verts_rphiz(:,1)),size(verts_rphiz(1,:))))
    if (coord_system.eq.2) allocate(verts(size(verts_sthetaphi(:,1)),size(verts_sthetaphi(1,:))))
    if (coord_system.eq.1) verts = verts_rphiz
    if (coord_system.eq.2) verts = verts_sthetaphi
!
    amin = minval(verts(ind_a,:))
    amax = maxval(verts(ind_a,:))
    cmin = minval(verts(ind_c,:))
    cmax = maxval(verts(ind_c,:))
!
    constant_part_of_weights = density*(amax-amin)*(cmax-cmin)*2*pi
!
    !compute starting conditions
    call RANDOM_NUMBER(rand_matrix)
    start_pos_pitch_mat(ind_a,:) = amin + (amax - amin)*rand_matrix(ind_a,:) !r in cylindrical, s in flux coordinates
    start_pos_pitch_mat(ind_b,:) = 2*pi*rand_matrix(ind_b,:) !phi in cylindrical and flux coordinates
    start_pos_pitch_mat(ind_c,:) = cmin + (cmax - cmin)*rand_matrix(ind_c,:) !z in cylindrical, theta in flux coordinates
    start_pos_pitch_mat(4,:) = 2*rand_matrix(4,:)-1 !pitch parameter
!
    if (boole_boltzmann_energies) then !compare with equation 133 of master thesis of Jonatan Schatzlmayr (remaining parts will be added later)
        start_pos_pitch_mat(5,:) = 5*energy_eV*rand_matrix(5,:) !boltzmann energy distribution
        constant_part_of_weights = constant_part_of_weights*10/sqrt(pi*energy_eV*ev2erg)
    endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! start antithetic variate !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (boole_antithetic_variate) then
        start_pos_pitch_mat(:,1:num_particles*n_simulations:2) = start_pos_pitch_mat(:,2:num_particles*n_simulations:2)
        start_pos_pitch_mat(4,1:num_particles*n_simulations:2) = -start_pos_pitch_mat(4,2:num_particles*n_simulations:2)
    endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! end antithetic variate !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    weights(:,1) = constant_part_of_weights
    if (boole_refined_sqrt_g.eqv..false.) weights(:,1) = constant_part_of_weights*start_pos_pitch_mat(ind_a,:)
!
end subroutine calc_starting_conditions
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine calc_boltzmann
!
    use orbit_timestep_gorilla_mod, only: initialize_gorilla
    use constants, only: ev2erg, pi
    use tetra_physics_mod, only: particle_mass,particle_charge,cm_over_e,mag_axis_R0, coord_system, tetra_physics
    use fluxtv_mod, only: load_flux_tube_volume,pos_fluxtv_mat
    use omp_lib, only: omp_get_thread_num, omp_get_num_threads
    use parmot_mod, only: rmu,ro0
    use velo_mod, only: isw_field_type
    use supporting_functions_mod, only: theta_sym_flux2theta_vmec,theta_vmec2theta_sym_flux
    use tetra_grid_settings_mod, only: n_field_periods
    use tetra_grid_mod, only: ntetr, nvert, verts_rphiz, tetra_grid, verts_sthetaphi
    use gorilla_settings_mod, only: boole_array_optional_quantities
!
    implicit none
!
    double precision, dimension(:,:), allocatable :: start_pos_pitch_mat
    double precision :: vmod,pitchpar,vpar,vperp,t_remain,t_confined,tau_out_can, velocity
    integer :: kpart,i,j,n,l,m,k,ind_tetr,iface,n_lost_particles,ierr
    integer :: n_start, n_end, i_part
    double precision, dimension(3) :: x_rand_beg,x
    logical :: boole_initialized,boole_particle_lost
    double precision :: dtau, dphi,dtaumin
    double precision, dimension(5) :: z
    Character(LEN=50) :: format_moments, format_fourier_moments
!
    ! open(35, file = 'outliers.dat')
    ! close(35,status='delete')
    !Load input for boltzmann computation
    call load_boltzmann_inp()
!
    !call read_in_starting_conditions(start_pos_pitch_mat, num_particles)
!
    num_particles = int(n_particles)
    n_start = 1
    n_end = num_particles
    n_species = 1
    n_fourier_modes = 5
!
    !prepare moment calculation
    n_moments = 0
    moments_selector = 0
    do i = 1,size(boole_array_optional_quantities)
        if (boole_array_optional_quantities(i).eqv..true.) then
            n_moments = n_moments + 1
            moments_selector(n_moments) = i
        endif
    enddo
!
    !Initialize GORILLA
    call initialize_gorilla()
!
    !compute maximum poloidal flux
    max_poloidal_flux = 0
    do i = 1, ntetr
        max_poloidal_flux = maxval((/max_poloidal_flux,tetra_physics(i)%Aphi1/))
    enddo
!
    n_prisms = ntetr/3
    allocate(tetr_moments(n_moments,ntetr,n_species))
    allocate(prism_moments(n_moments,n_prisms,n_species))
    allocate(prism_moments_statistics(n_moments,n_prisms,n_species,2)) !only sensible if n_simulations >> 1
    allocate(weights(num_particles*n_simulations,1))
    allocate(J_perp(num_particles*n_simulations))
    allocate(poloidal_flux(num_particles*n_simulations))
    allocate(elec_pot_vec(n_prisms))
    allocate(n_b(n_prisms))
!
    prism_moments_statistics = 0
    poloidal_flux = 0
    elec_pot_vec = 0
    n_b = 0 !boltzmann density for uniform spatial distribution
!
    !Load fluxtube volume for a starting position (File to be chosen in gorilla_applets.inp)
    !call load_flux_tube_volume()
!
    call calc_square_root_g
!
print*, 'start calc_volume_integrals'
    call calc_volume_integrals
print*, 'end calc_volume_integrals'
!
    !Compute velocity module from kinetic energy dependent on particle species
    vmod=sqrt(2.d0*energy_eV*ev2erg/particle_mass)
    !
print*, 'calc_starting_conditions started'
    call calc_starting_conditions(vmod,start_pos_pitch_mat)
print*, 'calc_starting_conditions finished'
    !
    !------------------------------------------------------------------------------------------------------------!
    !------------------------------------ Initialization of direct integrator -----------------------------------!
    !
                if(i_integrator_type.eq.2) then
    !
                    !inverse relativistic temperature
                    rmu=1d8
    !
                    !normalized larmor radius
                    ro0 = vmod*cm_over_e
    !
                    isw_field_type=1
    
                    !normalized slowing down time:
                    dtau = -1.d0*time_step*vmod
    !
                    !field line integration step step over phi (to check chamber wall crossing)
                    dphi=2.d0*pi/dble(100)
    !
                    !orbit integration time step (to check chamber wall crossing)
                    dtaumin=dphi*mag_axis_R0
    !
                endif
    !
    !------------------------------------------------------------------------------------------------------------!
    !
    !
    do j = 0,n_simulations-1
        kpart = 0
        n_lost_particles = 0
        n_pushings = 0
        counter_phi_0_mappings = 0
        tetr_moments = 0
        prism_moments = 0
!
        !$OMP PARALLEL DEFAULT(NONE) &
        !$OMP& SHARED(num_particles,pos_fluxtv_mat,kpart,vmod,time_step,i_integrator_type, &
        !$OMP& dtau,dtaumin,n_start,n_end, j, &
        !$OMP& start_pos_pitch_mat, particle_mass, boole_boltzmann_energies) &
        !$OMP& PRIVATE(n,boole_particle_lost,x_rand_beg,x,pitchpar,vpar,vperp,boole_initialized, &
        !$OMP& ind_tetr,iface,t_remain,t_confined,z,ierr,tau_out_can, velocity) &
        !$OMP& REDUCTION(+:n_lost_particles,tetr_moments, n_pushings, counter_phi_0_mappings)
print*, 'get number of threads', omp_get_num_threads()
        !$OMP DO
!
        !Loop over particles
        do n = n_start,n_end !1,num_particles
!
            !$omp critical
            !Counter for particles
            kpart = kpart+1 !in general not equal to n becuase of parallelisation
            !boole_particle_lost = .false.
if (n_end.gt.10) then
if (modulo(kpart,int(n_end/10)).eq.0) then
print *, kpart, ' / ', num_particles, 'particle: ', n, 'thread: ' !, omp_get_thread_num()
endif
else 
print *, kpart, ' / ', num_particles, 'particle: ', n, 'thread: ' !, omp_get_thread_num()
endif
            !$omp end critical
!print*, kpart
            !You need x_rand_beg(1,3), pitchpar(1) (between -1 and 1), energy is already given
            x_rand_beg = start_pos_pitch_mat(1:3,n+num_particles*j)
            pitchpar = start_pos_pitch_mat(4,n+num_particles*j)
!
            select case(i_integrator_type)
                case(1,0)
                    x = x_rand_beg
                case(2)
                    z(1) = x_rand_beg(1)
                    z(3) = x_rand_beg(3)
                    z(2) = theta_sym_flux2theta_vmec(z(1),x_rand_beg(2),z(3))  !Transform theta_symflux to theta_vmec
                    z(4) = 1.d0
            end select

            ! PRINT*, 'x_rand_beg = ', x_rand_beg
            ! PRINT*, 'pitchpar = ', pitchpar
!
            !Particle velocities in accordance with integrator
            select case(i_integrator_type)
                case(1,0)
                    vpar = pitchpar * vmod
                    vperp = sqrt(vmod**2-vpar**2)
!PRINT*, 'parallel velocity = ', vpar
                    if (boole_boltzmann_energies) then
                        velocity = sqrt(start_pos_pitch_mat(5,n+num_particles*j)*ev2erg*2/particle_mass)
                        vpar = pitchpar * velocity
                        vperp = sqrt(velocity**2-vpar**2)
                    endif
                case(2)
                    z(5) = pitchpar
            end select
!
            !Orbit integration
!
            select case(i_integrator_type)
!
                case(0)
                    !No orbit computation
!
                case(1)
                    boole_initialized = .false.
                    call orbit_timestep_gorilla_boltzmann(x,vpar,vperp,time_step,boole_initialized, &
                    & ind_tetr,iface,t_remain,n,vmod,start_pos_pitch_mat,j)
!
                    !Confinement time of particles
                    t_confined = time_step - t_remain
!
                    !Lost particle handling
                    if(ind_tetr.eq.-1) then
                        !$omp critical
                            n_lost_particles = n_lost_particles + 1
                            boole_particle_lost = .true.
                        !$omp end critical
                    endif
!
                    !Write results in file
                    ! !$omp critical
                    !     write(99,*) n, boole_particle_lost , x_rand_beg ,pitchpar,x(1),t_confined
                    ! !$omp end critical
!
                case(2)
                    ierr = 0
                    call orbit_timestep_can_boltzmann(z,dtau,dtaumin,ierr,tau_out_can)
!
                    if(ierr.eq.1) then
                        !$omp critical
                            n_lost_particles = n_lost_particles + 1
                            boole_particle_lost = .true.
                        !$omp end critical
                    endif
!
                    t_confined = -1.d0*tau_out_can/vmod
!
!                        !Write results in file
                    ! !$omp critical
                    !     write(99,*) n, boole_particle_lost , x_rand_beg ,pitchpar,z(1),t_confined
                    ! !$omp end critical
!
            end select
!
        enddo !n
        !$OMP END DO
        !$OMP END PARALLEL
!
print*, 'tetra moment', tetr_moments(:,1,1)
        prism_moments = (tetr_moments(:,tetra_indices_per_prism(:,1),:) + &
                        & tetr_moments(:,tetra_indices_per_prism(:,2),:) + &
                        & tetr_moments(:,tetra_indices_per_prism(:,3),:))
!
        do n = 1,n_moments
            do m = 1,n_species
                prism_moments(n,:,m) = prism_moments(n,:,m)/(prism_volumes*time_step*n_particles) !do normalisations
                if (boole_refined_sqrt_g) prism_moments(n,:,m) = prism_moments(n,:,m)*prism_volumes/refined_prism_volumes
            enddo
        enddo

        prism_moments_statistics(:,:,:,1) = prism_moments_statistics(:,:,:,1) + prism_moments
        prism_moments_statistics(:,:,:,2) = prism_moments_statistics(:,:,:,2) + prism_moments**2
print*, 'prism moments', prism_moments(4,1,1), prism_moments(4,1,1)**2
    enddo
!
    call fourier_transform_moments
!
    print*, 'Number of lost particles',n_lost_particles
    print*, 'max_poloidal_flux is', max_poloidal_flux
    print*, 'average number of pushings = ', n_pushings/n_particles
    print*, 'average number of toroidal revolutions = ', counter_phi_0_mappings/n_particles
    ! open(99,file='confined_fraction.dat')
    ! write(99,*) 1.d0-dble(n_lost_particles)/dble(num_particles)
!
    101 format(1000(e21.14,x))
    if (coord_system.eq.1) then
        ![R,phi,Z]: Write vertex coordinates to File
        open(55, file=filename_vertex_coordinates)
        do i=1, nvert
            write(55,101) verts_rphiz(1, i), verts_rphiz(2, i), verts_rphiz(3, i)
        end do
        close(55)
    elseif (coord_system.eq.2) then
        ![s,theta,phi]: Write vertex coordinates to File
        open(55, file=filename_vertex_coordinates)
        do i=1, nvert
            write(55,101) verts_sthetaphi(1, i), verts_sthetaphi(2, i), verts_sthetaphi(3, i)
        end do
        close(55)
    endif
!
    !Write vertex indices to File
    open(56, file=filename_vertex_indices)
    do i=1, ntetr
        write(56, *) tetra_grid(i)%ind_knot([1, 2, 3, 4])
    end do
    close(56)
!
    open(58, file = 'prism_volumes.dat')
    write(58,'(ES20.10E4)') prism_volumes
    close(58)
!
    if (n_moments.gt.0) then
        open(59, file = 'prism_moments.dat')
        do l = 1,n_prisms
            do i = 1,n_moments - 1
                write(59,'(2ES20.10E4)',advance="no") real(prism_moments(i,l,1)), aimag(prism_moments(i,l,1))
            enddo
                write(59,'(2ES20.10E4)') real(prism_moments(n_moments,l,1)), aimag(prism_moments(n_moments,l,1))
        enddo
        close(59)
        open(69, file = 'prism_moments_sums.dat')
        do l = 1,n_prisms
            do i = 1,n_moments - 1
                write(69,'(2ES20.10E4)',advance="no") real(prism_moments_statistics(i,l,1,1)), &
                                                    & aimag(prism_moments_statistics(i,l,1,1))
            enddo
                write(69,'(2ES20.10E4)') real(prism_moments_statistics(n_moments,l,1,1)), &
                                      & aimag(prism_moments_statistics(n_moments,l,1,1))
        enddo
        close(69)
        open(73, file = 'prism_moments_summed_squares.dat')
        do l = 1,n_prisms
            do i = 1,n_moments - 1
                write(73,'(2ES20.10E4)',advance="no") real(prism_moments_statistics(i,l,1,2)), &
                                                    & aimag(prism_moments_statistics(i,l,1,2))
            enddo
                write(73,'(2ES20.10E4)') real(prism_moments_statistics(n_moments,l,1,2)), &
                                      & aimag(prism_moments_statistics(n_moments,l,1,2))
        enddo
        close(73)
        open(60, file = 'tetr_moments.dat')
        do l = 1,ntetr
            do i = 1,n_moments - 1
                write(60,'(2ES20.10E4)',advance="no") real(tetr_moments(i,l,1)), aimag(tetr_moments(i,l,1))
            enddo
                write(60,'(2ES20.10E4)') real(tetr_moments(n_moments,l,1)), aimag(tetr_moments(n_moments,l,1))
        enddo
        close(60)
    endif
!
    open(61, file = 'tetra_indices_per_prism.dat')
    do i = 1,n_prisms
        write(61,*) tetra_indices_per_prism(i,:)
    enddo
    close(61)
!
    open(62, file = 'fourier_moments.dat')
    do l = 1,n_triangles
        do i = 1,n_fourier_modes-1
            write(62,'(2ES20.10E4)',advance="no") real(moments_in_frequency_space(1,l,i,1)), &
                                                    aimag(moments_in_frequency_space(1,l,i,1))!'(',real(moments_in_frequency_space(1,1,i,1)),',',aimag(moments_in_frequency_space(1,1,i,1)),')'
        enddo
            write(62,'(2ES20.10E4)') real(moments_in_frequency_space(1,l,n_fourier_modes,1)), &
                                       aimag(moments_in_frequency_space(1,l,n_fourier_modes,1))
    enddo   
    close(62)
!
    open(65, file = 'sqrt_g.dat')
    do i = 1,ntetr
        write(65,*) sqrt_g(i,:)
    enddo
    close(65)
!
    open(66, file = 'refined_prism_volumes.dat')
    write(66,'(ES20.10E4)') refined_prism_volumes
    close(66)
!
    open(67, file = 'r_integrand_constants.dat')
    do i = 1,n_prisms
        write(67,*) r_integrand_constants(i,:)
    enddo
    close(67)
!
    open(68, file = 'elec_pot_vec.dat')
    write(68,'(ES20.10E4)') elec_pot_vec
    close(68)
!
    open(70, file = 'boltzmann_densities.dat')
    write(70,'(ES20.10E4)') n_b
    close(70)
!
    open(72, file = 'h_phi.dat')
    do i = 1,n_prisms
        write(72,*) tetra_physics(i*3-2)%h2_1
    enddo
    close(72)
!
PRINT*, 'particle mass = ', particle_mass
! PRINT*, 'large radius = ', mag_axis_R0
! PRINT*, 'parallel velocity = ', vpar
PRINT*, 'absolute value of velocity = ', vmod
! PRINT*, 'perpendicular velocity = ', vperp
! PRINT*, 'pitch par =', pitchpar
PRINT*, 'particle charge = ', particle_charge
PRINT*, 'temperature = ', ev2erg*energy_eV
PRINT*, 'boole_refined_sqrt_g = ', boole_refined_sqrt_g
print*, prism_moments(4,1,1), real(prism_moments(4,1,1)), real(prism_moments(4,1,1))**2, 5.10106855d21**2, 1.0d21*1.0d21
!print*, tetra_physics(ntetr/2)%h1_1, tetra_physics(ntetr/2)%h2_1, tetra_physics(ntetr/2)%h3_1, tetra_physics(ntetr/2)%R1, &
    !sqrt(tetra_physics(ntetr/2)%h1_1**2+(tetra_physics(ntetr/2)%h2_1/tetra_physics(ntetr/2)%R1)**2+tetra_physics(ntetr/2)%h3_1**2)
!
deallocate(start_pos_pitch_mat, tetr_moments, prism_moments)
!
end subroutine calc_boltzmann
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine orbit_timestep_gorilla_boltzmann(x,vpar,vperp,t_step,boole_initialized,ind_tetr,iface, t_remain_out, n, &
                                            & vmod,start_pos_pitch_mat,j)
!
    use pusher_tetra_rk_mod, only: pusher_tetra_rk,initialize_const_motion_rk
    use pusher_tetra_poly_mod, only: pusher_tetra_poly,initialize_const_motion_poly
    use tetra_physics_poly_precomp_mod , only: make_precomp_poly_perpinv, initialize_boole_precomp_poly_perpinv, &
        & alloc_precomp_poly_perpinv
    use tetra_physics_mod, only: tetra_physics,particle_charge,particle_mass,cm_over_e
    use gorilla_settings_mod, only: ipusher, poly_order, optional_quantities_type, boole_array_optional_quantities
    use orbit_timestep_gorilla_mod, only: check_coordinate_domain
    use supporting_functions_mod, only: bmod_func, vperp_func
    use find_tetra_mod, only: find_tetra
    use constants, only: pi, ev2erg
    use tetra_grid_mod, only: tetra_grid
!
    implicit none
!
    double precision, dimension(3), intent(inout)   :: x
    double precision, intent(inout)                 :: vpar,vperp
    double precision, intent(in)                    :: t_step
    logical, intent(inout)                          :: boole_initialized
    integer, intent(inout)                          :: ind_tetr,iface
    double precision, intent(out), optional         :: t_remain_out
    double precision, dimension(3)                  :: z_save, x_save
    double precision                                :: vperp2,t_remain,t_pass,vpar_save, vmod
    logical                                         :: boole_t_finished
    integer                                         :: ind_tetr_save,iper_phi,k, m, n, j
    double precision                                :: perpinv,perpinv2, speed, r, z, phi, B, phi_elec_func
    double precision, dimension(:,:)                :: start_pos_pitch_mat
    type(optional_quantities_type)                  :: optional_quantities
!
    !If orbit_timestep is called for the first time without grid position
    if(.not.boole_initialized) then
!
        !Check coordinate domain (optionally perform modulo operation)
        call check_coordinate_domain(x)
!
        !Find tetrahedron index and face index for position x
        call find_tetra(x,vpar,vperp,ind_tetr,iface)
!
        !If particle doesn't lie inside any tetrahedron
        if(ind_tetr.eq.-1) then
            return
        endif
!
        r = x(1) - verts(1, tetra_grid(ind_tetr)%ind_knot(1))
        z = x(3) - verts(3,tetra_grid(ind_tetr)%ind_knot(1))
        phi = x(2) - verts(2,tetra_grid(ind_tetr)%ind_knot(1))
!
        if (boole_refined_sqrt_g) then
            weights(n+num_particles*j,1) = weights(n+num_particles*j,1)* &
                                    &  (sqrt_g(ind_tetr,1)+r*sqrt_g(ind_tetr,2)+z*sqrt_g(ind_tetr,3))/ &
                                    &  (sqrt_g(ind_tetr,4)+r*sqrt_g(ind_tetr,5)+z*sqrt_g(ind_tetr,6))
        endif
!       
        if (boole_boltzmann_energies) then
            !compare with equation 133 of master thesis of Jonatan Schatzlmayr (remaining parts have been added before)
            phi_elec_func = tetra_physics(ind_tetr)%Phi1 + sum(tetra_physics(ind_tetr)%gPhi*(/r,phi,z/))
            weights(n+num_particles*j,1) = weights(n+num_particles*j,1)*sqrt(start_pos_pitch_mat(5,n+num_particles*j)*ev2erg)* &
                        & exp(-(start_pos_pitch_mat(5,n+num_particles*j)*ev2erg+particle_charge*phi_elec_func)/(energy_eV*ev2erg))
        endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Linear density simulation !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
        if (boole_linear_density_simulation) then
            poloidal_flux(n+num_particles*j) = tetra_physics(ind_tetr)%Aphi1 + sum(tetra_physics(ind_tetr)%gAphi*(/r,phi,z/))
            weights(n+num_particles*j,1) = weights(n+num_particles*j,1)*(max_poloidal_flux*1.1-poloidal_flux(n))
        endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
        !compute J_perp for perpendicular pressure
        B = tetra_physics(ind_tetr)%bmod1+sum(tetra_physics(ind_tetr)%gb*(/r,phi,z/))
        J_perp(n+num_particles*j) = particle_mass*vmod**2*(1-start_pos_pitch_mat(4,n+num_particles*j)**2)*cm_over_e/(2*B)*(-1)
        !-1 because of negative gyrophase
!
        boole_initialized = .true.
    endif
!           
    !Exit the subroutine after initialization, if time step equals zero
    if(t_step.eq.0.d0) return
!
    !Squared perpendicular velocity
    vperp2 = vperp**2
!
    !Compute relative particle position
    z_save = x-tetra_physics(ind_tetr)%x1
!
    !Compute perpendicular invariant of particle
    perpinv=-0.5d0*vperp2/bmod_func(z_save,ind_tetr)
    perpinv2 = perpinv**2
!               
    !Initialize constants of motion in particle-private module
    select case(ipusher)
        case(1)
            call initialize_const_motion_rk(perpinv,perpinv2)
        case(2)
            call initialize_const_motion_poly(perpinv,perpinv2)
    end select        
!
!             !NOT FULLY IMPLEMENTED YET: Precompute quantities dependent on perpinv
!             call alloc_precomp_poly_perpinv(1,ntetr)
!             call initialize_boole_precomp_poly_perpinv()
!             call make_precomp_poly_perpinv(perpinv,perpinv2)
!
    !Integrate particle orbit for given time step
    t_remain = t_step
!
    !Logical for handling time integration
    boole_t_finished = .false.
!
    n_pushings = n_pushings-1 !set n_pushings to -1 because when entering the loop it wil go back to one without pushing
    !Loop for tetrahedron pushings until t_step is reached
    do
!
        n_pushings = n_pushings + 1
        !Domain Boundary
        if(ind_tetr.eq.-1) then
            !print *, 'WARNING: Particle lost.'
            if( present(t_remain_out)) then
                t_remain_out = t_remain
            endif
            exit
        endif
!                
        !Save the tetrahedron index for computation of vperp in the last step
        ind_tetr_save = ind_tetr
!
        !Save vpar for the computation of the parallel adiabatic invariant
        vpar_save = vpar
!  
        !t_remain (in) ... remaining time until t_step is finished
        !t_pass (out) ... time to pass the tetrahdron

        !Save x for the computation of the current
        x_save = x
!
        !Calculate trajectory
        select case(ipusher)
            case(1)
                call pusher_tetra_rk(ind_tetr,iface,x,vpar,z_save,t_remain,t_pass,boole_t_finished,iper_phi)
            case(2)
                call pusher_tetra_poly(poly_order,ind_tetr,iface,x,vpar,z_save,t_remain,&
                                                    & t_pass,boole_t_finished,iper_phi,optional_quantities)
        end select
!
        t_remain = t_remain - t_pass
        if (iper_phi.ne.0) counter_phi_0_mappings = counter_phi_0_mappings + iper_phi
!
        do m = 1,n_moments
        !print*, n_moments, moments_selector
            select case(moments_selector(m))
                case(1)
                    tetr_moments(m,ind_tetr_save,1) = tetr_moments(m,ind_tetr_save,1) + weights(n+num_particles*j,1)* &
                                                        & optional_quantities%t_hamiltonian!* &
                                                        !& (exp(2*(0,1)*x(2))+exp(3*(0,1)*x(2)))
                case(2)
                    tetr_moments(m,ind_tetr_save,1) = tetr_moments(m,ind_tetr_save,1) + weights(n+num_particles*j,1)* &
                                                                & optional_quantities%gyrophase*J_perp(n+num_particles*j)
                case(3)
                    tetr_moments(m,ind_tetr_save,1) = tetr_moments(m,ind_tetr_save,1) + weights(n+num_particles*j,1)* &
                                                                & optional_quantities%vpar_int
                case(4)
                    tetr_moments(m,ind_tetr_save,1) = tetr_moments(m,ind_tetr_save,1) + weights(n+num_particles*j,1)* &
                                                                & optional_quantities%vpar2_int
            end select
!
        enddo
!
        !Orbit stops within cell, because "flight"-time t_step has finished
        if(boole_t_finished) then
            if( present(t_remain_out)) then
                t_remain_out = t_remain
            endif
            exit
        endif
!
    enddo !Loop for tetrahedron pushings
!
    !Compute vperp from position
    vperp = vperp_func(z_save,perpinv,ind_tetr_save)
!            
!             !NOT FULLY IMPLEMENTED YET: Deallocate precomputed quantities dependent on perpinv
!             call alloc_precomp_poly_perpinv(2,ntetr)
!print*, 'number of pushings is', n_pushings
!         
end subroutine orbit_timestep_gorilla_boltzmann
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine calc_volume_integrals
!
    use constants, only: pi, ev2erg
    use tetra_grid_mod, only: ntetr, verts_rphiz, tetra_grid
    use tetra_grid_settings_mod, only: grid_size, n_field_periods
    use tetra_physics_mod, only: particle_mass,particle_charge, tetra_physics
!
    implicit none
!
    integer                                     :: i,j,k
    integer, dimension(2)                       :: triangle_indices
    double precision, dimension(3)              :: r_values, z_values, r_values_intermediate, z_values_intermediate, gradient
    double precision, dimension(2)              :: r, z, limits
    double precision                            :: z_star, alpha, beta, gamma, delta, epsilon, capital_gamma, capital_delta, &
                                                   a, a_dash, b, b_dash, c, c_dash, rmin, delta_r, delta_z, phi_0, eta
!
    allocate(tetra_indices_per_prism(n_prisms,3))    
    allocate(prism_volumes(n_prisms))
    allocate(refined_prism_volumes(n_prisms))
    allocate(r_integrand_constants(n_prisms,22)) !collects all constants used for integration in order to print them on a file
!
    refined_prism_volumes = 0
    r_integrand_constants = 0
!
    do i = 1,3
        tetra_indices_per_prism(:,i) = (/(i+3*k,k = 0,n_prisms-1)/)
    enddo
!
    !$OMP PARALLEL DEFAULT(NONE) &
    !$OMP& SHARED(n_prisms,verts_rphiz,tetra_grid,grid_size,tetra_indices_per_prism,prism_volumes, particle_charge,energy_ev, &
    !$OMP& refined_prism_volumes,sqrt_g,r_integrand_constants, elec_pot_vec,n_b,tetra_physics,boole_boltzmann_energies,density, &
    !$OMP& boole_refined_sqrt_g,n_field_periods) &
    !$OMP& PRIVATE(r_values,z_values,rmin,triangle_indices,r_values_intermediate,z_values_intermediate, delta_r, delta_z, eta, &
    !$OMP& r,z,gradient,z_star,alpha,beta,gamma,delta,epsilon,capital_gamma,capital_delta, limits,a,a_dash,b,b_dash,c,c_dash,phi_0)
    !$OMP DO
    do i = 1,n_prisms
        !calculate prism volumes using the basic approach (compare with chapter 4.2 of master thesis from Jonatan Schatzlmayr)
        r_values = verts_rphiz(1,tetra_grid(tetra_indices_per_prism(i,1))%ind_knot([1,2,3]))
        z_values = verts_rphiz(3,tetra_grid(tetra_indices_per_prism(i,1))%ind_knot([1,2,3]))
        rmin = minval(r_values)
        r_values = r_values - rmin
        z_values = z_values - z_values(minloc(r_values,1))
!
        triangle_indices(2) = maxloc(r_values,1)
        r_values_intermediate = r_values
        r_values_intermediate(triangle_indices(2)) = minval(r_values)
        if (sum(r_values_intermediate).eq.0) then
            z_values_intermediate = z_values
            z_values_intermediate(triangle_indices(2)) = minval(z_values)
            triangle_indices(1) = maxloc(abs(z_values_intermediate),1)
        else
            triangle_indices(1) = maxloc(r_values_intermediate,1)
        endif
!
        r = (/r_values(triangle_indices(1)),r_values(triangle_indices(2))/)
        z = (/z_values(triangle_indices(1)),z_values(triangle_indices(2))/)
!
        if (r(2).eq.0) then
            gradient(1) = 0
        else
            gradient(1) = z(2)/r(2)
        endif
        if (r(1).eq.0) then
            gradient(2) = 0
        else
            gradient(2) = z(1)/r(1)
        endif
        if ((r(2)-r(1)).eq.0) then
            gradient(3) = 0
        else
            gradient(3) = (z(2)-z(1))/(r(2)-r(1))
        endif
!
        r_integrand_constants(i,20:22) = (/minloc(r_values),triangle_indices(1),triangle_indices(2)/)
!
        z_star = z(1) - r(1)*gradient(3)
        alpha = abs(gradient(2)-gradient(1))
        beta = gradient(3) - gradient(1)
!
        r_integrand_constants(i,1:6) = (/r(1),r(2),alpha,rmin, z_star,beta/)
!
        prism_volumes(i) =  2*pi/(grid_size(2)*n_field_periods)*(alpha/3*r(1)**3+alpha*rmin/2*r(1)**2+ &
                            abs(z_star*rmin*(r(2)-r(1))+(z_star+beta*rmin)/2*(r(2)**2-r(1)**2)+beta/3*(r(2)**3-r(1)**3)))
!
        !calculate other volme integrals (compare with appendix B of master thesis of Jonatan Schatzlmayr)
        delta_r = verts_rphiz(1,tetra_grid(tetra_indices_per_prism(i,1))%ind_knot(1)) - &
                & verts_rphiz(1,tetra_grid(tetra_indices_per_prism(i,1))%ind_knot(int(r_integrand_constants(i,20))))
        delta_z = verts_rphiz(3,tetra_grid(tetra_indices_per_prism(i,1))%ind_knot(1)) - &
                & verts_rphiz(3,tetra_grid(tetra_indices_per_prism(i,1))%ind_knot(int(r_integrand_constants(i,20))))
!
        if (boole_refined_sqrt_g) then
            !Calculate prism volumes using the refined approach
            !(compare with appendix B (introductory pages + B1) of master thesis from Jonatan Schatzlmayr)
            a = sqrt_g(3*i-2,1) - sqrt_g(3*i-2,2)*delta_r - sqrt_g(3*i-2,3)*delta_z
            a_dash = sqrt_g(3*i-2,4) - sqrt_g(3*i-2,5)*delta_r - sqrt_g(3*i-2,6)*delta_z
            b = sqrt_g(3*i-2,2)
            b_dash = sqrt_g(3*i-2,5)
            c = sqrt_g(3*i-2,3)
            c_dash = sqrt_g(3*i-2,6)
    !
            !calculate the contribution from 0 to r(1) to the prism volumes
            limits = (/dble(0),r(1)/)
    !
            alpha = c/c_dash*(gradient(2)-gradient(1))
            beta = b_dash+c_dash*gradient(2)
            gamma = b_dash+c_dash*gradient(1)
            delta = (c_dash*a-c*a_dash)/c_dash**2
            epsilon = (c_dash*b-c*b_dash)/c_dash**2
    !
            r_integrand_constants(i,7:12) = (/alpha,beta,gamma,delta,epsilon,a_dash/)
    !
            refined_prism_volumes(i) = r(1)*epsilon*a_dash/(2*gamma*beta)*(gamma-beta) + &
                    & r(1)**2*alpha/2 + &
                    & log((a_dash+gamma*r(1))/a_dash)*(delta*a_dash/(gamma*beta)*(gamma-beta)+epsilon*a_dash**2/(2*gamma**2)) + &
                    & log((a_dash+beta*r(1))/(a_dash+gamma*r(1)))*(delta/beta*(a_dash+beta*r(1))+epsilon/2*r(1)**2) - &
                    & log((a_dash+beta*r(1))/a_dash)*epsilon*a_dash**2/(2*beta**2)
    !
            !calculate the contribution from r(1) to r(2) to the prism volumes
            limits = (/r(1),r(2)/)
    !
            alpha = c/c_dash*(gradient(3)-gradient(1))
            beta = b_dash+c_dash*gradient(3)
            gamma = b_dash+c_dash*gradient(1)
            delta = (c_dash*a-c*a_dash)/c_dash**2
            epsilon = (c_dash*b-c*b_dash)/c_dash**2
            capital_gamma = c/c_dash*z_star
            capital_delta = a_dash + c_dash*z_star
    !
            r_integrand_constants(i,13:19) = (/alpha,beta,gamma,delta,epsilon,capital_gamma,capital_delta/)
    !
            refined_prism_volumes(i) = refined_prism_volumes(i) + &
                    & (r(2)-r(1))*(capital_gamma+epsilon/(2*gamma*beta)*(gamma*capital_delta-beta*a_dash)) + &
                    & (r(2)**2-r(1)**2)*alpha/2 + &
                    & log((a_dash+gamma*r(2))/(a_dash+gamma*r(1)))*(delta/(gamma*beta)*(gamma*capital_delta-a_dash*beta) + &
                                                                                            & epsilon*a_dash**2/(2*gamma**2)) + &
                    & log((capital_delta+beta*r(2))/(a_dash+gamma*r(2))*(a_dash+gamma*r(1))/(capital_delta+beta*r(1))) * &
                                                                    & (delta/beta*(capital_delta+beta*r(2))+epsilon/2*r(2)**2) + &
                    & log((capital_delta+beta*r(1))/(a_dash+gamma*r(1)))*(delta*(r(2)-r(1))+epsilon/2*(r(2)**2-r(1)**2)) - &
                    & log((capital_delta+beta*r(2))/(capital_delta+beta*r(1)))*epsilon*capital_delta**2/(2*beta**2)
        endif
!
        if (boole_boltzmann_energies) then
!
            !calculate electric potential using the refined approach
            !(compare with appendix B (introductory pages + B2) of master thesis from Jonatan Schatzlmayr)
!
            a = tetra_physics(3*i-2)%gPhi(1)
            b = tetra_physics(3*i-2)%gPhi(3)
            phi_0 = tetra_physics(3*i-2)%Phi1 - a*delta_r - b*delta_z
!
            !calculate the contribution from 0 to r(1) to the electric potential
            alpha = phi_0*(gradient(2)-gradient(1))*rmin
            beta = (a*rmin+phi_0)*(gradient(2)-gradient(1))+b/2*(gradient(2)**2-gradient(1)**2)*rmin
            gamma = a*(gradient(2)-gradient(1))+b/2*(gradient(2)**2-gradient(1)**2)

            elec_pot_vec(i) = elec_pot_vec(i) + alpha/2*r(1)**2+beta/3*r(1)**3+gamma/4*r(1)**4
!
            !calculate the contribution from r(1) to r(2) to the electric potential
            alpha  = phi_0*z_star*rmin+b/2*z_star**2*rmin
            beta = a*z_star*rmin+phi_0*(gradient(3)-gradient(1))*rmin+phi_0*z_star+b/2*z_star**2+b*z_star*gradient(3)*rmin
            gamma = phi_0*(gradient(3)-gradient(1))+a*z_star+a*(gradient(3)-gradient(1))*rmin + &
                    & b/2*(gradient(3)**2-gradient(1)**2)*rmin+b*z_star*gradient(3)
            delta = a*(gradient(3)-gradient(1))+b/2*(gradient(3)**2-gradient(1)**2)

            elec_pot_vec(i) = elec_pot_vec(i) + alpha*(r(2)-r(1))+beta/2*(r(2)**2-r(1)**2)+&
                            & gamma/3*(r(2)**3-r(1)**3)+delta/4*(r(2)**4-r(1)**4)
!
            !calculate Boltzmann densiy using the refined approach
            !(compare with appendix B (introductory pages + B3) of master thesis from Jonatan Schatzlmayr)
!
            !calculate contribution from 0 to r(1) to the boltzmann density    
            alpha = density*exp(-particle_charge*phi_0/(energy_ev*ev2erg))*rmin*(gradient(2)-gradient(1))
!
            beta = density*exp(-particle_charge*phi_0/(energy_ev*ev2erg))*((gradient(2)-gradient(1))* &
                 & (1-rmin*particle_charge*a/(energy_ev*ev2erg))- &
                 & rmin*particle_charge*b/(2*energy_ev*ev2erg)*(gradient(2)**2-gradient(1)**2))
!
            gamma = density*exp(-particle_charge*phi_0/(energy_ev*ev2erg))* &
                  & (-particle_charge*a/(energy_ev*ev2erg)*(gradient(2)-gradient(1))- &
                  & particle_charge*b/(2*energy_ev*ev2erg)*(gradient(2)**2-gradient(1)**2))
!
            n_b(i) = n_b(i) + alpha/2*r(1)**2+beta/3*r(1)**3+gamma/4*r(1)**4
!
            !calculate contribution from r(1) to r(2) to the boltzmann density
            alpha = density*exp(-particle_charge*phi_0/(energy_ev*ev2erg))*rmin*&
                  & (z_star-particle_charge*b/(2*energy_ev*ev2erg)*z_star**2)
!
            beta = density*exp(-particle_charge*phi_0/(energy_ev*ev2erg))*(z_star- & 
                 & rmin*particle_charge*a/(energy_ev*ev2erg)*z_star+rmin*(gradient(3)-gradient(1))- &
                 & particle_charge*b/(2*energy_ev*ev2erg)*z_star**2-rmin*particle_charge*b/(energy_ev*ev2erg)*z_star*gradient(3))
!
            gamma = density*exp(-particle_charge*phi_0/(energy_ev*ev2erg))*((-rmin*particle_charge*a/(energy_eV*ev2erg)+1)* &
                  & (gradient(3)-gradient(1))-particle_charge*a/(energy_ev*ev2erg)*z_star- &
                  & rmin*particle_charge*b/(2*energy_ev*ev2erg)*(gradient(3)**2-gradient(1)**2)- &
                  & particle_charge*b/(energy_eV*ev2erg)*z_star*gradient(3))
!
            delta = density*exp(-particle_charge*phi_0/(energy_ev*ev2erg))*(-particle_charge*a/(energy_eV*ev2erg)* &
                  & (gradient(3)-gradient(1))-particle_charge*b/(2*energy_eV*ev2erg)*(gradient(3)**2-gradient(1)**2))
!
            n_b(i) = n_b(i) + alpha*(r(2)-r(1)) + beta/2*(r(2)**2-r(1)**2) + gamma/3*(r(2)**3-r(1)**3) + delta/4*(r(2)**4-r(1)**4)
!
        endif
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
!
    refined_prism_volumes = abs(refined_prism_volumes)*2*pi/(grid_size(2)*n_field_periods)
    elec_pot_vec = abs(elec_pot_vec)*2*pi/(grid_size(2)*n_field_periods*prism_volumes)
    n_b = abs(n_b)*2*pi/(grid_size(2)**n_field_periods*prism_volumes)
!
end subroutine calc_volume_integrals
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine fourier_transform_moments
!
    use constants, only: pi
    use tetra_grid_settings_mod, only: grid_size
!
    implicit none 
!
    integer                                                        :: n,m,j,k,p,q,l
    complex                                                        :: i
    complex, dimension(:,:,:,:), allocatable              :: prism_moments_ordered_for_ft
!
print*, 'fourier transform started'
!
    n_triangles = n_prisms/grid_size(ind_b)
    allocate(prism_moments_ordered_for_ft(n_moments,n_triangles,grid_size(ind_b),n_species))
    do q = 1,grid_size(ind_b)
        prism_moments_ordered_for_ft(:,:,q,:) = prism_moments(:,n_triangles*(q-1)+1:n_triangles*q,:)
    enddo
!
    if (n_fourier_modes.gt.grid_size(ind_b)) then
        print*, 'n_fourier_modes was chosen to be bigger than n_phi, it is therefore reduced to n_phi'
        n_fourier_modes = grid_size(ind_b)
    endif
!
    allocate(moments_in_frequency_space(n_moments,n_triangles,n_fourier_modes,n_species))
    moments_in_frequency_space = 0
    i = (0,1)
!
    do p = 1,n_triangles
        do n = 1,n_moments
            do m = 1,n_species
                do k = 0,n_fourier_modes-1
                    do j = 0,grid_size(ind_b)-1
                        moments_in_frequency_space(n,p,k+1,m) = moments_in_frequency_space(n,p,k+1,m) + 1/dble(grid_size(ind_b))* &
                                                        & exp(-2*pi*i*j*k/grid_size(ind_b))*prism_moments_ordered_for_ft(n,p,j+1,m)
                    enddo
                enddo
            enddo
        enddo
    enddo
!
    open(63, file = 'prism_moments_for_fourier.dat')
    do l = 1,n_triangles
        do m = 1,grid_size(ind_b)-1
            write(63,'(2ES20.10E4)',advance="no") real(prism_moments_ordered_for_ft(1,l,m,1)), &
                                                  aimag(prism_moments_ordered_for_ft(1,l,m,1))     
        enddo
            write(63,'(2ES20.10E4)') real(prism_moments_ordered_for_ft(1,l,grid_size(ind_b),1)), &
                                     aimag(prism_moments_ordered_for_ft(1,l,grid_size(ind_b),1))
    enddo
    close(63)
!
    deallocate(prism_moments_ordered_for_ft)
!
print*, 'fourier transform finished'
!
end subroutine fourier_transform_moments
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine calc_square_root_g
!
    use tetra_physics_mod, only: tetra_physics, hamiltonian_time
    use tetra_grid_mod, only: ntetr, tetra_grid, verts_rphiz
!
    implicit none
!
    integer            :: ind_tetr
!
    allocate(sqrt_g(ntetr,7))
    !compare first 6 entries with chapter 4.5 of master thesis from Jonatan Schatzlmayr, entry 7 is the radius (old metric determinant)
!
    do ind_tetr = 1, ntetr
        sqrt_g(ind_tetr,1) = hamiltonian_time(ind_tetr)%h1_in_curlA
!
        sqrt_g(ind_tetr,2) = tetra_physics(ind_tetr)%gh1(1)*tetra_physics(ind_tetr)%curlA(1) + &
                           & tetra_physics(ind_tetr)%gh2(1)*tetra_physics(ind_tetr)%curlA(2) + &
                           & tetra_physics(ind_tetr)%gh3(1)*tetra_physics(ind_tetr)%curlA(3)
!
        sqrt_g(ind_tetr,3) = tetra_physics(ind_tetr)%gh1(3)*tetra_physics(ind_tetr)%curlA(1) + &
                           & tetra_physics(ind_tetr)%gh2(3)*tetra_physics(ind_tetr)%curlA(2) + &
                           & tetra_physics(ind_tetr)%gh3(3)*tetra_physics(ind_tetr)%curlA(3)
!
        sqrt_g(ind_tetr,4) = tetra_physics(ind_tetr)%bmod1
!
        sqrt_g(ind_tetr,5) = tetra_physics(ind_tetr)%gB(1)
!
        sqrt_g(ind_tetr,6) = tetra_physics(ind_tetr)%gB(3)
!
        sqrt_g(ind_tetr,7) = verts_rphiz(1,tetra_grid(ind_tetr)%ind_knot(1))
    enddo
!
end subroutine calc_square_root_g
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!!!!!!!!!!!!!!!!! The following two subroutines are currently not made use of !!!!!!!!!!!!
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine orbit_timestep_can_boltzmann(z,dtau,dtaumin,ierr,tau_out)
    !
          use odeint_mod, only: odeint_allroutines
          use runge_kutta_mod, only: runge_kutta_allroutines
          use gorilla_settings_mod, only: rel_err_ode45
    !
          implicit none
    !
          integer, parameter          :: ndim=5, nstepmax=10000000
          double precision            :: relerr
    !
          integer :: ierr,j
          double precision :: dtau,dtaumin,phi,tau1,tau2,tau_out
    !
          double precision, dimension(2)    :: y
          double precision, dimension(ndim) :: z
    !
          external velo_can
    !
          !Use relative error definition from input file
          relerr = rel_err_ode45
    !
          if(dtaumin*nstepmax.le.abs(dtau)) then
            ierr=2
            print *,'orbit_timestep: number of steps exceeds nstepmax'
            return
          endif
    !
          ierr=0
          y(1)=z(1)
          y(2)=z(2)
          phi=z(3)
    !
          call chamb_can(y,phi,ierr)
    !
          tau_out = 0.d0
    !
          if(ierr.ne.0) return
          tau1=0.d0
          tau2=sign(dtaumin,dtau)      
    !
          do while(abs(tau2).lt.abs(dtau))
    !
    !print *, '(1) before RK4', z(1)
            call odeint_allroutines(z,ndim,tau1,tau2,relerr,velo_can)
    
    !        call runge_kutta_allroutines(4,z,ndim,tau1,tau2,velo_can)
    !print *, '(1) after RK4', z(1)
    !
            y(1)=z(1)
            y(2)=z(2)
            phi=z(3)
    !
            call chamb_can(y,phi,ierr)
    !
            tau_out = tau2
    !
            if(ierr.ne.0) return
            tau1=tau2
            tau2=tau2+sign(dtaumin,dtau)
          enddo
    !
          tau2=dtau
    !
    !print *, '(2) before RK4', z
          call odeint_allroutines(z,ndim,tau1,tau2,relerr,velo_can)
    
    !      call runge_kutta_allroutines(4,z,ndim,tau1,tau2,velo_can)
    !print *, '(2) after RK4', z
    !
          y(1)=z(1)
          y(2)=z(2)
          phi=z(3)
    !
          call chamb_can(y,phi,ierr)
    !
          tau_out = tau2
    !
          if(ierr.ne.0) return
    !
          end subroutine orbit_timestep_can_boltzmann
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine read_in_starting_conditions(start_pos_pitch_mat, n_start_pos)
    double precision, dimension(:,:), allocatable, intent(out) :: start_pos_pitch_mat
    integer, intent(out)                                       :: n_start_pos
    integer                                                    :: file_id_read_start, i_os, i
! 
    file_id_read_start = 100
    open(unit=file_id_read_start, file = filename_starting_conditions, iostat=i_os, status='old')
!
        !Error, if file does not exist.
        if ( i_os /= 0 ) then
                !Symmetry flux coordinates
                    print *, "Error opening file with starting positions and starting pitch parameter: ", &
                    & filename_starting_conditions
            stop
        endif
!
        !Count number of lines
            n_start_pos = 0
!
        do
            read(file_id_read_start, '(A)', iostat=i_os)
            if (i_os /= 0) exit
            n_start_pos = n_start_pos + 1
        end do
!
        print*, "File with starting positions and pitch parameter contains ", n_start_pos, "starting values."
!
        allocate(start_pos_pitch_mat(n_start_pos,4))
!
        rewind(file_id_read_start)
!
        do i = 1, n_start_pos
            read(file_id_read_start,*) start_pos_pitch_mat(i,:)
        end do
!
    close(file_id_read_start)

end subroutine read_in_starting_conditions
!
end module boltzmann_mod