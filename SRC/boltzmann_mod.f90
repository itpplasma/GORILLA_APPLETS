module boltzmann_mod
    implicit none
!
    private
!
    double precision :: time_step,energy_max_eV
    integer, dimension(:,:), allocatable :: tetra_indices_per_prism
    double precision, dimension(:), allocatable :: prism_volumes
    double precision, dimension(:,:,:), allocatable :: tetr_moments, prism_moments
    complex, dimension(:,:,:,:), allocatable :: moments_in_frequency_space
    integer :: i_integrator_type, seed_option, n_moments, n_species, n_prisms, n_t_measurements, num_particles, ind_a, ind_b, ind_c
    double precision :: n_particles
    integer, dimension(4) :: moments_selector
    logical :: boole_random_precalc
    character(1024) :: filename_dwell_times, filename_starting_conditions, filename_vertex_coordinates, &
    & filename_vertex_indices
    double precision, dimension(:), allocatable :: dwell_times
    double precision, dimension(:,:), allocatable :: currents !, time_resolved_energies
!
    !Namelist for boltzmann input
    NAMELIST /boltzmann_nml/ time_step,energy_max_eV,n_particles,i_integrator_type, &
    & seed_option,boole_random_precalc,filename_dwell_times,filename_starting_conditions, &
    & filename_vertex_coordinates, filename_vertex_indices!, n_t_measurements
!
    !boole_random_precalc,
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
    use constants, only: pi
    use tetra_grid_mod, only: verts_rphiz, verts_sthetaphi, ntetr
    use find_tetra_mod, only: find_tetra
    use tetra_physics_mod, only: coord_system
!
    implicit none
    double precision, intent(in)                                   :: vmod
    double precision, dimension(:,:), allocatable, intent(out)     :: start_pos_pitch_mat !dimension(4,n_particle)
    double precision                                               :: rand_scalar, vpar, vperp
    double precision                                               :: amin, amax, cmin, cmax
    double precision, dimension(:), allocatable                    :: rand_vector
    integer, dimension(:), allocatable                             :: starting_tetrahedra, starting_prisms
    double precision, dimension(:,:), allocatable                  :: rand_matrix
    integer                                                        :: i
    logical                                                        :: inside
    double precision, dimension(3)                                 :: x
    integer                                                        :: ind_tetr_out,iface, counter
    double precision, dimension(:,:), allocatable                  :: verts
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
    allocate(start_pos_pitch_mat(4,num_particles))
    allocate(rand_vector(num_particles))
    allocate(rand_matrix(4,num_particles))
    allocate(starting_tetrahedra(ntetr))
    allocate(starting_prisms(n_prisms))
!
    inside = .false.
    start_pos_pitch_mat = 0
    counter = 0
    starting_tetrahedra = 0
    starting_prisms = 0
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
    ! PRINT*, 'amin, amax, cmin and cmax are:', amin, amax, cmin, cmax
    ! PRINT*, 'numel(verts) = ', size(verts(1,:))
!
    ! open(55, file = 'vertices.dat')
    ! write(55,'(3ES15.3E4)') verts
!
    call RANDOM_NUMBER(rand_matrix)
    start_pos_pitch_mat(ind_a,:) = amin + (amax - amin)*rand_matrix(ind_a,:) !r in cylindrical, s in flux coordinates
    start_pos_pitch_mat(ind_b,:) = 2*pi*rand_matrix(ind_b,:) !phi in cylindrical and flux coordinates
    start_pos_pitch_mat(ind_c,:) = cmin + (cmax - cmin)*rand_matrix(ind_c,:) !z in cylindrical, theta in flux coordinates
    start_pos_pitch_mat(4,:) = 2*rand_matrix(4,:)-1 !pitch parameter
!
!
!
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  

    j = 0
    !$OMP PARALLEL DEFAULT(NONE) &
    !$OMP& SHARED(num_particles,start_pos_pitch_mat,amin,amax,cmin,cmax,vmod,ind_a,ind_b,ind_c,j) &
    !$OMP& PRIVATE(i,inside,rand_scalar,vpar,vperp,x,ind_tetr_out,iface) &
    !$OMP& REDUCTION(+:starting_tetrahedra)
    !$OMP DO
    do i = 1,num_particles
!         ! do while(inside.eqv..false.)
!             ! !$omp critical
!             ! counter = counter + 1
!             ! !$omp end critical
!             call RANDOM_NUMBER(rand_scalar)
!             start_pos_pitch_mat(ind_a,i) = amin + (amax - amin)*rand_scalar !r in cylindrical, s in flux coordinates
!             call RANDOM_NUMBER(rand_scalar)
!             start_pos_pitch_mat(ind_b,i) = 2*pi*rand_scalar !phi in cylindrical and flux coordinates
!             call RANDOM_NUMBER(rand_scalar)
!             start_pos_pitch_mat(ind_c,i) = cmin + (cmax - cmin)*rand_scalar !z in cylindrical, theta in flux coordinates
              x = start_pos_pitch_mat(1:3,i)
!             call RANDOM_NUMBER(rand_scalar)
!             start_pos_pitch_mat(4,i) = 2*rand_scalar-1 !pitch parameter
              vpar = start_pos_pitch_mat(4,i)*vmod
              vperp = sqrt(vmod**2-vpar**2)
              call find_tetra(x,vpar,vperp,ind_tetr_out,iface)
              if (ind_tetr_out.ne.-1) then
                starting_tetrahedra(ind_tetr_out) = starting_tetrahedra(ind_tetr_out) + 1!x(1)
              endif
!         !     if (ind_tetr_out.ne.-1) inside = .true.
!         ! enddo
!         ! inside = .false.
             !$omp critical
             j = j+1
             !$omp end critical
if (num_particles.ge.100) then
    if (modulo(j,int(num_particles/100)).eq.0) then
        print *, nint(dble(j)/dble(num_particles)*100), ' percent done'
    endif
endif
!         ! Print*, 'next'
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
!
! ! PRINT*, 'counter = ', counter
!
    ! call RANDOM_NUMBER(rand_vector)
    ! start_pos_pitch_mat(5,:) = energy_max_eV!*rand_vector !energy
!
    do i = 1,3
        starting_prisms = starting_prisms + starting_tetrahedra(tetra_indices_per_prism(:,i))
    enddo
    open(74, file = 'starting_positions.dat')
    write(74,'(2ES20.10E4)') start_pos_pitch_mat((/1,3/),:)
    close(74)
    open(75, file = 'starting_prisms.dat')
    do i=1, n_prisms
        write(75, *) starting_prisms(i)
    end do
    close(75)
!
end subroutine calc_starting_conditions
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine calc_boltzmann
!
    use orbit_timestep_gorilla_mod, only: initialize_gorilla !,orbit_timestep_gorilla
    use constants, only: ev2erg, pi
    use tetra_physics_mod, only: particle_mass,particle_charge,cm_over_e,mag_axis_R0, coord_system
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
    double precision, dimension(:), allocatable :: densities, single_particle_dwell_times!, measuring_times
    double precision, dimension(:,:), allocatable :: fluid_velocities, start_pos_pitch_mat
    double precision :: vmod,pitchpar,vpar,vperp,t_remain,t_confined,tau_out_can
    integer :: kpart,i,n,l,m,k,ind_tetr,iface,n_lost_particles,ierr
    integer :: n_start, n_end, i_part, n_fourier_modes
    double precision, dimension(3) :: x_rand_beg,x
    logical :: boole_initialized,boole_particle_lost
    double precision :: dtau, dphi,dtaumin
    double precision, dimension(5) :: z
    !Character(LEN=50) :: format_time_resolved_energies
    Character(LEN=50) :: format_moments
!
    open(35, file = 'outliers.dat')
    close(35,status='delete')
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
    n_moments = 0
    moments_selector = 0
    do i = 1,size(boole_array_optional_quantities)
        if (boole_array_optional_quantities(i).eqv..true.) then
            n_moments = n_moments + 1
            moments_selector(n_moments) = i
        endif
    enddo
!
    !measuring_times = (/(i,i=n_t_measurements-1,0,-1)/)
    !measuring_times = measuring_times/(n_t_measurements-1)*time_step
!
    !Initialize GORILLA
    call initialize_gorilla()
!
    n_prisms = ntetr/3
    allocate(dwell_times(1:ntetr))
    allocate(single_particle_dwell_times(1:ntetr))
    allocate(densities(1:ntetr))
    allocate(currents(1:3,1:ntetr))
    allocate(fluid_velocities(1:3,1:ntetr))
    !allocate(time_resolved_energies(1:n_t_measurements,1:num_particles))
    allocate(tetr_moments(n_moments,ntetr,n_species))
    allocate(prism_moments(n_moments,n_prisms,n_species))
! 
    !time_resolved_energies = 0
    dwell_times = 0
    currents = 0
    tetr_moments = 0
    prism_moments = 0
!
    !Load fluxtube volume for a starting position (File to be chosen in gorilla_applets.inp)
    !call load_flux_tube_volume()
!
    call calc_prism_volumes
!
    !Compute velocity module from kinetic energy dependent on particle species
    vmod=sqrt(2.d0*energy_max_eV*ev2erg/particle_mass)
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
                kpart = 0
                n_lost_particles = 0
    !
print*, 'now'
                !$OMP PARALLEL DEFAULT(NONE) &
                !$OMP& SHARED(num_particles,pos_fluxtv_mat,kpart,vmod,time_step,i_integrator_type, &
                !$OMP& dtau,dtaumin,boole_random_precalc,n_start,n_end, dwell_times, &
                !$OMP& currents, start_pos_pitch_mat, particle_mass) &
                !$OMP& PRIVATE(n,boole_particle_lost,x_rand_beg,x,pitchpar,vpar,vperp,boole_initialized, &
                !$OMP& ind_tetr,iface,t_remain,t_confined,z,ierr,tau_out_can, single_particle_dwell_times) &
                !$OMP& REDUCTION(+:n_lost_particles,tetr_moments)
!print*, 'get number of threads', omp_get_num_threads()
                !$OMP DO
    !
                !Loop over particles
                do n = n_start,n_end !1,num_particles
    !
                    !$omp critical
                    !Counter for particles
                    kpart = kpart+1
                    !boole_particle_lost = .false.
if (n_end.gt.10) then
    if (modulo(kpart,int(n_end/10)).eq.0) then
        print *, kpart, ' / ', num_particles, 'particle: ', n, 'thread: ' !, omp_get_thread_num()
    endif
else 
    print *, kpart, ' / ', num_particles, 'particle: ', n, 'thread: ' !, omp_get_thread_num()
endif
                    !$omp end critical
                    

                    !You need x_rand_beg(1,3), pitchpar(1) (between -1 and 1), energy is already given
                    x_rand_beg = start_pos_pitch_mat(1:3,n)
                    pitchpar = start_pos_pitch_mat(4,n)
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
                        case(2)
                            z(5) = pitchpar
                    end select

    !
                    !time_resolved_energies(1,n) = particle_mass*vmod**2/2
                    single_particle_dwell_times = 0
    !
                    !Orbit integration
                    select case(i_integrator_type)
    !
                        case(0)
                            !No orbit computation
    !
                        case(1)
                            boole_initialized = .false.
!print*, 'now'
                            call orbit_timestep_gorilla_boltzmann(x,vpar,vperp,time_step,boole_initialized, &
                            & ind_tetr,iface,t_remain,n,single_particle_dwell_times)
    !
                            !Confinement time of particles
                            t_confined = time_step - t_remain

                            !dwell times of particles

    !
                            !Lost particle handling
                            if(ind_tetr.eq.-1) then
                                !!$omp critical
                                    n_lost_particles = n_lost_particles + 1
                                    !boole_particle_lost = .true.
                                !!$omp end critical
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
                                !!$omp critical
                                    n_lost_particles = n_lost_particles + 1
                                    !boole_particle_lost = .true.
                                !!$omp end critical
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
                    !$omp critical
                    dwell_times = dwell_times + single_particle_dwell_times
                    !$omp end critical
    !
                enddo !n
                !$OMP END DO
                !$OMP END PARALLEL
print*, 'finished'
    !
    !            close(99)
    !            
                densities = dwell_times / time_step !compute densities
                do l = 1,ntetr !compute fluid velocity
                        if (densities(l) .EQ. 0) then
                            fluid_velocities(:,l) = 0
                        else
                            do m = 1,3
                                fluid_velocities(m,l) = currents(m,l)/(densities(l)*particle_charge)
                            enddo
                        endif
                enddo
!
                prism_moments = tetr_moments(:,tetra_indices_per_prism(:,1),:) + &
                              & tetr_moments(:,tetra_indices_per_prism(:,2),:) + &
                              & tetr_moments(:,tetra_indices_per_prism(:,3),:)
!
                ! do n = 1,n_moments
                !     do m = 1,n_species
                !         prism_moments(n,:,m) = prism_moments(n,:,m)/prism_volumes
                !     enddo
                ! enddo
!
    ! print *, 'Number of lost particles',n_lost_particles
    ! open(99,file='confined_fraction.dat')
    ! write(99,*) 1.d0-dble(n_lost_particles)/dble(num_particles)
!   
    open(50, file = filename_dwell_times)
    write(50,'(ES20.10E4)') dwell_times
    close(50)
    open(51, file = 'densities.dat')
    write(51,'(ES20.10E4)') densities
    close(51)
    open(52, file = 'currents.dat')
    write(52,'(3ES20.10E4)') currents
    close(52)
    open(53, file = 'fluid_velocities.dat')
    write(53,'(3ES20.10E4)') fluid_velocities
    close(53)
    ! if (num_particles.gt.0) then
    !     write(format_time_resolved_energies, *) '(',num_particles,'ES15.3E4)'
    !     open(54, file = 'time_resolved_energies.dat')
    !     write(54,format_time_resolved_energies) time_resolved_energies
    !     close(54)
    ! endif
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
                ![R,phi,Z]: Write vertex coordinates to File
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
        write(format_moments, *) '(',n_moments,'ES20.10E4)'
        write(59,format_moments) prism_moments(:,:,1)
        close(59)
        open(60, file = 'tetr_moments.dat')
        write(60,format_moments) tetr_moments(:,:,1)
        close(60)
    endif
    open(61, file = 'tetra_indices_per_prism.dat')
    do i = 1,n_prisms
        write(61,*) tetra_indices_per_prism(i,:)
    enddo
    close(61)
!
    call fourier_transform_moments(n_fourier_modes)
!
    deallocate(dwell_times, currents, fluid_velocities, densities, start_pos_pitch_mat, &
    & tetr_moments, prism_moments)
!
! PRINT*, 'particle mass = ', particle_mass
! PRINT*, 'large radius = ', mag_axis_R0
! PRINT*, 'parallel velocity = ', vpar
! PRINT*, 'absolute value of velocity = ', vmod
! PRINT*, 'perpendicular velocity = ', vperp
! PRINT*, 'pitch par =', pitchpar
! PRINT*, 'particle charge = ', particle_charge
!
end subroutine calc_boltzmann
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine orbit_timestep_gorilla_boltzmann(x,vpar,vperp,t_step,boole_initialized,ind_tetr,iface, t_remain_out, n, &
                                            & single_particle_dwell_times)
!
                use pusher_tetra_rk_mod, only: pusher_tetra_rk,initialize_const_motion_rk
                use pusher_tetra_poly_mod, only: pusher_tetra_poly,initialize_const_motion_poly
                use tetra_physics_poly_precomp_mod , only: make_precomp_poly_perpinv, initialize_boole_precomp_poly_perpinv, &
                    & alloc_precomp_poly_perpinv
                use tetra_physics_mod, only: tetra_physics,particle_charge,particle_mass
                use gorilla_settings_mod, only: ipusher, poly_order, optional_quantities_type, boole_array_optional_quantities
                use orbit_timestep_gorilla_mod, only:  check_coordinate_domain
                use supporting_functions_mod, only: bmod_func, vperp_func
                use find_tetra_mod, only: find_tetra
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
                double precision                                :: vperp2,t_remain,t_pass,vpar_save
                logical                                         :: boole_t_finished
                integer                                         :: ind_tetr_save,iper,k, i,m, n_pushings,n
                double precision                                :: perpinv,perpinv2, speed
                !double precision, dimension(:), intent(inout)   :: measuring_times
                type(optional_quantities_type)                  :: optional_quantities
                double precision, dimension(:), intent(inout)   :: single_particle_dwell_times
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
!             !NOT FULLY IMPLEMENTED YET: Precompute quatities dependent on perpinv
!             call alloc_precomp_poly_perpinv(1,ntetr)
!             call initialize_boole_precomp_poly_perpinv()
!             call make_precomp_poly_perpinv(perpinv,perpinv2)
!
                !Integrate particle orbit for given time step
                t_remain = t_step
!
                !Logical for handling time integration
                boole_t_finished = .false.

                !i checks when to do time measurements
                i = 2
!
                n_pushings = -1
                !Loop for tetrahedron pushings until t_step is reached
                do
!
                    n_pushings = n_pushings + 1
                    !Domain Boundary
                    if(ind_tetr.eq.-1) then
                        print *, 'WARNING: Particle lost.'
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
                            call pusher_tetra_rk(ind_tetr,iface,x,vpar,z_save,t_remain,t_pass,boole_t_finished,iper)
                        case(2)
                            call pusher_tetra_poly(poly_order,ind_tetr,iface,x,vpar,z_save,t_remain,&
                                                                & t_pass,boole_t_finished,iper,optional_quantities)
                    end select
!
                    t_remain = t_remain - t_pass
!
                    ! if (t_remain .LE. measuring_times(i)) then
                    !     speed = sqrt((x(1)-x_save(1))**2+(x(2)-x_save(2))**2+(x(3)-x_save(3))**2)/t_pass                        
                    !     time_resolved_energies(i,n) = particle_mass*speed**2/2
                    !     i = i+1
                    ! endif
!
                    !!$omp critical
                    single_particle_dwell_times(ind_tetr_save) = single_particle_dwell_times(ind_tetr_save) + &
                                                                & optional_quantities%t_hamiltonian
                    !dwell_times(ind_tetr_save) = dwell_times(ind_tetr_save) + optional_quantities%t_hamiltonian
                    currents(:,ind_tetr_save) = particle_charge*(x-x_save)/t_step ! /t_pass vor the speed but *t_pass to calculate the time averge together with /t_step
!
                    do m = 1,n_moments
                    !print*, n_moments, moments_selector
                        select case(moments_selector(m))
                            case(1)
                                tetr_moments(m,ind_tetr_save,1) = tetr_moments(m,ind_tetr_save,1) + &
                                                                            & optional_quantities%t_hamiltonian
if ((optional_quantities%t_hamiltonian.gt.1.4d-30).or.(optional_quantities%t_hamiltonian.lt.0.8d-30)) then
    print*, 'hamiltonian time and r are : ', optional_quantities%t_hamiltonian, x(1)
endif

                            case(2)
                                tetr_moments(m,ind_tetr_save,1) = tetr_moments(m,ind_tetr_save,1) + &
                                                                            & optional_quantities%gyrophase
                            case(3)
                                tetr_moments(m,ind_tetr_save,1) = tetr_moments(m,ind_tetr_save,1) + &
                                                                            & optional_quantities%vpar_int
                            case(4)
                                tetr_moments(m,ind_tetr_save,1) = tetr_moments(m,ind_tetr_save,1) + &
                                                                            & optional_quantities%vpar2_int
                        end select
!print*, 'hamiltonian time is ', optional_quantities%t_hamiltonian
                    enddo
                    !!$omp end critical
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
subroutine calc_prism_volumes
!
    use constants, only: pi
    use tetra_grid_mod, only: ntetr, verts_rphiz, tetra_grid
    use tetra_grid_settings_mod, only: grid_size
!
    implicit none
!
    integer                                     :: i,k
    integer, dimension(2)                       :: triangle_indices
    double precision, dimension(3)              :: r_values, z_values, r_values_intermediate, z_values_intermediate, gradient
    double precision, dimension(2)              :: r, z
    double precision                            :: z_star, alpha, beta, rmin
!
    print*, 'start prism volume computation'
    allocate(tetra_indices_per_prism(n_prisms,3))    
    allocate(prism_volumes(n_prisms))
!
    do i = 1,3
        tetra_indices_per_prism(:,i) = (/(i+3*k,k = 0,n_prisms-1)/)
    enddo
!
    !$OMP PARALLEL DEFAULT(NONE) &
    !$OMP& SHARED(n_prisms,verts_rphiz,tetra_grid,grid_size,tetra_indices_per_prism,prism_volumes) &
    !$OMP& PRIVATE(r_values,z_values,rmin,triangle_indices,r_values_intermediate,z_values_intermediate, &
    !$OMP& r,z,gradient,z_star,alpha,beta)
    !$OMP DO
    do i = 1,n_prisms
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
        gradient(1) = z(1)/r(1)
        gradient(2) = z(2)/r(2)
        gradient(3) = (z(2)-z(1))/(r(2)-r(1))
!
        if (r(1).eq.0) then
            gradient(1) = 0
        else
            gradient(1) = z(1)/r(1)
        endif
        if (r(2).eq.0) then
            gradient(2) = 0
        else
            gradient(2) = z(2)/r(2)
        endif
        if ((r(2)-r(1)).eq.0) then
            gradient(3) = 0
        else
            gradient(3) = (z(2)-z(1))/(r(2)-r(1))
        endif
!
        z_star = z(1) - r(1)*gradient(3)
        alpha = abs(gradient(1)-gradient(2))
        beta = gradient(3) - gradient(2)
!
        prism_volumes(i) =  2*pi/grid_size(2)*(alpha/3*r(1)**3+alpha*rmin/2*r(1)**2+ & 
                            abs(z_star*rmin*(r(2)-r(1))+(z_star+beta*rmin)/2*(r(2)**2-r(1)**2)+beta/3*(r(2)**3-r(1)**3)))
!
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    print*, 'end prism volume computation'
!
end subroutine calc_prism_volumes
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine fourier_transform_moments(n_fourier_modes)
!
    use constants, only: pi
    use tetra_grid_settings_mod, only: grid_size
!
    implicit none 
!
    integer, intent(inout)                                         :: n_fourier_modes
    integer                                                        :: n,m,j,k,p,q,n_triangles
    complex                                                        :: i
    double precision, dimension(:,:,:,:), allocatable              :: prism_moments_ordered_for_ft
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
    i = (0,1)
!
    do p = 1,n_triangles
        do n = 1,n_moments
            do m = 1,n_species
                do k = 0,n_fourier_modes-1
                    do j = 0,grid_size(ind_b)-1
                        moments_in_frequency_space(n,p,k+1,m) = moments_in_frequency_space(n,p,k+1,m) + &
                                                            & exp(-2*pi*i*j*k/n_prisms)*prism_moments_ordered_for_ft(n,p,j+1,m)
                    enddo
                enddo
            enddo
        enddo
    enddo
!
    deallocate(prism_moments_ordered_for_ft)
!
print*, 'fourier transform finished'
!
end subroutine fourier_transform_moments
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