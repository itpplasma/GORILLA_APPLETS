module boltzmann_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

    private

    real(dp) :: time_step,energy_eV, min_poloidal_flux, max_poloidal_flux, amax
    integer, dimension(:,:), allocatable :: tetra_indices_per_prism
    real(dp), dimension(:), allocatable :: prism_volumes, refined_prism_volumes, elec_pot_vec, n_b
    real(dp), dimension(:,:), allocatable :: verts, r_integrand_constants
    double complex, dimension(:,:), allocatable :: tetr_moments, prism_moments, &
                                                   & prism_moments_squared
    double complex, dimension(:,:,:), allocatable :: moments_in_frequency_space
    integer :: i_integrator_type, seed_option, n_moments, n_prisms, num_particles,&
               & ind_a, ind_b, ind_c, n_pushings, counter_phi_0_mappings, lost_outside, lost_inside
    real(dp) :: n_particles, density, constant_part_of_weights
    integer, dimension(4) :: moments_selector
    double complex, dimension(:,:), allocatable :: weights
    real(dp), dimension(:), allocatable :: J_perp, poloidal_flux, temperature_vector
    logical :: boole_refined_sqrt_g, boole_boltzmann_energies
    integer :: n_fourier_modes, n_triangles
    logical :: boole_linear_density_simulation, boole_antithetic_variate, boole_linear_temperature_simulation
    logical :: boole_collisions, boole_squared_moments, boole_point_source, boole_precalc_collisions
    real(dp), dimension(:,:,:), allocatable :: randcol
    integer :: randcoli = int(1.0d5)
    integer :: et_unit, rp_unit, pm_unit, di_unit
    type boole_writing_data_t
        logical :: vertex_indices = .true.
        logical :: vertex_coordinates = .true.
        logical :: prism_volumes = .false.
        logical :: refined_prism_volumes = .true.
        logical :: boltzmann_density = .true.
        logical :: electric_potential = .false.
        logical :: moments = .false.
        logical :: fourier_moments = .false.
    end type boole_writing_data_t
    type filenames_t
        character(len=100) :: exit_times
        character(len=100) :: remaining_particles
        character(len=100) :: poincare_maps
        character(len=100) :: prism_moments
        character(len=100) :: prism_moments_summed_squares
        character(len=100) :: vertex_coordinates
        character(len=100) :: vertex_indices
        character(len=100) :: prism_volumes
        character(len=100) :: fourier_moments
        character(len=100) :: refined_prism_volumes
        character(len=100) :: electric_potential
        character(len=100) :: boltzmann_densities
        character(len=100) :: divertor_intersections
        character(len=100) :: tetr_moments
    end type filenames_t

    !Namelist for boltzmann input
    NAMELIST /boltzmann_nml/ time_step,energy_eV,n_particles,boole_squared_moments,boole_point_source,boole_collisions, &
    & boole_precalc_collisions,density,boole_refined_sqrt_g,boole_boltzmann_energies, boole_linear_density_simulation, &
    & boole_antithetic_variate,boole_linear_temperature_simulation,i_integrator_type,seed_option

    public :: calc_boltzmann
   
contains

subroutine calc_boltzmann

    use orbit_timestep_gorilla_mod, only: initialize_gorilla
    use constants, only: ev2erg,pi,echarge,ame,amp,clight
    use tetra_physics_mod, only: particle_mass,particle_charge,cm_over_e,mag_axis_R0, coord_system, tetra_physics
    use omp_lib, only: omp_get_thread_num, omp_get_num_threads, omp_set_num_threads
    use parmot_mod, only: rmu,ro0
    use velo_mod, only: isw_field_type
    use supporting_functions_mod, only: theta_sym_flux2theta_vmec,theta_vmec2theta_sym_flux
    use tetra_grid_settings_mod, only: n_field_periods, grid_size
    use tetra_grid_mod, only: ntetr, nvert, verts_rphiz, tetra_grid, verts_sthetaphi
    use gorilla_settings_mod, only: boole_array_optional_quantities, ispecies
    use collis_ions, only: collis_init, stost
    use find_tetra_mod, only: find_tetra
    use gorilla_applets_settings_mod, only: i_option
    use field_mod, only: ipert
    use volume_integrals_and_sqrt_g_mod, only: calc_square_root_g, calc_volume_integrals

    real(dp), dimension(:,:), allocatable :: start_pos_pitch_mat, dens_mat, temp_mat, vpar_mat, efcolf_mat, &
                                                     velrat_mat, enrat_mat, dens_mat_tetr, temp_mat_tetr
    real(dp) :: v0,pitchpar,vpar,vperp,t_remain,t_confined, v, maxcol
    integer :: kpart,i,j,n,l,m,k,p,ind_tetr,iface,n_lost_particles,ierr,err,iantithetic, num_background_species, inorout
    integer :: n_start, n_end, i_part, count_integration_steps
    real(dp), dimension(3) :: x_rand_beg,x,randnum
    logical :: boole_initialized,boole_particle_lost
    real(dp) :: dtau, dphi,dtaumin, t_step
    real(dp), dimension(5) :: z, zet
    Character(LEN=50) :: format_moments, format_fourier_moments
    double complex, dimension(:,:), allocatable :: single_particle_tetr_moments
    real(dp) :: m0,z0,hamiltonian_time
    real(dp), dimension(:), allocatable :: efcolf,velrat,enrat,vpar_background,mass_num,charge_num,dens,temp
    integer :: Te_unit, Ti_unit, ne_unit, t_moments_unit, ipert_unit
    type(filenames_t) :: filenames
    real(dp), dimension(:,:), allocatable :: sqrt_g

    !Load input for boltzmann computation
    call read_boltzmann_inp()

    !call read_in_starting_conditions(start_pos_pitch_mat, num_particles)

    num_particles = int(n_particles)
    n_start = 1
    n_end = num_particles
    n_fourier_modes = 5

    !prepare moment calculation
    n_moments = 0
    moments_selector = 0
    do i = 1,size(boole_array_optional_quantities)
        if (boole_array_optional_quantities(i).eqv..true.) then
            n_moments = n_moments + 1
            moments_selector(n_moments) = i
        endif
    enddo

    open(newunit = ipert_unit, file='field_divB0.inp')
    read(ipert_unit,*) ipert        ! 0=eq only, 1=vac, 2=vac+plas no derivatives,
    close(ipert_unit)

    !Initialize GORILLA
    call initialize_gorilla(i_option,ipert)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! delete this again afterwards !!!!!!!!!!!!!!!!!!!!!!!
    if (ispecies.eq.4) particle_charge = 15*echarge
    print*, 'particle charge number = ', particle_charge/echarge
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    n_prisms = ntetr/3
    allocate(tetr_moments(n_moments,ntetr))
    allocate(single_particle_tetr_moments(n_moments,ntetr))
    if (boole_squared_moments) allocate(prism_moments_squared(n_moments,n_prisms))
    allocate(prism_moments(n_moments,n_prisms))
    allocate(weights(num_particles,1))
    allocate(J_perp(num_particles))
    allocate(poloidal_flux(num_particles))
    allocate(temperature_vector(num_particles))
    allocate(elec_pot_vec(n_prisms))
    allocate(prism_volumes(n_prisms))
    allocate(refined_prism_volumes(n_prisms))
    allocate(n_b(n_prisms))
    allocate(tetra_indices_per_prism(n_prisms,3))

    do i = 1,3
        tetra_indices_per_prism(:,i) = (/(i+3*k,k = 0,n_prisms-1)/)
    enddo

    poloidal_flux = 0
    temperature_vector = 0
    elec_pot_vec = 0
    n_b = 0 !boltzmann density for uniform spatial distribution

    call calc_square_root_g(sqrt_g)

print*, 'start calc_volume_integrals'
    call calc_volume_integrals(boole_boltzmann_energies,boole_refined_sqrt_g, density, energy_eV, &
                               prism_volumes, refined_prism_volumes, elec_pot_vec, n_b)
print*, 'end calc_volume_integrals'

    !Compute velocity module from kinetic energy dependent on particle species
    v0=sqrt(2.d0*energy_eV*ev2erg/particle_mass)

print*, 'calc_starting_conditions started'
    call calc_starting_conditions(v0,start_pos_pitch_mat)
print*, 'calc_starting_conditions finished'

    !compute maximum poloidal flux
    max_poloidal_flux = 0
    min_poloidal_flux = tetra_physics(1)%Aphi1
    do i = 1, ntetr
        max_poloidal_flux = max(max_poloidal_flux,tetra_physics(i)%Aphi1 + sum(tetra_physics(i)%gAphi* &
                            & (verts([1,2,3],tetra_grid(i)%ind_knot(4))-verts([1,2,3],tetra_grid(i)%ind_knot(1)))))
        min_poloidal_flux = min(min_poloidal_flux,tetra_physics(i)%Aphi1)
    enddo

    if (boole_collisions) then
        num_background_species = 2 
        allocate(dens_mat(num_background_species-1,size(verts_rphiz(1,:))))
        allocate(temp_mat(num_background_species,size(verts_rphiz(1,:))))
        allocate(vpar_mat(num_background_species,ntetr))
        allocate(efcolf_mat(num_background_species,ntetr))
        allocate(velrat_mat(num_background_species,ntetr))
        allocate(enrat_mat(num_background_species,ntetr))
        allocate(mass_num(num_background_species-1))
        allocate(charge_num(num_background_species-1))
        allocate(dens(num_background_species))
        allocate(temp(num_background_species))
        allocate(efcolf(num_background_species))
        allocate(velrat(num_background_species))
        allocate(enrat(num_background_species))
        mass_num = 0
        charge_num = 0
        mass_num(1) = 2
        !mass_num(2) = 3
        charge_num(1) = 1
        !charge_num(2) = 2
        dens_mat = 5.0d13
        temp_mat = energy_eV
        vpar_mat = 0 !ask Sergei when this will be needed!!!
        m0 = particle_mass/amp
        z0 = particle_charge/echarge
        print*, 'm0 = ', m0, 'z0 = ', z0

!!!!!!!!!!!!!!!!!!!! Alternative route is taken because data is not available per vertex but per tetrahedron !!!!!!!!!!!!!!!!!!!!!!!

        allocate(dens_mat_tetr(num_background_species-1,ntetr))
        allocate(temp_mat_tetr(num_background_species,ntetr))

        open(newunit = Te_unit, file = 'background/Te_d.dat')
        read(Te_unit,'(e16.9)') (temp_mat_tetr(2,i),i=1,ntetr/grid_size(2),3)
        close(Te_unit)

        open(newunit = Ti_unit, file = 'background/Ti_d.dat')
        read(Ti_unit,'(e16.9)') (temp_mat_tetr(1,i),i=1,ntetr/grid_size(2),3)
        close(Ti_unit)

        open(newunit = ne_unit, file = 'background/ne_d.dat')
        read(ne_unit,'(e16.9)') (dens_mat_tetr(1,i),i=1,ntetr/grid_size(2),3)
        close(ne_unit)

        do i = 1,grid_size(2)-1
            temp_mat_tetr(:,i*ntetr/grid_size(2)+1:(i+1)*ntetr/grid_size(2):3) = temp_mat_tetr(:,1:ntetr/grid_size(2):3)
            dens_mat_tetr(:,i*ntetr/grid_size(2)+1:(i+1)*ntetr/grid_size(2):3) = dens_mat_tetr(:,1:ntetr/grid_size(2):3)
        enddo
        do i = 1,2
            temp_mat_tetr(:,1+i:ntetr:3) = temp_mat_tetr(:,1:ntetr:3)
            dens_mat_tetr(:,1+i:ntetr:3) = dens_mat_tetr(:,1:ntetr:3)
        enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        do i = 1, ntetr
            do j = 1,num_background_species
                !following if statement because electron density will be calculated in collis init
                if (j.lt.num_background_species) dens(j) = sum(dens_mat_tetr(j,tetra_grid(i)%ind_knot([1,2,3,4])))/4
                temp(j) = sum(temp_mat_tetr(j,tetra_grid(i)%ind_knot([1,2,3,4])))/4
            enddo
            call collis_init(m0,z0,mass_num,charge_num,dens,temp,energy_eV,v0,efcolf,velrat,enrat)
            efcolf_mat(:,i) = efcolf
            velrat_mat(:,i) = velrat
            enrat_mat(:,i) = enrat
        enddo
    endif

        kpart = 0
        n_lost_particles = 0
        maxcol = 0
        lost_outside = 0
        lost_inside = 0
        n_pushings = 0
        counter_phi_0_mappings = 0
        tetr_moments = 0
        if (boole_squared_moments) prism_moments_squared = 0
        prism_moments = 0
        iantithetic = 1
        if (boole_antithetic_variate) iantithetic = 2
        count_integration_steps = 0

        call name_files(filenames)
        call unlink_files(filenames)
        call open_files_before_main_loop

        if (boole_collisions) deallocate(efcolf,velrat,enrat)

        !$OMP PARALLEL DEFAULT(NONE) &
        !$OMP& SHARED(num_particles,kpart,v0,time_step,i_integrator_type,boole_collisions,sqrt_g, &
        !$OMP& dtau,dtaumin,n_start,n_end,tetra_indices_per_prism, prism_moments_squared,boole_squared_moments, &
        !$OMP& start_pos_pitch_mat,boole_boltzmann_energies, prism_moments,count_integration_steps, et_unit, &
        !$OMP& density,energy_eV,dens_mat,temp_mat,vpar_mat,tetra_grid,iantithetic,tetra_physics, di_unit, pm_unit, rp_unit, &
        !$OMP& efcolf_mat,velrat_mat,enrat_mat,num_background_species,randcol,randcoli,maxcol,boole_precalc_collisions) &
        !$OMP& FIRSTPRIVATE(particle_mass, particle_charge) &
        !$OMP& PRIVATE(p,l,n,boole_particle_lost,x_rand_beg,x,pitchpar,vpar,vperp,boole_initialized,t_step,err,zet, &
        !$OMP& ind_tetr,iface,t_remain,t_confined,z,ierr, v, single_particle_tetr_moments,hamiltonian_time, &
        !$OMP& m0,z0,i,efcolf,velrat,enrat,vpar_background,inorout,randnum) &
        !$OMP& REDUCTION(+:n_lost_particles,tetr_moments, n_pushings, counter_phi_0_mappings)
        print*, 'get number of threads', omp_get_num_threads()
        if (boole_collisions) allocate(efcolf(num_background_species),velrat(num_background_species),enrat(num_background_species))
        !$OMP DO

        !Loop over particles
        do p = n_start,n_end/iantithetic !1,num_particles/iantithetic
            do l = 1,iantithetic

                n = (p-1)*iantithetic+l
                !$omp critical
                !Counter for particles
                kpart = kpart+1 !in general not equal to n becuase of parallelisation
                boole_particle_lost = .false.
if (n_end.gt.10) then
    if (modulo(kpart,int(n_end/10)).eq.0) then
        print *, kpart, ' / ', num_particles, 'particle: ', n, 'thread: ' !, omp_get_thread_num()
    endif
else 
    print *, kpart, ' / ', num_particles, 'particle: ', n, 'thread: ' !, omp_get_thread_num()
endif
                !$omp end critical

!if (abs(start_pos_pitch_mat(4,n)).lt.0.5) cycle !only trace particles for which vpar > vperp
                t_step = time_step
                t_confined = 0
                if (l.eq.1) single_particle_tetr_moments = 0

                !You need x_rand_beg(1,3), pitchpar(1) (between -1 and 1), energy is already given
                x_rand_beg = start_pos_pitch_mat(1:3,n)
                pitchpar = start_pos_pitch_mat(4,n)
                ! print*, start_pos_pitch_mat(5,n)
                ! if (n.eq.972) print*, x_rand_beg,start_pos_pitch_mat(5,n),pitchpar

                if (i_integrator_type.eq.2) print*, 'Error: i_integratpr_type set to 2, this module only works with &
                                                    & i_integrator_type set to 1'
                x = x_rand_beg
                vpar = pitchpar * v0
                vperp = sqrt(v0**2-vpar**2)
                if (boole_boltzmann_energies) then
                    v = sqrt(start_pos_pitch_mat(5,n)*ev2erg*2/particle_mass)
                    vpar = pitchpar * v
                    vperp = sqrt(v**2-vpar**2)
                endif

                i = 0
                do while (t_confined.lt.time_step)
                    i = i+1
                    !Orbit integration
                    if (i.eq.1) then
                        boole_initialized = .false.
                    endif
                    if (boole_collisions) then
                        if (i.eq.1) call find_tetra(x,vpar,vperp,ind_tetr,iface)
                        if (.not.(ind_tetr.eq.-1)) then
                            efcolf = efcolf_mat(:,ind_tetr)
                            velrat = velrat_mat(:,ind_tetr)
                            enrat = enrat_mat(:,ind_tetr)
                            vpar_background = vpar_mat(:,ind_tetr)
                            !print*, vpar_background

                            vpar = vpar - vpar_background(1)
                            !since vpar_background actually has num_background_particles entries, consider giving it as an extra
                            !optional input variable to stost, before randnum (maybe also check if radnum could then be set by 
                            !randnum = variable eve if vpar_background is not set and other variables are not set by name indexing)
                            !since it came up when writing these lines: replace expressions like
                            !"verts(size(verts_rphiz(:,1)),size(verts_rphiz(1,:)))" with "3,nvert"
                            zet(1:3) = x !spatial position
                            zet(4) = sqrt(vpar**2+vperp**2)/v0 !normalized velocity module 
                            zet(5) = vpar/sqrt(vpar**2+vperp**2) !pitch parameter

                            if (boole_precalc_collisions) then
                                randnum = randcol(n,mod(i-1,randcoli)+1,:) 
                                call stost(efcolf,velrat,enrat,zet,t_step,1,err,(time_step-t_confined)*v0,randnum)
                            else
                                call stost(efcolf,velrat,enrat,zet,t_step,1,err,(time_step-t_confined)*v0)
                            endif

                            t_step = t_step/v0
                            x = zet(1:3)
                            vpar = zet(5)*zet(4)*v0+vpar_background(1)
                            vperp = sqrt(1-zet(5)**2)*zet(4)*v0

                            !optionally still change particle_mass, particle_charge and cm_over_e, e.g.:
                            !particle_charge = particle_charge + echarge
                            !particle_mass = particle_mass + ame - amp 
                            !cm_over_e = clight*particle_mass/particle_charge
                        endif
                    endif

                    call orbit_timestep_gorilla_boltzmann(x,vpar,vperp,t_step,boole_initialized,ind_tetr,iface,n,v, &
                                    & sqrt_g,start_pos_pitch_mat,single_particle_tetr_moments, hamiltonian_time,inorout,t_remain)

                    t_confined = t_confined + t_step - t_remain
                    !Lost particle handling
                    if(ind_tetr.eq.-1) then
!write another if clause (if hole size = minimal .and. particle lost inside .and. boole_cut_out_hole = .true. 
!(use extra variable m in orbit routine (0 normally, 1 when lost outside, -1 when lost inside))),
!if m = 1 do as now, if m = -1 select arbitrary newp position and update x, vpar and vperp)
                        write(et_unit,*) t_confined, x, n
                        !print*, t_confined, x
                        !$omp critical
                            n_lost_particles = n_lost_particles + 1
                            boole_particle_lost = .true.
                        !$omp end critical
                        exit
                    endif

                    v = sqrt(vpar**2+vperp**2)
                enddo
                !$omp critical
                !print*, i
                count_integration_steps = count_integration_steps + i
                !print*, dble(i)/dble(randcoli)
                maxcol = max(dble(i)/dble(randcoli),maxcol)
                !Write results in file
                ! !$omp critical
                !     write(949,*) n, boole_particle_lost , x_rand_beg ,pitchpar,x(1),t_confined
                ! !$omp end critical
                !$omp end critical
                if (t_confined.eq.time_step) then
                    write(rp_unit,*) x, v, vpar, vperp, i, n
                endif
            enddo

                if (boole_squared_moments) then
                    !$omp critical
                    prism_moments_squared = prism_moments_squared + &
                                        & (single_particle_tetr_moments(:,tetra_indices_per_prism(:,1)) + &
                                        &  single_particle_tetr_moments(:,tetra_indices_per_prism(:,2)) + &
                                        &  single_particle_tetr_moments(:,tetra_indices_per_prism(:,3)))**2
                    prism_moments = prism_moments + &
                                & (single_particle_tetr_moments(:,tetra_indices_per_prism(:,1)) + &
                                &  single_particle_tetr_moments(:,tetra_indices_per_prism(:,2)) + &
                                &  single_particle_tetr_moments(:,tetra_indices_per_prism(:,3)))
                    !$omp end critical
                endif
        enddo !n
        !$OMP END DO
        !$OMP END PARALLEL


        if(.not.boole_squared_moments) then
            prism_moments = (tetr_moments(:,tetra_indices_per_prism(:,1)) + &
                        & tetr_moments(:,tetra_indices_per_prism(:,2)) + &
                        & tetr_moments(:,tetra_indices_per_prism(:,3)))
        endif

        do n = 1,n_moments
            prism_moments(n,:) = prism_moments(n,:)/(prism_volumes*time_step*n_particles) !do normalisations
            if (boole_squared_moments) then
                prism_moments_squared(n,:) = prism_moments_squared(n,:)/(prism_volumes**2*time_step**2*n_particles) !do normalisations
            endif
            if (boole_refined_sqrt_g) then
                    prism_moments(n,:) = prism_moments(n,:)*prism_volumes/refined_prism_volumes
                    if (boole_squared_moments) then
                    prism_moments_squared(n,:) = prism_moments_squared(n,:)*prism_volumes**2/refined_prism_volumes**2
                    endif
            endif
        enddo

    call fourier_transform_moments

    if (boole_precalc_collisions) print*, "maxcol = ", maxcol
    print*, 'Number of lost particles',n_lost_particles
    print*, 'max_poloidal_flux is', max_poloidal_flux
    print*, 'min_poloidal_flux is', min_poloidal_flux
    print*, 'average number of pushings = ', n_pushings/n_particles
    print*, 'average number of toroidal revolutions = ', counter_phi_0_mappings/n_particles
    print*, 'average number of integration steps = ', count_integration_steps/n_particles

    call close_files
    call write_data_to_files(filenames)


PRINT*, 'particle mass = ', particle_mass
PRINT*, 'absolute value of velocity = ', v0
PRINT*, 'particle charge = ', particle_charge
PRINT*, 'temperature = ', ev2erg*energy_eV
print*, 'energy in eV = ', energy_eV
print*, 'tracing time in seconds = ', time_step
print*, 'number of particles left through the outside = ', lost_outside
print*, 'number of particles left through the inside = ', lost_inside

deallocate(start_pos_pitch_mat, tetr_moments, prism_moments, single_particle_tetr_moments)
if (boole_squared_moments) deallocate(prism_moments_squared)

end subroutine calc_boltzmann

subroutine read_boltzmann_inp()

    integer :: b_inp_unit
    open(newunit = b_inp_unit, file='boltzmann.inp', status='unknown')
    read(b_inp_unit,nml=boltzmann_nml)
    close(b_inp_unit)
  
    print *,'GORILLA: Loaded input data from boltzmann.inp'
end subroutine read_boltzmann_inp

subroutine calc_starting_conditions(v0,start_pos_pitch_mat)
    
    use constants, only: pi, ev2erg
    use tetra_grid_mod, only: verts_rphiz, verts_sthetaphi, ntetr
    use find_tetra_mod, only: find_tetra
    use tetra_grid_settings_mod, only: grid_kind
    use tetra_physics_mod, only: coord_system
    use collis_ions, only: collis_init, stost

    real(dp), intent(in)                                   :: v0
    real(dp), dimension(:,:), allocatable, intent(out)     :: start_pos_pitch_mat
    real(dp)                                               :: rand_scalar, vpar, vperp
    real(dp)                                               :: amin, cmin, cmax !amax is set globally
    real(dp), dimension(:), allocatable                    :: rand_vector
    real(dp), dimension(:,:), allocatable                  :: rand_matrix1, rand_matrix2
    integer                                                        :: i
    real(dp), dimension(3)                                 :: x
    integer                                                        :: ind_tetr_out,iface,seed_inp_unit

!!!!comment out the following section to make starting conditions really random!!!

    integer,dimension(:), allocatable                              :: seed
    integer                                                        :: n

    open(newunit = seed_inp_unit, file='seed.inp', status='old',action = 'read')
    read(seed_inp_unit,*) n
    allocate(seed(n))
    read(seed_inp_unit,*) seed
    close(seed_inp_unit)
    CALL RANDOM_SEED (PUT=seed)
    deallocate(seed)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    allocate(start_pos_pitch_mat(5,num_particles))
    allocate(rand_vector(num_particles))
    allocate(rand_matrix1(3,num_particles))
    allocate(rand_matrix2(2,num_particles))

    start_pos_pitch_mat = 0

    ind_a = 1 !(R in cylindrical coordinates, s in flux coordinates)
    ind_b = 2 !(phi in cylindrical and flux coordinates)
    ind_c = 3 !(z in cylindrical coordinates, theta in flux coordinates)

    if (coord_system.eq.2) then
        ind_b = 3
        ind_c = 2
    endif

    if (coord_system.eq.1) allocate(verts(size(verts_rphiz(:,1)),size(verts_rphiz(1,:))))
    if (coord_system.eq.2) allocate(verts(size(verts_sthetaphi(:,1)),size(verts_sthetaphi(1,:))))
    if (coord_system.eq.1) verts = verts_rphiz
    if (coord_system.eq.2) verts = verts_sthetaphi

    amin = minval(verts(ind_a,:))
    amax = maxval(verts(ind_a,:))
    cmin = minval(verts(ind_c,:))
    cmax = maxval(verts(ind_c,:))

    constant_part_of_weights = density*(amax-amin)*(cmax-cmin)*2*pi

    !compute starting conditions
    if (boole_point_source) then
        if (grid_kind.eq.2) then
            start_pos_pitch_mat(1,:) = 160
            start_pos_pitch_mat(2,:) = 0
            start_pos_pitch_mat(3,:) = 70
        elseif (grid_kind.eq.4) then
            start_pos_pitch_mat(1,:) = 205
            start_pos_pitch_mat(2,:) = 0
            start_pos_pitch_mat(3,:) = 0
        endif
        if (coord_system.eq.2) print*, 'error: point source is only implemented for cylindrical coordinate system'
    else
        call RANDOM_NUMBER(rand_matrix1)
        start_pos_pitch_mat(ind_a,:) = amin + (amax - amin)*rand_matrix1(ind_a,:) !r in cylindrical, s in flux coordinates
        start_pos_pitch_mat(ind_b,:) = 2*pi*rand_matrix1(ind_b,:) !phi in cylindrical and flux coordinates
        start_pos_pitch_mat(ind_c,:) = cmin + (cmax - cmin)*rand_matrix1(ind_c,:) !z in cylindrical, theta in flux coordinates
        ! start_pos_pitch_mat(ind_a,:) = (/(214 + i*(216-214)/n_particles, i=1,num_particles)/)!r
        ! start_pos_pitch_mat(ind_b,:) = 0.0d0  !1d-1 !phi in cylindrical and flux coordinates
        ! start_pos_pitch_mat(ind_c,:) = 12d0 !z in cylindrical, theta in flux coordinates
    endif

    call RANDOM_NUMBER(rand_matrix2)
    !start_pos_pitch_mat(4,:) = 2*rand_matrix2(1,:)-1 !pitch parameter
    start_pos_pitch_mat(4,:) = 0.7d0 !pitch parameter

    if (boole_boltzmann_energies) then !compare with equation 133 of master thesis of Jonatan Schatzlmayr (remaining parts will be added later)
        start_pos_pitch_mat(5,:) = 5*energy_eV*rand_matrix2(2,:) !boltzmann energy distribution
        constant_part_of_weights = constant_part_of_weights*10/sqrt(pi)*energy_eV*ev2erg
    endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! start antithetic variate !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (boole_antithetic_variate) then
        start_pos_pitch_mat(:,1:num_particles:2) = start_pos_pitch_mat(:,2:num_particles:2)
        start_pos_pitch_mat(4,1:num_particles:2) = -start_pos_pitch_mat(4,2:num_particles:2)
    endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! end antithetic variate !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    weights(:,1) = constant_part_of_weights
    if (boole_refined_sqrt_g.eqv..false.) weights(:,1) = constant_part_of_weights*start_pos_pitch_mat(ind_a,:)

    if (boole_precalc_collisions) then
        allocate(randcol(num_particles,randcoli,3))
        call RANDOM_NUMBER(randcol)
        !3.464102 = sqrt(12), this creates a random number with zero average and unit variance
        randcol(:,:,1:2:3) =  3.464102*(randcol(:,:,1:2:3)-.5)
    endif
end subroutine calc_starting_conditions

subroutine orbit_timestep_gorilla_boltzmann(x,vpar,vperp,t_step,boole_initialized,ind_tetr,iface, n,v,sqrt_g,start_pos_pitch_mat, &
                                          & single_particle_tetr_moments,hamiltonian_time, inorout, t_remain_out)

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
    use tetra_grid_mod, only: tetra_grid, ntetr

    real(dp), dimension(3), intent(inout)   :: x
    real(dp), intent(inout)                 :: vpar,vperp
    real(dp), intent(in)                    :: t_step
    logical, intent(inout)                          :: boole_initialized
    integer, intent(inout)                          :: ind_tetr,iface
    real(dp), intent(out), optional         :: t_remain_out
    real(dp), dimension(3)                  :: z_save, x_save
    real(dp)                                :: vperp2,t_remain,t_pass,vpar_save, v, hamiltonian_time, aphi
    logical                                         :: boole_t_finished
    integer                                         :: ind_tetr_save,iper_phi,k, m, n, inorout
    real(dp)                                :: perpinv,perpinv2, speed, r, z, phi, B, phi_elec_func
    real(dp), dimension(:,:)                :: start_pos_pitch_mat
    type(optional_quantities_type)                  :: optional_quantities
    double complex, dimension(:,:)                  :: single_particle_tetr_moments
    integer                                         :: single_particle_counter_phi0_mappings
    real(dp), dimension(:,:)                :: sqrt_g
    

    !If orbit_timestep is called for the first time without grid position
    if(.not.boole_initialized) then

        !Check coordinate domain (optionally perform modulo operation)
        call check_coordinate_domain(x)

        !Find tetrahedron index and face index for position x
        call find_tetra(x,vpar,vperp,ind_tetr,iface)

        !If particle doesn't lie inside any tetrahedron
        if(ind_tetr.eq.-1) then
            t_remain_out = t_step
            return
        endif

        r = x(1) - verts(1, tetra_grid(ind_tetr)%ind_knot(1))
        z = x(3) - verts(3,tetra_grid(ind_tetr)%ind_knot(1))
        phi = x(2) - verts(2,tetra_grid(ind_tetr)%ind_knot(1))

        if (boole_refined_sqrt_g) then
            weights(n,1) = weights(n,1)* &
                                    &  (sqrt_g(ind_tetr,1)+r*sqrt_g(ind_tetr,2)+z*sqrt_g(ind_tetr,3))/ &
                                    &  (sqrt_g(ind_tetr,4)+r*sqrt_g(ind_tetr,5)+z*sqrt_g(ind_tetr,6))
        endif

        if (boole_linear_density_simulation.or.boole_linear_temperature_simulation) then
            poloidal_flux(n) = tetra_physics(ind_tetr)%Aphi1 + sum(tetra_physics(ind_tetr)%gAphi*(/r,phi,z/)) + &
                                             & cm_over_e*vpar*&
                                             & (tetra_physics(ind_tetr)%h2_1+sum(tetra_physics(ind_tetr)%gh2*(/r,phi,z/)))
        endif
        if (boole_linear_density_simulation) then
            weights(n,1) = weights(n,1)*(max_poloidal_flux*1.1-poloidal_flux(n))/(max_poloidal_flux*1.1)
        endif
      
        if (boole_boltzmann_energies) then
            !compare with equation 133 of master thesis of Jonatan Schatzlmayr (remaining parts have been added before)
            phi_elec_func = tetra_physics(ind_tetr)%Phi1 + sum(tetra_physics(ind_tetr)%gPhi*(/r,phi,z/))
            if (.not. boole_linear_temperature_simulation) then
                weights(n,1) = weights(n,1)*sqrt(start_pos_pitch_mat(5,n)*ev2erg)/(energy_eV*ev2erg)**1.5* &
                            & exp(-(start_pos_pitch_mat(5,n)*ev2erg+particle_charge*phi_elec_func)/(energy_eV*ev2erg))
            else
                temperature_vector(n) = energy_eV*ev2erg*(max_poloidal_flux*1.1-poloidal_flux(n))/(max_poloidal_flux*1.1)
                weights(n,1) = weights(n,1)*sqrt(start_pos_pitch_mat(5,n)*ev2erg)/temperature_vector(n)**1.5* &
                & exp(-(start_pos_pitch_mat(5,n)*ev2erg+particle_charge*phi_elec_func)/temperature_vector(n))
            endif
        endif

        !compute J_perp for perpendicular pressure
        B = tetra_physics(ind_tetr)%bmod1+sum(tetra_physics(ind_tetr)%gb*(/r,phi,z/))
        J_perp(n) = particle_mass*v**2*(1-start_pos_pitch_mat(4,n)**2)*cm_over_e/(2*B)*(-1)
        !-1 because of negative gyrophase

        boole_initialized = .true.
    endif
          
    !Exit the subroutine after initialization, if time step equals zero
    if(t_step.eq.0.d0) return

    inorout = 0

    !Squared perpendicular velocity
    vperp2 = vperp**2

    !Compute relative particle position
    z_save = x-tetra_physics(ind_tetr)%x1

    !Compute perpendicular invariant of particle
    perpinv=-0.5d0*vperp2/bmod_func(z_save,ind_tetr)
    perpinv2 = perpinv**2
             
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

    !Logical for handling time integration
    boole_t_finished = .false.

    single_particle_counter_phi0_mappings = 0

    n_pushings = n_pushings-1 !set n_pushings to -1 because when entering the loop it wil go back to one without pushing
    !Loop for tetrahedron pushings until t_step is reached
    do

        n_pushings = n_pushings + 1
        !Domain Boundary
        if(ind_tetr.eq.-1) then
            !print *, 'WARNING: Particle lost.'
            aphi = tetra_physics(ind_tetr_save)%Aphi1+sum(tetra_physics(ind_tetr_save)%gAphi*z_save)
            !$omp critical
            if (abs(aphi-max_poloidal_flux).lt.abs(aphi-min_poloidal_flux)) then
                lost_outside = lost_outside + 1
                inorout = 1
            endif
            if (abs(aphi-max_poloidal_flux).gt.abs(aphi-min_poloidal_flux)) then
                lost_inside = lost_inside + 1
                inorout = -1
            endif
            !$omp end critical
            if( present(t_remain_out)) then
                t_remain_out = t_remain
            endif
            exit
        endif
          
        !Save the tetrahedron index for computation of vperp in the last step
        ind_tetr_save = ind_tetr

        !Save vpar for the computation of the parallel adiabatic invariant
        vpar_save = vpar

        !t_remain (in) ... remaining time until t_step is finished
        !t_pass (out) ... time to pass the tetrahdron

        !Save x for the computation of the current
        x_save = x

        !Calculate trajectory
        select case(ipusher)
            case(1)
                call pusher_tetra_rk(ind_tetr,iface,x,vpar,z_save,t_remain,t_pass,boole_t_finished,iper_phi)
            case(2)
                call pusher_tetra_poly(poly_order,ind_tetr,iface,x,vpar,z_save,t_remain,&
                                                    & t_pass,boole_t_finished,iper_phi,optional_quantities)
        end select

        t_remain = t_remain - t_pass
        if (iper_phi.ne.0) then
            counter_phi_0_mappings = counter_phi_0_mappings + iper_phi
            single_particle_counter_phi0_mappings = single_particle_counter_phi0_mappings + 1
            if ((single_particle_counter_phi0_mappings.gt.10) &
                .and.(single_particle_counter_phi0_mappings.le.3000)) then
                !$omp critical
                    write(pm_unit,*) x
                !$omp end critical
            endif
        endif

        if (x(3).lt.-105d0) then
            boole_t_finished = .true.
           if (single_particle_counter_phi0_mappings.gt.10) then
            call calc_plane_intersection(x_save,x,-105d0)
            !$omp critical
                write(di_unit,*) x, n
            !$omp end critical
           endif
        endif

        hamiltonian_time = 0

        if (.not.boole_squared_moments) then
            do m = 1,n_moments
            !print*, n_moments, moments_selector
                select case(moments_selector(m))
                    case(1)
                        tetr_moments(m,ind_tetr_save) = tetr_moments(m,ind_tetr_save) + weights(n,1)* &
                                                        & optional_quantities%t_hamiltonian!* &
                                                        !& (exp(2*(0,1)*x(2))+exp(3*(0,1)*x(2)))
                        hamiltonian_time = optional_quantities%t_hamiltonian
                    case(2)
                        tetr_moments(m,ind_tetr_save) = tetr_moments(m,ind_tetr_save) + weights(n,1)* &
                                                        & optional_quantities%gyrophase*J_perp(n)
                    case(3)
                        tetr_moments(m,ind_tetr_save) = tetr_moments(m,ind_tetr_save) + weights(n,1)* &
                                                        & optional_quantities%vpar_int
                    case(4)
                        tetr_moments(m,ind_tetr_save) = tetr_moments(m,ind_tetr_save) + weights(n,1)* &
                                                        & optional_quantities%vpar2_int
                end select
            enddo
        else
            do m = 1,n_moments
            !print*, n_moments, moments_selector
                select case(moments_selector(m))
                    case(1)
                        single_particle_tetr_moments(m,ind_tetr_save) = single_particle_tetr_moments(m,ind_tetr_save) + &
                                                                        & weights(n,1)*optional_quantities%t_hamiltonian!* &
                                                                       !& (exp(2*(0,1)*x(2))+exp(3*(0,1)*x(2)))
                        hamiltonian_time = optional_quantities%t_hamiltonian
                    case(2)
                        single_particle_tetr_moments(m,ind_tetr_save) = single_particle_tetr_moments(m,ind_tetr_save) + &
                                                                        & weights(n,1)*optional_quantities%gyrophase*J_perp(n)
                    case(3)
                        single_particle_tetr_moments(m,ind_tetr_save) = single_particle_tetr_moments(m,ind_tetr_save) + &
                                                                        & weights(n,1)*optional_quantities%vpar_int
                    case(4)
                        single_particle_tetr_moments(m,ind_tetr_save) = single_particle_tetr_moments(m,ind_tetr_save) + &
                                                                        & weights(n,1)*optional_quantities%vpar2_int
                end select
            enddo
        endif

        !Orbit stops within cell, because "flight"-time t_step has finished
        if(boole_t_finished) then
            if( present(t_remain_out)) then
                t_remain_out = t_remain
            endif
            exit
        endif

    enddo !Loop for tetrahedron pushings

    !Compute vperp from position
    vperp = vperp_func(z_save,perpinv,ind_tetr_save)
        
!             !NOT FULLY IMPLEMENTED YET: Deallocate precomputed quantities dependent on perpinv
!             call alloc_precomp_poly_perpinv(2,ntetr)
!print*, 'number of pushings is', n_pushings

end subroutine orbit_timestep_gorilla_boltzmann

subroutine calc_plane_intersection(x_save,x,z_plane)

    use constants, only : pi

    real(dp), dimension(3), intent(in) :: x_save
    real(dp), dimension(3), intent(inout) :: x
    real(dp), intent(in) :: z_plane
    real(dp) :: rel_dist_z

    rel_dist_z = (z_plane-x_save(3))/(x(3)-x_save(3))
    x(1) = x_save(1) + rel_dist_z*(x(1)-x_save(1))
    if (abs(x(2)-x_save(2)).gt.pi) then
        x(2) = modulo(x_save(2) + 2*pi-abs(x(2)-x_save(2)),2*pi)
    else
        x(2) = x_save(2) + rel_dist_z*(x(2)-x_save(2))
    endif
    x(3) = z_plane

end subroutine calc_plane_intersection

subroutine fourier_transform_moments

    use constants, only: pi
    use tetra_grid_settings_mod, only: grid_size

    integer                                                        :: n,m,j,k,p,q,l
    integer                                                        :: pmff_unit
    complex                                                        :: i
    complex, dimension(:,:,:), allocatable              :: prism_moments_ordered_for_ft

print*, 'fourier transform started'

    n_triangles = n_prisms/grid_size(2)
    allocate(prism_moments_ordered_for_ft(n_moments,n_triangles,grid_size(2)))
    do q = 1,grid_size(2)
        prism_moments_ordered_for_ft(:,:,q) = prism_moments(:,n_triangles*(q-1)+1:n_triangles*q)
    enddo

    if (n_fourier_modes.gt.grid_size(2)) then
        print*, 'n_fourier_modes was chosen to be bigger than n_phi, it is therefore reduced to n_phi'
        n_fourier_modes = grid_size(2)
    endif

    allocate(moments_in_frequency_space(n_moments,n_triangles,n_fourier_modes))
    moments_in_frequency_space = 0
    i = (0,1)

    do p = 1,n_triangles
        do n = 1,n_moments
            do k = 0,n_fourier_modes-1
                do j = 0,grid_size(2)-1
                    moments_in_frequency_space(n,p,k+1) = moments_in_frequency_space(n,p,k+1) + 1/dble(grid_size(2))* &
                                                    & exp(-2*pi*i*j*k/grid_size(2))*prism_moments_ordered_for_ft(n,p,j+1)
                enddo
            enddo
        enddo
    enddo

    open(newunit = pmff_unit, file = 'prism_moments_for_fourier.dat')
    do l = 1,n_triangles
        do m = 1,grid_size(2)-1
            write(pmff_unit,'(2ES20.10E4)',advance="no") real(prism_moments_ordered_for_ft(1,l,m)), &
                                                  aimag(prism_moments_ordered_for_ft(1,l,m))     
        enddo
            write(pmff_unit,'(2ES20.10E4)') real(prism_moments_ordered_for_ft(1,l,grid_size(2))), &
                                     aimag(prism_moments_ordered_for_ft(1,l,grid_size(2)))
    enddo
    close(pmff_unit)

    deallocate(prism_moments_ordered_for_ft)

print*, 'fourier transform finished'

end subroutine fourier_transform_moments

subroutine write_data_to_files(filenames)

    type(filenames_t) :: filenames
    type(boole_writing_data_t) :: boole_writing_data

    if (boole_writing_data%vertex_indices) call write_vertex_indices(filenames%vertex_indices)
    if (boole_writing_data%vertex_coordinates) call write_vertex_coordinates(filenames%vertex_coordinates)
    if (boole_writing_data%prism_volumes) call write_prism_volumes(filenames%prism_volumes)
    if (boole_writing_data%refined_prism_volumes) call write_refined_prism_volumes(filenames%refined_prism_volumes)
    if (boole_writing_data%boltzmann_density) call write_boltzmann_densities(filenames%boltzmann_densities)
    if (boole_writing_data%electric_potential) call write_electric_potential(filenames%electric_potential)
    if (boole_writing_data%moments) call write_moments(filenames%prism_moments, &
                                                       filenames%prism_moments_summed_squares, &
                                                       filenames%tetr_moments)
    if (boole_writing_data%fourier_moments) call write_fourier_moments(filenames%fourier_moments)

end subroutine write_data_to_files

subroutine write_vertex_indices(filename_vertex_indices)

    use tetra_grid_mod, only: ntetr, tetra_grid

    character(len=100) :: filename_vertex_indices
    integer :: vi_unit
    integer :: i

    open(newunit = vi_unit, file = filename_vertex_indices)
    do i=1, ntetr
        write(vi_unit, *) tetra_grid(i)%ind_knot([1, 2, 3, 4])
    end do
    close(vi_unit)

end subroutine write_vertex_indices

subroutine write_vertex_coordinates(filename_vertex_coordinates)

    use tetra_physics_mod, only: coord_system
    use tetra_grid_mod, only: verts_rphiz, verts_sthetaphi, nvert

    character(len=100) :: filename_vertex_coordinates
    integer :: vc_unit
    integer :: i

    open(newunit = vc_unit, file = filename_vertex_coordinates)
    101 format(1000(e21.14,x))
    if (coord_system.eq.1) then
        ![R,phi,Z]: Write vertex coordinates to file
        do i=1, nvert
            write(vc_unit,101) verts_rphiz(1, i), verts_rphiz(2, i), verts_rphiz(3, i)
        end do
    elseif (coord_system.eq.2) then
        ![s,theta,phi]: Write vertex coordinates to file
        do i=1, nvert
            write(vc_unit,101) verts_sthetaphi(1, i), verts_sthetaphi(2, i), verts_sthetaphi(3, i)
        end do
    endif
    close(vc_unit)

end subroutine write_vertex_coordinates

subroutine write_prism_volumes(filename_prism_volumes)

    character(len=100) :: filename_prism_volumes
    integer :: pv_unit

    open(newunit = pv_unit, file = filename_prism_volumes)
    write(pv_unit,'(ES20.10E4)') prism_volumes
    close(pv_unit)

end subroutine write_prism_volumes

subroutine write_refined_prism_volumes(filename_refined_prism_volumes)

    character(len=100) :: filename_refined_prism_volumes
    integer :: rpv_unit

    open(newunit = rpv_unit, file = filename_refined_prism_volumes)
    write(rpv_unit,'(ES20.10E4)') refined_prism_volumes
    close(rpv_unit)

end subroutine write_refined_prism_volumes

subroutine write_boltzmann_densities(filename_boltzmann_densities)

    character(len=100) :: filename_boltzmann_densities
    integer :: bd_unit

    open(newunit = bd_unit, file = filename_boltzmann_densities)
    write(bd_unit,'(ES20.10E4)') n_b
    close(bd_unit)

end subroutine write_boltzmann_densities

subroutine write_electric_potential(filename_electric_potential)

    character(len=100) :: filename_electric_potential
    integer :: epv_unit

    open(newunit = epv_unit, file = filename_electric_potential)
    write(epv_unit,'(ES20.10E4)') elec_pot_vec
    close(epv_unit)

end subroutine write_electric_potential

subroutine write_moments(filename_prism_moments, filename_prism_moments_summed_squares, filename_tetr_moments)

    use tetra_grid_mod, only: ntetr

    character(len=100) :: filename_prism_moments, filename_prism_moments_summed_squares, filename_tetr_moments
    integer :: p_moments_unit, pmss_unit, t_moments_unit
    integer :: l, i

    open(newunit = p_moments_unit, file = filename_prism_moments)
    open(newunit = pmss_unit, file = filename_prism_moments_summed_squares)
    open(newunit = t_moments_unit, file = filename_tetr_moments)

    if (n_moments.gt.0) then
        do l = 1,n_prisms
            do i = 1,n_moments - 1
                write(p_moments_unit,'(2ES20.10E4)',advance="no") real(prism_moments(i,l)), aimag(prism_moments(i,l))
            enddo
                write(p_moments_unit,'(2ES20.10E4)') real(prism_moments(n_moments,l)), aimag(prism_moments(n_moments,l))
        enddo
        if (boole_squared_moments) then
            do l = 1,n_prisms
                do i = 1,n_moments - 1
                    write(pmss_unit,'(2ES20.10E4)',advance="no") real(prism_moments_squared(i,l)), &
                                                                    & aimag(prism_moments_squared(i,l))
                enddo
                    write(pmss_unit,'(2ES20.10E4)') real(prism_moments_squared(n_moments,l)), &
                                                       & aimag(prism_moments_squared(n_moments,l))
            enddo
        endif
        do l = 1,ntetr
            do i = 1,n_moments - 1
                write(t_moments_unit,'(2ES20.10E4)',advance="no") real(tetr_moments(i,l)), aimag(tetr_moments(i,l))
            enddo
                write(t_moments_unit,'(2ES20.10E4)') real(tetr_moments(n_moments,l)), aimag(tetr_moments(n_moments,l))
        enddo
    endif

    close(p_moments_unit)
    close(pmss_unit)
    close(t_moments_unit)

end subroutine write_moments

subroutine write_fourier_moments(filename_fourier_moments)

    character(len=100) :: filename_fourier_moments
    integer :: fm_unit, l, i

    open(newunit = fm_unit, file = filename_fourier_moments)
    do l = 1,n_triangles
        do i = 1,n_fourier_modes-1
            write(fm_unit,'(2ES20.10E4)',advance="no") real(moments_in_frequency_space(1,l,i)), &
                                                       aimag(moments_in_frequency_space(1,l,i))
        enddo
            write(fm_unit,'(2ES20.10E4)') real(moments_in_frequency_space(1,l,n_fourier_modes)), &
                                          aimag(moments_in_frequency_space(1,l,n_fourier_modes))
    enddo
    close(fm_unit)

end subroutine write_fourier_moments

subroutine name_files(filenames)

    type(filenames_t), intent(out) :: filenames

    filenames%exit_times = 'exit_times.dat'
    filenames%remaining_particles = 'remaining_particles.dat'
    filenames%poincare_maps = 'poincare_maps.dat'
    filenames%prism_moments = 'prism_moments.dat'
    filenames%prism_moments_summed_squares = 'prism_moments_summed_squares.dat'
    filenames%vertex_coordinates = 'vertex_coordinates.dat'
    filenames%vertex_indices = 'vertex_indices.dat'
    filenames%prism_volumes = 'prism_volumes.dat'
    filenames%fourier_moments = 'fourier_moments.dat'
    filenames%refined_prism_volumes = 'refined_prism_volumes.dat'
    filenames%electric_potential = 'elec_pot_vec.dat'
    filenames%boltzmann_densities = 'boltzmann_densities.dat'
    filenames%divertor_intersections = 'divertor_intersections.dat'
    filenames%tetr_moments = 'tetr_moments.dat'

end subroutine name_files

subroutine unlink_files(filenames)

    type(filenames_t) :: filenames

    call unlink(filenames%exit_times)
    call unlink(filenames%remaining_particles)
    call unlink(filenames%poincare_maps)
    call unlink(filenames%prism_moments)
    call unlink(filenames%prism_moments_summed_squares)
    call unlink(filenames%vertex_coordinates)
    call unlink(filenames%vertex_indices)
    call unlink(filenames%prism_volumes)
    call unlink(filenames%fourier_moments)
    call unlink(filenames%refined_prism_volumes)
    call unlink(filenames%electric_potential)
    call unlink(filenames%boltzmann_densities)
    call unlink(filenames%divertor_intersections)
    call unlink(filenames%tetr_moments)

end subroutine unlink_files

subroutine open_files_before_main_loop

    open(newunit = et_unit, file = 'exit_times.dat')
    open(newunit = rp_unit, file = 'remaining_particles.dat')
    open(newunit = pm_unit, file = 'poincare_maps.dat')
    open(newunit = di_unit, file = 'divertor_intersections.dat')

end subroutine open_files_before_main_loop

subroutine close_files

    close(et_unit)
    close(rp_unit)
    close(pm_unit)
    close(di_unit)

end subroutine close_files

end module boltzmann_mod