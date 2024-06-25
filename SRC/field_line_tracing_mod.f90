module field_line_tracing_mod
    implicit none
!
    private
!
    double precision :: time_step,energy_eV, min_poloidal_flux, max_poloidal_flux
    double precision, dimension(:,:), allocatable :: verts
    integer :: seed_option, num_particles,&
               & ind_a, ind_b, ind_c, lost_outside, lost_inside
    double precision :: n_particles, density, constant_part_of_weights, z_div_plate
    double complex, dimension(:,:), allocatable :: weights
    double precision, dimension(:,:), allocatable :: sqrt_g
    double precision, dimension(:), allocatable :: J_perp, poloidal_flux, temperature_vector
    logical :: boole_linear_density_simulation, boole_linear_temperature_simulation, &
             & boole_poincare_plot, boole_divertor_intersection, boole_collisions, boole_point_source, boole_precalc_collisions, &
             & boole_refined_sqrt_g, boole_boltzmann_energies
    double precision, dimension(:,:,:), allocatable :: randcol
    integer :: randcoli = int(1.0d5)
    integer :: n_poincare_mappings, n_mappings_ignored
    integer :: pm_unit, di_unit, check_unit !file units
    type counter_array
        integer :: lost_particles = 0
        integer :: tetr_pushings = 0
        integer :: phi_0_mappings = 0
        integer :: divertor_intersections = 0
        integer :: survived_particles = 0
        integer :: ignored_particles = 0
    end type counter_array
!
    !Namelist for field_line_tracing input
    NAMELIST /field_line_tracing_nml/ time_step,energy_eV,n_particles,boole_poincare_plot,n_poincare_mappings,n_mappings_ignored, &
    & boole_divertor_intersection, z_div_plate,boole_point_source,boole_collisions, &
    & boole_precalc_collisions,density,boole_refined_sqrt_g,boole_boltzmann_energies, boole_linear_density_simulation, &
    & boole_linear_temperature_simulation,seed_option
!
    public :: calc_field_lines
!
contains
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine field_line_tracing_inp()
!
    open(unit=71, file='field_line_tracing.inp', status='unknown')
    read(71,nml=field_line_tracing_nml)
    close(71)
!
    print *,'GORILLA: Loaded input data from field_line_tracing.inp'
end subroutine field_line_tracing_inp
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine calc_starting_conditions(v0,start_pos_pitch_mat)
!
    use constants, only: pi, ev2erg
    use tetra_grid_mod, only: verts_rphiz, verts_sthetaphi, ntetr
    use find_tetra_mod, only: find_tetra
    use tetra_grid_settings_mod, only: grid_kind
    use tetra_physics_mod, only: coord_system
    use collis_ions, only: collis_init, stost
!
    double precision, intent(in)                                   :: v0
    double precision, dimension(:,:), allocatable, intent(out)     :: start_pos_pitch_mat
    double precision                                               :: rand_scalar, vpar, vperp
    double precision                                               :: amin, cmin, cmax, amax !is set globally
    double precision, dimension(:), allocatable                    :: rand_vector
    double precision, dimension(:), allocatable                    :: rand_matrix2
    integer                                                        :: i
    double precision, dimension(3)                                 :: x
    integer                                                        :: ind_tetr_out,iface
!
!!!!comment out the following section to make starting conditions really random!!!
!
    integer,dimension(:), allocatable                              :: seed
    integer                                                        :: n
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
    allocate(start_pos_pitch_mat(5,num_particles))
    allocate(rand_vector(num_particles))
    allocate(rand_matrix2(num_particles))
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
        start_pos_pitch_mat(ind_a,:) = (/(214 + i*(216-214)/n_particles, i=1,num_particles)/)!r
        start_pos_pitch_mat(ind_b,:) = 0.0d0  !1d-1 !phi in cylindrical and flux coordinates
        start_pos_pitch_mat(ind_c,:) = 12d0 !z in cylindrical, theta in flux coordinates
    endif
!
    start_pos_pitch_mat(4,:) = 1 !delete this once i have a proper subroutine for field line tracing
!
    call RANDOM_NUMBER(rand_matrix2)
    if (boole_boltzmann_energies) then !compare with equation 133 of master thesis of Jonatan Schatzlmayr (remaining parts will be added later)
        start_pos_pitch_mat(5,:) = 5*energy_eV*rand_matrix2(:) !boltzmann energy distribution
        constant_part_of_weights = constant_part_of_weights*10/sqrt(pi)*energy_eV*ev2erg
    endif
!
    weights(:,1) = constant_part_of_weights
    if (boole_refined_sqrt_g.eqv..false.) weights(:,1) = constant_part_of_weights*start_pos_pitch_mat(ind_a,:)
!
    if (boole_precalc_collisions) then
        allocate(randcol(num_particles,randcoli,3))
        call RANDOM_NUMBER(randcol)
        !3.464102 = sqrt(12), this creates a random number with zero average and unit variance
        randcol(:,:,1:2:3) =  3.464102*(randcol(:,:,1:2:3)-.5)
    endif
end subroutine calc_starting_conditions
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine calc_field_lines
!
    use orbit_timestep_gorilla_mod, only: initialize_gorilla
    use constants, only: ev2erg,pi,echarge,ame,amp,clight
    use tetra_physics_mod, only: particle_mass,particle_charge,cm_over_e,mag_axis_R0, coord_system, tetra_physics
    use omp_lib, only: omp_get_thread_num, omp_get_num_threads, omp_set_num_threads
    use parmot_mod, only: rmu,ro0
    use velo_mod, only: isw_field_type
    use supporting_functions_mod, only: theta_sym_flux2theta_vmec,theta_vmec2theta_sym_flux
    use tetra_grid_settings_mod, only: n_field_periods, grid_size
    use tetra_grid_mod, only: ntetr, nvert, verts_rphiz, tetra_grid, verts_sthetaphi
    use gorilla_settings_mod, only: ispecies
    use collis_ions, only: collis_init, stost
    use find_tetra_mod, only: find_tetra
    use gorilla_applets_settings_mod, only: i_option
    use field_mod, only: ipert
!

!
    double precision, dimension(:,:), allocatable :: start_pos_pitch_mat, dens_mat, temp_mat, vpar_mat, efcolf_mat, &
                                                     velrat_mat, enrat_mat, dens_mat_tetr, temp_mat_tetr
    double precision :: v0,pitchpar,vpar,vperp,t_remain,t_confined, v, maxcol
    integer :: kpart,i,j,n,m,k,ind_tetr,iface,ierr,err,num_background_species, inorout
    integer :: n_start, n_end, i_part, count_integration_steps
    double precision, dimension(3) :: x_rand_beg,x,randnum
    logical :: boole_initialized,boole_particle_lost
    double precision :: dtau, dphi,dtaumin, t_step
    double precision, dimension(5) :: z, zet
    double precision :: m0,z0
    double precision, dimension(:), allocatable :: efcolf,velrat,enrat,vpar_background,mass_num,charge_num,dens,temp
    type(counter_array) counter, counter_loop
!
    !Load input for boltzmann computation
    call field_line_tracing_inp()
!
    num_particles = int(n_particles)
    n_start = 1
    n_end = num_particles
!

    open(83, file='field_divB0.inp')
    read(83,*) ipert        ! 0=eq only, 1=vac, 2=vac+plas no derivatives,
    close(83)

    !Initialize GORILLA
    call initialize_gorilla(i_option,ipert)

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! delete this again afterwards !!!!!!!!!!!!!!!!!!!!!!!
    if (ispecies.eq.4) particle_charge = 15*echarge
    print*, 'particle charge number = ', particle_charge/echarge
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    allocate(weights(num_particles,1))
    allocate(J_perp(num_particles))
    allocate(poloidal_flux(num_particles))
    allocate(temperature_vector(num_particles))
!
    poloidal_flux = 0
    temperature_vector = 0
!
    call calc_square_root_g
!
    !Compute velocity module from kinetic energy dependent on particle species
    v0=sqrt(2.d0*energy_eV*ev2erg/particle_mass)
!
print*, 'calc_starting_conditions started'
    call calc_starting_conditions(v0,start_pos_pitch_mat)
print*, 'calc_starting_conditions finished'
!
    max_poloidal_flux = 0
    min_poloidal_flux = tetra_physics(1)%Aphi1
    do i = 1, ntetr
        max_poloidal_flux = max(max_poloidal_flux,tetra_physics(i)%Aphi1 + sum(tetra_physics(i)%gAphi* &
                            & (verts([1,2,3],tetra_grid(i)%ind_knot(4))-verts([1,2,3],tetra_grid(i)%ind_knot(1)))))
        min_poloidal_flux = min(min_poloidal_flux,tetra_physics(i)%Aphi1)
    enddo
!
    !write subroutine collis_precomp
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
        vpar_mat = 0
        m0 = particle_mass/amp
        z0 = particle_charge/echarge
        print*, 'm0 = ', m0, 'z0 = ', z0
!
!!!!!!!!!!!!!!!!!!!! Alternative route is taken because data is not available per vertex but per tetrahedron !!!!!!!!!!!!!!!!!!!!!!!
!
        allocate(dens_mat_tetr(num_background_species-1,ntetr))
        allocate(temp_mat_tetr(num_background_species,ntetr))
!
        open(78, file = 'background/Te_d.dat')
        read(78,'(e16.9)') (temp_mat_tetr(2,i),i=1,ntetr/grid_size(2),3)
        close(78)
!
        open(79, file = 'background/Ti_d.dat')
        read(79,'(e16.9)') (temp_mat_tetr(1,i),i=1,ntetr/grid_size(2),3)
        close(79)
!
        open(80, file = 'background/ne_d.dat')
        read(80,'(e16.9)') (dens_mat_tetr(1,i),i=1,ntetr/grid_size(2),3)
        close(80)
!
        do i = 1,grid_size(2)-1
            temp_mat_tetr(:,i*ntetr/grid_size(2)+1:(i+1)*ntetr/grid_size(2):3) = temp_mat_tetr(:,1:ntetr/grid_size(2):3)
            dens_mat_tetr(:,i*ntetr/grid_size(2)+1:(i+1)*ntetr/grid_size(2):3) = dens_mat_tetr(:,1:ntetr/grid_size(2):3)
        enddo
        do i = 1,2
            temp_mat_tetr(:,1+i:ntetr:3) = temp_mat_tetr(:,1:ntetr:3)
            dens_mat_tetr(:,1+i:ntetr:3) = dens_mat_tetr(:,1:ntetr:3)
        enddo
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
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
!
        kpart = 0
        maxcol = 0
        lost_outside = 0
        lost_inside = 0
        count_integration_steps = 0
!
        call unlink_files
        call open_files
!
        if (boole_collisions) deallocate(efcolf,velrat,enrat)
!
        !$OMP PARALLEL DEFAULT(NONE) &
        !$OMP& SHARED(num_particles,kpart,v0,time_step,boole_collisions, &
        !$OMP& dtau,dtaumin,n_start,n_end, &
        !$OMP& start_pos_pitch_mat,boole_boltzmann_energies,count_integration_steps, &
        !$OMP& density,energy_eV,dens_mat,temp_mat,vpar_mat,tetra_grid,tetra_physics, &
        !$OMP& efcolf_mat,velrat_mat,enrat_mat,num_background_species,randcol,randcoli,maxcol,boole_precalc_collisions, counter) &
        !$OMP& FIRSTPRIVATE(particle_mass, particle_charge) &
        !$OMP& PRIVATE(n,boole_particle_lost,x_rand_beg,x,pitchpar,vpar,vperp,boole_initialized,t_step,err,zet, &
        !$OMP& ind_tetr,iface,t_remain,t_confined,z,ierr,v, &
        !$OMP& i,efcolf,velrat,enrat,vpar_background,inorout,randnum,j,counter_loop)
!
        print*, 'get number of threads', omp_get_num_threads()
        if (boole_collisions) allocate(efcolf(num_background_species),velrat(num_background_species),enrat(num_background_species))
        !$OMP DO
!
        !Loop over particles
        do n = n_start,n_end !1,num_particles
            !if (.not.any(n.eq.(/31997,8046,16148,35518,12921,16318,3807,652,15296,19990,16976,6843,2603/))) cycle
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

            !set counter variables to zero
            call set_counters_zero(counter_loop)

            t_step = time_step
            t_confined = 0
!
            !You need x_rand_beg(1,3), pitchpar(1) (between -1 and 1), energy is already given
            x_rand_beg = start_pos_pitch_mat(1:3,n)
            pitchpar = start_pos_pitch_mat(4,n)
!
            x = x_rand_beg
            vpar = pitchpar * v0
            vperp = sqrt(v0**2-vpar**2)
            if (boole_boltzmann_energies) then
                v = sqrt(start_pos_pitch_mat(5,n)*ev2erg*2/particle_mass)
                vpar = pitchpar * v
                vperp = sqrt(v**2-vpar**2)
            endif
!
            i = 0
            do while (t_confined.lt.time_step)
                i = i+1
                !Orbit integration
                if (i.eq.1) then
                    boole_initialized = .false.
                endif
!
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
!
                call orbit_timestep_gorilla_boltzmann(x,vpar,vperp,t_step,boole_initialized,ind_tetr,iface, &
                                & n,v,start_pos_pitch_mat,inorout,counter_loop,t_remain)
!
                t_confined = t_confined + t_step - t_remain
                !Lost particle handling
                if(ind_tetr.eq.-1) then
!write another if clause (if hole size = minimal .and. particle lost inside .and. boole_cut_out_hole = .true.
!(use extra variable m in orbit routine (0 normally, 1 when lost outside, -1 when lost inside))),
!if m = 1 do as now, if m = -1 select arbitrary newp position and update x, vpar and vperp)
                    write(75,*) t_confined, x, n
                    !print*, t_confined, x
                        counter_loop%lost_particles = 1
                        boole_particle_lost = .true.
                    exit
                endif
!
                v = sqrt(vpar**2+vperp**2)
            enddo
            !$omp critical
            count_integration_steps = count_integration_steps + i
            maxcol = max(dble(i)/dble(randcoli),maxcol)
            call add_counter_loop_to_counter(counter_loop,counter)
            !$omp end critical
        enddo !n
        !$OMP END DO
        !$OMP END PARALLEL
!
        call close_files
!
if (boole_precalc_collisions) print*, "maxcol = ", maxcol
print*, 'Number of lost particles',counter%lost_particles
print*, 'average number of pushings = ', counter%tetr_pushings/n_particles
print*, 'average number of toroidal revolutions = ', counter%phi_0_mappings/n_particles
print*, 'number of particles having hit the divertor = ', counter%divertor_intersections
print*, 'number of particles having hit the z-value of the divertor plate within the first 10 revolutions = ', &
         & counter%ignored_particles
print*, 'number of particles still in orbit after tracing time has passed = ', counter%survived_particles
print*, 'tracing time in seconds = ', time_step
!
deallocate(start_pos_pitch_mat)
!
end subroutine calc_field_lines
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine orbit_timestep_gorilla_boltzmann(x,vpar,vperp,t_step,boole_initialized,ind_tetr,iface, n,v,start_pos_pitch_mat, &
                                          & inorout, counter_loop,t_remain_out)
!
    use pusher_tetra_rk_mod, only: pusher_tetra_rk,initialize_const_motion_rk
    use pusher_tetra_poly_mod, only: pusher_tetra_poly,initialize_const_motion_poly
    use pusher_tetra_field_lines_mod, only: pusher_tetra_field_lines
    use tetra_physics_poly_precomp_mod , only: make_precomp_poly_perpinv, initialize_boole_precomp_poly_perpinv, &
        & alloc_precomp_poly_perpinv
    use tetra_physics_mod, only: tetra_physics,particle_charge,particle_mass,cm_over_e
    use gorilla_settings_mod, only: poly_order
    use orbit_timestep_gorilla_mod, only: check_coordinate_domain
    use supporting_functions_mod, only: bmod_func, vperp_func
    use find_tetra_mod, only: find_tetra
    use constants, only: pi, ev2erg
    use tetra_grid_mod, only: tetra_grid, ntetr
!

!
    double precision, dimension(3), intent(inout)   :: x
    double precision, intent(inout)                 :: vpar,vperp
    double precision, intent(in)                    :: t_step
    logical, intent(inout)                          :: boole_initialized
    integer, intent(inout)                          :: ind_tetr,iface
    double precision, intent(out), optional         :: t_remain_out
    double precision, dimension(3)                  :: z_save, x_save
    double precision                                :: vperp2,t_remain,t_pass,vpar_save, v, aphi
    logical                                         :: boole_t_finished
    integer                                         :: ind_tetr_save,iper_phi,k, m, n, inorout
    double precision                                :: perpinv,speed, r, z, phi, B, phi_elec_func, wlin
    double precision, dimension(:,:)                :: start_pos_pitch_mat
    double precision, dimension(3)                  :: xyz, xyz_save
    type(counter_array), intent(inout)              :: counter_loop

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
            t_remain_out = t_step
            return
        endif
!
        r = x(1) - verts(1, tetra_grid(ind_tetr)%ind_knot(1))
        z = x(3) - verts(3,tetra_grid(ind_tetr)%ind_knot(1))
        phi = x(2) - verts(2,tetra_grid(ind_tetr)%ind_knot(1))
!
        if (boole_refined_sqrt_g) then
            weights(n,1) = weights(n,1)* &
                                    &  (sqrt_g(ind_tetr,1)+r*sqrt_g(ind_tetr,2)+z*sqrt_g(ind_tetr,3))/ &
                                    &  (sqrt_g(ind_tetr,4)+r*sqrt_g(ind_tetr,5)+z*sqrt_g(ind_tetr,6))
        endif
!
        if (boole_linear_density_simulation.or.boole_linear_temperature_simulation) then
            poloidal_flux(n) = tetra_physics(ind_tetr)%Aphi1 + sum(tetra_physics(ind_tetr)%gAphi*(/r,phi,z/)) + &
                                             & cm_over_e*vpar*&
                                             & (tetra_physics(ind_tetr)%h2_1+sum(tetra_physics(ind_tetr)%gh2*(/r,phi,z/)))
        endif
        if (boole_linear_density_simulation) then
            weights(n,1) = weights(n,1)*(max_poloidal_flux*1.1-poloidal_flux(n))/(max_poloidal_flux*1.1)
        endif
!
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
!
        !compute J_perp for perpendicular pressure
        B = tetra_physics(ind_tetr)%bmod1+sum(tetra_physics(ind_tetr)%gb*(/r,phi,z/))
        J_perp(n) = particle_mass*v**2*(1-start_pos_pitch_mat(4,n)**2)*cm_over_e/(2*B)*(-1)
        !-1 because of negative gyrophase
!
        boole_initialized = .true.
    endif
!
    !Exit the subroutine after initialization, if time step equals zero
    if(t_step.eq.0.d0) return
!
    inorout = 0
!
    !Squared perpendicular velocity
    vperp2 = vperp**2
!
    !Compute relative particle position
    z_save = x-tetra_physics(ind_tetr)%x1
!
    !Compute perpendicular invariant of particle
    perpinv=-0.5d0*vperp2/bmod_func(z_save,ind_tetr)
!
    !Integrate particle orbit for given time step
    t_remain = t_step
!
    !Logical for handling time integration
    boole_t_finished = .false.
!
    counter_loop%tetr_pushings = counter_loop%tetr_pushings - 1 !set tetr_pushings to -1 because when entering the loop
    !it wil go back to one without pushing

    !Loop for tetrahedron pushings until t_step is reached
    do
!
        counter_loop%tetr_pushings = counter_loop%tetr_pushings + 1
        !Domain Boundary
        if(ind_tetr.eq.-1) then
            aphi = tetra_physics(ind_tetr_save)%Aphi1+sum(tetra_physics(ind_tetr_save)%gAphi*z_save)
            print *, 'WARNING: Particle lost.'
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
!
        ind_tetr_save = ind_tetr
        vpar_save = vpar
        x_save = x
!
        !Calculate trajectory
        call pusher_tetra_field_lines(ind_tetr,iface,x,vpar,z_save,t_remain,t_pass,boole_t_finished,iper_phi)
        if (boole_t_finished) counter_loop%survived_particles = 1

        t_remain = t_remain - t_pass
        if (iper_phi.ne.0) then
            counter_loop%phi_0_mappings = counter_loop%phi_0_mappings + 1!iper_phi
            ! if ((counter_loop%phi_0_mappings.gt.n_mappings_ignored).and. &
            !   & (any(n.eq.(/31997,8046,16148,35518,12921,16318,3807,652,15296,19990,16976,6843,2603/)))) then
            !     !$omp critical
            !         write(check_unit,*) n, x
            !     !$omp end critical
            ! endif
            if ((boole_poincare_plot).and.(counter_loop%phi_0_mappings.gt.n_mappings_ignored)) then
                !$omp critical
                    write(pm_unit,*) x
                !$omp end critical
            endif
            if ((boole_poincare_plot).and.(counter_loop%phi_0_mappings.eq.n_poincare_mappings)) then
                boole_t_finished = .true.
            endif
        endif
        if ((boole_divertor_intersection).and.(x(3).lt.z_div_plate)) then
            boole_t_finished = .true.
           if ((counter_loop%phi_0_mappings.gt.n_mappings_ignored)) then
            call calc_plane_intersection(x_save,x,z_div_plate)

            !$omp critical
                write(di_unit,*) x, n
            !$omp end critical
                counter_loop%divertor_intersections = 1
           else
            counter_loop%ignored_particles  = 1
           endif
        endif
!
        !Orbit stops within cell, because "flight"-time t_step has finished
        if(boole_t_finished) then
            if(present(t_remain_out)) then
                t_remain_out = 0
            endif
            if (x(3).gt.z_div_plate) print*, n, x, ind_tetr
            exit
        endif
!
    enddo !Loop for tetrahedron pushings
!
    !Compute vperp from position
    vperp = vperp_func(z_save,perpinv,ind_tetr_save)
!
end subroutine orbit_timestep_gorilla_boltzmann
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine calc_plane_intersection(x_save,x,z_plane)
!
    use constants, only : pi
!
    double precision, dimension(3), intent(in) :: x_save
    double precision, dimension(3), intent(inout) :: x
    double precision, intent(in) :: z_plane
    double precision :: rel_dist_z
!
    rel_dist_z = (z_plane-x_save(3))/(x(3)-x_save(3))
    x(1) = x_save(1) + rel_dist_z*(x(1)-x_save(1))
    if (abs(x(2)-x_save(2)).gt.pi) then
        x(2) = modulo(x_save(2) + 2*pi-abs(x(2)-x_save(2)),2*pi)
    else
        x(2) = x_save(2) + rel_dist_z*(x(2)-x_save(2))
    endif
    x(3) = z_plane
!
end subroutine calc_plane_intersection
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine cyl_to_cart(xcyl,xcart)
!
    double precision, dimension(3), intent(in) :: xcyl
    double precision, dimension(3), intent(out) :: xcart
!
    xcart(1) = xcyl(1)*cos(xcyl(2))
    xcart(2) = xcyl(1)*sin(xcyl(2))
    xcart(3) = xcyl(3)
!
end subroutine cyl_to_cart
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine cart_to_cyl(xcart,xcyl)
!
    double precision, dimension(3), intent(in) :: xcart
    double precision, dimension(3), intent(out) :: xcyl
!
    xcyl(1) = hypot(xcart(1),xcart(2))
    xcyl(2) = atan2(xcart(2),xcart(1))
    xcyl(3) = xcart(3)
end subroutine cart_to_cyl
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine set_counters_zero(counter_loop)
!
    type(counter_array), intent(inout) :: counter_loop

    counter_loop%lost_particles = 0
    counter_loop%tetr_pushings = 0
    counter_loop%phi_0_mappings = 0
    counter_loop%divertor_intersections = 0
    counter_loop%survived_particles = 0
    counter_loop%ignored_particles = 0

end subroutine set_counters_zero
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine add_counter_loop_to_counter(counter_loop,counter)
!
        type(counter_array), intent(in) :: counter_loop
        type(counter_array), intent(inout) :: counter
!
        counter%lost_particles = counter%lost_particles + counter_loop%lost_particles
        counter%tetr_pushings = counter%tetr_pushings + counter_loop%tetr_pushings
        counter%phi_0_mappings = counter%phi_0_mappings + counter_loop%phi_0_mappings
        counter%divertor_intersections = counter%divertor_intersections + counter_loop%divertor_intersections
        counter%survived_particles = counter%survived_particles + counter_loop%survived_particles
        counter%ignored_particles = counter%ignored_particles + counter_loop%ignored_particles
!
end subroutine add_counter_loop_to_counter
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine calc_square_root_g
!
    use tetra_physics_mod, only: tetra_physics, hamiltonian_time
    use tetra_grid_mod, only: ntetr, tetra_grid, verts_rphiz
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
subroutine unlink_files
!
    call unlink('poincare_maps.dat')
    call unlink('divertor_intersections.dat')
    call unlink('check.dat')
!
end subroutine unlink_files
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine open_files
!
    open(newunit = pm_unit, file = 'poincare_maps.dat')
    open(newunit = di_unit, file = 'divertor_intersections.dat')
    open(newunit = check_unit, file = 'check.dat')
!
end subroutine open_files
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine close_files
!
    close(pm_unit)
    close(di_unit)
    close(check_unit)
!
end subroutine close_files
!
end module field_line_tracing_mod