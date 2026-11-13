module field_line_tracing_mod

    use, intrinsic :: iso_fortran_env, only: int64
    use volume_integrals_and_sqrt_g_mod, only: sqrt_g

    implicit none

    private

    double precision :: time_step,energy_eV, min_poloidal_flux, max_poloidal_flux
    double precision, dimension(:,:), allocatable :: verts
    integer :: seed_option, num_particles,&
               & ind_a, ind_b, ind_c
    double precision :: n_particles, constant_part_of_weights, z_div_plate
    double complex, dimension(:,:), allocatable :: weights
    double precision, dimension(:), allocatable :: J_perp, poloidal_flux, temperature_vector
    logical :: boole_linear_density_simulation, boole_linear_temperature_simulation, &
             & boole_poincare_plot, boole_divertor_intersection, boole_point_source, &
             & boole_refined_sqrt_g, boole_monoenergetic
    integer :: n_poincare_mappings, n_mappings_ignored
    integer :: pm_unit, di_unit !file units
    !Start points from file, traced in both directions along B (see read_start_points)
    logical :: boole_start_points_file = .false.
    character(len=1024) :: filename_start_points = 'field_line_start_points.dat'
    character(len=1024) :: filename_exit_points = 'field_line_exits.dat'
    double precision :: max_arc_length = 1.0d4 !maximum arc length per direction in m (start-file mode only)
    logical :: boole_write_trajectories = .false.
    character(len=1024) :: filename_trajectories = 'field_line_trajectories.dat'
    integer, parameter :: n_traj_max = 1000
    integer, dimension(n_traj_max) :: traj_start_indices = 0 !start-point indices for which trajectories are written
    integer :: n_start_points = 0, traj_unit
    double precision, dimension(:,:), allocatable :: start_points_rphiz
    type counter_t
        integer :: lost_particles = 0
        integer(int64) :: tetr_pushings = 0
        integer :: phi_0_mappings = 0
        integer :: divertor_intersections = 0
        integer :: survived_particles = 0
        integer :: ignored_particles = 0
        integer(int64) :: integration_steps = 0
    end type counter_t

    !Namelist for field_line_tracing input
    NAMELIST /field_line_tracing_nml/ time_step,energy_eV,n_particles,boole_poincare_plot,n_poincare_mappings,n_mappings_ignored, &
    & boole_divertor_intersection, z_div_plate,boole_point_source, &
    & boole_refined_sqrt_g,boole_monoenergetic, boole_linear_density_simulation, &
    & boole_linear_temperature_simulation,seed_option, &
    & boole_start_points_file, filename_start_points, filename_exit_points, max_arc_length, &
    & boole_write_trajectories, filename_trajectories, traj_start_indices

    public :: calc_field_lines

contains

subroutine calc_field_lines

    use orbit_timestep_gorilla_mod, only: initialize_gorilla
    use constants, only: ev2erg,pi,echarge
    use tetra_physics_mod, only: particle_mass,particle_charge,cm_over_e,mag_axis_R0, coord_system, tetra_physics
    use omp_lib, only: omp_get_thread_num, omp_get_num_threads, omp_set_num_threads
    use parmot_mod, only: rmu,ro0
    use velo_mod, only: isw_field_type
    use supporting_functions_mod, only: theta_sym_flux2theta_vmec,theta_vmec2theta_sym_flux
    use tetra_grid_settings_mod, only: n_field_periods, grid_size
    use tetra_grid_mod, only: ntetr, nvert, verts_rphiz, tetra_grid, verts_sthetaphi
    use gorilla_settings_mod, only: ispecies
    use find_tetra_mod, only: find_tetra
    use gorilla_applets_settings_mod, only: i_option
    use field_mod, only: ipert
    use volume_integrals_and_sqrt_g_mod, only: calc_square_root_g
    use utils_data_pre_and_post_processing_mod, only: get_ipert

    double precision, dimension(:,:), allocatable :: start_pos_pitch_mat
    double precision :: v0,pitchpar,vpar,vperp,t_remain,t_confined, v
    integer :: kpart,i,n,ind_tetr,iface,ierr
    integer :: n_start, n_end
    double precision, dimension(3) :: x_rand_beg,x
    logical :: boole_initialized,boole_particle_lost
    double precision :: t_step
    type(counter_t) :: counter, local_counter
    !start-file mode: per (start point, direction) results
    integer, dimension(:), allocatable :: exit_status
    integer(int64), dimension(:), allocatable :: exit_pushings
    double precision, dimension(:), allocatable :: exit_s_geom, exit_s_time
    double precision, dimension(:,:), allocatable :: exit_x
    double precision :: s_geom
    integer :: ex_unit
    logical :: boole_traj_line, boole_stuck

    !Load input for boltzmann computation
    call read_field_line_tracing_inp()

    num_particles = int(n_particles)
    if (boole_start_points_file) then
        call read_start_points()
        num_particles = 2*n_start_points !each start point is traced in both directions (+1 and -1 along B)
        if (.not. boole_monoenergetic) then
            print *, 'start-points file mode: boole_monoenergetic is set to .true.'
            boole_monoenergetic = .true.
        endif
        allocate(exit_status(num_particles), exit_pushings(num_particles), exit_s_geom(num_particles), &
               & exit_s_time(num_particles), exit_x(3,num_particles))
        exit_status = 0
        exit_pushings = 0
        exit_s_geom = 0.d0
        exit_s_time = 0.d0
        exit_x = 0.d0
    endif
    n_start = 1
    n_end = num_particles

    call get_ipert()

    !Initialize GORILLA
    call initialize_gorilla(i_option,ipert)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! delete this again afterwards !!!!!!!!!!!!!!!!!!!!!!!
    if (ispecies.eq.4) particle_charge = 15*echarge
    print*, 'particle charge number = ', particle_charge/echarge
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    allocate(weights(num_particles,1))
    allocate(J_perp(num_particles))
    allocate(poloidal_flux(num_particles))
    allocate(temperature_vector(num_particles))

    poloidal_flux = 0
    temperature_vector = 0

    call calc_square_root_g

    !Compute velocity module from kinetic energy dependent on particle species
    v0=sqrt(2.d0*energy_eV*ev2erg/particle_mass)

    !Start-file mode: tracing time from maximum arc length. Field lines are traced with |vpar| = v0,
    !so the arc length along B is s = v0*t (cm), see orbit_timestep_gorilla_field_lines.
    if (boole_start_points_file) then
        time_step = max_arc_length*1.d2/v0
        print *, 'start-points file mode: ', n_start_points, ' start points, max. arc length (m) = ', max_arc_length, &
               & ', tracing time (s) = ', time_step
    endif

print*, 'calc_starting_conditions started'
    call calc_starting_conditions(v0,start_pos_pitch_mat)
print*, 'calc_starting_conditions finished'

    max_poloidal_flux = 0
    min_poloidal_flux = tetra_physics(1)%Aphi1
    do i = 1, ntetr
        max_poloidal_flux = max(max_poloidal_flux,tetra_physics(i)%Aphi1 + sum(tetra_physics(i)%gAphi* &
                            & (verts([1,2,3],tetra_grid(i)%ind_knot(4))-verts([1,2,3],tetra_grid(i)%ind_knot(1)))))
        min_poloidal_flux = min(min_poloidal_flux,tetra_physics(i)%Aphi1)
    enddo

        kpart = 0

        call unlink_files
        call open_files
        if (boole_start_points_file .and. boole_write_trajectories) then
            open(newunit = traj_unit, file = trim(filename_trajectories), status = 'replace')
            write(traj_unit,'(a)') '# start_index direction s(m) R(cm) phi(rad) Z(cm) |B|(G)'
        endif

        !$OMP PARALLEL DEFAULT(NONE) &
        !$OMP& SHARED(num_particles,kpart,v0,time_step, &
        !$OMP& n_start,n_end, &
        !$OMP& start_pos_pitch_mat,boole_monoenergetic, &
        !$OMP& energy_eV,tetra_grid,tetra_physics, counter, &
        !$OMP& exit_status,exit_pushings,exit_s_geom,exit_s_time,exit_x, &
        !$OMP& boole_start_points_file,boole_write_trajectories,traj_start_indices) &
        !$OMP& FIRSTPRIVATE(particle_mass, particle_charge) &
        !$OMP& PRIVATE(n,boole_particle_lost,x_rand_beg,x,pitchpar,vpar,vperp,boole_initialized,t_step, &
        !$OMP& ind_tetr,iface,t_remain,t_confined,ierr,v, &
        !$OMP& i,local_counter,s_geom,boole_traj_line,boole_stuck)

        print*, 'get number of threads', omp_get_num_threads()
        !$OMP DO

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
            call set_counters_zero(local_counter)

            t_step = time_step
            t_confined = 0

            !You need x_rand_beg(1,3), pitchpar(1) (between -1 and 1), energy is already given
            x_rand_beg = start_pos_pitch_mat(1:3,n)
            pitchpar = start_pos_pitch_mat(4,n)

            x = x_rand_beg
            v = v0
            vpar = pitchpar * v0
            vperp = sqrt(v0**2-vpar**2)
            if (.not. boole_monoenergetic) then
                v = sqrt(start_pos_pitch_mat(5,n)*ev2erg*2/particle_mass)
                vpar = pitchpar * v
                vperp = sqrt(v**2-vpar**2)
            endif

            s_geom = 0.d0
            boole_traj_line = .false.
            if (boole_start_points_file .and. boole_write_trajectories) &
                & boole_traj_line = any(traj_start_indices.eq.(n+1)/2)

            i = 0
            do while (t_confined.lt.time_step)
                i = i+1
                !Orbit integration
                if (i.eq.1) then
                    boole_initialized = .false.
                endif

                call orbit_timestep_gorilla_field_lines(x,vpar,vperp,t_step,boole_initialized,ind_tetr,iface, &
                                & n,v,start_pos_pitch_mat,local_counter,t_remain,s_geom,boole_traj_line,boole_stuck)

                t_confined = t_confined + t_step - t_remain
                !Start-file mode: record exit (status 1), start outside the mesh (status -1)
                if (boole_start_points_file) then
                    exit_s_geom(n) = s_geom
                    exit_s_time(n) = abs(vpar)*t_confined
                    exit_x(:,n) = x
                    exit_pushings(n) = local_counter%tetr_pushings
                    if (boole_stuck) exit_status(n) = -2
                    if (ind_tetr.eq.-1) then
                        exit_status(n) = 1
                        if (i.eq.1 .and. local_counter%tetr_pushings.le.0) exit_status(n) = -1
                    endif
                endif
                if (boole_stuck) then
                    print *, 'warning: field line ', n, ' does not advance, tracing stopped at x = ', x
                    exit
                endif
                !Lost particle handling
                if(ind_tetr.eq.-1) then
!write another if clause (if hole size = minimal .and. particle lost inside .and. boole_cut_out_hole = .true.
!(use extra variable m in orbit routine (0 normally, 1 when lost outside, -1 when lost inside))),
!if m = 1 do as now, if m = -1 select arbitrary newp position and update x, vpar and vperp)
                    write(75,*) t_confined, x, n
                    !print*, t_confined, x
                        local_counter%lost_particles = 1
                        boole_particle_lost = .true.
                    exit
                endif

                v = sqrt(vpar**2+vperp**2)
            enddo
            !$omp critical
            counter%integration_steps = counter%integration_steps + i
            call add_local_counter_to_counter(local_counter,counter)
            !$omp end critical
        enddo !n
        !$OMP END DO
        !$OMP END PARALLEL

        call close_files

    if (boole_start_points_file) then
        if (boole_write_trajectories) close(traj_unit)
        open(newunit = ex_unit, file = trim(filename_exit_points), status = 'replace')
        write(ex_unit,'(a)') '# start_index direction(+1 along B,-1 against B) status(1 left mesh,0 not left within '// &
            & 'max_arc_length,-1 start point not in mesh,-2 stopped: no progress) s_geom(m) s_vt(m) R_end(cm) phi_end(rad) Z_end(cm) n_pushings'
        do n = 1,num_particles
            write(ex_unit,'(i8,1x,i3,1x,i3,2(1x,es16.8),3(1x,es16.8),1x,i12)') (n+1)/2, merge(1,-1,mod(n,2).eq.1), &
                & exit_status(n), exit_s_geom(n)*1.d-2, exit_s_time(n)*1.d-2, exit_x(:,n), exit_pushings(n)
        enddo
        close(ex_unit)
        print *, 'start-points file mode: lines leaving the mesh = ', count(exit_status.eq.1), &
               & ', not leaving within max_arc_length = ', count(exit_status.eq.0), &
               & ', start point not in mesh = ', count(exit_status.eq.-1), &
               & ', stopped (no progress) = ', count(exit_status.eq.-2)
        print *, 'exit points written to ', trim(filename_exit_points)
    endif

print*, 'average number of pushings = ', counter%tetr_pushings/n_particles
print*, 'average number of toroidal revolutions = ', counter%phi_0_mappings/n_particles
if (boole_divertor_intersection) then
    print*, 'number of particles having hit the divertor = ', counter%divertor_intersections
    print*, 'number of particles having hit the z-value of the divertor plate within the first 10 revolutions = ', &
            & counter%ignored_particles
endif
print*, 'tracing time in seconds = ', time_step
print*, 'Number of lost particles',counter%lost_particles
print*, 'number of particles still in orbit after tracing time has passed = ', counter%survived_particles

deallocate(start_pos_pitch_mat)

end subroutine calc_field_lines

subroutine orbit_timestep_gorilla_field_lines(x,vpar,vperp,t_step,boole_initialized,ind_tetr,iface, n,v,start_pos_pitch_mat, &
                                          & local_counter,t_remain_out,arc_length,boole_traj,boole_stuck)

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
    use marker_distribution_mod, only: linear_profile_factor
    use tetra_grid_mod, only: tetra_grid, ntetr
    use tetra_grid_settings_mod, only: n_field_periods

    double precision, dimension(3), intent(inout)   :: x
    double precision, intent(inout)                 :: vpar,vperp
    double precision, intent(in)                    :: t_step
    logical, intent(inout)                          :: boole_initialized
    integer, intent(inout)                          :: ind_tetr,iface
    double precision, intent(out), optional         :: t_remain_out
    double precision, dimension(3)                  :: z_save, x_save
    double precision                                :: vperp2,t_remain,t_pass,vpar_save, v
    logical                                         :: boole_t_finished
    integer                                         :: ind_tetr_save,iper_phi,k, m, n
    double precision                                :: perpinv,speed, r, z, phi, B, phi_elec_func, wlin
    double precision, dimension(:,:)                :: start_pos_pitch_mat
    double precision, dimension(3)                  :: xyz, xyz_save
    type(counter_t), intent(inout)                  :: local_counter
    double precision, intent(inout), optional       :: arc_length !geometric arc length along the field line (cm)
    logical, intent(in), optional                   :: boole_traj   !write trajectory (s,R,phi,Z,|B|) to traj_unit
    logical, intent(out), optional                  :: boole_stuck  !no progress in many consecutive pushings
    logical                                         :: boole_arc, boole_write_traj
    integer                                         :: n_no_progress
    integer, parameter                              :: n_no_progress_max = 1000
    double precision                                :: dphi, phi_period, r_mid


    if (present(boole_stuck)) boole_stuck = .false.

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
            weights(n,1) = weights(n,1)*linear_profile_factor(poloidal_flux(n)/max_poloidal_flux)
        endif

        if (.not. boole_monoenergetic) then
            !compare with equation 133 of master thesis of Jonatan Schatzlmayr (remaining parts have been added before)
            phi_elec_func = tetra_physics(ind_tetr)%Phi1 + sum(tetra_physics(ind_tetr)%gPhi*(/r,phi,z/))
            if (.not. boole_linear_temperature_simulation) then
                weights(n,1) = weights(n,1)*sqrt(start_pos_pitch_mat(5,n)*ev2erg)/(energy_eV*ev2erg)**1.5* &
                            & exp(-(start_pos_pitch_mat(5,n)*ev2erg+particle_charge*phi_elec_func)/(energy_eV*ev2erg))
            else
                temperature_vector(n) = energy_eV*ev2erg*linear_profile_factor(poloidal_flux(n)/max_poloidal_flux)
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

    n_no_progress = 0
    boole_arc = present(arc_length)
    boole_write_traj = .false.
    if (present(boole_traj)) boole_write_traj = boole_traj .and. boole_arc
    phi_period = 2.d0*pi/dble(n_field_periods)
    if (boole_write_traj) then
        !$omp critical
        write(traj_unit,'(i8,1x,i3,5(1x,es16.8))') (n+1)/2, int(sign(1.d0,vpar)), arc_length*1.d-2, x, &
            & bmod_func(x-tetra_physics(ind_tetr)%x1,ind_tetr)
        !$omp end critical
    endif

    !Squared perpendicular velocity
    vperp2 = vperp**2

    !Compute relative particle position
    z_save = x-tetra_physics(ind_tetr)%x1

    !Compute perpendicular invariant of particle
    perpinv=-0.5d0*vperp2/bmod_func(z_save,ind_tetr)

    !Integrate particle orbit for given time step
    t_remain = t_step

    !Logical for handling time integration
    boole_t_finished = .false.

    local_counter%tetr_pushings = local_counter%tetr_pushings - 1 !set tetr_pushings to -1 because when entering the loop
    !it wil go back to one without pushing

    !Loop for tetrahedron pushings until t_step is reached
    do

        local_counter%tetr_pushings = local_counter%tetr_pushings + 1
        !Domain Boundary
        if(ind_tetr.eq.-1) then
            !print *, 'WARNING: Particle lost.'
            if( present(t_remain_out)) then
                t_remain_out = t_remain
            endif
            exit
        endif

        ind_tetr_save = ind_tetr
        vpar_save = vpar
        x_save = x

        !Calculate trajectory
        call pusher_tetra_field_lines(ind_tetr,iface,x,vpar,z_save,t_remain,t_pass,boole_t_finished,iper_phi)

        t_remain = t_remain - t_pass

        !Round-off safety: the pusher may report a passing time slightly beyond the remaining time
        if (t_remain.lt.0.d0) then
            t_remain = 0.d0
            boole_t_finished = .true.
        endif

        !Guard against a field line that does not advance (e.g. ping-pong between tetrahedra at an edge)
        if (t_pass.eq.0.d0 .and. .not.boole_t_finished) then
            n_no_progress = n_no_progress + 1
            if (n_no_progress.gt.n_no_progress_max .and. present(boole_stuck)) then
                boole_stuck = .true.
                if(present(t_remain_out)) t_remain_out = t_remain
                exit
            endif
        else
            n_no_progress = 0
        endif

        !Geometric arc length: inside a tetrahedron the field line is a straight line in (R,phi,Z),
        !its length is evaluated with the mid-point radius (cylindrical coordinates only)
        if (boole_arc) then
            dphi = x(2)-x_save(2)
            dphi = dphi - phi_period*anint(dphi/phi_period) !undo periodic relocation
            r_mid = 0.5d0*(x(1)+x_save(1))
            arc_length = arc_length + sqrt((x(1)-x_save(1))**2 + (r_mid*dphi)**2 + (x(3)-x_save(3))**2)
            if (boole_write_traj) then
                !$omp critical
                write(traj_unit,'(i8,1x,i3,5(1x,es16.8))') (n+1)/2, int(sign(1.d0,vpar)), arc_length*1.d-2, x, &
                    & bmod_func(x-tetra_physics(ind_tetr_save)%x1,ind_tetr_save)
                !$omp end critical
            endif
        endif

        if (iper_phi.ne.0) then
            local_counter%phi_0_mappings = local_counter%phi_0_mappings + 1!iper_phi
            if ((boole_poincare_plot).and.(local_counter%phi_0_mappings.gt.n_mappings_ignored) &
                                     .and.(local_counter%phi_0_mappings.le.n_poincare_mappings)) then
                !$omp critical
                    write(pm_unit,*) x
                !$omp end critical
            endif
            if ((boole_poincare_plot).and.(.not.boole_divertor_intersection) &
                 .and.(local_counter%phi_0_mappings.eq.n_poincare_mappings)) then
                boole_t_finished = .true.
            endif
        endif

        if (boole_t_finished) local_counter%survived_particles = 1

        if ((boole_divertor_intersection).and.(x(3).lt.z_div_plate)) then
            boole_t_finished = .true.
           if (local_counter%phi_0_mappings.gt.n_mappings_ignored) then
            call calc_plane_intersection(x_save,x,z_div_plate)

            !$omp critical
                write(di_unit,*) x, n
            !$omp end critical
                local_counter%divertor_intersections = 1
           else
            local_counter%ignored_particles  = 1
           endif
        endif

        !Orbit stops within cell, because "flight"-time t_step has finished
        if(boole_t_finished) then
            if(present(t_remain_out)) then
                t_remain_out = 0
            endif
            !if (x(3).gt.z_div_plate) print*, n, x, ind_tetr
            exit
        endif

    enddo !Loop for tetrahedron pushings

    !Compute vperp from position
    vperp = vperp_func(z_save,perpinv,ind_tetr_save)

end subroutine orbit_timestep_gorilla_field_lines

subroutine calc_plane_intersection(x_save,x,z_plane)

    use constants, only : pi

    double precision, dimension(3), intent(in) :: x_save
    double precision, dimension(3), intent(inout) :: x
    double precision, intent(in) :: z_plane
    double precision :: rel_dist_z

    rel_dist_z = (z_plane-x_save(3))/(x(3)-x_save(3))
    x(1) = x_save(1) + rel_dist_z*(x(1)-x_save(1))
    if (abs(x(2)-x_save(2)).gt.pi) then
        x(2) = modulo(x_save(2) + 2*pi-abs(x(2)-x_save(2)),2*pi)
    else
        x(2) = x_save(2) + rel_dist_z*(x(2)-x_save(2))
    endif
    x(3) = z_plane

end subroutine calc_plane_intersection

subroutine set_counters_zero(local_counter)

    type(counter_t), intent(inout) :: local_counter

    local_counter%lost_particles = 0
    local_counter%tetr_pushings = 0
    local_counter%phi_0_mappings = 0
    local_counter%divertor_intersections = 0
    local_counter%survived_particles = 0
    local_counter%ignored_particles = 0

end subroutine set_counters_zero

subroutine add_local_counter_to_counter(local_counter,counter)

        type(counter_t), intent(in) :: local_counter
        type(counter_t), intent(inout) :: counter

        counter%lost_particles = counter%lost_particles + local_counter%lost_particles
        counter%tetr_pushings = counter%tetr_pushings + local_counter%tetr_pushings
        counter%phi_0_mappings = counter%phi_0_mappings + local_counter%phi_0_mappings
        counter%divertor_intersections = counter%divertor_intersections + local_counter%divertor_intersections
        counter%survived_particles = counter%survived_particles + local_counter%survived_particles
        counter%ignored_particles = counter%ignored_particles + local_counter%ignored_particles

end subroutine add_local_counter_to_counter

subroutine unlink_files

    call unlink('poincare_maps.dat')
    call unlink('divertor_intersections.dat')
    call unlink('check.dat')

end subroutine unlink_files

subroutine open_files

    open(newunit = pm_unit, file = 'poincare_maps.dat')
    open(newunit = di_unit, file = 'divertor_intersections.dat')

end subroutine open_files

subroutine close_files

    close(pm_unit)
    close(di_unit)

end subroutine close_files

subroutine read_field_line_tracing_inp()

    integer :: flt_inp_unit

    open(newunit = flt_inp_unit, file='field_line_tracing.inp', status='unknown')
    read(flt_inp_unit,nml=field_line_tracing_nml)
    close(flt_inp_unit)

    print *,'GORILLA: Loaded input data from field_line_tracing.inp'
end subroutine read_field_line_tracing_inp

subroutine read_start_points()

    !Reads start points (R in cm, phi in rad, Z in cm; one per line, '#' or '!' lines are comments)

    integer :: sp_unit, ios, k
    character(len=1024) :: line
    double precision, dimension(3) :: xx

    open(newunit = sp_unit, file = trim(filename_start_points), status = 'old', action = 'read')
    do k = 1,2
        n_start_points = 0
        rewind(sp_unit)
        do
            read(sp_unit,'(a)',iostat=ios) line
            if (ios.ne.0) exit
            line = adjustl(line)
            if (len_trim(line).eq.0) cycle
            if (line(1:1).eq.'#' .or. line(1:1).eq.'!') cycle
            read(line,*,iostat=ios) xx
            if (ios.ne.0) then
                print *, 'error: cannot read start point from line: ', trim(line)
                stop 1
            endif
            n_start_points = n_start_points + 1
            if (k.eq.2) start_points_rphiz(:,n_start_points) = xx
        enddo
        if (k.eq.1) allocate(start_points_rphiz(3,n_start_points))
    enddo
    close(sp_unit)
    if (n_start_points.eq.0) then
        print *, 'error: no start points in ', trim(filename_start_points)
        stop 1
    endif
    print *, 'read ', n_start_points, ' start points from ', trim(filename_start_points)

end subroutine read_start_points

subroutine calc_starting_conditions(v0,start_pos_pitch_mat)

    use constants, only: pi, ev2erg
    use tetra_grid_mod, only: verts_rphiz, verts_sthetaphi, ntetr
    use find_tetra_mod, only: find_tetra
    use tetra_grid_settings_mod, only: grid_kind
    use tetra_physics_mod, only: coord_system
    use utils_data_pre_and_post_processing_mod, only: set_seed_for_random_numbers

    double precision, intent(in)                                   :: v0
    double precision, dimension(:,:), allocatable, intent(out)     :: start_pos_pitch_mat
    double precision                                               :: rand_scalar, vpar, vperp
    double precision                                               :: amin, cmin, cmax, amax !is set globally
    double precision, dimension(:), allocatable                    :: rand_vector
    double precision, dimension(:), allocatable                    :: rand_matrix2
    integer                                                        :: i
    double precision, dimension(3)                                 :: x
    integer                                                        :: ind_tetr_out,iface

!!!!comment out the following call to make starting conditions really random!!!
    call set_seed_for_random_numbers

    allocate(start_pos_pitch_mat(5,num_particles))
    allocate(rand_vector(num_particles))
    allocate(rand_matrix2(num_particles))

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

    constant_part_of_weights = (amax-amin)*(cmax-cmin)*2*pi

    !compute starting conditions
    if (boole_start_points_file) then
        !odd n: start point (n+1)/2 along B (pitch +1), even n: same start point against B (pitch -1)
        do i = 1,num_particles
            start_pos_pitch_mat(1:3,i) = start_points_rphiz(:,(i+1)/2)
        enddo
        if (coord_system.eq.2) print*, 'error: start-points file is only implemented for cylindrical coordinate system'
    elseif (boole_point_source) then
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

    start_pos_pitch_mat(4,:) = 1 !delete this once i have a proper subroutine for field line tracing
    if (boole_start_points_file) start_pos_pitch_mat(4,2:num_particles:2) = -1

    call RANDOM_NUMBER(rand_matrix2)
    if (.not. boole_monoenergetic) then !compare with equation 133 of master thesis of Jonatan Schatzlmayr (remaining parts will be added later)
        start_pos_pitch_mat(5,:) = 5*energy_eV*rand_matrix2(:) !boltzmann energy distribution
        constant_part_of_weights = constant_part_of_weights*10/sqrt(pi)*energy_eV*ev2erg
    endif

    weights(:,1) = constant_part_of_weights
    if (boole_refined_sqrt_g.eqv..false.) weights(:,1) = constant_part_of_weights*start_pos_pitch_mat(ind_a,:)

end subroutine calc_starting_conditions

end module field_line_tracing_mod