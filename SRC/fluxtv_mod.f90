!
  module fluxtv_mod
!
    implicit none
!
    private
!
    double precision, dimension(:,:), allocatable, public,protected :: pos_fluxtv_mat
    integer, dimension(2), public,protected :: pos_fluxtv_mat_shape
    double precision, dimension(5), public,protected :: grid_configuration
!
    public :: calc_flux_tube_volume,load_flux_tube_volume
!
  contains
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine calc_flux_tube_volume()
!
      use tetra_grid_settings_mod, only: grid_kind,grid_size
      use tetra_physics_mod, only: tetra_physics,coord_system,particle_mass
      use fluxtv_pusher_mod, only: pusher_tetr_flux_tube_volume,initialize_const_motion_fluxtv
      use pusher_tetra_rk_mod, only: pusher_tetra_rk,initialize_const_motion_rk
      use orbit_timestep_gorilla_mod, only: find_tetra, initialize_gorilla
      use supporting_functions_mod, only : sym_flux_in_cyl
      use constants, only: ev2erg
      use gorilla_applets_settings_mod, only: filename_fluxtv_precomp, start_pos_x1, start_pos_x2, start_pos_x3, &
        & t_step => t_step_fluxtv, energy_eV => energy_eV_fluxtv, nt_steps => nt_steps_fluxtv
!
      implicit none
!
      logical :: file_exist,boole_t_finished
      integer :: iface,i,j,k,l,m,n,ind_tetr,ind_tetr_save,iper,io_error
      integer :: ind_tetr_temp,ipert
      integer(kind=8) :: counter_mappings,counter_t_steps
      double precision :: phi_start, R_torus,R0,Z0,t_pass,t_remain, &
                          vmod,vpar,vperp,fluxtv_tetra,t_sum, &
                          vpar_temp,vperp_temp
      double precision :: bmod_multiplier
!      double precision :: rrr,ppp,zzz,B_r,B_p,B_z,dBrdR,dBrdp,dBrdZ    &
!                          ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ
      double precision, dimension(3) :: x,x_temp,z_final
!
      !Starting position for fluxtube volume computation
      x = [start_pos_x1, start_pos_x2, start_pos_x3]
!
      allocate(pos_fluxtv_mat(nt_steps,5)) !IMPORTANT nt_steps,4
!
      !Initialize grid
      call initialize_gorilla(1)
!
      !Option to calculate magnetic field line (0th order drift orbits)
      bmod_multiplier = 1.d5 !For normal conditions (drift orbits) 1.d0
      ipert = 0
!
      call initialize_gorilla(2,ipert,bmod_multiplier)
!
      !Save grid configuration
      grid_configuration = [grid_kind,coord_system,grid_size(1),grid_size(2),grid_size(3)]
!
      !Compute vmod from energy_eV
      vmod =   sqrt(2.d0*energy_eV*ev2erg/particle_mass)
      vpar = 1.d0*vmod
      vperp = 0.d0
      fluxtv_tetra = 0.d0 !Initialize effective volume (tube flux volume)
!
      counter_mappings = 0       ! Counter for mappings
      counter_t_steps = 0
!
      !Cumulated orbit flight time
      t_sum = 0.d0
!
      !Initialize total time of orbit flight
      t_remain = t_step
!
      call find_tetra(x,vpar,vperp,ind_tetr,iface)
!   
      !Set constants of motion in module (vperp = 0)
      call initialize_const_motion_rk(0.d0,0.d0)
      call initialize_const_motion_fluxtv(0.d0,0.d0)
!
      print *, 'Start flux tube volume calculation with configuration', grid_configuration
!
      !Loop over tetrahedron passes until t_step is reached
      do
!
        ind_tetr_save=ind_tetr
!
        !Domain Boundary
        if(ind_tetr.eq.-1) then
          print *,'Error: Particle has left the domain.'
          exit
        endif
!
        !t_remain (in) ... remaining time until t_step is finished
        !t_pass (out) ... time to pass the tetrahdron
!
        !Save certain variables
        x_temp = x
        vpar_temp = vpar
        ind_tetr_temp = ind_tetr
!
        !Calculate trajectory
        call pusher_tetra_rk(ind_tetr,iface,x,vpar,z_final,t_remain,t_pass,boole_t_finished,iper)
!        
        !Calculate flux tube volume with adaptive ODE integrator
        call pusher_tetr_flux_tube_volume(ind_tetr_temp,x_temp,z_final,vpar_temp,t_pass,fluxtv_tetra)
!
        t_remain = t_remain - t_pass
        t_sum = t_sum + t_pass
!
        !Orbit stops within cell, because "flight"-time t_step has finished
        if(boole_t_finished) then
          !print *, 'Stop within cell detected in main program'
          counter_t_steps = counter_t_steps +1
          t_remain = t_step
          pos_fluxtv_mat(counter_t_steps,4) = fluxtv_tetra
          pos_fluxtv_mat(counter_t_steps,1:3) = x
          pos_fluxtv_mat(counter_t_steps,5) = t_sum
        endif
!
        if(iper.ne.0) then
            counter_mappings = counter_mappings + iper
        endif
!
        if(counter_t_steps.gt.nt_steps-1) then
!           print *,"Exit (timesteps): Maximum number of time-steps reached."
          exit
        endif
!
      enddo !pusher_tetr
!
      print *,'number of mappings', counter_mappings
      print *,'number of time steps', counter_t_steps
!
      !Normalize effective volume
      pos_fluxtv_mat(:,4) = pos_fluxtv_mat(:,4)/pos_fluxtv_mat(nt_steps,4)
!
! ---- Write Orbit Position and Flux Tube Volume in file ----!
!
      inquire(file=filename_fluxtv_precomp ,exist = file_exist)
      if(file_exist) then
        open(unit=20,file=filename_fluxtv_precomp,status='old',iostat=io_error)
        if(io_error.eq.0) close(unit=20,status='delete')
      else
        !print *,'Error: Could not open and delete existing file'
        print *,'No existing file. Created new one to store data.'
      endif
      open(unit=20,file=filename_fluxtv_precomp,status='new',action='write', &
      & iostat=io_error)
!
      if ( io_error == 0) then
        write(20,*) shape(pos_fluxtv_mat(:,1:4))
        write(20,*) grid_configuration
        do i = 1,nt_steps
          write(20,*) pos_fluxtv_mat(i,1:4)
        enddo
      else
        write(*,*) 'Beim OEffenen der Datei ist ein Fehler Nr.', &
        io_error,' aufgetreten'
      end if
      close(unit=20)
!
!       !If symmetry flux coordinates are used, transform to cylindrical coordinates
!       if(coord_system.eq.2) then
!         call sym_flux_in_cyl('fluxtubevolume.dat','fluxtubevolume_rphiz.dat',1)
!       endif
! ----------------------------------------------------------------!
!
!      !Option to calculate magnetic field line (0th order drift orbits)
!      !Set back to standard bmod-value
!      bmod_multiplier = 1.d0
!
    end subroutine calc_flux_tube_volume
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine load_flux_tube_volume()
!
        use gorilla_applets_settings_mod, only: filename_fluxtv_load
        use gorilla_settings_mod, only: coord_system
        use tetra_grid_settings_mod, only: grid_kind,grid_size
!
        implicit none
!
        integer :: i
!
        if(.not.(allocated(pos_fluxtv_mat))) then
            !Read the coordinates data
            open(unit = 22, file=filename_fluxtv_load, status='old',action = 'read')
            read(22,*) pos_fluxtv_mat_shape(1),pos_fluxtv_mat_shape(2)
            read(22,*) grid_configuration
!
            !Check, if Fluxtube volume was correctly computed before diffusion coefficent computation
            if(any(grid_configuration.ne.[grid_kind,coord_system,grid_size(1),grid_size(2),grid_size(3)])) then
                print *, 'Properties of fluxtube data:', grid_configuration
                print *, 'Properites of requested grid:', [grid_kind,coord_system,grid_size(1),grid_size(2),grid_size(3)]
                print *, 'Error: Fluxtube volume must be computed before diffusion measurement'
                stop
            endif
!
            allocate(pos_fluxtv_mat(1:pos_fluxtv_mat_shape(1),1:pos_fluxtv_mat_shape(2)))
            do i=1,pos_fluxtv_mat_shape(1)
                read(22,*) pos_fluxtv_mat(i,1:pos_fluxtv_mat_shape(2))
            enddo
            close(unit = 22)
        endif
!
    end subroutine load_flux_tube_volume
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  end module fluxtv_mod
