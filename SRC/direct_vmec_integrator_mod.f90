!
module direct_vmec_integrator_mod
!
    contains
!
    subroutine direct_vmec_integrator()
!
        use orbit_timestep_gorilla_mod, only: initialize_gorilla
        use constants, only: ev2erg, pi
        use tetra_physics_mod, only: particle_mass,cm_over_e,mag_axis_R0
        use parmot_mod, only: rmu,ro0
        use velo_mod, only: isw_field_type
        use supporting_functions_mod, only: theta_sym_flux2theta_vmec,theta_vmec2theta_sym_flux
        use tetra_grid_settings_mod, only: n_field_periods
!
        implicit none
!
        integer :: i,n_timesteps,counter_mappings_phi_0,counter_mappings_vpar_0,ierr,n_mappings
        double precision :: vmod, energy_eV,time_step,tau_out_can,trace_time
        double precision :: dtau, dphi,dtaumin
        double precision, dimension(5) :: z
        double precision :: s,theta,phi,alam0
        double precision :: varphi_save1,varphi_save2,delta_phi,theta_symflux,phi_save,zerophi
        logical :: boole_poincare_phi0,boole_poincare_vpar0
!
!-----------------------------------------------------------------------------------------------------------!
! Prepare calculation of orbit tip by interpolation
!
      integer                                       :: nplagr,nder,itip,npl_half
      integer                                       :: ifp,npassing,ntr_regular,ntr_chaotic
      double precision                              :: alam_prev,zerolam,twopi,fraction
      double precision, dimension(6)                :: z_tip
      integer,          dimension(:),   allocatable :: ipoi
      double precision, dimension(:),   allocatable :: xp
      double precision, dimension(:,:), allocatable :: coef,orb_sten
!
      zerolam=0.d0
      zerophi=0.d0
      twopi=2.d0*pi
      nplagr=4
      nder=0
      npl_half=nplagr/2
      allocate(ipoi(nplagr),coef(0:nder,nplagr),orb_sten(5,nplagr),xp(nplagr))
      do i=1,nplagr
        ipoi(i)=i
      enddo
!
!------------------------------------------------------------------------------------------------------------!
! SETTINGS
!
        !Plotting options
        boole_poincare_phi0 = .true.
        boole_poincare_vpar0 = .false.
!
        if(boole_poincare_phi0) open(4002,file='poincare_plot_phi_0_sthetaphi.dat')
        if(boole_poincare_vpar0) open(4001,file='poincare_plot_vpar_0_sthetaphi.dat')
!
        !time for tracing alpha particle
!
        n_mappings = 10 
!
        trace_time = 1.0d-2 / 35.d0 *dble(n_mappings)
        n_timesteps = 30000 * n_mappings 
        time_step = trace_time/dble(n_timesteps)
!
print *, 'trace_time', trace_time
print *, 'n_timesteps',n_timesteps

        !Initialize GORILLA
        call initialize_gorilla()
!
        !Compute velocity module from kinetic energy dependent on particle species
        energy_eV = 3.5d6
        vmod=sqrt(2.d0*energy_eV*ev2erg/particle_mass)
!
print *, 'cm_over_e',cm_over_e
!
        !Guiding-center position and pitch parameter
        s = 0.62d0
        theta = 0.1d0
        phi = 0.1d0
        alam0 = 0.95d0
!
!------------------------------------------------------------------------------------------------------------!
! Initialization of direct integrator
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
        dphi=2.d0*pi/dble(30)/n_field_periods
!
        !orbit integration time step (to check chamber wall crossing)
        dtaumin=dphi*mag_axis_R0
!
!------------------------------------------------------------------------------------------------------------!
! Particle coordinates in VMEC
!
print *, 'theta_symflux_1',theta
        z(1) = s
        z(2) = theta_sym_flux2theta_vmec(s,theta,phi)
        z(3) = phi
        z(4) = 1.d0
        z(5) = alam0
!
print * ,'theta_symflux_2', theta_vmec2theta_sym_flux(z(1),mod(z(2),2.d0*pi),mod(z(3),2.d0*pi))

print *, 'z',z
print *, 'ro0',ro0


        counter_mappings_phi_0 = 0
        counter_mappings_vpar_0 = 0
!
!------------------------------------------------------------------------------------------------------------!
! Initialize tip detector

    ifp = 0
    itip=3
    alam_prev=z(5)
!
!------------------------------------------------------------------------------------------------------------!
!

!isw_field_type=0
!call spline_vmec_data
!isw_field_type=1



    do i=1,n_timesteps
        varphi_save1 = mod(z(3),2.d0*pi/n_field_periods)
!
        call orbit_timestep_can(z,dtau,dtaumin,ierr,tau_out_can)
!
        !Poincaré sections at v_par = 0
        if(boole_poincare_vpar0) then
!
            if(alam_prev.lt.0.d0.and.z(5).gt.0.d0) then
                itip=0   !<=tip has been passed
                counter_mappings_vpar_0 = counter_mappings_vpar_0 + 1
            endif
!
            itip=itip+1
            alam_prev=z(5)
    
            if(i.le.nplagr) then          !<=first nplagr points to initialize stencil
                orb_sten(1:5,i)=z
            else                          !<=normal case, shift stencil
                orb_sten(1:5,ipoi(1))=z
                ipoi=cshift(ipoi,1)
                if(itip.eq.npl_half) then   !<=stencil around tip is complete, interpolate
                    xp=orb_sten(5,ipoi)
!
                    call plag_coeff(nplagr,nder,zerolam,xp,coef)
                    z_tip=matmul(orb_sten(:,ipoi),coef(0,:))
!
                    theta_symflux = theta_vmec2theta_sym_flux(z_tip(1),mod(z_tip(2),2.d0*pi),mod(z_tip(3),2.d0*pi/n_field_periods))
                    write (4001,*) z_tip(1),theta_symflux,mod(z_tip(3),2.d0*pi/n_field_periods)
!
                endif
            endif
!
        endif !boole_poincare_vpar0
!
        !Poincaré plots for varphi = 0
        if(boole_poincare_phi0) then
            varphi_save2 = mod(z(3),2.d0*pi/n_field_periods)
            delta_phi = varphi_save2-varphi_save1
!
            if(abs(delta_phi).gt.(pi/n_field_periods)) then
!
                itip=0          !Poincaré plane has been passed
                phi_save = z(3) !Save toroidal phi angle after the Poincaré plane has been passed
!
                !Plotting of Poincaré cuts
                if(delta_phi.gt.0.d0) then
                    counter_mappings_phi_0 = counter_mappings_phi_0 - 1
                elseif(delta_phi.lt.0.d0) then
                    counter_mappings_phi_0 = counter_mappings_phi_0 + 1
                endif
            endif
!
            itip=itip+1
!
            if(i.le.nplagr) then          !<=first nplagr points to initialize stencil
                orb_sten(1:5,i)=z
            else                          !<=normal case, shift stencil
                orb_sten(1:5,ipoi(1))=z
                ipoi=cshift(ipoi,1)
                if(itip.eq.npl_half) then   !<=stencil around tip is complete, interpolate
                xp=orb_sten(3,ipoi)
!
                zerophi = floor(phi_save / (2.d0*pi/n_field_periods)) * 2.d0*pi/n_field_periods
!
!print *, 'zerophi',zerophi, 'phi_save', phi_save, '2.d0*pi/n_field_periods', zerophi * 2.d0*pi/n_field_periods
!
                call plag_coeff(nplagr,nder,zerophi,xp,coef)
                z_tip=matmul(orb_sten(:,ipoi),coef(0,:))
!
                theta_symflux = theta_vmec2theta_sym_flux(z_tip(1),mod(z_tip(2),2.d0*pi),mod(z_tip(3),2.d0*pi/n_field_periods))
                write (4002,*) z_tip(1),theta_symflux,mod(z_tip(3),2.d0*pi/n_field_periods)
!
                endif
            endif
        endif !boole_poincare_phi0
!
    enddo !
!
if(boole_poincare_phi0) then
    print *, 'Number of toroidal mappings:',counter_mappings_phi_0
    print *, 'Average number of DIRECT time steps per toroidal mapping:', n_timesteps/counter_mappings_phi_0
endif
!
if(boole_poincare_vpar0) then
    print *, 'Number of banana bounces:',counter_mappings_vpar_0
    print *, 'Average number of DIRECT time steps per banana bounce:', n_timesteps/counter_mappings_vpar_0
endif
!
if(boole_poincare_phi0) close(4002)
if(boole_poincare_vpar0) close(4001)
!
    end subroutine direct_vmec_integrator
!
end module direct_vmec_integrator_mod
