!
  module fluxtv_pusher_mod
!
    implicit none
!
    private
!
    integer                              :: ind_tetr
    double precision                     :: spamat,B0,dt_dtau_const
    double precision, dimension(3)       :: Bvec,gradB,x_init,x_ver1
    double precision, dimension(4)       :: b
    double precision, dimension(3,3)     :: amat
    double precision                     :: perpinv,perpinv2,tau
    double precision, dimension(5)       :: z_init
!
    public :: pusher_tetr_flux_tube_volume,initialize_const_motion_fluxtv
!
  contains
!
    subroutine initialize_const_motion_fluxtv(perpinv_in,perpinv2_in)
!
        implicit none
!            
        double precision, intent(in)    :: perpinv_in,perpinv2_in
!
        perpinv = perpinv_in
        perpinv2 = perpinv2_in
!            
    end subroutine initialize_const_motion_fluxtv
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine initialize_pusher_fluxtv(ind_tetr_in,x,vpar,z_final,t_pass,fluxtv_tetra)
!    
        use gorilla_settings_mod, only: boole_dt_dtau
        use constants, only: clight
        use tetra_physics_mod, only: tetra_physics, cm_over_e,dt_dtau
!
        implicit none
!
        integer                          :: ind_tetr_in
        double precision,dimension(3)    :: x,x_exit,z_final
        double precision                 :: vpar,t_pass,fluxtv_tetra,bmod,phi_elec,vpar2,vperp2
!
        ind_tetr = ind_tetr_in
!
        x_init = x
!        
        x_ver1 = tetra_physics(ind_tetr)%x1
!
        z_init(1:3) = x_init-x_ver1
        z_init(4) = vpar
        z_init(5) = fluxtv_tetra
!
        !Tetrahedron constants
        B0 = tetra_physics(ind_tetr)%bmod1
        gradB = tetra_physics(ind_tetr)%gb
        
        !Calculation of orbit parameter dependent on boole_dt_dtau
        if(boole_dt_dtau) then
            dt_dtau_const = tetra_physics(ind_tetr)%dt_dtau_const
        else
            x_exit = x_ver1 + z_final
            dt_dtau_const = dt_dtau(ind_tetr,x_init,x_exit)
        endif 
        tau = t_pass/dt_dtau_const
!    
        !Module of B at the entry point of the particle
        bmod=tetra_physics(ind_tetr)%bmod1+sum(tetra_physics(ind_tetr)%gb*z_init(1:3))
!
        !Phi at the entry point of the particle
        phi_elec=tetra_physics(ind_tetr)%Phi1+sum(tetra_physics(ind_tetr)%gPhi*z_init(1:3))
!
!       !Auxiliary quantities
        vperp2 = -2.d0*perpinv*bmod
        vpar2 = vpar**2
!
        !ODE coefficients
        b(1:3)=( tetra_physics(ind_tetr)%curlh*(vperp2+vpar2+2.d0*perpinv*B0) &
            + perpinv*tetra_physics(ind_tetr)%gBxh1 )*cm_over_e &
            - clight*(2.d0*(tetra_physics(ind_tetr)%Phi1-phi_elec)*tetra_physics(ind_tetr)%curlh &
            + tetra_physics(ind_tetr)%gPhixh1)
!
        b(4)=perpinv*tetra_physics(ind_tetr)%gBxcurlA-clight/cm_over_e*tetra_physics(ind_tetr)%gPhixcurlA
!
        Bvec=tetra_physics(ind_tetr)%curlA    !B-Vector
        amat=perpinv*cm_over_e*tetra_physics(ind_tetr)%alpmat &
            - clight* tetra_physics(ind_tetr)%betmat    ! a-Matrix (elements 1:3)
!                
        spamat=perpinv*cm_over_e*tetra_physics(ind_tetr)%spalpmat &
            - clight* tetra_physics(ind_tetr)%spbetmat    !tr(a-Matrix)
!    
    end subroutine initialize_pusher_fluxtv
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine pusher_tetr_flux_tube_volume(ind_tetr_in,x,z_final,vpar,t_pass,fluxtv_tetra)
!
      use odeint_mod, only: odeint_allroutines
!
      implicit none
!
      double precision, dimension(3)   :: x,z_final
      double precision, dimension(5) :: z
      double precision :: fluxtv_tetra,tau_old,relerr,vpar,t_pass
      integer :: ndim,ind_tetr_in
!
      call initialize_pusher_fluxtv(ind_tetr_in,x,vpar,z_final,t_pass,fluxtv_tetra)
!
      z = z_init
!
      ndim = 5             !Number of variables for ODE45-Integrator
      relerr = 1.d-8       !Relative accuracy of RK45 integration
!
      tau_old = 0.d0
!
      call odeint_allroutines(z,ndim,tau_old,tau,relerr,rhs_pusher_fluxtv)
!
      x = z(1:3)+x_ver1
      if(abs(z(4)-vpar).gt.1.d-12) print *, 'Error in pusher_tetr_flux_tube_volume: vpar changes during integration of orbit.'
      vpar = z(4)
      fluxtv_tetra = z(5)
!
    end subroutine pusher_tetr_flux_tube_volume
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine rhs_pusher_fluxtv(dummy,z,dzdtau)
!
      implicit none
!
      double precision, dimension(5) :: z,dzdtau
      double precision :: dummy
!
      dzdtau(1:3) = b(1:3)+matmul(amat,z(1:3))+Bvec*z(4)
      dzdtau(4) = b(4)+spamat*z(4)
      dzdtau(5) = z(4)*dt_dtau_const/( B0 +sum(gradB*z(1:3)) ) !Flux tube volume
!
    end subroutine rhs_pusher_fluxtv
!
  end module fluxtv_pusher_mod
!
