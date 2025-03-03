module sub_alpha_lifetime_can_mod

      implicit none

      contains 
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      subroutine elefie_can(x,derphi)
!
! Computes the derivatives of the electrostatic potential over
! coordinates (covariant electric field). Potential is normalized
! to reference temperature : phinor = e phi / T.
!
!   Input parameters:
!             formal:    x       -   array of coordinates
!
!   Output parameters:
!             formal:    derphi  -   array of derivatives
!
    use gorilla_settings_mod, only: eps_Phi
    use tetra_physics_mod, only: particle_charge,cm_over_e,particle_mass
    use parmot_mod, only: ro0
    use splint_vmec_data_mod, only: vmec_field
!
    implicit none
!
    integer :: ierr
    double precision, dimension(3) :: x,derphi
    double precision :: r,phi,z,psi,phi_el,phi_el_pr,phi_el_prpr,temp
    double precision :: s,theta,varphi,A_theta,A_phi,dA_theta_ds,dA_phi_ds,aiota,     &
                      sqg,alam,dl_ds,dl_dt,dl_dp,Bctrvr_vartheta,Bctrvr_varphi,     &
                      Bcovar_r,Bcovar_vartheta,Bcovar_varphi
!
   s = x(1)
   theta = x(2)
   varphi = x(3)
!
   call vmec_field(s,theta,varphi,A_theta,A_phi,dA_theta_ds,dA_phi_ds,aiota,     &
               sqg,alam,dl_ds,dl_dt,dl_dp,Bctrvr_vartheta,Bctrvr_varphi,     &
               Bcovar_r,Bcovar_vartheta,Bcovar_varphi)
!
   !Derivatives of Phi with respect to theta and varphi are set to zero. 
   derphi=0.d0
!
   !Temperature
   temp = (ro0/cm_over_e)**2*particle_mass/2.d0
!
   !Derivative of Phi with respect to s
   derphi(1) = dA_theta_ds * eps_Phi * & ! dPhi_ds
               & particle_charge/(temp)          ! e/T
!
      end subroutine elefie_can
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      subroutine velo_can(tau,z,vz)
!
!
!  Computes the components of the 5-D velocity          -  vz
!  for given set of phase space coordinates             -  z.
!
!  Warning !!!  The dimensionless time is used (see below)
!
!  Order of the coordinates is the following:
!  z(i) = x(i)  for i=1,2,3     - spatial coordinates with real
!                                 dimension of the general covariant
!                                 space coordinate system
!  z(4) = p                     - momentum  module normalized to
!                                 thermal momentum and sqrt(2);
!  z(5) = alambd                - cosine of the pitch-angle
!
!  Input parameters:
!            formal:  tau    -   dimensionless time: tau=sqrt(2*T/m)*t
!                     z      -   see above
!            common:  rmu    -   inverse relativistic temperature
!                     ro0    -   Larmor radius for the reference
!                                magnetic field and temperature:
!                                ro0=sqrt(2*T/m)/(e*B_ref/m*c)
!                     p0     -   dimensionless momentum module in the
!                                initial point
!                     alamb0 -   cos of pitch-angle in the initial point
!  Output parameters:
!            formal:  vz     -   see above
!
!  Called routines: magfie_can, magfie_vmec, elefie_can
!
      use parmot_mod, only : rmu,ro0
      use velo_mod,   only : isw_field_type
      use magfie_mod, only: magfie_vmec, magfie_can
!
      implicit none
!
      integer :: i
!
      double precision tau,z,vz
      double precision x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl
      double precision derphi
      double precision p,alambd,p2,ovmu,gamma2,gamma,ppar,vpa,coala
      double precision rmumag,rovsqg,rosqgb,rovbm
      double precision a_phi,a_b,a_c,hstar
      double precision s_hc,hpstar,phidot,blodot,bra
      double precision pardeb
!
      dimension z(5),vz(5)
      dimension x(3),bder(3),hcovar(3),hctrvr(3),hcurl(3)
      dimension derphi(3)
      dimension a_phi(3),a_b(3),a_c(3),hstar(3)
!
      do 1 i=1,3
        x(i)=z(i)
 1    continue
!
! in magfie: x(i)   - set of 3 curvilinear space coordinates (input)
!            bmod   - dimensionless magnetic field module: bmod=B/B_ref
!            sqrtg  - Jacobian of space coordinates (square root of
!                     metric tensor
!            bder   - derivatives of logarithm of bmod over space coords 
!                     (covariant vector)
!            hcovar - covariant components of the unit vector along
!                     the magnetic field
!            hctrvr - contravariant components of the unit vector along
!                     the magnetic field
!            hcurl  - contravariant components of the curl of this vector
!
      if(isw_field_type.eq.0) then
!
        call magfie_can(x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
!
      elseif(isw_field_type.eq.1) then
!
        call magfie_vmec(x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
!
      else
        print *,'velo: unknown field type'
        return
      endif
!
! in elefie: x(i)   - space coords (input, see above)
!            derphi - derivatives of the dimensionless electric potential
!                     phihat=e*phi/T over space coords (covar. vector)
!
      call elefie_can(x,derphi)
!
      p=z(4)
      alambd=z(5)
!
      p2=p*p
      ovmu=2.d0/rmu
      gamma2=p2*ovmu+1.d0
      gamma=dsqrt(gamma2)
      ppar=p*alambd
! vpa - dimensionless parallel velocity: vpa=v_parallel/sqrt(2*T/m)
      vpa=ppar/gamma
      coala=(1.d0-alambd**2)
! rmumag - magnetic moment
      rmumag=.5d0*p2*coala/bmod
!
      rovsqg=ro0/sqrtg
      rosqgb=.5d0*rovsqg/bmod
      rovbm=ro0/bmod
!
      a_phi(1)=(hcovar(2)*derphi(3)-hcovar(3)*derphi(2))*rosqgb
      a_b(1)=(hcovar(2)*bder(3)-hcovar(3)*bder(2))*rovsqg
      a_phi(2)=(hcovar(3)*derphi(1)-hcovar(1)*derphi(3))*rosqgb
      a_b(2)=(hcovar(3)*bder(1)-hcovar(1)*bder(3))*rovsqg
      a_phi(3)=(hcovar(1)*derphi(2)-hcovar(2)*derphi(1))*rosqgb
      a_b(3)=(hcovar(1)*bder(2)-hcovar(2)*bder(1))*rovsqg
!
      s_hc=0.d0
      do i=1,3
        a_c(i)=hcurl(i)*rovbm
        s_hc=s_hc+a_c(i)*hcovar(i)
        hstar(i)=hctrvr(i)+ppar*a_c(i)
      enddo
      hpstar=1.d0+ppar*s_hc
!
! velocities in the coordinate space
!
! phidot - derivative of the dmls el. potential over dmls time
! blodot - derivative of the logarith of the mag. field module over dmls time
      phidot=0.d0
      blodot=0.d0
      do i=1,3
        bra=vpa*hstar(i)+a_phi(i)+a_b(i)*rmumag/gamma
        vz(i)=bra/hpstar
        phidot=phidot+vz(i)*derphi(i)
        blodot=blodot+vz(i)*bder(i)
      enddo
!
! velocities in the phase space
!
      vz(4)=-0.5d0*gamma*phidot/p
!      vz(5)=coala/(alambd+dsign(1.d-32,alambd))*(vz(4)/p-0.5d0*blodot)
      vz(5)=-(0.5d0*coala/hpstar)*(sum(hstar*derphi)/p                 &
            + p*sum(hstar*bder)/gamma+alambd*sum(a_phi*bder))
!
      end subroutine velo_can
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      subroutine orbit_timestep_can(z,dtau,dtaumin,ierr,tau_out)
!
      use odeint_mod, only: odeint_allroutines
      use runge_kutta_mod, only: runge_kutta_allroutines
      use gorilla_settings_mod, only: rel_err_ode45
      use chamb_can_mod, only: chamb_can
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
      end subroutine orbit_timestep_can
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      subroutine rhs_mflint_can(phi,y,dery)
!
! Computes the right hand side of the magnetic field line equation for
! the integration over the toroidal angle, subintegrand for the flux tube 
! volume $1/B^\varphi$ and subintegrants for Boozer $B_{00}$ computation
! $B^2/B^\varphi$ and $B^3/B^\varphi$                       -   dery
!
! Oder of coordinates in magfie is the following: x(1)=R (big radius),
! x(2)=phi (toroidal angle), x(3)=z (altitude).
!
!  Input parameters:
!            formal:
!                    y(1:2) - coordinates in the poloidal plane (phi=const)
!                    y(3:5) - integrals
!  Output parameters:
!            formal:
!                 dery(1:5) - vector of the right-hand side
!  Called routines:  magfie_vmec, magfie_can
!
      use velo_mod,   only : isw_field_type
      use magfie_mod, only: magfie_vmec, magfie_can
!
      double precision :: phi
      double precision, dimension(5) :: y,dery
      double precision x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl
      dimension x(3),bder(3),hcovar(3),hctrvr(3),hcurl(3)
!
      x(1)=y(1)
      x(2)=y(2)
      x(3)=phi
!
      if(isw_field_type.eq.0) then
!
        call magfie_can(x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
!
      elseif(isw_field_type.eq.1) then
!
        call magfie_vmec(x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
!
      else
        print *,'rhs_mflint_can: unknown field type'
        return
      endif
!
      dery(1)=hctrvr(1)/hctrvr(3)
      dery(2)=hctrvr(2)/hctrvr(3)
      dery(3)=1.d0/(bmod*hctrvr(3))
      dery(4)=bmod/hctrvr(3)
      dery(5)=bmod**2/hctrvr(3)
!
      end subroutine rhs_mflint_can
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine integrate_mfl_can(npoi,dphi,rbeg,phibeg,zbeg,         &
                               xstart,bstart,volstart,bmod00,ierr)
!
      use velo_mod,   only : isw_field_type
      use odeint_mod, only: odeint_allroutines
      use chamb_can_mod, only: chamb_can
      use magfie_mod, only: magfie_vmec, magfie_can
!
      implicit none
!
      integer, parameter          :: ndim=5
      double precision, parameter :: relerr=1d-8
      integer :: npoi,i,ierr
      double precision :: dphi,rbeg,phibeg,zbeg,bmod00,phi,phiold
      double precision, dimension(3,npoi) :: xstart
      double precision, dimension(npoi)   :: bstart,volstart
      double precision, dimension(ndim)   :: y
!
      double precision x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl
      dimension x(3),bder(3),hcovar(3),hctrvr(3),hcurl(3)
!
      ierr=0
!
      phi=phibeg
      y(1)=rbeg
      y(2)=zbeg
      y(3:5)=0.d0
!
      call chamb_can(y(1:2),phi,ierr)
!
      if(ierr.ne.0) return
      x(1)=y(1)
      x(2)=y(2)
      x(3)=phi
!
      if(isw_field_type.eq.0) then
!
        call magfie_can(x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
!
      elseif(isw_field_type.eq.1) then
!
        call magfie_vmec(x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
!
      else
        print *,'rhs_mflint_can: unknown field type'
        return
      endif
!
      xstart(:,1)=x
      bstart(1)=bmod
      volstart(1)=y(3)
!
      do i=2,npoi
        phiold=phi
        phi=phiold+dphi
!
        call odeint_allroutines(y,ndim,phiold,phi,relerr,rhs_mflint_can)
!
        call chamb_can(y(1:2),phi,ierr)
!
        if(ierr.ne.0) return
        x(1)=y(1)
        x(2)=y(2)
        x(3)=phi
!
        if(isw_field_type.eq.0) then
!
          call magfie_can(x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
!
        elseif(isw_field_type.eq.1) then
!
          call magfie_vmec(x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
!
        else
          print *,'rhs_mflint_can: unknown field type'
          return
        endif
!
        xstart(:,i)=x
        bstart(i)=bmod
        volstart(i)=y(3)
      enddo
!
      volstart=volstart/volstart(npoi)
      bmod00=y(5)/y(4)
!
      end subroutine integrate_mfl_can

end module sub_alpha_lifetime_can_mod
