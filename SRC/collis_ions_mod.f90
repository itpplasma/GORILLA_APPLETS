module collis_ions

      use, intrinsic :: iso_fortran_env, only: dp => real64

implicit none

contains
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine coleff(efcolf,velrat,enrat,p,dpp_vec,dhh_vec,fpeff_vec)

!  Computes local values of dimensionless contravariant components
!  of collisional diffusion tensor and friction force for nonrelativistic
!  plasma. Backgound temperature is the same for all sorts.
!
!     Input variables:
!        formal: p      - dimensionless momentum module (p/(sqrt(2)*p_T) (=v0)
!        common: efcolf - dmls collision frequencies
!                velrat - ratio of test species thermal velocity to
!                         background species thermal velocity
!     Output variables:
!        formal: dpp    - dimensionless momentum module diffusion
!                         coefficient
!                dhh    - dimensionless pitch angle diffusion coeff.
!                fpeff  - effective dimensionless drag force (prop. to linear
!                         deviation in Fokker-Planck eq.)

  integer :: i, n
  real(dp), dimension(:), intent(in) :: efcolf,velrat,enrat
  real(dp), dimension(:), intent(out) :: dpp_vec,dhh_vec,fpeff_vec
  real(dp) :: p,plim,xbeta,dpd

  plim=max(p,1.d-8)
  n = size(efcolf)

  do i=1,n
    xbeta=p*velrat(i)

    call onseff(xbeta,dpp_vec(i),dhh_vec(i),dpd)

    fpeff_vec(i) = (dpd/plim-2.0*dpp_vec(i)*p*enrat(i))*efcolf(i)
    dpp_vec(i) = dpp_vec(i)*efcolf(i)
    dhh_vec(i) = dhh_vec(i)*efcolf(i)
  enddo

  dhh_vec = dhh_vec/plim**2

end subroutine coleff
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine onseff(v,dp,dh,dpd)

!  dp - dimensionless dpp
!  dh - dhh*p^2     (p - dmls)
!  dpd - (1/p)(d/dp)p^2*dp   (p - dmls)

! square root of pi
  double precision, parameter :: sqp=1.7724538d0
! cons=4./(3.*sqrt(pi))
  double precision, parameter :: cons=.75225278d0
  double precision :: v,dp,dh,dpd,v2,v3,ex,er

  v2=v**2
  v3=v2*v
  if(v.lt.0.01d0) then
    dp=cons*(1.d0-0.6d0*v2)
    dh=cons*(1.d0-0.2d0*v2)
    dpd=2.d0*cons*(1.d0-1.2d0*v2)
  elseif(v.gt.6.d0) then
    dp=1.d0/v3
    dh=(1.d0-0.5d0/v2)/v
    dpd=-1.d0/v3
  else
    ex=exp(-v2)/sqp
    er=erf(v)
    dp=er/v3-2.d0*ex/v2
    dh=er*(1.d0-0.5d0/v2)/v+ex/v2
    dpd=4.d0*ex-dp
  endif

end subroutine onseff

subroutine collis_init(m1,Z1,m,Z,dens,temp,e0,v0,efcolf,velrat,enrat)

!   Performs precomputation of the constants for Coulomb collision
!   operator for test particles colliding with n-1 sorts of ions and with electrons
!
!   Normalisation: test particle velocity is normalized by v0,
!   time is multiplied with v0 and has a meaning of free path of test
!   particle with v0. Respectively, collision frequencies have the meaning of inverse
!   mean free paths.
!
!   Input variables:
!        formal: m1,Z1         - mass and charge number of the colliding particle
!                m             - mass of the ion species (n entries) (where n is the number of particle species 
!                                with which the test particle collides)
!                Z             - charge numbers of these species (n entries)
!                dens          - densities of ion species and electrons (n entries), 1/cm**3
!                temp          - temperatures of ion species and electrons (n entries), eV
!                e0            - test particle energy used for normalisation, eV
!   Output variables:
!        formal: v0            - test particle velocity corresponding to e0, cm/s
!                efcolf        - normalized collision frequencies (n entries)
!                velrat        - ratio of v0 to the background particle thermal velocity $v_{t}=\sqrt(2T/m)$ (n entries)
!                enrat         - ratio of e0 to the background species energy (n entries)

  use constants, only: ev2erg, pi, echarge

  integer :: n,i, i_end
  real(dp), dimension(:) :: m,Z,dens,temp,efcolf,velrat,enrat
  real(dp), dimension(:), allocatable :: lambda
  real(dp) :: m1,Z1,e0, v0, frecol_base

  n = size(temp)
  allocate(lambda(n))

  v0=sqrt(2.d0*e0*ev2erg/m1)
  frecol_base = 2.d0*pi*echarge**4*Z1**2/(m1**2*v0**3)

  do i = 1,n
    call lambda_alpha_beta(Z1, Z(i), m1, m(i), temp(i), temp(i), dens(i), dens(i), lambda(i))
    enrat(i)=e0/temp(i)
    velrat(i)=v0/sqrt(2.d0*temp(i)*ev2erg/m(i))
    efcolf(i)=frecol_base*Z(i)**2*lambda(i)*dens(i)*velrat(i)/v0
  enddo

end subroutine collis_init
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine lambda_alpha_beta(z1, z2, m1, m2, T1, T2, n1, n2, lambda_ab)

  use constants, only: ame,amp

  real(dp), intent(in) :: z1, z2, m1, m2, T1, T2, n1, n2
  real(dp), intent(out) :: lambda_ab
  real(dp) :: ne, ni, zi, Te, Ti, mi
  integer :: species1, species2

  species1 = sign(1,int(z1))
  species2 = sign(1,int(z2))

  if ((species1.eq.1).and.(species2.eq.1)) then
    lambda_ab = 23.d0-log(max(epsilon(1.d0), z1*z2*(m1+m2)  /  (m1*T2+m2*T1) * sqrt(n1*z1**2/T1 + n2*z2**2/T2)))
  elseif ((species1.eq.2).and.(species2.eq.2)) then
    lambda_ab = 23.5d0 - log(max(epsilon(1.d0), sqrt(n2)*T2**(-5.d0/4.d0))) - sqrt(1.d-5 + (log(T2)-2.d0)**2/16.d0)
  else
    if (species1.eq.1) then
      ne = n2
      ni = n1
      zi = z1
      Te = T2
      Ti = T1
      mi = m1
    else
      ne = n1
      ni = n2
      zi = z2
      Te = T1
      Ti = T2
      mi = m2
    endif
    if (Te.lt.Ti*ame/(mi)) then
      lambda_ab = 16.d0-log(max(epsilon(1.d0), sqrt(ni)*Ti**(-1.5d0)*zi**2*mi/amp))
    elseif (Te.lt.10*zi**2) then
      lambda_ab = 23.d0-log(max(epsilon(1.d0), sqrt(ne)*zi*Te**(-1.5d0)))
    else
      lambda_ab = 24.d0-log(max(epsilon(1.d0), sqrt(ne)/Te))
    endif
  endif

end subroutine lambda_alpha_beta
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine stost(efcolf,velrat,enrat,z,dtau,iswmode,ierr,tau,randnum)
!
!  z(1:5)   - phase space coordinates: z(1:3) - spatial position,
!                                      z(4)   - normalized velocity module
!                                      z(5)   - pitch parameter
!  dtau     - normalized time step: dtau=dtau*v0 - has the dimension of length
!  iswmode  - switch of the collision mode:
!               1 - full operator (pitch-angle and energy scattering and drag)
!               2 - energy scattering and drag only
!               3 - drag only
!               4 - pitch-angle scattering only
!  ierr     - error code:
!               0 - good,
!               1 - bad argument (pitch |z(5)| > 1 ),
!               2 - step over pitch exceeds 1 (pitch was
!                   replaced by randomly distributed on [-1,1]),
!               3 - new pitch module exceeds 1, reflection from
!                   the boudary was performed,
!               10 or >10 - new momentum module is less then
!                   prescribed minimum, reflection was performed.
!
  integer :: iswmode,ierr,n
  real(dp), parameter :: pmin=5.e-2
  real(dp) :: dtau,p,dpp,dhh,fpeff,alam,dalam,coala, upper_limit
  real(dp), dimension(5) :: z
  real(dp) :: ur, epsilon, q
  real(dp), dimension(:), intent(in) :: efcolf,velrat,enrat
  real(dp), dimension(:), allocatable :: dpp_vec,dhh_vec,fpeff_vec
  real(dp), optional :: tau
  real(dp), dimension(3), intent(in), optional :: randnum
  real(dp) :: z4_save

  epsilon = 0.1
  q = 0.3
  upper_limit = 30
  n = size(efcolf)
  allocate(dpp_vec(n))
  allocate(dhh_vec(n))
  allocate(fpeff_vec(n))

  p=z(4)
  call coleff(efcolf,velrat,enrat,p,dpp_vec,dhh_vec,fpeff_vec)

  dpp = sum(dpp_vec)
  dhh = sum(dhh_vec)
  fpeff = sum(fpeff_vec)

  ierr=0

  if (present(tau)) then
    dtau = min(epsilon**2/(2*dhh),tau,upper_limit)
    if (z(4).lt.q) then
      !dtau = min(dtau*(q/z(4))**2,tau),upper_limit)
      dtau = min(1.0d-2/(2*dhh),tau,upper_limit)
    endif
  endif

  if(iswmode.eq.1.or.iswmode.eq.4) then
    alam=z(5)
    coala=1.d0-alam**2

    if(coala.lt.0.d0) then
      ierr=1
      return
    endif

    if (present(randnum)) ur = randnum(1)
    if (.not.present(randnum)) call getran(1,ur)

    dalam=sqrt(2.d0*dhh*coala*dtau)*dble(ur)-2.d0*alam*dhh*dtau

    if(abs(dalam).gt.1.d0) then
      ierr=2
      if (present(randnum)) ur = randnum(2)
      if (.not.present(randnum)) call random_number(ur)

      alam=2.d0*(dble(ur)-0.5d0)
    else
      alam=alam+dalam
      if(alam.gt.1.d0) then
        ierr=3
        alam=2.d0-alam
      elseif(alam.lt.-1.d0) then
        ierr=3
        alam=-2.d0-alam
      endif
    endif

    z(5)=alam
    if(iswmode.eq.4) return
  endif

  if(iswmode.lt.3) then

    if (present(randnum)) ur = randnum(3)
    if (.not.present(randnum)) call getran(0,ur)

    z4_save = z(4)
    z(4)=z(4)+sqrt(abs(2.d0*dpp*dtau))*dble(ur)+fpeff*dtau
    if (z(4)/z4_save.gt.10) then
        print*, 'v_old/v0 = ', z4_save, 'v_new/v0 = ', z(4), 'ratio =', z(4)/z4_save  !, 'v_old = ', z4_save*3.508831372d9
    endif
  else
    z(4)=z(4)+fpeff*dtau
  endif

  if(z(4).lt.pmin) then
    ierr=ierr+10
    z(4)=pmin+abs(pmin-z(4))
  endif

  return
end subroutine stost
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine getran(irand,ur)

  ! Produces a random number with zero average and unit variance

  ! Input parameters: irand - 0 for continuous, 1 for +1 -1,
  ! Output parameters: ur   - random number

    integer :: irand
    real(dp) :: ur

    call random_number(ur)

    if(irand.eq.0) then !continuous random number, constant is sqrt(12)
      ur=3.464102*(ur-.5)
    else !discrete random number
      if(ur.gt..5) then
        ur=1.
      else
        ur=-1.
      endif
    endif

end subroutine getran
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
end module collis_ions