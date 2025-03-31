module collis_ions

implicit none
!
contains
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine coleff(efcolf,velrat,enrat,p,dpp_vec,dhh_vec,fpeff_vec)
!
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
!
  integer :: i, n
  double precision, dimension(:), intent(in) :: efcolf,velrat,enrat
  double precision, dimension(3) :: dpp_vec,dhh_vec,fpeff_vec
  double precision :: p,plim,xbeta,dpd
!
  plim=max(p,1.d-8)
  n = size(efcolf)
!
  do i=1,n
    xbeta=p*velrat(i)
!
    call onseff(xbeta,dpp_vec(i),dhh_vec(i),dpd)
!
    fpeff_vec(i) = (dpd/plim-2.0*dpp_vec(i)*p*enrat(i))*efcolf(i)
    dpp_vec(i) = dpp_vec(i)*efcolf(i)
    dhh_vec(i) = dhh_vec(i)*efcolf(i)
  enddo
!
  dhh_vec = dhh_vec/plim**2
!
  return
end subroutine coleff
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine onseff(v,dp,dh,dpd)
!
!  dp - dimensionless dpp
!  dh - dhh*p^2     (p - dmls)
!  dpd - (1/p)(d/dp)p^2*dp   (p - dmls)
!
! square root of pi
  double precision, parameter :: sqp=1.7724538d0
! cons=4./(3.*sqrt(pi))
  double precision, parameter :: cons=.75225278d0
  double precision :: v,dp,dh,dpd,v2,v3,ex,er
!
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
!
  return
end subroutine onseff

subroutine collis_init(am0,Z0,m,Z,dens,temp,eion,v0,efcolf,velrat,enrat,boole_no_electrons)

!   Performs precomputation of the constants for Coulomb collision
!   operator for test particles colliding with n-1 sorts of ions and with electrons
!
!   Normalisation: test particle velocity is normalized by v0,
!   time is multiplied with v0 and has a meaning of free path of test
!   particle with v0. Respectively, collision frequencies have the meaning of inverse
!   mean free paths.
!
!   Input variables:
!        formal: am0,Z0        - mass number and charge number of the colliding particle
!                m             - mass numbers of the ion species (n-1 entries) (where n is the number of particle species 
!                                with which the test particle collides (n-1 ion species and electrons))
!                Z             - charge numbers of these species (n-1 entries)
!                dens          - densities of ion species and electrons (n entries), 1/cm**3
!                temp          - temperatures of ion species and electrons (n entries), eV
!                eion          - test particle energy used for normalisation, eV
!   Output variables:
!        formal: v0            - test particle velocity corresponding to eion, cm/s
!                efcolf        - normalized collision frequencies (n entries)
!                velrat        - ratio of v0 to the background particle thermal velocity $v_{t}=\sqrt(2T/m)$ (n entries)
!                enrat         - ratio of eion to the background species energy (n entries)

  integer :: n,i, i_end
  double precision, dimension(:) :: m,Z,dens,temp,efcolf,velrat,enrat
  double precision, dimension(:), allocatable :: lambda
  double precision :: am0,Z0,eion
  double precision :: v0
  double precision :: pi,pmass,emass,e,ev,frecol_base
  double precision :: k
  logical, intent(in), optional :: boole_no_electrons

  pi=3.14159265358979d0
  pmass=1.6726d-24
  emass=9.1094d-28
  e=4.8032d-10
  ev=1.6022d-12
  k=1.60d-12 !ev/erg

  n = size(temp)
  allocate(lambda(n))

  v0=sqrt(2.d0*eion*ev/(pmass*am0))

  i_end = n-1
  if (present(boole_no_electrons)) then
    if (boole_no_electrons) i_end = n
  endif

  do i = 1,i_end !go through all ion species in the loop and treat electrons afterwards
    enrat(i)=eion/temp(i)
    velrat(i)=v0/sqrt(2.d0*temp(i)*ev/(pmass*m(i)))
    lambda(i)=23.d0-log(max(epsilon(1.d0), &
          sqrt(dens(i)*Z(i)**2/temp(i))*Z0*Z(i)*(am0+m(i))/(am0*temp(i)+m(i)*eion)))
  enddo

  frecol_base=2.d0*pi*dens(i_end)*e**4*Z0**2/((am0*pmass)**2*v0**3) !usual
  frecol_base=frecol_base/v0                                  !normalized

  do i=1,i_end
    efcolf(i)=frecol_base*Z(i)**2*lambda(i)*dens(i)/dens(i_end)
  enddo

  if (i_end.eq.n-1) then
    enrat(n)=eion/temp(n)
    velrat(n)=v0/sqrt(2.d0*temp(n)*ev/emass)
    dens(n)=sum(dens(1:n-1)*Z)
    if (temp(n).lt.eion*emass/(am0*pmass)) then
      lambda(n)=16.d0-log(sqrt(dens(n))*eion**(-1.5)*Z0**2*am0)
    elseif (temp(n).lt.10*Z0**2) then
      lambda(n)=23.d0-log(sqrt(dens(n))*Z0*temp(n)**(-1.5))
    else
      lambda(n)=24.d0-log(sqrt(dens(n))/temp(n))
    endif

    efcolf(n)=frecol_base*lambda(n)
  endif

  efcolf=efcolf*velrat

end subroutine collis_init
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
  double precision, parameter :: pmin=1.e-8
  double precision :: dtau,p,dpp,dhh,fpeff,alam,dalam,coala, upper_limit
  double precision, dimension(5) :: z
  double precision :: ur, epsilon, q
  double precision, dimension(:), intent(in) :: efcolf,velrat,enrat
  double precision, dimension(:), allocatable :: dpp_vec,dhh_vec,fpeff_vec
  double precision, optional :: tau
  double precision, dimension(3), intent(in), optional :: randnum
!
  epsilon = 0.1
  q = 0.3
  upper_limit = 30
  n = size(efcolf)
  allocate(dpp_vec(n))
  allocate(dhh_vec(n))
  allocate(fpeff_vec(n))
!
  p=z(4)
  call coleff(efcolf,velrat,enrat,p,dpp_vec,dhh_vec,fpeff_vec)

  dpp = sum(dpp_vec)
  dhh = sum(dhh_vec)
  fpeff = sum(fpeff_vec)
!
  ierr=0
!
  if (present(tau)) then
    dtau = min(epsilon**2/(2*dhh),tau,upper_limit)
    if (z(4).lt.q) then
      dtau = min(dtau*(q/z(4))**2,tau)
    endif
  endif
!
  if(iswmode.eq.1.or.iswmode.eq.4) then
    alam=z(5)
    coala=1.d0-alam**2
!
    if(coala.lt.0.d0) then
      ierr=1
      return
    endif
!  
    if (present(randnum)) ur = randnum(1)
    if (.not.present(randnum)) call getran(1,ur)
!
    dalam=sqrt(2.d0*dhh*coala*dtau)*dble(ur)-2.d0*alam*dhh*dtau
!
    if(abs(dalam).gt.1.d0) then
      ierr=2
      if (present(randnum)) ur = randnum(2)
      if (.not.present(randnum)) call random_number(ur)
!
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
!
    z(5)=alam
    if(iswmode.eq.4) return
  endif
!
  if(iswmode.lt.3) then
!
  if (present(randnum)) ur = randnum(3)
  if (.not.present(randnum)) call getran(0,ur)
!
    z(4)=z(4)+sqrt(abs(2.d0*dpp*dtau))*dble(ur)+fpeff*dtau
  else
    z(4)=z(4)+fpeff*dtau
  endif
!
  if(z(4).lt.pmin) then
    ierr=ierr+10
    z(4)=pmin+abs(pmin-z(4))
  endif
!
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
    double precision :: ur

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
    return
end subroutine getran
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
end module collis_ions