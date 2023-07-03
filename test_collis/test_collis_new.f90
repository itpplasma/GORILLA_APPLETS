!to compile, write "gfortran test_collis_new.f90 -fopenmp -o new_test"
program test_collis
  use collis_ions
  use omp_lib, only: omp_get_thread_num, omp_get_num_threads, omp_set_num_threads

  implicit none

  integer, parameter :: n=3 !number of background particle species
  double precision, dimension(5) :: z
  double precision :: eion, eV, p_mass,nu_perp
  double precision, dimension(n) :: dens,temp,efcolf,velrat,enrat
  double precision, dimension(n-1) :: mass_num,charge_num
  double precision :: v0,dtau
  double precision :: m0=4.0d0, z0=2.0d0  ! Test particle mass and charge no.
  double precision :: trace_time = 1.0d-1   ! Total trace time in seconds

  integer, parameter :: num_t_step=100000 
  integer :: i,j,k,ierr,num_particles
  double precision, dimension(num_t_step) :: energy_of_t

  p_mass=1.6726d-24
  ev=1.6022d-12

  mass_num = (/2.0d0,3.0d0/)
  charge_num = (/1.0d0,1.0d0/)
  dens = (/0.5d14,0.5d14,0.d0/)
  temp = (/1.0d4,1.0d4,1.0d4/)
  eion = 3.5d6

  call collis_init(m0,z0,mass_num,charge_num,dens,temp,eion,v0,nu_perp,efcolf,velrat,enrat)

  trace_time=1.d0

  print *, 'v0 = ', v0, ' = ', sqrt(2.d0*eion*ev/(m0*p_mass))

  dtau = trace_time*v0/dble(num_t_step-1)

  print *, 'dtau = ', dtau

  energy_of_t=0.d0
  num_particles=10000
  k = 0

  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP& SHARED(k,num_particles,efcolf,velrat,enrat,dtau) &
  !$OMP& PRIVATE(i,z,j,ierr) &
  !$OMP& REDUCTION(+:energy_of_t)
  print*, 'get number of threads', omp_get_num_threads()
  !$OMP DO
  do i=1,num_particles
    z = 0.0d0
    z(4) = 1.0d0
    energy_of_t(1)=energy_of_t(1)+1.d0
    do j = 2, num_t_step
      call stost(efcolf,velrat,enrat,z,dtau,1,ierr)
  !    if(ierr.ne.0) print *,ierr
      energy_of_t(j)=energy_of_t(j)+z(4)**2
    enddo
    !$omp critical
      k = k+1
      if (modulo(k,int(num_particles/100)).eq.0) print *,k,' of ',num_particles,' completed'
    !$omp end critical
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

energy_of_t=energy_of_t/dble(num_particles)
  open(1, file='energy_average_new.dat')
do i = 1, num_t_step
    write (1,*) dble(i-1)*dtau/v0, energy_of_t(i)
enddo
  close(1)
  print*, 'End of test of new version of collision operator'
end program test_collis

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
    double precision :: p,plim,xbeta,dp,dh,dpd,dpp,dhh,fpeff
  !
    plim=max(p,1.d-8)
    n = size(efcolf)
  !
    dpp=0.0d0
    dhh=0.0d0
    fpeff=0.0d0
  !
    do i=1,n
      xbeta=p*velrat(i)
  !
      call onseff(xbeta,dp,dh,dpd)
  !
      dpp=dpp+dp*efcolf(i)
      dhh=dhh+dh*efcolf(i)
      fpeff=fpeff+(dpd/plim-2.0*dp*p*enrat(i))*efcolf(i)
  !
      call onseff(xbeta,dpp_vec(i),dhh_vec(i),dpd)
  !
      fpeff_vec(i) = (dpd/plim-2.0*dpp_vec(i)*p*enrat(i))*efcolf(i)
      dpp_vec(i) = dpp_vec(i)*efcolf(i)
      dhh_vec(i) = dhh_vec(i)*efcolf(i)
    enddo
  !
    dhh=dhh/plim**2
    dhh_vec = dhh_vec/plim**2
  !
    return
    end
  !
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
    subroutine onseff(v,dp,dh,dpd)
  !
  !  dp - dimensionless dpp
  !  dh - dhh*p^2     (p - dmls)
  !  dpd - (1/p)(d/dp)p^2*dp   (p - dmls)
  !
    implicit none
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
    end
  !
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !      FUNCTION ERF(X)
  !      PARAMETER  ( A1 = 0.07052 30784, A2 = 0.04228 20123,
  !     ,             A3 = 0.00927 05272, A4 = 0.00015 10143,
  !     ,             A5 = 0.00027 65672, A6 = 0.00004 30638 )
  !      F(T) = 1./((1.+T*(A1+T*(A2+T*(A3+T*(A4+T*(A5+T*A6))))))**4)**4
  !      W = 1. - F(ABS(X))
  !      ERF = SIGN(W,X)
  !      RETURN
  !      END
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
    subroutine collis_init(am0,Z0,m,Z,dens,temp,eion,v0,nu_perp0,efcolf,velrat,enrat)
  !
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
  !
    integer :: n,i
    double precision, dimension(:) :: m,Z,dens,temp,efcolf,velrat,enrat
    double precision, dimension(:), allocatable :: lambda
    double precision :: am0,Z0,eion
    double precision :: v0
    double precision :: pi,pmass,emass,e,ev,frecol_base
    double precision :: x1,x2,xe,k,nu_perp0
  !
    pi=3.14159265358979d0
    pmass=1.6726d-24
    emass=9.1094d-28
    e=4.8032d-10
    ev=1.6022d-12
    k=1.60d-12 !ev/erg
  !
    n = size(temp)
    allocate(lambda(n))
  !
    v0=sqrt(2.d0*eion*ev/(pmass*am0))
  !
    do i = 1,n-1 !go through all ion species in the loop and treat electrons afterwards
      enrat(i)=eion/temp(i)
      velrat(i)=v0/sqrt(2.d0*temp(i)*ev/(pmass*m(i)))
      lambda(i)=23.d0-log(max(epsilon(1.d0), &
            sqrt(dens(i)*Z(i)**2/temp(i))*Z0*Z(i)*(am0+m(i))/(am0*temp(i)+m(i)*eion)))
    enddo
  !
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
  !
    frecol_base=2.d0*pi*dens(n)*e**4*Z0**2/((am0*pmass)**2*v0**3) !usual
    frecol_base=frecol_base/v0                                  !normalized
  !
    do i=1,n-1
      efcolf(i)=frecol_base*Z(i)**2*lambda(i)*dens(i)/dens(n)
    enddo
    efcolf(n)=frecol_base*lambda(n)
  !
    efcolf=efcolf*velrat
  !
    x1 = m(1)*pmass*v0**2/(2*k*temp(1))
    x2 = m(2)*pmass*v0**2/(2*k*temp(2))
    xe =     emass*v0**2/(2*k*temp(3))
    !since x1, x2 and xe are all much smaller than one (verify this within the subroutine), the approximation for fast particles 
    !can be taken in nrl page 32
    !(2018 edition) to evaluate nu_perp (do this in a more refined way, e.g. by using the exact formula on page 31 later)
    nu_perp0 = Z0**2*sum(Z(1:n-1)**2*lambda(1:n-1)*dens(1:n-1))*1.8d-7*am0**(-0.5)*(ev*eion)**(-1.5)
  !
    end
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
    subroutine stost(efcolf,velrat,enrat,z,dtauc,iswmode,ierr)
  !
  !  z(1:5)   - phase space coordinates: z(1:3) - spatial position,
  !                                      z(4)   - normalized velocity module
  !                                      z(5)   - pitch parameter
  !  dtauc    - normalized time step: dtauc=dtau*v0 - has the dimension of length
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
    double precision :: dtauc,p,dpp,dhh,fpeff,alam,dalam,coala
    double precision, dimension(5) :: z
    real :: ur
    double precision, dimension(:) :: efcolf,velrat,enrat
    double precision, dimension(:), allocatable :: dpp_vec,dhh_vec,fpeff_vec
  !
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
    if(iswmode.eq.1.or.iswmode.eq.4) then
      alam=z(5)
      coala=1.d0-alam**2
  !
      if(coala.lt.0.d0) then
        ierr=1
        return
      endif
  !  
      call getran(1,ur)
  !
      dalam=sqrt(2.d0*dhh*coala*dtauc)*dble(ur)-2.d0*alam*dhh*dtauc
  !
      if(abs(dalam).gt.1.d0) then
        ierr=2
  !
      call random_number(ur)   
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
    call getran(0,ur)
  !
      z(4)=z(4)+sqrt(abs(2.d0*dpp*dtauc))*dble(ur)+fpeff*dtauc
    else
      z(4)=z(4)+fpeff*dtauc
    endif
  !
    if(z(4).lt.pmin) then
      ierr=ierr+10
      z(4)=pmin+abs(pmin-z(4))
    endif
  !
    return
    end
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
end module collis_ions
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine getran(irand,ur)
!
!  Produces the random number with zero average and unit variance
!
!  Input parameters: irand - 0 for continious, 1 for +1 -1,
!  Output parameters: ur   - random number
!
  call random_number(ur)
!
  if(irand.eq.0) then
!
!  continiuos random number, constant is sqrt(12)
!
    ur=3.464102*(ur-.5)
  else
!
!  discrete random number
!
    if(ur.gt..5) then
      ur=1.
    else
      ur=-1.
    endif
  endif
  return
  end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc