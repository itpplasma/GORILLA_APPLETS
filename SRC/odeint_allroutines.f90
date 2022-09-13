module odeint_mod
   
  private

  integer :: kmax=0, kount=0, kmaxx=200, ialloc
  double precision :: dxsav=0.d0
  double precision, dimension(:),   allocatable :: dydx,xp,y,yscal
  double precision, dimension(:,:), allocatable :: yp
  double precision, dimension(:),   allocatable :: ak2,ak3,ak4,ak5
  double precision, dimension(:),   allocatable :: ak6,ytemp
  double precision, dimension(:),   allocatable :: yerr,ytemp1
  
  !$OMP THREADPRIVATE(kmax,kount,kmaxx,ialloc,dxsav,dydx,xp,y,yscal,&
  !$OMP& yp,ak2,ak3,ak4,ak5,ak6,ytemp,yerr,ytemp1)
  
  public :: odeint_allroutines

contains

  subroutine odeint_allroutines(y,nvar,x1,x2,eps,derivs)
 
    implicit none

    external :: derivs
    integer :: nvar,nok,nbad
    double precision :: x1,x2,eps,h1,hmin
    double precision, dimension(nvar) :: y

    h1=x2-x1
    hmin=0.d0

    call odeint(y,nvar,x1,x2,eps,h1,hmin,nok,nbad,derivs)

  end subroutine odeint_allroutines

  subroutine alloc_odeint(nvar)

    implicit none 
    
    integer, intent(in) :: nvar
    
    if (ialloc.eq.1) then
      allocate(dydx(nvar),xp(kmaxx),y(nvar))
      allocate(yp(nvar,kmaxx),yscal(nvar))
      allocate(ak2(nvar),ak3(nvar),ak4(nvar),ak5(nvar))
      allocate(ak6(nvar),ytemp(nvar))
      allocate(yerr(nvar),ytemp1(nvar))
    else
      deallocate(dydx,xp,y,yp,yscal)
      deallocate(ak2,ak3,ak4,ak5,ak6,ytemp)
      deallocate(yerr,ytemp1)
    end if

  end subroutine alloc_odeint

  subroutine odeint(ystart, nvar, x1, x2, eps, h1, hmin, nok, nbad, &
    & derivs)

    implicit none

    INTEGER nbad,nok,nvar,KMAXX,MAXSTP
    double precision eps,h1,hmin,x1,x2,ystart(nvar),TINY
    EXTERNAL derivs
    PARAMETER (MAXSTP=1000000,TINY=1.e-30)
    INTEGER i,nstp
    double precision h,hdid,hnext,x,xsav

    xsav = 1.234e5

    ialloc = 1
    call alloc_odeint(nvar)

    x = x1
    h = sign(h1,x2-x1)
    nok = 0
    nbad = 0
    kount = 0

    do i=1,nvar
      y(i) = ystart(i)
    end do

    if (kmax.gt.0) xsav = x-2.*dxsav
    do nstp=1,MAXSTP
      call derivs(x,y,dydx)
      do i=1,nvar
        yscal(i) = abs(y(i))+abs(h*dydx(i))+TINY
      end do

      if (kmax .gt. 0) then
        if (abs(x-xsav) .gt. abs(dxsav)) then
          if (kount .lt. kmax-1) then
            kount = kount+1
            xp(kount) = x
            do i=1,nvar
              yp(i,kount)=y(i)
            end do
            xsav = x
          end if
        end if
      end if
      if ((x+h-x2)*(x+h-x1) .gt. 0.0) h=x2-x
      call rkqs(y,dydx,nvar,x,h,eps,yscal,hdid,hnext,derivs)
      if (hdid .eq. h) then
        nok = nok+1
      else
        nbad = nbad+1
      end if
      if ((x-x2)*(x2-x1) .ge. 0.) then
        do i=1,nvar
          ystart(i) = y(i)
        end do
        if (kmax .ne. 0) then
          kount = kount+1
          xp(kount) = x
          do i=1,nvar
            yp(i,kount) = y(i)
          end do
        end if
        ialloc=0
        call alloc_odeint(nvar)
        return
      end if
      if (abs(hnext) .lt. hmin) write(*,*) 'stepsize smaller than minimum in odeint'
      h=hnext
    end do

    write(*,*) 'too many steps in odeint'
    ialloc=0
    call alloc_odeint(nvar)
  end subroutine odeint

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine rkck(y, dydx, n, x, h, yout, yerr, derivs)
  
    implicit none

    EXTERNAL derivs

    INTEGER i, n
    double precision h, x, dydx(n), y(n), yerr(n), yout(n)

    double precision, parameter :: A2=0.2d0, A3=0.3d0, A4=0.6d0, A5=1.d0, &
      & A6 = 0.875d0, &
      & B21 = 0.2d0, B31 = 3.d0/40.d0, &
      & B32 = 9.d0/40.d0, B41 = 0.3d0, B42 = -0.9d0, B43 = 1.2d0, &
      & B51 = -11.d0/54.d0, B52 = 2.5d0, &
      & B53 = -70.d0/27.d0, B54 = 35.d0/27.d0, B61 = 1631.d0/55296.d0, &
      & B62 = 175.d0/512.d0, &
      & B63 = 575.d0/13824.d0, B64 = 44275.d0/110592.d0, B65 = 253.d0/4096.d0, &
      & C1 = 37.d0/378.d0, &
      & C3 = 250.d0/621.d0, C4 = 125.d0/594.d0, C6 = 512.d0/1771.d0, &
      & DC1 = C1-2825.d0/27648.d0, &
      & DC3 = C3-18575.d0/48384.d0, DC4 = C4-13525.d0/55296.d0, &
      & DC5 = -277.d0/14336.d0, &
      & DC6 = C6-0.25d0
    
    do i=1,n
      ytemp(i) = y(i) + B21*h*dydx(i)
    end do
   
    call derivs(x+A2*h, ytemp, ak2)
 
    do i=1,n
      ytemp(i) = y(i) + h*(B31*dydx(i) + B32*ak2(i))
    end do
   
    call derivs(x+A3*h, ytemp, ak3)

    do i=1,n
      ytemp(i) = y(i) + h*(B41*dydx(i) + B42*ak2(i) + B43*ak3(i))
    end do

    call derivs(x+A4*h, ytemp, ak4)

    do i=1,n
      ytemp(i) = y(i) + h*(B51*dydx(i) + B52*ak2(i) + B53*ak3(i) + B54*ak4(i))
    end do 
    
    call derivs(x+A5*h, ytemp, ak5)

    do i=1,n
      ytemp(i) = y(i) + h*(B61*dydx(i) + B62*ak2(i) + B63*ak3(i) &
        & + B64*ak4(i) + B65*ak5(i))
    end do

    call derivs(x+A6*h, ytemp, ak6)

    do i=1,n
        yout(i) = y(i) + h*(C1*dydx(i) + C3*ak3(i) + C4*ak4(i) + C6*ak6(i))
    end do

    do i=1,n
        yerr(i) = h*(DC1*dydx(i) + DC3*ak3(i) + DC4*ak4(i) + DC5*ak5(i) &
          & + DC6*ak6(i))
    end do

  end subroutine rkck

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine rkqs(y, dydx, n, x, htry, eps, yscal, hdid, hnext, derivs)

    implicit none

    external derivs

    integer i, n , k
    double precision eps,hdid,hnext,htry,x,dydx(n),y(n),yscal(n)
    double precision errmax,h,htemp,xnew

    double precision, parameter :: SAFETY=0.9d0, PGROW=-.2d0, &
      & PSHRNK=-.25d0, ERRCON=1.89d-4

    h = htry

    do

      call rkck(y, dydx, n, x, h, ytemp1, yerr, derivs)

      errmax=0.0
      do i=1,n
        errmax = max(errmax, abs(yerr(i)/yscal(i)))
      end do
      errmax=errmax/eps

      if (errmax .le. 1.0) exit

      htemp = SAFETY*h*(errmax**PSHRNK)
      h=sign(max(abs(htemp),0.1*abs(h)),h)
      xnew=x+h

      if (xnew.eq.x) then 
        print *,'stepsize underflow in rkqs, x,y = ', x, y
        stop
      end if

    end do

    if (errmax .gt. ERRCON) then
      hnext = SAFETY*h*(errmax**PGROW)
    else
      hnext = 5.0*h
    end if
    hdid = h
    x = x+h
    do i=1,n
      y(i)=ytemp1(i)
    end do

  end subroutine rkqs

end module odeint_mod
