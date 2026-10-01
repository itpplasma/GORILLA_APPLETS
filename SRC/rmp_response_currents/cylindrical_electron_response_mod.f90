module cylindrical_electron_response_mod
    ! Linear, zero-FLR electrons for the exact KIM cylinder (OU model 0).
    ! The perpendicular Maxwellian source energy is analytically averaged.
    ! Independent helical potential and radial magnetic sources are retained.
    use, intrinsic :: iso_fortran_env, only: real64
    use cylindrical_characteristics_mod, only: cylindrical_ou_resolvent
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: assemble_cylindrical_electrons
contains
    subroutine assemble_cylindrical_electrons(background,M,L,clight,blocks,refinement)
        real(real64),intent(in) :: background(:,:),L,clight
        integer,intent(in) :: M,refinement
        complex(real64),allocatable,intent(out) :: blocks(:,:,:)
        complex(real64),parameter :: ii=(0.0_real64,1.0_real64)
        real(real64),parameter :: pi=acos(-1.0_real64)
        real(real64) :: r,ks,kpar,omega_e,lambda_d,nu,vt,wc,a1,a2,pref,kr,krp
        complex(real64) :: moments(0:3,0:3),local(4),phase,phi_pref,b_pref
        integer :: j,row,col,N,dim,info
        N=size(background,1)
        dim=2*M+1
        if (size(background,2)/=13 .or. M<0 .or. N<dim .or. &
            L<=0.0_real64 .or. clight<=0.0_real64) error stop 'invalid cylindrical background'
        if (.not.all(ieee_is_finite(background))) error stop 'nonfinite cylindrical background'
        allocate(blocks(dim,dim,4))
        blocks=0.0_real64
        do j=1,N
            r=background(j,1)
            ks=background(j,2)
            kpar=background(j,3)
            omega_e=background(j,4)
            lambda_d=background(j,5)
            nu=background(j,6)
            vt=background(j,7)
            wc=background(j,8)
            a1=background(j,9)
            a2=background(j,10)
            if (lambda_d<=0.0_real64 .or. nu<=0.0_real64 .or. vt<=0.0_real64 .or. &
                wc==0.0_real64) error stop 'invalid electron thermal/collision parameters'
            call cylindrical_ou_resolvent(kpar*vt/nu,-omega_e/nu,moments,info,refinement)
            if (info/=0) error stop 'cylindrical characteristic resolvent failed or singular'
            ! q^2*n/T = 1/(4*pi*lambda_D^2), vT^2/omega_c = c*T/(q*B).
            ! Charge sources: q*n*(i*c*ks*Phi/B-vT*u*Br/B)*[A1+A2*(1+u^2/2)].
            ! The explicit adiabatic charge -q^2*n*Phi/T is included exactly once.
            pref=1.0_real64/(4.0_real64*pi*lambda_d**2)
            phi_pref=ii*pref*vt**2*ks/(wc*nu)
            b_pref=-pref*vt**3/(wc*nu*clight)
            local(1)=-pref+phi_pref*((a1+a2)*moments(0,0)+0.5_real64*a2*moments(0,2))
            local(2)=b_pref*((a1+a2)*moments(0,1)+0.5_real64*a2*moments(0,3))
            local(3)=vt*phi_pref*((a1+a2)*moments(1,0)+0.5_real64*a2*moments(1,2))
            local(4)=vt*b_pref*((a1+a2)*moments(1,1)+0.5_real64*a2*moments(1,3))
            do col=1,dim
                krp=2.0_real64*pi*real(col-M-1,real64)/L
                do row=1,dim
                    kr=2.0_real64*pi*real(row-M-1,real64)/L
                    phase=exp(ii*(krp-kr)*r)/real(N,real64)
                    blocks(row,col,:)=blocks(row,col,:)+phase*local
                end do
            end do
        end do
    end subroutine
end module
