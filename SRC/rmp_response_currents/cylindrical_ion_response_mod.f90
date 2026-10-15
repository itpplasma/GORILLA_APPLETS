module cylindrical_ion_response_mod
    ! Static m_phi=0 ions: equation (2.38) at source AND particle observation.
    ! w=v_perp**2/(2*vT**2); the perpendicular Maxwellian measure is exp(-w) dw.
    use, intrinsic :: iso_fortran_env, only: real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use cylindrical_characteristics_mod, only: cylindrical_ou_resolvent
    use cylindrical_electron_response_mod, only: restore_cylindrical_energy
    implicit none
    private
    public :: assemble_cylindrical_ions, ion_gyroharmonic, perpendicular_rule
contains
    pure real(real64) function ion_gyroharmonic(kperp, rho_thermal, w) result(harmonic)
        ! rho_thermal=vT/abs(omega_ci) [cm], kperp [cm^-1]; full integral FLR.
        real(real64), intent(in) :: kperp, rho_thermal, w
        harmonic=bessel_j0(sqrt(2.0_real64*w)*rho_thermal*kperp)
    end function

    subroutine perpendicular_rule(order, nodes, weights)
        ! Golub-Welsch Gauss-Laguerre rule for integral_0^infinity exp(-w) f(w) dw.
        integer, intent(in) :: order
        real(real64), allocatable, intent(out) :: nodes(:), weights(:)
        real(real64), allocatable :: off(:), eigenvectors(:,:), work(:)
        integer :: j, info
        if(order<2) error stop 'perpendicular quadrature needs at least two nodes'
        allocate(nodes(order),weights(order),off(order-1),eigenvectors(order,order),work(2*order-2))
        nodes=real([(2*j-1,j=1,order)],real64)
        off=real([(j,j=1,order-1)],real64)
        call dstev('V',order,nodes,off,eigenvectors,order,work,info)
        if(info/=0) error stop 'Gauss-Laguerre eigenproblem failed'
        weights=eigenvectors(1,:)**2
    end subroutine

    subroutine assemble_cylindrical_ions(background,M,L,clight,blocks,refinement,model,n_perp)
        real(real64), intent(in) :: background(:,:),L,clight
        integer, intent(in) :: M,refinement,model,n_perp
        complex(real64), allocatable, intent(out) :: blocks(:,:,:)
        real(real64), parameter :: pi=acos(-1.0_real64)
        complex(real64), parameter :: ii=(0.0_real64,1.0_real64)
        real(real64), allocatable :: nodes(:),weights(:),wave(:),gyro(:,:),s0(:,:),s1(:,:)
        real(real64) :: ks,kpar,omega,nu,vt,wc,a1,a2,pref,rho_thermal,force
        complex(real64) :: moments(0:3,0:3),phi_pref,b_pref,phase,local(4)
        integer :: N,dim,j,row,col,k,info
        N=size(background,1);dim=2*M+1
        if(size(background,2)/=14.or.M<0.or.N<dim.or.L<=0.or.clight<=0.or.refinement<1) &
            error stop 'invalid ion background/grid/refinement'
        if(.not.all(ieee_is_finite(background)).or..not.all(ieee_is_finite([L,clight]))) &
            error stop 'nonfinite ion background'
        if(model/=0.and.model/=1) error stop 'ion provider supports OU models 0/1'
        call perpendicular_rule(n_perp,nodes,weights)
        allocate(blocks(dim,dim,4),wave(dim),gyro(n_perp,dim),s0(dim,dim),s1(dim,dim))
        blocks=0
        wave=2*pi*real([(j-M-1,j=1,dim)],real64)/L
        do j=1,N
            ks=background(j,2);kpar=background(j,3);omega=background(j,4)
            nu=background(j,6);vt=background(j,7);wc=background(j,8)
            a1=background(j,9);a2=background(j,10);rho_thermal=background(j,14)
            if(background(j,5)<=0.or.nu<=0.or.vt<=0.or.wc==0.or.rho_thermal<0) &
                error stop 'invalid ion thermal/collision parameters'
            call cylindrical_ou_resolvent(kpar*vt/nu,-omega/nu,moments,info,refinement)
            if(info/=0) error stop 'ion characteristic resolvent failed or singular'
            if(model==1) then
                call restore_cylindrical_energy(moments,kpar*vt/nu,-omega/nu,info)
                if(info/=0) error stop 'ion energy-restoring resolvent failed or singular'
            end if
            do col=1,dim
                do k=1,n_perp
                    gyro(k,col)=ion_gyroharmonic(sqrt(ks**2+wave(col)**2),rho_thermal,nodes(k))
                end do
            end do
            ! Separate source-column and observation-row gyroaverages.
            ! S0=<J0_obs J0_src>, S1=<w J0_obs J0_src> over perpendicular speed.
            s0=matmul(transpose(gyro),gyro*spread(weights,2,dim))
            s1=matmul(transpose(gyro),gyro*spread(weights*nodes,2,dim))
            ! P=q_i^2*n_i/T_i [statC/(statV cm^3)]; omega_ci retains the charge sign.
            pref=1/(4*pi*background(j,5)**2)
            phi_pref=ii*pref*vt**2*ks/(wc*nu)
            b_pref=-pref*vt**3/(wc*nu*clight)
            do col=1,dim
                do row=1,dim
                    force=a1*s0(row,col)+a2*s1(row,col)
                    ! KIM's periodic model uses the local particle Boltzmann term -P.
                    ! Do not multiply that term by the source/observer Bessel product.
                    local(1)=-pref+phi_pref*(force*moments(0,0) &
                        +0.5_real64*a2*s0(row,col)*moments(0,2))
                    local(2)=b_pref*(force*moments(0,1) &
                        +0.5_real64*a2*s0(row,col)*moments(0,3))
                    local(3)=vt*phi_pref*(force*moments(1,0) &
                        +0.5_real64*a2*s0(row,col)*moments(1,2))
                    local(4)=vt*b_pref*(force*moments(1,1) &
                        +0.5_real64*a2*s0(row,col)*moments(1,3))
                    phase=exp(ii*(wave(col)-wave(row))*background(j,1))/N
                    blocks(row,col,:)=blocks(row,col,:)+phase*local
                end do
            end do
        end do
        if(.not.all(ieee_is_finite(real(blocks))).or..not.all(ieee_is_finite(aimag(blocks)))) &
            error stop 'nonfinite ion response'
    end subroutine
end module
