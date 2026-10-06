module cylindrical_electron_response_mod
    ! Linear, zero-FLR electrons for the exact KIM cylinder (OU models 0/1).
    ! The perpendicular Maxwellian source energy is analytically averaged.
    ! Independent helical potential and radial magnetic sources are retained.
    use, intrinsic :: iso_fortran_env, only: real64
    use cylindrical_characteristics_mod, only: cylindrical_ou_resolvent, restore_ou_energy
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: assemble_cylindrical_electrons, close_energy_response, restore_cylindrical_energy
    public :: damped_ou_moments,assemble_diffused_cylinder,assemble_adiabatic_charge
contains
    subroutine damped_ou_moments(x1,x2,damping,moments,info,energy_moments)
        ! Frozen Fourier diffusion adds damping*I to the OU Hermite operator.
        ! This bare reference is closed only by the existing model1 field-space
        ! energy correction; it introduces no additional collision closure.
        real(real64),intent(in) :: x1,x2,damping
        complex(real64),intent(out) :: moments(0:3,0:3)
        integer,intent(out) :: info
        complex(real64),intent(out),optional :: energy_moments(0:4)
        complex(real64),allocatable :: dl(:),du(:),diagonal(:),rhs(:,:)
        complex(real64) :: current(0:3,0:3),previous(0:3,0:3),energy_current(0:4),energy_previous(0:4)
        real(real64) :: coeff(4,0:3),err,scale
        complex(real64),parameter :: ii=(0.0_real64,1.0_real64)
        integer :: order,j,a,b,lapack_info
        info=1;moments=0
        if(.not.all(ieee_is_finite([x1,x2,damping])).or.damping<0) return
        coeff=0;coeff(1,0)=1;coeff(2,1)=1;coeff(1,2)=1
        coeff(3,2)=sqrt(2.0_real64);coeff(2,3)=3;coeff(4,3)=sqrt(6.0_real64)
        previous=0;energy_previous=0;order=64
        do while(order<=8192)
            allocate(dl(order-1),du(order-1),diagonal(order),rhs(order,5))
            diagonal=cmplx(real([(j-1,j=1,order)],real64)+damping,-x2,real64)
            dl=ii*x1*sqrt(real([(j,j=1,order-1)],real64));du=dl
            rhs=0;rhs(:4,:4)=coeff;rhs(3,5)=sqrt(2.0_real64)
            call zgtsv(order,5,dl,diagonal,du,rhs,order,lapack_info)
            if(lapack_info/=0) return
            do b=0,3
                do a=0,3
                    current(a,b)=sum(coeff(:,a)*rhs(:4,b+1))
                end do
            end do
            ! Solve E=u^2-1 directly instead of subtracting nearly equal raw
            ! moments. The Hermite resolvent is transpose-symmetric (no complex
            ! conjugation), so these also give the energy-observer moments.
            do a=0,3
                energy_current(a)=sum(coeff(:,a)*rhs(:4,5))
            end do
            energy_current(4)=sqrt(2.0_real64)*rhs(3,5)
            scale=sqrt(sum(abs(current)**2))
            err=sqrt(sum(abs(current-previous)**2))/max(scale,tiny(1.0_real64))
            err=max(err,sqrt(sum(abs(energy_current-energy_previous)**2)) &
                /max(sqrt(sum(abs(energy_current)**2)),tiny(1.0_real64)))
            deallocate(dl,du,diagonal,rhs)
            ! Existing model1 direct-resolvent accuracy target and cutoff cap.
            if(order>64.and.err<1.0e-12_real64) then
                moments=current
                if(present(energy_moments)) energy_moments=energy_current
                info=0;return
            end if
            previous=current;energy_previous=energy_current;order=2*order
        end do
    end subroutine
    subroutine assemble_adiabatic_charge(background,M,L,channels,adiabatic)
        ! The -1/(4*pi*lambda_D^2) charge term carried only by the analytical
        ! reference: a weighted control variate must restore it when beta<1.
        real(real64),intent(in) :: background(:,:),L
        integer,intent(in) :: M,channels
        complex(real64),allocatable,intent(out) :: adiabatic(:,:,:)
        complex(real64),parameter :: ii=(0.0_real64,1.0_real64)
        real(real64),parameter :: pi=acos(-1.0_real64)
        real(real64),allocatable :: wave(:)
        integer :: dim,N,j,col
        N=size(background,1);dim=2*M+1
        allocate(adiabatic(dim,dim,channels),wave(dim));adiabatic=0
        wave=2*pi*real([(j-M-1,j=1,dim)],real64)/L
        do j=1,N
            do col=1,dim
                adiabatic(:,col,1)=adiabatic(:,col,1)-exp(-ii*wave*background(j,1))/N &
                    *exp(ii*wave(col)*background(j,1))/(4*pi*background(j,5)**2)
            end do
        end do
    end subroutine
    subroutine assemble_diffused_cylinder(background,M,L,clight,diffusion,channels,reference)
        ! Analytical expectation of a frozen, flat cylindrical diffusion pair.
        ! gamma_p=D*(k_p^2+ks^2); radial/helical Brownian coordinates have zero
        ! covariance. Background is frozen at each observation radius, so this
        ! is a numerical control variate, not a new inhomogeneous field solver.
        real(real64),intent(in) :: background(:,:),L,clight,diffusion
        integer,intent(in) :: M,channels
        complex(real64),allocatable,intent(out) :: reference(:,:,:)
        complex(real64),parameter :: ii=(0.0_real64,1.0_real64)
        real(real64),parameter :: pi=acos(-1.0_real64)
        complex(real64) :: moments(0:3,0:3),energy_moments(0:4),aux(9),phi_pref,b_pref,source_phase
        complex(real64),allocatable :: projection(:)
        real(real64),allocatable :: wave(:)
        real(real64) :: r,ks,kpar,omega,nu,vt,wc,a1,a2,pref,qn,damping
        integer :: dim,N,j,col,k,info
        if(diffusion<=0.or.(channels/=4.and.channels/=9)) error stop 'invalid diffused reference'
        N=size(background,1);dim=2*M+1
        allocate(reference(dim,dim,channels),projection(dim),wave(dim));reference=0
        wave=2*pi*real([(j-M-1,j=1,dim)],real64)/L
        do j=1,N
            r=background(j,1);ks=background(j,2);kpar=background(j,3);omega=background(j,4)
            nu=background(j,6);vt=background(j,7);wc=background(j,8)
            a1=background(j,9);a2=background(j,10)
            pref=1/(4*pi*background(j,5)**2)
            phi_pref=ii*pref*vt**2*ks/(wc*nu);b_pref=-pref*vt**3/(wc*nu*clight)
            qn=pref*vt**2*background(j,13)/(wc*clight)
            projection=exp(-ii*wave*r)/N
            do col=1,dim
                damping=diffusion*(wave(col)**2+ks**2)/nu
                call damped_ou_moments(kpar*vt/nu,-omega/nu,damping,moments,info,energy_moments)
                if(info/=0) error stop 'diffused cylinder moment resolvent failed'
                aux(1)=-pref+phi_pref*((a1+a2)*moments(0,0)+0.5_real64*a2*moments(0,2))
                aux(2)=b_pref*((a1+a2)*moments(0,1)+0.5_real64*a2*moments(0,3))
                aux(3)=vt*phi_pref*((a1+a2)*moments(1,0)+0.5_real64*a2*moments(1,2))
                aux(4)=vt*b_pref*((a1+a2)*moments(1,1)+0.5_real64*a2*moments(1,3))
                if(channels==9) then
                    aux(5)=phi_pref/qn*((a1+a2)*energy_moments(0)+0.5_real64*a2*energy_moments(2))
                    aux(6)=b_pref/qn*((a1+a2)*energy_moments(1)+0.5_real64*a2*energy_moments(3))
                    aux(7)=qn*energy_moments(0)
                    aux(8)=qn*vt*energy_moments(1)
                    aux(9)=energy_moments(4)
                end if
                source_phase=exp(ii*wave(col)*r)
                do k=1,channels
                    reference(:,col,k)=reference(:,col,k)+projection*source_phase*aux(k)
                end do
            end do
        end do
    end subroutine
    subroutine assemble_cylindrical_electrons(background,M,L,clight,blocks,refinement,model,bare_blocks)
        real(real64),intent(in) :: background(:,:),L,clight
        integer,intent(in) :: M,refinement
        integer,intent(in),optional :: model
        complex(real64),allocatable,intent(out),optional :: bare_blocks(:,:,:)
        complex(real64),allocatable,intent(out) :: blocks(:,:,:)
        complex(real64),parameter :: ii=(0.0_real64,1.0_real64)
        real(real64),parameter :: pi=acos(-1.0_real64)
        real(real64) :: r,ks,kpar,omega_e,lambda_d,nu,vt,wc,a1,a2,pref,kr,krp
        complex(real64) :: moments(0:3,0:3),local(4),aux(9),phase,phi_pref,b_pref
        integer :: j,row,col,N,dim,info,selected_model
        real(real64) :: qn
        selected_model=0
        if(present(model)) selected_model=model
        if(selected_model/=0.and.selected_model/=1) error stop 'electron model must be 0 or 1'
        N=size(background,1)
        dim=2*M+1
        if (size(background,2)/=13 .or. M<0 .or. N<dim .or. &
            L<=0.0_real64 .or. clight<=0.0_real64) error stop 'invalid cylindrical background'
        if (.not.all(ieee_is_finite(background))) error stop 'nonfinite cylindrical background'
        allocate(blocks(dim,dim,4))
        blocks=0.0_real64
        if(present(bare_blocks)) then
            allocate(bare_blocks(dim,dim,9))
            bare_blocks=0.0_real64
        end if
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
            if(present(bare_blocks)) then
                ! Bare channels: rho-Phi,rho-B,j-Phi,j-B,E-Phi,E-B,rho-E,j-E,E-E.
                ! E=<p*g>, p=u^2-1, g=h/F. Reinjection source is nu*p*E,
                ! so the last three channels have no extra 1/nu prefactor.
                qn=pref*vt**2*background(j,13)/(wc*clight)
                aux(1:4)=local
                aux(5)=phi_pref/qn*((a1+a2)*(moments(2,0)-moments(0,0)) &
                    +0.5_real64*a2*(moments(2,2)-moments(0,2)))
                aux(6)=b_pref/qn*((a1+a2)*(moments(2,1)-moments(0,1)) &
                    +0.5_real64*a2*(moments(2,3)-moments(0,3)))
                aux(7)=qn*(moments(0,2)-moments(0,0))
                aux(8)=qn*vt*(moments(1,2)-moments(1,0))
                aux(9)=moments(2,2)-moments(2,0)-moments(0,2)+moments(0,0)
            end if
            if(selected_model==1) then
                call restore_cylindrical_energy(moments,kpar*vt/nu,-omega_e/nu,info)
                if(info/=0) error stop 'energy-restoring characteristic resolvent singular'
                local(1)=-pref+phi_pref*((a1+a2)*moments(0,0)+0.5_real64*a2*moments(0,2))
                local(2)=b_pref*((a1+a2)*moments(0,1)+0.5_real64*a2*moments(0,3))
                local(3)=vt*phi_pref*((a1+a2)*moments(1,0)+0.5_real64*a2*moments(1,2))
                local(4)=vt*b_pref*((a1+a2)*moments(1,1)+0.5_real64*a2*moments(1,3))
            end if
            do col=1,dim
                krp=2.0_real64*pi*real(col-M-1,real64)/L
                do row=1,dim
                    kr=2.0_real64*pi*real(row-M-1,real64)/L
                    phase=exp(ii*(krp-kr)*r)/real(N,real64)
                    blocks(row,col,:)=blocks(row,col,:)+phase*local
                    if(present(bare_blocks)) bare_blocks(row,col,:)=bare_blocks(row,col,:)+phase*aux
                end do
            end do
        end do
    end subroutine
    subroutine restore_cylindrical_energy(moments,x1,x2,info)
        ! Use rank-one correction normally. Close to an undamped resonance,
        ! cancellation in its denominator destroys relative accuracy; solve
        ! the SAME Hermite operator directly, with an adaptive velocity cutoff.
        complex(real64),intent(inout) :: moments(0:3,0:3)
        real(real64),intent(in) :: x1,x2
        integer,intent(out) :: info
        complex(real64),allocatable :: dl(:),du(:),diagonal(:),rhs(:,:)
        complex(real64) :: current(0:3,0:3),previous(0:3,0:3)
        real(real64) :: coeff(4,0:3),err,scale
        complex(real64),parameter :: ii=(0.0_real64,1.0_real64)
        integer :: order,j,m,n,lapack_info
        call restore_ou_energy(moments,info)
        if(info/=2) return
        coeff=0.0_real64
        coeff(1,0)=1.0_real64
        coeff(2,1)=1.0_real64
        coeff(1,2)=1.0_real64
        coeff(3,2)=sqrt(2.0_real64)
        coeff(2,3)=3.0_real64
        coeff(4,3)=sqrt(6.0_real64)
        previous=0.0_real64
        order=256
        do while(order<=8192)
            allocate(dl(order-1),du(order-1),diagonal(order),rhs(order,4))
            diagonal=cmplx(real([(j-1,j=1,order)],real64),-x2,real64)
            diagonal(3)=-ii*x2
            dl=ii*x1*sqrt(real([(j,j=1,order-1)],real64))
            du=dl
            rhs=0.0_real64
            rhs(:4,:)=coeff
            call zgtsv(order,4,dl,diagonal,du,rhs,order,lapack_info)
            if(lapack_info/=0) then
                info=1
                return
            end if
            do n=0,3
                do m=0,3
                    current(m,n)=sum(coeff(:,m)*rhs(:4,n+1))
                end do
            end do
            scale=sqrt(sum(abs(current)**2))
            err=sqrt(sum(abs(current-previous)**2))/max(scale,tiny(1.0_real64))
            deallocate(dl,du,diagonal,rhs)
            if(order>256.and.err<1.0e-12_real64) then
                moments=current
                info=0
                return
            end if
            previous=current
            order=2*order
        end do
        info=1
    end subroutine

    subroutine close_energy_response(bare,closed,backward_error,rcond)
        ! Fourier-space closure E = V_Phi*Phi+V_B*Br+A*E. It acts on
        ! the energy moment of the distribution, not on each particle energy.
        complex(real64),intent(in) :: bare(:,:,:)
        complex(real64),allocatable,intent(out) :: closed(:,:,:)
        real(real64),intent(out) :: backward_error,rcond
        complex(real64),allocatable :: matrix(:,:),lu(:,:),rhs(:,:),initial(:,:),work(:)
        real(real64),allocatable :: rwork(:)
        integer,allocatable :: pivots(:)
        real(real64) :: anorm,scale
        integer :: d,j,info
        d=size(bare,1)
        if(size(bare,2)/=d.or.size(bare,3)/=9) error stop 'invalid energy response channels'
        allocate(matrix(d,d),lu(d,d),rhs(d,2*d),initial(d,2*d),pivots(d),work(2*d),rwork(2*d))
        matrix=-bare(:,:,9)
        do j=1,d
            matrix(j,j)=matrix(j,j)+1.0_real64
        end do
        lu=matrix
        initial(:,:d)=bare(:,:,5)
        initial(:,d+1:)=bare(:,:,6)
        rhs=initial
        call zgesv(d,2*d,lu,d,pivots,rhs,d,info)
        if(info/=0) error stop 'energy closure solve failed'
        anorm=maxval(sum(abs(matrix),dim=1))
        call zgecon('1',d,lu,d,anorm,rcond,work,rwork,info)
        if(info/=0.or.rcond<=0.0_real64) error stop 'singular energy closure condition estimate'
        scale=sqrt(sum(abs(matrix)**2))*sqrt(sum(abs(rhs)**2))+sqrt(sum(abs(initial)**2))
        backward_error=0.0_real64
        if(scale>0) backward_error=sqrt(sum(abs(matmul(matrix,rhs)-initial)**2))/scale
        allocate(closed(d,d,4))
        closed=bare(:,:,1:4)
        closed(:,:,1)=closed(:,:,1)+matmul(bare(:,:,7),rhs(:,:d))
        closed(:,:,2)=closed(:,:,2)+matmul(bare(:,:,7),rhs(:,d+1:))
        closed(:,:,3)=closed(:,:,3)+matmul(bare(:,:,8),rhs(:,:d))
        closed(:,:,4)=closed(:,:,4)+matmul(bare(:,:,8),rhs(:,d+1:))
        if(.not.all(ieee_is_finite(real(closed))).or. &
            .not.all(ieee_is_finite(aimag(closed)))) error stop 'nonfinite closed electron response'
    end subroutine
end module
