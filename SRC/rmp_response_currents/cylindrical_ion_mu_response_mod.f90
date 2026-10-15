module cylindrical_ion_mu_response_mod
    ! Resolve the collision closure at fixed mu=v_perp^2/(2 B) [cm^2 s^-2 G^-1].
    ! lambda=mu/mu_ref; e^-lambda is the reference Gauss-Laguerre measure.
    use, intrinsic :: iso_fortran_env, only: dp=>real64
    use cylindrical_characteristics_mod, only: cylindrical_ou_resolvent
    use cylindrical_electron_response_mod, only: restore_cylindrical_energy
    implicit none
    private
    public :: ion_mu_geometry,ion_mu_source,ion_mu_deposit,assemble_cylindrical_ion_mu
contains
    subroutine ion_mu_geometry(bg,wave,lambda,mu_ref,field,gyro,measure,w)
        real(dp),intent(in) :: bg(14),wave(:),lambda,mu_ref,field
        real(dp),intent(out) :: gyro(size(wave)),measure,w
        real(dp) :: a,rho
        if(min(field,mu_ref,bg(7),bg(13))<=0.or.lambda<0) error stop 'invalid mu geometry'
        a=mu_ref*field/bg(7)**2
        w=lambda*a
        measure=a*exp(-(a-1)*lambda)
        ! omega_c scales with actual B. bg(14) also retains the explicit KIM FLR scale.
        rho=sqrt(2*w)*bg(14)*bg(13)/field
        gyro=bessel_j0(rho*sqrt(bg(2)**2+wave**2))
    end subroutine

    subroutine ion_mu_source(bg,wave,lambda,mu_ref,field,vpar,phase,clight,source)
        real(dp),intent(in) :: bg(14),wave(:),lambda,mu_ref,field,vpar,clight
        complex(dp),intent(in) :: phase(size(wave))
        complex(dp),intent(out) :: source(size(wave),3)
        complex(dp),parameter :: ii=(0.0_dp,1.0_dp)
        real(dp) :: gyro(size(wave)),measure,w,u,force
        call ion_mu_geometry(bg,wave,lambda,mu_ref,field,gyro,measure,w)
        u=vpar/bg(7)
        force=bg(9)+bg(10)*(w+u*u/2)
        source(:,1)=ii*clight*bg(2)/bg(13)*force*gyro*phase
        source(:,2)=-vpar/bg(13)*force*gyro*phase
        ! Energy reinjection acts on the gyrocentre distribution, not on a field ring.
        source(:,3)=bg(6)*(u*u-1)*phase
    end subroutine

    subroutine ion_mu_deposit(bg,wave,lambda,mu_ref,field,v0,clight,integrals,result)
        real(dp),intent(in) :: bg(14),wave(:),lambda,mu_ref,field,v0,clight
        complex(dp),intent(in) :: integrals(size(wave),3)
        complex(dp),intent(out) :: result(size(wave),size(wave),9)
        real(dp),parameter :: pi=acos(-1.0_dp)
        real(dp) :: gyro(size(wave)),measure,w,qn,p0
        integer :: row,col
        call ion_mu_geometry(bg,wave,lambda,mu_ref,field,gyro,measure,w)
        qn=bg(7)**2*bg(13)/(4*pi*bg(5)**2*bg(8)*clight)
        p0=(v0/bg(7))**2-1
        do col=1,size(wave)
            do row=1,size(wave)
                result(row,col,1:2)=qn*measure*gyro(row)*integrals(col,1:2)
                result(row,col,3:4)=v0*result(row,col,1:2)
                ! E=<p*g> at fixed mu: no perpendicular measure or observation J0.
                result(row,col,5:6)=p0*integrals(col,1:2)
                result(row,col,7)=qn*measure*gyro(row)*integrals(col,3)
                result(row,col,8)=v0*result(row,col,7)
                result(row,col,9)=p0*integrals(col,3)
            end do
        end do
    end subroutine

    subroutine assemble_cylindrical_ion_mu(background,M,L,clight,lambda,mu_ref,model,exact,bare,refinement)
        real(dp),intent(in) :: background(:,:),L,clight,lambda,mu_ref
        integer,intent(in) :: M,model,refinement
        complex(dp),allocatable,intent(out) :: exact(:,:,:),bare(:,:,:)
        real(dp),parameter :: pi=acos(-1.0_dp)
        complex(dp),parameter :: ii=(0.0_dp,1.0_dp)
        real(dp) :: wave(2*M+1),gyro(2*M+1),measure,w,vt,pref,qn,force
        complex(dp) :: moments(0:3,0:3),restored(0:3,0:3),local(9),full(4),phase,cp,cb
        integer :: N,d,j,row,col,info
        N=size(background,1);d=2*M+1
        if(size(background,2)/=14.or.N<d.or.M<0.or.L<=0.or.clight<=0) &
            error stop 'invalid fixed-mu background'
        if(model/=0.and.model/=1) error stop 'fixed-mu closure supports OU models 0/1'
        allocate(exact(d,d,4),bare(d,d,9));exact=0;bare=0
        wave=2*pi*real([(j-M-1,j=1,d)],dp)/L
        do j=1,N
            vt=background(j,7)
            call cylindrical_ou_resolvent(background(j,3)*vt/background(j,6), &
                -background(j,4)/background(j,6),moments,info,refinement)
            if(info/=0) error stop 'fixed-mu bare resolvent failed'
            restored=moments
            if(model==1) then
                call restore_cylindrical_energy(restored,background(j,3)*vt/background(j,6), &
                    -background(j,4)/background(j,6),info)
                if(info/=0) error stop 'fixed-mu conserving resolvent failed'
            end if
            call ion_mu_geometry(background(j,:),wave,lambda,mu_ref,background(j,13),gyro,measure,w)
            pref=1/(4*pi*background(j,5)**2)
            qn=pref*vt**2*background(j,13)/(background(j,8)*clight)
            cp=ii*pref*vt**2*background(j,2)/(background(j,8)*background(j,6))
            cb=-pref*vt**3/(background(j,8)*background(j,6)*clight)
            force=background(j,9)+background(j,10)*w
            do col=1,d
                do row=1,d
                    call particle_channels(moments,local(1:4))
                    call particle_channels(restored,full)
                    local(5)=cp/qn*gyro(col)*(force*(moments(2,0)-moments(0,0)) &
                        +background(j,10)/2*(moments(2,2)-moments(0,2)))
                    local(6)=cb/qn*gyro(col)*(force*(moments(2,1)-moments(0,1)) &
                        +background(j,10)/2*(moments(2,3)-moments(0,3)))
                    local(7)=qn*measure*gyro(row)*(moments(0,2)-moments(0,0))
                    local(8)=qn*measure*gyro(row)*vt*(moments(1,2)-moments(1,0))
                    local(9)=moments(2,2)-moments(2,0)-moments(0,2)+moments(0,0)
                    phase=exp(ii*(wave(col)-wave(row))*background(j,1))/N
                    bare(row,col,:)=bare(row,col,:)+phase*local
                    exact(row,col,:)=exact(row,col,:)+phase*full
                end do
            end do
        end do
    contains
        subroutine particle_channels(mom,result)
            complex(dp),intent(in) :: mom(0:3,0:3)
            complex(dp),intent(out) :: result(4)
            real(dp) :: product
            product=measure*gyro(row)*gyro(col)
            result(1)=-pref*measure+cp*product*(force*mom(0,0)+background(j,10)/2*mom(0,2))
            result(2)=cb*product*(force*mom(0,1)+background(j,10)/2*mom(0,3))
            result(3)=vt*cp*product*(force*mom(1,0)+background(j,10)/2*mom(1,2))
            result(4)=vt*cb*product*(force*mom(1,1)+background(j,10)/2*mom(1,3))
        end subroutine
    end subroutine
end module
