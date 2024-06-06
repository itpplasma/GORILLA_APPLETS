module poly_without_precomp_mod1
!
    implicit none
!
    double precision, dimension(4,4)    :: amat
    double precision, dimension(4)      :: amat_in_z
!    
    !$OMP THREADPRIVATE(amat,amat_in_z)
!    
end module poly_without_precomp_mod1
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
module pusher_tetra_field_lines_mod

    implicit none
!
    private
!    
    integer                             :: iface_init
    integer, public, protected          :: ind_tetr
    integer, public, protected          :: sign_rhs
    double precision                    :: vmod0
    double precision, public, protected :: perpinv,dt_dtau_const,bmod0
    double precision                    :: t_remain
    double precision, dimension(3)      :: x_init
    double precision, dimension(4), public, protected :: z_init
    double precision                    :: k1, k3
!
    !$OMP THREADPRIVATE(ind_tetr,iface_init,perpinv,dt_dtau_const,bmod0,t_remain,x_init,  &
    !$OMP& z_init,k1,k3,vmod0,sign_rhs)
! 
    public :: pusher_tetra_field_lines,analytic_integration
!  
    contains
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
        subroutine initialize_pusher_tetra_poly(ind_tetr_inout,x,iface,vpar,t_remain_in)
!
            use tetra_physics_mod, only: tetra_physics, sign_sqg
            use supporting_functions_mod, only: bmod_func, phi_elec_func, v2_E_mod_func
!
            implicit none
!
            integer, intent(in)                        :: ind_tetr_inout,iface
            double precision, intent(in)               :: vpar
            double precision, dimension(3), intent(in) :: x
            double precision                           :: vperp2,vpar2,t_remain_in,phi_elec
!
            t_remain = t_remain_in
!    
            ind_tetr=ind_tetr_inout           !Save the index of the tetrahedron locally
!
            !Sign of the right hand side of ODE - ensures that tau is ALWAYS positive inside the algorithm
            sign_rhs = sign_sqg * int(sign(1.d0,t_remain))
!
            z_init(1:3)=x-tetra_physics(ind_tetr)%x1       !x is the entry point of the particle into the tetrahedron in (R,phi,z)-coordinates
!
            z_init(4)=vpar                         !Transform to z_init: 1st vertex of tetrahdron is new origin
!
            !Save initial orbit parameters
            x_init = x
!             vperp_init = vperp
            iface_init = iface
!
            !Tetrahedron constants
            dt_dtau_const = tetra_physics(ind_tetr)%dt_dtau_const
!
            !Multiply with sign of rhs - ensures that tau is ALWAYS positive inside the algorithm
            dt_dtau_const = dt_dtau_const*dble(sign_rhs)
!    
            !Module of B at the entry point of the particle
            bmod0 = bmod_func(z_init(1:3),ind_tetr) !tetra_physics(ind_tetr)%bmod1+sum(tetra_physics(ind_tetr)%gb*z_init(1:3))
!
            !Phi at the entry point of the particle
            phi_elec = phi_elec_func(z_init(1:3),ind_tetr)   !tetra_physics(ind_tetr)%Phi1+sum(tetra_physics(ind_tetr)%gPhi*z_init(1:3))
!
            !Auxiliary quantities
            vperp2 = -2.d0*perpinv*bmod0
            vpar2 = vpar**2
            !This is the total speed viewed in the MOVING FRAME of ExB drift (here it only acts as a coefficient for the EOM-set)
            !For physical estimates v_E is considered seperately anyway
            vmod0 = sqrt(vpar2+vperp2)
!
            k1 = vperp2+vpar2+2.d0*perpinv*tetra_physics(ind_tetr)%bmod1
            k3 = tetra_physics(ind_tetr)%Phi1-phi_elec
!
        end subroutine initialize_pusher_tetra_poly
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
        subroutine pusher_tetra_field_lines(ind_tetr_inout,iface,x,vpar,z_save,t_remain_in,t_pass, &
                        & boole_t_finished,iper_phi)
!
            use tetra_physics_mod, only: tetra_physics,particle_charge,particle_mass
            use gorilla_diag_mod,  only: diag_pusher_tetry_poly
            use pusher_tetra_func_mod, only: pusher_handover2neighbour
            use gorilla_settings_mod, only: boole_guess, boole_adaptive_time_steps
!
            implicit none
!
            double precision, intent(in)                          :: t_remain_in
!
            integer,intent(inout)                                 :: ind_tetr_inout,iface
            double precision, dimension(3), intent(inout)         :: x
            double precision, intent(inout)                       :: vpar
!
            double precision, dimension(3), intent(out)           :: z_save
            double precision, intent(out)                         :: t_pass
            logical, intent(out)                                  :: boole_t_finished
            integer,intent(out)                                   :: iper_phi
!
            logical, dimension(4)                                 :: boole_faces
            integer                                               :: i,j,k
            double precision, dimension(4)                        :: z,z_dummy
            integer                                               :: iface_new
            double precision                                      :: tau,vperp2,tau_save,tau_max, energy_init,energy_current
            logical                                               :: boole_analytical_approx,boole_face_correct
            integer                                               :: i_step_root
            double precision                                      :: t_remain_new
!
            call initialize_pusher_tetra_poly(ind_tetr_inout,x,iface,vpar,t_remain_in)
!
            !Initial computation values
            z = z_init
! 
if(diag_pusher_tetry_poly) then
    print *, 'z_init', z_init
    print *, 'iface init', iface_init
    print *, 'norm at start'
    do i = 1,4
        print *,i, 'norm', normal_distance_func(z(1:3),i)
    enddo
!
endif
!
            !Initialize iper_phi (for the case that handover2neighbour does not set this value)
            iper_phi = 0
!
            !Initialize boole_t_finished
            boole_t_finished = .false.
!
            !Set iface_new start value for polynomial approximation
            iface_new = iface_init !instead of iface_new = iface
!
            !boole_faces ... Boolean array for switching on/off root computation for individual face
            !Initialize boole_faces, such that roots for all 4 faces are computed
            boole_faces = .true.
!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!FIRST ATTEMPT WITH SECOND ORDER GUESS!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
            !Analytical calculation of orbit parameter to guess exit face and estimate of tau
            call analytic_approx(boole_faces,z,iface_new,tau,boole_analytical_approx)
!
if(diag_pusher_tetry_poly) print *, 'boole_analytical_approx',boole_analytical_approx
!
            !Initialize face error recognition procedure
            boole_face_correct = .true.   
!
            !Analytical result does not exist.
            if(.not.boole_analytical_approx) then
                boole_face_correct = .false.
if(diag_pusher_tetry_poly) print *, 'Error in predicted integrator: Analytic approximation'
            endif
!
            !Integrate trajectory analytically, if root exists
            if(boole_face_correct) then
                call analytic_integration(z,tau)
!
                call check_three_planes(z,iface_new,boole_face_correct)
                call check_face_convergence(z,iface_new,boole_face_correct)        
!
            endif !boole face correct
!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!END OF FIRST ATTEMPT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
if(diag_pusher_tetry_poly) then
    print *, 'boole_face_correct',boole_face_correct
    print *, 'After orbit pushing'
    print *, 'iface_new',iface_new
    print *, 'norm'
    do i = 1,4
        print *,i, 'norm', normal_distance_func(z(1:3),i) 
    enddo
endif
!
            !Final processing
            x=z(1:3)+tetra_physics(ind_tetr)%x1
            vpar=z(4)
!
            !Compute passing time dependent
            t_pass = tau*dt_dtau_const
!
if(diag_pusher_tetry_poly) print *, 'tau total',tau
if(diag_pusher_tetry_poly) print *, 't_pass',t_pass
!if(diag_pusher_tetry_poly) then
!    print *, 't_remain',t_remain
!    if (t_remain .lt. 0) stop
!    if (t_pass .lt. 0) stop
!endif
!
            !Particle stops inside the tetrahedron - Absolute value is used, because negative time can be allowed
            if(abs(t_pass).ge.abs(t_remain)) then
!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!SECOND ATTEMPT IF PARTICLE DOES NOT LEAVE CELL IN REMAINING TIME!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
                !Set z back to z_init
                z = z_init
!
                !Compute orbit parameter tau from t_remain
                tau = t_remain/dt_dtau_const

                
if(diag_pusher_tetry_poly) print *, 'tau until t finished',tau
!
                !Integrate trajectory analytically from start until t_remain
                call analytic_integration(z,tau)
!
                ind_tetr_inout = ind_tetr
                iface = 0
!
                boole_face_correct = .true.
                !Particle must be inside the tetrahedron
                do i = 1,4
                    if(normal_distance_func(z(1:3),i).lt.0.d0) then
                        boole_face_correct = .false.
                    endif
                enddo
!
                if(boole_face_correct) then
                    boole_t_finished = .true.
                    
                    z_save = z(1:3)
                    x=z(1:3)+tetra_physics(ind_tetr)%x1
                    vpar=z(4)
!
                    t_pass = t_remain
!                    
                else 
                    print *, 'Error: Particle is outside the tetrahedron when time is finished.'
                endif
!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!END OF SECOND ATTEMPT IF PARTICLE DOES NOT LEAVE CELL!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
            !Normal orbit that passes the whole tetrahedron
            else
!            
                !Save relative coordinates after pushing
                z_save = z(1:3)
!
                iface = iface_new
                call pusher_handover2neighbour(ind_tetr,ind_tetr_inout,iface,x,iper_phi)
!                
            endif
!
        end subroutine pusher_tetra_field_lines
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
        subroutine check_three_planes(z,iface_new,boole_face_correct)

            use gorilla_diag_mod,only: diag_pusher_tetry_poly

            implicit none
            double precision, dimension(4), intent(in)         :: z
            integer, intent(in)                                :: iface_new
            logical                                            :: boole_face_correct
            integer                                            :: j,k

            !Validation loop ('3-planes'-control)
            do j=1,3    !Just consider the planes without the "exit-plane"
            k=modulo(iface_new+j-1,4)+1
                if(normal_distance_func(z(1:3),k).lt.0.d0) then     !If distance is negative, exitpoint of the considered plane is outside the tetrahedron
                    boole_face_correct = .false.
if(diag_pusher_tetry_poly) print *, 'Error: three planes'
if(diag_pusher_tetry_poly) print *, 'face', k,'normal_distance',normal_distance_func(z(1:3),k)
                endif        
            enddo

        end subroutine check_three_planes
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
        subroutine check_face_convergence(z,iface_new,boole_face_correct)

            use gorilla_diag_mod,only: diag_pusher_tetry_poly

            implicit none
            double precision, dimension(4), intent(in)         :: z
            integer, intent(in)                                :: iface_new
            logical                                            :: boole_face_correct

        if(abs(normal_distance_func(z(1:3),iface_new)).gt.1.d-11) then
            boole_face_correct = .false.
if(diag_pusher_tetry_poly) print *, 'Error: distance'
        endif

        
        end subroutine check_face_convergence
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine analytic_approx(boole_faces,z,iface_inout,dtau,boole_approx) 
!
        use gorilla_diag_mod, only: diag_pusher_tetry_poly
!
        implicit none
!
        integer,intent(inout)                               :: iface_inout
        logical, dimension(4), intent(in)                   :: boole_faces
        double precision,intent(out)                        :: dtau
        double precision, dimension(4),intent(in)           :: z
        logical,intent(out)                                 :: boole_approx
        integer                                             :: i,iface
        double precision, dimension(4)                      :: dtau_vec
        logical, dimension(4)                               :: boole_valid_dtau
        double precision, dimension(4,2)                    :: coef_mat
        double precision                                    :: lin_a, lin_b
!
        !Calculation of root of polynomial equation
        call analytic_coeff_without_precomp(boole_faces,z,coef_mat)
!     
        iface = iface_inout   
!
        dtau_vec = huge(0.d0) !Only results with positive times smaller than 'huge' will be evaluated
!
        !loop over faces
        faces_loop: do i = 1,4
            if(.not.boole_faces(i)) cycle
!
            !Classification of coefficients
            if( (i.eq.iface).or.(coef_mat(i,1).eq.0.d0) ) then
                dtau_vec(i) = 0.d0
                cycle !No further computations are needed below
            else
                lin_a = coef_mat(i,2)
                lin_b = coef_mat(i,1)
!
                !Handling, if coefficients are exactly ZERO
                if(lin_a.eq.0.d0) then
                    dtau_vec(i) = 0.d0
                    cycle
                endif ! lin_a.eq.0.d0
!
            endif
!
            !Computation of roots
            call Linear_Solver(lin_a,lin_b,dtau_vec(i))
!
        enddo faces_loop
!
        boole_valid_dtau = (dtau_vec.lt.huge(0.d0)).and.(dtau_vec.gt.0.d0)
        if( any( boole_valid_dtau) ) then
            boole_approx = .true.
!
            iface_inout = minloc(dtau_vec,1,boole_valid_dtau)
            dtau = dtau_vec(iface_inout)
        else
            boole_approx = .false.
        endif
! 
    end subroutine analytic_approx
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
        subroutine analytic_coeff_without_precomp(boole_faces,z,coef_mat)
!
            use tetra_physics_mod, only: tetra_physics,cm_over_e
            use constants, only: clight
            use poly_without_precomp_mod1
!
            implicit none
!
            logical, dimension(4), intent(in)   :: boole_faces
            logical, dimension(4)               :: boole_faces_not
            integer                             :: n
            double precision, dimension(4)      :: z
            double precision                    :: dist1
            double precision, dimension(4,2)      :: coef_mat
!
            amat = 0.d0
            amat(1:3,4) = tetra_physics(ind_tetr)%curlA
!
            !Multiply amat with appropriate sign (which ensures that tau remains positive inside the algorithm)
            amat = amat * dble(sign_rhs)
!
            dist1= -tetra_physics(ind_tetr)%dist_ref
!
            boole_faces_not = .not.boole_faces
!
                do n = 1,4
                    if(boole_faces_not(n)) cycle
                    coef_mat(n,1)=sum(tetra_physics(ind_tetr)%anorm(:,n)*z(1:3))
                enddo
                coef_mat(1,1)=coef_mat(1,1)-dist1

                amat_in_z = matmul(amat,z)
                do n = 1,4
                    if(boole_faces_not(n)) cycle
                    coef_mat(n,2) = sum(tetra_physics(ind_tetr)%anorm(:,n)*amat_in_z(1:3))
                enddo
!
        end subroutine analytic_coeff_without_precomp
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine Linear_Solver(a,b,dtau)
!
        !Find the root of the equation
        !f(tau) = a * tau + b
!
        implicit none
!
        double precision,intent(out)                :: dtau
        double precision,intent(in)                 :: a,b
!
        if(a.eq.0.d0) then
            dtau = huge(0.d0)
        else
            dtau = -b/a    
        endif        
!
    end subroutine Linear_Solver
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
        subroutine analytic_integration(z,tau)
!
            use poly_without_precomp_mod1
!
            implicit none
!
            double precision, intent(in)                 :: tau
            double precision, dimension(4),intent(inout) :: z
!
            z = z + tau*(amat_in_z)
!
        end subroutine analytic_integration
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!    
    function normal_distance_func(z123,iface)
!
        use tetra_physics_mod, only: tetra_physics
!
        implicit none
!
        integer,intent(in)                          :: iface
        double precision, dimension(3),intent(in)   :: z123
        double precision                            :: normal_distance_func,dist1
!
        dist1= -tetra_physics(ind_tetr)%dist_ref
!
        normal_distance_func=sum(z123*tetra_physics(ind_tetr)%anorm(:,iface))
        if(iface.eq.1) normal_distance_func=normal_distance_func-dist1   !Correction for the one plane that is not lying in the first vertex
!        
    end function normal_distance_func
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
end module pusher_tetra_field_lines_mod