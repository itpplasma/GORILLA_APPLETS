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
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
module field_line_tracing_mod
    implicit none
!
    private
!
    double precision :: time_step,energy_eV, min_poloidal_flux, max_poloidal_flux, amax
    integer, dimension(:,:), allocatable :: tetra_indices_per_prism
    double precision, dimension(:), allocatable :: prism_volumes, refined_prism_volumes, elec_pot_vec, n_b
    double precision, dimension(:,:), allocatable :: verts, sqrt_g, r_integrand_constants
    integer :: seed_option, n_prisms, num_particles,&
               & ind_a, ind_b, ind_c, n_pushings, counter_phi_0_mappings, lost_outside, lost_inside
    double precision :: n_particles, density, constant_part_of_weights
    double complex, dimension(:,:), allocatable :: weights
    double precision, dimension(:), allocatable :: J_perp, poloidal_flux, temperature_vector
    logical :: boole_refined_sqrt_g, boole_boltzmann_energies
    character(1024) :: filename_starting_conditions, filename_vertex_coordinates, &
    & filename_vertex_indices
    integer :: n_fourier_modes, n_triangles
    logical :: boole_linear_density_simulation, boole_antithetic_variate, boole_linear_temperature_simulation
    logical :: boole_collisions, boole_point_source, boole_precalc_collisions
    double precision, dimension(:,:,:), allocatable :: randcol
    integer :: randcoli = int(1.0d5)
!
    !Namelist for boltzmann input
    NAMELIST /field_line_tracing_nml/ time_step,energy_eV,n_particles,boole_point_source,boole_collisions, &
    & boole_precalc_collisions,density,boole_refined_sqrt_g,boole_boltzmann_energies, boole_linear_density_simulation, &
    & boole_antithetic_variate,boole_linear_temperature_simulation,seed_option, &
    & filename_starting_conditions,filename_vertex_coordinates, filename_vertex_indices
!
    public :: calc_field_lines
!    
contains
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine field_line_tracing_inp()
!
    open(unit=71, file='field_line_tracing.inp', status='unknown')
    read(71,nml=field_line_tracing_nml)
    close(71)
!    
    print *,'GORILLA: Loaded input data from field_line_tracing.inp'
end subroutine field_line_tracing_inp
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine calc_starting_conditions(v0,start_pos_pitch_mat)
!
    use constants, only: pi, ev2erg
    use tetra_grid_mod, only: verts_rphiz, verts_sthetaphi, ntetr
    use find_tetra_mod, only: find_tetra
    use tetra_grid_settings_mod, only: grid_kind
    use tetra_physics_mod, only: coord_system
    use collis_ions, only: collis_init, stost
!
    implicit none
    double precision, intent(in)                                   :: v0
    double precision, dimension(:,:), allocatable, intent(out)     :: start_pos_pitch_mat
    double precision                                               :: rand_scalar, vpar, vperp
    double precision                                               :: amin, cmin, cmax !amax is set globally
    double precision, dimension(:), allocatable                    :: rand_vector
    double precision, dimension(:), allocatable                    :: rand_matrix2
    integer                                                        :: i
    double precision, dimension(3)                                 :: x
    integer                                                        :: ind_tetr_out,iface
!
!!!!comment out the following section to make starting conditions really random!!!
!
    integer,dimension(:), allocatable                              :: seed
    integer                                                        :: n
!
    open(unit = 85, file='seed.inp', status='old',action = 'read')
    read(85,*) n
    allocate(seed(n))
    read(85,*) seed
    close(85)
    CALL RANDOM_SEED (PUT=seed)
    deallocate(seed)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    allocate(start_pos_pitch_mat(5,num_particles))
    allocate(rand_vector(num_particles))
    allocate(rand_matrix2(num_particles))
!
    start_pos_pitch_mat = 0
!
    ind_a = 1 !(R in cylindrical coordinates, s in flux coordinates)
    ind_b = 2 !(phi in cylindrical and flux coordinates)
    ind_c = 3 !(z in cylindrical coordinates, theta in flux coordinates)
!
    if (coord_system.eq.2) then
        ind_b = 3
        ind_c = 2
    endif
!
    if (coord_system.eq.1) allocate(verts(size(verts_rphiz(:,1)),size(verts_rphiz(1,:))))
    if (coord_system.eq.2) allocate(verts(size(verts_sthetaphi(:,1)),size(verts_sthetaphi(1,:))))
    if (coord_system.eq.1) verts = verts_rphiz
    if (coord_system.eq.2) verts = verts_sthetaphi
!
    amin = minval(verts(ind_a,:))
    amax = maxval(verts(ind_a,:))
    cmin = minval(verts(ind_c,:))
    cmax = maxval(verts(ind_c,:))
!
    constant_part_of_weights = density*(amax-amin)*(cmax-cmin)*2*pi
!
    !compute starting conditions
    if (boole_point_source) then
        if (grid_kind.eq.2) then
            start_pos_pitch_mat(1,:) = 160
            start_pos_pitch_mat(2,:) = 0
            start_pos_pitch_mat(3,:) = 70
        elseif (grid_kind.eq.4) then
            start_pos_pitch_mat(1,:) = 205
            start_pos_pitch_mat(2,:) = 0
            start_pos_pitch_mat(3,:) = 0
        endif
        if (coord_system.eq.2) print*, 'error: point source is only implemented for cylindrical coordinate system'
    else
        start_pos_pitch_mat(ind_a,:) = (/(214 + i*(216-214)/n_particles, i=1,num_particles)/)!r 
        start_pos_pitch_mat(ind_b,:) = 0.d0 !phi in cylindrical and flux coordinates
        start_pos_pitch_mat(ind_c,:) = 12d0 !z in cylindrical, theta in flux coordinates
        print*, maxval(start_pos_pitch_mat(ind_a,:))
    endif

!     nsurf=5000 !1000
!     nmap=3000 !300 !00
!     nskip=1
!     icount_begplot = 10
!     rbeg=214.d0
!   !  rend=216.d0 !216.d0
!     rend=216.d0
!   !
!   !
!     do isurf=1,nsurf
!       rrr=rbeg+(rend-rbeg)*isurf/nsurf
!       z=12d0
!       phi = 0.d0
!
    
    start_pos_pitch_mat(4,:) = 1 !delete this once i have a proper subroutine for field line tracing
!
    call RANDOM_NUMBER(rand_matrix2)
    if (boole_boltzmann_energies) then !compare with equation 133 of master thesis of Jonatan Schatzlmayr (remaining parts will be added later)
        start_pos_pitch_mat(5,:) = 5*energy_eV*rand_matrix2(:) !boltzmann energy distribution
        constant_part_of_weights = constant_part_of_weights*10/sqrt(pi)*energy_eV*ev2erg
    endif
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! start antithetic variate !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (boole_antithetic_variate) then
        start_pos_pitch_mat(:,1:num_particles:2) = start_pos_pitch_mat(:,2:num_particles:2)
        start_pos_pitch_mat(4,1:num_particles:2) = -start_pos_pitch_mat(4,2:num_particles:2)
    endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! end antithetic variate !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    weights(:,1) = constant_part_of_weights
    if (boole_refined_sqrt_g.eqv..false.) weights(:,1) = constant_part_of_weights*start_pos_pitch_mat(ind_a,:)
!
    if (boole_precalc_collisions) then
        allocate(randcol(num_particles,randcoli,3))
        call RANDOM_NUMBER(randcol)
        !3.464102 = sqrt(12), this creates a random number with zero average and unit variance
        randcol(:,:,1:2:3) =  3.464102*(randcol(:,:,1:2:3)-.5)
    endif
end subroutine calc_starting_conditions
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine calc_field_lines
!
    use orbit_timestep_gorilla_mod, only: initialize_gorilla
    use constants, only: ev2erg,pi,echarge,ame,amp,clight
    use tetra_physics_mod, only: particle_mass,particle_charge,cm_over_e,mag_axis_R0, coord_system, tetra_physics
    use omp_lib, only: omp_get_thread_num, omp_get_num_threads, omp_set_num_threads
    use parmot_mod, only: rmu,ro0
    use velo_mod, only: isw_field_type
    use supporting_functions_mod, only: theta_sym_flux2theta_vmec,theta_vmec2theta_sym_flux
    use tetra_grid_settings_mod, only: n_field_periods, grid_size
    use tetra_grid_mod, only: ntetr, nvert, verts_rphiz, tetra_grid, verts_sthetaphi
    use gorilla_settings_mod, only: ispecies
    use collis_ions, only: collis_init, stost
    use find_tetra_mod, only: find_tetra
!
    implicit none
!
    double precision, dimension(:,:), allocatable :: start_pos_pitch_mat, dens_mat, temp_mat, vpar_mat, efcolf_mat, &
                                                     velrat_mat, enrat_mat, dens_mat_tetr, temp_mat_tetr
    double precision :: v0,pitchpar,vpar,vperp,t_remain,t_confined, v, maxcol
    integer :: kpart,i,j,n,l,m,k,p,ind_tetr,iface,n_lost_particles,ierr,err,iantithetic, num_background_species, inorout
    integer :: n_start, n_end, i_part, count_integration_steps
    double precision, dimension(3) :: x_rand_beg,x,randnum
    logical :: boole_initialized,boole_particle_lost
    double precision :: dtau, dphi,dtaumin, t_step
    double precision, dimension(5) :: z, zet
    double precision :: m0,z0
    double precision, dimension(:), allocatable :: efcolf,velrat,enrat,vpar_background,mass_num,charge_num,dens,temp
!
    ! open(35, file = 'outliers.dat')
    ! close(35,status='delete')
    !Load input for boltzmann computation
    call field_line_tracing_inp()
!
    !call read_in_starting_conditions(start_pos_pitch_mat, num_particles)
!
    num_particles = int(n_particles)
    n_start = 1
    n_end = num_particles
    n_fourier_modes = 5
!
    !Initialize GORILLA
    call initialize_gorilla()
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! delete this again afterwards !!!!!!!!!!!!!!!!!!!!!!!
    if (ispecies.eq.4) particle_charge = 15*echarge
    print*, 'particle charge number = ', particle_charge/echarge
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    n_prisms = ntetr/3
    allocate(weights(num_particles,1))
    allocate(J_perp(num_particles))
    allocate(poloidal_flux(num_particles))
    allocate(temperature_vector(num_particles))
    allocate(elec_pot_vec(n_prisms))
    allocate(n_b(n_prisms))
!
    poloidal_flux = 0
    temperature_vector = 0
    elec_pot_vec = 0
    n_b = 0 !boltzmann density for uniform spatial distribution
!
    call calc_square_root_g
!
print*, 'start calc_volume_integrals'
    call calc_volume_integrals
print*, 'end calc_volume_integrals'
!
    !Compute velocity module from kinetic energy dependent on particle species
    v0=sqrt(2.d0*energy_eV*ev2erg/particle_mass)
!
print*, 'calc_starting_conditions started'
    call calc_starting_conditions(v0,start_pos_pitch_mat)
print*, 'calc_starting_conditions finished'
!
    !compute maximum poloidal flux
    call unlink('poloidal_flux.dat')
    open(77, file = 'poloidal_flux.dat')
    max_poloidal_flux = 0
    min_poloidal_flux = tetra_physics(1)%Aphi1
    do i = 1, ntetr
        max_poloidal_flux = max(max_poloidal_flux,tetra_physics(i)%Aphi1 + sum(tetra_physics(i)%gAphi* &
                            & (verts([1,2,3],tetra_grid(i)%ind_knot(4))-verts([1,2,3],tetra_grid(i)%ind_knot(1)))))
        min_poloidal_flux = min(min_poloidal_flux,tetra_physics(i)%Aphi1)
        write(77,*) tetra_physics(i)%Aphi1
    enddo
    write(77,*) max_poloidal_flux
    close(77)
    !tetra_physics(i)%gAphi, 
    !minval(verts(ind_a,tetra_grid(i)%ind_knot(4)))
!
    if (boole_collisions) then
        num_background_species = 2 
        allocate(dens_mat(num_background_species-1,size(verts_rphiz(1,:))))
        allocate(temp_mat(num_background_species,size(verts_rphiz(1,:))))
        allocate(vpar_mat(num_background_species,ntetr))
        allocate(efcolf_mat(num_background_species,ntetr))
        allocate(velrat_mat(num_background_species,ntetr))
        allocate(enrat_mat(num_background_species,ntetr))
        allocate(mass_num(num_background_species-1))
        allocate(charge_num(num_background_species-1))
        allocate(dens(num_background_species))
        allocate(temp(num_background_species))
        allocate(efcolf(num_background_species))
        allocate(velrat(num_background_species))
        allocate(enrat(num_background_species))
        mass_num = 0
        charge_num = 0
        mass_num(1) = 2
        !mass_num(2) = 3
        charge_num(1) = 1
        !charge_num(2) = 2
        dens_mat = 5.0d13
        temp_mat = energy_eV
        vpar_mat = 0 !ask Sergei when this will be needed!!!
        m0 = particle_mass/amp
        z0 = particle_charge/echarge
        print*, 'm0 = ', m0, 'z0 = ', z0
!
!!!!!!!!!!!!!!!!!!!! Alternative route is taken because data is not available per vertex but per tetrahedron !!!!!!!!!!!!!!!!!!!!!!!
!
        allocate(dens_mat_tetr(num_background_species-1,ntetr))
        allocate(temp_mat_tetr(num_background_species,ntetr))
!
        open(78, file = 'background/Te_d.dat')
        read(78,'(e16.9)') (temp_mat_tetr(2,i),i=1,ntetr/grid_size(2),3)
        close(78)
!
        open(79, file = 'background/Ti_d.dat')
        read(79,'(e16.9)') (temp_mat_tetr(1,i),i=1,ntetr/grid_size(2),3)
        close(79)
!
        open(80, file = 'background/ne_d.dat')
        read(80,'(e16.9)') (dens_mat_tetr(1,i),i=1,ntetr/grid_size(2),3)
        close(80)
!
        do i = 1,grid_size(2)-1
            temp_mat_tetr(:,i*ntetr/grid_size(2)+1:(i+1)*ntetr/grid_size(2):3) = temp_mat_tetr(:,1:ntetr/grid_size(2):3)
            dens_mat_tetr(:,i*ntetr/grid_size(2)+1:(i+1)*ntetr/grid_size(2):3) = dens_mat_tetr(:,1:ntetr/grid_size(2):3)
        enddo
        do i = 1,2
            temp_mat_tetr(:,1+i:ntetr:3) = temp_mat_tetr(:,1:ntetr:3)
            dens_mat_tetr(:,1+i:ntetr:3) = dens_mat_tetr(:,1:ntetr:3)
        enddo
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
        do i = 1, ntetr
            do j = 1,num_background_species
                !following if statement because electron density will be calculated in collis init
                if (j.lt.num_background_species) dens(j) = sum(dens_mat_tetr(j,tetra_grid(i)%ind_knot([1,2,3,4])))/4
                temp(j) = sum(temp_mat_tetr(j,tetra_grid(i)%ind_knot([1,2,3,4])))/4
            enddo
            call collis_init(m0,z0,mass_num,charge_num,dens,temp,energy_eV,v0,efcolf,velrat,enrat)
            efcolf_mat(:,i) = efcolf
            velrat_mat(:,i) = velrat
            enrat_mat(:,i) = enrat
        enddo
    endif
!
        kpart = 0
        n_lost_particles = 0
        maxcol = 0
        lost_outside = 0
        lost_inside = 0
        n_pushings = 0
        counter_phi_0_mappings = 0
        iantithetic = 1
        if (boole_antithetic_variate) iantithetic = 2
        count_integration_steps = 0
        call unlink('exit_times.dat')
        call unlink('remaining_particles.dat')
        call unlink('field_line_points.dat')
        open(75, file = 'exit_times.dat')
        open(76, file = 'remaining_particles.dat')
        open(81, file = 'field_line_points.dat')

        if (boole_collisions) deallocate(efcolf,velrat,enrat)
!
        !$OMP PARALLEL DEFAULT(NONE) &
        !$OMP& SHARED(num_particles,kpart,v0,time_step,boole_collisions, &
        !$OMP& dtau,dtaumin,n_start,n_end,tetra_indices_per_prism, &
        !$OMP& start_pos_pitch_mat,boole_boltzmann_energies,count_integration_steps, &
        !$OMP& density,energy_eV,dens_mat,temp_mat,vpar_mat,tetra_grid,iantithetic,tetra_physics, &
        !$OMP& efcolf_mat,velrat_mat,enrat_mat,num_background_species,randcol,randcoli,maxcol,boole_precalc_collisions) &
        !$OMP& FIRSTPRIVATE(particle_mass, particle_charge) &
        !$OMP& PRIVATE(p,l,n,boole_particle_lost,x_rand_beg,x,pitchpar,vpar,vperp,boole_initialized,t_step,err,zet, &
        !$OMP& ind_tetr,iface,t_remain,t_confined,z,ierr, v, &
        !$OMP& m0,z0,i,efcolf,velrat,enrat,vpar_background,inorout,randnum) &
        !$OMP& REDUCTION(+:n_lost_particles, n_pushings, counter_phi_0_mappings)
        print*, 'get number of threads', omp_get_num_threads()
        if (boole_collisions) allocate(efcolf(num_background_species),velrat(num_background_species),enrat(num_background_species))
        !$OMP DO
!
        !Loop over particles
        do p = n_start,n_end/iantithetic !1,num_particles/iantithetic
            do l = 1,iantithetic
!
                n = (p-1)*iantithetic+l
                !$omp critical
                !Counter for particles
                kpart = kpart+1 !in general not equal to n becuase of parallelisation
                boole_particle_lost = .false.
if (n_end.gt.10) then
    if (modulo(kpart,int(n_end/10)).eq.0) then
        print *, kpart, ' / ', num_particles, 'particle: ', n, 'thread: ' !, omp_get_thread_num()
    endif
else 
    print *, kpart, ' / ', num_particles, 'particle: ', n, 'thread: ' !, omp_get_thread_num()
endif
                !$omp end critical
!print*, kpart
!
!if (abs(start_pos_pitch_mat(4,n)).lt.0.5) cycle !only trace particles for which vpar > vperp
                t_step = time_step
                t_confined = 0
!
                !You need x_rand_beg(1,3), pitchpar(1) (between -1 and 1), energy is already given
                x_rand_beg = start_pos_pitch_mat(1:3,n)
                pitchpar = start_pos_pitch_mat(4,n)
                ! print*, start_pos_pitch_mat(5,n)
                ! if (n.eq.972) print*, x_rand_beg,start_pos_pitch_mat(5,n),pitchpar
!
                x = x_rand_beg
                vpar = pitchpar * v0
                vperp = sqrt(v0**2-vpar**2)
                if (boole_boltzmann_energies) then
                    v = sqrt(start_pos_pitch_mat(5,n)*ev2erg*2/particle_mass)
                    vpar = pitchpar * v
                    vperp = sqrt(v**2-vpar**2)
                endif
!
                i = 0
                do while (t_confined.lt.time_step)
                    i = i+1
                    !Orbit integration
                    if (i.eq.1) then
                        boole_initialized = .false.
                    endif
                    if (boole_collisions) then
                        if (i.eq.1) call find_tetra(x,vpar,vperp,ind_tetr,iface)
                        if (.not.(ind_tetr.eq.-1)) then
                            efcolf = efcolf_mat(:,ind_tetr)
                            velrat = velrat_mat(:,ind_tetr)
                            enrat = enrat_mat(:,ind_tetr)
                            vpar_background = vpar_mat(:,ind_tetr)
                            !print*, vpar_background
!
                            vpar = vpar - vpar_background(1)
                            !since vpar_background actually has num_background_particles entries, consider giving it as an extra
                            !optional input variable to stost, before randnum (maybe also check if radnum could then be set by 
                            !randnum = variable eve if vpar_background is not set and other variables are not set by name indexing)
                            !since it came up when writing these lines: replace expressions like
                            !"verts(size(verts_rphiz(:,1)),size(verts_rphiz(1,:)))" with "3,nvert"
                            zet(1:3) = x !spatial position
                            zet(4) = sqrt(vpar**2+vperp**2)/v0 !normalized velocity module 
                            zet(5) = vpar/sqrt(vpar**2+vperp**2) !pitch parameter
!
                            if (boole_precalc_collisions) then
                                randnum = randcol(n,mod(i-1,randcoli)+1,:) 
                                call stost(efcolf,velrat,enrat,zet,t_step,1,err,(time_step-t_confined)*v0,randnum)
                            else
                                call stost(efcolf,velrat,enrat,zet,t_step,1,err,(time_step-t_confined)*v0)
                            endif
!
                            t_step = t_step/v0
                            x = zet(1:3)
                            vpar = zet(5)*zet(4)*v0+vpar_background(1)
                            vperp = sqrt(1-zet(5)**2)*zet(4)*v0
!
                            !optionally still change particle_mass, particle_charge and cm_over_e, e.g.:
                            !particle_charge = particle_charge + echarge
                            !particle_mass = particle_mass + ame - amp 
                            !cm_over_e = clight*particle_mass/particle_charge
                        endif
                    endif
!
                    call orbit_timestep_gorilla_boltzmann(x,vpar,vperp,t_step,boole_initialized,ind_tetr,iface, &
                                    & n,v,start_pos_pitch_mat,inorout,t_remain)
!
                    t_confined = t_confined + t_step - t_remain
                    !Lost particle handling
                    if(ind_tetr.eq.-1) then
!write another if clause (if hole size = minimal .and. particle lost inside .and. boole_cut_out_hole = .true. 
!(use extra variable m in orbit routine (0 normally, 1 when lost outside, -1 when lost inside))),
!if m = 1 do as now, if m = -1 select arbitrary newp position and update x, vpar and vperp)
                        write(75,*) t_confined, x, n
                        !print*, t_confined, x
                        !$omp critical
                            n_lost_particles = n_lost_particles + 1
                            boole_particle_lost = .true.
                        !$omp end critical
                        exit
                    endif
!
                    v = sqrt(vpar**2+vperp**2)
                enddo
                !$omp critical
                !print*, i
                count_integration_steps = count_integration_steps + i
                !print*, dble(i)/dble(randcoli)
                maxcol = max(dble(i)/dble(randcoli),maxcol)
                !Write results in file
                ! !$omp critical
                !     write(99,*) n, boole_particle_lost , x_rand_beg ,pitchpar,x(1),t_confined
                ! !$omp end critical
                !$omp end critical
                if (t_confined.eq.time_step) then
                    write(76,*) x, v, vpar, vperp, i, n
                endif
            enddo
        enddo !n
        !$OMP END DO
        !$OMP END PARALLEL
!
        close(75)
        close(76)
        close(81)
!
    if (boole_precalc_collisions) print*, "maxcol = ", maxcol
    print*, 'Number of lost particles',n_lost_particles
    print*, 'max_poloidal_flux is', max_poloidal_flux
    print*, 'min_poloidal_flux is', min_poloidal_flux
    print*, 'average number of pushings = ', n_pushings/n_particles
    print*, 'average number of toroidal revolutions = ', counter_phi_0_mappings/n_particles
    print*, 'average number of integration steps = ', count_integration_steps/n_particles
!
    !delete all files before writing them to avoid confusion about files that are lying
    !around from previous runs and are commented out in the current run
    call unlink('data_mingle.dat')
    call unlink(filename_vertex_coordinates)
    call unlink(filename_vertex_indices)
    call unlink('prism_volumes.dat')
    call unlink('tetra_indices_per_prism.dat')
    call unlink('sqrt_g.dat')
    call unlink('refined_prism_volumes.dat')
    call unlink('r_integrand_constants.dat')
    call unlink('elec_pot_vec.dat')
    call unlink('boltzmann_densities.dat')
    call unlink('h_phi.dat')
    call unlink('tetrahedron_neighbours.dat')
!
    open(99,file='data_mingle.dat')
    write(99,*) 1.d0-dble(n_lost_particles)/dble(num_particles) !confined fraction
    write(99,*) grid_size(1) !grid size in r/s direction
    write(99,*) grid_size(2) !grid size in phi direction
    write(99,*) grid_size(3) !grid size in z/theta direction
    write(99,*) time_step !total tracing time
    close(99)
!
    101 format(1000(e21.14,x))
    if (coord_system.eq.1) then
        ![R,phi,Z]: Write vertex coordinates to File
        open(55, file=filename_vertex_coordinates)
        do i=1, nvert
            write(55,101) verts_rphiz(1, i), verts_rphiz(2, i), verts_rphiz(3, i)
        end do
        close(55)
    elseif (coord_system.eq.2) then
        ![s,theta,phi]: Write vertex coordinates to File
        open(55, file=filename_vertex_coordinates)
        do i=1, nvert
            write(55,101) verts_sthetaphi(1, i), verts_sthetaphi(2, i), verts_sthetaphi(3, i)
        end do
        close(55)
    endif
!
    !Write vertex indices to File
    open(56, file=filename_vertex_indices)
    do i=1, ntetr
        write(56, *) tetra_grid(i)%ind_knot([1, 2, 3, 4])
    end do
    close(56)
!
    open(58, file = 'prism_volumes.dat')
    write(58,'(ES20.10E4)') prism_volumes
    close(58)
!
    open(66, file = 'refined_prism_volumes.dat')
    write(66,'(ES20.10E4)') refined_prism_volumes
    close(66)
!
    open(68, file = 'elec_pot_vec.dat')
    write(68,'(ES20.10E4)') elec_pot_vec
    close(68)
! !
    open(70, file = 'boltzmann_densities.dat')
    write(70,'(ES20.10E4)') n_b
    close(70)
!
!
PRINT*, 'particle mass = ', particle_mass
PRINT*, 'absolute value of velocity = ', v0
PRINT*, 'particle charge = ', particle_charge
PRINT*, 'temperature = ', ev2erg*energy_eV
print*, 'energy in eV = ', energy_eV
print*, 'tracing time in seconds = ', time_step
print*, 'number of particles left through the outside = ', lost_outside
print*, 'number of particles left through the inside = ', lost_inside
!
deallocate(start_pos_pitch_mat)
!
end subroutine calc_field_lines
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine orbit_timestep_gorilla_boltzmann(x,vpar,vperp,t_step,boole_initialized,ind_tetr,iface, n,v,start_pos_pitch_mat, &
                                          & inorout, t_remain_out)
!
    use pusher_tetra_rk_mod, only: pusher_tetra_rk,initialize_const_motion_rk
    use pusher_tetra_poly_mod, only: pusher_tetra_poly,initialize_const_motion_poly
    use pusher_tetra_field_lines_mod, only: pusher_tetra_field_lines
    use tetra_physics_poly_precomp_mod , only: make_precomp_poly_perpinv, initialize_boole_precomp_poly_perpinv, &
        & alloc_precomp_poly_perpinv
    use tetra_physics_mod, only: tetra_physics,particle_charge,particle_mass,cm_over_e
    use gorilla_settings_mod, only: poly_order
    use orbit_timestep_gorilla_mod, only: check_coordinate_domain
    use supporting_functions_mod, only: bmod_func, vperp_func
    use find_tetra_mod, only: find_tetra
    use constants, only: pi, ev2erg
    use tetra_grid_mod, only: tetra_grid, ntetr
!
    implicit none
!
    double precision, dimension(3), intent(inout)   :: x
    double precision, intent(inout)                 :: vpar,vperp
    double precision, intent(in)                    :: t_step
    logical, intent(inout)                          :: boole_initialized
    integer, intent(inout)                          :: ind_tetr,iface
    double precision, intent(out), optional         :: t_remain_out
    double precision, dimension(3)                  :: z_save, x_save
    double precision                                :: vperp2,t_remain,t_pass,vpar_save, v, aphi
    logical                                         :: boole_t_finished
    integer                                         :: ind_tetr_save,iper_phi,k, m, n, inorout
    integer                                         :: single_particle_counter_phi_0_mappings
    double precision                                :: perpinv,speed, r, z, phi, B, phi_elec_func
    double precision, dimension(:,:)                :: start_pos_pitch_mat
!
    !If orbit_timestep is called for the first time without grid position
    if(.not.boole_initialized) then
!
        !Check coordinate domain (optionally perform modulo operation)
        call check_coordinate_domain(x)
!
        !Find tetrahedron index and face index for position x
        call find_tetra(x,vpar,vperp,ind_tetr,iface)
!
        !If particle doesn't lie inside any tetrahedron
        if(ind_tetr.eq.-1) then
            t_remain_out = t_step
            return
        endif
!
        r = x(1) - verts(1, tetra_grid(ind_tetr)%ind_knot(1))
        z = x(3) - verts(3,tetra_grid(ind_tetr)%ind_knot(1))
        phi = x(2) - verts(2,tetra_grid(ind_tetr)%ind_knot(1))
!
        if (boole_refined_sqrt_g) then
            weights(n,1) = weights(n,1)* &
                                    &  (sqrt_g(ind_tetr,1)+r*sqrt_g(ind_tetr,2)+z*sqrt_g(ind_tetr,3))/ &
                                    &  (sqrt_g(ind_tetr,4)+r*sqrt_g(ind_tetr,5)+z*sqrt_g(ind_tetr,6))
        endif
!
        if (boole_linear_density_simulation.or.boole_linear_temperature_simulation) then
            poloidal_flux(n) = tetra_physics(ind_tetr)%Aphi1 + sum(tetra_physics(ind_tetr)%gAphi*(/r,phi,z/)) + &
                                             & cm_over_e*vpar*&
                                             & (tetra_physics(ind_tetr)%h2_1+sum(tetra_physics(ind_tetr)%gh2*(/r,phi,z/)))
        endif
        if (boole_linear_density_simulation) then
            weights(n,1) = weights(n,1)*(max_poloidal_flux*1.1-poloidal_flux(n))/(max_poloidal_flux*1.1)
        endif
!       
        if (boole_boltzmann_energies) then
            !compare with equation 133 of master thesis of Jonatan Schatzlmayr (remaining parts have been added before)
            phi_elec_func = tetra_physics(ind_tetr)%Phi1 + sum(tetra_physics(ind_tetr)%gPhi*(/r,phi,z/))
            if (.not. boole_linear_temperature_simulation) then
                weights(n,1) = weights(n,1)*sqrt(start_pos_pitch_mat(5,n)*ev2erg)/(energy_eV*ev2erg)**1.5* &
                            & exp(-(start_pos_pitch_mat(5,n)*ev2erg+particle_charge*phi_elec_func)/(energy_eV*ev2erg))
            else
                temperature_vector(n) = energy_eV*ev2erg*(max_poloidal_flux*1.1-poloidal_flux(n))/(max_poloidal_flux*1.1)
                weights(n,1) = weights(n,1)*sqrt(start_pos_pitch_mat(5,n)*ev2erg)/temperature_vector(n)**1.5* &
                & exp(-(start_pos_pitch_mat(5,n)*ev2erg+particle_charge*phi_elec_func)/temperature_vector(n))
            endif
        endif
!
        !compute J_perp for perpendicular pressure
        B = tetra_physics(ind_tetr)%bmod1+sum(tetra_physics(ind_tetr)%gb*(/r,phi,z/))
        J_perp(n) = particle_mass*v**2*(1-start_pos_pitch_mat(4,n)**2)*cm_over_e/(2*B)*(-1)
        !-1 because of negative gyrophase
!
        boole_initialized = .true.
    endif
!           
    !Exit the subroutine after initialization, if time step equals zero
    if(t_step.eq.0.d0) return
!
    inorout = 0
!
    !Squared perpendicular velocity
    vperp2 = vperp**2
!
    !Compute relative particle position
    z_save = x-tetra_physics(ind_tetr)%x1
!
    !Compute perpendicular invariant of particle
    perpinv=-0.5d0*vperp2/bmod_func(z_save,ind_tetr)
!
    !Integrate particle orbit for given time step
    t_remain = t_step
!
    !Logical for handling time integration
    boole_t_finished = .false.

    single_particle_counter_phi_0_mappings = 0
!
    n_pushings = n_pushings-1 !set n_pushings to -1 because when entering the loop it wil go back to one without pushing
    !Loop for tetrahedron pushings until t_step is reached
    do
!
        n_pushings = n_pushings + 1
        !Domain Boundary
        if(ind_tetr.eq.-1) then
            !print *, 'WARNING: Particle lost.'
            aphi = tetra_physics(ind_tetr_save)%Aphi1+sum(tetra_physics(ind_tetr_save)%gAphi*z_save)
            !$omp critical
            if (abs(aphi-max_poloidal_flux).lt.abs(aphi-min_poloidal_flux)) then
                lost_outside = lost_outside + 1
                inorout = 1
            endif
            if (abs(aphi-max_poloidal_flux).gt.abs(aphi-min_poloidal_flux)) then
                lost_inside = lost_inside + 1
                inorout = -1
            endif
            !$omp end critical
            if( present(t_remain_out)) then
                t_remain_out = t_remain
            endif
            exit
        endif
!                
        !Save the tetrahedron index for computation of vperp in the last step
        ind_tetr_save = ind_tetr
!
        !Save vpar for the computation of the parallel adiabatic invariant
        vpar_save = vpar
!  
        !t_remain (in) ... remaining time until t_step is finished
        !t_pass (out) ... time to pass the tetrahdron

        !Save x for the computation of the current
        x_save = x
!
        !Calculate trajectory
        call pusher_tetra_field_lines(ind_tetr,iface,x,vpar,z_save,t_remain,t_pass,boole_t_finished,iper_phi)
!
        t_remain = t_remain - t_pass
        if (iper_phi.ne.0) then
            !$omp critical
            counter_phi_0_mappings = counter_phi_0_mappings + 1!iper_phi
            if (single_particle_counter_phi_0_mappings.gt.10) write(81,*) x
            !$omp end critical
            single_particle_counter_phi_0_mappings = single_particle_counter_phi_0_mappings + 1
            if (single_particle_counter_phi_0_mappings.eq.3000) then
                boole_t_finished = .true.
            endif
        endif
!
        !Orbit stops within cell, because "flight"-time t_step has finished
        if(boole_t_finished) then
            if( present(t_remain_out)) then
                t_remain_out = 0
            endif
            exit
        endif
!
    enddo !Loop for tetrahedron pushings
!
    !Compute vperp from position
    vperp = vperp_func(z_save,perpinv,ind_tetr_save)
!            
!             !NOT FULLY IMPLEMENTED YET: Deallocate precomputed quantities dependent on perpinv
!             call alloc_precomp_poly_perpinv(2,ntetr)
!print*, 'number of pushings is', n_pushings
!         
end subroutine orbit_timestep_gorilla_boltzmann
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine calc_volume_integrals
!
    use constants, only: pi, ev2erg
    use tetra_grid_mod, only: ntetr, verts_rphiz, tetra_grid
    use tetra_grid_settings_mod, only: grid_size, n_field_periods
    use tetra_physics_mod, only: particle_mass,particle_charge, tetra_physics
!
    implicit none
!
    integer                                     :: i,k
    integer, dimension(2)                       :: triangle_indices
    double precision, dimension(3)              :: r_values, z_values, r_values_intermediate, z_values_intermediate, gradient
    double precision, dimension(2)              :: r, z, limits
    double precision                            :: z_star, alpha, beta, gamma, delta, epsilon, capital_gamma, capital_delta, &
                                                   a, a_dash, b, b_dash, c, c_dash, rmin, delta_r, delta_z, phi_0, eta
!
    allocate(tetra_indices_per_prism(n_prisms,3))    
    allocate(prism_volumes(n_prisms))
    allocate(refined_prism_volumes(n_prisms))
    allocate(r_integrand_constants(n_prisms,22)) !collects all constants used for integration in order to print them on a file
!
    refined_prism_volumes = 0
    r_integrand_constants = 0
!
    do i = 1,3
        tetra_indices_per_prism(:,i) = (/(i+3*k,k = 0,n_prisms-1)/)
    enddo
!
    !$OMP PARALLEL DEFAULT(NONE) &
    !$OMP& SHARED(n_prisms,verts_rphiz,tetra_grid,grid_size,tetra_indices_per_prism,prism_volumes, particle_charge,energy_ev, &
    !$OMP& refined_prism_volumes,sqrt_g,r_integrand_constants, elec_pot_vec,n_b,tetra_physics,boole_boltzmann_energies,density, &
    !$OMP& boole_refined_sqrt_g,n_field_periods) &
    !$OMP& PRIVATE(r_values,z_values,rmin,triangle_indices,r_values_intermediate,z_values_intermediate, delta_r, delta_z, eta, &
    !$OMP& r,z,gradient,z_star,alpha,beta,gamma,delta,epsilon,capital_gamma,capital_delta, limits,a,a_dash,b,b_dash,c,c_dash,phi_0)
    !$OMP DO
    do i = 1,n_prisms
        !calculate prism volumes using the basic approach (compare with chapter 4.2 of master thesis from Jonatan Schatzlmayr)
        r_values = verts_rphiz(1,tetra_grid(tetra_indices_per_prism(i,1))%ind_knot([1,2,3]))
        z_values = verts_rphiz(3,tetra_grid(tetra_indices_per_prism(i,1))%ind_knot([1,2,3]))
        rmin = minval(r_values)
        r_values = r_values - rmin
        z_values = z_values - z_values(minloc(r_values,1))
!
        triangle_indices(2) = maxloc(r_values,1)
        r_values_intermediate = r_values
        r_values_intermediate(triangle_indices(2)) = minval(r_values)
        if (sum(r_values_intermediate).eq.0) then
            z_values_intermediate = z_values
            z_values_intermediate(triangle_indices(2)) = minval(z_values)
            triangle_indices(1) = maxloc(abs(z_values_intermediate),1)
        else
            triangle_indices(1) = maxloc(r_values_intermediate,1)
        endif
!
        r = (/r_values(triangle_indices(1)),r_values(triangle_indices(2))/)
        z = (/z_values(triangle_indices(1)),z_values(triangle_indices(2))/)
!
        if (r(2).eq.0) then
            gradient(1) = 0
        else
            gradient(1) = z(2)/r(2)
        endif
        if (r(1).eq.0) then
            gradient(2) = 0
        else
            gradient(2) = z(1)/r(1)
        endif
        if ((r(2)-r(1)).eq.0) then
            gradient(3) = 0
        else
            gradient(3) = (z(2)-z(1))/(r(2)-r(1))
        endif
!
        r_integrand_constants(i,20:22) = (/minloc(r_values),triangle_indices(1),triangle_indices(2)/)
!
        z_star = z(1) - r(1)*gradient(3)
        alpha = abs(gradient(2)-gradient(1))
        beta = gradient(3) - gradient(1)
!
        r_integrand_constants(i,1:6) = (/r(1),r(2),alpha,rmin, z_star,beta/)
!
        prism_volumes(i) =  2*pi/(grid_size(2)*n_field_periods)*(alpha/3*r(1)**3+alpha*rmin/2*r(1)**2+ &
                            abs(z_star*rmin*(r(2)-r(1))+(z_star+beta*rmin)/2*(r(2)**2-r(1)**2)+beta/3*(r(2)**3-r(1)**3)))
!
        !calculate other volme integrals (compare with appendix B of master thesis of Jonatan Schatzlmayr)
        delta_r = verts_rphiz(1,tetra_grid(tetra_indices_per_prism(i,1))%ind_knot(1)) - &
                & verts_rphiz(1,tetra_grid(tetra_indices_per_prism(i,1))%ind_knot(int(r_integrand_constants(i,20))))
        delta_z = verts_rphiz(3,tetra_grid(tetra_indices_per_prism(i,1))%ind_knot(1)) - &
                & verts_rphiz(3,tetra_grid(tetra_indices_per_prism(i,1))%ind_knot(int(r_integrand_constants(i,20))))
!
        if (boole_refined_sqrt_g) then
            !Calculate prism volumes using the refined approach
            !(compare with appendix B (introductory pages + B1) of master thesis from Jonatan Schatzlmayr)
            a = sqrt_g(3*i-2,1) - sqrt_g(3*i-2,2)*delta_r - sqrt_g(3*i-2,3)*delta_z
            a_dash = sqrt_g(3*i-2,4) - sqrt_g(3*i-2,5)*delta_r - sqrt_g(3*i-2,6)*delta_z
            b = sqrt_g(3*i-2,2)
            b_dash = sqrt_g(3*i-2,5)
            c = sqrt_g(3*i-2,3)
            c_dash = sqrt_g(3*i-2,6)
    !
            !calculate the contribution from 0 to r(1) to the prism volumes
            limits = (/dble(0),r(1)/)
    !
            alpha = c/c_dash*(gradient(2)-gradient(1))
            beta = b_dash+c_dash*gradient(2)
            gamma = b_dash+c_dash*gradient(1)
            delta = (c_dash*a-c*a_dash)/c_dash**2
            epsilon = (c_dash*b-c*b_dash)/c_dash**2
    !
            r_integrand_constants(i,7:12) = (/alpha,beta,gamma,delta,epsilon,a_dash/)
    !
            refined_prism_volumes(i) = r(1)*epsilon*a_dash/(2*gamma*beta)*(gamma-beta) + &
                    & r(1)**2*alpha/2 + &
                    & log((a_dash+gamma*r(1))/a_dash)*(delta*a_dash/(gamma*beta)*(gamma-beta)+epsilon*a_dash**2/(2*gamma**2)) + &
                    & log((a_dash+beta*r(1))/(a_dash+gamma*r(1)))*(delta/beta*(a_dash+beta*r(1))+epsilon/2*r(1)**2) - &
                    & log((a_dash+beta*r(1))/a_dash)*epsilon*a_dash**2/(2*beta**2)
    !
            !calculate the contribution from r(1) to r(2) to the prism volumes
            limits = (/r(1),r(2)/)
    !
            alpha = c/c_dash*(gradient(3)-gradient(1))
            beta = b_dash+c_dash*gradient(3)
            gamma = b_dash+c_dash*gradient(1)
            delta = (c_dash*a-c*a_dash)/c_dash**2
            epsilon = (c_dash*b-c*b_dash)/c_dash**2
            capital_gamma = c/c_dash*z_star
            capital_delta = a_dash + c_dash*z_star
    !
            r_integrand_constants(i,13:19) = (/alpha,beta,gamma,delta,epsilon,capital_gamma,capital_delta/)
    !
            refined_prism_volumes(i) = refined_prism_volumes(i) + &
                    & (r(2)-r(1))*(capital_gamma+epsilon/(2*gamma*beta)*(gamma*capital_delta-beta*a_dash)) + &
                    & (r(2)**2-r(1)**2)*alpha/2 + &
                    & log((a_dash+gamma*r(2))/(a_dash+gamma*r(1)))*(delta/(gamma*beta)*(gamma*capital_delta-a_dash*beta) + &
                                                                                            & epsilon*a_dash**2/(2*gamma**2)) + &
                    & log((capital_delta+beta*r(2))/(a_dash+gamma*r(2))*(a_dash+gamma*r(1))/(capital_delta+beta*r(1))) * &
                                                                    & (delta/beta*(capital_delta+beta*r(2))+epsilon/2*r(2)**2) + &
                    & log((capital_delta+beta*r(1))/(a_dash+gamma*r(1)))*(delta*(r(2)-r(1))+epsilon/2*(r(2)**2-r(1)**2)) - &
                    & log((capital_delta+beta*r(2))/(capital_delta+beta*r(1)))*epsilon*capital_delta**2/(2*beta**2)
        endif
!
        if (boole_boltzmann_energies) then
!
            !calculate electric potential using the refined approach
            !(compare with appendix B (introductory pages + B2) of master thesis from Jonatan Schatzlmayr)
!
            a = tetra_physics(3*i-2)%gPhi(1)
            b = tetra_physics(3*i-2)%gPhi(3)
            phi_0 = tetra_physics(3*i-2)%Phi1 - a*delta_r - b*delta_z
!
            !calculate the contribution from 0 to r(1) to the electric potential
            alpha = phi_0*(gradient(2)-gradient(1))*rmin
            beta = (a*rmin+phi_0)*(gradient(2)-gradient(1))+b/2*(gradient(2)**2-gradient(1)**2)*rmin
            gamma = a*(gradient(2)-gradient(1))+b/2*(gradient(2)**2-gradient(1)**2)

            elec_pot_vec(i) = elec_pot_vec(i) + alpha/2*r(1)**2+beta/3*r(1)**3+gamma/4*r(1)**4
!
            !calculate the contribution from r(1) to r(2) to the electric potential
            alpha  = phi_0*z_star*rmin+b/2*z_star**2*rmin
            beta = a*z_star*rmin+phi_0*(gradient(3)-gradient(1))*rmin+phi_0*z_star+b/2*z_star**2+b*z_star*gradient(3)*rmin
            gamma = phi_0*(gradient(3)-gradient(1))+a*z_star+a*(gradient(3)-gradient(1))*rmin + &
                    & b/2*(gradient(3)**2-gradient(1)**2)*rmin+b*z_star*gradient(3)
            delta = a*(gradient(3)-gradient(1))+b/2*(gradient(3)**2-gradient(1)**2)

            elec_pot_vec(i) = elec_pot_vec(i) + alpha*(r(2)-r(1))+beta/2*(r(2)**2-r(1)**2)+&
                            & gamma/3*(r(2)**3-r(1)**3)+delta/4*(r(2)**4-r(1)**4)
!
            !calculate Boltzmann densiy using the refined approach
            !(compare with appendix B (introductory pages + B3) of master thesis from Jonatan Schatzlmayr)
!
            !calculate contribution from 0 to r(1) to the boltzmann density    
            alpha = density*exp(-particle_charge*phi_0/(energy_ev*ev2erg))*rmin*(gradient(2)-gradient(1))
!
            beta = density*exp(-particle_charge*phi_0/(energy_ev*ev2erg))*((gradient(2)-gradient(1))* &
                 & (1-rmin*particle_charge*a/(energy_ev*ev2erg))- &
                 & rmin*particle_charge*b/(2*energy_ev*ev2erg)*(gradient(2)**2-gradient(1)**2))
!
            gamma = density*exp(-particle_charge*phi_0/(energy_ev*ev2erg))* &
                  & (-particle_charge*a/(energy_ev*ev2erg)*(gradient(2)-gradient(1))- &
                  & particle_charge*b/(2*energy_ev*ev2erg)*(gradient(2)**2-gradient(1)**2))
!
            n_b(i) = n_b(i) + alpha/2*r(1)**2+beta/3*r(1)**3+gamma/4*r(1)**4
!
            !calculate contribution from r(1) to r(2) to the boltzmann density
            alpha = density*exp(-particle_charge*phi_0/(energy_ev*ev2erg))*rmin*&
                  & (z_star-particle_charge*b/(2*energy_ev*ev2erg)*z_star**2)
!
            beta = density*exp(-particle_charge*phi_0/(energy_ev*ev2erg))*(z_star- & 
                 & rmin*particle_charge*a/(energy_ev*ev2erg)*z_star+rmin*(gradient(3)-gradient(1))- &
                 & particle_charge*b/(2*energy_ev*ev2erg)*z_star**2-rmin*particle_charge*b/(energy_ev*ev2erg)*z_star*gradient(3))
!
            gamma = density*exp(-particle_charge*phi_0/(energy_ev*ev2erg))*((-rmin*particle_charge*a/(energy_eV*ev2erg)+1)* &
                  & (gradient(3)-gradient(1))-particle_charge*a/(energy_ev*ev2erg)*z_star- &
                  & rmin*particle_charge*b/(2*energy_ev*ev2erg)*(gradient(3)**2-gradient(1)**2)- &
                  & particle_charge*b/(energy_eV*ev2erg)*z_star*gradient(3))
!
            delta = density*exp(-particle_charge*phi_0/(energy_ev*ev2erg))*(-particle_charge*a/(energy_eV*ev2erg)* &
                  & (gradient(3)-gradient(1))-particle_charge*b/(2*energy_eV*ev2erg)*(gradient(3)**2-gradient(1)**2))
!
            n_b(i) = n_b(i) + alpha*(r(2)-r(1)) + beta/2*(r(2)**2-r(1)**2) + gamma/3*(r(2)**3-r(1)**3) + delta/4*(r(2)**4-r(1)**4)
!
        endif
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
!
    refined_prism_volumes = abs(refined_prism_volumes)*2*pi/(grid_size(2)*n_field_periods)
    elec_pot_vec = abs(elec_pot_vec)*2*pi/(grid_size(2)*n_field_periods*prism_volumes)
    n_b = abs(n_b)*2*pi/(grid_size(2)**n_field_periods*prism_volumes)
!
end subroutine calc_volume_integrals
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine calc_square_root_g
!
    use tetra_physics_mod, only: tetra_physics, hamiltonian_time
    use tetra_grid_mod, only: ntetr, tetra_grid, verts_rphiz
!
    implicit none
!
    integer            :: ind_tetr
!
    allocate(sqrt_g(ntetr,7))
    !compare first 6 entries with chapter 4.5 of master thesis from Jonatan Schatzlmayr, entry 7 is the radius (old metric determinant)
!
    do ind_tetr = 1, ntetr
        sqrt_g(ind_tetr,1) = hamiltonian_time(ind_tetr)%h1_in_curlA
!
        sqrt_g(ind_tetr,2) = tetra_physics(ind_tetr)%gh1(1)*tetra_physics(ind_tetr)%curlA(1) + &
                           & tetra_physics(ind_tetr)%gh2(1)*tetra_physics(ind_tetr)%curlA(2) + &
                           & tetra_physics(ind_tetr)%gh3(1)*tetra_physics(ind_tetr)%curlA(3)
!
        sqrt_g(ind_tetr,3) = tetra_physics(ind_tetr)%gh1(3)*tetra_physics(ind_tetr)%curlA(1) + &
                           & tetra_physics(ind_tetr)%gh2(3)*tetra_physics(ind_tetr)%curlA(2) + &
                           & tetra_physics(ind_tetr)%gh3(3)*tetra_physics(ind_tetr)%curlA(3)
!
        sqrt_g(ind_tetr,4) = tetra_physics(ind_tetr)%bmod1
!
        sqrt_g(ind_tetr,5) = tetra_physics(ind_tetr)%gB(1)
!
        sqrt_g(ind_tetr,6) = tetra_physics(ind_tetr)%gB(3)
!
        sqrt_g(ind_tetr,7) = verts_rphiz(1,tetra_grid(ind_tetr)%ind_knot(1))
    enddo
!
end subroutine calc_square_root_g
!
end module field_line_tracing_mod