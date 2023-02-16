!
module reversibility_test_load_mod
!
    implicit none
!  
    private
!
    !Setting quantities from input file reversibility_test.inp
!
    double precision, public, protected    :: t_total
!
    double precision, public, protected    :: s_0
    double precision, public, protected    :: ds_0
    double precision, public, protected    :: theta_0
    double precision, public, protected    :: dtheta_0
    double precision, public, protected    :: R_0
    double precision, public, protected    :: dR_0
    double precision, public, protected    :: Z_0
    double precision, public, protected    :: dZ_0
    double precision, public, protected    :: phi_0
    double precision, public, protected    :: contour_fraction
    double precision, public, protected    :: pitchpar_0
    double precision, public, protected    :: dpitchpar
!
    double precision, public, protected    :: energy_eV_0
    double precision, public, protected    :: relative_bandwith
!
    integer, public, protected             :: n_steps
    integer, public, protected             :: n_orbits
    integer, public, protected             :: n_snapshots
!
    logical, public, protected             :: boole_apply_noise
    double precision, public, protected    :: noise_amplitude
    integer, public, protected             :: seed_option
!
    character(50), public, protected       :: filename_reversibility_test
    character(50), public, protected       :: filename_reversibility_test_back
!
    logical, public, protected             :: boole_diag_reversibility_test
!
    !Namelist for reversibility test input
    NAMELIST /REVERSIBILITY_TEST_NML/ t_total, s_0, ds_0, theta_0, dtheta_0, R_0, dR_0, Z_0, dZ_0, phi_0,&
                                        & contour_fraction, pitchpar_0, dpitchpar, &
                                        & energy_eV_0, relative_bandwith, n_steps, n_orbits, n_snapshots, &
                                        & boole_apply_noise, noise_amplitude, seed_option, &
                                        & filename_reversibility_test, filename_reversibility_test_back, &
                                        & boole_diag_reversibility_test
!
    public :: load_reversibility_test_inp
!
    contains

    subroutine load_reversibility_test_inp()
!
            open(unit=10, file='reversibility_test.inp', status='unknown')
            read(10,nml=reversibility_test_nml)
            close(10)
!
    end subroutine load_reversibility_test_inp

end module reversibility_test_load_mod
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
module reversibility_test_mod
!
    use reversibility_test_load_mod, only: load_reversibility_test_inp, &
                                        & boole_diag_reversibility_test
!
    implicit none
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    !starting values of additional quantites for backward integration
    type starting_values_type
        sequence
        double precision :: t_hamiltonian     !real time of tetrahedron passing
        double precision :: gyrophase         !gyrophase of particle
    end type starting_values_type
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    contains

    subroutine make_reversibility_test()
!
        use constants, only: pi
        use orbit_timestep_gorilla_mod, only: initialize_gorilla
        use reversibility_test_load_mod, only: t_total, contour_fraction, pitchpar_0, dpitchpar, &
                                            & n_steps, n_orbits, n_snapshots, energy_eV_0, relative_bandwith, boole_apply_noise, &
                                            & filename_reversibility_test, filename_reversibility_test_back
        use, intrinsic :: ieee_arithmetic, only: IEEE_VALUE, IEEE_QUIET_NAN
!
        implicit none
!
        double precision, dimension(3)                  :: x_0
        double precision, dimension(:), allocatable     :: alpha_0_vec
        double precision                                :: alpha_0, perpinv, perpinv2, energy_eV, lambda, vpar_0
        integer                                         :: t_total_sign, ind_tetr, iface
!
        double precision,dimension(:),allocatable       :: x1_vec, x2_vec, x3_vec, pitchpar_vec, &
                                                        & kin_energy_eV_vec, &
                                                        & hamiltonian_time_vec, gyrophase_vec, &
                                                        & rand_noise_vec
        double precision, dimension(:,:), allocatable   :: x1_mat, x2_mat, x3_mat, pitchpar_mat, & 
                                                        & kin_energy_eV_mat, &
                                                        & hamiltonian_time_mat, gyrophase_mat
        double precision, dimension(:,:), allocatable   :: x1_back_mat, x2_back_mat, x3_back_mat, pitchpar_back_mat, & 
                                                        & kin_energy_eV_back_mat, &
                                                        & hamiltonian_time_back_mat, gyrophase_back_mat
        type(starting_values_type)                      :: starting_values
!
        integer                                         :: i,j,k,l, counter_particles, delta_snapshot
        integer                                         :: file_id_reversibility_test,file_id_reversibility_test_back
        double precision, dimension(8)                  :: snapshot_end
!
!------------------------------------------------------------------------------------------------------------!
! Loading settings
!
        call load_reversibility_test_inp()    
!
!------------------------------------------------------------------------------------------------------------!
! Allocations and preparations
!
        !Allocate matrices to save forward data
        allocate(x1_mat(n_steps+1,n_orbits+1), x2_mat(n_steps+1,n_orbits+1), x3_mat(n_steps+1,n_orbits+1), & 
            & pitchpar_mat(n_steps+1,n_orbits+1), kin_energy_eV_mat(n_steps+1,n_orbits+1), &
            & hamiltonian_time_mat(n_steps+1,n_orbits+1), gyrophase_mat(n_steps+1,n_orbits+1))

        !Allocate matrices to save backward data
        allocate(x1_back_mat(n_steps+1,n_orbits+1), x2_back_mat(n_steps+1,n_orbits+1), x3_back_mat(n_steps+1,n_orbits+1), & 
        & pitchpar_back_mat(n_steps+1,n_orbits+1), kin_energy_eV_back_mat(n_steps+1,n_orbits+1), &
        & hamiltonian_time_back_mat(n_steps+1,n_orbits+1), gyrophase_back_mat(n_steps+1,n_orbits+1))
!
        !Allocate vectors to save individual orbit data
        allocate(x1_vec(n_steps+1), x2_vec(n_steps+1), x3_vec(n_steps+1), &
                & pitchpar_vec(n_steps+1),kin_energy_eV_vec(n_steps+1), &
                & hamiltonian_time_vec(n_steps+1), gyrophase_vec(n_steps+1))
!
        !Create vector with equidistant theta values
        allocate(alpha_0_vec(n_orbits+1))
        alpha_0_vec(:) = contour_fraction * [(2.d0*pi/dble(n_orbits)*dble(i), i = 0,n_orbits)]
!
        !Initialize GORILLA and data save matrices
        call initialize_gorilla()
        t_total_sign = int(sign(1.d0,t_total))
        if (boole_apply_noise) then
            allocate(rand_noise_vec(n_orbits+1))
            rand_noise_vec = get_rand_noise_vec()
        endif
!
        x1_mat = IEEE_VALUE(x1_mat, IEEE_QUIET_NAN)
        x2_mat = IEEE_VALUE(x2_mat, IEEE_QUIET_NAN)
        x3_mat = IEEE_VALUE(x3_mat, IEEE_QUIET_NAN)
        pitchpar_mat = IEEE_VALUE(pitchpar_mat, IEEE_QUIET_NAN)
        kin_energy_eV_mat = IEEE_VALUE(kin_energy_eV_mat, IEEE_QUIET_NAN)
        hamiltonian_time_mat = IEEE_VALUE(hamiltonian_time_mat, IEEE_QUIET_NAN)
        gyrophase_mat = IEEE_VALUE(gyrophase_mat, IEEE_QUIET_NAN)
!
        x1_back_mat = IEEE_VALUE(x1_back_mat, IEEE_QUIET_NAN)
        x2_back_mat = IEEE_VALUE(x2_back_mat, IEEE_QUIET_NAN)
        x3_back_mat = IEEE_VALUE(x3_back_mat, IEEE_QUIET_NAN)
        pitchpar_back_mat = IEEE_VALUE(pitchpar_back_mat, IEEE_QUIET_NAN)
        kin_energy_eV_back_mat = IEEE_VALUE(kin_energy_eV_back_mat, IEEE_QUIET_NAN)
        hamiltonian_time_back_mat = IEEE_VALUE(hamiltonian_time_back_mat, IEEE_QUIET_NAN)
        gyrophase_back_mat = IEEE_VALUE(gyrophase_back_mat, IEEE_QUIET_NAN)
!
        snapshot_end = IEEE_VALUE(snapshot_end, IEEE_QUIET_NAN)
!
!------------------------------------------------------------------------------------------------------------!
! Compute orbits
!
        counter_particles = 0
!
!
        !$OMP PARALLEL DEFAULT(NONE) &
        !$OMP& SHARED(alpha_0_vec,t_total,n_steps,n_orbits, &
        !$OMP& pitchpar_0,dpitchpar,energy_eV_0,relative_bandwith,boole_apply_noise,rand_noise_vec, &
        !$OMP& x1_mat,x2_mat,x3_mat,pitchpar_mat,kin_energy_eV_mat, &
        !$OMP& hamiltonian_time_mat,gyrophase_mat, &
        !$OMP& x1_back_mat,x2_back_mat,x3_back_mat,pitchpar_back_mat,kin_energy_eV_back_mat, &
        !$OMP& hamiltonian_time_back_mat,gyrophase_back_mat, &
        !$OMP& counter_particles,t_total_sign) &
        !$OMP& PRIVATE(k,x_0,alpha_0,lambda,perpinv,perpinv2, &
        !$OMP& vpar_0,ind_tetr,iface, &
        !$OMP& x1_vec,x2_vec,x3_vec,pitchpar_vec,kin_energy_eV_vec,energy_eV, &
        !$OMP& hamiltonian_time_vec,gyrophase_vec, starting_values)
        !$OMP DO
!
        !Loop over orbits and store data
        do k = 1,n_orbits+1
!
            !orbit starting positions
            alpha_0 = alpha_0_vec(k)
            x_0 = starting_position(alpha_0)
!
            !Define kinetic energy for individual orbit
            energy_eV = energy_eV_0 * (1.d0 + relative_bandwith * sin(alpha_0) )
!
            !Define starting pitch parameter for individual orbit
            lambda = pitchpar_0 + dpitchpar * cos(alpha_0)
!
            !Setup to get starting tetrahedra, velocities, ...
            !energy_eV, x_0 -> x_0,vpar_0,ind_tetr,iface, ...
            call gorilla_integration_setup(lambda,t_total_sign, &
                                        & x_0,vpar_0,energy_eV,ind_tetr,iface,perpinv,perpinv2)
            if(ind_tetr.eq.-1) cycle
!
            !Compute orbits
            call gorilla_integration(perpinv,perpinv2,t_total,n_steps, &
                                    & x_0,vpar_0,energy_eV,ind_tetr,iface, &
                                    & x1_vec,x2_vec,x3_vec,pitchpar_vec,kin_energy_eV_vec,hamiltonian_time_vec,gyrophase_vec)
            if(ind_tetr.eq.-1) cycle
!
            !Write values in matrix
            x1_mat(:,k) = x1_vec
            x2_mat(:,k) = x2_vec
            x3_mat(:,k) = x3_vec
            pitchpar_mat(:,k) = pitchpar_vec
            kin_energy_eV_mat(:,k) = kin_energy_eV_vec
            hamiltonian_time_mat(:,k) = hamiltonian_time_vec
            gyrophase_mat(:,k) = gyrophase_vec

            !Save starting values of additional quantities for backward integration
            starting_values%t_hamiltonian = 0.0d0
            starting_values%gyrophase = gyrophase_vec(n_steps+1)
            !Compute orbits backwards (-t_total) from last z = (x_0,vpar_0)
            !Apply optional noise on coordinate set
            if (boole_apply_noise) call apply_noise(rand_noise_vec(k),-t_total_sign,x_0,vpar_0, & 
                                                    & energy_eV,ind_tetr,iface,perpinv,perpinv2)

            !As have already correct coordinate set, do not need extra setup like in forward case
            call gorilla_integration(perpinv,perpinv2,-t_total,n_steps, &
                                    & x_0,vpar_0,energy_eV,ind_tetr,iface, &
                                    & x1_vec,x2_vec,x3_vec,pitchpar_vec,kin_energy_eV_vec,hamiltonian_time_vec,gyrophase_vec, &
                                    & starting_values=starting_values)
            if(ind_tetr.eq.-1) cycle
!
            !Write values in backwards matrix
            x1_back_mat(:,k) = x1_vec
            x2_back_mat(:,k) = x2_vec
            x3_back_mat(:,k) = x3_vec
            pitchpar_back_mat(:,k) = pitchpar_vec
            kin_energy_eV_back_mat(:,k) = kin_energy_eV_vec
            hamiltonian_time_back_mat(:,k) = hamiltonian_time_vec
            gyrophase_back_mat(:,k) = gyrophase_vec
!
            !$omp critical
                counter_particles = counter_particles +1
                print *, 'Counter orbits', counter_particles, '/',n_orbits+1
            !$omp end critical

        enddo
        !$OMP END DO
        !$OMP END PARALLEL
!
!------------------------------------------------------------------------------------------------------------!
!

        delta_snapshot = n_steps/n_snapshots
        !File ID for orbit position
        file_id_reversibility_test = 100
        open(file_id_reversibility_test,file=filename_reversibility_test)
        !Write orbit positions
        do l = 1, n_snapshots
            j = 1 + (l-1)*delta_snapshot
            do k = 1,n_orbits+1
                write(file_id_reversibility_test,*) j, hamiltonian_time_mat(j,k), &
                    & x1_mat(j,k), x2_mat(j,k), x3_mat(j,k), &
                    & pitchpar_mat(j,k), gyrophase_mat(j,k)/(2.0d0*pi), kin_energy_eV_mat(j,k)
            enddo
            write(file_id_reversibility_test,*) snapshot_end
        enddo
        !Final orbit position
        j = n_steps + 1
        do k = 1,n_orbits+1
            write(file_id_reversibility_test,*) j, hamiltonian_time_mat(j,k), &
                & x1_mat(j,k), x2_mat(j,k), x3_mat(j,k), &
                & pitchpar_mat(j,k), gyrophase_mat(j,k)/(2.0d0*pi), kin_energy_eV_mat(j,k)
        enddo
        write(file_id_reversibility_test,*) snapshot_end
        close(file_id_reversibility_test)
!
        !File ID for orbit position backwards
        file_id_reversibility_test_back = 200
        open(file_id_reversibility_test_back,file=filename_reversibility_test_back)
        !Write orbit positions
        do l = 1, n_snapshots
            !We save the backwards integrated quantities in reversed order to be comparable to forward integrated ones
            j = 1 + (n_steps - (l-1)*delta_snapshot)
            do k = 1,n_orbits+1
                write(file_id_reversibility_test_back,*) j, hamiltonian_time_back_mat(j,k), &
                    & x1_back_mat(j,k), x2_back_mat(j,k), x3_back_mat(j,k), &
                    & pitchpar_back_mat(j,k), gyrophase_back_mat(j,k)/(2.0d0*pi), &
                    kin_energy_eV_back_mat(j,k)
            enddo
            write(file_id_reversibility_test_back,*) snapshot_end
        enddo
        !The first orbit position is printed last in case of backward integration
        j = 1
        do k = 1,n_orbits+1
            write(file_id_reversibility_test_back,*) j, hamiltonian_time_back_mat(j,k), &
                & x1_back_mat(j,k), x2_back_mat(j,k), x3_back_mat(j,k), &
                & pitchpar_back_mat(j,k), gyrophase_back_mat(j,k)/(2.0d0*pi), &
                kin_energy_eV_back_mat(j,k)
        enddo
        write(file_id_reversibility_test_back,*) snapshot_end
        close(file_id_reversibility_test_back)
!
        deallocate(x1_mat,x2_mat,x3_mat,& 
                    & pitchpar_mat,kin_energy_eV_mat, & 
                    & hamiltonian_time_mat,gyrophase_mat)
        deallocate(x1_back_mat,x2_back_mat,x3_back_mat, & 
                    & pitchpar_back_mat,kin_energy_eV_back_mat, & 
                    & hamiltonian_time_back_mat,gyrophase_back_mat)
        deallocate(x1_vec,x2_vec,x3_vec,pitchpar_vec,kin_energy_eV_vec,hamiltonian_time_vec,gyrophase_vec)
        deallocate(alpha_0_vec)
!
    end subroutine make_reversibility_test
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine gorilla_integration_setup(lambda,t_total_sign, &
                                        & x,vpar,energy_eV,ind_tetr,iface,perpinv,perpinv2)
!
        use constants, only: ev2erg
        use tetra_physics_mod, only: tetra_physics, particle_mass, particle_charge
        use orbit_timestep_gorilla_mod, only: check_coordinate_domain
        use supporting_functions_mod, only: bmod_func
        use pusher_tetra_rk_mod, only: find_tetra
!
        implicit none
!
        double precision, intent(in)                    :: lambda
        integer, intent(in)                             :: t_total_sign
!
        double precision, dimension(3), intent(inout)   :: x
        double precision, intent(inout)                 :: energy_eV
!
        integer, intent(out)                            :: ind_tetr, iface
        double precision, intent(out)                   :: vpar, perpinv, perpinv2
!
        double precision                                :: vmod, vperp, vperp2 
        double precision, dimension(3)                  :: z_save
!
if(boole_diag_reversibility_test) print*, energy_eV
if(boole_diag_reversibility_test) print*, ev2erg
if(boole_diag_reversibility_test) print*, particle_mass
        !Compute velocity module from energy (ignoring potential energy -> only kinetic energy) dependent on particle species
        vmod=sqrt(2.d0*energy_eV*ev2erg/particle_mass)
!
        !--- Find tetrahedron for starting positions by neglecting electrostatic potential energy
        vpar = lambda*vmod
        vperp = sqrt((1.d0-lambda**2))*vmod
!
        !Check coordinate domain (optionally perform modulo operation)
        call check_coordinate_domain(x)
!
        !Find tetrahedron index and face index for position x
        call find_tetra(x,vpar,vperp,ind_tetr,iface,sign_t_step_in=t_total_sign)
!
        !If particle doesn't lie inside any tetrahedron
        if(ind_tetr.eq.-1) then
            print *, 'Particle position not found'
            return
        endif
!
        ! Compute perpinv
        vperp2 = vperp**2
        z_save = x-tetra_physics(ind_tetr)%x1
        perpinv= -0.5d0*vperp2/bmod_func(z_save,ind_tetr)
        perpinv2 = perpinv**2

        !Update Energy to include electrical potential
if(boole_diag_reversibility_test) print*, 'kinetic energy', energy_eV
        energy_eV = E_tot_eV_func(x,vpar,perpinv,ind_tetr)
!
if(boole_diag_reversibility_test) print*, 'total energy', energy_eV
if(boole_diag_reversibility_test) print *, 'vpar', vpar
if(boole_diag_reversibility_test) print *, 'vperp', vperp
if(boole_diag_reversibility_test) print *, 'vmod', vmod
if(boole_diag_reversibility_test) print *, 'lambda', vpar/vmod_func_new(energy_eV,x,ind_tetr)
!
    end subroutine gorilla_integration_setup
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine gorilla_integration(perpinv,perpinv2,t_total,n_steps, &
                                & x,vpar,energy_eV,ind_tetr,iface, &
                                & x1_vec,x2_vec,x3_vec,pitchpar_vec,kin_energy_eV_vec,hamiltonian_time_vec,gyrophase_vec, &
                                & starting_values)
!
        use constants, only: ev2erg
        use tetra_physics_mod, only: tetra_physics,particle_mass,particle_charge
        use pusher_tetra_rk_mod, only: pusher_tetra_rk,initialize_const_motion_rk
        use pusher_tetra_poly_mod, only: pusher_tetra_poly,initialize_const_motion_poly
        use gorilla_settings_mod, only: ipusher, poly_order, optional_quantities_type
!
        implicit none
!
        double precision, intent(in)                            :: perpinv, perpinv2, t_total
        integer, intent(in)                                     :: n_steps
        type(starting_values_type), optional, intent(in)        :: starting_values
!
        double precision, dimension(3), intent(inout)           :: x
        double precision, intent(inout)                         :: vpar, energy_eV
        integer, intent(inout)                                  :: ind_tetr, iface
!
        double precision,dimension(:),intent(out)               :: x1_vec, x2_vec, x3_vec, pitchpar_vec, &
                                                                & kin_energy_eV_vec, &
                                                                & hamiltonian_time_vec, gyrophase_vec
!
        double precision                                        :: t_step, vmod
        integer                                                 :: i, iper
        double precision, dimension(3)                          :: z_save
        double precision                                        :: t_remain, t_pass
        logical                                                 :: boole_t_finished
!
        double precision                                        :: t_hamiltonian, gyrophase
        type(optional_quantities_type)                          :: optional_quantities        
!
!------------------------------------------------------------------------------------------------------------!
! Precomputations and settings
!
        if(present(starting_values)) then
            t_hamiltonian = starting_values%t_hamiltonian
            gyrophase = starting_values%gyrophase
        else
            t_hamiltonian = 0.d0
            gyrophase = 0.d0
        endif
!
        !Define time step
        t_step = t_total/dble(n_steps)

        !Initialize constants of motion in particle-private module
        select case(ipusher)
            case(1)
                call initialize_const_motion_rk(perpinv,perpinv2)
            case(2)
                call initialize_const_motion_poly(perpinv,perpinv2)
        end select
!
!------------------------------------------------------------------------------------------------------------!
! Integrate orbit for time steps
!
        do i = 1, n_steps
!
            !Integrate particle orbit for given time step
            t_remain = t_step
!
            !Logical for handling time integration
            boole_t_finished = .false.
!
            !Output quantities (of the previous time step)
            x1_vec(i) = x(1)
            x2_vec(i) = x(2)
            x3_vec(i) = x(3)
            energy_eV = E_tot_eV_func(x,vpar,perpinv,ind_tetr)
            vmod=vmod_func_new(energy_eV,x,ind_tetr)
            pitchpar_vec(i) = vpar/vmod
            kin_energy_eV_vec(i) = particle_mass*vmod**2/2/ev2erg
            hamiltonian_time_vec(i) = t_hamiltonian
            gyrophase_vec(i) = gyrophase
!
            !Compute relative particle position
            z_save = x-tetra_physics(ind_tetr)%x1
!
            !Loop for tetrahedron pushings until t_step is reached
            do
!
                !Domain Boundary
                if(ind_tetr.eq.-1) then
                    print *, 'WARNING: Particle lost.'
if(boole_diag_reversibility_test) print *, x,vpar
if(boole_diag_reversibility_test) print *, 't_step',i
if(boole_diag_reversibility_test) stop
                    exit
                endif
!
                !t_remain (in) ... remaining time until t_step is finished
                !t_pass (out) ... time to pass the tetrahdron
!
                !Calculate trajectory
                select case(ipusher)
                    case(1)
                        call pusher_tetra_rk(ind_tetr,iface,x,vpar,z_save,t_remain,t_pass,boole_t_finished,iper)
                    case(2)
                        call pusher_tetra_poly(poly_order,ind_tetr,iface,x,vpar,z_save,t_remain,&
                                                         & t_pass,boole_t_finished,iper,optional_quantities)
                end select
!
                t_remain = t_remain - t_pass
!
                t_hamiltonian = t_hamiltonian + optional_quantities%t_hamiltonian
                gyrophase = gyrophase + optional_quantities%gyrophase
!
                !Orbit stops within cell, because "flight"-time t_step has finished
                if(boole_t_finished) then
                    exit
                endif
!
            enddo !Loop for tetrahedron pushings
!
        end do !i t_steps
!
        !Output quantities (of the last time step)
        x1_vec(i) = x(1)
        x2_vec(i) = x(2)
        x3_vec(i) = x(3)
        energy_eV = E_tot_eV_func(x,vpar,perpinv,ind_tetr)
        vmod = vmod_func_new(energy_eV,x,ind_tetr)
        pitchpar_vec(i) = vpar/vmod
        kin_energy_eV_vec(i) = particle_mass*vmod**2/2/ev2erg
        hamiltonian_time_vec(i) = t_hamiltonian
        gyrophase_vec(i) = gyrophase
!
    end subroutine gorilla_integration
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!     function phi_elec_func_new(x,ind_tetr)
! !
!         use tetra_physics_mod, only: tetra_physics
!         use supporting_functions_mod, only: phi_elec_func
! !
!         implicit none
! !
!         double precision :: phi_elec_func_new
!         integer, intent(in) :: ind_tetr
!         double precision, dimension(3),intent(in) :: x
!         double precision, dimension(3) :: z
! !
!         z = x-tetra_physics(ind_tetr)%x1
!         phi_elec_func_new = phi_elec_func(z,ind_tetr)
! !
!     end function phi_elec_func_new
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    function E_tot_eV_func(x,vpar,perpinv,ind_tetr)
!
        use constants, only: ev2erg
        use tetra_physics_mod, only: tetra_physics
        use supporting_functions_mod, only: energy_tot_func
!
        implicit none
!
        double precision                            :: E_tot_eV_func
!
        integer, intent(in)                         :: ind_tetr
        double precision, dimension(3),intent(in)   :: x
        double precision, intent(in)                :: perpinv,vpar
!
        double precision, dimension(3)              :: z
!
        z = x-tetra_physics(ind_tetr)%x1
!
        E_tot_eV_func = energy_tot_func([z,vpar],perpinv,ind_tetr)
        E_tot_eV_func = E_tot_eV_func/ev2erg
!
    end function E_tot_eV_func
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    function vmod_func_new(energy_eV,x,ind_tetr)
!
        use constants, only: ev2erg
        use tetra_physics_mod, only: tetra_physics
        use supporting_functions_mod, only: vmod_func
!
        implicit none
!
        double precision                            :: vmod_func_new
!
        double precision, intent(in)                :: energy_eV
        double precision, dimension(3),intent(in)   :: x
        integer, intent(in)                         :: ind_tetr
!
        double precision, dimension(3)              :: z123
        double precision                            :: energy
!
        z123 = x-tetra_physics(ind_tetr)%x1
        energy = energy_eV*ev2erg
        vmod_func_new = vmod_func(energy,z123,ind_tetr)
!
    end function vmod_func_new
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    function starting_position(alpha_0)
!
        use gorilla_settings_mod, only: coord_system
        use reversibility_test_load_mod, only: s_0, ds_0, theta_0, dtheta_0, R_0, dR_0, Z_0, dZ_0, phi_0
!
        implicit none
!
        double precision, dimension(3)  :: starting_position
!
        double precision                :: alpha_0
!
        select case(coord_system)
        case(1)
            starting_position(1) = R_0 + cos(alpha_0) * dR_0
            starting_position(2) = phi_0
            starting_position(3) = Z_0 + sin(alpha_0) * dZ_0
        case(2)
            starting_position(1) = s_0 + ds_0*cos(alpha_0)
            starting_position(2) = theta_0 + dtheta_0*sin(alpha_0)
            starting_position(3) = phi_0
        end select
!
    end function starting_position
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    function get_rand_noise_vec()
!
        use reversibility_test_load_mod, only: n_orbits, seed_option
!
        implicit none
!
        double precision, dimension(:), allocatable         :: get_rand_noise_vec
!
        double precision, dimension(:), allocatable         :: rand_noise_vec
        integer, dimension(:), allocatable                   :: seed
        double precision, dimension(:), allocatable         :: rd_seed
        integer                                             :: n
!
        allocate(rand_noise_vec(n_orbits+1),get_rand_noise_vec(n_orbits+1))
!
        !Seed options for setting the rand noises
        select case(seed_option)
            case(1) !No fixed seed
                n = 0
                call random_seed(size=n)
                allocate(seed(n))
                allocate(rd_seed(n))
                call random_number(rd_seed)
                seed = int(rd_seed*10.d0)
                deallocate(rd_seed)
                open(85,file='seed.inp')
                write(85,*) n
                write(85,*) seed
                close(85)
            case(2) !Load seed
                open(unit = 85, file='seed.inp', status='old',action = 'read')
                read(85,*) n
                allocate(seed(n))
                read(85,*) seed
                close(85)
        end select
!
        CALL RANDOM_SEED (PUT=seed)
        deallocate(seed)
!
        call random_number(rand_noise_vec)
        get_rand_noise_vec = (rand_noise_vec-0.5d0)
!
    end function get_rand_noise_vec
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine apply_noise(rand_noise,t_total_sign,x,vpar,energy_eV,ind_tetr,iface,perpinv,perpinv2)
!
        use tetra_physics_mod, only: tetra_physics
        use pusher_tetra_rk_mod, only: find_tetra
        use orbit_timestep_gorilla_mod, only: check_coordinate_domain
        use supporting_functions_mod, only: bmod_func
        use reversibility_test_load_mod, only: noise_amplitude
!
        implicit none
!
        integer, intent(in)                             :: t_total_sign
        double precision, intent(in)                    :: rand_noise
!
        double precision, dimension(:), intent(inout)   :: x
        double precision, intent(inout)                 :: vpar, energy_eV
        integer, intent(inout)                          :: ind_tetr
!
        integer, intent(out)                            :: iface
        double precision, intent(out)                   :: perpinv, perpinv2
!
        double precision                                :: vmod, vperp_disturbed, vperp2
        double precision, dimension(3)                  :: x_disturbed, z_save
        double precision                                :: vpar_disturbed
        integer                                         :: ind_tetr_disturbed, iface_disturbed
!
if(boole_diag_reversibility_test) print *, '---------------------------------------'     
if(boole_diag_reversibility_test) print *, 'energy', energy_eV
if(boole_diag_reversibility_test) print *, 'x', x
if(boole_diag_reversibility_test) print *, 'vpar', vpar
if(boole_diag_reversibility_test) print *, 'ind_tetr', ind_tetr
if(boole_diag_reversibility_test) print *, 'iface', iface
if(boole_diag_reversibility_test) print *, 'lambda', vpar/vmod_func_new(energy_eV,x,ind_tetr)
!
        !Get the UNTOUCHED vmod (aka before applying noise)
        vmod = vmod_func_new(energy_eV,x,ind_tetr)
!
        !Apply noise
        x_disturbed = x * (1 + noise_amplitude*rand_noise)
        vpar_disturbed = vpar * (1 + noise_amplitude*rand_noise)
        if (vpar_disturbed**2 .ge. vmod**2) vpar_disturbed = vmod !Safety boudary
!
        !Make other quantities consistent with disturbed coordinates (using the undisturbed vmod therefore leaving E_kin the same)
        vperp_disturbed = sqrt(vmod**2-vpar_disturbed**2)
!
        !Find tetrahedron index and face index for disturbed coordinate set
        call check_coordinate_domain(x_disturbed)
        call find_tetra(x_disturbed,vpar_disturbed,vperp_disturbed,ind_tetr_disturbed,iface_disturbed,sign_t_step_in=t_total_sign)
!
        !If particle doesn't lie inside any tetrahedron after applying noise
        if(ind_tetr.eq.-1) then
            print *, 'Noise could not be applied. Used undisturbed coordinates instead.'
            return
        endif
!
        !Save the disturbed quantities to the output
        x = x_disturbed
        vpar = vpar_disturbed
        ind_tetr = ind_tetr_disturbed
        iface = iface_disturbed
!
        ! Compute perpinv with saved disturbed quantities
        vperp2 = vperp_disturbed**2
        z_save = x-tetra_physics(ind_tetr)%x1
        perpinv= -0.5d0*vperp2/bmod_func(z_save,ind_tetr)
        perpinv2 = perpinv**2

        !Update Energy to include electrical potential (now also including disturbed quantities)
        energy_eV = E_tot_eV_func(x,vpar,perpinv,ind_tetr)
!
if(boole_diag_reversibility_test) print *, 'disturbed energy', energy_eV
if(boole_diag_reversibility_test) print *, 'disturbed x', x
if(boole_diag_reversibility_test) print *, 'disturbed vpar', vpar
if(boole_diag_reversibility_test) print *, 'disturbed ind_tetr', ind_tetr
if(boole_diag_reversibility_test) print *, 'disturbed iface', iface
if(boole_diag_reversibility_test) print *, 'disturbed lambda', vpar/vmod_func_new(energy_eV,x,ind_tetr)
!
    end subroutine apply_noise
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
end module reversibility_test_mod
