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
    double precision, public, protected    :: minor_0
    double precision, public, protected    :: dminor_0
    double precision, public, protected    :: R_center
    double precision, public, protected    :: Z_center
    double precision, public, protected    :: phi_0
    double precision, public, protected    :: pitchpar_0
!
    double precision, public, protected    :: energy_eV_0
    double precision, public, protected    :: relative_bandwith
!
    integer, public, protected             :: n_steps
    integer, public, protected             :: n_orbits
    integer, public, protected             :: n_snapshots
    character(50), public, protected       :: filename_reversibility_test
!
    logical, public, protected             :: boole_diag_reversibility_test
!
    !Namelist for reversibility test input
    NAMELIST /reversibility_test_nml/ t_total, minor_0, dminor_0, R_center, Z_center, phi_0, pitchpar_0, &
                                        & energy_eV_0, relative_bandwith, n_steps, n_orbits, n_snapshots, &
                                        & filename_reversibility_test, boole_diag_reversibility_test
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
    contains

    subroutine make_reversibility_test()
!
        use constants, only: pi
        use orbit_timestep_gorilla_mod, only: initialize_gorilla
        use gorilla_settings_mod, only: coord_system, boole_time_Hamiltonian, boole_gyrophase
        use reversibility_test_load_mod, only: t_total, minor_0, dminor_0, R_center, Z_center, phi_0, pitchpar_0, &
                                            & n_steps, n_orbits, n_snapshots, energy_eV_0, relative_bandwith, &
                                            & filename_reversibility_test
!
        implicit none
!
        double precision, dimension(3)                  :: x_0
        double precision, dimension(:), allocatable     :: alpha_0_vec
        double precision                                :: alpha_0, perpinv_0, perpinv, perpinv2, energy_eV
!
        double precision,dimension(:),allocatable       :: x1_vec, x2_vec, x3_vec, pitchpar_vec, &
                                                        & energy_eV_vec, perpinv_vec, &
                                                        & hamiltonian_time_vec, gyrophase_vec
        double precision, dimension(:,:), allocatable   :: x1_mat, x2_mat, x3_mat, pitchpar_mat, & 
                                                        & hamiltonian_time_mat, gyrophase_mat
!
        integer                                         :: i,j,k,l, counter_particles, delta_snapshot
        integer                                         :: file_id_reversibility_test
        double precision, dimension(7), parameter       :: snapshot_end = -1
!
!------------------------------------------------------------------------------------------------------------!
! Loading settings
!
        call load_reversibility_test_inp()    
!
!------------------------------------------------------------------------------------------------------------!
! Allocations and preparations
!
        !Allocate matrices to save complete data
        allocate(x1_mat(n_steps+1,n_orbits+1), x2_mat(n_steps+1,n_orbits+1), x3_mat(n_steps+1,n_orbits+1), & 
            & pitchpar_mat(n_steps+1,n_orbits+1), &
            & hamiltonian_time_mat(n_steps+1,n_orbits+1), gyrophase_mat(n_steps+1,n_orbits+1))
!
        !Allocate vectors to save individual orbit data
        allocate(x1_vec(n_steps+1), x2_vec(n_steps+1), x3_vec(n_steps+1), &
                & pitchpar_vec(n_steps+1), &
                & hamiltonian_time_vec(n_steps+1), gyrophase_vec(n_steps+1))
!
        !Create vector with equidistant theta values
        allocate(alpha_0_vec(n_orbits+1),energy_eV_vec(n_orbits+1),perpinv_vec(n_orbits+1))
        alpha_0_vec(:) = [(2.d0*pi/dble(n_orbits)*dble(i), i = 0,n_orbits)]
!
        !Initialize GORILLA
        call initialize_gorilla()
!
!------------------------------------------------------------------------------------------------------------!
! Compute orbits
!
        counter_particles = 0
!
        !Compute perpendicular invariant (constant for all orbits)
!
        !Orbit starting position for 1st orbit
        alpha_0 = alpha_0_vec(1)
        select case(coord_system)
            case(1)
                x_0(1) = R_center + cos(alpha_0)*(minor_0 + dminor_0*cos(alpha_0))
                x_0(2) = phi_0
                x_0(3) = Z_center + sin(alpha_0)*(minor_0 + dminor_0*cos(alpha_0))
            case(2)
                x_0(1) = minor_0 + dminor_0*cos(alpha_0)
                x_0(2) = sin(alpha_0)
                x_0(3) = phi_0
        end select
!
        call compute_perpinv_gorilla(x_0,energy_eV_0,pitchpar_0,perpinv_0)
!
!
        !$OMP PARALLEL DEFAULT(NONE) &
        !$OMP& SHARED(n_orbits,alpha_0_vec,t_total,n_steps, &
        !$OMP& minor_0,dminor_0,R_center,Z_center,phi_0,pitchpar_0,energy_eV_0,relative_bandwith, &
        !$OMP& x1_mat,x2_mat,x3_mat,pitchpar_mat,counter_particles, &
        !$OMP& perpinv_0,energy_eV_vec,perpinv_vec, &
        !$OMP& hamiltonian_time_mat,gyrophase_mat,coord_system) &
        !$OMP& PRIVATE(k,x_0,alpha_0,perpinv,perpinv2, &
        !$OMP& x1_vec,x2_vec,x3_vec,pitchpar_vec,energy_eV, &
        !$OMP& hamiltonian_time_vec,gyrophase_vec)
        !$OMP DO
!
        !Loop over orbits and store data
        do k = 1,n_orbits+1
!
            !orbit starting positions
            alpha_0 = alpha_0_vec(k)
            select case(coord_system)
                case(1)
                    x_0(1) = R_center + cos(alpha_0)*(minor_0 + dminor_0*cos(alpha_0))
                    x_0(2) = phi_0
                    x_0(3) = Z_center + sin(alpha_0)*(minor_0 + dminor_0*cos(alpha_0))
                case(2)
                    x_0(1) = minor_0 + dminor_0*cos(alpha_0)
                    x_0(2) = sin(alpha_0)
                    x_0(3) = phi_0
            end select
!
            !Define total energy (Hamiltonian) for individual orbit
            energy_eV = energy_eV_0 * (1.d0 + relative_bandwith * sin(alpha_0) )
            energy_eV_vec(k) = energy_eV
!
            !Define perpinv for individual orbit
            perpinv = perpinv_0 * (1.d0 + relative_bandwith * cos(alpha_0) )
            perpinv2 = perpinv**2
            perpinv_vec(k) = perpinv
!
            ! Compute orbits
            call gorilla_integration(x_0,pitchpar_0,perpinv,perpinv2,energy_eV,t_total,n_steps, &
                                & x1_vec,x2_vec,x3_vec,pitchpar_vec,hamiltonian_time_vec,gyrophase_vec)
!
            !Write values in matrix
            x1_mat(:,k) = x1_vec
            x2_mat(:,k) = x2_vec
            x3_mat(:,k) = x3_vec
            pitchpar_mat(:,k) = pitchpar_vec
            hamiltonian_time_mat(:,k) = hamiltonian_time_vec
            gyrophase_mat(:,k) = gyrophase_vec
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
            do k = 1,n_orbits
                write(file_id_reversibility_test,*) j, x1_mat(j,k), x2_mat(j,k), x3_mat(j,k), &
                    & pitchpar_mat(j,k), hamiltonian_time_mat(j,k), gyrophase_mat(j,k)
            enddo
            write(file_id_reversibility_test,*) snapshot_end
        enddo
        !Final orbit position
        do k = 1,n_orbits
            write(file_id_reversibility_test,*) j, x1_mat(n_steps+1,k), x2_mat(n_steps+1,k), x3_mat(n_steps+1,k), &
                & pitchpar_mat(n_steps+1,k), hamiltonian_time_mat(n_steps+1,k), gyrophase_mat(n_steps+1,k)
        enddo
        write(file_id_reversibility_test,*) snapshot_end
        close(file_id_reversibility_test)
!
        deallocate(x1_mat,x2_mat,x3_mat,pitchpar_mat,hamiltonian_time_mat,gyrophase_mat)
        deallocate(x1_vec,x2_vec,x3_vec,pitchpar_vec,hamiltonian_time_vec,gyrophase_vec)
        deallocate(alpha_0_vec,energy_eV_vec,perpinv_vec)
!
    end subroutine make_reversibility_test
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine gorilla_integration(x_0,pitchpar,perpinv,perpinv2,energy_eV,t_total,n_steps, &
                                & x1_vec,x2_vec,x3_vec,pitchpar_vec,hamiltonian_time_vec,gyrophase_vec)
!
        use constants, only: ev2erg
        use tetra_physics_mod, only: tetra_physics,particle_mass,particle_charge
        use orbit_timestep_gorilla_mod, only: bmod_func
        use pusher_tetra_rk_mod, only: pusher_tetra_rk,initialize_const_motion_rk
        use pusher_tetra_poly_mod, only: pusher_tetra_poly,initialize_const_motion_poly
        use gorilla_settings_mod, only: ipusher, poly_order, optional_quantities_type
!
        implicit none
!
        double precision, dimension(3), intent(in)  :: x_0
        double precision, intent(in)                :: pitchpar, perpinv, perpinv2, energy_eV, t_total
        integer, intent(in)                         :: n_steps
!
        double precision,dimension(:),intent(out)   :: x1_vec, x2_vec, x3_vec, pitchpar_vec, &
                                                    & hamiltonian_time_vec, gyrophase_vec
!
        double precision, dimension(3)              :: x
        double precision                            :: t_step, vmod, vpar, vperp
        integer                                     :: ind_tetr, iface,i, ind_tetr_save, iper
        double precision, dimension(3)              :: z_save
        double precision                            :: vperp2, t_remain, t_pass
        logical                                     :: boole_t_finished
        integer                                     :: counter_tor_mappings
!
        double precision                            :: t_hamiltonian, gyrophase
        type(optional_quantities_type)              :: optional_quantities        
!
!------------------------------------------------------------------------------------------------------------!
! Precomputations and settings
!
        !Hamiltonian time
        t_hamiltonian = 0.d0
!
        !Gyrophase
        gyrophase = 0.d0
!
        !Define time step
        t_step = t_total/dble(n_steps)
!
        !Define start position
        x = x_0
!
        call find_tetra_ePhi(x,energy_eV,pitchpar,ind_tetr,iface,vmod=vmod)
        if(ind_tetr.eq.-1) return
!
!------------------------------------------------------------------------------------------------------------!
! Use perpinv to re-compute vpar and vperp
!
        !Initialize constants of motion in particle-private module
        select case(ipusher)
            case(1)
                call initialize_const_motion_rk(perpinv,perpinv2)
            case(2)
                call initialize_const_motion_poly(perpinv,perpinv2)
        end select
!
        z_save = x-tetra_physics(ind_tetr)%x1
        vperp2 = -2.d0 * perpinv * bmod_func(z_save,ind_tetr)
!
        if(vperp2.ge.vmod**2) then
            print *, 'Please, choose a larger pitch parameter for starting orbit.'
            stop
        endif
!
        !Define vpar according to perpendicular invariant
        vpar = sqrt(vmod**2-vperp2)
!
if(boole_diag_reversibility_test) then
print *, 'vpar', vpar
print *, 'vperp', sqrt(vperp2)
print *, 'vmod', vmod
counter_tor_mappings = 0
endif
!
!------------------------------------------------------------------------------------------------------------!
! Integrate orbit for time steps
!
        !Save the tetrahedron index for computation of one-form at first step
        ind_tetr_save = ind_tetr
!
        do i = 1, n_steps+1
!
if(boole_diag_reversibility_test) then
    print *, 'step', i, '/', n_steps+1
endif
!
            !Integrate particle orbit for given time step
            t_remain = t_step
!
            !Logical for handling time integration
            boole_t_finished = .false.
!
            !Output quantities
            x1_vec(i) = x(1)
            x2_vec(i) = x(2)
            x3_vec(i) = x(3)
!
            !Compute relative particle position
            z_save = x-tetra_physics(ind_tetr_save)%x1
!
            vmod=sqrt(2.d0* (energy_eV*ev2erg - particle_charge * phi_elec_func(x,ind_tetr_save) ) / particle_mass)
            pitchpar_vec(i) = vpar/vmod
!
            hamiltonian_time_vec(i) = t_hamiltonian
            gyrophase_vec(i) = gyrophase
!
            !Loop for tetrahedron pushings until t_step is reached
            do
!
                !Domain Boundary
                if(ind_tetr.eq.-1) then
                    print *, 'WARNING: Particle lost.'
if(boole_diag_reversibility_test) then
print *, x,vpar
print *, 't_step',i
stop
endif
                    exit
                endif
!
                !Save the tetrahedron index
                ind_tetr_save = ind_tetr
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
if(boole_diag_reversibility_test) then
if(iper.ne.0) then
counter_tor_mappings = counter_tor_mappings + 1
endif
endif
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
if(boole_diag_reversibility_test) then
print *, 'counter_tor_mappings',counter_tor_mappings
print *, 't_hamiltonian',t_hamiltonian
endif
!
    end subroutine gorilla_integration
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine compute_perpinv_gorilla(x,energy_eV,pitchpar,perpinv)
!
        use tetra_physics_mod, only: tetra_physics
        use orbit_timestep_gorilla_mod, only: bmod_func
!
        implicit none
!
        double precision, intent(in)                    :: energy_eV,pitchpar
!
        double precision, dimension(3), intent(inout)   :: x
!
        double precision, intent(out)                   :: perpinv
!
        double precision                                :: vperp,vperp2
        integer                                         :: ind_tetr,iface
        double precision, dimension(3)                  :: z_save
!
        call find_tetra_ePhi(x,energy_eV,pitchpar,ind_tetr,iface,vperp=vperp)
        if(ind_tetr.eq.-1) return
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
    end subroutine compute_perpinv_gorilla
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    function phi_elec_func(x,ind_tetr)
!
        use tetra_physics_mod, only: tetra_physics
!
        implicit none
!
        double precision :: phi_elec_func
        integer, intent(in) :: ind_tetr
        double precision, dimension(3),intent(in) :: x
        double precision, dimension(3) :: z
!
        z = x-tetra_physics(ind_tetr)%x1
        phi_elec_func = tetra_physics(ind_tetr)%Phi1 + sum(tetra_physics(ind_tetr)%gPhi * z)
!
    end function phi_elec_func
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    function bmod_func_new(x,ind_tetr)
!
        use tetra_physics_mod, only: tetra_physics

        implicit none
!
        double precision :: bmod_func_new
        integer, intent(in) :: ind_tetr
        double precision, dimension(3),intent(in) :: x
        double precision, dimension(3) :: z
!
        z = x-tetra_physics(ind_tetr)%x1
        bmod_func_new = tetra_physics(ind_tetr)%bmod1+sum(tetra_physics(ind_tetr)%gb*z)
!
    end function bmod_func_new
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine find_tetra_ePhi(x,energy_eV,pitchpar,ind_tetr,iface,vmod,vperp)
!
        use constants, only: ev2erg
        use tetra_physics_mod, only: tetra_physics, particle_mass, particle_charge
        use orbit_timestep_gorilla_mod, only: check_coordinate_domain
        use pusher_tetra_rk_mod, only: find_tetra
!
        implicit none
!
        double precision, intent(in)                    :: energy_eV,pitchpar
!
        double precision, dimension(3), intent(inout)   :: x
!
        integer, intent(out)                            :: ind_tetr, iface
        double precision, intent(out), optional         :: vmod,vperp
!
        integer                                         :: ind_tetr_save
        double precision                                :: vpar
!
        !Compute velocity module from kinetic energy dependent on particle species
        vmod=sqrt(2.d0*energy_eV*ev2erg/particle_mass)
!
        !--- Find tetrahedron for starting positions by neglecting electrostatic potential energy
        vpar = pitchpar*vmod
        vperp = sqrt((1.d0-pitchpar**2))*vmod
!
        !Check coordinate domain (optionally perform modulo operation)
        call check_coordinate_domain(x)
!
        !Find tetrahedron index and face index for position x
        call find_tetra(x,vpar,vperp,ind_tetr,iface)
!
        !If particle doesn't lie inside any tetrahedron
        if(ind_tetr.eq.-1) then
            print *, 'Particle position not found'
            return
        endif
!
        !--- Take into account electrostatic potential energy
!
        !Compute velocity module from kinetic energy dependent on particle species
        vmod=sqrt(2.d0* (energy_eV*ev2erg - particle_charge * phi_elec_func(x,ind_tetr) ) / particle_mass)
!
        vpar = pitchpar*vmod
        vperp = sqrt((1.d0-pitchpar**2)*vmod**2)
!
        !Repeat find tetra (by taking into account electrostatic potential energy)
        ind_tetr_save = ind_tetr
        call find_tetra(x,vpar,vperp,ind_tetr,iface)
!
        if(ind_tetr.ne.ind_tetr_save) then
            print *, 'ERROR: Electrostatic potential energy term affects find_tetra'
            print *, 'E_PHI = ', particle_charge * phi_elec_func(x,ind_tetr_save) / ev2erg
            stop
        endif
!
    end subroutine find_tetra_ePhi
!
!------------------------------------------------------------------------------------------------------------!
!
end module reversibility_test_mod
