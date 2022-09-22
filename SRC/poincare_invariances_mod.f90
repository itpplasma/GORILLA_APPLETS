!
module poincare_invariances_mod
!
    !Energy electrons: 1.d-2
    !Energy D ions: 3.d3

    double precision :: energy_eV_0 = 3.d3
    double precision :: E_ePhi
    logical :: boole_diagnostic_pi = .false.
!
    contains
!
    subroutine compute_first_poincare_invariance()
!
        use constants, only: ev2erg, pi
        use tetra_physics_mod, only: particle_mass,cm_over_e,get_dt_dtau_const_1,set_all_dt_dtau_const
        use orbit_timestep_gorilla_mod, only: initialize_gorilla
        use tetra_grid_settings_mod, only: n_field_periods
        use gorilla_settings_mod, only: boole_time_Hamiltonian,boole_gyrophase
!
        implicit none
!
        integer :: n_orbits, i,j,k,l, i_integrator_type, n_steps, counter_particles, file_id_orbit_position,n_orbit_positions
        integer :: delta_orbit_position, file_id_orbit_cyl
        double precision :: s_0, theta_0, phi_0, pitchpar_0
        double precision, dimension(:), allocatable :: alpha_0_vec
        double precision :: J_perp, t_total, vmod, bmod, delta_theta, delta_phi,perpinv,perpinv2
        double precision :: A_theta,A_phi,dA_theta_ds,dA_phi_ds,aiota,     &
                            sqg,alam,dl_ds,dl_dt,dl_dp,Bctrvr_vartheta,Bctrvr_varphi,     &
                            Bcovar_r,Bcovar_vartheta,Bcovar_varphi
        double precision,dimension(:),allocatable :: s_vec,theta_vec,phi_vec,vpar_vec,A_theta_vec,A_phi_vec, &
            & h_s_vec, h_theta_vec, h_phi_vec, pitchpar_vec, hamiltonian_time_vec, energy_eV_vec, ePhi_vec, perpinv_vec, &
            & gyrophase_vec,bmod_vec
        double precision, dimension(:,:), allocatable :: s_mat,theta_mat,phi_mat,vpar_mat,A_theta_mat,A_phi_mat, &
            & h_s_mat, h_theta_mat, h_phi_mat, pitchpar_mat, hamiltonian_time_mat, ePhi_mat, gyrophase_mat, bmod_mat
        double precision, dimension(:), allocatable :: first_poincare_invariance, gyro_term
        double precision :: dt_dtau_const,dt_dtau_const_new
        double precision :: energy_eV
        logical :: boole_integration_over_tau
        double precision :: R_cyl, PHI_cyl, Z_cyl, delta_gyro_term
!
!------------------------------------------------------------------------------------------------------------!
! Settings
!
        !Total time evolution
        ! HYDRA ions (3keV);        t = 6.5d-5* 2.0d0
        ! TOk ions (3keV);          t = 6.5d-5* 2.0d0 / 4.d0
        ! TOK electrons (1E-2eV)    t = 6.5d-5* 2.0d0 / (4.d0*82.d0) * sqrt(3.d5)
        t_total = 6.5d-5* 2.0d0 / 4.d0 !save 2022 01 25
        !t_total = 2.5d-6
!
print *, ''
print *, 't_total (INPUT)',t_total
print *, 'boole_time_Hamiltonian',boole_time_Hamiltonian
print *, 'boole_gyrophase',boole_gyrophase
print *, ''
!
        !boole_integration_over_tau (ONLY for GORILLA)
        ! .true. ... tau is used as an integration variable
        ! .false. ... ordinary time is used as an integration variable
        boole_integration_over_tau = .true.
!
        !number of times steps
        n_steps = 100 * 2.d0
!
        !number of orbits
        n_orbits = 100000
!
        !Integrator type
        ! 1 ... GORILLA
        ! 2 ... DIRECT
        i_integrator_type = 1
!
        select case(i_integrator_type)
            case(1)
                open(2,file='results/poincare_invariance_gorilla_n_1E7_e_var_jperp_var_dtdtau_ham.dat')
            case(2)
                open(2,file='results/poincare_invariance_direct_1E-08_ns_s3.dat')
        end select
!
        !File ID for orbit position
        file_id_orbit_position = 100

        !File ID for cylindrical coordinates
        file_id_orbit_cyl = 200
!
        n_orbit_positions = 2
!
!------------------------------------------------------------------------------------------------------------!
! Allocations and preparations
!
        !Allocate matrices for computation of Poincaré invariance
        allocate(s_mat(n_steps+1,n_orbits+1),theta_mat(n_steps+1,n_orbits+1),phi_mat(n_steps+1,n_orbits+1), &
            & vpar_mat(n_steps+1,n_orbits+1), A_theta_mat(n_steps+1,n_orbits+1),A_phi_mat(n_steps+1,n_orbits+1), &
            & h_s_mat(n_steps+1,n_orbits+1), h_theta_mat(n_steps+1,n_orbits+1), h_phi_mat(n_steps+1,n_orbits+1), &
            & pitchpar_mat(n_steps+1,n_orbits+1),hamiltonian_time_mat(n_steps+1,n_orbits+1),ePhi_mat(n_steps+1,n_orbits+1), &
            & gyrophase_mat(n_steps+1,n_orbits+1),bmod_mat(n_steps+1,n_orbits+1))
!
        !Allocate output vectors for Poincaré invariance
        allocate(s_vec(n_steps+1),theta_vec(n_steps+1),phi_vec(n_steps+1), &
                & vpar_vec(n_steps+1),A_theta_vec(n_steps+1),A_phi_vec(n_steps+1), &
                & h_s_vec(n_steps+1), h_theta_vec(n_steps+1), h_phi_vec(n_steps+1), pitchpar_vec(n_steps+1), &
                & hamiltonian_time_vec(n_steps+1),ePhi_vec(n_steps+1), gyrophase_vec(n_steps+1), bmod_vec(n_steps+1))
!
        !Allocate first Poincaré invariance
        allocate(first_poincare_invariance(n_steps+1),gyro_term(n_steps+1))
!
        !Create vector with equidistant theta values
        allocate(alpha_0_vec(n_orbits+1),energy_eV_vec(n_orbits+1),perpinv_vec(n_orbits+1))
        alpha_0_vec(:) = [(2.d0*pi/dble(n_orbits)*dble(i), i = 0,n_orbits)]
!
        !Initialize GORILLA
        call initialize_gorilla()
!
        !Compute velocity module from kinetic energy dependent on particle species
        vmod=sqrt(2.d0*energy_eV_0*ev2erg/particle_mass)
!
!------------------------------------------------------------------------------------------------------------!
! Compute orbits and one-form elements for Poincaré invariances
!
        pitchpar_0 = 0.46d0 !0.6d0
        counter_particles = 0
!
        !Compute perpendicular invariant (constant for all orbits)
!
        !Orbit starting position for 1st orbit
        s_0 = 0.5d0+0.25d0*cos(alpha_0_vec(1))
        theta_0 = sin(alpha_0_vec(1))
        phi_0 = 0.d0
!
        select case(i_integrator_type)
!
            case(1) !GORILLA
!
                call compute_perpinv_gorilla(s_0, theta_0, phi_0, pitchpar_0,perpinv,perpinv2)
!
            case(2) !DIRECT
!
                call vmec_field(s_0,theta_0,phi_0,A_theta,A_phi,dA_theta_ds,dA_phi_ds,aiota,     &
                                sqg,alam,dl_ds,dl_dt,dl_dp,Bctrvr_vartheta,Bctrvr_varphi,     &
                                Bcovar_r,Bcovar_vartheta,Bcovar_varphi)
!
                bmod = sqrt(Bctrvr_vartheta*Bcovar_vartheta + Bctrvr_varphi*Bcovar_varphi)
!
                J_perp = (1.d0-pitchpar_0**2)*vmod**2 * particle_mass * cm_over_e / (2.d0*bmod)
!
        end select
!
!
        !Integration in orbit parameter - Manipulation of dt_dtau_const
        if(boole_integration_over_tau) then
        
            if(i_integrator_type.eq.1) then !ONLY for GORILLA
!
                !Get dt_dtau_const value of first vertex
                call get_dt_dtau_const_1(dt_dtau_const)
!
                !Set relation between orbit parameter and physical time to 1.d0
                t_total = t_total/dt_dtau_const
!
                dt_dtau_const_new = 1.d0
                call set_all_dt_dtau_const(dt_dtau_const_new)
!
            endif !ONLY for GORILLA
!
        endif !boole_integration_over_tau
!
        !$OMP PARALLEL DEFAULT(NONE) &
        !$OMP& SHARED(n_orbits,i_integrator_type,alpha_0_vec,particle_mass,cm_over_e,t_total,n_steps,ePhi_mat,perpinv_vec, &
        !$OMP& pitchpar_0,s_mat,theta_mat,phi_mat,vpar_mat,A_theta_mat,A_phi_mat, hamiltonian_time_mat,energy_eV_0, &
        !$OMP& h_s_mat, h_theta_mat, h_phi_mat,pitchpar_mat,counter_particles,J_perp,perpinv,perpinv2,energy_eV_vec, &
        !$OMP& gyrophase_mat,bmod_mat ) &
        !$OMP& PRIVATE(k,s_0,theta_0,phi_0,A_theta,A_phi,dA_theta_ds,dA_phi_ds,aiota, &
        !$OMP& sqg,alam,dl_ds,dl_dt,dl_dp,Bctrvr_vartheta,Bctrvr_varphi,Bcovar_r,Bcovar_vartheta,Bcovar_varphi, &
        !$OMP& bmod,s_vec,theta_vec,phi_vec,vpar_vec,A_theta_vec,A_phi_vec,h_s_vec, h_theta_vec, h_phi_vec,pitchpar_vec, &
        !$OMP& hamiltonian_time_vec,ePhi_vec,energy_eV,gyrophase_vec,bmod_vec)
        !$OMP DO
!
        !Loop over orbits and store data
        do k = 1,n_orbits+1
!
!if(k.eq.703) then
!    boole_diagnostic_pi = .true.
!endif
            !orbit starting positions
            s_0 = 0.5d0+0.25d0*cos(alpha_0_vec(k))
            theta_0 = sin(alpha_0_vec(k))
            phi_0 = 0.d0
!
            !Define total energy (Hamiltonian) for individual orbit
            energy_eV = energy_eV_0 * (1.d0 + 0.1d0 * sin(alpha_0_vec(k)) )
!            energy_eV = energy_eV_0 * (1.d0)
            energy_eV_vec(k) = energy_eV
!

            !Define perpinv for individual orbit
            perpinv_vec(k) = perpinv * (1.d0 + 0.1d0 * cos(alpha_0_vec(k)) )
!            perpinv_vec(k) = perpinv
            perpinv2 = perpinv_vec(k)**2

!
            ! Compute orbits
            select case(i_integrator_type)
!
                    case(1) !GORILLA
!
                        call gorilla_integration(s_0, theta_0, phi_0, pitchpar_0,perpinv_vec(k),perpinv2,energy_eV,t_total, &
                            & n_steps,s_vec,theta_vec,phi_vec,vpar_vec,A_theta_vec,A_phi_vec,h_s_vec, h_theta_vec, h_phi_vec, &
                            & pitchpar_vec, hamiltonian_time_vec,ePhi_vec,gyrophase_vec,bmod_vec)
!
                    case(2) !DIRECT
!
                        call direct_integration(s_0, theta_0, phi_0, J_perp, t_total, &
                            & n_steps,s_vec,theta_vec,phi_vec,vpar_vec,A_theta_vec,A_phi_vec,h_s_vec, h_theta_vec, h_phi_vec, &
                            & pitchpar_vec)
!
            end select
!
            !Write values in matrix
            s_mat(:,k) = s_vec
            theta_mat(:,k) = theta_vec
            phi_mat(:,k) = phi_vec
            vpar_mat(:,k) = vpar_vec
            A_theta_mat(:,k) = A_theta_vec
            A_phi_mat(:,k) = A_phi_vec
            h_s_mat(:,k) = h_s_vec
            h_theta_mat(:,k) = h_theta_vec
            h_phi_mat(:,k) = h_phi_vec
            pitchpar_mat(:,k) = pitchpar_vec
            hamiltonian_time_mat(:,k) = hamiltonian_time_vec
            ePhi_mat(:,k) = ePhi_vec
            gyrophase_mat(:,k) = gyrophase_vec
            bmod_mat(:,k) = bmod_vec
!
            !$omp critical
                counter_particles = counter_particles +1
                print *, 'Counter orbits', counter_particles, '/',n_orbits+1
            !$omp end critical

        enddo
!
        !$OMP END DO
        !$OMP END PARALLEL
!

!write(88,*) ''
!do k = 1,n_orbits+1
!
!write(88,*) hamiltonian_time_mat(n_steps+1,k)
!enddo
!stop
!------------------------------------------------------------------------------------------------------------!
! Compute Poincaré invariances via trapez formula for separate time steps
!
        counter_particles = 0
!
        !$OMP PARALLEL DEFAULT(NONE) &
        !$OMP& SHARED(n_orbits,particle_mass,cm_over_e,n_steps,boole_integration_over_tau,i_integrator_type, &
        !$OMP& s_mat,theta_mat,phi_mat,vpar_mat,A_theta_mat,A_phi_mat,hamiltonian_time_mat,n_field_periods, &
        !$OMP& h_s_mat, h_theta_mat, h_phi_mat,pitchpar_mat,counter_particles,energy_eV_vec,first_poincare_invariance, &
        !$OMP& gyrophase_mat,perpinv_vec,gyro_term) &
        !$OMP& PRIVATE(j,k,delta_theta,delta_phi,delta_gyro_term)
        !$OMP DO

        do j = 1,n_steps+1
            first_poincare_invariance(j) = 0.d0
            gyro_term(j) = 0.d0

            do k = 1,n_orbits
!
                ! Consider periodic boundaries
                delta_theta = theta_mat(j,k+1)-theta_mat(j,k)
                if(abs(delta_theta).gt.pi) then
                    if(delta_theta.gt.0.d0) then
!print *, 'delta theta is positive', delta_theta, 'j',j,'k',k
                        delta_theta = mod(theta_mat(j,k+1)-theta_mat(j,k)-2.d0*pi,2.d0*pi)
                    else
!print *, 'delta theta is negative', delta_theta, 'j',j,'k',k
                        delta_theta = mod(theta_mat(j,k+1)-theta_mat(j,k)+2.d0*pi,2.d0*pi)
                    endif
                endif
!                delta_theta = abs(delta_theta)
!
                delta_phi = (phi_mat(j,k+1)-phi_mat(j,k))
                if(abs(delta_phi).gt.(pi/n_field_periods)) then

                    if(delta_phi.gt.0.d0) then
!                    print *, 'delta phi is positive', delta_phi, 'j',j,'k',k
                        delta_phi = mod(phi_mat(j,k+1)-phi_mat(j,k)-2.d0*pi/n_field_periods,2.d0*pi/n_field_periods)
                    else
!                    print *, 'delta phi is negative', delta_phi, 'j',j,'k',k
                        delta_phi = mod(phi_mat(j,k+1)-phi_mat(j,k)+2.d0*pi/n_field_periods,2.d0*pi/n_field_periods)
                    endif
                    !delta_phi = mod(delta_phi+2.d0*pi/n_field_periods,2.d0*pi/n_field_periods)
                endif
!
!if(j.eq.120) then
!
!write(11,*) k,delta_theta,delta_phi
!write(12,*) k,theta_mat(j,k),phi_mat(j,k)
!
!endif
!
                ! Compute integral via trapez formula
                first_poincare_invariance(j) = first_poincare_invariance(j) + &
!
                    ! First term
                    & (vpar_mat(j,k)*cm_over_e*h_s_mat(j,k) + vpar_mat(j,k+1)*cm_over_e*h_s_mat(j,k+1))/2.d0 * &
                    & (s_mat(j,k+1)-s_mat(j,k)) * (particle_mass/cm_over_e) + &
!
                    ! Second term
                    & (A_theta_mat(j,k) + vpar_mat(j,k)*cm_over_e*h_theta_mat(j,k) + &
                    & A_theta_mat(j,k+1) + vpar_mat(j,k+1)*cm_over_e*h_theta_mat(j,k+1))/2.d0 * &
                    & delta_theta * (particle_mass/cm_over_e) + &
!
                    ! Third term
                    & (A_phi_mat(j,k) + vpar_mat(j,k)*cm_over_e*h_phi_mat(j,k) + &
                    & A_phi_mat(j,k+1) + vpar_mat(j,k+1)*cm_over_e*h_phi_mat(j,k+1))/2.d0 * &
                    & delta_phi * &
                    & (particle_mass/cm_over_e)
!
                    ! Term for gyrophase
                    delta_gyro_term =  - (perpinv_vec(k) + perpinv_vec(k+1))/2.d0 * (-particle_mass*cm_over_e) * &
                    & (gyrophase_mat(j,k+1) - gyrophase_mat(j,k) )
!
                    first_poincare_invariance(j) = first_poincare_invariance(j) + delta_gyro_term
                    gyro_term(j) = gyro_term(j) + delta_gyro_term
!
                    ! Extra term in the case of extended phase space
                    if(boole_integration_over_tau) then
                        if(i_integrator_type.eq.1) then !ONLY for GORILLA
                        first_poincare_invariance(j) = first_poincare_invariance(j) &
                            & - ev2erg * ( energy_eV_vec(k) + energy_eV_vec(k+1) )/2.d0 * &
                            & (hamiltonian_time_mat(j,k+1)-hamiltonian_time_mat(j,k))
                        endif
                    endif
!
            enddo
!
            !$omp critical
                counter_particles = counter_particles +1
                print *, 'Counter Poincaré invariant steps:', counter_particles, '/',n_steps+1
            !$omp end critical
!
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
!
        !Write first Poincaré invariance
        do j = 1,n_steps+1
            write(2,*) j,first_poincare_invariance(j),gyro_term(j)
        enddo
!
        close(2)
!
!
!

        delta_orbit_position = n_steps/n_orbit_positions
!
        j = 1
        !Write orbit positions
        do l = 1,n_orbit_positions
!
            do k = 1,n_orbits
                write(file_id_orbit_position,*) j,s_mat(j,k),theta_mat(j,k),phi_mat(j,k), &
                    & pitchpar_mat(j,k),energy_eV_vec(k),ePhi_mat(j,k),(-particle_mass*cm_over_e)*perpinv_vec(k), &
                    & gyrophase_mat(j,k), bmod_mat(j,k)

                call sym_flux_to_cyl(s_mat(j,k),theta_mat(j,k),phi_mat(j,k),R_cyl,PHI_cyl,Z_cyl)
                write(file_id_orbit_cyl,*) R_cyl,PHI_cyl,Z_cyl,pitchpar_mat(n_steps+1,k)
            enddo
!
            j = j+delta_orbit_position
            file_id_orbit_position = file_id_orbit_position + 1
            file_id_orbit_cyl = file_id_orbit_cyl + 1

        enddo
!
        !Final orbit position
        do k = 1,n_orbits
            write(file_id_orbit_position,*) j,s_mat(n_steps+1,k),theta_mat(n_steps+1,k),phi_mat(n_steps+1,k), &
                & pitchpar_mat(n_steps+1,k),energy_eV_vec(k),ePhi_mat(j,k),(-particle_mass*cm_over_e)*perpinv_vec(k), &
                & gyrophase_mat(n_steps+1,k), bmod_mat(j,k)

            call sym_flux_to_cyl(s_mat(n_steps+1,k),theta_mat(n_steps+1,k),phi_mat(n_steps+1,k),R_cyl,PHI_cyl,Z_cyl)
            write(file_id_orbit_cyl,*) R_cyl,PHI_cyl,Z_cyl,pitchpar_mat(n_steps+1,k)
        enddo

!
!Print electro-static potential energy
!print *, 'Electro-static potential energy'
!print *, E_ePhi


        
!        !Write out trajectory of a single orbit
!        l = 200
!        do k = 1,10
!            l = l+1
!            do j = 1,n_steps+1
!                write(l,*) hamiltonian_time_mat(j,k),gyrophase_mat(j,k),bmod_mat(j,k)
!            enddo
!        enddo



        deallocate(s_mat,theta_mat,phi_mat,vpar_mat,A_theta_mat,A_phi_mat, &
            & h_s_mat, h_theta_mat, h_phi_mat, pitchpar_mat)
!
        deallocate(s_vec,theta_vec,phi_vec,vpar_vec,A_theta_vec,A_phi_vec,h_s_vec, h_theta_vec, h_phi_vec, pitchpar_vec)
!
    end subroutine compute_first_poincare_invariance
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine direct_integration(s_0, theta_0, phi_0, J_perp, t_total, &
        & n_steps,s_vec,theta_vec,phi_vec,vpar_vec,A_theta_vec,A_phi_vec, &
        & h_s_vec, h_theta_vec, h_phi_vec, pitchpar_vec)
!
        use constants, only: ev2erg, pi
        use tetra_physics_mod, only: particle_mass,cm_over_e,mag_axis_R0
        use parmot_mod, only: rmu,ro0
        use velo_mod, only: isw_field_type
        use supporting_functions_mod, only: theta_sym_flux2theta_vmec,theta_vmec2theta_sym_flux
        use tetra_grid_settings_mod, only: n_field_periods
!
        implicit none
!
        double precision, intent(in) :: s_0, theta_0, phi_0, J_perp, t_total
        integer, intent(in) :: n_steps
        double precision,dimension(:),intent(out):: s_vec,theta_vec,phi_vec,vpar_vec,A_theta_vec,A_phi_vec, &
            & h_s_vec, h_theta_vec, h_phi_vec, pitchpar_vec
!
        integer :: i,ierr
        double precision :: dtau, dphi,dtaumin,time_step,tau_out_can
        double precision :: vmod
        double precision, dimension(5) :: z
!
        double precision :: bmod,vperp2,pitchpar
        double precision :: A_theta,A_phi,dA_theta_ds,dA_phi_ds,aiota,     &
                            sqg,alam,dl_ds,dl_dt,dl_dp,Bctrvr_vartheta,Bctrvr_varphi,     &
                            Bcovar_r,Bcovar_vartheta,Bcovar_varphi
!
        !Define time step
        time_step = t_total/dble(n_steps)
!
        !Compute velocity module from kinetic energy dependent on particle species
        vmod=sqrt(2.d0*energy_eV_0*ev2erg/particle_mass)
!
!------------------------------------------------------------------------------------------------------------!
! Initialization of direct integrator
!
        !inverse relativistic temperature
        rmu=1d18
!
        !normalized larmor radius
        ro0 = vmod*cm_over_e
!
        isw_field_type=1

        !normalized slowing down time:
        dtau = -1.d0*time_step*vmod
!
        !field line integration step step over phi (to check chamber wall crossing)
        dphi=2.d0*pi/n_field_periods
!
        !orbit integration time step (to check chamber wall crossing)
        dtaumin=2.d0*dphi*mag_axis_R0
!
!------------------------------------------------------------------------------------------------------------!
!
        !Guiding-center position and pitch parameter
        z(1) = s_0
        z(2) = theta_sym_flux2theta_vmec(s_0,theta_0,phi_0)
        !z(2) = theta_0
        z(3) = phi_0
        z(4) = 1.d0
!
        !Compute starting pitch parameter
!
        call vmec_field(z(1),z(2),z(3),A_theta,A_phi,dA_theta_ds,dA_phi_ds,aiota,     &
                        sqg,alam,dl_ds,dl_dt,dl_dp,Bctrvr_vartheta,Bctrvr_varphi,     &
                        Bcovar_r,Bcovar_vartheta,Bcovar_varphi)
!
        bmod = sqrt(Bctrvr_vartheta*Bcovar_vartheta + Bctrvr_varphi*Bcovar_varphi)
!
        vperp2 = 2.d0*J_perp*bmod/(cm_over_e*particle_mass)
        pitchpar = sqrt(1.d0-(vperp2/vmod**2))

        z(5) = pitchpar
!
        !Loop over time steps
!
        do i = 1,n_steps+1
!
            !Call field
            call vmec_field(z(1),z(2),z(3),A_theta,A_phi,dA_theta_ds,dA_phi_ds,aiota,     &
                            sqg,alam,dl_ds,dl_dt,dl_dp,Bctrvr_vartheta,Bctrvr_varphi,     &
                            Bcovar_r,Bcovar_vartheta,Bcovar_varphi)
!
            bmod = sqrt(Bctrvr_vartheta*Bcovar_vartheta + Bctrvr_varphi*Bcovar_varphi)
!
            !Output quantities
            s_vec(i) = z(1)
            theta_vec(i) = mod(theta_vmec2theta_sym_flux(z(1),mod(z(2),2.d0*pi),mod(z(3),2.d0*pi)),2.d0*pi)
            !theta_vec(i) = z(2)
            phi_vec(i) = mod(z(3),2.d0*pi/n_field_periods)
            vpar_vec(i) = z(4)*vmod
            A_theta_vec(i) = A_theta
            A_phi_vec(i) = A_phi
            h_s_vec(i) = Bcovar_r/bmod
            h_theta_vec(i) = Bcovar_vartheta/bmod
            h_phi_vec(i) = Bcovar_varphi/bmod
            pitchpar_vec(i) = z(5)
!
!write(3,*) i, z

            !Integrate guiding-center equations of motion
            call orbit_timestep_can(z,dtau,dtaumin,ierr,tau_out_can)
!
        enddo
!
    end subroutine direct_integration
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine gorilla_integration(s_0, theta_0, phi_0, pitchpar_0,perpinv,perpinv2,energy_eV,t_total, &
        & n_steps,s_vec,theta_vec,phi_vec,vpar_vec,A_theta_vec,A_phi_vec,h_s_vec, h_theta_vec, h_phi_vec, pitchpar_vec, &
        & hamiltonian_time_vec,ePhi_vec,gyrophase_vec,bmod_vec)
!
        use constants, only: ev2erg, pi
        use tetra_physics_mod, only: tetra_physics,particle_mass,hamiltonian_time,particle_charge,isinside
        use orbit_timestep_gorilla_mod, only: check_coordinate_domain,bmod_func
        use pusher_tetra_rk_mod, only: find_tetra,pusher_tetra_rk,initialize_const_motion_rk
        use pusher_tetra_poly_mod, only: pusher_tetra_poly,initialize_const_motion_poly
        use gorilla_settings_mod, only: ipusher, poly_order, optional_quantities_type
        use gorilla_diag_mod, only: diag_pusher_tetry_poly
!
        implicit none
!
        double precision, intent(in) :: s_0, theta_0, phi_0,pitchpar_0, perpinv,perpinv2,energy_eV,t_total
        integer, intent(in) :: n_steps
        double precision,dimension(:),intent(out) :: s_vec,theta_vec,phi_vec,vpar_vec,A_theta_vec,A_phi_vec, &
            & h_s_vec, h_theta_vec, h_phi_vec, pitchpar_vec,hamiltonian_time_vec,ePhi_vec,gyrophase_vec,bmod_vec
!
        double precision, dimension(3) :: x
        double precision :: t_step,vmod, vpar, vperp
        integer :: ind_tetr,iface,i,ind_tetr_save,iper
        double precision, dimension(3)                  :: z_save
        double precision                                :: vperp2,t_remain,t_pass,perpinv_temp,perpinv2_temp
        logical                                         :: boole_t_finished
        integer :: counter_tetra_pushings, counter_tor_mappings
!
        double precision :: t_hamiltonian,gyrophase
        integer :: l
        type(optional_quantities_type) :: optional_quantities        
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
        !Compute velocity module from kinetic energy dependent on particle species
        vmod=sqrt(2.d0*energy_eV*ev2erg/particle_mass)
!
        !Define start position
        x(1) = s_0
        x(2) = theta_0
        x(3) = phi_0
!
        !--- Find tetrahedron for starting positions by neglecting electrostatic potential energy
        vpar = pitchpar_0*vmod
        vperp = sqrt((1.d0-pitchpar_0**2))*vmod
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
        vpar = pitchpar_0*vmod
        vperp = sqrt((1.d0-pitchpar_0**2)*vmod**2)
!
        !Repeat find tetra (by taking into account electrostatic potential energy)
        ind_tetr_save = ind_tetr
        call find_tetra(x,vpar,vperp,ind_tetr,iface)
!
        if(ind_tetr.ne.ind_tetr_save) then
            print *, 'ERROR: Electrostatic potential energy term affects find_tetra'
            stop
        endif
!
!------------------------------------------------------------------------------------------------------------!
! Use perpinv to re-compute vpar and vperp
!
        !Compute relative particle position
        z_save = x-tetra_physics(ind_tetr)%x1

!vperp2 = vperp**2
!!
!!Compute perpendicular invariant of particle
!perpinv_temp=-0.5d0*vperp2/bmod_func(z_save,ind_tetr)
!perpinv2_temp = perpinv**2

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

if(boole_diagnostic_pi) then
print *, 'vpar', vpar
print *, 'vperp', sqrt(vperp2)
print *, 'vmod',vmod
endif

counter_tor_mappings = 0

!------------------------------------------------------------------------------------------------------------!
! Integrate orbit for time steps
!
        !Save the tetrahedron index for computation of one-form at first step
        ind_tetr_save = ind_tetr
!
        do i = 1,n_steps+1
!
if(boole_diagnostic_pi) then
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
            s_vec(i) = x(1)
            theta_vec(i) = x(2)
            phi_vec(i) = x(3)
            vpar_vec(i) = vpar

            !Compute relative particle position
            z_save = x-tetra_physics(ind_tetr_save)%x1

            !Use ind_tetr_save for obtaining appropriate quantities for Poincaré invariance
            A_theta_vec(i) = tetra_physics(ind_tetr_save)%Atheta1 + sum(tetra_physics(ind_tetr_save)%gAtheta*z_save)
            A_phi_vec(i) = tetra_physics(ind_tetr_save)%Aphi1 + sum(tetra_physics(ind_tetr_save)%gAphi*z_save)
            h_s_vec(i) = tetra_physics(ind_tetr_save)%h1_1 + sum(tetra_physics(ind_tetr_save)%gh1*z_save)
            h_theta_vec(i) = tetra_physics(ind_tetr_save)%h2_1 + sum(tetra_physics(ind_tetr_save)%gh2*z_save)
            h_phi_vec(i) = tetra_physics(ind_tetr_save)%h3_1 + sum(tetra_physics(ind_tetr_save)%gh3*z_save)
!
            vmod=sqrt(2.d0* (energy_eV*ev2erg - particle_charge * phi_elec_func(x,ind_tetr_save) ) / particle_mass)
            pitchpar_vec(i) = vpar/vmod
!
            hamiltonian_time_vec(i) = t_hamiltonian
            gyrophase_vec(i) = gyrophase
!
            ePhi_vec(i) = particle_charge * phi_elec_func(x,ind_tetr_save)/ev2erg

!if(.not.isinside(ind_tetr_save,x)) then
!    print *, 'NOT inside'
!    stop
!endif

            bmod_vec(i) = bmod_func_new(x,ind_tetr_save)
!
            !Loop for tetrahedron pushings until t_step is reached
counter_tetra_pushings = 0
!
            do
!
if(boole_diagnostic_pi) then
 print *, 'counter_pushings',counter_tetra_pushings
 diag_pusher_tetry_poly = .true.
endif
!
                !Domain Boundary
                if(ind_tetr.eq.-1) then
                    print *, 'WARNING: Particle lost.'
print *, x,vpar
print *, 't_step',i
print *, 'counter_tetra_pushings',counter_tetra_pushings
stop
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
!                t_hamiltonian = t_hamiltonian + t_pass * hamiltonian_time(ind_tetr_save)%dt_dtau_const_save
!
                gyrophase = gyrophase + optional_quantities%gyrophase
!
counter_tetra_pushings = counter_tetra_pushings + 1
if(iper.ne.0) then
    !print *, 'iper',iper
    counter_tor_mappings = counter_tor_mappings + 1
endif
!
                !Orbit stops within cell, because "flight"-time t_step has finished
                if(boole_t_finished) then
!
                    exit
                endif
!
            enddo !Loop for tetrahedron pushings
!
        end do !i t_steps
!
!print *, 'counter_tor_mappings',counter_tor_mappings
!print *, 't_hamiltonian',t_hamiltonian
!
    end subroutine gorilla_integration
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine compute_perpinv_gorilla(s_0, theta_0, phi_0, pitchpar_0,perpinv,perpinv2)
!
        use constants, only: ev2erg
        use tetra_physics_mod, only: tetra_physics,particle_mass,particle_charge
        use orbit_timestep_gorilla_mod, only: check_coordinate_domain,bmod_func
        use pusher_tetra_rk_mod, only: find_tetra
!
        implicit none
!
        double precision, intent(in) :: s_0, theta_0, phi_0,pitchpar_0
        double precision, intent(out) :: perpinv,perpinv2
        double precision, dimension(3) :: x
        double precision :: vmod, vpar, vperp,vperp2
        integer :: ind_tetr,iface,ind_tetr_save
        double precision, dimension(3)                  :: z_save
!
        !Compute velocity module from kinetic energy dependent on particle species
        vmod=sqrt(2.d0*energy_eV_0*ev2erg/particle_mass)
!
        !Define start position
        x(1) = s_0
        x(2) = theta_0
        x(3) = phi_0
!
        !--- Find tetrahedron for starting positions by neglecting electrostatic potential energy
!
        vpar = pitchpar_0*vmod
        vperp = sqrt((1.d0-pitchpar_0**2)*vmod**2)
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
        vmod=sqrt(2.d0* (energy_eV_0*ev2erg - particle_charge * phi_elec_func(x,ind_tetr) ) / particle_mass)
!
        vpar = pitchpar_0*vmod
        vperp = sqrt((1.d0-pitchpar_0**2)*vmod**2)
!
        !Repeat find tetra (by taking into account electrostatic potential energy)
        ind_tetr_save = ind_tetr
        call find_tetra(x,vpar,vperp,ind_tetr,iface)
!
        if(ind_tetr.ne.ind_tetr_save) then
            print *, 'ERROR: Electrostatic potential energy term affects find_tetra'
            stop
        endif
!
        !Squared perpendicular velocity
        vperp2 = vperp**2
!
        !Compute relative particle position
        z_save = x-tetra_physics(ind_tetr)%x1
!
        !Compute perpendicular invariant of particle
        perpinv=-0.5d0*vperp2/bmod_func(z_save,ind_tetr)
        perpinv2 = perpinv**2
!
E_ePhi = particle_charge * phi_elec_func(x,ind_tetr)/ev2erg
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
    subroutine sym_flux_to_cyl(x1,x2,x3,R,PHI,Z)
!
        use tetra_physics_mod, only: coord_system
        use tetra_grid_settings_mod, only: grid_kind
        use supporting_functions_mod, only: theta_sym_flux2theta_vmec
!
        implicit none
!
        double precision :: x1,x2,x3,R,PHI,Z
!
        !Variables for magdata in symflux coordinates
        integer :: inp_label
        double precision :: psi_pol,q,dq_ds,sqrtg,bmod,dbmod_dtheta,dR_ds,dR_dtheta
        double precision :: dZ_ds, dZ_dtheta
        double precision :: A_phi,A_theta,dA_phi_ds,dA_theta_ds,aiota,       &
                              alam,dR_dt,dR_dp,dZ_dt,dZ_dp,dl_ds,dl_dt,dl_dp
!
        select case(grid_kind)
            case(2) !EFIT field-aligned grid
                inp_label = 1
                call magdata_in_symfluxcoord_ext(inp_label,x1,psi_pol,x2,q,dq_ds, &
                                             sqrtg,bmod,dbmod_dtheta,R,dR_ds,dR_dtheta,       &
                                             Z,dZ_ds,dZ_dtheta)
            case(3) !VMEC field-aligned grid
!
                !Find VMEC theta for sym-flux-theta
                x2 = theta_sym_flux2theta_vmec(x1,x2,x3)
!
                call splint_vmec_data(x1,x2,x3,A_phi,A_theta,dA_phi_ds,dA_theta_ds,aiota,       &
                                    R,Z,alam,dR_ds,dR_dt,dR_dp,dZ_ds,dZ_dt,dZ_dp,dl_ds,dl_dt,dl_dp)
!
        end select !grid_kind
!
        PHI = x3
!
    end subroutine sym_flux_to_cyl
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
end module poincare_invariances_mod
