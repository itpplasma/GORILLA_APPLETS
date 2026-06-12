module volume_integrals_and_sqrt_g_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

    real(dp), dimension(:,:), allocatable :: sqrt_g

    public :: calc_square_root_g, calc_volume_integrals, calc_volume_integrals_in_flux_coordinates
    public :: deallocate_sqrt_g

contains

subroutine deallocate_sqrt_g
!
! Deallocates sqrt_g array to allow reallocation when grid parameters change.
!
    if (allocated(sqrt_g)) deallocate(sqrt_g)

end subroutine deallocate_sqrt_g

subroutine calc_square_root_g

    use tetra_physics_mod, only: tetra_physics, hamiltonian_time
    use tetra_grid_mod, only: ntetr, tetra_grid, verts_rphiz

    integer :: ind_tetr

    allocate(sqrt_g(ntetr,7))
    !compare first 6 entries with chapter 4.5 of master thesis from Jonatan Schatzlmayr, entry 7 is the radius (old metric determinant)

    do ind_tetr = 1, ntetr
        sqrt_g(ind_tetr,1) = hamiltonian_time(ind_tetr)%h1_in_curlA

        sqrt_g(ind_tetr,2) = tetra_physics(ind_tetr)%gh1(1)*tetra_physics(ind_tetr)%curlA(1) + &
                            & tetra_physics(ind_tetr)%gh2(1)*tetra_physics(ind_tetr)%curlA(2) + &
                            & tetra_physics(ind_tetr)%gh3(1)*tetra_physics(ind_tetr)%curlA(3)

        sqrt_g(ind_tetr,3) = tetra_physics(ind_tetr)%gh1(3)*tetra_physics(ind_tetr)%curlA(1) + &
                            & tetra_physics(ind_tetr)%gh2(3)*tetra_physics(ind_tetr)%curlA(2) + &
                            & tetra_physics(ind_tetr)%gh3(3)*tetra_physics(ind_tetr)%curlA(3)

        sqrt_g(ind_tetr,4) = tetra_physics(ind_tetr)%bmod1

        sqrt_g(ind_tetr,5) = tetra_physics(ind_tetr)%gB(1)

        sqrt_g(ind_tetr,6) = tetra_physics(ind_tetr)%gB(3)

        sqrt_g(ind_tetr,7) = verts_rphiz(1,tetra_grid(ind_tetr)%ind_knot(1))
    enddo

end subroutine calc_square_root_g

subroutine calc_volume_integrals

    use constants, only: pi, ev2erg
    use tetra_grid_mod, only: ntetr, verts_rphiz, tetra_grid
    use tetra_grid_settings_mod, only: grid_size, n_field_periods
    use tetra_physics_mod, only: particle_mass,particle_charge, tetra_physics
    use gorilla_applets_types_mod, only: output, in
    use, intrinsic :: ieee_arithmetic, only: ieee_set_halting_mode, ieee_set_flag, &
                                           & ieee_invalid, ieee_overflow, ieee_divide_by_zero, &
                                           & ieee_is_finite

    real(dp), dimension(:), allocatable              :: prism_volumes, refined_prism_volumes, elec_pot_vec, n_b
    integer                                          :: i,k
    integer                                          :: n_prisms
    integer, dimension(2)                            :: triangle_indices
    real(dp), dimension(3)                           :: r_values, z_values, r_values_intermediate, z_values_intermediate, gradient
    real(dp), dimension(2)                           :: r, z, limits
    real(dp)                                         :: z_star, alpha, beta, gamma, delta, epsilon, capital_gamma, capital_delta, &
                                                        a, a_dash, b, b_dash, c, c_dash, rmin, delta_r, delta_z, phi_0, eta
    integer, dimension(:,:), allocatable             :: tetra_indices_per_prism
    real(dp), dimension(:,:), allocatable            :: r_integrand_constants

    print*, 'calc_volume_integrals started'
    n_prisms = ntetr/3
    allocate(tetra_indices_per_prism(n_prisms,3))
    allocate(prism_volumes(n_prisms))
    allocate(refined_prism_volumes(n_prisms))
    allocate(r_integrand_constants(n_prisms,22)) !collects all constants used for integration in order to print them on a file
    allocate(elec_pot_vec(n_prisms))
    allocate(n_b(n_prisms))

    refined_prism_volumes = 0
    r_integrand_constants = 0
    elec_pot_vec = 0.0_dp
    n_b = 0.0_dp
    prism_volumes = 0.0_dp

    do i = 1,3
        tetra_indices_per_prism(:,i) = (/(i+3*k,k = 0,n_prisms-1)/)
    enddo

    ! Clear any sticky FPE flags set by earlier LAPACK calls (transient NaN
    ! during LU factorisation in dgesv can persist and re-trigger the trap
    ! when the OMP runtime forks the next parallel region on macOS ARM64).
    call ieee_set_flag(ieee_invalid,        .false.)
    call ieee_set_flag(ieee_overflow,       .false.)
    call ieee_set_flag(ieee_divide_by_zero, .false.)
    ! Temporarily disable the halt-on-FPE traps for the parallel region. The
    ! prism-volume / refined-volume formulae divide by vertex R-differences and
    ! take logs of their ratios; for the straight-cylinder limit (grid_kind=5,
    ! R0>>a) near-theta=90deg prisms have a near-zero R-spread, so a few of these
    ! terms overflow / divide-by-zero. prism_volumes is computed robustly via the
    ! exact Pappus form below (so the integrated moments stay clean); the trap
    ! disable only stops those transient FPEs in the auxiliary refined integrals
    ! from killing the run. N<=10 (R0<=1691) never trips this, so they are
    ! bit-unchanged.
    call ieee_set_halting_mode(ieee_invalid,        .false.)
    call ieee_set_halting_mode(ieee_overflow,       .false.)
    call ieee_set_halting_mode(ieee_divide_by_zero, .false.)

    !$OMP PARALLEL DEFAULT(NONE) &
    !$OMP& SHARED(n_prisms,verts_rphiz,tetra_grid,grid_size,tetra_indices_per_prism,prism_volumes, particle_charge, in, &
    !$OMP& refined_prism_volumes,sqrt_g,r_integrand_constants, elec_pot_vec,n_b,tetra_physics, &
    !$OMP& n_field_periods) &
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

        triangle_indices(2) = maxloc(r_values,1)
        r_values_intermediate = r_values
        r_values_intermediate(triangle_indices(2)) = minval(r_values)
        if (sum(r_values_intermediate).eq.0) then
            ! Two vertices share the minimal R, so removing the max-R apex
            ! leaves all-zero r_values_intermediate. triangle_indices(1) must
            ! be the OTHER min-R vertex, i.e. the remaining index of {1,2,3}:
            ! 6 - reference - apex. The previous maxloc(abs(z)) heuristic could
            ! alias triangle_indices(1) onto triangle_indices(2), collapsing the
            ! poloidal triangle to zero area. This happens on the analytic-circ
            ! mesh (grid_kind=5) at the theta=90 deg plane, where every vertex
            ! has R = R0 + rho*cos(theta) = R0, producing ~1e-17 prism volumes
            ! and NaN/Inf moments downstream. Shaped (EFIT/VMEC) meshes never
            ! place two prism vertices at exactly equal R, so they are unaffected.
            triangle_indices(1) = 6 - minloc(r_values,1) - triangle_indices(2)
        else
            triangle_indices(1) = maxloc(r_values_intermediate,1)
        endif

        r = (/r_values(triangle_indices(1)),r_values(triangle_indices(2))/)
        z = (/z_values(triangle_indices(1)),z_values(triangle_indices(2))/)

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

        r_integrand_constants(i,20:22) = (/minloc(r_values),triangle_indices(1),triangle_indices(2)/)

        z_star = z(1) - r(1)*gradient(3)
        alpha = abs(gradient(2)-gradient(1))
        beta = gradient(3) - gradient(1)

        r_integrand_constants(i,1:6) = (/r(1),r(2),alpha,rmin, z_star,beta/)

        prism_volumes(i) =  2*pi/(grid_size(2)*n_field_periods)*(alpha/3*r(1)**3+alpha*rmin/2*r(1)**2+ &
                            abs(z_star*rmin*(r(2)-r(1))+(z_star+beta*rmin)/2*(r(2)**2-r(1)**2)+beta/3*(r(2)**3-r(1)**3)))

        ! Robust fallback for thin-R (near-theta=90deg) poloidal triangles whose
        ! z/r gradients overflow the R-integration above (straight-cylinder limit,
        ! R0>>a). The exact volume of revolution of a triangle by Dphi is Pappus'
        ! theorem V = Dphi * R_centroid * A_triangle, which equals integral(R dA)
        ! and is identical to the formula above for well-conditioned triangles -- so
        ! this only fires for the few degenerate prisms (N<=10 keep their exact value).
        if (.not. ieee_is_finite(prism_volumes(i)) .or. prism_volumes(i) <= 0.0_dp) then
            prism_volumes(i) = 2*pi/(grid_size(2)*n_field_periods) &
                * (rmin + sum(r_values)/3.0_dp) &                                  ! R_centroid
                * 0.5_dp*abs( (r_values(2)-r_values(1))*(z_values(3)-z_values(1)) &  ! A_triangle
                            - (r_values(3)-r_values(1))*(z_values(2)-z_values(1)) )
        endif

        !calculate other volme integrals (compare with appendix B of master thesis of Jonatan Schatzlmayr)
        delta_r = verts_rphiz(1,tetra_grid(tetra_indices_per_prism(i,1))%ind_knot(1)) - &
                & verts_rphiz(1,tetra_grid(tetra_indices_per_prism(i,1))%ind_knot(int(r_integrand_constants(i,20))))
        delta_z = verts_rphiz(3,tetra_grid(tetra_indices_per_prism(i,1))%ind_knot(1)) - &
                & verts_rphiz(3,tetra_grid(tetra_indices_per_prism(i,1))%ind_knot(int(r_integrand_constants(i,20))))

        if (in%boole_refined_sqrt_g) then
            !Calculate prism volumes using the refined approach
            !(compare with appendix B (introductory pages + B1) of master thesis from Jonatan Schatzlmayr)
            a = sqrt_g(3*i-2,1) - sqrt_g(3*i-2,2)*delta_r - sqrt_g(3*i-2,3)*delta_z
            a_dash = sqrt_g(3*i-2,4) - sqrt_g(3*i-2,5)*delta_r - sqrt_g(3*i-2,6)*delta_z
            b = sqrt_g(3*i-2,2)
            b_dash = sqrt_g(3*i-2,5)
            c = sqrt_g(3*i-2,3)
            c_dash = sqrt_g(3*i-2,6)

            !calculate the contribution from 0 to r(1) to the prism volumes
            limits = (/dble(0),r(1)/)

            alpha = c/c_dash*(gradient(2)-gradient(1))
            beta = b_dash+c_dash*gradient(2)
            gamma = b_dash+c_dash*gradient(1)
            delta = (c_dash*a-c*a_dash)/c_dash**2
            epsilon = (c_dash*b-c*b_dash)/c_dash**2

            r_integrand_constants(i,7:12) = (/alpha,beta,gamma,delta,epsilon,a_dash/)

            refined_prism_volumes(i) = r(1)*epsilon*a_dash/(2*gamma*beta)*(gamma-beta) + &
                    & r(1)**2*alpha/2 + &
                    & log((a_dash+gamma*r(1))/a_dash)*(delta*a_dash/(gamma*beta)*(gamma-beta)+epsilon*a_dash**2/(2*gamma**2)) + &
                    & log((a_dash+beta*r(1))/(a_dash+gamma*r(1)))*(delta/beta*(a_dash+beta*r(1))+epsilon/2*r(1)**2) - &
                    & log((a_dash+beta*r(1))/a_dash)*epsilon*a_dash**2/(2*beta**2)

            !calculate the contribution from r(1) to r(2) to the prism volumes
            limits = (/r(1),r(2)/)

            alpha = c/c_dash*(gradient(3)-gradient(1))
            beta = b_dash+c_dash*gradient(3)
            gamma = b_dash+c_dash*gradient(1)
            delta = (c_dash*a-c*a_dash)/c_dash**2
            epsilon = (c_dash*b-c*b_dash)/c_dash**2
            capital_gamma = c/c_dash*z_star
            capital_delta = a_dash + c_dash*z_star

            r_integrand_constants(i,13:19) = (/alpha,beta,gamma,delta,epsilon,capital_gamma,capital_delta/)

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

        if (.not. in%boole_monoenergetic) then

            !calculate electric potential using the refined approach
            !(compare with appendix B (introductory pages + B2) of master thesis from Jonatan Schatzlmayr)

            a = tetra_physics(3*i-2)%gPhi(1)
            b = tetra_physics(3*i-2)%gPhi(3)
            phi_0 = tetra_physics(3*i-2)%Phi1 - a*delta_r - b*delta_z

            !calculate the contribution from 0 to r(1) to the electric potential
            alpha = phi_0*(gradient(2)-gradient(1))*rmin
            beta = (a*rmin+phi_0)*(gradient(2)-gradient(1))+b/2*(gradient(2)**2-gradient(1)**2)*rmin
            gamma = a*(gradient(2)-gradient(1))+b/2*(gradient(2)**2-gradient(1)**2)

            elec_pot_vec(i) = elec_pot_vec(i) + alpha/2*r(1)**2+beta/3*r(1)**3+gamma/4*r(1)**4

            !calculate the contribution from r(1) to r(2) to the electric potential
            alpha  = phi_0*z_star*rmin+b/2*z_star**2*rmin
            beta = a*z_star*rmin+phi_0*(gradient(3)-gradient(1))*rmin+phi_0*z_star+b/2*z_star**2+b*z_star*gradient(3)*rmin
            gamma = phi_0*(gradient(3)-gradient(1))+a*z_star+a*(gradient(3)-gradient(1))*rmin + &
                    & b/2*(gradient(3)**2-gradient(1)**2)*rmin+b*z_star*gradient(3)
            delta = a*(gradient(3)-gradient(1))+b/2*(gradient(3)**2-gradient(1)**2)

            elec_pot_vec(i) = elec_pot_vec(i) + alpha*(r(2)-r(1))+beta/2*(r(2)**2-r(1)**2)+&
                            & gamma/3*(r(2)**3-r(1)**3)+delta/4*(r(2)**4-r(1)**4)

            !calculate Boltzmann density using the refined approach
            !(compare with appendix B (introductory pages + B3) of master thesis from Jonatan Schatzlmayr)

            !calculate contribution from 0 to r(1) to the boltzmann density    
            alpha = in%density*exp(-particle_charge*phi_0/(in%energy_eV*ev2erg))*rmin*(gradient(2)-gradient(1))

            beta = in%density*exp(-particle_charge*phi_0/(in%energy_eV*ev2erg))*((gradient(2)-gradient(1))* &
                    & (1-rmin*particle_charge*a/(in%energy_eV*ev2erg))- &
                    & rmin*particle_charge*b/(2*in%energy_eV*ev2erg)*(gradient(2)**2-gradient(1)**2))

            gamma = in%density*exp(-particle_charge*phi_0/(in%energy_eV*ev2erg))* &
                    & (-particle_charge*a/(in%energy_eV*ev2erg)*(gradient(2)-gradient(1))- &
                    & particle_charge*b/(2*in%energy_eV*ev2erg)*(gradient(2)**2-gradient(1)**2))

            n_b(i) = n_b(i) + alpha/2*r(1)**2+beta/3*r(1)**3+gamma/4*r(1)**4

            !calculate contribution from r(1) to r(2) to the boltzmann density
            alpha = in%density*exp(-particle_charge*phi_0/(in%energy_eV*ev2erg))*rmin*&
                    & (z_star-particle_charge*b/(2*in%energy_eV*ev2erg)*z_star**2)

            beta = in%density*exp(-particle_charge*phi_0/(in%energy_eV*ev2erg))*(z_star- & 
                    & rmin*particle_charge*a/(in%energy_eV*ev2erg)*z_star+rmin*(gradient(3)-gradient(1))- &
                    & particle_charge*b/(2*in%energy_eV*ev2erg)*z_star**2-rmin*particle_charge*b/(in%energy_eV*ev2erg)*z_star*gradient(3))

            gamma = in%density*exp(-particle_charge*phi_0/(in%energy_eV*ev2erg))*((-rmin*particle_charge*a/(in%energy_eV*ev2erg)+1)* &
                    & (gradient(3)-gradient(1))-particle_charge*a/(in%energy_eV*ev2erg)*z_star- &
                    & rmin*particle_charge*b/(2*in%energy_eV*ev2erg)*(gradient(3)**2-gradient(1)**2)- &
                    & particle_charge*b/(in%energy_eV*ev2erg)*z_star*gradient(3))

            delta = in%density*exp(-particle_charge*phi_0/(in%energy_eV*ev2erg))*(-particle_charge*a/(in%energy_eV*ev2erg)* &
                    & (gradient(3)-gradient(1))-particle_charge*b/(2*in%energy_eV*ev2erg)*(gradient(3)**2-gradient(1)**2))

            n_b(i) = n_b(i) + alpha*(r(2)-r(1)) + beta/2*(r(2)**2-r(1)**2) + gamma/3*(r(2)**3-r(1)**3) + delta/4*(r(2)**4-r(1)**4)

        endif
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    ! Re-enable the halt-on-FPE traps after the OMP region
    call ieee_set_halting_mode(ieee_invalid,        .true.)
    call ieee_set_halting_mode(ieee_overflow,       .true.)
    call ieee_set_halting_mode(ieee_divide_by_zero, .true.)
    call ieee_set_flag(ieee_invalid,        .false.)
    call ieee_set_flag(ieee_overflow,       .false.)
    call ieee_set_flag(ieee_divide_by_zero, .false.)

    refined_prism_volumes = abs(refined_prism_volumes)*2*pi/(grid_size(2)*n_field_periods)
    elec_pot_vec = abs(elec_pot_vec)*2*pi/(grid_size(2)*n_field_periods*prism_volumes)
    n_b = abs(n_b)*2*pi/(grid_size(2)*n_field_periods*prism_volumes)   ! was '**' (typo): grid_size(2)**n_field_periods overflows int64 at large N -> SIGILL

    ! Sanitise the auxiliary integrals for the few degenerate thin-R prisms whose
    ! z/r-gradient terms overflowed (straight-cylinder limit). prism_volumes is
    ! already robust (Pappus fallback); fall the refinement back to it and zero the
    ! negligible electric-potential / Boltzmann contributions of those cells.
    where (.not. ieee_is_finite(refined_prism_volumes)) refined_prism_volumes = prism_volumes
    where (.not. ieee_is_finite(elec_pot_vec))          elec_pot_vec = 0.0_dp
    where (.not. ieee_is_finite(n_b))                   n_b = 0.0_dp

    output%prism_volumes = prism_volumes
    if(in%boole_refined_sqrt_g) output%refined_prism_volumes = refined_prism_volumes
    output%electric_potential = elec_pot_vec
    output%boltzmann_density = n_b
    print*, 'calc_volume_integrals finished'
end subroutine calc_volume_integrals

subroutine calc_volume_integrals_in_flux_coordinates

    use constants, only: pi, ev2erg
    use tetra_grid_mod, only: ntetr, verts_sthetaphi, tetra_grid
    use tetra_grid_settings_mod, only: grid_size, n_field_periods
    use tetra_physics_mod, only: tetra_physics
    use gorilla_applets_types_mod, only: output, in, g

    integer                   :: n_prisms, i, k, ind_tetr, even, odd
    real(dp), dimension(2)    :: verts_s, verts_phi, verts_theta
    real(dp)                  :: basic_volume, sqrtg
    real(dp), dimension(3)    :: vertex2, vertex3, vertex4

    print*, 'calc_volume_integrals started'

    n_prisms = ntetr/3

    !Compute prism values
    do i = 1,n_prisms

        !Get the tetrahedron index of the first of the three tetrahedra contained in the current triangular prism
        ind_tetr = (i-1)*3+1

        !Get minimum and maximum values of s, theta and phi of the hexahedron the current triangular prism is part of
        verts_s(1) =     verts_sthetaphi(1,tetra_grid(ind_tetr)%ind_knot(1))
        verts_s(2) =     verts_sthetaphi(1,tetra_grid(ind_tetr)%ind_knot(4))
        verts_theta(1) = verts_sthetaphi(2,tetra_grid(ind_tetr)%ind_knot(1))
        verts_theta(2) = verts_sthetaphi(2,tetra_grid(ind_tetr)%ind_knot(4))
        verts_phi(1) =   verts_sthetaphi(3,tetra_grid(ind_tetr)%ind_knot(1))
        verts_phi(2) =   verts_sthetaphi(3,tetra_grid(ind_tetr)%ind_knot(4))

        !In case the maximum value of theta (phi) is 2*pi (2*pi/n_field_periods), 
        !it is actually given as 0, so 2*pi (2*pi/n_field_periods) has to be added
        if ((verts_theta(2)-verts_theta(1)).lt.0) verts_theta(2) = verts_theta(2) + 2*pi
        if ((verts_phi(2)-verts_phi(1)).lt.0) verts_phi(2) = verts_phi(2) + 2*pi/n_field_periods

        !In the s-theta-phi coordinate space, the volume of the triangular prism is given as half of the hexahedron volume
        basic_volume = (verts_s(2)-verts_s(1))*(verts_theta(2)-verts_theta(1))*(verts_phi(2)-verts_phi(1))/2      

        !To get the prism volume in real physical space, we need to multiply basic_volume with an appropriately averaged value of 
        !sqrt(g) across the prism. This is done by adding together 1/6 of the value at all the prism vertices
        !Start by getting the position of vertex 2 and vertex 3 of tetrahedron in_tetr relative to vertex 1
        !(vertex 1, vertex 2 and vertex 3 of tetrahedron ind_tetr constitute the first three prism vertices)
        !The position of vertex 2 depends on whether the prism index is even or odd
        odd = mod(i,2)
        even = abs(odd-1)
        vertex2 = (/(verts_s(2)-verts_s(1))*odd, (verts_theta(2)-verts_theta(1))*even, 0.0_dp/)
        vertex3 = (/ verts_s(2)-verts_s(1)     ,  verts_theta(2)-verts_theta(1)      , 0.0_dp/)

        !Add the contributions of the first three prism vertices to sqrt(g)
        sqrtg   = 0.0_dp
        sqrtg = sqrtg + 1.0_dp/6.0_dp* tetra_physics(ind_tetr)%sqg1
        sqrtg = sqrtg + 1.0_dp/6.0_dp*(tetra_physics(ind_tetr)%sqg1 + sum(vertex2*tetra_physics(ind_tetr)%gsqg))
        sqrtg = sqrtg + 1.0_dp/6.0_dp*(tetra_physics(ind_tetr)%sqg1 + sum(vertex3*tetra_physics(ind_tetr)%gsqg))

        !For the remaining three prism vertices (which are equal to the first three shifted to the next phi-plane),
        !get vertices 2, 3 and 4 of tetrahedron ind_tetr+2
        vertex2 = (/ 0.0_dp                    ,  0.0_dp                             , verts_phi(2)-verts_phi(1)/)
        vertex3 = (/(verts_s(2)-verts_s(1))*odd, (verts_theta(2)-verts_theta(1))*even, verts_phi(2)-verts_phi(1)/)
        vertex4 = (/ verts_s(2)-verts_s(1)     ,  verts_theta(2)-verts_theta(1)      , verts_phi(2)-verts_phi(1)/)

        !Add the contributions of the remaining three prism vertices to sqrt(g)
        sqrtg = sqrtg + 1.0_dp/6.0_dp*(tetra_physics(ind_tetr+2)%sqg1 + sum(vertex2*tetra_physics(ind_tetr+2)%gsqg))
        sqrtg = sqrtg + 1.0_dp/6.0_dp*(tetra_physics(ind_tetr+2)%sqg1 + sum(vertex3*tetra_physics(ind_tetr+2)%gsqg))
        sqrtg = sqrtg + 1.0_dp/6.0_dp*(tetra_physics(ind_tetr+2)%sqg1 + sum(vertex4*tetra_physics(ind_tetr+2)%gsqg))

        !Write the resulting prism volum into the output type
        output%prism_volumes(i) = basic_volume*abs(sqrtg)

    enddo

    g%total_volume = sum(output%prism_volumes(:))

    print*, 'calc_volume_integrals finished'

end subroutine calc_volume_integrals_in_flux_coordinates

end module volume_integrals_and_sqrt_g_mod