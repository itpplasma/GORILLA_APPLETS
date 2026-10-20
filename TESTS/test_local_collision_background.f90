! Exercise the production centroid/background update, including its collision rates.
program test_local_collision_background
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use gorilla_settings_mod, only: load_gorilla_inp
    use tetra_grid_settings_mod, only: load_tetra_grid_inp
    use orbit_timestep_gorilla_mod, only: initialize_gorilla
    use tetra_grid_mod, only: ntetr, tetra_grid, verts_sthetaphi
    use tetra_physics_mod, only: particle_mass, particle_charge
    use gorilla_applets_types_mod, only: in, c
    use utils_data_pre_and_post_processing_mod, only: set_c
    use utils_rmp_response_currents_mod, only: set_local_collision_background
    use profile_data_mod, only: load_profiles, eval_profiles, profile_values_t
    use collis_ions, only: collis_init, coleff
    use constants, only: echarge
    implicit none
    character(len=64) :: arg
    integer :: i, j, a, u
    real(dp) :: s, v0, ef(2), vr(2), er(2), dens(2), temp(2)
    real(dp), allocatable :: dpp(:), dhh(:), drag(:)
    type(profile_values_t) :: pv

    call get_command_argument(1, arg)
    read(arg, *) in%n_background_species
    in%energy_eV = 3500.0_dp
    call load_tetra_grid_inp()
    call load_gorilla_inp()
    call initialize_gorilla(ipert_in=0)
    call load_profiles('profiles', 'equil.dat')
    call set_c()
    allocate(dpp(c%n), dhh(c%n), drag(c%n))
    open(newunit=u, file='rates.dat', status='replace')
    do a = 1, c%n
        ! First obtain the source identities, then select each as the test particle.
        call set_local_collision_background()
        particle_mass = c%mass(a)
        if (c%n == 2 .and. a == 1) particle_mass = 3.3436e-24_dp
        particle_charge = c%charge_num(a)*echarge
        call set_local_collision_background()
        do i = 1, ntetr
            s = sum(verts_sthetaphi(1, tetra_grid(i)%ind_knot))/4.0_dp
            call eval_profiles(s, pv)
            if (c%n == 2) then
                ! Frozen pre-fix algorithm, including its evaluation order.
                dens = pv%n_e
                temp = [pv%Ti, pv%Te]
                call collis_init(particle_mass, c%charge_num(a), c%mass, &
                    c%charge_num, dens, temp, in%energy_eV, v0, ef, vr, er)
                if (any(transfer(ef, [0_int64], 2) /= &
                    transfer(c%efcolf_mat(:,i), [0_int64], 2))) &
                    error stop 'Pure D/e collision frequencies changed bits'
                if (any(transfer(vr, [0_int64], 2) /= &
                    transfer(c%velrat_mat(:,i), [0_int64], 2))) &
                    error stop 'Pure D/e velocity ratios changed bits'
                if (any(transfer(er, [0_int64], 2) /= &
                    transfer(c%enrat_mat(:,i), [0_int64], 2))) &
                    error stop 'Pure D/e energy ratios changed bits'
            end if
            v0 = sqrt(2.0_dp*in%energy_eV*1.6022e-12_dp/particle_mass)
            call coleff(c%efcolf_mat(:,i), c%velrat_mat(:,i), &
                c%enrat_mat(:,i), 1.0_dp, dpp, dhh, drag)
            do j = 1, c%n
                write(u, '(3i8,9es25.16)') a, j, i, s, particle_mass, &
                    c%mass(j), c%charge_num(a), c%charge_num(j), &
                    c%dens_mat(j,i), c%temp_mat(j,i), 2.0_dp*v0*dhh(j), &
                    c%efcolf_mat(j,i)
            end do
        end do
    end do
    close(u)
    print *, 'PASS: production local collision background'
end program test_local_collision_background
