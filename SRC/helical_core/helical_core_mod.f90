module helical_core_mod
!
! Module for helical core particle tracing.
!
! This is a simplified version of anomalous_transport_mod that:
!   - Traces particles in a helical core geometry
!   - Does NOT include anomalous transport displacement (straight-line displacements)
!   - Does NOT include electric potential perturbation scans
!
! The module uses the same infrastructure as anomalous transport but removes
! all anomalous-specific features.
!

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

contains

subroutine calc_helical_core

    use orbit_timestep_gorilla_mod, only: initialize_gorilla
    use constants, only: ev2erg
    use tetra_grid_settings_mod, only: grid_kind
    use tetra_grid_mod, only: verts_rphiz
    use gorilla_applets_settings_mod, only: i_option
    use field_mod, only: ipert
    use volume_integrals_and_sqrt_g_mod, only: calc_square_root_g, calc_volume_integrals
    use gorilla_applets_types_mod, only: moment_specs, counter, c, in, start, output, weights
    use utils_write_data_to_files_mod, only: write_data_to_files, give_file_names, unlink_files
    use utils_data_pre_and_post_processing_mod, only: set_seed_for_random_numbers, &
        get_ipert, set_moment_specifications, initialise_output, initialize_exit_data, calc_poloidal_flux, &
        calc_collision_coefficients_for_all_tetrahedra, normalise_prism_moments_and_prism_moments_squared, &
        fourier_transform_moments, calc_starting_conditions
    use utils_helical_core_mod, only: read_helical_core_inp_into_type, &
        parallelised_particle_pushing_helical_core, eliminate_particles_outside_flux_threshold, &
        set_unidirectional_beam, write_collision_relaxation_diagnostics, write_collision_relaxation_profiles

    integer :: i_slice, relaxation_unit
    logical :: boole_relaxation_diagnostics
    real(dp), dimension(:), allocatable :: pitch_0, energy_0
    real(dp), dimension(:,:), allocatable :: x_0
    logical, dimension(:), allocatable :: lost_0
    real(dp), dimension(:), allocatable :: u_b_0, temp_b_0

    call set_seed_for_random_numbers
    call read_helical_core_inp_into_type
    call get_ipert()
    call initialize_gorilla(i_option,ipert)

    call set_moment_specifications
    call initialise_output
    call calc_square_root_g
    call calc_volume_integrals

    call initialize_exit_data
    call calc_poloidal_flux(verts_rphiz)
    call give_file_names
    call unlink_files

    if (in%boole_collisions) call calc_collision_coefficients_for_all_tetrahedra

    call calc_starting_conditions
    call eliminate_particles_outside_flux_threshold
    call set_unidirectional_beam

    ! Standard particle tracing (no scan, no diffusion coefficient calculation). With n_time_slices > 1,
    ! all markers are pushed slice by slice, so markers and the background reservoir of the conserving
    ! collisions evolve in common time.
    boole_relaxation_diagnostics = in%boole_collisions.and.in%boole_preserve_energy_and_momentum_during_collisions
    if (boole_relaxation_diagnostics) then
        x_0 = start%x(:,:,1)
        pitch_0 = start%pitch(:,1)
        energy_0 = start%energy(:,1)
        lost_0 = start%lost(:,1)
        u_b_0 = c%vpar_mat(1,:)
        temp_b_0 = c%temp_mat(1,:)
        open(newunit=relaxation_unit, file='collision_relaxation.dat', status='replace')
        write(relaxation_unit,'(A)') '#time <v_par>_m/v0 <u_par>_b/v0 P_m P_b P_lost E_m E_b E_lost E_tot ' // &
                                     'L_m L_b L_lost L_tot <T>_b T_m min_T_b max|u_b|/v0 n_confined'
    endif
    start%t = in%time_step/in%n_time_slices
    do i_slice = 1, in%n_time_slices
        call parallelised_particle_pushing_helical_core(species=1, boole_reuse_weights_in=(i_slice.gt.1), i_slice_in=i_slice)
        if (boole_relaxation_diagnostics) then
            !weights are known after the first slice and stay fixed, so the initial row can be written now
            if (i_slice.eq.1) call write_collision_relaxation_diagnostics(relaxation_unit, 0.0_dp, x_0, pitch_0, energy_0, &
                                                                          lost_0, u_b_0, temp_b_0)
            call write_collision_relaxation_diagnostics(relaxation_unit, i_slice*start%t(1), start%x(:,:,1), start%pitch(:,1), &
                                                        start%energy(:,1), start%lost(:,1), c%vpar_mat(1,:), c%temp_mat(1,:))
        endif
    enddo
    start%t = in%time_step
    if (boole_relaxation_diagnostics) then
        close(relaxation_unit)
        call write_collision_relaxation_profiles
    endif

    call normalise_prism_moments_and_prism_moments_squared(boole_skip_time_normalisation_in=in%boole_delta_f)

    if (moment_specs%n_moments.gt.0) call fourier_transform_moments
    call write_data_to_files

    if (in%boole_precalc_collisions) print*, "maxcol = ", c%maxcol
    print*, 'Number of lost ions',counter%lost_particles
    print*, 'average number of pushings = ', counter%tetr_pushings/in%n_particles
    print*, 'average number of toroidal revolutions = ', counter%phi_0_mappings/in%n_particles
    print*, 'average number of integration steps = ', counter%integration_steps/in%n_particles
    PRINT*, 'ion mass = ', start%particle_mass(1)
    PRINT*, 'absolute value of ion velocity = ', start%v0(1)
    PRINT*, 'ion charge = ', start%particle_charge(1)
    PRINT*, 'temperature = ', ev2erg*in%energy_eV
    print*, 'energy in eV = ', in%energy_eV
    print*, 'tracing time in seconds = ', in%time_step
    if((grid_kind.eq.2).or.(grid_kind.eq.3)) then
         print*, 'number of times that particles were pushed across the inside hole = ', counter%lost_inside
    endif

end subroutine calc_helical_core

end module helical_core_mod
