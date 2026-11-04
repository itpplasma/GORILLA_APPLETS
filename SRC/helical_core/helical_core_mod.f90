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

    use utils_data_pre_and_post_processing_mod, only: setup_applet, print_applet_summary
    use tetra_grid_mod, only: verts_rphiz
    use volume_integrals_and_sqrt_g_mod, only: calc_square_root_g, calc_volume_integrals
    use gorilla_applets_types_mod, only: moment_specs, in, start, output, weights
    use utils_write_data_to_files_mod, only: write_data_to_files, give_file_names, unlink_files
    use utils_data_pre_and_post_processing_mod, only: set_seed_for_random_numbers, initialize_exit_data, &
    calc_poloidal_flux, calc_collision_coefficients_for_all_tetrahedra, &
    normalise_prism_moments_and_prism_moments_squared, fourier_transform_moments, calc_starting_conditions
    use utils_helical_core_mod, only: read_helical_core_inp_into_type, &
        parallelised_particle_pushing_helical_core, eliminate_particles_outside_flux_threshold

    call set_seed_for_random_numbers
    call read_helical_core_inp_into_type
    call setup_applet

    call calc_square_root_g
    call calc_volume_integrals

    call initialize_exit_data
    call calc_poloidal_flux(verts_rphiz)
    call give_file_names
    call unlink_files

    if (in%boole_collisions) call calc_collision_coefficients_for_all_tetrahedra

    call calc_starting_conditions
    call eliminate_particles_outside_flux_threshold

    ! Standard particle tracing (no scan, no diffusion coefficient calculation)
    call parallelised_particle_pushing_helical_core(species=1)

    call normalise_prism_moments_and_prism_moments_squared(boole_skip_time_normalisation_in=in%boole_delta_f)

    if (moment_specs%n_moments.gt.0) call fourier_transform_moments
    call write_data_to_files

    call print_applet_summary('ion', start%particle_mass(1), start%v0(1), start%particle_charge(1))

end subroutine calc_helical_core

end module helical_core_mod
