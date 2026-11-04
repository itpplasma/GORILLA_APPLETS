module anomalous_transport_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

contains

subroutine calc_anomalous_transport

    use utils_data_pre_and_post_processing_mod, only: setup_applet, print_applet_summary
    use tetra_grid_mod, only: verts_rphiz
    use volume_integrals_and_sqrt_g_mod, only: calc_square_root_g, calc_volume_integrals
    use gorilla_applets_types_mod, only: moment_specs, in, start, output
    use utils_write_data_to_files_mod, only: write_data_to_files, give_file_names, unlink_files
    use utils_data_pre_and_post_processing_mod, only: set_seed_for_random_numbers, initialize_exit_data, &
    calc_poloidal_flux, calc_collision_coefficients_for_all_tetrahedra, &
    normalise_prism_moments_and_prism_moments_squared, fourier_transform_moments, calc_starting_conditions
    use utils_anomalous_transport_mod, only: read_anomalous_transport_inp_into_type, &
        parallelised_particle_pushing_anomalous_transport, calc_diffusion_coefficient, &
        scan_anomalous_transport

    call set_seed_for_random_numbers
    call read_anomalous_transport_inp_into_type
    call setup_applet

    call calc_square_root_g
    call calc_volume_integrals

    call initialize_exit_data
    call calc_poloidal_flux(verts_rphiz)
    call give_file_names
    call unlink_files

    if (in%boole_collisions) call calc_collision_coefficients_for_all_tetrahedra

    call calc_starting_conditions

    ! Scan over parameter, calculate diffusion coefficient, or run standard transport
    ! i_scan_option: 0=no scan, 1=scan over eps_Phi, 2=scan over n3, 3=scan over n2
    if (in%i_scan_option > 0) then
        call scan_anomalous_transport
    elseif (in%boole_calc_diffusion_coefficient) then
        call calc_diffusion_coefficient
    else
        call parallelised_particle_pushing_anomalous_transport(species=1)
    endif

    call normalise_prism_moments_and_prism_moments_squared


    if (moment_specs%n_moments.gt.0) call fourier_transform_moments
    call write_data_to_files

    call print_applet_summary('ion', start%particle_mass(1), start%v0(1), start%particle_charge(1))

end subroutine calc_anomalous_transport

end module anomalous_transport_mod
