module global_transport_fit_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

contains

subroutine calc_global_transport_fit()

    use global_transport_fit_core_mod, only: fit_global_transport
    use global_transport_fit_gorilla_mod, only: run_gorilla_mc_global_transport_fit
    use global_transport_fit_io_mod, only: load_global_transport_experiments, write_convergence_history, &
        write_global_transport_fit_outputs
    use global_transport_fit_settings_mod, only: control, filename_boundary_s, filename_convergence_history, &
        filename_density_1, filename_density_2, filename_density_variance_1, filename_density_variance_2, &
        filename_fit_profiles, filename_fit_summary, filename_flux_1, filename_flux_2, filename_flux_variance_1, &
        filename_flux_variance_2, filename_shell_volumes, filename_source_1, filename_source_2, &
        load_global_transport_fit_inp, run_mode
    use global_transport_fit_types_mod, only: global_transport_experiment_t, global_transport_fit_result_t

    type(global_transport_experiment_t), allocatable :: experiments(:)
    type(global_transport_fit_result_t) :: result
    character(len=256), dimension(2) :: density_files
    character(len=256), dimension(2) :: density_var_files
    character(len=256), dimension(2) :: flux_files
    character(len=256), dimension(2) :: flux_var_files
    character(len=256), dimension(2) :: source_files

    call load_global_transport_fit_inp()
    if (trim(adjustl(run_mode)) == 'gorilla_mc') then
        call run_gorilla_mc_global_transport_fit(control, result)
        print *, 'Global transport fit complete.'
        print *, 'objective = ', result%objective
        print *, 'converged = ', result%converged
        print *, 'iterations = ', result%n_iterations
        return
    end if

    source_files = [filename_source_1, filename_source_2]
    density_files = [filename_density_1, filename_density_2]
    density_var_files = [filename_density_variance_1, filename_density_variance_2]
    flux_files = [filename_flux_1, filename_flux_2]
    flux_var_files = [filename_flux_variance_1, filename_flux_variance_2]

    call load_global_transport_experiments(filename_boundary_s, filename_shell_volumes, source_files, density_files, &
        density_var_files, flux_files, flux_var_files, experiments)
    call fit_global_transport(experiments, control, result)
    call write_global_transport_fit_outputs(experiments(1)%boundary_s, result, filename_fit_summary, filename_fit_profiles)
    call write_convergence_history(result, filename_convergence_history)

    print *, 'Global transport fit complete.'
    print *, 'objective = ', result%objective
    print *, 'converged = ', result%converged
    print *, 'iterations = ', result%n_iterations

end subroutine calc_global_transport_fit

end module global_transport_fit_mod
