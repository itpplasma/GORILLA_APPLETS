program test_km_benchmark_settings

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use km_benchmark_settings_mod, only: boole_run_energy_scan, &
        boole_run_trace_scan, boole_write_surface_trace, collision_profile_file, &
        diagnostics_prefix, energy_scan_max_factor, energy_scan_min_factor, &
        energy_scan_points, filename_output, find_nearest_profile_surface, &
        fit_end_fraction, fit_start_fraction, &
        load_km_benchmark_inp, n_background_species, n_particles, &
        n_profile_surfaces, n_surfaces, n_trace_scan_multipliers, &
        profile_density, profile_energy_eV, profile_s, profile_temperature, &
        profile_tracer_energy_eV, surface_s_values, trace_scan_multipliers, &
        trace_time_multiplier

    implicit none

    call write_mock_benchmark_input()
    call load_km_benchmark_inp()
    call verify_loaded_values()
    call verify_profile_matching()
    call write_legacy_benchmark_input()
    call load_km_benchmark_inp()
    call verify_legacy_profile_values()
    print '(A)', 'PASS: test_km_benchmark_settings'

contains

subroutine write_mock_benchmark_input()

    integer :: input_unit, profile_unit

    open(newunit=input_unit, file='km_benchmark.inp', status='replace', action='write')
    write(input_unit, '(A)') '&km_benchmark_nml'
    write(input_unit, '(A)') '  n_particles = 321 ,'
    write(input_unit, '(A)') '  n_surfaces = 2 ,'
    write(input_unit, '(A)') '  n_background_species = 3 ,'
    write(input_unit, '(A)') '  boole_run_energy_scan = .true. ,'
    write(input_unit, '(A)') '  energy_scan_points = 5 ,'
    write(input_unit, '(A)') '  energy_scan_min_factor = 0.5d0 ,'
    write(input_unit, '(A)') '  energy_scan_max_factor = 2.0d0 ,'
    write(input_unit, '(A)') '  boole_run_trace_scan = .true. ,'
    write(input_unit, '(A)') '  n_trace_scan_multipliers = 3 ,'
    write(input_unit, '(A)') '  trace_scan_multipliers = 1.0d0, 2.0d0, 4.0d0 ,'
    write(input_unit, '(A)') '  boole_write_surface_trace = .true. ,'
    write(input_unit, '(A)') '  trace_time_multiplier = 3.5d0 ,'
    write(input_unit, '(A)') '  fit_start_fraction = 0.3d0 ,'
    write(input_unit, '(A)') '  fit_end_fraction = 0.9d0 ,'
    write(input_unit, '(A)') "  diagnostics_prefix = 'km_detail' ,"
    write(input_unit, '(A)') "  filename_output = 'km_out.csv' ,"
    write(input_unit, '(A)') '  surface_s_values = 0.2d0, 0.6d0 ,'
    write(input_unit, '(A)') "  collision_profile_file = 'mock_collision_profile.dat'"
    write(input_unit, '(A)') '/'
    close(input_unit)

    open(newunit=profile_unit, file='mock_collision_profile.dat', status='replace', &
        action='write')
    write(profile_unit, '(A)') '# mock collision profile'
    write(profile_unit, '(A)') '2 3'
    write(profile_unit, '(A)') '# s n_1 n_2 n_3 T_1 T_2 T_3 E_ref_eV'
    write(profile_unit, '(A)') '0.20000000  1.0e13  2.0e12  3.0e11  400.0  450.0  500.0  475.0'
    write(profile_unit, '(A)') '0.60000000  1.2e13  2.2e12  3.2e11  510.0  520.0  530.0  515.0'
    close(profile_unit)

end subroutine write_mock_benchmark_input

subroutine write_legacy_benchmark_input()

    integer :: input_unit, profile_unit

    open(newunit=input_unit, file='km_benchmark.inp', status='replace', action='write')
    write(input_unit, '(A)') '&km_benchmark_nml'
    write(input_unit, '(A)') '  n_particles = 111 ,'
    write(input_unit, '(A)') '  n_surfaces = 1 ,'
    write(input_unit, '(A)') '  n_background_species = 2 ,'
    write(input_unit, '(A)') "  collision_profile_file = 'legacy_collision_profile.dat'"
    write(input_unit, '(A)') '/'
    close(input_unit)

    open(newunit=profile_unit, file='legacy_collision_profile.dat', status='replace', &
        action='write')
    write(profile_unit, '(A)') '# legacy collision profile'
    write(profile_unit, '(A)') '1 2'
    write(profile_unit, '(A)') '# s n_1 n_2 E_ref_eV'
    write(profile_unit, '(A)') '0.35000000  4.0e13  4.2e13  650.0'
    close(profile_unit)

end subroutine write_legacy_benchmark_input

subroutine verify_loaded_values()

    call assert_true(boole_run_energy_scan, 'energy scan flag was not loaded')
    call assert_true(boole_run_trace_scan, 'trace scan flag was not loaded')
    call assert_true(boole_write_surface_trace, 'surface trace flag was not loaded')
    call assert_true(n_particles == 321, 'n_particles mismatch')
    call assert_true(n_surfaces == 2, 'n_surfaces mismatch')
    call assert_true(n_background_species == 3, 'n_background_species mismatch')
    call assert_true(energy_scan_points == 5, 'energy_scan_points mismatch')
    call assert_close(energy_scan_min_factor, 0.5_dp, 1.0d-12, &
        'energy_scan_min_factor mismatch')
    call assert_close(energy_scan_max_factor, 2.0_dp, 1.0d-12, &
        'energy_scan_max_factor mismatch')
    call assert_true(n_trace_scan_multipliers == 3, 'trace scan count mismatch')
    call assert_close(trace_scan_multipliers(1), 1.0_dp, 1.0d-12, 'trace multiplier 1 mismatch')
    call assert_close(trace_scan_multipliers(3), 4.0_dp, 1.0d-12, 'trace multiplier 3 mismatch')
    call assert_close(trace_time_multiplier, 3.5_dp, 1.0d-12, 'trace_time_multiplier mismatch')
    call assert_close(fit_start_fraction, 0.3_dp, 1.0d-12, 'fit_start_fraction mismatch')
    call assert_close(fit_end_fraction, 0.9_dp, 1.0d-12, 'fit_end_fraction mismatch')
    call assert_true(trim(diagnostics_prefix) == 'km_detail', 'diagnostics_prefix mismatch')
    call assert_true(trim(filename_output) == 'km_out.csv', 'filename_output mismatch')
    call assert_true(trim(collision_profile_file) == 'mock_collision_profile.dat', &
        'collision_profile_file mismatch')
    call assert_close(surface_s_values(1), 0.2_dp, 1.0d-12, 'surface_s_values(1) mismatch')
    call assert_close(surface_s_values(2), 0.6_dp, 1.0d-12, 'surface_s_values(2) mismatch')
    call assert_true(n_profile_surfaces == 2, 'n_profile_surfaces mismatch')
    call assert_close(profile_s(2), 0.6_dp, 1.0d-12, 'profile_s mismatch')
    call assert_close(profile_density(1, 2), 2.0d12, 1.0d2, 'profile_density mismatch')
    call assert_close(profile_temperature(1, 3), 500.0_dp, 1.0d-12, 'profile_temperature mismatch')
    call assert_close(profile_tracer_energy_eV(2), 515.0_dp, 1.0d-12, 'profile_tracer_energy mismatch')
    call assert_close(profile_energy_eV(2), 515.0_dp, 1.0d-12, 'profile_energy mismatch')

end subroutine verify_loaded_values

subroutine verify_profile_matching()

    call assert_true(find_nearest_profile_surface(0.21_dp) == 1, &
        'nearest profile match for first surface failed')
    call assert_true(find_nearest_profile_surface(0.58_dp) == 2, &
        'nearest profile match for second surface failed')

end subroutine verify_profile_matching

subroutine verify_legacy_profile_values()

    call assert_true(n_particles == 111, 'legacy n_particles mismatch')
    call assert_true(n_surfaces == 1, 'legacy n_surfaces mismatch')
    call assert_true(n_background_species == 2, 'legacy n_background_species mismatch')
    call assert_true(n_profile_surfaces == 1, 'legacy n_profile_surfaces mismatch')
    call assert_close(profile_s(1), 0.35_dp, 1.0d-12, 'legacy profile_s mismatch')
    call assert_close(profile_density(1, 1), 4.0d13, 1.0d2, 'legacy profile_density(1) mismatch')
    call assert_close(profile_density(1, 2), 4.2d13, 1.0d2, 'legacy profile_density(2) mismatch')
    call assert_close(profile_temperature(1, 1), 650.0_dp, 1.0d-12, &
        'legacy profile_temperature(1) mismatch')
    call assert_close(profile_temperature(1, 2), 650.0_dp, 1.0d-12, &
        'legacy profile_temperature(2) mismatch')
    call assert_close(profile_tracer_energy_eV(1), 650.0_dp, 1.0d-12, &
        'legacy profile_tracer_energy mismatch')

    call cleanup_file('legacy_collision_profile.dat')
    call cleanup_file('km_benchmark.inp')
    call cleanup_file('mock_collision_profile.dat')

end subroutine verify_legacy_profile_values

subroutine cleanup_file(path)

    character(len=*), intent(in) :: path

    integer :: file_unit
    logical :: exists

    inquire(file=path, exist=exists)
    if (.not. exists) return
    open(newunit=file_unit, file=path, status='old')
    close(file_unit, status='delete')

end subroutine cleanup_file

subroutine assert_close(actual, expected, tolerance, message)

    real(dp), intent(in) :: actual
    real(dp), intent(in) :: expected
    real(dp), intent(in) :: tolerance
    character(len=*), intent(in) :: message

    if (abs(actual - expected) > tolerance) then
        print '(A,3(1X,ES16.8))', trim(message), actual, expected, tolerance
        error stop 1
    end if

end subroutine assert_close

subroutine assert_true(condition, message)

    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) then
        print '(A)', trim(message)
        error stop 1
    end if

end subroutine assert_true

end program test_km_benchmark_settings
