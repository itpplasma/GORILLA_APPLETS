program test_km_benchmark_diagnostics

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use km_benchmark_diagnostics_mod, only: build_log_energy_grid, &
        collision_trace_t, compute_collision_trace, compute_energy_convolution, &
        compute_maxwellian_energy_weights, compute_relative_change, &
        compute_tracing_time_seconds, estimate_collision_time_seconds
    use constants, only: ame, amp, echarge
    use gorilla_applets_types_mod, only: start
    use kramers_moyal_transport_mod, only: km_d11_result_t, set_trace_time_for_species, &
        write_km_csv
    use utils_parallelised_particle_pushing_mod, only: is_valid_tetra_index

    implicit none

    call test_collision_trace_is_positive()
    call test_collision_time_scales_with_density()
    call test_log_energy_grid_is_geometric()
    call test_maxwellian_weights_normalize()
    call test_energy_convolution_preserves_constant_profile()
    call test_tracing_time_prefers_override()
    call test_relative_change_is_symmetric()
    call test_tetra_index_validation()
    call test_trace_time_is_written_to_start_state()
    call test_write_km_csv_emits_enriched_columns()

contains

subroutine test_collision_trace_is_positive()

    type(collision_trace_t) :: trace
    real(dp), dimension(2) :: charges, densities, masses, temperatures

    masses = (/2.0_dp * amp, ame/)
    charges = (/1.0_dp, -1.0_dp/)
    densities = (/3.0d13, 3.0d13/)
    temperatures = (/500.0_dp, 500.0_dp/)

    call compute_collision_trace(ame, -echarge, masses, charges, densities, &
        temperatures, 500.0_dp, trace)

    call assert_true(all(trace%lambda > 0.0_dp), 'collision trace lambda must be positive')
    call assert_true(all(trace%efcolf > 0.0_dp), 'collision trace efcolf must be positive')
    call assert_close(sum(trace%collision_frequency_hz), &
        trace%total_collision_frequency_hz, 1.0d-12, 'collision frequency sum mismatch')
    call assert_true(estimate_collision_time_seconds(trace) > 0.0_dp, &
        'collision time must be positive')
    print '(A)', 'PASS: test_collision_trace_is_positive'

end subroutine test_collision_trace_is_positive

subroutine test_collision_time_scales_with_density()

    type(collision_trace_t) :: trace_dense
    type(collision_trace_t) :: trace_sparse
    real(dp), dimension(2) :: charges, masses, temperatures
    real(dp), dimension(2) :: dense_densities, sparse_densities
    real(dp) :: ratio

    masses = (/2.0_dp * amp, ame/)
    charges = (/1.0_dp, -1.0_dp/)
    sparse_densities = (/2.0d13, 2.0d13/)
    dense_densities = 2.0_dp * sparse_densities
    temperatures = (/500.0_dp, 500.0_dp/)

    call compute_collision_trace(ame, -echarge, masses, charges, sparse_densities, &
        temperatures, 500.0_dp, trace_sparse)
    call compute_collision_trace(ame, -echarge, masses, charges, dense_densities, &
        temperatures, 500.0_dp, trace_dense)

    ratio = estimate_collision_time_seconds(trace_sparse) / &
        estimate_collision_time_seconds(trace_dense)
    call assert_close(ratio, 2.0_dp, 5.0d-2, 'collision time should scale inversely with density')
    print '(A)', 'PASS: test_collision_time_scales_with_density'

end subroutine test_collision_time_scales_with_density

subroutine test_log_energy_grid_is_geometric()

    real(dp), allocatable :: energy_grid(:)

    call build_log_energy_grid(500.0_dp, 0.5_dp, 2.0_dp, 5, energy_grid)

    call assert_close(energy_grid(1), 250.0_dp, 1.0d-10, 'log grid minimum mismatch')
    call assert_close(energy_grid(5), 1000.0_dp, 1.0d-10, 'log grid maximum mismatch')
    call assert_close(energy_grid(3)**2, energy_grid(2) * energy_grid(4), 1.0d-8, &
        'log grid should be geometric')
    print '(A)', 'PASS: test_log_energy_grid_is_geometric'

end subroutine test_log_energy_grid_is_geometric

subroutine test_maxwellian_weights_normalize()

    real(dp), dimension(5) :: energy_eV
    real(dp), allocatable :: weights(:)

    energy_eV = (/250.0_dp, 353.5533905932738_dp, 500.0_dp, 707.1067811865476_dp, 1000.0_dp/)
    call compute_maxwellian_energy_weights(energy_eV, 500.0_dp, weights)

    call assert_close(sum(weights), 1.0_dp, 1.0d-12, 'energy weights must normalize')
    call assert_true(all(weights > 0.0_dp), 'energy weights must be positive')
    call assert_true(weights(4) > weights(1), 'D11 energy weights must favor supra-thermal energies')
    print '(A)', 'PASS: test_maxwellian_weights_normalize'

end subroutine test_maxwellian_weights_normalize

subroutine test_energy_convolution_preserves_constant_profile()

    real(dp), dimension(5) :: energy_eV
    real(dp), dimension(5) :: transport_values
    real(dp), allocatable :: weights(:)
    real(dp) :: convolved_value

    energy_eV = (/250.0_dp, 353.5533905932738_dp, 500.0_dp, 707.1067811865476_dp, 1000.0_dp/)
    transport_values = 3.5_dp
    call compute_energy_convolution(energy_eV, transport_values, 500.0_dp, &
        convolved_value, weights)

    call assert_close(convolved_value, 3.5_dp, 1.0d-12, &
        'energy convolution must preserve constant profiles')
    call assert_true(all(weights > 0.0_dp), 'energy convolution weights must be positive')
    print '(A)', 'PASS: test_energy_convolution_preserves_constant_profile'

end subroutine test_energy_convolution_preserves_constant_profile

subroutine test_tracing_time_prefers_override()

    real(dp) :: trace_time

    trace_time = compute_tracing_time_seconds(1.2d-4, 4.0_dp, 3.0d-5)
    call assert_close(trace_time, 1.2d-4, 1.0d-15, 'trace time override must win')

    trace_time = compute_tracing_time_seconds(0.0_dp, 4.0_dp, 3.0d-5)
    call assert_close(trace_time, 1.2d-4, 1.0d-15, 'trace time multiplier mismatch')
    print '(A)', 'PASS: test_tracing_time_prefers_override'

end subroutine test_tracing_time_prefers_override

subroutine test_relative_change_is_symmetric()

    real(dp) :: change_a, change_b

    change_a = compute_relative_change(3.0_dp, 4.5_dp)
    change_b = compute_relative_change(4.5_dp, 3.0_dp)
    call assert_close(change_a, change_b, 1.0d-12, 'relative change must be symmetric')
    print '(A)', 'PASS: test_relative_change_is_symmetric'

end subroutine test_relative_change_is_symmetric

subroutine test_tetra_index_validation()

    call assert_true(is_valid_tetra_index(1, 8), 'first tetra index must be valid')
    call assert_true(is_valid_tetra_index(8, 8), 'last tetra index must be valid')
    call assert_true(.not. is_valid_tetra_index(0, 8), 'zero tetra index must be invalid')
    call assert_true(.not. is_valid_tetra_index(9, 8), 'out-of-range tetra index must be invalid')
    call assert_true(.not. is_valid_tetra_index(-1, 8), 'negative tetra index must be invalid')
    print '(A)', 'PASS: test_tetra_index_validation'

end subroutine test_tetra_index_validation

subroutine test_trace_time_is_written_to_start_state()

    allocate(start%t(2))
    start%t = 0.0_dp

    call set_trace_time_for_species(3.25d-4, 2)
    call assert_close(start%t(2), 3.25d-4, 1.0d-15, 'trace time handoff mismatch')

    deallocate(start%t)
    print '(A)', 'PASS: test_trace_time_is_written_to_start_state'

end subroutine test_trace_time_is_written_to_start_state

subroutine test_write_km_csv_emits_enriched_columns()

    type(km_d11_result_t) :: result
    character(len=256) :: header
    character(len=256) :: row
    integer :: file_unit

    allocate(result%boundary_s(3))
    allocate(result%surface_indices(1))
    allocate(result%d11(1))
    allocate(result%convection(1))
    allocate(result%energy_eV(1))
    allocate(result%collision_time_s(1))
    allocate(result%trace_time_s(1))
    allocate(result%total_collision_frequency_hz(1))
    allocate(result%convolved_d11(1))

    result%n_surfaces = 1
    result%boundary_s = (/0.1_dp, 0.4_dp, 0.8_dp/)
    result%surface_indices = (/2/)
    result%d11 = (/1.25_dp/)
    result%convection = (/-0.5_dp/)
    result%energy_eV = (/750.0_dp/)
    result%collision_time_s = (/2.0d-4/)
    result%trace_time_s = (/8.0d-4/)
    result%total_collision_frequency_hz = (/5.0d3/)
    result%convolved_d11 = (/1.75_dp/)

    call write_km_csv(result, 'km_test_profile.csv')

    open(newunit=file_unit, file='km_test_profile.csv', status='old', action='read')
    read(file_unit, '(A)') header
    read(file_unit, '(A)') row
    close(file_unit)

    call assert_true(trim(header) == &
        's,d11_s,convection_A,energy_eV,collision_time_s,trace_time_s,total_collision_frequency_hz,d11_convolved_s', &
        'km csv header mismatch')
    call assert_true(index(row, '4.0000000000000002E-01') > 0, 'km csv must contain surface s')
    call assert_true(index(row, '1.7500000000000000E+00') > 0, 'km csv must contain convolved D11')

    open(newunit=file_unit, file='km_test_profile.csv', status='old')
    close(file_unit, status='delete')
    print '(A)', 'PASS: test_write_km_csv_emits_enriched_columns'

end subroutine test_write_km_csv_emits_enriched_columns

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

end program test_km_benchmark_diagnostics
