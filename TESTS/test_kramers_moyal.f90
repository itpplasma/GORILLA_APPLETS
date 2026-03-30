program test_kramers_moyal
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use transport_benchmark_utils_mod, only: fit_transport_coefficients
    implicit none

    call test_fit_transport_coefficients_pure_diffusion()
    call test_fit_transport_coefficients_with_drift()
    call test_kramers_moyal_1d_random_walk()

contains

subroutine test_fit_transport_coefficients_pure_diffusion()
    real(dp), parameter :: D_true = 0.03_dp
    integer, parameter :: n_samples = 1000
    real(dp) :: time(n_samples), delta_s(n_samples), delta_s_squared(n_samples)
    real(dp) :: A_fit, B_fit
    integer :: i

    do i = 1, n_samples
        time(i) = real(i, dp) * 1.0d-4
        delta_s(i) = 0.0_dp
        delta_s_squared(i) = 2.0_dp * D_true * time(i)
    end do

    call fit_transport_coefficients(time, delta_s, delta_s_squared, A_fit, B_fit)

    call assert_close(A_fit, 0.0_dp, 1.0d-10, 'Pure diffusion: A must be zero')
    call assert_close(B_fit, D_true, 1.0d-10, 'Pure diffusion: B must equal D_true')
    print '(A)', 'PASS: test_fit_transport_coefficients_pure_diffusion'

end subroutine

subroutine test_fit_transport_coefficients_with_drift()
    real(dp), parameter :: D_true = 0.02_dp
    real(dp), parameter :: A_true = 0.5_dp
    integer, parameter :: n_samples = 1000
    real(dp) :: time(n_samples), delta_s(n_samples), delta_s_squared(n_samples)
    real(dp) :: A_fit, B_fit
    integer :: i

    do i = 1, n_samples
        time(i) = real(i, dp) * 1.0d-4
        delta_s(i) = A_true * time(i)
        delta_s_squared(i) = 2.0_dp * D_true * time(i) + A_true**2 * time(i)**2
    end do

    call fit_transport_coefficients(time, delta_s, delta_s_squared, A_fit, B_fit)

    call assert_close(A_fit, A_true, 1.0d-6, 'With drift: A recovery failed')
    call assert_close(B_fit, D_true, 1.0d-6, 'With drift: B recovery failed')
    print '(A)', 'PASS: test_fit_transport_coefficients_with_drift'

end subroutine

subroutine test_kramers_moyal_1d_random_walk()
    integer, parameter :: n_cells = 20
    integer, parameter :: n_particles = 50000
    integer, parameter :: n_steps = 1000
    real(dp), parameter :: dx = 1.0_dp / n_cells
    real(dp), parameter :: dt = 1.0d-4
    real(dp), parameter :: pi = 3.141592653589793d0

    real(dp) :: D_true(n_cells)
    real(dp) :: km_D(n_cells), km_A(n_cells), km_count(n_cells)
    real(dp) :: positions(n_particles)
    real(dp) :: x, new_x, delta_x, D_local, u1, u2, xi
    real(dp) :: max_err
    integer :: i, step, cell, n_valid
    logical :: alive(n_particles)

    ! True D(x) = 0.02 + 0.02 * sin(pi*x)
    do i = 1, n_cells
        D_true(i) = 0.02_dp + 0.02_dp * sin(pi * (real(i, dp) - 0.5_dp) * dx)
    end do

    km_D = 0.0_dp
    km_A = 0.0_dp
    km_count = 0.0_dp
    alive = .true.

    ! Initialize positions uniformly in [0.1, 0.9]
    do i = 1, n_particles
        call random_number(x)
        positions(i) = 0.1_dp + 0.8_dp * x
    end do

    ! Random walk with KM accumulation
    do step = 1, n_steps
        do i = 1, n_particles
            if (.not. alive(i)) cycle
            x = positions(i)
            cell = min(max(int(x / dx) + 1, 1), n_cells)
            D_local = D_true(cell)

            call random_number(u1)
            call random_number(u2)
            xi = sqrt(-2.0_dp * log(max(u1, 1.0d-30))) * cos(6.283185307179586d0 * u2)
            delta_x = sqrt(2.0_dp * D_local * dt) * xi

            ! KM accumulation at pre-step position
            km_A(cell) = km_A(cell) + delta_x
            km_D(cell) = km_D(cell) + delta_x**2
            km_count(cell) = km_count(cell) + 1.0_dp

            new_x = x + delta_x
            if (new_x < 0.0_dp) new_x = -new_x
            if (new_x >= 1.0_dp) then
                alive(i) = .false.
                cycle
            end if
            positions(i) = new_x
        end do
    end do

    ! Compute KM coefficients
    do i = 1, n_cells
        if (km_count(i) > 100.0_dp) then
            km_A(i) = km_A(i) / (km_count(i) * dt)
            km_D(i) = km_D(i) / (2.0_dp * km_count(i) * dt)
        else
            km_A(i) = -1.0_dp
            km_D(i) = -1.0_dp
        end if
    end do

    ! Check: D recovery within 10%
    max_err = 0.0_dp
    n_valid = 0
    do i = 2, n_cells - 1
        if (km_count(i) > 1000.0_dp) then
            max_err = max(max_err, abs(km_D(i) - D_true(i)) / D_true(i))
            n_valid = n_valid + 1
        end if
    end do

    if (max_err > 0.10_dp) then
        print '(A,ES12.4)', 'FAIL: KM 1D random walk max relative error = ', max_err
        stop 1
    end if
    if (n_valid < n_cells / 2) then
        print '(A,I0)', 'FAIL: KM too few valid cells = ', n_valid
        stop 1
    end if

    ! Check: A near zero (no drift in pure diffusion)
    max_err = 0.0_dp
    do i = 2, n_cells - 1
        if (km_count(i) > 1000.0_dp) then
            max_err = max(max_err, abs(km_A(i)))
        end if
    end do

    print '(A,ES10.2,A,I0,A)', 'PASS: test_kramers_moyal_1d_random_walk (D_err=', &
        max_err, ', valid=', n_valid, ' cells)'

end subroutine

subroutine assert_close(actual, expected, tolerance, message)
    real(dp), intent(in) :: actual, expected, tolerance
    character(len=*), intent(in) :: message
    if (abs(actual - expected) > tolerance) then
        print '(A,3(1X,ES16.8))', trim(message), actual, expected, tolerance
        stop 1
    end if
end subroutine

end program test_kramers_moyal
