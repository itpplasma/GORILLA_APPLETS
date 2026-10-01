program test_rmp_gc_measure
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use rmp_gc_measure_mod, only: gc_measure_correction, gc_log_measure_change
    implicit none
    real(dp), parameter :: beta = 0.08_dp, shift = 0.2_dp, pi = acos(-1.0_dp)
    real(dp) :: h(3), bs(3), u, proposed, chi, next_chi, rnd(3), log_accept
    real(dp) :: mass, momentum, energy, dx, weight, mean_u, var_u, sum_u, sum_u2
    integer :: i, nseed
    integer, allocatable :: seed(:)

    h = [0.0_dp, 2.0_dp, 0.0_dp]
    dx = 16.0_dp/20000
    mass = 0.0_dp
    momentum = 0.0_dp
    energy = 0.0_dp
    do i = 1, 20000
        u = -8.0_dp + (i-0.5_dp)*dx
        bs = [0.0_dp, 3.0_dp*(1.0_dp+beta*u), 0.0_dp]
        chi = gc_measure_correction(h, bs, 2.0_dp, 3.0_dp)
        weight = chi*exp(-0.5_dp*u*u)/sqrt(2.0_dp*pi)*dx
        mass = mass + weight
        momentum = momentum + u*weight
        energy = energy + u*u*weight
    end do
    if (abs(mass-1.0_dp) > 1.0e-12_dp) error stop 1
    if (abs(momentum-beta) > 1.0e-12_dp) error stop 2
    if (abs(energy-1.0_dp) > 1.0e-12_dp) error stop 3
    call random_seed(size=nseed)
    allocate(seed(nseed))
    seed = 8123
    call random_seed(put=seed)
    u = 0.0_dp
    sum_u = 0.0_dp
    sum_u2 = 0.0_dp
    do i = 1, 600000
        call random_number(rnd)
        proposed = 0.5_dp*u + sqrt(0.75_dp) &
            * sqrt(-2.0_dp*log(max(rnd(1), tiny(1.0_dp)))) &
            * cos(2.0_dp*pi*rnd(2))
        if (abs(proposed) <= 8.0_dp) then
            chi = 1.0_dp + beta*u
            next_chi = 1.0_dp + beta*proposed
            log_accept = shift*(proposed-u) + gc_log_measure_change(chi, next_chi)
            if (log(max(rnd(3), tiny(1.0_dp))) <= min(0.0_dp, log_accept)) &
                u = proposed
        end if
        if (i <= 10000) cycle
        sum_u = sum_u + u
        sum_u2 = sum_u2 + u*u
    end do
    mean_u = sum_u/590000
    var_u = sum_u2/590000 - mean_u*mean_u
    ! Analytic Gaussian moments of (1+beta*u)*exp(-(u-shift)^2/2).
    if (abs(mean_u-shift-beta/(1.0_dp+beta*shift)) > 0.008_dp) error stop 4
    if (abs(var_u-1.0_dp+(beta/(1.0_dp+beta*shift))**2) > 0.015_dp) error stop 5
    print *, 'PASS: invariant-measure birth moments and MH stationary moments'
end program test_rmp_gc_measure
