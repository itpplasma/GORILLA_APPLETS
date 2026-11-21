program test_marker_species

! Charge-state sampling of marker_species_mod (impurity tracing, i_option = 17) without orbits:
! markers change their charge state with the same calls as the applet (next_event_time,
! choose_transition), and the charge-state fractions at fixed times are compared with the exact
! solution of the rate equations dn_iz/dt = sum of rates in - sum of rates out.
!   1. ionisation chain Be1+ -> 2+ -> 3+ -> 4+, constant rates: closed form (Bateman 1910)
!   2. ionisation and recombination 1+ <-> 2+ <-> 3+ <-> 4+: rate equations integrated with RK4
!   3. branching of choose_transition at the threshold u = rate_ion/(rate_ion + rate_rec)
!   4. per-cell rates: the marker alternates between two cells (dwell times dwell(1), dwell(2)) with
!      different ionisation and recombination rates; transitions sampled as in the applet
!      (draw_optical_depth, time_to_event, sum of nu dt per piece, choose_transition with the cell);
!      exact: rate equations with the piecewise constant rates (RK4)
!   5. with neutrals: start as Be0, ionisation 0 -> 1 -> ... -> 4 and recombination down to 0 (RK4)
! Statistical tolerance: 4.5 binomial standard deviations sqrt(p (1 - p)/N) per fraction.

    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use marker_species_mod, only: species_set_t, init_species_set, element_data, total_rate, &
                                  next_event_time, choose_transition, set_cell_rates, cell_rates, &
                                  draw_optical_depth, time_to_event
    use thermal_force_rng_mod, only: rng_t, rng_seed, rng_uniform

    implicit none

    integer, parameter :: n_markers = 20000, n_times = 5, z_be = 4
    real(dp), parameter :: amu_g = 1.66053906660e-24_dp, n_sigma = 4.5_dp
    real(dp), parameter :: times(n_times) = [2.0e-6_dp, 5.0e-6_dp, 1.0e-5_dp, 2.0e-5_dp, 4.0e-5_dp]
    real(dp), parameter :: a = 2.0e5_dp, b = 1.0e5_dp, c = 5.0e4_dp   ! chain rates 1->2, 2->3, 3->4
    type(species_set_t) :: ss
    type(rng_t) :: rng
    real(dp), parameter :: dwell(2) = [3.0e-6_dp, 2.0e-6_dp]          ! case 4: time in cell 1, 2
    real(dp) :: rate_ion(0:z_be), rate_rec(0:z_be), frac(0:z_be, n_times), exact(0:z_be, n_times), r
    real(dp) :: cri(0:z_be, 2), crr(0:z_be, 2)
    integer :: z_nuc, n_fail

    n_fail = 0

    call element_data('Be', z_nuc, r)
    if (z_nuc /= z_be .or. abs(r - 9.0121831_dp) > 1.0e-6_dp) call fail('element_data Be')

    ! 1. pure ionisation chain
    rate_ion = 0.0_dp; rate_rec = 0.0_dp
    rate_ion(1:3) = [a, b, c]
    call make_set(rate_ion, rate_rec)
    call sample(1, 11_int64, frac)
    exact = bateman_chain(times)
    call compare('chain', frac, exact)
    call rng_seed(rng, 1_int64, 1_int64)
    if (next_event_time(ss, 4, 1.0_dp, rng) /= huge(1.0_dp)) call fail('next_event_time without rates')

    ! 2. ionisation and recombination among the ions
    rate_ion = 0.0_dp; rate_rec = 0.0_dp
    rate_ion(1:3) = [1.0e5_dp, 1.5e5_dp, 5.0e4_dp]
    rate_rec(2:4) = [3.0e5_dp, 1.0e5_dp, 2.0e5_dp]
    call make_set(rate_ion, rate_rec)
    call sample(2, 12_int64, frac)
    exact = rate_equation_rk4(2, times)
    call compare('ionisation+recombination', frac, exact)

    ! 3. branching threshold of choose_transition from charge state 2
    r = rate_ion(2)/(rate_ion(2) + rate_rec(2))
    if (choose_transition(ss, 2, r*(1.0_dp - 1.0e-12_dp)) /= 3) call fail('choose_transition below threshold')
    if (choose_transition(ss, 2, r*(1.0_dp + 1.0e-12_dp)) /= 1) call fail('choose_transition above threshold')
    if (choose_transition(ss, 1, 0.999_dp) /= 2) call fail('choose_transition, ionisation only')

    ! 4. per-cell rates (cell 1 ionising, cell 2 mainly recombining)
    cri = 0.0_dp; crr = 0.0_dp
    cri(1:3, 1) = [4.0e5_dp, 2.0e5_dp, 1.0e5_dp]
    crr(2:4, 1) = [2.0e4_dp, 1.0e4_dp, 1.0e4_dp]
    cri(1:3, 2) = [2.0e4_dp, 1.0e4_dp, 0.0_dp]
    crr(2:4, 2) = [3.0e5_dp, 4.0e5_dp, 2.0e5_dp]
    call set_cell_rates(ss, cri, crr)
    call sample_cells(1, 13_int64, frac)
    exact = rate_equation_rk4(1, times)
    call compare('per-cell rates', frac, exact)
    r = cri(2, 2)/(cri(2, 2) + crr(2, 2))
    if (choose_transition(ss, 2, r*(1.0_dp + 1.0e-12_dp), 2) /= 1) call fail('choose_transition with cell')
    if (choose_transition(ss, 2, r*(1.0_dp + 1.0e-12_dp)) /= 3) call fail('choose_transition, uniform rates')

    ! 5. with neutrals (uniform rates; a new set without the cell rates of case 4)
    rate_ion = 0.0_dp; rate_rec = 0.0_dp
    rate_ion(0:3) = [3.0e5_dp, 1.0e5_dp, 1.5e5_dp, 5.0e4_dp]
    rate_rec(1:4) = [2.0e5_dp, 3.0e5_dp, 1.0e5_dp, 2.0e5_dp]
    call make_set(rate_ion, rate_rec)
    call sample(0, 14_int64, frac)
    exact = rate_equation_rk4(0, times)
    call compare('with neutrals', frac, exact)
    r = rate_ion(1)/(rate_ion(1) + rate_rec(1))
    if (choose_transition(ss, 1, r*(1.0_dp + 1.0e-12_dp)) /= 0) call fail('choose_transition to the neutral')

    if (n_fail > 0) then
        print '(a,i0,a)', 'test_marker_species: ', n_fail, ' check(s) FAILED'
        error stop 1
    end if
    print '(a)', 'test_marker_species: all checks passed'

contains

subroutine make_set(rate_ion, rate_rec)

    real(dp), intent(in) :: rate_ion(0:), rate_rec(0:)

    ! D+ background (only enters the collision coefficients, not used here)
    call init_species_set(ss, 'Be', [2.01410177812_dp*amu_g], [1.0_dp], [1.0e13_dp], [20.0_dp], 20.0_dp, &
                          rate_ion, rate_rec, .true.)

end subroutine make_set

subroutine sample(iz_start, seed, frac)

    ! Fractions of markers in charge state 0 ... z_be at the given times; one random stream per
    ! marker as in the applet

    integer, intent(in) :: iz_start
    integer(int64), intent(in) :: seed
    real(dp), intent(out) :: frac(0:z_be, n_times)
    type(rng_t) :: rng
    real(dp) :: t, t_event
    integer :: ip, iz, k, count(0:z_be, n_times)

    count = 0
    do ip = 1, n_markers
        call rng_seed(rng, seed, int(ip, int64))
        iz = iz_start
        t = 0.0_dp
        k = 1
        do while (k <= n_times)
            t_event = next_event_time(ss, iz, t, rng)
            do while (k <= n_times)
                if (times(k) >= t_event) exit
                count(iz, k) = count(iz, k) + 1
                k = k + 1
            end do
            if (k > n_times) exit
            iz = choose_transition(ss, iz, rng_uniform(rng))
            if (iz < 0 .or. iz > z_be) call fail('charge state out of range')
            t = t_event
        end do
    end do
    frac = real(count, dp)/real(n_markers, dp)

end subroutine sample

pure function cell_at(t) result(icell)

    ! case 4: cell of the marker at time t (periodic dwell(1) in cell 1, then dwell(2) in cell 2)

    real(dp), intent(in) :: t
    integer :: icell

    icell = 2
    if (modulo(t, sum(dwell)) < dwell(1)) icell = 1

end function cell_at

pure function seg_start(iseg) result(t)

    ! case 4: start time of cell visit iseg (0, 1, 2, ...)

    integer, intent(in) :: iseg
    real(dp) :: t

    t = (iseg/2)*sum(dwell) + mod(iseg, 2)*dwell(1)

end function seg_start

subroutine sample_cells(iz_start, seed, frac)

    ! As sample, with per-cell rates: the time line is cut at the cell changes and the optical
    ! depth is accumulated piece by piece (the applet does the same per tetrahedron)

    integer, intent(in) :: iz_start
    integer(int64), intent(in) :: seed
    real(dp), intent(out) :: frac(0:z_be, n_times)
    type(rng_t) :: rng
    real(dp) :: t, t_end, t_push, tau, tau_acc, nu_ion, nu_rec, nu
    integer :: ip, iz, k, icell, iseg, count(0:z_be, n_times)

    count = 0
    do ip = 1, n_markers
        call rng_seed(rng, seed, int(ip, int64))
        iz = iz_start
        t = 0.0_dp
        k = 1
        iseg = 0                       ! cell visit: even = cell 1, odd = cell 2
        tau = draw_optical_depth(rng)
        tau_acc = 0.0_dp
        do while (k <= n_times)
            icell = 1 + mod(iseg, 2)
            t_end = min(seg_start(iseg + 1), times(k))
            call cell_rates(ss, iz, icell, nu_ion, nu_rec)
            nu = nu_ion + nu_rec
            t_push = min(t_end - t, time_to_event(nu, tau_acc, tau))
            t = t + t_push
            tau_acc = tau_acc + nu*t_push
            if (nu > 0.0_dp .and. tau_acc >= tau*(1.0_dp - 1.0e-12_dp)) then
                iz = choose_transition(ss, iz, rng_uniform(rng), icell)
                if (iz < 0 .or. iz > z_be) call fail('charge state out of range')
                tau = draw_optical_depth(rng)
                tau_acc = 0.0_dp
            end if
            if (t >= seg_start(iseg + 1)*(1.0_dp - 1.0e-12_dp)) iseg = iseg + 1
            do while (k <= n_times)
                if (t < times(k)*(1.0_dp - 1.0e-12_dp)) exit
                count(iz, k) = count(iz, k) + 1
                k = k + 1
            end do
        end do
    end do
    frac = real(count, dp)/real(n_markers, dp)

end subroutine sample_cells

function bateman_chain(t) result(n)

    ! Chain 1 -> 2 -> 3 -> 4 with distinct rates a, b, c, all markers in 1 at t = 0

    real(dp), intent(in) :: t(:)
    real(dp) :: n(0:z_be, size(t))

    n(0,:) = 0.0_dp
    n(1,:) = exp(-a*t)
    n(2,:) = a/(b - a)*(exp(-a*t) - exp(-b*t))
    n(3,:) = a*b*(exp(-a*t)/((b - a)*(c - a)) + exp(-b*t)/((a - b)*(c - b)) + exp(-c*t)/((a - c)*(b - c)))
    n(4,:) = 1.0_dp - n(1,:) - n(2,:) - n(3,:)

end function bateman_chain

function rate_equation_rk4(iz_start, t) result(n)

    ! Rate equations for the rates in ss, all markers in iz_start at t = 0, classical RK4

    integer, intent(in) :: iz_start
    real(dp), intent(in) :: t(:)
    real(dp) :: n(0:z_be, size(t)), y(0:z_be), k1(0:z_be), k2(0:z_be), k3(0:z_be), k4(0:z_be), h, tt
    integer :: k, nstep, i

    y = 0.0_dp
    y(iz_start) = 1.0_dp
    tt = 0.0_dp
    do k = 1, size(t)
        nstep = max(1, ceiling((t(k) - tt)/1.0e-8_dp))
        h = (t(k) - tt)/nstep
        do i = 1, nstep
            ! rates of the cell at the step midpoint (cell changes fall on step boundaries)
            k1 = rhs(y, tt + (i - 0.5_dp)*h)
            k2 = rhs(y + 0.5_dp*h*k1, tt + (i - 0.5_dp)*h)
            k3 = rhs(y + 0.5_dp*h*k2, tt + (i - 0.5_dp)*h)
            k4 = rhs(y + h*k3, tt + (i - 0.5_dp)*h)
            y = y + h/6.0_dp*(k1 + 2.0_dp*k2 + 2.0_dp*k3 + k4)
        end do
        tt = t(k)
        n(:,k) = y
    end do

end function rate_equation_rk4

function rhs(y, t) result(dy)

    ! uniform rates, or the cell rates at time t (case 4)

    real(dp), intent(in) :: y(0:z_be), t
    real(dp) :: dy(0:z_be), nu_ion, nu_rec
    integer :: iz

    dy = 0.0_dp
    do iz = 0, z_be
        call cell_rates(ss, iz, cell_at(t), nu_ion, nu_rec)
        dy(iz) = dy(iz) - (nu_ion + nu_rec)*y(iz)
        if (iz < z_be) dy(iz + 1) = dy(iz + 1) + nu_ion*y(iz)
        if (iz > 0) dy(iz - 1) = dy(iz - 1) + nu_rec*y(iz)
    end do

end function rhs

subroutine compare(label, frac, exact)

    character(len=*), intent(in) :: label
    real(dp), intent(in) :: frac(0:z_be, n_times), exact(0:z_be, n_times)
    real(dp) :: sigma
    integer :: k, iz

    print '(a,a)', 'case ', label
    do k = 1, n_times
        print '(a,es10.3,a,5f8.4,a,5f8.4)', '  t =', times(k), ' s  sampled', frac(:,k), '  exact', exact(:,k)
        do iz = 0, z_be
            sigma = sqrt(max(exact(iz,k)*(1.0_dp - exact(iz,k)), 1.0e-12_dp)/n_markers)
            if (abs(frac(iz,k) - exact(iz,k)) > n_sigma*sigma) call fail(label//': fraction outside tolerance')
        end do
    end do

end subroutine compare

subroutine fail(msg)

    character(len=*), intent(in) :: msg

    print '(a,a)', 'FAIL: ', msg
    n_fail = n_fail + 1

end subroutine fail

end program test_marker_species
