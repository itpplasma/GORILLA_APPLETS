!===============================================================================
! test_marker_transport.f90
!
! Unit / integration test for the reusable marker-transport library
! (SRC/transport).  Covers the verification points of issue #47:
!
!   1. Deterministic seeded result  - the same seed reproduces the same
!      ensemble trajectory.
!   2. No-collision analytical case conserves marker weight and energy via the
!      conservation ledger.
!   3. A synthetic process transfers a known amount into cell and surface
!      ledgers (exact, integer arithmetic).
!   4. Thread-count changes preserve statistical stream identity: the same
!      ensemble run with 1 and 8 OpenMP threads yields bit-identical marker
!      states because work is partitioned per marker stream.
!
! The test uses purely analytical steppers (no GORILLA grid is required).
!===============================================================================
program test_marker_transport

    use, intrinsic :: iso_fortran_env, only: int64, dp => real64, output_unit
    use rng_stream_mod, only: rng_stream_t
    use marker_transport_mod
    use conservative_tallies_mod, only: tallies_t, TALLY_WEIGHT, TALLY_MOMENTUM, &
                                        TALLY_ENERGY
    use omp_lib, only: omp_set_num_threads

    implicit none

    integer :: n_fail

    n_fail = 0

    !--- 1. Deterministic seeded streams -------------------------------------
    call test_stream_reproducibility(n_fail)

    !--- 2. No-collision analytical conservation ------------------------------
    call test_no_collision_conservation(n_fail)

    !--- 3. Synthetic process -> known cell / surface ledger transfers --------
    call test_synthetic_transfers(n_fail)

    !--- 4. Thread-count invariance of stream identity ------------------------
    call test_thread_invariance(n_fail)

    if (n_fail == 0) then
        write(output_unit, '(A)') 'PASS: all marker_transport library checks succeeded'
    else
        write(output_unit, '(A,I0,A)') 'FAIL: ', n_fail, ' check(s) failed'
        error stop 1
    endif

contains

!===============================================================================
! 1. Reproducibility of a single stream.
!===============================================================================
    subroutine test_stream_reproducibility(n_fail)
        integer, intent(inout) :: n_fail
        type(rng_stream_t) :: s1, s2
        real(dp) :: a(5), b(5)
        integer :: k

        call s1%init(42_int64, 1_int64)
        call s2%init(42_int64, 1_int64)
        do k = 1, 5
            a(k) = s1%next_u()
            b(k) = s2%next_u()
        enddo
        if (any(a /= b)) then
            n_fail = n_fail + 1
            write(output_unit, '(A)') '  FAIL: identical seed/index stream diverges'
        else
            write(output_unit, '(A)') '  ok: stream reproducible from (seed,index)'
        endif
        ! Streams with different indices must differ.
        call s2%init(42_int64, 2_int64)
        b(1) = s2%next_u()
        if (b(1) == a(1)) then
            n_fail = n_fail + 1
            write(output_unit, '(A)') '  FAIL: distinct indices produced identical draws'
        endif
    end subroutine test_stream_reproducibility

!===============================================================================
! 2. No-collision analytical case conserves weight & energy.
!===============================================================================
    subroutine test_no_collision_conservation(n_fail)
        integer, intent(inout) :: n_fail
        type(ensemble_t) :: ens
        type(tallies_t)  :: tallies
        type(transport_config_t) :: cfg
        type(termination_event_t), allocatable :: evts(:)
        real(dp) :: w0, p0, e0, w1, p1, e1, resid
        integer  :: i
        integer, parameter :: NP = 4, NSTEPS = 50
        real(dp), parameter :: DT = 0.1_dp

        call ens%init(NP, 7_int64)
        do i = 1, NP
            ens%markers(i)%x     = [0.1_dp*i, 0.2_dp*i, 0.0_dp]
            ens%markers(i)%vpar  = 0.3_dp*i
            ens%markers(i)%vperp = 0.25_dp*i
            ens%markers(i)%weight = 2.0_dp
            ens%markers(i)%cell  = 1
            ens%markers(i)%iface = 0
            ens%markers(i)%active = .true.
        enddo

        call tallies%cells%init(10)
        call tallies%surfaces%init(10)
        call tallies%cells%norm%init(1.0_dp, DT*NSTEPS, NP)
        call tallies%surfaces%norm%init(1.0_dp, DT*NSTEPS, NP)
        call ens%totals(w0, p0, e0)
        call tallies%ledger%begin_from_totals(w0, p0, e0)

        allocate(evts(NP))
        cfg%n_steps = NSTEPS
        cfg%dt      = DT
        call advance_ensemble(ens, straight_stepper, noop_process, tallies, cfg, evts)

        call ens%totals(w1, p1, e1)
        resid = tallies%ledger%residual()

        ! Weight and energy must be conserved exactly (no deposits, no losses).
        if (abs(w1 - w0) > 1.0e-14_dp) then
            n_fail = n_fail + 1
            write(output_unit, '(A,ES10.3)') '  FAIL: weight not conserved, dW=', w1-w0
        endif
        if (abs(e1 - e0) > 1.0e-14_dp * max(1.0_dp, e0)) then
            n_fail = n_fail + 1
            write(output_unit, '(A,ES10.3)') '  FAIL: energy not conserved, dE=', e1-e0
        endif
        if (resid > 1.0e-13_dp) then
            n_fail = n_fail + 1
            write(output_unit, '(A,ES10.3)') '  FAIL: conservation ledger residual =', resid
        else
            write(output_unit, '(A,ES10.3)') '  ok: no-collision conservation, ledger residual =', resid
        endif
        ! Normalization metadata must be exposed and correct.
        if (tallies%cells%norm%total_time /= DT*NSTEPS .or. &
            tallies%cells%norm%n_markers /= NP) then
            n_fail = n_fail + 1
            write(output_unit, '(A)') '  FAIL: cell tally normalization metadata wrong'
        endif
    end subroutine test_no_collision_conservation

!===============================================================================
! 3. Synthetic process transfers a known amount into cell and surface ledgers.
!===============================================================================
    subroutine test_synthetic_transfers(n_fail)
        integer, intent(inout) :: n_fail
        type(ensemble_t) :: ens
        type(tallies_t)  :: tallies
        type(transport_config_t) :: cfg
        type(termination_event_t), allocatable :: evts(:)
        real(dp) :: w0, p0, e0
        integer  :: i
        integer, parameter :: NSTEPS = 100
        real(dp), parameter :: DT = 0.01_dp
        ! Expected: synthetic process deposits 1.0 weight every step,
        ! alternating between cell 2 and surface 3.
        real(dp) :: exp_cell = 50.0_dp, exp_surf = 50.0_dp

        call ens%init(1, 123_int64)
        ens%markers(1)%x      = [0.0_dp, 0.0_dp, 0.0_dp]
        ens%markers(1)%vpar   = 1.0_dp
        ens%markers(1)%vperp  = 1.0_dp
        ens%markers(1)%weight = 1.0e6_dp
        ens%markers(1)%cell   = 1
        ens%markers(1)%iface  = 0
        ens%markers(1)%active = .true.

        call tallies%cells%init(10)
        call tallies%surfaces%init(10)
        call ens%totals(w0, p0, e0)
        call tallies%ledger%begin_from_totals(w0, p0, e0)

        allocate(evts(1))
        cfg%n_steps = NSTEPS
        cfg%dt      = DT
        call advance_ensemble(ens, cell_walk_stepper, synthetic_process, tallies, cfg, evts)

        if (tallies%ledger%cell_deposited(TALLY_WEIGHT) /= exp_cell) then
            n_fail = n_fail + 1
            write(output_unit, '(A,ES12.4,A,ES12.4)') '  FAIL: cell ledgers weight=', &
                tallies%ledger%cell_deposited(TALLY_WEIGHT), ' expected ', exp_cell
        endif
        if (tallies%ledger%surface_deposited(TALLY_WEIGHT) /= exp_surf) then
            n_fail = n_fail + 1
            write(output_unit, '(A,ES12.4,A,ES12.4)') '  FAIL: surface ledger weight=', &
                tallies%ledger%surface_deposited(TALLY_WEIGHT), ' expected ', exp_surf
        endif
        ! Conservation: marker weight dropped by NSTEPS (1.0 per step).
        if (ens%markers(1)%weight /= w0 - NSTEPS) then
            n_fail = n_fail + 1
            write(output_unit, '(A)') '  FAIL: marker weight did not drop by 1.0 per step'
        endif
        ! The whole transfer is accounted for in the ledger.
        if (tallies%ledger%residual() > 1.0e-13_dp) then
            n_fail = n_fail + 1
            write(output_unit, '(A,ES10.3)') '  FAIL: ledger residual after transfers =', &
                tallies%ledger%residual()
        else
            write(output_unit, '(A)') '  ok: synthetic process deposited 50u (cell) + 50u (surface), ledger balances'
        endif
        ! Cell-walk stepper changes the cell every step -> NSTEPS-1 total
        ! crossings (the first step is not a crossing).
        if (nint(sum(tallies%surfaces%crossing_count)) /= NSTEPS - 1) then
            n_fail = n_fail + 1
            write(output_unit, '(A,I0)') '  FAIL: expected NSTEPS-1 surface crossings, got ', &
                nint(sum(tallies%surfaces%crossing_count))
        else
            write(output_unit, '(A,I0)') '  ok: recorded ', NSTEPS-1, ' surface crossings'
        endif
    end subroutine test_synthetic_transfers

!===============================================================================
! 4. Thread-count changes preserve statistical stream identity.
!===============================================================================
    subroutine test_thread_invariance(n_fail)
        integer, intent(inout) :: n_fail
        type(ensemble_t) :: ens1, ens2
        type(tallies_t)  :: tallies1, tallies2
        type(transport_config_t) :: cfg
        type(termination_event_t), allocatable :: evts(:)
        integer :: i, nbad
        integer, parameter :: NP = 32, NSTEPS = 20
        real(dp), parameter :: DT = 0.05_dp

        call make_pitch_ensemble(ens1, NP)
        call tallies1%cells%init(10)
        call tallies1%surfaces%init(10)
        allocate(evts(NP))
        cfg%n_steps = NSTEPS
        cfg%dt      = DT

        ! Run with 1 thread.
        ens2 = ens1                      ! deep copy
        call tallies2%cells%init(10)
        call tallies2%surfaces%init(10)
        call omp_set_num_threads(1)
        call advance_ensemble(ens2, cell_walk_stepper, pitch_scatter_process, tallies2, cfg, evts)

        ! Run with 8 threads.
        call tallies1%cells%init(10)     ! zero tallies (already zero)
        call tallies1%surfaces%init(10)
        call omp_set_num_threads(8)
        call advance_ensemble(ens1, cell_walk_stepper, pitch_scatter_process, tallies1, cfg, evts)

        nbad = 0
        do i = 1, NP
            if (any(ens1%markers(i)%x /= ens2%markers(i)%x)) nbad = nbad + 1
            if (ens1%markers(i)%vpar /= ens2%markers(i)%vpar) nbad = nbad + 1
            if (ens1%markers(i)%vperp /= ens2%markers(i)%vperp) nbad = nbad + 1
            if (any(tallies1%cells%density_raw /= tallies2%cells%density_raw)) nbad = nbad + 1
        enddo
        if (nbad == 0) then
            write(output_unit, '(A)') '  ok: 1-thread and 8-thread runs are bit-identical (stream identity preserved)'
        else
            n_fail = n_fail + 1
            write(output_unit, '(A,I0,A)') '  FAIL: ', nbad, ' markers differ across thread counts'
        endif
    end subroutine test_thread_invariance

!===============================================================================
! Helpers
!===============================================================================
    subroutine make_pitch_ensemble(ens, n)
        type(ensemble_t), intent(inout) :: ens
        integer, intent(in) :: n
        integer :: i
        call ens%init(n, 999_int64)
        do i = 1, n
            ens%markers(i)%x      = [0.01_dp*i, 0.0_dp, 0.0_dp]
            ens%markers(i)%vpar   = 0.5_dp + 0.01_dp*i
            ens%markers(i)%vperp  = 1.5_dp
            ens%markers(i)%weight = 1.0_dp
            ens%markers(i)%cell   = 1
            ens%markers(i)%iface  = 0
            ens%markers(i)%active = .true.
        enddo
    end subroutine make_pitch_ensemble

    character(len=16) function int8str(v)
        integer(kind=8) :: tmp
        integer, intent(in) :: v
        tmp = v
        write(int8str, '(I0)') tmp
    end function int8str

!===============================================================================
! Test steppers / processes (module procedures for explicit interfaces).
!===============================================================================
    subroutine straight_stepper(x, vpar, vperp, dt, boole_initialized, cell, iface, t_remain, boole_lost)
        real(dp), intent(inout) :: x(3)
        real(dp), intent(inout) :: vpar, vperp
        real(dp), intent(in)    :: dt
        logical,  intent(inout) :: boole_initialized
        integer,  intent(inout) :: cell, iface
        real(dp), intent(out)   :: t_remain
        logical,  intent(out)   :: boole_lost
        x(3) = x(3) + vpar*dt
        boole_initialized = .true.
        cell = 1
        iface = 0
        t_remain = 0.0_dp
        boole_lost = .false.
    end subroutine straight_stepper

    subroutine cell_walk_stepper(x, vpar, vperp, dt, boole_initialized, cell, iface, t_remain, boole_lost)
        real(dp), intent(inout) :: x(3)
        real(dp), intent(inout) :: vpar, vperp
        real(dp), intent(in)    :: dt
        logical,  intent(inout) :: boole_initialized
        integer,  intent(inout) :: cell, iface
        real(dp), intent(out)   :: t_remain
        logical,  intent(out)   :: boole_lost
        boole_initialized = .true.
        if (cell < 1) cell = 1
        cell = cell + 1
        if (cell > 10) cell = 1
        iface = 0
        t_remain = 0.0_dp
        boole_lost = .false.
    end subroutine cell_walk_stepper

    subroutine synthetic_process(m, dt, rng, transfer)
        type(marker_t),       intent(inout) :: m
        real(dp),             intent(in)    :: dt
        type(rng_stream_t),   intent(inout) :: rng
        type(process_transfer_t), intent(out) :: transfer
        ! Transfer exactly 1.0 weight every step, alternating cell / surface.
        transfer%amount = 0.0_dp
        transfer%amount(TALLY_WEIGHT) = 1.0_dp
        if (mod(m%n_steps, 2) == 0) then
            transfer%target_cell = 2
            transfer%target_surface = 0
        else
            transfer%target_cell = 0
            transfer%target_surface = 3
        endif
        transfer%boole_terminate = .false.
        m%weight = m%weight - 1.0_dp
    end subroutine synthetic_process

    subroutine pitch_scatter_process(m, dt, rng, transfer)
        type(marker_t),       intent(inout) :: m
        real(dp),             intent(in)    :: dt
        type(rng_stream_t),   intent(inout) :: rng
        type(process_transfer_t), intent(out) :: transfer
        real(dp) :: u1, u2
        transfer%target_cell = 0
        transfer%target_surface = 0
        transfer%amount = 0.0_dp
        transfer%boole_terminate = .false.
        u1 = rng%next_u()
        u2 = rng%next_u()
        m%vpar  = m%vpar  * (0.5_dp + u1)
        m%vperp = m%vperp * (0.5_dp + u2)
    end subroutine pitch_scatter_process

end program test_marker_transport
