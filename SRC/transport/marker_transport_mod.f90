!===============================================================================
! marker_transport_mod.f90
!
! Reusable marker-ensemble transport services.
!
! The library advances a marker ensemble (positions, velocities, weights)
! through a *supplied* stepper.  A stepper is any routine with the
! step_signature_t interface -- in practice a thin wrapper around GORILLA's
! orbit_timestep_gorilla / characteristic-step contract.  An optional process
! callback (process_signature_t) couples collisions or other physics to the
! stepper and communicates any resulting source/sink/termination back to the
! caller through a process_transfer_t record.
!
! The transport loop is OpenMP-parallel over markers.  Each marker owns its own
! reproducible random stream (rng_stream_mod), so the collective result is
! independent of the thread count as long as concurrent work is partitioned per
! marker and never per draw.
!
! The module deliberately owns NO application globals: all configuration is
! passed in through transport_config_t, all state lives in the ensemble, and
! diagnostics go to a caller-supplied tallies_t plus a per-marker termination
! event array.
!===============================================================================
module marker_transport_mod

    use, intrinsic :: iso_fortran_env, only: int64, dp => real64
    use rng_stream_mod, only: rng_stream_t

    implicit none
    private

    public :: marker_t
    public :: ensemble_t
    public :: termination_event_t
    public :: process_transfer_t
    public :: transport_config_t
    public :: step_signature_t
    public :: process_signature_t
    public :: ensemble_init
    public :: ensemble_totals
    public :: advance_ensemble
    public :: noop_process
    public :: TERM_NONE, TERM_LOST, TERM_REACHED_TARGET, TERM_MAX_STEPS, TERM_PROCESS

    ! Termination event types.
    integer, parameter :: TERM_NONE           = 0
    integer, parameter :: TERM_LOST           = 1
    integer, parameter :: TERM_REACHED_TARGET = 2
    integer, parameter :: TERM_MAX_STEPS      = 3
    integer, parameter :: TERM_PROCESS        = 4

!-------------------------------------------------------------------------------
! A single computational marker.
!-------------------------------------------------------------------------------
    type :: marker_t
        real(dp) :: x(3) = 0.0_dp
        real(dp) :: vpar  = 0.0_dp
        real(dp) :: vperp = 0.0_dp
        real(dp) :: pitch = 0.0_dp
        real(dp) :: weight = 1.0_dp
        real(dp) :: time   = 0.0_dp
        integer  :: cell   = -1
        integer  :: iface  = -1
        integer  :: n_steps = 0
        logical  :: boole_initialized = .false.
        logical  :: active = .true.
        type(rng_stream_t) :: rng
    end type marker_t

!-------------------------------------------------------------------------------
! An ensemble of markers (the state advanced by the transport loop).
!-------------------------------------------------------------------------------
    type :: ensemble_t
        type(marker_t), allocatable :: markers(:)
        integer  :: n_markers = 0
        integer(int64) :: global_seed = 0_int64
        logical  :: boole_per_marker_streams = .true.
    contains
        procedure :: init => ensemble_init
        procedure :: totals => ensemble_totals
    end type ensemble_t

!-------------------------------------------------------------------------------
! Record describing why / how a marker left the ensemble (or finished).
! Slot i of the events array passed to advance_ensemble corresponds to
! marker i (markers are numbered 1..n_markers).
!-------------------------------------------------------------------------------
    type :: termination_event_t
        integer  :: marker_id     = 0
        integer  :: event_type    = TERM_NONE
        real(dp) :: time          = 0.0_dp
        real(dp) :: x(3)          = 0.0_dp
        integer  :: cell          = -1
        integer  :: reason        = 0   ! free-form app specific reason
        logical  :: occurred      = .false.
    end type termination_event_t

!-------------------------------------------------------------------------------
! Source / sink transfer reported by a process callback.
!
! A transfer describes a single deposit of `amount` (weight, momentum,
! energy) that is removed from the marker.  It is destined either for a cell
! ledger (target_cell > 0) or a surface ledger (target_surface > 0); the
! conservation ledger is debited exactly once per transfer.
!-------------------------------------------------------------------------------
    type :: process_transfer_t
        integer  :: target_cell     = 0      ! 0 = none
        integer  :: target_surface  = 0      ! 0 = none
        real(dp) :: amount(3)       = 0.0_dp ! weight, momentum, energy
        logical  :: boole_terminate = .false.
        integer  :: termination_type = TERM_NONE
        integer  :: reason          = 0
    end type process_transfer_t

!-------------------------------------------------------------------------------
! Run configuration for advance_ensemble.
!-------------------------------------------------------------------------------
    type :: transport_config_t
        integer(kind=8) :: n_steps   = 1
        real(dp)        :: dt        = 0.0_dp   ! step size (caller must set)
        integer         :: max_steps = huge(1)
        logical         :: boole_record_cell_tallies      = .true.
        logical         :: boole_record_surface_crossings = .true.
        real(dp)        :: marker_weight = 1.0_dp
    end type transport_config_t

!-------------------------------------------------------------------------------
! Stepper interface: advances one marker by dt.
!
! On entry boole_initialized is .false. for a fresh marker; the stepper sets
! it .true. after it has established the pusher's per-tetrahedron invariants.
! cell/iface identify the current tetrahedron / exit face.  The stepper
! returns the portion of dt for which the marker remained in the domain
! (t_remain); for a fully internal step t_remain == 0.  boole_lost signals
! that the marker escaped the computational domain.
!-------------------------------------------------------------------------------
    abstract interface
        subroutine step_signature_t(x, vpar, vperp, dt, boole_initialized, &
                                    cell, iface, t_remain, boole_lost)
            import :: dp
            real(dp), intent(inout) :: x(3)
            real(dp), intent(inout) :: vpar, vperp
            real(dp), intent(in)    :: dt
            logical,  intent(inout) :: boole_initialized
            integer,  intent(inout) :: cell, iface
            real(dp), intent(out)   :: t_remain
            logical,  intent(out)   :: boole_lost
        end subroutine step_signature_t
    end interface

!-------------------------------------------------------------------------------
! Process (collision / source-sink) callback interface.
!-------------------------------------------------------------------------------
    abstract interface
        subroutine process_signature_t(m, dt, rng, transfer)
            import :: marker_t, dp, rng_stream_t, process_transfer_t
            type(marker_t),       intent(inout) :: m
            real(dp),             intent(in)    :: dt
            type(rng_stream_t),   intent(inout) :: rng
            type(process_transfer_t), intent(out) :: transfer
        end subroutine process_signature_t
    end interface

contains

!===============================================================================
! ensemble_init
!===============================================================================
    subroutine ensemble_init(self, n_markers, global_seed, boole_per_marker_streams)
        class(ensemble_t), intent(inout) :: self
        integer, intent(in) :: n_markers
        integer(int64), intent(in) :: global_seed
        logical, intent(in), optional :: boole_per_marker_streams
        integer :: i

        if (allocated(self%markers)) deallocate(self%markers)
        allocate(self%markers(n_markers))
        self%n_markers = n_markers
        self%global_seed = global_seed
        self%boole_per_marker_streams = .true.
        if (present(boole_per_marker_streams)) &
            self%boole_per_marker_streams = boole_per_marker_streams

        ! Each marker owns a deterministic stream indexed by its marker number,
        ! so its random draws do not depend on thread scheduling.
        if (self%boole_per_marker_streams) then
            !$omp parallel do
            do i = 1, n_markers
                call self%markers(i)%rng%init(global_seed, int(i, int64))
            enddo
            !$omp end parallel do
        endif
    end subroutine ensemble_init

!===============================================================================
! ensemble_totals
!===============================================================================
    subroutine ensemble_totals(self, weight, momentum, energy)
        class(ensemble_t), intent(in) :: self
        real(dp), intent(out) :: weight, momentum, energy
        integer :: i
        weight = 0.0_dp
        momentum = 0.0_dp
        energy = 0.0_dp
        !$omp parallel do reduction(+:weight,momentum,energy)
        do i = 1, self%n_markers
            weight   = weight   + self%markers(i)%weight
            momentum = momentum + self%markers(i)%weight * self%markers(i)%vpar
            energy   = energy   + self%markers(i)%weight * &
                        0.5_dp * (self%markers(i)%vpar**2 + self%markers(i)%vperp**2)
        enddo
        !$omp end parallel do
    end subroutine ensemble_totals

!===============================================================================
! noop_process - a process callback that does nothing (pure transport case).
!===============================================================================
    subroutine noop_process(m, dt, rng, transfer)
        type(marker_t),       intent(inout) :: m
        real(dp),             intent(in)    :: dt
        type(rng_stream_t),   intent(inout) :: rng
        type(process_transfer_t), intent(out) :: transfer
        transfer%target_cell     = 0
        transfer%target_surface  = 0
        transfer%amount          = 0.0_dp
        transfer%boole_terminate = .false.
        transfer%termination_type = TERM_NONE
        transfer%reason           = 0
    end subroutine noop_process

!===============================================================================
! advance_ensemble
!
! Advances every active marker for cfg%n_steps steps of cfg%dt through the
! supplied stepper, applies the optional process callback after every step,
! fills the per-marker cell/surface tallies and the conservation ledger, and
! records termination events.  `events` must have at least ens%n_markers
! entries; slot i corresponds to marker i.
!===============================================================================
    subroutine advance_ensemble(ens, stepper, process, tallies, cfg, events)
        use conservative_tallies_mod, only: tallies_t
        type(ensemble_t),      intent(inout) :: ens
        procedure(step_signature_t)    :: stepper
        procedure(process_signature_t), optional :: process
        type(tallies_t),       intent(inout) :: tallies
        type(transport_config_t), intent(in) :: cfg
        type(termination_event_t), intent(inout) :: events(:)
        type(process_transfer_t) :: transfer
        real(dp) :: dt, t_remain
        integer :: i, k, prev_cell
        logical :: boole_lost

        if (size(events) < ens%n_markers) then
            error stop 'marker_transport_mod: events array too small'
        endif

        !$omp parallel default(shared) &
        !$omp& private(i, k, dt, t_remain, boole_lost, transfer, prev_cell)
        !$omp do
        do i = 1, ens%n_markers
            events(i)%occurred  = .false.
            events(i)%event_type = TERM_NONE
            events(i)%time       = 0.0_dp
            if (.not. ens%markers(i)%active) cycle

            dt = cfg%dt
            prev_cell = -1
            do k = 1, int(cfg%n_steps)
                if (.not. ens%markers(i)%active) exit
                if (ens%markers(i)%n_steps >= cfg%max_steps) then
                    events(i)%occurred   = .true.
                    events(i)%event_type = TERM_MAX_STEPS
                    events(i)%time       = ens%markers(i)%time
                    events(i)%x          = ens%markers(i)%x
                    events(i)%cell       = ens%markers(i)%cell
                    ens%markers(i)%active = .false.
                    exit
                endif

                ! --- characteristic step through the supplied stepper ---
                call stepper(ens%markers(i)%x, ens%markers(i)%vpar, &
                             ens%markers(i)%vperp, dt, &
                             ens%markers(i)%boole_initialized, &
                             ens%markers(i)%cell, ens%markers(i)%iface, &
                             t_remain, boole_lost)
                ens%markers(i)%n_steps = ens%markers(i)%n_steps + 1
                ens%markers(i)%time    = ens%markers(i)%time + (dt - t_remain)
                ens%markers(i)%pitch   = ens%markers(i)%vpar / &
                            sqrt(ens%markers(i)%vpar**2 + ens%markers(i)%vperp**2)

                ! --- surface-crossing detection (cell change across a step) ---
                if (cfg%boole_record_surface_crossings .and. prev_cell > 0 .and. &
                    prev_cell /= ens%markers(i)%cell) then
                    call tallies%surfaces%record_crossing(ens%markers(i)%cell, &
                                    ens%markers(i)%weight)
                endif

                ! --- cell tally for this step ---
                if (cfg%boole_record_cell_tallies) then
                    call tallies%cells%record_marker(ens%markers(i)%cell, &
                            ens%markers(i)%weight, ens%markers(i)%vpar, &
                            ens%markers(i)%vperp, dt - t_remain)
                endif

                ! --- loss termination ---
                if (boole_lost .or. ens%markers(i)%cell < 1) then
                    events(i)%occurred   = .true.
                    events(i)%event_type = TERM_LOST
                    events(i)%time       = ens%markers(i)%time
                    events(i)%x          = ens%markers(i)%x
                    events(i)%cell       = ens%markers(i)%cell
                    call tallies%ledger%record_lost( &
                        ens%markers(i)%weight, &
                        ens%markers(i)%weight * ens%markers(i)%vpar, &
                        ens%markers(i)%weight * 0.5_dp * &
                            (ens%markers(i)%vpar**2 + ens%markers(i)%vperp**2))
                    ens%markers(i)%active = .false.
                    exit
                endif

                ! --- process (collision / source-sink) callback ---
                if (present(process)) then
                    call process(ens%markers(i), dt - t_remain, &
                                 ens%markers(i)%rng, transfer)
                    ! Deposit into the cell tally / ledger (if destined for a cell).
                    if (transfer%target_cell > 0) then
                        if (cfg%boole_record_cell_tallies) then
                            call tallies%cells%record_deposit(transfer%target_cell, &
                                                              transfer%amount)
                        endif
                        call tallies%ledger%record_cell_deposit(transfer%amount)
                    endif
                    ! Deposit onto a surface ledger (if destined for a surface).
                    if (transfer%target_surface > 0) then
                        call tallies%surfaces%record_deposit(transfer%target_surface, &
                                                             transfer%amount)
                        if (transfer%target_cell <= 0) then
                            call tallies%ledger%record_surface_deposit(transfer%amount)
                        endif
                    endif
                    ! Termination requested by the process.
                    if (transfer%boole_terminate) then
                        events(i)%occurred   = .true.
                        events(i)%event_type = transfer%termination_type
                        events(i)%time       = ens%markers(i)%time
                        events(i)%x          = ens%markers(i)%x
                        events(i)%cell       = ens%markers(i)%cell
                        events(i)%reason     = transfer%reason
                        ens%markers(i)%active = .false.
                        exit
                    endif
                endif

                prev_cell = ens%markers(i)%cell
            enddo
        enddo
        !$omp end do
        !$omp end parallel
    end subroutine advance_ensemble

end module marker_transport_mod
