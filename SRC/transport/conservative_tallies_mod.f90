!===============================================================================
! conservative_tallies_mod.f90
!
! Cell and surface tallies plus a conservation ledger for marker ensembles.
!
! The tallies collect the physical moments (density, momentum, energy) and
! target-event information that edge applications need, decoupled from any
! particular applet's global state.  Every tally exposes its normalization
! metadata (reference marker weight, integrated time span, cell volumes) and
! talks to a shared conservation ledger, so that a client can verify that
! marker weight / momentum / energy are conserved up to the deposits that
! entered the cell and surface ledgers and the amount carried away by lost
! markers.
!
! Raw accumulations are time-weighted (accumulated value * dt), which is the
! natural form for orbit integration.  A call to normalize() converts them to
! time-averaged, optionally volume-normalized densities using the metadata.
!
! Process deposits (explicit source/sink amounts) are integrated quantities,
! not time-weighted residence, so they are accumulated in dedicated per-cell
! *_deposit ledgers and exposed separately from the time-averaged moments.
!===============================================================================
module conservative_tallies_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none
    private

    public :: tally_norm_t
    public :: cell_tally_t
    public :: surface_tally_t
    public :: conservation_ledger_t
    public :: tallies_t
    public :: TALLY_WEIGHT, TALLY_MOMENTUM, TALLY_ENERGY

    integer, parameter :: TALLY_WEIGHT   = 1
    integer, parameter :: TALLY_MOMENTUM = 2
    integer, parameter :: TALLY_ENERGY   = 3

!-------------------------------------------------------------------------------
! Normalization metadata attached to a tally.
!-------------------------------------------------------------------------------
    type :: tally_norm_t
        real(dp) :: marker_weight = 1.0_dp   ! physical weight one marker represents
        real(dp) :: total_time    = 0.0_dp   ! integrated time span of the run
        integer  :: n_markers     = 0        ! number of markers contributing
        logical  :: boole_volume_normalized  = .false.
        real(dp), allocatable :: cell_volume(:)
    contains
        procedure :: init => tally_norm_init
    end type tally_norm_t

!-------------------------------------------------------------------------------
! Per-cell conservative moments (density, momentum, energy along vpar).
!-------------------------------------------------------------------------------
    type :: cell_tally_t
        integer  :: n_cells = 0
        real(dp), allocatable :: density_raw(:)    ! weight * dt      (time-weighted)
        real(dp), allocatable :: momentum_raw(:)   ! weight * vpar * dt (time-weighted)
        real(dp), allocatable :: energy_raw(:)     ! weight * 0.5 v^2 * dt (time-weighted)
        real(dp), allocatable :: density_deposit(:)   ! integrated weight   deposited by processes
        real(dp), allocatable :: momentum_deposit(:)  ! integrated momentum deposited by processes
        real(dp), allocatable :: energy_deposit(:)    ! integrated energy   deposited by processes
        real(dp), allocatable :: density(:)        ! normalized residence
        real(dp), allocatable :: momentum(:)       ! normalized residence
        real(dp), allocatable :: energy(:)         ! normalized residence
        type(tally_norm_t) :: norm
        logical :: initialized = .false.
    contains
        procedure :: init         => cell_tally_init
        procedure :: reset        => cell_tally_reset
        procedure :: record_marker=> cell_tally_record_marker
        procedure :: record_deposit => cell_tally_record_deposit
        procedure :: normalize    => cell_tally_normalize
    end type cell_tally_t

!-------------------------------------------------------------------------------
! Target-event (surface) tallies: crossing counts and deposited quantities.
!-------------------------------------------------------------------------------
    type :: surface_tally_t
        integer  :: n_targets = 0
        real(dp), allocatable :: crossing_count(:)
        real(dp), allocatable :: deposited_weight(:)
        real(dp), allocatable :: deposited_momentum(:)
        real(dp), allocatable :: deposited_energy(:)
        type(tally_norm_t) :: norm
        logical :: initialized = .false.
    contains
        procedure :: init           => surface_tally_init
        procedure :: reset          => surface_tally_reset
        procedure :: record_crossing=> surface_tally_record_crossing
        procedure :: record_deposit => surface_tally_record_deposit
    end type surface_tally_t

!-------------------------------------------------------------------------------
! Conservation ledger.
!
! begin_from_ensemble() snapshots the initial ensemble totals.  Every deposit
! into a cell or surface ledger and every lost marker decrements the current
! totals.  A conservative process that neither deposits nor loses anything
! leaves the residual at machine precision.
!-------------------------------------------------------------------------------
    type :: conservation_ledger_t
        real(dp) :: init_weight   = 0.0_dp
        real(dp) :: init_momentum = 0.0_dp
        real(dp) :: init_energy   = 0.0_dp
        real(dp) :: weight        = 0.0_dp   ! current ensemble totals
        real(dp) :: momentum      = 0.0_dp
        real(dp) :: energy        = 0.0_dp
        real(dp) :: cell_deposited(3)   = 0.0_dp   ! weight, momentum, energy
        real(dp) :: surface_deposited(3)= 0.0_dp
        real(dp) :: lost(3)             = 0.0_dp
        logical  :: boole_initialized  = .false.
    contains
        procedure :: begin_from_totals   => ledger_begin_from_totals
        procedure :: record_cell_deposit => ledger_record_cell_deposit
        procedure :: record_surface_deposit => ledger_record_surface_deposit
        procedure :: record_lost         => ledger_record_lost
        procedure :: update_current      => ledger_update_current
        procedure :: residual            => ledger_residual
        procedure :: reset               => ledger_reset
    end type conservation_ledger_t

!-------------------------------------------------------------------------------
! Bundle of all tallies used by a transport run.
!-------------------------------------------------------------------------------
    type :: tallies_t
        type(cell_tally_t)        :: cells
        type(surface_tally_t)     :: surfaces
        type(conservation_ledger_t) :: ledger
    end type tallies_t

contains

!===============================================================================
! tally_norm_t
!===============================================================================
    subroutine tally_norm_init(self, marker_weight, total_time, n_markers, &
                                cell_volume)
        class(tally_norm_t), intent(inout) :: self
        real(dp), intent(in) :: marker_weight
        real(dp), intent(in) :: total_time
        integer,  intent(in) :: n_markers
        real(dp), intent(in), optional :: cell_volume(:)

        self%marker_weight = marker_weight
        self%total_time    = total_time
        self%n_markers     = n_markers
        if (present(cell_volume)) then
            allocate(self%cell_volume(size(cell_volume)))
            self%cell_volume = cell_volume
            self%boole_volume_normalized = .true.
        else
            self%boole_volume_normalized = .false.
        endif
    end subroutine tally_norm_init

!===============================================================================
! cell_tally_t
!===============================================================================
    subroutine cell_tally_init(self, n_cells)
        class(cell_tally_t), intent(inout) :: self
        integer, intent(in) :: n_cells

        if (allocated(self%density_raw))  deallocate(self%density_raw)
        if (allocated(self%momentum_raw)) deallocate(self%momentum_raw)
        if (allocated(self%energy_raw))   deallocate(self%energy_raw)
        if (allocated(self%density_deposit))  deallocate(self%density_deposit)
        if (allocated(self%momentum_deposit)) deallocate(self%momentum_deposit)
        if (allocated(self%energy_deposit))   deallocate(self%energy_deposit)
        if (allocated(self%density))      deallocate(self%density)
        if (allocated(self%momentum))     deallocate(self%momentum)
        if (allocated(self%energy))       deallocate(self%energy)
        self%n_cells = n_cells
        allocate(self%density_raw(n_cells))
        allocate(self%momentum_raw(n_cells))
        allocate(self%energy_raw(n_cells))
        allocate(self%density_deposit(n_cells))
        allocate(self%momentum_deposit(n_cells))
        allocate(self%energy_deposit(n_cells))
        allocate(self%density(n_cells))
        allocate(self%momentum(n_cells))
        allocate(self%energy(n_cells))
        self%density_raw  = 0.0_dp
        self%momentum_raw = 0.0_dp
        self%energy_raw   = 0.0_dp
        self%density_deposit  = 0.0_dp
        self%momentum_deposit = 0.0_dp
        self%energy_deposit   = 0.0_dp
        self%density      = 0.0_dp
        self%momentum     = 0.0_dp
        self%energy       = 0.0_dp
        self%initialized  = .true.
    end subroutine cell_tally_init

    subroutine cell_tally_reset(self)
        class(cell_tally_t), intent(inout) :: self
        if (.not. self%initialized) return
        self%density_raw  = 0.0_dp
        self%momentum_raw = 0.0_dp
        self%energy_raw   = 0.0_dp
        self%density_deposit  = 0.0_dp
        self%momentum_deposit = 0.0_dp
        self%energy_deposit   = 0.0_dp
        self%density      = 0.0_dp
        self%momentum     = 0.0_dp
        self%energy       = 0.0_dp
        call self%norm%init(1.0_dp, 0.0_dp, 0)
    end subroutine cell_tally_reset

!-------------------------------------------------------------------------------
! Record the time-weighted residence of one marker in a cell for a step dt.
! Thread-safe: individual array elements are updated atomically.
!-------------------------------------------------------------------------------
    subroutine cell_tally_record_marker(self, cell, weight, vpar, vperp, dt)
        use omp_lib, only: omp_in_parallel
        class(cell_tally_t), intent(inout) :: self
        integer,  intent(in) :: cell
        real(dp), intent(in) :: weight, vpar, vperp, dt
        real(dp) :: dv2, einc

        if (.not. self%initialized) then
            error stop 'conservative_tallies_mod: cell tally not initialized'
        endif
        if (cell < 1 .or. cell > self%n_cells) then
            ! Markers outside the computational domain are not tallied in cells;
            ! they are handled by the surface tally / lost ledger instead.
            return
        endif
        dv2  = vpar*vpar + vperp*vperp
        einc = 0.5_dp * weight * dv2 * dt
        !$omp atomic
        self%density_raw(cell) = self%density_raw(cell) + weight*dt
        !$omp atomic
        self%momentum_raw(cell) = self%momentum_raw(cell) + weight*vpar*dt
        !$omp atomic
        self%energy_raw(cell) = self%energy_raw(cell) + einc
    end subroutine cell_tally_record_marker

!-------------------------------------------------------------------------------
! Record an explicit deposit into a cell ledger (e.g. a process source/sink).
!-------------------------------------------------------------------------------
    subroutine cell_tally_record_deposit(self, cell, amount)
        class(cell_tally_t), intent(inout) :: self
        integer,  intent(in) :: cell
        real(dp), intent(in) :: amount(3)   ! weight, momentum, energy

        if (cell < 1 .or. cell > self%n_cells) return
        ! Process deposits are integrated amounts (not time-weighted), so they
        ! are accumulated in the dedicated deposit ledgers, never mixed into
        ! the time-weighted residence arrays that normalize() averages.
        !$omp atomic
        self%density_deposit(cell) = self%density_deposit(cell) + amount(TALLY_WEIGHT)
        !$omp atomic
        self%momentum_deposit(cell) = self%momentum_deposit(cell) + amount(TALLY_MOMENTUM)
        !$omp atomic
        self%energy_deposit(cell) = self%energy_deposit(cell) + amount(TALLY_ENERGY)
    end subroutine cell_tally_record_deposit

!-------------------------------------------------------------------------------
! Convert raw time-weighted accumulations into normalized moments.
!-------------------------------------------------------------------------------
    subroutine cell_tally_normalize(self)
        class(cell_tally_t), intent(inout) :: self
        integer :: i
        real(dp) :: norm_factor

        if (.not. self%initialized) return
        ! Time average of the residence raw arrays only: divide by
        ! (marker weight * total time).  Process deposits are integrated
        ! amounts stored separately in the *_deposit ledgers and are not
        ! folded into these time-averaged moments.
        norm_factor = self%norm%marker_weight * self%norm%total_time
        if (norm_factor <= 0.0_dp) norm_factor = 1.0_dp
        do i = 1, self%n_cells
            self%density(i)  = self%density_raw(i)  / norm_factor
            self%momentum(i) = self%momentum_raw(i) / norm_factor
            self%energy(i)   = self%energy_raw(i)   / norm_factor
            if (self%norm%boole_volume_normalized) then
                self%density(i)  = self%density(i)  / self%norm%cell_volume(i)
                self%momentum(i) = self%momentum(i) / self%norm%cell_volume(i)
                self%energy(i)   = self%energy(i)   / self%norm%cell_volume(i)
            endif
        enddo
    end subroutine cell_tally_normalize

!===============================================================================
! surface_tally_t
!===============================================================================
    subroutine surface_tally_init(self, n_targets)
        class(surface_tally_t), intent(inout) :: self
        integer, intent(in) :: n_targets

        if (allocated(self%crossing_count))     deallocate(self%crossing_count)
        if (allocated(self%deposited_weight))   deallocate(self%deposited_weight)
        if (allocated(self%deposited_momentum)) deallocate(self%deposited_momentum)
        if (allocated(self%deposited_energy))   deallocate(self%deposited_energy)
        self%n_targets = n_targets
        allocate(self%crossing_count(n_targets))
        allocate(self%deposited_weight(n_targets))
        allocate(self%deposited_momentum(n_targets))
        allocate(self%deposited_energy(n_targets))
        self%crossing_count     = 0.0_dp
        self%deposited_weight   = 0.0_dp
        self%deposited_momentum = 0.0_dp
        self%deposited_energy   = 0.0_dp
        self%initialized        = .true.
    end subroutine surface_tally_init

    subroutine surface_tally_reset(self)
        class(surface_tally_t), intent(inout) :: self
        integer :: i
        if (.not. self%initialized) return
        !$omp parallel do
        do i = 1, self%n_targets
            self%crossing_count(i)     = 0.0_dp
            self%deposited_weight(i)   = 0.0_dp
            self%deposited_momentum(i) = 0.0_dp
            self%deposited_energy(i)   = 0.0_dp
        enddo
        !$omp end parallel do
    end subroutine surface_tally_reset

!-------------------------------------------------------------------------------
! Record a surface crossing event for a target.
!-------------------------------------------------------------------------------
    subroutine surface_tally_record_crossing(self, target, weight)
        class(surface_tally_t), intent(inout) :: self
        integer,  intent(in) :: target
        real(dp), intent(in) :: weight

        if (target < 1 .or. target > self%n_targets) return
        !$omp atomic
        self%crossing_count(target) = self%crossing_count(target) + 1.0_dp
        !$omp atomic
        self%deposited_weight(target) = self%deposited_weight(target) + weight
    end subroutine surface_tally_record_crossing

!-------------------------------------------------------------------------------
! Record an explicit deposit onto a surface / target ledger.
!-------------------------------------------------------------------------------
    subroutine surface_tally_record_deposit(self, target, amount)
        class(surface_tally_t), intent(inout) :: self
        integer,  intent(in) :: target
        real(dp), intent(in) :: amount(3)   ! weight, momentum, energy

        if (target < 1 .or. target > self%n_targets) return
        !$omp atomic
        self%deposited_weight(target) = self%deposited_weight(target) + amount(TALLY_WEIGHT)
        !$omp atomic
        self%deposited_momentum(target) = self%deposited_momentum(target) + amount(TALLY_MOMENTUM)
        !$omp atomic
        self%deposited_energy(target) = self%deposited_energy(target) + amount(TALLY_ENERGY)
    end subroutine surface_tally_record_deposit

!===============================================================================
! conservation_ledger_t
!===============================================================================
    subroutine ledger_begin_from_totals(self, weight, momentum, energy)
        class(conservation_ledger_t), intent(inout) :: self
        real(dp), intent(in) :: weight, momentum, energy

        self%init_weight   = weight
        self%init_momentum = momentum
        self%init_energy   = energy
        self%weight        = weight
        self%momentum      = momentum
        self%energy        = energy
        self%cell_deposited    = 0.0_dp
        self%surface_deposited = 0.0_dp
        self%lost              = 0.0_dp
        self%boole_initialized = .true.
    end subroutine ledger_begin_from_totals

    subroutine ledger_record_cell_deposit(self, amount)
        class(conservation_ledger_t), intent(inout) :: self
        real(dp), intent(in) :: amount(3)
        !$omp critical
        self%cell_deposited(:) = self%cell_deposited(:) + amount(:)
        self%weight    = self%weight    - amount(TALLY_WEIGHT)
        self%momentum  = self%momentum  - amount(TALLY_MOMENTUM)
        self%energy    = self%energy    - amount(TALLY_ENERGY)
        !$omp end critical
    end subroutine ledger_record_cell_deposit

    subroutine ledger_record_surface_deposit(self, amount)
        class(conservation_ledger_t), intent(inout) :: self
        real(dp), intent(in) :: amount(3)
        !$omp critical
        self%surface_deposited(:) = self%surface_deposited(:) + amount(:)
        self%weight    = self%weight    - amount(TALLY_WEIGHT)
        self%momentum  = self%momentum  - amount(TALLY_MOMENTUM)
        self%energy    = self%energy    - amount(TALLY_ENERGY)
        !$omp end critical
    end subroutine ledger_record_surface_deposit

    subroutine ledger_record_lost(self, weight, momentum, energy)
        class(conservation_ledger_t), intent(inout) :: self
        real(dp), intent(in) :: weight, momentum, energy
        !$omp critical
        self%lost(TALLY_WEIGHT)   = self%lost(TALLY_WEIGHT)   + weight
        self%lost(TALLY_MOMENTUM) = self%lost(TALLY_MOMENTUM) + momentum
        self%lost(TALLY_ENERGY)   = self%lost(TALLY_ENERGY)   + energy
        self%weight    = self%weight    - weight
        self%momentum  = self%momentum  - momentum
        self%energy    = self%energy    - energy
        !$omp end critical
    end subroutine ledger_record_lost

!-------------------------------------------------------------------------------
! Update the current ensemble totals (used by non-conservative processes that
! directly alter marker weight / energy without going through a deposit).
!-------------------------------------------------------------------------------
    subroutine ledger_update_current(self, weight, momentum, energy)
        class(conservation_ledger_t), intent(inout) :: self
        real(dp), intent(in) :: weight, momentum, energy
        self%weight   = weight
        self%momentum = momentum
        self%energy   = energy
    end subroutine ledger_update_current

!-------------------------------------------------------------------------------
! Maximum relative conservation residual across the three channels.
!-------------------------------------------------------------------------------
    real(dp) function ledger_residual(self)
        class(conservation_ledger_t), intent(in) :: self
        real(dp) :: init(3), bal(3), res(3)
        integer :: i

        init(1) = self%init_weight
        init(2) = self%init_momentum
        init(3) = self%init_energy
        bal(1) = self%weight + self%cell_deposited(TALLY_WEIGHT) &
               + self%surface_deposited(TALLY_WEIGHT) + self%lost(TALLY_WEIGHT)
        bal(2) = self%momentum + self%cell_deposited(TALLY_MOMENTUM) &
               + self%surface_deposited(TALLY_MOMENTUM) + self%lost(TALLY_MOMENTUM)
        bal(3) = self%energy + self%cell_deposited(TALLY_ENERGY) &
               + self%surface_deposited(TALLY_ENERGY) + self%lost(TALLY_ENERGY)
        ledger_residual = 0.0_dp
        do i = 1, 3
            res(i) = abs(init(i) - bal(i))
            if (abs(init(i)) > 0.0_dp) res(i) = res(i) / abs(init(i))
            ledger_residual = max(ledger_residual, res(i))
        enddo
    end function ledger_residual

    subroutine ledger_reset(self)
        class(conservation_ledger_t), intent(inout) :: self
        self%init_weight   = 0.0_dp
        self%init_momentum = 0.0_dp
        self%init_energy   = 0.0_dp
        self%weight        = 0.0_dp
        self%momentum      = 0.0_dp
        self%energy        = 0.0_dp
        self%cell_deposited    = 0.0_dp
        self%surface_deposited = 0.0_dp
        self%lost              = 0.0_dp
        self%boole_initialized = .false.
    end subroutine ledger_reset

end module conservative_tallies_mod
