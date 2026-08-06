!===============================================================================
! rng_stream_mod.f90
!
! Reproducible, splittable pseudo-random number streams for marker transport.
!
! The library exposes an explicit 64-bit splitmix64 generator.  Every stream is
! identified by a (base_seed, stream_index) pair and derived deterministically
! from that pair.  A stream therefore produces the identical sequence of draws
! no matter which OpenMP thread consumes it and no matter how many threads the
! surrounding ensemble run happens to use.  This is the property that makes
! "statistical stream identity" invariant under thread-count changes: as long
! as concurrent worker threads partition *streams* (one per marker), and never
! partition individual draws, the collected Monte-Carlo result is bit-for-bit
! reproducible.
!
! The generator is splitmix64 (public domain, Steele et al. 2014), a fast
! non-cryptographic 64-bit PRNG with a full-period 2^64 state space and good
! statistical behaviour for the moderate correlation-length Monte-Carlo
! workloads used in the guiding-centre applets.
!===============================================================================
module rng_stream_mod

    use, intrinsic :: iso_fortran_env, only: int64, dp => real64

    implicit none
    private

    public :: rng_stream_t
    public :: rng_mix64
    public :: rng_seed_from_index

    ! splitmix64 constants.
    integer(int64), parameter :: GOLDEN_GAMMA = int(Z'9E3779B97F4A7C15', int64)
    integer(int64), parameter :: MIX_C1       = int(Z'BF58476D1CE4E5B9', int64)
    integer(int64), parameter :: MIX_C2       = int(Z'94D049BB133111EB', int64)
    integer(int64), parameter :: SEED_SALT    = int(Z'D1342543DE82EF95', int64)

    ! 53-bit fraction mask / conversion factor for uniform [0,1) draws.
    integer(int64), parameter :: FRAC53_MASK  = int(Z'001FFFFFFFFFFFFF', int64)
    real(dp),        parameter :: TWO_NEG_53   = 1.0_dp / 9007199254740992.0_dp

    type :: rng_stream_t
        private
        integer(int64) :: state        = 0_int64
        integer(int64) :: base_seed    = 0_int64
        integer(int64) :: stream_index = 0_int64
        logical        :: initialized  = .false.
    contains
        procedure :: init          => rng_init
        procedure :: reset         => rng_reset
        procedure :: next_u        => rng_next_u
        procedure :: next_int64    => rng_next_int64
        procedure :: index         => rng_index
        procedure :: base          => rng_base
        procedure :: is_initialized=> rng_is_initialized
    end type rng_stream_t

contains

!-------------------------------------------------------------------------------
! splitmix64 mixing function (pure, gives the next un-mixed 64-bit value).
!-------------------------------------------------------------------------------
    pure function rng_mix64(z) result(out)
        integer(int64), intent(in) :: z
        integer(int64) :: out
        integer(int64) :: x
        x = z
        x = ieor(x, ishft(x, -30)) * MIX_C1
        x = ieor(x, ishft(x, -27)) * MIX_C2
        x = ieor(x, ishft(x, -31))
        out = x
    end function rng_mix64

!-------------------------------------------------------------------------------
! Derive a well-mixed 64-bit seed from a (base_seed, stream_index) pair.
!-------------------------------------------------------------------------------
    pure function rng_seed_from_index(base_seed, stream_index) result(seed)
        integer(int64), intent(in) :: base_seed
        integer(int64), intent(in) :: stream_index
        integer(int64) :: seed
        seed = ieor(rng_mix64(ieor(base_seed, SEED_SALT)), &
                    rng_mix64(stream_index * GOLDEN_GAMMA))
    end function rng_seed_from_index

!-------------------------------------------------------------------------------
! Initialise a stream from a base seed and an integer stream index.
!-------------------------------------------------------------------------------
    subroutine rng_init(self, base_seed, stream_index)
        class(rng_stream_t), intent(inout) :: self
        integer(int64), intent(in) :: base_seed
        integer(int64), intent(in) :: stream_index

        self%base_seed    = base_seed
        self%stream_index = stream_index
        self%state        = rng_seed_from_index(base_seed, stream_index)
        self%initialized  = .true.
    end subroutine rng_init

!-------------------------------------------------------------------------------
! Reset a stream back to its initial state (reproducible re-run).
!-------------------------------------------------------------------------------
    subroutine rng_reset(self)
        class(rng_stream_t), intent(inout) :: self
        if (.not. self%initialized) then
            error stop 'rng_stream_mod: rng_reset called on uninitialized stream'
        endif
        self%state = rng_seed_from_index(self%base_seed, self%stream_index)
    end subroutine rng_reset

!-------------------------------------------------------------------------------
! Draw a uniform random number in [0,1).
!-------------------------------------------------------------------------------
    real(dp) function rng_next_u(self)
        class(rng_stream_t), intent(inout) :: self
        integer(int64) :: frac

        if (.not. self%initialized) then
            error stop 'rng_stream_mod: rng_next_u called on uninitialized stream'
        endif
        ! Advance the splitmix64 state.
        self%state = self%state + GOLDEN_GAMMA
        ! Use the 53 most-significant bits of the mixed state as the fraction.
        frac = iand(ishft(rng_mix64(self%state), -11), FRAC53_MASK)
        rng_next_u = real(frac, dp) * TWO_NEG_53
    end function rng_next_u

!-------------------------------------------------------------------------------
! Draw a raw 64-bit integer (useful for stream seeds / diagnostics).
!-------------------------------------------------------------------------------
    integer(int64) function rng_next_int64(self)
        class(rng_stream_t), intent(inout) :: self
        if (.not. self%initialized) then
            error stop 'rng_stream_mod: rng_next_int64 called on uninitialized stream'
        endif
        self%state = self%state + GOLDEN_GAMMA
        rng_next_int64 = rng_mix64(self%state)
    end function rng_next_int64

!-------------------------------------------------------------------------------
! Accessors.
!-------------------------------------------------------------------------------
    integer(int64) function rng_index(self)
        class(rng_stream_t), intent(in) :: self
        rng_index = self%stream_index
    end function rng_index

    integer(int64) function rng_base(self)
        class(rng_stream_t), intent(in) :: self
        rng_base = self%base_seed
    end function rng_base

    logical function rng_is_initialized(self)
        class(rng_stream_t), intent(in) :: self
        rng_is_initialized = self%initialized
    end function rng_is_initialized

end module rng_stream_mod
