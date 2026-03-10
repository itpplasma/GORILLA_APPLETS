module marker_distribution_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

    integer, parameter :: default_table_size = 1000

    !---------------------------------------------------------------------------
    ! Unified abstract interface for PDF functions
    ! Accepts a vector of any length: 3 elements for 3D, 1 for 1D
    !---------------------------------------------------------------------------
    abstract interface
        function pdf_interface(x) result(val)
            import :: dp
            real(dp), intent(in) :: x(:)
            real(dp) :: val
        end function pdf_interface
    end interface

    !---------------------------------------------------------------------------
    ! Distribution types
    !---------------------------------------------------------------------------

    ! 3D distribution (sampled via rejection sampling)
    type :: distribution_3d_t
        real(dp) :: xmin(3) = 0.0_dp
        real(dp) :: xmax(3) = 1.0_dp
        real(dp) :: pdf_max = 1.0_dp         ! maximum PDF value for rejection sampling
        real(dp) :: normalization = 1.0_dp   ! integral of PDF over domain (Monte Carlo estimate)
        procedure(pdf_interface), pointer, nopass :: pdf => null()
        logical :: initialized = .false.
        logical :: is_flat = .false.         ! true if using flat/uniform distribution
    end type distribution_3d_t

    ! 1D distribution (sampled via inverse CDF)
    type :: distribution_1d_t
        real(dp) :: xmin = 0.0_dp
        real(dp) :: xmax = 1.0_dp
        integer :: table_size = default_table_size
        real(dp), allocatable :: cdf(:)
        real(dp), allocatable :: x_values(:)
        real(dp) :: normalization = 1.0_dp
        procedure(pdf_interface), pointer, nopass :: pdf => null()
        logical :: initialized = .false.
        logical :: is_flat = .false.         ! true if using flat/uniform distribution
    end type distribution_1d_t

contains

!===============================================================================
! 3D DISTRIBUTION - uses rejection sampling
!===============================================================================

!-------------------------------------------------------------------------------
! Initialize a 3D distribution
!-------------------------------------------------------------------------------
subroutine init_distribution_3d(dist, pdf, xmin, xmax, n_probe)

    type(distribution_3d_t), intent(inout) :: dist
    procedure(pdf_interface) :: pdf
    real(dp), intent(in) :: xmin(3), xmax(3)
    integer, intent(in), optional :: n_probe

    integer :: i, n
    real(dp) :: x(3), val, volume, pdf_sum

    dist%pdf => pdf
    dist%xmin = xmin
    dist%xmax = xmax

    ! Check if this is a flat distribution - if so, skip Monte Carlo
    if (associated(dist%pdf, pdf_flat)) then
        dist%pdf_max = 1.0_dp
        dist%normalization = product(xmax - xmin)
        dist%initialized = .true.
        dist%is_flat = .true.
        return
    endif

    ! Non-flat distribution: do full Monte Carlo estimation
    dist%is_flat = .false.

    ! Number of probe points for Monte Carlo estimation
    n = 10000
    if (present(n_probe)) n = n_probe

    ! Estimate pdf_max and normalization by probing random points
    dist%pdf_max = 0.0_dp
    pdf_sum = 0.0_dp
    do i = 1, n
        call random_number(x)
        x = xmin + x * (xmax - xmin)
        val = pdf(x)
        if (val > dist%pdf_max) dist%pdf_max = val
        pdf_sum = pdf_sum + val
    enddo
    ! Add safety margin to pdf_max
    dist%pdf_max = dist%pdf_max * 1.1_dp

    ! Monte Carlo estimate of integral: (volume / n) * sum(pdf)
    volume = product(xmax - xmin)
    dist%normalization = volume * pdf_sum / dble(n)

    dist%initialized = .true.

end subroutine init_distribution_3d

!-------------------------------------------------------------------------------
! Cleanup 3D distribution
!-------------------------------------------------------------------------------
subroutine cleanup_distribution_3d(dist)

    type(distribution_3d_t), intent(inout) :: dist

    dist%pdf => null()
    dist%initialized = .false.

end subroutine cleanup_distribution_3d

!-------------------------------------------------------------------------------
! Evaluate normalized 3D PDF at a given point
!-------------------------------------------------------------------------------
function evaluate_distribution_3d(dist, x) result(val)

    type(distribution_3d_t), intent(in) :: dist
    real(dp), intent(in) :: x(3)
    real(dp) :: val

    if (.not. dist%initialized) then
        print *, 'Error: evaluate_distribution_3d called on uninitialized distribution'
        error stop 'Distribution not initialized'
    endif

    val = dist%pdf(x)
    if (dist%normalization > 0.0_dp) then
        val = val / dist%normalization
    endif

end function evaluate_distribution_3d

!-------------------------------------------------------------------------------
! Sample a single point from the 3D distribution using rejection sampling
!-------------------------------------------------------------------------------
function sample_distribution_3d(dist) result(x)

    type(distribution_3d_t), intent(in) :: dist
    real(dp) :: x(3)

    real(dp) :: u, val

    if (.not. dist%initialized) then
        print *, 'Error: sample_distribution_3d called on uninitialized distribution'
        error stop 'Distribution not initialized'
    endif

    do
        call random_number(x)
        x = dist%xmin + x * (dist%xmax - dist%xmin)
        val = dist%pdf(x)
        call random_number(u)
        if (u * dist%pdf_max <= val) exit
    enddo

end function sample_distribution_3d

!-------------------------------------------------------------------------------
! Fill arrays with 3D samples
!-------------------------------------------------------------------------------
subroutine sample_array_3d(dist, samples)

    type(distribution_3d_t), intent(in) :: dist
    real(dp), intent(out) :: samples(:,:)  ! shape (3, n_particles)

    integer :: i

    if (.not. dist%initialized) then
        print *, 'Error: sample_array_3d called on uninitialized distribution'
        error stop 'Distribution not initialized'
    endif

    do i = 1, size(samples, 2)
        samples(:, i) = sample_distribution_3d(dist)
    enddo

end subroutine sample_array_3d

!===============================================================================
! 1D DISTRIBUTION - uses inverse CDF sampling
!===============================================================================

!-------------------------------------------------------------------------------
! Initialize 1D distribution by building the CDF table
!-------------------------------------------------------------------------------
subroutine init_distribution_1d(dist, pdf, xmin, xmax, table_size)

    type(distribution_1d_t), intent(inout) :: dist
    procedure(pdf_interface) :: pdf
    real(dp), intent(in) :: xmin, xmax
    integer, intent(in), optional :: table_size

    integer :: i, n
    real(dp) :: x, hx

    dist%pdf => pdf
    dist%xmin = xmin
    dist%xmax = xmax

    ! Check if this is a flat distribution - if so, skip CDF table building
    if (associated(dist%pdf, pdf_flat)) then
        dist%normalization = xmax - xmin
        dist%initialized = .true.
        dist%is_flat = .true.
        return
    endif

    ! Non-flat distribution: build CDF table
    dist%is_flat = .false.

    if (present(table_size)) then
        dist%table_size = table_size
    else
        dist%table_size = default_table_size
    endif
    n = dist%table_size

    if (allocated(dist%cdf)) deallocate(dist%cdf)
    if (allocated(dist%x_values)) deallocate(dist%x_values)
    allocate(dist%cdf(0:n))
    allocate(dist%x_values(0:n))

    hx = (xmax - xmin) / dble(n)

    ! Build CDF using midpoint rule
    dist%cdf(0) = 0.0_dp
    dist%x_values(0) = xmin
    do i = 1, n
        x = xmin + (dble(i) - 0.5_dp) * hx
        dist%cdf(i) = dist%cdf(i-1) + pdf([x]) * hx
        dist%x_values(i) = xmin + dble(i) * hx
    enddo

    dist%normalization = dist%cdf(n)
    if (dist%normalization > 0.0_dp) then
        dist%cdf = dist%cdf / dist%normalization
    endif

    dist%initialized = .true.

end subroutine init_distribution_1d

!-------------------------------------------------------------------------------
! Cleanup 1D distribution
!-------------------------------------------------------------------------------
subroutine cleanup_distribution_1d(dist)

    type(distribution_1d_t), intent(inout) :: dist

    if (allocated(dist%cdf)) deallocate(dist%cdf)
    if (allocated(dist%x_values)) deallocate(dist%x_values)
    dist%pdf => null()
    dist%initialized = .false.

end subroutine cleanup_distribution_1d

!-------------------------------------------------------------------------------
! Evaluate normalized 1D PDF at a given point
!-------------------------------------------------------------------------------
function evaluate_distribution_1d(dist, x) result(val)

    type(distribution_1d_t), intent(in) :: dist
    real(dp), intent(in) :: x
    real(dp) :: val

    if (.not. dist%initialized) then
        print *, 'Error: evaluate_distribution_1d called on uninitialized distribution'
        error stop 'Distribution not initialized'
    endif

    val = dist%pdf([x])
    if (dist%normalization > 0.0_dp) then
        val = val / dist%normalization
    endif

end function evaluate_distribution_1d

!-------------------------------------------------------------------------------
! Sample a single value from the 1D distribution
!-------------------------------------------------------------------------------
function sample_distribution_1d(dist) result(x)

    type(distribution_1d_t), intent(in) :: dist
    real(dp) :: x

    real(dp) :: xi, t
    integer :: k

    if (.not. dist%initialized) then
        print *, 'Error: sample_distribution_1d called on uninitialized distribution'
        error stop 'Distribution not initialized'
    endif

    call random_number(xi)

    ! For flat distributions, use direct uniform sampling
    if (dist%is_flat) then
        x = dist%xmin + xi * (dist%xmax - dist%xmin)
        return
    endif

    call binsrc_local(dist%cdf, 0, dist%table_size, xi, k)

    ! Linear interpolation
    if (k > 0 .and. dist%cdf(k) > dist%cdf(k-1)) then
        t = (xi - dist%cdf(k-1)) / (dist%cdf(k) - dist%cdf(k-1))
        x = dist%x_values(k-1) + t * (dist%x_values(k) - dist%x_values(k-1))
    else
        x = dist%x_values(k)
    endif

end function sample_distribution_1d

!-------------------------------------------------------------------------------
! Fill array with 1D samples
!-------------------------------------------------------------------------------
subroutine sample_array_1d(dist, n, samples)

    type(distribution_1d_t), intent(in) :: dist
    integer, intent(in) :: n
    real(dp), intent(out) :: samples(n)

    integer :: i

    if (.not. dist%initialized) then
        print *, 'Error: sample_array_1d called on uninitialized distribution'
        error stop 'Distribution not initialized'
    endif

    do i = 1, n
        samples(i) = sample_distribution_1d(dist)
    enddo

end subroutine sample_array_1d

!===============================================================================
! UTILITY FUNCTIONS
!===============================================================================

!-------------------------------------------------------------------------------
! Binary search: finds k such that p(k-1) < xi <= p(k)
!-------------------------------------------------------------------------------
subroutine binsrc_local(p, nmin, nmax, xi, i)

    integer, intent(in) :: nmin, nmax
    real(dp), intent(in) :: p(nmin:nmax)
    real(dp), intent(in) :: xi
    integer, intent(out) :: i

    integer :: imin, imax, k, n

    imin = nmin
    imax = nmax
    n = nmax - nmin

    do k = 1, n
        i = (imax - imin) / 2 + imin
        if (p(i) > xi) then
            imax = i
        else
            imin = i
        endif
        if (imax == imin + 1) exit
    enddo

    i = imax

end subroutine binsrc_local

!===============================================================================
! FLAT/UNIFORM PDF FUNCTION
! Works for any dimension: 3D or 1D
!===============================================================================

function pdf_flat(x) result(val)
    real(dp), intent(in) :: x(:)
    real(dp) :: val
    val = 1.0_dp
end function pdf_flat

end module marker_distribution_mod
