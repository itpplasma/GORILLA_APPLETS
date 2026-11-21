module adf11_mod

! Reader for ADAS ADF11 files (unresolved, IGRD = 1, IPRT = 1), e.g. scd96_be.dat (effective
! ionisation) and acd96_be.dat (effective recombination) from OPEN-ADAS
! (https://open.adas.ac.uk/adf11). The files are not part of this repository (OPEN-ADAS terms);
! see EXAMPLES/impurity_tracing/README.md.
!
! File layout:
!   header line      IZMAX IDMAX ITMAX IZ1MIN IZ1MAX  /ELEMENT/ ...
!   separator line
!   IDMAX values     log10 n_e [cm^-3]
!   ITMAX values     log10 T_e [eV]
!   per Z1 = IZ1MIN ... IZ1MAX:
!     separator line containing 'Z1='
!     ITMAX*IDMAX values log10 <sigma v> [cm^3/s], n_e running fastest
!   comment lines (C...)
! Z1 is the charge of the recombining ion: block Z1 of an scd file is ionisation Z1-1 -> Z1, block
! Z1 of an acd file is recombination Z1 -> Z1-1.
!
! Interpolation: bilinear in (log10 n_e, log10 T_e) on log10 <sigma v>; outside the grid the
! values at the grid edge are used (no extrapolation).

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

    private

    type, public :: adf11_t
        character(len=256) :: file = ''
        integer :: z_nuc = 0, nd = 0, nt = 0, z1_min = 0, z1_max = 0
        real(dp), allocatable :: log_ne(:), log_te(:)
        real(dp), allocatable :: log_coef(:,:,:)        ! (nd, nt, z1_min:z1_max)
    end type adf11_t

    public :: read_adf11, adf11_coef

contains

subroutine read_adf11(file, tab)

    character(len=*), intent(in) :: file
    type(adf11_t), intent(out) :: tab
    character(len=256) :: line
    integer :: iu, ios, iz1, z1, k

    open(newunit=iu, file=file, status='old', action='read', iostat=ios)
    if (ios /= 0) then
        print *, 'adf11: cannot open ', trim(file)
        error stop
    end if
    tab%file = file

    read(iu, '(a)', iostat=ios) line
    if (ios /= 0) call fail('no header line')
    k = index(line, '/')
    if (k < 2) call fail('header line without /ELEMENT/')
    read(line(1:k-1), *, iostat=ios) tab%z_nuc, tab%nd, tab%nt, tab%z1_min, tab%z1_max
    if (ios /= 0) call fail('cannot read IZMAX IDMAX ITMAX IZ1MIN IZ1MAX from the header')
    if (tab%nd < 2 .or. tab%nt < 2 .or. tab%z1_min < 1 .or. tab%z1_max > tab%z_nuc .or. tab%z1_min > tab%z1_max) &
        call fail('inconsistent header')

    read(iu, '(a)', iostat=ios) line                     ! separator
    if (ios /= 0) call fail('file ends after the header')
    allocate(tab%log_ne(tab%nd), tab%log_te(tab%nt), tab%log_coef(tab%nd, tab%nt, tab%z1_min:tab%z1_max))
    read(iu, *, iostat=ios) tab%log_ne
    if (ios /= 0) call fail('cannot read the density grid')
    read(iu, *, iostat=ios) tab%log_te
    if (ios /= 0) call fail('cannot read the temperature grid')
    if (any(tab%log_ne(2:) <= tab%log_ne(:tab%nd-1)) .or. any(tab%log_te(2:) <= tab%log_te(:tab%nt-1))) &
        call fail('grids not strictly increasing')

    do iz1 = tab%z1_min, tab%z1_max
        read(iu, '(a)', iostat=ios) line
        if (ios /= 0) call fail('file ends before all Z1 blocks')
        k = index(line, 'Z1=')
        if (k == 0) call fail('expected a block separator with Z1= (metastable-resolved files are not supported)')
        read(line(k+3:), *, iostat=ios) z1
        if (ios /= 0 .or. z1 /= iz1) call fail('Z1 blocks not in order IZ1MIN ... IZ1MAX')
        read(iu, *, iostat=ios) tab%log_coef(:, :, iz1)
        if (ios /= 0) call fail('cannot read a Z1 data block')
    end do
    close(iu)

contains

    subroutine fail(msg)
        character(len=*), intent(in) :: msg
        print *, 'adf11: ', trim(file), ': ', msg
        error stop
    end subroutine fail

end subroutine read_adf11

function adf11_coef(tab, z1, ne_cm3, te_ev) result(coef)

    ! <sigma v> [cm^3/s] of block z1 at electron density ne_cm3 [cm^-3] and temperature te_ev [eV]

    type(adf11_t), intent(in) :: tab
    integer, intent(in) :: z1
    real(dp), intent(in) :: ne_cm3, te_ev
    real(dp) :: coef
    real(dp) :: x, y, wx, wy
    integer :: i, j

    if (z1 < tab%z1_min .or. z1 > tab%z1_max) then
        print *, 'adf11: ', trim(tab%file), ': no block Z1 = ', z1
        error stop
    end if
    if (ne_cm3 <= 0.0_dp .or. te_ev <= 0.0_dp) error stop 'adf11: n_e and T_e must be positive'

    x = min(max(log10(ne_cm3), tab%log_ne(1)), tab%log_ne(tab%nd))
    y = min(max(log10(te_ev), tab%log_te(1)), tab%log_te(tab%nt))
    i = interval(tab%log_ne, x)
    j = interval(tab%log_te, y)
    wx = (x - tab%log_ne(i))/(tab%log_ne(i+1) - tab%log_ne(i))
    wy = (y - tab%log_te(j))/(tab%log_te(j+1) - tab%log_te(j))
    coef = 10.0_dp**((1.0_dp - wx)*(1.0_dp - wy)*tab%log_coef(i, j, z1) + wx*(1.0_dp - wy)*tab%log_coef(i+1, j, z1) &
                     + (1.0_dp - wx)*wy*tab%log_coef(i, j+1, z1) + wx*wy*tab%log_coef(i+1, j+1, z1))

end function adf11_coef

pure function interval(grid, x) result(i)

    ! i with grid(i) <= x <= grid(i+1), 1 <= i <= size(grid)-1, for grid(1) <= x <= grid(n)

    real(dp), intent(in) :: grid(:), x
    integer :: i, lo, hi, mid

    lo = 1
    hi = size(grid)
    do while (hi - lo > 1)
        mid = (lo + hi)/2
        if (grid(mid) <= x) then
            lo = mid
        else
            hi = mid
        end if
    end do
    i = lo

end function interval

end module adf11_mod
