!> Optional axisymmetric first wall / limiter contour used as loss boundary.
!>
!> The contour is a closed polygon in the poloidal (R, Z) plane, read from a
!> text file with two columns R [cm], Z [cm] per vertex ('#' comments allowed;
!> the closing segment back to the first vertex is implied). A g-file
!> limiter (RLIM, ZLIM in m) converts by a factor 100. A marker whose
!> guiding centre leaves the polygon (even-odd rule) is lost there.
module wall_contour_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

    private
    public :: load_wall_contour, wall_contains, wall_active, clear_wall_contour

    real(dp), allocatable :: r_wall(:), z_wall(:)

contains

subroutine load_wall_contour(filename)
    character(len=*), intent(in) :: filename
    integer :: u, ios, n, i
    character(len=512) :: line

    call clear_wall_contour()
    open(newunit=u, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'wall contour: cannot open loss_wall_file'
    n = 0
    do
        read(u, '(a)', iostat=ios) line
        if (ios /= 0) exit
        line = adjustl(line)
        if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
        n = n + 1
    end do
    if (n < 3) error stop 'wall contour: need at least 3 vertices'
    allocate(r_wall(n), z_wall(n))
    rewind(u)
    i = 0
    do
        read(u, '(a)', iostat=ios) line
        if (ios /= 0) exit
        line = adjustl(line)
        if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
        i = i + 1
        read(line, *, iostat=ios) r_wall(i), z_wall(i)
        if (ios /= 0) error stop 'wall contour: invalid vertex line'
    end do
    close(u)
    print '(a, i0, a)', ' Loss wall contour: ', n, ' vertices'
end subroutine load_wall_contour

logical function wall_active()
    wall_active = allocated(r_wall)
end function wall_active

!> Even-odd point-in-polygon test of (R, Z); true when no wall is loaded.
logical function wall_contains(r, z) result(inside)
    real(dp), intent(in) :: r, z
    integer :: i, j, n

    inside = .true.
    if (.not. allocated(r_wall)) return
    inside = .false.
    n = size(r_wall)
    j = n
    do i = 1, n
        if ((z_wall(i) > z) .neqv. (z_wall(j) > z)) then
            if (r < (r_wall(j) - r_wall(i))*(z - z_wall(i))/(z_wall(j) - z_wall(i)) &
                + r_wall(i)) inside = .not. inside
        end if
        j = i
    end do
end function wall_contains

subroutine clear_wall_contour()
    if (allocated(r_wall)) deallocate(r_wall, z_wall)
end subroutine clear_wall_contour

end module wall_contour_mod
