module global_transport_fit_math_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

    private

    public :: build_piecewise_linear_basis
    public :: build_second_difference_matrix
    public :: compute_geometry_from_boundaries
    public :: solve_linear_system

contains

subroutine build_piecewise_linear_basis(eval_points, knot_points, basis)

    real(dp), intent(in) :: eval_points(:)
    real(dp), intent(in) :: knot_points(:)
    real(dp), allocatable, intent(out) :: basis(:, :)

    integer :: i
    integer :: j
    integer :: n_eval
    integer :: n_knots
    real(dp) :: left
    real(dp) :: right
    real(dp) :: x

    n_eval = size(eval_points)
    n_knots = size(knot_points)
    allocate(basis(n_eval, n_knots))
    basis = 0.0_dp

    do i = 1, n_eval
        x = eval_points(i)
        if (x <= knot_points(1)) then
            basis(i, 1) = 1.0_dp
        else if (x >= knot_points(n_knots)) then
            basis(i, n_knots) = 1.0_dp
        else
            do j = 1, n_knots - 1
                left = knot_points(j)
                right = knot_points(j + 1)
                if ((x >= left) .and. (x <= right)) then
                    basis(i, j) = (right - x) / (right - left)
                    basis(i, j + 1) = (x - left) / (right - left)
                    exit
                end if
            end do
        end if
    end do

end subroutine build_piecewise_linear_basis

subroutine build_second_difference_matrix(n_points, matrix)

    integer, intent(in) :: n_points
    real(dp), allocatable, intent(out) :: matrix(:, :)

    integer :: i

    if (n_points < 3) then
        allocate(matrix(0, n_points))
        return
    end if

    allocate(matrix(n_points - 2, n_points))
    matrix = 0.0_dp

    do i = 1, n_points - 2
        matrix(i, i) = 1.0_dp
        matrix(i, i + 1) = -2.0_dp
        matrix(i, i + 2) = 1.0_dp
    end do

end subroutine build_second_difference_matrix

subroutine compute_geometry_from_boundaries(boundary_s, shell_volumes, cell_centers, boundary_areas)

    real(dp), intent(in) :: boundary_s(:)
    real(dp), intent(in) :: shell_volumes(:)
    real(dp), allocatable, intent(out) :: cell_centers(:)
    real(dp), allocatable, intent(out) :: boundary_areas(:)

    real(dp), allocatable :: cell_volume_derivative(:)
    integer :: i
    integer :: n_cells
    real(dp) :: delta_s

    n_cells = size(shell_volumes)
    allocate(cell_centers(n_cells))
    allocate(boundary_areas(n_cells + 1))
    allocate(cell_volume_derivative(n_cells))

    do i = 1, n_cells
        cell_centers(i) = 0.5_dp * (boundary_s(i) + boundary_s(i + 1))
        delta_s = boundary_s(i + 1) - boundary_s(i)
        cell_volume_derivative(i) = shell_volumes(i) / max(delta_s, 1.0d-12)
    end do

    boundary_areas = 0.0_dp
    boundary_areas(1) = 0.0_dp
    if (n_cells == 1) then
        boundary_areas(2) = cell_volume_derivative(1)
        return
    end if
    do i = 2, n_cells
        boundary_areas(i) = 0.5_dp * (cell_volume_derivative(i - 1) + cell_volume_derivative(i))
    end do
    boundary_areas(n_cells + 1) = cell_volume_derivative(n_cells)

end subroutine compute_geometry_from_boundaries

subroutine solve_linear_system(matrix_in, rhs_in, solution)

    real(dp), intent(in) :: matrix_in(:, :)
    real(dp), intent(in) :: rhs_in(:, :)
    real(dp), allocatable, intent(out) :: solution(:, :)

    integer :: info
    integer :: n
    integer :: nrhs
    integer, allocatable :: ipiv(:)
    real(dp), allocatable :: matrix(:, :)

    interface
        subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            integer, intent(in) :: n
            integer, intent(in) :: nrhs
            integer, intent(in) :: lda
            integer, intent(in) :: ldb
            integer, intent(out) :: ipiv(*)
            integer, intent(out) :: info
            double precision, intent(inout) :: a(lda, *)
            double precision, intent(inout) :: b(ldb, *)
        end subroutine dgesv
    end interface

    n = size(matrix_in, 1)
    nrhs = size(rhs_in, 2)

    allocate(matrix(n, n))
    allocate(solution(n, nrhs))
    allocate(ipiv(n))

    matrix = matrix_in
    solution = rhs_in
    call dgesv(n, nrhs, matrix, n, ipiv, solution, n, info)
    if (info /= 0) then
        print *, 'Error: dgesv failed in solve_linear_system with info = ', info
        stop 1
    end if

end subroutine solve_linear_system

end module global_transport_fit_math_mod
