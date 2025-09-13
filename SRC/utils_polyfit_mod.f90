module utils_polyfit_mod

  implicit none

contains

  subroutine quadratic_fit(n, x, y, a0, a1, a2, success)
    ! Least-squares fit of y ≈ a2*x^2 + a1*x + a0
    integer, intent(in) :: n
    real   (kind=8), intent(in) :: x(n), y(n)
    real   (kind=8), intent(out) :: a0, a1, a2
    logical, intent(out) :: success

    integer :: degree
    real(kind=8) :: coeff(3)

    degree = 2
    call polyfit_ls(n, x, y, degree, coeff, success)

    if (success) then
      a0 = coeff(1)
      a1 = coeff(2)
      a2 = coeff(3)
    else
      a0 = 0.0d0
      a1 = 0.0d0
      a2 = 0.0d0
    end if
  end subroutine quadratic_fit

  subroutine cubic_fit(n, x, y, a0, a1, a2, a3, success)
    ! Least-squares fit of y ≈ a3*x^3 + a2*x^2 + a1*x + a0
    integer, intent(in) :: n
    real   (kind=8), intent(in) :: x(n), y(n)
    real   (kind=8), intent(out) :: a0, a1, a2, a3
    logical, intent(out) :: success

    integer :: degree
    real(kind=8) :: coeff(4)

    degree = 3
    call polyfit_ls(n, x, y, degree, coeff, success)

    if (success) then
      a0 = coeff(1)
      a1 = coeff(2)
      a2 = coeff(3)
      a3 = coeff(4)
    else
      a0 = 0.0d0
      a1 = 0.0d0
      a2 = 0.0d0
      a3 = 0.0d0
    end if
  end subroutine cubic_fit

  subroutine polyfit_ls(n, x, y, degree, coeff, success)
    ! Build normal equations for polynomial LSQ and solve
    integer, intent(in) :: n
    real   (kind=8), intent(in) :: x(n), y(n)
    integer, intent(in) :: degree              ! 2 (quadratic) or 3 (cubic)
    real   (kind=8), intent(out) :: coeff(*)   ! size degree+1, coeff(1)=a0
    logical, intent(out) :: success

    integer :: m, i, j, k
    integer :: nmat
    real(kind=8) :: eps
    real(kind=8), allocatable :: S(:)   ! sums of powers of x
    real(kind=8), allocatable :: T(:)   ! sums of x^k * y
    real(kind=8), allocatable :: A(:,:), b(:), c(:)
    real(kind=8) :: powx

    success = .false.

    if (degree < 1 .or. degree > 3) return
    nmat = degree + 1
    if (n < nmat) return

    allocate(S(0:2*degree))
    allocate(T(0:degree))
    allocate(A(nmat,nmat), b(nmat), c(nmat))

    S = 0.0d0
    T = 0.0d0

    do i = 1, n
      ! Accumulate powers efficiently
      ! powx(k) = x(i)^k up to 2*degree
      
      powx = 1.0d0
      S(0) = S(0) + 1.0d0
      T(0) = T(0) + y(i)
      do k = 1, 2*degree
        powx = powx * x(i)
        S(k) = S(k) + powx
        if (k <= degree) then
          T(k) = T(k) + powx * y(i)
        end if
      end do
    end do

    do i = 1, nmat
      do j = 1, nmat
        A(i,j) = S(i-1 + j-1)
      end do
      b(i) = T(i-1)
    end do

    eps = 1.0d-12
    call solve_linear_system(A, b, nmat, c, success, eps)
    if (success) then
      do i = 1, nmat
        coeff(i) = c(i)
      end do
    end if

    deallocate(S)
    deallocate(T)
    deallocate(A)
    deallocate(b)
    deallocate(c)
  end subroutine polyfit_ls

  subroutine solve_linear_system(A, b, nmat, x, success, eps)
    ! Gaussian elimination with partial pivoting for small dense systems
    integer, intent(in) :: nmat
    real(kind=8), intent(inout) :: A(nmat,nmat)
    real(kind=8), intent(inout) :: b(nmat)
    real(kind=8), intent(out) :: x(nmat)
    logical, intent(out) :: success
    real(kind=8), intent(in), optional :: eps

    integer :: i, j, k, p
    real(kind=8) :: pivot, factor, tol, maxval
    integer :: maxrow

    tol = 1.0d-12
    if (present(eps)) tol = eps

    success = .false.

    do k = 1, nmat
      ! Find pivot row
      maxval = 0.0d0
      maxrow = k
      do p = k, nmat
        if (abs(A(p,k)) > maxval) then
          maxval = abs(A(p,k))
          maxrow = p
        end if
      end do

      if (maxval <= tol) return

      if (maxrow /= k) then
        call swap_rows(A, b, nmat, k, maxrow)
      end if

      pivot = A(k,k)
      do i = k + 1, nmat
        factor = A(i,k) / pivot
        A(i,k) = 0.0d0
        do j = k + 1, nmat
          A(i,j) = A(i,j) - factor * A(k,j)
        end do
        b(i) = b(i) - factor * b(k)
      end do
    end do

    ! Back substitution
    do i = nmat, 1, -1
      if (abs(A(i,i)) <= tol) return
      x(i) = b(i)
      do j = i + 1, nmat
        x(i) = x(i) - A(i,j) * x(j)
      end do
      x(i) = x(i) / A(i,i)
    end do

    success = .true.
  end subroutine solve_linear_system

  subroutine swap_rows(A, b, nmat, r1, r2)
    integer, intent(in) :: nmat, r1, r2
    real(kind=8), intent(inout) :: A(nmat,nmat)
    real(kind=8), intent(inout) :: b(nmat)
    real(kind=8) :: tmp
    integer :: j

    do j = 1, nmat
      tmp = A(r1,j)
      A(r1,j) = A(r2,j)
      A(r2,j) = tmp
    end do
    tmp = b(r1)
    b(r1) = b(r2)
    b(r2) = tmp
  end subroutine swap_rows

end module utils_polyfit_mod

