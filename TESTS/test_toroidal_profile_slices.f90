program test_toroidal_profile_slices
    use, intrinsic :: iso_fortran_env, only: int64
    use toroidal_profile_slices_mod, only: toroidal_profile_bounds
    implicit none
    integer :: n1, n2, n3, ntetr, i, first, last
    integer :: case_index, previous, j
    integer :: shapes(3, 3)
    integer(int64) :: expected_first, expected_last, count
    real :: sample(18), target(18)

    shapes(:, 1) = [64, 96, 180]
    shapes(:, 2) = [64, 192, 180]
    shapes(:, 3) = [128, 192, 360]
    do case_index = 1, 3
        n1 = shapes(1, case_index)
        n2 = shapes(2, case_index)
        n3 = shapes(3, case_index)
        ! Independent mesh-count oracle in 64-bit arithmetic, including the
        ! multiplication-first formula that overflows the old implementation.
        count = 6_int64*n1*n2*n3
        ntetr = int(count)
        previous = 0
        do i = 0, n2 - 1
            call toroidal_profile_bounds(ntetr, n2, i, first, last)
            expected_first = int(i, int64)*count/n2 + 1_int64
            expected_last = int(i + 1, int64)*count/n2
            if (int(first, int64) /= expected_first) error stop 'wrong slice start'
            if (int(last, int64) /= expected_last) error stop 'wrong slice end'
            if (first /= previous + 1) error stop 'gap or overlap'
            previous = last
        end do
        if (previous /= ntetr) error stop 'incomplete mesh coverage'
    end do
    ! Small physical profile oracle: each phi slice has the same prism values.
    sample = 0.0
    sample(1:6:3) = [2.0, 5.0]
    do i = 1, 2
        call toroidal_profile_bounds(18, 3, i, first, last)
        sample(first:last:3) = sample(1:6:3)
    end do
    do j = 1, 18
        target(j) = 0.0
        if (mod(j - 1, 3) == 0) target(j) = 2.0 + 3.0*mod((j - 1)/3, 2)
    end do
    if (any(sample /= target)) error stop 'non-axisymmetric replicated profile'
    print *, 'Toroidal profile oracle: small/large mesh, coverage, replication PASS'
end program test_toroidal_profile_slices
