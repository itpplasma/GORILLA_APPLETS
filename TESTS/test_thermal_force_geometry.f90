program test_thermal_force_geometry

! cylindrical_step of thermal_force_mod: the coordinate step must move a point
! (R, phi, Z) by exactly the given physical displacement in Cartesian space.
! The linearised step (dx_R, dx_phi/R, dx_Z) misses dx_phi^2/(2R) in R; the
! test also checks that this is the size of the difference.

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use thermal_force_mod, only: cylindrical_step

    implicit none

    real(dp) :: x(3), dx(3), dq(3), u(7), err, max_err, max_lin
    real(dp) :: lab_dx(3), moved(3)
    integer :: i, n_seed
    integer, allocatable :: seed(:)
    logical :: all_passed

    all_passed = .true.
    max_err = 0.0_dp
    max_lin = 0.0_dp
    call random_seed(size=n_seed)
    allocate(seed(n_seed))
    seed = 20261116
    call random_seed(put=seed)
    do i = 1, 10000
        call random_number(u)
        x = [10.0_dp + 1.0e3_dp*u(1), 6.283185307179586_dp*u(2), -50.0_dp + 100.0_dp*u(3)]
        ! between 0.1 % and 10 % of R (smaller steps lose digits to rounding
        ! in the difference of two Cartesian positions of size R)
        dx = (u(4:6) - 0.5_dp)*0.2_dp*x(1)*(0.01_dp + 0.99_dp*u(7))
        dq = cylindrical_step(x, dx)
        ! physical displacement in Cartesian components
        lab_dx = dx(1)*[cos(x(2)), sin(x(2)), 0.0_dp] + dx(2)*[-sin(x(2)), cos(x(2)), 0.0_dp] &
               + dx(3)*[0.0_dp, 0.0_dp, 1.0_dp]
        moved = cartesian(x + dq) - cartesian(x)
        err = norm2(moved - lab_dx)/max(norm2(dx), tiny(1.0_dp))
        max_err = max(max_err, err)
        ! linearised step: R error dx_phi^2/(2R) to leading order
        max_lin = max(max_lin, abs((dq(1) - dx(1)) - dx(2)**2/(2.0_dp*x(1)))/max(dx(2)**2/(2.0_dp*x(1)), tiny(1.0_dp)))
    end do
    print '(a,es10.2)', 'cylindrical_step: max rel. error of the Cartesian displacement = ', max_err
    print '(a,es10.2)', 'linearised step misses dx_phi^2/(2R) in R; max rel. deviation from it = ', max_lin
    if (max_err > 1.0e-11_dp) all_passed = .false.
    if (max_lin > 0.2_dp) all_passed = .false.      ! higher orders for steps up to 10 % of R

    if (.not. all_passed) then
        print *, 'FAIL: cylindrical_step'
        error stop 1
    end if
    print *, 'PASS: cylindrical_step'

contains

    pure function cartesian(q) result(c)
        real(dp), intent(in) :: q(3)
        real(dp) :: c(3)
        c = [q(1)*cos(q(2)), q(1)*sin(q(2)), q(3)]
    end function cartesian

end program test_thermal_force_geometry
