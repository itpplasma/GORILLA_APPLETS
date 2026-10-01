program test_rmp_equilibrium_potential
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use profile_data_mod, only: load_profiles, eval_profiles, cleanup_profiles, &
        profile_values_t
    implicit none
    real(dp), parameter :: r(5) = [0.0_dp, 0.2_dp, 0.7_dp, 1.3_dp, 2.0_dp]
    character(len=6), parameter :: names(4) = &
        [character(len=6) :: 'n.dat', 'Te.dat', 'Ti.dat', 'Er.dat']
    type(profile_values_t) :: value
    real(dp) :: expected, field, slope
    integer :: iu, i, j, model

    ! Analytic antiderivative of E(r)=3+slope*r on an irregular radial grid.
    ! Check the production loader and spline values, not an isolated formula.
    do model = 1, 2
        slope = real(model-1, dp)*4.0_dp
        open(newunit=iu, file='potential_fixture/equil.dat', status='replace')
        do i = 1, size(r)
            write(iu, *) 10.0_dp+r(i), r(i), 2.0_dp, r(i), r(i)
        end do
        close(iu)
        do j = 1, size(names)
            open(newunit=iu, file='potential_fixture/'//trim(names(j)), &
                status='replace')
            do i = 1, size(r)
                field = 100.0_dp
                if (j == 4) field = 3.0_dp+slope*r(i)
                write(iu, *) r(i), field
            end do
            close(iu)
        end do
        call load_profiles('potential_fixture', 'potential_fixture/equil.dat')
        do i = 1, size(r)
            call eval_profiles(r(i)/r(size(r)), value)
            expected = -3.0_dp*r(i)-0.5_dp*slope*r(i)**2
            if (abs(value%Phi0-expected) > 1.0e-12_dp) error stop 1
        end do
        call cleanup_profiles()
    end do
    print *, 'PASS: equilibrium potential matches analytic constant/linear Er integral'
end program test_rmp_equilibrium_potential
