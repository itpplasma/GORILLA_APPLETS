program test_multi_species_collision

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use collis_ions, only: collis_init
    use constants, only: amp, ame

    implicit none

    logical :: all_passed

    all_passed = .true.
    call test_two_species(all_passed)
    call test_five_species(all_passed)

    if (.not. all_passed) then
        print *, 'FAIL: some tests failed'
        error stop 1
    end if
    print *, 'PASS: all multi-species collision tests passed'

contains

subroutine test_two_species(all_passed)

    logical, intent(inout) :: all_passed

    integer, parameter :: n = 2
    real(dp) :: m1, z1, e0, v0
    real(dp), dimension(n) :: masses, charges, densities, temperatures
    real(dp), dimension(n) :: efcolf, velrat, enrat
    integer :: i

    m1 = 2.0_dp * amp
    z1 = 1.0_dp
    e0 = 3500.0_dp

    masses = (/2.0_dp * amp, ame/)
    charges = (/1.0_dp, -1.0_dp/)
    densities = (/3.0d13, 3.0d13/)
    temperatures = (/3500.0_dp, 3500.0_dp/)

    call collis_init(m1, z1, masses, charges, densities, temperatures, &
                     e0, v0, efcolf, velrat, enrat)

    print *, '--- 2-species test (D+ and e-) ---'
    print *, '  v0 = ', v0
    do i = 1, n
        print *, '  efcolf(', i, ') = ', efcolf(i)
    end do

    do i = 1, n
        if (efcolf(i) <= 0.0_dp) then
            print *, 'FAIL: efcolf(', i, ') <= 0 in 2-species test'
            all_passed = .false.
        end if
    end do

    if (v0 <= 0.0_dp) then
        print *, 'FAIL: v0 <= 0 in 2-species test'
        all_passed = .false.
    end if

end subroutine test_two_species

subroutine test_five_species(all_passed)

    logical, intent(inout) :: all_passed

    integer, parameter :: n2 = 2
    integer, parameter :: n5 = 5
    real(dp) :: m1, z1, e0, v0_2, v0_5
    real(dp), dimension(n2) :: masses_2, charges_2, densities_2, temperatures_2
    real(dp), dimension(n2) :: efcolf_2, velrat_2, enrat_2
    real(dp), dimension(n5) :: masses_5, charges_5, densities_5, temperatures_5
    real(dp), dimension(n5) :: efcolf_5, velrat_5, enrat_5
    real(dp) :: total_freq_2, total_freq_5
    real(dp) :: ne10_z2n
    integer :: i
    integer :: max_impurity_idx
    real(dp) :: max_impurity_z2n

    m1 = 2.0_dp * amp
    z1 = 1.0_dp
    e0 = 3500.0_dp

    masses_2 = (/2.0_dp * amp, ame/)
    charges_2 = (/1.0_dp, -1.0_dp/)
    densities_2 = (/3.0d13, 3.0d13/)
    temperatures_2 = (/3500.0_dp, 3500.0_dp/)

    call collis_init(m1, z1, masses_2, charges_2, densities_2, &
                     temperatures_2, e0, v0_2, efcolf_2, velrat_2, enrat_2)

    total_freq_2 = sum(efcolf_2)

    ! 5-species: e-, D+, Ne10+, Ne9+, Ne8+
    ! Neon mass ~ 20 amu
    masses_5(1) = 2.0_dp * amp       ! D+
    masses_5(2) = 20.0_dp * amp      ! Ne10+
    masses_5(3) = 20.0_dp * amp      ! Ne9+
    masses_5(4) = 20.0_dp * amp      ! Ne8+
    masses_5(5) = ame                 ! e-

    charges_5(1) = 1.0_dp            ! D+
    charges_5(2) = 10.0_dp           ! Ne10+
    charges_5(3) = 9.0_dp            ! Ne9+
    charges_5(4) = 8.0_dp            ! Ne8+
    charges_5(5) = -1.0_dp           ! e-

    densities_5(1) = 2.0d13          ! D+
    densities_5(2) = 1.0d12          ! Ne10+
    densities_5(3) = 5.0d11          ! Ne9+
    densities_5(4) = 3.0d11          ! Ne8+
    densities_5(5) = 3.88d13         ! e- (quasi-neutrality)

    temperatures_5 = 3500.0_dp

    call collis_init(m1, z1, masses_5, charges_5, densities_5, &
                     temperatures_5, e0, v0_5, efcolf_5, velrat_5, enrat_5)

    total_freq_5 = sum(efcolf_5)

    print *, ''
    print *, '--- 5-species test (D+, Ne10+, Ne9+, Ne8+, e-) ---'
    print *, '  v0 = ', v0_5
    do i = 1, n5
        print *, '  efcolf(', i, ') = ', efcolf_5(i)
    end do
    print *, '  total_freq_2 = ', total_freq_2
    print *, '  total_freq_5 = ', total_freq_5

    ! All efcolf must be positive
    do i = 1, n5
        if (efcolf_5(i) <= 0.0_dp) then
            print *, 'FAIL: efcolf(', i, ') <= 0 in 5-species test'
            all_passed = .false.
        end if
    end do

    ! Total collision frequency must be higher with impurities
    if (total_freq_5 <= total_freq_2) then
        print *, 'FAIL: 5-species total freq should exceed 2-species'
        print *, '  total_freq_5 = ', total_freq_5
        print *, '  total_freq_2 = ', total_freq_2
        all_passed = .false.
    end if

    ! Ne10+ (index 2) should have the highest Z^2*n among impurities
    ! (indices 2, 3, 4)
    ne10_z2n = charges_5(2)**2 * densities_5(2)
    max_impurity_idx = 2
    max_impurity_z2n = ne10_z2n
    do i = 3, 4
        if (charges_5(i)**2 * densities_5(i) > max_impurity_z2n) then
            max_impurity_z2n = charges_5(i)**2 * densities_5(i)
            max_impurity_idx = i
        end if
    end do

    if (max_impurity_idx /= 2) then
        print *, 'FAIL: Ne10+ should dominate impurity Z^2*n'
        all_passed = .false.
    end if

    ! Ne10+ collision frequency should be highest among impurities
    ! (efcolf scales with Z^2*n for similar masses/temperatures)
    if (efcolf_5(2) <= efcolf_5(3) .or. efcolf_5(2) <= efcolf_5(4)) then
        print *, 'FAIL: Ne10+ efcolf should dominate among impurities'
        all_passed = .false.
    end if

end subroutine test_five_species

end program test_multi_species_collision
