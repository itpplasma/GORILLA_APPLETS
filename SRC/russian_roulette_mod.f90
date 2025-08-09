module russian_roulette_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

    type russian_roulette_t
    real(dp), dimension(:), allocatable :: velocity_bounds
    real(dp), dimension(:), allocatable :: roulette_numbers
    end type russian_roulette_t

    type(russian_roulette_t) :: russian_roulette

contains

subroutine prepare_russian_roulette(species)

    use gorilla_applets_types_mod, only: start

    integer, intent(in) :: species
    real(dp) :: v0, v_max
    integer :: n_domains, j
    real(dp) :: integral1, integral2, normalise_integral1, normalise_integral2
    real(dp) :: lower_limit, upper_limit
    integer :: n_intervals

    n_domains = 5
    v0 = start%v0(species)
    v_max = v0*sqrt(start%epsilon_max)

    allocate(russian_roulette%velocity_bounds(n_domains+1))
    russian_roulette%velocity_bounds = [(v_max/n_domains*(j-1), j=1,n_domains+1)]

    allocate(russian_roulette%roulette_numbers(n_domains+1))
    russian_roulette%roulette_numbers(1) = 0.0_dp
    russian_roulette%roulette_numbers(n_domains+1) = 0.0_dp

    do j = 2,n_domains

        lower_limit = russian_roulette%velocity_bounds(1)/v0
        upper_limit = russian_roulette%velocity_bounds(n_domains+1)/v0
        n_intervals = 1000
        normalise_integral1 = integrate_function(x2_exp_minus_x2, lower_limit, upper_limit, n_intervals)
        normalise_integral2 = integrate_function(x_power_5, lower_limit, upper_limit, n_intervals)

        lower_limit = russian_roulette%velocity_bounds(j-1)/v0
        upper_limit = russian_roulette%velocity_bounds(j)/v0
        n_intervals = 1000
        print*, lower_limit, upper_limit
        integral1 = integrate_function(x2_exp_minus_x2, lower_limit, upper_limit, n_intervals)/normalise_integral1
        integral2 = integrate_function(x_power_5, lower_limit, upper_limit, n_intervals)/normalise_integral2
        print*, integral1, integral2

        russian_roulette%roulette_numbers(j) = integral2/integral1

        lower_limit = russian_roulette%velocity_bounds(j)/v0
        upper_limit = russian_roulette%velocity_bounds(j+1)/v0
        n_intervals = 1000
        integral1 = integrate_function(x2_exp_minus_x2, lower_limit, upper_limit, n_intervals)/normalise_integral1
        integral2 = integrate_function(x_power_5, lower_limit, upper_limit, n_intervals)/normalise_integral2

        russian_roulette%roulette_numbers(j) = russian_roulette%roulette_numbers(j)*integral1/integral2
    enddo

    print*, russian_roulette%roulette_numbers
    stop

end subroutine prepare_russian_roulette

function integrate_function(func, a, b, n) result(integral)

    implicit none

    interface
        function func(x)
            use, intrinsic :: iso_fortran_env, only: dp => real64
            real(dp), intent(in) :: x
            real(dp) :: func
        end function func
    end interface

    real(dp), intent(in) :: a, b
    integer, intent(in) :: n
    real(dp) :: integral

    real(dp) :: h, x
    integer :: i

    if (n <= 0) then
        integral = 0.0_dp
        return
    end if

    h = (b - a) / real(n, dp)
    integral = 0.5_dp * (func(a) + func(b))

    do i = 1, n - 1
        x = a + real(i, dp) * h
        integral = integral + func(x)
    end do

    integral = integral * h

end function integrate_function

function x2_exp_minus_x2(x) result(f)
    implicit none
    real(dp), intent(in) :: x
    real(dp) :: f
    
    f = x**2 * exp(-x**2)
end function x2_exp_minus_x2

function x_power_5(x) result(f)
    implicit none
    real(dp), intent(in) :: x
    real(dp) :: f
    
    f = x**0.0_dp
end function x_power_5

end module russian_roulette_mod