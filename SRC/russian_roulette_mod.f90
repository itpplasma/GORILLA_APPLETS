module russian_roulette_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use gorilla_applets_types_mod, only: particle_status_t, time_t

    implicit none

    type russian_roulette_t
    real(dp), dimension(:), allocatable    :: velocity_bounds
    real(dp), dimension(:), allocatable    :: roulette_numbers
    integer(dp), dimension(:), allocatable :: boundary_fluxes
    logical                                :: boole_russian_roulette = .true.
    end type russian_roulette_t

    type(russian_roulette_t) :: rr

    type local_rr_t
    real(dp), dimension(:), allocatable    :: weight_deposits
    logical :: boole_eliminated = .false.
    real(dp), dimension(:,:), allocatable :: x
    real(dp), dimension(:),   allocatable :: vpar
    real(dp), dimension(:),   allocatable :: vperp
    real(dp), dimension(:),   allocatable :: weight
    integer,  dimension(:),   allocatable :: ind_tetr
    integer,  dimension(:),   allocatable :: iface
    integer,  dimension(:),   allocatable :: multiplicity
    !type(particle_status_t), dimension(:), allocatable :: particle_status
    type(time_t), dimension(:), allocatable :: t
    end type local_rr_t

contains

subroutine prepare_russian_roulette(species)

    use gorilla_applets_types_mod, only: start

    integer, intent(in) :: species
    real(dp) :: v0, v_max
    integer :: n_domains, j
    real(dp) :: integral1, integral2, normalise_integral1, normalise_integral2
    real(dp) :: lower_limit, upper_limit
    integer :: n_intervals

    n_domains = 10
    v0 = start%v0(species)
    v_max = v0*sqrt(start%epsilon_max)

    allocate(rr%boundary_fluxes(n_domains+1))
    rr%boundary_fluxes = 0

    allocate(rr%velocity_bounds(n_domains+1))
    rr%velocity_bounds = [(v_max/n_domains*(j-1), j=1,n_domains+1)]

    allocate(rr%roulette_numbers(n_domains+1))
    rr%roulette_numbers(1) = 1.0_dp
    rr%roulette_numbers(n_domains+1) = 1.0_dp

    do j = 2,n_domains

        lower_limit = rr%velocity_bounds(1)/v0
        upper_limit = rr%velocity_bounds(n_domains+1)/v0
        n_intervals = 1000
        normalise_integral1 = integrate_function(x2_exp_minus_x2, lower_limit, upper_limit, n_intervals)
        normalise_integral2 = integrate_function(x_power_1, lower_limit, upper_limit, n_intervals)

        lower_limit = rr%velocity_bounds(j-1)/v0
        upper_limit = rr%velocity_bounds(j)/v0
        n_intervals = 1000

        integral1 = integrate_function(x2_exp_minus_x2, lower_limit, upper_limit, n_intervals)/normalise_integral1
        integral2 = integrate_function(x_power_1, lower_limit, upper_limit, n_intervals)/normalise_integral2

        rr%roulette_numbers(j) = integral2/integral1

        lower_limit = rr%velocity_bounds(j)/v0
        upper_limit = rr%velocity_bounds(j+1)/v0
        n_intervals = 1000
        integral1 = integrate_function(x2_exp_minus_x2, lower_limit, upper_limit, n_intervals)/normalise_integral1
        integral2 = integrate_function(x_power_1, lower_limit, upper_limit, n_intervals)/normalise_integral2

        rr%roulette_numbers(j) = rr%roulette_numbers(j)*integral1/integral2
    enddo

    print*, rr%roulette_numbers


end subroutine prepare_russian_roulette

subroutine play_russian_roulette(vpar_save,vperp_save,vpar,vperp,t,x,particle_status,ind_tetr,iface,species,n,local_rr)

    use binsrc_mod, only: binsrc
    use gorilla_applets_types_mod, only: start, time_t

    real(dp), intent(in) :: vpar_save,vperp_save,vpar,vperp
    integer, intent(in) :: species, ind_tetr, iface, n
    type(time_t), intent(in) :: t
    type(particle_status_t), intent(in) :: particle_status
    real(dp), dimension(3) :: x
    type(local_rr_t) :: local_rr
    real(dp) :: v_old,v_new,roulette_number, xi, splitting_number, passing_probability
    integer :: ind_old,ind_new,n_domains,ind_boundary, j
    integer, dimension(2) :: nearest_ints



    v_old = sqrt(vpar_save**2+vperp_save**2)
    v_new = sqrt(vpar**2+vperp**2)

    n_domains = size(rr%velocity_bounds)-1
    call binsrc(rr%velocity_bounds,1,n_domains+1,v_old,ind_old)
    call binsrc(rr%velocity_bounds,1,n_domains+1,v_new,ind_new) 

    !first, check if boundary is hit, if boundary is not hit, return
    if (ind_old.eq.ind_new) return

    ind_boundary = min(ind_old,ind_new)
    roulette_number = rr%roulette_numbers(ind_boundary)

    !$omp critical
    if (ind_new.gt.ind_old) then
        rr%boundary_fluxes(ind_boundary) = rr%boundary_fluxes(ind_boundary) + 1
    else
        rr%boundary_fluxes(ind_boundary) = rr%boundary_fluxes(ind_boundary) - 1
    endif
    if (ind_new.lt.ind_old) roulette_number = 1.0_dp/roulette_number
    !$omp end critical

    if (roulette_number.gt.1.0_dp) then
        passing_probability = 1.0_dp/roulette_number
        call random_number(xi)
        if (xi.lt.passing_probability) then
            !$omp critical
            start%weight(n,species) =  start%weight(n,species) + local_rr%weight_deposits(ind_boundary)
            local_rr%weight_deposits(ind_boundary) = 0.0_dp
            !$omp end critical
        else
            local_rr%boole_eliminated = .true.
            local_rr%weight_deposits(ind_boundary) = local_rr%weight_deposits(ind_boundary) + start%weight(n,species)
        endif
    else
        splitting_number = 1.0_dp/roulette_number
        nearest_ints = (/floor(splitting_number),ceiling(splitting_number)/)
        call random_number(xi)
        if (xi.lt. (dble(nearest_ints(2))-splitting_number)) then
            j = 1
        else
            j = 2
        endif
        start%weight(n,species) = start%weight(n,species)/splitting_number
        !add nearest_ints(j)-1 particles with the same weight to waiting list, for that use a lot more input arguments
    endif

end subroutine play_russian_roulette

subroutine initiate_local_rr(local_rr, j)

    integer, intent(in) :: j
    type(local_rr_t) :: local_rr

    allocate(local_rr%weight_deposits(size(rr%boundary_fluxes)))
    local_rr%weight_deposits = 0.0_dp
    local_rr%boole_eliminated = .false.
    allocate(local_rr%x(3,j))
    allocate(local_rr%vpar(j))
    allocate(local_rr%vperp(j))
    allocate(local_rr%weight(j))
    allocate(local_rr%ind_tetr(j))
    allocate(local_rr%iface(j))
    allocate(local_rr%multiplicity(j))
    local_rr%multiplicity = 0
    !allocate(local_rr%particle_status(j))
    allocate(local_rr%t(j))

end subroutine initiate_local_rr

subroutine initiate_next_split_particle(local_rr)

    type(local_rr_t):: local_rr

    local_rr%boole_eliminated = .false.

end subroutine initiate_next_split_particle

function integrate_function(func, a, b, n) result(integral)

    implicit none

    interface
        function func(xx)
            use, intrinsic :: iso_fortran_env, only: dp => real64
            real(dp), intent(in) :: xx
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

    use gorilla_applets_types_mod, only: s, in

    implicit none
    real(dp), intent(in) :: x
    real(dp) :: f
    
    f = x**2 * exp(-x**2/(s%temperature/in%energy_eV))
end function x2_exp_minus_x2

function x_power_1(x) result(f)
    implicit none
    real(dp), intent(in) :: x
    real(dp) :: f
    
    f = x**1.0_dp
end function x_power_1

end module russian_roulette_mod