module russian_roulette_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

    type russian_roulette_t
    real(dp), dimension(:), allocatable    :: velocity_bounds
    real(dp), dimension(:), allocatable    :: roulette_numbers
    integer, dimension(:), allocatable     :: boundary_fluxes_plus
    integer, dimension(:), allocatable     :: boundary_fluxes_minus
    integer                                :: maximum_storage = 0
    logical                                :: boole_russian_roulette = .true.
    logical                                :: boole_weight_windows = .true.
    end type russian_roulette_t

    type(russian_roulette_t) :: rr

    type local_rr_t
    logical :: boole_eliminated = .false.
    real(dp), dimension(:), allocatable     :: weight_deposits
    integer,  dimension(:),   allocatable   :: multiplicity
    real(dp), dimension(:),   allocatable   :: v
    real(dp), dimension(:),   allocatable   :: weight
    real(dp), dimension(:,:), allocatable   :: particle_state_reals
    integer, dimension(:,:), allocatable    :: particle_state_integers
    end type local_rr_t

    integer :: n_reals, n_integers

contains

subroutine prepare_russian_roulette(v0,v_max,weights_before_redistribution, num_reals, num_integers)

    real(dp), intent(in) :: v0, v_max, weights_before_redistribution
    integer, intent(in) :: num_reals, num_integers
    integer :: n_domains

    n_domains = 10
    n_reals = num_reals
    n_integers = num_integers

    if (.not.rr%boole_weight_windows) call prepare_rr_with_boundary_fluxes(n_domains,v0,v_max)
    if (     rr%boole_weight_windows) call prepare_rr_with_weight_windows(n_domains,v0,v_max,weights_before_redistribution)

end subroutine prepare_russian_roulette

subroutine prepare_rr_with_boundary_fluxes(n_domains,v0,v_max)

    real(dp) :: v0, v_max
    integer :: n_domains, j
    real(dp) :: integral1, integral2, normalise_integral1, normalise_integral2
    real(dp) :: lower_limit, upper_limit
    integer :: n_intervals

    allocate(rr%boundary_fluxes_plus(n_domains+1))
    rr%boundary_fluxes_plus = 0

    allocate(rr%boundary_fluxes_minus(n_domains+1))
    rr%boundary_fluxes_minus = 0

    allocate(rr%velocity_bounds(n_domains+1))
    rr%velocity_bounds = [(v_max/n_domains*(j-1), j=1,n_domains+1)]

    allocate(rr%roulette_numbers(n_domains+1))
    rr%roulette_numbers(1) = 1.0_dp
    rr%roulette_numbers(n_domains+1) = 1.0_dp

    n_intervals = 1000

    do j = 2,n_domains

        lower_limit = rr%velocity_bounds(1)/v0
        upper_limit = rr%velocity_bounds(n_domains+1)/v0
       
        normalise_integral1 = integrate_function(x2_exp_minus_x2, lower_limit, upper_limit, n_intervals)
        normalise_integral2 = integrate_function(x3_exp_neg_x, lower_limit, upper_limit, n_intervals)

        lower_limit = rr%velocity_bounds(j-1)/v0
        upper_limit = rr%velocity_bounds(j)/v0


        integral1 = integrate_function(x2_exp_minus_x2, lower_limit, upper_limit, n_intervals)/normalise_integral1
        integral2 = integrate_function(x3_exp_neg_x, lower_limit, upper_limit, n_intervals)/normalise_integral2

        rr%roulette_numbers(j) = integral2/integral1

        lower_limit = rr%velocity_bounds(j)/v0
        upper_limit = rr%velocity_bounds(j+1)/v0

        integral1 = integrate_function(x2_exp_minus_x2, lower_limit, upper_limit, n_intervals)/normalise_integral1
        integral2 = integrate_function(x3_exp_neg_x, lower_limit, upper_limit, n_intervals)/normalise_integral2

        rr%roulette_numbers(j) = rr%roulette_numbers(j)*integral1/integral2
    enddo

end subroutine prepare_rr_with_boundary_fluxes

subroutine prepare_rr_with_weight_windows(n_domains,v0,v_max,weights_before_redistribution)

    real(dp) :: v0, v_max
    integer :: n_domains, j
    real(dp) :: weights_before_redistribution, factor
    real(dp), dimension(:), allocatable :: maxwellian_integrals
    real(dp) :: lower_limit, upper_limit
    integer :: n_intervals

    allocate(rr%velocity_bounds(n_domains+1))
    rr%velocity_bounds = [(v_max/n_domains*(j-1), j=1,n_domains+1)]
    allocate(rr%roulette_numbers(n_domains))
    allocate(maxwellian_integrals(n_domains))

    n_intervals = 1000

    do j = 1,n_domains
        lower_limit = rr%velocity_bounds(j)/v0
        upper_limit = rr%velocity_bounds(j+1)/v0

        rr%roulette_numbers(j) = 1.0_dp/(((upper_limit + lower_limit)/2)**7+1)
        maxwellian_integrals(j) = integrate_function(x2_exp_minus_x2, lower_limit, upper_limit, n_intervals)
    enddo

    factor = sum(weights_before_redistribution*maxwellian_integrals/rr%roulette_numbers)
    rr%roulette_numbers = rr%roulette_numbers*factor

end subroutine prepare_rr_with_weight_windows

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Preparation of russian roulette finished, start playing russian rouletette
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine play_russian_roulette(weight,v,v_old,reals,integers,local_rr)

    real(dp), intent(in) :: v,v_old
    real(dp),dimension(n_reals), intent(in) :: reals
    integer, dimension(n_integers), intent(in) :: integers
    real(dp), intent(inout) :: weight
    type(local_rr_t) :: local_rr

    if (rr%boole_weight_windows) then 
        call play_rr_with_weight_windows(weight,v,reals, integers,local_rr)
    else
        call play_rr_with_boundary_fluxes(weight,v,v_old,reals, integers,local_rr)
    endif

end subroutine play_russian_roulette

subroutine play_rr_with_boundary_fluxes(weight,v,v_old,reals, integers,local_rr)

    real(dp), intent(in) :: v,v_old
    real(dp),dimension(n_reals), intent(in) :: reals
    integer, dimension(n_integers), intent(in) :: integers
    real(dp), intent(inout) :: weight
    type(local_rr_t) :: local_rr
    real(dp) :: roulette_number, xi, splitting_number, passing_probability
    integer :: ind_old,ind_new,n_domains,ind_boundary

    n_domains = size(rr%velocity_bounds)-1
    call binsrc(rr%velocity_bounds,1,n_domains+1,v_old,ind_old)
    call binsrc(rr%velocity_bounds,1,n_domains+1,v,ind_new) 

    !first, check if boundary is hit, if boundary is not hit, return
    if (ind_old.eq.ind_new) return

    ind_boundary = min(ind_old,ind_new)
    roulette_number = rr%roulette_numbers(ind_boundary)
    !$omp critical
    if (ind_new.gt.ind_old) then
        rr%boundary_fluxes_plus(ind_boundary) = rr%boundary_fluxes_plus(ind_boundary) + 1
    else
        rr%boundary_fluxes_minus(ind_boundary) = rr%boundary_fluxes_minus(ind_boundary) + 1
    endif
    !$omp end critical

    if (ind_new.lt.ind_old) roulette_number = 1.0_dp/roulette_number

    if (roulette_number.gt.1.0_dp) then
        passing_probability = 1.0_dp/roulette_number
        call random_number(xi)
        if (xi.lt.passing_probability) then
            weight =  weight + local_rr%weight_deposits(ind_boundary)!*roulette_number
            local_rr%weight_deposits(ind_boundary) = 0.0_dp
        else
            local_rr%boole_eliminated = .true.
            local_rr%weight_deposits(ind_boundary) = local_rr%weight_deposits(ind_boundary) + weight
        endif
    else
        splitting_number = 1.0_dp/roulette_number
        call split_particle(weight,v,reals,integers,local_rr,splitting_number)
    endif

end subroutine play_rr_with_boundary_fluxes

subroutine play_rr_with_weight_windows(weight, v, reals, integers,local_rr)

    use gorilla_applets_types_mod, only: time_t

    real(dp), intent(in) :: v
    real(dp),dimension(n_reals), intent(in) :: reals
    integer, dimension(n_integers), intent(in) :: integers
    real(dp), intent(inout) :: weight
    type(local_rr_t) :: local_rr
    real(dp) :: factor, xi
    integer :: ind_v, n_domains

    n_domains = size(rr%velocity_bounds)-1
    call binsrc(rr%velocity_bounds,1,n_domains+1,v,ind_v)
    ind_v = ind_v - 1

    factor = weight/rr%roulette_numbers(ind_v)

    if (factor.lt.0.5_dp) then !play roulette and maybe eliminate particle
        call random_number(xi)
        if (xi.lt.factor) then
            weight =  weight/factor
        else
            local_rr%boole_eliminated = .true.
        endif
    elseif (factor.gt.2.0_dp) then
        call split_particle(weight,v,reals,integers,local_rr,factor)
    endif

end subroutine play_rr_with_weight_windows

subroutine split_particle(weight,v,reals,integers,local_rr,splitting_number)

    real(dp), intent(in) :: v,splitting_number
    real(dp),dimension(n_reals), intent(in) :: reals
    integer, dimension(n_integers), intent(in) :: integers
    real(dp), intent(inout) :: weight
    type(local_rr_t) :: local_rr
    real(dp) :: xi
    integer :: j, splitting_id, i
    integer, dimension(2) :: nearest_ints

    nearest_ints = (/floor(splitting_number),ceiling(splitting_number)/)
    call random_number(xi)
    if (xi.lt. (dble(nearest_ints(2))-splitting_number)) then
        j = 1
    else
        j = 2
    endif

    splitting_id = findloc(local_rr%multiplicity,0,dim=1)

    if (splitting_id.eq.0) then
        i = size(local_rr%v)
        call enlarge_local_rr(local_rr)
        splitting_id = i+1
    endif

    if (nearest_ints(j).gt.1) then
        local_rr%v(splitting_id) = v
        local_rr%weight(splitting_id) = weight/splitting_number
        local_rr%multiplicity(splitting_id) = nearest_ints(j)-1
        local_rr%particle_state_reals(splitting_id,:) = reals
        local_rr%particle_state_integers(splitting_id,:) = integers
        !$omp critical
        rr%maximum_storage = max(rr%maximum_storage,splitting_id)
        !$omp end critical
        weight = weight/nearest_ints(j)
    endif

end subroutine split_particle

subroutine initiate_local_rr(local_rr, i)

    integer, intent(in) :: i
    type(local_rr_t) :: local_rr

    if (.not.allocated(local_rr%weight_deposits)) allocate(local_rr%weight_deposits(size(rr%boundary_fluxes_plus)))
    local_rr%weight_deposits = 0.0_dp
    local_rr%boole_eliminated = .false.

    if (.not.allocated(local_rr%multiplicity))            allocate(local_rr%multiplicity(i))
    if (.not.allocated(local_rr%weight))                  allocate(local_rr%weight(i))
    local_rr%multiplicity = 0
    if (.not.allocated(local_rr%v))                       allocate(local_rr%v(i))
    local_rr%v = 0.0_dp
    
    if (.not.allocated(local_rr%particle_state_integers)) allocate(local_rr%particle_state_integers(i,n_integers))
    if (.not.allocated(local_rr%particle_state_reals))    allocate(local_rr%particle_state_reals(i,n_reals))

end subroutine initiate_local_rr

subroutine enlarge_local_rr(local_rr)

    type(local_rr_t) :: local_rr
    type(local_rr_t) :: local_rr_for_movealloc
    integer :: i

    i = size(local_rr%v)
    call initiate_local_rr(local_rr_for_movealloc,i+100)
    call copy_local_rr(local_rr,local_rr_for_movealloc)
    call move_allocation(local_rr,local_rr_for_movealloc)

end subroutine enlarge_local_rr

subroutine copy_local_rr(rr_small,rr_big)

    type(local_rr_t) :: rr_small,rr_big
    integer :: size_small, size_big

    size_small = size(rr_small%v)
    size_big = size(rr_big%v)
    
    rr_big%weight_deposits = rr_small%weight_deposits
    rr_big%boole_eliminated = rr_small%boole_eliminated
    rr_big%multiplicity(1:size_small) = rr_small%multiplicity
    rr_big%v(1:size_small) = rr_small%v
    rr_big%weight(1:size_small) = rr_small%weight

    rr_big%particle_state_integers(1:size_small,:) = rr_small%particle_state_integers
    rr_big%particle_state_reals(1:size_small,:) = rr_small%particle_state_reals

    rr_big%multiplicity(size_small+1:size_big) = 0
    rr_big%v(size_small+1:size_big) = 0.0_dp

end subroutine copy_local_rr

subroutine move_allocation(rr_small,rr_big)

    type(local_rr_t) :: rr_small,rr_big

    call move_alloc(rr_big%multiplicity, rr_small%multiplicity)
    call move_alloc(rr_big%v, rr_small%v)
    call move_alloc(rr_big%weight, rr_small%weight)
    call move_alloc(rr_big%particle_state_integers, rr_small%particle_state_integers)
    call move_alloc(rr_big%particle_state_reals, rr_small%particle_state_reals)

end subroutine move_allocation

subroutine prepare_next_split_particle(local_rr,id)

    type(local_rr_t):: local_rr
    integer, intent(out) :: id 

    id = maxloc(local_rr%v,dim=1)

    if (local_rr%v(id).eq.0.0_dp) then
        local_rr%boole_eliminated = .true.
    else
        local_rr%multiplicity(id) = local_rr%multiplicity(id) -1
        if (local_rr%multiplicity(id).eq.0) then
             local_rr%v(id) = 0.0_dp
        endif
        local_rr%boole_eliminated = .false.
    endif

end subroutine prepare_next_split_particle

subroutine binsrc(p,nmin,nmax,xi,i)

! Finds the index  i  of the array of increasing numbers   p  with dimension  n 
! which satisfies   p(i-1) <  xi  <  p(i) . Uses binary search algorithm.

  implicit none

  integer                                :: n,nmin,nmax,i,imin,imax,k
  double precision                       :: xi
  double precision, dimension(nmin:nmax) :: p

  imin=nmin
  imax=nmax
  n=nmax-nmin

  do k=1,n
    i=(imax-imin)/2+imin
    if(p(i).gt.xi) then
      imax=i
    else
      imin=i
    endif
    if(imax.eq.imin+1) exit
  enddo

  i=imax

  return
end subroutine binsrc

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

function x3_exp_neg_x(x) result(f)

    implicit none
    real(dp), intent(in) :: x
    real(dp) :: f
    
    f = x**7 * exp(-x**2)
end function x3_exp_neg_x

function x2_exp_minus_x2(x) result(f)

    implicit none
    real(dp), intent(in) :: x
    real(dp) :: f
    
    f = x**2 * exp(-x**2)
end function x2_exp_minus_x2

end module russian_roulette_mod