module russian_roulette_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

    type russian_roulette_t
    real(dp)                               :: weights_before_redistribution
    real(dp)                               :: v0
    integer                                :: maximum_storage = 0
    logical                                :: boole_russian_roulette = .false.
    end type russian_roulette_t

    type(russian_roulette_t) :: rr

    type local_rr_t
    logical :: boole_eliminated = .false.
    integer,  dimension(:),   allocatable   :: multiplicity
    real(dp), dimension(:),   allocatable   :: v
    real(dp), dimension(:),   allocatable   :: weight
    real(dp), dimension(:,:), allocatable   :: particle_state
    end type local_rr_t

    integer :: n_particle_state

contains

subroutine prepare_russian_roulette(v0,weights_before_redistribution,num_particle_state)

    !May be called repeatedly, e.g. once per species with a different v0
    real(dp), intent(in) :: v0, weights_before_redistribution
    integer, intent(in) :: num_particle_state

    n_particle_state = num_particle_state
    rr%weights_before_redistribution = weights_before_redistribution
    rr%v0 = v0

end subroutine prepare_russian_roulette

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Preparation of russian roulette finished, start playing russian rouletette
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine play_russian_roulette(weight,v,particle_state,local_rr)

    real(dp), intent(in) :: v
    real(dp),dimension(n_particle_state), intent(in) :: particle_state
    real(dp), intent(inout) :: weight
    type(local_rr_t) :: local_rr
    real(dp) :: factor, xi
    real(dp) :: pi = 3.141592653589793238462643383_dp

    factor = weight/(rr%weights_before_redistribution*14035.0_dp*2/sqrt(pi)/(1+(v/rr%v0)**7))

    if (factor.lt.0.5_dp) then !play roulette and maybe eliminate particle
        call random_number(xi)
        if (xi.lt.factor) then
            weight =  weight/factor
        else
            local_rr%boole_eliminated = .true.
        endif
    elseif (factor.gt.2.0_dp) then
        call split_particle(weight,v,particle_state,local_rr,factor)
    endif

end subroutine play_russian_roulette

subroutine split_particle(weight,v,particle_state,local_rr,splitting_number)

    real(dp), intent(in) :: v,splitting_number
    real(dp),dimension(n_particle_state), intent(in) :: particle_state
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
        local_rr%weight(splitting_id) = weight/nearest_ints(j)
        local_rr%multiplicity(splitting_id) = nearest_ints(j)-1
        local_rr%particle_state(splitting_id,:) = particle_state
        !$omp critical
        rr%maximum_storage = max(rr%maximum_storage,splitting_id)
        !$omp end critical
        weight = weight/nearest_ints(j)
    endif

end subroutine split_particle

subroutine initiate_local_rr(local_rr, i)

    integer, intent(in) :: i
    type(local_rr_t) :: local_rr

    local_rr%boole_eliminated = .false.

    if (.not.allocated(local_rr%multiplicity))            allocate(local_rr%multiplicity(i))
    if (.not.allocated(local_rr%weight))                  allocate(local_rr%weight(i))
    local_rr%multiplicity = 0
    if (.not.allocated(local_rr%v))                       allocate(local_rr%v(i))
    
    if (.not.allocated(local_rr%particle_state))    allocate(local_rr%particle_state(i,n_particle_state))

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
    
    rr_big%boole_eliminated = rr_small%boole_eliminated
    rr_big%multiplicity(1:size_small) = rr_small%multiplicity
    rr_big%v(1:size_small) = rr_small%v
    rr_big%weight(1:size_small) = rr_small%weight

    rr_big%particle_state(1:size_small,:) = rr_small%particle_state

    rr_big%multiplicity(size_small+1:size_big) = 0
    rr_big%v(size_small+1:size_big) = 0.0_dp

end subroutine copy_local_rr

subroutine move_allocation(rr_small,rr_big)

    type(local_rr_t) :: rr_small,rr_big

    call move_alloc(rr_big%multiplicity, rr_small%multiplicity)
    call move_alloc(rr_big%v, rr_small%v)
    call move_alloc(rr_big%weight, rr_small%weight)
    call move_alloc(rr_big%particle_state, rr_small%particle_state)

end subroutine move_allocation

subroutine prepare_next_split_particle(local_rr,id)

    type(local_rr_t):: local_rr
    integer, intent(out) :: id
    integer, allocatable :: nonzero_multiplicities(:)
    integer :: i

    nonzero_multiplicities = pack([(i,i=1,size(local_rr%v))],local_rr%multiplicity>0)

    if (size(nonzero_multiplicities).eq.0) then
        local_rr%boole_eliminated = .true.
    else
        i = maxloc(local_rr%v(nonzero_multiplicities),dim=1)
        id = nonzero_multiplicities(i)
        local_rr%multiplicity(id) = local_rr%multiplicity(id) -1
        local_rr%boole_eliminated = .false.
    endif

end subroutine prepare_next_split_particle

end module russian_roulette_mod