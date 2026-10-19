program test_thread_seed
    ! The per-thread random streams of the OpenMP particle loops must be a
    ! function of seed.inp: equal seeds give equal streams, different seeds
    ! or thread numbers give different streams.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use utils_data_pre_and_post_processing_mod, only: set_seed_for_random_numbers
    use utils_parallelised_particle_pushing_mod, only: &
        initialise_seed_for_random_numbers_for_each_thread
    implicit none
    integer, parameter :: n_draws = 16
    real(dp), dimension(n_draws) :: a0, a0_again, a1, b0

    call draw(11, 0, a0)
    call draw(11, 0, a0_again)
    call draw(11, 1, a1)
    call draw(12, 0, b0)

    if (any(a0 /= a0_again)) error stop 'same seed.inp, different thread stream'
    if (all(a0 == b0)) error stop 'thread stream ignores seed.inp'
    if (all(a0 == a1)) error stop 'threads share one stream'
    print *, 'PASS: per-thread streams follow seed.inp'

contains

    subroutine draw(seed, thread_num, values)
        integer, intent(in) :: seed, thread_num
        real(dp), dimension(:), intent(out) :: values
        integer :: unit, k

        open(newunit=unit, file='seed.inp', status='replace', action='write')
        write(unit, *) 8
        write(unit, *) (seed + 7*k, k = 1, 8)
        close(unit)
        call set_seed_for_random_numbers
        call initialise_seed_for_random_numbers_for_each_thread(thread_num, 1)
        call random_number(values)
    end subroutine draw

end program test_thread_seed
