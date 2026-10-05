program test_transport_loss_moments

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use flux_deviation_mod, only: accumulate_valid_moments

    implicit none

    real(dp) :: sum_second(3),sum_fourth(3)
    integer :: contributors(3)

    sum_second = 0.0_dp
    sum_fourth = 0.0_dp
    contributors = 0

    call accumulate_valid_moments([1.0_dp,4.0_dp,9.0_dp],3_8, &
        sum_second,sum_fourth,contributors)
    call accumulate_valid_moments([4.0_dp,16.0_dp,36.0_dp],2_8, &
        sum_second,sum_fourth,contributors)

    if(any(contributors /= [2,2,1])) error stop 'invalid contributor counts'
    if(any(abs(sum_second-[5.0_dp,20.0_dp,9.0_dp]) > 1.0e-14_dp)) then
        error stop 'lost trajectory prefix was not retained in second moment'
    endif
    if(any(abs(sum_fourth-[17.0_dp,272.0_dp,81.0_dp]) > 1.0e-14_dp)) then
        error stop 'lost trajectory prefix was not retained in fourth moment'
    endif

    print *, 'PASS: valid prefixes contribute until the recorded loss step'

end program test_transport_loss_moments
