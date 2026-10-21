program test_canonical_ntv_progress
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use canonical_ntv_mod, only: ntv_incremental, ntv_n_check, ntv_check_dt, &
        canonical_ntv_allocate, canonical_ntv_progress_open, canonical_ntv_end, &
        canonical_ntv_write
    implicit none
    real(dp) :: partial(21, 4), final(21, 4), row(21)
    character(len=2048) :: line
    integer :: n, u, count, ios

    ntv_incremental = .true.
    ntv_n_check = 2
    ntv_check_dt = 0.1_dp
    call canonical_ntv_allocate(4)
    call canonical_ntv_progress_open(1.0_dp, 0.2_dp)
    !$omp parallel do schedule(dynamic, 1)
    do n = 1, 4
        call canonical_ntv_end(n, mod(n, 2) == 0, real(n, dp)*0.03_dp)
    end do
    !$omp end parallel do
    ! Read while the writer is open: each finished record must already be visible.
    open(newunit=u, file='canonical_ntv_markers.dat.partial', status='old')
    count = 0
    do
        read(u, '(a)', iostat=ios) line
        if (ios /= 0) exit
        if (line(1:1) == '#') cycle
        read(line, *) row
        n = nint(row(1))
        partial(:, n) = row
        count = count + 1
    end do
    close(u)
    if (count /= 4) error stop 'progress stream did not flush all finished markers'
    call canonical_ntv_write(1.0_dp, 0.2_dp)
    open(newunit=u, file='canonical_ntv_markers.dat', status='old')
    do
        read(u, '(a)', iostat=ios) line
        if (ios /= 0) exit
        if (line(1:1) == '#') cycle
        read(line, *) row
        n = nint(row(1))
        final(:, n) = row
    end do
    close(u)
    if (any(partial /= final)) error stop 'progress and final marker records differ'
    print *, 'PASS: concurrent progress records equal final records'
end program test_canonical_ntv_progress
