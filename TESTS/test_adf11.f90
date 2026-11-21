program test_adf11

! adf11_mod and marker_species_mod adas_rates on a synthetic ADF11 file (Be layout: Z = 4,
! blocks Z1 = 1 ... 4, irregular grids). The tabulated log10 <sigma v> is bilinear in
! (log10 n_e, log10 T_e) with a different function per block, so the interpolation must
! reproduce it exactly (inside the grid) and give the edge values outside. Real ADAS files
! are not used: OPEN-ADAS data may not be distributed with the code.

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use adf11_mod, only: adf11_t, read_adf11, adf11_coef
    use marker_species_mod, only: adas_rates

    implicit none

    integer, parameter :: nd = 5, nt = 6, z_nuc = 4
    real(dp), parameter :: log_ne(nd) = [8.0_dp, 10.0_dp, 11.5_dp, 13.0_dp, 15.0_dp]
    real(dp), parameter :: log_te(nt) = [-0.5_dp, 0.0_dp, 0.7_dp, 1.3_dp, 2.0_dp, 3.5_dp]
    character(len=*), parameter :: fscd = 'test_adf11_scd.dat', facd = 'test_adf11_acd.dat'
    type(adf11_t) :: tab
    real(dp) :: ne, te, got, want, rate_ion(0:z_nuc), rate_rec(0:z_nuc)
    integer :: z1, k, nfail

    nfail = 0
    call write_file(fscd, 0.0_dp)
    call write_file(facd, 3.0_dp)

    call read_adf11(fscd, tab)
    if (tab%z_nuc /= z_nuc .or. tab%nd /= nd .or. tab%nt /= nt .or. tab%z1_min /= 1 .or. tab%z1_max /= z_nuc) then
        print *, 'header read wrong: ', tab%z_nuc, tab%nd, tab%nt, tab%z1_min, tab%z1_max
        nfail = nfail + 1
    end if
    if (maxval(abs(tab%log_ne - log_ne)) > 1e-12_dp .or. maxval(abs(tab%log_te - log_te)) > 1e-12_dp) then
        print *, 'grids read wrong'
        nfail = nfail + 1
    end if

    ! grid points, interior points, outside the grid (clamped)
    do z1 = 1, z_nuc
        do k = 1, 7
            select case (k)
            case (1); ne = 10.0_dp**log_ne(1);  te = 10.0_dp**log_te(1)
            case (2); ne = 10.0_dp**log_ne(3);  te = 10.0_dp**log_te(4)
            case (3); ne = 10.0_dp**log_ne(nd); te = 10.0_dp**log_te(nt)
            case (4); ne = 3.7e12_dp;           te = 17.0_dp
            case (5); ne = 2.0e9_dp;            te = 0.45_dp
            case (6); ne = 1.0e5_dp;            te = 1.0e5_dp     ! below n_e grid, above T_e grid
            case (7); ne = 1.0e17_dp;           te = 0.01_dp      ! above n_e grid, below T_e grid
            end select
            got = adf11_coef(tab, z1, ne, te)
            want = 10.0_dp**f(z1, 0.0_dp, min(max(log10(ne), log_ne(1)), log_ne(nd)), &
                              min(max(log10(te), log_te(1)), log_te(nt)))
            call check('adf11_coef', z1, k, got, want)
        end do
    end do

    ! adas_rates: ionisation iz -> iz+1 = block Z1 = iz+1 (scd), recombination iz -> iz-1 = block
    ! Z1 = iz (acd), rate = n_e <sigma v>; no recombination to neutrals, no ionisation of Be4+
    ne = 4.2e13_dp; te = 33.0_dp
    rate_ion = -1.0_dp; rate_rec = -1.0_dp
    call adas_rates('Be', fscd, facd, ne, te, rate_ion, rate_rec)
    do k = 0, z_nuc
        want = 0.0_dp
        if (k < z_nuc) want = ne*10.0_dp**f(k + 1, 0.0_dp, log10(ne), log10(te))
        call check('rate_ion', k, 0, rate_ion(k), want)
        want = 0.0_dp
        if (k >= 2) want = ne*10.0_dp**f(k, 3.0_dp, log10(ne), log10(te))
        call check('rate_rec', k, 0, rate_rec(k), want)
    end do

    ! empty file names leave the rates unchanged
    rate_ion = 7.0_dp; rate_rec = 9.0_dp
    call adas_rates('Be', '', '', ne, te, rate_ion, rate_rec)
    if (any(rate_ion /= 7.0_dp) .or. any(rate_rec /= 9.0_dp)) then
        print *, 'adas_rates changed rates without files'
        nfail = nfail + 1
    end if

    if (nfail > 0) then
        print '(a,i0,a)', 'test_adf11: ', nfail, ' failures'
        error stop 1
    end if
    print '(a)', 'test_adf11: all checks passed'

contains

    pure function f(z1, shift, x, y) result(v)
        ! synthetic log10 <sigma v>, bilinear in x = log10 n_e, y = log10 T_e
        integer, intent(in) :: z1
        real(dp), intent(in) :: shift, x, y
        real(dp) :: v
        v = -10.0_dp - shift - 0.5_dp*z1 + 0.1_dp*z1*(x - 12.0_dp) + (1.5_dp - 0.2_dp*z1)*y &
            + 0.03_dp*(x - 12.0_dp)*y
    end function f

    subroutine write_file(file, shift)
        character(len=*), intent(in) :: file
        real(dp), intent(in) :: shift
        integer :: iu, z1, it, id
        open(newunit=iu, file=file, status='replace', action='write')
        write(iu, '(5i5,a)') z_nuc, nd, nt, 1, z_nuc, '     /BERYLLIUM          /SYNTHETIC TEST'
        write(iu, '(a)') repeat('-', 80)
        write(iu, '(8f10.5)') log_ne
        write(iu, '(8f10.5)') log_te
        do z1 = 1, z_nuc
            write(iu, '(a,i2,a)') '---------------------/ IGRD= 1  / IPRT= 1  /--------/ Z1=', z1, '   / DATE= 21/11/26'
            write(iu, '(8f10.5)') ((f(z1, shift, log_ne(id), log_te(it)), id = 1, nd), it = 1, nt)
        end do
        write(iu, '(a)') 'C-------------------------------------------------------------------------------'
        write(iu, '(a)') 'C  synthetic file for test_adf11'
        close(iu)
    end subroutine write_file

    subroutine check(what, i, k, got, want)
        character(len=*), intent(in) :: what
        integer, intent(in) :: i, k
        real(dp), intent(in) :: got, want
        real(dp) :: tol
        ! file values have 5 decimals in log10: relative error <= ln(10)*0.5e-5 ~ 1.2e-5
        tol = 2.0e-5_dp*abs(want)
        if (abs(got - want) > tol .or. (want == 0.0_dp .and. got /= 0.0_dp)) then
            print '(a,a,2i3,2es16.8)', 'FAIL ', what, i, k, got, want
            nfail = nfail + 1
        end if
    end subroutine check

end program test_adf11
