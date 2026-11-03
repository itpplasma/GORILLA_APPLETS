!> Metropolis-Hastings chain built from ntv_lorentz_proposal (canonical NTV
!> Lorentz kicks): with target pi(lambda) ~ exp(a lambda) (the lambda dependence
!> of f0 through p_phi at fixed position and speed) the chain must sample pi;
!> with a = 0 it must keep the isotropic distribution (<lambda> = 0,
!> <lambda^2> = 1/3). Also checks the first two moments of one proposal.
program test_ntv_lorentz_proposal

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use canonical_ntv_mod, only: ntv_lorentz_proposal

    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    integer :: nfail, ns, k
    integer, allocatable :: seed(:)

    call random_seed(size=ns)
    allocate(seed(ns))
    seed = [(12345 + 17 * k, k = 1, ns)]
    call random_seed(put=seed)
    nfail = 0
    call check_chain(0.0_dp, nfail)
    call check_chain(0.5_dp, nfail)
    call check_single_step(nfail)
    if (nfail > 0) then
        print '(a, i0, a)', 'FAIL: ', nfail, ' Lorentz proposal checks failed'
        error stop 1
    end if
    print '(a)', 'PASS: Lorentz MH proposal'

contains

    function normal() result(xi)
        real(dp) :: xi, u(2)
        call random_number(u)
        u(1) = max(u(1), tiny(1.0_dp))
        xi = sqrt(-2.0_dp * log(u(1))) * cos(2.0_dp * pi * u(2))
    end function normal

    subroutine check_chain(a, nfail)
        real(dp), intent(in) :: a
        integer, intent(inout) :: nfail
        integer, parameter :: nchain = 400, nstep = 20000, nburn = 2000
        real(dp), parameter :: nu_dt = 0.02_dp
        real(dp) :: lam, lam_new, logq, u, s1, s2, cnt, m1, m2, e1, e2
        logical :: valid
        integer :: ic, k

        s1 = 0.0_dp
        s2 = 0.0_dp
        cnt = 0.0_dp
        do ic = 1, nchain
            call random_number(u)
            lam = 2.0_dp * u - 1.0_dp
            do k = 1, nstep
                call ntv_lorentz_proposal(lam, nu_dt, normal(), lam_new, logq, valid)
                if (valid) then
                    call random_number(u)
                    if (log(max(u, tiny(1.0_dp))) < a * (lam_new - lam) + logq) &
                        lam = lam_new
                end if
                if (k > nburn .and. mod(k, 10) == 0) then
                    s1 = s1 + lam
                    s2 = s2 + lam**2
                    cnt = cnt + 1.0_dp
                end if
            end do
        end do
        m1 = s1 / cnt
        m2 = s2 / cnt
        if (a == 0.0_dp) then
            e1 = 0.0_dp
            e2 = 1.0_dp / 3.0_dp
        else
            e1 = 1.0_dp / tanh(a) - 1.0_dp / a
            e2 = 1.0_dp - 2.0_dp * e1 / a
        end if
        print '(a, f5.2, a, 2f10.5, a, 2f10.5)', ' a = ', a, ': <lam>, <lam^2> =', &
            m1, m2, '  exact', e1, e2
        if (abs(m1 - e1) > 6.0e-3_dp .or. abs(m2 - e2) > 6.0e-3_dp) nfail = nfail + 1
    end subroutine check_chain

    subroutine check_single_step(nfail)
        integer, intent(inout) :: nfail
        integer, parameter :: n = 400000
        real(dp), parameter :: lam0 = 0.4_dp, nu_dt = 1.0e-3_dp
        real(dp) :: lam_new, logq, s1, s2, mean, var
        logical :: valid
        integer :: k

        s1 = 0.0_dp
        s2 = 0.0_dp
        do k = 1, n
            call ntv_lorentz_proposal(lam0, nu_dt, normal(), lam_new, logq, valid)
            s1 = s1 + (lam_new - lam0)
            s2 = s2 + (lam_new - lam0)**2
        end do
        mean = s1 / n
        var = s2 / n - mean**2
        print '(a, 2es12.4, a, 2es12.4)', ' one step: mean, var =', mean, var, &
            '  expected', -lam0 * nu_dt, (1.0_dp - lam0**2) * nu_dt
        if (abs(mean + lam0 * nu_dt) > 5.0_dp * sqrt(var / n)) nfail = nfail + 1
        if (abs(var / ((1.0_dp - lam0**2) * nu_dt) - 1.0_dp) > 1.0e-2_dp) &
            nfail = nfail + 1
    end subroutine check_single_step

end program test_ntv_lorentz_proposal
