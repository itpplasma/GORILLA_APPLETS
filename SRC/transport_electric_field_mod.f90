!Electrostatic potential normalisation for the mono-energetic transport applet.
!The transport applet prescribes the normalised ExB velocity
!    v_E = c*E_r/(B_00*v)
!at the flux tube, with E_r = dPhi/dr and r = sqrt(s*|psi_tor,a|/(pi*B_00)).
!GORILLA uses Phi = eps_Phi*A_theta(s), so eps_Phi follows from dA_theta/ds.
module transport_electric_field_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none

    private

    public :: eps_phi_from_v_e, chartmap_bmod00

contains

    pure function eps_phi_from_v_e(v_e, vmod, bmod00, da_theta_ds, ds_dr) result(eps_phi)

        use constants, only: clight

        real(dp), intent(in) :: v_e, vmod, bmod00, da_theta_ds, ds_dr
        real(dp) :: eps_phi

        eps_phi = v_e*vmod*bmod00/(da_theta_ds*ds_dr*clight)

    end function eps_phi_from_v_e

    !B_00 on flux surface s of a loaded Boozer chartmap: in Boozer angles the
    !(0,0) harmonic of |B| is its mean over a uniform (theta, zeta) grid.
    function chartmap_bmod00(s) result(bmod00)

        use constants, only: pi
        use boozer_chartmap_mod, only: evaluate_boozer_chartmap, chartmap_nfp

        real(dp), intent(in) :: s
        real(dp) :: bmod00

        integer, parameter :: ntheta = 64, nzeta = 64
        integer :: i, j
        real(dp) :: theta, zeta, a_s, a_theta, a_phi, h_s, h_theta, h_phi, bmod, sqrtg

        bmod00 = 0.0_dp
        do j = 1, nzeta
            zeta = 2.0_dp*pi/real(chartmap_nfp, dp)*real(j - 1, dp)/real(nzeta, dp)
            do i = 1, ntheta
                theta = 2.0_dp*pi*real(i - 1, dp)/real(ntheta, dp)
                call evaluate_boozer_chartmap(s, theta, zeta, a_s, a_theta, a_phi, &
                                              h_s, h_theta, h_phi, bmod, sqrtg)
                bmod00 = bmod00 + bmod
            end do
        end do
        bmod00 = bmod00/real(ntheta*nzeta, dp)

    end function chartmap_bmod00

end module transport_electric_field_mod
