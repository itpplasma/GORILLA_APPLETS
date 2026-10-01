module rmp_gc_measure_mod
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private
    public :: gc_measure_correction, gc_log_measure_change
contains
    pure function gc_measure_correction(h, symplectic_curl, metric, B) result(chi)
        real(dp), intent(in) :: h(3), symplectic_curl(3), metric, B
        real(dp) :: chi

        if (metric <= 0.0_dp .or. B <= 0.0_dp) error stop 'Invalid GC metric or B'
        chi = dot_product(h, symplectic_curl)/(metric*B)
        if (chi <= 0.0_dp) error stop 'Nonpositive GC phase-space measure'
    end function gc_measure_correction

    pure function gc_log_measure_change(chi_old, chi_new) result(dlog)
        real(dp), intent(in) :: chi_old, chi_new
        real(dp) :: dlog

        if (chi_old <= 0.0_dp .or. chi_new <= 0.0_dp) &
            error stop 'Nonpositive collision phase-space measure'
        dlog = log(chi_new) - log(chi_old)
    end function gc_log_measure_change
end module rmp_gc_measure_mod
