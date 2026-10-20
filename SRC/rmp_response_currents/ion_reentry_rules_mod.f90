module ion_reentry_rules_mod
    ! Acceptance rule for opt-in periodic radial re-entry. The pusher reports every
    ! removal as a loss (ind_tetr=-1): a real exit through the radial mesh boundary,
    ! but also troubleshooting/no-root failures, the stagnant-handover guard and the
    ! time-reversal guard. Only a real boundary exit may be re-entered:
    !  - the exit point lies on a boundary face (no neighbour) of the last cell, within
    !    the face roundoff bound (a negative bound flags "no boundary face found");
    !  - the remaining time has the sign of the step and does not exceed it.
    ! Anything else must fall through to the caller's abort with the saved state.
    use, intrinsic :: iso_fortran_env, only: dp=>real64
    implicit none
    private
    public :: reentry_exit_accepted
contains
    pure logical function reentry_exit_accepted(face_distance,bound,remaining,step) result(accepted)
        real(dp),intent(in) :: face_distance,bound,remaining,step
        accepted=bound>=0.and.abs(face_distance)<=bound.and.remaining*step>=0.and.abs(remaining)<=abs(step)
    end function
end module ion_reentry_rules_mod
