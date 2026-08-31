module helical_core_boundary_mod

    implicit none

    private
    public :: boundary_action
    public :: boundary_absorb, boundary_reflect

    integer, parameter :: boundary_reflect = 1
    integer, parameter :: boundary_absorb = 2

contains

pure integer function boundary_action(crossed_inner_hole)

    logical, intent(in) :: crossed_inner_hole

    if (crossed_inner_hole) then
        boundary_action = boundary_reflect
    else
        boundary_action = boundary_absorb
    end if

end function boundary_action

end module helical_core_boundary_mod
