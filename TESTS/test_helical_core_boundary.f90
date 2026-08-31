program test_helical_core_boundary

    use helical_core_boundary_mod, only: boundary_absorb, boundary_action, &
        boundary_reflect

    implicit none

    if (boundary_action(.true.) /= boundary_reflect) &
        error stop 'inner computational hole must reflect markers'
    if (boundary_action(.false.) /= boundary_absorb) &
        error stop 'outer SOL boundary must absorb markers'

    print *, 'PASS: helical-core boundary contract'

end program test_helical_core_boundary
