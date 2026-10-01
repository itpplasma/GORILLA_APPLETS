module toroidal_profile_slices_mod
    implicit none
    private
    public :: toroidal_profile_bounds
contains
    pure subroutine toroidal_profile_bounds(ntetr, nphi, iphi, first, last)
        integer, intent(in) :: ntetr, nphi, iphi
        integer, intent(out) :: first, last
        integer :: per_slice

        ! The mesh stores equally sized toroidal slices. Divide first so that
        ! every intermediate index remains bounded by the total mesh size.
        per_slice = ntetr/nphi
        first = iphi*per_slice + 1
        last = (iphi + 1)*per_slice
    end subroutine toroidal_profile_bounds
end module toroidal_profile_slices_mod
