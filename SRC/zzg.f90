module zzg_mod

  implicit none

  contains

function zzg()
    real :: zzg,xi
    call random_number(xi)
    zzg=xi
  end function zzg

end module zzg_mod