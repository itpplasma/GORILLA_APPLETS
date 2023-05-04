!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine random_num(ur)
    ur=zzg()
    return
    END
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine getran(irand,ur)
!
!  Produces the random number with zero average and unit variance
!
!  Input parameters: irand - 0 for continious, 1 for +1 -1,
!  Output parameters: ur   - random number
!
    call random_num(ur)
!
    if(irand.eq.0) then
!
!  continiuos random number, constant is sqrt(12)
!
      ur=3.464102*(ur-.5)
    else
!
!  discrete random number
!
      if(ur.gt..5) then
        ur=1.
      else
        ur=-1.
      endif
    endif
    return
    end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc