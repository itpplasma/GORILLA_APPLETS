program test
    implicit none

    integer :: b=5, i

    do i = -5,5
        print*, mod(i-1,b)+1
    enddo
    
end program test