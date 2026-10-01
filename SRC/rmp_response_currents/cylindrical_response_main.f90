program gorilla_cylindrical_response
    use, intrinsic :: iso_fortran_env, only: real64
    use cylindrical_electron_response_mod, only: assemble_cylindrical_electrons
    implicit none
    character(1024) :: input_path,output_path,header,arg
    integer :: unit,M,N,model,mmode,nmode,refinement=1,j,k,row,col
    real(real64) :: L,rm,clight
    real(real64),allocatable :: background(:,:)
    complex(real64),allocatable :: blocks(:,:,:)
    call get_command_argument(1,input_path)
    call get_command_argument(2,output_path)
    call get_command_argument(3,arg)
    if (len_trim(input_path)==0 .or. len_trim(output_path)==0) &
        error stop 'usage: gorilla_cylindrical_response.x background.dat response.dat [refinement]'
    if (len_trim(arg)>0) read(arg,*) refinement
    open(newunit=unit,file=trim(input_path),status='old',action='read')
    read(unit,'(a)') header
    if (trim(header)/='GK_BACKGROUND_V1') error stop 'wrong background format'
    read(unit,*) M,N,model,mmode,nmode
    read(unit,*) L,rm,clight
    if (model/=0) error stop 'cylindrical provider implements number-conserving OU model 0 only'
    if (N<1 .or. refinement<1) error stop 'invalid grid or quadrature refinement'
    allocate(background(N,13))
    do j=1,N
        read(unit,*) background(j,:)
    end do
    close(unit)
    call assemble_cylindrical_electrons(background,M,L,clight,blocks,refinement)
    open(newunit=unit,file=trim(output_path),status='new',action='write')
    write(unit,'(a)') 'GK_RESPONSE_V1'
    write(unit,*) M,N,model,mmode,nmode
    write(unit,'(3es26.17e3)') L,rm,clight
    do j=1,N
        write(unit,'(13es26.17e3)') background(j,:)
    end do
    do k=1,4
        do col=1,2*M+1
            do row=1,2*M+1
                write(unit,'(2es26.17e3)') real(blocks(row,col,k)),aimag(blocks(row,col,k))
            end do
        end do
    end do
    close(unit)
    print *, 'Wrote drift-kinetic cylindrical electron response; OU model 0, refinement=',refinement
end program
