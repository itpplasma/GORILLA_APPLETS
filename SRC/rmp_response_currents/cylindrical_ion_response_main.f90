program gorilla_cylindrical_ion_response
    use, intrinsic :: iso_fortran_env, only: real64
    use cylindrical_ion_response_mod, only: assemble_cylindrical_ions
    implicit none
    character(1024) :: input_path,output_path,header,arg
    integer :: unit,M,N,model,mmode,nmode,species,refinement=1,n_perp=64,j,k,row,col,ios
    real(real64) :: L,rm,clight
    real(real64),allocatable :: background(:,:)
    complex(real64),allocatable :: blocks(:,:,:)
    call get_command_argument(1,input_path)
    call get_command_argument(2,output_path)
    if(len_trim(input_path)==0.or.len_trim(output_path)==0) &
        error stop 'usage: gorilla_cylindrical_ion_response.x input output [n_perp] [refinement]'
    call get_command_argument(3,arg)
    if(len_trim(arg)>0) then
        read(arg,*,iostat=ios) n_perp
        if(ios/=0) error stop 'invalid n_perp'
    end if
    call get_command_argument(4,arg)
    if(len_trim(arg)>0) then
        read(arg,*,iostat=ios) refinement
        if(ios/=0) error stop 'invalid refinement'
    end if
    open(newunit=unit,file=trim(input_path),status='old',action='read')
    read(unit,'(a)') header
    if(trim(header)/='GK_ION_BACKGROUND_V1') error stop 'wrong ion background format'
    read(unit,*) M,N,model,mmode,nmode,species
    read(unit,*) L,rm,clight
    if(N<1.or.species<1) error stop 'invalid ion grid/species'
    allocate(background(N,14))
    do j=1,N
        read(unit,*) background(j,:)
    end do
    close(unit)
    call assemble_cylindrical_ions(background,M,L,clight,blocks,refinement,model,n_perp)
    open(newunit=unit,file=trim(output_path),status='new',action='write')
    write(unit,'(a)') 'GK_ION_RESPONSE_V1'
    write(unit,*) M,N,model,mmode,nmode,species
    write(unit,'(3es26.17e3)') L,rm,clight
    do j=1,N
        write(unit,'(14es26.17e3)') background(j,:)
    end do
    do k=1,4
        do col=1,2*M+1
            do row=1,2*M+1
                write(unit,'(2es26.17e3)') real(blocks(row,col,k)),aimag(blocks(row,col,k))
            end do
        end do
    end do
    close(unit)
    print *, 'Wrote integral-FLR ion response; species=',species,' model=',model,' n_perp=',n_perp
end program
