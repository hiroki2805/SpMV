subroutine mesh_gen(nne,nelem,nnode,PETOT,my_rank)
  implicit none

  integer :: i, j, k, l, n

  integer, intent(in) :: PETOT, my_rank

  integer, intent(in) :: nne, nelem, nnode

  integer :: eleno(8)
  integer :: count, elem
  real(kind(0d0)) :: xnode, ynode, znode

  integer, allocatable :: eptr(:), eind(:)
  integer :: ijk

  character :: filename*20

  ! ----------------------------------------------------------------------------

  open(20,file="mesh.msh")

  write(20,'(I8)') nelem

  do k = 1, nne-1
    do j = 1, nne-1
      do i = 1, nne-1
        eleno(1) = i+nne*(j-1)+nne*nne*(k-1)
        eleno(2) = eleno(1) + 1
        eleno(3) = eleno(2) + nne
        eleno(4) = eleno(3) - 1
        eleno(5) = eleno(1) + nne*nne
        eleno(6) = eleno(5) + 1
        eleno(7) = eleno(6) + nne
        eleno(8) = eleno(7) - 1
        write(20,'(8(I8,1X))')  (eleno(l),l=1,8)
      end do
    end do
  end do

  close(20)

  ! ----------------------------------------------------------------------------

  ! ここにmesh data の arrayを作成 eind eptr必要
  ! eind : 16,387,064 (127*127*127*8)のはず

  allocate(eptr(nelem+1),eind(nelem*8))

  ! make eptr array
  eptr(1) = 0
  do i = 1, nelem
    eptr(i+1) = 8*i
  end do


  ! make eind array
  do k = 1, nne-1
    do j = 1, nne-1
      do i = 1, nne-1
        ijk = (i-1)*8 + (j-1)*(nne-1)*8 + (k-1)*(nne-1)*(nne-1)*8
        eind(ijk+1) = i+nne*(j-1)+nne*nne*(k-1)
        eind(ijk+2) = eind(ijk+1) + 1
        eind(ijk+3) = eind(ijk+2) + nne
        eind(ijk+4) = eind(ijk+3) - 1
        eind(ijk+5) = eind(ijk+1) + nne*nne
        eind(ijk+6) = eind(ijk+5) + 1
        eind(ijk+7) = eind(ijk+6) + nne
        eind(ijk+8) = eind(ijk+7) - 1
      end do
    end do
  end do



  deallocate(eptr,eind)

  ! ----------------------------------------------------------------------------

  n = PETOT

  write(filename,*) n

  filename='mesh.msh.epart.'//trim(adjustl(filename))

  ! write(*,*) filename

  ! ----------------------------------------------------------------------------

  ! サブルーチン呼び出しにしていきたい

  ! ----------------------------------------------------------------------------

  open(21,file='part.inp')

  write(21,'(5(I8,1X))') nnode, nelem, 0, 1, 0

  count = 0

  do k = 1, nne
    do j = 1, nne
      do i = 1, nne
        count = count + 1
        xnode = 0.0d0 + 1.0d-3*dble(i-1)
        ynode = 0.0d0 + 1.0d-3*dble(j-1)
        znode = 0.0d0 + 1.0d-3*dble(k-1)
        write(21,'((I8,1X),3(E17.8,1X))') count, xnode, ynode, znode
      end do
    end do
  end do


  count = 0

  do k = 1, nne-1
    do j = 1, nne-1
      do i = 1, nne-1
        count = count + 1
        eleno(1) = i+nne*(j-1)+nne*nne*(k-1)
        eleno(2) = eleno(1) + 1
        eleno(3) = eleno(2) + nne
        eleno(4) = eleno(3) - 1
        eleno(5) = eleno(1) + nne*nne
        eleno(6) = eleno(5) + 1
        eleno(7) = eleno(6) + nne
        eleno(8) = eleno(7) - 1
        write(21,'(2(I8,1X),(A5,1X),8(I8,1X))') count, 1, " hex", (eleno(l),l=1,8)
      end do
    end do
  end do

  write(21,'(2(I8,1X))') 1, 1
  write(21,'(A)') "ELEMENT, "

  ! ----------------------------------------------------------------------------

  open(30,file=filename)

  count = 0
  do i = 1, nelem
    count = count + 1
    read(30,*) elem
    write(21,*) count, elem
  end do

  close(21)
  close(30)

  ! ----------------------------------------------------------------------------

end subroutine mesh_gen
