program mesh
  implicit none

  integer :: i, j, k, l, n
  integer :: nne, nelem, nnode
  integer :: eleno(8)
  integer :: count, elem
  real(kind(0d0)) :: xnode, ynode, znode

  character :: filename*20

  write(*,*) "enter number 'nne' :: "
  read(*,*) nne

  nelem = (nne-1)**3
  nnode = nne**3

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

  write(*,*) "enter division number (exit:0) :: "
  read(*,*) n
  if ( n==0 ) then
    stop
  end if
  write(filename,*) n

  filename='mesh.msh.epart.'//trim(adjustl(filename))

  ! write(*,*) filename

  ! ----------------------------------------------------------------------------

  ! 今後はサブルーチン呼び出しにしていきたいな

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

end program mesh
