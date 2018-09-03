program mesh
  implicit none

  integer :: i, j, k, l
  integer :: nne, nelem, nnode
  integer :: eleno(8)

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

end program mesh
