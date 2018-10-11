subroutine mesh_gen(nne,nelem,nnode,PETOT,my_rank,crs_col,crs_row,crs_count)
  implicit none

  integer :: i, j, k, l, m

  integer, intent(in) :: PETOT, my_rank

  integer, intent(in) :: nne, nelem, nnode

  integer :: eleno(8)
  integer :: count, elem
  real(kind(0d0)) :: xnode, ynode, znode

  ! for metis
  integer, allocatable :: eptr(:), eind(:)
  integer, allocatable :: epart(:), npart(:)
  integer :: objval
  integer, allocatable :: vwgt(:), vsize(:), tpwgts(:), options(:)
  integer :: ijk

  integer :: par_count
  integer, allocatable :: connect_par(:,:)

  integer, allocatable :: val_loc(:)

  integer, allocatable :: crs_col(:), crs_row(:)
  integer :: crs_count

  integer :: prev_no, prev_row

  character :: filename*20

  ! -------------------------- Mesh Gen ----------------------------------------

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

  ! -------------------------- Mesh Partitioning -------------------------------

  ! call Metis
  allocate(epart(nelem),npart(nnode))
  call METIS_PartMeshNodal(nelem,nnode,eptr,eind,vwgt,vsize,PETOT,tpwgts,options,objval,epart,npart)

  ! count number of element
  par_count = 0
  do i = 1, nelem
    if ( epart(i) == my_rank ) then
      par_count = par_count + 1
    end if
  end do

  ! make PartMesh
  allocate(connect_par(par_count,8))
  par_count = 0
  do i = 1, nelem
    if ( epart(i) == my_rank ) then
      par_count = par_count + 1
      do j = 1, 8
        connect_par(par_count,j) = eind((i-1)*8+j)
      end do
    end if
  end do

  deallocate(eptr,eind)
  deallocate(epart,npart)

  !! Make CRS
  ! 要素のある場所の配列
  allocate(val_loc(8*8*par_count))
  do i = 1, par_count
    do j = 1, 8
      do k = 1, 8
        ijk = (par_count-1)*i + (j-1)*8 + k
        val_loc(ijk) = (connect_par(i,j)-1)*nnode + connect_par(i,k)
      end do
    end do
  end do
  deallocate(connect_par)

  ! sort
  call heapsort(8*8*par_count,val_loc)

  ! Make CRS
  crs_count = 0
  prev_no = 0
  do i = 1, 8*8*par_count
    if ( val_loc(i) > prev_no ) then
      crs_count = crs_count + 1
      prev_no = val_loc(i)
    end if
  end do
  allocate(crs_col(crs_count),crs_row(nnode+1))

  crs_count = 0
  prev_no = 0
  prev_row = 1
  crs_row(1) = 1
  do i = 1, 8*8*par_count
    if ( val_loc(i) > prev_no ) then
      crs_count = crs_count + 1
      prev_no = val_loc(i)
      crs_col(crs_count) = mod(val_loc(i),nnode)
      if ( val_loc(i) > nnode*prev_row ) then
        crs_row(prev_row + 1) = crs_count
        prev_row = prev_row + 1
      end if
    end if
  end do
  crs_row(nnode) = crs_count





  ! reset element number
  !prev_no = 0
  !do i = 1, nnode
  !  elemdo : do j = 1, par_count
  !    do k = 1, 8
  !      if ( (i == prev_no + 1) .and. (connect_par(j,k) == i) ) then
  !        prev_no = i
  !        exit elemdo
  !      end if
  !      if (( i > prev_no + 1 ) .and. (connect_par(j,k) == i) )then
  !        do l = 1, par_count
  !          do m = 1, 8
  !            if ( connect_par(l,m) == i ) then
  !              connect_par(l,m) = prev_no + 1
  !            end if
  !          end do
  !        end do
  !        prev_no = prev_no + 1
  !      end if
  !    end do
  !  end do elemdo

  !end do

  ! ----------------------------------------------------------------------------

  !write(filename,*) PETOT

  !filename='mesh.msh.epart.'//trim(adjustl(filename))

  ! write(*,*) filename

  ! ----------------------------------------------------------------------------

  ! ----------------------------------------------------------------------------

  !open(21,file='part.inp')

  !write(21,'(5(I8,1X))') nnode, nelem, 0, 1, 0

  !count = 0

  !do k = 1, nne
  !  do j = 1, nne
  !    do i = 1, nne
  !      count = count + 1
  !      xnode = 0.0d0 + 1.0d-3*dble(i-1)
  !      ynode = 0.0d0 + 1.0d-3*dble(j-1)
  !      znode = 0.0d0 + 1.0d-3*dble(k-1)
  !      write(21,'((I8,1X),3(E17.8,1X))') count, xnode, ynode, znode
  !    end do
  !  end do
  !end do


  !count = 0

  !do k = 1, nne-1
  !  do j = 1, nne-1
  !    do i = 1, nne-1
  !      count = count + 1
  !      eleno(1) = i+nne*(j-1)+nne*nne*(k-1)
  !      eleno(2) = eleno(1) + 1
  !      eleno(3) = eleno(2) + nne
  !      eleno(4) = eleno(3) - 1
  !      eleno(5) = eleno(1) + nne*nne
  !      eleno(6) = eleno(5) + 1
  !      eleno(7) = eleno(6) + nne
  !      eleno(8) = eleno(7) - 1
  !      write(21,'(2(I8,1X),(A5,1X),8(I8,1X))') count, 1, " hex", (eleno(l),l=1,8)
  !    end do
  !  end do
  !end do

  !write(21,'(2(I8,1X))') 1, 1
  !write(21,'(A)') "ELEMENT, "

  ! ----------------------------------------------------------------------------

  !open(30,file=filename)

  !count = 0
  !do i = 1, nelem
  !  count = count + 1
  !  read(30,*) elem
  !  write(21,*) count, elem
  !end do

  !close(21)
  !close(30)

  ! ----------------------------------------------------------------------------

end subroutine mesh_gen
