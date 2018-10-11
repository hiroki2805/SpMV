program spmv_main
  implicit none

  ! --- for MPI ---
  include 'mpif.h'
  integer :: PETOT, my_rank, ierr
  ! ---------------

  integer :: itermax
  integer :: nnodes_edge, nne, nelem, nnode

  real(kind(0d0)) :: mev, vev

  ! ---------------------------- Run MPI -------------------------------

  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,PETOT,ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)

  ! --------------------------------------------------------------------



  ! ----------------------------- Setting -----------------------------

  ! SpMV iteration times
  itermax = 50

  ! number of nodes on a edge of cube
  nnodes_edge = 128

  ! value of matrix element
  mev = 1.2

  ! value of vector element
  vev = 10.2

  ! -----------------------------------------------------------------

  nne  = nnodes_edge
  nelem = (nne-1)**3
  nnode = nne**3

  ! -----------------------------------------------------------------

  ! ----------------------- Mesh Generator --------------------------

  call mesh_gen(nne,nelem,nnode,PETOT,my_rank)

  ! -----------------------------------------------------------------

  ! ------------------------ SpMV (MPI) -----------------------------

  call testSpMV_block11_mpi(itermax,nne,nelem,nnode,mev,vev,PETOT,my_rank,crs_col,crs_row,crs_count)

  ! -----------------------------------------------------------------

  ! ------------------------- Close MPI -----------------------------

  call MPI_FINALIZE(ierr)

  ! -----------------------------------------------------------------


end program spmv_main
