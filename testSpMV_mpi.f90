!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!/
!
! SpMV calculation sample code.
!
! This code run Sparse Matrix Vector product and measure the calculation time.
! The matrix data structure mimic the stiffness matrix used in Finite Element Method.
! Calculation model is a cube, divided by hexa linear elements.
! Each edge is divided as nnodes_edge nodes.
!
! Usage: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Just compile and execute.
!
! icc -fast testSpMV_block11.c && ./a.out
!
!
! Input: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Nothing. Calculation parameters are hard coded.
!
!
! Parameter variables: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! nnodes_edge: This variable define the size of problem.
!              nnodes_edge^3 is the total dimension of sparse matrix.
!              (number of non-zero elements on the matrix is given by mimic stiffness matrix)
!
!
! itermax:     This variable define how many times SpMV is calculated.
!
!
! Output: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!/
!
! Followings are outputted in stdout.
!
!                                                               ! comment
! nne   128                                                     ! number of nodes per edge
! nne2  16384                                                   ! nne^2
! nrows 2097152                                                 ! number of rows of sparse matrix (nne^3)
! ncols 2097152                                                 ! number of columns of sparse matrix (nne^3)
! nnz 55742968                                                  ! number of non-zero elements in sparse matrix
! CLOCKS_PER_SEC 1000000                                        ! constant for clock function
! user time for count number of non-zero element (sec) 0.050000
! user time for create irow jcol matrix (sec) 0.200000
! now making CRS format
! user time for create CRS matrix (sec) 0.310000
! user time for create vector B (sec) 0.000000
! now doing SpMV calculation 50 times
! user time for SpMV (sec) 4.170000                             ! This line gives SpMV calculation time in sec.
! test: now write y vector
! 4798.080000                                      ! y vector of Ab=y calculation (50th iteration) is shown.
! 7197.120000                                      ! It prevent compiler optimization of iter loop and
! 7197.120000                                      ! useful to check calculation result.
! 7197.120000
! 7197.120000
! ..
subroutine testSpMV_block11_mpi(itermax,nne,nelem,nnode,mev,vev,PETOT,my_rank,crs_col,crs_row,crs_count)
  implicit none

  integer :: x1, y1, z1
  integer :: x2, y2, z2

  integer, intent(in) :: itermax, nne, nelem, nnode

  integer, intent(in) :: PETOT, my_rank

  real(kind(0d0)), intent(in) :: mev, vev

  integer, intent(in) :: crs_col(crs_count), crs_row(nnode+1), crs_count

  integer :: i, j, k
  integer :: iter








end subroutine testSpMV_block11_mpi
