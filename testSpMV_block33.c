/////////////////////////////////////////////////////////////////////////////////
//
// SpMV calculation sample code. 
//
// This code run Sparse Matrix Vector product and measure the calculation time. 
// The matrix data structure mimic the stiffness matrix used in Finite Element Method. 
// Calculation model is a cube, divided by hexa linear elements. 
// Each edge is divided as nnodes_edge nodes. 
// Each non-zero matrix element have 3x3 blocked value.
// 
// Usage: //////////////////////////////////////////////////////////////////
// Just compile and execute. 
// 
// icc -fast testSpMV_block33.c && ./a.out
//
//
// Input: //////////////////////////////////////////////////////////////////
// Nothing. Calculation parameters are hard coded. 
// 
//
// Parameter variables: //////////////////////////////////////////////////////////////////
//
// nnodes_edge: This variable define the size of problem. 
//              nnodes_edge^3 is the total dimension of sparse matrix.
//              (number of non-zero elements on the matrix is given by mimic stiffness matrix)
//
//
// itermax:     This variable define how many times SpMV is calculated. 
// 
//
// Output: /////////////////////////////////////////////////////////////////////////////
//
// Followings are outputted in stdout.
// 
//                                                               // comment
// nne   128                                                     // number of nodes per edge
// nne2  16384                                                   // nne^2
// nrows 2097152                                                 // number of rows of sparse matrix (nne^3)
// ncols 2097152                                                 // number of columns of sparse matrix (nne^3)
// nnz 55742968                                                  // number of non-zero elements in sparse matrix
// CLOCKS_PER_SEC 1000000                                        // constant for clock function
// user time for count number of non-zero element (sec) 0.050000 
// user time for create irow jcol matrix (sec) 0.200000 
// now making CRS format
// user time for create CRS matrix (sec) 0.310000 
// user time for create vector B (sec) 0.000000 
// now doing SpMV calculation 50 times
// user time for SpMV (sec) 4.170000                             // This line gives SpMV calculation time in sec.
// test: now write y vector
// 14394.240000                                      // y vector of Ab=y calculation (50th iteration) is shown.
// 14394.240000                                      // It prevent compiler optimization of iter loop and 
// 14394.240000                                      // useful to check calculation result.
// 21591.360000
// 21591.360000
// ..
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// #include "omp.h"


main()
{
  int x1, y1, z1;
  int x2, y2, z2;

  int itermax;
  int nnodes_edge;

  double mev;
  double vev;

  int i,j,k;
  int iter;


// SpMV iteration times
  itermax = 50; // numbr of SpMV

// number of nodes on a edge of cube
  nnodes_edge = 128;

// value of matrix element
  mev = 1.2;

// value of vector element
  vev = 10.2;

// irow jcol matrix
  int nrows, ncols, nnz;
  int matdim;
  int nne, nne2;
  int *irow, *jcol;

// csr matrix
  int *csr_idx, *csr_item;
  int jsta, jend;
  double *csr_val;
  int row;

  int coef_block=9; //BCSR parameter : number of element in a block (fix to 9)
  int vec_coef_block=3; //BCSR parameter : number of element in a block (fix to 3)

// vector
  double *b_val;
  double *y;

  double b1, b2, b3, yv1, yv2, yv3;

  int inod;

// time
  clock_t t1, t2;
  double sec;


/////////////////////////////////////////////////////////
//
// create matrix A
// 

nne  = nnodes_edge;
nne2 = nnodes_edge * nnodes_edge;

matdim = nnodes_edge * nnodes_edge * nnodes_edge;

nrows = matdim;
ncols = matdim;

printf("nne   %d\n",nne   ); //DEBUG
printf("nne2  %d\n",nne2  ); //DEBUG
printf("nrows %d\n",nrows ); //DEBUG
printf("ncols %d\n",ncols ); //DEBUG

// count number of non zero elements
nnz=0; 

t1 = clock(); // TIMER

//printf("t1 %ld \n",t1);// DEBUG

for (z1=1; z1<=nne; z1++){
  for (y1=1; y1<=nne; y1++){
    for (x1=1; x1<=nne; x1++){
      i = nne2*(z1 - 1) + nne*(y1 - 1) + x1;
//printf("x1 %d y1 %d z1 %d i %d\n",x1,y1,z1,i ); //DEBUG

//////////////////////////////////////////////////////////////////// z=-1
x2=x1-1;
y2=y1-1;
z2=z1-1;
j=nne2*(z2-1)+nne*(y2-1)+x2;
if ((x2 >= 1) && (y2 >= 1) && (z2 >= 1) && (x2 <= nne) && (y2 <= nne) && (z2 <= nne)) {
  nnz++;
//printf("nnz++ x2 %d y2 %d z2 %d i %d\n",x2,y2,z2,i ); //DEBUG
}
//////////////////////////////
x2=x1;
y2=y1-1;
z2=z1-1;
j=nne2*(z2-1)+nne*(y2-1)+x2;
if ((x2 >= 1) && (y2 >= 1) && (z2 >= 1) && (x2 <= nne) && (y2 <= nne) && (z2 <= nne)) {
  nnz = nnz+1;
}
//////////////////////////////
x2=x1+1;
y2=y1-1;
z2=z1-1;
j=nne2*(z2-1)+nne*(y2-1)+x2;
if ((x2 >= 1) && (y2 >= 1) && (z2 >= 1) && (x2 <= nne) && (y2 <= nne) && (z2 <= nne)) {
  nnz = nnz+1;
}
////////////////////////////////////////////////////////////
x2=x1-1;
y2=y1;
z2=z1-1;
j=nne2*(z2-1)+nne*(y2-1)+x2;
if ((x2 >= 1) && (y2 >= 1) && (z2 >= 1) && (x2 <= nne) && (y2 <= nne) && (z2 <= nne)) {
  nnz = nnz+1;
}
//////////////////////////////
x2=x1;
y2=y1;
z2=z1-1;
j=nne2*(z2-1)+nne*(y2-1)+x2;
if ((x2 >= 1) && (y2 >= 1) && (z2 >= 1) && (x2 <= nne) && (y2 <= nne) && (z2 <= nne)) {
  nnz = nnz+1;
}
//////////////////////////////
x2=x1+1;
y2=y1;
z2=z1-1;
j=nne2*(z2-1)+nne*(y2-1)+x2;
if ((x2 >= 1) && (y2 >= 1) && (z2 >= 1) && (x2 <= nne) && (y2 <= nne) && (z2 <= nne)) {
  nnz = nnz+1;
}
////////////////////////////////////////////////////////////
x2=x1-1;
y2=y1+1;
z2=z1-1;
j=nne2*(z2-1)+nne*(y2-1)+x2;
if ((x2 >= 1) && (y2 >= 1) && (z2 >= 1) && (x2 <= nne) && (y2 <= nne) && (z2 <= nne)) {
  nnz = nnz+1;
}
//////////////////////////////
x2=x1;
y2=y1+1;
z2=z1-1;;
j=nne2*(z2-1)+nne*(y2-1)+x2;
if ((x2 >= 1) && (y2 >= 1) && (z2 >= 1) && (x2 <= nne) && (y2 <= nne) && (z2 <= nne)) {
  nnz = nnz+1;
}
//////////////////////////////
x2=x1+1;
y2=y1+1;
z2=z1-1;
j=nne2*(z2-1)+nne*(y2-1)+x2;
if ((x2 >= 1) && (y2 >= 1) && (z2 >= 1) && (x2 <= nne) && (y2 <= nne) && (z2 <= nne)) {
  nnz = nnz+1;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// z=0
x2=x1-1;
y2=y1-1;
z2=z1;
j=nne2*(z2-1)+nne*(y2-1)+x2;
if ((x2 >= 1) && (y2 >= 1) && (z2 >= 1) && (x2 <= nne) && (y2 <= nne) && (z2 <= nne)) {
  nnz = nnz+1;
}
//////////////////////////////
x2=x1;
y2=y1-1;
z2=z1;
j=nne2*(z2-1)+nne*(y2-1)+x2;
if ((x2 >= 1) && (y2 >= 1) && (z2 >= 1) && (x2 <= nne) && (y2 <= nne) && (z2 <= nne)) {
  nnz = nnz+1;
}
//////////////////////////////
x2=x1+1;
y2=y1-1;
z2=z1;
j=nne2*(z2-1)+nne*(y2-1)+x2;
if ((x2 >= 1) && (y2 >= 1) && (z2 >= 1) && (x2 <= nne) && (y2 <= nne) && (z2 <= nne)) {
  nnz = nnz+1;
}
////////////////////////////////////////////////////////////
x2=x1-1;
y2=y1;
z2=z1;
j=nne2*(z2-1)+nne*(y2-1)+x2;
if ((x2 >= 1) && (y2 >= 1) && (z2 >= 1) && (x2 <= nne) && (y2 <= nne) && (z2 <= nne)) {
  nnz = nnz+1;
}
//////////////////////////////
x2=x1;
y2=y1;
z2=z1;
j=nne2*(z2-1)+nne*(y2-1)+x2;
if ((x2 >= 1) && (y2 >= 1) && (z2 >= 1) && (x2 <= nne) && (y2 <= nne) && (z2 <= nne)) {
  nnz = nnz+1;
}
//////////////////////////////
x2=x1+1;
y2=y1;
z2=z1;
j=nne2*(z2-1)+nne*(y2-1)+x2;
if ((x2 >= 1) && (y2 >= 1) && (z2 >= 1) && (x2 <= nne) && (y2 <= nne) && (z2 <= nne)) {
  nnz = nnz+1;
}
////////////////////////////////////////////////////////////
x2=x1-1;
y2=y1+1;
z2=z1;
j=nne2*(z2-1)+nne*(y2-1)+x2;
if ((x2 >= 1) && (y2 >= 1) && (z2 >= 1) && (x2 <= nne) && (y2 <= nne) && (z2 <= nne)) {
  nnz = nnz+1;
}
//////////////////////////////
x2=x1;
y2=y1+1;
z2=z1;
j=nne2*(z2-1)+nne*(y2-1)+x2;
if ((x2 >= 1) && (y2 >= 1) && (z2 >= 1) && (x2 <= nne) && (y2 <= nne) && (z2 <= nne)) {
  nnz = nnz+1;
}
//////////////////////////////
x2=x1+1;
y2=y1+1;
z2=z1;
j=nne2*(z2-1)+nne*(y2-1)+x2;
if ((x2 >= 1) && (y2 >= 1) && (z2 >= 1) && (x2 <= nne) && (y2 <= nne) && (z2 <= nne)) {
  nnz = nnz+1;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// z=+1;
x2=x1-1;
y2=y1-1;
z2=z1+1;
j=nne2*(z2-1)+nne*(y2-1)+x2;
if ((x2 >= 1) && (y2 >= 1) && (z2 >= 1) && (x2 <= nne) && (y2 <= nne) && (z2 <= nne)) {
  nnz = nnz+1;
}
//////////////////////////////
x2=x1;
y2=y1-1;
z2=z1+1;
j=nne2*(z2-1)+nne*(y2-1)+x2;
if ((x2 >= 1) && (y2 >= 1) && (z2 >= 1) && (x2 <= nne) && (y2 <= nne) && (z2 <= nne)) {
  nnz = nnz+1;
}
//////////////////////////////
x2=x1+1;
y2=y1-1;
z2=z1+1;
j=nne2*(z2-1)+nne*(y2-1)+x2;
if ((x2 >= 1) && (y2 >= 1) && (z2 >= 1) && (x2 <= nne) && (y2 <= nne) && (z2 <= nne)) {
  nnz = nnz+1;
}
////////////////////////////////////////////////////////////
x2=x1-1;
y2=y1;
z2=z1+1;
j=nne2*(z2-1)+nne*(y2-1)+x2;
if ((x2 >= 1) && (y2 >= 1) && (z2 >= 1) && (x2 <= nne) && (y2 <= nne) && (z2 <= nne)) {
  nnz = nnz+1;
}
//////////////////////////////
x2=x1;
y2=y1;
z2=z1+1;
j=nne2*(z2-1)+nne*(y2-1)+x2;
if ((x2 >= 1) && (y2 >= 1) && (z2 >= 1) && (x2 <= nne) && (y2 <= nne) && (z2 <= nne)) {
  nnz = nnz+1;
}
//////////////////////////////
x2=x1+1;
y2=y1;
z2=z1+1;
j=nne2*(z2-1)+nne*(y2-1)+x2;
if ((x2 >= 1) && (y2 >= 1) && (z2 >= 1) && (x2 <= nne) && (y2 <= nne) && (z2 <= nne)) {
  nnz = nnz+1;
}
////////////////////////////////////////////////////////////
x2=x1-1;
y2=y1+1;
z2=z1+1;
j=nne2*(z2-1)+nne*(y2-1)+x2;
if ((x2 >= 1) && (y2 >= 1) && (z2 >= 1) && (x2 <= nne) && (y2 <= nne) && (z2 <= nne)) {
  nnz = nnz+1;
}
//////////////////////////////
x2=x1;
y2=y1+1;
z2=z1+1;
j=nne2*(z2-1)+nne*(y2-1)+x2;
if ((x2 >= 1) && (y2 >= 1) && (z2 >= 1) && (x2 <= nne) && (y2 <= nne) && (z2 <= nne)) {
  nnz = nnz+1;
}
//////////////////////////////
x2=x1+1;
y2=y1+1;
z2=z1+1;
j=nne2*(z2-1)+nne*(y2-1)+x2;
if ((x2 >= 1) && (y2 >= 1) && (z2 >= 1) && (x2 <= nne) && (y2 <= nne) && (z2 <= nne)) {
  nnz = nnz+1;
}


    }
  }
}
printf("nnz %d\n",nnz); //DEBUG

t2 = clock(); // TIMER

// printf("t2 %ld \n",t2);// DEBUG
printf("CLOCKS_PER_SEC %ld \n",CLOCKS_PER_SEC);// DEBUG
sec = (double)(t2 - t1) / CLOCKS_PER_SEC; // TIME
printf("user time for count number of non-zero element (sec) %f \n", sec ); // TIMER

// allocate irow, jcol matrix

t1 = clock(); // TIMER

irow = malloc(sizeof(int) * nnz);
jcol = malloc(sizeof(int) * nnz);

// set irow jcol matrix

k=0;

for (z1=1; z1<=nne; z1++){
  for (y1=1; y1<=nne; y1++){
    for (x1=1; x1<=nne; x1++){
      i = nne2*(z1 - 1) + nne*(y1 - 1) + x1;

//////////////////////////////////////////////////////////////////// z=-1
x2=x1-1;
y2=y1-1;
z2=z1-1;
j=nne2*(z2-1)+nne*(y2-1)+x2;
if ((x2 >= 1) && (y2 >= 1) && (z2 >= 1) && (x2 <= nne) && (y2 <= nne) && (z2 <= nne)) {
  irow[k]=i;
  jcol[k]=j;
  k++;
//printf("nnz++ x2 %d y2 %d z2 %d i %d\n",x2,y2,z2,i ); //DEBUG
}
//////////////////////////////
x2=x1;
y2=y1-1;
z2=z1-1;
j=nne2*(z2-1)+nne*(y2-1)+x2;
if ((x2 >= 1) && (y2 >= 1) && (z2 >= 1) && (x2 <= nne) && (y2 <= nne) && (z2 <= nne)) {
  irow[k]=i;
  jcol[k]=j;
  k++;
}
//////////////////////////////
x2=x1+1;
y2=y1-1;
z2=z1-1;
j=nne2*(z2-1)+nne*(y2-1)+x2;
if ((x2 >= 1) && (y2 >= 1) && (z2 >= 1) && (x2 <= nne) && (y2 <= nne) && (z2 <= nne)) {
  irow[k]=i;
  jcol[k]=j;
  k++;
}
////////////////////////////////////////////////////////////
x2=x1-1;
y2=y1;
z2=z1-1;
j=nne2*(z2-1)+nne*(y2-1)+x2;
if ((x2 >= 1) && (y2 >= 1) && (z2 >= 1) && (x2 <= nne) && (y2 <= nne) && (z2 <= nne)) {
  irow[k]=i;
  jcol[k]=j;
  k++;
}
//////////////////////////////
x2=x1;
y2=y1;
z2=z1-1;
j=nne2*(z2-1)+nne*(y2-1)+x2;
if ((x2 >= 1) && (y2 >= 1) && (z2 >= 1) && (x2 <= nne) && (y2 <= nne) && (z2 <= nne)) {
  irow[k]=i;
  jcol[k]=j;
  k++;
}
//////////////////////////////
x2=x1+1;
y2=y1;
z2=z1-1;
j=nne2*(z2-1)+nne*(y2-1)+x2;
if ((x2 >= 1) && (y2 >= 1) && (z2 >= 1) && (x2 <= nne) && (y2 <= nne) && (z2 <= nne)) {
  irow[k]=i;
  jcol[k]=j;
  k++;
}
////////////////////////////////////////////////////////////
x2=x1-1;
y2=y1+1;
z2=z1-1;
j=nne2*(z2-1)+nne*(y2-1)+x2;
if ((x2 >= 1) && (y2 >= 1) && (z2 >= 1) && (x2 <= nne) && (y2 <= nne) && (z2 <= nne)) {
  irow[k]=i;
  jcol[k]=j;
  k++;
}
//////////////////////////////
x2=x1;
y2=y1+1;
z2=z1-1;;
j=nne2*(z2-1)+nne*(y2-1)+x2;
if ((x2 >= 1) && (y2 >= 1) && (z2 >= 1) && (x2 <= nne) && (y2 <= nne) && (z2 <= nne)) {
  irow[k]=i;
  jcol[k]=j;
  k++;
}
//////////////////////////////
x2=x1+1;
y2=y1+1;
z2=z1-1;
j=nne2*(z2-1)+nne*(y2-1)+x2;
if ((x2 >= 1) && (y2 >= 1) && (z2 >= 1) && (x2 <= nne) && (y2 <= nne) && (z2 <= nne)) {
  irow[k]=i;
  jcol[k]=j;
  k++;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// z=0
x2=x1-1;
y2=y1-1;
z2=z1;
j=nne2*(z2-1)+nne*(y2-1)+x2;
if ((x2 >= 1) && (y2 >= 1) && (z2 >= 1) && (x2 <= nne) && (y2 <= nne) && (z2 <= nne)) {
  irow[k]=i;
  jcol[k]=j;
  k++;
}
//////////////////////////////
x2=x1;
y2=y1-1;
z2=z1;
j=nne2*(z2-1)+nne*(y2-1)+x2;
if ((x2 >= 1) && (y2 >= 1) && (z2 >= 1) && (x2 <= nne) && (y2 <= nne) && (z2 <= nne)) {
  irow[k]=i;
  jcol[k]=j;
  k++;
}
//////////////////////////////
x2=x1+1;
y2=y1-1;
z2=z1;
j=nne2*(z2-1)+nne*(y2-1)+x2;
if ((x2 >= 1) && (y2 >= 1) && (z2 >= 1) && (x2 <= nne) && (y2 <= nne) && (z2 <= nne)) {
  irow[k]=i;
  jcol[k]=j;
  k++;
}
////////////////////////////////////////////////////////////
x2=x1-1;
y2=y1;
z2=z1;
j=nne2*(z2-1)+nne*(y2-1)+x2;
if ((x2 >= 1) && (y2 >= 1) && (z2 >= 1) && (x2 <= nne) && (y2 <= nne) && (z2 <= nne)) {
  irow[k]=i;
  jcol[k]=j;
  k++;
}
//////////////////////////////
x2=x1;
y2=y1;
z2=z1;
j=nne2*(z2-1)+nne*(y2-1)+x2;
if ((x2 >= 1) && (y2 >= 1) && (z2 >= 1) && (x2 <= nne) && (y2 <= nne) && (z2 <= nne)) {
  irow[k]=i;
  jcol[k]=j;
  k++;
}
//////////////////////////////
x2=x1+1;
y2=y1;
z2=z1;
j=nne2*(z2-1)+nne*(y2-1)+x2;
if ((x2 >= 1) && (y2 >= 1) && (z2 >= 1) && (x2 <= nne) && (y2 <= nne) && (z2 <= nne)) {
  irow[k]=i;
  jcol[k]=j;
  k++;
}
////////////////////////////////////////////////////////////
x2=x1-1;
y2=y1+1;
z2=z1;
j=nne2*(z2-1)+nne*(y2-1)+x2;
if ((x2 >= 1) && (y2 >= 1) && (z2 >= 1) && (x2 <= nne) && (y2 <= nne) && (z2 <= nne)) {
  irow[k]=i;
  jcol[k]=j;
  k++;
}
//////////////////////////////
x2=x1;
y2=y1+1;
z2=z1;
j=nne2*(z2-1)+nne*(y2-1)+x2;
if ((x2 >= 1) && (y2 >= 1) && (z2 >= 1) && (x2 <= nne) && (y2 <= nne) && (z2 <= nne)) {
  irow[k]=i;
  jcol[k]=j;
  k++;
}
//////////////////////////////
x2=x1+1;
y2=y1+1;
z2=z1;
j=nne2*(z2-1)+nne*(y2-1)+x2;
if ((x2 >= 1) && (y2 >= 1) && (z2 >= 1) && (x2 <= nne) && (y2 <= nne) && (z2 <= nne)) {
  irow[k]=i;
  jcol[k]=j;
  k++;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// z=+1;
x2=x1-1;
y2=y1-1;
z2=z1+1;
j=nne2*(z2-1)+nne*(y2-1)+x2;
if ((x2 >= 1) && (y2 >= 1) && (z2 >= 1) && (x2 <= nne) && (y2 <= nne) && (z2 <= nne)) {
  irow[k]=i;
  jcol[k]=j;
  k++;
}
//////////////////////////////
x2=x1;
y2=y1-1;
z2=z1+1;
j=nne2*(z2-1)+nne*(y2-1)+x2;
if ((x2 >= 1) && (y2 >= 1) && (z2 >= 1) && (x2 <= nne) && (y2 <= nne) && (z2 <= nne)) {
  irow[k]=i;
  jcol[k]=j;
  k++;
}
//////////////////////////////
x2=x1+1;
y2=y1-1;
z2=z1+1;
j=nne2*(z2-1)+nne*(y2-1)+x2;
if ((x2 >= 1) && (y2 >= 1) && (z2 >= 1) && (x2 <= nne) && (y2 <= nne) && (z2 <= nne)) {
  irow[k]=i;
  jcol[k]=j;
  k++;
}
////////////////////////////////////////////////////////////
x2=x1-1;
y2=y1;
z2=z1+1;
j=nne2*(z2-1)+nne*(y2-1)+x2;
if ((x2 >= 1) && (y2 >= 1) && (z2 >= 1) && (x2 <= nne) && (y2 <= nne) && (z2 <= nne)) {
  irow[k]=i;
  jcol[k]=j;
  k++;
}
//////////////////////////////
x2=x1;
y2=y1;
z2=z1+1;
j=nne2*(z2-1)+nne*(y2-1)+x2;
if ((x2 >= 1) && (y2 >= 1) && (z2 >= 1) && (x2 <= nne) && (y2 <= nne) && (z2 <= nne)) {
  irow[k]=i;
  jcol[k]=j;
  k++;
}
//////////////////////////////
x2=x1+1;
y2=y1;
z2=z1+1;
j=nne2*(z2-1)+nne*(y2-1)+x2;
if ((x2 >= 1) && (y2 >= 1) && (z2 >= 1) && (x2 <= nne) && (y2 <= nne) && (z2 <= nne)) {
  irow[k]=i;
  jcol[k]=j;
  k++;
}
////////////////////////////////////////////////////////////
x2=x1-1;
y2=y1+1;
z2=z1+1;
j=nne2*(z2-1)+nne*(y2-1)+x2;
if ((x2 >= 1) && (y2 >= 1) && (z2 >= 1) && (x2 <= nne) && (y2 <= nne) && (z2 <= nne)) {
  irow[k]=i;
  jcol[k]=j;
  k++;
}
//////////////////////////////
x2=x1;
y2=y1+1;
z2=z1+1;
j=nne2*(z2-1)+nne*(y2-1)+x2;
if ((x2 >= 1) && (y2 >= 1) && (z2 >= 1) && (x2 <= nne) && (y2 <= nne) && (z2 <= nne)) {
  irow[k]=i;
  jcol[k]=j;
  k++;
}
//////////////////////////////
x2=x1+1;
y2=y1+1;
z2=z1+1;
j=nne2*(z2-1)+nne*(y2-1)+x2;
if ((x2 >= 1) && (y2 >= 1) && (z2 >= 1) && (x2 <= nne) && (y2 <= nne) && (z2 <= nne)) {
  irow[k]=i;
  jcol[k]=j;
  k++;
}

    }
  }
}
t2 = clock(); // TIMER
sec = (double)(t2 - t1) / CLOCKS_PER_SEC; // TIME
printf("user time for create irow jcol matrix (sec) %f \n", sec ); // TIMER

/*
// DEBUG
  for (i=0; i<k; i++) {
    printf("irow[%d] = %d \n", i, irow[i]);
    printf("jcol[%d] = %d \n", i, jcol[i]);
    printf("\n");
  }
*/

printf("now making CRS format\n");

t1 = clock(); //TIMER

csr_idx  = malloc(sizeof(int) * (nnz + 1));
csr_item = malloc(sizeof(int) * nnz);
csr_val  = malloc(sizeof(double) * nnz * coef_block);

csr_idx[0]=0;

row = 1;
for (i=0; i<nnz; i++){

// matrix format error check
  if(row > irow[i]) {
    printf("ERROR: row is not sorted\n");
    exit(1);
  }
  if(row < irow[i]-1) {
    printf("ERROR: vacant row\n");
    exit(1);
  }

// new row?
  if (row == irow[i]-1) {
    csr_idx[row] = i;
  }

  row = irow[i];
  csr_item[i] = jcol[i]; // array element i 1-offset
  csr_val[9*i+0]  = mev;
  csr_val[9*i+1]  = mev;
  csr_val[9*i+2]  = mev;
  csr_val[9*i+3]  = mev;
  csr_val[9*i+4]  = mev;
  csr_val[9*i+5]  = mev;
  csr_val[9*i+6]  = mev;
  csr_val[9*i+7]  = mev;
  csr_val[9*i+8]  = mev;
}
csr_idx[nrows] = nnz;


/*
// DEBUG
printf("test: now write CSR format\n");
for (i=1; i<=nrows; i++){
  jsta = csr_idx[i-1];
  jend = csr_idx[i];
  for (j=jsta; j<jend; j++){
    printf("j: %d\n",j);
    printf("row: %d column: %d \n", i, csr_item[j]);
    printf("value: %f %f %f %f %f %f %f %f %f \n", csr_val[9*j+0], csr_val[9*j+1], csr_val[9*j+2], csr_val[9*j+3], csr_val[9*j+4], csr_val[9*j+5], csr_val[9*j+6], csr_val[9*j+7], csr_val[9*j+8]); // DEBUG
  }
}
*/


t2 = clock(); // TIMER
sec = (double)(t2 - t1) / CLOCKS_PER_SEC; // TIME
printf("user time for create CRS matrix (sec) %f \n", sec ); // TIMER

/////////////////////////////////////////////////////////////////
//
// set vector b for Ab=x
//

t1 = clock(); // TIMER

b_val  = malloc(sizeof(double) * (nrows*vec_coef_block + 1)); // b_val[1] .. b_val[nrows*vec_coef_block] will be used.

for (i=1; i<=nrows*vec_coef_block; i++){
  b_val[i] = vev;
}

/*
//DEBUG
printf("test: now write b vector\n");
for (i=1; i<=nrows; i++){
  printf("%f\n",b_val[i]);
}
*/

t2 = clock(); // TIMER
sec = (double)(t2 - t1) / CLOCKS_PER_SEC; // TIME
printf("user time for create vector B (sec) %f \n", sec ); // TIMER

/////////////////////////////////////////////////////////////////
//
// matrix vector product
//

y = malloc(sizeof(double) * (nrows*vec_coef_block + 1)); // y[1] .. y[nrows*vec_coef_block] will be used.

printf("now doing SpMV calculation %d times\n", itermax);

sec=0; // TIMER

unsigned long int count = 0; //DEBUG

for (iter=0; iter<itermax; iter++){

// printf("iter %d\n",iter);

// clean up y
  for (i=1; i<=nrows*vec_coef_block; i++){
    y[i]=0.0;
  }

// change b vector for each iteration // this code is inserted to prevent compiler optimization for iter loop.
for (i=1; i<=nrows*vec_coef_block; i++){    //DEBUG
  b_val[i] = vev*iter;    //DEBUG
}    //DEBUG

// SpMV
t1 = clock(); // TIME
// #pragma omp parallel for private(i, j, jsta, jend, inod, yv1, yv2, yv3, b1, b2, b3)
  for (i=1; i<=nrows; i++){
    jsta = csr_idx[i-1];
    jend = csr_idx[i];

    yv1 = 0;
    yv2 = 0;
    yv3 = 0;

    for (j=jsta; j<jend; j++){
      inod = csr_item[j];

      b1 = b_val[3*inod-2];
      b2 = b_val[3*inod-1];
      b3 = b_val[3*inod  ];

      yv1 += csr_val[9*j+0]*b1 + csr_val[9*j+1]*b2 + csr_val[9*j+2]*b3;
      yv2 += csr_val[9*j+3]*b1 + csr_val[9*j+4]*b2 + csr_val[9*j+5]*b3;
      yv3 += csr_val[9*j+6]*b1 + csr_val[9*j+7]*b2 + csr_val[9*j+8]*b3;

count++; //DEBUG

    }
    y[3*i-2] = yv1;
    y[3*i-1] = yv2;
    y[3*i  ] = yv3;
  }
t2 = clock(); // TIME
sec += (double)(t2 - t1) / CLOCKS_PER_SEC ; // TIME


}
printf("user time for SpMV (sec) %f \n", sec); // TIME
printf("count %ld\n",count);

printf("test: now write y vector\n"); // DEBUG
for (i=1; i<=nrows*vec_coef_block; i++){ // DEBUG
  printf("%f\n",y[i]); // DEBUG
} // DEBUG


}
