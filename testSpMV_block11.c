/////////////////////////////////////////////////////////////////////////////////
//
// SpMV calculation sample code. 
//
// This code run Sparse Matrix Vector product and measure the calculation time. 
// The matrix data structure mimic the stiffness matrix used in Finite Element Method. 
// Calculation model is a cube, divided by hexa linear elements. 
// Each edge is divided as nnodes_edge nodes. 
// 
// Usage: //////////////////////////////////////////////////////////////////
// Just compile and execute. 
// 
// icc -fast testSpMV_block11.c && ./a.out
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
// 4798.080000                                      // y vector of Ab=y calculation (50th iteration) is shown.
// 7197.120000                                      // It prevent compiler optimization of iter loop and 
// 7197.120000                                      // useful to check calculation result.
// 7197.120000
// 7197.120000
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

// vector
  double *b_val;
  double *y;

  int inod;

// time
  clock_t t1, t2;
  double sec;


/////////////////////////////////////////////////////////////////////////////////////////////////////////
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

// printf("t1 %ld \n",t1);// DEBUG

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

//// DEBUG
//  for (i=0; i<k; i++) {
//    printf("irow[%d] = %d \n", i, irow[i]);
//    printf("jcol[%d] = %d \n", i, jcol[i]);
//    printf("\n");
//  }

printf("now making CRS format\n");

t1 = clock(); //TIMER

csr_idx  = malloc(sizeof(int) * (nnz + 1));
csr_item = malloc(sizeof(int) * nnz);
csr_val  = malloc(sizeof(double) * nnz);

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
  csr_val[i]  = mev;     // matrix element
}
csr_idx[nrows] = nnz;

//// DEBUG
//printf("test: now write CSR format\n");
//for (i=1; i<=nrows; i++){
//  jsta = csr_idx[i-1];
//  jend = csr_idx[i];
//  for (j=jsta; j<jend; j++){
//    printf("row: %d column: %d value: %f \n", i, csr_item[j], csr_val[j]);
//  }
//}

t2 = clock(); // TIMER
sec = (double)(t2 - t1) / CLOCKS_PER_SEC; // TIME
printf("user time for create CRS matrix (sec) %f \n", sec ); // TIMER

/////////////////////////////////////////////////////////////////
//
// set vector b for Ab=y
//

t1 = clock(); // TIMER

b_val  = malloc(sizeof(double) * (nrows+1)); // b_val[1] .. b_val[nrows] will be used.

for (i=1; i<=nrows; i++){
  b_val[i] = vev;
}

////DEBUG
//printf("test: now write X vector\n");
//for (i=1; i<=nrows; i++){
//  printf("%f\n",b_val[i]);
//}

t2 = clock(); // TIMER
sec = (double)(t2 - t1) / CLOCKS_PER_SEC; // TIME
printf("user time for create vector B (sec) %f \n", sec ); // TIMER

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// matrix vector product
//

y = malloc(sizeof(double) * (nrows+1)); // y[1] .. y[nrows] will be used.

printf("now doing SpMV calculation %d times\n", itermax);

sec=0; // TIMER


for (iter=0; iter<itermax; iter++){

// clean up y
  for (i=1; i<=nrows; i++){
    y[i]=0.0;
  }

// change b value // this code is inserted to prevent compiler optimization for iter loop.
for (i=1; i<=nrows; i++){
  b_val[i] = vev * iter;
}

// SpMV
t1 = clock(); // TIME
// #pragma omp parallel for private(i, j, jsta, jend, inod)
  for (i=1; i<=nrows; i++){
    jsta = csr_idx[i-1];
    jend = csr_idx[i];
    for (j=jsta; j<jend; j++){
      inod = csr_item[j];
      y[i] += csr_val[j]*b_val[inod];
    }
  }
t2 = clock(); // TIME
sec += (double)(t2 - t1) / CLOCKS_PER_SEC ; // TIME
}
printf("user time for SpMV (sec) %f \n", sec); // TIME

printf("test: now write y vector\n"); // DEBUG
for (i=1; i<=nrows; i++){ // DEBUG
  printf("%f\n",y[i]); // DEBUG
} // DEBUG

}
