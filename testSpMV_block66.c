/////////////////////////////////////////////////////////////////////////////////
//
// SpMV calculation sample code. 
//
// This code run Sparse Matrix Vector product and measure the calculation time. 
// The matrix data structure mimic the stiffness matrix used in Finite Element Method. 
// Calculation model is a cube, divided by hexa linear elements. 
// Each edge is divided as nnodes_edge nodes. 
// Each non-zero matrix element have 6x6 blocked value.
// 
// Usage: //////////////////////////////////////////////////////////////////
// Just compile and execute. 
// 
// icc -fast testSpMV_block66.c && ./a.out
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
//  28788.480000                                     // y vector of Ab=y calculation (50th iteration) is shown.
//  28788.480000                                     // It prevent compiler optimization of iter loop and 
//  28788.480000                                     // useful to check calculation result.
//  28788.480000
// 
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

  int coef_block=36; //BCSR parameter : number of element in a block (fix to 36)
  int vec_coef_block=6; //BCSR parameter : number of element in a block (fix to 6)

// vector
  double *b_val;
  double *y;

  double b1, b2, b3, b4, b5, b6, yv1, yv2, yv3, yv4, yv5, yv6;

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

//printf("t2 %ld \n",t2);// DEBUG
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
  csr_val[36*i+ 0]  = mev ;
  csr_val[36*i+ 1]  = mev ;
  csr_val[36*i+ 2]  = mev ;
  csr_val[36*i+ 3]  = mev ;
  csr_val[36*i+ 4]  = mev ;
  csr_val[36*i+ 5]  = mev ;
  csr_val[36*i+ 6]  = mev ;
  csr_val[36*i+ 7]  = mev ;
  csr_val[36*i+ 8]  = mev ;
  csr_val[36*i+ 9]  = mev ;
  csr_val[36*i+10]  = mev ;
  csr_val[36*i+11]  = mev ;
  csr_val[36*i+12]  = mev ;
  csr_val[36*i+13]  = mev ;
  csr_val[36*i+14]  = mev ;
  csr_val[36*i+15]  = mev ;
  csr_val[36*i+16]  = mev ;
  csr_val[36*i+17]  = mev ;
  csr_val[36*i+18]  = mev ;
  csr_val[36*i+19]  = mev ;
  csr_val[36*i+20]  = mev ;
  csr_val[36*i+21]  = mev ;
  csr_val[36*i+22]  = mev ;
  csr_val[36*i+23]  = mev ;
  csr_val[36*i+24]  = mev ;
  csr_val[36*i+25]  = mev ;
  csr_val[36*i+26]  = mev ;
  csr_val[36*i+27]  = mev ;
  csr_val[36*i+28]  = mev ;
  csr_val[36*i+29]  = mev ;
  csr_val[36*i+30]  = mev ;
  csr_val[36*i+31]  = mev ;
  csr_val[36*i+32]  = mev ;
  csr_val[36*i+33]  = mev ;
  csr_val[36*i+34]  = mev ;
  csr_val[36*i+35]  = mev ;
}
csr_idx[nrows] = nnz;

//// DEBUG
//printf("test: now write CSR format\n");
//for (i=1; i<=nrows; i++){
//  jsta = csr_idx[i-1];
//  jend = csr_idx[i];
//  for (j=jsta; j<jend; j++){
//    printf("row: %d column: %d \n", i, csr_item[j]);
////    printf("value: %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f \n", 
/*
  csr_val[36*i+ 0], 
  csr_val[36*i+ 1], 
  csr_val[36*i+ 2], 
  csr_val[36*i+ 3], 
  csr_val[36*i+ 4], 
  csr_val[36*i+ 5], 
  csr_val[36*i+ 6], 
  csr_val[36*i+ 7], 
  csr_val[36*i+ 8], 
  csr_val[36*i+ 9], 
  csr_val[36*i+10], 
  csr_val[36*i+11], 
  csr_val[36*i+12], 
  csr_val[36*i+13], 
  csr_val[36*i+14], 
  csr_val[36*i+15], 
  csr_val[36*i+16], 
  csr_val[36*i+17], 
  csr_val[36*i+18], 
  csr_val[36*i+19], 
  csr_val[36*i+20], 
  csr_val[36*i+21], 
  csr_val[36*i+22], 
  csr_val[36*i+23], 
  csr_val[36*i+24], 
  csr_val[36*i+25], 
  csr_val[36*i+26], 
  csr_val[36*i+27], 
  csr_val[36*i+28], 
  csr_val[36*i+29], 
  csr_val[36*i+30], 
  csr_val[36*i+31], 
  csr_val[36*i+32], 
  csr_val[36*i+33], 
  csr_val[36*i+34], 
  csr_val[36*i+35]
    ); // DEBUG
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

////DEBUG
//printf("test: now write b vector\n");
//for (i=1; i<=nrows; i++){
//  printf("%f\n",b_val[i]);
//}

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


for (iter=0; iter<itermax; iter++){

// clean up y
  for (i=1; i<=nrows*vec_coef_block; i++){
    y[i]=0.0;
  }

// change b value // this code is inserted to prevent compiler optimization for iter loop.
for (i=1; i<=nrows*vec_coef_block; i++){
  b_val[i] = vev * iter;
}

// SpMV
t1 = clock(); // TIME
// #pragma omp parallel for private(i, j, jsta, jend, inod, yv1, yv2, yv3, yv4, yv5, yv6, b1, b2, b3, b4, b5, b6)
  for (i=1; i<=nrows; i++){
    jsta = csr_idx[i-1];
    jend = csr_idx[i];

    yv1 = 0;
    yv2 = 0;
    yv3 = 0;
    yv4 = 0;
    yv5 = 0;
    yv6 = 0;

    for (j=jsta; j<jend; j++){
      inod = csr_item[j];

      b1 = b_val[6*inod-5];
      b2 = b_val[6*inod-4];
      b3 = b_val[6*inod-3];
      b4 = b_val[6*inod-2];
      b5 = b_val[6*inod-1];
      b6 = b_val[6*inod  ];

      yv1 += csr_val[36*j+ 0]*b1 + csr_val[36*j+ 1]*b2 + csr_val[36*j+ 2]*b3 + csr_val[36*j+ 3]*b4 + csr_val[36*j+ 4]*b5 + csr_val[36*j+ 5]*b6;
      yv2 += csr_val[36*j+ 6]*b1 + csr_val[36*j+ 7]*b2 + csr_val[36*j+ 8]*b3 + csr_val[36*j+ 9]*b4 + csr_val[36*j+10]*b5 + csr_val[36*j+11]*b6;
      yv3 += csr_val[36*j+12]*b1 + csr_val[36*j+13]*b2 + csr_val[36*j+14]*b3 + csr_val[36*j+15]*b4 + csr_val[36*j+16]*b5 + csr_val[36*j+17]*b6;
      yv4 += csr_val[36*j+18]*b1 + csr_val[36*j+19]*b2 + csr_val[36*j+20]*b3 + csr_val[36*j+21]*b4 + csr_val[36*j+22]*b5 + csr_val[36*j+23]*b6;
      yv5 += csr_val[36*j+24]*b1 + csr_val[36*j+25]*b2 + csr_val[36*j+26]*b3 + csr_val[36*j+27]*b4 + csr_val[36*j+28]*b5 + csr_val[36*j+29]*b6;
      yv6 += csr_val[36*j+30]*b1 + csr_val[36*j+31]*b2 + csr_val[36*j+32]*b3 + csr_val[36*j+33]*b4 + csr_val[36*j+34]*b5 + csr_val[36*j+35]*b6;

    }
    y[6*i-5] = yv1;
    y[6*i-4] = yv2;
    y[6*i-3] = yv3;
    y[6*i-2] = yv4;
    y[6*i-1] = yv5;
    y[6*i  ] = yv6;

  }
t2 = clock(); // TIME
sec += (double)(t2 - t1) / CLOCKS_PER_SEC ; // TIME
}
printf("user time for SpMV (sec) %f \n", sec); // TIME


printf("test: now write y vector\n"); // DEBUG
for (i=1; i<=nrows*vec_coef_block; i++){ // DEBUG
  printf("%f\n",y[i]); // DEBUG
} // DEBUG


}
