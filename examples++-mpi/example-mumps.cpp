// g++ `../examples++-load/ff-get-dep -ff example-mumps.cpp` -o example-mumps
//ff-c++-LIBRARY-dep:  mumps parmetis ptscotch  scalapack blacs  mpifc  fc   mpi blas 
/*
 *
 *  This file is part of MUMPS 4.9.2, built on Thu Nov  5 07:05:08 UTC 2009
 *
 */
/* Example program using the C interface to the 
 * double real arithmetic version of MUMPS, dmumps_c.
 * We solve the system A x = RHS with
 *   A = diag(1 2) and RHS = [1 4]^T
 * Solution is [1 2]^T */
#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include "dmumps_c.h"
#define JOB_INIT -1
#define JOB_END -2
#define USE_COMM_WORLD -987654

#if defined(MAIN_COMP)
/*
 * Some Fortran compilers (COMPAQ fort) define main inside
 * their runtime library while a Fortran program translates
 * to MAIN_ or MAIN__ which is then called from "main". This
 * is annoying because MAIN__ has no arguments and we must
 * define argc/argv arbitrarily !!
 */
int MAIN__();
int MAIN_()
  {
    return MAIN__();
  }

int MAIN__()
{
  int argc=1;
  char * name = "c_example";
  char ** argv ;
#else
int main(int argc, char ** argv)
{
#endif
  DMUMPS_STRUC_C id;
  int n = 2;
  int nz = 2;
  int irn[] = {1,2};
  int jcn[] = {1,2};
  double a[2];
  double rhs[2];

  int myid, ierr;
#if defined(MAIN_COMP)
  argv = &name;
#endif
  ierr = MPI_Init(&argc, &argv);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  /* Define A and rhs */
  rhs[0]=1.0;rhs[1]=4.0;
  a[0]=1.0;a[1]=2.0;

  /* Initialize a MUMPS instance. Use MPI_COMM_WORLD */
  id.job=JOB_INIT; id.par=1; id.sym=0;id.comm_fortran=USE_COMM_WORLD;
  dmumps_c(&id);
  /* Define the problem on the host */
  if (myid == 0) {
    id.n = n; id.nz =nz; id.irn=irn; id.jcn=jcn;
    id.a = a; id.rhs = rhs;
  }
#define ICNTL(I) icntl[(I)-1] /* macro s.t. indices match documentation */
/* No outputs */
  id.ICNTL(1)=-1; id.ICNTL(2)=-1; id.ICNTL(3)=-1; id.ICNTL(4)=0;
/* Call the MUMPS package. */
  id.job=6;
  dmumps_c(&id);
  id.job=JOB_END; dmumps_c(&id); /* Terminate instance */
  if (myid == 0) {
    printf("Solution is : (%8.2f  %8.2f)\n", rhs[0],rhs[1]);
  }
  ierr = MPI_Finalize();
  return 0;
}
