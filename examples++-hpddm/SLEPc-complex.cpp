//ff-c++-LIBRARY-dep: cxx11 hpddm petsccomplex slepccomplex [mumps parmetis ptscotch scotch scalapack|umfpack] [mkl|blas] mpi
//ff-c++-cpp-dep:

#define WITH_slepccomplex
#define  PETScandSLEPc 1
#include "PETSc-code.hpp"
#include "SLEPc-code.hpp"
LOADFUNC(Init)
