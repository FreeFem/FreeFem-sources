//ff-c++-LIBRARY-dep: cxx11 hpddm petsc slepc [mumps parmetis ptscotch scotch scalapack|umfpack] [mkl|blas] mpi
//ff-c++-cpp-dep:
#ifdef WITH_mkl
#define HPDDM_MKL 1
#define MKL_PARDISOSUB
#elif defined(WITH_mumps)
#define MUMPSSUB
#else
#define SUITESPARSESUB
#endif
#define WITH_slepc

#define  PETScandSLEPc 1
#include "PETSc-code.hpp"
#include "SLEPc-code.hpp"
LOADFUNC(Init)
