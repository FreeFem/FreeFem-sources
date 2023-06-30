//ff-c++-LIBRARY-dep: [slepc petsc|petsc] [mkl|blas] hpddm mpi
//ff-c++-cpp-dep:
#define  PETScandSLEPc 1
#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wunused-result"
#elif defined(__GNUC__) || defined(__GNUG__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wunused-result"
#endif
#include "PETSc-code.hpp"
#include "SLEPc-code.hpp"
#if defined(__clang__)
  #pragma clang diagnostic pop
#elif defined(__GNUC__) || defined(__GNUG__)
  #pragma GCC diagnostic pop
#endif
LOADFUNC(Init)
