//ff-c++-LIBRARY-dep: cxx11 [slepc|petsc] [mkl|blas] hpddm mpi
//ff-c++-cpp-dep:
#define  PETScandSLEPc 1
#include "PETSc-code.hpp"
#include "SLEPc-code.hpp"
LOADFUNC(Init)
