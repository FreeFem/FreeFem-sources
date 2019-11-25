//ff-c++-LIBRARY-dep: cxx11 slepccomplex [mkl|blas] hpddm mpi
//ff-c++-cpp-dep:

#define WITH_slepccomplex
#define PETScandSLEPc 1
#include "PETSc-code.hpp"
#include "SLEPc-code.hpp"
LOADFUNC(Init)
