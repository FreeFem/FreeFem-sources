//ff-c++-LIBRARY-dep: cxx11 [slepccomplex petsccomplex|petsccomplex] [mkl|blas] [bemtool boost] hpddm mpi
//ff-c++-cpp-dep:

#define  PETScandSLEPc 1
#include "PETSc-code.hpp"
#include "SLEPc-code.hpp"
LOADFUNC(Init)
