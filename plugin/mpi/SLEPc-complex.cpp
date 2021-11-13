//ff-c++-LIBRARY-dep: cxx11 slepccomplex petsccomplex [mkl|blas] hpddm mpi
//ff-c++-cpp-dep:

#include "ff++.hpp"

static void Load_Init() {
    if (mpirank == 0)
        cerr << " ++ WARNING: SLEPc-complex has been superseded by PETSc-complex" << endl;
}

LOADFUNC(Load_Init)
