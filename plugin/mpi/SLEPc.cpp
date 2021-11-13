//ff-c++-LIBRARY-dep: cxx11 slepc petsc [mkl|blas] hpddm mpi
//ff-c++-cpp-dep:

#include "ff++.hpp"

static void Load_Init() {
    if (mpirank == 0)
        cerr << " ++ WARNING: SLEPc has been superseded by PETSc" << endl;
}

LOADFUNC(Load_Init)
