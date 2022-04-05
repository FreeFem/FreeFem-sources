//ff-c++-LIBRARY-dep: cxx11
//ff-c++-cpp-dep:

#include "ff++.hpp"

static void Load_Init() {
    if (mpirank == 0)
        cerr << " ++ WARNING: SLEPc has been superseded by PETSc" << endl;
}

LOADFUNC(Load_Init)
