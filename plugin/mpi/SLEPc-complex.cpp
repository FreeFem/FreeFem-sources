/* clang-format off */
//ff-c++-LIBRARY-dep: cxx11
//ff-c++-cpp-dep:
/* clang-format on */

#include "ff++.hpp"

static void Load_Init() {
    if (mpirank == 0)
        cerr << " ++ WARNING: SLEPc-complex has been superseded by PETSc-complex" << endl;
}

LOADFUNC(Load_Init)
