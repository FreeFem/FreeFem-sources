//ff-c++-LIBRARY-dep: petsc mpi

#include <ff++.hpp>
#include <petsc.h>

long initialized() {
    PetscBool isInitialized;
    PetscInitialized(&isInitialized);
    return static_cast<long>(isInitialized);
}

static void Init_function( ) {
  Global.Add("PetscInitialized", "(", new OneOperator0< long >(initialized));
}
LOADFUNC(Init_function)
