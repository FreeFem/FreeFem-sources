//ff-c++-LIBRARY-dep: petsc mpi

#include <ff++.hpp>
#include <petsc.h>

long initialized() {
    PetscBool isInitialized;
    PetscErrorCode ierr = PetscInitialized(&isInitialized);
#ifdef _WIN32
    if(!isInitialized) {
        ierr = PetscInitializeNoArguments();
        ierr = PetscInitialized(&isInitialized);
    }
#endif
    (void)ierr;
    return static_cast<long>(isInitialized);
}

static void Init_function( ) {
  Global.Add("PetscInitialized", "(", new OneOperator0< long >(initialized));
}
LOADFUNC(Init_function)
