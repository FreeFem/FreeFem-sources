#ifndef PETSC_HPP_
#define PETSC_HPP_

#if PETSC_VERSION_LT(3,7,0)
#define FFPetscOptionsInsert(a,b,c) PetscOptionsInsert(a,b,c)
#else
#define FFPetscOptionsInsert(a,b,c) PetscOptionsInsert(NULL,a,b,c)
#endif

#if PETSC_VERSION_LT(3,6,0)
#define MatCreateVecs MatGetVecs
#endif

#include "common.hpp"

namespace PETSc {
template<class HpddmType>
class DistributedCSR {
    public:
        HpddmType*                  _A;
        Mat                     _petsc;
        Vec                         _x;
        ISLocalToGlobalMapping   _rmap;
        VecScatter            _scatter;
        Vec                     _isVec;
        KSP                       _ksp;
        unsigned int*             _num;
        unsigned int            _first;
        unsigned int             _last;
        DistributedCSR() : _A(), _petsc(), _x(), _ksp(), _num(), _first(), _last() { };
        ~DistributedCSR() {
            MatDestroy(&_petsc);
            VecDestroy(&_x);
            KSPDestroy(&_ksp);
            if(_A) {
                if(!std::is_same<HpddmType, HpSchwarz<PetscScalar>>::value) {
                    ISLocalToGlobalMappingDestroy(&_rmap);
                    VecDestroy(&_isVec);
                    VecScatterDestroy(&_scatter);
                }
                else
                    _A->clearBuffer();
                delete _A;
                _A = nullptr;
            }
            delete [] _num;
            _num = nullptr;
        }
};
}
#endif
