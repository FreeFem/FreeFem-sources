#ifndef _ALL_IN_ONE_
#define _ALL_IN_ONE_
#endif
//ff-c++-LIBRARY-dep: cxx11   hpddm [petsc|mumps parmetis  ptscotch scotch] umfpack amd  scalapack blas [mkl]   mpifc  fc mpi  pthread
//ff-c++-cpp-dep:
// mumps est avec petsc ..
#define SCHWARZ
//#define BDD

#ifdef WITH_mumps
#define MUMPSSUB
#define DMUMPS
#endif

#ifdef WITH_mkl
#define HPDDM_MKL 1
#endif


#ifndef WITH_petsc
#pragma message("schwarz plugin compile without PETSc")
#else
#define WITH_PETSC
#define MUMPSSUB
#define DMUMPS
#endif

#include <mpi.h>
#include <HPDDM.hpp>
#ifndef __powerpc__
#if defined(__clang__) && !defined(__INTEL_COMPILER)
    #pragma clang diagnostic push
    #pragma clang diagnostic ignored "-Wdeprecated"
    #pragma clang diagnostic ignored "-Wunsequenced"
    #pragma clang diagnostic ignored "-Wdangling-else"
    #pragma clang diagnostic ignored "-Wcomment"
    #pragma clang diagnostic ignored "-Wunused-member-function"
    #pragma clang diagnostic ignored "-Wmismatched-tags"
    #pragma clang diagnostic ignored "-Woverloaded-virtual"
    #pragma clang diagnostic ignored "-Wreorder"
    #pragma clang diagnostic ignored "-Wunused-variable"
    #pragma clang diagnostic ignored "-Wunused-parameter"
#if __has_warning("-Wundefined-bool-conversion")
    #pragma clang diagnostic ignored "-Wundefined-bool-conversion"
#endif
#elif !defined(__INTEL_COMPILER)
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wdeprecated"
    #pragma GCC diagnostic ignored "-Wcomment"
    #pragma GCC diagnostic ignored "-Wunused-local-typedefs"
    #pragma GCC diagnostic ignored "-Woverloaded-virtual"
    #pragma GCC diagnostic ignored "-Wreorder"
    #pragma GCC diagnostic ignored "-Wunused-variable"
    #pragma GCC diagnostic ignored "-Wunused-parameter"
    #pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#endif
#include <ff++.hpp>
#include <AFunction_ext.hpp>
#if defined(__clang__) && !defined(__INTEL_COMPILER)
    #pragma clang diagnostic pop
#elif !defined(__INTEL_COMPILER)
    #pragma GCC diagnostic pop
#endif
#endif
#include <vector>
#include <cmath>

template<class T>
class STL {
    T* const _it;
    const int _size;
    public:
        STL(const KN<T>& v) : _it(v), _size(v.size()) { };
        int size() const {
            return _size;
        }
        T* begin() const {
            return _it;
        }
        T* end() const {
            return _it + _size;
        }
};
template<class K>
class Pair {
    public:
        Pair() : p() { };
        std::pair<MPI_Request, const K*>* p;
        void init() {
        }
        void destroy() {
        }
};

#include "skeleton.cpp"
#ifdef SCHWARZ
#include "schwarz.hpp"
#endif
#if defined(BDD) || defined(FETI)
//#include "substructuring.cpp"
#endif
//#include "scotch.cpp"
#include "removeDOF.cpp"
#include "mymsh3.cpp"
//#include "brainmodel.cpp"
/*
#ifdef __powerpc__
#include "../../ff++/examples++-load/metis.cpp"
#include "../../ff++/examples++-load/msh3.cpp"
#include "../../ff++/examples++-load/Element_P3.cpp"
#include "../../ff++/examples++-load/Element_P4.cpp"
#include "../../ff++/examples++-load/gmsh.cpp"
#endif
*/
#ifdef WITH_PETSC
#include "PETSc.hpp"
#endif

static void Init_hpddm() {
  /* 
#ifdef __powerpc__
    Dcl_Type<listMesh3>();
    typedef Mesh *pmesh;
    typedef Mesh3 *pmesh3;

    if (verbosity && mpirank == 0)
        cout << " load: msh3  " << endl;

    TheOperators->Add("+",new OneBinaryOperator_st< Op3_addmesh<listMesh3,pmesh3,pmesh3>  >      );
    TheOperators->Add("+",new OneBinaryOperator_st< Op3_addmesh<listMesh3,listMesh3,pmesh3>  >      );

    TheOperators->Add("=",new OneBinaryOperator_st< Op3_setmesh<false,pmesh3*,pmesh3*,listMesh3>  >     );
    TheOperators->Add("<-",new OneBinaryOperator_st< Op3_setmesh<true,pmesh3*,pmesh3*,listMesh3>  >     );

    Global.Add("change","(",new SetMesh3D);
    Global.Add("movemesh3","(",new Movemesh3D);
    Global.Add("buildlayers","(",new  BuildLayerMesh);
    Global.Add("trunc","(", new Op_trunc_mesh3);

    Global.Add("AddLayers","(",new OneOperator4_<bool, Mesh3 * , KN<double> *,long, KN<double> * >(AddLayers));
  if (verbosity && mpirank == 0)
    cout << " load: mymsh3  " << endl;

  Global.Add("gluemesh3","(",new OneOperator2_<pmesh3, KN<pmesh3>*, long> (GluMesh3tab)); 
  if (verbosity && mpirank == 0)
    cout << " load: brainmodel  " << endl;
    
  Dcl_Type<Brainmodel *>(InitP<Brainmodel>,Destroy<Brainmodel>);
  zzzfff->Add("Brainmodel",atype<Brainmodel*>());
  TheOperators->Add("<-", new OneOperator3_<Brainmodel *,Brainmodel* ,string*, Complex>(&init_Brainmodel));
  atype< Brainmodel * >()->Add("(","",new OneOperator4_<Complex,Brainmodel *,double,double,double>(Brainmodel_eval));
    if(verbosity && mpirank == 0)
        cout << " lood: init metis  " << endl;
    Global.Add("metisdual","(",new OneOperator3_<KN<double> *,KN<double> *,Mesh *,long , E_F_stackF0F0F0_<KN<double> *,KN<double> *,Mesh *,long> >(&partmetis<Mesh,1>));
    Global.Add("metisdual","(",new OneOperator3_<KN<double> *,KN<double> *,Mesh3 *,long , E_F_stackF0F0F0_<KN<double> *,KN<double> *,Mesh3 *,long> >(&partmetis<Mesh3,1>));
    Global.Add("gmshload3","(",new GMSH_LoadMesh3);
    Global.Add("gmshload","(",new GMSH_LoadMesh);
#endif
  
#ifndef __powerpc__
    Global.Add("scotch", "(", new SCOTCH<Mesh, pmesh, long>);
    Global.Add("scotch", "(", new SCOTCH<Mesh3, pmesh3, long>);
#endif
  */
    Global.Add("mytrunc","(", new Op_mytrunc_mesh3);
    Global.Add("gluemesh3","(",new Op_GluMesh3tab); 
    // Global.Add("scotch", "(", new SCOTCH<Mesh, pmesh, double>);
    //Global.Add("scotch", "(", new SCOTCH<Mesh3, pmesh3, double>);
    Global.Add("removeInteraction", "(", new removeInteraction<double>);
    Global.Add("removeInteraction", "(", new removeInteraction<std::complex<double> >);
    Global.Add("buildSkeleton", "(", new Skeleton);
#include "init.hpp"
}

LOADFUNC(Init_hpddm)
