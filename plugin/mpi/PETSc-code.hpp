#include "petsc.h"

#include "PETSc.hpp"
#include "compositeFESpace.hpp"

typedef PETSc::DistributedCSR< HpSchwarz< PetscScalar > > Dmat;
typedef PETSc::DistributedCSR< HpSchwarz< PetscReal > > DmatR;
typedef PETSc::DistributedCSR< HpSchwarz< PetscComplex > > DmatC;
typedef PETSc::DistributedCSR< HpSchur< PetscScalar > > Dbddc;
typedef PETSc::DistributedCSR< HpSchur< PetscReal > > DbddcR;
typedef PETSc::DistributedCSR< HpSchur< PetscComplex > > DbddcC;

#if defined(WITH_bemtool) && defined(WITH_htool) && defined(PETSC_HAVE_HTOOL)
namespace PETSc {

template<typename P, typename MeshBemtool>
struct HtoolCtx {
    htool::VirtualGenerator<PetscScalar>* generator;
    bemtool::Dof<P>* dof;
    MeshBemtool* mesh;
    bemtool::Geometry* node;
    HtoolCtx() : generator(), dof(), mesh(), node() { }
    ~HtoolCtx() {
        delete generator;
        delete dof;
        delete mesh;
        delete node;
    }
};

template<typename P, typename MeshBemtool>
static PetscErrorCode GenEntriesFromCtx(PetscInt sdim,PetscInt M,PetscInt N, const PetscInt *const J, const PetscInt *const K, PetscScalar* ptr,void *ctx) {
    PetscFunctionBeginUser;
    HtoolCtx<P,MeshBemtool>* user = reinterpret_cast<HtoolCtx<P,MeshBemtool>*>(ctx);
    user->generator->copy_submatrix(M,N,J,K,ptr);
    PetscFunctionReturn(PETSC_SUCCESS);
}

template<typename P, typename MeshBemtool>
static PetscErrorCode DestroyHtoolCtx(void *ctx) {
    HtoolCtx<P,MeshBemtool>* user = (HtoolCtx<P,MeshBemtool>*)ctx;

    PetscFunctionBeginUser;
    delete user;
    PetscFunctionReturn(PETSC_SUCCESS);
}

template<template<typename P, typename MeshBemtool> class Gen, typename P, typename MeshBemtool, typename std::enable_if< std::is_same< HtoolCtx<P, MeshBemtool>, Gen<P, MeshBemtool> >::value >::type* = nullptr>
htool::VirtualGenerator<PetscScalar>* get_gen(Gen<P, MeshBemtool>* generator) {
    return generator->generator;
}

template<template<typename P, typename MeshBemtool> class Gen, typename P, typename MeshBemtool, typename std::enable_if< !std::is_same< HtoolCtx<P, MeshBemtool>, Gen<P, MeshBemtool> >::value >::type* = nullptr>
htool::VirtualGenerator<PetscScalar>* get_gen(Gen<P, MeshBemtool>* generator) {
    return generator;
}

template<class Matrix, template<typename P, typename MeshBemtool> class Gen, typename P, typename MeshBemtool, class R = PetscScalar, typename std::enable_if< std::is_same< Matrix, Dmat >::value >::type* = nullptr>
void Assembly(Matrix* A, Gen<P, MeshBemtool>* generator, string compressor,vector<double> &p1,vector<double> &p2,MPI_Comm comm,int dim,bool sym = false) {
    PetscInt m, M;
    KSPDestroy(&A->_ksp);
    if(A->_vS) {
        for(int i = 0; i < A->_vS->size(); ++i)
            MatDestroy(&(*A->_vS)[i]);
        delete A->_vS;
        A->_vS = nullptr;
    }
    if(!A->_petsc) {
        MatCreateHtoolFromKernel(PETSC_COMM_SELF,p2.size()/3,p1.size()/3,p2.size()/3,p1.size()/3,3,p2.data(),p1.data(),nullptr,get_gen(generator),&A->_petsc);
    }
    else {
        Mat B;
        PetscInt m, n, M, N, rbegin, cbegin;
        MatGetLocalSize(A->_petsc, &m, &n);
        MatGetSize(A->_petsc, &M, &N);
        MatGetOwnershipRange(A->_petsc, &rbegin, NULL);
        MatGetOwnershipRangeColumn(A->_petsc, &cbegin, NULL);
        ffassert(N == p1.size()/3 && M == p2.size()/3);
        MatCreateHtoolFromKernel(PetscObjectComm((PetscObject)A->_petsc),m,n,M,N,3,p2.data()+rbegin*3,p1.data()+cbegin*3,nullptr,get_gen(generator),&B);
        MatHeaderReplace(A->_petsc, &B);
    }
    if(compressor.size()) {
        PetscOptionsInsertString(NULL, compressor.c_str());
    }
    if(sym) {
        MatSetOption(A->_petsc, MAT_SYMMETRIC, PETSC_TRUE);
    }
    MatSetFromOptions(A->_petsc);
    MatAssemblyBegin(A->_petsc, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A->_petsc, MAT_FINAL_ASSEMBLY);
    if(std::is_same<HtoolCtx<P, MeshBemtool>, Gen<P, MeshBemtool>>::value) {
        MatHtoolSetKernel(A->_petsc, GenEntriesFromCtx<P, MeshBemtool>, generator);
    }
}

template<class fes1, class fes2, typename std::enable_if< (fes1::FESpace::Mesh::RdHat::d >= 3) || std::is_same<typename fes1::FESpace::Mesh, Mesh>::value >::type* = nullptr >
void varfBem(const typename fes1::FESpace*& PUh, const typename fes2::FESpace*& PVh, bool same, int VFBEM, Stack stack, const list<C_F0>& bargs, const Data_Sparse_Solver& ds, Dmat* B) {
    ffassert(VFBEM == -1);
}

template<typename P, typename MeshBemtool, typename std::enable_if< std::is_same<P, bemtool::RT0_2D>::value || std::is_same<P, bemtool::P0_1D>::value || std::is_same<P, bemtool::P0_2D>::value >::type* = nullptr >
bemtool::Dof<P>* new_dof(MeshBemtool *mesh) {
    return new bemtool::Dof<P>(*mesh);
}

template<typename P, typename MeshBemtool, typename std::enable_if< !std::is_same<P, bemtool::RT0_2D>::value && !std::is_same<P, bemtool::P0_1D>::value && !std::is_same<P, bemtool::P0_2D>::value >::type* = nullptr >
bemtool::Dof<P>* new_dof(MeshBemtool *mesh) {
    return new bemtool::Dof<P>(*mesh, true);
}

template<typename P, typename MeshBemtool, typename Mesh1, bool T, typename Mesh2, typename FESpace2, typename std::enable_if< !T >::type* = nullptr >
void dispatch(MeshBemtool *mesh, bemtool::Geometry *node, int VFBEM, Stack stack, const list<C_F0>& bargs, const Data_Sparse_Solver& ds, Dmat* B, const Mesh2& ThV, std::vector<double>& p1, std::vector<double>& p2, bool same, const FESpace2& Vh) {
    HtoolCtx<P,MeshBemtool>* ctx = new HtoolCtx<P,MeshBemtool>;
    ctx->dof = new_dof<P>(mesh);
    ctx->mesh = mesh;
    ctx->node = node;
    bemtool::Geometry node_output;
    if (VFBEM == 1) {
        pair<BemKernel*, std::complex<double>> kernel = getBemKernel(stack, bargs);
        BemKernel *Ker = kernel.first;
        std::complex<double> alpha = kernel.second;
        ff_BIO_Generator_Maxwell<PetscScalar>(ctx->generator,Ker,*ctx->dof,alpha);
    }
    else if (VFBEM == 2) {
        BemPotential *Pot = getBemPotential(stack, bargs);
        if (Vh.MaxNbNodePerElement == Mesh2::RdHat::d+1)
            Mesh2Bemtool(ThV,node_output);
        else if (Vh.MaxNbNodePerElement == 1) {
            int m = Vh.NbOfDF;
            typename Mesh2::RdHat pbt(1./(Mesh2::RdHat::d+1),1./(Mesh2::RdHat::d+1));
            for (int i=0; i<m; i++) {
                Fem2D::R3 p = ThV[i](pbt);
                bemtool::R3 q;
                q[0]=p.x; q[1]=p.y; q[2]=p.z;
                node_output.setnodes(q);
            }
        }
        else
            ffassert(0);
        ff_POT_Generator_Maxwell<PetscScalar,P>(ctx->generator,Pot,*ctx->dof,*mesh,node_output);
    }
    if (same) {
        PetscContainer ptr;
        PetscObjectQuery((PetscObject)B->_petsc, "HtoolCtx", (PetscObject*)&ptr);
        PetscContainerDestroy(&ptr);
        Assembly(B,ctx,ds.sparams,p1,p1,MPI_COMM_NULL,Mesh1::RdHat::d+1,ds.sym);
        PetscContainerCreate(PetscObjectComm((PetscObject)B->_petsc), &ptr);
        PetscContainerSetPointer(ptr, ctx);
        PetscContainerSetUserDestroy(ptr, DestroyHtoolCtx<P,MeshBemtool>);
        PetscObjectCompose((PetscObject)B->_petsc, "HtoolCtx", (PetscObject)ptr);
    }
    else {
        Assembly(B,ctx,ds.sparams,p1,p2,MPI_COMM_NULL,Mesh1::RdHat::d+1,ds.sym);
        delete ctx;
    }
}

template<typename P, typename MeshBemtool, typename Mesh1, bool T, typename Mesh2, typename FESpace2, typename std::enable_if< T >::type* = nullptr >
void dispatch(MeshBemtool *mesh, bemtool::Geometry *node, int VFBEM, Stack stack, const list<C_F0>& bargs, const Data_Sparse_Solver& ds, Dmat* B, const Mesh2& ThV, std::vector<double>& p1, std::vector<double>& p2, bool same, const FESpace2& Vh) {
    HtoolCtx<P,MeshBemtool>* ctx = new HtoolCtx<P,MeshBemtool>;
    ctx->dof = new_dof<P>(mesh);
    ctx->mesh = mesh;
    ctx->node = node;
    bemtool::Geometry node_output;
    if (VFBEM == 1) {
        pair<BemKernel*, std::complex<double>> kernel = getBemKernel(stack, bargs);
        BemKernel *Ker = kernel.first;
        std::complex<double> alpha = kernel.second;
        ff_BIO_Generator<PetscScalar,P,Mesh1>(ctx->generator,Ker,*ctx->dof,alpha);
    }
    else if (VFBEM == 2) {
        BemPotential *Pot = getBemPotential(stack, bargs);
        if (Vh.MaxNbNodePerElement == Mesh2::RdHat::d+1)
            Mesh2Bemtool(ThV,node_output);
        else if (Vh.MaxNbNodePerElement == 1) {
            int m = Vh.NbOfDF;
            typename Mesh2::RdHat pbt(1./(Mesh2::RdHat::d+1),1./(Mesh2::RdHat::d+1));
            for (int i=0; i<m; i++) {
                Fem2D::R3 p = ThV[i](pbt);
                bemtool::R3 q;
                q[0]=p.x; q[1]=p.y; q[2]=p.z;
                node_output.setnodes(q);
            }
        }
        else
            ffassert(0);
        ff_POT_Generator<PetscScalar,P,MeshBemtool,Mesh1>(ctx->generator,Pot,*ctx->dof,*mesh,node_output);
    }
    if (same) {
        PetscContainer ptr;
        PetscObjectQuery((PetscObject)B->_petsc, "HtoolCtx", (PetscObject*)&ptr);
        PetscContainerDestroy(&ptr);
        Assembly(B,ctx,ds.sparams,p1,p1,MPI_COMM_NULL,Mesh1::RdHat::d+1,ds.sym);
        PetscContainerCreate(PetscObjectComm((PetscObject)B->_petsc), &ptr);
        PetscContainerSetPointer(ptr, ctx);
        PetscContainerSetUserDestroy(ptr, DestroyHtoolCtx<P,MeshBemtool>);
        PetscObjectCompose((PetscObject)B->_petsc, "HtoolCtx", (PetscObject)ptr);
    }
    else {
        Assembly(B,ctx,ds.sparams,p1,p2,MPI_COMM_NULL,Mesh1::RdHat::d+1,ds.sym);
        delete ctx;
    }
}

template<class fes1, class fes2, typename std::enable_if< (fes1::FESpace::Mesh::RdHat::d < 3) && !std::is_same<typename fes1::FESpace::Mesh, Mesh>::value >::type* = nullptr >
void varfBem(const typename fes1::FESpace*& PUh, const typename fes2::FESpace*& PVh, bool same, int VFBEM, Stack stack, const list<C_F0>& bargs, const Data_Sparse_Solver& ds, Dmat* B) {
    if (VFBEM == 1)
        ffassert(same);
    typedef typename fes1::pfes pfes1;
    typedef typename fes1::FESpace FESpace1;
    typedef typename FESpace1::Mesh Mesh1;

    typedef typename fes2::pfes pfes2;
    typedef typename fes2::FESpace FESpace2;
    typedef typename FESpace2::Mesh Mesh2;

    typedef typename std::conditional<Mesh1::RdHat::d==1, bemtool::Mesh1D, bemtool::Mesh2D>::type MeshBemtool;
    typedef typename std::conditional<Mesh1::RdHat::d==1, bemtool::P0_1D, bemtool::P0_2D>::type P0;
    typedef typename std::conditional<Mesh1::RdHat::d==1, bemtool::P1_1D, bemtool::P1_2D>::type P1;
    typedef typename std::conditional<Mesh1::RdHat::d==1, bemtool::P2_1D, bemtool::P2_2D>::type P2;

    const FESpace1& Uh = *PUh;
    const FESpace2& Vh = *PVh;
    const Mesh1& ThU = Uh.Th; // line
    const Mesh2& ThV = Vh.Th; // colunm
    int n = Uh.NbOfDF;
    int m = Vh.NbOfDF;
    bemtool::Geometry *node = new bemtool::Geometry;
    MeshBemtool *mesh = new MeshBemtool;
    Mesh2Bemtool(ThU, *node, *mesh);
    bemtool::Dof<P1> *dof = new bemtool::Dof<P1>(*mesh,true);
    vector<double> p1;
    p1.reserve(3*n);
    vector<double> p2;
    typename Mesh1::RdHat pbs(1./(Mesh1::RdHat::d+1),1./(Mesh1::RdHat::d+1));
    int Snbv = Uh.TFE[0]->ndfonVertex;
    int Snbe = Uh.TFE[0]->ndfonEdge;
    int Snbt = Uh.TFE[0]->ndfonFace;
    bool SP0 = Mesh1::RdHat::d == 1 ? (Snbv == 0) && (Snbe == 1) && (Snbt == 0) : (Snbv == 0) && (Snbe == 0) && (Snbt == 1);
    bool SP1 = (Snbv == 1) && (Snbe == 0) && (Snbt == 0);
    bool SP2 = (Snbv == 1) && (Snbe == 1) && (Snbt == 0);
    bool SRT0 = (Mesh1::RdHat::d == 2) && (Snbv == 0) && (Snbe == 1) && (Snbt == 0);
    if (SP2) {
        bemtool::Dof<P2> dof(*mesh,true);
        for (int i=0; i<n; i++) {
            const std::vector<bemtool::N2>& j = dof.ToElt(i);
            bemtool::R3 p = dof(j[0][0])[j[0][1]];
            p1.emplace_back(p[0]);
            p1.emplace_back(p[1]);
            p1.emplace_back(p[2]);
        }
    }
    else if (SRT0) {
        bemtool::Dof<bemtool::RT0_2D> dof(*mesh);
        for (int i=0; i<n; i++) {
            const std::vector<bemtool::N2>& j = dof.ToElt(i);
            bemtool::R3 p = dof(j[0][0])[j[0][1]];
            p1.emplace_back(p[0]);
            p1.emplace_back(p[1]);
            p1.emplace_back(p[2]);
        }
    }
    else {
        for (int i=0; i<n; i++) {
            Fem2D::R3 p;
            if (SP1)
                p = ThU.vertices[i];
            else if (SP0)
                p = ThU[i](pbs);
            else {
                if (mpirank == 0) std::cerr << "ff-BemTool error: only P0, P1 and P2 discretizations are available for now." << std::endl;
                ffassert(0);
            }
            p1.emplace_back(p.x);
            p1.emplace_back(p.y);
            p1.emplace_back(p.z);
        }
    }
    if (!same) {
        if(Vh.TFE[0]->N == 1) {
            typename Mesh2::RdHat pbt(1./(Mesh2::RdHat::d+1),1./(Mesh2::RdHat::d+1));
            p2.reserve(3*m);
            for (int i=0; i<m; i++) {
                Fem2D::R3 p;
                if (Vh.MaxNbNodePerElement == Mesh2::RdHat::d+1)
                    p = ThV.vertices[i];
                else if (Vh.MaxNbNodePerElement == 1)
                    p = ThV[i](pbt);
                else {
                    if (mpirank == 0) std::cerr << "ff-BemTool error: only P0 and P1 FEspaces are available for reconstructions." << std::endl;
                    ffassert(0);
                }
                p2.emplace_back(p.x);
                p2.emplace_back(p.y);
                p2.emplace_back(p.z);
            }
        }
        else {
            ffassert(SRT0 && Mesh1::RdHat::d == 2 && VFBEM == 2);
            int nnn = Vh.TFE[0]->N;

            typename Mesh2::RdHat pbt(1./(Mesh2::RdHat::d+1),1./(Mesh2::RdHat::d+1));
            p2.reserve(3*m);

            int mDofScalar = m/nnn; // computation of the dof of one component 

            for (int i=0; i<mDofScalar; i++) {
                Fem2D::R3 p;
                if (Vh.MaxNbNodePerElement == Mesh2::RdHat::d + 1)
                    p = ThV.vertices[i];
                else if (Vh.MaxNbNodePerElement == 1)
                    p = ThV[i](pbt);
                else {
                    if (mpirank == 0) std::cerr << "ff-BemTool error: only P0 and P1 FEspaces are available for reconstructions." << std::endl;
                    ffassert(0);
                }
                for(int iii=0; iii<nnn; iii++){
                    ffassert( nnn*3*i+3*iii+2 < nnn*3*m );
                    p2.emplace_back(p.x);
                    p2.emplace_back(p.y);
                    p2.emplace_back(p.z);
                }
            }
        }
    }
    if (SP1)
        dispatch<P1,MeshBemtool,Mesh1,true>(mesh, node, VFBEM, stack, bargs, ds, B, ThV, p1, p2, same, Vh);
    else if (SP0)
        dispatch<P0,MeshBemtool,Mesh1,true>(mesh, node, VFBEM, stack, bargs, ds, B, ThV, p1, p2, same, Vh);
    else if (SP2)
        dispatch<P2,MeshBemtool,Mesh1,true>(mesh, node, VFBEM, stack, bargs, ds, B, ThV, p1, p2, same, Vh);
    else if (SRT0 && Mesh1::RdHat::d == 2) {
        dispatch<bemtool::RT0_2D,MeshBemtool,Mesh1,false>(mesh, node, VFBEM, stack, bargs, ds, B, ThV, p1, p2, same, Vh);
    }
    else
        ffassert(0);
}
} // namespace PETSc
#endif

namespace PETSc {
  template<class K, class MMesh, class fes1, class fes2 >
  struct varfToMat : public OneOperator {
    class Op : public E_F0mps {
    public:
      Call_FormBilinear<fes1, fes2>* b;
      Expression a;
      AnyType operator()(Stack s) const;

      Op(Expression x, Expression y) : b(new Call_FormBilinear<fes1, fes2>(*dynamic_cast<const Call_FormBilinear<fes1, fes2>*>(y))), a(x) {
          assert(b && b->nargs);
          ffassert(FieldOfForm(b->largs, IsComplexType<upscaled_type<K>>::value) == IsComplexType<upscaled_type<K>>::value);
          
          #if defined(WITH_bemtool) && defined(WITH_htool) && defined(PETSC_HAVE_HTOOL)
          // Check the nbitem of inconnu and test in BemFormBilinear
          checkNbItemFEspacesInconnuAndTest(b->largs,b->N,b->M);
          #endif
      }
      operator aType () const { return atype<Dmat*>(); }
    };
    E_F0* code(const basicAC_F0& args) const {
        return new Op(to<Dmat*>(args[0]), args[1]);
    }
    varfToMat() : OneOperator(atype<Dmat*>(), atype<Dmat*>(), atype<const Call_FormBilinear<fes1, fes2>*>()) {}
  };
  template<class fes1, class fes2, typename std::enable_if< std::is_same< fes1, fes2 >::value >::type* = nullptr >
  void assert_ptr(fes1* pUh, fes2* pVh) {
    ffassert(pUh == pVh);
  }
  template<class fes1, class fes2, typename std::enable_if< !std::is_same< fes1, fes2 >::value >::type* = nullptr >
  void assert_ptr(fes1* pUh, fes2* pVh) { }
  template<class K, class MMesh, class fes1, class fes2>
  AnyType varfToMat<K, MMesh, fes1, fes2>::Op::operator()(Stack stack) const {
    typedef typename fes1::pfes pfes1;
    typedef typename fes1::FESpace FESpace1;
    typedef typename FESpace1::Mesh Mesh1;

    typedef typename fes2::pfes pfes2;
    typedef typename fes2::FESpace FESpace2;
    typedef typename FESpace2::Mesh Mesh2;

    assert(b && b->nargs);
    pfes1* pUh = GetAny<pfes1*>((*b->euh)(stack));
    pfes2* pVh = GetAny<pfes2*>((*b->evh)(stack));
    const FESpace1* PUh = (FESpace1*)**pUh;
    const FESpace2* PVh = (FESpace2*)**pVh;

    Data_Sparse_Solver ds;
    ds.factorize = 0;
    ds.initmat = true;
    int np = OpCall_FormBilinear_np::n_name_param - NB_NAME_PARM_HMAT;
    SetEnd_Data_Sparse_Solver<upscaled_type<K>>(stack, ds, b->nargs, np);

    WhereStackOfPtr2Free(stack) = new StackOfPtr2Free(stack);

    Dmat& B(*GetAny<Dmat*>((*a)(stack)));
    if(!PUh || !PVh)
      return SetAny<Dmat*>(&B);
    const FESpace1& Uh = *PUh;
    const FESpace2& Vh = *PVh;
    const Mesh1& Th = Uh.Th;
    bool same = isSameMesh(b->largs, &Uh.Th, &Vh.Th, stack);
#if defined(WITH_bemtool) && defined(WITH_htool) && defined(PETSC_HAVE_HTOOL)
    int VFBEM = typeVFBEM(b->largs, stack);
#else
    int VFBEM = -1;
#endif
    if (VFBEM == -1) {
      ffassert((std::is_same< fes1, fes2 >::value));
      Matrice_Creuse<upscaled_type<K>> A;
      A.init();
      if(same) {
        if(A.Uh != Uh || A.Vh != Vh) {
          A.Uh = Uh;
          A.Vh = Vh;
          if(ds.sym) {
            A.A.master(new MatriceMorse<upscaled_type<K>>(Vh.NbOfDF, Vh.NbOfDF, 0, ds.sym));
            assert_ptr(&Uh, &Vh);
          }
          else
            A.A.master(new MatriceMorse<upscaled_type<K>>(Vh.NbOfDF, Uh.NbOfDF, 2 * Vh.NbOfDF, 0));
        }
        if(AssembleVarForm<upscaled_type<K>, MatriceCreuse<upscaled_type<K>>, MMesh, FESpace1,FESpace2>(stack, Th, Uh, Vh, ds.sym, A.A, 0, b->largs))
          AssembleBC<upscaled_type<K>, MMesh,FESpace1, FESpace2>(stack, Th, Uh, Vh, ds.sym, A.A, 0, 0, b->largs, ds.tgv);
      }
      else {
        MatriceMorse<upscaled_type<K>> *pMA = new MatriceMorse<upscaled_type<K>>(Vh.NbOfDF, Uh.NbOfDF, 0, ds.sym);
        MatriceMap<upscaled_type<K>>& D = *pMA;
        bool bc = AssembleVarForm<upscaled_type<K>, MatriceMap<upscaled_type<K>>, MMesh, FESpace1,FESpace2>(stack, Th, Uh, Vh, ds.sym, &D, 0, b->largs);
        A.A.master(pMA);
        if(bc)
          AssembleBC<upscaled_type<K>>(stack, Th, Uh, Vh, ds.sym, A.A, 0, 0, b->largs, ds.tgv);
      }
      changeOperatorSimple(&B, &A);
      if(B._A)
          B._A->setMatrix(nullptr);
    }
#if defined(WITH_bemtool) && defined(WITH_htool) && defined(PETSC_HAVE_HTOOL)
    else {
        varfBem<fes1, fes2>(PUh, PVh, same, VFBEM, stack, b->largs, ds, &B);
    }
#endif
    return SetAny<Dmat*>(&B);
  }
} // namespace PETSc

namespace PETSc {
  template< class Type >
  struct _n_User;
  template< class Type >
  using User = _n_User< Type >*;
  template< bool C, class HpddmType,
            typename std::enable_if< std::is_same< HpddmType, Dmat >::value >::type* = nullptr >
  void initPETScStructure(
    HpddmType* ptA, PetscInt bs, PetscBool symmetric,
    KN< typename std::conditional< std::is_same< HpddmType, Dmat >::value, double, long >::type >* ptD) {
    double timing = MPI_Wtime( );
    long long global;
    if (ptD) {
      if(!ptA->_D) {
        PetscReal* d = reinterpret_cast<PetscReal*>(ptD->operator double*());
        if(!std::is_same<upscaled_type<PetscReal>, PetscReal>::value) {
            for(int i = 0; i < ptD->n; ++i)
                d[i] = ptD->operator[](i);
        }
        ptA->_A->restriction(d);
        if (!C) ptA->_A->initialize(d);
        else {
          ptA->_D = new KN<PetscReal>(ptD->n);
          for(int i = 0; i < ptD->n; ++i)
              ptA->_D->operator[](i) = d[i];
          ptA->_A->initialize(*ptA->_D);
        }
      }
      else {
        ptA->_A->initialize(*ptA->_D);
      }
    }
    if (!C && !ptD) global = PETSC_DECIDE;
    else ptA->_A->distributedNumbering(ptA->_num, ptA->_first, ptA->_last, global);
    if (verbosity > 0 && mpirank == 0)
      cout << " --- global numbering created (in " << MPI_Wtime( ) - timing << ")" << endl;
    timing = MPI_Wtime( );
    PetscInt* ia = nullptr;
    PetscInt* ja = nullptr;
    PetscScalar* c = nullptr;
    bool free = ptA->_A->getMatrix( )->HPDDM_ia
                  ? ptA->_A->distributedCSR(ptA->_num, ptA->_first, ptA->_last, ia, ja, c)
                  : false;
    MatCreate(ptA->_A->getCommunicator(), &ptA->_petsc);
    if (bs > 1) MatSetBlockSize(ptA->_petsc, bs);
    MatSetSizes(ptA->_petsc, ptA->_last - ptA->_first, ptA->_last - ptA->_first, global, global);
    bool sym = ptA->_A->getMatrix( )->HPDDM_sym;
    if (ia) {
      if (sym) {
        MatSetType(ptA->_petsc, MATSBAIJ);
        MatSeqSBAIJSetPreallocationCSR(ptA->_petsc, 1, ia, ja, c);
        MatMPISBAIJSetPreallocationCSR(ptA->_petsc, 1, ia, ja, c);
      } else {
        MatSetType(ptA->_petsc, MATAIJ);
        MatSeqAIJSetPreallocationCSR(ptA->_petsc, ia, ja, c);
        MatMPIAIJSetPreallocationCSR(ptA->_petsc, ia, ja, c);
        MatSetOption(ptA->_petsc, MAT_SYMMETRIC, symmetric);
      }
    } else {
      MatSetType(ptA->_petsc, MATAIJ);
      MatSetUp(ptA->_petsc);
    }
    if (free) {
      delete[] ia;
      delete[] ja;
      delete[] c;
    }
    ptA->_A->setBuffer( );
    if (verbosity > 0 && mpirank == 0)
      cout << " --- global CSR created (in " << MPI_Wtime( ) - timing << ")" << endl;
  }
  template< bool C, class HpddmType,
            typename std::enable_if< !std::is_same< HpddmType, Dmat >::value >::type* = nullptr >
  void initPETScStructure(
    HpddmType* ptA, PetscInt& bs, PetscBool symmetric,
    KN< typename std::conditional< std::is_same< HpddmType, Dmat >::value, double, long >::type >* ptD) {
    const HPDDM::MatrixCSR< PetscScalar >* M = ptA->_A->getMatrix( );
    if (!M->HPDDM_sym) cout << "Please assemble a symmetric CSR" << endl;
    double timing = MPI_Wtime( );
    ptA->_A->template renumber< false >(STL< long >(*ptD), nullptr);
    long long global;
    ptA->_A->distributedNumbering(ptA->_num, ptA->_first, ptA->_last, global);
    if (verbosity > 0 && mpirank == 0)
      cout << " --- global numbering created (in " << MPI_Wtime( ) - timing << ")" << endl;
    timing = MPI_Wtime( );
    PetscInt* indices;
    PetscMalloc(sizeof(PetscInt) * M->HPDDM_n / bs, &indices);
    for (unsigned int i = 0; i < M->HPDDM_n; i += bs) indices[i / bs] = ptA->_num[i] / bs;
    ISLocalToGlobalMapping rmap;
    ISLocalToGlobalMappingCreate(ptA->_A->getCommunicator(), bs, M->HPDDM_n / bs, indices, PETSC_OWN_POINTER,
                                 &rmap);
    MatCreateIS(ptA->_A->getCommunicator(), bs, PETSC_DECIDE, PETSC_DECIDE, global, global, rmap, NULL,
                &ptA->_petsc);
    Mat local;
    MatISGetLocalMat(ptA->_petsc, &local);
    MatSetType(local, MATSEQSBAIJ);
    std::vector< std::vector< std::pair< int, PetscScalar > > > transpose(M->HPDDM_n);
    for (int i = 0; i < transpose.size( ); ++i)
      for (int j = M->HPDDM_ia[i]; j < M->HPDDM_ia[i + 1]; ++j) {
        transpose[M->HPDDM_ja[j]].emplace_back(i, M->HPDDM_a[j]);
        if (bs > 1 && (i - M->HPDDM_ja[j] <= (i % bs)) && M->HPDDM_ja[j] != i)
          transpose[i].emplace_back(M->HPDDM_ja[j], M->HPDDM_a[j]);
      }
    int nnz = 0;
    for (int i = 0; i < transpose.size( ); ++i) {
      std::sort(transpose[i].begin( ), transpose[i].end( ),
                [](const std::pair< int, PetscScalar >& lhs,
                   const std::pair< int, PetscScalar >& rhs) { return lhs.first < rhs.first; });
      nnz += transpose[i].size( );
    }
    PetscInt* ia = new PetscInt[M->HPDDM_n / bs + 1];
    PetscInt* ja = new PetscInt[nnz / (bs * bs)];
    PetscScalar* a = new PetscScalar[nnz];
    ia[0] = 0;
    for (int i = 0; i < transpose.size( ); ++i) {
      for (int j = 0; j < transpose[i].size( ); ++j) {
        if (i % bs == 0 && j % bs == 0) ja[ia[i / bs] + j / bs] = transpose[i][j].first / bs;
        a[ia[i / bs] * (bs * bs) + j % bs + (j / bs) * (bs * bs) + (i % bs) * bs] =
          transpose[i][j].second;
      }
      if (i % bs == 0) ia[i / bs + 1] = ia[i / bs] + transpose[i].size( ) / bs;
    }
    MatSeqSBAIJSetPreallocationCSR(local, bs, ia, ja, a);
    MatSetOption(ptA->_petsc, MAT_NO_OFF_PROC_ENTRIES, PETSC_TRUE);
    MatSetOption(ptA->_petsc, MAT_SYMMETRIC, symmetric);
    MatAssemblyBegin(ptA->_petsc, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(ptA->_petsc, MAT_FINAL_ASSEMBLY);
    delete[] a;
    delete[] ja;
    delete[] ia;
    IS to, from;
    PetscInt nr;
    Vec rglobal;
    ISLocalToGlobalMappingGetSize(rmap, &nr);
    ISCreateStride(PETSC_COMM_SELF, nr, 0, 1, &to);
    ISLocalToGlobalMappingApplyIS(rmap, to, &from);
    MatCreateVecs(ptA->_petsc, &rglobal, NULL);
    Vec isVec;
    VecCreate(PETSC_COMM_SELF, &isVec);
    VecSetType(isVec, VECSTANDARD);
    VecSetSizes(isVec, PETSC_DECIDE, nr);
    VecScatterCreate(rglobal, from, isVec, to, &ptA->_scatter);
    VecDestroy(&isVec);
    VecDestroy(&rglobal);
    ISDestroy(&from);
    ISDestroy(&to);
    // ISLocalToGlobalMappingView(rmap, PETSC_VIEWER_STDOUT_WORLD);
    ISLocalToGlobalMappingDestroy(&rmap);
    if (verbosity > 0 && mpirank == 0)
      cout << " --- global CSR created (in " << MPI_Wtime( ) - timing << ")" << endl;
  }
  template< class Type, class K >
  long globalNumbering(Type* const& A, KN< K >* const& numbering) {
    if (A) {
      numbering->resize(A->_A->getDof( ));
      if (A->_num)
        for (int i = 0; i < numbering->n; ++i) numbering->operator[](i) = A->_num[i];
    }
    return 0L;
  }
  static Mat ff_to_PETSc(const HPDDM::MatrixCSR< PetscScalar >* const A) {
    Mat aux;
    if (A->HPDDM_sym) {
      std::vector< std::pair< int, int > >* transpose =
        new std::vector< std::pair< int, int > >[A->HPDDM_n]( );
      for (int i = 0; i < A->HPDDM_n; ++i)
        for (int j = A->HPDDM_ia[i] - (HPDDM_NUMBERING == 'F');
             j < A->HPDDM_ia[i + 1] - (HPDDM_NUMBERING == 'F'); ++j)
          transpose[A->HPDDM_ja[j] - (HPDDM_NUMBERING == 'F')].emplace_back(i, j);
      for (int i = 0; i < A->HPDDM_n; ++i)
        std::sort(transpose[i].begin( ), transpose[i].end( ));
      PetscInt* ia = new PetscInt[A->HPDDM_n + 1];
      PetscInt* ja = new PetscInt[A->HPDDM_nnz];
      PetscScalar* c = new PetscScalar[A->HPDDM_nnz];
      ia[0] = 0;
      for (int i = 0; i < A->HPDDM_n; ++i) {
        for (int j = 0; j < transpose[i].size( ); ++j) {
          c[ia[i] + j] = A->HPDDM_a[transpose[i][j].second];
          ja[ia[i] + j] = transpose[i][j].first;
        }
        ia[i + 1] = ia[i] + transpose[i].size( );
      }
      delete[] transpose;
      MatCreate(PETSC_COMM_SELF, &aux);
      MatSetSizes(aux, A->HPDDM_n, A->HPDDM_n, A->HPDDM_n, A->HPDDM_n);
      MatSetType(aux, MATSEQSBAIJ);
      MatSeqSBAIJSetPreallocationCSR(aux, 1, ia, ja, c);
      delete[] c;
      delete[] ja;
      delete[] ia;
    } else if(std::is_same<decltype(A->HPDDM_ia), PetscInt>::value)
      MatCreateSeqAIJWithArrays(PETSC_COMM_SELF, A->HPDDM_n, A->HPDDM_m, reinterpret_cast<PetscInt*>(A->HPDDM_ia), reinterpret_cast<PetscInt*>(A->HPDDM_ja), A->HPDDM_a, &aux);
    else {
      MatCreate(PETSC_COMM_SELF, &aux);
      MatSetSizes(aux, A->HPDDM_n, A->HPDDM_n, A->HPDDM_n, A->HPDDM_n);
      MatSetType(aux, MATSEQAIJ);
#if defined(PETSC_USE_64BIT_INDICES)
      PetscInt* ia = new PetscInt[A->HPDDM_n + 1];
      std::copy_n(A->HPDDM_ia, A->HPDDM_n + 1, ia);
      PetscInt* ja = new PetscInt[A->HPDDM_nnz];
      std::copy_n(A->HPDDM_ja, A->HPDDM_nnz, ja);
      PetscScalar* c = new PetscScalar[A->HPDDM_nnz];
      std::copy_n(A->HPDDM_a, A->HPDDM_nnz, c);
      MatSeqAIJSetPreallocationCSR(aux, ia, ja, c);
      delete[] c;
      delete[] ja;
      delete[] ia;
#else
      MatSeqAIJSetPreallocationCSR(aux, A->HPDDM_ia, A->HPDDM_ja, A->HPDDM_a);
#endif
    }
    return aux;
  }
  template< class Type >
  long ParMmgCommunicators(Type* const& A, KN< double >* const& gamma, KN< long >* const& rest, KN< KN< long >>* const& communicators) {
    if (A && A->_petsc && A->_A) {
      std::unordered_map<int, std::pair<int, PetscInt>> map;
      PetscScalar* val = new PetscScalar[A->_A->getDof()]();
      for(int i = 0; i < rest->n; ++i) {
        if(std::abs(gamma->operator[](i) - 1.0) < 1.0e-6) {
          map[rest->operator[](i)] = std::make_pair(i, A->_num[rest->operator[](i)]);
          val[rest->operator[](i)] = HPDDM::Wrapper<PetscScalar>::d__1;
        }
      }
      A->_A->recvBuffer(val);
      delete [] val;
      PetscScalar** const buffer = A->_A->getBuffer();
      const HPDDM::vectorNeighbor& neighbors = A->_A->getMap();
      communicators->resize(2 * neighbors.size() + 1);
      communicators->operator[](0).resize(neighbors.size());
      unsigned short k = 0;
      for(unsigned short i = 0; i < neighbors.size(); ++i) {
          communicators->operator[](0)[k] = neighbors[i].first;
          communicators->operator[](2 * k + 1).resize(neighbors[i].second.size());
          communicators->operator[](2 * k + 2).resize(neighbors[i].second.size());
          int m = 0;
          for(unsigned int j = 0; j < neighbors[i].second.size(); ++j) {
              if(std::abs(buffer[i][j] - HPDDM::Wrapper<PetscScalar>::d__1) < 1.0e-6) {
                  std::unordered_map<int, std::pair<int, PetscInt>>::const_iterator it = map.find(neighbors[i].second[j]);
                  if(it != map.cend()) {
                      communicators->operator[](2 * k + 1)[m] = it->second.first + 1;
                      communicators->operator[](2 * k + 2)[m++] = it->second.second + 1;
                  }
              }
          }
          if(m > 0) {
              communicators->operator[](2 * k + 1).resize(m);
              communicators->operator[](2 * k + 2).resize(m);
              ++k;
          }
      }
      int size;
      MPI_Comm_size(A->_A->getCommunicator(), &size);
      assert(size == 1 || k);
      communicators->operator[](0).resize(k);
      communicators->resize(2 * k + 1);
    }
    return 0L;
  }
  template< class Type >
  class changeOperator : public OneOperator {
   public:
    const int c;
    class changeOperator_Op : public E_F0mps {
     public:
      Expression A;
      Expression B;
      const int c;
      static const int n_name_param = 2;
      static basicAC_F0::name_and_type name_param[];
      Expression nargs[n_name_param];
      changeOperator_Op(const basicAC_F0& args, int d) : A(0), B(0), c(d) {
        args.SetNameParam(n_name_param, name_param, nargs);
        A = to< Type* >(args[0]);
        if (c == 0)
          B = to< Matrice_Creuse< upscaled_type<PetscScalar> >* >(args[1]);
        else
          B = to< Type* >(args[1]);
      }

      AnyType operator( )(Stack stack) const;
      operator aType( ) const { return atype< long >( ); }
    };
    E_F0* code(const basicAC_F0& args) const { return new changeOperator_Op(args, c); }
    changeOperator( )
      : OneOperator(atype< long >( ), atype< Type* >( ),
                    atype< Matrice_Creuse< upscaled_type<PetscScalar> >* >( )),
        c(0) {}
    changeOperator(int)
      : OneOperator(atype< long >( ), atype< Type* >( ), atype< Type* >( )), c(1) {}
  };
  template< class Type >
  basicAC_F0::name_and_type changeOperator< Type >::changeOperator_Op::name_param[] = {
    {"restriction", &typeid(Matrice_Creuse< double >*)}, {"parent", &typeid(Type*)}};
  template< class Type >
  void change(Type* const& ptA, Matrice_Creuse< upscaled_type<PetscScalar> >* const& mat, Type* const& ptB,
              Matrice_Creuse< double >* const& pList, Type* const& ptParent) {
    if (mat) {
      if (ptA) {
        MatriceMorse< upscaled_type<PetscScalar> >* mN = nullptr;
        if (mat->A) mN = static_cast< MatriceMorse< upscaled_type<PetscScalar> >* >(&*(mat->A));
        PetscBool assembled = PETSC_FALSE;
        if (ptA->_petsc) MatAssembled(ptA->_petsc, &assembled);
        if (mN) {
          HPDDM::MatrixCSR< void >* dL = nullptr;
          if (pList && pList->A) {
            MatriceMorse< double >* mList = static_cast< MatriceMorse< double >* >(&*(pList->A));
            ffassert(mList->n == mList->nnz);
            ffassert(mList->m == mN->n);
            dL = new_HPDDM_MatrixCSRvoid(mList, false);
          }
          HPDDM::MatrixCSR< PetscScalar >* dM = new_HPDDM_MatrixCSR< PetscScalar >(mN);
          HPDDM::MatrixCSR< PetscScalar >* dN;
          if (!dL)
            dN = dM;
          else {
            unsigned int* perm = new unsigned int[dM->HPDDM_n]( );
            for (unsigned int i = 0; i < dL->HPDDM_n; ++i) perm[dL->HPDDM_ja[i]] = i + 1;
            dN = new HPDDM::MatrixCSR< PetscScalar >(dM, dL, perm);
            delete[] perm;
            delete dM;
            dM = nullptr;
          }
          if (ptA->_A) {
              ptA->_A->setMatrix(dN);
#if defined(PCHPDDM) && defined(PETSC_USE_SHARED_LIBRARIES)
              PC pc = nullptr;
              if(ptParent) {
                  PetscInt M, N;
                  Mat** mat;
                  MatNestGetSubMats(ptParent->_petsc, &M, &N, &mat);
                  PetscInt i;
                  for(i = 0; i < std::min(N, M); ++i) {
                      if(mat[i][i] == ptA->_petsc)
                          break;
                  }
                  if(i < std::min(N, M)) {
                      PC parent;
                      KSPGetPC(ptParent->_ksp, &parent);
                      KSP *subksp;
                      PetscInt nsplits;
                      PCFieldSplitGetSubKSP(parent, &nsplits, &subksp);
                      KSPGetPC(subksp[i], &pc);
                      PetscFree(subksp);
                  }
              }
              else if(ptA->_ksp) {
                  KSPGetPC(ptA->_ksp, &pc);
              }
              if(pc) {
                  PCType type;
                  PCGetType(pc, &type);
                  PetscBool isType;
                  PetscStrcmp(type, PCHPDDM, &isType);
                  if(isType) {
                      const HPDDM::MatrixCSR<PetscScalar>* const A = ptA->_A->getMatrix();
                      Mat aux = ff_to_PETSc(A);
                      Mat N;
                      PetscObjectQuery((PetscObject)pc, "_PCHPDDM_Neumann_Mat", (PetscObject*)&N);
                      if(!N) {
                          PetscInt* idx;
                          PetscMalloc1(dN->HPDDM_n, &idx);
                          std::copy_n(ptA->_num, dN->HPDDM_n, idx);
                          IS is;
                          ISCreateGeneral(PETSC_COMM_SELF, ptA->_A->getMatrix()->HPDDM_n, idx, PETSC_OWN_POINTER, &is);
                          PetscObjectCompose((PetscObject)pc, "_PCHPDDM_Neumann_Mat", (PetscObject)aux);
                          PCHPDDMSetAuxiliaryMat(pc, is, aux, NULL, NULL);
                          PCSetFromOptions(pc);
                          MatDestroy(&aux);
                          ISDestroy(&is);
                      }
                      else {
                          PetscObjectCompose((PetscObject)pc, "_PCHPDDM_Neumann_Mat", (PetscObject)aux);
                          PCHPDDMSetAuxiliaryMat(pc, NULL, aux, NULL, NULL);
                          MatDestroy(&aux);
                      }
                  }
              }
#endif
          }
          ffassert(ptA->_num);
          PetscInt* ia = nullptr;
          PetscInt* ja = nullptr;
          PetscScalar* c = nullptr;
          bool free = true;
          if (ptA->_cnum) {
            PetscInt* cnum = ptA->_num + dN->HPDDM_n;
            if (ptA->_petsc) {
                int rank, size, N;
                MPI_Comm_size(PetscObjectComm((PetscObject)ptA->_petsc), &size);
                if (size > 1) {
                  const PetscInt* ranges;
                  MatGetOwnershipRangesColumn(ptA->_petsc, &ranges);
                  MPI_Comm_rank(PetscObjectComm((PetscObject)ptA->_petsc), &rank);
                  MatGetSize(ptA->_petsc, NULL, &N);
                  if (ranges[1] == N) {
                      if (rank != 0)
                          cnum = new PetscInt[N];
                      MPI_Bcast(cnum, N, MPIU_INT, 0, PetscObjectComm((PetscObject)ptA->_petsc));
                  }
                }
            }
            free = HPDDM::template Subdomain< PetscScalar >::distributedCSR(
              ptA->_num, ptA->_first, ptA->_last, ia, ja, c, dN, cnum);
            if (cnum != ptA->_num + dN->HPDDM_n)
                delete [] cnum;
          }
          else
            free = ptA->_A->distributedCSR(ptA->_num, ptA->_first, ptA->_last, ia, ja, c);
          if (assembled) {
            MatZeroEntries(ptA->_petsc);
            for (PetscInt i = 0; i < ptA->_last - ptA->_first; ++i) {
              PetscInt row = ptA->_first + i;
              MatSetValues(ptA->_petsc, 1, &row, ia[i + 1] - ia[i],
                           ja + ia[i], c + ia[i], INSERT_VALUES);
            }
            MatAssemblyBegin(ptA->_petsc, MAT_FINAL_ASSEMBLY);
            MatAssemblyEnd(ptA->_petsc, MAT_FINAL_ASSEMBLY);
          } else {
            if (!ptA->_petsc) {
              MatCreate(PETSC_COMM_WORLD, &ptA->_petsc);
              MatSetSizes(ptA->_petsc, ptA->_last - ptA->_first, ptA->_last - ptA->_first,
                          PETSC_DECIDE, PETSC_DECIDE);
            }
            if (ptA->_A && ptA->_A->getMatrix( )->HPDDM_sym) {
              MatSetType(ptA->_petsc, MATSBAIJ);
              MatSeqSBAIJSetPreallocationCSR(ptA->_petsc, 1, ia, ja, c);
              MatMPISBAIJSetPreallocationCSR(ptA->_petsc, 1, ia, ja, c);
            } else {
              MatSetType(ptA->_petsc, MATAIJ);
              MatSeqAIJSetPreallocationCSR(ptA->_petsc, ia, ja, c);
              MatMPIAIJSetPreallocationCSR(ptA->_petsc, ia, ja, c);
            }
          }
          if (free) {
            delete[] ia;
            delete[] ja;
            delete[] c;
          }
        } else if (!assembled) {
          PetscInt m;
          MatGetLocalSize(ptA->_petsc, &m, NULL);
          PetscInt* ia = new PetscInt[m + 1]( );
          MatSeqAIJSetPreallocationCSR(ptA->_petsc, ia, NULL, NULL);
          MatMPIAIJSetPreallocationCSR(ptA->_petsc, ia, NULL, NULL);
          delete[] ia;
        } else {
          MatAssemblyBegin(ptA->_petsc, MAT_FINAL_ASSEMBLY);
          MatAssemblyEnd(ptA->_petsc, MAT_FINAL_ASSEMBLY);
        }
        if (ptParent) {
          PetscBool assembled;
          MatAssembled(ptParent->_petsc, &assembled);
          if (!assembled) {
            PetscInt M, N;
            Mat** mat;
            MatNestGetSubMats(ptParent->_petsc, &M, &N, &mat);
            PetscBool assemble = PETSC_TRUE;
            for (PetscInt i = 0; i < M && assemble; ++i) {
              for (PetscInt j = 0; j < N && assemble; ++j) {
                if (mat[i][j]) {
                  PetscBool assembled;
                  MatAssembled(mat[i][j], &assembled);
                  if (!assembled) assemble = PETSC_FALSE;
                }
              }
            }
            if (assemble) {
              MatAssemblyBegin(ptParent->_petsc, MAT_FINAL_ASSEMBLY);
              MatAssemblyEnd(ptParent->_petsc, MAT_FINAL_ASSEMBLY);
            }
          }
        }
        if (ptA->_ksp) {
          KSPSetOperators(ptA->_ksp, ptA->_petsc, ptA->_petsc);
          if (std::is_same< Type, Dmat >::value && ptA->_vS && !ptA->_vS->empty( )) {
            KSPSetFromOptions(ptA->_ksp);
            PetscBool assembled;
            MatAssembled(ptA->_petsc, &assembled);
            if (assembled) {
              if (ptParent) MatAssembled(ptParent->_petsc, &assembled);
              if (assembled) KSPSetUp(ptA->_ksp);
            }
            PC pc;
            KSPGetPC(ptA->_ksp, &pc);
            PCType type;
            PCGetType(pc, &type);
            PetscBool isFieldSplit;
            PetscStrcmp(type, PCFIELDSPLIT, &isFieldSplit);
            if (isFieldSplit) {
              PetscInt nsplits;
              KSP* subksp;
              PCFieldSplitGetSubKSP(pc, &nsplits, &subksp);
              for (int i = 0; i < nsplits; ++i) {
                PC subpc;
                KSPGetPC(subksp[i], &subpc);
                PCGetType(subpc, &type);
                PetscStrcmp(type, PCFIELDSPLIT, &isFieldSplit);
                if (isFieldSplit) {
                  pc = subpc;
                  break;
                }
              }
              setCompositePC(pc, ptA->_vS);
            }
          }
        }
      }
    } else {
      MatType type;
      PetscBool isType = PETSC_FALSE;
      if (ptB->_petsc) {
        MatGetType(ptB->_petsc, &type);
        PetscStrcmp(type, MATNEST, &isType);
      }
      if (!isType) {
        Mat backup = ptA->_petsc;
        ptA->_petsc = nullptr;
        ptA->dtor( );
        ptA->_petsc = backup;
        PetscMPIInt flag;
        if (ptA->_petsc)
          MPI_Comm_compare(PetscObjectComm((PetscObject)ptA->_petsc),
                           PetscObjectComm((PetscObject)ptB->_petsc), &flag);
        if (!ptA->_petsc || (flag != MPI_CONGRUENT && flag != MPI_IDENT)) {
          MatDestroy(&ptA->_petsc);
          if (ptB->_petsc) {
            PetscBool assembled;
            MatAssembled(ptB->_petsc, &assembled);
            MatDuplicate(ptB->_petsc, assembled ? MAT_COPY_VALUES : MAT_DO_NOT_COPY_VALUES, &ptA->_petsc);
            MatDestroy(&ptB->_petsc);
          }
        } else
          MatHeaderReplace(ptA->_petsc, &ptB->_petsc);
        KSPDestroy(&ptB->_ksp);
        ptA->_A = ptB->_A;
        ptB->_A = nullptr;
        ptA->_num = ptB->_num;
        ptA->_cnum = ptB->_cnum;
        ptA->_first = ptB->_first;
        ptA->_cfirst = ptB->_cfirst;
        ptA->_last = ptB->_last;
        ptA->_clast = ptB->_clast;
        ptB->_num = nullptr;
        ptA->_exchange = ptB->_exchange;
        if (ptB->_exchange) {
          ptA->_exchange[1] = ptB->_exchange[1];
          ptB->_exchange = nullptr;
        }
        ptA->_D = ptB->_D;
        ptB->_D = nullptr;
      }
    }
  }
  template< class Type >
  AnyType changeOperator< Type >::changeOperator_Op::operator( )(Stack stack) const {
    Type* ptA = GetAny< Type* >((*A)(stack));
    Matrice_Creuse< upscaled_type<PetscScalar> >* mat =
      c == 0 ? GetAny< Matrice_Creuse< upscaled_type<PetscScalar> >* >((*B)(stack)) : nullptr;
    Type* ptB = c != 0 ? GetAny< Type* >((*B)(stack)) : nullptr;
    Matrice_Creuse< double >* pList =
      nargs[0] ? GetAny< Matrice_Creuse< double >* >((*nargs[0])(stack)) : nullptr;
    Type* ptParent = nargs[1] ? GetAny< Type* >((*nargs[1])(stack)) : nullptr;
    change(ptA, mat, ptB, pList, ptParent);
    return 0L;
  }
  Dmat* changeOperatorSimple(Dmat* const& dA, Dmat* const& A) {
    Dmat* const null = nullptr;
    change(dA, nullptr, A, nullptr, null);
    return dA;
  }
  Dmat* changeOperatorSimple(Dmat* const& dA, Matrice_Creuse< upscaled_type<PetscScalar> >* const& A) {
    Dmat* const null = nullptr;
    change(dA, A, null, nullptr, null);
    return dA;
  }
  template<class A, class B>
  struct scale {
    using first_argument_type  = A;
    using second_argument_type = B;
    using result_type          = long;
    static long f(A const& a, const B& b) {
      MatScale(a->_petsc, PetscScalar(b));
      return 0L;
    }
  };
  template<class A, class B>
  struct AXPY {
    using first_argument_type  = A;
    using second_argument_type = B;
    using result_type          = long;
    static long f(A const& a, const B& b) {
      MatAXPY(a->_petsc, 1.0, b->_petsc, DIFFERENT_NONZERO_PATTERN);
      return 0L;
    }
  };
  std::pair<PetscScalar, Dmat*>* toP(Dmat* const& M) {
    return new std::pair<PetscScalar, Dmat*>{1.0, M};
  }
  AnyType M2P(Stack, const AnyType& p) {
    return toP(GetAny<Dmat*>(p));
  }
  template<class K>
  AnyType AddCombDmat(Stack stack, Expression emat, Expression combMat) {
    Dmat* A = GetAny<Dmat*>((*emat)(stack));
    ffassert(A && A->_petsc);
    std::pair<upscaled_type<PetscScalar>, Dmat*>* p = GetAny<std::pair<upscaled_type<PetscScalar>, Dmat*>*>((*combMat)(stack));
    ffassert(p && p->second && p->second->_petsc);
    MatAXPY(A->_petsc, PetscScalar(p->first), p->second->_petsc, DIFFERENT_NONZERO_PATTERN);
    delete p;
    return A;
  }
  template<class K>
  struct Op2 {
    using first_argument_type  = K;
    using second_argument_type = Dmat*;
    using result_type          = std::pair<K, Dmat*>*;
    typedef std::pair<K, Dmat*> P;
    static P* f(K const& a,Dmat* const& b) {
      return new P{a, b};
    }
  };
  long changeSchur(Dmat* const& dA, KN< Matrice_Creuse< upscaled_type<PetscScalar> > >* const& schurPreconditioner,
                   KN< double >* const& schurList) {
    setVectorSchur(dA, schurPreconditioner, schurList);
    return 0L;
  }
  template< class Type, class K >
  long originalNumbering(Type* const& A, KN< K >* const& in, KN< long >* const& interface) {
    if (A) A->_A->originalNumbering(STL< long >(*interface), *in);
    return 0L;
  }
  template< class K >
  long renumber(KN< K >* const& in, KN< long >* const& numbering, KN< K >* const& out) {
    PetscInt low;
    out->resize(in->n);
    Vec x, y;
    VecCreateMPIWithArray(PETSC_COMM_WORLD, 1, in->n, PETSC_DECIDE,
                          in->operator PetscScalar*(), &x);
    VecGetOwnershipRange(x, &low, NULL);
    VecCreateMPIWithArray(PETSC_COMM_WORLD, 1, in->n, PETSC_DECIDE,
                          out->operator PetscScalar*(), &y);

    IS to;
    PetscInt* idx = new PetscInt[in->n];
    for (int i = 0; i < in->n; ++i) idx[i] = numbering->operator[](i);
    ISCreateGeneral(PETSC_COMM_SELF, in->n, idx, PETSC_USE_POINTER, &to);
    VecScatter scatter;
    VecScatterCreate(x, NULL, y, to, &scatter);
    VecScatterBegin(scatter, x, y, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(scatter, x, y, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterDestroy(&scatter);
    ISDestroy(&to);
    delete[] idx;
    VecDestroy(&x);
    VecDestroy(&y);
    return 0L;
  }
  long stagePush(string* const& str) {
#if defined(PETSC_USE_LOG)
      PetscLogStage stage;
      PetscLogStageGetId(str->c_str(), &stage);
      if(stage == -1) {
        PetscLogStageRegister(str->c_str(), &stage);
      }
      PetscLogStagePush(stage);
#endif
      return 0L;
  }
  long stagePop() {
#if defined(PETSC_USE_LOG)
      PetscLogStagePop();
#endif
      return 0L;
  }
  double memoryGetCurrentUsage() {
      PetscLogDouble mem;
      PetscMemoryGetCurrentUsage(&mem);
      return mem;
  }
  long hasType(string* const& obj, string* const& type) {
      if (obj->size() > 0 && type->size() > 0) {
          std::string o(*obj);
          std::transform(o.begin(), o.end(), o.begin(), ::toupper);
          std::string t(*type);
          std::transform(t.begin(), t.end(), t.begin(), ::tolower);
          if (o.compare("PC") == 0) {
              PetscErrorCode (*pc)(PC*);
              PetscFunctionListFind(PCList, t.c_str(), &pc);
              if(pc)
                  return 1L;
          }
          else if (o.compare("KSP") == 0) {
              PetscErrorCode (*ksp)(KSP*);
              PetscFunctionListFind(KSPList, t.c_str(), &ksp);
              if(ksp)
                  return 1L;
          }
          else if (o.compare("MAT") == 0) {
              PetscErrorCode (*mat)(Mat*);
              PetscFunctionListFind(MatList, t.c_str(), &mat);
              if(mat)
                  return 1L;
          }
          else if (o.compare("MATSOLVER") == 0) {
              PetscBool foundtype;
              MatSolverTypeGet(t.c_str(), MATAIJ, MAT_FACTOR_LU, &foundtype, NULL, NULL);
              if(foundtype)
                  return 1L;
          }
      }
      return 0L;
  }

  template< class HpddmType, int C >
  class initCSRfromDMatrix : public OneOperator {
   public:
    const int c;
    class initCSRfromDMatrix_Op : public E_F0mps {
     public:
      Expression A;
      Expression B;
      Expression K;
      const OneOperator *codeA;
      const int c;
      static const int n_name_param = 3;
      static basicAC_F0::name_and_type name_param[];
      Expression nargs[n_name_param];
      initCSRfromDMatrix_Op(const basicAC_F0& args, int d)
        : A(0), B(0), K(0), codeA(nullptr), c(d) {
        args.SetNameParam(n_name_param, name_param, nargs);
        A = to< DistributedCSR< HpddmType >* >(args[0]);
        B = to< DistributedCSR< HpddmType >* >(args[1]);
        if (c == 0) K = to< Matrice_Creuse< upscaled_type<PetscScalar> >* >(args[2]);
        else if (c == 1) {
          const Polymorphic* op = dynamic_cast< const Polymorphic* >(args[2].LeftValue( ));
          ffassert(op);
          codeA = op->Find("(", ArrayOfaType(atype< KN< upscaled_type<PetscScalar> >* >( ), false));
        }
      }

      AnyType operator( )(Stack stack) const;
      operator aType( ) const { return atype< DistributedCSR< HpddmType >* >( ); }
    };
    E_F0* code(const basicAC_F0& args) const {
      return new initCSRfromDMatrix_Op(args, c);
    }
    initCSRfromDMatrix( )
      : OneOperator(
          atype< DistributedCSR< HpddmType >* >( ), atype< DistributedCSR< HpddmType >* >( ),
          atype< DistributedCSR< HpddmType >* >( ), atype< typename std::conditional<C == 0, Matrice_Creuse< upscaled_type<PetscScalar> >*, Polymorphic*>::type >( )),
        c(C) {}
    initCSRfromDMatrix(int)
      : OneOperator(
          atype< DistributedCSR< HpddmType >* >( ), atype< DistributedCSR< HpddmType >* >( ),
          atype< DistributedCSR< HpddmType >* >( )),
        c(2) {}
  };
  template< class HpddmType, int C >
  basicAC_F0::name_and_type
    initCSRfromDMatrix< HpddmType, C >::initCSRfromDMatrix_Op::name_param[] = {
      {C == 1 ? "transpose" : "clean", C == 1 ? &typeid(Polymorphic*) : &typeid(bool)},
      {"symmetric", &typeid(bool)},
      {"restriction", &typeid(Matrice_Creuse< double >*)}};

  template< class HpddmType, int D >
  class initRectangularCSRfromDMatrix : public OneOperator {
   public:
    const int c;
    class initRectangularCSRfromDMatrix_Op : public E_F0mps {
     public:
      Expression A;
      Expression B;
      Expression C;
      Expression K;
      const OneOperator *codeA;
      const int c;
      static const int n_name_param = 2;
      static basicAC_F0::name_and_type name_param[];
      Expression nargs[n_name_param];
      initRectangularCSRfromDMatrix_Op(const basicAC_F0& args, int d)
        : A(0), B(0), C(0), K(0), codeA(nullptr), c(d) {
        args.SetNameParam(n_name_param, name_param, nargs);
        A = to< DistributedCSR< HpddmType >* >(args[0]);
        B = to< DistributedCSR< HpddmType >* >(args[1]);
        C = to< DistributedCSR< HpddmType >* >(args[2]);
        if (c == 0 || c == 3) K = to< Matrice_Creuse< upscaled_type<PetscScalar> >* >(args[3]);
        if (c == 2) {
          const Polymorphic* op = dynamic_cast< const Polymorphic* >(args[3].LeftValue( ));
          ffassert(op);
          codeA = op->Find("(", ArrayOfaType(atype< KN< upscaled_type<PetscScalar> >* >( ), false));
        }
      }

      AnyType operator( )(Stack stack) const;
      operator aType( ) const { return c != 3 ? atype< DistributedCSR< HpddmType >* >( ) : atype< long >( ); }
    };
    E_F0* code(const basicAC_F0& args) const {
      return new initRectangularCSRfromDMatrix_Op(args, c);
    }
    initRectangularCSRfromDMatrix( )
      : OneOperator(
          atype< DistributedCSR< HpddmType >* >( ), atype< DistributedCSR< HpddmType >* >( ),
          atype< DistributedCSR< HpddmType >* >( ), atype< DistributedCSR< HpddmType >* >( ),
          atype< typename std::conditional<D == 0, Matrice_Creuse< upscaled_type<PetscScalar> >*, Polymorphic* >::type >( )),
        c(D ? 2 : 0) {}
    initRectangularCSRfromDMatrix(int)
      : OneOperator(
          atype< DistributedCSR< HpddmType >* >( ), atype< DistributedCSR< HpddmType >* >( ),
          atype< DistributedCSR< HpddmType >* >( ), atype< DistributedCSR< HpddmType >* >( )),
        c(1) {}
    initRectangularCSRfromDMatrix(int, int)
      : OneOperator(
          atype< long >( ), atype< DistributedCSR< HpddmType >* >( ),
          atype< DistributedCSR< HpddmType >* >( ), atype< DistributedCSR< HpddmType >* >( ),
          atype< Matrice_Creuse< upscaled_type<PetscScalar> >* >( )),
        c(3) {}
  };
  template< class HpddmType, int D >
  basicAC_F0::name_and_type
    initRectangularCSRfromDMatrix< HpddmType, D >::initRectangularCSRfromDMatrix_Op::name_param[] = {
      {D == 1 ? "transpose" : "clean", D == 1 ? &typeid(Polymorphic*) : &typeid(bool)},
      {"numbering", &typeid(KN<upscaled_type<PetscScalar>>*)}};

#if !defined(PETSC_USE_REAL_SINGLE)
  template< class HpddmType >
  class initCSRfromMatrix_Op : public E_F0mps {
   public:
    Expression A;
    Expression K;
    Expression size;
    static const int n_name_param = 6;
    static basicAC_F0::name_and_type name_param[];
    Expression nargs[n_name_param];
    initCSRfromMatrix_Op(const basicAC_F0& args, Expression param1, Expression param2,
                         Expression param3)
      : A(param1), K(param2), size(param3) {
      args.SetNameParam(n_name_param, name_param, nargs);
    }

    AnyType operator( )(Stack stack) const;
  };
  template< class HpddmType >
  basicAC_F0::name_and_type initCSRfromMatrix_Op< HpddmType >::name_param[] = {
    {"communicator", &typeid(pcommworld)},
    {"bs", &typeid(long)},
    {"symmetric", &typeid(bool)},
    {"clean", &typeid(bool)},
    {"bsr", &typeid(bool)},
    {"prune", &typeid(bool)}};
  template< class HpddmType >
  class initCSRfromMatrix : public OneOperator {
   public:
    initCSRfromMatrix( )
      : OneOperator(atype< DistributedCSR< HpddmType >* >( ),
                    atype< DistributedCSR< HpddmType >* >( ),
                    atype< Matrice_Creuse< PetscScalar >* >( ), atype< KN< long >* >( )) {}

    E_F0* code(const basicAC_F0& args) const {
      return new initCSRfromMatrix_Op< HpddmType >(args, t[0]->CastTo(args[0]),
                                                   t[1]->CastTo(args[1]), t[2]->CastTo(args[2]));
    }
  };
  template< class HpddmType >
  AnyType initCSRfromMatrix_Op< HpddmType >::operator( )(Stack stack) const {
    DistributedCSR< HpddmType >* ptA = GetAny< DistributedCSR< HpddmType >* >((*A)(stack));
    Matrice_Creuse< PetscScalar >* ptK = GetAny< Matrice_Creuse< PetscScalar >* >((*K)(stack));
    KN< long >* ptSize = GetAny< KN< long >* >((*size)(stack));
    MatriceMorse< PetscScalar >* mK = static_cast< MatriceMorse< PetscScalar >* >(&(*(ptK->A)));
    PetscInt bs = nargs[1] ? GetAny< long >((*nargs[1])(stack)) : 1;
    MatCreate(nargs[0] ? *static_cast< MPI_Comm* >(GetAny< pcommworld >((*nargs[0])(stack))) : PETSC_COMM_WORLD, &ptA->_petsc);
    if (bs > 1) MatSetBlockSize(ptA->_petsc, bs);
    bool bsr = nargs[4] && GetAny< bool >((*nargs[4])(stack));
    bool prune = bsr ? (nargs[5] && GetAny< bool >((*nargs[5])(stack))) : false;
    if (prune) {
      ffassert(0);
    }
    if (mpisize > 1) {
      ffassert(ptSize->n >= 3 + mK->n);
      ptA->_first = ptSize->operator( )(0);
      ptA->_last = ptSize->operator( )(1);
      MatSetSizes(ptA->_petsc, mK->n * (bsr ? bs : 1), mK->n * (bsr ? bs : 1),
                  ptSize->operator( )(2) * (bsr ? bs : 1), ptSize->operator( )(2) * (bsr ? bs : 1));
      if (prune) {
        ffassert(0);
      }
      ptA->_num = new PetscInt[ptSize->n - 3];
      for (int i = 3; i < ptSize->n; ++i) ptA->_num[i - 3] = ptSize->operator( )(i);
    } else {
      ptA->_first = 0;
      ptA->_last = mK->n;
      ptA->_num = nullptr;
      MatSetSizes(ptA->_petsc, mK->n * (bsr ? bs : 1), mK->n * (bsr ? bs : 1),
                  mK->n * (bsr ? bs : 1), mK->n * (bsr ? bs : 1));
      if (prune) {
        ffassert(0);
      }
    }
    bool clean = nargs[3] && GetAny< bool >((*nargs[3])(stack));
    if (clean) ptSize->resize(0);
    MatSetType(ptA->_petsc, bsr ? MATBAIJ : MATAIJ);
    mK->CSR( );
    if (bsr) {
      MatSeqBAIJSetPreallocationCSR(ptA->_petsc, bs, reinterpret_cast< PetscInt* >(mK->p),
                                    reinterpret_cast< PetscInt* >(mK->j), mK->aij);
      MatMPIBAIJSetPreallocationCSR(ptA->_petsc, bs, reinterpret_cast< PetscInt* >(mK->p),
                                    reinterpret_cast< PetscInt* >(mK->j), mK->aij);
    }
    else {
      MatSeqAIJSetPreallocationCSR(ptA->_petsc, reinterpret_cast< PetscInt* >(mK->p),
                                   reinterpret_cast< PetscInt* >(mK->j), mK->aij);
      MatMPIAIJSetPreallocationCSR(ptA->_petsc, reinterpret_cast< PetscInt* >(mK->p),
                                   reinterpret_cast< PetscInt* >(mK->j), mK->aij);
    }
    if (prune) {
      ffassert(0);
    }
    MatSetOption(
      ptA->_petsc, MAT_SYMMETRIC,
      nargs[2] && GetAny< bool >((*nargs[2])(stack)) ? PETSC_TRUE : PETSC_FALSE);
    if (clean) ptK->destroy( );
    if (prune) {
      ffassert(0);
    }
    return ptA;
  }
  template< class HpddmType >
  class initCSRfromArray_Op : public E_F0mps {
   public:
    Expression A;
    Expression K;
    static const int n_name_param = 5;
    static basicAC_F0::name_and_type name_param[];
    Expression nargs[n_name_param];
    initCSRfromArray_Op(const basicAC_F0& args, Expression param1, Expression param2)
      : A(param1), K(param2) {
      args.SetNameParam(n_name_param, name_param, nargs);
    }

    AnyType operator( )(Stack stack) const;
  };
  template< class HpddmType >
  basicAC_F0::name_and_type initCSRfromArray_Op< HpddmType >::name_param[] = {
    {"columns", &typeid(KN< long >*)},
    {"communicator", &typeid(pcommworld)},
    {"bs", &typeid(long)},
    {"symmetric", &typeid(bool)},
    {"clean", &typeid(bool)}};
  template< class HpddmType >
  class initCSRfromArray : public OneOperator {
   public:
    initCSRfromArray( )
      : OneOperator(atype< DistributedCSR< HpddmType >* >( ),
                    atype< DistributedCSR< HpddmType >* >( ),
                    atype< KN< Matrice_Creuse< PetscScalar > >* >( )) {}

    E_F0* code(const basicAC_F0& args) const {
      return new initCSRfromArray_Op< HpddmType >(args, t[0]->CastTo(args[0]),
                                                  t[1]->CastTo(args[1]));
    }
  };
  template< class HpddmType >
  AnyType initCSRfromArray_Op< HpddmType >::operator( )(Stack stack) const {
    MPI_Comm comm = nargs[1] ? *static_cast< MPI_Comm* >(GetAny< pcommworld >((*nargs[1])(stack))) : PETSC_COMM_WORLD;
    int size;
    MPI_Comm_size(comm, &size);
    DistributedCSR< HpddmType >* ptA = GetAny< DistributedCSR< HpddmType >* >((*A)(stack));
    KN< Matrice_Creuse< PetscScalar > >* ptK =
      GetAny< KN< Matrice_Creuse< PetscScalar > >* >((*K)(stack));
    KN< long >* ptJ = nargs[0] ? GetAny< KN< long >* >((*nargs[0])(stack)) : nullptr;
    if (!ptJ || ptJ->n == ptK->n) {
      double timing = MPI_Wtime( );
      std::vector< std::pair< int, int > > v;
      if (ptJ) {
        v.reserve(ptJ->n);
        for (int i = 0; i < ptJ->n; ++i) v.emplace_back(std::make_pair(ptJ->operator[](i), i));
        std::sort(v.begin( ), v.end( ));
      }
      PetscInt bs = nargs[2] ? GetAny< long >((*nargs[2])(stack)) : 1;
      int rank;
      MPI_Comm_rank(comm, &rank);
      int* dims = new int[size]( );
      std::vector< std::pair< int, int > >::const_iterator it =
        std::lower_bound(v.cbegin( ), v.cend( ), std::make_pair(rank, 0),
                         [](const std::pair< int, int >& lhs, const std::pair< int, int >& rhs) {
                           return lhs.first < rhs.first;
                         });
      int n = 0;
      if (!ptJ) {
        dims[rank] = ptK->operator[](mpisize / ptK->n > 1 ? 0 : rank).M( );
        n = ptK->operator[](mpisize / ptK->n > 1 ? 0 : rank).N( );
      } else if (it->first == rank) {
        dims[rank] = ptK->operator[](it->second).M( );
        n = ptK->operator[](it->second).N( );
      }
      MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, dims, 1, MPI_INT, comm);
      if (ptJ)
        size = v.size( );
      else
        size = ptK->n;
      PetscInt* ia = new PetscInt[n + 1];
      for (int i = 0; i < n + 1; ++i) {
        ia[i] = 0;
        for (int k = 0; k < size; ++k) {
          MatriceMorse< PetscScalar >* mA =
            static_cast< MatriceMorse< PetscScalar >* >(&(*(ptK->operator[](k)).A));
          mA->CSR( );
          if (i < mA->n + 1) ia[i] += mA->p[i];
        }
      }
      PetscInt* ja = new PetscInt[ia[n]];
      PetscScalar* c = new PetscScalar[ia[n]];
      int nnz = 0;
      int offset = (ptJ ? std::accumulate(dims, dims + v.begin( )->first, 0) : 0);
      for (int i = 0; i < n; ++i) {
        int backup = offset;
        for (int k = 0; k < size; ++k) {
          MatriceMorse< PetscScalar >* mA = static_cast< MatriceMorse< PetscScalar >* >(
            &(*(ptK->operator[](ptJ ? v[k].second : k)).A));
          mA->CSR( );
          if (i < mA->n) {
            int j = mA->p[i];
            for (; j < mA->p[i + 1] && mA->p[j] < dims[ptJ ? v[k].first : k]; ++j)
              ja[nnz + j - mA->p[i]] = mA->j[j] + offset;
            std::copy_n(mA->aij + mA->p[i], j - mA->p[i], c + nnz);
            nnz += j - mA->p[i];
            if (k < (ptJ ? v.size( ) : size) - 1)
              offset +=
                (ptJ ? std::accumulate(dims + v[k].first, dims + v[k + 1].first, 0) : dims[k]);
          }
        }
        offset = backup;
      }
      delete[] dims;
      MatCreate(comm, &ptA->_petsc);
      if (bs > 1) MatSetBlockSize(ptA->_petsc, bs);
      {
        int size;
        MPI_Comm_size(comm, &size);
        if (!ptJ && (size / ptK->n > 1))
          MatSetSizes(ptA->_petsc, n, n, PETSC_DECIDE, PETSC_DECIDE);
        else
          MatSetSizes(ptA->_petsc, n, ptK->operator[](ptJ ? it->second : rank).M( ), PETSC_DECIDE,
                      PETSC_DECIDE);
      }
      ptA->_first = 0;
      ptA->_last = n;
      if (nargs[4] && GetAny< bool >((*nargs[4])(stack))) {
        int* cl = nullptr;
        int* lg = nullptr;
        PetscScalar* a = nullptr;
        for (int k = 0; k < size; ++k) {
          MatriceMorse< PetscScalar >* mA = static_cast< MatriceMorse< PetscScalar >* >(
            &(*(ptK->operator[](ptJ ? v[k].second : k)).A));
          if (k == 0) {
            mA->CSR( );
            cl = mA->j;
            lg = mA->p;
            a = mA->aij;
          } else {
            if (mA->j == cl) mA->j = nullptr;
            if (mA->p == lg) mA->p = nullptr;
            if (mA->aij == a) mA->aij = nullptr;
          }
        }
        ptK->resize(0);
      }
      MatSetType(ptA->_petsc, MATAIJ);
      MatSeqAIJSetPreallocationCSR(ptA->_petsc, ia, ja, c);
      MatMPIAIJSetPreallocationCSR(ptA->_petsc, ia, ja, c);
      MatSetOption(
        ptA->_petsc, MAT_SYMMETRIC,
        nargs[3] && GetAny< bool >((*nargs[3])(stack)) ? PETSC_TRUE : PETSC_FALSE);
      delete[] ia;
      delete[] ja;
      delete[] c;
      if (verbosity > 0 && mpirank == 0)
        cout << " --- global CSR created (in " << MPI_Wtime( ) - timing << ")" << endl;
      timing = MPI_Wtime( );
    }
    return ptA;
  }
#endif

  template< class HpddmType, bool C = false >
  class initCSR : public OneOperator {
   public:
    const int c;
    class E_initCSR : public E_F0mps {
     public:
      Expression A;
      Expression K;
      Expression R;
      Expression D;
      const int c;
      static const int n_name_param = 6;
      static basicAC_F0::name_and_type name_param[];
      Expression nargs[n_name_param];
      E_initCSR(const basicAC_F0& args, int d) : A(0), K(0), R(0), D(0), c(d) {
        args.SetNameParam(n_name_param, name_param, nargs);
        A = to< DistributedCSR< HpddmType >* >(args[0]);
        if (c == 1 || c == 3 || c == 4)
          K = to< long >(args[1]);
        else
          K = to< Matrice_Creuse< upscaled_type<PetscScalar> >* >(args[1]);
        if (c == 0 || c == 1) {
          R = to< KN< KN< long > >* >(args[2]);
          D = to< KN< typename std::conditional<
            std::is_same< HpddmType, HpSchwarz< PetscScalar > >::value, double, long >::type >* >(
            args[3]);
        }
        else if (c == 4)
          R = to< long >(args[2]);
      }

      AnyType operator( )(Stack stack) const;
      operator aType( ) const { return atype< DistributedCSR< HpddmType >* >( ); }
    };
    E_F0* code(const basicAC_F0& args) const { return new E_initCSR(args, c); }
    initCSR( )
      : OneOperator(
          atype< DistributedCSR< HpddmType >* >( ), atype< DistributedCSR< HpddmType >* >( ),
          atype< Matrice_Creuse< upscaled_type<PetscScalar> >* >( ), atype< KN< KN< long > >* >( ),
          atype< KN<
            typename std::conditional< std::is_same< HpddmType, HpSchwarz< PetscScalar > >::value,
                                       double, long >::type >* >( )),
        c(0) {}
    initCSR(int)
      : OneOperator(
          atype< DistributedCSR< HpddmType >* >( ), atype< DistributedCSR< HpddmType >* >( ),
          atype< long >( ), atype< KN< KN< long > >* >( ),
          atype< KN<
            typename std::conditional< std::is_same< HpddmType, HpSchwarz< PetscScalar > >::value,
                                       double, long >::type >* >( )),
        c(1) {}
    initCSR(int, int)
      : OneOperator(atype< DistributedCSR< HpddmType >* >( ),
                    atype< DistributedCSR< HpddmType >* >( ),
                    atype< Matrice_Creuse< upscaled_type<PetscScalar> >* >( )),
        c(2) {}
    initCSR(int, int, int)
      : OneOperator(atype< DistributedCSR< HpddmType >* >( ),
                    atype< DistributedCSR< HpddmType >* >( ), atype< long >( )),
        c(3) {}
    initCSR(int, int, int, int)
      : OneOperator(atype< DistributedCSR< HpddmType >* >( ),
                    atype< DistributedCSR< HpddmType >* >( ), atype< long >( ), atype< long >( )),
        c(4) {}
  };
  template< class HpddmType, bool C >
  basicAC_F0::name_and_type initCSR< HpddmType, C >::E_initCSR::name_param[] = {
    {"communicator", &typeid(pcommworld)},
    {"bs", &typeid(long)},
    {"clean", &typeid(bool)},
    {"symmetric", &typeid(bool)},
    {"restriction", &typeid(Matrice_Creuse< double >*)},
    {"level", &typeid(long)}};
  template< class HpddmType, bool C >
  AnyType initCSR< HpddmType, C >::E_initCSR::operator( )(Stack stack) const {
    DistributedCSR< HpddmType >* ptA = GetAny< DistributedCSR< HpddmType >* >((*A)(stack));
    KN< KN< long > >* ptR = (c == 0 || c == 1 ? GetAny< KN< KN< long > >* >((*R)(stack)) : nullptr);
    KN< typename std::conditional< std::is_same< HpddmType, HpSchwarz< PetscScalar > >::value,
                                   double, long >::type >* ptD =
      (c == 0 || c == 1
         ? GetAny< KN< typename std::conditional<
             std::is_same< HpddmType, HpSchwarz< PetscScalar > >::value, double, long >::type >* >(
             (*D)(stack))
         : nullptr), *empty = nullptr;
    PetscInt bs = nargs[1] ? GetAny< long >((*nargs[1])(stack)) : 1;
    int dof = 0;
    MatriceMorse< upscaled_type<PetscScalar> >* mA = nullptr;
    if (c == 0 || c == 2)
      mA = static_cast< MatriceMorse< upscaled_type<PetscScalar> >* >(
        &(*GetAny< Matrice_Creuse< upscaled_type<PetscScalar> >* >((*K)(stack))->A));
    else
      dof = GetAny< long >((*K)(stack));
    MPI_Comm* comm = nargs[0] ? (MPI_Comm*)GetAny< pcommworld >((*nargs[0])(stack)) : 0;
    if ((c == 2 && mA->n != mA->m) || c == 4) {
      ptA->_first = 0;
      ptA->_cfirst = 0;
      if (c != 4) {
        ptA->_last = mA->n;
        ptA->_clast = mA->m;
      }
      else {
        ptA->_last = dof;
        ptA->_clast = GetAny< long >((*R)(stack));
      }
      MatCreate(comm ? *comm : PETSC_COMM_WORLD, &ptA->_petsc);
      MatSetSizes(ptA->_petsc, ptA->_last, ptA->_clast, PETSC_DECIDE, PETSC_DECIDE);
      MatSetType(ptA->_petsc, MATAIJ);
      MatSetUp(ptA->_petsc);
      ptA->_num = new PetscInt[ptA->_last + ptA->_clast];
      ptA->_cnum = ptA->_num + ptA->_last;
      PetscInt rbegin, cbegin;
      MatGetOwnershipRange(ptA->_petsc, &rbegin, NULL);
      MatGetOwnershipRangeColumn(ptA->_petsc, &cbegin, NULL);
      std::iota(ptA->_num, ptA->_cnum, rbegin);
      std::iota(ptA->_cnum, ptA->_cnum + ptA->_clast, cbegin);
      ptA->_first += rbegin;
      ptA->_last += rbegin;
      ptA->_cfirst += cbegin;
      ptA->_clast += cbegin;
      if (c != 4) {
        ff_HPDDM_MatrixCSR< PetscScalar > dA(mA);
        if(cbegin)
          for(int i = 0; i < dA.HPDDM_nnz; ++i)
            dA.HPDDM_ja[i] += cbegin;
#if defined(PETSC_USE_64BIT_INDICES)
        PetscInt* ia = new PetscInt[dA.HPDDM_n + 1];
        std::copy_n(dA.HPDDM_ia, dA.HPDDM_n + 1, ia);
        PetscInt* ja = new PetscInt[dA.HPDDM_nnz];
        std::copy_n(dA.HPDDM_ja, dA.HPDDM_nnz, ja);
        PetscScalar* c = new PetscScalar[dA.HPDDM_nnz];
        std::copy_n(dA.HPDDM_a, dA.HPDDM_nnz, c);
        MatSeqAIJSetPreallocationCSR(ptA->_petsc, ia, ja, c);
        MatMPIAIJSetPreallocationCSR(ptA->_petsc, ia, ja, c);
        delete[] c;
        delete[] ja;
        delete[] ia;
#else
        MatSeqAIJSetPreallocationCSR(ptA->_petsc, dA.HPDDM_ia, dA.HPDDM_ja, dA.HPDDM_a);
        MatMPIAIJSetPreallocationCSR(ptA->_petsc, dA.HPDDM_ia, dA.HPDDM_ja, dA.HPDDM_a);
#endif
        if(cbegin)
          for(int i = 0; i < dA.HPDDM_nnz; ++i)
            dA.HPDDM_ja[i] -= cbegin;
      }
    }
    else {
      ptA->_A = new HpddmType;
      if ((ptR || c == 2 || c == 3) && (mA || (c != 0 && c != 2))) {
        HPDDM::MatrixCSR< PetscScalar >* dA;
        if (mA)
          dA = new_HPDDM_MatrixCSR< PetscScalar >(mA);
        else
          dA = new HPDDM::MatrixCSR< PetscScalar >(dof, dof, 0, nullptr, nullptr, nullptr, false);
        Matrice_Creuse< double >* pList =
          nargs[4] ? GetAny< Matrice_Creuse< double >* >((*nargs[4])(stack)) : 0;
        HPDDM::MatrixCSR< void >* dL = nullptr;
        int level = nargs[5] ? std::abs(GetAny< long >((*nargs[5])(stack))) : 0;
        KN_< KN< long > > sub((c == 0 || c == 1) && ptR->n > 0 && ptR->operator[](0).n > 0
                                ? (*ptR)(FromTo(1 + level * ptR->operator[](0).n,
                                                1 + (level + 1) * ptR->operator[](0).n - 1))
                                : KN< KN< long > >( ));
        if (std::is_same< HpddmType, HpSchwarz< PetscScalar > >::value && pList &&
            (mA || (c != 0 && c != 2))) {
          int n = 0;
          ptA->_exchange = new HPDDM::template Subdomain< PetscScalar >*[2]( );
          ptA->_exchange[0] = new HPDDM::template Subdomain< PetscScalar >( );
          ptA->_exchange[0]->initialize(
            dA, STL< long >((c == 0 || c == 1) && ptR->n > 0 ? ptR->operator[](0) : KN< long >( )),
            sub, comm);
          ptA->_exchange[0]->setBuffer( );
          if (pList->A) {
            MatriceMorse< double >* mList = static_cast< MatriceMorse< double >* >(&*(pList->A));
            mList->CSR( );
            ffassert(mList->n == mList->nnz);
            ffassert(mList->m == (mA ? mA->n : dof));
            n = mList->n;
            dL = new HPDDM::MatrixCSR< void >(n, mA ? mA->n : dof, n, mList->p, mList->j, false);
            ptA->_D = new KN<PetscReal>(n);
            for (int i = 0; i < n; ++i) ptA->_D->operator[](i) = ptD->operator[](mList->j[i]);
          } else {
            dL = new HPDDM::MatrixCSR< void >(0, mA ? mA->n : dof, 0, nullptr, nullptr, false);
            empty = new KN< typename std::conditional< std::is_same< HpddmType, HpSchwarz< PetscScalar > >::value,
                                     double, long >::type >(0);
          }
        }
        ptA->_A->HPDDM::template Subdomain< PetscScalar >::initialize(
          dA, STL< long >((c == 0 || c == 1) && ptR->n > 0 ? ptR->operator[](0) : KN< long >( )), sub,
          comm, dL);
        delete dL;
        ptA->_num = new PetscInt[ptA->_A->getMatrix( )->HPDDM_n];
        if (!C && (c == 2 || c == 3)) {
            long long global;
            ptA->_A->distributedNumbering(ptA->_num, ptA->_first, ptA->_last, global);
        }
        initPETScStructure<C>(ptA, bs,
          nargs[3] && GetAny< bool >((*nargs[3])(stack)) ? PETSC_TRUE : PETSC_FALSE, empty ? empty : ptD);
        delete empty;
        if (!std::is_same< HpddmType, HpSchwarz< PetscScalar > >::value) {
          mA->p = ptA->_A->getMatrix( )->HPDDM_ia;
        }
        if (nargs[2] && GetAny< bool >((*nargs[2])(stack))) {
          if (ptR) ptR->resize(0);
          if (c == 0 || c == 2) GetAny< Matrice_Creuse< PetscScalar >* >((*K)(stack))->destroy( );
        }
      }
    }
    return ptA;
  }
  template< class HpddmType >
  class initCSRfromBlockMatrix : public E_F0 {
   public:
    typedef DistributedCSR< HpddmType >* Result;
    PetscInt N, M;
    int init;
    Expression emat;
    Expression** e_M;
    int** t_M;
    initCSRfromBlockMatrix(const basicAC_F0& args, int i = 0) {
      init = i;
      M = 0;
      args.SetNameParam( );
      emat = args[0];
      const E_Array& eM = *dynamic_cast< const E_Array* >((Expression)args[1]);
      N = eM.size( );
      int err = 0;
      for (int i = 0; i < N; ++i) {
        const E_Array* e = dynamic_cast< const E_Array* >((Expression)eM[i]);
        if (!e)
          ++err;
        else {
          if (i == 0)
            M = e->size( );
          else if (M != e->size( ))
            ++err;
        }
      }
      if (err) CompileError("Wrong format of block matrix");
      assert(N && M);
      e_M = new Expression*[N];
      t_M = new int*[N];
      for (int i = 0; i < N; ++i) {
        const E_Array l = *dynamic_cast< const E_Array* >((Expression)eM[i]);
        e_M[i] = new Expression[M];
        t_M[i] = new int[M];
        for (int j = 0; j < M; ++j) {
          C_F0 c_M(l[j]);
          Expression e = c_M.LeftValue( );
          aType r = c_M.left( );
          if (r == atype< long >( ) && e->EvaluableWithOutStack( )) {
            long c = GetAny< long >((*e)(NullStack));
            if (c == 0) {
              e_M[i][j] = 0;
              t_M[i][j] = 0;
            } else if (atype< PetscScalar >( )->CastingFrom(r)) {
              e_M[i][j] = to< PetscScalar >(c_M);
              t_M[i][j] = 7;
            } else
              CompileError("Empty block matrix");
          } else if (r == atype< DistributedCSR< HpddmType >* >( )) {
            e_M[i][j] = e;
            t_M[i][j] = 1;
          } else if (r == atype< OpTrans< DistributedCSR< HpddmType > > >( )) {
            e_M[i][j] = e;
            t_M[i][j] = 2;
          } else if (atype< KN_< PetscScalar > >( )->CastingFrom(r)) {
            e_M[i][j] = to< KN_< PetscScalar > >(c_M);
            t_M[i][j] = 3;
          } else if (atype< Transpose< KN_< PetscScalar > > >( )->CastingFrom(r)) {
            e_M[i][j] = to< Transpose< KN_< PetscScalar > > >(c_M);
            t_M[i][j] = 4;
          } else if (atype< PetscScalar >( )->CastingFrom(r)) {
            e_M[i][j] = to< PetscScalar >(c_M);
            t_M[i][j] = 7;
          } else if (atype< std::pair<Dmat*, Dmat*>* >( )->CastingFrom(r)) {
            e_M[i][j] = to< std::pair<Dmat*, Dmat*>* >(c_M);
            t_M[i][j] = 8;
          } else
            CompileError("Unsupported type in submatrix");
        }
      }
    }
    ~initCSRfromBlockMatrix( ) {
      if (e_M) {
        for (int i = 0; i < N; ++i) {
          delete[] e_M[i];
          delete[] t_M[i];
        }
        delete[] e_M;
        delete[] t_M;
      }
    }
    static ArrayOfaType typeargs( ) {
      return ArrayOfaType(atype< Result >( ), atype< E_Array >( ));
    }
    static E_F0* f(const basicAC_F0& args) { return new initCSRfromBlockMatrix(args, 0); }
    AnyType operator( )(Stack s) const {
      Dmat** exchange = new Dmat*[N * M]();
      Mat* a = new Mat[N * M]( );
      std::vector< int > zeros;
      zeros.reserve(std::min(N, M));
      std::vector<std::pair<int, int>> destroy;
      destroy.reserve(N * M);
      MPI_Comm comm1 = MPI_COMM_NULL, comm2 = MPI_COMM_NULL;
      for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M; ++j) {
          Expression e = e_M[i][j];
          int t = t_M[i][j];
          if (e) {
            AnyType e_ij = (*e)(s);
            if (t == 1) {
              DistributedCSR< HpddmType >* pt = GetAny< DistributedCSR< HpddmType >* >(e_ij);
              a[i * M + j] = pt->_petsc;
              PetscObjectGetComm((PetscObject)pt->_petsc, &comm2);
              if (comm1 != MPI_COMM_NULL) {
                int flag;
                MPI_Comm_compare(comm1, comm2, &flag);
                ffassert(flag == MPI_CONGRUENT || flag == MPI_IDENT);
              }
              comm1 = comm2;
              exchange[i * M + j] = pt;
            } else if (t == 2) {
              DistributedCSR< HpddmType >* pt = GetAny< DistributedCSR< HpddmType >* >(e_ij);
              Mat B;
              MatCreateHermitianTranspose(pt->_petsc, &B);
              a[i * M + j] = B;
              destroy.emplace_back(i, j);
              exchange[i * M + j] = pt;
              PetscObjectGetComm((PetscObject)pt->_petsc, &comm2);
              if (comm1 != MPI_COMM_NULL) {
                int flag;
                MPI_Comm_compare(comm1, comm2, &flag);
                ffassert(flag == MPI_CONGRUENT || flag == MPI_IDENT);
              }
              comm1 = comm2;
            } else if (t == 7) {
              PetscScalar r = GetAny< PetscScalar >(e_ij);
              if (std::abs(r) > 1.0e-16) {
                Mat B;
                MatCreateDense(comm1 != MPI_COMM_NULL ? comm1 : PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, 1, 1, NULL, &B);
                PetscScalar* array;
                MatDenseGetArray(B, &array);
                if (array) array[0] = r;
                MatDenseRestoreArray(B, &array);
                a[i * M + j] = B;
                destroy.emplace_back(i, j);
              } else if (i == j)
                zeros.emplace_back(i);
            } else if (t == 3 || t == 4) {
              KN< PetscScalar > x;
              Mat B;
              if (t == 3) x = GetAny< KN_< PetscScalar > >(e_ij);
              else x = GetAny< Transpose< KN_< PetscScalar > > >(e_ij);
              MatCreateDense(comm1 != MPI_COMM_NULL ? comm1 : PETSC_COMM_WORLD, x.n, PETSC_DECIDE, PETSC_DECIDE, 1, NULL, &B);
              PetscScalar* array;
              MatDenseGetArray(B, &array);
              if (array) std::copy_n(static_cast< PetscScalar* >(x), x.n, array);
              MatDenseRestoreArray(B, &array);
              if (t == 3) {
                a[i * M + j] = B;
              }
              else {
                Mat C;
                MatCreateHermitianTranspose(B, &C);
                MatDestroy(&B);
                a[i * M + j] = C;
              }
              destroy.emplace_back(i, j);
            } else if (t == 8) {
              std::pair<Dmat*, Dmat*>* p = GetAny<std::pair<Dmat*, Dmat*>*>(e_ij);
              ffassert(p && p->first->_petsc && p->second->_petsc);
              Mat mats[2] = { p->second->_petsc, p->first->_petsc };
              Mat C;
              MatCreateComposite(PetscObjectComm((PetscObject)p->first->_petsc), 2, mats, &C);
              MatCompositeSetType(C, MAT_COMPOSITE_MULTIPLICATIVE);
              a[i * M + j] = C;
              delete p;
              destroy.emplace_back(i, j);
            } else {
              ExecError("Unknown type in submatrix");
            }
          } else if (i == j)
            zeros.emplace_back(i);
        }
      }
      for (int i = 0; i < zeros.size( ); ++i) {
        int posX = -1, posY = -1;
        for (int j = 0; j < M && posX == -1; ++j) {
          if (j != zeros[i] && a[zeros[i] * M + j]) posX = j;
        }
        for (int j = 0; j < N && posY == -1; ++j) {
          if (j != zeros[i] && a[zeros[i] + j * M]) posY = j;
        }
        if (posX == -1 && posY == -1)
          ExecError("Zero row and zero column");
        else {
          PetscInt x, X, y, Y;
          if (posX != -1) {
            MatGetSize(a[zeros[i] * M + posX], &X, &Y);
            MatGetLocalSize(a[zeros[i] * M + posX], &x, &y);
          }
          if (posY != -1) {
            MatGetSize(a[zeros[i] + posY * M], posX == -1 ? &X : NULL, &Y);
            MatGetLocalSize(a[zeros[i] + posY * M], posX == -1 ? &x : NULL, &y);
            if (posX == -1) {
              x = y;
              X = Y;
            }
          } else {
            y = x;
            Y = X;
          }
          MatCreate(comm1 != MPI_COMM_NULL ? comm1 : PETSC_COMM_WORLD, a + zeros[i] * M + zeros[i]);
          MatSetSizes(a[zeros[i] * M + zeros[i]], x, y, X, Y);
          MatSetType(a[zeros[i] * M + zeros[i]], MATAIJ);
          MatSeqAIJSetPreallocation(a[zeros[i] * M + zeros[i]], 0, NULL);
          MatMPIAIJSetPreallocation(a[zeros[i] * M + zeros[i]], 0, NULL, 0, NULL);
          MatAssemblyBegin(a[zeros[i] * M + zeros[i]], MAT_FINAL_ASSEMBLY);
          MatAssemblyEnd(a[zeros[i] * M + zeros[i]], MAT_FINAL_ASSEMBLY);
          destroy.emplace_back(zeros[i], zeros[i]);
        }
      }
      Result sparse_mat = GetAny< Result >((*emat)(s));
      if (sparse_mat->_petsc) sparse_mat->dtor( );
      MatCreateNest(comm1 != MPI_COMM_NULL ? comm1 : PETSC_COMM_WORLD, N, NULL, M, NULL, a, &sparse_mat->_petsc);
      for(std::pair<int, int> p : destroy)
          MatDestroy(a + p.first * M + p.second);
      sparse_mat->_exchange = reinterpret_cast<HPDDM::Subdomain<PetscScalar>**>(exchange);
      delete[] a;
      return sparse_mat;
    }
  };
  template< class HpddmType >
  class assignBlockMatrix : public initCSRfromBlockMatrix< HpddmType > {
   public:
    assignBlockMatrix(const basicAC_F0& args) : initCSRfromBlockMatrix< HpddmType >(args, 1) {}
    static E_F0* f(const basicAC_F0& args) {
      return new initCSRfromBlockMatrix< HpddmType >(args, 1);
    }
  };

  template< class Type >
  class setOptions : public OneOperator {
   public:
    const int c;
    class setOptions_Op : public E_F0mps {
     public:
      Expression A;
      Expression P;
      const int c;
      static const int n_name_param = 21;
      static basicAC_F0::name_and_type name_param[];
      Expression nargs[n_name_param];
      setOptions_Op(const basicAC_F0& args, int d) : A(0), P(0), c(d) {
        args.SetNameParam(n_name_param, name_param, nargs);
        if (c == 0)
          A = to< Type* >(args[0]);
        else {
          A = to< KN< Dmat >* >(args[0]);
          if (c == 1)
            P = to< KN< Matrice_Creuse< double > >* >(args[1]);
          else if (c == 2)
            P = to< long >(args[1]);
          else
            P = to< KN< Dmat >* >(args[1]);
        }
      }

      AnyType operator( )(Stack stack) const;
      operator aType( ) const { return atype< long >( ); }
    };
    E_F0* code(const basicAC_F0& args) const { return new setOptions_Op(args, c); }
    setOptions( ) : OneOperator(atype< long >( ), atype< Type* >( )), c(0) {}
    setOptions(int)
      : OneOperator(atype< long >( ), atype< KN< Dmat >* >( ),
                    atype< KN< Matrice_Creuse< double > >* >( )),
        c(1) {}
    setOptions(int, int)
      : OneOperator(atype< long >( ), atype< KN< Dmat >* >( ), atype< long >( )), c(2) {}
    setOptions(int, int, int)
      : OneOperator(atype< long >( ), atype< KN< Dmat >* >( ), atype< KN< Dmat >* >( )), c(3) {}
  };
  template< class Type >
  basicAC_F0::name_and_type setOptions< Type >::setOptions_Op::name_param[] = {
    {"sparams", &typeid(std::string*)},                                                           // 0
    {"nearnullspace", &typeid(FEbaseArrayKn< upscaled_type<PetscScalar> >*)},                     // 1
    {"fields", &typeid(KN< double >*)},                                                           // 2
    {"names", &typeid(KN< String >*)},                                                            // 3
    {"prefix", &typeid(std::string*)},                                                            // 4
    {"schurPreconditioner", &typeid(KN< Matrice_Creuse< upscaled_type<PetscScalar> > >*)},        // 5
    {"schurList", &typeid(KN< double >*)},                                                        // 6
    {"parent", &typeid(Type*)},                                                                   // 7
    {"MatNullSpace", &typeid(KNM< upscaled_type<PetscScalar> >*)},                                // 8
    {"fieldsplit", &typeid(long)},                                                                // 9
    {"schurComplement", &typeid(KNM< upscaled_type<PetscScalar> >*)},                             // 10
    {"schur", &typeid(KN< Dmat >*)},                                                              // 11
    {"aux", &typeid(Matrice_Creuse< upscaled_type<PetscScalar> >*)},                              // 12
    {"coordinates", &typeid(KNM< double >*)},                                                     // 13
    {"gradient", &typeid(Dmat*)},                                                                 // 14
    {"O", &typeid(Matrice_Creuse< upscaled_type<PetscScalar> >*)},                                // 15
    {"bs", &typeid(long)},                                                                        // 16
    {"precon", &typeid(Polymorphic*)},                                                            // 17
    {"setup", &typeid(bool)},                                                                     // 18
    {"monitor", &typeid(Polymorphic*)},                                                           // 19
    {"deflation", &typeid(FEbaseArrayKn<upscaled_type<PetscScalar>>*)}                            // 20
  };
  class ShellInjection;
  template< class Type >
  class LinearSolver;
  template< class Type >
  class NonlinearSolver;
  template< class Type >
  class TimeStepper;
  template< class Type >
  class TaoSolver;
  template< class Type >
  struct _n_User {
    typename Type::MatF_O* op;
  };
  template<>
  struct _n_User< ShellInjection > {
    Matrice_Creuse< double >* P;
    Dmat* f;
    Dmat* C;
  };
  template< class Type >
  struct _n_User< NonlinearSolver< Type > > {
    typename NonlinearSolver< Type >::VecF_O* op;
    typename LinearSolver< Type >::MatF_O* r;
    typename NonlinearSolver< Type >::IConvF_O* conv;
  };
  template< class Type >
  struct _n_User< TimeStepper< Type > > {
    typename NonlinearSolver< Type >::IVecF_O* op;
    typename NonlinearSolver< Type >::IMatF_O* r;
    typename NonlinearSolver< Type >::IMatF_O* rhs;
    typename NonlinearSolver< Type >::IMonF_O* mon;
  };
  template< class Type >
  struct _n_User< TaoSolver< Type > > {
    typename NonlinearSolver< Type >::VecF_O* J;
    typename LinearSolver< Type >::MatF_O* r;
    typename NonlinearSolver< Type >::VecF_O* op;
    typename LinearSolver< Type >::MatF_O* ic;
    typename NonlinearSolver< Type >::VecF_O* icJ;
    typename LinearSolver< Type >::MatF_O* ec;
    typename NonlinearSolver< Type >::VecF_O* ecJ;
  };
  static Mat* O;
  PetscErrorCode CustomCreateSubMatrices(Mat,PetscInt,const IS*,const IS*,MatReuse scall,Mat *submat[]) {
    PetscErrorCode ierr;

    PetscFunctionBeginUser;
    if (scall == MAT_INITIAL_MATRIX) {
      ierr = PetscCalloc1(1,submat);CHKERRQ(ierr);
    }
    (*submat)[0] = *O;
    PetscFunctionReturn(PETSC_SUCCESS);
  }
  template< class T, typename std::enable_if< std::is_same< T, KN< PetscScalar > >::value || std::is_same< T, KN< long > >::value >::type* =
                       nullptr >
  void resize(T* v, int n, int m) {
    v->resize(n);
  }
  template< class T, typename std::enable_if< std::is_same< T, KNM< PetscScalar > >::value || std::is_same< T, KNM< long > >::value >::type* =
                       nullptr >
  void resize(T* v, int n, int m) {
    v->resize(n, m);
  }
  template< class T, typename std::enable_if<
                       !std::is_same< T, KN< PetscScalar > >::value && !std::is_same< T, KNM< PetscScalar > >::value &&
                       !std::is_same< T, KN< long > >::value && !std::is_same< T, KNM< long > >::value>::type* = nullptr >
  void resize(T* v, int n, int m) {}
  template< class T, class U >
  void changeNumbering_func(PetscInt* const num, PetscInt first, PetscInt last,
                            PetscInt m, PetscInt n, PetscInt bs, T* ptIn, U* ptOut, bool inverse) {
    PetscScalar* out;
    if (!inverse) {
      resize(ptOut, m ? m : n, ptIn->M());
      out = ptOut->operator PetscScalar*();
      if (last - first) {
        for(int i = 0; i < ptIn->M(); ++i) {
          HPDDM::Subdomain< PetscScalar >::template distributedVec< 0 >(
            num, first, last, ptIn->operator PetscScalar*() + i * ptIn->N(), out, static_cast<PetscInt>(ptIn->N() / bs), bs);
          out += last - first;
        }
      }
      else
        std::copy_n(ptIn->operator PetscScalar*(), ptIn->N() * ptOut->M(), out);
    } else {
      resize(ptIn, n, ptOut->M());
      if(ptIn->N())
          *ptIn = PetscScalar( );
      out = ptOut->operator PetscScalar*();
      if (num) {
        for(int i = 0; i < ptIn->M(); ++i) {
          HPDDM::Subdomain< PetscScalar >::template distributedVec< 1 >(
            num, first, last, ptIn->operator PetscScalar*() + i * ptIn->N(), out, static_cast<PetscInt>(ptIn->N() / bs), bs);
          out += last - first;
        }
      }
      else
        std::copy_n(out, ptIn->N() * ptOut->M(), ptIn->operator PetscScalar*());
    }
  }
  template< bool T, class K, class U, typename std::enable_if< std::is_same< K, upscaled_type<PetscScalar> >::value && std::is_same< U, upscaled_type<U>>::value >::type* = nullptr >
  void MatMult(MatriceMorse<K>* const& A, KN_<U>& in, KN_<U>& out) {
    A->addMatMul(in, out, T, in.step, out.step);
  }
  template< bool T, class K, class U, typename std::enable_if< !std::is_same< K, upscaled_type<PetscScalar> >::value && std::is_same< U, upscaled_type<U>>::value >::type* = nullptr >
  void MatMult(MatriceMorse<K>* const& A, KN_<U>& in, KN_<U>& out) {
    upscaled_type<PetscScalar>* pc = in;
    double* pr = reinterpret_cast<double*>(pc);
    KN_< K > realIn(pr + 0, in.N(), 2 * in.step);
    KN_< K > imagIn(pr + 1, in.N(), 2 * in.step);
    pc = out;
    pr = reinterpret_cast<double*>(pc);
    KN_< K > realOut(pr + 0, out.N(), 2 * out.step);
    KN_< K > imagOut(pr + 1, out.N(), 2 * out.step);
    A->addMatMul(realIn, realOut, T, 2 * in.step, 2 * out.step);
    A->addMatMul(imagIn, imagOut, T, 2 * in.step, 2 * out.step);
  }
  template< bool T, class K, class U, typename std::enable_if< !std::is_same< U, upscaled_type<U>>::value >::type* = nullptr >
  void MatMult(MatriceMorse<K>* const& A, KN_<U>& in, KN_<U>& out) {
    KN< upscaled_type<PetscScalar> > inUp(in.N());
    KN< upscaled_type<PetscScalar> > outUp(out.N());
    for(int j = 0; j < in.N(); ++j)
        inUp[j] = in[j];
    MatMult<T>(A, inUp, outUp);
    for(int j = 0; j < out.N(); ++j)
        out[j] = outUp[j];
  }
  template< bool T >
  static PetscErrorCode ShellInjectionOp(Mat A, Vec x, Vec y) {
    User< ShellInjection > user;
    const PetscScalar* in;
    PetscScalar* out;
    PetscErrorCode ierr;

    PetscFunctionBeginUser;
    ierr = MatShellGetContext(A, &user);CHKERRQ(ierr);
    MatriceMorse< double >* mP = user->P->A ? static_cast< MatriceMorse< double >* >(&(*user->P->A)) : nullptr;
    if (mP) {
      ierr = VecGetArrayRead(x, &in);CHKERRQ(ierr);
      ierr = VecGetArray(y, &out);CHKERRQ(ierr);
      if (!T) {
        PetscInt mFine, nCoarse;
        MatGetLocalSize(A, &mFine, nullptr);
        MatGetLocalSize(A, nullptr, &nCoarse);
        KN< PetscScalar > coarse(user->C->_A->getDof( ));
        KN_< PetscScalar > coarseIn(const_cast< PetscScalar* >(in), nCoarse);
        changeNumbering_func(user->C->_num, user->C->_first, user->C->_last, nCoarse,
                             user->C->_A->getDof(), 1, &coarse, &coarseIn, true);
        KN< PetscScalar > fine(user->f->_A->getDof( ));
        fine = PetscScalar( );
        user->C->_A->exchange(coarse);
        MatMult<false>(mP, coarse, fine);
        KN_< PetscScalar > fineOut(out, mFine);
        fineOut = PetscScalar( );
        changeNumbering_func(user->f->_num, user->f->_first, user->f->_last, mFine,
                             user->f->_A->getDof(), 1, &fine, &fineOut, false);
      } else {
        PetscInt nFine, mCoarse;
        MatGetLocalSize(A, nullptr, &nFine);
        MatGetLocalSize(A, &mCoarse, nullptr);
        KN< PetscScalar > fine(user->f->_A->getDof( ));
        KN_< PetscScalar > fineIn(const_cast< PetscScalar* >(in), nFine);
        changeNumbering_func(user->f->_num, user->f->_first, user->f->_last, nFine,
                             user->f->_A->getDof(), 1, &fine, &fineIn, true);
        KN< PetscScalar > coarse(user->C->_A->getDof( ));
        coarse = PetscScalar( );
        user->f->_A->exchange(fine);
        MatMult<true>(mP, fine, coarse);
        KN_< PetscScalar > coarseOut(out, mCoarse);
        changeNumbering_func(user->C->_num, user->C->_first, user->C->_last, mCoarse,
                             user->C->_A->getDof(), 1, &coarse, &coarseOut, false);
      }
      ierr = VecRestoreArray(y, &out);CHKERRQ(ierr);
      ierr = VecRestoreArrayRead(x, &in);CHKERRQ(ierr);
    } else VecCopy(x, y);
    PetscFunctionReturn(PETSC_SUCCESS);
  }
  template< class Type, typename std::enable_if<!std::is_same<Type, ShellInjection>::value>::type* = nullptr  >
  static PetscErrorCode ShellDestroy(Mat A) {
    User< Type > user;
    PetscErrorCode ierr;

    PetscFunctionBeginUser;
    ierr = MatShellGetContext(A, &user);CHKERRQ(ierr);
    delete user->op;
    ierr = PetscFree(user);CHKERRQ(ierr);
    PetscFunctionReturn(PETSC_SUCCESS);
  }
  template< class Type, typename std::enable_if<std::is_same<Type, ShellInjection>::value>::type* = nullptr  >
  static PetscErrorCode ShellDestroy(Mat A) {
    User< Type > user;
    PetscErrorCode ierr;

    PetscFunctionBeginUser;
    ierr = MatShellGetContext(A, &user);CHKERRQ(ierr);
    ierr = PetscFree(user);CHKERRQ(ierr);
    PetscFunctionReturn(PETSC_SUCCESS);
  }
  template< class Type >
  static PetscErrorCode PCShellDestroy(PC pc) {
    User< Type > user;
    PetscErrorCode ierr;

    PetscFunctionBeginUser;
    ierr = PCShellGetContext(pc, (void**)&user);CHKERRQ(ierr);
    delete user->op;
    ierr = PetscFree(user);CHKERRQ(ierr);
    PetscFunctionReturn(PETSC_SUCCESS);
  }
  template< class Type >
  static PetscErrorCode MonitorDestroy(void** ctx) {
    PetscFunctionBeginUser;
    delete reinterpret_cast< typename LinearSolver< Type >::MonF_O* >(*ctx);
    PetscFunctionReturn(PETSC_SUCCESS);
  }
  struct WrapperSubKSP {
    IS  is;
    KSP *ksp;
  };
  static PetscErrorCode WrapperApply(PC pc, Vec in, Vec out) {
    WrapperSubKSP *ctx;
    PetscFunctionBeginUser;
    VecCopy(in, out);
    PCShellGetContext(pc, (void**)&ctx);
    VecPermute(out, ctx->is, PETSC_FALSE);
    KSPSolve(*ctx->ksp, out, out);
    VecPermute(out, ctx->is, PETSC_TRUE);
    PetscFunctionReturn(PETSC_SUCCESS);
  }
  static PetscErrorCode WrapperDestroy(PC pc) {
    WrapperSubKSP *ctx;
    PetscFunctionBeginUser;
    PCShellGetContext(pc, (void**)&ctx);
    ISDestroy(&ctx->is);
    PetscFree(ctx);
    PetscFunctionReturn(PETSC_SUCCESS);
  }
  static PetscErrorCode WrapperView(PC pc, PetscViewer viewer) {
    WrapperSubKSP *ctx;
    PetscFunctionBeginUser;
    PCShellGetContext(pc, (void**)&ctx);
    KSPView(*ctx->ksp, viewer);
    PetscFunctionReturn(PETSC_SUCCESS);
  }
  template< class, class Container, char N = 'N' >
  static PetscErrorCode Op_User(Container, Vec, Vec);
  template< class Type >
  PetscErrorCode Monitor(KSP ksp, PetscInt it, PetscReal rnorm, void* ctx) {
    PetscFunctionBeginUser;
    typename LinearSolver< Type >::MonF_O* mat =
      reinterpret_cast< typename LinearSolver< Type >::MonF_O* >(ctx);
    mat->apply(it, rnorm);
    PetscFunctionReturn(PETSC_SUCCESS);
  }
  template< class Type >
  AnyType setOptions< Type >::setOptions_Op::operator( )(Stack stack) const {
    Type *ptA, *ptParent;
    KN< Dmat >* tabA;
    int level = 0;
    if (c != 0) {
      tabA = GetAny< KN< Dmat >* >((*A)(stack));
      assert(tabA->N( ) > 0);
      if (c == 2) {
        level = GetAny< long >((*P)(stack));
        assert(0 <= level && level < tabA->N( ));
        level = tabA->N( ) - level - 1;
      }
      ptA = reinterpret_cast< Type* >(&tabA->operator[](level));
      level = tabA->N( ) - level - 1;
    }
    else
      ptA = GetAny< Type* >((*A)(stack));
    ptParent = nargs[7] ? GetAny< Type* >((*nargs[7])(stack)) : 0;
    if (ptParent && ptParent->_ksp) {
      KSP* subksp;
      PetscInt n;
      PC pc;
      KSPGetPC(ptParent->_ksp, &pc);
      PCType type;
      PCGetType(pc, &type);
      PetscBool isType;
      PetscStrcmp(type, PCFIELDSPLIT, &isType);
      if (!isType) {
        PetscStrcmp(type, PCASM, &isType);
        if (!isType)
          return 0L;
        else {
          PCASMGetSubKSP(pc, &n, NULL, &subksp);
          ffassert(n == 1);
          KSPSetOptionsPrefix(subksp[0], NULL);
          KSPSetType(subksp[0], KSPPREONLY);
          PC subpc;
          KSPGetPC(subksp[0], &subpc);
          PCSetType(subpc, PCSHELL);
          WrapperSubKSP* wrapper = nullptr;
          PetscNew(&wrapper);
          PCShellSetContext(subpc, wrapper);
          PCShellSetApply(subpc, WrapperApply);
          PCShellSetDestroy(subpc, WrapperDestroy);
          PCShellSetView(subpc, WrapperView);
          KSPSetSkipPCSetFromOptions(subksp[0], PETSC_TRUE);
          PetscInt* idx;
          PetscMalloc1(ptA->_A->getDof(), &idx);
          std::map<int, int> map;
          ffassert(ptA->_A->getDof() == ptParent->_A->getDof());
          for(int i = 0; i < ptA->_A->getDof(); ++i)
              map[ptParent->_num[i]] = i;
          int i = 0;
          for(const std::pair<const int, int>& v : map)
              idx[i++] = v.second;
          IS perm;
          ISCreateGeneral(PETSC_COMM_SELF, ptA->_A->getDof(), idx, PETSC_OWN_POINTER, &perm);
          ISInvertPermutation(perm, PETSC_DECIDE, &wrapper->is);
          ISDestroy(&perm);
          wrapper->ksp = &ptA->_ksp;
          PetscBool assembled = PETSC_FALSE;
          if (ptA->_petsc)
            MatAssembled(ptA->_petsc, &assembled);
          if (!assembled) {
            Mat A;
            MatType type;
            KSPGetOperators(subksp[0], &A, NULL);
            MatGetType(A, &type);
            PetscBool convert;
            PetscStrcmp(type, MATSEQSBAIJ, &convert);
            if(convert)
                MatConvert(A, MATSEQAIJ, MAT_INPLACE_MATRIX, &A);
            MatDestroy(&ptA->_petsc);
            MatPermute(A, wrapper->is, wrapper->is, &ptA->_petsc);
            if(convert)
                MatConvert(ptA->_petsc, MATSBAIJ, MAT_INPLACE_MATRIX, &ptA->_petsc);
          }
        }
      }
      else {
        PCFieldSplitGetSubKSP(pc, &n, &subksp);
        Mat** mat;
        PetscInt M, N;
        MatNestGetSubMats(ptParent->_petsc, &M, &N, &mat);
        ffassert(M == N);
        for (int i = 0; i < M; ++i) {
          if (mat[i][i] == ptA->_petsc) {
            KSPDestroy(&ptA->_ksp);
            ptA->_ksp = subksp[i];
            PetscObjectReference((PetscObject)subksp[i]);
            KSPSetOperators(subksp[i], ptA->_petsc, ptA->_petsc);
            break;
          }
        }
        PetscFree(subksp);
      }
    }
    if (nargs[0]) {
      std::string* options = GetAny< std::string* >((*nargs[0])(stack));
      PetscOptionsInsertString(NULL, options->c_str());
    }
    long FS = nargs[9] ? GetAny< long >((*nargs[9])(stack)) : -1;
    KSP ksp = nullptr;
    if (ptA && ptA->_petsc) {
      {
        long bs = nargs[16] ? GetAny< long >((*nargs[16])(stack)) : -1;
        if(bs >= 1)
            MatSetBlockSize(ptA->_petsc, bs);
      }
      PetscBool assembled;
      MatAssembled(ptA->_petsc, &assembled);
      if (!ptA->_ksp && c != 2) {
        KSPCreate(PetscObjectComm((PetscObject)ptA->_petsc), &ptA->_ksp);
        KSPSetOperators(ptA->_ksp, ptA->_petsc, ptA->_petsc);
      }
      if (nargs[4] && c != 2)
        KSPSetOptionsPrefix(ptA->_ksp, GetAny< std::string* >((*nargs[4])(stack))->c_str( ));
      if (c == 1 || c == 3) {
        KN< Matrice_Creuse< double > >* mP = (c == 1 ? GetAny< KN< Matrice_Creuse< double > >* >((*P)(stack)) : nullptr);
        KN< Dmat >* mD = (c == 3 ? GetAny< KN< Dmat >* >((*P)(stack)) : nullptr);
        ffassert((c == 1 && mP->N( ) + 1 == tabA->N( )) || (c == 3 && (mD->N( ) + 1 == tabA->N( ) || tabA->N( ) == 1)));
        ksp = ptA->_ksp;
        PC pc;
        KSPGetPC(ksp, &pc);
        PCSetType(pc, PCMG);
        const PetscInt level = (c == 3 ? mD->N( ) + 1 : tabA->N( ));
        PCMGSetLevels(pc, level, NULL);
        if (c == 3 && tabA->N( ) == 1)
          PCMGSetGalerkin(pc, PC_MG_GALERKIN_BOTH);
        else
          PCMGSetGalerkin(pc, PC_MG_GALERKIN_NONE);
        PCSetFromOptions(pc);
        for (int i = 0; i < level; ++i) {
          KSP smoother;
          PCMGGetSmoother(pc, level - i - 1, &smoother);
          if (c != 3 || tabA->N( ) > 1)
            KSPSetOperators(smoother, tabA->operator[](i)._petsc,
                            tabA->operator[](i)._petsc);
          if (i < level - 1) {
            if (c == 1) {
                User< ShellInjection > user = nullptr;
                PetscNew(&user);
                user->P = &mP->operator[](i);
                user->f = &tabA->operator[](i);
                user->C = &tabA->operator[](i + 1);
                Mat P;
                PetscInt mFine, MFine, mCoarse, MCoarse;
                MatGetLocalSize(tabA->operator[](i)._petsc, &mFine, nullptr);
                MatGetSize(tabA->operator[](i)._petsc, &MFine, nullptr);
                MatGetLocalSize(tabA->operator[](i + 1)._petsc, &mCoarse, nullptr);
                MatGetSize(tabA->operator[](i + 1)._petsc, &MCoarse, nullptr);
                MatCreateShell(PetscObjectComm((PetscObject)tabA->operator[](i + 1)._petsc), mFine, mCoarse, MFine, MCoarse, user, &P);
                MatShellSetOperation(P, MATOP_MULT, (void (*)(void))ShellInjectionOp< false >);
                MatShellSetOperation(P, MATOP_MULT_TRANSPOSE, (void (*)(void))ShellInjectionOp< true >);
                MatShellSetOperation(P, MATOP_DESTROY, (void (*)(void))ShellDestroy< ShellInjection >);
                PCMGSetInterpolation(pc, level - i - 1, P);
                MatDestroy(&P);
            } else {
                PCMGSetInterpolation(pc, level - i - 1, mD->operator[](i)._petsc);
                MatDestroy(&(mD->operator[](i)._petsc));
            }
          }
        }
      } else if (c == 2) {
        PC pc;
        KSPGetPC(tabA->operator[](0)._ksp, &pc);
        PCMGGetSmoother(pc, level, &ksp);
      } else if (std::is_same< Type, Dmat >::value) {
        PC pc;
        if (FS < 0)
          ksp = ptA->_ksp;
        else {
          KSPGetPC(ptA->_ksp, &pc);
          PetscInt nsplits;
          KSP* subksp;
          PCFieldSplitGetSubKSP(pc, &nsplits, &subksp);
          if (FS < nsplits) ksp = subksp[FS];
          else ffassert(0);
          PetscFree(subksp);
        }
        KSPGetPC(ksp, &pc);
        PCSetFromOptions(pc);
        PCType type;
        PCGetType(pc, &type);
        PetscBool isFieldSplit;
        PetscStrcmp(type, PCFIELDSPLIT, &isFieldSplit);
        if (isFieldSplit) {
          KN< double >* fields = nargs[2] ? GetAny< KN< double >* >((*nargs[2])(stack)) : 0;
          KN< String >* names = nargs[3] ? GetAny< KN< String >* >((*nargs[3])(stack)) : 0;
          KN< Matrice_Creuse< upscaled_type<PetscScalar> > >* mS =
            nargs[5] ? GetAny< KN< Matrice_Creuse< upscaled_type<PetscScalar> > >* >((*nargs[5])(stack)) : 0;
          KN< double >* pL = nargs[6] ? GetAny< KN< double >* >((*nargs[6])(stack)) : 0;
          KN< Dmat >* mdS = nargs[11] ? GetAny< KN< Dmat >* >((*nargs[11])(stack)) : 0;
          if (mdS)
            setFieldSplitPC(ptA, ksp, fields, names, mdS);
          else
            setFieldSplitPC(ptA, ksp, fields, names, mS, pL);
        }
      }
      else ksp = ptA->_ksp;
      KSPSetFromOptions(ksp);
      if (c != 1) {
        if (std::is_same< Type, Dmat >::value) {
          FEbaseArrayKn< upscaled_type<PetscScalar> >* ptNS =
            nargs[1] ? GetAny< FEbaseArrayKn< upscaled_type<PetscScalar> >* >((*nargs[1])(stack)) : 0;
          KNM< upscaled_type<PetscScalar> >* ptPETScNS =
            nargs[8] ? GetAny< KNM< upscaled_type<PetscScalar> >* >((*nargs[8])(stack)) : 0;
          int dim = ptNS ? ptNS->N : 0;
          int dimPETSc = ptPETScNS ? ptPETScNS->M( ) : 0;
          if (dim || dimPETSc) {
            Vec v;
            MatCreateVecs(ptA->_petsc, &v, NULL);
            Vec* ns;
            VecDuplicateVecs(v, std::max(dim, dimPETSc), &ns);
            VecDestroy(&v);
            if (std::max(dim, dimPETSc) == dimPETSc) {
              PetscInt m;
              VecGetLocalSize(ns[0], &m);
              for (unsigned short i = 0; i < dimPETSc; ++i) {
                PetscScalar* x;
                VecGetArray(ns[i], &x);
                for (int j = 0; j < m; ++j) x[j] = (*ptPETScNS)(j, i);
                VecRestoreArray(ns[i], &x);
              }
            } else
              for (unsigned short i = 0; i < dim; ++i) {
                PetscScalar* x;
                VecGetArray(ns[i], &x);
                upscaled_type<PetscScalar>* get = *ptNS->get(i);
                PetscScalar* base = reinterpret_cast<PetscScalar*>(get);
                if(!std::is_same<upscaled_type<PetscReal>, PetscReal>::value) {
                    for(int j = 0; j < ptNS->get(i)->n; ++j)
                        base[j] = get[j];
                }
                HPDDM::Subdomain< PetscScalar >::template distributedVec< 0 >(
                  ptA->_num, ptA->_first, ptA->_last, base,
                  x, static_cast<PetscInt>(ptNS->get(i)->n));
                VecRestoreArray(ns[i], &x);
              }
            PetscScalar* dots = new PetscScalar[std::max(dim, dimPETSc)];
            for (unsigned short i = 0; i < std::max(dim, dimPETSc); ++i) {
              if (i > 0) {
                VecMDot(ns[i], i, ns, dots);
                for (int j = 0; j < i; ++j) dots[j] *= -1.0;
                VecMAXPY(ns[i], i, dots, ns);
              }
              VecNormalize(ns[i], NULL);
            }
            delete[] dots;
            MatNullSpace sp;
            MatNullSpaceCreate(PetscObjectComm((PetscObject)ptA->_petsc), PETSC_FALSE, std::max(dim, dimPETSc), ns, &sp);
            MatSetNearNullSpace(ptA->_petsc, sp);
            MatNullSpaceDestroy(&sp);
            VecDestroyVecs(std::max(dim, dimPETSc), &ns);
          }
          PC pc;
          KSPGetPC(ksp, &pc);
          PCType type;
          PCGetType(pc, &type);
          PetscBool isType;
          KNM< upscaled_type<PetscReal> >* coordinates =
            nargs[13] ? GetAny< KNM< upscaled_type<PetscReal> >* >((*nargs[13])(stack)) : nullptr;
          Dmat* G = nargs[14] ? GetAny< Dmat* >((*nargs[14])(stack)) : nullptr;
          if (coordinates && std::is_same<PetscReal, upscaled_type<PetscReal>>::value) PCSetCoordinates(pc, coordinates->N( ), coordinates->M( ), reinterpret_cast<PetscReal*>(coordinates->operator double*()));
#if defined(PETSC_HAVE_HYPRE)
          if (G) {
            PetscStrcmp(type, PCHYPRE, &isType);
            if (isType) PCHYPRESetDiscreteGradient(pc, G->_petsc);
          }
#endif
          if (assembled && ptA->_A && ptA->_num) {
            PetscInt* idx;
            PetscMalloc1(ptA->_A->getDof(), &idx);
            std::copy_n(ptA->_num, ptA->_A->getDof(), idx);
            IS is;
            ISCreateGeneral(PETSC_COMM_SELF, ptA->_A->getDof(), idx, PETSC_OWN_POINTER,
                            &is);
#if defined(PETSC_HAVE_HPDDM) && defined(PETSC_USE_SHARED_LIBRARIES)
            PetscStrcmp(type, PCHPDDM, &isType);
            const HPDDM::MatrixCSR< PetscScalar >* const A = ptA->_A->getMatrix( );
            if (isType) {
              if(nargs[20]) {
#if PETSC_VERSION_GE(3, 18, 0)
                  Mat Z;
                  FEbaseArrayKn<upscaled_type<PetscScalar>>* deflation = GetAny<FEbaseArrayKn<upscaled_type<PetscScalar>>*>((*nargs[20])(stack));
                  MatCreateSeqDense(PETSC_COMM_SELF, ptA->_A->getDof(), deflation->N, NULL, &Z);
                  PetscScalar* data;
                  MatDenseGetArrayWrite(Z, &data);
                  for(int i = 0; i < deflation->N; ++i)
                    std::copy_n(&(*deflation->get(i))[0], deflation->get(0)->n, data + i * deflation->get(0)->n);
                  MatDenseRestoreArrayWrite(Z, &data);
                  PCHPDDMSetDeflationMat(pc, is, Z);
                  MatDestroy(&Z);
#else
                  ffassert(0);
#endif
              }
              else if(A) {
                if (!A->HPDDM_ia) {
                  Matrice_Creuse< upscaled_type<PetscScalar> >* ptK =
                    nargs[15] ? GetAny< Matrice_Creuse< upscaled_type<PetscScalar> >* >((*nargs[15])(stack)) : nullptr;
                  if (ptK && ptK->A) {
                    MatriceMorse< upscaled_type<PetscScalar> >* mA =
                      static_cast< MatriceMorse< upscaled_type<PetscScalar> >* >(&(*ptK->A));
                    HPDDM::MatrixCSR< PetscScalar >* B = new_HPDDM_MatrixCSR< PetscScalar >(mA);
                    Mat aux = ff_to_PETSc(B);
                    Mat N;
                    PetscObjectQuery((PetscObject)pc, "_PCHPDDM_Neumann_Mat", (PetscObject*)&N);
                    if(!N) {
                        PCHPDDMSetAuxiliaryMat(pc, is, aux, NULL, NULL);
                        PCSetFromOptions(pc);
                        MatDestroy(&aux);
                    }
                    else {
                        PetscObjectCompose((PetscObject)pc, "_PCHPDDM_Neumann_Mat", (PetscObject)aux);
                        PCHPDDMSetAuxiliaryMat(pc, NULL, aux, NULL, NULL);
                        MatDestroy(&aux);
                    }
                    delete B;
                  }
                }
                else {
                  Mat aux = ff_to_PETSc(A);
                  Mat N;
                  PetscObjectQuery((PetscObject)pc, "_PCHPDDM_Neumann_Mat", (PetscObject*)&N);
                  if(!N) {
                      PCHPDDMSetAuxiliaryMat(pc, is, aux, NULL, NULL);
                      PCSetFromOptions(pc);
                      MatDestroy(&aux);
                  }
                  else {
                      PetscObjectCompose((PetscObject)pc, "_PCHPDDM_Neumann_Mat", (PetscObject)aux);
                      PCHPDDMSetAuxiliaryMat(pc, NULL, aux, NULL, NULL);
                      MatDestroy(&aux);
                  }
                  Matrice_Creuse< upscaled_type<PetscScalar> >* ptK =
                    nargs[12] ? GetAny< Matrice_Creuse< upscaled_type<PetscScalar> >* >((*nargs[12])(stack)) : nullptr;
                  if (ptK && ptK->A) {
                    MatriceMorse< upscaled_type<PetscScalar> >* mA =
                      static_cast< MatriceMorse< upscaled_type<PetscScalar> >* >(&(*ptK->A));
                    HPDDM::MatrixCSR< PetscScalar >* B = new_HPDDM_MatrixCSR< PetscScalar >(mA);
                    aux = ff_to_PETSc(B);
                    PCHPDDMSetRHSMat(pc, aux);
                    MatDestroy(&aux);
                    delete B;
                  }
                }
              }
            }
            else {
#endif
              Matrice_Creuse< upscaled_type<PetscScalar> >* ptO =
                nargs[15] ? GetAny< Matrice_Creuse< upscaled_type<PetscScalar> >* >((*nargs[15])(stack)) : nullptr;
              if (ptO && ptO->A) {
                MatriceMorse< upscaled_type<PetscScalar> >* mO =
                  static_cast< MatriceMorse< upscaled_type<PetscScalar> >* >(&(*ptO->A));
                ff_HPDDM_MatrixCSR< PetscScalar > dO(mO);
                PCSetType(pc, PCASM);
                IS loc;
                PetscInt n, first;
                MatGetOwnershipRange(ptA->_petsc, &first, &n);
                n -= first;
                ISCreateStride(PETSC_COMM_SELF, n, first, 1, &loc);
                PCASMSetLocalSubdomains(pc, 1, &is, &loc);
                int nnz = dO.HPDDM_nnz;
                MPI_Allreduce(MPI_IN_PLACE, &nnz, 1, MPI_INT, MPI_MAX, PetscObjectComm((PetscObject)ptA->_petsc));
                if (nnz) {
                  MatSetOperation(ptA->_petsc, MATOP_CREATE_SUBMATRICES, (void(*)(void))CustomCreateSubMatrices);
                  Mat aux = ff_to_PETSc(&dO);
                  IS perm;
                  ISSortPermutation(is, PETSC_TRUE, &perm);
                  O = new Mat;
                  MatPermute(aux, perm, perm, O);
                  ISDestroy(&perm);
                  MatDestroy(&aux);
                }
                PCSetUp(pc);
                ISDestroy(&loc);
                if (nnz) {
                  MatSetOperation(ptA->_petsc, MATOP_CREATE_SUBMATRICES, (void(*)(void))MatCreateSubMatrices);
                  delete O;
                }
                O = nullptr;
              }
#if defined(PETSC_HAVE_HPDDM) && defined(PETSC_USE_SHARED_LIBRARIES)
            }
#endif
            ISDestroy(&is);
          }
        }
        if (std::is_same< Type, Dmat >::value && (nargs[6] || nargs[11])) {
          if (nargs[2] && (nargs[5] || nargs[11])) {
            if (assembled) {
              if (ptParent) MatAssembled(ptParent->_petsc, &assembled);
              if (assembled) KSPSetUp(ksp);
            }
            PC pc;
            KSPGetPC(ksp, &pc);
            setCompositePC(pc, ptA->_vS);
          } else if (mpisize == 1 && nargs[10]) {
            IS is;
            std::vector< PetscInt > idx;
            KN< double >* pL = GetAny< KN< double >* >((*nargs[6])(stack));
            idx.reserve(pL->n);
            KNM< PetscScalar >* pS = GetAny< KNM< PetscScalar >* >((*nargs[10])(stack));
            for (int i = 0; i < pL->n; ++i)
              if (std::abs((*pL)[i]) > 1.0e-12) idx.emplace_back(i);
            ISCreateGeneral(PetscObjectComm((PetscObject)ptA->_ksp), idx.size( ), idx.data( ), PETSC_COPY_VALUES, &is);
            PC pc;
            KSPGetPC(ptA->_ksp, &pc);
            PCFactorSetUpMatSolverType(pc);
            Mat F, S;
            PCFactorGetMatrix(pc, &F);
            MatFactorSetSchurIS(F, is);
            KSPSetUp(ptA->_ksp);
            MatFactorSchurStatus status;
            MatFactorGetSchurComplement(F, &S, &status);
            PetscInt m, n;
            MatGetSize(S, &m, &n);
            pS->resize(m, n);
            PetscScalar* data;
            MatDenseGetArray(S, &data);
            std::copy_n(data, m * n, pS->operator PetscScalar*());
            MatDenseRestoreArray(S, &data);
            MatFactorRestoreSchurComplement(F, &S, status);
            ISDestroy(&is);
          }
        }
        if (std::is_same< Type, Dmat >::value && ptA->_petsc && nargs[17]) {
            PC pc;
            KSPGetPC(ksp, &pc);
            PetscBool isType;
            PetscObjectTypeCompare((PetscObject)pc, PCSHELL, &isType);
            User< LinearSolver< Type > > userPC = nullptr;
            if(isType) {
                PCShellGetContext(pc, &userPC);
                if(userPC)
                    delete userPC->op;
            }
            else {
                PCSetType(pc, PCSHELL);
            }
            if(!userPC)
                PetscNew(&userPC);
            const Polymorphic* op = dynamic_cast< const Polymorphic* >(nargs[17]);
            ffassert(op);
            const OneOperator* codeA = op->Find("(", ArrayOfaType(atype< KN< PetscScalar >* >( ), false));
            PetscInt n;
            MatGetLocalSize(ptA->_petsc, &n, NULL);
            userPC->op = new typename LinearSolver< Type >::MatF_O(n, stack, codeA);
            PCShellSetContext(pc, userPC);
            PCShellSetApply(pc, Op_User< LinearSolver< Type >, PC >);
            PCShellSetDestroy(pc, PCShellDestroy< LinearSolver< Dmat >  >);
        }
        const Polymorphic* op = nargs[19] ? dynamic_cast< const Polymorphic* >(nargs[19]) : nullptr;
        if (op) {
          ffassert(op);
          const OneOperator* codeM = op->Find(
            "(", ArrayOfaType(atype< long >( ), atype< double >( ), false));
          typename LinearSolver< Type >::MonF_O* mon = new typename LinearSolver< Type >::MonF_O(stack, codeM);
          KSPMonitorSet(ksp, Monitor< LinearSolver< Type > >, mon, MonitorDestroy< Type >);
        }
      }
    }
    if(nargs[18] && GetAny< bool >((*nargs[18])(stack)))
        KSPSetUp(ksp);
    return 0L;
  }

  template< class Type, unsigned short O >
  class view_Op : public E_F0mps {
   public:
    Expression A;
    static const int n_name_param = 3;
    static basicAC_F0::name_and_type name_param[];
    Expression nargs[n_name_param];
    view_Op(const basicAC_F0& args, Expression param1) : A(param1) {
      args.SetNameParam(n_name_param, name_param, nargs);
    }

    AnyType operator( )(Stack stack) const;
  };
  template< class Type, unsigned short O >
  basicAC_F0::name_and_type view_Op< Type, O >::name_param[] = {
    {!std::is_same<Type, KNM<PetscScalar>>::value && !O ? "object" : "communicator", !std::is_same<Type, KNM<PetscScalar>>::value && !O ? &typeid(std::string*) : &typeid(pcommworld)},
    {"format", &typeid(std::string*)},
    {"name", &typeid(std::string*)}
  };
  template< class Type, unsigned short O >
  class view : public OneOperator {
   public:
    view( ) : OneOperator(atype< long >( ), atype< Type* >( )) {}

    E_F0* code(const basicAC_F0& args) const {
      return new view_Op< Type, O >(args, t[0]->CastTo(args[0]));
    }
  };
  template< class Type, unsigned short O, typename std::enable_if<!std::is_same<Type, KNM<PetscScalar>>::value>::type* = nullptr >
  AnyType view_dispatched(Type* const& ptA, std::string const& o, std::string* const& type, std::string* const& name, MPI_Comm const& comm) {
    if (o.size() == 0 || o.compare("MAT") == 0) {
      bool pop = false;
      PetscViewer viewer = NULL;
      if(type) {
          if(type->compare("matlab") == 0) {
              if(name) PetscViewerASCIIOpen(!O ? PetscObjectComm((PetscObject)ptA->_petsc) : comm, name->c_str(), &viewer);
              else viewer = PETSC_VIEWER_STDOUT_(PetscObjectComm((PetscObject)ptA->_petsc));
              PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
              pop = true;
          }
          else if(type->compare("info") == 0) {
              ffassert(!O);
              if(name) PetscViewerASCIIOpen(!O ? PetscObjectComm((PetscObject)ptA->_petsc) : comm, name->c_str(), &viewer);
              else viewer = PETSC_VIEWER_STDOUT_(PetscObjectComm((PetscObject)ptA->_petsc));
              PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_INFO);
              pop = true;
          }
          else if(type->compare("binary") == 0) {
              if(name) PetscViewerBinaryOpen(!O ? PetscObjectComm((PetscObject)ptA->_petsc) : comm, name->c_str(), O ? FILE_MODE_READ : FILE_MODE_WRITE, &viewer);
              else viewer = PETSC_VIEWER_BINARY_(PetscObjectComm((PetscObject)ptA->_petsc));
              PetscViewerPushFormat(viewer, PETSC_VIEWER_NATIVE);
              pop = true;
          }
          else if(type->compare("draw") == 0) {
              ffassert(!O);
              MatView(ptA->_petsc, PETSC_VIEWER_DRAW_WORLD);
              return 0L;
          }
      }
      if(!viewer && name) PetscViewerASCIIOpen(PetscObjectComm((PetscObject)ptA->_petsc), name->c_str(), &viewer);
      if(!O)
          MatView(ptA->_petsc, viewer ? viewer : PETSC_VIEWER_STDOUT_(PetscObjectComm((PetscObject)ptA->_petsc)));
      else {
          ffassert(viewer);
          ptA->dtor( );
          MatCreate(comm, &ptA->_petsc);
          MatLoad(ptA->_petsc, viewer);
      }
      if(pop)
          PetscViewerPopFormat(PETSC_VIEWER_STDOUT_(PetscObjectComm((PetscObject)ptA->_petsc)));
      if(name)
          PetscViewerDestroy(&viewer);
    }
    else {
      ffassert(!O);
      if (ptA->_ksp) {
        if (o.compare("KSP") == 0)
          KSPView(ptA->_ksp, PETSC_VIEWER_STDOUT_(PetscObjectComm((PetscObject)ptA->_ksp)));
        else if (o.compare("PC") == 0) {
          PC pc;
          KSPGetPC(ptA->_ksp, &pc);
          PCView(pc, PETSC_VIEWER_STDOUT_(PetscObjectComm((PetscObject)ptA->_ksp)));
        }
      }
    }
    return 0L;
  }
  template< class Type, unsigned short O, typename std::enable_if<std::is_same<Type, KNM<PetscScalar>>::value>::type* = nullptr >
  AnyType view_dispatched(Type* const& ptA, std::string const& o, std::string* const& type, std::string* const& name, MPI_Comm const& comm) {
    bool pop = false;
    PetscViewer viewer = NULL;
    Mat A;
    if(!O)
        MatCreateDense(comm, ptA->N( ), PETSC_DECIDE, PETSC_DECIDE, ptA->M(), &ptA->operator( )(0, 0), &A);
    else
        MatCreate(comm, &A);
    if(type) {
        if(type->compare("matlab") == 0) {
            if(name) PetscViewerASCIIOpen(PetscObjectComm((PetscObject)A), name->c_str(), &viewer);
            else viewer = PETSC_VIEWER_STDOUT_(PetscObjectComm((PetscObject)A));
            PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
            pop = true;
        }
        else if(type->compare("info") == 0) {
            ffassert(!O);
            if(name) PetscViewerASCIIOpen(PetscObjectComm((PetscObject)A), name->c_str(), &viewer);
            else viewer = PETSC_VIEWER_STDOUT_(PetscObjectComm((PetscObject)A));
            PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_INFO);
            pop = true;
        }
        else if(type->compare("binary") == 0) {
            if(name) PetscViewerBinaryOpen(PetscObjectComm((PetscObject)A), name->c_str(), O ? FILE_MODE_READ : FILE_MODE_WRITE, &viewer);
            else viewer = PETSC_VIEWER_BINARY_(PetscObjectComm((PetscObject)A));
            PetscViewerPushFormat(viewer, PETSC_VIEWER_NATIVE);
            pop = true;
        }
        else if(type->compare("draw") == 0) {
            ffassert(!O);
            MatView(A, PETSC_VIEWER_DRAW_WORLD);
            MatDestroy(&A);
            return 0L;
        }
    }
    if(!viewer && name) PetscViewerASCIIOpen(PetscObjectComm((PetscObject)A), name->c_str(), &viewer);
    if(!O)
        MatView(A, viewer ? viewer : PETSC_VIEWER_STDOUT_(PetscObjectComm((PetscObject)A)));
    else {
        ffassert(viewer);
        MatSetType(A, MATDENSE);
        MatLoad(A, viewer);
        PetscInt m, N;
        MatGetLocalSize(A, &m, NULL);
        MatGetSize(A, NULL, &N);
        ptA->resize(m, N);
        PetscScalar* array;
        MatDenseGetArray(A, &array);
        if (array) std::copy_n(array, m * N, &(ptA->operator( )(0, 0)));
        MatDenseRestoreArray(A, &array);
    }
    if(pop)
        PetscViewerPopFormat(PETSC_VIEWER_STDOUT_(PetscObjectComm((PetscObject)A)));
    if(name)
        PetscViewerDestroy(&viewer);
    MatDestroy(&A);
    return 0L;
  }
  template< class Type, unsigned short O >
  AnyType view_Op< Type, O >::operator( )(Stack stack) const {
    Type* ptA = GetAny< Type* >((*A)(stack));
    std::string* object = !std::is_same<Type, KNM<PetscScalar>>::value && !O && nargs[0] ? GetAny< std::string* >((*nargs[0])(stack)) : NULL;
    MPI_Comm comm = (std::is_same<Type, KNM<PetscScalar>>::value || O) && nargs[0] ? *static_cast< MPI_Comm* >(GetAny< pcommworld >((*nargs[0])(stack))) : PETSC_COMM_WORLD;
    std::string o(object ? *object : "");
    std::transform(o.begin(), o.end(), o.begin(), ::toupper);
    std::string* type = nargs[1] ? GetAny< std::string* >((*nargs[1])(stack)) : NULL;
    std::string* name = nargs[2] ? GetAny< std::string* >((*nargs[2])(stack)) : NULL;
    return view_dispatched<Type, O>(ptA, o, type, name, comm);
  }

  template< class Type, template<class> class Storage >
  class changeNumbering : public OneOperator {
   public:
    const int c;
    class changeNumbering_Op : public E_F0mps {
     public:
      std::vector< std::pair< Expression, Expression > > E;
      std::vector< Expression > v;
      Expression A;
      Expression in;
      Expression out;
      const int c;
      static const int n_name_param = 2;
      static basicAC_F0::name_and_type name_param[];
      Expression nargs[n_name_param];
      changeNumbering_Op(const basicAC_F0& args, int d) : A(0), in(0), out(0), c(d) {
        args.SetNameParam(n_name_param, name_param, nargs);
        if (c != 2) {
          if (c == 1) {
            A = to< KN< long >* >(args[0]);
            in = to< Storage< PetscScalar >* >(args[1]);
          } else if (c != 4) {
            E.reserve(1);
            E.emplace_back(to< Type* >(args[0]), to< Storage< PetscScalar >* >(args[1]));
          } else {
            E.reserve(1);
            E.emplace_back(to< Type* >(args[0]), to< FEbaseArrayKn< PetscScalar >* >(args[1]));
          }
        } else {
          const E_Array* EA = dynamic_cast< const E_Array* >(args[0].LeftValue( ));
          const E_Array* Ex = dynamic_cast< const E_Array* >(args[1].LeftValue( ));
          ffassert(EA->size( ) == Ex->size( ) && EA->size( ));
          E.reserve(EA->size( ));
          v.reserve(EA->size( ));
          for (int i = 0; i < EA->size( ); ++i) {
            E.emplace_back(to< Type* >((*EA)[i]), to< Storage< PetscScalar >* >((*Ex)[i]));
            v.emplace_back(CastTo< KN_< PetscScalar > >((*Ex)[i]));
          }
        }
        if(c == 3)
            out = to< Storage< long >* >(args[2]);
        else
            out = to< Storage< PetscScalar >* >(args[2]);
      }

      AnyType operator( )(Stack stack) const;
      operator aType( ) const { return atype< long >( ); }
    };
    E_F0* code(const basicAC_F0& args) const { return new changeNumbering_Op(args, c); }
    changeNumbering( )
      : OneOperator(atype< long >( ), atype< Type* >( ), atype< Storage< PetscScalar >* >( ),
                    atype< Storage< PetscScalar >* >( )),
        c(0) {}
    changeNumbering(int)
      : OneOperator(atype< long >( ), atype< KN< long >* >( ), atype< Storage< PetscScalar >* >( ),
                    atype< Storage< PetscScalar >* >( )),
        c(1) {}
    changeNumbering(int, int)
      : OneOperator(atype< long >( ), atype< E_Array >( ), atype< E_Array >( ),
                    atype< Storage< PetscScalar >* >( )),
        c(2) {}
    changeNumbering(int, int, int)
      : OneOperator(atype< long >( ), atype< Type* >( ), atype< Storage< PetscScalar >* >( ),
                    atype< Storage< long >* >( )),
        c(3) {}
    changeNumbering(int, int, int, int)
      : OneOperator(atype< long >( ), atype< Type* >( ), atype< FEbaseArrayKn< PetscScalar >* >( ),
                    atype< Storage< PetscScalar >* >( )),
        c(4) {}
  };
  template< class Type, template<class> class Storage >
  basicAC_F0::name_and_type changeNumbering< Type, Storage >::changeNumbering_Op::name_param[] = {
    {"inverse", &typeid(bool)}, {"exchange", &typeid(bool)}};
  template< class Type, template<class> class Storage >
  AnyType changeNumbering< Type, Storage >::changeNumbering_Op::operator( )(Stack stack) const {
    Storage< long >* ptOutCast = (c == 3 ? GetAny< Storage< long >* >((*out)(stack)) : nullptr);
    Storage< PetscScalar >* ptOut = (c != 3 ? GetAny< Storage< PetscScalar >* >((*out)(stack)) : nullptr);
    if(c == 3) {
        ptOut = new Storage< PetscScalar >;
        PETSc::resize(ptOut, ptOutCast->N(), 1);
    }
    bool inverse = nargs[0] && GetAny< bool >((*nargs[0])(stack));
    if(inverse)
        ffassert(ptOut->operator PetscScalar*());
    int sum = 0;
    if (c == 0 || c == 2 || c == 3 || c == 4) {
      PetscScalar* pt = *ptOut;
      for (int j = 0; j < E.size( ); ++j) {
        Type* ptA = GetAny< Type* >((*(E[j].first))(stack));
        if (ptA && (ptA->_last - ptA->_first || (inverse && ptA->_num))) {
          PetscInt m;
          ffassert(ptA->_petsc);
          MatGetLocalSize(ptA->_petsc, &m, NULL);
          if (c != 4) {
            Storage< PetscScalar >* ptIn = GetAny< Storage< PetscScalar >* >((*(E[j].second))(stack));
            if (!inverse) ffassert(ptIn->N() == ptA->_A->getDof());
            if (c != 2) {
              if (inverse) ffassert(ptOut->N() == m);
              changeNumbering_func(ptA->_num, ptA->_first, ptA->_last, m, ptA->_A->getDof(),
                                   1, ptIn, ptOut, inverse);
            } else {
              sum += m;
              if (ptOut->N() < sum && !inverse) {
                ffassert(ptIn->M() == 1);
                ffassert(ptOut->M() == 1);
                PETSc::resize(ptOut, sum, ptIn->M());
                pt = *ptOut + sum - m;
              }
              KN_< PetscScalar > ptOutShift(pt, m);
              changeNumbering_func(ptA->_num, ptA->_first, ptA->_last, m, ptA->_A->getDof(),
                                   1, ptIn, &ptOutShift, inverse);
            }
            if (inverse && nargs[1] && GetAny< bool >((*nargs[1])(stack))) {
              for(int i = 0; i < ptIn->M(); ++i)
                ptA->_A->exchange(*ptIn + i * ptIn->N());
            }
            if (c == 2) {
              pt += m;
              if (j == E.size( ) - 1)
                ffassert(ptOut->N() == std::distance(static_cast< PetscScalar* >(*ptOut), pt) ||
                         (!inverse && ptOut->N() == sum));
              KN_<PetscScalar> view = GetAny<KN_<PetscScalar>>((*(v[j]))(stack));
              if(view.operator PetscScalar*() != ptIn->operator PetscScalar*())
                view = *ptIn;
            }
          } else {
            FEbaseArrayKn< PetscScalar >* ptIn = GetAny< FEbaseArrayKn< PetscScalar >* >((*(E[j].second))(stack));
            ffassert(ptIn && ptIn->N == ptOut->M());
            ffassert(ptIn->get(0) && ptIn->get(0)->n == ptA->_A->getDof());
            ffassert(ptOut->N() == m);
            for(int i = 0; i < ptIn->N; ++i) {
                KN_< PetscScalar > out(pt + i * m, m);
                changeNumbering_func(ptA->_num, ptA->_first, ptA->_last, m, ptA->_A->getDof(),
                                     1, ptIn->get(i), &out, inverse);
            }
            if (inverse && nargs[1] && GetAny< bool >((*nargs[1])(stack))) {
              for(int i = 0; i < ptIn->N; ++i)
                ptA->_A->exchange(&(*ptIn->get(i))[0]);
            }
          }
        }
      }
      if(c == 3) {
          PETSc::resize(ptOutCast, ptOut->N(), 1);
          PetscScalar* pt = *ptOut;
          long* ptCast = *ptOutCast;
          for(int i = 0; i < ptOut->N(); ++i)
              ptCast[i] = std::lround(std::real(pt[i]));
          delete ptOut;
      }
    } else {
      Storage< PetscScalar >* ptIn = GetAny< Storage< PetscScalar >* >((*in)(stack));
      ffassert(ptIn->M() == 1);
      ffassert(ptOut->M() == 1);
      KN< long >* ptA = GetAny< KN< long >* >((*A)(stack));
      PetscInt* num = reinterpret_cast< PetscInt* >(ptA->operator long*());
      changeNumbering_func(num + 2, num[0], num[1], num[1] - num[0], ptA->n - 2, 1, ptIn, ptOut,
                           inverse);
    }
    return 0L;
  }

  template< char N, class Type >
  long MatMult(Type* const& A, KN< PetscScalar >* const& in, KN< PetscScalar >* const& out) {
    if (A) {
      Vec x, y;
      PetscInt size;
      if (N == 'T' || N == 'H') {
        MatCreateVecs(A->_petsc, &y, &x);
        VecGetLocalSize(x, &size);
        ffassert(in->n == size);
        VecGetLocalSize(y, &size);
        VecPlaceArray(x, *in);
        out->resize(size);
        VecPlaceArray(y, *out);
        if (N == 'T')
          MatMultTranspose(A->_petsc, x, y);
        else
          MatMultHermitianTranspose(A->_petsc, x, y);
      } else {
        MatCreateVecs(A->_petsc, &x, &y);
        VecGetLocalSize(x, &size);
        ffassert(in->n == size);
        VecGetLocalSize(y, &size);
        VecPlaceArray(x, *in);
        out->resize(size);
        VecPlaceArray(y, *out);
        MatMult(A->_petsc, x, y);
      }
      VecResetArray(y);
      VecResetArray(x);
      VecDestroy(&y);
      VecDestroy(&x);
    }
    return 0L;
  }
  template< char P, class Type, class Storage, typename std::enable_if< std::is_same< Storage, KNM< PetscScalar > >::value >::type* = nullptr  >
  long MatMatMult(Type* const& A, Storage* const& in, Storage* const& out) {
    if (A->_petsc) {
      Mat x, y;
      PetscInt n, m, N, M;
      MatGetLocalSize(A->_petsc, &n, &m);
      MatGetSize(A->_petsc, &N, &M);
      MatCreateDense(PetscObjectComm((PetscObject)A->_petsc), P == 'T' || P == 'H' ? n : m, PETSC_DECIDE,
                     P == 'T' || P == 'H' ? N : M, in->M( ), &(in->operator( )(0, 0)), &x);
      if (P == 'T' || P == 'H') {
        ffassert(in->N( ) == n);
        out->resize(m, in->M( ));
        if (P == 'T')
          MatTransposeMatMult(A->_petsc, x, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &y);
        else {
          Mat w;
          MatDuplicate(x, MAT_COPY_VALUES, &w);
          MatConjugate(w);
          MatTransposeMatMult(A->_petsc, w, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &y);
          MatDestroy(&w);
          MatConjugate(y);
        }
        PetscScalar* array;
        MatDenseGetArray(y, &array);
        if (array) std::copy_n(array, m * in->M( ), &(out->operator( )(0, 0)));
        MatDenseRestoreArray(y, &array);
      } else {
        ffassert(in->N( ) == m);
        out->resize(n, in->M( ));
        MatCreateDense(PetscObjectComm((PetscObject)A->_petsc), P == 'T' || P == 'H' ? m : n, PETSC_DECIDE,
                       P == 'T' || P == 'H' ? M : N, in->M( ), &(out->operator( )(0, 0)), &y);
        MatMatMult(A->_petsc, x, MAT_REUSE_MATRIX, PETSC_DEFAULT, &y);
      }
      MatDestroy(&y);
      MatDestroy(&x);
    }
    return 0L;
  }
  template< char P, class Type, class Storage, typename std::enable_if< std::is_same< Storage, Dmat >::value >::type* = nullptr  >
  long MatMatMult(Type* const& A, Storage* const& in, Storage* const& out) {
    ffassert(out != in && out != A);
    if (out->_petsc) out->dtor( );
    PetscBool isDense;
    PetscObjectTypeCompareAny((PetscObject)in->_petsc, &isDense, MATMPIDENSE, MATSEQDENSE, "");
    if (isDense) {
        PetscInt m, M, n;
        MatGetLocalSize(A->_petsc, &m, NULL);
        MatGetSize(A->_petsc, &M, NULL);
        MatGetLocalSize(in->_petsc, &n, NULL);
        MatCreateDense(PetscObjectComm((PetscObject)A->_petsc), m, n, M, PETSC_DETERMINE, NULL, &out->_petsc);
        MatMatMult(A->_petsc, in->_petsc, MAT_REUSE_MATRIX, PETSC_DEFAULT, &out->_petsc);
    }
    else
        MatMatMult(A->_petsc, in->_petsc, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &out->_petsc);
    return 0L;
  }
  template< class Type >
  long MatPtAP(Type* const& A, Type* const& P, Type* const& B) {
    if (A->_petsc && P->_petsc) {
      B->dtor( );
      MatPtAP(A->_petsc, P->_petsc, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &B->_petsc);
    }
    return 0L;
  }
  void prepareConvert(Mat A, Mat* B) {
    MatType type;
    PetscBool isType;
    Mat** mat;
    PetscInt M, N;
    MatNestGetSubMats(A, &M, &N, &mat);
    std::vector< std::pair< std::pair< PetscInt, PetscInt >, Mat > > b;
    b.reserve(M * N);
    for (PetscInt i = 0; i < M; ++i) {
      for (PetscInt j = 0; j < N; ++j) {
        if (mat[i][j]) {
          MatGetType(mat[i][j], &type);
          PetscStrcmp(type, MATHERMITIANTRANSPOSEVIRTUAL, &isType);
          if (isType) {
            b.emplace_back(std::make_pair(std::make_pair(i, j), Mat( )));
            Mat D = mat[i][j];
            Mat C;
            MatHermitianTransposeGetMat(D, &b.back( ).second);
            MatHermitianTranspose(b.back( ).second, MAT_INITIAL_MATRIX, &C);
            PetscObjectReference((PetscObject)b.back().second);
            MatNestSetSubMat(A, i, j, C);
            MatDestroy(&C);
          }
        }
      }
    }
    MatConvert(A, MATAIJ, MAT_INITIAL_MATRIX, B);
    for (std::pair< std::pair< PetscInt, PetscInt >, Mat > p : b) {
      Mat C;
      MatCreateHermitianTranspose(p.second, &C);
      MatDestroy(&p.second);
      MatNestSetSubMat(A, p.first.first, p.first.second, C);
      MatDestroy(&C);
    }
  }
  template< class Type >
  class convert_Op : public E_F0mps {
   public:
    Expression A;
    Expression B;
    static const int n_name_param = 1;
    static basicAC_F0::name_and_type name_param[];
    Expression nargs[n_name_param];
    convert_Op(const basicAC_F0& args, Expression param1, Expression param2) : A(param1), B(param2) {
      args.SetNameParam(n_name_param, name_param, nargs);
    }

    AnyType operator( )(Stack stack) const;
  };
  template< class Type >
  basicAC_F0::name_and_type convert_Op< Type >::name_param[] = {
    {"type", &typeid(std::string*)}
  };
  template< class Type >
  class convert : public OneOperator {
   public:
    convert( ) : OneOperator(atype< long >( ), atype< Type* >( ), atype< Type* >( )) {}

    E_F0* code(const basicAC_F0& args) const {
      return new convert_Op< Type >(args, t[0]->CastTo(args[0]), t[1]->CastTo(args[1]));
    }
  };
  template< class Type >
  AnyType convert_Op< Type >::operator( )(Stack stack) const {
    Type* ptA = GetAny< Type* >((*A)(stack));
    Type* ptB = GetAny< Type* >((*B)(stack));
    std::string stype = nargs[0] ? *GetAny< std::string* >((*nargs[0])(stack)) : "aij";
    if (nargs[0]) std::transform(stype.begin(), stype.end(), stype.begin(), ::tolower);
    if (ptA->_petsc) {
      if (ptB->_petsc) ptB->dtor( );
      MatType type;
      PetscBool isType;
      MatGetType(ptA->_petsc, &type);
      PetscStrcmp(type, MATNEST, &isType);
      if (isType) {
        prepareConvert(ptA->_petsc, &ptB->_petsc);
      }
      else {
        Mat C = ptA->_petsc;
#if defined(PETSC_HAVE_HTOOL)
        PetscStrcmp(type, MATHTOOL, &isType);
        if (isType) {
          MatConvert(ptA->_petsc, MATDENSE, MAT_INITIAL_MATRIX, &ptB->_petsc);
          if (nargs[0] && stype.compare("dense") != 0) {
            C = ptB->_petsc;
            isType = PETSC_FALSE;
          }
        }
#endif
        if (!isType) {
          MatType mtype(stype.c_str());
          MatConvert(C, mtype, MAT_INITIAL_MATRIX, &ptB->_petsc);
          if (C != ptA->_petsc)
              MatDestroy(&C);
        }
      }
    }
    return 0L;
  }
  template< class Type >
  long GetConvergedReason(Type* const& A) {
    if (A->_ksp) {
      KSPConvergedReason reason;
      KSPGetConvergedReason(A->_ksp, &reason);
      return static_cast<long>(reason);
    }
    return 0L;
  }
  template< class Type >
  long GetIterationNumber(Type* const& A) {
    if (A->_ksp) {
      PetscInt its;
      KSPGetIterationNumber(A->_ksp, &its);
      return static_cast<long>(its);
    }
    return 0L;
  }
  template< class Type >
  long SetResidualHistory(Type* const& A, KN< double >* const& hist) {
    if (A->_ksp) {
      KSPSetResidualHistory(A->_ksp, hist->operator double*(), hist->n, PETSC_TRUE);
    }
    return 0L;
  }
  template< class Type >
  long MatZeroRows(Type* const& A, KN< double >* const& ptRows) {
    if (A->_petsc) {
      PetscInt m;
      MatGetLocalSize(A->_petsc, &m, NULL);
      ffassert(ptRows->n == m);
      std::vector< PetscInt > rows;
      rows.reserve(m);
      PetscInt start;
      MatGetOwnershipRange(A->_petsc, &start, NULL);
      for (int i = 0; i < m; ++i) {
        if (std::abs(ptRows->operator[](i)) > 1.0e-12) rows.emplace_back(start + i);
      }
      MatSetOption(A->_petsc, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
      MatSetOption(A->_petsc, MAT_NO_OFF_PROC_ZERO_ROWS, PETSC_TRUE);
      MatZeroRows(A->_petsc, rows.size( ), rows.data( ), 0.0, NULL, NULL);
    }
    return 0L;
  }
  template< class Type >
  class LinearSolver : public OneOperator {
   public:
    typedef KN< PetscScalar > Kn;
    typedef KN_< PetscScalar > Kn_;
    class MatF_O : public RNM_VirtualMatrix< PetscScalar > {
     public:
      Stack stack;
      mutable Kn x;
      C_F0 c_x;
      Expression mat, matT;
      MatF_O(int n, Stack stk, const OneOperator* op, int m = -1, const OneOperator* opT = nullptr)
        : RNM_VirtualMatrix< PetscScalar >(n, m), stack(stk), x(n), c_x(CPValue(x)),
          mat(op ? CastTo< Kn_ >(C_F0(op->code(basicAC_F0_wa(c_x)), (aType)*op)) : 0), matT(opT ? CastTo< Kn_ >(C_F0(opT->code(basicAC_F0_wa(c_x)), (aType)*opT)) : 0) {}
      ~MatF_O( ) {
        delete matT;
        delete mat;
        Expression zzz = c_x;
        delete zzz;
      }
      void addMatMul(const Kn_& xx, Kn_& Ax) const {
        if(x.N( ) != xx.N( ))
            x.resize(xx.N( ));
        x = xx;
        ffassert(mat);
        Ax += GetAny< Kn_ >((*mat)(stack));
        WhereStackOfPtr2Free(stack)->clean( );
      }
      void addMatTransMul(const Kn_& xx, Kn_& Ax) const {
        if(x.N( ) != xx.N( ))
            x.resize(xx.N( ));
        x = xx;
        ffassert(matT);
        Ax += GetAny< Kn_ >((*matT)(stack));
        WhereStackOfPtr2Free(stack)->clean( );
      }
    };
    class MonF_O {
     public:
      Stack stack;
      mutable long s;
      C_F0 c_s;
      mutable double t;
      C_F0 c_t;
      Expression mat;
      MonF_O(Stack stk, const OneOperator* op)
        : stack(stk), s(0), c_s(CPValue(s)), t(0), c_t(CPValue(t)),
          mat(op ? CastTo< long >(C_F0(op->code(basicAC_F0_wa({c_s, c_t})), (aType)*op)) : 0) {
      }
      ~MonF_O( ) {
        delete mat;
        Expression zzz = c_s;
        delete zzz;
        zzz = c_t;
        delete zzz;
      }
      long apply(const long& ss, const double& tt) const {
        s = ss;
        t = tt;
        ffassert(mat);
        long ret = GetAny< long >((*mat)(stack));
        WhereStackOfPtr2Free(stack)->clean( );
        return ret;
      }
    };
    const int c;
    class E_LinearSolver : public E_F0mps {
     public:
      Expression A;
      Expression x;
      Expression y;
      const OneOperator *codeA, *codeC;
      const int c;
      static const int n_name_param = 2;
      static basicAC_F0::name_and_type name_param[];
      Expression nargs[n_name_param];
      E_LinearSolver(const basicAC_F0& args, int d) : A(0), x(0), y(0), codeA(0), codeC(0), c(d) {
        args.SetNameParam(n_name_param, name_param, nargs);
        if (c != 2) {
          if (c == 1) {
            const Polymorphic* op = dynamic_cast< const Polymorphic* >(args[0].LeftValue( ));
            ffassert(op);
            codeA = op->Find("(", ArrayOfaType(atype< KN< PetscScalar >* >( ), false));
            if (nargs[0]) {
              op = dynamic_cast< const Polymorphic* >(nargs[0]);
              ffassert(op);
              codeC = op->Find("(", ArrayOfaType(atype< Kn* >( ), false));
            }
          } else {
            A = to< Type* >(args[0]);
          }
          if (c != 4) {
            x = to< KN< PetscScalar >* >(args[1]);
            y = to< KN< PetscScalar >* >(args[2]);
          } else {
            x = to< Type* >(args[1]);
            y = to< Type* >(args[2]);
          }
        } else {
          A = to< Type* >(args[0]);
          x = to< KNM< PetscScalar >* >(args[1]);
          y = to< KNM< PetscScalar >* >(args[2]);
        }
      }

      AnyType operator( )(Stack stack) const;
      operator aType( ) const { return atype< long >( ); }
    };
    E_F0* code(const basicAC_F0& args) const { return new E_LinearSolver(args, c); }
    LinearSolver( )
      : OneOperator(atype< long >( ), atype< Type* >( ), atype< KN< PetscScalar >* >( ),
                    atype< KN< PetscScalar >* >( )),
        c(0) {}
    LinearSolver(int)
      : OneOperator(atype< long >( ), atype< Polymorphic* >( ), atype< KN< PetscScalar >* >( ),
                    atype< KN< PetscScalar >* >( )),
        c(1) {}
    LinearSolver(int, int)
      : OneOperator(atype< long >( ), atype< Type* >( ), atype< KNM< PetscScalar >* >( ),
                    atype< KNM< PetscScalar >* >( )),
        c(2) {}
    LinearSolver(int, int, int)
      : OneOperator(atype< long >( ), atype< Type* >( ), atype< KN< PetscScalar >* >( ),
                    atype< KN< PetscScalar >* >( )),
        c(3) {}
    LinearSolver(int, int, int, int)
      : OneOperator(atype< long >( ), atype< Type* >( ), atype< Type* >( ),
                    atype< Type* >( )),
        c(4) {}
    LinearSolver(int, int, int, int, int)
      : OneOperator(atype< long >( ), atype< Type* >( ), atype< KN< PetscScalar >* >( ),
                    atype< KN< PetscScalar >* >( )),
        c(5) {}
  };
  template< class Type >
  basicAC_F0::name_and_type LinearSolver< Type >::E_LinearSolver::name_param[] = {
    {"precon", &typeid(Polymorphic*)}, {"sparams", &typeid(std::string*)}};
  template< class Type >
  AnyType LinearSolver< Type >::E_LinearSolver::operator( )(Stack stack) const {
    if (c != 2) {
      KN< PetscScalar >* in, *out;
      if (c != 4) {
        in = GetAny< KN< PetscScalar >* >((*x)(stack));
        out = GetAny< KN< PetscScalar >* >((*y)(stack));
      }
      if (A) {
        Type* ptA = GetAny< Type* >((*A)(stack));
        if (!ptA->_ksp) {
          ffassert(ptA->_petsc);
          KSPCreate(PetscObjectComm((PetscObject)ptA->_petsc), &ptA->_ksp);
          KSPSetOperators(ptA->_ksp, ptA->_petsc, ptA->_petsc);
          KSPSetFromOptions(ptA->_ksp);
        }
        if (c != 4) {
          Vec x = NULL, y;
          MatCreateVecs(ptA->_petsc, &y, in == out ? NULL : &x);
          PetscInt size;
          if (x) {
            VecGetLocalSize(x, &size);
            ffassert(in->n == size);
          }
          VecGetLocalSize(y, &size);
          if (!x)
            ffassert(in->n == size);
          if (out->n != size) {
            out->resize(size);
            *out = PetscScalar( );
          }
          if (x)
            VecPlaceArray(x, *in);
          else {
            x = y;
            PetscObjectReference((PetscObject)y);
          }
          VecPlaceArray(y, *out);
          PetscInt N, rbegin;
          PetscScalar* tmpIn, *tmpOut;
          if (c != 3 && c != 5)
            KSPSolve(ptA->_ksp, x, y);
          else {
            ffassert(in != out);
            if (c == 3) VecConjugate(x);
            KSPSolveTranspose(ptA->_ksp, x, y);
            if (c == 3) VecConjugate(y);
          }
          if (x != y) {
            if (c == 3) VecConjugate(x);
            VecResetArray(x);
            VecDestroy(&x);
          }
          VecResetArray(y);
          VecDestroy(&y);
        } else {
          Dmat* X = GetAny< Dmat* >((*x)(stack));
          Dmat* Y = GetAny< Dmat* >((*y)(stack));
          ffassert(X && Y && X != Y);
          PetscBool isDense;
          PetscObjectTypeCompareAny((PetscObject)X->_petsc, &isDense, MATMPIDENSE, MATSEQDENSE, "");
          Mat XD;
          if (!isDense) {
            MatConvert(X->_petsc, MATDENSE, MAT_INITIAL_MATRIX, &XD);
          } else XD = X->_petsc;
          if (!Y->_petsc) MatDuplicate(XD, MAT_DO_NOT_COPY_VALUES, &Y->_petsc);
          KSPMatSolve(ptA->_ksp, XD, Y->_petsc);
          if (!isDense) MatDestroy(&XD);
        }
      } else {
        User< LinearSolver< Type > > user = nullptr;
        PetscNew(&user);
        user->op = new LinearSolver< Type >::MatF_O(in->n, stack, codeA);
        Mat S;
        MatCreateShell(PETSC_COMM_WORLD, in->n, in->n, PETSC_DECIDE, PETSC_DECIDE, user, &S);
        MatShellSetOperation(S, MATOP_MULT,
                             (void (*)(void))Op_User< LinearSolver< Type >, Mat >);
        Vec x, y;
        MatCreateVecs(S, &x, &y);
        if (out->n != in->n) {
          out->resize(in->n);
          *out = PetscScalar( );
        }
        VecPlaceArray(x, *in);
        VecPlaceArray(y, *out);
        KSP ksp;
        KSPCreate(PETSC_COMM_WORLD, &ksp);
        KSPSetOperators(ksp, S, S);
        if (nargs[1]) {
          std::string* options = GetAny< std::string* >((*nargs[1])(stack));
          PetscOptionsInsertString(NULL, options->c_str());
        }
        PC pc;
        KSPGetPC(ksp, &pc);
        User< LinearSolver< Type > > userPC = nullptr;
        if (codeC) {
          PCSetType(pc, PCSHELL);
          PetscNew(&userPC);
          userPC->op = new LinearSolver< Type >::MatF_O(in->n, stack, codeC);
          PCShellSetContext(pc, userPC);
          PCShellSetApply(pc, Op_User< LinearSolver< Type >, PC >);
        }
        KSPSetFromOptions(ksp);
        KSPSolve(ksp, x, y);
        VecResetArray(y);
        VecResetArray(x);
        VecDestroy(&y);
        VecDestroy(&x);
        MatDestroy(&S);
        KSPDestroy(&ksp);
        if (codeC) {
          delete userPC->op;
          PetscFree(userPC);
        }
        delete user->op;
        PetscFree(user);
      }
    } else {
      KNM< PetscScalar >* in = GetAny< KNM< PetscScalar >* >((*x)(stack));
      KNM< PetscScalar >* out = GetAny< KNM< PetscScalar >* >((*y)(stack));
      if (A) {
        Type* ptA = GetAny< Type* >((*A)(stack));
        PetscBool isType;
        PetscInt m, M;
        MatGetLocalSize(ptA->_petsc, &m, NULL);
        ffassert(in->N( ) == m);
        out->resize(m, in->M( ));
        if (!ptA->_ksp) {
          KSPCreate(PetscObjectComm((PetscObject)ptA->_petsc), &ptA->_ksp);
          KSPSetOperators(ptA->_ksp, ptA->_petsc, ptA->_petsc);
          isType = PETSC_FALSE;
        }
        else {
          KSPType type;
          KSPGetType(ptA->_ksp, &type);
#if defined(KSPHPDDM)
          PetscStrcmp(type, KSPHPDDM, &isType);
#else
          isType = PETSC_FALSE;
#endif
        }
        isType = PETSC_TRUE;
#if defined(PETSC_HAVE_HPDDM)
        if(isType) {
          MatGetSize(ptA->_petsc, &M, NULL);
          Mat B, C;
          MatCreateDense(PetscObjectComm((PetscObject)ptA->_ksp), in->N( ), PETSC_DECIDE, M, in->M(), &in->operator( )(0, 0), &B);
          MatCreateDense(PetscObjectComm((PetscObject)ptA->_ksp), out->N( ), PETSC_DECIDE, M, out->M(), &out->operator( )(0, 0), &C);
          KSPMatSolve(ptA->_ksp, B, C);
          MatDestroy(&C);
          MatDestroy(&B);
        }
        else
            ffassert(0);
#else
        HPDDM::PETScOperator op(ptA->_ksp, m);
        op.apply(&in->operator( )(0, 0), &out->operator( )(0, 0), in->M( ));
#endif
      }
    }
    return 0L;
  }
  template< class HpddmType, int C >
  AnyType initCSRfromDMatrix< HpddmType, C >::initCSRfromDMatrix_Op::operator( )(
    Stack stack) const {
    ffassert((C == 0 && (c == 0 || c == 2)) || (C == c && C == 1));
    DistributedCSR< HpddmType >* ptA = GetAny< DistributedCSR< HpddmType >* >((*A)(stack));
    DistributedCSR< HpddmType >* ptB = GetAny< DistributedCSR< HpddmType >* >((*B)(stack));
    Matrice_Creuse< upscaled_type<PetscScalar> >* ptK = (c == 0 ? GetAny< Matrice_Creuse< upscaled_type<PetscScalar> >* >((*K)(stack)) : nullptr);
    if (ptB->_A) {
      HPDDM::MatrixCSR< PetscScalar >* dA;
      if (c == 0 && ptK->A) {
        MatriceMorse< upscaled_type<PetscScalar> >* mA = static_cast< MatriceMorse< upscaled_type<PetscScalar> >* >(&(*ptK->A));
        dA = new_HPDDM_MatrixCSR< PetscScalar >(mA);
      } else {
        if (c == 0) {
          int* ia = new int[1]( );
          dA = new HPDDM::MatrixCSR< PetscScalar >(0, 0, 0, nullptr, ia, nullptr, false, true);
        } else {
          int m = ptB->_A->getDof();
          int* ia = (c == 1 ? new int[m + 1]( ) : nullptr);
          dA = new HPDDM::MatrixCSR< PetscScalar >(m, m, 0, nullptr, ia, nullptr, false, true);
        }
      }
      Matrice_Creuse< double >* pList =
        nargs[2] && c != 1 ? GetAny< Matrice_Creuse< double >* >((*nargs[2])(stack)) : nullptr;
      PetscInt bs;
      MatGetBlockSize(ptB->_petsc, &bs);
      if (!pList) {
        ptA->_A = new HpddmType(static_cast< const HPDDM::Subdomain< PetscScalar >& >(*ptB->_A));
        ptA->_A->setMatrix(dA);
        ptA->_num = new PetscInt[dA->HPDDM_n];
        std::copy_n(ptB->_num, dA->HPDDM_n, ptA->_num);
        ptA->_first = ptB->_first;
        ptA->_last = ptB->_last;
        initPETScStructure<false>(ptA, bs,
          nargs[1] && GetAny< bool >((*nargs[1])(stack)) ? PETSC_TRUE : PETSC_FALSE,
          static_cast< KN< double >* >(nullptr));
      } else {
        int n = ptB->_A->getDof();
        ffassert(dA->HPDDM_n == n);
        HPDDM::MatrixCSR< void >* L;
        KN< double >* empty = nullptr;
        if (pList->A) {
          MatriceMorse< double >* mList = static_cast< MatriceMorse< double >* >(&*(pList->A));
          mList->CSR( );
          ffassert(mList->n == mList->nnz);
          ffassert(mList->m == n);
          L = new HPDDM::MatrixCSR< void >(mList->n, n, mList->n, mList->p, mList->j, false);
          ptA->_D = new KN<PetscReal>(mList->n);
          for (int i = 0; i < mList->n; ++i) ptA->_D->operator[](i) = ptB->_A->getScaling()[mList->j[i]];
          empty = new KN< double >(n, (double*)(ptB->_A->getScaling()));
        } else {
          L = new HPDDM::MatrixCSR< void >(0, n, 0, nullptr, nullptr, false);
          empty = new KN< double >(0);
        }
        const HPDDM::vectorNeighbor& map = ptB->_A->getMap();
        std::vector<int> o;
        std::vector<std::vector<int>> r;
        o.reserve(map.size());
        r.reserve(map.size());
        for(const auto& i : map) {
            o.emplace_back(i.first);
            r.emplace_back(i.second);
        }
        const MPI_Comm& comm = ptB->_A->getCommunicator();
        ptA->_A = new HpddmType;
        ptA->_A->HPDDM::template Subdomain< PetscScalar >::initialize(
          dA, o, r, const_cast<MPI_Comm*>(&comm), L);
        delete L;
        ptA->_num = new PetscInt[ptA->_A->getDof()];
        initPETScStructure<false>(ptA, bs,
          nargs[1] && GetAny< bool >((*nargs[1])(stack)) ? PETSC_TRUE : PETSC_FALSE, empty);
        delete empty;
      }
      if (c == 1) {
        MatSetType(ptA->_petsc, MATSHELL);
        User< LinearSolver< Dmat > > user = nullptr;
        PetscNew(&user);
        user->op = nullptr;
        const Polymorphic* op = nargs[0] ? dynamic_cast< const Polymorphic* >(nargs[0]) : nullptr;
        if(op) {
            const OneOperator* codeAt = op->Find("(", ArrayOfaType(atype< KN< PetscScalar >* >( ), false));
            if (codeAt) {
              user->op = new LinearSolver< Dmat >::MatF_O(ptB->_last - ptB->_first, stack, codeA, -1, codeAt);
              MatShellSetOperation(ptA->_petsc, MATOP_MULT_TRANSPOSE,
                      (void (*)(void))Op_User< LinearSolver< Mat >, Mat, 'T' >);
            }
        }
        if(!user->op) user->op = new LinearSolver< Dmat >::MatF_O(ptB->_last - ptB->_first, stack, codeA);
        MatShellSetContext(ptA->_petsc, user);
        MatShellSetOperation(ptA->_petsc, MATOP_MULT,
                             (void (*)(void))Op_User< LinearSolver< Mat >, Mat >);
        MatShellSetOperation(ptA->_petsc, MATOP_DESTROY, (void (*)(void))ShellDestroy< LinearSolver< Dmat >  >);
        MatSetUp(ptA->_petsc);
      }
    } else if(ptB->_petsc) {
      MatDuplicate(ptB->_petsc, MAT_COPY_VALUES, &ptA->_petsc);
    }
    if (c == 0 && nargs[0] && GetAny< bool >((*nargs[0])(stack))) ptK->destroy( );
    return ptA;
  }
  template< class HpddmType, int D >
  AnyType initRectangularCSRfromDMatrix< HpddmType, D >::initRectangularCSRfromDMatrix_Op::operator( )(
    Stack stack) const {
    ffassert((D == 0 && (c == 0 || c == 1 || c == 3)) || (c == 2 && D == 1));
    DistributedCSR< HpddmType >* ptA = GetAny< DistributedCSR< HpddmType >* >((*A)(stack));
    DistributedCSR< HpddmType >* ptB = GetAny< DistributedCSR< HpddmType >* >((*B)(stack));
    DistributedCSR< HpddmType >* ptC = GetAny< DistributedCSR< HpddmType >* >((*C)(stack));
    Matrice_Creuse< upscaled_type<PetscScalar> >* ptK =
      (c == 0 || c == 3 ? GetAny< Matrice_Creuse< upscaled_type<PetscScalar> >* >((*K)(stack)) : nullptr);
    if (ptB->_A && ptC->_A) {
      ptA->_first = ptB->_first;
      ptA->_last = ptB->_last;
      ptA->_cfirst = ptC->_first;
      ptA->_clast = ptC->_last;
      PetscInt bsB, bsC;
      MatGetBlockSize(ptB->_petsc, &bsB);
      MatGetBlockSize(ptC->_petsc, &bsC);
      MatCreate(PetscObjectComm((PetscObject)ptB->_petsc), &ptA->_petsc);
      if (bsB == bsC && bsB > 1) MatSetBlockSize(ptA->_petsc, bsB);
      MatSetSizes(ptA->_petsc, ptB->_last - ptB->_first, ptC->_last - ptC->_first, PETSC_DECIDE, PETSC_DECIDE);
      MatSetType(ptA->_petsc, MATAIJ);
      if (c == 0 || c == 3) {
        PetscInt* ia = nullptr;
        PetscInt* ja = nullptr;
        PetscScalar* a = nullptr;
        bool free = true;
        if (ptK->A) {
          MatriceMorse< upscaled_type<PetscScalar> >* mA = static_cast< MatriceMorse< upscaled_type<PetscScalar> >* >(&(*ptK->A));
          ff_HPDDM_MatrixCSR< PetscScalar > dA(mA);
          ptA->_num = new PetscInt[mA->n + (ptC->_A && ptC->_A->getMatrix() ? ptC->_A->getMatrix()->HPDDM_m : mA->m)];
          ptA->_cnum = ptA->_num + mA->n;
          std::copy_n(ptB->_num, mA->n, ptA->_num);
          std::copy_n(ptC->_num, (ptC->_A && ptC->_A->getMatrix() ? ptC->_A->getMatrix()->HPDDM_m : mA->m), ptA->_cnum);
          KN<PetscScalar>* numbering = nargs[1] ? GetAny< KN<PetscScalar>* >((*nargs[1])(stack)) : NULL;
          if (c == 0 || !numbering)
            free = HPDDM::template Subdomain< PetscScalar >::distributedCSR(
              ptA->_num, ptA->_first, ptA->_last, ia, ja, a, &dA, ptA->_num + mA->n);
          else {
              KN<PetscInt> col(mA->m);
              ffassert(mA->m == numbering->N());
              for(int i = 0; i < col.N(); ++i)
                  col[i] = std::lround(std::real(numbering->operator[](i)));
              free = HPDDM::template Subdomain< PetscScalar >::distributedCSR(
                ptA->_num, ptA->_first, ptA->_last, ia, ja, a, &dA, col.operator PetscInt*());
          }
        } else
          ia = new PetscInt[ptB->_last - ptB->_first + 1]( );
        MatSeqAIJSetPreallocationCSR(ptA->_petsc, ia, ja, a);
        MatMPIAIJSetPreallocationCSR(ptA->_petsc, ia, ja, a);
        if (free) {
          delete[] ia;
          delete[] ja;
          delete[] a;
        }
      } else {
        if (c == 2) {
          MatSetType(ptA->_petsc, MATSHELL);
          User< LinearSolver< Dmat > > user = nullptr;
          PetscNew(&user);
          user->op = nullptr;
          const Polymorphic* op = nargs[0] ? dynamic_cast< const Polymorphic* >(nargs[0]) : nullptr;
          if(op) {
              const OneOperator* codeAt = op->Find("(", ArrayOfaType(atype< KN< PetscScalar >* >( ), false));
              if (codeAt) {
                  user->op = new LinearSolver< Dmat >::MatF_O(ptC->_last - ptC->_first, stack, codeA, ptB->_last - ptB->_first, codeAt);
                  MatShellSetOperation(ptA->_petsc, MATOP_MULT_TRANSPOSE,
                          (void (*)(void))Op_User< LinearSolver< Mat >, Mat, 'T' >);
              }
          }
          if(!user->op) user->op = new LinearSolver< Dmat >::MatF_O(ptC->_last - ptC->_first, stack, codeA, ptB->_last - ptB->_first);
          MatShellSetContext(ptA->_petsc, user);
          MatShellSetOperation(ptA->_petsc, MATOP_MULT,
                               (void (*)(void))Op_User< LinearSolver< Mat >, Mat >);
          MatShellSetOperation(ptA->_petsc, MATOP_DESTROY, (void (*)(void))ShellDestroy< LinearSolver< Dmat >  >);
        }
        MatSetUp(ptA->_petsc);
        MatSetOption(ptA->_petsc, MAT_NO_OFF_PROC_ENTRIES, PETSC_TRUE);
        if(ptB->_A->getMatrix() && ptC->_A->getMatrix()) {
          ptA->_num = new PetscInt[ptB->_A->getMatrix( )->HPDDM_n + ptC->_A->getMatrix( )->HPDDM_m];
          ptA->_cnum = ptA->_num + ptB->_A->getMatrix( )->HPDDM_n;
          std::copy_n(ptB->_num, ptB->_A->getMatrix( )->HPDDM_n, ptA->_num);
          std::copy_n(ptC->_num, ptC->_A->getMatrix( )->HPDDM_m, ptA->_cnum);
        } else {
          ptA->_num = new PetscInt[ptB->_A->getDof() + ptC->_A->getDof()];
          ptA->_cnum = ptA->_num + ptB->_A->getDof();
          std::copy_n(ptB->_num, ptB->_A->getDof(), ptA->_num);
          std::copy_n(ptC->_num, ptC->_A->getDof(), ptA->_cnum);
        }
      }
      ptA->_exchange = new HPDDM::template Subdomain< PetscScalar >*[2];
      ptA->_exchange[0] = new HPDDM::template Subdomain< PetscScalar >(*ptB->_A);
      ptA->_exchange[0]->setBuffer( );
      ptA->_exchange[1] = new HPDDM::template Subdomain< PetscScalar >(*ptC->_A);
      ptA->_exchange[1]->setBuffer( );
    }
    if (c == 0 && nargs[0] && GetAny< bool >((*nargs[0])(stack))) ptK->destroy( );
    if (c != 3)
        return ptA;
    else
        return 0L;
  }
  template< class Type >
  class NonlinearSolver : public OneOperator {
   public:
    typedef KN< PetscScalar > Kn;
    typedef KN_< PetscScalar > Kn_;
    class VecF_O {
     public:
      Stack stack;
      mutable Kn x;
      C_F0 c_x;
      mutable Kn x_e;
      C_F0 c_x_e;
      mutable Kn x_i;
      C_F0 c_x_i;
      Expression mat;
      VecF_O(int n, Stack stk, const OneOperator* op)
        : stack(stk), x(n), c_x(CPValue(x)),
          mat(op ? CastTo< long >(C_F0(op->code(basicAC_F0_wa(c_x)), (aType)*op)) : 0) {}
      VecF_O(int n, Stack stk, const OneOperator* op, int)
        : stack(stk), x(n), c_x(CPValue(x)),
          mat(op ? CastTo< PetscReal >(C_F0(op->code(basicAC_F0_wa(c_x)), (aType)*op)) : 0) {}
      VecF_O(int n, Stack stk, const OneOperator* op, int, int ni, int ne)
        : stack(stk), x(n), c_x(CPValue(x)), x_e(std::max(0, ne)), c_x_e( ), x_i(std::max(0, ni)),
          c_x_i( ), mat( ) {
        if (ne >= 0) c_x_e = CPValue(x_e);
        if (ni >= 0) c_x_i = CPValue(x_i);
        if (op) {
          if (ni >= 0 && ne >= 0)
            mat = CastTo< long >(C_F0(op->code(basicAC_F0_wa({c_x, c_x_e, c_x_i})), (aType)*op));
          else if (ni >= 0)
            mat = CastTo< long >(C_F0(op->code(basicAC_F0_wa({c_x, c_x_i})), (aType)*op));
          else
            mat = CastTo< long >(C_F0(op->code(basicAC_F0_wa({c_x, c_x_e})), (aType)*op));
        }
      }
      ~VecF_O( ) {
        delete mat;
        Expression zzz = c_x;
        delete zzz;
        zzz = c_x_e;
        delete zzz;
        zzz = c_x_i;
        delete zzz;
      }
      long apply(const Kn_& xx) const {
        x = xx;
        ffassert(mat);
        long ret = GetAny< long >((*mat)(stack));
        WhereStackOfPtr2Free(stack)->clean( );
        return ret;
      }
      long apply(const Kn_& xx, const Kn_& xx_e, const Kn_& xx_i) const {
        x = xx;
        x_e = xx_e;
        x_i = xx_i;
        ffassert(mat);
        long ret = GetAny< long >((*mat)(stack));
        WhereStackOfPtr2Free(stack)->clean( );
        return ret;
      }
      long apply(const Kn_& xx, const Kn_& xx_d) const {
        x = xx;
        Expression zzz = c_x_e;
        if (zzz)
          x_e = xx_d;
        else
          x_i = xx_d;
        ffassert(mat);
        long ret = GetAny< long >((*mat)(stack));
        WhereStackOfPtr2Free(stack)->clean( );
        return ret;
      }
      long apply(const Kn_& xx, PetscReal* f) const {
        x = xx;
        ffassert(mat);
        double ret = GetAny< double >((*mat)(stack));
        *f = ret;
        WhereStackOfPtr2Free(stack)->clean( );
        return 0;
      }
    };
    class IVecF_O {
     public:
      Stack stack;
      mutable double t;
      C_F0 c_t;
      mutable Kn x;
      C_F0 c_x;
      mutable Kn x_t;
      C_F0 c_x_t;
      mutable double a;
      C_F0 c_a;
      Expression mat;
      IVecF_O(int n, Stack stk, const OneOperator* op)
        : stack(stk), t(0), c_t(CPValue(t)), x(n), c_x(CPValue(x)), x_t(n), c_x_t(CPValue(x_t)),
          a(0), c_a(CPValue(a)),
          mat(op ? CastTo< long >(C_F0(op->code(basicAC_F0_wa({c_t, c_x, c_x_t, c_a})), (aType)*op))
                 : 0) {}
      ~IVecF_O( ) {
        delete mat;
        Expression zzz = c_t;
        delete zzz;
        zzz = c_x;
        delete zzz;
        zzz = c_x_t;
        delete zzz;
        zzz = c_a;
        delete zzz;
      }
      long apply(const double& tt, const Kn_& xx, const Kn_& xx_t, const double& aa) const {
        t = tt;
        x = xx;
        x_t = xx_t;
        a = aa;
        long ret = mat ? GetAny< long >((*mat)(stack)) : 1;
        WhereStackOfPtr2Free(stack)->clean( );
        return ret;
      }
    };
    class IMatF_O {
     public:
      Stack stack;
      mutable double t;
      C_F0 c_t;
      mutable Kn x;
      C_F0 c_x;
      mutable Kn x_t;
      C_F0 c_x_t;
      Expression mat;
      IMatF_O(int n, Stack stk, const OneOperator* op)
        : stack(stk), t(0), c_t(CPValue(t)), x(n), c_x(CPValue(x)), x_t(n), c_x_t(CPValue(x_t)),
          mat(op ? CastTo< Kn_ >(C_F0(op->code(basicAC_F0_wa({c_t, c_x, c_x_t})), (aType)*op))
                 : 0) {}
      IMatF_O(int n, Stack stk, const OneOperator* op, int)
        : stack(stk), t(0), c_t(CPValue(t)), x(n), c_x(CPValue(x)),
          mat(op ? CastTo< Kn_ >(C_F0(op->code(basicAC_F0_wa({c_t, c_x})), (aType)*op)) : 0) {}
      ~IMatF_O( ) {
        delete mat;
        Expression zzz = c_t;
        delete zzz;
        zzz = c_x;
        delete zzz;
        zzz = c_x_t;
        delete zzz;
      }
      void apply(const double& tt, const Kn_& xx, const Kn_& xx_t, Kn_& res) const {
        t = tt;
        x = xx;
        x_t = xx_t;
        ffassert(mat);
        res = GetAny< Kn_ >((*mat)(stack));
        WhereStackOfPtr2Free(stack)->clean( );
      }
      void apply(const double& tt, const Kn_& xx, Kn_& res) const {
        t = tt;
        x = xx;
        ffassert(mat);
        res = GetAny< Kn_ >((*mat)(stack));
        WhereStackOfPtr2Free(stack)->clean( );
      }
    };
    class IMonF_O {
     public:
      Stack stack;
      mutable long s;
      C_F0 c_s;
      mutable double t;
      C_F0 c_t;
      mutable Kn x;
      C_F0 c_x;
      Expression mat;
      IMonF_O(int n, Stack stk, const OneOperator* op)
        : stack(stk), s(0), c_s(CPValue(s)), t(0), c_t(CPValue(t)), x(n), c_x(CPValue(x)),
          mat(op ? CastTo< long >(C_F0(op->code(basicAC_F0_wa({c_s, c_t, c_x})), (aType)*op)) : 0) {
      }
      ~IMonF_O( ) {
        delete mat;
        Expression zzz = c_s;
        delete zzz;
        zzz = c_t;
        delete zzz;
        zzz = c_x;
        delete zzz;
      }
      long apply(const long& ss, const double& tt, const Kn_& xx) const {
        s = ss;
        t = tt;
        x = xx;
        ffassert(mat);
        long ret = GetAny< long >((*mat)(stack));
        WhereStackOfPtr2Free(stack)->clean( );
        return ret;
      }
    };
    class IConvF_O {
     public:
      Stack stack;
      mutable long it;
      C_F0 c_it;
      mutable double xnorm;
      C_F0 c_xnorm;
      mutable double gnorm;
      C_F0 c_gnorm;
      mutable double f;
      C_F0 c_f;
      mutable Kn x;
      C_F0 c_x;
      mutable Kn dx;
      C_F0 c_dx;
      Expression mat;
      IConvF_O(int n, Stack stk, const OneOperator* op)
        : stack(stk), it(0), c_it(CPValue(it)), xnorm(0), c_xnorm(CPValue(xnorm)), gnorm(0),
          c_gnorm(CPValue(gnorm)), f(0), c_f(CPValue(f)), x(n), c_x(CPValue(x)), dx(n),
          c_dx(CPValue(dx)),
          mat(op ? CastTo< long >(C_F0(
                     op->code(basicAC_F0_wa({c_it, c_xnorm, c_gnorm, c_f, c_x, c_dx})), (aType)*op))
                 : 0) {}
      ~IConvF_O( ) {
        delete mat;
        Expression zzz = c_it;
        delete zzz;
        zzz = c_xnorm;
        delete zzz;
        zzz = c_gnorm;
        delete zzz;
        zzz = c_f;
        delete zzz;
        zzz = c_x;
        delete zzz;
        zzz = c_dx;
        delete zzz;
      }
      long apply(const long& iit, const double& ixnorm, const double& ignorm, const double& iff,
                 const Kn_& ix, const Kn_& idx) const {
        it = iit;
        xnorm = ixnorm;
        gnorm = ignorm;
        f = iff;
        x = ix;
        dx = idx;
        ffassert(mat);
        long ret = GetAny< long >((*mat)(stack));
        WhereStackOfPtr2Free(stack)->clean( );
        return ret;
      }
    };
    const int c;
    class E_NonlinearSolver : public E_F0mps {
     public:
      Expression A;
      Expression J;
      Expression r;
      Expression b;
      Expression x;
      const OneOperator *codeJ, *codeR, *codeRHS;
      const int c;
      static const int n_name_param = 12;
      static basicAC_F0::name_and_type name_param[];
      Expression nargs[n_name_param];
      E_NonlinearSolver(const basicAC_F0& args, int d)
        : A(0), J(0), r(0), x(0), codeJ(0), codeR(0), codeRHS(0), c(d) {
        args.SetNameParam(n_name_param, name_param, nargs);
        A = to< Type* >(args[0]);
        const Polymorphic* op = dynamic_cast< const Polymorphic* >(args[1].LeftValue( ));
        ffassert(op);
        if (c == 0 || c == 1 || c == 4) {
          codeJ = op->Find("(", ArrayOfaType(atype< KN< PetscScalar >* >( ), false));
          op = dynamic_cast< const Polymorphic* >(args[2].LeftValue( ));
          ffassert(op);
          codeR = op->Find("(", ArrayOfaType(atype< KN< PetscScalar >* >( ), false));
        } else {
          codeJ =
            op->Find("(", ArrayOfaType(atype< double >( ), atype< KN< PetscScalar >* >( ),
                                       atype< KN< PetscScalar >* >( ), atype< double >( ), false));
          op = dynamic_cast< const Polymorphic* >(args[2].LeftValue( ));
          ffassert(op);
          codeR = op->Find("(", ArrayOfaType(atype< double >( ), atype< KN< PetscScalar >* >( ),
                                             atype< KN< PetscScalar >* >( ), false));
          if (!codeR)
            codeRHS = op->Find(
              "(", ArrayOfaType(atype< double >( ), atype< KN< PetscScalar >* >( ), false));
        }
        if (c == 0 || c == 2 || c == 4)
          x = to< KN< PetscScalar >* >(args[3]);
        else {
          if (c == 1)
            b = to< KN< PetscScalar >* >(args[3]);
          else {
            op = dynamic_cast< const Polymorphic* >(args[3].LeftValue( ));
            ffassert(op);
            codeRHS = op->Find(
              "(", ArrayOfaType(atype< double >( ), atype< KN< PetscScalar >* >( ), false));
          }
          x = to< KN< PetscScalar >* >(args[4]);
        }
      }

      AnyType operator( )(Stack stack) const;
      operator aType( ) const { return atype< long >( ); }
    };
    E_F0* code(const basicAC_F0& args) const { return new E_NonlinearSolver(args, c); }
    NonlinearSolver(int I)
      : OneOperator(atype< long >( ), atype< Type* >( ), atype< Polymorphic* >( ),
                    atype< Polymorphic* >( ), atype< KN< PetscScalar >* >( )),
        c(I == 1 ? 0 : 4) {}
    NonlinearSolver( )
      : OneOperator(atype< long >( ), atype< Type* >( ), atype< Polymorphic* >( ),
                    atype< Polymorphic* >( ), atype< KN< PetscScalar >* >( ),
                    atype< KN< PetscScalar >* >( )),
        c(1) {}
    NonlinearSolver(int, int)
      : OneOperator(atype< long >( ), atype< Type* >( ), atype< Polymorphic* >( ),
                    atype< Polymorphic* >( ), atype< KN< PetscScalar >* >( )),
        c(2) {}
    NonlinearSolver(int, int, int)
      : OneOperator(atype< long >( ), atype< Type* >( ), atype< Polymorphic* >( ),
                    atype< Polymorphic* >( ), atype< Polymorphic* >( ),
                    atype< KN< PetscScalar >* >( )),
        c(3) {}
  };
  template< class Type >
  basicAC_F0::name_and_type NonlinearSolver< Type >::E_NonlinearSolver::name_param[] = {
    {"sparams", &typeid(std::string*)},
    {"xl", &typeid(KN< PetscScalar >*)},
    {"xu", &typeid(KN< PetscScalar >*)},
    {"monitor", &typeid(Polymorphic*)},
    {"HessianRoutine", &typeid(Polymorphic*)},
    {"InequalityConstraints", &typeid(Polymorphic*)},
    {"EqualityConstraints", &typeid(Polymorphic*)},
    {"JacobianInequality", &typeid(Polymorphic*)},
    {"JacobianEquality", &typeid(Polymorphic*)},
    {"JI", &typeid(Type*)},
    {"JE", &typeid(Type*)},
    {"reason", &typeid(long*)}};
  template< class Type >
  PetscErrorCode FormJacobian(SNES snes, Vec x, Mat J, Mat B, void* ctx) {
    User< Type >* user;
    const PetscScalar* in;

    PetscFunctionBeginUser;
    user = reinterpret_cast< User< Type >* >(ctx);
    typename NonlinearSolver< Type >::VecF_O* mat =
      reinterpret_cast< typename NonlinearSolver< Type >::VecF_O* >((*user)->op);
    VecGetArrayRead(x, &in);
    KN_< PetscScalar > xx(const_cast< PetscScalar* >(in), mat->x.n);
    long ret = mat->apply(xx);
    VecRestoreArrayRead(x, &in);
    PetscFunctionReturn(PetscErrorCode(ret));
  }
  template< class Type >
  PetscErrorCode Convergence(SNES snes, PetscInt it, PetscReal xnorm, PetscReal gnorm, PetscReal f,
                             SNESConvergedReason* reason, void* ctx) {
    User< Type >* user;
    const PetscScalar* in;
    const PetscScalar* din;
    Vec u, du;

    PetscFunctionBeginUser;
    user = reinterpret_cast< User< Type >* >(ctx);
    typename NonlinearSolver< Type >::IConvF_O* mat =
      reinterpret_cast< typename NonlinearSolver< Type >::IConvF_O* >((*user)->conv);
    SNESGetSolution(snes, &u);
    SNESGetSolutionUpdate(snes, &du);
    VecGetArrayRead(u, &in);
    VecGetArrayRead(du, &din);
    KN_< PetscScalar > xx(const_cast< PetscScalar* >(in), mat->x.n);
    KN_< PetscScalar > ww(const_cast< PetscScalar* >(din), mat->x.n);
    SNESConvergedDefault(snes, it, xnorm, gnorm, f, reason, ctx);
    if (*reason == SNES_CONVERGED_ITERATING)
      *reason = SNESConvergedReason(mat->apply(it, xnorm, gnorm, f, xx, ww));
    VecRestoreArrayRead(du, &din);
    VecRestoreArrayRead(u, &in);
    PetscFunctionReturn(PETSC_SUCCESS);
  }
  template< class Type, int I >
  PetscErrorCode FormJacobianTao(Tao tao, Vec x, Mat J, Mat B, void* ctx) {
    User< Type >* user;
    const PetscScalar* in;

    PetscFunctionBeginUser;
    user = reinterpret_cast< User< Type >* >(ctx);
    typename NonlinearSolver< Type >::VecF_O* mat =
      reinterpret_cast< typename NonlinearSolver< Type >::VecF_O* >(
        I == 0 ? (*user)->icJ : (I == 1 ? (*user)->ecJ : (*user)->op));
    VecGetArrayRead(x, &in);
    KN_< PetscScalar > xx(const_cast< PetscScalar* >(in), mat->x.n);
    long ret;
    MatSetOption(J, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    if (I < 2 || (!(*user)->icJ && !(*user)->ecJ)) {
      ret = mat->apply(xx);
    } else {
      Vec DE = nullptr, DI = nullptr;
      const PetscScalar* in_e;
      const PetscScalar* in_i;
      TaoGetDualVariables(tao, (*user)->ecJ ? &DE : nullptr, (*user)->icJ ? &DI : nullptr);
      if (DE) VecGetArrayRead(DE, &in_e);
      if (DI) VecGetArrayRead(DI, &in_i);
      if (DE && DI) {
        KN_< PetscScalar > xx_e(const_cast< PetscScalar* >(in_e), mat->x_e.n);
        KN_< PetscScalar > xx_i(const_cast< PetscScalar* >(in_i), mat->x_i.n);
        ret = mat->apply(xx, xx_e, xx_i);
      } else if (DE) {
        KN_< PetscScalar > xx_e(const_cast< PetscScalar* >(in_e), mat->x_e.n);
        ret = mat->apply(xx, xx_e);
      } else {
        KN_< PetscScalar > xx_i(const_cast< PetscScalar* >(in_i), mat->x_i.n);
        ret = mat->apply(xx, xx_i);
      }
      if (DI) VecRestoreArrayRead(DI, &in_i);
      if (DE) VecRestoreArrayRead(DE, &in_e);
    }
    VecRestoreArrayRead(x, &in);
    PetscFunctionReturn(PetscErrorCode(ret));
  }
  template< class Type, int I >
  PetscErrorCode FormConstraintsTao(Tao obj, Vec x, Vec c, void* ctx) {
    User< Type >* user;
    const PetscScalar* in;
    PetscScalar* out;

    PetscFunctionBeginUser;
    user = reinterpret_cast< User< Type >* >(ctx);
    typename LinearSolver< Type >::MatF_O* mat =
      reinterpret_cast< typename LinearSolver< Type >::MatF_O* >(I == 0 ? (*user)->ic
                                                                             : (*user)->ec);
    VecGetArrayRead(x, &in);
    VecGetArray(c, &out);
    KN_< PetscScalar > xx(const_cast< PetscScalar* >(in), mat->N);
    KN_< PetscScalar > yy(out, mat->M);
    yy = *mat * xx;
    VecRestoreArray(c, &out);
    VecRestoreArrayRead(x, &in);
    PetscFunctionReturn(PETSC_SUCCESS);
  }
  template< class PType, class Type >
  PetscErrorCode FormFunction(PType obj, Vec x, Vec f, void* ctx) {
    User< Type >* user;
    const PetscScalar* in;
    PetscScalar* out;

    PetscFunctionBeginUser;
    user = reinterpret_cast< User< Type >* >(ctx);
    typename LinearSolver< Type >::MatF_O* mat =
      reinterpret_cast< typename LinearSolver< Type >::MatF_O* >((*user)->r);
    VecGetArrayRead(x, &in);
    VecGetArray(f, &out);
    KN_< PetscScalar > xx(const_cast< PetscScalar* >(in), mat->N);
    KN_< PetscScalar > yy(out, mat->N);
    yy = *mat * xx;
    VecRestoreArray(f, &out);
    VecRestoreArrayRead(x, &in);
    PetscFunctionReturn(PETSC_SUCCESS);
  }
  template< class Type >
  PetscErrorCode FormObjectiveRoutine(Tao tao, Vec x, PetscReal* f, void* ctx) {
    User< Type >* user;
    const PetscScalar* in;

    PetscFunctionBeginUser;
    user = reinterpret_cast< User< Type >* >(ctx);
    typename NonlinearSolver< Type >::VecF_O* J =
      reinterpret_cast< typename NonlinearSolver< Type >::VecF_O* >((*user)->J);
    VecGetArrayRead(x, &in);
    KN_< PetscScalar > xx(const_cast< PetscScalar* >(in), J->x.n);
    J->apply(xx, f);
    VecRestoreArrayRead(x, &in);
    PetscFunctionReturn(PETSC_SUCCESS);
  }
  template< class Type >
  PetscErrorCode FormIJacobian(TS ts, PetscReal t, Vec u, Vec u_t, PetscReal a, Mat J, Mat B,
                               void* ctx) {
    User< Type >* user;
    const PetscScalar *in, *in_t;

    PetscFunctionBeginUser;
    user = reinterpret_cast< User< Type >* >(ctx);
    typename NonlinearSolver< Type >::IVecF_O* mat =
      reinterpret_cast< typename NonlinearSolver< Type >::IVecF_O* >((*user)->op);
    VecGetArrayRead(u, &in);
    VecGetArrayRead(u_t, &in_t);
    KN_< PetscScalar > xx(const_cast< PetscScalar* >(in), mat->x.n);
    KN_< PetscScalar > xx_t(const_cast< PetscScalar* >(in_t), mat->x.n);
    long ret = mat->apply(t, xx, xx_t, a);
    VecRestoreArrayRead(u_t, &in);
    VecRestoreArrayRead(u, &in);
    PetscFunctionReturn(PetscErrorCode(ret));
  }
  template< class Type >
  PetscErrorCode FormIFunction(TS ts, PetscReal t, Vec u, Vec u_t, Vec F, void* ctx) {
    User< Type >* user;
    const PetscScalar *in, *in_t;
    PetscScalar* out;

    PetscFunctionBeginUser;
    user = reinterpret_cast< User< Type >* >(ctx);
    typename NonlinearSolver< Type >::IMatF_O* mat =
      reinterpret_cast< typename NonlinearSolver< Type >::IMatF_O* >((*user)->r);
    VecGetArrayRead(u, &in);
    VecGetArrayRead(u_t, &in_t);
    VecGetArray(F, &out);
    KN_< PetscScalar > xx(const_cast< PetscScalar* >(in), mat->x.n);
    KN_< PetscScalar > xx_t(const_cast< PetscScalar* >(in_t), mat->x.n);
    KN_< PetscScalar > yy(out, mat->x.n);
    mat->apply(t, xx, xx_t, yy);
    VecRestoreArray(F, &out);
    VecRestoreArrayRead(u_t, &in);
    VecRestoreArrayRead(u, &in);
    PetscFunctionReturn(PETSC_SUCCESS);
  }
  template< class Type >
  PetscErrorCode FormRHSFunction(TS ts, PetscReal t, Vec u, Vec F, void* ctx) {
    User< Type >* user;
    const PetscScalar* in;
    PetscScalar* out;

    PetscFunctionBeginUser;
    user = reinterpret_cast< User< Type >* >(ctx);
    typename NonlinearSolver< Type >::IMatF_O* mat =
      reinterpret_cast< typename NonlinearSolver< Type >::IMatF_O* >((*user)->rhs);
    VecGetArrayRead(u, &in);
    VecGetArray(F, &out);
    KN_< PetscScalar > xx(const_cast< PetscScalar* >(in), mat->x.n);
    KN_< PetscScalar > yy(out, mat->x.n);
    mat->apply(t, xx, yy);
    VecRestoreArray(F, &out);
    VecRestoreArrayRead(u, &in);
    PetscFunctionReturn(PETSC_SUCCESS);
  }
  template< class Type >
  PetscErrorCode Monitor(TS ts, PetscInt step, PetscReal time, Vec u, void* ctx) {
    User< Type >* user;
    const PetscScalar* in;

    PetscFunctionBeginUser;
    user = reinterpret_cast< User< Type >* >(ctx);
    typename NonlinearSolver< Type >::IMonF_O* mat =
      reinterpret_cast< typename NonlinearSolver< Type >::IMonF_O* >((*user)->mon);
    VecGetArrayRead(u, &in);
    KN_< PetscScalar > xx(const_cast< PetscScalar* >(in), mat->x.n);
    mat->apply(step, time, xx);
    VecRestoreArrayRead(u, &in);
    PetscFunctionReturn(PETSC_SUCCESS);
  }
  template< class Type >
  AnyType NonlinearSolver< Type >::E_NonlinearSolver::operator( )(Stack stack) const {
    Type* ptA = GetAny< Type* >((*A)(stack));
    if (ptA->_petsc) {
      KN< PetscScalar >* in = GetAny< KN< PetscScalar >* >((*x)(stack));
      PetscInt first = ptA->_first;
      PetscInt last = ptA->_last;
      if (!ptA->_num) MatGetOwnershipRange(ptA->_petsc, &first, &last);
      ffassert(in->n == 0 || in->n == last - first);
      if (in->n == 0) in->resize(last - first);
      if(nargs[0]) {
        std::string* options = GetAny< std::string* >((*nargs[0])(stack));
        PetscOptionsInsertString(NULL, options->c_str());
      }
      Vec r, x;
      MatCreateVecs(ptA->_petsc, &r, NULL);
      {
        PetscInt n;
        VecGetSize(r, &n);
        VecCreateMPIWithArray(PetscObjectComm((PetscObject)ptA->_petsc), 1, in->n, n, static_cast< PetscScalar* >(*in), &x);
      }
      KSP ksp;
      if (c < 2 || c == 4) {
        KN< PetscScalar >* fl =
          nargs[1] ? GetAny< KN< PetscScalar >* >((*nargs[1])(stack)) : nullptr;
        bool setXl = (fl && fl->n == in->n);
        KN< PetscScalar >* fu =
          nargs[2] ? GetAny< KN< PetscScalar >* >((*nargs[2])(stack)) : nullptr;
        bool setXu = (fu && fu->n == in->n);
        Vec xu, xl = nullptr;
        if (setXl || setXu) {
          MatCreateVecs(ptA->_petsc, &xu, &xl);
          PetscScalar* ptr;
          if (setXl) {
            VecGetArray(xl, &ptr);
            std::copy_n(static_cast< PetscScalar* >(*fl), in->n, ptr);
            VecRestoreArray(xl, &ptr);
          } else
            VecSet(xl, PETSC_NINFINITY);
          if (setXu) {
            VecGetArray(xu, &ptr);
            std::copy_n(static_cast< PetscScalar* >(*fu), in->n, ptr);
            VecRestoreArray(xu, &ptr);
          } else
            VecSet(xu, PETSC_INFINITY);
        }
        if (c == 4) {
          User< TaoSolver< Type > > user = nullptr;
          PetscNew(&user);
          user->J = new NonlinearSolver< Type >::VecF_O(in->n, stack, codeJ, 1);
          user->r = new typename LinearSolver< Type >::MatF_O(in->n, stack, codeR);
          Tao tao;
          TaoCreate(PETSC_COMM_WORLD, &tao);
          TaoSetObjective(tao, FormObjectiveRoutine< TaoSolver< Type > >, &user);
          TaoSetGradient(tao, NULL, FormFunction< Tao, TaoSolver< Type > >, &user);
          if (setXl || setXu) TaoSetVariableBounds(tao, xl, xu);
          TaoSetSolution(tao, x);
          PetscInt sizeE, sizeI;
          Vec ce = nullptr, ci = nullptr;
          const Polymorphic* op = nargs[5] ? dynamic_cast< const Polymorphic* >(nargs[5]) : nullptr;
          Type* pt = nargs[9] ? GetAny< Type* >((*nargs[9])(stack)) : nullptr;
          if (pt && pt->_petsc) {
            MatGetLocalSize(pt->_petsc, &sizeI, NULL);
            if (op) {
              const OneOperator* code = op->Find("(", ArrayOfaType(atype< Kn* >( ), false));
              user->ic = new typename LinearSolver< Type >::MatF_O(in->n, stack, code, sizeI);
              VecCreateMPI(PETSC_COMM_WORLD, sizeI, PETSC_DECIDE, &ci);
              TaoSetInequalityConstraintsRoutine(tao, ci,
                                                 FormConstraintsTao< TaoSolver< Type >, 0 >, &user);
            }
            if (user->ic) {
              op = nargs[7] ? dynamic_cast< const Polymorphic* >(nargs[7]) : nullptr;
              if (op && pt && pt->_petsc) {
                const OneOperator* code = op->Find("(", ArrayOfaType(atype< Kn* >( ), false));
                user->icJ = new NonlinearSolver< Type >::VecF_O(in->n, stack, code);
                TaoSetJacobianInequalityRoutine(tao, pt->_petsc, pt->_petsc,
                                                FormJacobianTao< TaoSolver< Type >, 0 >, &user);
              }
            }
          }
          pt = nargs[10] ? GetAny< Type* >((*nargs[10])(stack)) : nullptr;
          if (pt && pt->_petsc) {
            MatGetLocalSize(pt->_petsc, &sizeE, NULL);
            op = nargs[6] ? dynamic_cast< const Polymorphic* >(nargs[6]) : nullptr;
            if (op) {
              const OneOperator* code = op->Find("(", ArrayOfaType(atype< Kn* >( ), false));
              user->ec = new typename LinearSolver< Type >::MatF_O(in->n, stack, code, sizeE);
              VecCreateMPI(PETSC_COMM_WORLD, sizeE, PETSC_DECIDE, &ce);
              TaoSetEqualityConstraintsRoutine(tao, ce, FormConstraintsTao< TaoSolver< Type >, 1 >,
                                               &user);
            }
            if (user->ec) {
              op = nargs[8] ? dynamic_cast< const Polymorphic* >(nargs[8]) : nullptr;
              if (op) {
                const OneOperator* code = op->Find("(", ArrayOfaType(atype< Kn* >( ), false));
                user->ecJ = new NonlinearSolver< Type >::VecF_O(in->n, stack, code);
                TaoSetJacobianEqualityRoutine(tao, pt->_petsc, pt->_petsc,
                                              FormJacobianTao< TaoSolver< Type >, 1 >, &user);
              }
            }
          }
          op = nargs[4] ? dynamic_cast< const Polymorphic* >(nargs[4]) : nullptr;
          if (op) {
            if (user->icJ && user->ecJ) {
              const OneOperator* codeH = op->Find(
                "(", ArrayOfaType(atype< Kn* >( ), atype< Kn* >( ), atype< Kn* >( ), false));
              user->op = new NonlinearSolver< Type >::VecF_O(in->n, stack, codeH, 2, sizeI, sizeE);
              TaoSetHessian(tao, ptA->_petsc, ptA->_petsc,
                            FormJacobianTao< TaoSolver< Type >, 5 >, &user);
            } else if (user->icJ) {
              const OneOperator* codeH =
                op->Find("(", ArrayOfaType(atype< Kn* >( ), atype< Kn* >( ), false));
              user->op = new NonlinearSolver< Type >::VecF_O(in->n, stack, codeH, 2, sizeI, -1);
              TaoSetHessian(tao, ptA->_petsc, ptA->_petsc,
                            FormJacobianTao< TaoSolver< Type >, 3 >, &user);
            } else if (user->ecJ) {
              const OneOperator* codeH =
                op->Find("(", ArrayOfaType(atype< Kn* >( ), atype< Kn* >( ), false));
              user->op = new NonlinearSolver< Type >::VecF_O(in->n, stack, codeH, 2, -1, sizeE);
              TaoSetHessian(tao, ptA->_petsc, ptA->_petsc,
                            FormJacobianTao< TaoSolver< Type >, 4 >, &user);
            } else {
              const OneOperator* codeH = op->Find("(", ArrayOfaType(atype< Kn* >( ), false));
              user->op = new NonlinearSolver< Type >::VecF_O(in->n, stack, codeH);
              TaoSetHessian(tao, ptA->_petsc, ptA->_petsc,
                            FormJacobianTao< TaoSolver< Type >, 2 >, &user);
            }
          }
          TaoSetFromOptions(tao);
          TaoSolve(tao);
          if (ce) VecDestroy(&ce);
          if (ci) VecDestroy(&ci);
          TaoDestroy(&tao);
          delete user->ecJ;
          delete user->ec;
          delete user->icJ;
          delete user->ic;
          delete user->op;
          delete user->r;
          delete user->J;
          PetscFree(user);
        } else {
          User< NonlinearSolver< Type > > user = nullptr;
          PetscNew(&user);
          user->op = new NonlinearSolver< Type >::VecF_O(in->n, stack, codeJ);
          user->r = new typename LinearSolver< Type >::MatF_O(in->n, stack, codeR);
          user->conv = nullptr;
          SNES snes;
          SNESCreate(PETSC_COMM_WORLD, &snes);
          Vec f;
          SNESSetFunction(snes, r, FormFunction< SNES, NonlinearSolver< Type > >, &user);
          SNESSetJacobian(snes, ptA->_petsc, ptA->_petsc, FormJacobian< NonlinearSolver< Type > >,
                          &user);
          if (setXl || setXu) SNESVISetVariableBounds(snes, xl, xu);
          SNESSetFromOptions(snes);
          if (ptA->_ksp) {
            SNESGetKSP(snes, &ksp);
            PetscObjectReference((PetscObject)ksp);
            SNESSetKSP(snes, ptA->_ksp);
          }
          if (c == 1) {
            PetscInt n;
            VecGetSize(r, &n);
            KN< PetscScalar >* fb = GetAny< KN< PetscScalar >* >((*b)(stack));
            ffassert(fb->n == last - first);
            VecCreateMPIWithArray(PETSC_COMM_WORLD, 1, fb->n, n, static_cast< PetscScalar* >(*fb),
                                  &f);
          }
          const Polymorphic* op = nargs[3] ? dynamic_cast< const Polymorphic* >(nargs[3]) : nullptr;
          if (op) {
            ffassert(op);
            const OneOperator* codeM =
              op->Find("(", ArrayOfaType(atype< long >( ), atype< double >( ), atype< double >( ),
                                         atype< double >( ), atype< KN< PetscScalar >* >( ),
                                         atype< KN< PetscScalar >* >( ), false));
            user->conv = new NonlinearSolver< Type >::IConvF_O(in->n, stack, codeM);
            SNESSetConvergenceTest(snes, Convergence< NonlinearSolver< Type > >, &user, NULL);
          }
          SNESSolve(snes, c == 1 ? f : NULL, x);
          if (c == 1) VecDestroy(&f);
          if (ptA->_ksp) {
            SNESSetKSP(snes, ksp);
            KSPDestroy(&ksp);
          }
          if(nargs[11]) {
            SNESConvergedReason reason;
            SNESGetConvergedReason(snes, &reason);
            long* ret = GetAny< long* >((*nargs[11])(stack));
            *ret = static_cast<long>(reason);
          }
          SNESDestroy(&snes);
          delete user->conv;
          delete user->r;
          delete user->op;
          PetscFree(user);
        }
        if (xl) {
          VecDestroy(&xl);
          VecDestroy(&xu);
        }
      } else {
        User< TimeStepper< Type > > user = nullptr;
        PetscNew(&user);
        user->op = new NonlinearSolver< Type >::IVecF_O(in->n, stack, codeJ);
        user->r = (codeR ? new NonlinearSolver< Type >::IMatF_O(in->n, stack, codeR) : nullptr);
        user->rhs =
          (c == 3 || !codeR ? new NonlinearSolver< Type >::IMatF_O(in->n, stack, codeRHS, 1)
                            : nullptr);
        user->mon = nullptr;
        TS ts;
        TSCreate(PetscObjectComm((PetscObject)ptA->_petsc), &ts);
        TSSetIJacobian(ts, ptA->_petsc, ptA->_petsc, FormIJacobian< TimeStepper< Type > >, &user);
        if (user->r) TSSetIFunction(ts, r, FormIFunction< TimeStepper< Type > >, &user);
        if (user->rhs) TSSetRHSFunction(ts, NULL, FormRHSFunction< TimeStepper< Type > >, &user);
        TSSetFromOptions(ts);
        SNES snes;
        if (ptA->_ksp) {
          TSGetSNES(ts, &snes);
          SNESGetKSP(snes, &ksp);
          PetscObjectReference((PetscObject)ksp);
          SNESSetKSP(snes, ptA->_ksp);
        }
        const Polymorphic* op = nargs[3] ? dynamic_cast< const Polymorphic* >(nargs[3]) : nullptr;
        if (op) {
          ffassert(op);
          const OneOperator* codeM = op->Find(
            "(", ArrayOfaType(atype< long >( ), atype< double >( ), atype< Kn* >( ), false));
          user->mon = new NonlinearSolver< Type >::IMonF_O(in->n, stack, codeM);
          TSMonitorSet(ts, Monitor< TimeStepper< Type > >, &user, NULL);
        }
        TSSetSolution(ts, x);
        TSSolve(ts, x);
        if (ptA->_ksp) {
          SNESSetKSP(snes, ksp);
          KSPDestroy(&ksp);
        }
        TSDestroy(&ts);
        delete user->mon;
        delete user->r;
        delete user->op;
        PetscFree(user);
      }
      VecDestroy(&x);
      VecDestroy(&r);
    }
    return 0L;
  }
  template< class Type, class Container,
            typename std::enable_if< std::is_same< Container, Mat >::value >::type* = nullptr >
  static PetscErrorCode ContainerGetContext(Container A, Type& user) {
    return MatShellGetContext(A, &user);
  }
  template< class Type, class Container,
            typename std::enable_if< std::is_same< Container, PC >::value >::type* = nullptr >
  static PetscErrorCode ContainerGetContext(Container A, Type& user) {
    return PCShellGetContext(A, (void**)&user);
  }
  template< class Type, class Container, char N >
  static PetscErrorCode Op_User(Container A, Vec x, Vec y) {
    static_assert(std::is_same<Container, Mat>::value || N == 'N', "Wrong types");
    User< Type > user;
    const PetscScalar* in;
    PetscScalar* out;
    PetscErrorCode ierr;

    PetscFunctionBeginUser;
    ierr = ContainerGetContext(A, user);CHKERRQ(ierr);
    typename LinearSolver< Type >::MatF_O* mat =
      reinterpret_cast< typename LinearSolver< Type >::MatF_O* >(user->op);
    VecGetArrayRead(x, &in);
    VecGetArray(y, &out);
    if(N == 'N') {
        KN_< PetscScalar > xx(const_cast< PetscScalar* >(in), mat->N);
        KN_< PetscScalar > yy(out, mat->M);
        yy = *mat * xx;
    } else {
        KN_< PetscScalar > xx(const_cast< PetscScalar* >(in), mat->M);
        KN_< PetscScalar > yy(out, mat->N);
        typename LinearSolver< Type >::MatF_O::plusAtx Atxx(mat, xx);
        yy = Atxx;
    }
    VecRestoreArray(y, &out);
    VecRestoreArrayRead(x, &in);
    PetscFunctionReturn(PETSC_SUCCESS);
  }

  template< bool U, char T, class K >
  void loopDistributedVec(Mat nest, HPDDM::Subdomain<PetscScalar>** exchange, KN_<upscaled_type<K>>* const& in, K* out) {
      Mat** mat;
      PetscInt M, N;
      MatNestGetSubMats(nest, &M, &N, &mat);
      Dmat** cast = reinterpret_cast<Dmat**>(exchange);
      PetscScalar* ptr = reinterpret_cast<PetscScalar*>(in->operator upscaled_type<PetscScalar>*());
      MatType type;
      PetscBool isType;
      if(T == 'N') {
          for(PetscInt i = 0; i < M; ++i) {
              for(PetscInt j = 0; j < N; ++j) {
                  if(mat[i][j]) {
                      MatGetType(mat[i][j], &type);
                      PetscObjectTypeCompareAny((PetscObject)mat[i][j], &isType, MATMPIDENSE, MATSEQDENSE, "");
                      PetscInt n = 0, m = 0;
                      if(!isType) {
                          PetscStrcmp(type, MATHERMITIANTRANSPOSEVIRTUAL, &isType);
                          if(isType) {
                              Mat C;
                              MatHermitianTransposeGetMat(mat[i][j], &C);
                              PetscObjectTypeCompareAny((PetscObject)C, &isType, MATMPIDENSE, MATSEQDENSE, "");
                              if(isType) MatGetSize(C, &m, &n);
                              type = MATHERMITIANTRANSPOSEVIRTUAL;
                          }
                      } else MatGetSize(mat[i][j], &n, &m);
                      if(isType && (m > 1 || n == 1)) {
                          if(mpirank == 0) {
                              if(U)
                                  *ptr++ = *out++;
                              else
                                  *out++ = *ptr++;
                          }
                          break;
                      }
                      else if(cast[i * N + j]) {
                          PetscStrcmp(type, MATHERMITIANTRANSPOSEVIRTUAL, &isType);
                          PetscInt* num = nullptr, first, last, n;
                          const HPDDM::Subdomain<PetscScalar>* A = nullptr;
                          if(cast[i * N + j]->_cnum) {
                              if(isType) {
                                  num = cast[i * N + j]->_cnum;
                                  first = cast[i * N + j]->_cfirst;
                                  last = cast[i * N + j]->_clast;
                              }
                              else {
                                  num = cast[i * N + j]->_num;
                                  first = cast[i * N + j]->_first;
                                  last = cast[i * N + j]->_last;
                              }
                              A = cast[i * N + j]->_exchange[isType ? 1 : 0];
                              n = A ? A->getDof() : 0;
                          }
                          else {
                              num = cast[i * N + j]->_num;
                              first = cast[i * N + j]->_first;
                              last = cast[i * N + j]->_last;
                              A = cast[i * N + j]->_A;
                              n = A ? A->getDof() : 0;
                          }
                          if(num) {
                              HPDDM::Subdomain< K >::template distributedVec< U >(num, first, last, ptr, out, n, 1);
                              if (U && A)
                                  A->exchange(ptr);
                              ptr += n;
                              out += last - first;
                              break;
                          }
                      }
                  }
              }
          }
      }
      else {
          ffassert(0);
      }
  }

  template< class T, class U, class K, char trans >
  class InvPETSc {
    static_assert(std::is_same< K, PetscScalar >::value, "Wrong types");

   public:
    const T t;
    const U u;
    InvPETSc(T v, U w) : t(v), u(w) {}
    void solve(U out) const {
      if ((*t)._petsc) {
        Vec x, y;
        double timing = MPI_Wtime( );
        MatCreateVecs((*t)._petsc, &x, &y);
        PetscScalar* ptr;
        MatType type;
        PetscBool isType;
        MatGetType((*t)._petsc, &type);
        PetscStrcmp(type, MATNEST, &isType);
        PetscScalar* p = reinterpret_cast<PetscScalar*>(u->operator upscaled_type<PetscScalar>*());

        if( (*t)._vector_global ){
          for(int i = 0; i < u->n; ++i)
              p[i] = u->operator[](i);
          MPI_Allreduce(MPI_IN_PLACE, p, u->n, HPDDM::Wrapper<K>::mpi_type(), MPI_SUM, PETSC_COMM_WORLD);
          VecGetArray(x, &ptr);
          if(isType) { // case nested matrix
            Mat** mat;
            PetscInt M, N;
            PetscInt Mi = 0, Ni = 0;
            MatNestGetSubMats((*t)._petsc, &M, &N, &mat);

            PetscInt offset_cols=0;
            PetscInt offset_local_cols=0;
            for(PetscInt j = 0; j < N; ++j) { // cols
              for(PetscInt i = 0; i < M; ++i) { // rows
                Ni = 0, Mi = 0;
                if(mat[i][j]) {
                  MatGetSize(mat[i][j], &Mi, &Ni);

                  PetscInt cbegin;
                  PetscInt ni,mi;
                  MatGetLocalSize( mat[i][j], &mi, &ni);
                  MatGetOwnershipRangeColumn( mat[i][j], &cbegin, NULL);

                  for(PetscInt i = 0; i < ni; ++i)
                    ptr[i+offset_local_cols] = p[i+cbegin+offset_cols];

                  offset_cols += Ni;
                  offset_local_cols += ni;
                  break;
                } // if mat[i][j]

              } // loop rows
            } // loop cols
          } // if nested matrix
          else{
            // case not nested matrix
            PetscInt cbegin;
            PetscInt n,m,N,M;
            MatGetLocalSize( (*t)._petsc, &m, &n);
            MatGetOwnershipRangeColumn( (*t)._petsc, &cbegin, NULL);
            for(PetscInt i = 0; i < n; ++i)
              ptr[i]  = p[i+cbegin];
          }
          VecRestoreArray(x, &ptr);
        }
        else if ( std::is_same< typename std::remove_reference< decltype(*t.A->_A) >::type,
                          HpSchwarz< PetscScalar > >::value) {
          VecGetArray(x, &ptr);
          for(int i = 0; i < u->n; ++i)
              p[i] = u->operator[](i);
          if(isType) {
            ffassert((std::is_same<PetscReal, upscaled_type<PetscReal>>::value));
            loopDistributedVec<0, 'N'>((*t)._petsc, (*t)._exchange, u, ptr);
          }
          else
            HPDDM::Subdomain< K >::template distributedVec< 0 >((*t)._num, (*t)._first, (*t)._last,
                                                                p, ptr,
                                                                static_cast<PetscInt>(u->n), 1);
          VecRestoreArray(x, &ptr);
          if ((*t)._ksp) {
            PetscBool nonZero;
            KSPGetInitialGuessNonzero((*t)._ksp, &nonZero);
            if (nonZero) {
              VecGetArray(y, &ptr);
              p = reinterpret_cast<PetscScalar*>(out->operator upscaled_type<PetscScalar>*());
              for(int i = 0; i < out->n; ++i)
                  p[i] = out->operator[](i);
              if(isType)
                loopDistributedVec<0, 'N'>((*t)._petsc, (*t)._exchange, out, ptr);
              else
                HPDDM::Subdomain< K >::template distributedVec< 0 >(
                  (*t)._num, (*t)._first, (*t)._last, p, ptr,
                  static_cast<PetscInt>(out->n), 1);
              VecRestoreArray(y, &ptr);
            }
          }
          if (t.A->_A) std::fill_n(out->operator upscaled_type<PetscScalar>*(), out->n, 0.0);
        } else {
          VecSet(x, PetscScalar( ));
          Vec isVec;
          VecCreateMPIWithArray(PETSC_COMM_SELF, 1, (*t)._A->getMatrix( )->HPDDM_n,
                                (*t)._A->getMatrix( )->HPDDM_n, p, &isVec);
          VecScatterBegin((*t)._scatter, isVec, x, ADD_VALUES, SCATTER_REVERSE);
          VecScatterEnd((*t)._scatter, isVec, x, ADD_VALUES, SCATTER_REVERSE);
          VecDestroy(&isVec);
        }
        timing = MPI_Wtime( );
        if (!(*t)._ksp) {
          KSPCreate(PetscObjectComm((PetscObject)(*t)._petsc), &(*t)._ksp);
          KSPSetOperators((*t)._ksp, (*t)._petsc, (*t)._petsc);
          KSPSetFromOptions((*t)._ksp);
        }
        if (isType) {
          PC pc;
          KSPGetPC((*t)._ksp, &pc);
          PCType type;
          PCGetType(pc, &type);
          PetscStrcmp(type, PCFIELDSPLIT, &isType);
          if (!isType) {
            Mat C;
            prepareConvert((*t)._petsc, &C);
            KSPSetOperators((*t)._ksp, (*t)._petsc, C);
            MatDestroy(&C);
            isType = PETSC_TRUE;
          }
        }
        if (trans == 'N')
          KSPSolve((*t)._ksp, x, y);
        else {
          if (!std::is_same< PetscScalar, PetscReal >::value && t.conjugate) VecConjugate(x);
          KSPSolveTranspose((*t)._ksp, x, y);
          if (!std::is_same< PetscScalar, PetscReal >::value && t.conjugate) VecConjugate(y);
        }
        if (verbosity > 0 && mpirank == 0)
          cout << " --- system solved with PETSc (in " << MPI_Wtime( ) - timing << ")" << endl;
        p = reinterpret_cast<PetscScalar*>(out->operator upscaled_type<PetscScalar>*());
        if((*t)._vector_global){
          VecGetArray(y, &ptr);
          if(isType) {
            Mat** mat;
            PetscInt M, N;
            PetscInt Mi = 0, Ni = 0;
            MatNestGetSubMats((*t)._petsc, &M, &N, &mat);

            PetscInt offset_rows=0;
            PetscInt offset_local_rows=0;
            for(PetscInt i = 0; i < M; ++i) { // rows
              for(PetscInt j = 0; j < N; ++j) { // cols

                Ni = 0, Mi = 0;
                if(mat[i][j]) {
                  MatGetSize(mat[i][j], &Mi, &Ni);

                  PetscInt rbegin;
                  PetscInt ni,mi;
                  MatGetLocalSize( mat[i][j], &mi, &ni);
                  MatGetOwnershipRange( mat[i][j], &rbegin, NULL);

                  for(PetscInt i = 0; i < mi; ++i)
                    p[i+rbegin+offset_rows] = ptr[i+offset_local_rows];

                  offset_rows += Mi;
                  offset_local_rows += mi;
                  break;
                } // if mat[i][j]

              } // loop rows
            } // loop cols
          } // if nested matrix
          else{
            // case not nested matrix
            PetscInt rbegin;
            PetscInt n,m,N,M;
            MatGetOwnershipRange( (*t)._petsc, &rbegin, NULL);
            MatGetLocalSize( (*t)._petsc, &m, &n);
            std::fill(p,p+u->n,0.);
            
            for(PetscInt i = 0; i < m; ++i)
              p[i+rbegin] = ptr[i];
          }

          VecRestoreArray(y, &ptr);
          MPI_Allreduce(MPI_IN_PLACE, p, u->n, HPDDM::Wrapper<K>::mpi_type(), MPI_SUM, PETSC_COMM_WORLD);
        }
        else if (std::is_same< typename std::remove_reference< decltype(*t.A->_A) >::type,
                          HpSchwarz< PetscScalar > >::value) {
          VecGetArray(y, &ptr);
          if(isType)
            loopDistributedVec<1, 'N'>((*t)._petsc, (*t)._exchange, out, ptr);
          else
            HPDDM::Subdomain< K >::template distributedVec< 1 >((*t)._num, (*t)._first, (*t)._last, p,
                                                                ptr, static_cast<PetscInt>(out->n), 1);
          VecRestoreArray(y, &ptr);
        } else {
          Vec isVec;
          VecCreateMPIWithArray(PETSC_COMM_SELF, 1, (*t)._A->getMatrix( )->HPDDM_n,
                                (*t)._A->getMatrix( )->HPDDM_n, p, &isVec);
          VecScatterBegin((*t)._scatter, y, isVec, INSERT_VALUES, SCATTER_FORWARD);
          VecScatterEnd((*t)._scatter, y, isVec, INSERT_VALUES, SCATTER_FORWARD);
          VecDestroy(&isVec);
        }
        VecDestroy(&x);
        VecDestroy(&y);
        if (std::is_same< typename std::remove_reference< decltype(*t.A->_A) >::type,
                          HpSchwarz< PetscScalar > >::value &&
            t.A->_A)
          (*t)._A->exchange(p);
        if(!std::is_same<PetscReal, upscaled_type<PetscReal>>::value) {
          for(int i = out->n - 1; i >= 0; --i)
            out->operator[](i) = p[i];
          p = reinterpret_cast<PetscScalar*>(u->operator upscaled_type<PetscScalar>*());
          for(int i = u->n - 1; i >= 0; --i)
            u->operator[](i) = p[i];
        }
      }
    };
    static U inv(U Ax, InvPETSc< T, U, K, trans > A) {
      A.solve(Ax);
      return Ax;
    }
    static U init(U Ax, InvPETSc< T, U, K, trans > A) {
      PetscInt n, m;
      MatGetSize(A.t.A->_petsc, &n, &m);
      ffassert(n == m);
      Ax->init(A.u->n);
      return inv(Ax, A);
    }
  };

  template< class T, class U, class K, char N >
  class ProdPETSc {
    static_assert(std::is_same< K, PetscScalar >::value, "Wrong types");

   public:
    const T t;
    const U u;
    ProdPETSc(T v, U w) : t(v), u(w) {}
    void prod(U out) const {
      if ((*t)._petsc) {
        PetscBool assembled;
        MatAssembled((*t)._petsc, &assembled);
        if (assembled) {
          Vec x, y;
          PetscScalar* ptr;
          if (N == 'N')
            MatCreateVecs((*t)._petsc, &x, &y);
          else
            MatCreateVecs((*t)._petsc, &y, &x);
          VecGetArray(x, &ptr);
          MatType type;
          MatGetType((*t)._petsc, &type);
          PetscBool isType;
          PetscScalar* p = reinterpret_cast<PetscScalar*>(u->operator upscaled_type<PetscScalar>*());
          PetscStrcmp(type, MATNEST, &isType);
          
          if((*t)._vector_global){
            int nb_comm_size;
            MPI_Comm_size(PETSC_COMM_WORLD,&nb_comm_size);
            for(int i = 0; i < u->n; ++i)
              p[i] = u->operator[](i)/nb_comm_size; // Is there I divided by mpisize
            MPI_Allreduce(MPI_IN_PLACE, p, u->n, HPDDM::Wrapper<K>::mpi_type(), MPI_SUM, PETSC_COMM_WORLD);
            // VecGetArray(x, &ptr) is defined before
            if(isType) { // case nested matrix
              Mat** mat;
              PetscInt Mb, Nb;
              PetscInt Mi = 0, Ni = 0;
              MatNestGetSubMats((*t)._petsc, &Mb, &Nb, &mat);

              PetscInt offset_cols=0;
              PetscInt offset_local_cols=0;
              for(PetscInt j = 0; j < Nb; ++j) { // cols
                for(PetscInt i = 0; i < Mb; ++i) { // rows
                  Ni = 0, Mi = 0;
                  if(mat[i][j]) {
                    MatGetSize(mat[i][j], &Mi, &Ni);

                    PetscInt cbegin;
                    PetscInt ni,mi;
                    MatGetLocalSize( mat[i][j], &mi, &ni);
                    MatGetOwnershipRangeColumn( mat[i][j], &cbegin, NULL);

                    for(PetscInt i = 0; i < ni; ++i)
                      ptr[i+offset_local_cols] = p[i+cbegin+offset_cols];

                    offset_cols += Ni;
                    offset_local_cols += ni;
                    break;
                  } // if mat[i][j]

                } // loop rows
              } // loop cols
            } // if nested matrix
            else{
              // case not nested matrix
              PetscInt cbegin;
              PetscInt n,m;
              MatGetLocalSize( (*t)._petsc, &m, &n);
              MatGetOwnershipRangeColumn( (*t)._petsc, &cbegin, NULL);
              for(PetscInt i = 0; i < n; ++i)
                ptr[i]  = p[i+cbegin];
            }
          }
          else if (isType) {
            ffassert((std::is_same<PetscReal, upscaled_type<PetscReal>>::value));
            loopDistributedVec<0, N>(t->_petsc, t->_exchange, u, ptr);
          } else {
            if(t->_A)
              assert(u->n == t->_A->getDof() && out->n == t->_A->getDof());
            else if (t->_exchange) {
              if (N == 'T') {
                assert(u->n == t->_exchange[0]->getDof() && out->n == t->_exchange[1]->getDof());
              } else {
                assert(u->n == t->_exchange[1]->getDof() && out->n == t->_exchange[0]->getDof());
              }
            }
            for(int i = 0; i < u->n; ++i)
                p[i] = u->operator[](i);
            if (!t->_cnum || N == 'T') {
              if (t->_num)
                HPDDM::Subdomain< K >::template distributedVec< 0 >(
                  t->_num, t->_first, t->_last, p, ptr, static_cast<PetscInt>(u->n), 1);
            } else
              HPDDM::Subdomain< K >::template distributedVec< 0 >(
                t->_cnum, t->_cfirst, t->_clast, p, ptr, static_cast<PetscInt>(u->n), 1);
          }


          VecRestoreArray(x, &ptr);
          if (N == 'T')
            MatMultTranspose(t->_petsc, x, y);
          else
            MatMult(t->_petsc, x, y);
          VecDestroy(&x);
          VecGetArray(y, &ptr);

          if((*t)._vector_global){
            p = reinterpret_cast<PetscScalar*>(out->operator upscaled_type<PetscScalar>*());
            if(isType) {
              Mat** mat;
              PetscInt Mb, Nb;
              PetscInt Mi = 0, Ni = 0;
              MatNestGetSubMats((*t)._petsc, &Mb, &Nb, &mat);

              PetscInt offset_rows=0;
              PetscInt offset_local_rows=0;
              for(PetscInt i = 0; i < Mb; ++i) { // rows
                for(PetscInt j = 0; j < Nb; ++j) { // cols

                  Ni = 0, Mi = 0;
                  if(mat[i][j]) {
                    MatGetSize(mat[i][j], &Mi, &Ni);

                    PetscInt rbegin;
                    PetscInt ni,mi;
                    MatGetLocalSize( mat[i][j], &mi, &ni);
                    MatGetOwnershipRange( mat[i][j], &rbegin, NULL);

                    for(PetscInt i = 0; i < mi; ++i)
                      p[i+rbegin+offset_rows] = ptr[i+offset_local_rows];

                    offset_rows += Mi;
                    offset_local_rows += mi;
                    break;
                  } // if mat[i][j]

                } // loop rows
              } // loop cols
            } // if nested matrix
            else{
              // case not nested matrix
              PetscInt rbegin;
              PetscInt n,m;
              MatGetOwnershipRange( (*t)._petsc, &rbegin, NULL);
              MatGetLocalSize( (*t)._petsc, &m, &n);
              std::fill(p,p+u->n,0.);
              
              for(PetscInt i = 0; i < m; ++i)
                p[i+rbegin] = ptr[i];
            }

            VecRestoreArray(y, &ptr);
            MPI_Allreduce(MPI_IN_PLACE, p, u->n, HPDDM::Wrapper<K>::mpi_type(), MPI_SUM, PETSC_COMM_WORLD);
            for(int i = out->n - 1; i >= 0; --i)
                out->operator[](i) = p[i];
            if(!std::is_same<PetscReal, upscaled_type<PetscReal>>::value) {
              p = reinterpret_cast<PetscScalar*>(u->operator upscaled_type<PetscScalar>*());
              for(int i = u->n - 1; i >= 0; --i)
                u->operator[](i) = p[i];
            }
          }
          else{
          if (!t->_A) std::fill_n(out->operator upscaled_type<PetscScalar>*(), out->n, 0.0);
          if (isType) {
            loopDistributedVec<1, N>(t->_petsc, t->_exchange, out, ptr);
          } else {
            p = reinterpret_cast<PetscScalar*>(out->operator upscaled_type<PetscScalar>*());
            if (!t->_cnum || N == 'N') {
              if (t->_num)
                HPDDM::Subdomain< K >::template distributedVec< 1 >(t->_num, t->_first, t->_last, p,
                                                                    ptr, static_cast<PetscInt>(out->n), 1);
            } else
              HPDDM::Subdomain< K >::template distributedVec< 1 >(t->_cnum, t->_cfirst, t->_clast, p,
                                                                  ptr, static_cast<PetscInt>(out->n), 1);
            if (t->_A)
              (*t)._A->exchange(p);
            else if (t->_exchange) {
              if (N == 'N')
                (*t)._exchange[0]->exchange(p);
              else
                (*t)._exchange[1]->exchange(p);
            }
            for(int i = out->n - 1; i >= 0; --i)
                out->operator[](i) = p[i];
            if(!std::is_same<PetscReal, upscaled_type<PetscReal>>::value) {
                p = reinterpret_cast<PetscScalar*>(u->operator upscaled_type<PetscScalar>*());
                for(int i = u->n - 1; i >= 0; --i)
                    u->operator[](i) = p[i];
            }
          }
          } 

          VecRestoreArray(y, &ptr);
          VecDestroy(&y);
        }
      }
    }
    static U mv(U Ax, ProdPETSc< T, U, K, N > A) {
      *Ax = K( );
      A.prod(Ax);
      return Ax;
    }
    static U init(U Ax, ProdPETSc< T, U, K, N > A) {
      if( A.t->_vector_global ){
        // case rhs and sol are global in freefem
        Ax->init(A.u->n);
      }
      else{
        if (A.t->_A)
          Ax->init(A.u->n);
        else {
          if (A.t->_exchange) {
            if (N == 'T')
              Ax->init(A.t->_exchange[1]->getDof());
            else
              Ax->init(A.t->_exchange[0]->getDof());
          }
        }
      }
      return mv(Ax, A);
    }
  };
  long destroyCSR(Dmat* p) {
    ffassert(p);
    p->dtor();
    return 0L;
  }
  KN_<double> Dmat_D(Dmat* p) {
    throwassert(p && p->_A);
    KN_<double> D(reinterpret_cast<upscaled_type<PetscReal>*>(const_cast<PetscReal*>(p->_A->getScaling())), p->_A->getDof());
    return D;
  }
  long Dmat_n(Dmat* p) {
    throwassert(p);
    PetscInt n = 0;
    if(p->_petsc)
        MatGetLocalSize(p->_petsc, &n, NULL);
    return static_cast<long>(n);
  }
  KN_<long> Dmat_range(Stack stack, Dmat* const& p) {
    throwassert(p);
    PetscInt m, n;
    long* ptr = Add2StackOfPtr2FreeA<long>(stack, new long[2]);
    KN_<long> ret(ptr, 2);
    if(p->_petsc) {
        MatGetOwnershipRange(p->_petsc, &m, &n);
        ptr[0] = m;
        ptr[1] = n;
    }
    return ret;
  }
  template<class K, typename std::enable_if<std::is_same<K, upscaled_type<K>>::value>::type* = nullptr>
  static void init() {
#if !defined(PETSC_USE_REAL_SINGLE)
    if (std::is_same< PetscInt, int >::value) {
      TheOperators->Add("<-", new PETSc::initCSRfromMatrix< HpSchwarz< PetscScalar > >);
    }
#endif
    addScalarProduct< Dmat, PetscScalar >( );

    TheOperators->Add("<-", new PETSc::initCSR< HpSchur< PetscScalar > >);

    Global.Add("ChangeNumbering", "(", new PETSc::changeNumbering< Dmat, KN >( ));
    Global.Add("ChangeNumbering", "(", new PETSc::changeNumbering< Dmat, KN >(1));
    Global.Add("ChangeNumbering", "(", new PETSc::changeNumbering< Dmat, KN >(1, 1));
    Global.Add("ChangeNumbering", "(", new PETSc::changeNumbering< Dmat, KN >(1, 1, 1));
    Global.Add("ChangeNumbering", "(", new PETSc::changeNumbering< Dmat, KNM >( ));
    Global.Add("ChangeNumbering", "(", new PETSc::changeNumbering< Dmat, KNM >(1, 1, 1, 1));
    Global.Add("MatMult", "(",
               new OneOperator3_< long, Dmat*, KN< PetscScalar >*, KN< PetscScalar >* >(
                 PETSc::MatMult< 'N' >));
    Global.Add("MatMatMult", "(",
               new OneOperator3_< long, Dmat*, KNM< PetscScalar >*, KNM< PetscScalar >* >(
                 PETSc::MatMatMult< 'N' >));
    Global.Add("MatMatMult", "(",
               new OneOperator3_< long, Dmat*, Dmat*, Dmat* >(
                 PETSc::MatMatMult< 'N' >));
    Global.Add("MatMultTranspose", "(",
               new OneOperator3_< long, Dmat*, KN< PetscScalar >*, KN< PetscScalar >* >(
                 PETSc::MatMult< 'T' >));
    Global.Add("MatTransposeMatMult", "(",
               new OneOperator3_< long, Dmat*, KNM< PetscScalar >*, KNM< PetscScalar >* >(
                 PETSc::MatMatMult< 'T' >));
    if (!std::is_same< PetscScalar, PetscReal >::value) {
      Global.Add("MatMultHermitianTranspose", "(",
                 new OneOperator3_< long, Dmat*, KN< PetscScalar >*, KN< PetscScalar >* >(
                   PETSc::MatMult< 'H' >));
      Global.Add("MatHermitianTransposeMatMult", "(",
                 new OneOperator3_< long, Dmat*, KNM< PetscScalar >*, KNM< PetscScalar >* >(
                   PETSc::MatMatMult< 'H' >));
    }
    Global.Add("MatPtAP", "(",
               new OneOperator3_< long, Dmat*, Dmat*, Dmat* >(PETSc::MatPtAP));
    Global.Add("MatConvert", "(", new PETSc::convert< Dmat >);
    Global.Add("MatZeroRows", "(",
               new OneOperator2_< long, Dmat*, KN< double >* >(PETSc::MatZeroRows));
    Global.Add("KSPSolve", "(", new PETSc::LinearSolver< Dmat >( ));
    Global.Add("KSPSolve", "(", new PETSc::LinearSolver< Dmat >(1));
    Global.Add("KSPSolve", "(", new PETSc::LinearSolver< Dmat >(1, 1));
    if (!std::is_same< PetscScalar, PetscReal >::value)
      Global.Add("KSPSolveHermitianTranspose", "(", new PETSc::LinearSolver< Dmat >(1, 1, 1));
    Global.Add("KSPSolve", "(", new PETSc::LinearSolver< Dmat >(1, 1, 1, 1));
    Global.Add("KSPSolveTranspose", "(", new PETSc::LinearSolver< Dmat >(1, 1, 1, 1, 1));
    Global.Add("KSPGetConvergedReason", "(", new OneOperator1_< long, Dmat* >(PETSc::GetConvergedReason< Dmat >));
    Global.Add("KSPGetIterationNumber", "(", new OneOperator1_< long, Dmat* >(PETSc::GetIterationNumber< Dmat >));
    Global.Add("KSPSetResidualHistory", "(", new OneOperator2_< long, Dmat*, KN< double >* >(PETSc::SetResidualHistory< Dmat >));
    Global.Add("SNESSolve", "(", new PETSc::NonlinearSolver< Dmat >(1));
    Global.Add("SNESSolve", "(", new PETSc::NonlinearSolver< Dmat >( ));
    Global.Add("TSSolve", "(", new PETSc::NonlinearSolver< Dmat >(1, 1));
    Global.Add("TSSolve", "(", new PETSc::NonlinearSolver< Dmat >(1, 1, 1));
    Global.Add("TaoSolve", "(", new PETSc::NonlinearSolver< Dmat >(4));
    Global.Add("GlobalNumbering", "(",
               new OneOperator2_< long, Dmat*, KN< long >* >(PETSc::globalNumbering< Dmat >));
    Global.Add("GlobalNumbering", "(",
               new OneOperator2_< long, Dmat*, KN< double >* >(PETSc::globalNumbering< Dmat >));
    Global.Add("GlobalNumbering", "(",
               new OneOperator2_< long, Dbddc*, KN< long >* >(PETSc::globalNumbering< Dbddc >));
    Global.Add("ParMmgCommunicators", "(",
               new OneOperator4_< long, Dmat*, KN< double >*, KN< long >*, KN< KN< long >>*>(PETSc::ParMmgCommunicators< Dmat >));
    Global.Add("ChangeSchur", "(",
               new OneOperator3_< long, Dmat*, KN< Matrice_Creuse< upscaled_type<PetscScalar> > >*, KN< double >* >(
                 PETSc::changeSchur));
    Global.Add("ObjectView", "(", new PETSc::view< Dmat, 0 >);
    Global.Add("ObjectView", "(", new PETSc::view< KNM<PetscScalar>, 0 >);
    Global.Add("MatLoad", "(", new PETSc::view< Dmat, 1 >);
    Global.Add("MatLoad", "(", new PETSc::view< KNM<PetscScalar>, 1 >);
    Global.Add(
      "OriginalNumbering", "(",
      new OneOperator3_< long, Dbddc*, KN< PetscScalar >*, KN< long >* >(PETSc::originalNumbering));
    Global.Add("renumber", "(",
               new OneOperator3_< long, KN< PetscScalar >*, KN< long >*, KN< PetscScalar >* >(
                 PETSc::renumber));
    Global.Add("set", "(", new PETSc::setOptions< Dbddc >( ));
    addInv< Dbddc, PETSc::InvPETSc, KN< PetscScalar >, PetscScalar >( );
#ifdef GENERATE_DEPRECATED_FUNCTIONS
    Global.Add("changeNumbering", "(", new PETSc::changeNumbering< Dmat, KN >( ));
    Global.Add("changeNumbering", "(", new PETSc::changeNumbering< Dmat, KN >(1));
    Global.Add("changeNumbering", "(", new PETSc::changeNumbering< Dmat, KN >(1, 1));
    Global.Add("changeNumbering", "(", new PETSc::changeNumbering< Dmat, KN >(1, 1, 1));
    Global.Add("changeNumbering", "(", new PETSc::changeNumbering< Dmat, KNM >( ));
    Global.Add("globalNumbering", "(",
               new OneOperator2_< long, Dmat*, KN< long >* >(PETSc::globalNumbering< Dmat >));
    Global.Add("globalNumbering", "(",
               new OneOperator2_< long, Dmat*, KN< double >* >(PETSc::globalNumbering< Dmat >));
    Global.Add("globalNumbering", "(",
               new OneOperator2_< long, Dbddc*, KN< long >* >(PETSc::globalNumbering< Dbddc >));
    Global.Add(
      "originalNumbering", "(",
      new OneOperator3_< long, Dbddc*, KN< PetscScalar >*, KN< long >* >(PETSc::originalNumbering));
#endif
  }
  template<class K, typename std::enable_if<!std::is_same<K, upscaled_type<K>>::value>::type* = nullptr>
  static void init() { }
  class initDM_Op : public E_F0mps {
      public:
          Expression A;
          Expression B;
          static const int n_name_param = 6;
          static basicAC_F0::name_and_type name_param[];
          Expression nargs[n_name_param];
          initDM_Op(const basicAC_F0& args, Expression param1, Expression param2) : A(param1), B(param2) {
              args.SetNameParam(n_name_param, name_param, nargs);
          }

          AnyType operator()(Stack stack) const;
  };
  basicAC_F0::name_and_type initDM_Op::name_param[] = {
      {"communicator", &typeid(pcommworld)},
      {"overlap", &typeid(long)},
      {"neighbors", &typeid(KN<long>*)},
      {"sparams", &typeid(string*)},
      {"partition", &typeid(KN<long>*)},
      {"prefix", &typeid(string*)}
  };
  class initDM : public OneOperator {
      public:
          initDM() : OneOperator(atype<DMPlex*>(), atype<DMPlex*>(), atype<string*>()) { }

          E_F0* code(const basicAC_F0& args) const {
              return new initDM_Op(args, t[0]->CastTo(args[0]), t[1]->CastTo(args[1]));
          }
  };
  AnyType initDM_Op::operator()(Stack stack) const {
      DMPlex* pA = GetAny<DMPlex*>((*A)(stack));
      std::string* pB = GetAny<std::string*>((*B)(stack));
      std::string* prefix = nargs[5] ? GetAny< std::string* >((*nargs[5])(stack)) : NULL;
      KN<long>* part = nargs[4] ? GetAny< KN<long>* >((*nargs[4])(stack)) : NULL;
      KN<long>* neighbors = nargs[2] ? GetAny< KN<long>* >((*nargs[2])(stack)) : NULL;
      if(nargs[3]) {
          std::string* options = GetAny< std::string* >((*nargs[3])(stack));
          PetscOptionsInsertString(NULL, options->c_str());
      }
      PetscInt overlap = nargs[1] ? GetAny< long >((*nargs[1])(stack)) : 0;
      MPI_Comm comm = nargs[0] ? *static_cast< MPI_Comm* >(GetAny< pcommworld >((*nargs[0])(stack))) : PETSC_COMM_WORLD;
      int size;
      MPI_Comm_size(comm, &size);
      DMPlexCreateFromFile(comm, pB->c_str(), NULL, PETSC_TRUE, &pA->_dm);
      if(prefix)
          DMSetOptionsPrefix(pA->_dm, prefix->c_str());
      DMSetFromOptions(pA->_dm);
      PetscPartitioner partitioner;
      DMPlexGetPartitioner(pA->_dm, &partitioner);
      PetscPartitionerSetFromOptions(partitioner);
      DM pdm;
      DMPlexDistribute(pA->_dm, 0, NULL, &pdm);
      if (pdm) {
          if(overlap == 0) {
              DM idm;
              DMPlexInterpolate(pdm, &idm);
              DMDestroy(&pdm);
              pdm = idm;
          }
          DMLabel label;
          DMGetLabel(pdm, "Face Sets", &label);
          if (!label) {
              DMCreateLabel(pdm, "Face Sets");
              DMGetLabel(pdm, "Face Sets", &label);
          }
#if PETSC_VERSION_GE(3,19,0)
          PetscSF sf;
          DMGetPointSF(pdm, &sf);
          DMSetPointSF(pdm, NULL);
#endif
          DMPlexMarkBoundaryFaces(pdm, 111111, label);
#if PETSC_VERSION_GE(3,19,0)
          DMSetPointSF(pdm, sf);
#endif
          DMDestroy(&pA->_dm);
          pA->_dm = pdm;
          if(neighbors) {
              PetscInt nranks;
              const PetscMPIInt *ranks;
              DMGetNeighbors(pA->_dm, &nranks, &ranks);
              neighbors->resize(nranks);
              std::copy_n(ranks, nranks, neighbors->operator long*());
          }
      } else if(neighbors) neighbors->resize(0);
      if(part) {
          PetscSection rootSection, leafSection;
          PetscSectionCreate(comm, &rootSection);
          PetscSectionCreate(comm, &leafSection);
          {
              PetscInt pStart, pEnd, cStart, cEnd;
              DMPlexGetChart(pA->_dm, &pStart, &pEnd);
              DMPlexGetHeightStratum(pA->_dm, 0, &cStart, &cEnd);
              PetscSectionSetChart(rootSection, pStart, pEnd);
              for(PetscInt c = cStart; c < cEnd; ++c)
                  PetscSectionSetDof(rootSection, c, 1);
              PetscSectionSetUp(rootSection);
              PetscSectionCopy(rootSection, leafSection);
              DMSetLocalSection(pA->_dm, rootSection);
              DMSetLocalSection(pA->_dm, leafSection);
          }
          Vec ranks, local;
          DMPlexCreateRankField(pA->_dm, &ranks);
          DMCreateLocalVector(pA->_dm, &local);
          DMGlobalToLocal(pA->_dm, ranks, INSERT_VALUES, local);
          const PetscScalar* val;
          VecGetArrayRead(local, &val);
          PetscInt n;
          VecGetLocalSize(local, &n);
          part->resize(n);
          for(PetscInt i = 0; i < n; ++i)
              part->operator[](i) = PetscRealPart(val[i]);
          VecRestoreArrayRead(local, &val);
          VecDestroy(&ranks);
          VecDestroy(&local);
          PetscSectionDestroy(&rootSection);
          PetscSectionDestroy(&leafSection);
      }
      return pA;
  }
  static void findPerm(int *ivt, int *pivt, Vertex3 *v) {
    std::copy(ivt, ivt + 4, pivt);
    R3 AB(v[ivt[0]], v[ivt[1]]);
    R3 AC(v[ivt[0]], v[ivt[2]]);
    R3 AD(v[ivt[0]], v[ivt[3]]);
    if (det(AB, AC, AD) > 0) {
      return;
    } else if (det(AB, AD, AC) > 0) {
      std::swap(pivt[2], pivt[3]);
    } else if (det(AC, AB, AD) > 0) {
      std::swap(pivt[1], pivt[2]);
    }
  }
  static void findPerm(int *ivt, int *pivt, Vertex *v) {
    std::copy(ivt, ivt + 3, pivt);
    R2 A(v[ivt[0]]);
    R2 B(v[ivt[1]]);
    R2 C(v[ivt[2]]);
    if (((B-A)^(C-A)) < 0)
      std::swap(pivt[1], pivt[2]);
  }
  class DMPlexToFF : public OneOperator {
   public:
    const int c;
    class DMPlexToFF_Op : public E_F0mps {
     public:
      Expression A;
      Expression B;
      const int c;
      DMPlexToFF_Op(const basicAC_F0& args, int d) : A(0), B(0), c(d) {
        B = to<DMPlex*>(args[1]);
        if(c == 0)
            A = to< Mesh const** >(args[0]);
        else
            A = to< Mesh3 const** >(args[0]);
      }

      AnyType operator( )(Stack stack) const;
      operator aType( ) const { return c == 0 ? atype< Mesh const** >( ) : atype < Mesh3 const** >( ); }
    };
    E_F0* code(const basicAC_F0& args) const { return new DMPlexToFF_Op(args, c); }
    DMPlexToFF( )
      : OneOperator(atype< Mesh const** >( ), atype< Mesh const** >( ), atype< DMPlex* >( )), c(0) {}
    DMPlexToFF(int)
      : OneOperator(atype< Mesh3 const** >( ), atype< Mesh3 const** >( ), atype< DMPlex* >( )), c(1) {}
  };
  AnyType DMPlexToFF::DMPlexToFF_Op::operator( )(Stack stack) const {
      DMPlex* p = GetAny<DMPlex*>((*B)(stack));
      PetscSection      coordSection;
      Vec               coordinates;
      DMLabel           label;
      const PetscScalar *a;
      PetscInt          dim, pStart, pEnd, cStart, cEnd, vStart, vEnd, depth, numValues;
      IS                valueIS;
      const PetscInt    *values;
      DMGetDimension(p->_dm, &dim);
      DMGetCoordinatesLocal(p->_dm, &coordinates);
      DMGetCoordinateSection(p->_dm, &coordSection);
      DMPlexGetDepth(p->_dm, &depth);
      DMPlexGetHeightStratum(p->_dm, depth, &vStart, &vEnd);
      DMPlexGetHeightStratum(p->_dm, 0, &cStart, &cEnd);
      PetscSectionGetChart(coordSection, &pStart, &pEnd);
      VecGetArrayRead(coordinates, &a);
      if(c == 1) {
          Mesh3** pA = GetAny<Mesh3**>((*A)(stack));
          (*pA)->destroy();
          if(!p->_dm) return pA;
          ffassert(dim == 3);
          Vertex3 *v = new Vertex3[vEnd - vStart];
          Tet *t = new Tet[cEnd - cStart];
          Tet *tt = t;
          for (PetscInt i = 0; i < vEnd - vStart; ++i) {
              v[i].x = PetscRealPart(a[3 * i + 0]);
              v[i].y = PetscRealPart(a[3 * i + 1]);
              v[i].z = PetscRealPart(a[3 * i + 2]);
              v[i].lab = 0;
          }
          VecRestoreArrayRead(coordinates, &a);
          DMGetLabel(p->_dm, "Face Sets", &label);
          DMLabelGetNumValues(label, &numValues);
          DMLabelGetValueIS(label, &valueIS);
          ISGetIndices(valueIS, &values);
          PetscInt interface = -1;
          std::unordered_set<PetscInt> set;
          for (PetscInt v = 0; v < numValues; ++v) {
              if (values[v] == 111111) {
                  interface = v;
                  continue;
              }
              IS face;
              PetscInt size;
              const PetscInt* indices;

              DMLabelGetStratumIS(label, values[v], &face);
              ISGetLocalSize(face, &size);
              ISGetIndices(face, &indices);
              for(PetscInt i = 0; i < size;++i)
                  set.insert(indices[i]);
              ISRestoreIndices(face, &indices);
              ISDestroy(&face);
          }
          PetscInt size = set.size();
          if (interface != -1) {
              IS face;
              PetscInt iss;
              const PetscInt* indices;

              DMLabelGetStratumIS(label, values[interface], &face);
              ISGetLocalSize(face, &iss);
              ISGetIndices(face, &indices);
              for(PetscInt i = 0; i < iss;++i) {
                  if (set.find(indices[i]) == set.end())
                      ++size;
              }
              ISRestoreIndices(face, &indices);
              ISDestroy(&face);
          }
          Triangle3 *b = new Triangle3[size];
          Triangle3 *bb = b;
          for (PetscInt j = 0; j < numValues; ++j) {
              IS face;
              PetscInt size;
              const PetscInt* indices;

              DMLabelGetStratumIS(label, values[j], &face);
              ISGetLocalSize(face, &size);
              ISGetIndices(face, &indices);
              for (PetscInt c = 0; c < size; ++c) {
                  if (set.find(indices[c]) != set.end() && j == interface) continue;
                  const PetscInt *points, *orientations;
                  PetscInt       size, i;

                  DMPlexGetConeSize(p->_dm, indices[c], &size);
                  DMPlexGetCone(p->_dm, indices[c], &points);
                  std::set<PetscInt> set;
                  for (PetscInt j = 0; j < size; ++j) {
                      const PetscInt *vertex;
                      PetscInt       size, i;

                      DMPlexGetConeSize(p->_dm, points[j], &size);
                      DMPlexGetCone(p->_dm, points[j], &vertex);
                      for (PetscInt w = 0; w < size; ++w)
                          set.insert(vertex[w]);
                  }
                  std::vector<PetscInt> conv(set.begin(), set.end());
                  ffassert(conv.size() == 3);
                  int iv[3];
                  for (PetscInt j = 0; j < conv.size(); ++j)
                      iv[j] = conv[j] - vStart;
                  int lab = values[j] == 111111 ? -111111 : values[j];
                  bb++->set(v, iv, lab);
              }
              ISRestoreIndices(face, &indices);
              ISDestroy(&face);
          }
          ISRestoreIndices(valueIS, &values);
          ISDestroy(&valueIS);
          for (PetscInt c = cStart; c < cEnd; ++c) {
              PetscInt *closure = NULL;
              PetscInt  closureSize, cl;
              DMPlexGetTransitiveClosure(p->_dm, c, PETSC_TRUE, &closureSize, &closure);
              int iv[4];
              PetscInt lab;
              DMGetLabel(p->_dm, "Cell Sets", &label);
              if (label) DMGetLabelValue(p->_dm, "Cell Sets", c, &lab);
              else {
                DMGetLabel(p->_dm, "celltype", &label);
                if (label) DMGetLabelValue(p->_dm, "celltype", c, &lab);
              }
              if (lab == -1) lab = 0;
              int* ivv = iv;
              for (cl = 0; cl < 2 * closureSize; cl += 2) {
                  PetscInt point = closure[cl], dof, off, d, p;

                  if ((point < pStart) || (point >= pEnd)) continue;
                  PetscSectionGetDof(coordSection, point, &dof);
                  if (!dof) continue;
                  PetscSectionGetOffset(coordSection, point, &off);
                  *ivv++ = point - cEnd;
              }
              ffassert(ivv - iv == 4);
              int pivt[4];
              findPerm(iv, pivt, v);
              tt++->set(v, pivt, lab);
              DMPlexRestoreTransitiveClosure(p->_dm, c, PETSC_TRUE, &closureSize, &closure);
          }
          *pA = new Mesh3(vEnd - vStart, cEnd - cStart, size, v, t, b);
          (*pA)->BuildGTree();
          return pA;
      }
      else {
          Mesh const** pA = GetAny<Mesh const**>((*A)(stack));
          (*pA)->destroy();
          if(!p->_dm) return pA;
          ffassert(dim == 2);
          Vertex *v = new Vertex[vEnd - vStart];
          Triangle *t = new Triangle[cEnd - cStart];
          Triangle *tt = t;
          for (PetscInt i = 0; i < vEnd - vStart; ++i) {
              v[i].x = PetscRealPart(a[2 * i + 0]);
              v[i].y = PetscRealPart(a[2 * i + 1]);
              v[i].lab = 0;
          }
          VecRestoreArrayRead(coordinates, &a);
          DMGetLabel(p->_dm, "Face Sets", &label);
          DMLabelGetNumValues(label, &numValues);
          DMLabelGetValueIS(label, &valueIS);
          ISGetIndices(valueIS, &values);
          PetscInt interface = -1;
          std::unordered_set<PetscInt> set;
          for (PetscInt v = 0; v < numValues; ++v) {
              if (values[v] == 111111) {
                  interface = v;
                  continue;
              }
              IS face;
              PetscInt size;
              const PetscInt* indices;

              DMLabelGetStratumIS(label, values[v], &face);
              ISGetLocalSize(face, &size);
              ISGetIndices(face, &indices);
              for(PetscInt i = 0; i < size;++i)
                  set.insert(indices[i]);
              ISRestoreIndices(face, &indices);
              ISDestroy(&face);
          }
          PetscInt size = set.size();
          if (interface != -1) {
              IS face;
              PetscInt iss;
              const PetscInt* indices;

              DMLabelGetStratumIS(label, values[interface], &face);
              ISGetLocalSize(face, &iss);
              ISGetIndices(face, &indices);
              for(PetscInt i = 0; i < iss;++i) {
                  if (set.find(indices[i]) == set.end())
                      ++size;
              }
              ISRestoreIndices(face, &indices);
              ISDestroy(&face);
          }
          BoundaryEdge *b = new BoundaryEdge[size];
          BoundaryEdge *bb = b;
          for (PetscInt j = 0; j < numValues; ++j) {
              IS face;
              PetscInt size;
              const PetscInt* indices;

              DMLabelGetStratumIS(label, values[j], &face);
              ISGetLocalSize(face, &size);
              ISGetIndices(face, &indices);
              for (PetscInt c = 0; c < size; ++c) {
                  if (set.find(indices[c]) != set.end() && j == interface) continue;
                  const PetscInt *points, *orientations;
                  PetscInt       size, i;

                  DMPlexGetConeSize(p->_dm, indices[c], &size);
                  DMPlexGetCone(p->_dm, indices[c], &points);
                  ffassert(size == 2);
                  int iv[2];
                  for (PetscInt j = 0; j < size; ++j)
                      iv[j] = points[j] - vStart;
                  int lab = values[j] == 111111 ? -111111 : values[j];
                  *bb++ = BoundaryEdge(v, iv[0], iv[1], lab);
              }
              ISRestoreIndices(face, &indices);
              ISDestroy(&face);
          }
          ISRestoreIndices(valueIS, &values);
          ISDestroy(&valueIS);
          DMGetLabel(p->_dm, "Cell Sets", &label);
          if (label) DMLabelSetDefaultValue(label, 0);
          else {
            DMGetLabel(p->_dm, "celltype", &label);
            if (label) DMLabelSetDefaultValue(label, 0);
          }
          for (PetscInt c = cStart; c < cEnd; ++c) {
              PetscInt *closure = NULL;
              PetscInt  closureSize, cl;
              DMPlexGetTransitiveClosure(p->_dm, c, PETSC_TRUE, &closureSize, &closure);
              int iv[3], lab = 0;
              if (label) {
                DMLabelGetValue(label, c, &cl);
                lab = cl;
              }
              int* ivv = iv;
              for (cl = 0; cl < 2 * closureSize; cl += 2) {
                  PetscInt point = closure[cl], dof, off, d, p;

                  if ((point < pStart) || (point >= pEnd)) continue;
                  PetscSectionGetDof(coordSection, point, &dof);
                  if (!dof) continue;
                  PetscSectionGetOffset(coordSection, point, &off);
                  *ivv++ = point - cEnd;
              }
              ffassert(ivv - iv == 3);
              int pivt[3];
              findPerm(iv, pivt, v);
              tt++->set(v, pivt[0], pivt[1], pivt[2], lab);
              DMPlexRestoreTransitiveClosure(p->_dm, c, PETSC_TRUE, &closureSize, &closure);
          }
          Mesh* m = new Mesh(vEnd - vStart, cEnd - cStart, size, v, t, b);
          R2 Pn, Px;
          m->BoundingBox(Pn, Px);
          m->quadtree = new Fem2D::FQuadTree(m, Pn, Px, m->nv);
          *pA = m;
          return pA;
      }
  }
  class buildSolution : public OneOperator {
   public:
    class buildSolution_Op : public E_F0mps {
     public:
      Expression A;
      Expression B;
      buildSolution_Op(const basicAC_F0& args) : A(0), B(0) {
        A = to<Dmat*>(args[0]);
        B = to<KN<upscaled_type<PetscScalar>>*>(args[1]);
      }

      AnyType operator( )(Stack stack) const;
      operator aType( ) const { return atype< long >( ); }
    };
    E_F0* code(const basicAC_F0& args) const { return new buildSolution_Op(args); }
    buildSolution( )
      : OneOperator(atype<long>( ), atype<Dmat*>( ), atype<KN<upscaled_type<PetscScalar>>*>( )) {}
  };
  AnyType buildSolution::buildSolution_Op::operator( )(Stack stack) const {
      Dmat* ptA = GetAny<Dmat*>((*A)(stack));
      KN<upscaled_type<PetscScalar>>* ptKN = GetAny<KN<upscaled_type<PetscScalar>>*>((*B)(stack));
      ffassert(ptA->_ksp);
      Mat A;
      KSPGetOperators(ptA->_ksp, &A, NULL);
      PetscInt n, N;
      MatGetLocalSize(A, &n, NULL);
      MatGetSize(A, &N, NULL);
      ptKN->resize(n);
      Vec v;
      PetscScalar* p = reinterpret_cast<PetscScalar*>(ptKN->operator upscaled_type<PetscScalar>*());
      VecCreateMPIWithArray(PetscObjectComm((PetscObject)A), 1, n, N, p, &v);
      KSPBuildSolution(ptA->_ksp, v, NULL);
      if(!std::is_same<upscaled_type<PetscReal>, PetscReal>::value) {
        for(int i = n - 1; i >= 0; --i)
          ptKN->operator[](i) = p[i];
      }
      VecDestroy(&v);
      return 0L;
  }
} // namespace PETSc

template< class R, class FESpaceT1, class FESpaceT2 >
Matrice_Creuse<R> *  PETSC_buildMatrixInterpolationForCompositeFESpace(const FESpaceT1 * Uh ,const FESpaceT2 * Vh, bool transpose=false){
  ffassert(Uh);
  ffassert(Vh);
  int NUh = Uh->N;
  int NVh = Vh->N;

  if(verbosity>3) cout << "NUh=" << NUh << ", NVh=" << NVh << endl;
  Matrice_Creuse<R> * sparse_mat= new Matrice_Creuse<R>();

  // Remarque pas de U2Vc pour l'instant
  int* data = new int[4 + NUh];
  // default value for the interpolation matrix
  data[0]=transpose;         // transpose not
  data[1]=(long) op_id;  // get just value
  data[2]=false;         // get just value
  data[3]=0L;            // get just value

  for(int i=0;i<NUh;++i) data[4+i]=i;//

  if(verbosity>3){
    for(int i=0;i<NUh;++i)
    {
      cout << "The Uh componante " << i << " -> " << data[4+i] << "  Componante of Vh  " <<endl;
    }
  }
  for(int i=0;i<NUh;++i){
    if(data[4+i]>=NVh)
    {
      cout << "The Uh componante " << i << " -> " << data[4+i] << " >= " << NVh << " number of Vh Componante " <<endl;
      ExecError("Interpolation incompability between componante ");
    }
  }
  const FESpaceT1 &rUh = *Uh;
  const FESpaceT2 &rVh = *Vh;

  MatriceMorse<R>* titi=buildInterpolationMatrixT<FESpaceT1,FESpaceT2>(rUh,rVh,data);

  sparse_mat->init();
  sparse_mat->typemat=0;//(TypeSolveMat::NONESQUARE); //  none square matrice (morse)
  sparse_mat->A.master( titi );	  //  sparse_mat->A.master(new MatriceMorse<R>(*Uh,*Vh,buildInterpolationMatrix,data));
  if(verbosity>3){
    cout << "sparse_mat->typemat=" << sparse_mat->typemat << endl;
    cout << "N=" << sparse_mat->A->n << endl;
    cout << "M=" << sparse_mat->A->m << endl;
  }
  delete [] data;

  return sparse_mat;
}

template<class HpddmType>  //  to make   A=linearform(x)
struct OpMatrixtoBilinearFormVGPETSc
  : public OneOperator
{
  typedef typename Call_CompositeFormBilinear<vect_generic_v_fes,vect_generic_v_fes>::const_iterator const_iterator;
  int init;
  
  class Op : public E_F0mps {
    public:
      Call_CompositeFormBilinear<vect_generic_v_fes,vect_generic_v_fes> *b;
      Expression a;
      int init;
      //AnyType operator()(Stack s)  const;
      
      Op(Expression aa,Expression  bb,int initt)
        : b(new Call_CompositeFormBilinear<vect_generic_v_fes,vect_generic_v_fes>(* dynamic_cast<const Call_CompositeFormBilinear<vect_generic_v_fes,vect_generic_v_fes> *>(bb))),a(aa),init(initt)
    { 
      assert(b && b->nargs);
      int NN = (int) b->euh->componentNbitem().size();
      int MM = (int) b->evh->componentNbitem().size();

      bool total_iscmplx=false;
      // loop over block
      for(int i=0; i<NN; i++){
        for(int j=0; j<MM; j++){
          // FieldOfForm : optimize the terms (flags -O3) of the variational form and verifies the type of the variational form
          bool iscmplx=FieldOfForm(b->block_largs(i,j),IsComplexType<PetscScalar>::value)  ;
          // cout<< "FieldOfForm:iscmplx " << iscmplx << " " << IsComplexType<R>::value << " " << ((iscmplx) == IsComplexType<R>::value) << endl;
          ffassert( (iscmplx) == IsComplexType<PetscScalar>::value);
          if( !total_iscmplx ) total_iscmplx=iscmplx;
        }
      }
    }
    operator aType () const { return atype<PETSc::DistributedCSR< HpddmType >*>();}

    AnyType operator()(Stack s)  const;
  };

  E_F0 * code(const basicAC_F0 & args) const
  { return  new Op(to<PETSc::DistributedCSR< HpddmType >*>(args[0]),args[1],init); }
  OpMatrixtoBilinearFormVGPETSc(int initt=0) :
    OneOperator(atype<PETSc::DistributedCSR< HpddmType >*>(),atype<PETSc::DistributedCSR< HpddmType >*>(),atype<const Call_CompositeFormBilinear<vect_generic_v_fes,vect_generic_v_fes>*>()),
    init(initt){};

};

// function to transform freefem matrix in matIS.
void ff_createMatIS( MatriceMorse<PetscScalar> &ff_mat, Mat &matIS){
    std::set<PetscInt> irows;
    std::set<PetscInt> jcols;

    ff_mat.CSR(); // transform the matrix to CSR format
    
    std::vector<PetscInt> perm_row(ff_mat.n,-1); 
    std::vector<PetscInt> perm_col(ff_mat.m,-1);
      
    for (int ii=0; ii < ff_mat.n; ii++) {
      for (int la = ff_mat.p[ii]; la < ff_mat.p[ii+1]; la++) {
        perm_row[ii] = 1; 
        perm_col[ff_mat.j[la]] = 1;
        if( la > ff_mat.p[ii]){
          //
          if( !( ff_mat.j[la] > ff_mat.j[ff_mat.p[ii]]) ){
            cerr << " The column index must be croissant :: Error " << ff_mat.j[la] << " " << ff_mat.j[ff_mat.p[ii]]  <<endl;
          }
          ffassert( ff_mat.j[la] > ff_mat.j[ff_mat.p[ii]] ); 
        }
      }
    }
  
  if(verbosity> 2)
    for(int ii=0; ii<10; ii++)
      cout << "perm_row["<<ii<<"]=" << perm_row[ii] << endl;

  // construction of irows
  for (int ii=0; ii < ff_mat.n; ii++) 
    if (perm_row[ii] == 1) {
      auto it = irows.insert(ii);
      perm_row[ii] = std::distance(irows.begin(),it.first);
    }
  // construction of jcols
  for (int ii=0; ii < ff_mat.m; ii++) 
    if (perm_col[ii] == 1) {
      auto it = jcols.insert(ii);
      perm_col[ii] = std::distance(jcols.begin(),it.first);
    }

  if(verbosity>2){
    MPI_Barrier(PETSC_COMM_WORLD);
    if(mpirank ==1){
      cout << "nnz=" << ff_mat.nnz << endl;
      cout << "n=" << ff_mat.n << endl;
      //for (int ii=0; ii < ff_mat.n; ii++) 
      //  cout << "ii=" << ii << " " << perm[ii] << endl;   
    }
    MPI_Barrier(PETSC_COMM_WORLD);
  }
  PetscInt *IA = new PetscInt[irows.size()+1];
  PetscInt *JA = new PetscInt[ff_mat.nnz];
  PetscScalar* aa = new PetscScalar[ff_mat.nnz];

  int cpt = 0;
  IA[0] = 0;
  for (int ii=0; ii < ff_mat.n; ii++)
    if (ff_mat.p[ii] != ff_mat.p[ii+1]){
      IA[cpt] = ff_mat.p[ii];
      ffassert( perm_row[ii] == cpt);
      cpt++;
    }
  IA[ cpt ] = ff_mat.nnz; 
  ffassert( IA[cpt] == ff_mat.nnz);
  ffassert(cpt==irows.size());

  for (int ii=0; ii < ff_mat.nnz; ii++) {
    JA[ii] = perm_col[ff_mat.j[ii]];
  }

  std::copy_n(ff_mat.aij, ff_mat.nnz, aa);

  std::vector<PetscInt> indices_row; indices_row.reserve(irows.size());
  for(const auto& p : irows) indices_row.emplace_back(p);

  std::vector<PetscInt> indices_col; indices_col.reserve(jcols.size());
  for(const auto& p : jcols) indices_col.emplace_back(p);

  ISLocalToGlobalMapping mapping_row;
  ISLocalToGlobalMappingCreate(PETSC_COMM_WORLD, 1, indices_row.size(), indices_row.data(), PETSC_COPY_VALUES, &mapping_row);

  ISLocalToGlobalMapping mapping_col;
  ISLocalToGlobalMappingCreate(PETSC_COMM_WORLD, 1, indices_col.size(), indices_col.data(), PETSC_COPY_VALUES, &mapping_col);

  Mat matISlocal, matIJ;

  // Remark: If we put the true size of the global matrix ==> Error : Abort trap
  MatCreateIS(PETSC_COMM_WORLD, 1, PETSC_DECIDE, PETSC_DECIDE, ff_mat.n, ff_mat.m, mapping_row, mapping_col, &matIS);

  // This call supposed that local and global indices are the same
  // MatCreateSeqAIJ(PETSC_COMM_SELF, irows.size(), irows.size(), 0, nullptr, &matISlocal);
  // MatCreateSeqAIJ(PETSC_COMM_SELF, A.pHM()->n, A.pHM()->m, 0, nullptr, &matISlocal);

  ISLocalToGlobalMappingDestroy(&mapping_row);
  ISLocalToGlobalMappingDestroy(&mapping_col);

  // This 4 lines are equivalent to MatCreateSeqAIJ
  MatCreate(PETSC_COMM_SELF, &matISlocal);
  MatSetType(matISlocal,MATSEQAIJ);
  MatSetSizes(matISlocal, irows.size(), jcols.size(), irows.size(), jcols.size());
  // nullptr can be replaced by the vector of nnz
  MatSeqAIJSetPreallocation(matISlocal, 0, nullptr);

  if(verbosity>2){
    MPI_Barrier(PETSC_COMM_WORLD);
    if(mpirank==0) cout << " matISlocal" <<  mpirank << endl;
    if(mpirank==0) cout << " irows.size()=" << irows.size() << endl;
    if(mpirank==0) cout << " jcols.size()=" << jcols.size() << endl;
    MPI_Barrier(PETSC_COMM_WORLD);
  }

  MatSeqAIJSetPreallocationCSR(matISlocal, IA, JA, aa);
  //MatMPIAIJSetPreallocationCSR(matISlocal, IA, JA, aa);

  MatISSetLocalMat(matIS, matISlocal);

  MatAssemblyBegin(matIS, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(matIS, MAT_FINAL_ASSEMBLY);

  MatConvert(matIS, MATAIJ, MAT_INPLACE_MATRIX, &matIS);

  // delete [] IA;
  // delete [] JA;
  // delete [] aa;

  delete IA;
  delete JA;
  delete aa;
  MatDestroy(&matISlocal);
}

template<class HpddmType>
AnyType OpMatrixtoBilinearFormVGPETSc<HpddmType>::Op::operator()(Stack stack) const
{
  typedef PetscScalar R;
  assert(b && b->nargs);

  pvectgenericfes  * pUh= GetAny<pvectgenericfes *>((*b->euh)(stack));
  pvectgenericfes  * pVh= GetAny<pvectgenericfes *>((*b->evh)(stack));

  ffassert( *pUh && *pVh ); 
  // Update is necessary when we get "pvectgenericfes" to take account a new mesh.
  (*pUh)->update();
  (*pVh)->update();

  if( verbosity > 5){
    (*pUh)->printPointer();
    (*pVh)->printPointer();
  }
  int NpUh = (*pUh)->N; // number of fespace in pUh
  int NpVh = (*pVh)->N; // number of fespace in pVh

  KN<int> UhNbOfDf = (*pUh)->vectOfNbOfDF(); // A changer en long
  KN<int> VhNbOfDf = (*pVh)->vectOfNbOfDF();

  KN<int> UhNbItem = (*pUh)->vectOfNbitem();
  KN<int> VhNbItem = (*pVh)->vectOfNbitem();

  //
  const KNM<list<C_F0>> & block_largs=b->block_largs; 

  // check if we have a square matrix
  bool A_is_square= (void*)pUh == (void*)pVh || ((*pUh)->totalNbOfDF()) == ( (*pVh)->totalNbOfDF()) ;

  // === simple check if A is symetrical === // 
  // voir avec les autres.
  bool A_is_maybe_sym = (void*)pUh == (void*)pVh; 

  // VF == true => VF type of Matrix
  //bool VF=isVF(b->block_largs);    //=== used to set the solver ??? block matrix ??? ===/
  bool VF = 0;

  // set parameteer of the matrix :: 
  Data_Sparse_Solver ds;
  ds.factorize=0;
  ds.initmat=true;
  int np = OpCall_FormBilinear_np::n_name_param - NB_NAME_PARM_HMAT;
  SetEnd_Data_Sparse_Solver<R>(stack,ds, b->nargs,np);

  // J'ai repris ce qu'il y avait. 
  // PAC(e)     :: Attention peut être pas compatible avec les matrices bloques.
  // A repenser :: surtout pour le parametre symetrique? on le met ce parametre à zéro pour l'instant.
  // set ds.sym = 0 

  ds.sym = 0;
  if(verbosity>3)
    cout << " === we consider the block matrix as a non symetric matrix === (to be change in the future)" << endl; 

  if (! A_is_square )
   {
     if(verbosity>3) cout << " -- the solver  is un set  on rectangular matrix  " << endl;
    }

  // A quoi cela correspond?? Gestion du stack + autre
  WhereStackOfPtr2Free(stack)=new StackOfPtr2Free(stack);// FH aout 2007
  
  PETSc::DistributedCSR< HpddmType > * Ares( GetAny<PETSc::DistributedCSR< HpddmType >*>((*a)(stack)));

  // test function (Vh) are the line
  // inconnu function (Uh) are the column

  // Assemble the variationnal form
  int maxJVh=NpVh;

  Mat* a = new Mat[NpUh * maxJVh]();
  std::vector<std::pair<int, int>> destroy;
  destroy.reserve(NpUh * maxJVh);

  int offsetMatrixUh = 0;
  // loop over the block
  for( int i=0; i<NpUh; i++){
    int offsetMatrixVh = 0;
    if( ds.sym > 0 ){ maxJVh=(i+1); ffassert(maxJVh<NpVh);}
    for( int j=0; j<maxJVh; j++){
      if(mpirank ==0 && verbosity>3){
      cout << "offsetMatrixUh= " << offsetMatrixUh << ", offsetMatrixVh= " << offsetMatrixVh << endl;
      cout << "construction of block i,j=" << i << ","<< j << endl;
      cout << "                          " << endl;
      }
      MPI_Barrier(PETSC_COMM_WORLD);
      // construction du block (i,j)
      const list<C_F0> & b_largs=block_largs(i,j); 

      // size of the block
      int N_block = UhNbOfDf[i];
      int M_block = VhNbOfDf[j];

      // initialise the bem block and fem block 
      Matrice_Creuse<R> Afem;
      Afem.init();
      Afem.resize(M_block,N_block);
      int nsparseblocks = 0;
      Mat Abem = PETSC_NULL;

      if(verbosity>2) cout << "size_block =" << b_largs.size() << endl; 
      if( b_largs.size()> 0){
        
        // compute largs due BEM and FEM part
        list<C_F0> largs_FEM;
        list<C_F0> largs_BEM;
        largs_FEM.clear();
        largs_BEM.clear();
        separateFEMpartBemPart( b_largs, largs_FEM, largs_BEM );

        if(verbosity>2){
          cout << " FEM.size()=" << largs_FEM.size() << endl;
          cout << " BEM.size()=" << largs_BEM.size() << endl;
        }

        if( largs_BEM.size() >0 ){
          // compute of BEM par
#if defined(WITH_bemtool) && defined(WITH_htool) && defined(PETSC_HAVE_HTOOL)
          // need bemtool, htool
          const list<C_F0> & b_largs_zz = largs_BEM;

          int VFBEM = typeVFBEM(b_largs_zz,stack);
          if(VFBEM == 2){ cerr << " not implemented with BEM POTENTIAL" << endl; ffassert(0);}
          Data_Bem_Solver dsbem;
          dsbem.factorize=0;
          dsbem.initmat=true;
          SetEnd_Data_Bem_Solver<R>(stack, dsbem, b->nargs,OpCall_FormBilinear_np::n_name_param);  // LIST_NAME_PARM_HMAT

          HMatrixVirt<R> ** Hmat = new HMatrixVirt<R> *();
         
          //
          // avoir dans le futur si la difference entre bloc diagonal et bloc non diagonal a un sens.
          //
          if( i==j ){

            bool samemesh = (void*) (*pUh)->vect[i]->getppTh() == (void*) (*pVh)->vect[j]->getppTh();  // same Fem2D::Mesh     +++ pot or kernel
            if (VFBEM==1)
              ffassert (samemesh);

            PETSc::DistributedCSR< HpddmType > * Abemblock = new PETSc::DistributedCSR< HpddmType >;

            // block diagonal matrix
            if( (*pUh)->typeFE[i] == 4 && (*pVh)->typeFE[j] == 4 ){
              ffassert( i==j ); // If not a block diagonal not coded yet
              // MeshS --- MeshS
              // ==== FESpace 3d Surf: inconnue et test ===
              const FESpaceS * PUh = (FESpaceS *) (*pUh)->vect[i]->getpVh();
              const FESpaceS * PVh = (FESpaceS *) (*pVh)->vect[j]->getpVh();

              MatCreate(PETSC_COMM_WORLD, &Abemblock->_petsc);
              MatSetSizes(Abemblock->_petsc, PETSC_DECIDE, PETSC_DECIDE, PUh->NbOfDF, PUh->NbOfDF);
              varfBem<v_fesS, v_fesS>(PUh, PUh, 1, VFBEM, stack, b_largs_zz, dsbem, Abemblock);


            }
            else if( (*pUh)->typeFE[i] == 5 && (*pVh)->typeFE[j] == 5 ){
              ffassert( i==j ); // If not a block diagonal not coded yet
              // MeshL --- MeshL
              // ==== FESpace 3d Curve: inconnue et test ===
              const FESpaceL * PUh = (FESpaceL *) (*pUh)->vect[i]->getpVh();
              const FESpaceL * PVh = (FESpaceL *) (*pVh)->vect[j]->getpVh();

              MatCreate(PETSC_COMM_WORLD, &Abemblock->_petsc);
              MatSetSizes(Abemblock->_petsc, PETSC_DECIDE, PETSC_DECIDE, PUh->NbOfDF, PUh->NbOfDF);
              varfBem<v_fesL, v_fesL>(PUh, PUh, 1, VFBEM, stack, b_largs_zz, dsbem, Abemblock);

            }
            else{
              cerr << " BEM bilinear form " << endl;
              cerr << " Block ("<< i <<" ,"<< j << ")" << endl;
              cerr << " =: Pas prise en compte des FESpace inconnue de type := "<< typeFEtoString( (*pUh)->typeFE[i] ) << endl;
              cerr << " =:                 avec des FESpace test de type    := "<< typeFEtoString( (*pVh)->typeFE[j] ) << endl;
              ffassert(0);
            }
            Abem = Abemblock->_petsc;
          }
          else{

            bool samemesh = (void*) (*pUh)->vect[i]->getppTh() == (void*) (*pVh)->vect[j]->getppTh();  // same Fem2D::Mesh     +++ pot or kernel
          
            if(init)
              *Hmat =0;
            //*Hmat =0;
            if( *Hmat)
              delete *Hmat;
            *Hmat =0;
            
            PETSc::DistributedCSR< HpddmType > * Abemblock = new PETSc::DistributedCSR< HpddmType >;

            // block non diagonal matrix        
            if( (*pUh)->typeFE[i] == 5 && (*pVh)->typeFE[j] == 2 ){
              // case Uh[i] == MeshL et Vh[j] = Mesh2  // Est ce que cela a un sens?
              
              if(verbosity >3) cout << " === creation de la matrice BEM pour un bloc non diagonaux === " << endl;
              const FESpaceL * PUh = (FESpaceL *) (*pUh)->vect[i]->getpVh();

              MatCreate(PETSC_COMM_WORLD, &Abemblock->_petsc);
              MatSetSizes(Abemblock->_petsc, PETSC_DECIDE, PETSC_DECIDE, PUh->NbOfDF, PUh->NbOfDF);
              varfBem<v_fesL, v_fesL>(PUh, PUh, 1, VFBEM, stack, b_largs_zz, dsbem, Abemblock);
            }

            else{
              cerr << " BEM bilinear form " << endl;
              cerr << " Block ("<< i <<" ,"<< j << ")" << endl;
              cerr << " =: Pas prise en compte des FESpace inconnue de type := "<< typeFEtoString( (*pUh)->typeFE[i] ) << endl;
              cerr << " =:                 avec des FESpace test de type    := "<< typeFEtoString( (*pVh)->typeFE[j] ) << endl;
              ffassert(0);
            }
            
            // creation de la matrice dense 
            
            // BEM matrix is constructed with different FESpace
            ffassert( (*pUh)->vect[i]->getpVh() != (*pVh)->vect[j]->getpVh() ) ;
            
            if( (*pUh)->typeFE[i] == 5 && (*pVh)->typeFE[j] == 2 ){
              // case Uh[i] == MeshL et Vh[j] = Mesh2 
              const FESpaceL * PUh = (FESpaceL *) (*pUh)->vect[i]->getpVh();
              const FESpace * PVh = (FESpace *) (*pVh)->vect[j]->getpVh();
              // construction of the matrix of interpolation
              
              // The transpose is in the build matrix now
              bool transpose_MI = false;
              Matrice_Creuse<double> *  MI_BBB = PETSC_buildMatrixInterpolationForCompositeFESpace<double,FESpaceL,FESpace>( PUh, PVh, transpose_MI );
              MatriceMorse<double> * mr =MI_BBB->pHM();
              
              // Transform Real Matrix in Complex Matrix 
              MatriceMorse<R> * mA = new MatriceMorse<R>(mr->n,mr->m,0,0);
              // we divide by mpisize because the interpolation matrix is computed in sequential
              *mr *= 1.0/mpisize;
              *mA = *mr;
              
              // transform intrerpolation matrix to matrix IS
              Mat matIS_MI;
              ff_createMatIS( *mA, matIS_MI);

              if(verbosity>2){
                MPI_Barrier(MPI_COMM_WORLD);
                cout <<"compute the transpose matrix" << endl;
                MPI_Barrier(MPI_COMM_WORLD);
              }

              if( transpose_MI) ffassert(0); // transpose ==true doesn't work. Error Matrice Interpolation
              Mat mAAT;
              if (std::is_same< PetscScalar, PetscReal >::value) MatCreateTranspose(matIS_MI, &mAAT);
              else MatCreateHermitianTranspose(matIS_MI, &mAAT);
              MatDestroy(&matIS_MI);
              
              //  create composite 
              Mat mats[2] = { Abemblock->_petsc , mAAT};
              Mat C;
              MatCreateComposite(PetscObjectComm((PetscObject)Abemblock->_petsc), 2, mats, &C);
              MatCompositeSetType(C, MAT_COMPOSITE_MULTIPLICATIVE);
              Abem = C;

              // Abemblock->dtor();  // ???
              MatDestroy(&mAAT);  // ???
              // we need to do that because R=Complex with BEM
              MI_BBB->destroy();
              delete MI_BBB;
            }
            else{
              cerr << "==== to do ==== " << endl;
              ffassert(0);
            }

            
          }
#endif
        }
        if( largs_FEM.size() >0){
          const list<C_F0> & b_largs_zz = largs_FEM;

          varfToCompositeBlockLinearSystemALLCASE_pfes<R>( i, j, (*pUh)->typeFE[i], (*pVh)->typeFE[j], 
                                                        0, 0, (*pUh)->vect[i], (*pVh)->vect[j],
                                                        true, false, ds.sym, ds.tgv, 
                                                        b_largs_zz, stack, 
                                                        0, 0, Afem.pHM());
         
          nsparseblocks++;
        }
      }
    offsetMatrixVh += VhNbOfDf[j];

    if(b_largs.size()> 0){
      Afem.pHM()->half = ds.sym;
  
      Mat matIS;
      if(nsparseblocks){
        ff_createMatIS( *(Afem.pHM()), matIS);
        Afem.destroy();
      }
    
      if (Abem != PETSC_NULL) {
        if (!nsparseblocks) {
          a[j * maxJVh + i] = Abem;
        }
        else {
          Mat mats[2] = { matIS, Abem };
          Mat C;
          MatCreateComposite(PetscObjectComm((PetscObject)matIS), 2, mats, &C);
          MatCompositeSetType(C, MAT_COMPOSITE_ADDITIVE);
          a[j * maxJVh + i] = C;
        }
      }
      else {
        a[j * maxJVh + i] = matIS; 
      }
      destroy.emplace_back(j, i);
    }
      
    } // end loop j
    offsetMatrixUh += UhNbOfDf[i];
  } // end loop i
  
  // 
  if(verbosity>3) cout << "Ares->_vector_global=" << Ares->_vector_global << endl; 

  if( NpUh==1 && maxJVh==1 ){
    Ares->_petsc = a[0];
    Ares->_vector_global = (PetscBool) 1;
    a = PETSC_NULL; // ???
  }else{
    MatCreateNest(PETSC_COMM_WORLD, NpUh, NULL, maxJVh, NULL, a, &Ares->_petsc);
    Ares->_vector_global = (PetscBool) 1;
    for(std::pair<int, int> p : destroy)
        MatDestroy(a + p.first * maxJVh + p.second);
    delete[] a;
  }
  return SetAny<PETSc::DistributedCSR< HpddmType >*>(Ares);
}

static void Init_PETSc( ) {
  if (verbosity > 1 && mpirank == 0)
    cout << " PETSc (" << typeid(PetscScalar).name( ) << ")" << endl;
  if (exist_type< DmatC* >( )) {
    if(mpirank == 0) cout << "Cannot load both \"PETSc\" and \"PETSc-complex\", please pick a single one" << endl;
    ffassert(0);
  }
  int argc = pkarg->n;
  char** argv = new char*[argc];
  for (int i = 0; i < argc; ++i) argv[i] = const_cast< char* >((*(*pkarg)[i].getap( ))->c_str( ));
  PetscInitialize(&argc, &argv, 0, "");
  PetscSysInitializePackage( );
  MatInitializePackage( );
#ifndef HPDDM_SLEPC
  KSPRegister("hpddm", HPDDM::KSPCreate_HPDDM);
  if (argc > 1) {
    HPDDM::Option& opt = *HPDDM::Option::get( );
    opt.parse(argc - 1, argv + 1, mpirank == 0);
    if (mpirank != 0) opt.remove("verbosity");
  }
#endif
  delete[] argv;
  ff_atend(PETSc::finalizePETSc);
  Dcl_Type< DmatR* >(Initialize< DmatR >, DeleteDTOR< DmatR >);
  zzzfff->Add("Mat", atype< DmatR* >( ));
  if (!exist_type< DmatC* >( )) Dcl_Type< DmatC* >(Initialize< DmatC >, DeleteDTOR< DmatC >);
  Dcl_Type< DbddcR* >(Initialize< DbddcR >, DeleteDTOR< DbddcR >);
  zzzfff->Add("MatIS", atype< DbddcR* >( ));
  if (!exist_type< DbddcC* >( )) Dcl_Type< DbddcC* >(Initialize< DbddcC >, DeleteDTOR< DbddcC >);
  map_type_of_map[make_pair(atype< DmatR* >( ), atype< Complex* >( ))] = atype< DmatC* >( );
  map_type_of_map[make_pair(atype< DmatR* >( ), atype< double* >( ))] = atype< DmatR* >( );
  map_type_of_map[make_pair(atype< DbddcR* >( ), atype< Complex* >( ))] = atype< DbddcC* >( );
  map_type_of_map[make_pair(atype< DbddcR* >( ), atype< double* >( ))] = atype< DbddcR* >( );

  addArray< Dmat >( );

  TheOperators->Add("<-", new PETSc::initCSR< HpSchwarz< PetscScalar > >);
  TheOperators->Add("<-", new PETSc::initCSR< HpSchwarz< PetscScalar > >(1));
  TheOperators->Add("<-", new PETSc::initCSR< HpSchwarz< PetscScalar > >(1, 1));
  TheOperators->Add("<-", new PETSc::initCSR< HpSchwarz< PetscScalar > >(1, 1, 1));
  TheOperators->Add("<-", new PETSc::initCSR< HpSchwarz< PetscScalar > >(1, 1, 1, 1));
  Global.Add("constructor", "(", new PETSc::initCSR< HpSchwarz< PetscScalar >, true >);
  Global.Add("constructor", "(", new PETSc::initCSR< HpSchwarz< PetscScalar >, true >(1));
  Global.Add("constructor", "(", new PETSc::initCSR< HpSchwarz< PetscScalar >, true >(1, 1, 1));
  Global.Add("MatDestroy", "(", new OneOperator1< long, Dmat* >(PETSc::destroyCSR));
  zzzfff->Add("PetscScalar", atype<std::conditional<std::is_same<PetscScalar, PetscReal>::value, double, Complex>::type*>());
  Add< Dmat* >("D", ".", new OneOperator1< KN_<double>, Dmat* >(PETSc::Dmat_D));
  Add< Dmat* >("n", ".", new OneOperator1< long, Dmat* >(PETSc::Dmat_n));
  Add< Dmat* >("range", ".", new OneOperator1s_< KN_<long>, Dmat* >(PETSc::Dmat_range));
#if !defined(PETSC_USE_REAL_SINGLE)
  TheOperators->Add("<-", new PETSc::initCSRfromArray< HpSchwarz< PetscScalar > >);
#endif
  TheOperators->Add("<-", new PETSc::initCSRfromDMatrix< HpSchwarz< PetscScalar >, 0 >);
  TheOperators->Add("<-", new PETSc::initCSRfromDMatrix< HpSchwarz< PetscScalar >, 1 >);
  TheOperators->Add("<-", new PETSc::initCSRfromDMatrix< HpSchwarz< PetscScalar >, 0 >(1));
  TheOperators->Add("<-", new PETSc::initRectangularCSRfromDMatrix< HpSchwarz< PetscScalar >, 0 >);
  TheOperators->Add("<-", new PETSc::initRectangularCSRfromDMatrix< HpSchwarz< PetscScalar >, 1 >);
  TheOperators->Add("<-", new PETSc::initRectangularCSRfromDMatrix< HpSchwarz< PetscScalar >, 0 >(1));
  Global.Add("constructor", "(", new PETSc::initRectangularCSRfromDMatrix< HpSchwarz< PetscScalar >, 0 >(1, 1));
  TheOperators->Add(
    "<-", new OneOperatorCode< PETSc::initCSRfromBlockMatrix< HpSchwarz< PetscScalar > > >( ));
  TheOperators->Add("<-", new OpMatrixtoBilinearFormVGPETSc< HpSchwarz< PetscScalar> >(1));
  TheOperators->Add(
    "=", new OneOperatorCode< PETSc::assignBlockMatrix< HpSchwarz< PetscScalar > > >( ),
         new PETSc::varfToMat< PetscScalar, Mesh, v_fes, v_fes >,
         new PETSc::varfToMat< PetscScalar, Mesh3, v_fes3, v_fes3 >,
         new PETSc::varfToMat< PetscScalar, MeshS, v_fesS, v_fesS >,
         new PETSc::varfToMat< PetscScalar, MeshL, v_fesL, v_fesL >
#if defined(WITH_bemtool) && defined(WITH_htool)
                                                                   ,
         new PETSc::varfToMat< PetscScalar, MeshL, v_fesL, v_fes  >,
         new PETSc::varfToMat< PetscScalar, MeshL, v_fesL, v_fesS >,
         new PETSc::varfToMat< PetscScalar, MeshS, v_fesS, v_fes  >
#endif
                                                                   );
#if defined(WITH_bemtool) && defined(WITH_htool) && defined(PETSC_HAVE_HTOOL)
  typedef const BemKernel fkernel;
  if (!exist_type< fkernel* >( )) map_type[typeid(const BemFormBilinear *).name( )] = new TypeFormBEM;
#endif
  addProd< Dmat, PETSc::ProdPETSc, KN< upscaled_type<PetscScalar> >, PetscScalar, 'N' >( );
  addProd< Dmat, PETSc::ProdPETSc, KN< upscaled_type<PetscScalar> >, PetscScalar, 'T' >( );
  addInv< Dmat, PETSc::InvPETSc, KN< upscaled_type<PetscScalar> >, PetscScalar, 'N' >( );
  addInv< Dmat, PETSc::InvPETSc, KN< upscaled_type<PetscScalar> >, PetscScalar, 'T' >( );
  Global.Add("set", "(", new PETSc::setOptions< Dmat >( ));
  Global.Add("set", "(", new PETSc::setOptions< Dmat >(1));
  Global.Add("set", "(", new PETSc::setOptions< Dmat >(1, 1));
  Global.Add("set", "(", new PETSc::setOptions< Dmat >(1, 1, 1));
  Global.Add("exchange", "(", new exchangeIn< Dmat, PetscScalar >);
  Global.Add("exchange", "(", new exchangeInOut< Dmat, PetscScalar >);
  PETSc::init<PetscReal>();
  Global.Add("ChangeOperator", "(", new PETSc::changeOperator< Dmat >( ));
  Global.Add("ChangeOperator", "(", new PETSc::changeOperator< Dmat >(1));
  TheOperators->Add("=", new OneOperator2_< Dmat*, Dmat*, Matrice_Creuse< upscaled_type<PetscScalar> >* >(PETSc::changeOperatorSimple));
  TheOperators->Add("=", new OneOperator2_< Dmat*, Dmat*, Dmat* >(PETSc::changeOperatorSimple));
  TheOperators->Add("*=", new OneBinaryOperator< PETSc::scale<Dmat*, upscaled_type<PetscScalar>> >);
  Dcl_Type< std::pair<upscaled_type<PetscScalar>, Dmat*>* >();
  Dcl_Type< std::pair<Dmat*, Dmat*>* >();
  atype<std::pair<upscaled_type<PetscScalar>, Dmat*>*>()->AddCast(new E_F1_funcT<std::pair<upscaled_type<PetscScalar>, Dmat*>*, Dmat*>(PETSc::M2P));
  TheOperators->Add("+=", new OneBinaryOperator< PETSc::AXPY<Dmat*, Dmat*> >,
                          new OneOperator2_< Dmat*, Dmat*, std::pair<upscaled_type<PetscScalar>, Dmat*>*, E_F_StackF0F0>(PETSc::AddCombDmat<PetscScalar>));
  TheOperators->Add("*", new OneBinaryOperator<PETSc::Op2<upscaled_type<PetscScalar>>>);
  TheOperators->Add("*", new OneBinaryOperator<PETSc::Op2<Dmat*>>);
  Global.Add("PetscLogStagePush", "(", new OneOperator1_< long, string* >(PETSc::stagePush));
  Global.Add("PetscLogStagePop", "(", new OneOperator0< long >(PETSc::stagePop));
  Global.Add("PetscMemoryGetCurrentUsage", "(", new OneOperator0< double >(PETSc::memoryGetCurrentUsage));
  Global.Add("HasType", "(", new OneOperator2_< long, string*, string* >(PETSc::hasType));
  Global.Add("KSPBuildSolution", "(", new PETSc::buildSolution());
  Init_Common( );
  Dcl_Type< PETSc::DMPlex* >(Initialize< PETSc::DMPlex >, DeleteDTOR< PETSc::DMPlex >);
  zzzfff->Add("DM", atype< PETSc::DMPlex* >( ));
  TheOperators->Add("<-", new PETSc::initDM);
  TheOperators->Add("<-", new PETSc::DMPlexToFF);
  TheOperators->Add("=", new PETSc::DMPlexToFF);
  TheOperators->Add("<-", new PETSc::DMPlexToFF(1));
  TheOperators->Add("=", new PETSc::DMPlexToFF(1));
#ifdef GENERATE_DEPRECATE_FUNCTIONS
  Global.Add("changeOperator", "(", new PETSc::changeOperator< Dmat >( ));
  Global.Add("changeOperator", "(", new PETSc::changeOperator< Dmat >(1));
#endif
}
#ifndef PETScandSLEPc
LOADFUNC(Init_PETSc)
#endif
