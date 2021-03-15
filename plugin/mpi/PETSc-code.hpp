#include "petsc.h"

#include "PETSc.hpp"

typedef PETSc::DistributedCSR< HpSchwarz< PetscScalar > > Dmat;
typedef PETSc::DistributedCSR< HpSchwarz< PetscReal > > DmatR;
typedef PETSc::DistributedCSR< HpSchwarz< PetscComplex > > DmatC;
typedef PETSc::DistributedCSR< HpSchur< PetscScalar > > Dbddc;
typedef PETSc::DistributedCSR< HpSchur< PetscReal > > DbddcR;
typedef PETSc::DistributedCSR< HpSchur< PetscComplex > > DbddcC;

namespace PETSc {
  template<class K, class MMesh, class fes1, class fes2 >
  struct varfToMat : public OneOperator {
    class Op : public E_F0mps {
    public:
      Call_FormBilinear<fes1,fes2>* b;
      Expression a;
      AnyType operator()(Stack s) const;

      Op(Expression x, Expression  y) : b(new Call_FormBilinear<fes1,fes2>(*dynamic_cast<const Call_FormBilinear<fes1,fes2>*>(y))), a(x) {
          assert(b && b->nargs);
          ffassert(FieldOfForm(b->largs, IsComplexType<HPDDM::upscaled_type<K>>::value) == IsComplexType<HPDDM::upscaled_type<K>>::value);
      }
      operator aType () const { return atype<Dmat*>(); }
    };
    E_F0* code(const basicAC_F0& args) const {
        return new Op(to<Dmat*>(args[0]), args[1]);
    }
    varfToMat() : OneOperator(atype<Dmat*>(), atype<Dmat*>(), atype<const Call_FormBilinear<fes1,fes2>*>()) {}
  };
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
    int np = OpCall_FormBilinear_np::n_name_param - 6;
    SetEnd_Data_Sparse_Solver<HPDDM::upscaled_type<K>>(stack, ds, b->nargs, np);

    WhereStackOfPtr2Free(stack) = new StackOfPtr2Free(stack);

    Dmat& B(*GetAny<Dmat*>((*a)(stack)));
    Matrice_Creuse<HPDDM::upscaled_type<K>> A;
    A.init();
    if(!PUh || !PVh)
      return SetAny<Dmat*>(&B);
    const FESpace1& Uh = *PUh;
    const FESpace2& Vh = *PVh;
    const Mesh1& Th = Uh.Th;
    bool same = isSameMesh(b->largs, &Uh.Th, &Vh.Th, stack);
    if(same) {
      if(A.Uh != Uh || A.Vh != Vh) {
        A.Uh = Uh;
        A.Vh = Vh;
        if(ds.sym) {
          A.A.master(new MatriceMorse<HPDDM::upscaled_type<K>>(Vh.NbOfDF, Vh.NbOfDF, 0, ds.sym));
          ffassert(&Uh == &Vh);
        }
        else
          A.A.master(new MatriceMorse<HPDDM::upscaled_type<K>>(Vh.NbOfDF, Uh.NbOfDF, 2 * Vh.NbOfDF, 0));
      }
      if(AssembleVarForm<HPDDM::upscaled_type<K>, MatriceCreuse<HPDDM::upscaled_type<K>>, MMesh, FESpace1,FESpace2>(stack, Th, Uh, Vh, ds.sym, A.A, 0, b->largs))
        AssembleBC<HPDDM::upscaled_type<K>, MMesh,FESpace1, FESpace2>(stack, Th, Uh, Vh, ds.sym, A.A, 0, 0, b->largs, ds.tgv);
    }
    else {
      MatriceMorse<HPDDM::upscaled_type<K>> *pMA = new MatriceMorse<HPDDM::upscaled_type<K>>(Vh.NbOfDF, Uh.NbOfDF, 0, ds.sym);
      MatriceMap<HPDDM::upscaled_type<K>>& D = *pMA;
      bool bc = AssembleVarForm<HPDDM::upscaled_type<K>, MatriceMap<HPDDM::upscaled_type<K>>, MMesh, FESpace1,FESpace2>(stack, Th, Uh, Vh, ds.sym, &D, 0, b->largs);
      A.A.master(pMA);
      if(bc)
        AssembleBC<HPDDM::upscaled_type<K>>(stack, Th, Uh, Vh, ds.sym, A.A, 0, 0, b->largs, ds.tgv);
    }
    changeOperatorSimple(&B, &A);
    B._A->setMatrix(nullptr);
    return SetAny<Dmat*>(&B);
  }

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
    int size;
    MPI_Comm_size(ptA->_A->getCommunicator(), &size);
    if (ptD || size == 1) {
      if (ptD) {
        if(!ptA->_D) {
          PetscReal* d = reinterpret_cast<PetscReal*>(ptD->operator double*());
          if(!std::is_same<HPDDM::upscaled_type<PetscReal>, PetscReal>::value) {
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
      ptA->_A->distributedNumbering(ptA->_num, ptA->_first, ptA->_last, global);
      if (verbosity > 0 && mpirank == 0)
        cout << " --- global numbering created (in " << MPI_Wtime( ) - timing << ")" << endl;
    } else
      global = PETSC_DECIDE;
    timing = MPI_Wtime( );
    PetscInt* ia = nullptr;
    PetscInt* ja = nullptr;
    PetscScalar* c = nullptr;
    bool free = ptA->_A->getMatrix( )->_ia
                  ? ptA->_A->distributedCSR(ptA->_num, ptA->_first, ptA->_last, ia, ja, c)
                  : false;
    MatCreate(ptA->_A->getCommunicator(), &ptA->_petsc);
    if (bs > 1) MatSetBlockSize(ptA->_petsc, bs);
    MatSetSizes(ptA->_petsc, ptA->_last - ptA->_first, ptA->_last - ptA->_first, global, global);
    bool sym = ptA->_A->getMatrix( )->_sym;
    if (ia) {
      if (sym) {
        MatSetType(ptA->_petsc, MATMPISBAIJ);
        MatMPISBAIJSetPreallocationCSR(ptA->_petsc, 1, ia, ja, c);
      } else {
        MatSetType(ptA->_petsc, MATMPIAIJ);
        MatMPIAIJSetPreallocationCSR(ptA->_petsc, ia, ja, c);
        MatSetOption(ptA->_petsc, MAT_SYMMETRIC, symmetric);
      }
    } else {
      MatSetType(ptA->_petsc, MATMPIAIJ);
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
    if (!M->_sym) cout << "Please assemble a symmetric CSR" << endl;
    double timing = MPI_Wtime( );
    ptA->_A->template renumber< false >(STL< long >(*ptD), nullptr);
    long long global;
    ptA->_A->distributedNumbering(ptA->_num, ptA->_first, ptA->_last, global);
    if (verbosity > 0 && mpirank == 0)
      cout << " --- global numbering created (in " << MPI_Wtime( ) - timing << ")" << endl;
    timing = MPI_Wtime( );
    PetscInt* indices;
    PetscMalloc(sizeof(PetscInt) * M->_n / bs, &indices);
    for (unsigned int i = 0; i < M->_n; i += bs) indices[i / bs] = ptA->_num[i] / bs;
    ISLocalToGlobalMapping rmap;
    ISLocalToGlobalMappingCreate(ptA->_A->getCommunicator(), bs, M->_n / bs, indices, PETSC_OWN_POINTER,
                                 &rmap);
    MatCreateIS(ptA->_A->getCommunicator(), bs, PETSC_DECIDE, PETSC_DECIDE, global, global, rmap, NULL,
                &ptA->_petsc);
    Mat local;
    MatISGetLocalMat(ptA->_petsc, &local);
    MatSetType(local, MATSEQSBAIJ);
    std::vector< std::vector< std::pair< int, PetscScalar > > > transpose(M->_n);
    for (int i = 0; i < transpose.size( ); ++i)
      for (int j = M->_ia[i]; j < M->_ia[i + 1]; ++j) {
        transpose[M->_ja[j]].emplace_back(i, M->_a[j]);
        if (bs > 1 && (i - M->_ja[j] <= (i % bs)) && M->_ja[j] != i)
          transpose[i].emplace_back(M->_ja[j], M->_a[j]);
      }
    int nnz = 0;
    for (int i = 0; i < transpose.size( ); ++i) {
      std::sort(transpose[i].begin( ), transpose[i].end( ),
                [](const std::pair< int, PetscScalar >& lhs,
                   const std::pair< int, PetscScalar >& rhs) { return lhs.first < rhs.first; });
      nnz += transpose[i].size( );
    }
    PetscInt* ia = new PetscInt[M->_n / bs + 1];
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
    VecSetType(isVec, VECMPI);
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
    if (A->_sym) {
      std::vector< std::pair< int, int > >* transpose =
        new std::vector< std::pair< int, int > >[A->_n]( );
      for (int i = 0; i < A->_n; ++i)
        for (int j = A->_ia[i] - (HPDDM_NUMBERING == 'F');
             j < A->_ia[i + 1] - (HPDDM_NUMBERING == 'F'); ++j)
          transpose[A->_ja[j] - (HPDDM_NUMBERING == 'F')].emplace_back(i, j);
      for (int i = 0; i < A->_n; ++i)
        std::sort(transpose[i].begin( ), transpose[i].end( ));
      PetscInt* ia = new PetscInt[A->_n + 1];
      PetscInt* ja = new PetscInt[A->_nnz];
      PetscScalar* c = new PetscScalar[A->_nnz];
      ia[0] = 0;
      for (int i = 0; i < A->_n; ++i) {
        for (int j = 0; j < transpose[i].size( ); ++j) {
          c[ia[i] + j] = A->_a[transpose[i][j].second];
          ja[ia[i] + j] = transpose[i][j].first;
        }
        ia[i + 1] = ia[i] + transpose[i].size( );
      }
      delete[] transpose;
      MatCreate(PETSC_COMM_SELF, &aux);
      MatSetSizes(aux, A->_n, A->_n, A->_n, A->_n);
      MatSetType(aux, MATSEQSBAIJ);
      MatSeqSBAIJSetPreallocationCSR(aux, 1, ia, ja, c);
      delete[] c;
      delete[] ja;
      delete[] ia;
    } else if(std::is_same<decltype(A->_ia), PetscInt>::value)
      MatCreateSeqAIJWithArrays(PETSC_COMM_SELF, A->_n, A->_m, reinterpret_cast<PetscInt*>(A->_ia), reinterpret_cast<PetscInt*>(A->_ja), A->_a, &aux);
    else {
      PetscInt* ia = new PetscInt[A->_n + 1];
      std::copy_n(A->_ia, A->_n + 1, ia);
      PetscInt* ja = new PetscInt[A->_nnz];
      std::copy_n(A->_ja, A->_nnz, ja);
      PetscScalar* c = new PetscScalar[A->_nnz];
      std::copy_n(A->_a, A->_nnz, c);
      MatCreate(PETSC_COMM_SELF, &aux);
      MatSetSizes(aux, A->_n, A->_n, A->_n, A->_n);
      MatSetType(aux, MATSEQAIJ);
      MatSeqAIJSetPreallocationCSR(aux, ia, ja, c);
      delete[] c;
      delete[] ja;
      delete[] ia;
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
          B = to< Matrice_Creuse< HPDDM::upscaled_type<PetscScalar> >* >(args[1]);
        else
          B = to< Type* >(args[1]);
      }

      AnyType operator( )(Stack stack) const;
      operator aType( ) const { return atype< long >( ); }
    };
    E_F0* code(const basicAC_F0& args) const { return new changeOperator_Op(args, c); }
    changeOperator( )
      : OneOperator(atype< long >( ), atype< Type* >( ),
                    atype< Matrice_Creuse< HPDDM::upscaled_type<PetscScalar> >* >( )),
        c(0) {}
    changeOperator(int)
      : OneOperator(atype< long >( ), atype< Type* >( ), atype< Type* >( )), c(1) {}
  };
  template< class Type >
  basicAC_F0::name_and_type changeOperator< Type >::changeOperator_Op::name_param[] = {
    {"restriction", &typeid(Matrice_Creuse< double >*)}, {"parent", &typeid(Type*)}};
  template< class Type >
  void change(Type* const& ptA, Matrice_Creuse< HPDDM::upscaled_type<PetscScalar> >* const& mat, Type* const& ptB,
              Matrice_Creuse< double >* const& pList, Type* const& ptParent) {
    if (mat) {
      if (ptA) {
        MatriceMorse< HPDDM::upscaled_type<PetscScalar> >* mN = nullptr;
        if (mat->A) mN = static_cast< MatriceMorse< HPDDM::upscaled_type<PetscScalar> >* >(&*(mat->A));
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
            unsigned int* perm = new unsigned int[dM->_n]( );
            for (unsigned int i = 0; i < dL->_n; ++i) perm[dL->_ja[i]] = i + 1;
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
                          PetscMalloc1(dN->_n, &idx);
                          std::copy_n(ptA->_num, dN->_n, idx);
                          IS is;
                          ISCreateGeneral(PETSC_COMM_SELF, ptA->_A->getMatrix()->_n, idx, PETSC_OWN_POINTER, &is);
                          PetscObjectCompose((PetscObject)pc, "_PCHPDDM_Neumann_IS", (PetscObject)is);
                          PetscObjectCompose((PetscObject)pc, "_PCHPDDM_Neumann_Mat", (PetscObject)aux);
                          ISDestroy(&is);
                          MatDestroy(&aux);
                      }
                      else {
                          MatHeaderReplace(N, &aux);
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
          if (ptA->_cnum)
            free = HPDDM::template Subdomain< PetscScalar >::distributedCSR(
              ptA->_num, ptA->_first, ptA->_last, ia, ja, c, dN, ptA->_num + dN->_n);
          else
            free = ptA->_A->distributedCSR(ptA->_num, ptA->_first, ptA->_last, ia, ja, c);
          if (assembled) {
            MatZeroEntries(ptA->_petsc);
            for (PetscInt i = 0; i < ptA->_last - ptA->_first; ++i) {
              PetscInt row = ptA->_first + i;
              MatSetValues(ptA->_petsc, 1, &row, ia[i + 1] - ia[i],
                           ja + ia[i], c + ia[i], INSERT_VALUES);
            }
          } else {
            if (!ptA->_petsc) {
              MatCreate(PETSC_COMM_WORLD, &ptA->_petsc);
              MatSetSizes(ptA->_petsc, ptA->_last - ptA->_first, ptA->_last - ptA->_first,
                          PETSC_DECIDE, PETSC_DECIDE);
            }
            if (ptA->_A && ptA->_A->getMatrix( )->_sym) {
              MatSetType(ptA->_petsc, MATMPISBAIJ);
              MatMPISBAIJSetPreallocationCSR(ptA->_petsc, 1, ia, ja, c);
            } else {
              MatSetType(ptA->_petsc, MATMPIAIJ);
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
          MatMPIAIJSetPreallocationCSR(ptA->_petsc, ia, NULL, NULL);
          delete[] ia;
        }
        MatAssemblyBegin(ptA->_petsc, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(ptA->_petsc, MAT_FINAL_ASSEMBLY);
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
            if (assembled) MatConvert(ptB->_petsc, MATSAME, MAT_INITIAL_MATRIX, &ptA->_petsc);
            MatDestroy(&ptB->_petsc);
          }
        } else {
          MatHeaderReplace(ptA->_petsc, &ptB->_petsc);
        }
        if (ptB->_ksp) KSPDestroy(&ptB->_ksp);
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
    Matrice_Creuse< HPDDM::upscaled_type<PetscScalar> >* mat =
      c == 0 ? GetAny< Matrice_Creuse< HPDDM::upscaled_type<PetscScalar> >* >((*B)(stack)) : nullptr;
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
  Dmat* changeOperatorSimple(Dmat* const& dA, Matrice_Creuse< HPDDM::upscaled_type<PetscScalar> >* const& A) {
    Dmat* const null = nullptr;
    change(dA, A, null, nullptr, null);
    return dA;
  }
  long changeSchur(Dmat* const& dA, KN< Matrice_Creuse< HPDDM::upscaled_type<PetscScalar> > >* const& schurPreconditioner,
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

    IS from, to;
    PetscInt* idx = new PetscInt[in->n];
    for (int i = 0; i < in->n; ++i) idx[i] = numbering->operator[](i);
    ISCreateGeneral(PETSC_COMM_SELF, in->n, idx, PETSC_USE_POINTER, &to);
    ISCreateStride(PETSC_COMM_SELF, in->n, low, 1, &from);
    VecScatter scatter;
    VecScatterCreate(x, NULL, y, to, &scatter);
    VecScatterBegin(scatter, x, y, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(scatter, x, y, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterDestroy(&scatter);
    ISDestroy(&to);
    ISDestroy(&from);
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
        if (c == 0) K = to< Matrice_Creuse< PetscScalar >* >(args[2]);
        else if (c == 1) {
          const Polymorphic* op = dynamic_cast< const Polymorphic* >(args[2].LeftValue( ));
          ffassert(op);
          codeA = op->Find("(", ArrayOfaType(atype< KN< PetscScalar >* >( ), false));
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
          atype< DistributedCSR< HpddmType >* >( ), atype< typename std::conditional<C == 0, Matrice_Creuse< PetscScalar >*, Polymorphic*>::type >( )),
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
        if (c == 0 || c == 3) K = to< Matrice_Creuse< PetscScalar >* >(args[3]);
        if (c == 2) {
          const Polymorphic* op = dynamic_cast< const Polymorphic* >(args[3].LeftValue( ));
          ffassert(op);
          codeA = op->Find("(", ArrayOfaType(atype< KN< PetscScalar >* >( ), false));
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
          atype< typename std::conditional<D == 0, Matrice_Creuse< PetscScalar >*, Polymorphic* >::type >( )),
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
          atype< Matrice_Creuse< PetscScalar >* >( )),
        c(3) {}
  };
  template< class HpddmType, int D >
  basicAC_F0::name_and_type
    initRectangularCSRfromDMatrix< HpddmType, D >::initRectangularCSRfromDMatrix_Op::name_param[] = {
      {D == 1 ? "transpose" : "clean", D == 1 ? &typeid(Polymorphic*) : &typeid(bool)},
      {"numbering", &typeid(KN<PetscScalar>*)}};

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
    MatSetType(ptA->_petsc, bsr ? MATMPIBAIJ : MATMPIAIJ);
    mK->CSR( );
    if (bsr)
      MatMPIBAIJSetPreallocationCSR(ptA->_petsc, bs, reinterpret_cast< PetscInt* >(mK->p),
                                    reinterpret_cast< PetscInt* >(mK->j), mK->aij);
    else
      MatMPIAIJSetPreallocationCSR(ptA->_petsc, reinterpret_cast< PetscInt* >(mK->p),
                                   reinterpret_cast< PetscInt* >(mK->j), mK->aij);
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
      MatSetType(ptA->_petsc, MATMPIAIJ);
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
        if (c == 1 || c == 3)
          K = to< long >(args[1]);
        else
          K = to< Matrice_Creuse< HPDDM::upscaled_type<PetscScalar> >* >(args[1]);
        if (c == 0 || c == 1) {
          R = to< KN< KN< long > >* >(args[2]);
          D = to< KN< typename std::conditional<
            std::is_same< HpddmType, HpSchwarz< PetscScalar > >::value, double, long >::type >* >(
            args[3]);
        }
      }

      AnyType operator( )(Stack stack) const;
      operator aType( ) const { return atype< DistributedCSR< HpddmType >* >( ); }
    };
    E_F0* code(const basicAC_F0& args) const { return new E_initCSR(args, c); }
    initCSR( )
      : OneOperator(
          atype< DistributedCSR< HpddmType >* >( ), atype< DistributedCSR< HpddmType >* >( ),
          atype< Matrice_Creuse< HPDDM::upscaled_type<PetscScalar> >* >( ), atype< KN< KN< long > >* >( ),
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
                    atype< Matrice_Creuse< HPDDM::upscaled_type<PetscScalar> >* >( )),
        c(2) {}
    initCSR(int, int, int)
      : OneOperator(atype< DistributedCSR< HpddmType >* >( ),
                    atype< DistributedCSR< HpddmType >* >( ), atype< long >( )),
        c(3) {}
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
    MatriceMorse< HPDDM::upscaled_type<PetscScalar> >* mA = nullptr;
    if (c == 0 || c == 2)
      mA = static_cast< MatriceMorse< HPDDM::upscaled_type<PetscScalar> >* >(
        &(*GetAny< Matrice_Creuse< HPDDM::upscaled_type<PetscScalar> >* >((*K)(stack))->A));
    else
      dof = GetAny< long >((*K)(stack));
    MPI_Comm* comm = nargs[0] ? (MPI_Comm*)GetAny< pcommworld >((*nargs[0])(stack)) : 0;
    if(c == 2 || c == 3) {
        static MPI_Comm self = MPI_COMM_SELF;
        MPI_Comm& rval = self;
        comm = &rval;
    }
    if (c == 2 && mA->n != mA->m) {
      ptA->_first = 0;
      ptA->_last = mA->n;
      ptA->_cfirst = 0;
      ptA->_clast = mA->m;
      MatCreate(*comm, &ptA->_petsc);
      MatSetSizes(ptA->_petsc, mA->n, mA->m, PETSC_DECIDE, PETSC_DECIDE);
      MatSetType(ptA->_petsc, MATMPIAIJ);
      ff_HPDDM_MatrixCSR< PetscScalar > dA(mA);
      ptA->_num = new PetscInt[mA->n + mA->m];
      ptA->_cnum = ptA->_num + mA->n;
      std::iota(ptA->_num, ptA->_num + mA->n, 0);
      std::iota(ptA->_cnum, ptA->_cnum + mA->m, 0);
      MatMPIAIJSetPreallocationCSR(ptA->_petsc, dA._ia, dA._ja, dA._a);
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
        ptA->_num = new PetscInt[ptA->_A->getMatrix( )->_n];
        initPETScStructure<C>(ptA, bs,
          nargs[3] && GetAny< bool >((*nargs[3])(stack)) ? PETSC_TRUE : PETSC_FALSE, empty ? empty : ptD);
        delete empty;
        if (!std::is_same< HpddmType, HpSchwarz< PetscScalar > >::value) {
          mA->p = ptA->_A->getMatrix( )->_ia;
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
        N = 0;
        M = 0;
        e_M = nullptr;
        t_M = nullptr;
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
      for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M; ++j) {
          Expression e = e_M[i][j];
          int t = t_M[i][j];
          if (e) {
            AnyType e_ij = (*e)(s);
            if (t == 1) {
              DistributedCSR< HpddmType >* pt = GetAny< DistributedCSR< HpddmType >* >(e_ij);
              a[i * M + j] = pt->_petsc;
              exchange[i * M + j] = pt;
            } else if (t == 2) {
              DistributedCSR< HpddmType >* pt = GetAny< DistributedCSR< HpddmType >* >(e_ij);
              Mat B;
              if (std::is_same< PetscScalar, PetscReal >::value)
                MatCreateTranspose(pt->_petsc, &B);
              else
                MatCreateHermitianTranspose(pt->_petsc, &B);
              a[i * M + j] = B;
              destroy.emplace_back(i, j);
              exchange[i * M + j] = pt;
            } else if (t == 7) {
              PetscScalar r = GetAny< PetscScalar >(e_ij);
              if (std::abs(r) > 1.0e-16) {
                Mat B;
                MatCreate(PETSC_COMM_WORLD, &B);
                MatSetSizes(B, PETSC_DECIDE, PETSC_DECIDE, 1, 1);
                MatSetType(B, MATMPIDENSE);
                MatMPIDenseSetPreallocation(B, PETSC_NULL);
                PetscScalar* array;
                MatDenseGetArray(B, &array);
                if (array) array[0] = r;
                MatDenseRestoreArray(B, &array);
                MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY);
                MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY);
                a[i * M + j] = B;
                destroy.emplace_back(i, j);
              } else if (i == j)
                zeros.emplace_back(i);
            } else if (t == 3 || t == 4) {
              KN< PetscScalar > x;
              Mat B;
              MatCreate(PETSC_COMM_WORLD, &B);
              if (t == 3) x = GetAny< KN_< PetscScalar > >(e_ij);
              else x = GetAny< Transpose< KN_< PetscScalar > > >(e_ij);
              MatSetSizes(B, x.n, PETSC_DECIDE, PETSC_DECIDE, 1);
              MatSetType(B, MATMPIDENSE);
              MatMPIDenseSetPreallocation(B, PETSC_NULL);
              PetscScalar* array;
              MatDenseGetArray(B, &array);
              if (array) std::copy_n(static_cast< PetscScalar* >(x), x.n, array);
              MatDenseRestoreArray(B, &array);
              MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY);
              MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY);
              if (t == 3) {
                a[i * M + j] = B;
              }
              else {
                Mat C;
                if (std::is_same< PetscScalar, PetscReal >::value) MatCreateTranspose(B, &C);
                else MatCreateHermitianTranspose(B, &C);
                MatDestroy(&B);
                a[i * M + j] = C;
              }
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
          MatCreate(PETSC_COMM_WORLD, a + zeros[i] * M + zeros[i]);
          MatSetSizes(a[zeros[i] * M + zeros[i]], x, y, X, Y);
          MatSetType(a[zeros[i] * M + zeros[i]], MATMPIAIJ);
          MatMPIAIJSetPreallocation(a[zeros[i] * M + zeros[i]], 0, NULL, 0, NULL);
          MatAssemblyBegin(a[zeros[i] * M + zeros[i]], MAT_FINAL_ASSEMBLY);
          MatAssemblyEnd(a[zeros[i] * M + zeros[i]], MAT_FINAL_ASSEMBLY);
          destroy.emplace_back(zeros[i], zeros[i]);
        }
      }
      Result sparse_mat = GetAny< Result >((*emat)(s));
      if (sparse_mat->_petsc) sparse_mat->dtor( );
      MatCreateNest(PETSC_COMM_WORLD, N, NULL, M, NULL, a, &sparse_mat->_petsc);
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
      static const int n_name_param = 19;
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
    {"nearnullspace", &typeid(FEbaseArrayKn< HPDDM::upscaled_type<PetscScalar> >*)},              // 1
    {"fields", &typeid(KN< double >*)},                                                           // 2
    {"names", &typeid(KN< String >*)},                                                            // 3
    {"prefix", &typeid(std::string*)},                                                            // 4
    {"schurPreconditioner", &typeid(KN< Matrice_Creuse< HPDDM::upscaled_type<PetscScalar> > >*)}, // 5
    {"schurList", &typeid(KN< double >*)},                                                        // 6
    {"parent", &typeid(Type*)},                                                                   // 7
    {"MatNullSpace", &typeid(KNM< HPDDM::upscaled_type<PetscScalar> >*)},                         // 8
    {"fieldsplit", &typeid(long)},                                                                // 9
    {"schurComplement", &typeid(KNM< HPDDM::upscaled_type<PetscScalar> >*)},                      // 10
    {"schur", &typeid(KN< Dmat >*)},                                                              // 11
    {"aux", &typeid(Matrice_Creuse< HPDDM::upscaled_type<PetscScalar> >*)},                       // 12
    {"coordinates", &typeid(KNM< double >*)},                                                     // 13
    {"gradient", &typeid(Dmat*)},                                                                 // 14
    {"O", &typeid(Matrice_Creuse< HPDDM::upscaled_type<PetscScalar> >*)},                         // 15
    {"bs", &typeid(long)},                                                                        // 16
    {"precon", &typeid(Polymorphic*)},                                                            // 17
    {"setup", &typeid(bool)}                                                                      // 18
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

    PetscFunctionBegin;
    if (scall == MAT_INITIAL_MATRIX) {
      ierr = PetscCalloc1(1,submat);CHKERRQ(ierr);
    }
    (*submat)[0] = *O;
    PetscFunctionReturn(0);
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
  template< bool T, class K, class U, typename std::enable_if< std::is_same< K, HPDDM::upscaled_type<PetscScalar> >::value && std::is_same< U, HPDDM::upscaled_type<U>>::value >::type* = nullptr >
  void MatMult(MatriceMorse<K>* const& A, KN_<U>& in, KN_<U>& out) {
    A->addMatMul(in, out, T, in.step, out.step);
  }
  template< bool T, class K, class U, typename std::enable_if< !std::is_same< K, HPDDM::upscaled_type<PetscScalar> >::value && std::is_same< U, HPDDM::upscaled_type<U>>::value >::type* = nullptr >
  void MatMult(MatriceMorse<K>* const& A, KN_<U>& in, KN_<U>& out) {
    HPDDM::upscaled_type<PetscScalar>* pc = in;
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
  template< bool T, class K, class U, typename std::enable_if< !std::is_same< U, HPDDM::upscaled_type<U>>::value >::type* = nullptr >
  void MatMult(MatriceMorse<K>* const& A, KN_<U>& in, KN_<U>& out) {
    KN< HPDDM::upscaled_type<PetscScalar> > inUp(in.N());
    KN< HPDDM::upscaled_type<PetscScalar> > outUp(out.N());
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

    PetscFunctionBegin;
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
    PetscFunctionReturn(0);
  }
  template< class Type, typename std::enable_if<!std::is_same<Type, ShellInjection>::value>::type* = nullptr  >
  static PetscErrorCode ShellDestroy(Mat A) {
    User< Type > user;
    PetscErrorCode ierr;

    PetscFunctionBegin;
    ierr = MatShellGetContext(A, &user);CHKERRQ(ierr);
    delete user->op;
    ierr = PetscFree(user);CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
  template< class Type, typename std::enable_if<std::is_same<Type, ShellInjection>::value>::type* = nullptr  >
  static PetscErrorCode ShellDestroy(Mat A) {
    User< Type > user;
    PetscErrorCode ierr;

    PetscFunctionBegin;
    ierr = MatShellGetContext(A, &user);CHKERRQ(ierr);
    ierr = PetscFree(user);CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
  template< class Type >
  static PetscErrorCode PCShellDestroy(PC pc) {
    User< Type > user;
    PetscErrorCode ierr;

    PetscFunctionBegin;
    ierr = PCShellGetContext(pc, (void**)&user);CHKERRQ(ierr);
    delete user->op;
    ierr = PetscFree(user);CHKERRQ(ierr);
    PetscFunctionReturn(0);
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
    PetscFunctionReturn(0);
  }
  static PetscErrorCode WrapperDestroy(PC pc) {
    WrapperSubKSP *ctx;
    PetscFunctionBeginUser;
    PCShellGetContext(pc, (void**)&ctx);
    ISDestroy(&ctx->is);
    PetscFree(ctx);
    PetscFunctionReturn(0);
  }
  static PetscErrorCode WrapperView(PC pc, PetscViewer viewer) {
    WrapperSubKSP *ctx;
    PetscFunctionBeginUser;
    PCShellGetContext(pc, (void**)&ctx);
    KSPView(*ctx->ksp, viewer);
    PetscFunctionReturn(0);
  }
  template< class, class Container, char N = 'N' >
  static PetscErrorCode Op_User(Container, Vec, Vec);
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
                MatConvert(ptA->_petsc, MATMPISBAIJ, MAT_INPLACE_MATRIX, &ptA->_petsc);
            else
                MatConvert(ptA->_petsc, MATMPIAIJ, MAT_INPLACE_MATRIX, &ptA->_petsc);
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
      if (c == 1 || c == 3) {
        if (nargs[4])
            KSPSetOptionsPrefix(ptA->_ksp, GetAny< std::string* >((*nargs[4])(stack))->c_str( ));
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
          KN< Matrice_Creuse< HPDDM::upscaled_type<PetscScalar> > >* mS =
            nargs[5] ? GetAny< KN< Matrice_Creuse< HPDDM::upscaled_type<PetscScalar> > >* >((*nargs[5])(stack)) : 0;
          KN< double >* pL = nargs[6] ? GetAny< KN< double >* >((*nargs[6])(stack)) : 0;
          KN< Dmat >* mdS = nargs[11] ? GetAny< KN< Dmat >* >((*nargs[11])(stack)) : 0;
          if (mdS)
            setFieldSplitPC(ptA, ksp, fields, names, mdS);
          else
            setFieldSplitPC(ptA, ksp, fields, names, mS, pL);
        }
      }
      else ksp = ptA->_ksp;
      if (nargs[4] && c != 1 && c != 2 || c == 3)
        KSPSetOptionsPrefix(ptA->_ksp, GetAny< std::string* >((*nargs[4])(stack))->c_str( ));
      KSPSetFromOptions(ksp);
      if (c != 1) {
        if (std::is_same< Type, Dmat >::value) {
          FEbaseArrayKn< HPDDM::upscaled_type<PetscScalar> >* ptNS =
            nargs[1] ? GetAny< FEbaseArrayKn< HPDDM::upscaled_type<PetscScalar> >* >((*nargs[1])(stack)) : 0;
          KNM< HPDDM::upscaled_type<PetscScalar> >* ptPETScNS =
            nargs[8] ? GetAny< KNM< HPDDM::upscaled_type<PetscScalar> >* >((*nargs[8])(stack)) : 0;
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
                HPDDM::upscaled_type<PetscScalar>* get = *ptNS->get(i);
                PetscScalar* base = reinterpret_cast<PetscScalar*>(get);
                if(!std::is_same<HPDDM::upscaled_type<PetscReal>, PetscReal>::value) {
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
          KNM< HPDDM::upscaled_type<PetscReal> >* coordinates =
            nargs[13] ? GetAny< KNM< HPDDM::upscaled_type<PetscReal> >* >((*nargs[13])(stack)) : nullptr;
          Dmat* G = nargs[14] ? GetAny< Dmat* >((*nargs[14])(stack)) : nullptr;
          if (coordinates && std::is_same<PetscReal, HPDDM::upscaled_type<PetscReal>>::value) PCSetCoordinates(pc, coordinates->N( ), coordinates->M( ), reinterpret_cast<PetscReal*>(coordinates->operator double*()));
#ifdef PETSC_HAVE_HYPRE
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
#if defined(PCHPDDM) && defined(PETSC_USE_SHARED_LIBRARIES)
            PetscStrcmp(type, PCHPDDM, &isType);
            const HPDDM::MatrixCSR< PetscScalar >* const A = ptA->_A->getMatrix( );
            if (isType && A) {
              Mat aux = ff_to_PETSc(A);
              PCHPDDMSetAuxiliaryMat(pc, is, aux, NULL, NULL);
              MatDestroy(&aux);
              Matrice_Creuse< HPDDM::upscaled_type<PetscScalar> >* ptK =
                nargs[12] ? GetAny< Matrice_Creuse< HPDDM::upscaled_type<PetscScalar> >* >((*nargs[12])(stack)) : nullptr;
#if PETSC_VERSION_GE(3, 13, 0)
              if (ptK && ptK->A) {
                MatriceMorse< HPDDM::upscaled_type<PetscScalar> >* mA =
                  static_cast< MatriceMorse< HPDDM::upscaled_type<PetscScalar> >* >(&(*ptK->A));
                HPDDM::MatrixCSR< PetscScalar >* B = new_HPDDM_MatrixCSR< PetscScalar >(mA);
                aux = ff_to_PETSc(B);
                PCHPDDMSetRHSMat(pc, aux);
                MatDestroy(&aux);
                delete B;
              }
#endif
            }
#endif
            Matrice_Creuse< HPDDM::upscaled_type<PetscScalar> >* ptO =
              nargs[15] ? GetAny< Matrice_Creuse< HPDDM::upscaled_type<PetscScalar> >* >((*nargs[15])(stack)) : nullptr;
            if (ptO && ptO->A) {
              MatriceMorse< HPDDM::upscaled_type<PetscScalar> >* mO =
                static_cast< MatriceMorse< HPDDM::upscaled_type<PetscScalar> >* >(&(*ptO->A));
              ff_HPDDM_MatrixCSR< PetscScalar > dO(mO);
              PCSetType(pc, PCASM);
              IS loc;
              PetscInt n, first;
              MatGetOwnershipRange(ptA->_petsc, &first, &n);
              n -= first;
              ISCreateStride(PETSC_COMM_SELF, n, first, 1, &loc);
              PCASMSetLocalSubdomains(pc, 1, &is, &loc);
              int nnz = dO._nnz;
              MPI_Allreduce(MPI_IN_PLACE, &nnz, 1, MPI_INT, MPI_MAX, PetscObjectComm((PetscObject)ptA->_petsc));
              PetscErrorCode (*CreateSubMatrices)(Mat,PetscInt,const IS*,const IS*,MatReuse,Mat**);
              if (nnz) {
                MatGetOperation(ptA->_petsc, MATOP_CREATE_SUBMATRICES, (void(**)(void))&CreateSubMatrices);
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
                MatSetOperation(ptA->_petsc, MATOP_CREATE_SUBMATRICES, (void(*)(void))CreateSubMatrices);
                delete O;
              }
              O = nullptr;
            }
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
            const Polymorphic* op = dynamic_cast< const Polymorphic* >(nargs[17]);
            ffassert(op);
            const OneOperator* codeA = op->Find("(", ArrayOfaType(atype< KN< PetscScalar >* >( ), false));
            PCSetType(pc, PCSHELL);
            User< LinearSolver< Type > > userPC = nullptr;
            PetscNew(&userPC);
            PetscInt n;
            MatGetLocalSize(ptA->_petsc, &n, NULL);
            userPC->op = new typename LinearSolver< Type >::MatF_O(n, stack, codeA);
            PCShellSetContext(pc, userPC);
            PCShellSetApply(pc, Op_User< LinearSolver< Type >, PC >);
            PCShellSetDestroy(pc, PCShellDestroy< LinearSolver< Dmat >  >);
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
          } else {
            E.reserve(1);
            E.emplace_back(to< Type* >(args[0]), to< Storage< PetscScalar >* >(args[1]));
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
    if (c == 0 || c == 2 || c == 3) {
      PetscScalar* pt = *ptOut;
      for (int j = 0; j < E.size( ); ++j) {
        Type* ptA = GetAny< Type* >((*(E[j].first))(stack));
        if (ptA && (ptA->_last - ptA->_first || (inverse && ptA->_num))) {
          PetscInt m;
          MatGetLocalSize(ptA->_petsc, &m, NULL);
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
  template< char P, class Type >
  long MatMatMult(Type* const& A, KNM< PetscScalar >* const& in, KNM< PetscScalar >* const& out) {
    if (A->_petsc) {
      Mat C;
      MatType type;
      MatGetType(A->_petsc, &type);
      PetscBool isType;
      PetscStrcmp(type, MATMPIAIJ, &isType);
      int size;
      MPI_Comm_size(PetscObjectComm((PetscObject)A->_petsc), &size);
      if (size == 1 && isType)
        MatMPIAIJGetLocalMat(A->_petsc, MAT_INITIAL_MATRIX, &C);
      else
        C = A->_petsc;
      Mat x, y;
      PetscInt n, m, N, M;
      MatGetLocalSize(A->_petsc, &n, &m);
      MatGetSize(A->_petsc, &N, &M);
      MatCreate(PetscObjectComm((PetscObject)A->_petsc), &x);
      MatSetSizes(x, P == 'T' || P == 'H' ? n : m, PETSC_DECIDE,
                  P == 'T' || P == 'H' ? N : M, in->M( ));
      MatSetType(x, MATDENSE);
      {
        int size;
        MPI_Comm_size(PetscObjectComm((PetscObject)A->_petsc), &size);
        if (size == 1)
          MatSeqDenseSetPreallocation(x, &(in->operator( )(0, 0)));
        else
          MatMPIDenseSetPreallocation(x, &(in->operator( )(0, 0)));
      }
      MatAssemblyBegin(x, MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd(x, MAT_FINAL_ASSEMBLY);
      if (P == 'T' || P == 'H') {
        ffassert(in->N( ) == n);
        out->resize(m, in->M( ));
        if (P == 'T')
          MatTransposeMatMult(C, x, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &y);
        else {
          Mat w;
          MatDuplicate(x, MAT_COPY_VALUES, &w);
          MatConjugate(w);
          MatTransposeMatMult(C, w, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &y);
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
        MatCreate(PetscObjectComm((PetscObject)A->_petsc), &y);
        MatSetSizes(y, P == 'T' || P == 'H' ? m : n, PETSC_DECIDE,
                    P == 'T' || P == 'H' ? M : N, in->M( ));
        MatSetType(y, MATDENSE);
        int size;
        MPI_Comm_size(PetscObjectComm((PetscObject)A->_petsc), &size);
        if (size == 1)
          MatSeqDenseSetPreallocation(y, &(out->operator( )(0, 0)));
        else
          MatMPIDenseSetPreallocation(y, &(out->operator( )(0, 0)));
        MatAssemblyBegin(y, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(y, MAT_FINAL_ASSEMBLY);
        MatMatMult(C, x, MAT_REUSE_MATRIX, PETSC_DEFAULT, &y);
      }
      MatDestroy(&y);
      MatDestroy(&x);
      if (C != A->_petsc) MatDestroy(&C);
    }
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
          PetscStrcmp(type, MATTRANSPOSEMAT, &isType);
          if (isType) {
            b.emplace_back(std::make_pair(std::make_pair(i, j), Mat( )));
            Mat D = mat[i][j];
            Mat C;
            if (std::is_same< PetscScalar, PetscReal >::value) {
              MatTransposeGetMat(D, &b.back( ).second);
              MatTranspose(b.back( ).second, MAT_INITIAL_MATRIX, &C);
            } else {
              MatHermitianTransposeGetMat(D, &b.back( ).second);
              MatHermitianTranspose(b.back( ).second, MAT_INITIAL_MATRIX, &C);
            }
            PetscObjectReference((PetscObject)b.back().second);
            MatNestSetSubMat(A, i, j, C);
            MatDestroy(&C);
          }
        }
      }
    }
    MatConvert(A, MATMPIAIJ, MAT_INITIAL_MATRIX, B);
    for (std::pair< std::pair< PetscInt, PetscInt >, Mat > p : b) {
      Mat C;
      if (std::is_same< PetscScalar, PetscReal >::value)
        MatCreateTranspose(p.second, &C);
      else
        MatCreateHermitianTranspose(p.second, &C);
      MatDestroy(&p.second);
      MatNestSetSubMat(A, p.first.first, p.first.second, C);
      MatDestroy(&C);
    }
  }
  template< class Type >
  long MatConvert(Type* const& A, Type* const& B) {
    if (A->_petsc) {
      MatType type;
      PetscBool isType;
      MatGetType(A->_petsc, &type);
      PetscStrcmp(type, MATNEST, &isType);
      if (isType) {
        if (B->_petsc) B->dtor( );
        prepareConvert(A->_petsc, &B->_petsc);
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
          x = to< KN< PetscScalar >* >(args[1]);
          y = to< KN< PetscScalar >* >(args[2]);
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
  };
  template< class Type >
  basicAC_F0::name_and_type LinearSolver< Type >::E_LinearSolver::name_param[] = {
    {"precon", &typeid(Polymorphic*)}, {"sparams", &typeid(std::string*)}};
  template< class Type >
  AnyType LinearSolver< Type >::E_LinearSolver::operator( )(Stack stack) const {
    if (c != 2) {
      KN< PetscScalar >* in = GetAny< KN< PetscScalar >* >((*x)(stack));
      KN< PetscScalar >* out = GetAny< KN< PetscScalar >* >((*y)(stack));
      if (A) {
        Type* ptA = GetAny< Type* >((*A)(stack));
        Vec x, y;
        MatCreateVecs(ptA->_petsc, &x, &y);
        PetscInt size;
        VecGetLocalSize(y, &size);
        ffassert(in->n == size);
        if (out->n != size) {
          out->resize(size);
          *out = PetscScalar( );
        }
        VecPlaceArray(x, *in);
        VecPlaceArray(y, *out);
        if (!ptA->_ksp) {
          KSPCreate(PetscObjectComm((PetscObject)ptA->_petsc), &ptA->_ksp);
          KSPSetOperators(ptA->_ksp, ptA->_petsc, ptA->_petsc);
        }
        PetscInt N, rbegin;
        PetscScalar* tmpIn, *tmpOut;
        PetscContainer ctx;
        PetscErrorCode (*s2c)(PetscContainer, PetscScalar*, PetscScalar*), (*c2s)(PetscContainer, PetscScalar*, PetscScalar*);
        PetscObjectQuery((PetscObject)ptA->_petsc, "Hmat", (PetscObject*)&ctx);
        if(ctx) {
            PetscObjectQueryFunction((PetscObject)ctx, "s2c_C", &s2c);
            PetscObjectQueryFunction((PetscObject)ctx, "c2s_C", &c2s);
            MatGetSize(ptA->_petsc, &N, nullptr);
            MatGetOwnershipRange(ptA->_petsc, &rbegin, nullptr);
            tmpIn = new PetscScalar[N]();
            tmpOut = new PetscScalar[N]();
            std::copy_n(in->operator PetscScalar*(), in->n, tmpIn + rbegin);
            MPI_Allreduce(MPI_IN_PLACE, tmpIn, N, HPDDM::Wrapper<PetscScalar>::mpi_type(), MPI_SUM, PetscObjectComm((PetscObject)ptA->_petsc));
            (*s2c)(ctx, tmpIn, tmpOut);
            std::copy_n(tmpOut + rbegin, in->n, in->operator PetscScalar*());
        }
        if (c != 3)
          KSPSolve(ptA->_ksp, x, y);
        else {
          VecConjugate(x);
          KSPSolveTranspose(ptA->_ksp, x, y);
          VecConjugate(y);
        }
        if(ctx) {
            std::fill_n(tmpIn, N, PetscScalar());
            std::fill_n(tmpOut, N, PetscScalar());
            std::copy_n(out->operator PetscScalar*(), out->n, tmpIn + rbegin);
            MPI_Allreduce(MPI_IN_PLACE, tmpIn, N, HPDDM::Wrapper<PetscScalar>::mpi_type(), MPI_SUM, PetscObjectComm((PetscObject)ptA->_petsc));
            (*c2s)(ctx, tmpIn, tmpOut);
            std::copy_n(tmpOut + rbegin, out->n, out->operator PetscScalar*());
            delete [] tmpIn;
            delete [] tmpOut;
        }
        VecResetArray(y);
        VecResetArray(x);
        VecDestroy(&y);
        VecDestroy(&x);
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
#if PETSC_VERSION_GE(3, 14, 0)
        isType = PETSC_TRUE;
#endif
#if PETSC_VERSION_GT(3, 13, 1) && defined(PETSC_HAVE_HPDDM)
        if(isType) {
          MatGetSize(ptA->_petsc, &M, NULL);
          Mat B, C;
          MatCreateDense(PetscObjectComm((PetscObject)ptA->_ksp), in->N( ), PETSC_DECIDE, M, in->M(), &in->operator( )(0, 0), &B);
          MatCreateDense(PetscObjectComm((PetscObject)ptA->_ksp), out->N( ), PETSC_DECIDE, M, out->M(), &out->operator( )(0, 0), &C);
          int size;
          MPI_Comm_size(PetscObjectComm((PetscObject)ptA->_ksp), &size);
          if(size == 1) {
              MatType type;
              MatGetType(ptA->_petsc, &type);
              PetscStrcmp(type, MATMPIAIJ, &isType);
              if(isType) {
                  Mat Ad, Ac;
                  MatMPIAIJGetSeqAIJ(ptA->_petsc, &Ad, NULL, NULL);
                  MatConvert(Ad, MATSAME, MAT_INITIAL_MATRIX, &Ac);
                  MatHeaderReplace(ptA->_petsc, &Ac);
              }
          }
#if PETSC_VERSION_GE(3, 14, 0)
          KSPMatSolve(ptA->_ksp, B, C);
#else
          KSPHPDDMMatSolve(ptA->_ksp, B, C);
#endif
          MatDestroy(&C);
          MatDestroy(&B);
        }
        else
            ffassert(0);
#else
        HPDDM::PETScOperator op(ptA->_ksp, m
#if PETSC_VERSION_LT(3, 13, 2)
                                                 , 1
#endif
                                                      );
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
    Matrice_Creuse< HPDDM::upscaled_type<PetscScalar> >* ptK = (c == 0 ? GetAny< Matrice_Creuse< HPDDM::upscaled_type<PetscScalar> >* >((*K)(stack)) : nullptr);
    if (ptB->_A) {
      HPDDM::MatrixCSR< PetscScalar >* dA;
      if (c == 0 && ptK->A) {
        MatriceMorse< HPDDM::upscaled_type<PetscScalar> >* mA = static_cast< MatriceMorse< HPDDM::upscaled_type<PetscScalar> >* >(&(*ptK->A));
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
      if (!pList || !pList->A) {
        ptA->_A = new HpddmType(static_cast< const HPDDM::Subdomain< PetscScalar >& >(*ptB->_A));
        ptA->_A->setMatrix(dA);
        ptA->_num = new PetscInt[dA->_n];
        std::copy_n(ptB->_num, dA->_n, ptA->_num);
        ptA->_first = ptB->_first;
        ptA->_last = ptB->_last;
        initPETScStructure<false>(ptA, bs,
          nargs[1] && GetAny< bool >((*nargs[1])(stack)) ? PETSC_TRUE : PETSC_FALSE,
          static_cast< KN< double >* >(nullptr));
      } else {
        MatriceMorse< double >* mList = static_cast< MatriceMorse< double >* >(&*(pList->A));
        mList->CSR( );
        ffassert(mList->n == mList->nnz);
        int n = ptB->_A->getDof();
        ffassert(mList->m == n);
        ffassert(dA->_n == n);
        HPDDM::MatrixCSR< void > L(mList->n, n, mList->n, mList->p, mList->j, false);
        ptA->_D = new KN<PetscReal>(mList->n);
        for (int i = 0; i < mList->n; ++i) ptA->_D->operator[](i) = ptB->_A->getScaling()[mList->j[i]];
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
          dA, o, r, const_cast<MPI_Comm*>(&comm), &L);
        ptA->_num = new PetscInt[ptA->_A->getDof()];
        initPETScStructure<false>(ptA, bs,
          nargs[1] && GetAny< bool >((*nargs[1])(stack)) ? PETSC_TRUE : PETSC_FALSE, ptA->_D);
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
    Matrice_Creuse< HPDDM::upscaled_type<PetscScalar> >* ptK =
      (c == 0 || c == 3 ? GetAny< Matrice_Creuse< HPDDM::upscaled_type<PetscScalar> >* >((*K)(stack)) : nullptr);
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
      MatSetType(ptA->_petsc, MATMPIAIJ);
      if (c == 0 || c == 3) {
        PetscInt* ia = nullptr;
        PetscInt* ja = nullptr;
        PetscScalar* a = nullptr;
        bool free = true;
        if (ptK->A) {
          MatriceMorse< HPDDM::upscaled_type<PetscScalar> >* mA = static_cast< MatriceMorse< HPDDM::upscaled_type<PetscScalar> >* >(&(*ptK->A));
          ff_HPDDM_MatrixCSR< PetscScalar > dA(mA);
          ptA->_num = new PetscInt[mA->n + (ptC->_A && ptC->_A->getMatrix() ? ptC->_A->getMatrix()->_m : mA->m)];
          ptA->_cnum = ptA->_num + mA->n;
          std::copy_n(ptB->_num, mA->n, ptA->_num);
          std::copy_n(ptC->_num, (ptC->_A && ptC->_A->getMatrix() ? ptC->_A->getMatrix()->_m : mA->m), ptA->_cnum);
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
          ptA->_num = new PetscInt[ptB->_A->getMatrix( )->_n + ptC->_A->getMatrix( )->_m];
          ptA->_cnum = ptA->_num + ptB->_A->getMatrix( )->_n;
          std::copy_n(ptB->_num, ptB->_A->getMatrix( )->_n, ptA->_num);
          std::copy_n(ptC->_num, ptC->_A->getMatrix( )->_m, ptA->_cnum);
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
    PetscFunctionReturn(ret);
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
    PetscFunctionReturn(0);
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
    PetscFunctionReturn(ret);
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
    PetscFunctionReturn(0);
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
    PetscFunctionReturn(0);
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
    PetscFunctionReturn(0);
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
    PetscFunctionReturn(ret);
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
    PetscFunctionReturn(0);
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
    PetscFunctionReturn(0);
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
    PetscFunctionReturn(0);
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
          TaoSetObjectiveRoutine(tao, FormObjectiveRoutine< TaoSolver< Type > >, &user);
          TaoSetGradientRoutine(tao, FormFunction< Tao, TaoSolver< Type > >, &user);
          if (setXl || setXu) TaoSetVariableBounds(tao, xl, xu);
          TaoSetInitialVector(tao, x);
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
              TaoSetHessianRoutine(tao, ptA->_petsc, ptA->_petsc,
                                   FormJacobianTao< TaoSolver< Type >, 5 >, &user);
            } else if (user->icJ) {
              const OneOperator* codeH =
                op->Find("(", ArrayOfaType(atype< Kn* >( ), atype< Kn* >( ), false));
              user->op = new NonlinearSolver< Type >::VecF_O(in->n, stack, codeH, 2, sizeI, -1);
              TaoSetHessianRoutine(tao, ptA->_petsc, ptA->_petsc,
                                   FormJacobianTao< TaoSolver< Type >, 3 >, &user);
            } else if (user->ecJ) {
              const OneOperator* codeH =
                op->Find("(", ArrayOfaType(atype< Kn* >( ), atype< Kn* >( ), false));
              user->op = new NonlinearSolver< Type >::VecF_O(in->n, stack, codeH, 2, -1, sizeE);
              TaoSetHessianRoutine(tao, ptA->_petsc, ptA->_petsc,
                                   FormJacobianTao< TaoSolver< Type >, 4 >, &user);
            } else {
              const OneOperator* codeH = op->Find("(", ArrayOfaType(atype< Kn* >( ), false));
              user->op = new NonlinearSolver< Type >::VecF_O(in->n, stack, codeH);
              TaoSetHessianRoutine(tao, ptA->_petsc, ptA->_petsc,
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

    PetscFunctionBegin;
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
    PetscFunctionReturn(0);
  }

  template< bool U, char T, class K >
  void loopDistributedVec(Mat nest, HPDDM::Subdomain<PetscScalar>** exchange, KN_<HPDDM::upscaled_type<K>>* const& in, K* out) {
      Mat** mat;
      PetscInt M, N;
      MatNestGetSubMats(nest, &M, &N, &mat);
      Dmat** cast = reinterpret_cast<Dmat**>(exchange);
      PetscScalar* ptr = reinterpret_cast<PetscScalar*>(in->operator HPDDM::upscaled_type<PetscScalar>*());
      MatType type;
      PetscBool isType;
      if(T == 'N') {
          for(PetscInt i = 0; i < M; ++i) {
              for(PetscInt j = 0; j < N; ++j) {
                  if(mat[i][j]) {
                      MatGetType(mat[i][j], &type);
                      PetscStrcmp(type, MATMPIDENSE, &isType);
                      PetscInt n = 0, m = 0;
                      if(!isType) {
                          PetscStrcmp(type, MATTRANSPOSEMAT, &isType);
                          if(isType) {
                              Mat C;
                              if (std::is_same< PetscScalar, PetscReal >::value) MatTransposeGetMat(mat[i][j], &C);
                              else MatHermitianTransposeGetMat(mat[i][j], &C);
                              MatGetType(C, &type);
                              PetscStrcmp(type, MATMPIDENSE, &isType);
                              if(isType) MatGetSize(C, &m, &n);
                              type = MATTRANSPOSEMAT;
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
                          PetscStrcmp(type, MATTRANSPOSEMAT, &isType);
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
                              n = A->getDof();
                          }
                          else {
                              num = cast[i * N + j]->_num;
                              first = cast[i * N + j]->_first;
                              last = cast[i * N + j]->_last;
                              A = cast[i * N + j]->_A;
                              n = A->getDof();
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
        PetscScalar* p = reinterpret_cast<PetscScalar*>(u->operator HPDDM::upscaled_type<PetscScalar>*());
        if (std::is_same< typename std::remove_reference< decltype(*t.A->_A) >::type,
                          HpSchwarz< PetscScalar > >::value) {
          VecGetArray(x, &ptr);
          for(int i = 0; i < u->n; ++i)
              p[i] = u->operator[](i);
          if(isType) {
            ffassert((std::is_same<PetscReal, HPDDM::upscaled_type<PetscReal>>::value));
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
              p = reinterpret_cast<PetscScalar*>(out->operator HPDDM::upscaled_type<PetscScalar>*());
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
          if (t.A->_A) std::fill_n(out->operator HPDDM::upscaled_type<PetscScalar>*(), out->n, 0.0);
        } else {
          VecSet(x, PetscScalar( ));
          Vec isVec;
          VecCreateMPIWithArray(PETSC_COMM_SELF, 1, (*t)._A->getMatrix( )->_n,
                                (*t)._A->getMatrix( )->_n, p, &isVec);
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
        p = reinterpret_cast<PetscScalar*>(out->operator HPDDM::upscaled_type<PetscScalar>*());
        if (std::is_same< typename std::remove_reference< decltype(*t.A->_A) >::type,
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
          VecCreateMPIWithArray(PETSC_COMM_SELF, 1, (*t)._A->getMatrix( )->_n,
                                (*t)._A->getMatrix( )->_n, p, &isVec);
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
          PetscScalar* p = reinterpret_cast<PetscScalar*>(u->operator HPDDM::upscaled_type<PetscScalar>*());
          PetscStrcmp(type, MATNEST, &isType);
          if (isType) {
            ffassert((std::is_same<PetscReal, HPDDM::upscaled_type<PetscReal>>::value));
            loopDistributedVec<0, N>(t->_petsc, t->_exchange, u, ptr);
          } else {
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
          if (!t->_A) std::fill_n(out->operator HPDDM::upscaled_type<PetscScalar>*(), out->n, 0.0);
          if (isType) {
            loopDistributedVec<1, N>(t->_petsc, t->_exchange, out, ptr);
          } else {
            p = reinterpret_cast<PetscScalar*>(out->operator HPDDM::upscaled_type<PetscScalar>*());
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
      PetscInt n, m;
      MatGetSize(A.t->_petsc, &n, &m);
      ffassert(n == m);
      Ax->init(A.u->n);
      return mv(Ax, A);
    }
  };
  long destroyCSR(Dmat* p) {
    ffassert(p);
    p->dtor();
    return 0L;
  }
  template<class K, typename std::enable_if<std::is_same<K, HPDDM::upscaled_type<K>>::value>::type* = nullptr>
  KN_<K> Dmat_D(Dmat* p) {
    throwassert(p && p->_A);
    KN_<K> D(const_cast<K*>(p->_A->getScaling()), p->_A->getDof());
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
  template<class K, typename std::enable_if<std::is_same<K, HPDDM::upscaled_type<K>>::value>::type* = nullptr>
  static void init() {
    if (std::is_same< PetscInt, int >::value) {
      TheOperators->Add("<-", new PETSc::initCSRfromMatrix< HpSchwarz< PetscScalar > >);
    }
    addScalarProduct< Dmat, PetscScalar >( );

    TheOperators->Add("<-", new PETSc::initCSR< HpSchur< PetscScalar > >);

    Global.Add("ChangeNumbering", "(", new PETSc::changeNumbering< Dmat, KN >( ));
    Global.Add("ChangeNumbering", "(", new PETSc::changeNumbering< Dmat, KN >(1));
    Global.Add("ChangeNumbering", "(", new PETSc::changeNumbering< Dmat, KN >(1, 1));
    Global.Add("ChangeNumbering", "(", new PETSc::changeNumbering< Dmat, KN >(1, 1, 1));
    Global.Add("ChangeNumbering", "(", new PETSc::changeNumbering< Dmat, KNM >( ));
    Global.Add("MatMult", "(",
               new OneOperator3_< long, Dmat*, KN< PetscScalar >*, KN< PetscScalar >* >(
                 PETSc::MatMult< 'N' >));
    Global.Add("MatMatMult", "(",
               new OneOperator3_< long, Dmat*, KNM< PetscScalar >*, KNM< PetscScalar >* >(
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
    Global.Add("MatConvert", "(", new OneOperator2_< long, Dmat*, Dmat* >(PETSc::MatConvert));
    Global.Add("MatZeroRows", "(",
               new OneOperator2_< long, Dmat*, KN< double >* >(PETSc::MatZeroRows));
    Global.Add("KSPSolve", "(", new PETSc::LinearSolver< Dmat >( ));
    Global.Add("KSPSolve", "(", new PETSc::LinearSolver< Dmat >(1));
    Global.Add("KSPSolve", "(", new PETSc::LinearSolver< Dmat >(1, 1));
    if (!std::is_same< PetscScalar, PetscReal >::value)
      Global.Add("KSPSolveHermitianTranspose", "(", new PETSc::LinearSolver< Dmat >(1, 1, 1));
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
               new OneOperator3_< long, Dmat*, KN< Matrice_Creuse< HPDDM::upscaled_type<PetscScalar> > >*, KN< double >* >(
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
  template<class K, typename std::enable_if<!std::is_same<K, HPDDM::upscaled_type<K>>::value>::type* = nullptr>
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
      DMPlexCreateFromFile(comm, pB->c_str(), overlap > 0 || size == 1 ? PETSC_TRUE : PETSC_FALSE, &pA->_dm);
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
          DMPlexMarkBoundaryFaces(pdm, 111111, label);
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
#if PETSC_VERSION_GE(3, 12, 0)
              DMSetLocalSection(pA->_dm, rootSection);
              DMSetLocalSection(pA->_dm, leafSection);
#else
              DMLocalSection(pA->_dm, rootSection);
              DMLocalSection(pA->_dm, leafSection);
#endif
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
      static const int n_name_param = 0;
      static basicAC_F0::name_and_type name_param[];
      Expression nargs[n_name_param];
      DMPlexToFF_Op(const basicAC_F0& args, int d) : A(0), B(0), c(d) {
        args.SetNameParam(n_name_param, name_param, nargs);
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
  basicAC_F0::name_and_type DMPlexToFF::DMPlexToFF_Op::name_param[] = { };
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
          int nt = 0;
          for (PetscInt v = 0; v < numValues; ++v) {
              IS face;
              PetscInt size;
              const PetscInt* indices;

              DMLabelGetStratumIS(label, values[v], &face);
              ISGetLocalSize(face, &size);
              ISGetIndices(face, &indices);
              ISDestroy(&face);
              nt += size;
          }
          Triangle3 *b = new Triangle3[nt];
          Triangle3 *bb = b;
          for (PetscInt j = numValues - 1; j >= 0; --j) {
              IS face;
              PetscInt size;
              const PetscInt* indices;

              DMLabelGetStratumIS(label, values[j], &face);
              ISGetLocalSize(face, &size);
              ISGetIndices(face, &indices);
              for (PetscInt c = 0; c < size; ++c) {
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
              int iv[4], lab = 0;
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
          *pA = new Mesh3(vEnd - vStart, cEnd - cStart, nt, v, t, b);
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
          int nt = 0;
          for (PetscInt v = 0; v < numValues; ++v) {
              IS face;
              PetscInt size;
              const PetscInt* indices;

              DMLabelGetStratumIS(label, values[v], &face);
              ISGetLocalSize(face, &size);
              ISGetIndices(face, &indices);
              ISDestroy(&face);
              nt += size;
          }
          BoundaryEdge *b = new BoundaryEdge[nt];
          BoundaryEdge *bb = b;
          for (PetscInt j = numValues - 1; j >= 0; --j) {
              IS face;
              PetscInt size;
              const PetscInt* indices;

              DMLabelGetStratumIS(label, values[j], &face);
              ISGetLocalSize(face, &size);
              ISGetIndices(face, &indices);
              for (PetscInt c = 0; c < size; ++c) {
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
          for (PetscInt c = cStart; c < cEnd; ++c) {
              PetscInt *closure = NULL;
              PetscInt  closureSize, cl;
              DMPlexGetTransitiveClosure(p->_dm, c, PETSC_TRUE, &closureSize, &closure);
              int iv[3], lab = 0;
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
          Mesh* m = new Mesh(vEnd - vStart, cEnd - cStart, nt, v, t, b);
          R2 Pn, Px;
          m->BoundingBox(Pn, Px);
          m->quadtree = new Fem2D::FQuadTree(m, Pn, Px, m->nv);
          *pA = m;
          return pA;
      }
  }
} // namespace PETSc

static void Init_PETSc( ) {
  if (verbosity > 1 && mpirank == 0)
    cout << " PETSc (" << typeid(PetscScalar).name( ) << ")" << endl;
  if (exist_type< DmatC* >( )) {
    if(mpirank == 0) cout << "Cannot load both \"PETSc\"/\"SLEPc\" and \"PETSc-complex\"/\"SLEPc-complex\", please pick a single one" << endl;
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
  Global.Add("constructor", "(", new PETSc::initCSR< HpSchwarz< PetscScalar >, true >);
  Global.Add("constructor", "(", new PETSc::initCSR< HpSchwarz< PetscScalar >, true >(1));
  Global.Add("constructor", "(", new PETSc::initCSR< HpSchwarz< PetscScalar >, true >(1, 1));
  Global.Add("constructor", "(", new PETSc::initCSR< HpSchwarz< PetscScalar >, true >(1, 1, 1));
  Global.Add("MatDestroy", "(", new OneOperator1< long, Dmat* >(PETSc::destroyCSR));
  Add< Dmat* >("D", ".", new OneOperator1< KN_<double>, Dmat* >(PETSc::Dmat_D));
  Add< Dmat* >("n", ".", new OneOperator1< long, Dmat* >(PETSc::Dmat_n));
  Add< Dmat* >("range", ".", new OneOperator1s_< KN_<long>, Dmat* >(PETSc::Dmat_range));
  TheOperators->Add("<-", new PETSc::initCSRfromArray< HpSchwarz< PetscScalar > >);
  TheOperators->Add("<-", new PETSc::initCSRfromDMatrix< HpSchwarz< PetscScalar >, 0 >);
  TheOperators->Add("<-", new PETSc::initCSRfromDMatrix< HpSchwarz< PetscScalar >, 1 >);
  TheOperators->Add("<-", new PETSc::initCSRfromDMatrix< HpSchwarz< PetscScalar >, 0 >(1));
  TheOperators->Add("<-", new PETSc::initRectangularCSRfromDMatrix< HpSchwarz< PetscScalar >, 0 >);
  TheOperators->Add("<-", new PETSc::initRectangularCSRfromDMatrix< HpSchwarz< PetscScalar >, 1 >);
  TheOperators->Add("<-", new PETSc::initRectangularCSRfromDMatrix< HpSchwarz< PetscScalar >, 0 >(1));
  Global.Add("constructor", "(", new PETSc::initRectangularCSRfromDMatrix< HpSchwarz< PetscScalar >, 0 >(1, 1));
  TheOperators->Add(
    "<-", new OneOperatorCode< PETSc::initCSRfromBlockMatrix< HpSchwarz< PetscScalar > > >( ));
  TheOperators->Add(
    "=", new OneOperatorCode< PETSc::assignBlockMatrix< HpSchwarz< PetscScalar > > >( ),
         new PETSc::varfToMat< PetscScalar, Mesh, v_fes, v_fes >,
         new PETSc::varfToMat< PetscScalar, Mesh3, v_fes3, v_fes3 >,
         new PETSc::varfToMat< PetscScalar, MeshS, v_fesS, v_fesS >);
  if (std::is_same< PetscInt, int >::value) {
    TheOperators->Add("<-", new PETSc::initCSRfromMatrix< HpSchwarz< PetscScalar > >);
  }
  addProd< Dmat, PETSc::ProdPETSc, KN< HPDDM::upscaled_type<PetscScalar> >, PetscScalar, 'N' >( );
  addProd< Dmat, PETSc::ProdPETSc, KN< HPDDM::upscaled_type<PetscScalar> >, PetscScalar, 'T' >( );
  addInv< Dmat, PETSc::InvPETSc, KN< HPDDM::upscaled_type<PetscScalar> >, PetscScalar, 'N' >( );
  addInv< Dmat, PETSc::InvPETSc, KN< HPDDM::upscaled_type<PetscScalar> >, PetscScalar, 'T' >( );
  Global.Add("set", "(", new PETSc::setOptions< Dmat >( ));
  Global.Add("set", "(", new PETSc::setOptions< Dmat >(1));
  Global.Add("set", "(", new PETSc::setOptions< Dmat >(1, 1));
  Global.Add("set", "(", new PETSc::setOptions< Dmat >(1, 1, 1));
  Global.Add("exchange", "(", new exchangeIn< Dmat, PetscScalar >);
  Global.Add("exchange", "(", new exchangeInOut< Dmat, PetscScalar >);
  PETSc::init<PetscReal>();
  Global.Add("ChangeOperator", "(", new PETSc::changeOperator< Dmat >( ));
  Global.Add("ChangeOperator", "(", new PETSc::changeOperator< Dmat >(1));
  TheOperators->Add("=", new OneOperator2_< Dmat*, Dmat*, Matrice_Creuse< HPDDM::upscaled_type<PetscScalar> >* >(PETSc::changeOperatorSimple));
  TheOperators->Add("=", new OneOperator2_< Dmat*, Dmat*, Dmat* >(PETSc::changeOperatorSimple));
  Global.Add("PetscLogStagePush", "(", new OneOperator1_< long, string* >(PETSc::stagePush));
  Global.Add("PetscLogStagePop", "(", new OneOperator0< long >(PETSc::stagePop));
  Global.Add("PetscMemoryGetCurrentUsage", "(", new OneOperator0< double >(PETSc::memoryGetCurrentUsage));
  Global.Add("HasType", "(", new OneOperator2_< long, string*, string* >(PETSc::hasType));
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
