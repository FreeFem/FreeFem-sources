/****************************************************************************/
/* This file is part of FreeFEM.                                            */
/*                                                                          */
/* FreeFEM is free software: you can redistribute it and/or modify          */
/* it under the terms of the GNU Lesser General Public License as           */
/* published by the Free Software Foundation, either version 3 of           */
/* the License, or (at your option) any later version.                      */
/*                                                                          */
/* FreeFEM is distributed in the hope that it will be useful,               */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of           */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            */
/* GNU Lesser General Public License for more details.                      */
/*                                                                          */
/* You should have received a copy of the GNU Lesser General Public License */
/* along with FreeFEM. If not, see <http://www.gnu.org/licenses/>.          */
/****************************************************************************/
// SUMMARY : ...
// LICENSE : LGPLv3
// ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
// AUTHORS : ...
// E-MAIL  : ...

// *INDENT-OFF* //
//ff-c++-LIBRARY-dep: superlu blas
//ff-c++-cpp-dep:
// *INDENT-ON* //

#include "ff++.hpp"
#include "slu_ddefs.h"
#include "superlu_enum_consts.h"

#include "slu_zdefs.h"

template< class R >
struct SuperLUDriver {};

template<>
struct SuperLUDriver< double > {
  /* Driver routines */
  static Dtype_t R_SLU_T( ) { return SLU_D; }
  static trans_t trans( ) { return TRANS; }

  static void gssv(superlu_options_t *p1, SuperMatrix *p2, int *p3, int *p4, SuperMatrix *p5,
                   SuperMatrix *p6, SuperMatrix *p7, SuperLUStat_t *p8, int *p9) {
    dgssv(p1, p2, p3, p4, p5, p6, p7, p8, p9);
  }

  static void gssvx(superlu_options_t *p1, SuperMatrix *p2, int *p3, int *p4, int *p5, char *p6,
                    double *p7, double *p8, SuperMatrix *p9, SuperMatrix *p10, void *p11, int p12,
                    SuperMatrix *p13, SuperMatrix *p14, double *p15, double *p16, double *p17,
                    double *p18, GlobalLU_t *pGlu, mem_usage_t *p19, SuperLUStat_t *p20, int *p21) {
    dgssvx(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, pGlu,
           p19, p20, p21);
  }

  /* Supernodal LU factor related */
  static void Create_CompCol_Matrix(SuperMatrix *p1, int p2, int p3, int p4, double *p5, int *p6,
                                    int *p7, Stype_t p8, Dtype_t p9, Mtype_t p10) {
    dCreate_CompCol_Matrix(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10);
  }

  static void Create_CompRow_Matrix(SuperMatrix *p1, int p2, int p3, int p4, double *p5, int *p6,
                                    int *p7, Stype_t p8, Dtype_t p9, Mtype_t p10) {
    dCreate_CompRow_Matrix(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10);
  }

  static void Create_Dense_Matrix(SuperMatrix *p1, int p2, int p3, double *p4, int p5, Stype_t p6,
                                  Dtype_t p7, Mtype_t p8) {
    dCreate_Dense_Matrix(p1, p2, p3, p4, p5, p6, p7, p8);
  }

  static void Create_SuperNode_Matrix(SuperMatrix *p1, int p2, int p3, int p4, double *p5, int *p6,
                                      int *p7, int *p8, int *p9, int *p10, Stype_t p11, Dtype_t p12,
                                      Mtype_t p13) {
    dCreate_SuperNode_Matrix(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13);
  }

  static void CompRow_to_CompCol(int p1, int p2, int p3, double *p4, int *p5, int *p6, double **p7,
                                 int **p8, int **p9) {
    dCompRow_to_CompCol(p1, p2, p3, p4, p5, p6, p7, p8, p9);
  }
};

template<>
struct SuperLUDriver< Complex > {
  /* Driver routines */
  static Dtype_t R_SLU_T( ) { return SLU_Z; }
  static trans_t trans( ) { return CONJ; }

  static doublecomplex *dc(Complex *p) { return (doublecomplex *)(void *)p; }

  static doublecomplex **dc(Complex **p) { return (doublecomplex **)(void *)p; }

  static void gssv(superlu_options_t *p1, SuperMatrix *p2, int *p3, int *p4, SuperMatrix *p5,
                   SuperMatrix *p6, SuperMatrix *p7, SuperLUStat_t *p8, int *p9) {
    zgssv(p1, p2, p3, p4, p5, p6, p7, p8, p9);
  }

  static void gssvx(superlu_options_t *p1, SuperMatrix *p2, int *p3, int *p4, int *p5, char *p6,
                    double *p7, double *p8, SuperMatrix *p9, SuperMatrix *p10, void *p11, int p12,
                    SuperMatrix *p13, SuperMatrix *p14, double *p15, double *p16, double *p17,
                    double *p18, GlobalLU_t *pGlu, mem_usage_t *p19, SuperLUStat_t *p20, int *p21) {
    zgssvx(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, pGlu,
           p19, p20, p21);
  }

  /* Supernodal LU factor related */
  static void Create_CompCol_Matrix(SuperMatrix *p1, int p2, int p3, int p4, Complex *p5, int *p6,
                                    int *p7, Stype_t p8, Dtype_t p9, Mtype_t p10) {
    zCreate_CompCol_Matrix(p1, p2, p3, p4, dc(p5), p6, p7, p8, p9, p10);
  }

  static void Create_CompRow_Matrix(SuperMatrix *p1, int p2, int p3, int p4, Complex *p5, int *p6,
                                    int *p7, Stype_t p8, Dtype_t p9, Mtype_t p10) {
    zCreate_CompRow_Matrix(p1, p2, p3, p4, dc(p5), p6, p7, p8, p9, p10);
  }

  static void Create_Dense_Matrix(SuperMatrix *p1, int p2, int p3, Complex *p4, int p5, Stype_t p6,
                                  Dtype_t p7, Mtype_t p8) {
    zCreate_Dense_Matrix(p1, p2, p3, dc(p4), p5, p6, p7, p8);
  }

  static void Create_SuperNode_Matrix(SuperMatrix *p1, int p2, int p3, int p4, Complex *p5, int *p6,
                                      int *p7, int *p8, int *p9, int *p10, Stype_t p11, Dtype_t p12,
                                      Mtype_t p13) {
    zCreate_SuperNode_Matrix(p1, p2, p3, p4, dc(p5), p6, p7, p8, p9, p10, p11, p12, p13);
  }

  static void CompRow_to_CompCol(int p1, int p2, int p3, Complex *p4, int *p5, int *p6,
                                 Complex **p7, int **p8, int **p9) {
    zCompRow_to_CompCol(p1, p2, p3, dc(p4), p5, p6, dc(p7), p8, p9);
  }
};

int s_(char *str, const char *cmp[]) {
  int i = 0;

  while (cmp[i] != 0) {
    if (strcmp(str, cmp[i]) == 0) {
      return i + 1;
    }

    i++;
  }

  return 0;
}

void read_options_freefem(string string_option, superlu_options_t *options) {
  static const yes_no_t enumyes_no_t[2] = {NO, YES};
  static const fact_t enumfact_t[4] = {DOFACT, SamePattern, SamePattern_SameRowPerm, FACTORED};
  static const colperm_t enumcolperm_t[5] = {NATURAL, MMD_ATA, MMD_AT_PLUS_A, COLAMD, MY_PERMC};
  static const trans_t enumtrans_t[3] = {NOTRANS, TRANS, CONJ};
  static const IterRefine_t enumIterRefine_t[4] = {NOREFINE, SLU_SINGLE, SLU_DOUBLE, SLU_EXTRA};
  static const char *compyes_no_t[] = {"NO", "YES", 0};
  static const char *compfact_t[] = {"DOFACT", "SamePattern", "SamePattern_SameRowPerm", "FACTORED",
                                     0};
  static const char *compcolperm_t[] = {"NATURAL", "MMD_ATA",  "MMD_AT_PLUS_A",
                                        "COLAMD",  "MY_PERMC", 0};
  static const char *comptrans_t[] = {"NOTRANS", "TRANS", "CONJ", 0};
  static const char *compIterRefine_t[] = {"NOREFINE", "SINGLE", "DOUBLE", "EXTRA", 0};
  static const char *comp[] = {
    "Fact",          "Equil",       "ColPerm",         "DiagPivotThresh", "Trans", "IterRefine",
    "SymmetricMode", "PivotGrowth", "ConditionNumber", "PrintStat",       0};

  /* Set the default values for options argument:
   *   options.Fact = DOFACT;
   *   options.Equil = YES;
   *   options.ColPerm = COLAMD;
   *   options.DiagPivotThresh = 1.0;
   *   options.Trans = NOTRANS;
   *   options.IterRefine = NOREFINE;
   *   options.SymmetricMode = NO;
   *   options.PivotGrowth = NO;
   *   options.ConditionNumber = NO;
   *   options.PrintStat = YES;
   */
  // cout << "string_option" <<  *string_option << endl;
  KN< char > kdata(string_option.size( ) + 1);

  char *data = kdata;
  strcpy(data, string_option.c_str( ));
  cout << "data=" << data << endl;
  char *tictac;
  tictac = strtok(data, " =,\t\n");
  cout << "tictac=" << data << endl;

  while (tictac != NULL) {
    int id_option = s_(tictac, comp);
    tictac = strtok(NULL, " =,\t\n");
    int val_options;

    switch (id_option) {
      case 1:    // Fact
        val_options = s_(tictac, compfact_t);
        if (val_options == 0) {
          printf("value given for SuperLU for options %s is not correct\n", "Fact");
          exit(1);
        }

        options->Fact = enumfact_t[val_options - 1];
        break;
      case 2:    // Equil
        val_options = s_(tictac, compyes_no_t);
        if (val_options == 0) {
          printf("value given for SuperLU for options %s is not correct\n", "Equil");
          exit(1);
        }

        options->Equil = enumyes_no_t[val_options - 1];
        break;
      case 3:    // ColPerm
        val_options = s_(tictac, compcolperm_t);
        if (val_options == 0) {
          printf("value given for SuperLU for options %s is not correct\n", "ColPerm");
          exit(1);
        }

        options->ColPerm = enumcolperm_t[val_options - 1];
        break;
      case 4:    // DiagPivotThresh
        options->DiagPivotThresh = strtod(tictac, &tictac);
        break;
      case 5:    // Trans
        val_options = s_(tictac, comptrans_t);
        if (val_options == 0) {
          printf("value given for SuperLU for options %s is not correct\n", "Trans");
          exit(1);
        }

        options->Trans = enumtrans_t[val_options - 1];
        break;
      case 6:    // IterRefine
        val_options = s_(tictac, compIterRefine_t);
        if (val_options == 0) {
          printf("value given for SuperLU for options %s is not correct\n", "IterRefine");
          exit(1);
        }

        options->IterRefine = enumIterRefine_t[val_options - 1];
        break;
      case 7:    // SymmetricMode
        val_options = s_(tictac, compyes_no_t);
        if (val_options == 0) {
          printf("value given for SuperLU for options %s is not correct\n", "SymmetricMode");
          exit(1);
        }

        options->SymmetricMode = enumyes_no_t[val_options - 1];
        break;
      case 8:    // PivotGrowth
        val_options = s_(tictac, compyes_no_t);
        if (val_options == 0) {
          printf("value given for SuperLU for options %s is not correct\n", "PivotGrowth");
          exit(1);
        }

        options->PivotGrowth = enumyes_no_t[val_options - 1];
        break;
      case 9:    // ConditionNumber
        val_options = s_(tictac, compyes_no_t);
        if (val_options == 0) {
          printf("value given for SuperLU for options %s is not correct\n", "ConditionNumber");
          exit(1);
        }

        options->ConditionNumber = enumyes_no_t[val_options - 1];
        break;
      case 10:    // PrintStat
        val_options = s_(tictac, compyes_no_t);
        if (val_options == 0) {
          printf("value given for SuperLU for options %s is not correct\n", "PrintStat");
          exit(1);
        }

        options->PrintStat = enumyes_no_t[val_options - 1];
        break;
      case 0:    // Equivalent of case default
        break;
    }

    tictac = strtok(NULL, " =,\t\n");
  }

  // #endif
}

template< typename R >
class VirtualSolverSuperLU : public VirtualSolver< int, R >, public SuperLUDriver< R > {
 public:
  //  1 unsym , 2 sym, 4 pos , 8 nopos, 16  seq, 32  ompi, 64 mpi ,
  static const int orTypeSol = 1 & 8 & 16;

  typedef R K;
  typedef int Z;
  typedef HashMatrix< Z, K > HMat;
  HMat *AH;

  double tol_pivot_sym, tol_pivot;    // Add 31 oct 2005
  mutable char equed[1];
  yes_no_t equil;
  mutable SuperMatrix A, L, U;
  mutable GlobalLU_t Glu;
  R *a;
  int *asub, *xa;
  KN< int > perm_c; /* column permutation vector */
  KN< int > perm_r; /* row permutations from partial pivoting */
  string string_option;
  KN< int > etree;
  R *rhsb, *rhsx, *xact;
  KN< double > RR, CC;
  int m, n, nnz;

  mutable superlu_options_t options;
  mutable mem_usage_t mem_usage;
  long verb;
  int cs, cn;
  SuperLUStat_t stat;
  VirtualSolverSuperLU(HMat &AA, const Data_Sparse_Solver &ds, Stack stack)
    : AH(&AA), etree(0), string_option(ds.sparams), perm_r(ds.perm_r), perm_c(ds.perm_c), RR(0),
      CC(0), tol_pivot_sym(ds.tol_pivot_sym), tol_pivot(ds.tol_pivot), verb(ds.verb), cn(0), cs(0) {

    A.Store = 0;
    L.Store = 0;
    U.Store = 0;

    set_default_options(&options);
    if (AH->half) {
      cerr << " Sorry SUPERLU need a no symmetric matrix " << endl;
      cerr << " bug in choose Solver " << endl;
      ExecError("SuperLU solver");
    }
    options.SymmetricMode = AH->half ? YES : NO;
    StatInit(&stat);
  }
  void dosolver(R *x, R *b, int N, int trans) {
    if (verb > 2 || verbosity > 9)
      cout << "dosolver SuperLU double/int  " << N << " " << trans << endl;
    ffassert(trans == 0);
    options.Trans = trans ? SuperLUDriver< R >::trans( ) : NOTRANS;
    int info = 0, lwork = 0;
    void *work = 0;
    double ferr[1], berr[1];
    double rpg, rcond;

    SuperMatrix B, X;

    Dtype_t R_SLU = SuperLUDriver< R >::R_SLU_T( );
    this->Create_Dense_Matrix(&B, m, N, b, m, SLU_DN, R_SLU, SLU_GE);
    this->Create_Dense_Matrix(&X, m, N, x, m, SLU_DN, R_SLU, SLU_GE);

    SuperLUDriver< R >::gssvx(&options, &A, perm_c, perm_r, etree, equed, RR, CC, &L, &U, work,
                              lwork, &B, &X, &rpg, &rcond, ferr, berr, &Glu, &mem_usage, &stat,
                              &info);
    if (verbosity > 2) printf("Triangular solve: dgssvx() returns info %d\n", info);
    if (verbosity > 3) {
      if (info == 0 || info == n + 1) {

        if (options.IterRefine) {
          int i = 0;
          printf("Iterative Refinement:\n");
          printf("%8s%8s%16s%16s\n", "rhs", "Steps", "FERR", "BERR");
          printf("%8d%8d%16e%16e\n", i + 1, stat.RefineSteps, ferr[0], berr[0]);
        }

        fflush(stdout);
      } else if (info > 0 && lwork == -1) {
        printf("** Estimated memory: %d bytes\n", info - n);
      }
    }

    if (B.Store) Destroy_SuperMatrix_Store(&B);
    if (X.Store) Destroy_SuperMatrix_Store(&X);
  }

  void UpdateState( ) {
    if (verb > 2 || verbosity > 9)
      std::cout << " UpdateState " << AH->re_do_numerics << " " << AH->re_do_symbolic << std::endl;
    if (AH->GetReDoNumerics( )) cn++;
    if (AH->GetReDoSymbolic( )) cs++;
    this->ChangeCodeState(AH->n, cs, cn);
  }
  void fac_init( ) {
    n = AH->n;
    m = AH->m;
    nnz = AH->nnz;
    if (RR.size( ) != n) {
      RR.resize(n);
    }
    if (CC.size( ) != n) {
      CC.resize(n);
    }
    if (etree.size( ) != n) {
      etree.resize(n);
    }
    if (perm_r.size( ) != n) {
      perm_r.resize(n);
    }
    if (perm_c.size( ) != n) {
      perm_c.resize(n);
    }
    options.Fact = DOFACT;

    RR = 1.;
    CC = 1.;
  }
  void fac_symbolic( ) {
    if (verb > 2 || verbosity > 9)
      cout << "fac_symbolic SuperLU R: nnz U "
           << " nnz= " << AH->nnz << endl;
  }
  void fac_numeric( ) {
    Dtype_t R_SLU = SuperLUDriver< R >::R_SLU_T( );
    int info = 0, lwork = 0;
    void *work = 0;
    double ferr[1], berr[1];
    double rpg, rcond;

    SuperMatrix B, X;
    if (A.Store) Destroy_SuperMatrix_Store(&A);
    if (L.Store) Destroy_SuperNode_Matrix(&L);
    if (U.Store) Destroy_CompCol_Matrix(&U);
    AH->CSC(xa, asub, a);
    this->Create_CompCol_Matrix(&A, m, n, nnz, a, asub, xa, SLU_NC, R_SLU, SLU_GE);
    /* Indicate not to solve the system  ncol = 0 */
    // no X and B ..
    this->Create_Dense_Matrix(&B, m, 0, 0, m, SLU_DN, R_SLU, SLU_GE);
    this->Create_Dense_Matrix(&X, m, 0, 0, m, SLU_DN, R_SLU, SLU_GE);
    B.ncol = 0;
    options.Fact = DOFACT;
    SuperLUDriver< R >::gssvx(&options, &A, perm_c, perm_r, etree, equed, RR, CC, &L, &U, work,
                              lwork, &B, &X, &rpg, &rcond, ferr, berr, &Glu, &mem_usage, &stat,
                              &info);
    options.Fact = FACTORED;
    if (B.Store) Destroy_SuperMatrix_Store(&B);
    if (X.Store) Destroy_SuperMatrix_Store(&X);
  }
  ~VirtualSolverSuperLU( ) {
    if (A.Store) Destroy_SuperMatrix_Store(&A);
    if (L.Store) Destroy_SuperNode_Matrix(&L);
    if (U.Store) Destroy_CompCol_Matrix(&U);
    A.Store = 0;
    U.Store = 0;
    L.Store = 0;
  }
};

static void Load_Init( ) {
  addsolver< VirtualSolverSuperLU< double > >("SuperLU", 50, 1);
  addsolver< VirtualSolverSuperLU< Complex > >("SuperLU", 50, 1);
  setptrstring(def_solver, "SuperLU");
}
LOADFUNC(Load_Init)
