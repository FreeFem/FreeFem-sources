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
// SUMMARY : PARDISO interface
// LICENSE : LGPLv3
// ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
// AUTHORS : Pierre Jolivet
// E-MAIL  : pierre.joliver@enseeiht.fr

/* clang-format off */
//ff-c++-LIBRARY-dep: mkl
//ff-c++-cpp-dep:
/* clang-format on */

#include <mkl_pardiso.h>
#include <mkl_spblas.h>
#include <mkl_types.h>
#if 0
#include <omp.h>
#else

extern "C" {
extern int omp_get_max_threads(void);
extern int omp_get_num_threads(void);
extern void omp_set_num_threads(int);
}
#endif
#include "ff++.hpp"

template< typename RR >
struct PARDISO_STRUC_TRAIT {
  typedef void R;
  static const int unSYM = 0;
  static const int SYM = 0;
};
template<>
struct PARDISO_STRUC_TRAIT< double > {
  typedef double R;
  static const int unSYM = 11;
  static const int SYM = 2;
};
template<>
struct PARDISO_STRUC_TRAIT< Complex > {
  typedef Complex R;
  static const int unSYM = 13;
  static const int SYM = -4;
};

void mkl_csrcsc(MKL_INT *job, MKL_INT *n, Complex *Acsr, MKL_INT *AJ0, MKL_INT *AI0, Complex *Acsc,
                MKL_INT *AJ1, MKL_INT *AI1, MKL_INT *info) {
  mkl_zcsrcsc(job, n, reinterpret_cast< MKL_Complex16 * >(Acsr), AJ0, AI0,
              reinterpret_cast< MKL_Complex16 * >(Acsc), AJ1, AI1, info);
}

void mkl_csrcsc(MKL_INT *job, MKL_INT *n, double *Acsr, MKL_INT *AJ0, MKL_INT *AI0, double *Acsc,
                MKL_INT *AJ1, MKL_INT *AI1, MKL_INT *info) {
  mkl_dcsrcsc(job, n, Acsr, AJ0, AI0, Acsc, AJ1, AI1, info);
}

template< class R >
class SolverPardiso : public VirtualSolver< int, R > {
 private:
  typedef HashMatrix< int, R > HMat;
  mutable void *pt[64];
  MKL_INT iparm[64];
  MKL_INT maxfct, mnum, phase, error, msglvl;
  /* Auxiliary variables. */
  double ddum;  /* Double dummy */
  MKL_INT idum; /* Integer dummy. */

  mutable HMat *ptA;
  MKL_INT n;
  MKL_INT *ia;
  MKL_INT *ja;
  R *a;    // Coef
  long verb;
  long cn, cs;
  MKL_INT pmtype, mtype, nrhs;

 public:
  static const int orTypeSol = 1 & 2 & 8 & 16;    // Do all
  static const MKL_INT pmtype_unset = -1000000000;
  typedef typename PARDISO_STRUC_TRAIT< R >::R MR;

  SolverPardiso(HMat &AH, const Data_Sparse_Solver &ds, Stack stack)
    : ptA(&AH), ia(0), ja(0), a(0), verb(ds.verb), cn(0), cs(0), pmtype(pmtype_unset) {
    if (verb > 2)
      cout << "   SolverPardiso " << this << " mat: " << ptA << "sym " << ds.sym << " half "
           << ptA->half << " spd " << ds.positive << endl;

    if (ds.lparams.N( ) > 1) pmtype = ds.lparams[0];    // bof bof ...
    fill(iparm, iparm + 64, 0);
    fill(pt, pt + 64, (void *)0);
    iparm[0] = 1;   /* No solver default */
    iparm[1] = 2;   /* Fill-in reordering from METIS */
    iparm[3] = 0;   /* No iterative-direct algorithm */
    iparm[4] = 0;   /* No user fill-in reducing permutation */
    iparm[5] = 0;   /* Write solution into x */
    iparm[6] = 0;   /* Not in use */
    iparm[7] = 2;   /* Max numbers of iterative refinement steps */
    iparm[8] = 0;   /* Not in use */
    iparm[9] = 13;  /* Perturb the pivot elements with 1E-13 */
    iparm[10] = 1;  /* Use nonsymmetric permutation and scaling MPS */
    iparm[11] = 0;  /* Not in use */
    iparm[12] = 0;  /* Maximum weighted matching algorithm is switched-off (default for symmetric).
                       Try iparm[12] = 1 in case of inappropriate accuracy */
    iparm[13] = 0;  /* Output: Number of perturbed pivots */
    iparm[14] = 0;  /* Not in use */
    iparm[15] = 0;  /* Not in use */
    iparm[16] = 0;  /* Not in use */
    iparm[17] = -1; /* Output: Number of nonzeros in the factor LU */
    iparm[18] = -1; /* Output: Mflops for LU factorization */
    iparm[19] = 0;  /* Output: Numbers of CG Iterations */
    maxfct = 1;     /* Maximum number of numerical factorizations. */
    mnum = 1;       /* Which factorization to use. */
    msglvl = verb > 4; /* Print statistical information in file */
    error = 0;         //(const MatriceMorse<R> &A, KN<long> &param_int, KN<double> &param_double) {
    if (ptA->half)
      mtype = -2; /* Real symmetric matrix */
    else
      mtype = 11;    // CRS
    nrhs = 0;
  }

  void UpdateState( ) {
    if (verb > 2 || verbosity > 9)
      std::cout << " UpdateState " << ptA->re_do_numerics << " " << ptA->re_do_symbolic
                << std::endl;
    if (ptA->GetReDoNumerics( )) cn++;
    if (ptA->GetReDoSymbolic( )) cs++;
    this->ChangeCodeState(ptA->n, cs, cn);
  }
  void fac_init( ) {
    ptA->setfortran(true);
    n = ptA->N;
    if (ptA->half) {
      ptA->CSR(ia, ja, a);
    } else {
      ptA->CSR(ia, ja, a);
    }
    if (verb > 99)
      for (int i = 0; i < n; ++i)
        for (int k = ia[i] - 1; k < ia[i + 1]; ++k)
          cout << i + 1 << " " << ja[k] << " " << a[k] << endl;
  }
  void fac_symbolic( ) {    // phase 11

    phase = 11;
    PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum,
            &ddum, &error);
    if (error != 0) {
      printf("\nERROR during symbolic factorization: %d\n", error);
      exit(1);
    }
    if (verb > 3) {
      printf("\n Pardiso: Reordering completed ... ");
      printf("\n      Number of nonzeros in factors = %d", iparm[17]);
      printf("\n      Number of factorization MFLOPS = %d", iparm[18]);
    }
  }
  void fac_numeric( ) {
    phase = 22;
    PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum,
            &ddum, &error);
    if (error != 0) {
      printf("\nERROR during numerical factorization: %d", error);
      exit(2);
    }
    if (verb > 3) printf("\n    Pardiso: Factorization completed ... ");
  }
  void dosolver(R *x, R *b, int N, int trans) {
    phase = 33;
    iparm[7] = 2; /* Max numbers of iterative refinement steps. */
    nrhs = N;
    iparm[11] = trans;    //
    /* Set right hand side to one. */

    PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, b, x,
            &error);
    if (error != 0) {
      printf("\nERROR during solution: %d", error);
      exit(3);
    }
    if (verbosity > 0) printf("\nSolve completed ... ");
  }

  ~SolverPardiso( ) {
    //   Warning the solver is del afer the Matrix ptA;
    if (verb > 99)
      cout << " ~SolverPardiso" << this << " " << ptA << " " << ptA->type_state << endl;
    if (ptA->type_state != HMat::type_isdeleted) {
      ptA->setfortran(false);
      ptA->HM( );
    }
    nrhs = 0;
    phase = -1; /* Release internal memory. */

    PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, &ddum, ia, ja, &idum, &nrhs, iparm, &msglvl,
            &ddum, &ddum, &error);
  }
};

static long ffompgetnumthreads( ) { return omp_get_num_threads( ); }

static long ffompgetmaxthreads( ) { return omp_get_max_threads( ); }

static long ffompsetnumthreads(long n) {
  omp_set_num_threads(n);
  return n;
}

static void Load_Init( ) {
  addsolver< SolverPardiso< double > >("PARDISO", 50, 3);
  addsolver< SolverPardiso< Complex > >("PARDISO", 50, 3);

  if (!Global.Find("ompsetnumthreads").NotNull( )) {
    Global.Add("ompsetnumthreads", "(", new OneOperator1< long, long >(ffompsetnumthreads));
  }

  if (!Global.Find("ompgetnumthreads").NotNull( )) {
    Global.Add("ompgetnumthreads", "(", new OneOperator0< long >(ffompgetnumthreads));
  }

  if (!Global.Find("ompgetmaxthreads").NotNull( )) {
    Global.Add("ompgetmaxthreads", "(", new OneOperator0< long >(ffompgetmaxthreads));
  }
}

// LOADFUNC(initPARDISO);
LOADFUNC(Load_Init)
