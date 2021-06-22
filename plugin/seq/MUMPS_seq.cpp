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
// AUTHORS : Frederic Hecht
// E-MAIL  : frederic.hecht@sorbonne-unviersite.fr

/* clang-format off */
//ff-c++-LIBRARY-dep: mumps_seq blas libseq fc pthread
//ff-c++-cpp-dep:
/* clang-format on */

// F. Hecht  december 2011
// ----------------------------
// file to add MUMPS sequentiel interface for sparce linear solver with dynamic load.

#include <iostream>
using namespace std;

#include "ff++.hpp"

#include "mumps_seq/mpi.h"
#include "dmumps_c.h"
#include "zmumps_c.h"

const int JOB_INIT = -1;
const int JOB_END = -2;
const int JOB_ANA = 1;
const int JOB_FAC = 2;
const int JOB_ANA_FAC = 4;
const int JOB_SOLVE = 3;
const int USE_COMM_WORLD = -987654;

template< typename RR >
struct MUMPS_STRUC_TRAIT {
  typedef void MUMPS;
  typedef void R;
};

template<>
struct MUMPS_STRUC_TRAIT< double > {
  typedef DMUMPS_STRUC_C MUMPS;
  typedef double R;
};

template<>
struct MUMPS_STRUC_TRAIT< Complex > {
  typedef ZMUMPS_STRUC_C MUMPS;
  typedef ZMUMPS_COMPLEX R;
};

void mumps_c(DMUMPS_STRUC_C *id) { dmumps_c(id); }

void mumps_c(ZMUMPS_STRUC_C *id) { zmumps_c(id); }

// template<typename R>
template< class R = double >
class SolveMUMPS_seq : public VirtualSolver< int, R > {
 public:
  static const int orTypeSol;
  typedef HashMatrix< int, R > HMat;
  typedef R K;    //
  HMat &A;

  // typedef double R;
  long verb;
  double eps;
  mutable double epsr;
  double tgv;
  int cn, cs;
  typedef typename MUMPS_STRUC_TRAIT< R >::R MR;
  mutable typename MUMPS_STRUC_TRAIT< R >::MUMPS id;
  KN< double > *rinfog;
  KN< long > *infog;

  int &ICNTL(int i) const { return id.icntl[i - 1]; }
  double &CNTL(int i) const { return id.cntl[i - 1]; }
  int &INFO(int i) const { return id.info[i - 1]; }
  double &RINFO(int i) const { return id.rinfo[i - 1]; }
  int &INFOG(int i) const { return id.infog[i - 1]; }
  double &RINFOG(int i) const { return id.rinfog[i - 1]; }

  void SetVerb( ) const {
    ICNTL(1) = 6;    //   output stream for error messages.
    ICNTL(2) = 6;    //  stream for diagnostic printing, statistics, and warning messages.
    ICNTL(3) = 6;    //  output stream global information, collected on the host.
    ICNTL(4) = min(max(verb - 2, 1L), 4L);    // the level of printing for error, warning, and diag
    if (verb == 0) ICNTL(4) = 0;
    ICNTL(11) = 0;    // noerroranalysisisperformed(nostatistics).
    if (id.job == JOB_SOLVE &&
        verb > 99) {    // computes statistics related to an error analysis of the linear system
      if (verb > 999)
        ICNTL(11) = 1;    // All Stat (veryexpensive)
      else
        ICNTL(11) = 2;    // compute main statistics
    }
  }
  void Clean( ) {
    delete[] id.irn;
    delete[] id.jcn;
    delete[] id.a;
    id.irn = 0;
    id.jcn = 0;
    id.a = 0;
  }
  void to_mumps_mat( ) {
    Clean( );

    id.nrhs = 0;    //
    int n = A.n;
    int nz = A.nnz;
    ffassert(A.n == A.m);

    int *irn = new int[nz];
    int *jcn = new int[nz];
    R *a = new R[nz];
    A.CSR( );

    for (int i = 0; i < n; ++i) {
      for (int k = A.p[i]; k < A.p[i + 1]; ++k) {
        irn[k] = i + 1;
        jcn[k] = A.j[k] + 1;
        a[k] = A.aij[k];
      }
    }

    id.n = n;
    id.nz = nz;
    id.irn = irn;
    id.jcn = jcn;
    id.a = (MR *)(void *)a;
    id.rhs = 0;
    ffassert(A.half == (id.sym != 0));    //
    ICNTL(5) = 0;                         // input matrix type
    ICNTL(7) = 7;                         // NUMBERING ...

    ICNTL(9) = 1;    // 1: A x = b, !1 : tA x = b  during slove phase
    ICNTL(18) = 0;
  }
  void Check(const char *msg = "mumps_seq") {
    if (INFO(1) != 0) {
      cout << " Erreur Mumps seq: number " << INFO(1) << endl;
      cout << " Fatal Erreur  " << msg << endl;
      Clean( );
      id.job = JOB_END;
      mumps_c(&id); /* Terminate instance */
      ErrorExec(msg, INFO(1));
    }
  }
  void CopyInfo( ) {
    if (rinfog) {
      // copy rinfog
      if (rinfog->N( ) < 40) {
        rinfog->resize(40);
      }

      for (int i = 0; i < 40; ++i) {
        (*rinfog)[i] = RINFOG(i + 1);
      }
    }

    if (infog) {
      // copy ginfo
      if (infog->N( ) < 40) {
        infog->resize(40);
      }

      for (int i = 0; i < 40; ++i) {
        (*infog)[i] = INFOG(i + 1);
      }
    }
  }
  SolveMUMPS_seq(HMat &AA, const Data_Sparse_Solver &ds, Stack stack)
    : A(AA), verb(ds.verb), eps(ds.epsilon), epsr(0), tgv(ds.tgv), cn(0), cs(0), rinfog(ds.rinfo),
      infog(ds.info) {

    int myid = 0;

    id.irn = 0;
    id.jcn = 0;
    id.a = 0;

    id.job = JOB_INIT;
    id.par = 1;
    id.sym = A.half;
    id.comm_fortran = USE_COMM_WORLD;
    SetVerb( );
    mumps_c(&id);

    Check("MUMPS_seq build/init");
    if (verbosity > 3) {
      cout << "  -- MUMPS   n=  " << id.n << ", peak Mem: " << INFOG(22) << " Mb"
           << " sym: " << id.sym << endl;
    }
  }

  ~SolveMUMPS_seq( ) {
    Clean( );
    id.job = JOB_END;
    SetVerb( );
    mumps_c(&id); /* Terminate instance */
                  /*int ierr = */
  }

  void dosolver(K *x, K *b, int N, int trans) {
    if (verbosity > 1) {
      cout << " -- MUMPS solve,  peak Mem : " << INFOG(22) << " Mb,   n = " << id.n
           << " sym =" << id.sym << " trans = " << trans << endl;
    }
    ICNTL(9) = trans == 0;    // 1: A x = b, !1 : tA x = b  during slove phase
    id.nrhs = N;
    id.lrhs = id.n;
    myscopy(id.n, b, x);
    id.rhs = (MR *)(void *)(R *)x;
    id.job = JOB_SOLVE;    // performs the analysis. and performs the factorization.
    SetVerb( );
    mumps_c(&id);
    Check("MUMPS_seq dosolver");

    if (verb > 9) {

      for (int j = 0; j < N; ++j) {
        KN_< R > B(b + j * id.n, id.n);
        cout << j << "   b linfty " << B.linfty( ) << endl;
      }
    }

    if (verb > 2) {

      for (int j = 0; j < N; ++j) {
        KN_< R > B(x + j * id.n, id.n);
        cout << "   x  " << j << "  linfty " << B.linfty( ) << endl;
      }
    }
    CopyInfo( );
  }

  void fac_init( ) { to_mumps_mat( ); }    // n, nzz fixe
  void fac_symbolic( ) {
    id.job = JOB_ANA;
    SetVerb( );
    mumps_c(&id);
    Check("MUMPS_seq Analyse");
    CopyInfo( );
  }
  void fac_numeric( ) {
    id.job = JOB_FAC;
    SetVerb( );
    mumps_c(&id);
    Check("MUMPS_seq Factorize");
    CopyInfo( );
  }
  void UpdateState( ) {
    if (A.GetReDoNumerics( )) cn++;
    if (A.GetReDoSymbolic( )) cs++;
    this->ChangeCodeState(A.n, cs, cn);
  }
};
struct InitEnd {
  InitEnd( ) {
    cout << "init MUMPS_SEQ: MPI_Init" << endl;
    int argc = 0;
    char **argv = 0;
    // BOF BOF
    MPI_Init(&argc, &argv);
  }
  ~InitEnd( ) {
    cout << "close  MUMPS_SEQ: MPI_Finalize" << endl;
    MPI_Finalize( );
  }
};
static InitEnd global;    // To init SEQ MPI ????

// 1 unsym , 2 herm, 4 sym, 8 pos , 16 nopos, 32  seq, 64  ompi, 128 mpi
template<> const int SolveMUMPS_seq<double>::orTypeSol = 1|2|4|8|16|32;
template<> const int SolveMUMPS_seq<std::complex<double>>::orTypeSol = 1|4|8|16|32;

static void Load_Init( ) {
  addsolver< SolveMUMPS_seq< double > >("MUMPS", 50, 1);
  addsolver< SolveMUMPS_seq< Complex > >("MUMPS", 50, 1);
  addsolver< SolveMUMPS_seq< double > >("MUMPSSEQ", 50, 1);
  addsolver< SolveMUMPS_seq< Complex > >("MUMPSSEQ", 50, 1);
  setptrstring(def_solver, "MUMPSSEQ");
}
LOADFUNC(Load_Init)
