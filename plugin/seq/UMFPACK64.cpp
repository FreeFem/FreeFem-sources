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
// ff-c++-LIBRARY-dep: umfpack amd blas
// ff-c++-cpp-dep:
// *INDENT-ON* //

// file to add UMFPACK solver with dynamic load.

#include <iostream>
using namespace std;
#include "ff++.hpp"

template< class K = double >
class VirtualSolverUMFPACK64 : public VirtualSolver< int, K > {
 public:
  //  1 unsym , 2 sym, 4 pos , 8 nopos, 16  seq, 32  ompi, 64 mpi ,
  static const int orTypeSol = 1 & 8 & 16;
  typedef HashMatrix< int, K > HMat;
  typedef HashMatrix< SuiteSparse_long, K > HMat64;
  HMat *pA;
  HMat64 *pA64;
  VirtualSolverUMFPACK< SuiteSparse_long, K > v64;

  VirtualSolverUMFPACK64(HMat &AA, const Data_Sparse_Solver &ds, Stack stack)
    : pA(&AA), pA64(new HashMatrix< SuiteSparse_long, K >(AA)), v64(*pA64, ds, stack) {}
  void dosolver(K *x, K *b, int N, int trans) { return v64.dosolver(x, b, N, trans); }

  void fac_init( ) { v64.fac_init( ); }    // n, nzz fixe
  void fac_symbolic( ) { v64.fac_symbolic( ); }
  void fac_numeric( ) { v64.fac_numeric( ); }
  void UpdateState( ) { v64.UpdateState( ); }

  ~VirtualSolverUMFPACK64( ) { delete pA64; }
};

static void Load_Init( ) {
  addsolver< VirtualSolverUMFPACK64< double > >("UMFPACK64", 50, 1);
  addsolver< VirtualSolverUMFPACK64< Complex > >("UMFPACK64", 50, 1);
  setptrstring(def_solver, "UMFPACK64");
}
LOADFUNC(Load_Init)
