/*! \file   TridiagMatrix.hpp
    \brief  management of threads for factorization and Fw/Bw substitution
    \author Atsushi Suzuki, Laboratoire Jacques-Louis Lions
    \date   Jun. 20th 2014
    \date   Jul. 12th 2015
    \date   Nov. 30th 2016
*/

// This file is part of Dissection
// 
// Dissection is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Linking Dissection statically or dynamically with other modules is making
// a combined work based on Disssection. Thus, the terms and conditions of 
// the GNU General Public License cover the whole combination.
//
// As a special exception, the copyright holders of Dissection give you 
// permission to combine Dissection program with free software programs or 
// libraries that are released under the GNU LGPL and with independent modules 
// that communicate with Dissection solely through the Dissection-fortran 
// interface. You may copy and distribute such a system following the terms of 
// the GNU GPL for Dissection and the licenses of the other code concerned, 
// provided that you include the source code of that other code when and as
// the GNU GPL requires distribution of source code and provided that you do 
// not modify the Dissection-fortran interface.
//
// Note that people who make modified versions of Dissection are not obligated 
// to grant this special exception for their modified versions; it is their
// choice whether to do so. The GNU General Public License gives permission to 
// release a modified version without this exception; this exception also makes
// it possible to release a modified version which carries forward this
// exception. If you modify the Dissection-fortran interface, this exception 
// does not apply to your modified version of Dissection, and you must remove 
// this exception when you distribute your modified version.
//
// This exception is an additional permission under section 7 of the GNU 
// General Public License, version 3 ("GPLv3")
//
// Dissection is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Dissection.  If not, see <http://www.gnu.org/licenses/>.
//

# ifndef _DRIVER_TRIDIAGMATRIX_
# define _DRIVER_TRIDIAGMATRIX_
#include <complex>
#include <vector>
#include "Driver/C_threads_tasks.hpp"
#include "Algebra/SquareBlockMatrix.hpp"

using std::vector;

template<typename T>
class TridiagMatrix
{
public:
  TridiagMatrix(int dim, 
		int nnz, int *ptrows, int *indcols, int *indvals,
		int *remap_eqn,
		bool isSym,
		bool isMapped) :
    _nrow(dim), _nnz(nnz), _isSym(isSym), _isMapped(isMapped)
  {
    _ptrows = new int[dim + 1];
    _indcols = new int[nnz];
    _indvals = new int[nnz];
    _remap_eqn = new int[dim];
    
    for (int i = 0; i < (dim + 1); i++) {
      _ptrows[i] = ptrows[i] + 1;  // C to Fortran
    }
    for (int i = 0; i < nnz; i++) {
      _indcols[i] = indcols[i] + 1; // C to Fortran
      _indvals[i] = indvals[i] + 1; // C to Fortran
    }
    for (int i = 0; i < dim; i++) {
      _remap_eqn[i] = remap_eqn[i];
    }
    _diag =  new SquareBlockMatrix<T>; // to keep pivot information
  }
  /** Destructor.
   */
  ~TridiagMatrix()
  {
    //    FORTRAN_DECL(tridiag_free)(_tridiag_sparse);
    delete [] _ptrows;
    delete [] _indcols;
    delete [] _indvals;
    delete [] _remap_eqn;
    delete _diag;
  }

  int nrow() {
    return _nrow;
  }
  int nnz() {
    return _nnz;
  }
  void setGlobalNonzero(int nnz) {
    _nnz_global = nnz;
  }
  int nnz_global() {
    return _nnz_global;
  }
  int *ptRows() {
    return _ptrows;
  }
  int *indCols() {
    return _indcols;
  } 
  int *indVals() {
    return _indvals;
  }
  int *remap_eqn() {
    return _remap_eqn;
  }
  void setCoef(T *coef) {
    _coef = coef;
  }
  T *getCoef() {
    return _coef;
  }
  SquareBlockMatrix<T>* Diag() {
    return _diag;
  }
  bool isSym() const { return _isSym; }

  bool isMapped() const { return _isMapped; }
  
  void setPtr(void *tridiag_sparse)
  {
    _tridiag_sparse = tridiag_sparse;
  }
  void* &getPtr() {return _tridiag_sparse; }
  const void* &getPtr() const {return _tridiag_sparse; }
  
private:
  // Attributs :
  //  Dissection::Tree _btree;
  int _nrow;
  int _nnz;
  int _nnz_global;
  int *_ptrows;
  int *_indcols;
  int *_indvals;
  int *_remap_eqn;
  T* _coef;
  void *_tridiag_sparse;
  bool _isSym;
  bool _isMapped;
  SquareBlockMatrix<T> *_diag;
};

#endif
