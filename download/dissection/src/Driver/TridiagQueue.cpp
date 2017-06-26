/*! \file TridiagQueue.cpp
    \brief task mangemanet of tridiagonal factorization algorithm 
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

#include "Driver/TridiagQueue.hpp"
#include "Algebra/VectorArray.hpp"

#include <float.h>

template<typename T, typename U>
void TridiagQueue<T, U>::
generate_queue(TridiagBlockMatrix<T, U> *tridiag,
	       const int dim,
	       const int nnz,
	       const bool isMapped,
	       int *remap_eqn,
	       int *ptUnsymRows,
	       int *indUnsymCol,
	       int *indVals,
	       T *coef)
{
  _dim = dim;
  _nnz = nnz;
  _tridiag =tridiag;
  _isMapped = isMapped;
  _remap_eqn = new int[_dim];
  for (int i = 0; i < dim; i++) {
    _remap_eqn[i] = remap_eqn[i];
  }
  _ptRows= new int[_dim + 1];
  for (int i = 0; i < (dim + 1); i++) {
    _ptRows[i] = ptUnsymRows[i];
  }
  _indCols = new int[nnz];
  _indVals = new int[nnz];
  for (int i = 0; i < nnz; i++) {
    _indCols[i] = indUnsymCol[i];
    _indVals[i] = indVals[i];
  }
  _coef = coef;

  if (_tridiag_solver == false) {
    fprintf(stderr, "%s %d : tridiga_solver is not defined\n",
	    __FILE__, __LINE__);
  }
  _allocated = true;
}

template
void TridiagQueue<double>::
generate_queue(TridiagBlockMatrix<double> *tridiag,
	       const int dim,
	       const int nnz,
	       const bool isMapped,
	       int *remap_eqn,
	       int *ptUnsymRows,
	       int *indUnsymCol,
	       int *indVals,
	       double *coef);

template
void TridiagQueue<quadruple>::
generate_queue(TridiagBlockMatrix<quadruple> *tridiag,
	       const int dim,
	       const int nnz,	       
	       const bool isMapped,
	       int *remap_eqn,
	       int *ptUnsymRows,
	       int *indUnsymCol,
	       int *indVals,
	       quadruple *coef);

template
void TridiagQueue<complex<double>, double>::
generate_queue(TridiagBlockMatrix<complex<double>, double> *tridiag,
	       const int dim,
	       const int nnz,
	       const bool isMapped,
	       int *remap_eqn,
	       int *ptUnsymRows,
	       int *indUnsymCol,
	       int *indVals,
	       complex<double> *coef);


template
void TridiagQueue<complex<quadruple>, quadruple>::
generate_queue(TridiagBlockMatrix<complex<quadruple>, quadruple> *tridiag,
	       const int dim,
	       const int nnz,
	       const bool isMapped,
	       int *remap_eqn,
	       int *ptUnsymRows,
	       int *indUnsymCol,
	       int *indVals,
	       complex<quadruple> *coef);
//

template<typename T, typename U>
void TridiagQueue<T, U>::generate_queue_fwbw() {} // dummy
template
void TridiagQueue<double>::generate_queue_fwbw();

template
void TridiagQueue<quadruple>::generate_queue_fwbw();

template
void TridiagQueue<complex<double>, double>::generate_queue_fwbw();


template
void TridiagQueue<complex<quadruple>, quadruple>::generate_queue_fwbw();
//

template<typename T, typename U>
void TridiagQueue<T, U>::exec_symb_fact() 
{
  if (_tridiag_solver == false) {
    fprintf(stderr, "%s %d : tridiga_solver is not defined\n",
	    __FILE__, __LINE__);
  }
  vector<int> color_mask(_dim, 1);
  _tridiag->SymbolicFact(1, 1, &color_mask[0], _dim, // color = color_max = 1
			 _nnz, _ptRows, _indCols, _indVals);
  color_mask.clear();
}

template
void TridiagQueue<double>::exec_symb_fact();

template
void TridiagQueue<quadruple>::exec_symb_fact();

template
void TridiagQueue<complex<double>, double>::exec_symb_fact();

template
void TridiagQueue<complex<quadruple>, quadruple>::exec_symb_fact();
//

template<typename T, typename U>
void TridiagQueue<T, U>::exec_num_fact(const int called,
				       const double eps_pivot,
				       const bool kernel_detection,
				       const int aug_dim,
				       const U eps_machine)
{
  double pivot;
  vector<int> list_sing;
  double nopd;
  if (_tridiag_solver == false) {
    fprintf(stderr, "%s %d : tridiga_solver is not defined\n",
	    __FILE__, __LINE__);
  }

  _tridiag->NumericFact(_coef,
			eps_pivot,
			&pivot,
			kernel_detection,
			aug_dim,
			eps_machine,
			&nopd);
}

template
void TridiagQueue<double>::exec_num_fact(const int called,
					 const double eps_pivot,
					 const bool kernel_detection,
					 const int aug_dim,
					 const double eps_machine);

template
void TridiagQueue<quadruple>::
exec_num_fact(const int called,
	      const double eps_pivot,
	      const bool kernel_detection,
	      const int aug_dim,
	      const quadruple eps_machine);

template
void TridiagQueue<complex<double>, double>::
exec_num_fact(const int called,
	      const double eps_pivot,
	      const bool kernel_detection,
	      const int aug_dim,
	      const double eps_machine);

template
void TridiagQueue<complex<quadruple>, quadruple>::
exec_num_fact(const int called,
	      const double eps_pivot,
	      const bool kernel_detection,
	      const int aug_dim,
	      const quadruple eps_machine);
//

template<typename T, typename U>
void TridiagQueue<T, U>::exec_fwbw(T *x, const int nrhs, bool isTrans)
{
  if (_tridiag_solver == false) {
    fprintf(stderr, "%s %d : tridiga_solver is not defined\n",
	    __FILE__, __LINE__);
  }

  const int nrow = _tridiag->nrow();
  if (_isMapped) {
    if (nrhs == 1) {
      VectorArray<T> xx(nrow);
      for (int i = 0; i < nrow; i++) {
	xx[i] = x[_remap_eqn[i]];
      }
      _tridiag->SolveSingle(true, isTrans, xx.addrCoefs());
      for (int i = 0; i < nrow; i++) {
	x[_remap_eqn[i]] = xx[i];
      }
      xx.free();
    }
    else {
      ColumnMatrix<T> xx(nrow, nrhs);
      for (int n = 0; n < nrhs; n++) {
	for (int i = 0; i < nrow; i++) {
	  xx(i, n) = x[_remap_eqn[i] + n * nrow];
	}
      }
      _tridiag->SolveMulti(true, isTrans, nrhs, xx);
      for (int n = 0; n < nrhs; n++) {
	for (int i = 0; i < nrow; i++) {
	  x[_remap_eqn[i] + n * nrow] = xx(i, n);
	}
      }
      xx.free();
    }
  }
  else {
    if (nrhs == 1) {
      _tridiag->SolveSingle(true, isTrans, x);
    }
    else {
      ColumnMatrix<T> xx(nrow, nrhs, x, false);
      _tridiag->SolveMulti(true, isTrans, nrhs, xx);
    }
  }
}

template
void TridiagQueue<double>::exec_fwbw(double *x, const int nrhs, bool isTrans);

template
void TridiagQueue<quadruple>::exec_fwbw(quadruple *x, const int nrhs,
					bool isTrans);

template
void TridiagQueue<complex<double>, double>::
exec_fwbw(complex<double> *x, const int nrhs, bool isTrans);

template
void TridiagQueue<complex<quadruple>, quadruple>::
exec_fwbw(complex<quadruple> *x, const int nrhs, bool isTrans);
//
