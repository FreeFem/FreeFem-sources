/*! \file   TridiagBlockMatrix.cpp
    \brief  tridiagonal factorization algorithm with Cuthill-McKee
    \author François-Xavier Roux, ONERA, Laboratoire Jacques-Louis Lions
    \author Atsushi Suzuki, Laboratoire Jacques-Louis Lions
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

#include <cstdlib>
#include "Driver/TridiagBlockMatrix.hpp"
#include "Driver/C_BlasRoutines.hpp"
#include "Driver/C_KernDetect.hpp"
#include "Driver/DissectionDefault.hpp"
#include "Compiler/OptionCompiler.hpp"
#include "Compiler/arithmetic.hpp"
#include "Compiler/elapsed_time.hpp"
#include "Algebra/SparseRenumbering.hpp"
#include "Algebra/VectorArray.hpp"

template<typename T, typename U>
const T TridiagBlockMatrix<T, U>::_one = T(1.0);
template<typename T, typename U>
const T TridiagBlockMatrix<T, U>::_none = T(-1.0);
template<typename T, typename U>
const T TridiagBlockMatrix<T, U>::_zero = T(0.0);

template<typename T, typename U>
void TridiagBlockMatrix<T, U>::SymbolicFact(const int color,
					    const int color_max,
					    int *color_mask,
					    const int dim_,
					    const int nnz_,
					    const int *prow_,
					    const int *indcols_,
					    const int *indvals_)
{
  long long n1, n2;
  //  const void *fp = (void *)fp_;
  int dim1, dim2, dim3;
  _color = color;
  _color_max = color_max;
  // need count _dim and _nnz
  _nnz = 0;
  dim1 = 0;
  for (int i = 0; i < dim_; i++) {
    if (color_mask[i] == color) {
      _nnz += prow_[i + 1] - prow_[i];
      dim1++;
    }
  }
  // isolated entries are dealt by color 1
  dim3 = 0;
  if (color == 1) {
    for (int i = 0; i < dim_; i++) {
      if (color_mask[i] == 0) {
	_nnz++;
	dim3++;
      }
    }
  }
  dim2 = 0;
  for (int i = 0; i < dim_; i++) { // treatment fused block whose size is
                                   // less than DIM_AUG_KERN
    if (color_mask[i] == (-color)) {
      _nnz += prow_[i + 1] - prow_[i];
      dim2++;
    }
  }
  _dim = dim1 + dim2 + dim3;
  _new2old.resize(_dim);

  vector<int> remap_eqn, map_eqn, prow1, indcols1;
  remap_eqn.resize(dim_);
  map_eqn.resize(dim_);
  prow1.resize(dim_ + 1);
  indcols1.resize(nnz_);
  {
    int i0, i1, i2, i3;
    i0 = 0;
    i1 = dim1;
    i2 = dim1 + dim2;
    i3 = _dim; // dim1 + dim2 + dim3
    for (int i = 0; i < dim_; i++) {
      if (color_mask[i] == color) {
	remap_eqn[i0++] = i;
      }
      else if (color_mask[i] == (-color)) {
	remap_eqn[i1++] = i;
      }
      else if ((color == 1) && (color_mask[i] == 0)) {
	remap_eqn[i2++] = i;
      }
      else {
	remap_eqn[i3++] = i;
      }
    }
  }  
  for (int i = 0; i < dim_; i++) {
    map_eqn[remap_eqn[i]] = i; // map : old(original) to new (mapped) index
  }
  RenumberCSR(dim_, remap_eqn, map_eqn, 
	      prow_, indcols_, prow1, indcols1);

  vector<int> new2old, old2new;
  new2old.resize(dim_);
  old2new.resize(dim_);
    
  if (dim1 > MIN_TRIDIAG_SIZE) {
    CMK_number(dim1, &prow1[0], &indcols1[0], new2old, _verbose, _fp);
    
    for (int i = 0; i < dim3; i++) {
      _new2old[i] = remap_eqn[i + dim1 + dim2];
    }
    for (int i = 0; i < dim1; i++) {
      _new2old[i + dim3] = remap_eqn[new2old[i]];
    }
    for (int i = 0; i < dim2; i++) {
      _new2old[i + dim1 + dim3] = remap_eqn[i + dim1];
    }
    _nfront = point_front(dim1, &prow1[0], &indcols1[0], new2old,
			 _p_front);
  }
  else {
    _nfront = 1;
    _p_front.resize(2);
    _p_front[0] = 0;
    _p_front[1] = dim1 + dim2;
    for (int i = 0; i < _dim; i++) {
      _new2old[i] = remap_eqn[i];
    }
  }
  if (dim3 > 0) {
    for (int i = 1; i < _nfront; i++) {
      _p_front[i] += dim3;
    }
  }
  _p_front[_nfront] = _dim; // replace the end of block
    
  prow1.clear();
  indcols1.clear();
  //  remap_eqn.clear();
  map_eqn.clear();

  if (_verbose) {
    fprintf(_fp,
	    "%s %d : SymbolicFact : %d ", __FILE__, __LINE__, _nb);
    fprintf(_fp,
	    "color = %d dim = %d = %d + %d + %d nfront = %d : ", 
	  color, _dim, dim1, dim2, dim3, _nfront);
    for (int i = 0; i <= _nfront; i++) {
      fprintf(_fp, "%d ", _p_front[i]);
    }
    fprintf(_fp, "\n");
  } // if (_verbose)
  _ptRows.resize(_dim + 1);
  _indCols.resize(_nnz);
  _indVals.resize(_nnz);

  for (int i = 0; i < _dim; i++) {
      new2old[i] = _new2old[i];
  }
  for (int i = _dim; i < dim_; i++) {
    new2old[i] = remap_eqn[i];
  }
  for (int i = 0; i < dim_; i++) {
    old2new[new2old[i]] = i;
  }
  {
    bool shrink_flag = (_dim < dim_);
    RenumberCSR(shrink_flag,
		_dim, nnz_, new2old, old2new, prow_, indcols_, indvals_, 
	      _ptRows, _indCols, _indVals);
  }
  new2old.clear();
  old2new.clear();
  
  _p_diag.resize(_dim);
  _p_upper.resize(_dim);
  
  TridiagStruct(_ptRows, _indCols, _nfront, _p_front, _p_diag, _p_upper);

  if (_isSymmetric) {
    _nop = 0.0;
    n2 = (double)(_p_front[1] - _p_front[0]);
    _nop += (n2 * n2 * n2) / 3.0;
    for (int n = 1; n < _nfront; n++) {
      n1 = n2;
      n2 = (double)(_p_front[n + 1] - _p_front[n]);
      _nop += (n1 * n1 * n2) / 3.0;
      _nop += (n1 * n2 * n2) / 3.0;
      _nop += (n2 * n2 * n2) / 3.0;
    }
  }
  else {
    _nop = 0.0;
    n2 = (double)(_p_front[1] - _p_front[0]);
    _nop += (2.0 * n2 * n2 * n2) / 3.0;
    for (int n = 1; n < _nfront; n++) {
      n1 = n2;
      n2 = (double)(_p_front[n + 1] - _p_front[n]);
      _nop += (2.0 * n1 * n1 * n2) / 3.0;
      _nop += (2.0 * n1 * n2 * n2) / 3.0;
      _nop += (2.0 * n2 * n2 * n2) / 3.0;
    }
  }
  _diag_blocks = new ColumnMatrix<T>[_nfront];
#ifndef SPARSE_OFFDIAG
  _upper_blocks = new ColumnMatrix<T>[_nfront];  // 27 Nov.2015
  _lower_blocks = new ColumnMatrix<T>[_nfront];  // 27 Nov.2015
#endif
#ifdef STORE_WHOLE_FACTORIZED
  _factorized_whole.init(_dim, _dim);
  _factorized_whole.ZeroClear();
#endif
  _diag_block_alloc_status = true;
  _maxdim = 0;
  for (int n = 0; n < _nfront; n++) {
    const int itmp = _p_front[n + 1] - _p_front[n];
    _maxdim = _maxdim > itmp ? _maxdim : itmp;
  }
  _n0 = 0; // initialization
  _nscol = 0;
  remap_eqn.clear();
}

template
void TridiagBlockMatrix<double>::SymbolicFact(const int color,
					      const int color_max,
					      int *color_mask,
					      const int dim_,
					      const int nnz_,
					      const int *prow_,
					      const int *indcols_,
					      const int *indvals_);
template
void TridiagBlockMatrix<quadruple>::SymbolicFact(const int color,
						 const int color_max,
						 int *color_mask,
						 const int dim_,
						 const int nnz_,
						 const int *prow_,
						 const int *indcols_,
						 const int *indvals_);
template
void TridiagBlockMatrix<complex<double>, double>::
SymbolicFact(const int color,
	     const int color_max,
	     int *color_mask,
	     const int dim_,
	     const int nnz_,
	     const int *prow_,
	     const int *indcols_,
	     const int *indvals_);

template
void TridiagBlockMatrix<complex<quadruple>, quadruple>::
SymbolicFact(const int color,
	     const int color_max,
	     int *color_mask,
	     const int dim_,
	     const int nnz_,
	     const int *prow_,
	     const int *indcols_,
	     const int *indvals_);
//

template<typename T, typename U>
void TridiagBlockMatrix<T, U>::NumericFact(const T* coef,
					   const double eps_pivot,
					   double *pivot,
					   const bool kernel_detection,
					   const int dim_aug_kern,
					   const U eps_machine,
					   double *nopd)
{
  int nsing, n0, nn0;
  ColumnMatrix<T> a12, a21, a22, s22, b12, b21, b22;
  double eps_piv, pivot_val;
  bool flag_repeat_piv;
  vector<int> perm;
  ColumnMatrix<T>* diag_block_save;
  vector<int> num_null_aug;
  vector<int> list_sing;
  const double machine_eps = todouble<U>(eps_machine);
  
  _coef = coef;
  diag_block_save = new ColumnMatrix<T>[_nfront];
  TridiagNumericFact(pivot, eps_pivot,
		     dim_aug_kern,
		     diag_block_save,
		     num_null_aug, nopd);
  nsing = 0;
  for (int n = 0; n < _nfront; n++) {
    nsing += _num_null[n];
  }

  if (nsing == 0) {
    _nscol = 0;
    _detected = true;   // without using kernel detection
    _n0 = 0;
    _list_elim.clear();
    _list_schur.clear();
    _a12.free();
    _a21.free();
    _s22.free();
  }
  else { // if (nsing > 0)
    _detected = false; // initialization : detection is not activated
    if (_verbose) {
      fprintf(_fp, "%s %d nsing  = %d nfront=%d : ",
	      __FILE__, __LINE__, nsing, _nfront);
      for (int n = 0; n < _nfront; n++) {
	fprintf(_fp, "%d/%d ", _num_null[n], (_p_front[n + 1] - _p_front[n]));
	if (num_null_aug[n] > 0) {
	  fprintf(_fp, "(+%d)", num_null_aug[n]);
	}
      }
      fprintf(_fp, "\n");
    } // if (_verbose)
    s22.init(nsing, nsing);
    perm.resize(nsing);
    list_sing.resize(nsing);
    nsing = 0;
    for (int n = 0; n < _nfront; n++) {
      const int offset = _p_front[n];
      const int dim1 = _p_front[n + 1] - offset;
      const int num_nonsing = dim1 - _num_null[n];
      for (int i = num_nonsing; i < dim1; i++, nsing++) { // singular entries
	list_sing[nsing] = _permute[i + offset] + offset;  
      }
    }
#if 0
    if (_verbose) {
      fprintf(_fp, "%s %d : list_sing : %d\n", __FILE__, __LINE__, nsing);
      for (int i = 0; i < nsing; i++) {
	fprintf(_fp, "%d ", list_sing[i]);
      }
      fprintf(_fp, "\n");
    }
#endif
    a12.init(_dim, nsing);
    a21.init(_dim, nsing);
    a22.init(nsing, nsing);
    ComputeSchurComplement(nsing, list_sing, a12, a21, a22);

    eps_piv = eps_pivot; // initial is given by user and then repeated
    flag_repeat_piv = true;
    int count_repeat = 0;
    while (flag_repeat_piv) {
      s22.copy(a22);
      // find maximum diagonal from s22[]
      pivot_val = 0.0;
      for (int i = 0; i < nsing; i++) {
  		                                     // pivot value in double
	const double tmp = blas_abs<T, double>(s22(i, i));
	pivot_val = tmp > pivot_val ? tmp : pivot_val;
      }
      n0 = 0; // by assuming the matrix is invertible
      if (_isSymmetric) {
	double fop;
	full_ldlt_permute<T, U>(&nn0, n0, nsing, s22.addrCoefs(), nsing,
				&pivot_val, &perm[0], eps_piv, &fop);
      }
      else {
	double fop;
	full_ldu_permute<T, U>(&nn0, n0, nsing, s22.addrCoefs(), nsing,
			       &pivot_val, &perm[0], eps_piv, &fop);
      }
      if (_verbose) {
	fprintf(_fp, "%s %d factorize dim = %d eps_piv = %g -> sing = %d\n",
		__FILE__, __LINE__, nsing, eps_piv, nn0);
      }
      //
      if (((nsing - dim_aug_kern) >= nn0) || (eps_piv < TOL_PIVOT)) {
	flag_repeat_piv = false;
	break;
      }
      else {
	eps_piv /= 10.0;
	count_repeat++;
      }
    } // while (flag_repeat_piv)
    n0 = nn0;
    if (_verbose) {
      fprintf(_fp, "%s %d : n0 = %d count_repeat = %d\n", __FILE__, __LINE__,
	    n0, count_repeat);
    }
#if 0
    if (_verbose) {
      fprintf(_fp, "%s %d : s22() dim = %d\n", __FILE__, __LINE__, nsing);
      for (int i = 0; i < nsing; i++) {
	fprintf(_fp, "%d ", i);
	for (int j = 0; j < nsing; j++) {
	  printscalar<T>(_fp, s22(i, j));
	}
	fprintf(_fp, "\n");
      }
    }
#endif
    if ((n0 == 0) && (count_repeat == 0)) {
      _n0 = 0;
      _detected = true;  // without kernel check
    }
    else if (!kernel_detection) {
      _n0 = (-1);  // for safety
      _detected = false;
    }
    else { // if ((n0 == 0) && (count_repeat == 0)) && kernel_detection
      int nsing1, nsing2;
      if (n0 == 0) {
	nsing1 = nsing;
	b22.init(nsing1, nsing1);
	b22.copy(a22); //
      }
      else { // if (n0 > 0) { // compute Schur complement again
	nsing1 = n0 + dim_aug_kern;
	nsing2 = nsing - nsing1;
	if (nsing2 > 0) {
	  b21.init(nsing2, nsing1);
	  b12.init(nsing2, nsing1);
	  b22.init(nsing1, nsing1);
	  if (_isSymmetric) {
	    for (int j = 0; j < nsing1; j++) {
	      const int jj = perm[nsing2 + j];
	      for (int i = 0; i <= j; i++) {
		const int ii = perm[nsing2 + i];
		b22(i, j) = a22(ii, jj);
	      }
	    }
	    // symmetrize
	    for (int j = 0 ; j < nsing1; j++) {
	      for (int i = 0; i < j; i++) {
		b22(j, i) = b22(i, j);
	      }
	    }
	    for (int i = 0; i < nsing2; i++) {
	      const int ii = perm[i];
	      for (int j = 0; j < nsing1; j++) {
		const int jj = perm[nsing2 + j];
		b12(i, j) = a22(ii, jj);
	      }
	    }
	    // alpha = 1
	    blas_trsm<T>(CblasLeft, CblasLower, CblasNoTrans, CblasUnit,
			 nsing2, nsing1, _one, s22.addrCoefs(), nsing,
			 b12.addrCoefs(), nsing2);
	    for (int i = 0; i < nsing2; i++) {
	      for (int j = 0; j < nsing1; j++) {
		b21(i, j) = b12(i, j) * s22(i, i);
	      }
	    }
	  } // _isSymmetric
	  else {
	    for (int j = 0; j < nsing1; j++) {
	      const int jj = perm[nsing2 + j];
	      for (int i = 0; i < nsing1; i++) {
		const int ii = perm[nsing2 + i];
		b22(i, j) = a22(ii, jj);
	      }
	    }
	    for (int i = 0; i < nsing2; i++) {
	      const int ii = perm[i];
	      for (int j = 0; j < nsing1; j++) {
		const int jj = perm[nsing2 + j];
		b12(i, j) = a22(ii, jj);
		b21(i, j) = a22(jj, ii);
	      }
	    }
	  // alpha = T(1)
	    blas_trsm<T>(CblasLeft, CblasLower, CblasNoTrans, CblasUnit,
			 nsing2, nsing1, _one, s22.addrCoefs(), nsing,
			 b12.addrCoefs(), nsing2);
	    blas_trsm<T>(CblasLeft, CblasLower, CblasTrans, CblasUnit,
			 nsing2, nsing1, _one, s22.addrCoefs(), nsing,
			 b21.addrCoefs(), nsing2);
	    for (int i = 0; i < nsing2; i++) {
	      for (int j = 0; j < nsing1; j++) {
		b12(i, j) *= s22(i, i);
	      }
	    }
	  } // _isSymmetric
	  // alpha = -1
	  // beta = 1
	  blas_gemm<T>(CblasTrans, CblasNoTrans, nsing1, nsing1, nsing2, _none,
		     b21.addrCoefs(), nsing2, b12.addrCoefs(), nsing2, _one,
		     b22.addrCoefs(), nsing1);
	  b21.free();
	  b12.free();
	} // if (nsing2 > 0)
	else {
	  b22.init(nsing1, nsing1);
	  b22.copy(a22); //
	}
      } // if (n0 > 0)
      bool flag, flag_2x2;
      flag = ComputeDimKernel<T, U>(&nn0, &flag_2x2, b22.addrCoefs(),
				    nsing1, _isSymmetric, 
				    dim_aug_kern, eps_machine, eps_piv,
				    _verbose, _fp);
      if (!flag) {
	_n0 = (-1);     // for safety
	_detected = false;
      }
      else {
	_n0 = nn0;
	_detected = true;
      }
    } // if ((n0 == 0) && (count_repeat == 0)) && kernel_detection
    for (int n = 0; n < _nfront; n++) {
      if (diag_block_save[n].size() > 0) {
	_diag_blocks[n].copy(diag_block_save[n]);
	_num_null[n] -= num_null_aug[n];
	num_null_aug[n] = 0; // for safety
      }
      diag_block_save[n].free(); // release memory
    }
    // recounting "nsing" after determinated singularity with restoring aug_dim
    nsing = 0;
    for (int n = 0; n < _nfront; n++) {
      const int offset = _p_front[n];
      const int dim1 = _p_front[n + 1] - offset;
      const int num_nonsing = dim1 - _num_null[n];
      for (int i = num_nonsing; i < dim1; i++, nsing++) { 
	list_sing[nsing] = _permute[i + offset] + offset;
      }
    }
#if 0
    if (_verbose) {
      fprintf(_fp, "%s %d : list_sing in tridiag: nsing = %d : ",
	      __FILE__, __LINE__, nsing);
      for (int i = 0; i < nsing; i++) {
	fprintf(_fp, "%d ", list_sing[i]);
      }
      fprintf(_fp, "\n");
      fprintf(_fp, "%s %d : singular entries : nsing = %d : ",
	      __FILE__, __LINE__, nsing);
      for (int i = 0; i < nsing; i++) {
	fprintf(_fp, "%d ", _new2old[list_sing[i]]);
      }
      fprintf(_fp, "\n");
    } // if (_verbose)
#endif
    a12.free();
    a21.free();
    a21.free();
    s22.free();
    if ((nsing == _n0) || (!_detected)) { // nsing >= _n0
      _nscol = 0;
      _list_schur.clear();
      _list_elim.resize(nsing);
      for (int i = 0; i < nsing; i++){
	_list_elim[i] = list_sing[i];
      }
      _a12.free();
      _a21.free();
      _s22.free();
      if (!_detected) {
	_n0 = nsing;
      }
    }
    else { // if ((nsing < _n0)  && _detected)
      a12.init(_dim, nsing);
      a21.init(_dim, nsing);
      a22.init(nsing, nsing);
      ComputeSchurComplement(nsing, list_sing, a12, a21, a22);
#if 0
      // debug : begin
      if (_verbose) {
	fprintf(_fp, "%s %d : a22() dim = %d\n", __FILE__, __LINE__, nsing);
	for (int i = 0; i < nsing; i++) {
	  fprintf(_fp, "%d : ", i);
	  for (int j = 0; j < nsing; j++) {
	    printscalar<T>(_fp, a22(i, j));
	  }
	  fprintf(_fp, "\n");
	}
	
	fprintf(_fp, "%s %d : a12() dim = %d nsing = %d \n",
		__FILE__, __LINE__, _dim, nsing);
	for (int j = 0; j < nsing; j++) {
	  fprintf(_fp, "%d : ", j);
	  for (int  i = 0; i < _dim; i++) {
	    printscalar<T>(_fp, a12(i, j));
	  }
	  fprintf(_fp, "\n");
	}
      } // if (_verbose)
#endif
    // debug : end
      pivot_val = 0.0;
      for (int i = 0; i < nsing; i++) {
	                                           // pivot value in double
	const double tmp = blas_abs<T, double>(a22(i, i));  
	pivot_val = tmp > pivot_val ? tmp : pivot_val;
      }
      if (_verbose) {
	fprintf(_fp, "%s %d : _n0 = %d nsing = %d\n",
		__FILE__, __LINE__, _n0, nsing);
      }
      // factorize Schur complement by knowing the dimension of the kernle : _n0
      if (_isSymmetric) {
	double fop;
	bool flag;
	flag = full_ldlt_permute<T, U>(&nn0, _n0, nsing, a22.addrCoefs(), nsing,
				       &pivot_val, &perm[0], machine_eps, &fop);
	if (_verbose && !flag) {
	  fprintf(_fp, "%s %d : full_ldlt_permute fails : %d != %d\n",
		  __FILE__, __LINE__, _n0, nn0);
	}
      }
      else {
	double fop;
	bool flag;
	flag = full_ldu_permute<T, U>(&nn0, _n0, nsing, a22.addrCoefs(), nsing,
				      &pivot_val, &perm[0], machine_eps, &fop);
	if (_verbose && !flag) {
	fprintf(_fp, "%s %d : full_ldu_permute fails : %d != %d\n",
		__FILE__, __LINE__, _n0, nn0);
	}
      }
      _nscol = nsing - _n0;
      if (_verbose) {
	fprintf(_fp, "%s %d : detected = %s : _nscol = %d _n0 = %d : perm[] = ",
		__FILE__, __LINE__, (_detected ? "true" : "false"),
		_nscol, _n0);
	for (int i = 0; i < nsing; i++) {
	  fprintf(_fp, "%d ", perm[i]);
	}
	fprintf(_fp, "\n");
      } // if (_verbose)
      _list_schur.resize(_nscol);
      if (_n0 > 0) {
	_list_elim.resize(_n0);
      }
      for (int i = 0; i < _nscol; i++){
	_list_schur[i] = _permute_ginv[list_sing[perm[i]]]; // regular part
      }
      for (int i = 0; i < _n0; i++){
	_list_elim[i] = list_sing[perm[i + _nscol]];
      }

      if (_verbose) {
	fprintf(_fp, "%s %d : _nscol = %d : _list_schur[] = ",
		__FILE__, __LINE__, _nscol);
	for (int i = 0; i < _nscol; i++) {
	  fprintf(_fp, "%d ", _list_schur[i]);
	}
	fprintf(_fp, "\n");
	if (_n0 > 0) {
	  fprintf(_fp, "%s %d : n0 = %d : _list_elim[] = ",
		  __FILE__, __LINE__, _n0);
	  for (int i = 0; i < _n0; i++) {
	    fprintf(_fp, "%d ", _list_elim[i]);
	  }
	  fprintf(_fp, "\n");
	  
	  fprintf(_fp, "%s %d : n0 = %d : _new2old[_list_elim[]] = ",
		  __FILE__, __LINE__, _n0);
	  for (int i = 0; i < _n0; i++) {
	    fprintf(_fp, "%d ", _new2old[_list_elim[i]]);
	  }
	  fprintf(_fp, "\n");
	} // if (_n0 > 0)
      }   // if (_verbose)
      
      _a12.init(_dim, _nscol);
      for (int j = 0; j < _nscol; j++) {
	for (int i = 0; i < _dim; i++) {
	  _a12(_permute_ginv[i], j) = a12(i, perm[j]);
	}
      }
      if (!_isSymmetric) {
	_a21.init(_dim, _nscol);
	for (int j = 0; j < _nscol; j++) {
	  for (int i = 0; i < _dim; i++) {
	    _a21(_permute_ginv[i], j) = a21(i, perm[j]);
	  }
	}
      }
      else {
	_a21.free();
      }
      _s22.init(_nscol, _nscol);
      // a22 is already permuted by perm[] during LDLt/LDU,
      // _s22 is smaller than a22
      for (int i = 0; i < _nscol; i++) {
	for (int j = 0; j < _nscol; j++) {
	  _s22(i, j) = a22(i, j);
	}
      }
      a12.free();
      a21.free();
      a22.free();
    }
    perm.clear();
  }   // if (nsing > 0) 
  delete [] diag_block_save;
}

template
void TridiagBlockMatrix<double>::
NumericFact(const double* coef,
	    const double eps_pivot,
	    double *pivot,
	    const bool kernel_detection,
	    const int dim_aug_kern,
	    const double eps_machine,
	    double *nopd);

template
void TridiagBlockMatrix<quadruple>::
NumericFact(const quadruple* coef,
	    const double eps_pivot,
	    double *pivot,
	    const bool kernel_detection,
	    const int dim_aug_kern,
	    const quadruple eps_machine,
	    double *nopd);

template
void TridiagBlockMatrix<complex<double>, double>::
NumericFact(const complex<double>* coef,
	    const double eps_pivot,
	    double *pivot,
	    const bool kernel_detection,
	    const int dim_aug_kern,
	    const double eps_machine,
	    double *nopd);

template
void TridiagBlockMatrix<complex<quadruple>, quadruple>::
NumericFact(const complex<quadruple>* coef,
	    const double eps_pivot,
	    double *pivot,
	    const bool kernel_detection,
	    const int dim_aug_kern,
	    const quadruple eps_machine,
	    double *nopd);
//

template<typename T, typename U>
void TridiagBlockMatrix<T, U>::ComputeSchurComplement(const int nsing,
						      vector<int> &list_sing,
						      ColumnMatrix<T>& a12,
						      ColumnMatrix<T>& a21,
						      ColumnMatrix<T>& a22)
{
  // _nscol should be set 0 before calling
  a12.ZeroClear();
  a21.ZeroClear();
  a22.ZeroClear();
  for (int j = 0; j < nsing; j++) {
    extract_column(list_sing[j], a12.addrCoefs() + (j * _dim));	
  }
  for (int j = 0; j < nsing; j++) {
    for (int i = 0; i < nsing; i++) {
      a22(i, j) = a12(list_sing[i], j);
    }
  }
  for (int j = 0; j < nsing; j++) {
    for (int i = 0; i < nsing; i++) {
      a12(list_sing[i], j) = _zero;
    }
  }
  if (_isSymmetric) {
    a21.copy(a12);
  }
  else {
    for (int j = 0; j < nsing; j++) {
      extract_row(list_sing[j], a21.addrCoefs() + (j * _dim));
    }
    for (int j = 0; j < nsing; j++) {
      for (int i = 0; i < nsing; i++) {
	a21(list_sing[i], j) = _zero;
      }
    }
  }
  SolveMulti(false, false, nsing, a12, 0);
  // alpha = -1
  // beta = 1
  blas_gemm<T>(CblasTrans, CblasNoTrans, nsing, nsing, _dim, _none,
	       a21.addrCoefs(), _dim,  a12.addrCoefs(), _dim, _one,
	       a22.addrCoefs(), nsing);
  if (_isSymmetric) {
  // symmetrize
    for (int i = 0; i < nsing; i++) {
      for (int j = i + 1; j < nsing; j++) {
	a22(j, i) = a22(i, j);  // acess lower
      }
    }
  }
  else {
    SolveMulti(false, true, nsing, a21, 0);
  }
}

template
void TridiagBlockMatrix<double>::
ComputeSchurComplement(const int nsing, vector<int> &list_sing,
		       ColumnMatrix<double>& a12,
		       ColumnMatrix<double>& a21,
		       ColumnMatrix<double>& a22);

template
void TridiagBlockMatrix<quadruple>::
ComputeSchurComplement(const int nsing, vector<int> &list_sing,
		       ColumnMatrix<quadruple>& a12,
		       ColumnMatrix<quadruple>& a21,
		       ColumnMatrix<quadruple>& a22);

template
void TridiagBlockMatrix<complex<double>, double>::
ComputeSchurComplement(const int nsing, vector<int> &list_sing,
		       ColumnMatrix<complex<double> >& a12,
		       ColumnMatrix<complex<double> >& a21,
		       ColumnMatrix<complex<double> >& a22);


template
void TridiagBlockMatrix<complex<quadruple>, quadruple>::
ComputeSchurComplement(const int nsing, vector<int> &list_sing,
		       ColumnMatrix<complex<quadruple> >& a12,
		       ColumnMatrix<complex<quadruple> >& a21,
		       ColumnMatrix<complex<quadruple> >& a22);
//

template<typename T, typename U>
void TridiagBlockMatrix<T, U>::SingularNode(vector<int> &list_sing)
{
  if (_n0 >= 0) {
    list_sing.resize(_n0);
    for (int i = 0; i < _n0; i++) {
      list_sing[i] = _new2old[_list_elim[i]];
    }
  }
  else {
    fprintf(_fp, "%s %d : _n0 = %d\n", __FILE__, __LINE__, _n0);
  }
}

template
void TridiagBlockMatrix<double>::
SingularNode(vector<int> &list_sing);

template
void TridiagBlockMatrix<quadruple>::
SingularNode(vector<int> &list_sing);

template
void TridiagBlockMatrix<complex<double>, double>::
SingularNode(vector<int> &list_sing);

template
void TridiagBlockMatrix<complex<quadruple>, quadruple>::
SingularNode(vector<int> &list_sing);
//

//#define DEBUG_FACTORIZATION
template<typename T, typename U>
void TridiagBlockMatrix<T, U>::
TridiagNumericFact(double *pivot,
		   const double eps_pivot,
		   const int dim_aug_kern,
		   ColumnMatrix<T> *diag_block_save,
		   vector<int> &num_null_aug,
		   double *nopd)
{
  int dim1, dim2, nn0;
  const int n0 = 0;
  int n0_total = 0;
  ColumnMatrix<T> diag_works, diag_origs; // supposing "dim2 <= size_b2"
#ifdef SPARSE_OFFDIAG
  ColumnMatrix<T> upper, lower;
#endif
  int upper_size_max;
  int diag_size_max;
  int *permute_diag, *permute_offdiag, *permute_work;
  int *permute_diag_inv, *permute_offdiag_inv;
  double fop;
  //  const T zero(0.0);
  
  _num_null.resize(_nfront);
  //  _null_lists_local.resize(_nfront);
  _permute.resize(_dim);
  *pivot = 1.0;
  upper_size_max = 0;
  for (int n = 1; n < _nfront; n++) {
    const int itmp = ((_p_front[n] - _p_front[n - 1]) *
		      (_p_front[n + 1] - _p_front[n]));
    upper_size_max = itmp > upper_size_max ? itmp : upper_size_max;
  }
#ifdef SPARSE_OFFDIAG
  VectorArray<T> upper_work(upper_size_max);
  VectorArray<T> lower_work(upper_size_max);
#else
  VectorArray<T> permute_vals(upper_size_max);
#endif
  diag_size_max = 0;
  for (int n = 0; n < _nfront; n++) {
    const int itmp = _p_front[n + 1] - _p_front[n];
    _diag_blocks[n].init(itmp, itmp);
    diag_size_max  = itmp > diag_size_max ? itmp : diag_size_max;
  }
  permute_offdiag = new int[diag_size_max];
  permute_offdiag_inv = new int[diag_size_max];
  permute_diag_inv = new int[diag_size_max];
  permute_work = new int[diag_size_max];
#ifdef DEBUG_FACTORIZATION
  ColumnMatrix<T> diag_orig(diag_size_max, diag_size_max);
#endif
  { // n = 0
    dim2 = _p_front[1] - _p_front[0];
    ColumnMatrix<T> &diag = _diag_blocks[0];
    diag.ZeroClear();
    permute_diag = &_permute[0] + _p_front[0];
    
    for (int i = _p_front[0]; i < _p_front[1]; i++) {
      for (int k = _p_diag[i]; k < _p_upper[i]; k++) {
	const int ii = i - _p_front[0];
	const int jj = _indCols[k] - _p_front[0];
	diag(ii, jj) = _coef[_indVals[k]];
      }
    }
#ifdef DEBUG_FACTORIZATION
    diag_origs.init(dim2, dim2, diag_orig.addrCoefs(), false);
    for (int j = 0; j < dim2; j++) {
      for (int i = 0; i < dim2; i++) {
	diag_origs(i, j) = diag(i, j);
      }
    }
#endif
    if (_isSymmetric) {
      bool flag;
      flag = full_ldlt_permute<T, U>(&nn0, n0, dim2,
				     _diag_blocks[0].addrCoefs(), dim2, pivot,
				     permute_diag, eps_pivot, &fop);
    }
    else {
      bool flag;
      flag = full_ldu_permute<T, U>(&nn0, n0, dim2,
				    _diag_blocks[0].addrCoefs(), dim2, pivot,
				    permute_diag, eps_pivot, &fop);
    }
#ifdef STORE_WHOLE_FACTORIZED
    for (int j = 0; j < dim2; j++) {
      for (int i = 0; i < dim2; i++) {
	_factorized_whole(i, j) = diag(i, j);
      }
    }
#endif
#ifdef DEBUG_FACTORIZATION
    if (_verbose) {
      fprintf(_fp, "%s %d matrix dump residual dim2 = %d\n",
	      __FILE__, __LINE__,
	      dim2);
      fprintf(_fp, "permute_diag[] = ");
      for (int i = 0; i < dim2; i++) {
	fprintf(_fp, "%d ", permute_diag[i]);
      }
      fprintf(_fp, "\n");
      fprintf(stderr, "%s %d : n = %d\n", __FILE__, __LINE__, 0);
      for (int i = 0; i < dim2; i++) {
	fprintf(stderr, "%d : ", i);
	for (int j = 0; j < dim2; j++) {
	  printscalar<T>(stderr, diag(i, j));
	  fprintf(stderr, " ");
	}
	fprintf(stderr, "\n");
      }
    }
    for (int i = 0; i < dim2; i++) {
      if (_verbose) {
	fprintf(_fp, "%d : ", i);
	printscalar<T>(_fp, diag_origs(i, i));
      	fprintf(_fp, " : ");
	printscalar<T>(_fp, diag(i, i));
	fprintf(_fp, " : ");
      } // if (_verbose)
      for (int j = 0; j <= i; j++) { // lower
	T tmpx(0.0);
	if (_isSymmetric) {
	  for (int k = 0; k < j; k++) {
	    tmpx += diag(k, i) * diag(k, j) / diag(k, k);
	  }
	  if (j < i) {
	    tmpx +=  diag(j, i) / diag(j, j);
	  }
	  else {
	    tmpx +=  _one / diag(j, j);
	  }
	}
	else {
	  for (int k = 0; k < j; k++) {
	    tmpx += diag(i, k) * diag(k, j) / diag(k, k);
	  }
	  if (j < i) {
	    tmpx +=  diag(i, j) / diag(j, j);
	  }
	  else {
	    tmpx +=  _one / diag(j, j);
	  }
	}
	const T tmpy = diag_origs(permute_diag[i], permute_diag[j]) - tmpx;
	const double tmpabs = blas_abs<T, double>(tmpy);
	const double eps_machine = machine_epsilon<double, double>();
	if (_verbose && (tmpabs > eps_machine)) {
	    fprintf(_fp, "%d : ", j);
	    printscalar<T>(_fp, tmpy);
	}
      }
      if (_verbose) {
	fprintf(_fp, "\n");
      }
    }
#endif
    _num_null[0] = nn0;
    *nopd += fop;
    for (int i = 0; i < nn0; i++) {
      //      _null_lists_local[0].push_back(i + (dim2 - nn0));
      // nullifying rows and columns
      for (int i = (dim2 - nn0); i < dim2; i++) {
	for (int j = 0; j < dim2; j++) {
	  diag(i, j) = _zero;
	  diag(j, i) = _zero;
	}
      }
    }
    n0_total += nn0;
  }
  for (int n = 1; n < _nfront; n++) {
    ColumnMatrix<T> &diag1 = _diag_blocks[n - 1];
    vector<int> i0;
    dim1 = dim2;
    dim2 = _p_front[n + 1] - _p_front[n];
#ifdef SPARSE_OFFDIAG
    upper.init(dim1, dim2, upper_work.addrCoefs(), false);
    lower.init(dim1, dim2, lower_work.addrCoefs(), false);
#else
    ColumnMatrix<T> &upper = _upper_blocks[n]; // 27 Nov.2015
    ColumnMatrix<T> &lower = _lower_blocks[n]; // 27 Nov.2015
    upper.init(dim1, dim2);
    lower.init(dim1, dim2);
#endif
    //  ColumnMatrix<T> upper, lower
    // permute_diag[] is defined in the previous step n-1
    // generate permute_offdiag[] from previous permute_diag[]
    GenPermuteOffdiag(_p_front[n - 1], _p_front[n], _p_front[n + 1],
		      _ptRows, _p_diag, _indCols,
		      permute_diag, permute_diag_inv,
		      permute_offdiag, permute_offdiag_inv,
		      i0);
    upper.ZeroClear();
    if (_isSymmetric) {
      for (int i = _p_front[n]; i < _p_front[n + 1]; i++) {
	const int jj = permute_offdiag_inv[i - _p_front[n]];
	  for (int k = _ptRows[i]; k < _p_diag[i]; k++) {
	    // access lower block to create upper block
	  const int ii = permute_diag_inv[_indCols[k] - _p_front[n - 1]];
	  upper(ii, jj) = _coef[_indVals[k]];
	}
      }
    }
    else {
      lower.ZeroClear();
      for (int i = _p_front[n - 1]; i < _p_front[n]; i++) {
	const int ii = permute_diag_inv[i - _p_front[n - 1]];
	for (int k = _p_upper[i]; k < _ptRows[i + 1]; k++) {
	  const int jj = permute_offdiag_inv[_indCols[k] - _p_front[n]];
	  upper(ii, jj) = _coef[_indVals[k]];
	}
      }
      for (int i = _p_front[n]; i < _p_front[n + 1]; i++) {
	const int jj = permute_offdiag_inv[i - _p_front[n]];
	for (int k = _ptRows[i]; k < _p_diag[i]; k++) {
	  // same as upper() of _isSymmetric
	  const int ii = permute_diag_inv[_indCols[k] - _p_front[n - 1]];
	  lower(ii, jj) = _coef[_indVals[k]];
	}
      }
    }
    full_fw_multiprofile<T>(false, dim1, _num_null[n - 1], dim2,
			    diag1.addrCoefs(), dim1,
			    upper.addrCoefs(), dim1, i0, &fop);
    if (_isSymmetric) { // reduce arithmeic by setting T(0) for 0 <=i< i0[j]
      for (int j = 0; j < dim2; j++) {
	for (int i = 0; i < dim1; i++) {
	  lower(i, j) = upper(i, j) * diag1(i, i);
	}
      }
    }
    else {
      full_fw_multiprofile<T>(true, dim1, _num_null[n - 1], dim2,
			      diag1.addrCoefs(), dim1,
			      lower.addrCoefs(), dim1, i0, &fop);
      for (int j = 0; j < dim2; j++) {
	for (int i = 0; i < dim1; i++) {
	  lower(i, j) *= diag1(i, i);
	}
      }
    }
    ColumnMatrix<T> &diag = _diag_blocks[n];
    // alpha = -1, beta = 0
    SparseSchur<T>(_isSymmetric, dim2, dim1, i0, upper, lower, diag, &fop);
    for (int i = _p_front[n]; i < _p_front[n + 1]; i++) {
      for (int k = _p_diag[i]; k < _p_upper[i]; k++) { 
	const int ii = permute_offdiag_inv[i - _p_front[n]];
	const int jj = permute_offdiag_inv[_indCols[k] - _p_front[n]];
	diag(ii, jj) += _coef[_indVals[k]];
      }
    }
#ifndef SPARSE_OFFDIAG
    for (int j = 0; j < dim2; j++) {
      for (int i = 0; i < dim1; i++) {
	upper(i, j) *= diag1(i, i);
      }
    }
#endif
#ifdef DEBUG_FACTORIZATION
    diag_origs.init(dim2, dim2, diag_orig, false);
    for (int j = 0; j < dim2; j++) {
      for (int i = 0; i < dim2; i++) {
	diag_origs(i, j) = diag(i, j);
      }
    }
#endif
    if (_isSymmetric) {
      bool flag;
      flag = full_ldlt_permute<T, U>(&nn0, n0, dim2, diag.addrCoefs(), dim2,
				     pivot, permute_work, eps_pivot, &fop);
    }
    else {
      bool flag;
      flag = full_ldu_permute<T, U>(&nn0, n0, dim2, diag.addrCoefs(), dim2,
				    pivot, permute_work, eps_pivot, &fop);
    }

    _num_null[n] = nn0;
    permute_diag = &_permute[0] + _p_front[n];
        // generate permute_diag from pemrute_work + perumte_upper
    for (int i = 0; i < dim2; i++) {
      permute_diag[i] = permute_offdiag[permute_work[i]];
    }
#ifndef SPARSE_OFFDIAG
    for (int j = 0; j < dim2; j++) {
      for (int i = 0; i < dim1; i++) {
	permute_vals[i + j * dim1] = upper(i, permute_work[j]);
      }
    }
    for (int j = 0; j < dim2; j++) {
      for (int i = 0; i < dim1; i++) {
	upper(i, j) = permute_vals[i + j * dim1];
      }
    }
    for (int j = 0; j < dim2; j++) {
      for (int i = 0; i < dim1; i++) {
	permute_vals[i + j * dim1] = lower(i, permute_work[j]);
      }
    }
    for (int j = 0; j < dim2; j++) {
      for (int i = 0; i < dim1; i++) {
	lower(i, j) = permute_vals[i + j * dim1];
      }
    }
#endif
#ifdef STORE_WHOLE_FACTORIZED
    for (int j = 0; j < dim2; j++) {
      for (int i = 0; i < dim2; i++) {
	_factorized_whole(_p_front[n] + i, _p_front[n] + j) = diag(i, j);
      }
    }
    if (_isSymmetric) {
      for (int j = 0; j < dim2; j++) {
	for (int i = 0; i < dim1; i++) {
	  _factorized_whole(_p_front[n] + j, _p_front[n - 1] + i) =
	    lower(i, permute_work[j]);
	}
      }
    }
    else {
      for (int j = 0; j < dim2; j++) {
	for (int i = 0; i < dim1; i++) {
	  _factorized_whole(_p_front[n - 1] + i, _p_front[n] + j) =
	    upper(i, permute_work[j]) * diag1(i, i);
	  _factorized_whole(_p_front[n] + j, _p_front[n - 1] + i) =
	    lower(i, permute_work[j]);
	}
      }
    }
#endif
    // debug
#ifdef DEBUG_FACTORIZATION
    if (_verbose) {
      fprintf(_fp, "%s %d matrix dump residual dim2 = %d\n",
	      __FILE__, __LINE__,
	      dim2);
      fprintf(_fp, "permute_diag[] = ");
      for (int i = 0; i < dim2; i++) {
	fprintf(_fp, "%d ", permute_diag[i]);
      }
      fprintf(_fp, "\n");
      fprintf(stderr, "%s %d : n = %d\n", __FILE__, __LINE__, n);
      for (int i = 0; i < dim2; i++) {
	fprintf(stderr, "%d : ", i);
	for (int j = 0; j < dim2; j++) {
	  printscalar<T>(stderr, diag(i, j));
	  fprintf(stderr, " ");
	}
	fprintf(stderr, "\n");
      }
    } // if (_verbose)
     for (int i = 0; i < dim2; i++) {
       if (_verbose) {
	 fprintf(_fp, "%d : ", i);
	 printscalar<T>(_fp, diag_origs(i, i));
	 fprintf(_fp, " : ");
	 printscalar<T>(_fp, diag(i, i));
	 fprintf(_fp, " : ");
       } // if (_verbose)
       for (int j = 0; j <= i; j++) { // lower
	T tmpx(0.0);
	if (_isSymmetric) {
	  for (int k = 0; k < j; k++) {
	    tmpx += diag(k, i) * diag(k, j) / diag(k, k);
	  }
	  if (j < i) {
	    tmpx += diag(j, i) / diag(j, j);
	  }
	  else {
	    tmpx += _one / diag(j, j);
	  }
	}
	else {
	  for (int k = 0; k < j; k++) {
	    tmpx += diag(i, k) * diag(k, j) / diag(k, k);
	  }
	  if (j < i) {
	    tmpx += diag(i, j) / diag(j, j);
	  }
	  else {
	    tmpx += _one / diag(j, j);
	  }
	}
	if (blas_abs<T, U>(diag_origs(permute_work[i], permute_work[j]) - tmpx) >
	    machine_epsilon<T, U>()) {
	  if (_verbose) {
	    fprintf(_fp, "%d : ", j);
	    printscalar<T>(_fp,
			   (diag_origs(permute_work[i], permute_work[j])
			    - tmpx));
	  } // if (_verbose)
	}
      }
      if (_verbose) {
	fprintf(_fp, "\n");
      }
    }
#endif
    // debug
    *nopd += fop;
    for (int i = 0; i < nn0; i++) {
      // nullifying rows and columns
      for (int i = (dim2 - nn0); i < dim2; i++) {
	for (int j = 0; j < dim2; j++) {
	  diag(i, j) = _zero;
	  diag(j, i) = _zero;
	}
      }
    }
    n0_total += nn0;
  } // loop : n
  if ((n0_total > 0) && (dim_aug_kern > 0)) {
    int count = 0;
    bool flag = false;
    num_null_aug.clear();
    num_null_aug.resize(_nfront, 0);
    for (int n = (_nfront - 1); n >= 0; n--) {
      dim1 = _p_front[n + 1] - _p_front[n];
      for (int k = (dim1 - _num_null[n] - 1); k >= 0; k--) {
	num_null_aug[n]++;
	count++;
	if (count >= dim_aug_kern) {
	  flag = true;
	  break;
	}
      }
      if (flag) {
	break;
      }
    }
    if (flag == false) {
      if (_verbose) {
	fprintf(_fp, "%s %d : dim_aug_kern = %d > dim = %d?\n",
		__FILE__, __LINE__, dim_aug_kern, _dim);
      }
    }
    for (int n = (_nfront - 1); n >= 0; n--) {
      diag_block_save[n].free();
      if (num_null_aug[n] > 0) {
	dim1 = _p_front[n + 1] - _p_front[n];
	diag_block_save[n].init(dim1, dim1);
	diag_block_save[n].copy(_diag_blocks[n]);
	_num_null[n] += num_null_aug[n]; //
	for (int j = 0; j < num_null_aug[n]; j++) {
	  const int jj = dim1 - _num_null[n] + j;
	  for (int i = 0; i < dim1; i++) {
	    _diag_blocks[n](i, jj) = _zero;
	    _diag_blocks[n](jj, i) = _zero;
	  }
	}
      }
    }  // loop : n
  } // if ((n0_total > 0) && kernel_detection) {
  _permute_ginv.resize(_dim);
  for (int n = 0; n < _nfront; n++) {
    const int offset = _p_front[n];
    for (int i = offset; i < _p_front[n + 1]; i++) {
      _permute_ginv[_permute[i] + offset] = i;
    }
  }
  //#ifdef SPARSE_OFFDIAG   
  //  delete [] upper_work;
  //  delete [] lower_work;
  //#else
  //  delete [] permute_vals;
  //#endif
  delete [] permute_offdiag;
  delete [] permute_diag_inv;
  delete [] permute_offdiag_inv;
  delete [] permute_work;
  //#ifdef DEBUG_FACTORIZATION
  //  delete [] diag_orig;
  //#endif

}

template
void TridiagBlockMatrix<double>::
TridiagNumericFact(double *pivot,
		   const double eps_pivot,
		   const int dim_aug_kern,
		   ColumnMatrix<double> *diag_block_save,
		   vector<int> &num_null_aug,
		   double *nopd);

template
void TridiagBlockMatrix<quadruple>::
TridiagNumericFact(double *pivot,
		   const double eps_pivot,
		   const int dim_aug_kern,
		   ColumnMatrix<quadruple> *diag_block_save,
		   vector<int> &num_null_aug,
		   double *nopd);

template
void TridiagBlockMatrix<complex<double>, double>::
TridiagNumericFact(double *pivot,
		   const double eps_pivot,
		   const int dim_aug_kern,
		   ColumnMatrix<complex<double> > *diag_block_save,
		   vector<int> &num_null_aug,
		   double *nopd);

template
void TridiagBlockMatrix<complex<quadruple>, quadruple>::
TridiagNumericFact(double *pivot,
		   const double eps_pivot,
		   const int dim_aug_kern,
		   ColumnMatrix<complex<quadruple> > *diag_block_save,
		   vector<int> &num_null_aug,
		   double *nopd);
//

template<typename T, typename U>
void TridiagBlockMatrix<T, U>::ComputeSchur(const int dim_,
					    int* color_mask,
					    const int ncol_,
					    const int *ptrow1,
					    const int *indcols1,
					    const int *indvals1,
					    const int *indvals2, 
					    const T *coef,
					    const int size_b1,
					    SquareBlockMatrix<T> &local_s,
					    double *nopd,
					    elapsed_t *tt)
{
  //input : CSR data of upper block : ptrow1, indcols1, indvals1, indvals2, coef
  //output: Schur complment : -c^T a^-1 b
  //      : permutation of upper block : old2new_j  
  ColumnMatrix<T> b, c, y, z;
  vector<T> dscale;
  //  const int nrow = _dim;
  vector<int> i0;
  vector<int> permute_diag_inv, permute_upper, permute_upper_inv;
  //  T *s;
  int ncol;
  vector<int> map_eqn, remap_eqn;

  get_realtime(&tt[0]); // entering
  map_eqn.resize(dim_);
  remap_eqn.resize(dim_);

  int jtmp = _dim;
  int itmp = 0;
  //  fprintf(_fp, "%s %d : graph_colors = %d %d x %d _nscol = %d\n",
  //	  __FILE__, __LINE__, _color, _dim, ncol_, _nscol);
  for (int i = 0; i < dim_; i++) {
    if ((color_mask[i] == _color) || (color_mask[i] == (-_color)) ||
	((_color == 1) && color_mask[i] == 0)) {
      map_eqn[itmp++] = i;
    }
    else {
      map_eqn[jtmp++] = i;
    }
  }
  for (int i = 0; i < dim_; i++) {
    remap_eqn[map_eqn[i]] = i;
  }
  const int nrow = _dim;
  GenPermuteUpper(nrow, remap_eqn, ncol_, ptrow1, indcols1, _nfront, _p_front,
		  _new2old, _permute, permute_diag_inv,
		  permute_upper,
		  permute_upper_inv, i0, _verbose, _fp);
  ncol = 0;
  for (vector<int>::const_iterator it = i0.begin(); it != i0.end(); ++it) {
    ncol++;
    if (*it >= nrow) {
      break;
    }
  }
#ifdef PERMUTE_UPPER
  ColumnMatrix<T> s(ncol, ncol); 
#endif
  b.init(nrow, ncol);
  c.init(nrow, ncol);
  dscale.resize(nrow);
  b.ZeroClear();
  FillBlockSparse<T>(coef, nrow, map_eqn, ptrow1, indcols1, indvals1,
		     permute_diag_inv, permute_upper_inv, b);
  if (!_isSymmetric) {
    c.ZeroClear();
    FillBlockSparse<T>(coef, nrow, map_eqn, ptrow1, indcols1, indvals2,
		       permute_diag_inv, permute_upper_inv, c);
  }
  get_realtime(&tt[1]); // after filling enteries
// A_12' = A_11^-1 A_12   A_21'^T = A_11^-T A_21^T 
//[L_11      ][D_11     ][U_11 A_12'] forward : x_1 = L_11^-1 b_1
//[A_21' L_22][     D_22][     U_22 ]           x_2 = L_22^-1 (b_2 - A_21' x_1)
//                                    backward: x_1 = U_11^-T b_1
//                                            x_2 = U_22^-1 (b_2 - A_12'^T x_1)
// symmetric => A_21' = A_12'^T 
// unsymmetric A_21'^T is stored ==> DGEMM('T','N',...,(tri%a12) or (tri%a21),
  if (_nscol > 0) {
    y.init(_nscol, ncol);
    z.init(_nscol, ncol);
    for (int j = 0; j < ncol; j++) {
      for (int i = 0; i < _nscol; i++) {
	y(i, j) = b(_list_schur[i], j);
	// nullifying
	b(_list_schur[i], j) = _zero;
      }
    }
#if 0
    fprintf(_fp, "y nscol = %d ncol = %d\n", _nscol, ncol);
    for (int i = 0; i < _nscol; i++) {
      for (int j = 0; j < ncol; j++) {
	printscalar<T>(stderr, y(i, j));
      }
      fprintf(_fp, "\n");
    }
#endif
    // alpha = -1; beta = 1;
    if (_isSymmetric) {
      // y -= A_12 ^T b
      blas_gemm<T>(CblasTrans, CblasNoTrans, _nscol, ncol, nrow, _none,
		   _a12.addrCoefs(), nrow, b.addrCoefs(), nrow, _one,
		   y.addrCoefs(), _nscol);
    }
    else {
      // y -= A_21 ^T b
      blas_gemm<T>(CblasTrans, CblasNoTrans, _nscol, ncol, nrow, _none,
		   _a21.addrCoefs(), nrow, b.addrCoefs(), nrow, _one,
		   y.addrCoefs(), _nscol);
    }
#if 0
    fprintf(_fp, "y nscol = %d ncol = %d\n", _nscol, ncol);
    for (int i = 0; i < _nscol; i++) {
      for (int j = 0; j < ncol; j++) {
	printscalar<T>(_fp, y(i, j));
      }
      fprintf(_fp, "\n");
    }
#endif
  } // if (_nscol > 0)
  ForwardUpper(false, ncol, b, i0, dscale); // normal
#if 0
  if (_verbose) {
    fprintf(_fp, "%s %d singular part nscol = %d\n",
	    __FILE__, __LINE__, _nscol);
  }
#endif
  if (_nscol > 0) {
    // alpha = 1.0
    blas_trsm<T>(CblasLeft, CblasLower, CblasNoTrans, CblasUnit, _nscol, ncol,
		 _one, _s22.addrCoefs(), _nscol, y.addrCoefs(), _nscol);
    if (_isSymmetric) {
      for (int i = 0; i < _nscol; i++) {         // exclude _2 entries
	dscale[_list_schur[i]] = _zero;
      }
      for (int j = 0; j < ncol; j++) {
	for (int i = 0; i < _nscol; i++) {
	  z(i, j) = y(i, j) * _s22(i, i);
	}
      }
#if 0
      if (_verbose) {
	fprintf(_fp, "%s %d : s22() dim = %d\n", __FILE__, __LINE__, _nscol);
	for (int i = 0; i < _nscol; i++) {
	  fprintf(_fp, "%d ", i);
	  for (int j = 0; j < _nscol; j++) {
	    printscalar<T>(_fp, _s22(i, j));
	  }
	  fprintf(_fp, "\n");
	}
	fprintf(_fp, "y nscol = %d ncol = %d\n", _nscol, ncol);
	for (int i = 0; i < _nscol; i++) {
	  for (int j = 0; j < ncol; j++) {
	    printscalar<T>(_fp, y(i, j));
	  }
	  fprintf(_fp, "\n");
	}
	fprintf(_fp, "z nscol = %d ncol = %d\n", _nscol, ncol);
	for (int i = 0; i < _nscol; i++) {
	  for (int j = 0; j < ncol; j++) {
	    printscalar<T>(_fp, z(i, j));
	  }
	  fprintf(_fp, "\n");
	}
      } // if (_verbose)
#endif
    }
  }
  if (_isSymmetric) {
    // ddtimesu
    for (int j = 0; j < ncol; j++) {
      for (int i = 0; i < i0[j]; i++) {
	c(i, j) = _zero;
      }
      for (int i = i0[j]; i < nrow; i++) {
       	c(i, j) = dscale[i] * b(i, j);
      }
      // fop += (nrow-j0(j)-1)
    }
  }
  else {   // if (isSymmetric)
    if (_nscol > 0) {
      for (int j = 0; j < ncol; j++) {
	for (int i = 0; i < _nscol; i++) {
	  z(i, j) = c(_list_schur[i], j);
	  // nullifying : easier than changing entries in dscle[]
	  c(_list_schur[i], j) = _zero;
	}
      }
      // alpha = -1; beta = 1;
      // z -= A_12 ^T c
      blas_gemm<T>(CblasTrans, CblasNoTrans, _nscol, ncol, nrow, _none,
		   _a12.addrCoefs(), nrow, c.addrCoefs(), nrow, _one,
		   z.addrCoefs(), _nscol);
    }
    ForwardUpper(true, ncol, c, i0, dscale); // Transposed
    if (_nscol > 0) {
      // alpha = 1.0
      blas_trsm<T>(CblasLeft, CblasUpper, CblasTrans, CblasUnit,
		   _nscol, ncol, _one,
		   _s22.addrCoefs(), _nscol, z.addrCoefs(), _nscol);
      for (int i = 0; i < _nscol; i++) {
	dscale[_list_schur[i]] = _zero;
      }
      for (int j = 0; j < ncol; j++) {
	for (int i = 0; i < _nscol; i++) {
	  z(i, j) *= _s22(i, i);
	}
      }
    }
    // ddtimesu
    for (int j = 0; j < ncol; j++) {
      for (int i = 0; i < i0[j]; i++) {
	c(i, j) = _zero;
      }
      for (int i = i0[j]; i < nrow; i++) {
	c(i, j) *= dscale[i];
      }
      // fop += (nrow-j0(j)-1)
    }
  } // if (_isSymmetric)
  //nops += fop;
  double fop;
  get_realtime(&tt[2]); // after tridiag forward
#ifdef PERMUTE_UPPER
  // s = b * c : alpha = 1.0; beta = 0.0
  if (_isSymmetric) {
    SchurProfileSym<T>(nrow, ncol, i0, b, c, s.addrCoefs(), size_b1, &fop);
  }
  else {
    SchurProfileUnSym<T>(nrow, ncol, i0, b, c, s.addrCoefs(), size_b1, &fop);
  }
  if (_nscol > 0) {
    // s += z y : alpha = 1.0; beta = 1.0;
    blas_gemm<T>(CblasTrans, CblasNoTrans, ncol, ncol, _nscol, _one,
		 z.addrCoefs(), _nscol, y.addrCoefs(), _nscol, _one,
		 s.addrCoefs(), ncol);
  }
  get_realtime(&tt[3]); // after sparse dgemm
  // nops = tri%nop
  if (_color == 1) {
    local_s.ZeroClear();
    if (_isSymmetric) {
      for (int i = 0; i < ncol; i++) {
	const int ii = permute_upper[i];
	for (int j = i; j < ncol; j++) {
	  const int jj = permute_upper[j];
	  const int iii = (ii <= jj ? ii : jj);
	  const int jjj = (ii <= jj ? jj : ii);
	  local_s(iii, jjj) = s(i, j);
     	}
      }
    }
    else {
      for (int i = 0; i < ncol; i++) {
     	const int ii = permute_upper[i];
     	for (int j = 0; j < ncol; j++) {
     	  const int jj = permute_upper[j];
     	  local_s(ii, jj) = s(i, j);
	}
      }
    }
    //    local_s.copyFromArrayPermute(s, ncol, &permute_upper_inv[0]);
  }
  else {
    if (_isSymmetric) {
      for (int i = 0; i < ncol; i++) {
	const int ii = permute_upper[i];
	for (int j = i; j < ncol; j++) {
	  const int jj = permute_upper[j];
	  const int iii = (ii <= jj ? ii : jj);
	  const int jjj = (ii <= jj ? jj : ii);
	  local_s(iii, jjj) += s(i, +j);
	}
      }
    }
    else {
      for (int i = 0; i < ncol; i++) {
	const int ii = permute_upper[i];
	for (int j = 0; j < ncol; j++) {
	  const int jj = permute_upper[j];
	  local_s(ii, jj) += s(i, j);
	}
      }
    }
    //    local_s.addFromArrayPermute(s, ncol, &permute_upper_inv[0]);
  } // if (_color == 1)
#else
 // s = b * c : alpha = 1.0; beta = 0.0
  if (_isSymmetric) {
    for (int j = 0; j < local_s.num_blocks(); j++) {
      for (int i = 0; i < j; i++) {
	const int nnrow = local_s.nrowBlock(i, j);
	const int nncol = local_s.ncolBlock(i, j);
 	const int ishift = local_s.IndexBlock(i) * nrow;
	const int jshift = local_s.IndexBlock(j) * nrow;
	blas_gemm<T>(CblasTrans, CblasNoTrans,
		     nnrow, nncol, nrow, _one,
		     c.addrCoefs() + ishift, nrow,
		     b.addrCoefs() + jshift, nrow,
		     (_color == 1 ? _zero : _one),
		     local_s.addrCoefBlock(i, j), nnrow);
	if (_nscol > 0) {
 	  const int ishift1 = local_s.IndexBlock(i) * _nscol;
	  const int jshift1 = local_s.IndexBlock(j) * _nscol;
	  // s += z y : alpha = 1.0; beta = 1.0;
	  blas_gemm<T>(CblasTrans, CblasNoTrans, nnrow, nncol, _nscol, _one,
		       z.addrCoefs() + ishift1, _nscol,
		       y.addrCoefs() + jshift1, _nscol, _one,
		       local_s.addrCoefBlock(i, j), nnrow);
	}
      }
      {
      	const int nnrow = local_s.nrowBlock(j, j);
	const int jshift = local_s.IndexBlock(j) * nrow;
	C_gemm_symm<T>(nnrow, nrow, _one,
		       c.addrCoefs() + jshift, nrow,
		       b.addrCoefs() + jshift, nrow,
		       (_color == 1 ? _zero : _one),
		       local_s.addrCoefBlock(j, j), nnrow);
	if (_nscol > 0) {
	  const int jshift1 = local_s.IndexBlock(j) * _nscol;
	  // s += z y : alpha = 1.0; beta = 1.0;
	  C_gemm_symm<T>(nnrow,  _nscol, _one,
			 z.addrCoefs() + jshift1, _nscol,
			 y.addrCoefs() + jshift1, _nscol, _one,
			 local_s.addrCoefBlock(j, j), nnrow);
	}
      }	
    } // loop : j
  } // if (_isSymmetric)
  else {
    for (int j = 0; j < local_s.num_blocks(); j++) {
      for (int i = 0; i <= j; i++) {
	const int nnrow = local_s.nrowBlock(i, j);
	const int nncol = local_s.ncolBlock(i, j);
	const int ishift = local_s.IndexBlock(i) * nrow;
	const int jshift = local_s.IndexBlock(j) * nrow;
	blas_gemm<T>(CblasTrans, CblasNoTrans,
		     nnrow, nncol, nrow, _one,
		     c.addrCoefs() + ishift, nrow,
		     b.addrCoefs() + jshift, nrow,
		     (_color == 1 ? _zero : _one),
		     local_s.addrCoefBlock(i, j), nnrow);
	if (_nscol > 0) {
	  const int ishift1 = local_s.IndexBlock(i) * _nscol;
	  const int jshift1 = local_s.IndexBlock(j) * _nscol;
	  // s += z y : alpha = 1.0; beta = 1.0;
	  blas_gemm<T>(CblasTrans, CblasNoTrans, nnrow, nncol, _nscol, _one,
		       z.addrCoefs() + ishift1, _nscol,
		       y.addrCoefs() + jshift1, _nscol, _one,
		       local_s.addrCoefBlock(i, j), nnrow);
	}
      } // lower part of local_s stores data in transposed way
      for (int i = (j + 1); i < local_s.num_blocks(); i++) {
	const int nnrow = local_s.nrowBlock(i, j);
	const int nncol = local_s.ncolBlock(i, j);
	const int ishift = local_s.IndexBlock(i) * nrow;
	const int jshift = local_s.IndexBlock(j) * nrow;
	blas_gemm<T>(CblasTrans, CblasNoTrans,
		     nnrow, nncol, nrow, _one,
		     b.addrCoefs() + jshift, nrow,
		     c.addrCoefs() + ishift, nrow,
		     (_color == 1 ? _zero : _one),
		     local_s.addrCoefBlock(i, j), nnrow);
	if (_nscol > 0) {
	  const int ishift1 = local_s.IndexBlock(i) * _nscol;
	  const int jshift1 = local_s.IndexBlock(j) * _nscol;
	  // s += z y : alpha = 1.0; beta = 1.0;
	  blas_gemm<T>(CblasTrans, CblasNoTrans, nnrow, nncol, _nscol, _one,
		       y.addrCoefs() + jshift1, _nscol, 
		       z.addrCoefs() + ishift1, _nscol, _one,
		       local_s.addrCoefBlock(i, j), nnrow);
	}
      }
    } 
  } // if (_isSymmetric)
  get_realtime(&tt[3]); // after sparse dgemm
#endif
#if 0
  if (_verbose) {
    fprintf(_fp, "%s %d : local_s : ncol_ = %d color_max = %d\n",
	    __FILE__, __LINE__, ncol_, _color_max);
    
    fprintf(_fp, "permute_diag_inv : ");
    vector<int>::const_iterator it, jt;
    it = permute_diag_inv.begin();
    for ( ; it != permute_diag_inv.end(); ++it) {
      fprintf(_fp, "%d ", *it);
    }
    fprintf(_fp, "\n");
    fprintf(_fp, "permute_upper_inv : ");
    it = permute_upper_inv.begin();
    jt = i0.begin();
    for ( ; it != permute_upper_inv.end(); ++it, ++jt) {
      fprintf(_fp, "%d:%d ", *it, *jt);
    }
    fprintf(_fp, "\n");
    
    if (_color == _color_max) {
      for (int i = 0; i < ncol_; i++) {
	fprintf(_fp, "%d : ", i);
	for (int j = i; j < ncol_; j++) {
	  printscalar<T>(_fp, local_s(i, j));
	}
	fprintf(_fp, "\n");
      }
    }
  } // if (_verbose)
#endif
  b.free();
  c.free();
  y.free();
  z.free();
  dscale.clear();
  permute_diag_inv.clear();
  permute_upper_inv.clear();
#ifdef PERMUTE_UPPER
  s.free(); //  delete [] s;
#endif
  get_realtime(&tt[4]); // end of the routine
}

template
void TridiagBlockMatrix<double>::
ComputeSchur(const int dim_,
	     int* color_mask,
	     const int ncol,
	     const int *ptrow1,
	     const int *indcols1,
	     const int *indvals1,
	     const int *indvals2, 
	     const double *coef,
	     const int size_b1,
	     SquareBlockMatrix<double> &local_s,
	     double *nopd,
	     elapsed_t *tt);

template
void TridiagBlockMatrix<quadruple>::
ComputeSchur(const int dim_,
	     int* color_mask,
	     const int ncol,
	     const int *ptrow1,
	     const int *indcols1,
	     const int *indvals1,
	     const int *indvals2, 
	     const quadruple *coef,
	     const int size_b1,
	     SquareBlockMatrix<quadruple> &local_s,
	     double *nopd,
	     elapsed_t *tt);

template
void TridiagBlockMatrix<complex<double>, double>::
ComputeSchur(const int dim_,
	     int* color_mask,
	     const int ncol,
	     const int *ptrow1, const int *indcols1,
	     const int *indvals1, const int *indvals2, 
	     const complex<double> *coef,
	     const int size_b1,
	     SquareBlockMatrix<complex<double> > &local_s,
	     double *nopd,
	     elapsed_t *tt);

template
void TridiagBlockMatrix<complex<quadruple>, quadruple>::
ComputeSchur(const int dim_,
	     int* color_mask,
	     const int ncol,
	     const int *ptrow1, const int *indcols1,
	     const int *indvals1, const int *indvals2, 
	     const complex<quadruple> *coef,
	     const int size_b1,
	     SquareBlockMatrix<complex<quadruple> > &local_s,
	     double *nopd,
	     elapsed_t *tt);
//

template<typename T, typename U>
void TridiagBlockMatrix<T, U>::SolveMulti(const bool flag_new2old,
					  const bool isTrans,
					  const int nrhs, ColumnMatrix<T>& x,
					  const int nscol_)
{
  ColumnMatrix<T> xx, y, zz;
  double fop;
  const int nscol = (nscol_ == (-1) ? _nscol : nscol_);
  xx.init(_dim, nrhs);  
  if (flag_new2old) {
    for (int m = 0; m < nrhs; m++) {
      for (int i = 0; i < _dim; i++) {
	xx(_permute_ginv[i], m) = x(_new2old[i], m);
      }
    }
  }
  else {
    for (int m = 0; m < nrhs; m++) {
      for (int i = 0; i < _dim; i++) {
	xx(_permute_ginv[i], m) = x(i, m);
      }
    }
  }
  // x(, m) is reorderd to follow the permutation generated by
  // TridiagBlockFactorization
  if (nscol > 0) {
    y.init(nscol, nrhs);
    for (int m = 0; m < nrhs; m++) {
      for (int i = 0; i < nscol; i++) {
	y(i, m) = xx(_list_schur[i], m);
	xx(_list_schur[i], m) = _zero;     		// nullifying
      }
    }
    // alpha = (-1.0); beta = 1.0
    if (_isSymmetric) {
      blas_gemm<T>(CblasTrans, CblasNoTrans, nscol, nrhs, _dim, _none,
		   _a12.addrCoefs(), _dim,
		   xx.addrCoefs(), _dim, _one, y.addrCoefs(), nscol);
    }
    else {
      if (isTrans) {
	blas_gemm<T>(CblasTrans, CblasNoTrans, nscol, nrhs, _dim, _none,
		     _a12.addrCoefs(), _dim,
		     xx.addrCoefs(), _dim, _one, y.addrCoefs(), nscol);
      }
      else {
	blas_gemm<T>(CblasTrans, CblasNoTrans, nscol, nrhs, _dim, _none,
		     _a21.addrCoefs(), _dim,
		     xx.addrCoefs(), _dim, _one, y.addrCoefs(), nscol);
      }
    }
    blas_trsm<T>(CblasLeft,
		 isTrans ? CblasUpper : CblasLower,
		 isTrans ? CblasTrans : CblasNoTrans,
		 CblasUnit, nscol, nrhs,
		 _one, _s22.addrCoefs(), nscol, y.addrCoefs(), nscol);
    for (int m = 0; m < nrhs; m++) {
      for (int i = 0; i < nscol; i++) {
	y(i, m) *= _s22(i, i);
      }
    }
    blas_trsm<T>(CblasLeft,
		 isTrans ? CblasLower : CblasUpper,
		 isTrans ? CblasTrans : CblasNoTrans,
		 CblasUnit, nscol, nrhs,
		 _one, _s22.addrCoefs(), nscol, y.addrCoefs(), nscol);
  } // if (nscol > 0)
  zz.init(_maxdim, nrhs);
  // forward substitution
  for (int n = 0; n < (_nfront - 1); n++) {
    // copy to working array zz
    const int dim1 = _p_front[n + 1] - _p_front[n];
    for (int m = 0; m < nrhs; m++) {
      for (int i = 0; i < dim1; i++) {
	zz(i, m) = xx(_p_front[n] + i, m);
      }
    }
#ifdef SPARSE_OFFDIAG
    full_fwbw_multi<T>(isTrans, dim1, _num_null[n],
		       _diag_blocks[n].addrCoefs(), dim1,
		       nrhs,
		       zz.addrCoefs(), _maxdim);
    // sparse matrix * dense matrix
    if (isTrans) {
      for (int m = 0; m < nrhs; m++) {
	for (int i = _p_front[n]; i < _p_front[n + 1]; i++) {
	  const int ii = _permute_ginv[i] - _p_front[n];
	  for (int k = _p_upper[i];  k < _ptRows[i + 1]; k++) {
	    const int jj = _permute_ginv[_indCols[k]];
	    xx(jj, m) -= _coef[_indVals[k]] * zz(ii, m);
	  }
	}
      }
    }
    else {
      for (int m = 0; m < nrhs; m++) {
	for (int i = _p_front[n + 1]; i < _p_front[n + 2]; i++) {
	  const int ii = _permute_ginv[i];
	  for (int k = _ptRows[i];  k < _p_diag[i]; k++) {
	    const int jj = _permute_ginv[_indCols[k]] - _p_front[n];
	    xx(ii, m) -= _coef[_indVals[k]] * zz(jj, m);
	  }
	}
      }
    } // if (isTrans)
#else
    full_fw_multi<T>(isTrans, dim1, _num_null[n],
		     _diag_blocks[n].addrCoefs(), dim1,
		     nrhs,
		     zz.addrCoefs(), _maxdim, &fop);
    int dim2 = _p_front[n + 2] - _p_front[n + 1];
    if (isTrans) {
      ColumnMatrix<T> &upper = _upper_blocks[n + 1];
      blas_gemm<T>(CblasTrans, CblasNoTrans, dim2, nrhs, dim1,
		   _none, upper.addrCoefs(), upper.nbRows(),
		   zz.addrCoefs(), _maxdim,
		   _one,
		   xx.addrCoefs() + _p_front[n + 1], _dim);
    }
    else {
      ColumnMatrix<T> &lower = _lower_blocks[n + 1];
      blas_gemm<T>(CblasTrans, CblasNoTrans, dim2, nrhs, dim1,
		   _none, lower.addrCoefs(), lower.nbRows(),
		   zz.addrCoefs(), _maxdim,
		   _one,
		   xx.addrCoefs() + _p_front[n + 1], _dim);
    }
    for (int m = 0; m < nrhs; m++) {
      for (int i = 0; i < dim1; i++) {
	xx(_p_front[n] + i, m) = zz(i, m) * _diag_blocks[n](i, i);
      }
    }
#endif
  }
#ifndef SPARSE_OFFDIAG
  {
    const int n = _nfront - 1;
    const int dim1 = _p_front[n + 1] - _p_front[n];
    for (int m = 0; m < nrhs; m++) {
      for (int i = 0; i < dim1; i++) {
	zz(i, m) = xx(_p_front[n] + i, m);
      }
    }
    full_fw_multi<T>(isTrans, dim1, _num_null[n],
		     _diag_blocks[n].addrCoefs(), dim1,
		     nrhs,
		     zz.addrCoefs(), _maxdim, &fop);
    for (int m = 0; m < nrhs; m++) {
      for (int i = 0; i < dim1; i++) {
	xx(_p_front[n] + i, m) = zz(i, m) * _diag_blocks[n](i, i);
      }
    }
  }
#endif
  // backward substitution
  for (int n = (_nfront - 1); n > 0; n--) {
    const int dim1 = _p_front[n + 1] - _p_front[n];
#ifdef SPARSE_OFFDIAG
    full_fwbw_multi<T>(isTrans, dim1, _num_null[n],
		       _diag_blocks[n].addrCoefs(), dim1,
		       nrhs,
		       xx.addrCoefs() + _p_front[n], _dim);
    if (_isSymmetric) {
      // Transposed sparse matrix * dense matrix
      for (int m = 0; m < nrhs; m++) {
	for (int i = _p_front[n]; i < _p_front[n + 1]; i++) {
	  const int ii = _permute_ginv[i];
	  for (int k = _ptRows[i];  k < _p_diag[i]; k++) {
	    const int jj = _permute_ginv[_indCols[k]];
	    xx(jj, m) -= _coef[_indVals[k]] * xx(ii, m);
	  }
	}
      }
    }
    else {
      if (isTrans) {
	for (int m = 0; m < nrhs; m++) {
	  for (int i = _p_front[n]; i < _p_front[n + 1]; i++) {
	    const int ii = _permute_ginv[i];
	    for (int k = _ptRows[i];  k < _p_diag[i]; k++) {
	      const int jj = _permute_ginv[_indCols[k]];
	      xx(jj, m) -= _coef[_indVals[k]] * xx(ii, m);
	    }
	  }
	}
      }
      else {
	for (int m = 0; m < nrhs; m++) {
	  for (int i = _p_front[n - 1]; i < _p_front[n]; i++) {
	    const int ii = _permute_ginv[i];
	    for (int k = _p_upper[i];  k < _ptRows[i + 1]; k++) {
	      const int jj = _permute_ginv[_indCols[k]];
	      xx(ii, m) -= _coef[_indVals[k]] * xx(jj, m);
	    }
	  }
	}
      } // if (isTrans)
    }
#else
    full_bw_multi<T>(isTrans, dim1, _num_null[n],
		     _diag_blocks[n].addrCoefs(), dim1,
		     nrhs,
		     xx.addrCoefs() + _p_front[n], _dim, &fop);
    int dim2 = _p_front[n] - _p_front[n - 1];
    if (isTrans) {
      ColumnMatrix<T> &lower = _lower_blocks[n];
      blas_gemm<T>(CblasNoTrans, CblasNoTrans, dim2, nrhs, dim1,
		   _none, lower.addrCoefs(), lower.nbRows(),
		   xx.addrCoefs() + _p_front[n], _dim,
		   _one,
		   xx.addrCoefs() + _p_front[n - 1], _dim);
    }
    else {
      ColumnMatrix<T> &upper = _upper_blocks[n];
      blas_gemm<T>(CblasNoTrans, CblasNoTrans, dim2, nrhs, dim1,
		   _none, upper.addrCoefs(), upper.nbRows(),
		   xx.addrCoefs() + _p_front[n], _dim,
		   _one,
		   xx.addrCoefs() + _p_front[n - 1], _dim);
    }
#endif
  } // loop : n
  {
    const int dim1 = _p_front[1] - _p_front[0];
#ifdef SPARSE_OFFDIAG
    full_fwbw_multi<T>(isTrans, dim1, _num_null[0],
		       _diag_blocks[0].addrCoefs(), dim1,
		       nrhs,
		       xx.addrCoefs() + _p_front[0], _dim);
#else
    full_bw_multi<T>(isTrans, dim1, _num_null[0],
		     _diag_blocks[0].addrCoefs(), dim1,
		     nrhs,
		     xx.addrCoefs() + _p_front[0], _dim, &fop);
#endif
  }
  if (nscol > 0) {
    // alpha = -1.0; beta = 1.0
    if (isTrans) {
      blas_gemm<T>(CblasNoTrans, CblasNoTrans, _dim, nrhs, nscol, _none,
		   _a21.addrCoefs(), _dim, y.addrCoefs(), nscol, _one,
		   xx.addrCoefs(), _dim);
    } 
    else {
      blas_gemm<T>(CblasNoTrans, CblasNoTrans, _dim, nrhs, nscol, _none,
		   _a12.addrCoefs(), _dim, y.addrCoefs(), nscol, _one,
		   xx.addrCoefs(), _dim);
    }
    for (int m = 0; m < nrhs; m++) {
      for (int i = 0; i < nscol; i++) {
	xx(_list_schur[i], m) = y(i, m);
      }
    }
    y.free();
  }
  if (flag_new2old) {
    for (int j = 0; j < nrhs; j++) {
      for (int i = 0; i < _dim; i++) {
	x(_new2old[i], j) = xx(_permute_ginv[i], j);
      }
    }
  }
  else {
    for (int j = 0; j < nrhs; j++) {
      for (int i = 0; i < _dim; i++) {
	x(i, j) = xx(_permute_ginv[i], j);
      }
    }
  }
}

template
void TridiagBlockMatrix<double>::SolveMulti(const bool flag_new2old,
					    const bool isTrans,
					    const int nhrs,
					    ColumnMatrix<double>& x,
					    const int nscol_);
template
void TridiagBlockMatrix<quadruple>::
SolveMulti(const bool flag_new2old,
	   const bool isTrans,
	   const int nhrs,
	   ColumnMatrix<quadruple>& x,
	   const int nscol_);
template
void TridiagBlockMatrix<complex<double>, double>::
SolveMulti(const bool flag_new2old, const bool isTrans,
	   const int nhrs, ColumnMatrix<complex<double> >& x,
	   const int nscol_);

template
void TridiagBlockMatrix<complex<quadruple>, quadruple>::
SolveMulti(const bool flag_new2old, const bool isTrans,
	   const int nhrs, ColumnMatrix<complex<quadruple> >& x,
	   const int nscol_);

//

template<typename T, typename U>
void TridiagBlockMatrix<T, U>::SolveSingle(const bool flag_new2old,
					   const bool isTrans, T* x_,
					   const int nscol_)
{
  const int nscol = (nscol_ == (-1) ? _nscol : nscol_);
  double fop;
  
  //  zz.resize(_maxdim);
  VectorArray<T> xx(_dim);
  VectorArray<T> zz(_maxdim);
  VectorArray<T> y;        // allocated only for nscol > 0
  if (flag_new2old) {
    for (int i = 0; i < _dim; i++) {
      xx[_permute_ginv[i]] = x_[_new2old[i]];
    }
  }
  else {
    for (int i = 0; i < _dim; i++) {
      xx[_permute_ginv[i]] = x_[i];
    }
  }
  // xx[ ] is reorderd to follow the permutation generated by
  // TridiagBlockFactorization
  if (nscol > 0) {
    y.init(nscol);
    for (int i = 0; i < nscol; i++) {
      y[i] = xx[_list_schur[i]];
      xx[_list_schur[i]] = _zero;
    }
    if (_isSymmetric) {
      blas_gemv<T>(CblasTrans, _dim, nscol, _none, _a12.addrCoefs(), _dim,
		   xx.addrCoefs(), 1, _one, y.addrCoefs(), 1);
    }
    else {
      if (isTrans) {
	blas_gemv<T>(CblasTrans, _dim, nscol, _none, _a12.addrCoefs(), _dim,
		     xx.addrCoefs(), 1, _one, y.addrCoefs(), 1);
      }
      else {
	blas_gemv<T>(CblasTrans, _dim, nscol, _none, _a21.addrCoefs(), _dim,
		     xx.addrCoefs(), 1, _one, y.addrCoefs(), 1);
      }
    }

    blas_trsv<T>(isTrans ? CblasUpper : CblasLower,
		 isTrans ? CblasTrans : CblasNoTrans,
		 CblasUnit, nscol,
		 _s22.addrCoefs(), nscol, y.addrCoefs(), 1);
    for (int i = 0; i < nscol; i++) {
      y[i] *= _s22(i, i);
    }
    blas_trsv<T>(isTrans ? CblasLower : CblasUpper,
		 isTrans ? CblasTrans : CblasNoTrans,
		 CblasUnit, nscol,
		 _s22.addrCoefs(), nscol, y.addrCoefs(), 1);
  } // if (nscol > 0)
  // forward substitution
  for (int n = 0; n < (_nfront - 1); n++) {
    // copy to working array zz
    const int dim1 = _p_front[n + 1] - _p_front[n];
    for (int i = 0; i < dim1; i++) {
      zz[i] = xx[_p_front[n] + i];
    }
#ifdef SPARSE_OFFDIAG
    full_fwbw_single<T>(isTrans, dim1, _num_null[n],
			_diag_blocks[n].addrCoefs(), dim1,
			&zz[0]);
    // sparse matrix * vector
    if (isTrans) {
      for (int i = _p_front[n]; i < _p_front[n + 1]; i++) {
	const int ii = _permute_ginv[i] - _p_front[n];
	for (int k = _p_upper[i];  k < _ptRows[i + 1]; k++) {
	  const int jj = _permute_ginv[_indCols[k]];
	  xx[jj] -= _coef[_indVals[k]] * zz[ii];
	}
      }
    }
    else {
      for (int i = _p_front[n + 1]; i < _p_front[n + 2]; i++) {
	const int ii = _permute_ginv[i];
	for (int k = _ptRows[i];  k < _p_diag[i]; k++) {
	  const int jj = _permute_ginv[_indCols[k]] - _p_front[n];
	  xx[ii] -= _coef[_indVals[k]] * zz[jj];
	}
      }
    }
#else
    full_fw_single<T>(isTrans, dim1, _num_null[n],
		      _diag_blocks[n].addrCoefs(), dim1,
		      &zz[0], &fop);
    int dim2 = _p_front[n + 2] - _p_front[n + 1];
    if (isTrans) {
      ColumnMatrix<T> &upper = _upper_blocks[n + 1];
      blas_gemv<T>(CblasTrans, dim1, dim2,
		   _none, upper.addrCoefs(), upper.nbRows(),
		   &zz[0], 1,
		   _one,
		   xx.addrCoefs() + _p_front[n + 1], 1);
    }
    else {
      ColumnMatrix<T> &lower = _lower_blocks[n + 1];
      blas_gemv<T>(CblasTrans, dim1, dim2,
		   _none, lower.addrCoefs(), lower.nbRows(),
		   &zz[0], 1,
		   _one,
		   xx.addrCoefs() + _p_front[n + 1], 1);
    }
    for (int i = 0; i < dim1; i++) {
      xx[_p_front[n] + i] = zz[i] * _diag_blocks[n](i, i);
    }
#endif
  }
#ifndef SPARSE_OFFDIAG
  {
    const int n = _nfront - 1;
    const int dim1 = _p_front[n + 1] - _p_front[n];
    for (int i = 0; i < dim1; i++) {
      zz[i] = xx[_p_front[n] + i];
    }
    full_fw_single<T>(isTrans, dim1, _num_null[n],
		      _diag_blocks[n].addrCoefs(), dim1,
		      &zz[0], &fop);
    for (int i = 0; i < dim1; i++) {
      xx[_p_front[n] + i] = zz[i] * _diag_blocks[n](i, i);
    }
  }
#endif
  // backward substitution
  for (int n = (_nfront - 1); n > 0; n--) {
    const int dim1 = _p_front[n + 1] - _p_front[n];
#ifdef SPARSE_OFFDIAG
    full_fwbw_single<T>(isTrans, dim1, _num_null[n],
			_diag_blocks[n].addrCoefs(), dim1,
			xx.addrCoefs() + _p_front[n]);
    if (_isSymmetric) {
      // Transposed sparse matrix * vector
      for (int i = _p_front[n]; i < _p_front[n + 1]; i++) {
	const int ii = _permute_ginv[i];
	for (int k = _ptRows[i];  k < _p_diag[i]; k++) {
	  const int jj = _permute_ginv[_indCols[k]];
	  xx[jj] -= _coef[_indVals[k]] * xx[ii];
	}
      }
    }
    else {
      if (isTrans) {
	for (int i = _p_front[n]; i < _p_front[n + 1]; i++) {
	  const int ii = _permute_ginv[i];
	  for (int k = _ptRows[i];  k < _p_diag[i]; k++) {
	    const int jj = _permute_ginv[_indCols[k]];
	    xx[jj] -= _coef[_indVals[k]] * xx[ii];
	  }
	}
      }
      else {
	for (int i = _p_front[n - 1]; i < _p_front[n]; i++) {
	  const int ii = _permute_ginv[i];
	  for (int k = _p_upper[i];  k < _ptRows[i + 1]; k++) {
	    const int jj = _permute_ginv[_indCols[k]];
	    xx[ii] -= _coef[_indVals[k]] * xx[jj];
	  }
	}
      }
    }
#else
    full_bw_single<T>(isTrans, dim1, _num_null[n],
		      _diag_blocks[n].addrCoefs(), dim1,
		      xx.addrCoefs() + _p_front[n], &fop);
    int dim2 = _p_front[n] - _p_front[n - 1];
    if (isTrans) {
      ColumnMatrix<T> &lower = _lower_blocks[n];
      blas_gemv<T>(CblasNoTrans, dim2, dim1,
		   _none, lower.addrCoefs(), lower.nbRows(),
		   xx.addrCoefs() + _p_front[n], 1,
		   _one,
		   xx.addrCoefs() + _p_front[n - 1], 1);
    }
    else {
      ColumnMatrix<T> &upper = _upper_blocks[n];
      blas_gemv<T>(CblasNoTrans, dim2, dim1,
		   _none, upper.addrCoefs(), upper.nbRows(),
		   xx.addrCoefs() + _p_front[n], 1,
		   _one,
		   xx.addrCoefs() + _p_front[n - 1], 1);
    }
#endif
  } // loop : n
  {
    const int dim1 = _p_front[1] - _p_front[0];
#ifdef SPARSE_OFFDIAG
    full_fwbw_single<T>(isTrans, dim1, _num_null[0],
			_diag_blocks[0].addrCoefs(), dim1,
			xx.addrCoefs() + _p_front[0]);
#else
    full_bw_single<T>(isTrans, dim1, _num_null[0],
		      _diag_blocks[0].addrCoefs(), dim1,
		      xx.addrCoefs() + _p_front[0], &fop);
#endif    
  }
  if (nscol > 0) {
    // alpha = -1.0; beta = 1.0
    if (isTrans) {
      blas_gemv<T>(CblasNoTrans, _dim, nscol, _none, _a21.addrCoefs(), _dim,
		   y.addrCoefs(), 1, _one, xx.addrCoefs(), 1);
    }
    else {
      blas_gemv<T>(CblasNoTrans, _dim, nscol, _none, _a12.addrCoefs(), _dim,
		   y.addrCoefs(), 1, _one, xx.addrCoefs(), 1);
    }
    for (int i = 0; i < nscol; i++) {
      xx[_list_schur[i]] = y[i];
    }
  }
  if (flag_new2old) {
    for (int i = 0; i < _dim; i++) {
      x_[_new2old[i]] = xx[_permute_ginv[i]];
    }
  }
  else {
    for (int i = 0; i < _dim; i++) {
      x_[i] = xx[_permute_ginv[i]];
    }
  }
  //  delete [] xx;
}

template
void TridiagBlockMatrix<double>::SolveSingle(const bool flag_new2old,
					     const bool isTrans,
					     double *x, const int nscol_);

template
void TridiagBlockMatrix<quadruple>::
SolveSingle(const bool flag_new2old,
	    const bool isTrans,
	    quadruple *x, const int nscol_);

template
void TridiagBlockMatrix<complex<double>, double>::
SolveSingle(const bool flag_new2old,
	    const bool isTrans,
	    complex<double> *x, const int nscol_);

template
void TridiagBlockMatrix<complex<quadruple>, quadruple>::
SolveSingle(const bool flag_new2old,
	    const bool isTrans,
	    complex<quadruple> *x, const int nscol_);
//

template<typename T, typename U>
void TridiagBlockMatrix<T, U>::ForwardUpper(bool isTransposed,
					    int ncol, ColumnMatrix<T> &b,
					    vector<int>& i0,
					    vector<T> &dscale)
{
  int ncol0;
  double fop;
#ifdef STORE_WHOLE_FACTORIZED
  blas_trsm<T>(CblasLeft,
	       (isTransposed ? CblasUpper : CblasLower),
	       (isTransposed ? CblasTrans : CblasNoTrans),
	       CblasUnit,
	       _dim, ncol, _one, _factorized_whole.addrCoefs(), _dim,
	       b.addrCoefs(), b.nbRows());
  for (int i = 0; i < _dim; i++) {
    dscale[i] = _factorized_whole(i, i);
  }
#else // #ifdef STORE_WHOLE_FACTORIZED
  ColumnMatrix<T> zz;
  vector<int> i0_local;
  i0_local.resize(ncol);

  zz.init(_maxdim, ncol);
  //  zz.ZeroClear();
  for (int n = 0; n < (_nfront - 1); n++) {
    const int dim1 = _p_front[n + 1] - _p_front[n];
    ColumnMatrix<T> &diag = _diag_blocks[n];
    ncol0 = ncol;
    for (int j = 0; j < ncol; j++) {
      int itmp = std::max((i0[j] - _p_front[n]), 0);
      if (i0[j] >= _p_front[n + 1]) {
	ncol0 = j;
	break;
      }
      i0_local[j] = itmp;
    }
    if (ncol0 > 0) {
      full_fw_multiprofile<T>(isTransposed,
			      dim1, _num_null[n], ncol0, diag.addrCoefs(), dim1,
			      b.addrCoefs() + _p_front[n], _dim,
			      i0_local, &fop);
    }
#ifdef SPARSE_OFFDIAG
    // dddiag_times_s1
    for (int j = 0; j < ncol0; j++) {
      for (int i = 0; i < i0_local[j]; i++) {
	zz(i, j) = _zero;
      }
      for (int i = i0_local[j]; i < dim1; i++) {
	zz(i, j) = b(i + _p_front[n], j) * diag(i, i);
      }
      fop += (double)(dim1 - i0[j]);
    }
    if (ncol0 > 0) {
      full_bw_multi<T>(isTransposed,
		       dim1, _num_null[n], diag.addrCoefs(), dim1, ncol0,
		       zz.addrCoefs(), _maxdim,
		       &fop);
    }
    // sparse matrix * dense matrix
    if (isTransposed) {
       for (int j = 0; j < ncol0; j++) {
	 for (int i = _p_front[n]; i < _p_front[n + 1]; i++) {
	   const int i1 = _permute_ginv[i] - _p_front[n];
	   for (int k = _p_upper[i];  k < _ptRows[i + 1]; k++) {
	     const int ii = _permute_ginv[_indCols[k]];
	     b(ii, j) -= _coef[_indVals[k]] * zz(i1, j);
	  }
	}
      }
    }
    else {
      for (int j = 0; j < ncol0; j++) {
	for (int i = _p_front[n + 1]; i < _p_front[n + 2]; i++) {
	  const int ii = _permute_ginv[i];
	  for (int k = _ptRows[i];  k < _p_diag[i]; k++) {
	    const int i1 = _permute_ginv[_indCols[k]] - _p_front[n];
	    b(ii, j) -= _coef[_indVals[k]] * zz(i1, j);
	  }
	}
      }
    }
#else // #ifdef SPARSE_OFFDIAG
    int dim2 = _p_front[n + 2] - _p_front[n + 1];
    if (isTransposed) {
      ColumnMatrix<T> &upper = _upper_blocks[n + 1];
	blas_gemm<T>(CblasTrans, CblasNoTrans, dim2, ncol0, dim1,
		     _none, upper.addrCoefs(), upper.nbRows(),
		     b.addrCoefs() + _p_front[n], b.nbRows(),
		     _one,
		     b.addrCoefs() + _p_front[n + 1], b.nbRows());
    }
    else {
      ColumnMatrix<T> &lower = _lower_blocks[n + 1];
      blas_gemm<T>(CblasTrans, CblasNoTrans, dim2, ncol0, dim1,
		   _none, lower.addrCoefs(), lower.nbRows(),
		   b.addrCoefs() + _p_front[n], b.nbRows(),
		   _one,
		   b.addrCoefs() + _p_front[n + 1], b.nbRows());
    }
#endif // #ifdef SPARSE_OFFDIAG
  } // loop : n
  {
    int n = (_nfront - 1);
    const int dim1 = _p_front[n + 1] - _p_front[n];
    ColumnMatrix<T> &diag = _diag_blocks[n];
    ncol0 = ncol;
    for (int j = 0; j < ncol; j++) {
      int itmp = std::max((i0[j] - _p_front[n]), 0);
      if (i0[j] >= _p_front[n + 1]) {
	ncol0 = j;
	break;
      }
      i0_local[j] = itmp;
    }
    if (ncol0 > 0) {
      full_fw_multiprofile<T>(isTransposed,
			      dim1, _num_null[n], ncol0, diag.addrCoefs(), dim1,
			      b.addrCoefs() + _p_front[n], _dim,
			      i0_local, &fop);
    }
  }
  for (int n = 0; n < _nfront; n++) {
    const int dim1 = _p_front[n + 1] - _p_front[n];
     ColumnMatrix<T> &diag = _diag_blocks[n];
     for (int i = 0; i < dim1; i++) {
       dscale[_p_front[n] + i] = diag(i, i);
    }
  }
  zz.free();
#endif // #ifdef STORE_WHOLE_FACTORIZED
}

template
void TridiagBlockMatrix<double>::ForwardUpper(bool isTransposed,
					      int ncol, ColumnMatrix<double> &b,
					      vector<int>& i0,
					      vector<double> &dscale);
template
void TridiagBlockMatrix<quadruple>::ForwardUpper(bool isTransposed,
						 int ncol,
						 ColumnMatrix<quadruple> &b,
						 vector<int>& i0,
						 vector<quadruple> &dscale);
template
void TridiagBlockMatrix<complex<double>, double>::
ForwardUpper(bool isTransposed,
	     int ncol, ColumnMatrix<complex<double> > &b,
	     vector<int>& i0,
	     vector<complex<double> > &dscale);

template
void TridiagBlockMatrix<complex<quadruple>, quadruple>::
ForwardUpper(bool isTransposed,
	     int ncol, ColumnMatrix<complex<quadruple> > &b,
	     vector<int>& i0,
	     vector<complex<quadruple> > &dscale);
//

void RenumberCSR(const bool shrink_flag,
		 const int dim,
		 const int nnz,
		 vector<int> &b2a,
		 vector<int> &a2b,
		 const int *aptrows,
		 const int *aindcols,
		 const int *aindvals,
		 vector<int> &bptrows,
		 vector<int> &bindcols, vector<int> &bindvals)
{
  vector<int> *bbindcols, *bbindvals;
  if (shrink_flag) {
    bbindcols = new vector<int>;
    bbindvals = new vector<int>;
    (*bbindcols).resize(nnz);
    (*bbindvals).resize(nnz);
  }
  else {
    bbindcols = &bindcols;
    bbindvals = &bindvals;
  }
  bptrows[0] = 0;
  for (int i = 0; i < dim; i++) {
    const int ii = b2a[i];
    bptrows[i + 1] = bptrows[i] + (aptrows[ii + 1] - aptrows[ii]);
  }
  for (int i = 0; i < dim; i++) {
    const int ii = b2a[i];
    list<int> local_indcols;
    list<int> local_indvals;
    for (int kk = aptrows[ii]; kk < aptrows[ii + 1]; kk++) {
      const int jj = aindcols[kk];
      const int j = a2b[jj];
      const int ll = aindvals[kk];
      list<int>::iterator ic = local_indcols.begin();
      list<int>::iterator iv = local_indvals.begin();
      bool flag = false;
      for ( ; ic != local_indcols.end(); ++ic, ++iv) {
	if (*ic > j) {
	  local_indcols.insert(ic, j);
	  local_indvals.insert(iv, ll);
	  flag = true;
	  break;
	}
      }
      if (flag == false) {
	local_indcols.push_back(j);
	local_indvals.push_back(ll);
      }
    } // loop : kk
    {
      list<int>::const_iterator ic = local_indcols.begin();
      list<int>::const_iterator iv = local_indvals.begin();
      int k = bptrows[i];
      for ( ; k < bptrows[i + 1]; k++, ++ic, ++iv) {
	(*bbindcols)[k] = *ic;
	(*bbindvals)[k] = *iv;
      }
    }
    local_indcols.clear();
    local_indvals.clear();
  } // loop : i
  if (shrink_flag) {
    int nnz1 = bptrows[dim];
    if (nnz1 != bindcols.size()) {
      fprintf(stderr, "%s %d : %d != %d\n",
	      __FILE__, __LINE__, nnz1, (int)bindcols.size());
    }
    for (int i = 0; i < nnz1; i++) {
      bindcols[i] = (*bbindcols)[i];
      bindvals[i] = (*bbindvals)[i];
    }
    (*bbindcols).clear();
    (*bbindvals).clear();
    delete bbindcols;
    delete bbindvals;
  }
}

void RenumberCSR(const int dim,
		 vector<int> &b2a,
		 vector<int> &a2b,
		 const int *aptrows,
		 const int *aindcols,
		 vector<int> &bptrows,
		 vector<int> &bindcols)
{
  //  const int dim = b2a.size();
  bptrows[0] = 0;
  for (int i = 0; i < dim; i++) {
    const int ii = b2a[i];
    bptrows[i + 1] = bptrows[i] + (aptrows[ii + 1] - aptrows[ii]);
  }
  for (int i = 0; i < dim; i++) {
    const int ii = b2a[i];
    list<int> local_indcols;
    for (int kk = aptrows[ii]; kk < aptrows[ii + 1]; kk++) {
      const int jj = aindcols[kk];
      const int j = a2b[jj];
      list<int>::iterator ic = local_indcols.begin();
      bool flag = false;
      for ( ; ic != local_indcols.end(); ++ic) {
	if (*ic > j) {
	  local_indcols.insert(ic, j);
	  flag = true;
	  break;
	}
      }
      if (flag == false) {
	local_indcols.push_back(j);
      }
    } // loop : kk
    {
      list<int>::const_iterator ic = local_indcols.begin();
      int k = bptrows[i];
      for ( ; k < bptrows[i + 1]; k++, ++ic) {
	bindcols[k] = *ic;
      }
    }
    local_indcols.clear();
  } // loop : i
}

void TridiagStruct(vector<int> &ptrow, vector<int> &indcols,
		   const int nfront,
		   vector<int> &pfront, vector<int> &p_diag,
		   vector<int> &p_upp)
//, vector<int> &p_dia_blo)
{
  for (int n = 0; n < nfront; n++) {
    const int i1 = pfront[n];
    const int i2 = pfront[n + 1];
    for (int i = i1; i < i2; i ++) {
      for (int k = ptrow[i]; k < ptrow[i + 1]; k++) {
	if (indcols[k] >= i1) {
	  p_diag[i] = k;
	  break;
	}
      }
    } // loop : i
  }
  for (int n = 0; n < nfront; n++) {
    const int i1 = pfront[n];
    const int i2 = pfront[n + 1];
    for (int i = i1; i < i2; i ++) {
      bool flag = false;
      for (int k = ptrow[i]; k < ptrow[i + 1]; k++) {
	if (indcols[k] >= i2) {
	  p_upp[i] = k;
	  flag = true;
	  break;
	}
      }
      if (flag == false) {
	p_upp[i] = ptrow[i + 1];
      }
    } // loop : i
  }
}

void GenPermuteOffdiag(const int pfront0, const int pfront1, const int pfront2,
		       vector<int> &ptrow, vector<int> &ptdiag,
		       vector<int> &indcols,
		       const int *permute_diag, int *permute_diag_inv,
		       int *permute_offdiag, int *permute_offdiag_inv,
		       vector<int> &i0)
{
  const int dim1 = pfront1 - pfront0;
  const int dim2 = pfront2 - pfront1;
  vector<int> i0_loc;
  i0_loc.resize(dim2);
  i0.resize(dim2);
  for (int i = 0; i < dim1; i++) {
    permute_diag_inv[permute_diag[i]] = i;
  }
  list<int> perm_offdg;
  for (int j = 0; j < dim2; j++) {
    const int jj = j + pfront1;
    int itmp = dim1;
    for (int k = ptrow[jj]; k < ptdiag[jj]; k++) {
      const int i = permute_diag_inv[indcols[k] - pfront0];
      itmp = itmp > i ? i : itmp;
    }
    bool flag = false;
    i0_loc[j] = itmp;
    for (list<int>::iterator it = perm_offdg.begin(); it != perm_offdg.end();
	 ++it) {
      if (itmp < i0_loc[(*it)]) {
	perm_offdg.insert(it, j);
	flag = true;
	break;
      }
    }
    if (flag == false) {
      perm_offdg.push_back(j);
    }
  } // loop : j
  { // copy list<int> to vector<int>
    list<int>::const_iterator jt = perm_offdg.begin();
    int j = 0;
    for (; j < dim2; j++, ++jt) {
      permute_offdiag[j] = (*jt);
      i0[j] = i0_loc[(*jt)];
    }
    for (int i = 0 ; i < dim2; i++) {
      permute_offdiag_inv[permute_offdiag[i]] = i;
    }
  }
  perm_offdg.clear();
  i0_loc.clear();
}

void GenPermuteUpper(const int nrow,
		     vector<int> &remap_eqn,
		     const int ncol,
		     const int *ptrow, const int *indcols,
		     const int nfront,
		     vector<int> &p_fronts,
		     vector<int> &new2old,
		     vector<int> &permute_diag,
		     vector<int> &permute_diag_inv,
		     vector<int> &permute_upper,
		     vector<int> &permute_upper_inv,
		     vector<int> &i0,
		     const bool verbose,
		     FILE *fp)
{
  vector<bool> colflag;
  vector<int> permute_tmp;
  colflag.resize(ncol, false);
  permute_diag_inv.resize(nrow);
  permute_tmp.resize(nrow);
  permute_upper.resize(ncol);
  permute_upper_inv.resize(ncol);
  i0.resize(ncol);
  for (int n = 0; n < nfront; n++) {
    for (int i = p_fronts[n]; i < p_fronts[n + 1]; i++) {
      permute_tmp[i] = new2old[permute_diag[i] + p_fronts[n]];
    }
  }
#ifdef PERMUTE_UPPER
  int l = 0;
  for (int i = 0; i < nrow; i++) {
    const int ii = permute_tmp[i];
    for (int k = ptrow[ii]; k < ptrow[ii + 1]; k++) {
      const int jj = indcols[k];
      if (colflag[jj] == false) {
	colflag[jj] = true;
	i0[l] = i;
	permute_upper[l++] = jj;
      }
    }
  }
  for(int j = 0; j < ncol; j++) {
    if (colflag[j] == false) {
      i0[l] = nrow;
      permute_upper[l++] = j;
    }
  }
  for (int j = 0; j < ncol; j++) {
    permute_upper_inv[permute_upper[j]] = j;
  }
#else
  for(int j = 0; j < ncol; j++) {
    permute_upper[j] = j;
    permute_upper_inv[j] = j;
  }
#endif
  for (int i = 0; i < nrow; i++) {
    permute_diag_inv[remap_eqn[permute_tmp[i]]] = i;
  }
  colflag.clear();
  permute_tmp.clear();
}

template<typename T, typename U>
void TridiagBlockMatrix<T, U>::extract_column(const int jcol, T *dcol)
{
  for (int i = 0; i < _dim; i++) {
    for (int k = _ptRows[i]; k < _ptRows[i + 1]; k++) {
      if (_indCols[k] == jcol) {
	dcol[i] = _coef[_indVals[k]];
	break;
      }
    }
  }
}

template
void TridiagBlockMatrix<double>::extract_column(const int jcol,
							double *dcol);

template
void TridiagBlockMatrix<quadruple>::extract_column(const int jcol,
							      quadruple *dcol);

template
void TridiagBlockMatrix<complex<double>, double>::
extract_column(const int jcol, 
	       complex<double> *dcol);

template
void TridiagBlockMatrix<complex<quadruple>, quadruple>::
extract_column(const int jcol, 
	       complex<quadruple> *dcol);
//

template<typename T, typename U>
void TridiagBlockMatrix<T, U>::extract_row(const int irow, T *drow)
{
  for (int k = _ptRows[irow]; k < _ptRows[irow + 1]; k++) {
    const int jj = _indCols[k];
    drow[jj] = _coef[_indVals[k]];
  }
}

template
void TridiagBlockMatrix<double>::extract_row(const int irow,
					     double *drow);

template
void TridiagBlockMatrix<quadruple>::extract_row(const int irow,
						quadruple *drow);

template
void TridiagBlockMatrix<complex<double>, double>::extract_row(const int irow,
						       complex<double> *drow);
			
template
void TridiagBlockMatrix<complex<quadruple>, quadruple>::
extract_row(const int irow,
	    complex<quadruple> *drow);
//

template<typename T>
void FillBlockSparse(const T *coef,
		     const int dim,
		     vector<int>& map_eqn,
		     const int *ptrow,
		     const int *indcols,
		     const int *indvals,
		     vector<int> &old2new_i,
		     vector<int> &old2new_j,
		     ColumnMatrix<T> &b)
{
  for (int i = 0; i < dim; i++) {
    const int ii = old2new_i[i];
    const int i0 = map_eqn[i];
    for (int k = ptrow[i0]; k < ptrow[i0 + 1]; k++) {
      const int jj = old2new_j[indcols[k]];
      b(ii, jj) = coef[indvals[k]];
    }
  } // loop : i
}

template
void FillBlockSparse<double>(const double *coef,
			     const int dim,
			     vector<int>& map_eqn,
			     const int *ptrow,
			     const int *indcols,
			     const int *indvals,
			     vector<int> &old2new_i,
			     vector<int> &old2new_j,
			     ColumnMatrix<double> &b);

template
void FillBlockSparse<quadruple>(const quadruple *coef,
				const int dim,
				vector<int>& map_eqn,
				const int *ptrow,
				const int *indcols,
				const int *indvals,
				vector<int> &old2new_i,
				vector<int> &old2new_j,
				ColumnMatrix<quadruple> &b);

template
void FillBlockSparse<complex<double> >(const complex<double> *coef,
				       const int dim,
				       vector<int>& map_eqn,
				       const int *ptrow,
				       const int *indcols,
				       const int *indvals,
				       vector<int> &old2new_i,
				       vector<int> &old2new_j,
				       ColumnMatrix<complex<double> > &b);


template
void FillBlockSparse<complex<quadruple> >(const complex<quadruple> *coef,
					  const int dim,
					  vector<int>& map_eqn,
					  const int *ptrow,
					  const int *indcols,
					  const int *indvals,
					  vector<int> &old2new_i,
					  vector<int> &old2new_j,
					  ColumnMatrix<complex<quadruple> > &b);
//

bool isnegative(const double x) {
  return (x < 0.0);
}

bool isnegative(const quadruple x) {
  return (x < quadruple(0.0));
}

bool isnegative(const complex<double> &x) {
  return ((x.real() < 0.0) && (x.imag() == 0.0));
}

bool isnegative(const complex<quadruple> &x) {
  return (x.real() < quadruple(0.0)) && (x.imag() == quadruple(0.0));
}

template<typename T, typename U>
int TridiagBlockMatrix<T, U>::NumNegativeDiags()
{
  int count = 0 ;
  for (int n = 0; n < _nfront; n++) {
    ColumnMatrix<T> &diag = _diag_blocks[n];
    for (int i = _p_front[n]; i < _p_front[n + 1]; i++) {
      const int ii = i - _p_front[n];
      if (isnegative(diag(ii, ii))) {
	count++;
      }
    }
  }
  for (int i = 0; i < _nscol; i++) {
    if (isnegative(_s22(i, i))) {
      count++;
    }
  }	
  return count;
}

template
int TridiagBlockMatrix<double>::NumNegativeDiags();

template
int TridiagBlockMatrix<quadruple>::NumNegativeDiags();

template
int TridiagBlockMatrix<complex<double>, double>::NumNegativeDiags();

template
int TridiagBlockMatrix<complex<quadruple>, quadruple>::NumNegativeDiags();
//

template<typename T, typename U>
void TridiagBlockMatrix<T, U>::KernelBasis(const bool isTrans,
					   ColumnMatrix<T> &a12)
{
  //  const T zero(0.0);
  //  const T none(-1.0);
  
  vector<T> v;
  v.resize(_dim);
  a12.free();
  a12.init(_dim, _n0);
  a12.ZeroClear();
  for (int j = 0; j < _n0; j++) {
    if (isTrans) {
      extract_row(_list_elim[j], a12.addrCoefs() + (j * _dim));
    }
    else {
      extract_column(_list_elim[j], a12.addrCoefs() + (j * _dim));
    }
    for (int i = 0; i < _n0; i++) {
      a12(_list_elim[i], j) = _zero;
    }
  }
  if (_n0 > 1) {
    SolveMulti(false, isTrans, _n0, a12);
  }
  else {
    SolveSingle(false, isTrans, a12.addrCoefs());  
  }
  for (int j = 0; j < _n0; j++) {
    for (int i = 0; i < _dim; i++) {
      v[i] = a12(i, j);
    }
    for (int i = 0; i < _n0; i++) {
      v[_list_elim[i]] = _zero;
    }
    v[_list_elim[j]] = _none;
    for (int i = 0; i < _dim; i++) {
      a12(_new2old[i], j) = v[i];
    }
  }

  v.clear();
}

template
void TridiagBlockMatrix<double>::KernelBasis(const bool isTrans,
					     ColumnMatrix<double> &a12);

template
void TridiagBlockMatrix<quadruple>::KernelBasis(const bool isTrans,
						ColumnMatrix<quadruple> &a12);

template
void TridiagBlockMatrix<complex<double>, double>::
KernelBasis(const bool isTrans,
	    ColumnMatrix<complex<double> > & a12);

template
void TridiagBlockMatrix<complex<quadruple>, quadruple >::
KernelBasis(const bool isTrans,
	    ColumnMatrix<complex<quadruple> > & a12);
//
