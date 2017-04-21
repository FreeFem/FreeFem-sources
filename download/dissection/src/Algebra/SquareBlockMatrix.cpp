/*! \file   SquareBlockMatrix.cpp
    \brief  Block storage for symmetric/unsymmetric Square matrix
    \author Atsushi Suzuki, Laboratoire Jacques-Louis Lions
    \date   Jun. 11th 2013
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

#include "Algebra/SquareBlockMatrix.hpp"
#include "Compiler/blas.hpp"                   // #include<complex> is inside

template<typename T>
void SquareBlockMatrix<T>::init(int dim, int block_size, bool isSym,
				int first_block)
{
  int itmp;
  _pivrelaxed = false;
  _coefs_alloc_status = false; // for safty
  _dim = dim;
  _block_size = block_size;
  _isSym = isSym;
  _block_size2 = _block_size * _block_size;

  if (dim == 0) {
    _isblocked = true;
    _kernelDetected = false;
    _lower_allocated = false;
    _coefs = NULL;
    _coefs_alloc_status = false;
    _allocation_status = NULL;
    return;
  }
  if (first_block > 0) {
    _isdecomposed = true;
    _dim0 = first_block;
    _dim1 = _dim - _dim0;
      
    _num_blocks0 = (_dim0 + _block_size - 1) / _block_size;
    _num_blocks1 = (_dim1 + _block_size - 1) / _block_size;
    _num_blocks = _num_blocks0 + _num_blocks1;
      
    _block_sizes.resize(_num_blocks);
    itmp = _dim0 % _block_size;
    _block_size_last0 = (itmp == 0) ? _block_size : itmp;
    itmp = _dim1 % _block_size;
    _block_size_last1 = (itmp == 0) ? _block_size : itmp;
    for (int i = 0; i < (_num_blocks0 - 1); i++) {
      _block_sizes[i] = _block_size;
    }
    _block_sizes[_num_blocks0 - 1] = _block_size_last0;
    if (_num_blocks1 > 0) {
      for (int i = _num_blocks0; i < (_num_blocks - 1); i++) {
	_block_sizes[i] = _block_size;
      }
      _block_sizes[_num_blocks - 1] = _block_size_last1;
    }
    else {
      _block_size_last1 = 0;
    }
  }
  else {
    _isdecomposed = false;
    _num_blocks = (_dim + _block_size - 1) / _block_size;
    _block_sizes.resize(_num_blocks);
    itmp = _dim % _block_size;
    _block_size_last = (itmp == 0) ? _block_size : itmp;
    for (int i = 0; i < (_num_blocks - 1); i++) {
      _block_sizes[i] = _block_size;
    }
    _block_sizes[_num_blocks - 1] = _block_size_last;
  }
  _nsing_block.resize(_num_blocks);
  if (_num_blocks > 0) {
    for (int i = 0; i < _num_blocks; i++) {
      _nsing_block[i] = 0;
    }
  }
  _isblocked = true;
  _kernelDetected = false;
  _lower_allocated = false;

  if (_num_blocks > 0) {
    try {
      _coefs = new T*[_num_blocks * _num_blocks];  // need to be managed
    } catch (const std::bad_alloc& e) {
      fprintf(stderr, "%s %d : allocation failed : %s", __FILE__, __LINE__,
	      e.what());
    }      
    try {
      _allocation_status = new bool[_num_blocks * _num_blocks];
    } catch (const std::bad_alloc& e) {
      fprintf(stderr, "%s %d : allocation failed : %s", __FILE__, __LINE__,
	      e.what());
    }      
    _coefs_alloc_status = true;
    for (int i = 0; i < (_num_blocks * _num_blocks); i++) {
      _allocation_status[i] = false;
    }
    if (_isSym) {
      _size = 0;
      for (int i = 0; i < _num_blocks; i++) {
	for (int j = i; j < _num_blocks; j++) {
	  _size += _block_sizes[i] * _block_sizes[j];
	}
      }
    }
    else {
      _size = 0;
      for (int i = 0; i < _num_blocks; i++) {
	for (int j = 0; j < _num_blocks; j++) {
	  _size += _block_sizes[i] * _block_sizes[j];
	}
      }
    }
  }
  else {
    _coefs = NULL;
    _coefs_alloc_status = false;
    _allocation_status = NULL;
  }
}

template
void SquareBlockMatrix<double>::init(int dim, int block_size, bool isSym,
				     int first_block);
template
void SquareBlockMatrix<quadruple>::init(int dim, int block_size, bool isSym,
					int first_block);
template
void SquareBlockMatrix<complex<double> >::init(int dim, int block_size,
					       bool isSym, int first_block);
template
void SquareBlockMatrix<complex<quadruple> >::init(int dim, int block_size,
						  bool isSym, int first_block);
//

template<typename T>
void SquareBlockMatrix<T>::allocateBlock(int i, int j)
{
  const int nrow = _block_sizes[i];
  const int ncol = _block_sizes[j];
  //    fprintf(stderr, "%s %d : allocateBlock %d %d : %d x %d\n", 
  //	    __FILE__, __LINE__, i, j, nrow, ncol);
  const int itmp = i + j * _num_blocks;
  if (!_allocation_status[itmp]) {
    try {
      _coefs[itmp] = new T[nrow * ncol];
    } catch (const std::bad_alloc& e) {
      fprintf(stderr, "%s %d : allocation failed : %s", __FILE__, __LINE__,
	      e.what());
    }      
    _allocation_status[itmp] = true;
  }
}

template
void SquareBlockMatrix<double>::allocateBlock(int i, int j);
template
void SquareBlockMatrix<quadruple>::allocateBlock(int i, int j);
template
void SquareBlockMatrix<complex<double> >::allocateBlock(int i, int j);
template
void SquareBlockMatrix<complex<quadruple> >::allocateBlock(int i, int j);
//

template<typename T>
void SquareBlockMatrix<T>::allocate()
{
  if (_isSym) {
    //    int ioffset = 0;
    for (int i = 0; i < _num_blocks; i++) {
      for (int j = i; j < _num_blocks; j++) {
	const int nrow = _block_sizes[i];
	const int ncol = _block_sizes[j];
	const int itmp = i + j * _num_blocks;
	if (!_allocation_status[itmp]) {
	  try {
	    _coefs[itmp] = new T[nrow * ncol];
	  } catch (const std::bad_alloc& e) {
	    fprintf(stderr, "%s %d : allocation failed : %s",
		    __FILE__, __LINE__,
		    e.what());
	  }      
	}
	_allocation_status[itmp] = true;
      }
    }
  }
  else {
    //    int ioffset = 0;
    for (int i = 0; i < _num_blocks; i++) {
      for (int j = 0; j < _num_blocks; j++) {
	const int nrow = _block_sizes[i];
	const int ncol = _block_sizes[j];
	const int itmp = i + j * _num_blocks;
	if (!_allocation_status[itmp]) {
	  try {
	    _coefs[itmp] = new T[nrow * ncol];
	  } catch (const std::bad_alloc& e) {
	    fprintf(stderr, "%s %d : allocation failed : %s",
		    __FILE__, __LINE__,
		    e.what());
	  }      
	}
	_allocation_status[itmp] = true;
      }
    }
  }
}

template
void SquareBlockMatrix<double>::allocate();
template
void SquareBlockMatrix<quadruple>::allocate();
template
void SquareBlockMatrix<complex<double> >::allocate();
template
void SquareBlockMatrix<complex<quadruple> >::allocate();
//

template<typename T>
void SquareBlockMatrix<T>::free(int i, int j)
{
  const int itmp = i + j * _num_blocks;
  if (_allocation_status[itmp]) {
    delete [] _coefs[itmp];  // need to be mananged
  }
  _allocation_status[itmp] = false;
}

template
void SquareBlockMatrix<double>::free(int i, int j);
template
void SquareBlockMatrix<quadruple>::free(int i, int j);
template
void SquareBlockMatrix<complex<double> >::free(int i, int j);
template
void SquareBlockMatrix<complex<quadruple> >::free(int i, int j);
//

template<typename T>
void SquareBlockMatrix<T>::free()
{
  if (_isSym) {
    for (int i = 0; i < _num_blocks; i++) {
      for (int j = i; j < _num_blocks; j++) {
	const int itmp = i + j * _num_blocks;
	if (_allocation_status[itmp]) {
	  delete [] _coefs[itmp];  // need to be mananged
	}
	_allocation_status[itmp]= false;
      }
    }
  }
  else {
    for (int i = 0; i < _num_blocks; i++) {
      for (int j = 0; j < _num_blocks; j++) {
	const int itmp = i + j * _num_blocks;
	if (_allocation_status[itmp]) {
	  delete [] _coefs[itmp];  // need to be managed
	}
	_allocation_status[itmp] = false;
      }
    }
  }
}
template
void SquareBlockMatrix<double>::free();
template
void SquareBlockMatrix<quadruple>::free();
template
void SquareBlockMatrix<complex<double> >::free();
template
void SquareBlockMatrix<complex<quadruple> >::free();
//

template<typename T>
T* SquareBlockMatrix<T>::addrCoefBlock(int i, int j) 
{
  if (_isdecomposed) {                        // is necessary?
    if (_isSym && (i > j)) {
      return _coefs[j + i * _num_blocks];
    }
  }
  return _coefs[i + j * _num_blocks];
}

template
double* SquareBlockMatrix<double>::addrCoefBlock(int i, int j);

template
quadruple* SquareBlockMatrix<quadruple>::addrCoefBlock(int i, int j);

template
complex<double>* SquareBlockMatrix<complex<double> >::
addrCoefBlock(int i, int j);

template
complex<quadruple>* SquareBlockMatrix<complex<quadruple> >::
addrCoefBlock(int i, int j);
//

template<typename T>
void SquareBlockMatrix<T>::ZeroClear()
{
  const T zero(0.0);
  if (_isSym) {
    for (int i = 0; i < _num_blocks; i++) {
      for (int j = i; j < _num_blocks; j++) {
	const int nrow = _block_sizes[i];
	const int ncol = _block_sizes[j];
	const int nsize = nrow * ncol;
	const int itmp = i + j * _num_blocks;
	for (int k = 0; k < nsize; k++) {
	  _coefs[itmp][k] = zero;   // needs to be replaced by memcpy
	}
      }
    }
  }
  else {
    for (int i = 0; i < _num_blocks; i++) {
      for (int j = 0; j < _num_blocks; j++) {
	const int nrow = _block_sizes[i];
	const int ncol = _block_sizes[j];
	const int nsize = nrow * ncol;
	const int itmp = i + j * _num_blocks;
	for (int k = 0; k < nsize; k++) {
	  _coefs[itmp][k] = zero;   // needs to be replaced by memcpy
	}
      }
    }
  }
}

template
void SquareBlockMatrix<double>::ZeroClear();

template
void SquareBlockMatrix<quadruple>::ZeroClear();

template
void SquareBlockMatrix<complex<double> >::ZeroClear();

template
void SquareBlockMatrix<complex<quadruple> >::ZeroClear();
//

template<typename T>
T& SquareBlockMatrix<T>::diag(int i)
{
  int i0, i1;
  if(_isdecomposed) {
    fprintf(stderr, "%s %d : operator() is not suppoesd to be used\n",
	    __FILE__, __LINE__);
    i0 = BlockIndex(i);
    i1 = BlockOffset(i);
  }
  else {
    i0 = i / _block_size;
    i1 = i % _block_size;
  }
  //    const int itmp = offsetBlock(i0, i0) + i1 * (nrowBlock(i0, i0) + 1);
  //    return coefs()[itmp];
  return _coefs[i0 * (_num_blocks + 1)][i1 * (nrowBlock(i0, i0) + 1)];
}

template
double& SquareBlockMatrix<double>::diag(int i);

template
quadruple& SquareBlockMatrix<quadruple>::diag(int i);

template
complex<double> & SquareBlockMatrix<complex<double> >::diag(int i);

template
complex<quadruple> & SquareBlockMatrix<complex<quadruple> >::diag(int i);
//

template<typename T>
const T& SquareBlockMatrix<T>::diag(int i) const
{
  int i0, i1;
  if(_isdecomposed) {
    fprintf(stderr, "%s %d : operator() is not suppoesd to be used\n",
	    __FILE__, __LINE__);
    i0 = BlockIndex(i);
    i1 = BlockOffset(i);
  }
  else {
    i0 = i / _block_size;
    i1 = i % _block_size;
  }
  return _coefs[i0 * (_num_blocks + 1)][i1 * (nrowBlock(i0, i0) + 1)];
}

template
const double & SquareBlockMatrix<double>::diag(int i) const;

template
const quadruple & SquareBlockMatrix<quadruple>::diag(int i) const;

template
const complex<double> & SquareBlockMatrix<complex<double> >::diag(int i) const;

template
const complex<quadruple> & SquareBlockMatrix<complex<quadruple> >::
diag(int i) const;
//

template<typename T>
T& SquareBlockMatrix<T>::operator () (int i, int j)
{ 
  int i0, i1, j0, j1;
  if(_isdecomposed) {
    fprintf(stderr, "%s %d : operator() is not suppoesd to be used\n",
	    __FILE__, __LINE__);
    if (_isSym) {
      if (i < j) {
	i0 = BlockIndex(i);
	j0 = BlockIndex(j);
	i1 = BlockOffset(i);
	j1 = BlockOffset(j);
      }
      else {
	i0 = BlockIndex(j);
	j0 = BlockIndex(i);
	i1 = BlockOffset(j);
	j1 = BlockOffset(i);
      }
    }
    else {
      i0 = BlockIndex(i);
      j0 = BlockIndex(j);
      if(i0 > j0) {
	i1 = BlockOffset(j);
	j1 = BlockOffset(i);
      }
      else {
	i1 = BlockOffset(i);
	j1 = BlockOffset(j);
      }
    }
  }  //     if(_isdecomposed) {
  else {
    if (_isSym) {
      if (i < j) {
	i0 = i / _block_size;
	i1 = i % _block_size;
	j0 = j / _block_size;
	j1 = j % _block_size;
	}
      else {
	i0 = j / _block_size;
	i1 = j % _block_size;
	j0 = i / _block_size;
	j1 = i % _block_size;
      }
    }
    else {
      i0 = i / _block_size;
      j0 = j / _block_size;
      if (i0 > j0) {
	i1 = j % _block_size;
	j1 = i % _block_size;
      }
      else {
	i1 = i % _block_size;
	j1 = j % _block_size;
      }
    }
  }  //  if(_isdecomposed) else
    return _coefs[i0 + j0 * _num_blocks][i1 + j1 * nrowBlock(i0, j0)];
}

template
double& SquareBlockMatrix<double>::operator () (int i, int j);

template
quadruple& SquareBlockMatrix<quadruple>::operator () (int i, int j);

template
complex<double> & SquareBlockMatrix<complex<double> >::
operator () (int i, int j);

template
complex<quadruple> & SquareBlockMatrix<complex<quadruple> >::
operator () (int i, int j);
//

template<typename T>
const T& SquareBlockMatrix<T>::operator () (int i, int j) const
{ 
  int i0, i1, j0, j1;
  if(_isdecomposed) {
    fprintf(stderr, "%s %d : operator() is not suppoesd to be used\n",
	    __FILE__, __LINE__);
    if (_isSym) {
      if (i < j) {
	i0 = BlockIndex(i);
	j0 = BlockIndex(j);
	i1 = BlockOffset(i);
	j1 = BlockOffset(j);
      }
      else {
	i0 = BlockIndex(j);
	j0 = BlockIndex(i);
	i1 = BlockOffset(j);
	  j1 = BlockOffset(i);
      }
    }
    else {
      i0 = BlockIndex(i);
      j0 = BlockIndex(j);
      if(i0 > j0) {
	i1 = BlockOffset(j);
	j1 = BlockOffset(i);
      }
      else {
	i1 = BlockOffset(i);
	  j1 = BlockOffset(j);
      }
    }
  }  //     if(_isdecomposed) {
  else {
    if (_isSym) {
      if (i < j) {
	i0 = i / _block_size;
	i1 = i % _block_size;
	j0 = j / _block_size;
	j1 = j % _block_size;
      }
      else {
	i0 = j / _block_size;
	i1 = j % _block_size;
	j0 = i / _block_size;
	  j1 = i % _block_size;
      }
    }
    else {
      i0 = i / _block_size;
      j0 = j / _block_size;
      if (i0 > j0) {
	i1 = j % _block_size;
	j1 = i % _block_size;
      }
      else {
	i1 = i % _block_size;
	j1 = j % _block_size;
	}
    }
  }  //  if(_isdecomposed) else
  return _coefs[i0 + j0 * _num_blocks][i1 + j1 * nrowBlock(i0, j0)];
}

template
const double& SquareBlockMatrix<double>::operator () (int i, int j) const;

template
const quadruple& SquareBlockMatrix<quadruple>::operator () (int i, int j) const;

template
const complex<double>& SquareBlockMatrix<complex<double> >::
operator () (int i, int j) const;

template
const complex<quadruple>& SquareBlockMatrix<complex<quadruple> >::
operator () (int i, int j) const;
//

template<typename T>
void SquareBlockMatrix<T>::copyFromArray(const T *a, int dim)
{
  //  fprintf(stderr, "%s %d dim = %d block_size = %d\n", __FILE__, __LINE__,
  //	  _dim, _block_size);
  if (_isdecomposed) {
    fprintf(stderr, "%s %d : copyFromArray() cannot copy decompesd data\n",
	    __FILE__, __LINE__);
    return;
  }
  if (dim != _dim) {
    fprintf(stderr, "%s %d : mismatched dim %d != %d\n", 
	    __FILE__, __LINE__, dim, _dim);
    return;
  }
  for (int ki = 0; ki < _num_blocks; ki++) {
    const int n_row = ((ki == (_num_blocks - 1)) ?
		       _block_size_last : _block_size);
    const int dst_row = n_row;
    for (int kj = ki; kj < _num_blocks; kj++) {
      const int n_col = ((kj == (_num_blocks - 1)) ?
			 _block_size_last : _block_size);
      T *a_upper = addrCoefBlock(ki, kj);
      int j, j0, j1;
      for (j1 = (ki + kj * _dim) * _block_size, j0 = 0, j = 0; 
	   j < n_col; j++, j1 += _dim, j0 += dst_row) {
	blas_copy(n_row, (T *)a + j1, 1, a_upper + j0, 1); 
      }
    } // loop : kj
  } // loop : ki
  if (!_isSym) {
    for (int ki = 0; ki < _num_blocks; ki++) {
      const int n_row = (ki == (_num_blocks - 1)) ?
	_block_size_last : _block_size;
      for (int kj = 0; kj < ki; kj++) {
	const int n_col = _block_size;
	T *a_lower = addrCoefBlock(ki, kj);
	int dst_row = nrowBlock(ki, kj);
	int j, j0, j1;
	for (j1 = (ki + kj * _dim) * _block_size, j0 = 0, j = 0; 
	     j < n_col; j++, j1 += _dim, j0++) {
	  blas_copy(n_row, (T *)a + j1, 1, a_lower + j0, dst_row); 
	}
      } // loop : kj
    } // loop : ki
  }  // if (!isSym)
}

template
void SquareBlockMatrix<double>::copyFromArray(const double *a, int dim);

template
void SquareBlockMatrix<quadruple>::copyFromArray(const quadruple *a, int dim);
template
void SquareBlockMatrix<complex<double> >::
copyFromArray(const complex<double> *a, int dim);

template
void SquareBlockMatrix<complex<quadruple> >::
copyFromArray(const complex<quadruple> *a, int dim);
//

template<typename T>
void SquareBlockMatrix<T>::copyFromArrayPermute(const T *a, int dim,
						int *permute)
{
  //  fprintf(stderr, "%s %d dim = %d block_size = %d\n", __FILE__, __LINE__,
  //	  _dim, _block_size);
  if (_isdecomposed) {
    fprintf(stderr, "%s %d : copyFromArray() cannot copy decompesd data\n",
	    __FILE__, __LINE__);
    return;
  }
  if (dim != _dim) {
    fprintf(stderr, "%s %d : mismatched dim %d != %d\n", 
	    __FILE__, __LINE__, dim, _dim);
    return;
  }
  for (int ki = 0; ki < _num_blocks; ki++) {
    const int n_row = ((ki == (_num_blocks - 1)) ?
		       _block_size_last : _block_size);
    const int dst_row = n_row;
    for (int kj = ki; kj < _num_blocks; kj++) {
      const int n_col =	((kj == (_num_blocks - 1)) ?
			 _block_size_last : _block_size);
      T *a_upper = addrCoefBlock(ki, kj);
      int j, j0, j1;
      for (j1 = (ki + kj * _dim) * _block_size, j0 = 0, j = 0; 
	   j < n_col; j++, j1 += _dim, j0 += dst_row) {
	const int jj = permute[kj * _block_size + j];
	int ii = ki * _block_size;
	int kk = j0;
	if (!_isSym) {
	  for (int i = 0; i < n_row; i++, kk++, ii++) {
	    a_upper[kk] = a[permute[ii] + jj * _dim];
	  }
	}
	else {
	  for (int i = 0; i < n_row; i++, kk++, ii++) {
	    const int itmp = permute[ii];
	    a_upper[kk] = a[itmp <= jj ? itmp + jj * _dim : jj + itmp * _dim];
	  }
	} // if (!isSym) {
      //	blas_copy(n_row, (T *)a + j1, 1, a_upper + j0, 1);
      }
    } // loop : kj
  } // loop : ki
  if (!_isSym) {
    for (int ki = 0; ki < _num_blocks; ki++) {
      const int n_row = ((ki == (_num_blocks - 1)) ?
			 _block_size_last : _block_size);
      for (int kj = 0; kj < ki; kj++) {
	const int n_col = _block_size;
	T *a_lower = addrCoefBlock(ki, kj);
	const int dst_row = nrowBlock(ki, kj);
	int j, j0, j1;
	for (j1 = (ki + kj * _dim) * _block_size, j0 = 0, j = 0; 
	     j < n_col; j++, j1 += _dim, j0++) {
	  const int jj = permute[kj * _block_size + j];
	  int ii = ki * _block_size;
	  int kk = j0;
	  for (int i = 0; i < n_row; i++, ii++, kk += dst_row) {
	    a_lower[kk] = a[permute[ii] + jj * _dim];
	  }
	  // 	  blas_copy(n_row, (T *)a + j1, 1, a_lower + j0, dst_row);
	}
      } // loop : kj
    } // loop : ki
  }  // if (!isSym)
}

template
void SquareBlockMatrix<double>::copyFromArrayPermute(const double *a, int dim,
						     int *permute);

template
void SquareBlockMatrix<quadruple>::copyFromArrayPermute(const quadruple *a,
							int dim,
							int *permute);

template
void SquareBlockMatrix<complex<double> >::
copyFromArrayPermute(const complex<double> *a, int dim,
		     int *permute);

template
void SquareBlockMatrix<complex<quadruple> >::
copyFromArrayPermute(const complex<quadruple> *a, int dim,
		     int *permute);
//

template<typename T> 
void SquareBlockMatrix<T>::copyToArrayFull(const T *a, int dim)
{
  if (dim != _dim) {
    fprintf(stderr, "%s %d : mismatched dim %d != %d\n", 
	    __FILE__, __LINE__, dim, _dim);
    return;
  }
  // copy upper part
  for (int ki = 0; ki < _num_blocks; ki++) {
    const int n_row = ((ki == (_num_blocks - 1)) ?
		       _block_size_last : _block_size);
    const int dst_row = n_row;
    for (int kj = ki; kj < _num_blocks; kj++) {
      const int n_col = ((kj == (_num_blocks - 1)) ?
			 _block_size_last : _block_size);
      T *a_upper = addrCoefBlock(ki, kj);
      int j, j0, j1;
      for (j1 = (ki + kj * _dim) * _block_size, j0 = 0, j = 0; 
	   j < n_col; j++, j1 += _dim, j0 += dst_row) {
	blas_copy(n_row, a_upper + j0, 1, (T *)a + j1, 1); 
      }
    } // loop : kj
  } // loop : ki
  // copy lower part
  if (_isSym) {
    for (int ki = 0; ki < _num_blocks; ki++) {
      const int n_row = ((ki == (_num_blocks - 1)) ?
			 _block_size_last : _block_size);
      for (int kj = 0; kj < ki; kj++) {
	const int n_col = _block_size;
	const int dst_col = n_col;
	T *a_upper = addrCoefBlock(kj, ki); // transposed access
	int j, j0, j1;
	for (j1 = (ki + kj * _dim) * _block_size, j0 = 0, j = 0; 
	     j < n_row; j++, j1++, j0 += dst_col) {
	  blas_copy(n_col, a_upper + j0, 1, (T *)a + j1, _dim);
	}
      } // loop : kj
    } // loop : ki
  }
  else {
    for (int ki = 0; ki < _num_blocks; ki++) {
      const int n_row = ((ki == (_num_blocks - 1)) ?
			 _block_size_last : _block_size);
      for (int kj = 0; kj < ki; kj++) {
	const int n_col = _block_size;
	T *a_upper = addrCoefBlock(ki, kj);
	const int src_row = nrowBlock(ki, kj);
	int j, j0, j1;
	for (j1 = (ki + kj * _dim) * _block_size, j0 = 0, j = 0; 
	     j < n_col; j++, j1 += _dim, j0++) {
	  blas_copy(n_row, a_upper + j0, src_row, (T *)a + j1, 1);
	}
      } // loop : kj
    } // loop : ki
  }
}

template
void SquareBlockMatrix<double>::copyToArrayFull(const double *a, int dim);

template
void SquareBlockMatrix<quadruple>::copyToArrayFull(const quadruple *a, int dim);
template
void SquareBlockMatrix<complex<double> >::
copyToArrayFull(const complex<double> *a, int dim);


template
void SquareBlockMatrix<complex<quadruple> >::
copyToArrayFull(const complex<quadruple> *a, int dim);
//

template<typename T>
void SquareBlockMatrix<T>::copyBlockToArray(int i, int j, const T *a, int dim)
{
  if (dim != _dim) {
  fprintf(stderr, "%s %d : mismatched dim %d != %d\n", 
	    __FILE__, __LINE__, dim, _dim);
    return;
  }
  if (_isSym && i > j) {
    fprintf(stderr, 
	    "%s %d : symmetric matrix has upper block only : %d > %d\n", 
	    __FILE__, __LINE__, i, j);
    return;
  }
  T *a_src = addrCoefBlock(i, j);
  int j0, j1, j2;
  if (i <= j) {
    const int n_row = ((i == (_num_blocks - 1)) ?
		       _block_size_last : _block_size);
    const int dst_row = n_row;
    const int n_col = ((j == (_num_blocks - 1)) ?
		       _block_size_last : _block_size);
    for (j2 = (i + j * _dim) * _block_size, j1 = 0, j0 = 0;
	 j0 < n_col; j0++, j1 += dst_row, j2 += _dim) {
      blas_copy(n_row, a_src + j1, 1, (T *)a + j2, 1);
    }
  }
  else {
    const int n_row = ((i == (_num_blocks - 1)) ?
		       _block_size_last : _block_size);
    const int n_col = ((j == (_num_blocks - 1)) ?
		       _block_size_last : _block_size);
    const int dst_col = n_col;
    for (j2 = (i + j * _dim) * _block_size, j1 = 0, j0 = 0;
	 j0 < n_col; j0++, j1++, j2 += _dim) {
      blas_copy(n_row, a_src + j1, dst_col, (T *)a + j2, 1);
    }
  }
}

template
void SquareBlockMatrix<double>::copyBlockToArray(int i, int j, 
						 const double *a, int dim);

template
void SquareBlockMatrix<quadruple>::copyBlockToArray(int i, int j, 
						    const quadruple *a,
						    int dim);
template
void SquareBlockMatrix<complex<double> >::
copyBlockToArray(int i, int j, const complex<double> *a, int dim);

template
void SquareBlockMatrix<complex<quadruple> >::
copyBlockToArray(int i, int j, const complex<quadruple> *a, int dim);
