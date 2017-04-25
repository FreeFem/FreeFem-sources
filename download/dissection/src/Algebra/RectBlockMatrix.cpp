/*! \file   RectBlockMatrix.cpp
    \brief  Block storage for Rectangular matrix
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

#include "Algebra/RectBlockMatrix.hpp"
#include "Compiler/blas.hpp"
#include <cstdio>

template<typename T>
void RectBlockMatrix<T>::init(int dim_r, int dim_c, int block_size, int first_block)
{
  int itmp;

  _dim_r = dim_r;
  _dim_c = dim_c;
  _block_size = block_size;

  _block_size2 = _block_size * _block_size;

  if ((dim_r * dim_c) == 0) {
    _isdecomposed_c = false;
    _num_blocks_r = 0;
    _num_blocks_c = 0;
    _block_size_last_r = 0;
    _block_size_last_c = 0;
    _coefs = NULL; 
    _allocation_status = NULL; 
    _deallocation_status = false;
    return;
  }
  if (first_block > 0) {
    _isdecomposed_c = true;
    _dim0_c = first_block;
    _dim1_c = _dim_c - _dim0_c;
    _num_blocks_r = (_dim_r + _block_size - 1) / _block_size;
    _num_blocks0_c = (_dim0_c + _block_size - 1) / _block_size;
    _num_blocks1_c = (_dim1_c + _block_size - 1) / _block_size;
    _num_blocks_c = _num_blocks0_c + _num_blocks1_c;
    _block_sizes_r.resize(_num_blocks_r);
    _block_sizes_c.resize(_num_blocks_c);
    _indexblock_r.resize(_num_blocks_r + 1);
    _indexblock_c.resize(_num_blocks_c + 1);

    itmp = _dim_r % _block_size;
    _block_size_last_r = (itmp == 0) ? _block_size : itmp;
    int ktmp;
    ktmp = 0;
    if (_num_blocks_r > 0) {
      for (int i = 0; i < (_num_blocks_r - 1); i++, ktmp += _block_size) {
	_block_sizes_r[i] = _block_size;
	_indexblock_r[i] = ktmp;
      }
      _block_sizes_r[_num_blocks_r - 1] = _block_size_last_r;
      _indexblock_r[_num_blocks_r - 1] = ktmp;
      _indexblock_r[_num_blocks_r] = _dim_r;
    }
    itmp = _dim0_c % _block_size;
    _block_size_last0_c = (itmp == 0) ? _block_size : itmp;
    ktmp = 0;
    if (_num_blocks0_c > 0) {
      for (int i = 0; i < (_num_blocks0_c - 1); i++, ktmp += _block_size) {
	_block_sizes_c[i] = _block_size;
	_indexblock_c[i] = ktmp;
      }
      _block_sizes_c[_num_blocks0_c - 1] = _block_size_last0_c;
      _indexblock_c[_num_blocks0_c - 1] = ktmp;
    }
    itmp = _dim1_c % _block_size;
    _block_size_last1_c = (itmp == 0) ? _block_size : itmp;
    ktmp = _dim0_c;
    if (_num_blocks1_c > 0) {
      for (int i = _num_blocks0_c; i < (_num_blocks_c - 1); 
	   i++, ktmp += _block_size) {
	_block_sizes_c[i] = _block_size;
	_indexblock_c[i] = ktmp;
      }
      _block_sizes_c[_num_blocks_c - 1] = _block_size_last1_c;
      _indexblock_c[_num_blocks_c - 1] = ktmp;
    }
    _indexblock_c[_num_blocks_c] = _dim_c;
  }
  else {
    _isdecomposed_c = false;
    _num_blocks_r = (_dim_r + _block_size - 1) / _block_size;
    _num_blocks_c = (_dim_c + _block_size - 1) / _block_size;
    _block_sizes_r.resize(_num_blocks_r);
    _block_sizes_c.resize(_num_blocks_c);
    _indexblock_r.resize(_num_blocks_r + 1);
    _indexblock_c.resize(_num_blocks_c + 1);

    itmp = _dim_r % _block_size;
    _block_size_last_r = (itmp == 0) ? _block_size : itmp;
    int ktmp;
    ktmp = 0;
    if (_num_blocks_r > 0) {
      for (int i = 0; i < (_num_blocks_r - 1); i++, ktmp += _block_size) {
	_block_sizes_r[i] = _block_size;
	_indexblock_r[i] = ktmp;
      }
      _block_sizes_r[_num_blocks_r - 1] = _block_size_last_r;
      _indexblock_r[_num_blocks_r - 1] = ktmp;
      _indexblock_r[_num_blocks_r] = _dim_r;
    }
    itmp = _dim_c % _block_size;
    ktmp = 0;
    if (_num_blocks_c > 0) {
      _block_size_last_c = (itmp == 0) ? _block_size : itmp;
      for (int i = 0; i < (_num_blocks_c - 1); i++, ktmp += _block_size) {
	_block_sizes_c[i] = _block_size;
	_indexblock_c[i] = ktmp;
      }
      _block_sizes_c[_num_blocks_c - 1] = _block_size_last_c;
      _indexblock_c[_num_blocks_c - 1] = ktmp;
      _indexblock_c[_num_blocks_c] = _dim_c;
    }
  }
  if (_num_blocks_r * _num_blocks_c == 0) {
    _coefs = NULL; // new T*[1]; // dummy
    _allocation_status = NULL; // new bool[1]; // dummy
    _deallocation_status = false;
    //      _deallocation_status = true;
  }
  else {
    try {
      _coefs = new T*[_num_blocks_r * _num_blocks_c];  // need to be managed
    } catch (const std::bad_alloc& e) {
      fprintf(stderr, "%s %d : allocation failed : %s", __FILE__, __LINE__,
	      e.what());
    }
    try {
      _allocation_status = new bool[_num_blocks_r * _num_blocks_c];
    } catch (const std::bad_alloc& e) {
      fprintf(stderr, "%s %d : allocation failed : %s", __FILE__, __LINE__,
	      e.what());
    }      
    for (int i = 0; i < (_num_blocks_r * _num_blocks_c); i++) {
      _allocation_status[i] = false;
    }
    _deallocation_status = false;
  }
}

template
void RectBlockMatrix<double>::init(int dim_r, int dim_c, int block_size,
				   int first_block);
template
void RectBlockMatrix<quadruple>::init(int dim_r, int dim_c, int block_size,
				      int first_block);
template
void RectBlockMatrix<complex<double> >::init(int dim_r, int dim_c,
					     int block_size, int first_block);
template
void RectBlockMatrix<complex<quadruple> >::init(int dim_r, int dim_c,
						int block_size,
						int first_block);
//
  
template<typename T>
void RectBlockMatrix<T>::allocateBlock(int i, int j)
{
  const int nrow = _block_sizes_r[i];
  const int ncol = _block_sizes_c[j];
  if (!_allocation_status[i + j * _num_blocks_r]) {   //
    try {
      _coefs[i + j * _num_blocks_r] = new T[nrow * ncol];
    } catch (const std::bad_alloc& e) {
      fprintf(stderr, "%s %d : allocation failed : %s", __FILE__, __LINE__,
	      e.what());
    }      
  }
  _allocation_status[i + j * _num_blocks_r] = true;
}

template
void RectBlockMatrix<double>::allocateBlock(int i, int j);
template
void RectBlockMatrix<quadruple>::allocateBlock(int i, int j);
template
void RectBlockMatrix<complex<double> >::allocateBlock(int i, int j);
template
void RectBlockMatrix<complex<quadruple> >::allocateBlock(int i, int j);
//

template<typename T>
void RectBlockMatrix<T>::allocate()
{
  for (int i = 0; i < _num_blocks_r; i++) {
    for (int j = 0; j < _num_blocks_c; j++) {
      const int nrow = _block_sizes_r[i];
      const int ncol = _block_sizes_c[j];
      if (!_allocation_status[i + j * _num_blocks_r]) {
	try {
	  _coefs[i + j * _num_blocks_r] = new T[nrow * ncol];
	} catch (const std::bad_alloc& e) {
	  fprintf(stderr, "%s %d : allocation failed : %s", __FILE__, __LINE__,
		  e.what());
	}      
      }
      _allocation_status[i + j * _num_blocks_r] = true;
    }
  }
  _deallocation_status = false;
}

template
void RectBlockMatrix<double>::allocate();
template
void RectBlockMatrix<quadruple>::allocate();
template
void RectBlockMatrix<complex<double> >::allocate();
template
void RectBlockMatrix<complex<quadruple> >::allocate();
//

template<typename T>
void RectBlockMatrix<T>::free(int i, int j)
{
    //    fprintf(stderr, "%s %d : free()\n", __FILE__, __LINE__);
  if (_allocation_status[i + j * _num_blocks_r]) {
    delete[] _coefs[i + j * _num_blocks_r];  // need to be mananged
  }
  _allocation_status[i + j * _num_blocks_r] = false;
}

template
void RectBlockMatrix<double>::free(int i, int j);
template
void RectBlockMatrix<quadruple>::free(int i, int j);
template
void RectBlockMatrix<complex<double> >::free(int i, int j);
template
void RectBlockMatrix<complex<quadruple> >::free(int i, int j);
//

template<typename T>
void RectBlockMatrix<T>::free()
{
  //    fprintf(stderr, "%s %d : free()\n", __FILE__, __LINE__);
  for (int i = 0; i < _num_blocks_r; i++) {
    for (int j = 0; j < _num_blocks_c; j++) {
      if (_allocation_status[i + j * _num_blocks_r]) {
	delete[] _coefs[i + j * _num_blocks_r];  // need to be managed
	_allocation_status[i + j * _num_blocks_r] = false;
      }
    }
  }
}

template
void RectBlockMatrix<double>::free();
template
void RectBlockMatrix<quadruple>::free();
template
void RectBlockMatrix<complex<double> >::free();
template
void RectBlockMatrix<complex<quadruple> >::free();
//

template<typename T>
T* RectBlockMatrix<T>::addrCoefBlock(int i, int j) 
{
  return _coefs[i + j * _num_blocks_r];
}

template
double* RectBlockMatrix<double>::addrCoefBlock(int i, int j);

template
quadruple* RectBlockMatrix<quadruple>::addrCoefBlock(int i, int j);

template
complex<double>* RectBlockMatrix<complex<double> >::
addrCoefBlock(int i, int j);

template
complex<quadruple>* RectBlockMatrix<complex<quadruple> >::
addrCoefBlock(int i, int j);
//

template<typename T>
void RectBlockMatrix<T>::ZeroClear()
{
  const T zero(0.0);
  for (int i = 0; i < _num_blocks_r; i++) {
    for (int j = 0; j < _num_blocks_c; j++) {
      const int nrow = _block_sizes_r[i];
      const int ncol = _block_sizes_c[j];
      const int nsize = nrow * ncol;
      const int itmp = i + j * _num_blocks_r;
      for (int k = 0; k < nsize; k++) {
	_coefs[itmp][k] = zero;   // need to be replaced by memcpy
      }
    }
  }
}

template
void RectBlockMatrix<double>::ZeroClear();

template
void RectBlockMatrix<quadruple>::ZeroClear();

template
void RectBlockMatrix<complex<double> >::ZeroClear();

template
void RectBlockMatrix<complex<quadruple> >::ZeroClear();
//

template<typename T>
T& RectBlockMatrix<T>::operator () (int i, int j)
{ 
  int i0, i1, j0, j1;

  if (_isdecomposed_c) {
    i0 = i / _block_size;
    i1 = i % _block_size;
    if (j < _dim0_c) {
      j0 = j / _block_size;
      j1 = j % _block_size;
    }
    else {
      j0 = (j - _dim0_c) / _block_size + _num_blocks0_c;
      j1 = (j - _dim0_c) % _block_size;
    }
  }
  else {
    i0 = i / _block_size;
    j0 = j / _block_size;
    i1 = i % _block_size;
    j1 = j % _block_size;
  }
  return _coefs[i0 + j0 * _num_blocks_r][i1 + j1 * _block_sizes_r[i0]];
}

template
double& RectBlockMatrix<double>::operator () (int i, int j);

template
quadruple& RectBlockMatrix<quadruple>::operator () (int i, int j);

template
complex<double> & RectBlockMatrix<complex<double> >::
operator () (int i, int j);

template
complex<quadruple> & RectBlockMatrix<complex<quadruple> >::
operator () (int i, int j);
//

template<typename T>
const T& RectBlockMatrix<T>::operator () (int i, int j) const
{ 
  int i0, i1, j0, j1;

  i0 = i / _block_size;
  j0 = j / _block_size;
  i1 = i % _block_size;
  j1 = j % _block_size;
  return _coefs[i0 + j0 * _num_blocks_r][i1 + j1 * _block_sizes_r[i0]];
}

template
const double& RectBlockMatrix<double>::operator () (int i, int j) const;

template
const quadruple& RectBlockMatrix<quadruple>::operator () (int i, int j) const;

template
const complex<double>& RectBlockMatrix<complex<double> >::
operator () (int i, int j) const;

template
const complex<quadruple>& RectBlockMatrix<complex<quadruple> >::
operator () (int i, int j) const;
//
