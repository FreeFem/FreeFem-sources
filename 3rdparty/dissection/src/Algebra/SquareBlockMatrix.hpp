/*! \file   SquareBlockMatrix.hpp
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

#ifndef _ALGEBRA_SQUAREBLOCKMATRIX_HPP
#define _ALGEBRA_SQUAREBLOCKMATRIX_HPP
#include <cstdio>
#include <cmath>
#include <vector>

using std::vector;

template<typename T>
class SquareBlockMatrix 
{
public: 

  SquareBlockMatrix() : _dim(0), _size(0), _block_size(0), _block_size2(0),
			_block_size_last(0), _isSym(true)
  { 
    _isblocked = false;
    _isdecomposed = false;
    _num_blocks = 0;
    _num_blocks0 = 0;
    _num_blocks1 = 0;
    _dim0 = 0;
    _dim1 = 0;
    _block_size_last0 = 0;
    _block_size_last1 = 0;
    _lower_allocated = false;
    _coefs_alloc_status = false;
    _pivrelaxed = false;
  }

  SquareBlockMatrix(int dim, int block_size, bool isSym) 
  {
#ifdef DEBUG_SQUAREBLOCKMATRIX
    fprintf(stderr, "%s %d : constructor with args %d %d\n", 
	    __FILE__, __LINE__, dim, block_size);
#endif
    if (dim > 0) {
      init(dim, block_size, isSym, 0);
      allocate();
      _permute.resize(dim);
    }
  }

  void init(int dim, int block_size, bool isSym, int first_block = 0);
  void allocateBlock(int i, int j);
  void allocate();
  void free(int i, int j);
  void free();

  ~SquareBlockMatrix() {
    free();
    if (_coefs_alloc_status) {
      delete [] _coefs;
      delete [] _allocation_status;
      _coefs_alloc_status = false;
    }
    _block_sizes.clear();
    _nsing_block.clear();
    _num_blocks = 0; // ? to avoid double free
    _dim = 0;
    _permute.clear();
    _singIdx.clear();
    _singIdx0.clear();
  }

  int dimension() const { return (int)_dim; }
  int dimension0() const { return (int)_dim0; }
  int dimension1() const { return (int)_dim1; }
  int storage_size() const { return (int)_size; }
  int block_size() const { return (int)_block_size; }
  int num_blocks() const { return _num_blocks; }
  int block_size_last() const { return (int)_block_size_last; }
  bool isSym() const { return _isSym; }
  bool isDecomposed() const { return _isdecomposed; }
  int num_blocks0() const { return _num_blocks0; }
  int num_blocks1() const { return _num_blocks1; }
  int block_size_last0() const { return (int)_block_size_last0; }
  int block_size_last1() const { return (int)_block_size_last1; }
  int nrowBlock(int i) const { return _block_sizes[i]; }

  //  int offsetBlock(int i, int j) const;

  T** coefs() { return _coefs; }
  T* addrCoefBlock(int i, int j);

  void allocateLowerBlocks() {
    if ((!_lower_allocated) && _isSym) {
      for (int i = 0; i < _num_blocks; i++) {
	for (int j = 0; j < i; j++) {
	  const int nrow = _block_sizes[i];
	  const int ncol = _block_sizes[j];
	  const int itmp = i + j * _num_blocks;
	  _coefs[itmp] = new T[nrow * ncol];
	  _allocation_status[itmp] = true; 
	}
      }
      _lower_allocated = true;
    }
  }
  
  void freeLowerBlocks() {
#ifdef ALLOCATE_LOWER_BLOCKS
    if (_lower_allocated && _isSym) {
#else
    if (_isSym) {
#endif
      for (int i = 0; i < _num_blocks; i++) {
	for (int j = 0; j < i; j++) {
	  const int itmp = i + j * _num_blocks;
	  if (_allocation_status[itmp]) {
	    delete [] _coefs[itmp];
	  }
	  _allocation_status[itmp] = false;
	}
      }
      _lower_allocated = false;
    }
    else {
      fprintf(stderr, "double free _lower %s %d\n", __FILE__, __LINE__);
    }
  }

  int nrowBlock(int i, int j) const
  {
    return (i <= j ? _block_sizes[i] : _block_sizes[j]);
  }

  int ncolBlock(int i, int j) const
  {
    return (i <= j ? _block_sizes[j] : _block_sizes[i]);
  }

  bool isTransposed(int i, int j);

  int BlockIndex(int i) const {
    int iblock;
    if (_isdecomposed) {
      if (i < _dim0) {
	iblock = i / _block_size;
      }
      else {
	iblock = (i - _dim0) / _block_size + _num_blocks0;
      }
    }
    else {
      iblock = i / _block_size;
    }
    return iblock;
  }

  int IndexBlock(int iblock) {
    int index;
    if (_isdecomposed) {
      if (iblock < _num_blocks0) {
	index = iblock * _block_size;
      }
      else {
	index = _dim0 + (iblock - _num_blocks0) * _block_size;
      }
    }
    else {
      index = iblock * _block_size;
    }
    return index;
  }

  int BlockOffset(int i) const {
    int ires;
    if (_isdecomposed) {
      if (i < _dim0) {
	ires = i % _block_size;
      }
      else {
	ires = (i - _dim0) % _block_size;
      }
    }
    else {
      ires = i % _block_size;
    }
    return ires;
  }
  void copyFromArrayPermute(const T *a, int dim, int *permute);
  void copyFromArray(const T *a, int dim);
  void copyToArrayFull(const T *a,int dim);
  void copyBlockToArray(int i, int j, const T *a,int dim);

  void ZeroClear();

  T& diag(int i);

  const T& diag(int i) const;

  T& operator () (int i, int j);

  const T& operator () (int i, int j) const;

  SquareBlockMatrix<T>* clone() const
  {
    SquareBlockMatrix<T> *ret = new SquareBlockMatrix<T>;
    ret->copy(*this);
    return(ret);
  }

  void copy ( const SquareBlockMatrix& A )
  {
    //    PlainMatrix<T>::copy(A); //_loc2glob=A._loc2glob;
    _block_size = A._block_size;
    _block_size2 = A._block_size2;
    _num_blocks = A._num_blocks;
    _block_size_last = A._block_size_last;
    _isSym = A._isSym;
    _rank = A._rank;
    _nsing = A._nsing;
    _nsing_block = A._nsing_block;
    _dim = A._dim;
    _size = A._size;
    _permute = A._permute;
    _kernelDetected = A._kernelDetected;
    _singIdx = A._singIdx;
    _singIdx0 = A._singIdx0;
    _isblocked = A._isblocked;
  }

  bool isBlocked()
  {
    return _isblocked;
  }
  void unsetBlocked()
  {
    _isblocked = false;
  }
  void setBlocked()
  {
    _isblocked = true;
  }

  void set_rank(int rank) { 
    _rank = rank; 
    _nsing = _dim - _rank;
  }

  void set_KernelDetected(int state) {
    _kernelDetected = state;
  }

  int KernelDetected() const {
    return _kernelDetected;
  }

  vector<int> &getNsingBlock() {
    return _nsing_block;
  }

  vector<int> &getSingIdx() {
    return _singIdx;
  }

  vector<int> &getSingIdx0() {
    return _singIdx0;
  }

  int rank() const
  {
    return _rank;
  }

  int dim_kern() const 
  {
    return _nsing;
  }
  void set_dim_kern(int nsing) 
  { _nsing = nsing; }

  int dim_kern_block(int k) const
  {
    return _nsing_block[k];
  }
  void set_dim_kern_block(int k, int nsing)
  {
    _nsing_block[k] = nsing;
  }

  vector<int>& getPermute()
  {
    return _permute;
  }
  const vector<int>& getPermute() const
  {
    return _permute;
  }

 void set_lastPivot(double lastpiv)
  {
    _lastpiv = lastpiv;
  }

  double lastPivot() const
  {
# ifdef DISDEBUG
    assert( isFactorized() );
# endif
    return _lastpiv;
  }

  void unset_pivrelaxed() { _pivrelaxed = false; }
  void set_pivrelaxed() { _pivrelaxed = true; }
  bool is_pivrelaxed() { return _pivrelaxed; }
private:
  int _dim;   //<! Dimension of the matrix
  int _size;
  int _block_size;
  int _block_size2;
  int _num_blocks;
  int _block_size_last;
  bool _isSym;
  vector<int> _block_sizes;
  int _dim0;
  int _dim1;
  int _block_size_last0;
  int _block_size_last1;
  int _num_blocks0;
  int _num_blocks1;
  bool _isdecomposed;
  //
  bool _coefs_alloc_status;
  bool *_allocation_status;
  bool _lower_allocated;
  //
  int  _rank;     //!< Rank of the matrix (initialized after factorization)
  int _nsing;
  double _lastpiv; // in case of quadruple ?
  T**   _coefs;
  //
  bool _isblocked;
  vector<int> _nsing_block;
  vector<int> _permute;
  bool        _kernelDetected;  // true : done, false : unknown
  vector<int> _singIdx;  // with permutation
  vector<int> _singIdx0; // without permutation
  bool _pivrelaxed;
  // 
}; // End class SquareMatrix

#endif
