/*! \file   RectBlockMatrix.hpp
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

#ifndef _ALGEBRA_RECTBLOCKMATRIX_HPP
#define _ALGEBRA_RECTBLOCKMATRIX_HPP
#include <cstdio>
#include <cassert>
#include <cmath>
#include <vector>

using std::vector;

template<typename T>
class RectBlockMatrix 
{
public: 

  RectBlockMatrix()
  { 
    _dim_r = 0;
    _dim_c = 0;
    _block_size = 0; 
    _block_size_last_r = 0;
    _block_size_last_c = 0;
    _num_blocks_r = 0;
    _num_blocks_c = 0;
    _deallocation_status = true; 
  }

  RectBlockMatrix(int dim_r, int dim_c, int block_size)
  {
#ifdef DEBUG_SQUAREBLOCKMATRIX
    fprintf(stderr, "%s %d : constructor with args %d %d\n", 
	    __FILE__, __LINE__, dim, block_size);
#endif
    init(dim_r, dim_c, block_size, 0);
    allocate();
    _deallocation_status = false;
  }

  RectBlockMatrix(int dim_r, int dim_c, double *aa)
  {
#ifdef DEBUG_SQUAREBLOCKMATRIX
    fprintf(stderr, "%s %d : constructor with args %d %d\n", 
	    __FILE__, __LINE__, dim, block_size);
#endif
    _block_size = dim_r > dim_c ? dim_r : dim_c;
    _dim_c = dim_c;
    _dim_r = dim_r;
    _num_blocks_r = 1;
    _num_blocks_c = 1;

    _coefs = new T*[1];
    _coefs[0] = aa;
    _block_sizes_r.resize(1);
    _block_sizes_c.resize(1);
    _indexblock_r.resize(2);
    _indexblock_c.resize(2);
    _block_sizes_r[0] = _dim_r;
    _block_sizes_c[0] = _dim_c;
    _indexblock_r[0] = 0;
    _indexblock_r[1] = _dim_r;
    _indexblock_c[0] = 0;
    _indexblock_c[1] = _dim_c;
    _allocation_status = new bool[1];
    _allocation_status[0] = false;   // _coefs[0] is not allowed to free
    _isdecomposed_c = false;
    _deallocation_status = true;
  }

  void init(int dim_r, int dim_c, int block_size, int first_block = 0);
  void allocateBlock(int i, int j);
  void allocate();
  void free(int i, int j);
  void free();
  
  ~RectBlockMatrix() { 
    free();
    if (!_deallocation_status) {
      delete [] _coefs;
      delete [] _allocation_status;
    }
    _deallocation_status = true;
    _block_sizes_r.clear();
    _block_sizes_c.clear();
    _num_blocks_r = 0; // ? to avoid double free
    _num_blocks_c = 0; // ? to avoid double free
    _dim_r = 0;
    _dim_c = 0;
  }

  int dimension_r() const { return (int)_dim_r; }
  int dimension_c() const { return (int)_dim_c; }
  int block_size() const { return (int)_block_size; }
  int num_blocks_r() const { return _num_blocks_r; }
  int num_blocks_c() const { return _num_blocks_c; }
  int nbRows() const { return _dim_r; }
  int nbColumns() const { return _dim_c; }

  int block_size_last_r() const { return (int)_block_size_last_r; }
  int block_size_last_c() const { return (int)_block_size_last_c; }

  int nrowBlock(int i) const { return _block_sizes_r[i]; }
  int ncolBlock(int i) const { return _block_sizes_c[i]; }

  int IndexBlock_r(int i) const { return _indexblock_r[i]; }
  int IndexBlock_c(int i) const { return _indexblock_c[i]; }
  T* addrCoefBlock(int i, int j);

  void ZeroClear();

  T& operator () (int i, int j);

  const T& operator () (int i, int j) const;

  RectBlockMatrix<T>* clone() const
  {
    RectBlockMatrix<T> *ret = new RectBlockMatrix<T>;
    ret->copy(*this);
    return(ret);
  }

  void copy ( const RectBlockMatrix& A )
  {
    _block_size = A._block_size;
    _block_size2 = A._block_size2;
    _num_blocks = A._num_blocks;
    _block_sizes_r = A._block_sizes_r;
    _block_sizes_c = A._block_sizes_c;
    _dim_r = A._dim_r;
    _dim_c = A._dim_c;
    _block_size_last_r = A._block_size_last_r;
    _block_size_last_c = A._block_size_last_c;
    _num_blocks_r = A._num_blocks_r;
    _num_blocks_c = A._num_blocks_c;
    _coefs = A._coefs;
    _allocation_status = A._allocation_status;
    _deallocation_status = A._deallocation_status;
  }

private:
  int _block_size;
  int _block_size2;
  int _num_blocks;
  vector<int> _block_sizes_r;
  vector<int> _block_sizes_c;
  vector<int> _indexblock_r;
  vector<int> _indexblock_c;
  int _dim_r;
  int _dim_c;
  int _block_size_last_r;
  int _block_size_last_c;
  int _block_size_last0_c;
  int _block_size_last1_c;
  int _num_blocks_r;
  int _num_blocks_c;
  int _num_blocks0_c;
  int _num_blocks1_c;
  int _dim0_c;
  int _dim1_c;
  T**   _coefs;
  bool *_allocation_status;
  bool _deallocation_status;
  bool _isdecomposed_c;
}; // End class RectMatrix

#endif
