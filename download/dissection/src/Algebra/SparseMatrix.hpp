/*! \file   SparseMatrix.hpp
    \brief  Sparse matrix definition
    \author Xavier Juvigny, ONERA
    \date   Jan. 25th 2005
    \modification allocation of array by STL vector class
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

#ifndef _ALGEBRA_SPARSEMATRIX_HPP_
#define _ALGEBRA_SPARSEMATRIX_HPP_

#include <vector>
#include "Compiler/blas.hpp"

using std::vector;

template<typename T, typename W = T, typename Z = W>
class SparseMatrix
{
public:

  SparseMatrix() : _isSymmetric(false), _isUpper(true), _isWhole(false) {}

  SparseMatrix(bool isSym, bool isUpper = true, bool isWhole = false) : 
    _isSymmetric(isSym), _isUpper(isUpper), _isWhole(isWhole)
  {}

  SparseMatrix( int n, int nnz, const T* coefs,
		bool isSym = false, 
		bool isUpper = true,
		bool isWhole = false) : 
    _isSymmetric(isSym), _isUpper(isUpper), _isWhole(isWhole)
  {
    _ptRows.resize(n+1);
    _indCols.resize(nnz);
    _coefs.resize(nnz);
    for (int i = 0 ; i < nnz; i++) {
      _coefs[i] = coefs[i];
    }
  }

  // added by Atsushi, 01 Feb.2013 : symbolic data are given by integer arrays
  SparseMatrix( int n, int nnz, const int *ptRows, const int *indCols,
		bool isSym = false, bool isUpper = true,
		bool isWhole = false) : 
    _isSymmetric(isSym),
    _isUpper(isUpper),
    _isWhole(isWhole)
  {
    _ptRows.resize(n+1);
    for (int i = 0; i < (n + 1); i++) {
      _ptRows[i] = ptRows[i];
    }
    _indCols.resize(nnz);
    for (int i = 0; i < nnz; i++) {
      _indCols[i] = indCols[i];
    }
    _coefs.resize(nnz);
  }

  SparseMatrix( int n, int nnz, const vector<int>& ptRows, 
		bool isSym = false, 
		bool isUpper = true,
		bool isWhole = false) : 
    _isSymmetric(isSym), _isUpper(isUpper), _isWhole(isWhole)
  {
    _indCols.resize(nnz);
    _ptRows.resize(n + 1);
    for (int i = 0; i < (n + 1); i++) {
      _ptRows[i] = ptRows[i];
    }
    _coefs.resize(nnz);
  }

  SparseMatrix( const SparseMatrix& A ) :
    _ptRows(A._ptRows), _indCols(A._indCols), _coefs(A._coefs), 
    _isSymmetric(A._isSymmetric), _isUpper(A._isUpper), _isWhole(A._isWhole)
  { }

  ~SparseMatrix() { free(); }

  void free() {
    _ptRows.clear();
    _indCols.clear();
    _coefs.clear();
  }

  inline int dimension() const { return (_ptRows.size() - 1); }

  inline int nnz() const { return _indCols.size(); }

  inline bool isSymmetric() const { return _isSymmetric; }

  inline bool isUpper() const { return _isUpper; }

  inline bool isWhole() const { return _isWhole; }
  
  //
  inline int *getRows() { return &_ptRows[0]; }
  //
  inline int *getIndCols() { return &_indCols[0]; }

  inline T* getCoef() { return &_coefs[0]; } 

  inline T* getCoef() const { return &_coefs[0]; } 

  inline int ptRow(int i) const { return _ptRows[i]; }

  inline int indCol(int i) const { return _indCols[i]; }

  inline int& indCol(int i) { return _indCols[i]; }

  inline const T Coef(int i) const { return _coefs[i]; }

  inline T& Coef(int i) { return _coefs[i]; }

  inline vector<int>& ptRows() { return _ptRows; }
  inline vector<int>& indCols() { return _indCols; }
  inline vector<T> & coefs() { return _coefs; }
  
  void prod(const T *u, T *v) const;
  void prodt(const T *u, T *v ) const;

  void normalize(const int type, const W* coefs0, Z* precDiag);

  void extractSquareMatrix(T *DSsingCofes, vector<int> &singVal);

 
  SparseMatrix<T>* PartialCopyCSR(vector<int> &permute, const int n,
				  bool transposed);

  int CSR_sym2unsym(int *ptRows, int *indCols, int *toSym, 
		    const int *ptSymRows, const int *indSymCols, 
		    const int dim, const bool upper_flag);
   
private:
  SparseMatrix& operator = ( const SparseMatrix& A );
  vector<int> _ptRows;  
  vector<int> _indCols; 
  vector<T> _coefs; 
  bool        _isSymmetric;  
  bool        _isUpper;
  bool        _isWhole;
}; // End class SparseMatrix

#endif
