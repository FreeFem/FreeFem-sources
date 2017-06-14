/*! \file   SquareMatrix.hpp
    \brief  Definition of a template of square plain matrix
    \author X. Juvigny
    \date   January 19th 2005
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


#ifndef _ALGEBRA_SQUAREMATRIX_HPP
#define _ALGEBRA_SQUAREMATRIX_HPP
#include <cassert>
#include <cmath>
#include <vector>
#include "Algebra/PlainMatrix.hpp"
#include "Algebra/ColumnMatrix.hpp"

using std::vector;

/** \brief Square plain matrix

   This matrix is a square plain matrix with a global indices array
   to have correspondence with local numerotation.
*/
template<typename T>
class SquareMatrix : public PlainMatrix<T>
{
public:

  using PlainMatrix<T>::coefs;
  using PlainMatrix<T>::addrCoefs;
  using PlainMatrix<T>::addrCoefs_pt;

  SquareMatrix() : PlainMatrix<T>(), _rank(-1), _permute(),
                   _isFactorized(false), _kernelDetected(-1),
		   _singIdx(), 
		   _singIdx0()
  {
    _dim = 0;
    _num_blocks = 0;
    _isblocked = false;
  }

  SquareMatrix(int dim, int block_size) : 
    PlainMatrix<T>(dim * dim), _dim(dim), _rank(-1), 
    _permute(dim), _singIdx(), _singIdx0(),
    _kernelDetected(-1), _isFactorized(false)
  {
    _block_size = block_size;
    _num_blocks = (dim + _block_size - 1) / _block_size;
    if (_num_blocks > 0) {
      _nsing_block = new int[_num_blocks];
      for (int i = 0; i < _num_blocks; i++) {
	_nsing_block[i] = 0;
      }
    }
    _isblocked = true;
  } 

  SquareMatrix(int dim, int block_size, T* coefs,
	       int nsing, int *singidx0, int *singidx, 
	       int *permute,
	       bool isblocked,
	       bool isOwner) :
    PlainMatrix<T>(isOwner)
  { 
    _dim = dim;
    _block_size = block_size;
    _num_blocks = 0;
    _isblocked = isblocked;
    PlainMatrix<T>::init(dim * dim, coefs, isOwner);
    _singIdx.resize(nsing);
    for (int i = 0; i < nsing; i++) {
      _singIdx[i] = singidx[i];
    }
    _singIdx0.resize(nsing);
    for (int i = 0; i < nsing; i++) {
      _singIdx0[i] = singidx0[i];
    }
    _permute.resize(dim);
    for (int i = 0; i < dim; i++) {
      _permute[i] = permute[i];
    }
  }
  
  SquareMatrix(const SquareMatrix& A) : PlainMatrix<T>(A), 
					  _rank(A._rank), _permute(A._permute), 
					  _eps(A._eps), _singIdx(), _singIdx0(),
					  _dim(A._dim), 
					  _kernelDetected(-1),
					  _isFactorized(A._isFactorized),
					  _isblocked(A._isblocked)
  {}

  ~SquareMatrix() { 
    if (_num_blocks > 0) {
      delete [] _nsing_block;
    }
  }

  int dimension() const { return _dim; }

  virtual T& operator () (int i, int j)
  {
# ifdef DISDEBUG
    assert(int(i)<dimension());
    assert(int(j)<dimension());
#endif
    //  return coefs()[j+i*dimension()];
    //  05 Jul.2012 : Atsushi  
    return coefs()[i + j * _dim];
  }
  /** \brief Return \f$A_{ij}\f$ as constant reference

      \param i The ith row of the matrix
      \param j The jth column of the matrix
      \pre   0 <= i < dimension()
      \pre   0 <= j < dimension()
   */
  virtual const T& operator () (int i, int j) const
  {
# ifdef DISDEBUG
    assert(int(i)<dimension());
    assert(int(j)<dimension());
#endif
    //  return coefs()[j+i*dimension()];
    //  05 Jul.2012 : Atsushi  
    return coefs()[i + j * _dim];
  }
  /** \brief Tell if the squarematrix is factorized or not.

      \return true if factorized, else false
   */
  bool isFactorized() const { return _isFactorized; }
  /** \brief Return the rank of the factorized matrix

      \pre The squarematrix must be factorized or return -1
  */
  int rank() const
  {
# ifdef DISDEBUG
    assert(isFactorized());
# endif
    return _rank;
  }
  int* rank_val() { return &_rank; }
  /** \brief Return the last non zero pivot of the matrix

      \pre The squarematrix must be factorized else return 0
  */
  int dim_kern() const 
  {
    return _nsing;
  }
  void set_dim_kern(int nsing) 
  { _nsing = nsing; }

  int dim_kern_block(int k) const
  {
#ifdef DISDEBUG
    assert(k >= 0);
    assert(k < _num_blocks);
#endif
    return _nsing_block[k];
  }
  void set_dim_kern_block(int k, int nsing)
  {
#ifdef DISDEBUG
    assert(k >= 0);
    assert(k < _num_blocks);
#endif
    _nsing_block[k] = nsing;
  }
   
  void set_lastPivot(T lastpiv)
  {
    _lastpiv = lastpiv;
  }

  T lastPivot() const
  {
# ifdef DISDEBUG
    assert(isFactorized());
# endif
    return _lastpiv;
  }
  T *lastPivot_val() { return &_lastpiv; }

  /** \brief Return the pivot obtained after factorization

      \pre The squarematrix must be factorized
  */
  vector<int>& getPermute()
  {
#ifdef DISDEBUG
    assert(isFactorized());
#endif
    return _permute;
  }
  /** \brief Return the pivot obtained after factorization
      \pre The squarematrix must be factorized
  */
  const vector<int>& getPermute() const
  {
#ifdef DISDEBUG
    assert(isFactorized());
#endif
    return _permute;
  }
  /** \name Operations on the Squarematrix
   */
  //@{
  /** \brief Copy operator 

      \param A The squareMatrix to copy
      \return A reference on the current squareMatrix
  */
  SquareMatrix& operator = (const SquareMatrix& A)
  {
    if (this != &A)
    {
      PlainMatrix<T>::operator = (A); // _loc2glob=A._loc2glob;
      _rank = A._rank; 
      _permute = A._permute; 
      _eps = A._eps; 
      _lastpiv = A._lastpiv;
      _singIdx = A._singIdx; 
      _singIdx0 = A._singIdx0; 
      _isFactorized = A._isFactorized;
      _kernelDetected = A._kernelDetected;
    }
    return(*this);
  }
  /** \brief  Clone operator

      \return a new copy of the current squareMatrix
  */

  virtual SquareMatrix<T>* clone() const
  {
    //return new SquareMatrix<T>(*this);
    SquareMatrix<T> *ret = new SquareMatrix<T>;
    ret->copy(*this);
    return(ret);
  }

  /** \brief Deep copy of a square matrix in the current SquareMatrix

      \param A The matrix to copy
  */
  void copy (const SquareMatrix& A)
  {
    PlainMatrix<T>::copy(A); //_loc2glob=A._loc2glob;
    _rank = A._rank; 
    _permute = A._permute; 
    _eps = A._eps; 
    _lastpiv = A._lastpiv;
    _singIdx = A._singIdx;
    _singIdx0 = A._singIdx0;
    _isFactorized = A._isFactorized;
    _kernelDetected = A._kernelDetected;
  }

  void init(const int dim)
  {
    _dim = dim;
    PlainMatrix<T>::init(dim * dim);
    _rank = -1; 
    _permute.resize(dim);
    _singIdx.resize(0);
    _singIdx0.resize(0);
    _isFactorized = false;
  } 

  void free()
  {
    PlainMatrix<T>::free();
    _dim = 0;
    _permute.clear();
     _singIdx.clear();
    _singIdx0.clear();
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

  void setFactorizedState_pub() { _isFactorized = true; }
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
  void set_epspiv(T eps) { _eps = eps; }

  vector<int> &getSingIdx() {
    return _singIdx;
  }
 vector<int> &getSingIdx0() {
    return _singIdx0;
  }
  //  void set_lastpiv(int lastpiv) { _lastpiv = lastpiv; }
  
protected:
  // Internal methods
  /// \brief Set the matrix as factorized
  void setFactorizedState() { _isFactorized = true; }
  /// \brief Set dimension of the matrix
  void setDimension(int n) { _dim = n; }
  /// \brief Return the diagonal issued from the factorization
  // Attributs
  int  _rank;     //!< Rank of the matrix (initialized after factorization)
  int _nsing;
  int *_nsing_block;
  int _block_size;
  int _num_blocks;
  T    _lastpiv;  //!< Last non null pivot found during the factorization
  T    _eps;
private:
  // Private attributs
  int _dim;   //<! Dimension of the matrix
  vector<int> _permute; //<! 
  bool        _isFactorized; //<! Tell if the matrix is factorized or not
  int        _kernelDetected;  // 0 : without, 1 : with kernel (-1) : unknown
  //  IndiceArray _loc2glob; //!< Indice array describing local to global indices for columns
  // int _sing_max;
  vector<int> _singIdx;  // with permutation
  vector<int> _singIdx0; // without permutation 
  // 
  bool _isblocked;
}; // End class SquareMatrix

// # include "Algebra/SquareMatrix.tpp"
// ------------------------------------------------------------------------
/** \brief Output stream operator 

    Print the square matrix in human readable form.

    \param out The stream outpu
    \param A   The square matrix to print on the out stream
    \return The modified output stream operator
*/

template<typename T>
class SubSquareMatrix : public SquareMatrix<T>
{
public:
  SubSquareMatrix() : SquareMatrix<T>(), _loc2glob()
  { 
    _offdiag_2x2_status = false;
  }

  SubSquareMatrix(const vector<int>& l2g) : SquareMatrix<T>()
  { 
    init(l2g); 
  }

  void init(const vector<int>& l2g)
  { 
    int n = l2g.size();
    _loc2glob = l2g; // copy constructor
    SquareMatrix<T>::setDimension(n); 
    PlainMatrix<T>::init(n * n); 
    SquareMatrix<T>::getPermute().resize(n);
    _pivot_width.resize(n);
    _pivot_2x2.resize(0);
    _offdiag_2x2 = new T[n];
    _offdiag_2x2_status = true;
  }
  vector<int>& loc2glob()
  {
    return _loc2glob;
  }
  const vector<int>& loc2glob() const
  { 
    return _loc2glob; 
  }

  virtual ~SubSquareMatrix() {}

  void free()
  {
    SquareMatrix<T>::free();
    _loc2glob.clear();
    _pivot_width.clear();
    _pivot_2x2.clear();
    if (_offdiag_2x2_status) {
      delete [] _offdiag_2x2;
      _offdiag_2x2_status = false;
    }
  }

  vector<int>& getPivotWidth()
  {
    return _pivot_width;
  }
  /** \brief Return the pivot obtained after factorization
      \pre The squarematrix must be factorized
  */
  const vector<int>& getPivotWidth() const
  {
    return _pivot_width;
  }
  vector<int>& getPivot2x2()
  {
    return _pivot_2x2;
  }
  /** \brief Return the pivot obtained after factorization
      \pre The squarematrix must be factorized
  */
  const vector<int>& getPivot2x2() const
  {
    return _pivot_2x2;
  }

  T *addr2x2()
  {
    return _offdiag_2x2;
  }
private:
  vector<int> _loc2glob;
  vector<int> _pivot_width;
  vector<int> _pivot_2x2;
  T*          _offdiag_2x2;
  bool        _offdiag_2x2_status;
};
#endif
