/*! \file   ColumnMatrix.hpp
    \brief  Rectangular matrix view as a set of column vectors
    \author Xavier Juvigny, ONERA
    \date   Jan. 19th 2005
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

#ifndef _ALGEBRA_COLUMNMATRIX_HPP
# define _ALGEBRA_COLUMNMATRIX_HPP

# include "Algebra/PlainMatrix.hpp"

template<typename T>
class ColumnMatrix : public PlainMatrix<T>
{
public:

  using PlainMatrix<T>::coefs;
  using PlainMatrix<T>::addrCoefs;
  using PlainMatrix<T>::addrCoefs_pt;

  ColumnMatrix() : PlainMatrix<T>(), _nrows(0), _ncols(0)
  {}

  ColumnMatrix(int nrows, int ncols) :
    PlainMatrix<T>(), _nrows(0), _ncols(0)
  { init(nrows, ncols); }

  ColumnMatrix(int nrows, int ncols, T* coefs, bool isOwner) :
    PlainMatrix<T>(isOwner), _nrows(0), _ncols(0)
  { init(nrows, ncols, coefs, isOwner); }

  void init(int nrows, int ncols, bool later_allocation = false)
  {
    _nrows = nrows; 
    _ncols = ncols;
    if (!later_allocation) {
      PlainMatrix<T>::init((_nrows * _ncols));
    }
  }

  void init(int nrows, int ncols, T* cfs, bool isOwner)
  {
    _nrows = nrows; 
    _ncols = ncols;
    PlainMatrix<T>::init((_nrows * _ncols), cfs, isOwner);
  }

  void allocate()
  {
    PlainMatrix<T>::init((_nrows * _ncols));
  }

  ~ColumnMatrix() { }

  int nbColumns() const { return _ncols; }
  int nbRows() const { return _nrows; }
  int size() const { return _ncols * _nrows; }
  
  virtual T& operator () (int i, int j)
  {
# ifdef DISDEBUG
    assert(int(i) < nbRows());
    assert(int(j) < nbColumns());
#endif
    return coefs()[i + j * _nrows];
  }

  virtual const T& operator () (int i, int j) const
  {
# ifdef DISDEBUG
    assert(int(i) < nbRows());
    assert(int(j) < nbColumns());
#endif
    return coefs()[i + j * _nrows];
  }

  void ZeroClear()
  {
    PlainMatrix<T>::ZeroClear();
  }
 
  virtual ColumnMatrix<T>* clone() const
  {
    ColumnMatrix<T> *ret=new ColumnMatrix<T>;
    ret->copy(*this);
    return(ret);
  }
  /// \brief Deep copy of B
  void copy(const ColumnMatrix<T>& B)
  {
    _nrows=B._nrows;
    _ncols=B._ncols;
    PlainMatrix<T>::copy(B);
  }

  void free()
  {
    PlainMatrix<T>::free();
    _nrows = 0;
    _ncols = 0;
  }

private:
  ColumnMatrix(const ColumnMatrix& A);
  ColumnMatrix& operator = (const ColumnMatrix& A);
  int _nrows;
  int _ncols; 
};

#endif
