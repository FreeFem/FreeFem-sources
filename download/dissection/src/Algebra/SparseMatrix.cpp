/*! \file   SparseMatrix.cpp
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

#include "Algebra/SparseMatrix.hpp"
#include "Algebra/VectorArray.hpp"
#include "Compiler/OptionCompiler.hpp"
#include "Driver/DissectionDefault.hpp" // for definition of scaling option
#include <cstring>
#include <cstdio>
#include <cstdlib>

inline
double SpMAX(double a, double b)
{
  return (a > b ? a : b);
}

inline
quadruple SpMAX(quadruple a, quadruple b)
{
  return (a > b ? a : b);
}


// T : higher precision to keep scaled matrix, "complex/real" value
// Z : precision to keep diagonal scaling matrix, "real" value
// W : precision to keep input entries, "complex/real" value
template<typename T, typename W, typename Z> void
SparseMatrix<T, W, Z>::normalize(const int type, const W* coefs0, Z* u)
{
  const Z zero(0.0);
  const Z one(1.0);
  const W Wzero(0.0);
  const int n = dimension();
  VectorArray<Z> v(n);    //  Z* v = new Z[n];
  VectorArray<Z> d(n);     //  Z* d = new Z[n];
  int *ptUnsymRows, *indUnsymCols, *indVals;
  if (_isSymmetric && (type == KKT_SCALING)) {
    if (_isWhole) {
      ptUnsymRows = &_ptRows[0];
      indUnsymCols = &_indCols[0];
      indVals = new int[nnz()];
      for (int i = 0; i < nnz(); i++) {
	indVals[i] = i;
      }
    }
    else {
      const int nnz0 = nnz() * 2 - n;
      ptUnsymRows = new int[n + 1];
      indUnsymCols = new int[nnz0]; // diagonal entries exits though
      indVals = new int[nnz0];      // coeficient equals to zero
      // only used as extend symmetric symbolic structure to unsymmetric
      int nnz1;
      nnz1 = CSR_sym2unsym(ptUnsymRows, indUnsymCols, indVals, getRows(), 
			   getIndCols(), 
			   n, isUpper());
      if (nnz1 != nnz0) {
	fprintf(stderr, 
		"%s %d : symmetric matrix has no diagonal entry %d != %d\n",
		__FILE__, __LINE__, nnz0, nnz1);
	exit(-1);
      }
    }
  }

  for (int i = 0; i < n; i++) {
    d[i] = zero;
    u[i] = zero;
    v[i] = zero;
  }
  for (int i = 0; i < n; i++) {
    for (unsigned k = (_ptRows[i] + 1); k < _ptRows[i + 1]; k++) {
      int j = _indCols[k];
      if (i == j) {
	d[i] = blas_abs<W, Z>(coefs0[k]);
      }
    }
  }
  // 23 Jun.2011 Atsushi
  // u[i] = |a(i, i)|, v[i] = max_{k}{|a(i,k)|,|a(k,j)|} 
  if (_isSymmetric && (!_isWhole)) {
    for (int i = 0; i < n; i++) {
      u[i] = blas_abs<W, Z>(coefs0[_ptRows[i]]);
      v[i] = SpMAX(v[i], u[i]);
      for (int k = (_ptRows[i] + 1); k < _ptRows[i + 1]; k++) {
	int j = _indCols[k];
	Z acoef = blas_abs<W, Z>(coefs0[k]);
	v[i] = SpMAX(v[i], acoef);
	v[j] = SpMAX(v[j], acoef);
      }
    }
  }
  else {
    for (int i = 0; i < n; i++) {
      for (int k = _ptRows[i]; k < _ptRows[i + 1]; k++) {
	int j = _indCols[k];
	if (j == i) {
	  u[i] = blas_abs<W, Z>(coefs0[k]);
	  v[i] = SpMAX(v[i],u[i]);
	}
	else {
	  Z acoef = blas_abs<W, Z>(coefs0[k]);
	  v[i] = SpMAX(v[i], acoef);
	  v[j] = SpMAX(v[j], acoef);
	}
      }
    }
  }
  W alower, aupper, adiag;
  switch(type) {
  case DIAGONAL_SCALING:
    for (int i = 0; i < n; i++) {
      u[i] = ((u[i] != zero) ? sqrt(one / u[i]) : 
	      ((v[i] != zero) ? sqrt(one / v[i]) : one));
    }
    break;
  case KKT_SCALING:
    if (_isSymmetric && (!_isWhole)) {
      for (int i = 0 ; i < n; i++) {
	if (d[i] != zero) {
	  u[i] = sqrt(one / u[i]);
	}
	else {
	  W xtmp = Wzero;
	  for (int m = ptUnsymRows[i]; m < ptUnsymRows[i + 1]; m++) {
	    const int j = indUnsymCols[m];
	    alower = coefs0[indVals[m]];
	    int flag = 0;
	    for (int n = ptUnsymRows[j]; n < ptUnsymRows[j + 1]; n++) {
	      const int k = indUnsymCols[n];
	      if (k == j) {
		adiag = coefs0[indVals[n]];
		if (adiag != Wzero) {
		  flag++;
		}
		else {
		  continue;
		}
	      }
	      if (k == i) {
		aupper = coefs0[indVals[n]];
		flag++;
	      }
	      if (flag == 2) {
		break;
	      }
	    } // loop : n
	    if (flag == 2) {
	      xtmp += alower * aupper / adiag;
	    }
	  } // loop : m
	  u[i] = one / sqrt(blas_abs<W, Z>(xtmp));
	} // if (d[i] != zero)
      } // loop : i
    } // if (_isSymmetric && (!_isWhole))
    else {
      for (int i = 0 ; i < n; i++) {
	if (d[i] != zero) {
	  u[i] = sqrt(one / u[i]);
	}
	else {
	  W xtmp = Wzero;
	  for (int m = _ptRows[i]; m < _ptRows[i + 1]; m++) {
	    const int j = _indCols[m];
	    alower = coefs0[m];
	    int flag = 0;
	    for (int n = _ptRows[j]; n < _ptRows[j + 1]; n++) {
	      const int k = _indCols[n];
	      if (k == j) {
		adiag = coefs0[n];
		if (adiag != Wzero) {
		  flag++;
		}
		else {
		  continue;
		}
	      }
	      if (k == i) {
		aupper = coefs0[n];
		flag++;
	      }
	      if (flag == 2) {
		break;
	      }
	    } // loop : n
	    if (flag == 2) {
	      xtmp += alower * aupper / adiag;
	    }
	  } // loop : m
	  u[i] = one / sqrt(blas_abs<W, Z>(xtmp));
	} // if (d[i] != zero)
      }   // loop : i
    }  // if (_isSymmetric && (!_isWhole))
    break;
  default: 
    for (int i = 0; i < n; i++) {
      u[i] = one;
    }
    break;
  } // witch (type)
  // Scaling
  if (type == NO_SCALING) {
    for (int i = 0; i < _ptRows[n]; i++) {
      _coefs[i] = coefs0[i];
    }
    // blas_copy<T>(_ptRows[n], coefs0, 1, &_coefs[0], 1);
    // const W* coefs0 to T* _coefs is not supported becuase of type conversion
  }
  else {
    for ( int i = 0; i < n; i++) {
      for (int k = _ptRows[i]; k < _ptRows[i+1]; k++) {
	int j = _indCols[k];
	//      _coefs[k] = tohigher((1), coefs0[k] * W(u[i]) * W(u[j]));
	_coefs[k] = coefs0[k] * W(u[i]) * W(u[j]);
      }
    }
  }
  //  delete [] v;
  //  delete [] d;
  if (_isSymmetric && (type == KKT_SCALING)) {
    if (!_isWhole) {
      delete [] ptUnsymRows;
      delete [] indUnsymCols;
    }
    delete [] indVals;
  }
}

template void
SparseMatrix<double>::normalize(const int type, 
				const double *coefs0,
				double *u);

template void
SparseMatrix<quadruple>::normalize(const int type, 
				   const quadruple *coefs0,
				   quadruple *u);
template void
SparseMatrix<quadruple, double>::normalize(const int type, 
					   const double *coefs0,
					   double *u);


template void
SparseMatrix<complex<double>, complex<double>,
	     double>::normalize(const int type, 
				const complex<double> *coefs0,
				double *u);

template void
SparseMatrix<complex<quadruple>, complex<quadruple>,
	     quadruple>::normalize(const int type, 
				   const complex<quadruple> *coefs0,
				   quadruple *u);

template void
SparseMatrix<complex<quadruple>, complex<double>,
	     double>::normalize(const int type, 
				const complex<double> *coefs0,
				double *u);

//

template<typename T, typename W, typename Z> void
SparseMatrix<T, W, Z>::extractSquareMatrix(T* DSsingCoefs, vector<int> &singVal)
{
  const T zero(0.0);
  const int nsing = singVal.size();
  if (_isSymmetric && (!_isWhole)) {
    if (_isUpper) {
      for (int i = 0; i < nsing; i++) {
	for (int j = i; j < nsing; j++) {
	  DSsingCoefs[i + j * nsing] = 0.0;
	  const int icol = singVal[i];
	  int itmp = icol;
	  for (int it = ptRow(icol); it < ptRow(icol + 1); it++) {
	    itmp = indCol(it);
	    if (itmp == singVal[j]) {
	      DSsingCoefs[i + j * nsing] = Coef(it);
	    }
	  }
	}
	//symmetrize
	for (int j = i + 1; j < nsing; j++) {
	  DSsingCoefs[j + i * nsing] = DSsingCoefs[i + j * nsing];
	}
      }
    }  //  if (isUpper()) 
    else {
      for (int i = 0; i < nsing; i++) {
	for (int j = 0; j <= i; j++) {
	  DSsingCoefs[i + j * nsing] = zero;
	  const int icol = singVal[i];
	  int itmp = icol;
	  for (int it = ptRow(icol); it < ptRow(icol + 1); it++) {
	    itmp = indCol(it);
	    if (itmp == singVal[j]) {
	      DSsingCoefs[i + j * nsing] = Coef(it);
	    }
	  }
	}
	//symmetrize
	for (int j = 0; j < i; j++) {
	  DSsingCoefs[j + i * nsing] = DSsingCoefs[i + j * nsing];
	}
      }
    } //  if (isUpper()) 
  }
  else {
    // 
    for (int i = 0; i < nsing; i++) {
      for (int j = 0; j < nsing; j++) {
	DSsingCoefs[i + j * nsing] = zero;
	const int icol = singVal[i];
	for (int it = ptRow(icol); it != ptRow(icol + 1); it++) {
	  if (indCol(it) == singVal[j]) {
	    DSsingCoefs[i + j * nsing] = Coef(it);
	    break;
	  }
	}
      }
    }
  }
}

template void
SparseMatrix<double>::extractSquareMatrix(double* DSsingCoefs, 
					  vector<int> &singVal);
template void
SparseMatrix<quadruple>::extractSquareMatrix(quadruple* DSsingCoefs, 
					     vector<int> &singVal);

template void
SparseMatrix<quadruple, double>::extractSquareMatrix(quadruple* DSsingCoefs, 
					     vector<int> &singVal);
template void
SparseMatrix<complex<double>, complex<double>,
	     double >::extractSquareMatrix(complex<double>* DSsingCoefs, 
					   vector<int> &singVal);

template void
SparseMatrix<complex<quadruple>, complex<quadruple>,
	     quadruple>::extractSquareMatrix(complex<quadruple>* DSsingCoefs, 
					     vector<int> &singVal);

template void
SparseMatrix<complex<quadruple>, complex<double>,
	     double>::extractSquareMatrix(complex<quadruple>* DSsingCoefs, 
					     vector<int> &singVal);
//

template<typename T, typename W, typename Z>
void SparseMatrix<T, W, Z>::prod(const T *u, T *v) const
{
  const T zero(0.0);
  for (int i = 0; i < dimension(); i++) {
    v[i] = zero;
  }
  if (_isSymmetric && (!_isWhole)) {
    for(int i = 0; i < dimension(); i++) {
      // assumption data structure
      // _indCols[k] < _indCols[k + 1]
      // _indCols[_ptRows[i]] = i          for _isUpper == true
      // _indCols[_ptRows[i + 1] - 1] = i  for _isUpper == false
      int ibegin, iend;
      if (_isUpper) {
	ibegin = _ptRows[i] + 1;
	iend = _ptRows[i + 1];
      }
      else {
	ibegin = _ptRows[i];
	iend = _ptRows[i + 1] - 1;
      }
      for (int k = ibegin; k < iend; k++) {
	const int icol = _indCols[k];
	v[i] += _coefs[k] * u[icol];
	v[icol] += _coefs[k] * u[i];
      }
      const int k = _isUpper ? _ptRows[i] : (_ptRows[i + 1] - 1); // diagonal
      v[i] += _coefs[k] * u[_indCols[k]];
    } // loop : i
  }
  else {
    for(int i = 0; i < dimension(); i++) {
      for (int k = _ptRows[i]; k < _ptRows[i + 1]; k++) {
	v[i] += _coefs[k] * u[_indCols[k]];
      }
    }
  }
}

template
void SparseMatrix<double>::prod(const double *u, double *v) const;

template
void SparseMatrix<quadruple>::prod(const quadruple *u, quadruple *v) const;

template
void SparseMatrix<quadruple, double>::prod(const quadruple *u,
					   quadruple *v) const;

template
void SparseMatrix<complex<double>, complex<double>,
		  double>::prod(const complex<double> *u, 
				complex<double> *v) const;

template
void SparseMatrix<complex<quadruple>, complex<quadruple>,
		  quadruple>::prod(const complex<quadruple> *u, 
				   complex<quadruple> *v) const;
template
void SparseMatrix<complex<quadruple>, complex<double>,
		  double>::prod(const complex<quadruple> *u, 
				   complex<quadruple> *v) const;
//

template<typename T, typename W, typename Z>
void SparseMatrix<T, W, Z>::prodt(const T *u, T *v ) const
{
  const T zero(0.0);
  for (int i = 0; i < dimension(); i++) {
    v[i] = zero;
  }
  for(int i = 0; i < dimension(); i++) {
    for (int k = _ptRows[i]; k < _ptRows[i + 1]; k++) {
      v[_indCols[k]] += _coefs[k]*u[i];
    }
  }
}

template
void SparseMatrix<double>::prodt(const double *u, double *v) const;

template
void SparseMatrix<quadruple>::prodt(const quadruple *u, quadruple *v) const;

template
void SparseMatrix<quadruple, double>::prodt(const quadruple *u,
					    quadruple *v) const;

template
void SparseMatrix<complex<double>, complex<double>,
		  double>::prodt(const complex<double> *u, 
				 complex<double> *v) const;

template
void SparseMatrix<complex<quadruple>, complex<quadruple>,
		  quadruple>::prodt(const complex<quadruple> *u, 
				    complex<quadruple> *v) const;
template
void SparseMatrix<complex<quadruple>, complex<double>,
		  double>::prodt(const complex<quadruple> *u, 
				    complex<quadruple> *v) const;
//

template<typename T, typename W, typename Z>
SparseMatrix<T>* SparseMatrix<T, W, Z>::PartialCopyCSR(vector<int> &permute,
							const int n,
							bool transposed)
{  
  vector<int> dist_ptRows;
  dist_ptRows.resize(n + 1);
  int *source_ptRows = getRows(); //&source->beginRows()[0];
  SparseMatrix<T> *dist;
  if (_isSymmetric && !(_isWhole)) {
    int *width = new int[n];
    for (int i = 0; i < n; i++) {
      width[i] = source_ptRows[permute[i] + 1] - source_ptRows[permute[i]];
    }
    for (int i = 0; i < dimension(); i++) {
      for (int it = ptRow(i); it < ptRow(i + 1); it++) {
	const int icol = indCol(it);
	// does not count diagonal
	if (icol != i) {
	  for (int k = 0; k < n; k++) {
	    if (permute[k] == icol) {
	      width[k]++;
	      break;
	    }
	  }
	} // if (icol != i)
      }
    } // loop : i
    dist_ptRows[0] = 0;
    for (int i = 0; i < n; i++) {
      dist_ptRows[i + 1] = dist_ptRows[i] + width[i];
    }
    delete [] width;
    dist = new SparseMatrix<T>(n, dist_ptRows[n], dist_ptRows, 
			       false, false, false);
    if (isUpper()) {
      for (int k = 0; k < n; k++) {
	const int irow = permute[k];
	int jt = dist->ptRow(k);
	for (int i = 0; i < irow; i++) {
	  for (int it = ptRow(i); it < ptRow(i + 1); it++) {
	    const int icol = indCol(it);
	    // strictly lower == not touch the diagonal
	    if (icol != i) {
	      if (irow == icol) {
		dist->indCol(jt) = i;
		dist->Coef(jt) = Coef(it);
		jt++;
		break;
	      }
	    }
	  }
	} // loop : i
	// copy upper to upper
	for (int it = ptRow(irow); it < ptRow(irow + 1); it++) {
	  dist->indCol(jt) = indCol(it);
	  dist->Coef(jt) = Coef(it);
	  jt++;
	}
      }  // loop : k
    }       // order of filling nonzero elements is different from the upper 
    else {  // case to keep increaing indcols[]
           
      for (int k = 0; k < n; k++) {
	const int irow = permute[k];
	int jt = dist->ptRow(k);
	// copy lower to lower
	for (int it = ptRow(irow); it < ptRow(irow + 1); it++) {
	  dist->indCol(jt) = indCol(it);
	  dist->Coef(jt) = Coef(it);
	  jt++;
	}
	for (int i = irow; i < dimension(); i++) {
	  for (int it = ptRow(i); it < ptRow(i + 1); it++) {
	    const int icol = indCol(it);
	    // strictly lower == not touch the diagonal
	    if (icol != i) {
	      if (irow == icol) {
		dist->indCol(jt) = i;
		dist->Coef(jt) = Coef(it);
		jt++;
		break;
	      }
	    }
	  }
	} // loop : i
      }  // loop : k
    }
  } //   if (_isSymmetric)
  else {
    dist_ptRows[0] = 0;
    // non-zero pattern is symmetric
    for (int i = 0; i < n; i++) {
      const int m = permute[i];
      dist_ptRows[i + 1] = (dist_ptRows[i] + 
			    (source_ptRows[m + 1] - source_ptRows[m]));
    }
				// isOwner = true, isSym = false
    dist = new SparseMatrix<T>(n, dist_ptRows[n], dist_ptRows, 
			       false, false, false);
    if (transposed) {
      for (int k = 0; k < n; k++) {
	const int irow = permute[k];
	int jt = dist->ptRow(k);
	for (int i = 0; i < dimension(); i++) {
	  for (int it = ptRow(i); it < ptRow(i + 1); it++){
	    const int icol = indCol(it);
	    if (irow == icol) {
	      dist->indCol(jt) = i;
	      dist->Coef(jt) = Coef(it);
	      jt++;
	      break;
	    }
	  }
	}  // loop : i
      }
    }
    else {
      for (int k = 0; k < n; k++) {
	const int irow = permute[k];
	int it = ptRow(irow);
	int jt = dist->ptRow(k);
	// copy lower to lower
	for (; it < ptRow(irow + 1); it++, jt++) {;
	  dist->indCol(jt) = indCol(it);
	  dist->Coef(jt) = Coef(it);
	}
      }
    }
  } //   if (_isSymmetric)
  return dist;
}

template
SparseMatrix<double>* SparseMatrix<double>::
PartialCopyCSR(vector<int> &permute, const int n,
	       bool transposed);

template
SparseMatrix<quadruple>* SparseMatrix<quadruple>::
PartialCopyCSR(vector<int> &permute, const int n,
	       bool transposed);

template
SparseMatrix<quadruple>* SparseMatrix<quadruple, double>::
PartialCopyCSR(vector<int> &permute, const int n,
	       bool transposed);

template
SparseMatrix<complex<double> >*
SparseMatrix<complex<double>, complex<double>, double>::
PartialCopyCSR(vector<int> &permute, const int n,
	       bool transposed);

template
SparseMatrix<complex<quadruple> >*
SparseMatrix<complex<quadruple>, complex<quadruple>, quadruple>::
PartialCopyCSR(vector<int> &permute, const int n,
	       bool transposed);

template
SparseMatrix<complex<quadruple> >*
SparseMatrix<complex<quadruple>, complex<double>, double>::
PartialCopyCSR(vector<int> &permute, const int n,
	       bool transposed);
//

template<typename T, typename W, typename Z> int 
SparseMatrix<T, W, Z>::CSR_sym2unsym(int *ptRows, int *indCols, int *toSym, 
				     const int *ptSymRows,
				     const int *indSymCols, 
				     const int dim, const bool upper_flag)
{
  int* nbIndPerRow = new int[dim];

  memset(nbIndPerRow, 0, dim*sizeof(int));

  for (int i = 0; i < dim; i++) {
    nbIndPerRow[i] += ptSymRows[i + 1] - ptSymRows[i];
    int ibegin = ptSymRows[i] + (upper_flag ? 1 : 0);
    int iend = ptSymRows[i + 1] + (upper_flag ? 0 : (-1));
    for (int k = ibegin; k < iend; k++) {
      nbIndPerRow[indSymCols[k]]++;
    }
  }
  // Build ptRows array :
  // ...................
  ptRows[0] = 0;
  for (int i = 0; i < dim; i++) {
    ptRows[i + 1] = ptRows[i] + nbIndPerRow[i];
  }
  //  CHECK(ptRows[dim] == (2 * nz - dim),
  //	"error in sym2unsym() : Wrong number of non zeros elemnts in ptRows !");
	    // Allocate and fill indices columns :
  memset(nbIndPerRow, 0, (dim * sizeof(int)));
  // for upper case, nbIndPerRow[i] keeps entries added by transposed operation
  // but for lower case counts all nonzero entries in progress 
  for (int i = 0; i < dim; i++) {
    int itmp = ptRows[i] + nbIndPerRow[i];
    for (int k = ptSymRows[i]; k < ptSymRows[i + 1]; k++) {
      indCols[itmp] = indSymCols[k];
      toSym[itmp] = k;
      itmp++;
    } // loop : k
    if (!upper_flag) {
      nbIndPerRow[i] = itmp - ptRows[i];
    }
    //    memcpy(indCols + (ptRows[i] + nbIndPerRow[i]), 
    //	   indSymCols + ptSymRows[i],
    //	   (ptSymRows[i + 1] - ptSymRows[i]) * sizeof(int));
    int ibegin = ptSymRows[i] + (upper_flag ? 1 : 0);
    int iend = ptSymRows[i + 1] + (upper_flag ? 0 : (-1));
    for (int k = ibegin; k < iend; k++) {
      const int j = indSymCols[k];
      const int jtmp = ptRows[j] + nbIndPerRow[j];
      indCols[jtmp] = i;
      toSym[jtmp] = k;
      nbIndPerRow[j]++;
    } // loop : k
  }   // loop : i
  // extract other toSym[]
  delete [] nbIndPerRow;
  return ptRows[dim];
}
