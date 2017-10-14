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
#include "Compiler/arithmetic.hpp"
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


template<typename T, typename U> void
normalize(const int type, const T* coefs0, SparseMatrix<T> *ptDA, U* u)
{
  const T zero(0.0);
  const T one(1.0);
  const U Uone(1.0);
  const U Uzero(0.0);
  const int n = ptDA->dimension();
  VectorArray<U> v(n);    //  Z* v = new Z[n];
  VectorArray<U> d(n);     //  Z* d = new Z[n];
  int *ptUnsymRows, *indUnsymCols, *indVals;
  if (ptDA->isSymmetric() && (type == KKT_SCALING)) {
    const int nnz = ptDA->nnz();
    if (ptDA->isWhole()) {
      ptUnsymRows = ptDA->getRows();
      indUnsymCols = ptDA->getIndCols();
      indVals = new int[nnz];
      for (int i = 0; i < nnz; i++) {
	indVals[i] = i;
      }
    }
    else {
      const int nnz0 = nnz * 2 - n;
      ptUnsymRows = new int[n + 1];
      indUnsymCols = new int[nnz0]; // diagonal entries exits though
      indVals = new int[nnz0];      // coeficient equals to zero
      // only used as extend symmetric symbolic structure to unsymmetric
      int nnz1;
      nnz1 = CSR_sym2unsym(ptUnsymRows, indUnsymCols, indVals, ptDA->getRows(), 
			   ptDA->getIndCols(), 
			   n, ptDA->isUpper());
      if (nnz1 != nnz0) {
	fprintf(stderr, 
		"%s %d : symmetric matrix has no diagonal entry %d != %d\n",
		__FILE__, __LINE__, nnz0, nnz1);
	exit(-1);
      }
    }
  }

  for (int i = 0; i < n; i++) {
    d[i] = Uzero;
    u[i] = Uzero;
    v[i] = Uzero;
  }
  for (int i = 0; i < n; i++) {
    for (unsigned k = (ptDA->getRows()[i] + 1); k < ptDA->getRows()[i + 1]; k++) {
      int j = ptDA->getIndCols()[k];
      if (i == j) {
	d[i] = blas_abs<T, U>(coefs0[k]);
      }
    }
  }
  // 23 Jun.2011 Atsushi
  // u[i] = |a(i, i)|, v[i] = max_{k}{|a(i,k)|,|a(k,j)|} 
  if (ptDA->isSymmetric() && (!ptDA->isWhole())) {
    for (int i = 0; i < n; i++) {
      u[i] = blas_abs<T, U>(coefs0[ptDA->getRows()[i]]);
      v[i] = SpMAX(v[i], u[i]);
      for (int k = (ptDA->getRows()[i] + 1); k < ptDA->getRows()[i + 1]; k++) {
	int j = ptDA->getIndCols()[k];
	U acoef = blas_abs<T, U>(coefs0[k]);
	v[i] = SpMAX(v[i], acoef);
	v[j] = SpMAX(v[j], acoef);
      }
    }
  }
  else {
    for (int i = 0; i < n; i++) {
      for (int k = ptDA->getRows()[i]; k < ptDA->getRows()[i + 1]; k++) {
	int j = ptDA->getIndCols()[k];
	if (j == i) {
	  u[i] = blas_abs<T, U>(coefs0[k]);
	  v[i] = SpMAX(v[i],u[i]);
	}
	else {
	  U acoef = blas_abs<T, U>(coefs0[k]);
	  v[i] = SpMAX(v[i], acoef);
	  v[j] = SpMAX(v[j], acoef);
	}
      }
    }
  }
  T alower, aupper, adiag;
  switch(type) {
  case DIAGONAL_SCALING:
    for (int i = 0; i < n; i++) {
      u[i] = ((u[i] != Uzero) ? sqrt<U>(Uone / u[i]) : 
	      ((v[i] != Uzero) ? sqrt<U>(Uone / v[i]) : Uone));
    }
    break;
  case KKT_SCALING:
    if (ptDA->isSymmetric() && (!ptDA->isWhole())) {
      for (int i = 0 ; i < n; i++) {
	if (d[i] != Uzero) {
	  u[i] = sqrt<U>(Uone / u[i]);
	}
	else {
	  T xtmp = zero;
	  for (int m = ptUnsymRows[i]; m < ptUnsymRows[i + 1]; m++) {
	    const int j = indUnsymCols[m];
	    alower = coefs0[indVals[m]];
	    int flag = 0;
	    for (int n = ptUnsymRows[j]; n < ptUnsymRows[j + 1]; n++) {
	      const int k = indUnsymCols[n];
	      if (k == j) {
		adiag = coefs0[indVals[n]];
		if (adiag != zero) {
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
	  u[i] = Uone / sqrt<U>(blas_abs<T, U>(xtmp));
	} // if (d[i] != zero)
      } // loop : i
    } // if (_isSymmetric && (!_isWhole))
    else {
      for (int i = 0 ; i < n; i++) {
	if (d[i] != Uzero) {
	  u[i] = sqrt<U>(Uone / u[i]);
	}
	else {
	  T xtmp = zero;
	  for (int m = ptDA->getRows()[i]; m < ptDA->getRows()[i + 1]; m++) {
	    const int j = ptDA->getIndCols()[m];
	    alower = coefs0[m];
	    int flag = 0;
	    for (int n = ptDA->getRows()[j]; n < ptDA->getRows()[j + 1]; n++) {
	      const int k = ptDA->getIndCols()[n];
	      if (k == j) {
		adiag = coefs0[n];
		if (adiag != zero) {
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
	  u[i] = Uone / sqrt<U>(blas_abs<T, U>(xtmp));
	} // if (d[i] != zero)
      }   // loop : i
    }  // if (_isSymmetric && (!_isWhole))
    break;
  default: 
    for (int i = 0; i < n; i++) {
      u[i] = Uone;
    }
    break;
  } // witch (type)
  // Scaling
  switch(type) {
  case NO_SCALING:
    for (int i = 0; i < ptDA->getRows()[n]; i++) {
      ptDA->getCoef()[i] = coefs0[i];
    }
    break;
  case DIAGONAL_SCALING:
  case KKT_SCALING:
    for ( int i = 0; i < n; i++) {
      for (int k = ptDA->getRows()[i]; k < ptDA->getRows()[i+1]; k++) {
	int j = ptDA->getIndCols()[k];
		ptDA->getCoef()[k] = coefs0[k] * u[i] *u[j];
      }
    }
    break;
  }
  
  //  delete [] v;
  //  delete [] d;
  if (ptDA->isSymmetric() && (type == KKT_SCALING)) {
    if (!ptDA->isWhole()) {
      delete [] ptUnsymRows;
      delete [] indUnsymCols;
    }
    delete [] indVals;
  }
}

template void
normalize<double, double>(const int type, const double *coefs0,
				  SparseMatrix<double> *ptDA, double *u);

template void
normalize<quadruple, quadruple>(const int type, 
					   const quadruple *coefs0,
					   SparseMatrix<quadruple> *ptDA,
					   quadruple *u);

template void
normalize<complex<double>, double>
   (const int type, const complex<double> *coefs0,
    SparseMatrix<complex<double> > *ptDA, double *u);

template void
normalize<complex<quadruple>, quadruple>
    (const int type, const complex<quadruple> *coefs0,
     SparseMatrix<complex<quadruple> > *ptDA, quadruple *u);

template void
normalize<float, float>(const int type, const float *coefs0,
				  SparseMatrix<float> *ptDA, float *u);

template void
normalize<complex<float>, float>
   (const int type, const complex<float> *coefs0,
    SparseMatrix<complex<float> > *ptDA, float *u);

//

template<typename T> void
SparseMatrix<T>::extractSquareMatrix(T* DSsingCoefs, vector<int> &singVal)
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
SparseMatrix<complex<double> >::extractSquareMatrix(complex<double>* DSsingCoefs, 
					   vector<int> &singVal);

template void
SparseMatrix<complex<quadruple> >::extractSquareMatrix(complex<quadruple>* DSsingCoefs, 
					     vector<int> &singVal);

template void
SparseMatrix<float>::extractSquareMatrix(float* DSsingCoefs, 
					  vector<int> &singVal);

template void
SparseMatrix<complex<float> >::extractSquareMatrix(complex<float>* DSsingCoefs, 
					   vector<int> &singVal);
//

template<typename T>
void SparseMatrix<T>::prod(const T *u, T *v) const
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
void SparseMatrix<complex<double> >::prod(const complex<double> *u, 
					  complex<double> *v) const;

template
void SparseMatrix<complex<quadruple> >::prod(const complex<quadruple> *u, 
					     complex<quadruple> *v) const;

template
void SparseMatrix<float>::prod(const float *u, float *v) const;

template
void SparseMatrix<complex<float> >::prod(const complex<float> *u, 
					  complex<float> *v) const;
//

template<typename T>
void SparseMatrix<T>::prodt(const T *u, T *v ) const
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
void SparseMatrix<complex<double> >::prodt(const complex<double> *u, 
				 complex<double> *v) const;

template
void SparseMatrix<complex<quadruple> >::prodt(const complex<quadruple> *u, 
				    complex<quadruple> *v) const;

template
void SparseMatrix<float>::prodt(const float *u, float *v) const;

template
void SparseMatrix<complex<float> >::prodt(const complex<float> *u, 
				 complex<float> *v) const;
//

template<typename T>
SparseMatrix<T>* SparseMatrix<T>::PartialCopyCSR(vector<int> &permute,
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
SparseMatrix<complex<double> >* SparseMatrix<complex<double> >::
PartialCopyCSR(vector<int> &permute, const int n,
	       bool transposed);

template
SparseMatrix<complex<quadruple> >* SparseMatrix<complex<quadruple> >::
PartialCopyCSR(vector<int> &permute, const int n,
	       bool transposed);

template
SparseMatrix<float>* SparseMatrix<float>::
PartialCopyCSR(vector<int> &permute, const int n,
	       bool transposed);

template
SparseMatrix<complex<float> >* SparseMatrix<complex<float> >::
PartialCopyCSR(vector<int> &permute, const int n,
	       bool transposed);

//
int CSR_sym2unsym(int *ptRows, int *indCols, int *toSym, 
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
