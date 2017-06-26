/*! \file CopyMatrix.cpp
    \brief task mangemanet of dissection algorithm
    \author Atsushi Suzuki, Laboratoire Jacques-Louis Lions
    \date   Aug. 09th 2015
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

#include "Driver/CopyMatrix.hpp"

void CopySparseMatrix(SparseMatrix<double> *b,
		      SparseMatrix<quadruple, double, double> *a)
{
  b->ptRows() = a->ptRows();
  b->indCols() = a->indCols();
  const int nnz = a->nnz();
  b->coefs().resize(nnz);
  for (int i = 0; i < nnz; i++) {
    b->coefs()[i] = tolower<quadruple, double>(a->coefs()[i]);
  }
}

void CopySparseMatrix(SparseMatrix<complex<double>, complex<double>, double> *b,
		      SparseMatrix<complex<quadruple>, complex<double>, double> *a)
{
  b->ptRows() = a->ptRows();
  b->indCols() = a->indCols();
  const int nnz = a->nnz();
  b->coefs().resize(nnz);
  for (int i = 0; i < nnz; i++) {
    b->coefs()[i] = tolower<complex<quadruple>, complex<double> >(a->coefs()[i]);
  }
}
		      
template<typename T, typename W>
void CopySquareBlockMatrix(SquareBlockMatrix<W> &b,
			   SquareBlockMatrix<T> &a)
{
  // b.init(a.dimension(), a.block_size(), a.isSym(), 0);   // initialize b
  b.allocate();
  b.getPermute() = a.getPermute();          // copy vector<int> _permute
  b.getNsingBlock() = a.getNsingBlock();    // copy vector<int> _nsing_block
  b.getSingIdx() = a.getSingIdx();          // copy vector<int> _singIdx
  b.getSingIdx0() = a.getSingIdx0();        // copy vector<int> _singIdx0
  b.set_KernelDetected(a.KernelDetected()); // copy bool _kernelDetected
  b.set_rank(a.rank());                     // copy int _rank, int _nsing
  b.set_lastPivot(a.lastPivot());           // copy double _lastpiv
  if (a.isBlocked()) {                      // copy bool _isblocked
    b.setBlocked();
  }
  else {
    b.unsetBlocked();
  }
  if (a.isSym()) {
    const int num_blocks = a.num_blocks();
    for (int i = 0 ; i < num_blocks; i++) {
      for (int j = i ; j < num_blocks; j++) {
	const int nrow = a.nrowBlock(i);
	const int ncol = a.nrowBlock(j);
	for (int k = 0; k < (nrow * ncol); k++) {
	  b.addrCoefBlock(i, j)[k] = tolower<T, W>(a.addrCoefBlock(i, j)[k]);
	}
      } // loop : j
    }   // loop : i
  }
  else {
    const int num_blocks = a.num_blocks();
    for (int i = 0 ; i < num_blocks; i++) {
      for (int j = 0 ; j < num_blocks; j++) {
	const int nrow = a.nrowBlock(i);
	const int ncol = a.nrowBlock(j);
	for (int k = 0; k < (nrow * ncol); k++) {
	  b.addrCoefBlock(i, j)[k] = tolower<T, W>(a.addrCoefBlock(i, j)[k]);
	}
      } // loop : j
    }   // loop : i
  }
}

template
void CopySquareBlockMatrix<quadruple, double>(SquareBlockMatrix<double> &b,
					      SquareBlockMatrix<quadruple> &a);

template
void CopySquareBlockMatrix<complex<quadruple>,
			   complex<double> >(SquareBlockMatrix<complex<double> > &b,
			   SquareBlockMatrix<complex<quadruple> > &a);

//

template<typename T, typename W>
void CopyRectBlockMatrix(RectBlockMatrix<W> &b,
			 RectBlockMatrix<T> &a)
{
  for (int i = 0 ; i < a.num_blocks_r(); i++) {
    for (int j = 0 ; j < a.num_blocks_c(); j++) {
      const int nrow = a.nrowBlock(i);
      const int ncol = a.ncolBlock(j);
      for (int k = 0; k < (nrow * ncol); k++) {
	b.addrCoefBlock(i, j)[k] = tolower<T, W>(a.addrCoefBlock(i, j)[k]);
      }
    } // loop : j
  }   // loop : i
}

template
void CopyRectBlockMatrix<quadruple, double>(RectBlockMatrix<double> &b,
					    RectBlockMatrix<quadruple> &a);

template
void CopyRectBlockMatrix<complex<quadruple>, complex<double> >(RectBlockMatrix<complex<double> > &b,
			 RectBlockMatrix<complex<quadruple> > &a);
//

template<typename T, typename U, typename W, typename Z>
void CopyTridiagBlockMatrix(TridiagBlockMatrix<W, Z> &b,
			    TridiagBlockMatrix<T, U> &a,
			    W *coef)
{
  b.init(a.dimension(), a.block_size(), a.isSym());

  b.setNfront(a.Nfront());
  b.setMaxdim(a.maxdim());
  b.setNnz(a.nnz());
  b.setNop(a.nop());
  b.setDiag_block_alloc_status(a.diag_block_alloc_status());
  b.setNscol(a.nscol());
  b.setNsing(a.nsing());
  b.setDetected(a.detected());

  b.getPtRows() = a.getPtRows();
  b.getIndCols() = a.getIndCols();
  b.getIndVals() = a.getIndVals();
  b.getNew2old() = a.getNew2old();
  b.getP_front() = a.getP_front();
  b.getPermute() = a.getPermute();
  b.getPermute_ginv() = a.getPermute_ginv();
  b.getP_diag() = a.getP_diag();
  b.getP_upper() =  a.getP_upper();
  b.getList_schur() = a.getList_schur();
  b.getList_elim() = a.getList_elim();
  b.getNum_null() = a.getNum_null();
  b.setCoef(coef);

  const int nfront = a.Nfront();
  b.getaddrDiagMatrix() = new ColumnMatrix<W>[nfront];
#ifndef SPARSE_OFFDIAG
  b.getaddrLowerMatrix() = new ColumnMatrix<W>[nfront];
  b.getaddrUpperMatrix() = new ColumnMatrix<W>[nfront];
#endif
  for (int n = 0; n < nfront; n++) {
    {
      const int nbCols = a.getaddrDiagMatrix()[n].nbColumns();
      const int nbRows = a.getaddrDiagMatrix()[n].nbRows();
      const int size = a.getaddrDiagMatrix()[n].size();
      b.getaddrDiagMatrix()[n].init(nbRows, nbCols);  // allocation
      for (int i = 0; i < size; i++) {
	b.getaddrDiagMatrix()[n].coefs()[i] =
	  tolower<T, W>(a.getaddrDiagMatrix()[n].coefs()[i]);
      }
    }
#ifndef SPARSE_OFFDIAG
    {
      const int nbCols = a.getaddrLowerMatrix()[n].nbColumns();
      const int nbRows = a.getaddrLowerMatrix()[n].nbRows();
      const int size = a.getaddrLowerMatrix()[n].size();
      b.getaddrLowerMatrix()[n].init(nbRows, nbCols);  // allocation
      b.getaddrUpperMatrix()[n].init(nbRows, nbCols);  // allocation
      for (int i = 0; i < size; i++) {
	b.getaddrLowerMatrix()[n].coefs()[i] =
	  tolower<T, W>(a.getaddrLowerMatrix()[n].coefs()[i]);
	b.getaddrUpperMatrix()[n].coefs()[i] =
	  tolower<T, W>(a.getaddrUpperMatrix()[n].coefs()[i]);
      }
    }
    #endif
  }
  b.getA12().init(a.getA12().nbRows(), a.getA12().nbColumns()); // allocation
  for (int i = 0; i < a.getA12().size(); i++) {
    b.getA12().coefs()[i] = tolower<T, W>(a.getA12().coefs()[i]);
  }
  b.getA21().init(a.getA21().nbRows(), a.getA21().nbColumns()); // allocation
  for (int i = 0; i < a.getA21().size(); i++) {
    b.getA21().coefs()[i] = tolower<T, W>(a.getA21().coefs()[i]);
  }
  b.getS22().init(a.getS22().nbRows(), a.getS22().nbColumns()); // allocation
  for (int i = 0; i < a.getS22().size(); i++) {
    b.getS22().coefs()[i] = tolower<T, W>(a.getS22().coefs()[i]);
  }
}

template
void CopyTridiagBlockMatrix(TridiagBlockMatrix<double, double> &b,
			    TridiagBlockMatrix<quadruple, quadruple> &a,
			    double *coef);
template
void CopyTridiagBlockMatrix(TridiagBlockMatrix<complex<double>, double> &b,
			    TridiagBlockMatrix<complex<quadruple>,
			                      quadruple> &a,
			    complex<double> *coef);
//

template<typename T, typename U, typename W, typename Z>
void CopyDissectionMatrix(DissectionMatrix<W, Z> *b,
			  DissectionMatrix<T, U> *a,
			  SquareBlockMatrix<W> *diag,           // pointers
			  RectBlockMatrix<W> *lower,
			  RectBlockMatrix<W> *upper)
{
  b->setNb(a->nb());
  b->setLevel(a->level());
  b->setNrow(a->nrow());
  b->setNcol_offdiag(a->ncol_offdiag());
  b->setIsSym(a->isSym());
  b->setIslast(a->islast());
  b->setAlignedFather(a->alignedFather());
  // copy pointer to matrices that have lower accuracy <W, Z>
  b->setColorTridiagBlockMatrix(a->ColorTridiagBlockMatrix());
  b->paddrdiagBlock() = diag;
  //  b->paddrtridiagBlock() = tridiag;
  b->paddrupperBlock() = lower;
  b->paddrupperBlock() = upper;
}

template
void CopyDissectionMatrix<quadruple, quadruple,
			  double, double>(DissectionMatrix<double, double> *b,
			  DissectionMatrix<quadruple, quadruple> *a,
			  SquareBlockMatrix<double> *diag,
			  RectBlockMatrix<double> *lower,
			  RectBlockMatrix<double> *upper);

template
void CopyDissectionMatrix<complex<quadruple>, quadruple,
			  complex<double>, double>(DissectionMatrix<complex<double>, double> *b,
			  DissectionMatrix<complex<quadruple>, quadruple> *a,
			  SquareBlockMatrix<complex<double> > *diag,
			  RectBlockMatrix<complex<double> > *lower,
			  RectBlockMatrix<complex<double> > *upper);
//

template<typename T, typename W>
void CopySchurMatrix(SchurMatrix<W> &b,
		     SchurMatrix<T> &a)
{
  b.getSlduList() = a.getSlduList();
  {
    const int dim = a.getSldu().dimension();
    b.getSldu().init(a.getSldu().loc2glob());
    for (int i = 0; i < (dim * dim); i++) {
      b.getSldu().addrCoefs()[i] = tolower<T, W>(a.getSldu().addrCoefs()[i]);
    }
    for (int i = 0; i < dim; i++) {
      b.getSldu().addr2x2()[i] = tolower<T, W>(a.getSldu().addr2x2()[i]);
    }
    b.getSldu().getPivotWidth() = a.getSldu().getPivotWidth();
    b.getSldu().getPivot2x2() = a.getSldu().getPivot2x2();
    b.getSldu().getPermute() = a.getSldu().getPermute();
  }
#if 0
  {
    const int dim = a.getSchur().dimension();
    for (int i = 0; i < (dim * dim); i++) {
      b.getSchur().addrCoefs()[i] = tolower<T, W>(a.getSchur().addrCoefs()[i]);
    }
    b.getSchur().getPermute() = a.getSchur().getPermute();
  }
#endif
  {
    const bool isUpper = a.getArow()->isUpper();
    const bool isSym = a.getArow()->isSymmetric();
    const int dim = a.getArow()->dimension();
    const int nnz = a.getArow()->nnz();
    b.getArow() = new SparseMatrix<W>(dim, nnz,
				      &(a.getArow()->getRows()[0]),
				      &(a.getArow()->getIndCols()[0]),
				      isSym, isUpper);
    for (int i = 0; i < nnz; i++) {
      b.getArow()->Coef(i) = tolower<T, W>(a.getArow()->Coef(i));
    }
  }
  {
    const bool isUpper = a.getAcol()->isUpper();
    const bool isSym = a.getAcol()->isSymmetric();
    const int dim = a.getAcol()->dimension();
    const int nnz = a.getAcol()->nnz();
    b.getAcol() = new SparseMatrix<W>(dim, nnz,
				      &(a.getAcol()->getRows()[0]),
				      &(a.getAcol()->getIndCols()[0]),
				      isSym, isUpper);
    for (int i = 0; i < nnz; i++) {
      b.getAcol()->Coef(i) = tolower<T, W>(a.getAcol()->Coef(i));
    }
  }
  {
    const int size = a.getScol().size();
    for (int i = 0; i < size; i++) {
      b.getScol().addrCoefs()[i] = tolower<T, W>(a.getScol().addrCoefs()[i]);
    }
  }
}

template
void CopySchurMatrix(SchurMatrix<double> &b,
		     SchurMatrix<quadruple> &a);

template
void CopySchurMatrix(SchurMatrix<complex<double> > &b,
		     SchurMatrix<complex<quadruple> > &a);

//

template<typename T, typename W>
void CopyKernelMatrix(KernelMatrix<W> &b,
		      KernelMatrix<T> &a)
{
  b.set_dimension(a.dimension());
  b.getSingIdx() = a.getSingIdx();
  b.getKernListEq() = a.getKernListEq();
  {
    b.getKernBasis().init(a.getKernBasis().nbRows(),
			  a.getKernBasis().nbColumns());
    const int size = a.getKernBasis().size();
    for (int i = 0; i < size; i++) {
      b.getKernBasis().addrCoefs()[i] =
	tolower<T, W>(a.getKernBasis().addrCoefs()[i]);
    }
  }
  {
    b.getTKernBasis().init(a.getKernBasis().nbRows(),
			  a.getKernBasis().nbColumns());
	
    const int size = a.getTKernBasis().size();
    for (int i = 0; i < size; i++) {
      b.getTKernBasis().addrCoefs()[i] =
	tolower<T, W>(a.getTKernBasis().addrCoefs()[i]);
    }
  }
  {
    const int dim = a.getKernProj().dimension();
    b.getKernProj().init(dim);
    for (int i = 0; i < (dim * dim); i++) {
      b.getKernProj().addrCoefs()[i] = tolower<T, W>(a.getKernProj().addrCoefs()[i]);
    }
    b.getKernProj().getPermute() = a.getKernProj().getPermute();
  }
  {

    const int dim = a.getTKernProj().dimension();
    b.getTKernProj().init(dim);
    for (int i = 0; i < (dim * dim); i++) {
      b.getTKernProj().addrCoefs()[i] =
	tolower<T, W>(a.getTKernProj().addrCoefs()[i]);
    }
    b.getTKernProj().getPermute() = a.getTKernProj().getPermute();
  }
  {
    const int dim = a.getNTKernProj().dimension();
    b.getNTKernProj().init(dim);
    for (int i = 0; i < (dim * dim); i++) {
      b.getNTKernProj().addrCoefs()[i] =
	tolower<T, W>(a.getNTKernProj().addrCoefs()[i]);
    }
    b.getNTKernProj().getPermute() = a.getNTKernProj().getPermute();
  }

}

template
void CopyKernelMatrix<quadruple, double>(KernelMatrix<double> &b,
					 KernelMatrix<quadruple> &a);

template
void CopyKernelMatrix<complex<quadruple>,
		      complex<double> >(KernelMatrix<complex<double> > &b,
					KernelMatrix<complex<quadruple> > &a);
