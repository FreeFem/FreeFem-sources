/*! \file   C_Dsub.cpp
    \brief  routines for substiution of off-diagonal matrix with strip
    \author Atsushi Suzuki, Laboratoire Jacques-Louis Lions
    \date   Jun. 20th 2014
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

#include "Compiler/blas.hpp"
#include "Compiler/OptionLibrary.h"
#include "Driver/C_threads_tasks.hpp"
#include "Driver/C_Dsub.hpp"
#include <string>

template<typename T>
void C_Dsub_task_exec(void *arg_)
{
  list<C_Dsub_task<T>* >*arg = (list<C_Dsub_task<T>* >*)arg_;
  int k = 0;
  for (typename list<C_Dsub_task<T>* >::const_iterator it = arg->begin();
       it != arg->end(); ++it, k++) {
    (*it)->func(*it);
  }
}
template
void C_Dsub_task_exec<double>(void *arg_);

template
void C_Dsub_task_exec<quadruple>(void *arg_);

template
void C_Dsub_task_exec<complex<double> >(void *arg_);

template
void C_Dsub_task_exec<complex<quadruple> >(void *arg_);


// #define DEBUG_BLOCKSUBT

template<typename T>
void dsub_sym2sym_diag(C_Dsub_task<T> *arg)
{
  const int ir_bgn = arg->ir_bgn;
  const int ir_end = arg->ir_end;
  const int ir_bgn_src = arg->ir_bgn_src;
  T *dst_pt = arg->dst_mtrx->addrCoefBlock(arg->ir_block, arg->jc_block);
  const int dst_row = arg->dst_mtrx->nrowBlock(arg->ir_block, arg->jc_block);
  SquareBlockMatrix<T> *src_mtrx = arg->src_pt;
  const int nrow = ir_end - ir_bgn;
  const int iblock0 = src_mtrx->BlockIndex(ir_bgn_src);
  const int iblock1 = src_mtrx->BlockIndex(ir_bgn_src + nrow - 1);
  const int nrow_src_offset = src_mtrx->BlockOffset(ir_bgn_src);
  const T none(-1.0);

  if ((iblock0 + 1) == iblock1) {
    const int ir_mid_src = src_mtrx->IndexBlock(iblock1); // iblock1 * SIZE_B1;
    const int nrow0 = ir_mid_src - ir_bgn_src;
    const int nrow1 = nrow - nrow0;
    const int ir_mid = ir_bgn + nrow0;
    {
      int jj, j0, j1;
      const int nrow_src = src_mtrx->nrowBlock(iblock0, iblock0);
      T *src_pt = src_mtrx->addrCoefBlock(iblock0, iblock0);
      for (j1 = nrow_src_offset * (nrow_src + 1), 
	   j0 = ir_bgn * (dst_row + 1),    //D(ir_bgn,ir_bgn)-(ir_mid,ir_mid)
	   jj = ir_bgn; jj < ir_mid; jj++, // driving index
	   j0 += dst_row, j1 += nrow_src) { 
	blas_axpy<T>(nrow0, none, src_pt + j1, 1, dst_pt + j0, 1);
      }
    }
    {
      int jj, j0, j1, j2, j3, j4;
      const int nrow_src1 = src_mtrx->nrowBlock(iblock1, iblock1);
      const int nrow_src0 = src_mtrx->nrowBlock(iblock0, iblock1);
      T *src_pt1 = src_mtrx->addrCoefBlock(iblock1, iblock1);
      T *src_pt0 = src_mtrx->addrCoefBlock(iblock0, iblock1);
      for (j4 = 0, j3 = nrow_src_offset,
	   j0 = ir_bgn + ir_mid * dst_row, // (ir_bgn,ir_mid)-(ir_mid,ir_end)
	   j1 = ir_mid * (dst_row + 1),    //D(ir_mid,ir_mid)-(ir_end,ir_end)
	   j2 = ir_mid + ir_bgn * dst_row, //T(ir_bgn,ir_mid)-(ir_mid,ir_end)
	   jj = ir_mid; jj < ir_end; jj++, // driving index
	   j0 += dst_row, j1 += dst_row, j2++,
	   j3 += nrow_src0, j4 += nrow_src1) { 
	blas_axpy<T>(nrow0, none, src_pt0 + j3, 1, dst_pt + j0, 1);
	blas_axpy<T>(nrow1, none, src_pt1 + j4, 1, dst_pt + j1, 1);
	// lower part of Schur complement matrix
	// blas_axpy<T>(nrow0, none, src_pt0 + j3, 1, dst_pt + j2, dst_row);
      }
    }
  }
  else if (iblock0 == iblock1) {
    int jj, j0, j1;
    const int nrow_src = src_mtrx->nrowBlock(iblock0, iblock0);
    T *src_pt = src_mtrx->addrCoefBlock(iblock0, iblock0);
    for (j1 = nrow_src_offset * (nrow_src + 1), 
	 j0 = ir_bgn * (dst_row + 1),      //D(ir_bgn,ir_bgn)-(ir_end,ir_end)
	 jj = ir_bgn; jj < ir_end; jj++,   // driving index
	 j0 += dst_row, j1 += nrow_src) { 
      blas_axpy<T>(nrow, none, src_pt + j1, 1, dst_pt + j0, 1);
    } // loop : j0, j1
  }
  else {
    if (arg->verbose) {
      fprintf(arg->fp, "%s %d : bad block indices in diagonal %d %d\n", 
	      __FILE__, __LINE__, iblock0, iblock1);
    }
  }
}

template
void dsub_sym2sym_diag<double>(C_Dsub_task<double> *arg);

template
void dsub_sym2sym_diag<quadruple>(C_Dsub_task<quadruple> *arg);

template
void dsub_sym2sym_diag<complex<double> >(C_Dsub_task<complex<double> > *arg);

template
void dsub_sym2sym_diag<complex<quadruple> >(C_Dsub_task<complex<quadruple> > *arg);
// end of template function dsub_sym2sym_diag<T>
 
template<typename T>
void dsub_unsym2unsym_diag(C_Dsub_task<T> *arg)
{
  const int ir_bgn = arg->ir_bgn;
  const int ir_end = arg->ir_end;
  const int ir_bgn_src = arg->ir_bgn_src;
  T *dst_pt = arg->dst_mtrx->addrCoefBlock(arg->ir_block, arg->jc_block);
  const int dst_row = arg->dst_mtrx->nrowBlock(arg->ir_block, arg->jc_block);
  SquareBlockMatrix<T> *src_mtrx = arg->src_pt;
  const int nrow = ir_end - ir_bgn;
  const int iblock0 = src_mtrx->BlockIndex(ir_bgn_src);
  const int iblock1 = src_mtrx->BlockIndex(ir_bgn_src + nrow - 1);
  const int nrow_src_offset = src_mtrx->BlockOffset(ir_bgn_src);
  // block copy : lower triangle to lower is enough
  const T none(-1.0);

  if ((iblock0 + 1) == iblock1) {
    const int ir_mid_src = src_mtrx->IndexBlock(iblock1);
    const int nrow0 = ir_mid_src - ir_bgn_src;
    const int nrow1 = nrow - nrow0;
    const int ir_mid = ir_bgn + nrow0;
    {
      int jj, j0, j1;
      const int nrow_src = src_mtrx->nrowBlock(iblock0, iblock0);
      T *src_pt = src_mtrx->addrCoefBlock(iblock0, iblock0);
      for (j1 = nrow_src_offset * (nrow_src + 1), 
	   j0 = ir_bgn * (dst_row + 1),    //D(ir_bgn,ir_bgn)-(ir_mid,ir_mid)
	   jj = ir_bgn; jj < ir_mid; jj++, // driving index
	   j0 += dst_row, j1 += nrow_src) { 
	blas_axpy<T>(nrow0, none, src_pt + j1, 1, dst_pt + j0, 1);
      }
    }
    {
      int jj, j0, j1, j2, j3, j4;
      const int nrow_src1 = src_mtrx->nrowBlock(iblock1, iblock1);
      const int nrow_src0 = src_mtrx->nrowBlock(iblock0, iblock1);
      T *src_pt1 = src_mtrx->addrCoefBlock(iblock1, iblock1);
      T *src_pt0 = src_mtrx->addrCoefBlock(iblock0, iblock1);
      // lower part of Schur complement matrix
      T *src_pt2 = src_mtrx->addrCoefBlock(iblock1, iblock0);
      for (j4 = 0, j3 = nrow_src_offset, 
  	   j0 = ir_bgn + ir_mid * dst_row, // (ir_bgn,ir_mid)-(ir_mid,ir_end) 
	   j1 = ir_mid * (dst_row + 1),    //D(ir_mid,ir_mid)-(ir_end,ir_end) 
	   j2 = ir_mid + ir_bgn * dst_row, //T(ir_bgn,ir_mid)-(ir_mid,ir_end)
	   jj = ir_mid; jj < ir_end; jj++, // driving index
	   j0 += dst_row, j1 += dst_row, j2++, 
	   j3 += nrow_src0, j4 += nrow_src1) { 
	blas_axpy<T>(nrow0, none, src_pt0 + j3, 1, dst_pt + j0, 1);
	blas_axpy<T>(nrow1, none, src_pt1 + j4, 1, dst_pt + j1, 1);
      // lower part of Schur complement matrix
	blas_axpy<T>(nrow0, none, src_pt2 + j3, 1, dst_pt + j2, dst_row);
      }
    }
  }
  else if (iblock0 == iblock1) {
    int jj, j0, j1;
    const int nrow_src = src_mtrx->nrowBlock(iblock0, iblock0);
    T *src_pt = src_mtrx->addrCoefBlock(iblock0, iblock0);
    for (j1 = nrow_src_offset * (nrow_src + 1), 
	 j0 = ir_bgn * (dst_row + 1),     //D(ir_bgn,ir_bgn)-(ir_end,ir_end)
	 jj = ir_bgn; jj < ir_end; jj++,  // driving index
	 j0 += dst_row, j1 += nrow_src) { 
      blas_axpy<T>(nrow, none, src_pt + j1, 1, dst_pt + j0, 1);
    } // loop : j0, j1
  }
  else {
    if (arg->verbose) {
      fprintf(arg->fp, "%s %d : bad block indices in diagonal %d %d\n", 
	      __FILE__, __LINE__, iblock0, iblock1);
    }
  }
}

template
void dsub_unsym2unsym_diag<double>(C_Dsub_task<double> *arg);

template
void dsub_unsym2unsym_diag<quadruple>(C_Dsub_task<quadruple> *arg);

template
void dsub_unsym2unsym_diag<complex<double> >(C_Dsub_task<complex<double> >*arg);

template
void dsub_unsym2unsym_diag<complex<quadruple> >(C_Dsub_task<complex<quadruple> >*arg);
// end of template function dsub_unsym2unsym_diag<T>
 
template<typename T>
void dsub_sym2sym(C_Dsub_task<T> *arg)
{
  const int ir_bgn = arg->ir_bgn;
  const int ir_end = arg->ir_end;
  const int jc_bgn = arg->jc_bgn;
  const int jc_end = arg->jc_end;
  const int ir_bgn_src = arg->ir_bgn_src;
  const int jc_bgn_src = arg->jc_bgn_src;
  T *dst_pt = arg->dst_mtrx->addrCoefBlock(arg->ir_block, arg->jc_block);
  const int dst_row = arg->dst_mtrx->nrowBlock(arg->ir_block, arg->jc_block);
  SquareBlockMatrix<T> *src_mtrx = arg->src_pt;
  const int nrow = ir_end - ir_bgn;
  const int ncol = jc_end - jc_bgn;
  const int iblock0 = src_mtrx->BlockIndex(ir_bgn_src);
  const int iblock1 = src_mtrx->BlockIndex(ir_bgn_src + nrow - 1);
  const int nrow_src_offset = src_mtrx->BlockOffset(ir_bgn_src);
  const int jblock0 = src_mtrx->BlockIndex(jc_bgn_src);
  const int jblock1 = src_mtrx->BlockIndex(jc_bgn_src + ncol - 1);
  const int ncol_src_offset = src_mtrx->BlockOffset(jc_bgn_src);
  const T none(-1.0);

  if ((iblock0 + 1) == iblock1) {
    // decomposition is defined by src block with SIZE_B1
    const int ir_mid_src = src_mtrx->IndexBlock(iblock1); // * SIZE_B1;
    const int nrow0 = ir_mid_src - ir_bgn_src;
    const int nrow1 = nrow - nrow0;
    const int ir_mid = ir_bgn + nrow0;  
    if ((jblock0 + 1) == jblock1) {
      const int jc_mid_src = src_mtrx->IndexBlock(jblock1); // * SIZE_B1;
      const int ncol0 = jc_mid_src - jc_bgn_src;
      const int jc_mid = jc_bgn + ncol0;
      {
	int jj, j0, j1, j2, j3;
	const int nrow_src0 = src_mtrx->nrowBlock(iblock0, jblock0);
	const int nrow_src1 = src_mtrx->nrowBlock(iblock1, jblock0);
	T *src_pt0 = src_mtrx->addrCoefBlock(iblock0, jblock0);
	T *src_pt1 = src_mtrx->addrCoefBlock(iblock1, jblock0);
	// distnation : transpose of 
	// distnation : transpose of 
	for (j3 = ncol_src_offset * nrow_src1, 
	     j2 = nrow_src_offset + ncol_src_offset * nrow_src0,
	     j0 = ir_bgn + jc_bgn * dst_row, // (ir_bgn,jc_bgn)-(ir_mid,jc_mid)
	     j1 = ir_mid + jc_bgn * dst_row, // (ir_mid,jc_bgn)-(ir_end,jc_mid)
	     jj = jc_bgn; jj < jc_mid; jj++, // drinving index
	     j0 += dst_row, j1 += dst_row, j2 += nrow_src0, j3 += nrow_src1) { 
	  blas_axpy<T>(nrow0, none, src_pt0 + j2, 1, dst_pt + j0, 1);
	  blas_axpy<T>(nrow1, none, src_pt1 + j3, 1, dst_pt + j1, 1);
	}
      }
      {
	int jj, j0, j1, j2, j3;
	const int nrow_src0 = src_mtrx->nrowBlock(iblock0, jblock1);
	const int nrow_src1 = src_mtrx->nrowBlock(iblock1, jblock1);
	T *src_pt0 = src_mtrx->addrCoefBlock(iblock0, jblock1);
	T *src_pt1 = src_mtrx->addrCoefBlock(iblock1, jblock1);
	// distnation : transpose of 
	// distnation : transpose of 
	for (j3 = 0, j2 = nrow_src_offset, 	      
	     j0 = ir_bgn + jc_mid * dst_row, // (ir_bgn,jc_mid)-(ir_mid,jc_end)
	     j1 = ir_mid + jc_mid * dst_row, // (ir_mid,jc_mid)-(ir_end,jc_end)
	     jj = jc_mid; jj < jc_end; jj++, // driving index
	     j0 += dst_row, j1 += dst_row, j2 += nrow_src0, j3 += nrow_src1) { 
	  blas_axpy<T>(nrow0, none, src_pt0 + j2, 1, dst_pt + j0, 1);
	  blas_axpy<T>(nrow1, none, src_pt1 + j3, 1, dst_pt + j1, 1);
	}
      }
    } // if ((jblock0 + 1) == jblock1) {
    else if (jblock0 == jblock1) {
      int jj, j0, j1, j2, j3;
      const int nrow_src0 = src_mtrx->nrowBlock(iblock0, jblock0);
      const int nrow_src1 = src_mtrx->nrowBlock(iblock1, jblock0);
      T *src_pt0 = src_mtrx->addrCoefBlock(iblock0, jblock0);
      T *src_pt1 = src_mtrx->addrCoefBlock(iblock1, jblock0);
	// distnation : transpose of 
	// distnation : transpose of 
      for (j3 = ncol_src_offset * nrow_src1, 
	   j2 = nrow_src_offset + ncol_src_offset * nrow_src0,
	   j0 = ir_bgn + jc_bgn * dst_row,  // (ir_bgn,jc_bgn)-(ir_mid,jc_end)
	   j1 = ir_mid + jc_bgn * dst_row,  // (ir_mid,jc_bgn)-(ir_end,jc_end)
	   jj = jc_bgn; jj < jc_end; jj++,  // driving index
	   j0 += dst_row, j1 += dst_row, j2 += nrow_src0, j3 += nrow_src1) { 
	blas_axpy<T>(nrow0, none, src_pt0 + j2, 1, dst_pt + j0, 1);
	blas_axpy<T>(nrow1, none, src_pt1 + j3, 1, dst_pt + j1, 1);
      }
    }
    else {
      if (arg->verbose) {
	fprintf(arg->fp, "%s %d : bad block indices in column : %d %d\n", 
		__FILE__, __LINE__, jblock0, jblock1);
      }
    }
  }
  else if (iblock0 == iblock1) {
    if ((jblock0 + 1) == jblock1) {
      const int jc_mid_src = src_mtrx->IndexBlock(jblock1);// * SIZE_B1;
      const int ncol0 = jc_mid_src - jc_bgn_src;
      const int jc_mid = jc_bgn + ncol0;
      {
	int jj, j0, j1;
	const int nrow_src = src_mtrx->nrowBlock(iblock0, jblock0);
	T *src_pt = src_mtrx->addrCoefBlock(iblock0, jblock0);
	for (j1 = nrow_src_offset + ncol_src_offset * nrow_src,
	     j0 = ir_bgn + jc_bgn * dst_row, // (ir_bgn,jc_bgn)-(ir_end,jc_mid)
	     jj = jc_bgn; jj < jc_mid; jj++, // driving index
	     j0 += dst_row, j1 += nrow_src) { 
	  blas_axpy<T>(nrow, none, src_pt + j1, 1, dst_pt + j0, 1);
	}
      }
      {
	int jj, j0, j1;
	const int nrow_src = src_mtrx->nrowBlock(iblock0, jblock1);
	T *src_pt = src_mtrx->addrCoefBlock(iblock0, jblock1);
	// distnation : transpose of 
	for (j1 = nrow_src_offset,
	     j0 = ir_bgn + jc_mid * dst_row, // (ir_bgn,jc_mid)-(ir_end,jc_end)
	     jj = jc_mid; jj < jc_end; jj++, // driving index
	     j0 += dst_row, j1 += nrow_src) {
	  blas_axpy<T>(nrow, none, src_pt + j1, 1, dst_pt + j0, 1);
	}
      }
    }
    else {
      int jj, j0, j1;
      const int nrow_src = src_mtrx->nrowBlock(iblock0, jblock0);
      T *src_pt = src_mtrx->addrCoefBlock(iblock0, jblock0);
      // distnation : transpose of 
      for (j1 = nrow_src_offset + ncol_src_offset * nrow_src,
	   j0 = ir_bgn + jc_bgn * dst_row, // (ir_bgn,jc_bgn)-(ir_end,jc_end)
	   jj = jc_bgn; jj < jc_end; jj++, // driving index
	   j0 += dst_row, j1 += nrow_src) { 
	blas_axpy<T>(nrow, none, src_pt + j1, 1, dst_pt + j0, 1);
      }
    }
  } // if ((jblock0 + 1) == jblock1) 
  else {
    if (arg->verbose) {
      fprintf(arg->fp, "%s %d : bad block indices in row : %d %d\n",
	      __FILE__, __LINE__, iblock0, iblock1);
    }
  }
}
template
void dsub_sym2sym<double>(C_Dsub_task<double> *arg);

template
void dsub_sym2sym<quadruple>(C_Dsub_task<quadruple> *arg);

template
void dsub_sym2sym<complex<double> >(C_Dsub_task<complex<double> > *arg);

template
void dsub_sym2sym<complex<quadruple> >(C_Dsub_task<complex<quadruple> > *arg);
// end of template function dsub_sym2sym<T>

template<typename T>
void dsub_unsym2unsym(C_Dsub_task<T> *arg)
{
  const int ir_bgn = arg->ir_bgn;
  const int ir_end = arg->ir_end;
  const int jc_bgn = arg->jc_bgn;
  const int jc_end = arg->jc_end;
  const int ir_bgn_src = arg->ir_bgn_src;
  const int jc_bgn_src = arg->jc_bgn_src;
  T *dst_pt =  arg->dst_mtrx->addrCoefBlock(arg->ir_block, arg->jc_block);
  const int dst_row =   arg->dst_mtrx->nrowBlock(arg->ir_block, arg->jc_block);
  T *dst_pt2 = arg->dst_mtrx->addrCoefBlock(arg->jc_block, arg->ir_block);
  SquareBlockMatrix<T> *src_mtrx = arg->src_pt;
  const int nrow = ir_end - ir_bgn;
  const int ncol = jc_end - jc_bgn;
  const int iblock0 = src_mtrx->BlockIndex(ir_bgn_src);
  const int iblock1 = src_mtrx->BlockIndex(ir_bgn_src + nrow - 1);
  const int nrow_src_offset = src_mtrx->BlockOffset(ir_bgn_src);
  const int jblock0 = src_mtrx->BlockIndex(jc_bgn_src);
  const int jblock1 = src_mtrx->BlockIndex(jc_bgn_src + ncol - 1);
  const int ncol_src_offset = src_mtrx->BlockOffset(jc_bgn_src);
  const T none(-1.0);

  if ((iblock0 + 1) == iblock1) {
    // decomposition is defined by src block with SIZE_B1
    const int ir_mid_src = src_mtrx->IndexBlock(iblock1);// * SIZE_B1;
    const int nrow0 = ir_mid_src - ir_bgn_src;
    const int nrow1 = nrow - nrow0;
    const int ir_mid = ir_bgn + nrow0;  
    if ((jblock0 + 1) == jblock1) {
      const int jc_mid_src = src_mtrx->IndexBlock(jblock1);// * SIZE_B1;
      const int ncol0 = jc_mid_src - jc_bgn_src;
      const int jc_mid = jc_bgn + ncol0;
      {
	int jj, j0, j1, j2, j3;
	const int nrow_src0 = src_mtrx->nrowBlock(iblock0, jblock0);
	const int nrow_src1 = src_mtrx->nrowBlock(iblock1, jblock0);
	T *src_pt0 = src_mtrx->addrCoefBlock(iblock0, jblock0);
	T *src_pt1 = src_mtrx->addrCoefBlock(iblock1, jblock0);
	for (j3 = ncol_src_offset * nrow_src1, 
	     j2 = nrow_src_offset + ncol_src_offset * nrow_src0,
	     j0 = ir_bgn + jc_bgn * dst_row, //(ir_bgn,jc_bgn)-(ir_mid,jc_mid)
	     j1 = ir_mid + jc_bgn * dst_row, //(ir_mid,jc_bgn)-(ir_end,jc_mid)
             jj = jc_bgn; jj < jc_mid; jj++, // drinving index
	     j0 += dst_row, j1 += dst_row,
             j2 += nrow_src0, j3 += nrow_src1) { 
	  blas_axpy<T>(nrow0, none, src_pt0 + j2, 1, dst_pt + j0, 1);
	  blas_axpy<T>(nrow1, none, src_pt1 + j3, 1, dst_pt + j1, 1);
	}
      }
      // lower part of Schur complement matrix
      {  	
	int jj, j0, j2;
	const int nrow_src0 = src_mtrx->nrowBlock(iblock0, jblock0);
	T *src_pt2 = src_mtrx->addrCoefBlock(jblock0, iblock0);
	{
	  for (j2 = nrow_src_offset + ncol_src_offset * nrow_src0,
	       j0 = ir_bgn + jc_bgn * dst_row, //(ir_bgn,jc_bgn)-(ir_mid,jc_mid)
               jj = jc_bgn; jj < jc_mid; jj++, // drinving index
	       j0 += dst_row,
	       j2 += nrow_src0) {
	    blas_axpy<T>(nrow0, none, src_pt2 + j2, 1, dst_pt2 + j0, 1);
	  }
	}
      }
      {
	int jj, j1, j3;
	const int nrow_src1 = src_mtrx->nrowBlock(iblock1, jblock0);
	T *src_pt3 = src_mtrx->addrCoefBlock(jblock0, iblock1);
	if (iblock1 == jblock0) {
	  for (j3 = ncol_src_offset,
	       j1 = ir_mid + jc_bgn * dst_row, //(ir_mid,jc_bgn)-(ir_end,jc_mid)
               jj = jc_bgn; jj < jc_mid; jj++, // drinving index
	       j1 += dst_row, 
               j3++) { 
	    blas_axpy<T>(nrow1, none, src_pt3 + j3, nrow_src1, dst_pt2 + j1, 1);
	  }
	}
	else {
	  for (j3 = ncol_src_offset * nrow_src1, 
	       j1 = ir_mid + jc_bgn * dst_row, //(ir_mid,jc_bgn)-(ir_end,jc_mid)
               jj = jc_bgn; jj < jc_mid; jj++, // drinving index
	       j1 += dst_row, 
               j3 += nrow_src1) { 
	    blas_axpy<T>(nrow1, none, src_pt3 + j3, 1, dst_pt2 + j1, 1);
	  }
	}
      }
      {
	int jj, j0, j1, j2, j3;
        const int nrow_src0 = src_mtrx->nrowBlock(iblock0, jblock1);
        const int nrow_src1 = src_mtrx->nrowBlock(iblock1, jblock1);
	T *src_pt0 = src_mtrx->addrCoefBlock(iblock0, jblock1);
	T *src_pt1 = src_mtrx->addrCoefBlock(iblock1, jblock1);
	for (j3 = 0, j2 = nrow_src_offset, 
	     j0 = ir_bgn + jc_mid * dst_row, //(ir_bgn,jc_mid)-(ir_mid,jc_end)
	     j1 = ir_mid + jc_mid * dst_row, //(ir_mid,jc_mid)-(ir_end,jc_end)
	     jj = jc_mid; jj < jc_end; jj++, //driving index
             j0 += dst_row, j1 += dst_row, 
	     j2 += nrow_src0, j3 += nrow_src1) { 
	  blas_axpy<T>(nrow0, none, src_pt0 + j2, 1, dst_pt + j0, 1);
	  blas_axpy<T>(nrow1, none, src_pt1 + j3, 1, dst_pt + j1, 1);
	}
      }
      {
	int jj, j0, j2;
        const int nrow_src0 = src_mtrx->nrowBlock(iblock0, jblock1);
	T *src_pt2 = src_mtrx->addrCoefBlock(jblock1, iblock0);
	{
	  for (j2 = nrow_src_offset, 
	       j0 = ir_bgn + jc_mid * dst_row, //(ir_bgn,jc_mid)-(ir_mid,jc_end)
	       jj = jc_mid; jj < jc_end; jj++, //driving index
               j0 += dst_row,
	       j2 += nrow_src0) { 
	    blas_axpy<T>(nrow0, none, src_pt2 + j2, 1, dst_pt2 + j0, 1);
	  }
	}
      }
      {
	int jj, j1, j3;
        const int nrow_src1 = src_mtrx->nrowBlock(iblock1, jblock1);
	T *src_pt3 = src_mtrx->addrCoefBlock(jblock1, iblock1);
	{
	  for (j3 = 0, 
	       j1 = ir_mid + jc_mid * dst_row, //(ir_mid,jc_mid)-(ir_end,jc_end)
	       jj = jc_mid; jj < jc_end; jj++, //driving index
	       j1 += dst_row, 
	       j3 += nrow_src1) { 
	    blas_axpy<T>(nrow1, none, src_pt3 + j3, 1, dst_pt2 + j1, 1);
	  }
	}
      }
    } // if ((jblock0 + 1) == jblock1) {
    else if (jblock0 == jblock1) {
      {
	int jj, j0, j1, j2, j3;
	const int nrow_src0 = src_mtrx->nrowBlock(iblock0, jblock0);
	const int nrow_src1 = src_mtrx->nrowBlock(iblock1, jblock0);
	T *src_pt0 = src_mtrx->addrCoefBlock(iblock0, jblock0);
	T *src_pt1 = src_mtrx->addrCoefBlock(iblock1, jblock0);
	for (j3 = ncol_src_offset * nrow_src1, 
	     j2 = nrow_src_offset + ncol_src_offset * nrow_src0,
  	     j0 = ir_bgn + jc_bgn * dst_row, // (ir_bgn,jc_bgn)-(ir_mid,jc_end)
             j1 = ir_mid + jc_bgn * dst_row, // (ir_mid,jc_bgn)-(ir_end,jc_end)
             jj = jc_bgn; jj < jc_end; jj++, // driving index
	     j0 += dst_row, j1 += dst_row, 
             j2 += nrow_src0, j3 += nrow_src1) { 
	  blas_axpy<T>(nrow0, none, src_pt0 + j2, 1, dst_pt + j0, 1);
	  blas_axpy<T>(nrow1, none, src_pt1 + j3, 1, dst_pt + j1, 1);
	}
      }
      {
	int jj, j0, j2;
	const int nrow_src0 = src_mtrx->nrowBlock(iblock0, jblock0);
	T *src_pt2 = src_mtrx->addrCoefBlock(jblock0, iblock0);
	{
	  for (j2 = nrow_src_offset + ncol_src_offset * nrow_src0,
	       j0 = ir_bgn + jc_bgn * dst_row, //(ir_bgn,jc_bgn)-(ir_mid,jc_end)
               jj = jc_bgn; jj < jc_end; jj++, // driving index
	       j0 += dst_row, 
	       j2 += nrow_src0) {
	    blas_axpy<T>(nrow0, none, src_pt2 + j2, 1, dst_pt2 + j0, 1);
	  }
	}
      }
      {
	int jj, j1, j3;
	const int nrow_src1 = src_mtrx->nrowBlock(iblock1, jblock0);
	T *src_pt3 = src_mtrx->addrCoefBlock(jblock0, iblock1);
	if (iblock1 == jblock0) {
	  for (j3 = ncol_src_offset,           // transposed
               j1 = ir_mid + jc_bgn * dst_row, //(ir_mid,jc_bgn)-(ir_end,jc_end)
               jj = jc_bgn; jj < jc_end; jj++, // driving index
	       j1 += dst_row, 
	       j3++) {
	    blas_axpy<T>(nrow1, none, src_pt3 + j3, nrow_src1, dst_pt2 + j1, 1);
	  }
	}
	else {
	  for (j3 = ncol_src_offset * nrow_src1, 
               j1 = ir_mid + jc_bgn * dst_row, //(ir_mid,jc_bgn)-(ir_end,jc_end)
               jj = jc_bgn; jj < jc_end; jj++, // driving index
	       j1 += dst_row, 
               j3 += nrow_src1) { 
	    blas_axpy<T>(nrow1, none, src_pt3 + j3, 1, dst_pt2 + j1, 1);
	  }
	}
      }
    }
    else {
      if (arg->verbose) {
	fprintf(arg->fp, "%s %d : bad block indices in column : %d %d\n", 
		__FILE__, __LINE__, jblock0, jblock1);
      }
    }
  }
  else if (iblock0 == iblock1) {
    if ((jblock0 + 1) == jblock1) {
      const int jc_mid_src = src_mtrx->IndexBlock(jblock1);// * SIZE_B1;
      const int ncol0 = jc_mid_src - jc_bgn_src;
      const int jc_mid = jc_bgn + ncol0;
      {
	int jj, j0, j1;
	const int nrow_src = src_mtrx->nrowBlock(iblock0, jblock0);
	T *src_pt = src_mtrx->addrCoefBlock(iblock0, jblock0);
	for (j1 = nrow_src_offset + ncol_src_offset * nrow_src,
	     j0 = ir_bgn + jc_bgn * dst_row, //(ir_bgn,jc_bgn)-(ir_end,jc_mid)
	     jj = jc_bgn; jj < jc_mid; jj++, // driving index
	     j0 += dst_row, j1 += nrow_src) { 
	  blas_axpy<T>(nrow, none, src_pt + j1, 1, dst_pt + j0, 1);
	}
      }
      {
	int jj, j0, j1;
	const int nrow_src = src_mtrx->nrowBlock(iblock0, jblock0);
	T *src_pt2 = src_mtrx->addrCoefBlock(jblock0, iblock0);
	if (iblock0 == jblock0) {
	  for (j1 = ncol_src_offset + nrow_src_offset * nrow_src, // transposed
	       j0 = ir_bgn + jc_bgn * dst_row, //(ir_bgn,jc_bgn)-(ir_end,jc_mid)
	       jj = jc_bgn; jj < jc_mid; jj++, // driving index
	       j0 += dst_row, j1++) { 
	    blas_axpy<T>(nrow, none, src_pt2 + j1, nrow_src, dst_pt2 + j0, 1);
	  }
	}
	else {
	  for (j1 = nrow_src_offset + ncol_src_offset * nrow_src,
	       j0 = ir_bgn + jc_bgn * dst_row, //(ir_bgn,jc_bgn)-(ir_end,jc_mid)
	       jj = jc_bgn; jj < jc_mid; jj++, // driving index
	       j0 += dst_row, j1 += nrow_src) { 
	    blas_axpy<T>(nrow, none, src_pt2 + j1, 1, dst_pt2 + j0, 1);
	  }
	}
      }
      {
	int jj, j0, j1;
	const int nrow_src = src_mtrx->nrowBlock(iblock0, jblock1);
	T *src_pt = src_mtrx->addrCoefBlock(iblock0, jblock1);
	for (j1 = nrow_src_offset,
	     j0 = ir_bgn + jc_mid * dst_row, // (ir_bgn,jc_mid)-(ir_end,jc_end)
	     jj = jc_mid; jj < jc_end; jj++, // driving index
	     j0 += dst_row, j1 += nrow_src) {
	  blas_axpy<T>(nrow, none, src_pt + j1, 1, dst_pt + j0, 1);
	}
      }
      {
	int jj, j0, j1;
	const int nrow_src = src_mtrx->nrowBlock(iblock0, jblock1);
	T *src_pt2 = src_mtrx->addrCoefBlock(jblock1, iblock0);
	{
	  for (j1 = nrow_src_offset,
	       j0 = ir_bgn + jc_mid * dst_row, //(ir_bgn,jc_mid)-(ir_end,jc_end)
	       jj = jc_mid; jj < jc_end; jj++, // driving index
	       j0 += dst_row, j1 += nrow_src) {
	    blas_axpy<T>(nrow, none, src_pt2 + j1, 1, dst_pt2 + j0, 1);
	  }
	}
      }
    }
    else {
      {
	int jj, j0, j1;
	const int nrow_src = src_mtrx->nrowBlock(iblock0, jblock0);
	T *src_pt = src_mtrx->addrCoefBlock(iblock0, jblock0);
	for (j1 = nrow_src_offset + ncol_src_offset * nrow_src,
	   j0 = ir_bgn + jc_bgn * dst_row, // (ir_bgn,jc_bgn)-(ir_end,jc_end)
	   jj = jc_bgn; jj < jc_end; jj++, // driving index
	   j0 += dst_row, j1 += nrow_src) { 
	  blas_axpy<T>(nrow, none, src_pt + j1, 1, dst_pt + j0, 1);
	}
      }
      {
	int jj, j0, j1;
	const int nrow_src = src_mtrx->nrowBlock(iblock0, jblock0);
	T *src_pt2 = src_mtrx->addrCoefBlock(jblock0, iblock0);
	if (iblock0 == jblock0) {
	  for (j1 = ncol_src_offset + nrow_src_offset * nrow_src, // transposed
	       j0 = ir_bgn + jc_bgn * dst_row, //(ir_bgn,jc_bgn)-(ir_end,jc_end)
	       jj = jc_bgn; jj < jc_end; jj++, // driving index
	       j0 += dst_row, j1++) { 
	    blas_axpy<T>(nrow, none, src_pt2 + j1, nrow_src, dst_pt2 + j0, 1);
	  }
	}
	else {
	  for (j1 = nrow_src_offset + ncol_src_offset * nrow_src,
	       j0 = ir_bgn + jc_bgn * dst_row, //(ir_bgn,jc_bgn)-(ir_end,jc_end)
	       jj = jc_bgn; jj < jc_end; jj++, // driving index
	       j0 += dst_row, j1 += nrow_src) { 
	    blas_axpy<T>(nrow, none, src_pt2 + j1, 1, dst_pt2 + j0, 1);
	  }
	}
      }
    }
  } // if ((jblock0 + 1) == jblock1) 
  else {
    if (arg->verbose) {
      fprintf(arg->fp, "%s %d : bad block indices in row : %d %d\n",
	      __FILE__, __LINE__, iblock0, iblock1);
    }
  }
}
template
void dsub_unsym2unsym<double>(C_Dsub_task<double> *arg);

template
void dsub_unsym2unsym<quadruple>(C_Dsub_task<quadruple> *arg);

template
void dsub_unsym2unsym<complex<double> >(C_Dsub_task<complex<double> > *arg);

template
void dsub_unsym2unsym<complex<quadruple> >(C_Dsub_task<complex<quadruple> > *arg);
// end of template function dsub_unsym2unsym<T>

template<typename T>
void dsub_unsym2diag(C_Dsub_task<T> *arg)
{
  // destination is inside of a diagonal block : ir_block == jc_block
  const int ir_bgn = arg->ir_bgn;
  const int ir_end = arg->ir_end;
  const int jc_bgn = arg->jc_bgn;
  const int jc_end = arg->jc_end;
  const int ir_bgn_src = arg->ir_bgn_src;
  const int jc_bgn_src = arg->jc_bgn_src;
  T *dst_pt =  arg->dst_mtrx->addrCoefBlock(arg->ir_block, arg->jc_block);
  const int dst_row =   arg->dst_mtrx->nrowBlock(arg->ir_block, arg->jc_block);
  SquareBlockMatrix<T> *src_mtrx = arg->src_pt;
  const int nrow = ir_end - ir_bgn;
  const int ncol = jc_end - jc_bgn;
  const int iblock0 = src_mtrx->BlockIndex(ir_bgn_src);
  const int iblock1 = src_mtrx->BlockIndex(ir_bgn_src + nrow - 1);
  const int nrow_src_offset = src_mtrx->BlockOffset(ir_bgn_src);
  const int jblock0 = src_mtrx->BlockIndex(jc_bgn_src);
  const int jblock1 = src_mtrx->BlockIndex(jc_bgn_src + ncol - 1);
  const int ncol_src_offset = src_mtrx->BlockOffset(jc_bgn_src);
  const T none(-1.0);

  if ((iblock0 + 1) == iblock1) {
    // decomposition is defined by src block with SIZE_B1
    const int ir_mid_src = src_mtrx->IndexBlock(iblock1);
    const int nrow0 = ir_mid_src - ir_bgn_src;
    const int nrow1 = nrow - nrow0;
    const int ir_mid = ir_bgn + nrow0;  
    if ((jblock0 + 1) == jblock1) {
      const int jc_mid_src = src_mtrx->IndexBlock(jblock1);
      const int ncol0 = jc_mid_src - jc_bgn_src;
      const int jc_mid = jc_bgn + ncol0;
      {
	int jj, j0, j1, j2, j3, j4, j5;
	const int nrow_src0 = src_mtrx->nrowBlock(iblock0, jblock0);
	const int nrow_src1 = src_mtrx->nrowBlock(iblock1, jblock0);
	T *src_pt0 = src_mtrx->addrCoefBlock(iblock0, jblock0);
	T *src_pt1 = src_mtrx->addrCoefBlock(iblock1, jblock0);
	// lower part of Schur complement matrix
	for (j5 = ncol_src_offset * nrow_src1, 
	     j4 = nrow_src_offset + ncol_src_offset * nrow_src0,
	     j0 = ir_bgn + jc_bgn * dst_row, // (ir_bgn,jc_bgn)-(ir_mid,jc_mid)
	     j1 = ir_mid + jc_bgn * dst_row, // (ir_mid,jc_bgn)-(ir_end,jc_mid)
	     jj = jc_bgn; jj < jc_mid; jj++, // drinving index
	     j0 += dst_row, j1 += dst_row, 
	     j4 += nrow_src0, j5 += nrow_src1) { 
	    blas_axpy<T>(nrow0, none, src_pt0 + j4, 1, dst_pt + j0, 1);
	    blas_axpy<T>(nrow1, none, src_pt1 + j5, 1, dst_pt + j1, 1);
	}
	// lower part of Schur complement matrix
	T *src_pt2 = src_mtrx->addrCoefBlock(jblock0, iblock0);
	T *src_pt3 = src_mtrx->addrCoefBlock(jblock0, iblock1);
	if (iblock1 == jblock0) {
	  for (j5 = ncol_src_offset,   // transposted
	       j4 = nrow_src_offset + ncol_src_offset * nrow_src0,
	       j2 = jc_bgn + ir_bgn * dst_row,//(ir_bgn,jc_bgn)-(ir_mid,jc_mid)
	       j3 = jc_bgn + ir_mid * dst_row,//(ir_mid,jc_bgn)-(ir_end,jc_mid)
	       jj = jc_bgn; jj < jc_mid; jj++, // drinving index
               j2++, j3++, 
	       j4 += nrow_src0, j5++) { 
	      blas_axpy<T>(nrow0, none, src_pt2 + j4, 1,         // transposed
			                dst_pt + j2, dst_row);
	      blas_axpy<T>(nrow1, none, src_pt3 + j5, nrow_src1, // both in the
			                dst_pt + j3, dst_row); // diagonal
	  }
	}
	else {
	  for (j5 = ncol_src_offset * nrow_src1, 
	       j4 = nrow_src_offset + ncol_src_offset * nrow_src0,
	       j2 = jc_bgn + ir_bgn * dst_row,//(ir_bgn,jc_bgn)-(ir_mid,jc_mid)
	       j3 = jc_bgn + ir_mid * dst_row,//(ir_mid,jc_bgn)-(ir_end,jc_mid)
	       jj = jc_bgn; jj < jc_mid; jj++, // drinving index
	       j2++, j3++, 
	       j4 += nrow_src0, j5 += nrow_src1) { 
	      blas_axpy<T>(nrow0, none, src_pt2 + j4, 1,         // transposed
                                        dst_pt + j2, dst_row);
	      blas_axpy<T>(nrow1, none, src_pt3 + j5, 1,         // transposed
			                dst_pt + j3, dst_row);
	  }
	}
      }
      {
	int jj, j0, j1, j2, j3, j4, j5;
	const int nrow_src0 = src_mtrx->nrowBlock(iblock0, jblock1);
	const int nrow_src1 = src_mtrx->nrowBlock(iblock1, jblock1);
	T *src_pt0 = src_mtrx->addrCoefBlock(iblock0, jblock1);
	T *src_pt1 = src_mtrx->addrCoefBlock(iblock1, jblock1);
	for (j5 = 0, j4 = nrow_src_offset, 
	     j0 = ir_bgn + jc_mid * dst_row, // (ir_bgn,jc_mid)-(ir_mid,jc_end)
	     j1 = ir_mid + jc_mid * dst_row, // (ir_mid,jc_mid)-(ir_end,jc_end)
	     jj = jc_mid; jj < jc_end; jj++, // driving index
	     j0 += dst_row, j1 += dst_row , 
	     j4 += nrow_src0, j5 += nrow_src1) { 
	  blas_axpy<T>(nrow0, none, src_pt0 + j4, 1, dst_pt + j0, 1);
	  blas_axpy<T>(nrow1, none, src_pt1 + j5, 1, dst_pt + j1, 1);
	}
	// lower part of Schur complement matrix
        T *src_pt2 = src_mtrx->addrCoefBlock(jblock1, iblock0);
	T *src_pt3 = src_mtrx->addrCoefBlock(jblock1, iblock1);
        {
	  for (j5 = 0, j4 = nrow_src_offset, 
	       j2 = jc_mid + ir_bgn * dst_row,//(ir_bgn,jc_mid)-(ir_mid,jc_end)
	       j3 = jc_mid + ir_mid * dst_row,//(ir_mid,jc_mid)-(ir_end,jc_end)
	       jj = jc_mid; jj < jc_end; jj++, // driving index
	       j2++, j3++, 
	       j4 += nrow_src0, j5 += nrow_src1) { 
	    blas_axpy<T>(nrow0, none, src_pt2 + j4, 1, dst_pt + j2, dst_row);
	    blas_axpy<T>(nrow1, none, src_pt3 + j5, 1, dst_pt + j3, dst_row);
	  }
	}
      }
    } // if ((jblock0 + 1) == jblock1) {
    else if (jblock0 == jblock1) {
      int jj, j0, j1, j2, j3, j4, j5;
      const int nrow_src0 = src_mtrx->nrowBlock(iblock0, jblock0);
      const int nrow_src1 = src_mtrx->nrowBlock(iblock1, jblock0);
      T *src_pt0 = src_mtrx->addrCoefBlock(iblock0, jblock0);
      T *src_pt1 = src_mtrx->addrCoefBlock(iblock1, jblock0);
      for (j5 = ncol_src_offset * nrow_src1, 
	   j4 = nrow_src_offset + ncol_src_offset * nrow_src0,
	   j0 = ir_bgn + jc_bgn * dst_row, // (ir_bgn,jc_bgn)-(ir_mid,jc_end)
	   j1 = ir_mid + jc_bgn * dst_row, // (ir_mid,jc_bgn)-(ir_end,jc_end)
	   jj = jc_bgn; jj < jc_end; jj++, // driving index
	   j0 += dst_row, j1 += dst_row, 
	   j4 += nrow_src0, j5 += nrow_src1) { 
	blas_axpy<T>(nrow0, none, src_pt0 + j4, 1, dst_pt + j0, 1);
	blas_axpy<T>(nrow1, none, src_pt1 + j5, 1, dst_pt + j1, 1);
      }
      // lower part of Schur complement matrix
      T *src_pt2 = src_mtrx->addrCoefBlock(jblock0, iblock0);
      T *src_pt3 = src_mtrx->addrCoefBlock(jblock0, iblock1);
      if (iblock1 == jblock0) {
	for (j5 = ncol_src_offset, 
	     j4 = nrow_src_offset + ncol_src_offset * nrow_src0,
	     j2 = jc_bgn + ir_bgn * dst_row,//(ir_bgn,jc_bgn)-(ir_mid,jc_end)
	     j3 = jc_bgn + ir_mid * dst_row,//(ir_mid,jc_bgn)-(ir_end,jc_end)
	     jj = jc_bgn; jj < jc_end; jj++, // driving index
	     j2++, j3++, 
	     j4 += nrow_src0, j5++) {
	  blas_axpy<T>(nrow0, none, src_pt2 + j4, 1,         // transposed
                                    dst_pt + j2, dst_row);
	  blas_axpy<T>(nrow1, none, src_pt3 + j5, nrow_src1, // both in the 
		                    dst_pt + j3, dst_row); // diagonal
	}
      }
      else {
	for (j5 = ncol_src_offset * nrow_src1, 
	     j4 = nrow_src_offset + ncol_src_offset * nrow_src0,
	     j2 = jc_bgn + ir_bgn * dst_row,//(ir_bgn,jc_bgn)-(ir_mid,jc_end)
	     j3 = jc_bgn + ir_mid * dst_row,//(ir_mid,jc_bgn)-(ir_end,jc_end)
	     jj = jc_bgn; jj < jc_end; jj++, // driving index
	     j2++, j3++, 
	     j4 += nrow_src0, j5 += nrow_src1) { 
	  blas_axpy<T>(nrow0, none, src_pt2 + j4, 1, dst_pt + j2, dst_row);
	  blas_axpy<T>(nrow1, none, src_pt3 + j5, 1, dst_pt + j3, dst_row);
	}
      }
    }
    else {
      if (arg->verbose) {
	fprintf(arg->fp, "%s %d : bad block indices in column : %d %d\n", 
		__FILE__, __LINE__, jblock0, jblock1);
      }
    }
  }
  else if (iblock0 == iblock1) {
    if ((jblock0 + 1) == jblock1) {
      const int jc_mid_src = src_mtrx->IndexBlock(jblock1);
      const int ncol0 = jc_mid_src - jc_bgn_src;
      const int jc_mid = jc_bgn + ncol0;
      {
	int jj, j0, j1, j2;
	const int nrow_src = src_mtrx->nrowBlock(iblock0, jblock0);
	T *src_pt = src_mtrx->addrCoefBlock(iblock0, jblock0);
	for (j2 = nrow_src_offset + ncol_src_offset * nrow_src,
	     j0 = ir_bgn + jc_bgn * dst_row, // (ir_bgn,jc_bgn)-(ir_end,jc_mid)
	     jj = jc_bgn; jj < jc_mid; jj++, // driving index
	     j0 += dst_row, j2 += nrow_src) { 
	  blas_axpy<T>(nrow, none, src_pt + j2, 1, dst_pt + j0, 1);
	}
      // lower part of Schur complement matrix
        T *src_pt2 = src_mtrx->addrCoefBlock(jblock0, iblock0);
	if (iblock0 == jblock0) {
	  for (j2 = ncol_src_offset + nrow_src_offset * nrow_src,
	       j1 = jc_bgn + ir_bgn * dst_row,//(ir_bgn,jc_bgn)-(ir_end,jc_mid)
	       jj = jc_bgn; jj < jc_mid; jj++, // driving index
	       j1++, j2++) { 
	    blas_axpy<T>(nrow, none, src_pt2 + j2, nrow_src,
                                     dst_pt + j1, dst_row);
	  }
	}
	else {
	  for (j2 = nrow_src_offset + ncol_src_offset * nrow_src,
	       j1 = jc_bgn + ir_bgn * dst_row,//(ir_bgn,jc_bgn)-(ir_end,jc_mid)
	       jj = jc_bgn; jj < jc_mid; jj++, // driving index
	       j1++, j2 += nrow_src) { 
	    blas_axpy<T>(nrow, none, src_pt2 + j2, 1,
                                     dst_pt + j1, dst_row);
	  }
	}
      }
      {
	int jj, j0, j1, j2;
	const int nrow_src = src_mtrx->nrowBlock(iblock0, jblock1);
	T *src_pt = src_mtrx->addrCoefBlock(iblock0, jblock1);
	T *src_pt2 = src_mtrx->addrCoefBlock(jblock1, iblock0);
	for (j2 = nrow_src_offset,
	     j0 = ir_bgn + jc_mid * dst_row, // (ir_bgn,jc_mid)-(ir_end,jc_end)
	     j1 = jc_mid + ir_bgn * dst_row,//T(ir_bgn,jc_mid)-(ir_end,jc_end)
	     jj = jc_mid; jj < jc_end; jj++, // driving index
	     j0 += dst_row, j1++, j2 += nrow_src) {
	  blas_axpy<T>(nrow, none, src_pt + j2, 1, dst_pt + j0, 1);
      // lower part of Schur complement matrix
	  blas_axpy<T>(nrow, none, src_pt2 + j2, 1, dst_pt + j1, dst_row);
	}
      }
    }
    else {
      int jj, j0, j1, j2;
      const int nrow_src = src_mtrx->nrowBlock(iblock0, jblock0);
      T *src_pt = src_mtrx->addrCoefBlock(iblock0, jblock0);
      for (j2 = nrow_src_offset + ncol_src_offset * nrow_src,
	   j0 = ir_bgn + jc_bgn * dst_row, // (ir_bgn,jc_bgn)-(ir_end,jc_end)
	   jj = jc_bgn; jj < jc_end; jj++, // driving index
	   j0 += dst_row, j2 += nrow_src) { 
	blas_axpy<T>(nrow, none, src_pt + j2, 1, dst_pt + j0, 1);
      }
      T *src_pt2 = src_mtrx->addrCoefBlock(jblock0, iblock0);
      if (iblock0 == jblock0) {
	for (j2 = ncol_src_offset + nrow_src_offset * nrow_src, //trans
	     j1 = jc_bgn + ir_bgn * dst_row, //T(ir_bgn,jc_bgn)-(ir_end,jc_end)
	     jj = jc_bgn; jj < jc_end; jj++, // driving index
	     j1++, j2++) { 
	  blas_axpy<T>(nrow, none, src_pt2 + j2, nrow_src, 
                                   dst_pt + j1, dst_row);
	}
      }
      else {
	for (j2 = nrow_src_offset + ncol_src_offset * nrow_src,
  	     j1 = jc_bgn + ir_bgn * dst_row, //T(ir_bgn,jc_bgn)-(ir_end,jc_end)
	     jj = jc_bgn; jj < jc_end; jj++, // driving index
	     j1++, j2 += nrow_src) { 
	  blas_axpy<T>(nrow, none, src_pt2 + j2, 1, dst_pt + j1, dst_row);
	}
      }
    }
  } // if ((jblock0 + 1) == jblock1) 
  else {
    if (arg->verbose) {
      fprintf(arg->fp, "%s %d : bad block indices in row : %d %d\n",
	      __FILE__, __LINE__, iblock0, iblock1);
    }
  }
}

template
void dsub_unsym2diag<double>(C_Dsub_task<double> *arg);

template
void dsub_unsym2diag<quadruple>(C_Dsub_task<quadruple> *arg);

template
void dsub_unsym2diag<complex<double> >(C_Dsub_task<complex<double> > *arg);

template
void dsub_unsym2diag<complex<quadruple> >(C_Dsub_task<complex<quadruple> > *arg);
// end of template function dsub_unsym2diag<T>

template<typename T>
void dsub_sym2rct(C_Dsub_task<T> *arg)
{
  const int ir_bgn = arg->ir_bgn;
  const int ir_end = arg->ir_end;
  const int jc_bgn = arg->jc_bgn;
  const int jc_end = arg->jc_end;
  const int ir_bgn_src = arg->ir_bgn_src;
  const int jc_bgn_src = arg->jc_bgn_src;
  //  const int dst_row = arg->dst_row;
  T *dst_pt = arg->dst_pt->addrCoefBlock(arg->ir_block, arg->jc_block);
  const int dst_row = arg->dst_pt->nrowBlock(arg->ir_block);
  SquareBlockMatrix<T> *src_mtrx = arg->src_pt;
  const int nrow = ir_end - ir_bgn;
  const int ncol = jc_end - jc_bgn;
  const int iblock0 = src_mtrx->BlockIndex(ir_bgn_src);
  const int iblock1 = src_mtrx->BlockIndex(ir_bgn_src + nrow - 1);
  const int nrow_src_offset = src_mtrx->BlockOffset(ir_bgn_src);
  const int jblock0 = src_mtrx->BlockIndex(jc_bgn_src);
  const int jblock1 = src_mtrx->BlockIndex(jc_bgn_src + ncol - 1);
  const int ncol_src_offset = src_mtrx->BlockOffset(jc_bgn_src);
  const T none(-1.0);
  if (iblock0 < iblock1) {
    // decomposition is defined by src block with SIZE_B1
    if ((jblock0 + 1) == jblock1) {
      const int jc_mid = jc_bgn + (src_mtrx->IndexBlock(jblock1) - jc_bgn_src);
      int ir_start = ir_bgn;
      {
	int jj, j0, j1;
	const int nrow_src = src_mtrx->nrowBlock(iblock0, jblock0); // SIZE_B1?
	const int nnrow = nrow_src - nrow_src_offset;
	T *src_pt;
	src_pt = src_mtrx->addrCoefBlock(iblock0, jblock0);
	for (j1 = nrow_src_offset + ncol_src_offset * nrow_src,
	     j0 = ir_start + jc_bgn * dst_row, // (   ,jc_bgn)-(  +nnrow,jc_mid)
	     jj = jc_bgn; jj < jc_mid; jj++,   // drinving index
	     j0 += dst_row, j1 += nrow_src) { 
	  blas_axpy<T>(nnrow, none, src_pt + j1, 1, dst_pt + j0, 1);
	}
	src_pt = src_mtrx->addrCoefBlock(iblock0, jblock1);
	for (j1 = nrow_src_offset,
	     j0 = ir_start + jc_mid * dst_row, // (   ,jc_mid)-(  +nnrow,jc_end)
	     jj = jc_mid; jj < jc_end; jj++, // drinving index
	     j0 += dst_row, j1 += nrow_src) { 
	  blas_axpy<T>(nnrow, none, src_pt + j1, 1, dst_pt + j0, 1);
	}
	ir_start += nnrow;
      }
      for (int iblock = iblock0 + 1; iblock < iblock1; iblock++) {
	int jj, j0, j1;
	const int nrow_src = src_mtrx->nrowBlock(iblock, jblock0); // SIZE_B1?
	T *src_pt;
	src_pt = src_mtrx->addrCoefBlock(iblock, jblock0);
	for (j1 = ncol_src_offset * nrow_src,
	     j0 = ir_start + jc_bgn * dst_row, // (  ,jc_bgn)-( +SIZE_B1,jc_mid)
	     jj = jc_bgn; jj < jc_mid; jj++, // drinving index
	     j0 += dst_row, j1 += nrow_src) { 
	  blas_axpy<T>(nrow_src, none, src_pt + j1, 1, dst_pt + j0, 1);
	}
	src_pt = src_mtrx->addrCoefBlock(iblock, jblock1);
	for (j1 = 0,
	     j0 = ir_start + jc_mid * dst_row, // (  ,jc_mid)-( +SIZE_B1,jc_end)
	     jj = jc_mid; jj < jc_end; jj++, // drinving index
	     j0 += dst_row, j1 += nrow_src) { 
	  blas_axpy<T>(nrow_src, none, src_pt + j1, 1, dst_pt + j0, 1);
	}
	ir_start += nrow_src;
      }
      {
	int jj, j0, j1;
	const int nrow_src = src_mtrx->nrowBlock(iblock1, jblock0);
	const int nnrow = ir_bgn_src + nrow - src_mtrx->IndexBlock(iblock1);
	T *src_pt;
	src_pt = src_mtrx->addrCoefBlock(iblock1, jblock0);
	for (j1 = ncol_src_offset * nrow_src,
	     j0 = ir_start + jc_bgn * dst_row,  // (  ,jc_bgn)-( +nnrow,jc_mid)
	     jj = jc_bgn; jj < jc_mid; jj++, // drinving index
	     j0 += dst_row, j1 += nrow_src) { 
	  blas_axpy<T>(nnrow, none, src_pt + j1, 1, dst_pt + j0, 1);
	}
	src_pt = src_mtrx->addrCoefBlock(iblock1, jblock1);
	for (j1 = 0,
	     j0 = ir_start + jc_mid * dst_row,  // (  ,jc_bgn)-( +nnrow,jc_mid)
	     jj = jc_mid; jj < jc_end; jj++, // drinving index
	     j0 += dst_row, j1 += nrow_src) { 
	  blas_axpy<T>(nnrow, none, src_pt + j1, 1, dst_pt + j0, 1);
	}
      }
    }
    else if (jblock0 == jblock1) {
      int ir_start = ir_bgn;
      {
	int jj, j0, j1;
	const int nrow_src = src_mtrx->nrowBlock(iblock0, jblock0);
	const int nnrow = nrow_src - nrow_src_offset;
	T *src_pt = src_mtrx->addrCoefBlock(iblock0, jblock0);
	for (j1 = nrow_src_offset + ncol_src_offset * nrow_src,
	     j0 = ir_start + jc_bgn * dst_row, // (  ,jc_bgn)-( +nnrow,jc_end)
	     jj = jc_bgn; jj < jc_end; jj++,   // drinving index
	     j0 += dst_row, j1 += nrow_src) { 
	  blas_axpy<T>(nnrow, none, src_pt + j1, 1, dst_pt + j0, 1);
	}
	ir_start += nnrow;
      }
      for (int iblock = iblock0 + 1; iblock < iblock1; iblock++) {
	int jj, j0, j1;
	const int nrow_src = src_mtrx->nrowBlock(iblock, jblock0);
	T *src_pt = src_mtrx->addrCoefBlock(iblock, jblock0);
	for (j1 = ncol_src_offset * nrow_src,
	     j0 = ir_start + jc_bgn * dst_row, // (  ,jc_bgn)-( +SIZE_B1,jc_end)
	     jj = jc_bgn; jj < jc_end; jj++, // drinving index
	     j0 += dst_row, j1 += nrow_src) { 
	  blas_axpy<T>(nrow_src, none, src_pt + j1, 1, dst_pt + j0, 1);
	}
	ir_start += nrow_src;
      }
      {
	int jj, j0, j1;
	const int nrow_src = src_mtrx->nrowBlock(iblock1, jblock0);
	//	const int itmp = (ir_bgn_src + nrow) % SIZE_B1;
	//	const int nnrow = (itmp == 0) ? SIZE_B1 : itmp;
	const int nnrow = ir_bgn_src + nrow - src_mtrx->IndexBlock(iblock1);
	T *src_pt = src_mtrx->addrCoefBlock(iblock1, jblock0);
	for (j1 = ncol_src_offset * nrow_src,
	     j0 = ir_start + jc_bgn * dst_row, // (  ,jc_bgn)-( +nnrow,jc_end)
	     jj = jc_bgn; jj < jc_end; jj++, // drinving index
	     j0 += dst_row, j1 += nrow_src) { 
	  blas_axpy<T>(nnrow, none, src_pt + j1, 1, dst_pt + j0, 1);
	}
      }
    }
    else {
      if (arg->verbose) {
	fprintf(arg->fp, "%s %d : bad block indices in column : %d %d\n", 
		__FILE__, __LINE__, jblock0, jblock1);
      }
    }
  }   // if (iblock0 < iblock1)
  else {
    if ((jblock0 + 1) == jblock1) {
      const int jc_mid_src = src_mtrx->IndexBlock(jblock1);
      const int ncol0 = jc_mid_src - jc_bgn_src;
      const int jc_mid = jc_bgn + ncol0;
      {
	int jj, j0, j1;
	const int nrow_src = src_mtrx->nrowBlock(iblock0, jblock0);
	T *src_pt = src_mtrx->addrCoefBlock(iblock0, jblock0);
	for (j1 = nrow_src_offset + ncol_src_offset * nrow_src,
	     j0 = ir_bgn + jc_bgn * dst_row, // (ir_bgn,jc_bgn)-(ir_end,jc_mid)
	     jj = jc_bgn; jj < jc_mid; jj++, // driving index
	     j0 += dst_row, j1 += nrow_src) { 
	  blas_axpy<T>(nrow, none, src_pt + j1, 1, dst_pt + j0, 1);
	}
      }
      {
	int jj, j0, j1;
	const int nrow_src = src_mtrx->nrowBlock(iblock0, jblock1);
	T *src_pt = src_mtrx->addrCoefBlock(iblock0, jblock1);
	for (j1 = nrow_src_offset,
	     j0 = ir_bgn + jc_mid * dst_row, // (ir_bgn,jc_mid)-(ir_end,jc_end)
	     jj = jc_mid; jj < jc_end; jj++, // driving index
	     j0 += dst_row, j1 += nrow_src) {
	  blas_axpy<T>(nrow, none, src_pt + j1, 1, dst_pt + j0, 1);
	}
      }
    }
    else if (jblock0 == jblock1) {
      int jj, j0, j1;
      const int nrow_src = src_mtrx->nrowBlock(iblock0, jblock0);
      T *src_pt = src_mtrx->addrCoefBlock(iblock0, jblock0);
      for (j1 = nrow_src_offset + ncol_src_offset * nrow_src,
	   j0 = ir_bgn + jc_bgn * dst_row,  // (ir_bgn,jc_bgn)-(ir_end,jc_end)
	   jj = jc_bgn; jj < jc_end; jj++,  // driving index
	   j0 += dst_row, j1 += nrow_src) {
	blas_axpy<T>(nrow, none, src_pt + j1, 1, dst_pt + j0, 1);
      }
    }
    else {
      if (arg->verbose) {
	fprintf(arg->fp, "%s %d : bad block indices in column : %d %d\n", 
		__FILE__, __LINE__, jblock0, jblock1);
      }
    }
  } // if (iblock0 < iblock1) 
}

template
void dsub_sym2rct<double>(C_Dsub_task<double> *arg);

template
void dsub_sym2rct<quadruple>(C_Dsub_task<quadruple> *arg);

template
void dsub_sym2rct<complex<double> >(C_Dsub_task<complex<double> > *arg);

template
void dsub_sym2rct<complex<quadruple> >(C_Dsub_task<complex<quadruple> > *arg);
// end of template function dsub_sym2rct<T>

template<typename T>
void dsub_unsym2rct(C_Dsub_task<T> *arg)
{
  const int ir_bgn = arg->ir_bgn;
  const int ir_end = arg->ir_end;
  const int jc_bgn = arg->jc_bgn;
  const int jc_end = arg->jc_end;
  const int ir_bgn_src = arg->ir_bgn_src;
  const int jc_bgn_src = arg->jc_bgn_src;
  //  const int dst_row = arg->dst_row;
  T *dst_pt = arg->dst_pt->addrCoefBlock(arg->ir_block, arg->jc_block);
  T *dst_pt2 = arg->dst_pt2->addrCoefBlock(arg->ir_block, arg->jc_block);
  const int dst_row = arg->dst_pt->nrowBlock(arg->ir_block);
  SquareBlockMatrix<T> *src_mtrx = arg->src_pt;
  const int nrow = ir_end - ir_bgn;
  const int ncol = jc_end - jc_bgn;
  const int iblock0 = src_mtrx->BlockIndex(ir_bgn_src);
  const int iblock1 = src_mtrx->BlockIndex(ir_bgn_src + nrow - 1);
  const int nrow_src_offset = src_mtrx->BlockOffset(ir_bgn_src);
  const int jblock0 = src_mtrx->BlockIndex(jc_bgn_src);
  const int jblock1 = src_mtrx->BlockIndex(jc_bgn_src + ncol - 1);
  const int ncol_src_offset = src_mtrx->BlockOffset(jc_bgn_src);
  const T none(-1.0);

  if (iblock0 < iblock1) {
    // decomposition is defined by src block with SIZE_B1
    if ((jblock0 + 1) == jblock1) {
      const int jc_mid = jc_bgn + (src_mtrx->IndexBlock(jblock1) - jc_bgn_src);
      int ir_start = ir_bgn;
      {
	int jj, j0, j1, j2;
	const int nrow_src = src_mtrx->nrowBlock(iblock0, jblock0); // SIZE_B1?
	const int nnrow = nrow_src - nrow_src_offset;
	T *src_pt0, *src_pt1;
	src_pt0 = src_mtrx->addrCoefBlock(iblock0, jblock0);
	src_pt1 = src_mtrx->addrCoefBlock(jblock0, iblock0);
	if (jblock0 == iblock0) {
	  for (j1 = nrow_src_offset + ncol_src_offset * nrow_src,
	       j2 = ncol_src_offset + nrow_src_offset * nrow_src, // transposed
	       j0 = ir_start + jc_bgn * dst_row,//(   ,jc_bgn)-(  +nnrow,jc_mid)
	       jj = jc_bgn; jj < jc_mid; jj++,   // drinving index
	       j0 += dst_row, j1 += nrow_src, j2++) { 
	    blas_axpy<T>(nnrow, none, src_pt0 + j1, 1, dst_pt + j0, 1);
	    blas_axpy<T>(nnrow, none, src_pt1 + j2, nrow_src, dst_pt2 + j0, 1);
	  }
	}
	else {
	  for (j1 = nrow_src_offset + ncol_src_offset * nrow_src,
	       j0 = ir_start + jc_bgn * dst_row,//(   ,jc_bgn)-(  +nnrow,jc_mid)
	       jj = jc_bgn; jj < jc_mid; jj++,   // drinving index
	       j0 += dst_row, j1 += nrow_src) { 
	    blas_axpy<T>(nnrow, none, src_pt0 + j1, 1, dst_pt + j0, 1);
	    blas_axpy<T>(nnrow, none, src_pt1 + j1, 1, dst_pt2 + j0, 1);
	  }
	}
	src_pt0 = src_mtrx->addrCoefBlock(iblock0, jblock1);
	src_pt1 = src_mtrx->addrCoefBlock(jblock1, iblock0); 
	if (jblock1 == iblock0) {
	  for (j1 = nrow_src_offset,
	       j2 = nrow_src_offset * nrow_src, // transposed
	       j0 = ir_start + jc_mid * dst_row,//(   ,jc_mid)-(  +nnrow,jc_end)
	       jj = jc_mid; jj < jc_end; jj++, // drinving index
	       j0 += dst_row, j1 += nrow_src, j2++) { 
	    blas_axpy<T>(nnrow, none, src_pt0 + j1, 1, dst_pt + j0, 1);
	    blas_axpy<T>(nnrow, none, src_pt1 + j2, nrow_src, dst_pt2 + j0, 1);
	  }
	}
	else {
	  for (j1 = nrow_src_offset,
	       j0 = ir_start + jc_mid * dst_row,//(   ,jc_mid)-(  +nnrow,jc_end)
	       jj = jc_mid; jj < jc_end; jj++, // drinving index
	       j0 += dst_row, j1 += nrow_src) { 
	    blas_axpy<T>(nnrow, none, src_pt0 + j1, 1, dst_pt + j0, 1);
	    blas_axpy<T>(nnrow, none, src_pt1 + j1, 1, dst_pt2 + j0, 1);
	  }
	}
	ir_start += nnrow;
      }
      for (int iblock = iblock0 + 1; iblock < iblock1; iblock++) {
	int jj, j0, j1, j2;
	const int nrow_src = src_mtrx->nrowBlock(iblock, jblock0); // SIZE_B1?
	T *src_pt0, *src_pt1;
	src_pt0 = src_mtrx->addrCoefBlock(iblock, jblock0);
	src_pt1 = src_mtrx->addrCoefBlock(jblock0, iblock); 
	if (jblock0 == iblock) {
	  for (j1 = ncol_src_offset * nrow_src,
	       j2 = ncol_src_offset,             // transposed
	       j0 = ir_start + jc_bgn * dst_row,//(  ,jc_bgn)-( +SIZE_B1,jc_mid)
	       jj = jc_bgn; jj < jc_mid; jj++, // drinving index
	       j0 += dst_row, j1 += nrow_src, j2++) { 
	    blas_axpy<T>(nrow_src, none, src_pt0 + j1, 1, 
                                         dst_pt + j0, 1);
	    blas_axpy<T>(nrow_src, none, src_pt1 + j2, nrow_src, 
			                 dst_pt2 + j0, 1);
	  }
	}
	else {
	  for (j1 = ncol_src_offset * nrow_src,
	       j0 = ir_start + jc_bgn * dst_row,//(  ,jc_bgn)-( +SIZE_B1,jc_mid)
	       jj = jc_bgn; jj < jc_mid; jj++, // drinving index
	       j0 += dst_row, j1 += nrow_src) { 
	    blas_axpy<T>(nrow_src, none, src_pt0 + j1, 1, dst_pt + j0, 1);
	    blas_axpy<T>(nrow_src, none, src_pt1 + j1, 1, dst_pt2 + j0, 1);
	  }
	}
	src_pt0 = src_mtrx->addrCoefBlock(iblock, jblock1);
	src_pt1 = src_mtrx->addrCoefBlock(jblock1, iblock); 
	if (jblock1 == iblock) {
	  for (j1 = 0,
	       j2 = 0,
       	       j0 = ir_start + jc_mid * dst_row,//(  ,jc_mid)-( +SIZE_B1,jc_end)
	       jj = jc_mid; jj < jc_end; jj++, // drinving index
	       j0 += dst_row, j1 += nrow_src, j2++) { 
	    blas_axpy<T>(nrow_src, none, src_pt0 + j1, 1, 
                                         dst_pt + j0, 1);
	    blas_axpy<T>(nrow_src, none, src_pt1 + j2, nrow_src, 
                                         dst_pt2 + j0, 1);
	  }
	}
	else {
	  for (j1 = 0,
       	       j0 = ir_start + jc_mid * dst_row,//(  ,jc_mid)-( +SIZE_B1,jc_end)
	       jj = jc_mid; jj < jc_end; jj++, // drinving index
	       j0 += dst_row, j1 += nrow_src) { 
	    blas_axpy<T>(nrow_src, none, src_pt0 + j1, 1, dst_pt + j0, 1);
	    blas_axpy<T>(nrow_src, none, src_pt1 + j1, 1, dst_pt2 + j0, 1);
	  }
	}
	ir_start += nrow_src;
      }
      {
	int jj, j0, j1, j2;
	const int nrow_src = src_mtrx->nrowBlock(iblock1, jblock0);
	//	const int itmp = (ir_bgn_src + nrow) % SIZE_B1;
	//	const int nnrow = (itmp == 0) ? SIZE_B1 : itmp;
	const int nnrow = ir_bgn_src + nrow - src_mtrx->IndexBlock(iblock1);
	T *src_pt0, *src_pt1;
	src_pt0 = src_mtrx->addrCoefBlock(iblock1, jblock0);
	src_pt1 = src_mtrx->addrCoefBlock(jblock0, iblock1);
	if (jblock0 == iblock1) {
	  for (j1 = ncol_src_offset * nrow_src,
	       j2 = ncol_src_offset,              // transposed
	       j0 = ir_start + jc_bgn * dst_row, // (  ,jc_bgn)-( +nnrow,jc_mid)
	       jj = jc_bgn; jj < jc_mid; jj++, // drinving index
	       j0 += dst_row, j1 += nrow_src, j2++) { 
	    blas_axpy<T>(nnrow, none, src_pt0 + j1, 1, dst_pt + j0, 1);
	    blas_axpy<T>(nnrow, none, src_pt1 + j2, nrow_src, dst_pt2 + j0, 1);
	  }
	}
	else {
	  for (j1 = ncol_src_offset * nrow_src,
	       j0 = ir_start + jc_bgn * dst_row, // (  ,jc_bgn)-( +nnrow,jc_mid)
	       jj = jc_bgn; jj < jc_mid; jj++, // drinving index
	       j0 += dst_row, j1 += nrow_src) { 
	    blas_axpy<T>(nnrow, none, src_pt0 + j1, 1, dst_pt + j0, 1);
	    blas_axpy<T>(nnrow, none, src_pt1 + j1, 1, dst_pt2 + j0, 1);
	  }
	}
	src_pt0 = src_mtrx->addrCoefBlock(iblock1, jblock1);
	src_pt1 = src_mtrx->addrCoefBlock(jblock1, iblock1); 
	if (jblock1 == iblock1) {
	  for (j2 = 0,
	       j1 = 0, //
	       j0 = ir_start + jc_mid * dst_row, // (  ,jc_bgn)-( +nnrow,jc_mid)
	       jj = jc_mid; jj < jc_end; jj++, // drinving index
	       j0 += dst_row, j1 += nrow_src, j2++) { 
	    blas_axpy<T>(nnrow, none, src_pt0 + j1, 1, dst_pt + j0, 1);
	    blas_axpy<T>(nnrow, none, src_pt1 + j2, nrow_src, dst_pt2 + j0, 1);
	  }
	}
	else {
	  for (j1 = 0,
	     j0 = ir_start + jc_mid * dst_row,  // (  ,jc_bgn)-( +nnrow,jc_mid)
	     jj = jc_mid; jj < jc_end; jj++, // drinving index
	     j0 += dst_row, j1 += nrow_src) { 
	    blas_axpy<T>(nnrow, none, src_pt0 + j1, 1, dst_pt + j0, 1);
	    blas_axpy<T>(nnrow, none, src_pt1 + j1, 1, dst_pt2 + j0, 1);
	  }
	}
      }
    }
    else if (jblock0 == jblock1) {
      int ir_start = ir_bgn;
      {
	int jj, j0, j1, j2;
	const int nrow_src = src_mtrx->nrowBlock(iblock0, jblock0);
	const int nnrow = nrow_src - nrow_src_offset;
	T *src_pt0 = src_mtrx->addrCoefBlock(iblock0, jblock0);
	T *src_pt1 = src_mtrx->addrCoefBlock(jblock0, iblock0);
	if (jblock0 == iblock0) {
	  for (j1 = nrow_src_offset + ncol_src_offset * nrow_src,
	       j2 = ncol_src_offset + nrow_src_offset * nrow_src, // transpopsed
  	       j0 = ir_start + jc_bgn * dst_row, // (  ,jc_bgn)-( +nnrow,jc_end)
	       jj = jc_bgn; jj < jc_end; jj++,   // drinving index
	       j0 += dst_row, j1 += nrow_src, j2++) { 
	    blas_axpy<T>(nnrow, none, src_pt0 + j1, 1, dst_pt + j0, 1);
	    blas_axpy<T>(nnrow, none, src_pt1 + j2, nrow_src, dst_pt2 + j0, 1);
	  }
	}
	else {
	  for (j1 = nrow_src_offset + ncol_src_offset * nrow_src,
	     j0 = ir_start + jc_bgn * dst_row, // (  ,jc_bgn)-( +nnrow,jc_end)
	     jj = jc_bgn; jj < jc_end; jj++,   // drinving index
	     j0 += dst_row, j1 += nrow_src) { 
	    blas_axpy<T>(nnrow, none, src_pt0 + j1, 1, dst_pt + j0, 1);
	    blas_axpy<T>(nnrow, none, src_pt1 + j1, 1, dst_pt2 + j0, 1);
	  }
	}
	ir_start += nnrow;
      }
      for (int iblock = iblock0 + 1; iblock < iblock1; iblock++) {
	int jj, j0, j1, j2;
	const int nrow_src = src_mtrx->nrowBlock(iblock, jblock0);
	T *src_pt0 = src_mtrx->addrCoefBlock(iblock, jblock0);
	T *src_pt1 = src_mtrx->addrCoefBlock(jblock0, iblock); 
	if (jblock0 == iblock) {
	  for (j1 = ncol_src_offset * nrow_src,
	       j2 = ncol_src_offset,            // transposed
	       j0 = ir_start + jc_bgn * dst_row,//(  ,jc_bgn)-( +SIZE_B1,jc_end)
	       jj = jc_bgn; jj < jc_end; jj++, // drinving index
	       j0 += dst_row, j1 += nrow_src, j2++) { 
	    blas_axpy<T>(nrow_src, none, src_pt0 + j1, 1, dst_pt + j0, 1);
	    blas_axpy<T>(nrow_src, none, src_pt1 + j2, nrow_src, 
                                         dst_pt2 + j0, 1);
	  }
	}
	else {
	  for (j1 = ncol_src_offset * nrow_src,
	       j0 = ir_start + jc_bgn * dst_row, // (  ,jc_bgn)-( +SIZE_B1,jc_end)
	       jj = jc_bgn; jj < jc_end; jj++, // drinving index
	       j0 += dst_row, j1 += nrow_src) { 
	    blas_axpy<T>(nrow_src, none, src_pt0 + j1, 1, dst_pt + j0, 1);
	    blas_axpy<T>(nrow_src, none, src_pt1 + j1, 1, dst_pt2 + j0, 1);
	  }
	}
	ir_start += nrow_src;
      }
      {
	int jj, j0, j1, j2;
	const int nrow_src = src_mtrx->nrowBlock(iblock1, jblock0);
	//	const int itmp = (ir_bgn_src + nrow) % SIZE_B1;
	//	const int nnrow = (itmp == 0) ? SIZE_B1 : itmp;
	const int nnrow = ir_bgn_src + nrow - src_mtrx->IndexBlock(iblock1);
	T *src_pt0 = src_mtrx->addrCoefBlock(iblock1, jblock0);
        T *src_pt1 = src_mtrx->addrCoefBlock(jblock0, iblock1); 
	if (jblock0 == iblock1) {
	  for (j1 = ncol_src_offset * nrow_src,
	       j2 = ncol_src_offset,             // transposed
	       j0 = ir_start + jc_bgn * dst_row, // (  ,jc_bgn)-( +nnrow,jc_end)
	       jj = jc_bgn; jj < jc_end; jj++, // drinving index
  	       j0 += dst_row, j1 += nrow_src, j2++) { 
	    blas_axpy<T>(nnrow, none, src_pt0 + j1, 1, dst_pt + j0, 1);
	    blas_axpy<T>(nnrow, none, src_pt1 + j2, nrow_src, dst_pt2 + j0, 1);
	  }
	}
	else {
	  for (j1 = ncol_src_offset * nrow_src,
	       j0 = ir_start + jc_bgn * dst_row, // (  ,jc_bgn)-( +nnrow,jc_end)
	       jj = jc_bgn; jj < jc_end; jj++, // drinving index
	       j0 += dst_row, j1 += nrow_src) { 
	    blas_axpy<T>(nnrow, none, src_pt0 + j1, 1, dst_pt + j0, 1);
	    blas_axpy<T>(nnrow, none, src_pt1 + j1, 1, dst_pt2 + j0, 1);
	  }
	}
      }
    }
    else {
      if (arg->verbose) {
	fprintf(arg->fp, "%s %d : bad block indices in column : %d %d\n", 
	      __FILE__, __LINE__, jblock0, jblock1);
      }
    }
  }
  else {
    if ((jblock0 + 1) == jblock1) {
      const int jc_mid_src = src_mtrx->IndexBlock(jblock1);
      const int ncol0 = jc_mid_src - jc_bgn_src;
      const int jc_mid = jc_bgn + ncol0;
      {
	int jj, j0, j1, j2;
	const int nrow_src = src_mtrx->nrowBlock(iblock0, jblock0);
	T *src_pt0 = src_mtrx->addrCoefBlock(iblock0, jblock0);
	T *src_pt1 = src_mtrx->addrCoefBlock(jblock0, iblock0);
	if (jblock0 == iblock0) {
	  for (j1 = nrow_src_offset + ncol_src_offset * nrow_src,
	       j2 = ncol_src_offset + nrow_src_offset * nrow_src, //transposed
	       j0 = ir_bgn + jc_bgn * dst_row, //(ir_bgn,jc_bgn)-(ir_end,jc_mid)
	       jj = jc_bgn; jj < jc_mid; jj++, // driving index
	       j0 += dst_row, j1 += nrow_src, j2++) { 
	    blas_axpy<T>(nrow, none, src_pt0 + j1, 1, dst_pt + j0, 1);
	    blas_axpy<T>(nrow, none, src_pt1 + j2, nrow_src, dst_pt2 + j0, 1);
	  }
	}
	else {
	  for (j1 = nrow_src_offset + ncol_src_offset * nrow_src,
	       j0 = ir_bgn + jc_bgn * dst_row, //(ir_bgn,jc_bgn)-(ir_end,jc_mid)
	       jj = jc_bgn; jj < jc_mid; jj++, // driving index
	       j0 += dst_row, j1 += nrow_src) { 
	    blas_axpy<T>(nrow, none, src_pt0 + j1, 1, dst_pt + j0, 1);
	    blas_axpy<T>(nrow, none, src_pt1 + j1, 1, dst_pt2 + j0, 1);
	  }
	}
      }
      {
	int jj, j0, j1, j2;
	const int nrow_src = src_mtrx->nrowBlock(iblock0, jblock1);
	T *src_pt0 = src_mtrx->addrCoefBlock(iblock0, jblock1);
	T *src_pt1 = src_mtrx->addrCoefBlock(jblock1, iblock0); 
	if (jblock1 == iblock0) {
	  for (j1 = nrow_src_offset,
               j2 = nrow_src_offset * nrow_src, //
  	       j0 = ir_bgn + jc_mid * dst_row, //(ir_bgn,jc_mid)-(ir_end,jc_end)
	       jj = jc_mid; jj < jc_end; jj++, // driving index
	       j0 += dst_row, j1 += nrow_src, j2++) {
	    blas_axpy<T>(nrow, none, src_pt0 + j1, 1, dst_pt + j0, 1);
	    blas_axpy<T>(nrow, none, src_pt1 + j2, nrow_src, dst_pt2 + j0, 1);
	  }
	}
	else {
	// distnation : 
	  for (j1 = nrow_src_offset,
  	       j0 = ir_bgn + jc_mid * dst_row, //(ir_bgn,jc_mid)-(ir_end,jc_end)
	       jj = jc_mid; jj < jc_end; jj++, // driving index
	       j0 += dst_row, j1 += nrow_src) {
	    blas_axpy<T>(nrow, none, src_pt0 + j1, 1, dst_pt + j0, 1);
	    blas_axpy<T>(nrow, none, src_pt1 + j1, 1, dst_pt2 + j0, 1);
	  }
	}
      }
    }
    else if (jblock0 == jblock1) {
      int jj, j0, j1, j2;
      const int nrow_src = src_mtrx->nrowBlock(iblock0, jblock0);
      T *src_pt0 = src_mtrx->addrCoefBlock(iblock0, jblock0);
      T *src_pt1 = src_mtrx->addrCoefBlock(jblock0, iblock0);
      if (jblock0 == iblock0) {
	for (j1 = nrow_src_offset + ncol_src_offset * nrow_src,
	     j2 = ncol_src_offset + nrow_src_offset * nrow_src,  
	     j0 = ir_bgn + jc_bgn * dst_row,  // (ir_bgn,jc_bgn)-(ir_end,jc_end)
	     jj = jc_bgn; jj < jc_end; jj++,  // driving index
	     j0 += dst_row, j1 += nrow_src, j2++) {
	  blas_axpy<T>(nrow, none, src_pt0 + j1, 1, dst_pt + j0, 1);
	  blas_axpy<T>(nrow, none, src_pt1 + j2, nrow_src, dst_pt2 + j0, 1);
	}
      }
      else {
	for (j1 = nrow_src_offset + ncol_src_offset * nrow_src,
	     j0 = ir_bgn + jc_bgn * dst_row,  // (ir_bgn,jc_bgn)-(ir_end,jc_end)
	     jj = jc_bgn; jj < jc_end; jj++,  // driving index
	     j0 += dst_row, j1 += nrow_src) {
	  blas_axpy<T>(nrow, none, src_pt0 + j1, 1, dst_pt + j0, 1);
	  blas_axpy<T>(nrow, none, src_pt1 + j1, 1, dst_pt2 + j0, 1);
	}
      }
    }
    else {
      if (arg->verbose) {
	fprintf(arg->fp, "%s %d : bad block indices in column : %d %d\n", 
		__FILE__, __LINE__, jblock0, jblock1);
      }
    }
  } // if (iblock0 < iblock1)
}

template
void dsub_unsym2rct<double>(C_Dsub_task<double> *arg);

template
void dsub_unsym2rct<quadruple>(C_Dsub_task<quadruple> *arg);

template
void dsub_unsym2rct<complex<double> >(C_Dsub_task<complex<double> > *arg);

template
void dsub_unsym2rct<complex<quadruple> >(C_Dsub_task<complex<quadruple> > *arg);
// end of template function dsub_unsym2rct<T>

template<typename T>
void dsub_sym2sym_diag_two(C_Dsub_task<T> *arg)
{
  if (arg->isSkip) {
    return;
  }  
  C_Dsub_task<T> *tmp;
  tmp = new C_Dsub_task<T>(arg->atomic_size,
			   arg->atomic_id,
			   arg->ir_bgn, 
			   arg->ir_end,
			   (-1), // jc_bgn
			   (-1), // jc_end
			   arg->ir_bgn_src,
			   (-1), // ir_bgn_src2
			   (-1), // jc_bgn_src
			   (-1), // jc_bgn_src2
			   arg->dst_row,
			   arg->dst_mtrx,
			   (RectBlockMatrix<T> *)NULL, // dst_pt
			   arg->ir_block,
			   arg->jc_block,
			   (RectBlockMatrix<T> *)NULL, // dst_pt2
			   (-1),
			   (-1),
			   arg->src_pt,  
			   (SquareBlockMatrix<T>*)NULL, // src_pt2
			   dsub_sym2sym_diag,
			   false,  // dummy
			   *(arg->ops_complexity),
			   arg->father_id,
			   arg->level,
			   arg->verbose,
			   arg->fp);
  dsub_sym2sym_diag<T>(tmp);
  delete tmp;
  tmp = new C_Dsub_task<T>(arg->atomic_size,
			   arg->atomic_id,
			   arg->ir_bgn, 
			   arg->ir_end,
			   (-1), // jc_bgn
			   (-1), // jc_end
			   arg->ir_bgn_src2,
			   (-1), // ir_bgn_src2
			   (-1), // jc_bgn_src
			   (-1), // jc_bgn_src2
			   arg->dst_row,
			   arg->dst_mtrx,
			   (RectBlockMatrix<T> *)NULL, // dst_pt2
			   arg->ir_block,
			   arg->jc_block,
			   (RectBlockMatrix<T> *)NULL, // dst_pt2
			   (-1),
			   (-1),
			   arg->src_pt2,  
			   (SquareBlockMatrix<T>*)NULL, // src_pt2
			   dsub_sym2sym_diag,
			   false,  // dummy
			   *(arg->ops_complexity),
			   arg->father_id,
			   arg->level,
			   arg->verbose,
			   arg->fp);
  dsub_sym2sym_diag<T>(tmp);
  delete tmp;
}

template
void dsub_sym2sym_diag_two<double>(C_Dsub_task<double> *arg);

template
void dsub_sym2sym_diag_two<quadruple>(C_Dsub_task<quadruple> *arg);

template
void dsub_sym2sym_diag_two<complex<double> >(C_Dsub_task<complex<double> > *arg);

template
void dsub_sym2sym_diag_two<complex<quadruple> >(C_Dsub_task<complex<quadruple> > *arg);
// end of template function dsub_sym2sym_diag_two<T> 

template<typename T>
void dsub_unsym2unsym_diag_two(C_Dsub_task<T> *arg)
{
  if (arg->isSkip) {
    return;
  }
  C_Dsub_task<T> *tmp;

  tmp = new C_Dsub_task<T>(arg->atomic_size,
			   arg->atomic_id,
			   arg->ir_bgn, 
			   arg->ir_end,
			   (-1), // jc_bgn
			   (-1), // jc_end
			   arg->ir_bgn_src,
			   (-1), // ir_bgn_src2
			   (-1), // jc_bgn_src
			   (-1), // jc_bgn_src2
			   arg->dst_row,
			   arg->dst_mtrx,
			   (RectBlockMatrix<T> *)NULL, // dst_pt2
			   arg->ir_block,
			   arg->jc_block,
			   (RectBlockMatrix<T> *)NULL, // dst_pt2
			   (-1),
			   (-1),
			   arg->src_pt,  
			   (SquareBlockMatrix<T>*)NULL, // src_pt2
			   dsub_unsym2unsym_diag,
			   false,  // dummy
			   *(arg->ops_complexity),
			   arg->father_id,
			   arg->level,
			   arg->verbose,
			   arg->fp);
  dsub_unsym2unsym_diag<T>(tmp);
  delete tmp;
  tmp = new C_Dsub_task<T>(arg->atomic_size,
			   arg->atomic_id,
			   arg->ir_bgn, 
			   arg->ir_end,
			   (-1), // jc_bgn
			   (-1), // jc_end
			   arg->ir_bgn_src2,
			   (-1), // ir_bgn_src2
			   (-1), // jc_bgn_src
			   (-1), // jc_bgn_src2
			   arg->dst_row,
			   arg->dst_mtrx,
			   (RectBlockMatrix<T> *)NULL, // dst_pt2
			   arg->ir_block,
			   arg->jc_block,
			   (RectBlockMatrix<T> *)NULL, // dst_pt2
			   (-1),
			   (-1),
			   arg->src_pt2,  
			   (SquareBlockMatrix<T>*)NULL, // src_pt2
			   dsub_unsym2unsym_diag,
			   false,  // dummy
			   *(arg->ops_complexity),
			   arg->father_id,
			   arg->level,
			   arg->verbose,
			   arg->fp);
  dsub_unsym2unsym_diag<T>(tmp);
  delete tmp;
}

template
void dsub_unsym2unsym_diag_two<double>(C_Dsub_task<double> *arg);

template
void dsub_unsym2unsym_diag_two<quadruple>(C_Dsub_task<quadruple> *arg);

template
void dsub_unsym2unsym_diag_two<complex<double> >(C_Dsub_task<complex<double> > *arg);

template
void dsub_unsym2unsym_diag_two<complex<quadruple> >(C_Dsub_task<complex<quadruple> > *arg);
// end of template function dsub_unsym2unsym_diag_two<T>

template<typename T>
void dsub_sym2sym_two(C_Dsub_task<T> *arg)
{
  if (arg->isSkip) {
    return;
  }    
  C_Dsub_task<T> *tmp;
  tmp = new C_Dsub_task<T>(arg->atomic_size,
			   arg->atomic_id,
			   arg->ir_bgn, 
			   arg->ir_end,
			   arg->jc_bgn,
			   arg->jc_end,
			   arg->ir_bgn_src,
			   (-1), // ir_bgn_src2
			   arg->jc_bgn_src,
			   (-1), // jc_bgn_src2
			   arg->dst_row,
			   arg->dst_mtrx,
			   (RectBlockMatrix<T> *)NULL, // dst_pt2
			   arg->ir_block,
			   arg->jc_block,
			   (RectBlockMatrix<T> *)NULL, // dst_pt2
			   (-1),
			   (-1),
			   arg->src_pt,  
			   (SquareBlockMatrix<T>*)NULL, // src_pt2
			   dsub_sym2sym,
			   false,  // dummy
			   *(arg->ops_complexity),
			   arg->father_id,
			   arg->level,
			   arg->verbose,
			   arg->fp);
  dsub_sym2sym<T>(tmp);
  delete tmp;
  tmp = new C_Dsub_task<T>(arg->atomic_size,
			   arg->atomic_id,
			   arg->ir_bgn, 
			   arg->ir_end,
			   arg->jc_bgn,
			   arg->jc_end,
			   arg->ir_bgn_src2,
			   (-1), // ir_bgn_src2
			   arg->jc_bgn_src2,
			   (-1), // jc_bgn_src2
			   arg->dst_row,
			   arg->dst_mtrx,
			   (RectBlockMatrix<T> *)NULL, // dst_pt2
			   arg->ir_block,
			   arg->jc_block,
			   (RectBlockMatrix<T> *)NULL, // dst_pt2
			   (-1),
			   (-1),
			   arg->src_pt2,  			
			   (SquareBlockMatrix<T>*)NULL, // src_pt2
			   dsub_sym2sym,
			   false,  // dummy
			   *(arg->ops_complexity),
			   arg->father_id,
			   arg->level,
			   arg->verbose,
			   arg->fp);
  dsub_sym2sym<T>(tmp);
  delete tmp;
}

template
void dsub_sym2sym_two<double>(C_Dsub_task<double> *arg);

template
void dsub_sym2sym_two<quadruple>(C_Dsub_task<quadruple> *arg);

template
void dsub_sym2sym_two<complex<double> >(C_Dsub_task<complex<double> > *arg);

template
void dsub_sym2sym_two<complex<quadruple> >(C_Dsub_task<complex<quadruple> > *arg);
// end of template function dsub_sym2sym_two<T>


template<typename T>
void dsub_unsym2unsym_two(C_Dsub_task<T> *arg)
{
  if (arg->isSkip) {
    return;
  }  
  C_Dsub_task<T> *tmp;
  tmp = new C_Dsub_task<T>(arg->atomic_size,
			   arg->atomic_id,
			   arg->ir_bgn, 
			   arg->ir_end,
			   arg->jc_bgn,
			   arg->jc_end,
			   arg->ir_bgn_src,
			   (-1), // ir_bgn_src2
			   arg->jc_bgn_src,
			   (-1), // jc_bgn_src2
			   arg->dst_row,
			   arg->dst_mtrx,
			   (RectBlockMatrix<T> *)NULL, // dst_pt2
			   arg->ir_block,
			   arg->jc_block,
			   (RectBlockMatrix<T> *)NULL, // dst_pt2
			   (-1),
			   (-1),
			   arg->src_pt,  
			   (SquareBlockMatrix<T>*)NULL, // src_pt2
			   dsub_unsym2unsym,
			   false,  // dummy
			   *(arg->ops_complexity),
			   arg->father_id,
			   arg->level,
			   arg->verbose,
			   arg->fp);
  dsub_unsym2unsym<T>(tmp);
  delete tmp;
  tmp = new C_Dsub_task<T>(arg->atomic_size,
			   arg->atomic_id,
			   arg->ir_bgn, 
			   arg->ir_end,
			   arg->jc_bgn,
			   arg->jc_end,
			   arg->ir_bgn_src2,
			   (-1), // ir_bgn_src2
			   arg->jc_bgn_src2,
			   (-1), // jc_bgn_src2
			   arg->dst_row,
			   arg->dst_mtrx,
			   (RectBlockMatrix<T> *)NULL, // dst_pt2
			   arg->ir_block,
			   arg->jc_block,
			   (RectBlockMatrix<T> *)NULL, // dst_pt2
			   (-1),
			   (-1),
			   arg->src_pt2,  
			   (SquareBlockMatrix<T>*)NULL, // src_pt2
			   dsub_unsym2unsym,
			   false,  // dummy
			   *(arg->ops_complexity),
			   arg->father_id,
			   arg->level,
			   arg->verbose,
			   arg->fp);
  dsub_unsym2unsym<T>(tmp);
  delete tmp;
}
template
void dsub_unsym2unsym_two<double>(C_Dsub_task<double> *arg);

template
void dsub_unsym2unsym_two<quadruple>(C_Dsub_task<quadruple> *arg);

template
void dsub_unsym2unsym_two<complex<double> >(C_Dsub_task<complex<double> > *arg);

template
void dsub_unsym2unsym_two<complex<quadruple> >(C_Dsub_task<complex<quadruple> > *arg);
// end of template function dsub_unsym2unsym_two<T>

template<typename T>
void dsub_unsym2diag_two(C_Dsub_task<T> *arg)
{
  if (arg->isSkip) {
    return;
  }  
  C_Dsub_task<T> *tmp;
  tmp = new C_Dsub_task<T>(arg->atomic_size,
			   arg->atomic_id,
			   arg->ir_bgn, 
			   arg->ir_end,
			   arg->jc_bgn,
			   arg->jc_end,
			   arg->ir_bgn_src,
			   (-1), // ir_bgn_src2
			   arg->jc_bgn_src,
			   (-1), // jc_bgn_src2
			   arg->dst_row,
			   arg->dst_mtrx,
			   (RectBlockMatrix<T> *)NULL, // dst_pt2
			   arg->ir_block,
			   arg->jc_block,
			   (RectBlockMatrix<T> *)NULL, // dst_pt2
			   (-1),
			   (-1),
			   arg->src_pt,  
			   (SquareBlockMatrix<T>*)NULL, // src_pt2
			   dsub_unsym2unsym,
			   false,  // dummy
			   *(arg->ops_complexity),
			   arg->father_id,
			   arg->level,
			   arg->verbose,
			   arg->fp);
  dsub_unsym2diag<T>(tmp);
  delete tmp;
  tmp = new C_Dsub_task<T>(arg->atomic_size,
			   arg->atomic_id,
			   arg->ir_bgn, 
			   arg->ir_end,
			   arg->jc_bgn,
			   arg->jc_end,
			   arg->ir_bgn_src2,
			   (-1), // ir_bgn_src2
			   arg->jc_bgn_src2,
			   (-1), // jc_bgn_src2
			   arg->dst_row,
			   arg->dst_mtrx,
			   (RectBlockMatrix<T> *)NULL, // dst_pt2
			   arg->ir_block,
			   arg->jc_block,
			   (RectBlockMatrix<T> *)NULL, // dst_pt2
			   (-1),
			   (-1),
			   arg->src_pt2,  
			   (SquareBlockMatrix<T>*)NULL, // src_pt2
			   dsub_unsym2unsym,
			   false, // dummy
			   *(arg->ops_complexity),
			   arg->father_id,
			   arg->level,
			   arg->verbose,
			   arg->fp);
  dsub_unsym2diag<T>(tmp);
  delete tmp;
}

template
void dsub_unsym2diag_two<double>(C_Dsub_task<double> *arg);

template
void dsub_unsym2diag_two<quadruple>(C_Dsub_task<quadruple> *arg);

template
void dsub_unsym2diag_two<complex<double> >(C_Dsub_task<complex<double> > *arg);


template
void dsub_unsym2diag_two<complex<quadruple> >(C_Dsub_task<complex<quadruple> > *arg);
// end of template function dsub_unsym2diag_two<T>

template<typename T>
void dsub_sym2rct_two(C_Dsub_task<T> *arg) 
{
  C_Dsub_task<T> *tmp;
  tmp = new C_Dsub_task<T>(arg->atomic_size,
			   arg->atomic_id,
			   arg->ir_bgn, 
			   arg->ir_end,
			   arg->jc_bgn,
			   arg->jc_end,
			   arg->ir_bgn_src,
			   (-1), // ir_bgn_src2
			   arg->jc_bgn_src,
			   (-1), // jc_bgn_src2
			   arg->dst_row,
			   (SquareBlockMatrix<T>*)NULL, // src_pt2
			   arg->dst_pt,
			   arg->ir_block, // 0,
			   arg->jc_block, // 0,
			   (RectBlockMatrix<T> *)NULL, // dst_pt2
			   (-1),
			   (-1),
			   arg->src_pt,  
			   (SquareBlockMatrix<T>*)NULL, // src_pt2
			   dsub_sym2rct,
			   false,  // dummy
			   *(arg->ops_complexity),
			   arg->father_id,
			   arg->level,
			   arg->verbose,
			   arg->fp);
  //  tmp->verbose = false;
  dsub_sym2rct<T>(tmp);
  delete tmp;
  tmp = new C_Dsub_task<T>(arg->atomic_size,
			   arg->atomic_id,
			   arg->ir_bgn, 
			   arg->ir_end,
			   arg->jc_bgn,
			   arg->jc_end,
			   arg->ir_bgn_src2,
			   (-1), // ir_bgn_src2
			   arg->jc_bgn_src2,
			   (-1), // jc_bgn_src2
			   arg->dst_row,
			   (SquareBlockMatrix<T>*)NULL, // src_pt2
			   arg->dst_pt,
			   arg->ir_block, //0,
			   arg->jc_block, //0,
			   (RectBlockMatrix<T> *)NULL, // dst_pt2
			   (-1),
			   (-1),
			   arg->src_pt2,  
			   (SquareBlockMatrix<T>*)NULL, // src_pt2
			   dsub_sym2rct,
			   false,  // dummy
			   *(arg->ops_complexity),
			   arg->father_id,
			   arg->level,
			   arg->verbose,
			   arg->fp);
  //  tmp->verbose = false;
  dsub_sym2rct<T>(tmp);
  delete tmp;
}

template
void dsub_sym2rct_two<double>(C_Dsub_task<double> *arg);

template
void dsub_sym2rct_two<quadruple>(C_Dsub_task<quadruple> *arg);
template
void dsub_sym2rct_two<complex<double> >(C_Dsub_task<complex<double> > *arg);

template
void dsub_sym2rct_two<complex<quadruple> >(C_Dsub_task<complex<quadruple> > *arg);
// end of template function dsub_sym2rct_two<T>

template<typename T>
void dsub_unsym2rct_two(C_Dsub_task<T> *arg) 
{
  C_Dsub_task<T> *tmp;
  tmp = new C_Dsub_task<T>(arg->atomic_size,
			   arg->atomic_id,
			   arg->ir_bgn, 
			   arg->ir_end,
			   arg->jc_bgn,
			   arg->jc_end,
			   arg->ir_bgn_src,
			   (-1), // ir_bgn_src2
			   arg->jc_bgn_src,
			   (-1), // jc_bgn_src2
			   arg->dst_row,
			   (SquareBlockMatrix<T>*)NULL, // src_pt2
			   arg->dst_pt,
			   arg->ir_block,//0,
			   arg->jc_block, //0,
			   arg->dst_pt2,
			   0,
			   0,
			   arg->src_pt,  
			   (SquareBlockMatrix<T>*)NULL, // src_pt2
			   dsub_unsym2rct,
			   false,  // dummy
			   *(arg->ops_complexity),
			   arg->father_id,
			   arg->level,
			   arg->verbose,
			   arg->fp);
  dsub_unsym2rct<T>(tmp);
  delete tmp;
  tmp = new C_Dsub_task<T>(arg->atomic_size,
			   arg->atomic_id,
			   arg->ir_bgn, 
			   arg->ir_end,
			   arg->jc_bgn,
			   arg->jc_end,
			   arg->ir_bgn_src2,
			   (-1), // ir_bgn_src2
			   arg->jc_bgn_src2,
			   (-1), // jc_bgn_src2
			   arg->dst_row,
			   (SquareBlockMatrix<T>*)NULL, // src_pt2
			   arg->dst_pt,
			   arg->ir_block, //0,
			   arg->jc_block,//0,
			   arg->dst_pt2,
			   0,
			   0,
			   arg->src_pt2,  
			   (SquareBlockMatrix<T>*)NULL, // src_pt2
			   dsub_unsym2rct,
			   false,  // dummy
			   *(arg->ops_complexity),
			   arg->father_id,
			   arg->level,
			   arg->verbose,
			   arg->fp);
  dsub_unsym2rct<T>(tmp);
  delete tmp;
}

template
void dsub_unsym2rct_two<double>(C_Dsub_task<double> *arg);

template
void dsub_unsym2rct_two<quadruple>(C_Dsub_task<quadruple> *arg);

template
void dsub_unsym2rct_two<complex<double> >(C_Dsub_task<complex<double> > *arg);

template
void dsub_unsym2rct_two<complex<quadruple> >(C_Dsub_task<complex<quadruple> > *arg);
// end of template function dsub_unsym2rct_two<T>

//#define DEBUG_QUEUE_GENERATION2
//#define DEBUG_QUEUE_GENERATION

//#define  DEBUG_PREPARE_THREAD
//#define  DEBUG_PREPARE_THREAD_DEBUG

template<typename T>
void C_Dsub_queue(bool isSym, 
		  int father_id,
		  bool skip_flag,
		  vector<C_task *>& queue,
		  list <child_contribution<T> > &child_contrib,
		  vector<C_task *>* tasks_p, // _tasks_DSymmGEMM
		  vector<int>* tasks_p_indcol,
		  const bool tasks_p_flag,
		  vector<C_task *>* tasks_q, // _tasks_DfillSymm
		  vector<C_task *>* tasks_r, // _tasks_SparseLocalSchur
		  vector<C_task *>* tasks_s, // _tasks_DSub[level + 1][(*it)]
		  vector<C_task *>* tasks_d, // _tasks_deallocateLocalSchur
		  vector<int>* tasks_d_indcol,
		  int level,
		  const bool verbose,
		  FILE *fp)
{
  //    list <child_contribution>& child_contrib = child_contribs[*it];
  const int diag_size =    child_contrib.front().diag_size;
  const int father_row =   child_contrib.front().father_row;

  const int size_res = diag_size % SIZE_B1;
  //  const int size_res2 = offdiag_size % SIZE_B1;
  int num_block;
  num_block = diag_size / SIZE_B1 + (size_res != 0);
  //  const int num_block2 = (offdiag_size / SIZE_B1 + (size_res2 != 0));
  const int num_block2 = child_contrib.front().father_offdiag_pt->num_blocks_c();
  const int block_diag_size = (num_block * (num_block + 1)) / 2;
  const int block_offdiag_size = num_block * num_block2;

  //  vector<C_task *>* tasks_p = ((tasks_p_ == NULL) ? NULL :
  //			       (tasks_p_->size() > 0 ? tasks_p_ : NULL));
  if (diag_size == 0) {
    queue.resize(1);
    int nb = father_id + 1;
    string task_name = ("i dummy : " +
			to_string(level) + " : " +to_string(nb));
    C_dummy_arg *arg = new C_dummy_arg(verbose, &fp, nb);
    //    *(arg->ops_complexity) = (-1L);
    queue[0] = new C_task(C_DUMMY,
			  task_name,
			  (void *)arg,
			  C_dummy,
			  1, // atomic_size,
			  0, // atomic_id,
			  arg->ops_complexity);
    queue[0]->parents->clear();
    fprintf(fp, "%s %d : %s\n", __FILE__, __LINE__, task_name.c_str());
    return;
  }
#ifdef DEBUG_PREPARE_THREAD
  cout << "father = " << (father_id + 1) // selfIndex^-1
       << " # of child = " << child_contrib.size() << " [ ";
#endif      
  for (typename list<child_contribution<T> >::const_iterator jt = child_contrib.begin();
       jt != child_contrib.end(); ++jt) {
#ifdef DEBUG_PREPARE_THREAD
    cout << (*jt).child_id << " ";
#endif
  }
#ifdef DEBUG_PREPARE_THREAD
  cout << "]" << endl;
#endif
#ifdef DEBUG_STRIPES
  for (typename list<child_contribution<T> >::const_iterator jt = child_contrib.begin();
       jt != child_contrib.end(); ++jt) {
    cout << "child = " << (*jt).child_id << endl;
    cout << "diag = ";
    for (list <index_strip>::const_iterator kt = (*jt).diag_strip.begin();
	 kt != (*jt).diag_strip.end();
	 ++kt) {
      cout << "[ " << (*kt).begin_dst << " , " 
	   << (*kt).begin_src << " , "
	   << (*kt).width << " ] ";
    }
    cout << endl;
    cout << "offdiag = ";
    for (list <index_strip>::const_iterator kt = (*jt).offdiag_strip.begin();
	 kt != (*jt).offdiag_strip.end();
	 ++kt) {
      cout << "[ " << (*kt).begin_dst << " , " 
	   << (*kt).begin_src << " , "
	   << (*kt).width << " ] ";
    }
    cout << endl;
  } // loop: jt
#endif
  //
  if (child_contrib.size() == 2) {
    //      vector<C_task *>&queue = tasks_C_DSUB[*it];
    //     queue.resize(block_diag_size + num_block2);
    //      vector<list<C_Dsub_task *>*>&C_task_arg = lists_C_DSUB[*it];
    //      C_task_arg.resize(block_diag_size + num_block2);
#ifdef  DEBUG_PREPARE_THREAD_DEBUG
    cout << "2 children resize = " << block_diag_size + num_block2 << " ";
#endif
    child_contribution<T> child0 = child_contrib.front();
    child_contribution<T> child1 = child_contrib.back();
#if 0
    for (typename list<child_contribution<T> >::const_iterator jt = child_contrib.begin();
	 jt != child_contrib.end(); ++jt) {
      cout << "child = " << (*jt).child_id 
	   << " offdiag = ";
      for (list <index_strip>::const_iterator kt = (*jt).offdiag_strip.begin();
	   kt != (*jt).offdiag_strip.end();
	   ++kt) {
	cout << "[ " << (*kt).begin_dst << " , " 
	     << (*kt).begin_src << " , "
	     << (*kt).width << " ] ";
      }
      cout << endl;
    }
#endif
    if ((child0.child_pt->dimension() == 0) ||
	(child1.child_pt->dimension() == 0)) {
      fprintf(fp, "%s %d : %d %d %d : %d %d %d\n", __FILE__, __LINE__,
	      child0.child_id,
	      (int)child0.diag_strip.size(),
	      (int)child0.offdiag_strip.size(),
	      child1.child_id,
	      (int)child1.diag_strip.size(),
	      (int)child1.offdiag_strip.size());
    }

    list <index_strip> *strips_r, *strips_c;
    strips_r = new list <index_strip>[2];
    strips_c = new list <index_strip>[2];
    list <index_strip2> strips_r01, strips_c01;
    int *child_id = new int[2];
    child_id[0] = child0.child_id;
    child_id[1] = child1.child_id;
    //    double ***child_pt = new double**[2];
    SquareBlockMatrix<T>** child_pt = new SquareBlockMatrix<T>*[2];
    child_pt[0] = child0.child_pt;
    child_pt[1] = child1.child_pt;
    SquareBlockMatrix<T> *father_diag_pt = child0.father_diag_pt;
    RectBlockMatrix<T> *father_offdiag_pt = child0.father_offdiag_pt;
    RectBlockMatrix<T> *father_offdia2_pt = child0.father_offdiag_unsym_pt;

    if(verbose) {
      fprintf(fp, "%s %d : father id = %d\n",
	      __FILE__, __LINE__, (father_id + 1));
    }
    if (child0.diag_strip.size() > 1 || child1.diag_strip.size() > 1) {
#ifdef DEBUG_PREPARE_THREAD
      fprintf(fp, "%s %d : row block is not continous\n",
	      __FILE__, __LINE__);
#endif
      combine_two_strips(strips_r[0], strips_r[1], strips_r01, 
			 child0.diag_strip,
			 child1.diag_strip,
			 child0.diag_size);
    }
    else {
      if ((child0.diag_strip.size() == 0) &&
	  (child1.diag_strip.size() == 0)) {
	if (verbose) {
	  fprintf(fp, "%s %d : both rows are null : father = %d\n",
		  __FILE__, __LINE__, (father_id + 1));
	}
	strips_r[0].clear();
	strips_r[1].clear();
	strips_r01.clear();
      }
      if ((child0.diag_strip.size() == 0) || 
	  (child1.diag_strip.size() == 0)) {
	if (verbose) {
	  fprintf(fp, "%s %d : one of rows is null / ",
		  __FILE__, __LINE__);
	}
	if (child0.diag_strip.size() == 0) {
	  if (verbose) {
	    fprintf(fp, " r0 null / ");
	  }
	  strips_r[0].clear();
	  strips_r[1] = child1.diag_strip;
	}
	else {
	  if (verbose) {
	    fprintf(fp, "%d : %d / ",  
		    child0.diag_strip.front().begin_dst,
		    child0.diag_strip.front().width);
	  }
	}
	if (child1.diag_strip.size() == 0) {
	  if (verbose) {
	    fprintf(fp, " r1 null\n");
	  }
	  strips_r[0] = child0.diag_strip;
	  strips_r[1].clear();
	}
	else {
	  if (verbose) {
	    fprintf(fp, "%d : %d\n",  
		    child1.diag_strip.front().begin_dst,
		    child1.diag_strip.front().width);
	  }
	} // if (verbose)
	strips_r01.clear();
      } //  (child0.diag_strip.size() == 0) || (child1.diag_strip.size() == 0)
      else {
	if ((child0.diag_strip.front().begin_dst != 
	     child1.diag_strip.front().begin_dst) ||
	    (child0.diag_strip.front().width != 
	     child1.diag_strip.front().width)) {
	  if (verbose) {
	    fprintf(fp, 
		    "%s %d : row blocks are not same / ",
		    __FILE__, __LINE__);
	    fprintf(fp, "%d : %d : %d / %d : %d : %d\n",
		    (int)child0.diag_strip.size(),
		    child0.diag_strip.front().begin_dst,
		    child0.diag_strip.front().width,
		    (int)child1.diag_strip.size(),
		    child1.diag_strip.front().begin_dst,
		    child1.diag_strip.front().width);
	  } // if (verbose)
	  split_two_strips(strips_r[0], strips_r[1], strips_r01, 
			   child0.diag_strip.front(),
			   child1.diag_strip.front());
	  if (verbose) {
	    fprintf(fp, "%s %d : split of strips r0 ",
		    __FILE__, __LINE__);
	    for (list<index_strip>::const_iterator kt = strips_r[0].begin(); 
		 kt != strips_r[0].end(); ++kt) {
	      fprintf(fp, "%d : %d : %d ", 
		      (*kt).begin_dst, (*kt).begin_src, (*kt).width);
	    }
	    fprintf(fp, " r1 ");
	    for (list<index_strip>::const_iterator kt = strips_r[1].begin(); 
		 kt != strips_r[1].end(); ++kt) {
	      fprintf(fp, "%d : %d : %d ", 
		      (*kt).begin_dst, (*kt).begin_src, (*kt).width);
	    }
	    fprintf(fp, "r01 ");
	    for (list<index_strip2>::const_iterator kt = strips_r01.begin(); 
		 kt != strips_r01.end(); ++kt) {
	      fprintf(fp, "%d : %d : %d : %d ", 
		    (*kt).begin_dst, (*kt).begin_src0, (*kt).begin_src1,
		      (*kt).width);
	    }
	    fprintf(fp, "\n");
	  } // if (verbose)
	}
	else {
	  strips_r[0].clear();
	  strips_r[1].clear();
	  copy_two_strips(strips_r01, child0.diag_strip, child1.diag_strip);
	  if (verbose) {
	    fprintf(fp, "%s %d : r01 ", __FILE__, __LINE__);
	    for (list<index_strip2>::const_iterator kt = strips_r01.begin(); 
		 kt != strips_r01.end(); ++kt) {
	      fprintf(fp, "%d : %d : %d : %d ", 
		      (*kt).begin_dst, (*kt).begin_src0,
		      (*kt).begin_src1, (*kt).width);
	    }
	    fprintf(fp, "\n");
	  } // if (verbose)
	}
      }//  (child0.diag_strip.size() == 0) || (child1.diag_strip.size() == 0)
    }  //  (child0.diag_strip.size() > 1 || child1.diag_strip.size() > 1)
    if ((child0.offdiag_strip.size() > 0) &&
	(child1.offdiag_strip.size() > 0)) {
      combine_two_strips(strips_c[0], strips_c[1], strips_c01, 
			 child0.offdiag_strip,
			 child1.offdiag_strip,
			 child0.offdiag_size);
    }
    else if (child1.offdiag_strip.size() == 0) {
      strips_c[0]= child0.offdiag_strip;
      strips_c[1].clear();
      strips_c01.clear();
    }
    else { // (child0.offdiag_strip.size() == 0)
      strips_c[1]= child1.offdiag_strip;
      strips_c[0].clear();
      strips_c01.clear();
    }
    if ((child0.child_pt->dimension() == 0) ||
	(child1.child_pt->dimension() == 0)) {
      fprintf(fp, "%s %d : %d %d %d\n", __FILE__, __LINE__,
	      (int)strips_c[0].size(), (int)strips_c[1].size(),
	      (int)strips_c01.size());
    }
    queue.resize(block_diag_size + num_block2);
    //    queue.resize(block_diag_size + block_offdiag_size);
    if (verbose) {
#if 0
      print_strips("child0.diag_strip", child0.diag_strip, fp);
      print_strips("child1.diag_strip", child1.diag_strip, fp);
      print_strips("strips_r0", strips_r[0], fp);
      print_strips("strips_r1", strips_r[1], fp);
      print_strips("strips_r01", strips_r01, fp);
      print_strips("child0.offdiag_strip", child0.offdiag_strip, fp);
      print_strips("child1.offdiag_strip", child1.offdiag_strip, fp);
      print_strips("strips_c0", strips_c[0], fp);
      print_strips("strips_c1", strips_c[1], fp);
      print_strips("strips_c01", strips_c01, fp);
#endif
      fprintf(fp, "%s %d : C_Dsub with father id = %d %d %d ", 
	      __FILE__, __LINE__, (father_id + 1), num_block, num_block2);
      if (skip_flag) {
	fprintf(fp, "skipped\n");
      }
      else {
	fprintf(fp, "\n");
      }
    } // if (verbose)
    int kk = 0;
    for (int kc0 = 0; kc0 < num_block; kc0++) {
      const int kc = kc0 * SIZE_B1;
      const int kc_end = (kc0 == num_block - 1) ? diag_size : (kc + SIZE_B1);
      for (int kr0 = 0; kr0 <= kc0; kr0++) { // running over upper blocks
	const int kr = kr0 * SIZE_B1;
	const int kr_end = 
	  (kr0 == num_block - 1) ? diag_size : (kr + SIZE_B1);

	list<C_Dsub_task<T> *> *C_task_arg = new list<C_Dsub_task<T> *>;
	//	  C_task_arg[kk] = new list<C_Dsub_task *>;
	long *ops_sum = new long;
	*ops_sum = 0L;
	list<int>* parents_r = new list<int>[2];
	list<int>* parents_c = new list<int>[2];

	// (child0 + child01 + child1) * (child0 + child01 + child1)
	// = child0  * child0  + child0  * child01
	// + child01 * child0  + child01 * child01 + child01 * child1
	//                     + child1  * child01 + child1  * child1
	// 
	// diagonal contribution
	// tasks are stored in each queue _c_dsub_arg[*it][ ]
	// ---- child_ll * (child_ll + child01) ---- ll = 0, 1
	for (int ll = 0; ll < 2; ll++) { 
	  for (list <index_strip>::const_iterator mt = strips_r[ll].begin();
	       mt != strips_r[ll].end(); ++mt) {
	    const int r_bgn_dst = (*mt).begin_dst;
	    const int r_end_dst = r_bgn_dst + (*mt).width;
	    if (r_end_dst < kr) { 
	      continue;   // continue loop : mt : strips_r[ll]
	    }
	    if (r_bgn_dst >= kr_end) { 
	      break;      // break loop : mt : strips_r[ll]
	    }
	    const int ir_bgn = imax(kr, r_bgn_dst); 
	    const int ir_end = imin(kr_end, r_end_dst); 
	    const int ic_bgn = imax(kc, r_bgn_dst);    
	    const int ic_end = imin(kc_end, r_end_dst);
	    // -- child_ll * child_ll
	    if (kr == kc) {		
	      const long ops = ((long)(ir_end - ir_bgn) *
				((long)(ir_end - ir_bgn) + 1L) / 2L);
	      // diagonal : half size
	      const int ir_bgn_src = 
		((*mt).begin_src + (ir_bgn - r_bgn_dst));
	      const int ir_end_src = ir_bgn_src + (ir_end - ir_bgn);
	      if (ir_bgn < ir_end) {
		C_Dsub_task<T> *tmp = 
		  new C_Dsub_task<T>(1,
				     0,
				     (ir_bgn - kr), 
				     (ir_end - kr), 
				     (-1), // jc_bgn
				     (-1), // jc_end
				     ir_bgn_src,
				     (-1), // ir_bgn_src2
				     (-1), // jc_bgn_src
				     (-1), // jc_bgn_src2
				     father_row, 
				     father_diag_pt, 
				     (RectBlockMatrix<T> *)NULL, // dst_pt
				     kr0,
				     kr0,
				     (RectBlockMatrix<T> *)NULL, // dist_pt2
				     (-1),
				     (-1),
				     child_pt[ll],  
				     (SquareBlockMatrix<T>*)NULL, // src_pt2
				     (isSym ? dsub_sym2sym_diag<T> : 
				      dsub_unsym2unsym_diag<T>), 
				     false,  // skip_flag
				     ops,
				     father_id,
				     level,
				     verbose,
				     fp);
		tmp->child0_id = child_id[0];
		tmp->child1_id = child_id[1];
		C_task_arg->push_back(tmp);
		*ops_sum += ops;
		update_parents_list(parents_r[ll], 
				    ir_bgn_src, ir_end_src, child_pt[ll]);
	      } // if (ir_bgn < ir_end)
	    }      // if (kr == kc)
	    else { // kr < kc
	      const long ops = ((long)(ir_end - ir_bgn) *
				(long)(ic_end - ic_bgn));
	      const int ir_bgn_src = 
		((*mt).begin_src + (ir_bgn - r_bgn_dst));
	      const int ic_bgn_src = 
		((*mt).begin_src + (ic_bgn - r_bgn_dst));
	      const int ir_end_src = ir_bgn_src + (ir_end - ir_bgn);
	      const int ic_end_src = ic_bgn_src + (ic_end - ic_bgn);
	      if ((ir_bgn < ir_end) && (ic_bgn < ic_end)) {
		C_Dsub_task<T> *tmp = 
		  new C_Dsub_task<T>(1,
				     0,
				     (ir_bgn - kr), 
				     (ir_end - kr), 
				     (ic_bgn - kc), 
				     (ic_end - kc), 
				     ir_bgn_src,
				     (-1),  // ir_bgn_src2
				     ic_bgn_src,
				     (-1),  // jc_bgn_src2
				     father_row, 
				     father_diag_pt, 
				     (RectBlockMatrix<T> *)NULL, // dst_pt
				     kr0,
				     kc0,
				     (RectBlockMatrix<T> *)NULL,  // dist_pt2
				     (-1),
				     (-1),
				     child_pt[ll],
				     (SquareBlockMatrix<T>*)NULL, // src_pt2
				     isSym ? dsub_sym2sym<T> : dsub_unsym2unsym<T>,
				     false, // skip_flag
				     ops,
				     father_id,
				     level,
				     verbose,
				     fp);
		C_task_arg->push_back(tmp);
		*ops_sum += ops;
		  
		update_parents_list(parents_r[ll], 
				    ir_bgn_src, ir_end_src, child_pt[ll]);
		update_parents_list(parents_c[ll], 
				    ic_bgn_src, ic_end_src, child_pt[ll]);
#if 0
		if (father_id == 0 && ll == 1) {
		  fprintf(stderr, "%s : %d [%d %d] %d\n",
			  __FILE__, __LINE__,
			  ic_bgn_src, ic_end_src, parents_c[ll].back());
		}
#endif
		//SIZE_B1);
	      } // if ((ir_bgn < ir_end) && (ic_bgn < ic_end)) 
	    } // if (kc == kr)
	      // off-diagonal blocks
	    list <index_strip>::const_iterator nt = mt;
	    ++nt;
	    for (; nt != strips_r[ll].end(); ++nt) {
	      const int c_bgn_dst = (*nt).begin_dst;
	      const int c_end_dst = c_bgn_dst + (*nt).width;
	      if (c_end_dst < kc) {
		continue;
	      }
	      if (c_bgn_dst >= kc_end) {
		break;
	      }
	      const int ir_bgn = imax(kr, r_bgn_dst); 
	      const int ir_end = imin(kr_end, r_end_dst); 
	      const int jc_bgn = imax(kc, c_bgn_dst); 
	      const int jc_end = imin(kc_end, c_end_dst); 
	      const long ops = ((long)(ir_end - ir_bgn) * 
				(long)(jc_end - jc_bgn));

	      const int ir_bgn_src = ((*mt).begin_src + (ir_bgn - r_bgn_dst));
	      const int jc_bgn_src = ((*nt).begin_src + (jc_bgn - c_bgn_dst));
	      const int ir_end_src = ir_bgn_src + (ir_end - ir_bgn);
	      const int jc_end_src = jc_bgn_src + (jc_end - jc_bgn);
	      if ((ir_bgn < ir_end) && (jc_bgn < jc_end)) {
		C_Dsub_task<T> *tmp = 
		  new C_Dsub_task<T>(1,
				     0,
				     (ir_bgn - kr), 
				     (ir_end - kr), 
				     (jc_bgn - kc), 
				     (jc_end - kc), 
				     ir_bgn_src,
				     (-1), // ir_bgn_src2
				     jc_bgn_src,
				     (-1), // jc_bgn_src2
				     father_row, 
				     father_diag_pt, 
				     (RectBlockMatrix<T> *)NULL, // dst_pt
				     kr0,
				     kc0,
				     (RectBlockMatrix<T> *)NULL, // dst_pt2
				     (-1),
				     (-1),
				     child_pt[ll],
				     (SquareBlockMatrix<T>*)NULL, // src_pt2
				     isSym ? dsub_sym2sym<T> : ((kr0 == kc0) ? dsub_unsym2diag<T> :dsub_unsym2unsym<T>),
				     false, // skip_flag
				     ops,
				     father_id,
				     level,
				     verbose,
				     fp);
		C_task_arg->push_back(tmp);
		*ops_sum += ops;
		  
		update_parents_list(parents_r[ll], 
				    ir_bgn_src, ir_end_src, child_pt[ll]);
		update_parents_list(parents_c[ll], 
				    jc_bgn_src, jc_end_src, child_pt[ll]);
#if 0
		if (father_id == 0 && ll == 1) {
		  fprintf(stderr, "%s : %d [%d %d] %d\n",
			  __FILE__, __LINE__,
			  jc_bgn_src, jc_end_src, parents_c[ll].back());
		}
#endif
	      } // if ((ir_bgn < ir_end) && (jc_bgn < jc_end))
	    } // loop : nt
	    // -- child_ll * child01
	    for (list <index_strip2>::const_iterator nt = strips_r01.begin();
		 nt != strips_r01.end(); ++nt) {
	      const int c_bgn_dst = (*nt).begin_dst;
	      const int c_end_dst = c_bgn_dst + (*nt).width;
	      if (c_bgn_dst < r_bgn_dst) {
		continue;
	      }
	      if (c_end_dst < kc) {
		continue;
	      }
	      if (c_bgn_dst >= kc_end) {
		break;
	      }
	      const int ir_bgn = imax(kr, r_bgn_dst); 
	      const int ir_end = imin(kr_end, r_end_dst);
	      const int jc_bgn = imax(kc, c_bgn_dst); 
	      const int jc_end = imin(kc_end, c_end_dst); 
	      const long ops = ((long)(ir_end - ir_bgn) *
				(long)(jc_end - jc_bgn));
	      const int ir_bgn_src = 
		((*mt).begin_src + (ir_bgn - r_bgn_dst));
	      const int jc_bgn_src = 
		((ll == 0 ? (*nt).begin_src0 : (*nt).begin_src1) +
		 (jc_bgn - c_bgn_dst));
	      const int ir_end_src = ir_bgn_src + (ir_end - ir_bgn);
	      const int jc_end_src = jc_bgn_src + (jc_end - jc_bgn);
	      if ((ir_bgn < ir_end) && (jc_bgn < jc_end)) {
		C_Dsub_task<T> *tmp = 
		  new C_Dsub_task<T>(1,
				     0,
				     (ir_bgn - kr), 
				     (ir_end - kr), 
				     (jc_bgn - kc), 
				     (jc_end - kc), 
				     ir_bgn_src,
				     (-1),  // ir_bgn_src2
				     jc_bgn_src, // (*nt).begin_src[ll] + ... 
				     (-1),  // jc_bgn_src2
				     father_row, 
				     father_diag_pt, 
				     (RectBlockMatrix<T> *)NULL, // dst_pt
				     kr0,
				     kc0,
				     (RectBlockMatrix<T> *)NULL,  // dst_pt2
				     (-1),
				     (-1),
				     child_pt[ll],
				     (SquareBlockMatrix<T>*)NULL, // src_pt2
				     isSym ? dsub_sym2sym<T> : ((kr0 == kc0) ? dsub_unsym2diag<T> :dsub_unsym2unsym<T>),
			     // dsub_unsym2unsym<T>, bug : 17 Jan.2016 found
				     false, // skip_flag
				     ops,
				     father_id,
				     level,
				     verbose,
				     fp);
		C_task_arg->push_back(tmp);
		*ops_sum += ops;
		  
		update_parents_list(parents_r[ll], 
				    ir_bgn_src, ir_end_src, child_pt[ll]);
		update_parents_list(parents_c[ll], 
				    jc_bgn_src, jc_end_src, child_pt[ll]);
#if 0
		if (father_id == 0 && ll == 1) {
		  fprintf(stderr, "%s : %d [%d %d] %d\n",
			  __FILE__, __LINE__,
			  jc_bgn_src, jc_end_src, parents_c[ll].back());
		}
#endif
	      } //  if ((ir_bgn < ir_end) && (jc_bgn < jc_end))v
	    } // loop : nt
	  }   // loop : mt
	}     // loop : ll
	// ---- child01 * (child0 + child01 + child1) ----
	// -- child01 * child01 
	for (list <index_strip2>::const_iterator mt = strips_r01.begin();
	     mt != strips_r01.end(); ++mt) {
	  const int r_bgn_dst = (*mt).begin_dst;
	  const int r_end_dst = r_bgn_dst + (*mt).width;
	  if (r_end_dst < kr) {
	    continue;
	  }
	  if (r_bgn_dst >= kr_end) {
	    break;
	  }
	  const int ir_bgn = imax(kr, r_bgn_dst); 
	  const int ir_end = imin(kr_end, r_end_dst); 
	  const int ic_bgn = imax(kc, r_bgn_dst);
	  const int ic_end = imin(kc_end, r_end_dst); 
	  if (kr == kc) {
	    const long ops = ((long)(ir_end - ir_bgn) *
			      ((long)(ir_end - ir_bgn) + 1L));
	    // diagonal : half size
	    int ir_bgn_srcs[2], ir_end_srcs[2];
	    ir_bgn_srcs[0] = (*mt).begin_src0 + (ir_bgn - r_bgn_dst);
	    ir_bgn_srcs[1] = (*mt).begin_src1 + (ir_bgn - r_bgn_dst);
	    ir_end_srcs[0] = ir_bgn_srcs[0] + (ir_end - ir_bgn);
	    ir_end_srcs[1] = ir_bgn_srcs[1] + (ir_end - ir_bgn);
	    if (ir_bgn < ir_end) {
	      C_Dsub_task<T> *tmp = 
		new C_Dsub_task<T>(1,
				   0,
				   (ir_bgn - kr), 
				   (ir_end - kr), 
				   (-1), // jc_bgn
				   (-1), // jc_end
				   ir_bgn_srcs[0],
				   ir_bgn_srcs[1],
				   (-1), // jc_bgn_src
				   (-1), // jc_bgn_src2
				   father_row, 
				   father_diag_pt, 
				   (RectBlockMatrix<T> *)NULL, // dst_pt
				   kr0,
				   kr0,
				   (RectBlockMatrix<T> *)NULL,  // dst_pt2
				   (-1),
				   (-1),
				   child_pt[0],  
				   child_pt[1],  
				   (isSym ? dsub_sym2sym_diag_two<T> : 
				    dsub_unsym2unsym_diag_two<T>), 
				   skip_flag,
				   ops,
				   father_id,
				   level,
				   verbose,
				   fp);
	      C_task_arg->push_back(tmp);
	      *ops_sum += skip_flag ? 0L : ops;
	      for (int ll = 0; ll < 2; ll++) {
		update_parents_list(parents_r[ll], 
				    ir_bgn_srcs[ll], ir_end_srcs[ll], 
				    child_pt[ll]);
	      }
	    } // if (ir_bgn < ir_end)
	  }  // kr != kc
	  else {
	    const long ops = ((long)(ir_end - ir_bgn) *
			      (long)(ic_end - ic_bgn) * 2L);
	    int ir_bgn_srcs[2], ir_end_srcs[2];
	    int ic_bgn_srcs[2], ic_end_srcs[2];
	    ir_bgn_srcs[0] = (*mt).begin_src0 + (ir_bgn - r_bgn_dst);
	    ir_bgn_srcs[1] = (*mt).begin_src1 + (ir_bgn - r_bgn_dst);
	    ic_bgn_srcs[0] = (*mt).begin_src0 + (ic_bgn - r_bgn_dst);
	    ic_bgn_srcs[1] = (*mt).begin_src1 + (ic_bgn - r_bgn_dst);
	      
	    ir_end_srcs[0] = ir_bgn_srcs[0] + (ir_end - ir_bgn);
	    ir_end_srcs[1] = ir_bgn_srcs[1] + (ir_end - ir_bgn);
	    ic_end_srcs[0] = ic_bgn_srcs[0] + (ic_end - ic_bgn);
	    ic_end_srcs[1] = ic_bgn_srcs[1] + (ic_end - ic_bgn);
	    if ((ir_bgn < ir_end) && (ic_bgn < ic_end)) {
	      C_Dsub_task<T> *tmp = 
		new C_Dsub_task<T>(1,
				   0,
				   (ir_bgn - kr), 
				   (ir_end - kr), 
				   (ic_bgn - kc), 
				   (ic_end - kc), 
				   ir_bgn_srcs[0],
				   ir_bgn_srcs[1],
				   ic_bgn_srcs[0],
				   ic_bgn_srcs[1],
				   father_row, 
				   father_diag_pt, 
				   (RectBlockMatrix<T> *)NULL, // dst_pt
				   kr0,
				   kc0,
				   (RectBlockMatrix<T> *)NULL,  // dst_pt2
				   (-1),
				   (-1),
				   child_pt[0],
				   child_pt[1],
				   isSym ? dsub_sym2sym_two<T> : dsub_unsym2unsym_two<T>,
				   skip_flag,
				   ops,
				   father_id,
				   level,
				   verbose,
				   fp);
	      C_task_arg->push_back(tmp);
	      *ops_sum += skip_flag ? 0L : ops;
	      for (int ll = 0; ll < 2; ll++) {
		update_parents_list(parents_r[ll], 
				    ir_bgn_srcs[ll], ir_end_srcs[ll], 
				    child_pt[ll]);
		update_parents_list(parents_c[ll], 
				    ic_bgn_srcs[ll], ic_end_srcs[ll], 
				    child_pt[ll]);
	      }
	    } // if ((ir_bgn < ir_end) && (ic_bgn < ic_end)) 
	  } // if (kr == kc)
	  // off-diagonal blocks
	  list <index_strip2>::const_iterator nt = mt;
	  ++nt;
	  for (; nt != strips_r01.end(); ++nt) {
	    const int c_bgn_dst = (*nt).begin_dst;
	    const int c_end_dst = c_bgn_dst + (*nt).width;
	    if (c_end_dst < kc) {
	      continue;
	    }
	    if (c_bgn_dst >= kc_end) {
	      break;
	    }
	      
	    const int ir_bgn = imax(kr, r_bgn_dst); 
	    const int ir_end = imin(kr_end, r_end_dst); 
	    const int jc_bgn = imax(kc, c_bgn_dst); 
	    const int jc_end = imin(kc_end, c_end_dst); 
	    const long ops = ((long)(ir_end - ir_bgn) * 
			      (long)(jc_end - jc_bgn) * 2L);
	    int ir_bgn_srcs[2], ir_end_srcs[2];
	    int jc_bgn_srcs[2], jc_end_srcs[2];
	    ir_bgn_srcs[0] = (*mt).begin_src0 + (ir_bgn - r_bgn_dst);
	    ir_bgn_srcs[1] = (*mt).begin_src1 + (ir_bgn - r_bgn_dst);
	    jc_bgn_srcs[0] = (*nt).begin_src0 + (jc_bgn - c_bgn_dst);
	    jc_bgn_srcs[1] = (*nt).begin_src1 + (jc_bgn - c_bgn_dst);
	      
	    ir_end_srcs[0] = ir_bgn_srcs[0] + (ir_end - ir_bgn);
	    ir_end_srcs[1] = ir_bgn_srcs[1] + (ir_end - ir_bgn);
	    jc_end_srcs[0] = jc_bgn_srcs[0] + (jc_end - jc_bgn);
	    jc_end_srcs[1] = jc_bgn_srcs[1] + (jc_end - jc_bgn);

	    if ((ir_bgn < ir_end) && (jc_bgn < jc_end)){ 
	      C_Dsub_task<T> *tmp = 
		new C_Dsub_task<T>(1,
				   0,
				   (ir_bgn - kr), 
				   (ir_end - kr), 
				   (jc_bgn - kc), 
				   (jc_end - kc), 
				   ir_bgn_srcs[0],
				   ir_bgn_srcs[1],
				   jc_bgn_srcs[0],
				   jc_bgn_srcs[1],
				   father_row, 
				   father_diag_pt, 
				   (RectBlockMatrix<T> *)NULL, // dst_pt
				   kr0,
				   kc0,
				   (RectBlockMatrix<T> *)NULL,  // dst_pt2
				   (-1),
				   (-1),
				   child_pt[0],
				   child_pt[1],
				   isSym ? dsub_sym2sym_two<T> : ((kr0 == kc0) ? dsub_unsym2diag_two<T> : dsub_unsym2unsym_two<T>),
				   skip_flag,
				   ops,
				   father_id,
				   level,
				   verbose,
				   fp);
	      C_task_arg->push_back(tmp);
	      *ops_sum += skip_flag ? 0L : ops;
		
	      for (int ll = 0; ll < 2; ll++) {
		update_parents_list(parents_r[ll], 
				    ir_bgn_srcs[ll], ir_end_srcs[ll], 
				    child_pt[ll]);
		update_parents_list(parents_c[ll], 
				    jc_bgn_srcs[ll], jc_end_srcs[ll], 
				    child_pt[ll]);
#if 0
		if (father_id == 0 && ll == 1) {
		  fprintf(stderr, "%s : %d [%d %d] %d\n",
			  __FILE__, __LINE__,
			  jc_bgn_srcs[ll], jc_end_srcs[ll], parents_c[ll].back());
		}
#endif
	      }
	    } //	if ((ir_bgn < ir_end) && (jc_bgn < jc_end)) {
	  } // loop : nt
	  // -- child01 * child_ll
	  for (int ll = 0; ll < 2; ll++) { 
	    for (list <index_strip>::const_iterator nt = strips_r[ll].begin();
		 nt != strips_r[ll].end(); ++nt) {
	      const int c_bgn_dst = (*nt).begin_dst;
	      const int c_end_dst = c_bgn_dst + (*nt).width;
	      if (c_bgn_dst < r_bgn_dst) {
		continue;
	      }
	      if (c_end_dst < kc) {
		continue;
	      }
	      if (c_bgn_dst >= kc_end) {
		break;
	      }
	      const int ir_bgn = imax(kr, r_bgn_dst);
	      const int ir_end = imin(kr_end, r_end_dst); 
	      const int jc_bgn = imax(kc, c_bgn_dst); 
	      const int jc_end = imin(kc_end, c_end_dst); 
	      const long ops = ((long)(ir_end - ir_bgn) * 
				(long)(jc_end - jc_bgn));
	      const int ir_bgn_src = 
		((ll == 0 ? (*mt).begin_src0 : (*mt).begin_src1) +
		 (ir_bgn - r_bgn_dst));
	      const int jc_bgn_src = 
		// ((*nt).begin_dst + (jc_bgn - c_bgn_dst));
		(*nt).begin_src + (jc_bgn - c_bgn_dst);
	      const int ir_end_src = ir_bgn_src + (ir_end - ir_bgn);
	      const int jc_end_src = jc_bgn_src + (jc_end - jc_bgn);

	      if ((ir_bgn < ir_end) && (jc_bgn < jc_end)) {
		C_Dsub_task<T> *tmp = 
		  new C_Dsub_task<T>(1,
				     0,
				     (ir_bgn - kr), 
				     (ir_end - kr), 
				     (jc_bgn - kc), 
				     (jc_end - kc), 
				     ir_bgn_src, // ((*mt).begin_dst[ll] + ...
				     (-1), // ir_bgn_src2
				     jc_bgn_src,
				     (-1), // jc_bgn_src2
				     father_row, 
				     father_diag_pt, 
				     (RectBlockMatrix<T> *)NULL, // dst_pt
				     kr0,
				     kc0,
				     (RectBlockMatrix<T> *)NULL,  // dst_pt2
				     (-1),
				     (-1),
				     child_pt[ll],
				     (SquareBlockMatrix<T>*)NULL, // src_pt2
				     isSym ? dsub_sym2sym<T> : ((kr0 == kc0) ? dsub_unsym2diag<T> : dsub_unsym2unsym<T>),
				     false,
				     ops,
				     father_id,
				     level,
				     verbose,
				     fp);
		C_task_arg->push_back(tmp);
		*ops_sum += ops;
		  
		update_parents_list(parents_r[ll], 
				    ir_bgn_src, ir_end_src, child_pt[ll]);
		update_parents_list(parents_c[ll], 
				    jc_bgn_src, jc_end_src, child_pt[ll]);
#if 0
		if (father_id == 0 && ll == 1) {
		  fprintf(stderr, "%s : %d [%d %d] %d : %d %d %d\n",
			  __FILE__, __LINE__,
			  jc_bgn_src, jc_end_src, parents_c[ll].back(),
			  (*nt).begin_src, (*nt).begin_dst, (*nt).width);
		}
#endif
	      } 	//      if ((ir_bgn < ir_end) && (jc_bgn < jc_end)) 
	    } // loop : nt
	  }   // loop : ll
	}    // loop : mt
	string task_name = ("g " + to_string(kr0) + " " + to_string(kc0)
			    + " : " + to_string(level) + " : "
			    + to_string(father_id + 1)); // selfIndex^-1
	queue[kk] = new C_task(C_DSUB,
			       task_name, //task_name.str(),
			       //			 task_name_cstr,
			       (void *)C_task_arg,
			       C_Dsub_task_exec<T>,
			       1,
			       0,
			       ops_sum);
	queue[kk]->parallel_max = block_diag_size;
	queue[kk]->parallel_id = kk;
	// added to manage parents of tasks_p and tasks_q
	if (tasks_r != NULL) {
	  for (typename list<child_contribution<T> >::const_iterator jt =
		 child_contrib.begin();
	       jt != child_contrib.end(); ++jt) {
	    queue[kk]->parents->push_back(tasks_r[(*jt).child_id][0]); 
	  }
	}
	if (tasks_q != NULL) {
	  queue[kk]->parents->push_back(tasks_q[father_id][0]); // diag
	}
	if (tasks_s != NULL) {
	  const int itmp = (kc0 * (kc0 + 1)) / 2 + kr0;
	  queue[kk]->parents->push_back((*tasks_s)[itmp]);     
	}
	//	if (tasks_p != NULL) {
#ifdef DEBUG_QUEUE_GENERATION
	cout << "++ diag two children " << task_name.str().c_str() << " ";
	list<C_Dsub_task<T> *> *task_tmp =
	  (list<C_Dsub_task<T> *> *)queue[kk]->func_arg;
	cout << "task size = " << task_tmp->size() << " : " << endl;
	for (list<C_Dsub_task<T> *>::const_iterator jt = task_tmp->begin();
	     jt != task_tmp->end(); jt++) {
	    cout << (*jt)->atomic_id << " / " << (*jt)->atomic_size << " : "
		 << (*jt)->ir_bgn << " / " << (*jt)->ir_end << " : "
		 << (*jt)->jc_bgn << " / " << (*jt)->jc_end
		 << endl;
	}
#endif
	if (!skip_flag) {
	  for (int ll = 0; ll < 2; ll++) {
	    for (list<int>::const_iterator mt = parents_r[ll].begin(); 
		 mt != parents_r[ll].end(); ++mt) {
	      for (list<int>::const_iterator nt = mt;
		   nt != parents_r[ll].end(); ++nt) { // running upper
		const int idx = (isSym ? 
				 (((*nt) * ((*nt) + 1)) / 2 + (*mt)) :
				 ((*nt) * (*nt) + 2 * (*mt)));
		vector<C_task *> &tasks_tmq = tasks_d[child_id[ll]];
		vector<int> &indcolq = tasks_d_indcol[child_id[ll]];
		if (tasks_p_flag) {
		  vector<C_task *> &tasks_tmp = tasks_p[child_id[ll]];
		  vector<int> &indcolp = tasks_p_indcol[child_id[ll]];
		  queue[kk]->parents->push_back(tasks_tmp[indcolp[idx]]);
		}
		if (tasks_tmq.size() > 0) {
		  tasks_tmq[indcolq[idx]]->parents->push_back(queue[kk]);
		}
		if (!isSym && ((*nt) > (*mt))) {
		  if (tasks_p_flag) {
		    vector<C_task *> &tasks_tmp = tasks_p[child_id[ll]];
		    vector<int> &indcolp = tasks_p_indcol[child_id[ll]];
		    queue[kk]->parents->push_back(tasks_tmp[indcolp[idx + 1]]);
		  }
		  if (tasks_tmq.size() > 0) {
		    tasks_tmq[indcolq[idx + 1]]->parents->push_back(queue[kk]);
		  }
		}
	      }
	      for (list<int>::const_iterator nt = parents_c[ll].begin(); 
		   nt != parents_c[ll].end(); ++nt) {
		if ((*nt) > (*mt)) {
		  const int idx = (isSym ? 
				   (((*nt) * ((*nt) + 1)) / 2 + (*mt)) :
				   ((*nt) * (*nt) + 2 * (*mt)));
		  vector<C_task *> &tasks_tmq = tasks_d[child_id[ll]];
		  vector<int> &indcolq = tasks_d_indcol[child_id[ll]];
		  if (tasks_p_flag) {
		    vector<C_task *> &tasks_tmp = tasks_p[child_id[ll]];
		    vector<int> &indcolp = tasks_p_indcol[child_id[ll]];
		    //		    if (indcolp.size() > idx) {
		    queue[kk]->parents->push_back(tasks_tmp[indcolp[idx]]);
#if 0
		    }
		    else {
			fprintf(stderr,
				"%s %d : %d < %d : prnts_c[%d].size=%d %d %d kr0=%d kc0=%d\n",
				__FILE__, __LINE__,
				(int)indcolp.size(), idx,
				ll, parents_c[ll].size(), *mt, *nt, kr0, kc0);
		    }
#endif
		  }
		  if (tasks_tmq.size() > 0) {
		    //		    if (indcolq.size() > idx) {
		      tasks_tmq[indcolq[idx]]->parents->push_back(queue[kk]);
#if 0
		    }
		    else {
			fprintf(stderr, "%s %d : % d : (%d %d) %d < %d\n",
				__FILE__, __LINE__, ll, (*nt), (*mt),
				(int)indcolq.size(), idx);
		    }
		  #endif
		  if (!isSym) {
		    if (tasks_p_flag) {
		      vector<C_task *> &tasks_tmp = tasks_p[child_id[ll]];
		      vector<int> &indcolp = tasks_p_indcol[child_id[ll]];
		      //		      if (indcolp.size() > (idx + 1)) {
			queue[kk]->parents->push_back(tasks_tmp[indcolp[idx + 1]]);		      }
#if 0
		      else {
			fprintf(stderr, "%s %d : %d < %d\n", __FILE__, __LINE__,
				(int)indcolp.size(), (idx + 1));
		      }
#endif
		    }
		    if (tasks_tmq.size() > 0) {
		      //		      if (indcolq.size() > (idx + 1)) {
			tasks_tmq[indcolq[idx + 1]]->parents->push_back(queue[kk]);
#if 0
		      }
		      else {
			fprintf(stderr, "%s %d : %d < %d\n", __FILE__, __LINE__,
				(int)indcolq.size(), (idx + 1));
		      }
#endif
		    }
		  }
		}
	      }  // loop : nt
	    } // loop :
	  }   // loop :  ll
	}   //  if (!skip_flag)
	else {
	  // verify skip
	  if (verbose ) {
	    for (int ll = 0; ll < 2; ll++) {
	      for (list<int>::const_iterator mt = parents_r[ll].begin(); 
		   mt != parents_r[ll].end(); ++mt) {
		if ((*mt) > child_pt[ll]->num_blocks0()) {
		  fprintf(stderr, "%s %d : %d incorrect skip : %d %d\n", 
			  __FILE__, __LINE__, ll, (*mt), 
			  child_pt[ll]->num_blocks0());
		}
	      }
	      for (list<int>::const_iterator mt = parents_c[ll].begin(); 
		   mt != parents_c[ll].end(); ++mt) {
		if ((*mt) > child_pt[ll]->num_blocks0()) {
		  fprintf(stderr, "%s %d : %d incorrect skip : %d %d\n", 
			  __FILE__, __LINE__, ll, (*mt), 
			  child_pt[ll]->num_blocks0());
		}
	      }
	    }
	  }  // if (verbose)
	}
	// }   // 	if (tasks_p != NULL)
	queue[kk]->parents->sort(compare_task_name);
	queue[kk]->parents->unique();
	EraseNullParents(queue[kk]);
	for (int m = 0; m < 2; m++) {
	  parents_r[m].clear();
	  parents_c[m].clear();
	}
	delete [] parents_r;
	delete [] parents_c;
	kk++;
      }     // loop : kr0
    }      // loop : kc0
    // offdiagonal contribution
   
    for (int kc0 = 0; kc0 < num_block2; kc0++) {
      //      const int kc = kc0 * SIZE_B1;
      //      const int kc_end = ((kc0 == num_block2 - 1) ? offdiag_size : 
      //			  (kc + SIZE_B1));
      const int kc = father_offdiag_pt->IndexBlock_c(kc0);
      const int kc_end = father_offdiag_pt->IndexBlock_c(kc0 + 1);

      list<C_Dsub_task<T> *> *C_task_arg = new list<C_Dsub_task<T> *>;
      list<int>* parents_r = new list<int>[2];
      list<int>* parents_c = new list<int>[2];
      long *ops_sum = new long;
      *ops_sum = 0L;
	
      // ---- child_ll * (child_ll + child01) ----
      for (int ll = 0; ll < 2; ll++) { 
	for (int kr0 = 0; kr0 < num_block; kr0++) { // 02 Jul.2014 : Atsushi
	  const int kr = kr0 * SIZE_B1;
	  const int kr_end = (kr0 == num_block - 1) ? diag_size : (kr + SIZE_B1);
	  for (list <index_strip>::const_iterator mt = strips_r[ll].begin();
	       mt != strips_r[ll].end(); ++mt) {

	    const int r_bgn_dst = (*mt).begin_dst;
	    const int r_end_dst = r_bgn_dst + (*mt).width;
	    if (r_end_dst < kr) {
	      continue;  // of loop mt 
	    }
	    if (r_bgn_dst >= kr_end ) {
	      break;     // of loop mt
	    }
	    const int ir_bgn = imax(kr, r_bgn_dst);
	    const int ir_end = imin(kr_end, r_end_dst);

	  // -- child_ll * child_ll	      
	    for (list <index_strip>::const_iterator nt = strips_c[ll].begin();
		 nt != strips_c[ll].end(); ++nt) {
	      const int c_bgn_dst = (*nt).begin_dst;
	      const int c_end_dst = c_bgn_dst + (*nt).width;
	      if (c_end_dst < kc) {
		continue;
	      }
	      if (c_bgn_dst >= kc_end) {
		break;
	      }
	      const int jc_bgn = imax(kc, c_bgn_dst); 
	      const int jc_end = imin(kc_end, c_end_dst); 
	      const long ops = ((long)(*mt).width * 
				(long)(jc_end - jc_bgn));
	      const int ir_bgn_src = (*mt).begin_src + (ir_bgn - r_bgn_dst);
	      const int jc_bgn_src = (*nt).begin_src + (jc_bgn - c_bgn_dst);
	      const int ir_end_src = ir_bgn_src + (ir_end - ir_bgn);
	      const int jc_end_src = jc_bgn_src + (jc_end - jc_bgn);
	      if ((ir_bgn < ir_end) && (jc_bgn < jc_end)) {
		C_Dsub_task<T> *tmp = 
		  new C_Dsub_task<T>(1,
				     0,
				     (ir_bgn - kr), 
				     (ir_end - kr), 
				     (jc_bgn - kc), 
				     (jc_end - kc), 
				     ir_bgn_src, // (*mt).begin_src,
				     (-1), // ir_bgn_src2
				     jc_bgn_src,
				     (-1), // jc_bgn_src2
				     father_row, 
				     (SquareBlockMatrix<T>*)NULL, // father_pt
				     father_offdiag_pt,
				     kr0, // ir_block
				     kc0, // jc_block
				     isSym ? (RectBlockMatrix<T> *)NULL : father_offdia2_pt,
				     0, // withou block
				     0, // withou block
				     child_pt[ll],
				     (SquareBlockMatrix<T>*)NULL, // src_pt2
				     isSym ? dsub_sym2rct<T> : dsub_unsym2rct<T>,
				     false,
				     ops,
				     father_id,
				     level,
				     verbose,
				     fp);
		C_task_arg->push_back(tmp);
		*ops_sum += ops;
		update_parents_list(parents_r[ll], 
				    ir_bgn_src, ir_end_src, child_pt[ll]);
		update_parents_list(parents_c[ll], 
				    jc_bgn_src, jc_end_src, child_pt[ll]);
	      } // if ((ir_bgn < ir_end) && (jc_bgn < jc_end))
	    } // loop : nt
	    // -- child_ll * child01	      
	    for (list <index_strip2>::const_iterator nt = strips_c01.begin();
		 nt != strips_c01.end(); ++nt) {
	      const int c_bgn_dst= (*nt).begin_dst;
	      const int c_end_dst = c_bgn_dst + (*nt).width;
	      if (c_end_dst < kc) {
		continue;
	      }
	      if (c_bgn_dst >= kc_end) {
		break;
	      }
	      const int jc_bgn = imax(kc, c_bgn_dst); 
	      const int jc_end = imin(kc_end, c_end_dst); 
	      const long ops = ((long)(*mt).width * 
				(long)(jc_end - jc_bgn));
	      const int ir_bgn_src = (*mt).begin_src + (ir_bgn - r_bgn_dst);
	      const int jc_bgn_src = 
		((ll == 0 ? (*nt).begin_src0 : (*nt).begin_src1) +
		 (jc_bgn - c_bgn_dst));
	      const int ir_end_src = ir_bgn_src + (ir_end - ir_bgn);
	      const int jc_end_src = jc_bgn_src + (jc_end - jc_bgn);
	      
	      if ((ir_bgn < ir_end) && (jc_bgn < jc_end)) {
		C_Dsub_task<T> *tmp = 
		  new C_Dsub_task<T>(1,
				     0,
				     (ir_bgn - kr), 
				     (ir_end - kr), 
				     (jc_bgn - kc), 
				     (jc_end - kc),
				     ir_bgn_src,
				     (-1), // ir_bgn_src2
				     jc_bgn_src,  // (*nt)->begin_src[ll] + ...
				     (-1), // jc_bgn_src2
				     father_row, 
				     (SquareBlockMatrix<T>*)NULL, // father_pt
				     father_offdiag_pt,
				     kr0, // ir_block
				     kc0, // jc_block
				     isSym ? (RectBlockMatrix<T> *)NULL : father_offdia2_pt,
				     0, // withou block
				     0, // withou block
				     child_pt[ll],
				     (SquareBlockMatrix<T>*)NULL, // src_pt2
				     isSym ? dsub_sym2rct<T> : dsub_unsym2rct<T>,
				     false,
				     ops,
				     father_id,
				     level,
				     verbose,
				     fp);
		C_task_arg->push_back(tmp);
		*ops_sum += ops;
		
		update_parents_list(parents_r[ll], 
				    ir_bgn_src, ir_end_src, child_pt[ll]);
		update_parents_list(parents_c[ll], 
				    jc_bgn_src, jc_end_src, child_pt[ll]);
	      } // if ((ir_bgn < ir_end) && (jc_bgn < jc_end)) 
	    } // loop : nt
	  }   // loop : mt
	} // loop kr0
      }    // loop : ll
      // ---- child01 * (child0 + child01 + child1) ----
      for (int kr0 = 0; kr0 < num_block; kr0++) { // 02 Jul.2014 : Atsushi
	const int kr = kr0 * SIZE_B1;
	const int kr_end = (kr0 == num_block - 1) ? diag_size : (kr + SIZE_B1);
	  
	for (list <index_strip2>::const_iterator mt = strips_r01.begin();
	     mt != strips_r01.end(); ++mt) {
	  
	  const int r_bgn_dst = (*mt).begin_dst;
	  const int r_end_dst = r_bgn_dst + (*mt).width;
	  if (r_end_dst < kr) {
	    continue;
	  }
	  if (r_bgn_dst >= kr_end) {
	    break;
	  }
	  const int ir_bgn = imax(kr, r_bgn_dst);
	  const int ir_end = imin(kr_end, r_end_dst);
	  
	  // -- child01 * child_ll
	  for (int ll = 0; ll < 2; ll++) {
	    for (list <index_strip>::const_iterator nt = strips_c[ll].begin();
		 nt != strips_c[ll].end(); ++nt) {
	      const int c_bgn_dst = (*nt).begin_dst;
	      const int c_end_dst = c_bgn_dst + (*nt).width;
	      if (c_end_dst < kc) {
		continue;
	      }
	      if (c_bgn_dst >= kc_end) {
		break;
	      }
	      const int jc_bgn = imax(kc, c_bgn_dst); 
	      const int jc_end = imin(kc_end, c_end_dst);
	      const long ops = ((long)(*mt).width * 
				(long)(jc_end - jc_bgn));
	      const int ir_bgn_src = 
		(ll == 0 ? (*mt).begin_src0 : (*mt).begin_src1) + 
		(ir_bgn - r_bgn_dst);
	      const int jc_bgn_src = ((*nt).begin_src + (jc_bgn - c_bgn_dst));
	      const int ir_end_src = ir_bgn_src + (ir_end - ir_bgn);
	      const int jc_end_src = jc_bgn_src + (jc_end - jc_bgn);
	      if ((ir_bgn < ir_end) && (jc_bgn < jc_end)) {
		C_Dsub_task<T> *tmp = 
		  new C_Dsub_task<T>(1,
				     0,
				     (ir_bgn - kr), 
				     (ir_end - kr), 
				     (jc_bgn - kc), 
				     (jc_end - kc), 
				     ir_bgn_src, // ((*mt).begin_src[ll] + ...
				     (-1), // ir_bgn_src2
				     jc_bgn_src,
				     (-1), // jc_bgn_src2
				     father_row, 
				     (SquareBlockMatrix<T>*)NULL, // father_pt
				     father_offdiag_pt, 
				     kr0, // ir_block
				     kc0, // jc_block
				     isSym ? (RectBlockMatrix<T> *)NULL : father_offdia2_pt,
				     0, // withou block
				     0, // withou block
				     child_pt[ll],
				     (SquareBlockMatrix<T>*)NULL, // src_pt2
				     isSym ? dsub_sym2rct<T>: dsub_unsym2rct<T>,
				     false,
				     ops,
				     father_id,
				     level,
				     verbose,
				     fp);
		C_task_arg->push_back(tmp);
		*ops_sum += ops;
		
		update_parents_list(parents_r[ll], 
				    ir_bgn_src, ir_end_src, child_pt[ll]);
		update_parents_list(parents_c[ll], 
				    jc_bgn_src, jc_end_src, child_pt[ll]);
	      } // if ((ir_bgn < ir_end) && (jc_bgn < jc_end)) 
	    } // loop : nt
	  }   // loop : ll
	  // -- child01 * child01
	  for (list <index_strip2>::const_iterator nt = strips_c01.begin();
	       nt != strips_c01.end(); ++nt) {
	    const int c_bgn_dst = (*nt).begin_dst;
	    const int c_end_dst = c_bgn_dst + (*nt).width;
	    if (c_end_dst < kc) {
	      continue;
	    }
	    if (c_bgn_dst >= kc_end) {
	      break;
	    }
	    const int jc_bgn = imax(kc, c_bgn_dst);
	    const int jc_end = imin(kc_end, c_end_dst); 
	    const long ops = ((long)(*mt).width * 
			      (long)(jc_end - jc_bgn) * 2L);
	    
	    int ir_bgn_srcs[2], ir_end_srcs[2];
	    int jc_bgn_srcs[2], jc_end_srcs[2];
	    ir_bgn_srcs[0] = (*mt).begin_src0 + (ir_bgn - r_bgn_dst);
	    ir_bgn_srcs[1] = (*mt).begin_src1 + (ir_bgn - r_bgn_dst);
	    jc_bgn_srcs[0] = (*nt).begin_src0 + (jc_bgn - c_bgn_dst);
	    jc_bgn_srcs[1] = (*nt).begin_src1 + (jc_bgn - c_bgn_dst);
	    
	    ir_end_srcs[0] = ir_bgn_srcs[0] + (ir_end - ir_bgn);
	    ir_end_srcs[1] = ir_bgn_srcs[1] + (ir_end - ir_bgn);
	    jc_end_srcs[0] = jc_bgn_srcs[0] + (jc_end - jc_bgn);
	    jc_end_srcs[1] = jc_bgn_srcs[1] + (jc_end - jc_bgn);
	    
	    if ((ir_bgn < ir_end) && (jc_bgn < jc_end)) {
	      C_Dsub_task<T> *tmp = 
		  new C_Dsub_task<T>(1,
				   0,
				   (ir_bgn - kr), 
				   (ir_end - kr), 
				   (jc_bgn - kc), 
				   (jc_end - kc), 
				   ir_bgn_srcs[0],
				   ir_bgn_srcs[1],
				   jc_bgn_srcs[0],
				   jc_bgn_srcs[1],
				   father_row, 
				   (SquareBlockMatrix<T>*)NULL, // father_pt
				   father_offdiag_pt,
				   kr0, // ir_block
				   kc0, // jc_block
				   isSym ? (RectBlockMatrix<T> *)NULL : father_offdia2_pt,
				   0, // withou block
				   0, // withou block
				   child_pt[0],
				   child_pt[1],
				   isSym ? dsub_sym2rct_two<T> : dsub_unsym2rct_two<T>,
				   false,
				   ops,
				   father_id,
				   level,
				   verbose,
				   fp);
	      C_task_arg->push_back(tmp);
	      *ops_sum += ops;
	      for (int ll = 0; ll < 2; ll++) {
		update_parents_list(parents_r[ll], 
				    ir_bgn_srcs[ll], ir_end_srcs[ll], 
				    child_pt[ll]);
		update_parents_list(parents_c[ll], 
				    jc_bgn_srcs[ll], jc_end_srcs[ll], 
				    child_pt[ll]);
	      } // loop : ll
	    } //   if ((ir_bgn < ir_end) && (jc_bgn < jc_end))
	  } // loop : nt
	}    // loop : mt
      } // loop : kr0
      string task_name = ("h    " + to_string(kc0) + " : "
			  + to_string(level) + " : "
			  + to_string(father_id + 1)); // selfIndex^-1
      queue[kk] = new C_task(C_DSUB,
			     task_name, //task_name.str(),
			     (void *)C_task_arg,
			     C_Dsub_task_exec<T>,
			     1,
			     0,
			     ops_sum);
      queue[kk]->parallel_max = num_block2;
      queue[kk]->parallel_id = kc0;
      // added to manage parents of tasks_p and tasks_q
      if (tasks_r != NULL) {
	for (typename list<child_contribution<T> >::const_iterator jt =
	       child_contrib.begin();
	     jt != child_contrib.end(); ++jt) {
	  queue[kk]->parents->push_back(tasks_r[(*jt).child_id][0]); 
	} 
      }
      if (tasks_q != NULL) {
	queue[kk]->parents->push_back(tasks_q[father_id][1]); //  offdiag
      }
      if (tasks_s != NULL) {   // offiag
	for (int i = 0; i < num_block; i++) {
	  const int itmp = kc0 * num_block + i + (*tasks_s)[0]->parallel_max; 
	  queue[kk]->parents->push_back((*tasks_s)[itmp]);     
	}
      }
      //      if (tasks_p != NULL) {
      for (int ll = 0; ll < 2; ll++) {
	for (list<int>::const_iterator mt = parents_r[ll].begin(); 
	     mt != parents_r[ll].end(); ++mt) {
	  for (list<int>::const_iterator nt = parents_c[ll].begin(); 
	       nt != parents_c[ll].end(); ++nt) {
	    if ((*nt) >= (*mt)) {  // upper block
	      const int idx = (isSym ? 
			       (((*nt) * ((*nt) + 1)) / 2 + (*mt)) :
			       ((*nt) * (*nt) + 2 * (*mt)));

	      vector<C_task *> &tasks_tmq = tasks_d[child_id[ll]];
	      vector<int> &indcolq = tasks_d_indcol[child_id[ll]];
	      if (tasks_p_flag) {
		vector<C_task *> &tasks_tmp = tasks_p[child_id[ll]];
		vector<int> &indcolp = tasks_p_indcol[child_id[ll]];
		queue[kk]->parents->push_back(tasks_tmp[indcolp[idx]]);
	      }
	      if (tasks_tmq.size() > 0) {
		tasks_tmq[indcolq[idx]]->parents->push_back(queue[kk]);
	      }
	      if (!isSym && ((*nt) > (*mt))) {
		if (tasks_p_flag) {
		  vector<C_task *> &tasks_tmp = tasks_p[child_id[ll]];
		  vector<int> &indcolp = tasks_p_indcol[child_id[ll]];
		  queue[kk]->parents->push_back(tasks_tmp[indcolp[idx + 1]]);
		}
		if (tasks_tmq.size() > 0) {
		  tasks_tmq[indcolq[idx + 1]]->parents->push_back(queue[kk]);
		}
	      }
	    }
	  }  // loop : nt
	}    // loop : mt
      } // loop : ll
	//      }//     if (tasks_p != NULL) {
      queue[kk]->parents->sort(compare_task_name);
      queue[kk]->parents->unique();
      EraseNullParents(queue[kk]);
      for (int m = 0; m < 2; m++) {
	parents_r[m].clear();
	parents_c[m].clear();
      }
      delete [] parents_r;
      delete [] parents_c;
      kk++;
    } // loop : kc0
    delete [] strips_r;
    delete [] strips_c;
    delete [] child_pt;
    delete [] child_id;
#ifdef DEBUG_QUEUE_GENERATION2
    cout << "father = " << (father_id + 1) << " contrib = 2 : queue size = " 
	 << queue.size() << endl;
    for (vector<C_task *>::const_iterator kt = queue.begin();
	 kt != queue.end(); ++kt) {
      cout << (*kt)->task_name << " :: " << (*kt)->parents->size() << " : " ;
      for (list<C_task *>::const_iterator jt = (*kt)->parents->begin();
	   jt != (*kt)->parents->end(); ++jt) {
	cout << (*jt)->task_name << " / ";
      }
      cout << endl;
    } // loop : kt
#endif
#ifdef DEBUG_PREPARE_THREAD_DEBUG
    cout << "kk = " << kk << endl;
#endif
  } // if (child_contrib.size() == 2) {
  else {
    queue.resize(block_diag_size + block_offdiag_size);
#ifdef  DEBUG_PREPARE_THREAD_DEBUG      
    cout << "resize = " << block_diag_size + block_offdiag_size << " ";
#endif
    int kk = 0;
    // -- diagonal 
    for (int kc0 = 0; kc0 < num_block; kc0++) {
      const int kc = kc0 * SIZE_B1;
      const int kc_end = (kc0 == (num_block - 1)) ? diag_size : (kc + SIZE_B1);
      for (int kr0 = 0; kr0 <= kc0; kr0++) {
	const int kr = kr0 * SIZE_B1;
	const int kr_end = 
	  (kr0 == (num_block - 1)) ? diag_size : (kr + SIZE_B1);
	list<C_Dsub_task<T> *> *C_task_arg = new list<C_Dsub_task<T> *>;
	// C_task_arg[kk] = new list<C_Dsub_task *>;
	list<int>* parents_r = new list<int>[child_contrib.size()];
	list<int>* parents_c = new list<int>[child_contrib.size()];
	long *ops_sum = new long;
	*ops_sum = 0L;
	  
	int count_atomics = 0;
	for (typename list<child_contribution<T> >::const_iterator jt = child_contrib.begin();
	     jt != child_contrib.end(); ++jt) {
	  for (list <index_strip>::const_iterator mt = (*jt).diag_strip.begin();
	       mt != (*jt).diag_strip.end(); ++mt) {
	    const int r_bgn_dst = (*mt).begin_dst;
	    const int r_end_dst = r_bgn_dst + (*mt).width;
	    if (r_end_dst < kr || r_end_dst < kc) {
	      continue;
	    }
	    if (r_bgn_dst >= kr_end || r_bgn_dst >= kc_end) {
	      break;
	    }	    
	    count_atomics++;
	  }
	  for (list <index_strip>::const_iterator mt = (*jt).diag_strip.begin();
	       mt != (*jt).diag_strip.end(); ++mt) {
	    const int r_bgn_dst = (*mt).begin_dst;
	    const int r_end_dst = r_bgn_dst + (*mt).width;
	    if (r_end_dst < kr) {
	      continue;
	    }
	    if (r_bgn_dst >= kr_end ) {
	      break;
	    }
	    list <index_strip>::const_iterator nt = mt;
	    ++nt;
	    for (; nt != (*jt).diag_strip.end(); ++nt) {
	      const int c_bgn_dst = (*nt).begin_dst;
	      const int c_end_dst = c_bgn_dst + (*nt).width;
	      if (c_end_dst < kc) {
		continue;
	      }
	      if (c_bgn_dst >= kc_end) {
		break;
	      }
	      count_atomics++;
	    }  // loop : nt
	  } // loop  : mt
	} // loop : jt
	  
	int atomic_id = 0;
	int ll = 0;
	for (typename list<child_contribution<T> >::const_iterator jt = child_contrib.begin();
	     jt != child_contrib.end(); ++jt, ll++) {
	  SquareBlockMatrix<T>* child_pt = (*jt).child_pt;
	  SquareBlockMatrix<T>* father_diag_pt = (*jt).father_diag_pt;
	  for (list <index_strip>::const_iterator mt = (*jt).diag_strip.begin();
	       mt != (*jt).diag_strip.end(); ++mt) {
	    const int r_bgn_dst = (*mt).begin_dst;
	    const int r_end_dst = r_bgn_dst + (*mt).width;
	    if (r_end_dst < kr) {
	      continue;
	    }
	    if (r_bgn_dst >= kr_end) {
	      break;
	    }
	    const int ir_bgn = imax(kr, r_bgn_dst); 
	    const int ir_end = imin(kr_end, r_end_dst); 
	    const int ic_bgn = imax(kc, r_bgn_dst); 
	    const int ic_end = imin(kc_end, r_end_dst); 
	    if (kr == kc) {		
	      const long ops = ((long)(ir_end - ir_bgn) *
				((long)(ir_end - ir_bgn) + 1L) / 2L);
	      const int ir_bgn_src = 
		((*mt).begin_src + (ir_bgn - r_bgn_dst));
	      const int ir_end_src = ir_bgn_src + (ir_end - ir_bgn);
	      if (ir_bgn < ir_end) { 
		C_Dsub_task<T> *tmp = 
		  new C_Dsub_task<T>(count_atomics,
				     atomic_id++,
				     (ir_bgn - kr), 
				     (ir_end - kr), 
				     (-1), // jc_bgn
				     (-1), // jc_end
				     ir_bgn_src,
				     (-1), // ir_bgn_src2
				     (-1), // jc_bgn_src
				     (-1), // jc_bgn_src2
				     father_row, 
				     father_diag_pt, 
				     (RectBlockMatrix<T> *)NULL, // dst_pt
				     kr0,
				     kr0,
				     (RectBlockMatrix<T> *)NULL,  // dst_pt2
				     (-1),
				     (-1),
				     child_pt,  
				     (SquareBlockMatrix<T>*)NULL, // src_pt2
				     (isSym ? dsub_sym2sym_diag<T> : 
				      dsub_unsym2unsym_diag<T>), 
				     false,
				     ops,
				     father_id,
				     level,
				     verbose,
				     fp);
		C_task_arg->push_back(tmp);
		*ops_sum += ops;
		update_parents_list(parents_r[ll], 
				    ir_bgn_src, ir_end_src, child_pt);
	      } //  if (ir_bgn < ir_end)
	    }
	    else { // kc > kr
	      const long ops = ((long)(ir_end - ir_bgn) *
				(long)(ic_end - ic_bgn));
	      const int ir_bgn_src = 
		((*mt).begin_src + (ir_bgn - r_bgn_dst));
	      const int ic_bgn_src = 
		((*mt).begin_src + (ic_bgn - r_bgn_dst));
	      const int ir_end_src = ir_bgn_src + (ir_end - ir_bgn);
	      const int ic_end_src = ic_bgn_src + (ic_end - ic_bgn);
	      if ((ir_bgn < ir_end) && (ic_bgn < ic_end)) { 
		C_Dsub_task<T> *tmp = 
		  new C_Dsub_task<T>(count_atomics,
				     atomic_id++,
				     (ir_bgn - kr), 
				     (ir_end - kr), 
				     (ic_bgn - kc), 
				     (ic_end - kc), 
				     ir_bgn_src,
				     (-1),  // ir_bgn_src2
				     ic_bgn_src,
				     (-1),  // jc_bgn_src2
				     father_row, 
				     father_diag_pt,
				     (RectBlockMatrix<T> *)NULL, // dst_pt
				     kr0,
				     kc0,
				     (RectBlockMatrix<T> *)NULL,  // dst_pt2
				     (-1),
				     (-1),
				     child_pt,
				     (SquareBlockMatrix<T>*)NULL, // src_pt2
				     isSym ? dsub_sym2sym<T> : dsub_unsym2unsym<T>,
				     false,
				     ops,
				     father_id,
				     level,
				     verbose,
				     fp);
		C_task_arg->push_back(tmp);
		*ops_sum += ops;
		update_parents_list(parents_r[ll], 
				    ir_bgn_src, ir_end_src, child_pt);
		update_parents_list(parents_c[ll], 
				    ic_bgn_src, ic_end_src, child_pt);
	      } // if ((ir_bgn < ir_end) && (ic_bgn < ic_end)) 
	    }   // if (kr == kc)
	  }     // loop mt
	  //
	  for (list <index_strip>::const_iterator mt = (*jt).diag_strip.begin();
	       mt != (*jt).diag_strip.end(); ++mt) {
	    const int r_bgn_dst = (*mt).begin_dst;
	    const int r_end_dst = r_bgn_dst + (*mt).width;
	    if (r_end_dst < kr) {
	      continue;
	    }
	    if (r_bgn_dst >= kr_end ) {
	      break;
	    }
	    list <index_strip>::const_iterator nt = mt;
	    ++nt;
	    for (; nt != (*jt).diag_strip.end(); ++nt) {
	      const int c_bgn_dst = (*nt).begin_dst;
	      const int c_end_dst = c_bgn_dst + (*nt).width;
	      if (c_end_dst < kc) {
		continue;
	      }
	      if (c_bgn_dst >= kc_end) {
		break;
	      }
	      const int ir_bgn = imax(kr, r_bgn_dst);
	      const int ir_end = imin(kr_end, r_end_dst); 
	      const int jc_bgn = imax(kc, c_bgn_dst); 
	      const int jc_end = imin(kc_end, c_end_dst);
	      // index runs upper diagonal block ( mt < nt ) => (i < j)
	      const long ops = ((long)(ir_end - ir_bgn) * 
				(long)(jc_end - jc_bgn));
	      const int ir_bgn_src = ((*mt).begin_src + (ir_bgn - r_bgn_dst));
	      const int jc_bgn_src = ((*nt).begin_src + (jc_bgn - c_bgn_dst));
	      const int ir_end_src = ir_bgn_src + (ir_end - ir_bgn);
	      const int jc_end_src = jc_bgn_src + (jc_end - jc_bgn);
	      if ((ir_bgn < ir_end) && (jc_bgn < jc_end)) { 
		C_Dsub_task<T> *tmp = 
		  new C_Dsub_task<T>(count_atomics,
				     atomic_id++,
				     (ir_bgn - kr), 
				     (ir_end - kr), 
				     (jc_bgn - kc), 
				     (jc_end - kc), 
				     ir_bgn_src,
				     (-1), // ir_bgn_src2
				     jc_bgn_src,
				     (-1), // jc_bgn_src2
				     father_row, 
				     father_diag_pt, 
				     (RectBlockMatrix<T> *)NULL, // dst_pt
				     kr0,
				     kc0,
				     (RectBlockMatrix<T> *)NULL,  // dst_pt2
				     (-1),
				     (-1),
				     child_pt,
				     (SquareBlockMatrix<T>*)NULL, // src_pt2
				     isSym ? dsub_sym2sym<T> : ((kr0 == kc0) ? dsub_unsym2diag<T> : dsub_unsym2unsym<T>),
				     false,
				     ops,
				     father_id,
				     level,
				     verbose,
				     fp);
		//
		C_task_arg->push_back(tmp);
		*ops_sum += ops;
		update_parents_list(parents_r[ll], 
				    ir_bgn_src, ir_end_src, child_pt);
		update_parents_list(parents_c[ll], 
				    jc_bgn_src, jc_end_src, child_pt);
	      } // if ((ir_bgn < ir_end) && (jc_bgn < jc_end)) 
	    }     // loop : nt
	  }       // loop : mt
	}         // loop : jt
	string task_name = ("i " + to_string(kr0) + " " + to_string(kc0)
			    + " : " + to_string(level) + " : "
			    + to_string(father_id + 1)); // selfIndex^-1
	queue[kk] = new C_task(C_DSUB,
			       task_name,
			       (void *)C_task_arg,
			       C_Dsub_task_exec<T>,
			       1,
			       0,
			       ops_sum);
	queue[kk]->parallel_max = block_diag_size; // block_total_size;
	queue[kk]->parallel_id = kk;
	//  added to manage parents of tasks_p and tasks_q
	if (tasks_r != NULL) {
	  for (typename list<child_contribution<T> >::const_iterator jt = child_contrib.begin();
	       jt != child_contrib.end(); ++jt) {
	    queue[kk]->parents->push_back(tasks_r[(*jt).child_id][0]); 
	  }
	}
	if (tasks_q != NULL) {
	  queue[kk]->parents->push_back(tasks_q[father_id][0]); // diag
	}
	if (tasks_s != NULL) {
	  const int itmp = ((kc0 * (kc0 + 1)) / 2) + kr0;  // diag
	  queue[kk]->parents->push_back((*tasks_s)[itmp]);     
	}
	//	if (tasks_p != NULL) {
	ll = 0;
	for (typename list<child_contribution<T> >::const_iterator jt = child_contrib.begin();
	     jt != child_contrib.end(); ++jt, ll++) {
	  for (list<int>::const_iterator mt = parents_r[ll].begin(); 
	       mt != parents_r[ll].end(); ++mt) {
	    for (list<int>::const_iterator nt = mt;
		 nt != parents_r[ll].end(); ++nt) {
	      const int idx = (isSym ? 
			       (((*nt) * ((*nt) + 1)) / 2 + (*mt)) :
			       ((*nt) * (*nt) + 2 * (*mt)));
	      
	      vector<C_task *> &tasks_tmq = tasks_d[(*jt).child_id];
	      vector<int> &indcolq = tasks_d_indcol[(*jt).child_id];
	      if (tasks_p_flag) {
		vector<C_task *> &tasks_tmp = tasks_p[(*jt).child_id];
		vector<int> &indcolp = tasks_p_indcol[(*jt).child_id];
		queue[kk]->parents->push_back(tasks_tmp[indcolp[idx]]);
	      }
	      if (tasks_tmq.size() > 0) {
		tasks_tmq[indcolq[idx]]->parents->push_back(queue[kk]);
	      }
	      if (!isSym && ((*nt) > (*mt))) { 
		if (tasks_p_flag) {
		  vector<C_task *> &tasks_tmp = tasks_p[(*jt).child_id];
		  vector<int> &indcolp = tasks_p_indcol[(*jt).child_id];
		  queue[kk]->parents->push_back(tasks_tmp[indcolp[idx + 1]]);
		}
		if (tasks_tmq.size() > 0) {
		  tasks_tmq[indcolq[idx + 1]]->parents->push_back(queue[kk]);
		}
	      }
	    }
	    for (list<int>::const_iterator nt = parents_c[ll].begin(); 
		 nt != parents_c[ll].end(); ++nt) {
	      if ((*nt) > (*mt)) {  // upper block
		const int idx = (isSym ? 
				 (((*nt) * ((*nt) + 1)) / 2 + (*mt)) :
				 ((*nt) * (*nt) + 2 * (*mt)));
		vector<C_task *> &tasks_tmq = tasks_d[(*jt).child_id];
		vector<int> &indcolq = tasks_d_indcol[(*jt).child_id];
		if (tasks_p_flag) {
		  vector<C_task *> &tasks_tmp = tasks_p[(*jt).child_id];
		  vector<int> &indcolp = tasks_p_indcol[(*jt).child_id];
		  queue[kk]->parents->push_back(tasks_tmp[indcolp[idx]]);
		}
		if (tasks_tmq.size() > 0) {
		  tasks_tmq[indcolq[idx]]->parents->push_back(queue[kk]);
		}
		if (!isSym) {
		  if (tasks_p_flag) {
		    vector<C_task *> &tasks_tmp = tasks_p[(*jt).child_id];
		    vector<int> &indcolp = tasks_p_indcol[(*jt).child_id];
		    queue[kk]->parents->push_back(tasks_tmp[indcolp[idx + 1]]);
		  }
		  if (tasks_tmq.size() > 0) {
		    tasks_tmq[indcolq[idx + 1]]->parents->push_back(queue[kk]);
		  }
		}
	      }
	    } // loop : nt
	  }   // loop : mt
	}    // loop : jt
	queue[kk]->parents->sort(compare_task_name);
	queue[kk]->parents->unique();
	EraseNullParents(queue[kk]);
	for (int m = 0; m < child_contrib.size(); m++) {
	  parents_r[m].clear();
	  parents_c[m].clear();
	}
	delete [] parents_r;
	delete [] parents_c;
	kk++;
      }          // loop : kc
    }            // loop : kr
    // -- offdiagonal
    for (int kc0 = 0; kc0 < num_block2; kc0++) {
      const int kc = child_contrib.front().father_offdiag_pt->IndexBlock_c(kc0);
      const int kc_end = child_contrib.front().father_offdiag_pt->IndexBlock_c(kc0 + 1);

      for (int kr0 = 0; kr0 < num_block; kr0++) {
	const int kr = kr0 * SIZE_B1;
	const int kr_end = 
	  (kr0 == (num_block - 1)) ? diag_size : (kr + SIZE_B1);

	list<C_Dsub_task<T> *> *C_task_arg = new list<C_Dsub_task<T> *>;
	// C_task_arg[kk] = new list<C_Dsub_task *>;
	list<int>* parents_r = new list<int>[child_contrib.size()];
	list<int>* parents_c = new list<int>[child_contrib.size()];
	long *ops_sum = new long;
	*ops_sum = 0L;
	int count_atomics = 0;
	  
	for (typename list<child_contribution<T> >::const_iterator jt = child_contrib.begin();
	     jt != child_contrib.end(); ++jt) {
	    
	  for (list <index_strip>::const_iterator mt = (*jt).diag_strip.begin();
	       mt != (*jt).diag_strip.end(); ++mt) {
	    const int r_bgn_dst = (*mt).begin_dst;
	    const int r_end_dst = r_bgn_dst + (*mt).width;
	    if (r_end_dst < kr) {
	      continue;  // of loop mt 
	    }
	    if (r_bgn_dst >= kr_end ) {
	      break;     // of loop mt
	    }
	    for (list <index_strip>::const_iterator nt = 
		   (*jt).offdiag_strip.begin();
		 nt != (*jt).offdiag_strip.end(); ++nt) {
	      const int c_bgn_dst = (*nt).begin_dst;
	      const int c_end_dst = c_bgn_dst + (*nt).width;
	      if (c_end_dst < kc) {
		continue;
	      }
	      if (c_bgn_dst >= kc_end) {
		break;
	      }
	      count_atomics++;
	    }
	  }
	} // loop : jt
	  
	int atomic_id = 0;
	int ll = 0;
	for (typename list<child_contribution<T> >::const_iterator jt = child_contrib.begin();
	     jt != child_contrib.end(); ++jt, ll++) {
	    
	  const int father_row = (*jt).father_row;
	  SquareBlockMatrix<T> *child_pt = (*jt).child_pt;
	  RectBlockMatrix<T> *father_offdiag_pt = (*jt).father_offdiag_pt;
	  RectBlockMatrix<T> *father_offdia2_pt = (*jt).father_offdiag_unsym_pt;
	    
	  for (list <index_strip>::const_iterator mt = (*jt).diag_strip.begin();
	       mt != (*jt).diag_strip.end(); ++mt) {
	    const int r_bgn_dst = (*mt).begin_dst;
	    const int r_end_dst = r_bgn_dst + (*mt).width;
	    if (r_end_dst < kr) {
	      continue;  // of loop mt 
	    }
	    if (r_bgn_dst >= kr_end ) {
	      break;     // of loop mt
	    }
	    const int ir_bgn = imax(kr, r_bgn_dst);
	    const int ir_end = imin(kr_end, r_end_dst);
	    for (list <index_strip>::const_iterator nt = (*jt).offdiag_strip.begin();
		 nt != (*jt).offdiag_strip.end(); ++nt) {
	      const int c_bgn_dst = (*nt).begin_dst;
	      const int c_end_dst = c_bgn_dst + (*nt).width;
	      if (c_end_dst < kc) {
		continue;
	      }
	      if (c_bgn_dst >= kc_end) {
		break;
	      }
	      const int jc_bgn = imax(kc, c_bgn_dst); 
	      const int jc_end = imin(kc_end, c_end_dst); 
	      //
	      // dsub_sym2rct(ir_bgn, ir_end, jc_bgn, jc_end, 
	      //               ir_bgn_src, jc_bgn_src, 
	      //               src_pt, src_pt) 
	      const long ops = ((long)(ir_end - ir_bgn) *
				(long)(jc_end - jc_bgn));
	      const int ir_bgn_src =(*mt).begin_src + (ir_bgn - r_bgn_dst);
	      const int jc_bgn_src =(*nt).begin_src + (jc_bgn - c_bgn_dst);
	      const int ir_end_src = ir_bgn_src + (ir_end - ir_bgn);
	      const int jc_end_src = jc_bgn_src + (jc_end - jc_bgn);

	      if ((ir_bgn < ir_end) && (jc_bgn < jc_end)) { 
		C_Dsub_task<T> *tmp = 
		  new C_Dsub_task<T>(count_atomics,
				     atomic_id++,
				     (ir_bgn - kr), 
				     (ir_end - kr), 
				     (jc_bgn - kc), 
				     (jc_end - kc), 
				     ir_bgn_src,
				     (-1), // ir_bgn_src2
				     jc_bgn_src,
				     (-1), // jc_bgn_src2
				     father_row, 
				     (SquareBlockMatrix<T>*)NULL, // father_pt
				     father_offdiag_pt, 
				     kr0, // ir_block
				     kc0, // jc_block
				     isSym ? (RectBlockMatrix<T> *)NULL : father_offdia2_pt,
				     0, // without block
				     0, // without block
				     child_pt,
				     (SquareBlockMatrix<T>*)NULL, // src_pt2
				     isSym ? dsub_sym2rct<T> : dsub_unsym2rct<T>,
				     false,
				     ops,
				     father_id,
				     level,
				     verbose,
				     fp);
		C_task_arg->push_back(tmp);
		*ops_sum += ops;
		update_parents_list(parents_r[ll], 
				    ir_bgn_src, ir_end_src, child_pt);
		update_parents_list(parents_c[ll], 
				    jc_bgn_src, jc_end_src, child_pt);
	      } // if ((ir_bgn < ir_end) && (jc_bgn < jc_end)) 
	    }     // loop : nt
	  }       // loop : mt
	}         // loop : jt
	string task_name = ("j " + to_string(kr0) + " " + to_string(kc0)
			    + " : " + to_string(level) + " : " +
			    to_string(father_id + 1)); // selfIndex^-1
	queue[kk] = new C_task(C_DSUB,
			       task_name, //  task_name.str(),
			       (void *)C_task_arg,
			       C_Dsub_task_exec<T>,
			       1,
			       0,
			       ops_sum);
	queue[kk]->parallel_max = block_offdiag_size; // block_total_size;
	queue[kk]->parallel_id = kk - block_diag_size;
	// added to manage parents of tasks_p and tasks_q
	if (tasks_r != NULL) {
	  for (typename list<child_contribution<T> >::const_iterator jt = child_contrib.begin();
	       jt != child_contrib.end();
	       ++jt) {
	    queue[kk]->parents->push_back(tasks_r[(*jt).child_id][0]); 
	  }
	}
	if (tasks_q != NULL) {
	  queue[kk]->parents->push_back(tasks_q[father_id][1]); //  offdiag
	}
	if (tasks_s != NULL) {   // offiag
	  const int itmp = kc0 * num_block + kr0 + (*tasks_s)[0]->parallel_max;
	  queue[kk]->parents->push_back((*tasks_s)[itmp]);     
	}
	//	if (tasks_p != NULL) {
	ll = 0;
	for (typename list<child_contribution<T> >::const_iterator jt = child_contrib.begin();
	     jt != child_contrib.end(); ++jt, ll++) {
	  for (list<int>::const_iterator mt = parents_r[ll].begin(); 
	       mt != parents_r[ll].end(); ++mt) {
	    for (list<int>::const_iterator nt = parents_c[ll].begin(); 
		 nt != parents_c[ll].end(); ++nt) {
	      // ?? offdiag : to be debugged : 27 Mar.2012 Atsushi
	      if ((*nt) >= (*mt)) {  // upper block
		const int idx = (isSym? 
				 (((*nt) * ((*nt) + 1)) / 2 + (*mt)) :
				 ((*nt) * (*nt) + 2 * (*mt)));
		vector<C_task *>&tasks_tmq = tasks_d[(*jt).child_id];
		vector<int> &indcolq = tasks_d_indcol[(*jt).child_id];
		if (tasks_p_flag) {
		  vector<C_task *>&tasks_tmp = tasks_p[(*jt).child_id];
		  vector<int> &indcolp = tasks_p_indcol[(*jt).child_id];
		  queue[kk]->parents->push_back(tasks_tmp[indcolp[idx]]);
		}
		if (tasks_tmq.size() > 0) {
		  tasks_tmq[indcolq[idx]]->parents->push_back(queue[kk]);
		}
		if (!isSym && ((*nt) > (*mt))) {
		  if (tasks_p_flag) {
		    vector<C_task *>&tasks_tmp = tasks_p[(*jt).child_id];
		    vector<int> &indcolp = tasks_p_indcol[(*jt).child_id];
		    queue[kk]->parents->push_back(tasks_tmp[indcolp[idx + 1]]);
		  }
		  if (tasks_tmq.size() > 0) {
		    tasks_tmq[indcolq[idx + 1]]->parents->push_back(queue[kk]);
		  }
		}
	      }  
	    } // loop : nt
	  }   // loop : mt
	}   // loop : jt
	//	} // if (tasks_p != NULL) {
	queue[kk]->parents->sort(compare_task_name);
	queue[kk]->parents->unique();
	EraseNullParents(queue[kk]);
	for (int m = 0; m < child_contrib.size(); m++) {
	  parents_r[m].clear();
	  parents_c[m].clear();
	}
	delete [] parents_r;
	delete [] parents_c;
	kk++;
      }           // loop : kr0
    }             // loop : kc0
  } // else (child_contrib.size() == 2)
#if 0
  if (verbose) {
    for (vector<C_task *>::const_iterator kt = queue.begin();
	 kt != queue.end(); ++kt) {
      fprintf(fp, "%s :: %d : parents = %d ",
	      (*kt)->task_name, (int)*((*kt)->ops_complexity),
	      (int)(*kt)->parents->size());
      for (list<C_task *>::const_iterator jt = (*kt)->parents->begin();
	   jt != (*kt)->parents->end(); ++jt) {
	fprintf(fp, "%s /", (*jt)->task_name);
      }
      fprintf(fp, "\n");
    }
  }
#endif
}

template
void C_Dsub_queue<double>(bool isSym, 
			  int father_id,
			  bool skip_flag,
			  vector<C_task *>& queue,
			  list <child_contribution<double> > &child_contrib,
			  vector<C_task *>* tasks_p, 
			  vector<int>* tasks_p_indcol,
			  const bool tasks_p_flag,
			  vector<C_task *>* tasks_q,
			  vector<C_task *>* tasks_r,
			  vector<C_task *>* tasks_s,
			  vector<C_task *>* tasks_d, 
			  vector<int>* tasks_d_indcol,
			  int level,
			  const bool verbose,
			  FILE *fp);

template
void C_Dsub_queue<quadruple>(bool isSym, 
			     int father_id,
			     bool skip_flag,
			     vector<C_task *>& queue,
			     list <child_contribution<quadruple> > &child_contrib,
			     vector<C_task *>* tasks_p, 
			     vector<int>* tasks_p_indcol,
			     const bool tasks_p_flag,
			     vector<C_task *>* tasks_q,
			     vector<C_task *>* tasks_r,
			     vector<C_task *>* tasks_s,
			     vector<C_task *>* tasks_d, 
			     vector<int>* tasks_d_indcol,
			     int level,
			     const bool verbose,
			     FILE *fp);
 
template
void C_Dsub_queue<complex<double> >(bool isSym, 
				    int father_id,
				    bool skip_flag,
				    vector<C_task *>& queue,
				    list <child_contribution<complex<double> > > &child_contrib,
				    vector<C_task *>* tasks_p, 
				    vector<int>* tasks_p_indcol,
				    const bool tasks_p_flag,
				    vector<C_task *>* tasks_q,
				    vector<C_task *>* tasks_r,
				    vector<C_task *>* tasks_s,
				    vector<C_task *>* tasks_d, 
				    vector<int>* tasks_d_indcol,
				    int level,
				    const bool verbose,
				    FILE *fp);

template
void C_Dsub_queue<complex<quadruple> >(bool isSym, 
				       int father_id,
				       bool skip_flag,
				       vector<C_task *>& queue,
				       list <child_contribution<complex<quadruple> > > &child_contrib,
				       vector<C_task *>* tasks_p, 
				       vector<int>* tasks_p_indcol,
				       const bool tasks_p_flag,
				       vector<C_task *>* tasks_q,
				       vector<C_task *>* tasks_r,
				       vector<C_task *>* tasks_s,
				       vector<C_task *>* tasks_d, 
				       vector<int>* tasks_d_indcol,
				       int level,
				       const bool verbose,
				       FILE *fp);
//
  
template<typename T>
void update_parents_list(list <int>& parents,
			 const int begin, const int end, 
			 SquareBlockMatrix<T>* mtrx)
			 //			 const int size_block)
{
  const int ibgn = mtrx->BlockIndex(begin);
  const int iend = mtrx->BlockIndex(end - 1);
  for (int m = ibgn; m <= iend; m++) {
    bool flag = true;
    for (list<int>::const_iterator it = parents.begin();
	 it != parents.end(); ++it) {
      if (*it == m) {
	flag = false;
	break; 
      }
    }
    if (flag) {
      parents.push_back(m);
    }
  }
}

template
void update_parents_list<double>(list <int>& parents,
				 const int begin, const int end, 
				 SquareBlockMatrix<double>* mtrx);

template
void update_parents_list<quadruple>(list <int>& parents,
				 const int begin, const int end, 
				 SquareBlockMatrix<quadruple>* mtrx);
template
void update_parents_list<complex<double> >(list <int>& parents,
					   const int begin, const int end, 
					   SquareBlockMatrix<complex<double> >* mtrx);

template
void update_parents_list<complex<quadruple> >(list <int>& parents,
					   const int begin, const int end, 
					   SquareBlockMatrix<complex<quadruple> >* mtrx);
//
