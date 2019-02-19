/*! \file   C_DFullLDLt.cpp
    \brief  block factorization routines LDL^t, LDU
    \author Atsushi Suzuki, Laboratoire Jacques-Louis Lions
    \date   Apr. 22th 2013            
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

#include "Driver/C_threads_tasks.hpp"
#include "Driver/C_KernDetect.hpp"
#include "Driver/C_BlasRoutines.hpp"
#include "Driver/DissectionDefault.hpp"
#include "Compiler/arithmetic.hpp"
#include "Algebra/VectorArray.hpp"
#include "Compiler/DissectionIO.hpp"
#include <time.h>

#ifdef BLAS_MKL
#define MKL_DOMATCOPY
#endif

template<typename T>
void C_dupdateb_Schur_diag(void *arg_)
{
  C_dupdateb_Schur_arg<T> *arg = (C_dupdateb_Schur_arg<T> *)arg_;
  const int task_position = arg->task_position;
  const int id_block = arg->id_block;
  const int id_level = arg->id_level;
  const int n = arg->n;
  const int nrow = arg->nrow; // size of diag. block whose matrix is factorized
  const int ncol = arg->ncol; // 
  const int i1_block = arg->i1_block;
  const int ii_block = arg->ii_block;
  SquareBlockMatrix<T> &D = *(arg->D);
  const int n0 = D.dim_kern_block(id_block);

  const T none(-1.0);
  const T one(1.0);

  if ((task_position % 2 == 1) && (id_level == 0)) {
    T *a = arg->a->addrCoefs();    //T *a = *(arg->a);
    D.copyBlockToArray(arg->ii_block, arg->ii_block, a, n);
  }

  if (arg->isSym) {
    // alpha=-1
    // beta=1
    C_gemm_symm<T>(ncol,
		   (nrow - n0), // nrow : zero padding is better ?
		   none, // alpha
		   D.addrCoefBlock(ii_block, i1_block), // lower transposed
		   D.nrowBlock(ii_block, i1_block),
		   D.addrCoefBlock(i1_block, ii_block), // upper
		   D.nrowBlock(i1_block, ii_block),
		   one, // beta,
		   D.addrCoefBlock(ii_block, ii_block), // Schur
		   D.nrowBlock(ii_block, ii_block));
  }
  else {
    blas_gemm<T>(CblasTrans, CblasNoTrans,
		 ncol,
		 ncol,
		 (nrow - n0), // nrow : zero padding is better ?
		 none, // alpha,
		 D.addrCoefBlock(ii_block, i1_block), // lower transposed
		 D.nrowBlock(ii_block, i1_block),
		 D.addrCoefBlock(i1_block, ii_block), // upper
		 D.nrowBlock(i1_block, ii_block),
		 one, // beta,
		 D.addrCoefBlock(ii_block, ii_block), // Schur
		 D.nrowBlock(ii_block, ii_block));
  }
}

template
void C_dupdateb_Schur_diag<double>(void *arg_);

template
void C_dupdateb_Schur_diag<quadruple>(void *arg_);

template
void C_dupdateb_Schur_diag<complex<double> >(void *arg_);

template
void C_dupdateb_Schur_diag<complex<quadruple> >(void *arg_);

template
void C_dupdateb_Schur_diag<float>(void *arg_);

template
void C_dupdateb_Schur_diag<complex<float> >(void *arg_);
//

template<typename T>
void C_dupdateb_Schur_offdiag(void *arg_)
{
  C_dupdateb_Schur_arg<T> *arg =  (C_dupdateb_Schur_arg<T> *)arg_;
  const int task_position = arg->task_position;
  const int id_block = arg->id_block;
  const int id_level = arg->id_level;
  const int n = arg->n;
  const int nrow = arg->nrow;
  const int ncol = arg->ncol;
  const int i1_block = arg->i1_block;
  const int ii_block = arg->ii_block;
  const int jj_block = arg->jj_block;

  SquareBlockMatrix<T> &D = *(arg->D);
  const int n0 = D.dim_kern_block(id_block);

  const T none(-1.0);
  const T one(1.0);

  if ((task_position % 2 == 1) && (id_level == 0)) {
    T *a = arg->a->addrCoefs(); //T *a = *(arg->a);
    D.copyBlockToArray(arg->ii_block, arg->jj_block, a, n);
    if (!arg->isSym) {
      D.copyBlockToArray(arg->jj_block, arg->ii_block, a, n);
    }
  }

  blas_gemm<T>(CblasTrans, CblasNoTrans,
	       nrow, 
	       ncol,
	       (nrow - n0), //
	       none, // alpha,
	       D.addrCoefBlock(ii_block, i1_block), // lower
	       D.nrowBlock(ii_block, i1_block), 
	       D.addrCoefBlock(i1_block, jj_block), // upper
	       D.nrowBlock(i1_block, jj_block), 
	       one, // beta,
	       D.addrCoefBlock(ii_block, jj_block), // Schur
	       D.nrowBlock(ii_block, jj_block));
  if (!arg->isSym) {
    // S(jj,ii)-=L(jj)*U(ii) <==> S(jj,ii)^T-=L(jj)^T*U(ii)^T
    blas_gemm<T>(CblasTrans, CblasNoTrans,
		 nrow,        // row of S^T
		 ncol,
		 (nrow - n0), //
		 none, // alpha,
		 D.addrCoefBlock(i1_block, ii_block),  // upper^T
		 D.nrowBlock(i1_block, ii_block),
		 D.addrCoefBlock(jj_block, i1_block),  // lower^T
		 D.nrowBlock(jj_block, i1_block),  
		 one, // beta,
		 D.addrCoefBlock(jj_block, ii_block), // Schur
		 D.nrowBlock(jj_block, ii_block));
  }
}
template
void C_dupdateb_Schur_offdiag<double>(void *arg_);

template
void C_dupdateb_Schur_offdiag<quadruple>(void *arg_);

template
void C_dupdateb_Schur_offdiag<complex<double> >(void *arg_);

template
void C_dupdateb_Schur_offdiag<complex<quadruple> >(void *arg_);

template
void C_dupdateb_Schur_offdiag<float>(void *arg_);

template
void C_dupdateb_Schur_offdiag<complex<float> >(void *arg_);
//

template<typename T, typename U>
void C_dfull_gauss_b(void *arg_)
{
  const T zero(0.0);
  C_dfull_gauss_arg<T, U> *arg = (C_dfull_gauss_arg<T, U> *)arg_;
  const int task_position = arg->task_position;
  const int id_block = arg->id_block;
  const int id_level = arg->id_level;
  const int n = arg->n;
  const bool verbose = arg->verbose;
  FILE *fp = *(arg->fp);

  //  int *nnn0 = arg->n0;
  int nn0;
  const int nrow = arg->nrow;
  const int i1_block = arg->i1_block;

  SquareBlockMatrix<T> &D = *(arg->D);
  const int i1 = D.IndexBlock(i1_block);

  int *permute_block = new int[D.nrowBlock(i1_block)]; // block_size()] ; 
  const int aug_dim = *(arg->aug_dim);
  const U eps_machine = *(arg->eps_machine);
  const double eps_piv = *(arg->eps_piv);
  
  // get smaller pivot of two children

  T *a_diag = D.addrCoefBlock(arg->i1_block, arg->i1_block);
  int nrow_block = D.nrowBlock(arg->i1_block, arg->i1_block);
  vector<int>& permute = D.getPermute();
  double pivot, pivot0, pivot1;
  vector<int>& singIdx0 = D.getSingIdx0();

  if (task_position % 2 == 1) {
    // copy the whole matrix for kernel detection and refactorization
    //    
    if (*(arg->kernel_detection) || (id_level == 0)) {  // 30 Jul.2014
      ColumnMatrix<T> &a = *(arg->a);
      a.init(n, n);
      // symmetrizie the first diagonal block : copy upper to lower
      if (arg->isSym) {
	for (int i = 0; i < nrow_block; i++) {
	  for (int j = 0; j < i; j++) { // access lower
	    a_diag[i + j * nrow_block] = a_diag[j + i * nrow_block];
	  }
	} 
      }
      D.copyBlockToArray(arg->i1_block, arg->i1_block, a.addrCoefs(), n);
      //      D.copyToArrayFull(a, n);
    } // if (id_level == 0)
    singIdx0.clear();
    singIdx0.reserve(n);
    //    singIdx0.resize(n, -1);
    D.set_dim_kern(0);    // counting total number of null pivots among blocks
    // initialize permute[]
    for (int i = 0; i < n; i++) {
      permute[i] = i;
    }
    pivot0 = *(arg->pivot0);
    pivot1 = *(arg->pivot1);
    pivot = (pivot0 < pivot1 ? pivot0 : pivot1);
    *(arg->quit) = false;  // default is to continue queue for refactorization
  } //  if (task_position % 2 == 1) 
  else {
    pivot = *(arg->pivot);
  }

  const int n0 = 0;
  double fop;

  if ((task_position % 2 == 1) && (nrow > aug_dim)) {
    bool flag_repeat_piv = true;
    //    bool repeat_first = true;
    double eps_piv1 = eps_piv;
    int count_repeat = 0;
    ColumnMatrix<T> a_diag1(nrow_block, nrow_block);
    while (flag_repeat_piv) {
      if (count_repeat == 0) {
	blas_copy<T>((nrow_block * nrow_block),
		     a_diag, 1, a_diag1.addrCoefs(), 1);
      }
      else {
	blas_copy<T>((nrow_block * nrow_block),
		     a_diag1.addrCoefs(), 1, a_diag, 1);
      } // if (count_repeat > 0)
      if (arg->isSym) {
	full_ldlt_permute<T, U>(&nn0, n0, nrow, a_diag, nrow_block, &pivot, 
				permute_block, 
				eps_piv1, &fop);
      } 
      else {
	full_ldu_permute<T, U>(&nn0, n0, nrow, a_diag, nrow_block, &pivot, 
			       permute_block, 
			       eps_piv1, &fop);
      } //   if (arg->isSym) 
      if (((nrow - nn0) >= aug_dim) || (eps_piv1 < TOL_PIVOT)) {
	flag_repeat_piv = false;
      }
      else {
	eps_piv1 /= 10.0;
	count_repeat++;
      }
    } // while (flag_repeat_piv)
    if (eps_piv1 < TOL_PIVOT) {
      D.set_pivrelaxed();
    }
    a_diag1.free();
    if (count_repeat > 0) {
      diss_printf(verbose, fp,
	      "%s %d : eps_piv = %g pivot = %g n0 = %d count_repeat = %d\n",
	      __FILE__, __LINE__, eps_piv1, pivot, n0, count_repeat);
    }    
  } //   if ((id_level == 0) && (task_position % 2 == 1)) 
  else {
    if (arg->isSym) {
      full_ldlt_permute<T, U>(&nn0, n0, nrow, a_diag, nrow_block, &pivot, 
			      permute_block, 
			      eps_piv, &fop);
    } //   if (arg->isSym) 
    else {
      full_ldu_permute<T, U>(&nn0, n0, nrow, a_diag, nrow_block, &pivot, 
			     permute_block, 
			     eps_piv, &fop);
    }
  }
  if (nn0 > 0) {
    diss_printf(verbose, fp,
		"%s %d : nd = %d : level = %d block = %d null = %d / %d\n",
		__FILE__, __LINE__, arg->nb,  arg->id_level, arg->id_block,
		nn0, nrow);
  }
  // permute_block is defined by Fortran array, i.e. takes index starting 1
  for (int i = i1; i < i1 + nrow; i++) {
    permute[i] = permute_block[i - i1] + i1;
  }
  // store singular nodes
  int sing_max = D.dim_kern(); // D.sing_max();
  {
    int itmp = i1 + nrow - nn0;
    for (int i = 0; i < nn0; i++, itmp++, sing_max++) {
      //      singIdx0[sing_max] = itmp;
      singIdx0.push_back(itmp);
    }
    D.set_dim_kern(sing_max);  // total number of null pivots
    D.set_dim_kern_block(id_block, nn0); // #null pivots candidate in the block
  }

  if (task_position / 2 == 1) {
    { // scope for list<int> tmp_idx
      list<int> tmp_idx;
      for (int i = 0; i < sing_max; i++) {
	for (int j = -1; j <= 1; j++) {
	  const int itmp = singIdx0[i] + j;
	  if (itmp >= 0 && itmp < n) {
	    tmp_idx.push_back(itmp);
	  }
	}
      }
      tmp_idx.sort();
      tmp_idx.unique();

      if (sing_max > 0) {
	diss_printf(verbose, fp,
		    "%s %d : %d : %g\n", __FILE__, __LINE__, arg->nb, pivot);
      }
      for (list<int>::const_iterator it = tmp_idx.begin(); 
	   it != tmp_idx.end(); ++it) {
	bool dispflag = false;
	for (int j = 0; j < sing_max; j++) {
	  if ((*it) == singIdx0[j]) {
	    dispflag = true;
	  }
	}
	diss_printf(verbose, fp, "%d %s %s\n", (*it), (dispflag ? "*" :":"),
		    tostring<T>(D.diag(*it)).c_str());
      }
    }
    // nullifying diagonal entries
    for (int m = 0; m < sing_max; m++) {
      const int ii = singIdx0[m];
      D.diag(ii) = zero;
    }
    bool refactorize = false;   // default value
    // check dim of factorized matrix
    int dim_kern = 0;
    bool flagKernelDetect = true;
    if (sing_max > 0) {
      diss_printf(verbose, fp,
		  "%s %d : C_dfull_gauss_b : id_level %d : nb = %d : %d / %d ",
		  __FILE__, __LINE__,
		  id_level, arg->nb, arg->id_block, arg->num_block);
      diss_printf(verbose, fp, "task_position = %d, i1 = %d , nrow = %d\n",
		  task_position, i1, nrow);
      diss_printf(verbose, fp, "%s %d : sing_max = %d : ",
		  __FILE__, __LINE__, sing_max);
      for (int i = 0; i < sing_max; i++) {
	diss_printf(verbose, fp, "%d ", singIdx0[i]);
      }
      diss_printf(verbose, fp, "\n");

      if ((task_position == 3) && (sing_max + aug_dim) > nrow) {
	// if (task_position == 3) {
	*(arg->quit) = true;  // refactorize is only applied for the root matrix
	diss_printf(verbose, fp,
		    "%s %d : nonsingular part is too small %d : %d - %d : %d\n",
		    __FILE__, __LINE__,
		    aug_dim, nrow, nn0, sing_max);
	flagKernelDetect = false;
	dim_kern = sing_max;       // nullify rows with all suspicious pivots
	D.set_KernelDetected(false);  // kernel is unknown
	singIdx0.resize(sing_max);
      }
    } // if (sing_max > 0)
    else {
      flagKernelDetect = false;
      singIdx0.resize(0);
      dim_kern = 0;
      D.set_KernelDetected(true);
    }
    if (flagKernelDetect) {
      if (*(arg->kernel_detection) || (id_level == 0)) {
	singIdx0.resize(sing_max);
	dim_kern = dimKernDense<T, U>(singIdx0, n,
				      aug_dim,
				      eps_machine,
				      eps_piv,
				      D,
				      arg->a->addrCoefs(), //*(arg->a), 
				      refactorize,
				      false,  // isFullPermute
				      arg->isSym,
				      verbose,
				      fp);
	diss_printf(verbose, fp, "%s %d : dim_kern = %d : singIdx0 = [",
		    __FILE__, __LINE__, dim_kern);
	for (int i = 0; i < singIdx0.size(); i++) {
	  diss_printf(verbose, fp, "%d ", singIdx0[i]);
	}
	diss_printf(verbose, fp, "]\n");
	if (dim_kern == (-1)) {
	  // refactorization happens in only the last level
	  refactorize = true;
	  D.set_KernelDetected(false);
	}
	else if (dim_kern != singIdx0.size()) {  // 20 Feb.2017
	  dim_kern = sing_max;       // nullify rows with all suspicious pivots
	  D.set_KernelDetected(false);  // kernel is unknown
	  singIdx0.resize(sing_max);
	}
	else {
	  D.set_KernelDetected(true);
	}
      }  //else if (*(arg->kernel_detection) || (id_level == 0)) {
      else {
	diss_printf(verbose, fp,
		    "%s %d : kernel detection will be done in the last\n",
		    __FILE__, __LINE__);
	dim_kern = sing_max;       // nullify rows with all suspicious pivots
	D.set_KernelDetected(false);  // kernel is unknown
	singIdx0.resize(sing_max);
      }  //else if (*(arg->kernel_detection) || (id_level == 0)) {
    } // if (flagKernelDetect)
     // if (sing_max > 0)
    if(arg->isSym) {
      D.freeLowerBlocks();
    }
    if (refactorize){
      //reset permutation for refactorization
      for (int i = 0; i < n; i++) {
	permute[i] = i;
      }
      diss_printf(verbose, fp,
		  "%s %d : refactorization with whole diagonal pivots starts\n",
		  __FILE__, __LINE__);
      *(arg->quit) = false;
    }
    else { // if (!refactorize)
      if (*(arg->kernel_detection) || (id_level == 0)) {  // 30 Jul.2014
	arg->a->free();
      }
      if(id_level == 0) {
	*(arg->quit) = true;  // refactorize is only applied for the root matrix
      }
      D.getSingIdx().resize(sing_max);
      for (int i = 0; i < sing_max; i++) {
	D.getSingIdx()[i] = permute[singIdx0[i]];
      }
      D.set_lastPivot(pivot);
      D.set_dim_kern(dim_kern);
      D.set_rank(n - dim_kern); // 
    } // if (refactorize)
  } // if (task_position / 2 == 1)
  else {
    if (nn0 > 0) {
      // nullifying diangol entries
      for (int ii = (i1 + nrow - nn0); ii < (i1 + nrow); ii++) {
	D.diag(ii) = zero;	//  a[ii * (n + 1)] = 0.0;
      }
    } // if (nn0 > 0)
  }  // if (task_position / 2 == 1)
  *(arg->pivot) = pivot;  // for passing information to other task
  delete [] permute_block;
}

template
void C_dfull_gauss_b<double, double>(void *arg_);

template
void C_dfull_gauss_b<complex<double>, double>(void *arg_);

template
void C_dfull_gauss_b<quadruple, quadruple>(void *arg_);

template
void C_dfull_gauss_b<complex<quadruple>, quadruple>(void *arg_);

template
void C_dfull_gauss_b<float, float>(void *arg_);

template
void C_dfull_gauss_b<complex<float>, float>(void *arg_);

//

template<typename T>
void C_dinvDL_timesU(void *arg_)
{
  const T one(1.0);
  C_dinvDL_timesU_arg<T> *arg = (C_dinvDL_timesU_arg<T> *)arg_;
  const int task_position = arg->task_position;
  const int id_level = arg->id_level;
  const int id_block = arg->id_block;
  const int n = arg->n;
  const int nrow = arg->nrow;
  const int ncol = arg->ncol;

  SquareBlockMatrix<T> &D = *(arg->D);
  const int i1 = D.IndexBlock(arg->i1_block);
  
  T *a_diag  = D.addrCoefBlock(arg->i1_block, arg->i1_block);
  T *a_upper = D.addrCoefBlock(arg->i1_block, arg->jj_block);
  const int nrow_block_diag = D.nrowBlock(arg->i1_block, arg->i1_block);
  const int nrow_block_upper= D.nrowBlock(arg->i1_block, arg->jj_block);
  if (arg->isSym) {
    D.allocateBlock(arg->jj_block, arg->i1_block);
  }

  T *a_lower = D.addrCoefBlock(arg->jj_block, arg->i1_block);
  const int nrow_block_lower = D.nrowBlock(arg->jj_block, arg->i1_block);
  
  vector<int>& permute = D.getPermute();
  const int n0 = D.dim_kern_block(id_block);
  // assuming to be allocated in the cache memory
  VectorArray<T> a_tmp(nrow);
  if ((task_position % 2 == 1) && (id_level == 0)) {
    //    T *a = *(arg->a);
    ColumnMatrix<T> &a = *(arg->a);
    // copy a_upper
    D.copyBlockToArray(arg->i1_block, arg->jj_block, a.addrCoefs(), n);  
    if (!arg->isSym) {
      // copy a_lower
      D.copyBlockToArray(arg->jj_block, arg->i1_block, a.addrCoefs(), n);  
    }
  }
  {
    // copy upper to upper with permutation
    for (int j = 0; j < ncol; j++) {
      const int jn = j * nrow;
      for (int i = 0; i < nrow; i++) {
	const int ip = permute[i1 + i] - i1; 
	a_tmp[i] = a_upper[ip + jn];
      }
      blas_copy<T>(nrow, a_tmp.addrCoefs(), 1, a_upper + jn, 1);
    }
  }
  if (!arg->isSym) {
    // copy transposed lower to transposed lower with permutation
    for (int j = 0; j < ncol; j++) {
      const int jn = j * nrow;
      for (int i = 0; i < nrow; i++) {
	const int ip = permute[i1 + i] - i1; 
	a_tmp[i] = a_lower[ip + jn];
      }
      blas_copy<T>(nrow, a_tmp.addrCoefs(), 1, a_lower + jn, 1);
    }
  }

  blas_trsm<T>(CblasLeft, CblasLower, CblasNoTrans, CblasUnit,
	       (nrow - n0),    // n0 : dim of suspicious pivots
	       ncol,   
	       one, //  alpha,
	       a_diag,
	       nrow_block_diag, // 3 Dec. 2015
	       a_upper, 
	       nrow_block_upper); // Dec. 2015
  if (arg->isSym) {
    for (int j = 0; j < ncol; j++) {
      const int jn = j * nrow;
      blas_copy<T>(nrow, a_upper + jn, 1, a_lower + jn, 1);
    }
  } 
  else {
    // (A_{jj i1} U_{i1 i1}^-1)^T = U_{i1 i1}^{-T} A_{jj i1}^{T}
    blas_trsm<T>(CblasLeft, CblasUpper, CblasTrans, CblasUnit,
		 (nrow - n0),      // n0 : dim of suspicious pivots
		 ncol,         
		 one, // alpha, 
		 a_diag,
		 nrow_block_diag, // 3 Dec.2015
		 a_lower,
		 nrow_block_lower); // 3 Dec.2015
  }
  for (int i = 0; i < nrow; i++) {
    const T aa = D.diag(i1 + i);
    for (int j = 0; j < ncol; j++) {
      a_upper[i + j * nrow] *= aa;
    }
  }
}

template
void C_dinvDL_timesU<double>(void *arg_);

template
void C_dinvDL_timesU<quadruple>(void *arg_);

template
void C_dinvDL_timesU<complex<double> >(void *arg_);

template
void C_dinvDL_timesU<complex<quadruple> >(void *arg_);

template
void C_dinvDL_timesU<float>(void *arg_);

template
void C_dinvDL_timesU<complex<float> >(void *arg_);

//

template<typename T, typename U>
void C_gauss_whole_pivot(void *arg_)
{
  const T zero(0.0);
  const T one(1.0);
  C_dfull_gauss_arg<T, U> *arg = (C_dfull_gauss_arg<T, U> *)arg_;

  const int task_position = arg->task_position;
  const int n = arg->n;
  SquareBlockMatrix<T> &D = *(arg->D);
  ColumnMatrix<T> &a = *(arg->a);
  const double eps_piv = *(arg->eps_piv);

  double pivot = *(arg->pivot);
  
  //  double *a_sym = D.addrCoefs();
  vector<int>& permute = D.getPermute();
  vector<int>& singIdx0 = D.getSingIdx0();
  const int aug_dim = *(arg->aug_dim);
  const U eps_machine = *(arg->eps_machine);
  const bool verbose = arg->verbose;
  FILE *fp = *(arg->fp);
  int n1, dim_kern;
 

  if (task_position % 2 == 1) {
    ColumnMatrix<T> aa(n, n);
    aa.copy(a);          //  blas_copy<T>((n * n), a.addrCoefs(), 1, aa, 1);
    if (arg->isSym) {
    // copy upper to lower for rank-1 update procedure by dsyr('L', ) 
      for (int j = 0; j < n; j++) {
	for (int i = 0; i < j; i++) {
	  //	  aa[j + i * n] = aa[i + j * n];
	  aa(j, i) = aa(i, j);
	}
      }
    }
    const double pivot0 = *(arg->pivot0);
    const double pivot1 = *(arg->pivot1);
    pivot = (pivot0 < pivot1 ? pivot0 : pivot1);
    int n0 = 0;
    diss_printf(verbose, fp,
		"%s %d : C_gauss_whole_pivot : serial factroization : n = %d\n",
	      __FILE__, __LINE__, n);
    int count_repeat = 0;
    bool flag_repeat_piv = true;
    double eps_piv1 = eps_piv;
    while (flag_repeat_piv) {
      if (count_repeat > 0) {
	aa.copy(a);
	//	blas_copy<T>((n * n), a.addrCoefs(), 1, aa.addrCoefs(), 1);
	if (arg->isSym) {
	  // copy upper to lower for rank-1 update procedure by dsyr('L', ) 
	  for (int j = 0; j < n; j++) {
	    for (int i = 0; i < j; i++) {
	      aa[j + i * n] = aa[i + j * n];
	    }
	  }
	}
      } // if (count_repeat > 0)
      n0 = 0;
      full_gauss3<T, U>(&n0, aa.addrCoefs(), n, &pivot, &permute[0], 
			arg->isSym,
			eps_piv1,
			verbose,
			fp);
      diss_printf(verbose, fp,
		  "%s %d : C_gauss_whole_pivot : pivot = %g n0 = %d\n",
		  __FILE__, __LINE__, pivot, n0);
      if (((n - n0) >= aug_dim) || (eps_piv1 < TOL_PIVOT)) {
	flag_repeat_piv = false;
      }
      else {
	eps_piv1 /= 10.0;
	count_repeat++;
      }
    } // while (flag_repeat_piv)
    diss_printf(verbose, fp,
		"%s %d : eps_pvi = %g pivot = %g n0 = %d count_repeat = %d\n",
		__FILE__, __LINE__, eps_piv1, pivot, n0, count_repeat);
    //    D.unsetBlocked(); // 
    D.copyFromArray(aa.addrCoefs(), n);
    // modification of lower blocks by removing D^-1 in the same mannar of
    // factorization by blocks
    if (!arg->isSym) {
      const int num_block = D.num_blocks();
      for (int k = 0; k < num_block; k++) {
	for (int m = (k + 1); m < num_block; m++) {
	  T *lower = D.addrCoefBlock(m, k);
	  const int nrow = D.nrowBlock(m, k);
	  const int ncol = D.ncolBlock(m, k);
	  for (int i = 0; i < nrow; i++) {
	    const T aa = one / D.addrCoefBlock(k, k)[i * (nrow + 1)];
	    for (int j = 0; j < ncol; j++) {
	      lower[i + j * nrow] *= aa;
	    }
	  }
	}
      }
    }
    //
    dim_kern = 0;
    if (n0 > 0) {
      singIdx0.resize(n0);
		      // singular nodes are continuously located from the last
      for (int i = 0; i < n0; i++) {  
	singIdx0[i] = n - n0 + i;
      }
      dim_kern = dimKernDense<T, U>(singIdx0, n,
				    aug_dim,
				    eps_machine,
				    eps_piv,
				    D, 
				    a.addrCoefs(), 
				    false, // refactorize
				    true, // isFullPermute
				    arg->isSym,
				    verbose,
				    fp);
      dim_kern = n0;
      diss_printf(verbose, fp,
		  "%s %d : C_gauss_whole_pivot : dim_kern = %d\n",
		  __FILE__, __LINE__, dim_kern);
      if (dim_kern == (-1)) {
	diss_printf(verbose, fp,
		    "%s %d : strict diagonal pivot is not enough!\n",
		    __FILE__, __LINE__);
	dim_kern = n0;
	D.set_KernelDetected(false);
      }
      else {
	n0 = dim_kern;
	D.set_KernelDetected(true);
      }
    } // if (n0 > 0)
    aa.free();   // need to move after copy to SquareBlockMatrix D(,)
    D.unsetBlocked(); //     
    if (n0 > 0) {
      n1 = n - n0; // n - dim_kern; ==> regular part needs to be restored
  // nullify lower part corresponding to the suspicious pivots
      for (int j = 0; j < n; j++) {
	for (int i = n1; i < n; i++){
	  a[i + j * n] = zero;
	}
      }
      if (!arg->isSym) {
	for (int j = 0; j < n; j++) {
	  for (int i = n1; i < n; i++){
	    a[j + i * n] = zero;
	  }
	}
      }
      D.getSingIdx().resize(n0);
      for (int i = 0; i < n0; i++) {
	D.getSingIdx()[i] = permute[singIdx0[i]];
      }
      //
      D.set_lastPivot(pivot);
      D.set_dim_kern(dim_kern);
      D.set_rank(n - dim_kern);
    } //     if (n0 > 0) 
    else {
      singIdx0.resize(0);
      D.getSingIdx().resize(0);
      D.set_KernelDetected(true);
      D.set_lastPivot(pivot);
      D.set_dim_kern(0); 
      D.set_rank(n);
    } //    if (n0 > 0)
    {
      int dim_kern1 = dim_kern;
      for (int k = (D.num_blocks() - 1); k >= 0 ; k--) {
	int nrow_local = D.nrowBlock(k);
	int dim_kern2 = nrow_local < dim_kern1 ? nrow_local : dim_kern1;
	D.set_dim_kern_block(k, dim_kern2);
	if (dim_kern1 > nrow_local) {
	  dim_kern1 -= nrow_local;
	}
	else {
	  dim_kern1 = 0;
	}
      }
    }
#ifdef DEBUG_MATRIX_DFULLLDLT
    cout << arg->task_name << "n= " << n << endl;
    for (int i = 0; i < n; i++) {
      cout << " " << permute[i];
    }
    cout << endl;
    for (int i = 0; i < ((n * (n + 1)) / 2); i++) {
      cout << " " << a_sym[i];
    }
    cout << endl;
#endif
  } // if (task_position % 2 == 1) 
  *(arg->quit) = true;
  arg->a->free();
#ifdef DEBUG_MEMORY_ALLOC
  cerr << "C_gauss_whole_pivot : memory deallocate." << endl;
  //
  cout << "pivot " << pivot << endl;
#endif
  *(arg->pivot) = pivot;  // for passing information to other task

}    

template
void C_gauss_whole_pivot<double, double>(void *arg_);

template
void C_gauss_whole_pivot<complex<double>, double>(void *arg_);

template
void C_gauss_whole_pivot<quadruple, quadruple>(void *arg_);

template
void C_gauss_whole_pivot<complex<quadruple>, quadruple>(void *arg_);

template
void C_gauss_whole_pivot<float, float>(void *arg_);

template
void C_gauss_whole_pivot<complex<float>, float>(void *arg_);

//

template<typename T>
void C_dupdateb_Schur_offdiag_t(void *arg_)
{
  C_dupdateb_Schur_arg<T> *arg =  (C_dupdateb_Schur_arg<T> *)arg_;

  const int nrow = arg->nrow;
  const int ncol = arg->ncol;
  const int b_size = arg->b_size;
  SquareBlockMatrix<T> &D = *(arg->D);
  
  const int i1 = D.IndexBlock(arg->i1_block); // arg->i1_block * SIZE_B1;
  const int ii = D.IndexBlock(arg->ii_block); // arg->ii_block * SIZE_B1;
  const int jj = D.IndexBlock(arg->jj_block); // arg->jj_block * SIZE_B1;
  FILE *fp = *(arg->fp);
  const bool verbose = arg->verbose;
  diss_printf(verbose, fp, 
	    "C_dupdate_Schur_offdiag_t : i1=%d ii=%d jj=%d nrow=%d ncol=%d b_size=%d\n",
	    i1, ii, jj, nrow, ncol, b_size);
  diss_printf(verbose, fp, "no computation\n");
}

template
void C_dupdateb_Schur_offdiag_t<double>(void *arg_);

template
void C_dupdateb_Schur_offdiag_t<quadruple>(void *arg_);

template
void C_dupdateb_Schur_offdiag_t<complex<double> >(void *arg_);

template
void C_dupdateb_Schur_offdiag_t<complex<quadruple> >(void *arg_);

template
void C_dupdateb_Schur_offdiag_t<float>(void *arg_);

template
void C_dupdateb_Schur_offdiag_t<complex<float> >(void *arg_);
//
