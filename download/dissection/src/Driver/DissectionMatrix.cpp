/*! \file   DissectionMatrix.hpp
    \brief  management of threads for factorization and Fw/Bw substitution
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

#include "Driver/DissectionMatrix.hpp"
#include "Algebra/SparseRenumbering.hpp"

// to_stirng is available in C++11 but NEC SXC++ does not support it
#include "Compiler/OptionLibrary.h"

#include <string>

template<typename T, typename U>
DissectionMatrix<T, U>::DissectionMatrix(Dissection::Tree *btree, 
					 const int nb, 
					 bool isSym,
					 const bool verbose,
					 FILE *fp) :
  _nb(nb), _diag(NULL), _lower(NULL), _upper(NULL), _isSym(isSym)
{
  const int level_last = btree->NumberOfLevels() - 1;
  _level = btree->nodeLayer(_nb);
  _nrow = btree->sizeOfDomain(_nb);
  _ncol_offdiag = btree->sizeOfFathersStrips(_nb);
  if (_nrow == 0) {
    _ncol_offdiag = 0;  // for safety
  }
  // no need to be allocated through the whole process
  //  
  //  _loc2glob_diag = btree->getDiagLoc2Glob(_nb);
  //  _loc2glob_offdiag = btree->getOffdiagLoc2Glob(_nb);

  _csr_diag = &(btree->getDiagCSR(_nb));
  _csr_offdiag = &(btree->getOffdiagCSR(_nb));
  _upper = new RectBlockMatrix<T>;
  _lower = new RectBlockMatrix<T>;
  // for unsymmetric matrix _lower block is stored in transposed way to use
  // the same strips as _upper block
                                  // bool later_allocation
   _localSchur = new SquareBlockMatrix<T>;

   if (_level == level_last) {
     _islast = true;
    _color_mask = new int[_nrow];
    _colors = getColorMaskCSR(_color_mask, _csr_diag, verbose,fp);
    _tridiag = new TridiagBlockMatrix<T,U>*[_colors];
    for (int i = 0; i < _colors; i++) {
      _tridiag[i] = new TridiagBlockMatrix<T,U>(_nrow, SIZE_B1, _isSym, nb,
						verbose, fp);
    }
    // the value of the last pivot and singluar nodes to be stored in _diag
    _diag = new SquareBlockMatrix<T>;
    _upper->init(_nrow, _ncol_offdiag, SIZE_B1, 0);
    _upper->allocate();      
    _lower->init(_nrow, _ncol_offdiag, SIZE_B1, 0);
    if (!_isSym) {
      _lower->allocate();
    }
    _localSchur->init(_ncol_offdiag, SIZE_B1, _isSym, 0); // 02 Feb.2014
    _alignedFather = false;
  }
  else {
    _islast = false;
    if (_nrow > 0) {
      _diag = new SquareBlockMatrix<T>(_nrow, SIZE_B1, _isSym);
    }
    else {
      _diag = new SquareBlockMatrix<T>(); // dummy
    }
    _tridiag = (TridiagBlockMatrix<T,U> **)NULL;
    if (_nb > 1) {
      int ntmp = 0;
      if (_nrow > 0) {
	const int b_id = btree->brotherIndex(_nb);
	const int ll = _level - 1;
	const Dissection::SetOfStrips& f0 =
	  btree->getFathersStrips(_nb)[ll];
	const Dissection::SetOfStrips& f1 =
	  btree->getFathersStrips(b_id)[ll];
	_alignedFather = false;
	if ((btree->sizeOfDomain(b_id) > 0) &&
	    (f0.numberOfStrips() == 1) && (f0.numberOfStrips() == 1)) {
	  Dissection::SetOfStrips::const_iterator it0 = f0.begin();
	  Dissection::SetOfStrips::const_iterator it1 = f1.begin();
	  if (((*it0).begin_src == (*it1).begin_src) &&
	      ((*it0).width == (*it1).width)) {
	    _alignedFather = true;
	    ntmp = f0.numberOfIndices();
	  }
	}
      }
      else {
	_alignedFather = false; // the other brother is treaed as nonaligned
      }
      if (!_alignedFather) {
      	fprintf(fp, "%s %d : %d : %s\n", __FILE__, __LINE__,
      		_nb, (_alignedFather ? "aligned" : "nonaligned"));
      }
      _upper->init(_nrow, _ncol_offdiag, SIZE_B1, ntmp);
      _upper->allocate();      
      _lower->init(_nrow, _ncol_offdiag, SIZE_B1, ntmp);
      if (!_isSym) {
	_lower->allocate();
      }
      _localSchur->init(_ncol_offdiag, SIZE_B1, _isSym, ntmp); // 02 Feb.2014
    }
  }
            // pointer to working array which is allocated and deallocated in
            // C_dfull_sym_gauss_b() 
  _factorize_LDLt = new ColumnMatrix<T>; //new T*;     
  //  _factorize_LDLt_diag = new ColumnMatrix<T>; // new T*;
  if (verbose) {
    fprintf(fp,
	    "%s %d : %d : nrow : %d : %d = %d + %d num_blocks : %d = %d + %d\n",
	    __FILE__,
	    __LINE__,
	    _nb,
	    _nrow,
	    _ncol_offdiag, 
	    _localSchur->dimension0(),
	    _localSchur->dimension1(),
	    _localSchur->num_blocks(),
	    _localSchur->num_blocks0(),
	    _localSchur->num_blocks1());
  }
}

template
DissectionMatrix<double, double>::
DissectionMatrix(Dissection::Tree *btree, 
		 const int nb, 
		 const bool isSym,
		 const bool verbose,
		 FILE *fp);

template
DissectionMatrix<complex<double>, double>::
DissectionMatrix(Dissection::Tree *btree, 
		 const int nb, 
		 const bool isSym,
		 const bool verbose,
		 FILE *fp);

template
DissectionMatrix<quadruple, quadruple>::
DissectionMatrix(Dissection::Tree *btree, 
		 const int nb, 
		 const bool isSym,
		 const bool verbose,
		 FILE *fp);

template
DissectionMatrix<complex<quadruple>, quadruple>::
DissectionMatrix(Dissection::Tree *btree, 
		 const int nb, 
		 const bool isSym,
		 const bool verbose,
		 FILE *fp);
//

template<typename T, typename U>
void DissectionMatrix<T, U>::C_SparseSymbFact_queue(vector<C_task*>& queue,
						    Dissection::Tree *btree,
						    const bool verbose,
						    FILE **fp)
{
  string task_name = "sparse_s_fact : " + to_string(_nb);
  long opsl = _csr_diag->n; // since factarization is based on skyline
                                      // complexity is a linear function of size
   C_SparseSymbFact_arg<T, U> *arg = 
     new C_SparseSymbFact_arg<T, U>(_tridiag,
				    _colors,
				    _color_mask,
				    btree->sizeOfDomain(_nb),
				    _csr_diag,
				    verbose,
				    fp);
   //				 
   *(arg->nopd) = opsl;  // very rough estimate, will be updated after run
   *(arg->ops_complexity) = opsl;   
   
   queue[0] = new C_task(C_SPARSESYMBFACT,
			 task_name,
			 (void *)arg,
			 C_SparseSymbFact<T, U>,
			 1, // atomic_size
			 0, // atomic_id
			 arg->ops_complexity    // ops_complexity
			 );
}

template
void DissectionMatrix<double, double>::
C_SparseSymbFact_queue(vector<C_task*>& queue,
		       Dissection::Tree *btree,
		       const bool verbose,
		       FILE **fp);


template
void DissectionMatrix<complex<double>, double>::
C_SparseSymbFact_queue(vector<C_task*>& queue,
		       Dissection::Tree *btree,
		       const bool verbose,
		       FILE **fp);

template
void DissectionMatrix<quadruple, quadruple>::
C_SparseSymbFact_queue(vector<C_task*>& queue,
		       Dissection::Tree *btree,
		       const bool verbose,
		       FILE **fp);

template
void DissectionMatrix<complex<quadruple>, quadruple>::
C_SparseSymbFact_queue(vector<C_task*>& queue,
		       Dissection::Tree *btree,
		       const bool verbose,
		       FILE **fp);
//

template<typename T, typename U>
void DissectionMatrix<T, U>::C_SparseNumFact_queue(vector<C_task*>& queue,
						   Dissection::Tree *btree,
						   int nnz,
						   T *coefs,
						   double *eps_pivot,
						   double *pivot,
						   bool *kernel_detection,
						   int *aug_dim,
						   U *eps_machine,
						   vector<C_task*>& task_q,
						   const bool verbose,
						   FILE **fp)
{
  string task_name = "n : " + to_string(_nb);

   //   char *task_name_cstr = new char[task_name.str().size() + 1];
   //   strcpy(task_name_cstr, task_name.str().c_str());

                                     // complexity is a linear function of size
   const int nrow = btree->sizeOfDomain(_nb);
   const int ncol = btree->sizeOfFathersStrips(_nb);
   C_SparseNumFact_arg<T, U> *arg = 
     new C_SparseNumFact_arg<T, U>(_tridiag, //_num_Sys,
				   _isSym,
				   _colors,
				   _color_mask,
				   nnz, //ptDA->nz(), // nnz
				   coefs, //ptDA->getCoef(), // *coefs
				   nrow,
				   ncol,
				   _csr_diag,
				   _csr_offdiag,
				   _diag,
				   eps_pivot,
				   pivot,
				   kernel_detection,
				   aug_dim,
				   eps_machine,
				   _localSchur,
				   verbose,
				   fp,
				   _nb);
   const long ncoll = ncol;
   const long nrowl = nrow;
   *(arg->nopd) = ncoll * ncoll * nrowl; // very rough estimate ,
   // will be updated after run
   queue[0] = new C_task(C_SPARSENUMFACT,
			 task_name,
			 (void *)arg,
			 C_SparseNumFact<T, U>,
			 1, // atomic_size
			 0, // atomic_id
			 ((C_SparseSymbFact_arg<T, U> *)(task_q[0]->func_arg))->nopd);
}

template
void DissectionMatrix<double, double>::
C_SparseNumFact_queue(vector<C_task*>& queue,
		      Dissection::Tree *btree,
		      int nnz,
		      double *coefs,
		      double *eps_pivot,
		      double *pivot,
		      bool *kernel_detection,
		      int *aug_dim,
		      double *eps_machine,
		      vector<C_task*>& task_q,
		      const bool verbose,
		      FILE **fp);

template
void DissectionMatrix<complex<double>, double>::
C_SparseNumFact_queue(vector<C_task*>& queue,
		      Dissection::Tree *btree,
		      int nnz,
		      complex<double> *coefs,
		      double *eps_pivot,
		      double *pivot,
		      bool *kernel_detection,
		      int *aug_dim,
		      double *eps_machine,
		      vector<C_task*>& task_q,
		      const bool verbose,
		      FILE **fp);

template
void DissectionMatrix<quadruple, quadruple>::
C_SparseNumFact_queue(vector<C_task*>& queue,
		      Dissection::Tree *btree,
		      int nnz,
		      quadruple *coefs,
		      double *eps_pivot,
		      double *pivot,
		      bool *kernel_detection,
		      int *aug_dim,
		      quadruple *eps_machine,
		      vector<C_task*>& task_q,
		      const bool verbose,
		      FILE **fp);

template
void DissectionMatrix<complex<quadruple>, quadruple>::
C_SparseNumFact_queue(vector<C_task*>& queue,
		      Dissection::Tree *btree,
		      int nnz,
		      complex<quadruple> *coefs,
		      double *eps_pivot,
		      double *pivot,
		      bool *kernel_detection,
		      int *aug_dim,
		      quadruple *eps_machine,
		      vector<C_task*>& task_q,
		      const bool verbose,
		      FILE **fp);
//

template<typename T, typename U>
void DissectionMatrix<T, U>::C_SparseLocalSchur_queue(vector<C_task*>& queue,
						      Dissection::Tree *btree,
						      int nnz,
						      T *coefs,
						      vector<C_task*>& task_p,
						      const bool verbose,
						      FILE **fp)
{
  string task_name = "o : " + to_string(_nb);
   //   char *task_name_cstr = new char[task_name.str().size() + 1];
   //   strcpy(task_name_cstr, task_name.str().c_str());

                                     // complexity is a linear function of size
   C_SparseNumFact_arg<T, U> *arg = 
     new C_SparseNumFact_arg<T, U>(_tridiag,
				   _isSym,
				   _colors,
				   _color_mask,
				   nnz, //	   ptDA->nz(), // nnz
				   coefs, // ptDA->getCoef(), // *coefs
				   btree->sizeOfDomain(_nb),
				   btree->sizeOfFathersStrips(_nb),
				   _csr_diag,
				   _csr_offdiag,
				   _diag,
				   (double *)NULL, // eps_pivot
				   (double *)NULL, // pivot
				   (bool *)NULL, // dummy
				   (int *)NULL, // dummy 
				   (U *)NULL, // dummy
				   _localSchur,
				   verbose,
				   fp,
				   _nb);
   
   queue[0] = new C_task(C_SPARSESCHUR,
			 task_name,
			 (void *)arg,
			 C_SparseLocalSchur<T, U>,
			 1, // atomic_size
			 0, // atomic_id
			 ((C_SparseNumFact_arg<T, U> *)(task_p[0]->func_arg))->nopd
			 );
   queue[0]->fp = &(arg->fp);
   queue[0]->parents->push_back(task_p[0]);
}

template
void DissectionMatrix<double, double>::
C_SparseLocalSchur_queue(vector<C_task*>& queue,
			 Dissection::Tree *btree,
			 int nnz,
			 double *coefs,
			 vector<C_task*>& task_p,
			 const bool verbose,
			 FILE **fp);

template
void DissectionMatrix<complex<double>,  double>::
C_SparseLocalSchur_queue(vector<C_task*>& queue,
			 Dissection::Tree *btree,
			 int nnz,
			 complex<double> *coefs,
			 vector<C_task*>& task_p,
			 const bool verbose,
			 FILE **fp);

template
void DissectionMatrix<quadruple, quadruple>::
C_SparseLocalSchur_queue(vector<C_task*>& queue,
			 Dissection::Tree *btree,
			 int nnz,
			 quadruple *coefs,
			 vector<C_task*>& task_p,
			 const bool verbose,
			 FILE **fp);

template
void DissectionMatrix<complex<quadruple>, quadruple>::
C_SparseLocalSchur_queue(vector<C_task*>& queue,
			 Dissection::Tree *btree,
			 int nnz,
			 complex<quadruple> *coefs,
			 vector<C_task*>& task_p,
			 const bool verbose,
			 FILE **fp);
//

template<typename T, typename U>
void DissectionMatrix<T, U>::C_FillMatrix_queue(vector<C_task*>& queue,
						int nnz,
						T *coefs,
						const bool verbose,
						FILE **fp)
{
  {
    string task_name = "l : " + to_string(_nb);
    // 16 Sep.2014 : Atsushi
    // -1L is used for initialization of the array <=> dependecy 
    long ops;
    if (_level > 0) {
      // replace -1L by 1L : 29 Nov.2016 Atsushi
      //      ops = _csr_offdiag->nnz == 0 ? (-1L) : (long)_csr_offdiag->nnz;
      ops = _csr_offdiag->nnz == 0 ? 1L : (long)_csr_offdiag->nnz;
    }
    else {
      ops = 0L;
    }
    //  if (ops > 0) need to be skippid : 29 Aug.2014 Atsushi
    C_FillMatrix_arg<T> *arg = 
      new C_FillMatrix_arg<T>(_isSym,
			      (SquareBlockMatrix<T>*)NULL,
			      _upper,
			      (_isSym ? (RectBlockMatrix<T>*)NULL : _lower),
			      (CSR_indirect*)NULL,
			      _csr_offdiag,
			      coefs,
			      verbose,
			      fp,
			      _nb);
    *(arg->ops_complexity) = ops;  
 
    queue[1] = new C_task(C_FILLMATRIX,
			  task_name,
			  (void *)arg,
			  C_FillMatrix_offdiag<T>,
			  1, // atomic_size
			  0, // atomic_id
			  arg->ops_complexity // ops_complexity 
			  );
  }
  { // begin scope : task_name

    string task_name = "k : " + to_string(_nb);

    long ops = (long)(_csr_diag->nnz + _csr_diag->n) / 2L;
    C_FillMatrix_arg<T> *arg = 
      new C_FillMatrix_arg<T>(_isSym,
			      _diag,
			      (RectBlockMatrix<T>*)NULL, // upper
			      (RectBlockMatrix<T>*)NULL, // lower
			      _csr_diag,
			      (CSR_indirect*)NULL,
			      coefs,
			      verbose,
			      fp,
			   _nb);
    // ptDA->getCoef());
    *(arg->ops_complexity) = ops;   

    queue[0] = new C_task(C_FILLMATRIX,
			  task_name,
			  (void *)arg,
			  C_FillMatrix_diag<T>,
			  1, // atomic_size
			  0, // atomic_id
			  arg->ops_complexity // ops_complexity
			  );
  }
}

template
void DissectionMatrix<double, double>::
C_FillMatrix_queue(vector<C_task*>& queue,
		   int nnz,
		   double *coefs,
		   const bool verbose,
		   FILE **fp);

template
void DissectionMatrix<complex<double>, double>::
C_FillMatrix_queue(vector<C_task*>& queue,
		   int nnz,
		   complex<double> *coefs,
		   const bool verbose,
		   FILE **fp);

template
void DissectionMatrix<quadruple, quadruple>::
C_FillMatrix_queue(vector<C_task*>& queue,
		   int nnz,
		   quadruple *coefs,
		   const bool verbose,
		   FILE **fp);

template
void DissectionMatrix<complex<quadruple>, quadruple>::
C_FillMatrix_queue(vector<C_task*>& queue,
		   int nnz,
		   complex<quadruple> *coefs,
		   const bool verbose,
		   FILE **fp);
//

template<typename T, typename U>
int DissectionMatrix<T, U>::C_DFullLDLt_queue(vector<C_task*>& queue,
					      vector<int>& task_indcol,
					      vector<int>& task_ptr,
					      double *eps_piv,
					      bool *kernel_detection,
					      int *aug_dim,
					      U *eps_machine,
					      double *pivot,
					      double *pivot0,
					      double *pivot1,
					      vector<C_task*>& task_p,
					      const bool isChldrnAlgnd,
					      const bool verbose,
					      FILE **fp)
{
  const int n = _diag->dimension();       // _nrow
  const int upper_ncol = getUpperNCol();  // _ncol_offdiag

  const int num_block = _diag->num_blocks();
  int num_tasks1;
  int ipos, jpos;
  int parent_id;

  if (_level == 0) {
    num_tasks1 = ((num_block - 1) * num_block * (num_block + 1) / 6 +
		  num_block);
  }
  else {
    num_tasks1 = 0;
  }

  if (n == 0) {
    queue.resize(1);
    string task_name = ("a dummy : " +
			to_string(_level) + " : " +to_string(_nb));
    C_dummy_arg *arg = new C_dummy_arg(verbose, fp, _nb);
    //    *(arg->ops_complexity) = (-1L);
    queue[0] = new C_task(C_DUMMY,
			  task_name,
			  (void *)arg,
			  C_dummy,
			  1, // atomic_size,
			  0, // atomic_id,
			  arg->ops_complexity);
    queue[0]->parents->clear();
    task_ptr.clear();
    return 1;
  }
  int atomic_size, atomic_id;
  //  int *ptr;
  if (_level == 0) {
    task_ptr.resize(2 * num_block + 1);
    //    ptr = new int[2 * num_block + 1];
  }
  else {
    task_ptr.resize(num_block + 1);
  }
#if 0 // verbose
  if (queue.size() != (num_tasks0 + num_tasks1)) {
    cout << "queue.size() is incorrect : " << queue.size()
	 << " : " << (num_tasks0 + num_tasks1) << endl;
  }
#endif
#ifdef DEBUG_MUTEX_C_DFULLLDLT
  cout << "**** C_DFullLDLt " 
       << " starts with " << n << " sized-matrix by " << SIZE_B1 
       << " ****" << endl;
#endif  
  // keeping starting index of triangular factorizations at k-th level
  task_ptr[0] = 0;
  for (int k = 0; k < num_block; k++) {
    task_ptr[k + 1] = task_ptr[k] + ((num_block - k) * (num_block - k + 1)) / 2;
  }
  if (_level == 0) {
    int kk = num_block;
    for (int k = 0; k < num_block; k++, kk++) {
      task_ptr[kk + 1] = (task_ptr[kk] + 
			  ((num_block - k - 1) * (num_block - k)) / 2 + 1);
    }
  }
#ifdef DEBUG_C_DFULLLDLT
  cout << "group_id = " << group_id << " : ";
  for (int k = 0; k < num_block + 1; k++) {
    cout <<  task_ptr[k] << " ";
  }
  cout << endl;
#endif
  // queue with "blocked" pivot strategy
  // generating tasks 
  // a^0, {b^0_0, c^0_00, a^1}, b^0_1,....,b^0_m, c^0_01,c^0_11,...,c^0_mm, 
  // {b^1_0,c^1_00 a^2}, b^1_1,....,b^1_{m-1}, c^0_01,...,c^0_{m-1}{m-1}, 
  ipos = (-1);
  //
  jpos = 0;
  for (int k = 0; k < num_block; k++) {
    const int nrow = _diag->nrowBlock(k);
    int task_position;
    task_position = 0;
    
    if (k == 0) {
      task_position = 1;
    }
    if (k == (num_block - 1)) {
      task_position += 2;
    }

    // task a^k
    { // scope for char *task_name_cstr
      if (k == 0) {
	ipos = 0;
	atomic_size = 1; 		
	atomic_id = 0;
      }
      else {
	ipos = task_ptr[k - 1] + 3;
	atomic_size = 3;
	atomic_id = 2;
      }
      string task_name = ("a " + to_string(k) + " : " + to_string(_level) +
			  " : " + to_string(_nb));
      const long nrowl = (long)nrow;
      const long ops = nrowl * (nrowl + 1L) * (nrowl + 2L) + nrowl;
      C_dfull_gauss_arg<T, U> *arg = 
	new C_dfull_gauss_arg<T, U>(_isSym,
				    task_position,   // int task_position
				    _level,
				    num_block,
				    k,
				    _diag, 
				    _localSchur,
				    _factorize_LDLt, // pointer for T[n * n]
			    //	    _factorize_LDLt_diag,  // diag
				    (int *)NULL,     // int *permute_block
				    n,               // int n
				    upper_ncol,      // lower matrix size
				    nrow,            // int nrow
				    k,              // int i1
				    eps_piv,
				    pivot,
				    pivot0,
				    pivot1,
				    kernel_detection,
				    aug_dim,
				    eps_machine,
				    verbose,
				    fp,
				    _nb);
      *(arg->ops_complexity) = ops;   
      queue[ipos] = new C_task(C_DFULL_SYM_GAUSS,
			      task_name,
			       (void *)arg,
			       C_dfull_gauss_b<T, U>,
			       atomic_size,
			       atomic_id,
			       arg->ops_complexity // ops_complexity
			       );  
      arg->quit = &queue[ipos]->quit_queue;
      queue[ipos]->fp = &(arg->fp);
      task_indcol[jpos++] = ipos;
      if (k == 0) {
	// a^0 is the root
	queue[ipos]->parallel_max = 1;
	queue[ipos]->parallel_id = 0;
	// depndency on dgemm/dsub of previous factorization level
	queue[ipos]->parents->push_back(task_p[0]);
      }
      else {
	queue[ipos]->parallel_max = (-1);
	queue[ipos]->parallel_id = 0;
	// no need to write dependency : sequential operation inside of atom
	//                               resolves dependency
	//	parent_id = task_ptr[k - 1] + 2;    // c^{k-1}_{0,0}
	//	queue[ipos]->parents->push_back(queue[parent_id]);
      }
      if (task_position / 2 == 1) {
	queue[ipos]->to_next_task = num_tasks1;
      }
    }  // scope for char *task_namefp_cstr
    for (int j = 0; j < (num_block - 1 - k); j++) {
      const int jj = k + j + 1;
      // task b^k_j
      {  // scope for char *task_name_cstr
	                                     // lower part stored in transposed
	const int ncol = _diag->ncolBlock(jj, k); 
	if (j == 0) {
	  // a^{k}, b^{k}_0 
	  ipos = task_ptr[k] + 1; 
	  atomic_size = 3;
	  atomic_id = 0;
	}
	else {
	  // a^{k}, b^{k}_0 c^{k}_00 a^{k+1}, b^{k}_1, ...  
	  //                   three tasks a^{k} c^{k}_00 a^{k+1} prior {b^{k}}
	  ipos = task_ptr[k] + j + 3; 
	  atomic_size = 1;
	  atomic_id = 0;
	}
	string task_name = ("b " + to_string(k) + " " + to_string(j) 
			    + " : " + to_string(_level) + " : "
			    + to_string(_nb));
	const long nrowl = (long)nrow;
	const long ncoll = (long)ncol;
	const long ops = (nrowl - 1L) * nrowl * ncoll + nrowl * ncoll;
	// dtrsm for (nrow x ncol matrix) + D^-1 (nrow x ncol matrix)
	C_dinvDL_timesU_arg<T> *arg = 
	  new C_dinvDL_timesU_arg<T>(_isSym,
				     task_position,  
				     _level,
				     num_block,
				     k,
				     _diag, //D,
				     _factorize_LDLt,
				     n,                   // int n
				     nrow,                // int nrow
				     ncol,                // int ncol
				     k,                  // int i1
				     jj,                  // k + j + 1,
				     verbose,
				     fp,
				     _nb);
	*(arg->ops_complexity) = ops;   
	queue[ipos] = new C_task(C_DINV_DL_TIMESU,
				 task_name,
				 (void *)arg,
				 C_dinvDL_timesU<T>,
				 atomic_size,   
				 atomic_id,  
				 arg->ops_complexity // ops_complexity
				 );
	task_indcol[jpos++] = ipos;	
	queue[ipos]->parallel_max = (num_block - k - 1);
	queue[ipos]->parallel_id = j;
	
	if (k == 0) {
	 	// depndency on dgemm/dsub of previous factorization level
	  if (_isSym || !isChldrnAlgnd) {
	    const int jtmp = ((j + 1) * (j + 2)) / 2;
	    queue[ipos]->parents->push_back(task_p[jtmp]);
	  }
	  else {
	    const int jtmp = (j + 1) * (j + 1);
	    queue[ipos]->parents->push_back(task_p[jtmp]);
	    queue[ipos]->parents->push_back(task_p[jtmp + 1]);
	  }
	  parent_id = task_ptr[k];                    // a^0
	  queue[ipos]->parents->push_back(queue[parent_id]);        
	}
	else {
	  parent_id = task_ptr[k - 1] + 3;             // a^k
	  queue[ipos]->parents->push_back(queue[parent_id]); 
				                             // c^{k-1}_{0,j+1}
	  parent_id = (task_ptr[k - 1] + (num_block - k + 2) +
		       ((j + 1) * (j + 2)) / 2);
	  queue[ipos]->parents->push_back(queue[parent_id]);
	}  // if (k == 0)
      }    // scope for char *task_name_cstr : task b^k_j
    }
    for (int j = 0; j < (num_block - 1 - k); j++) {
      const int jj = k + j + 1;
      const int ncol = _diag->ncolBlock(jj, k); 
      for (int i = 0; i < j; i++) {
	// task c^k_{i,j}
	{
	  ipos = (task_ptr[k] + 
		  2 +                        // a^0 a^1
		  (num_block - k - 1) +      // b^0_0 ,... b^0_{num_blok-k-2}
		  (j * (j + 1)) / 2 +        // size of upper triangle 
		  i);                        //
	  string task_name = ("c " + to_string(k) + " " + to_string(i)
			      + " " + to_string(j) + " : " + to_string(_level)
			      + " : " + to_string(_nb));
	  const long nrowl = (long)nrow;
	  const long ncoll = (long)ncol;
	  const long ops = nrowl * nrowl * ncoll * 2L;
	  // optimized BLAS3 gemm : opt. factor = 1/0.5
	  C_dupdateb_Schur_arg<T> *arg = 
	    new  C_dupdateb_Schur_arg<T>(_isSym,
					 task_position,
					 _level,
					 num_block,
					 k,				
					 _diag, //D,
					 _factorize_LDLt,
					 // _factorize_LDLt_work,
					 n, 
					 nrow, 
					 ncol,
					 SIZE_B1,
					 k,                 // i1
					 k + i + 1, 
					 jj, // k + j + 1,
					 verbose,
					 fp,
					 _nb);
	  *(arg->ops_complexity) = ops;   
	  queue[ipos] = new C_task(C_DHALF_SCHUR_B,
				   task_name,
				   //			   task_name_cstr,
				   (void *)arg,
				   C_dupdateb_Schur_offdiag<T>,
				   1, // atomic_size
				   0, // atomic_id
				   arg->ops_complexity // ops_complexity
				   );
	  task_indcol[jpos++] = ipos;
	  //
	  // c^k_00 is computed with b^k_0 and a^{k+1}
	  queue[ipos]->parallel_max = ((num_block - k - 1) *
				       (num_block - k)) / 2 - 1;
	  queue[ipos]->parallel_id = (j * (j + 1)) / 2 + i - 1 ;
#if 0
	  fprintf(stderr, "%s %d : %d/", __FILE__, __LINE__, (int)task_p.size());
	  for (vector<C_task *>::iterator jt = task_p.begin();
	       jt != task_p.end(); ++jt) {
	    fprintf(stderr, "%s/", (*jt)->task_name);
	  }
	  fprintf(stderr, "\n");
#endif
	  if (k == 0) {
	    // depndency on dgemm/dsub of previous factorization level
	    // (i + 1)-th row, (j + 1)-th column
	    if (_isSym || !isChldrnAlgnd) {
	      const int jtmp = ((j + 1) * (j + 2)) / 2 + i + 1;
	      queue[ipos]->parents->push_back(task_p[jtmp]);
	    } 
	    else {
	      const int jtmp = ((j + 1) * (j + 1)) + 2 * i + 2;
	      queue[ipos]->parents->push_back(task_p[jtmp]);
	      queue[ipos]->parents->push_back(task_p[jtmp + 1]);
	    }
	  }
	  else {
	    //	  if (k > 0) { 
	    // c^{k-1}_{i+1,j+1}
	    parent_id = (task_ptr[k - 1] + (num_block - k + 2) +
			 ((j + 1) * (j + 2)) / 2 + (i + 1));
	    queue[ipos]->parents->push_back(queue[parent_id]);
	  }
	  if (i == 0) {
	    parent_id = task_ptr[k] + 1;              // b^k_i
	    queue[ipos]->parents->push_back(queue[parent_id]); 
	  }
	  else {
	    parent_id = task_ptr[k] + 3 + i;
	    queue[ipos]->parents->push_back(queue[parent_id]);
	  } // if (i == 0)
	    // by definition of the loop, i < j
	  parent_id = task_ptr[k] + 3 + j;       	  // b^k_j
	  queue[ipos]->parents->push_back(queue[parent_id]); 
	}  //  scope for char *task_name_cstr :: task c^k_{i,j}
      } // loop : i

      // task c^k_{j,j}      
      {
	if (j == 0) {
	  ipos = task_ptr[k] + 2; 
	  atomic_size = 3;
	  atomic_id = 1;
	}
	else {
	  ipos = (task_ptr[k] + 
		  2 +                        // a^0 a^1
		  (num_block - k -1) +       // b^0_0 ,... b^0_{num_blok-k-2}
		  (j * (j + 1)) / 2 +        // size of upper triangle 
		  j);        
	  atomic_size = 1;
	  atomic_id = 0;
	}
	string task_name = ("c " + to_string(k) + " " + to_string(j) + " "
			    + to_string(j) +  " : " + to_string(_level)
			    + " : " + to_string(_nb));
	//	char *task_name_cstr = new char[task_name.str().size() + 1];
	//	strcpy(task_name_cstr, task_name.str().c_str());
	const long nrowl = (long)nrow;
	const long ncoll = (long)ncol;
	const long ops = nrowl * ncoll * (ncoll + 1L);
	// non-optimized gemm' : opt. factor = 1 
	C_dupdateb_Schur_arg<T> *arg = 
	  new C_dupdateb_Schur_arg<T>(_isSym,
				      task_position,
				      _level,
				      num_block,
				      k,
				      _diag, // D,
				      _factorize_LDLt,
				      // _factorize_LDLt_work,
				      n, 
				      nrow,                // int nrow
				      ncol,                // int ncol
				      SIZE_B1,
				      k, 
				      k + j + 1, 
				      (-1),
				      verbose,
				      fp,
				      _nb);               // jj
	*(arg->ops_complexity) = ops;   
	
	queue[ipos] = new C_task(C_DHALF_SCHUR_B,
				 task_name,
				 (void *)arg,
				 C_dupdateb_Schur_diag<T>,
				 atomic_size,
				 atomic_id,
				 arg->ops_complexity // ops_complexity
				 );
	task_indcol[jpos++] = ipos;
	if (j > 0) {
	  // c^k_00 is computed with b^k_0 and a^{k+1}
	  queue[ipos]->parallel_max = ((num_block - k - 1) * 
				       (num_block - k)) / 2 - 1;
	  queue[ipos]->parallel_id = (j * (j + 1)) / 2 + j - 1 ;
	}
	else {
	  queue[ipos]->parallel_max = (-1);
	  queue[ipos]->parallel_id = 0;
	}
	if (k == 0) {
	  // depndency on dgemm/dsub of previous factorization level
	  // (j + 1)-th row, (j + 1)-th column
	  if (_isSym || !isChldrnAlgnd) {
	    const int jtmp = ((j + 1) * (j + 2)) / 2 + j + 1;
	    queue[ipos]->parents->push_back(task_p[jtmp]);
	  }
	  else {
	    const int jtmp = (j + 1) * (j + 1) + 2 * j + 2;
	    queue[ipos]->parents->push_back(task_p[jtmp]);
	  }
	}
	else {
	  //	if (k > 0) {
	  // c^{k-1}_{j+1,j+1}
	  parent_id = (task_ptr[k - 1] + (num_block - k + 2) +      
		       ((j + 1) * (j + 2)) / 2 + (j + 1));
	  queue[ipos]->parents->push_back(queue[parent_id]);
	}
	// else {  depndency to dgemm of previous factorization level }
	if (j > 0) {
	  parent_id = task_ptr[k] + 3 + j;           // b^k_j
	  queue[ipos]->parents->push_back(queue[parent_id]);
	}
	// no need to write dependency : sequential operation inside of atom
	//                               resolves dependency
	// else { 
        //   parent_id = task_ptr[k] + 1;               // b^k_0
	//  queue[ipos]->parents->push_back(queue[parent_id]);  	
	//}
      }    // scope for char *task_name_cstr : c^k_{j,j}
    } // loop : j
  } // loop : k

  if (_level == 0) { // 22 Jan.2013 : Atsushi : without refact. for level > 0
  // queue for refactorization with full pivot strategy
    ipos = task_ptr[num_block];
    for (int k = 0; k < num_block; k++) {
      const int nrow = _diag->nrowBlock(k);
      int task_position;
      task_position = 0;
      
      if (k == 0) {
	task_position = 1;
      }
      if (k == (num_block - 1)) {
	task_position += 2;
      }
      // task a^k
      { // scope for char *task_name_cstr
	//  pv_param.pivot is updating during factorizaion
	string task_name = ("A " + to_string(k) + " : " +to_string(_level)
			    + " : " + to_string(_nb));
	const long nrowl = (long)nrow;
	const long ops = nrowl * (nrowl + 1L) * (nrowl + 2L) + nrowl;
	C_dfull_gauss_arg<T, U> *arg = 
	  new C_dfull_gauss_arg<T, U>(_isSym, 
				      task_position,   // int task_position
				      _level,
				      num_block,
				      k,
				      _diag, // D, 
				      //		   _lower,
				      _localSchur,
				      _factorize_LDLt,
				      //_factorize_LDLt_diag,  // diag
				      (int *)NULL,     // int *permute_block
				      n,               // int n
				      upper_ncol,      // lower matrix size
				      nrow,            // int nrow
				      k, //kk,              // int i1
				      eps_piv,
				      pivot,
				      pivot0,
				      pivot1,
				      kernel_detection,
				      aug_dim,
				      eps_machine,
				      verbose,
				      fp,
				      _nb);
	*(arg->ops_complexity) = ops;   
	queue[ipos] = new C_task(C_DFULL_SYM_GAUSS,
				 task_name,
				 //				 task_name_cstr,
				 (void *)arg,
				 C_gauss_whole_pivot<T, U>,
				 atomic_size,
				 atomic_id,
				 arg->ops_complexity // ops_complexity
				 );  
	arg->quit = &queue[ipos]->quit_queue;
	queue[ipos]->fp = &(arg->fp);
	task_indcol[jpos++] = ipos;
	queue[ipos]->parallel_max = 1;
	queue[ipos]->parallel_id = 0;
	queue[ipos]->to_next_task = (((num_block - 1 - k) * (num_block - k) * 
				      (num_block + 1 - k)) / 6 +
				     num_block - k - 1);
	if (k == 0) {
	  const int jj = task_ptr[num_block - 1];
	  queue[ipos]->parents->push_back(queue[jj]); // the last 
	}
	else {
         // depending on all C_dupdateb_Schur_offdiagt()s in the previous level
	  const int jbegin = task_ptr[(k - 1) + num_block] + 1;
	  const int jend = task_ptr[k + num_block];
	  for (int j = jbegin; j < jend; j++) {
	    queue[ipos]->parents->push_back(queue[j]);
	  }
	}
	ipos++;
      }  // scope for char *task_name_cstr
      int itmp = 0;
      for (int j = 0; j < (num_block - 1 - k); j++) {
	const int jj = k + j + 1;
	//   const int nrow1 = (j == (num_block - 2 - k)) ? size_res : SIZE_B1;
	const int nrow1 = _diag->nrowBlock(jj);
	for (int i = j; i < (num_block - 1 - k); i++) {
	  const int ii = k + i + 1;
	  // const int ncol1 = (i == (num_block - 2 - k)) ? size_res : SIZE_B1;
	  const int ncol1 = _diag->nrowBlock(ii);
	  // task c^k_{i,j}
	  {
	    string task_name = ("C " + to_string(k) + " " + to_string(i)
	      + " " + to_string(j) + " : " + to_string(_level)
	      + " : " + to_string(_nb));
	    const long nrowl = (long)nrow1;
	    const long ncoll = (long)ncol1;
	    const long ops = nrowl * nrowl * ncoll * 2L;
	    // optimized BLAS3 gemm : opt. factor = 1/0.5
	    C_dupdateb_Schur_arg<T> *arg = 
	      new C_dupdateb_Schur_arg<T>(_isSym,
					  task_position,
					  _level,
					  num_block,
					  k,
					  _diag, // D,
					  _factorize_LDLt,
					  n, 
					  nrow1, 
					  ncol1,
					  SIZE_B1,
					  k,  //kk, 
					  ii, // k + i + 1
					  jj, // k + j + 1,
					  verbose,
					  fp,
					  _nb);
	    *(arg->ops_complexity) = ops;   
	    queue[ipos] = new C_task(C_DHALF_SCHUR_BT,
				     task_name,
				     (void *)arg,
				     C_dupdateb_Schur_offdiag_t<T>,
				     1, // atomic_size
				     0, // atomic_id
				     arg->ops_complexity // ops_complexity
				     );
	    task_indcol[jpos++] = ipos;
	    //
	    // c^k_00 is computed with b^k_0 and a^{k+1}
	    queue[ipos]->parallel_max = ((num_block - k - 1) *
					 (num_block - k)) / 2;
	    queue[ipos]->parallel_id = itmp++;   // itmp = i + j * n_block - ...
	    const int jj = task_ptr[k + num_block];
	    queue[ipos]->parents->push_back(queue[jj]);
	    ipos++;
	  }
	}   // loop : j
      } // loop : i
    } // loop : k
    //  nop += (double)n * (double)n * (double)n / 3.0; 
  } // if (_level == 0)
  //
  return ipos;
}

template
int DissectionMatrix<double, double>::
C_DFullLDLt_queue(vector<C_task*>& queue,
		  vector<int>& task_indcol,
		  vector<int>& task_ptr,
		  double *eps_piv,
		  bool *kernel_detection,
		  int *aug_dim,
		  double *eps_machine,
		  double *pivot,
		  double *pivot0,
		  double *pivot1,
		  vector<C_task*>& task_p,
	          const bool isChldrnAlgnd,
		  const bool verbose,
		  FILE **fp);

template
int DissectionMatrix<complex<double>, double>::
C_DFullLDLt_queue(vector<C_task*>& queue,
		  vector<int>& task_indcol,
		  vector<int>& task_ptr,
		  double *eps_piv,
		  bool *kernel_detection,
		  int *aug_dim,
		  double *eps_machine,
		  double *pivot,
		  double *pivot0,
		  double *pivot1,
		  vector<C_task*>& task_p,
	          const bool isChldrnAlgnd,
		  const bool verbose,
		  FILE **fp);

template
int DissectionMatrix<quadruple, quadruple>::
C_DFullLDLt_queue(vector<C_task*>& queue,
		  vector<int>& task_indcol,
		  vector<int>& task_ptr,
		  double *eps_piv,
		  bool *kernel_detection,
		  int *aug_dim,
		  quadruple *eps_machine,
		  double *pivot,
		  double *pivot0,
		  double *pivot1,
		  vector<C_task*>& task_p,
	          const bool isChldrnAlgnd,
		  const bool verbose,
		  FILE **fp);

template
int DissectionMatrix<complex<quadruple>, quadruple >::
C_DFullLDLt_queue(vector<C_task*>& queue,
		  vector<int>& task_indcol,
		  vector<int>& task_ptr,
		  double *eps_piv,
		  bool *kernel_detection,
		  int *aug_dim,
		  quadruple *eps_machine,
		  double *pivot,
		  double *pivot0,
		  double *pivot1,
		  vector<C_task*>& task_p,
                  const bool isChldrnAlgnd,
		  const bool verbose,
		  FILE **fp);
//

template<typename T, typename U>
int DissectionMatrix<T, U>::
C_DTRSMScale_queue(vector<C_task*> &queue,
		   vector<list<int> > &qparents_index,
		   vector<int> &queue_index,
		   vector<int> &indcol,
		   vector<int> &ptr,
		   Dissection::Tree *btree,
		   vector<C_task*> &task_o,
		   vector<int> &task_o_indcol,
		   vector<int> &task_o_ptr,
		   vector<C_task*> &task_p,
		   const bool verbose,
		   FILE **fp)
{
  int itmp, jtmp;
  const int nrow = getUpperNRow();  // _nrow
  const int num_block = _diag->num_blocks();
  const int task_p_offset = (num_block * (num_block + 1)) / 2;
  const int num_block_col = _upper->num_blocks_c(); //
  vector<int> &singidx = singIdx();
  if (nrow == 0) {
    queue.resize(1);
    string task_name = ("du dummy : " +
			to_string(_level) + " : " +to_string(_nb));
    C_dummy_arg *arg = new C_dummy_arg(verbose, fp, _nb);
    //    *(arg->ops_complexity) = (-1L);
    queue[0] = new C_task(C_DUMMY,
			  task_name,
			  (void *)arg,
			  C_dummy,
			  1, // atomic_size,
			  0, // atomic_id,
			  arg->ops_complexity);
    queue[0]->parents->clear();
    return 1;
  }
  itmp = ((num_block * (num_block + 1)) / 2) * num_block_col;
  if (!_isSym) {
    itmp *= 2;
  }
  queue.resize(itmp);
  qparents_index.resize(itmp);
  queue_index.resize(itmp);
  // initialize
  for (int i = 0; i < itmp; i++) {
    queue_index[i] = i;
  }
  //  vector<int> ptr, indcol;
  ptr.resize(num_block + 1);
  jtmp = (num_block + 1) * num_block;
  itmp = _isSym ? (jtmp / 2) : jtmp;
  indcol.resize(itmp);             // upper and lower are stored sequentially
  ptr[0] = 0;
  for (int k = 0; k < num_block; k++) {
    ptr[k + 1] = ptr[k] + (num_block - k);  // based on the symmetric case
  }
  int ipos, jpos, atomic_size, atomic_id;
  jpos = 0;
  for (int k = 0; k < num_block; k++) {
    if (k == 0) {
      ipos = 0;
      atomic_size = 1;
      atomic_id = 0;
    }
    else {
      ipos = _isSym ? (ptr[k - 1] + 2) : (ptr[k - 1] + 2) * 2;
      atomic_size = 2;
      atomic_id = 1;
    }
    int nrow_block = _diag->nrowBlock(k, k);
    const long nrow_blockl = (long)nrow_block;
    int iipos = ipos * num_block_col;
    for (int l = 0; l < num_block_col; l++) {
      int nrhs_block = _upper->ncolBlock(l);
      const long nrhs_blockl = (long)nrhs_block;
      const long ops = nrow_blockl * nrow_blockl *nrhs_blockl;
      DTRSMScale_arg<T> *arg = 
	new DTRSMScale_arg<T>(_isSym,
			      _diag,  // ldlt,
			      _upper, //
			      _lower,  // not yet allocated
			      //     pt_offset,
			      nrow,
			      nrhs_block,
			      k,    //
			      l,
			      (-1), // mblock
			      &singidx,
			      true, // localPermute
			      verbose,
			      fp,
			      _nb);
      *(arg->ops_complexity) = ops;
      string  task_name = ("du " + to_string(k) + "   " + to_string(l) 
			   + " : " + to_string(_level) + " : "
			   + to_string(_nb));
      queue[iipos] = new C_task(C_DTRSMSCALE,
				task_name,
				(void *)arg,
				C_DTRSMScale_diag_upper<T>,
				atomic_size,
				atomic_id,
				arg->ops_complexity); 
      queue[iipos]->fp = &(arg->fp);
      // dependency on LDLt : alpha^(k)
      jtmp = task_o_indcol[task_o_ptr[k]];
      queue[iipos]->parents->push_back(task_o[jtmp]);
      if (k == 0) {
	queue[iipos]->parents->push_back(task_p[task_p_offset + l]);
      }
      else {
	if (_isSym) {
	  jtmp = indcol[ptr[k - 1]] * num_block_col + l;
	}
	else {
	  jtmp = indcol[2 * ptr[k - 1]] * num_block_col + l;
	}
	//	queue[iipos]->parents->push_back(queue[jtmp]);
	qparents_index[iipos].push_back(jtmp);
      }
      indcol[jpos] = ipos;
      iipos++;
    } //loop : l
    if (!_isSym) {
      for (int l = 0; l < num_block_col; l++) {
	int nrhs_block = _lower->ncolBlock(l);
	const long nrhs_blockl = (long)nrhs_block;
	// lower DTRSM without scaling
	// to prevent elimination from the task queue to keep task-dependency
	const long ops = nrow_blockl == 1L ? 1L : (nrow_blockl * (nrow_blockl - 1L) * nrhs_blockl); 
	DTRSMScale_arg<T> *arg = 
	  new DTRSMScale_arg<T>(_isSym,
				_diag, //ldlt,
				(RectBlockMatrix<T> *)NULL,
				_lower,   // not yet allocated
				//	pt_offset,
				nrow,
				nrhs_block,
				k,    //
				l,
				(-1), // mblock
				&singidx,
				true, // localPermute
				verbose,
				fp,
				_nb);
	*(arg->ops_complexity) = ops;

	string task_name = ("dl " + to_string(k) + "   " + to_string(l) 
			    + " : " + to_string(_level) + " : "
			    + to_string(_nb));
	queue[iipos] = new C_task(C_DTRSMSCALE,
				  task_name,
				  (void *)arg,
				  C_DTRSMScale_diag_lower<T>,
				  atomic_size,
				  atomic_id,
				  arg->ops_complexity);  
	queue[iipos]->fp = &(arg->fp);
	// dependency on LDLt : alpha^(k)
	jtmp = task_o_indcol[task_o_ptr[k]];
	queue[iipos]->parents->push_back(task_o[jtmp]);
	if (k == 0) {
	  queue[iipos]->parents->push_back(task_p[task_p_offset + l]);
	}
	else {
	  jtmp = indcol[2 * ptr[k - 1] + 1] * num_block_col + l;
	  //	  queue[iipos]->parents->push_back(queue[jtmp]);
	  qparents_index[iipos].push_back(jtmp);
	}
	iipos++;
	indcol[jpos + 1] = ipos + 1;
      }
    } // loop : l
    jpos += _isSym ? 1 : 2;
    for (int m = (k + 1); m < num_block; m++) {
      if (m == (k + 1)) {
	ipos = _isSym ? (ptr[k] + 1) : 
                       ((ptr[k] + 1) * 2);
	atomic_size = 2;
	atomic_id = 0;
      }
      else {
	ipos = _isSym ? (ptr[k] + (m - k + 1)) : 
                       ((ptr[k] + (m - k + 1)) * 2);
	atomic_size = 1;
	atomic_id = 0;
      }
      int ncol_block = _diag->ncolBlock(k, m); 
      const long ncol_blockl = (long)ncol_block;
      int iipos = ipos * num_block_col;
      for (int l = 0; l < num_block_col; l++) {
	int nrhs_block = _upper->ncolBlock(l);
	const long nrhs_blockl = (long)nrhs_block;
	const long ops = (2L * nrow_blockl * ncol_blockl * nrhs_blockl);
	DTRSMScale_arg<T> *arg = 
	  new DTRSMScale_arg<T>(_isSym,
				_diag, //ldlt,
				_upper,  //
				_lower,   // not yet allocated
				//	pt_offset, 
				nrow,
				nrhs_block,
				k,    //
				l,
				m, // mblock
				&singidx,
				true, // localPermute
				verbose,
				fp,
				_nb);

	string task_name = ("dU " + to_string(k) + " " + to_string(m)
	                   + " " + to_string(l) + " : " + to_string(_level)
	                   + " : " + to_string(_nb));

	*(arg->ops_complexity) = ops;
	queue[iipos] = new C_task(C_DTRSMSCALE,
				  task_name,
				  (void *)arg,
				  C_DTRSMScale_offdiag_upper<T>,
				  atomic_size,
				  atomic_id,
				  arg->ops_complexity);  
	queue[iipos]->fp = &(arg->fp);
	// Dependency on LDLt : beta^(k)_(m)
	jtmp = task_o_indcol[task_o_ptr[k] + m - k];
	queue[iipos]->parents->push_back(task_o[jtmp]);
	jtmp = (indcol[_isSym ? ptr[k] : ptr[k] * 2] * num_block_col + l);
	qparents_index[iipos].push_back(jtmp); //
	if (k > 0) {
  	  const int ktmp = (_isSym ? (ptr[k - 1] + m - k + 1) : 
		  	            ((ptr[k - 1] + m - k + 1) * 2));
	  jtmp = (indcol[ktmp] * num_block_col + l);
	  qparents_index[iipos].push_back(jtmp); //
	}
	indcol[jpos] = ipos;
	iipos++;
      }
      if (!_isSym) {
	for (int l = 0; l < num_block_col; l++) {
	  int nrhs_block = _lower->ncolBlock(l);
	  const long nrhs_blockl = (long)nrhs_block;
	  const long ops = (2L * nrow_blockl * ncol_blockl * nrhs_blockl);
	  DTRSMScale_arg<T> *arg = 
	    new DTRSMScale_arg<T>(_isSym,
				  _diag, //ldlt,
				  (RectBlockMatrix<T> *)NULL,
				  _lower, // not yet allocated
				  //			  pt_offset,
				  nrow,
				  nrhs_block,
				  k,    //
				  l,
				  m, // mblock
				  &singidx,
				  true, // localPermute
				  verbose,
				  fp,
				  _nb);
	  *(arg->ops_complexity) = ops;

	  string task_name = ("dL " + to_string(k) + " " + to_string(m)
			      + " " + to_string(l) + " : " + to_string(_level)
			      + " : " + to_string(_nb));

	  queue[iipos] = new C_task(C_DTRSMSCALE,
				    task_name,
				    (void *)arg,
				    C_DTRSMScale_offdiag_lower<T>,
				    atomic_size,
				    atomic_id,
				    arg->ops_complexity);  
	  queue[iipos]->fp = &(arg->fp);
	  // dependency on LDLt : beta^(k)_(m)
	  jtmp = task_o_indcol[task_o_ptr[k] + m - k];
	  queue[iipos]->parents->push_back(task_o[jtmp]);
	  jtmp = indcol[2 * ptr[k] + 1] * num_block_col + l;
	  qparents_index[iipos].push_back(jtmp); // k-th beginning block
	  if (k > 0) {
	    const int ktmp = (ptr[k - 1] + m - k + 1) * 2 + 1;
	    jtmp = (indcol[ktmp] * num_block_col + l);
	    qparents_index[iipos].push_back(jtmp); //
	  }
	  iipos++;
	  indcol[jpos + 1] = ipos + 1;
	} // if (!isSym) 
      } // loop : l
      jpos += _isSym ? 1 : 2;
    } // loop : m
  } // loop k
#if 0
  fprintf(*fp, "%s %d : C_DTRSMScale_queue()\n", __FILE__, __LINE__);
  fprintf(*fp, "ptr[] = ");
  for (vector<int>::const_iterator it = ptr.begin(); it != ptr.end(); ++it) {
    fprintf(*fp, "%d ", (*it));
  }
  fprintf(*fp, " ");
  fprintf(*fp, "indcol[] = ");
  for (vector<int>::const_iterator it = indcol.begin(); it != indcol.end(); ++it) {
    fprintf(*fp, "%d ", (*it));
  }
  fprintf(*fp, "\nnum_block_col = %d ", num_block_col);
  fprintf(*fp, ": num_block0 = %d num_block1 = %d\n", 
	  _localSchur->num_blocks0(), _localSchur->num_blocks1());
  {
    int j = 0;
    for (vector<C_task *>::const_iterator it = queue.begin(); 
	 it != queue.end(); ++it, j++) {
      fprintf(*fp, "%s : %d / %d : %d / %d : %d parents : ", 
	      (*it)->task_name, 
	      (*it)->atomic_id, (*it)->atomic_size,
	      (*it)->parallel_id, (*it)->parallel_max,
	    (*it)->parents->size());
      for (list<int>::const_iterator jt = qparents_index[j].begin();
	   jt != qparents_index[j].end(); ++jt) {
      fprintf(*fp, "%s : ", queue[(*jt)]->task_name);
      }
      fprintf(*fp, "\n");
    }
  }
#endif
  return queue.size();
}

template
int DissectionMatrix<double, double>::
C_DTRSMScale_queue(vector<C_task*> &queue,
		   vector<list<int> > &qparents_index,
		   vector<int> &queue_index,
		   vector<int> &indcol,
		   vector<int> &ptr,
		   Dissection::Tree *btree,
		   vector<C_task*> &task_o,
		   vector<int> &task_o_indcol,
		   vector<int> &task_o_ptr,
		   vector<C_task*> &task_p,
		   const bool verbose,
		   FILE **fp);

template
int DissectionMatrix<complex<double>, double>::
C_DTRSMScale_queue(vector<C_task*> &queue,
		   vector<list<int> > &qparents_index,
		   vector<int> &queue_index,
		   vector<int> &indcol,
		   vector<int> &ptr,
		   Dissection::Tree *btree,
		   vector<C_task*> &task_o,
		   vector<int> &task_o_indcol,
		   vector<int> &task_o_ptr,
		   vector<C_task*> &task_p,
		   const bool verbose,
		   FILE **fp);

template
int DissectionMatrix<quadruple, quadruple>::
C_DTRSMScale_queue(vector<C_task*> &queue,
		   vector<list<int> > &qparents_index,
		   vector<int> &queue_index,
		   vector<int> &indcol,
		   vector<int> &ptr,
		   Dissection::Tree *btree,
		   vector<C_task*> &task_o,
		   vector<int> &task_o_indcol,
		   vector<int> &task_o_ptr,
		   vector<C_task*> &task_p,
		   const bool verbose,
		   FILE **fp);

template
int DissectionMatrix<complex<quadruple>, quadruple>::
C_DTRSMScale_queue(vector<C_task*> &queue,
		   vector<list<int> > &qparents_index,
		   vector<int> &queue_index,
		   vector<int> &indcol,
		   vector<int> &ptr,
		   Dissection::Tree *btree,
		   vector<C_task*> &task_o,
		   vector<int> &task_o_indcol,
		   vector<int> &task_o_ptr,
		   vector<C_task*> &task_p,
		   const bool verbose,
		   FILE **fp);
//


template<typename T, typename U>
void DissectionMatrix<T, U>::
C_DTRSMScale_rearrange(vector<C_task*> &queue,
		       vector<list<int> > &qparents_index,
		       vector<int> &queue_index,
		       vector<int>& indcol,
		       vector<int>& ptr,
		       const bool verbose,
		       FILE *fp)
{
  int itmp, jtmp;
  const int num_block = _diag->num_blocks();
  const int num_block_col = _upper->num_blocks_c();
 
  vector<C_task *> tmp;
  vector<int> tmp_index;

  if (num_block == 1) {
    if (!_isSym) {
      tmp.resize(2 * num_block_col);
      tmp_index.resize(2 * num_block_col);
      itmp = 0;
      for (int l = 0; l < _localSchur->num_blocks0(); l++, itmp += 2) {
	tmp[itmp] = queue[l];                     // upper
	tmp_index[itmp] = queue_index[l];
	tmp[itmp + 1] = queue[l + num_block_col]; // lower
	tmp_index[itmp + 1] = queue_index[l + num_block_col];
      }
      for (int l = _localSchur->num_blocks0(); l < num_block_col; l++, 
	                                                   itmp += 2) {
	tmp[itmp] = queue[l];                     // upper
	tmp_index[itmp] = queue_index[l];  
	tmp[itmp + 1] = queue[l + num_block_col]; // lower
	tmp_index[itmp + 1] = queue_index[l + num_block_col];
      }
      for (int i = 0; i < (2 * num_block_col); i++) {
	queue[i] = tmp[i];
	queue_index[i] = tmp_index[i];
      }
    }
  }
  else if (num_block == 2) {
    tmp.resize(3 * num_block_col);
    tmp_index.resize(3 * num_block_col);
    if (!_isSym) {
      swap_queues_n(queue, queue_index, 1, 2, num_block_col, tmp, tmp_index);
      swap_queues_n(queue, queue_index, 3, 4, num_block_col, tmp, tmp_index);
      swap_queues_n(queue, queue_index, 2, 3, num_block_col, tmp, tmp_index);
    }
    for (int l = 0; l < num_block_col; l++) {
      tmp[3 * l] = queue[l];
      tmp_index[3 * l] = queue_index[l];
      tmp[3 * l + 1] = queue[num_block_col + l];
      tmp_index[3 * l + 1] = queue_index[num_block_col + l];
      tmp[3 * l + 2] = queue[2 * num_block_col + l];
      tmp_index[3 * l + 2] = queue_index[2 * num_block_col + l];
      for (int m = 0; m < 3; m++) {
	tmp[3 * l + m]->atomic_size = 3;
	tmp[3 * l + m]->atomic_id = m;
      }
    }
    for (int l = 0; l < (3 * num_block_col); l++) {
      queue[l] = tmp[l];
      queue_index[l] = tmp_index[l];
    } // loop : l
    if (!_isSym) {
      for (int l = 0; l < num_block_col; l++) {
	tmp[3 * l] = queue[3 * num_block_col + l];
	tmp_index[3 * l] = queue_index[3 * num_block_col + l];
	tmp[3 * l + 1] = queue[4 * num_block_col + l];
	tmp_index[3 * l + 1] = queue_index[4 * num_block_col + l];
	tmp[3 * l + 2] = queue[5 * num_block_col + l];
	tmp_index[3 * l + 2] = queue_index[5 * num_block_col + l];
	for (int m = 0; m < 3; m++) {
	  tmp[3 * l + m]->atomic_size = 3;
	  tmp[3 * l + m]->atomic_id = m;
	}
      }
      for (int l = 0; l < (3 * num_block_col); l++) {
	queue[3 * num_block_col + l] = tmp[l];
	queue_index[3 * num_block_col + l] = tmp_index[l];
      } // loop : l
    } 
    if (!_isSym) {
      tmp.resize(6 * num_block_col);
      tmp_index.resize(6 * num_block_col);
      itmp = 0;
      for (int l = 0; l < _localSchur->num_blocks0(); l++, itmp += 6) {
	for (int k = 0; k < 3; k++) {
	  tmp[itmp + k] = queue[3 * l + k];                       // upper
	  tmp_index[itmp + k] = queue_index[3 * l + k];
	  tmp[itmp + 3 + k] = queue[3 * (l + num_block_col) + k]; // lower
	  tmp_index[itmp + 3 + k] = queue_index[3 * (l + num_block_col) + k];
	}
      }
      for (int l = _localSchur->num_blocks0(); l < num_block_col; l++, 
	                                                   itmp += 6) {
	for (int k = 0; k < 3; k++) {
	  tmp[itmp + k] = queue[3 * l + k];                       // upper
	  tmp_index[itmp + k] = queue_index[3 * l + k];
	  tmp[itmp + 3 + k] = queue[3 * (l + num_block_col) + k]; // lower
	  tmp_index[itmp + 3 + k] = queue_index[3 * (l + num_block_col) + k];
	}
      }
      for (int i = 0; i < (6 * num_block_col); i++) {
	queue[i] = tmp[i];
	queue_index[i] = tmp_index[i];
      }
    }

  }   // if (num_block == 1) else if (num_block == 2) 
  else {
    tmp.resize(2 * num_block_col);
    tmp_index.resize(2 * num_block_col);
    if (!_isSym) {
      for (int i = 0; i < indcol.size(); i++) {
	const int ii = i * num_block_col;
	if ((i % 2 == 1) && (queue[ii]->atomic_size == 2) &&
	    (queue[ii]->atomic_id == 0)) {
	  swap_queues_n(queue, queue_index, i, (i + 1), num_block_col, tmp, 
			tmp_index);
	  i++;
	}
      } // loop : i
    }   //  if (!isSym) 
    for (int i = 0; i < indcol.size(); i++) {
      const int ii = i * num_block_col;
      if ((queue[ii]->atomic_size == 2) && (queue[ii]->atomic_id == 0)) {
	for (int l = 0; l < num_block_col; l++) {
	  tmp[2 * l] = queue[ii + l];
	  tmp_index[2 * l] = queue_index[ii + l];
	  tmp[2 * l + 1] = queue[(i + 1) * num_block_col + l];
	  tmp_index[2 * l + 1] = queue_index[(i + 1) * num_block_col + l];
	}
	for (int l = 0; l < (2 * num_block_col); l++) {
	  queue[ii + l] = tmp[l];
	  queue_index[ii + l] = tmp_index[l];
	}
	i++;
      }
    } // loop : it
  }  // if (num_block == 2) 
  
  if (num_block < 3) {
    if (_isSym) {
      int parallel_max = _localSchur->num_blocks0();
      int nsize = ((num_block == 1) ? 1 : 3);
      for (int l = 0; l < parallel_max; l++) {
	itmp = l * nsize;
	for (int k = 0; k < nsize; k++, itmp++) {
	  queue[itmp]->parallel_max = parallel_max;
	  queue[itmp]->parallel_id = l;
	}
      }
      parallel_max = num_block_col - _localSchur->num_blocks0();
      for (int l = 0; l < parallel_max; l++) { 
	itmp = (l + _localSchur->num_blocks0()) * nsize;
	for (int k = 0; k < nsize; k++, itmp++) {
	  queue[itmp]->parallel_max = parallel_max;
	  queue[itmp]->parallel_id = l;
	}
      }
    }
    else {
      int parallel_max = _localSchur->num_blocks0() * 2;
      int nsize = ((num_block == 1) ? 1 : 3);
      for (int l = 0; l < parallel_max; l++) {
	itmp = l * nsize;
	for (int k = 0; k < nsize; k++, itmp++) {
	  queue[itmp]->parallel_max = parallel_max;
	  queue[itmp]->parallel_id = l;
	}
      }
      parallel_max = (num_block_col - _localSchur->num_blocks0()) * 2;
      for (int l = 0; l < parallel_max; l++) { 
	itmp = (l + 2 * _localSchur->num_blocks0()) * nsize;
	for (int k = 0; k < nsize; k++, itmp++) {
	  queue[itmp]->parallel_max = parallel_max;
	  queue[itmp]->parallel_id = l;
	}
      }
    }
  }
  else { // if (num_block < 3) 
    int parallel_max = (_isSym ? 1 : 2) * num_block_col;
    for (int l = 0; l < num_block_col; l++) {
      queue[l]->parallel_max = parallel_max;
      queue[l]->parallel_id = l;
    }
    if (!_isSym) {
      for (int l = 0; l < num_block_col; l++) {
	queue[num_block_col + l]->parallel_max = parallel_max;
	queue[num_block_col + l]->parallel_id = num_block_col + l;
      }
    }
    for (int k = 0; k < (num_block - 1); k++) {
      if (_isSym) {
	parallel_max = (ptr[k + 1] - ptr[k] - 1) * num_block_col;
	jtmp = 0;
	for (int i = (ptr[k] + 1); i < (ptr[k + 1] + 1); i++) {
	  itmp = i * num_block_col;
	  for (int l = 0; l < num_block_col; l++, itmp++) {
	    queue[itmp]->parallel_max = parallel_max;
	    queue[itmp]->parallel_id = jtmp;
	    if (((queue[itmp]->atomic_size == 2) && 
		 (queue[itmp]->atomic_id == 1)) || 
		(queue[itmp]->atomic_size == 1)) {
	      jtmp++;
	    }
	  } // loop : l
	} // loop : i
      }   // _isSym
      else {
	parallel_max = (2 * (ptr[k + 1] - ptr[k]) - 2) * num_block_col;
	jtmp = 0;
	for (int i = (ptr[k] * 2 + 2); i < (ptr[k + 1] * 2 + 2); i++) {
	  itmp = i * num_block_col;
	  for (int l = 0; l < num_block_col; l++, itmp++) {
	    queue[itmp]->parallel_max = parallel_max;
	    queue[itmp]->parallel_id = jtmp;
	    if (((queue[itmp]->atomic_size == 2) && 
		 (queue[itmp]->atomic_id == 1)) || 
		(queue[itmp]->atomic_size == 1)) {
	      jtmp++;
	    }
	  } // loop : l
	} // loop : i
      }
    } // loop : k
  } // if (num_block < 2)

  tmp_index.resize(queue_index.size());
  vector<int> tmp_old(queue_index);
  for (int i = 0; i < queue_index.size(); i++) {
    tmp_index[queue_index[i]] = i;
  }
  for (int i = 0; i < queue_index.size(); i++) {
    queue_index[i] = tmp_index[i];
  }
  { // begin : scope j
    int j = 0;
    for (vector<C_task *>::const_iterator it = queue.begin(); 
	 it != queue.end(); ++it, j++) {
      const int jj = tmp_old[j];
      for (list<int>::const_iterator kt = qparents_index[jj].begin(); 
	   kt != qparents_index[jj].end(); ++kt) {
	(*it)->parents->push_back(queue[queue_index[(*kt)]]);
      } // loop : kt
    } // loop : it
 }  // end : scope j

  // merging all dependency inside of atomic operation : not optimal but 
  // minimal modification of cheking of task dependency during excution
  for (vector<C_task *>::const_iterator it = queue.begin(); 
                                  it != queue.end(); ++it) {
    switch((*it)->atomic_size) {
    case 2:
      {
	vector<C_task *>::const_iterator jt = it;
	C_task *task_tmp = (*jt); 
	++it;
	(*jt)->parents->merge(*((*it)->parents));
	(*jt)->parents->unique();
	// remove dependency inside of atomic operation
	for (list<C_task*>::iterator mt = (*jt)->parents->begin();
	     mt != (*jt)->parents->end(); ++mt) {
	  if ((*mt) == (*jt)) {
	    mt = (*jt)->parents->erase(mt);
	  }
	}
	*((*it)->parents) = *((*jt)->parents);
	(*it)->parents->push_back(task_tmp);
      }
      break;
    case 3:
      {
	vector<C_task *>::const_iterator jt = it;
	C_task *task_tmp0 = (*jt);
	++it;
	(*jt)->parents->merge(*((*it)->parents));
	(*jt)->parents->unique();
	// remove dependency inside of atomic operation
	for (list<C_task*>::iterator mt = (*jt)->parents->begin();
	     mt != (*jt)->parents->end(); ++mt) {
	  if ((*mt) == (*jt)) {
	    task_tmp0 = (*mt);
	    mt = (*jt)->parents->erase(mt);
	  }
	}
	vector<C_task *>::const_iterator kt = it;
	C_task *task_tmp1 = (*kt);
	++it;
	(*jt)->parents->merge(*((*it)->parents));
	(*jt)->parents->unique();
	// remove dependency inside of atomic operation
	for (list<C_task*>::iterator mt = (*jt)->parents->begin();
	     mt != (*jt)->parents->end(); ++mt) {
	  if (((*mt) == (*jt)) || ((*mt) == (*kt))) {
	    mt = (*jt)->parents->erase(mt);
	  }
	}
	*((*kt)->parents) = *((*jt)->parents);
	*((*it)->parents) = *((*jt)->parents);
	(*kt)->parents->push_back(task_tmp0);
	(*it)->parents->push_back(task_tmp0);
	(*it)->parents->push_back(task_tmp1);
      }
      break;
    default:
      break;
    }
  }

  for (vector<C_task *>::const_iterator it = queue.begin(); 
                                  it != queue.end(); ++it) {
    (*it)->parents->sort(compare_task_name);
    (*it)->parents->unique();
  }
#if 0
  fprintf(fp, "%s %d : C_DTRSMScale_rearrange()\n", __FILE__, __LINE__);
  for (vector<C_task *>::const_iterator it = queue.begin(); 
                                  it != queue.end(); ++it) {
    fprintf(fp, "%s : %d / %d : %d / %d : %d parents : ", 
	    (*it)->task_name, 
	    (*it)->atomic_id, (*it)->atomic_size,
	    (*it)->parallel_id, (*it)->parallel_max,
	    (*it)->parents->size());
    for (list<C_task*>::const_iterator jt = (*it)->parents->begin(); 
	 jt != (*it)->parents->end(); ++jt) {
      fprintf(fp, "%s : ", (*jt)->task_name);
    }
    fprintf(fp, "\n");
  }

  fprintf(fp, "%s %d C_DTRSMScale_rearrange() end : ", 
	  __FILE__, __LINE__);
  for (vector<int>::const_iterator it = queue_index.begin(); it != queue_index.end();
       ++it) {
    fprintf(fp, "%d ", *it);
  }
  fprintf(fp, "\n");
#endif
}

template
void DissectionMatrix<double, double>::
C_DTRSMScale_rearrange(vector<C_task*> &queue,
		       vector<list<int> > &qparents_index,
		       vector<int> &queue_index,
		       vector<int>& indcol,
		       vector<int>& ptr,
		       const bool verbose,
		       FILE *fp);

template
void DissectionMatrix<complex<double>, double>::
C_DTRSMScale_rearrange(vector<C_task*> &queue,
		       vector<list<int> > &qparents_index,
		       vector<int> &queue_index,
		       vector<int>& indcol,
		       vector<int>& ptr,
		       const bool verbose,
		       FILE *fp);

template
void DissectionMatrix<quadruple, quadruple>::
C_DTRSMScale_rearrange(vector<C_task*> &queue,
		       vector<list<int> > &qparents_index,
		       vector<int> &queue_index,
		       vector<int>& indcol,
		       vector<int>& ptr,
		       const bool verbose,
		       FILE *fp);

template
void DissectionMatrix<complex<quadruple>, quadruple>::
C_DTRSMScale_rearrange(vector<C_task*> &queue,
		       vector<list<int> > &qparents_index,
		       vector<int> &queue_index,
		       vector<int>& indcol,
		       vector<int>& ptr,
		       const bool verbose,
		       FILE *fp);
//

template<typename T, typename U>
int DissectionMatrix<T, U>::C_DGEMM_local_queue(vector<C_task*> &queue,
						 vector<int> &indcol,
						 RectBlockMatrix<T> *upper1,
						 RectBlockMatrix<T> *lower1,
						 bool isSkip,
						 bool isDirect,
						 SquareBlockMatrix<T>* fdiag,
						 vector<C_task *> &task_p,
						 vector<int> &task_p_index,
						 vector<int> &task_p_indcol,
						 vector<int> &task_p_ptr,
						 vector<C_task *> &task_q,
						 vector<int> &task_q_index,
						 vector<int> &task_q_indcol,
						 vector<int> &task_q_ptr,
						 vector<C_task *> *task_s,
						 const bool verbose,
						 FILE **fp)
{
  const int nrow = _upper->nbRows();     // _nrow
  const int nrow1 = isSkip ? (-1) : upper1->nbRows(); // dummy for isSkip

  const int num_block = _localSchur->num_blocks(); 
  const int num_block0 = _localSchur->num_blocks0(); 
  const int num_block1 = num_block - num_block0;
  int queue_size;
  if (nrow == 0 && _level != 0) {
    queue.resize(1);
    string task_name = ("f dummy : " +
			to_string(_level) + " : " +to_string(_nb));
    fprintf(*fp, "%s %d : %s\n", __FILE__, __LINE__, task_name.c_str());
    C_dummy_arg *arg = new C_dummy_arg(verbose, fp, _nb);
    // *(arg->ops_complexity) = (-1L);
    queue[0] = new C_task(C_DUMMY,
			  task_name,
			  (void *)arg,
			  C_dummy,
			  1,  // atomic_size
			  0,  // atomic_id
			  arg->ops_complexity);
    queue[0]->parents->clear();
    return 1;
  }
  queue_size = (_isSym ? ((num_block * (num_block + 1)) / 2) : 
              		  (num_block * num_block));

  queue.resize(queue_size);
  indcol.resize(queue_size);
  int pncol_block = _upper->num_blocks_c(); //
  int qncol_block = isSkip ? (-1) : upper1->num_blocks_c(); //
 
  int ipos = 0;
  int parallel_max0, parallel_max1, parallel_max2;
  parallel_max0 = (_isSym ? (num_block0 * (num_block0 + 1)) / 2 : 
		            (num_block0 * num_block0));
  parallel_max2 = (_isSym ? (num_block1 * (num_block1 + 1)) / 2 : 
		            (num_block1 * num_block1));
  parallel_max1 = (_isSym ? (num_block * (num_block + 1)) / 2 : 
		  (num_block * num_block)) - parallel_max0 - parallel_max2;
  for (int j = 0; j < num_block0; j++) {
    for (int i = 0; i < j; i++) {
      const long block_nrowl = (long)fdiag->nrowBlock(i, j);
      const long block_ncoll = (long)fdiag->ncolBlock(i, j);
      const long nrowl = (long)nrow;
      const long nrowwl = isSkip ? nrowl : nrowl + (long)nrow1;
      {
	if (isDirect) {
	  const long ops = (isSkip ? 0L :   // flag to skip the dependency
			    (block_nrowl * block_ncoll * nrowwl * 2L));
	  string task_name = ("F " + to_string(i) + " " + to_string(j)
			      + " : " + to_string(_level) + " : "
			      + to_string(_nb));
	  DSchurGEMM_two_arg<T> *arg = 
	    new DSchurGEMM_two_arg<T>(_isSym,
				      false, // isTrans
				      _lower,
				      _upper, 
				      nrow,
				      lower1,
				      upper1, 
				      nrow1,
				      i,
				      j, 
				      fdiag,
				      isSkip,
				      verbose,
				      fp,
				      _nb);
	  *(arg->ops_complexity) = ops;
	  queue[ipos] = new C_task(C_DGEMM_DIRECT_TWO,
				   task_name,
				   (void *)arg,
				   DSchurGEMM_offdiag_two<T>,
				   1,              // atomic_size
				   0,              // atomic_id
				   arg->ops_complexity);
	  queue[ipos]->fp = &(arg->fp);
	  if(!isSkip) {
	    const int mmp = task_p_ptr.size() - 1;
	    const int mmq = task_q_ptr.size() - 1;
	    if (_isSym) {
	      for (int m = 0; m < task_p_ptr[mmp]; m++) {
		int nn;
		nn = task_p_indcol[m] * pncol_block + i;
		queue[ipos]->parents->push_back(task_p[task_p_index[nn]]);	
		nn = task_p_indcol[m] * pncol_block + j;
		queue[ipos]->parents->push_back(task_p[task_p_index[nn]]);
	      }
	      for (int m = 0; m < task_q_ptr[mmq]; m++) {
		int nn;
		nn = task_q_indcol[m] * qncol_block + i;
		queue[ipos]->parents->push_back(task_q[task_q_index[nn]]);	
		nn = task_q_indcol[m] * qncol_block + j;
		queue[ipos]->parents->push_back(task_q[task_q_index[nn]]);
	      }
	    }
	    else {
	      for (int m = 0; m < task_p_ptr[mmp]; m++) {
		int nn;
		nn = task_p_indcol[2 * m + 1] * pncol_block + i;   // lower
		queue[ipos]->parents->push_back(task_p[task_p_index[nn]]);	
		nn = task_p_indcol[2 * m] * pncol_block + j;       // upper
		queue[ipos]->parents->push_back(task_p[task_p_index[nn]]);
	      }
 	      for (int m = 0; m < task_q_ptr[mmq]; m++) {
		int nn;
		nn = task_q_indcol[2 * m + 1] * qncol_block + i;   // lower
		queue[ipos]->parents->push_back(task_q[task_q_index[nn]]);	
		nn = task_q_indcol[2 * m] * qncol_block + j;       // upper
		queue[ipos]->parents->push_back(task_q[task_q_index[nn]]);
	      }
	    }
	    if (task_s != NULL) {
	      const int ktmp = (j * (j + 1)) / 2 + i;  // diagonal block DSUB
	      queue[ipos]->parents->push_back((*task_s)[ktmp]);
	    }
	  } // if (!isSkip)
	} // if (isDirect)
	else {
	  const long ops = (isSkip ? 0L :   // flag to skip the dependency
			    (block_nrowl * block_ncoll * nrowl * 2L));
	  string task_name = ("f " + to_string(i) + " " + to_string(j)
  	                     + " : " + to_string(_level) + " : "
	                     + to_string(_nb));
	  DSchurGEMM_arg<T> *arg = 
	    new DSchurGEMM_arg<T>(_isSym,
				  false, // isTrans
				  _lower,
				  _upper, 
				  nrow,
				  i,
				  j, 
				  _localSchur,
				  verbose,
				  fp,
				  _nb);
	  *(arg->ops_complexity) = ops;
	  queue[ipos] = new C_task(C_DGEMM_LOCAL_TWO,
				   task_name,
				   (void *)arg,
				   DSchurGEMM_offdiag<T>,
				   1,              // atomic_size
				   0,              // atomic_id
				   arg->ops_complexity);
	  queue[ipos]->fp = &(arg->fp);
	  const int mmp = task_p_ptr.size() - 1;
	  if (_isSym) {
	    for (int m = 0; m < task_p_ptr[mmp]; m++) {
	      int nn;
	      nn = task_p_indcol[m] * pncol_block + i;       // lower
	      queue[ipos]->parents->push_back(task_p[task_p_index[nn]]);	
	      nn = task_p_indcol[m] * pncol_block + j;       // upper
	      queue[ipos]->parents->push_back(task_p[task_p_index[nn]]);
	    }
	  }
	  else {
	    for (int m = 0; m < task_p_ptr[mmp]; m++) {
	      int nn;
	      nn = task_p_indcol[2 * m + 1] * pncol_block + i; // lower
	      queue[ipos]->parents->push_back(task_p[task_p_index[nn]]);	
	      nn = task_p_indcol[2 * m] * pncol_block + j;     // upper
	      queue[ipos]->parents->push_back(task_p[task_p_index[nn]]);
	    }
	  }
	} // else if (isDirect)
	queue[ipos]->parents->unique();
	queue[ipos]->parallel_max = parallel_max0;
	queue[ipos]->parallel_id = ipos;
	const int itmp = _isSym ? ((j * (j + 1)) / 2 + i) : (j * j + 2 * i);
	indcol[itmp] = ipos;
	ipos++;      // lower block
      }
      if (!_isSym) {
	if (isDirect) {
	  const long ops = (isSkip ? 0L :   // flag to skip the dependency
			    (block_nrowl * block_ncoll * nrowwl * 2L));
	  string task_name = ("F " + to_string(j) + " " + to_string(i)
			      + " : " + to_string(_level) + " : "
			      + to_string(_nb));
	  DSchurGEMM_two_arg<T> *arg = 
	    new DSchurGEMM_two_arg<T>(_isSym,
				      true, // isTrans
				      _upper, 
				      _lower,
				      nrow,
				      upper1, 
				      lower1,
				      nrow1,
				      j,  // transposed
				      i, 
				      fdiag,
				      isSkip, //
				      verbose,
				      fp,
				      _nb);
	  *(arg->ops_complexity) = ops;
	  queue[ipos] = new C_task(C_DGEMM_DIRECT_TWO,
				   task_name,
				   (void *)arg,
				   DSchurGEMM_offdiag_two<T>,
				   1, 
				   0, 
				   arg->ops_complexity);
	  queue[ipos]->fp = &(arg->fp);
    	  if(!isSkip) {
	    for (int m = 0; m < (task_p_ptr.size() - 1); m++) {
	      int mm, nn;
	      mm = 2 * task_p_ptr[m]; // upper
	      nn = task_p_indcol[mm] * pncol_block + i;
	      queue[ipos]->parents->push_back(task_p[task_p_index[nn]]);
	      mm = 2 * task_p_ptr[m] + 1;  // lower
	      nn = task_p_indcol[mm] * pncol_block + j;
	      queue[ipos]->parents->push_back(task_p[task_p_index[nn]]);
	    }
	    for (int m = 0; m < (task_q_ptr.size() - 1); m++) {
	      int mm, nn;
	      mm = 2 * task_q_ptr[m]; // upper
	      nn = task_q_indcol[mm] * qncol_block + i;
	      queue[ipos]->parents->push_back(task_q[task_q_index[nn]]);
	      mm = 2 * task_q_ptr[m] + 1;  // lower
	      nn = task_q_indcol[mm] * qncol_block + j;
	      queue[ipos]->parents->push_back(task_q[task_q_index[nn]]);
	    }
	    if (task_s != NULL) {
	      const int ktmp = (j * (j + 1)) / 2 + i;   // diagonal block DSUB
	      queue[ipos]->parents->push_back((*task_s)[ktmp]);
	    }
	  }
	} // if (isDirect) 
	else {
	  const long ops = (isSkip ? 0L :   // flag to skip the dependency
			    (block_nrowl * block_ncoll * nrowl * 2L));
	  string task_name = ("f " + to_string(j) + " " + to_string(i)
	                     + " : " + to_string(_level) + " : "
	                     + to_string(_nb));
	  DSchurGEMM_arg<T> *arg = 
	    new DSchurGEMM_arg<T>(_isSym,
				  true, // isTrans
				  _upper, 
				  _lower,
				  nrow,
				  j,   // trasnposed
				  i, 
				  _localSchur,
				  verbose,
				  fp,
				  _nb);
	  *(arg->ops_complexity) = ops;
	  queue[ipos] = new C_task(C_DGEMM_LOCAL_TWO,
				   task_name,
				   (void *)arg,
				   DSchurGEMM_offdiag<T>,
				   1, 
				   0, 
				   arg->ops_complexity); 
	  queue[ipos]->fp = &(arg->fp);
	  const int mmp = task_p_ptr.size() - 1;
	  for (int m = 0; m < task_p_ptr[mmp]; m++) {
	    int nn;
	    nn = task_p_indcol[2 * m + 1] * pncol_block + j;   // lower
	    queue[ipos]->parents->push_back(task_p[task_p_index[nn]]);	
	    nn = task_p_indcol[2 * m] * pncol_block + i;       // upper
	    queue[ipos]->parents->push_back(task_p[task_p_index[nn]]);
	  }
	} // else if (isDirect)
	queue[ipos]->parents->unique();
	queue[ipos]->parallel_max = parallel_max0;
	queue[ipos]->parallel_id = ipos;
	const int itmp = j * j + 2 * i + 1; // lower part of unsymmetric matrix
	indcol[itmp] = ipos;
	ipos++;
      } // (!_isSym)
    } // loop : i
      // diagonal
    {
      const long block_ncoll = (long)fdiag->ncolBlock(j, j);
      const long nrowl = (long)nrow;
      const long nrowwl = isSkip ? nrowl : nrowl + (long)nrow1;
      if (isDirect) {
	const long ops = (isSkip ? 0L :
 			 (_isSym ? (block_ncoll * (block_ncoll + 1L) * nrowwl) :
			           (block_ncoll * block_ncoll * nrowwl * 2L)));
	string task_name = ("F " + to_string(j) + " " + to_string(j)
			    + " : " + to_string(_level) + " : "
			    + to_string(_nb));
	DSchurGEMM_two_arg<T> *arg = 
	  new DSchurGEMM_two_arg<T>(_isSym,
				    false, //
				    _lower,  // 
				    _upper, //
				    nrow,
				    lower1,  // 
				    upper1, //
				    nrow1,
				    j,
				    (-1),			       
				    fdiag,
				    isSkip,
				    verbose,
				    fp,
				    _nb);
	*(arg->ops_complexity) = ops;
	queue[ipos] = new C_task(C_DGEMM_DIRECT_TWO,
				 task_name, // task_name_cstr,
				 (void *)arg,
				 DSchurGEMM_diag_two<T>,
				 1,              // atomic_size
				 0,              // atomic_id
				 arg->ops_complexity);
	queue[ipos]->fp = &(arg->fp);
	if(!isSkip) {
	  const int mmp = task_p_ptr.size() - 1;
	  const int mmq = task_q_ptr.size() - 1;
	  if (_isSym) {
	    for (int m = 0; m < task_p_ptr[mmp]; m++) {
	      int nn;
	      nn = task_p_indcol[m] * pncol_block + j;
	      queue[ipos]->parents->push_back(task_p[task_p_index[nn]]);	
	    }
	    for (int m = 0; m < task_q_ptr[mmq]; m++) {
	      int nn;
	      nn = task_q_indcol[m] * qncol_block + j;
	      queue[ipos]->parents->push_back(task_q[task_q_index[nn]]);	
	    }
	  }
	  else {
	    for (int m = 0; m < task_p_ptr[mmp]; m++) {
	      int nn;
	      nn = task_p_indcol[2 * m + 1] * pncol_block + j;   // lower
	      queue[ipos]->parents->push_back(task_p[task_p_index[nn]]);	
	      nn = task_p_indcol[2 * m] * pncol_block + j;       // upper
	      queue[ipos]->parents->push_back(task_p[task_p_index[nn]]);
	    }
	    for (int m = 0; m < task_q_ptr[mmq]; m++) {
	      int nn;
	      nn = task_q_indcol[2 * m + 1] * qncol_block + j;   // lower
	      queue[ipos]->parents->push_back(task_q[task_q_index[nn]]);	
	      nn = task_q_indcol[2 * m] * qncol_block + j;       // upper
	      queue[ipos]->parents->push_back(task_q[task_q_index[nn]]);
	    }
	  }
	  if (task_s != NULL) {
	    const int ktmp = (j * (j + 1)) / 2 + j;   // diagonal block DSUB
	    queue[ipos]->parents->push_back((*task_s)[ktmp]);
	  }
	}
      }  // if (isDirect) 
      else {
	const long ops = (isSkip ? 0L :
			 (_isSym ? (block_ncoll * (block_ncoll + 1L) * nrowl) : 
			           (block_ncoll * block_ncoll * nrowl * 2L)));
	string task_name = ("f " + to_string(j) + " " + to_string(j)
			    + " : " + to_string(_level) + " : "
			    + to_string(_nb));
	DSchurGEMM_arg<T> *arg = 
	  new DSchurGEMM_arg<T>(_isSym,
				false, // isTrans
				_lower,  // 
				_upper, //
				nrow,
				j,
				(-1),			       
				_localSchur,
				verbose,
				fp,
				_nb);
	*(arg->ops_complexity) = ops;
	queue[ipos] = new C_task(C_DGEMM_LOCAL_TWO,
				 task_name, // task_name_cstr,
				 (void *)arg,
				 DSchurGEMM_diag<T>,
				 1,              // atomic_size
				 0,              // atomic_id
				 arg->ops_complexity);  
	queue[ipos]->fp = &(arg->fp);
	const int mmp = task_p_ptr.size() - 1;
	if (_isSym) {
	  for (int m = 0; m < task_p_ptr[mmp]; m++) {
	    int nn;
	    nn = task_p_indcol[m] * pncol_block + j;
	    queue[ipos]->parents->push_back(task_p[task_p_index[nn]]);	
	  }
	}
	else {
	  for (int m = 0; m < task_p_ptr[mmp]; m++) {
	    int nn;
	    nn = task_p_indcol[2 * m + 1] * pncol_block + j; // lower
	    queue[ipos]->parents->push_back(task_p[task_p_index[nn]]);	
	    nn = task_p_indcol[2 * m] * pncol_block + j;       // upper
	    queue[ipos]->parents->push_back(task_p[task_p_index[nn]]);
	  }
	}
      } // else if (isDirect)
      queue[ipos]->parents->unique();
      queue[ipos]->parallel_max = parallel_max0;
      queue[ipos]->parallel_id = ipos;
      const int itmp = _isSym ? ((j * (j + 1)) / 2 + j) : (j * j + 2 * j);
      indcol[itmp] = ipos;
      ipos++;
    } 
  } // loop : j
  for (int j = num_block0; j < num_block; j++) {
    for (int i = 0; i < num_block0; i++) {
      const long block_nrowl = (long)_localSchur->nrowBlock(i, j); 
      const long block_ncoll = (long)_localSchur->ncolBlock(i, j); 
      const long nrowl = (long)nrow;
      {
	const long ops = block_nrowl * block_ncoll * nrowl * 2L;
	string task_name = ("f " + to_string(i) + " " + to_string(j)
	                   + " : " + to_string(_level) + " : "
	                   + to_string(_nb));
	DSchurGEMM_arg<T> *arg = 
	  new DSchurGEMM_arg<T>(_isSym,
				false, // isTrans
				_lower,
				_upper, 
				nrow,
				i, //(i * SIZE_B1),
				j, //(j * SIZE_B1),
				_localSchur,
				verbose,
				fp,
				_nb);
	*(arg->ops_complexity) = ops;
	queue[ipos] = new C_task(C_DGEMM_LOCAL_TWO,
				 task_name,
				 (void *)arg,
				 DSchurGEMM_offdiag<T>,
				 1,              // atomic_size
				 0,              // atomic_id
				 arg->ops_complexity); 
	queue[ipos]->fp = &(arg->fp);
	const int mmp = task_p_ptr.size() - 1;
	if (_isSym) {
	  for (int m = 0; m < task_p_ptr[mmp]; m++) {
	    int nn;
	    nn = task_p_indcol[m] * pncol_block + i;
	    queue[ipos]->parents->push_back(task_p[task_p_index[nn]]);
	    nn = task_p_indcol[m] * pncol_block + j;
	    queue[ipos]->parents->push_back(task_p[task_p_index[nn]]);	
	  }
	}
	else {
	  for (int m = 0; m < task_p_ptr[mmp]; m++) {
	    int nn;
	    nn = task_p_indcol[2 * m + 1] * pncol_block + i;   // lower
	    queue[ipos]->parents->push_back(task_p[task_p_index[nn]]);	
	    nn = task_p_indcol[2 * m] * pncol_block + j;       // upper
	    queue[ipos]->parents->push_back(task_p[task_p_index[nn]]);
	  }
	}
	queue[ipos]->parents->sort(compare_task_name);
	queue[ipos]->parents->unique();
	queue[ipos]->parallel_max = parallel_max1;
	queue[ipos]->parallel_id = ipos - parallel_max0;
	const int itmp = _isSym ? ((j * (j + 1)) / 2 + i) : (j * j + 2 * i);
	indcol[itmp] = ipos;
	ipos++;      // lower block
      }
      if (!_isSym) {
	const long ops = block_nrowl * block_ncoll * nrowl * 2L;
	string task_name = ("f " + to_string(j) + " " + to_string(i)
		            + " : " + to_string(_level) + " : "
			    + to_string(_nb));
      	DSchurGEMM_arg<T> *arg = 
	  new DSchurGEMM_arg<T>(_isSym,
				true, // isTrans
				_upper, 
				_lower,
				nrow,
				j, // transposed
				i, 
				_localSchur,
				verbose,
				fp,
				_nb);
        *(arg->ops_complexity) = ops;
        queue[ipos] = new C_task(C_DGEMM_LOCAL_TWO,
				 task_name,
				 (void *)arg,
				 DSchurGEMM_offdiag<T>,
				 1,              // atomic_size
				 0,              // atomic_id
				 arg->ops_complexity);
	queue[ipos]->fp = &(arg->fp);
	const int mmp = task_p_ptr.size() - 1;
	for (int m = 0; m < task_p_ptr[mmp]; m++) {
	  int nn;
	  nn = task_p_indcol[2 * m + 1] * pncol_block + j;   // lower
	  queue[ipos]->parents->push_back(task_p[task_p_index[nn]]);	
	  nn = task_p_indcol[2 * m] * pncol_block + i;       // upper
	  queue[ipos]->parents->push_back(task_p[task_p_index[nn]]);
	}
	queue[ipos]->parents->sort(compare_task_name);
	queue[ipos]->parents->unique();
	queue[ipos]->parallel_max = parallel_max1;
	queue[ipos]->parallel_id = ipos - parallel_max0;	
	const int itmp = j * j + 2 * i + 1;
	indcol[itmp] = ipos;
	ipos++;
      } // (!_isSym)
    } // loop : i
  }
  for (int j = num_block0; j < num_block; j++) {
    for (int i = num_block0; i < j; i++) {
      const long block_nrowl = (long)_localSchur->nrowBlock(i, j); 
      const long block_ncoll = (long)_localSchur->ncolBlock(i, j); 
      const long nrowl = (long)nrow;
      {
	const long ops = block_nrowl * block_ncoll * nrowl * 2L;
	string task_name = ("f " + to_string(i) + " " + to_string(j)
			    + " : " + to_string(_level) + " : "
			    + to_string(_nb));
	DSchurGEMM_arg<T> *arg = 
	  new DSchurGEMM_arg<T>(_isSym,
				false, // isTrans
				_lower,
				_upper, 
				nrow,
				i, //(i * SIZE_B1),
				j, //(j * SIZE_B1),
				_localSchur,
				verbose,
				fp,
				_nb);
	*(arg->ops_complexity) = ops;
	queue[ipos] = new C_task(C_DGEMM_LOCAL_MULT,
				 task_name,
				 (void *)arg,
				 DSchurGEMM_offdiag<T>,
				 1,              // atomic_size
				 0,              // atomic_id
				 arg->ops_complexity);
	queue[ipos]->fp = &(arg->fp); 
	const int mmp = task_p_ptr.size() - 1;
	if (_isSym) {
	  for (int m = 0; m < task_p_ptr[mmp]; m++) {
	    int nn;
	    nn = task_p_indcol[m] * pncol_block + i;
	    queue[ipos]->parents->push_back(task_p[task_p_index[nn]]);
	    nn = task_p_indcol[m] * pncol_block + j;
	    queue[ipos]->parents->push_back(task_p[task_p_index[nn]]);	
	  }
	}
	else {
	  for (int m = 0; m < task_p_ptr[mmp]; m++) {
	    int nn;
	    nn = task_p_indcol[2 * m + 1] * pncol_block + i;   // lower
	    queue[ipos]->parents->push_back(task_p[task_p_index[nn]]);	
	    nn = task_p_indcol[2 * m] * pncol_block + j;       // upper
	    queue[ipos]->parents->push_back(task_p[task_p_index[nn]]);
	  }
	}
	queue[ipos]->parents->sort(compare_task_name);
	queue[ipos]->parents->unique();
	queue[ipos]->parallel_max = parallel_max2;
	queue[ipos]->parallel_id = ipos - parallel_max0 - parallel_max1;
	const int itmp = _isSym ? ((j * (j + 1)) / 2 + i) : (j * j + 2 * i);
	indcol[itmp] = ipos;
	ipos++;      // lower block
      }
      if (!_isSym) {
	const long ops = block_nrowl * block_ncoll * nrowl * 2L;
	string task_name = ("f " + to_string(j) + " " + to_string(i)
	                    + " : " + to_string(_level) + " : "
	                    + to_string(_nb));
      	DSchurGEMM_arg<T> *arg = 
	  new DSchurGEMM_arg<T>(_isSym,
				true, // isTrans
				_upper, 
				_lower,
				nrow,
				j, // transposed
				i, 
				_localSchur,
				verbose,
				fp,
				_nb);
        *(arg->ops_complexity) = ops;
        queue[ipos] = new C_task(C_DGEMM_LOCAL_MULT,
				 task_name,
				 (void *)arg,
				 DSchurGEMM_offdiag<T>,
				 1,              // atomic_size
				 0,              // atomic_id
				 arg->ops_complexity);
	queue[ipos]->fp = &(arg->fp);
	const int mmp = task_p_ptr.size() - 1;
	for (int m = 0; m < task_p_ptr[mmp]; m++) {
	  int nn;
	  nn = task_p_indcol[2 * m + 1] * pncol_block + j;   // lower
	  queue[ipos]->parents->push_back(task_p[task_p_index[nn]]);	
	  nn = task_p_indcol[2 * m] * pncol_block + i;       // upper
	  queue[ipos]->parents->push_back(task_p[task_p_index[nn]]);
	}
	queue[ipos]->parents->sort(compare_task_name);
	queue[ipos]->parents->unique();
	queue[ipos]->parallel_max = parallel_max2;
	queue[ipos]->parallel_id = ipos - parallel_max0 - parallel_max1;
	const int itmp = j * j + 2 * i + 1;
	indcol[itmp] = ipos;
	ipos++;
      } // (!_isSym)
    } // loop : i

    // diagonal
    {
      string task_name = ("f " + to_string(j) + " " + to_string(j)
			  + " : " + to_string(_level) + " : "
			  + to_string(_nb));
      const long nrowl = (long)nrow;
      const long block_ncoll = (long)_localSchur->ncolBlock(j, j);
      const long ops = (_isSym ? 
			(block_ncoll * (block_ncoll + 1L) * nrowl) : 
			(block_ncoll * block_ncoll * nrowl * 2L));
      DSchurGEMM_arg<T> *arg = 
	new DSchurGEMM_arg<T>(_isSym,
			      false, // isTrans
			      _lower,  // 
			      _upper, //
			      nrow,
			      j, //(j * SIZE_B1), 
			      (-1),
			      _localSchur,
			      verbose,
			      fp,
			      _nb);
      *(arg->ops_complexity) = ops;
      queue[ipos] = new C_task(C_DGEMM_LOCAL_MULT,
			      task_name, // task_name_cstr,
			       (void *)arg,
			       DSchurGEMM_diag<T>,
			       1,              // atomic_size
			       0,              // atomic_id
			       arg->ops_complexity);
      queue[ipos]->fp = &(arg->fp); 
      const int mmp = task_p_ptr.size() - 1;
      if (_isSym) {
	for (int m = 0; m < task_p_ptr[mmp]; m++) {
	  int nn;
	  nn = task_p_indcol[m] * pncol_block + j;
	  queue[ipos]->parents->push_back(task_p[task_p_index[nn]]);
	}
      }
      else {
	for (int m = 0; m < task_p_ptr[mmp]; m++) {
	  int nn;
	  nn = task_p_indcol[2 * m + 1] * pncol_block + j;   // lower
	  queue[ipos]->parents->push_back(task_p[task_p_index[nn]]);	
	  nn = task_p_indcol[2 * m] * pncol_block + j;       // upper
	  queue[ipos]->parents->push_back(task_p[task_p_index[nn]]);
	}
      }
      queue[ipos]->parents->sort(compare_task_name);
      queue[ipos]->parents->unique();
      queue[ipos]->parallel_max = parallel_max2;
      queue[ipos]->parallel_id = ipos - parallel_max0 - parallel_max1;	
      const int itmp = _isSym ? ((j * (j + 1)) / 2 + j) : (j * j + 2 * j);
      indcol[itmp] = ipos;
      ipos++;
    } 
  } // loop : j
  return queue.size();
}

template
int DissectionMatrix<double, double>::
C_DGEMM_local_queue(vector<C_task*> &queue,
		    vector<int> &indcol,
		    RectBlockMatrix<double> *upper1,
		    RectBlockMatrix<double> *lower1,
		    bool isSkip,
		    bool isDirect,
		    SquareBlockMatrix<double>* fdiag,
		    vector<C_task *> &task_p,
		    vector<int> &task_p_index,
		    vector<int> &task_p_indcol,
		    vector<int> &task_p_ptr,
		    vector<C_task *> &task_q,
		    vector<int> &task_q_index,
		    vector<int> &task_q_indcol,
		    vector<int> &task_q_ptr,
		    vector<C_task *> *task_s,
		    const bool verbose,
		    FILE **fp);

template
int DissectionMatrix<complex<double>, double>::
C_DGEMM_local_queue(vector<C_task*> &queue,
		    vector<int> &indcol,
		    RectBlockMatrix<complex<double> > *upper1,
		    RectBlockMatrix<complex<double> > *lower1,
		    bool isSkip,
		    bool isDirect,
		    SquareBlockMatrix<complex<double> > *fdiag,
		    vector<C_task *> &task_p,
		    vector<int> &task_p_index,
		    vector<int> &task_p_indcol,
		    vector<int> &task_p_ptr,
		    vector<C_task *> &task_q,
		    vector<int> &task_q_index,
		    vector<int> &task_q_indcol,
		    vector<int> &task_q_ptr,
		    vector<C_task *> *task_s,
		    const bool verbose,
		    FILE **fp);

template
int DissectionMatrix<quadruple, quadruple>::
C_DGEMM_local_queue(vector<C_task*> &queue,
		    vector<int> &indcol,
		    RectBlockMatrix<quadruple> *upper1,
		    RectBlockMatrix<quadruple> *lower1,
		    bool isSkip,
		    bool isDirect,
		    SquareBlockMatrix<quadruple>* fdiag,
		    vector<C_task *> &task_p,
		    vector<int> &task_p_index,
		    vector<int> &task_p_indcol,
		    vector<int> &task_p_ptr,
		    vector<C_task *> &task_q,
		    vector<int> &task_q_index,
		    vector<int> &task_q_indcol,
		    vector<int> &task_q_ptr,
		    vector<C_task *> *task_s,
		    const bool verbose,
		    FILE **fp);

template
int DissectionMatrix<complex<quadruple>, quadruple>::
C_DGEMM_local_queue(vector<C_task*> &queue,
		    vector<int> &indcol,
		    RectBlockMatrix<complex<quadruple> > *upper1,
		    RectBlockMatrix<complex<quadruple> > *lower1,
		    bool isSkip,
		    bool isDirect,
		    SquareBlockMatrix<complex<quadruple> > *fdiag,
		    vector<C_task *> &task_p,
		    vector<int> &task_p_index,
		    vector<int> &task_p_indcol,
		    vector<int> &task_p_ptr,
		    vector<C_task *> &task_q,
		    vector<int> &task_q_index,
		    vector<int> &task_q_indcol,
		    vector<int> &task_q_ptr,
		    vector<C_task *> *task_s,
		    const bool verbose,
		    FILE **fp);
//

template<typename T, typename U>
void DissectionMatrix<T, U>::C_deallocLocalSchur_queue(vector<C_task*> &queue,
						       vector<int> &indcol,
						       const bool verbose,
						       FILE **fp)
{
  const int num_block = _localSchur->num_blocks(); 
  const int num_block0 = _localSchur->num_blocks0(); 
  const int num_block1 = num_block - num_block0;
  int sizeq;
  sizeq = (_isSym ? 
	   (num_block0 * num_block1 + (num_block1 * (num_block1 + 1)) / 2) : 
	   (num_block0 * num_block1 * 2 + num_block1 * num_block1));

  queue.resize(sizeq);

  sizeq = (_isSym ? 
	   (num_block * (num_block + 1) / 2) : (num_block * num_block));
	   
  indcol.resize(sizeq);

  int ipos = 0;
  for (int j = num_block0; j < num_block; j++) {
    for (int i = 0; i < num_block0; i++) {
      {
	string task_name = ("m " + to_string(i) + " " + to_string(j)
			    + " : " + to_string(_level) + " : "
			    + to_string(_nb));
	const long nrow = (long)_localSchur->nrowBlock(i, j);
	const long ncol = (long)_localSchur->ncolBlock(i, j);
	const long ops = nrow * ncol;
	C_deallocLocalSchur_arg<T> *arg =
	  new C_deallocLocalSchur_arg<T>(_isSym,
					 _localSchur,
					 i,
					 j,
					 verbose,
					 fp,
					 _nb);
	*(arg->ops_complexity) = ops;
	queue[ipos] = new C_task(C_DEALLOCLOCALSCHUR,
				 task_name,
				 (void *)arg,
				 C_deallocLocalSchur<T>,
				 1,              // atomic_size
				 0,              // atomic_id
				 arg->ops_complexity); 
	queue[ipos]->fp = &(arg->fp);
	const int itmp = _isSym ? ((j * (j + 1)) / 2 + i) : (j * j + 2 * i);
	indcol[itmp] = ipos;
	ipos++;      // lower block
      }
      if (!_isSym) {
	string task_name = ("m " + to_string(j) + " " + to_string(i)
			    + " : " + to_string(_level) + " : "
			    + to_string(_nb));
	const long ops = 1L; // dummy
	C_deallocLocalSchur_arg<T> *arg =
	  new C_deallocLocalSchur_arg<T>(_isSym,
					 _localSchur,
					 j,
					 i,
					 verbose,
					 fp,
					 _nb);
	*(arg->ops_complexity) = ops;
	queue[ipos] = new C_task(C_DEALLOCLOCALSCHUR,
				 task_name,
				 (void *)arg,
				 C_deallocLocalSchur<T>,
				 1,              // atomic_size
				 0,              // atomic_id
				 arg->ops_complexity); 
	queue[ipos]->fp = &(arg->fp);
	const int itmp = j * j + 2 * i + 1;
	indcol[itmp] = ipos;
	ipos++;
      } // (!_isSym)
    } // loop : i
  }
  for (int j = num_block0; j < num_block; j++) {
    for (int i = num_block0; i < j; i++) {
      {
	string task_name = ("m " + to_string(i) + " " + to_string(j)
			    + " : " + to_string(_level) + " : "
			    + to_string(_nb));
	const long nrow = (long)_localSchur->nrowBlock(i, j);
	const long ncol = (long)_localSchur->ncolBlock(i, j);	
	const long ops = nrow * ncol;
	C_deallocLocalSchur_arg<T> *arg =
	  new C_deallocLocalSchur_arg<T>(_isSym,
					 _localSchur,
					 i,
					 j,
					 verbose,
					 fp,
					 _nb);
	*(arg->ops_complexity) = ops;
	queue[ipos] = new C_task(C_DEALLOCLOCALSCHUR,
				 task_name,
				 (void *)arg,
				 C_deallocLocalSchur<T>,
				 1,              // atomic_size
				 0,              // atomic_id
				 arg->ops_complexity); 
	queue[ipos]->fp = &(arg->fp);
	const int itmp = _isSym ? ((j * (j + 1)) / 2 + i) : (j * j + 2 * i);
	indcol[itmp] = ipos;
	ipos++;      // lower block
      }
      if (!_isSym) {
	string task_name = ("m " + to_string(j) + " " + to_string(i)
	                   + " : " + to_string(_level) + " : "
	                   + to_string(_nb));
	const long nrow = (long)_localSchur->nrowBlock(i, j);
	const long ncol = (long)_localSchur->ncolBlock(i, j);
	const long ops = nrow * ncol; // dummy
	C_deallocLocalSchur_arg<T> *arg =
	  new C_deallocLocalSchur_arg<T>(_isSym,
					 _localSchur,
					 j,
					 i,
					 verbose,
					 fp,
					 _nb);
	*(arg->ops_complexity) = ops;
	queue[ipos] = new C_task(C_DEALLOCLOCALSCHUR,
				 task_name,
				 (void *)arg,
				 C_deallocLocalSchur<T>,
				 1,              // atomic_size
				 0,              // atomic_id
				 arg->ops_complexity); 
	queue[ipos]->fp = &(arg->fp);
	const int itmp = j * j + 2 * i + 1;
	indcol[itmp] = ipos;
	ipos++;
      } // (!_isSym)
    } // loop : i

    // diagonal
    {
      string task_name = ("m " + to_string(j) + " " + to_string(j)
			   + " : " + to_string(_level) + " : "
			   + to_string(_nb));
      const long nrow = (long)_localSchur->nrowBlock(j, j);
      const long ops = nrow * nrow; // dummy
      C_deallocLocalSchur_arg<T> *arg =
	new C_deallocLocalSchur_arg<T>(_isSym,
				       _localSchur,
				       j,
				       j,
				       verbose,
				       fp,
				       _nb);
      *(arg->ops_complexity) = ops;
      queue[ipos] = new C_task(C_DEALLOCLOCALSCHUR,
			       task_name,
			       (void *)arg,
			       C_deallocLocalSchur<T>,
			       1,              // atomic_size
			       0,              // atomic_id
			       arg->ops_complexity); 
      queue[ipos]->fp = &(arg->fp);
      const int itmp = _isSym ? ((j * (j + 1)) / 2 + j) : (j * j + 2 * j);
      indcol[itmp] = ipos;
      ipos++;
    } 
  } // loop : j
#if 0
  for (vector<C_task *>::const_iterator it = queue.begin(); it != queue.end(); ++it) {
    fprintf(stderr, "%s : %d / %d : %d parents : ", 
	    (*it)->task_name, 
	    (*it)->parallel_id,
	    (*it)->parallel_max,
	    (*it)->parents->size());
    for (list<C_task*>::const_iterator jt = (*it)->parents->begin(); 
	 jt != (*it)->parents->end(); ++jt) {
      fprintf(stderr, "%s : ", (*jt)->task_name);
    }
    fprintf(stderr, "\n");
  }
#endif
}

template
void DissectionMatrix<double, double>::
C_deallocLocalSchur_queue(vector<C_task*> &queue,
			  vector<int> &indcol,
			  const bool verbose,
			  FILE **fp);

template
void DissectionMatrix<complex<double>, double>::
C_deallocLocalSchur_queue(vector<C_task*> &queue,
			  vector<int> &indcol,
			  const bool verbose,
			  FILE **fp);

template
void DissectionMatrix<quadruple, quadruple>::
C_deallocLocalSchur_queue(vector<C_task*> &queue,
			  vector<int> &indcol,
			  const bool verbose,
			  FILE **fp);

template
void DissectionMatrix<complex<quadruple>, quadruple>::
C_deallocLocalSchur_queue(vector<C_task*> &queue,
			  vector<int> &indcol,
			  const bool verbose,
			  FILE **fp);

//

template<typename T, typename U> 
void DissectionMatrix<T, U>::
ChildContrib(list<child_contribution<T> > *child_contribs,
	     Dissection::Tree *btree,
	     vector<DissectionMatrix<T, U>* >& dM,
	     FILE **fp)
{
#if 0
  cout << "** strips : _nb = " << _nb << endl;
#endif
  int offset_diag_src = 0;
  for (int ll = (_level - 1); ll >= 0; ll--) {
    // copy of strips with shifting position from inside of each block to 
    // continous block
    list<index_strip> diag;
    list<index_strip> offdiag;
    const int father_id = btree->nthfatherIndex(_nb, (_level - 1 - ll) + 1);
    const Dissection::SetOfStrips &diag_xj = btree->getFathersStrips(_nb)[ll];
    int father_id0 = btree->selfIndex(father_id);
    DissectionMatrix &fatherM = *dM[father_id0];

    for (Dissection::SetOfStrips::const_iterator it = diag_xj.begin();
	 it != diag_xj.end(); ++it) {
      // destination is in dense matrix, then without offest 
      // offset = (*it).begin_dst - (*itf).begin_src = (*it).begin_dst;
      // offset_offdiagf_src = 0;
      // (*itf).begin_src + offset_offdiagf_src + offset 
      //       == (*it).begin_dst
      diag.push_back(index_strip((*it).begin_dst,
				 (*it).begin_src + offset_diag_src, 
				 (*it).width));
    }
    offset_diag_src += diag_xj.numberOfIndices();
    int offset_offdiagc_src = offset_diag_src;    
    for (int m = (ll - 1); m >= 0; m--) {
      const Dissection::SetOfStrips &offdiag_xjc = 
	btree->getFathersStrips(_nb)[m];
      const Dissection::SetOfStrips &offdiag_xjf = 
	btree->getFathersStrips(father_id)[m];
#if 0
      cout << "chlid = " << _nb << " level = " << m << endl;
      for (Dissection::SetOfStrips::const_iterator it = offdiag_xjc.begin();
	   it != offdiag_xjc.end(); ++it) {
	cout << " / " << (*it).begin_dst << " : " << (*it).begin_src
	     << " : " << (*it).width;
      }
      cout << endl;
      cout << "father= " << father_id << " level = " << m << endl;
      for (Dissection::SetOfStrips::const_iterator it = offdiag_xjf.begin();
	   it != offdiag_xjf.end(); ++it) {
	cout << " / " << (*it).begin_dst << " : " << (*it).begin_src
	     << " : " << (*it).width;
      }
      cout << endl;
#endif
      //
      int offset_offdiagf_src = 0;
      for (int k = ll - 1; k > m; k--) {
	offset_offdiagf_src += 
	  btree->getFathersStrips(father_id)[k].numberOfIndices();
      }
      // strip is in increasing order with both begin_src (within condensed 
      // adderss of strips) and begin_dst (with column address in the block)
      //      Dissection::SetOfStrips::const_iterator itf = offdiag_xjf.begin();
      Dissection::SetOfStrips::const_iterator itf;
      for (Dissection::SetOfStrips::const_iterator itc = offdiag_xjc.begin();
	 itc != offdiag_xjc.end(); ++itc) {
	// find a strip in offdiag strips with father_id
	//	bool found = false;
	int offset = 0;
	for (itf = offdiag_xjf.begin(); itf != offdiag_xjf.end(); ++itf) {
	//	for (; itf != offdiag_xjf.end(); ++itf) {
	  if (((*itf).begin_dst <= (*itc).begin_dst) &&
	      ((*itf).begin_dst + (*itf).width) >= 
	      ((*itc).begin_dst + (*itc).width)) {
	    offset = (*itc).begin_dst - (*itf).begin_dst;
	    //	    found = true;
	    break;
	  }
	}
	offdiag.push_back(index_strip(((*itf).begin_src + 
				       offset_offdiagf_src +
				       offset),
				      (*itc).begin_src + 
				      offset_offdiagc_src, 
				      (*itc).width));
      }
      offset_offdiagc_src += offdiag_xjc.numberOfIndices();
    } //loop : ll

    //    cout << "** father_id = " << father_id << endl;
#if 0
    cout << "diag strips = " << diag.size();
    for (list<index_strip>::const_iterator it = diag.begin();
	 it != diag.end(); ++it) {
      cout << " / " << (*it).begin_dst << " : " << (*it).begin_src
	   << " : " << (*it).width;
    }
    cout << endl;
    cout << "offdiag strips = " << offdiag.size() << " = ";
    for (int m = (ll - 1); m >= 0; m--) {
      const Dissection::SetOfStrips &offdiag_xj = 
	btree->getFathersStrips(_nb)[m];
      cout << offdiag_xj.numberOfStrips() << " ";
    }

    int begin_src0 = (-1);
    int begin_dst0 = (-1);
    //    int begin_dst0 = (-1);
    for (list<index_strip>::const_iterator it = offdiag.begin();
	 it != offdiag.end(); ++it) {
      if ((begin_src0 > (*it).begin_src) || 
	  (begin_dst0 > (*it).begin_dst)) {
	cerr << "**ERROR**";
      }
      cout << " / " << (*it).begin_dst << " : " << (*it).begin_src
	   << " : " << (*it).width;
      begin_src0 = (*it).begin_src;
      begin_dst0 = (*it).begin_dst;
    }
    cout << endl;
#endif
    
    list<child_contribution<T> > &tmp = child_contribs[father_id0];
    if (fatherM.nrow() == 0 || _nrow == 0) {
      fprintf(*fp, "%s %d : %d/%d %d/%d %d %d:%d\n",
	      __FILE__, __LINE__,
	      _nrow, _nb, fatherM.nrow(), father_id, ll,
	      (int)diag.size(), (int)offdiag.size());
      diag.clear();
      if (_nrow == 0) {
	offdiag.clear();
      }
    }
    int child_id = btree->selfIndex(_nb);
    tmp.push_back(child_contribution<T>(child_id,
					fatherM.nrow(),
					fatherM.getUpperNCol(),
					diag, 
					offdiag,
					fatherM.addrdiagBlock(),
					fatherM.addrupperBlock(),
					(!_isSym) ? fatherM.addrlowerBlock() : (RectBlockMatrix<T> *)NULL,
					_ncol_offdiag,
					_localSchur));
  } // loop : ll
}

template
void DissectionMatrix<double, double>::
ChildContrib(list<child_contribution<double> > *child_contribs,
	     Dissection::Tree *btree,
	     vector<DissectionMatrix<double>* >& dM,
	     FILE **fp);

template
void DissectionMatrix<complex<double>, double>::
ChildContrib(list<child_contribution<complex<double> > > *child_contribs,
	     Dissection::Tree *btree,
	     vector<DissectionMatrix<complex<double>, double>* >& dM,
	     FILE **fp);

template
void DissectionMatrix<quadruple, quadruple>::
ChildContrib(list<child_contribution<quadruple> > *child_contribs,
	     Dissection::Tree *btree,
	     vector<DissectionMatrix<quadruple>* >& dM,
	     FILE **fp);

template
void DissectionMatrix<complex<quadruple>, quadruple>::
ChildContrib(list<child_contribution<complex<quadruple> > > *child_contribs,
	     Dissection::Tree *btree,
	     vector<DissectionMatrix<complex<quadruple>, quadruple>* >& dM,
	     FILE **fp);
//

template<typename T, typename U> 
void DissectionMatrix<T, U>::deallocLower_queue(C_task*& queue,
						bool isSym,
						vector<C_task*> &task_p,
						bool isDirect,
						vector<C_task*> &task_q)
{
  string task_name = "M  : " + to_string(_nb);
  
  //  char *task_name_cstr = new char[task_name.str().size() + 1];
  //  strcpy(task_name_cstr, task_name.str().c_str());
  
  const int ops_complexity = _ncol_offdiag * _nrow;
  C_deallocLower_arg<T> *arg = new C_deallocLower_arg<T>(isSym,
							 _lower,
							 ops_complexity);
					       
  queue = new C_task(C_DEALLOCLOWER,
		     task_name,
		     (void *)arg,
		     C_deallocLower<T>,
		     1,  // atomic_size
		     0,  // atomic_id
		     arg->ops_complexity);
  for (vector<C_task *>::const_iterator it = task_p.begin(); 
       it != task_p.end(); ++it) {
    queue->parents->push_back(*it);
  }
  if (isDirect) {
    for (int i = 0; i < task_q[0]->parallel_max; i++) {
      queue->parents->push_back(task_q[i]);
    }
    queue->parents->sort(compare_task_name);
    queue->parents->unique();  
  }
}

template
void DissectionMatrix<double, double>::
deallocLower_queue(C_task*& queue,
		   bool isSym,
		   vector<C_task*> &task_p,
		   bool isDirect,
		   vector<C_task*> &task_q);

template
void DissectionMatrix<complex<double>, double>::
deallocLower_queue(C_task*& queue,
		   bool isSym,
		   vector<C_task*> &task_p,
		   bool isDirect,
		   vector<C_task*> &task_q);

template
void DissectionMatrix<quadruple, quadruple>::
deallocLower_queue(C_task*& queue,
		   bool isSym,
		   vector<C_task*> &task_p,
		   bool isDirect,
		   vector<C_task*> &task_q);

template
void DissectionMatrix<complex<quadruple>, quadruple>::
deallocLower_queue(C_task*& queue,
		   bool isSym,
		   vector<C_task*> &task_p,
		   bool isDirect,
		   vector<C_task*> &task_q);
//
