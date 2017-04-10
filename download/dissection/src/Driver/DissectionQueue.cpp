/*! \file   DissectionQueue.hpp
    \brief  management of threads for factorization and Fw/Bw substitution
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
#include "Driver/C_Dsub.hpp"
#include "Driver/DissectionQueue.hpp"
#include "Driver/QueueRuntime.hpp"

#include <string>

#define RATIO_DTRSM_MERGED 0.4

template<typename T, typename U>
const T DissectionQueue<T, U>::_one = T(1.0);
template<typename T, typename U>
const T DissectionQueue<T, U>::_none = T(-1.0);
template<typename T, typename U>
const T DissectionQueue<T, U>::_zero = T(0.0);


template<typename T, typename U>
DissectionQueue<T, U>::DissectionQueue(Dissection::Tree *btree,
				       vector<DissectionMatrix<T, U>* >& dM,
				       const int num_threads,
				       const bool isSym,
				       const bool verbose,
				       FILE *fp) :
  _btree(btree), _num_threads(num_threads), 
  _verbose(verbose), _fp(fp), _isSym(isSym),
  _queue_symb_allocated(false), _queue_numrc_allocated(false),
  _queue_fwbw_allocated(false)
{
  const int nb_level = _btree->NumberOfLevels();
  _nb_level = nb_level;
  const int level_last = nb_level - 1;
  const int nb_doms = _btree->NumberOfSubdomains(); //(1U<<(level_last+1))-1
  _nb_doms = nb_doms;
  const int nb_doms_dense = (1U << level_last) - 1;
  const int nb_doms_sparse = (1U << level_last);

  //  _queue_dynamic->reserve(10);

  for (int d = 1; d <= nb_doms; d++) {
    const int d1 = _btree->selfIndex(d);
    dM[d1] = new DissectionMatrix<T, U>(_btree, d, isSym, verbose, fp); 
  } 

  // set up working arrays for C_Dsub
  _children = new vector<int>[nb_doms];
  _tasks_SparseSymb = new vector<C_task *>[nb_doms_sparse];
  _tasks_SparseNum = new vector<C_task *>[nb_doms];        // nd_doms_sparse is 
  _tasks_SparseLocalSchur = new vector<C_task *>[nb_doms]; // enough but
  _tasks_DFillSym = new vector<C_task *>[nb_doms];         // for convinience of
                                                           // dependency
  _tasks_DFullLDLt = new vector<C_task *>[nb_doms_dense];
  _tasks_DTRSMScale = new vector<C_task *>[nb_doms_dense];
  _tasks_DSymmGEMM = new vector<C_task *>[nb_doms_dense];
  _tasks_deallocLocalSchur = 
    new vector<C_task *>[nb_doms_dense + nb_doms_sparse];
  _tasks_Dsub = new vector<C_task *>*[nb_level];
  _tasks_Dsub[0] = new vector<C_task *>[nb_level * nb_doms_dense];
  for (int i = 1; i < nb_level; i++) {
    _tasks_Dsub[i] = _tasks_Dsub[0] + i * nb_doms_dense;
  }
  _tasks_deallocLower = new vector<C_task *>[nb_level];
  _child_contribs = new list<child_contribution<T> >[nb_doms];
 
  _pivots = new double[nb_doms];  // pass pivot value from a DissectionMatrix to
                                  // the other is not easy in thread execution
  // initialize state to controle tasks

  _dissectionRuntime = new QueueRuntime(nb_doms, num_threads, verbose, fp);

#if 0 // move to QueueRuntime
  _zone_entered = new int[DIST_TASK_CRITICAL];
  _zone_finished = new int[DIST_TASK_CRITICAL];
  _zone_static_assigned = new int[DIST_TASK_CRITICAL];

  allocate_int2d(_begins, num_threads);
  allocate_int2d(_ends, num_threads);

  allocate_int3d(_group_entered, num_threads);
  allocate_int3d(_group_finished, num_threads);
  allocate_int3d(_group_task_ends, num_threads);
  allocate_int3d(_group_static_assigned, num_threads);

  _group_task_id = new int[num_threads];

  allocate_int3d(_begins_group, num_threads);
  allocate_int3d(_ends_group, num_threads);
  allocate_unsigned2d(_group_nops, num_threads);

  _mutex_group = new pthread_mutex_t[num_threads];
#endif
}

template
DissectionQueue<double, double>::
DissectionQueue(Dissection::Tree *btree,
		vector<DissectionMatrix<double>* >& dM,
		const int num_threads,
		const bool isSym,
		const bool verbose,
		FILE *fp);
template
DissectionQueue<quadruple, quadruple>::
DissectionQueue(Dissection::Tree *btree,
		vector<DissectionMatrix<quadruple>* >& dM,
		const int num_threads,
		const bool isSym,
		const bool verbose,
		FILE *fp);
template
DissectionQueue<complex<double>, double>::
DissectionQueue(Dissection::Tree *btree,
		vector<DissectionMatrix<complex<double>, double>* >& dM,
		const int num_threads,
		const bool isSym,
		const bool verbose,
		FILE *fp);

template
DissectionQueue<complex<quadruple>, quadruple>::
DissectionQueue(Dissection::Tree *btree,
		vector<DissectionMatrix<complex<quadruple>, quadruple>* >& dM,
		const int num_threads,
		const bool isSym,
		const bool verbose,
		FILE *fp);
//

template<typename T, typename U>
DissectionQueue<T, U>::~DissectionQueue() 
{
  //  const int nb_level = _btree->NumberOfLevels();
  //  const int num_threads = _num_threads;
#if 0
  for (list<C_task *>::iterator it = _queue_dummy.begin();
       it != _queue_dummy.end(); ++it) {
    fprintf(stderr, "%s %d : %s\n", __FILE__, __LINE__, (*it)->task_name);
    C_dummy_arg *arg = (C_dummy_arg *)(*it)->func_arg;
    delete arg;
    (*it)->func_arg = NULL;
    delete (*it);
    (*it) = NULL;
  }
  _queue_dummy.clear();
#endif
  erase_queue();
  erase_queue_fwbw();
  // set up working arrays for C_Dsub
  delete [] _children;
  delete [] _tasks_SparseSymb;
  delete [] _tasks_SparseNum;
  delete [] _tasks_SparseLocalSchur;
  delete [] _tasks_DFillSym; // _tasks_DFillSym = new vector<C_task *>[nb_doms];
  delete [] _tasks_DFullLDLt;
  delete [] _tasks_DTRSMScale;
  delete [] _tasks_DSymmGEMM;
  delete [] _tasks_Dsub[0];
  delete [] _tasks_Dsub;
  delete [] _tasks_deallocLocalSchur;
  delete [] _tasks_deallocLower;
  //  delete [] _tasks_fillSym;
  delete [] _child_contribs;
  delete [] _pivots;
  // nullify the pointer
  _btree = NULL;
  delete _dissectionRuntime;
}

template
DissectionQueue<double, double>::~DissectionQueue();

template
DissectionQueue<quadruple, quadruple>::~DissectionQueue();

template
DissectionQueue<complex<double>, double>::~DissectionQueue();

template
DissectionQueue<complex<quadruple>, quadruple>::~DissectionQueue();
//

template<typename T, typename U> 
void DissectionQueue<T, U>::
generate_queue(vector<DissectionMatrix<T, U>* >& dM, 
	       int nnz,
	       T *coefs)
{
  const int num_threads = _num_threads;
  const int nb_level = _btree->NumberOfLevels();
  const int level_last = nb_level - 1;
  
  const int nb_doms = _btree->NumberOfSubdomains();
  const int nb_doms_dense = (1U << level_last) - 1;
  const int nb_doms_sparse = (1U << level_last);
  // set up working arrays for C_Dsub

  for (int level = (level_last - 1); level >= 0; level--){
    const int begdom = 1U << level;
    const int enddom = begdom * 2;
    for (int d = begdom; d < enddom; d++) {
      const int j = _btree->selfIndex(d);
      const int mm = _btree->childIndex(d);
      const int nn = _btree->brotherIndex(mm);
      _children[j].push_back(_btree->selfIndex(mm));
      _children[j].push_back(_btree->selfIndex(nn));
    }
  }

  vector<vector<int> > all_fathersIndex(nb_level);

  long **nops_queue = new long*[num_threads + 1];
  nops_queue[0] = new long[(num_threads + 1) * nb_level];
  for (int p = 1; p <= num_threads; p++) {
    nops_queue[p] = nops_queue[0] + p * nb_level;
    for (int j = 0; j <= level_last; j++) {
      nops_queue[p][j] = 0L;
    }
  }
  vector<int>* tasks_deallocLocalSchur_indcol = 
    new vector<int>[nb_doms_dense + nb_doms_sparse];
  // sprase subdomains : begin : level == level_last
  {
    all_fathersIndex[level_last].clear();
    for (int i = 0; i < nb_doms; i++) {
      _child_contribs[i].clear();
    }

    const int begdom = 1U << level_last;
    const int enddom = begdom * 2;

    int k = 0;
    for (int d = begdom; d < enddom; d++, k++) {
      const int j = _btree->selfIndex(d);
      _tasks_SparseSymb[k].resize(1);
      _tasks_SparseNum[j].resize(1);
      _tasks_SparseLocalSchur[j].resize(1);
      //      _tasks_DFillSym[j].resize(2);

      dM[j]->C_SparseSymbFact_queue(_tasks_SparseSymb[k],
				    _btree,
				    _verbose,
				    &_fp);
      dM[j]->C_SparseNumFact_queue(_tasks_SparseNum[j], 
				   _btree,
				   nnz, 
				   coefs, 
				   &_eps_piv,
				   &_pivots[j],
				   &_kernel_detection,
				   &_aug_dim,
				   &_eps_machine,
				   _tasks_SparseSymb[k],
				   _verbose,
				   &_fp);
      dM[j]->C_SparseLocalSchur_queue(_tasks_SparseLocalSchur[j], 
				      _btree,
				      nnz, 
				      coefs,
				      _tasks_SparseNum[j],
				      _verbose,
				      &_fp);

      dM[j]->C_deallocLocalSchur_queue(_tasks_deallocLocalSchur[j],
				       tasks_deallocLocalSchur_indcol[j],
				       _verbose,
				       &_fp);
      
      dM[j]->ChildContrib(_child_contribs, 
			  _btree,
			  dM,
			  &_fp);
    } // loop : d
    
    for (int d = 1; d < begdom; d++) {
      const int j = _btree->selfIndex(d);
      _tasks_DFillSym[j].resize(2);
      dM[j]->C_FillMatrix_queue(_tasks_DFillSym[j],
				nnz, 
				coefs,
				_verbose,
				&_fp); 
    }
    all_fathersIndex[level_last].resize((1U << level_last) - 1);
    for (int i = 1; i < (1U << level_last); i++) {
      const int j = _btree->selfIndex(i);
      all_fathersIndex[level_last][j] = j; 
    }
    for (vector<int>::const_iterator it = all_fathersIndex[level_last].begin();
	 it != all_fathersIndex[level_last].end(); ++it) {
      
      C_Dsub_queue<T>(_isSym, (*it), 
		      false, // all sparse matrix parts need to be done
		      _tasks_Dsub[level_last][(*it)], 
		      _child_contribs[(*it)],
		      (vector<C_task*> *)NULL, // _tasks_DSymmGEMM
		      (vector<int> *)NULL,
		      false,
		      _tasks_DFillSym,
		      _tasks_SparseLocalSchur, 
		      (vector<C_task*> *)NULL,
		      _tasks_deallocLocalSchur,
		      tasks_deallocLocalSchur_indcol,
		      level_last, _verbose, _fp);
    }
    // tasks to deallocate local Schur complement
    _tasks_deallocLower[level_last].resize(0); // no C- working array of DTSRM
  }  // end : level == level_last 
  vector<int> nrow_DFullLDLt;
  nrow_DFullLDLt.resize(nb_doms_dense);
  vector<bool> isMergedDTRSM(nb_doms_dense, false);
  vector<bool> isDividedDTRSM(nb_doms_dense, false);
  vector<int>* tasks_DTRSMScale_indcol = new vector<int>[nb_doms_dense];
  vector<int>* tasks_DTRSMScale_index = new vector<int>[nb_doms_dense];
  vector<list<int> >* tasks_DTRSMScale_parents_index = new vector<list<int> >[nb_doms_dense];
  vector<int>* tasks_DTRSMScale_ptr = new vector<int>[nb_doms_dense];
  vector<int>* tasks_DSymmGEMM_indcol = new vector<int>[nb_doms_dense];

  for (int level = (level_last - 1); level >= 0; level--) {
    int itmp;
    all_fathersIndex[level].clear();
    for (int i = 0; i < nb_doms; i++) {
      _child_contribs[i].clear();
    }
    const int begdom = 1U << level;
    const int enddom = begdom * 2;
    for (int d = begdom; d < enddom; d++) {
      const int j = _btree->selfIndex(d);

      const int num_row_b = dM[j]->diagBlock().num_blocks();
      int num_tasks_ldlt;
      num_tasks_ldlt = (num_row_b * (num_row_b + 1) * (num_row_b + 2)) / 6;
      if (level == 0) {
	num_tasks_ldlt += (((num_row_b - 1) * num_row_b * (num_row_b + 1)) / 6 +
			   num_row_b);
      }
      nrow_DFullLDLt[j] = num_row_b;
      // factorization by block + whole pivotting
      
      _tasks_DFullLDLt[j].resize(num_tasks_ldlt); 
      vector<int> task_indcol;
      vector<int> task_ptr;
      task_indcol.resize(num_tasks_ldlt);
      //      _tasks_DTRSMScale[j].resize(num_col_b2);
      //      _tasks_DSymmGEMM[j].resize(num_dgemm);
    
      bool isChldrnAlgnd = false;
      int evnchld_id;
      if (level < (level_last - 1)) {
	evnchld_id = _btree->childIndex(d); // "selfIndex(j) + 1 == j" is used
	isChldrnAlgnd = dM[evnchld_id]->isAlignedFather();
      }
      vector<C_task *> &task_p = (isChldrnAlgnd ? _tasks_DSymmGEMM[evnchld_id] :
				  _tasks_Dsub[level + 1][j]);
      
      itmp = dM[j]->C_DFullLDLt_queue(_tasks_DFullLDLt[j],
				      task_indcol,
				      task_ptr,
				      // value not yet defined
				      &_eps_piv,   
				      &_kernel_detection,
				      &_aug_dim,
				      &_eps_machine,
				      &_pivots[j],
				      &_pivots[_children[j][0]],
				      &_pivots[_children[j][1]],
				      task_p,
				      isChldrnAlgnd,
				      _verbose,
				      &_fp);
      if (itmp > 1) {
	itmp = EraseNullParents(_tasks_DFullLDLt[j]);
      }
#ifdef DEBUG_PREPARE_THREAD
      cout << j << " : "
	   << " generated queue : C_DFullLDLt = " << j
	   << " num_block = " << num_row_b
	   << " num_tasks = " << num_tasks_ldlt
	   << " / " << itmp << endl;
#endif
      if (level > 0) {
	vector<C_task *> &task_pp = _tasks_Dsub[level + 1][j];
	
	itmp = dM[j]->C_DTRSMScale_queue(_tasks_DTRSMScale[j],
				  tasks_DTRSMScale_parents_index[j],
				  tasks_DTRSMScale_index[j],
				  tasks_DTRSMScale_indcol[j],
				  tasks_DTRSMScale_ptr[j],
				  _btree,
				  _tasks_DFullLDLt[j],
				  task_indcol,
				  task_ptr,
				  task_pp,
				  _verbose,
				  &_fp);
	if (itmp > 1) {
	  itmp = EraseNullParents(_tasks_DTRSMScale[j]);
	  dM[j]->C_DTRSMScale_rearrange(_tasks_DTRSMScale[j],
					tasks_DTRSMScale_parents_index[j],
					tasks_DTRSMScale_index[j],
					tasks_DTRSMScale_indcol[j],
					tasks_DTRSMScale_ptr[j],
					_verbose,
					_fp);
	}
      }
#ifdef DEBUG_PREPARE_THREAD
      cout << j << " : "
	   << " generated queue : C_DTRSMScale_queue " << j
	   << " num_tasks = " << num_col_b2 << " " 
	   << _tasks_DTRSMScale[j].size() << endl;
#endif

      SquareBlockMatrix<T>* f_diag;
      int jf = (-1);
      if (level > 0) { // there is no father for level == 0 
	jf = _btree->selfIndex(_btree->fatherIndex(d));
	f_diag = dM[jf]->addrdiagBlock();
      }
      else {
	f_diag = (SquareBlockMatrix<T>*)NULL;
      }
      bool isDirect = dM[j]->isAlignedFather();
      bool isEven = (level > 0) ? (j % 2 == 0) : false;
      RectBlockMatrix<T> *upper1, *lower1;
      upper1 = (isEven ? dM[j - 1]->addrupperBlock() :
		(RectBlockMatrix<T> *) NULL);
      lower1 = (isEven ? dM[j - 1]->addrlowerBlock() :
		(RectBlockMatrix<T> *) NULL);
      vector<C_task *> null_task;
      vector<C_task *> &tasks_q = (isEven ? _tasks_DTRSMScale[j - 1] :
				   null_task);
      vector<int> null_idx;
      vector<int> &tasks_q_index = (isEven ? tasks_DTRSMScale_index[j - 1] :
				   null_idx);
      
      vector<int> &tasks_q_indcol = (isEven ? tasks_DTRSMScale_indcol[j - 1] :
				     null_idx);
      vector<int> &tasks_q_ptr = (isEven ? tasks_DTRSMScale_ptr[j - 1] :
				     null_idx);

      itmp = dM[j]->C_DGEMM_local_queue(_tasks_DSymmGEMM[j],
				 tasks_DSymmGEMM_indcol[j],
				 upper1,
				 lower1,
				 (!isEven),
				 isDirect,
				 f_diag,
				 _tasks_DTRSMScale[j],
				 tasks_DTRSMScale_index[j],
				 tasks_DTRSMScale_indcol[j],
				 tasks_DTRSMScale_ptr[j],
				 tasks_q,
				 tasks_q_index,
				 tasks_q_indcol,
				 tasks_q_ptr,
				 (jf >= 0) ? &_tasks_Dsub[level + 1][jf] : 
				 (vector<C_task *>*)NULL,
				 _verbose,
				 &_fp);
      if (itmp > 1 ) {
	itmp = EraseNullParents(_tasks_DSymmGEMM[j]);
      }
      dM[j]->C_deallocLocalSchur_queue(_tasks_deallocLocalSchur[j],
				       tasks_deallocLocalSchur_indcol[j],
				       _verbose,
				       &_fp);
				    

#ifdef DEBUG_PREPARE_THREAD
      cout << j << " : "
	   << " generated queue : C_DGEMM_local_queue " << j
	   << " num_tasks = " << num_dgemm << " " 
	   << _tasks_DGEMM[j].size() << endl;
#endif
      
      dM[j]->ChildContrib(_child_contribs, 
			  _btree,
			  dM,
			  &_fp);
	
    } // loop : d(j)
    
    all_fathersIndex[level].resize((1U << level) - 1);
    for (int i = 1; i < (1U << level); i++) {
      const int j = _btree->selfIndex(i);
      all_fathersIndex[level][j] = j;
    }
    for (vector<int>::const_iterator it = all_fathersIndex[level].begin();
	 it != all_fathersIndex[level].end(); ++it) {
      bool skip_flag = false;
      const int father_id = _btree->Index2Node((*it));
      if (_btree->nodeLayer(father_id) == (level - 1)) {
	if (_child_contribs[(*it)].size() == 2) {
	  int child_id = _btree->childIndex(father_id);
	  int child_id0 = _btree->selfIndex(child_id);
	  int child_id1 = _btree->selfIndex(child_id + 1);
	  skip_flag = ((dM[child_id0]->isAlignedFather() &&
			dM[child_id1]->isAlignedFather()) ? 
		       true : false);
	  if (_verbose) {
	    fprintf(_fp,
		    "%s %d : father_id = %d children_id = %d %d skip_flag = %s\n",
		    __FILE__, __LINE__, 
		    ((*it) + 1), (child_id0 + 1), (child_id1 + 1), 
		    skip_flag ? "true" : "false");
	  }
	}
      }
      C_Dsub_queue<T>(_isSym, (*it),
		      skip_flag,
		      _tasks_Dsub[level][(*it)], 
		      _child_contribs[(*it)],
		      _tasks_DSymmGEMM,
		      tasks_DSymmGEMM_indcol,
		      true,
		      (vector<C_task*> *)NULL, // _tasks_DFillSym,
		      (vector<C_task*>*)NULL,
		      &_tasks_Dsub[level + 1][(*it)],
		      _tasks_deallocLocalSchur,
		      tasks_deallocLocalSchur_indcol,
		      level,
		      _verbose,
		      _fp);
      // itmp = EraseNullParents(_tasks_Dsub[level][(*it)]);
    }
    // tasks to deallocate lower column matrix and local Schur complement
    if (level > 0) {
      _tasks_deallocLower[level].resize(enddom - begdom);
      {
	int k = 0;
	for (int d = begdom; d < enddom; d++, k++) {
	  const int j = _btree->selfIndex(d);
	  bool isDirect = dM[j]->isAlignedFather();
	  const int jj = (j % 2 == 1) ? (j + 1) : j;
	  dM[j]->deallocLower_queue(_tasks_deallocLower[level][k],
				    _isSym,
				    _tasks_DSymmGEMM[j],
				    isDirect,
				    _tasks_DSymmGEMM[jj]);

	} // loop : d, k
      } // scope for k
      itmp = EraseNullParents(_tasks_deallocLower[level]);
    }  // if (level > 0)
  } // loop : level
#if 0
  if (_verbose) {
    for (int level = level_last; level > 0; level--) {
      for (vector<C_task *>::const_iterator it = _tasks_deallocLower[level].begin();
	   it != _tasks_deallocLower[level].end(); ++it) {
	fprintf(stderr, "%s : parents = %d : ", 
		(*it)->task_name, (*it)->parents->size());
	for (list<C_task *>::const_iterator jt = (*it)->parents->begin();
	     jt != (*it)->parents->end(); ++jt) {
	  fprintf(stderr, "%s / ", (*jt)->task_name);
	}
	fprintf(stderr, "\n");
      }
    }
    
    for (int level = level_last; level > 0; level--) {
      const int begdom = 1U << level;
      const int enddom = begdom * 2;
      for (int d = begdom; d < enddom; d++) {
	const int j = _btree->selfIndex(d);
	fprintf(stderr, "nb = %d\n", j);
	for (vector<C_task *>::const_iterator it = _tasks_deallocLocalSchur[j].begin();
	     it != _tasks_deallocLocalSchur[j].end(); ++it) {
	  fprintf(stderr, "%s : parents = %d : ", 
		  (*it)->task_name, (*it)->parents->size());
	  for (list<C_task *>::const_iterator jt = (*it)->parents->begin();
	       jt != (*it)->parents->end(); ++jt) {
	    fprintf(stderr, "%s / ", (*jt)->task_name);
	  }
	  fprintf(stderr, "\n");
	}
      } // loop : d
    } // loop : level
  }
#endif
  for (int level = (level_last - 1); level > 0; level--) {
    const int begdom = 1U << level;
    const int enddom = begdom * 2;
    for (int d = begdom; d < enddom; d++) {
      const int j = _btree->selfIndex(d);
      if(nrow_DFullLDLt[j] > 2) {
	list <C_task *> tmp1;
	list <C_task *> tmp2;
	const int ncol = dM[j]->upperBlock().num_blocks_c();
	const int ncol1 = (int)((double)ncol * RATIO_DTRSM_MERGED);
	int i0, i1, k1, kk1, kk2;	//	int i0, i1, k1, k2, kk1, kk2;
	tmp1.push_back(_tasks_DFullLDLt[j][0]);
	k1 = 1;
	i0 = 1;
	i1 = 0;
	bool first_loop_done = false;
	kk2 = 0;
	while (i0 < _tasks_DFullLDLt[j].size()) {
	  const int parallel_max0 =_tasks_DFullLDLt[j][i0]->parallel_max;
	  const int parallel_max1 =_tasks_DTRSMScale[j][i1]->parallel_max;
	  kk1 = 0;
	  for (int m = 0; m < parallel_max0; m++) {
	    const int atomic_size = _tasks_DFullLDLt[j][i0]->atomic_size;
	    for (int n = 0; n < atomic_size; n++) {
	      _tasks_DFullLDLt[j][i0 + n]->parallel_id = kk1;
	      _tasks_DFullLDLt[j][i0 + n]->parallel_max = 0; // for conuting
	      tmp1.push_back(_tasks_DFullLDLt[j][i0 + n]);
	    }
	    kk1++;
	    i0 += atomic_size;
	  } // loop : m
 	  int ll = 0;
	  for (int m = 0; m < parallel_max1; m++) {
	    const int atomic_size = _tasks_DTRSMScale[j][i1]->atomic_size;
	    for (int n = 0; n < atomic_size; n++) {
	      if (first_loop_done) {
		if (ll < ncol1) {
		  _tasks_DTRSMScale[j][i1 + n]->parallel_id = kk1;
		  _tasks_DTRSMScale[j][i1 + n]->parallel_max = 0;
		  tmp1.push_back(_tasks_DTRSMScale[j][i1 + n]);
		}
		else {
		  _tasks_DTRSMScale[j][i1 + n]->parallel_id = kk2;
		  _tasks_DTRSMScale[j][i1 + n]->parallel_max = 0;
		  tmp2.push_back(_tasks_DTRSMScale[j][i1 + n]);
		}
	      }
	      else {
		_tasks_DTRSMScale[j][i1 + n]->parallel_id = kk1;
		_tasks_DTRSMScale[j][i1 + n]->parallel_max = 0;
		tmp1.push_back(_tasks_DTRSMScale[j][i1 + n]);
	      }
	    } // loop : n
	    if (first_loop_done) {
	      if (ll < ncol1) {
		kk1++;
	      }
	      else {
		kk2++;
	      }
	    }
	    else {
	      kk1++;
	    }
 	    ll++;
	    if (ll == ncol) {
	      ll = 0;
	    }
	    i1 += atomic_size;
	  } // loop : m
	  for (list<C_task *>::const_iterator it = tmp1.begin(); it != tmp1.end(); 
	       ++it) {
	    if ((*it)->parallel_max == 0) { 
	      (*it)->parallel_max = kk1;
	    }
	  }
#if 0
	  for (list<C_task *>::const_iterator it = tmp2.begin(); it != tmp2.end(); 
	       ++it) {
	    if ((*it)->parallel_max == 0) { 
	      (*it)->parallel_max = kk2;
	    }
	  }
#endif
	  if (i0 >= _tasks_DFullLDLt[j].size()) {
	    break;
	  }
	  // number of C_dupdateb_Schur_{diag,offdiag} may be large enough
	  const int parallel_max2 =_tasks_DFullLDLt[j][i0]->parallel_max;
	  kk1 = 0;
	  for (int m = 0; m < parallel_max2; m++) {
	    const int atomic_size = _tasks_DFullLDLt[j][i0]->atomic_size;
	    for (int n = 0; n < atomic_size; n++) {
	      _tasks_DFullLDLt[j][i0 + n]->parallel_id = kk1;
	      _tasks_DFullLDLt[j][i0 + n]->parallel_max = 0;
	      tmp1.push_back(_tasks_DFullLDLt[j][i0 + n]);
	    }
	    i0 += atomic_size;
	    kk1++;
	  } // loop : m 
	  for (list<C_task *>::const_iterator it = tmp1.begin(); it != tmp1.end(); 
	       ++it) {
	    if ((*it)->parallel_max == 0) { 
	      (*it)->parallel_max = kk1;
	    }
	  }
	  first_loop_done = true;
	} // while (i0 < _tasks_DFullLDLt[j].size())
	int ll = 0;
	while (i1 < _tasks_DTRSMScale[j].size()) {
	  kk1 = 0;
	  const int parallel_max =_tasks_DTRSMScale[j][i1]->parallel_max;
	  for (int m = 0; m < parallel_max; m++) {
	    const int atomic_size = _tasks_DTRSMScale[j][i1]->atomic_size;
	    for (int n = 0; n < atomic_size; n++) {
	      if (ll < ncol1) {
		_tasks_DTRSMScale[j][i1 + n]->parallel_id = kk1;
		_tasks_DTRSMScale[j][i1 + n]->parallel_max = 0;
		tmp1.push_back(_tasks_DTRSMScale[j][i1 + n]);
	      }
	      else {
		_tasks_DTRSMScale[j][i1 + n]->parallel_id = kk2;
		_tasks_DTRSMScale[j][i1 + n]->parallel_max = 0;
		tmp2.push_back(_tasks_DTRSMScale[j][i1 + n]);
	      }
	    }
	    if (ll < ncol1) {
	      kk1++;
	    }
	    else {
	      kk2++;
	    }
	    i1 += atomic_size;
	    ll++;
	    if (ll == ncol) {
	      ll = 0;
	    }
	  } // loop : m
	  for (list<C_task *>::const_iterator it = tmp1.begin(); it != tmp1.end(); 
	       ++it) {
	    if ((*it)->parallel_max == 0) { 
	      (*it)->parallel_max = kk1;
	    }
	  }
	  for (list<C_task *>::const_iterator it = tmp2.begin(); it != tmp2.end(); 
	       ++it) {
	    if ((*it)->parallel_max == 0) { 
	      (*it)->parallel_max = kk2;
	    }
	  }
	} // while (i1 < _tasks_DTRSMScale[j].size())
	_tasks_DFullLDLt[j].resize(tmp1.size());
	vector<C_task *>::iterator jt;
	jt = _tasks_DFullLDLt[j].begin();
	for (list<C_task *>::const_iterator it = tmp1.begin(); it != tmp1.end(); 
	     ++it, ++jt) {
	  (*jt) = (*it);
	}
	tmp1.clear();
	_tasks_DTRSMScale[j].clear();
	_tasks_DTRSMScale[j].resize(tmp2.size());
	jt = _tasks_DTRSMScale[j].begin();
	for (list<C_task *>::const_iterator it = tmp2.begin(); it != tmp2.end(); 
	     ++it, ++jt) {
	  (*jt) = (*it);
	}
#if 0
	if (_verbose) {
	  fprintf(stderr, "%s %d\n", __FILE__, __LINE__);
	  for (vector<C_task *>::const_iterator it = _tasks_DFullLDLt[j].begin();
	       it != _tasks_DFullLDLt[j].end(); 
	       ++it) {
	    fprintf(stderr, "%s : %d / %d : %d / %d\n", 
		    (*it)->task_name, (*it)->atomic_id, (*it)->atomic_size,
		    (*it)->parallel_id, (*it)->parallel_max);
	  }
	  fprintf(stderr, "%s %d\n", __FILE__, __LINE__);
	  for (vector<C_task *>::const_iterator it = _tasks_DTRSMScale[j].begin();
	       it != _tasks_DTRSMScale[j].end(); 
	       ++it) {
	    fprintf(stderr, "%s : %d / %d : %d / %d\n", 
		    (*it)->task_name, (*it)->atomic_id, (*it)->atomic_size,
		    (*it)->parallel_id, (*it)->parallel_max);
	  }
	}
#endif
	
	isMergedDTRSM[j] = true;
	isDividedDTRSM[j] = true;
      } //   if(nrow_DFullLDLt[j] > 2)
      else {
	if (nrow_DFullLDLt[j + ((j % 2 == 0 ) ? (-1) : 1)] > 2) {
	  // brother has merged LDLt and DTRSM
	  int size0 = _tasks_DFullLDLt[j].size();
	  int size1 = _tasks_DTRSMScale[j].size();
	  vector <C_task *> tmp;
	  tmp.resize(size0 + size1);
	  int ii = 0;
	  for (int k = 0; k < _tasks_DFullLDLt[j].size(); k++, ii++) {
	    tmp[ii] = _tasks_DFullLDLt[j][k];
	    tmp[ii]->parallel_max += size1;
	  }
	  for (int k = 0; k <  _tasks_DTRSMScale[j].size(); k++, ii++) {
	    tmp[ii] = _tasks_DTRSMScale[j][k];
	    tmp[ii]->parallel_max += size0;
	    tmp[ii]->parallel_id += size0;
	  }
#if 0
	  if (_verbose) {
	    fprintf(stderr, "%s %d\n", __FILE__, __LINE__);
	    for (vector<C_task *>::const_iterator it = tmp.begin(); it != tmp.end(); 
		 ++it) {
	      fprintf(stderr, "%s : %d / %d : %d / %d\n", 
		      (*it)->task_name, (*it)->atomic_id, (*it)->atomic_size,
		      (*it)->parallel_id, (*it)->parallel_max);
	    }
	  }
#endif
	  _tasks_DFullLDLt[j].resize(tmp.size());
	  vector<C_task *>::iterator jt = _tasks_DFullLDLt[j].begin();
	  for (vector<C_task *>::const_iterator it = tmp.begin(); it != tmp.end(); 
	       ++it, ++jt) {
	    (*jt) = (*it);
	  }
	  _tasks_DTRSMScale[j].clear();
	  isMergedDTRSM[j] = true;
	} // if 
      }// else  if(nrow_DFullLDLt[j] > 2) 
    } // loop : d
  }

  list<C_task *> queue_null;
  _dissectionRuntime->generate_queue(_queue_symb,
				     _queue_static,
				     _queue_dynamic,
				     queue_null,
     //				     _queue_dummy,
				     _btree,
				     _children,
				     _tasks_SparseSymb,
				     _tasks_SparseNum,
				     _tasks_SparseLocalSchur,
				     _tasks_DFillSym,
				     _tasks_DFullLDLt,
				     _tasks_DTRSMScale,
				     _tasks_DSymmGEMM,
				     _tasks_Dsub,
				     _tasks_deallocLower,
				     _tasks_deallocLocalSchur,
				     nops_queue,
				     all_fathersIndex,
				     nrow_DFullLDLt,
				     isMergedDTRSM,
				     isDividedDTRSM,
				     level_last);

  queue_null.sort(compare_task_name);
  queue_null.unique();
#if 0
  _queue_dummy.sort(compare_task_name); // to avoid multiple copies for safty
  _queue_dummy.unique();
#endif
// #define DEBUG_QUEUE_NULL
#ifdef DEBUG_QUEUE_NULL
  fprintf(_fp, "%s %d : NULL queues : ", __FILE__, __LINE__);
#endif
  for (list<C_task *>::iterator it = queue_null.begin(); 
       it != queue_null.end(); ++it) {
#ifdef DEBUG_QUEUE_NULL
    fprintf(_fp, "%s %d : %s : %ld\n",
	    __FILE__, __LINE__, (*it)->task_name, (*(*it)->ops_complexity));
#endif
    erase_task<T, U>(*it);
  }
  queue_null.clear();

  delete [] tasks_deallocLocalSchur_indcol;
  delete [] tasks_DTRSMScale_indcol;
  delete [] tasks_DTRSMScale_index;
  delete [] tasks_DTRSMScale_parents_index;
  delete [] tasks_DTRSMScale_ptr;
  delete [] tasks_DSymmGEMM_indcol;
  delete [] nops_queue[0];
  delete [] nops_queue;  
  
#if 0
  for (int level = level_last; level > 0; level--) {
    cout << "_tasks_Dsub: level = " << level << endl;
    const int begdom = 1U << (level - 1);
    const int enddom = begdom * 2;
    for (int d = 1; d < enddom; d++) {
      const int j = _btree->selfIndex(d);
      cout << "id = " << j
	   << " size = " << _tasks_Dsub[level][j].size() << endl;
      for (vector<C_task *>::const_iterator it = _tasks_Dsub[level][j].begin();
	   it != _tasks_Dsub[level][j].end(); ++it) {
	cout << (*it)->task_name << " :: ";
	list<C_Dsub_task *>*arg = (list<C_Dsub_task*>*)(*it)->func_arg;
	cout << arg->size() << " | ";
	int k = 0;
	for (list<C_Dsub_task *>::const_iterator it = arg->begin();
	     it != arg->end(); ++it, k++) {
	  cerr << k ;
	  cerr << "[ " << (*it)->atomic_id << " " << (*it)->atomic_size 
	       << (*it)->ir_bgn << " " << (*it)->ir_end << " ] ";
	}
      }
      cout << endl;
    } // loop : d
  }
  //  exit(-1);
#endif
  _queue_symb_allocated = true;
  _queue_numrc_allocated = true;
}

template 
void DissectionQueue<double, double>::
generate_queue(vector<DissectionMatrix<double>*>& dM,   
	       int nnz, double *coefs);

template 
void DissectionQueue<quadruple, quadruple>::
generate_queue(vector<DissectionMatrix<quadruple>*>& dM,   
	       int nnz, quadruple *coefs);

template 
void DissectionQueue<complex<double>, double>::
generate_queue(vector<DissectionMatrix<complex<double>, double>*>& dM,   
	       int nnz, complex<double> *coefs);

template 
void DissectionQueue<complex<quadruple>, quadruple>::
generate_queue(vector<DissectionMatrix<complex<quadruple>, quadruple>*>& dM,   
	       int nnz, complex<quadruple> *coefs);
//

template<typename T, typename U>
void DissectionQueue<T, U>::
generate_queue_fwbw(vector<DissectionMatrix<T, U>*>& dM,
		    int dim,
		    int nnz,
		    T *coefs)
{
  //  const T zero(0.0);
  //  const T none(-1.0);
  //  const T one(1.0);
  const int nb_level = _btree->NumberOfLevels();
  const int level_last = nb_level - 1;
  
  const int nb_doms = _btree->NumberOfSubdomains();

  // allocation of working array for sereial execution
  // find maximum size of submatrix among dense domains
  //
  _dim = dim;
  _nnz = nnz;
  _nrhs = new int*;
  _isTrans = new bool*;
  _x = new T*;
  _xi = new T**[nb_doms];
  _yi = new T**[nb_doms];
  _zi = new T**[nb_doms];
  _wi = new T**[nb_doms];

  for (int d = 1; d <= nb_doms; d++) {
    const int d0 = _btree->selfIndex(d);
    _xi[d0] = new T*;
    _yi[d0] = new T*;
    _wi[d0] = new T*;
    _zi[d0] = new T*;
  }

  _diag_contribs = new list<diag_contribution>[nb_doms];
  for (int level = level_last; level >= 0; level--) {
    const unsigned begdom = 1U << level;
    const unsigned enddom = begdom * 2;
    for (int d = begdom; d < enddom; d++) {
      int offset_src = 0;
      for (int ll = (level - 1); ll >=0; ll--) {
	const int father_id = _btree->nthfatherIndex(d, (level - ll));
	const int father_id0 = _btree->selfIndex(father_id);
	// copy data from Xavier's SetOfStrips to list<index_strip>
	list<index_strip> diag;
	const Dissection::SetOfStrips &diag_x = _btree->getFathersStrips(d)[ll];
	for (Dissection::SetOfStrips::const_iterator it = diag_x.begin();
	     it != diag_x.end(); ++it) {
	  diag.push_back(index_strip((*it).begin_dst,
				     (*it).begin_src + offset_src,
				     (*it).width));
	}
	list<diag_contribution> &tmp = _diag_contribs[father_id0];
	tmp.push_back(diag_contribution(d,
					_btree->sizeOfDomain(d),
					_btree->sizeOfFathersStrips(d),
					_btree->sizeOfDomain(father_id),
					diag));

	offset_src += diag_x.numberOfIndices();
      }
    }
  }
  vector<C_task*> tasks_a(1U << level_last);
  // sparse solver forward substitution phase
  {  // begin : sparse
    const unsigned begdom = 1U << level_last;
    const unsigned enddom = begdom * 2;
    int ipos = 0;
    for (int d = begdom; d < enddom; d++, ipos++) {  // to be parallelized 
      const int dd0 = d - begdom;
      const int d0 = _btree->selfIndex(d);
      const CSR_indirect &offdiag = _btree->getOffdiagCSR(d);
      int *indVals = offdiag.indVals;
      // get from the global array
      //      const int *loc2glob_diag = _btree->getDiagLoc2Glob(d);
      const int colors = dM[d0]->ColorTridiagBlockMatrix();
      C_SparseFw_arg<T, U> *arg = 
	new C_SparseFw_arg<T, U>(colors,
				 d,
				 _isSym,
				 dim,
				 _isTrans, // 
				 _nrhs,
				 dM[d0]->addrtridiagBlock(),
				 _x,
				 _yi[d0],
				 _zi[d0],
				 coefs,
				 (int)_btree->sizeOfDomain(d),
				 (int)_btree->sizeOfFathersStrips(d),
				 offdiag.ptRows,
				 offdiag.indCols,
				 (_isSym ? indVals : offdiag.indVals_unsym),
				 indVals,
				 _btree->getDiagLoc2Glob(d));
      string task_name = "a : " + to_string(d);
      *(arg->ops_complexity) = (long)dim;
      C_task *task = new C_task(C_SPARSESYMFW,
				task_name,
				(void *)arg,
				C_SparseFw<T, U>,
				1, 
				0,
				arg->ops_complexity);
      task->parallel_max = begdom;
      task->parallel_id = dd0;
      task->parents->clear(); // no dependency
      tasks_a[dd0] = task;
    } // loop : d
  }   // end : sparse
// updating RHS vector from level_last to all levels less than level_last
  vector<C_task*> tasks_i((1U << level_last) - 1);
  {
    const unsigned begdom = 1U;
    const unsigned enddom = 1U << level_last;
    for (int d = begdom; d < enddom; d++) {    // to be parallelized 
      const int dd0 = d - begdom;
      const int d0 = _btree->selfIndex(d);
      const int n_diag = _btree->sizeOfDomain(d);
      if (n_diag == 0) {
	string task_name = ("i dummy : " + to_string(d));
	C_dummy_arg *arg = new C_dummy_arg(_verbose, &_fp, d);
	tasks_i[dd0] = new C_task(C_DUMMY,
				  task_name,
				  (void *)arg,
				  C_dummy,
				  1, // atomic_size,
				  0, // atomic_id,
				  arg->ops_complexity);
      }
      else {
	C_Dsub_FwBw_arg<T> *arg = 
	  new C_Dsub_FwBw_arg<T>(dim,
				 _nrhs,
				 n_diag,
				 true,
				 level_last,
				 _btree,
				 &_diag_contribs[d0],
				 _x,
				 _yi[d0],
				 _zi,
				 _btree->getDiagLoc2Glob(d));
	string task_name = "i : " + to_string(d);
	*(arg->ops_complexity) = (long)n_diag;
	C_task *task = new C_task(C_DSUB_FWBW,
				  task_name,
				  (void *)arg,
				  C_Dsub_FwBw<T>,
				  1, 
				  0,
				  arg->ops_complexity);
	//task->func(task->func_arg);
	task->parallel_max = enddom;
	task->parallel_id = dd0;
	const int begdom0 = 1U << level_last;
	for (list<diag_contribution>::const_iterator it = _diag_contribs[d0].begin();
	     it != _diag_contribs[d0].end(); ++it) {
	  const int child_id = (*it).child_id;
	  if (_btree->nodeLayer(child_id) == level_last) {
	    task->parents->push_back(tasks_a[child_id - begdom0]); 
	  }
	}
	tasks_i[dd0] = task;
      }
    }         // loop : d
  }
// forward substitution of dense part
  vector<C_task*>** tasks_de = new vector<C_task*>*[level_last];
  tasks_de[0] = new vector<C_task*>[(1U << (level_last + 1)) - 1];
  vector<C_task*>** tasks_k = new vector<C_task*>*[level_last];
  tasks_k[0] = new vector<C_task*>[(1U << (level_last + 1)) - 1];
  vector<C_task*>* tasks_j = new vector<C_task *>[level_last];
  for (int level = (level_last - 1); level >= 0; level--) {
    const unsigned begdom = 1U << level;
    const unsigned enddom = begdom * 2;
    if (level > 0) {
      tasks_de[level] = tasks_de[0] + ((1U << level) - 1);
      tasks_k[level] = tasks_k[0] + ((1U << level) - 1);
    }
    for (int d = begdom; d < enddom; d++) {   // to be parallelized 
      const int dd0 = d - begdom;
      const int d0 = _btree->selfIndex(d);
      const int n_diag = _btree->sizeOfDomain(d);
      //      fprintf(stderr, "%s %d : %d %d\n", __FILE__, __LINE__, d, n_diag);
      const int n_offdiag = _btree->sizeOfFathersStrips(d); 
      SquareBlockMatrix<T> &Diag = dM[d0]->diagBlock();
      //      const int *permute = Diag.getPermute().getAddr();
      { // closure of tasks_de
	const int num_block = Diag.num_blocks();
	if (n_diag == 0) {
	  tasks_de[level][dd0].resize(1);
	  string task_name = ("d dummy : " + to_string(d));
	  C_dummy_arg *arg = new C_dummy_arg(_verbose, &_fp, d);
	  tasks_de[level][dd0][0] = new C_task(C_DUMMY,
					       task_name,
					       (void *)arg,
					       C_dummy,
					       1, // atomic_size,
					       0, // atomic_id,
					       arg->ops_complexity);
	}
	else {
	  tasks_de[level][dd0].resize((num_block * (num_block + 1)) / 2);
	  for (int k = 0; k < num_block; k++) {
	    const int kk = Diag.IndexBlock(k);
	    const int ncol = Diag.nrowBlock(k);
	    C_DenseFwBw_arg<T> *arg = 
	      new C_DenseFwBw_arg<T>(_isSym,
				     false,  // isBackward
				     (int)dim,
				     _isTrans,
				     _nrhs,
				     (int)n_diag,
				     (int)ncol, // 
				     (int)k,
				     (int)kk,
				     _xi[d0],
				     _wi[d0],
				     _yi[d0],
				     (T **)NULL,
				     dM[d0]->addrdiagBlock(),
				     (k == 0),  // isFirstBlock
				     (k == (num_block - 1)), // isLastBlock
				     (int *)NULL,
				     _verbose,
				     &_fp);
	    string task_name = "d : " + to_string(d) + " : " + to_string(k);
	    *(arg->ops_complexity) = (long)n_diag;
	    C_task *task = new C_task(C_DENSE_SYMFW_DIAG,
				      task_name,
				      (void *)arg,
				      C_DenseFwBw_diag<T>,
				      1, 
				      0,
				      arg->ops_complexity);
	    if (k == 0) {
	      task->parallel_max = 1;
	      task->parallel_id = 0;
	      tasks_de[level][dd0][0] = task;
	      if (level == (level_last - 1)) {
		task->parents->push_back(tasks_i[d0]); //
	      }
	      else {
		task->parents->push_back(tasks_j[level + 1][d0]);
	      }
	    }
	    else {
	      task->parallel_max = num_block - k;
	      task->parallel_id = 0;
	      task->atomic_size = 2;
	      task->atomic_id = 1;
	      int ipos;
	      // depending on e_{k k-1} 
	      ipos = 2 + num_block * (k - 1) - ((k - 1) * k) / 2 + k - 1;
	      task->parents->push_back(tasks_de[level][dd0][ipos - 1]); //
	      tasks_de[level][dd0][ipos] = task;
	      for (int m = 0; m < (k - 1); m++) {
		ipos = num_block * m + k + 1 - (m * (m + 1)) / 2;
		task->parents->push_back(tasks_de[level][dd0][ipos]); //
	      }
	    }
// updating of lower blocks xi[d0] + ii are independent among i
	    for (int i = (k + 1); i < num_block; i++) {
	      const int ii = Diag.IndexBlock(i);
	      const int nrow = Diag.nrowBlock(i);
	      C_DenseFwBwOffdiag_arg<T> *arg = 
		new C_DenseFwBwOffdiag_arg<T>(_isSym ? true : false,  //trans
					      true,  //isLower
					      (int)dim,
					      _isTrans,
					      _nrhs,
					      (int)n_diag,
					      (int)n_diag,
					      (int)n_diag, //ldc
					      (int)nrow,
					      (int)ncol,
					      _wi[d0],
					      kk, // ii
					      _yi[d0], //
					      _xi[d0], //
					      ii, // jj
					      dM[d0]->addrdiagBlock(),
					      i,  // i_block
					      k,  // j_block
					      _none, // alpha,
					      _one); // beta
	      string task_name = ("e : " + to_string(d) + " : " + to_string(k)
				  + " " + to_string(i));
	      *(arg->ops_complexity) = (long)dim;
	      C_task *task = new C_task(C_DENSE_SYMFW_OFFDIAG,
					task_name,
					(void *)arg,
					C_DenseFwBw_offdiag<T>,
					1, 
					0,
					arg->ops_complexity);
	      task->parallel_max = num_block - (k + 1);
	      task->parallel_id = i - (k + 1);
	      int ipos;
	      // depending on d_{k-1}
	      if (k > 0) {
		ipos = num_block * (k - 1) + i + 1 - ((k - 1) * k) / 2;
		task->parents->push_back(tasks_de[level][dd0][ipos]);
	      }
	      if (k > 0) {
		ipos = 2 + num_block * (k - 1) - ((k - 1) * k) / 2 + k - 1;
	      }
	      else {
		ipos = 0;
	      }
	      task->parents->push_back(tasks_de[level][dd0][ipos]);
	      if (i == (k + 1)) {
		task->atomic_size = 2; // 
		task->atomic_id = 0;   //
		ipos = 2 + num_block * k - (k * (k + 1)) / 2 + k - 1;
		tasks_de[level][dd0][ipos] = task;
	      }
	      else { // i > (k + 1)
		ipos = num_block * k + i + 1 - (k * (k + 1)) / 2;	      
		tasks_de[level][dd0][ipos] = task;
	      }  // if (i == (k + 1))
	    } // loop : i
	  } // loop : k
	} // if (n_diag == 0) 
      }   // closure of tasks_de
      if (level > 0) {
	//	const int num_block = Diag.num_blocks(); //
	const int num_block_r = Diag.num_blocks(); 

	RectBlockMatrix<T> *upperblock = 
	    (_isSym ? dM[d0]->addrupperBlock() : dM[d0]->addrlowerBlock());
	RectBlockMatrix<T> *lowerblock = dM[d0]->addrupperBlock();
	const int num_block_c = upperblock->num_blocks_c();
	if (n_diag == 0) {
	  tasks_k[level][dd0].resize(1);
	  string task_name = ("k dummy : " + to_string(d));
	  C_dummy_arg *arg = new C_dummy_arg(_verbose, &_fp, d);
	  tasks_k[level][dd0][0] = new C_task(C_DUMMY,
					      task_name,
					      (void *)arg,
					      C_dummy,
					      1, // atomic_size,
					      0, // atomic_id,
					      arg->ops_complexity);
	}
	else {
	  tasks_k[level][dd0].resize(num_block_r * num_block_c);

	//	const double alpha = 1.0;
	  for (int i = 0; i < num_block_r; i++) {
	    const int ii = Diag.IndexBlock(i);
	    const int nrow = Diag.nrowBlock(i);
	    const T beta = (i == 0) ? _zero : _one;
	    
	    for (int j = 0; j < num_block_c; j++) {
	      const int jj = upperblock->IndexBlock_c(j); // j * SIZE_B1;
	      const int ncol = upperblock->ncolBlock(j);
	      C_StripsFwBwOffdiag_arg<T> *arg = 
		new C_StripsFwBwOffdiag_arg<T>(true,  // isLower
					       (int)dim,
					       _isTrans,
					       _nrhs,
					       (int)n_diag,
					       (int)n_diag,
					       (int)n_offdiag,
					       (int)ncol, // refering to upper
					       (int)nrow,
					       _wi[d0],
					       ii, // ii,
					       _zi[d0],
					       jj, // jj,
					       upperblock,
					       lowerblock, 
					       i, 
					       j,
					       _one, // alpha,
					       beta);
	      string task_name = ("k : " + to_string(d) + " : " + to_string(i)
				  +  " " + to_string(j));
	      *(arg->ops_complexity) = (long)n_diag;
	      C_task *task = new C_task(C_STRIPS_SYMFW_OFFDIAG,
					task_name,
					(void *)arg,
				    C_StripsFwBw_offdiag<T>,
					1, 
					0,
					arg->ops_complexity);
	      task->parallel_max = num_block_c;  // sequential among i
	      task->parallel_id = j;             // parallel among j
	      int ipos;
	      if (i == 0) {
		task->parents->push_back(tasks_de[level][dd0][0]);
	      }
	      else {
		ipos = 2 + num_block_r * (i - 1) - ((i - 1) * i) / 2 + i - 1;
		task->parents->push_back(tasks_de[level][dd0][ipos]);
		ipos = (i - 1) * num_block_c + j;
		task->parents->push_back(tasks_k[level][dd0][ipos]);
	      }
	      ipos = i * num_block_c + j;
	      tasks_k[level][dd0][ipos] = task;
	    } // loop : j
	  }   // loop : i
	}
      }     // if (level > 0)
      // solving D y = x
    } // loop : d
// updating RHS vector from level to all levels less than level
    {
      const unsigned begdom = 1U;
      const unsigned enddom = 1U << level;
      tasks_j[level].resize(enddom - 1);
      for (int d = begdom; d < enddom; d++) {   // to be parallelized 
	const int dd0 = d - begdom;
	const int d0 = _btree->selfIndex(d);
	const int n_diag = _btree->sizeOfDomain(d);
	if (n_diag == 0) {
	  string task_name = ("j dummy : " + to_string(d));
	  C_dummy_arg *arg = new C_dummy_arg(_verbose, &_fp, d);
	  tasks_j[level][dd0] = new C_task(C_DUMMY,
					       task_name,
					       (void *)arg,
					       C_dummy,
					       1, // atomic_size,
					       0, // atomic_id,
					       arg->ops_complexity);
	}
	else {
	  C_Dsub_FwBw_arg<T> *arg = 
	    new C_Dsub_FwBw_arg<T>(dim,
				   _nrhs,
				   -1, // n_diag
				   false,                // access_global
				   level,
				   _btree,
				   &_diag_contribs[d0],
				   (T **)NULL,      // x
				   _yi[d0],
				   _zi,
				   _btree->getDiagLoc2Glob(d));
	  string task_name = "j : " + to_string(level) + " : " + to_string(d);
	  *(arg->ops_complexity) = (long)dim;
	  C_task *task = new C_task(C_DSUB_FWBW,
				    task_name,
				    (void *)arg,
				    C_Dsub_FwBw<T>,
				    1, 
				    0,
				    arg->ops_complexity);
	  //      task->func(task->func_arg);
	  task->parallel_max = enddom;
	  task->parallel_id = d0;
	  task->parents->push_back(tasks_i[dd0]);
	  for (list<diag_contribution>::const_iterator it = _diag_contribs[d0].begin();
	       it != _diag_contribs[d0].end(); ++it) {
	    if ((*it).diag_strip.size() > 0) {      // to avoid null contribution
	      const int child_id = (*it).child_id;
	      if (_btree->nodeLayer(child_id) == level) {
		const int child_id0 = _btree->selfIndex(child_id);
		RectBlockMatrix<T> *upperblock = dM[child_id0]->addrupperBlock();
		const int num_block_r = upperblock->num_blocks_r();
		const int num_block_c = upperblock->num_blocks_c();
		for (int j = 0; j < num_block_c; j++) {
		  int ipos = (num_block_r - 1) * num_block_c + j;
		  task->parents->push_back(tasks_k[0][child_id0][ipos]);
		}
	      } // if (_btree->nodeLayer(child_id) == level)
	    }  // (if (*it).diag_strip.size() > 0) 
	  }         // loop : it
	  tasks_j[level][dd0] = task;
	}
      }  // loop : d
    }
  }   // loop : level
  vector<C_task*>* tasks_h = new vector<C_task *>[level_last];
  vector<C_task*>** tasks_l = new vector<C_task*>*[level_last];
  vector<C_task*>** tasks_fg = new vector<C_task*>*[level_last];
  tasks_l[0] = new vector<C_task*>[(1U << (level_last + 1)) - 1];
  tasks_fg[0] = new vector<C_task*>[(1U << (level_last + 1)) - 1];
  for (int level = 0; level < level_last; level++) {
    const unsigned begdom = 1U << level;
    const unsigned enddom = begdom * 2;
    tasks_h[level].resize(enddom - 1);
    if (level > 0) {
      tasks_fg[level] = tasks_fg[0] + ((1U << level) - 1);
      tasks_l[level] = tasks_l[0] + ((1U << level) - 1);
    }
    for (int d = begdom; d < enddom; d++) {    // to be parallelized 
      const int dd0 = d - begdom;
      const int d0 = _btree->selfIndex(d);
      const int n_diag = _btree->sizeOfDomain(d);
      const int n_offdiag = _btree->sizeOfFathersStrips(d); 
      SquareBlockMatrix<T> &Diag = dM[d0]->diagBlock();
      if (level > 0) {
	if (n_diag == 0) {
	  string task_name = ("h dummy : " + to_string(d));
	  C_dummy_arg *arg = new C_dummy_arg(_verbose, &_fp, d);
	  tasks_h[level][dd0] = new C_task(C_DUMMY,
					       task_name,
					       (void *)arg,
					       C_dummy,
					       1, // atomic_size,
					       0, // atomic_id,
					       arg->ops_complexity);
	}
	else {
	  C_Dfill_FwBw_arg<T> *arg = 
	    new C_Dfill_FwBw_arg<T>(_nrhs,
				    d,
				    level,
				    _btree,
				    (int)n_offdiag,
				    _yi,
				    _zi[d0]);
	  string task_name = "h : " + to_string(d);
	  *(arg->ops_complexity) = (long)n_offdiag;
	  C_task *task = new C_task(C_DENSE_SYMFILL,
				    task_name,
				    (void *)arg,
				    C_Dfill_FwBw<T>,
				    1, 
				    0,
				    arg->ops_complexity);
	  task->parallel_max = enddom - 1;
	  task->parallel_id = dd0;
	  task->parents->clear();
	  for (int ll = (level - 1); ll >= 0; ll--) {
	    const int father_id = _btree->nthfatherIndex(d, (level - ll));
	    const int father_id0 = _btree->selfIndex(father_id);
	    const Dissection::SetOfStrips &diag = _btree->getFathersStrips(d)[ll];
	    if (diag.numberOfIndices() > 0) {
	      const int ipos = tasks_fg[0][father_id0].size() - 1; // the last
	      C_task *task_parent = tasks_fg[0][father_id0][ipos];
	      bool flag = false;
	      for (list<C_task *>::const_iterator it = task->parents->begin();
		   it != task->parents->end(); ++it) {
		if ((*it) == task_parent) {
		  flag = true;
		  break;
		}
	      } // loop : it
	      if (flag == false) {
		task->parents->push_back(task_parent);
	      }
	    }
	  }
	  tasks_h[level][dd0] = task;
	}
	const int num_block_r = Diag.num_blocks();
	if (n_diag == 0) {
	  tasks_l[level][dd0].resize(1);
	  string task_name = ("l dummy : " + to_string(d));
	  C_dummy_arg *arg = new C_dummy_arg(_verbose, &_fp, d);
	  tasks_l[level][dd0][0] = new C_task(C_DUMMY,
					       task_name,
					       (void *)arg,
					       C_dummy,
					       1, // atomic_size,
					       0, // atomic_id,
					       arg->ops_complexity);
	}
	else {
	  tasks_l[level][dd0].resize(num_block_r);
	  for (int i = 0; i < num_block_r; i++) {
	    const int ii = Diag.IndexBlock(i);
	    const int nrow = Diag.nrowBlock(i);
	    C_StripsFwBwOffdiag_arg<T> *arg = 
	      new C_StripsFwBwOffdiag_arg<T>(false, // isLower
					     (int)dim,
					     _isTrans,
					     _nrhs,
					     (int)n_diag,
					     (int)n_offdiag,
					     (int)n_diag,
					     (int)nrow,
					     (int)n_offdiag,
					     _zi[d0],
					     0, // ii,
					     _xi[d0],
					     ii, // jj,
					     dM[d0]->addrupperBlock(),
					     dM[d0]->addrlowerBlock(),
					     0, 
					     i,
					     _none, //alpha,
					     _one); // beta);
	    string task_name = "l : " + to_string(d) + " : " + to_string(i);
	    *(arg->ops_complexity) = (long)n_diag;
	    C_task *task = new C_task(C_STRIPS_SYMFW_OFFDIAG,
				      task_name,
				      (void *)arg,
				      C_StripsFwBw_offdiag<T>,
				      1, 
				      0,
				      arg->ops_complexity);
	    int ipos;
	    task->parents->push_back(tasks_h[level][dd0]);
	    const int father_id = _btree->fatherIndex(d);
	    const int father_id0 = _btree->selfIndex(father_id);
	    ipos = tasks_fg[0][father_id0].size() - 1; // the last of tasks_fg
	    if (ipos >= 0) {
	      task->parents->push_back(tasks_fg[0][father_id0][ipos]);
	    }
	    tasks_l[level][dd0][i] = task;
	  } // loop : i
	} // if (n_diag == 0)
      } // if (level > 0)
      {
	SquareBlockMatrix<T> &Diag = dM[d0]->diagBlock();
	const int num_block = Diag.num_blocks();
	if (n_diag == 0) {
	  tasks_fg[level][dd0].resize(1);
	  string task_name = ("d dummy : " + to_string(d));
	  C_dummy_arg *arg = new C_dummy_arg(_verbose, &_fp, d);
	  tasks_fg[level][dd0][0] = new C_task(C_DUMMY,
					       task_name,
					       (void *)arg,
					       C_dummy,
					       1, // atomic_size,
					       0, // atomic_id,
					       arg->ops_complexity);

	}
	else {
	  tasks_fg[level][dd0].resize((num_block * (num_block + 1)) / 2);
	  int k0 = 0;
	  for (int k = (num_block - 1); k >= 0; k--, k0++) {
	    const int kk = Diag.IndexBlock(k);
	    const int ncol = Diag.nrowBlock(k);
	    C_DenseFwBw_arg<T> *arg = 
	      new C_DenseFwBw_arg<T>(_isSym, 
				     true,  // isBackward
				     (int)dim,
				     _isTrans,
				     _nrhs,
				     (int)n_diag,
				     (int)ncol, // 
				     (int)k,
				     (int)kk,
				     _xi[d0],
				     (T **)NULL,
				     _yi[d0], //
				     (k == 0) ? _x : (T **)NULL,
				     dM[d0]->addrdiagBlock(),
				     (k == 0), // isFirstBlock
				     (k == (num_block - 1)), // isLastBlock
				     (k == 0) ? _btree->getDiagLoc2Glob(d) : (int *)NULL,
				     _verbose,
				     &_fp);
	    string task_name = "f : "  + to_string(d) + " : " + to_string(k);
	    *(arg->ops_complexity) = (long)dim;
	    C_task *task = new C_task(C_DENSE_SYMFW_DIAG,
				      task_name,
				      (void *)arg,
				      C_DenseFwBw_diag<T>,
				      1, 
				      0,
				      arg->ops_complexity);
	    if (k0 == 0) {
	      task->parallel_max = 1;
	      task->parallel_id = 0;
	      tasks_fg[level][dd0][0] = task;
	      if (level == 0){ 
		const int ipos = tasks_de[0][0].size() - 1; //the last of the fw
		if (ipos >= 0) {
		  task->parents->push_back(tasks_de[0][0][ipos]);
		}
	      }
	      else {
		task->parents->push_back(tasks_l[0][d0][k]);
	      }
	    }
	    else {
	      task->parallel_max = num_block - k0;
	      task->parallel_id = 0;
	      task->atomic_size = 2;
	      task->atomic_id = 1;
	      int ipos;
	      ipos = 2 + num_block * (k0 - 1) - ((k0 - 1) * k0) / 2 + k0 - 1;
	      task->parents->push_back(tasks_fg[level][dd0][ipos - 1]); //
	      tasks_fg[level][dd0][ipos] = task;
	      for (int m = 0; m < (k0 - 1); m++) {
		ipos = num_block * m + k0 + 1 - (m * (m + 1)) / 2;	 
		task->parents->push_back(tasks_fg[level][dd0][ipos]); //
	      }
	    }
	    int i0 = k0 + 1; // a trick to use the same formula to count ipos
	    for (int i = (k - 1); i >= 0; i--, i0++) {
	      const int ii = Diag.IndexBlock(i);
	      const int nrow = Diag.nrowBlock(i); // SIZE_B1; // ?
	      C_DenseFwBwOffdiag_arg<T> *arg = 
		new C_DenseFwBwOffdiag_arg<T>(false, // trans
					      false, // isLower
					      (int)dim,
					      _isTrans,
					      _nrhs,
					      (int)n_diag,
					      (int)n_diag,
					      (int)n_diag,
					      (int)nrow, 
					      (int)ncol, 
					      _xi[d0],
					      kk, // ii
					      _xi[d0],
					      _xi[d0], // (double **)NULL,
					      ii, // jj
					      dM[d0]->addrdiagBlock(),
					      i,  // i_block
					      k,  // j_block
					      _none, // alpha,
					      _one); //beta);
	      string task_name = ("g : " + to_string(d) + " : " + to_string(k)
				  + " " + to_string(i));
	      *(arg->ops_complexity) = (long)dim;
	      C_task *task = new C_task(C_DENSE_SYMFW_OFFDIAG,
					task_name,
					(void *)arg,
					C_DenseFwBw_offdiag<T>,
					1, 
					0,
					arg->ops_complexity);
	      if (level > 0) {
		task->parents->push_back(tasks_l[0][d0][i]);
	      }
	      task->parallel_max = num_block - (k0 + 1);
	      int ipos;
	      if (k0 > 0) {
		ipos = num_block * (k0 - 1) + i0 + 1 - ((k0 - 1) * k0) / 2; 
		task->parents->push_back(tasks_fg[level][dd0][ipos]);
	      }
	      if (k0 > 0) {
		ipos = 2 + num_block * (k0 - 1) - ((k0 - 1) * k0) / 2 + k0 - 1;
	      }
	      else {
		ipos = 0;
	      }
	      task->parents->push_back(tasks_fg[level][dd0][ipos]);
	      if (i0 == (k0 + 1)) {
		task->parallel_id = 0;
		task->atomic_size = 2; // 
		task->atomic_id = 0;   //
		ipos = 2 + num_block * k0 - (k0 * (k0 + 1)) / 2 + k0 - 1;
	      }
	      else { // i0 > (k0 + 1)
		task->parallel_id = i - (k0 + 1);
		ipos = num_block * k0 + i0 + 1 - (k0 * (k0 + 1)) / 2;	      
	      }
	      tasks_fg[level][dd0][ipos] = task;
	    } // loop : i
	  } // loop : k
	} // if (n_diag == 0)
      }
    } // loop : d
  }   // loop : level
  vector<C_task*> tasks_b(1U << level_last);
  {   // begin : sparse
    const unsigned begdom = 1U << level_last;
    const unsigned enddom = begdom * 2;
    for (int d = begdom; d < enddom; d++) {     // to be parallelized 
      const int dd0 = d - begdom;
      const int d0 = _btree->selfIndex(d);
      const CSR_indirect &offdiag = _btree->getOffdiagCSR(d);
      const int colors = dM[d0]->ColorTridiagBlockMatrix();
      C_SparseBw_arg<T, U> *arg = 
	new C_SparseBw_arg<T, U>(colors,
				 d,
				 _isSym,
				 (int)dim,
				 _isTrans,
				 _nrhs,
				 _btree,
				 (int)level_last,
				 dM[d0]->addrtridiagBlock(),
				 _x,
				 _yi,
				 _xi[d0],
				 _yi[d0],
				 _zi[d0],
				 coefs, 
				 offdiag.ptRows,
				 offdiag.indCols,
				 offdiag.indVals,
				 offdiag.indVals_unsym);
      string task_name = "b : " + to_string(d);
      *(arg->ops_complexity) = (long)dim;
      C_task *task = new C_task(C_SPARSESYMBW,
				task_name,
				(void *)arg,
				C_SparseBw<T, U>,
				1, 
				0,
				arg->ops_complexity);
      task->parallel_max = begdom;
      task->parallel_id = dd0;
      for (int ll = (level_last - 1); ll >= 0; ll--) {
	const int father_id = _btree->nthfatherIndex(d, (level_last - ll));
	const int father_id0 = _btree->selfIndex(father_id);
	const Dissection::SetOfStrips &diag = _btree->getFathersStrips(d)[ll];
	const int ipos = tasks_fg[0][father_id0].size() - 1;
	if ((diag.numberOfIndices() > 0) && (ipos >= 0)) {
	  task->parents->push_back(tasks_fg[0][father_id0][ipos]);
	}
      }
      tasks_b[dd0] = task;
    } // loop : d
  }   // end : sparse

  int count_tasks = 0;
  for (int i = 0; i < (1U << level_last); i++) {
    count_tasks++;
  }
  for (int i = 0; i < ((1U << level_last) - 1); i++) {
    count_tasks++;
  }
  for (int level = (level_last - 1); level >= 0; level--) {
    const int begdom = (int)(1U << level);
    const int enddom = begdom * 2;
    for (int d = begdom; d < enddom; d++) {
      const int dd0 = d - begdom;
      for (int i = 0; i < tasks_de[level][dd0].size(); i++) {
	count_tasks++;
      }
      if (level > 0) {
	for (int i = 0; i < tasks_k[level][dd0].size(); i++) {
	  count_tasks++;
	}
      }
    } // loop : d 
   for (int d = 1; d < begdom; d++) {   // to be parallelized 
	  count_tasks++;	  
    }
  }
  for (int level = 0; level < level_last; level++) {
    const int begdom = (int)1U << level;
    const int enddom = begdom * 2;
    for (int d = begdom; d < enddom; d++) {
      const int dd0 = d - begdom;
      if (level > 0) {
	count_tasks++;
	for (int i = 0; i < tasks_l[level][dd0].size(); i++) {
	  count_tasks++;
	}
      }
      for (int i = 0; i < tasks_fg[level][dd0].size(); i++) {
	count_tasks++;
      }
    } // loop : d
  }
  for (int i = 0; i < (1U << level_last); i++) {
    count_tasks++;
  }

  vector<C_task *> *tasks_tmp = new vector<C_task *>;
  tasks_tmp->resize(count_tasks);
  count_tasks = 0;
  for (int i = 0; i < (1U << level_last); i++) {
    (*tasks_tmp)[count_tasks++] = tasks_a[i];
  }
  for (int i = 0; i < ((1U << level_last) - 1); i++) {
    (*tasks_tmp)[count_tasks++] = tasks_i[i];
  }
  for (int level = (level_last - 1); level >= 0; level--) {
    const int begdom = (int)(1U << level);
    const int enddom = begdom * 2;
    for (int d = begdom; d < enddom; d++) {
      const int dd0 = d - begdom;
      for (int i = 0; i < tasks_de[level][dd0].size(); i++) {
	(*tasks_tmp)[count_tasks++] = tasks_de[level][dd0][i];
      }
      if (level > 0) {
	for (int i = 0; i < tasks_k[level][dd0].size(); i++) {
	  (*tasks_tmp)[count_tasks++] = tasks_k[level][dd0][i];
	}
      }
    } // loop : d 
   for (int d = 1; d < begdom; d++) {   // to be parallelized 
      const int dd0 = d - 1;
      (*tasks_tmp)[count_tasks++] = tasks_j[level][dd0];
    }
  }
  for (int level = 0; level < level_last; level++) {
    const int begdom = (int)(1U << level);
    const int enddom = begdom * 2;
    for (int d = begdom; d < enddom; d++) {
      const int dd0 = d - begdom;
      if (level > 0) {
	(*tasks_tmp)[count_tasks++] = tasks_h[level][dd0];
	for (int i = 0; i < tasks_l[level][dd0].size(); i++) {
	  (*tasks_tmp)[count_tasks++] = tasks_l[level][dd0][i];
	}
      }
      for (int i = 0; i < tasks_fg[level][dd0].size(); i++) {
	(*tasks_tmp)[count_tasks++] = tasks_fg[level][dd0][i];
      }
    } // loop : d
  }
  for (int i = 0; i < (1U << level_last); i++) {
    (*tasks_tmp)[count_tasks++] = tasks_b[i];
  }
  string task_name = "fwbw whole"; 
  _queue_fwbw =  new C_task_seq(C_FWBW,
				task_name, // dummy
				(-1), // mutex_id
				TASK_SINGLE, //TASK_SINGLE,
				_num_threads,
				-1, // level
				-1, // phase
				tasks_tmp,
				0,
				//itmp,
				count_tasks,
				0); // dummy

  _dissectionRuntime->set_queue_fwbw(_queue_fwbw);
  if (_verbose) {
    fprintf(_fp, 
	    "%s %d : void DissectionQueue::generate_queue_fwbw\n",
	  __FILE__, __LINE__);
  }
#if 0
  fprintf(_fp, "-- queue begin --\n");
  for (vector<C_task*>::const_iterator it = _queue_fwbw->queue->begin(); 
       it != _queue_fwbw->queue->end(); ++it) {
    fprintf(_fp, "%s : %d :: ", (*it)->task_name, (*it)->parents->size());
    for (list<C_task*>::const_iterator jt = (*it)->parents->begin();
	 jt != (*it)->parents->end(); ++jt) {
      fprintf(_fp, "%s : ", (*jt)->task_name);
    }
    fprintf(_fp, "\n");
  }
  fprintf(_fp, "-- queue end --\n");
#endif

  delete [] tasks_de[0]; // 
  delete [] tasks_de;    // vector<C_task*>** 
  delete [] tasks_k[0];  //
  delete [] tasks_k;     // vector<C_task*>** 
  delete [] tasks_j;     // vector<C_task*>*
  delete [] tasks_h;     // vector<C_task*>*
  delete [] tasks_l[0];  //
  delete [] tasks_l;     // vector<C_task*>**
  delete [] tasks_fg[0];
  delete [] tasks_fg;    // vector<C_task*>**
  // vector<C_task*> tasks_a, tasks_b, and tasks_i are deallocated automatically
  _queue_fwbw_allocated = true;
}

template
void DissectionQueue<double, double>::
generate_queue_fwbw(vector<DissectionMatrix<double>*>& dissectionMatrix,
		    int dim, int nnz, double *coefs);

template
void DissectionQueue<quadruple, quadruple>::
generate_queue_fwbw(vector<DissectionMatrix<quadruple>*>& dissectionMatrix,
		    int dim, int nnz, quadruple *coefs);

template
void DissectionQueue<complex<double>, double>::
generate_queue_fwbw(vector<DissectionMatrix<complex<double>, double>*>& dissectionMatrix,
		    int dim, int nnz, complex<double> *coefs);

template
void DissectionQueue<complex<quadruple>, quadruple>::
generate_queue_fwbw(vector<DissectionMatrix<complex<quadruple>, quadruple>*>& dissectionMatrix,
		    int dim, int nnz, complex<quadruple> *coefs);
//

template<typename T, typename U>
void DissectionQueue<T, U>::exec_symb_fact(void)
{
  _dissectionRuntime->exec_symb_fact();
}

template
void DissectionQueue<double, double>::exec_symb_fact();

template
void DissectionQueue<quadruple, quadruple>::exec_symb_fact();

template
void DissectionQueue<complex<double>, double>::exec_symb_fact();

template
void DissectionQueue<complex<quadruple>, quadruple>::exec_symb_fact();

//

template<typename T, typename U>
void DissectionQueue<T, U>::exec_num_fact(const int called,
					  const double eps_piv, 
					  const bool kernel_detection,
					  const int aug_dim,
					  const U eps_machine)
{
  _eps_piv = eps_piv;
  _kernel_detection = kernel_detection;
  _aug_dim = aug_dim;
  _eps_machine = eps_machine;
  _dissectionRuntime->exec_num_fact(called);
}

template
void DissectionQueue<double, double>::
exec_num_fact(const int called,
	      const double eps_piv,
	      const bool kernel_detection,
	      const int aug_dim,
	      const double eps_machine);
template
void DissectionQueue<quadruple, quadruple>::
exec_num_fact(const int called,
	      const double eps_piv,
	      const bool kernel_detection,
	      const int aug_dim,
	      const quadruple eps_machine);
template
void DissectionQueue<complex<double>, double>::
exec_num_fact(const int called,
	      const double eps_piv, 
	      const bool kernel_detection,
	      const int aug_dim,
	      const double eps_machine);

template
void DissectionQueue<complex<quadruple>, quadruple>::
exec_num_fact(const int called,
	      const double eps_piv, 
	      const bool kernel_detection,
	      const int aug_dim,
	      const quadruple eps_machine);
//

template<typename T, typename U>
void DissectionQueue<T, U>::exec_fwbw_seq(T *x, const int nrhs, bool isTrans)
{
  const int nb_doms = _btree->NumberOfSubdomains();
  int ntmp;

  *_x = x;
  *_nrhs = (int *)&nrhs;
  *_isTrans = &isTrans;

  ntmp = 0;
  for (int d = 1; d <= nb_doms; d++) {
    ntmp += _btree->sizeOfDomain(d);
  }
  ColumnMatrix<T> xxi(ntmp, nrhs);
  ColumnMatrix<T> yyi(ntmp, nrhs);
  ColumnMatrix<T> wwi(ntmp, nrhs);
  ntmp = 0;
  for (int d = 1; d <= nb_doms; d++) {
    const int d0 = _btree->selfIndex(d);
    *_xi[d0] = xxi.addrCoefs() + ntmp * nrhs;
    *_yi[d0] = yyi.addrCoefs() + ntmp * nrhs;
    *_wi[d0] = wwi.addrCoefs() + ntmp * nrhs;
    ntmp += _btree->sizeOfDomain(d);
  }
  ntmp = 0;
  for (int d = 1; d <= nb_doms; d++) {
    ntmp += _btree->sizeOfFathersStrips(d);
  }
  ColumnMatrix<T> zzi(ntmp, nrhs);
  ntmp = 0;
  for (int d = 1; d <= nb_doms; d++) {
    const int d0 = _btree->selfIndex(d);
    *_zi[d0] = zzi.addrCoefs() + ntmp * nrhs;
    ntmp += _btree->sizeOfFathersStrips(d);
  }
  if(_verbose) {
    fprintf(_fp, "exec_fwbw_seq for nrhs = %d \n", nrhs);
  }
  _dissectionRuntime->exec_fwbw_seq();

  xxi.free();
  yyi.free();
  wwi.free();
  zzi.free();
}

template 
void DissectionQueue<double, double>::exec_fwbw_seq(double *x, const int nrhs, 
						    bool isTrans);

template 
void DissectionQueue<quadruple, quadruple>::exec_fwbw_seq(quadruple *x,
							  const int nrhs, 
							  bool isTrans);
template 
void DissectionQueue<complex<double>, double>::
exec_fwbw_seq(complex<double> *x, 
	      const int nrhs, 
	      bool isTrans);

template 
void DissectionQueue<complex<quadruple>, quadruple>::
exec_fwbw_seq(complex<quadruple> *x, 
	      const int nrhs, 
	      bool isTrans);

//

template<typename T, typename U>
void DissectionQueue<T, U>::exec_fwbw(T *x, const int nrhs, bool isTrans)
{
  const int num_threads = _num_threads; //_num_threads_symb;

  //  struct timespec ts0, ts1;
  const int nb_doms = _btree->NumberOfSubdomains();
  int ntmp;
  if (_verbose) {
    fprintf(_fp, "fwbw for nrhs = %d with %d threads\n", nrhs, num_threads);
  }
  *_x = x;
  *_nrhs = (int *)&nrhs;
  *_isTrans = &isTrans;

  ntmp = 0;
  for (int d = 1; d <= nb_doms; d++) {
    ntmp += _btree->sizeOfDomain(d);
  }
  ColumnMatrix<T> xxi(ntmp, nrhs);
  ColumnMatrix<T> yyi(ntmp, nrhs);
  ColumnMatrix<T> wwi(ntmp, nrhs);
  ntmp = 0;
  for (int d = 1; d <= nb_doms; d++) {
    const int d0 = _btree->selfIndex(d);
    *_xi[d0] = xxi.addrCoefs() + ntmp * nrhs;
    *_yi[d0] = yyi.addrCoefs() + ntmp * nrhs;
    *_wi[d0] = wwi.addrCoefs() + ntmp * nrhs;
    ntmp += _btree->sizeOfDomain(d);
  }
  ntmp = 0;
  for (int d = 1; d <= nb_doms; d++) {
    ntmp += _btree->sizeOfFathersStrips(d);
  }
  ColumnMatrix<T> zzi(ntmp, nrhs);
  ntmp = 0;
  for (int d = 1; d <= nb_doms; d++) {
    const int d0 = _btree->selfIndex(d);
    //    *_zi[d0] = &zzi[ntmp * nrhs];
    *_zi[d0] = zzi.addrCoefs()+ ntmp * nrhs;
    ntmp += _btree->sizeOfFathersStrips(d);
  }

  _dissectionRuntime->exec_fwbw();

  xxi.free();
  yyi.free();
  wwi.free();
  zzi.free();
}

template
void DissectionQueue<double, double>::exec_fwbw(double *x, const int nrhs,
						bool isTrans); 
template
void DissectionQueue<quadruple, quadruple>::exec_fwbw(quadruple *x,
						      const int nrhs,
						      bool isTrans); 

template
void DissectionQueue<complex<double>, double>::exec_fwbw(complex<double> *x, 
							 const int nrhs,
							 bool isTrans);

template
void DissectionQueue<complex<quadruple>, quadruple>::
exec_fwbw(complex<quadruple> *x, 
	   const int nrhs, bool isTrans);
//
// #define DEBUG_ERASE

template<typename T, typename U>
void DissectionQueue<T, U>::erase_queue(void)
{
  const int num_threads = _num_threads;
  if (_verbose) {
    fprintf(_fp, "%s %d : void QueueRuntime::erase_queue(void)",
	    __FILE__, __LINE__);
  }
  if (_queue_symb_allocated) {
#ifdef DEBUG_ERASE
  cerr << "symbolic ";
#endif
 vector<C_task *> &queue = *(_queue_symb->queue);
  for (int j = _queue_symb->begin; j < _queue_symb->end; j++) {
    C_task *task = queue[j];
#ifdef DEBUG_ERASE
    cerr << task->task_name << " ";
#endif
    //    delete [] task->task_name; // char *task_name_cstr = new char []
    erase_task<T, U>(task);
  }  // loop : j
#ifdef DEBUG_ERASE
  cerr << _queue_symb->task_name << " " << endl;
#endif
  //  delete [] _queue_symb->task_name;
  delete _queue_symb->queue;
  delete _queue_symb;
  }

#ifdef DEBUG_ERASE
  cerr << "zero reset" << endl;
#endif
  if (_queue_numrc_allocated) {
  // zero reset could be done redundantly
  for (int p = 0; p < num_threads; p++) {
     for (list<C_task_seq *>::const_iterator it = _queue_static[p].begin();
	  it != _queue_static[p].end(); ++it) {
       (*it)->referred = 0;
       vector<C_task *> &queue = *((*it)->queue);
       for (int j = (*it)->begin; j < (*it)->end; j++) {
	 queue[j]->referred = 0;
       }  // loop : j
     }    // loop : it
  } // loop : p
  for (vector<C_task_seq *>::const_iterator it = _queue_dynamic->begin();
       it != _queue_dynamic->end(); ++it) {
    (*it)->referred = 0;
    for (int j = (*it)->begin; j < (*it)->end; j++) {
      (*(*it)->queue)[j]->referred = 0;
    }
  }
#ifdef DEBUG_ERASE
  cerr << "increment referred counter" << endl;
#endif
  // increment referred counter
  for (int p = 0; p < num_threads; p++) {
     for (list<C_task_seq *>::const_iterator it = _queue_static[p].begin();
	  it != _queue_static[p].end(); ++it) {
       (*it)->referred++;
       vector<C_task *> &queue = *((*it)->queue);
       for (int j = (*it)->begin; j < (*it)->end; j++) {
	 queue[j]->referred++;
       }  // loop : j
     }    // loop : it
  } // loop : p
  for (vector<C_task_seq *>::const_iterator it = _queue_dynamic->begin();
       it != _queue_dynamic->end(); ++it) {
    (*it)->referred++;
    for (int j = (*it)->begin; j < (*it)->end; j++) {
      (*(*it)->queue)[j]->referred++;
    }
  }
  // delete task
#ifdef DEBUG_ERASE
  cerr << "delete task" << endl;
#endif
  for (int p = 0; p < num_threads; p++) {
     for (list<C_task_seq *>::const_iterator it = _queue_static[p].begin();
	  it != _queue_static[p].end(); ++it) {
       vector<C_task *> &queue = *((*it)->queue);
       for (int j = (*it)->begin; j < (*it)->end; j++) {
	 C_task *task = queue[j];
	 task->referred--;
	 if (task->referred == 0) {
	   erase_task<T, U>(task); // 
	 }
       }  // loop : j
     }    // loop : it
  } // loop : p

  for (vector<C_task_seq *>::const_iterator it = _queue_dynamic->begin();
       it != _queue_dynamic->end(); ++it) {
    for (int j = (*it)->begin; j < (*it)->end; j++) {
      C_task *task = (*(*it)->queue)[j];
      task->referred--;
      if (task->referred == 0) {
	erase_task<T, U>(task); // 
      }
    }
  }
  // delete task_seq
#ifdef DEBUG_ERASE
  cerr << "static : C_task_seq" << endl;
#endif
  // _queue_static[] share tasks in _queue_dynamic
  for (int p = 0; p < num_threads; p++) {
     for (list<C_task_seq *>::iterator it = _queue_static[p].begin();
	  it != _queue_static[p].end(); ++it) {
       (*it)->referred--;
       if ((*it)->referred == 0) {
#ifdef DEBUG_ERASE
	 fprintf(stderr, "%x %s\n", (*it), (*it)->task_name);
#endif
	 delete (*it)->queue;
	 (*it)->queue = NULL;
	 delete (*it);
	 (*it) = NULL;
       }
     } // loop : it
  } // loop : p
#ifdef DEBUG_ERASE
  cerr << "dynamic : C_task_seq" << endl;
#endif
  for (vector<C_task_seq *>::iterator it = _queue_dynamic->begin();
       it != _queue_dynamic->end(); ++it) {
    (*it)->referred--;
    if ((*it)->referred == 0) {
#ifdef DEBUG_ERASE
      fprintf(stderr, "%x %s\n", (*it), (*it)->task_name);
#endif
      delete (*it)->queue;
      (*it)->queue = NULL;
      delete (*it);
      (*it) = NULL;
    }
  } // loop : it

  delete _queue_dynamic;
  // Here all queues are deallocated.
  for (int p = 0; p < num_threads; p++) {
    _queue_static[p].clear();
  }
  delete [] _queue_static;
#ifdef DEBUG_ERASE
  cerr << "queues are deallocated" << endl;
#endif
  }
  _queue_symb_allocated = false;
  _queue_numrc_allocated = false;
}

template
void DissectionQueue<double, double>::erase_queue(void);

template
void DissectionQueue<quadruple, quadruple>::erase_queue(void);

template
void DissectionQueue<complex<double>, double>::erase_queue(void);

template
void DissectionQueue<complex<quadruple>, quadruple>::erase_queue(void);

//

template<typename T, typename U>
void DissectionQueue<T, U>::erase_queue_fwbw(void)
{
  if (_verbose) {
    fprintf(_fp, "%s %d : void QueueRuntime::erase_queue_fwbw(void)",
	    __FILE__, __LINE__);
  }
  if (_queue_fwbw_allocated) {
  delete [] _diag_contribs;
  delete _x;
  for (int d = 1; d <= _nb_doms; d++) {
    const int d0 = _btree->selfIndex(d);
    delete _xi[d0];
    delete _yi[d0];
    delete _wi[d0];
    delete _zi[d0];
  }

  delete [] _xi;
  delete [] _yi;
  delete [] _zi;
  delete [] _wi;
  
  for (vector<C_task *>::iterator it = _queue_fwbw->queue->begin();
       it != _queue_fwbw->queue->end(); ++it) {
    erase_task<T, U>(*it);
  }
  _queue_fwbw->queue->clear();
  delete _queue_fwbw->queue;
  delete _queue_fwbw;

  delete _nrhs; // added 01 Oct.2013 Atsushi
  delete _isTrans;
  }
  _queue_fwbw_allocated = false;
}

template
void DissectionQueue<double, double>::erase_queue_fwbw(void);

template
void DissectionQueue<quadruple, quadruple>::erase_queue_fwbw(void);

template
void DissectionQueue<complex<double>, double>::erase_queue_fwbw(void);


template
void DissectionQueue<complex<quadruple>, quadruple>::erase_queue_fwbw(void);

