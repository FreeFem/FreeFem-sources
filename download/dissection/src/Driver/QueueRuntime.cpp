/*! \file QueueRuntime.cpp
    \brief  management of threads for Dissection Matrix
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

#include "Driver/QueueRuntime.hpp"

#include <string>

#if 0
#define DEBUG_MUTEX
#define DEBUG_MUTEX1
#define DEBUG_MUTEX2
#define DEBUG_MUTEX3
#define DEBUG_STATIC_QUEUE1
#define DEBUG_STATIC_QUEUE2
#define DEBUG_MUTEX_C_BLAS_BLOCK1
#define DEBUG_MUTEX_C_DFULLLDLT
#endif

#define RATIO_STATIC 0.7       // for slow DGEMM'
#define RATIO_QUEUE_GREEDY 0.8

// #define DEBUG_PRINT_TASK
// #define DEBUG_EXEC_THREAD_21Jun2012
// #define DEBUG_THREAD_DONE_PRINT
// #define DEBUG_CHECKPARENTS_DONE 
// #define DEBUG_EXEC_THREAD
// #define DEBUG_EXEC_THREAD_IDLE
// #define DEBUG_PREPARE_THREAD
#if 0
#define DEBUG_EXEC_THREAD    
#define DEBUG_THREAD_LOOP
#define DEBUG_THREAD_DONE_PRINT
#define DEBUG_EXEC_THREAD_IDLE
#define DEBUG_EXEC_THREAD_FILE
#define DEBUG_CHECKPARENTS_DONE 
#endif

//#define DEBUG_THREAD_TIME
//#define DEBUG_EXEC_THREAD
#define DEBUG_DEADLOCK

// DEBUG_DEADLOCK only works without setting DEBUG_EXEC_THREAD_IDLE
#ifdef DEBUG_DEADLOCK
#ifdef DEBUG_EXEC_THREAD_IDLE
#undef DEBUG_EXEC_THREAD_IDLE
#endif // DEBUG_EXEC_THREAD_IDLE
#endif //  DEBUG_DEADLOCK

extern pthread_mutex_t _mutex_debug; 

// constructor
QueueRuntime::QueueRuntime(int nb_doms, int num_threads,
			   const bool verbose, FILE *fp)
{

  _nb_doms = nb_doms;
  _num_threads = num_threads;
  _verbose = verbose;
  _fp = fp;

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
}


QueueRuntime::~QueueRuntime()
{
  // initialize state to controle tasks
  delete [] _zone_entered;
  delete [] _zone_finished;
  delete [] _zone_static_assigned;

  deallocate_int2d(_begins);
  deallocate_int2d(_ends);
  deallocate_int3d(_group_entered);
  deallocate_int3d(_group_finished);
  deallocate_int3d(_group_task_ends);
  deallocate_int3d(_group_static_assigned);
  deallocate_int3d(_begins_group);
  deallocate_int3d(_ends_group);
  deallocate_unsigned2d(_group_nops);

  delete [] _group_task_id;
  delete [] _mutex_group;
}


void allocate_int2d(int **&array, const int num_threads)
{
  array = new int*[num_threads];
  array[0] = new int[num_threads * DIST_TASK_CRITICAL];
  int p, pp;
  p = 1;
  pp = DIST_TASK_CRITICAL;
  for ( ; p < num_threads; p++, pp += DIST_TASK_CRITICAL) {
    array[p] = &array[0][pp];
  }
}

void allocate_unsigned2d(long **&array, const int num_threads)
{
  array = new long*[num_threads];
  array[0] = new long[num_threads * DIST_TASK_CRITICAL];
  int p, pp;
  p = 1;
  pp = DIST_TASK_CRITICAL;
  for ( ; p < num_threads; p++, pp += DIST_TASK_CRITICAL) {
    array[p] = &array[0][pp];
  }
}

void allocate_int3d(int ***&array, const int num_threads)
{
  array = new int**[num_threads];
  array[0] = new int*[num_threads * DIST_TASK_CRITICAL];
  array[0][0] = new int[(num_threads * 
			 DIST_TASK_CRITICAL * DIST_TASK_CRITICAL)];
  int p, q, pp;
  p = 1;
  pp = DIST_TASK_CRITICAL;
  for ( ; p < num_threads; p++, pp += DIST_TASK_CRITICAL) {
    array[p] = &array[0][pp];
  }
  p = 0;
  pp = 0;
  for ( ; p < num_threads; p++, pp += DIST_TASK_CRITICAL) {
    for (q = 0; q < DIST_TASK_CRITICAL; q++) {
      if (!((p == 0) && (q == 0))) {
	  array[p][q] = &array[0][0][(q + pp) * DIST_TASK_CRITICAL];
      }
    }
  }
}

void deallocate_int2d(int **&array)
{
  delete [] array[0];
  delete [] array;
}

void deallocate_unsigned2d(long **&array)
{
  delete [] array[0];
  delete [] array;
}

void deallocate_int3d(int ***&array)
{
  delete [] array[0][0];
  delete [] array[0];
  delete [] array;
}

int IndexTask(vector<C_task *> &queue, int pos)
{
  int j = 0;
  int pos1 = pos > queue.size() ? queue.size() : pos;
  for (int i = 0; i < pos1; i++) {
    j += queue[j]->atomic_size;
  }
  return j;
}

// #define DEBUG_NULL_TASK
void copytask_list2seq(list<C_task_seq*> &queue_static,
		       list<C_task_seq*> queue_lists,
		       list<C_task *> &queue_null,
//		       list<C_task *> &queue_dummy,
		       string task_name,
		       int task_id, 
		       int mutex_id,
		       int parallel_single,
		       int num_threads,
		       int level,
		       int phase)
{
  // #define DEBUG_ASSIGN_TASK
#ifdef DEBUG_ASSIGN_TASK
  cerr << endl << "** " << task_name << " ** " << queue_lists.size() << endl;
#endif
  int num_tasks = 0;
  for (list<C_task_seq*>::const_iterator it = queue_lists.begin();
       it != queue_lists.end(); ++it) {
    //    int itmp = 0;
    for (int j = (*it)->begin; j < (*it)->end; j++) {
      if (*((*(*it)->queue)[j]->ops_complexity) != 0L) {//cf. -1L : 16 Sep.2014
	num_tasks++;
//    	itmp++;
      }
#if 0
      else if (*((*(*it)->queue)[j]->ops_complexity) == (-1L)) {
	fprintf(stderr, "%s %d : copytask_list2seq() dummy task = %s\n",
		__FILE__, __LINE__, (*(*it)->queue)[j]->task_name);
	queue_dummy.push_back((*(*it)->queue)[j]);
      }
#endif
      else {
	(*(*it)->queue)[j]->status = TASK_DONE;
#ifdef DEBUG_NULL_TASK
	fprintf(stderr, "%s %d : copytask_list2seq() null task = %s\n",
		__FILE__, __LINE__, (*(*it)->queue)[j]->task_name);
#endif
	queue_null.push_back((*(*it)->queue)[j]);
      }
    }
#ifdef DEBUG_ASSIGN_TASK
    cerr << (*it)->task_id
	 << " [ " << (*it)->end << " " << (*it)->begin << " ] ";
#endif
  }
#ifdef DEBUG_ASSIGN_TASK
  cerr << endl;
#endif
  vector<C_task *> *tasks_tmp = new vector<C_task *>;
  long nops = 0L;
  tasks_tmp->resize(num_tasks);	
  // tasks_tmp->reserve(num_tasks);
  int k = 0;
  for (list<C_task_seq*>::const_iterator it = queue_lists.begin();
       it != queue_lists.end(); ++it) {
    for (int j = (*it)->begin; j < (*it)->end; j++) {
      if (*((*(*it)->queue)[j]->ops_complexity) != 0L) { //cf. -1L : 16 Sep.2014
	(*tasks_tmp)[k++] = (*(*it)->queue)[j];
	if (*((*(*it)->queue)[j]->ops_complexity) > 0L) {
	  nops += *((*(*it)->queue)[j]->ops_complexity);
	}
      }
#if 0
      else if (*((*(*it)->queue)[j]->ops_complexity) == 0L) {
	(*(*it)->queue)[j]->status = TASK_DONE;
#ifdef DEBUG_NULL_TASK
	fprintf(stderr, "%s %d : copytask_list2seq() null task = %s\n",
		__FILE__, __LINE__, (*(*it)->queue)[j]->task_name);
#endif
	queue_null.push_back((*(*it)->queue)[j]);
      }
#endif
    }
  }
  C_task_seq* seq_tmp = new C_task_seq(task_id,
				       task_name,
				       mutex_id,
				       parallel_single,
				       num_threads,
				       level,
				       phase,
				       tasks_tmp,
				       0,                  // begin
				       num_tasks,  // end
				       nops);
  queue_static.push_back(seq_tmp);
}

void task_assign_diag1(list<C_task_seq*> * &queue_static,
		       vector<C_task *>* &tasks_queue,
		       list<C_task *> &queue_null,
     //		       list<C_task *> &queue_dummy,
		       vector<int> &nrow_DFullLDLt,
		       vector<bool> &isMergedDTRSM,
		       string queue_symbol,
		       int queue_id,
		       int level_id,
		       int num_threads,
		       int level_last,
		       int level,
		       int begdom,
		       int enddom,
		       long *nops_sum,
		       vector<int> starts,
		       list<C_task_seq*>* &queue_lists,
		       bool queue_lists_clear,
		       Dissection::Tree* btree,
		       vector<int>* children)
{
  int num_tmp_shared;
  list<C_task_seq *> tmp_dynamic;
  list<C_task_seq *> task_seq_tmp;
  long nops_block_total;
  vector<C_task *> *tmp_shared = new vector<C_task *>;

  if (queue_lists_clear) {
    for (int p = 0; p < num_threads; p++) {
#if 1
      for (list<C_task_seq *>::iterator it = queue_lists[p].begin();
	   it != queue_lists[p].end(); ++it) {
	delete (*it);
	(*it) = NULL;
      }
#endif
      queue_lists[p].clear();
    }
  }
  task_seq_tmp.clear();
  nops_block_total = 0L;
  for (int d = begdom; d < enddom; d++) {
    const int j = btree->selfIndex(d);
    for (int ic = 0; ic < 2; ic++) {
      const int jc = children[j][ic];
      if ((queue_id == C_DTRSM) && (nrow_DFullLDLt[jc] > 2)) {
	continue;
      }
      long nops = 0L;
      const int begin = starts[jc];
      int end;
      if (queue_id == C_DTRSM) {
	if (!isMergedDTRSM[j]) { // looking the father node
	  end = IndexTask(tasks_queue[jc], tasks_queue[jc][0]->parallel_max);
	}
	else {
	  end = tasks_queue[jc].size();
	}
      }
      else {
	int itmp = tasks_queue[jc][0]->parallel_max;
	int jtmp = tasks_queue[jc].size();
	//	if (itmp < (_isSym ? 3 : 4)) {
	if (!isMergedDTRSM[j]) { 
	  end = itmp;
	}
	else {
	  if (itmp == jtmp) {
	    end = jtmp;
	  }
	  else {
	    end = itmp + tasks_queue[jc][itmp]->parallel_max;
	  }
	}
      }
      if (level == 0) {
	end = (int)((double)end * RATIO_STATIC);
	const int ktmp = tasks_queue[jc][end]->atomic_id;
	end -= ktmp;
      }
      for (int i = begin; i < end; i++) {
	nops += *(tasks_queue[jc][i]->ops_complexity);
      }
      string task_name = to_string(jc);
      if (nops > 0L) {
	C_task_seq* tmp = 
	  new C_task_seq(queue_id,
			 task_name,  // for debugging
			 (-1), // mutex_id
			 TASK_SINGLE,
			 1,
			 (level + 1),
			 (level_last - 1 - level) * 6 + level_id,
			 &tasks_queue[jc],
			 begin,
			 end,
			 nops);
	task_seq_tmp.push_back(tmp);
	nops_block_total += nops;
      }
      else {
	for (int i = begin; i < end; i++) {
	  //	  fprintf(stderr, "%s %d : %s\n",
	  //		  __FILE__, __LINE__, tasks_queue[jc][i]->task_name);
	  queue_null.push_back(tasks_queue[jc][i]);
	}
      }
    }  // loop : ic
  }    // loop : j
  num_tmp_shared = 0;
  tmp_dynamic.clear();

  if (task_seq_tmp.size() > 0) { 
    task_seq_tmp.sort(C_task_seq_complexity_greater);
    string queue_name;
    if (level == 0) {
      queue_name = queue_symbol + "% ";
    }
    else {
      queue_name = queue_symbol + "* ";
    }
    assign_tasks_statically(queue_lists,
			    tmp_dynamic,
			    nops_sum,
			    task_seq_tmp,
			    queue_id,
			    queue_name.c_str(),
			    (level + 1),
			    (level_last - 1 - level) * 6 + level_id,
			    nops_block_total,
			    num_threads);
    //    task_seq_tmp.clear(); moved in  assign_tasks_statically
    if (tmp_dynamic.size() > 0) {
      num_tmp_shared += tmp_dynamic.size();
    }
  }
  if (level == 0) {
    for (int d = begdom; d < enddom; d++) {
      const int j = btree->selfIndex(d);
      for (int ic = 0; ic < 2; ic++) {
	const int jc = children[j][ic];
	int parallel_max;
	if ((queue_id == C_DTRSM) && (nrow_DFullLDLt[jc] > 2)) {
	  continue;
	}
   	if (queue_id == C_DTRSM) {
	  if (!isMergedDTRSM[j]) { // looking father node
	    parallel_max = IndexTask(tasks_queue[jc],
				     tasks_queue[jc][0]->parallel_max);
	  }
	  else {
	    parallel_max = tasks_queue[jc].size();
	  }
	}
	else {
	  int itmp = tasks_queue[jc][0]->parallel_max;
	  int jtmp = tasks_queue[jc].size();
	  //	  if (itmp < _isSym ? 3 : 4) {
	if (!isMergedDTRSM[j]) { 
	    parallel_max = itmp;
	  }
	  else {
	    parallel_max = ((itmp == jtmp) ? jtmp : 
			    (itmp + tasks_queue[jc][itmp]->parallel_max));
	  }
	}
	const int end = (int)((double)parallel_max * RATIO_STATIC);
	const int ktmp = tasks_queue[jc][end]->atomic_id;
	num_tmp_shared += parallel_max - end + ktmp;
      }
    }
  }
  tmp_shared->reserve(num_tmp_shared);
  if (tmp_dynamic.size() > 0) {
    long nops = 0L;
    for (list<C_task_seq *>::iterator mt = tmp_dynamic.begin();
	 mt != tmp_dynamic.end(); ++mt) {
      nops += (*mt)->ops_complexity;
    }
#ifdef DEBUG_PREPARE_THREAD
    cout << "static -> dynamic : " << tmp_dynamic.size() 
	 << " nops_sum = " << nops << " " ;
#endif
    for (list<C_task_seq *>::iterator mt = tmp_dynamic.begin();
	 mt != tmp_dynamic.end(); ++mt) {
      for (int i = (*mt)->begin; i < (*mt)->end; i++) {
	tmp_shared->push_back((*(*mt)->queue)[i]);
      }
      delete (*mt);
      (*mt) = NULL;
    }
#ifdef DEBUG_PREPARE_THREAD
    cout << endl;
#endif
  } // if (tmp_dynamic.size() > 0)
  if (level == 0) {
    for (int d = begdom; d < enddom; d++) {
      const int j = btree->selfIndex(d);
      for (int ic = 0; ic < 2; ic++) {
	const int jc = children[j][ic];
	int parallel_max;
	if ((queue_id == C_DTRSM) && (nrow_DFullLDLt[jc] > 2)) {
	  continue;
	}
	if (queue_id == C_DTRSM) {
	  if (!isMergedDTRSM[j]) { // looking father node
	    parallel_max = IndexTask(tasks_queue[jc],
				     tasks_queue[jc][0]->parallel_max);
	  }
	  else {
	    parallel_max = tasks_queue[jc].size();
	  }
	}
	else {
	  int itmp = tasks_queue[jc][0]->parallel_max;
	  int jtmp = tasks_queue[jc].size();
	  //	  if (itmp <= (_isSym ? 3 : 4)) {
	  if (!isMergedDTRSM[j]) { 
	    parallel_max = itmp;
	  }
	  else {
	    parallel_max = ((itmp == jtmp) ? jtmp : 
			    (itmp + tasks_queue[jc][itmp]->parallel_max));
	  }
	}
	int end = (int)((double)parallel_max * RATIO_STATIC);
	const int ktmp = tasks_queue[jc][end]->atomic_id;
	end -= ktmp;
	for (int i = end; i < parallel_max; i++) {
	                                                 // exclude null tasks
	  if (*(tasks_queue[jc][i]->ops_complexity) > 0L) { 
	    tmp_shared->push_back(tasks_queue[jc][i]);
	  }
	  else {
	    queue_null.push_back(tasks_queue[jc][i]);
	  }
	  // cyclic by index jc is better ?
	}
      } // loop : ic
    }
  }
  {
    string task_name = queue_symbol + " " + to_string(level + 1);
    for (int p = 0; p < num_threads; p++) {
      copytask_list2seq(queue_static[p],
			queue_lists[p],
			queue_null,
			//			queue_dummy,
			task_name,
			queue_id,
			p,
			TASK_SINGLE,
			1,
			(level + 1),
			(level_last - 1 - level) * 6 + level_id);
      for (list<C_task_seq *>::iterator it = queue_lists[p].begin();
	   it != queue_lists[p].end(); ++it) {
	delete (*it);
	(*it) = NULL;
      }
    }
  }
#ifdef DEBUG_PREPARE_THREAD
  cout << queue_symbol << num_tmp_shared << " : " << tmp_shared->size() 
       << endl;
#endif
  if (tmp_shared->size() > 0) {
    string task_name = queue_symbol + "+ " + to_string(level + 1);
    long nops = 0L;
    for (vector<C_task *>::const_iterator it = tmp_shared->begin();
	 it != tmp_shared->end(); ++it) {
      nops += *((*it)->ops_complexity);
    }
    C_task_seq* tmp = 
      new C_task_seq(queue_id,
		     task_name,
		     0,
		     TASK_PARALLEL,
		     num_threads,
		     (level + 1),
		     (level_last - 1 - level) * 6 + level_id,
		     tmp_shared,
		     0,
		     tmp_shared->size(),
		     nops);
    for (int p = 0; p < num_threads; p++) {
      queue_static[p].push_back(tmp);
    }
  }  // if (num_tmp_shared > 0)
  else {
    delete tmp_shared;
  }
  //  delete tmp_shared;
}

void task_assign_diag2(list<C_task_seq*> *&queue_static,
		       vector<C_task *>* &tasks_queue, 
		       list<C_task *> &queue_null,
		       //		       list<C_task *> &queue_dummy,
		       vector<bool> &isMergedDTRSM,
		       string queue_symbol,
		       int queue_id,
		       int level_id,
		       int num_threads,
		       int level_last,
		       int level,
		       int begdom,
		       int enddom,
		       long *nops_sum,
		       vector<int> starts,
		       vector<int> starts_sub,
		       list<C_task_seq*>* &queue_lists,
		       bool queue_lists_clear,
		       Dissection::Tree* btree,
		       vector<int>* children)

{
  int num_tmp_shared;
  list<C_task_seq*> tmp_dynamic;
  list<C_task_seq *>task_seq_tmp;
  long nops_block_total;
  vector<C_task *> *tmp_shared = new vector<C_task *>;

  task_seq_tmp.clear();
  nops_block_total = 0L;
  if (queue_lists_clear) {
    for (int p = 0; p < num_threads; p++) {
#if 1
      for (list<C_task_seq *>::iterator it = queue_lists[p].begin();
	   it != queue_lists[p].end(); ++it) {
	delete (*it);
	(*it) = NULL;
      }
#endif
      queue_lists[p].clear();
    }
  }
  for (int d = begdom; d < enddom; d++) {
    const int j = btree->selfIndex(d);
    long nops = 0L;
    const int begin = starts_sub[j];
    int end;
    //   if (tasks_queue0[j].size() <= 4) {
    if (!isMergedDTRSM[j]) {
      end = tasks_queue[j][0]->parallel_max;
    }
    else {
      end = tasks_queue[j].size();
    }
    if (level == 0) {
      end = (int)((double)end * RATIO_STATIC);
    }
    for (int i = begin; i < end; i++) {
      nops += *(tasks_queue[j][i]->ops_complexity);
    }
    string task_name = to_string(j);
    if (nops > 0L) {
      C_task_seq* tmp = 
	new C_task_seq(queue_id,
		       task_name, // for debugging
		       (-1), // mutex_id
		       TASK_SINGLE,
		       1,
		       (level + 1),
		       (level_last - 1 - level) * 6 + level_id,
		       &tasks_queue[j],
		       begin,
		       end,
		       nops);
      task_seq_tmp.push_back(tmp);
      nops_block_total += nops;
    }
    else {
      for (int i = begin; i < end; i++) {
	queue_null.push_back(tasks_queue[j][i]);
      }
    }
  } // loop : j
  tmp_dynamic.clear();
  num_tmp_shared = 0;
  if (task_seq_tmp.size() > 0) {
    task_seq_tmp.sort(C_task_seq_complexity_greater);
    string queue_name;
    if (level == 0) {
      queue_name = queue_symbol + "% " ;
    }
    else {
      queue_name = queue_symbol + "* " ;
    }
    assign_tasks_statically(queue_lists,
			    tmp_dynamic,
			    nops_sum,
			    task_seq_tmp,
			    queue_id,
			    queue_name.c_str(), 
			    (level + 1),
			    (level_last - 1 - level) * 6 + level_id,
			    nops_block_total,
			    num_threads);
    if (tmp_dynamic.size() > 0) {
      num_tmp_shared += tmp_dynamic.size();
    }
  }
  if (level == 0) {
    for (int d = begdom; d < enddom; d++) {
      const int j = btree->selfIndex(d);
      const int parallel_max = tasks_queue[j].size();
      const int end = (int)((double)parallel_max * RATIO_STATIC);
      num_tmp_shared += parallel_max - end;
    }
  }
  tmp_shared->reserve(num_tmp_shared);
  if (tmp_dynamic.size() > 0) {
    long nops = 0L;
    for (list<C_task_seq *>::const_iterator mt = tmp_dynamic.begin();
	 mt != tmp_dynamic.end(); ++mt) {
      nops += (*mt)->ops_complexity;
    }
    for (list<C_task_seq *>::iterator mt = tmp_dynamic.begin();
	 mt != tmp_dynamic.end(); ++mt) {
      for (int i = (*mt)->begin; i < (*mt)->end; i++) {
	tmp_shared->push_back((*(*mt)->queue)[i]);
      }
      delete (*mt);
      (*mt) = NULL;
    }
  }
  if (level == 0) {
    for (int d = begdom; d < enddom; d++) {
      const int j = btree->selfIndex(d);
      int parallel_max;
      if (!isMergedDTRSM[j]) {
	parallel_max = tasks_queue[j][0]->parallel_max;
      }
      else {
	parallel_max = tasks_queue[j].size();
      }
      const int end = (int)((double)parallel_max * RATIO_STATIC);
      for (int i = end; i < parallel_max; i++) {
	if (*(tasks_queue[j][i]->ops_complexity) > 0L) {  // exclude null tasks
	  tmp_shared->push_back(tasks_queue[j][i]);
	}
	else {
	  queue_null.push_back(tasks_queue[j][i]);
	}
	// cyclic by index jc is better ?
      }
    }
  }
  {
    string task_name = queue_symbol + " " + to_string(level + 1);
    for (int p = 0; p < num_threads; p++) {
      copytask_list2seq(queue_static[p],
			queue_lists[p],
			queue_null,
//			queue_dummy,
			task_name,
			queue_id,
			p,
			TASK_SINGLE,
			1,
			(level + 1),
			(level_last - 1 - level) * 6 + level_id);
      // 3 Jul 2016
      for (list<C_task_seq *>::iterator it = queue_lists[p].begin();
	   it != queue_lists[p].end(); ++it) {
	delete (*it);
	(*it) = NULL;
      }
    }
  }
#ifdef DEBUG_PREPARE_THREAD
  cout << queue_symbol << num_tmp_shared << " : " << tmp_shared->size() 
       << endl;
#endif
  if (tmp_shared->size() > 0) {
    string task_name = queue_symbol + "+ " + to_string(level + 1);
    long nops = 0L;
    for (vector<C_task *>::const_iterator it = tmp_shared->begin();
	 it != tmp_shared->end(); ++it) {
      nops += *((*it)->ops_complexity);
    }
    C_task_seq* tmp = 
      new C_task_seq(queue_id,
		     task_name,
		     0,
		     TASK_PARALLEL,
		     num_threads,
		     (level + 1),
		     (level_last - 1 - level) * 6 + level_id,
		     tmp_shared,
		     0,
		     tmp_shared->size(),
		     nops);
    for (int p = 0; p < num_threads; p++) {
      queue_static[p].push_back(tmp);
    }
  } //  if (tmp_shared->size() > 0)
  else {
    delete tmp_shared;
  }
}

void QueueRuntime::generate_queue(C_task_seq* &queue_symb_,
				  list<C_task_seq*>* & queue_static_,
				  vector<C_task_seq *>* &queue_dynamic_,
				  list<C_task *> &queue_null_,
  //				  list<C_task *> &queue_dummy_,
				  Dissection::Tree* btree,
				  vector<int>* children,
				  vector<C_task *>* tasks_SparseSymb,
				  vector<C_task *>* tasks_SparseNum,
				  vector<C_task *>* tasks_SparseLocalSchur,
				  vector<C_task *>* tasks_DFillSym,
				  vector<C_task *>* tasks_DFullLDLt,
				  vector<C_task *>* tasks_DTRSMScale,
				  vector<C_task *>* tasks_DGEMM,
				  vector<C_task *>** tasks_Dsub,
				  vector<C_task *>* tasks_deallocLower,
				  vector<C_task *>* tasks_deallocLocalSchur,
				  long **nops_queue,
				  vector<vector<int> > all_fathersIndex,
				  vector<int> nrow_DFullLDLt,
				  vector<bool> isMergedDTRSM,
				  vector<bool> isDividedDTRSM,
				  const int level_last)
{  
  // begin : task squence for symmbolic factorization
  const int num_threads = _num_threads;
  const int nb_doms_dense = (1U << level_last) - 1;

  {  //begin scope : task_name
    list<C_task_seq*> queue_tmp;
    string task_name = "l " + to_string(level_last);

    const int begdom = 1U << level_last; 
    const int enddom = begdom * 2;

    for (int d = begdom; d < enddom; d++) {
      const int j = d - begdom;
      C_task_seq* tmp = 
	new C_task_seq(C_SPARSESYMBFACT,
		       _null_name, // dummy
		       (-1), // mutex_id
		       TASK_SINGLE, 
		       1,
		       level_last,
		       (-1), // phase
		       &tasks_SparseSymb[j],
		       0,
		       1,
		       *(tasks_SparseSymb[j][0]->ops_complexity));
      queue_tmp.push_back(tmp);
    }
    //
    queue_tmp.sort(C_task_seq_complexity_greater); 
    vector<C_task *> *task_tmp = new vector<C_task *>;
    task_tmp->resize(begdom);
    int k = 0;
    for (list<C_task_seq*>::iterator it = queue_tmp.begin();
	 it != queue_tmp.end(); ++it, k++) {
      (*task_tmp)[k] = (*it)->queue[0][0];
    }
    queue_symb_ = new C_task_seq(C_SPARSESYMBFACT,
				 task_name, // dummy
				 (-1), // mutex_id
				 TASK_SINGLE,
				 1,
				 level_last,
				 (-1), // phase
				 task_tmp,
				 0,
				 begdom,
				 0); // dummy
    _queue_symb = queue_symb_;
    for (list<C_task_seq*>::iterator it = queue_tmp.begin();
	 it != queue_tmp.end(); ++it) {
      delete (*it);
      (*it) = NULL;
    }
  } //end scope : task_name
  queue_static_ = new list<C_task_seq*>[num_threads];
  _queue_static = queue_static_;
  queue_dynamic_ = new vector<C_task_seq*>;
  _queue_dynamic = queue_dynamic_;

// end : task squence for symmbolic factorization
  list<C_task_seq*>* queue_dynamic = new list<C_task_seq*>[level_last + 1];
  list<C_task_seq*>* queue_dynamic0 = new list<C_task_seq*>;
  list<C_task_seq*>* queue_dynamic1 = new list<C_task_seq*>[level_last + 1];
  list<C_task_seq*>* queue_dynamic2 = new list<C_task_seq*>[level_last + 1];
  //  list<C_task_seq*>* queue_dynamic3 = new list<C_task_seq*>[level_last + 1];
  list<C_task_seq*>* queue_lists = new list<C_task_seq*>[num_threads];
  list<C_task_seq*> queue_tmp0, queue_tmp1;
  // assuming that num_leaves(level_last) > num_threads
  // num_leaves(level_last) = _begLevel[level_last] - _begLevel[level_last + 1]
  // 
// Sparse numeric factorization : tasks_SparseNum[] -> _queue_static[]
  queue_tmp0.clear();
#define SPARSE_MERGED
#ifdef SPARSE_MERGED
  {
  // sparse subdomains    : tasks_SparseNum[] 
  // all dense subdomains : tasks_DFillSym[]
  // the begining of the dense subdomains
    const int begdom = 1U << (level_last - 1); 
    const int enddom = begdom * 2;
    for (int d = begdom; d < enddom; d++) {
      int j;
      vector<C_task *>* tasks = new vector<C_task *>;
      j = btree->selfIndex(2 * d);
      long nops = 0L;
      tasks->push_back(tasks_SparseNum[j][0]);
      nops += *(tasks_SparseNum[j][0]->ops_complexity);
      tasks->push_back(tasks_SparseLocalSchur[j][0]);
      nops += *(tasks_SparseLocalSchur[j][0]->ops_complexity);
      j = btree->selfIndex(2 * d + 1);
      tasks->push_back(tasks_SparseNum[j][0]);
      nops += *(tasks_SparseNum[j][0]->ops_complexity);
      tasks->push_back(tasks_SparseLocalSchur[j][0]);
      nops += *(tasks_SparseLocalSchur[j][0]->ops_complexity);
      j = btree->selfIndex(d);
      for (vector<C_task *>::const_iterator it = tasks_DFillSym[j].begin();
	   it != tasks_DFillSym[j].end(); ++it) {
	tasks->push_back(*it);
	nops += *((*it)->ops_complexity);
      }
      // task complexity of DFillSym[] in higher level is very large
      int dd;
      if (d % 2 == 0) {
	dd = (d - begdom) / 2;
      }
      else {
      	dd = begdom - (d - begdom + 1) / 2;
      }
      if (dd > 0) {
	j = btree->selfIndex(dd);
	for (vector<C_task *>::const_iterator it = tasks_DFillSym[j].begin();
	     it != tasks_DFillSym[j].end(); ++it) {
	  tasks->push_back(*it);
	  nops += *((*it)->ops_complexity);
	}
      }
      C_task_seq* tmp = 
	new C_task_seq(C_SPARSESOLVER,
		       _null_name,
		       (-1), // mutex_id
		       (-1), //TASK_SINGLE,
		       (-1), // 1,
		       (-1), // level_last,
		       (-1), // phase
		       tasks,
		       0,
		       tasks->size(),
		       nops);
      queue_tmp0.push_back(tmp);
    } // loop : d(j)
    string task_name = "z " + to_string(level_last);
    copytask_list2seq(queue_dynamic[level_last],
		      queue_tmp0,
		      queue_null_,
     //		      queue_dummy_,
		      task_name,
		      C_SPARSESOLVER, // C_FILLMATRIX1,
		      0,               // shared by all threads
		      TASK_PARALLEL,
		      num_threads,
		      level_last,
		      2);
  // erase temporary C_task_seq whose elements are copied to queue_dynamic[]
    for (list<C_task_seq*>::iterator it = queue_tmp0.begin();
	 it != queue_tmp0.end(); ++it) {
      delete (*it)->queue; // deallocate a container : tasks in L880
      delete (*it);
      (*it) = NULL;
    } 
  }
#else
  {
  // sparse subdomains    : tasks_SparseNum[] 
  // all dense subdomains : tasks_DFillSym[]
  // the begining of the dense subdomains
    const int begdom = 1U << level_last; 
    const int enddom = begdom * 2;
    queue_tmp0.clear();
    for (int d = begdom; d < enddom; d++) {
      const int j = btree->selfIndex(d);
      C_task_seq* tmp = 
	new C_task_seq(C_SPARSENUMFACT,
		       _null_name,
		       (-1), // mutex_id
		       (-1), //TASK_SINGLE,
		       (-1), // 1,
		       (-1), // level_last,
		       (-1), // phase
		       &tasks_SparseNum[j],
		       0,
		       1,
		       *(tasks_SparseNum[j][0]->ops_complexity));
      queue_tmp0.push_back(tmp);
    } // loop : d(j)
    queue_tmp0.sort(C_task_seq_complexity_greater);
    {  
      int itmp = 0;
      for (list<C_task_seq*>::const_iterator it = queue_tmp0.begin();
	   it != queue_tmp0.end(); ++it, itmp++) {
	queue_lists[itmp % num_threads].push_back(*it);
      }
    } // loop : it, itmp
    {
      string task_name = "u " + to_string(level_last);
      for (int p = 0; p < num_threads; p++) {
	copytask_list2seq(_queue_static[p],
			  queue_lists[p],
			  queue_null_,
  //			  queue_dummy_,
			  task_name,
			  C_SPARSENUMFACT,
			  p,
			  TASK_SINGLE,
			  1,
			  level_last,
			  0);
      }
      for (list<C_task_seq*>::iterator it = queue_tmp0.begin();
	   it != queue_tmp0.end(); ++it) {
	delete (*it);
	(*it) = NULL;
      } 
    }
// Sparse local Schur complement : tasks_SparseLocalSchur[] -> queue_dynamic[]
    queue_tmp0.clear();
    for (int d = begdom; d < enddom; d++) {
      const int j = btree->selfIndex(d);
      C_task_seq* tmp = 
	new C_task_seq(C_SPARSELOCALSCHUR,
		       _null_name, //(char *)NULL, // dummy
		       (-1), // mutex_id
		       (-1), //TASK_SINGLE,
		       (-1), // 1,
		       (-1), // level_last,
		       (-1), // phase
		       &tasks_SparseLocalSchur[j],
		       0,
		       1,
		       *(tasks_SparseLocalSchur[j][0]->ops_complexity));
      queue_tmp0.push_back(tmp);
    } // loop : d(j)
    queue_tmp0.sort(C_task_seq_complexity_greater);
    {
      string task_name = "w " + to_string(level_last);

      // C_SPARSELOCALSCHUR is simple greedy estimated complexity by 
      //                    assuming blocks are dense
      // C_SPARSELOCALSCHUR1 uses estimated complexity by C_SPARSENUMFACT
      copytask_list2seq(*queue_dynamic0,
			queue_tmp0,
			queue_null_,
//			queue_dummy_,
			task_name,
			C_SPARSELOCALSCHUR, // C_SPARSELOCALSCHUR1 
			0,               // shared by all threads
			TASK_PARALLEL,
			num_threads,
			level_last,
			1);              // SparseLDLt : 0
      // erase temporary C_task_seq whose elements are copied to queue_dynamic[]
      for (list<C_task_seq*>::iterator it = queue_tmp0.begin();
	   it != queue_tmp0.end(); ++it) {
	delete (*it);
	(*it) = NULL;
      } 
    }
    for (int p = 0; p < num_threads; p++) {
      for(list<C_task_seq *>::const_iterator it = queue_dynamic0->begin();
	  it != queue_dynamic0->end(); ++it) {
	    _queue_static[p].push_back(*it);
      }
    }
  }
    // -- C_FILLMATRIX diagonal -
    // without using children
  queue_tmp0.clear();
  {
    // sparse subdomains    : tasks_SparseNum[] 
    // all dense subdomains : tasks_DFillSym[]
    // the begining of the dense subdomains
    const int begdom = 1U << (level_last  - 1); 
    const int enddom = begdom * 2;
    
    for (int p = 0; p < num_threads; p++) {
      queue_lists[p].clear();
    }
#if 0  // treat tasks of C_FILLMATRIX in greedy manner
    for (int d = begdom; d < enddom; d++) {
      const int j = btree->selfIndex(d);
//
#ifdef DEBUG_PREPARE_THREAD
      cout << "* level = " << (level + 1) << " j = " << j << endl;
#endif
      C_task_seq* tmp = 
	new C_task_seq(C_FILLMATRIX,
		       _null_name, // dummy
		       (-1), // mutex_id
		       (-1), // TASK_SINGLE,
		       (-1), // 1,
		       (-1), // (level + 1),
		       (-1), // (level_last - 1 - level) * 5 + 1,
		       &tasks_DFillSym[j],
		       0,
		       1,
		       *(tasks_DFillSym[j][0]->ops_complexity));
      
      queue_tmp0.push_back(tmp);
    }
    {
      int itmp = 0;
      queue_tmp0.sort(C_task_seq_complexity_greater);
      for (list<C_task_seq*>::const_iterator it = queue_tmp0.begin();
	   it != queue_tmp0.end(); ++it, itmp++) {
	queue_lists[itmp % num_threads].push_back(*it);
      }
    } // scope of itmp
    {  // begin : scope of task_name
      string task_name = "v " + to_string(level_last);
      for (int p = 0; p < num_threads; p++) {
	copytask_list2seq(_queue_static[p],
			  queue_lists[p],
			  queue_null_,
  //			  queue_dummy_,
			  task_name,
			  C_FILLMATRIX,
			  p,
			  TASK_SINGLE,
			  1,
			  level_last,
			  2);
      } // loop : p
// erase temporary C_task_seq whose elements are copied to queue_dynamic[]
      for (list<C_task_seq*>::iterator it = queue_tmp0.begin();
	   it != queue_tmp0.end(); ++it) {
	delete (*it);
	(*it) = NULL;
      }
    }  // end : scope of task_name
#endif
    {
      string task_name = "v " + to_string(level_last) + " : ";
      queue_tmp0.clear();
      for (int d = begdom; d < enddom; d++) {
	const int j = btree->selfIndex(d);
	const int begin = 0;  // better that all (diag/offdiag) in greedy
	//	const int begin = tasks_DFillSym[j][0]->parallel_max;
	// -- C_FILLMATRIX offdiagonal --
	const int end = tasks_DFillSym[j].size();
	if ((end - begin) > 0) {   // always true
	  task_name += to_string(d) + " ";
	  C_task_seq* tmp = 
	    new C_task_seq(C_FILLMATRIX,
			   _null_name,      //	     task_name_cstr,
			   (-1),           // mutex_id
			   (-1),           // TASK_PARALLEL,
			   (-1),           // num_threads,
			   (-1),           //        level
			   (-1),           //         phase
			   &tasks_DFillSym[j],
			   begin,
			   end,
			   0L // nops;
			   );
	  queue_tmp0.push_back(tmp);
	}
      } // loop : d(j)
	// -- C_FILLMATRIX diag/offdiagonal for grand fathers --
      for (int d = 1; d < begdom; d++) {
	const int j = btree->selfIndex(d);
	task_name += to_string(d) + " ";
	C_task_seq* tmp = 
	  new C_task_seq(C_FILLMATRIX,
			   _null_name,      //     task_name_cstr,
			 (-1),           // mutex_id
			 (-1),           // TASK_PARALLEL,
			 (-1),           // num_threads,
			 (-1),           //        level
			 (-1),           //         phase
			 &tasks_DFillSym[j],
			 0,
			 tasks_DFillSym[j].size(),
			 0L // nops;
			 );
	queue_tmp0.push_back(tmp);
      } // loop : d(j)
      if (queue_tmp0.size() > 0) {              // only the last block
	copytask_list2seq(queue_dynamic[level_last],
			  queue_tmp0,
			  queue_null_,
  //			  queue_dummy_,
			  task_name,
			  C_FILLMATRIX, // C_FILLMATRIX1,
			  0,               // shared by all threads
			  TASK_PARALLEL,
			  num_threads,
			  level_last,
			  2);

  // erase temporary C_task_seq whose elements are copied to queue_dynamic[]
	for (list<C_task_seq*>::iterator it = queue_tmp0.begin();
	     it != queue_tmp0.end(); ++it) {
	  delete (*it);
	  (*it) = NULL;
	} 
      }
    }
  }
#endif
  {
    for (int p = 0; p < num_threads; p++) {
      for(list<C_task_seq *>::const_iterator it = queue_dynamic[level_last].begin();
	  it != queue_dynamic[level_last].end(); ++it) {
	_queue_static[p].push_back(*it);
      }
    }
  }
    // for dense subdomains
  for (int level = (level_last - 1); level >=0; level--) {
    const int begdom = 1U << level;
    const int enddom = begdom * 2;
    const int num_leaves = 1 << (level + 1); // 1 << level
    if (level < (level_last - 1)) {
    // better management by mixed Assignment of blocks of DTRSM/DGEMM 
    // -- DTRSM offdiagonal --
      {
	string task_name = "q offdiag " + to_string(level + 1) + " : ";
	queue_tmp0.clear();
    // off-diagonal parallel block between children and grandfather, ancenstors
	for (int d = begdom; d < enddom; d++) {
	  const int j = btree->selfIndex(d);

	  for (int ic = 0; ic < 2; ic++) {
	    const int jc = children[j][ic];
	    const int end = tasks_DTRSMScale[jc].size();
	    int begin = tasks_DTRSMScale[jc].size();
	    if (!isDividedDTRSM[jc] && (end > 0) && !isMergedDTRSM[j]) {
	      begin = IndexTask(tasks_DTRSMScale[jc],
				  tasks_DTRSMScale[jc][0]->parallel_max);
	    }
	    if ((end - begin) > 0) { 
	      task_name += to_string(jc) + " ";
	      C_task_seq* tmp = 
		new C_task_seq(C_DTRSM,
			       _null_name,   //	     task_name_cstr,
			       (-1),           // mutex_id
			       (-1),           // TASK_PARALLEL,
			       (-1),           // num_threads,
			       (-1),           //        level
			       (-1),           //         phase
			       &tasks_DTRSMScale[jc],
			       begin,
			       end,
			       0L // nops;
			       );
	      queue_tmp0.push_back(tmp);
	    } //  if ((end - begin) > 0) 
	  } // loop : ic
	} // loop : j
	if (queue_tmp0.size() > 0) {              // only the last block
	  copytask_list2seq(queue_dynamic[level],
			    queue_tmp0,
			    queue_null_,
 //			    queue_dummy_,
			    task_name,
			    C_DTRSM1,
			    0,               // shared by all threads
			    TASK_PARALLEL,
			    num_threads,
			    (level + 1),
			    (level_last - 1 - level) * 6 + 1);
 // erase temporary C_task_seq whose elements are copied to queue_dynamic[]
	  for (list<C_task_seq*>::iterator it = queue_tmp0.begin();
	       it != queue_tmp0.end(); ++it) {
	    delete (*it);
	    (*it) = NULL;
	  } 
	}
      }
    // dependency
    // -- D_GEMM offdiagonal --
      {
	string task_name = "r " + to_string(level + 1) + " : ";
	queue_tmp0.clear();
	for (int d = begdom; d < enddom; d++) {
	  const int j = btree->selfIndex(d);
	  for (int ic = 0; ic < 2; ic++) {
	    const int jc = children[j][ic];
	    const int end = tasks_DGEMM[jc].size();
	    if (end > 0) {
	      int begin;
	      const int itmp = tasks_DGEMM[jc][0]->parallel_max;
	      if (!isMergedDTRSM[j]) {
		begin = itmp;
	      }
	      else {
		begin = ((itmp == end) ? end : 
			 (itmp + tasks_DGEMM[jc][itmp]->parallel_max));
	      }
	      if ((end - begin) > 0) {
		task_name += to_string(jc) + " ";
		C_task_seq* tmp = 
		  new C_task_seq(C_DGEMM,
				 _null_name, // task_name_cstr,
				 (-1), // mutex_id
				 (-1), // TASK_PARALLEL,
				 (-1), // num_threads,
				 (-1), // (level + 1),
				 (-1), // (level_last - 1 - level) * 5 + 2,
				 &tasks_DGEMM[jc],
				 begin, 
				 end,
				 0L// nops);
				 );
		queue_tmp0.push_back(tmp);
	      }
	    } // if (tasks_DGEMM[jc].size() > 0)
	  } // loop : ic
	}  // loop j
	if (queue_tmp0.size() > 0) {              // only the last part
	  copytask_list2seq(queue_dynamic[level],
			    queue_tmp0,
			    queue_null_,
 //			    queue_dummy_,
			    task_name,
			    C_DGEMM1,
			    0,                // shared by all threads
			    TASK_PARALLEL,
			    num_threads,
			    (level + 1),
			    (level_last - 1 - level) * 6 + 2);
     // erase temporary C_task_seq whose elements are copied to queue_dynamic[]
	  for (list<C_task_seq*>::iterator it = queue_tmp0.begin();
	       it != queue_tmp0.end(); ++it) {
	    delete (*it);
	    (*it) = NULL;
	  } 
	} //  if (queue_tmp0.size() > 0)
      }  // scope of task_name
    } //  if (level < (level_last - 1)) 
    // -- DEALLACATE
    {
      string task_name = "X " + to_string(level + 1) + " : ";
      queue_tmp0.clear();
      C_task_seq* tmp = 
	new C_task_seq(C_DEALLOCATE,
		       _null_name,
		       (-1), // mutex_id
		       (-1), //TASK_PARALLEL,
		       (-1),//num_threads,
		       (-1),//(level + 1),
		       (-1),//(level_last - 1 - level) * 5 + 3,
		       &tasks_deallocLower[level + 1],
		       0,
		       tasks_deallocLower[level + 1].size(),
		       0L// nops);
		       );
      queue_tmp0.push_back(tmp);
      copytask_list2seq(queue_dynamic[level],
			queue_tmp0,
			queue_null_,
//			queue_dummy_,
			task_name,
			C_DEALLOCATE1,
			0,               // shared by all threads
			TASK_PARALLEL,
			num_threads,
			(level + 1),
			(level_last - 1 - level) * 6 + 3);
 // erase temporary C_task_seq whose elements are copied to queue_dynamic[]
      for (list<C_task_seq*>::iterator it = queue_tmp0.begin();
	   it != queue_tmp0.end(); ++it) {
	delete (*it);
	(*it) = NULL;
      } 
    }

    // -- DSUB offdiagonal --
    {
      string task_name = "s " + to_string(level + 1) + " : ";
      queue_tmp0.clear();
      for (int d = begdom; d < enddom; d++) {
	const int j = btree->selfIndex(d);
	if (!isMergedDTRSM[j]) {
	  const int end = tasks_Dsub[level + 1][j].size();
	  if (end > 0) {  // to skip non allocated tasks_Dsub[level + 1][j][0]
	    const int begin = tasks_Dsub[level + 1][j][0]->parallel_max;
	    if ((end - begin) > 0) {
	      task_name += to_string(j) + " ";
	      C_task_seq* tmp = 
		new C_task_seq(C_DSUB,
			       _null_name,
			       (-1), // mutex_id
			       (-1), //TASK_PARALLEL,
			       (-1),//num_threads,
			       (-1),//(level + 1),
			       (-1),//(level_last - 1 - level) * 5 + 3,
			       &tasks_Dsub[level + 1][j],
			       begin,
			       end,
			       0L// nops);
			       );
	      queue_tmp0.push_back(tmp);
	    } // if ((end - begin) > 0)
	  }   // if (end > 0)
	}     // if (!isMergedDTRSM[j])
      }  // loop : d / j
      for (vector<int>::const_iterator it = all_fathersIndex[level].begin();
	   it != all_fathersIndex[level].end(); ++it) {
	const int size = tasks_Dsub[level + 1][*it].size();
	if (size > 0) {
	  task_name += to_string(*it) + " ";
	  C_task_seq* tmp = 
	    new C_task_seq(C_DSUB,
			   _null_name,
			   (-1), // mutex_id
			   (-1), //TASK_PARALLEL,
			   (-1), //num_threads,
			   (-1), //(level + 1),
			   (-1), //(level_last - 1 - level) * 5 + 3,
			   &tasks_Dsub[level + 1][*it],
			   0,
			   size,
			   0L// nops);
			   );
	  queue_tmp0.push_back(tmp);
	}
      } // loop : it
      if (queue_tmp0.size() > 0) {              // only the last part
	copytask_list2seq(queue_dynamic[level],
			  queue_tmp0,
			  queue_null_,
  //			  queue_dummy_,
			  task_name,
			  C_DSUB1,
			  0,               // shared by all threads
			  TASK_PARALLEL,
			  num_threads,
			  (level + 1),
			  (level_last - 1 - level) * 6 + 4);
	// erase temporary C_task_seq whose elements are copied to 
	// queue_dynamic[]
	for (list<C_task_seq*>::iterator it = queue_tmp0.begin();
	     it != queue_tmp0.end(); ++it) {
	  delete (*it);
	  (*it) = NULL;
	} 
      }
    }
#if 0
    // -- DEALLACATE
    {
      string task_name = "X " + to_string(level + 1) + " : ";
      queue_tmp0.clear();
      C_task_seq* tmp = 
	new C_task_seq(C_DEALLOCATE,
		       _null_name,
		       (-1), // mutex_id
		       (-1), //TASK_PARALLEL,
		       (-1),//num_threads,
		       (-1),//(level + 1),
		       (-1),//(level_last - 1 - level) * 5 + 3,
		       &tasks_deallocLower[level + 1],
		       0,
		       tasks_deallocLower[level + 1].size(),
		       0L// nops);
		       );
      queue_tmp0.push_back(tmp);
      copytask_list2seq(queue_dynamic2[level],
			queue_tmp0,
			queue_null_,
//			queue_dummy_,
			task_name,
			C_DEALLOCATE1,
			0,               // shared by all threads
			TASK_PARALLEL,
			num_threads,
			(level + 1),
			(level_last - 1 - level) * 6 + 4);
 // erase temporary C_task_seq whose elements are copied to queue_dynamic[]
      for (list<C_task_seq*>::iterator it = queue_tmp0.begin();
	   it != queue_tmp0.end(); ++it) {
	delete (*it);
	(*it) = NULL;
      } 
    }
#endif
    // -- DEALLACATE
    {
      string task_name = "x " + to_string(level + 1)  + " : ";
      queue_tmp0.clear();
      const int begdom1 = 1U << (level + 1);
      const int enddom1 = begdom1 * 2;
      for (int d = begdom1; d < enddom1; d++) {
	const int j = btree->selfIndex(d);
	C_task_seq* tmp = 
	  new C_task_seq(C_DEALLOCATE,
			 _null_name,
			 (-1), // mutex_id
			 (-1), //TASK_PARALLEL,
			 (-1),//num_threads,
			 (-1),//(level + 1),
			 (-1),//(level_last - 1 - level) * 5 + 3,
			 &tasks_deallocLocalSchur[j],
			 0,
			 tasks_deallocLocalSchur[j].size(),
			 0L// nops);
			 );
	queue_tmp0.push_back(tmp);
      }
      copytask_list2seq(queue_dynamic2[level],
			queue_tmp0,
			queue_null_,
//			queue_dummy_,
			task_name,
			C_DEALLOCATE1,
			0,               // shared by all threads
			TASK_PARALLEL,
			num_threads,
			(level + 1),
			(level_last - 1 - level) * 6 + 5);
 // erase temporary C_task_seq whose elements are copied to queue_dynamic[]
      for (list<C_task_seq*>::iterator it = queue_tmp0.begin();
	   it != queue_tmp0.end(); ++it) {
	delete (*it);
	(*it) = NULL;
      } 
    }

    // level == 0 <==> num_leaves = 2
    if (num_leaves <= num_threads) { 
      // trick to merge DFullLDLt into previous stage of DTRSM to follow
      // the critical path
      long *nops_sum = new long[num_threads];
      for (int i = 0; i < num_threads; i++) {
	nops_sum[i] = 0L;
      }
      vector<int> starts(nb_doms_dense, 0);
      vector<int> starts_dtrsm(nb_doms_dense, 0);
      vector<int> starts_sub(nb_doms_dense, 0);

      long nops_block_total = 0L;
      // -- STARTING BLOCK : DTSRM + DGEMM + DSUB + DFULLLDLT --
      list<vector<C_task *> *> diag_starting;
      list<string> task_name_strs;
      for (int p = 0; p < num_threads; p++) {
	queue_lists[p].clear();
      }
      for (int d = begdom; d < enddom; d++) {
	const int j = btree->selfIndex(d);
#ifdef DEBUG_PREPARE_THREAD
	cout << "+ level = " << (level + 1) << " j = " << j << " children = " 
	     << children[j][0] << " " << children[j][1] << " : " 
	   << "# DFullLDLt @ " << level << " : " << j << " / " 
	     << tasks_DFullLDLt[j].size() << endl;
#endif
	// 3 * 4 * 5 / 6 + 1 = 11
	// 3 * 4 * 5 / 6 + 2 * 3 * 4 / 6 + 3 = 17
	//	if (tasks_DFullLDLt[j].size() > 4) { // with parallelizing
	if(nrow_DFullLDLt[j] > 2) {
	  vector<C_task *> *tmp = new vector<C_task *>;
	  tmp->reserve(4 + 2 * (_isSym ? 3 : 6));  // 3 = 2 + 1
	  for (int ic = 0; ic < 2; ic++) {
	    const int jj = children[j][ic];
	    if ((nrow_DFullLDLt[jj] <= 2) && !isMergedDTRSM[jj]) {
// the case DFullLDLt and DTRSMScale are concatenated can be skipped
	      const int itmp = tasks_DTRSMScale[jj][0]->atomic_size;
	      for (int k = 0; k < itmp; k++) {
		tmp->push_back(tasks_DTRSMScale[jj][k]);
	      }
	      if (_isSym) {
		starts_dtrsm[jj] = tasks_DTRSMScale[jj][0]->atomic_size;
	      }
	      else {
		for (int k = 0; k < tasks_DTRSMScale[jj][itmp]->atomic_size; 
		     k++) {
		  tmp->push_back(tasks_DTRSMScale[jj][itmp + k]);
		}
		starts_dtrsm[jj] = (itmp + 
				      tasks_DTRSMScale[jj][itmp]->atomic_size);
	      } // if (_isSym)
	    }   // if ((nrow_DFullLDLt[jj] <= 2) && !isMergedDTRSM[jj]) 
	    tmp->push_back(tasks_DGEMM[jj][0]);
	    starts[jj] = 1;
	  }     // loop : ic
	  tmp->push_back(tasks_Dsub[level + 1][j][0]); // null for many cases 
	  tmp->push_back(tasks_DFullLDLt[j][0]);
	  diag_starting.push_back(tmp);

	  string task_name = "p " + to_string(level + 1) + " " + to_string(j);
	  task_name_strs.push_back(task_name);
	  starts_sub[j] = 1;
	} //  if(nrow_DFullLDLt[j] > 2)
      }   // loop : d
      list<string>::const_iterator kt = task_name_strs.begin();
      {
	int itmp = 0;
	for (list<vector<C_task *> *>::const_iterator it = diag_starting.begin();
	     it != diag_starting.end(); ++it, ++kt, itmp++){
	  long nops = 0L;
	  for (vector<C_task *>::const_iterator jt = (*it)->begin();
	       jt != (*it)->end(); ++jt) {
	    nops += *((*jt)->ops_complexity);
	  }
	  C_task_seq* tmp = 
	    new C_task_seq(C_DIAG_START,
			   (*kt),
			   TASK_SINGLE,
			   (-1), // mutex_id
			   1,
			   (level + 1),
			   (level_last - 1 - level) * 6 + 1,
			   (*it),
			   0,
			   (*it)->size(),
			   nops);
	  queue_lists[itmp % num_threads].push_back(tmp);
	  nops_sum[itmp % num_threads] = nops;
	  nops_block_total += nops;
#ifdef DEBUG_PREPARE_THREAD
	  cout << "ll = " << itmp % num_threads 
	       << " critical path " << *kt << " nops = " << nops << endl;
#endif
	} // loop : it, kt, itmp
      } // scope of itmp
// -- DTRSM -- diangonal
      task_assign_diag1(_queue_static,
			tasks_DTRSMScale,
			queue_null_,
//			queue_dummy_,
			nrow_DFullLDLt,
			isMergedDTRSM,
			"q",
			C_DTRSM,
			1, // int level_id,
			num_threads,
			level_last,
			level,
			begdom,
			enddom,
			nops_sum,
			starts_dtrsm,
			queue_lists,
			false,
			btree,
			children);
      // erase STL container defined in line 1499
      for (list<vector<C_task *> *>::iterator it = diag_starting.begin();
	   it != diag_starting.end(); ++it) {
	delete (*it);
	(*it) = NULL;
      }   
      { // begin : scope queue_tmp2
	list<C_task_seq *> queue_tmp2;
	queue_tmp2.clear();	
	for (int d = begdom; d < enddom; d++) {
	  const int j = btree->selfIndex(d);
	  for (int ic = 0; ic < 2; ic++) {
	    const int jc = children[j][ic];
	    long nops = 0L;
	    if (isDividedDTRSM[jc]) {
	      const int end = tasks_DTRSMScale[jc].size();
	      for (int i = 0; i < end; i++) {
		nops += *(tasks_DTRSMScale[jc][i]->ops_complexity);
	      }
	      C_task_seq* tmp = 
		new C_task_seq(C_DTRSM,
			       _null_name,
			       (-1), // mutex_id
			       (-1), // TASK_SINGLE,
			       (-1), // 1,
			       (-1), // (level + 1),
			       (-1), // (level_last - 1 - level) * 5 + 1,
			       &tasks_DTRSMScale[jc],
			       0,
			       end,
			       nops);
	      queue_tmp2.push_back(tmp);
	    } //  if (isDividedDTRSM[jc])
	  } // loop : ic
	} // loop : d
	{ // begin : scope task_name
	  string task_name = "q Offdiag " + to_string(level + 1);
	  if (queue_tmp2.size() > 0) {              // only the last block
	    copytask_list2seq(queue_dynamic1[level],
			      queue_tmp2,
			      queue_null_,
     //			      queue_dummy_,
			      task_name,
			      C_DTRSM,
			      0,               // shared by all threads
			      TASK_PARALLEL,
			      num_threads,
			      (level + 1),
			      (level_last - 1 - level) * 6 + 1);
	    for (list<C_task_seq*>::iterator it = queue_tmp2.begin();
		 it != queue_tmp2.end(); ++it) {
	      delete (*it);
	      (*it) = NULL;
	    } // loop : it
	  } // if (queue_tmp2.size() > 0) 
	}  // end : scope task_name 
      }   // end : scope queue_tmp2
      if (level >= 0) {
	for (int p = 0; p < num_threads; p++) {
	  for(list<C_task_seq *>::const_iterator it = queue_dynamic1[level].begin();
	      it != queue_dynamic1[level].end(); ++it) {
	      _queue_static[p].push_back(*it);
	  }
	} // loop : p
// the case queue_dynamic consists of only C_task_seq == 'x', copied here
	for (int p = 0; p < num_threads; p++) {
	  for(list<C_task_seq *>::const_iterator it = queue_dynamic2[level + 1].begin();
	      it != queue_dynamic2[level + 1].end(); ++it) {
	    _queue_static[p].push_back(*it);
	  }
	} // loop : p
      }

// -- DGEMM -- diagonal
      vector<int> null_idx;
      task_assign_diag1(_queue_static,
			tasks_DGEMM,
			queue_null_,
//			queue_dummy_,
			null_idx, // not used
			isMergedDTRSM,
			"r",
			C_DGEMM,
			2, // int level_id,
			num_threads,
			level_last,
			level,
			begdom,
			enddom,
			nops_sum,
			starts,
			queue_lists,
			true,
			btree,
			children);

// -- DSUB -- diagonal
      task_assign_diag2(_queue_static,
			tasks_Dsub[level + 1],
			queue_null_,
//			queue_dummy_,
			isMergedDTRSM,
			//tasks_DFullLDLt,
			"s",
			C_DSUB,
			3, // int level_id,
			num_threads,
			level_last,
			level,
			begdom,
			enddom,
			nops_sum,
			starts,
			starts_sub,
			queue_lists,
			true,
			btree,
			children);

      if (level == 0) {
// the case queue_dynamic consists of only C_task_seq == 'x', copied here
	for (int p = 0; p < num_threads; p++) {
	  for(list<C_task_seq *>::const_iterator it = queue_dynamic2[level].begin();
	      it != queue_dynamic2[level].end(); ++it) {
	    _queue_static[p].push_back(*it);
	  }
	} // loop : p
      }
      // -- C_DFULL -- diagonal
      queue_tmp0.clear(); 
      queue_tmp1.clear(); 
      for (int p = 0; p < num_threads; p++) {
	queue_lists[p].clear();
      }
      for (int d = begdom; d < enddom; d++) {
	const int j = btree->selfIndex(d);
	long nops = 0L;
	for (int i = starts_sub[j]; i < tasks_DFullLDLt[j].size(); i++) {
	  nops +=  *(tasks_DFullLDLt[j][i]->ops_complexity);
	}
	string task_name = "t& " + to_string(level) + " : " + to_string(j);
	C_task_seq* tmp = new C_task_seq(C_DFULL,
					 task_name,
					 0, // mutex_id
					 TASK_SINGLE,
					 1, 
					 level,
					 (level_last - level) * 6,
					 &tasks_DFullLDLt[j],
					 starts_sub[j],
					 tasks_DFullLDLt[j].size(),
					 nops);
	if (starts_sub[j] > 0) {
	  queue_tmp0.push_back(tmp);
	}
	else {
	  queue_tmp1.push_back(tmp);
	}
      }
      const int size_task0 = queue_tmp0.size();
      const int size_threads = num_threads - queue_tmp1.size();
#ifdef DEBUG_PREPARE_THREAD
      cout << "*** C_DFULL : parallel = " << size_task0 << " : serial " 
	   << queue_tmp1.size() << " / ";
#endif
      if (size_task0 > 0 && size_threads > 0) {
	queue_tmp0.sort(C_task_seq_complexity_greater);
	int mm0 = size_threads % size_task0;
	const int mm = (size_threads / size_task0) + (mm0 != 0);
	if (mm0 == 0) {
	  mm0 = size_task0;
	}
	int jj = 0;
	for (list<C_task_seq *>::const_iterator it = queue_tmp0.begin();
	     it != queue_tmp0.end(); ++it, jj++) {
	  
	  const int multplcty = (jj < mm0) ? mm : (mm - 1);
#ifdef DEBUG_PREPARE_THREAD
	  cout << (*it)->task_id << " : " << (*it)->task_name 
	       << " : < " << multplcty << " > ";
#endif
	  const int parallel_single = 
	    (multplcty == 1) ? TASK_SINGLE : TASK_PARALLEL;
	  // copy raw task into a task sequence whose length is one
	  list<C_task_seq *> queue_tmp2, queue_tmp3;
	  queue_tmp2.push_back(*it);
	  copytask_list2seq(queue_tmp3,
			    queue_tmp2,
			    queue_null_,
    //			    queue_dummy_,
			    (*it)->task_name,
			    (*it)->task_id,
			    jj,
			    parallel_single,
			    multplcty,
			    (*it)->level,
			    (*it)->phase);
	  // threads share one task (*it)
	  for (int m = 0; m < mm; m++)  {
	    const int n = jj + m * size_task0;
	    if (n < size_threads) {
	      _queue_static[n % num_threads].push_back(queue_tmp3.front());
	      //  _queue_static[n % num_threads].push_back(seq_tmp);
#ifdef DEBUG_PREPARE_THREAD
	      cout << "[ " << jj << " @ " << n << " : " 
		   << n % num_threads << " ]";
#endif
	    }
	  }
	} // loop : j
#ifdef DEBUG_PREPARE_THREAD
	cout << "/ ";
#endif
	jj = size_threads;
	for (list<C_task_seq *>::const_iterator it = queue_tmp1.begin();
	     it != queue_tmp1.end(); ++it, jj++) {
	  // copy for each thread
	    list<C_task_seq *> queue_tmp2, queue_tmp3;
	    queue_tmp2.push_back(*it);
	    copytask_list2seq(queue_tmp3,
			      queue_tmp2,
			      queue_null_,
     //			      queue_dummy_,
			      (*it)->task_name,
			      (*it)->task_id,
			      0,
			      TASK_SINGLE,
			      1,
			      (*it)->level,
			      (*it)->phase);

	    _queue_static[jj % num_threads].push_back(queue_tmp3.front());
#ifdef DEBUG_PREPARE_THREAD
	  cout << jj << " ; " << jj % num_threads << " / ";
#endif
	}
#ifdef DEBUG_PREPARE_THREAD
	cout << endl;
#endif
	// delete original C_task_seq after copying to all threads
	for (list<C_task_seq*>::iterator it = queue_tmp0.begin();
	     it != queue_tmp0.end(); ++it) {
	  delete (*it);
	  (*it) = NULL;
	}
	for (list<C_task_seq*>::iterator it = queue_tmp1.begin();
	       it != queue_tmp1.end(); ++it) {
	  delete (*it);
	  (*it) = NULL;
	}  
      }   // if (size_task0 > 0)
      else { // (size_task0 == 0) || (size_threads <= 0)
#ifdef DEBUG_PREPARE_THREAD
	cout << endl;
#endif
	{
	  int itmp = 0;
	  queue_tmp1.sort(C_task_seq_complexity_greater);
	  for (list<C_task_seq*>::const_iterator it = queue_tmp1.begin();
	       it != queue_tmp1.end(); ++it, itmp++) {
	    (*it)->parallel_single = TASK_SINGLE;
	    queue_lists[itmp % num_threads].push_back(*it);
	  } // loop : it, itmp
	}
	{
	  string task_name = "t " + to_string(level);
	  for (int p = 0; p < num_threads; p++) {
	    copytask_list2seq(_queue_static[p],
			      queue_lists[p],
			      queue_null_,
     //			      queue_dummy_,
			      task_name,
			      C_DFULL,
			      p,
			      TASK_SINGLE,
			      1,
			      level,
			      (level_last - level) * 6);
	  }  // loop : p
	  // (queue_tmp0.size() == 0) ? for safty of non-memory-leak
	  for (list<C_task_seq*>::iterator it = queue_tmp0.begin();
	       it != queue_tmp0.end(); ++it) {
	    delete (*it);
	    (*it) = NULL;
	  }
	  for (list<C_task_seq*>::iterator it = queue_tmp1.begin();
	       it != queue_tmp1.end(); ++it) {
	    delete (*it);
	    (*it) = NULL;
	  }  
	}
      }  // if (size_task > 0)
      if (level > 0) {
	for (int p = 0; p < num_threads; p++) {
	  for(list<C_task_seq *>::const_iterator it = queue_dynamic[level].begin();
	      it != queue_dynamic[level].end(); ++it) {
	    _queue_static[p].push_back(*it);
	  }
	}
      }
      delete [] nops_sum;
    } //  
    else { // if (num_leaves <= num_threads) 
      if (level < (level_last - 1)) {
      // -- DTRSM -- diagonal
	queue_tmp0.clear();
	for (int p = 0; p < num_threads; p++) {
	  queue_lists[p].clear();
	}
	list<C_task_seq *> queue_tmp2;
	queue_tmp2.clear();
	for (int d = begdom; d < enddom; d++) {
	  const int j = btree->selfIndex(d);
#ifdef DEBUG_PREPARE_THREAD
	  cout << "* level = " << (level + 1) << " j = " << j << " children = " 
	       << children[j][0] << " " << children[j][1] << " " << endl;
#endif
	  for (int ic = 0; ic < 2; ic++) {
	    const int jc = children[j][ic];
	    if (tasks_DTRSMScale[jc].size() > 0) {
	      long nops = 0L;
	      if (isDividedDTRSM[jc]) {
		const int end = tasks_DTRSMScale[jc].size();
		for (int i = 0; i < end; i++) {
		  nops += *(tasks_DTRSMScale[jc][i]->ops_complexity);
		}
		C_task_seq* tmp = 
		  new C_task_seq(C_DTRSM,
				 _null_name,
				 (-1), // mutex_id
				 (-1), // TASK_SINGLE,
				 (-1), // 1,
				 (-1), // (level + 1),
				 (-1), // (level_last - 1 - level) * 5 + 1,
				 &tasks_DTRSMScale[jc],
				 0,
				 end,
				 nops);
		queue_tmp2.push_back(tmp);
	      } //  if (isDividedDTRSM[jc])
	      else {
		int end;
		if (!isMergedDTRSM[j]) {  // looking the father node
		  end = IndexTask(tasks_DTRSMScale[jc],
				  tasks_DTRSMScale[jc][0]->parallel_max);
		}
		else {
		  end = tasks_DTRSMScale[jc].size();
		}
		for (int i = 0; i < end; i++) {
		  nops += *(tasks_DTRSMScale[jc][i]->ops_complexity);
		}
		C_task_seq* tmp = 
		  new C_task_seq(C_DTRSM,
				 _null_name,
				 (-1), // mutex_id
				 (-1), // TASK_SINGLE,
				 (-1), // 1,
				 (-1), // (level + 1),
				 (-1), // (level_last - 1 - level) * 5 + 1,
				 &tasks_DTRSMScale[jc],
				 0,
				 end,
				 nops);
		queue_tmp0.push_back(tmp);
	      } // if (isDividedDTRSM[jc])
	    } // if (tasks_DTRSMScale[jc].size() > 0)
	  } // loop : ic
	}  // loop : d
	queue_tmp0.sort(C_task_seq_complexity_smaller);
	{ // begin : scope itmp
	  int itmp = 0;
	  for (list<C_task_seq*>::const_iterator it = queue_tmp0.begin();
	       it != queue_tmp0.end(); ++it, itmp++) {
	    queue_lists[itmp % num_threads].push_back(*it);
	  } // loop : it
	} // end : scope itmp
	{ // begin : scope task_name
	  string task_name = "q diag " + to_string(level + 1);
	  for (int p = 0; p < num_threads; p++) {
	    copytask_list2seq(_queue_static[p],
			      queue_lists[p],
			      queue_null_,
     //			      queue_dummy_,
			      task_name,
			      C_DTRSM,
			      p,
			      TASK_SINGLE,
			      1,
			      (level + 1),
			      (level_last - 1 - level) * 6 + 1);
	  }
	  // erase temporary C_task_seq whose elements are copied to 
	  // queue_static[]
	  for (list<C_task_seq*>::iterator it = queue_tmp0.begin();
	       it != queue_tmp0.end(); ++it) {
	    delete (*it);
	    (*it) = NULL;
	  } 
	}  // end : scope task_name 
	{ // begin : scope task_name
	  string task_name = "q Offdiag " + to_string(level + 1);
	  if (queue_tmp2.size() > 0) {              // only the last block
	    copytask_list2seq(queue_dynamic1[level],
			      queue_tmp2,
			      queue_null_,
     //			      queue_dummy_,
			      task_name,
			      C_DTRSM,
			      0,               // shared by all threads
			      TASK_PARALLEL,
			      num_threads,
			      (level + 1),
			      (level_last - 1 - level) * 6 + 1);
    // erase temporary C_task_seq whose elements are copied to queue_dynamic[]
	    for (list<C_task_seq*>::iterator it = queue_tmp2.begin();
		 it != queue_tmp2.end(); ++it) {
	      delete (*it);
	      (*it) = NULL;
	    }
	  } // if (queue_tmp2.size() > 0) 
	}  // end : scope task_name 
	if (level >= 0) {
	  for (int p = 0; p < num_threads; p++) {
	    for(list<C_task_seq *>::const_iterator it = queue_dynamic1[level].begin();
		it != queue_dynamic1[level].end(); ++it) {
	      _queue_static[p].push_back(*it);
	    }
	  } // loop : p
	  // C_task_seq == 'x', copied here
	  for (int p = 0; p < num_threads; p++) {
	    for(list<C_task_seq *>::const_iterator it = queue_dynamic2[level + 1].begin();
		it != queue_dynamic2[level + 1].end(); ++it) {
	      _queue_static[p].push_back(*it);
	    }
	  } // loop : p
	} //  if (level >= 0) 
      // -- DGEMM -- diagonal
	queue_tmp0.clear();
	for (int p = 0; p < num_threads; p++) {
	  queue_lists[p].clear();
	}
	for (int d = begdom; d < enddom; d++) {
	  const int j = btree->selfIndex(d);
	  for (int ic = 0; ic < 2; ic++) {
	    const int jc = children[j][ic];
	    const int tasks_size = tasks_DGEMM[jc].size();
	    if (tasks_size > 0) {
	      long nops = 0L;
	      int end;
	      int itmp = tasks_DGEMM[jc][0]->parallel_max;
	      //	      if (itmp <= (_isSym ? 3 : 4)) {
#if 1
	      if (!isMergedDTRSM[j]) {
		end = itmp;
	      }
	      else {
		end = ((itmp == tasks_size) ? tasks_size : 
		       itmp + tasks_DGEMM[jc][itmp]->parallel_max);
	      }
#else
	      end = ((itmp == tasks_size) ? tasks_size : 
		     itmp + tasks_DGEMM[jc][itmp]->parallel_max);
#endif
	      for (int i = 0; i < end; i++) {
		nops += *(tasks_DGEMM[jc][i]->ops_complexity);
	      }
	      C_task_seq* tmp = 
		new C_task_seq(C_DGEMM,
			       _null_name,
			       (-1), // mutex_id
			       (-1), // TASK_SINGLE,
			       (-1), //1,
			       (-1), //(level + 1),
			       (-1), //(level_last - 1 - level) * 5 + 2,
			       &tasks_DGEMM[jc],
			       0,
			       end,
			       nops);
	      queue_tmp0.push_back(tmp);
	    }  //    if (tasks_DGEMM[jc].size() > 0) 
	  } //  loop : ic
	} // loop : d
	queue_tmp0.sort(C_task_seq_complexity_greater);
	{
	  int itmp = 0;
	  for (list<C_task_seq*>::const_iterator it = queue_tmp0.begin();
	       it != queue_tmp0.end(); ++it, itmp++) {
	    queue_lists[itmp % num_threads].push_back(*it);
	  } // loop : it
	}
	{
	  string task_name = "r " + to_string(level + 1);
	  for (int p = 0; p < num_threads; p++) {
	    copytask_list2seq(_queue_static[p],
			      queue_lists[p],
			      queue_null_,
     //			      queue_dummy_,
			      task_name,
			      C_DGEMM,
			      p,
			      TASK_SINGLE,
			      1,
			      (level + 1),
			      (level_last - 1 - level) * 6 + 2);
	  }
	  // erase temporary C_task_seq whose elements are copied to 
	  // queue_static[]
	  for (list<C_task_seq*>::iterator it = queue_tmp0.begin();
	       it != queue_tmp0.end(); ++it) {
	    delete (*it);
	    (*it) = NULL;
	  } 
	}
      } // if (level == (level_last - 1))
      // -- DSUB -- diagonal
      queue_tmp0.clear();
      for (int p = 0; p < num_threads; p++) {
	queue_lists[p].clear();
      } 
      for (int d = begdom; d < enddom; d++) {
	const int j = btree->selfIndex(d);
	long nops = 0L;
	if (tasks_Dsub[level + 1][j].size() > 0) {
	  int end;
	  //	  if (tasks_DFullLDLt[j].size() <= 4) {
	  //	  if(nrow_DFullLDLt[j] <= 2) {
	  if (!isMergedDTRSM[j]) {
	    end = tasks_Dsub[level + 1][j][0]->parallel_max;
	  }
	  else {
	    end = tasks_Dsub[level + 1][j].size();
	  }
	  for (int i = 0; i < end; i++) {
	    nops += *(tasks_Dsub[level + 1][j][i]->ops_complexity);
	  }
	  C_task_seq* tmp = 
	    new C_task_seq(C_DSUB,
			   _null_name,
			   (-1), // mutex_id
			   (-1), // TASK_SINGLE,
			   (-1), //1,
			   (-1), //(level + 1),
			   (-1), //(level_last - 1 - level) * 5 + 3,
			   &tasks_Dsub[level + 1][j],
			   0,
			   end,
			   //	   tasks_Dsub[level + 1][j][0]->parallel_max,
			   nops);
	  queue_tmp0.push_back(tmp);
	}
      } // loop : d
      queue_tmp0.sort(C_task_seq_complexity_greater);
      {
	int itmp = 0;
	for (list<C_task_seq*>::const_iterator it = queue_tmp0.begin();
	     it != queue_tmp0.end(); ++it, itmp++) {
	  queue_lists[itmp % num_threads].push_back(*it);
	} // loop : it
      }
      {  // begin : scope for task_name;
	string task_name = "s " + to_string(level + 1);
	for (int p = 0; p < num_threads; p++) {
	  copytask_list2seq(_queue_static[p],
			    queue_lists[p],
			    queue_null_,
   //			    queue_dummy_,
			    task_name, 
			    C_DSUB,
			    p,
			    TASK_SINGLE,
			    1,
			    (level + 1),
			    (level_last - 1 - level) * 6 + 3);
	}
	// erase temporary C_task_seq whose elements are copied to 
	// queue_static[]
	for (list<C_task_seq*>::iterator it = queue_tmp0.begin();
	     it != queue_tmp0.end(); ++it) {
	  delete (*it);
	  (*it) = NULL;
	} 
      } // end : scope for task_name;
      if (level == 0) {
// the case queue_dynamic consists of only C_task_seq == 'x', copied here
	for (int p = 0; p < num_threads; p++) {
	  for(list<C_task_seq *>::const_iterator it = queue_dynamic2[level].begin();
	      it != queue_dynamic2[level].end(); ++it) {
	    _queue_static[p].push_back(*it);
	  }
	} // loop : p
      }
      // -- C_DFULL --
      queue_tmp0.clear();
      for (int p = 0; p < num_threads; p++) {
	queue_lists[p].clear();
      }
      for (int d = begdom; d < enddom; d++) {
	const int j = btree->selfIndex(d);
	long nops = 0L;
	for (int i = 0; i < tasks_DFullLDLt[j].size(); i++) {
	  nops += *(tasks_DFullLDLt[j][i]->ops_complexity);
	}
	string task_name = "t " + to_string(level) + " : " + to_string(j);
	C_task_seq* tmp = new C_task_seq(C_DFULL,
					 task_name,
					 (-1), // mutex_id
					 TASK_SINGLE,
					 1,
					 level,
					 (level_last - level) * 6,
					 &tasks_DFullLDLt[j],
					 0,
					 tasks_DFullLDLt[j].size(),
					 nops);
	queue_tmp0.push_back(tmp);
      } // loop : d
      queue_tmp0.sort(C_task_seq_complexity_greater);
      {
	int itmp = 0;
	for (list<C_task_seq*>::const_iterator it = queue_tmp0.begin();
	     it != queue_tmp0.end(); ++it, itmp++) {
	  queue_lists[itmp % num_threads].push_back(*it);
	} // loop : it
      }
      // begin : scope for task_name;
      {
	string task_name = "t " + to_string(level);
	for (int p = 0; p < num_threads; p++) {
#if 0
	  cerr << __FILE__ << " " << __LINE__ << " " << endl;
	  for (list<C_task_seq*>::const_iterator it = queue_lists[p].begin();
	       it != queue_lists[p].end(); ++it) {
	    cerr << (*it)->task_id
		 << " [ " << (*it)->end << " " << (*it)->begin << " ] ";
	  }
	  cerr << endl;	  
#endif
	  copytask_list2seq(_queue_static[p],
			    queue_lists[p],
			    queue_null_,
    //			    queue_dummy_,
			    task_name, 
			    C_DFULL,
			    p,
			    TASK_SINGLE,
			    1,
			    level,
			    (level_last - level) * 6);
	}
	// erase temporary C_task_seq whose elements are copied to 
	// queue_static[]
	for (list<C_task_seq*>::iterator it = queue_tmp0.begin();
	     it != queue_tmp0.end(); ++it) {
	  delete (*it);
	  (*it) = NULL;
	} 
      } // end : scope for task_name
      if ((level > 0)) { // && (level < (level_last - 1))) {
	for (int p = 0; p < num_threads; p++) {
	  for(list<C_task_seq *>::const_iterator it = queue_dynamic[level].begin();
	      it != queue_dynamic[level].end(); ++it) {
	    _queue_static[p].push_back(*it);
	  }
	} // loop : p
      }
      // dependency
    } // if (num_leaves < num_threads)
  }   // loop : level
// begin : counting task size
  {
    for (int p = 0; p < num_threads; p++) {
#ifdef DEBUG_PREPARE_THREAD
      cout << "thread = " << p << " @ " << _queue_static[p].size() << " / ";
#endif
      for (list<C_task_seq*>::const_iterator it = _queue_static[p].begin();
	   it != _queue_static[p].end(); ++it) {
	if ((*it)->parallel_single == TASK_SINGLE) { 
#ifdef DEBUG_PREPARE_THREAD
	  cout << (*it)->task_name << " [ "
	       << (*it)->begin << " : " << (*it)->end << " ] ; "
	       << (*it)->ops_complexity << " / ";
#endif
	  nops_queue[p][(*it)->level] += (*it)->ops_complexity;
	}
	if ((*it)->parallel_single == TASK_PARALLEL) {
#ifdef DEBUG_PREPARE_THREAD
	  cout << (*it)->task_name << " [[ "
	       << (*it)->begin << " : " << (*it)->end << " ]] ; "
	       << (*it)->ops_complexity << " / ";
#endif
	  if (((*it)->num_threads > 0) && ((*it)->task_id == C_DFULL)) {
	    //	    cout << "mutex_id = " << (*it)->mutex_id;
	    nops_queue[p][(*it)->level] += 
	      (*it)->ops_complexity / (long)(*it)->num_threads;
	  }
	  if ((*it)->level == 1) {
	    nops_queue[num_threads][(*it)->level] += (*it)->ops_complexity;
	  }
	}
      }

#ifdef DEBUG_PREPARE_THREAD
      cout << endl;
#endif
    }
#ifdef DEBUG_PREPARE_THREAD
    cout << "shared / ";
#endif
    for (int level = level_last; level >=0; level--) {
#ifdef DEBUG_PREPARE_THREAD
      cout << "level = " << level << endl;
#endif
      for (list<C_task_seq*>::const_iterator it = queue_dynamic[level].begin();
	   it != queue_dynamic[level].end(); ++it) {
#ifdef DEBUG_PREPARE_THREAD
	cout << (*it)->task_name << " [ "
	     << (*it)->begin << " : " << (*it)->end << " ] ; "
	     << (*it)->ops_complexity << " / ";
#endif
	nops_queue[num_threads][(*it)->level] += (*it)->ops_complexity;
      }
#ifdef DEBUG_PREPARE_THREAD
      cout << endl;
#endif
    }
#ifdef DEBUG_PREPARE_THREAD
    for (int level = (level_last + 1); level >= 0; level--) {
      cout << "level = " << level << " : ";
      for (int j = 0; j <= num_threads; j++) { 
	cout << nops_queue[j][level] << " ";
      }
      cout << endl;
    }
#endif
  }
// end : counting task size
// copy tasks in reverse order
  {
    int itmp = 0;  
#if 0
    itmp += queue_dynamic0->size();          // C_SPARSELOCALSCHUR
#endif
    for (int i = 0; i <= level_last; i++) {
      itmp += queue_dynamic[i].size();
      itmp += queue_dynamic1[i].size();
      itmp += queue_dynamic2[i].size();     // C_DEALLOCATE
      //      itmp += queue_dynamic3[i].size();     // C_DEALLOCATE
    }
    _queue_dynamic->reserve(itmp);
  }
#if 0
  for (list<C_task_seq*>::const_iterator it = queue_dynamic0->begin();
       it != queue_dynamic0->end(); ++it) {
      _queue_dynamic->push_back(*it);
  }
#endif
  for (int i = 0; i <= level_last; i++) {
    const int level = level_last - i;
    for (list<C_task_seq*>::const_iterator it = queue_dynamic[level].begin();
	 it != queue_dynamic[level].end(); ++it) {
      _queue_dynamic->push_back(*it);
    }
    for (list<C_task_seq*>::const_iterator it = queue_dynamic1[level].begin();
	 it != queue_dynamic1[level].end(); ++it) {
      _queue_dynamic->push_back(*it);
    }
    for (list<C_task_seq*>::const_iterator it = queue_dynamic2[level].begin();
	 it != queue_dynamic2[level].end(); ++it) {
      _queue_dynamic->push_back(*it);
    }
#if 0
    for (list<C_task_seq*>::const_iterator it = queue_dynamic3[level].begin();
	 it != queue_dynamic3[level].end(); ++it) {
      _queue_dynamic->push_back(*it);
    }
#endif
  }
#if 0
 for (int p = 0; p < num_threads; p++) {
   C_task_seq * queue_tmp = _queue_static[p].back();
   cerr << p << " : " << queue_tmp->task_id << " : " 
	<< queue_tmp->task_name << endl;
 }
#endif
// begin : scope for fout
 if (0) {
   char fname[256];
   int pid = get_process_id();
   sprintf(fname, "tasks-created.%d.data", pid);
   FILE *fp;
   if ((fp = fopen(fname, "a")) != NULL) {
     write_dependency(fp);
   }
   else {
    fprintf(stderr,
	    "%s %d : fail to open %s\n",
	    __FILE__, __LINE__, fname);
    exit(-1);
   }
 }
  queue_dynamic0->clear();
  for (int i = 0; i < (level_last + 1); i++) {
    queue_dynamic[i].clear();
    queue_dynamic1[i].clear();
    queue_dynamic2[i].clear();
    //    queue_dynamic3[i].clear();
  }
  delete [] queue_dynamic;
  delete [] queue_dynamic1;
  delete [] queue_dynamic2;
  //  delete [] queue_dynamic3;
  delete queue_dynamic0;
  for (int i = 0; i < num_threads; i++) {
    queue_lists[i].clear();
  }
  delete [] queue_lists;
}

void QueueRuntime::write_dependency(FILE *fp)
{
  for (int p = 0; p < _num_threads; p++) {
    fprintf(fp, "*** thread = %d @ %d ***\n",
	    p, (int)_queue_static[p].size());
    for (list<C_task_seq*>::const_iterator it = _queue_static[p].begin();
	 it != _queue_static[p].end(); ++it) {
      fprintf(fp, "** %d C_task_seq = %d %s %s **\n",
	      p, (*it)->task_id, (*it)->task_name,
	      (((*it)->parallel_single == TASK_SINGLE) ? " single " : " shared "));
      for (int j = (*it)->begin; j < (*it)->end; j++) {
	fprintf(fp, "%s : ", (*(*it)->queue)[j]->task_name);
	if ((*(*it)->queue)[j]->atomic_size > 1) {
	  fprintf(fp, "atomic %d / %d ",
		  (*(*it)->queue)[j]->atomic_id,
		  (*(*it)->queue)[j]->atomic_size);
	}
//	list<C_task *>& parents = *(*(*it)->queue)[j]->parents_work;
	list<C_task *>& parents = *(*(*it)->queue)[j]->parents;
	fprintf(fp, "parents = %d / ", (int)parents.size());
	for (list<C_task *>::const_iterator jt = parents.begin();
	     jt != parents.end(); ++jt) {
	  fprintf(fp, "[%p] %s / ", (*jt), (*jt)->task_name);
	}
	fprintf(fp, "\n");
      }
    } // loop : it
  }   // loop : p
  fprintf(fp, "*** dynamic : %d ***\n", (int)_queue_dynamic->size());
  for (vector<C_task_seq*>::const_iterator it = _queue_dynamic->begin();
       it != _queue_dynamic->end(); ++it) {
    fprintf(fp, "** C_task_seq = %d : %s [ %d %d] %d **\n",
	    (*it)->task_id,
	    (*it)->task_name,
	    (*it)->begin, (*it)->end,
	    (*it)->pos);
    for (int j = (*it)->begin; j < (*it)->end; j++) {
      fprintf(fp, "%s : ", (*(*it)->queue)[j]->task_name);
      if ((*(*it)->queue)[j]->atomic_size > 1) {
	fprintf(fp, " atomic %d / %d ",
		(*(*it)->queue)[j]->atomic_id,
		(*(*it)->queue)[j]->atomic_size);
      }
//    list<C_task *>& parents = *(*(*it)->queue)[j]->parents_work;
      list<C_task *>& parents = *(*(*it)->queue)[j]->parents;
      
      fprintf(fp, "parents = %d / ", (int)parents.size());
      for (list<C_task *>::const_iterator jt = parents.begin();
	   jt != parents.end(); ++jt) {
	fprintf(fp, "[%p] %s / ", (*jt), (*jt)->task_name);
      }
      fprintf(fp, "\n");
    }
  } // loop : it
}  // end : scope for fout

void QueueRuntime::exec_symb_fact()
{
  const int num_threads = _num_threads;
  void* results;
  pthread_attr_t th_attr;
  pthread_t *threads;

  clock_t t0_cpu, t1_cpu;
  elapsed_t t0_elapsed, t1_elapsed;
  //  struct timespec ts0, ts1;
  int ierr;
  if (_verbose) {
    fprintf(_fp, "symbolic factorization of sparse matrices with %d threads\n",
	    num_threads);
  }
  threads = new pthread_t[num_threads]; 

  ierr = pthread_mutex_init(&_mutex_root, NULL);
  if (_verbose && (ierr != 0)) {
    fprintf(_fp, "%s %d : pthread_mutex_init(&_mutex_root, NULL) : %d\n",
	    __FILE__, __LINE__, ierr);
  }
#ifdef DEBUG_EXEC_THREAD
  ierr = pthread_mutex_init(&_mutex_debug, NULL);
  if (_verbose && (ierr != 0)) {
    fprintf(_fp, "%s %d : pthread_mutex_init(&_mutex_debug, NULL) %d\n", 
	    __FILE__, __LINE__, ierr);
  }
#endif
  pthread_attr_init(&th_attr);
  pthread_attr_setdetachstate(&th_attr, PTHREAD_CREATE_JOINABLE);       

#ifdef DEBUG_EXEC_THREAD_FILE
  {
    int pid = get_process_id();
    char fname[256];
    sprintf(fname, "task-s.%d.data", pid);
    _fout.open(fname);
  }
#endif

  t0_cpu = clock();
  get_realtime(&t0_elapsed);

  THREAD_QUEUE_EXEC **params = new THREAD_QUEUE_EXEC*[num_threads];
  for (int p = 0; p < num_threads; p++) {
    params[p] = new THREAD_QUEUE_EXEC(p, num_threads, this);
    int pid = pthread_create(&threads[p], &th_attr, 
			     &thread_queue_symb_factorize_,
			     (void *)params[p]);
    if (pid != 0) {
      if (_verbose) {
	fprintf(_fp, "bad thread creation ? : %d\n", pid);
      }
      exit(0);
    }
  }
  pthread_attr_destroy(&th_attr);
  for (int p = 0; p < num_threads; p++) {
    int pid = pthread_join(threads[p], &results);
    if (pid != 0) {
      if (_verbose) {
	fprintf(_fp, "bad thread creation ? : %d\n", pid);
      }
      exit(0);
    }
    delete params[p];
  }
  delete [] params;

  t1_cpu = clock();
  get_realtime(&t1_elapsed);
  if (_verbose) {
    fprintf(_fp, 
	    "execution of symb queue : cpu time = %.4e elapsed time = %.4e\n", 
	    (double)(t1_cpu - t0_cpu) / (double)CLOCKS_PER_SEC,
	    convert_time(t1_elapsed, t0_elapsed));
  }
#ifdef DEBUG_EXEC_THREAD_FILE
  _fout.close();
#endif

#ifdef DEBUG_THREAD_TIME
  {
    int pid = get_process_id();
    char filename[256];		  
    FILE *fp;
    sprintf(filename, "threadtime-symb.%d.data", pid);
    if ((fp = fopen(filename, "w")) == NULL) {
      fprintf(stderr, "%s %d : fail to open %s\n",
	      __FILE__, __LINE__, filename);
      exit(-1);
    }
    else {
      for (int j = _queue_symb->begin; j < _queue_symb->end; j++) {
	fprintf(fp, "%s\t%d\t%d\t%d\t%d\t%d\t%ld\n",
		(*_queue_symb->queue)[j]->task_name,
		(*_queue_symb->queue)[j]->thread_id,
		convert_sec((*_queue_symb->queue)[j]->t0),
		convert_microsec((*_queue_symb->queue)[j]->t0),
		convert_sec((*_queue_symb->queue)[j]->t1),
		convert_microsec((*_queue_symb->queue)[j]->t1),
		*((*_queue_symb->queue)[j]->ops_complexity));
      }
      fclose(fp);
    } // if fopen()
  }
#endif

  pthread_mutex_destroy(&_mutex_root);
#ifdef DEBUG_EXEC_THREAD
  pthread_mutex_destroy(&_mutex_debug);
#endif

  delete [] threads;
}

void QueueRuntime::exec_num_fact_debug()
{
  if (_verbose) {
    fprintf(_fp, "numerical factorization with single threads :: DEBUG\n");
  }
  int *permute_block = new int[SIZE_B1];
  list<C_task_seq*>::const_iterator *its = new list<C_task_seq*>::const_iterator[_num_threads];
  if (_verbose) {
    fprintf(_fp, "*** thread = 0 @ %d ***\n", (int)_queue_static[0].size());
  }
  for (int p = 0; p < _num_threads; p++) {
    its[p] = _queue_static[p].begin();
  }
  for (list<C_task_seq*>::const_iterator it = _queue_static[0].begin();
       it != _queue_static[0].end(); ++it) {
    if (_verbose) {
      fprintf(_fp, "** C_task_seq = %d : %s, %s **\n",
	      (*it)->task_id, (*it)->task_name , 
	      (((*it)->parallel_single == TASK_SINGLE) ? " single " : " shared "));
    } 
    for (int j = (*it)->begin; j < (*it)->end; j++) {
      if (_verbose) {
	fprintf(_fp, "%s : ", (*(*it)->queue)[j]->task_name);
	if ((*(*it)->queue)[j]->atomic_size > 1) {
	  fprintf(_fp, "atomic %d / %d ",
		  (*(*it)->queue)[j]->atomic_id,
		  (*(*it)->queue)[j]->atomic_size);
	}
	else {
	  fprintf(_fp, "\n");
	}
      }
      execute_task_debug(*it, j, permute_block, 0);
      if ((*(*it)->queue)[j]->quit_queue) {
	return;
      }
    }
    for (int p = 1 ; p < _num_threads; p++) {
      fprintf(_fp, "** C_task_seq = %d : %s %s **\n",
	      (*its[p])->task_id, (*its[p])->task_name,
	      (((*its[p])->parallel_single == TASK_SINGLE) ? " single " : " shared "));
      if ((*its[p])->parallel_single == TASK_SINGLE) {
	for (int j = (*its[p])->begin; j < (*its[p])->end; j++) {
	  if (_verbose) {
	    fprintf(_fp, "%s : ", (*(*its[p])->queue)[j]->task_name);
	    if ((*(*its[p])->queue)[j]->atomic_size > 1) {
	      fprintf(_fp, " atomic %d / %d\n",
		      (*(*its[p])->queue)[j]->atomic_id,
		      (*(*its[p])->queue)[j]->atomic_size);
	    }
	    else {
	      fprintf(_fp, "\n");
	    }
	  }
	  execute_task_debug(*its[p], j, permute_block, 0);
	  if ((*(*its[p])->queue)[j]->quit_queue) {
	    return;
	  }
	}
      }
      ++its[p];
    } // loop : p > 0
  } // loop : it
}

void QueueRuntime::execute_task_debug(C_task_seq *seq, int pos, 
				      int *permute_block, 
				      int pid)
{
  C_task *task = (*seq->queue)[pos];
#ifdef DEBUG_CHECKPARENTS_DONE
 // debugging
  int waiting;
  waiting = check_parents_done(task);
  if (waiting > 0) {
    {
      cerr << pid << " parents of task " << task->task_name 
	   << " not finished : ";
      for(list<C_task *>::const_iterator nt = task->parents_work->begin();
	  nt != task->parents_work->end(); ++nt) {
	cerr << (*nt)->task_name << " ";
      }
      cerr << endl;
    }
  }
  //
#endif
#if 0
  switch (task->task_id) {
  case C_DFULL_SYM_GAUSS:
    ((C_dfull_gauss_arg*)task->func_arg)->permute_block = 
      permute_block;
    break;
  }
#endif
  // accuessing to unsigned char does not need mutex
  task->status = TASK_WORKING;
  //
#ifdef DEBUG_THREAD_TIME
  get_realtime(&(task->t0));
  //  clock_gettime(CLOCK_MONOTONIC, &(task->t0));
#endif
  // debugging : 13 :Apr.2012 Atsushi
  task->func(task->func_arg);
  if (task->task_id == C_DFULL_SYM_GAUSS) {
    //    task->quit_queue = ((C_dfull_gauss_arg*)task->func_arg)->quit;
    for (int i = (pos + 1); i <= (pos + task->to_next_task); i++) {
      (*seq->queue)[i]->status = TASK_DONE;
      (*seq->queue)[i]->quit_queue = true;
    }
  }
  // accuessing to unsigned char does not neet mutex
  task->status = TASK_DONE;
  //
}

void QueueRuntime::exec_num_fact(const int called)
{
  //   unsigned int ui;
  const int num_threads = _num_threads;
  void* results;
  pthread_attr_t th_attr;
  pthread_t *threads;

  clock_t t0_cpu, t1_cpu;
  elapsed_t t0_elapsed, t1_elapsed;
  //  struct timespec ts0, ts1;
#ifdef DEBUG_PRINT_TASK
  _fps = new FILE*[num_threads];

  char fname[256];
  for (int p = 0; p < num_threads; p++) {
    sprintf(fname, "log.%04d.data", p);
    _fps[p] = fopen(fname, "w");
  }
#endif
  t0_cpu = clock(); 
  get_realtime(&t0_elapsed);
  //  clock_gettime(CLOCK_REALTIME, &ts0);
  if (_verbose) {
    fprintf(_fp, "numerical factorization with %d threads\n", num_threads);
  }
  _queue_dynamic_pos_start = (int)((double)(_queue_dynamic->size())
				   * RATIO_QUEUE_GREEDY);
  _queue_dynamic_pos = _queue_dynamic_pos_start;
  _queue_dynamic_notcopied = num_threads;

  // reset status for numeric queue execution
  _waiting_root = 0;
  _phase_dynamic = 0;

  for (int i = 0; i < DIST_TASK_CRITICAL; i++) {
    _zone_entered[i] = 0;
    _zone_finished[i] = 0;
  }

  for (int i = 0; i < DIST_TASK_CRITICAL; i++) {
    _zone_static_assigned[i] = 0;
  }

  for (int p = 0; p < num_threads; p++) {
    for (int i = 0; i < DIST_TASK_CRITICAL; i++) {
      _begins[p][i] = 0;
      _ends[p][i] = 0;
    }
    for (int q = 0; q < DIST_TASK_CRITICAL; q++) {
       for (int i = 0; i < DIST_TASK_CRITICAL; i++) {
	_group_entered[p][q][i] = 0;
	_group_finished[p][q][i] = 0;
	_group_task_ends[p][q][i] = (-1);
	_group_static_assigned[p][q][i] = (-1);
      }
    } // loop : q
  } // loop : p

  pthread_mutex_init(&_mutex_file, NULL);
  pthread_mutex_init(&_mutex_dependency, NULL);

  pthread_attr_init(&th_attr);
  pthread_attr_setdetachstate(&th_attr, PTHREAD_CREATE_JOINABLE);       

  for (int p = 0; p < num_threads; p++) {
    pthread_mutex_init(&_mutex_group[p], NULL);
  }
  threads = new pthread_t[num_threads]; 

  pthread_mutex_init(&_mutex_root, NULL);
  pthread_mutex_init(&_mutex_debug, NULL);
  pthread_cond_init(&_cond_root, NULL);

#ifdef DEBUG_EXEC_THREAD_FILE
  {
    int pid = get_process_id();
    char fname[256];
    sprintf(fname, "task-n.%d.data", pid);
    cerr << "== task log == " << fname << endl;
    _fout.open(fname);
  }
#endif
#if 1
 {
  // copy dependency data : parents -> parents_work 
   for (int p = 0; p < _num_threads; p++) {
     for (list<C_task_seq *>::const_iterator it = _queue_static[p].begin();
	  it != _queue_static[p].end(); ++it) {
       vector<C_task *> &queue = *(*it)->queue;
       for (int j = (*it)->begin; j < (*it)->end; j++) {
	 list<C_task *>& parents_work = *(queue[j]->parents_work);
#if 0
	 cerr << queue[j]->task_name 
	      << " parents_work.size() = " << parents_work.size() << endl;
#endif
	 if (parents_work.size() > 0) { // to avoid double free
	   parents_work.clear();  
	 }
       }  // loop : j
     }    // loop : it
   }
   for (vector<C_task_seq *>::const_iterator it = _queue_dynamic->begin();
	it != _queue_dynamic->end(); ++it) {
     vector<C_task *> &queue = *((*it)->queue);
     for (int j = (*it)->begin; j < (*it)->end; j++) {
       list<C_task *>& parents_work = *(queue[j]->parents_work);
#if 0
       cerr << queue[j]->task_name 
	    << " parents_work.size() = " << parents_work.size() << endl;
#endif
       if (parents_work.size() > 0) {  // to avoid double free
	 parents_work.clear(); 
       }
     }  // loop : j
   }    // loop : it


   for (int p = 0; p < _num_threads; p++) {
     for (list<C_task_seq *>::const_iterator it = _queue_static[p].begin();
	  it != _queue_static[p].end(); ++it) {
       vector<C_task *> &queue = *(*it)->queue;
       for (int j = (*it)->begin; j < (*it)->end; j++) {
	 queue[j]->quit_queue = false;
	 list<C_task *>& parents = *(queue[j]->parents);
	 list<C_task *>& parents_work = *(queue[j]->parents_work);
	 if (parents_work.size() == 0) {
#ifdef SX_ACE 
	   for (list<C_task *>::const_iterator lt = parents.begin();
		lt != parents.end(); ++lt) {
	     parents_work.push_back(*lt);
	   }
#else // SX-ACE C++ rev 110 with C++98/03 does not understand back_insertera
	   std::copy(parents.begin(), parents.end(), 
		     back_inserter(parents_work));
#endif
	 }
#if 0   // verify task name with printing single/parallel task sequence
	 else {
	   cerr << queue[j]->task_name << " " 
		<< (*it)->task_name << " "
		<< ((*it)->parallel_single == TASK_SINGLE ? "serl" : "para")
	      << " parents_work.size() = " << parents_work.size() << endl;
	 }
#endif
	 queue[j]->status = TASK_WAITING;  // reset status
       }  // loop : j
       (*it)->pos = (*it)->begin;
     }    // loop : it
   }

   for (vector<C_task_seq *>::const_iterator it = _queue_dynamic->begin();
	it != _queue_dynamic->end(); ++it) {
     vector<C_task *> &queue = *((*it)->queue);
     for (int j = (*it)->begin; j < (*it)->end; j++) {
       queue[j]->quit_queue = false;
       list<C_task *>& parents = *(queue[j]->parents);
       list<C_task *>& parents_work = *(queue[j]->parents_work);
       if (parents_work.size() == 0) {
#ifdef SX_ACE
	   for (list<C_task *>::const_iterator lt = parents.begin();
		lt != parents.end(); ++lt) {
	     parents_work.push_back(*lt);
	   }
#else // SX-ACE C++ rev 110 with C++98/03 does not understand back_inserter
	 std::copy(parents.begin(), parents.end(), back_inserter(parents_work));
#endif     
       }
#if 0    // verify task name with printing single/parallel task sequence
       else {
	 cerr << queue[j]->task_name << " " 
	      << (*it)->task_name << " "
	      << ((*it)->parallel_single == TASK_SINGLE ? "serl" : "para")
	      << " parents_work.size() = " << parents_work.size() << endl;
       }
#endif
       queue[j]->status = TASK_WAITING;  // reset status
     }  // loop : j
     (*it)->pos = (*it)->begin;
   }    // loop : it
 }
 #endif
// begin : scope for fout
  if (0) {
    char fname[256];
    int pid = get_process_id();
    sprintf(fname, "tasks-before.%d.%d.data", pid, called);
    FILE *fp;
    if ((fp = fopen(fname, "a")) != NULL) {
      write_dependency(fp);
      fclose(fp);
    }
    else {
      fprintf(stderr,
	      "%s %d : fail to open %s\n",
	      __FILE__, __LINE__, fname);
      exit(-1);
    }
  }  // end : scope for fout

  t0_cpu = clock();
  get_realtime(&t0_elapsed);

  THREAD_QUEUE_EXEC **params = new THREAD_QUEUE_EXEC*[num_threads];
  for (int p = 0; p < num_threads; p++) {
    params[p] = new THREAD_QUEUE_EXEC(p, num_threads, this);
    int pid = pthread_create(&threads[p], &th_attr, 
			     &thread_queue_num_factorize_,
			     (void *)params[p]);
    if (pid != 0) {
      if (_verbose) {
	fprintf(_fp, "bad thread creation ? : %d\n", pid);
      }
      exit(0);
    }
  }
  pthread_attr_destroy(&th_attr);
  for (int p = 0; p < num_threads; p++) {
    int pid = pthread_join(threads[p], &results);
    if (pid != 0) {
      if (_verbose) {
	fprintf(_fp, "bad thread creation ? : %d\n", pid);
      }
      exit(0);
    }
    delete params[p];
  }
  delete [] params;

  t1_cpu = clock();
  get_realtime(&t1_elapsed);
  if (_verbose) {
    fprintf(_fp,
	    "execution of numr queue : cpu time = %.4e elapsed time = %.4e\n", 
	    (double)(t1_cpu - t0_cpu) / (double)CLOCKS_PER_SEC,
	    convert_time(t1_elapsed, t0_elapsed));
  }
#ifdef DEBUG_EXEC_THREAD_FILE
  _fout.close();
#endif

  // begin : scope for fout
  if (0) {
    char fname[256];
    int ppid = get_process_id();
    sprintf(fname, "tasks-copied.%d.data", ppid);
    FILE *fp;
    if ((fp = fopen(fname, "a")) != NULL) {
      write_dependency(fp);
      fclose(fp);
    }
    else {
      fprintf(stderr,
	      "%s %d : fail to open %s\n",
	      __FILE__, __LINE__, fname);
      exit(-1);
    }
  }

#ifdef DEBUG_THREAD_TIME
  {
    int pid = get_process_id();
    char filename[256];		  
    FILE *fp;

    for (int p = 0; p < num_threads; p++) {
      for (list<C_task_seq*>::const_iterator it = _queue_static[p].begin();
	   it != _queue_static[p].end(); ++it) {
	(*it)->status = TASK_WAITING;
      }
    }
    for (vector<C_task_seq*>::const_iterator it = _queue_dynamic->begin();
	 it != _queue_dynamic->end(); ++it) {
      (*it)->status = TASK_WAITING;
    }
    sprintf(filename, "threadtime-num.%d.data", pid);
    if ((fp = fopen(filename, "w")) == NULL) {
      fprintf(stderr, "%s %d : fail to open %s\n",
	      __FILE__, __LINE__, filename);
      exit(-1);
    }
    else {
      fprintf(fp, "**** numeric factorization start ****\n");
      for (int p = 0; p < num_threads; p++) {
	fprintf(fp, "*** thread = %d @ %d ** ", p,
		(int)_queue_static[p].size());
	for (list<C_task_seq*>::const_iterator it = _queue_static[p].begin();
	     it != _queue_static[p].end(); ++it) {
	  fprintf(fp, "** %d  C_task_seq = %d : %s %s **\n",
		  p, (*it)->task_id,
		  ((*it)->parallel_single == TASK_SINGLE ? "distributed " : "shared "),
		(*it)->task_name);
	  if ((*it)->status == TASK_DONE) {
	    continue;
	  }
	  for (int j = (*it)->begin; j < (*it)->end; j++) {
	    
	    if ((*(*it)->queue)[j]->task_id == C_SPARSESCHUR) {
	      elapsed_t t0, t1;
	      t0 = (*(*it)->queue)[j]->t0;
	      t1 = ((C_SparseNumFact_arg<double, double> *)(*(*it)->queue)[j]->func_arg)->tt[0];
	      string task_name = "o0 : ";
	      fprintf(fp, "%s\t%d\t%d\t%d\t%d\t%d\t%ld\n",
		      task_name.c_str(),
		      //			(*(*it)->queue)[j]->task_name,
		      (*(*it)->queue)[j]->thread_id,
		      convert_sec(t0),
		      convert_microsec(t0),
		      convert_sec(t1),
		      convert_microsec(t1), 0L);
	      for (int i = 0; i < 4; i++) {
		string task_name = "o" + to_string(i + 1) + " : ";

		t0 = ((C_SparseNumFact_arg<double, double> *)(*(*it)->queue)[j]->func_arg)->tt[i];
		t1 = ((C_SparseNumFact_arg<double, double> *)(*(*it)->queue)[j]->func_arg)->tt[i + 1];		      
		fprintf(fp, "%s\t%d\t%d\t%d\t%d\t%d\t%ld\n",
			task_name.c_str(),
			//			(*(*it)->queue)[j]->task_name,
			(*(*it)->queue)[j]->thread_id,
			convert_sec(t0),
			convert_microsec(t0),
			convert_sec(t1),
			convert_microsec(t1), 0L);
	      }
	    }
	    else {
	      fprintf(fp, "%s\t%d\t%d\t%d\t%d\t%d\t%ld\n",
		      (*(*it)->queue)[j]->task_name,
		      (*(*it)->queue)[j]->thread_id,
		      convert_sec((*(*it)->queue)[j]->t0),
		      convert_microsec((*(*it)->queue)[j]->t0),
		      convert_sec((*(*it)->queue)[j]->t1),
		      convert_microsec((*(*it)->queue)[j]->t1),
		      *((*(*it)->queue)[j]->ops_complexity));
	    }
	    //*((C_SparseNumFact_arg<double> *)(*(*it)->queue)[j]->func_arg)->nops :
	  }
	  (*it)->status = TASK_DONE;
	}
      }
      fprintf(fp, "*** dynamic : %d ***", (int)_queue_dynamic->size());
      for (vector<C_task_seq*>::const_iterator it = _queue_dynamic->begin();
	   it != _queue_dynamic->end(); ++it) {
	fprintf(fp, "** C_task_seq = %s : [%d, %d] %d **\n",
	      (*it)->task_name,
		(*it)->begin, (*it)->end,(*it)->pos);
	if ((*it)->status == TASK_DONE) {
	  continue;
	}
	for (int j = (*it)->begin; j < (*it)->end; j++) {
	  fprintf(fp, "%s\t%d\t%d\t%d\t%d\t%d\t%ld\n",
		  (*(*it)->queue)[j]->task_name,
		  (*(*it)->queue)[j]->thread_id,
		  convert_sec((*(*it)->queue)[j]->t0),
		  convert_microsec((*(*it)->queue)[j]->t0),
		  convert_sec((*(*it)->queue)[j]->t1),
		  convert_microsec((*(*it)->queue)[j]->t1),
		  *((*(*it)->queue)[j]->ops_complexity));
	}
      (*it)->status = TASK_DONE;
      }
      fprintf(fp, "**** numeric factorization end ****\n");
      fclose(fp);
    }
  } // if fopen()
#endif // DEBUG_THREAD_TIME
#ifdef DEBUG_PRINT_TASK
  for (int p = 0; p < num_threads; p++) {
    fclose(_fps[p]);
  }
  delete [] _fps;
#endif
  delete [] threads;
}

void *thread_queue_symb_factorize_(void *arg)
{
  THREAD_QUEUE_EXEC *params = (THREAD_QUEUE_EXEC *)arg;
  params->dissectionRuntime->thread_queue_symb_factorize(params->id, 
							 params->num_threads);
  pthread_exit(arg);

  return (void *)NULL;
}

void QueueRuntime::thread_queue_symb_factorize(const int pid,
					       const int num_threads)
{
  //  const int pid = params->id;
  //  const int num_threads = params->num_threads;
  // greedy -- better to be in seperated function
  int pos, end;
  C_task_seq* it = _queue_symb;

  while(1) {
    pthread_mutex_lock(&_mutex_root);
    {
      pos = it->pos;
      end = it->end;
      if (pos < end) {
	it->pos++;
      }
    }
    pthread_mutex_unlock(&_mutex_root);
    if (pos >= end) {
      break;
    }
    C_task *task = (*it->queue)[pos];
    task->status = TASK_WORKING;
#ifdef DEBUG_THREAD_TIME
    get_realtime(&(task->t0));
    //    clock_gettime(CLOCK_MONOTONIC, &(task->t0));
#endif
    task->func(task->func_arg);
#ifdef DEBUG_THREAD_TIME
    get_realtime(&(task->t1));
    //    clock_gettime(CLOCK_MONOTONIC, &(task->t1));
    task->thread_id = pid;
#endif
  // accuessing to unsigned char does not neet mutex
    task->status = TASK_DONE;

#ifdef DEBUG_EXEC_THREAD
    pthread_mutex_lock(&_mutex_debug);
    {
      cerr << "( " << pid << " " << pos << " ) " ;
    }
    pthread_mutex_unlock(&_mutex_debug);
#endif
  } // while (1)
#ifdef DEBUG_EXEC_THREAD
  pthread_mutex_lock(&_mutex_debug);
  {
    cerr << pid << " " << it->task_name 
	 << " @ " << it->task_id << " greedy end." << endl;
  }
  pthread_mutex_unlock(&_mutex_debug);
#endif     
}

void *thread_queue_num_factorize_(void *arg)
{
  THREAD_QUEUE_EXEC *params = (THREAD_QUEUE_EXEC *)arg;
  params->dissectionRuntime->thread_queue_num_factorize(params->id, 
							params->num_threads);
  pthread_exit(arg);

  return (void *)NULL;
}

void QueueRuntime::thread_queue_num_factorize(const int pid, 
					      const int num_threads)
{
  int *permute_block = new int[SIZE_B1];

  list<C_task_seq*>::const_iterator it = _queue_static[pid].begin();
  int zone, zone_idxn, zone_idxp;
  bool zone_flag;
  int zone_first_entered = 0;
  int zone_last_entered = 0;
  int zone_first_finished = 0;
  int zone_last_finished = 0;
  int cnt_cdfull = 0;
  // 
  zone = 0;
  while(it != _queue_static[pid].end()) {
#ifdef DEBUG_THREAD_LOOP
    pthread_mutex_lock(&_mutex_debug);
    {
      cerr << pid << " : " << (*it)->task_name 
	   << " _phase_dynamic = " << _phase_dynamic << endl;
    }
    pthread_mutex_unlock(&_mutex_debug);
#endif
    //
    zone_idxp = (zone + DIST_TASK_CRITICAL - 1) % DIST_TASK_CRITICAL;
    zone_idxn = zone % DIST_TASK_CRITICAL;
    zone_first_entered = 0;
    zone_last_entered = 0;
    pthread_mutex_lock(&_mutex_root);
    {
      if (_zone_entered[zone_idxn] == 0) {
	zone_first_entered = 1;
	_zone_static_assigned[zone_idxn] = 0;
      }
      _zone_entered[zone_idxn]++;
      zone_flag = (_zone_entered[zone_idxn] == num_threads);
      if (zone_flag) {
	// clear cyclic buffer
 	_zone_entered[zone_idxp] = 0;
	_zone_static_assigned[zone_idxp] = 0;  // only fol safety
	zone_last_entered = 1;
      }
    }
    pthread_mutex_unlock(&_mutex_root);
    zone++;
    //
    if ((*it)->parallel_single == TASK_PARALLEL) {
      switch((*it)->task_id) {
      case C_DFULL :
#ifdef DEBUG_EXEC_THREAD
	pthread_mutex_lock(&_mutex_debug);
	{
	  cerr << pid << " C_DFULL " 
	       << (*it)->task_name << " : " << (*it)->task_id 
	       << " @ " << zone << " : " << zone_flag << " : "
	       << " zone ( " << zone_first_entered << " : "
	       << zone_last_entered  << " )"
	       << " : " << (*it)->pos << " " << (*it)->end << endl;
	}
	pthread_mutex_unlock(&_mutex_debug);
#endif
	thread_queue_C_DFULL(pid, num_threads, *it, 
			     permute_block, 
			     cnt_cdfull,
			     zone_first_entered,
			     zone_idxn);
#ifdef DEBUG_THREAD_LOOP
	pthread_mutex_lock(&_mutex_debug);
	{
	  cerr << pid << " C_DFULL " 
	       << (*it)->task_name << " end." << endl;
	}
	pthread_mutex_unlock(&_mutex_debug);
#endif
	break;
      case C_SPARSELOCALSCHUR:
      case C_DTRSM:
      case C_DGEMM:
      case C_DSUB:
      case C_FILLMATRIX:
      case C_DEALLOCATE:
      case C_SPARSESOLVER: 
#ifdef DEBUG_EXEC_THREAD
	pthread_mutex_lock(&_mutex_debug);
	{
	  cerr << pid << " remaining diagonal " 
	       << (*it)->task_name << " : " << (*it)->task_id 
	       << " @ " << zone << " : " << zone_flag << " : "
	       << " zone ( " << zone_first_entered << " : "
	       << zone_last_entered  << " )"
	       << " : " << (*it)->pos << " " << (*it)->end << endl;
	}
	pthread_mutex_unlock(&_mutex_debug);
#endif
	thread_queue_parallel_dynamic(pid, num_threads, *it,
				      permute_block, 
				      zone_idxn);
#ifdef DEBUG_THREAD_LOOP
	pthread_mutex_lock(&_mutex_debug);
	{
	  cerr << pid << " remaining diagonal " 
	       << (*it)->task_name << " end." << endl;
	}
	pthread_mutex_unlock(&_mutex_debug);
#endif
	break;
      case C_SPARSELOCALSCHUR1:
      case C_DTRSM1:
      case C_DGEMM1:
      case C_DSUB1:
      case C_FILLMATRIX1:
      case C_DEALLOCATE1:
#ifdef DEBUG_EXEC_THREAD
	pthread_mutex_lock(&_mutex_debug);
	{
	  cerr << pid << " remaining off-diagonal " 
	       << (*it)->task_name << " : " << (*it)->task_id 
	       << " @ " << zone << " : " << zone_flag << " : "
	       << " zone ( " << zone_first_entered << " : "
	       << zone_last_entered  << " )"
	       << " : " << (*it)->pos << " " << (*it)->end << endl;
	}
	pthread_mutex_unlock(&_mutex_debug);
#endif
	thread_queue_parallel_static(pid, num_threads, *it, 
				     permute_block, 
				     zone_idxp,
				     zone_idxn,
				     zone_first_entered,
				     zone_last_entered);
#ifdef DEBUG_THREADS_LOOP
	pthread_mutex_lock(&_mutex_debug);
	{
	  cerr << pid << " remaining off-diagonal " 
	       << (*it)->task_name << " end." << endl;
	}
	pthread_mutex_unlock(&_mutex_debug);
#endif
	break;
      }
    }
    else {  //     if ((*it)->parallel_single == TASK_SINGLE)
#ifdef DEBUG_EXEC_THREAD
      pthread_mutex_lock(&_mutex_debug);
      {
	cerr << pid << " single " 
	     << (*it)->task_name << " : " << (*it)->task_id 
	     << " @ " << zone << " : " << zone_flag << " : "
	     << " zone ( " << zone_first_entered << " : "
	     << zone_last_entered  << " )"
	     << " : " << (*it)->pos << " " << (*it)->end << endl;
      }
      pthread_mutex_unlock(&_mutex_debug);
#endif
      thread_queue_single(pid, num_threads, *it,
			  permute_block, 
			  zone_idxn);
#ifdef DEBUG_THREAD_LOOP
	pthread_mutex_lock(&_mutex_debug);
	{
	  cerr << pid << " single " 
	       << (*it)->task_name << " end." << endl;
	}
	pthread_mutex_unlock(&_mutex_debug);
#endif
    } //  if ((*it)->parallel_single == TASK_SINGLE) 
    int clear_flag = 0;
    pthread_mutex_lock(&_mutex_root);
    {
      if (_zone_finished[zone_idxn] == 0) {
	zone_first_finished = 1;
      }
      _zone_finished[zone_idxn]++;
      if (_zone_finished[zone_idxn] == num_threads) {
	zone_last_finished = 1;
	// clear status of cyclic buffer
 	_zone_finished[zone_idxp] = 0;
	// the last task
	(*it)->status = TASK_DONE;
	if((*it)->task_id == C_DFULL) {
	  clear_flag = 1;
	}
#ifdef DEBUG_THREAD_LOOP
	pthread_mutex_lock(&_mutex_debug);
	{
	  cerr << pid << " " 
	       << (*it)->task_name << " end." << endl;
	}
	pthread_mutex_unlock(&_mutex_debug);
#endif
      }
    }
    pthread_mutex_unlock(&_mutex_root);
    // without mutex
    if((*it)->task_id == C_DFULL) {
#if 0
      pthread_mutex_lock(&_mutex_debug);
      {
	cerr << pid << " cnt_cdfull = " << cnt_cdfull << " : "
	     << (*it)->task_id << " : " << (*it)->task_name << endl;
      }
      pthread_mutex_unlock(&_mutex_debug);
#endif
      cnt_cdfull = (cnt_cdfull + 1) % DIST_TASK_CRITICAL;
      if (clear_flag) {
	const int cnt_cdfull_p = 
	  ((cnt_cdfull + DIST_TASK_CRITICAL - 1) % DIST_TASK_CRITICAL);
	for (int p = 0; p < num_threads; p++) {
	  for (int i = 0; i < DIST_TASK_CRITICAL; i++) {
	    _group_entered[p][cnt_cdfull_p][i] = 0;
	    _group_finished[p][cnt_cdfull_p][i] = 0;
	    _group_task_ends[p][cnt_cdfull_p][i] = (-1);
	  }
	}
      }
    } //     if((*it)->task_id == C_DFULL) {
    ++it;
  } // while (it != _queue_static[pid].end()) {
  delete [] permute_block;
}


void QueueRuntime::thread_queue_C_DFULL(const int pid,
					const int num_threads,
					C_task_seq *it,
					int *permute_block, 
					const int cnt_cdfull,
					const int zone_first_entered,
					const int zone_idxn)
{
  const int pid_g = it->mutex_id;
  int group = 0;
  int group_first_entered = 0, group_last_entered = 0;
  int group_first_finished = 0;
  int pos;
  int pid0;
  int *_grp_entered = _group_entered[pid_g][cnt_cdfull];
  int *_grp_finished = _group_finished[pid_g][cnt_cdfull];
  int *_grp_task_ends = _group_task_ends[pid_g][cnt_cdfull];
  int *_grp_static_assigned = _group_static_assigned[pid_g][cnt_cdfull];
  bool zone_flag;
  
  // get status other threads entered the same zone
  pthread_mutex_lock(&_mutex_root);
  {
    zone_flag = (_zone_entered[zone_idxn] == num_threads);
  }
  pthread_mutex_unlock(&_mutex_root);
  
#ifdef DEBUG_EXEC_THREAD
  pthread_mutex_lock(&_mutex_debug);
  {
    cerr << pid << " cnt_cdfull = " << cnt_cdfull << " parallel : c_dfull " 
	 << it->task_name << " pos = " << it->pos << endl;
  }
  pthread_mutex_unlock(&_mutex_debug);
#endif

  pthread_mutex_lock(&_mutex_group[pid_g]);
  {
    if (_grp_entered[0] == 0) {
      // reset flags to keep status of entered/finshed tasks
      for (int i = 1; i < DIST_TASK_CRITICAL; i++) {
	_grp_entered[i] = 0;
	_grp_finished[i] = 0;
      }
      for (int i = 0; i < it->num_threads; i++) {
	_group_nops[pid_g][i] = 0L;
      }
      group_first_entered = 1;
      _grp_static_assigned[0] = 0; // 12 Feb.2013 Atsusti : for safety?
#ifdef DEBUG_EXEC_THREAD
      pthread_mutex_lock(&_mutex_debug);	      
      {
	cerr << pid << " : " << it->pos << endl;
      }
      pthread_mutex_unlock(&_mutex_debug);
#endif
    }
    // pid0 is decided by arrived order
    pid0 = _grp_entered[0]++;
    if (_grp_entered[0] == it->num_threads) {
      group_last_entered = 1;
    }
  }
  pthread_mutex_unlock(&_mutex_group[pid_g]);


  //
  while (1) {    // global loop of task queue
    //
    const int gidxp = (group + DIST_TASK_CRITICAL - 1) % DIST_TASK_CRITICAL;
    const int gidxn = group % DIST_TASK_CRITICAL;
    int finish_group = 0;
    int skip_flag = 0;
    if (group > 0) {
      pthread_mutex_lock(&_mutex_group[pid_g]);
      {
	_begins_group[pid_g][pid0][gidxn] = (-1);
	_ends_group[pid_g][pid0][gidxn] = (-1);
	group_first_entered = 0;
	group_last_entered = 0;
	if (_grp_entered[gidxn] == 0) {
	  group_first_entered = 1;
	  // assuming other jobs are assiged by static
	  for (int i = 0; i < it->num_threads; i++) {
	    _group_nops[pid_g][i] = 0L;
	  }
#if 1
	  pos = it->pos; 
	  int jtmp = pos; 
	  for (int i = (*it->queue)[pos]->parallel_id;
	       i < (*it->queue)[pos]->parallel_max; i++) {
	    jtmp += (*it->queue)[jtmp]->atomic_size;
	  }
	  _grp_task_ends[gidxn] = jtmp;
	  _grp_task_ends[gidxp] = (-1);
	  _grp_static_assigned[gidxn] = 0;
#else
	  _group_task_ends[pid_g][cnt_cdfull][gidxn] = (-1);
#endif
	}
	_grp_entered[gidxn]++;
	if (_grp_entered[gidxn] == it->num_threads) {
	  // clear status of cyclic buffer
	  _grp_entered[gidxp] = 0;
	  group_last_entered = 1;
	}
	if ((_grp_finished[gidxn] > 0) && (group_first_finished == 0)) {
	  _grp_finished[gidxn]++;
	  if (_grp_finished[gidxn] == it->num_threads) {
	    // clear status of cyclic buffer
	    _grp_finished[gidxp] = 0;
	    _group_task_ends[pid_g][cnt_cdfull][gidxp] = (-1);
	  }
	  skip_flag = 1;
	} // if ((_grp_finished[gidxn] > 0) && (group_first_finished == 0))
      }
      pthread_mutex_unlock(&_mutex_group[pid_g]);
    }
    else {
      pthread_mutex_lock(&_mutex_group[pid_g]);
      {
	if ((_grp_finished[gidxn] > 0) && 
	    (group_first_finished == 0)) {
	  _grp_finished[gidxn]++;
	  if (_grp_finished[gidxn] == it->num_threads) {
	// clear status of cyclic buffer
	    _grp_finished[gidxp] = 0;
	    _group_task_ends[pid_g][cnt_cdfull][gidxp] = (-1);
	  }
	  skip_flag = 1;
	}
      }
      pthread_mutex_unlock(&_mutex_group[pid_g]);
    }
#ifdef DEBUG_EXEC_THREAD
    pthread_mutex_lock(&_mutex_debug);
    {
      cerr << pid << " pid_g = " << pid_g << " pid0 " << pid0
	   << " group " << group << " entered " 
	   << group_first_entered << " / " << group_last_entered 
	   << " (" << _grp_entered[gidxp] << "@" << gidxp 
	   << ") / (" << _grp_entered[gidxn] << "@" << gidxn 
	   << ") finished "
	   << " (" << _grp_finished[gidxp] << "@" << gidxp 
	   << ") / ( " << _grp_finished[gidxn] << "@" << gidxn 
	   << ") " << it->pos << " : " << it->end
	   << " skip_flag = " << skip_flag;
      if (skip_flag) {
	cerr << " task is already finished : my job is skipped";
      }
      cerr << " zone_flag = " << zone_flag
	   << " zone_idxn = " << zone_idxn
	   << endl;
    }
    pthread_mutex_unlock(&_mutex_debug);
#endif
    group++; 
    //
    if (skip_flag) {
      continue;
    }
#if 0  // this has no effect ? 20 Jun.2012 Atsushi
    if (group_first_entered) {
      pthread_mutex_lock(&_mutex_group[pid_g]);
      { 
	const int group_static_assigned = _grp_static_assigned[gidxn];
	if (group_static_assigned == 0) {
	  for (int p = 0; p < it->num_threads; p++) {
	    _begins_group[pid_g][p][gidxn] = (-1);
	    _ends_group[pid_g][p][gidxn] = (-1);
	  }
	}
      }
      pthread_mutex_unlock(&_mutex_group[pid_g]);
    }
#endif
#if 0
    if (group_first_entered) {
      pthread_mutex_lock(&_mutex_group[pid_g]);
      {
	pos = it->pos; 
	int jtmp = pos; 
	for (int i = (*it->queue)[pos]->parallel_id;
	     i < (*it->queue)[pos]->parallel_max; i++) {
	  jtmp += (*it->queue)[jtmp]->atomic_size;
	}
	const int end_queue = jtmp;
	_grp_task_ends[gidxn] = jtmp;
	_grp_task_ends[gidxp] = (-1);
	_grp_static_assigned[gidxn] = 0;
      }
      pthread_mutex_unlock(&_mutex_group[pid_g]);
    }
#endif
#if 1
    // wait until other threads enter the same zone
    while (!zone_flag) {
      // looking for dynamic queue
      int flag = 1;
      while ((flag == 1) && (!zone_flag)) {
	flag = execute_task_dynamic_buffer(permute_block, 
					   pid,
					   &_mutex_dependency);
	pthread_mutex_lock(&_mutex_root);
	{
	  zone_flag = (_zone_entered[zone_idxn] == num_threads);
	}
	pthread_mutex_unlock(&_mutex_root);
      }
      if (zone_flag) {
#ifdef DEBUG_EXEC_THREAD
	pthread_mutex_lock(&_mutex_debug);
	{
	  cerr << pid << " pid_g = " << pid_g << " pid0 " << pid0
	       << " zone_idxn = " << zone_idxn
	       << " zone_flag = " << zone_flag 
	       << endl;
	}
	pthread_mutex_unlock(&_mutex_debug);
#endif
	break;
      }
// get tasks with atomic_size from the queue : *it
      int itmp, jtmp;
      int exit_flag = 0;
      pthread_mutex_lock(&_mutex_group[pid_g]);
      {
	pos = it->pos;
	// do not increase queue more than its size
	if (pos == it->end) {
	  exit_flag = 1;
	}
	else {
	  // there might be segments with more than 1 automic size
	  // between parallel_id and parallel_max
	  itmp = pos + (*it->queue)[pos]->atomic_size;
	  jtmp = (itmp - 1 + ((*it->queue)[pos]->parallel_max -
			      (*it->queue)[pos]->parallel_id));
#ifdef DEBUG_EXEC_THREAD
	  pthread_mutex_lock(&_mutex_debug);
	  {
	    cerr << pid << " pos = " << pos << " +atomic_size = " << itmp
		 << " parallel_max " << jtmp << " : " << it->end 
		 << endl;
	  }
	  pthread_mutex_unlock(&_mutex_debug);
#endif
	  if (itmp <= jtmp) {
	    it->pos = itmp;
	  }
	}
      } 
      pthread_mutex_unlock(&_mutex_group[pid_g]);
      if (exit_flag) {
#ifdef DEBUG_EXEC_THREAD
	pthread_mutex_lock(&_mutex_debug);
	{
	  cerr << pid << " pod = " << pos << "it->end = " << it->end
	       << " : " << it->task_name << " already finished.." << endl;
	}
	pthread_mutex_unlock(&_mutex_debug);
#endif
	return;
      }
      if (itmp == jtmp) {
#ifdef DEBUG_EXEC_THREAD
	pthread_mutex_lock(&_mutex_debug);
	{
	  cerr << pid << " pos = " << pos << " jtmp = " << jtmp << endl;
	}
	pthread_mutex_unlock(&_mutex_debug);
#endif
	finish_group = 1;
      }
      // reserve tasks and conunt estimated job size
      pthread_mutex_lock(&_mutex_group[pid_g]);
      {
	for (int i = 0; i < (*it->queue)[pos]->atomic_size; i++) {
	  _group_nops[pid_g][pid0] += 
	    *((*it->queue)[pos + i]->ops_complexity);
	}
      }
      pthread_mutex_unlock(&_mutex_group[pid_g]);
	
#ifdef DEBUG_EXEC_THREAD
      pthread_mutex_lock(&_mutex_debug);
      {
	if ((*it->queue)[pos]->atomic_size > 1) {
	  C_task *kt = (*it->queue)[pos];
	  cerr << pid << " b-c-a " << pos << " : " 
	       << kt->task_name << " : "
	       << kt->parents_work->size() << " / ";
	  for (list<C_task *>::const_iterator mt = kt->parents_work->begin();
	       mt != kt->parents_work->end(); ++mt) {
	    cerr << (int)(*mt)->status << "@" << (*mt)->task_name << " / ";
	  }
	  cerr << "before statically assigned task : with status check." 
	       << endl;
	  }
      }
      pthread_mutex_unlock(&_mutex_debug);
#endif
      for (int i = 0; i < (*it->queue)[pos]->atomic_size; i++) {
	C_task *kt = (*it->queue)[pos + i];
	int waiting1 = check_parents_done(kt, &_mutex_dependency);
	
	while(waiting1 > 0) {
	  //
	  pthread_mutex_lock(&_mutex_root);
	  {
	    _waiting_root++;
#ifdef DEBUG_DEADLOCK
	    if (_waiting_root > num_threads) {
	      fprintf(stderr, "dead lock occured : %s %d\n",
		      __FILE__, __LINE__);
	      exit(-1);
	    }
#endif
#ifdef DEBUG_EXEC_THREAD_IDLE
	    // pthread_mutex_lock(&_mutex_debug);
	    {
	      cerr << pid << " " << kt->task_name << " "
		   << "waiting = "<< waiting1 << " : " 
		   << (int)_waiting_root << " "
		   << __FILE__ << " " << __LINE__ << " ";
	      for (list<C_task *>::const_iterator jt = kt->parents_work->begin();
		   jt != kt->parents_work->end(); ++jt) {
		cerr << (*jt)->task_name << " ";
	      }
	      cerr << " : sleeping" << endl;
	    }
	    // pthread_mutex_unlock(&_mutex_debug);
#endif
	    pthread_cond_wait(&_cond_root, &_mutex_root);
#ifdef DEBUG_EXEC_THREAD_IDLE
	    // pthread_mutex_lock(&_mutex_debug);
	    {
	      cerr << "\t" << pid 
		   << " waked up : _waiting_root = " << (int)_waiting_root 
		     << endl;
	    }
	    // pthread_mutex_unlock(&_mutex_debug);
#endif
	    _waiting_root--;
	    // other thread broacast wake up and decreasing _waiting_root
	    waiting1 = check_parents_done(kt, &_mutex_dependency);
	  }
	  pthread_mutex_unlock(&_mutex_root);
	} // while (waiting > 0)
	if (waiting1 == (-1)) {
#ifdef DEBUG_EXEC_THREAD
	  pthread_mutex_lock(&_mutex_debug);
	  {
	    cerr << pid << " pod = " << pos << " : "
		 << (*it->queue)[pos + i]->task_name		
		 << " : " << it->task_name << " quitted..." << endl;
	  }
	  pthread_mutex_unlock(&_mutex_debug);
#endif
	  return;
	}
	execute_task(it, (pos + i), 
		     permute_block, 
		     pid, &_mutex_dependency);
	if ((*it->queue)[pos + i]->quit_queue) {
#ifdef DEBUG_EXEC_THREAD
	  pthread_mutex_lock(&_mutex_debug);
	  {
	    cerr << pid << " pos = " << pos << " atomic_size= " 
		 << (*it->queue)[pos]->atomic_size
		   << " : " << it->task_name << " finished..." 
		 << "skip " << (*it->queue)[pos + i]->to_next_task
		 << " : " << __FILE__ << " : " << __LINE__
		 << endl;
	    }
	  pthread_mutex_unlock(&_mutex_debug);
#endif
	  pos += (*it->queue)[pos + i]->to_next_task;
	  return;
	  //	  break;
	}
#ifdef DEBUG_EXEC_THREAD
	pthread_mutex_lock(&_mutex_debug);
	{
 	  cerr << "{ " 
	       << (*it->queue)[pos + i]->task_name << " @ "
	       << pid << " " << (pos + i) << " }";
	}
	  pthread_mutex_unlock(&_mutex_debug);
#endif
      } // loop : i
	// static tasks are assigned by other thread
 	
      if (finish_group) {
#ifdef DEBUG_EXEC_THREAD
	pthread_mutex_lock(&_mutex_debug);
	{
	  cerr << pid << " C_FULL : the first greedy end the queue." 
	       << it->pos << endl;
	}
	pthread_mutex_unlock(&_mutex_debug);
#endif
	pthread_mutex_lock(&_mutex_group[pid_g]);
	{
	  pos = it->pos; 
	  for (int p = 0; p < it->num_threads; p++) {
	    _begins_group[pid_g][p][gidxn] = pos;
	    _ends_group[pid_g][p][gidxn] = pos;
	  }
	  _grp_static_assigned[gidxn] = 1;
	}
	pthread_mutex_unlock(&_mutex_group[pid_g]);
	continue;
      }
      pthread_mutex_lock(&_mutex_root);
      {
	zone_flag = (_zone_entered[zone_idxn] == num_threads);
      }
      pthread_mutex_unlock(&_mutex_root);
    } //  while (zone_flag)
#ifdef DEBUG_EXEC_THREAD
    pthread_mutex_lock(&_mutex_debug);
    {
     cerr << pid << " pid_g = " << pid_g << " pid0 " << pid0
	  << " group " << (group - 1)   // after increment of group
	  << " zone_flag = " << zone_flag
	  << endl;
    }
    pthread_mutex_unlock(&_mutex_debug);
#endif
#endif
    if (group_last_entered) {
      // 
      long nops, nops_total, nops_static, nops_per_thread; 
      nops = 0L;
      nops_total = 0L;
      pthread_mutex_lock(&_mutex_group[pid_g]);
      {
	pos = it->pos;   
      }
      pthread_mutex_unlock(&_mutex_group[pid_g]);
      // this value might be changed but
      // no difference on condition queue_size < num_threads
      // pos becomes larger then queue_size smaller 
      if (pos == it->end) {
	break;
      }
      const int queue_size = ((*it->queue)[pos]->parallel_max - 
			      (*it->queue)[pos]->parallel_id);
      if (queue_size < it->num_threads) {
	pthread_mutex_lock(&_mutex_group[pid_g]);
	{
	  pos = it->pos; 
	  for (int p = 0; p < it->num_threads; p++) {
	    _begins_group[pid_g][p][gidxn] = pos;
	    _ends_group[pid_g][p][gidxn] = pos;
	  }
	  _grp_static_assigned[gidxn] = 1;
	}
	pthread_mutex_unlock(&_mutex_group[pid_g]);
      }
      else {
	pthread_mutex_lock(&_mutex_group[pid_g]);
	{
	  pos = it->pos;
	  int jj = pos; 
#if 0
	  for (int i = 0; i < it->num_threads; i++) {
	    nops_total += _group_nops[pid_g][i];
	  }
#endif
	  for (int i = (*it->queue)[pos]->parallel_id;
	       i < (*it->queue)[pos]->parallel_max; i++) {
	    // atomic_size
	    const int itmp = jj;
	    for (int j = 0; j < (*it->queue)[itmp]->atomic_size; 
		 j++, jj++) {
	      nops_total += *((*it->queue)[jj]->ops_complexity);
	    }
	  }
	  nops_static = (long)((double)nops_total * RATIO_QUEUE_GREEDY);
	  nops_per_thread = (nops_static / (long)it->num_threads);
#ifdef DEBUG_EXEC_THREAD
	  pthread_mutex_lock(&_mutex_debug);
	  {
	    cerr << pos << " " <<  _grp_task_ends[gidxn] << " " 
		 << nops_total << " " << nops_per_thread << endl;
	  }
	  pthread_mutex_unlock(&_mutex_debug);
#endif
	  long ntmp = 0L;
	  for (int p = 0; p < it->num_threads; p++) {
	    if (ntmp < _group_nops[pid_g][p]) {
	      ntmp = _group_nops[pid_g][p];  
	    }
	  }
	  if (ntmp >= nops_per_thread) { // greedy
#ifdef DEBUG_EXEC_THREAD
	    pthread_mutex_lock(&_mutex_debug);
	    {
	      cerr << pos << " max _group_nops[pid_g][] =  " << ntmp
		   << " nops_per_threads = " << nops_per_thread << endl;
	    }
	    pthread_mutex_unlock(&_mutex_debug);
#endif
	    for (int p = 0; p < it->num_threads; p++) {
	      _begins_group[pid_g][p][gidxn] = pos;
	      _ends_group[pid_g][p][gidxn] = pos;
	    }
	    //  _group_task_ends[pid_g][cnt_cdfull][gidxn] = end_queue;
	  }
	  else {   // if (ntmp < nops_per_thread)
	    _begins_group[pid_g][0][gidxn] = pos;
	    // for fail-safe
	    //_ends_group[pid_g][it->num_threads - 1][gidxn] = end_queue;
	    _ends_group[pid_g][it->num_threads - 1][gidxn] = 
	      _grp_task_ends[gidxn];
	    nops = 0L;
	    int p = 0;
	    for (int j = pos; j < _grp_task_ends[gidxn]; ) {
	      for (int k = 0; k < (*it->queue)[j]->atomic_size; 
		   k++) {
		nops += *((*it->queue)[j + k]-> ops_complexity);
	      }
	      j += (*it->queue)[j]->atomic_size;
	      // do not devide atomic operation
	      if (nops > nops_per_thread) {
		_ends_group[pid_g][p][gidxn] = j;
		p++;
		if (p == it->num_threads) {
		  break;
		}
		_begins_group[pid_g][p][gidxn] = j;
		nops = 0L;
	      } // if
	    } // loop : j
	    if (p < it->num_threads) {
	      // static assignment is failed.
#ifdef DEBUG_EXEC_THREAD
	      pthread_mutex_lock(&_mutex_debug);
	      {
		cerr << pid << "pid_g = " << pid_g << " gidxn = " << gidxn
		     << " : pos = " << pos << " static assignment is failed"
		     << endl;
	      }
	      pthread_mutex_unlock(&_mutex_debug);
#endif
	      for (int j = 0; j < it->num_threads; j++) {
		_begins_group[pid_g][j][gidxn] = pos;
		_ends_group[pid_g][j][gidxn] = pos;
	      }
	    }
	    it->pos = _ends_group[pid_g][it->num_threads - 1][gidxn];
	    //    cerr << pid << " pos = " << it->pos << endl;
	    //	    _group_task_ends[pid_g][cnt_cdfull][gidxn] = end_queue;
	  } // if (ntmp < nops_per_thread)
	  _grp_static_assigned[gidxn] = 1;
	}
	pthread_mutex_unlock(&_mutex_group[pid_g]);
      } //  if (queue_size < it->num_threads) 
#ifdef DEBUG_EXEC_THREAD
      pthread_mutex_lock(&_mutex_debug);
      {
	cerr << pid << " : ";
	for (int p = 0; p < it->num_threads; p++) {
	  cerr << "[ " << _begins_group[pid_g][p][gidxn] 
	       << " : "  << _ends_group[pid_g][p][gidxn] << " ] ";
	}
	cerr << _grp_task_ends[gidxn] << endl;
      }
      pthread_mutex_unlock(&_mutex_debug);
#endif
    }  // if (group_last_entered)
    else {
      int group_task_end, group_static_assigned;
      pthread_mutex_lock(&_mutex_group[pid_g]);
      {
	group_task_end = _grp_task_ends[gidxn];
	group_static_assigned = _grp_static_assigned[gidxn];
	pos = it->pos;
      }
      pthread_mutex_unlock(&_mutex_group[pid_g]);
#ifdef DEBUG_EXEC_THREAD
      pthread_mutex_lock(&_mutex_debug);
      {
	cerr << pid << " pos = " << pos << " gidxn = " << gidxn
	     << " group_task_end " << group_task_end 
	     << " group_static_assigned " << group_static_assigned
	     << endl;
      }
      pthread_mutex_unlock(&_mutex_debug);
#endif
      if (pos == it->end && (group_static_assigned == 0)) {
	// not yet static assigned but queue is exhausted
#ifdef DEBUG_EXEC_THREAD
	pthread_mutex_lock(&_mutex_debug);
	{
	  cerr << pid << " pod = " << pos << "it->end = " << it->end
	       << " : " << it->task_name << " finished." << endl;
	}
	pthread_mutex_unlock(&_mutex_debug);
#endif
	return;
      }
      // include group_static_assigned == (-1) : 21 Jun.2012, Atsushi
      while (group_static_assigned != 1) { 	      // greedy
	int exit_flag = 0;
	// to skip un-necessary access to _dynamic_queue
	pthread_mutex_lock(&_mutex_group[pid_g]);
	{
	  pos = it->pos;
	  if (pos == it->end) {
	    exit_flag = 1;
	  }
	} 
	pthread_mutex_unlock(&_mutex_group[pid_g]);
	if (exit_flag) {
#ifdef DEBUG_EXEC_THREAD
	  pthread_mutex_lock(&_mutex_debug);
	  {
	    cerr << pid << " pod = " << pos << "it->end = " << it->end
		 << " : " << it->task_name << " already finished." << endl;
	  }
	  pthread_mutex_unlock(&_mutex_debug);
#endif
	  return;
	}
	// looking for dynamic queue
	int flag = 1;
	// include group_static_assigned == (-1) : 21 Jun.2012, Atsushi
	while ((flag == 1) && (group_static_assigned != 1)) {
	  flag = execute_task_dynamic_buffer(permute_block, 
					     pid,
					     &_mutex_dependency);
	  
	  pthread_mutex_lock(&_mutex_group[pid_g]);
	  { 
	    group_static_assigned = _grp_static_assigned[gidxn];
	  }
	  pthread_mutex_unlock(&_mutex_group[pid_g]);
	} // while ((flag == 1) && (group_static_assigned != 0))
	if (group_static_assigned) {
	  break;
	}
	// get tasks with atomic_size from the queue : *it
	int itmp, jtmp;
	pthread_mutex_lock(&_mutex_group[pid_g]);
	{
	  pos = it->pos;
	  // do not increase queue more than its size
	  if (pos == it->end) {
	    exit_flag = 1;
	  }
	  else {
	    // there might be segments with more than 1 automic size
	    // between parallel_id and parallel_max
	    itmp = pos + (*it->queue)[pos]->atomic_size;
	    jtmp = (itmp - 1 + ((*it->queue)[pos]->parallel_max -
				(*it->queue)[pos]->parallel_id));
#ifdef DEBUG_EXEC_THREAD
	    pthread_mutex_lock(&_mutex_debug);
	    {
	      cerr << pid << " pos = " << pos << " +atomic_size = " << itmp
		   << " parallel_max " << jtmp << " : " << it->end 
		   << endl;
	    }
	    pthread_mutex_unlock(&_mutex_debug);
#endif
	    if (itmp <= jtmp) {
	      it->pos = itmp;
	    }
	  }
	} 
	pthread_mutex_unlock(&_mutex_group[pid_g]);
	if (exit_flag) {
#ifdef DEBUG_EXEC_THREAD
	  pthread_mutex_lock(&_mutex_debug);
	  {
	    cerr << pid << " pod = " << pos << "it->end = " << it->end
		 << " : " << it->task_name << " already finished.." << endl;
	  }
	  pthread_mutex_unlock(&_mutex_debug);
#endif
	  return;
	}
	if (itmp == jtmp) {
#ifdef DEBUG_EXEC_THREAD
	  pthread_mutex_lock(&_mutex_debug);
	  {
	    cerr << pid << " pos = " << pos << " jtmp = " << jtmp << endl;
	  }
	  pthread_mutex_unlock(&_mutex_debug);
#endif
	  finish_group = 1;
	}
	// reserve tasks and conunt estimated job size
	pthread_mutex_lock(&_mutex_group[pid_g]);
	{
	  for (int i = 0; i < (*it->queue)[pos]->atomic_size; i++) {
	    _group_nops[pid_g][pid0] += 
	      *((*it->queue)[pos + i]->ops_complexity);
	  }
	}
	pthread_mutex_unlock(&_mutex_group[pid_g]);
	
#ifdef DEBUG_EXEC_THREAD
	pthread_mutex_lock(&_mutex_debug);
	{
	  if ((*it->queue)[pos]->atomic_size > 1) {
	    C_task *kt = (*it->queue)[pos];
	    cerr << pid << " b-c-a " << pos << " : " 
		 << kt->task_name << " : "
		 << kt->parents_work->size() << " / ";
	    for (list<C_task *>::const_iterator mt = kt->parents_work->begin();
		 mt != kt->parents_work->end(); ++mt) {
	      cerr << (int)(*mt)->status << "@" << (*mt)->task_name << " / ";
	    }
	    cerr << "before statically assigned task : with status check." 
		 << endl;
	  }
	}
	pthread_mutex_unlock(&_mutex_debug);
#endif

	for (int i = 0; i < (*it->queue)[pos]->atomic_size; i++) {
	  C_task *kt = (*it->queue)[pos + i];
	  int waiting = check_parents_done(kt, &_mutex_dependency);
	  
	  while(waiting > 0) {
	    //
	    pthread_mutex_lock(&_mutex_root);
	    {
	      _waiting_root++;
#ifdef DEBUG_DEADLOCK
	      if (_waiting_root > num_threads) {
		fprintf(stderr, "dead lock occured : %s %d\n",
			__FILE__, __LINE__);
		exit(-1);
	      }
#endif
#ifdef DEBUG_EXEC_THREAD_IDLE
	      // pthread_mutex_lock(&_mutex_debug);
	      {
		cerr << pid << " " << kt->task_name << " "
		     << "waiting = "<< waiting << " : " 
		     << (int)_waiting_root << " "
		     << __FILE__ << " " << __LINE__ << " ";
		for (list<C_task *>::const_iterator jt = kt->parents_work->begin();
		     jt != kt->parents_work->end(); ++jt) {
		  cerr << (*jt)->task_name << " ";
		}
		cerr << " : sleeping" << endl;
	      }
	      // pthread_mutex_unlock(&_mutex_debug);
#endif
	      pthread_cond_wait(&_cond_root, &_mutex_root);
#ifdef DEBUG_EXEC_THREAD_IDLE
	      // pthread_mutex_lock(&_mutex_debug);
	      {
		cerr << "\t" << pid 
		     << " waked up : _waiting_root = " << (int)_waiting_root 
		     << endl;
	      }
	      // pthread_mutex_unlock(&_mutex_debug);
#endif
	      _waiting_root--;
	      // other thread broacast wake up and decreasing _waiting_root
	      waiting = check_parents_done(kt, &_mutex_dependency);
	    }
	    pthread_mutex_unlock(&_mutex_root);
	  } // while (waiting > 0)
	  if (waiting == (-1)) {
#ifdef DEBUG_EXEC_THREAD
	    pthread_mutex_lock(&_mutex_debug);
	    {
	      cerr << pid << " pod = " << pos + i << " : " 
		   << (*it->queue)[pos + i]->task_name
		   << " : " << it->task_name << " quitted..." << endl;
	    }
	    pthread_mutex_unlock(&_mutex_debug);
#endif
	    return;
	  }
	  execute_task(it, (pos + i), 
		       permute_block, 
		       pid, &_mutex_dependency);
	  if ((*it->queue)[pos + i]->quit_queue) {
#ifdef DEBUG_EXEC_THREAD
	    pthread_mutex_lock(&_mutex_debug);
	    {
	      cerr << pid << " pos = " << pos << " atomic_size= " 
		   << (*it->queue)[pos]->atomic_size
		   << " : " << it->task_name << " finished..." 
		   << "skip " << (*it->queue)[pos + i]->to_next_task
		   << " : " << __FILE__ << " : " << __LINE__
		   << endl;
	    }
	    pthread_mutex_unlock(&_mutex_debug);
#endif
	    pos += (*it->queue)[pos + i]->to_next_task;
	    return;
	    //	    break;
	  }
#ifdef DEBUG_EXEC_THREAD
	  pthread_mutex_lock(&_mutex_debug);
	  {
	    cerr << "{ " 
		 << (*it->queue)[pos + i]->task_name << " @ "
		 << pid << " " << (pos + i) << " }";
	  }
	  pthread_mutex_unlock(&_mutex_debug);
#endif
	} // loop : i
	// static tasks are assigned by other thread
	pthread_mutex_lock(&_mutex_group[pid_g]);
	{ 
	  group_static_assigned = _grp_static_assigned[gidxn];
	  if ((finish_group == 1) && (group_static_assigned == 0)) {
	    // static task is not assigned but queue becomes empty
	    for (int p = 0; p < it->num_threads; p++) {
	      _begins_group[pid_g][p][gidxn] = (-1);
	      _ends_group[pid_g][p][gidxn] = (-1);
	    }
	    _grp_task_ends[gidxn] = jtmp;
	  }
	}
	pthread_mutex_unlock(&_mutex_group[pid_g]);
	
	if (finish_group) {
#ifdef DEBUG_EXEC_THREAD
	  pthread_mutex_lock(&_mutex_debug);
	  {
	    cerr << pid << " C_FULL : the first greedy end the queue." 
		 << it->pos << endl;
	  }
	  pthread_mutex_unlock(&_mutex_debug);
#endif
	  pthread_mutex_lock(&_mutex_group[pid_g]);
	  {
	    pos = it->pos; 
	    for (int p = 0; p < it->num_threads; p++) {
	      _begins_group[pid_g][p][gidxn] = pos;
	      _ends_group[pid_g][p][gidxn] = pos;
	    }
	    _grp_static_assigned[gidxn] = 1;
	  }
	  pthread_mutex_unlock(&_mutex_group[pid_g]);
	  continue;
	}
      // just for debugging : group_task_end >=0 <=> braek while ( < 0 )
#ifdef DEBUG_EXEC_THREAD
	if (group_static_assigned) {
	  pthread_mutex_lock(&_mutex_debug);
	  {
	    cerr << pid << " end the first greedy " 
		 << it->pos << endl;
	  }
	  pthread_mutex_unlock(&_mutex_debug);
	}
#endif
      } // while (group_static_assigned != 1) 
    } //   if (group_last_entered)
    // static assignment
    int _begin_grp, _end_grp;
        
    if ( _grp_static_assigned[gidxn] == (-1)) {
      // this never happen?
      _begin_grp = _end_grp = (-1);
#ifdef DEBUG_EXEC_THREAD_21Jun2012
      pthread_mutex_lock(&_mutex_debug);
      {
	cerr << pid << ": pid0 = " << pid0 << " gidxn = " << gidxn
	     << " static assignment is not yet done! "
	     << endl;
      }
      pthread_mutex_unlock(&_mutex_debug);
#endif
    }
    else {
      _begin_grp = _begins_group[pid_g][pid0][gidxn];
      _end_grp = _ends_group[pid_g][pid0][gidxn];
    }

#ifdef DEBUG_EXEC_THREAD
    pthread_mutex_lock(&_mutex_debug);
    {
      cerr << pid << ": pid0 = " << pid0 << " gidxn = " << gidxn
	   << " < " << _begin_grp << " " << _end_grp << " > assigned" 
	   << endl;
    }
    pthread_mutex_unlock(&_mutex_debug);
#endif
    for (int i = _begin_grp; i < _end_grp; i++) {
      // no need to check dependency for static assignement 
      execute_task(it, i, permute_block, 
		   pid);
      if ((*it->queue)[i]->quit_queue) {
	fprintf(stderr, "%d i = %d : %s finished ... skip %d\n",
		pid, i, it->task_name,
		(*it->queue)[i]->to_next_task);
	i += (*it->queue)[i]->to_next_task;
      }
#ifdef DEBUG_EXEC_THREAD
      pthread_mutex_lock(&_mutex_debug);
      {
	cerr << "( "
	     << (*it->queue)[i]->task_name << " @ "
	     << pid << " " << i << " )";
      }
      pthread_mutex_unlock(&_mutex_debug);
#endif
    } // loop : i
    // greedy
    while(1) {
      int group_task_end;
      pthread_mutex_lock(&_mutex_group[pid_g]);
      {
	pos = it->pos;
	group_task_end = _grp_task_ends[gidxn];
	if (pos < group_task_end) {
	  const int itmp = (it->pos + 
			    (*it->queue)[pos]->atomic_size);
	  if (itmp <= group_task_end) {
	    it->pos = itmp;
	  }
	}
      }
      pthread_mutex_unlock(&_mutex_group[pid_g]);
      if (pos >= group_task_end) {
	break;
      }
      if (_begin_grp < _end_grp) {
      // <==> statically assigned block is computed
#ifdef DEBUG_EXEC_THREAD
	pthread_mutex_lock(&_mutex_debug);
	{
	  if ((*it->queue)[pos]->atomic_size > 1) {
	    C_task *kt = (*it->queue)[pos];
	    cerr << pid << " b-c-a " << pos << " : " << kt->task_name << " : "
		 << kt->parents_work->size() << " / ";
	    
	    for (list<C_task *>::const_iterator mt = kt->parents_work->begin();
		 mt != kt->parents_work->end(); ++mt) {
	      cerr << (int)(*mt)->status << "@" << (*mt)->task_name << " / ";
	    }
	    cerr << "within statically assigned task : without status check." 
		 << endl;
	  }
	}
	pthread_mutex_unlock(&_mutex_debug);
#endif
	for (int i = 0; i < (*it->queue)[pos]->atomic_size; i++) {
	  execute_task(it, (pos + i), permute_block, 
		       pid,
		       &_mutex_group[pid_g]);
	  if ((*it->queue)[pos + i]->quit_queue) {
#ifdef DEBUG_EXEC_THREAD
	    pthread_mutex_lock(&_mutex_debug);
	    {
	      cerr << pid << " pos = " << pos << " atomic_size= " 
		   << (*it->queue)[pos]->atomic_size
		   << " : " << it->task_name << " finished..." 
		   << "skip " << (*it->queue)[pos + i]->to_next_task
		   << endl;
	    }
	    pthread_mutex_unlock(&_mutex_debug);
#endif
	    //	    const int ibegin = pos + i + 1;
	    pos += (*it->queue)[pos + i]->to_next_task;
	    return;
	    // break;
	  }
#ifdef DEBUG_EXEC_THREAD
	  pthread_mutex_lock(&_mutex_debug);
	  {
	    cerr << "( "
		 << (*it->queue)[pos + i]->task_name << " @ "
	       << pid << " " << (pos + i) << " )";
	  }
	  pthread_mutex_unlock(&_mutex_debug);
#endif
	} // loop : i
      } // if (_begins_grp == _ends_grp)
      else {
	// <==> no statically assigned block 
#ifdef DEBUG_EXEC_THREAD
	pthread_mutex_lock(&_mutex_debug);
	{
	  if ((*it->queue)[pos]->atomic_size > 1) {
	    C_task *kt = (*it->queue)[pos];
	    cerr << pid << " b-c-a " << pos << " : " << kt->task_name << " : "
		 << kt->parents_work->size() << " / ";
	    
	    for (list<C_task *>::const_iterator mt = kt->parents_work->begin();
		 mt != kt->parents_work->end(); ++mt) {
	      cerr << (int)(*mt)->status << "@" << (*mt)->task_name << " / ";
	    }
	    cerr << "no statically assigned task : with status check." << endl;
	  }
	}
	pthread_mutex_unlock(&_mutex_debug);
#endif
	for (int i = 0; i < (*it->queue)[pos]->atomic_size; i++) {
	  C_task *kt = (*it->queue)[pos + i];
	  int waiting = check_parents_done(kt, &_mutex_dependency);
	  while(waiting > 0) {
	    pthread_mutex_lock(&_mutex_root);
	    {
	      _waiting_root++;
#ifdef DEBUG_DEADLOCK
	      if (_waiting_root > num_threads) {
		fprintf(stderr, "dead lock occured : %s %d\n",
			__FILE__,  __LINE__);
		exit(-1);
	      }
#endif
#ifdef DEBUG_EXEC_THREAD_IDLE
	      // pthread_mutex_lock(&_mutex_debug);
	      {
		cerr << pid << " " << kt->task_name << " "
		     << "waiting = " << waiting 
		     << " : " << (int)_waiting_root << " : "
		     << __FILE__ << " " << __LINE__ << " ";
		for (list<C_task *>::const_iterator jt = kt->parents_work->begin();
		     jt != kt->parents_work->end(); ++jt) {
		  cerr << (*jt)->task_name << " ";
		}
		cerr << " : sleeping" << endl;
	      }
	      // pthread_mutex_unlock(&_mutex_debug);
#endif
	      pthread_cond_wait(&_cond_root, &_mutex_root);
	      _waiting_root--;
#ifdef DEBUG_EXEC_THREAD_IDLE
	      // pthread_mutex_lock(&_mutex_debug);
	      {
		cerr << endl
		     << pid << " " << kt->task_name << " "
		     << " waked up : _waiting_root = " << (int)_waiting_root 
		     << endl;
	      }
	      // pthread_mutex_unlock(&_mutex_debug);
#endif
	      waiting = check_parents_done(kt, &_mutex_dependency);
	    }
	    pthread_mutex_unlock(&_mutex_root);
	  }
	  if (waiting == (-1)) {
#ifdef DEBUG_EXEC_THREAD
	    pthread_mutex_lock(&_mutex_debug);
	    {
	      cerr << pid << " pod = " << pos + i << " : " 
		   << (*it->queue)[pos + i]->task_name
		   << " : " << it->task_name << " quitted..." << endl;
	    }
	    pthread_mutex_unlock(&_mutex_debug);
#endif
	    return;
	  } // if (wating == (-1))
	  execute_task(it, (pos + i), permute_block, 
		       pid,
		       &_mutex_group[pid_g]);
	  if ((*it->queue)[pos + i]->quit_queue) {
#ifdef DEBUG_EXEC_THREAD
	    pthread_mutex_lock(&_mutex_debug);
	    {
	      cerr << pid << " i = " << i
		   << " : " << it->task_name << " finished... skip " 
		   << (*it->queue)[pos + i]->to_next_task << endl;
	    }
	    pthread_mutex_unlock(&_mutex_debug);
#endif
	    //	    const int ibegin = pos + i + 1;
	    i += (*it->queue)[pos + i]->to_next_task;
	    return;
	  } // if (quit_queue)
#ifdef DEBUG_EXEC_THREAD
	  pthread_mutex_lock(&_mutex_debug);
	  {
	    cerr << "( "
		 << (*it->queue)[pos + i]->task_name << " @ "
		 << pid << " " << (pos + i) << " )";
	  }
	  pthread_mutex_unlock(&_mutex_debug);
#endif
	} // loop : i
      } // // if (_begin_grp <_end_grp)
    } // while (1)
#ifdef DEBUG_EXEC_THREAD
    pthread_mutex_lock(&_mutex_debug);
    {
      cerr << pid << " finish " << group 
	   << " pos = " << it->pos << endl;
    }
    pthread_mutex_unlock(&_mutex_debug);
#endif
    if (pos >= it->end) {
#ifdef DEBUG_EXEC_THREAD
      pthread_mutex_lock(&_mutex_debug);
      {
	cerr << pid << " pod = " << pos << "it->end = " << it->end
	     << " : " << it->task_name << " finished..." << endl;
      }
      pthread_mutex_unlock(&_mutex_debug);
#endif
      return;
    }
    pthread_mutex_lock(&_mutex_group[pid_g]);
    {
      group_first_finished = 0;
      if (_grp_finished[gidxn] == 0) {
	group_first_finished = 1;
      }
      _grp_finished[gidxn]++;
      if (_grp_finished[gidxn] == it->num_threads) {
      // clear status of cyclic buffer
	_grp_finished[gidxp] = 0;
      }
    }
    pthread_mutex_unlock(&_mutex_group[pid_g]);
  } // while (1)
}

void QueueRuntime::thread_queue_parallel_dynamic(const int pid, 
						 const int num_threads,
						 C_task_seq *it,
						 int *permute_block, 
						 const int zone_idxn)
{
  // completely greedy
  int pos, end, atomic_size, atomic_id, update_flag;
  while(1) {
    pthread_mutex_lock(&_mutex_root);
    {
      pos = it->pos;
      end = it->end;
      if (pos < end) {
	atomic_size = (*it->queue)[pos]->atomic_size;
	atomic_id = (*it->queue)[pos]->atomic_id;
	it->pos += atomic_size;  // it->pos++;
	if (it->pos >= end) {
	  update_flag = 1;
	}
      }
    }
    pthread_mutex_unlock(&_mutex_root);
    if (pos >= end) {
      break;
    }
    C_task *kt = (*it->queue)[pos]; // [pos + atomic_size - 1]; 
    int waiting = check_parents_done(kt, &_mutex_dependency);
    while (waiting > 0) {
      // looking for dynamic queue
      int flag = execute_task_dynamic_buffer(permute_block, 
					     pid,
					     &_mutex_dependency);
      if (flag != 1) {
	//
	while(waiting > 0) {
	  pthread_mutex_lock(&_mutex_root);
	  {
	    _waiting_root++;
#ifdef DEBUG_DEADLOCK
	      if (_waiting_root > num_threads) {
		fprintf(stderr, "dead lock occured : %s %d\n",
			__FILE__, __LINE__);
		exit(-1);
	      }
#endif
#ifdef DEBUG_EXEC_THREAD_IDLE
	    // pthread_mutex_lock(&_mutex_debug);
	    {
	      cerr << pid << " " << kt->task_name << " "
		   << "waiting = "<< waiting << " : " 
		   << (int)_waiting_root << " : "
		   << __FILE__ << " " << __LINE__ << " ";
	      for (list<C_task *>::const_iterator jt = kt->parents_work->begin();
		   jt != kt->parents_work->end(); ++jt) {
		cerr << (*jt)->task_name << " ";
	      }
	      cerr << " : sleeping" << endl;
	    }
	    // pthread_mutex_unlock(&_mutex_debug);
#endif
	    pthread_cond_wait(&_cond_root, &_mutex_root);
#ifdef DEBUG_EXEC_THREAD_IDLE
	    // pthread_mutex_lock(&_mutex_debug);
	    {
	      cerr << endl 
		   << pid << " " << kt->task_name << " "
		   << " waked up : _waiting_root = " << (int)_waiting_root 
		   << endl;
	    }
	    // pthread_mutex_unlock(&_mutex_debug);
#endif
	    _waiting_root--;
	    // other thread broacast wake up and decreasing _waiting_root
	  }
	  pthread_mutex_unlock(&_mutex_root);
	  waiting = check_parents_done(kt, &_mutex_dependency);
	} // while (waiting > 0)
      } // if (flag != 1)
      waiting = check_parents_done(kt, &_mutex_dependency);
    } // while (waiting > 0)
    for (int m = 0; m < atomic_size; m++) {
      execute_task(it, (pos + m), permute_block,  // loop with atomic_size
		   pid,
		   &_mutex_dependency);
    }
    if ((*it->queue)[pos]->quit_queue) {
#ifdef DEBUG_EXEC_THREAD
      pthread_mutex_lock(&_mutex_debug);
      {
	cerr << pid << " pos = " << pos
		   << " : " << it->task_name << " finished..." 
	     << "skip " << (*it->queue)[pos]->to_next_task
		   << endl;
      }
      pthread_mutex_unlock(&_mutex_debug);
#endif
      //      const int ibegin = pos + 1;
      pos += (*it->queue)[pos]->to_next_task;
    } // if (quit_queue)
#ifdef DEBUG_EXEC_THREAD
    pthread_mutex_lock(&_mutex_debug);
    {
      cerr << "( " << pid << " " << pos << " ) " ;
    }
    pthread_mutex_unlock(&_mutex_debug);
#endif
  } // while (1)
#ifdef DEBUG_EXEC_THREAD
  pthread_mutex_lock(&_mutex_debug);
  {
    cerr << pid << " " << it->task_name 
	 << " @ " << it->task_id << " greedy end." << endl;
  }
  pthread_mutex_unlock(&_mutex_debug);
#endif
}

void QueueRuntime::thread_queue_parallel_static(const int pid, 
						const int num_threads,
						C_task_seq* it, 
						int *permute_block, 
						const int zone_idxp,
						const int zone_idxn,
						const int zone_first_entered,
						const int zone_last_entered)
{
  // the last thread which turns into the zone assigns jobs statiscally
  int pos, end, atomic_size, atomic_id;
  int finished_status = 0;
#if 0 // 12 Feb.2013 Atsushi --> move to just after zone_first_entered = 1
  if (zone_first_entered) {
    pthread_mutex_lock(&_mutex_root);
    {
      _zone_static_assigned[zone_idxn] = 0;
    }
    pthread_mutex_unlock(&_mutex_root);
  } 
#endif
  pthread_mutex_lock(&_mutex_root);
  {
    pos = it->pos;
  }
  pthread_mutex_unlock(&_mutex_root);

#ifdef DEBUG_EXEC_THREAD
  pthread_mutex_lock(&_mutex_debug);
  {
    cerr << pid << "_phase_dynamic = " << _phase_dynamic << " @ "
	 << _queue_dynamic->size() << endl;
  }
  pthread_mutex_unlock(&_mutex_debug);
#endif
  // reach the end before all threads entered
  if (pos == it->end) {
    return;
  }
  if (zone_last_entered) {
    long nops, nops_total, nops_static, nops_per_thread; 
    int begin, end;
    // this mutex is relatively larage
    //  : lock out other threads until assingnment is finished
    
    pthread_mutex_lock(&_mutex_root);
    {
#if 0 // 12 Feb.2013 Atsushi --> move to just after zone_last_entered = 1 
      _zone_static_assigned[zone_idxp] = 0;
#endif
      begin = it->pos;
      end = it->end;
      if ((end - begin) < num_threads) {
	// execute tasks by greedy way
	for (int q = 0; q < num_threads; q++) {
	    _begins[q][zone_idxn] = _ends[q][zone_idxn] = pos;
	}
	_zone_static_assigned[zone_idxn] = 1;
      }
      else {
      //	      nops_total = it->ops_complexity;
	nops_total = 0L;
	for (int j = begin; j < end; j++) {
	  nops_total += *((*it->queue)[j]->ops_complexity);
	  // 28 May 2012 : Atsushi
	  //	  cerr << (*it->queue)[j]->task_name
	  //   << *((*it->queue)[j]->ops_complexity) << " ";
	}
	//cerr << endl;

	nops_static = (long)((double)nops_total * RATIO_QUEUE_GREEDY);
	nops_per_thread = nops_static / (long)num_threads;
#ifdef DEBUG_EXEC_THREAD
	pthread_mutex_lock(&_mutex_debug);
	{
	  cerr << pid << " begin = " << begin << " end = " 
	       << end << " "
	       << nops_total << " " << nops_static << " " 
	       << nops_per_thread << endl;
	}
	pthread_mutex_unlock(&_mutex_debug);
#endif
	_begins[0][zone_idxn] = begin;
	_ends[num_threads - 1][zone_idxn] = end;
	nops = 0L;
	int p = 0;
	for (int j = begin; j < end; j++) {
	  nops += *((*it->queue)[j]->ops_complexity);
	  if (nops > nops_per_thread) {
	    //	     &&  ((*it->queue)[j]->atomic_size ==   // 13 May 2014
	    //	       ((*it->queue)[j]->atomic_id + 1))) {
	    const int jj = j - (*it->queue)[j]->atomic_id;
	    _ends[p][zone_idxn] = jj;
	    p++;
	    if (p == num_threads) {
	      break;
	    }
	    _begins[p][zone_idxn] = jj;
	    nops = 0L;
	  }
	} // loop : j
	if (p < num_threads) { // equally assignment is failed
	  // execute tasks by greedy way
	  for (int q = 0; q < num_threads; q++) {
	    _begins[q][zone_idxn] = _ends[q][zone_idxn] = pos;
	  }
	}
	else {
	  it->pos = _ends[num_threads - 1][zone_idxn];
	}
	_zone_static_assigned[zone_idxn] = 1;
      }
    } //  if ((end - begin) < num_threads) 
    pthread_mutex_unlock(&_mutex_root);
#ifdef DEBUG_EXEC_THREAD
    pthread_mutex_lock(&_mutex_debug);
    {
      for (int p = 0; p < num_threads; p++) {
	cerr << "[ " << _begins[p][zone_idxn] 
	     << " : "  << _ends[p][zone_idxn] << " ] ";
      }
      cerr << it->pos << endl;
    }
    pthread_mutex_unlock(&_mutex_debug);
#endif
    //
  }  // if (zone_last_entered)
  else { // if (zone_last_entered)
    int zone_static_assigned = 0;
#if 0
    pthread_mutex_lock(&_mutex_root);
    {
      zone_static_assigned = _zone_static_assigned[zone_idxn];
    }
    pthread_mutex_unlock(&_mutex_root);
#ifdef DEBUG_EXEC_THREAD
    pthread_mutex_lock(&_mutex_debug);
    {
      cerr << pid << " pos = " << pos << " zone_idxn = " << zone_idxn
	   << " zone_static_assigned " << zone_static_assigned << " : "
	   << _zone_static_assigned[zone_idxn]
	   << endl;
    }
    pthread_mutex_unlock(&_mutex_debug);
#endif
#endif
    while(zone_static_assigned == 0) {
      // similar to the routine : execute_task_dynamic_buffer()
      pthread_mutex_lock(&_mutex_root);
      {
	zone_static_assigned = _zone_static_assigned[zone_idxn];
	if (zone_static_assigned == 0) {
	  pos = it->pos;
	  end = it->end;
	  if (pos < end) {
	    atomic_size = (*it->queue)[pos]->atomic_size;
	    atomic_id = (*it->queue)[pos]->atomic_id;
	    it->pos += atomic_size;                      // pos++;
	    if (it->pos >= end) {
	      finished_status = 1;
	      //
	      _phase_dynamic++;
	      if (_phase_dynamic == _queue_dynamic->size()) {
		_phase_dynamic = (-1);
	      } 
	// queue is exhaused by this step and enforcing other tasks to quit
	      for (int p = 0; p < num_threads; p++) {
		_begins[p][zone_idxn] = _ends[p][zone_idxn] = end;
	      }
	    } // if (it->pos == end) {
	  }
	} // if (zone_static_assigned == 0) 
      }
#ifdef DEBUG_EXEC_THREAD
      pthread_mutex_lock(&_mutex_debug);
      {
	cerr << pid << " pos = " << pos << " zone_idxn = " << zone_idxn
	     << " zone_static_assigned " << zone_static_assigned << " : "
	     << _zone_static_assigned[zone_idxn]
	     << endl;
      }
      pthread_mutex_unlock(&_mutex_debug);
#endif
      pthread_mutex_unlock(&_mutex_root);
      if (zone_static_assigned == 1) {
	break;
      }
      if (pos >= end) { 
	// skip my task
	_begins[pid][zone_idxn] = _ends[pid][zone_idxn] = end;

	//	pthread_mutex_unlock(&_mutex_root);  // bug : 25 Jun.2012
#ifdef DEBUG_EXEC_THREAD
	pthread_mutex_lock(&_mutex_debug);
	{
	  cerr << pid << " queue is already exhausted " 
	       << it->pos << endl;
	}
	pthread_mutex_unlock(&_mutex_debug);
#endif
	break;
      }
      C_task *kt = (*it->queue)[pos]; // [pos + atomic_size - 1];   
      int waiting = check_parents_done(kt, &_mutex_dependency);
      if (waiting > 0) {
	pthread_mutex_lock(&_mutex_root);
	{
	  while (waiting > 0) { // real sleeping
	    _waiting_root++;
#ifdef DEBUG_DEADLOCK
	      if (_waiting_root > num_threads) {
		fprintf(stderr, "dead lock occured : %s %d\n",
			__FILE__,  __LINE__);
		exit(-1);
	      }
#endif
#ifdef DEBUG_EXEC_THREAD_IDLE
	      // pthread_mutex_lock(&_mutex_debug);
	    {
	      cerr << pid << " " << kt->task_name << " "
		   << "waiting = "<< waiting << " : " 
		   << (int)_waiting_root << " : " 
		   << __FILE__ << " " << __LINE__ << " ";
	      for (list<C_task *>::const_iterator jt = kt->parents_work->begin();
		   jt != kt->parents_work->end(); ++jt) {
		cerr << (*jt)->task_name << " ";
	      }
	      cerr << " : sleeping" << endl;
	    }
	    // pthread_mutex_unlock(&_mutex_debug);
#endif
	    //
	    pthread_cond_wait(&_cond_root, &_mutex_root);
#ifdef DEBUG_EXEC_THREAD_IDLE
	    // pthread_mutex_lock(&_mutex_debug);
	    {
	      cerr << endl
		   << pid << " " << kt->task_name << " "
		   << " waked up : _waiting_root = " << (int)_waiting_root 
		   << endl;
	    }
	    // pthread_mutex_unlock(&_mutex_debug);
#endif
	    _waiting_root--;
	    // other thread broacast wake up and decreasing _waiting_root
	    waiting = check_parents_done(kt, &_mutex_dependency);
	  } // while (waiting > 0)
	}
	pthread_mutex_unlock(&_mutex_root);
      } // if (waiting > 0)
      atomic_size = (*it->queue)[pos]->atomic_size;
      for (int m = 0; m < atomic_size; m++) {
	execute_task(it, (pos + m), permute_block,  // loop with atomic_size
		     pid, &_mutex_dependency);
      }
#ifdef DEBUG_EXEC_THREAD
      pthread_mutex_lock(&_mutex_debug);
      {
	cerr << "{ " << pid << " " << pos << " }";
      }
      pthread_mutex_unlock(&_mutex_debug);
#endif
      if (finished_status) { //  && (zone_static_assigned == 0)) {
#if 0
	// enforcing other tasks to quit
	pthread_mutex_lock(&_mutex_root);
	{
	  for (int p = 0; p < num_threads; p++) {
	    _begins[p][zone_idxn] = _ends[p][zone_idxn] = end;
	  }
	}
	pthread_mutex_unlock(&_mutex_root);
#endif
#ifdef DEBUG_EXEC_THREAD
	pthread_mutex_lock(&_mutex_debug);
	{
	  cerr << pid << " queue is exhausted in the first greedy" 
	       << it->pos << endl;
	}
	pthread_mutex_unlock(&_mutex_debug);
#endif
	break;
      }
      pthread_mutex_lock(&_mutex_root);
      {
	it->ops_complexity -= *((*it->queue)[pos]->ops_complexity);
	zone_static_assigned = _zone_static_assigned[zone_idxn];
      }
      pthread_mutex_unlock(&_mutex_root);
      if (zone_static_assigned == 1) {
#ifdef DEBUG_EXEC_THREAD
	pthread_mutex_lock(&_mutex_debug);
	{
	  cerr << pid << " end the first greedy " 
	       << it->pos << endl;
	}
	pthread_mutex_unlock(&_mutex_debug);
#endif
      } 
    } // while (zone_static_assigned == 0)
  } //    if (zone_last_entered)
  // static
#ifdef DEBUG_EXEC_THREAD
  pthread_mutex_lock(&_mutex_debug);
  {
    cerr << "pid = " << pid << " zone_idxn = " << zone_idxn 
	 << "< " << _begins[pid][zone_idxn]
	 << " " << _ends[pid][zone_idxn] << " > assigned" << endl;
  }
  pthread_mutex_unlock(&_mutex_debug);
#endif

  for (int i = _begins[pid][zone_idxn]; i < _ends[pid][zone_idxn]; i++) {
    execute_task(it, i, permute_block, 
		 pid);
  }  
#ifdef DEBUG_EXEC_THREAD
  pthread_mutex_lock(&_mutex_debug);
  {
    cerr << pid << ": " << "< " << _begins[pid][zone_idxn]
	 << " " << _ends[pid][zone_idxn] << " > done." << endl;
  }
  pthread_mutex_unlock(&_mutex_debug);
#endif
 // greedy
#ifdef DEBUG_EXEC_THREAD
  pthread_mutex_lock(&_mutex_debug);
  {
    cerr << pid << " " 
	 << it->task_name << " @ " << it->task_id 
	 << " greedy begin: " << endl;
  }
  pthread_mutex_unlock(&_mutex_debug);
#endif
  while(1) {
    pthread_mutex_lock(&_mutex_root);
    {
      pos = it->pos;
      end = it->end;
      if (pos < end) {
	atomic_size = (*it->queue)[pos]->atomic_size;
	atomic_id = (*it->queue)[pos]->atomic_id;
	it->pos += atomic_size;                      // pos++
	if (it->pos >= end) {
	  finished_status = 1;
	  _phase_dynamic++;
	  if (_phase_dynamic == _queue_dynamic->size()) {
	    _phase_dynamic = (-1);
	  } // queue is exhausted
#ifdef DEBUG_EXEC_THREAD
	  pthread_mutex_lock(&_mutex_debug);
	  {
	    cerr << pid << " _phase_dynamic = " 
		 << _phase_dynamic << " @ " 
		 << _queue_dynamic->size() << endl;
	    if (_phase_dynamic >= 0) {
	      cerr << pid << " next =  " 
		   << (*_queue_dynamic)[_phase_dynamic]->task_name 
		   << endl;
	    }
	  }
	  pthread_mutex_unlock(&_mutex_debug);
#endif
	}
      } // if (pos < end)
    }
    pthread_mutex_unlock(&_mutex_root);
    if (pos >= end) {
#ifdef DEBUG_EXEC_THREAD
      pthread_mutex_lock(&_mutex_debug);
      {
	cerr << pid << " queue is already exhausted " << endl;
      }
      pthread_mutex_unlock(&_mutex_debug);
#endif
      break;
    }
    atomic_size = (*it->queue)[pos]->atomic_size;
    for (int m = 0; m < atomic_size; m++) {
      execute_task(it, (pos + m), permute_block,    // loop with atomic_size
		   pid,
		   &_mutex_dependency);
    }
#ifdef DEBUG_EXEC_THREAD1
    pthread_mutex_lock(&_mutex_debug);
    {
      cerr << "( " << pid << " " << pos << " : " 
	   << (*it->queue)[pos]->task_name << " ) ";
    }
    pthread_mutex_unlock(&_mutex_debug);
#endif	      

    if (finished_status) {
      break;
    }
  } // while(1)
  //
#ifdef DEBUG_EXEC_THREAD
  pthread_mutex_lock(&_mutex_debug);
  {
    cerr << pid << " greedy end : " << endl;
  }
  pthread_mutex_unlock(&_mutex_debug);
#endif
}

void QueueRuntime::thread_queue_single(const int pid, 
				       const int num_threads,
				       C_task_seq *it, 
				       int *permute_block, 
				       const int zone_idxn)
{
// tasks are excuted by a single thread : atomic operation will not be divided.
  int start = it->begin;
  bool zone_flag;
  pthread_mutex_lock(&_mutex_root);
  {
    zone_flag = (_zone_entered[zone_idxn] == num_threads);
  }
  pthread_mutex_unlock(&_mutex_root);

#ifdef DEBUG_EXEC_THREAD
  pthread_mutex_lock(&_mutex_debug);
  {
    cerr << pid << " @ " << num_threads << " single " 
       << it->task_name << " : " << it->task_id 
       << " zone_flag = " << zone_flag << endl;
  }
  pthread_mutex_unlock(&_mutex_debug);
#endif

  while ((start < it->end) && (!zone_flag)) {
    C_task *kt = (*it->queue)[start];
    int waiting = check_parents_done(kt, &_mutex_dependency);
    while(waiting > 0) {
      // looking for dynamic queue
      int flag = execute_task_dynamic_buffer(permute_block, 
					     pid,
					     &_mutex_dependency);
      if (flag != 1) {
	//
	int waiting2;
	waiting2 = check_parents_done(kt, &_mutex_dependency);
	while(waiting2 > 0) {
	  pthread_mutex_lock(&_mutex_root);
	  {
	    _waiting_root++;
#ifdef DEBUG_DEADLOCK
	    if (_waiting_root > num_threads) {
	      fprintf(stderr, "dead lock occured : %s %d\n",
		      __FILE__, __LINE__);
	      exit(-1);
	    }
#endif
#ifdef DEBUG_EXEC_THREAD_IDLE
	    //	    pthread_mutex_lock(&_mutex_debug);
	    {
	      cerr << pid << " " << kt->task_name << " "
		   << "waiting = "<< waiting << " : " 
		   << (int)_waiting_root << " : "
		   << __FILE__ << " " << __LINE__ << " ";
	      for (list<C_task *>::const_iterator jt = kt->parents_work->begin();
		   jt != kt->parents_work->end(); ++jt) {
		cerr << (*jt)->task_name << " ";
	      }
	      cerr << " : sleeping" << endl;
	    }
	    //pthread_mutex_unlock(&_mutex_debug);
#endif
	    pthread_cond_wait(&_cond_root, &_mutex_root);
#ifdef DEBUG_EXEC_THREAD_IDLE
	    //	    pthread_mutex_lock(&_mutex_debug);
	    {
	      cerr << endl
		   << pid << " " << kt->task_name << " "
		   << " waked up : _waiting_root = " << (int)_waiting_root
		   << endl;
	    }
	    //  pthread_mutex_unlock(&_mutex_debug);
#endif
	    _waiting_root--;
	    // other thread broacast wake up and decreasing _waiting_root
	    waiting2 = check_parents_done(kt, &_mutex_dependency);
	  }
	  pthread_mutex_unlock(&_mutex_root);
	} // while (waiting2)
      } // if (flag != 1)
      waiting = check_parents_done(kt, &_mutex_dependency);
      pthread_mutex_lock(&_mutex_root);
      {
	zone_flag = (_zone_entered[zone_idxn] == num_threads); 
      }
      pthread_mutex_unlock(&_mutex_root);
      if (zone_flag) {
	break;
      }
    } // while (waiting > 0)
    execute_task(it, start, permute_block, 
		 pid,
		 &_mutex_dependency);
    if ((*it->queue)[start]->quit_queue) {
#ifdef DEBUG_EXEC_THREAD
      pthread_mutex_lock(&_mutex_debug);
      {
	cerr << pid << " start = " << start
	     << " : " << it->task_name << " finished... skip " 
	     << (*it->queue)[start]->to_next_task << endl;
      }
      pthread_mutex_unlock(&_mutex_debug);
#endif
      //      const int ibegin = start + 1;
      start += (*it->queue)[start]->to_next_task;
   }
#ifdef DEBUG_EXEC_THREAD      
    pthread_mutex_lock(&_mutex_debug);
    {
      cerr << "( " << pid << " : " 
	   << (*it->queue)[start]->task_name << " ) ";
    }
    pthread_mutex_unlock(&_mutex_debug);
#endif
    start++;
    pthread_mutex_lock(&_mutex_root);
    {
      zone_flag = (_zone_entered[zone_idxn] == num_threads);     
    }
    pthread_mutex_unlock(&_mutex_root);
  } // while ((start < it->end) || (!zone_flag))
  for (int j = start; j < it->end; j++) {
    execute_task(it, j, permute_block, 
		 pid);
    if ((*it->queue)[j]->quit_queue) {
#ifdef DEBUG_EXEC_THREAD
      pthread_mutex_lock(&_mutex_debug);
      {
	cerr << pid << " j = " << j
	     << " : " << it->task_name << " finished... skip " 
	     << (*it->queue)[j]->to_next_task << endl;
      }
      pthread_mutex_unlock(&_mutex_debug);
#endif
      //      const int ibegin = j + 1;
      j += (*it->queue)[j]->to_next_task;
    }
  } // loop : j

#ifdef DEBUG_EXEC_THREAD2
  pthread_mutex_lock(&_mutex_debug);
  {
    for (int j = start; j < it->end; j++) {
      cerr << "( " << pid << " : " 
	   << (*it->queue)[j]->task_name << " ) ";
    }
  }
  pthread_mutex_unlock(&_mutex_debug);
#endif
  //
}
// #define TASKSTATUS_MUTEX
#ifdef TASKSTATUS_MUTEX
inline void QueueRuntime::taskstatus_mutex_lock(pthread_mutex_t *mutex_dependency)
{
  pthread_mutex_lock(mutex_dependency);
}

inline void QueueRuntime::taskstatus_mutex_unlock(pthread_mutex_t *mutex_dependency)
{
  pthread_mutex_unlock(mutex_dependency);
}
#else
inline void QueueRuntime::taskstatus_mutex_lock(pthread_mutex_t *mutex_dependency)
{
  //dummy
}

inline void QueueRuntime::taskstatus_mutex_unlock(pthread_mutex_t *mutex_dependency)
{
  //dummy
}
#endif

int QueueRuntime::execute_task_dynamic_buffer(int *permute_block, 
						 int pid,
						 pthread_mutex_t *mutex_dependency)
{
  int ph_dynmc, pos, end, atomic_size, atomic_id;
  int waiting;
  int update_status = 0;
  C_task *kt;

  pthread_mutex_lock(&_mutex_root);
  {
    ph_dynmc = _phase_dynamic;
    if (ph_dynmc >= 0) {
      pos = (*_queue_dynamic)[ph_dynmc]->pos;
      end = (*_queue_dynamic)[ph_dynmc]->end;
      if (pos < end) {
	atomic_size = (*(*_queue_dynamic)[ph_dynmc]->queue)[pos]->atomic_size;
	atomic_id = (*(*_queue_dynamic)[ph_dynmc]->queue)[pos]->atomic_id;
#if 0
	if (atomic_id == (atomic_size - 1)) {
	  atomic_size = 1;
	}
#endif
	kt = (*(*_queue_dynamic)[ph_dynmc]->queue)[pos];
	// to avoid other thread takes the same pos
	waiting = check_parents_done(kt, mutex_dependency);
	if (waiting <= 0) {
	  (*_queue_dynamic)[ph_dynmc]->pos += atomic_size;    // pos++
	  if ((*_queue_dynamic)[ph_dynmc]->pos == end) {
	    update_status = 1;
	    _phase_dynamic++;
	    if (_phase_dynamic == _queue_dynamic->size()) {
	      _phase_dynamic = (-1);
	    }
	  } // if ((*_queue_dynamic)[ph_dynmc]->pos == end)
	}
      }
      else {
	ph_dynmc = (-1);
      }
    } // if (ph_dynmc > 0)
  }
  pthread_mutex_unlock(&_mutex_root);

  if (ph_dynmc < 0) {
    return (-1);  // queue exhausted
  }
  if (waiting > 0) {  // waiting == (-1) <=> quit_queue == true
    return 0;     // the first task in the queue is not ready
  }

#ifdef DEBUG_EXEC_THREAD
    pthread_mutex_lock(&_mutex_debug);
    {
      cerr << pid << " from dynamic " 
	   << (*_queue_dynamic)[ph_dynmc]->task_name << "  "
	   << pos << " @ " << ph_dynmc
	   << " " << (*_queue_dynamic)[ph_dynmc]->end
	   << " " << kt->task_name << " : ";
    }
    pthread_mutex_unlock(&_mutex_debug);
#endif
    for (int m = 0; m < atomic_size; m++) {
      execute_task((*_queue_dynamic)[ph_dynmc], (pos + m),  // pos
		   permute_block, 
		   pid, mutex_dependency);
    }

#ifdef DEBUG_EXEC_THREAD
  pthread_mutex_lock(&_mutex_debug);
  {
    cerr << kt->task_name << " : " 
	 << pos << " done" << endl;
  }
  pthread_mutex_unlock(&_mutex_debug);
#endif
  // ops_complexity is recalculated in each routine :-
  pthread_mutex_lock(&_mutex_root);
  {
    (*_queue_dynamic)[ph_dynmc]->ops_complexity -= *(kt->ops_complexity);
    if (update_status) {
      (*_queue_dynamic)[ph_dynmc]->status = TASK_DONE;
    }
  }
  pthread_mutex_unlock(&_mutex_root);

  return 1;              // success 
}

int QueueRuntime::check_parents_done(C_task *it, 
				     pthread_mutex_t *mutex_dependency)
{
  unsigned char status;
  int dependency = 0;
  bool quit_queue = false;
  //  cerr << it->task_name << " size = " << it->parents_work->size() << " || ";
  for(list<C_task *>::iterator nt = it->parents_work->begin();
      nt != it->parents_work->end(); ) {
    //    cerr << (*nt)->task_name << " ";
    taskstatus_mutex_lock(mutex_dependency); 
    {
      status = (*nt)->status;
    }
    taskstatus_mutex_unlock(mutex_dependency); 
    if (status == TASK_DONE) {
      if ((*nt)->quit_queue) {
	quit_queue = true;
      }
      // 5 Jul.2015 Atsushi
      nt = it->parents_work->erase(nt);
      // ++nt;
    }
    else {
      dependency++;
      ++nt;
    }
  }
  //  cerr << " >> " << it->parents_work->size() << " . " << endl;
#if 0 // 8 Jul.2015 Atsushi
  if (it->parents_work->size() > 0) {
    cerr << it->task_name << " size = " << it->parents_work->size() << " | ";
    for(list<C_task *>::const_iterator nt = it->parents_work->begin();
	nt != it->parents_work->end(); ++nt) {
      cerr << (*nt)->task_name << " ";
    }
    cerr << "\n";
  }
#endif
  if (quit_queue && (dependency == 0)) {
    dependency = (-1);
  }
  return dependency;
}

int QueueRuntime::check_parents_done(C_task *it)
{
  int dependency = 0;
  bool quit_queue = false;
  //  cerr << it->task_name << " size = " << it->parents_work->size() << " | ";
  for(list<C_task *>::iterator nt = it->parents_work->begin();
      nt != it->parents_work->end(); ) {
    //    cerr << (*nt)->task_name << " ";
    if ((*nt)->status == TASK_DONE) {
      if ((*nt)->quit_queue) {
	quit_queue = true;
      }
      // 5 Jul.2015 Atsushi
      nt = it->parents_work->erase(nt);
      // ++nt;
    }
    else {
      dependency++;
      ++nt;
    }
  }
#if 0 // 8 Jul.2015 Atsushi
  //  cerr << " >> " << it->parents_work->size() << " . " << endl;
  if (it->parents_work->size() > 0) {
    cerr << it->task_name << " size = " << it->parents_work->size() << " | ";

    for(list<C_task *>::const_iterator nt = it->parents->begin();
	nt != it->parents->end(); ++nt) {
      cerr << (*nt)->task_name << " ";
    }
    cerr << " || ";
    
    for(list<C_task *>::const_iterator nt = it->parents_work->begin();
	nt != it->parents_work->end(); ++nt) {
      cerr << (*nt)->task_name << " ";
    }
    cerr << "\n";
  }
#endif
  if (quit_queue && (dependency == 0)) {
    dependency = (-1);
  }
  return dependency;
}

void QueueRuntime::execute_task(C_task_seq *seq, int pos, 
				int *permute_block, 
				int pid,
				pthread_mutex_t *mutex_dependency)
{
  C_task *task = (*seq->queue)[pos];
  // debugging
#ifdef DEBUG_CHECKPARENTS_DONE
  int waiting;
  waiting = check_parents_done(task, mutex_dependency);
  if (waiting > 0) {

    pthread_mutex_lock(&_mutex_debug);
    {
      cerr << pid << " parents of task " << task->task_name 
	   << " not finished : ";
      for(list<C_task *>::const_iterator nt = task->parents_work->begin();
	  nt != task->parents_work->end(); ++nt) {
	cerr << (*nt)->task_name << " ";
      }
      cerr << endl;
    }
    pthread_mutex_unlock(&_mutex_debug);
  }
#endif
  //
#ifdef DEBUG_PRINT_TASK
  switch (task->task_id) {
  case C_FILLMATRIX:
  case C_DFULL_SYM_GAUSS:
  case C_DTRSMSCALE:
    //  case C_SPARSESCHUR:
  case C_DGEMM_LOCAL_MULT:
  case C_DGEMM_LOCAL_TWO:
  case C_DGEMM_DIRECT_TWO:
    fprintf(_fps[pid], "%s\n", task->task_name);
    *(task->fp) = &_fps[pid];
    break;
  }
#endif
  taskstatus_mutex_lock(mutex_dependency);
  {
    task->status = TASK_WORKING;
  }
  taskstatus_mutex_unlock(mutex_dependency);
#ifdef DEBUG_THREAD_TIME
  get_realtime(&(task->t0));
  //  clock_gettime(CLOCK_MONOTONIC, &(task->t0));
#endif
  task->func(task->func_arg);
#ifdef DEBUG_THREAD_DONE_PRINT
  pthread_mutex_lock(&_mutex_debug);
  {
     cerr << "pid = " << pid << " : " << task->task_name << endl;
  }
  pthread_mutex_unlock(&_mutex_debug);
#endif
  if (task->task_id == C_DFULL_SYM_GAUSS) {
    // task->quit_queue = ((C_dfull_gauss_arg*)task->func_arg)->quit;
    taskstatus_mutex_lock(mutex_dependency);
    {
      for (int i = (pos + 1); i <= (pos + task->to_next_task); i++) {
	(*seq->queue)[i]->status = TASK_DONE;
	(*seq->queue)[i]->quit_queue = true;
#ifdef DEBUG_EXEC_THREAD2
	pthread_mutex_lock(&_mutex_debug);
	{
	  cerr << (*seq->queue)[i]->task_name << " ";
	}
	pthread_mutex_unlock(&_mutex_debug);
#endif
      }
    }
    taskstatus_mutex_unlock(mutex_dependency);
    //    task->to_next_task = ((C_dfull_gauss_arg*)task->func_arg)->to_next_task;
  }
#ifdef DEBUG_THREAD_TIME
    get_realtime(&(task->t1));
  //  clock_gettime(CLOCK_MONOTONIC, &(task->t1));
  task->thread_id = pid;
#endif
#ifdef DEBUG_EXEC_THREAD_FILE
  pthread_mutex_lock(&_mutex_file);
  {
    const int sec_n0 = convert_sec(task->t0); 
    const int sec_m0 = convert_microsec(task->t0);
    _fout << task->task_name << " / ";
    _fout << sec_n0 << " " << sec_m0 << "\t";
    for(list<C_task *>::const_iterator nt = task->parents->begin();
	nt != task->parents->end(); ++nt) {
      int sec_n1 = convert_sec((*nt)->t1);
      int sec_m1 = convert_microsec((*nt)->t1);
      if ((sec_n0 < sec_n1) || ((sec_n0 == sec_n1) && (sec_m0 < sec_m1))) {
	_fout << "*";
      }
      _fout << (*nt)->task_name << " / ";
      _fout << sec_n1 << " " << sec_m1 << "\t";
    }
    _fout << endl; // " : " << pid << endl;
  }
  pthread_mutex_unlock(&_mutex_file);
#endif
  taskstatus_mutex_lock(mutex_dependency);
  {
    task->status = TASK_DONE;
  }
  taskstatus_mutex_unlock(mutex_dependency);
  //
  if (_waiting_root > 0) {
    pthread_mutex_lock(&_mutex_root);
    {
      pthread_cond_broadcast(&_cond_root);
    }
    pthread_mutex_unlock(&_mutex_root);
#ifdef DEBUG_EXEC_THREAD_IDLE
    //    pthread_mutex_lock(&_mutex_debug);
    {
      cerr << pid << " " << task->task_name 
	   << "_waiting_root = " << (int)_waiting_root 
	   << " broadcast." << endl;
    }
    //    pthread_mutex_unlock(&_mutex_debug);
#endif
  }

}

void QueueRuntime::execute_task(C_task_seq *seq, int pos, 
				int *permute_block, 
				int pid)
{
  C_task *task = (*seq->queue)[pos];
#ifdef DEBUG_CHECKPARENTS_DONE
 // debugging
  int waiting;
  waiting = check_parents_done(task);
  if (waiting > 0) {
    pthread_mutex_lock(&_mutex_debug);
    {
      cerr << pid << " parents of task " << task->task_name 
	   << " not finished : " << task->parents_work->size() << " : ";
      // debug : 5 Jul. 2015, Atsushi
      for(list<C_task *>::const_iterator nt = task->parents_work->begin();
	  nt != task->parents_work->end(); ++nt) {
	cerr << (*nt)->task_name << " / ";
      }
      cerr << endl;
    }
    pthread_mutex_unlock(&_mutex_debug);
  }
  //
#endif
#ifdef DEBUG_PRINT_TASK
  switch (task->task_id) {
  case C_FILLMATRIX:
  case C_DFULL_SYM_GAUSS:
  case C_DTRSMSCALE:
    //  case C_SPARSESCHUR:
  case C_DGEMM_LOCAL_MULT:
  case C_DGEMM_LOCAL_TWO:
  case C_DGEMM_DIRECT_TWO:
    fprintf(_fps[pid], "%s\n", task->task_name);
    *(task->fp) = &_fps[pid];
    break;
  }
#endif
  // accuessing to unsigned char does not need mutex
  task->status = TASK_WORKING;
  //
#ifdef DEBUG_THREAD_TIME
  get_realtime(&(task->t0));
  //  clock_gettime(CLOCK_MONOTONIC, &(task->t0));
#endif
  // debugging : 13 :Apr.2012 Atsushi
  task->func(task->func_arg);
#ifdef DEBUG_THREAD_DONE_PRINT
  pthread_mutex_lock(&_mutex_debug);
  {
     cerr << "pid = " << pid << " : " << task->task_name << endl;
  }
  pthread_mutex_unlock(&_mutex_debug);
#endif
  if (task->task_id == C_DFULL_SYM_GAUSS) {
    // task->quit_queue = ((C_dfull_gauss_arg*)task->func_arg)->quit;
    // task->to_next_task = ((C_dfull_gauss_arg*)task->func_arg)->to_next_task;
    for (int i = (pos + 1); i <= (pos + task->to_next_task); i++) {
      (*seq->queue)[i]->status = TASK_DONE;
      (*seq->queue)[i]->quit_queue = true;
#ifdef DEBUG_EXEC_THREAD2
    pthread_mutex_lock(&_mutex_debug);
    {
      cerr << (*seq->queue)[i]->task_name << " ";
    }
    pthread_mutex_unlock(&_mutex_debug);
#endif
    }
    //    cout << task->task_name << " / " << task->to_next_task << endl;
  }

#ifdef DEBUG_THREAD_TIME
  get_realtime(&(task->t1));
  //  clock_gettime(CLOCK_MONOTONIC, &(task->t1));
  task->thread_id = pid;
#endif
  // accuessing to unsigned char does not neet mutex
  task->status = TASK_DONE;
  //
  if (_waiting_root > 0) {
    pthread_mutex_lock(&_mutex_root);
    {
      pthread_cond_broadcast(&_cond_root);
    }
    pthread_mutex_unlock(&_mutex_root);
#ifdef DEBUG_EXEC_THREAD_IDLE
    //    pthread_mutex_lock(&_mutex_debug);
    {
      cerr << pid << " " << task->task_name
	   << "_waiting_root = " << (int)_waiting_root 
	   << " broadcast." << endl;
    }
    //    pthread_mutex_unlock(&_mutex_debug);
#endif
  }
  //
#ifdef DEBUG_EXEC_THREAD_FILE
  pthread_mutex_lock(&_mutex_file);
  {
    const int sec_n0 = convert_sec(task->t1); 
    const int sec_m0 = convert_microsec(task->t1);
    _fout << task->task_name << " / ";
    _fout << sec_n0 << " " << sec_m0 << "\t";
    for(list<C_task *>::const_iterator nt = task->parents->begin();
	nt != task->parents->end(); ++nt) {
      int sec_n1 = convert_sec((*nt)->t1);
      int sec_m1 = convert_microsec((*nt)->t1);
      if ((sec_n0 < sec_n1) || ((sec_n0 == sec_n1) && (sec_m0 < sec_m1))) {
	_fout << "*";
      }
      _fout << (*nt)->task_name << " / ";
      _fout << sec_n1 << " " << sec_m1 << "\t";
    }
    _fout << endl; // " : " << pid << endl;
  }
  pthread_mutex_unlock(&_mutex_file);
#endif
  task->status = TASK_DONE;
}

void QueueRuntime::exec_fwbw_seq()
{
  clock_t t0_cpu, t1_cpu;
  elapsed_t t0_elapsed, t1_elapsed;

  t0_cpu = clock();
  get_realtime(&t0_elapsed);
  for (vector<C_task*>::const_iterator it = _queue_fwbw->queue->begin();  
      it != _queue_fwbw->queue->end(); ++it) {
    (*it)->func((*it)->func_arg);
  }
  t1_cpu = clock();
  get_realtime(&t1_elapsed);
  if (_verbose) {
    fprintf(_fp, 
	    "execution of fw/bw : cpu time = %.4e elapsed time = %.4e\n", 
	    (double)(t1_cpu - t0_cpu) / (double)CLOCKS_PER_SEC,
	    convert_time(t1_elapsed, t0_elapsed));
  }
}

void QueueRuntime::exec_fwbw()
{
  const int num_threads = _num_threads; //_num_threads_symb;
  void* results;
  pthread_attr_t th_attr;
  pthread_t *threads;

  clock_t t0_cpu, t1_cpu;
  elapsed_t t0_elapsed, t1_elapsed;
  //  struct timespec ts0, ts1;
  int ierr;
  if (_verbose) {
    fprintf(_fp, "fwbw with %d threads\n", num_threads);
  }

 // copy dependency data : parents -> parents_work 
  {
    for (vector<C_task *>::const_iterator it = _queue_fwbw->queue->begin();
	 it != _queue_fwbw->queue->end(); ++it) {
      list<C_task *>& parents_work = *((*it)->parents_work);
      if (parents_work.size() > 0) { // to avoid double free
	parents_work.clear();  
      }
    }    // loop : it
    for (vector<C_task *>::const_iterator it = _queue_fwbw->queue->begin();
	 it != _queue_fwbw->queue->end(); ++it) {
      list<C_task *>& parents = *((*it)->parents);
      list<C_task *>& parents_work = *((*it)->parents_work);
      (*it)->status = TASK_WAITING;  // reset status
      (*it)->quit_queue = false;
      if (parents_work.size() == 0) {
#ifdef SX_ACE
	   for (list<C_task *>::const_iterator lt = parents.begin();
		lt != parents.end(); ++lt) {
	     parents_work.push_back(*lt);
	   }
#else // SX-ACE C++ rev 110 with C++98/03 does not understand back_inserter
	std::copy(parents.begin(), parents.end(), back_inserter(parents_work));
#endif
      }
      (*it)->broadcast_deadlock = 0;
    }    // loop : it
    _queue_fwbw->pos = _queue_fwbw->begin;
  }
#if 0  // for debugging : execution in sereial
  pos = _queue_fwbw->begin;
  cerr << "pos = " << pos << endl;
  for (int i = 0; i < pos; i++) {
    C_task * task = (*_queue_fwbw->queue)[i];
    task->func(task->func_arg);
    task->status = TASK_DONE;
  }
#endif
  threads = new pthread_t[num_threads]; 

#ifdef DEBUG_EXEC_THREAD
  pthread_mutex_init(&_mutex_debug, NULL);
#endif
  pthread_mutex_init(&_mutex_dependency, NULL);
  pthread_cond_init(&_cond_root, NULL);

  ierr = pthread_mutex_init(&_mutex_root, NULL);
  if (ierr != 0) {
    fprintf(stderr, " pthread_mutex_init(&_mutex_root, NULL) %s %d : %d\n",
	    __FILE__, __LINE__, ierr);
  }
#ifdef DEBUG_EXEC_THREAD
  ierr = pthread_mutex_init(&_mutex_debug, NULL);
  if (verbose && (ierr != 0)) {
    fprintf(_fp, "%s %d : pthread_mutex_init(&_mutex_debug, NULL) %d\n", 
	    __FILE__, __LINE__, ierr);
  }
#endif
  pthread_attr_init(&th_attr);
  pthread_attr_setdetachstate(&th_attr, PTHREAD_CREATE_JOINABLE);       

#ifdef DEBUG_EXEC_THREAD_FILE
  {
    int pid = get_process_id();
    char fname[256];
    sprintf(fname, "task-s.%d.data", pid);
    _fout.open(fname);
  }
#endif

  t0_cpu = clock();
  get_realtime(&t0_elapsed);
  //  clock_gettime(CLOCK_REALTIME, &ts0);

  THREAD_QUEUE_EXEC **params = new THREAD_QUEUE_EXEC*[num_threads];
  for (int p = 0; p < num_threads; p++) {
    params[p] = new THREAD_QUEUE_EXEC(p, num_threads, this);
    int pid = pthread_create(&threads[p], &th_attr, 
			     &thread_queue_fwbw_,
			     (void *)params[p]);
    if (pid != 0) {
      if (_verbose) {
	fprintf(_fp, "bad thread creation ? : %d\n", pid);
      }
      exit(0);
    }
  }
  pthread_attr_destroy(&th_attr);
  for (int p = 0; p < num_threads; p++) {
    int pid = pthread_join(threads[p], &results);
    if (pid != 0) {
      if (_verbose) {
	fprintf(_fp, "bad thread creation ? : %d\n", pid);
      }
      exit(0);
    }
    delete params[p];
  }
  delete [] params;
#if 0 // for debugging : execution in sereial
  pos = _queue_fwbw->pos;
  cerr << "pos = " << pos << endl;
  for (int i = pos; i < _queue_fwbw->queue->size(); i++) {
    C_task * task = (*_queue_fwbw->queue)[i];
    task->func(task->func_arg);
    task->status = TASK_DONE;
  }
#endif
  t1_cpu = clock();
  get_realtime(&t1_elapsed);
  // clock_gettime(CLOCK_REALTIME, &ts1);
  if (_verbose) {
    fprintf(_fp, 
	    "execution of fw/bw : cpu time = %.4e elapsed time = %.4e\n", 
	    (double)(t1_cpu - t0_cpu) / (double)CLOCKS_PER_SEC,
	    convert_time(t1_elapsed, t0_elapsed));
  }

#ifdef DEBUG_EXEC_THREAD_FILE
  _fout.close();
#endif

#ifdef DEBUG_THREAD_TIME
  {
    int pid = get_process_id();
    char filename[256];		  
    FILE *fp;
    sprintf(filename, "threadtime-fwbw.%d.data", pid);
    if ((fp = fopen(filename, "w")) == NULL) {
      fprintf(stderr,
	      "%s %d : fail to open %s\n",
	      __FILE__, __LINE__, filename);
      exit(-1);
    }
    else {
      fprintf(fp, "queue_fwbw\n");
      for (int j = _queue_fwbw->begin; j < _queue_fwbw->end; j++) {
	double t0 = (convert_sec((*_queue_fwbw->queue)[j]->t0) + 
		     convert_microsec((*_queue_fwbw->queue)[j]->t0) * 1.0e-6);
	
	fprintf(fp, "%2d:\t%-16s\t%2d\t%4d\t%7d\t%7d\t%7d\t%7d",
		(*_queue_fwbw->queue)[j]->thread_id,
		(*_queue_fwbw->queue)[j]->task_name,
		(int)(*_queue_fwbw->queue)[j]->parents_work->size(),
		(*_queue_fwbw->queue)[j]->broadcast_deadlock,
		convert_sec((*_queue_fwbw->queue)[j]->t0),
		convert_microsec((*_queue_fwbw->queue)[j]->t0),
		convert_sec((*_queue_fwbw->queue)[j]->t1),
		convert_microsec((*_queue_fwbw->queue)[j]->t1));
	for (list<C_task *>::const_iterator it = (*_queue_fwbw->queue)[j]->parents->begin();
	     it != (*_queue_fwbw->queue)[j]->parents->end(); ++it) {
	  double tt1 = (convert_sec((*it)->t1) + 
			convert_microsec((*it)->t1) * 1.0e-6);
	  if (t0 < tt1) {
	    fprintf(fp, " [ %10s ]", (*it)->task_name);
	  }
	}
	fprintf(fp, "\n");
      }
      fclose(fp);
    } // if ((fp = fopen(filename, "a")) == NULL)
  }
#endif

  pthread_mutex_destroy(&_mutex_root);
#ifdef DEBUG_EXEC_THREAD
  pthread_mutex_destroy(&_mutex_debug);
#endif

  delete [] threads;
}


void *thread_queue_fwbw_(void *arg)
{
  THREAD_QUEUE_EXEC *params = (THREAD_QUEUE_EXEC *)arg;
  params->dissectionRuntime->thread_queue_fwbw(params->id, 
					       params->num_threads);
  pthread_exit(arg);

  return (void *)NULL;
}

// #define DEBUG_EXEC_THREAD_IDLEa

void QueueRuntime::thread_queue_fwbw(const int pid, const int num_threads)
{
  //  const int pid = params->id;
  //  const int num_threads = params->num_threads;
  // greedy -- better to be in seperated function
  int pos, end;
  C_task_seq* it = _queue_fwbw;
  //  cerr << "pid = " << pid 
  //     << " _queue_symb : pos = " << it->pos << " end = " << it->end << endl;
  while(1) {
    pthread_mutex_lock(&_mutex_root);
    {
      pos = it->pos;
      end = it->end;
      if (pos < end) {
	it->pos = it->pos + (*it->queue)[pos]->atomic_size;
      }
    }
    pthread_mutex_unlock(&_mutex_root);
    if (pos >= end) {
      break;
    }
    C_task *task_p = (*it->queue)[pos];
    for (int i = 0; i < task_p->atomic_size; i++) {
      C_task *task = (*it->queue)[pos + i];
      int waiting = check_parents_done(task, &_mutex_dependency);
      while (waiting > 0) {
	pthread_mutex_lock(&_mutex_root);
	{
	  _waiting_root++;
	  if (_waiting_root == num_threads) {
	    // wake up other threads
#ifdef DEBUG_EXEC_THREAD_IDLEa
	    pthread_mutex_lock(&_mutex_debug);
	    {
	      cerr << pid << " " << task->task_name 
		   << " _waiting_root = " << (int)_waiting_root 
		   << " dead locked! : emergency broadcast." << endl;
	    }
	    pthread_mutex_unlock(&_mutex_debug);
#endif // DEBUG_EXEC_THREAD_IDLEa
	    pthread_cond_broadcast(&_cond_root); 
	    task->broadcast_deadlock++;
	  }
#ifdef DEBUG_EXEC_THREAD_IDLEa
	  pthread_mutex_lock(&_mutex_debug);
	  {
	    cerr << pid << " " << task->task_name << " "
		 << "waiting = "<< waiting << " : " 
		 << (int)_waiting_root << " "
		 << __FILE__ << " " << __LINE__ << " ";
	    for (list<C_task *>::const_iterator jt = task->parents_work->begin();
		 jt != task->parents_work->end(); ++jt) {
	      cerr << (*jt)->task_name << " ";
	    }
	    cerr << " : sleeping" << endl;
	  }
	  pthread_mutex_unlock(&_mutex_debug);
#endif // DEBUG_EXEC_THREAD_IDLEa
	  pthread_cond_wait(&_cond_root, &_mutex_root);
	  _waiting_root--;
#ifdef DEBUG_EXEC_THREAD_IDLEa
	  pthread_mutex_lock(&_mutex_debug);
	  {
	    cerr << pid << " " << task->task_name 
		 << " waked up : _waiting_root = " << (int)_waiting_root 
		 << endl;
	  }
	  pthread_mutex_unlock(&_mutex_debug);
#endif // DEBUG_EXEC_THREAD_IDLEa
	  // other thread broacast wake up and decreasing _waiting_root
	} 
	waiting = check_parents_done(task, &_mutex_dependency);
	pthread_mutex_unlock(&_mutex_root);
      } // while (waiting > 0)
      pthread_mutex_lock(&_mutex_dependency);
      task->status = TASK_WORKING;
      pthread_mutex_unlock(&_mutex_dependency);
#ifdef DEBUG_THREAD_TIME
      get_realtime(&(task->t0));
#endif
      task->func(task->func_arg);
#ifdef DEBUG_THREAD_TIME
      get_realtime(&(task->t1));
      task->thread_id = pid;
#endif
      pthread_mutex_lock(&_mutex_dependency);  // for completely greedy algorithm
      task->status = TASK_DONE;
      pthread_mutex_unlock(&_mutex_dependency);
      if (_waiting_root > 0) {
	pthread_mutex_lock(&_mutex_root);
	{
#ifdef DEBUG_EXEC_THREAD_IDLEa
	  pthread_mutex_lock(&_mutex_debug);
	  {
	    cerr << pid << " " << task->task_name 
		 << " _waiting_root = " << (int)_waiting_root 
		 << " broadcast." << endl;
	  }
	  pthread_mutex_unlock(&_mutex_debug);
#endif // DEBUG_EXEC_THREAD_IDLEa
	  pthread_cond_broadcast(&_cond_root);
	}
	pthread_mutex_unlock(&_mutex_root);
      }
#ifdef DEBUG_EXEC_THREAD_IDLE
#if 0 // better with sprintf
      string stmp = (to_string(pid) + " " + to_string(pos) + " "
	+ task->task_name + "\t" + convert_sec(task->t0) << "." << convert_microsec(task->t0) << " : "
	   << convert_sec(task->t1) << "." << convert_microsec(task->t1) << endl;
      pthread_mutex_lock(&_mutex_debug);
      {
	cerr << stmp.str();
      }
      pthread_mutex_unlock(&_mutex_debug);
#endif
#endif // DEBUG_EXEC_THREAD_IDLEa
    } // loop : i
  } // while (1)
#ifdef DEBUG_EXEC_THREAD
  pthread_mutex_lock(&_mutex_debug);
  {
    cerr << pid << " " << it->task_name 
	 << " @ " << it->task_id << " greedy end." << endl;
  }
  pthread_mutex_unlock(&_mutex_debug);
#endif     
}
