/*! \file QueueRuntime.hpp
    \brief management of threads for Dissection Matrix
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

#ifndef _QUEUE_RUNTIME_
# define _QUEUE_RUNTIME_
#ifdef _MSC_VER
#include <iterator>
#endif
#ifdef SX_ACE
#include <iterator>
#endif
#include "Driver/C_threads_tasks.hpp"
#include "Compiler/OptionLibrary.h"
#include "Compiler/elapsed_time.hpp"
#include <time.h>
#include <cstdlib>
#include <string>
#include <pthread.h>

class QueueRuntime
{
public:
  QueueRuntime(int nb_doms, int num_threads, const bool verbose, FILE *fp);

  ~QueueRuntime(); 

  void generate_queue(C_task_seq* &_queue_symb,
		      list<C_task_seq*>* &_queue_static,
		      vector<C_task_seq *>* &_queue_dynamic,
		      list<C_task *> &queue_null,
     //		      list<C_task *> &queue_dummy,
		      Dissection::Tree* btree,
		      vector<int>* children,
		      vector<C_task *>* tasks_SparseSymb,
		      vector<C_task *>* tasks_SparseNum,
		      vector<C_task *>* tasks_SparseLocalSchur,
		      vector<C_task *>* tasks_DFillSym,
		      vector<C_task *>* tasks_DFullLDLt,
		      vector<C_task *>* tasks_DTRSMScale,
		      vector<C_task *>* tasks_DSymmGEMM,
		      vector<C_task *>** tasks_Dsub,
		      vector<C_task *>* tasks_deallocLower,
		      vector<C_task *>* tasks_deallocLocalSchur,
		      long **nops_queue,
		      vector<vector<int> > all_fathersIndex,
		      vector<int> nrow_DFullLDLt,
		      vector<bool> isMergedDTRSM,
		      vector<bool> isDividedDTRSM,
		      const int level_last);

  void write_dependency(FILE *fp);
  void exec_symb_fact();
  void exec_num_fact(const int called);
  void exec_num_fact_debug();
  void execute_task_debug(C_task_seq *seq, int pos, 
			  int *permute_block, 
			  int pid);
  void exec_fwbw();
  void exec_fwbw_seq();
  void thread_queue_symb_factorize(const int pid, const int num_threads);
  void thread_queue_num_factorize(const int pid, const int num_threads);
  void thread_queue_fwbw(const int pid, const int num_threads);

  void thread_queue_C_DFULL(const int pid, const int num_thraeds,
			    C_task_seq *it,
			    int *permute_block, 
			    const int cnt_cdfull,
			    const int zone_first_entered,
			    const int zone_idxn);

  void thread_queue_parallel_dynamic(const int pid, const int num_threads,
				     C_task_seq *it,
				     int *permute_block, 
				     const int zone_idxn);

  void thread_queue_parallel_static(const int pid, const int num_threads,
				    C_task_seq *it, 
				    int *permute_block, 
				    const int zone_idxp,
				    const int zone_idxn,
				    const int zone_first_entered,
				    const int zone_last_entered);

  void thread_queue_single(const int pid, const int num_threads,
			   C_task_seq *it,
			   int *permute_block, 
			   const int zone_idxn);

  void taskstatus_mutex_lock(pthread_mutex_t *mutex_dependency);
  void taskstatus_mutex_unlock(pthread_mutex_t *mutex_dependency);

  int execute_task_dynamic_buffer(int *permute_block, 
				  int pid,
				  pthread_mutex_t *mutex_dependency);
  
  void execute_task(C_task_seq *seq, int pos, 
		   int *permute_block, 
		   int pid,
		   pthread_mutex_t *mutex_dependency);
  void execute_task(C_task_seq *seq, int pos, 
		    int *permute_block,
		    int pid);
  
  int check_parents_done(C_task *it, pthread_mutex_t *mutex_dependency);
  int check_parents_done(C_task *it);

  void set_queue_fwbw(C_task_seq* &queue_fwbw)
  {
    _queue_fwbw = queue_fwbw;
  }

  QueueRuntime(const QueueRuntime &s) { } // copy constrcutor : dummy

private:
  int _nb_doms;
  int _num_threads;

  C_task_seq* _queue_symb;     
  C_task_seq* _queue_fwbw;
  list<C_task_seq*> *_queue_static;
  vector<C_task_seq *> *_queue_dynamic;
  //  list<C_task *>* _queue_dummy;
  pthread_mutex_t _mutex_root; 
  pthread_mutex_t _mutex_dependency; 
  pthread_mutex_t *_mutex_group;
  pthread_mutex_t _mutex_debug;
  pthread_mutex_t _mutex_file;
  
  pthread_cond_t _cond_root;
  int *_zone_entered;
  int *_zone_finished;

  int ***_group_entered;
  int ***_group_finished;
  int ***_group_task_ends;
  int ***_group_static_assigned;

  int *_group_task_id;
  int *_zone_static_assigned;

  int **_begins;
  int **_ends;
  int ***_begins_group;
  int ***_ends_group;
  long **_group_nops;

  unsigned char _waiting_root; // num. of threads should be less than 255
  int _phase_dynamic;
  int _queue_dynamic_pos_start;
  int _queue_dynamic_pos;
  int _queue_dynamic_notcopied;
  //
  //  ofstream _fout;

  // forward-backward substitutions

  bool _verbose;
  FILE *_fp;
  bool _isSym;

}; // End class DissictionQueue

struct THREAD_QUEUE_EXEC {
  int id;
  int num_threads;
  QueueRuntime* dissectionRuntime;
  THREAD_QUEUE_EXEC(int id_,
		    int num_threads_,
		    QueueRuntime* dissectionRuntime_) :
  id(id_), 
  num_threads(num_threads_),
  dissectionRuntime(dissectionRuntime_) {}
};

void copytask_list2seq(list<C_task_seq*> &queue_static,
		       list<C_task_seq*> queue_lists,
		       list<C_task *> &queue_null,
		       string task_name,
		       int task_id, 
		       int mutex_id,
		       int parallel_single,
		       int num_threads,
		       int level,
		       int phase);

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
		       vector<int>* children);

void task_assign_diag2(list<C_task_seq*> *&queue_static,
		       vector<C_task *>* &tasks_queue, 
		       vector<C_task *> &queue_null,
    //		       vector<C_task *> &queue_dummy,
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
		       vector<int>* children);

void allocate_int2d(int **&array, const int num_threads);
void allocate_unsigned2d(long **&array, const int num_threads);
void allocate_int3d(int ***&array, const int num_threads);
void deallocate_int2d(int **&array);
void deallocate_unsigned2d(long **&array);
void deallocate_int3d(int ***&array);

void *thread_queue_num_factorize_(void *arg);
void *thread_queue_symb_factorize_(void *arg);
void *thread_queue_fwbw_(void *arg);

#endif
