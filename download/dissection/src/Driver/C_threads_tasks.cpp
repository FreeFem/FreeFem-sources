/*! \file C_threads_tasks.hpp
    \brief tasks executed asynchronously with threads
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

#include "Compiler/OptionLibrary.hpp"
#include <iostream>
#include <float.h>
#include <cmath> // for isnan()

#include "Compiler/blas.hpp"
#include "Driver/C_threads_tasks.hpp"
#include "Driver/C_KernDetect.hpp"
#include "Algebra/VectorArray.hpp"
#include "Compiler/DissectionIO.hpp"

// #define DEBUG_SVD
#ifdef DEBUG_SVD
#include "mkl_lapack.h"
#endif

void C_dummy(void *arg_) { };

//#define DEBUG_ERASENULLPARENTS
int EraseNullParents(vector<C_task *> &queue)
{
  int count = 0;
  for (vector<C_task *>::const_iterator it = queue.begin(); it != queue.end();
       ++it) {
    for (list<C_task *>::iterator jt = (*it)->parents->begin(); 
	 jt != (*it)->parents->end(); ) {
      if (*(*jt)->ops_complexity == 0L) {
	count++;
	jt = (*it)->parents->erase(jt);
      }
      else {
	++jt;
      }
    }
  }
  return count;
}

int EraseNullParents(C_task * task)
{
  int count = 0;
  for (list<C_task *>::iterator jt = task->parents->begin(); 
	 jt != task->parents->end(); ) {
    if (*(*jt)->ops_complexity == 0L) {
      count++;
      jt = task->parents->erase(jt);
    }
    else {
      ++jt;
    }
  }
  return count;
}

template<typename T, typename U>
void C_SparseSymbFact(void *arg_)
{
  C_SparseSymbFact_arg<T, U> *arg = (C_SparseSymbFact_arg<T, U> *)arg_;
  TridiagBlockMatrix<T, U>** tridiag = arg->tridiag;
  const int colors = arg->colors;
  int *color_mask = arg->color_mask;
  const int nrow =  arg->nrow;
  const int nnz0 =  arg->csr_diag->nnz;
  int* prow0 =      arg->csr_diag->ptRows;
  int* indcols0 =   arg->csr_diag->indCols;
  int* indvals0 =   arg->csr_diag->indVals;
  //  bool verbose =    arg->verbose;
  //  FILE *fp =        *(arg->fp);
  for (int i = 0; i < colors; i++) {
    tridiag[i]->SymbolicFact((i + 1), colors, color_mask, nrow,
			     nnz0, prow0, indcols0, indvals0); //
  }  
}

template
void C_SparseSymbFact<double, double>(void *arg_);

template
void C_SparseSymbFact<quadruple, quadruple>(void *arg_);

template
void C_SparseSymbFact<complex<double>, double>(void *arg_);

template
void C_SparseSymbFact<complex<quadruple>, quadruple>(void *arg_);

template
void C_SparseSymbFact<float, float>(void *arg_);

template
void C_SparseSymbFact<complex<float>, float>(void *arg_);

//

void c_getrealtime_(uint64_t &tmprofiles, const int &m)
{
  elapsed_t *tmprfs = (elapsed_t *)tmprofiles;
  get_realtime(&tmprfs[m]);
}

void c_fileout_(uint64_t &fp_prt, char *strgs, const int &force_stderr)
{
  if (force_stderr == 1) {
    fprintf(stderr, "%s\n", strgs);
  }
  else {
    FILE *fp = (FILE *)fp_prt;
    fprintf(fp, "%s\n", strgs);
  }
}

int compare_source_dist_index(const void *_a, const void *_b)
{
  source_dist_index *a = (source_dist_index *)_a;
  source_dist_index *b = (source_dist_index *)_b;
  if (a->global_i == b->global_i) {
    return (int)(a->global_j - b->global_j);
  }
  else {
    return (int)(a->global_i - b->global_i);
  }
}

int convert_array2strip(list<index_strip> &strips, 
			const vector<int>& array)
{
  const int size = array.size();
  strips.clear();
  int in_strip = 0;
  int begin_dst, width;

  for (int i = 0; i < size; i++) {
    if ((in_strip == 0) && (array[i] != (-1))) {
      in_strip = 1;
      begin_dst = i;
    }
    if (i < (size - 1)) {
      if ((in_strip == 1) && 
	  ((array[i + 1] == (-1)) || (array[i + 1] != (array[i] + 1)))) {
	in_strip = 0;
	width = (i + 1) - begin_dst;
	// constructor
	strips.push_back(index_strip(begin_dst, 
				     array[begin_dst], // begin_src
				     width));
	if (array[i + 1] != (array[i] + 1)) {
	  begin_dst = i + 1;
	}
      }
    }
  }
  if (in_strip == 1) {
    width = size - begin_dst;
    strips.push_back(index_strip(begin_dst, 
				 array[begin_dst], // begin_src
				 width));
  }

  return strips.size();
}

int combine_two_strips(list<index_strip> &stripsa,
		       list<index_strip> &stripsb,
		       list<index_strip2> &stripsc,
		       list<index_strip> &strips0, 
		       list<index_strip> &strips1,
		       const int size) 
{
  stripsa.clear();
  stripsb.clear();
  stripsc.clear();

  if (strips0.size() == 0 && strips1.size() == 0) {
    return 0;
  }

  if (strips0.size() == 0) {
    stripsa.clear();
    stripsb.clear();
    for (list<index_strip>::const_iterator it = strips1.begin();
	 it != strips1.end();
	 ++it) {
      stripsc.push_back(index_strip2((*it).begin_dst,
				     -1,
				     (*it).begin_src, 
				     (*it).width));
    }
    return strips1.size();
  }

  if (strips1.size() == 0) {
    stripsa.clear();
    stripsb.clear();
    for (list<index_strip>::const_iterator it = strips0.begin();
	 it != strips0.end();
	 ++it) {
      	  stripsc.push_back(index_strip2((*it).begin_dst, 
					 (*it).begin_src, 
					 -1,
					 (*it).width));
    }
    return strips0.size();
  }

  vector<int> ary0(size, (-1));
  vector<int> ary1(size, (-1));
  int in_strip0, in_strip;
  int begin_dst = (-1), width;

  for (list<index_strip>::const_iterator it = strips0.begin();
       it != strips0.end();
       ++it) {
    int i0 = (*it).begin_dst;
    int i1 = (*it).begin_src;
    for (int i = 0; i < (*it).width; i++) {
      ary0[i0++] = i1++;
    }
  }
  for (list<index_strip>::const_iterator it = strips1.begin();
       it != strips1.end();
       ++it) {
    int i0 = (*it).begin_dst;
    int i1 = (*it).begin_src;
    for (int i = 0; i < (*it).width; i++) {
      ary1[i0++] = i1++;
    }
  }

  in_strip0 = in_strip = 0;

  for (int i = 0; i < size; i++) {
    // save current status on strip
    in_strip0 = in_strip;
    if ((in_strip == 0 || in_strip == 2) && (ary0[i] != (-1))) {
      in_strip += 1;
    }
    if ((in_strip == 0 || in_strip == 1) && (ary1[i] != (-1))) {
      in_strip += 2;
    }
    if (in_strip > in_strip0) {
      // strip0 or strip1 -> union of strips0 and strips1 
      if (in_strip0 > 0) {
	width = i - begin_dst;
	switch (in_strip0) {
	case 1:
	  stripsa.push_back(index_strip(begin_dst, ary0[begin_dst], width));
	  break;
	case 2:
	  stripsb.push_back(index_strip(begin_dst, ary1[begin_dst], width));
	  break;
	}
      }
      begin_dst = i;
    }
    // save current status on strip
    in_strip0 = in_strip;
    //
    if (i < (size - 1)) {
      // decreasing into strip0 or strip1 or null
      if ((in_strip == 1 || in_strip == 3) && 
	  ((ary0[i + 1] == (-1)) || (ary0[i + 1] != (ary0[i] + 1)))) {
	in_strip -= 1;
      }
      if ((in_strip == 2 || in_strip == 3) && 
	  ((ary1[i + 1] == (-1)) || (ary1[i + 1] != (ary1[i] + 1)))) {
	in_strip -= 2;
      }
      if (in_strip0 > in_strip) {
	// output union of strip0 and strip1, strip0 or strip1
	width = (i + 1) - begin_dst;
	switch (in_strip0) {
	case 1:
	  stripsa.push_back(index_strip(begin_dst, ary0[begin_dst], width));
	  break;
	case 2:
	  stripsb.push_back(index_strip(begin_dst, ary1[begin_dst], width));
	  break;
	case 3:
	  stripsc.push_back(index_strip2(begin_dst, 
					 ary0[begin_dst],
					 ary1[begin_dst], width));
	  break;
	}
	// strip0 or strip1
	if (in_strip > 0) {
	  begin_dst = i + 1;
	}
	in_strip0 = in_strip;
      } // if (in_strip0 > in_strip) {
    }   // if (i < (size - 1))
  }     // loop : i
  if (in_strip > 0) {
    width = size - begin_dst;

    switch (in_strip0) {
    case 1:
      stripsa.push_back(index_strip(begin_dst, ary0[begin_dst], width));
      break;
    case 2:
      stripsb.push_back(index_strip(begin_dst, ary1[begin_dst], width));
      break;
    case 3:
      stripsc.push_back(index_strip2(begin_dst, 
				     ary0[begin_dst], 
				     ary1[begin_dst], width));
      break;
    }
  }
  return (stripsa.size() + stripsb.size() + stripsc.size());
}
void copy_two_strips(list<index_strip2> &strips2,
		     list<index_strip> &strips0, 
		     list<index_strip> &strips1) 
{
  strips2.clear();

  list <index_strip>::const_iterator mt0 = strips0.begin();
  list <index_strip>::const_iterator mt1 = strips1.begin();
  for ( ; ((mt0 != strips0.end()) && (mt1 != strips1.end())); ++mt0, ++mt1) {
    strips2.push_back(index_strip2((*mt0).begin_dst,
				   (*mt0).begin_src,
				   (*mt1).begin_src,
				   (*mt0).width));
  }
}

void split_two_strips(list<index_strip> &strips0,
		      list<index_strip> &strips1,
		      list<index_strip2> &strips2,
		      index_strip strip0,
		      index_strip strip1) 
{
  // assumption : union of strip0 and strip1 is not empty
  // strips0 = strip0 \setminus strip1
  // strips1 = strip1 \setminus strip0
  // strips01 = strip0 \cap strip1
  strips0.clear();
  strips1.clear();
  strips2.clear();
  int begin_dst = 0;
  int width = 0;
  int begin_src0 = strip0.begin_src;
  int begin_src1 = strip1.begin_src;
  int end_dst0 = strip0.begin_dst + strip0.width;
  int end_dst1 = strip1.begin_dst + strip1.width;

  if (strip0.begin_dst < strip1.begin_dst) {
    width = strip1.begin_dst - strip0.begin_dst;
    strips0.push_back(index_strip(strip0.begin_dst,
				  strip0.begin_src,
				  width));
    begin_src0 += width;
  }
  else if (strip0.begin_dst > strip1.begin_dst) {
    width = strip0.begin_dst - strip1.begin_dst;
    strips1.push_back(index_strip(strip1.begin_dst,
				  strip1.begin_src,
				  width));
    begin_src1 += width;
  }
  begin_dst += width;
  const int end_dst = end_dst0 < end_dst1 ? end_dst0 : end_dst1;
  const int width2 = end_dst - begin_dst;
  strips2.push_back(index_strip2(begin_dst,
				 begin_src0,
				 begin_src1,
				 width2));
  if (end_dst0 < end_dst1) {
    strips1.push_back(index_strip(end_dst,              // 27 Nov.2016 bug found
				  begin_src1 + width2,  // begin_src0 + width
				  end_dst1 - end_dst0));
  }
  else if (end_dst0 > end_dst1) {
    strips0.push_back(index_strip(end_dst,
				  begin_src0 + width2,  // begin_src1 + width
				  end_dst0 - end_dst1));
  }
}

void print_strips(const char *name, list<index_strip> &strips,
		  FILE *fp)
{
  fprintf(fp, "%s : ", name);
  for (list <index_strip>::const_iterator mt = strips.begin();
       mt != strips.end(); ++mt) {
    fprintf(fp, "[ %d %d %d ] ",  
	    (*mt).begin_dst, (*mt).begin_src, (*mt).width);
  }
  fprintf(fp, "\n");
}

void print_strips(const char *name, list<index_strip2> &strips, FILE *fp)
{
  fprintf(fp, "%s : ", name);
  for (list <index_strip2>::const_iterator mt = strips.begin();
       mt != strips.end(); ++mt) {
    fprintf(fp, "[ %d %d %d %d ] ",  
	    (*mt).begin_dst, (*mt).begin_src0, (*mt).begin_src1, (*mt).width);
  }
  fprintf(fp, "\n");
}
  
		     
bool C_task_seq_complexity_smaller(C_task_seq *a, C_task_seq *b)
{
  if (a->ops_complexity < b->ops_complexity) {
    return true;
  }
  else {
    return false;
  }
}
bool C_task_seq_complexity_greater(C_task_seq *a, C_task_seq *b)
{
  if (a->ops_complexity > b->ops_complexity) {
    return true;
  }
  else {
    return false;
  }
}

bool C_task_seq_beginidx_smaller(C_task_seq *a, C_task_seq *b)
{
  if (a->begin < b->begin) {
    return true;
  }
  else {
    return false;
  }
}
   

void assign_tasks_statically(list<C_task_seq*> *queue_static,
			     list<C_task_seq*> &queue_dynamic,
			     long *nops_sum,
			     list<C_task_seq *> &task_seq_tmp,
			     const int task_id,
			     const char *task_name_,
			     const int level,
			     const int phase,
			     const long nops_block_total,
			     int num_threads)
{
  const long nops_per_thread = nops_block_total / (long)num_threads;
#ifdef DEBUG_PREPARE_THREAD
  cout << "nops_bock_total = " << nops_block_total 
       << " nops_per_thread = " << nops_per_thread << endl;
#endif
  list<C_task_seq *>::const_iterator it;
  int task_begin;
  vector<list<C_task_seq *> > queue_divided(num_threads);

  vector<int> excluding_list(num_threads, 0);
  long *nops_sum0 = new long[num_threads]; // for debugging
#ifdef DEBUG_PREPARE_THREAD
  cout << "excluded indices ";
#endif
  int excluded = 0;
  for (int i = 0; i < num_threads; i++) {
    nops_sum0[i] = nops_sum[i];                    // for debugging
    if (nops_sum[i] >= nops_per_thread) {
      excluding_list[i] = 1;
      excluded++;
#ifdef DEBUG_PREPARE_THREAD
      cout << i << " ";
#endif
    }
  }
#ifdef DEBUG_PREPARE_THREAD
  cout << endl;
#endif
  if (excluded == num_threads) {
    for (int i = 0; i < num_threads; i++) {
      excluding_list[i] = 0;
    }
  }
  int ll = (-1);
  while(excluding_list[++ll] == 1) {
#ifdef DEBUG_PREPARE_THREAD
    if (nops_sum[ll] > 0) {
      cout << "+ll = " << ll << " start =  " << nops_sum[ll] << " " ;
    }
#endif
  }

  it = task_seq_tmp.begin();
  task_begin = (*it)->begin;
  int flag = 0;
  for ( ; it != task_seq_tmp.end(); ++it) {
    long nops = 0L;
    for (int i = (*it)->begin; i < (*it)->end; i++) {
      flag = 0;
      vector<C_task* >& jt = *((*it)->queue);
      const long ltmp = *(jt[i]->ops_complexity);
      int ii = 0; // initializing
      if ((nops_sum[ll] < nops_per_thread) &&
	  ((nops_sum[ll] + (ltmp / 2UL)) >= nops_per_thread)) {
	if (ll == (num_threads - 1)) { // include blocks into the last thread
	  nops_sum[ll] += ltmp;
	  nops += ltmp;
#if 1
	  flag = (-2);
	  ii = i;
#else
	  flag = 0;
#endif
	}
	else {
	  ii = i;
	  flag = 1;
	}
      }
      else {
	nops_sum[ll] += ltmp;
	nops += ltmp;
	if (nops_sum[ll] > nops_per_thread) {
	  if (ll == (num_threads - 1)) { // include blocks into the last thread
#if 1
	    flag = (-2);
	    ii = i + 1;
#else
	    flag = 0;
#endif
	  }
	  else {
	    ii = i + 1;
	    flag = (-1);
	  }
	}
      }
      // 7 May 2014 : to accpet atomic_size > 1
      if (ii < (*it)->queue->size()) {
	ii -= (*(*it)->queue)[ii]->atomic_id;
      }
      if (flag != 0) {
	if (nops > 0L) { // to avoid zero addition when nops_sum[ll] exceed
	                  // the nops_per_thread at the beginning
	  string task_name = (task_name_ + to_string(level) + " : " 
			      + (*it)->task_name + " :: "
			      + to_string(task_begin) + " : "
			      + to_string(ii));
	  //	  char *task_name_cstr = new char[task_name.str().size() + 1];
	  //	  strcpy(task_name_cstr, task_name.str().c_str());
	  //	 cout << "-ll = " << ll << " nops_sum = " << nops_sum[ll] << " "
	  //     << task_name_cstr << endl;
	  C_task_seq* tmp = 
	    new C_task_seq(task_id,
			   task_name,
		 	   (-1),            // mutex_id
			   TASK_SINGLE,
			   1,
			   level,
			   phase,
			   (*it)->queue,
			   task_begin,
			   ii,
			   nops);
	  queue_divided[ll].push_back(tmp);
	}
#if 1
	if (flag == (-2)) {
	  task_begin = ii;
	  break;  // loop : i
	}
#endif
	while(excluding_list[++ll] == 1) {
	  //  cout << "+ll = " << ll << " start =  " << nops_sum[ll] << " " ;
	}
	if (ii < (*it)->end) { // flag == 1
	  task_begin = ii;
	}
	else {
	  list<C_task_seq *>::const_iterator kt = it;
	  ++kt;
	  if (kt != task_seq_tmp.end()) { // for the next queue
	    task_begin = (*kt)->begin;
	  }
	}
	if (flag == 1) {
	  nops = ltmp;
	  nops_sum[ll] += ltmp;
	}
	else {
	  nops = 0L;
	}
      }
    } // loop : i
#if 1              
    if (flag == (-2)) {  // special case moving static to dynamic
      break; // loop : it
    }
#endif
    if (flag >= 0) {
      string task_name = (task_name_ + to_string(level) + " : " 
			  + (*it)->task_name + " :: "
			  + to_string(task_begin) + " : "
			  + to_string((*it)->end));
      C_task_seq* tmp = 
	new C_task_seq(task_id,
		       task_name,
		       //		       task_name_cstr,
		       (-1),            // mutex_id
		       TASK_SINGLE,
		       1,
		       level,
		       phase,
		       (*it)->queue,
		       task_begin,
		       (*it)->end,
		       nops);
      queue_divided[ll].push_back(tmp);
      list<C_task_seq *>::const_iterator jt = it;
      ++jt;
      if (jt != task_seq_tmp.end()) { // for the next queue
	task_begin = (*jt)->begin;
      }
    }
  }     // loop : it
#if 1  
  if (flag == (-2)) {
    //    task_begin -= (*(*it)->queue)[task_begin]->atomic_id;
    long nops = 0L;
    for (int i = task_begin; i < (*it)->end; i++) {
      nops += *((*(*it)->queue)[i]->ops_complexity);
    }
    string task_name = (task_name_ + to_string(level) + " : " 
			+ (*it)->task_name + " :: "
			+ to_string(task_begin) + " : "
			+ to_string((*it)->end));
    C_task_seq* tmp = 
      new C_task_seq(task_id,
		     task_name,
		     0,            // mutex_id
		     TASK_PARALLEL,
		     num_threads,
		     level,
		     phase,
		     (*it)->queue,
		     task_begin,
		     (*it)->end,
		     nops);
    queue_dynamic.push_back(tmp);
    ++it;
    for ( ; it != task_seq_tmp.end(); ++it) {
      for (int i = (*it)->begin; i < (*it)->end; i++) {
	nops += *((*(*it)->queue)[i]->ops_complexity);
      }
      string task_name = (task_name_ + to_string(level) + " : " 
			  + (*it)->task_name + " :: "
			  + to_string((*it)->begin) + " : "
			  + to_string((*it)->end));
      C_task_seq* tmp = 
	new C_task_seq(task_id,
		       task_name,
		       0,            // mutex_id
		       TASK_PARALLEL,
		       num_threads,
		       level,
		       phase,
		       (*it)->queue,
		       (*it)->begin,
		       (*it)->end,
		       nops);
      queue_dynamic.push_back(tmp);
    } // loop : it
  }   // flag == (-2)
#endif
#ifdef DEBUG_PREPARE_THREAD
  cout << "-- assigned task -- per thread = "
       << nops_per_thread << " -- ";
  for (list<C_task_seq *>::const_iterator it = task_seq_tmp.begin();
       it != task_seq_tmp.end(); ++it) {
    cout << (*it)->task_name << " ; " 
	 << (*it)->begin << " : " 
	 << (*it)->end << " / " ;
  }
  cout << "-- " << endl;
#endif
  for (ll = 0; ll < num_threads; ll++) {
#ifdef DEBUG_PREPARE_THREAD  
    cout << "thread_id " << ll 
	 << " nops_sum0 = " << nops_sum0[ll]
	 << " : nops_sum = " << nops_sum[ll] << " : ";
#endif
    queue_divided[ll].sort(C_task_seq_beginidx_smaller);
    for (list<C_task_seq *>::const_iterator it = queue_divided[ll].begin();
	 it != queue_divided[ll].end(); ++it) {
#ifdef DEBUG_PREPARE_THREAD
      cout << (*it)->task_name << " ";
#endif
      queue_static[ll % num_threads].push_back(*it);
    }
#ifdef DEBUG_PREPARE_THREAD
    cout << endl;
#endif
  } // loop : ll

  for (list<C_task_seq *>::iterator it = task_seq_tmp.begin();
       it != task_seq_tmp.end(); ++it) {
    delete (*it);    //    delete [] (*it)->task_name;
    (*it) = NULL;
  }
  task_seq_tmp.clear();
  
  // unbalance in assigned tasks is to be reduced by greedy execution but
  // unbalance from the beginning is passed to the next block
  for (int i = 0; i < num_threads; i++) {
    if (excluding_list[i]) {
      nops_sum[i] -= nops_per_thread;
    }
    else {
      nops_sum[i] = 0L;
    }
  }
  delete [] nops_sum0;
}

template<typename T, typename U>
int dimKernDense(vector<int> &singIdx, const int n,
		 const int aug_dim,
		 const U eps_machine,
		 const double eps_piv,
		 SquareBlockMatrix<T> &D,
		 T *a,
		 const bool refactorize,
		 const bool isFullPermute,
		 const bool isSym,
		 const bool verbose,
		 FILE *fp)
{
  // output == (-1) : refactorized / >= 0 : dim of the kernel
 // check dim of factorized matrix
  int sing_max = singIdx.size();
  int aug_max = sing_max + aug_dim;
  // map<int, int> aug_ind;
  list<int> aug_ind;
  vector<int> aug_ind0, aug_ind1;
  aug_ind0.resize(aug_max);
  aug_ind1.resize(aug_max);
  VectorArray<T> a_diag(aug_max);
  const T zero(0.0);
  const T one(1.0);
  const T none(-1.0);
  
  diss_printf(verbose, fp,
	      "%s %d : Schur complement form %d x %d : block size = %d\n", 
	      __FILE__, __LINE__, n, aug_max, D.block_size());
  RectBlockMatrix<T> b(n, aug_max, D.block_size());
  RectBlockMatrix<T> c(n, aug_max, D.block_size());

  ColumnMatrix<T> s;
  s.init(aug_max, aug_max);
  int *permute_d = new int[sing_max];
  vector<int> &permute = D.getPermute();
  for (int k = 0; k < sing_max; k++) {
    aug_ind.push_back(singIdx[k]);
  }
  int k0 = (n - 1);
  int k1 = sing_max;
  while ((k1 < aug_max) && k0 >= 0) {
    bool flag = true;
    for (int i = 0; i < sing_max; i++) {
      if (singIdx[i] == k0) {
	flag = false;
	break;
      }
    }
    if (flag) {
      aug_ind.push_back(k0);
      k1++; 
    }
    k0--;
  } // while
  aug_ind.sort();
  {  // for scope of int i
    int i = 0;
    for (list<int>::const_iterator it = aug_ind.begin(); it != aug_ind.end();
	 ++it, i++) {
      aug_ind0[i] = (*it);
      aug_ind1[i] = permute[(*it)];
    }
  }
  aug_ind.clear();
  diss_printf(verbose, fp, "%s %d : suspicious dimension of pivots %d + %d\n",
	      __FILE__, __LINE__, sing_max, aug_dim);
  for (int i = 0; i < aug_max; i++) {
    diss_printf(verbose, fp, "%d : %d %d %s\n", i, aug_ind0[i], aug_ind1[i],
		tostring<T>(D.diag(aug_ind0[i])).c_str());
  }
  // save diagonal entries coressponding to nullification
  for (int i = 0; i < aug_max; i++) {
    a_diag[i] = D.diag(aug_ind0[i]); // a[aug_ind0[j] * (n + 1)];
  }

  if (isSym) {
    for (int j = 0; j < aug_max; j++) {
      D.diag(aug_ind0[j]) = zero; // nullifying diagonals
      for (int i = 0; i < n; i++) {
	// access upper for symmetric matrix
	const int ii = i > aug_ind1[j] ? aug_ind1[j] : i;
	const int jj = i > aug_ind1[j] ? i : aug_ind1[j];
	b(i, j) = a[ii + jj * n];
      } // loop : i
    }   // loop : j
  }
  else {
    for (int j = 0; j < aug_max; j++) {
      D.diag(aug_ind0[j]) = zero;
      for (int i = 0; i < n; i++) {
	const int jj = aug_ind1[j];
	b(i, j) = a[i + jj * n];
	c(i, j) = a[jj + i * n];
      }
    }
  }
  for (int j = 0; j < aug_max; j++) {
    for (int i = 0; i < aug_max; i++) {
      s(i, j) = b(aug_ind1[i], j);
    }
  }

  DTRSMScale_arg<T> *tmp_arg = 
    new DTRSMScale_arg<T>(isSym,
			  &D,
			  &b, 
			  &c,
			  n,       // nrow
			  aug_max, // ncol
			  (-1),    // kblock
			  0,       // lblock
			  (-1),    // mblock
			  &aug_ind0, //singval,
			  (!isFullPermute),
			  verbose,
			  &fp,
			  -1); //dummy
  C_DTRSMScale_solve<T>(tmp_arg);
  delete tmp_arg;
  
  for (int i = 0; i < b.num_blocks_c(); i++) {
    for (int j = 0; j < b.num_blocks_c(); j++) {    
      for (int k = 0; k < b.num_blocks_r(); k++) {
	int nrow = c.nrowBlock(k);
	blas_gemm<T>(CblasTrans, CblasNoTrans,
		     c.ncolBlock(i), b.ncolBlock(j), nrow,
		     none,  // alpha = -1
		     c.addrCoefBlock(k, i), nrow,
		     b.addrCoefBlock(k, j), nrow,
		     one,   // beta = 1
		     s.addrCoefs() + (c.IndexBlock_r(i) + b.IndexBlock_c(j) * aug_max),
		     aug_max);
      }
    } // loop : i
  }   // loop : j
      // symmetrize
  if (isSym) {
    for (int j = 0; j < aug_max; j++) {
      for (int i = 0; i < j; i++) {
	s[j + i * aug_max] = s[i + j * aug_max];
      }
    }
  }
  ColumnMatrix<T> ss(aug_max, aug_max);
  ss.copy(s);
  bool flagShrinkSchur = false;
  int nn0 = sing_max;
  int *permute1 = new int[aug_max];
#if 0
  {
    double pivot = 1.0;
    double fop;
    int n0;
    n0 = 0;
    if (isSym) {
      full_ldlt_permute<T, U>(&nn0, n0, aug_max, ss.addrCoefs(), aug_max,
			      &pivot, permute1, eps_piv, &fop);
    } //   if (arg->isSym) 
    else {
      full_ldu_permute<T, U>(&nn0, n0, aug_max, ss.addrCoefs(), aug_max,
			     &pivot, permute1, eps_piv, &fop);
    }
    diss_printf(verbose, fp,
		"%s %d : %d -> %d ", __FILE__, __LINE__, sing_max, nn0);

    if (nn0 < sing_max) {
      flagShrinkSchur = true;
      diss_printf(verbose, fp,
		  "schrink the Schur complement by recursive computing\n");
      for (int i = 0; i < aug_max; i++) {
	diss_printf(verbose, fp, "%d ", permute1[i]);
      }
      diss_printf(verbose, fp, "\n");
    }
    else {
      diss_printf(verbose, fp, "\n");
    }
  }
#endif
  int aug_max1 = aug_max;
  if (flagShrinkSchur) {
    aug_max1 = aug_dim + nn0;
    int aug_max0 = aug_max - aug_max1;
    ColumnMatrix<T> sc(aug_max1, aug_max1);
    for (int j = 0; j < aug_max1; j++) {
      const int jj = permute1[aug_max - aug_max1 + j];
      for (int i = 0; i < aug_max1; i++) {
	const int ii = permute1[aug_max - aug_max1 + i];
	sc(i, j) = s(ii, jj);
      }
    }
    ColumnMatrix<T> upper(aug_max0, aug_max1);
    ColumnMatrix<T> lower(aug_max0, aug_max1);
    if (!isSym) {
      for (int j = 0; j < aug_max1; j++) {
	const int jj = permute1[aug_max - aug_max1 + j];
	for (int i = 0; i < aug_max0; i++) {
	  const int ii = permute1[i];
	  upper(i, j) = s(ii, jj);
	  lower(i, j) = s(jj, ii);
	}
      }
      blas_trsm<T>(CblasLeft, CblasLower, CblasNoTrans, CblasUnit,
		   aug_max0, aug_max1, one, ss.addrCoefs(), aug_max,
		   upper.addrCoefs(), aug_max0);
      for (int i = 0; i < aug_max0; i++) {
	const T stmp = ss(i, i);
	for (int j = 0; j < aug_max1; j++) {
	  upper(i, j) *= stmp;
	}
      }
      blas_trsm<T>(CblasLeft, CblasUpper, CblasTrans, CblasUnit,
		   aug_max0, aug_max1, one, ss.addrCoefs(), aug_max,
		   lower.addrCoefs(), aug_max0);
    }
    else {
      for (int j = 0; j < aug_max1; j++) {
	const int jj = permute1[aug_max - aug_max1 + j];
	for (int i = 0; i < aug_max0; i++) {
	  const int ii = permute1[i];
	  lower(i, j) = s(ii, jj); // copy from upper
	}
      }
      blas_trsm<T>(CblasLeft, CblasLower, CblasNoTrans, CblasUnit,
		   aug_max0, aug_max1, one, ss.addrCoefs(), aug_max,
		   lower.addrCoefs(), aug_max0);
      for (int i = 0; i < aug_max0; i++) {
	const T stmp = ss(i, i);
	  for (int j = 0; j < aug_max1; j++) {
	  upper(i, j) = lower(i, j) * stmp;
	}
      }
    }
    blas_gemm<T>(CblasTrans, CblasNoTrans, aug_max1, aug_max1, aug_max0,
		 none, lower.addrCoefs(), aug_max0, upper.addrCoefs(),
		 aug_max0, one,
		 sc.addrCoefs(), aug_max1);

    upper.free();
    lower.free();
    s.free();
    s.init(aug_max1, aug_max1);
    s.copy(sc);
    //    blas_copy<T>((aug_max1 * aug_max1), sc, 1, s, 1);
    sc.free(); 
    sing_max = nn0;
  } //   if (flagShrinkSchur) {
  ss.free();
  delete [] permute1;

  int n2;
  bool flag, flag_unsym_pivot;

  flag = ComputeDimKernel<T, U>(&n2, &flag_unsym_pivot, s.addrCoefs(), aug_max1, isSym,
				aug_dim, eps_machine, eps_piv, verbose, fp);
  diss_printf(verbose, fp, "%s %d : detected dim. of the kernel = %d %s\n",
	      __FILE__, __LINE__, n2, flag_unsym_pivot ? "full" : "symmetric");
  // 
  // restore stored diagonal entries nullified for augmented dimenison
  for (int i = 0; i < aug_max; i++) {
    D.diag(aug_ind0[i]) = a_diag[i];
  }
  if (flag_unsym_pivot) {
    return 0;     // force kernel detection failure : 15 May 2018
  }
  if ((n2 != sing_max)) {
    if (refactorize) {
      diss_printf(verbose, fp,
		  "%s %d : sing_max = %d n2 = %d -> refactorized\n",
		  __FILE__, __LINE__, sing_max, n2);
      return (-1);
    } 
  } // if (n2 != sing_max)
  
  aug_ind0.clear();
  aug_ind1.clear();
  a_diag.free();
  b.free();
  c.free(); 
  s.free();
  delete [] permute_d;

  return n2;
}

template
int dimKernDense<double, double>(vector<int> &singIdx,
				 const int n,
				 const int aug_dim,
				 const double eps_machine,
				 const double eps_piv,
				 SquareBlockMatrix<double> &D,
				 double *a,
				 const bool refactorize,
				 const bool isFullPermute,
				 const bool isSym,
				 const bool verbose,
				 FILE *fp);

template
int dimKernDense<quadruple,
		 quadruple>(vector<int> &singIdx,
			    const int n,
			    const int aug_dim,
			    const quadruple eps_machine,
			    const double eps_piv,
			    SquareBlockMatrix<quadruple> &D,
			    quadruple *a,
			    const bool refactorize,
			    const bool isFullPermute,
			    const bool isSym,
			    const bool verbose,
			    FILE *fp);

template
int dimKernDense<complex<double>,
		 double>(vector<int> &singIdx,
			 const int n,
			 const int aug_dim,
			 const double eps_mahcine,
			 const double eps_piv,
			 SquareBlockMatrix<complex<double> > &D,
			 complex<double> *a,
			 const bool refactorize,
			 const bool isFullPermute,
			 const bool isSym,
			 const bool verbose,
			 FILE *fp);

template
int dimKernDense<complex<quadruple>,
		 quadruple>(vector<int> &singIdx,
			    const int n,
			    const int aug_dim,
			    const quadruple eps_machine,
			    const double eps_piv,
			    SquareBlockMatrix<complex<quadruple> > &D,
			    complex<quadruple> *a,
			    const bool refactorize,
			    const bool isFullPermute,
			    const bool isSym,
			    const bool verbose,
			    FILE *fp);

template
int dimKernDense<float, float>(vector<int> &singIdx,
				 const int n,
				 const int aug_dim,
				 const float eps_machine,
				 const double eps_piv,
				 SquareBlockMatrix<float> &D,
				 float *a,
				 const bool refactorize,
				 const bool isFullPermute,
				 const bool isSym,
				 const bool verbose,
				 FILE *fp);

template
int dimKernDense<complex<float>,
		 float>(vector<int> &singIdx,
			 const int n,
			 const int aug_dim,
			 const float eps_mahcine,
			 const double eps_piv,
			 SquareBlockMatrix<complex<float> > &D,
			 complex<float> *a,
			 const bool refactorize,
			 const bool isFullPermute,
			 const bool isSym,
			 const bool verbose,
			 FILE *fp);
//

template<typename T>
void calc_relative_norm(double *norm_l2,  double *norm_infty, 
			const T *v, const T *u, const int dim)
{
  fprintf(stderr, "%s %d : specialized template is not yet defined.\n",
	  __FILE__, __LINE__);
}

template<>
void calc_relative_norm<double>(double *norm_l2,  double *norm_infty, 
				const double *v, const double *u, const int dim)
{
  double xtmp1 = 0.0, xtmp0 = 0.0;
  double ytmp1, ytmp0;
  double ztmp1 = 0.0, ztmp0 = 0.0;
  for (int i = 0; i < dim; i++) {
    ytmp1 = v[i];
    ytmp0 = u[i];
    if (ytmp1 < 0.0) {
      ytmp1 *= (-1.0);
    }
    if (ytmp0 < 0.0) {
      ytmp0 *= (-1.0);
    }
    xtmp0 = (ytmp0 > xtmp0) ? ytmp0 : xtmp0;
    ztmp0 += ytmp0 * ytmp0;
    
    xtmp1 = (ytmp1 > xtmp1) ? ytmp1 : xtmp1;
    ztmp1 += ytmp1 * ytmp1;
  }
  *norm_l2 = sqrt(ztmp1 / ztmp0);
  *norm_infty = xtmp1 / xtmp0;
}

template<>
void calc_relative_norm<quadruple>(double *norm_l2,  double *norm_infty, 
				   const quadruple *v,
				   const quadruple *u, const int dim)
{
  double xtmp1 = 0.0, xtmp0 = 0.0;
  double ytmp1, ytmp0;
  double ztmp1 = 0.0, ztmp0 = 0.0;
  for (int i = 0; i < dim; i++) {
    ytmp1 = quad2double(v[i]);
    ytmp0 = quad2double(u[i]);
    if (ytmp1 < 0.0) {
      ytmp1 *= (-1.0);
    }
    if (ytmp0 < 0.0) {
      ytmp0 *= (-1.0);
    }
    xtmp0 = (ytmp0 > xtmp0) ? ytmp0 : xtmp0;
    ztmp0 += ytmp0 * ytmp0;
    
    xtmp1 = (ytmp1 > xtmp1) ? ytmp1 : xtmp1;
    ztmp1 += ytmp1 * ytmp1;
  }
  *norm_l2 = sqrt(ztmp1 / ztmp0);
  *norm_infty = xtmp1 / xtmp0;
}

template<>
void calc_relative_norm<complex<double> >(double *norm_l2,  double *norm_infty, 
					  const complex<double> *v,
					  const complex<double> *u,
					  const int dim)
{
  double xtmp1 = 0.0, xtmp0 = 0.0;
  double ztmp1 = 0.0, ztmp0 = 0.0;
  for (int i = 0; i < dim; i++) {
    const double ytmp1 = std::abs(v[i]);
    const double ytmp0 = std::abs(u[i]);

    xtmp0 = ytmp0 > xtmp0 ? ytmp0 : xtmp0;
    ztmp0 += std::real(u[i] * std::conj(u[i]));
    xtmp1 = ytmp1 > xtmp1 ? ytmp1 : xtmp1;
    ztmp1 += std::real(v[i] * std::conj(v[i]));
  }
  *norm_l2 = sqrt(ztmp1 / ztmp0);
  *norm_infty = xtmp1 / xtmp0;
}

template<>
void calc_relative_norm<complex<quadruple> >(double *norm_l2,
					     double *norm_infty, 
					     const complex<quadruple> *v,
					     const complex<quadruple> *u,
					     const int dim)
{
  double xtmp1 = 0.0, xtmp0 = 0.0;
  double ztmp1 = 0.0, ztmp0 = 0.0;
  for (int i = 0; i < dim; i++) {
    quadruple vr = v[i].real();
    quadruple vi = v[i].imag();
    quadruple ur = u[i].real();
    quadruple ui = u[i].imag();
    
    const double ytmp1 = quad2double(sqrt(vr * vr + vi * vi));
    const double ytmp0 = quad2double(sqrt(ur * ur + ui * ui));

    xtmp0 = ytmp0 > xtmp0 ? ytmp0 : xtmp0;
    ztmp0 += quad2double(ur * ur + ui * ui);
    xtmp1 = ytmp1 > xtmp1 ? ytmp1 : xtmp1;
    ztmp1 += quad2double(vr * vr + vi * vi);
  }
  *norm_l2 = sqrt(ztmp1 / ztmp0);
  *norm_infty = xtmp1 / xtmp0;
}

template
void calc_relative_norm<float>(double *norm_l2, double *norm_infty, 
			       const float *v, const float *u, const int dim);

template
void calc_relative_norm<complex<float> >(double *norm_l2, double *norm_infty, 
					 const complex<float> *v,
					 const complex<float> *u,
					 const int dim);

//
template<typename T, typename Z>
void calc_relative_normscaled(double *norm_l2,  double *norm_infty, 
			       const T *v, const T *u, const Z *w,
			       const int dim)
{
  fprintf(stderr, "%s %d : specialized template is not yet defined.\n",
	  __FILE__, __LINE__);
}

template<>
void calc_relative_normscaled<double,
			      double>(double *norm_l2,  double *norm_infty, 
				      const double *v, const double *u,
				      const double *w, const int dim)
{
  double xtmp1 = 0.0, xtmp0 = 0.0;
  double ytmp1, ytmp0;
  double ztmp1 = 0.0, ztmp0 = 0.0;
  for (int i = 0; i < dim; i++) {
    ytmp0 = u[i] / w[i];
    ytmp1 = v[i] * w[i];
    if (ytmp1 < 0.0) {
      ytmp1 *= (-1.0);
    }
    if (ytmp0 < 0.0) {
      ytmp0 *= (-1.0);
    }
    xtmp0 = (ytmp0 > xtmp0) ? ytmp0 : xtmp0;
    ztmp0 += ytmp0 * ytmp0;
    
    xtmp1 = (ytmp1 > xtmp1) ? ytmp1 : xtmp1;
    ztmp1 += ytmp1 * ytmp1;
  }
  *norm_l2 = sqrt(ztmp1 / ztmp0);
  *norm_infty = xtmp1 / xtmp0;
}

template<>
void calc_relative_normscaled<quadruple,
			      quadruple>(double *norm_l2,
					 double *norm_infty, 
					 const quadruple *v,
					 const quadruple *u,
					 const quadruple *w, const int dim)
{
  double xtmp1 = 0.0, xtmp0 = 0.0;
  double ytmp1, ytmp0;
  double ztmp1 = 0.0, ztmp0 = 0.0;
  for (int i = 0; i < dim; i++) {
    ytmp0 = quad2double(u[i] / w[i]);
    ytmp1 = quad2double(v[i] * w[i]);
    if (ytmp1 < 0.0) {
      ytmp1 *= (-1.0);
    }
    if (ytmp0 < 0.0) {
      ytmp0 *= (-1.0);
    }
    xtmp0 = (ytmp0 > xtmp0) ? ytmp0 : xtmp0;
    ztmp0 += ytmp0 * ytmp0;
    
    xtmp1 = (ytmp1 > xtmp1) ? ytmp1 : xtmp1;
    ztmp1 += ytmp1 * ytmp1;
  }
  *norm_l2 = sqrt(ztmp1 / ztmp0);
  *norm_infty = xtmp1 / xtmp0;
}
template<>
void calc_relative_normscaled<complex<double>, double>(double *norm_l2,
						       double *norm_infty, 
						       const complex<double> *v,
						       const complex<double> *u,
						       const double *w,
						       const int dim)
{
  double xtmp1 = 0.0, xtmp0 = 0.0;
  double ztmp1 = 0.0, ztmp0 = 0.0;
  for (int i = 0; i < dim; i++) {
    const double ytmp0 = std::abs(u[i]) / w[i];
    const double ytmp1 = std::abs(v[i]) * w[i];
    const double ww = w[i] * w[i];
    xtmp0 = ytmp0 > xtmp0 ? ytmp0 : xtmp0;
    ztmp0 += std::real(u[i] * std::conj(u[i])) / ww;
    xtmp1 = ytmp1 > xtmp1 ? ytmp1 : xtmp1;
    ztmp1 += std::real(v[i] * std::conj(v[i])) * ww;
  }
  *norm_l2 = sqrt(ztmp1 / ztmp0);
  *norm_infty = xtmp1 / xtmp0;
}

template<>
void calc_relative_normscaled<complex<quadruple>,
			      quadruple>(double *norm_l2,
					 double *norm_infty, 
					 const complex<quadruple> *v,
					 const complex<quadruple> *u,
					 const quadruple *w,
					 const int dim)
{
  double xtmp1 = 0.0, xtmp0 = 0.0;
  double ztmp1 = 0.0, ztmp0 = 0.0;
  for (int i = 0; i < dim; i++) {
    quadruple vr = v[i].real();
    quadruple vi = v[i].imag();
    quadruple ur = u[i].real();
    quadruple ui = u[i].imag();

    const double ytmp0 = quad2double(sqrt(ur * ur + ui * ui) / w[i]);
    const double ytmp1 = quad2double(sqrt(vr * vr + vi * vi) * w[i]);
    const double ww = quad2double(w[i] * w[i]);
    xtmp0 = ytmp0 > xtmp0 ? ytmp0 : xtmp0;
    ztmp0 += quad2double(ur * ur + ui * ui) / ww;
    xtmp1 = ytmp1 > xtmp1 ? ytmp1 : xtmp1;
    ztmp1 += quad2double(vr * vr + vi * vi) * ww;
  }
  *norm_l2 = sqrt(ztmp1 / ztmp0);
  *norm_infty = xtmp1 / xtmp0;
}
//
template
void calc_relative_normscaled<float,
			      float>(double *norm_l2,  double *norm_infty, 
				     const float *v, const float *u,
				     const float *w, const int dim);

template
void calc_relative_normscaled<complex<float>,
			      float>(double *norm_l2,  double *norm_infty, 
				     const complex<float> *v,
				     const complex<float> *u,
				     const float *w, const int dim);
//

int CSR_sym2unsym(CSR_indirect *unsym,
		  const int *ptSymRows, const int *indSymCols, 
		  const int *map_eqn, const int *remap_eqn, 
		  const int dim, const bool upper_flag, bool verbose, FILE *fp)
{
  int* nbIndPerRow = new int[dim];

  //  memset(nbIndPerRow, 0, dim*sizeof(int));
  for (int i = 0; i < dim; i++) {
    nbIndPerRow[i] = 0;
  }

  for (int i = 0; i < dim; i++) {
    const int ii = remap_eqn[i];
    nbIndPerRow[i] += ptSymRows[ii + 1] - ptSymRows[ii];
    int ibegin = ptSymRows[ii] + (upper_flag ? 1 : 0);
    int iend = ptSymRows[ii + 1] + (upper_flag ? 0 : (-1));
    for (int kk = ibegin; kk < iend; kk++) {
      nbIndPerRow[map_eqn[indSymCols[kk]]]++; // 3 Jan.2014
    }
  }
  // Build unsym->ptRows array :
  // ...................
  unsym->ptRows[0] = 0;
  for (int i = 0; i < dim; i++) {
    unsym->ptRows[i + 1] = unsym->ptRows[i] + nbIndPerRow[i];
  }
  //  CHECK(unsym->ptRows[dim] == (2 * nz - dim),
  //	"error in sym2unsym() : Wrong number of non zeros elemnts in ptRows !");
	    // Allocate and fill indices columns :
  //  memset(nbIndPerRow, 0, (dim * sizeof(int)));
    for (int i = 0; i < dim; i++) {
    nbIndPerRow[i] = 0;
  }

  // for upper case, nbIndPerRow[i] keeps entries added by transposed operation
  // but for lower case counts all nonzero entries in progress 
  for (int i = 0; i < dim; i++) {
    int itmp = unsym->ptRows[i] + nbIndPerRow[i];
    const int ii = remap_eqn[i];
    for (int kk = ptSymRows[ii]; kk < ptSymRows[ii + 1]; kk++) {
      unsym->indCols[itmp] = map_eqn[indSymCols[kk]];
      unsym->indVals[itmp] = kk; //map_indcols[kk];
      itmp++;
    } // loop : kk
    if (!upper_flag) {
      nbIndPerRow[i] = itmp - unsym->ptRows[i];
    }
    //    memcpy(indCols + (ptRows[i] + nbIndPerRow[i]), 
    //	   indSymCols + ptSymRows[i],
    //	   (ptSymRows[i + 1] - ptSymRows[i]) * sizeof(int));
    int ibegin = ptSymRows[ii] + (upper_flag ? 1 : 0);
    int iend = ptSymRows[ii + 1] + (upper_flag ? 0 : (-1));
    for (int kk = ibegin; kk < iend; kk++) {
      const int j = map_eqn[indSymCols[kk]];
      const int jtmp = unsym->ptRows[j] + nbIndPerRow[j];
      unsym->indCols[jtmp] = i;
      unsym->indVals[jtmp] = kk; //map_indcols[kk];
      nbIndPerRow[j]++;
    } // loop : kk
  }   // loop : i
  delete [] nbIndPerRow;
  return unsym->ptRows[dim];
}

bool CSR_unsym2unsym(CSR_indirect *unsym,
		     const int *ptUnSymRows, const int *indUnSymCols, 
		     const int *map_eqn, const int *remap_eqn, 
		     //const int *map_indcols,
		     const int dim, const bool verbose, FILE *fp)
{
  int* nbIndPerRow = new int[dim];
  for (int i = 0; i < dim; i++) {
    const int ii = remap_eqn[i];
    nbIndPerRow[i] = ptUnSymRows[ii + 1] - ptUnSymRows[ii];
  }
  unsym->ptRows[0] = 0;
  for (int i = 0; i < dim; i++) {
    unsym->ptRows[i + 1] = unsym->ptRows[i] + nbIndPerRow[i];
  }

  for (int i = 0; i < dim; i++) { // running over new index
    const int ii = remap_eqn[i];  // access to the original data of node
    int itmp = unsym->ptRows[i];
    for (int kk = ptUnSymRows[ii]; kk < ptUnSymRows[ii + 1]; kk++, itmp++) {
      const int jj = indUnSymCols[kk];
      const int j = map_eqn[jj];
      unsym->indCols[itmp] = j;    // access to the new data
      unsym->indVals[itmp] = kk; // access to the orignal data of unsym->indCols
      bool found = false;
      for (int mm = ptUnSymRows[jj]; mm < ptUnSymRows[jj + 1]; mm++) {
	if (indUnSymCols[mm] == ii) { // based on the original CSR data
	  unsym->indVals_unsym[itmp] = mm; // map_indcols[mm];
	  found = true;
	  break;
	}
      }
      if (!found) {
	diss_printf(verbose, fp,
		    "%s %d : symmetric place of (%d, %d) -> (%d %d) not found\n",
		  __FILE__, __LINE__, ii, jj, i, j);
	return false;
      }
    }
  }
  delete [] nbIndPerRow;
  return true;
}



void swap_queues_n(vector <C_task *> &queue,
		   vector <int> &queue_index,
		   const int ii, 
		   const int jj,
		   const int n,
		   vector <C_task *> &tmp,
		   vector <int> &tmp_index)
{
  for (int l = 0; l < n; l++) {
    tmp[l] = queue[ii * n + l];
    tmp_index[l] = queue_index[ii * n + l];
  }
  for (int l = 0; l < n; l++) {
    queue[ii * n + l] = queue[(ii + 1) * n + l];
    queue_index[ii * n + l] = queue_index[(ii + 1) * n + l];
  }
  for (int l = 0; l < n; l++) {
    queue[jj * n + l] = tmp[l];
    queue_index[jj * n + l] = tmp_index[l];
  }
}

bool compare_task_name(C_task *first, C_task *second) {
  // 'g', 'h', 'i', 'j' should be in higher order
  if (first->task_name[0] >= 'h') {
    if (second->task_name[0] >= 'h') {
      return (strcmp(first->task_name, second->task_name) >=0 ? false : true);
    }
    else {
      return true;
    }
  }
  if (second->task_name[0] >= 'h') {
    return false;
  }
  return (strcmp(first->task_name, second->task_name) >=0 ? false : true);
}

// =====================================================================
template<typename T>
int count_diag_negative_real(SquareBlockMatrix<T>& Diag)
{
  int count = 0;
  const int dim_diag = Diag.dimension();
  for (int i = 0; i < dim_diag; i++) {
    const T xtmp = Diag.diag(i);
    if (xtmp < T(0.0)) {
      count++;
    }
  }
  return count;
}

template<typename T>
int count_diag_negative_complex(SquareBlockMatrix<complex<T> >& Diag)
{
  int count = 0;
  const int dim_diag = Diag.dimension();
  for (int i = 0; i < dim_diag; i++) {
    const complex<T> xtmp = Diag.diag(i);
    if ((xtmp.real() < T(0.0)) && xtmp.imag() == T(0.0)) {
      count++;
    }
  }
  return count;
}

template<typename T>
int count_diag_negative(SquareBlockMatrix<T>& Diag)
{
   fprintf(stderr, "%s %d : specialized template is not yet defined.\n",
          __FILE__, __LINE__);
   return (-1);
}

template<>
int count_diag_negative<double>(SquareBlockMatrix<double>& Diag)
{
  return count_diag_negative_real<double>(Diag);
}

template<>
int count_diag_negative<quadruple>(SquareBlockMatrix<quadruple>& Diag)
{
  return count_diag_negative_real<quadruple>(Diag);
}

template<>
int count_diag_negative<complex<double> >(SquareBlockMatrix<complex<double> >& Diag)
{
  return count_diag_negative_complex<double>(Diag);
}

template<>
int count_diag_negative<complex<quadruple> >(SquareBlockMatrix<complex<quadruple> >& Diag)
{
  return count_diag_negative_complex<quadruple>(Diag);
}

template<>
int count_diag_negative<float>(SquareBlockMatrix<float>& Diag)
{
  return count_diag_negative_real<float>(Diag);
}


template<>
int count_diag_negative<complex<float> >(SquareBlockMatrix<complex<float> >& Diag)
{
  return count_diag_negative_complex<float>(Diag);
}
//
template<typename T>
int count_diag_negative(SubSquareMatrix<T>& Diag)
{
   fprintf(stderr, "%s %d : specialized template is not yet defined.\n",
	  __FILE__, __LINE__);
   return (-1);
}

template<>
int count_diag_negative<double>(SubSquareMatrix<double>& Diag)
{
  int count = 0;
  const int dim_diag = Diag.dimension();
  for (int i = 0; i < dim_diag; i++) {
    const double xtmp = Diag(i, i);
    if (xtmp < 0.0) {
      count++;
    }
  }
  return count;
}

template<>
int count_diag_negative<quadruple>(SubSquareMatrix<quadruple>& Diag)
{
  int count = 0;
  const int dim_diag = Diag.dimension();
  for (int i = 0; i < dim_diag; i++) {
    const double xtmp = quad2double(Diag(i, i));
    if (xtmp < 0.0) {
      count++;
    }
  }
  return count;
}

template<>
int count_diag_negative<complex<quadruple> >(SubSquareMatrix<complex<quadruple> >& Diag)
{
  int count = 0;
  const int dim_diag = Diag.dimension();
  for (int i = 0; i < dim_diag; i++) {
    const complex<quadruple> &xtmp = Diag(i, i);
    if (quad2double(xtmp.imag()) < 0.0 && quad2double(xtmp.real()) == 0.0) {
      count++;
    }
  }
  return count;
}

template<>
int count_diag_negative<complex<double> >(SubSquareMatrix<complex<double> >& Diag)
{
  int count = 0;
  const int dim_diag = Diag.dimension();
  for (int i = 0; i < dim_diag; i++) {
    const complex<double> xtmp = Diag(i, i);
    if (xtmp.imag() < 0.0 && xtmp.real() == 0.0) {
      count++;
    }
  }
  return count;
}

template<>
int count_diag_negative<float>(SubSquareMatrix<float>& Diag)
{
  int count = 0;
  const int dim_diag = Diag.dimension();
  for (int i = 0; i < dim_diag; i++) {
    const float xtmp = Diag(i, i);
    if (xtmp < 0.0) {
      count++;
    }
  }
  return count;
}

template<>
int count_diag_negative<complex<float> >(SubSquareMatrix<complex<float> >& Diag)
{
  int count = 0;
  const int dim_diag = Diag.dimension();
  for (int i = 0; i < dim_diag; i++) {
    const complex<float> xtmp = Diag(i, i);
    if (xtmp.imag() < 0.0 && xtmp.real() == 0.0) {
      count++;
    }
  }
  return count;
}

//
template<typename T, typename U>
void C_SparseNumFact(void *arg_)
{
  C_SparseNumFact_arg<T, U> *arg = (C_SparseNumFact_arg<T, U> *)arg_;
  TridiagBlockMatrix<T, U> **tridiag = arg->tridiag;
  int colors = arg->colors;
  int nrow = arg->nrow;
  double eps_pivot = *(arg->eps_pivot);
  const bool kernel_detection = *(arg->kernel_detection);
  const bool higher_precision = *(arg->higher_precision);
  const int dim_aug_kern = *(arg->dim_aug_kern);
  const U eps_machine = *(arg->eps_machine);
  T *coefs = arg->coefs;

  double *pivot = arg->pivot;  // all diagonal entries of sparse matrices 
			        // equal to 1, set at the initialization

  SquareBlockMatrix<T>& D = *(arg->D);
  vector<int> &list_sing = D.getSingIdx();
  double nopd;
  const bool verbose = arg->verbose;
  FILE *fp = *(arg->fp);
  elapsed_t t0, t1;
  get_realtime(&t0);
  int nsing = 0;
  bool detected = true;
  nopd = 0.0;
  *pivot = 1.0; // matrix is scaled so that the diagonal entries take 1
  for (int i = 0; i < colors; i++) {
    double pivot1;
    double nopd1;
    tridiag[i]->NumericFact(coefs, eps_pivot,
			    &pivot1,
			    kernel_detection,
			    higher_precision,
			    dim_aug_kern,
			    eps_machine,
			    &nopd1);
    nopd += nopd1;
    //    tridiag[i0->SingularNode(list_sing);
    if (!tridiag[i]->detected()) {
      detected = false; // tridiag[i]->nsing() == (-1);     
    }
    nsing += tridiag[i]->nsing();
    *pivot = *pivot < pivot1 ? *pivot : pivot1;
  }
  list_sing.resize(nsing);
  {
    int j = 0;
    for (int i = 0; i < colors; i++) {
      vector<int> list_sing_tmp;
      tridiag[i]->SingularNode(list_sing_tmp);
      for (vector<int>::const_iterator it = list_sing_tmp.begin();
	   it != list_sing_tmp.end(); ++it, j++) {
	list_sing[j] = (*it);
      }
    }
  }
	
  get_realtime(&t1);
  diss_printf(verbose, fp, "%s %d : %d : pivot = %g n0 = %d detected = %s\n",
	      __FILE__, __LINE__,
	      arg->nb, *pivot, nsing, detected ? "true" : "false");
  
  D.set_lastPivot(*pivot);
  D.set_KernelDetected(detected);
  D.set_rank(nrow - nsing);

  //*(arg->pivot) = pivot;
  *(arg->nopd) = (long)nopd; // ?? uninitialized ??
}

template
void C_SparseNumFact<double, double>(void *arg_);

template
void C_SparseNumFact<quadruple, quadruple>(void *arg_);

template
void C_SparseNumFact<complex<double>, double>(void *arg_);

template
void C_SparseNumFact<complex<quadruple>,  quadruple>(void *arg_);

template
void C_SparseNumFact<float, float>(void *arg_);

template
void C_SparseNumFact<complex<float>, float>(void *arg_);
//


template<typename T, typename U>
void C_SparseLocalSchur(void *arg_)
{
  C_SparseNumFact_arg<T, U> *arg = (C_SparseNumFact_arg<T, U> *)arg_;

  TridiagBlockMatrix<T, U> **tridiag = arg->tridiag;
  int colors = arg->colors;
  int *color_mask = arg->color_mask;
  int nrow = arg->nrow;
  int ncol = arg->ncol;
  T *coefs = arg->coefs;
  int *prow1 = arg->csr_offdiag->ptRows;
  int *indcols1 = arg->csr_offdiag->indCols;
  int *indvals1 = arg->csr_offdiag->indVals;
  int *indvals2 = arg->isSym ? (int *)NULL : arg->csr_offdiag->indVals_unsym;

  elapsed_t *tt = arg->tt;
  double *nops = new double[3];
  //
  SquareBlockMatrix<T> &localSchur = *(arg->localSchur);
  localSchur.allocate();
  double nopsum = 0.0;
  for (int i = 0; i < colors; i++) {
    tridiag[i]->ComputeSchur(nrow, color_mask,
			     ncol, prow1, indcols1, indvals1, indvals2,
			     coefs, SIZE_B1, localSchur, nops, tt);
                                      // get the excat flop by direct counting
    nopsum += nops[0] + nops[2];
  }
  *(arg->nops) = (long)nopsum;
  delete [] nops;
}

template
void C_SparseLocalSchur<double, double>(void *arg_);

template
void C_SparseLocalSchur<quadruple, quadruple>(void *arg_);

template
void C_SparseLocalSchur<complex<double>, double>(void *arg_);

template
void C_SparseLocalSchur<complex<quadruple>, quadruple>(void *arg_);

template
void C_SparseLocalSchur<float, float>(void *arg_);

template
void C_SparseLocalSchur<complex<float>, float>(void *arg_);
//


template<>
void dump_matrix<double>(FILE *fp, const int nrow, double *a)
{
  for (int i = 0; i < nrow; i++) {
    fprintf(fp, "%d : ", i);
    for (int j = i; j < nrow; j++) {
      fprintf(fp, "%g ", a[i + nrow * j]);
    }
    fprintf(fp, "\n");
  }
}

template<typename T>
void dump_matrix(FILE *fp, const int nrow, T *a)
{
  fprintf(stderr, "%s %d : general case is not defined\n", __FILE__, __LINE__);
}

template
void dump_matrix<quadruple>(FILE *fp, const int nrow, 
			    quadruple *a);

template
void dump_matrix<complex<double> >(FILE *fp, const int nrow, 
				    complex<double> *a);

template
void dump_matrix<complex<quadruple> >(FILE *fp, const int nrow, 
				      complex<quadruple> *a);

template
void dump_matrix<complex<float> >(FILE *fp, const int nrow, 
				    complex<float> *a);

template<>
void dump_matrix<double>(FILE *fp, const int nrow, const int ncol, double *a)
{
  for (int i = 0; i < nrow; i++) {
    fprintf(fp, "%d : ", i);
    for (int j = 0; j < ncol; j++) {
      fprintf(fp, "%g ", a[i + nrow * j]);
    }
    fprintf(fp, "\n");
  }
}

template<>
void dump_matrix<float>(FILE *fp, const int nrow, const int ncol, float *a)
{
  for (int i = 0; i < nrow; i++) {
    fprintf(fp, "%d : ", i);
    for (int j = 0; j < ncol; j++) {
      fprintf(fp, "%g ", a[i + nrow * j]);
    }
    fprintf(fp, "\n");
  }
}

template<typename T>
void dump_matrix(FILE *fp, const int nrow, const int ncol, T *a)
{
  fprintf(stderr, "%s %d : general case is not defined\n", __FILE__, __LINE__);
}

template
void dump_matrix<quadruple>(FILE *fp, const int nrow, const int ncol, 
			    quadruple *a);

template
void dump_matrix<complex<double> >(FILE *fp, const int nrow, const int ncol, 
				    complex<double> *a);
template
void dump_matrix<complex<quadruple> >(FILE *fp, const int nrow, const int ncol, 
				      complex<quadruple> *a);

template
void dump_matrix<complex<float> >(FILE *fp, const int nrow, const int ncol, 
				    complex<float> *a);

template<>
void dump_matrix<double>(FILE *fp, const int kk, 
			 const int nrow, const int ncol, const int nn,
			 double *a)
{
  for (int i = 0; i < nrow; i++) {
    fprintf(fp, "%d : ", i);
    for (int j = 0; j < ncol; j++) {
      fprintf(fp, "%g ", a[kk + i + nn * j]);
    }
    fprintf(fp, "\n");
  }
}

template<>
void dump_matrix<float>(FILE *fp, const int kk, 
			 const int nrow, const int ncol, const int nn,
			 float *a)
{
  for (int i = 0; i < nrow; i++) {
    fprintf(fp, "%d : ", i);
    for (int j = 0; j < ncol; j++) {
      fprintf(fp, "%g ", a[kk + i + nn * j]);
    }
    fprintf(fp, "\n");
  }
}

template<typename T>
void dump_matrix(FILE *fp, const int kk, 
		 const int nrow, const int ncol, const int nn, T *a)
{
  fprintf(stderr, "%s %d : general case is not defined\n", __FILE__, __LINE__);
}

template
void dump_matrix<quadruple>(FILE *fp, const int kk, 
			    const int nrow, const int ncol, const int nn,
			    quadruple *a);
template
void dump_matrix<complex<double> >(FILE *fp, const int kk, 
				   const int nrow, const int ncol, const int nn,
				    complex<double> *a);

template
void dump_matrix<complex<quadruple> >(FILE *fp, const int kk, 
				      const int nrow, const int ncol,
				      const int nn,
				      complex<quadruple> *a);
template
void dump_matrix<complex<float> >(FILE *fp, const int kk, 
				   const int nrow, const int ncol, const int nn,
				    complex<float> *a);

template<>
void dump_matrix<double>(FILE *fp, 
			 RectBlockMatrix<double> &a)
{
  for (int i = 0; i < a.dimension_r(); i++) {
    fprintf(fp, "%d : ", i);
    for (int j = 0; j < a.dimension_c(); j++) {
      fprintf(fp, "%g ", a(i, j));
    }
    fprintf(fp, "\n");
  }
}

template<>
void dump_matrix<float>(FILE *fp, 
			 RectBlockMatrix<float> &a)
{
  for (int i = 0; i < a.dimension_r(); i++) {
    fprintf(fp, "%d : ", i);
    for (int j = 0; j < a.dimension_c(); j++) {
      fprintf(fp, "%g ", a(i, j));
    }
    fprintf(fp, "\n");
  }
}

template<typename T>
void dump_matrix(FILE *fp, 
		 RectBlockMatrix<T> &a)
{
  fprintf(stderr, "%s %d : general case is not defined\n", __FILE__, __LINE__);
}

template<>
void dump_matrix<double>(FILE *fp, 
			 SquareBlockMatrix<double> &a)
{
  // if (a.isSym()) {
  if (1) { // 16 Jan.debug
    for (int i = 0; i < a.dimension(); i++) {
      fprintf(fp, "%d : ", i);
      for (int j = i; j < a.dimension() ; j++) {
	fprintf(fp, "%16.8e ", a(i, j));
      }
      fprintf(fp, "\n");
    }
  }
  else {
    for (int i = 0; i < a.dimension(); i++) {
      fprintf(fp, "%d : ", i);
      for (int j = 0; j < a.dimension(); j++) {
	fprintf(fp, "%16.8e ", a(i, j));
      }
      fprintf(fp, "\n");
    }
  }
}
template
void dump_matrix<quadruple>(FILE *fp, RectBlockMatrix<quadruple> &a);

template
void dump_matrix<complex<double> >(FILE *fp, 
				   RectBlockMatrix<complex<double> > &a);

template
void dump_matrix<complex<quadruple> >(FILE *fp, 
				      RectBlockMatrix<complex<quadruple> > &a);

template
void dump_matrix<complex<float> >(FILE *fp, 
				   RectBlockMatrix<complex<float> > &a);
//

template<typename T>
void dump_matrix(FILE *fp, 
		 SquareBlockMatrix<T> &a)
{
  fprintf(stderr, "%s %d : general case is not defined\n", __FILE__, __LINE__);
}

template
void dump_matrix<quadruple>(FILE *fp, SquareBlockMatrix<quadruple> &a);

template
void dump_matrix<complex<double> >(FILE *fp, 
				   SquareBlockMatrix<complex<double> > &a);

template
void
dump_matrix<complex<quadruple> >(FILE *fp, 
				 SquareBlockMatrix<complex<quadruple> > &a);

template
void dump_matrix<complex<float> >(FILE *fp, 
				   SquareBlockMatrix<complex<float> > &a);

//
template<>
void dump_matrix<double>(FILE *fp, const int nrow, const int nnz, int *prow,
			 int *indcols, int *indvals, double *a)
{
  //  verify_nan<double>(fp, nnz, a);

  for (int i = 0; i < nrow; i++) {
    if (prow[i + 1] > prow[i]) {
      fprintf(fp, "%d : ", i);
      for (int j = prow[i]; j < prow[i + 1]; j++) {
	fprintf(fp, "%d:%d:%g ", indcols[j], indvals[j], a[indvals[j]]);
      }
      fprintf(fp, "\n");
    }
  }
}

template<typename T>
void dump_matrix(FILE *fp, const int nrow, const int nnz, int *prow,
		 int *indcols, int *indvals, T *a)
{
  fprintf(stderr, "%s %d : general case is not defined\n", __FILE__, __LINE__);
}

template
void dump_matrix<complex<double> >(FILE *fp, const int nrow, const int nnz, 
				   int *prow, 
				   int *indcols, int *indvals, 
				   complex<double> *a);
//

template<typename T>
void C_FillMatrix_diag(void *arg_)
{
  C_FillMatrix_arg<T> *arg = (C_FillMatrix_arg<T> *)arg_;

  SquareBlockMatrix<T>& D = *arg->D;
  CSR_indirect csr_diag = *arg->csr_diag;
  T *coefs = arg->coefs;

  D.ZeroClear();
  for (int i = 0; i < csr_diag.n; i++) {
    for (int k = csr_diag.ptRows[i]; k < csr_diag.ptRows[i + 1]; k++) {
      const int j = csr_diag.indCols[k];
      D(i, j) = coefs[csr_diag.indVals[k]];
    }
  }
  // 
}

template
void C_FillMatrix_diag<double>(void *arg_);

template
void C_FillMatrix_diag<quadruple>(void *arg_);

template
void C_FillMatrix_diag<complex<double> >(void *arg_);

template
void C_FillMatrix_diag<complex<quadruple> >(void *arg_);

template
void C_FillMatrix_diag<float>(void *arg_);

template
void C_FillMatrix_diag<complex<float> >(void *arg_);

//

template<typename T>
void C_FillMatrix_offdiag(void *arg_)
{
  C_FillMatrix_arg<T> *arg = (C_FillMatrix_arg<T> *)arg_;

  RectBlockMatrix<T>& upper = *arg->upper;
  CSR_indirect csr_offdiag = *arg->csr_offdiag;
  T *coefs = arg->coefs;

  upper.ZeroClear();
  for (int i = 0; i < csr_offdiag.n; i++) {
    for (int k = csr_offdiag.ptRows[i]; k < csr_offdiag.ptRows[i + 1]; k++) {
      const int j = csr_offdiag.indCols[k];
      // access to (i, j) = i + j * n is not efficient
      upper(i, j) = coefs[csr_offdiag.indVals[k]];
    }
  }
  if (!arg->isSym) {
    RectBlockMatrix<T>& lower = *arg->lower;
    lower.ZeroClear();
    for (int i = 0; i < csr_offdiag.n; i++) {
      for (int k = csr_offdiag.ptRows[i]; k < csr_offdiag.ptRows[i + 1]; k++) {
	const int j = csr_offdiag.indCols[k];
	// access to (i, j) = i + j * n is not efficient
	lower(i, j) = coefs[csr_offdiag.indVals_unsym[k]];
      }
    }
  }
}

template
void C_FillMatrix_offdiag<double>(void *arg_);

template
void C_FillMatrix_offdiag<quadruple>(void *arg_);

template
void C_FillMatrix_offdiag<complex<double> >(void *arg_);

template
void C_FillMatrix_offdiag<complex<quadruple> >(void *arg_);

template
void C_FillMatrix_offdiag<float>(void *arg_);

template
void C_FillMatrix_offdiag<complex<float> >(void *arg_);
//

template<typename T>
void DSchurGEMM_diag(void *arg_)
{
  const T zero(0.0);
  const T one(1.0);
  DSchurGEMM_arg<T> *arg = (DSchurGEMM_arg<T> *)arg_;
  SquareBlockMatrix<T> *localSchur = arg->localSchur;
  localSchur->allocateBlock(arg->i_block, arg->i_block);
  T *s = localSchur->addrCoefBlock(arg->i_block, arg->i_block);
  const int nrow = localSchur->nrowBlock(arg->i_block, arg->i_block);

  //  alpha = 1;
  // to replace by symmetric dgemm which compute only upper part of the matrix 
  if (arg->isSym) {
    for (int k = 0; k < arg->lower->num_blocks_r(); k++) {
      const T beta = ((k == 0) ? zero : one);
      const T *ptLt = arg->lower->addrCoefBlock(k, arg->i_block);
      const T *ptU = arg->upper->addrCoefBlock(k, arg->i_block);
      int nnrow = arg->lower->nrowBlock(k);
      C_gemm_symm(nrow, // arg->block_nrow, 
		  nnrow, // arg->nrow, 
		  one,   // alpha
		  ptLt, 
		  nnrow, // arg->nrow,
		  ptU, // ptL, // arg->ptL, 
		  nnrow, // arg->nrow,
		  beta,
		  s, // arg->s + (arg->i_row + arg->i_row * arg->ncol),
		  nrow);
    } // loop : k
  }
  else {
    for (int k = 0; k < arg->lower->num_blocks_r(); k++) {
      const T beta = ((k == 0) ? zero : one);
      const T *ptLt = arg->lower->addrCoefBlock(k, arg->i_block);
      const T *ptU = arg->upper->addrCoefBlock(k, arg->i_block);
      int nnrow = arg->lower->nrowBlock(k);
      blas_gemm<T>(CblasTrans, CblasNoTrans, 
		   //	arg->block_nrow, arg->block_nrow, 
		   nrow, nrow,
		   nnrow, // arg->nrow,
		   one, // alpha
		   ptLt, // ptUt, // arg->ptUt, 
		   nnrow, // arg->nrow,
		   ptU, // ptL, // arg->ptL, 
		   nnrow, // arg->nrow,
		   beta,
		   s, // arg->s + (arg->i_row + arg->i_row * arg->ncol),
		   nrow);
    } // loop : k
  }
}

template
void DSchurGEMM_diag<double>(void *arg_);

template
void DSchurGEMM_diag<quadruple>(void *arg_);

template
void DSchurGEMM_diag<complex<double> >(void *arg_);

template
void DSchurGEMM_diag<complex<quadruple> >(void *arg_);

template
void DSchurGEMM_diag<float>(void *arg_);

template
void DSchurGEMM_diag<complex<float> >(void *arg_);
//

template<typename T>
void DSchurGEMM_diag_two(void *arg_)
{
  DSchurGEMM_two_arg<T> *arg = (DSchurGEMM_two_arg<T> *)arg_;
  if (arg->isSkip) {
    return;
  }

  SquareBlockMatrix<T> *fatherDiag = arg->localSchur;
  T *s = fatherDiag->addrCoefBlock(arg->i_block, arg->i_block);
  const int nrow = fatherDiag->nrowBlock(arg->i_block, arg->i_block);
  const T none(-1.0);
  const T one(1.0);
  // to replace by symmetric dgemm which compute only upper part of the matrix 
  if (arg->isSym) {
    for (int k = 0; k < arg->lower0->num_blocks_r(); k++) {
      const T *ptLt0 = arg->lower0->addrCoefBlock(k, arg->i_block);
      const T *ptU0 = arg->upper0->addrCoefBlock(k, arg->i_block);
      int nnrow = arg->lower0->nrowBlock(k);
      C_gemm_symm(nrow, // arg->block_nrow, 
		  nnrow, // arg->nrow0, 
		  none,  // alpha
		  ptLt0, nnrow, // arg->nrow0,
		  ptU0, nnrow, // arg->nrow0,
		  one,   // beta
		  s, nrow);
    }
    for (int k = 0; k < arg->lower1->num_blocks_r(); k++) {
      const T *ptLt1 = arg->lower1->addrCoefBlock(k, arg->i_block);
      const T *ptU1 = arg->upper1->addrCoefBlock(k, arg->i_block);
      int nnrow = arg->lower1->nrowBlock(k);
      C_gemm_symm(nrow, // arg->block_nrow, 
		  nnrow, // arg->nrow1, 
		  none, // alpha
		  ptLt1, nnrow, //arg->nrow1,
		  ptU1, nnrow, // arg->nrow1,
		  one, // beta
		  s, nrow);
    }
  }
  else {
    for (int k = 0; k < arg->lower0->num_blocks_r(); k++) {
      const T *ptLt0 = arg->lower0->addrCoefBlock(k, arg->i_block);
      const T *ptU0 = arg->upper0->addrCoefBlock(k, arg->i_block);
      int nnrow = arg->lower0->nrowBlock(k);
      blas_gemm<T>(CblasTrans, CblasNoTrans, 
		   // arg->block_nrow, arg->block_nrow, 
		   nrow, nrow,
		   nnrow, // arg->nrow0,
		   none,  // alpha
		   ptLt0, nnrow, // arg->nrow0, 
		   ptU0, nnrow, // arg->nrow0,
		   one, // beta,
		   s, nrow);
    }
    for (int k = 0; k < arg->lower1->num_blocks_r(); k++) {
      const T *ptLt1 = arg->lower1->addrCoefBlock(k, arg->i_block);
      const T *ptU1 = arg->upper1->addrCoefBlock(k, arg->i_block);
      int nnrow = arg->lower1->nrowBlock(k);
      blas_gemm<T>(CblasTrans, CblasNoTrans, 
		   //arg->block_nrow, arg->block_nrow, 
		   nrow, nrow,
		   nnrow, // arg->nrow1,
		   none, // alpha
		   ptLt1, nnrow,  // arg->nrow1, 
		   ptU1, nnrow, //arg->nrow1,
		   one, // beta,
		   s, nrow);
    }
  }
}

template
void DSchurGEMM_diag_two<double>(void *arg_);

template
void DSchurGEMM_diag_two<quadruple>(void *arg_);

template
void DSchurGEMM_diag_two<complex<double> >(void *arg_);

template
void DSchurGEMM_diag_two<complex<quadruple> >(void *arg_);

template
void DSchurGEMM_diag_two<float>(void *arg_);

template
void DSchurGEMM_diag_two<complex<float> >(void *arg_);

//

template<typename T>
void DSchurGEMM_offdiag(void *arg_)
{
  const T zero(0.0);
  const T one(1.0);
  DSchurGEMM_arg<T> *arg = (DSchurGEMM_arg<T> *)arg_;
  SquareBlockMatrix<T> *localSchur = arg->localSchur;
  localSchur->allocateBlock(arg->i_block, arg->j_block);
  T *s = localSchur->addrCoefBlock(arg->i_block, arg->j_block);
  const int nrow = localSchur->nrowBlock(arg->i_block, arg->j_block);
  const int ncol = localSchur->ncolBlock(arg->i_block, arg->j_block);
  // alpha = 1
  // block_nrow * block_ncol : S(i,j) = L(i)^T U(j)
  if (arg->isTrans) {
    for (int k = 0; k < arg->lower->num_blocks_r(); k++) {
      const T beta = ((k == 0) ? zero : one);
      const T *ptLt = arg->lower->addrCoefBlock(k, arg->j_block);
      const T *ptU = arg->upper->addrCoefBlock(k, arg->i_block);
      int nnrow = arg->lower->nrowBlock(k);
      blas_gemm<T>(CblasTrans, CblasNoTrans, 
		   nrow, ncol,
		   nnrow, //arg->nrow,
		   one, // alpha,
		   ptLt, 
		   nnrow, // arg->nrow,
		   ptU, 
		   nnrow, // arg->nrow,
		   beta,
		   s, 
		   nrow);
    }
  }
  else {
    for (int k = 0; k < arg->lower->num_blocks_r(); k++) {
      const T beta = ((k == 0) ? zero : one);
      const T *ptLt = arg->lower->addrCoefBlock(k, arg->i_block);
      const T *ptU = arg->upper->addrCoefBlock(k, arg->j_block);
      int nnrow = arg->lower->nrowBlock(k);
      blas_gemm<T>(CblasTrans, CblasNoTrans, 
		   nrow, ncol,
		   nnrow, // arg->nrow,
		   one, // alpha
		   ptLt, 
		   nnrow, // arg->nrow,
		   ptU, 
		   nnrow, // arg->nrow,
		   beta,
		   s, 
		   nrow);
    }
  }
}

template
void DSchurGEMM_offdiag<double>(void *arg_);

template
void DSchurGEMM_offdiag<quadruple>(void *arg_);

template
void DSchurGEMM_offdiag<complex<double> >(void *arg_);

template
void DSchurGEMM_offdiag<complex<quadruple> >(void *arg_);

template
void DSchurGEMM_offdiag<float>(void *arg_);

template
void DSchurGEMM_offdiag<complex<float> >(void *arg_);
//

template<typename T>
void DSchurGEMM_offdiag_two(void *arg_)
{
  DSchurGEMM_two_arg<T> *arg = (DSchurGEMM_two_arg<T> *)arg_;
  if (arg->isSkip) {
    return;
  }
  SquareBlockMatrix<T> *fatherDiag = arg->localSchur;  
  T *s = fatherDiag->addrCoefBlock(arg->i_block, arg->j_block);
  const int nrow = fatherDiag->nrowBlock(arg->i_block, arg->j_block);
  const int ncol = fatherDiag->ncolBlock(arg->i_block, arg->j_block);
  const T none(-1.0);
  const T one(1.0);
  // block_nrow * block_ncol : S(i,j) = L(i)^T U(j)
  if (arg->isTrans) {
    for (int k = 0; k < arg->lower0->num_blocks_r(); k++) {
      const T *ptLt0 = arg->lower0->addrCoefBlock(k, arg->j_block);
      const T *ptU0 = arg->upper0->addrCoefBlock(k, arg->i_block);
      int nnrow = arg->lower0->nrowBlock(k);    
      blas_gemm<T>(CblasTrans, CblasNoTrans, 
		   nrow, ncol,
		   nnrow, // arg->nrow0,
		   none, //alpha,
		   ptLt0, nnrow, // arg->nrow0,
		   ptU0, nnrow, // arg->nrow0,
		   one, //beta,
		   s, nrow);
    }
    for (int k = 0; k < arg->lower1->num_blocks_r(); k++) {
      const T *ptLt1 = arg->lower1->addrCoefBlock(k, arg->j_block);
      const T *ptU1 = arg->upper1->addrCoefBlock(k, arg->i_block);
      int nnrow = arg->lower1->nrowBlock(k);
      blas_gemm<T>(CblasTrans, CblasNoTrans, 
		   nrow, ncol,
		   nnrow, // arg->nrow1,
		   none, // alpha,
		   ptLt1, nnrow, // arg->nrow1,
		   ptU1, nnrow, // arg->nrow1,
		   one, // beta,
		   s, nrow);
    }
  }
  else {
    for (int k = 0; k < arg->lower0->num_blocks_r(); k++) {
      const T *ptLt0 = arg->lower0->addrCoefBlock(k, arg->i_block);
      const T *ptU0 = arg->upper0->addrCoefBlock(k, arg->j_block);
      int nnrow = arg->lower0->nrowBlock(k);    
      blas_gemm<T>(CblasTrans, CblasNoTrans, 
		   nrow, ncol,
		   nnrow, // arg->nrow0,
		   none, //alpha,
		   ptLt0, nnrow, // arg->nrow0,
		   ptU0, nnrow, // arg->nrow0,
		   one, //beta,
		   s, nrow);
    }
    for (int k = 0; k < arg->lower1->num_blocks_r(); k++) {
      const T *ptLt1 = arg->lower1->addrCoefBlock(k, arg->i_block);
      const T *ptU1 = arg->upper1->addrCoefBlock(k, arg->j_block);
      int nnrow = arg->lower1->nrowBlock(k);
      blas_gemm<T>(CblasTrans, CblasNoTrans, 
		   nrow, ncol,
		   nnrow, // arg->nrow1,
		   none, // alpha,
		   ptLt1, nnrow, // arg->nrow1,
		   ptU1, nnrow, // arg->nrow1,
		   one, //beta,
		   s, nrow);
    }
  }
}

template
void DSchurGEMM_offdiag_two<double>(void *arg_);

template
void DSchurGEMM_offdiag_two<quadruple>(void *arg_);

template
void DSchurGEMM_offdiag_two<complex<double> >(void *arg_);

template
void DSchurGEMM_offdiag_two<complex<quadruple> >(void *arg_);

template
void DSchurGEMM_offdiag_two<float>(void *arg_);

template
void DSchurGEMM_offdiag_two<complex<float> >(void *arg_);
//

template<typename T>
void C_DTRSMScale_diag_upper(void *arg_)
{
  DTRSMScale_arg<T> *arg = (DTRSMScale_arg<T> *)arg_;
  int ncol = arg->ncol;
  vector<int> &permute = arg->LDLt->getPermute();
  const int k = arg->kblock;
  const int kk = arg->LDLt->IndexBlock(k);
  vector<int>* singLocNodes0 = arg->singLocNodes0;
  const int n0 = singLocNodes0->size();
  T *LDLt = arg->LDLt->addrCoefBlock(k, k);
  T *upper = arg->upper->addrCoefBlock(k, arg->lblock);
  //  T *lower; 
  const int nrow = arg->LDLt->nrowBlock(k, k);
  const T one(1.0);
  const T zero(0.0);
  int nn0 = 0;
  if (arg->isSym) {
    arg->lower->allocateBlock(k, arg->lblock);
    T *lower = arg->lower->addrCoefBlock(k, arg->lblock);
    if (arg->localPermute) {
      VectorArray<T> xtmp(arg->upper->nrowBlock(k));
      for (int j = 0; j < ncol; j++) {
	int i, ii;
	int jnrow = j * nrow;
	for (ii = kk, i = 0; i < nrow; i++, ii++) {
	  const int ip = permute[ii] - kk;  
	  xtmp[i] = upper[ip + jnrow];  
	}
	blas_copy<T>(nrow, xtmp.addrCoefs(), 1, lower + jnrow, 1);
      } // loop : j
      xtmp.free();
    }
    if (n0 > 0) {
	for (vector<int>::const_iterator it = singLocNodes0->begin(); 
	   it != singLocNodes0->end();
	   ++it) {
	if (((*it) >= kk) && ((*it) < (kk + nrow))) {
	  nn0++;
	  int itmp = (*it) - kk;
	  for (int j = 0; j < ncol; j++, itmp += nrow) {
	    lower[itmp] = zero;
	  }
	}
      } // loop : it
    }
    blas_trsm<T>(CblasLeft, CblasLower, CblasNoTrans, CblasUnit,
		 (nrow - nn0),  // skip computation nullified entries
		 ncol,
		 one, //alpha, 
		 LDLt,
		 nrow, //
		 lower, nrow); //lower + kk,  n
    for (int j = 0; j < ncol; j++) {   
      const int jnrow = j * nrow;
      for (int i = 0; i < nrow; i++) {
	upper[i + jnrow] = lower[i + jnrow] * LDLt[i * (nrow + 1)];
      }
    }
  }  //   if (arg->isSym)
  else {
    //    lower = arg->lower->addrCoefBlock(k, arg->lblock);
    if (arg->localPermute) {
      VectorArray<T> xtmp(arg->upper->nrowBlock(k));
      for (int j = 0; j < ncol; j++) {
	int i, ii;
	int jnrow = j * nrow;
	for (ii = kk, i = 0; i < nrow; i++, ii++) {
	  const int ip = permute[ii] - kk;  
	  xtmp[i] = upper[ip + jnrow];  
	}
	blas_copy<T>(nrow, xtmp.addrCoefs(), 1, upper + jnrow, 1);
      } // loop : j
      xtmp.free();
    }
    if (n0 > 0) {
      for (vector<int>::const_iterator it = singLocNodes0->begin(); 
	   it != singLocNodes0->end();
	   ++it) {
	if (((*it) >= kk) && ((*it) < (kk + nrow))) {
	  int itmp = (*it) - kk;
	  nn0++;
	  for (int j = 0; j < ncol; j++, itmp += nrow) {
	    upper[itmp] = zero;
	  }
	}
      } // loop : it
    }
    blas_trsm<T>(CblasLeft, CblasLower, CblasNoTrans, CblasUnit,
		 (nrow - nn0),  // skip computation nullified entries
		 ncol,
		 one, //alpha, 
		 LDLt, nrow,
		 upper, nrow);  // upper + kk, n);
    for (int j = 0; j < ncol; j++) {   
      const int jnrow = j * nrow;
      for (int i = 0; i < nrow; i++) {
	upper[i + jnrow] *= LDLt[i * (nrow + 1)];
      }
    }
  }  // else if (arg->isSym)
}

template
void C_DTRSMScale_diag_upper<double>(void *arg_);

template
void C_DTRSMScale_diag_upper<quadruple>(void *arg_);

template
void C_DTRSMScale_diag_upper<complex<double> >(void *arg_);

template
void C_DTRSMScale_diag_upper<complex<quadruple> >(void *arg_);

template
void C_DTRSMScale_diag_upper<float>(void *arg_);

template
void C_DTRSMScale_diag_upper<complex<float> >(void *arg_);

//

template<typename T>
void C_DTRSMScale_offdiag_upper(void *arg_)
{
  const T none(-1.0);
  const T one(1.0);
  
  DTRSMScale_arg<T> *arg = (DTRSMScale_arg<T> *)arg_;
  int ncol = arg->ncol;
  const int k = arg->kblock;
  const int m = arg->mblock;
  vector<int>* singLocNodes0 = arg->singLocNodes0;
  const int n0 = singLocNodes0->size();
  const int kk = arg->LDLt->IndexBlock(k);
  const int nrow = arg->LDLt->nrowBlock(k, k);
  const int nrow1 = arg->LDLt->ncolBlock(k, m); //

  int nn0 = 0;
  if (n0 > 0) {
    for (vector<int>::const_iterator it = singLocNodes0->begin(); 
	 it != singLocNodes0->end();
	 ++it) {
      if (((*it) >= kk) && ((*it) < (kk + nrow))) {
	nn0++;
      }
    } // loop : it
  }
  if (arg->isSym) {
    blas_gemm<T>(CblasTrans, CblasNoTrans,
		 nrow1, ncol, // based on row/column sizes of "upper"
		 (nrow - nn0), // skip computation nullified entries
		 none, // alpha, 
		 arg->LDLt->addrCoefBlock(k, m), // transposed upper
		 arg->LDLt->nrowBlock(k, m),
		 arg->lower->addrCoefBlock(k, arg->lblock),
		 arg->lower->nrowBlock(k),
		 one, // beta, 
		 arg->upper->addrCoefBlock(m, arg->lblock),
		 arg->upper->nrowBlock(m));
  }
  else {
    blas_gemm<T>(CblasTrans, CblasNoTrans,
		 nrow1, ncol, // based on row/column sizes of "upper"
		 (nrow - nn0), // skip computation nullified entries
		 none, // alpha, 
		 arg->LDLt->addrCoefBlock(m, k), // transposed lower
		 arg->LDLt->nrowBlock(m, k),
		 arg->upper->addrCoefBlock(k, arg->lblock),
		 arg->upper->nrowBlock(k), 
		 one, // beta, 
		 arg->upper->addrCoefBlock(m, arg->lblock),
		 arg->upper->nrowBlock(m));
  }
}

template
void C_DTRSMScale_offdiag_upper<double>(void *arg_);

template
void C_DTRSMScale_offdiag_upper<quadruple>(void *arg_);

template
void C_DTRSMScale_offdiag_upper<complex<double> >(void *arg_);

template
void C_DTRSMScale_offdiag_upper<complex<quadruple> >(void *arg_);

template
void C_DTRSMScale_offdiag_upper<float>(void *arg_);

template
void C_DTRSMScale_offdiag_upper<complex<float> >(void *arg_);
//

template<typename T>
void C_DTRSMScale_diag_lower(void *arg_)
{
  const T one(1.0);
  const T zero(0.0);
  DTRSMScale_arg<T> *arg = (DTRSMScale_arg<T> *)arg_;
  int ncol = arg->ncol;
  vector<int> &permute = arg->LDLt->getPermute();
  const int k = arg->kblock;
  T *lower = arg->lower->addrCoefBlock(k, arg->lblock);
  const int nrow = arg->LDLt->nrowBlock(k, k);
  const int kk = arg->LDLt->IndexBlock(k);
  vector<int>* singLocNodes0 = arg->singLocNodes0;
  const int n0 = singLocNodes0->size();
      // direct use of lower block matrix
  if (arg->localPermute) {
    VectorArray<T> xtmp(arg->lower->nrowBlock(k));   // SIZE_B1];
    for (int j = 0; j < ncol; j++) {
      int i, ii;
      int jnrow = j * nrow;
      for (ii = kk, i = 0; i < nrow; i++, ii++) {
	const int ip = permute[ii] - kk;  
	xtmp[i] = lower[ip + jnrow];  
      }
      blas_copy<T>(nrow, xtmp.addrCoefs(), 1, lower + jnrow, 1);
    } // loop : j
    xtmp.free();
  }
  // nullifying rows corresponding to singular nodes
  int nn0 = 0;
  if (n0 > 0) {
    for (vector<int>::const_iterator it = singLocNodes0->begin(); 
	 it != singLocNodes0->end();
	 ++it) {
      if (((*it) >= kk) && ((*it) < (kk + nrow))) {
	nn0++;
	int itmp = (*it) - kk;
	for (int j = 0; j < ncol; j++, itmp += nrow) {
	  lower[itmp] = zero;
	}
      }
    } // loop : it
  } // if (n0 > 0)
  blas_trsm<T>(CblasLeft, CblasUpper, CblasTrans, CblasUnit,
	       (nrow - nn0),   // skip computation nullified entries
	       ncol,
	       one, // alpha, 
	       arg->LDLt->addrCoefBlock(k, k), 
	       arg->LDLt->nrowBlock(k, k),     
	       lower,
	       nrow);
}

template
void C_DTRSMScale_diag_lower<double>(void *arg_);

template
void C_DTRSMScale_diag_lower<complex<double> >(void *arg_);

template
void C_DTRSMScale_diag_lower<quadruple>(void *arg_);

template
void C_DTRSMScale_diag_lower<complex<quadruple> >(void *arg_);

template
void C_DTRSMScale_diag_lower<float>(void *arg_);

template
void C_DTRSMScale_diag_lower<complex<float> >(void *arg_);
//

template<typename T>
void C_DTRSMScale_offdiag_lower(void *arg_)
{
  DTRSMScale_arg<T> *arg = (DTRSMScale_arg<T> *)arg_;
  int ncol = arg->ncol;
  const int k = arg->kblock;
  const int m = arg->mblock;
  const int nrow = arg->LDLt->nrowBlock(k, k);
  const int kk = arg->LDLt->IndexBlock(k);
  const int nrow1 = arg->LDLt->ncolBlock(k, m);
  vector<int>* singLocNodes0 = arg->singLocNodes0;
  const int n0 = singLocNodes0->size();
  int nn0 = 0;

  if (n0 > 0) {
    for (vector<int>::const_iterator it = singLocNodes0->begin(); 
	 it != singLocNodes0->end();
	 ++it) {
      if (((*it) >= kk) && ((*it) < (kk + nrow))) {
	nn0++;
      }
    } // loop : it
  }
  // x -= A_21 P_1^T L_11^-1 D_11^-1 zz
  const T none(-1.0);
  const T one(1.0);
  blas_gemm<T>(CblasTrans, CblasNoTrans,
	       nrow1, ncol, 
	       (nrow - nn0),  // skip computation nullified entries
	       none, // alpha, 
	       arg->LDLt->addrCoefBlock(k, m), // transposed upper
	       arg->LDLt->nrowBlock(k, m), 
	       arg->lower->addrCoefBlock(k, arg->lblock), 
	       arg->lower->nrowBlock(k),
	       one, // beta, 
	       arg->lower->addrCoefBlock(m, arg->lblock), 
	       arg->lower->nrowBlock(m));
}
template
void C_DTRSMScale_offdiag_lower<double>(void *arg_);

template
void C_DTRSMScale_offdiag_lower<quadruple>(void *arg_);

template
void C_DTRSMScale_offdiag_lower<complex<double> >(void *arg_);

template
void C_DTRSMScale_offdiag_lower<complex<quadruple> >(void *arg_);

template
void C_DTRSMScale_offdiag_lower<float>(void *arg_);

template
void C_DTRSMScale_offdiag_lower<complex<float> >(void *arg_);
//

template<typename T>
void C_DTRSMScale_solve(void *arg_)
{
  DTRSMScale_arg<T> *arg = (DTRSMScale_arg<T> *)arg_;
                         // should be a local array belonging to the same CPU
  bool verbose = arg->verbose;
  FILE *fp = *(arg->fp);
  const int num_block = arg->LDLt->num_blocks();

   if (!arg->localPermute) {
    VectorArray<T> xtmp(arg->nrow);
    vector<int> &permute = arg->LDLt->getPermute();
    diss_printf(verbose, fp, "%s %d : non blocked forward substituion\n",
	      __FILE__, __LINE__);
    //    const int num_block = arg->LDLt->num_blocks(); 
    if (arg->isSym) {
      for (int k = 0; k < num_block; k++) {
	arg->lower->allocateBlock(k, arg->lblock);
      }
    }
    for (int j = 0; j < arg->ncol; j++) {
      for (int i = 0; i < arg->nrow; i++) {
	xtmp[i] = (*arg->upper)(permute[i], j);
      }
      // blas_copy<T>(arg->nrow, xtmp, 1, arg->upper, 1); // block copy
      // arg->nrow need to be composed into blocks  : 22 Feb.2016
      for (int i = 0; i < arg->nrow; i++) {
	(*arg->upper)(i, j) = xtmp[i];
      }
    }
    if (arg->isSym) {
      for (int j = 0; j < arg->ncol; j++) {
	for (int i = 0; i < arg->nrow; i++) {
	  (*arg->lower)(i, j) = (*arg->upper)(i, j);
	}
      }      
    }
    else {
      for (int j = 0; j < arg->ncol; j++) {
	for (int i = 0; i < arg->nrow; i++) {
	  xtmp[i] = (*arg->lower)(permute[i], j);
	}
	for (int i = 0; i < arg->nrow; i++) {
	  (*arg->lower)(i, j) = xtmp[i];
	}
      }
    }
    xtmp.free();
  } // if (!arg->localPermute) 
  //    const int num_block = arg->LDLt->num_blocks();
  for (int l = 0; l < (*arg->upper).num_blocks_c(); l++) {
    arg->lblock = l;
    arg->ncol = (*arg->upper).ncolBlock(l);
    for (int k = 0; k < num_block; k++) {
      arg->kblock = k;
      C_DTRSMScale_diag_upper<T>(arg_);
      for (int m = (k + 1); m < num_block; m++) {
	arg->mblock = m;
	C_DTRSMScale_offdiag_upper<T>(arg_);
      }
      if (!arg->isSym) {
	C_DTRSMScale_diag_lower<T>(arg_);
	for (int m = (k + 1); m < num_block; m++) {
	  arg->mblock = m;
	  C_DTRSMScale_offdiag_lower<T>(arg_);
	}
      } // if (!arg->isSym)
    } // loop : k
  }
}

template
void C_DTRSMScale_solve<double>(void *arg_);

template
void C_DTRSMScale_solve<quadruple>(void *arg_);

template
void C_DTRSMScale_solve<complex<double> >(void *arg_);

template
void C_DTRSMScale_solve<complex<quadruple> >(void *arg_);

template
void C_DTRSMScale_solve<float>(void *arg_);

template
void C_DTRSMScale_solve<complex<float> >(void *arg_);

//

template<typename T>
void C_deallocLower(void *arg_)
{
  C_deallocLower_arg<T> *arg = (C_deallocLower_arg<T> *)arg_;
  if (arg->isSym) {
    arg->lower->free();  
  }
}
template
void C_deallocLower<double>(void *arg_);

template
void C_deallocLower<quadruple>(void *arg_);

template
void C_deallocLower<complex<double> >(void *arg_);

template
void C_deallocLower<complex<quadruple> >(void *arg_);

template
void C_deallocLower<float>(void *arg_);

template
void C_deallocLower<complex<float> >(void *arg_);
//

template<typename T>
void C_deallocLocalSchur(void *arg_)
{
  C_deallocLocalSchur_arg<T> *arg = (C_deallocLocalSchur_arg<T> *)arg_;
  arg->localSchur->free(arg->i_block, arg->j_block); 
}
template
void C_deallocLocalSchur<double>(void *arg_);

template
void C_deallocLocalSchur<quadruple>(void *arg_);

template
void C_deallocLocalSchur<complex<double> >(void *arg_);

template
void C_deallocLocalSchur<complex<quadruple> >(void *arg_);

template
void C_deallocLocalSchur<float>(void *arg_);

template
void C_deallocLocalSchur<complex<float> >(void *arg_);

//

template<typename T, typename U>
void C_SparseFw(void *arg_)
{
  const T zero(0.0);
  C_SparseFw_arg<T, U> *arg = (C_SparseFw_arg<T, U> *)arg_;
  int dim = arg->dim;
  const bool isTrans = **(arg->isTrans);
  const int nrhs = **(arg->nrhs);
  //  const void* diag_sparse = arg->diag_sparse;
  int colors = arg->colors;
  TridiagBlockMatrix<T, U> **tridiag = arg->tridiag;
  T *x = *(arg->x);
  T *yi = *(arg->yi);
  T *zi = *(arg->zi);
  T *coef = arg->coef; // arg->ptDA->getCoef();
  const int n_diag = arg->n_diag;
  const int n_offdiag = arg->n_offdiag;
  const int *ptRows = arg->ptRows;
  const int *indCols = arg->indCols;
  // arg->inTrans is selected as runing time
  const int *indVals =  **(arg->isTrans) ? arg->indVals_unsym : arg->indVals;
  const int *loc2glob_diag = arg->loc2glob_diag;
  
  // get from the global array

  for (int m = 0; m < nrhs; m++) {
    const int mn_diag = m * n_diag;
    const int mdim = m * dim;
    for (int i = 0; i < n_diag; i++) {
      const int ii = loc2glob_diag[i]; 
      yi[i + mn_diag] = x[ii + mdim];
    }
  } // loop : m
  if (nrhs == 1) {
    for (int i = 0; i < colors; i++) {
      tridiag[i]->SolveSingle(true, isTrans, yi);
    }
  }
  else {
    ColumnMatrix<T> rhs(n_diag, nrhs, yi, false);
    for (int i = 0; i < colors; i++) {
      tridiag[i]->SolveMulti(true, isTrans, nrhs, rhs);
    }
  }
  // zero clean
  for (int i = 0; i < n_offdiag * nrhs; i++) {
    zi[i] = zero;
  }
  // transposed Block SpMV : unsymmetric indVals is given as indVals_unsym
  for (int i = 0; i < n_diag; i++) {
    for (int k = ptRows[i]; k < ptRows[i + 1]; k++) {
      const int j = indCols[k];
      T tmp = coef[indVals[k]];
      for (int m = 0; m < nrhs; m++) {
	const int mn_diag = m * n_diag;
	const int mn_offdiag = m * n_offdiag;
	zi[j + mn_offdiag] += tmp * yi[i + mn_diag];
      }
    }
  }
}
template
void C_SparseFw<double, double>(void *arg_);

template
void C_SparseFw<quadruple, quadruple>(void *arg_);

template
void C_SparseFw<complex<double>, double>(void *arg_);

template
void C_SparseFw<complex<quadruple>, quadruple>(void *arg_);

template
void C_SparseFw<float, float>(void *arg_);

template
void C_SparseFw<complex<float>, float>(void *arg_);
//

template<typename T, typename U>
void C_SparseBw(void *arg_)
{
  const T zero(0.0);
  C_SparseBw_arg<T, U> *arg = (C_SparseBw_arg<T, U> *)arg_;
  int d = arg->nb;
  int dim = arg->dim;
  const bool isTrans = **(arg->isTrans);
  const int nrhs = **(arg->nrhs);
  Dissection::Tree *btree = arg->btree;
  const int level_last = arg->level_last;
  //  const void* diag_sparse = arg->diag_sparse;
  int colors = arg->colors;
  TridiagBlockMatrix<T, U> **tridiag = arg->tridiag;
  T *x = *(arg->x);
  T *xi = *(arg->xi);
  T ***yy = arg->yy;
  T *yi = *(arg->yi);
  T *zi = *(arg->zi);
  T *coef = arg->coef; // arg->ptDA->getCoef();
  const int *ptRows = arg->ptRows;
  const int *indCols = arg->indCols;
  const int *indVals = isTrans ? arg->indVals_unsym : arg->indVals;

  const int n_diag = btree->sizeOfDomain(d);
  const int n_offdiag = btree->sizeOfFathersStrips(d); 

  for (int i = 0; i < (n_diag * nrhs); i++) {
    xi[i] = zero;
  }
// reading data does not cause conflict between processors
  int offset_src = 0;
  for (int ll = (level_last - 1); ll >= 0; ll--) {
    const int father_id = btree->nthfatherIndex(d, (level_last - ll));
    const int father_id0 = btree->selfIndex(father_id);
    const Dissection::SetOfStrips &diag = btree->getFathersStrips(d)[ll];
    for (Dissection::SetOfStrips::const_iterator it = diag.begin();
	 it != diag.end(); ++it) {
      for (int m = 0; m < nrhs; m++) {
	const int mn_offdiag = m * n_offdiag;
	const int mn_diag = m * btree->sizeOfDomain(father_id);
	int i, i0, i1;
	for (i = 0, i0 = (*it).begin_src + offset_src, i1 = (*it).begin_dst;
	     i < (*it).width; i++, i0++, i1++) {
	  zi[i0 + mn_offdiag] = (*yy[father_id0])[i1 + mn_diag];
	} // loop : i
      }   // loop : m
    }     // loop : it
    offset_src += diag.numberOfIndices();
  }       // loop : ll
  // SpMV operations : y_d = A_ds y_s
  for (int i = 0; i < n_diag; i++) {
    for (int k = ptRows[i]; k < ptRows[i + 1]; k++) {
      const int jj = indCols[k];
      const T tmp = coef[indVals[k]];
      for (int m = 0; m < nrhs; m++) {
	const int mn_diag = m * n_diag;
	const int mn_offdiag = m * n_offdiag;
	xi[i + mn_diag] += tmp * zi[jj + mn_offdiag];
      }
    }
  }
  //      const int idom = _btree->selfIndex(d);
  if (nrhs == 1) {
    for (int i = 0; i < colors; i++) {
      tridiag[i]->SolveSingle(true, isTrans, xi);
    }
  }
  else {
    ColumnMatrix<T> rhs(n_diag, nrhs, xi, false);
    for (int i = 0; i < colors; i++) {
      tridiag[i]->SolveMulti(true, isTrans, nrhs, rhs);
    }
  }
  for (int i = 0; i < (n_diag * nrhs); i++) {
    yi[i] -= xi[i]; // without permutation
  }
  // write back to the global array
  const int *loc2glob_diag = btree->getDiagLoc2Glob(d);
  for (int m = 0; m < nrhs; m++) {
    const int mn_diag = m * n_diag;
    const int mdim = m * dim;
    for (int i = 0; i < n_diag; i++) {
      const int ii = loc2glob_diag[i];
      x[ii + mdim] = yi[i + mn_diag];
    }
  } // loop : m
}

template
void C_SparseBw<double, double>(void *arg_);

template
void C_SparseBw<quadruple, quadruple>(void *arg_);

template
void C_SparseBw<complex<double>, double>(void *arg_);

template
void C_SparseBw<complex<quadruple>, quadruple>(void *arg_);

template
void C_SparseBw<float, float>(void *arg_);

template
void C_SparseBw<complex<float>, float>(void *arg_);
//

template<typename T>
void C_Dsub_FwBw(void *arg_)
{
  C_Dsub_FwBw_arg<T> *arg = (C_Dsub_FwBw_arg<T> *)arg_;

  const int dim = arg->dim;
  const int nrhs = *(arg->nrhs)[0];
  const int level = arg->level;
  Dissection::Tree *btree = arg->btree;
  list<diag_contribution>* diag_contribs = arg->diag_contribs;
  T *yi = *(arg->yi);
  T ***zi = arg->zi;
  bool access_global = arg->access_global;
  if (access_global) {
    T *x = *(arg->x);
  // get from the global array
    const int n_diag = arg->n_diag;
    const int *loc2glob_diag = arg->loc2glob_diag;
    for (int m = 0; m < nrhs; m++) {
      const int mn_diag = m * n_diag;
      const int mdim = m * dim;
      for (int i = 0; i < n_diag; i++) {
	const int ii = loc2glob_diag[i]; 
	yi[i + mn_diag] = x[ii + mdim];
      }
    } // loop : m
  }
  for (list<diag_contribution>::const_iterator it = diag_contribs->begin();
       it != diag_contribs->end(); ++it) {
    const int child_id = (*it).child_id;
    const int child_id0 = btree->selfIndex(child_id);
    if ((btree->nodeLayer(child_id) == level) &&
	(btree->sizeOfDomain(child_id) > 0)) {   // 08 Nov.2016 Atsushi
      for (list<index_strip>::const_iterator mt = (*it).diag_strip.begin();
	   mt != (*it).diag_strip.end(); ++mt) {
	for (int m = 0; m < nrhs; m++) {
	  const int mn_offdiag = m * (*it).child_column;
	  const int mn_diag = m * (*it).father_row;
	  int i, i0, i1;
	  for (i = 0, i0 = (*mt).begin_src, i1 = (*mt).begin_dst;
	       i < (*mt).width; i++, i0++, i1++) {
	    yi[i1 + mn_diag] -= (*zi[child_id0])[i0 + mn_offdiag];
	  } // loop : i
	}   // loop : m
      }     // loop : mt
    }       // if (_btree->nodeLayer((*it).child_id) == level_last) {
  }       // loop : it 
} 

template
void C_Dsub_FwBw<double>(void *arg_);

template
void C_Dsub_FwBw<quadruple>(void *arg_);

template
void C_Dsub_FwBw<complex<double> >(void *arg_);

template
void C_Dsub_FwBw<complex<quadruple> >(void *arg_);

template
void C_Dsub_FwBw<float>(void *arg_);

template
void C_Dsub_FwBw<complex<float> >(void *arg_);
//

template<typename T>
void C_Dfill_FwBw(void *arg_)
{
  C_Dfill_FwBw_arg<T> *arg = (C_Dfill_FwBw_arg<T> *)arg_;

  const int nrhs = *(arg->nrhs)[0];
  const int d = arg->d;
  const int level = arg->level;
  Dissection::Tree *btree = arg->btree;
  const int n_offdiag = arg->n_offdiag;
  T ***yi = arg->yi;
  T *zi = *(arg->zi);

  int offset_src = 0;
  for (int ll = (level - 1); ll >= 0; ll--) {
    const int father_id = btree->nthfatherIndex(d, (level - ll));
    const int father_id0 = btree->selfIndex(father_id);
    const Dissection::SetOfStrips &diag = btree->getFathersStrips(d)[ll];
    for (Dissection::SetOfStrips::const_iterator it = diag.begin();
	 it != diag.end(); ++it) {
      for (int m = 0; m < nrhs; m++) {
	const int mn_offdiag = m * n_offdiag;
	const int mn_diag = m * btree->sizeOfDomain(father_id);
	int i, i0, i1;
	for (i = 0, i0 = (*it).begin_src + offset_src, 
	       i1 = (*it).begin_dst;
	     i < (*it).width; i++, i0++, i1++) {
	  zi[i0 + mn_offdiag] = (*yi[father_id0])[i1 + mn_diag];
	} // loop : i
      }   // loop : m
    }     // loop : it
    offset_src += diag.numberOfIndices();
  }       // loop : ll

}
template
void C_Dfill_FwBw<double>(void *arg_);

template
void C_Dfill_FwBw<quadruple>(void *arg_);

template
void C_Dfill_FwBw<complex<double> >(void *arg_);

template
void C_Dfill_FwBw<complex<quadruple> >(void *arg_);

template
void C_Dfill_FwBw<float>(void *arg_);

template
void C_Dfill_FwBw<complex<float> >(void *arg_);
//

template<typename T>
void C_DenseFwBw_diag(void *arg_)
{
  const T one(1.0);
  const T zero(0.0);
  C_DenseFwBw_arg<T> *arg = (C_DenseFwBw_arg<T> *)arg_;
  const int dim = arg->dim;   // global dimension, used when nrhs > 1
  const bool isTrans = **(arg->isTrans);
  const int nrhs = **(arg->nrhs);
  const int n_diag = arg->n_diag;
  const int nrow = arg->nrow;
  const int k_block = arg->k_block;
  T *xi = *(arg->xi);
  const bool isSym = arg->isSym;
  const bool isBackward = arg->isBackward;
  const bool isFirstBlock = arg->isFirstBlock;
  SquareBlockMatrix<T> &Diag = *(arg->LDLt);
  const bool isBlocked = Diag.isBlocked();
  T *LDLt = Diag.addrCoefBlock(k_block, k_block);
  const int LDLt_nrow = Diag.nrowBlock(k_block, k_block);
  const int kk = Diag.IndexBlock(k_block);
  vector<int> singIdx0 = Diag.getSingIdx0();
  const bool verbose = arg->verbose;
  FILE *fp = *(arg->fp);
  if (isBlocked) {
    if (!isBackward) {
      T *yi = *(arg->yi);
      vector<int> &permute = Diag.getPermute();
      // permutation : original index -> factorized index
      for (int m = 0; m < nrhs; m++) {
	const int mn_diag = m * n_diag;
	for (int i = kk; i < kk + nrow; i++) {
	  xi[i + mn_diag] = yi[permute[i] + mn_diag];
	}
      }
    }
  }
  else { //   if (isBlocked)
    // permute the whole RHS by the k == 0 block
    if (isFirstBlock && (isBackward == false)) {
      T *yi = *(arg->yi);
      vector<int> &permute = Diag.getPermute();
      // permutation : original index -> factorized index
      for (int m = 0; m < nrhs; m++) {
	const int mn_diag = m * n_diag;
	for (int i = 0; i < n_diag; i++) {
	  xi[i + mn_diag] = yi[permute[i] + mn_diag];
	}
      }
    }
  }      //  if (isBlocked) 
  if (nrhs == 1) {
    int nn0 = 0;
    if (singIdx0.size() > 0) {
      for (vector<int>::const_iterator it = singIdx0.begin(); 
	   it != singIdx0.end(); ++it) {
	// inside of the block
	if (((*it) >= kk) && ((*it) < (kk + nrow))) {
	  xi[(*it)] = zero;
	  nn0++;
	}
      }
      // if (d0 == 0) degmv for the regular part of singIdx0
    }  // if (singIdx.size() > 0)
    const int nn1 = nrow - nn0; // invertible part
    if (isSym) {
      blas_trsv<T>(CblasLower, (isBackward ? CblasTrans: CblasNoTrans), 
		   CblasUnit,
		   nn1, 
		   LDLt, LDLt_nrow, //n_diag, 
		   xi + kk, 1);
    }
    else {
      if (isTrans) {
	if (isBackward) { //Lower bocks of matrix are not scaled then scale here
	  int itmp = 0;
	  for (int i = kk; i < (kk + nrow); i++, itmp += (LDLt_nrow + 1)) {
	    xi[i] *= LDLt[itmp];
	  }
	}
	blas_trsv<T>((isBackward ? CblasLower : CblasUpper), CblasTrans, 
		     CblasUnit,
		     nn1, 
		     LDLt, LDLt_nrow, //n_diag,  
		     xi + kk, 1);
      }
      else {
	blas_trsv<T>((isBackward ? CblasUpper : CblasLower), CblasNoTrans, 
		     CblasUnit,
		     nn1, 
		     LDLt, LDLt_nrow, //n_diag,  
		     xi + kk, 1);
      }
    }
      // if (d0 == 0) degmv for the regular part of singIdx0
      // if (singIdx.size() > 0) 
    if (!isBackward) {
      T *wi = *(arg->wi);
      if (isSym) {
	int itmp = 0;
	for (int i = kk; i < (kk + nrow); i++, itmp += (LDLt_nrow + 1)) {
	  wi[i] = xi[i];
	  xi[i] *= LDLt[itmp];
	}
      }
      else {
	if(isTrans) {
	  for (int i = kk; i < (kk + nrow); i++) {
	    wi[i] = xi[i];      // Upper bocks of matrix are scaled
	  }
	}
	else {
	  int itmp = 0;
	  for (int i = kk; i < (kk + nrow); i++, itmp += (LDLt_nrow + 1)) {
	    xi[i] *= LDLt[itmp]; // Lower bocks of matrix are not scaled
	    wi[i] = xi[i];
	  }
	}
      }
    }
  }  
  else { //   if (nrhs != 1) 
    int nn0 = 0;
    if (singIdx0.size() > 0) {  
      for (vector<int>::const_iterator it = singIdx0.begin(); 
	   it != singIdx0.end(); ++it) {
	// inside of th block
	if (((*it) >= kk) && ((*it) < (kk + nrow))) {
          nn0++;
	  int itmp = (*it); 
	  for (int m = 0; m < nrhs; m++, itmp += n_diag) {
	    xi[itmp] = zero;
	  }
	}
      }
    } // if (singIdx.size() > 0)
    const int nn1 = nrow - nn0; // invertible part
    if (isSym) {
      blas_trsm<T>(CblasLeft, 
		   CblasLower, (isBackward ? CblasTrans : CblasNoTrans), 
		   CblasUnit,
		   nn1, 
		   nrhs,
		   one, //alpha, 
		   LDLt, LDLt_nrow, //n_diag,  
		   xi + kk, n_diag);
    }
    else {
      if (isTrans) {
	if (isBackward) {
	  int itmp = 0;
	  for (int i = kk; i < (kk + nrow); i++, itmp += (LDLt_nrow + 1)) {
	    for (int m = 0; m < nrhs; m++) {
	      const int mn_diag = m * n_diag;
	      xi[i + mn_diag] *= LDLt[itmp];
	    }
	  }
	}
	blas_trsm<T>(CblasLeft, 
		     (isBackward ? CblasLower : CblasUpper), CblasTrans, 
		     CblasUnit,
		     nn1, 
		     nrhs,
		     one, // alpha, 
		     LDLt, LDLt_nrow, //n_diag,  
		     xi + kk, n_diag);
      }
      else {
	blas_trsm<T>(CblasLeft, 
		     (isBackward ? CblasUpper : CblasLower), CblasNoTrans, 
		     CblasUnit,
		     nn1, 
		     nrhs,
		     one, //alpha, 
		     LDLt, LDLt_nrow, //n_diag,  
		     xi + kk, n_diag);
      }
    }
    if (!isBackward) {
      T *wi = *(arg->wi);
      if(isSym) {
	int itmp = 0;
	for (int i = kk; i < (kk + nrow); i++, itmp += (LDLt_nrow + 1)) {
	  for (int m = 0; m < nrhs; m++) {
	    const int mn_diag = m * n_diag;
	    wi[i + mn_diag] = xi[i + mn_diag];  // for yl[] - A_li xi[] : l>i
	    xi[i + mn_diag] *= LDLt[itmp];
	  }
	}
      }
      else {
	if (isTrans) {
	  for (int i = kk; i < (kk + nrow); i++) {
	    for (int m = 0; m < nrhs; m++) {
	      const int mn_diag = m * n_diag;
	      wi[i + mn_diag] = xi[i + mn_diag];  // for yl[] - A_li xi[] : l>i
	    }
	  }
	}
	else {
	  int itmp = 0;
	  for (int i = kk; i < (kk + nrow); i++, itmp += (LDLt_nrow + 1)) {
	    for (int m = 0; m < nrhs; m++) {
	      const int mn_diag = m * n_diag;
	      xi[i + mn_diag] *= LDLt[itmp];
	      wi[i + mn_diag] = xi[i + mn_diag];  // for yl[] - A_li xi[] : l>i
	    }
	  }
	}
      }
    } // if (!isBackward)
  } //   if (nrhs == 1) 
  if (isBlocked) {
    if (isBackward) {
      T *yi = *(arg->yi);
      vector<int> &permute = Diag.getPermute();
      for (int m = 0; m < nrhs; m++) {
	const int mn_diag = m * n_diag;
	for (int i = kk; i < kk + nrow; i++) {
	  yi[permute[i] + mn_diag] = xi[i + mn_diag];
	}
	for (int i = kk; i < kk + nrow; i++) {
	  xi[i + mn_diag] = yi[i + mn_diag];
	}
      }
      if(isFirstBlock) {
	int *loc2glob_diag = arg->loc2glob;
	T *x = *(arg->x);
	for (int m = 0; m < nrhs; m++) {
	  const int mn_diag = m * n_diag;
	  const int mdim = m * dim;
	  for (int i = 0; i < n_diag; i++) {
	    const int ii = loc2glob_diag[i];
	    x[ii + mdim] = xi[i + mn_diag];
	  }
	} // loop : m
      } // if(isFirstBlock) 
    } // if (isBackward)
    // forward substitution with block permutation : internally no permutation
  }
  else {
    if (isFirstBlock && isBackward) {
      T *yi = *(arg->yi);
      vector<int> &permute = Diag.getPermute();
      for (int m = 0; m < nrhs; m++) {
	const int mn_diag = m * n_diag;
	for (int i = 0; i < n_diag; i++) {
	  yi[permute[i] + mn_diag] = xi[i + mn_diag];
	}
      }
      int *loc2glob_diag = arg->loc2glob;
      T *x = *(arg->x);
      for (int m = 0; m < nrhs; m++) {
	const int mn_diag = m * n_diag;
	const int mdim = m * dim;
	for (int i = 0; i < n_diag; i++) {
	  const int ii = loc2glob_diag[i];
	  x[ii + mdim] = yi[i + mn_diag];
	}
      } // loop : m
    }   //  if (isFirstBlock && isBackward) {
  }     //  if (isBlocked)
}
template
void C_DenseFwBw_diag<double>(void *arg_);

template
void C_DenseFwBw_diag<quadruple>(void *arg_);

template
void C_DenseFwBw_diag<complex<double> >(void *arg_);

template
void C_DenseFwBw_diag<complex<quadruple> >(void *arg_);

template
void C_DenseFwBw_diag<float>(void *arg_);

template
void C_DenseFwBw_diag<complex<float> >(void *arg_);

//

template<typename T>
void C_DenseFwBw_offdiag(void *arg_)
{
  C_DenseFwBwOffdiag_arg<T> *arg = (C_DenseFwBwOffdiag_arg<T> *)arg_;
  const bool isTrans = **(arg->isTrans);
  const int nrhs = **(arg->nrhs);
  const int ldb = arg->ldb;
  const int ldc = arg->ldc;
  const int nrow = arg->nrow;
  const int ncol = arg->ncol;
  T *xi = *(arg->xi);
  T *yi;
  const int ii = arg->ii;
  const int jj = arg->jj;
  T *a; 
  const T alpha = arg->alpha;
  const T beta = arg->beta;
  const bool trans = arg->trans;
  const bool isLower = arg->isLower;
  int i_block, j_block;
  int lda;
  // #define DEBUG_INDEX
#ifdef DEBUG_INDEX
  fprintf(stderr, "%s %d : nrow = %d ncol = %d %s [%d %d] [%d %d] %d : %s %s ", 
	  __FILE__, __LINE__, nrow, ncol,
	  isLower ? "Lower" : "Upper", arg->i_block, arg->j_block, ii, jj, nrhs,
	  trans ? "T" : "N", isTrans ? "trans" : "normal");
#endif
  {
    if (arg->LDLt->isBlocked()) {
      yi = *(arg->yi);
    }
    else {
      yi = *(arg->zi);
    }
    if (trans) {  // noly for _isSym == true
      i_block = arg->j_block;
      j_block = arg->i_block;
    }
    else {
      if (isTrans) {
	i_block = arg->j_block;
	j_block = arg->i_block;
      }
      else {
	i_block = arg->i_block;
	j_block = arg->j_block;
      }
    }
    a = arg->LDLt->addrCoefBlock(i_block, j_block);
    lda = arg->LDLt->nrowBlock(i_block, j_block);
#ifdef DEBUG_INDEX
	fprintf(stderr, "%d x %d\n",
		arg->LDLt->nrowBlock(arg->j_block, arg->i_block),
		arg->LDLt->ncolBlock(arg->j_block, arg->i_block));
#endif
    //
    if (isLower) {
      if (nrhs == 1) {
	// nrow and ncol are based on geometrical size of the block
	// data storage is transposed in the case of the lower block
	blas_gemv<T>(CblasTrans,
		     ncol, nrow, // based on a
		     alpha, 
		     a, lda,
		     xi + ii, 1, 
		     beta, 
		     yi + jj, 1);
      }
      else {
	blas_gemm<T>(CblasTrans,
		     CblasNoTrans,
		     nrow, nrhs, ncol,  // based on yi[]
		     alpha, 
		     a, lda,
		     xi + ii, ldb,
		     beta, 
		     yi + jj, ldc);
      }
    }
    else {
      if (nrhs == 1) {
	blas_gemv<T>(CblasNoTrans,
		     nrow, ncol, // based on a
		     alpha, 
		     a, lda,
		     xi + ii, 1, 
		     beta, 
		     yi + jj, 1);
      }
      else {
	blas_gemm<T>(CblasNoTrans,
		     CblasNoTrans,
		     nrow, nrhs, ncol, // based on yi[]
		     alpha, 
		     a, lda,
		     xi + ii, ldb, 
		     beta, 
		     yi + jj, ldc);
      }
    }
  } // if (isLocalDiag)
}
template
void C_DenseFwBw_offdiag<double>(void *arg_);

template
void C_DenseFwBw_offdiag<quadruple>(void *arg_);

template
void C_DenseFwBw_offdiag<complex<double> >(void *arg_);

template
void C_DenseFwBw_offdiag<complex<quadruple> >(void *arg_);

template
void C_DenseFwBw_offdiag<float>(void *arg_);

template
void C_DenseFwBw_offdiag<complex<float> >(void *arg_);
//

template<typename T>
void C_StripsFwBw_offdiag(void *arg_)
{
  C_StripsFwBwOffdiag_arg<T> *arg = (C_StripsFwBwOffdiag_arg<T> *)arg_;
  const bool isTrans = **(arg->isTrans);
  const int nrhs = **(arg->nrhs);
  const int ldb = arg->ldb;
  const int ldc = arg->ldc;
  const int ii = arg->ii;
  const int jj = arg->jj;
  const T alpha = arg->alpha;
  const T beta = arg->beta;
  //  const bool trans = arg->trans;
  const bool isLower = arg->isLower;

  T *xi = *(arg->xi);
  T *yi = *(arg->yi);

  RectBlockMatrix<T> *aa = isTrans ? arg->lower : arg->upper;
  const int i_block = isLower ? arg->i_block : arg->j_block;
  const int j_block = isLower ? arg->j_block : arg->i_block;
  
  if (isLower) {
    if (nrhs == 1) {
	// nrow and ncol are based on geometrical size of the block
	// data storage is transposed in the case of the lower block
      blas_gemv<T>(CblasTrans,
		   aa->nrowBlock(i_block), // based on matrix A
		   aa->ncolBlock(j_block), // 
		   alpha,
		   aa->addrCoefBlock(i_block, j_block), 
		   aa->nrowBlock(i_block),	  //  lda,
		   xi + ii, 1, 
		   beta, 
		   yi + jj, 1);
    }
    else {
      blas_gemm<T>(CblasTrans,
		   CblasNoTrans,
		   aa->ncolBlock(j_block), // based on yi[]
		   nrhs, 
		   aa->nrowBlock(i_block), // 
		   alpha, 
		   aa->addrCoefBlock(i_block, j_block), 
		   aa->nrowBlock(i_block), // lda,
		   xi + ii, ldb,
		   beta, 
		   yi + jj, ldc);
    }
  }
  else {
    const int nblocks_col = aa->num_blocks_c();
    if (nrhs == 1) {
      for (int j = 0; j < nblocks_col; j++) { // sequential operation.
	int jjj = aa->IndexBlock_c(j);
	blas_gemv<T>(CblasNoTrans,
		     // nrow, ncol, // based on a
		     aa->nrowBlock(i_block), // nrow, 
		     aa->ncolBlock(j), // ncol, 
		     alpha,
		     aa->addrCoefBlock(i_block, j),
		     aa->nrowBlock(i_block), // lda
		     xi + jjj, 1, 
		     beta, 
		     yi + jj, 1);
      }
    }
    else {
      for (int j = 0; j < nblocks_col; j++) { // sequential operation.
	int jjj = aa->IndexBlock_c(j);
	blas_gemm<T>(CblasNoTrans,
		     CblasNoTrans,
		     aa->nrowBlock(i_block), // nrow, 
		     nrhs, 
		     aa->ncolBlock(j), // ncol, // based on yi[]
		     alpha, 
		     aa->addrCoefBlock(i_block, j),
		     aa->nrowBlock(i_block), // lda,
		     xi + jjj, ldb, 
		     beta, 
		     yi + jj, ldc);
      }
    } // if (nrhs == 1)
  }   // if (isLower)
}  
template
void C_StripsFwBw_offdiag<double>(void *arg_);

template
void C_StripsFwBw_offdiag<quadruple>(void *arg_);

template
void C_StripsFwBw_offdiag<complex<double> >(void *arg_);

template
void C_StripsFwBw_offdiag<complex<quadruple> >(void *arg_);

template
void C_StripsFwBw_offdiag<float>(void *arg_);

template
void C_StripsFwBw_offdiag<complex<float> >(void *arg_);
//
template<typename T, typename U>
void erase_task(C_task *& task)
{
  bool flag = false;
  switch(task->task_id) {
  case C_DSUB:
    {
      list<C_Dsub_task<T>* > *tt = (list<C_Dsub_task<T>* > *)task->func_arg;
      for (typename list<C_Dsub_task<T>* >::iterator kt = tt->begin(); 
	   kt != tt->end(); ++kt) {
	delete (*kt);  // C_Dsub_task *tmp = new C_Dsub_task( ... )
	(*kt) = NULL;
      }
      flag = true;
      //  fprintf(stderr, "%s %d : %s\n", __FILE__, __LINE__, task->task_name);
      delete tt;
      tt = NULL;
    }
    break;
  case C_SPARSESYMBFACT:
    {
      C_SparseSymbFact_arg<T, U> *tt = (C_SparseSymbFact_arg<T, U> *)task->func_arg;
      delete tt;
      tt = NULL;
    }
    break;
  case C_SPARSENUMFACT:
    {
      C_SparseNumFact_arg<T, U> *tt = (C_SparseNumFact_arg<T, U> *)task->func_arg;
      delete tt;
      tt = NULL;
    }
    break;
  case C_SPARSESCHUR:
    {
      C_SparseNumFact_arg<T, U> *tt = (C_SparseNumFact_arg<T, U> *)task->func_arg;
      delete tt;
      tt = NULL;
    }
    break;
  case C_FILLMATRIX:
    {
      C_FillMatrix_arg<T> *tt = (C_FillMatrix_arg<T> *)task->func_arg;
      delete tt;
      tt = NULL;
    }
    break;
  case C_DFULL_SYM_GAUSS:
    {
      C_dfull_gauss_arg<T, U> *tt = (C_dfull_gauss_arg<T, U> *)task->func_arg;
      delete tt;
      tt = NULL;
    }
    break;
  case C_DINV_DL_TIMESU:
    {
      C_dinvDL_timesU_arg<T> *tt = (C_dinvDL_timesU_arg<T> *)task->func_arg;
      delete tt;
      tt = NULL;
    }
    break;
  case C_DHALF_SCHUR_B:
    {
      C_dupdateb_Schur_arg<T> *tt = (C_dupdateb_Schur_arg<T> *)task->func_arg;
      delete tt;
      tt = NULL;
    }
    break;
  case C_DHALF_SCHUR_BT:
    {
      C_dupdateb_Schur_arg<T> *tt = (C_dupdateb_Schur_arg<T> *)task->func_arg;
      delete tt;
      tt = NULL;
    }
    break;
  case C_DTRSMSCALE:
    {
      DTRSMScale_arg<T> *tt = (DTRSMScale_arg<T> *)task->func_arg;
      delete tt;
      tt = NULL;
    }
    break;
  case C_DGEMM_LOCAL_MULT:
    {
      DSchurGEMM_arg<T> *tt = (DSchurGEMM_arg<T> *)task->func_arg;
      delete tt;
      tt = NULL;
    }
    break;
  case C_DGEMM_LOCAL_TWO:
    {
      DSchurGEMM_arg<T> *tt = (DSchurGEMM_arg<T> *)task->func_arg;
      delete tt;
      tt = NULL;
    }
    break;
  case C_DGEMM_DIRECT_TWO:
    {
      DSchurGEMM_two_arg<T> *tt = (DSchurGEMM_two_arg<T> *)task->func_arg;
      delete tt;
      tt = NULL;
    }
    break;
  case C_DEALLOCLOWER:
    {
      C_deallocLower_arg<T> *tt = (C_deallocLower_arg<T> *)task->func_arg;
      delete tt;
      tt = NULL;
    }
    break;
  case C_DEALLOCLOCALSCHUR:
    {
      C_deallocLocalSchur_arg<T> *tt = (C_deallocLocalSchur_arg<T> *)task->func_arg;
      delete tt;
      tt = NULL;
    }
    break;

  case C_SPARSESYMFW:
    {
      C_SparseFw_arg<T, U> *tt = (C_SparseFw_arg<T, U> *)task->func_arg;
      delete tt;
      tt = NULL;
    }
    break;
  case C_DSUB_FWBW:
    {
      C_Dsub_FwBw_arg<T> *tt = (C_Dsub_FwBw_arg<T> *)task->func_arg;
      delete tt;
      tt = NULL;
    }
    break;
  case C_DENSE_SYMFW_DIAG:
    {
      C_DenseFwBw_arg<T> *tt = (C_DenseFwBw_arg<T> *)task->func_arg;
      delete tt;
      tt = NULL;
    }
    break;
  case C_DENSE_SYMFW_OFFDIAG:
    {
      C_DenseFwBwOffdiag_arg<T> *tt = (C_DenseFwBwOffdiag_arg<T> *)task->func_arg;
      delete tt;
      tt = NULL;
    }
    break;
    case C_STRIPS_SYMFW_OFFDIAG:
    {
      C_StripsFwBwOffdiag_arg<T> *tt = (C_StripsFwBwOffdiag_arg<T> *)task->func_arg;
      delete tt;
      tt = NULL;
    }
    break;
  case C_DENSE_SYMFILL:
    {
      C_Dfill_FwBw_arg<T> *tt = (C_Dfill_FwBw_arg<T> *)task->func_arg;
      delete tt;
      tt = NULL;
    }
    break;
  case C_SPARSESYMBW:
    {
      C_SparseBw_arg<T, U> *tt = (C_SparseBw_arg<T, U> *)task->func_arg;
      delete tt;
      tt = NULL;
    }
    break;
  case C_DUMMY:
    {
      C_dummy_arg *tt = (C_dummy_arg *)task->func_arg;
      delete tt;
      tt = NULL;
    }
    break;
  }
  if (flag) {
    delete task->ops_complexity; // allocated in C_Dsub_quee()
  }
  //  cout << task->task_name << endl;
  delete task;
  task = NULL;
}
template
void erase_task<double, double>(C_task *& task);

template
void erase_task<quadruple, quadruple>(C_task *& task);

template
void erase_task<complex<double>, double>(C_task *& task);

template
void erase_task<complex<quadruple>, quadruple>(C_task *& task);

template
void erase_task<float, float>(C_task *& task);

template
void erase_task<complex<float>, float>(C_task *& task);
//

template<typename T, typename U>
void full_gauss3(int *n0,
		 T *a,
		 const int n,
		 double *pivot,
		 int *permute,
		 const bool isSym,
		 const double eps,
		 const bool verbose,
		 FILE *fp)
{
  bool flag;
  int nn0, nn1;
  double fop;

  nn1 = *n0;
  diss_printf(verbose, fp,
	      "%s %d : full_sym_gauss3 is not yet implemented %d ",
	      __FILE__, __LINE__, nn1);
  if (isSym) {
    diss_printf(verbose, fp,"ldlt_permute\n");
    flag = full_ldlt_permute<T, U>(&nn0, nn1, n, a, n, pivot, permute, eps,
				   &fop);
  }
  else {
    diss_printf(verbose, fp,"ldu_permute\n");
    flag = full_ldu_permute<T, U>(&nn0, nn1, n, a, n, pivot, permute, eps,
				  &fop);
  }
  *n0 = nn0;
}

template
void full_gauss3<double, double>(int *n0,
				 double *a,
				 const int n,
				 double *pivot,
				 int *permute,
				 const bool isSym,
				 const double eps,
				 const bool verbose,
				 FILE *fp);

template
void full_gauss3<quadruple, quadruple>(int *n0,
				       quadruple *a,
				       const int n,
				       double *pivot,
				       int *permute,
				       const bool isSym,
				       const double eps,
				       const bool verbose,
				       FILE *fp);

template
void full_gauss3<complex<double>, double >(int *n0,
					   complex<double> *a,
					   const int n,
					   double *pivot,
					   int *permute,
					   const bool isSym,
					   const double eps,
					   const bool verbose,
					   FILE *fp);

template
void full_gauss3<complex<quadruple>, quadruple >(int *n0,
						 complex<quadruple> *a,
						 const int n,
						 double *pivot,
						 int *permute,
						 const bool isSym,
						 const double eps,
						 const bool verbose,
						 FILE *fp);

template
void full_gauss3<float, float>(int *n0,
				 float *a,
				 const int n,
				 double *pivot,
				 int *permute,
				 const bool isSym,
				 const double eps,
				 const bool verbose,
				 FILE *fp);

template
void full_gauss3<complex<float>, float >(int *n0,
					   complex<float> *a,
					   const int n,
					   double *pivot,
					   int *permute,
					   const bool isSym,
					   const double eps,
					   const bool verbose,
					   FILE *fp);
//

template<typename T>
void dump_vectors(int nrow, int nn0, T *v, string fname)
{
  fprintf(stderr, "%s %d : speci1alized template is not yet defined.\n",
	   __FILE__, __LINE__); 
}

template<>
void dump_vectors<double>(int nrow, int nn0, double *v, string fname)
{
  FILE *fp;
  if ((fp = fopen(fname.c_str(), "w")) != NULL) {
    for (int i = 0 ; i < nrow; i++) {
      fprintf(fp, "%d ", i);
      for (int j = 0; j < nn0; j++) {
	fprintf(fp, "%g ", v[i + j * nrow]);
      }
      fprintf(fp, "\n");
    }
    fclose(fp);
  }
  else {
    fprintf(stderr,
	    "%s %d : fail to open %s\n",
	    __FILE__, __LINE__, fname.c_str());
  }
}

template
void dump_vectors<complex<double> >(int nrow, int nn0, complex<double> *v,
				    string fname);
template
void dump_vectors<quadruple>(int nrow, int nn0, quadruple *v,
			     string fname);

template
void dump_vectors<complex<quadruple> >(int nrow, int nn0, complex<quadruple> *v,
				       string fname);

template
void dump_vectors<float>(int nrow, int nn0, float *v,
			     string fname);

template
void dump_vectors<complex<float> >(int nrow, int nn0, complex<float> *v,
				    string fname);

//
