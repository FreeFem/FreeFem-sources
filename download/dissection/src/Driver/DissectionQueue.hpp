/*! \file   DissectionQueue.hpp
    \brief  management of threads for factorization and Fw/Bw substitution
    \author Atsushi Suzuki, Laboratoire Jacques-Louis Lions
    \date   Mar. 30th 2012
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

#ifndef _DISSECTION_QUEUE_
#define _DISSECTION_QUEUE_
#include <cassert>
#include "Compiler/OptionLibrary.h"
#include "Splitters/BisectionTree.hpp"
#include "Driver/DissectionMatrix.hpp"
#include "Driver/QueueRuntime.hpp"
#include "Compiler/elapsed_time.hpp"
#include <time.h>
#include <string>
#include <pthread.h>

template<typename T, typename U = T>
class DissectionQueue
{
public:
  DissectionQueue(Dissection::Tree *btree,
		  vector<DissectionMatrix<T, U>*>& dissectionMatrix,
		  const int num_threads,
		  const bool isSym,
		  const bool verbose,
		  FILE *fp);

  ~DissectionQueue(); 

  void generate_queue(vector<DissectionMatrix<T, U>*>& dM,
		      int nnz, T *coefs);

  void generate_queue_fwbw(vector<DissectionMatrix<T, U>*>& dissectionMatrix,
			   int dim, int nnz, T *coefs);

  void exec_symb_fact();

  void exec_num_fact(const int called,
		     const double eps_piv, 
		     const bool kernel_detection,
		     const int aug_dim,
		     const U eps_machine);

  void exec_fwbw(T *x, const int nrhs, bool isTrans);
  void exec_fwbw_seq(T *x, const int nrhs, bool isTrans);

  void erase_queue(void);
  void erase_queue_fwbw(void);

  DissectionQueue(const DissectionQueue &s) 
  {
    _queue_symb = s._queue_symb;
    _queue_static = s._queue_static;
    _queue_dynamic = s._queue_dynamic;
    _pivots = s._pivots;
  }

  list<child_contribution<T> >& ChildContribs(int nb) {
    return _child_contribs[nb];
  }

  int dimension() const { return _dim; }
  int nnz() const { return _nnz; }
  
private:
  Dissection::Tree* _btree;

  C_task_seq* _queue_symb;     // single set of tasks among sparse subdomains
  C_task_seq* _queue_fwbw;
  list<C_task_seq*> *_queue_static;
  vector<C_task_seq *> *_queue_dynamic;
  list<C_task *> _queue_dummy;
  vector<int>* _children;
  vector<C_task *>* _tasks_SparseSymb;
  vector<C_task *>* _tasks_SparseNum;
  vector<C_task *>* _tasks_SparseLocalSchur;
  vector<C_task *>* _tasks_DFillSym;
  vector<C_task *>* _tasks_DFullLDLt;
  vector<C_task *>* _tasks_DTRSMScale;
  vector<C_task *>* _tasks_DSymmGEMM;
  vector<C_task *>* _tasks_deallocLocalSchur;
  vector<C_task *>** _tasks_Dsub;
  vector<C_task *>* _tasks_deallocLower;
  list<child_contribution<T> > *_child_contribs;

  int _num_threads;
  int _num_threads_symb;
  int _dim;
  int _nnz;
  int _nb_doms;
  int _nb_level;
 
  double _eps_piv;           // used in selecting pivot
  bool _kernel_detection;
  int _aug_dim;
  U _eps_machine;            // magnitude of numerical perturbation
  double *_pivots;
  int _max_size_work_y;
  list<diag_contribution>* _diag_contribs;
  T ***_xi;
  T ***_yi;
  T ***_zi;
  T ***_wi;
  T **_x;   // a pointer to keep rhs and solution vectors

  int **_nrhs;    // these values contain only address information during
  bool **_isTrans;
  // queue generation of fw/bw subtitutions
  bool _verbose;
  FILE *_fp;
  FILE **_fps;
  bool _isSym;

  QueueRuntime* _dissectionRuntime;
  bool _queue_symb_allocated;
  bool _queue_numrc_allocated;
  bool _queue_fwbw_allocated;
    
  static const T _one;  // (1.0);
  static const T _zero; // (0.0);
  static const T _none; // (-1.0);
}; // End class DissictionQueue

#endif
