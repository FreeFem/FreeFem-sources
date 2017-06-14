/*! \file   C_Dsub.hpp
    \brief  routines for substiution of off-diagonal matrix with strips
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

template<typename T>
void C_Dsub_task_exec(void *arg_); 

template<typename T>
void dsub_sym2sym_diag(C_Dsub_task<T> *arg);

template<typename T>
void dsub_unsym2unsym_diag(C_Dsub_task<T> *arg);

template<typename T>
void dsub_sym2sym(C_Dsub_task<T> *arg);

template<typename T>
void dsub_unsym2unsym(C_Dsub_task<T> *arg);

template<typename T>
void dsub_unsym2diag(C_Dsub_task<T> *arg);

template<typename T>
void dsub_sym2rct(C_Dsub_task<T> *arg);

template<typename T>
void dsub_unsym2rct(C_Dsub_task<T> *arg);

template<typename T>
void dsub_sym2sym_diag_two(C_Dsub_task<T> *arg);

template<typename T>
void dsub_unsym2unsym_diag_two(C_Dsub_task<T> *arg);

template<typename T>
void dsub_sym2sym_two(C_Dsub_task<T> *arg);

template<typename T>
void dsub_unsym2unsym_two(C_Dsub_task<T> *arg);

template<typename T>
void dsub_unsym2diag_two(C_Dsub_task<T> *arg);

template<typename T>
void dsub_sym2rct_two(C_Dsub_task<T> *arg);

template<typename T>
void dsub_unsym2rct_two(C_Dsub_task<T> *arg);

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
		  FILE *fp);

template<typename T>
void update_parents_list(list <int>& parents,
			 const int begin, const int end, 
			 SquareBlockMatrix<T>* mtrx);
