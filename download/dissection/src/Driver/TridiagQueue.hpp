/*! \file TridiagQueue.hpp
    \brief task mangemanet of tridiagonal factorization algorithm 
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

# ifndef _DRIVER_TRIDIAGQUEUE_
# define _DRIVER_TRIDIAGQUEUE_
#include <complex>
#include <vector>
#include "Driver/C_threads_tasks.hpp"
// #include "Driver/TridiagMatrix.hpp"

using std::vector;

template<typename T, typename U = T>
class TridiagQueue
{
public:
  TridiagQueue(bool tridiag_solver, bool verbose, FILE *fp) :
    _tridiag_solver(tridiag_solver), _verbose(verbose), _fp(fp)
  {
    _allocated = false;
  }

  void generate_queue(TridiagBlockMatrix<T, U> *tridiag,
		      const int dim,
		      const int nnz,
		      const bool isMapped,
		      int *remap_eqn,
		      int *ptUnsymRows,
		      int *indUnsymCol,
		      int *indVals,
		      T *coef);

  void exec_symb_fact();

  void generate_queue_fwbw();
 
  void exec_num_fact(const int called,
		     const double eps_pivot,
		     const bool kernel_detection_all,
		     const int aug_dim,
		     const U eps_machine);

  void exec_fwbw(T *x, const int nrhs, bool isTrans);

  bool tridiagSolver() {
    return _tridiag_solver;
  }
  
  ~TridiagQueue() {
    if (_allocated) {
      delete [] _remap_eqn;
      delete [] _ptRows;
      delete [] _indCols;
      delete [] _indVals;
    }
  }

  bool tridiag_solver() const { return _tridiag_solver; }
  bool isMapped() const {return _isMapped; }
  int dimension() const {return _dim; }
  int nnz() const {return _nnz; }
  int *remap_eqn() {return _remap_eqn; }
  int *ptRows() { return _ptRows; }
  int *indCols() { return _indCols; }
  int *indVals() { return _indVals; }
  
private:
  bool _tridiag_solver;
  TridiagBlockMatrix<T, U> *_tridiag;
  bool _isMapped;
  int _dim;
  int _nnz;
  int *_remap_eqn;
  int *_ptRows;
  int *_indCols;
  int *_indVals;
  T *_coef;
  bool _verbose;
  FILE *_fp;
  bool _allocated;
};

#endif
