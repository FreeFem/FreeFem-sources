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

# ifndef _DRIVER_DISSECTIONMATRIX_
# define _DRIVER_DISSECTIONMATRIX_
#include <complex>
#include <vector>
#include "Splitters/BisectionTree.hpp"
#include "Algebra/PlainMatrix.hpp"
#include "Algebra/SquareMatrix.hpp"
#include "Algebra/SquareBlockMatrix.hpp"
#include "Algebra/RectBlockMatrix.hpp"
#include "Algebra/SparseMatrix.hpp"
#include "Algebra/ColumnMatrix.hpp"
#include "Driver/TridiagBlockMatrix.hpp"
#include "Driver/C_threads_tasks.hpp"

using std::vector;

template<typename T, typename U = T>
class DissectionMatrix
{
public:

  DissectionMatrix(Dissection::Tree *btree,
		   const int nb, const bool isSym, const bool verbose,
		   FILE *fp);

  DissectionMatrix(const DissectionMatrix<T, U> & im)
  {
    (*this) = im;
  }

  ~DissectionMatrix()
  {
    if (_islast) {
      delete [] _color_mask;
      for (int i = 0; i < _colors; i++) {
	_tridiag[i]->free();
	delete _tridiag[i];
      }
      delete [] _tridiag;
    }
    _diag->free(); // free() is compatible to zero-sized matrix
    _lower->free();
    _upper->free();
    _localSchur->free();  // this may be redundant : 26 Jun.2014
    delete _factorize_LDLt;
    delete _diag;
    delete _lower;
    delete _upper;
    delete _localSchur;
  }

  DissectionMatrix& operator = ( const DissectionMatrix<T, U> &im)
  {
    if ( &im != this) {
      _nb = im._nb;
      _level = im._level;
      _nrow = im._nrow;
      _ncol_offdiag = im._ncol_offdiag;
      _nop = im._nop;
      _diag = im._diag->clone();
      if (!_isSym) {
	_lower->copy(*im._lower);
      }
      _upper->copy(*im._upper);
      _isSym = im._isSym; 
      _localSchur = &im.localSchurBlock();
      _factorize_LDLt = im._factorize_LDLt;
      //      _factorize_LDLt_diag = im._factorize_LDLt_diag;
      _csr_diag = im._csr_diag;
      _csr_offdiag = im._csr_offdiag;
      _alignedFather = im._alignedFather;
    }
    return *this;
  }
  int ColorTridiagBlockMatrix() { return _colors ; }
  void setColorTridiagBlockMatrix(int colors) { _colors = colors; }
  int nrow() { return _nrow; }
  SquareBlockMatrix<T>* addrdiagBlock() { return _diag; }
  SquareBlockMatrix<T>*& paddrdiagBlock() { return _diag; }
  SquareBlockMatrix<T>& diagBlock() { return *_diag; }
  const SquareBlockMatrix<T>& diagBlock() const { return *_diag;}
  TridiagBlockMatrix<T, U> **addrtridiagBlock() { return _tridiag; }
  TridiagBlockMatrix<T, U> **&paddrtridiagBlock() { return _tridiag; }
  RectBlockMatrix<T>* addrupperBlock() { return _upper; }
  RectBlockMatrix<T>*& paddrupperBlock() { return _upper; }
  RectBlockMatrix<T>& upperBlock() { return *_upper; }
  const RectBlockMatrix<T>& upperBlock() const { return *_upper; }
  RectBlockMatrix<T>* addrlowerBlock() { return _lower; }
  RectBlockMatrix<T>*& paddrlowerBlock() { return _lower; }
  RectBlockMatrix<T>& lowerBlock() { return *_lower; }
  const RectBlockMatrix<T>& lowerBlock() const { return *_lower; }
  bool isSym() const { return _isSym; }
  void setAlignedFather() { _alignedFather = true; }
  bool isAlignedFather() const { return _alignedFather; }
  bool isFactorized() const  { return _diag->isFactorized(); }
  PlainMatrix<T>& loclSchurBlock();
  const vector<int>& singIdxPermute() const { return _diag->getSingIdx(); }
  vector<int>& singIdxPermute() { return _diag->getSingIdx(); }
  const vector<int>& singIdx() const { return _diag->getSingIdx0(); }
  vector<int>& singIdx() { return _diag->getSingIdx0(); }
  int KernelDetected() const { return _diag->KernelDetected(); }
  int KernelDim() const { return _diag->dim_kern(); }
  void SetKernelDim(int nsing) { _diag->set_dim_kernel(nsing); }
  SquareMatrix<T>& localSchurBlock() { return *_localSchur; }
  const SquareMatrix<T>& localSchurBlock() const { return *_localSchur; }
  T getLastPivot() const { return _diag->lastPivot(); }
  void setLastPivot(T pivot_val) { _diag->set_lastPivot(pivot_val); }
  
  void C_SparseSymbFact_queue(vector<C_task*>& queue,
			      Dissection::Tree *btree,
			      const bool verbose,
			      FILE **fp);
  
  void C_SparseNumFact_queue(vector<C_task*>& queue,
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
			     FILE **fp);

  void C_SparseLocalSchur_queue(vector<C_task*>& queue,
				Dissection::Tree *btree,
				int nnz,
				T *cofes,
				vector<C_task*>& task_p,
				const bool verbose,
				FILE **fp);

  void C_FillMatrix_queue(vector<C_task*>& queue,
			  int nnz,
			  T *cofes,
			  const bool verbose,
			  FILE **fp);

  int C_DFullLDLt_queue(vector<C_task*>& queue,
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
			FILE **fp);

  int C_DTRSMScale_queue(vector<C_task*> &queue,
			  vector<list<int> > &queue_parents_index,
			  vector<int> &queue_index,
			  vector<int>& indcol,
			  vector<int>& ptr,
			  Dissection::Tree *btree,
			  vector<C_task*>& task_o,
			  vector<int>& task_indcol,
			  vector<int>& task_ptr,
			  vector<C_task*>& task_p,
			  const bool verbose,
			  FILE **fp);

  void C_DTRSMScale_rearrange(vector<C_task*> &queue,
			      vector<list<int> > &queue_parents_index,
			      vector<int> &queue_index,
			      vector<int>& indcol,
			      vector<int>& ptr,
			      const bool verbose,
			      FILE *fp);
  
  int C_DGEMM_local_queue(vector<C_task*> &queue,
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
			   FILE **fp);

  void C_deallocLocalSchur_queue(vector<C_task*> &queue,
				 vector<int> &indcol,
				 const bool verbose,
				 FILE **fp);

  void ChildContrib(list<child_contribution<T> > *child_contribs,
		    Dissection::Tree *btree,
		    vector<DissectionMatrix<T, U>* >& dissectionMatrix,
		    FILE **fp);

  void deallocLower_queue(C_task*& queue,
			  bool isSym,
			  vector<C_task*> &task_p,
			  bool isDirect,
			  vector<C_task*> &task_q);

  int getUpperNRow() const { return upperBlock().nbRows(); }
  int getUpperNCol() const { return upperBlock().nbColumns(); }

  unsigned nb() const { return _nb; }
  unsigned level() const { return _level; }  
  int ncol_offdiag() const {return _ncol_offdiag; }
  int nop() const { return _nop; }
  bool islast() const { return _islast; }
  bool alignedFather() const { return _alignedFather; }

  void setNb(unsigned nb) { _nb = nb; }
  void setLevel(unsigned level) { _level = level; }
  void setNrow(int nrow) { _nrow = nrow; }
  void setNcol_offdiag(int ncol_offdiag) { _ncol_offdiag = ncol_offdiag; }
  void setNop(int nop) { _nop = nop; }
  void setIsSym(bool isSym) { _isSym = isSym; }
  void setIslast(bool islast) { _islast = islast; }
  void setAlignedFather(bool alignedFather) { _alignedFather = alignedFather; }

private:
  unsigned _nb;
  unsigned _level;
  int _nrow;
  int _ncol_offdiag;
  int _nop;                            // complexity of factorization
  SquareBlockMatrix<T>* _diag;         // pointer to the diagonal block
  TridiagBlockMatrix<T, U>** _tridiag; // pointer to the tridiagonal blocks
  RectBlockMatrix<T>* _lower;          // pointer to lower off-diagonal block 
  RectBlockMatrix<T>* _upper;          // pointer to upper off-diagonal block
  bool _isSym;             
  bool _islast;
  bool _alignedFather;
  int _colors;
  int *_color_mask;
  
  SquareBlockMatrix<T> *_localSchur;
  ColumnMatrix<T> *_factorize_LDLt;
  //  ColumnMatrix<T> *_factorize_LDLt_diag;
  const CSR_indirect *_csr_diag;
  const CSR_indirect *_csr_offdiag;
};

template<typename T>
class SchurMatrix
{
public:
  SchurMatrix() { }
  ~SchurMatrix() { }
  void free() {
    _sldu.free();
    delete _arow;
    delete _acol;
    _scol.free();
    _schur.free();
  }
  SubSquareMatrix<T>& getSldu() { return _sldu; }
  vector<int>& getSlduList() { return _sldu_list; }
  SparseMatrix<T>*& getArow() { return _arow; }
  SparseMatrix<T>*& getAcol() { return _acol; }
  ColumnMatrix<T>& getScol() { return _scol; }
  ColumnMatrix<T>& getSchur() { return _schur; } // for debugging
private:
  SubSquareMatrix<T>  _sldu;
  vector<int>         _sldu_list;
  SparseMatrix<T>*    _arow;
  SparseMatrix<T>*    _acol;
  ColumnMatrix<T>     _scol;
  ColumnMatrix<T>     _schur;
};
  
template<typename T>
class KernelMatrix
{
public:
  KernelMatrix(int dimension = 0) : _dimension(dimension) { }
  ~KernelMatrix() { }
  void free()
  {
    _singIdx.clear();
    _kern_list_eq.clear();
    _kern_basis.free();
    _tkern_basis.free();
    _kern_proj.free();
    _tkern_proj.free();
    _ntkern_proj.free();
  }
  int dimension() const { return _dimension; }
  void set_dimension(int dimension){ _dimension = dimension; }
  vector<int>& getSingIdx() { return _singIdx; }
  vector<int>& getKernListEq() { return _kern_list_eq; }
  ColumnMatrix<T>& getKernBasis() {return _kern_basis; }
  ColumnMatrix<T>& getTKernBasis() { return _tkern_basis; }
  SquareMatrix<T>& getKernProj() { return _kern_proj; }
  SquareMatrix<T>& getTKernProj()  { return _tkern_proj; }
  SquareMatrix<T>& getNTKernProj() { return _ntkern_proj; }

private:
  int _dimension;
  vector<int>      _singIdx;
  vector<int>      _kern_list_eq;
  ColumnMatrix<T>  _kern_basis;
  ColumnMatrix<T>  _tkern_basis;
  SquareMatrix<T>  _kern_proj;
  SquareMatrix<T>  _tkern_proj; 
  SquareMatrix<T>  _ntkern_proj;
};

#endif
