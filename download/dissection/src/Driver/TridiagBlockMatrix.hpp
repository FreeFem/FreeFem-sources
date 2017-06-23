/*! \file TridiagBlockMatrix.hpp
    \brief tridiagonal factorization algorithm with Cuthill-McKee
    \author François-Xavier Roux, ONERA, Laboratoire Jacques-Louis Lions
    \author Atsushi Suzuki, Laboratoire Jacques-Louis Lions
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

#ifndef _ALGEBRA_TRIDIAGBLOCKMATRIX_HPP_
#define _ALGEBRA_TRIDIAGBLOCKMATRIX_HPP_

#include <vector>
#include <list>
#include "Compiler/blas.hpp"
#include "Compiler/elapsed_time.hpp"
#include "Algebra/ColumnMatrix.hpp"
#include "Algebra/SquareBlockMatrix.hpp"
#include "Algebra/CSR_matrix.hpp"

// setting for SX-ACE to reduce small tridiagonal blocks during DTRSM
// #define SPARSE_OFFDIAG   // CSR data is used in offdiag blocks of tridiag
// #define STORE_WHOLE_FACTORIZED
// setting for superscalar CPUs
#define PERMUTE_UPPER

using std::vector;
using std::list;

template<typename T, typename U = T>
class TridiagBlockMatrix
{
public:
  TridiagBlockMatrix(int dim, int size_b1, bool isSymmetric,
		     const int nb, const bool verbose, FILE *fp) :
    _diag_block_alloc_status(false), _nb(nb), _verbose(verbose), _fp(fp) {
    init(dim, size_b1, isSymmetric);
  }
  
  void init(int dim, int size_b1, bool isSymmetric) {
    _dim = dim;
    _size_b1 = size_b1;
    _isSymmetric = isSymmetric;
  }
    
  ~TridiagBlockMatrix() {
    free();
    if(_diag_block_alloc_status) {
      delete [] _diag_blocks;
#ifndef SPARSE_OFFDIAG
      delete [] _upper_blocks; // 27 Nov.2015
      delete [] _lower_blocks; // 27 Nov.2015
#endif
#ifdef STORE_WHOLE_FACTORIZED
      _factorized_whole.free();
#endif
      _diag_block_alloc_status = false;
    }
  }

  void free() {
      for (int n = 0; n < _nfront; n++) {
      _diag_blocks[n].free();
#ifndef SPARSE_OFFDIAG
      _upper_blocks[n].free(); // 27 Nov.2015
      _lower_blocks[n].free(); // 27 Nov.2015
#endif
    }
    _a21.free();
    _a12.free();
    _s22.free();
  }

  bool isSym() const {
    return _isSymmetric;
  }    
  void SymbolicFact(const int color,
		    const int color_max,
		    int *color_mask,
		    const int dim_,
		    const int nnz,
		    const int *prow_,
		    const int *indcols_,
		    const int *indvals_);

  void NumericFact(const T* coef,
		   const double eps_pivot,
		   double *pivot,
		   const bool kernel_detection,
		   const int dim_aug_kern,
		   const U eps_machine,
		   double *nopd);

  void ComputeSchurComplement(const int nsing,
			      vector<int> &list_sing,
			      ColumnMatrix<T>& a12,
			      ColumnMatrix<T>& a21,
			      ColumnMatrix<T>& a22);

  void TridiagNumericFact(double *pivot,
			  const double eps_pivot,
			  const int dim_aug_kern,
			  ColumnMatrix<T> *diag_block_save,
			  vector<int> &num_null_aug,
			  double *nopd);
  
  void ComputeSchur(const int dim_,
		    int *color_mask,
		    const int ncol,
		    const int *ptrow1, const int *indcols1,
		    const int *indvals1, const int *indvals2, 
		    const T *coef,
		    const int size_b1,
		    SquareBlockMatrix<T> & local_s,
		    double *nopd,
		    elapsed_t *tt);

  void SolveMulti(const bool flag_new2old, const bool isTrans,
		  const int nhrs, ColumnMatrix<T>& x, const int nscol_ = (-1));
  
  void SolveSingle(const bool flag_new2old,
		   const bool isTrans, T* x, const int nscol_ = (-1));

  void ForwardUpper(bool isTransposed,
		    int ncol, ColumnMatrix<T> &b,
		    vector<int>& i0,
		    vector<T> &dscale);

  void extract_column(const int jcol, T *dcol);

  void extract_row(const int irow, T *drow);

  void SingularNode(vector<int> &list_sing);

  int NumNegativeDiags();

  // void KernelIndex(vector<int> &list);

  void KernelBasis(const bool isTrans, ColumnMatrix<T> &a12);

  int dimension() const {return _dim; }
  int block_size() const {return _size_b1; }
  int nrow() const { return _dim; }
  int rank() const { return (_dim - _n0); }
  int nsing() const { return _n0; }
  bool detected() const { return _detected; }
  double nop() const {return _nop; }
  bool diag_block_alloc_status() const {return _diag_block_alloc_status; }
  int maxdim() const { return _maxdim; }
  int nnz() const { return _nnz; }
  int nscol() const {return _nscol; }
  int Nfront() const {return _nfront; }
  
  vector<int> &getPtRows()  { return _ptRows; }
  vector<int> &getIndCols() { return _indCols; }
  vector<int> &getIndVals() { return _indVals; }
  vector<int> &getNew2old() { return _new2old; }
  vector<int> &getP_front() { return _p_front; }
  vector<int> &getPermute() { return _permute; }
  vector<int> &getPermute_ginv() { return _permute_ginv; }
  vector<int> &getP_diag() { return _p_diag; }
  vector<int> &getP_upper() { return _p_upper; }
  vector<int> &getList_schur() { return _list_schur; }
  vector<int> &getList_elim() { return _list_elim; }
  vector<int> &getNum_null() { return _num_null; }
  void setCoef(const T *coef) { _coef = coef; }
  ColumnMatrix<T> &getA12() { return _a12; }
  ColumnMatrix<T> &getA21() { return _a21; }
  ColumnMatrix<T> &getS22() { return _s22; }

  void setNfront(int nfront) { _nfront = nfront; }
  void setMaxdim(int maxdim) { _maxdim = maxdim; }
  void setNnz(int nnz) { _nnz = nnz; }
  void setNop(double nop) { _nop = nop; }
  void setDiag_block_alloc_status(int diag_block_alloc_status) {
    _diag_block_alloc_status = diag_block_alloc_status;
  }
  void setNscol(int nscol) { _nscol = nscol; }
  void setNsing(int n0) { _n0 = n0; }
  void setDetected(bool detected) { _detected = detected; }

  //  void setNb(int nb) { _nb = nb; } // for debug
  int getNb() { return _nb; } // for debug
  ColumnMatrix<T>*& getaddrDiagMatrix() { return _diag_blocks; }
#ifndef SPARSE_OFFDIAG
  ColumnMatrix<T>*& getaddrLowerMatrix() { return _lower_blocks; }
  ColumnMatrix<T>*& getaddrUpperMatrix() { return _upper_blocks; }
#endif
private:
  int _nb;
  bool _verbose;
  FILE *_fp;
  bool _isSymmetric;
  const T* _coef;
  int _color;
  int _color_max;
  int _dim;
  int _size_b1;
  int _nfront;
  int _maxdim;
  int _nnz;
  double _nop;
  vector<int> _ptRows;  // prow
  vector<int> _indCols; // column_numb
  vector<int> _indVals; // column_numb
  vector<int> _new2old;
  vector<int> _p_front;
  vector<int> _permute;
  vector<int> _permute_ginv;
  vector<int> _p_diag;
  vector<int> _p_upper;
  ColumnMatrix<T>* _diag_blocks;
#ifndef SPARSE_OFFDIAG
  ColumnMatrix<T>* _upper_blocks;  // 27 Nov.2015
  ColumnMatrix<T>* _lower_blocks;  // 27 Nov.2015
#endif
#ifdef STORE_WHOLE_FACTORIZED
  ColumnMatrix<T> _factorized_whole;
#endif
  bool _diag_block_alloc_status;
  int _nscol;   // number of postponed entries
  int _n0;      // kernel dimension
  bool _detected;
  ColumnMatrix<T> _a12;
  ColumnMatrix<T> _a21;
  ColumnMatrix<T> _s22;
  vector<int> _list_schur;
  vector<int> _list_elim; // ?
  vector<int> _num_null;

  static const T _one;  // (1.0);
  static const T _zero; // (0.0);
  static const T _none; // (-1.0);
};

void RenumberCSR(const bool shrink_flag,
		 const int dim,
		 const int nnz,
		 vector<int> &b2a,
		 vector<int> &a2b,
		 const int *aptrows,
		 const int *aindcols,
		 const int *aindvals,
		 vector<int> &bptrows,
		 vector<int> &bindcols, vector<int> &bindvals);

void RenumberCSR(const int dim,
		 vector<int> &b2a,
		 vector<int> &a2b,
		 const int *aptrows,
		 const int *aindcols,
		 vector<int> &bptrows,
		 vector<int> &bindcols);

void TridiagStruct(vector<int> &ptrow, vector<int> &indcols,
		   const int nfront,
		   vector<int> &pfront, vector<int> &p_diag,
		   vector<int> &p_upp);

void GenPermuteOffdiag(const int pfront0, const int pfront1, const int pfront2,
		       vector<int> &ptrow, vector<int> &ptdiag,
		       vector<int> &indcols,
		       const int *permute_diag, int *permute_diag_inv,
		       int *permute_offdiag, int *permute_offdiag_inv,
		       vector<int> &i0);

void GenPermuteUpper(const int nrow,
		     vector<int> &remap_eqn,
		     const int ncol,
		     const int *ptrow, const int *indcols,
		     const int n_front,
		     vector<int> &p_fronts,
		     vector<int> &old2new,
		     vector<int> &permute_diag,
		     vector<int> &permute_diag_inv,
		     vector<int> &permute_upper,
		     vector<int> &permute_upper_inv,
		     vector<int> &i0,
		     const bool verbose,
		     FILE *fp);

template<typename T>
void FillBlockSparse(const T *coef,
		     const int dim,
		     vector<int>& map_eqn,
		     const int *ptrow,
		     const int *indcols, const int *indvals,
		     vector<int> &new2old_j,
		     vector<int> &old2new_j,
		     ColumnMatrix<T> &b);

#endif
