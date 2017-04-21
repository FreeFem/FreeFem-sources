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

#ifndef _C_THREADS_TASKS_
#define _C_THREADS_TASKS_
#include <sys/types.h>
#include <stdint.h>
# include "Compiler/blas.hpp"
# include "Compiler/OptionCompiler.hpp"
#include <vector>
#include <list>
#include <algorithm>
#include <string>
#include <cstring>
#include "Compiler/elapsed_time.hpp"
# include "Algebra/SquareMatrix.hpp"
# include "Algebra/SquareBlockMatrix.hpp"
# include "Algebra/RectBlockMatrix.hpp"
# include "Algebra/ColumnMatrix.hpp"
# include "Splitters/BisectionTree.hpp"
# include "Algebra/CSR_matrix.hpp"
# include "Driver/TridiagBlockMatrix.hpp"

using namespace std;

// set manually the same value defined in Moudels/blak_blocksize.f90

#define SIZE_B1     240
#define SIZE_DGEMM_SYMM_DTRSV   40
//#define SIZE_B1     600
//#define SIZE_DGEMM_SYMM_DTRSV   100

#define C_DHALF_SCHUR_BLOCK_GEMM 1.0
#define C_DHALF_SCHUR_BLOCK_GEMV 2.0
#define FW_SOLVESCLAE_DTRSM      3.0

#define BLAS3_GEMM_OPT 0.5

// #define TASK_NAME_SIZE 256

const string _null_name;

struct source_dist_index {
  int source;
  int dist;
  int global_i;
  int global_j;
  ~source_dist_index() {}
  source_dist_index() {}
  source_dist_index(int source_, int dist_, int global_i_, int global_j_) :
    source(source_), dist(dist_), global_i(global_i_), global_j(global_j_) {}
  source_dist_index(const source_dist_index &im)
  {
    //    if (&im != this) { : no need to be declared
      source = im.source;
      dist = im.dist;
      global_i = im.global_i;
      global_j = im.global_j;
      //    }
  }  
};

template<typename T>
struct C_Dsub_task {
  int ir_bgn;
  int ir_end;
  int jc_bgn; 
  int jc_end;
  int ir_bgn_src;
  int ir_bgn_src2;
  int jc_bgn_src;
  int jc_bgn_src2;
  int dst_row;
  SquareBlockMatrix<T> *dst_mtrx;
  RectBlockMatrix<T> *dst_pt;
  int ir_block;
  int jc_block;
  RectBlockMatrix<T> *dst_pt2;
  int ir2_block;
  int jc2_block;
  SquareBlockMatrix<T> *src_pt;
  SquareBlockMatrix<T> *src_pt2;
  void (*func)(C_Dsub_task<T> *);
  bool isSkip;
  int atomic_size;
  int atomic_id;
  int parallel_max;
  int parallel_id;
  long *ops_complexity;
  int father_id;
  int child0_id;
  int child1_id;
  int level;
  bool verbose;
  bool debug;
  int child_id;
  FILE *fp;
  
  ~C_Dsub_task() {
    delete ops_complexity;
  }
  C_Dsub_task() {}
  C_Dsub_task(int atomic_size_,
	      int atomic_id_,
	      int ir_bgn_,
	      int ir_end_,
	      int jc_bgn_, 
	      int jc_end_,
	      int ir_bgn_src_,
	      int ir_bgn_src2_,
	      int jc_bgn_src_,
	      int jc_bgn_src2_,
	      int dst_row_,
	      SquareBlockMatrix<T>* dst_mtrx_,
	      RectBlockMatrix<T> *dst_pt_,
	      int ir_block_,
	      int jc_block_,
	      RectBlockMatrix<T> *dst_pt2_,
	      int ir2_block_,
	      int jc2_block_,
	      SquareBlockMatrix<T> *src_pt_,
	      SquareBlockMatrix<T> *src_pt2_,
	      void (*func_)(C_Dsub_task<T> *),
	      bool isSkip_,
	      long ops_complexity_,
	      int father_id_,
	      int level_,
	      bool verbose_,
	      FILE *fp_) :
    ir_bgn(ir_bgn_), 
    ir_end(ir_end_),
    jc_bgn(jc_bgn_),
    jc_end(jc_end_),
    ir_bgn_src(ir_bgn_src_),
    ir_bgn_src2(ir_bgn_src2_),
    jc_bgn_src(jc_bgn_src_),
    jc_bgn_src2(jc_bgn_src2_),
    dst_row(dst_row_),
    dst_mtrx(dst_mtrx_),
    dst_pt(dst_pt_),
    ir_block(ir_block_),
    jc_block(jc_block_),
    dst_pt2(dst_pt2_),
    ir2_block(ir2_block_),
    jc2_block(jc2_block_),
    src_pt(src_pt_),
    src_pt2(src_pt2_),
    func(func_),
    isSkip(isSkip_),
    atomic_size(atomic_size_),
    atomic_id(atomic_id_),
    father_id(father_id_),
    level(level_),
    verbose(verbose_),
    fp(fp_)
  {
    debug = false;
    child_id = 0;
    ops_complexity = new long;
    *ops_complexity = ops_complexity_;
  }  

  C_Dsub_task(const C_Dsub_task &im)
  {
    atomic_size = im.atomic_size;
    atomic_id = im.atomic_id;
    ir_bgn = im.ir_bgn; 
    ir_end = im.ir_end;
    jc_bgn = im.jc_bgn;
    jc_end = im.jc_end;
    ir_bgn_src = im.ir_bgn_src;
    ir_bgn_src2 = im.ir_bgn_src2;
    jc_bgn_src = im.jc_bgn_src;
    jc_bgn_src2 = im.jc_bgn_src2;
    dst_row = im.dst_row;
    dst_mtrx = im.dst_mtrx;
    dst_pt = im.dst_pt;
    ir_block = im.ir_block;
    jc_block = im.jc_block;
    dst_pt2 = im.dst_pt2;
    ir2_block = im.ir2_block;
    jc2_block = im.jc2_block;
    src_pt = im.src_pt;
    src_pt2 = im.src_pt2;
    func = im.func;
    isSkip = im.isSkip;
    father_id = im.father_id;
    level = im.level;
    verbose = im.verbose;
    debug = im.debug;
    child_id = im.child_id;
    fp = im.fp;
    ops_complexity = new long;
    *ops_complexity = *(im.ops_complexity); // copy value
  }
};

typedef struct {
  double pivot;
  double w;
  double eps;
  double delta;
} pivot_param;

#define DIST_TASK_CRITICAL 20

#define TASK_WAITING 0
#define TASK_WORKING 1
#define TASK_DONE    2

#define TASK_SINGLE   8
#define TASK_PARALLEL 9

#define C_DUMMY           1
#define C_DIAG_START      2

#define C_DFULL           4
#define C_DFULL_SYM_GAUSS 5
#define C_DINV_DL_TIMESU  6
#define C_DHALF_SCHUR_B   70
#define C_DHALF_SCHUR_BT  71

#define C_DTRSM           8
#define C_DTRSM1          9
#define C_DTRSMSCALE      10

#define C_DGEMM            16
#define C_DGEMM1           17
#define C_DGEMM_LOCAL_MULT 18
#define C_DGEMM_LOCAL_TWO  19
#define C_DGEMM_DIRECT_TWO 20


#define C_DSUB            32
#define C_DSUB1           33
#define C_DEALLOCATE      34
#define C_DEALLOCATE1     35
#define C_DEALLOCLOCALSCHUR 36
#define C_DEALLOCLOWER    40

#define C_FILLMATRIX      64
#define C_FILLMATRIX1     65

#define C_SPARSESYMBFACT 128
#define C_SPARSENUMFACT  256
#define C_SPARSESCHUR    257
#define C_SPARSELOCALSCHUR  512
#define C_SPARSELOCALSCHUR1  513
#define C_SPARSESOLVER   384

#define C_FWBW              1024
#define C_SPARSESYMFW       1025
#define C_SPARSESYMBW       1026
#define C_DSUB_FWBW         1027
#define C_DENSE_SYMFW_DIAG  1028
#define C_DENSE_SYMFW_OFFDIAG  1029
#define C_STRIPS_SYMFW_OFFDIAG  1030
#define C_DENSE_SYMSCALE    1031
#define C_DENSE_SYMFILL     1032

template<typename T, typename U>
struct C_SparseSymbFact_arg {
  TridiagBlockMatrix<T, U> **tridiag;
  int nrow;
  int colors;
  int *color_mask;
  const CSR_indirect* csr_diag;
  long *ops_complexity;
  long *nopd;
  long *nops;
  bool verbose;
  FILE **fp;
 

  ~C_SparseSymbFact_arg() {
    delete ops_complexity;
    delete nopd;
    delete nops;
  }
  C_SparseSymbFact_arg() {}
  C_SparseSymbFact_arg (TridiagBlockMatrix<T, U> **tridiag_,
			int colors_,
			int *color_mask_,
			int nrow_,
			const CSR_indirect* csr_diag_,
			bool verbose_,
			FILE **fp_
			) :
    tridiag(tridiag_),
    nrow(nrow_),
    colors(colors_),
    color_mask(color_mask_),
    csr_diag(csr_diag_),
    verbose(verbose_),
    fp(fp_) {
    ops_complexity = new long;
    nopd = new long;
    nops = new long;
  }

  C_SparseSymbFact_arg(const C_SparseSymbFact_arg &im)
  {
    tridiag = im.tridiag;
    colors = im.colors;
    color_mask = im.color_mask;
    nrow = im.nrow;
    csr_diag = im.csr_diag;
    ops_complexity = new long;
    *ops_complexity = *(im.ops_complexity);
    nopd = new long;
    *nopd = *(im.nopd);
    nops = new long;
    *nops = *(im.nops);
    verbose = im.verbose;
    fp = im.fp;
  }
};

template<typename T, typename U>
struct C_SparseNumFact_arg {
  TridiagBlockMatrix<T, U> **tridiag;
  bool isSym;
  int colors;
  int *color_mask;
  int nnz;
  T *coefs; 
  int nrow;
  int ncol;
  const CSR_indirect* csr_diag;
  const CSR_indirect* csr_offdiag;
  SquareBlockMatrix<T>* D;  // to store pivot information
  double *eps_pivot;
  double *pivot;
  bool *kernel_detection;
  int *dim_aug_kern;
  U *eps_machine;
  SquareBlockMatrix<T>* localSchur;
  long *ops_complexity;
  long *nopd;
  long *nops;
  bool verbose;
  elapsed_t *tt;
  FILE **fp;
  int nb;

  ~C_SparseNumFact_arg() {
    delete ops_complexity;
    delete nopd;
    delete nops;
    delete [] tt;
  }
  C_SparseNumFact_arg() {}
  C_SparseNumFact_arg (TridiagBlockMatrix<T, U> **tridiag_,
		       bool isSym_,
		       int colors_,
		       int *color_mask_,
		       int nnz_,
		       T *coefs_,
		       int nrow_,
		       int ncol_,
		       const CSR_indirect* csr_diag_,
		       const CSR_indirect* csr_offdiag_,
		       SquareBlockMatrix<T>* D_,
		       double *eps_pivot_,
		       double *pivot_,
		       bool *kernel_detection_,
		       int *dim_aug_kern_,
		       U *eps_machine_,
		       SquareBlockMatrix<T>* localSchur_,
		       bool verbose_,
		       FILE **fp_,
		       int nb_) :
    tridiag(tridiag_),
    isSym(isSym_),
    colors(colors_),
    color_mask(color_mask_),
    nnz(nnz_),
    coefs(coefs_),
    nrow(nrow_),
    ncol(ncol_),
    csr_diag(csr_diag_),
    csr_offdiag(csr_offdiag_),
    D(D_),
    eps_pivot(eps_pivot_),
    pivot(pivot_),
    kernel_detection(kernel_detection_),
    dim_aug_kern(dim_aug_kern_),
    eps_machine(eps_machine_),
    localSchur(localSchur_),
    verbose(verbose_),
    fp(fp_),
    nb(nb_)
  {
    ops_complexity = new long;
    nopd = new long;
    nops = new long;
    tt = new elapsed_t[5];
  }

  C_SparseNumFact_arg(const C_SparseNumFact_arg &im)
  {
    tridiag = im.tridiag;
    isSym = im.isSym;
    colors = im.colors;
    color_mask = im.color_mask;
    nnz = im.nnz;
    coefs = im.coefs;
    nrow = im.nrow;
    ncol = im.ncol;
    csr_diag = im.csr_diag;
    csr_offdiag = im.csr_offdiag;
    D = im.D;
    eps_pivot = im.eps_pivot;
    pivot = im.pivot;
    kernel_detection = im.kernel_detection;
    dim_aug_kern = im.dim_aug_kern;
    eps_machine = im.eps_machine;
    localSchur = im.localSchur;
    ops_complexity = new long;
    *ops_complexity = *(im.ops_complexity);
    nopd = new long;
    *nopd = *(im.nopd);
    nops = new long;
    *nops = *(im.nops);
    verbose = im.verbose;
    tt = im.tt;
    fp = im.fp;
    nb = im.nb;
  }
};

template<typename T>
struct C_FillMatrix_arg {
  bool isSym;
  SquareBlockMatrix<T>* D;
  RectBlockMatrix<T>* upper;
  RectBlockMatrix<T>* lower;
  const CSR_indirect *csr_diag;
  const CSR_indirect *csr_offdiag;
  T *coefs;
  long *ops_complexity;
  bool verbose;
  FILE **fp;
  int nb;
  
  ~C_FillMatrix_arg() {
    delete ops_complexity;
  }
  C_FillMatrix_arg() {}
  C_FillMatrix_arg (bool isSym_,
		    SquareBlockMatrix<T>* D_,
		    RectBlockMatrix<T>* upper_,
		    RectBlockMatrix<T>* lower_,
		    const CSR_indirect *csr_diag_,
		    const CSR_indirect *csr_offdiag_,
		    T *coefs_,
		    bool verbose_,
		    FILE **fp_,
		    int nb_) : 
    isSym(isSym_),
    D(D_),
    upper(upper_),
    lower(lower_),
    csr_diag(csr_diag_),
    csr_offdiag(csr_offdiag_),
    coefs(coefs_),
    verbose(verbose_),
    fp(fp_),
    nb(nb_)
  {
    ops_complexity = new long;
  }

  C_FillMatrix_arg(const C_FillMatrix_arg &im)
  {
    isSym = im.isSym;
    D = im.D;
    upper = im.upper;
    lower = im.lower;
    csr_diag = im.csr_diag;
    csr_offdiag = im.csr_offdiag;
    coefs =im.coefs;
    ops_complexity = new long;
    *ops_complexity = *(im.ops_complexity);
    verbose = im.verbose;
    fp = im.fp;
    nb = im.nb;
  }
};

template<typename T, typename U>
struct C_dfull_gauss_arg {
  int task_position;
  int id_level;
  int num_block;
  int id_block;
  SquareBlockMatrix<T>* D;
  SquareBlockMatrix<T>* localSchur;
  ColumnMatrix<T> *a;
  //  ColumnMatrix<T> *diag;    //  T **diag;
  int *permute_block;
  int n;
  int upper_ncol;
  int nrow;
  int i1_block;
  double *eps_piv;
  double *pivot;
  double *pivot0;
  double *pivot1;
  bool isSym;
  long *ops_complexity;
  bool *kernel_detection;
  int *aug_dim;
  U *eps_machine;
  bool *quit;
  int to_next_task;
  bool verbose;
  FILE **fp;
  int nb;

  ~C_dfull_gauss_arg() {
    delete ops_complexity;
  }
  C_dfull_gauss_arg() {}
  C_dfull_gauss_arg(bool isSym_,
		    int task_position_,
		    int id_level_,
		    int num_block_,
		    int id_block_,
		    SquareBlockMatrix<T>* D_,
		    SquareBlockMatrix<T>* localSchur_,
		    ColumnMatrix<T> *a_, //T **a_,
		    //    		    ColumnMatrix<T> *diag_, //T **diag_,
		    int *permute_block_,
		    int n_,
		    int upper_ncol_,
		    int nrow_,
		    int i1_block_,
		    double *eps_piv_,
		    double *pivot_,
		    double *pivot0_,
		    double *pivot1_,
		    bool *kernel_detection_,
		    int *aug_dim_,
		    U *eps_machine_,
		    bool verbose_,
		    FILE **fp_,
		    int nb_) :
    task_position(task_position_),
    id_level(id_level_),
    num_block(num_block_),
    id_block(id_block_),
    D(D_),
    //    lower(lower_),
    localSchur(localSchur_),
    a(a_),
    //    diag(diag_),
    permute_block(permute_block_),
    n(n_),
    upper_ncol(upper_ncol_),
    nrow(nrow_),
    i1_block(i1_block_),
    eps_piv(eps_piv_),
    pivot(pivot_),
    pivot0(pivot0_),
    pivot1(pivot1_),
    isSym(isSym_),
    kernel_detection(kernel_detection_),
    aug_dim(aug_dim_),
    eps_machine(eps_machine_),
    to_next_task(0),
    verbose(verbose_),
    fp(fp_),
    nb(nb_) {
    ops_complexity = new long;
  }

  C_dfull_gauss_arg(const C_dfull_gauss_arg &im)
  {
    isSym = im.isSym;
    task_position = im.task_position;
    id_level = im.id_level;
    num_block = im.num_block;
    id_block = im.id_block;
    D = im.D;
    localSchur = im.localSchur;
    a = im.a;
    //    diag = im.diag;
    permute_block = im.permute_block;
    n = im.n;
    upper_ncol = im.upper_ncol;
    nrow = im.nrow;
    i1_block = im.i1_block;
    eps_piv = im.eps_piv;
    pivot = im.pivot;
    pivot0 = im.pivot0;
    pivot1 = im.pivot1;
    ops_complexity = new long;
    *ops_complexity = *(im.ops_complexity);
    kernel_detection = im.kernel_detection;
    aug_dim = im.aug_dim;
    eps_machine = im.eps_machine;
    fp = im.fp;
    nb = im.nb;
    quit = im.quit;
    to_next_task = im.to_next_task;
  }
};

template<typename T>
struct C_dinvDL_timesU_arg{
  bool isSym;
  int task_position;
  int id_level;
  int num_block;
  int id_block;
  SquareBlockMatrix<T>* D;
  ColumnMatrix<T>* a;  //  T **a;
  int n;
  int nrow;
  int ncol;
  int i1_block;
  int jj_block;
  long *ops_complexity;
  bool verbose;
  FILE **fp;
  int nb;

  ~C_dinvDL_timesU_arg() {
    delete ops_complexity;
  }
  C_dinvDL_timesU_arg() {}
  C_dinvDL_timesU_arg(bool isSym_,
		      int task_position_,
		      int id_level_,
		      int num_block_,
		      int id_block_,
		      SquareBlockMatrix<T>* D_,
		      ColumnMatrix<T> *a_,
		      int n_,
		      int nrow_,
		      int ncol_,
		      int i1_block_,
		      int jj_block_,
		      bool verbose_,
		      FILE **fp_,
		      int nb_) :
    isSym(isSym_),
    task_position(task_position_),
    id_level(id_level_),
    num_block(num_block_),
    id_block(id_block_),
    D(D_),
    a(a_),
    n(n_),
    nrow(nrow_),
    ncol(ncol_),
    i1_block(i1_block_),
    jj_block(jj_block_),
    verbose(verbose_),
    fp(fp_),
    nb(nb_) {
    ops_complexity = new long;
  }

  C_dinvDL_timesU_arg(const C_dinvDL_timesU_arg &im)
  {
    isSym = im.isSym;
    task_position = im.task_position;
    id_level = im.id_level;
    num_block = im.num_block;
    id_block = im.id_block;
    D = im.D;
    a = im.a;
    n = im.n;
    nrow = im.nrow;
    ncol = im.ncol;
    i1_block = im.i1_block;
    jj_block = im.jj_block;
    ops_complexity = new long;
    *ops_complexity = *(im.ops_complexity);
    verbose = im.verbose;
    fp = im.fp;
    nb = im.nb;
  }
};

template<typename T>
struct C_dupdateb_Schur_arg{ 
  bool isSym;
  int task_position;
  int id_level;
  int num_block;
  int id_block;
  SquareBlockMatrix<T>* D;
  ColumnMatrix<T> *a;  //  T **a;
  int n;
  int nrow;
  int ncol;
  int b_size;
  int i1_block;
  int ii_block;
  int jj_block;
  long *ops_complexity;
  bool verbose;
  FILE **fp;
  int nb;

  ~C_dupdateb_Schur_arg() {
    delete ops_complexity;
  }
  C_dupdateb_Schur_arg() {}
  C_dupdateb_Schur_arg(int isSym_,
		       int task_position_,
		       int id_level_,
		       int num_block_,
		       int id_block_,		      
		       SquareBlockMatrix<T>* D_,
		       ColumnMatrix<T> *a_, // T **a_,
		       int n_,
		       int nrow_,
		       int ncol_,
		       int b_size_,
		       int i1_block_,
		       int ii_block_,
		       int jj_block_,
		       bool verbose_,
		       FILE **fp_,
		       int nb_) :
    isSym(isSym_),
    task_position(task_position_),
    id_level(id_level_),
    num_block(num_block_),
    id_block(id_block_),
    D(D_),
    a(a_),
    n(n_),
    nrow(nrow_),
    ncol(ncol_),
    b_size(b_size_),
    i1_block(i1_block_),
    ii_block(ii_block_),
    jj_block(jj_block_),
    verbose(verbose_),
    fp(fp_),
    nb(nb_) {
    ops_complexity = new long;
  }

  C_dupdateb_Schur_arg(const C_dupdateb_Schur_arg &im)
  {
    isSym = im.isSym;
    task_position = im.task_position;
    id_level = im.id_level;
    num_block = im.num_block;
    id_block = im.id_block;
    D = im.D;
    a = im.a;
    n = im.n;
    nrow = im.nrow;
    ncol = im.ncol;
    b_size = im.b_size;
    i1_block = im.i1_block;
    ii_block = im.ii_block;
    jj_block = im.jj_block;
    ops_complexity = new long;
    *ops_complexity = *(im.ops_complexity);
    verbose = im.verbose;
    fp = im.fp;
    nb = im.nb;
  }
};

template<typename T>
struct DTRSMScale_arg {
  SquareBlockMatrix<T>* LDLt;
  RectBlockMatrix<T>* upper; //  double *x;
  RectBlockMatrix<T>* lower; // double *z;
  int nrow;
  int ncol;
  int kblock;
  int lblock;
  int mblock;
  vector<int>* singLocNodes0;
  bool localPermute;
  bool isSym;
  long *ops_complexity;
  bool verbose;
  FILE **fp;
  int nb;
  bool debug;
  
  ~DTRSMScale_arg() {
    delete ops_complexity;
  }
  DTRSMScale_arg() {}
  DTRSMScale_arg(bool isSym_,
		 SquareBlockMatrix<T> *LDLt_,
		 RectBlockMatrix<T>* upper_,  //  double *x_;
		 RectBlockMatrix<T>* lower_,  // double *z_;
		 int nrow_,
		 int ncol_,
		 int kblock_,
		 int lblock_,
		 int mblock_,
		 vector<int>* singLocNodes0_,
		 bool localPermute_,
		 bool verbose_,
		 FILE **fp_,
		 int nb_) :
    LDLt(LDLt_),
    upper(upper_),
    lower(lower_), 
    nrow(nrow_),
    ncol(ncol_),
    kblock(kblock_),
    lblock(lblock_),
    mblock(mblock_),
    singLocNodes0(singLocNodes0_),
    localPermute(localPermute_),
    isSym(isSym_),
    verbose(verbose_),
    fp(fp_),
    nb(nb_)
  {
    debug = false;
    ops_complexity = new long;
  }

  DTRSMScale_arg(const DTRSMScale_arg &im)
  {
    LDLt = im.LDLt;
    upper = im.upper;
    lower = im.lower;
    //    offset = im.offset;
    nrow = im.nrow;
    ncol = im.ncol;
    kblock = im.kblock;
    lblock = im.lblock;
    mblock = im.mblock;
    singLocNodes0 = im.singLocNodes0;
    ops_complexity = new long;
    *ops_complexity = *(im.ops_complexity);
    localPermute = im.localPermute;
    isSym = im.isSym;
    verbose = im.verbose;
    fp = im.fp;
    nb = im.nb;
    debug = im.debug;
  }
};

template<typename T>
struct DSchurGEMM_arg {
  bool isSym;
  bool isTrans;
  RectBlockMatrix<T>* lower; 
  RectBlockMatrix<T>* upper;
  int nrow;
  int i_block; // row block index
  int j_block; // column block index
  SquareBlockMatrix<T> *localSchur;
  long *ops_complexity;
  bool verbose;
  FILE **fp;
  int nb;

  ~DSchurGEMM_arg() {
    delete ops_complexity;
  }
  DSchurGEMM_arg() {}
  DSchurGEMM_arg(bool isSym_,
		 bool isTrans_,
		 RectBlockMatrix<T>* lower_, 
		 RectBlockMatrix<T>* upper_,
		 int nrow_,
		 int i_block_,
		 int j_block_,
		 SquareBlockMatrix<T> *localSchur_,
		 bool verbose_,
		 FILE **fp_,
		 int nb_) :
    isSym(isSym_),
    isTrans(isTrans_),
    lower(lower_),
    upper(upper_),
    nrow(nrow_),
    i_block(i_block_),
    j_block(j_block_),
    localSchur(localSchur_),
    verbose(verbose_),
    fp(fp_),
    nb(nb_)
  {
    ops_complexity = new long;
  }

  DSchurGEMM_arg(const DSchurGEMM_arg &im)
  {
    isSym = im.isSym;
    isTrans = im.isTrans;
    lower = im.lower;
    upper = im.upper;
    nrow = im.nrow;
    i_block = im.i_block;
    j_block = im.j_block;
    localSchur = im.localSchur;
    *ops_complexity = *(im.ops_complexity);
    verbose = im.verbose;
    fp = im.fp;
    nb = im.nb;
  }
};

template<typename T>
struct DSchurGEMM_two_arg {
  bool isSym;
  bool isTrans;
  RectBlockMatrix<T>* lower0; 
  RectBlockMatrix<T>* upper0;
  int nrow0;
  RectBlockMatrix<T>* lower1; 
  RectBlockMatrix<T>* upper1;
  int nrow1;
  int i_block; // row block index
  int j_block; // column block index
  SquareBlockMatrix<T> *localSchur;
  long *ops_complexity;
  bool isSkip;
  bool verbose;
  FILE **fp;
  int nb;

  ~DSchurGEMM_two_arg() {
    delete ops_complexity;
  }
  DSchurGEMM_two_arg() {}
  DSchurGEMM_two_arg(bool isSym_,
		     bool isTrans_,
		     RectBlockMatrix<T>* lower0_, 
		     RectBlockMatrix<T>* upper0_,
		     int nrow0_,
		     RectBlockMatrix<T>* lower1_, 
		     RectBlockMatrix<T>* upper1_,
		     int nrow1_,
		     int i_block_,
		     int j_block_,
		     SquareBlockMatrix<T> *localSchur_,
		     bool isSkip_,
		     bool verbose_,
		     FILE **fp_,
		     int nb_) :
    isSym(isSym_),
    isTrans(isTrans_),
    lower0(lower0_),
    upper0(upper0_),
    nrow0(nrow0_),
    lower1(lower1_),
    upper1(upper1_),
    nrow1(nrow1_),
    i_block(i_block_),
    j_block(j_block_),
    localSchur(localSchur_),
    isSkip(isSkip_),
    verbose(verbose_),
    fp(fp_),
    nb(nb_)
  {
    ops_complexity = new long;
  }

  DSchurGEMM_two_arg(const DSchurGEMM_two_arg &im)
  {
    isSym = im.isSym;
    isTrans = im.isTrans;
    lower0 = im.lower0;
    upper0 = im.upper0;
    nrow0 = im.nrow0;
    lower1 = im.lower1;
    upper1 = im.upper1;
    nrow1 = im.nrow1;
    i_block = im.i_block;
    j_block = im.j_block;
    localSchur = im.localSchur;
    *ops_complexity = *(im.ops_complexity);
    isSkip = im.isSkip;
    verbose = im.verbose;
    fp = im.fp;
    nb = im.nb;
  }
};

template<typename T>
struct C_deallocLower_arg {
  bool isSym;
  RectBlockMatrix<T>* lower; 
  long *ops_complexity;
  ~C_deallocLower_arg() {
    delete ops_complexity;
  }
  C_deallocLower_arg() {}
  C_deallocLower_arg(bool isSym_,
		     RectBlockMatrix<T>* lower_, 
		     long ops_complexity_) :
    isSym(isSym_),
    lower(lower_)
  {
    ops_complexity = new long;
    *ops_complexity = ops_complexity_;
  }
  C_deallocLower_arg(const C_deallocLower_arg &im)
  {
    isSym = im.isSym;
    lower = im.lower;
    ops_complexity = new long;
    *ops_complexity = *(im.ops_complexity);
  }
};

template<typename T>
struct C_deallocLocalSchur_arg {
  bool isSym;
  SquareBlockMatrix<T>* localSchur;
  int i_block;
  int j_block;
  long *ops_complexity;
  bool verbose;
  FILE **fp;
  int nb;
  ~C_deallocLocalSchur_arg() {
    delete ops_complexity;
  }
  C_deallocLocalSchur_arg() {}
  C_deallocLocalSchur_arg(bool isSym_,
			  SquareBlockMatrix<T> *localSchur_,
			  int i_block_,
			  int j_block_,
			  bool verbose_,
			  FILE **fp_,
			  int nb_) :
    isSym(isSym_),
    localSchur(localSchur_),
    i_block(i_block_),
    j_block(j_block_),
    verbose(verbose_),
    fp(fp_),
    nb(nb_) {
    ops_complexity = new long;
  } 
  C_deallocLocalSchur_arg(const C_deallocLocalSchur_arg &im)
  {
    isSym = im.isSym;
    i_block = im.i_block;
    j_block = im.j_block;
    localSchur = im.localSchur;
    ops_complexity = im.ops_complexity;
    verbose = im.verbose;
    fp = im.fp;
    nb = im.nb;
  }
};

struct C_dummy_arg {
  long *ops_complexity;
  bool verbose;
  FILE **fp;
  int nb;
  ~C_dummy_arg() {
    delete ops_complexity;
  }
  C_dummy_arg() {}
  C_dummy_arg(bool verbose_,
	      FILE **fp_,
	      int nb_) :
    verbose(verbose_),
    fp(fp_),
    nb(nb_) {
    ops_complexity = new long;
    *ops_complexity = (-1L);
  }
  C_dummy_arg(const C_dummy_arg & im)
  {
    ops_complexity = im.ops_complexity;
    fp = im.fp;
    nb = im.nb;
  }
};

template<typename T, typename U>
struct C_SparseFw_arg {
  int colors;
  int nb;
  bool isSym;
  int dim;
  bool **isTrans;
  int **nrhs;
  //  void* diag_sparse;
  TridiagBlockMatrix<T, U> **tridiag;
  T **x;
  T **yi;
  T **zi;
  T *coef; // SparseMatrix<double>* ptDA;  //  
  int n_diag;
  int n_offdiag;
  int *ptRows;
  int *indCols;
  int *indVals;
  int *indVals_unsym;
  int *loc2glob_diag;
  long *ops_complexity;
  ~C_SparseFw_arg() {
    delete ops_complexity;
  }
  C_SparseFw_arg() {}
  C_SparseFw_arg(int colors_,
		 int nb_,
		 bool isSym_,
		 int dim_,
		 bool **isTrans_,
		 int **nrhs_,
		 TridiagBlockMatrix<T, U> **tridiag_,
		 T **x_,
		 T **yi_,
		 T **zi_,
		 T *coef_, // SparseMatrix<double>* ptDA_, //
		 int n_diag_,
		 int n_offdiag_,
		 int *ptRows_,
		 int *indCols_,
		 int *indVals_,
		 int *indVals_unsym_,
		 int *loc2glob_diag_) :
    colors(colors_),
    nb(nb_),
    isSym(isSym_),
    dim(dim_),
    isTrans(isTrans_),
    nrhs(nrhs_),
    tridiag(tridiag_),
    x(x_),
    yi(yi_),
    zi(zi_),
    coef(coef_), // ptDA(ptDA_),   //    
    n_diag(n_diag_),
    n_offdiag(n_offdiag_),
    ptRows(ptRows_),
    indCols(indCols_),
    indVals(indVals_),
    indVals_unsym(indVals_unsym_),
    loc2glob_diag(loc2glob_diag_) { 
    ops_complexity = new long;
  }
  C_SparseFw_arg(const C_SparseFw_arg &im)
  {
    colors = im.colors;
    nb = im.nb;
    isSym = im.isSym;
    dim = im.dim;
    isTrans = im.isTrans;
    nrhs = im.nrhs;
    tridiag = im.tridiag;
    x = im.x;
    yi = im.yi;
    zi = im.zi;
    coef = im.coef; // ptDA = im.ptDA; // 
    n_diag = im.n_diag;
    n_offdiag = im.n_offdiag;
    ptRows = im.ptRows;
    indCols = im.indCols;
    indVals = im.indVals;
    indVals_unsym = im.indVals_unsym;
    loc2glob_diag = im.loc2glob_diag;
  }
};

template<typename T, typename U>
struct C_SparseBw_arg {
  int colors;
  int nb;
  bool isSym; // not used, for debugging
  int dim;
  bool **isTrans;
  int **nrhs;
  Dissection::Tree *btree; 
  int level_last;
  TridiagBlockMatrix<T, U> **tridiag;
  T **x;
  T ***yy;
  T **xi;
  T **yi;
  T **zi;
  T *coef; // SparseMatrix<double>* ptDA;  // 
  int *ptRows;
  int *indCols;
  int *indVals;
  int *indVals_unsym;
  long *ops_complexity;
  ~C_SparseBw_arg() {
    delete ops_complexity;
  }
  C_SparseBw_arg() {}
  C_SparseBw_arg(int colors_,
		 int nb_,
		 bool isSym_,
		 int dim_,
		 bool **isTrans_,
		 int **nrhs_,
		 Dissection::Tree *btree_,
		 int level_last_,
		 TridiagBlockMatrix<T, U> **tridiag_,
		 T **x_,
		 T ***yy_,
		 T **xi_,
		 T **yi_,
		 T **zi_,
		 T *coef_, // SparseMatrix<double>* ptDA_, //    
		 int *ptRows_,
		 int *indCols_,
		 int *indVals_,
		 int *indVals_unsym_) :
    colors(colors_),
    nb(nb_),
    isSym(isSym_),
    dim(dim_),
    isTrans(isTrans_),
    nrhs(nrhs_),
    btree(btree_),
    level_last(level_last_),
    //    diag_sparse(diag_sparse_),
    tridiag(tridiag_),
    x(x_),
    yy(yy_),
    xi(xi_),
    yi(yi_),
    zi(zi_),
    coef(coef_), // ptDA(ptDA_), //    
    ptRows(ptRows_),
    indCols(indCols_),
    indVals(indVals_),
    indVals_unsym(indVals_unsym_) { 
    ops_complexity = new long;
  }

  C_SparseBw_arg(const C_SparseBw_arg &im)
  {
    colors = im.colors;
    nb = im.nb;
    isSym = im.isSym;
    dim = im.dim;
    isTrans = im.isTrans;
    nrhs = im.nrhs;
    btree = im.btree;
    level_last = im.level_last;
    tridiag = im.tridiag;
    x = im.x;
    yy = im.yy;
    xi = im.xi;
    yi = im.yi;
    zi = im.zi;
    coef = im.coef; // ptDA = im.ptDA; //
    ptRows = im.ptRows;
    indCols = im.indCols;
    indVals = im.indVals;
    indVals_unsym = im.indVals_unsym;
  }
};

struct C_task {
  int task_id;
  int referred;
  char *task_name;   // is allocated by char [] when structure is constructed
  bool task_name_allocated;
  void *func_arg;
  void (*func)(void *);
  int atomic_size;
  int atomic_id;
  unsigned char status;
  list<C_task*>* parents;      
  list<C_task*>* parents_work;
  int parallel_max;
  int parallel_id;
  long *ops_complexity;  // value is assigned in func_arg->ops_complexity
  long flop;
  elapsed_t t0, t1;
  int thread_id;
  int broadcast_deadlock;
  bool quit_queue;
  int to_next_task;
  bool verbose;
  FILE ***fp;

  ~C_task() {
    if (task_name_allocated) {
      delete [] task_name;
    }
    delete parents;
    delete parents_work;
  }
  C_task() {}
  C_task(int task_id_,
	 string task_name_,
	 void *func_arg_,
	 void (*func_)(void *),
	 int atomic_size_,
	 int atomic_id_,
	 long *ops_complexity_ ) :
    task_id(task_id_),
    referred(0),
    func_arg(func_arg_),
    func(func_),
    atomic_size(atomic_size_),
    atomic_id(atomic_id_),
    status(TASK_WAITING),
    parallel_max(1),
    parallel_id(0),
    ops_complexity(ops_complexity_),
    broadcast_deadlock(0),
    quit_queue(false),
    to_next_task(0) {
    if (!task_name_.empty()) { 
      task_name = new char[task_name_.length() + 1];
      strcpy(task_name, task_name_.c_str());
      task_name_allocated = true;
    }
    else {
      task_name_allocated = false;
    }
    parents = new list<C_task*>;
    parents_work = new list<C_task*>;
    get_realtime(&t0); // reset time as created one
    COPYTIME(t1, t0);
    thread_id = (-1); // not executed yet
  }

  C_task(const C_task &im)
  {
    task_id = im.task_id;
    referred = im.referred;
    if (im.task_name_allocated) { 
      task_name = new char[strlen(im.task_name) + 1];
      strcpy(task_name, im.task_name);
      task_name_allocated = true;
    }
    else {
      task_name_allocated = false;
    }
    func_arg = im.func_arg;
    func = im.func;
    atomic_size = im.atomic_size;
    atomic_id = im.atomic_id;
    status = im.status;
    parents = im.parents;
    parents_work = im.parents_work;
    parallel_max = im.parallel_max;
    parallel_id = im.parallel_id;
    ops_complexity = im.ops_complexity; // copy pointer
    //    nops = im.nops;
    flop = im.flop;
    t0 = im.t0;
    t1 = im.t1;
    thread_id = im.thread_id;
    quit_queue = im.quit_queue;
    to_next_task = im.to_next_task;
    verbose = im.verbose;
    fp = im.fp;
  }
};

struct C_task_seq {
  int task_id;
  int referred;
  char *task_name; // is allocated by char [] when structure is constructed
  bool task_name_allocated;
  int mutex_id;
  int parallel_single; // TASK_PARALLEL/TASK_SINGLE
  int num_threads;
  //  how to describe some (DFullLDLt)s which are shared
  //  assignment of processor group for DFullLDLt is done by statically
  //  5 processors is assigned to 2 + 1 + 1 + 1 / 3 + 2 / 5
  //  with concerning task size
  int level;
  int phase;
  vector <C_task *> *queue;
  int begin;
  int end;
  unsigned char status;
  int pos;
  long ops_complexity;
  
  ~C_task_seq() {
    if (task_name_allocated) { 
      delete [] task_name;
    }
  }
  C_task_seq() {}
  C_task_seq(int task_id_,
	     string task_name_,
	     int mutex_id_,
	     int parallel_single_,
	     int num_threads_,
	     int level_,
	     int phase_,
	     vector <C_task *> *queue_,
	     int begin_,
	     int end_,
	     long ops_complexity_) :
    task_id(task_id_),
    referred(0),
    //    task_name(task_name_),
    mutex_id(mutex_id_),
    parallel_single(parallel_single_),
    num_threads(num_threads_),
    level(level_),
    phase(phase_),
    queue(queue_),
    begin(begin_),
    end(end_),
    status(TASK_WAITING),
    pos(begin_),
    ops_complexity(ops_complexity_) {
    if (!task_name_.empty()) {
      task_name = new char[task_name_.length() + 1];
      strcpy(task_name, task_name_.c_str());
      task_name_allocated = true;
    }
    else {
      task_name_allocated = false;
    }
  } ;

  C_task_seq(const C_task_seq &im)
  {
    task_id = im.task_id;
    referred = im.referred;
    if (im.task_name_allocated) { 
      task_name = new char[strlen(im.task_name) + 1];
      strcpy(task_name, im.task_name);
      task_name_allocated = true;
    }
    else {
      task_name_allocated = false;
    }
    mutex_id = im.mutex_id;
    parallel_single = im.parallel_single;
    num_threads = im.num_threads;
    level = im.level;
    phase = im.phase;
    queue = im.queue; // copy pointer
    begin = im.begin;
    end = im.end;
    status = im.status;
    pos = im.pos;
    ops_complexity = im.ops_complexity;
  }
};



bool C_task_seq_complexity_smaller(C_task_seq *a, C_task_seq *b);
bool C_task_seq_complexity_greater(C_task_seq *a, C_task_seq *b);
bool C_task_seq_beginidx_smaller(C_task_seq *a, C_task_seq *b);
    
struct index_strip {
  int begin_dst;
  int begin_src;
  int width;
  ~index_strip() {}
  index_strip() {}
  index_strip(int begin_dst_, int begin_src_, int width_) :
    begin_dst(begin_dst_),
    begin_src(begin_src_),
    width(width_) {}
  index_strip(const index_strip &im)
  {
    begin_dst = im.begin_dst;
    begin_src = im.begin_src;
    width = im.width;
  }
};

struct index_strip2 {
  int begin_dst;
  int begin_src0;
  int begin_src1;
  int width;
  ~index_strip2() {}
  index_strip2() {}
  index_strip2(int begin_dst_, int begin_src0_, int begin_src1_, int width_) :
    begin_dst(begin_dst_),
    begin_src0(begin_src0_),
    begin_src1(begin_src1_),
    width(width_) {}
  index_strip2(const index_strip2 &im)
  {
    begin_dst = im.begin_dst;
    begin_src0 = im.begin_src0;
    begin_src1 = im.begin_src1;
    width = im.width;
  }
};
  
template<typename T>
struct child_contribution {
  int child_id;
  int diag_size;
  int offdiag_size;
  list <index_strip> diag_strip;
  list <index_strip> offdiag_strip;
  int father_row;
  SquareBlockMatrix<T> *father_diag_pt;
  RectBlockMatrix<T> *father_offdiag_pt;
  RectBlockMatrix<T> *father_offdiag_unsym_pt;
  int child_row;
  SquareBlockMatrix<T>* child_pt;
  //  double **child_pt;
  ~child_contribution() {}
  child_contribution() {}
  child_contribution(int child_id_,
		     int diag_size_,
		     int offdiag_size_,
		     list <index_strip> diag_strip_,
		     list <index_strip> offdiag_strip_,
		     SquareBlockMatrix<T> *father_diag_pt_,
		     RectBlockMatrix<T> *father_offdiag_pt_,
		     RectBlockMatrix<T> *father_offdiag_unsym_pt_,
		     int child_row_,
		     SquareBlockMatrix<T>* child_pt_ ) :
    child_id(child_id_),
    diag_size(diag_size_),
    offdiag_size(offdiag_size_),
    diag_strip(diag_strip_),
    offdiag_strip(offdiag_strip_),
    father_row(diag_size_), //  father_row(father_row_),
    father_diag_pt(father_diag_pt_),
    father_offdiag_pt(father_offdiag_pt_),
    father_offdiag_unsym_pt(father_offdiag_unsym_pt_),
    child_row(child_row_),
    child_pt(child_pt_) {}

  child_contribution(const child_contribution &im)
  {
    child_id = im.child_id;
    diag_size = im.diag_size;
    offdiag_size = im.offdiag_size;
    diag_strip = im.diag_strip;
    offdiag_strip = im.offdiag_strip;
    father_row = im.father_row;
    father_diag_pt = im.father_diag_pt;
    father_offdiag_pt = im.father_offdiag_pt;
    father_offdiag_unsym_pt = im.father_offdiag_unsym_pt;
    child_row = im.child_row;
    child_pt = im.child_pt;
  }
};

struct diag_contribution {
  int child_id;
  int child_row;
  int child_column;
  int father_row;
  list <index_strip> diag_strip;
  ~diag_contribution() {}
  diag_contribution() {}
  diag_contribution(int child_id_,
		    int child_row_,
		    int child_column_,
		    int father_row_,
		    list <index_strip> diag_strip_) :
    child_id(child_id_), 
    child_row(child_row_),
    child_column(child_column_),
    father_row(father_row_),
    diag_strip(diag_strip_) {}

  diag_contribution(const diag_contribution &im)
  {
    child_id = im.child_id;
    child_row = im.child_row;
    child_column = im.child_column;
    father_row = im.father_row;
    diag_strip = im.diag_strip;
  }
};

template<typename T>
struct C_Dsub_FwBw_arg {
  int dim;
  int **nrhs;
  int n_diag;
  bool access_global;
  int level;
  Dissection::Tree *btree;
  list<diag_contribution>* diag_contribs;
  T **x;
  T **yi;
  T ***zi;
  int *loc2glob_diag;
  long *ops_complexity;
  ~C_Dsub_FwBw_arg() {
    delete ops_complexity;
  }
  C_Dsub_FwBw_arg() {}
  C_Dsub_FwBw_arg(int dim_,
		  int **nrhs_,
		  int n_diag_,
		  bool access_global_,
		  int level_,
		  Dissection::Tree *btree_,
		  list<diag_contribution>* diag_contribs_,
		  T **x_,
		  T **yi_,
		  T ***zi_,
		  int *loc2glob_diag_) : 
    dim(dim_),
    nrhs(nrhs_),
    n_diag(n_diag_),
    access_global(access_global_),
    level(level_),
    btree(btree_),
    diag_contribs(diag_contribs_),
    x(x_),
    yi(yi_),
    zi(zi_),
    loc2glob_diag(loc2glob_diag_) {
    ops_complexity = new long;
  }
  C_Dsub_FwBw_arg(const C_Dsub_FwBw_arg &im)
  {
    dim = im.dim;
    nrhs = im.nrhs;
    n_diag = im.n_diag;
    level = im.level;
    access_global = im.access_global;
    btree = im.btree;
    diag_contribs = im.diag_contribs;
    x = im.x;
    yi = im.yi;
    zi = im.zi;
    loc2glob_diag = im.loc2glob_diag;
    ops_complexity = im.ops_complexity;
  }    
};

template<typename T>
struct C_Dfill_FwBw_arg {
  int **nrhs;
  int d;
  int level;
  Dissection::Tree *btree;
  int n_offdiag;
  T ***yi;
  T **zi;
  long *ops_complexity;
  ~C_Dfill_FwBw_arg() {
    delete ops_complexity;
  }
  C_Dfill_FwBw_arg() {}
  C_Dfill_FwBw_arg(int **nrhs_,
		   int d_,
		   int level_,
		   Dissection::Tree *btree_,
		   int n_offdiag_,
		   T ***yi_,
		   T **zi_) :
    nrhs(nrhs_),
    d(d_),
    level(level_),
    btree(btree_),
    n_offdiag(n_offdiag_),
    yi(yi_),
    zi(zi_) {
    ops_complexity = new long;
  }
  C_Dfill_FwBw_arg(const C_Dfill_FwBw_arg &im)
  {
    nrhs = im.nrhs;
    d = im.d;
    level = im.level;
    btree = im.btree;
    n_offdiag = im.n_offdiag;
    yi = im.yi;
    zi = im.zi;
    ops_complexity = im.ops_complexity;
  }    
};

template<typename T>
struct C_DenseFwBw_arg
{
  bool isSym;
  bool isBackward;
  int dim;
  bool **isTrans;
  int **nrhs;
  int n_diag;
  int nrow;
  int k_block;
  int kk;
  T **xi;
  T **wi;
  T **yi;
  T **x;
  SquareBlockMatrix<T> *LDLt;
  bool isFirstBlock;
  bool isLastBlock;
  int *loc2glob;
  long *ops_complexity;
  bool verbose;
  FILE **fp;
  ~C_DenseFwBw_arg() {
    delete ops_complexity;
  }
  C_DenseFwBw_arg() {}

  C_DenseFwBw_arg(bool isSym_,
		  bool isBackward_,
		  int dim_,
		  bool **isTrans_,
		  int **nrhs_,
		  int n_diag_,
		  int nrow_,
		  int k_block_,
		  int kk_,
		  T **xi_,
		  T **wi_,
		  T **yi_,
		  T **x_,
		  SquareBlockMatrix<T> *LDLt_,
		  bool isFirstBlock_,
		  bool isLastBlock_,
		  int *loc2glob_,
		  bool verbose_,
		  FILE **fp_) :
    isSym(isSym_),
    isBackward(isBackward_),
    dim(dim_),
    isTrans(isTrans_),
    nrhs(nrhs_),
    n_diag(n_diag_),
    nrow(nrow_),
    k_block(k_block_),
    kk(kk_),
    xi(xi_),
    wi(wi_),
    yi(yi_),
    x(x_),
    LDLt(LDLt_),
    isFirstBlock(isFirstBlock_),
    isLastBlock(isLastBlock_),
    loc2glob(loc2glob_),
    verbose(verbose_),
    fp(fp_)
  {
    ops_complexity = new long;
  }
  C_DenseFwBw_arg(C_DenseFwBw_arg &im)
  {
    isSym = im.isSym;
    isBackward = im.isBackward;
    dim = im.dim;
    isTrans = im.isTrans;
    nrhs = im.nrhs;
    n_diag = im.n_diag;
    nrow = im.nrow;
    k_block = im.k_block;
    kk = im.kk;
    xi = im.xi;
    wi = im.wi;
    yi = im.yi; 
    x = im.x;
    LDLt = im.LDLt;
    isFirstBlock = im.isFirstBlock;
    isLastBlock = im.isLastBlock;
    loc2glob = im.loc2glob;
    verbose = im.verbose;
    fp = im.fp;
  }    
};

template<typename T>
struct C_DenseFwBwOffdiag_arg
{
  bool trans;
  bool isLower;
  int dim;
  bool **isTrans;
  int **nrhs;
  int lda;
  int ldb;
  int ldc;
  int nrow;
  int ncol;
  T **xi;
  int ii;
  T **yi;
  T **zi;
  int jj;
  SquareBlockMatrix<T> *LDLt;
  int i_block;
  int j_block;
  T alpha;
  T beta;
  long *ops_complexity;

  ~C_DenseFwBwOffdiag_arg() {
    delete ops_complexity;
  }
  C_DenseFwBwOffdiag_arg() {}

  C_DenseFwBwOffdiag_arg(bool trans_,
			 bool isLower_,
			 int dim_,
			 bool **isTrans_,
			 int **nrhs_,
			 int lda_,
			 int ldb_,
			 int ldc_,
			 int nrow_,
			 int ncol_,
			 T **xi_,
			 int ii_,
			 T **yi_,
			 T **zi_,
			 int jj_,
			 SquareBlockMatrix<T> *LDLt_,
			 int i_block_,
			 int j_block_,
			 T alpha_,
			 T beta_) :
    trans(trans_),
    isLower(isLower_),
    dim(dim_),
    isTrans(isTrans_),
    nrhs(nrhs_),
    lda(lda_),
    ldb(ldb_),
    ldc(ldc_),
    nrow(nrow_),
    ncol(ncol_),
    xi(xi_),
    ii(ii_),
    yi(yi_),
    zi(zi_),
    jj(jj_),
    LDLt(LDLt_),
    i_block(i_block_),
    j_block(j_block_),
    alpha(alpha_),
    beta(beta_)
  {
    ops_complexity = new long;
  }
  C_DenseFwBwOffdiag_arg(C_DenseFwBwOffdiag_arg &im)
  {
    trans = im.trans;
    isLower = im.isLower;
    dim = im.dim;
    isTrans = im.isTrans;
    nrhs = im.nrhs;
    lda = im.lda;
    ldb = im.ldb;
    ldc = im.ldc;
    nrow = im.nrow;
    ncol = im.ncol;
    xi = im.xi;
    ii = im.ii;
    yi = im.yi; 
    zi = im.zi; 
    jj = im.jj;
    LDLt = im.LDLt;
    i_block = im.i_block;
    j_block = im.j_block;
    alpha = im.alpha;
    beta = im.beta;
    ops_complexity = im.ops_complexity;
  }    
};

template<typename T>
struct C_StripsFwBwOffdiag_arg
{
  bool isLower;
  int dim;
  bool **isTrans;
  int **nrhs;
  int lda;
  int ldb;
  int ldc;
  int nrow;
  int ncol;
  T **xi;
  int ii;
  T **yi;
  int jj;
  RectBlockMatrix<T> *upper;
  RectBlockMatrix<T> *lower;
  int i_block;
  int j_block;
  T alpha;
  T beta;
  long *ops_complexity;

  ~C_StripsFwBwOffdiag_arg() {
    delete ops_complexity;
  }
  C_StripsFwBwOffdiag_arg() {}

  C_StripsFwBwOffdiag_arg(bool isLower_,
			  int dim_,
			  bool **isTrans_,
			  int **nrhs_,
			  int lda_,
			  int ldb_,
			  int ldc_,
			  int nrow_,
			  int ncol_,
			  T **xi_,
			  int ii_,
			  T **yi_,
			  int jj_,
			  RectBlockMatrix<T> *upper_,
			  RectBlockMatrix<T> *lower_,
			  int i_block_,
			  int j_block_,
			  T alpha_,
			  T beta_) :
    isLower(isLower_),
    dim(dim_),
    isTrans(isTrans_),
    nrhs(nrhs_),
    lda(lda_),
    ldb(ldb_),
    ldc(ldc_),
    nrow(nrow_),
    ncol(ncol_),
    xi(xi_),
    ii(ii_),
    yi(yi_),
    jj(jj_),
    upper(upper_),
    lower(lower_),
    i_block(i_block_),
    j_block(j_block_),
    alpha(alpha_),
    beta(beta_)
  {
    ops_complexity = new long;
  }
  C_StripsFwBwOffdiag_arg(C_StripsFwBwOffdiag_arg &im)
  {
    isLower = im.isLower;
    dim = im.dim;
    isTrans = im.isTrans;
    nrhs = im.nrhs;
    lda = im.lda;
    ldb = im.ldb;
    ldc = im.ldc;
    nrow = im.nrow;
    ncol = im.ncol;
    xi = im.xi;
    ii = im.ii;
    yi = im.yi; 
    jj = im.jj;
    upper = im.upper;
    lower = im.lower;
    i_block = im.i_block;
    j_block = im.j_block;
    alpha = im.alpha;
    beta = im.beta;
    ops_complexity = im.ops_complexity;
  }    
};


#define RESTRICT 

void print_strips(const char *name, list<index_strip> &strips, FILE *fp);
void print_strips(const char *name, list<index_strip2> &strips, FILE *fp);

void assign_tasks_statically(list<C_task_seq*> *queue_static,
			     list<C_task_seq*> &queue_dynamic,
			     long *nops_sum,
			     list<C_task_seq *> &task_seq_tmp,
			     const int task_id,
			     const char *task_name_,
			     const int level,
			     const int phase,
			     const long nops_block_total,
			     int num_threads);

template<typename T>
int count_diag_negative(SquareBlockMatrix<T>& Diag);
template<>
int count_diag_negative<double>(SquareBlockMatrix<double>& Diag);
template<>
int count_diag_negative<complex<double> >(SquareBlockMatrix<complex<double> >& Diag);
template<>
int count_diag_negative<quadruple>(SquareBlockMatrix<quadruple>& Diag);
template<>
int count_diag_negative<complex<quadruple> >(SquareBlockMatrix<complex<quadruple> >& Diag);

template<typename T>
int count_diag_negative(SubSquareMatrix<T>& Diag);
template<>
int count_diag_negative<double>(SubSquareMatrix<double>& Diag);
template<>
int count_diag_negative<complex<double> >(SubSquareMatrix<complex<double> >& Diag);
template<>
int count_diag_negative<quadruple>(SubSquareMatrix<quadruple>& Diag);
template<>
int count_diag_negative<complex<quadruple> >(SubSquareMatrix<complex<quadruple> >& Diag);

template<typename T, typename U>
void full_gauss3(int *n0,
		 T *a,
		 const int n,
		 double *pivot,
		 int *permute,
		 const bool isSym,
		 const double eps,
		 const bool verbose,
		 FILE *fp);

void CopyUpper2LowerSymm(double RESTRICT *s, const double RESTRICT *s_t, 
			 const int size_b2,
			 int i, int block_nrow, int ncol);

void CopyUpper2LowerSquare(double RESTRICT *s, const double RESTRICT *s_t, 
			   const int size_b2,
			   int i, int j, int block_nrow, int block_ncol,
			   int ncol);

template<typename T, typename U>
void C_SparseSymbFact(void *arg_); 

template<typename T, typename U>
void C_SparseNumFact(void *arg_);

template<typename T, typename U>
void C_SparseLocalSchur(void *arg_);

template<typename T>
void dump_matrix(FILE *fp, const int nrow, T *a);
template<>
void dump_matrix<double>(FILE *fp, const int nrow, double *a);

template<typename T>
void dump_matrix(FILE *fp, const int nrow, const int ncol, T *a);
template<>
void dump_matrix<double>(FILE *fp, const int nrow, const int ncol, double *a);

template<typename T>
void dump_matrix(FILE *fp, const int kk, 
		 const int nrow, const int ncol, const int nn, T *a);

template<typename T>
void dump_matrix(FILE *fp, 
		 RectBlockMatrix<T> &a);

template<typename T>
void dump_matrix(FILE *fp, 
		 SquareBlockMatrix<T> &a);

template<typename T>
void dump_matrix(FILE *fp, const int nrow, const int nnz, int *prow,
		 int *indcols, int *indvals, T *a);
#if 0
template<typename T>
void verify_nan(FILE *fp, const int nnz, T *a);
#endif
template<typename T, typename U>
void C_dfull_gauss_b(void *arg_);

template<typename T>
void C_dinvDL_timesU(void *arg_);

template<typename T>
void C_dupdateb_Schur_diag(void *arg_);

template<typename T>
void C_dupdateb_Schur_offdiag(void *arg_);

template<typename T, typename U>
void C_gauss_whole_pivot(void *arg_);

template<typename T>
void C_dupdateb_Schur_offdiag_t(void *arg_);

void C_CopyLowerMatrix(int n, const double *A, double *LDlt);

void C_DTRSMScale_solve_seq(void *arg_);

template<typename T>
void C_FillSymMatrix(void *arg_);

template<typename T>
void C_FillMatrix_diag(void *arg_);

template<typename T>
void C_FillMatrix_offdiag(void *arg_);

template<typename T>
void DSchurGEMM_diag(void *arg_);

template<typename T>
void DSchurGEMM_diag_two(void *arg_);


template<typename T>
void DSchurGEMM_offdiag(void *arg_);

template<typename T>
void DSchurGEMM_offdiag_two(void *arg_);

template<typename T>
void C_DTRSMScale_diag_upper(void *arg_);

template<typename T>
void C_DTRSMScale_offdiag_upper(void *arg_);

template<typename T>
void C_DTRSMScale_diag_lower(void *arg_);

template<typename T>
void C_DTRSMScale_offdiag_lower(void *arg_);

template<typename T>
void C_DTRSMScale_solve(void *arg_);

template<typename T>
void C_deallocLower(void *arg_);

template<typename T>
void C_deallocLocalSchur(void *arg_);

template<typename T, typename U>
void C_SparseFw(void *arg_);

template<typename T, typename U>
void C_SparseBw(void *arg_);

template<typename T>
void C_Dsub_FwBw(void *arg_);

template<typename T>
void C_DenseFwBw_diag(void *arg_);

template<typename T>
void C_DenseFwBw_offdiag(void *arg_);

template<typename T>
void C_StripsFwBw_offdiag(void *arg_);

template<typename T>
void C_Dfill_FwBw(void *arg_);

template<typename T, typename U>
void erase_task(C_task *& task);

void C_dummy(void *arg_);

#define imin(a, b) ((a) < (b) ? (a) : (b))
#define imax(a, b) ((a) > (b) ? (a) : (b))

int compare_DTRSMScale_task(const void *_a, const void *_b);
int compare_DSchurGEMM_task(const void* _a, const void* _b);
int compare_source_dist_index(const void *_a, const void *_b);

int convert_array2strip(list<index_strip> &strips, 
			vector<int>& array);

int combine_two_strips(list<index_strip> &stripsa,
		       list<index_strip> &stripsb,
		       list<index_strip2> &stripsc,
		       list<index_strip> &strips0, 
		       list<index_strip> &strips1,
		       const int size);
#if 0
void copy_one_strip(list<index_strip> &strips_dst, 
		    list<index_strip> &strips_src);
#endif
void copy_two_strips(list<index_strip2> &strips2,
		     list<index_strip> &strips0, 
		     list<index_strip> &strips1);

void split_two_strips(list<index_strip> &strips0,
		      list<index_strip> &strips1,
		      list<index_strip2> &strips2,
		      index_strip strip0,
		      index_strip strip1);

template<typename T, typename U>
int dimKernDense(vector<int> &singIdx,
		 const int n,
		 const int aug_dim,
		 const U eps_machine,
		 const double eps_piv,
		 SquareBlockMatrix<T> &D,
		 T *a,
		 const bool refactorize,
		 const bool isBlocked,
		 const bool isSym,
		 const bool verbose,
		 FILE *fp);

template<typename T>
void calc_relative_norm(double *norm_l2,  double *norm_infty, 
			const T *v, const T *u, const int dim);
template<>
void calc_relative_norm<double>(double *norm_l2,  double *norm_infty, 
				const double *v, const double *u,
				const int dim);
template<>
void calc_relative_norm<quadruple>(double *norm_l2,  double *norm_infty, 
				   const quadruple *v, const quadruple *u,
				   const int dim);
template<>
void calc_relative_norm<complex<double> >(double *norm_l2,  double *norm_infty, 
					  const complex<double> *v,
					  const complex<double> *u,
					  const int dim);
template<>
void calc_relative_norm<complex<quadruple> >(double *norm_l2,
					     double *norm_infty, 
					     const complex<quadruple> *v,
					     const complex<quadruple> *u,
					     const int dim);

template<typename T, typename U>
void calc_relative_normscaled(double *norm_l2,  double *norm_infty, 
			      const T *v, const T *u, const U *w,
			      const int dim);
template<>
void calc_relative_normscaled<double,
			      double>(double *norm_l2,  double *norm_infty, 
				      const double *v, const double *u,
				      const double *w,
				      const int dim);
template<>
void calc_relative_normscaled<quadruple,
			      quadruple>(double *norm_l2,  double *norm_infty, 
					 const quadruple *v, const quadruple *u,
					 const quadruple *w,
					 const int dim);

template<>
void calc_relative_normscaled<quadruple,
			      double>(double *norm_l2,  double *norm_infty, 
				      const quadruple *v, const quadruple *u,
				      const double *w,
				      const int dim);

template<>
void calc_relative_normscaled<complex<double>,
			      double>(double *norm_l2,  double *norm_infty, 
				      const complex<double> *v,
				      const complex<double> *u,
				      const double *w,
				      const int dim);
template<>
void calc_relative_normscaled<complex<quadruple>,
			      quadruple>(double *norm_l2,
					 double *norm_infty, 
					 const complex<quadruple> *v,
					 const complex<quadruple> *u,
					 const quadruple *w,
					 const int dim);

template<>
void calc_relative_normscaled<complex<quadruple>,
			      double>(double *norm_l2,
				      double *norm_infty, 
				      const complex<quadruple> *v,
				      const complex<quadruple> *u,
				      const double *w,
				      const int dim);

int CSR_sym2unsym(CSR_indirect *unsym, 
		  const int *ptSymRows, const int *indSymCols, 
		  const int *map_eqn, const int *remap_eqn,
		  const int dim, 
		  const bool upper_flag = true);

bool CSR_unsym2unsym(CSR_indirect *unsym,
		     const int *ptUnSymRows, const int *indUnSymCols, 
		     const int *map_eqn, const int *remap_eqn, 
		     const int dim);

void swap_queues_n(vector <C_task *> &queue,
		   vector<int> &queue_index,
		   const int ii, 
		   const int jj,
		   const int n,
		   vector <C_task *> &tmp,
		   vector<int> &tmp_index);

int EraseNullParents(vector<C_task *> &queue);
int EraseNullParents(C_task *task);

bool compare_task_name(C_task *first, C_task *second);

extern "C" {
  void c_getrealtime_(uint64_t &tmprofiles, const int &m);
  void c_fileout_(uint64_t &fp_prt, char *strgs, const int &force_stderr);
}

void swap_2x2pivots(const int way, 
		      int *pivot_width, int *permute_q, 
		      const int dim_augkern, const int nn0, 
		      const int n_dim, double *a1, long double *aq,
		      double *d1, long double *d1q, double *a_fact);
#if 0
  void test_2x2(int *print_cntrl);
#endif

template<typename T> void dump_vectors(int nrow, int nn0, T *v,
				       string fname);

void ComputeSVD(double *b, const double *a_, const int n);
#endif
