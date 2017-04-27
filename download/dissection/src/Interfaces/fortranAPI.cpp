/*! \file   fortranAPI.cpp
    \brief  Fortran style interface named as Dissecotion-fortran interface
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
#include <sys/types.h>
#include <stdint.h>
#include <cstdio>
#include <cstdlib>

#ifdef DISSECTION_FORTRAN
#include "Interfaces/fortranAPI.h"
#else
#include "Interfaces/Dissection.hpp"
#endif
#include "Driver/DissectionSolver.hpp"

struct dissection_solver_ptr
{
  int real_or_complex;
  bool quad_fact;
  DissectionSolver<double> *rptr;
  DissectionSolver<complex<double>, double> *cptr;
  DissectionSolver<quadruple, quadruple, double, double> *rqtr;
  DissectionSolver<double, double, double, double, quadruple, quadruple> *rqtr_fwbw;
  DissectionSolver<complex<quadruple>, quadruple,
		   complex<double>, double,
		   complex<quadruple>, quadruple> *cqtr;
  DissectionSolver<complex<double>, double,
		   complex<double>, double,
		   complex<quadruple>, quadruple> *cqtr_fwbw;
  FILE *fp;
  bool verbose;
  int called;
#ifdef BLAS_MKL
  int mkl_num_threads;
#endif
  int symbolic;
  int numeric;
};

DISSECTION_API void DISS_VERSION(int *versn,
				 int *reles,
				 int *patch)
{
  *versn = DISSECTION_VERSION;
  *reles = DISSECTION_RELEASE;
  *patch = DISSECTION_PATCHLEVEL;
}

DISSECTION_API void DISS_INIT(uint64_t &dslv_,
			      const int &called,
			      const int &real_or_complex,
			      const int &quad_fact,
			      const int &nthreads,
			      const int &verbose)
{
  int num_threads;
  dissection_solver_ptr *dslv;
  dslv_ = (uint64_t)new dissection_solver_ptr;
  dslv = (dissection_solver_ptr *)dslv_;
  dslv->real_or_complex = real_or_complex;
  dslv->quad_fact = (quad_fact == 1);
  dslv->called = called;
  dslv->symbolic = 0;
  dslv->numeric = 0;
  {
    int pid = (int)getpid();
    char fname[256];
    if (verbose > 0) {
      dslv->verbose = true;
    }
    else {
      dslv->verbose = false;
    }
#if 1 
    if (dslv->verbose > 0) {
      fprintf(stderr, "pid = %d\n", pid);
      sprintf(fname, "dissection.%04d.%04d.log", pid, called);
      //      sprintf(fname, "dissection.%04d.log", pid);
      dslv->fp = fopen(fname, "a");
    }
    else {
      dslv->fp = stderr;
    }
#else
    dslv->fp = stderr;
#endif
  }
  if (dslv->verbose > 0) {
    fprintf(dslv->fp, "%s %d : diss_init : called = %d\n", 
	    __FILE__, __LINE__,  dslv->called);
  }
  
  //  _called++;                   // counter for dumping matrix data to debug
#ifdef BLAS_MKL
  if (getenv("MKL_NUM_THREADS")) {
    sscanf(getenv("MKL_NUM_THREADS"), "%d", &dslv->mkl_num_threads);
    if (dslv->verbose > 0) {
      fprintf(dslv->fp,
	      "environmental variable MKL_NUM_THREADS = %d\n",
	      dslv->mkl_num_threads);
    }
  }
  else {
    dslv->mkl_num_threads = mkl_get_max_threads();
  }
  if (dslv->verbose > 0) {
    fprintf(dslv->fp,
	    "MKL_NUM_THREADS = %d\n", dslv->mkl_num_threads);
  }
#endif
  if (nthreads == (-1)) {
    if (getenv("NTHREADS")) {
      sscanf(getenv("NTHREADS"), "%d", &num_threads);
    }
    else {
      num_threads = 1;
    }
  }
  if (nthreads > 0) {
    num_threads = nthreads;
  }
  if (dslv->quad_fact) {
    switch(real_or_complex) {
    case DISSECTION_REAL_MATRIX:
      dslv->rqtr = new DissectionSolver<quadruple, quadruple,
				    double, double>(num_threads, 
						    (verbose != 0 ? true : false), 
						    dslv->called, dslv->fp);
      dslv->rqtr_fwbw = new DissectionSolver<double, double,
					     double, double,
					     quadruple, quadruple>(num_threads, 
								   (verbose != 0 ? true : false), 
								   dslv->called, dslv->fp);

      break;
    case DISSECTION_COMPLEX_MATRIX:
      dslv->cqtr = new DissectionSolver<complex<quadruple>, quadruple,
					complex<double>, double,
					complex<quadruple>, quadruple>(num_threads, 
						  (verbose != 0 ? true : false), 
						  dslv->called, dslv->fp);
      dslv->cqtr_fwbw = new DissectionSolver<complex<double>, double,
					     complex<double>, double,
					     complex<quadruple>, quadruple>(num_threads, 
						  (verbose != 0 ? true : false), 
						  dslv->called, dslv->fp);

      break;
    }
  }
  else {
    switch(real_or_complex) {
    case DISSECTION_REAL_MATRIX:
      dslv->rptr = new DissectionSolver<double>(num_threads, 
						(verbose != 0 ? true : false), 
						dslv->called, dslv->fp);
      break;
    case DISSECTION_COMPLEX_MATRIX:
      dslv->cptr = new DissectionSolver<complex<double>, double>(num_threads, 
								 (verbose != 0 ? true : false), 
								 dslv->called, dslv->fp);
      break;
    default:
      if (dslv->verbose > 0) {
	fprintf(dslv->fp, "%s %d : unknown matrix data type : %d\n", 
		__FILE__, __LINE__, dslv->real_or_complex);
      }
    }
  }
}

DISSECTION_API void DISS_FREE(uint64_t &dslv_) 
{
  dissection_solver_ptr *dslv = (dissection_solver_ptr *)dslv_;

  if (dslv->quad_fact) {
    switch(dslv->real_or_complex) {
    case DISSECTION_REAL_MATRIX:
      delete dslv->rqtr;
      //      delete dslv->rqtr_fwbw; 
      break;
    case DISSECTION_COMPLEX_MATRIX:
      delete dslv->cqtr;
      //      delete dslv->cqtr_fwbw; 
      break;
    default:
      if (dslv->verbose > 0) {
	fprintf(dslv->fp, "%s %d : unknown matrix data type : %d\n", 
		__FILE__, __LINE__, dslv->real_or_complex);
      }
    }
  }
  else {
    switch(dslv->real_or_complex) {
    case DISSECTION_REAL_MATRIX:
      delete dslv->rptr;
      break;
    case DISSECTION_COMPLEX_MATRIX:
      delete dslv->cptr;
      break;
    default:
      if (dslv->verbose > 0) {
	fprintf(dslv->fp, "%s %d : unknown matrix data type : %d\n", 
		__FILE__, __LINE__, dslv->real_or_complex);
      }
    }
  }
  delete dslv;
  dslv_ = (uint64_t)NULL;
  //  _called--;
  //  if ((_called == 0) && (dslv->fp != stderr)) {
  if (dslv->fp != stderr) {
    fclose(dslv->fp);
  }
#ifdef VECLIB
  unsetenv("VECLIB_MAXIMUM_THREADS");
#endif

}

DISSECTION_API void DISS_NUMERIC_FREE(uint64_t &dslv_)
{
  dissection_solver_ptr *dslv = (dissection_solver_ptr *)dslv_;
  if (dslv->quad_fact) {
    switch(dslv->real_or_complex) {
    case DISSECTION_REAL_MATRIX:
      dslv->rqtr->NumericFree();
      break;
    case DISSECTION_COMPLEX_MATRIX:
      dslv->cqtr->NumericFree();
      break;
    default:
      if (dslv->verbose > 0) {
	fprintf(dslv->fp, "%s %d unknown matrix data type : %d\n", 
		__FILE__, __LINE__, dslv->real_or_complex);
      }
    }
  }
  else {
    switch(dslv->real_or_complex) {
    case DISSECTION_REAL_MATRIX:
      dslv->rptr->NumericFree();
      break;
    case DISSECTION_COMPLEX_MATRIX:
      dslv->cptr->NumericFree();
      break;
    default:
      if (dslv->verbose > 0) {
	fprintf(dslv->fp, "%s %d unknown matrix data type : %d\n", 
		__FILE__, __LINE__, dslv->real_or_complex);
      }
    }
  }
}
  
DISSECTION_API void DISS_S_FACT(uint64_t &dslv_,
				const int &dim,
				const int *ptRows,
				const int *indCols,
				const int &sym,
				const int &decomposer)
{
  // sym = 1 : symmetric with upper
  //     = 0 : unsymmetric,
  //     = 3 : symmetric with lower
  // decomposer = 0 : SCOTCH 
  //            = 1 : METIS
  //            = 2 : TRIDAIG without nested bisection
  dissection_solver_ptr *dslv = (dissection_solver_ptr *)dslv_;
  //  int num_levels;
#if 0
  double *coefs;
  const int nnz = ptRows[dim];
  coefs = new double[nnz];
  memset(coefs, 0, sizeof(double) * nnz);
  SaveMMMatrix_(dim, nnz, false, false, ptRows, indCols, dslv->called, coefs);
  delete [] coefs;
#endif
  
#if 0
  switch(dslv->real_or_complex) {
  case DISSECTION_REAL_MATRIX:
    fout = dslv->rptr->get_filedescriptor();
    break;
  case DISSECTION_COMPLEX_MATRIX:
    fout = dslv->cptr->get_filedescriptor();
    break;
  default:
    if (dslv->verbose > 0) {
      fprintf(dslv->fp, "%s %d unknown matrix data type : %d\n", 
	      __FILE__, __LINE__, dslv->real_or_complex);
    }
  }
#endif
  if (dslv->quad_fact) {
    switch(dslv->real_or_complex) {
    case DISSECTION_REAL_MATRIX:
      dslv->rqtr->SymbolicFact(dim, (int *)ptRows, (int *)indCols,
			       (bool)(sym % 2 == 1), (bool)((sym / 2) == 0),
			       (bool)((sym / 4) == 1),
			       decomposer); // using default parameter
      break;
    case DISSECTION_COMPLEX_MATRIX:
      dslv->cqtr->SymbolicFact(dim, (int *)ptRows, (int *)indCols,
			       (bool)(sym % 2 == 1), (bool)((sym / 2) == 0),
			       (bool)((sym / 4) == 1),
			       decomposer); // using default parameter
      break;
    default:
      if (dslv->verbose > 0) {
	fprintf(dslv->fp, "%s %d unknown matrix data type : %d\n", 
		__FILE__, __LINE__, dslv->real_or_complex);
      }
    }
    if (dslv->verbose > 0) {
      fprintf(dslv->fp, "%s:%d Dissection::SymbolicFact done\n",
	      __FILE__, __LINE__);
    }
  }
  else {
    switch(dslv->real_or_complex) {
    case DISSECTION_REAL_MATRIX:
      dslv->rptr->SymbolicFact(dim, (int *)ptRows, (int *)indCols,
			       (bool)(sym % 2 == 1), (bool)((sym / 2) == 0),
			       (bool)((sym / 4) == 1),
			       decomposer); // using default parameter
      break;
    case DISSECTION_COMPLEX_MATRIX:
      dslv->cptr->SymbolicFact(dim, (int *)ptRows, (int *)indCols,
			       (bool)(sym % 2 == 1), (bool)((sym / 2) == 0),
			       (bool)((sym / 4) == 1),
			       decomposer); // using default parameter
      break;
    default:
      if (dslv->verbose > 0) {
	fprintf(dslv->fp, "%s %d unknown matrix data type : %d\n", 
		__FILE__, __LINE__, dslv->real_or_complex);
      }
    }
    if (dslv->verbose > 0) {
      fprintf(dslv->fp, "%s:%d Dissection::SymbolicFact done\n",
	      __FILE__, __LINE__);
    }
  }
  dslv->symbolic++;
}

DISSECTION_API void DISS_N_FACT(uint64_t &dslv_,
				const double *coefs,
				const int &scaling,
				const double &eps_pivot,
				const int &indefinite_flag)
{
  // scaling = 0 : without scaling
  //           1 : 1/sqrt(a_ii) or 1/sqrt(max|a_ij|)
  //           2 : 1/sqrt(a_ii) or Schur complement corresponding to diagonal
  //               kernel_detection_all (for KKT type)
  // eps_pivot = 1.0e-2 : threshold of pivot, ratio of contiguous diagonal
  //                      entries with absolute value
  // indefinite_flag = 1 : indefinite -> kernel_detection_all = false
  // indefinite_flag = 0 : semi-definite -> kernel_detection_all = true

  bool kernel_detection_all = indefinite_flag == 0 ? true : false;
    
  //  FILE *fout;
  dissection_solver_ptr *dslv = (dissection_solver_ptr *)dslv_;

#ifdef BLAS_MKL
  mkl_set_num_threads(1);
#endif
  #ifdef VECLIB
    setenv("VECLIB_MAXIMUM_THREADS", "1", true);
#endif
  //  dslv->NumericFree(); // for debugging : 20 Nov.2013
    if (dslv->quad_fact) {
      switch(dslv->real_or_complex) {
	case DISSECTION_REAL_MATRIX:
	  dslv->rqtr->NumericFact(dslv->numeric,
				  (double *)coefs, scaling, 
				  eps_pivot, 
				  kernel_detection_all);
	  dslv->rqtr_fwbw->CopyQueueFwBw(*(dslv->rqtr));
	  break;
      case DISSECTION_COMPLEX_MATRIX:
	dslv->cqtr->NumericFact(dslv->numeric,
				(complex<double> *)coefs, scaling, 
				eps_pivot, 
				kernel_detection_all);
	dslv->cqtr_fwbw->CopyQueueFwBw(*(dslv->cqtr));
	break;
      }
    }
    else {
      switch(dslv->real_or_complex) {
	case DISSECTION_REAL_MATRIX:
	  dslv->rptr->NumericFact(dslv->numeric,
				  (double *)coefs, scaling, 
				  eps_pivot, 
				  kernel_detection_all);
	  if (dslv->rptr->getFactorized() == false) {
	    dslv->rptr->SaveMMMatrix(dslv->called, coefs);
	  }
	  break;
      case DISSECTION_COMPLEX_MATRIX:
	dslv->cptr->NumericFact(dslv->numeric,
				(complex<double> *)coefs, scaling, 
				eps_pivot, 
				kernel_detection_all);
	break;
      }
    }
#ifdef BLAS_MKL
  mkl_set_num_threads(dslv->mkl_num_threads);
#endif
  if (dslv->verbose > 0) {
    fprintf(dslv->fp,
	    "%s %d : Dissection::NumericFact done : %d\n",
	    __FILE__, __LINE__, dslv->numeric);
  }
  dslv->numeric++;
}

DISSECTION_API void DISS_GET_COLORS(uint64_t &dslv_, 
				    int *n)
{
  dissection_solver_ptr *dslv = (dissection_solver_ptr *)dslv_;
  if (dslv->quad_fact) {
    switch(dslv->real_or_complex) {
    case DISSECTION_REAL_MATRIX:
      *n = dslv->rqtr->GetMaxColors();
      break;
    case DISSECTION_COMPLEX_MATRIX:
      *n =  dslv->cqtr->GetMaxColors();
      break;
    default:
      if (dslv->verbose > 0) {
	fprintf(dslv->fp, "%s %d unknown matrix data type : %d\n", 
		__FILE__, __LINE__, dslv->real_or_complex);
      }
    }
  }
  else {
    switch(dslv->real_or_complex) {
    case DISSECTION_REAL_MATRIX:
      *n = dslv->rptr->GetMaxColors();
      break;
    case DISSECTION_COMPLEX_MATRIX:
      *n =  dslv->cptr->GetMaxColors();
      break;
    default:
      if (dslv->verbose > 0) {
	fprintf(dslv->fp, "%s %d unknown matrix data type : %d\n", 
		__FILE__, __LINE__, dslv->real_or_complex);
      }
    }
  }
}

DISSECTION_API void DISS_GET_KERN_DIM(uint64_t &dslv_, 
				      int *n0)
{
  dissection_solver_ptr *dslv = (dissection_solver_ptr *)dslv_;

  if (dslv->quad_fact) {
    switch(dslv->real_or_complex) {
    case DISSECTION_REAL_MATRIX:
      *n0 = dslv->rqtr->kern_dimension();
      break;
    case DISSECTION_COMPLEX_MATRIX:
      if (dslv->verbose > 0) {
	fprintf(dslv->fp, 
		"%s %d diss_get_kern_dim() for complex is not yet implemented\n",
		__FILE__, __LINE__);
      }
      *n0 = 0;
      break;
    default:
      if (dslv->verbose > 0) {
	fprintf(dslv->fp, "%s %d unknown matrix data type : %d\n", 
		__FILE__, __LINE__, dslv->real_or_complex);
      }
    }
  }
  else {
    switch(dslv->real_or_complex) {
    case DISSECTION_REAL_MATRIX:
      *n0 = dslv->rptr->kern_dimension();
      break;
    case DISSECTION_COMPLEX_MATRIX:
      if (dslv->verbose > 0) {
	fprintf(dslv->fp, 
		"%s %d diss_get_kern_dim() for complex is not yet implemented\n",
		__FILE__, __LINE__);
      }
      *n0 = 0;
      break;
    default:
      if (dslv->verbose > 0) {
	fprintf(dslv->fp, "%s %d unknown matrix data type : %d\n", 
		__FILE__, __LINE__, dslv->real_or_complex);
      }
    }
  }
}

DISSECTION_API void DISS_GET_NULLPIVOTS(uint64_t &dslv_, 
					int *pivots)
{
  dissection_solver_ptr *dslv = (dissection_solver_ptr *)dslv_;
  if (dslv->quad_fact) {
    switch(dslv->real_or_complex) {
    case DISSECTION_REAL_MATRIX:
      dslv->rqtr->GetNullPivotIndices(pivots);
      break;
    case DISSECTION_COMPLEX_MATRIX:
      dslv->cqtr->GetNullPivotIndices(pivots);
      break;
    default:
      if (dslv->verbose > 0) {
	fprintf(dslv->fp, "%s %d unknown matrix data type : %d\n", 
		__FILE__, __LINE__, dslv->real_or_complex);
      }
    }
  }
  else {
    switch(dslv->real_or_complex) {
    case DISSECTION_REAL_MATRIX:
      dslv->rptr->GetNullPivotIndices(pivots);
      break;
    case DISSECTION_COMPLEX_MATRIX:
      dslv->cptr->GetNullPivotIndices(pivots);
      break;
    default:
      if (dslv->verbose > 0) {
	fprintf(dslv->fp, "%s %d unknown matrix data type : %d\n", 
		__FILE__, __LINE__, dslv->real_or_complex);
      }
    }
  }
}

DISSECTION_API void DISS_GET_SMALLPIVOTS(uint64_t &dslv_,
					 const int &n,
					 int *pivots)
{
  dissection_solver_ptr *dslv = (dissection_solver_ptr *)dslv_;
  if (dslv->quad_fact) {
    switch(dslv->real_or_complex) {
    case DISSECTION_REAL_MATRIX:
      dslv->rqtr->GetSmallestPivotIndices(n, pivots);
      break;
    case DISSECTION_COMPLEX_MATRIX:
      dslv->cqtr->GetSmallestPivotIndices(n, pivots);
      break;
    default:
      if (dslv->verbose > 0) {
      fprintf(dslv->fp, "%s %d unknown matrix data type : %d\n", 
	      __FILE__, __LINE__, dslv->real_or_complex);
      }
    }
  }
  else {
    switch(dslv->real_or_complex) {
    case DISSECTION_REAL_MATRIX:
      dslv->rptr->GetSmallestPivotIndices(n, pivots);
      break;
    case DISSECTION_COMPLEX_MATRIX:
      dslv->cptr->GetSmallestPivotIndices(n, pivots);
      break;
    default:
      if (dslv->verbose > 0) {
      fprintf(dslv->fp, "%s %d unknown matrix data type : %d\n", 
	      __FILE__, __LINE__, dslv->real_or_complex);
      }
    }
  }
}

DISSECTION_API void DISS_GET_KERN_VECS(uint64_t &dslv_, 
				       double *vec)
{
  dissection_solver_ptr *dslv = (dissection_solver_ptr *)dslv_;
  if (dslv->quad_fact) {
    switch(dslv->real_or_complex) {
    case DISSECTION_REAL_MATRIX:
      dslv->rqtr_fwbw->GetKernelVectors(vec);
      break;
    case DISSECTION_COMPLEX_MATRIX:
      if (dslv->verbose > 0) {
	fprintf(dslv->fp, 
		"%s %d diss_get_kern_vecs() for complex is not yet implemented\n",
		__FILE__, __LINE__);
      }
      break;
    default:
      if (dslv->verbose > 0) {
	fprintf(dslv->fp, "%s %d unknown matrix data type : %d\n", 
		__FILE__, __LINE__, dslv->real_or_complex);
      }
    }
  }
  else {
    switch(dslv->real_or_complex) {
    case DISSECTION_REAL_MATRIX:
      dslv->rptr->GetKernelVectors(vec);
      break;
    case DISSECTION_COMPLEX_MATRIX:
      if (dslv->verbose > 0) {
	fprintf(dslv->fp, 
		"%s %d diss_get_kern_vecs() for complex is not yet implemented\n",
		__FILE__, __LINE__);
      }
      break;
    default:
      if (dslv->verbose > 0) {
	fprintf(dslv->fp, "%s %d unknown matrix data type : %d\n", 
		__FILE__, __LINE__, dslv->real_or_complex);
      }
    }
  }
}

DISSECTION_API void DISS_GET_KERNT_VECS(uint64_t &dslv_, 
				       double *vec)
{
  dissection_solver_ptr *dslv = (dissection_solver_ptr *)dslv_;
  if (dslv->quad_fact) {
    switch(dslv->real_or_complex) {
    case DISSECTION_REAL_MATRIX:
      dslv->rqtr_fwbw->GetTransKernelVectors(vec);
      break;
    case DISSECTION_COMPLEX_MATRIX:
      if (dslv->verbose > 0) {
	fprintf(dslv->fp, 
		"%s %d diss_get_kern_vecs() for complex is not yet implemented\n",
		__FILE__, __LINE__);
      }
      break;
    default:
      if (dslv->verbose > 0) {
	fprintf(dslv->fp, "%s %d unknown matrix data type : %d\n", 
		__FILE__, __LINE__, dslv->real_or_complex);
      }
    }
  }
  else {
    switch(dslv->real_or_complex) {
    case DISSECTION_REAL_MATRIX:
      dslv->rptr->GetTransKernelVectors(vec);
      break;
    case DISSECTION_COMPLEX_MATRIX:
      if (dslv->verbose > 0) {
	fprintf(dslv->fp, 
		"%s %d diss_get_kern_vecs() for complex is not yet implemented\n",
		__FILE__, __LINE__);
      }
      break;
    default:
      if (dslv->verbose > 0) {
	fprintf(dslv->fp, "%s %d unknown matrix data type : %d\n", 
		__FILE__, __LINE__, dslv->real_or_complex);
      }
    }
  }
}

DISSECTION_API void DISS_PROJECT(uint64_t &dslv_, 
				 double *x)
{
  dissection_solver_ptr *dslv = (dissection_solver_ptr *)dslv_;
  int n0;
  if (dslv->quad_fact) {
    switch(dslv->real_or_complex) {
    case DISSECTION_REAL_MATRIX:
      n0 = dslv->rqtr_fwbw->kern_dimension();
      if (n0 > 0) {
	dslv->rqtr_fwbw->ProjectionImageSingle(x);
      }
      break;
    case DISSECTION_COMPLEX_MATRIX:
      if (dslv->verbose > 0) {
	fprintf(dslv->fp, 
		"%s %d diss_project() for complex is not yet implemented\n",
		__FILE__, __LINE__);
      }
      break;
    default:
      if (dslv->verbose > 0) {
	fprintf(dslv->fp, "%s %d unknown matrix data type : %d\n", 
		__FILE__, __LINE__, dslv->real_or_complex);
      }
    }
  }
  else {
    switch(dslv->real_or_complex) {
    case DISSECTION_REAL_MATRIX:
      n0 = dslv->rptr->kern_dimension();
      if (n0 > 0) {
	dslv->rptr->ProjectionImageSingle(x);
      }
      break;
    case DISSECTION_COMPLEX_MATRIX:
      if (dslv->verbose > 0) {
	fprintf(dslv->fp, 
		"%s %d diss_project() for complex is not yet implemented\n",
		__FILE__, __LINE__);
      }
      break;
    default:
      if (dslv->verbose > 0) {
	fprintf(dslv->fp, "%s %d unknown matrix data type : %d\n", 
		__FILE__, __LINE__, dslv->real_or_complex);
      }
    }
  }
}

DISSECTION_API void DISS_SOLVE_1(uint64_t &dslv_, 
				 double *x, 
				 const int &projection,
				 const int &trans)
{
  const bool isProj = (bool)(projection == 1);
  const bool isTrans = (bool)(trans == 1);
  
  dissection_solver_ptr *dslv = (dissection_solver_ptr *)dslv_;
  if (dslv->quad_fact) {
    switch(dslv->real_or_complex) {
    case DISSECTION_REAL_MATRIX:
      dslv->rqtr_fwbw->SolveSingle(x, isProj, isTrans, true); 
      // isTrans
      break;
    case DISSECTION_COMPLEX_MATRIX:    
      fprintf(dslv->cptr->get_filedescriptor(), 
	      "Dissection::SolveSingle : %p\n", dslv->cptr);
      dslv->cqtr_fwbw->SolveSingle((complex<double> *)x,
				   isProj, isTrans, true); 
      // isTrans
      break;
    default:
      if (dslv->verbose > 0) {
	fprintf(dslv->fp, "%s %d unknown matrix data type : %d\n", 
		__FILE__, __LINE__, dslv->real_or_complex);
      }
    }

  }
  else {
    switch(dslv->real_or_complex) {
    case DISSECTION_REAL_MATRIX:
      dslv->rptr->SolveSingle(x, isProj, isTrans, true);
      break;
    case DISSECTION_COMPLEX_MATRIX:    
      fprintf(dslv->cptr->get_filedescriptor(), 
	      "Dissection::SolveSingle : %p\n", dslv->cptr);
      dslv->cptr->SolveSingle((complex<double> *)x, isProj, isTrans, true);
      // isTrans
      break;
    default:
      if (dslv->verbose > 0) {
	fprintf(dslv->fp, "%s %d unknown matrix data type : %d\n", 
		__FILE__, __LINE__, dslv->real_or_complex);
      }
    }
  }
}

DISSECTION_API void DISS_SOLVE_N(uint64_t &dslv_, 
				 double *x, 
				 const int &nrhs, const int &projection,
				 const int &trans)
{
  const bool isProj = (bool)(projection == 1);
  const bool isTrans = (bool)(trans == 1);
  
  dissection_solver_ptr *dslv = (dissection_solver_ptr *)dslv_;
  if (dslv->quad_fact) {
    switch(dslv->real_or_complex) {
    case DISSECTION_REAL_MATRIX:
      dslv->rqtr_fwbw->SolveMulti(x, nrhs, isProj, isTrans, true);
      break;
    case DISSECTION_COMPLEX_MATRIX:
      dslv->cqtr_fwbw->SolveMulti((complex<double> *)x, nrhs,
				  isProj, isTrans, true);
      break;
    default:
      if (dslv->verbose > 0) {
	fprintf(dslv->fp, "%s %d unknown matrix data type : %d\n", 
		__FILE__, __LINE__, dslv->real_or_complex);
      }
    }
  }
  else {
    switch(dslv->real_or_complex) {
    case DISSECTION_REAL_MATRIX:
      dslv->rptr->SolveMulti(x, nrhs, isProj, isTrans, true);
      break;
    case DISSECTION_COMPLEX_MATRIX:
      dslv->cptr->SolveMulti((complex<double> *)x, nrhs, isProj, isTrans, true);
      break;
    default:
      if (dslv->verbose > 0) {
	fprintf(dslv->fp, "%s %d unknown matrix data type : %d\n", 
		__FILE__, __LINE__, dslv->real_or_complex);
      }
    }
  }
}

DISSECTION_API void DISS_MATRIX_PRODUCT(uint64_t &dslv_,
					const double* x, double* y)
{
  dissection_solver_ptr *dslv = (dissection_solver_ptr *)dslv_;
  if (dslv->quad_fact) {
    switch(dslv->real_or_complex) {
    case DISSECTION_REAL_MATRIX:
      dslv->rqtr_fwbw->SpMV(x, y);
      break;
    case DISSECTION_COMPLEX_MATRIX:
      dslv->cqtr_fwbw->SpMV((complex<double> *)x, (complex<double> *)y);    
      break;
    default:
      if (dslv->verbose > 0) {
	fprintf(dslv->fp, "%s %d unknown matrix data type : %d\n", 
		__FILE__, __LINE__, dslv->real_or_complex);
      }
    }
  }
  else {
    switch(dslv->real_or_complex) {
    case DISSECTION_REAL_MATRIX:
      dslv->rptr->SpMV(x, y);
      break;
    case DISSECTION_COMPLEX_MATRIX:
      dslv->cptr->SpMV((complex<double> *)x, (complex<double> *)y);    
      break;
    default:
      if (dslv->verbose > 0) {
	fprintf(dslv->fp, "%s %d unknown matrix data type : %d\n", 
		__FILE__, __LINE__, dslv->real_or_complex);
      }
    }
  }
}
