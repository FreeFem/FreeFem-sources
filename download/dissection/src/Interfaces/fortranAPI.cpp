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
  DissectionSolver<double, double, quadruple, quadruple> *rqtr_fwbw;
  DissectionSolver<complex<quadruple>, quadruple> *cqtr;
  DissectionSolver<complex<double>, double,
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
