
// ORIG-DATE: 02/2009
// -*- Mode : c++ -*-
//
// SUMMARY  :
// USAGE    : LGPL
// ORG      : INRIA Saclay
// AUTHOR   : Guy Atenekeng
// E-MAIL   : Guy_Antoine_Atenekeng_Kahou@lri.fr
//
//ff-c++-LIBRARY-dep: hypre  metis  blas  mpi
//ff-c++-cpp-dep:

/*
 This file is part of Freefem++

 Freefem++ is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation; either version 2.1 of the License, or
 (at your option) any later version.

 Freefem++  is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public License
 along with Freefem++; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

 Thank to the ARN ()  FF2A3 grant
 ref:ANR-07-CIS7-002-01
 */

// FFCS: add requirement for MPI
//ff-c++-LIBRARY-dep: hypre mpi
//ff-c++-cpp-dep:

// add F.Hecht ...  oct 2010
#define HYPRE_TIMING
// .. end add
#include <mpi.h>

#include <iostream>
using namespace std;

#include "rgraph.hpp"
#include "error.hpp"
#include "AFunction.hpp"

//#include "lex.hpp"
#include "MatriceCreuse_tpl.hpp"

#ifdef __cplusplus
extern "C" {
#include "metis.h"
#endif
#ifdef __cplusplus
}
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef MPI_WTIME_IS_GLOBAL
#undef MPI_WTIME_IS_GLOBAL
#endif
#define MPI_WTIME_IS_GLOBAL 1
#define STATS

#include "_hypre_utilities.h"
#include "HYPRE.h"
#include "HYPRE_krylov.h"
#include "HYPRE_IJ_mv.h"
#include "HYPRE_parcsr_ls.h"
#include "_hypre_parcsr_mv.h"
#include "fortran_matrix.h"
#include "HYPRE_lobpcg.h"

#include "interpreter.h"
#include "multivector.h"
#include "HYPRE_MatvecFunctions.h"

#include "HYPRE_parcsr_int.h"

#define CHECK_ZERO
#define MCW MPI_COMM_WORLD
#define BUFLEN 100
#define SCALE 0

#ifdef SUN
/*
 * 	It uses the system call gethrtime(3C), which is accurate to
 *	nanoseconds.
 */
#include <sys/time.h>

double dwalltime( ) { return ((double)gethrtime( ) / 1e9); }

#else

#ifndef NO_TIMER
#include <sys/types.h>
#include <sys/times.h>
#include <time.h>
#include <sys/time.h>
#endif

#ifndef CLK_TCK
#define CLK_TCK 60
#endif

double dwalltime( ) {
#ifdef NO_TIMER
  /* no sys/times.h on _WIN32 */
  double tmp;
  tmp = 0.0;
#else
  struct tms use;
  double tmp;
  times(&use);
  tmp = use.tms_utime;
  tmp += use.tms_stime;
#endif
  return (double)(tmp) / CLK_TCK;
}

#endif

int roscal(int n, int job, int nrm, double *AAv, int *p, int *pr, double *scaletmpr, int *ierr) {
  /*---------------------------------------------------------------------
  |
  | This routine scales each row of mata so that the norm is 1.
  |
  |----------------------------------------------------------------------
  | on entry:
  | mata  = the matrix (in SparRow form)
  | nrm   = type of norm
  |          0 (\infty),  1 or 2
  |
  | on return
  | diag  = diag[j] = 1/norm(row[j])
  |
  |     0 --> normal return
  |     j --> row j is a zero row
  |--------------------------------------------------------------------*/
  /*   local variables    */
  int i, k;
  double scal;

  for (i = 0; i < n; i++) {
    scal = 0.0;
    // kr = &(AAv[pr[i]]);
    if (nrm == 0) {
      for (k = pr[i]; k < pr[i + 1]; k++)
        if (fabs(AAv[k]) > scal) scal = fabs(AAv[k]);
    } else if (nrm == 1) {
      for (k = pr[i]; k < pr[i + 1]; k++) scal += fabs(AAv[k]);
    } else { /* nrm = 2 */
      for (k = pr[i]; k < (pr[i + 1]); k++) scal += AAv[k] * AAv[k];
    }
    if (nrm == 2) scal = sqrt(scal);
    if (scal == 0.0) {
      *ierr = i;
      return i + 1;
    } else
      scal = 1.0 / scal;
    scaletmpr[i] = scal;
    for (k = pr[i]; k < (pr[i + 1]); k++) AAv[k] = AAv[k] * scal;
  }
  *ierr = 0;
  return 0;
}
/*---------------end of roscalC-----------------------------------------
----------------------------------------------------------------------*/
int coscal(int n, int job, int nrm, double *AAv, int *p, int *pr, double *scaletmpc, int *ierr) {
  /*---------------------------------------------------------------------
  |
  | This routine scales each column of mata so that the norm is 1.
  |
  |----------------------------------------------------------------------
  | on entry:
  | mata  = the matrix (in SparRow form)
  | nrm   = type of norm
  |          0 (\infty),  1 or 2
  |
  | on return
  | diag  = diag[j] = 1/norm(row[j])
  |
  |     0 --> normal return
  |     j --> column j is a zero column
  |--------------------------------------------------------------------*/
  /*   local variables    */
  int i, j, k;
  double *kr;
  int *ki;
  for (i = 0; i < n; i++) scaletmpc[i] = 0.0;
  /*---------------------------------------
  |   compute the norm of each column
  |--------------------------------------*/
  for (i = 0; i < n; i++) {
    kr = &(AAv[pr[i]]);
    ki = &(pr[i]);
    if (nrm == 0) {
      for (k = pr[i]; k < pr[i + 1]; k++) {
        j = pr[i];
        if (fabs(AAv[k]) > scaletmpc[p[k]]) scaletmpc[p[k]] = fabs(AAv[k]);
      }
    } else if (nrm == 1) {
      for (k = pr[i]; k < pr[i + 1]; k++) scaletmpc[p[k]] += fabs(AAv[k]);
    } else {
      for (k = pr[i]; k < pr[i + 1]; k++) scaletmpc[p[k]] += fabs(AAv[k]) * fabs(AAv[k]);
    }
  }
  if (nrm == 2) {
    for (i = 0; i < n; i++) scaletmpc[i] = sqrt(scaletmpc[i]);
  }
  /*---------------------------------------
  |   invert
  |--------------------------------------*/
  for (i = 0; i < n; i++) {
    if (scaletmpc[i] == 0.0) {
      *ierr = i + 1;
      return i + 1;
    } else
      scaletmpc[i] = 1.0 / scaletmpc[i];
  }
  /*---------------------------------------
  |   C = A * D
  |--------------------------------------*/
  for (i = 0; i < n; i++) {

    for (k = pr[i]; k < pr[i + 1]; k++) AAv[k] = AAv[k] * scaletmpc[p[k]];
  }
  *ierr = 0;
  return 0;
}
/*---------------end of coscalC-----------------------------------------
----------------------------------------------------------------------*/

typedef struct sparse_mat_loc {
  int *ptr;     // index of the beginning of rows in id_cols and vals
  int *rows;    // index of non empty rows
  int *id_cols;
  double *vals;
  int ilower;    // lower index limit
  int iupper;    // upper index limit
  int n_loc;     /*number of rows*/
  int *ncols;    /*number of columns in each row*/
} sparse_mat_loc;

/*
  **Function to distribute a sparse matrix as blocks rows on several processes**

  A: (input) sparse matrix.
    (matrix A as input is available on each process
  type :(input) 0=CSR format, any other value=CSC
  size: (input) size of the communicator
  rk: (input) rank of the process in the communicator
  A_loc: (output) sparse matrix in CSR reduced on local process
  */

int dist_matrix(int n, int *ptr, int *id_rows, double *vals, int type, int size, int rk,
                int *mapptr, int *maptmp, int *iwork1, sparse_mat_loc *A_loc) {
  int i, j, ilower, iupper;
  int n_loc, nnz_loc, relpos;
  int *marker;

  /* Compute the number of rows to distribute to local process */
  n_loc = n / size;
  if (rk == size - 1) n_loc = n_loc + n % size;

  int i1, i2, ncols;
  // save the lower index (ilower) and upper (iupper) for each process
  (A_loc)->ilower = rk * (n / size);
  ilower = (A_loc)->ilower;
  (A_loc)->iupper = (rk + 1) * (n / size) - 1;
  if (rk == size - 1) (A_loc)->iupper = (A_loc)->iupper + n % size;
  iupper = (A_loc)->iupper;
  (A_loc)->n_loc = n_loc;

  if (!((A_loc)->ptr = (int *)malloc((n_loc + 1) * sizeof(int)))) {
    printf("%s", "Malloc fails for ptr \n");
    exit(1);
  }
  if (!((A_loc)->rows = (int *)malloc((n_loc) * sizeof(int)))) {
    printf("%s", "Malloc fails for rows \n");
    exit(1);
  }
  if (!((A_loc)->ncols = (int *)malloc((n_loc) * sizeof(int)))) {
    printf("%s", "Malloc fails for ncols \n");
    exit(1);
  }

  // Change global Input matrix (A) to local (A_loc) on each process
  if (type == 0) {    // Matrix A is in CSR format

    // Gets local nnz
    i1 = (A_loc)->ilower;
    i2 = (A_loc)->iupper;
    (A_loc)->ilower = i1;
    (A_loc)->iupper = i2;
    nnz_loc = 0;
    for (i = i1; i <= i2; i++) {
      nnz_loc += ptr[i + 1] - ptr[i];
    }
    // Allocate memory for local matrix

    if (!((A_loc)->id_cols = (int *)malloc(nnz_loc * sizeof(int)))) {
      printf("%s", "Malloc fails for id_cols \n");
      exit(1);
    }
    if (!((A_loc)->vals = (double *)malloc(nnz_loc * sizeof(double)))) {
      printf("%s", "Malloc fails for vals");
      exit(1);
    }

    // Transfer the corresponding values from global to local
    relpos = 0;
    // int ncols; //count number of elements in each row
    for (i = i1; i <= i2; i++) {
      (A_loc)->rows[i - i1] = i;
      (A_loc)->ptr[i - i1] = relpos;
      ncols = relpos;
      for (j = ptr[i]; j < ptr[i + 1]; j++) {
        (A_loc)->id_cols[relpos] = id_rows[j];
        (A_loc)->vals[relpos] = vals[j];
        relpos++;
      }
      (A_loc)->ncols[i - i1] = relpos - ncols;
    }
    // cout << "taille des sous domaines" << nnz_loc << endl;
  } else {                                     // matrix A is in CSC format
    marker = (int *)calloc(n, sizeof(int));    // count number of elements in each row
    for (i = 0; i < n; i++)
      for (j = ptr[i]; j < ptr[i + 1]; j++) marker[id_rows[j]]++;

    (A_loc)->ptr[0] = 0;    // set up the beginning of each row
    for (i = 0; i < n_loc; i++) {
      (A_loc)->ptr[i + 1] = (A_loc)->ptr[i] + marker[i + ilower];
      (A_loc)->id_cols[relpos] = id_rows[j];
      (A_loc)->vals[relpos] = vals[j];
      relpos++;
    }
    (A_loc)->ncols[i - ilower] = relpos - ncols;
  }
  return 0;
}

class hypreParam {
  // Solveur and preconditionner
 public:
  char solver[BUFLEN];
  char precon[BUFLEN];
  // BoomerAMG parameter
  int amg_coarsentype; /* Falgout coarsening */
  int amg_relaxtype;   /* hybrid Gauss-Seidel or SOR */
  int amg_interptype;  /* default*/
  int amg_maxlevels;
  int amg_numsweeps;          /*default*/
  double amg_strongthreshold; /*suitable for 3D Laplace Operator*/
  double amg_truncfactor;
  int amg_prntlevel; /* print setup info */
  double amg_tol;    // BoomerAMG Tolerance
  int amg_maxiter;
  int scale;
  int gsmg_sample, measure_type, cycle_type;
  int solv_stopcrit;
  double trunc_factor;

  // More complex smoothers (Schwarz methods, Pilut, Parasails, Euclid)
  int smooth_type;
  int smooth_numlevels;
  int smooth_numsweeps;
  double pilut_droptol;
  double pilut_maxnz;
  int schwarz_overlap;
  int schwarz_variant;
  int schwarz_domaintype;

  // parasails parameter
  int sai_max_levels;
  double sai_threshold;

  double sai_filter;
  int sai_sym;
  int sai_log;
  int VERBOSE;

  /***************************HYPRE_BOOMERAMG***********************/
  double strong_threshold;

  int *num_grid_sweeps;
  int *grid_relax_type;
  int *grid_relax_points;
  double *relax_weight;
  double *omega;
  // Solver parameter (used for GMRES , PCG or BiCGStab)
  double solv_tol;
  int solv_maxiter;
  int solv_kdim;
  int solv_log;
  int solv_prntlevel;
  int precond_id, solver_id, matrix_id, solver_type;
  int smooth_num_levels, smooth_num_sweeps, max_levels, Two_Norm;
  int domain_type, num_functions, variant, overlap, nonzeros_to_keep;
  double max_row_sum, drop_tol;
  int *dof_func;
  int pcg_max_its, rrow_size, Rel_change;
  int dscg_max_its, coarsen_type, hybrid, k_dim, num_sweep;
  int pmax_iter;
  double cf_tol, tol, pc_tol;
  double schwarz_rlx_weight;

  /*For timing*/
  int timing;

 public:
  hypreParam(const KN< long > &param_int, const KN< double > &param_double) {
    amg_coarsentype = 6; /* Falgout coarsening */
    amg_relaxtype = 3;   /* hybrid Gauss-Seidel or SOR */
    amg_interptype = 0;  /* default*/
    amg_maxlevels = 25;
    amg_numsweeps = 1;          /*default*/
    amg_strongthreshold = 0.25; /*suitable for 3D Laplace Operator*/
    amg_truncfactor = 0.3;
    amg_prntlevel = 1; /* print setup info */
    amg_tol = 0.0;     // BoomerAMG Tolerance
    amg_maxiter = 20;
    gsmg_sample = 1;
    // More complex smoothers (Schwarz methods, Pilut, Parasails, Euclid)
    smooth_type = 6;
    smooth_numlevels = 3;
    smooth_numsweeps = 1;
    pilut_droptol = 1.0e-4;
    pilut_maxnz = 100;
    schwarz_overlap = 10;
    schwarz_variant = 2;
    schwarz_domaintype = 2;
    // parasails parameter
    sai_max_levels = 1;
    sai_threshold = 0.1;
    sai_filter = 0.1;
    sai_sym = 0;
    sai_log = 1;
    int i;

    // Solver parameter (used for GMRES or BiCGStab)
    solv_tol = 1.0e-11;
    solv_maxiter = 1000;
    solv_kdim = 40;

    solv_log = 0;
    solv_prntlevel = 2;
    precond_id = 0;    // BOOMER AMG
    solver_id = 1;     // GMRES as solver
    VERBOSE = 0;
    scale = 1;
    pmax_iter = 30;
    rrow_size = 1000;
    solv_stopcrit = 1;

    amg_interptype = 6;
    gsmg_sample = 5;
    amg_coarsentype = 6;

    measure_type = 1;
    amg_strongthreshold = 0.25;
    trunc_factor = 1e-2;

    amg_maxiter = 20;
    cycle_type = 1;
    smooth_num_levels = 3;
    smooth_num_sweeps = 1;
    max_levels = 25;
    hybrid = 1;
    k_dim = 5;
    smooth_type = 6;
    num_functions = 1;
    smooth_num_levels = 3;
    smooth_num_sweeps = 2;
    num_sweep = 1;
    variant = 0;
    overlap = 10;
    domain_type = 2;
    nonzeros_to_keep = 1000;
    tol = 1.e-8;
    pc_tol = 0.;
    drop_tol = -1.;
    max_row_sum = 0.9;
    schwarz_rlx_weight = 1.;
    sai_threshold = 0.1;
    sai_filter = 0.1;

    relax_weight = hypre_CTAlloc(double, max_levels);
    omega = hypre_CTAlloc(double, max_levels);
    for (i = 0; i < max_levels; i++) {
      relax_weight[i] = 1.;
      omega[i] = 1.;
    }
    max_row_sum = 0.9;
    schwarz_rlx_weight = 1.;
    variant = 0;
    num_functions = 1;
    overlap = 10;
    domain_type = 2;

    if (param_int.N( ) > 0) {
      if ((param_int[0] >= 0) && (param_int[0] <= 9))
        solver_id = param_int[0];
      else
        solver_id = 1;
    }    // GMRES as solver

    if ((solver_id != 4) && (solver_id != 5)) {
      if (param_int.N( ) > 1) {
        if ((param_int[1] >= 0) && (param_int[1] <= 9))
          precond_id = param_int[1];
        else
          precond_id = 0;
      }    // BOOMER AMG }
      if (param_int.N( ) > 2) {
        if (param_int[2] > 0)
          solv_maxiter = param_int[2];
        else
          solv_maxiter = 1000;
      }
      if (param_int.N( ) > 3) {
        if (param_int[3] > 0)
          solv_kdim = param_int[3];
        else
          solv_kdim = 40;
      }
      if (param_int.N( ) > 4) {
        if (param_int[4] >= 0) solv_prntlevel = param_int[4];
      }
      if (param_int.N( ) > 5) {
        if (param_int[5] >= 0) solv_log = param_int[5];
      }
      if (param_int.N( ) > 6) {
        if (param_int[6] >= 0) solv_stopcrit = param_int[6];
      }

      if (param_double.N( ) > 0) {
        if (param_double[0] > 0) solv_tol = param_double[0];
      }
      switch (precond_id) {
        case 0:    // Preconditionner AMG
          if (param_int.N( ) > 7) {
            if (param_int[7] >= 0) amg_interptype = param_int[7];
          }
          if (param_int.N( ) > 8) {
            if (param_int[8] >= 0) gsmg_sample = param_int[8];
          }
          if (param_int.N( ) > 9) {
            if (param_int[9] >= 0) amg_coarsentype = param_int[9];
          }
          if (param_int.N( ) > 10) {
            if (param_int[10] >= 0) measure_type = param_int[10];
          }
          if (param_double.N( ) > 1) {
            if (param_double[1] > 0) amg_strongthreshold = param_double[1];
          }
          if (param_double.N( ) > 2) {
            if (param_double[2] > 0) trunc_factor = param_double[2];
          }
          // if(param_int.N()>11)  {if(param_int[11]>0) amg_maxiter=param_int[11];}
          if (param_int.N( ) > 11) {
            if (param_int[11] > 0) cycle_type = param_int[11];
          }
          if (param_int.N( ) > 12) {
            if (param_int[12] > 0) smooth_type = param_int[12];
          }
          if (param_int.N( ) > 13) {
            if (param_int[13] > 0) smooth_num_levels = param_int[13];
          }
          if (param_int.N( ) > 14) {
            if (param_int[14] > 0) smooth_num_sweeps = param_int[14];
          }
          if (param_int.N( ) > 15) {
            if (param_int[15] > 0) max_levels = param_int[15];
          }
          relax_weight = hypre_CTAlloc(double, max_levels);
          omega = hypre_CTAlloc(double, max_levels);
          for (i = 0; i < max_levels; i++) {
            relax_weight[i] = 1.;
            omega[i] = 1.;
          }
          if (param_double.N( ) > 3) {
            if (param_double[3] > 0)
              max_row_sum = param_double[3];
            else
              max_row_sum = 0.9;
          } else
            max_row_sum = 0.9;
          if (param_double.N( ) > 4) {
            if (param_double[4] > 0)
              schwarz_rlx_weight = param_double[4];
            else
              schwarz_rlx_weight = 1.;
          } else
            schwarz_rlx_weight = 1.;
          if (param_int.N( ) > 16) {
            if (param_int[16] > 0)
              variant = param_int[16];
            else
              variant = 3;
          } else
            variant = 3;
          if (param_int.N( ) > 17) {
            if (param_int[17] > 0)
              num_functions = param_int[17];
            else
              num_functions = 1;
          } else
            num_functions = 1;
          if (param_int.N( ) > 18) {
            if (param_int[18] > 0)
              overlap = param_int[18];
            else
              overlap = 10;
          } else
            overlap = 10;
          if (param_int.N( ) > 19) {
            if (param_int[19] > 0)
              domain_type = param_int[19];
            else
              domain_type = 2;
          } else
            domain_type = 2;
          break;

        case 1:    // Preconditionner PILUT

          if (param_double.N( ) > 1) {
            if (param_double[1] > 0)
              drop_tol = param_double[1];
            else
              drop_tol = 1e-5;
          } else
            drop_tol = 1e-5;
          if (param_int.N( ) > 7) {
            if (param_int[7] > 0)
              nonzeros_to_keep = param_int[7];
            else
              nonzeros_to_keep = 1000;
          } else
            nonzeros_to_keep = 1000;
          if (param_int.N( ) > 8) {
            if (param_int[8] > 0) pmax_iter = param_int[8];
          }
          if (param_int.N( ) > 9) {
            if (param_int[9] > 0) rrow_size = param_int[9];
          }

          break;

        case 2:    // Preconditionner ParaSails
          if (param_double.N( ) > 1) {
            if (param_double[1] > 0)
              sai_filter = param_double[1];
            else
              sai_filter = 0.1;
          } else
            sai_filter = 0.1;
          if (param_double.N( ) > 2) {
            if (param_double[2] > 0)
              sai_threshold = param_double[2];
            else
              sai_threshold = 0.1;
          } else
            sai_threshold = 0.1;
          if (param_int.N( ) > 7) {
            if (param_int[7] > 0)
              max_levels = param_int[7];
            else
              max_levels = 1;
          } else
            max_levels = 1;
          if (param_int.N( ) > 8) {
            if (param_int[8] > 0)
              sai_sym = param_int[8];
            else
              sai_sym = 0;
          } else
            sai_sym = 0;
          break;

        case 3:    // Preconditionner Schwarz
          if (param_double.N( ) > 1) {
            if (param_double[1] > 0) schwarz_rlx_weight = param_double[1];
            schwarz_rlx_weight = 1;
          } else
            schwarz_rlx_weight = 1.;
          if (param_int.N( ) > 7) {
            if (param_int[7] > 0)
              variant = param_int[7];
            else
              variant = 1;
          } else
            variant = 1;
          if (param_int.N( ) > 8) {
            if (param_int[8] > 0)
              overlap = param_int[8];
            else
              overlap = 1;
          } else
            overlap = 1;
          if (param_int.N( ) > 9) {
            if (param_int[9] > 0)
              domain_type = param_int[9];
            else
              domain_type = 3;
          } else
            domain_type = 3;
          break;

        default:
          break;
      }
    }

    if (solver_id == 4)    // Solver AMG_HYBRID
    {

      if (param_double.N( ) > 1) {
        if (param_double[1] >= 0)
          amg_tol = param_double[1];
        else
          amg_tol = 1e-9;
      } else
        amg_tol = 1e-9;
      if (param_double.N( ) > 2) {
        if (param_double[2] >= 0)
          cf_tol = param_double[2];
        else
          cf_tol = 1e-3;
      } else
        cf_tol = 1e-3;
      if (param_int.N( ) > 1) {
        if (param_int[1] >= 0)
          solver_type = param_int[1];
        else
          solver_type = 1;
      } else
        solver_type = 1;
      if (param_int.N( ) > 2) {
        if (param_int[2] > 0)
          dscg_max_its = param_int[2];
        else
          dscg_max_its = 1000;
      } else
        dscg_max_its = 1000;
      if (param_int.N( ) > 3) {
        if (param_int[3] > 0)
          pcg_max_its = param_int[3];
        else
          pcg_max_its = 200;
      } else
        pcg_max_its = 200;
      if (param_int.N( ) > 4) {
        if (param_int[4] > 0)
          coarsen_type = param_int[4];
        else
          coarsen_type = 6;
      } else
        coarsen_type = 6;
      if (param_double.N( ) > 3) {
        if (param_double[3] > 0)
          strong_threshold = param_double[3];
        else
          strong_threshold = 1e-3;
      } else
        strong_threshold = 1e-3;
      if (param_double.N( ) > 4) {
        if (param_double[4] > 0)
          trunc_factor = param_double[4];
        else
          trunc_factor = 1e-2;
      } else
        trunc_factor = 1e-2;
      if (param_int.N( ) > 5) {
        if (param_int[5] > 0)
          max_levels = param_int[5];
        else
          max_levels = 25;
      } else
        max_levels = 25;
      if (param_double.N( ) > 5) {
        if (param_double[5] > 0)
          max_row_sum = param_double[5];
        else
          max_row_sum = 0.9;
      } else
        max_row_sum = 0.9;

      relax_weight = hypre_CTAlloc(double, max_levels);
      omega = hypre_CTAlloc(double, max_levels);
      for (i = 0; i < max_levels; i++) {
        relax_weight[i] = 1.;
        omega[i] = 1.;
      }
    }
    if (solver_id == 3) {
      if (param_int.N( ) > 7) {
        if (param_int[7] >= 0)
          Two_Norm = param_int[7];
        else
          Two_Norm = 2;
      } else
        Two_Norm = 2;
      if (param_int.N( ) > 8) {
        if (param_int[8] >= 0)
          Rel_change = param_int[8];
        else
          Rel_change = 1;
      } else
        Rel_change = 1;
    }
    if (solver_id == 5)    // Solver AMG
    {
      if (param_int.N( ) > 7) {
        if (param_int[7] >= 0) amg_interptype = param_int[7];
      }
      if (param_int.N( ) > 8) {
        if (param_int[8] >= 0) gsmg_sample = param_int[8];
      }
      if (param_int.N( ) > 9) {
        if (param_int[9] >= 0) amg_coarsentype = param_int[9];
      }
      if (param_int.N( ) > 10) {
        if (param_int[10] >= 0) measure_type = param_int[10];
      }
      if (param_double.N( ) > 2) {
        if (param_double[2] > 0) amg_strongthreshold = param_double[2];
      }
      if (param_double.N( ) > 3) {
        if (param_double[3] > 0) trunc_factor = param_double[3];
      }
      if (param_int.N( ) > 11) {
        if (param_int[11] > 0) amg_maxiter = param_int[11];
      }
      if (param_int.N( ) > 12) {
        if (param_int[12] > 0) cycle_type = param_int[12];
      }
      if (param_int.N( ) > 13) {
        if (param_int[13] > 0) smooth_type = param_int[13];
      }
      if (param_int.N( ) > 14) {
        if (param_int[14] > 0) smooth_num_levels = param_int[14];
      }
      if (param_int.N( ) > 15) {
        if (param_int[15] > 0) smooth_num_sweeps = param_int[15];
      }
      if (param_int.N( ) > 16) {
        if (param_int[16] > 0) max_levels = param_int[16];
      }
      relax_weight = hypre_CTAlloc(double, max_levels);
      omega = hypre_CTAlloc(double, max_levels);
      for (i = 0; i < max_levels; i++) {
        relax_weight[i] = 1.;
        omega[i] = 1.;
      }
      if (param_double.N( ) > 4) {
        if (param_double[4] > 0)
          max_row_sum = param_double[4];
        else
          max_row_sum = 1e-1;
      } else
        max_row_sum = 1e-1;
      if (param_double.N( ) > 5) {
        if (param_double[5] > 0)
          schwarz_rlx_weight = param_double[5];
        else
          schwarz_rlx_weight = 1.;
      } else
        schwarz_rlx_weight = 1.;
      if (param_int.N( ) > 17) {
        if (param_int[17] > 0)
          variant = param_int[17];
        else
          variant = 1;
      } else
        variant = 1;
      if (param_int.N( ) > 19) {
        if (param_int[18] > 0)
          num_functions = param_int[18];
        else
          num_functions = 5;
      } else
        num_functions = 5;
      if (param_int.N( ) > 20) {
        if (param_int[19] > 0)
          overlap = param_int[19];
        else
          overlap = 1;
      } else
        overlap = 1;
      if (param_int.N( ) > 21) {
        if (param_int[20] > 0)
          domain_type = param_int[20];
        else
          domain_type = 1;
      } else
        domain_type = 1;
    }
    if (param_int.N( ) > 22) {
      if (param_int[22] > 0) VERBOSE = param_int[22];
    }
    if (param_int.N( ) > 23) {
      if (param_int[23] > 0)
        scale = param_int[23];
      else
        scale = 1;
    } else
      scale = 1;
    if (param_int.N( ) > 24) {
      if (param_int[24] > 0)
        timing = param_int[24];
      else
        timing = 1;
    } else
      timing = 1;
  }

 public:
  hypreParam( ) {
    int i;
    amg_coarsentype = 6; /* Falgout coarsening */
    amg_relaxtype = 3;   /* hybrid Gauss-Seidel or SOR */
    amg_interptype = 0;  /* default*/
    amg_maxlevels = 25;
    amg_numsweeps = 1;          /*default*/
    amg_strongthreshold = 0.25; /*suitable for 3D Laplace Operator*/
    amg_truncfactor = 0.3;
    amg_prntlevel = 1; /* print setup info */
    amg_tol = 1e-7;    // BoomerAMG Tolerance
    amg_maxiter = 1;
    gsmg_sample = 1;

    // Solver parameter (used for GMRES or BiCGStab)
    solv_tol = 1.0e-11;
    solv_maxiter = 1000;
    solv_kdim = 40;

    solv_log = 0;
    solv_prntlevel = 2;
    precond_id = 0;    // BOOMER AMG
    solver_id = 1;     // GMRES as solver
    VERBOSE = 0;
    scale = 1;
    pmax_iter = 30;
    rrow_size = 1000;

    amg_interptype = 0;
    gsmg_sample = 1;
    amg_coarsentype = 6;
    measure_type = 1;
    amg_strongthreshold = 0.25;
    trunc_factor = 1e-2;
    amg_maxiter = 20;
    cycle_type = 1;
    smooth_type = 6;
    smooth_num_levels = 0;
    smooth_num_sweeps = 2;
    max_levels = 25;

    relax_weight = hypre_CTAlloc(double, max_levels);
    omega = hypre_CTAlloc(double, max_levels);
    for (i = 0; i < max_levels; i++) {
      relax_weight[i] = 1.;
      omega[i] = 1.;
    }
    max_row_sum = 0.9;
    schwarz_rlx_weight = 1.;
    variant = 0;
    num_functions = 1;
    overlap = 10;
    domain_type = 0;
  }

 public:
  hypreParam(char *fileparameter, MPI_Comm comm) {
    FILE *f;
    char buf[BUFLEN];
    int num;
    int rk, size;
    MPI_Comm_rank(comm, &rk);
    MPI_Comm_size(comm, &size);

    amg_coarsentype = 6; /* Falgout coarsening */
    amg_relaxtype = 3;   /* hybrid Gauss-Seidel or SOR */
    amg_interptype = 0;  /* default*/
    amg_maxlevels = 25;
    amg_numsweeps = 1;          /*default*/
    amg_strongthreshold = 0.25; /*suitable for 3D Laplace Operator*/
    amg_truncfactor = 0.3;
    amg_prntlevel = 1; /* print setup info */
    amg_tol = 0.0;     // BoomerAMG Tolerance
    amg_maxiter = 20;
    // More complex smoothers (Schwarz methods, Pilut, Parasails, Euclid)
    smooth_type = 6;
    smooth_numlevels = 3;
    smooth_numsweeps = 1;
    pilut_droptol = 1.0e-4;
    pilut_maxnz = 100;
    schwarz_overlap = 10;
    schwarz_variant = 2;
    schwarz_domaintype = 2;
    // parasails parameter
    sai_max_levels = 1;
    sai_threshold = 0.1;
    sai_filter = 0.05;
    sai_sym = 0;
    sai_log = 1;

    // Solver parameter (used for GMRES or BiCGStab)
    solv_tol = 1.0e-30;
    solv_maxiter = 80;
    solv_kdim = 40;
    //  int                  solv_stopcrit = 1; //only for BiCGSTAB
    solv_log = 0;
    solv_prntlevel = 0;
    precond_id = 0;    // BOOMER AMG
    solver_id = 1;     // GMRES as solver
    VERBOSE = 0;
    scale = 1;

    if (fileparameter == NULL) {
      if (rk == 0)
        printf("%s", "Set default parameter because you not precise the parameter file \n \n");
      solver_id = 1;     // GMRES as solver
      precond_id = 0;    // BOOMER AMG

    } else if ((f = fopen(fileparameter, "r")) == NULL) {
      if (rk == 0)
        printf("%s", "Set default parameter because your parameter file not exist \n \n");
      solver_id = 1;     // GMRES as solver
      precond_id = 0;    // BOOMER AMG
    } else {
      if (rk == 0) printf("%s%s%s", "Read parameter from file ", fileparameter, "\n \n");
      num = 0;
      while (fgets(buf, BUFLEN, f) != NULL) {
        switch (num) {
          case 0:
            sscanf(buf, "%s", solver);
            break;
          case 1:
            sscanf(buf, "%d", &solver_id);
            break;
          case 2:
            sscanf(buf, "%lf", &solv_tol);
            break;
          case 3:
            sscanf(buf, "%d", &solv_maxiter);
            break;
          case 4:
            sscanf(buf, "%d", &solv_prntlevel);
            break;
          case 5:
            sscanf(buf, "%d", &solv_log);
            if (solver_id != 1) {
              fgets(buf, BUFLEN, f);
              num++;
            }
            break;
          case 6:
            sscanf(buf, "%d", &solv_kdim);
            break;
          case 7:
            sscanf(buf, "%s", precon);
            break;
          case 8:
            sscanf(buf, "%d", &precond_id);

            if (precond_id == 2)    // The parameter of preconditionner is inside file
            {
              fclose(f);
              exit(1);
            }
            break;
          case 9:
            if (precond_id == 0) sscanf(buf, "%d", &amg_coarsentype);
            if (precond_id == 1) sscanf(buf, "%lf", &sai_threshold);
            break;
          case 10:
            if (precond_id == 0) sscanf(buf, "%d", &amg_prntlevel);
            if (precond_id == 1) sscanf(buf, "%d", &sai_max_levels);
            break;
          case 11:
            if (precond_id == 0) sscanf(buf, "%d", &amg_interptype);
            if (precond_id == 1) sscanf(buf, "%lf", &sai_filter);
            break;
          case 12:
            if (precond_id == 0) sscanf(buf, "%d", &amg_maxlevels);
            if (precond_id == 1) case 17:
              sscanf(buf, "%d", &amg_prntlevel);
            break;
          case 18:
            sscanf(buf, "%lf", &amg_tol);
            break;
          case 19:
            sscanf(buf, "%d", &amg_maxiter);
            break;
          case 20:
            sscanf(buf, "%d", &scale);
            break;
          default:
            break;
        }
        num++;
      }
      if (fclose(f) == EOF)
        printf("%s", "Error while closing the file \n");
      else
        printf("%s", "File is well close \n");
    }
  }
};

class HypreSolver : public MatriceMorse< double >::VirtualSolver {
  mutable HYPRE_IJMatrix ij_A;
  mutable HYPRE_IJVector ij_B;
  mutable HYPRE_IJVector ij_X;
  mutable HYPRE_ParCSRMatrix par_A;
  mutable HYPRE_ParVector par_B;
  mutable HYPRE_ParVector par_X;
  void *object;
  HYPRE_ParCSRMatrix parcsr_A;

  string data_option;
  mutable int time_index, time_index1;
  mutable HYPRE_Solver solver;
  mutable HYPRE_Solver precond;
  mutable hypreParam *param;
  mutable double *scaletmpr, *scaletmpc;
  mutable int rk, size;
  int jlower, jupper;
  int ilower, iupper;
  mutable int num_iter;
  mutable double final_res_norm;
  mutable double tgv, eps, tol_pivot, tol_pivot_sym, epsr;
  sparse_mat_loc A_loc;
  mutable int *iwork, *maptmp, *mapptr, *iwork1, *riord;
  mutable int n_loc, n, VERBOSE;
  mutable MPI_Comm comm;
  mutable int pcg_num_its, dscg_num_its;

 public:
  HypreSolver(const MatriceMorse< double > &AA, string datafile, KN< long > &param_int,
              KN< double > &param_double, MPI_Comm *mpicommw) {
    int i, j, ierrr;
    if (mpicommw == NULL) {
      comm = MPI_COMM_WORLD;
    } else
      comm = *mpicommw;

    /*****INITIALIZE MPI ENVIRONMENT*****/

    ierrr = MPI_Comm_rank(comm, &rk);
    ierrr = MPI_Comm_size(comm, &size);

    if (((param_double != NULL) || (param_int != NULL)) && (datafile.empty( ))) {
      if (rk == 0) cout << "#########PARAMETERS ARE SET INSIDE A VECTOR###" << endl;
      param = new hypreParam(param_int, param_double);
    }
    if ((datafile.empty( )) && ((param_int == NULL) && (param_double == NULL))) {
      if (rk == 0) cout << "###########DEFAULT PARAMETERS WILL BE SET#######" << endl;
      param = new hypreParam( );
    }
    if ((!datafile.empty( )) && ((param_int == NULL) && (param_double == NULL))) {
      if (rk == 0) cout << "#########PARAMETERS ARE INSIDE A FILE#########" << endl;
      char *p;
      p = new char[datafile.length( ) + 1];
      strcpy(p, datafile.c_str( ));
      param = new hypreParam(p, comm);
    }

    /*###################################################
                         USING HYPRE
    ####################################################*/
    int n, nnz, *pr, *p, ierr;
    double *AAv;
    n = AA.n;
    nnz = AA.nbcoef;
    /*##################################################
                      COPY ENTRY MATRIX
     ##################################################*/
    if (param->timing) {
      time_index = hypre_InitializeTiming("MATRIX DISTRIBUTION");
      hypre_BeginTiming(time_index);
    }

    pr = new int[n + 1];
    p = new int[nnz];
    AAv = new double[nnz];
    for (i = 0; i < nnz; i++) {
      AAv[i] = AA.a[i];
      p[i] = AA.cl[i];
      if (i <= n) pr[i] = AA.lg[i];
    }
    int job, tmp;
    if (param->scale) {
      job = 1; /*-- compute  1-norm */
      tmp = 2; /*-- compute 2-norm  */
      scaletmpr = new double[n];
      scaletmpc = new double[n];
      roscal(n, job, tmp, AAv, p, pr, scaletmpr, &ierr);
      if (ierr) fprintf(stderr, "ERROR: IN ROSCAL, IERR = %d\n", ierr);
      /*------- scale the RHS according to row scaling coefficients */
      coscal(n, job, tmp, AAv, p, pr, scaletmpc, &ierr);
      if (ierr) fprintf(stderr, "ERROR : IN COSCAL, IERR = %d\n", ierr);
    } /*--- end of branch on scaling */

    /*************************************************************
    Distribute input matrix into local structures
    *************************************************************/
    int type = 0, wgtflag = 0, numflag = 0, volume;    // 0=CSR; 1=CSC
    int option[5];
    option[0] = 0;
    riord = (int *)malloc(sizeof(int) * n);
    if (riord == NULL) {
      if (rk == 0) printf("%s", "IN PARTITION , MEMORY ALLOCATION FAILED \n");
      exit(1);
    }
    /************************USE METIS FOR DATA DISTRIBUTION**************************/
    if (size > 1)
      METIS_PartGraphKway(&n, pr, p, NULL, NULL, &wgtflag, &numflag, &size, option, &volume, riord);
    else if (size == 1) {
      for (i = 0; i < n; i++) riord[i] = 0;
    }
    iwork = new int[n];
    maptmp = new int[n];
    mapptr = new int[size + 1];
    iwork1 = new int[size + 1];
    for (i = 0; i < size; i++) {
      iwork[i] = 0;
      iwork1[i] = 0;
    }
    for (j = 0; j < n; j++) {
      iwork[riord[j]]++;
      iwork1[riord[j]]++;
    }
    numflag = 0;
    for (i = 0; i < size; i++) {
      mapptr[i] = numflag;
      numflag += iwork[i];
    }
    mapptr[size] = numflag;

    for (i = 0; i < size; i++) {
      iwork[i] = mapptr[i];
    }
    for (i = 0; i < n; i++) {
      maptmp[iwork[riord[i]]] = i;
      iwork[riord[i]]++;
    }

    dist_matrix(AA.n, pr, p, AAv, type, size, rk, mapptr, maptmp, iwork1, &A_loc);
    /***Distribute vector***/
    n_loc = A_loc.n_loc;
    /**** Preparing Matrix *****/
    ierr = 0;
    ilower = A_loc.ilower;
    jlower = ilower;
    iupper = A_loc.iupper;
    jupper = iupper;
    ierr += HYPRE_IJMatrixCreate(comm, ilower, iupper, jlower, jupper, &ij_A);

    if (ierr) {
      printf("Error in driver building IJMatrix from parcsr matrix. \n");
      exit(-1);
    }

    ierr += HYPRE_IJMatrixSetObjectType(ij_A, HYPRE_PARCSR);
    ierr += HYPRE_IJMatrixInitialize(ij_A);

    int taille;

    for (i = ilower; i <= iupper; i++) {
      taille = pr[i + 1] - pr[i];

      ierr += HYPRE_IJMatrixSetValues(ij_A, 1, &taille, &i, &(p[pr[i]]), &(AAv[pr[i]]));
    }

    ierr += HYPRE_IJMatrixAssemble(ij_A);
    /*Extract HYPRE Object data structures suitable for the solvers*/
    ierr += HYPRE_IJMatrixGetObject(ij_A, (void **)&par_A);

    /*Free local matrix in A_loc, B and X*/
    free(A_loc.rows);
    free(A_loc.id_cols);
    free(A_loc.vals);
    free(A_loc.ncols);
    if (param->timing) {
      hypre_EndTiming(time_index);
      hypre_PrintTiming("IJ Matrix Setup", comm);
      hypre_FinalizeTiming(time_index);
      hypre_ClearTiming( );
    }

    /*************************************************************
                Create preconditioner

                  Setup and Use solver
   **************************************************************/

    param->timing = 1;

    switch (param->solver_id) {
      case 0:    // BiCGStab solver
        HYPRE_ParCSRBiCGSTABCreate(comm, &solver);
        HYPRE_ParCSRBiCGSTABSetTol(solver, param->solv_tol);
        HYPRE_ParCSRBiCGSTABSetMaxIter(solver, param->solv_maxiter);
        // HYPRE_ParCSRBiCGSTABSetStopCrit(solver, solv_stopcrit);
        HYPRE_ParCSRBiCGSTABSetPrintLevel(solver, param->solv_prntlevel);
        HYPRE_ParCSRBiCGSTABSetLogging(solver, param->solv_log);

        // Set Preconditioner
        switch (param->precond_id) {
          case 0:    // set BoomerAMG as preconditioner
            if (rk == 0) printf("SOLVER: BOOMERAMG-BiCGSTAB\n");
            HYPRE_BoomerAMGCreate(&precond);

            HYPRE_BoomerAMGSetInterpType(precond, param->amg_interptype);

            HYPRE_BoomerAMGSetTol(precond, param->amg_tol);
            HYPRE_BoomerAMGSetCoarsenType(precond, param->amg_coarsentype);
            HYPRE_BoomerAMGSetMeasureType(precond, param->measure_type);
            HYPRE_BoomerAMGSetStrongThreshold(precond, param->amg_strongthreshold);

            HYPRE_BoomerAMGSetTruncFactor(precond, param->trunc_factor);

            HYPRE_BoomerAMGSetMaxIter(precond, param->amg_maxiter);
            HYPRE_BoomerAMGSetCycleType(precond, param->cycle_type);

            HYPRE_BoomerAMGSetSmoothType(precond, param->smooth_type);
            HYPRE_BoomerAMGSetSmoothNumLevels(precond, param->smooth_num_levels);
            HYPRE_BoomerAMGSetSmoothNumSweeps(precond, param->smooth_num_sweeps);
            HYPRE_BoomerAMGSetMaxLevels(precond, param->max_levels);
            HYPRE_BoomerAMGSetMaxRowSum(precond, param->max_row_sum);

            HYPRE_BoomerAMGSetOverlap(precond, param->overlap);
            HYPRE_BoomerAMGSetVariant(precond, param->variant);
            HYPRE_BoomerAMGSetDomainType(precond, param->domain_type);

            HYPRE_BiCGSTABSetPrecond(solver, (HYPRE_PtrToSolverFcn)HYPRE_BoomerAMGSolve,
                                     (HYPRE_PtrToSolverFcn)HYPRE_BoomerAMGSetup, precond);
            break;
          case 4:
            /*use diagonal scaling as preconditioner*/
            if (rk == 0) printf("SOLVER: DS-BiCGSTAB\n");
            precond = NULL;
            HYPRE_ParCSRBiCGSTABSetPrecond(solver, (HYPRE_PtrToParSolverFcn)HYPRE_ParCSRDiagScale,
                                           (HYPRE_PtrToParSolverFcn)HYPRE_ParCSRDiagScaleSetup,
                                           precond);
            break;
          case 1:
            /*Use PILUT as preconditioner*/
            if (rk == 0) printf("SOLVER: PILUT-BiCGSTAB\n");
            ierr = HYPRE_ParCSRPilutCreate(comm, &precond);
            if (ierr) printf("ERROR: PILUT-BiCGSTAB\n");

            if (param->drop_tol >= 0) HYPRE_ParCSRPilutSetDropTolerance(precond, param->drop_tol);
            if (param->nonzeros_to_keep >= 0)
              HYPRE_ParCSRPilutSetFactorRowSize(precond, param->nonzeros_to_keep);

            HYPRE_ParCSRBiCGSTABSetPrecond(solver, (HYPRE_PtrToParSolverFcn)HYPRE_ParCSRPilutSolve,
                                           (HYPRE_PtrToParSolverFcn)HYPRE_ParCSRPilutSetup,
                                           precond);
            break;
          case 2:
            if (rk == 0) printf("SOLVER: ParaSails-BicGSTAB\n");
            ierr = HYPRE_ParaSailsCreate(comm, &precond);
            if (ierr) printf("ERROR: ParaSails-BicGSTAB\n");
            HYPRE_ParaSailsSetParams(precond, param->sai_threshold, param->max_levels);
            HYPRE_ParaSailsSetFilter(precond, param->sai_filter);
            HYPRE_ParaSailsSetSym(precond, param->sai_sym);

            HYPRE_ParCSRBiCGSTABSetPrecond(solver, (HYPRE_PtrToParSolverFcn)HYPRE_ParaSailsSolve,
                                           (HYPRE_PtrToParSolverFcn)HYPRE_ParaSailsSetup, precond);
            break;
          case 3:
            if (rk == 0) printf("SOLVER: Schwarz-PCG \n");
            HYPRE_SchwarzCreate(&precond);
            HYPRE_SchwarzSetVariant(precond, param->variant);
            HYPRE_SchwarzSetOverlap(precond, param->overlap);
            HYPRE_SchwarzSetDomainType(precond, param->domain_type);
            // HYPRE_SchwarzSetRelaxWeight(precond,param->schwarz_rlx_weight);
            /*HYPRE_BoomerAMGSetOverlap(precond, param->overlap);
            HYPRE_BoomerAMGSetVariant(precond, param->variant);
            HYPRE_BoomerAMGSetDomainType(precond, param->domain_type);*/

            HYPRE_BiCGSTABSetPrecond(solver, (HYPRE_PtrToSolverFcn)HYPRE_SchwarzSolve,
                                     (HYPRE_PtrToSolverFcn)HYPRE_SchwarzSetup, precond);

            break;

          default:
            break;
        }

        break;

      case 2: /*USE PCG AS SOLVER*/
        HYPRE_ParCSRPCGCreate(comm, &solver);
        HYPRE_ParCSRPCGSetTol(solver, param->solv_tol);
        HYPRE_PCGSetMaxIter(solver, param->solv_maxiter);
        HYPRE_PCGSetPrintLevel(solver, param->solv_prntlevel);

        HYPRE_PCGSetTwoNorm(solver, param->Two_Norm);
        HYPRE_PCGSetRelChange(solver, param->Rel_change);
        HYPRE_PCGGetPrecond(solver, &precond);
        switch (param->precond_id) {
          case 0:    // set BoomerAMG as preconditioner in PCG
            if (rk == 0) printf("SOLVER: AMG-PCG\n");
            HYPRE_BoomerAMGCreate(&precond);

            HYPRE_BoomerAMGSetInterpType(precond, param->amg_interptype);

            HYPRE_BoomerAMGSetTol(precond, param->amg_tol);
            HYPRE_BoomerAMGSetCoarsenType(precond, param->amg_coarsentype);
            HYPRE_BoomerAMGSetMeasureType(precond, param->measure_type);
            HYPRE_BoomerAMGSetStrongThreshold(precond, param->amg_strongthreshold);

            HYPRE_BoomerAMGSetTruncFactor(precond, param->trunc_factor);

            HYPRE_BoomerAMGSetMaxIter(precond, param->amg_maxiter);
            HYPRE_BoomerAMGSetCycleType(precond, param->cycle_type);

            HYPRE_BoomerAMGSetSmoothType(precond, param->smooth_type);
            HYPRE_BoomerAMGSetSmoothNumLevels(precond, param->smooth_num_levels);
            HYPRE_BoomerAMGSetSmoothNumSweeps(precond, param->smooth_num_sweeps);
            HYPRE_BoomerAMGSetMaxLevels(precond, param->max_levels);
            HYPRE_BoomerAMGSetMaxRowSum(precond, param->max_row_sum);

            HYPRE_BoomerAMGSetOverlap(precond, param->overlap);
            HYPRE_BoomerAMGSetVariant(precond, param->variant);
            HYPRE_BoomerAMGSetDomainType(precond, param->domain_type);

            HYPRE_PCGSetPrecond(solver, (HYPRE_PtrToSolverFcn)HYPRE_BoomerAMGSolve,
                                (HYPRE_PtrToSolverFcn)HYPRE_BoomerAMGSetup, precond);

            break;

          case 2:
            if (rk == 0) printf("SOLVER: ParaSails-PCG\n");
            ierr = HYPRE_ParaSailsCreate(comm, &precond);
            if (ierr) printf("ERROR: ParaSails-PCG\n");
            HYPRE_ParaSailsSetParams(precond, param->sai_threshold, param->max_levels);
            HYPRE_ParaSailsSetFilter(precond, param->sai_filter);
            HYPRE_ParaSailsSetSym(precond, param->sai_sym);
            HYPRE_PCGSetPrecond(solver, (HYPRE_PtrToSolverFcn)HYPRE_ParaSailsSolve,
                                (HYPRE_PtrToSolverFcn)HYPRE_ParaSailsSetup, precond);
            break;
          case 3:
            if (rk == 0) printf("SOLVER: Schwarz-PCG \n");
            HYPRE_SchwarzCreate(&precond);
            HYPRE_SchwarzSetVariant(precond, param->variant);
            HYPRE_SchwarzSetOverlap(precond, param->overlap);
            HYPRE_SchwarzSetDomainType(precond, param->domain_type);
            // HYPRE_SchwarzSetRelaxWeight(precond,param->schwarz_rlx_weight);
            /*HYPRE_BoomerAMGSetOverlap(precond, param->overlap);
            HYPRE_BoomerAMGSetVariant(precond, param->variant);
            HYPRE_BoomerAMGSetDomainType(precond, param->domain_type);*/

            HYPRE_PCGSetPrecond(solver, (HYPRE_PtrToSolverFcn)HYPRE_SchwarzSolve,
                                (HYPRE_PtrToSolverFcn)HYPRE_SchwarzSetup, precond);

            break;

          default:
            break;
        }

        break;
      case 5:
        if (rk == 0) printf("SOLVER: AMG \n");

        HYPRE_BoomerAMGCreate(&solver);
        HYPRE_BoomerAMGSetInterpType(solver, param->amg_interptype);
        HYPRE_BoomerAMGSetNumSamples(solver, param->gsmg_sample);
        HYPRE_BoomerAMGSetTol(solver, param->amg_tol);
        HYPRE_BoomerAMGSetCoarsenType(solver, param->amg_coarsentype);
        HYPRE_BoomerAMGSetMeasureType(solver, param->measure_type);
        HYPRE_BoomerAMGSetStrongThreshold(solver, param->amg_strongthreshold);
        HYPRE_BoomerAMGSetTruncFactor(solver, param->trunc_factor);
        HYPRE_BoomerAMGSetMaxIter(solver, param->amg_maxiter);
        HYPRE_BoomerAMGSetCycleType(solver, param->cycle_type);
        HYPRE_BoomerAMGSetRelaxWeight(solver, param->relax_weight);
        HYPRE_BoomerAMGSetOmega(solver, param->omega);

        HYPRE_BoomerAMGSetSmoothType(solver, param->smooth_type);
        HYPRE_BoomerAMGSetSmoothNumLevels(solver, param->smooth_num_levels);
        HYPRE_BoomerAMGSetSmoothNumSweeps(solver, param->smooth_num_sweeps);
        HYPRE_BoomerAMGSetMaxLevels(solver, param->max_levels);
        HYPRE_BoomerAMGSetMaxRowSum(solver, param->max_row_sum);

        HYPRE_BoomerAMGSetVariant(solver, param->variant);
        // HYPRE_BoomerAMGSetNumFunctions(solver, param->num_functions);
        // HYPRE_BoomerAMGSetSchwarzRlxWeight(solver,param->schwarz_rlx_weight);
        HYPRE_BoomerAMGSetOverlap(solver, param->overlap);
        HYPRE_BoomerAMGSetVariant(solver, param->domain_type);
        /* if(param->num_functions>1)
          HYPRE_BoomerAMGSetDofFunc(solver,param->dof_func);*/

        break;
      case 4:
        if (rk == 0) printf("SOLVER: AMG_HYBRID \n");

        HYPRE_ParCSRHybridCreate(&solver);
        HYPRE_ParCSRHybridSetTol(solver, param->amg_tol);

        HYPRE_ParCSRHybridSetCoarsenType(solver, param->coarsen_type);
        HYPRE_ParCSRHybridSetStrongThreshold(solver, param->strong_threshold);
        HYPRE_ParCSRHybridSetTruncFactor(solver, param->trunc_factor);
        HYPRE_ParCSRHybridSetDSCGMaxIter(solver, param->dscg_max_its);
        HYPRE_ParCSRHybridSetPCGMaxIter(solver, param->pcg_max_its);
        HYPRE_ParCSRHybridSetConvergenceTol(solver, param->cf_tol);
        HYPRE_ParCSRHybridSetSolverType(solver, param->solver_type);

        HYPRE_ParCSRHybridSetRelaxWeight(solver, param->relax_weight);
        HYPRE_ParCSRHybridSetOmega(solver, param->omega);
        HYPRE_ParCSRHybridSetMaxLevels(solver, param->max_levels);
        HYPRE_ParCSRHybridSetMaxRowSum(solver, param->max_row_sum);

        break;

      case 1: /*GMRES AS SOLVER*/
        HYPRE_ParCSRFlexGMRESCreate(comm, &solver);

        HYPRE_ParCSRFlexGMRESSetPrintLevel(solver, VERBOSE);

        HYPRE_ParCSRFlexGMRESSetKDim(solver, param->solv_kdim);
        HYPRE_ParCSRFlexGMRESSetMaxIter(solver, param->solv_maxiter);
        HYPRE_ParCSRFlexGMRESSetTol(solver, param->solv_tol);
        HYPRE_FlexGMRESSetLogging(solver, param->solv_log);
        HYPRE_ParCSRGMRESSetPrintLevel(solver, param->solv_prntlevel);

        // Set Preconditioner
        switch (param->precond_id) {
          case 0:    // set BoomerAMG as preconditioner
            if (rk == 0) printf("SOLVER: AMG-GMRES\n");
            HYPRE_BoomerAMGCreate(&precond);

            HYPRE_BoomerAMGSetInterpType(precond, param->amg_interptype);

            //  HYPRE_BoomerAMGSetTol(precond, param->amg_tol);
            HYPRE_BoomerAMGSetCoarsenType(precond, param->amg_coarsentype);
            HYPRE_BoomerAMGSetMeasureType(precond, param->measure_type);
            HYPRE_BoomerAMGSetStrongThreshold(precond, param->amg_strongthreshold);

            HYPRE_BoomerAMGSetTruncFactor(precond, param->trunc_factor);

            HYPRE_BoomerAMGSetMaxIter(precond, param->amg_maxiter);
            HYPRE_BoomerAMGSetCycleType(precond, param->cycle_type);

            HYPRE_BoomerAMGSetSmoothType(precond, param->smooth_type);
            HYPRE_BoomerAMGSetSmoothNumLevels(precond, param->smooth_num_levels);
            HYPRE_BoomerAMGSetSmoothNumSweeps(precond, param->smooth_num_sweeps);
            HYPRE_BoomerAMGSetMaxLevels(precond, param->max_levels);
            HYPRE_BoomerAMGSetMaxRowSum(precond, param->max_row_sum);

            HYPRE_BoomerAMGSetOverlap(precond, param->overlap);
            HYPRE_BoomerAMGSetVariant(precond, param->variant);
            HYPRE_BoomerAMGSetDomainType(precond, param->domain_type);

            HYPRE_ParCSRFlexGMRESSetPrecond(solver, (HYPRE_PtrToParSolverFcn)HYPRE_BoomerAMGSolve,
                                            (HYPRE_PtrToParSolverFcn)HYPRE_BoomerAMGSetup, precond);

            break;
          case 4:
            /*use diagonal scaling as preconditioner*/
            if (rk == 0) printf("SOLVER: DS-GMRES\n");
            precond = NULL;
            HYPRE_ParCSRFlexGMRESSetPrecond(solver, (HYPRE_PtrToParSolverFcn)HYPRE_ParCSRDiagScale,
                                            (HYPRE_PtrToParSolverFcn)HYPRE_ParCSRDiagScaleSetup,
                                            precond);

            break;
          case 1:
            /*Use PILUT as preconditioner*/
            if (rk == 0) printf("SOLVER: PILUT-GMRES\n");
            ierr = HYPRE_ParCSRPilutCreate(comm, &precond);

            HYPRE_ParCSRPilutSetMaxIter(precond, param->pmax_iter);
            if (param->drop_tol >= 0) HYPRE_ParCSRPilutSetDropTolerance(precond, param->drop_tol);
            if (param->nonzeros_to_keep >= 0)
              HYPRE_ParCSRPilutSetFactorRowSize(precond, param->nonzeros_to_keep);
            if (ierr) printf("ERROR: PILUT-GMRES \n");
            HYPRE_ParCSRFlexGMRESSetPrecond(solver, (HYPRE_PtrToParSolverFcn)HYPRE_ParCSRPilutSolve,
                                            (HYPRE_PtrToParSolverFcn)HYPRE_ParCSRPilutSetup,
                                            precond);
            break;

          case 3:
            if (rk == 0) printf("SOLVER: Schwarz-GMRES \n");
            HYPRE_SchwarzCreate(&precond);
            HYPRE_SchwarzSetVariant(precond, param->variant);
            HYPRE_SchwarzSetOverlap(precond, param->overlap);
            HYPRE_SchwarzSetDomainType(precond, param->domain_type);
            // HYPRE_SchwarzSetRelaxWeight(precond,param->schwarz_rlx_weight);
            /*HYPRE_BoomerAMGSetOverlap(precond, param->overlap);
            HYPRE_BoomerAMGSetVariant(precond, param->variant);
            HYPRE_BoomerAMGSetDomainType(precond, param->domain_type);*/

            HYPRE_ParCSRFlexGMRESSetPrecond(solver, (HYPRE_PtrToParSolverFcn)HYPRE_SchwarzSolve,
                                            (HYPRE_PtrToParSolverFcn)HYPRE_SchwarzSetup, precond);

            break;
          case 2:
            if (rk == 0) printf("SOLVER: ParaSails-GMRES\n");
            ierr = HYPRE_ParaSailsCreate(comm, &precond);

            if (ierr) printf("ERROR: ParaSails-GMRES\n");
            HYPRE_ParaSailsSetParams(precond, param->sai_threshold, param->max_levels);
            HYPRE_ParaSailsSetFilter(precond, param->sai_filter);
            HYPRE_ParaSailsSetSym(precond, param->sai_sym);

            HYPRE_ParCSRFlexGMRESSetPrecond(solver, (HYPRE_PtrToParSolverFcn)HYPRE_ParaSailsSolve,
                                            (HYPRE_PtrToParSolverFcn)HYPRE_ParaSailsSetup, precond);

            break;
        }

      default:

        break;
    }
  }

  void Solver(const MatriceMorse< double > &AA, KN_< double > &x, const KN_< double > &b) const {
    ffassert(&x[0] != &b[0]);
    epsr = (eps < 0) ? (epsr > 0 ? -epsr : -eps) : eps;
    int i, i1, i2;
    int *row, *row1;
    double *b_loc, *X_loc, *rhs, *xx;
    n = AA.n;

    rhs = (double *)malloc(sizeof(double) * n);
    xx = (double *)malloc(sizeof(double) * n);
    x = (double *)malloc(sizeof(double) * n);
    i1 = ilower;
    i2 = iupper;
    double t2, t1;

    for (i = 0; i < n; i++) rhs[i] = b[i];

    if (param->scale)
      for (i = 0; i < n; i++) rhs[i] *= scaletmpr[i];

    b_loc = (double *)malloc(n_loc * sizeof(double));
    X_loc = (double *)malloc(n_loc * sizeof(double));
    row = (int *)malloc(n_loc * sizeof(int));
    row1 = (int *)malloc(n_loc * sizeof(int));

    for (i = i1; i <= i2; i++) {
      // node = maptmp[i];
      b_loc[i - i1] = rhs[i];
      X_loc[i - i1] = 0.0;    // Initial Guest
      row[i - i1] = i;        // used to get results later
      row1[i - i1] = i;
    }
    int ierr = 0;
    HYPRE_IJVectorCreate(comm, ilower, iupper, &ij_B);
    HYPRE_IJVectorCreate(comm, ilower, iupper, &ij_X);

    HYPRE_IJVectorSetObjectType(ij_B, HYPRE_PARCSR);
    HYPRE_IJVectorSetObjectType(ij_X, HYPRE_PARCSR);

    HYPRE_IJVectorInitialize(ij_B);
    HYPRE_IJVectorInitialize(ij_X);

    HYPRE_IJVectorSetValues(ij_B, n_loc, row, b_loc);

    HYPRE_IJVectorSetValues(ij_X, n_loc, row, X_loc);

    HYPRE_IJVectorAssemble(ij_B);
    HYPRE_IJVectorAssemble(ij_X);

    ierr = HYPRE_IJVectorGetObject(ij_B, (void **)&par_B);

    ierr = HYPRE_IJVectorGetObject(ij_X, (void **)&par_X);

    switch (param->solver_id) {
      case 0:    // BICGSTAB solver
        if (param->timing) {
          time_index1 = hypre_InitializeTiming("BiCGSTAB SETUP");
          hypre_BeginTiming(time_index1);
        }
        HYPRE_ParCSRBiCGSTABSetup(solver, par_A, par_B, par_X);
        if (param->timing) {
          hypre_EndTiming(time_index1);
          hypre_PrintTiming("BiCGSTAB SETUP TIME", comm);
          hypre_FinalizeTiming(time_index1);
          hypre_ClearTiming( );
        }

        if (param->timing) {
          time_index1 = hypre_InitializeTiming("ParCSR BICGSTAB Solver");
          hypre_BeginTiming(time_index1);
        }

        HYPRE_ParCSRBiCGSTABSolve(solver, par_A, par_B, par_X);

        if (param->timing) {
          hypre_EndTiming(time_index1);
          hypre_PrintTiming("SOLVE PHASE TIMES", comm);
          hypre_FinalizeTiming(time_index1);
          hypre_ClearTiming( );
        }

        HYPRE_ParCSRBiCGSTABGetNumIterations(solver, &num_iter);
        HYPRE_ParCSRBiCGSTABGetFinalRelativeResidualNorm(solver, &final_res_norm);

        // HYPRE_ParCSRBiCGSTABDestroy(solver);
        break;
      case 1:    // GMRES Solver
        if (param->timing) {
          time_index1 = hypre_InitializeTiming("GMRES SETUP");
          hypre_BeginTiming(time_index1);
        }

        HYPRE_FlexGMRESSetup(solver, (HYPRE_Matrix)par_A, (HYPRE_Vector)par_B, (HYPRE_Vector)par_X);
        if (param->timing) {
          hypre_EndTiming(time_index);
          hypre_PrintTiming("SETUP PHASE TIME", comm);
          hypre_FinalizeTiming(time_index);
          hypre_ClearTiming( );
        }

        if (param->timing) {
          time_index1 = hypre_InitializeTiming("ParCSR GMRES Solver");
          hypre_BeginTiming(time_index1);
        }

        HYPRE_FlexGMRESSolve(solver, (HYPRE_Matrix)par_A, (HYPRE_Vector)par_B, (HYPRE_Vector)par_X);
        if (param->timing) {
          hypre_EndTiming(time_index1);
          hypre_PrintTiming("Solve phase times", comm);
          hypre_FinalizeTiming(time_index1);
          hypre_ClearTiming( );
        }

        HYPRE_GMRESGetNumIterations(solver, &num_iter);
        HYPRE_GMRESGetFinalRelativeResidualNorm(solver, &final_res_norm);

        // HYPRE_ParCSRGMRESDestroy(solver);
        break;
      case 2:    // PCG
        if (param->timing) {
          time_index1 = hypre_InitializeTiming("PCG SETUP");
          hypre_BeginTiming(time_index1);
        }

        HYPRE_PCGSetup(solver, (HYPRE_Matrix)par_A, (HYPRE_Vector)par_B, (HYPRE_Vector)par_X);
        if (param->timing) {
          hypre_EndTiming(time_index);
          hypre_PrintTiming("SETUP PHASE TIME", comm);
          hypre_FinalizeTiming(time_index);
          hypre_ClearTiming( );
        }

        if (param->timing) {
          time_index1 = hypre_InitializeTiming("ParCSR PCG Solve");
          hypre_BeginTiming(time_index1);
        }

        HYPRE_PCGSolve(solver, (HYPRE_Matrix)par_A, (HYPRE_Vector)par_B, (HYPRE_Vector)par_X);

        if (param->timing) {
          hypre_EndTiming(time_index1);
          hypre_PrintTiming("Solve phase times", comm);
          hypre_FinalizeTiming(time_index1);
          hypre_ClearTiming( );
        }

        HYPRE_ParCSRPCGGetNumIterations(solver, &dscg_num_its);
        HYPRE_ParCSRPCGGetFinalRelativeResidualNorm(solver, &final_res_norm);

        break;
      case 4:
        if (param->timing) {
          time_index1 = hypre_InitializeTiming("HYBRID SETUP");
          hypre_BeginTiming(time_index1);
        }
        HYPRE_ParCSRHybridSetup(solver, par_A, par_B, par_X);
        if (param->timing) {
          hypre_EndTiming(time_index);
          hypre_PrintTiming("SETUP PHASE TIME", comm);
          hypre_FinalizeTiming(time_index);
          hypre_ClearTiming( );
        }

        if (param->timing) {
          time_index1 = hypre_InitializeTiming("ParCSR Hybrid Solve");
          hypre_BeginTiming(time_index1);
        }

        HYPRE_ParCSRHybridSolve(solver, par_A, par_B, par_X);
        if (param->timing) {
          hypre_EndTiming(time_index1);
          hypre_PrintTiming("Solve phase times", comm);
          hypre_FinalizeTiming(time_index1);
          hypre_ClearTiming( );
        }

        HYPRE_ParCSRHybridGetNumIterations(solver, &num_iter);
        HYPRE_ParCSRHybridGetPCGNumIterations(solver, &pcg_num_its);
        HYPRE_ParCSRHybridGetDSCGNumIterations(solver, &dscg_num_its);
        HYPRE_ParCSRHybridGetFinalRelativeResidualNorm(solver, &final_res_norm);
        if (rk == 0) {
          printf("\n");
          printf("Iterations = %d\n", num_iter);
          printf("PCG_Iterations = %d\n", pcg_num_its);
          printf("DSCG_Iterations = %d\n", dscg_num_its);
          printf("Final Relative Residual Norm = %e\n", final_res_norm);
          printf("\n");
        }
        break;
      case 5:
        if (param->timing) {
          time_index1 = hypre_InitializeTiming("BoomerAMG  SETUP");
          hypre_BeginTiming(time_index1);
        }
        HYPRE_BoomerAMGSetup(solver, par_A, par_B, par_X);
        if (param->timing) {
          time_index1 = hypre_InitializeTiming("BoomerAMG SETUP Solve");
          hypre_BeginTiming(time_index1);
        }
        HYPRE_BoomerAMGSolve(solver, par_A, par_B, par_X);
        if (param->timing) {
          hypre_EndTiming(time_index1);
          hypre_PrintTiming("Solve phase times", comm);
          hypre_FinalizeTiming(time_index1);
          hypre_ClearTiming( );
        }
        break;

      default:
        break;
    }

    /*Reconstruction of vector*/
    delete[] iwork1;
    delete[] mapptr;
    iwork1 = new int[size + 1];
    mapptr = new int[size + 1];
    int disp = 0;
    for (i = 0; i < size; i++) {
      iwork1[i] = n / size;
      if (i == size - 1) iwork1[i] = (n / size) + n % size;
      mapptr[i] = disp;
      disp += iwork1[i];
    }
    mapptr[i] = disp;

    HYPRE_IJVectorGetValues(ij_X, n_loc, row1, X_loc);
    MPI_Gatherv(X_loc, n_loc, MPI_DOUBLE, xx, iwork1, mapptr, MPI_DOUBLE, 0, comm);
    MPI_Bcast(xx, n, MPI_DOUBLE, 0, comm);

    for (i = 0; i < n; i++) {
      x[i] = xx[i];
      if (param->scale) x[i] = x[i] * scaletmpc[i];
    }

    if (verbosity == 0) {
      MPI_Barrier(comm);
      t2 = dwalltime( );
      t2 = t2 - t1;
      MPI_Reduce(&t2, &t1, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
      MPI_Bcast(&t1, 1, MPI_DOUBLE, 0, comm);

      if (rk == 0) {
        printf("%s%18.6g%s", "TIME FOR SOLVING ", fabs(t1), "\n \n");
        printf("%s%18.6g%s", "RELATIVE RESIDU  ", final_res_norm, "\n \n");
        printf("%s%d%s", "NUMBER OF ITERATION ", num_iter, " \n \n");
      }
    }

    if (X_loc != NULL) free(X_loc);
    if (rhs != NULL) free(rhs);
    if (xx != NULL) free(xx);
    if (b_loc != NULL) free(b_loc);
    if (row != NULL) free(row);
    if (iwork1 != NULL) delete[] iwork1;
    if (mapptr != NULL) delete[] mapptr;
  }
  ~HypreSolver( ) {
    if (verbosity) {
      cout << "~HypreSolver():" << endl;
      HYPRE_IJMatrixDestroy(ij_A);
      HYPRE_IJVectorDestroy(ij_B);
      HYPRE_IJVectorDestroy(ij_X);
      switch (param->solver_id) {
        case 0:    // BICGSTAB solver
          HYPRE_ParCSRBiCGSTABDestroy(solver);
          switch (param->precond_id) {
            case 0:
              HYPRE_BoomerAMGDestroy(precond);
              break;
            case 5:
              HYPRE_ParaSailsDestroy(precond);
              break;
            case 6:
              HYPRE_SchwarzDestroy(precond);
              break;
            case 4:
              HYPRE_ParCSRPilutDestroy(precond);
              break;
            default:
              break;
          }
          break;

        case 1:    // GMRES Solver
          HYPRE_ParCSRFlexGMRESDestroy(solver);
          switch (param->precond_id) {
            case 0:
              HYPRE_BoomerAMGDestroy(precond);
              break;
            case 5:
              HYPRE_ParaSailsDestroy(precond);
              break;
            case 6:
              HYPRE_SchwarzDestroy(precond);
              break;
            case 4:
              HYPRE_ParCSRPilutDestroy(precond);
              break;
            default:
              break;
          }
          break;

        case 3:    // PCG

          HYPRE_ParCSRPCGDestroy(solver);
          switch (param->precond_id) {
            case 0:
              HYPRE_BoomerAMGDestroy(precond);
              break;
            case 5:
              HYPRE_ParaSailsDestroy(precond);
              break;
            case 6:
              HYPRE_SchwarzDestroy(precond);
              break;
            default:
              break;
          }

          break;
        case 4:    // AMG-Hybrid
          HYPRE_ParCSRHybridDestroy(solver);
          break;
        case 5:    // AMG
          HYPRE_BoomerAMGDestroy(solver);
          break;
        default:
          break;
      }
    }
  }
  void addMatMul(const KN_< Complex > &x, KN_< Complex > &Ax) const {
    ffassert(x.N( ) == Ax.N( ));
    Ax += (const MatriceMorse< Complex > &)(*this) * x;
  }
};

inline MatriceMorse< double >::VirtualSolver *BuildHypreSolver(DCL_ARG_SPARSE_SOLVER(double, A)) {
  if (verbosity > 9) cout << " BuildSolverHypre>" << endl;

  return new HypreSolver(*A, ds.data_filename, ds.lparams, ds.dparams, (MPI_Comm *)ds.commworld);
}

/* --FH:   class Init { public:
    Init();
    };*/

//  the 2 default sparse solver double and complex
DefSparseSolver< double >::SparseMatSolver SparseMatSolver_R;
// DefSparseSolver<Complex>::SparseMatSolver SparseMatSolver_C;
// the default probleme solver
TypeSolveMat::TSolveMat TypeSolveMatdefaultvalue = TypeSolveMat::defaultvalue;

bool SetDefault( ) {
  if (verbosity > 1) cout << " SetDefault sparse to default" << endl;
  DefSparseSolver< double >::solver = SparseMatSolver_R;

  TypeSolveMat::defaultvalue = TypeSolveMat::SparseSolver;
  return 1;
}

bool SetHypreSolver( ) {
  if (verbosity > 1) cout << " SetDefault sparse solver to Hyprempi" << endl;
  DefSparseSolver< double >::solver = BuildHypreSolver;
  TypeSolveMat::defaultvalue = TypeSolveMatdefaultvalue;
  return 1;
}
static void Load_Init( ) {

  SparseMatSolver_R = DefSparseSolver< double >::solver;
  if (verbosity > 1) cout << "\n Add: Hyprempi,  defaultsolver defaultsolverHyprempi" << endl;
  TypeSolveMat::defaultvalue = TypeSolveMat::SparseSolver;
  DefSparseSolver< double >::solver = BuildHypreSolver;
  if (!Global.Find("defaultsolver").NotNull( ))
    Global.Add("defaultsolver", "(", new OneOperator0< bool >(SetDefault));
  Global.Add("defaulttoHyprempi", "(", new OneOperator0< bool >(SetHypreSolver));
}

LOADFUNC(Load_Init)
