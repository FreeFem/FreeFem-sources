/*! \file   MM-DissectionSolver.cpp
    \brief  test rouinte of dissection solver reading Matrix Market format
    \author Atsushi Suzuki, Laboratoire Jacques-Louis Lions
    \date   Jun. 20th 2014
    \date   Jul. 12th 2015
    \date   Nov. 30th 2016
    \dahe   Apr. 24th 2018
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

#include "Compiler/OptionCompiler.hpp"
#include "Compiler/OptionLibrary.hpp"
#include "Driver/DissectionSolver.hpp"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <math.h>

#ifdef BLAS_MKL
#include <mkl_service.h>
#endif
#ifdef POSIX_THREADS
#include <pthread.h>
#endif
#include <iostream>
#include <sstream>
#include <fstream>
#include <complex>
using namespace std;

#ifdef POSIX_THREADS
static int _stat = (-1);

void *thread_child(void *arg)
{
  char buf[256];
  int *pid = (int *)arg;
  unsigned int mem_tmp, mem_min, mem_max;
  double avg_mem;
  avg_mem = 0.0;
  mem_min = (1U << 31) - 1;
  mem_max = 0U;
  int stat0, stat1;
  stat0 = _stat;
  unsigned int count = 0U;
//  fprintf(stderr, "thread_child forked\n");
  while(_stat != 0) {
    stat1 = _stat;
     if (stat1 == 1) {
       sprintf(buf, "/proc/%d/statm", *pid);
       ifstream fin(buf);
       fin >> mem_tmp;
       fin.close();
       if (mem_tmp > mem_max) {
	 mem_max = mem_tmp;
       }
       if (mem_tmp < mem_min) {
	 mem_min = mem_tmp;
       }
       avg_mem += (double)mem_tmp;
       count++;
     }
     if ((stat1 == (-1)) && (stat0 == 1)) {
       fprintf(stderr, 
	       "used memory :min: %14.8e  max: %14.8e avg: %14.8e count: %d\n", 
	       (double)mem_min * 4.0 / (1024.0 * 1024.0),
	       (double)mem_max * 4.0 / (1024.0 * 1024.0),
	       (avg_mem / (double)count) * 4.0 / (1024.0 * 1024.0),
	       count);
       count = 0U;
       avg_mem = 0.0;
       mem_min = (1U << 31) - 1;
       mem_max = 0U;
     }
     stat0 = stat1;
#ifdef POSIX_THREADS    
     usleep(1000);
#endif
  }
// fprintf(stderr, "thread_child join\n count = %ld\n", count);
  pthread_exit(arg);

  return (void *)NULL;
}
#endif
template<typename T>
bool generate_CSR(std::list<int>* ind_cols_tmp, std::list<T>* val_tmp, 
		  int nrow, int *nnz, int *mask, int *old2new,
		  int *irow, int *jcol, T* val, bool symmetrize)
{
  bool flag_modified = false;
  const T zero(0.0);
  //  ind_cols_tmp = new std::list<int>[nrow];
  //  val_tmp = new std::list<T>[nrow];
  int nnz1 = *nnz;
  for (int i = 0; i < *nnz; i++) {
    const int i0 = irow[i];
    const int j0 = jcol[i];
    const int ii = old2new[i0];
    const int jj = old2new[j0];
    if ((mask[i0] != 1) || (mask[j0] != 1)) {
      //      fprintf(stderr, "%d %d\n", i0, j0);
      nnz1--;
      continue;
    }
    //    fprintf(stderr, "%d %d -> %d %d \n", i0, j0, ii, jj);
    if (ind_cols_tmp[ii].empty()) {
      ind_cols_tmp[ii].push_back(jj);
      val_tmp[ii].push_back(val[i]);
    }
    else {
      if (ind_cols_tmp[ii].back() < jj) {
	ind_cols_tmp[ii].push_back(jj);
	val_tmp[ii].push_back(val[i]);
      }
      else {
	typename std::list<T>::iterator iv = val_tmp[ii].begin();
	std::list<int>::iterator it = ind_cols_tmp[ii].begin();
	for ( ; it != ind_cols_tmp[ii].end(); ++it, ++iv) {
	  if (*it == jj) {
	    fprintf(stderr, "already exits? (%d %d)\n", ii, jj);
	      break;
	  }
	  if (*it > jj) {
	    ind_cols_tmp[ii].insert(it, jj);
	    val_tmp[ii].insert(iv, val[i]);
	    break;
	  }
	}
      }
    }
  }
  // symmetrize
  if (symmetrize) {
    for (int i = 0; i < nrow; i++) {
      for (std::list<int>::iterator jt = ind_cols_tmp[i].begin();
	   jt != ind_cols_tmp[i].end(); ++jt) {
	const int jj = (*jt);
	bool flag = false;
	for (std::list<int>::iterator it = ind_cols_tmp[jj].begin();
	      it != ind_cols_tmp[jj].end(); ++it) {
	  if ((*it) == i) {
	    flag = true;
	    //	    fprintf(stderr, "%d %d symmetric position found\n", i, jj);
	    break;
	  }
	}
	if (!flag) {
	  flag_modified = true;
	  //	  fprintf(stderr, "%d %d need to be added\n", i, jj);
	  if (ind_cols_tmp[jj].back() < i) {
	    ind_cols_tmp[jj].push_back(i);
	    val_tmp[jj].push_back(zero);
	    //	    fprintf(stderr, "%d %d append\n", i, jj);
	    nnz1++;
	  }
	  else {
	    typename std::list<T>::iterator iv = val_tmp[jj].begin();
	    std::list<int>::iterator it = ind_cols_tmp[jj].begin();
	    for (; it != ind_cols_tmp[jj].end(); ++it, ++iv) {
	      if ((*it) > i) {
		ind_cols_tmp[jj].insert(it, i);
		val_tmp[jj].insert(iv, zero);
		nnz1++;
		//		fprintf(stderr, "%d %d inserted\n", i, jj);
		break;
	      }
	    }
	  }
	} // if (!flag);
      }
    }
  }
  if (symmetrize) {
    for (int i = 0; i < nrow; i++) {
      std::list<int>::iterator jt = ind_cols_tmp[i].begin();
      for ( ; jt != ind_cols_tmp[i].end(); ++jt) {
	const int jj = (*jt);
	bool flag = false;
	for (std::list<int>::iterator it = ind_cols_tmp[jj].begin();
	     it != ind_cols_tmp[jj].end(); ++it) {
	  if ((*it) == i) {
	    flag = true;
	    break;
	  }
	}
	if (!flag) {
	  fprintf(stderr, "%d %d position is not symmetric\n", i, jj);
	}
      }
    }
  }
  *nnz = nnz1;
  return flag_modified;
}

template
bool generate_CSR<double>(std::list<int>* ind_cols_tmp, std::list<double>* val_tmp, 
			  int nrow, int *nnz, int *mask, int *old2new,
			  int *irow, int *jcol, double* val, bool symmetrize);
template
bool generate_CSR<complex<double> >(std::list<int>* ind_cols_tmp,
				    std::list<complex<double> >* val_tmp, 
				    int nrow, int *nnz, int *mask, int *old2new,
				    int *irow, int *jcol, complex<double>* val,
				    bool symmetrize);

template<typename T>
void copy_CSR(int *indcols, int *ptrows, T* coefs, int nrow, 
	      bool upper_flag, bool isSym,
	      std::list<int>* ind_cols_tmp, std::list<T>* val_tmp)
{
  const T zero(0.0);
  ptrows[0] = 0;
  for (int i = 0; i < nrow; i++) {
    int k;
    int itmp = ind_cols_tmp[i].size();
    if (upper_flag) {
      if (ind_cols_tmp[i].front() == i) {
	ptrows[i + 1] = ptrows[i] + itmp;
	k = ptrows[i];
      }
      else {
	fprintf(stderr, "zero is added to diagonal : %d\n", i);
	ptrows[i + 1] = ptrows[i] + itmp + 1;
	indcols[ptrows[i]] = i;
	coefs[ptrows[i]] = zero;
	k = ptrows[i] + 1;
      }
    }
    else {
      k = ptrows[i];
      if (ind_cols_tmp[i].back() == i || (!isSym)) {
	ptrows[i + 1] = ptrows[i] + itmp;
      }
      else {
	fprintf(stderr, "zero is added to diagonal : %d\n", i);
	ptrows[i + 1] = ptrows[i] + itmp + 1;
	indcols[ptrows[i + 1] - 1] = i;
	coefs[ptrows[i + 1] - 1] = zero;
      }
    }
    std::list<int>::iterator it = ind_cols_tmp[i].begin();
    typename std::list<T>::iterator iv = val_tmp[i].begin();
    for ( ; it != ind_cols_tmp[i].end(); ++it, ++iv, k++) {
      indcols[k] = *it;
      coefs[k] = *iv;
    }
  } // loop : i
}

template
void copy_CSR<double>(int *indcols, int *ptrows, double* coefs, int nrow, 
	      bool upper_flag, bool isSym,
	      std::list<int>* ind_cols_tmp, std::list<double>* val_tmp);

template
void copy_CSR<complex<double> >(int *indcols, int *ptrows, complex<double>* coefs, int nrow, 
	      bool upper_flag, bool isSym,
	      std::list<int>* ind_cols_tmp, std::list<complex<double> >* val_tmp);

#if 0
inline
double SpMAX(double a, double b)
{
  return (a > b ? a : b);
}

template<typename T>
void normalize(const int type, const int n, 
	       const int *ptRows, const int *indCols,
	       const T* coefs0, T* coefs, T* u)
{
  const T zero(0.0);
  const T one(1.0);
  const T Wzero(0.0);
  T *v, *d;
  v = new T[n];
  d = new T[n];

  T alower, aupper, adiag;

  for (int i = 0; i < n; i++) {
    d[i] = zero;
    u[i] = zero;
    v[i] = zero;
  }

  for (int i = 0; i < n; i++) {
    for (int k = ptRows[i]; k < ptRows[i + 1]; k++) {
      int j = indCols[k];
      if (j == i) {
	u[i] = coefs0[k] < zero ? -coefs0[k] : coefs0[k];
	d[i] = u[i];
	v[i] = SpMAX(v[i],u[i]);
      }
      else {
	T acoef = coefs0[k] < zero ? -coefs0[k] : coefs0[k];
	v[i] = SpMAX(v[i], acoef);
	v[j] = SpMAX(v[j], acoef);
      }
    }
  }
  
  switch(type) {
  case DIAGONAL_SCALING:
    for (int i = 0; i < n; i++) {
      u[i] = ((u[i] != zero) ? sqrt(one / u[i]) : 
	      ((v[i] != zero) ? sqrt(one / v[i]) : one));
    }
    break;
  case KKT_SCALING:
    for (int i = 0 ; i < n; i++) {
      if (d[i] != zero) {
	u[i] = sqrt(one / u[i]);
      }
      else {
	T xtmp = Wzero;
	for (int m = ptRows[i]; m < ptRows[i + 1]; m++) {
	  const int j = indCols[m];
	  alower = coefs0[m];
	  int flag = 0;
	  for (int n = ptRows[j]; n < ptRows[j + 1]; n++) {
	    const int k = indCols[n];
	    if (k == j) {
	      adiag = coefs0[n];
	      if (adiag != Wzero) {
		flag++;
	      }
	      else {
		continue;
	      }
	    }
	    if (k == i) {
	      aupper = coefs0[n];
	      flag++;
	    }
	    if (flag == 2) {
	      break;
	    }
	  } // loop : n
	  if (flag == 2) {
	    xtmp += alower * aupper / adiag;
	  }
	} // loop : m
	if (xtmp == zero) {
	  u[i] = one;
	}
	else {
	  u[i] = one / sqrt(xtmp < zero ? -xtmp : xtmp);
	}
      } // if (d[i] != zero)
    }   // loop : i
    break;
  default: 
    for (int i = 0; i < n; i++) {
      u[i] = one;
    }
    break;
  } // witch (type)
  // Scaling
  if (type == NO_SCALING) {
    for (int i = 0; i < ptRows[n]; i++) {
      coefs[i] = coefs0[i];
    }
    // blas_copy<T>(_ptRows[n], coefs0, 1, &_coefs[0], 1);
    // const W* coefs0 to T* _coefs is not supported becuase of type conversion
  }
  else {
    for ( int i = 0; i < n; i++) {
      for (int k = ptRows[i]; k < ptRows[i+1]; k++) {
	int j = indCols[k];
	//      _coefs[k] = tohigher((1), coefs0[k] * W(u[i]) * W(u[j]));
	coefs[k] = coefs0[k] * u[i] * u[j];
      }
    }
  }
  delete [] d;
  delete [] v;
}

template void
normalize<double>(const int type, const int n, 
		  const int *ptRows, const int *indCols,
		  const double* coefs0, double* coefs, double* u);
#endif
int main(int argc, char **argv)
{
  int n, itmp, jtmp;
  char fname[256], fname1[256];
  char buf[1024];
  int nrow, nnz, flag, nnz1;
  int *ptrows, *indcols;
  int *irow, *jcol;
  double *val, *coefs;
  complex<double> *valc, *ccoefs;
  int decomposer;
  int num_threads;
  int scaling = 1;
  double eps_pivot;
  int numlevels = -1;
  int minNodes = 128;
  std::list<int>* ind_cols_tmp;
  std::list<double>* val_tmp;
  std::list<complex<double> >* val_tmpc;
  FILE *fp;
  bool isSym, isComplex;
  bool upper_flag = true;
  bool isWhole = false;
  bool kernel_detection_all = false;
  int *indx_excl;
  int nexcl = 0;
  bool excl_flag = false;
  bool flag_modified = false;
  bool assume_invertible = false; 
  if (argc < 6) {
    fprintf(stderr, "MM-dissection [data file] [decomposer] [num_threads] [eps_pivot] [num_levels] [scaling] [kerner_detection_all] [upper_flag] [minNodes]\n");
    exit(-1);
  }    
  strcpy(fname, argv[1]);
  decomposer = atoi(argv[2]);
  num_threads = atoi(argv[3]);
  eps_pivot = atof(argv[4]);
  numlevels = atof(argv[5]);
  if (argc >= 7) {
    scaling = atoi(argv[6]);
  }
  if (argc >= 8) {
    assume_invertible = atoi(argv[7]) > 0 ? true : false;
  }
  if (argc >= 9) {
    upper_flag = (atoi(argv[8]) == 1);
    isWhole = (atoi(argv[8]) == (-1));
  }
  if (argc >= 10) {
    strcpy(fname1, argv[9]);
    excl_flag = true;
  }
  if (argc >= 11) {
    minNodes = atoi(argv[10]);
  }

  // read from the file
  if ((fp = fopen(fname, "r")) == NULL) {
    fprintf(stderr, "fail to open %s\n", fname);
  }
  fgets(buf, 256, fp);
  //
  if (strstr(buf, "symmetric") != NULL) {
   isSym = true;
  }
  else {
    isSym = false;
    upper_flag = false;
  }

  if (strstr(buf, "complex") != NULL) {
   isComplex = true;
  }
  else {
    isComplex = false;
  }

  fprintf(stderr, "symmetric = %s\n", isSym ? "true " : "false");
  fprintf(stderr, "scaling = %d\n", scaling);
  fprintf(stderr, "upper = %s\n", upper_flag ? "true" : "false");
  if (kernel_detection_all) {
    fprintf(stderr, "kernel detection is activated for all submatrices\n");
  }
  if (excl_flag) {
    fprintf(stderr, "list of singular nodes %s\n", fname1);
  }
  while (1) {
    fgets(buf, 256, fp);
    if (buf[0] != '%') {
      sscanf(buf, "%d %d %d", &nrow, &itmp, &nnz);
      break;
    }
  }
  irow = new int[nnz];
  jcol = new int[nnz];
  nnz1 = 0;
  if (isComplex) {
    double xreal, ximag;
    valc = new complex<double>[nnz];
    if (upper_flag) {
      for (int i = 0; i < nnz; i++) {
	fscanf(fp, "%d %d %lf %lf", &jcol[i], &irow[i], &xreal, &ximag);
	valc[i] = complex<double>(xreal, ximag);
	irow[i]--;
	jcol[i]--;
	if (isSym && irow[i] > jcol[i]) {
	  fprintf(stderr, "exchanged : %d > %d\n", irow[i], jcol[i]);
	  itmp = irow[i];
	  irow[i] = jcol[i];
	  jcol[i] = itmp;
	}
      }
    }
    else {
      for (int i = 0; i < nnz; i++) {
	fscanf(fp, "%d %d %lf %lf", &irow[i], &jcol[i], &xreal, &ximag);
	valc[i] = complex<double>(xreal, ximag);
	irow[i]--;
	jcol[i]--;
	if (isSym && irow[i] < jcol[i]) {
	  fprintf(stderr, "exchanged : %d > %d\n", irow[i], jcol[i]);
	  itmp = irow[i];
	  irow[i] = jcol[i];
	  jcol[i] = itmp;
	}
      }
    }
  }
  else { // if (isComplex)
    val = new double[nnz];
    if (upper_flag) {
      for (int i = 0; i < nnz; i++) {
	fscanf(fp, "%d %d %lf", &jcol[i], &irow[i], &val[i]); // read lower
	irow[i]--;
	jcol[i]--;
	if (isSym && irow[i] > jcol[i]) {
	  fprintf(stderr, "exchanged : %d > %d\n", irow[i], jcol[i]);
	  itmp = irow[i];
	  irow[i] = jcol[i];
	  jcol[i] = itmp;
	}
      }
    }
    else {
      int ii = 0;
      int itmp, jtmp;
      double vtmp;
      for (int i = 0; i < nnz; i++) {
	fscanf(fp, "%d %d %lf", &itmp, &jtmp, &vtmp);
	if (vtmp != 0.0 || (itmp == jtmp)) {
	  // if (true) {   // 04 Apr.2018
	  irow[ii] = itmp - 1;
	  jcol[ii] = jtmp - 1;
	  val[ii] = vtmp;
	  ii++;
	}
	else {
	  nnz1++;
	}
      }
    } // else 
  } // if (isComplex)
  fclose (fp);
  if (nnz1 > 0) {
    fprintf(stderr, "%s %d : %d zero entries excluded %d -> %d\n",
	    __FILE__, __LINE__, nnz1, nnz, (nnz - nnz1));
    nnz = nnz - nnz1;
  }
  if (excl_flag) {
    if ((fp = fopen(fname1, "r")) == NULL) {
      fprintf(stderr, "fail to open %s\n", fname1);
    }
    fgets(buf, 256, fp);
    sscanf(buf, "# %d", &nexcl);
    indx_excl = new int[nexcl];
    for (int i = 0; i < nexcl; i++) {
      fgets(buf, 256, fp);
      sscanf(buf, "%d", &itmp);
      indx_excl[i] = itmp;
    }
    fclose(fp);
  }

  int *mask = new int[nrow];
  int *old2new = new int[nrow];
  for (int i = 0; i < nrow; i++) {
    mask[i] = 1;
  }
  for (int i = 0; i < nexcl; i++) {
    mask[indx_excl[i]] = 0;
  }
  itmp = 0;
  jtmp = nrow - nexcl;
  for (int i = 0; i < nrow; i++) {
    if (mask[i] == 1) {
      old2new[i] = itmp++;
    }
    else {
      old2new[i] = jtmp++;
    }
  }
  nrow -= nexcl;

  ind_cols_tmp = new std::list<int>[nrow];
  fprintf(stderr, "%s %d : getnerate_CSR\n", __FILE__, __LINE__);
  if (isComplex) {
    val_tmpc = new std::list<complex<double> >[nrow];
    nnz1 = nnz;
    flag_modified = generate_CSR<complex<double> >(ind_cols_tmp, val_tmpc, 
				   nrow, &nnz1, mask, old2new,
				   irow, jcol, valc, (!isSym));
  }
  else {
    val_tmp = new std::list<double>[nrow];
    nnz1 = nnz;
    flag_modified = generate_CSR<double>(ind_cols_tmp, val_tmp, 
			 nrow, &nnz1, mask, old2new,
			 irow, jcol, val, (!isSym));
  }
  if (flag_modified) {
    fprintf(stderr, "%s %d : matrix is not structual symmetric %d ->%d\n",
	    __FILE__, __LINE__, nnz, nnz1);
  }
  nnz = nnz1;
  delete [] irow;
  delete [] jcol;
  delete [] mask;
  delete [] old2new;
  if (isComplex) {
    delete [] valc;
  }
  else {
    delete [] val;
  }

  if (upper_flag) {
    for (int i = 0; i < nrow; i++) {
      if (ind_cols_tmp[i].front() != i) {
	nnz++;
      }
    }
  }
  else {
    for (int i = 0; i < nrow; i++) {
      if (ind_cols_tmp[i].back() != i) {
	nnz++;
      }
    }
  }
  fprintf(stderr, "%s %d : copy_CSR\n", __FILE__, __LINE__);
  ptrows = new int[nrow + 1];
  indcols = new int[nnz];
  if (isComplex) {
    ccoefs = new complex<double>[nnz];
    copy_CSR<complex<double> >(indcols, ptrows, ccoefs, 
			       nrow, upper_flag, isSym, 
			       ind_cols_tmp, val_tmpc);
  }
  else {
    coefs = new double[nnz];
    copy_CSR<double>(indcols, ptrows, coefs, 
		     nrow, upper_flag, isSym,
		     ind_cols_tmp, val_tmp);
  }
  delete [] ind_cols_tmp;

  if (isComplex) {
    delete [] val_tmpc;
  }
  else {
    delete [] val_tmp;
  }
  int pid = get_process_id();
#if 1 
  fprintf(stderr, "pid = %d\n", pid);
  sprintf(fname, "dissection.%04d.log", pid);
  fp = fopen(fname, "a");
#else
  fp = stderr;
#endif
  for (int i = 0; i < argc; i++) {
    fprintf(fp, "%s ", argv[i]);
  }
  fprintf(fp, "\n");
  fprintf(stderr, "%s %d : before pthread_create\n", __FILE__, __LINE__);
  void* results;
#ifdef POSIX_THREADS
  pthread_attr_t th_attr;
  pthread_t thread;
  pthread_attr_init(&th_attr);
  pthread_attr_setdetachstate(&th_attr, PTHREAD_CREATE_JOINABLE);       
  int pthid = pthread_create(&thread, &th_attr, 
			     &thread_child,
			     (void *)&pid);
  if (pthid != 0) {
    cout << "bad thread creation ? " << pid << endl;
    exit(0);
  }
  fprintf(stderr, "%s %d : after pthread_create\n", __FILE__, __LINE__);
#endif
  if (isWhole) {
    isSym = true;
    upper_flag = false;
  }
    
  if (isComplex) {
    DissectionSolver<complex<double>, double>*dslv = 
      new DissectionSolver<complex<double>, double>(num_threads, true, 0, fp);
    
    int called = 0;
    clock_t t0_cpu, t1_cpu, t2_cpu, t3_cpu, t4_cpu, t5_cpu;
    elapsed_t t0_elapsed, t1_elapsed, t2_elapsed, t3_elapsed, 
      t4_elapsed, t5_elapsed;
#ifdef BLAS_MKL
    mkl_set_num_threads(1);
#endif
#ifdef VECLIB
    setenv("VECLIB_MAXIMUM_THREADS", "1", true);
#endif
    t0_cpu = clock();
    get_realtime(&t0_elapsed);

    dslv->SymbolicFact(nrow, (int *)ptrows, (int *)indcols,
		       isSym,
		       upper_flag,
		       isWhole,
		       decomposer, numlevels, minNodes); 
    
    t1_cpu = clock();
    get_realtime(&t1_elapsed);
    
    t2_cpu = clock();
    get_realtime(&t2_elapsed);
#ifdef POSIX_THREADS
    _stat = 1;
    usleep(5000);
#endif
    dslv->NumericFact(0, (complex<double> *)ccoefs, scaling, 
		      eps_pivot, kernel_detection_all, 4, -1.0,
		      assume_invertible);
#ifdef POSIX_THREADS
    _stat = (-1);
#endif
    t3_cpu = clock();
    get_realtime(&t3_elapsed);
#ifdef POSIX_THREADS
    usleep(5000);
#endif
    fprintf(stderr, "%s %d : NumericFact() done\n", __FILE__, __LINE__);

    int n0;
    n0 = dslv->kern_dimension();
    fprintf(fp, "%s %d : ## kernel dimension = %d\n", __FILE__, __LINE__, n0);
    
    complex<double> *x = new complex<double>[nrow];
    complex<double> *y = new complex<double>[nrow];
    complex<double> *z = new complex<double>[nrow];
    for (int i = 0; i < nrow; i++) {
      y[i] = complex<double>((double)(i % 11), 0.0);
    }
    dslv->SpMV(y, x);
    dslv->SpMV(x, y);
    for (int i = 0; i < nrow; i++) {
      z[i] = y[i];
    }
    
    t4_cpu = clock();
    get_realtime(&t4_elapsed);  
#ifdef POSIX_THREADS
    _stat = 1;
    usleep(5000);
#endif
    dslv->SolveSingle(y, true, false, true); // with projection + scaling
#ifdef POSIX_THREADS
    _stat = (-1);
    usleep(5000);
#endif
    fprintf(stderr, "%s %d : SolveSingle() done\n", __FILE__, __LINE__);
    t5_cpu = clock();
    get_realtime(&t5_elapsed);
    double norm0, norm1;
    norm0 = 0.0;
    norm1 = 0.0;
    for (int i = 0; i < nrow; i++) {
      norm0 += x[i].real() * x[i].real() + x[i].imag() * x[i].imag();
      complex<double> ztmp = y[i] - x[i];
      norm1 += ztmp.real() * ztmp.real() + ztmp.imag() * ztmp.imag();
    }
    fprintf(fp, "%s %d : ## error    = %18.7e\n",
	    __FILE__, __LINE__,
	    sqrt(norm1 / norm0));
    
    dslv->SpMV(y, x);
  
    norm0 = 0.0;
    norm1 = 0.0;
    for (int i = 0; i < nrow; i++) {
      norm0 += z[i].real() * z[i].real() + z[i].imag() * z[i].imag();
      complex<double> ztmp = z[i] - x[i];
      norm1 += ztmp.real() * ztmp.real() + ztmp.imag() * ztmp.imag();
    }
    fprintf(fp, "%s %d : ## residual = %18.7e\n",
	    __FILE__, __LINE__,
	    sqrt(norm1 / norm0));
#ifdef POSIX_THREADS
    _stat = 0;
    pthread_attr_destroy(&th_attr);
    pthid = pthread_join(thread, &results);  
    if (pthid != 0) {
      cout << "bad thread join ? " << pthid << endl;
      exit(0);
    }
#endif
    fprintf(fp, "%s %d : ## symbolic fact    : cpu time = %.4e elapsed time = %.4e\n",
	    __FILE__, __LINE__,
	    (double)(t1_cpu - t0_cpu) / (double)CLOCKS_PER_SEC,
	    convert_time(t1_elapsed, t0_elapsed));
    
    fprintf(fp, "%s %d : ## numeric fact     : cpu time = %.4e elapsed time = %.4e\n",
	    __FILE__, __LINE__,
	    (double)(t3_cpu - t2_cpu) / (double)CLOCKS_PER_SEC,
	    convert_time(t3_elapsed, t2_elapsed));
    
    fprintf(fp, "%s %d : ## solve single RHS : cpu time = %.4e elapsed time = %.4e\n",
	    __FILE__, __LINE__,
	    (double)(t5_cpu - t4_cpu) / (double)CLOCKS_PER_SEC,
	    convert_time(t5_elapsed, t4_elapsed));

    delete dslv;
    delete [] ptrows;
    delete [] indcols;
    delete [] ccoefs;
    delete [] x;
    delete [] y;
    delete [] z;
  }  // if (isComplex)
  else { 
    DissectionSolver<double> *dslv = 
      new DissectionSolver<double>(num_threads, true, 0, fp);
    
    int called = 0;
    clock_t t0_cpu, t1_cpu, t2_cpu, t3_cpu, t4_cpu, t5_cpu;
    elapsed_t t0_elapsed, t1_elapsed, t2_elapsed, t3_elapsed, 
      t4_elapsed, t5_elapsed;
#ifdef BLAS_MKL
    mkl_set_num_threads(1);
#endif
#ifdef VECLIB
   setenv("VECLIB_MAXIMUM_THREADS", "1", true);
#endif
    t0_cpu = clock();
    get_realtime(&t0_elapsed);

    dslv->SymbolicFact(nrow, (int *)ptrows, (int *)indcols,
		       isSym,
		       upper_flag,
		       isWhole,
		       decomposer, numlevels, minNodes); 
    //                  sym, upper
    //    dslv->SaveMMMatrix(0, coefs);
    //    exit(-1);
    t1_cpu = clock();
    get_realtime(&t1_elapsed);
#ifdef POSIX_THREADS    
    _stat = 1;
    usleep(5000);
#endif
    t2_cpu = clock();
    get_realtime(&t2_elapsed);
    dslv->NumericFact(0, (double *)coefs, scaling, 
		      eps_pivot, kernel_detection_all,
		      4, -1.0,
		      assume_invertible);
    t3_cpu = clock();
    get_realtime(&t3_elapsed);
#ifdef POSIX_THREADS    
    _stat = (-1);
    usleep(5000);
#endif
    fprintf(stderr, "%s %d : NumericFact() done\n", __FILE__, __LINE__);
    int n0;
    n0 = dslv->kern_dimension();
    fprintf(fp, "%s %d : ## kernel dimension = %d\n", __FILE__, __LINE__, n0);
    
    double *x = new double[nrow];
    double *y = new double[nrow];
    double *z = new double[nrow];
    for (int i = 0; i < nrow; i++) {
      y[i] = (double)(i % 11);
    }
    dslv->SpMV(y, x);
    if (n0 > 0) {
      dslv->ProjectionKernelOrthSingle(y, "given", false);
    }
    dslv->SpMV(x, y);
    for (int i = 0; i < nrow; i++) {
      z[i] = y[i];
    }
#ifdef POSIX_THREADS    
    _stat = 1;
    usleep(5000);
#endif
    t4_cpu = clock();
    get_realtime(&t4_elapsed);  
    dslv->SolveSingle(y, true, false, true); // with projection + scaling
#ifdef POSIX_THREADS    
    _stat = (-1);
    usleep(5000);
#endif
    fprintf(stderr, "%s %d : SolveSingle() done\n", __FILE__, __LINE__);
    if (n0 > 0) {
      dslv->ProjectionImageSingle(y);
    }
    t5_cpu = clock();
    get_realtime(&t5_elapsed);

    double norm0, norm1;
    norm0 = 0.0;
    norm1 = 0.0;
    for (int i = 0; i < nrow; i++) {
      norm0 += x[i] * x[i];
      norm1 += (y[i] - x[i]) * (y[i] - x[i]);
    }
    fprintf(fp, "%s %d : ## error    = %18.7e\n",
	    __FILE__, __LINE__, sqrt(norm1 / norm0));

    dslv->SpMV(y, x);
    
    norm0 = 0.0;
    norm1 = 0.0;
    for (int i = 0; i < nrow; i++) {
      norm0 += z[i] * z[i];
      norm1 += (z[i] - x[i]) * (z[i] - x[i]);
    }
    fprintf(fp, "%s %d : ## residual = %18.7e\n",
	    __FILE__, __LINE__, sqrt(norm1 / norm0));
#ifdef POSIX_THREADS
    _stat = 0;
    pthread_attr_destroy(&th_attr);
    pthid = pthread_join(thread, &results);  
    if (pthid != 0) {
      cout << "bad thread join ? " << pthid << endl;
      exit(0);
    }
#endif
    fprintf(fp, "%s %d : ## symbolic fact    : cpu time = %.4e elapsed time = %.4e\n",
	    __FILE__, __LINE__, 
	    (double)(t1_cpu - t0_cpu) / (double)CLOCKS_PER_SEC,
	    convert_time(t1_elapsed, t0_elapsed));
    
    fprintf(fp, "%s %d : ## numeric fact     : cpu time = %.4e elapsed time = %.4e\n",
	    __FILE__, __LINE__, 
	    (double)(t3_cpu - t2_cpu) / (double)CLOCKS_PER_SEC,
	    convert_time(t3_elapsed, t2_elapsed));

    fprintf(fp, "%s %d : ## solve single RHS : cpu time = %.4e elapsed time = %.4e\n",
	    __FILE__, __LINE__, 
	    (double)(t5_cpu - t4_cpu) / (double)CLOCKS_PER_SEC,
	    convert_time(t5_elapsed, t4_elapsed));
    
    delete dslv;
    delete [] ptrows;
    delete [] indcols;
    delete [] coefs;
    delete [] x;
    delete [] y;
    delete [] z;
  }
  fclose(fp);

}
