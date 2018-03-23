/*! \file   C_KernDetect.cpp
    \brief  Kernel detection algorithm : symm <= DOI: 10.1002/nme.4729 / unsymm
    \author Atsushi. Suzuki, Laboratoire Jacques-Louis Lions
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

#include "Driver/C_KernDetect.hpp"
#include "Algebra/ColumnMatrix.hpp"
#include "Algebra/VectorArray.hpp"
#include "Compiler/OptionLibrary.h"
#include "Compiler/DissectionIO.hpp"

template<typename T, typename U>
void copy_matrix_permute_(const int lda, int n,
			  U *b, T *a, int *permute)
{
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      const int ij0 = permute[i] + permute[j] * lda;
      const int ij1 = i + j * n;
      b[ij1] = U(a[ij0]);
    }
  }
}

template
void copy_matrix_permute_<double, quadruple>(const int lda, int n,
					     quadruple *b, double *a,
					     int *permute);

template
void copy_matrix_permute_<float, double>(const int lda, int n,
					double *b, float *a,
					int *permute);

#ifndef NO_OCTRUPLE
template
void copy_matrix_permute_<quadruple, octruple>(const int lda, int n,
					       octruple *b, quadruple *a,
					       int *permute);
#endif
//

template<typename T, typename U>
void copy_matrix_permute_complex_(const int lda, int n,
				  complex<U> *b, complex<T> *a, int *permute)
{
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      const int ij0 = permute[i] + permute[j] * lda;
      const int ij1 = i + j * n;
      b[ij1] = complex<U>(U(a[ij0].real()),    // complex value is defined
			  U(a[ij0].imag()));   // by std::complex, then no
    }                                              // direct cast
  }
}

template
void copy_matrix_permute_complex_<double, quadruple>(const int lda, int n,
						     complex<quadruple> *b,
						     complex<double> *a,
						     int *permute);

template
void copy_matrix_permute_complex_<float, double>(const int lda, int n,
						 complex<double> *b,
						 complex<float> *a,
						 int *permute);

#ifndef NO_OCTRUPLE
template
void copy_matrix_permute_complex_<quadruple, octruple>(const int lda, int n,
						       complex<octruple> *b,
						       complex<quadruple> *a,
						       int *permute);
#endif
//

inline
void copy_matrix_permute(const int lda, int n,
			 quadruple *b, double *a, int *permute)
{
  copy_matrix_permute_<double, quadruple>(lda, n, b, a, permute);
}

inline
void copy_matrix_permute(const int lda, int n,
			 complex<quadruple> *b, complex<double> *a,
			 int *permute)
{
  copy_matrix_permute_complex_<double, quadruple>(lda, n, b, a, permute);
}

inline
void copy_matrix_permute(const int lda, int n,
			 double *b, float *a, int *permute)
{
  copy_matrix_permute_<float, double>(lda, n, b, a, permute);
}

inline
void copy_matrix_permute(const int lda, int n,
			 complex<double> *b, complex<float> *a,
			 int *permute)
{
  copy_matrix_permute_complex_<float, double>(lda, n, b, a, permute);
}

#ifndef NO_OCTRUPLE
inline
void copy_matrix_permute(const int lda, int n,
			 octruple *b, quadruple *a, int *permute)
{
  copy_matrix_permute_<quadruple, octruple>(lda, n, b, a, permute);
}

inline
void copy_matrix_permute(const int lda, int n,
			 complex<octruple> *b, complex<quadruple> *a,
			 int *permute)
{
  copy_matrix_permute_complex_<quadruple, octruple>(lda, n, b, a, permute);
}
#endif

template<typename T, typename U>
bool check_kern(const int n0, const int lda, const int n, T *a_ini,
		int *permute,
		const int dim_augkern, const U &eps,
		const U &eps_param, const bool flag_sym, U *errors)
{
  fprintf(stderr, "%s %d : specialized template not implemented\n",
	  __FILE__, __LINE__);
  return false;
}

template<>
bool check_kern<double, double>(const int n0, const int lda, const int n,
				double *a_ini, int *permute,
				const int dim_augkern, const double &eps,
				const double &eps_param, const bool flag_sym,
				double *errors)
{
  return check_kern_<double,
		     double,
		     quadruple,
		     quadruple>(n0, lda, n, a_ini, permute,
				dim_augkern, eps, eps_param,
				flag_sym, errors);
}

template<>
bool check_kern<complex<double>, double >(const int n0,
					  const int lda, const int n,
					  complex<double> *a_ini,
					  int *permute,
					  const int dim_augkern,
					  const double &eps,
					  const double &eps_param,
					  const bool flag_sym,
					  double *errors)
{
  return check_kern_<complex<double>,
		     double,
		     complex<quadruple>,
		     quadruple>(n0, lda, n, a_ini, permute,
				dim_augkern, eps, eps_param,
				flag_sym, errors);
}

template<>
bool check_kern<float, float>(const int n0, const int lda, const int n,
			      float *a_ini, int *permute,
			      const int dim_augkern, const float &eps,
			      const float &eps_param, const bool flag_sym,
			      float *errors)
{
  return check_kern_<float,
		     float,
		     double,
		     double>(n0, lda, n, a_ini, permute,
			     dim_augkern, eps, eps_param,
			     flag_sym, errors);
}

template<>
bool check_kern<complex<float>, float >(const int n0,
					  const int lda, const int n,
					  complex<float> *a_ini,
					  int *permute,
					  const int dim_augkern,
					  const float &eps,
					  const float &eps_param,
					  const bool flag_sym,
					  float *errors)
{
  return check_kern_<complex<float>,
		     float,
		     complex<double>,
		     double>(n0, lda, n, a_ini, permute,
				dim_augkern, eps, eps_param,
				flag_sym, errors);
}


#ifndef NO_OCTRUPLE
template<>
bool check_kern<quadruple, quadruple>(const int n0, const int lda, const int n,
				      quadruple *a_ini, int *permute,
				      const int dim_augkern,
				      const quadruple &eps,
				      const quadruple &eps_param,
				      const bool flag_sym,
				      quadruple *errors)
{
  return check_kern_<quadruple,
		     quadruple,
		     octruple,
		     octruple>(n0, lda, n, a_ini, permute,
			       dim_augkern, eps,
			       eps_param, flag_sym,
			       errors);
}

template<>
bool check_kern<complex<quadruple>, quadruple>(const int n0,
					       const int lda, const int n,
					       complex<quadruple> *a_ini,
					       int *permute,
					       const int dim_augkern,
					       const quadruple &eps,
					       const quadruple &eps_param,
					       const bool flag_sym,
					       quadruple *errors)
{
  return check_kern_<complex<quadruple>,
		     quadruple,
		     complex<octruple>,
		     octruple>(n0, lda, n, a_ini, permute,
			       dim_augkern, eps, eps_param,
			       flag_sym, errors);
}
#endif

// T may be std::complex of U and W is in higher precision than T
template<typename T, typename U, typename W, typename Y>
bool check_kern_(const int n0, const int lda, const int n, T *a_ini,
		 int *permute,
		 const int dim_augkern, const U &eps,
		 const U &eps_param, const bool flag_sym, U *errors)
{
  bool flag;
  //  W *a_q, *a_fq, *proj, *nsp, *nsp2;
  //  W *v, *alpha;
  const W zero(0.0);
  const W one(1.0);
  const W none(-1.0);
  
  ColumnMatrix<W> a_q(n, n); 
  ColumnMatrix<W> a_fq(n, n); 
  ColumnMatrix<W> proj(n, n); 
  ColumnMatrix<W> nsp(n, n); 
  ColumnMatrix<W> nsp2(n,n);
  VectorArray<W> v(n); 
  VectorArray<W> alpha(n);

  copy_matrix_permute(lda, n, a_fq.addrCoefs(), a_ini, permute);
  // duplicate of a_fq
  a_q.copy(a_fq);
  if (flag_sym) {
    full_ldlt<W>(n, a_fq.addrCoefs(), n);
  }
  else {
    full_ldu<W>(n, a_fq.addrCoefs(), n);
  }
  int n1 = n - n0;
  for (int j = 0; j < n0; j++) {
    for (int i = 0; i < n1; i++) {
      nsp(i, j) = a_q(i, (j + n1));  // nsp[i + j * n] = a_q[i + (j + n1) * n];

    }
    for (int i = n1; i < n; i++) {
      nsp(i, j) = zero;              // nsp[i + j * n] = zero;
    }
    nsp((n1 + j), j) = none;
  }
  full_fwbw_perturb_multi<W, U>(n1, n0, a_q.addrCoefs(), n,
				a_fq.addrCoefs(), nsp.addrCoefs(), dim_augkern,
				eps, flag_sym);
  // compute projection matrix
  for (int i = 0; i < n0; i++) {
    for (int j = 0; j <= i; j++) {
      proj(i, j) = blas_dot<W>(n,
			       nsp.addrCoefs() + (i * n), 1,
			       nsp.addrCoefs() + (j * n), 1); // lower
      proj(j, i) = blas_conj(proj(i, j));                  // upper
    }
  }
                         // Hermite symmteric with complex inner product
  full_ldlh<W>(n0, proj.addrCoefs(), n); 
  for (int m = (n0 - 1); m <= (n0 + 1); m++) {
    int k = n - m;
    U res_err = U(0.0);
    for (int j = 0; j < n; j++) {
      for (int i = 0; i < k; i++) {
	nsp2(i, j) = a_q(i, j);  // nsp2[i + j * n] = a_q[i + j * n];
      }
    }
    full_fwbw_perturb_multi<W, U>(k, n, a_q.addrCoefs(), n,
				a_fq.addrCoefs(), nsp2.addrCoefs(),
				dim_augkern, eps,
				  flag_sym);
    for (int j = 0; j < n; j++) {
      for (int i = k; i < n; i++) {
	nsp2(i, j) = zero;        //  nsp2[i + j * n] = zero;
      }
      nsp2(j, j) -= one;          //  nsp2[j + j * n] -= one;
    }
    for (int j = 0; j < k; j++) {
      for (int i = 0; i < n; i++) {
	v[i] = nsp2(i, j);      //  v[i] = nsp2[i + j * n];
      }
      U res = conv_prec<U, Y>(blas_l2norm<W, Y>(n, v.addrCoefs(), 1)); 
      res_err = res_err > res ? res_err : res;
    }
    for (int j = k; j < n; j++) {
      for (int i = 0; i < n; i++) {
	v[i] = nsp2(i, j);     	// v[i] = nsp2[i + j * n];
      }
      //      alpha = 1; beta = 0;
      blas_gemv<W>(CblasTrans, n, n0, one,
		   nsp.addrCoefs(), n, v.addrCoefs(), 1,
		   zero, alpha.addrCoefs(), 1);
      full_fwbw_part<W>(n0, proj.addrCoefs(), n, alpha.addrCoefs());
      //      alpha = -1; beta = 1;
      blas_gemv<W>(CblasNoTrans, n, n0, none,
		   nsp.addrCoefs(), n, alpha.addrCoefs(), 1,
		   one, v.addrCoefs(), 1);
      U res = conv_prec<U, Y>(blas_l2norm<W, Y>(n, v.addrCoefs(), 1)); 
      res_err = res_err > res ? res_err : res;
    }
    errors[m - n0 + 1] = res_err;
  }
  flag = false;
  if ((errors[0] > eps_param) && (errors[1] < eps_param) &&
      (errors[2] > eps_param)) {
    flag = true;
  }

  return flag;
}

template
bool check_kern_<double, double,
		 quadruple, quadruple>(const int n0, const int lda, const int n,
				       double *a_ini, int *permute,
				       const int dim_augkern, const double &eps,
				       const double &eps_param,
				       const bool flag_sym,
				       double *errors);
template
bool check_kern_<complex<double>,
		 double,
		 complex<quadruple>,
		 quadruple>(const int n0,
			    const int lda,
			    const int n,
			    complex<double> *a_ini,
			    int *permute,
			    const int dim_augkern,
			    const double &eps,
			    const double &eps_param,
			    const bool flag_sym,
			    double *errors);

template
bool check_kern_<float, float,
		 double, double>(const int n0, const int lda, const int n,
				 float *a_ini, int *permute,
				 const int dim_augkern, const float &eps,
				 const float &eps_param,
				 const bool flag_sym,
				 float *errors);
template
bool check_kern_<complex<float>,
		 float,
		 complex<double>,
		 double>(const int n0,
			 const int lda,
			 const int n,
			 complex<float> *a_ini,
			 int *permute,
			 const int dim_augkern,
			 const float &eps,
			 const float &eps_param,
			 const bool flag_sym,
			 float *errors);

#ifndef NO_OCTRUPLE
template
bool check_kern_<quadruple, quadruple,
		 octruple, octruple>(const int n0, const int lda,
				     const int n,
				     quadruple *a_ini, int *permute,
				     const int dim_augkern,
				     const quadruple &eps,
				     const quadruple &eps_param,
				     const bool flag_sym,
				     quadruple *errors);
template
bool
check_kern_<complex<quadruple>,
	    quadruple,
	    complex<octruple>,
	    octruple>(const int n0,
		      const int lda,
		      const int n,
		      complex<quadruple> *a_ini,
		      int *permute,
		      const int dim_augkern,
		      const quadruple &eps,
		      const quadruple &eps_param,
		      const bool flag_sym,
		      quadruple *errors);
#endif

template<typename T, typename U>
U check_matrixerr(const int lda, const int n,
		  T *a,
		  const int dim_augkern, const int k,
		  int *permute,
		  const U &eps,
		  const bool flag_sym)
{
  fprintf(stderr, "%s %d : specialized template not implemented\n",
	  __FILE__, __LINE__);
  return U(0.0);
};

template<>
double check_matrixerr<double, double>(const int lda, const int n,
				       double *a, const int dim_augkern,
				       const int k, int *permute,
				       const double &eps,
				       const bool flag_sym)
{
  return check_matrixerr_<double, double, quadruple>(lda, n, a, dim_augkern, k,
						     permute, eps, flag_sym);
}

template<>
double check_matrixerr<complex<double>, double>(const int lda,
						const int n,
						complex<double> *a,
						const int dim_augkern,
						const int k,
						int *permute,
						const double &eps,
						const bool flag_sym)
{
  return check_matrixerr_<complex<double>,
			  double,
			  complex<quadruple> >(lda, n, a, dim_augkern, k,
					       permute,
					       eps, flag_sym);
}

template<>
float check_matrixerr<float, float>(const int lda, const int n,
				       float *a, const int dim_augkern,
				       const int k, int *permute,
				       const float &eps,
				       const bool flag_sym)
{
  return check_matrixerr_<float, float, double>(lda, n, a, dim_augkern, k,
						     permute, eps, flag_sym);
}

template<>
float check_matrixerr<complex<float>, float>(const int lda,
						const int n,
						complex<float> *a,
						const int dim_augkern,
						const int k,
						int *permute,
						const float &eps,
						const bool flag_sym)
{
  return check_matrixerr_<complex<float>,
			  float,
			  complex<double> >(lda, n, a, dim_augkern, k,
					    permute,
					    eps, flag_sym);
}

#ifndef NO_OCTRUPLE
template<>
quadruple check_matrixerr<quadruple, quadruple>(const int lda, const int n,
						quadruple *a,
						const int dim_augkern,
						const int k, int *permute,
						const quadruple &eps,
						const bool flag_sym)
{
  return check_matrixerr_<quadruple,
			  quadruple,
			  octruple>(lda, n, a, dim_augkern,
				    k, permute, eps, flag_sym);
}

template<>
quadruple check_matrixerr<complex<quadruple>,
			  quadruple>(const int lda,
				     const int n,
				     complex<quadruple> *a,
				     const int dim_augkern,
				     const int k,
				     int *permute,
				     const quadruple &eps,
				     const bool flag_sym)
{
  return check_matrixerr_<complex<quadruple>,
			  quadruple, 
			  complex<octruple> >(lda, n, a, dim_augkern,
					      k, permute, eps, flag_sym);
}
#endif
//
// T may be std::complex of U and W is in higher precision than T
template<typename T, typename U, typename W>
U check_matrixerr_(const int lda, const int n,
		   T *a,
		   const int dim_augkern, const int k,
		   int *permute,
		   const U &eps,
		   const bool flag_sym)
  
{
  U error;
  const W one(1.0);
  ColumnMatrix<W> nsp2(n, n);
  ColumnMatrix<W> nsp3(n, n);
  ColumnMatrix<W> nsp4(n, n);
  VectorArray<W> v(n);       

  for (int i = 0; i < (n * n); i++) {
    nsp2.addrCoefs()[i] = one;
  }
  // permutation is given
  for (int j = 0; j < k; j++) {
    for (int i = 0; i < k; i++) {
      const int ij0 = permute[i] + permute[j] * lda;
      nsp2(i, j) = W(a[ij0]);        //      const int ij1 = i + j * n;
    }
  }
  nsp3.copy(nsp2);
  nsp4.copy(nsp2);
  if (flag_sym) {
    full_ldlt<W>(k, nsp3.addrCoefs(), n);      // factorization is done in
  }                                            // higher accurary
  else {
    full_ldu<W>(k, nsp3.addrCoefs(), n);
  }
  full_fwbw_perturb_multi<W, U>(k, k,
				nsp4.addrCoefs(), n,
				nsp3.addrCoefs(), nsp2.addrCoefs(),
				dim_augkern, eps,
				flag_sym);
  for (int i = 0; i < k; i++) {
    nsp2(i, i) -= one;        //    nsp2[i + i * n] -= one;
  }
  error = matrix_infty_norm<W, U>(k, nsp2.addrCoefs(), n); 

  return error;
}

template
double check_matrixerr_<double, double, quadruple>(const int lda, const int n,
						   double *a,
						   const int dim_augkern,
						   const int k, int *permute,
						   const double &eps,
						   const bool flag_sym);

template
double check_matrixerr_<complex<double>,
			double,
			complex<quadruple> >(const int lda,
					     const int n,
					     complex<double> *a,
					     const int dim_augkern,
					     const int k,
					     int *permute,
					     const double &eps,
					     const bool flag_sym);

template
float check_matrixerr_<float, float, double>(const int lda, const int n,
					      float *a,
					      const int dim_augkern,
					      const int k, int *permute,
					      const float &eps,
					      const bool flag_sym);

template
float check_matrixerr_<complex<float>,
			float,
			complex<double> >(const int lda,
					  const int n,
					  complex<float> *a,
					  const int dim_augkern,
					  const int k,
					  int *permute,
					  const float &eps,
					  const bool flag_sym);

#ifndef NO_OCTRUPLE
template
quadruple check_matrixerr_<quadruple,
			   quadruple,
			   octruple>(const int lda, const int n,
				     quadruple *a,
				     const int dim_augkern,
				     const int k, int *permute,
				     const quadruple &eps,
				     const bool flag_sym);

template
quadruple
check_matrixerr_<complex<quadruple>,
		 quadruple,
		 complex<octruple> >(const int lda,
				     const int n,
				     complex<quadruple> *a,
				     const int dim_augkern,
				     const int k,
				     int *permute,
				     const quadruple &eps,
				     const bool flag_sym);
#endif
//

template<typename T, typename U>
void verify_kernels(const int n0, const int n, const T *a_ini,
		    const bool isSym, const double eps, U *errors)
{
  const T one(1.0);
  const T zero(0.0);
  const T none(-1.0);

  ColumnMatrix<T> a_fact(n, n);
  ColumnMatrix<T> a_p(n, n);
  int *permute = new int[n];
  double fop;
  const int n1 = n - n0;
  blas_copy<T>((n * n), a_ini, 1, a_fact.addrCoefs(), 1);
  double pivot_ref = 0.0;
  for (int i = 0; i < n ; i++) {
    double tmp = blas_abs<T, double>(a_ini[i + i * n]); 
    pivot_ref = pivot_ref > tmp ? pivot_ref : tmp;
  }
  int nn0;
  if (isSym) {
    full_ldlt_permute<T, U>(&nn0, n0, n, a_fact.addrCoefs(), n,
			    &pivot_ref, permute, eps,
			    &fop);
  }
  else {
    full_ldu_permute<T, U>(&nn0, n0, n, a_fact.addrCoefs(), n,
			   &pivot_ref, permute, eps,
			   &fop);
  }
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n ; j++) {
      const int ij0 = permute[i] + permute[j] * n;
      a_p(i, j) = a_ini[ij0];         //      const int ij1 = i + j * n;
    }
  }
  //  compute [A_11^-1 A_12]
  //          [     -I     ]
  VectorArray<T> v(n);
  VectorArray<T> w(n);
  for (int j = n1; j < n; j++) {
    for (int i = 0; i < n1; i++) {
      const int ij0 = permute[i] + permute[j] * n;
      v[i] = a_ini[ij0];
    }
    full_fwbw_part<T>(n1, a_fact.addrCoefs(), n, v.addrCoefs());
    for (int i = n1; i < n; i++) {
      v[i] = zero;
    }
    v[j] = none;
    blas_gemv<T>(CblasNoTrans, n, n, one, a_p.addrCoefs(), n,
		 v.addrCoefs(), 1, zero,
		 w.addrCoefs(), 1);
    errors[j - n1] = blas_l2norm<T, U>(n, w.addrCoefs(), 1);
  }

  delete [] permute;
}

template
void verify_kernels<double, double>(const int n0, const int n,
				    const double *a_ini,
				    const bool isSym, const double eps,
				    double *errors);

template
void verify_kernels<complex<double>, double>(const int n0, const int n,
					     const complex<double> *a_ini,
					     const bool isSym, const double eps,
					     double *errors);

template
void verify_kernels<float, float>(const int n0, const int n,
				  const float *a_ini,
				  const bool isSym, const double eps,
				  float *errors);

template
void verify_kernels<complex<float>, float>(const int n0, const int n,
					   const complex<float> *a_ini,
					   const bool isSym, const double eps,
					   float *errors);

template
void verify_kernels<quadruple, quadruple>(const int n0, const int n,
					  const quadruple *a_ini,
					  const bool isSym, const double eps,
					  quadruple *errors);

template
void verify_kernels<complex<quadruple>,
		    quadruple>(const int n0, const int n,
			       const complex<quadruple> *a_ini,
			       const bool isSym,
			       const double eps,
			       quadruple *errors);
//


template<typename T, typename U>
void HouseholderVector_complex(int n, T *x, T *v, T *gamma)
{
  const U zero(0.0);
  const U one(1.0);
  const U two(2.0);
  const U onehalf(1.5);
  
  const T czero(zero, zero);
  const T cone(one, zero);
  const T ctwo(two, zero);
  const U pi(M_PI);
  U s = blas_l2norm2<T, U>((n - 1), &x[1], 1);
  
  v[0] = cone;
  for (int i = 1; i < n; i++) {
    v[i] = x[i];
  }
  if (s == zero) {
    *gamma = czero;
  }
  else {
//    U x0arg = std::arg(x[0]);
    U x0arg = atan2(x[0].imag(), x[0].real());
    const T alpha = complex<U>(cos(x0arg), sin(x0arg));
    U xabs = sqrt<U>(x[0].real() * x[0].real() +
		     x[0].imag() * x[0].imag());
    
    if ((x0arg >= pi / two) && (x0arg <  pi * onehalf)) {
      v[0] = alpha * (xabs + sqrt<U>(s));
    }
    else {
      v[0] = alpha * (-s) / (xabs + sqrt<U>(s));
    }
    const U v0r(v[0].real());
    const U v0i(v[0].imag());
    const U v0sq = v0r * v0r + v0i * v0i;
    *gamma = ctwo * v0sq / (s + v0sq); 
    T z = one / v[0];
    for (int i = 0; i < n; i++) {
      v[i] *= z;
    }
  }
}

template
void HouseholderVector_complex<complex<double>, double>(int n,
							complex<double> *x,
							complex<double> *v,
							complex<double> *gamma);

template
void HouseholderVector_complex<complex<quadruple>,
			       quadruple>(int n,
					  complex<quadruple> *x,
					  complex<quadruple> *v,
					  complex<quadruple> *gamma);

template
void HouseholderVector_complex<complex<float>, float>(int n,
							complex<float> *x,
							complex<float> *v,
							complex<float> *gamma);

//

template<typename T>
void HouseholderVector(int n, T *x, T *v, T *gamma)
{
  const T one(1.0);
  const T zero(0.0);
  const T two(2.0);
  const T s = blas_l2norm2<T, T>((n - 1), &x[1], 1);
  v[0] = one;
  for (int i = 1; i < n; i++) {
    v[i] = x[i];
  }
  if (s == zero) {
    *gamma = zero;
  }
  else {
    T z = sqrt<T>(x[0] * x[0] + s);
    if (x[0] <= zero) {
      v[0] = x[0] - z;
    }
    else {
      v[0] = (-s) / (x[0] + z);
    }
    *gamma = two * v[0] * v[0] / (s + v[0] * v[0]);
    z = one / v[0];
    for (int i = 0; i < n; i++) {
      v[i] *= z;
    }
  }
}

template<>
void HouseholderVector<complex<double> >(int n,
					 complex<double> *x,
					 complex<double> *v,
					 complex<double> *gamma)
{
  HouseholderVector_complex<complex<double>, double>(n, x, v, gamma);
}

template<>
void HouseholderVector<complex<quadruple> >(int n,
					    complex<quadruple> *x,
					    complex<quadruple> *v,
					    complex<quadruple> *gamma)
{
  HouseholderVector_complex<complex<quadruple>, quadruple>(n, x, v, gamma);
}

template<>
void HouseholderVector<complex<float> >(int n,
					 complex<float> *x,
					 complex<float> *v,
					 complex<float> *gamma)
{
  HouseholderVector_complex<complex<float>, float>(n, x, v, gamma);
}


template
void HouseholderVector<double>(int n, double *x, double *v, double *gamma);

template
void HouseholderVector<quadruple>(int n, quadruple *x, quadruple *v,
				  quadruple *gamma);

template
void HouseholderVector<float>(int n, float *x, float *v, float *gamma);

#if 0
template
void HouseholderVector<complex<double> >(int n,
					 complex<double> *x,
					 complex<double> *v,
					 complex<double> *gamma);

template
void HouseholderVector<complex<quadruple> >(int n,
					    complex<quadruple> *x,
					    complex<quadruple> *v,
					    complex<quadruple> *gamma);
#endif
//

template<typename T>
void HouseholderReflection(int n, T *a, int lda, T *v, T *w, const T &gamma)
{
  const T one(1.0);
  const T zero(0.0);
  blas_gemv<T>(CblasConjTrans, n, n, one, a, lda, v, 1, zero, w, 1);
  T ngamma = (-gamma);
  blas_gerc<T>(n, n, ngamma, v, 1, w, 1, a, lda);
}

template
void HouseholderReflection<double>(int n, double *a, int lda, double *v,
				   double *w, const double &gamma);

template
void HouseholderReflection<complex<double> >(int n, complex<double> *a, int lda,
					     complex<double> *v,
					     complex<double> *w,
					     const complex<double> &gamma);

template
void HouseholderReflection<quadruple>(int n, quadruple *a, int lda,
				      quadruple *v,
				      quadruple *w, const quadruple &gamma);

template
void
HouseholderReflection<complex<quadruple> >(int n,
					   complex<quadruple> *a, int lda,
					   complex<quadruple> *v,
					   complex<quadruple> *w,
					   const complex<quadruple> &gamma);
template
void HouseholderReflection<float>(int n, float *a, int lda, float *v,
				   float *w, const float &gamma);

template
void HouseholderReflection<complex<float> >(int n, complex<float> *a, int lda,
					     complex<float> *v,
					     complex<float> *w,
					     const complex<float> &gamma);

//

// T may be complex of U
template<typename T, typename U>
int hqr_pivot(const int n, T *a, int *permute)
{
  const T Tzero(0.0);
  const U Uzero(0.0);
  int n0, k;

  VectorArray<U> cc(n);
  VectorArray<T> col(n);
  VectorArray<T> v(n);
  VectorArray<T> w(n);

  for (int i = 0 ; i < n; i++) {
    cc[i] = blas_l2norm2<T, U>(n, &a[i * n], 1); // norm2 returns double value
    permute[i] = i;
  }
  k = 0;
  {
    U tmp(0.0);
    for (int i = 0; i < n; i++) {
      if (cc[i] > tmp) {
	tmp = cc[i];
	k = i;         // find the first entry that attains the maximum value
      }
    } // loop : i
  }
  n0 = 0;
  for (int m = 0; m < n; m++) {
    if (k > m) {   // swap k-th and m-th columns of A[]
      int kk = permute[m];
      permute[m] = permute[k];
      permute[k] = kk;
      for (int i = 0; i < n; i++) {
	col[i] = a[i + m * n];
      }
      for (int i = 0; i < n; i++) {
	a[i + m * n] = a[i + k * n];
      }
      for (int i = 0; i < n; i++) {
	a[i + k * n] = col[i];
      }
      U c = cc[m];
      cc[m] = cc[k];
      cc[k] = c;
    } // if (k > m) 
    int nm = n - m;
    T gamma;
    HouseholderVector<T>(nm, &a[m + m * n], v.addrCoefs(), &gamma);

    HouseholderReflection<T>(nm, &a[m + m * n], n,
			     v.addrCoefs(), w.addrCoefs(), gamma);
    for (int i = (m + 1); i < n; i++) {
      a[i + m * n] = Tzero; //  = v[i] : to keep Householder matrix
    }
    for (int i = (m + 1); i < n; i++) {
      cc[i] = blas_l2norm2<T, U>((n - m - 1), &a[m + 1 + i * n], 1);
    }
    U tt(0.0);
    for (int i = (m + 1); i < n; i++) {
      if (cc[i] > tt) {
	tt = cc[i];
	k = i;       // find the first entry that attains the maximum value
      }
    } // loop : i
    if (tt == Uzero) {
      n0 = n - (m + 1);
      break;
    }
  } // loop : m

  return (n - n0);
}

template
int hqr_pivot<double, double>(const int n, double *a, int *permute);

template
int hqr_pivot<complex<double>, double>(const int n, complex<double> *a,
				       int *permute);

template
int hqr_pivot<quadruple, quadruple>(const int n, quadruple *a, int *permute);

template
int hqr_pivot<complex<quadruple>, quadruple>(const int n,
					     complex<quadruple> *a,
					     int *permute);
template
int hqr_pivot<float, float>(const int n, float *a, int *permute);

template
int hqr_pivot<complex<float>, float>(const int n, complex<float> *a,
				       int *permute);

//


template<typename T, typename U>
bool ComputeDimKernel(int *n0, bool *flag_2x2, const T *a_, const int n, 
		      const bool sym_flag,
		      const int dim_augkern,
		      const U eps_machine, // for perturbation
		      const double eps_piv,
		      const bool verbose,
		      FILE *fp)
{
  const T zero(0.0);
  const U Uzero(0.0);
  // dimension of the image of the matrix a is at least one
  int nn0, n1, n2;
  bool flag;
  int n_dim = n + 1;
  int *permute = new int[n_dim];

  ColumnMatrix<T> a0(n_dim, n_dim);
                            // sizeof(long double) = sizeof(double) in Windows
  ColumnMatrix<T> a1(n_dim, n_dim);
  VectorArray<T> aa_diag(n_dim);
  VectorArray<U> rr(n_dim - 1);

  int *permute_d = new int[n_dim];
  int *permute_q = new int[n_dim];
  
  a1(n, n) = zero;
  //  a1.ZeroClear();
  // emulate numerical error from floating point operations
  if (sym_flag) {
    for (int i = 0; i < n; i++) {
      T tmp = zero;
      for (int j = 0; j < n; j++) {
	a1(i, j) = a_[i + j * n];
	tmp += a_[i + j * n] + T(random_bool() ? eps_machine : Uzero);
      }
      a1(i, n) = tmp;
      a1(n, i) = tmp; // keep symmetry of the matrix
      a1(n, n) += tmp + T(random_bool() ? eps_machine : Uzero);
    }
  }
  else {
    for (int i = 0; i < n; i++) {
      T tmp = zero;
      for (int j = 0; j < n; j++) {
	a1(i, j) = a_[i + j * n];
	tmp += a_[i + j * n] + T(random_bool() ? eps_machine : Uzero);
      }
      a1(i, n) = tmp;
      a1(n, n) += tmp + T(random_bool() ? eps_machine : Uzero);
    }
    for (int j = 0; j < n; j++) {
      T tmp = zero;
      for (int i = 0; i < n; i++) {
	tmp += a_[i + j * n] + T(random_bool() ? eps_machine : Uzero);
      }
      a1(n, j) = tmp;
    }
  }
  diss_printf(verbose, fp, "%s %d : machine eps = %s\n", __FILE__, __LINE__,
	      tostring<U>(eps_machine).c_str());
  a0.copy(a1);
  //
  n1 = hqr_pivot<T, U>(n_dim, a0.addrCoefs(), permute);
  diss_printf(verbose, fp,
	      "%s %d dimension of the image deteced by d_hqr_pivot() is %d\n",
	    __FILE__, __LINE__, n1);
  for (int i = 0; i < n1; i++) {
    aa_diag[i] = a0(i, i);
  }
#ifdef DEBUG_QR
  diss_printf(verbose, fp, "matrix\n");
  for (int i = 0; i < n_dim; i++) {
    diss_printf(verbose, fp, "%d : ", i); 
    for (int j = 0; j < n_dim; j++) {
      diss_printf(verbose, fp, "%s ", tostring<T>(aa_diag[i]).c_str());
    }
    diss_printf(verbose, fp, "\n");
  }
#else
  diss_printf(verbose, fp, "%s %d : diagonal entries of QR factorization\n",
	      __FILE__, __LINE__);
  for (int i = 0; i < n_dim; i++) {
    diss_printf(verbose, fp, "%d : %d %s\n", i, permute[i],
		tostring<T>(aa_diag[i]).c_str());
  }
#endif
  list<int> pos_gap;
  vector<int> kernel_dim;
  pos_gap.push_back(dim_augkern);

  flag = false;
  n2 = n1;
  for (int i = 0; i < (n1 - 1); i++) {
    if (blas_abs<T, U>(aa_diag[i]) > (eps_machine /
				      conv_prec<double, U>(sqrt<double>(eps_piv)))) {
      rr[i] = blas_abs<T, U>(aa_diag[i + 1] / aa_diag[i]);
    }
    else {
      n2 = i + 1;
      flag = true;
      break;
    }
  }

  // find largest gap inside of invertible part
  U aug_diff = U(1.0);
  for (int i = 0; i < (dim_augkern - 1); i++) {
    if (aug_diff > rr[i]) {
      aug_diff = rr[i];
    }
  }

  diss_printf(verbose, fp, "%s %d : aug_diff = %s esp = %.12e\n",
	      __FILE__, __LINE__, tostring<U>(aug_diff).c_str(), eps_piv);
  for (int i = (dim_augkern - 1); i < (n2 - 1); i++) {
    if (rr[i] < (aug_diff * U(eps_piv))) {
      diss_printf(verbose, fp, "rr[%d] : %s\n", i, tostring<U>(rr[i]).c_str());
    }
  }

  pos_gap.sort();
  pos_gap.unique();
  if (pos_gap.size() == 1 && (!flag)) {
    pos_gap.push_back(n1 - 1);
  }
  
  int n_dim1 = flag ? n2 : n_dim;

  for (list<int>::const_iterator it = pos_gap.begin(); it != pos_gap.end();
       ++it) {
    kernel_dim.push_back(n_dim1 - (*it));
  }

  diss_printf(verbose, fp, "%s %d : kernel_dim = %d : ",
	      __FILE__, __LINE__, (int)kernel_dim.size());
  for (vector<int>::const_iterator it = kernel_dim.begin();
       it != kernel_dim.end();
       ++it) {
    diss_printf(verbose, fp, "%d ", *it);
  }
  diss_printf(verbose, fp, "\n");
  ColumnMatrix<T> a_fact(n_dim1, n_dim1);
  for (int i = 0; i < n_dim1; i++) {
    for (int j = 0; j < n_dim1; j++) {
      //      a_fact[i + j * n_dim1] = a1[permute[i] + permute[j] * n_dim];
      a_fact(i, j) = a1(permute[i], permute[j]);
    }
  }
  
  //  if (kernel_dim.size() > 1) {
  flag = VerifyDimKernel(&nn0, n_dim1, a_fact.addrCoefs(), kernel_dim,
			 sym_flag, dim_augkern, eps_machine, verbose, fp);

  nn0 += (n_dim - n2) - 1;
  //  nn0--;   // detection of kernel is done for dim_augkern > 0 and at least 
  if (nn0 > 0) {
    U *errors_k = new U[nn0 + 1];
    vector<int> nnn;
    if (nn0 > 2) {
      nnn.push_back(nn0 - 1);
    }
    nnn.push_back(nn0);
    nnn.push_back(nn0 + 1);
    for (vector<int>::const_iterator it = nnn.begin(); it != nnn.end();
	 ++it) {
      const int nn = (*it);
      const double eps_machine_double = todouble<U>(eps_machine);
      verify_kernels(nn, n, a_, sym_flag, eps_machine_double, errors_k);
      diss_printf(verbose, fp, "errors of %d kernels\n", nn);
      for (int i = 0; i < nn; i++) {
	diss_printf(verbose, fp,
		    "%d : %s\n", i, tostring<U>(errors_k[i]).c_str());
      }
    }
    delete [] errors_k;
  } // if (nn0 > 0)
  delete [] permute;
  delete [] permute_d;
  delete [] permute_q;
  diss_printf(verbose, fp,
	      "%s %d : dimension of the kernel is %d\n",
	      __FILE__, __LINE__, nn0);
  *n0 = nn0;
  *flag_2x2 = false;

  return flag;
}

template
bool ComputeDimKernel<double, double>(int *n0, bool *flag_2x2,
				      const double *a_,
				      const int n, 
				      const bool sym_flag,
				      const int dim_augkern,
				      const double eps_machine,
				      const double eps_piv,
				      const bool verbose,
				      FILE *fp);

template
bool ComputeDimKernel<complex<double>, 
		      double>(int *n0, bool *flag_2x2,
			      const complex<double> *a_,
			      const int n, 
			      const bool sym_flag,
			      const int dim_augkern,
			      const double eps_machine,
			      const double eps_piv,
			      const bool verbose,
			      FILE *fp);

template
bool ComputeDimKernel<quadruple, 
		      quadruple>(int *n0, bool *flag_2x2,
				 const quadruple *a_,
				 const int n, 
				 const bool sym_flag,
				 const int dim_augkern,
				 const quadruple eps_machine,
				 const double eps_piv,
				 const bool verbose,
				 FILE *fp);

template
bool ComputeDimKernel<complex<quadruple>,
		      quadruple>(int *n0, bool *flag_2x2,
				 const complex<quadruple> *a_,
				 const int n, 
				 const bool sym_flag,
				 const int dim_augkern,
				 const quadruple eps_machine,
				 const double eps_piv,
				 const bool verbose,
				 FILE *fp);

template
bool ComputeDimKernel<float, float>(int *n0, bool *flag_2x2,
				    const float *a_,
				    const int n, 
				    const bool sym_flag,
				    const int dim_augkern,
				    const float eps_machine,
				    const double eps_piv,
				    const bool verbose,
				    FILE *fp);

template
bool ComputeDimKernel<complex<float>, 
		      float>(int *n0, bool *flag_2x2,
			      const complex<float> *a_,
			      const int n, 
			      const bool sym_flag,
			      const int dim_augkern,
			      const float eps_machine,
			      const double eps_piv,
			      const bool verbose,
			      FILE *fp);


//

template<typename T, typename U>
bool VerifyDimKernel(int *nn0_,
		     int n_dim, T* a_fact,
		     vector<int> &kernel_dim,
		     const bool sym_flag,
		     const int dim_augkern,
		     const U eps_machine,
		     const bool verbose,
		     FILE *fp)
{
  int nn, nn0;
  int k;
  double pivot_ref, pivot_ref_q;
  int *permute_q = new int[n_dim];
  ColumnMatrix<T> a2(n_dim, n_dim);
  U *errors = new U[6];

  blas_copy<T>((n_dim * n_dim), a_fact, 1, a2.addrCoefs(), 1);
  //  for (int i = 0; i < (n_dim * n_dim); i++) {
  //    a2[i] = a_fact[i];
  //  }
  pivot_ref = 0.0;
  k = 0;
  for (int i = 0; i < n_dim; i++) {
    double dtmp = blas_abs<T, double>(a_fact[k]); // accuracy? : 14 Jul.2015 Atsushi
    k += (n_dim + 1);
    pivot_ref = (pivot_ref < dtmp ? dtmp : pivot_ref);
  } 
  pivot_ref_q = pivot_ref;
  diss_printf(verbose, fp,
	      "pivot_ref_q = %.12e\n", pivot_ref_q);
  const int nn1 = 1; // matrix has at least one dimensional kernel
          // qfull_sym_gauss_part is only used to get permute_q[]

  if (sym_flag) {
    double fop;
    const double eps0 = todouble<U>(eps_machine);
    full_ldlt_permute<T, U>(&nn, nn1, n_dim, a_fact, n_dim, &pivot_ref_q,
			    permute_q, eps0, &fop);
  }
  else {
    double fop;
    const double eps0 = todouble<U>(eps_machine);
    full_ldu_permute<T, U>(&nn, nn1, n_dim, a_fact, n_dim, &pivot_ref_q,
			   permute_q, eps0, &fop);
  }
  // question: eps_machine is ok?  : 06 Jan.2013 
  // => jump between diagonal is within double precision : 27 Jan.2013
  diss_printf(verbose, fp,
	      "%s %d : permutation : ", __FILE__, __LINE__);
  for (int i = 0; i < n_dim; i++) {
    diss_printf(verbose, fp, "%d ", permute_q[i]);
  }
  diss_printf(verbose, fp, "\n");
  vector<int> dims(n_dim);
  for (int i = 0; i < n_dim; i++) {
    dims[i] = i + 1;
  }
  vector<U> errors_image(dims.size());
  for (int i = 0; i < dims.size(); i++) {
    errors_image[i] = check_matrixerr<T, U>(n_dim, n_dim,
					    a2.addrCoefs(), dim_augkern,
					    dims[i],
					    permute_q,
					    eps_machine,
					    sym_flag);
  }
  for (int i = 0; i < dims.size(); i++) {
    diss_printf(verbose, fp,
		"%d : %s\n", dims[i], tostring<U>(errors_image[i]).c_str());
  }
  U err_image = errors_image[0];
  for (int i = 1; i < dims.size(); i++) {
    if (dims[i] > dim_augkern) {
      break;
    }
    err_image = err_image < errors_image[i] ? errors_image[i] : err_image;
  }
  int count_err_image_updated = 0;
  int dim_image = dim_augkern;
  diss_printf(verbose, fp,
	      "err_image = %s ", tostring<U>(err_image).c_str());
  U eps_param0 = sqrt<U>(err_image * errors_image.back());
  bool flag = false;
  bool flag0 = false;
    diss_printf(verbose, fp,
		"eps_param0 = %s\n", tostring<U>(eps_param0).c_str());
  for (vector<int>::const_iterator it = kernel_dim.begin(); 
       it != kernel_dim.end(); ++it) {
    nn0 = *it;

    flag0 = check_kern<T, U>(nn0, n_dim, n_dim, a2.addrCoefs(),
			     permute_q,
			     dim_augkern, eps_machine,
			     eps_param0, sym_flag, errors);
      diss_printf(verbose, fp, "%d : %s / %s / %s\n",
	      nn0,
	      tostring<U>(errors[0]).c_str(), tostring<U>(errors[1]).c_str(),
	      tostring<U>(errors[2]).c_str());
    //    if (nn0 > 1) {
    if (!flag0 && (nn0 > 1)) {
      diss_printf(verbose, fp, 
		  "first trial by error from image %d %s -> %s fails\n",
		  n_dim, tostring<U>(errors_image.back()).c_str(),
		  tostring<U>(eps_param0).c_str());
      flag0 = check_kern<T, U>((nn0 - 1), n_dim, n_dim, a2.addrCoefs(),
			       permute_q,
			       dim_augkern, eps_machine,
			       eps_param0, sym_flag,
			       &errors[3]);
      diss_printf(verbose, fp, "%d : %s / %s / %s\n",
		  (nn0 - 1), tostring<U>(errors[3]).c_str(),
		  tostring<U>(errors[4]).c_str(),
		  tostring<U>(errors[5]).c_str());
      diss_printf(verbose, fp, "diff = %s\n",
		  tostring<U>(errors[0] - errors[4]).c_str());
      // criteria of equality
      U xtmp = ((fabs(errors[0]) > fabs(errors[4])) ? 
		     fabs(errors[0]) : fabs(errors[4]));
      xtmp = sqrt<U>(xtmp * eps_machine);
      U eps_param1 = (errors[3] + errors[4]) / U(2);
      eps_param1 = sqrt<U>(eps_param1 * err_image);
      eps_param0 = sqrt<U>(err_image * errors_image.back());
      diss_printf(verbose, fp, "%s + %s = %s / %s\n", 
		  tostring<U>(errors[3]).c_str(),
		  tostring<U>(errors[4]).c_str(),
		  tostring<U>(eps_param1).c_str(),
		  tostring<U>(eps_param0).c_str());
      if (errors[3] < eps_param0) {
	diss_printf(verbose, fp, "not satisfies the condition %s < %s\n", 
		    tostring<U>(eps_param0).c_str(),
		    tostring<U>(errors[3]).c_str());
	int itmp = n_dim - nn0;
	U err_image_tmp;

	err_image_tmp = check_matrixerr<T, U>(n_dim, n_dim,
					      a2.addrCoefs(), dim_augkern,
					      itmp,
					      permute_q,
					      eps_machine,
					      sym_flag);
	diss_printf(verbose, fp,
		    "part of %d (%d - %d) is regular, err_image=%s/%s\n",
		    itmp, n_dim, nn0,
		    tostring<U>(err_image_tmp).c_str(),
		    tostring<U>(err_image).c_str());
	if (err_image < err_image_tmp) {
	  dim_image = itmp;
	  count_err_image_updated++;
	  err_image = err_image_tmp;
	}
	continue;
      }
      if ((errors[0] > eps_param1) && (errors[1] < eps_param1)) {
	if ((fabs(errors[0] - errors[4]) < xtmp)) {
	  diss_printf(verbose, fp, "found with %d dim. %d updated\n", 
		      dim_image, count_err_image_updated);
	}
	else{
	  diss_printf(verbose, fp, 
		      "found with %d dim. %d updated, needs be refactorized\n",
		      dim_image, count_err_image_updated);
	}
	flag = true;
	break;
      }
    } // if (nn0 > 0)
    else {
      if (flag0) {
	diss_printf(verbose, fp, "found\n");
	flag = true;
	break;
      }
    }
  } // loop : it
  if (!flag) {
    diss_printf(verbose, fp, "kernel detection routine does not work : ");
    if (kernel_dim.size() == 1) {
      diss_printf(verbose, fp, "detection by Householder = %d\n", nn0);
      nn0 = 0;
      // n0 = kernel_dim.front();
    }
    else {
      diss_printf(verbose, fp, "unclear : ");
      nn0 = kernel_dim.front();
      for (vector<int>::const_iterator it = kernel_dim.begin();
	   it != kernel_dim.end(); it++) {
	diss_printf(verbose, fp, "%d ", (*it));
	if (((*it) + dim_augkern) != n_dim) {
	  nn0 = (*it);
	  diss_printf(verbose, fp, " / ");
	}
	else {
	  diss_printf(verbose, fp, "+ %d = %d /", dim_augkern, n_dim);
	}
      }
      diss_printf(verbose, fp, "\n");
      diss_printf(verbose, fp, "detection by Householder = %d\n", nn0);
    }
  }
  else {
    if ((kernel_dim.size() == 1) && (kernel_dim.front() != nn0)) {
      if (nn0 != 1) {
	diss_printf(verbose, fp, "strange_matrix\n");
	nn0 = 0;
	flag = false;
      }
    }
  }
  //  delete [] a2;
  delete [] errors;
  delete [] permute_q;
  *nn0_ = nn0; 
  return flag;
}

template
bool VerifyDimKernel<double, double>(int *nn0_, int n_dim, double* a_fact,
				     vector<int>& kernel_dim,
				     const bool sym_flag,
				     const int dim_augkern,
				     const double eps_machine,
				     const bool verbose,
				     FILE *fp);
template
bool VerifyDimKernel<complex<double>, double>(int *nn0_,
					      int n_dim, 
					      complex<double>* a_fact,
					      vector<int>& kernel_dim,
					      const bool sym_flag,
					      const int dim_augkern,
					      const double eps_machine,
					      const bool verbose,
					      FILE *fp);

template
bool VerifyDimKernel<float, float>(int *nn0_, int n_dim, float* a_fact,
				     vector<int>& kernel_dim,
				     const bool sym_flag,
				     const int dim_augkern,
				     const float eps_machine,
				     const bool verbose,
				     FILE *fp);
template
bool VerifyDimKernel<complex<float>, float>(int *nn0_,
					      int n_dim, 
					      complex<float>* a_fact,
					      vector<int>& kernel_dim,
					      const bool sym_flag,
					      const int dim_augkern,
					      const float eps_machine,
					      const bool verbose,
					      FILE *fp);

#ifndef NO_OCTRUPLE
template
bool VerifyDimKernel<quadruple, quadruple>(int *nn0_, int n_dim,
					   quadruple* a_fact,
					   vector<int>& kernel_dim,
					   const bool sym_flag,
					   const int dim_augkern,
					   const quadruple eps_machine,
					   const bool verbose,
					   FILE *fp);
template
bool VerifyDimKernel<complex<quadruple>, quadruple>(int *nn0_,
						    int n_dim, 
						    complex<quadruple>* a_fact,
						    vector<int>& kernel_dim,
						    const bool sym_flag,
						    const int dim_augkern,
						    const quadruple eps_machine,
						    const bool verbose,
						    FILE *fp);

#endif

