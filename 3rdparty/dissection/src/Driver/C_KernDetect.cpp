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

#include "Compiler/OptionLibrary.hpp"
#include "Driver/C_KernDetect.hpp"
#include "Algebra/ColumnMatrix.hpp"
#include "Algebra/VectorArray.hpp"
#include "Compiler/DissectionIO.hpp"

// T may be std::complex of U and W is in higher precision than T
template<typename T, typename U, typename W, typename Y>
bool check_kern(const int n0, const int lda, const int n, W *a_ini,
		 int *permute,
		 const int dim_augkern, const U &eps,
		 const U &eps_param, const bool flag_sym, U *errors,
		 const bool verbose, FILE *fp)
{
  bool flag;
  //  W *a_q, *a_fq, *proj, *nsp, *nsp2;
  //  W *v, *alpha;
  const W zero(0.0);
  const W one(1.0);
  const W none(-1.0);
  const U Uzero(0.0);
  
  ColumnMatrix<W> a_q(n, n); 
  ColumnMatrix<W> a_fq(n, n); 
  ColumnMatrix<W> proj(n, n); 
  ColumnMatrix<W> nsp(n, n); 
  ColumnMatrix<W> nsp2(n,n);
  VectorArray<W> v(n); 
  VectorArray<W> alpha(n);

  for (int j = 0; j < n; j++) {
    for (int i = 0; i < n; i++) {
      const int ij = permute[i] + permute[j] * lda;
      a_fq(i, j) = a_ini[ij];
    }                                             
  }
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
			        eps,
				flag_sym);
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
				  dim_augkern,
				  eps,
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
bool check_kern<double, double,
		quadruple, quadruple>(const int n0, const int lda, const int n,
				      quadruple *a_ini, int *permute,
				      const int dim_augkern, const double &eps,
				      const double &eps_param,
				      const bool flag_sym,
				      double *errors,
				      const bool verbose, FILE *fp);
template
bool check_kern<complex<double>,
		double,
		complex<quadruple>,
		quadruple>(const int n0,
			   const int lda,
			   const int n,
			   complex<quadruple> *a_ini,
			   int *permute,
			   const int dim_augkern,
			   const double &eps,
			   const double &eps_param,
			   const bool flag_sym,
			   double *errors,
			   const bool verbose, FILE *fp);

#ifndef NO_OCTRUPLE
template
bool check_kern<quadruple, quadruple,
		octruple, octruple>(const int n0, const int lda,
				    const int n,
				    octruple *a_ini, int *permute,
				    const int dim_augkern,
				    const quadruple &eps,
				    const quadruple &eps_param,
				    const bool flag_sym,
				    quadruple *errors,
				    const bool verbose, FILE *fp);
template
bool
check_kern<complex<quadruple>,
	   quadruple,
	   complex<octruple>,
	   octruple>(const int n0,
		     const int lda,
		     const int n,
		     complex<octruple> *a_ini,
		     int *permute,
		     const int dim_augkern,
		     const quadruple &eps,
		     const quadruple &eps_param,
		     const bool flag_sym,
		     quadruple *errors,
		     const bool verbose, FILE *fp);
#endif
//
// T may be std::complex of U and W is in higher precision than T
template<typename T, typename U, typename W, typename Y>
U check_matrixerr(const int lda, const int n,
		   W *a,
		   const int dim_augkern, const int k,
		   int *permute,
		   const U &eps,
		   const bool flag_sym)
  
{
  U error;
  const U Uzero(0.0);
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
double check_matrixerr<double, double,
		       quadruple, quadruple>(const int lda, const int n,
					     quadruple *a,
					     const int dim_augkern,
					     const int k, int *permute,
					     const double &eps,
					     const bool flag_sym);

template
double check_matrixerr<complex<double>,
		       double,
		       complex<quadruple>,
		       quadruple>(const int lda,
				  const int n,
				  complex<quadruple> *a,
				  const int dim_augkern,
				  const int k,
				  int *permute,
				  const double &eps,
				  const bool flag_sym);
#ifndef NO_OCTRUPLE
template
quadruple check_matrixerr<quadruple,
			  quadruple,
			  octruple,
			  octruple>(const int lda, const int n,
				    octruple *a,
				    const int dim_augkern,
				    const int k, int *permute,
				    const quadruple &eps,
				    const bool flag_sym);

template
quadruple
check_matrixerr<complex<quadruple>,
		quadruple,
		complex<octruple>,
		octruple>(const int lda,
			  const int n,
			  complex<octruple> *a,
			  const int dim_augkern,
			  const int k,
			  int *permute,
			  const quadruple &eps,
			  const bool flag_sym);
#endif
//
 
template<typename T, typename U>
void verify_kernels(const int n0, const int n, const T *a_ini,
		    const T *a_fact, const int lda,
		    const double eps, U *errors,
		    const bool verbose, FILE *fp)
{
  const T one(1.0);
  const T zero(0.0);
  const T none(-1.0);

  ColumnMatrix<T> a_p(n, n);
  const int n1 = n - n0;
  
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n ; j++) {
      // const int ij0 = permute[i] + permute[j] * n;
      const int ij0 = i + j * lda;
      a_p(i, j) = a_ini[ij0];         //      const int ij1 = i + j * n;
    }
  }
  //  compute [A_11^-1 A_12]
  //          [     -I     ]
  VectorArray<T> v(n);
  VectorArray<T> w(n);
  for (int j = n1; j < n; j++) {
    for (int i = 0; i < n1; i++) {
      //      const int ij0 = permute[i] + permute[j] * n;
      const int ij0 = i + j * lda;
      v[i] = a_ini[ij0];
    }
    full_fwbw_part<T>(n1, (T *)a_fact, lda, v.addrCoefs()); //cast : 22 May 2018
    for (int i = n1; i < n; i++) {
      v[i] = zero;
    }
    v[j] = none;
    blas_gemv<T>(CblasNoTrans, n, n, one, a_p.addrCoefs(), n,
		 v.addrCoefs(), 1, zero,
		 w.addrCoefs(), 1);
#if 0
    diss_printf(verbose, fp, "%s %d verify_kernels %d\n",
		__FILE__, __LINE__, n);
    for (int i = 0; i < n; i++) {
      diss_printf(verbose, fp, "%d %s %s\n", i, tostring<T>(v[i]).c_str(),
		  tostring<T>(w[i]).c_str());
    }
#endif
    errors[j - n1] = blas_l2norm<T, U>(n, w.addrCoefs(), 1);
  }
}

template
void verify_kernels<double, double>(const int n0, const int n,
				    const double *a_ini,
				    const double *a_fact,
				    const int lda,
				    const double eps,
				    double *errors,
				    const bool verbose, FILE *fp);

template
void verify_kernels<complex<double>, double>(const int n0, const int n,
					     const complex<double> *a_ini,
					     const complex<double> *a_fact,
					     const int lda,
					     const double eps,
					     double *errors,
					     const bool verbose, FILE *fp);
template
void verify_kernels<quadruple, quadruple>(const int n0, const int n,
					  const quadruple *a_ini,
					  const quadruple *a_fact,
					  const int lda,
					  const double eps,
					  quadruple *errors,
					  const bool verbose, FILE *fp);

template
void verify_kernels<complex<quadruple>,
		    quadruple>(const int n0, const int n,
			       const complex<quadruple> *a_ini,
			       const complex<quadruple> *a_fact,
			       const int lda,
			       const double eps,
			       quadruple *errors,
			       const bool verbose, FILE *fp);
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
void HouseholderVector<complex<octruple> >(int n,
					    complex<octruple> *x,
					    complex<octruple> *v,
					    complex<octruple> *gamma)
{
  HouseholderVector_complex<complex<octruple>, octruple>(n, x, v, gamma);
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


template<typename T, typename U, typename W, typename Y>
bool ComputeDimKernel_(int *n0, bool *flag_unsym_permute,
		       const T *a_, const int n, 
		      const bool sym_flag,
		      const int dim_augkern,
		      const U eps_machine, // for perturbation
		      const double eps_piv,
		      const bool verbose,
		      FILE *fp)
{
  const W zero(0.0);
  const Y Yzero(0.0);
  const Y Yone(1.0);
  const U Uzero(0.0);
  Y Yeps_machine;
  Yeps_machine = conv_prec<Y, U>(eps_machine);
  // dimension of the image of the matrix a is at least one
  int nn0, n1, n2;
  int n_dim = n + 1;
  int *permute = new int[n_dim];

  ColumnMatrix<W> aq0(n_dim, n_dim);
  ColumnMatrix<T> ad0(n, n);
  ColumnMatrix<W> a1(n_dim, n_dim);
  ColumnMatrix<W> aq_fact(n_dim, n_dim);
  ColumnMatrix<T> ad_fact(n, n);
  VectorArray<W> aa_diag(n_dim);
  //  VectorArray<Y> diag_scale(n_dim);
  VectorArray<double> rr(n_dim - 1);

  int *permute_d = new int[n_dim];
  int *permute_q = new int[n_dim];

  int *permute_right = new int[n_dim];
  int *permute_left = new int[n_dim];
  bool flag_tmp;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      a1(i, j) = conv_prec<W, T>(a_[i + j * n]);
    }
  }
  if (sym_flag) {
    flag_tmp = false;
    for (int i = 0; i < n; i++) {
      permute_left[i] = permute_right[i] = i;
    }
  }
  else {
    double pivot = 1.0;
    double fop;
    double epspiv = todouble<U>(eps_machine);
    ldu_full_permute<W, Y>(&nn0, 0, n, a1.addrCoefs(), n_dim, &pivot,
			   permute_right, permute_left,
			   epspiv, &fop);
    flag_tmp = false;
    for (int i = 0; i < n; i++) {
      if (permute_left[i] != permute_right[i]) {
	flag_tmp = true;
	break;
      }
    }
  }
  diss_printf(verbose, fp, "%s %d : %s pivots\n",
	      __FILE__, __LINE__,
	      flag_tmp ? "full" : "symmetric ");
  for (int i = 0; i < n; i++) {
    diss_printf(verbose, fp, "%d : %d %d %s\n",
		i, permute_left[i], permute_right[i],
		tostring<W>(a1(i, i)).c_str());
    //	      tostring<T>(tolower<W, T>(diag_scale[i])).c_str());
  }
  for (int j = 0; j < n; j++) {
    for (int i = 0; i < n; i++) {
      a1(i, j) = conv_prec<W, T>(a_[permute_left[i] + permute_right[j] * n]);
    }
  }
  diss_printf(verbose, fp, "%s %d : permuted Schur complement %d\n", 
	      __FILE__, __LINE__, n);
  for (int i = 0; i < n; i++) {
    diss_printf(verbose, fp, "%d ", i);
    for (int j = 0; j < n; j++) {
      diss_printf(verbose, fp, "%s ",
		  tostring<T>(conv_prec<T, W>(a1(i, j))).c_str());
    }
    diss_printf(verbose, fp, "\n");
  }
  for (int i = 0; i < n; i++) {
    W tmp = zero;
    for (int j = 0; j < n; j++) {
      tmp += a1(i, j) + (random_bool() ? Yeps_machine : Yzero);
    }
    a1(i, n) = tmp;
    a1(n, n) += tmp + (random_bool() ? Yeps_machine : Yzero);
  }
  for (int j = 0; j < n; j++) {
    W tmp = zero;
    for (int i = 0; i < n; i++) {
      tmp += a1(i, j) + (random_bool() ? Yeps_machine : Yzero);
    }
    a1(n, j) = tmp;
  }
  // column and row of inflated matrix are scaled
  Y diag_scalen = sqrt(Yone / blas_abs<W, Y>(a1(n, n)));
  for (int i = 0; i < n; i++) {
    a1(i, n) *= diag_scalen;
    a1(n, i) *= diag_scalen;
  }
  a1(n, n) *= (diag_scalen *  diag_scalen);
	
  diss_printf(verbose, fp, "%s %d : the original Schur complement %d\n", 
	      __FILE__, __LINE__, n);
  for (int i = 0; i < n; i++) {
    diss_printf(verbose, fp, "%d ", i);
    for (int j = 0; j < n; j++) {
      diss_printf(verbose, fp, "%s ",
		  tostring<T>(conv_prec<T, W>(a_[i + j * n])).c_str());
    }
    diss_printf(verbose, fp, "\n");
  }

  diss_printf(verbose, fp, "%s %d : the scaled last Schur complement %d\n", 
	      __FILE__, __LINE__, n_dim);
  for (int i = 0; i < n_dim; i++) {
    diss_printf(verbose, fp, "%d ", i);
    for (int j = 0; j < n_dim; j++) {
      diss_printf(verbose, fp, "%s ",
		  tostring<T>(conv_prec<T, W>(a1(i,j))).c_str());
    }
    diss_printf(verbose, fp, "\n");
  }
  aq0.copy(a1);
  aq_fact.copy(a1);
  n1 = hqr_pivot<W, Y>(n_dim, aq0.addrCoefs(), permute);
  diss_printf(verbose, fp,
	      "%s %d dimension of the image deteced by d_hqr_pivot() is %d\n",
	      __FILE__, __LINE__, n1);
  diss_printf(verbose, fp, "%s %d QR\n", __FILE__, __LINE__);
  for (int i = 0; i < n_dim; i++) {
    diss_printf(verbose, fp, "%d : %d ", i, permute[i]); 
    for (int j = 0; j < n_dim; j++) {
      diss_printf(verbose, fp, "%s ",
		  tostring<T>(conv_prec<T, W>(aq0(i, j))).c_str());
    }
    diss_printf(verbose, fp, "\n");
  }
  for (int i = 0; i < n1; i++) {
    aa_diag[i] = aq0(i, i);
  }
    
  list<int> pos_gap;
  vector<int> kernel_dim;
  pos_gap.push_back(dim_augkern);

  for (int i = 0; i < (n1 - 1); i++) {
    rr[i] = todouble<Y>(blas_abs<W, Y>(aa_diag[i + 1] / aa_diag[i]));
  }
  // find largest gap inside of invertible part
  double aug_diff = 1.0;
  for (int i = 0; i < (dim_augkern - 1); i++) {
    if (aug_diff > rr[i]) {
      aug_diff = rr[i];
    }
  }
  
  diss_printf(verbose, fp, "%s %d : aug_diff = %s esp = %.12e\n",
	      __FILE__, __LINE__,
	      tostring<double>(aug_diff).c_str(), eps_piv);
  for (int i = (dim_augkern - 1); i < (n1 - 1); i++) {
    if (rr[i] < (aug_diff * eps_piv)) {
      diss_printf(verbose, fp, "rr[%d] : %s\n",
		  i, tostring<double>(rr[i]).c_str());
      pos_gap.push_back(i + 1);
    }
  }
  const bool flag_perturb =
    ((rr[dim_augkern - 1] < aug_diff)
     && (blas_abs<W, Y>(aa_diag[n1 - 1]) > Yeps_machine * Y(eps_piv)));
    
  const U eps_perturb =  flag_perturb ? eps_machine : Uzero;
  diss_printf(verbose, fp, "%s %d : rr[%d] = %s => perturb = %s %s\n",
	      __FILE__, __LINE__,
	      (dim_augkern - 1),
	      tostring<double>(rr[dim_augkern - 1]).c_str(),
	      tostring<U>(eps_perturb).c_str(),
	      tostring<Y>(blas_abs<W, Y>(aa_diag[n1 - 1])).c_str()
	      );

  int n_dim1;
  if (flag_perturb) {  // condition number of the invertible part is moderate
    n2 = n1;           // and the diagonal entry becomes etreamly small
    // n1 may be less than n_dim when a diagnal of QR == 0
    for (int i = 0; i < (n1 - 1); i++) {
      if (rr[i] < (todouble<U>(eps_machine) / sqrt(eps_piv))) {
	n2 = i + 1;
	break;
      }
    }
    n_dim1 = n2 > dim_augkern ? n2 : n_dim;
    pos_gap.push_back(n2 - 1); // matrix being inflated 
  }
  else {
    n_dim1 = n_dim;
    pos_gap.push_back(n); // matrix being inflated : n = n_dim - 1
  }
  rr.free();
  pos_gap.sort();
  pos_gap.unique();

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

  int flag;
  flag = VerifyDimKernel<T, U, W, Y>(&nn0, permute_q,
				     n_dim1, aq_fact.addrCoefs(), kernel_dim,
				     sym_flag, dim_augkern, eps_perturb,
				     verbose, fp);
  if (flag == false) {
    nn0 = 1;           // iflated matrix is generated as singular
  }
  //  nn0 += (n_dim - n2) - 1;
  //  nn0--;   // detection of kernel is done for dim_augkern > 0 and at least 
  //  if (nn0 > 0) {
  {
    const double eps_machine_double = todouble<U>(eps_machine);
    U *errors_d = new U[nn0 + 1];
    list<int> nnn;
    if (nn0 > 2) {
      nnn.push_back(nn0 - 2);
    }
    nnn.push_back(nn0 - 1); //
    nnn.push_back(nn0);

    for (int j = 0; j < n; j++) {
      for (int i = 0; i < n; i++) {
	ad_fact(i, j) = a_[i + j * n];
      }
    }
    {
      double pivot_ref = 0.0;
      for (int i = 0; i < n ; i++) {
	double tmp = blas_abs<T, double>(a_[i + i * n]); 
	pivot_ref = pivot_ref > tmp ? pivot_ref : tmp;
      }
      int nnn0;
      double fop;
      if (sym_flag) {
	full_ldlt_permute<T, U>(&nnn0, 0, n, ad_fact.addrCoefs(), n,
				&pivot_ref, permute_d, eps_machine_double,
				&fop);
      }
      else {
	full_ldu_permute<T, U>(&nnn0, 0, n, ad_fact.addrCoefs(), n,
			       &pivot_ref, permute_d, eps_machine_double,
			       &fop);
      }
    }
    for (int j = 0; j < n; j++) {
      for (int i = 0; i < n; i++) {
	ad0(i, j) = a_[permute_d[i] + permute_d[j] * n];
      }
    }
    diss_printf(verbose, fp, "%s %d : verify errors : size = %d : ",
		__FILE__, __LINE__, nnn.size());
    for (list<int>::const_iterator it = nnn.begin(); it != nnn.end();
	 ++it) {
      diss_printf(verbose, fp, "%d ", (*it));
    }
    diss_printf(verbose, fp, "\n");
    for (list<int>::const_iterator it = nnn.begin(); it != nnn.end();
	 ++it) {
      const int nn = (*it);
      if (nn > 0) {
	verify_kernels<T, U>(nn, n, ad0.addrCoefs(), ad_fact.addrCoefs(),
			     n,
			     eps_machine_double, errors_d,
			     verbose, fp);
	diss_printf(verbose, fp, "errors of %d kernels\n", nn);
	for (int i = 0; i < nn; i++) {
	  diss_printf(verbose, fp,
		      "%d : %s\n", i, tostring<U>(errors_d[i]).c_str());
	}
      }
    }
    if (flag_tmp) {
      for (int j = 0; j < n; j++) {
	for (int i = 0; i < n; i++) {
	  ad_fact(i, j) = a_[permute_left[i] + permute_right[j] * n];
	}
      }
      {
	double pivot_ref = 0.0;
	for (int i = 0; i < n ; i++) {
	  double tmp = blas_abs<T, double>(a_[i + i * n]); 
	  pivot_ref = pivot_ref > tmp ? pivot_ref : tmp;
	}
	int nnn0;
	double fop;
	if (sym_flag) {
	  full_ldlt_permute<T, U>(&nnn0, 0, n, ad_fact.addrCoefs(), n,
				  &pivot_ref, permute_d, eps_machine_double,
				  &fop);
	}
	else {
	  full_ldu_permute<T, U>(&nnn0, 0, n, ad_fact.addrCoefs(), n,
				 &pivot_ref, permute_d, eps_machine_double,
				 &fop);
	}
      } 
      for (int j = 0; j < n; j++) {
	for (int i = 0; i < n; i++) {
	  ad0(i, j) = a_[permute_d[permute_left[i]] + permute_d[permute_right[j]] * n];
	}
      }
      diss_printf(verbose, fp, "%s %d : verify errors : size = %d : ",
		__FILE__, __LINE__, nnn.size());
      for (list<int>::const_iterator it = nnn.begin(); it != nnn.end();
	   ++it) {
	diss_printf(verbose, fp, "%d ", (*it));
      }
      diss_printf(verbose, fp, "\n");
      for (list<int>::const_iterator it = nnn.begin(); it != nnn.end();
	   ++it) {
	const int nn = (*it);
	if (nn > 0) {
	  verify_kernels<T, U>(nn, n, ad0.addrCoefs(), ad_fact.addrCoefs(), n,
			       eps_machine_double, errors_d,
			     verbose, fp);
	  diss_printf(verbose, fp, "errors of %d kernels\n", nn);
	  for (int i = 0; i < nn; i++) {
	    diss_printf(verbose, fp,
			"%d : %s\n", i, tostring<U>(errors_d[i]).c_str());
	  }
	}
      }
    } // if (flag_tmp) 
    //    nnn.clear();
    Y *errors_q = new Y[nn0 + 1];

    for (int j = 0; j < n_dim; j++) {
      for (int i = 0; i < n_dim; i++) {
	aq0(i, j) = a1(permute_q[i], permute_q[j]);
      }
    }

    for (list<int>::const_iterator it = nnn.begin(); it != nnn.end();
	 ++it) {
      const int nn = (*it) + 1;
      verify_kernels<W, Y>(nn, n_dim,
			   aq0.addrCoefs(), aq_fact.addrCoefs(), n_dim,
			   eps_machine_double, errors_q,
			   verbose, fp);
      diss_printf(verbose, fp, "errors of %d kernels\n", nn);
      for (int i = 0; i < nn; i++) {
	diss_printf(verbose, fp,
		    "%d :  %s\n", i,
		    tostring<Y>(errors_q[i]).c_str());
      }
    }

    
    for (int j = 0; j < n; j++) {
      for (int i = 0; i < n; i++) {
	aq_fact(i, j) = conv_prec<W, T>(a_[i + j * n]);
      }
    }
    {
      double pivot_ref = 0.0;
      for (int i = 0; i < n ; i++) {
	double tmp = blas_abs<T, double>(a_[i + i * n]); 
	pivot_ref = pivot_ref > tmp ? pivot_ref : tmp;
      }
      int nnn0;
      double fop;
      if (sym_flag) {
	full_ldlt_permute<W, Y>(&nnn0, 0, n, aq_fact.addrCoefs(), n_dim,
				&pivot_ref, permute_q, eps_machine_double,
				&fop);
      }
      else {
	full_ldu_permute<W, Y>(&nnn0, 0, n, aq_fact.addrCoefs(), n_dim,
			       &pivot_ref, permute_q, eps_machine_double,
			       &fop);
      }
    }
    for (int j = 0; j < n; j++) {
      for (int i = 0; i < n; i++) {
	aq0(i, j) = conv_prec<W, T>(a_[permute_q[i] + permute_q[j] * n]);
      }
    }
    diss_printf(verbose, fp, "%s %d : verify errors : size = %d : ",
		__FILE__, __LINE__, nnn.size());
    for (list<int>::const_iterator it = nnn.begin(); it != nnn.end();
	 ++it) {
      diss_printf(verbose, fp, "%d ", (*it));
    }
    diss_printf(verbose, fp, "\n");
    for (list<int>::const_iterator it = nnn.begin(); it != nnn.end();
	 ++it) {
      const int nn = (*it);
      if (nn > 0) {
	verify_kernels<W, Y>(nn, n,
			     aq0.addrCoefs(), aq_fact.addrCoefs(), n_dim,
			     eps_machine_double, errors_q,
			     verbose, fp);
	diss_printf(verbose, fp, "errors of %d kernels\n", nn);
	for (int i = 0; i < nn; i++) {
	  diss_printf(verbose, fp,
		      "%d : %s\n", i, tostring<Y>(errors_q[i]).c_str());
	}
      }
    }

    nnn.clear();
    delete [] errors_q;
    delete [] errors_d;
  } // if (nn0 > 0)
  delete [] permute;
  delete [] permute_d;
  delete [] permute_q;
  delete [] permute_left;
  delete [] permute_right;
  *n0 = nn0 + (n_dim - n_dim1) - 1;
  diss_printf(verbose, fp,
	      "%s %d : dimension of the kernel is %d\n",
	      __FILE__, __LINE__, *n0);


  
  *flag_unsym_permute = flag_tmp;

  return flag;
}

template
bool ComputeDimKernel_<double, double,
		       quadruple, quadruple>(int *n0, bool *flag_unsym_permute,
					     const double *a_,
					     const int n, 
					     const bool sym_flag,
					     const int dim_augkern,
					     const double eps_machine,
					     const double eps_piv,
					     const bool verbose,
					     FILE *fp);

template
bool ComputeDimKernel_<complex<double>, 
		       double,
		       complex<quadruple>,
		       quadruple>(int *n0, bool *flag_unsym_permute,
			      const complex<double> *a_,
			      const int n, 
			      const bool sym_flag,
			      const int dim_augkern,
			      const double eps_machine,
			      const double eps_piv,
			      const bool verbose,
			      FILE *fp);

template
bool ComputeDimKernel_<quadruple, 
		       quadruple,
		       octruple,
		       octruple>(int *n0, bool *flag_unsym_permute,
				 const quadruple *a_,
				 const int n, 
				 const bool sym_flag,
				 const int dim_augkern,
				 const quadruple eps_machine,
				 const double eps_piv,
				 const bool verbose,
				 FILE *fp);

template
bool ComputeDimKernel_<complex<quadruple>,
		       quadruple,
		       complex<octruple>,
		       octruple>(int *n0, bool *flag_unsym_permute,
				 const complex<quadruple> *a_,
				 const int n, 
				 const bool sym_flag,
				 const int dim_augkern,
				 const quadruple eps_machine,
				 const double eps_piv,
				 const bool verbose,
				 FILE *fp);


template<typename T, typename U>
bool ComputeDimKernel(int *n0, bool *flag_unsym_permute,
		      const T *a_, const int n, 
		      const bool sym_flag,
		      const int dim_augkern,
		      const U eps_machine, // for perturbation
		      const double eps_piv,
		      const bool verbose,
		      FILE *fp)
{
  fprintf(stderr, "%s %d : specialized template not implemented\n",
	  __FILE__, __LINE__);
}

template<>
bool ComputeDimKernel<double, double>(int *n0, bool *flag_unsym_permute,
				      const double *a_,
				      const int n, 
				      const bool sym_flag,
				      const int dim_augkern,
				      const double eps_machine,
				      const double eps_piv,
				      const bool verbose,
				      FILE *fp)
{
  return ComputeDimKernel_<double,
			   double,
			   quadruple,
			   quadruple>(n0, flag_unsym_permute, a_, n,
				      sym_flag, dim_augkern,
				      eps_machine, eps_piv,
				      verbose, fp);
}

template<>
bool ComputeDimKernel<complex<double>, 
		      double>(int *n0, bool *flag_unsym_permute,
			      const complex<double> *a_,
			      const int n, 
			      const bool sym_flag,
			      const int dim_augkern,
			      const double eps_machine,
			      const double eps_piv,
			      const bool verbose,
			      FILE *fp)
{
  return ComputeDimKernel_<complex<double>,
			   double,
			   complex<quadruple>,
			   quadruple>(n0, flag_unsym_permute, a_, n,
				      sym_flag, dim_augkern,
				      eps_machine, eps_piv,
				      verbose, fp);
  
}

template<>
bool ComputeDimKernel<quadruple, 
		      quadruple>(int *n0, bool *flag_unsym_permute,
				 const quadruple *a_,
				 const int n, 
				 const bool sym_flag,
				 const int dim_augkern,
				 const quadruple eps_machine,
				 const double eps_piv,
				 const bool verbose,
				 FILE *fp)
{
  return ComputeDimKernel_<quadruple,
			   quadruple,
			   octruple,
			   octruple>(n0, flag_unsym_permute, a_, n,
				     sym_flag, dim_augkern,
				     eps_machine, eps_piv,
				     verbose, fp);
}

template<>
bool ComputeDimKernel<complex<quadruple>,
		      quadruple>(int *n0, bool *flag_unsym_permute,
				 const complex<quadruple> *a_,
				 const int n, 
				 const bool sym_flag,
				 const int dim_augkern,
				 const quadruple eps_machine,
				 const double eps_piv,
				 const bool verbose,
				 FILE *fp)
{
  return ComputeDimKernel_<complex<quadruple>,
			   quadruple,
			   complex<octruple>,
			   octruple>(n0, flag_unsym_permute, a_, n,
				     sym_flag, dim_augkern,
				     eps_machine, eps_piv,
				     verbose, fp);
}

template<>
bool ComputeDimKernel<float, float>(int *n0, bool *flag_unsym_permute,
				      const float *a_,
				      const int n, 
				      const bool sym_flag,
				      const int dim_augkern,
				      const float eps_machine,
				      const double eps_piv,
				      const bool verbose,
				      FILE *fp)
{
  return ComputeDimKernel_<float,
			   float,
			   double,
			   double>(n0, flag_unsym_permute, a_, n,
				      sym_flag, dim_augkern,
				      eps_machine, eps_piv,
				      verbose, fp);
}

template<>
bool ComputeDimKernel<complex<float>, 
		      float>(int *n0, bool *flag_unsym_permute,
			      const complex<float> *a_,
			      const int n, 
			      const bool sym_flag,
			      const int dim_augkern,
			      const float eps_machine,
			      const double eps_piv,
			      const bool verbose,
			      FILE *fp)
{
  return ComputeDimKernel_<complex<float>,
			   float,
			   complex<double>,
			   double>(n0, flag_unsym_permute, a_, n,
				      sym_flag, dim_augkern,
				      eps_machine, eps_piv,
				      verbose, fp);
  
}


//

template<typename T, typename U, typename W, typename Y>
bool VerifyDimKernel(int *nn0_,
		     int *permute_q,
		     int n_dim, W* a_fact,
		     vector<int> &kernel_dim,
		     const bool sym_flag,
		     const int dim_augkern,
		     const U eps_machine,
		     const bool verbose,
		     FILE *fp)
{
  U Uzero(0.0);
  U Uone(1.0);
  int nn, nn0;
  int k;
  double pivot_ref, pivot_ref_q;
  //  int *permute_q = new int[n_dim];
  ColumnMatrix<W> a2(n_dim, n_dim);
  U *errors = new U[6];

  blas_copy<W>((n_dim * n_dim), a_fact, 1, a2.addrCoefs(), 1);
  //  for (int i = 0; i < (n_dim * n_dim); i++) {
  //    a2[i] = a_fact[i];
  //  }
  pivot_ref = 0.0;
  k = 0;
  for (int i = 0; i < n_dim; i++) {
    double dtmp = blas_abs<W, double>(a_fact[k]); // accuracy? : 14 Jul.2015 Atsushi
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
    full_ldlt_permute<W, Y>(&nn, nn1, n_dim, a_fact, n_dim, &pivot_ref_q,
			    permute_q, eps0, &fop);
  }
  else {
    double fop;
    const double eps0 = todouble<U>(eps_machine);
    full_ldu_permute<W, Y>(&nn, nn1, n_dim, a_fact, n_dim, &pivot_ref_q,
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
  U error_max, error_min, Ueps;
  error_min = Uone / eps_machine;
  error_max = Uzero;
  Ueps = machine_epsilon<U, U>();
  for (int i = 0; i < dims.size(); i++) {
    errors_image[i] = check_matrixerr<T, U, W, Y>(n_dim, n_dim,
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

    flag0 = check_kern<T, U, W, Y>(nn0, n_dim, n_dim, a2.addrCoefs(),
				   permute_q,
				   dim_augkern, eps_machine,
				   eps_param0, sym_flag, errors,
				   verbose, fp);
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
      flag0 = check_kern<T, U, W, Y>((nn0 - 1), n_dim, n_dim, a2.addrCoefs(),
				     permute_q,
				     dim_augkern, eps_machine,
				     eps_param0, sym_flag,
				     &errors[3], verbose, fp);
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
	
	err_image_tmp = check_matrixerr<T, U, W, Y>(n_dim, n_dim,
						    a2.addrCoefs(),
						    dim_augkern,
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
  //  delete [] permute_q;
  *nn0_ = nn0; 
  return flag;
}

template
bool VerifyDimKernel<double, double,
		     quadruple, quadruple>(int *nn0_,
					   int *permute_q,
					   int n_dim,
					   quadruple* a_fact,
					   vector<int>& kernel_dim,
					   const bool sym_flag,
					   const int dim_augkern,
					   const double eps_machine,
					   const bool verbose,
					   FILE *fp);
template
bool VerifyDimKernel<complex<double>, double,
		     complex<quadruple>, quadruple>(int *nn0_,
						    int *permute_q,
						    int n_dim, 
						    complex<quadruple>* a_fact,
						    vector<int>& kernel_dim,
						    const bool sym_flag,
						    const int dim_augkern,
						    const double eps_machine,
						    const bool verbose,
						    FILE *fp);
#ifndef NO_OCTRUPLE
template
bool VerifyDimKernel<quadruple, quadruple,
		     octruple, octruple>(int *nn0_,
					 int *permute_q,
					 int n_dim,
					 octruple* a_fact,
					 vector<int>& kernel_dim,
					 const bool sym_flag,
					 const int dim_augkern,
					 const quadruple eps_machine,
					 const bool verbose,
					 FILE *fp);
template
bool VerifyDimKernel<complex<quadruple>, quadruple,
		     complex<octruple>, octruple>(int *nn0_,
						  int *permute_q,
						  int n_dim, 
						  complex<octruple>* a_fact,
						  vector<int>& kernel_dim,
						  const bool sym_flag,
						  const int dim_augkern,
						  const quadruple eps_machine,
						  const bool verbose,
						  FILE *fp);

#endif

