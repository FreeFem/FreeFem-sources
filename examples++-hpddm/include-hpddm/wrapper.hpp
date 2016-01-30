/*
   This file is part of HPDDM.

   Author(s): Pierre Jolivet <pierre.jolivet@enseeiht.fr>
        Date: 2014-08-04

   Copyright (C) 2011-2014 Université de Grenoble
                 2015      Eidgenössische Technische Hochschule Zürich
                 2016-     Centre National de la Recherche Scientifique

   HPDDM is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published
   by the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   HPDDM is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with HPDDM.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _HPDDM_WRAPPER_
#define _HPDDM_WRAPPER_

#define HPDDM_GENERATE_EXTERN_MKL(C, T)                                                                      \
void cblas_ ## C ## gthr(const int, const T*, T*, const int*);                                               \
void cblas_ ## C ## sctr(const int, const T*, const int*, T*);

#if HPDDM_MKL && !defined(INTEL_MKL_VERSION)
extern "C" {
HPDDM_GENERATE_EXTERN_MKL(s, float)
HPDDM_GENERATE_EXTERN_MKL(d, double)
HPDDM_GENERATE_EXTERN_MKL(c, std::complex<float>)
HPDDM_GENERATE_EXTERN_MKL(z, std::complex<double>)
}
#endif // HPDDM_MKL && !defined(INTEL_MKL_VERSION)

namespace HPDDM {
/* Class: Wrapper
 *
 *  A class for handling dense and sparse linear algebra.
 *
 * Template Parameter:
 *    K              - Scalar type. */
template<class K>
struct Wrapper {
    /* Function: mpi_type
     *  Returns the MPI datatype of the template parameter of <Wrapper>. */
    static MPI_Datatype mpi_type();
    /* Function: mpi_underlying_type
     *  Returns the MPI datatype of the underlying type of the template parameter of <Wrapper>. */
    static MPI_Datatype mpi_underlying_type() {
        return Wrapper<underlying_type<K>>::mpi_type();
    }
    static constexpr bool is_complex = !std::is_same<K, underlying_type<K>>::value;
    /* Variable: transc
     *  Transposed real operators or Hermitian transposed complex operators. */
    static constexpr char transc = is_complex ? 'C' : 'T';
    /* Variable: I
     *  Numbering of a sparse <MatrixCSR>. */
#if HPDDM_MKL
    static constexpr char I = 'F';
#else
    static constexpr char I = HPDDM_NUMBERING;
#endif
    /* Variable: d__0
     *  Zero. */
    static constexpr K d__0 = K(0.0);
    /* Variable: d__1
     *  One. */
    static constexpr K d__1 = K(1.0);
    /* Variable: d__2
     *  Minus one. */
    static constexpr K d__2 = K(-1.0);

    /* Function: csrmv(square)
     *  Computes a sparse square matrix-vector product. */
    template<char N = HPDDM_NUMBERING>
    static void csrmv(bool, const int* const, const K* const, const int* const, const int* const, const K* const, K* const);
    /* Function: csrmv
     *  Computes a scalar-sparse matrix-vector product. */
    template<char N = HPDDM_NUMBERING>
    static void csrmv(const char* const, const int* const, const int* const, const K* const, bool,
                      const K* const, const int* const, const int* const, const K* const, const K* const, K* const);
    /* Function: csrmm(square)
     *  Computes a sparse square matrix-matrix product. */
    template<char N = HPDDM_NUMBERING>
    static void csrmm(bool, const int* const, const int* const, const K* const, const int* const, const int* const, const K* const, K* const);
    /* Function: csrmm
     *  Computes a scalar-sparse matrix-matrix product. */
    template<char N = HPDDM_NUMBERING>
    static void csrmm(const char* const, const int* const, const int* const, const int* const, const K* const, bool,
                      const K* const, const int* const, const int* const, const K* const, const int* const,
                      const K* const, K* const, const int* const);

    /* Function: csrcsc
     *  Converts a matrix stored in Compressed Sparse Row format into Compressed Sparse Column format. */
    template<char, char>
    static void csrcsc(const int* const, const K* const, const int* const, const int* const, K* const, int* const, int* const);
    /* Function: gthr
     *  Gathers the elements of a full-storage sparse vector into compressed form. */
    static void gthr(const int&, const K* const, K* const, const int* const);
    /* Function: sctr
     *  Scatters the elements of a compressed sparse vector into full-storage form. */
    static void sctr(const int&, const K* const, const int* const, K* const);
    /* Function: diag(in-place)
     *  Computes a vector-vector element-wise multiplication. */
    static void diag(const int&, const underlying_type<K>* const, K* const);
    /* Function: diag
     *  Computes a vector-vector element-wise multiplication. */
    static void diag(const int&, const underlying_type<K>* const, const K* const, K* const);
    /* Function: diag(in-place)
     *  Computes a vector-matrix element-wise multiplication. */
    static void diag(const int&, const int&, const underlying_type<K>* const, K* const);
    /* Function: diag
     *  Computes a vector-matrix element-wise multiplication. */
    static void diag(const int&, const int&, const underlying_type<K>* const, const K* const, K* const);
    /* Function: conj
     *  Conjugates a real or complex number. */
    template<class T, typename std::enable_if<!Wrapper<T>::is_complex>::type* = nullptr>
    static T conj(T& x) { return x; }
    template<class T, typename std::enable_if<Wrapper<T>::is_complex>::type* = nullptr>
    static T conj(T& x) { return std::conj(x); }
    /* Function: imatcopy
     *  Transforms (copy, transpose, conjugate transpose, conjugate) a dense matrix in-place. */
    template<char O>
    static void imatcopy(const int, const int, K* const, const int, const int);
    /* Function: omatcopy
     *  Transforms (copy, transpose, conjugate transpose, conjugate) a dense matrix out-place. */
    template<char O>
    static void omatcopy(const int, const int, const K* const, const int, K* const, const int);
};

template<>
inline MPI_Datatype Wrapper<float>::mpi_type() { return MPI_FLOAT; }
template<>
inline MPI_Datatype Wrapper<double>::mpi_type() { return MPI_DOUBLE; }
template<>
inline MPI_Datatype Wrapper<std::complex<float>>::mpi_type() { return MPI_COMPLEX; }
template<>
inline MPI_Datatype Wrapper<std::complex<double>>::mpi_type() { return MPI_DOUBLE_COMPLEX; }

template<class K>
constexpr char Wrapper<K>::transc;

template<class K>
constexpr K Wrapper<K>::d__0;
template<class K>
constexpr K Wrapper<K>::d__1;
template<class K>
constexpr K Wrapper<K>::d__2;

template<class K>
inline void Wrapper<K>::diag(const int& n, const underlying_type<K>* const d, K* const in) {
    diag(n, d, nullptr, in);
}
template<class K>
inline void Wrapper<K>::diag(const int& m, const int& n, const underlying_type<K>* const d, K* const in) {
    diag(m, n, d, nullptr, in);
}

#if HPDDM_MKL
template<char N>
struct matdescr {
    static const char a[];
    static const char b[];
};

template<char N>
const char matdescr<N>::a[6] { 'G', '0', '0', N, '0', '0' };
template<char N>
const char matdescr<N>::b[6] { 'S', 'L', 'N', N, '0', '0' };

#define HPDDM_GENERATE_MKL(C, T)                                                                             \
template<>                                                                                                   \
template<char N>                                                                                             \
inline void Wrapper<T>::csrmv(bool sym, const int* const n, const T* const a, const int* const ia,           \
                              const int* const ja, const T* const x, T* const y) {                           \
    if(N == 'C') {                                                                                           \
        if(sym)                                                                                              \
            mkl_cspblas_ ## C ## csrsymv("L", HPDDM_CONST(int, n), HPDDM_CONST(T, a), HPDDM_CONST(int, ia),  \
                                         HPDDM_CONST(int, ja), HPDDM_CONST(T, x), y);                        \
        else                                                                                                 \
            mkl_cspblas_ ## C ## csrgemv("N", HPDDM_CONST(int, n), HPDDM_CONST(T, a), HPDDM_CONST(int, ia),  \
                                         HPDDM_CONST(int, ja), HPDDM_CONST(T, x), y);                        \
    }                                                                                                        \
    else {                                                                                                   \
        if(sym)                                                                                              \
            mkl_ ## C ## csrsymv("L", HPDDM_CONST(int, n), HPDDM_CONST(T, a), HPDDM_CONST(int, ia),          \
                                 HPDDM_CONST(int, ja), HPDDM_CONST(T, x), y);                                \
        else                                                                                                 \
            mkl_ ## C ## csrgemv("N", HPDDM_CONST(int, n), HPDDM_CONST(T, a), HPDDM_CONST(int, ia),          \
                                 HPDDM_CONST(int, ja), HPDDM_CONST(T, x), y);                                \
    }                                                                                                        \
}                                                                                                            \
template<>                                                                                                   \
template<char N>                                                                                             \
inline void Wrapper<T>::csrmv(const char* const trans, const int* const m, const int* const k,               \
                              const T* const alpha, bool sym, const T* const a, const int* const ia,         \
                              const int* const ja, const T* const x, const T* const beta, T* const y) {      \
    mkl_ ## C ## csrmv(HPDDM_CONST(char, trans), HPDDM_CONST(int, m), HPDDM_CONST(int, k),                   \
                       HPDDM_CONST(T, alpha), HPDDM_CONST(char, sym ? matdescr<N>::b : matdescr<N>::a),      \
                       HPDDM_CONST(T, a), HPDDM_CONST(int, ja), HPDDM_CONST(int, ia),                        \
                       HPDDM_CONST(int, ia) + 1, HPDDM_CONST(T, x), HPDDM_CONST(T, beta), y);                \
}                                                                                                            \
template<>                                                                                                   \
template<char N>                                                                                             \
inline void Wrapper<T>::csrmm(const char* const trans, const int* const m, const int* const n,               \
                              const int* const k, const T* const alpha, bool sym,                            \
                              const T* const a, const int* const ia, const int* const ja,                    \
                              const T* const x, const int* const ldb, const T* const beta,                   \
                              T* const y, const int* const ldc) {                                            \
    if(*n != 1) {                                                                                            \
        if(N != 'F') {                                                                                       \
            std::for_each(const_cast<int*>(ja), const_cast<int*>(ja) + ia[*m], [](int& i) { ++i; });         \
            std::for_each(const_cast<int*>(ia), const_cast<int*>(ia) + *m + 1, [](int& i) { ++i; });         \
        }                                                                                                    \
        mkl_ ## C ## csrmm(HPDDM_CONST(char, trans), HPDDM_CONST(int, m), HPDDM_CONST(int, n),               \
                           HPDDM_CONST(int, k), HPDDM_CONST(T, alpha),                                       \
                           HPDDM_CONST(char, sym ? matdescr<'F'>::b : matdescr<'F'>::a), HPDDM_CONST(T, a),  \
                           HPDDM_CONST(int, ja), HPDDM_CONST(int, ia), HPDDM_CONST(int, ia) + 1,             \
                           HPDDM_CONST(T, x), HPDDM_CONST(int, ldb), HPDDM_CONST(T, beta),  y,               \
                           HPDDM_CONST(int, ldc));                                                           \
        if(N != 'F') {                                                                                       \
            std::for_each(const_cast<int*>(ia), const_cast<int*>(ia) + *m + 1, [](int& i) { --i; });         \
            std::for_each(const_cast<int*>(ja), const_cast<int*>(ja) + ia[*m], [](int& i) { --i; });         \
        }                                                                                                    \
    }                                                                                                        \
    else                                                                                                     \
        csrmv<N>(trans, m, k, alpha, sym, a, ia, ja, x, beta, y);                                            \
}                                                                                                            \
                                                                                                             \
template<>                                                                                                   \
template<char N, char M>                                                                                     \
inline void Wrapper<T>::csrcsc(const int* const n, const T* const a, const int* const ja,                    \
                               const int* const ia, T* const b, int* const jb, int* const ib) {              \
    int job[6] { 0, N == 'F', M == 'F', 0, 0, 1 };                                                           \
    int error;                                                                                               \
    mkl_ ## C ## csrcsc(job, HPDDM_CONST(int, n), const_cast<T*>(a), const_cast<int*>(ja),                   \
                        const_cast<int*>(ia), b, jb, ib, &error);                                            \
}                                                                                                            \
template<>                                                                                                   \
inline void Wrapper<T>::gthr(const int& n, const T* const y, T* const x, const int* const indx) {            \
    cblas_ ## C ## gthr(n, y, x, indx);                                                                      \
}                                                                                                            \
template<>                                                                                                   \
inline void Wrapper<T>::sctr(const int& n, const T* const x, const int* const indx, T* const y) {            \
    cblas_ ## C ## sctr(n, x, indx, y);                                                                      \
}                                                                                                            \
template<>                                                                                                   \
template<char O>                                                                                             \
inline void Wrapper<T>::imatcopy(const int n, const int m, T* const ab, const int lda, const int ldb) {      \
    static_assert(O == 'N' || O == 'R' || O == 'T' || O == 'C', "Unknown operation");                        \
    mkl_ ## C ## imatcopy('C', O, m, n, d__1, ab, lda, ldb);                                                 \
}                                                                                                            \
template<>                                                                                                   \
template<char O>                                                                                             \
inline void Wrapper<T>::omatcopy(const int n, const int m, const T* const a, const int lda,                  \
                                 T* const b, const int ldb) {                                                \
    static_assert(O == 'N' || O == 'R' || O == 'T' || O == 'C', "Unknown operation");                        \
    mkl_ ## C ## omatcopy('C', O, m, n, d__1, a, lda, b, ldb);                                               \
}
#define HPDDM_GENERATE_MKL_VML(C, T)                                                                         \
template<>                                                                                                   \
inline void Wrapper<T>::diag(const int& m, const int& n, const T* const d,                                   \
                             const T* const in, T* const out) {                                              \
    if(in)                                                                                                   \
        for(int i = 0; i < n; ++i)                                                                           \
            v ## C ## Mul(m, d, in + i * m, out + i * m);                                                    \
    else                                                                                                     \
        for(int i = 0; i < n; ++i)                                                                           \
            v ## C ## Mul(m, d, out + i * m, out + i * m);                                                   \
}                                                                                                            \
template<>                                                                                                   \
inline void Wrapper<T>::diag(const int& n, const T* const d, const T* const in, T* const out) {              \
    diag(n, i__1, d, in, out);                                                                               \
}
HPDDM_GENERATE_MKL(s, float)
HPDDM_GENERATE_MKL(d, double)
HPDDM_GENERATE_MKL(c, std::complex<float>)
HPDDM_GENERATE_MKL(z, std::complex<double>)
HPDDM_GENERATE_MKL_VML(s, float)
HPDDM_GENERATE_MKL_VML(d, double)
#else
template<class K>
template<char N>
inline void Wrapper<K>::csrmv(bool sym, const int* const n, const K* const a, const int* const ia, const int* const ja, const K* const x, K* const y) {
    csrmv<N>("N", n, n, &d__1, sym, a, ia, ja, x, &d__0, y);
}
template<class K>
template<char N>
inline void Wrapper<K>::csrmv(const char* const trans, const int* const m, const int* const k, const K* const alpha, bool sym,
                              const K* const a, const int* const ia, const int* const ja, const K* const x, const K* const beta, K* const y) {
    if(*trans == 'N') {
        if(sym) {
            if(beta == &d__0)
                std::fill_n(y, *m, K());
            else if(beta != &d__1)
                Blas<K>::scal(m, beta, y, &i__1);
            for(int i = 0; i < *m; ++i) {
                if(ia[i + 1] != ia[i]) {
                    K res = K();
                    int l = ia[i] - (N == 'F');
                    int j = ja[l] - (N == 'F');
                    while(l < ia[i + 1] - 1 - (N == 'F')) {
                        res += a[l] * x[j];
                        y[j] += *alpha * a[l] * x[i];
                        j = ja[++l] - (N == 'F');
                    }
                    if(i != j) {
                        res += a[l] * x[j];
                        y[j] += *alpha * a[l] * x[i];
                        y[i] += *alpha * res;
                    }
                    else
                        y[i] += *alpha * (res + a[l] * x[i]);
                }
            }
        }
        else {
            if(beta == &d__0)
                std::fill_n(y, *m, K());
#pragma omp parallel for schedule(static, HPDDM_GRANULARITY)
            for(int i = 0; i < *m; ++i) {
                K res = K();
                for(int l = ia[i] - (N == 'F'); l < ia[i + 1] - (N == 'F'); ++l)
                    res += a[l] * x[ja[l] - (N == 'F')];
                y[i] = *alpha * res + *beta * y[i];
            }
        }
    }
    else {
        if(beta == &d__0)
            std::fill_n(y, *k, K());
        else if(beta != &d__1)
            Blas<K>::scal(k, beta, y, &i__1);
        if(sym) {
            for(int i = 0; i < *m; ++i) {
                K res = K();
                for(int l = ia[i] - (N == 'F'); l < ia[i + 1] - (N == 'F'); ++l) {
                    int j = ja[l] - (N == 'F');
                    y[j] += *alpha * a[l] * x[i];
                    if(i != j)
                        res += a[l] * x[j];
                }
                y[i] += *alpha * res;
            }
        }
        else {
            for(int i = 0; i < *m; ++i)
                for(int l = ia[i] - (N == 'F'); l < ia[i + 1] - (N == 'F'); ++l)
                    y[ja[l] - (N == 'F')] += *alpha * a[l] * x[i];
        }
    }
}
template<class K>
template<char N>
inline void Wrapper<K>::csrmm(const char* const trans, const int* const m, const int* const n, const int* const k, const K* const alpha, bool sym,
                              const K* const a, const int* const ia, const int* const ja, const K* const x, const int* const ldb, const K* const beta, K* const y, const int* const ldc) {
    if(*trans == 'N') {
        if(*ldb != *k || *ldc != *m)
            return;
        int dimY = *m;
        K* res;
        if(sym) {
            int j;
            int dimX = *k;
            int dimNY = dimY * *n;
            if(beta == &d__0)
                std::fill_n(y, dimNY, K());
            else if(beta != &d__1)
                Blas<K>::scal(&dimNY, beta, y, &i__1);
            res = new K[*n];
            for(int i = 0; i < dimY; ++i) {
                std::fill_n(res, *n, K());
                for(int l = ia[i] - (N == 'F'); l < ia[i + 1] - (N == 'F'); ++l) {
                    j = ja[l] - (N == 'F');
                    if(i != j)
                        for(int r = 0; r < *n; ++r) {
                            res[r] += a[l] * x[j + r * dimX];
                            y[j + r * *ldb] += *alpha * a[l] * x[i + r * *ldc];
                        }
                    else
                        Blas<K>::axpy(n, a + l, x + j, k, res, &i__1);
                }
                Blas<K>::axpy(n, alpha, res, &i__1, y + i, m);
            }
            delete [] res;
        }
        else {
#pragma omp parallel private(res)
            {
                res = new K[*n];
#pragma omp for schedule(static, HPDDM_GRANULARITY)
                for(int i = 0; i < dimY; ++i) {
                    std::fill_n(res, *n, K());
                    for(int l = ia[i] - (N == 'F'); l < ia[i + 1] - (N == 'F'); ++l)
                        Blas<K>::axpy(n, a + l, x + ja[l] - (N == 'F'), k, res, &i__1);
                    Blas<K>::axpby(*n, *alpha, res, 1, *beta, y + i, dimY);
                }
                delete [] res;
            }
        }
    }
    else {
        if(*ldb != *m || *ldc != *k)
            return;
        int dimX = *m;
        int dimY = *k;
        int dimNY = dimY * *n;
        if(beta == &d__0)
            std::fill_n(y, dimNY, K());
        else if(beta != &d__1)
            Blas<K>::scal(&dimNY, beta, y, &i__1);
        if(sym) {
            K* res = new K[*n];
            for(int i = 0; i < *m; ++i) {
                std::fill_n(res, *n, K());
                for(int l = ia[i] - (N == 'F'); l < ia[i + 1] - (N == 'F'); ++l) {
                    int j = ja[l] - (N == 'F');
                    if(i != j)
                        for(int r = 0; r < *n; ++r) {
                            y[j + r * dimY] += *alpha * a[l] * x[i + r * dimX];
                            res[r] += a[l] * x[j + r * dimX];
                        }
                    else {
                        const K scal = *alpha * a[l];
                        Blas<K>::axpy(n, &scal, x + i, m, y + j, k);
                    }
                }
                Blas<K>::axpy(n, alpha, res, &i__1, y + i, k);
            }
            delete [] res;
        }
        else {
            for(int i = 0; i < *m; ++i)
                for(int l = ia[i] - (N == 'F'); l < ia[i + 1] - (N == 'F'); ++l) {
                    const K scal = *alpha * a[l];
                    Blas<K>::axpy(n, &scal, x + i, m, y + ja[l] - (N == 'F'), k);
                }
        }
    }
}

template<class K>
template<char N, char M>
inline void Wrapper<K>::csrcsc(const int* const n, const K* const a, const int* const ja, const int* const ia, K* const b, int* const jb, int* const ib) {
    unsigned int nnz = ia[*n] - (N == 'F');
    std::fill_n(ib, *n + 1, 0);
    for(unsigned int i = 0; i < nnz; ++i)
        ib[ja[i] + (N == 'C')]++;
    std::partial_sum(ib, ib + *n + 1, ib);
    for(unsigned int i = 0; i < *n; ++i)
        for(unsigned int j = ia[i] - (N == 'F'); j < ia[i + 1] - (N == 'F'); ++j) {
            unsigned int k = ib[ja[j] - (N == 'F')]++;
            jb[k] = i + (M == 'F');
            b[k] = a[j];
        }
    for(unsigned int i = *n; i > 0; --i)
        ib[i] = ib[i - 1] + (M == 'F');
    ib[0] = (M == 'F');
}
template<class K>
inline void Wrapper<K>::gthr(const int& n, const K* const y, K* const x, const int* const indx) {
    for(int i = 0; i < n; ++i)
        x[i] = y[indx[i]];
}
template<class K>
inline void Wrapper<K>::sctr(const int& n, const K* const x, const int* const indx, K* const y) {
    for(int i = 0; i < n; ++i)
        y[indx[i]] = x[i];
}
template<class K>
template<char O>
inline void Wrapper<K>::omatcopy(const int n, const int m, const K* const a, const int lda, K* const b, const int ldb) {
    static_assert(O == 'N' || O == 'R' || O == 'T' || O == 'C', "Unknown operation");
    if(O == 'T' || O == 'C')
        for(int i = 0; i < n; ++i)
            for(int j = 0; j < m; ++j) {
                if(O == 'T')
                    b[j * ldb + i] = a[i * lda + j];
                else
                    b[j * ldb + i] = conj(a[i * lda + j]);
            }
    if(O == 'R' && is_complex)
        for(int i = 0; i < n; ++i)
            std::transform(a + i * lda, a + i * lda + m, b + i * ldb, [](const K& z) { return conj(z); });
    else
        for(int i = 0; i < n; ++i)
            std::copy_n(a + i * lda, m, b + i * ldb);
}
template<class K>
template<char O>
inline void Wrapper<K>::imatcopy(const int n, const int m, K* const ab, const int lda, const int ldb) {
    static_assert(O == 'N' || O == 'R' || O == 'T' || O == 'C', "Unknown operation");
    if(O == 'T' || O == 'C') {
        if(n != 1 || m != 1) {
            if(lda == m && ldb == n) {
                if(n != m) {
                    const int size = n * m - 1;
                    std::bitset<1024> b;
                    b[0] = b[size] = 1;
                    int i = 1;
                    while(i < size) {
                        int it = i;
                        K t = ab[i];
                        do {
                            int next = (i * n) % size;
                            std::swap(ab[next], t);
                            b[i] = 1;
                            i = next;
                        } while(i != it);
                        if(O == 'C' && is_complex)
                            ab[i] = conj(ab[i]);

                        for(i = 1; i < size && b[i]; ++i);
                    }
                }
                else {
                    for(int i = 0; i < n - 1; ++i)
                        for(int j = i + 1; j < n; ++j) {
                            if(O == 'C' && is_complex) {
                                ab[i * n + j] = conj(ab[i * n + j]);
                                ab[j * n + i] = conj(ab[j * n + i]);
                                std::swap(ab[i * n + j], ab[j * n + i]);
                            }
                            else
                                std::swap(ab[i * n + j], ab[j * n + i]);
                        }
                }
            }
            else {
                K* tmp = new K[n * m];
                omatcopy<O>(n, m, ab, lda, tmp, n);
                Blas<K>::lacpy("A", &n, &m, tmp, &n, ab, &ldb);
                delete [] tmp;
            }
        }
    }
    else if(O == 'R' && is_complex) {
        if(lda == ldb) {
            for(int i = 0; i < n; ++i)
                std::for_each(ab + i * lda, ab + i * lda + m, [](K& z) { z = conj(z); });
        }
        else if (lda < ldb) {
            for(int i = n; i-- > 0; )
                for(int j = m; j-- > 0; )
                    ab[i * ldb + j] = conj(ab[i * lda + j]);
        }
        else {
            for(int i = 0; i < n; ++i)
                for(int j = 0; j < m; ++j)
                    ab[i * ldb + j] = conj(ab[i * lda + j]);
        }
    }
    else {
        if(lda < ldb)
            for(int i = n; i > 0; --i)
                std::copy_backward(ab + (i - 1) * lda, ab + (i - 1) * lda + m, ab + (i - 1) * ldb + m);
        else if(lda > ldb)
            for(int i = 1; i < n; ++i)
                std::copy_n(ab + i * ldb, m, ab + i * lda);
    }
}
#endif // HPDDM_MKL

template<class K>
inline void Wrapper<K>::diag(const int& n, const underlying_type<K>* const d, const K* const in, K* const out) {
    diag(n, i__1, d, in, out);
}
template<class K>
inline void Wrapper<K>::diag(const int& m, const int& n, const underlying_type<K>* const d, const K* const in, K* const out) {
    if(in)
        for(int i = 0; i < n; ++i)
            for(int j = 0; j < m; ++j)
                out[j + i * m] = d[j] * in[j + i * m];
    else
        for(int i = 0; i < n; ++i)
            for(int j = 0; j < m; ++j)
                out[j + i * m] *= d[j];
}
template<class K>
template<char N>
inline void Wrapper<K>::csrmm(bool sym, const int* const n, const int* const m, const K* const a, const int* const ia, const int* const ja, const K* const x, K* const y) {
    csrmm<N>("N", n, m, n, &d__1, sym, a, ia, ja, x, n, &d__0, y, n);
}

template<class Idx, class T>
inline void reorder(const Idx& i, const Idx& j, const T& v) {
    std::swap(v[i], v[j]);
}
template<class Idx, class First, class... Rest>
inline void reorder(const Idx& i, const Idx& j, const First& first, const Rest&... rest) {
    std::swap(first[i], first[j]);
    reorder(i, j, rest...);
}
/* Function: reorder
 *  Rearranges an arbitrary number of containers based on the permutation defined by the first argument. */
template<class T, class... Args>
inline void reorder(std::vector<T>& order, const Args&... args) {
    static_assert(sizeof...(args) > 0, "Nothing to reorder");
    for(T i = 0; i < order.size() - 1; ++i) {
        T j = order[i];
        if(j != i) {
            T k = i + 1;
            while(order[k] != i)
                ++k;
            std::swap(order[i], order[k]);
            reorder(i, j, args...);
        }
    }
}
} // HPDDM
#endif // _HPDDM_WRAPPER_
