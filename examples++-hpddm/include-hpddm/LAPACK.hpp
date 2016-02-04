/*
   This file is part of HPDDM.

   Author(s): Pierre Jolivet <pierre.jolivet@enseeiht.fr>
        Date: 2014-03-16

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

#ifndef _HPDDM_LAPACK_
#define _HPDDM_LAPACK_

#define HPDDM_GENERATE_EXTERN_LAPACK(C, T, U, SYM, ORT)                                                      \
void HPDDM_F77(C ## lapmt)(const int*, const int*, const int*, T*, const int*, int*);                        \
U    HPDDM_F77(C ## lange)(const char*, const int*, const int*, const T*, const int*, U*);                   \
U    HPDDM_F77(C ## lan ## SYM)(const char*, const char*, const int*, const T*, const int*, U*);             \
void HPDDM_F77(C ## SYM ## gst)(const int*, const char*, const int*, T*, const int*,                         \
                                const T*, const int*, int*);                                                 \
void HPDDM_F77(C ## SYM ## trd)(const char*, const int*, T*, const int*, U*, U*, T*, T*, const int*, int*);  \
void HPDDM_F77(C ## stein)(const int*, const U*, const U*, const int*, const U*, const int*,                 \
                           const int*, T*, const int*, U*, int*, int*, int*);                                \
void HPDDM_F77(C ## ORT ## mtr)(const char*, const char*, const char*, const int*, const int*,               \
                                const T*, const int*, const T*, T*, const int*, T*, const int*, int*);       \
void HPDDM_F77(C ## potrf)(const char*, const int*, T*, const int*, int*);                                   \
void HPDDM_F77(C ## potrs)(const char*, const int*, const int*, const T*, const int*, T*, const int*, int*); \
void HPDDM_F77(C ## pstrf)(const char*, const int*, T*, const int*, int*, int*, const U*, U*, int*);         \
void HPDDM_F77(C ## trtrs)(const char*, const char*, const char*, const int*, const int*, const T*,          \
                           const int*, T*, const int*, int*);                                                \
void HPDDM_F77(C ## geqrf)(const int*, const int*, T*, const int*, T*, T*, const int*, int*);                \
void HPDDM_F77(C ## geqrt)(const int*, const int*, const int*, T*, const int*, T*, const int*, T*, int*);    \
void HPDDM_F77(C ## gemqrt)(const char*, const char*, const int*, const int*, const int*, const int*,        \
                            const T*, const int*, const T*, const int*, T*, const int*, T*, int*);
#define HPDDM_GENERATE_EXTERN_LAPACK_COMPLEX(C, T, B, U)                                                     \
HPDDM_GENERATE_EXTERN_LAPACK(B, U, U, sy, or)                                                                \
HPDDM_GENERATE_EXTERN_LAPACK(C, T, U, he, un)                                                                \
void HPDDM_F77(B ## stebz)(const char*, const char*, const int*, const U*, const U*, const int*, const int*, \
                           const U*, const U*, const U*, int*, int*, U*, int*, int*, U*, int*, int*);        \
void HPDDM_F77(B ## pocon)(const char*, const int*, const U*, const int*, U*, U*, U*, int*, int*);           \
void HPDDM_F77(C ## pocon)(const char*, const int*, const T*, const int*, U*, U*, T*, U*, int*);             \
void HPDDM_F77(B ## geqp3)(const int*, const int*, U*, const int*, const int*, U*, U*, const int*, int*);    \
void HPDDM_F77(C ## geqp3)(const int*, const int*, T*, const int*, const int*, T*, T*, const int*, U*, int*);\
void HPDDM_F77(B ## ormqr)(const char*, const char*, const int*, const int*, const int*, const U*,           \
                           const int*, const U*, U*, const int*, U*, const int*, int*);                      \
void HPDDM_F77(C ## unmqr)(const char*, const char*, const int*, const int*, const int*, const T*,           \
                           const int*, const T*, T*, const int*, T*, const int*, int*);                      \
void HPDDM_F77(B ## hseqr)(const char*, const char*, const int*, const int*, const int*, U*, const int*, U*, \
                           U*, U*, const int*, U*, const int*, int*);                                        \
void HPDDM_F77(C ## hseqr)(const char*, const char*, const int*, const int*, const int*, T*, const int*, T*, \
                           T*, const int*, T*, const int*, int*);                                            \
void HPDDM_F77(B ## hsein)(const char*, const char*, const char*, int*, const int*, U*, const int*, U*,      \
                           const U*, U*, const int*, U*, const int*, const int*, int*, U*, int*, int*, int*);\
void HPDDM_F77(C ## hsein)(const char*, const char*, const char*, int*, const int*, T*, const int*, T*, T*,  \
                           const int*, T*, const int*, const int*, int*, T*, U*, int*, int*, int*);          \
void HPDDM_F77(B ## geev)(const char*, const char*, const int*, U*, const int*, U*, U*, U*, const int*, U*,  \
                          const int*, U*, const int*, int*);                                                 \
void HPDDM_F77(C ## geev)(const char*, const char*, const int*, T*, const int*, T*, T*, const int*, T*,      \
                          const int*, T*, const int*, U*, int*);                                             \
void HPDDM_F77(B ## ggev)(const char*, const char*, const int*, U*, const int*, U*, const int*, U*, U*, U*,  \
                          U*, const int*, U*, const int*, U*, const int*, int*);                             \
void HPDDM_F77(C ## ggev)(const char*, const char*, const int*, T*, const int*, T*, const int*, T*, T*,      \
                          T*, const int*, T*, const int*, T*, const int*, U*, int*);                         \
void HPDDM_F77(B ## gesdd)(const char*, const int*, const int*, U*, const int*, U*, U*, const int*, U*,      \
                           const int*, U*, const int*, int*, int*);                                          \
void HPDDM_F77(C ## gesdd)(const char*, const int*, const int*, T*, const int*, U*, T*, const int*, T*,      \
                           const int*, T*, const int*, U*, int*, int*);

#ifndef INTEL_MKL_VERSION
# ifdef __cplusplus
extern "C" {
HPDDM_GENERATE_EXTERN_LAPACK_COMPLEX(c, std::complex<float>, s, float)
HPDDM_GENERATE_EXTERN_LAPACK_COMPLEX(z, std::complex<double>, d, double)
}
# else
HPDDM_GENERATE_EXTERN_LAPACK_COMPLEX(c, void, s, float)
HPDDM_GENERATE_EXTERN_LAPACK_COMPLEX(z, void, d, double)
# endif // __cplusplus
#endif // INTEL_MKL_VERSION

#ifdef __cplusplus
namespace HPDDM {
/* Class: Lapack
 *
 *  A class inheriting from <Eigensolver> to use <Lapack> for dense eigenvalue problems.
 *
 * Template Parameter:
 *    K              - Scalar type. */
template<class K>
class Lapack : public Eigensolver<K> {
    private:
        /* Function: gst
         *  Reduces a symmetric or Hermitian definite generalized eigenvalue problem to a standard form. */
        static void gst(const int*, const char*, const int*, K*, const int*, K*, const int*, int*);
        /* Function: trd
         *  Reduces a symmetric or Hermitian matrix to a tridiagonal form. */
        static void trd(const char*, const int*, K*, const int*, underlying_type<K>*, underlying_type<K>*, K*, K*, const int*, int*);
        /* Function: stein
         *  Computes the eigenvectors corresponding to specified eigenvalues of a symmetric tridiagonal matrix. */
        static void stein(const int*, const underlying_type<K>*, const underlying_type<K>*, const int*, const underlying_type<K>*, const int*, const int*, K*, const int*, underlying_type<K>*, int*, int*, int*);
        /* Function: stebz
         *  Computes selected eigenvalues of a symmetric tridiagonal matrix by bisection. */
        static void stebz(const char*, const char*, const int*, const underlying_type<K>*, const underlying_type<K>*, const int*, const int*, const underlying_type<K>*, const underlying_type<K>*, const underlying_type<K>*, int*, int*, underlying_type<K>*, int*, int*, underlying_type<K>*, int*, int*);
        /* Function: mtr
         *  Multiplies a matrix by an orthogonal or unitary matrix obtained with <Lapack::trd>. */
        static void mtr(const char*, const char*, const char*, const int*, const int*, const K*, const int*, const K*, K*, const int*, K*, const int*, int*);
        /* Function: gesdd
         *  Computes the singular value decomposition of a rectangular matrix, and optionally the left and/or right singular vectors, using a divide and conquer algorithm. */
        static void gesdd(const char*, const int*, const int*, K*, const int*, underlying_type<K>*, K*, const int*, K*, const int*, K*, const int*, underlying_type<K>*, int*, int*);
    public:
        Lapack(int n)                                                               : Eigensolver<K>(n) { }
        Lapack(int n, int nu)                                                       : Eigensolver<K>(n, nu) { }
        Lapack(underlying_type<K> threshold, int n, int nu)                         : Eigensolver<K>(threshold, n, nu) { }
        Lapack(underlying_type<K> tol, underlying_type<K> threshold, int n, int nu) : Eigensolver<K>(tol, threshold, n, nu) { }
        /* Function: lapmt
         *  Performs a forward or backward permutation of the columns of a matrix. */
        static void lapmt(const int*, const int*, const int*, K*, const int*, int*);
        /* Function: lange
         *  Computes the norm of a general rectangular matrix. */
        static underlying_type<K> lange(const char*, const int*, const int*, const K*, const int*, underlying_type<K>*);
        /* Function: lan
         *  Computes the norm of a symmetric or Hermitian matrix. */
        static underlying_type<K> lan(const char*, const char*, const int*, const K*, const int*, underlying_type<K>*);
        /* Function: potrf
         *  Computes the Cholesky factorization of a symmetric or Hermitian positive definite matrix. */
        static void potrf(const char*, const int*, K*, const int*, int*);
        /* Function: pocon
         *  Estimates the reciprocal of the condition number of a symmetric or Hermitian positive definite matrix. */
        static void pocon(const char*, const int*, const K*, const int*, underlying_type<K>*, underlying_type<K>*, K*, typename std::conditional<Wrapper<K>::is_complex, underlying_type<K>*, int*>::type, int*);
        /* Function: potrs
         *  Solves a system of linear equations with a Cholesky-factored matrix. */
        static void potrs(const char*, const int*, const int*, const K*, const int*, K*, const int*, int*);
        /* Function: pstrf
         *  Computes the Cholesky factorization of a symmetric or Hermitian positive semidefinite matrix with pivoting. */
        static void pstrf(const char*, const int*, K*, const int*, int*, int*, const underlying_type<K>*, underlying_type<K>*, int*);
        /* Function: trtrs
         *  Solves a system of linear equations with a triangular matrix. */
        static void trtrs(const char*, const char*, const char*, const int*, const int*, const K*, const int*, K*, const int*, int*);
        /* Function: geqp3
         *  Computes a QR decomposition of a rectangular matrix with column pivoting. */
        static void geqp3(const int*, const int*, K*, const int*, int*, K*, K*, const int*, underlying_type<K>*, int*);
        /* Function: geqrf
         *  Computes a QR decomposition of a rectangular matrix. */
        static void geqrf(const int*, const int*, K*, const int*, K*, K*, const int*, int*);
        /* Function: geqrt
         *  Computes a blocked QR decomposition of a rectangular matrix using the compact WY representation of Q. */
        static void geqrt(const int*, const int*, const int*, K*, const int*, K*, const int*, K*, int*);
        /* Function: gemqrt
         *  Multiplies a matrix by an orthogonal or unitary matrix obtained with <Lapack::geqrt>. */
        static void gemqrt(const char*, const char*, const int*, const int*, const int*, const int*, const K*, const int*, const K*, const int*, K*, const int*, K*, int*);
        /* Function: mqr
         *  Multiplies a matrix by an orthogonal or unitary matrix obtained with <Lapack::geq>. */
        static void mqr(const char*, const char*, const int*, const int*, const int*, const K*, const int*, const K*, K*, const int*, K*, const int*, int*);
        /* Function: hseqr
         *  Computes all eigenvalues and (optionally) the Schur factorization of an upper Hessenberg matrix. */
        static void hseqr(const char*, const char*, const int*, const int*, const int*, K*, const int*, K*, K*, K*, const int*, K*, const int*, int*);
        /* Function: hsein
         *  Computes selected eigenvectors of an upper Hessenberg matrix that correspond to specified eigenvalues. */
        static void hsein(const char*, const char*, const char*, int*, const int*, K*, const int*, K*, const K*, K*, const int*, K*, const int*, const int*, int*, K*, underlying_type<K>*, int*, int*, int*);
        /* Function: geev
         *  Computes the eigenvalues and the eigenvectors of a nonsymmetric eigenvalue problem. */
        static void geev(const char*, const char*, const int*, K*, const int*, K*, K*, K*, const int*, K*, const int*, K*, const int*, underlying_type<K>*, int*);
        /* Function: ggev
         *  Computes the eigenvalues and the eigenvectors of a nonsymmetric generalized eigenvalue problem. */
        static void ggev(const char*, const char*, const int*, K*, const int*, K*, const int*, K*, K*, K*, K*, const int*, K*, const int*, K*, const int*, underlying_type<K>*, int*);

        /* Function: workspace
         *  Returns the optimal size of the workspace array. */
        int workspace() const {
            int info;
            int lwork = -1;
            K wkopt;
            trd("L", &(Eigensolver<K>::_n), nullptr, &(Eigensolver<K>::_n), nullptr, nullptr, nullptr, &wkopt, &lwork, &info);
            return static_cast<int>(std::real(wkopt));
        }
        /* Function: reduce
         *
         *  Reduces a symmetric or Hermitian definite generalized eigenvalue problem to a standard problem after factorizing the right-hand side matrix.
         *
         * Parameters:
         *    A              - Left-hand side matrix.
         *    B              - Right-hand side matrix. */
        void reduce(K* const& A, K* const& B) const {
            int info;
            potrf("L", &(Eigensolver<K>::_n), B, &(Eigensolver<K>::_n), &info);
            gst(&i__1, "L", &(Eigensolver<K>::_n), A, &(Eigensolver<K>::_n), B, &(Eigensolver<K>::_n), &info);
        }
        /* Function: expand
         *
         *  Computes the eigenvectors of a generalized eigenvalue problem after completion of <Lapack::solve>.
         *
         * Parameters:
         *    B              - Right-hand side matrix.
         *    ev             - Array of eigenvectors. */
        void expand(K* const& B, K* const* const ev) const {
            int info;
            trtrs("L", "T", "N", &(Eigensolver<K>::_n), &(Eigensolver<K>::_nu), B, &(Eigensolver<K>::_n), *ev, &(Eigensolver<K>::_n), &info);
        }
        /* Function: solve
         *
         *  Computes eigenvectors of the standard eigenvalue problem Ax = l x.
         *
         * Parameters:
         *    A              - Left-hand side matrix.
         *    ev             - Array of eigenvectors.
         *    work           - Workspace array.
         *    lwork          - Size of the input workspace array.
         *    communicator   - MPI communicator for selecting the threshold criterion. */
        void solve(K* const& A, K**& ev, K* const& work, int& lwork, const MPI_Comm& communicator) {
            int info;
            K* tau = work + lwork;
            underlying_type<K>* d = reinterpret_cast<underlying_type<K>*>(tau + Eigensolver<K>::_n);
            underlying_type<K>* e = d + Eigensolver<K>::_n;
            trd("L", &(Eigensolver<K>::_n), A, &(Eigensolver<K>::_n), d, e, tau, work, &lwork, &info);
            underlying_type<K> vl = -1.0 / HPDDM_EPS;
            underlying_type<K> vu = Eigensolver<K>::_threshold;
            int iu = Eigensolver<K>::_nu;
            int nsplit;
            underlying_type<K>* evr = e + Eigensolver<K>::_n - 1;
            int* iblock = new int[5 * Eigensolver<K>::_n];
            int* isplit = iblock + Eigensolver<K>::_n;
            int* iwork = isplit + Eigensolver<K>::_n;
            char range = Eigensolver<K>::_threshold > 0.0 ? 'V' : 'I';
            stebz(&range, "B", &(Eigensolver<K>::_n), &vl, &vu, &i__1, &iu, &(Eigensolver<K>::_tol), d, e, &(Eigensolver<K>::_nu), &nsplit, evr, iblock, isplit, reinterpret_cast<underlying_type<K>*>(work), iwork, &info);
            if(Eigensolver<K>::_nu) {
                ev = new K*[Eigensolver<K>::_nu];
                *ev = new K[Eigensolver<K>::_n * Eigensolver<K>::_nu];
                for(unsigned short i = 1; i < Eigensolver<K>::_nu; ++i)
                    ev[i] = *ev + i * Eigensolver<K>::_n;
                int* ifailv = new int[Eigensolver<K>::_nu];
                stein(&(Eigensolver<K>::_n), d, e, &(Eigensolver<K>::_nu), evr, iblock, isplit, *ev, &(Eigensolver<K>::_n), reinterpret_cast<underlying_type<K>*>(work), iwork, ifailv, &info);
                delete [] ifailv;
                mtr("L", "L", "N", &(Eigensolver<K>::_n), &(Eigensolver<K>::_nu), A, &(Eigensolver<K>::_n), tau, *ev, &(Eigensolver<K>::_n), work, &lwork, &info);
                if(!Wrapper<K>::is_complex)
                    lwork += 3 * Eigensolver<K>::_n - 1;
                else
                    lwork += 4 * Eigensolver<K>::_n - 1;
            }
            delete [] iblock;
        }
        int workspace(const char* jobz, const int* const m) const {
            int info;
            int lwork = -1;
            K wkopt;
            gesdd(jobz, &(Eigensolver<K>::_n), m, nullptr, &(Eigensolver<K>::_n), nullptr, nullptr, &(Eigensolver<K>::_n), nullptr, m, &wkopt, &lwork, nullptr, nullptr, &info);
            return static_cast<int>(std::real(wkopt));
        }
        void svd(const char* jobz, const int* m, K* a, underlying_type<K>* s, K* u, K* vt, K* work, const int* lwork, int* iwork, underlying_type<K>* rwork = nullptr) const {
            int info;
            gesdd(jobz, &(Eigensolver<K>::_n), m, a, &(Eigensolver<K>::_n), s, u, &(Eigensolver<K>::_n), vt, m, work, lwork, rwork, iwork, &info);
        }
        void purify(K* ev, const underlying_type<K>* const d = nullptr) {
            int lwork = workspace("N", &(Eigensolver<K>::_nu));
            K* a, *work;
            underlying_type<K>* rwork, *s;
            if(!Wrapper<K>::is_complex) {
                a = new K[&(Eigensolver<K>::_n) * Eigensolver<K>::_nu + lwork + Eigensolver<K>::_nu];
                work = a + &(Eigensolver<K>::_n) * Eigensolver<K>::_nu;
                s = reinterpret_cast<underlying_type<K>*>(work) + lwork;
                rwork = nullptr;
            }
            else {
                a = new K[&(Eigensolver<K>::_n) * Eigensolver<K>::_nu + lwork];
                work = a + &(Eigensolver<K>::_n) * Eigensolver<K>::_nu;
                s = new underlying_type<K>[Eigensolver<K>::_nu + std::max(1, Eigensolver<K>::_nu * std::max(5 * Eigensolver<K>::_nu + 7, 2 * &(Eigensolver<K>::_n) + 2 * Eigensolver<K>::_nu + 1))];
                rwork = s + Eigensolver<K>::_nu;
            }
            if(d)
                Wrapper<K>::diag(Eigensolver<K>::_n, d, ev, a, Eigensolver<K>::_nu);
            else
                std::copy_n(ev, &(Eigensolver<K>::_n) * Eigensolver<K>::_nu, a);
            int* iwork = new int[8 * Eigensolver<K>::_n];
            int info;
            gesdd("N", &(Eigensolver<K>::_n), &(Eigensolver<K>::_nu), a, &(Eigensolver<K>::_n), s, nullptr, &(Eigensolver<K>::_n), nullptr, &(Eigensolver<K>::_nu), work, &lwork, rwork, iwork, &info);
            delete [] iwork;
            if(Wrapper<K>::is_complex)
                delete [] s;
            delete [] a;
        }
};

/* Class: QR
 *
 *  A class to use LAPACK for computing QR decompositions.
 *
 * Template Parameter:
 *    K              - Scalar type. */
template<class K>
class QR {
    private:
        int                               _n;
        int                           _lwork;
        K* const                          _a;
        K* const                        _tau;
        K* const                       _work;
# if HPDDM_QR == 1
        std::vector<int>               _jpvt;
        int                            _rank;
# endif
        /* Function: workspace
         *  Returns the optimal size of the workspace array. */
        int workspace() const {
            int info;
            int lwork[2] { -1, -1 };
            K wkopt;
# if HPDDM_QR == 1
            Lapack<K>::geqp3(&_n, &_n, nullptr, &_n, nullptr, nullptr, &wkopt, lwork, nullptr, &info);
# else
            Lapack<K>::geqrf(&_n, &_n, nullptr, &_n, nullptr, &wkopt, lwork, &info);
# endif
            lwork[0] = static_cast<int>(std::real(wkopt));
            Lapack<K>::mqr("L", "T", &_n, &i__1, &_n, nullptr, &_n, nullptr, nullptr, &_n, &wkopt, lwork + 1, &info);
            lwork[1] = static_cast<int>(std::real(wkopt));
            return *std::max_element(lwork, lwork + 1);
        }
    public:
        QR(int n, const K* const cpy = nullptr) : _n(n), _lwork(workspace()), _a(new K[_n * (_n + 1) + _lwork]), _tau(_a + _n * _n), _work(_tau + _n) {
# if HPDDM_QR == 1
            _jpvt.resize(_n);
            _rank = n;
# endif
            if(cpy)
                for(unsigned int i = 0; i < _n; ++i) {
                    _a[i * (_n + 1)] = cpy[i * (_n + 1)];
                    for(unsigned int j = i + 1; j < _n; ++j)
                        _a[j * _n + i] = _a[i * _n + j] = cpy[i * _n + j];
                }
        }
        ~QR() {
            delete [] _a;
        }
        /* Function: getPointer
         *  Returns the pointer <QR::a>. */
        K* getPointer() const { return _a; }
        void decompose() {
            int info;
# if HPDDM_QR == 1
            underlying_type<K>* rwork = Wrapper<K>::is_complex ? new underlying_type<K>[2 * _n] : nullptr;
            Lapack<K>::geqp3(&_n, &_n, _a, &_n, _jpvt.data(), _tau, _work, &_lwork, rwork, &info);
            delete [] rwork;
            while(std::abs(_a[(_rank - 1) * (_n + 1)]) < HPDDM_EPS * std::abs(_a[0]) && _rank-- > 0);
# else
            Lapack<K>::geqrf(&_n, &_n, _a, &_n, _tau, _work, &_lwork, &info);
            std::vector<int> jpvt;
            jpvt.reserve(6);
            underlying_type<K> max = std::abs(_a[0]);
            for(unsigned int i = 1; i < _n; ++i) {
                if(std::abs(_a[(_n + 1) * i]) < max * 1.0e-6)
                    jpvt.emplace_back(i);
                else if(std::abs(_a[(_n + 1) * i]) > max / 1.0e-6) {
                    jpvt.clear();
                    max = std::abs(_a[(_n + 1) * i]);
                    i = 0;
                }
                else
                    max = std::max(std::abs(_a[(i + 1) * _n]), max);
            }
            std::for_each(jpvt.cbegin(), jpvt.cend(), [&](const int i) { std::fill_n(_a + _n * i, i, K()); _a[(_n + 1) * i] = Wrapper<K>::d__1; });
# endif
        }
        /* Function: solve
         *  Computes the solution of a least squares problem. */
        void solve(K* const x) const {
            int info;
            Lapack<K>::mqr("L", "T", &_n, &i__1, &_n, _a, &_n, _tau, x, &_n, _work, &_lwork, &info);
# if HPDDM_QR == 1
            Lapack<K>::trtrs("U", "N", "N", &_rank, &i__1, _a, &_n, x, &_n, &info);
            Lapack<K>::lapmt(&i__0, &i__1, &_n, x, &i__1, const_cast<int*>(_jpvt.data()));
# else
            Lapack<K>::trtrs("U", "N", "N", &_n, &i__1, _a, &_n, x, &_n, &info);
# endif
        }
};

# define HPDDM_GENERATE_LAPACK(C, T, B, U, SYM, ORT)                                                         \
template<>                                                                                                   \
inline void Lapack<T>::lapmt(const int* forwrd, const int* m, const int* n, T* x, const int* ldx, int* k) {  \
    HPDDM_F77(C ## lapmt)(forwrd, m, n, x, ldx, k);                                                          \
}                                                                                                            \
template<>                                                                                                   \
inline U Lapack<T>::lange(const char* norm, const int* m, const int* n, const T* a, const int* lda,          \
                          U* work) {                                                                         \
    return HPDDM_F77(C ## lange)(norm, m, n, a, lda, work);                                                  \
}                                                                                                            \
template<>                                                                                                   \
inline U Lapack<T>::lan(const char* norm, const char* uplo, const int* m, const T* a, const int* lda,        \
                        U* work) {                                                                           \
    return HPDDM_F77(C ## lan ## SYM)(norm, uplo, m, a, lda, work);                                          \
}                                                                                                            \
template<>                                                                                                   \
inline void Lapack<T>::gst(const int* itype, const char* uplo, const int* n,                                 \
                           T* a, const int* lda, T* b, const int* ldb, int* info) {                          \
    HPDDM_F77(C ## SYM ## gst)(itype, uplo, n, a, lda, b, ldb, info);                                        \
}                                                                                                            \
template<>                                                                                                   \
inline void Lapack<T>::trd(const char* uplo, const int* n, T* a, const int* lda,                             \
                           U* d, U* e, T* tau, T* work, const int* lwork, int* info) {                       \
    HPDDM_F77(C ## SYM ## trd)(uplo, n, a, lda, d, e, tau, work, lwork, info);                               \
}                                                                                                            \
template<>                                                                                                   \
inline void Lapack<T>::stein(const int* n, const U* d, const U* e, const int* m, const U* w,                 \
                             const int* iblock, const int* isplit, T* z, const int* ldz,                     \
                             U* work, int* iwork, int* ifailv, int* info) {                                  \
    HPDDM_F77(C ## stein)(n, d, e, m, w, iblock, isplit, z, ldz, work, iwork, ifailv, info);                 \
}                                                                                                            \
template<>                                                                                                   \
inline void Lapack<T>::stebz(const char* range, const char* order, const int* n, const U* vl, const U* vu,   \
                             const int* il, const int* iu, const U* abstol, const U* d, const U* e, int* m,  \
                             int* nsplit, U* w, int* iblock, int* isplit, U* work, int* iwork, int* info) {  \
    HPDDM_F77(B ## stebz)(range, order, n, vl, vu, il, iu, abstol, d, e, m, nsplit, w, iblock, isplit,       \
                          work, iwork, info);                                                                \
}                                                                                                            \
template<>                                                                                                   \
inline void Lapack<T>::mtr(const char* side, const char* uplo, const char* trans, const int* m,              \
                           const int* n, const T* a, const int* lda, const T* tau, T* c, const int* ldc,     \
                           T* work, const int* lwork, int* info) {                                           \
    HPDDM_F77(C ## ORT ## mtr)(side, uplo, trans, m, n, a, lda, tau, c, ldc, work, lwork, info);             \
}                                                                                                            \
template<>                                                                                                   \
inline void Lapack<T>::potrf(const char* uplo, const int* n, T* a, const int* lda, int* info) {              \
    HPDDM_F77(C ## potrf)(uplo, n, a, lda, info);                                                            \
}                                                                                                            \
template<>                                                                                                   \
inline void Lapack<T>::potrs(const char* uplo, const int* n, const int* nrhs, const T* a, const int* lda,    \
                             T* b, const int* ldb, int* info) {                                              \
    HPDDM_F77(C ## potrs)(uplo, n, nrhs, a, lda, b, ldb, info);                                              \
}                                                                                                            \
template<>                                                                                                   \
inline void Lapack<T>::pstrf(const char* uplo, const int* n, T* a, const int* lda, int* piv, int* rank,      \
                             const U* tol, U* work, int* info) {                                             \
    HPDDM_F77(C ## pstrf)(uplo, n, a, lda, piv, rank, tol, work, info);                                      \
}                                                                                                            \
template<>                                                                                                   \
inline void Lapack<T>::trtrs(const char* uplo, const char* trans, const char* diag, const int* n,            \
                             const int* nrhs, const T* a, const int* lda, T* b, const int* ldb, int* info) { \
    HPDDM_F77(C ## trtrs)(uplo, trans, diag, n, nrhs, a, lda, b, ldb, info);                                 \
}                                                                                                            \
template<>                                                                                                   \
inline void Lapack<T>::geqrf(const int* m, const int* n, T* a, const int* lda, T* tau, T* work,              \
                             const int* lwork, int* info) {                                                  \
    HPDDM_F77(C ## geqrf)(m, n, a, lda, tau, work, lwork, info);                                             \
}                                                                                                            \
template<>                                                                                                   \
inline void Lapack<T>::geqrt(const int* m, const int* n, const int* nb, T* a, const int* lda, T* t,          \
                             const int* ldt, T* work, int* info) {                                           \
    HPDDM_F77(C ## geqrt)(m, n, nb, a, lda, t, ldt, work, info);                                             \
}                                                                                                            \
template<>                                                                                                   \
inline void Lapack<T>::gemqrt(const char* side, const char* trans, const int* m, const int* n, const int* k, \
                              const int* nb, const T* v, const int* ldv, const T* t, const int* ldt, T* c,   \
                              const int* ldc, T* work, int* info) {                                          \
    HPDDM_F77(C ## gemqrt)(side, trans, m, n, k, nb, v, ldv, t, ldt, c, ldc, work, info);                    \
}
# define HPDDM_GENERATE_LAPACK_COMPLEX(C, T, B, U)                                                           \
HPDDM_GENERATE_LAPACK(B, U, B, U, sy, or)                                                                    \
HPDDM_GENERATE_LAPACK(C, T, B, U, he, un)                                                                    \
template<>                                                                                                   \
inline void Lapack<U>::pocon(const char* uplo, const int* n, const U* a, const int* lda, U* anorm, U* rcond, \
                             U* work, int* iwork, int* info) {                                               \
    HPDDM_F77(B ## pocon)(uplo, n, a, lda, anorm, rcond, work, iwork, info);                                 \
}                                                                                                            \
template<>                                                                                                   \
inline void Lapack<T>::pocon(const char* uplo, const int* n, const T* a, const int* lda, U* anorm, U* rcond, \
                             T* work, U* rwork, int* info) {                                                 \
    HPDDM_F77(C ## pocon)(uplo, n, a, lda, anorm, rcond, work, rwork, info);                                 \
}                                                                                                            \
template<>                                                                                                   \
inline void Lapack<U>::geqp3(const int* m, const int* n, U* a, const int* lda, int* jpvt, U* tau, U* work,   \
                             const int* lwork, U*, int* info) {                                              \
    HPDDM_F77(B ## geqp3)(m, n, a, lda, jpvt, tau, work, lwork, info);                                       \
}                                                                                                            \
template<>                                                                                                   \
inline void Lapack<T>::geqp3(const int* m, const int* n, T* a, const int* lda, int* jpvt, T* tau, T* work,   \
                             const int* lwork, U* rwork, int* info) {                                        \
    HPDDM_F77(C ## geqp3)(m, n, a, lda, jpvt, tau, work, lwork, rwork, info);                                \
}                                                                                                            \
template<>                                                                                                   \
inline void Lapack<U>::mqr(const char* side, const char* trans, const int* m, const int* n, const int* k,    \
                           const U* a, const int* lda, const U* tau, U* c, const int* ldc, U* work,          \
                           const int* lwork, int* info) {                                                    \
    HPDDM_F77(B ## ormqr)(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);                     \
}                                                                                                            \
template<>                                                                                                   \
inline void Lapack<T>::mqr(const char* side, const char* trans, const int* m, const int* n, const int* k,    \
                           const T* a, const int* lda, const T* tau, T* c, const int* ldc, T* work,          \
                           const int* lwork, int* info) {                                                    \
    HPDDM_F77(C ## unmqr)(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);                     \
}                                                                                                            \
template<>                                                                                                   \
inline void Lapack<U>::hseqr(const char* job, const char* compz, const int* n, const int* ilo,               \
                             const int* ihi, U* h, const int* ldh, U* wr, U* wi, U* z, const int* ldz,       \
                             U* work, const int* lwork, int* info) {                                         \
    HPDDM_F77(B ## hseqr)(job, compz, n, ilo, ihi, h, ldh, wr, wi, z, ldz, work, lwork, info);               \
}                                                                                                            \
template<>                                                                                                   \
inline void Lapack<T>::hseqr(const char* job, const char* compz, const int* n, const int* ilo,               \
                             const int* ihi, T* h, const int* ldh, T* w, T*, T* z, const int* ldz,           \
                             T* work, const int* lwork, int* info) {                                         \
    HPDDM_F77(C ## hseqr)(job, compz, n, ilo, ihi, h, ldh, w, z, ldz, work, lwork, info);                    \
}                                                                                                            \
template<>                                                                                                   \
inline void Lapack<U>::hsein(const char* side, const char* eigsrc, const char* initv, int* select,           \
                             const int* n, U* h, const int* ldh, U* wr, const U* wi, U* vl, const int* ldvl, \
                             U* vr, const int* ldvr, const int* mm, int* m, U* work, U*, int* ifaill,        \
                             int* ifailr, int* info) {                                                       \
    HPDDM_F77(B ## hsein)(side, eigsrc, initv, select, n, h, ldh, wr, wi, vl, ldvl, vr, ldvr, mm, m, work,   \
                          ifaill, ifailr, info);                                                             \
}                                                                                                            \
template<>                                                                                                   \
inline void Lapack<T>::hsein(const char* side, const char* eigsrc, const char* initv, int* select,           \
                             const int* n, T* h, const int* ldh, T* w, const T*, T* vl, const int* ldvl,     \
                             T* vr, const int* ldvr, const int* mm, int* m, T* work, U* rwork, int* ifaill,  \
                             int* ifailr, int* info) {                                                       \
    HPDDM_F77(C ## hsein)(side, eigsrc, initv, select, n, h, ldh, w, vl, ldvl, vr, ldvr, mm, m, work, rwork, \
                          ifaill, ifailr, info);                                                             \
}                                                                                                            \
template<>                                                                                                   \
inline void Lapack<U>::geev(const char* jobvl, const char* jobvr, const int* n, U* a, const int* lda,        \
                            U* wr, U* wi, U* vl, const int* ldvl, U* vr, const int* ldvr, U* work,           \
                            const int* lwork, U*, int* info) {                                               \
    HPDDM_F77(B ## geev)(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info);            \
}                                                                                                            \
template<>                                                                                                   \
inline void Lapack<T>::geev(const char* jobvl, const char* jobvr, const int* n, T* a, const int* lda, T* w,  \
                            T*, T* vl, const int* ldvl, T* vr, const int* ldvr, T* work, const int* lwork,   \
                            U* rwork, int* info) {                                                           \
    HPDDM_F77(C ## geev)(jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork, info);          \
}                                                                                                            \
template<>                                                                                                   \
inline void Lapack<U>::ggev(const char* jobvl, const char* jobvr, const int* n, U* a, const int* lda, U* b,  \
                            const int* ldb, U* alphar, U* alphai, U* beta, U* vl, const int* ldvl, U* vr,    \
                            const int* ldvr, U* work, const int* lwork, U*, int* info) {                     \
    HPDDM_F77(B ## ggev)(jobvl, jobvr, n, a, lda, b, ldb, alphar, alphai, beta, vl, ldvl, vr, ldvr, work,    \
                         lwork, info);                                                                       \
}                                                                                                            \
template<>                                                                                                   \
inline void Lapack<T>::ggev(const char* jobvl, const char* jobvr, const int* n, T* a, const int* lda, T* b,  \
                            const int* ldb, T* alpha, T*, T* beta, T* vl, const int* ldvl, T* vr,            \
                            const int* ldvr, T* work, const int* lwork, U* rwork, int* info) {               \
    HPDDM_F77(C ## ggev)(jobvl, jobvr, n, a, lda, b, ldb, alpha, beta, vl, ldvl, vr, ldvr, work, lwork,      \
                         rwork, info);                                                                       \
}                                                                                                            \
template<>                                                                                                   \
inline void Lapack<U>::gesdd(const char* jobz, const int* m, const int* n, U* a, const int* lda, U* s,       \
                             U* u, const int* ldu, U* vt, const int* ldvt, U* work, const int* lwork,        \
                             U*, int* iwork, int* info) {                                                    \
    HPDDM_F77(B ## gesdd)(jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, iwork, info);                \
}                                                                                                            \
template<>                                                                                                   \
inline void Lapack<T>::gesdd(const char* jobz, const int* m, const int* n, T* a, const int* lda, U* s,       \
                             T* u, const int* ldu, T* vt, const int* ldvt, T* work, const int* lwork,        \
                             U* rwork, int* iwork, int* info) {                                              \
    HPDDM_F77(C ## gesdd)(jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, iwork, info);         \
}
HPDDM_GENERATE_LAPACK_COMPLEX(c, std::complex<float>, s, float)
HPDDM_GENERATE_LAPACK_COMPLEX(z, std::complex<double>, d, double)
} // HPDDM
#endif // __cplusplus
#endif // _HPDDM_LAPACK_
