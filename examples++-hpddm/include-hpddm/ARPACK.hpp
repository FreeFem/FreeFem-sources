/*
   This file is part of HPDDM.

   Author(s): Pierre Jolivet <pierre.jolivet@enseeiht.fr>
        Date: 2012-12-15

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

#ifndef _HPDDM_ARPACK_
#define _HPDDM_ARPACK_

#define EIGENSOLVER HPDDM::Arpack

#define HPDDM_GENERATE_ARPACK_EXTERN(C, T, B, U)                                                              \
void HPDDM_F77(B ## saupd)(int*, const char*, const int*, const char*, const int*, const U*, U*, int*, U*,    \
                           const int*, int*, int*, U*, U*, int*, int*, int, int);                             \
void HPDDM_F77(B ## seupd)(const int*, const char*, int*, U*, U*, const int*, const U*, const char*,          \
                           const int*, const char*, const int*, const U*, U*, int*, U*, const int*, int*,     \
                           int*, U*, U*, int*, int*, int, int, int);                                          \
void HPDDM_F77(C ## naupd)(int*, const char*, const int*, const char*, const int*, const U*, T*, int*,        \
                           T*, const int*, int*, int*, T*, T*, int*, U*, int*, int, int);                     \
void HPDDM_F77(C ## neupd)(const int*, const char*, int*, T*, T*, const int*, const T*, T*, const char*,      \
                           const int*, const char*, const int*, const U*, T*, int*, T*, const int*, int*,     \
                           int*, T*, T*, int*, U*, int*, int, int, int);

extern "C" {
HPDDM_GENERATE_ARPACK_EXTERN(c, std::complex<float>, s, float)
HPDDM_GENERATE_ARPACK_EXTERN(z, std::complex<double>, d, double)
}

namespace HPDDM {
/* Class: Arpack
 *
 *  A class inheriting from <Eigensolver> to use <Arpack> for sparse eigenvalue problems.
 *
 * Template Parameter:
 *    K              - Scalar type. */
template<class K>
class Arpack : public Eigensolver<K> {
    private:
        /* Variable: it
         *  Maximum number of iterations of the IRAM. */
        unsigned short                        _it;
        /* Variable: which
         *  Eigenpairs to retrieve. */
        static constexpr const char* const _which = Wrapper<K>::is_complex ? "LM" : "LM";
        /* Function: aupd
         *  Iterates the implicitly restarted Arnoldi method. */
        static void aupd(int*, const char*, const int*, const char*, const int*, const underlying_type<K>*, K*, int*, K*, int*, int*, K*, K*, int*, underlying_type<K>*, int*);
        /* Function: eupd
         *  Post-processes the eigenpairs computed with <Arpack::aupd>. */
        static void eupd(const int*, const char*, int*, K*, K*, const int*, const K*, K*, const char*, const char*, const int*, const underlying_type<K>*, K*, int*, K*, int*, int*, K*, K*, int*, underlying_type<K>*, int*);
    public:
        Arpack(int n, int nu)                                                                          : Eigensolver<K>(n, nu), _it(100) { }
        Arpack(underlying_type<K> threshold, int n, int nu)                                            : Eigensolver<K>(threshold, n, nu), _it(100) { }
        Arpack(underlying_type<K> tol, underlying_type<K> threshold, int n, int nu, unsigned short it) : Eigensolver<K>(tol, threshold, n, nu), _it(it) { }
        /* Function: solve
         *
         *  Computes eigenvectors of the generalized eigenvalue problem Ax = l Bx.
         *
         * Parameters:
         *    A              - Left-hand side matrix.
         *    B              - Right-hand side matrix.
         *    ev             - Array of eigenvectors.
         *    communicator   - MPI communicator for selecting the threshold criterion. */
        template<template<class> class Solver>
        void solve(MatrixCSR<K>* const& A, MatrixCSR<K>* const& B, K**& ev, const MPI_Comm& communicator, Solver<K>* const& s = nullptr) {
            int iparam[11] { 1, 0, _it, 1, 0, 0, 3, 0, 0, 0, 0 };
            int ipntr[!Wrapper<K>::is_complex ? 11 : 14] { };
            if(4 * Eigensolver<K>::_nu > Eigensolver<K>::_n)
                Eigensolver<K>::_nu = std::max(1, Eigensolver<K>::_n / 4);
            int ncv = 2 * Eigensolver<K>::_nu + 1;
            int lworkl = ncv * (ncv + 8);
            K* workd;
            K* workev;
            if(!Wrapper<K>::is_complex) {
                workd = new K[lworkl + (ncv + 4) * Eigensolver<K>::_n];
                workev = nullptr;
            }
            else {
                lworkl += ncv * (2 * ncv - 3);
                workd = new K[lworkl + (ncv + 4) * Eigensolver<K>::_n + 2 * ncv];
                workev = workd + lworkl + (ncv + 4) * Eigensolver<K>::_n;
            }
            K* workl = workd + 3 * Eigensolver<K>::_n;
            K* vp = workl + lworkl;
            underlying_type<K>* rwork = nullptr;
            if(Wrapper<K>::is_complex)
                rwork = new underlying_type<K>[Eigensolver<K>::_n];
            K* vresid = vp + ncv * Eigensolver<K>::_n;
            Solver<K>* const prec = s ? s : new Solver<K>;
#ifdef MUMPSSUB
            prec->numfact(A, false);
#else
            prec->numfact(A, true);
#endif
            int info = 0;
            do {
                const int* const n = &(Eigensolver<K>::_n), *const nu = &(Eigensolver<K>::_nu);
                const underlying_type<K>* const tol = &(Eigensolver<K>::_tol);
                auto loop = [&]() {
                    int ido = 0;
                    while(ido != 99) {
                        aupd(&ido, "G", n, _which, nu, tol, vresid, &ncv,
                             vp, iparam, ipntr, workd, workl, &lworkl, rwork, &info);
                        if(ido == -1) {
                            Wrapper<K>::csrmv(B->_sym, n, B->_a, B->_ia, B->_ja, workd + ipntr[0] - 1, workd + ipntr[1] - 1);
                            prec->solve(workd + ipntr[1] - 1);
                        }
                        else if(ido == 1)
                            prec->solve(workd + ipntr[2] - 1, workd + ipntr[1] - 1);
                        else
                            Wrapper<K>::csrmv(B->_sym, n, B->_a, B->_ia, B->_ja, workd + ipntr[0] - 1, workd + ipntr[1] - 1);
                    }
                };
                loop();
                if(info == -9999) {
                    Eigensolver<K>::_nu = std::ceil(Eigensolver<K>::_nu / 3);
                    std::fill_n(iparam + 4, 7, 0);
                    iparam[2] = _it, iparam[6] = 3;
                    ncv = 2 * Eigensolver<K>::_nu + 1;
                }
            } while(info == -9999 && Eigensolver<K>::_nu > 1 && (info = 0));
            if(s == nullptr)
                delete prec;
            Eigensolver<K>::_nu = iparam[4];
            if(Eigensolver<K>::_nu) {
                K* evr = new K[Eigensolver<K>::_nu];
                ev = new K*[Eigensolver<K>::_nu];
                *ev = new K[Eigensolver<K>::_n * Eigensolver<K>::_nu];
                for(unsigned short i = 1; i < Eigensolver<K>::_nu; ++i)
                    ev[i] = *ev + i * Eigensolver<K>::_n;
                int* select = new int[ncv];
                eupd(&i__1, "A", select, evr, *ev, &(Eigensolver<K>::_n), &(Wrapper<K>::d__0), workev, "G",
                     _which, &(Eigensolver<K>::_nu), &(Eigensolver<K>::_tol), vresid, &ncv, vp, iparam,
                     ipntr, workd, workl, &lworkl, rwork, &info);
                delete [] select;
                if(Eigensolver<K>::_threshold > 0.0)
                    Eigensolver<K>::selectNu(evr, communicator);
                delete [] evr;
            }
            else {
                ev = new K*[0];
                *ev = nullptr;
            }
            delete [] workd;
            delete [] rwork;
        }
};

#define HPDDM_GENERATE_ARPACK(C, T, B, U)                                                                    \
template<>                                                                                                   \
inline void Arpack<U>::aupd(int* ido, const char* bmat, const int* n, const char* which, const int* nu,      \
                            const U* tol, U* vresid, int* ncv, U* vp, int* iparam, int* ipntr, U* workd,     \
                            U* workl, int* lworkl, U*, int* info) {                                          \
    HPDDM_F77(B ## saupd)(ido, bmat, n, which, nu, tol, vresid, ncv, vp, n, iparam,                          \
                          ipntr, workd, workl, lworkl, info, 1, 2);                                          \
}                                                                                                            \
template<>                                                                                                   \
inline void Arpack<U>::eupd(const int* rvec, const char* HowMny, int* select, U* evr, U* ev, const int* n,   \
                            const U* sigma, U*, const char* bmat, const char* which, const int* nu,          \
                            const U* tol, U* vresid, int* necv, U* vp, int* iparam, int* ipntr,              \
                            U* workd, U* workl, int* lworkl, U*, int* info) {                                \
    HPDDM_F77(B ## seupd)(rvec, HowMny, select, evr, ev, n, sigma, bmat,                                     \
                          n, which, nu, tol, vresid, necv, vp, n, iparam,                                    \
                          ipntr, workd, workl, lworkl, info, 1, 1, 2);                                       \
}                                                                                                            \
template<>                                                                                                   \
inline void Arpack<T>::aupd(int* ido, const char* bmat, const int* n, const char* which, const int* nu,      \
                            const U* tol, T* vresid, int* ncv, T* vp, int* iparam, int* ipntr, T* workd,     \
                            T* workl, int* lworkl, U* rwork, int* info) {                                    \
    HPDDM_F77(C ## naupd)(ido, bmat, n, which, nu, tol, vresid, ncv, vp, n, iparam,                          \
                          ipntr, workd, workl, lworkl, rwork, info, 1, 2);                                   \
}                                                                                                            \
template<>                                                                                                   \
inline void Arpack<T>::eupd(const int* rvec, const char* HowMny, int* select, T* evr, T* ev, const int* n,   \
                            const T* sigma, T* workev, const char* bmat, const char* which, const int* nu,   \
                            const U* tol, T* vresid, int* necv, T* vp, int* iparam, int* ipntr,              \
                            T* workd, T* workl, int* lworkl, U* rwork, int* info) {                          \
    HPDDM_F77(C ## neupd)(rvec, HowMny, select, evr, ev, n, sigma, workev, bmat,                             \
                          n, which, nu, tol, vresid, necv, vp, n, iparam,                                    \
                          ipntr, workd, workl, lworkl, rwork, info, 1, 1, 2);                                \
}
HPDDM_GENERATE_ARPACK(c, std::complex<float>, s, float)
HPDDM_GENERATE_ARPACK(z, std::complex<double>, d, double)
} // HPDDM
#endif // _HPDDM_ARPACK_
