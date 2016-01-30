/*
   This file is part of HPDDM.

   Author(s): Pierre Jolivet <pierre.jolivet@enseeiht.fr>
        Date: 2015-08-26

   Copyright (C) 2015      Eidgenössische Technische Hochschule Zürich
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

#ifndef _HPDDM_SCALAPACK_
#define _HPDDM_SCALAPACK_

#define HPDDM_GENERATE_EXTERN_SCALAPACK_COMPLEX(C, T, B, U)                                                  \
void HPDDM_F77(p ## B ## geqpf)(const int*, const int*, U*, const int*, const int*, const int*, int*, U*,    \
                                U*, const int*, int*);                                                       \
void HPDDM_F77(p ## C ## geqpf)(const int*, const int*, T*, const int*, const int*, const int*, int*, T*,    \
                                T*, const int*, U*, const int*, int*);                                       \
void HPDDM_F77(p ## B ## geqrf)(const int*, const int*, U*, const int*, const int*, const int*, U*, U*,      \
                                const int*, int*);                                                           \
void HPDDM_F77(p ## C ## geqrf)(const int*, const int*, T*, const int*, const int*, const int*, T*, T*,      \
                                const int*, U*, const int*, int*);                                           \
void HPDDM_F77(p ## B ## orgqr)(const int*, const int*, const int*, U*, const int*, const int*, const int*,  \
                                const U*, U*, const int*, int*);                                             \
void HPDDM_F77(p ## C ## ungqr)(const int*, const int*, const int*, T*, const int*, const int*, const int*,  \
                                const T*, T*, const int*, int*);

extern "C" {
void Cblacs_get(int, int, int*);
void Cblacs_gridinit(int*, const char*, int, int);
void Cblacs_gridexit(int);
void descinit_(int*, const int*, const int*, const int*, const int*, const int*, const int*, const int*, const int*, int*);
HPDDM_GENERATE_EXTERN_SCALAPACK_COMPLEX(c, std::complex<float>, s, float)
HPDDM_GENERATE_EXTERN_SCALAPACK_COMPLEX(z, std::complex<double>, d, double)
}

namespace HPDDM {
/* Class: ScaLapack
 *
 * Template Parameter:
 *    K              - Scalar type. */
template<class K>
class ScaLapack {
    public:
        ScaLapack() = delete;
        ScaLapack(const ScaLapack&) = delete;
        template<bool pivoting = true>
        static void workspace(const int* m, const int* n, const int* ia, const int* ja, const int* desca, int* lwork, int* lrwork, int* info) {
            *lwork = -1;
            K wkopt;
            if(!Wrapper<K>::is_complex) {
                if(pivoting)
                    geqpf(m, n, nullptr, ia, ja, desca, nullptr, nullptr, &wkopt, lwork, nullptr, nullptr, info);
                else
                    geqrf(m, n, nullptr, ia, ja, desca, nullptr, &wkopt, lwork, nullptr, nullptr, info);
                *lwork = static_cast<int>(std::real(wkopt));
                if(lrwork)
                    *lrwork = 0;
            }
            else {
                *lrwork = -1;
                underlying_type<K> rwkopt;
                if(pivoting)
                    geqpf(m, n, nullptr, ia, ja, desca, nullptr, nullptr, &wkopt, lwork, &rwkopt, lrwork, info);
                else
                    geqrf(m, n, nullptr, ia, ja, desca, nullptr, &wkopt, lwork, &rwkopt, lrwork, info);
                *lwork = static_cast<int>(std::real(wkopt));
                *lrwork = static_cast<int>(rwkopt);
            }
        }
        static void geqpf(const int*, const int*, K*, const int*, const int*, const int*, int*, K*, K*, const int*, underlying_type<K>*, const int*, int*);
        static void geqrf(const int*, const int*, K*, const int*, const int*, const int*, K*, K*, const int*, underlying_type<K>*, const int*, int*);
        static void gqr(const int*, const int*, const int*, K*, const int*, const int*, const int*, const K*, K*, const int*, int*);
};

#define HPDDM_GENERATE_SCALAPACK_COMPLEX(C, T, B, U)                                                         \
template<>                                                                                                   \
inline void ScaLapack<U>::geqpf(const int* m, const int* n, U* a, const int* ia, const int* ja,              \
                                const int* desca, int* jpvt, U* tau, U* work, const int* lwork, U* rwork,    \
                                const int* lrwork, int* info) {                                              \
    HPDDM_F77(p ## B ## geqpf)(m, n, a, ia, ja, desca, jpvt, tau, work, lwork, info);                        \
}                                                                                                            \
template<>                                                                                                   \
inline void ScaLapack<T>::geqpf(const int* m, const int* n, T* a, const int* ia, const int* ja,              \
                                const int* desca, int* jpvt, T* tau, T* work, const int* lwork, U* rwork,    \
                                const int* lrwork, int* info) {                                              \
    HPDDM_F77(p ## C ## geqpf)(m, n, a, ia, ja, desca, jpvt, tau, work, lwork, rwork, lrwork, info);         \
}                                                                                                            \
template<>                                                                                                   \
inline void ScaLapack<U>::geqrf(const int* m, const int* n, U* a, const int* ia, const int* ja,              \
                                const int* desca, U* tau, U* work, const int* lwork, U* rwork,               \
                                const int* lrwork, int* info) {                                              \
    HPDDM_F77(p ## B ## geqrf)(m, n, a, ia, ja, desca, tau, work, lwork, info);                              \
}                                                                                                            \
template<>                                                                                                   \
inline void ScaLapack<T>::geqrf(const int* m, const int* n, T* a, const int* ia, const int* ja,              \
                                const int* desca, T* tau, T* work, const int* lwork, U* rwork,               \
                                const int* lrwork, int* info) {                                              \
    HPDDM_F77(p ## C ## geqrf)(m, n, a, ia, ja, desca, tau, work, lwork, rwork, lrwork, info);               \
}                                                                                                            \
template<>                                                                                                   \
inline void ScaLapack<U>::gqr(const int* m, const int* n, const int* k, U* a, const int* ia, const int* ja,  \
                              const int* desca, const U* tau, U* work, const int* lwork, int* info) {        \
    HPDDM_F77(p ## B ## orgqr)(m, n, k, a, ia, ja, desca, tau, work, lwork, info);                           \
}                                                                                                            \
template<>                                                                                                   \
inline void ScaLapack<T>::gqr(const int* m, const int* n, const int* k, T* a, const int* ia, const int* ja,  \
                              const int* desca, const T* tau, T* work, const int* lwork, int* info) {        \
    HPDDM_F77(p ## C ## ungqr)(m, n, k, a, ia, ja, desca, tau, work, lwork, info);                           \
}
HPDDM_GENERATE_SCALAPACK_COMPLEX(c, std::complex<float>, s, float)
HPDDM_GENERATE_SCALAPACK_COMPLEX(z, std::complex<double>, d, double)
} // HPDDM
#endif // _HPDDM_SCALAPACK_
