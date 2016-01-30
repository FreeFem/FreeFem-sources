 /*
   This file is part of HPDDM.

   Author(s): Pierre Jolivet <pierre.jolivet@enseeiht.fr>
        Date: 2014-11-05

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

#ifndef _HPDDM_ITERATIVE_
#define _HPDDM_ITERATIVE_

#include "ScaLAPACK.hpp"

namespace HPDDM {
/* Class: Iterative method
 *  A class that implements various iterative methods. */
class IterativeMethod {
    private:
        /* Function: allocate
         *  Allocates workspace arrays for <Iterative method::CG>. */
        template<class K, typename std::enable_if<!Wrapper<K>::is_complex>::type* = nullptr>
        static void allocate(K*& dir, K*& p, const int& n, const unsigned short extra = 0, const unsigned short it = 1) {
            if(extra == 0) {
                dir = new K[2 + std::max(1, 4 * n)];
                p = dir + 2;
            }
            else {
                dir = new K[1 + 2 * it + std::max(1, (4 + extra * it) * n)];
                p = dir + 1 + 2 * it;
            }
        }
        template<class K, typename std::enable_if<Wrapper<K>::is_complex>::type* = nullptr>
        static void allocate(underlying_type<K>*& dir, K*& p, const int& n, const unsigned short extra = 0, const unsigned short it = 1) {
            if(extra == 0) {
                dir = new underlying_type<K>[2];
                p = new K[std::max(1, 4 * n)];
            }
            else {
                dir = new underlying_type<K>[1 + 2 * it];
                p = new K[std::max(1, (4 + extra * it) * n)];
            }
        }
        /* Function: depenalize
         *  Divides a scalar by <HPDDM_PEN>. */
        template<class K, typename std::enable_if<!Wrapper<K>::is_complex>::type* = nullptr>
        static void depenalize(const K& b, K& x) { x = b / HPDDM_PEN; }
        template<class K, typename std::enable_if<Wrapper<K>::is_complex>::type* = nullptr>
        static void depenalize(const K& b, K& x) { x = b / std::complex<underlying_type<K>>(HPDDM_PEN, HPDDM_PEN); }
        /* Function: update
         *
         *  Updates a solution vector after convergence of <Iterative method::GMRES>.
         *
         * Template Parameter:
         *    K              - Scalar type.
         *
         * Parameters:
         *    variant        - Type of preconditioning.
         *    n              - Size of the vector.
         *    x              - Solution vector.
         *    k              - Dimension of the Hessenberg matrix.
         *    h              - Hessenberg matrix.
         *    s              - Coefficients in the Krylov subspace.
         *    v              - Basis of the Krylov subspace. */
        template<class Operator, class K>
        static void update(const Operator& A, char variant, const int& n, K* const x, const K* const* const h, K* const s, const K* const* const v, const short* const hasConverged, const int& mu, K* const work, const int& deflated = -1) {
            int tmp = std::distance(h[0], h[1]);
            if(mu == 1 || deflated != -1) {
                int dim = std::abs(*hasConverged);
                int info;
                if(deflated != -1)
                    tmp /= deflated;
                Lapack<K>::trtrs("U", "N", "N", &dim, deflated != -1 ? &deflated : &mu, *h, &tmp, s, &tmp, &info);
            }
            else
                for(unsigned short nu = 0; nu < mu; ++nu) {
                    for(int i = std::abs(hasConverged[nu]); i-- > 0; ) {
                        K alpha = -(s[i * mu + nu] /= h[i][i * mu + nu]);
                        Blas<K>::axpy(&i, &alpha, h[i] + nu, &mu, s + nu, &mu);
                    }
                }
            K* const correction = (variant == 'R' ? const_cast<K*>(v[tmp / (deflated == -1 ? mu : deflated) - 1]) : work);
            if(deflated == -1) {
                tmp = mu * n;
                if(variant == 'L') {
                    for(unsigned short nu = 0; nu < mu; ++nu)
                        if(hasConverged[nu] != 0) {
                            int dim = std::abs(hasConverged[nu]);
                            Blas<K>::gemv("N", &n, &dim, &(Wrapper<K>::d__1), *v + nu * n, &tmp, s + nu, &mu, &(Wrapper<K>::d__1), x + nu * n, &i__1);
                        }
                }
                else {
                    for(unsigned short nu = 0; nu < mu; ++nu) {
                        int dim = std::abs(hasConverged[nu]);
                        Blas<K>::gemv("N", &n, &dim, &(Wrapper<K>::d__1), *v + nu * n, &tmp, s + nu, &mu, &(Wrapper<K>::d__0), work + nu * n, &i__1);
                    }
                    if(variant == 'R')
                        A.apply(work, correction, mu);
                    for(unsigned short nu = 0; nu < mu; ++nu)
                        if(hasConverged[nu] != 0)
                            Blas<K>::axpy(&n, &(Wrapper<K>::d__1), correction + nu * n, &i__1, x + nu * n, &i__1);
                }
            }
            else if(deflated == mu) {
                int dim = *hasConverged;
                if(variant == 'L')
                    Blas<K>::gemm("N", "N", &n, &mu, &dim, &(Wrapper<K>::d__1), *v, &n, s, &tmp, &(Wrapper<K>::d__1), x, &n);
                else {
                    Blas<K>::gemm("N", "N", &n, &mu, &dim, &(Wrapper<K>::d__1), *v, &n, s, &tmp, &(Wrapper<K>::d__0), work, &n);
                    if(variant == 'R')
                        A.apply(work, correction, mu);
                    Blas<K>::axpy(&(tmp = mu * n), &(Wrapper<K>::d__1), correction, &i__1, x, &i__1);
                }
            }
            else {
                int dim = *hasConverged;
                Blas<K>::gemm("N", "N", &n, &deflated, &dim, &(Wrapper<K>::d__1), *v, &n, s, &tmp, &(Wrapper<K>::d__0), work, &n);
                if(variant == 'R')
                    A.apply(work, correction, deflated);
                Blas<K>::gemm("N", "N", &n, &(dim = mu - deflated), &deflated, &(Wrapper<K>::d__1), correction, &n, s + deflated * tmp, &tmp, &(Wrapper<K>::d__1), x + deflated * n, &n);
                Blas<K>::axpy(&(tmp = deflated * n), &(Wrapper<K>::d__1), correction, &i__1, x, &i__1);
            }
        }
        template<class T, typename std::enable_if<std::is_pointer<T>::value>::type* = nullptr>
        static void clean(T* const& pt) {
            delete [] *pt;
            delete []  pt;
        }
        template<class T, typename std::enable_if<!std::is_pointer<T>::value>::type* = nullptr>
        static void clean(T* const& pt) {
            delete [] pt;
        }
        template<class K, class T, typename std::enable_if<std::is_pointer<T>::value>::type* = nullptr>
        static void axpy(const int* const n, const K* const a, const T* const x, const int* const incx, T* const y, const int* const incy) {
            static_assert(std::is_same<typename std::remove_pointer<T>::type, K>::value, "Wrong types");
            Blas<typename std::remove_pointer<T>::type>::axpy(n, a, *x, incx, *y, incy);
        }
        template<class K, class T, typename std::enable_if<std::is_pointer<T>::value>::type* = nullptr>
        static void axpy(const int* const, const K* const, const T* const, const int* const, T const, const int* const) { }
        template<class K, class T, typename std::enable_if<!std::is_pointer<T>::value>::type* = nullptr>
        static void axpy(const int* const n, const K* const a, const T* const x, const int* const incx, T* const y, const int* const incy) {
            static_assert(std::is_same<T, K>::value, "Wrong types");
            Blas<T>::axpy(n, a, x, incx, y, incy);
        }
        template<class T, typename std::enable_if<std::is_pointer<T>::value>::type* = nullptr>
        static typename std::remove_pointer<T>::type dot(const int* const n, const T* const x, const int* const incx, const T* const y, const int* const incy) {
            return Blas<typename std::remove_pointer<T>::type>::dot(n, *x, incx, *y, incy) / 2.0;
        }
        template<class T, typename std::enable_if<!std::is_pointer<T>::value>::type* = nullptr>
        static T dot(const int* const n, const T* const x, const int* const incx, const T* const y, const int* const incy) {
            return Blas<T>::dot(n, x, incx, y, incy);
        }
        template<class T, class U, typename std::enable_if<std::is_pointer<T>::value>::type* = nullptr>
        static void diag(const int&, const U* const* const, T* const, T* const = nullptr) { }
        template<class T, typename std::enable_if<!std::is_pointer<T>::value>::type* = nullptr>
        static void diag(const int& n, const underlying_type<T>* const d, T* const in, T* const out = nullptr) {
            if(out)
                Wrapper<T>::diag(n, d, in, out);
            else
                Wrapper<T>::diag(n, d, in);
        }
    public:
        /* Function: GMRES
         *
         *  Implements the GMRES.
         *
         * Template Parameters:
         *    Type           - See <Gmres>.
         *    excluded       - True if the master processes are excluded from the domain decomposition, false otherwise.
         *    K              - Scalar type.
         *
         * Parameters:
         *    A              - Global operator.
         *    x              - Solution vector(s).
         *    b              - Right-hand side(s).
         *    mu             - Number of right-hand sides.
         *    comm           - Global MPI communicator. */
        template<bool excluded = false, class Operator = void, class K = double>
        static int GMRES(const Operator& A, K* const x, const K* const b, const int& mu, const MPI_Comm& comm) {
            const Option& opt = *Option::get();
            const int n = excluded ? 0 : A.getDof();
            const unsigned short it = opt["max_it"];
            underlying_type<K> tol = opt["tol"];
            const char verbosity = opt.val<char>("verbosity");
            if(std::abs(tol) < std::numeric_limits<underlying_type<K>>::epsilon()) {
                if(verbosity > 0)
                    std::cout << "WARNING -- the tolerance of the iterative method was set to " << tol << " which is lower than the machine epsilon for type " << demangle(typeid(underlying_type<K>).name()) << ", forcing the tolerance to " << 2 * std::numeric_limits<underlying_type<K>>::epsilon() << std::endl;
                tol = 2 * std::numeric_limits<underlying_type<K>>::epsilon();
            }
            const unsigned short m = std::min(static_cast<unsigned short>(std::numeric_limits<short>::max()), std::min(static_cast<unsigned short>(opt["gmres_restart"]), it));
            const char variant = (opt["variant"] == 0 ? 'L' : opt["variant"] == 1 ? 'R' : 'F');

            K** const H = new K*[m * (2 + (variant == 'F')) + 1];
            K** const v = H + m;
            K* const s = new K[mu * ((m + 1) * (m + 1) + n * (2 + m * (1 + (variant == 'F'))) + (!Wrapper<K>::is_complex ? m + 1 : (m + 2) / 2))];
            K* const Ax = s + mu * (m + 1);
            *H = Ax + mu * n;
            for(unsigned short i = 1; i < m; ++i)
                H[i] = *H + i * mu * (m + 1);
            *v = *H + m * mu * (m + 1);
            for(unsigned short i = 1; i < m * (1 + (variant == 'F')) + 1; ++i)
                v[i] = *v + i * mu * n;
            underlying_type<K>* const norm = reinterpret_cast<underlying_type<K>*>(*v + (m * (1 + (variant == 'F')) + 1) * mu * n);
            underlying_type<K>* const sn = norm + mu;
            bool alloc = A.setBuffer(mu);
            short* const hasConverged = new short[mu];
            std::fill_n(hasConverged, mu, -m);

            A.template start<excluded>(b, x, mu);
            A.template apply<excluded>(b, *v, mu, Ax);
            for(unsigned short nu = 0; nu < mu; ++nu)
                norm[nu] = std::real(Blas<K>::dot(&n, *v + nu * n, &i__1, *v + nu * n, &i__1));

            if(!excluded)
                for(unsigned short nu = 0; nu < mu; ++nu)
                    for(unsigned int i = 0; i < n; ++i)
                        if(std::abs(b[nu * n + i]) > HPDDM_PEN * HPDDM_EPS)
                            depenalize(b[nu * n + i], x[nu * n + i]);

            unsigned short j = 1;
            while(j <= it) {
                if(!excluded)
                    A.GMV(x, variant == 'L' ? Ax : *v, mu);
                Blas<K>::axpby(mu * n, 1.0, b, 1, -1.0, variant == 'L' ? Ax : *v, 1);
                if(variant == 'L')
                    A.template apply<excluded>(Ax, *v, mu);
                for(unsigned short nu = 0; nu < mu; ++nu)
                    sn[nu] = std::real(Blas<K>::dot(&n, *v + nu * n, &i__1, *v + nu * n, &i__1));
                if(j == 1) {
                    MPI_Allreduce(MPI_IN_PLACE, norm, 2 * mu, Wrapper<K>::mpi_underlying_type(), MPI_SUM, comm);
                    for(unsigned short nu = 0; nu < mu; ++nu) {
                        norm[nu] = std::sqrt(norm[nu]);
                        if(norm[nu] < HPDDM_EPS)
                            norm[nu] = 1.0;
                        if(tol > 0.0 && sn[nu] / norm[nu] < tol) {
                            if(norm[nu] > 1.0 / HPDDM_EPS)
                                norm[nu] = 1.0;
                            else
                                hasConverged[nu] = 0;
                        }
                        else if(sn[nu] < -tol)
                            hasConverged[nu] = 0;
                    }
                }
                else
                    MPI_Allreduce(MPI_IN_PLACE, sn, mu, Wrapper<K>::mpi_underlying_type(), MPI_SUM, comm);
                for(unsigned short nu = 0; nu < mu; ++nu) {
                    if(hasConverged[nu] > 0)
                        hasConverged[nu] = 0;
                    s[nu] = std::sqrt(sn[nu]);
                    std::for_each(*v + nu * n, *v + (nu + 1) * n, [&](K& y) { y /= s[nu]; });
                }
                int i = 0;
                while(i < m && j <= it) {
                    if(variant == 'L') {
                        if(!excluded)
                            A.GMV(v[i], Ax, mu);
                        A.template apply<excluded>(Ax, v[i + 1], mu);
                    }
                    else {
                        A.template apply<excluded>(v[i], variant == 'F' ? v[i + m + 1] : Ax, mu, v[i + 1]);
                        if(!excluded)
                            A.GMV(variant == 'F' ? v[i + m + 1] : Ax, v[i + 1], mu);
                    }
                    if(excluded) {
                        std::fill_n(H[i], mu * (i + 1), K());
                        if(opt["gs"] == 1)
                            for(int k = 0; k < i + 1; ++k)
                                MPI_Allreduce(MPI_IN_PLACE, H[i] + mu * k, mu, Wrapper<K>::mpi_type(), MPI_SUM, comm);
                        else
                            MPI_Allreduce(MPI_IN_PLACE, H[i], mu * (i + 1), Wrapper<K>::mpi_type(), MPI_SUM, comm);
                        std::fill_n(sn + i * mu, mu, underlying_type<K>());
                        MPI_Allreduce(MPI_IN_PLACE, sn + i * mu, mu, Wrapper<K>::mpi_underlying_type(), MPI_SUM, comm);
                        std::transform(sn + i * mu, sn + (i + 1) * mu, H[i] + (i + 1) * mu, [](const underlying_type<K>& b) { return std::sqrt(b); });
                    }
                    else {
                        if(opt["gs"] == 1)
                            for(unsigned short k = 0; k < i + 1; ++k) {
                                for(unsigned short nu = 0; nu < mu; ++nu)
                                    H[i][k * mu + nu] = Blas<K>::dot(&n, v[k] + nu * n, &i__1, v[i + 1] + nu * n, &i__1);
                                MPI_Allreduce(MPI_IN_PLACE, H[i] + k * mu, mu, Wrapper<K>::mpi_type(), MPI_SUM, comm);
                                std::transform(H[i] + k * mu, H[i] + (k + 1) * mu, H[i] + (i + 1) * mu, [](const K& h) { return -h; });
                                for(unsigned short nu = 0; nu < mu; ++nu)
                                    Blas<K>::axpy(&n, H[i] + (i + 1) * mu + nu, v[k] + nu * n, &i__1, v[i + 1] + nu * n, &i__1);
                            }
                        else {
                            int tmp[2] { i + 1, mu * n };
                            for(unsigned short nu = 0; nu < mu; ++nu)
                                Blas<K>::gemv(&(Wrapper<K>::transc), &n, tmp, &(Wrapper<K>::d__1), *v + nu * n, tmp + 1, v[i + 1] + nu * n, &i__1, &(Wrapper<K>::d__0), H[i] + nu, &mu);
                            MPI_Allreduce(MPI_IN_PLACE, H[i], (i + 1) * mu, Wrapper<K>::mpi_type(), MPI_SUM, comm);
                            if(opt["gs"] == 0)
                                for(unsigned short nu = 0; nu < mu; ++nu)
                                    Blas<K>::gemv("N", &n, tmp, &(Wrapper<K>::d__2), *v + nu * n, tmp + 1, H[i] + nu, &mu, &(Wrapper<K>::d__1), v[i + 1] + nu * n, &i__1);
                            else
                                for(unsigned short nu = 0; nu < mu; ++nu)
                                    Blas<K>::axpby(n, -H[i][i * mu + nu], v[i] + nu * n, 1, 1.0, v[i + 1] + nu * n, 1);
                        }
                        for(unsigned short nu = 0; nu < mu; ++nu)
                            sn[i * mu + nu] = std::real(Blas<K>::dot(&n, v[i + 1] + nu * n, &i__1, v[i + 1] + nu * n, &i__1));
                        MPI_Allreduce(MPI_IN_PLACE, sn + i * mu, mu, Wrapper<K>::mpi_underlying_type(), MPI_SUM, comm);
                        for(unsigned short nu = 0; nu < mu; ++nu) {
                            H[i][(i + 1) * mu + nu] = std::sqrt(sn[i * mu + nu]);
                            if(i < m - 1)
                                std::for_each(v[i + 1] + nu * n, v[i + 1] + (nu + 1) * n, [&](K& y) { y /= H[i][(i + 1) * mu + nu]; });
                        }
                    }
                    for(unsigned short k = 0; k < i; ++k) {
                        for(unsigned short nu = 0; nu < mu; ++nu) {
                            K gamma = Wrapper<K>::conj(H[k][(k + 1) * mu + nu]) * H[i][k * mu + nu] + sn[k * mu + nu] * H[i][(k + 1) * mu + nu];
                            H[i][(k + 1) * mu + nu] = -sn[k * mu + nu] * H[i][k * mu + nu] + H[k][(k + 1) * mu + nu] * H[i][(k + 1) * mu + nu];
                            H[i][k * mu + nu] = gamma;
                        }
                    }
                    for(unsigned short nu = 0; nu < mu; ++nu) {
                        const int tmp = 2;
                        underlying_type<K> delta = Blas<K>::nrm2(&tmp, H[i] + i * mu + nu, &mu); // std::sqrt(H[i][i] * H[i][i] + H[i][i + 1] * H[i][i + 1]);
                        sn[i * mu + nu] = std::real(H[i][(i + 1) * mu + nu]) / delta;
                        H[i][(i + 1) * mu + nu] = H[i][i * mu + nu] / delta;
                        H[i][i * mu + nu] = delta;
                        s[(i + 1) * mu + nu] = -sn[i * mu + nu] * s[i * mu + nu];
                        s[i * mu + nu] *= Wrapper<K>::conj(H[i][(i + 1) * mu + nu]);
                        if(hasConverged[nu] == -m && ((tol > 0 && std::abs(s[(i + 1) * mu + nu]) / norm[nu] <= tol) || (tol < 0 && std::abs(s[(i + 1) * mu + nu]) <= -tol)))
                            hasConverged[nu] = i + 1;
                    }
                    if(verbosity > 0) {
                        int tmp[2] { 0, 0 };
                        underlying_type<K> beta = std::abs(s[(i + 1) * mu]);
                        for(unsigned short nu = 0; nu < mu; ++nu) {
                            if(hasConverged[nu] != -m)
                                ++tmp[0];
                            else if(std::abs(s[(i + 1) * mu + nu]) > beta) {
                                beta = std::abs(s[(i + 1) * mu + nu]);
                                tmp[1] = nu;
                            }
                        }
                        if(tol > 0)
                            std::cout << "GMRES: " << std::setw(3) << j << " " << std::scientific << beta << " " <<  norm[tmp[1]] << " " <<  beta / norm[tmp[1]] << " < " << tol;
                        else
                            std::cout << "GMRES: " << std::setw(3) << j << " " << std::scientific << beta << " < " << -tol;
                        if(mu > 1) {
                            std::cout << " (rhs #" << tmp[1] + 1;
                            if(tmp[0] > 0)
                                std::cout << ", " << tmp[0] << " converged rhs";
                            std::cout << ")";
                        }
                        std::cout << std::endl;
                    }
                    if(std::find(hasConverged, hasConverged + mu, -m) == hasConverged + mu)
                        break;
                    else
                        ++i, ++j;
                }
                if(j != it + 1 && i == m) {
                    if(!excluded)
                        update(A, variant, n, x, H, s, v + (m + 1) * (variant == 'F'), hasConverged, mu, Ax);
                    std::fill_n(s + mu, mu * m, K());
                    if(verbosity > 0)
                        std::cout << "GMRES restart(" << m << ")" << std::endl;
                }
                else
                    break;
            }
            if(!excluded)
                update(A, variant, n, x, H, s, v + (m + 1) * (variant == 'F'), hasConverged, mu, Ax);
            if(verbosity > 0) {
                if(j != it + 1)
                    std::cout << "GMRES converges after " << j << " iteration" << (j > 1 ? "s" : "") << std::endl;
                else
                    std::cout << "GMRES does not converges after " << it << " iteration" << (it > 1 ? "s" : "") << std::endl;
            }
            delete [] hasConverged;
            A.clearBuffer(alloc);
            delete [] s;
            delete [] H;
            return std::min(j, it);
        }
        /* Function: CG
         *
         *  Implements the CG method.
         *
         * Template Parameters:
         *    excluded       - True if the master processes are excluded from the domain decomposition, false otherwise.
         *    K              - Scalar type.
         *
         * Parameters:
         *    A              - Global operator.
         *    x              - Solution vector.
         *    b              - Right-hand side.
         *    comm           - Global MPI communicator. */
        template<bool excluded = false, class Operator, class K>
        static int CG(const Operator& A, K* const x, const K* const b, const MPI_Comm& comm) {
            const Option& opt = *Option::get();
            if(opt.any_of("schwarz_method", { 0, 1, 4 }) || opt.any_of("schwarz_coarse_correction", { 0 }))
                return GMRES(A, x, b, 1, comm);
            const int n = A.getDof();
            const unsigned short it = opt["max_it"];
            underlying_type<K> tol = opt["tol"];
            const char verbosity = opt.val<char>("verbosity");
            if(std::abs(tol) < std::numeric_limits<underlying_type<K>>::epsilon()) {
                if(verbosity > 0)
                    std::cout << "WARNING -- the tolerance of the iterative method was set to " << tol << " which is lower than the machine epsilon for type " << demangle(typeid(underlying_type<K>).name()) << ", forcing the tolerance to " << 2 * std::numeric_limits<underlying_type<K>>::epsilon() << std::endl;
                tol = 2 * std::numeric_limits<underlying_type<K>>::epsilon();
            }
            underlying_type<K>* dir;
            K* trash;
            allocate(dir, trash, n, opt["variant"] == 2 ? 2 : (opt["gs"] != 2 ? 1 : 0), it);
            bool alloc = A.setBuffer(1);
            K* z = trash + n;
            K* r = z + n;
            K* p = r + n;
            const underlying_type<K>* const d = A.getScaling();

            A.template start<excluded>(b, x);
            for(unsigned int i = 0; i < n; ++i)
                if(std::abs(b[i]) > HPDDM_PEN * HPDDM_EPS)
                    depenalize(b[i], x[i]);
            A.GMV(x, z);
            std::copy_n(b, n, r);
            Blas<K>::axpy(&n, &(Wrapper<K>::d__2), z, &i__1, r, &i__1);

            A.apply(r, p, 1, z);

            Wrapper<K>::diag(n, d, p, trash);
            dir[0] = std::real(Blas<K>::dot(&n, trash, &i__1, p, &i__1));
            MPI_Allreduce(MPI_IN_PLACE, dir, 1, Wrapper<K>::mpi_underlying_type(), MPI_SUM, comm);
            underlying_type<K> resInit = std::sqrt(dir[0]);

            unsigned short i = 1;
            while(i <= it) {
                dir[0] = std::real(Blas<K>::dot(&n, r, &i__1, trash, &i__1));
                if(opt["variant"] == 2 && i > 1) {
                    for(unsigned short k = 0; k < i - 1; ++k)
                        dir[1 + k] = -std::real(Blas<K>::dot(&n, trash, &i__1, p + (1 + it + k) * n, &i__1)) / dir[1 + it + k];
                    MPI_Allreduce(MPI_IN_PLACE, dir + 1, i - 1, Wrapper<K>::mpi_underlying_type(), MPI_SUM, comm);
                    std::copy_n(z, n, p);
                    for(unsigned short k = 0; k < i - 1; ++k) {
                        trash[0] = dir[1 + k];
                        Blas<K>::axpy(&n, trash, p + (1 + k) * n, &i__1, p, &i__1);
                    }
                }
                A.GMV(p, z);
                if(opt["gs"] != 2 && i > 1) {
                    Wrapper<K>::diag(n, d, z, trash);
                    for(unsigned short k = 0; k < i - 1; ++k)
                        dir[1 + k] = -std::real(Blas<K>::dot(&n, trash, &i__1, p + (1 + k) * n, &i__1)) / dir[1 + it + k];
                    MPI_Allreduce(MPI_IN_PLACE, dir + 1, i - 1, Wrapper<K>::mpi_underlying_type(), MPI_SUM, comm);
                    for(unsigned short k = 0; k < i - 1; ++k) {
                        trash[0] = dir[1 + k];
                        Blas<K>::axpy(&n, trash, p + (1 + k) * n, &i__1, p, &i__1);
                    }
                    A.GMV(p, z);
                }
                Wrapper<K>::diag(n, d, p, trash);
                dir[1] = std::real(Blas<K>::dot(&n, z, &i__1, trash, &i__1));
                MPI_Allreduce(MPI_IN_PLACE, dir, 2, Wrapper<K>::mpi_underlying_type(), MPI_SUM, comm);
                if(opt["gs"] != 2 || opt["variant"] == 2) {
                    dir[it + i] = dir[1];
                    std::copy_n(p, n, p + i * n);
                    if(opt["variant"] == 2)
                        std::copy_n(z, n, p + (it + i) * n);
                }
                trash[0] = dir[0] / dir[1];
                Blas<K>::axpy(&n, trash, p, &i__1, x, &i__1);
                trash[0] = -trash[0];
                Blas<K>::axpy(&n, trash, z, &i__1, r, &i__1);

                A.apply(r, z, 1, trash);
                Wrapper<K>::diag(n, d, z, trash);
                dir[1] = std::real(Blas<K>::dot(&n, r, &i__1, trash, &i__1)) / dir[0];
                dir[0] = std::real(Blas<K>::dot(&n, z, &i__1, trash, &i__1));
                MPI_Allreduce(MPI_IN_PLACE, dir, 2, Wrapper<K>::mpi_underlying_type(), MPI_SUM, comm);
                if(opt["variant"] != 2)
                    Blas<K>::axpby(n, 1.0, z, 1, dir[1], p, 1);
                dir[0] = std::sqrt(dir[0]);
                if(verbosity > 0) {
                    if(tol > 0)
                        std::cout << "CG: " << std::setw(3) << i << " " << std::scientific << dir[0] << " " << resInit << " " << dir[0] / resInit << " < " << tol << std::endl;
                    else
                        std::cout << "CG: " << std::setw(3) << i << " " << std::scientific << dir[0] << " < " << -tol << std::endl;
                }
                if((tol > 0.0 && dir[0] / resInit <= tol) || (tol < 0.0 && dir[0] <= -tol))
                    break;
                else
                    ++i;
            }
            if(verbosity > 0) {
                if(i != it + 1)
                    std::cout << "CG converges after " << i << " iteration" << (i > 1 ? "s" : "") << std::endl;
                else
                    std::cout << "CG does not converges after " << it << " iteration" << (it > 1 ? "s" : "") << std::endl;
            }
            delete [] dir;
            if(Wrapper<K>::is_complex)
                delete [] trash;
            A.clearBuffer(alloc);
            return std::min(i, it);
        }
        /* Function: PCG
         *
         *  Implements the projected CG method.
         *
         * Template Parameters:
         *    excluded       - True if the master processes are excluded from the domain decomposition, false otherwise.
         *    K              - Scalar type.
         *
         * Parameters:
         *    A              - Global operator.
         *    x              - Solution vector.
         *    f              - Right-hand side.
         *    comm           - Global MPI communicator. */
        template<bool excluded = false, class Operator, class K>
        static int PCG(const Operator& A, K* const x, const K* const f, const MPI_Comm& comm) {
            typedef typename std::conditional<std::is_pointer<typename std::remove_reference<decltype(*A.getScaling())>::type>::value, K**, K*>::type ptr_type;
            const Option& opt = *Option::get();
            const int n = std::is_same<ptr_type, K*>::value ? A.getDof() : A.getMult();
            const unsigned short it = opt["max_it"];
            underlying_type<K> tol = opt["tol"];
            const char verbosity = opt.val<char>("verbosity");
            if(std::abs(tol) < std::numeric_limits<underlying_type<K>>::epsilon()) {
                if(verbosity > 0)
                    std::cout << "WARNING -- the tolerance of the iterative method was set to " << tol << " which is lower than the machine epsilon for type " << demangle(typeid(underlying_type<K>).name()) << ", forcing the tolerance to " << 2 * std::numeric_limits<underlying_type<K>>::epsilon() << std::endl;
                tol = 2 * std::numeric_limits<underlying_type<K>>::epsilon();
            }
            const int offset = std::is_same<ptr_type, K*>::value ? A.getEliminated() : 0;
            ptr_type storage[std::is_same<ptr_type, K*>::value ? 1 : 2];
            // storage[0] = r
            // storage[1] = lambda
            A.allocateArray(storage);
            bool alloc = A.setBuffer(1);
            auto m = A.getScaling();
            if(std::is_same<ptr_type, K*>::value)
                A.template start<excluded>(f, x + offset, nullptr, storage[0]);
            else
                A.template start<excluded>(f, x, storage[1], storage[0]);

            std::vector<ptr_type> z;
            z.reserve(it);
            ptr_type zCurr;
            A.allocateSingle(zCurr);
            z.emplace_back(zCurr);
            if(!excluded)
                A.precond(storage[0], zCurr);                                                              //     z_0 = M r_0

            underlying_type<K> resInit;
            A.template computeDot<excluded>(&resInit, zCurr, zCurr, comm);
            resInit = std::sqrt(resInit);

            std::vector<ptr_type> p;
            p.reserve(it);
            ptr_type pCurr;
            A.allocateSingle(pCurr);
            p.emplace_back(pCurr);

            K* alpha = new K[excluded ? std::max(static_cast<unsigned short>(2), it) : 2 * it];
            underlying_type<K> resRel = std::numeric_limits<underlying_type<K>>::max();
            unsigned short i = 1;
            while(i <= it) {
                if(!excluded) {
                    A.template project<excluded, 'N'>(zCurr, pCurr);                                       //     p_i = P z_i
                    for(unsigned short k = 0; k < i - 1; ++k)
                        alpha[it + k] = dot(&n, z[k], &i__1, pCurr, &i__1);
                    MPI_Allreduce(MPI_IN_PLACE, alpha + it, i - 1, Wrapper<K>::mpi_type(), MPI_SUM, comm); // alpha_k = < z_k, p_i >
                    for(unsigned short k = 0; k < i - 1; ++k) {
                        alpha[it + k] /= -alpha[k];
                        axpy(&n, alpha + it + k, p[k], &i__1, pCurr, &i__1);                               //     p_i = p_i - sum < z_k, p_i > / < z_k, p_k > p_k
                    }
                    A.apply(pCurr, zCurr);                                                                 //     z_i = F p_i

                    A.allocateSingle(zCurr);
                    if(std::is_same<ptr_type, K*>::value) {
                        diag(n, m, pCurr, zCurr);
                        alpha[i - 1] = dot(&n, z.back(), &i__1, zCurr, &i__1);
                        alpha[i]     = dot(&n, storage[0], &i__1, zCurr, &i__1);
                    }
                    else {
                        alpha[i - 1] = dot(&n, z.back(), &i__1, pCurr, &i__1);
                        alpha[i]     = dot(&n, storage[0], &i__1, pCurr, &i__1);
                    }
                    MPI_Allreduce(MPI_IN_PLACE, alpha + i - 1, 2, Wrapper<K>::mpi_type(), MPI_SUM, comm);
                    alpha[it] = alpha[i] / alpha[i - 1];
                    if(std::is_same<ptr_type, K*>::value)
                        axpy(&n, alpha + it, pCurr, &i__1, x + offset, &i__1);
                    else
                        axpy(&n, alpha + it, pCurr, &i__1, storage[1], &i__1);                             // l_i + 1 = l_i + < r_i, p_i > / < z_i, p_i > p_i
                    alpha[it] = -alpha[it];
                    axpy(&n, alpha + it, z.back(), &i__1, storage[0], &i__1);                              // r_i + 1 = r_i - < r_i, p_i > / < z_i, p_i > z_i
                    A.template project<excluded, 'T'>(storage[0]);                                         // r_i + 1 = P^T r_i + 1

                    z.emplace_back(zCurr);
                    A.precond(storage[0], zCurr);                                                          // z_i + 1 = M r_i
                }
                else {
                    A.template project<excluded, 'N'>(zCurr, pCurr);
                    std::fill_n(alpha, i - 1, K());
                    MPI_Allreduce(MPI_IN_PLACE, alpha, i - 1, Wrapper<K>::mpi_type(), MPI_SUM, comm);
                    std::fill_n(alpha, 2, K());
                    MPI_Allreduce(MPI_IN_PLACE, alpha, 2, Wrapper<K>::mpi_type(), MPI_SUM, comm);
                    A.template project<excluded, 'T'>(storage[0]);
                }
                A.template computeDot<excluded>(&resRel, zCurr, zCurr, comm);
                resRel = std::sqrt(resRel);
                if(verbosity > 0)
                    std::cout << "CG: " << std::setw(3) << i << " " << std::scientific << resRel << " " << resInit << " " << resRel / resInit << " < " << tol << std::endl;
                if(resRel / resInit <= tol)
                    break;
                else
                    ++i;
                if(!excluded) {
                    A.allocateSingle(pCurr);
                    p.emplace_back(pCurr);
                    diag(n, m, z[i - 2]);
                }
            }
            if(verbosity > 0) {
                if(i != it + 1)
                    std::cout << "CG converges after " << i << " iteration" << (i > 1 ? "s" : "") << std::endl;
                else
                    std::cout << "CG does not converges after " << it << " iteration" << (it > 1 ? "s" : "") << std::endl;
            }
            if(std::is_same<ptr_type, K*>::value)
                A.template computeSolution<excluded>(f, x);
            else
                A.template computeSolution<excluded>(storage[1], x);
            delete [] alpha;
            for(auto zCurr : z)
                clean(zCurr);
            for(auto pCurr : p)
                clean(pCurr);
            clean(storage[0]);
            A.clearBuffer(alloc);
            return std::min(i, it);
        }
};
} // HPDDM
#endif // _HPDDM_ITERATIVE_
