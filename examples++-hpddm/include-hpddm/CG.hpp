 /*
   This file is part of HPDDM.

   Author(s): Pierre Jolivet <pierre.jolivet@enseeiht.fr>
        Date: 2015-12-21

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

#ifndef _HPDDM_CG_
#define _HPDDM_CG_

#include "iterative.hpp"

namespace HPDDM {
template<bool excluded, class Operator, class K>
inline int IterativeMethod::CG(const Operator& A, const K* const b, K* const x, const MPI_Comm& comm) {
    const Option& opt = *Option::get();
    if(opt.any_of("schwarz_method", { 0, 1, 4 }) || opt.any_of("schwarz_coarse_correction", { 0 }))
        return GMRES(A, b, x, 1, comm);
    const int n = A.getDof();
    const unsigned short it = opt["max_it"];
    underlying_type<K> tol = opt["tol"];
    const char verbosity = opt.val<char>("verbosity");
    std::cout << std::scientific;
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
                std::cout << "CG: " << std::setw(3) << i << " " << dir[0] << " " << resInit << " " << dir[0] / resInit << " < " << tol << std::endl;
            else
                std::cout << "CG: " << std::setw(3) << i << " " << dir[0] << " < " << -tol << std::endl;
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
    std::cout.unsetf(std::ios_base::scientific);
    return std::min(i, it);
}
template<bool excluded, class Operator, class K>
inline int IterativeMethod::PCG(const Operator& A, const K* const f, K* const x, const MPI_Comm& comm) {
    typedef typename std::conditional<std::is_pointer<typename std::remove_reference<decltype(*A.getScaling())>::type>::value, K**, K*>::type ptr_type;
    const Option& opt = *Option::get();
    const int n = std::is_same<ptr_type, K*>::value ? A.getDof() : A.getMult();
    const unsigned short it = opt["max_it"];
    underlying_type<K> tol = opt["tol"];
    const char verbosity = opt.val<char>("verbosity");
    std::cout << std::scientific;
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
            std::cout << "PCG: " << std::setw(3) << i << " " << resRel << " " << resInit << " " << resRel / resInit << " < " << tol << std::endl;
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
            std::cout << "PCG converges after " << i << " iteration" << (i > 1 ? "s" : "") << std::endl;
        else
            std::cout << "PCG does not converges after " << it << " iteration" << (it > 1 ? "s" : "") << std::endl;
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
    std::cout.unsetf(std::ios_base::scientific);
    return std::min(i, it);
}
} // HPDDM
#endif // _HPDDM_CG_
