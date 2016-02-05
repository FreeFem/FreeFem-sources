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

#ifndef _HPDDM_GMRES_
#define _HPDDM_GMRES_

#include "iterative.hpp"

namespace HPDDM {
template<bool excluded, class Operator, class K>
inline int IterativeMethod::GMRES(const Operator& A, const K* const b, K* const x, const int& mu, const MPI_Comm& comm) {
    const Option& opt = *Option::get();
    const int n = excluded ? 0 : A.getDof();
    const unsigned short it = opt["max_it"];
    underlying_type<K> tol = opt["tol"];
    const char verbosity = opt.val<char>("verbosity");
    std::cout << std::scientific;
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
    if(variant == 'L') {
        A.template apply<excluded>(b, *v, mu, Ax);
        for(unsigned short nu = 0; nu < mu; ++nu)
            norm[nu] = std::real(Blas<K>::dot(&n, *v + nu * n, &i__1, *v + nu * n, &i__1));
    }
    else
        for(unsigned short nu = 0; nu < mu; ++nu) {
            norm[nu] = 0.0;
            for(unsigned int i = 0; i < n; ++i)
                if(std::abs(b[nu * n + i]) <= HPDDM_PEN * HPDDM_EPS)
                    norm[nu] += std::norm(b[nu * n + i]);
        }

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
        unsigned short i = 0;
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
            Arnoldi<excluded>(A, static_cast<char>(opt["gs"]), m, H, v, s, sn, n, i++, mu, Ax, comm);
            for(unsigned short nu = 0; nu < mu; ++nu) {
                if(hasConverged[nu] == -m && ((tol > 0 && std::abs(s[i * mu + nu]) / norm[nu] <= tol) || (tol < 0 && std::abs(s[i * mu + nu]) <= -tol)))
                    hasConverged[nu] = i;
            }
            if(verbosity > 0) {
                int tmp[2] { 0, 0 };
                underlying_type<K> beta = std::abs(s[i * mu]);
                for(unsigned short nu = 0; nu < mu; ++nu) {
                    if(hasConverged[nu] != -m)
                        ++tmp[0];
                    else if(std::abs(s[i * mu + nu]) > beta) {
                        beta = std::abs(s[i * mu + nu]);
                        tmp[1] = nu;
                    }
                }
                if(tol > 0)
                    std::cout << "GMRES: " << std::setw(3) << j << " " << beta << " " <<  norm[tmp[1]] << " " <<  beta / norm[tmp[1]] << " < " << tol;
                else
                    std::cout << "GMRES: " << std::setw(3) << j << " " << beta << " < " << -tol;
                if(mu > 1) {
                    std::cout << " (rhs #" << tmp[1] + 1;
                    if(tmp[0] > 0)
                        std::cout << ", " << tmp[0] << " converged rhs";
                    std::cout << ")";
                }
                std::cout << std::endl;
            }
            if(std::find(hasConverged, hasConverged + mu, -m) == hasConverged + mu) {
                i = 0;
                break;
            }
            else
                ++j;
        }
        if(j != it + 1 && i == m) {
            if(!excluded) {
                if(mu > 1) {
                    for(i = 0; i < m; ++i)
                        Wrapper<K>::template imatcopy<'T'>(i + 1, mu, H[i], mu, m + 1);
                }
                updateSol(A, variant, n, x, H, s, v + (m + 1) * (variant == 'F'), hasConverged, mu, Ax);
            }
            if(verbosity > 0)
                std::cout << "GMRES restart(" << m << ")" << std::endl;
        }
        else
            break;
    }
    if(!excluded) {
        const int rem = it % m;
        std::for_each(hasConverged, hasConverged + mu, [&](short& d) { if(d < 0) d = rem > 0 ? rem : -d; });
        if(mu > 1) {
            unsigned short dim = *std::max_element(hasConverged, hasConverged + mu);
            for(unsigned short i = 0; i < dim; ++i)
                Wrapper<K>::template imatcopy<'T'>(i + 1, mu, H[i], mu, m + 1);
        }
        updateSol(A, variant, n, x, H, s, v + (m + 1) * (variant == 'F'), hasConverged, mu, Ax);
    }
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
    std::cout.unsetf(std::ios_base::scientific);
    return std::min(j, it);
}
template<bool excluded, class Operator, class K>
inline int IterativeMethod::BGMRES(const Operator& A, const K* const b, K* const x, const int& mu, const MPI_Comm& comm) {
    const Option& opt = *Option::get();
    const int n = excluded ? 0 : A.getDof();
    const unsigned short it = opt["max_it"];
    underlying_type<K> tol = opt["tol"];
    const char verbosity = opt.val<char>("verbosity");
    std::cout << std::scientific;
    if(std::abs(tol) < std::numeric_limits<underlying_type<K>>::epsilon()) {
        if(verbosity > 0)
            std::cout << "WARNING -- the tolerance of the iterative method was set to " << tol << " which is lower than the machine epsilon for type " << demangle(typeid(underlying_type<K>).name()) << ", forcing the tolerance to " << 2 * std::numeric_limits<underlying_type<K>>::epsilon() << std::endl;
        tol = 2 * std::numeric_limits<underlying_type<K>>::epsilon();
    }
    const unsigned short m = std::min(static_cast<unsigned short>(std::numeric_limits<short>::max()), std::min(static_cast<unsigned short>(opt["gmres_restart"]), it));
    const char variant = (opt["variant"] == 0 ? 'L' : opt["variant"] == 1 ? 'R' : 'F');

    K** const H = new K*[m * (2 + (variant == 'F')) + 1];
    K** const v = H + m;
    int ldh = mu * (m + 1);
    int info;
    int N = 2 * mu;
    int lwork = mu * std::max(n, opt["gs"] != 1 ? ldh : mu);
    *H = new K[lwork + mu * ((m + 1) * ldh + n * (m * (1 + (variant == 'F')) + 1) + 2 * m) + (Wrapper<K>::is_complex ? (mu + 1) / 2 : mu)];
    *v = *H + m * mu * ldh;
    K* const s = *v + mu * n * (m * (1 + (variant == 'F')) + 1);
    K* const tau = s + mu * ldh;
    K* const Ax = tau + m * N;
    underlying_type<K>* const norm = reinterpret_cast<underlying_type<K>*>(Ax + lwork);
    underlying_type<K>* const beta = norm - mu;
    bool alloc = A.setBuffer(mu);

    A.template start<excluded>(b, x, mu);
    if(variant == 'L') {
        A.template apply<excluded>(b, *v, mu, Ax);
        for(unsigned short nu = 0; nu < mu; ++nu)
            norm[nu] = std::real(Blas<K>::dot(&n, *v + nu * n, &i__1, *v + nu * n, &i__1));
    }
    else
        for(unsigned short nu = 0; nu < mu; ++nu) {
            norm[nu] = 0.0;
            for(unsigned int i = 0; i < n; ++i)
                if(std::abs(b[nu * n + i]) <= HPDDM_PEN * HPDDM_EPS)
                    norm[nu] += std::norm(b[nu * n + i]);
        }

    if(!excluded)
        for(unsigned short nu = 0; nu < mu; ++nu)
            for(unsigned int i = 0; i < n; ++i)
                if(std::abs(b[nu * n + i]) > HPDDM_PEN * HPDDM_EPS)
                    depenalize(b[nu * n + i], x[nu * n + i]);

    unsigned short j = 1;
    short dim = mu * m;
    int* const piv = new int[mu];
    underlying_type<K>* workpiv = norm - 2 * mu;
    int deflated = -1;
    while(j <= it) {
        if(!excluded)
            A.GMV(x, variant == 'L' ? Ax : *v, mu);
        Blas<K>::axpby(mu * n, 1.0, b, 1, -1.0, variant == 'L' ? Ax : *v, 1);
        if(variant == 'L')
            A.template apply<excluded>(Ax, *v, mu);
        if(j == 1) {
            for(unsigned short nu = 0; nu < mu; ++nu)
                beta[nu] = std::real(Blas<K>::dot(&n, *v + nu * n, &i__1, *v + nu * n, &i__1));
            MPI_Allreduce(MPI_IN_PLACE, beta, 2 * mu, Wrapper<K>::mpi_underlying_type(), MPI_SUM, comm);
            for(unsigned short nu = 0; nu < mu; ++nu) {
                norm[nu] = std::sqrt(norm[nu]);
                if(norm[nu] < HPDDM_EPS)
                    norm[nu] = 1.0;
                beta[nu] = std::sqrt(beta[nu]);
                if(tol > 0.0 && beta[nu] / norm[nu] < tol && norm[nu] > 1.0 / HPDDM_EPS)
                    norm[nu] = 1.0;
            }
        }
        Blas<K>::herk("U", "C", &mu, &n, &(Wrapper<underlying_type<K>>::d__1), *v, &n, &(Wrapper<underlying_type<K>>::d__0), Ax, &mu);
        for(unsigned short nu = 1; nu < mu; ++nu)
            std::copy_n(Ax + nu * mu, nu + 1, Ax + (nu * (nu + 1)) / 2);
        MPI_Allreduce(MPI_IN_PLACE, Ax, (mu * (mu + 1)) / 2, Wrapper<K>::mpi_type(), MPI_SUM, comm);
        for(unsigned short nu = mu; nu-- > 0; )
            std::copy_n(Ax + (nu * (nu + 1)) / 2, nu + 1, s + nu * mu);
        if(!opt.set("initial_deflation_tol")) {
            Lapack<K>::potrf("U", &mu, s, &mu, &info);
            if(verbosity > 3) {
                std::cout << "BGMRES diag(R), QR = block residual: ";
                std::cout << s[0];
                if(mu > 1) {
                    if(mu > 2)
                        std::cout << "\t...";
                    std::cout << "\t" << s[(mu - 1) * (mu + 1)];
                }
                std::cout << std::endl;
            }
            N = (info > 0 ? info - 1 : mu);
        }
        else {
            Lapack<K>::pstrf("U", &mu, s, &mu, piv, &N, &(Wrapper<underlying_type<K>>::d__0), workpiv, &info);
            if(verbosity > 3) {
                std::cout << "BGMRES diag(R), QR = block residual, with pivoting: ";
                std::cout << s[0] << " (" << piv[0] << ")";
                if(mu > 1) {
                    if(mu > 2)
                        std::cout << "\t...";
                    std::cout << "\t" << s[(mu - 1) * (mu + 1)] << " (" << piv[mu - 1] << ")";
                }
                std::cout << std::endl;
            }
            if(info == 0) {
                N = mu;
                while(N > 1 && std::abs(s[(N - 1) * (mu + 1)] / s[0]) <= opt.val("initial_deflation_tol"))
                    --N;
            }
            Lapack<K>::lapmt(&i__1, &n, &mu, *v, &n, piv);
            Lapack<underlying_type<K>>::lapmt(&i__1, &i__1, &mu, norm, &i__1, piv);
        }
        if(N != mu) {
            int nrhs = mu - N;
            Lapack<K>::trtrs("U", "N", "N", &N, &nrhs, s, &mu, s + N * mu, &mu, &info);
        }
        if(N != deflated) {
            deflated = N;
            dim = deflated * (j - 1 + m > it ? it - j + 1 : m);
            ldh = deflated * (m + 1);
            for(unsigned short i = 1; i < m; ++i)
                H[i] = *H + i * deflated * ldh;
            for(unsigned short i = 1; i < m * (1 + (variant == 'F')) + 1; ++i)
                v[i] = *v + i * deflated * n;
        }
        N *= 2;
        std::fill_n(tau, m * N, K());
        Wrapper<K>::template imatcopy<'N'>(mu, mu, s, mu, ldh);
        Blas<K>::trsm("R", "U", "N", "N", &n, &deflated, &(Wrapper<K>::d__1), s, &ldh, *v, &n);
        for(unsigned short nu = 0; nu < deflated; ++nu)
            std::fill(s + nu * (ldh + 1) + 1, s + (nu + 1) * ldh, K());
        std::fill(*H, *v, K());
        unsigned short i = 0;
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
            if(BlockArnoldi<excluded>(A, static_cast<char>(opt["gs"]), m, H, v, tau, s, lwork, n, i++, deflated, Ax, comm)) {
                dim = deflated * (i - 1);
                i = m;
                j = it + 1;
                break;
            }
            unsigned short converged = 0;
            for(unsigned short nu = 0; nu < deflated; ++nu) {
                beta[nu] = Blas<K>::nrm2(&deflated, s + deflated * i + nu * ldh, &i__1);
                if(((tol > 0 && beta[nu] / norm[nu] <= tol) || (tol < 0 && beta[nu] <= -tol)))
                    ++converged;
            }
            if(verbosity > 0) {
                underlying_type<K>* max = std::max_element(beta, beta + deflated);
                if(tol > 0)
                    std::cout << "BGMRES: " << std::setw(3) << j << " " << *max << " " <<  norm[std::distance(beta, max)] << " " <<  *max / norm[std::distance(beta, max)] << " < " << tol;
                else
                    std::cout << "BGMRES: " << std::setw(3) << j << " " << *max << " < " << -tol;
                std::cout << " (rhs #" << std::distance(beta, max) + 1;
                if(converged > 0)
                    std::cout << ", " << converged << " converged rhs";
                if(deflated != mu)
                    std::cout << ", " << mu - deflated << " deflated rhs";
                std::cout << ")" << std::endl;
            }
            if(converged == deflated) {
                dim = deflated * i;
                i = 0;
                break;
            }
            else
                ++j;
        }
        if(j != it + 1 && i == m) {
            if(opt.set("initial_deflation_tol"))
                Lapack<K>::lapmt(&i__1, &n, &mu, x, &n, piv);
            if(!excluded)
                updateSol(A, variant, n, x, H, s, v + (m + 1) * (variant == 'F'), &dim, mu, Ax, deflated);
            if(opt.set("initial_deflation_tol")) {
                Lapack<K>::lapmt(&i__0, &n, &mu, x, &n, piv);
                Lapack<underlying_type<K>>::lapmt(&i__0, &i__1, &mu, norm, &i__1, piv);
            }
            if(verbosity > 0)
                std::cout << "BGMRES restart(" << m << ")" << std::endl;
        }
        else
            break;
    }
    if(opt.set("initial_deflation_tol"))
        Lapack<K>::lapmt(&i__1, &n, &mu, x, &n, piv);
    if(!excluded) {
        const int rem = it % m;
        if(rem != 0)
            dim = deflated * rem;
        updateSol(A, variant, n, x, H, s, v + (m + 1) * (variant == 'F'), &dim, mu, Ax, deflated);
    }
    if(opt.set("initial_deflation_tol"))
        Lapack<K>::lapmt(&i__0, &n, &mu, x, &n, piv);
    delete [] piv;

    if(verbosity > 0) {
        if(j != it + 1)
            std::cout << "BGMRES converges after " << j << " iteration" << (j > 1 ? "s" : "") << std::endl;
        else
            std::cout << "BGMRES does not converges after " << it << " iteration" << (it > 1 ? "s" : "") << std::endl;
    }
    A.clearBuffer(alloc);
    delete [] *H;
    delete [] H;
    std::cout.unsetf(std::ios_base::scientific);
    return std::min(j, it);
}
} // HPDDM
#endif // _HPDDM_GMRES_
