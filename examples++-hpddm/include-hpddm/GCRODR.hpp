 /*
   This file is part of HPDDM.

   Author(s): Pierre Jolivet <pierre.jolivet@enseeiht.fr>
        Date: 2016-01-06

   Copyright (C) 2016-     Centre National de la Recherche Scientifique

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

#ifndef _HPDDM_GCRODR_
#define _HPDDM_GCRODR_

#include "GMRES.hpp"

namespace HPDDM {
template<class K>
class Recycling : private Singleton {
    private:
        K*              _storage;
        unsigned short       _mu;
        unsigned short        _k;
    public:
        bool               _same;
        template<int N>
        Recycling(Singleton::construct_key<N>, unsigned short mu) : _storage(), _mu(mu), _same(false) { }
        ~Recycling() {
            delete [] _storage;
        }
        bool recycling() const {
            return _storage != nullptr;
        }
        void allocate(int n, unsigned short k) {
            _k = k;
            _storage = new K[2 * _mu * _k * n];
        }
        K* storage() const {
            return _storage;
        }
        unsigned short k() const {
            return _k;
        }
        template<int N = 0>
        static std::shared_ptr<Recycling> get(unsigned short mu) {
            return Singleton::get<Recycling, N>(mu);
        }
};

template<bool excluded, class Operator, class K>
inline int IterativeMethod::GCRODR(const Operator& A, const K* const b, K* const x, const int& mu, const MPI_Comm& comm) {
    const Option& opt = *Option::get();
    int k = opt.val<int>("gmres_recycle", 0);
    if(k == 0)
        return GMRES(A, b, x, mu, comm);
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
    const int m = std::min(static_cast<unsigned short>(std::numeric_limits<short>::max()), std::min(static_cast<unsigned short>(opt["gmres_restart"]), it));
    k = std::min(m - 1, k);
    const char variant = (opt["variant"] == 0 ? 'L' : opt["variant"] == 1 ? 'R' : 'F');

    const int ldh = mu * (m + 1);
    K** const H = new K*[m * (3 + (variant == 'F')) + 1];
    K** const save = H + m;
    *save = new K[ldh * m];
    K** const v = save + m;
    K* const s = new K[mu * ((m + 1) * (m + 1) + n * (2 + (variant == 'R') + m * (1 + (variant == 'F'))) + (!Wrapper<K>::is_complex ? m + 1 : (m + 2) / 2))];
    K* const Ax = s + ldh;
    const int ldv = mu * n;
    *H = Ax + (1 + (variant == 'R')) * ldv;
    for(unsigned short i = 1; i < m; ++i) {
        H[i] = *H + i * ldh;
        save[i] = *save + i * ldh;
    }
    *v = *H + m * ldh;
    for(unsigned short i = 1; i < m * (1 + (variant == 'F')) + 1; ++i)
        v[i] = *v + i * ldv;
    underlying_type<K>* const norm = reinterpret_cast<underlying_type<K>*>(*v + (m * (1 + (variant == 'F')) + 1) * ldv);
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

    int info;
    unsigned short j = 1;
    bool recycling;
    K* U, *C;
    Recycling<K>& recycled = *Recycling<K>::get(mu);
    if(recycled.recycling()) {
        recycling = true;
        k = recycled.k();
        U = recycled.storage();
        C = U + k * ldv;
    }
    else
        recycling = false;
    while(j <= it) {
        unsigned short i = (recycling ? k : 0);
        if(!excluded)
            A.GMV(x, variant == 'L' ? Ax : v[i], mu);
        Blas<K>::axpby(ldv, 1.0, b, 1, -1.0, variant == 'L' ? Ax : v[i], 1);
        if(variant == 'L')
            A.template apply<excluded>(Ax, v[i], mu);
        if(j == 1 && recycling) {
            if(!recycled._same) {
                for(unsigned short nu = 0; nu < k; ++nu) {
                    if(variant == 'L') {
                        A.GMV(U + nu * ldv, Ax, mu);
                        A.template apply<excluded>(Ax, C + nu * ldv, mu);
                    }
                    else {
                        A.template apply<excluded>(U + nu * ldv, variant == 'F' ? v[m + 1] + nu * n : Ax, mu, v[i + 1]);
                        A.GMV(variant == 'F' ? v[m + 1] + nu * n : Ax, C + nu * ldv, mu);
                    }
                }
                K* work = new K[k * k * mu];
                for(unsigned short nu = 0; nu < mu; ++nu) {
                    Blas<K>::herk("U", "C", &k, &n, &(Wrapper<underlying_type<K>>::d__1), C + nu * n, &ldv, &(Wrapper<underlying_type<K>>::d__0), work + nu * (k * (k + 1)) / 2, &k);
                    for(unsigned short xi = 1; xi < k; ++xi)
                        std::copy_n(work + nu * (k * (k + 1)) / 2 + xi * k, xi + 1, work + nu * (k * (k + 1)) / 2 + (xi * (xi + 1)) / 2);
                }
                MPI_Allreduce(MPI_IN_PLACE, work, mu * (k * (k + 1)) / 2, Wrapper<K>::mpi_type(), MPI_SUM, comm);
                for(unsigned short nu = mu; nu-- > 0; ) {
                    for(unsigned short xi = k; xi > 0; --xi)
                        std::copy_backward(work + nu * (k * (k + 1)) / 2 + (xi * (xi - 1)) / 2, work + nu * (k * (k + 1)) / 2 + (xi * (xi + 1)) / 2, work + k * k * nu + k * xi - (k - xi));
                    Lapack<K>::potrf("U", &k, work + k * k * nu, &k, &info);
                    Blas<K>::trsm("R", "U", "N", "N", &n, &k, &(Wrapper<K>::d__1), work + k * k * nu, &k, U + nu * n, &ldv);
                }
                delete [] work;
            }
            orthonormalization<excluded>(0, n, k, mu, C, v[i], H[i], comm);
            if(variant == 'L')
                for(unsigned short nu = 0; nu < mu; ++nu)
                    Blas<K>::gemv("N", &n, &k, &(Wrapper<K>::d__1), U + nu * n, &ldv, H[i] + nu, &mu, &(Wrapper<K>::d__1), x + nu * n, &i__1);
            else {
                for(unsigned short nu = 0; nu < k; ++nu)
                    A.template apply<excluded>(U + nu * ldv, v[(m + 1) * (variant == 'F')] + nu * ldv, mu, Ax);
                for(unsigned short nu = 0; nu < mu; ++nu)
                    Blas<K>::gemv("N", &n, &k, &(Wrapper<K>::d__1), v[(m + 1) * (variant == 'F')] + nu * n, &ldv, H[i] + nu, &mu, &(Wrapper<K>::d__1), x + nu * n, &i__1);
            }
            std::copy_n(C, k * ldv, *v);
        }
        for(unsigned short nu = 0; nu < mu; ++nu)
            sn[i * mu + nu] = std::real(Blas<K>::dot(&n, v[i] + nu * n, &i__1, v[i] + nu * n, &i__1));
        if(j == 1) {
            MPI_Allreduce(MPI_IN_PLACE, norm, 2 * mu, Wrapper<K>::mpi_underlying_type(), MPI_SUM, comm);
            for(unsigned short nu = 0; nu < mu; ++nu) {
                norm[nu] = std::sqrt(norm[nu]);
                if(norm[nu] < HPDDM_EPS)
                    norm[nu] = 1.0;
                if(tol > 0.0 && sn[i * mu + nu] / norm[nu] < tol) {
                    if(norm[nu] > 1.0 / HPDDM_EPS)
                        norm[nu] = 1.0;
                    else
                        hasConverged[nu] = 0;
                }
                else if(sn[i * mu + nu] < -tol)
                    hasConverged[nu] = 0;
            }
        }
        else
            MPI_Allreduce(MPI_IN_PLACE, sn + i * mu, mu, Wrapper<K>::mpi_underlying_type(), MPI_SUM, comm);
        if(recycling) {
            *v += k * ldv;
            if(variant == 'F')
                v[m + 1] += k * ldv;
        }
        for(unsigned short nu = 0; nu < mu; ++nu) {
            if(hasConverged[nu] > 0)
                hasConverged[nu] = 0;
            s[mu * i + nu] = std::sqrt(sn[i * mu + nu]);
            if(recycling)
                sn[nu] = std::real(s[mu * i + nu]);
            std::for_each(*v + nu * n, *v + (nu + 1) * n, [&](K& y) { y /= s[mu * i + nu]; });
        }
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
            if(recycling) {
                const char gs = static_cast<char>(opt["gs"]);
                orthonormalization<excluded>(gs, n, k, mu, C, v[i + 1], H[i], comm);
                H[i - k] += k * (mu + ldh);
                v[i - k + 1] += k * ldv;
                if(variant == 'F')
                    v[i - k + m + 2] += k * ldv;
                Arnoldi<excluded>(A, gs, m - k, H, v, s + k * mu, sn + k * mu, n, i++ - k, mu, Ax, comm, save);
            }
            else
                Arnoldi<excluded>(A, static_cast<char>(opt["gs"]), m, H, v, s, sn, n, i++, mu, Ax, comm, save);
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
                    std::cout << "GCRODR: " << std::setw(3) << j << " " << beta << " " <<  norm[tmp[1]] << " " <<  beta / norm[tmp[1]] << " < " << tol;
                else
                    std::cout << "GCRODR: " << std::setw(3) << j << " " << beta << " < " << -tol;
                if(mu > 1) {
                    std::cout << " (rhs #" << tmp[1] + 1;
                    if(tmp[0] > 0)
                        std::cout << ", " << tmp[0] << " converged rhs";
                    std::cout << ")";
                }
                std::cout << std::endl;
            }
            if(std::find(hasConverged, hasConverged + mu, -m) == hasConverged + mu) {
                i += m - k;
                break;
            }
            else
                ++j;
        }
        if(recycling) {
            for(unsigned short j = 0, shift = (i <= m ? i - k : i - m); j < shift; ++j) {
                if(variant == 'F')
                    v[j + m + 2] -= k * ldv;
                v[j + 1] -= k * ldv;
                H[j] -= k * (mu + ldh);
            }
            if(variant == 'F')
                v[m + 1] -= k * ldv;
            *v -= k * ldv;
        }
        if(j != it + 1 && i == m) {
            if(!excluded) {
                if(mu > 1) {
                    for(i = (recycling ? k : 0); i < m; ++i) {
                        Wrapper<K>::template imatcopy<'T'>(i + 2, mu, H[i], mu, m + 1);
                        for(unsigned short nu = 0; nu < mu; ++nu)
                            std::fill_n(H[i] + i + 2 + nu * (m + 1), m - i - 1, K());
                    }
                }
                if(!recycling)
                    updateSol(A, variant, n, x, H, s, static_cast<const K* const* const>(v + (m + 1) * (variant == 'F')), hasConverged, mu, Ax);
                else {
                    computeMin(H, s + k * mu, hasConverged, mu, -1, k);
                    int inc = mu;
                    int diff = m - k;
                    for(unsigned short nu = 0; nu < mu; ++nu) {
                        **v = sn[nu];
                        Blas<K>::gemv(&(Wrapper<K>::transc), &n, &k, *v, C + nu * n, &ldv, *v + (k * mu + nu) * n , &i__1, &(Wrapper<K>::d__0), s + nu, &inc);
                    }
                    MPI_Allreduce(MPI_IN_PLACE, s, k * mu, Wrapper<K>::mpi_type(), MPI_SUM, comm);
                    for(unsigned short nu = 0; nu < mu; ++nu) {
                        Blas<K>::gemv("N", &k, &diff, &(Wrapper<K>::d__2), H[k] + nu * (m + 1), &ldh, s + k * mu + nu, &inc, &(Wrapper<K>::d__1), s + nu, &inc);
                    }
                    std::copy_n(U, k * ldv, v[(m + 1) * (variant == 'F')]);
                    addSol(A, variant, n, x, std::distance(H[0], H[1]), s, static_cast<const K* const* const>(v + (m + 1) * (variant == 'F')), hasConverged, mu, Ax);
                }
            }
            else {
                if(!recycling)
                    updateSol(A, variant, n, x, H, s, v + (m + 1) * (variant == 'F'), hasConverged, mu, Ax);
                else
                    addSol(A, variant, n, x, std::distance(H[0], H[1]), s, v + (m + 1) * (variant == 'F'), hasConverged, mu, Ax);
            }
            if(verbosity > 0)
                std::cout << "GCRODR restart(" << m << ", " << k << ")" << std::endl;
            if(mu > 1) {
                for(i = 0; i < m - (recycling ? k : 0); ++i)
                    Wrapper<K>::template imatcopy<'T'>(i + 2, mu, save[i], mu, m + 1);
            }
            if(!recycling) {
                recycling = true;
                recycled.allocate(n, k);
                U = recycled.storage();
                C = U + k * ldv;
                std::fill_n(s, m * mu, K());
                for(unsigned short nu = 0; nu < mu; ++nu) {
                    if(hasConverged[nu] == -m) {
                        K h = H[m - 1][(m + 1) * nu + m] / H[m - 1][(m + 1) * nu + m - 1];
                        for(i = m; i-- > 1; ) {
                            s[i + m * nu] = H[i - 1][(m + 1) * nu + i] * h;
                            h *= -sn[(i - 1) * mu + nu];
                        }
                        s[m * nu] = h;
                        for(i = 0; i < m; ++i) {
                            std::fill_n(save[i] + i + 2 + nu * (m + 1), m - i - 1, K());
                            std::copy_n(save[i] + nu * (m + 1), m + 1, H[i] + nu * (m + 1));
                        }
                        h = save[m - 1][m + nu * (m + 1)] * save[m - 1][m + nu * (m + 1)];
                        Blas<K>::axpy(&m, &h, s + m * nu, &i__1, H[m - 1] + nu * (m + 1), &i__1);
                        int* select = new int[m]();
                        int row = m + 1;
                        int lwork = -1;
                        Lapack<K>::hseqr("E", "N", &m, &i__1, &m, nullptr, &ldh, nullptr, nullptr, nullptr, &i__1, &h, &lwork, &info);
                        lwork = std::max(static_cast<int>(std::real(h)), Wrapper<K>::is_complex ? m * m : (m * (m + 2)));
                        *select = -1;
                        Lapack<K>::geqrf(&row, &k, nullptr, &ldh, nullptr, &h, select, &info);
                        lwork = std::max(static_cast<int>(std::real(h)), lwork);
                        Lapack<K>::mqr("R", "N", &n, &row, &k, nullptr, &ldh, nullptr, nullptr, &ldv, &h, select, &info);
                        *select = 0;
                        lwork = std::max(static_cast<int>(std::real(h)), lwork);
                        K* work = new K[lwork];
                        K* w = new K[Wrapper<K>::is_complex ? m : (2 * m)];
                        K* backup = new K[m * m];
                        for(i = 0; i < m; ++i)
                            std::copy_n(*H + nu * (m + 1) + i * ldh, m, backup + i * m);
                        Lapack<K>::hseqr("E", "N", &m, &i__1, &m, backup, &m, w, w + m, nullptr, &i__1, work, &lwork, &info);
                        delete [] backup;
                        std::vector<std::pair<unsigned short, underlying_type<K>>> p;
                        p.reserve(k + 1);
                        for(i = 0; i < m; ++i) {
                            underlying_type<K> magnitude = Wrapper<K>::is_complex ? std::norm(w[i]) : std::real(w[i] * w[i] + w[m + i] * w[m + i]);
                            typename decltype(p)::const_iterator it = std::lower_bound(p.cbegin(), p.cend(), std::make_pair(i, magnitude), [](const std::pair<unsigned short, underlying_type<K>>& lhs, const std::pair<unsigned short, underlying_type<K>>& rhs) { return lhs.second < rhs.second; });
                            if(p.size() < k || it != p.cend())
                                p.insert(it, std::make_pair(i, magnitude));
                            if(p.size() == k + 1)
                                p.pop_back();
                        }
                        int mm = Wrapper<K>::is_complex ? k : 0;
                        for(typename decltype(p)::const_iterator it = p.cbegin(); it != p.cend(); ++it) {
                            if(Wrapper<K>::is_complex)
                                select[it->first] = 1;
                            else {
                                if(std::abs(w[m + it->first]) < HPDDM_EPS) {
                                    select[it->first] = 1;
                                    ++mm;
                                }
                                else if(mm < k + 1) {
                                    select[it->first] = 1;
                                    mm += 2;
                                    ++it;
                                }
                                else
                                    break;
                            }
                        }
                        underlying_type<K>* rwork = Wrapper<K>::is_complex ? new underlying_type<K>[k] : nullptr;
                        K* vr = new K[mm * m];
                        int* ifailr = new int[mm];
                        int col;
                        Lapack<K>::hsein("R", "Q", "N", select, &m, *H + nu * (m + 1), &ldh, w, w + m, nullptr, &i__1, vr, &m, &mm, &col, work, rwork, nullptr, ifailr, &info);
                        delete [] ifailr;
                        delete [] select;
                        delete [] rwork;
                        delete [] w;
                        Blas<K>::gemm("N", "N", &n, &k, &m, &(Wrapper<K>::d__1), v[(m + 1) * (variant == 'F')] + nu * n, &ldv, vr, &m, &(Wrapper<K>::d__0), U + nu * n, &ldv);
                        Blas<K>::gemm("N", "N", &row, &k, &m, &(Wrapper<K>::d__1), *save, &ldh, vr, &m, &(Wrapper<K>::d__0), *H + nu * (m + 1), &ldh);
                        delete [] vr;
                        K* tau = new K[k];
                        std::for_each(v[m] + nu * n, v[m] + (nu + 1) * n, [&](K& y) { y /= save[m - 1][m + nu * (m + 1)]; });
                        Lapack<K>::geqrf(&row, &k, *H + nu * (m + 1), &ldh, tau, work, &lwork, &info);
                        Lapack<K>::mqr("R", "N", &n, &row, &k, *H + nu * (m + 1), &ldh, tau, *v + nu * n, &ldv, work, &lwork, &info);
                        Wrapper<K>::template omatcopy<'N'>(k, n, *v + nu * n, ldv, C + nu * n, ldv);
                        Blas<K>::trsm("R", "U", "N", "N", &n, &k, &(Wrapper<K>::d__1), *H + nu * (m + 1), &ldh, U + nu * n, &ldv);
                        delete [] tau;
                        delete [] work;
                    }
                }
            }
            else if(!recycled._same) {
                std::copy_n(C, k * ldv, *v);
                K* prod = new K[k * mu * (m + 2)];
                if(variant == 'F')
                    std::copy_n(v[m + 1], k * ldv, U);
                for(unsigned short nu = 0; nu < mu; ++nu) {
                    std::for_each(v[m] + nu * n, v[m] + (nu + 1) * n, [&](K& y) { y /= save[m - k - 1][m - k + nu * (m + 1)]; });
                    int row = m + 1;
                    Blas<K>::gemm(&(Wrapper<K>::transc), "N", &row, &k, &n, &(Wrapper<K>::d__1), *v + nu * n, &ldv, U + nu * n, &ldv, &(Wrapper<K>::d__0), prod + k * nu * row, &row);
                    for(i = 0; i < k; ++i)
                        prod[k * mu * row + k * nu + i] = Blas<K>::dot(&n, U + nu * n + i * ldv, &i__1, U + nu * n + i * ldv, &i__1);
                }
                MPI_Allreduce(MPI_IN_PLACE, prod, k * mu * (m + 2), Wrapper<K>::mpi_type(), MPI_SUM, comm);
                std::for_each(prod + k * mu * (m + 1), prod + k * mu * (m + 2), [](K& u) { u = 1.0 / std::sqrt(u); });
                for(unsigned short nu = 0; nu < mu; ++nu) {
                    if(hasConverged[nu] == -m) {
                        for(i = 0; i < k; ++i)
                            Blas<K>::scal(&n, prod + k * mu * (m + 1) + k * nu + i, U + nu * n + i * ldv, &i__1);
                        for(i = 0; i < m; ++i)
                            std::fill_n(save[i] + i + 2 + nu * (m + 1), m - i - 1, K());
                        int lda = m;
                        K* A = new K[m * m];
                        for(i = 0; i < k; ++i)
                            for(unsigned short j = 0; j < k; ++j)
                                A[j + i * m] = (i == j ? prod[k * mu * (m + 1) + k * nu + i] * prod[k * mu * (m + 1) + k * nu + i] : Wrapper<K>::d__0);
                        int diff = m - k;
                        Wrapper<K>::template omatcopy<'N'>(m - k, k, H[k] + nu * (m + 1), ldh, A + k * m, m);
                        for(i = 0; i < k; ++i)
                            Blas<K>::scal(&diff, prod + k * mu * (m + 1) + k * nu + i, A + k * m + i, &lda);
                        Wrapper<K>::template omatcopy<'C'>(m - k, k, A + k * m, m, A + k, m);
                        int row = diff + 1;
                        Blas<K>::gemm(&(Wrapper<K>::transc), "N", &diff, &diff, &k, &(Wrapper<K>::d__1), H[k] + nu * (m + 1), &ldh, H[k] + nu * (m + 1), &ldh, &(Wrapper<K>::d__0), A + k * m + k, &lda);
                        Blas<K>::gemm(&(Wrapper<K>::transc), "N", &diff, &diff, &row, &(Wrapper<K>::d__1), *save + nu * (m + 1), &ldh, *save + nu * (m + 1), &ldh, &(Wrapper<K>::d__1), A + k * m + k, &lda);
                        K* B = new K[m * (m + 1)]();
                        row = m + 1;
                         for(i = 0; i < k; ++i)
                             std::transform(prod + k * nu * (m + 1) + i * (m + 1), prod + k * nu * (m + 1) + (i + 1) * (m + 1), B + i * (m + 1), [&](const K& u) { return prod[k * mu * (m + 1) + k * nu + i] * u; });
                        Wrapper<K>::template omatcopy<'C'>(m - k, m - k, *save + nu * (m + 1), ldh, B + k + k * (m + 1), m + 1);
                        Blas<K>::gemm(&(Wrapper<K>::transc), "N", &diff, &k, &row, &(Wrapper<K>::d__1), *save + nu * (m + 1), &ldh, B + k, &row, &(Wrapper<K>::d__0), *H + k + 1 + nu * (m + 1), &ldh);
                        Wrapper<K>::template omatcopy<'N'>(k, m - k, *H + k + 1 + nu * (m + 1), ldh, B + k, m + 1);
                        Blas<K>::gemm(&(Wrapper<K>::transc), "N", &diff, &k, &k, &(Wrapper<K>::d__1), H[k] + nu * (m + 1), &ldh, B, &row, &(Wrapper<K>::d__1), B + k, &row);
                        for(i = 0; i < k; ++i)
                            Blas<K>::scal(&k, prod + k * mu * (m + 1) + k * nu + i, B + i, &row);
                        K* alpha = new K[(2 + !Wrapper<K>::is_complex) * lda];
                        int lwork = -1;
                        K* vr = new K[m * m];
                        Lapack<K>::ggev("N", "V", &lda, A, &lda, B, &row, alpha, alpha + 2 * lda, alpha + lda, nullptr, &i__1, nullptr, &m, *H + k + 1, &lwork, nullptr, &info);
                        lwork = std::real(H[0][k + 1]);
                        K* work = new K[Wrapper<K>::is_complex ? (lwork + lda / 4 + 2) : lwork];
                        underlying_type<K>* rwork = reinterpret_cast<underlying_type<K>*>(work + lwork);
                        Lapack<K>::ggev("N", "V", &lda, A, &lda, B, &row, alpha, alpha + 2 * lda, alpha + lda, nullptr, &i__1, vr, &m, work, &lwork, rwork, &info);
                        std::vector<std::pair<unsigned short, underlying_type<K>>> q;
                        q.reserve(m);
                        for(i = 0; i < m; ++i) {
                            underlying_type<K> magnitude = Wrapper<K>::is_complex ? std::norm(alpha[i] / alpha[lda + i]) : std::real((alpha[i] * alpha[i] + alpha[2 * lda + i] * alpha[2 * lda + i]) / (alpha[lda + i] * alpha[lda + i]));
                            q.emplace_back(i, magnitude);
                        }
                        std::sort(q.begin(), q.end(), [](const std::pair<unsigned short, underlying_type<K>>& lhs, const std::pair<unsigned short, underlying_type<K>>& rhs) { return lhs.second < rhs.second; });
                        info = std::accumulate(q.cbegin(), q.cbegin() + k, 0, [](int a, const std::pair<unsigned short, underlying_type<K>>& b) { return a + b.first; });
                        for(i = k; info != (k * (k - 1)) / 2 && i < m; ++i)
                            info += q[i].first;
                        int* perm = new int[i];
                        for(unsigned short j = 0; j < i; ++j)
                            perm[j] = q[j].first + 1;
                        Lapack<K>::lapmt(&i__1, &m, &(info = i), vr, &m, perm);
                        row = diff + 1;
                        Blas<K>::gemm("N", "N", &row, &k, &diff, &(Wrapper<K>::d__1), *save + nu * (m + 1), &ldh, vr + k, &m, &(Wrapper<K>::d__0), *H + k + nu * (m + 1), &ldh);
                        Wrapper<K>::template omatcopy<'N'>(k, k, vr, m, *H + nu * (m + 1), ldh);
                        for(i = 0; i < k; ++i)
                            Blas<K>::scal(&k, prod + k * mu * (m + 1) + k * nu + i, *H + nu * (m + 1) + i, &ldh);
                        Blas<K>::gemm("N", "N", &k, &k, &diff, &(Wrapper<K>::d__1), H[k] + nu * (m + 1), &ldh, vr + k, &m, &(Wrapper<K>::d__1), *H + nu * (m + 1), &ldh);
                        row = m + 1;
                        K* tau = new K[k];
                        *perm = -1;
                        Lapack<K>::geqrf(&row, &k, nullptr, &ldh, nullptr, work, perm, &info);
                        Lapack<K>::mqr("R", "N", &n, &row, &k, nullptr, &ldh, nullptr, nullptr, &ldv, work + 1, perm, &info);
                        delete [] perm;
                        if(std::real(work[0]) > lwork || std::real(work[1]) > lwork) {
                            lwork = std::max(std::real(work[0]), std::real(work[1]));
                            delete [] work;
                            work = new K[lwork];
                        }
                        Lapack<K>::geqrf(&row, &k, *H + nu * (m + 1), &ldh, tau, work, &lwork, &info);
                        Wrapper<K>::template omatcopy<'N'>(k, n, U + nu * n, ldv, v[(m + 1) * (variant == 'F')] + nu * n, ldv);
                        Blas<K>::gemm("N", "N", &n, &k, &m, &(Wrapper<K>::d__1), v[(m + 1) * (variant == 'F')] + nu * n, &ldv, vr, &m, &(Wrapper<K>::d__0), U + nu * n, &ldv);
                        Blas<K>::trsm("R", "U", "N", "N", &n, &k, &(Wrapper<K>::d__1), *H + nu * (m + 1), &ldh, U + nu * n, &ldv);
                        Wrapper<K>::template omatcopy<'N'>(k, n, C + nu * n, ldv, *v + nu * n, ldv);
                        Lapack<K>::mqr("R", "N", &n, &row, &k, *H + nu * (m + 1), &ldh, tau, *v + nu * n, &ldv, work, &lwork, &info);
                        Wrapper<K>::template omatcopy<'N'>(k, n, *v + nu * n, ldv, C + nu * n, ldv);
                        delete [] tau;
                        delete [] work;
                        delete [] vr;
                        delete [] alpha;
                        delete [] B;
                        delete [] A;
                    }
                }
                delete [] prod;
            }
        }
        else
            break;
    }
    if(!excluded) {
        int rem = (recycling ? (it - m) % (m - k) : it % m);
        if(rem != 0) {
            if(recycling)
                rem += k;
            std::for_each(hasConverged, hasConverged + mu, [&](short& dim) { if(dim < 0) dim = rem; });
        }
        if(mu > 1) {
            unsigned short dim = std::abs(*std::max_element(hasConverged, hasConverged + mu, [](const short& lhs, const short& rhs) { return std::abs(lhs) < std::abs(rhs); }));
            for(unsigned short i = (recycling ? k : 0); i < dim; ++i) {
                Wrapper<K>::template imatcopy<'T'>(i + 1, mu, H[i], mu, m + 1);
                for(unsigned short nu = 0; nu < mu; ++nu)
                    std::fill_n(H[i] + i + 2 + nu * (m + 1), m - i - 1, K());
            }
        }
        if(!recycling)
            updateSol(A, variant, n, x, H, s, static_cast<const K* const* const>(v + (m + 1) * (variant == 'F')), hasConverged, mu, Ax);
        else {
            computeMin(H, s + k * mu, hasConverged, mu, -1, k);
            int inc = mu;
            for(unsigned short nu = 0; nu < mu; ++nu) {
                int diff = std::abs(hasConverged[nu]) - k;
                if(diff != -k) {
                    **v = sn[nu];
                    Blas<K>::gemv(&(Wrapper<K>::transc), &n, &k, *v, C + nu * n, &ldv, *v + (k * mu + nu) * n , &i__1, &(Wrapper<K>::d__0), s + nu, &inc);
                }
            }
            MPI_Allreduce(MPI_IN_PLACE, s, k * mu, Wrapper<K>::mpi_type(), MPI_SUM, comm);
            for(unsigned short nu = 0; nu < mu; ++nu) {
                int diff = std::abs(hasConverged[nu]) - k;
                if(diff != -k) {
                    Blas<K>::gemv("N", &k, &diff, &(Wrapper<K>::d__2), H[k] + nu * (m + 1), &ldh, s + k * mu + nu, &inc, &(Wrapper<K>::d__0), s + nu, &inc);
                }
            }
            std::copy_n(U, k * ldv, v[(m + 1) * (variant == 'F')]);
            addSol(A, variant, n, x, std::distance(H[0], H[1]), s, static_cast<const K* const* const>(v + (m + 1) * (variant == 'F')), hasConverged, mu, Ax);
        }
    }
    if(verbosity > 0) {
        if(j != it + 1)
            std::cout << "GCRODR converges after " << j << " iteration" << (j > 1 ? "s" : "") << std::endl;
        else
            std::cout << "GCRODR does not converges after " << it << " iteration" << (it > 1 ? "s" : "") << std::endl;
    }
    delete [] hasConverged;
    A.clearBuffer(alloc);
    delete [] s;
    delete [] *save;
    delete [] H;
    std::cout.unsetf(std::ios_base::scientific);
    return std::min(j, it);
}
} // HPDDM
#endif // _HPDDM_GCRODR__
