/*
   This file is part of HPDDM.

   Author(s): Pierre Jolivet <pierre.jolivet@enseeiht.fr>
              Frédéric Nataf <nataf@ann.jussieu.fr>
        Date: 2013-03-10

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

#ifndef _HPDDM_SCHWARZ_
#define _HPDDM_SCHWARZ_

#include <set>
#include "preconditioner.hpp"

namespace HPDDM {
/* Class: Schwarz
 *
 *  A class for solving problems using Schwarz methods that inherits from <Preconditioner>.
 *
 * Template Parameters:
 *    Solver         - Solver used for the factorization of local matrices.
 *    CoarseOperator - Class of the coarse operator.
 *    S              - 'S'ymmetric or 'G'eneral coarse operator.
 *    K              - Scalar type. */
template<template<class> class Solver, template<class> class CoarseSolver, char S, class K>
class Schwarz : public Preconditioner<Solver, CoarseOperator<CoarseSolver, S, K>, K> {
    public:
        /* Enum: Prcndtnr
         *
         *  Defines the Schwarz method used as a preconditioner.
         *
         * NO           - No preconditioner.
         * SY           - Symmetric preconditioner, e.g. Additive Schwarz method.
         * GE           - Nonsymmetric preconditioner, e.g. Restricted Additive Schwarz method.
         * OS           - Optimized symmetric preconditioner, e.g. Optimized Schwarz method.
         * OG           - Optimized nonsymmetric preconditioner, e.g. Optimized Restricted Additive Schwarz method. */
        enum class Prcndtnr : char {
            NO, SY, GE, OS, OG
        };
    private:
        /* Variable: d
         *  Local partition of unity. */
        const underlying_type<K>* _d;
        /* Variable: type
         *  Type of <Prcndtnr> used in <Schwarz::apply> and <Schwarz::deflation>. */
        Prcndtnr               _type;
    public:
        Schwarz() : _d() { }
        ~Schwarz() { _d = nullptr; }
        /* Typedef: super
         *  Type of the immediate parent class <Preconditioner>. */
        typedef Preconditioner<Solver, CoarseOperator<CoarseSolver, S, K>, K> super;
        /* Function: initialize
         *  Sets <Schwarz::d>. */
        template<class Container = std::vector<int>>
        void initialize(underlying_type<K>* const& d) {
            _d = d;
        }
        /* Function: callNumfact
         *  Factorizes <Subdomain::a> or another user-supplied matrix, useful for <Prcndtnr::OS> and <Prcndtnr::OG>. */
        template<char N = HPDDM_NUMBERING>
        void callNumfact(MatrixCSR<K>* const& A = nullptr) {
            Option& opt = *Option::get();
            if(A != nullptr) {
                if(opt["schwarz_method"] == 1)
                    _type = Prcndtnr::OS;
                else
                    _type = Prcndtnr::OG;
            }
            else {
                if(opt["schwarz_method"] == 3)
                    _type = Prcndtnr::SY;
                else if(opt["schwarz_method"] == 5)
                    _type = Prcndtnr::NO;
                else {
                    _type = Prcndtnr::GE;
                    opt["schwarz_method"] = 0;
                }
            }
            unsigned short reuse = opt.val<unsigned short>("reuse_preconditioner", 0);
            if(reuse <= 1)
                super::_s.template numfact<N>(_type == Prcndtnr::OS || _type == Prcndtnr::OG ? A : Subdomain<K>::_a, _type == Prcndtnr::OS ? true : false);
            if(reuse >= 1)
                opt["reuse_preconditioner"] += 1;
        }
        void setMatrix(MatrixCSR<K>* const& a) {
            bool fact = super::setMatrix(a) && _type != Prcndtnr::OS && _type != Prcndtnr::OG;
            if(fact) {
                using type = alias<Solver<K>>;
                super::_s.~type();
                super::_s.numfact(a);
                Option& opt = *Option::get();
                if(opt.val<unsigned short>("reuse_preconditioner", 0) >= 1)
                    opt["reuse_preconditioner"] = 1;
            }
        }
        /* Function: multiplicityScaling
         *
         *  Builds the multiplicity scaling.
         *
         * Parameter:
         *    d              - Array of values. */
        void multiplicityScaling(underlying_type<K>* const d) const {
            Subdomain<K>::setBuffer(1);
            for(unsigned short i = 0, size = Subdomain<K>::_map.size(); i < size; ++i) {
                underlying_type<K>* const recv = reinterpret_cast<underlying_type<K>*>(Subdomain<K>::_buff[i]);
                underlying_type<K>* const send = reinterpret_cast<underlying_type<K>*>(Subdomain<K>::_buff[size + i]);
                MPI_Irecv(recv, Subdomain<K>::_map[i].second.size(), Wrapper<K>::mpi_underlying_type(), Subdomain<K>::_map[i].first, 0, Subdomain<K>::_communicator, Subdomain<K>::_rq + i);
                Wrapper<underlying_type<K>>::gthr(Subdomain<K>::_map[i].second.size(), d, send, Subdomain<K>::_map[i].second.data());
                MPI_Isend(send, Subdomain<K>::_map[i].second.size(), Wrapper<K>::mpi_underlying_type(), Subdomain<K>::_map[i].first, 0, Subdomain<K>::_communicator, Subdomain<K>::_rq + size + i);
            }
            std::fill_n(d, Subdomain<K>::_dof, 1.0);
            for(unsigned short i = 0, size = Subdomain<K>::_map.size(); i < size; ++i) {
                int index;
                MPI_Waitany(size, Subdomain<K>::_rq, &index, MPI_STATUS_IGNORE);
                underlying_type<K>* const recv = reinterpret_cast<underlying_type<K>*>(Subdomain<K>::_buff[index]);
                underlying_type<K>* const send = reinterpret_cast<underlying_type<K>*>(Subdomain<K>::_buff[size + index]);
                for(unsigned int j = 0; j < Subdomain<K>::_map[index].second.size(); ++j) {
                    if(std::abs(send[j]) < HPDDM_EPS)
                        d[Subdomain<K>::_map[index].second[j]] = 0.0;
                    else
                        d[Subdomain<K>::_map[index].second[j]] /= 1.0 + d[Subdomain<K>::_map[index].second[j]] * recv[j] / send[j];
                }
            }
            MPI_Waitall(Subdomain<K>::_map.size(), Subdomain<K>::_rq + Subdomain<K>::_map.size(), MPI_STATUSES_IGNORE);
            delete [] *Subdomain<K>::_buff;
        }
        /* Function: getScaling
         *  Returns a constant pointer to <Schwarz::d>. */
        const underlying_type<K>* getScaling() const { return _d; }
        /* Function: scaledExchange */
        template<bool allocate = false>
        void scaledExchange(K* const x, const unsigned short& mu = 1) const {
            bool free = false;
            if(allocate)
                free = Subdomain<K>::setBuffer(mu);
            Wrapper<K>::diag(Subdomain<K>::_dof, _d, x, mu);
            Subdomain<K>::exchange(x, mu);
            if(allocate)
                Subdomain<K>::clearBuffer(free);
        }
        /* Function: deflation
         *
         *  Computes a coarse correction.
         *
         * Template parameter:
         *    excluded       - True if the master processes are excluded from the domain decomposition, false otherwise. 
         *
         * Parameters:
         *    in             - Input vectors.
         *    out            - Output vectors.
         *    mu             - Number of vectors. */
        template<bool excluded>
        void deflation(const K* const in, K* const out, const unsigned short& mu) const {
            if(excluded)
                super::_co->template callSolver<excluded>(super::_uc, mu);
            else {
                Wrapper<K>::diag(Subdomain<K>::_dof, _d, in, out, mu);                                                                                                                                                                                      // out = D in
                int tmp = mu;
                Blas<K>::gemm(&(Wrapper<K>::transc), "N", super::getAddrLocal(), &tmp, &(Subdomain<K>::_dof), &(Wrapper<K>::d__1), *super::_ev, &(Subdomain<K>::_dof), out, &(Subdomain<K>::_dof), &(Wrapper<K>::d__0), super::_uc, super::getAddrLocal()); // _uc = _ev^T D in
                super::_co->template callSolver<excluded>(super::_uc, mu);                                                                                                                                                                                  // _uc = E \ _ev^T D in
                Blas<K>::gemm("N", "N", &(Subdomain<K>::_dof), &tmp, super::getAddrLocal(), &(Wrapper<K>::d__1), *super::_ev, &(Subdomain<K>::_dof), super::_uc, super::getAddrLocal(), &(Wrapper<K>::d__0), out, &(Subdomain<K>::_dof));                   // out = _ev E \ _ev^T D in
                scaledExchange(out, mu);
            }
        }
#if HPDDM_ICOLLECTIVE
        /* Function: Ideflation
         *
         *  Computes the first part of a coarse correction asynchronously.
         *
         * Template parameter:
         *    excluded       - True if the master processes are excluded from the domain decomposition, false otherwise. 
         *
         * Parameters:
         *    in             - Input vector.
         *    out            - Output vector.
         *    rq             - MPI request to check completion of the MPI transfers. */
        template<bool excluded>
        void Ideflation(const K* const in, K* const out, MPI_Request* rq) const {
            if(excluded)
                super::_co->template IcallSolver<excluded>(super::_uc, rq);
            else {
                Wrapper<K>::diag(Subdomain<K>::_dof, _d, in, out, mu);
                Blas<K>::gemv(&(Wrapper<K>::transc), &(Subdomain<K>::_dof), super::getAddrLocal(), &(Wrapper<K>::d__1), *super::_ev, &(Subdomain<K>::_dof), out, &i__1, &(Wrapper<K>::d__0), super::_uc, &i__1);
                super::_co->template IcallSolver<excluded>(super::_uc, rq);
            }
        }
#endif // HPDDM_ICOLLECTIVE
        /* Function: buildTwo
         *
         *  Assembles and factorizes the coarse operator by calling <Preconditioner::buildTwo>.
         *
         * Template Parameter:
         *    excluded       - Greater than 0 if the master processes are excluded from the domain decomposition, equal to 0 otherwise.
         *
         * Parameter:
         *    comm           - Global MPI communicator.
         *
         * See also: <Bdd::buildTwo>, <Feti::buildTwo>. */
        template<unsigned short excluded = 0>
        std::pair<MPI_Request, const K*>* buildTwo(const MPI_Comm& comm) {
            return super::template buildTwo<excluded, MatrixMultiplication<Schwarz<Solver, CoarseSolver, S, K>, K>>(this, comm);
        }
        template<bool excluded = false>
        void start(const K* const b, K* const x, const unsigned short& mu = 1) const {
            if(super::_co) {
                super::start(mu);
                if(Option::get()->val("schwarz_coarse_correction", -1) == 2)
                    deflation<excluded>(b, x, mu);
            }
        }
        /* Function: apply
         *
         *  Applies the global Schwarz preconditioner.
         *
         * Template Parameter:
         *    excluded       - Greater than 0 if the master processes are excluded from the domain decomposition, equal to 0 otherwise.
         *
         * Parameters:
         *    in             - Input vectors, modified internally if no workspace array is specified !
         *    out            - Output vectors.
         *    mu             - Number of vectors.
         *    work           - Workspace array. */
        template<bool excluded = false>
        void apply(const K* const in, K* const out, const unsigned short& mu = 1, K* work = nullptr) const {
            const int correction = Option::get()->val("schwarz_coarse_correction", -1);
            if(!super::_co || correction == -1) {
                if(_type == Prcndtnr::NO)
                    std::copy_n(in, mu * Subdomain<K>::_dof, out);
                else if(_type == Prcndtnr::GE || _type == Prcndtnr::OG) {
                    if(!excluded) {
                        super::_s.solve(in, out, mu);
                        scaledExchange(out, mu);         // out = D A \ in
                    }
                }
                else {
                    if(!excluded) {
                        if(_type == Prcndtnr::OS) {
                            Wrapper<K>::diag(Subdomain<K>::_dof, _d, in, out, mu);
                            super::_s.solve(out, mu);
                            Wrapper<K>::diag(Subdomain<K>::_dof, _d, out, mu);
                        }
                        else
                            super::_s.solve(in, out, mu);
                        Subdomain<K>::exchange(out, mu); // out = A \ in
                    }
                }
            }
            else {
                int tmp = mu * Subdomain<K>::_dof;
                if(!work)
                    work = const_cast<K*>(in);
                else
                    std::copy_n(in, tmp, work);
                if(correction == 1) {
#if HPDDM_ICOLLECTIVE
                    MPI_Request rq[2];
                    Ideflation<excluded>(in, out, mu, rq);
                    if(!excluded) {
                        super::_s.solve(work, mu);                                                                                                                                                         // out = A \ in
                        MPI_Waitall(2, rq, MPI_STATUSES_IGNORE);
                        for(unsigned short nu = 0; nu < mu; ++nu)
                            Blas<K>::gemv("N", &(Subdomain<K>::_dof), super::getAddrLocal(), &(Wrapper<K>::d__1), *super::_ev, &(Subdomain<K>::_dof), super::_uc + nu * super::getLocal(), &i__1, &(Wrapper<K>::d__0), out + nu * Subdomain<K>::_dof, &i__1); // out = Z E \ Z^T in
                        Blas<K>::axpy(&tmp, &(Wrapper<K>::d__1), work, &i__1, out, &i__1);
                        scaledExchange(out, mu);                                                                                                                                                                                                              // out = Z E \ Z^T in + A \ in
                    }
                    else
                        MPI_Wait(rq + 1, MPI_STATUS_IGNORE);
#else
                    deflation<excluded>(in, out, mu);
                    if(!excluded) {
                        super::_s.solve(work, mu);
                        Blas<K>::axpy(&tmp, &(Wrapper<K>::d__1), work, &i__1, out, &i__1);
                        scaledExchange(out, mu);
                    }
#endif // HPDDM_ICOLLECTIVE
                }
                else if(correction == 2) {
                    if(_type == Prcndtnr::OS)
                        Wrapper<K>::diag(Subdomain<K>::_dof, _d, work, mu);
                    super::_s.solve(work, out, mu);
                    scaledExchange(out, mu);
                    GMV(out, work, mu);
                    deflation<excluded>(nullptr, work, mu);
                    Blas<K>::axpy(&tmp, &(Wrapper<K>::d__2), work, &i__1, out, &i__1);
                }
                else {
                    deflation<excluded>(in, out, mu);                                                                  // out = Z E \ Z^T in
                    if(!excluded) {
                        Wrapper<K>::csrmm("N", &(Subdomain<K>::_dof), &(tmp = mu), &(Subdomain<K>::_dof), &(Wrapper<K>::d__2), Subdomain<K>::_a->_sym, Subdomain<K>::_a->_a, Subdomain<K>::_a->_ia, Subdomain<K>::_a->_ja, out, &(Subdomain<K>::_dof), &(Wrapper<K>::d__1), work, &(Subdomain<K>::_dof));
                        scaledExchange(work, mu);                                                                      //  in = (I - A Z E \ Z^T) in
                        if(_type == Prcndtnr::OS)
                            Wrapper<K>::diag(Subdomain<K>::_dof, _d, work, mu);
                        super::_s.solve(work, mu);
                        scaledExchange(work, mu);                                                                      //  in = D A \ (I - A Z E \ Z^T) in
                        Blas<K>::axpy(&(tmp = mu * Subdomain<K>::_dof), &(Wrapper<K>::d__1), work, &i__1, out, &i__1); // out = D A \ (I - A Z E \ Z^T) in + Z E \ Z^T in
                    }
                }
            }
        }
        /* Function: scaleIntoOverlap
         *
         *  Scales the input matrix using <Schwarz::d> on the overlap and sets the output matrix to zero elsewhere.
         *
         * Parameters:
         *    A              - Input matrix.
         *    B              - Output matrix used in GenEO.
         *
         * See also: <Schwarz::solveGEVP>. */
        template<char N = HPDDM_NUMBERING>
        void scaleIntoOverlap(const MatrixCSR<K>* const& A, MatrixCSR<K>*& B) const {
            std::set<unsigned int> intoOverlap;
            for(const pairNeighbor& neighbor : Subdomain<K>::_map)
                for(unsigned int i : neighbor.second)
                    if(_d[i] > HPDDM_EPS)
                        intoOverlap.insert(i);
            std::vector<std::vector<std::pair<unsigned int, K>>> tmp(intoOverlap.size());
            unsigned int k, iPrev = 0;
#pragma omp parallel for schedule(static, HPDDM_GRANULARITY) reduction(+ : iPrev)
            for(k = 0; k < intoOverlap.size(); ++k) {
                auto it = std::next(intoOverlap.cbegin(), k);
                tmp[k].reserve(A->_ia[*it + 1] - A->_ia[*it]);
                for(unsigned int j = A->_ia[*it] - (N == 'F'); j < A->_ia[*it + 1] - (N == 'F'); ++j) {
                    K value = _d[*it] * _d[A->_ja[j] - (N == 'F')] * A->_a[j];
                    if(std::abs(value) > HPDDM_EPS && intoOverlap.find(A->_ja[j] - (N == 'F')) != intoOverlap.cend())
                        tmp[k].emplace_back(A->_ja[j], value);
                }
                iPrev += tmp[k].size();
            }
            int nnz = iPrev;
            if(B != nullptr)
                delete B;
            B = new MatrixCSR<K>(Subdomain<K>::_dof, Subdomain<K>::_dof, nnz, A->_sym);
            nnz = iPrev = k = 0;
            for(unsigned int i : intoOverlap) {
                std::fill(B->_ia + iPrev, B->_ia + i + 1, nnz + (N == 'F'));
                for(const std::pair<unsigned int, K>& p : tmp[k]) {
                    B->_ja[nnz] = p.first;
                    B->_a[nnz++] = p.second;
                }
                ++k;
                iPrev = i + 1;
            }
            std::fill(B->_ia + iPrev, B->_ia + Subdomain<K>::_dof + 1, nnz + (N == 'F'));
        }
        /* Function: solveGEVP
         *
         *  Solves the generalized eigenvalue problem Ax = l Bx.
         *
         * Parameters:
         *    A              - Left-hand side matrix.
         *    B              - Right-hand side matrix (optional).
         *    nu             - Number of eigenvectors requested.
         *    threshold      - Precision of the eigensolver. */
        template<template<class> class Eps>
        void solveGEVP(MatrixCSR<K>* const& A, unsigned short& nu, const underlying_type<K>& threshold, MatrixCSR<K>* const& B = nullptr, const MatrixCSR<K>* const& pattern = nullptr) {
            Eps<K> evp(threshold, Subdomain<K>::_dof, nu);
#ifndef PY_MAJOR_VERSION
            bool free = pattern ? pattern->sameSparsity(A) : Subdomain<K>::_a->sameSparsity(A);
#else
            constexpr bool free = false;
#endif
            MatrixCSR<K>* rhs = nullptr;
            if(B)
                rhs = B;
            else
                scaleIntoOverlap(A, rhs);
            evp.template solve<Solver>(A, rhs, super::_ev, Subdomain<K>::_communicator, free ? &(super::_s) : nullptr);
            if(rhs != B)
                delete rhs;
            if(free) {
                A->_ia = nullptr;
                A->_ja = nullptr;
            }
            (*Option::get())["geneo_nu"] = nu = evp.getNu();
            const int n = Subdomain<K>::_dof;
            std::for_each(super::_ev, super::_ev + nu, [&](K* const v) { std::replace_if(v, v + n, [](K x) { return std::abs(x) < 1.0 / (HPDDM_EPS * HPDDM_PEN); }, K()); });
        }
        template<bool sorted = true, bool scale = false>
        void interaction(std::vector<const MatrixCSR<K>*>& blocks) const {
            Subdomain<K>::template interaction<HPDDM_NUMBERING, sorted, scale>(blocks, _d);
        }
        /* Function: GMV
         *
         *  Computes a global sparse matrix-vector product.
         *
         * Parameters:
         *    in             - Input vector.
         *    out            - Output vector. */
        void GMV(const K* const in, K* const out, const int& mu = 1) const {
#if 0
            K* tmp = new K[mu * Subdomain<K>::_dof];
            Wrapper<K>::diag(Subdomain<K>::_dof, _d, in, tmp, mu);
            if(HPDDM_NUMBERING == Wrapper<K>::I)
                Wrapper<K>::csrmm(Subdomain<K>::_a->_sym, &(Subdomain<K>::_dof), &mu, Subdomain<K>::_a->_a, Subdomain<K>::_a->_ia, Subdomain<K>::_a->_ja, tmp, out);
            else if(Subdomain<K>::_a->_ia[Subdomain<K>::_dof] == Subdomain<K>::_a->_nnz)
                Wrapper<K>::template csrmm<'C'>(Subdomain<K>::_a->_sym, &(Subdomain<K>::_dof), &mu, Subdomain<K>::_a->_a, Subdomain<K>::_a->_ia, Subdomain<K>::_a->_ja, tmp, out);
            else
                Wrapper<K>::template csrmm<'F'>(Subdomain<K>::_a->_sym, &(Subdomain<K>::_dof), &mu, Subdomain<K>::_a->_a, Subdomain<K>::_a->_ia, Subdomain<K>::_a->_ja, tmp, out);
            delete [] tmp;
            Subdomain<K>::exchange(out, mu);
#else
            if(HPDDM_NUMBERING == Wrapper<K>::I)
                Wrapper<K>::csrmm(Subdomain<K>::_a->_sym, &(Subdomain<K>::_dof), &mu, Subdomain<K>::_a->_a, Subdomain<K>::_a->_ia, Subdomain<K>::_a->_ja, in, out);
            else if(Subdomain<K>::_a->_ia[Subdomain<K>::_dof] == Subdomain<K>::_a->_nnz)
                Wrapper<K>::template csrmm<'C'>(Subdomain<K>::_a->_sym, &(Subdomain<K>::_dof), &mu, Subdomain<K>::_a->_a, Subdomain<K>::_a->_ia, Subdomain<K>::_a->_ja, in, out);
            else
                Wrapper<K>::template csrmm<'F'>(Subdomain<K>::_a->_sym, &(Subdomain<K>::_dof), &mu, Subdomain<K>::_a->_a, Subdomain<K>::_a->_ia, Subdomain<K>::_a->_ja, in, out);
            scaledExchange(out, mu);
#endif
        }
        /* Function: computeError
         *
         *  Computes the Euclidean norm of a right-hand side and of the difference between a solution vector and a right-hand side.
         *
         * Parameters:
         *    x              - Solution vector.
         *    f              - Right-hand side.
         *    storage        - Array to store both values.
         *
         * See also: <Schur::computeError>. */
        void computeError(const K* const x, const K* const f, underlying_type<K>* const storage, const unsigned short& mu = 1) const {
            int dim = mu * Subdomain<K>::_dof;
            K* tmp = new K[dim];
            bool alloc = Subdomain<K>::setBuffer(mu);
            GMV(x, tmp, mu);
            Subdomain<K>::clearBuffer(alloc);
            Blas<K>::axpy(&dim, &(Wrapper<K>::d__2), f, &i__1, tmp, &i__1);
            std::fill_n(storage, 2 * mu, 0.0);
            for(unsigned int i = 0; i < Subdomain<K>::_dof; ++i) {
                bool isBoundaryCond = true;
                unsigned int stop;
                if(!Subdomain<K>::_a->_sym)
                    stop = std::distance(Subdomain<K>::_a->_ja, std::upper_bound(Subdomain<K>::_a->_ja + Subdomain<K>::_a->_ia[i], Subdomain<K>::_a->_ja + Subdomain<K>::_a->_ia[i + 1], i));
                else
                    stop = Subdomain<K>::_a->_ia[i + 1];
                if(std::abs(Subdomain<K>::_a->_a[stop - 1]) > HPDDM_EPS * HPDDM_PEN)
                    continue;
                for(unsigned int j = Subdomain<K>::_a->_ia[i]; j < stop && isBoundaryCond; ++j) {
                    if(i != Subdomain<K>::_a->_ja[j] && std::abs(Subdomain<K>::_a->_a[j]) > HPDDM_EPS)
                        isBoundaryCond = false;
                    else if(i == Subdomain<K>::_a->_ja[j] && std::abs(Subdomain<K>::_a->_a[j] - K(1.0)) > HPDDM_EPS)
                        isBoundaryCond = false;
                }
                for(unsigned short nu = 0; nu < mu; ++nu) {
                    if(!isBoundaryCond)
                        storage[2 * nu + 1] += _d[i] * std::norm(tmp[nu * Subdomain<K>::_dof + i]);
                    if(std::abs(f[nu * Subdomain<K>::_dof + i]) > HPDDM_EPS * HPDDM_PEN)
                        storage[2 * nu] += _d[i] * std::norm(f[nu * Subdomain<K>::_dof + i] / K(HPDDM_PEN));
                    else
                        storage[2 * nu] += _d[i] * std::norm(f[nu * Subdomain<K>::_dof + i]);
                }
            }
            delete [] tmp;
            MPI_Allreduce(MPI_IN_PLACE, storage, 2 * mu, Wrapper<K>::mpi_underlying_type(), MPI_SUM, Subdomain<K>::_communicator);
            std::for_each(storage, storage + 2 * mu, [](underlying_type<K>& b) { b = std::sqrt(b); });
        }
        template<char N = HPDDM_NUMBERING>
        void distributedNumbering(unsigned int* const in, unsigned int& first, unsigned int& last, unsigned int& global) const {
            Subdomain<K>::template globalMapping<N>(in, in + Subdomain<K>::_dof, first, last, global, _d);
        }
        bool distributedCSR(unsigned int* const num, unsigned int first, unsigned int last, int*& ia, int*& ja, K*& c) const {
            return Subdomain<K>::distributedCSR(num, first, last, ia, ja, c, Subdomain<K>::_a);
        }
};
} // HPDDM
#endif // _HPDDM_SCHWARZ_
