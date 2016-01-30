/*
   This file is part of HPDDM.

   Author(s): Pierre Jolivet <pierre.jolivet@enseeiht.fr>
        Date: 2013-06-03

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

#ifndef _HPDDM_BDD_
#define _HPDDM_BDD_

#include "schur.hpp"

namespace HPDDM {
/* Class: Bdd
 *
 *  A class for solving problems using the BDD method.
 *
 * Template Parameters:
 *    Solver         - Solver used for the factorization of local matrices.
 *    CoarseOperator - Class of the coarse operator.
 *    S              - 'S'ymmetric or 'G'eneral coarse operator.
 *    K              - Scalar type. */
template<template<class> class Solver, template<class> class CoarseSolver, char S, class K>
class Bdd : public Schur<Solver, CoarseOperator<CoarseSolver, S, K>, K> {
    private:
        /* Variable: m
         *  Local partition of unity. */
        underlying_type<K>* _m;
    public:
        Bdd() : _m() { }
        ~Bdd() {
            delete []  _m;
        }
        /* Typedef: super
         *  Type of the immediate parent class <Schur>. */
        typedef Schur<Solver, CoarseOperator<CoarseSolver, S, K>, K> super;
        /* Function: initialize
         *  Allocates <Bdd::m> and calls <Schur::initialize>. */
        void initialize() {
            super::template initialize<false>();
            _m = new underlying_type<K>[Subdomain<K>::_dof];
        }
        void allocateSingle(K*& primal) const {
            primal = new K[Subdomain<K>::_dof];
        }
        template<unsigned short N>
        void allocateArray(K* (&array)[N]) const {
            *array = new K[N * Subdomain<K>::_dof];
            for(unsigned short i = 1; i < N; ++i)
                array[i] = *array + i * Subdomain<K>::_dof;
        }
        /* Function: buildScaling
         *
         *  Builds the local partition of unity <Bdd::m>.
         *
         * Parameters:
         *    scaling        - Type of scaling (multiplicity, stiffness or coefficient scaling).
         *    rho            - Physical local coefficients (optional). */
        void buildScaling(unsigned short& scaling, const K* const& rho = nullptr) {
            initialize();
            std::fill_n(_m, Subdomain<K>::_dof, 1.0);
            if((scaling == 2 && rho) || scaling == 1) {
                if(scaling == 1)
                    super::stiffnessScaling(super::_work);
                else
                    std::copy_n(rho + super::_bi->_m, Subdomain<K>::_dof, super::_work);
                bool alloc = Subdomain<K>::setBuffer(1, super::_structure, Subdomain<K>::_a->_n);
                Subdomain<K>::recvBuffer(super::_work);
                for(unsigned short i = 0, size = Subdomain<K>::_map.size(); i < size; ++i)
                    for(unsigned int j = 0; j < Subdomain<K>::_map[i].second.size(); ++j)
                        _m[Subdomain<K>::_map[i].second[j]] *= std::real(Subdomain<K>::_buff[size + i][j]) / std::real(Subdomain<K>::_buff[size + i][j] + _m[Subdomain<K>::_map[i].second[j]] * Subdomain<K>::_buff[i][j]);
                Subdomain<K>::clearBuffer(alloc);
            }
            else {
                scaling = 0;
                for(const pairNeighbor& neighbor : Subdomain<K>::_map)
                    for(pairNeighbor::second_type::const_reference p : neighbor.second)
                        _m[p] /= 1.0 + _m[p];
            }
        }
        /* Function: start
         *
         *  Projected Conjugate Gradient initialization.
         *
         * Template Parameter:
         *    excluded       - True if the master processes are excluded from the domain decomposition, false otherwise.
         *
         * Parameters:
         *    f              - Right-hand side.
         *    x              - Solution vector.
         *    b              - Condensed right-hand side.
         *    r              - First residual. */
        template<bool excluded>
        void start(const K* const f, K* const x, K* const b, K* r) const {
            if(super::_co) {
                if(!excluded) {
                    super::condensateEffort(f, b);
                    Subdomain<K>::exchange(b ? b : super::_structure + super::_bi->_m);
                    if(super::_ev) {
                        std::copy_n(b ? b : super::_structure + super::_bi->_m, Subdomain<K>::_dof, x);
                        Wrapper<K>::diag(Subdomain<K>::_dof, _m, x);
                        if(super::_schur) {
                            Blas<K>::gemv(&(Wrapper<K>::transc), &(Subdomain<K>::_dof), super::_co->getAddrLocal(), &(Wrapper<K>::d__1), *super::_ev, &(Subdomain<K>::_dof), x, &i__1, &(Wrapper<K>::d__0), super::_uc, &i__1);
                            super::_co->template callSolver<excluded>(super::_uc);
                            Blas<K>::gemv("N", &(Subdomain<K>::_dof), super::_co->getAddrLocal(), &(Wrapper<K>::d__1), *super::_ev, &(Subdomain<K>::_dof), super::_uc, &i__1, &(Wrapper<K>::d__0), x, &i__1);
                        }
                        else {
                            Blas<K>::gemv(&(Wrapper<K>::transc), &(Subdomain<K>::_dof), super::_co->getAddrLocal(), &(Wrapper<K>::d__1), *super::_ev + super::_bi->_m, &(Subdomain<K>::_a->_n), x, &i__1, &(Wrapper<K>::d__0), super::_uc, &i__1);
                            super::_co->template callSolver<excluded>(super::_uc);
                            Blas<K>::gemv("N", &(Subdomain<K>::_dof), super::_co->getAddrLocal(), &(Wrapper<K>::d__1), *super::_ev + super::_bi->_m, &(Subdomain<K>::_a->_n), super::_uc, &i__1, &(Wrapper<K>::d__0), x, &i__1);
                        }
                        Wrapper<K>::diag(Subdomain<K>::_dof, _m, x);
                    }
                    else {
                        std::fill_n(x, Subdomain<K>::_dof, K());
                        super::_co->template callSolver<excluded>(super::_uc);
                    }
                    Subdomain<K>::exchange(x);
                    super::applyLocalSchurComplement(x, r);
                    Subdomain<K>::exchange(r);
                    Blas<K>::axpby(Subdomain<K>::_dof, Wrapper<K>::d__1, b ? b : super::_structure + super::_bi->_m, 1, Wrapper<K>::d__2, r, 1);
                }
                else
                    super::_co->template callSolver<excluded>(super::_uc);
            }
            else if(!excluded) {
                super::condensateEffort(f, r);
                Subdomain<K>::exchange(r);
                std::fill_n(x, Subdomain<K>::_dof, K());
            }
        }
        /* Function: apply
         *
         *  Applies the global Schur complement to a single right-hand side.
         *
         * Parameters:
         *    in             - Input vector.
         *    out            - Output vector (optional). */
        void apply(K* const in, K* const out = nullptr) const {
            if(out) {
                super::applyLocalSchurComplement(in, out);
                Subdomain<K>::exchange(out);
            }
            else {
                super::applyLocalSchurComplement(in);
                Subdomain<K>::exchange(in);
            }
        }
        /* Function: precond
         *
         *  Applies the global preconditioner to a single right-hand side.
         *
         * Parameters:
         *    in             - Input vector.
         *    out            - Output vector (optional). */
        void precond(K* const in, K* const out = nullptr) const {
            Wrapper<K>::diag(Subdomain<K>::_dof, _m, in, super::_work + super::_bi->_m);
            if(!HPDDM_QR || !super::_schur) {
                std::fill_n(super::_work, super::_bi->_m, K());
                static_cast<Solver<K>*>(super::_pinv)->solve(super::_work);
            }
            else {
                if(super::_deficiency)
                    static_cast<QR<K>*>(super::_pinv)->solve(super::_work + super::_bi->_m);
                else {
                    int info;
                    Lapack<K>::potrs("L", &(Subdomain<K>::_dof), &i__1, static_cast<const K*>(super::_pinv), &(Subdomain<K>::_dof), super::_work + super::_bi->_m, &(Subdomain<K>::_dof), &info);
                }
            }
            if(out) {
                Wrapper<K>::diag(Subdomain<K>::_dof, _m, super::_work + super::_bi->_m, out);
                Subdomain<K>::exchange(out);
            }
            else {
                Wrapper<K>::diag(Subdomain<K>::_dof, _m, super::_work + super::_bi->_m, in);
                Subdomain<K>::exchange(in);
            }
        }
        /* Function: callNumfact
         *  Factorizes <Subdomain::a> or <Schur::schur> if available. */
        void callNumfact() {
            if(HPDDM_QR && super::_schur) {
                delete super::_bb;
                super::_bb = nullptr;
                if(super::_deficiency) {
                    super::_pinv = new QR<K>(Subdomain<K>::_dof, super::_schur);
                    QR<K>* qr = static_cast<QR<K>*>(super::_pinv);
                    qr->decompose();
                }
                else {
                    super::_pinv = new K[Subdomain<K>::_dof * Subdomain<K>::_dof];
                    Blas<K>::lacpy("L", &(Subdomain<K>::_dof), &(Subdomain<K>::_dof), super::_schur, &(Subdomain<K>::_dof), static_cast<K*>(super::_pinv), &(Subdomain<K>::_dof));
                    int info;
                    Lapack<K>::potrf("L", &(Subdomain<K>::_dof), static_cast<K*>(super::_pinv), &(Subdomain<K>::_dof), &info);
                }
            }
            else
                super::callNumfact();
        }
        /* Function: project
         *
         *  Projects into the coarse space.
         *
         * Template Parameters:
         *    excluded       - True if the master processes are excluded from the domain decomposition, false otherwise.
         *    trans          - 'T' if the transposed projection should be applied, 'N' otherwise.
         *
         * Parameters:
         *    in             - Input vector.
         *    out            - Output vector (optional). */
        template<bool excluded, char trans>
        void project(K* const in, K* const out = nullptr) const {
            static_assert(trans == 'T' || trans == 'N', "Unsupported value for argument 'trans'");
            if(super::_co) {
                if(!excluded) {
                    if(trans == 'N')
                        apply(in, super::_structure + super::_bi->_m);
                    if(super::_ev) {
                        if(trans == 'N')
                            Wrapper<K>::diag(Subdomain<K>::_dof, _m, super::_structure + super::_bi->_m);
                        else
                            Wrapper<K>::diag(Subdomain<K>::_dof, _m, in, super::_structure + super::_bi->_m);
                        if(super::_schur) {
                            Blas<K>::gemv(&(Wrapper<K>::transc), &(Subdomain<K>::_dof), super::_co->getAddrLocal(), &(Wrapper<K>::d__1), *super::_ev, &(Subdomain<K>::_dof), super::_structure + super::_bi->_m, &i__1, &(Wrapper<K>::d__0), super::_uc, &i__1);
                            super::_co->template callSolver<excluded>(super::_uc);
                            Blas<K>::gemv("N", &(Subdomain<K>::_dof), super::_co->getAddrLocal(), &(Wrapper<K>::d__1), *super::_ev, &(Subdomain<K>::_dof), super::_uc, &i__1, &(Wrapper<K>::d__0), super::_structure + super::_bi->_m, &i__1);
                        }
                        else {
                            Blas<K>::gemv(&(Wrapper<K>::transc), &(Subdomain<K>::_dof), super::_co->getAddrLocal(), &(Wrapper<K>::d__1), *super::_ev + super::_bi->_m, &(Subdomain<K>::_a->_n), super::_structure + super::_bi->_m, &i__1, &(Wrapper<K>::d__0), super::_uc, &i__1);
                            super::_co->callSolver(super::_uc);
                            Blas<K>::gemv("N", &(Subdomain<K>::_dof), super::_co->getAddrLocal(), &(Wrapper<K>::d__1), *super::_ev + super::_bi->_m, &(Subdomain<K>::_a->_n), super::_uc, &i__1, &(Wrapper<K>::d__0), super::_structure + super::_bi->_m, &i__1);
                        }
                    }
                    else {
                        super::_co->callSolver(super::_uc);
                        std::fill_n(super::_structure + super::_bi->_m, Subdomain<K>::_dof, K());
                    }
                    Wrapper<K>::diag(Subdomain<K>::_dof, _m, super::_structure + super::_bi->_m);
                    Subdomain<K>::exchange(super::_structure + super::_bi->_m);
                    if(trans == 'T')
                        apply(super::_structure + super::_bi->_m);
                    if(out)
                        for(unsigned int i = 0; i < Subdomain<K>::_dof; ++i)
                            out[i] = in[i] - *(super::_structure + super::_bi->_m + i);
                    else
                        Blas<K>::axpy(&(Subdomain<K>::_dof), &(Wrapper<K>::d__2), super::_structure + super::_bi->_m, &i__1, in, &i__1);
                }
                else
                    super::_co->template callSolver<excluded>(super::_uc);
            }
            else if(!excluded && out)
                std::copy_n(in, Subdomain<K>::_dof, out);
        }
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
         * See also: <Feti::buildTwo>, <Schwarz::buildTwo>.*/
        template<unsigned short excluded = 0>
        std::pair<MPI_Request, const K*>* buildTwo(const MPI_Comm& comm) {
            const Option& opt = *Option::get();
            if(!super::_schur && opt["geneo_nu"])
                super::_deficiency = opt["geneo_nu"];
            return super::template buildTwo<excluded, BddProjection<Bdd<Solver, CoarseSolver, S, K>, K>>(this, comm);
        }
        /* Function: computeSolution
         *
         *  Computes the solution after convergence of the Projected Conjugate Gradient.
         *
         * Template Parameter:
         *    excluded       - True if the master processes are excluded from the domain decomposition, false otherwise.
         *
         * Parameters:
         *    f              - Right-hand side.
         *    x              - Solution vector. */
        template<bool excluded>
        void computeSolution(const K* const f, K* const x) const {
            if(!excluded && super::_bi->_m) {
                std::copy_n(f, super::_bi->_m, x);
                Wrapper<K>::template csrmv<Wrapper<K>::I>(&(Wrapper<K>::transc), &(Subdomain<K>::_dof), &(super::_bi->_m), &(Wrapper<K>::d__2), false, super::_bi->_a, super::_bi->_ia, super::_bi->_ja, x + super::_bi->_m, &(Wrapper<K>::d__1), x);
                if(!super::_schur)
                    super::_s.solve(x);
                else {
                    std::copy_n(x, super::_bi->_m, super::_structure);
                    super::_s.solve(super::_structure);
                    std::copy_n(super::_structure, super::_bi->_m, x);
                }
            }
        }
        template<bool>
        void computeSolution(K* const* const, K* const) const { }
        /* Function: computeDot
         *
         *  Computes the dot product of two vectors.
         *
         * Template Parameter:
         *    excluded       - True if the master processes are excluded from the domain decomposition, false otherwise.
         *
         * Parameters:
         *    a              - Left-hand side.
         *    b              - Right-hand side. */
        template<bool excluded>
        void computeDot(underlying_type<K>* const val, const K* const a, const K* const b, const MPI_Comm& comm) const {
            if(!excluded) {
                Wrapper<K>::diag(Subdomain<K>::_dof, _m, a, super::_work);
                *val = std::real(Blas<K>::dot(&(Subdomain<K>::_dof), super::_work, &i__1, b, &i__1));
            }
            else
                *val = 0.0;
            MPI_Allreduce(MPI_IN_PLACE, val, 1, Wrapper<K>::mpi_underlying_type(), MPI_SUM, comm);
        }
        /* Function: getScaling
         *  Returns a constant pointer to <Bdd::m>. */
        const underlying_type<K>* getScaling() const { return _m; }
        /* Function: solveGEVP
         *
         *  Solves the GenEO problem.
         *
         * Template Parameter:
         *    L              - 'S'ymmetric or 'G'eneral transfer of the local Schur complements.
         *
         * Parameters:
         *    nu             - Number of eigenvectors requested.
         *    threshold      - Criterion for selecting the eigenpairs (optional). */
        template<char L = 'S'>
        void solveGEVP(unsigned short& nu, const underlying_type<K>& threshold = 0.0) {
            super::template solveGEVP<L>(_m, nu, threshold);
        }
};
} // HPDDM
#endif // _HPDDM_BDD_
