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

#ifndef _HPDDM_FETI_
#define _HPDDM_FETI_

#include "schur.hpp"

namespace HPDDM {
/* Class: Feti
 *
 *  A class for solving problems using the FETI method.
 *
 * Template Parameters:
 *    Solver         - Solver used for the factorization of local matrices.
 *    CoarseOperator - Class of the coarse operator.
 *    S              - 'S'ymmetric or 'G'eneral coarse operator.
 *    K              - Scalar type. */
template<template<class> class Solver, template<class> class CoarseSolver, char S, class K, FetiPrcndtnr P>
class Feti : public Schur<Solver, CoarseOperator<CoarseSolver, S, K>, K> {
    private:
        /* Variable: primal
         *  Storage for local primal unknowns. */
        K*              _primal;
        /* Variable: dual
         *  Storage for local dual unknowns. */
        K**               _dual;
        /* Variable: m
         *  Local partition of unity. */
        underlying_type<K>** _m;
        /* Function: A
         *
         *  Jump operator.
         *
         * Template Parameters:
         *    trans          - 'T' if the transposed jump operator should be applied, 'N' otherwise.
         *    scale          - True if the unknowns should be scale by <Feti::m>, false otherwise.
         *
         * Parameters:
         *    primal         - Primal unknowns.
         *    dual           - Dual unknowns. */
        template<char trans, bool scale>
        void A(K* const primal, K* const* const dual) const {
            static_assert(trans == 'T' || trans == 'N', "Unsupported value for argument 'trans'");
            if(trans == 'T') {
                std::fill_n(primal, Subdomain<K>::_dof, K());
                for(unsigned short i = 0; i < super::_signed; ++i)
                    for(unsigned int j = 0; j < Subdomain<K>::_map[i].second.size(); ++j)
                        primal[Subdomain<K>::_map[i].second[j]] -= scale ? _m[i][j] * dual[i][j] : dual[i][j];
                for(unsigned short i = super::_signed; i < Subdomain<K>::_map.size(); ++i)
                    for(unsigned int j = 0; j < Subdomain<K>::_map[i].second.size(); ++j)
                        primal[Subdomain<K>::_map[i].second[j]] += scale ? _m[i][j] * dual[i][j] : dual[i][j];
            }
            else {
                for(unsigned short i = 0; i < super::_signed; ++i) {
                    MPI_Irecv(Subdomain<K>::_buff[i], Subdomain<K>::_map[i].second.size(), Wrapper<K>::mpi_type(), Subdomain<K>::_map[i].first, 0, Subdomain<K>::_communicator, Subdomain<K>::_rq + i);
                    for(unsigned int j = 0; j < Subdomain<K>::_map[i].second.size(); ++j)
                        dual[i][j] = -(scale ? _m[i][j] * primal[Subdomain<K>::_map[i].second[j]] : primal[Subdomain<K>::_map[i].second[j]]);
                    MPI_Isend(dual[i], Subdomain<K>::_map[i].second.size(), Wrapper<K>::mpi_type(), Subdomain<K>::_map[i].first, 0, Subdomain<K>::_communicator, Subdomain<K>::_rq + Subdomain<K>::_map.size() + i);
                }
                for(unsigned short i = super::_signed; i < Subdomain<K>::_map.size(); ++i) {
                    MPI_Irecv(Subdomain<K>::_buff[i], Subdomain<K>::_map[i].second.size(), Wrapper<K>::mpi_type(), Subdomain<K>::_map[i].first, 0, Subdomain<K>::_communicator, Subdomain<K>::_rq + i);
                    for(unsigned int j = 0; j < Subdomain<K>::_map[i].second.size(); ++j)
                        dual[i][j] =  (scale ? _m[i][j] * primal[Subdomain<K>::_map[i].second[j]] : primal[Subdomain<K>::_map[i].second[j]]);
                    MPI_Isend(dual[i], Subdomain<K>::_map[i].second.size(), Wrapper<K>::mpi_type(), Subdomain<K>::_map[i].first, 0, Subdomain<K>::_communicator, Subdomain<K>::_rq + Subdomain<K>::_map.size() + i);
                }
                MPI_Waitall(2 * Subdomain<K>::_map.size(), Subdomain<K>::_rq, MPI_STATUSES_IGNORE);
                Blas<K>::axpy(&(super::_mult), &(Wrapper<K>::d__1), Subdomain<K>::_buff[0], &i__1, *dual, &i__1);
            }
        }
        template<class U, typename std::enable_if<!Wrapper<U>::is_complex>::type* = nullptr>
        void allocate(U**& dual, underlying_type<U>**& m) {
            static_assert(std::is_same<U, K>::value, "Wrong types");
            dual = new U*[2 * Subdomain<K>::_map.size()];
            m    = dual + Subdomain<K>::_map.size();
        }
        template<class U, typename std::enable_if<Wrapper<U>::is_complex>::type* = nullptr>
        void allocate(U**& dual, underlying_type<U>**& m) {
            static_assert(std::is_same<U, K>::value, "Wrong types");
            dual = new U*[Subdomain<K>::_map.size()];
            m    = new underlying_type<U>*[Subdomain<K>::_map.size()];
        }
    public:
        Feti() : _primal(), _dual(), _m() { }
        ~Feti() {
            delete [] super::_schur;
            super::_schur = nullptr;
            if(_m)
                delete [] *_m;
            if(Wrapper<K>::is_complex)
                delete []  _m;
            delete [] _dual;
        }
        /* Typedef: super
         *  Type of the immediate parent class <Schur>. */
        typedef Schur<Solver, CoarseOperator<CoarseSolver, S, K>, K> super;
        /* Function: initialize
         *  Allocates <Feti::primal>, <Feti::dual>, and <Feti::m> and calls <Schur::initialize>. */
        void initialize() {
            super::template initialize<true>();
            _primal = super::_structure + super::_bi->_m;
            allocate(_dual, _m);
            *_dual  = super::_work;
            *_m     = new underlying_type<K>[super::_mult];
            for(unsigned short i = 1; i < Subdomain<K>::_map.size(); ++i) {
                _dual[i] = _dual[i - 1] + Subdomain<K>::_map[i - 1].second.size();
                _m   [i] = _m[i - 1]    + Subdomain<K>::_map[i - 1].second.size();
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
        void start(const K* const f, K* const x, K* const* const l, K* const* const r) const {
            Solver<K>* p = static_cast<Solver<K>*>(super::_pinv);
            if(super::_co) {
                if(!excluded) {
                    if(super::_ev) {
                        if(super::_schur) {
                            super::condensateEffort(f, nullptr);
                            Blas<K>::gemv(&(Wrapper<K>::transc), &(Subdomain<K>::_dof), super::_co->getAddrLocal(), &(Wrapper<K>::d__1), *super::_ev, &(Subdomain<K>::_dof), super::_structure + super::_bi->_m, &i__1, &(Wrapper<K>::d__0), super::_uc, &i__1); //     _uc = R_b g
                            super::_co->template callSolver<excluded>(super::_uc);                                                                                                                                                                               //     _uc = (G Q G^T) \ R_b g
                            Blas<K>::gemv("N", &(Subdomain<K>::_dof), super::_co->getAddrLocal(), &(Wrapper<K>::d__1), *super::_ev, &(Subdomain<K>::_dof), super::_uc, &i__1, &(Wrapper<K>::d__0), _primal, &i__1);                                              // _primal = R_b (G Q G^T) \ R f
                        }
                        else {
                            Blas<K>::gemv(&(Wrapper<K>::transc), &(Subdomain<K>::_a->_n), super::_co->getAddrLocal(), &(Wrapper<K>::d__1), *super::_ev, &(Subdomain<K>::_a->_n), f, &i__1, &(Wrapper<K>::d__0), super::_uc, &i__1);                              //     _uc = R f
                            super::_co->template callSolver<excluded>(super::_uc);                                                                                                                                                                               //     _uc = (G Q G^T) \ R f
                            Blas<K>::gemv("N", &(Subdomain<K>::_dof), super::_co->getAddrLocal(), &(Wrapper<K>::d__1), *super::_ev + super::_bi->_m, &(Subdomain<K>::_a->_n), super::_uc, &i__1, &(Wrapper<K>::d__0), _primal, &i__1);                           // _primal = R_b (G Q G^T) \ R f
                        }
                    }
                    else {
                        super::_co->template callSolver<excluded>(super::_uc);
                        std::fill_n(_primal, Subdomain<K>::_dof, K());
                    }
                    A<'N', 0>(_primal, l);                                                            //       l = A R_b (G Q G^T) \ R f
                    precond(l);                                                                       //       l = Q A R_b (G Q G^T) \ R f
                    A<'T', 0>(_primal, l);                                                            // _primal = A^T Q A R_b (G Q G^T) \ R f
                    std::fill_n(super::_structure, super::_bi->_m, K());
                    p->solve(super::_structure);                                                      // _primal = S \ A^T Q A R_b (G Q G^T) \ R f
                }
                else
                    super::_co->template callSolver<excluded>(super::_uc);
            }
            if(!excluded) {
                p->solve(f, x);                                                                       //       x = S \ f
                if(!super::_co) {
                    A<'N', 0>(x + super::_bi->_m, r);                                                 //       r = A S \ f
                    std::fill_n(*l, super::_mult, K());                                               //       l = 0
                }
                else {
                    Blas<K>::axpby(Subdomain<K>::_dof, 1.0, x + super::_bi->_m, 1, -1.0, _primal, 1); // _primal = S \ (f - A^T Q A R_b (G Q G^T) \ R f)
                    A<'N', 0>(_primal, r);                                                            //       r = A S \ (f - A^T Q A R_b (G Q G^T) \ R f)
                    project<excluded, 'T'>(r);                                                        //       r = P^T r
                }
            }
            else if(super::_co)
                project<excluded, 'T'>(r);
        }
        /* Function: allocateSingle
         *
         *  Allocates a single Lagrange multiplier.
         *
         * Parameter:
         *    mult           - Reference to a Lagrange multiplier. */
        void allocateSingle(K**& mult) const {
            mult  = new K*[Subdomain<K>::_map.size()];
            *mult = new K[super::_mult];
            for(unsigned short i = 1; i < Subdomain<K>::_map.size(); ++i)
                mult[i] = mult[i - 1] + Subdomain<K>::_map[i - 1].second.size();
        }
        /* Function: allocateArray
         *
         *  Allocates an array of multiple Lagrange multipliers.
         *
         * Template Parameter:
         *    N              - Size of the array.
         *
         * Parameter:
         *    array          - Reference to an array of Lagrange multipliers. */
        template<unsigned short N>
        void allocateArray(K** (&array)[N]) const {
            *array  = new K*[N * Subdomain<K>::_map.size()];
            **array = new K[N * super::_mult];
            for(unsigned short i = 0; i < N; ++i) {
                array[i]  = *array + i * Subdomain<K>::_map.size();
                *array[i] = **array + i * super::_mult;
                for(unsigned short j = 1; j < Subdomain<K>::_map.size(); ++j)
                    array[i][j] = array[i][j - 1] + Subdomain<K>::_map[j - 1].second.size();
            }
        }
        /* Function: buildScaling
         *
         *  Builds the local partition of unity <Feti::m>.
         *
         * Parameters:
         *    scaling        - Type of scaling (multiplicity, stiffness or coefficient scaling).
         *    rho            - Physical local coefficients (optional). */
        void buildScaling(unsigned short& scaling, const K* const& rho = nullptr) {
            initialize();
            std::vector<std::pair<unsigned short, unsigned int>>* array = new std::vector<std::pair<unsigned short, unsigned int>>[Subdomain<K>::_dof];
            for(const pairNeighbor& neighbor: Subdomain<K>::_map)
                for(unsigned int j = 0; j < neighbor.second.size(); ++j)
                    array[neighbor.second[j]].emplace_back(neighbor.first, j);
            if((scaling == 2 && rho) || scaling == 1) {
                if(scaling == 1)
                    super::stiffnessScaling(_primal);
                else
                    std::copy_n(rho + super::_bi->_m, Subdomain<K>::_dof, _primal);
                bool alloc = Subdomain<K>::setBuffer(1, super::_work, super::_mult + super::_bi->_m);
                Subdomain<K>::exchange(_primal);
                for(unsigned short i = 0; i < Subdomain<K>::_map.size(); ++i)
                    for(unsigned int j = 0; j < Subdomain<K>::_map[i].second.size(); ++j)
                        _m[i][j] = std::real(Subdomain<K>::_buff[i][j] / _primal[Subdomain<K>::_map[i].second[j]]);
                Subdomain<K>::clearBuffer(alloc);
            }
            else {
                scaling = 0;
                for(unsigned short i = 0; i < Subdomain<K>::_map.size(); ++i)
                    for(unsigned int j = 0; j < Subdomain<K>::_map[i].second.size(); ++j)
                        _m[i][j] = 1.0 / (1.0 + array[Subdomain<K>::_map[i].second[j]].size());
            }
            delete [] array;
        }
        /* Function: apply
         *
         *  Applies the global FETI operator.
         *
         * Parameters:
         *    in             - Input vector.
         *    out            - Output vector (optional). */
        void apply(K* const* const in, K* const* const out = nullptr) const {
            A<'T', 0>(_primal, in);
            std::fill_n(super::_structure, super::_bi->_m, K());
            static_cast<Solver<K>*>(super::_pinv)->solve(super::_structure);
            A<'N', 0>(_primal, out ? out : in);
        }
        /* Function: applyLocalPreconditioner(n)
         *
         *  Applies the local preconditioner to multiple right-hand sides.
         *
         * Template Parameter:
         *    q              - Type of <FetiPrcndtnr> to apply.
         *
         * Parameters:
         *    u              - Input vectors.
         *    n              - Number of input vectors. */
        template<FetiPrcndtnr q = P>
        void applyLocalPreconditioner(K*& u, unsigned short n) const {
            switch(q) {
                case FetiPrcndtnr::DIRICHLET:   super::applyLocalSchurComplement(u, n); break;
                case FetiPrcndtnr::LUMPED:      super::applyLocalLumpedMatrix(u, n); break;
                case FetiPrcndtnr::SUPERLUMPED: super::applyLocalSuperlumpedMatrix(u, n); break;
                case FetiPrcndtnr::NONE:        break;
            }
        }
        /* Function: applyLocalPreconditioner
         *
         *  Applies the local preconditioner to a single right-hand side.
         *
         * Template Parameter:
         *    q              - Type of <FetiPrcndtnr> to apply.
         *
         * Parameter:
         *    u              - Input vector. */
        template<FetiPrcndtnr q = P>
        void applyLocalPreconditioner(K* const u) const {
            switch(q) {
                case FetiPrcndtnr::DIRICHLET:   super::applyLocalSchurComplement(u); break;
                case FetiPrcndtnr::LUMPED:      super::applyLocalLumpedMatrix(u); break;
                case FetiPrcndtnr::SUPERLUMPED: super::applyLocalSuperlumpedMatrix(u); break;
                case FetiPrcndtnr::NONE:        break;
            }
        }
        /* Function: precond
         *
         *  Applies the global preconditioner to a single right-hand side.
         *
         * Parameters:
         *    in             - Input vector.
         *    out            - Output vector (optional). */
        template<FetiPrcndtnr q = P>
        void precond(K* const* const in, K* const* const out = nullptr) const {
            A<'T', 1>(_primal, in);
            applyLocalPreconditioner<q>(_primal);
            A<'N', 1>(_primal, out ? out : in);
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
        void project(K* const* const in, K* const* const out = nullptr) const {
            static_assert(trans == 'T' || trans == 'N', "Unsupported value for argument 'trans'");
            if(super::_co) {
                if(!excluded) {
                    if(trans == 'T')
                        precond(in, _dual);
                    if(super::_ev) {
                        if(trans == 'T')
                            A<'T', 0>(_primal, _dual);
                        else
                            A<'T', 0>(_primal, in);
                        if(super::_schur) {
                            Blas<K>::gemv(&(Wrapper<K>::transc), &(Subdomain<K>::_dof), super::_co->getAddrLocal(), &(Wrapper<K>::d__1), *super::_ev, &(Subdomain<K>::_dof), _primal, &i__1, &(Wrapper<K>::d__0), super::_uc, &i__1);
                            super::_co->template callSolver<excluded>(super::_uc);
                            Blas<K>::gemv("N", &(Subdomain<K>::_dof), super::_co->getAddrLocal(), &(Wrapper<K>::d__1), *super::_ev, &(Subdomain<K>::_dof), super::_uc, &i__1, &(Wrapper<K>::d__0), _primal, &i__1);
                        }
                        else {
                            Blas<K>::gemv(&(Wrapper<K>::transc), &(Subdomain<K>::_dof), super::_co->getAddrLocal(), &(Wrapper<K>::d__1), *super::_ev + super::_bi->_m, &(Subdomain<K>::_a->_n), _primal, &i__1, &(Wrapper<K>::d__0), super::_uc, &i__1);
                            super::_co->template callSolver<excluded>(super::_uc);
                            Blas<K>::gemv("N", &(Subdomain<K>::_dof), super::_co->getAddrLocal(), &(Wrapper<K>::d__1), *super::_ev + super::_bi->_m, &(Subdomain<K>::_a->_n), super::_uc, &i__1, &(Wrapper<K>::d__0), _primal, &i__1);
                        }
                    }
                    else {
                        super::_co->template callSolver<excluded>(super::_uc);
                        std::fill_n(_primal, Subdomain<K>::_dof, K());
                    }
                    A<'N', 0>(_primal, _dual);
                    if(trans == 'N')
                        precond(_dual);
                    if(out)
                        for(unsigned int i = 0; i < super::_mult; ++i)
                            (*out)[i] = (*in)[i] - (*_dual)[i];
                    else
                        Blas<K>::axpy(&(super::_mult), &(Wrapper<K>::d__2), *_dual, &i__1, *in, &i__1);
                }
                else
                    super::_co->template callSolver<excluded>(super::_uc);
            }
            else if(!excluded && out)
                std::copy_n(*in, super::_mult, *out);
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
         * See also: <Bdd::buildTwo>, <Schwarz::buildTwo>.*/
        template<unsigned short excluded = 0>
        std::pair<MPI_Request, const K*>* buildTwo(const MPI_Comm& comm) {
            const Option& opt = *Option::get();
            if(!super::_schur && opt["geneo_nu"])
                super::_deficiency = opt["geneo_nu"];
            return super::template buildTwo<excluded, FetiProjection<decltype(*this), P, K>>(this, comm);
        }
        /* Function: computeSolution
         *
         *  Computes the solution after convergence of the Projected Conjugate Gradient.
         *
         * Template Parameter:
         *    excluded       - True if the master processes are excluded from the domain decomposition, false otherwise.
         *
         * Parameters:
         *    l              - Last iterate of the Lagrange multiplier.
         *    x              - Solution vector. */
        template<bool excluded>
        void computeSolution(K* const* const l, K* const x) const {
            if(!excluded) {
                A<'T', 0>(_primal, l);                                                                                                                                                                                                                   //    _primal = A^T l
                std::fill_n(super::_structure, super::_bi->_m, K());
                static_cast<Solver<K>*>(super::_pinv)->solve(super::_structure);                                                                                                                                                                         // _structure = S \ A^T l
                Blas<K>::axpy(&(Subdomain<K>::_a->_n), &(Wrapper<K>::d__2), super::_structure, &i__1, x, &i__1);                                                                                                                                         //          x = x - S \ A^T l
                if(super::_co) {
                    A<'N', 0>(x + super::_bi->_m, _dual);                                                                                                                                                                                                //      _dual = A (x - S \ A^T l)
                    precond(_dual);                                                                                                                                                                                                                      //      _dual = Q A (x - S \ A^T l)
                    if(!super::_ev)
                        super::_co->template callSolver<excluded>(super::_uc);
                    else {
                        A<'T', 0>(_primal, _dual);                                                                                                                                                                                                       //    _primal = A^T Q A (x - S \ A^T l)
                        if(super::_schur) {
                            Blas<K>::gemv(&(Wrapper<K>::transc), &(Subdomain<K>::_dof), super::_co->getAddrLocal(), &(Wrapper<K>::d__1), *super::_ev, &(Subdomain<K>::_dof), _primal, &i__1, &(Wrapper<K>::d__0), super::_uc, &i__1);                    //        _uc = R_b^T A^T Q A (x - S \ A^T l)
                            super::_co->template callSolver<excluded>(super::_uc);                                                                                                                                                                       //        _uc = (G Q G^T) \ R_b^T A^T Q A (x - S \ A^T l)
                            Blas<K>::gemv("N", &(Subdomain<K>::_dof), super::_co->getAddrLocal(), &(Wrapper<K>::d__1), *super::_ev, &(Subdomain<K>::_dof), super::_uc, &i__1, &(Wrapper<K>::d__0), _primal, &i__1);                                      //        x_b = x_b - R_b^T (G Q G^T) \ R_b^T A^T Q A (x - S \ A^T l)
                            Wrapper<K>::template csrmv<Wrapper<K>::I>(&(Wrapper<K>::transc), &(Subdomain<K>::_dof), &(super::_bi->_m), &(Wrapper<K>::d__2), false, super::_bi->_a, super::_bi->_ia, super::_bi->_ja, _primal, &(Wrapper<K>::d__0), super::_work);
                            if(super::_bi->_m)
                                super::_s.solve(super::_work);
                            Blas<K>::axpy(&(super::_bi->_m), &(Wrapper<K>::d__2), super::_work, &i__1, x, &i__1);
                            Blas<K>::axpy(&(Subdomain<K>::_dof), &(Wrapper<K>::d__2), _primal, &i__1, x + super::_bi->_m, &i__1);
                        }
                        else {
                            Blas<K>::gemv(&(Wrapper<K>::transc), &(Subdomain<K>::_dof), super::_co->getAddrLocal(), &(Wrapper<K>::d__1), *super::_ev + super::_bi->_m, &(Subdomain<K>::_a->_n), _primal, &i__1, &(Wrapper<K>::d__0), super::_uc, &i__1); //       _uc = R A^T Q A (x - S \ A^T l)
                            super::_co->template callSolver<excluded>(super::_uc);                                                                                                                                                                       //       _uc = (G Q G^T) \ R A^T Q A (x - S \ A^T l)
                            Blas<K>::gemv("N", &(Subdomain<K>::_a->_n), super::_co->getAddrLocal(), &(Wrapper<K>::d__2), *super::_ev, &(Subdomain<K>::_a->_n), super::_uc, &i__1, &(Wrapper<K>::d__1), x, &i__1);                                        //         x = x - R^T (G Q G^T) \ R A^T Q A (x - S \ A^T l)
                        }
                    }
                }
            }
            else if(super::_co)
                super::_co->template callSolver<excluded>(super::_uc);
        }
        template<bool>
        void computeSolution(const K* const, K* const) const { }
        /* Function: computeDot
         *
         *  Computes the dot product of two Lagrange multipliers.
         *
         * Template Parameter:
         *    excluded       - True if the master processes are excluded from the domain decomposition, false otherwise.
         *
         * Parameters:
         *    a              - Left-hand side.
         *    b              - Right-hand side. */
        template<bool excluded>
        void computeDot(underlying_type<K>* const val, const K* const* const a, const K* const* const b, const MPI_Comm& comm) const {
            if(!excluded)
                *val = std::real(Blas<K>::dot(&(super::_mult), *a, &i__1, *b, &i__1)) / 2.0;
            else
                *val = 0.0;
            MPI_Allreduce(MPI_IN_PLACE, val, 1, Wrapper<K>::mpi_underlying_type(), MPI_SUM, comm);
        }
        /* Function: getScaling
         *  Returns a constant pointer to <Feti::m>. */
        const underlying_type<K>* const* getScaling() const { return _m; }
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
        void solveGEVP(unsigned short& nu, const underlying_type<K>& threshold) {
            underlying_type<K>* const pt = reinterpret_cast<underlying_type<K>*>(_primal);
            for(unsigned short i = 0; i < Subdomain<K>::_map.size(); ++i)
                for(unsigned int j = 0; j < Subdomain<K>::_map[i].second.size(); ++j)
                    pt[Subdomain<K>::_map[i].second[j]] = _m[i][j];
            super::template solveGEVP<L>(pt, nu, threshold);
            nu = super::_deficiency;
            if(nu == 0 && super::_ev) {
                delete [] *super::_ev;
                delete []  super::_ev;
                super::_ev = nullptr;
            }
        }
};
} // HPDDM
#endif // _HPDDM_FETI_
