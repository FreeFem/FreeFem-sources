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

#ifndef _HPDDM_PRECONDITIONER_
#define _HPDDM_PRECONDITIONER_

#define HPDDM_LAMBDA_F(in, input, inout, output, N)                                                          \
    unsigned short* input = static_cast<unsigned short*>(in);                                                \
    unsigned short* output = static_cast<unsigned short*>(inout);                                            \
    output[0] = std::max(output[0], input[0]);                                                               \
    output[1] = output[1] & input[1];                                                                        \
    output[2] = output[2] & input[2];                                                                        \
    if(N == 3)                                                                                               \
        output[3] = output[3] & input[3];

#include "subdomain.hpp"
#include "coarse_operator_impl.hpp"
#include "operator.hpp"

namespace HPDDM {
/* Class: Preconditioner
 *
 *  A base class from which <Schwarz> and <Schur> inherit.
 *
 * Template Parameters:
 *    Solver         - Solver used for the factorization of local matrices.
 *    CoarseOperator - Class of the coarse operator.
 *    K              - Scalar type. */
template<template<class> class Solver, class CoarseOperator, class K>
class Preconditioner : public Subdomain<K> {
#ifdef __MINGW32__
    private:
        template<unsigned short N>
        static void __stdcall f(void* in, void* inout, int*, MPI_Datatype*) {
            HPDDM_LAMBDA_F(in, input, inout, output, N)
        }
#endif
    protected:
        /* Variable: s
         *  Solver used in <Schwarz::callNumfact> and <Schur::callNumfactPreconditioner> or <Schur::computeSchurComplement>. */
        Solver<K>           _s;
        /* Variable: co
         *  Pointer to a <Coarse operator>. */
        CoarseOperator*    _co;
        /* Variable: ev
         *  Array of deflation vectors as needed by <Preconditioner::co>. */
        K**                _ev;
        /* Variable: uc
         *  Workspace array of size <Coarse operator::local>. */
        K*                 _uc;
        /* Function: buildTwo
         *
         *  Assembles and factorizes the coarse operator.
         *
         * Template Parameter:
         *    excluded       - Greater than 0 if the master processes are excluded from the domain decomposition, equal to 0 otherwise.
         *
         * Parameters:
         *    A              - Operator used in the definition of the Galerkin matrix.
         *    comm           - Global MPI communicator. */
        template<unsigned short excluded, class Operator, class Prcndtnr>
        std::pair<MPI_Request, const K*>* buildTwo(Prcndtnr* B, const MPI_Comm& comm) {
            static_assert(std::is_same<typename Prcndtnr::super&, decltype(*this)>::value || std::is_same<typename Prcndtnr::super::super&, decltype(*this)>::value, "Wrong preconditioner");
            std::pair<MPI_Request, const K*>* ret = nullptr;
            constexpr unsigned short N = std::is_same<typename Prcndtnr::super&, decltype(*this)>::value ? 2 : 3;
            unsigned short allUniform[N + 1];
            allUniform[0] = Subdomain<K>::_map.size();
            const Option& opt = *Option::get();
            unsigned short nu = allUniform[1] = (_co ? _co->getLocal() : static_cast<unsigned short>(opt["geneo_nu"]));
            allUniform[2] = static_cast<unsigned short>(~nu);
            if(N == 3)
                allUniform[3] = nu > 0 ? nu : std::numeric_limits<unsigned short>::max();
            {
                MPI_Op op;
#ifdef __MINGW32__
                MPI_Op_create(&f<N>, 1, &op);
#else
                auto f = [](void* in, void* inout, int*, MPI_Datatype*) -> void {
                    HPDDM_LAMBDA_F(in, input, inout, output, N)
                };
                MPI_Op_create(f, 1, &op);
#endif
                MPI_Allreduce(MPI_IN_PLACE, allUniform, N + 1, MPI_UNSIGNED_SHORT, op, comm);
                MPI_Op_free(&op);
            }
            if(nu > 0 || allUniform[1] != 0 || allUniform[2] != std::numeric_limits<unsigned short>::max()) {
                if(!_co) {
                    _co = new CoarseOperator;
                    _co->setLocal(nu);
                }
                double construction = MPI_Wtime();
                if(allUniform[1] == nu && allUniform[2] == static_cast<unsigned short>(~nu))
                    ret = _co->template construction<1, excluded>(std::move(Operator(*B, allUniform[0])), comm);
                else if(N == 3 && allUniform[1] == 0 && allUniform[2] == static_cast<unsigned short>(~allUniform[3]))
                    ret = _co->template construction<2, excluded>(std::move(Operator(*B, allUniform[0])), comm);
                else
                    ret = _co->template construction<0, excluded>(std::move(Operator(*B, allUniform[0])), comm);
                construction = MPI_Wtime() - construction;
                if(_co->getRank() == 0 && opt.val<int>("verbosity") > 0) {
                    std::stringstream ss;
                    ss << std::setprecision(2) << construction;
                    std::string line = " --- coarse operator transferred and factorized by " + to_string(static_cast<int>(opt["master_p"])) + " process" + (static_cast<int>(opt["master_p"]) == 1 ? "" : "es") + " (in " + ss.str() + "s)";
                    std::cout << line << std::endl;
                    std::cout << std::right << std::setw(line.size()) << "(criterion = " + to_string(allUniform[1] == nu && allUniform[2] == static_cast<unsigned short>(~nu) ? nu : (N == 3 && allUniform[2] == static_cast<unsigned short>(~allUniform[3]) ? -_co->getLocal() : 0)) + " -- topology = " + to_string(static_cast<int>(opt["master_topology"])) + " -- distribution = " + to_string(static_cast<int>(opt["master_distribution"])) + ")" << std::endl;
                }
                _uc = new K[_co->getSizeRHS()];
            }
            else {
                delete _co;
                _co = nullptr;
            }
            return ret;
        }
    public:
        Preconditioner() : _co(), _ev(), _uc() { }
        Preconditioner(const Preconditioner&) = delete;
        ~Preconditioner() {
            delete _co;
            if(_ev && *_ev)
                delete [] *_ev;
            delete [] _ev;
            delete [] _uc;
        }
        /* Typedef: super
         *  Type of the immediate parent class <Subdomain>. */
        typedef Subdomain<K> super;
        /* Function: initialize
         *
         *  Initializes a two-level preconditioner.
         *
         * Parameter:
         *    deflation      - Number of local deflation vectors. */
        void initialize(const unsigned short& deflation) {
            if(!_co) {
                _co = new CoarseOperator;
                _co->setLocal(deflation);
            }
        }
        /* Function: callSolve
         *
         *  Applies <Preconditioner::s> to multiple right-hand sides in-place.
         *
         * Parameters:
         *    x              - Input right-hand sides, solution vectors are stored in-place.
         *    n              - Number of input right-hand sides. */
        void callSolve(K* const x, const unsigned short& n = 1) const { _s.solve(x, n); }
        /* Function: getVectors
         *  Returns a constant pointer to <Preconditioner::ev>. */
        K** getVectors() const { return _ev; }
        /* Function: setVectors
         *  Sets the pointer <Preconditioner::ev>. */
        void setVectors(K** const& ev) { _ev = ev; }
        /* Function: getLocal
         *  Returns the value of <Coarse operator::local>. */
        unsigned short getLocal() const { return _co ? _co->getLocal() : 0; }
        /* Function: getAddrLocal
         *  Returns the address of <Coarse operator::local> or <i__0> if <Preconditioner::co> is not allocated. */
        const int* getAddrLocal() const { return _co ? _co->getAddrLocal() : &i__0; }
};
} // HPDDM
#endif // _HPDDM_PRECONDITIONER_
