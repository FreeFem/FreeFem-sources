/****************************************************************************/
/* This file is part of FreeFEM.                                            */
/*                                                                          */
/* FreeFEM is free software: you can redistribute it and/or modify          */
/* it under the terms of the GNU Lesser General Public License as           */
/* published by the Free Software Foundation, either version 3 of           */
/* the License, or (at your option) any later version.                      */
/*                                                                          */
/* FreeFEM is distributed in the hope that it will be useful,               */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of           */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            */
/* GNU Lesser General Public License for more details.                      */
/*                                                                          */
/* You should have received a copy of the GNU Lesser General Public License */
/* along with FreeFEM. If not, see <http://www.gnu.org/licenses/>.          */
/****************************************************************************/
// SUMMARY : ...
// LICENSE : LGPLv3
// ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
// AUTHORS : ...
// E-MAIL  : ...

#ifndef NEWTONRAPHSON_HH
#define NEWTONRAPHSON_HH

#include "Optima.hpp"
#include "defs.hpp"

/* At the present version, you must use the CubicLS procedure */

// Need a "real" Matrix classe
template< class LS >
class Newt : public Optima< LS > {
  typedef typename LS::Real Real;
  typedef typename LS::Param Param;
  typedef typename LS::Vect Vect;
  typedef typename LS::VMat VMat;
  typedef typename LS::Mat Mat;
  typedef list< Real > mlist;

 public:
  Newt(LS* ls,           // pointer to the line-search object.
       int iter,         // Maximum number of iterations
       Real tol,         // minimum accepted gradient at optimum solution
       int verb = 0);    // verbose or quiet
  ~Newt( ) { ; }

  // the Newt search starting from model0, returns an optimum model
  Param optimizer(Param& model0) override;
};

template< class LS >
Newt< LS >::Newt(LS* p, int it, Real eps, int verb) : Optima< LS >(verb) {
  this->ls = p;
  this->iterMax = it;
  this->tol = eps;
  this->iterNum = 0;
}

// MAT class needs:
// - Ax=b resolution method
template< class LS >
typename Newt< LS >::Param Newt< LS >::optimizer(Param& model0) {
  // reset the this->residue history for every new optimizer
  this->iterNum = 0;
  if (this->residue != NULL) {
    delete this->residue;
    this->residue = new mlist;
  }

  // Initial settings for some parameters
  int n = model0.size( );
  Vect g0(n);
  double lambda = 0.025;

  g0 = *(this->ls->gradient(model0));

  // check the gradient, in case the initial model is the optimal
  double err = (Real)sqrt((g0, g0));

  if (this->isVerbose) cerr << "Initial this->residue : " << err << endl;

  this->appendResidue(err);    // residual

  if (err < this->tol) {
    if (this->isVerbose) cerr << "Initial guess was great!" << endl;
    this->isSuccess = 1;
    return model0;
  }

  Vect g1(n);
  Vect s(n);
  Param model1(model0);

  while (this->finalResidue( ) > this->tol && this->iterNum < this->iterMax) {
    // searching directions
    // s = g0/(*this->ls->hessian(model0)); // LU reinitialization at each step
    this->ls->hessian(model0)->Solve(s, g0);
    s = -1. * s;

    double descent = (s, g0);
    // Cubic Line Search
    model1 = this->ls->search(model0, s, descent, lambda);
    g1 = *(this->ls->gradient(model1));
    err = (Real)sqrt((g1, g1));

    if (this->isVerbose)
      cerr << "Iteration (" << this->iterNum << ") : "
           << "current value of the objective function: " << this->ls->currentValue( )
           << "\t current this->residue: " << err << endl;

    this->appendResidue(err);    // residual

    g0 = g1;
    model0 = model1;

    this->iterNum++;
  }

  return (model1);
}

#endif
