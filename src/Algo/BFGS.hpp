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
// SUMMARY : BFGS Quasi-Newton Optimization Method
// LICENSE : LGPLv3
// ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
// AUTHORS : ...
// E-MAIL  : ...

/* At the present version, you must use the CubicLineSearch procedure */

#ifndef _BFGS_HH_
#define _BFGS_HH_

#include "Optima.hpp"
#include "defs.hpp"

// Need a "true" Matrix class
template< class LS >
class BFGS : public Optima< LS > {
  typedef typename LS::Real Real;
  typedef typename LS::Param Param;
  typedef typename LS::Vect Vect;
  typedef typename LS::VMat VMat;
  typedef typename LS::Mat Mat;
  typedef list< Real > mlist;

 public:
  BFGS(LS* ls,         // pointer to the line-search object
       int iter,       // Maximum number of iterations
       Real tol,       // Minimum accepted gradient at optimum solution
       int verb = 0    // vebose or quiet
  );
  ~BFGS( ) { ; }

  // the BFGS search starting from model0, returns an optimum model
  Param optimizer(Param& model0) override;
};

template< class LS >
BFGS< LS >::BFGS(LS* p, int it, Real eps, int verb) : Optima< LS >(verb) {
  this->ls = p;
  this->iterMax = it;
  this->tol = eps;
  this->iterNum = 0;
}

// Matrix class requirements:
// - A matrix constructor with n_rows, n_cols
// - A matrix constructor from diagonal
// - A matrix-vector product
// - An outProduct method: from two vector u, v, build the matrix A_{i,j}=u(i)v(j)
// - A t() method
// - An operator +
// - An operator = Real
// - An operator = Vect
template< class LS >
typename BFGS< LS >::Param BFGS< LS >::optimizer(Param& model0) {
  // reset the residue history for every new optimizer
  this->iterNum = 0;
  if (this->residue != NULL) {
    delete this->residue;
    this->residue = new mlist;
  }

  // Initial settings for some parameters
  int n = model0.size( );
  Vect g0(n);
  double lambda = 0.025;
  double descent = 0.;

  g0 = *(this->ls->gradient(model0));

  // Check the gradient, in case the initial model is the optimal
  double err = (Real)sqrt((g0, g0));

  if (this->isVerbose) cerr << "Initial residue : " << err << endl;

  this->appendResidue(err);    // residual

  if (err < this->tol) {
    if (this->isVerbose) cerr << "Initial guess was great! \n";
    this->isSuccess = 1;
    return model0;
  }

  // Initial identical matrix for estimating inverse of the Hessian
  // Vect diag(n, 1.);

  Mat H(n, n);
  H = 0;
  diagonal(H) = 1;
  Real d, dd, scale;
  Param model1(model0);
  Vect s(n), gamma(n), delta(n), g1(n);

  // Searching directions
  s = Real( );
  s -= H * g0;
  descent = (s, g0);
  assert(s.max( ) < 1e100);

  // Cubic line search for a new model
  model1 = this->ls->search(model0, s, descent, lambda);
  g1 = *(this->ls->gradient(model1));
  err = (Real)sqrt((g1, g1));
  if (this->isVerbose)
    cerr << "Iteration (" << this->iterNum << ") : "
         << "current value of the objective function: " << this->ls->currentValue( )
         << "\t current residue: " << err << endl;

  this->appendResidue(err);    // residual
  this->iterNum++;

  Mat B(n, n);

  while (this->finalResidue( ) > this->tol && this->iterNum < this->iterMax) {
    gamma = g1 - g0;
    delta = model1 - model0;

    // Replace the searching direction with temporal storage
    s = H * gamma;

    // Factor of the denominator
    dd = (delta, gamma);

    // check denominator
    if (Abs(dd) < 1e-20) {
      // Re-initialize the Hessian Matrix
      // It must be zeroed first (cf. Matrix.hpp)
      H = 0.;
      diagonal(H) = 1.;
    } else {
      assert(dd);
      d = 1. / dd;

      scale = d * ((gamma, s));
      scale += 1;
      scale *= d;

      // Update the first term
      H += scale * delta * delta.t( );

      // Update the second term
      H -= d * s * delta.t( );
      H -= d * delta * s.t( );

      // Store the current model and gradient
      g0 = g1;
      model0 = model1;
    }

    s = Real( );
    s -= H * g0;
    descent = (s, g0);
    model1 = this->ls->search(model0, s, descent, lambda);
    g1 = *(this->ls->gradient(model1));
    err = (Real)sqrt((g1, g1));

    if (this->isVerbose)
      cerr << "Iteration (" << this->iterNum << ") : "
           << "current value of the objective function: " << this->ls->currentValue( )
           << "\t current residue: " << err << endl;

    this->appendResidue(err);    // residual
    this->iterNum++;
  }

  return (model1);
}

#endif    //_BFGS_HH_
