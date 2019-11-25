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
// SUMMARY : Conjugate gradient class
// LICENSE : LGPLv3
// ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
// AUTHORS : ...
// E-MAIL  : ...

// Definition of the conjugate gradient class
// Non-linear conjugate gradient algorithm
// author:  Wenceslau Gouveia
// modified:  H. Lydia Deng, 02/23/94, 03/14/94

// .NAME ConjugateGradient class
// .LIBRARY Base
// .HEADER Optimization Algorithms
// .INCLUDE defs.hpp
// .FILE CG.hpp

// .SECTION Description
// .B ConjugateGradient()
// The conjugate gradient procedure implemented in this object is an extension of
// the conjugate gradient used in linear system solving to handle non quadratic
// objective functions. Such extensions amount basically to the inclusion of a
// line search step and a modification in the computation of the conjugate directions.
// Details can be found in the classic (and tough!) book Practical Optimization by
// Powell.
//
// .SECTION Description
// Public Operations
// Constructors:
//		ConjugateGradient(LineSearch* ls, int iter, double tol)
// 		Here:
//			ls: Defines the line search used. At the present version you should use the
//			    CubicLineSearch procedure.
//			iter: Maximum number of iterations
//		 	tol: Minimum accepted module of the gradient at optimum solution
// Methods:
//		Model<>ConjugateGradient::optimizer(Model<double>&)
//		Here:
//			model0:  Initial model for the conjugate gradient procedure.
//
//			The optimum model is returned by the function.
//
// .SECTION Caveats
// This procedure requires the derivatives of the objective function. At the present version the
// COOOL library cannot handle analytic derivatives provided by the user. Here the derivatives
// are computed numerically by a centered finite differencing scheme in a automatic fashion.
// The hability to support user-provided derivatives will hopefully be
// implemented in the near future.

#ifndef CONJUGATE_GRADIENT_HH
#define CONJUGATE_GRADIENT_HH

#include "Optima.hpp"
#include "defs.hpp"
#include "Param.hpp"

// Nonlinear Conjugate Gradient
/*
  This conjugate gradient procedure implemented in this
  object is an extension of the conjugate gradient
  used in linear system solving to handle non quadratic
  objective functions. Such extensions amount basically
  to the inclusion of a line search step and a
  modification in the computation of the conjugate directions.
  Details can be found in the classic (and tough!) book
  Practical Optimization by Powell.
*/

template< class LS >
class ConjugateGradient : public Optima< LS > {
  typedef typename LS::Real Real;
  typedef typename LS::Param Param;
  typedef typename LS::Vect Vect;
  typedef typename LS::VMat VMat;
  typedef LS LineSearch;

 public:
  // a constructor
  ConjugateGradient(LS* ls,         // Pointer to the line-search object
                    int iter,       // Maximum number of iterations
                    Real tol,       // Minimum accepted gradient at optimum solution
                    int verb = 0    // verbose or quiet
  );
  /* At the present version, you must use the CubicLineSearch procedure */
  // a constructor
  ~ConjugateGradient( ) { ; }

  // Conjugate gradient search starting from m0, returns an optimum Model
  Param optimizer(Param& m0);
};

template< class LS >
ConjugateGradient< LS >::ConjugateGradient(LS* p, int it, Real eps, int verb) : Optima< LS >(verb) {
  ls = p;
  iterMax = it;
  tol = eps;
  iterNum = 0;
}

// Needs :
// A scalar product on Vect: (u,v)
// A Real Vect product: operator*
// A Real Vect product: operator*=
// A Vect difference: operator-
template< class LS >
ConjugateGradient< LS >::Param ConjugateGradient< LS >::optimizer(Param& model0) {
  // Reset the residue history for every new optimizer
  iterNum = 0;
  isSuccess = 0;
  if (residue != NULL) {
    delete residue;
    residue = new list< Real >;
  }

  int n = model0.size( );
  Param model1(model0);    // new model
  Vect search(n);          // search direction
  Vect g0(n);              // old gradient vector
  Vect g1(n);              // new gradient vector
  double beta;             // beta parameter
  double lambda = .025;    // line search parameter
  double descent = 0.;     // descent direction

  // Beginning iterations
  g0 = *(ls->gradient(model0));

  // Check the gradient, in case the initial model is the optimal
  Real err = (Real)sqrt((g0, g0));
  if (isVerbose) cout << "Initial residue : " << err << endl;
  appendResidue(err);    // residual
  if (err < tol) {
    if (isVerbose) cout << "Initial guess was great!" << endl;
    isSuccess = 1;
    return model0;
  }

  // Considering first iteration
  search = -1. * g0;
  descent = (search, g0);

  // We use CubicLineSearch
  //    cerr << "Line Search" << endl;
  //    cerr << "model0 " << model0;
  //    cerr << "search " << (Param)search; // we cast only for display
  //    cerr << "descent " << descent << endl;
  //    cerr << "lambda " << lambda << endl;
  //    cerr << endl;

  model1 = ls->search(model0, search, descent, lambda);
  g1 = *(ls->gradient(model1));    // Gradient at new model
  err = (Real)sqrt((g1, g1));
  if (isVerbose)
    cout << "Iteration (0) : "
         << "current value of the objective function: " << ls->currentValue( )
         << "\t current residue: " << err << endl;
  appendResidue(err);    // residual

  iterNum = 0;
  Real temp;
  do {
    iterNum++;

    temp = 1. / ((g0, g0));
    beta = ((g1 - g0), g1);
    beta *= temp;    // computation Polak & Ribiere

    search = beta * search - g1;    // search direction

    descent = (search, g1);    // descent
    if (descent > 0.) {
      if (isVerbose) cout << "Reset search direction to gradient vector!" << endl;
      search = -1. * g1;
      descent = (search, g1);
    }

    model0 = model1;
    g0 = g1;    // save the old model and gradient before new search

    // We use CubicLineSearch
    //  	  cerr << "Line Search" << endl;
    //  	  cerr << "model0 " << model0;
    //  	  cerr << "search " << (Param)search;
    //  	  cerr << "descent " << descent << endl;
    //  	  cerr << "lambda " << lambda << endl;
    //  	  cerr << endl;

    model1 = ls->search(model0, search, descent, lambda);    // line search
    g1 = *(ls->gradient(model1));

    err = (Real)sqrt((g1, g1));
    if (isVerbose)
      cout << "Iteration (" << iterNum << ") : "
           << "current value of the nrj : " << ls->currentValue( ) << "\t current residue : " << err
           << endl;
    appendResidue(err);    // residual

  } while (finalResidue( ) > tol && iterNum < iterMax);    // stopping criterion

  if (finalResidue( ) <= tol) isSuccess = 1;

  return (model1);    // hopefully answer
}

#endif
