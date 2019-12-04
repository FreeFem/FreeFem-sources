// Tony : lambda=0.25 ou lambda=0.025 ???
// cf. CG.hpp
//
//			The minimizer along the search direction is returned by the function.
//
/*
  CubicLineSearch() class implements the efficient line search procedure
  described in Dennis and Schbabel's book entitled "Numerical
  Methods for Unconstrained and Nonlinear Equations. The
  objective is to perform a unidimensional search for a minimum
  point along a specified direction in a multidimensional
  space.
*/
// This procedure seems to be fairly robust. It has worked for a fairly broad class of
// problems from optimization of standard test functons in optimization theory and
// to hard geophysical problems as stacking power optimization and amplitude
// seismogram inversion

#ifndef CUBIC_LINE_SEARCH_HH
#define CUBIC_LINE_SEARCH_HH

#include "LineSearch.hpp"
#include <limits.h>
#include <float.h>

template< class LS >
class CubicLineSearch : public LS {
  typedef typename LS::Real Real;
  typedef typename LS::Param Param;
  typedef typename LS::Vect Vect;
  typedef typename LS::VMat VMat;
  typedef LS LineSearch;
  typedef tNRJ< Param, Vect, VMat, Real > NRJ;

 public:
  // a constructor with the default delta
  CubicLineSearch(NRJ* f, int iter);
  // a constructor with the specified delta
  CubicLineSearch(NRJ* f, int iter, Vect* delta);
  /*    The parameter $delta$ is not used by the line search
        itself. Rather it is used in the numerical computation
        of the derivatives using centered differences. For
        example the derivative of f(x) at the point x0 would be
        given by
        \[f(x0 - delta) - f(x0 + delta) / 2 * delta\]
        */
  ~CubicLineSearch( );

  // search for the minimum model along a 1-D direction
  Param search(    /// initial model for the line search
    const Param& m0,
    /// search direction
    Vect& direction,
    /// product of search direction and gradient
    Real descent,
    /// a parameter
    double lambda);
  /* The parameter $lambda$ controls the accuraccy of the
    line search. $lambda = .25$ is a good choice.
   The minimizer along the search direction is returned
   by the function.  */

  // lambda=0.25 ou lambda=0.025 ???
  // cf. CG.hpp
};

template< class LS >
CubicLineSearch< LS >::CubicLineSearch(NRJ* f, int it) : LS(f) {
  this->iterMax = it;
}

template< class LS >
CubicLineSearch< LS >::CubicLineSearch(NRJ* f, int it, Vect* interval) : LS(f, interval) {
  this->iterMax = it;
}

template< class LS >
CubicLineSearch< LS >::~CubicLineSearch( ) {}

// Code for the Cubic Line Search
template< class LS >
typename CubicLineSearch< LS >::Param CubicLineSearch< LS >::search(const Param& current_solution,
                                                                    Vect& p, Real slope,
                                                                    double lambda) {

  int tst = 0;    // useful vars
  Real alpha = 0., alpha_prev = 0;
  Real fprev = 0;
  Real new_m = 0, old_m = 0;
  Param new_solution(current_solution);
  cout << " search " << p.max( ) << endl;
  assert(p.max( ) < 1e100);
  old_m = this->nrj->getVal(current_solution);    // Evaluation at the current solution
  this->iterNum = 0;
  this->iterNum++;    // iteration counter
  alpha = 1.;         // updating step

  new_solution = this->update(current_solution, 1, alpha, p);    // updating
  new_m = this->nrj->getVal(new_solution);                       // Evaluation at the
                                                                 // new solution
  this->iterNum++;

  // Implementing Goldstein's test for alpha too small
  while (new_m < old_m + (1. - lambda) * alpha * slope && this->iterNum < this->iterMax) {
    alpha *= 3;
    new_solution = this->update(current_solution, 1, alpha, p);
    new_m = this->nrj->getVal(new_solution);
    this->iterNum++;
  }
  if (this->iterNum == this->iterMax) cerr << "Alpha overflowed! \n";

  // Armijo's test for alpha too large
  alpha_prev = alpha;    // H.L. Deng, 6/13/95
  while (new_m > old_m + lambda * alpha * slope && this->iterNum < this->iterMax) {
    Real alpha2 = alpha * alpha, alpha_tmp = 0.;
    Real f1 = new_m - old_m - slope * alpha;
    Real a = 0, b = 0;
    Real c = 0, cm11 = 0, cm12 = 0, cm21 = 0, cm22 = 0;
    Real disc = 0;

    if (tst == 0) {
      alpha_tmp = -slope * alpha2 / (f1 * 2.);
      // tentative alpha

      tst = 1;
    } else {
      Real alpha_prev2 = alpha_prev * alpha_prev;
      Real f2 = fprev - old_m - alpha_prev * slope;

      c = 1. / (alpha - alpha_prev);
      cm11 = 1. / alpha2;
      cm12 = -1. / alpha_prev2;
      cm21 = -alpha_prev / alpha2;
      cm22 = alpha / alpha_prev2;

      a = c * (cm11 * f1 + cm12 * f2);
      b = c * (cm21 * f1 + cm22 * f2);
      disc = b * b - 3. * a * slope;

      if ((Abs(a) > FLT_MIN) && (disc > FLT_MIN))
        alpha_tmp = (-b + sqrt(disc)) / (3. * a);
      else
        alpha_tmp = slope * alpha2 / (2. * f1);

      if (alpha_tmp >= .5 * alpha) alpha_tmp = .5 * alpha;
    }
    alpha_prev = alpha;
    fprev = new_m;

    if (alpha_tmp < .1 * alpha)
      alpha *= .1;
    else
      alpha = alpha_tmp;

    new_solution = this->update(current_solution, 1, alpha, p);
    new_m = this->nrj->getVal(new_solution);
    this->iterNum++;
  }
  if (this->iterNum == this->iterMax) {
    cerr << "Alpha underflowed! \n";
    cerr << this->iterMax;
  }

  this->value = new_m;
  this->appendSearchNumber( );
  return (new_solution);    // # of iterations
}

#endif
