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

#ifndef LINE_SEARCH_HH
#define LINE_SEARCH_HH

#include <list>
#include "defs.hpp"
#include "NRJ.hpp"
#include <limits.h>
#include <float.h>

#define MAX_IT_LS 1000

/*
  This is the base class for one dimensional
  line search algorithms. Some optimization
  algorithms need to find the optimimal point
  on one searching directions.
  Classes derived from this class inherit features of
  LineSearch. This class cannot be used directly.
*/

// Warning:
// no copy operator
// no = operator

template< class P, class V, class M, class VM, class R >
class LineSearch {
 public:
  typedef R Real;
  typedef VM VMat;
  typedef M Mat;
  typedef P Param;
  typedef V Vect;
  typedef tNRJ< Param, Vect, VMat, Real > NRJ;

 private:
  // Used only if we calculate numerically the gradient
  // in this case, we allocate memory
  // else, this is only a pointer to the NRJ gradient
  Vect *grad;
  int ref;    // ref equals 1 if we calcualte numerically the gradient

  void numericalGradient(const Param &);

 protected:
  // maxinum number of iterations allowed
  int iterMax;
  // the number iterations so far
  // to reinitialize in inherited classes
  int iterNum;
  // the history of number of iteration of iteration
  list< int > *iterHistory;
  // the value at the minimum along the direction
  // to update in inherited classes
  Real value;
  // the delta used for calculating numerical gradient
  // same type as gradient
  Vect *step;
  // pointer to the nrj
  NRJ *nrj;

  void appendSearchNumber( );

 public:
  // a constructor with pointer to the nrj and
  // to the step if numerical gradient is computed
  LineSearch(NRJ *f, Vect *interval = NULL);
  virtual ~LineSearch( );

  // Start point and descent direction + parameters
  virtual Param search(const Param &, Vect &, Real, double);

  // compute the gradient Vector
  // linesearch gradient
  Vect *gradient(const Param &m);
  // linesearch hessian
  VMat *hessian(const Param &m);
  // evaluate the nrj
  Real evaluate(const Param &m);

  // number of iterations for all line searches
  list< int > allSearchIterations( );

  // number of iterations in the current line search
  int searchIterations( );
  // value of the objective function for the current Model
  Real currentValue( );

  Param update(const Param &, Real, Real, const Vect &) const;
};

template< class P, class V, class M, class VM, class R >
LineSearch< P, V, M, VM, R >::LineSearch(NRJ *p, Vect *interval) {
  iterMax = MAX_IT_LS;
  iterNum = 0;
  nrj = p;
  step = interval;
  iterHistory = new list< int >;
  grad = NULL;
  ref = 0;
  value = 0.;
}

template< class P, class V, class M, class VM, class R >
LineSearch< P, V, M, VM, R >::~LineSearch( ) {
  if (iterHistory != NULL) delete iterHistory;
  if (ref == 1 && grad != NULL) delete grad;
  if (step != NULL) delete step;
}

template< class P, class V, class M, class VM, class R >
void LineSearch< P, V, M, VM, R >::appendSearchNumber( ) {
  iterHistory[0].push_back(iterNum);
}

template< class P, class V, class M, class VM, class R >
list< int > LineSearch< P, V, M, VM, R >::allSearchIterations( ) {
  return iterHistory[0];
}

template< class P, class V, class M, class VM, class R >
int LineSearch< P, V, M, VM, R >::searchIterations( ) {
  return iterNum;
}

template< class P, class V, class M, class VM, class R >
typename LineSearch< P, V, M, VM, R >::Real LineSearch< P, V, M, VM, R >::currentValue( ) {
  return value;
}

template< class P, class V, class M, class VM, class R >
typename LineSearch< P, V, M, VM, R >::Param LineSearch< P, V, M, VM, R >::search(const Param &m,
                                                                                  Vect &v,
                                                                                  Real alpha,
                                                                                  double beta) {
  cerr << "You need to specify the LineSearch method!" << endl;
  exit(1);
  return 0;
}

// Centered finite difference
template< class P, class V, class M, class VM, class R >
void LineSearch< P, V, M, VM, R >::numericalGradient(const Param &m) {
  int n = m.size( );

  // At first call: gradient initialization
  if (grad == NULL) {
    grad = new Vect(n);
    ref = 1;
  }

  Real diffValue, newValue, off;
  Param m1(m);

  if (step != NULL) {
    for (int i = 0; i < n; i++) {
      off = step[0][i];
      m1[i] -= off;
      newValue = nrj->getVal(m1);
      m1[i] += 2 * off;
      diffValue = nrj->getVal(m1);
      diffValue -= newValue;
      (*grad)[i] = diffValue / (2 * off);
      m1[i] = m[i];
    }
  } else {    // Default: step of relative size 10^-3
    for (int i = 0; i < n; i++) {
      off = m[i] / 1000;
      m1[i] -= off;
      newValue = nrj->getVal(m1);
      m1[i] += 2 * off;
      diffValue = nrj->getVal(m1);
      diffValue -= newValue;
      (*grad)[i] = diffValue / (2 * off);
      m1[i] = m[i];
    }
  }
}

template< class P, class V, class M, class VM, class R >
typename LineSearch< P, V, M, VM, R >::Vect *LineSearch< P, V, M, VM, R >::gradient(
  const Param &m) {
  Vect *g;

  g = nrj->getGradient(m);

  if (g == NULL) {
    cerr << "Undefined gradient for this NRJ function!" << endl;
    cerr << "will use a numerical gradient..." << endl;
    numericalGradient(m);
  } else
    grad = g;

  return grad;
}

template< class P, class V, class M, class VM, class R >
typename LineSearch< P, V, M, VM, R >::VMat *LineSearch< P, V, M, VM, R >::hessian(const Param &m) {
  VMat *h;

  h = nrj->getHessian(m);

  if (h == NULL) {
    cerr << "Hessien is undefined for this NRJ function!" << endl;
    exit(1);
  }

  return h;
}

template< class P, class V, class M, class VM, class R >
typename LineSearch< P, V, M, VM, R >::Real LineSearch< P, V, M, VM, R >::evaluate(const Param &m) {
  return nrj->getVal(m);
}

// We moving forward from alpha*m, m is the start point in the direction dmod, of beta
// We check we are in the bounds
template< class P, class V, class M, class VM, class R >
typename LineSearch< P, V, M, VM, R >::Param LineSearch< P, V, M, VM, R >::update(
  const Param &m, Real alpha, Real beta, const Vect &dmod) const {
  Param newparam(m);
  long ndim = m.size( );

  Vect direction(dmod);

  if (direction.size( ) != ndim) {
    cerr << "Size of direction is different from dimensions" << endl;
    exit(1);
  }

  // check to see if the new model is out of the upper bound
  cout << "update " << m.min( ) << " " << m.max( ) << " " << direction.max( ) << " " << dmod.max( );
  direction = alpha * m + beta * direction;
  newparam = direction;
  cout << " ; " << newparam.min( ) << " " << newparam.max( ) << endl;
  if ((m.modMax( ) != NULL) && (m.modMin( ) != NULL)) {
    for (long i = 0; i < ndim; i++) {
      if ((newparam[i] - (*(m.modMax( )))[i]) > 0) newparam[i] = (*(m.modMax( )))[i];

      if ((newparam[i] - (*(m.modMin( )))[i]) < 0) newparam[i] = (*(m.modMin( )))[i];
    }
  }

  return newparam;
}

#endif
