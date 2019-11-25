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

#ifndef NRJ_HH
#define NRJ_HH

#include "defs.hpp"

// Warning :
// - no copy operator
// - no = operator

template< class Param, class Vect, class Mat, class Real >
class tNRJ {
 protected:
  int nappel_val;
  int nappel_grad;
  int nappel_hess;
  int nparam;

  Real val;
  Vect *grad;
  Mat *hess;

 public:
  // Warning:
  // In the tNRJ constructor, you must initialize grad and/or hess if you use it
  // That is not done automatically as it depends on type
  tNRJ(int);
  virtual ~tNRJ( );

  // To precise for each function of tNRJ
  virtual Real Val(const Param &) = 0;

  Real getVal(const Param &);

  // To precise if needed, else return null
  virtual Vect *Gradient(const Param &);

  Vect *getGradient(const Param &);

  // To precise if needed, else return null
  virtual Mat *Hessian(const Param &);

  Mat *getHessian(const Param &);

  int Appel_Val( ) const { return nappel_val; };
  int Appel_Grad( ) const { return nappel_grad; };
  int Appel_Hess( ) const { return nappel_hess; };
};

template< class Param, class Vect, class Mat, class Real >
tNRJ< Param, Vect, Mat, Real >::~tNRJ( ) {
  if (grad != NULL) delete grad;
  if (hess != NULL) delete hess;
}

template< class Param, class Vect, class Mat, class Real >
tNRJ< Param, Vect, Mat, Real >::tNRJ(int n) {
  nparam = n;
  nappel_val = 0;
  nappel_grad = 0;
  nappel_hess = 0;
  grad = NULL;
  hess = NULL;
  val = 0.;
}

template< class Param, class Vect, class Mat, class Real >
Real tNRJ< Param, Vect, Mat, Real >::getVal(const Param &p) {
  nappel_val++;
  return Val(p);
}

template< class Param, class Vect, class Mat, class Real >
Vect *tNRJ< Param, Vect, Mat, Real >::getGradient(const Param &p) {
  nappel_grad++;
  return Gradient(p);
}

template< class Param, class Vect, class Mat, class Real >
Vect *tNRJ< Param, Vect, Mat, Real >::Gradient(const Param &) {
  return NULL;
}

template< class Param, class Vect, class Mat, class Real >
Mat *tNRJ< Param, Vect, Mat, Real >::getHessian(const Param &p) {
  nappel_hess++;
  return Hessian(p);
}

template< class Param, class Vect, class Mat, class Real >
Mat *tNRJ< Param, Vect, Mat, Real >::Hessian(const Param &) {
  return NULL;
}

#endif
