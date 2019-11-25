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

// Param class, inherited from KN
// This is the class of parameters against which we optimize
// We can possibly add constraints

#ifndef PARAM_HH
#define PARAM_HH

#include <iostream>
using namespace std;

#include "defs.hpp"
//#include "mvvtp.h"
//#include "mvblas.h"

// We nees a Vect type for update function
template< class Real >
class Param : public KN< Real > {
 private:
  int ndim;
  // These bound can be used in the update method of LineSearch.hpp
  // This is a way to implement a projected gradient
  KN< Real >* maxParam;
  KN< Real >* minParam;

 public:
  Param( );
  // construct a continuous model with size n
  Param(int n);
  // construct a continuous model space with boundary and initial values
  Param(const KN< Real >& maxp, const KN< Real >& minp, const KN< Real >& initmod);
  // construct a continuous model space with boundary
  Param(const KN< Real >& maxp, const KN< Real >& minp);
  // construct a continuous model space with initial values
  Param(const KN< Real >& initmod);
  // copy operator
  Param(const Param< Real >& p);
  virtual ~Param( );

  Param< Real >& operator=(const Param< Real >&);

  KN< Real >* modMax( ) const;
  KN< Real >* modMin( ) const;
  void setModMax(const KN< Real >& v);
  void setModMin(const KN< Real >& v);

  // cf. Stroustrup page 612-613 for operator <<
  virtual ostream& toto(ostream&) const;

 private:    // no copy
};

// Constructors
template< class Real >
Param< Real >::Param(void) : KN< Real >( ) {
  // cerr << "Param default constructor" << endl;
  maxParam = NULL;
  minParam = NULL;
}

template< class Real >
Param< Real >::Param(int n) : KN< Real >(n) {
  // cerr << "Param constructor 0" << endl;
  maxParam = NULL;
  minParam = NULL;
}

template< class Real >
Param< Real >::Param(const KN< Real >& maxp, const KN< Real >& minp, const KN< Real >& initmod)
  : KN< Real >(initmod) {
  int ndim = initmod.size( );

  // cerr << "Param constructor 1" << endl;

  maxParam = new KN< Real >(ndim);
  minParam = new KN< Real >(ndim);
  maxParam[0] = maxp;
  minParam[0] = minp;
}

template< class Real >
Param< Real >::Param(const KN< Real >& maxp, const KN< Real >& minp)
  : KN< Real >(Min(minp.size( ), maxp.size( ))) {
  int ndim = Min(minp.size( ), maxp.size( ));

  // cerr << "Param constructor 2" << endl;

  maxParam = new KN< Real >(ndim);
  minParam = new KN< Real >(ndim);
  maxParam[0] = maxp;
  minParam[0] = minp;
}

template< class Real >
Param< Real >::Param(const KN< Real >& initmod) : KN< Real >(initmod) {
  // cerr << "Param constructor 3" << endl;

  maxParam = NULL;
  minParam = NULL;
}

// copy operator
template< class Real >
Param< Real >::Param(const Param& p) : KN< Real >(p) {
  // cerr << "Param copy operator" << endl;

  if ((p.maxParam) == NULL)
    maxParam = NULL;
  else {
    maxParam = new KN< Real >(ndim);
    *maxParam = *(p.maxParam);
  }

  if ((p.minParam) == NULL)
    minParam = NULL;
  else {
    minParam = new KN< Real >(ndim);
    *minParam = *(p.minParam);
  }
}

template< class Real >
Param< Real >& Param< Real >::operator=(const Param< Real >& p) {
  // cerr << "Param = operator" << endl;
  // modif FH 042005 for gcc4.0
  KN< Real >& a1 = *this;
  const KN< Real >& a2 = p;
  a1 = a2;    // KN copy operator

  int ndim = p.size( );

  if ((p.maxParam) == NULL)
    maxParam = NULL;
  else {
    if (maxParam) delete maxParam;
    maxParam = new KN< Real >(ndim);
    *maxParam = *(p.maxParam);
  }

  if ((p.minParam) == NULL)
    minParam = NULL;
  else {
    if (minParam) delete minParam;
    minParam = new KN< Real >(ndim);
    *minParam = *(p.minParam);
  }

  return (*this);
}

template< class Real >
Param< Real >::~Param( ) {
  // cerr << "Param destructor" << endl;

  if (maxParam != NULL) delete maxParam;
  if (minParam != NULL) delete minParam;
}

template< class Real >
KN< Real >* Param< Real >::modMax( ) const {
  return maxParam;
}

template< class Real >
KN< Real >* Param< Real >::modMin( ) const {
  return minParam;
}

template< class Real >
void Param< Real >::setModMax(const KN< Real >& v) {
  if (maxParam != NULL) delete maxParam;
  maxParam = new KN< Real >(v);
}

template< class Real >
void Param< Real >::setModMin(const KN< Real >& v) {
  if (minParam != NULL) delete minParam;
  minParam = new KN< Real >(v);
}

template< class Real >
ostream& Param< Real >::toto(ostream& os) const {
  for (long i = 0; i < (*this).size( ); i++) os << (*this)[i] << " ";
  os << endl;

  return os;
}

template< class Real >
ostream& operator<<(ostream& os, const Param< Real >& d) {
  return d.toto(os);
}

#endif
