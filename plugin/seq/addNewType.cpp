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
// SUMMARY : Example C++ function "CppModTemplate" dynamically loaded into "load.edp"
// LICENSE : LGPLv3
// ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
// AUTHORS : ...
// E-MAIL  : ...

#include "ff++.hpp"
using namespace Fem2D;    // see src/femlib/RNM.hpp

/*!
 * \class myType
 * \brief Example class
 */
class myType {
 public:
  string *nom;
  myType(char *nn) { cout << " nn = " << nn << endl; }    // initialization

  double x(double u, double v) const { return u + v; }

  void init( ) {
    cout << " init myTpe " << endl;
    nom = 0;
  }    // pointer initialization

  void destroy( ) {
    cout << " associated variable destruction " << endl;
    delete nom;
    nom = 0;
  }    // pointer delete
};

/*!
 * \class myType_uv
 * \brief Example class
 */
class myType_uv {
 public:
  myType *mt;
  double u, v;
  myType_uv(myType *mmt, double uu, double vv) : mt(mmt), u(uu), v(vv) {}
};

// The real constructor is here
/*!
 * \brief Constructor init_MyType
 * \param a myType *const &
 * \param s string *const &
 * \return NULL
 */
myType *init_MyType(myType *const &a, string *const &s) {
  a->nom = new string(*s);
  cout << " build MyType " << *a->nom << endl;
  return NULL;    // returned value never used for now (13.1)
}

/*!
 * \brief Set a type
 * \param mt myType *const &
 * \param u const double &
 * \param v const double &
 * \return myType_uv(mt, u, v);
 */
myType_uv set_myType_uv(myType *const &mt, const double &u, const double &v) {
  return myType_uv(mt, u, v);
}

/*!
 * \brief Set a type
 * \param muv const myType_uv &
 * \return muv.mt->x(muv.u, muv.v);
 */
double get_myType_uv_x(const myType_uv &muv) { return muv.mt->x(muv.u, muv.v); }

/*!
 * \brief Get a type
 * \param muv const myType_uv &
 * \return r static R3 &
 */
R3 *get_myType_uv_N(const myType_uv &muv) {
  static R3 r;

  r = R3(muv.mt->x(muv.u, muv.v), 0., 0.);
  return &r;
}

// Add the function name to the FreeFEM table
/*!
 * \brief Dynamic load function
 */
static void Load_Init( ) {
  Dcl_Type< myType * >(InitP< myType >, Destroy< myType >);    // declare two new types for FreeFEM
  Dcl_Type< myType_uv >( );
  // Dcl_Type<R3>();
  // cast of ** to *
  // atype<myType**>()->AddCast(new E_F1_funcT<myType*, myType **>(UnRef<myType*>));

  zzzfff->Add("myType", atype< myType * >( ));    // add type myType to FreeFEM
  // constructeur d'un type myType  dans freefem
  TheOperators->Add("<-", new OneOperator2_< myType *, myType *, string * >(&init_MyType));
  // in FreeFEM
  // myType ff("thisisastring");
  // add the function myType* (u,v) & create the type myType_uv

  // ff(0.1, 0.6).x
  // two steps:
  // add the method x on myType_uv
  // add function on myType_uv
  atype< myType * >( )->Add(
    "(", "", new OneOperator3_< myType_uv, myType *, double, double >(set_myType_uv));

  Add< myType_uv >("x", ".", new OneOperator1_< double, myType_uv >(get_myType_uv_x));
  Add< myType_uv >("N", ".", new OneOperator1_< R3 *, myType_uv >(get_myType_uv_N));
}

LOADFUNC(Load_Init)
