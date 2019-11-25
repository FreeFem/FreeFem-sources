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

// *INDENT-OFF* //
//ff-c++-LIBRARY-dep:
//ff-c++-cpp-dep:
// *INDENT-ON* //

// Example C++ function "myfunction", dynamically loaded into "load.edp"
#define _USE_MATH_DEFINES
#include <ff++.hpp>
using namespace Fem2D;
double myfunction(Stack stack) {
  // to get FreeFem++  data
  MeshPoint &mp = *MeshPointStack(stack);    // the struct to get x,y, normal , value
  double x = mp.P.x;                         // get the current x value
  double y = mp.P.y;                         // get the current y value

  return sin(x) * cos(y);
}

// A class build the link with freefem++
// generaly this class are already in AFunction.hpp
// but unfortunatly, I have no simple function with no parameter
// in freefem++ depending of the mesh,
template< class R >
class OneOperator0s : public OneOperator {
  // the class to defined a evaluated a new function
  // It  must devive from  E_F0 if it is mesh independent
  // or from E_F0mps if it is mesh dependent
  class E_F0_F : public E_F0mps {
   public:
    typedef R (*func)(Stack stack);
    func f;    // the pointeur to the fnction myfunction
    E_F0_F(func ff) : f(ff) {}

    // the operator evaluation in freefem++
    AnyType operator( )(Stack stack) const { return SetAny< R >(f(stack)); }
  };

  typedef R (*func)(Stack);
  func f;

 public:
  // the function which build the freefem++ byte code
  E_F0 *code(const basicAC_F0 &) const { return new E_F0_F(f); }

  // the constructor to say ff is a function without parameter
  // and returning a R
  OneOperator0s(func ff) : OneOperator(map_type[typeid(R).name( )]), f(ff) {}
};

double testio(Stack stack) {
  double x = M_PI;
  long l = (1 << 9);

  cout << " test cout " << x << endl;
  cout << " test cout " << l << endl;

  cerr << " test cerr " << x << endl;
  cerr << " test cerr " << l << endl;
  return 0.;
}

static void init( ) {
  Global.Add("myfunction", "(", new OneOperator0s< double >(myfunction));
  Global.Add("testio", "(", new OneOperator0s< double >(testio));
}

LOADFUNC(init);
