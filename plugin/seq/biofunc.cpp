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

#include <ff++.hpp>
#include <AFunction_ext.hpp>
#include <cstdlib>
using namespace Fem2D;
typedef double R;
double fmonod(const double &xx,const double &k,const double &kb)
{double x=max(xx,0.); return (k*x)/(kb+x);}
double dfmonod(const double &xx,const double& k,const double &kb)
{double x = max(xx,0.), a= kb+x;return k/a - k*x/(a*a);}
double fmonod(const double &xx,const double &kb)
{double x=max(xx,0.); return (x)/(kb+x);}
double dfmonod(const double &xx,const double &kb)
{double x = max(xx,0.), a= kb+x;return 1./a - x/(a*a);}

static void init () {
    Global.Add("fmonod", "(", new OneOperator3_<R,R,R,R >(fmonod));
    Global.Add("dfmonod", "(", new OneOperator3_<R,R,R,R>(dfmonod));
    Global.Add("fmonod", "(", new OneOperator2_<R,R,R >(fmonod));
    Global.Add("dfmonod", "(", new OneOperator2_<R,R,R>(dfmonod));
}

LOADFUNC(init);
