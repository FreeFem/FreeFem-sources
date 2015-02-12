// Example C++ function "myfunction", dynamically loaded into "load.edp"
// ---------------------------------------------------------------------
// $Id$
#include <ff++.hpp>
using namespace Fem2D;
double myf(string * s) {
  cout << *s << endl;
  return 0.;
}
double f(const double& x) { return x*x+1;} 
//  Hack to do something at initialisation time
//   to add the name myfunction to the freefem++ table 
/*  class Init { public:
  Init();
};
$1 */
static void Load_Init(){
  Global.Add("Why","(",new OneOperator1<double,string*>(myf));
  Global.Add("f","(",new OneOperator1_<double,double>(f));
}
LOADFUNC(Load_Init)
