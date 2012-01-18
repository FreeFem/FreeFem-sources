// Example C++ function "myfunction", dynamically loaded into "load.edp"
// ---------------------------------------------------------------------
// $Id$
#include <ff++.hpp>
using namespace Fem2D;
double myf(string * s) {
  cout << *s << endl;
  return 0.;
} 
//  Hack to do something at initialisation time
//   to add the name myfunction to the freefem++ table 
class Init { public:
  Init();
};
LOADINIT(Init);
Init::Init(){
  Global.Add("myfunction","(",new OneOperator1<double,string*>(myf));
}
