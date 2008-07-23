// Example C++ function "add random function ", dynamically loaded into "load.edp"                                                     
// ---------------------------------------------------------------------                                                     
// $Id$                                                              

#include "config.h"
#include  <iostream>
#include  <cfloat>
using namespace std;
#include "error.hpp"
#include "AFunction.hpp"
#include <cstdlib>

#ifdef HAVE_SRANDOMDEV
template<class R>
class  OneOperator_0 : public OneOperator {
  class E_F0_F :public  E_F0mps { public:
    typedef  R (*func)( ) ; 
    func f;
    E_F0_F(func ff)  : f(ff) {}
    AnyType operator()(Stack )  const {return SetAny<R>( f()) ;}  
    operator aType () const { return atype<R>();} 

  };

  typedef  R (*func)() ; 
  func  f;
public: 
  E_F0 * code(const basicAC_F0 & ) const 
  { return  new E_F0_F(f);} 
  OneOperator_0(func  ff): OneOperator(map_type[typeid(R).name()]),f(ff){}
};

long ffsrandom(long  s) { srandom( (unsigned int ) s); return 0;}
long ffsrandomdev() { srandomdev(); return 0;}

class Init { public:
  Init();
};
Init init;



Init::Init(){
  Global.Add("srandomdev","(",new OneOperator_0<long>(ffsrandomdev));
  Global.Add("srandom","(",new OneOperator1<long>(ffsrandom));
  Global.Add("random","(",new OneOperator_0<long>(random));
}
 
/* These real versions are due to Isaku Wada, 2002/01/09 added */
#endif
