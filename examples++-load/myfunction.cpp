#include  <iostream>
#include  <cfloat>
using namespace std;
#include "error.hpp"
#include "AFunction.hpp"
using namespace std;  
#include "rgraph.hpp"
#include "RNM.hpp"
#include "fem.hpp"
#include "FESpace.hpp" 
#include "MeshPoint.hpp"

using namespace Fem2D;
double myfunction(Stack stack)
{
  //  to get FreeFem++  data
  MeshPoint &mp= *MeshPointStack(stack);
  double x= mp.P.x;
  double y= mp.P.y;
  //  cout << "x = " << x << " y=" << y << endl;
  return sin(x)*cos(y);
}

//  a class ------
template<class R>
class  OneOperator0s : public OneOperator {
  virtual bool MeshIndependent() const {return false;} // 
  class E_F0_F :public  E_F0mps { public:
    typedef  R (*func)(Stack stack) ; 
    func f;
    E_F0_F(func ff)  : f(ff) {}
    AnyType operator()(Stack stack)  const {return SetAny<R>( f(stack)) ;}  
    operator aType () const { return atype<R>();} 
  };

  //  aType r; //  return type
  typedef  R (*func)(Stack ) ; 
  func  f;
public: 
  E_F0 * code(const basicAC_F0 & ) const 
  { return  new E_F0_F(f);} 
   OneOperator0s(func  ff): OneOperator(map_type[typeid(R).name()]),f(ff){}
};

//  hack to do something at in

class Init { public:
  Init();
};
Init init;
Init::Init(){
  cout << " lood: myfunction " << endl;
  Global.Add("myfunction","(",new OneOperator0s<double>(myfunction));

}


