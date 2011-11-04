// Example C++ function "myfunction", dynamically loaded into "load.edp"
// ---------------------------------------------------------------------
// $Id$
#include <ff++.hpp>
using namespace Fem2D;
double myfunction(Stack stack)
{
  //  to get FreeFem++  data
  MeshPoint &mp= *MeshPointStack(stack); // the struct to get x,y, normal , value 
  double x= mp.P.x; // get the current x value
  double y= mp.P.y; // get the current y value
  // cout << "x = " << x << " y=" << y << " " << sin(x)*cos(y) <<  endl;
  return sin(x)*cos(y);
}

//  A class build the link with freefem++ 
// generaly this class are already in AFunction.hpp
// but unfortunatly, I have no simple function with no parameter
// in freefem++ depending of the mesh,  
template<class R>
class  OneOperator0s : public OneOperator {

  // the class to defined a evaluated a new function 
  //  It  must devive from  E_F0 if it is mesh independent
  //   or from E_F0mps if it is mesh dependent 
  class E_F0_F :public  E_F0mps { public:   
    typedef  R (*func)(Stack stack) ; 
    func f; // the pointeur to the fnction myfunction 
    E_F0_F(func ff)  : f(ff) {}
    // the operator evaluation in freefem++ 
    AnyType operator()(Stack stack)  const {return SetAny<R>( f(stack)) ;}  

  };

  typedef  R (*func)(Stack ) ; 
  func  f; 
public: 
  // the function which build the freefem++ byte code 
  E_F0 * code(const basicAC_F0 & ) const { return  new E_F0_F(f);} 
  // the constructor to say ff is a function without parameter
  // and returning a R
  OneOperator0s(func  ff): OneOperator(map_type[typeid(R).name()]),f(ff){}
};

//  Hack to do something at initialisation time
//   to add the name myfunction to the freefem++ table 
class Init { public:
  Init();
};
LOADINIT(Init);
Init::Init(){
  Global.Add("myfunction","(",new OneOperator0s<double>(myfunction));
}
