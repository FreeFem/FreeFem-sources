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

double testio(Stack stack)
{
  double x=M_PI;
  long l= (1<<9) ;
  cout << " test cout " << x << endl;
  cout << " test cout " << l << endl;

  cerr << " test cerr " << x << endl;
  cerr << " test cerr " << l << endl;
  return 0.; 
}

static void init(){
  Global.Add("myfunction","(",new OneOperator0s<double>(myfunction));
  Global.Add("testio","(",new OneOperator0s<double>(testio));
}

LOADFUNC(init);
