
// Example C++ function "CppModTemplate" dynamically loaded into "load.edp"
// ------------------------------------------------------------------------

#include  <iostream>
#include  <cfloat>
using namespace std;
#include "error.hpp"
#include "AFunction.hpp"
#include "lex.hpp"
#include "rgraph.hpp"
#include "RNM.hpp"
#include "fem.hpp"
#include "FESpace.hpp" 
#include "MeshPoint.hpp"

using namespace Fem2D;

// see src/femlib/RNM.hpp

class myType { public:
  string * nom;
  myType(char * nn) { cout << " nn = " << nn << endl;  }
  double x(double u,double v) const { return u+v;}
  void init() { cout << " init myTpe \n" ;  nom =0;} // init des pointeur
  void destroy() { cout << " destroy de la variable associe \n";  delete  nom;nom=0; } 
};

class myType_uv { public:
  myType * mt;
  double u,v;
  myType_uv(myType * mmt,double uu,double vv): mt(mmt),u(uu),v(vv) {}
};

// le vrai constructeur est la
myType * init_MyType(myType * const &a, string * const & s)
{
  a->nom = new string(* s);
  cout << " build MyType " << *a->nom << endl;
} 


myType_uv set_myType_uv( myType * const & mt,const double & u,const double & v)
{  return myType_uv(mt,u,v);}

double get_myType_uv_x(const myType_uv & muv)
{
  return muv.mt->x(muv.u,muv.v);
}

R3  get_myType_uv_N(const myType_uv & muv)
{
  return R3(muv.mt->x(muv.u,muv.v),0.,0.);
}
//   Add the function name to the freefem++ table 
class Init { public:
  Init();
};
Init init;
Init::Init(){

  Dcl_Type<myType*>(InitP<myType >,Destroy<myType>); // declare deux nouveau type pour freefem++  un pointeur et 
  Dcl_Type<myType_uv>();
  //  cast d'un ** en * 
  //  atype<myType**>()->AddCast( new E_F1_funcT<myType*,myType **>(UnRef<myType*>)); 

  zzzfff->Add("myType",atype<myType*>()); // ajoute le type myType a freefem++ 
  // constructeur  d'un type myType  dans freefem 
  TheOperators->Add("<-", 
		    new OneOperator2_<myType *,myType* ,string*>(&init_MyType)); 
  // dans ff++
  //    myType ff("qsdlqdjlqsjdlkq");
  // ajoute la fonction  myType* (u,v) cree le type myType_uv
  //   ff(0.1,0.6).x
  //   deux etapes
  //    1)  ff(u,v) -> myType_uv
  //  ajoute la methode x sur myType_uv   ff(u,v).x 
  // ajoute des fonction sur myType_uv
  // 1)

  atype< myType * >()->Add("(","",new OneOperator3_<myType_uv,myType *,double,double  >(set_myType_uv));  

   Add<myType_uv>("x",".",new OneOperator1_<double,myType_uv>(get_myType_uv_x) );
   Add<myType_uv>("N",".",new OneOperator1_<R3,myType_uv>(get_myType_uv_N) );
}



