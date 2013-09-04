//  // ---------------------------------------------------------------------
// $Id$
// compile and link with ff-c++ metricPk.cpp 
#include  <iostream>
#include  <cfloat>
#include  <cmath>
using namespace std;
#include "ff++.hpp" 

using namespace std;








class MetricPk :  public E_F0mps 
{
public:
    static basicAC_F0::name_and_type name_param[] ;
    static const int n_name_param =1;
    
    Expression nargs[n_name_param];// stocker les argunments nommes

  typedef KN_<double>  Result;
  Expression expTh;
  Expression expu;
    
  MetricPk(const basicAC_F0 & args)
  {

    args.SetNameParam();// les arguments nommes
    expTh= to<pmesh>(args[0]);  // a the expression to get the mesh
    expu= to<long>(args[1]);  // a the expression to get the mesh
    /*
    exphmin= to<double>(args[2]);  // a the expression to get the mesh
    exphmax= to<double>(args[3]);  // a the expression to get the mesh
    experr= to<double>(args[4]);  // a the expression to get the mesh
    //  a array expression [ a, b]
    const E_Array * ma= dynamic_cast<const E_Array*>((Expression) args[5]);
    const E_Array * mp= dynamic_cast<const E_Array*>((Expression) args[6]);
    if (ma->size() != 3) CompileError("syntax: MetricKuate(Th,np,o,err,[m11,m12,m22],[xx,yy])");
    if (mp->size() != 2) CompileError("syntax: MetricKuate(Th,np,o,err,[m11,m12,m22],[xx,yy])");
    int err =0;
    m11= CastTo<KN<double> * >((*ma)[0]); // fist exp of the array (must be a  double)
    m12= CastTo<KN<double> * >((*ma)[1]); // second exp of the array (must be a  double)
    m22= CastTo<KN<double> * >((*ma)[2]); // second exp of the array (must be a  double)
    px= CastTo<double * >((*mp)[0]); // fist exp of the array (must be a  double)
    py= CastTo<double * >((*mp)[1]); // second exp of the array (must be a  double)
     */
  }
    double arg(int i,Stack stack,double a) const { return nargs[i] ? GetAny<double>( (*nargs[i])(stack) ): a;}
    long arg(int i,Stack stack,long a) const{ return nargs[i] ? GetAny<long>( (*nargs[i])(stack) ): a;}
    bool arg(int i,Stack stack,bool a) const{ return nargs[i] ? GetAny<bool>( (*nargs[i])(stack) ): a;}
 
    
  MetricPk~(){}

  static ArrayOfaType  typeargs()
  { return  ArrayOfaType(
			 atype<pmesh>(),
			 atype<double>());;
  }
  static  E_F0 * f(const basicAC_F0 & args){ return new MetricPk(args);}
  AnyType operator()(Stack s) const ;// la vrai fonction que fait fair le boulot 

};

basicAC_F0::name_and_type Adaptation::name_param[Adaptation::n_name_param] = {
    { 
	"order",             &typeid(long)}  // à 
  };

AnyType MetricPk:: operator()(Stack s) const 
{
    long order=  arg(0,stack,2);    // coef in the metric
    Mesh * Thh = GetAny<pmesh>((*getmesh)(stack));
    ffassert(Thh);
    Mesh & Th= *TTh; 


}
class Init { public:
  Init();
};
LOADINIT(Init);
Init::Init()
{
  cout << "\n  -- lood: init MetricPk\n";
  Global.Add("MetricPk","(", new OneOperatorCode<MetricKuate >( ));
}
