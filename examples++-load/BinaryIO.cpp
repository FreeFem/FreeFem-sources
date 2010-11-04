
#include  <iostream>
#include  <cfloat>
using namespace std;
#include "error.hpp"
#include "AFunction.hpp"
#include "rgraph.hpp"
#include "RNM.hpp"
#include "fem.hpp"
#include "FESpace.hpp" 
#include "MeshPoint.hpp"
#include "AFunction_ext.hpp" // Extension of "AFunction.hpp" to deal with more than 3 parameters function
using namespace Fem2D;

 
double SaveVec(KN<double> *const & f, string *const & nome)   
{
  std::ofstream outfile (nome->data(),ios_base::binary);
  // To access value at node i of vector N, do as follow: *(N[0]+i)
  // Explanation (C++ for dummies as I am ;-):
  //   N         is an alias to the KN object.
  //   N[0]      is a pointer to the first element of the vector.
  //   N[0]+i    is a pointer to the ith element of the vector.
  //   *(N[0]+i) is the value of the ith element of the vector.
 
  long int nn = f->N(); // get number of nodes
  long int dim=nn;
  outfile.write ((char*) &dim, sizeof(long int));//write the dimension of the vector
  double ftemp ;
  for(long int i=0; i<nn; i++) {
    
    ftemp = *(f[0]+i) ;
    outfile.write ((char*) &ftemp, sizeof(double));
  }
  outfile.close();
  return 0.0;  // dummy return value.
}

 
double LoadVec(KN<double> *const & ww, string *const & nome)   
{
  std::ifstream infile (nome->data(),ios_base::binary);
  long int dim;
  infile.read((char *) &dim, sizeof(long int));
  double dtemp;
  for(long int i=0; i<dim; i++)
   {
   infile.read((char *) &dtemp, sizeof(double));
    *(ww[0]+i)=dtemp ;
    }
  return 0.0;  // dummy return value.
}

 
double LoadFlag(long int *const & ww, string *const & nome)   
{
  std::ifstream infile (nome->data(),ios_base::binary);
  long int flag;
  infile.read((char *) &flag, sizeof(long int));
  *ww=flag;
  return 0.0;  // dummy return value.
}

double flag(long int *const & FLAG,string *const &nome)   
{
  std::ofstream outfile (nome->data(),ios_base::binary);
long int Flag;
 Flag= *FLAG;   
 outfile.write ((char*) &Flag, sizeof(long int));
outfile.close();
}



//   add the function name to the freefem++ table 
class Init { public:
  Init();
};	
Init init;
Init::Init(){	
		
  Global.Add("LoadVec","(",new OneOperator2_<double,  KN<double>*, string* >(LoadVec)); 
Global.Add("LoadFlag","(",new OneOperator2_<double,long int*, string* >(LoadFlag));
  Global.Add("SaveVec","(",new OneOperator2_<double,KN<double>*, string* >(SaveVec));
 Global.Add("flag","(",new OneOperator2_<double,long int*,string* >(flag));  
}

