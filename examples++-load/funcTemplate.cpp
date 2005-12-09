// Example C++ function "CppModTemplate" dynamically loaded into "load.edp"
// ------------------------------------------------------------------------

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

// see src/femlib/RNM.hpp

// dummy routine to understand how to use vector
double CppModTemplate3(KN<double> *const & A,                            // OUTPUT
		       KN<double> *const & B, KN<double> *const & C)     // INPUTS
{
  
  // Remarque:
  // It might prove usefull to have a look in the cpp file where KN is defined: src/femlib/RNM.hpp
  //
  // To access value at node i of vector N, do as follow: *(N[0]+i)
  // Explanation (C++ for dummies as I am ;-):
  //   N         is an alias to the KN object.
  //   N[0]      is a pointer to the first element of the vector.
  //   N[0]+i    is a pointer to the ith element of the vector.
  //   *(N[0]+i) is the value of the ith element of the vector.
  
  int nn = A->N(); // get number of nodes

  cout << "nn: " << nn << endl;
  
  for(int i=0; i<nn; i++) {
    (*(A[0]+i)) = (*(B[0]+i)) *  (*(C[0]+i));
    cout << (*(A[0]+i)) << endl;
  }

  return 0.0;  // dummy return value.
}

double CppModTemplate4(KN<double> *const & A,                            // OUTPUT
		       KN<double> *const & B, KN<double> *const & C,     // INPUTS
		       KN<double> *const & D)   
{
  int nn = A->N(); // get number of nodes
  cout << "nn: " << nn << endl;
  for(int i=0; i<nn; i++) {
    (*(A[0]+i)) = (*(B[0]+i)) *  (*(C[0]+i)) * (*(D[0]+i));
    cout << (*(A[0]+i)) << endl;
  }
  return 0.0;  // dummy return value.
}

double CppModTemplate5(KN<double> *const & A,                            // OUTPUT
		       KN<double> *const & B, KN<double> *const & C,     // INPUTS
		       KN<double> *const & D, KN<double> *const & E)   
{
  int nn = A->N(); // get number of nodes
  cout << "nn: " << nn << endl;
  for(int i=0; i<nn; i++) {
    (*(A[0]+i)) = (*(B[0]+i)) *  (*(C[0]+i)) * (*(D[0]+i)) * (*(E[0]+i));
    cout << (*(A[0]+i)) << endl;
  }
  return 0.0;  // dummy return value.

}


double CppModTemplate6(KN<double> *const & A,                            // OUTPUT
		       KN<double> *const & B, KN<double> *const & C,     // INPUTS
		       KN<double> *const & D, KN<double> *const & E,
		       KN<double> *const & F)   
{
  int nn = A->N(); // get number of nodes
  cout << "nn: " << nn << endl;
  for(int i=0; i<nn; i++) {
    (*(A[0]+i)) = (*(B[0]+i)) *  (*(C[0]+i)) * (*(D[0]+i)) * (*(E[0]+i)) * (*(F[0]+i));
    cout << (*(A[0]+i)) << endl;
  }
  return 0.0;  // dummy return value.
}

double CppModTemplate7(KN<double> *const & A,                            // OUTPUT
		       KN<double> *const & B, KN<double> *const & C,     // INPUTS
		       KN<double> *const & D, KN<double> *const & E,
		       KN<double> *const & F, KN<double> *const & G)   
{
  int nn = A->N(); // get number of nodes
  cout << "nn: " << nn << endl;
  for(int i=0; i<nn; i++) {
    (*(A[0]+i)) = (*(B[0]+i)) *  (*(C[0]+i)) * (*(D[0]+i)) * (*(E[0]+i)) * (*(F[0]+i)) * (*(G[0]+i));
    cout << (*(A[0]+i)) << endl;
  }
  return 0.0;  // dummy return value.
}

double CppModTemplate8(KN<double> *const & A,                            // OUTPUT
		       KN<double> *const & B, KN<double> *const & C,     // INPUTS
		       KN<double> *const & D, KN<double> *const & E,
		       KN<double> *const & F, KN<double> *const & G,
		       KN<double> *const & H)   
{
  int nn = A->N(); // get number of nodes
  cout << "nn: " << nn << endl;
  for(int i=0; i<nn; i++) {
    (*(A[0]+i)) = (*(B[0]+i)) *  (*(C[0]+i)) * (*(D[0]+i)) * (*(E[0]+i)) * (*(F[0]+i)) * (*(G[0]+i)) * (*(H[0]+i)) ;
    cout << (*(A[0]+i)) << endl;
  }
  return 0.0;  // dummy return value.
}



//   add the function name to the freefem++ table 
class Init { public:
  Init();
};
Init init;
Init::Init(){
  // Add function with 3 arguments
  Global.Add("CppModTemplate3","(",new OneOperator3_<double, KN<double>*, KN<double>*, KN<double>*>(CppModTemplate3)); 
  Global.Add("CppModTemplate4","(",new OneOperator4_<double, KN<double>*, KN<double>*, KN<double>*, KN<double>*>(CppModTemplate4)); 
  Global.Add("CppModTemplate5","(",new OneOperator5_<double, KN<double>*, KN<double>*, KN<double>*, KN<double>*, KN<double>*>(CppModTemplate5)); 
  Global.Add("CppModTemplate6","(",new OneOperator6_<double, KN<double>*, KN<double>*, KN<double>*, KN<double>*, KN<double>*, KN<double>*>(CppModTemplate6)); 
  Global.Add("CppModTemplate7","(",new OneOperator7_<double, KN<double>*, KN<double>*, KN<double>*, KN<double>*, KN<double>*, KN<double>*, KN<double>*>(CppModTemplate7)); 
  Global.Add("CppModTemplate8","(",new OneOperator8_<double, KN<double>*, KN<double>*, KN<double>*, KN<double>*, KN<double>*, KN<double>*, KN<double>*, KN<double>*>(CppModTemplate8)); 
}



