// Example C++ function "myfunction", dynamically loaded into "load.edp"
// ---------------------------------------------------------------------
// $Id$
#include <ff++.hpp>
using namespace Fem2D;
double BECtrap(Stack stack,KN<double> * const &  pd)
{
  MeshPoint &mp= *MeshPointStack(stack); // the struct to get x,y, normal , value 
  double *d = *pd;
  double x= mp.P.x; // get the current x value
  double y= mp.P.y; // get the current y value
  double z= mp.P.y; // get the current y value
  double x2=x*x, y2=y*y,z2=z*z, r2 = x2+y2, r4 = r2*r2;
  long  n = pd->N(); 

  // cout << "x = " << x << " y=" << y << " " << sin(x)*cos(y) <<  endl;
  double ret ;
  if( n ==4)
    ret = x2*d[0] + y2*d[1] + z2*d[2] + r4*d[3] ; 
  else  if(n==6)
    {
      double s = sin(d[5]*z);
      ret = x2*d[0] + y2*d[1] + z2*d[2] + r4*d[3] + d[4]*s*s; 
    }
  else ffassert(0); // 
  return ret; 
}


void init(){
  Global.Add("BECtrap","(",new OneOperator1s_<double,KN<double> * >(BECtrap));
}

LOADFUNC(init);
