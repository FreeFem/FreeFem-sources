
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cassert>
#include <cmath>
using namespace std; 

#include "RNM.hpp"    
#include "bmo.hpp"    

class BMO_user : public BijanMO  
{ 
public:
  void init(Vect & v) ;    
  double J(Vect & x);
  double * DJ(Vect &x,Vect & fpx);
  BMO_user(): BijanMO(5,5,5,5,5)
  {
      /* 
       BijanMO(
       ndim, 
       nbrestart=1,
       nbext1=1,
       nbbvp=5,
       nbgrad=5,
       epsfd=1e-5,
       rho000=100,
       epsloc=1e-4,
       epsij=1e-6,
       n100=100)
       
       */
  }
};

/* cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* cc               provided by the user */
/* cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */

void  BMO_user::init(Vect & v) 
{

    xmin=-5;
    xmax= 5.;

    for (int ii = 0; ii < ndim; ++ii)
      {
	double	xrd=abs(sin((ii+1)/(ndim*2.)));
	v[ii] = xrd * xmax[ii] + (1 - xrd) * xmin[ii];
      }

    return ;
} 


/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* functional definition */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */

double BMO_user::J(Vect & x)
{

  double f =ndim;
  assert(x.min()> -6 && x.max() < 6);
  for (int ii = 0; ii < ndim; ++ii)
    f=f+(x[ii]*x[ii]-cos(18*x[ii]));

  if(debug>5) 
    {
      cout << "     "  << f << " " ;
      for(int i=0;i<ndim;++i) cout << x[i] << " ";
      cout  << nbeval+1 <<" " <<  nbevalp+1 <<  endl; 
    }
  return f;
} /* func_ */


/* cccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* gradient exact,   no defini => DF */  
/* cccccccccccccccccccccccccccccccccccccccccccccccccccccccc */

double * BMO_user::DJ(Vect & x, Vect & fpx)
{

    for (int ii = 0; ii < ndim; ++ii) 
      fpx[ii] = x[ii] * 2. + sin(x[ii] * 18.) * 18.;
    return fpx;//fpx;
}



int main()
{
  BMO_user u;
  u.main();
  return 0;
}
