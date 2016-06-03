// Example C++ function "myfunction", dynamically loaded into "load.edp"
// ---------------------------------------------------------------------
// $Id$
#include <ff++.hpp>
#include <AFunction_ext.hpp>
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

Complex GPvortex(Stack stack,const double & x0,const double &y0,const double &kappa)
{
    MeshPoint &mp= *MeshPointStack(stack); // the struct to get x,y, normal , value
   // double *d = *pd;
    double x= mp.P.x; // get the current x value
    double y= mp.P.y; // get the current y value
    Complex p((x-x0)*kappa,(y-y0)*kappa);
    double r = abs(p);// theta=arg(p);
    double tanhr=tanh(r);
    double tr = tanhr/r ;
    return (r>1e-20) ? tr*p : p ;
    
}

Complex dxGPvortex(Stack stack,const double & x0,const double &y0,const double &kappa)
{
    MeshPoint &mp= *MeshPointStack(stack); // the struct to get x,y, normal , value
    // double *d = *pd;
    double x= mp.P.x; // get the current x value
    double dx=1.;
    double y= mp.P.y; // get the current y value
    Complex p((x-x0)*kappa,(y-y0)*kappa);
    Complex dp((dx)*kappa,0);
    
    double r = abs(p),dr=p.real()*dx/r;
    double tanhr=tanh(r),dtanhr=(1-tanhr*tanhr)*dr;
    double tr = tanhr/r,dtr = (dtanhr/r)- dr*tanhr/(r*r) ;
    
    return (r>1e-20) ? ( dtr*p + tr*dp )  : dp ;
    

}

Complex dyGPvortex(Stack stack,const double & x0,const double &y0,const double &kappa)
{
    MeshPoint &mp= *MeshPointStack(stack); // the struct to get x,y, normal , value
    // double *d = *pd;
    double x= mp.P.x; // get the current x value
    double dy=1.;
    double y= mp.P.y; // get the current y value
    Complex p((x-x0)*kappa,(y-y0)*kappa);
    Complex dp(0,(dy)*kappa);
    
    double r = abs(p),dr=p.imag()*dy/r;
    double tanhr=tanh(r),dtanhr=(1-tanhr*tanhr)*dr;
    double tr = tanhr/r,dtr = (dtanhr/r)- dr*tanhr/(r*r) ;
    
    return (r>1e-20) ? ( dtr*p + tr*dp )  : dp ;
    
}

Complex GPvortices(Stack stack,const KNM_<double> & ps)
{
    MeshPoint &mp= *MeshPointStack(stack); // the struct to get x,y, normal , value
    // double *d = *pd;
    double x= mp.P.x; // get the current x value
    double y= mp.P.y; // get the current y value
    Complex vs=1.;
    ffassert(ps.M()>=3);
    for(int i=0; i<ps.N();i++)
    {
    double x0=ps(i,0), y0=ps(i,1),kappa = ps(i,2);
    Complex p((x-x0)*kappa,(y-y0)*kappa);
    double r = abs(p);// theta=arg(p);
    double tanhr=tanh(r);
    double tr = tanhr/r ;
    vs*= (r>1e-20) ? tr*p : p ;
    }
    return vs;
}

static void init(){
  Global.Add("BECtrap","(",new OneOperator1s_<double,KN<double> * ,E_F_F0s_<double,KN<double> * ,E_F0mps > >(BECtrap));
  Global.Add("GPvortex","(",
              new OneOperator3s_<Complex,double,double,double,
                                 E_F_F0F0F0s_<Complex,double,double,double,E_F0mps> > (GPvortex));
     Global.Add("GPvortices","(",
               new OneOperator1s_<Complex,KNM_<double> ,
                                  E_F_F0s_<Complex,KNM_<double>,E_F0mps > >(GPvortices));
    Global.Add("dxGPvortex","(",new OneOperator3s_<Complex,double,double,double,E_F_F0F0F0s_<Complex,double,double,double,E_F0mps> >
               (dxGPvortex));
    Global.Add("dyGPvortex","(",new OneOperator3s_<Complex,double,double,double,E_F_F0F0F0s_<Complex,double,double,double,E_F0mps> >
               (dyGPvortex));

    
}

LOADFUNC(init);
