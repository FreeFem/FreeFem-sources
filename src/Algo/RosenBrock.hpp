#ifndef ROSENBROCK_HH
#define ROSENBROCK_HH

#include "NRJ.hpp"
#include "Param.hpp"
//#include "mvvtp.h"
//#include "mvblas.h"

template <class Real,class Mat>
class RosenBrock : public NRJ<Param<Real>,KN<Real>,Mat,Real> {
private:
protected:
  
public:
  //Constructors and Destructors
  RosenBrock(int);
  ~RosenBrock();
  
  
  Real Val(const Param<Real>&);
  KN<Real>* Gradient(const Param<Real>&);
  Mat* Hessian(const Param<Real>&);
  
};

template <class Real,class Mat>
RosenBrock<Real,Mat>::RosenBrock(int n)
  :NRJ<Param<Real>,KN<Real>,Mat,Real>(n)
{  // On initialise le gradient
  grad= new KN<Real>(n);
  // On initialise le hessien
  // Si on utlise pas le hessien, il suffit qu'il existe un constructeur
  // pour le type Mat qui prend un parametre entier. C'est vrai pour les
  // doubles par exemple...
  hess= new Mat(n); 
  (*hess)=0;
}

template <class Real,class Mat>
RosenBrock<Real,Mat>::~RosenBrock()
{ ; }

template <class Real,class Mat>
Real RosenBrock<Real,Mat>::Val(const Param<Real>& param)
{	
  Real ti_1, ti, tt;
  Real d=0;

  ti_1=param[0];
  
  for (int i=1; i < nparam; i++) 
	{
      ti = param[i];
      tt = 1-ti_1;
      d += tt*tt;
	  tt=ti-ti_1*ti_1;
	  d += 100*tt*tt;
	  
	  ti_1=ti;
	}

  val=d;
  
  return d;
}

template <class Real,class Mat>
KN<Real>* RosenBrock<Real,Mat>::Gradient(const Param<Real>& param)
{
   Real ti_1,ti;
   ti_1=param[0];
   ti=param[1];

   (*grad)[0]=-200*2*ti_1*(ti-ti_1*ti_1)-2*(1-ti_1);
   
   for (int i=1; i<nparam-1; i++){
	 (*grad)[i]=200*(ti-ti_1*ti_1);
	 ti_1=ti;
	 ti=param[i+1];
	 (*grad)[i] += -200*2*ti_1*(ti-ti_1*ti_1)-2*(1-ti_1);
   }
   
   (*grad)[nparam-1]=200*(ti-ti_1*ti_1);

   
   return	grad;
}


template <class Real,class Mat>
Mat* RosenBrock<Real,Mat>::Hessian(const Param<Real>& param)
{
  Real ti_1,ti;
  ti_1=param[0];
  ti=param[1];

  (*hess)(0,0)=-200*2*(ti-ti_1*ti_1)+200*2*2*ti_1*ti_1+2;
  (*hess)(0,1)=-200*2*ti_1;
  
  for (int i=1; i<nparam-1; i++){
	(*hess)(i,i)=200;
	(*hess)(i,i-1)=-200*2*ti_1;
	ti_1=ti;
	ti=param[i+1];
	(*hess)(i,i) += -200*2*(ti-ti_1*ti_1)+200*2*2*ti_1*ti_1+2;
	(*hess)(i,i+1)=-200*2*ti_1;
  }
  
  (*hess)(nparam-1,nparam-1)=200;
  (*hess)(nparam-1,nparam-2)=-200*2*ti_1;
  
  return hess;
}



#endif



