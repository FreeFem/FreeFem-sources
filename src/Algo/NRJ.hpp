#ifndef NRJ_HH
#define NRJ_HH

#include "defs.hpp"

// Attention :
// pas de copie
// pas d'opérateur =

template <class Param, class Vect, class Mat,class Real>
class NRJ{
protected:

  int				nappel_val;
  int				nappel_grad;
  int				nappel_hess;
  int               nparam;
  
  Real              val;
  Vect*             grad;
  Mat*              hess;
  
public:

  // Attention : dans le constructeur de l'NRJ, il faut vraiment
  // initialiser grad ou hess, si on les utilise... Ce n'est pas
  // fait ici (parce que l'initialisation dépend trop du type)
  NRJ(int);
  virtual ~NRJ();
  
  // à preciser pour chaque fonction NRJ
  virtual Real Val(const Param&) = 0;
  
  Real getVal(const Param&);

  // à preciser éventuellement, sinon ça rend nul
  virtual Vect* Gradient(const Param&);
  
  Vect* getGradient(const Param&);

  // à preciser éventuellement, sinon ça rend nul
  virtual Mat* Hessian(const Param&);

  Mat* getHessian(const Param&);
  

  int Appel_Val() const {return nappel_val;};
  int Appel_Grad() const {return nappel_grad;};
  int Appel_Hess() const {return nappel_hess;};
  

};

template <class Param, class Vect, class Mat,class Real>
NRJ<Param,Vect, Mat, Real>::~NRJ()
{
  if (grad!=NULL) delete grad;
  if (hess!=NULL) delete hess;
}

template <class Param, class Vect, class Mat,class Real>
NRJ<Param,Vect, Mat, Real>::NRJ(int n)
{
  nparam=n;
  nappel_val=0;
  nappel_grad=0;
  nappel_hess=0;
  grad = NULL;
  hess = NULL;
  val=0.;
}

template <class Param, class Vect, class Mat,class Real>
Real NRJ<Param,Vect, Mat, Real>::getVal(const Param& p)
{
	nappel_val++;
	return Val(p);
}

template <class Param, class Vect, class Mat,class Real>
Vect* NRJ<Param,Vect, Mat, Real>::getGradient(const Param& p)
{
	nappel_grad++;
	return Gradient(p);
}

template <class Param, class Vect, class Mat,class Real>
Vect* NRJ<Param,Vect, Mat, Real>::Gradient(const Param&)
{
  return NULL;
}

template <class Param, class Vect, class Mat,class Real>
Mat* NRJ<Param,Vect, Mat, Real>::getHessian(const Param& p)
{
	nappel_hess++;
	return Hessian(p);
}

template <class Param, class Vect, class Mat,class Real>
Mat* NRJ<Param,Vect, Mat, Real>::Hessian(const Param&)
{
  return NULL;
}


#endif
