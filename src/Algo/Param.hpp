// Class Param
// C'est la classe des paramètres par rapport auxquels on optimise
// Derive de la classe KN
// On pourra éventuellement lui ajouter des contraintes

#ifndef PARAM_HH
#define PARAM_HH

#include <iostream>
using namespace std;

#include "defs.hpp"
//#include "mvvtp.h"
//#include "mvblas.h"

// on a besoin du type Vect pour la fonction update
template <class Real>
class Param: public KN<Real> {
  
private:
  int	ndim;
  // Ces bornes servent éventuellement dans la méthode update de LineSearch.hpp
  // C'est une manière d'implémenter le gradient projeté...
  KN<Real>*	maxParam;
  KN<Real>*	minParam;
  
public:
  Param();
  //construct a continuous model with size n
  Param(int n);
  //construct a continuous model space with boundary and initial values
  Param(const KN<Real>& maxp, const KN<Real>& minp, const KN<Real>& initmod);
  //construct a continuous model space with boundary
  Param(const KN<Real>& maxp, const KN<Real>& minp);
  //construct a continuous model space with initial values 
  Param(const KN<Real>& initmod);
  //operateur de copie
  Param(const Param<Real>& p);
  virtual ~Param();

  Param<Real>& operator=(const Param<Real>&);
  
  KN<Real>* modMax() const;
  KN<Real>* modMin() const;
  void setModMax(const KN<Real>& v);
  void setModMin(const KN<Real>& v);
  
  //cf. Stroustrup page 612-613 pour l'implémentation de operator<<
  virtual ostream& toto(ostream&) const;
  
};

// Constructors

template <class Real>
Param<Real>::Param(void)
  :KN<Real>()
{
  //  cerr<<"Constructeur par défaut Param"<<endl;
  
   maxParam = NULL;
   minParam = NULL;
}


template <class Real>
Param<Real>::Param(int n)
  :KN<Real>(n)
{
  //  cerr<<"Constructeur 0 Param"<<endl;
  
   maxParam = NULL;
   minParam = NULL;
}


template <class Real>
Param<Real>::Param(const KN<Real>& maxp, const KN<Real>& minp, const KN<Real>& initmod)
  :KN<Real>(initmod)
{
  int ndim=initmod.size();

  //  cerr<<"Constructeur 1 Param"<<endl;
  
  maxParam = new KN<Real>(ndim);
  minParam = new KN<Real>(ndim);
  maxParam[0] = maxp;
  minParam[0] = minp;
}

template <class Real>
Param<Real>::Param(const KN<Real>& maxp, const KN<Real>& minp)
  :KN<Real>(Min(minp.size(),maxp.size()))
{
  int ndim=Min(minp.size(),maxp.size());

  //  cerr<<"Constructeur 2 Param"<<endl;
	
  maxParam = new KN<Real>(ndim);
  minParam = new KN<Real>(ndim);
  maxParam[0] = maxp;
  minParam[0] = minp;
}


template <class Real>
Param<Real>::Param(const KN<Real>& initmod)
  :KN<Real>(initmod)
{
  //  cerr<<"Constructeur 3 Param"<<endl;
  
  maxParam = NULL;
  minParam = NULL;
  
}

// opérateur de copie
template <class Real>
Param<Real>::Param(const Param& p)
  :KN<Real>(p)
{
  //  cerr<<"Operateur de copie de Param"<<endl;
  
  if ((p.maxParam)==NULL)
	maxParam=NULL;
  else{
	maxParam = new KN<Real>(ndim);
	*maxParam=*(p.maxParam);
  }
  
  if ((p.minParam)==NULL)
	minParam=NULL;
  else{
	minParam = new KN<Real>(ndim);
	*minParam=*(p.minParam);
  }
}

template<class Real>
Param<Real>& Param<Real>::operator=(const Param<Real>& p)
{
  //  cerr<<"Operateur = de Param"<<endl;

  KN<Real>* a1;
  KN<Real>* a2;
  a1=this;
  a2=(KN<Real> *)(&p);
  *a1=*a2; // Operateur de copie de KN
  
  int ndim=p.size();
  
  if ((p.maxParam)==NULL)
	maxParam=NULL;
  else{
	maxParam = new KN<Real>(ndim);
	*maxParam=*(p.maxParam);
  }
  
  if ((p.minParam)==NULL)
	minParam=NULL;
  else{
	minParam = new KN<Real>(ndim);
	*minParam=*(p.minParam);
  }
 
  return (*this);
  
}
 
template<class Real>
Param<Real>::~Param()
{
  //  cerr<<"Destructeur de Param"<<endl;

  if (maxParam != NULL)  delete maxParam;
  if (minParam != NULL)  delete minParam;

}



template<class Real>
KN<Real>* Param<Real>::modMax() const
{
  return maxParam;
}
 
template<class Real>
KN<Real>* Param<Real>::modMin() const
{
  return minParam;
}
 
template<class Real>
void Param<Real>::setModMax(const KN<Real>& v)
{ 
  if (maxParam!=NULL) delete maxParam;
  maxParam= new KN<Real>(v);
	}
 
template<class Real>
void Param<Real>::setModMin(const KN<Real>& v)
{
  if (minParam!=NULL) delete minParam;
  minParam= new KN<Real>(v);
	}



//cf. Stroustrup page 612-613 pour l'implémentation de operator<<

template<class Real>
ostream& Param<Real>::toto (ostream& os) const
{
  for (long i=0;i<(*this).size();i++)
	os << (*this)[i]<<" ";
  os<<endl;
  
  return os;
}

template<class Real>
ostream& operator <<(ostream& os, const Param<Real>& d)
{
  return d.toto(os);
}


#endif



