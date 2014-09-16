#ifndef OPTIMA_HH
#define OPTIMA_HH

#include "NRJ.hpp"
#include <limits.h>
#include "defs.hpp"

#define MAX_IT_OPT 1000

//==================
//define the base class for  optimization classes
//	H. Lydia Deng, 03/14/94
//======================

/*
  This is the base class for general optimization classes.
  Classes derived from this class inherit features of 
  Optima. This class cannot be instantiated directly.
*/

// Attention :
// pas de copie
// pas d'opérateur =

template <class LS>
class Optima {
  typedef typename LS::Real Real;
  typedef typename LS::Param Param;
  typedef typename LS::Vect Vect;
  typedef typename LS::VMat VMat;
  typedef LS  LineSearch;
  typedef list<Real> mlist;

protected:
  // maximum number of iterations
  int iterMax;
  int iterNum;
  // tolerance error
  // c'est une tolérance sur la norme euclidienne du gradient
  Real  tol;
  // mlist of residue 
  mlist *  residue;
  // pointer to LS
  LS * ls;
  // verbose or quiet
  int  isVerbose;
  // the flag indicating if the search was a success
  int  isSuccess; 
  // append the new residue to the mlist
  void appendResidue(Real res); 
  
public:
  
  Optima(int verb=0);
  virtual ~Optima();
  
  virtual Param optimizer(Param&) = 0;
  
  //Output residues

  int   ifSuccess();
  //  residue of the last iteration
  Real finalResidue();
  // residue of the first iteration
  Real firstResidue();
  //  residues of all iterations
  mlist allResidue();
  //  normalized residues of all iterations
  mlist normResidue();

};

template <class LS>
Optima<LS>::Optima(int verbose){ 
  residue 	= 	new mlist; 
  isVerbose 	= 	verbose;
  isSuccess       =       0;
  ls		=	NULL;
  tol		=	0.;
  iterMax		=	MAX_IT_OPT;
}

template <class LS>
Optima<LS>::~Optima() { 
  if (residue!=NULL) delete residue;
}

//Output residues
template <class LS>
int Optima<LS>::ifSuccess()
{return isSuccess;}

template <class LS>
typename Optima<LS>::Real Optima<LS>::finalResidue()
{return residue->back();}

template <class LS>
typename Optima<LS>::Real Optima<LS>::firstResidue()
{return residue->front();}

template <class LS>
typename Optima<LS>::mlist Optima<LS>::allResidue()
{return *residue;}

template <class LS>
typename Optima<LS>::mlist Optima<LS>::normResidue()
{return normalize(*residue);} //cf. defs.hpp

template <class LS>
void Optima<LS>::appendResidue(Real res)
{
  (residue[0]).push_back(res);
}      

#endif

