#ifndef NEWTONRAPHSON_HH
#define NEWTONRAPHSON_HH

#include "Optima.hpp"
#include "defs.hpp"

/* At the present version, you must use the CubicLS procedure */

// Necessite une "vraie" classe Matrice
template <class LS>
class Newt : public Optima<LS>
{
  typedef typename LS::Real Real;
  typedef typename LS::Param Param;
  typedef typename LS::Vect Vect;
  typedef typename LS::VMat VMat;
  typedef typename LS::Mat Mat;
  typedef list<Real> mlist;
public:
 
  Newt(///pointer to the line-search object.
	   LS * ls, 
	   ///Maximum number of iterations
	   int iter, 
	   ///minimum accepted gradient at optimum solution
	   Real tol,
	   ///vebose or quiet
	   int verb=0);
  ~Newt(){;}
  
  // the Newt search starting from model0, returns an optimum model 
  Param optimizer(Param& model0);

};


template <class LS>
Newt<LS>::Newt(LS* p, int it, Real eps, int verb) 
  :Optima<LS>(verb)
{
  ls=p;
  iterMax 	= 	it;
  tol 		= 	eps;
  iterNum 	= 	0;
}

// Nécessite pour la classe MAT:
// - une méthode de résolution de Ax=b
template <class LS>
typename Newt<LS>::Param Newt<LS>::optimizer( Param& model0)
{
  //reset the residue history for every new optimizer
  iterNum = 0;
  if (residue != NULL)
	{
      delete residue;
      residue = new mlist;
	}
  
  // Initial settings for some parameters
  int n = model0.size();
  Vect g0(n);
  double lambda = 0.025;
  double descent = 0.;

  g0= *(ls->gradient(model0));
  
  // check the gradient, in case the initial model is the optimal
  double err = (Real)sqrt((g0,g0));
  
  if (isVerbose) cerr << "Initial residue : " << err << endl;

  appendResidue(err);	// residual

  if (err < tol)
	{
      if (isVerbose) cerr << "Initial guess was great! \n";
      isSuccess = 1;
      return model0;
	}
   
  Vect g1(n);
  Vect s(n);
  Param model1(model0);
  
  while (finalResidue() > tol && iterNum < iterMax)
	{
	  //searching directions
	 // s = g0/(*ls->hessian(model0)); //on réinitialise LU a chaque fois
	  ls->hessian(model0)->Solve(s,g0);
	  s = -1.*s;

	  descent = (s,g0);
	  // Cubic Line Search
	  model1 = ls->search(model0, s, descent, lambda);
	  g1 = *(ls->gradient(model1));
	  err = (Real)sqrt((g1,g1));
	  
	  if (isVerbose)
		cerr << "Iteration (" << iterNum << ") : "<<"current value of the objective function: "
			 <<ls->currentValue() << "\t current residue: "<< err << endl;
	  
	  appendResidue(err);	// residual
	  
	  g0=g1;
	  model0=model1;
	  
	  iterNum ++;
	}
  
  return(model1);
}

#endif

