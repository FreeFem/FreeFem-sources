#ifndef BFGS_HH
#define BFGS_HH

#include "Optima.hpp"
#include "defs.hpp"

// BFGS Quasi-Newton Optimization Method
/* At the present version, you must use the CubicLineSearch 
   procedure */

// Necessite une "vraie" classe Matrice
template <class LS>
class BFGS : public Optima<LS>
{
  typedef typename LS::Real Real;
  typedef typename LS::Param Param;
  typedef typename LS::Vect Vect;
  typedef typename LS::VMat VMat;
  typedef typename LS::Mat Mat;
  typedef list<Real> mlist;

public:
  
  BFGS(///pointer to the line-search object.
	   LS* ls, 
	   ///Maximum number of iterations
	   int iter, 
	   ///minimum accepted gradient at optimum solution
	   Real tol,
	   ///vebose or quiet
	   int verb=0);
  ~BFGS(){;}
  
  // the BFGS search starting from model0, returns an optimum model 
  Param optimizer(Param& model0);

};

template <class LS>
BFGS<LS>::BFGS(LS* p, int it, Real eps, int verb) 
  :Optima<LS>(verb)
{
  ls=p;
  iterMax 	= 	it;
  tol 		= 	eps;
  iterNum 	= 	0;
}

// Nécessite pour la classe MAT:
// - un constructeur de matrice avec nbre de lig, nbre de col
// - un constructeur de matrice avec la diagonale comme paramètre
// - le produit matrice vecteur
// - une méthode outProduct qui construit à partir de 2 vecteurs u et v la matrice A_{i,j}=u(i)v(j)
// - une méthode t()
// - un opérateur +
// - un opérateur = Real
// - un opérateur = Vect
template <class LS>
typename BFGS<LS>::Param BFGS<LS>::optimizer(Param& model0)
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
   
   //initial identical matrix for estimating inverse of the Hessian
   //Vect diag(n,1.);

   Mat H(n,n);H=0;diagonal(H)=1;
   Real d, dd, scale;
   Param model1(model0);
   Vect s(n),gamma(n),delta(n),g1(n);

   //searching directions
   s=0;
   s -= H*g0;
   descent = (s,g0);
   assert(s.max() <1e100);
   
   //cubic line search for a new model
   model1 = ls->search(model0, s, descent, lambda);
   g1 = *(ls->gradient(model1));
   err = (Real)sqrt((g1,g1));
   if (isVerbose)
	 cerr << "Iteration (" << iterNum << ") : "
		  <<"current value of the objective function: "
		  <<ls->currentValue() << "\t current residue: "<< err << endl;

   appendResidue(err);	// residual
   iterNum ++;

   Mat B(n,n);

   
   while (finalResidue() > tol && iterNum < iterMax)
   {
	 gamma = g1 - g0;
	 delta = model1 - model0;

	 //replace the searching direction with temporal storage
	 s = H*gamma;

	 //factor of the denominator
	 dd = (delta,gamma);

	 // Modif TONY
	 // Au départ, il y avait : if (dd < 0.00000001)
	 // Mais Adel ne fait pas ce test non plus...
	 if (Abs(dd)<1e-20)
	   {
		 // re-initialize the Hessian Matrix
		 // Il faut d'abord le mettre a zero : cf. Matrix.hpp
		 H = 0.; diagonal(H)=1.;
		
	   }
	 else
	   {
		 assert(dd);
		 d = 1./dd;
		 
		 scale = d*((gamma,s));
		 scale += 1;
		 scale *= d;
		 
		 // update the first term
		 H += scale*delta*delta.t();
		 
		 //update the second term
		 H -= d*s*delta.t();
		 H -= d*delta*s.t();
		 
		 //store the current model and gradient
		 g0 = g1;
		 model0 = model1;
	   }
	 s=0;
	 s -= H*g0;
	 descent = (s,g0);
	 model1 = ls->search(model0, s, descent, lambda);
	 g1 = *(ls->gradient(model1));
	 err = (Real)sqrt((g1,g1));
	 
	 if (isVerbose)
	   cerr << "Iteration (" << iterNum << ") : "<<"current value of the objective function: "
			<<ls->currentValue() << "\t current residue: "<< err << endl;

	 appendResidue(err);	// residual
	 iterNum ++;
   }
   
   return(model1);
}

#endif

