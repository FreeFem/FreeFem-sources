//======================================================================
// Definition of the BrentLineSearch class
// Brent's line search algorithm
// author: Wenceslau Gouveia, Adapted from Numerical Recipes in C
// Modified to fit into new classes.  H. Lydia Deng, 02/21/94, 03/15/94
// Tony : cf. page 301 section 10.2
//========================================================================
// .B BrentLineSearch()
// This routine, inpired by the Numerical Recipes book, performs a unidimensional search
// for the minimum of the objective function along a specified direction. The minimum is
// at first bracket using the Golden search procedure. After bracketing the Brent's 
// algorithm is used to isolate the minimum to a fractional precision of about the specified 
// tolerance.
// 
// .SECTION Description
// Public Operations
// Constructors: 
//		BrentLineSearch (ObjectiveFunction *f, int iter);
// 		Here:
//			f: Define the objective function
//			iter: Maximum number of iterations
// Methods:
//		model <> BrentLineSearch::search(Model<double>& model0, Vector<double>& direction, 
//						  	                 double tol, double delta)
//		Here:
//			model0:  Initial model to initiate the bracketing procedure
// 			direction: Vector that defines the direction of the line search
//			tol: The minimum is within the returned value +/- tol
//			delta: Used in the bracketing procedure. The initial interval
//				  for the bracketing is from 0 to delta * STEP_MAX, where
//				  STEP_MAX is hard coded to 5.
//
//			The sought minimum is returned by the function.
//
// .SECTION Caveats
// This line search was not thoroughly tested. The CubicLineSearch procedure, that requires
// certain derivative information on the objective function (that can be provided by numerical
// methods) has demonstrated to be a  more efficient line search procedure.

// Cela me semble louche... cf. les commentaires BIZARRE !
// Ceci dit, le premiers tests semblent OK...

#include "defs.hpp"
#include <stdio.h>


#ifndef BRENT_LINE_SEARCH_HH
#define BRENT_LINE_SEARCH_HH

#include "LineSearch.hpp"

#define CGOLD 		.3819660
#define ZEPS 		1.e-10
#define GOLD 		1.618034
#define STEP_MAX	5
#define GLIMIT 		100.

/*
   BrentLineSearch class was inpired by the Numerical Recipes book, 
   performs an unidimensional search for the minimum of the objective 
   function along a specified direction. The minimum is at first
   bracket using the Golden search procedure. After bracketing 
   the Brent's algorithm is used to isolate the minimum to a 
   fractional precision of about the specified tolerance.

   This line search was not thoroughly tested. The 
   CubicLineSearch procedure, that requires certain 
   derivative information on the objective function
   (that can be provided by numerical methods) has
   demonstrated to be a  more efficient line search procedure.

*/

template <class LS>
class BrentLineSearch : public LS
{
  typedef typename LS::Real Real;
  typedef typename LS::Param Param;
  typedef typename LS::Vect Vect;
  typedef typename LS::VMat VMat;
  typedef typename LS::Mat Mat;
  typedef typename LS::NRJ NRJ;
   public:				

   BrentLineSearch(NRJ* f, int iter);
   ~BrentLineSearch();

  // Implementation of the Brent Search
  // search for minimum model along a 1-D direction
  Param search(///initial model to initiate the bracketing procedure
			Param& m0, 
			///the direction of the line search
			Vect& direction, 
			///the minimum is within the returned value +/- tol
			Real tol, 
			///a parameter used in the bracketing procedure
			double delta);
   /*The initial interval for the bracketing is from $0$ to 
   $delta \times STEP\_MAX$, where $STEP\_MAX$ is hard coded to $5$.
   The sought minimum is returned by the function.  */
   //@ManMemo: search for minimum model along a 1-D direction
};

template <class LS>
BrentLineSearch<LS>::BrentLineSearch(NRJ * f, int it) 
: LS(f)
{
  iterMax 	= 	it;
}

template <class LS>
BrentLineSearch<LS>::~BrentLineSearch()
{;}

// Code for the BrentBrent line search
template <class LS>
typename BrentLineSearch<LS>::Param BrentLineSearch<LS>::search(Param& model0,Vect& direction, Real tol, double delta)
{
    iterNum	=	0;
	KN<double> steps(3);			// brackets
	Vect of_values(3);			// OF evaluated inside bracket
	Real dum;					// auxiliary quantity
	Real r, q, ulim;				// auxiliary quantity
	Real u, fu;					// define the new bracket limit
	
	// Variables related to Brent's algorithm
	Real d, fv, fw, fx, p;			// auxiliary variables
	double tol1, tol2, v, w, xm;			// auxiliary variables
	double x, e = 0.,etemp;				// auxiliary variables
	int i;						// counter

	// Ajout Tony
	d=0; // Sinon, il n'est pas initialise : BIZARRE !
	
	
	// Beggining of the bracketing stage

	steps[0] = 0.;
	steps[1] = STEP_MAX * delta;

	of_values[0]= nrj->getVal(model0);
	of_values[1] =nrj->getVal(update(model0,1,steps[1],direction));
	iterNum += 2;

	if (of_values[1] > of_values[0])
	{
	   dum = steps[0]; steps[0] = steps[1]; steps[1] = dum;
	   dum = of_values[0]; of_values[0] = of_values[1]; of_values[1] = dum;
	}	

	steps[2] = steps[1] + GOLD * (steps[1] - steps[0]);
	of_values[2] = nrj->getVal(update(model0,1,steps[2],direction));
	iterNum++;
	while (of_values[1] > of_values[2])
	{
		r = (steps[1] - steps[0]) * (of_values[1] - of_values[2]);
		q = (steps[1] - steps[2]) * (of_values[1] - of_values[0]);
		u = steps[1] - ((steps[1] - steps[2]) * q - 
				(steps[1] - steps[0]) * r) / 
		    (2. * Abs(Max(Abs(q - r), TINY)) / Sgn(q - r));
		ulim = steps[1] + GLIMIT * (steps[2] - steps[1]);

		if ((steps[1] - u) * (u - steps[2]) > 0.)
		{
			fu = nrj->getVal(update(model0,1,u,direction));
			iterNum++;
			if (fu < of_values[2])
			{
				steps[0] = steps[1];
				steps[1] = u;
				of_values[0] = of_values[1];
				of_values[1] = fu;
				break;
			}
			else if (fu > of_values[1])
			{
				steps[2] = u;
				of_values[2] = fu;
				break;
			}

			u = steps[2] + GOLD * (steps[2] - steps[1]);
			fu = nrj->getVal(update(model0,1,u,direction));
			iterNum++;

		}
		else if ((steps[2] - u) * (u - ulim) > 0.)
		{
			fu = nrj->getVal(update(model0,1,u,direction));
			iterNum++;
			if (fu < of_values[2])
			{
				steps[1] = steps[2];
				steps[2] = u;
				u = steps[2] + GOLD * (steps[2] - steps[1]);

				of_values[1] = of_values[2];
				of_values[2] = fu;
				fu = nrj->getVal(update(model0,1,u,direction));
				iterNum ++;
			}
		}
		else if ((u - ulim) * (ulim * steps[2]) >= 0.)
		{
			u = ulim;
			fu = nrj->getVal(update(model0,1,u,direction));
			iterNum ++;
		}
		else
		{
			u = steps[2] + GOLD * (steps[2] - steps[1]);
			fu = nrj->getVal(update(model0,1,u,direction));
			iterNum ++;
		}

		steps[0] = steps[1]; steps[1] = steps[2]; steps[2] = u;
		of_values[0] = of_values[1]; of_values[1] = of_values[2]; 
		of_values[2] = fu;

	}
	
	// Sorting STEPS in ascending order
	for (i = 0; i < 2; i++)
	{
		if (steps[0] > steps[i+1])
		{
			dum = steps[0];
			steps[0] = steps[i+1];
			steps[i+1] = dum;

			dum = of_values[0];
			of_values[0] = of_values[i+1];
			of_values[i+1] = dum;
		}
	}

	if (steps[1] > steps[2])
	{
		dum = steps[1];
		steps[1] = steps[2];
		steps[2] = dum;

		dum = of_values[1];
		of_values[1] = of_values[2];
		of_values[2] = dum;
	}

	// The line minimization will be performed now using as
	// bracket the 3 first steps given in vector steps 
	// The algorithm is due to Brent

	// initializations
	x = w = v = steps[1];
	fw = fv = fx = of_values[1];

	for (; iterNum <= iterMax; iterNum++)
	{
		xm = .5 * (steps[0] + steps[2]);
		tol1 = tol * Abs(x) + ZEPS;
		tol2 = 2.0 * tol1;
		if (Abs(x - xm) <= (tol2 - .5 * (steps[2] - steps[0])))
		{
			value = fx;
			Param new_model(update(model0,1,x,direction));
			
			return new_model;
		}
						// minimum along a line
		if (Abs(e) > tol1)
		{
			r = (x - w) * (fx - fv);
			q = (x - v) * (fx - fw);
			p = (x - v) * q - (x - w) * r;
			q = 2.0 * (q - r);
			// Tony
			// ??? ERREUR ???
			// dans le numerical recipe, c'est p=-p...
			// BIZARRE !
			if (q > 0.){
			  cerr<<"BIZARRE BrentLS"<<endl;
			  p = -q;
			}
			
			q = Abs(q);
			etemp = e;
			e = d;
			if (Abs(p) >= Abs(.5 * q * etemp) ||
			    p <= q * (steps[0] - x) ||
			    p >= q * (steps[2] - x))
			{	// parabolic fit
				if (x >= xm)
					e = steps[0] - x;
				else
					e = steps[2] - x;
				d = CGOLD * e;
			}
			else
			{
				d = p / q;
				u = x + d;
				if (u - steps[0] < tol2 ||
				    steps[2] - u < tol2)
					d = Abs(tol1) / Sgn(xm-x);
			}
		}
		else
		{
			if (x >= xm)
				e = steps[0] - x;
			else
				e = steps[2] - x;

			d = CGOLD * e;
		}

		if (Abs(d) >= tol1)
			u = x + d;
		else
			u = x + Abs(tol1) / Sgn(d);

		fu  = nrj->getVal(update(model0,1,u,direction));

		if (fu <= fx)
		{
			if (u >= x)
				steps[0] = x;
			else
				steps[2] = x;
			v = w; w = x; x = u;
			fv = fw; fw = fx; fx = fu;
		}
		else
		{
			if (u < x)
				steps[0] = u;
			else
				steps[2] = u;

			if (fu <= fw || w == x)
			{
				v = w; w = u; fv = fw; fw = fu;
			}

			else if (fu <= fw || v == x || v == w)

			{
				v = u; fv = fu;
			}
		}
	}
	
    appendSearchNumber();
	//    cout << " Maximum number of iterations reached " << endl;


	Param new_model(update(model0,1,x,direction));
	value = nrj->getVal(new_model);			
	return new_model;
}

#endif
