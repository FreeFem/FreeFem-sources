//ff-c++-cpp-dep: newuoa.f
// -*- Mode : c++ -*-
//
// SUMMARY  :      
// USAGE    :        
// ORG      : 
// AUTHOR   : Frederic Hecht
// E-MAIL   : hecht@ann.jussieu.fr
//

/*
 
 This file is part of Freefem++
 
 Freefem++ is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation; either version 2.1 of the License, or
 (at your option) any later version.
 
 Freefem++  is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.
 
 You should have received a copy of the GNU Lesser General Public License
 along with Freefem++; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */


#include  <iostream>
#include  <cfloat>
using namespace std;
#include "error.hpp"
#include "AFunction.hpp"
#include "rgraph.hpp"
#include "RNM.hpp"
#include "MatriceCreuse_tpl.hpp"
#include "Mesh3dn.hpp"
#include "MeshPoint.hpp"
#include "lgfem.hpp"
#include "lgmesh3.hpp"
#include "lgsolver.hpp"
#include "problem.hpp"

typedef int integer;
typedef int logical;


typedef void  (*typecalfunc)( integer *, double *, double *f, void * );

#define F77newuoa newuoa_

extern "C" {
double F77newuoa(integer *N, integer *NPT, double *x , double * rhob, double *rhog,
		     integer *iprint, integer *maxfun,
		     double *w, void * iwf,  typecalfunc  calfun);
}


void   calfun( integer * n, double * x, double *f, void * t);

 
//template<class R>
extern Block *currentblock;

typedef double R;
void   calfun( integer * n, double * x, double *f, void * t);
class OptimNewoa : public OneOperator 
{
public:
  typedef KN<R> Kn;
  typedef KN_<R> Kn_;
   const int cas;


class ffcalfunc {  //   to call the freefem function .. J 
    public:
	Stack stack;
	Expression JJ,theparame;
    
      ffcalfunc(Stack s,Expression JJJ,Expression epar)
        : stack(s),JJ(JJJ), theparame(epar) {}
    
	double J(Kn_  x) const 
	{
	     KN<double> *p=GetAny<KN<double> *>( (*theparame)(stack) );
	    *p=x;
	    double  ret= GetAny<R>( (*JJ)(stack));
	    WhereStackOfPtr2Free(stack)->clean();
	return  ret; }
	
    }; 
  

  class E_newoa: public E_F0mps { public:
    const int cas;
    static basicAC_F0::name_and_type name_param[] ;
    static const int n_name_param =4;
    Expression nargs[n_name_param];
    Expression X;
    C_F0 inittheparam,theparam,closetheparam; 
    Expression JJ;
    long arg(int i,Stack stack,long a) const{ return nargs[i] ? GetAny<long>( (*nargs[i])(stack) ): a;}
    R arg(int i,Stack stack,R a) const{ return nargs[i] ? GetAny<R>( (*nargs[i])(stack) ): a;}
      
 
      
    E_newoa(const basicAC_F0 & args,int cc) :
      cas(cc)
    {
      int nbj= args.size()-1;
      Block::open(currentblock); // make a new block to 
      X = to<Kn*>(args[nbj]);
      C_F0 X_n(args[nbj],"n");
      //  the expression to init the theparam of all 
      inittheparam = currentblock->NewVar<LocalVariable>("the parameter",atype<KN<R> *>(),X_n);
      theparam = currentblock->Find("the parameter"); //  the expression for the parameter
      args.SetNameParam(n_name_param,name_param,nargs);
      const  Polymorphic * opJ=0;
      if (nbj>0)
	{
	  opJ=  dynamic_cast<const  Polymorphic *>(args[0].LeftValue());
	  assert(opJ);
	}      
      JJ= to<R>(C_F0(opJ,"(",theparam));
      closetheparam=currentblock->close(currentblock);   // the cleanning block expression 
    }
    
    virtual AnyType operator()(Stack stack)  const
    {
      double cost = 1e100;
      WhereStackOfPtr2Free(stack)=new StackOfPtr2Free(stack);// FH mars 2005 
      Kn &x = *GetAny<Kn *>((*X)(stack));	
      long n=x.N();	
      double rhobeg=arg(0,stack,1E-6); // not used ....
      double rhoend=arg(1,stack,2.); // not used ....
      long maxfun=arg(2,stack,1000L); // bof bof 
      long npt=arg(3,stack,n*2L+1L); // bof bof 
      long iprint = verbosity;	
      ffcalfunc ffJ(stack,JJ,theparam);
      int lw =   (npt+13)*(npt+n)+3*n*(n+3)/2;
      KN<double> w(lw);
      integer N=n,NPT=npt,IPRINT=iprint,MAXFUN=maxfun;
      cost= F77newuoa(&N,&NPT,(double *)x,&rhobeg,&rhoend,&IPRINT,&MAXFUN,(double *)w,(void *) &ffJ, calfun);
      closetheparam.eval(stack); // clean memory 
      WhereStackOfPtr2Free(stack)->clean(); // FH mars 2005 
      return cost; //SetAny<long>(0);  Modif FH  july 2005       
    }    
    
    
    operator aType () const { return atype<double>();}         
    
  };
  
  
  
  E_F0 * code(const basicAC_F0 & args) const {
    return new E_newoa(args,cas);}
  
  OptimNewoa(int c) :   OneOperator(atype<double>(),
				    atype<Polymorphic*>(),
				    atype<KN<R> *>()),cas(c){}
      
};

basicAC_F0::name_and_type  OptimNewoa::E_newoa::name_param[]= 
  {
    {  "rhobeg", &typeid(double)  },
    {  "rhoend", &typeid(double)  },
    { "maxfun",&typeid(long) },
    { "npt",&typeid(long) }
  };

void   calfun( integer * n, double * x, double *f, void * t)
{
  OptimNewoa::ffcalfunc * tt=static_cast<OptimNewoa::ffcalfunc *>(t); 
  *f=tt->J(KN_<double>(x,*n));
  if(verbosity>20)  cout << " F= " << * f << endl;
}

class Init { public:
  Init();
};

LOADINIT(Init);  //  une variable globale qui serat construite  au chargement dynamique 

Init::Init()  // le constructeur qui ajoute la fonction "splitmesh3"  a freefem++ 
{
  Global.Add("newuoa","(",new OptimNewoa(1));  //  j + dJ

}
