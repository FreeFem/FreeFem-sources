// ORIG-DATE: 02/2009
// -*- Mode : c++ -*-
//
// SUMMARY  :  
// USAGE    : LGPL      
// ORG      : LJLL Universite Pierre et Marie Curie, Paris,  FRANCE 
// AUTHOR   : Jacques Morice
// E-MAIL   : jacques.morice@ann.jussieu.fr
//
//ff-c++-LIBRARY-dep:  mpi
//ff-c++-cpp-dep: 

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

 Thank to the ARN ()  FF2A3 grant
 ref:ANR-07-CIS7-002-01 
 */



#include  <iostream>
using namespace std;

#include "ff++.hpp"
#include "CGNL.hpp"
#include "mpi.h"

template<class R,class DJ,class P> 
int NLCG(const DJ & dJ,const P & C,KN_<R> &x,const int nbitermax, double &eps,long kprint,MPI_Comm * )
{
    //  -------------
    assert(&x && &dJ && &C);
    typedef KN<R> Rn;
    int n=x.N();
    
    R ro=1;
    Rn g(n),h(n),Ah(n), & Cg(Ah);  // on utilise Ah pour stocke Cg  
    g=dJ*x;// dJ(x,g);  
    Cg = C*g; // gradient preconditionne 
    h =-Cg; 
    R g2 = (Cg,g);
    if (g2 < 1e-30) 
      { if(kprint>1)
	  cout << "GCNL  g^2 =" << g2 << " < 1.e-30  Nothing to do " << endl;
	  return 2;  }
    if (kprint>5 ) 
	cout << " 0 GCNL  g^2 =" << g2 << endl;
    R reps2 =eps >0 ?  eps*eps*g2 : -eps; // epsilon relatif 
    eps = reps2;
    for (int iter=0;iter<=nbitermax;iter++)
      { 
	  ro = argmin(ro,dJ,x,h,g,Ah);
	  
	  Cg = C*g;
	  R g2p=g2; 
	  g2 = (Cg,g);
	  if (  kprint >1 )
	      cout << "CGNL:" <<iter <<  ",  ro = " << ro << " ||g||^2 = " << g2 << endl; 
	  if (g2 < reps2) { 
	      if (kprint )
		  cout << "CGNL converge: " << iter <<",  ro = " << ro << " ||g||^2 = " << g2 << endl; 
	      return 1;// ok 
	  }
	  R gamma = g2/g2p;       
	  h *= gamma;
	  h -= Cg;  //  h = -Cg * gamma* h       
      }
    if(verbosity)
	cout << " Non convergence of the NL cojugate gradient " <<endl;
    return 0; 
}
template<class T> struct MPI_TYPE {};
template<> struct MPI_TYPE<long>      {static const MPI_Datatype  TYPE(){return MPI_LONG;}};
template<> struct MPI_TYPE<int>      {static const MPI_Datatype TYPE(){return MPI_INT;}};
template<> struct MPI_TYPE<double>    {static const MPI_Datatype TYPE(){return MPI_DOUBLE;}};
template<> struct MPI_TYPE<char>    {static const MPI_Datatype TYPE(){return MPI_BYTE;}};

template<class R>
R ReduceSum1(R s,MPI_Comm * comm)
{
  R r=0;
  //  nt MPI_Allreduce( void *sendbuf, void *recbuf,  int count,
  //	    MPI_Datatype datatype, MPI_Op op, MPI_Comm comm )
  MPI_Allreduce( &s, &r, 1 ,MPI_TYPE<R>::TYPE(),   MPI_SUM,  *comm );
  return r; 
}

template<class R,class M,class P> 
int ConjuguedGradient2(const M & A,const P & C,KN_<R> &x,const KN_<R> &b,const int nbitermax, double &eps,long kprint,MPI_Comm * commworld)
{
    //  ConjuguedGradient2 affine A*x = 0 est toujours appele avec les condition aux limites 
    //  -------------
    throwassert(&x  && &A && &C);
    typedef KN<R> Rn;
    int n=x.N();
    if (verbosity>99) kprint=1;
    R ro=1;
    Rn g(n),h(n),Ah(n), & Cg(Ah);  // on utilise Ah pour stocke Cg  
    g = A*x;
    g -= b;  
    Cg = C*g; // gradient preconditionne 
    h =-Cg; 
    R g2 = ReduceSum1((Cg,g),commworld);
    if (g2 < 1e-30) 
      { if(verbosity>1)
	  cout << "GC  g^2 =" << g2 << " < 1.e-30  Nothing to do " << endl;
	  return 2;  }
    if (verbosity>5 ) 
	cout << " 0 GC  g^2 =" << g2 << endl;
    R reps2 =eps >0 ?  eps*eps*g2 : -eps; // epsilon relatif 
    eps = reps2;
    for (int iter=0;iter<=nbitermax;iter++)
      { 
	  R rop = ro; 
	  x += rop*h;      //   x+ rop*h  , g=Ax   (x old)
	  //       ((Ah = A*x - b) - g);
	  // Ah -= b;        //   Ax + rop*Ah = rop*Ah + g  =
	  // Ah -= g;         //   Ah*rop  
	  Ah = A*x;
	  Ah -= b;        //   Ax + rop*Ah = rop*Ah + g  =
	  Ah -= g;         //   Ah*rop  
	 
	  R hAh =ReduceSum1((h,Ah),commworld);
	  R gh = ReduceSum1((g,h),commworld);
	  if (norm(hAh)<1e-60) ExecError("CG2: Matrix is not defined (/0), sorry ");
	  ro =  -gh*rop/hAh ; // ro optimal (produit scalaire usuel)
	  x += (ro-rop) *h;
	  g += (ro/rop) *Ah; // plus besoin de Ah, on utilise avec Cg optimisation
	  Cg = C*g;
	  R g2p=g2; 
	  g2 = ReduceSum1((Cg,g),commworld);
	  if ( ( (iter%kprint) == kprint-1)  &&  verbosity >1 )
	      cout << "CG:" <<iter <<  "  ro = " << ro << " ||g||^2 = " << g2 << endl; 
	  if (g2 < reps2) { 
	      if (verbosity )
		  cout << "CG converges " << iter <<  "  ro = " << ro << " ||g||^2 = " << g2 << endl; 
	      return 1;// ok 
          }
	  R gamma = g2/g2p;       
	  h *= gamma;
	  h -= Cg;  //  h = -Cg * gamma* h       
      }
    if (verbosity )
	cout << "CG does'nt converge: " << nbitermax <<   " ||g||^2 = " << g2 << " reps2= " << reps2 << endl; 
    return 0; 
}

template<class R>
class MPILinearCG : public OneOperator 
{ 
public:
    typedef KN<R> Kn;
    typedef KN_<R> Kn_;
    const int cas;
    
    class MatF_O: VirtualMatrice<R> { public:
	Stack stack;
	mutable  Kn x;
	C_F0 c_x;
	
	Expression  mat1,mat;
	typedef  typename VirtualMatrice<R>::plusAx plusAx;
	MatF_O(int n,Stack stk,const OneOperator * op) 
	: VirtualMatrice<R>(n),stack(stk),
	x(n),c_x(CPValue(x)),
	mat1(op->code(basicAC_F0_wa(c_x))),
	mat( CastTo<Kn_>(C_F0(mat1,(aType)*op))) {
	    //ffassert(atype<Kn_ >() ==(aType) *op);
	    // WhereStackOfPtr2Free(stack)=new StackOfPtr2Free(stack);// FH mars 2005   
	    
	}
	~MatF_O() { 
	    // cout << " del MatF_O mat " << endl;
	    if(mat1 != mat) 
		delete mat;
	    delete mat1;
	    // cout << " del MatF_Ocx ..." <<  endl;
	    Expression zzz = c_x;
	    // cout << " zzz "<< zzz << endl;
	    delete zzz;
	    // WhereStackOfPtr2Free(stack)->clean(); // FH mars 2005 
	    
	}
	void addMatMul(const  Kn_  & xx, Kn_ & Ax) const { 
	    ffassert(xx.N()==Ax.N());
	    x =xx;
	    Ax  += GetAny<Kn_>((*mat)(stack));
	    WhereStackOfPtr2Free(stack)->clean();
	} 
	plusAx operator*(const Kn &  x) const {return plusAx(this,x);} 
	virtual bool ChecknbLine(int n) const { return true;}  
	virtual bool ChecknbColumn(int m) const { return true;} 
	
    };  
    
    
    class E_LCG: public E_F0mps { public:
	const int cas;// <0 => Nolinear
	static const int n_name_param=5;
	
	static basicAC_F0::name_and_type name_param[] ;
	
	
	Expression nargs[n_name_param];
	
	const OneOperator *A, *C; 
	Expression X,B;
	
	E_LCG(const basicAC_F0 & args,int cc) :cas(cc)
	{
	  args.SetNameParam(n_name_param,name_param,nargs);
	  {  const  Polymorphic * op=  dynamic_cast<const  Polymorphic *>(args[0].LeftValue());
	      ffassert(op);
	      A = op->Find("(",ArrayOfaType(atype<Kn* >(),false)); }
	  if (nargs[2]) 
	    {  const  Polymorphic * op=  dynamic_cast<const  Polymorphic *>(nargs[2]);
		ffassert(op); 
		C = op->Find("(",ArrayOfaType(atype<Kn* >(),false)); }
	  else  C =0;
	  X = to<Kn*>(args[1]);
	  if (args.size()>2)
	      B = to<Kn*>(args[2]);
	  else 
	      B=0;
	}
	
	virtual AnyType operator()(Stack stack)  const {
	    int ret=-1;
	    
	    // WhereStackOfPtr2Free(stack)=new StackOfPtr2Free(stack);// FH mars 2005   
	    try {
		Kn &x = *GetAny<Kn *>((*X)(stack));
		int n=x.N();
		MatF_O AA(n,stack,A);
		double eps = 1.0e-6;
		int nbitermax=  100;
		pcommworld vcommworld=0;
		if (nargs[0]) eps= GetAny<double>((*nargs[0])(stack));
		if (nargs[1]) nbitermax = GetAny<long>((*nargs[1])(stack));
		if (nargs[3]) eps= *GetAny<double*>((*nargs[3])(stack));
		if (nargs[4]) vcommworld = GetAny<pcommworld>((*nargs[4])(stack));
		MPI_Comm mpiCommWorld = MPI_COMM_WORLD;
		MPI_Comm * commworld= vcommworld ? (MPI_Comm *) vcommworld: & mpiCommWorld ;
		KN<R>  bzero(B?1:n); // const array zero
		bzero=R(); 
		KN<R> *bb=&bzero; 
		if (B) {
		    Kn &b = *GetAny<Kn *>((*B)(stack));
		    R p = (b,b);
		    if (p) 
		      {
			// ExecError("Sorry MPILinearCG work only with nul right hand side, so put the right hand in the function");
		      }
		    bb = &b;
		}
		if (cas<0) {
		    if (C) 
		      { MatF_O CC(n,stack,C);
			ret = NLCG(AA,CC,x,nbitermax,eps, 51L-Min(Abs(verbosity),50L),commworld );}
		    else 
		      ret = NLCG(AA,MatriceIdentite<R>(n),x,nbitermax,eps, 51L-Min(Abs(verbosity),50L),commworld);
		}
		else 
		    if (C) 
		      { MatF_O CC(n,stack,C);
			ret = ConjuguedGradient2(AA,CC,x,*bb,nbitermax,eps, 51L-Min(Abs(verbosity),50L),commworld);}
		    else 
		      ret = ConjuguedGradient2(AA,MatriceIdentite<R>(n),x,*bb,nbitermax,eps, 51L-Min(Abs(verbosity),50L),commworld);
		if( nargs[3]) *GetAny<double*>((*nargs[3])(stack)) = -(eps);
	    }
	    catch(...)
	  {
	    // WhereStackOfPtr2Free(stack)->clean(); // FH mars 2005 
	    throw;
	  }
	    // WhereStackOfPtr2Free(stack)->clean(); // FH mars 2005 
	    
	    return SetAny<long>(ret);
	    
	}  
	operator aType () const { return atype<long>();}         
	
    };
    
    E_F0 * code(const basicAC_F0 & args) const {
	return new E_LCG(args,cas);}
    MPILinearCG() :   OneOperator(atype<long>(),
			       atype<Polymorphic*>(),
			       atype<KN<R> *>(),atype<KN<R> *>()),cas(2){}
    MPILinearCG(int cc) :   OneOperator(atype<long>(),
				     atype<Polymorphic*>(),
				     atype<KN<R> *>()),cas(cc){}
    
};


template<class R>
basicAC_F0::name_and_type  MPILinearCG<R>::E_LCG::name_param[]= {
    {   "eps", &typeid(double)  },
    {   "nbiter",&typeid(long) },
    {   "precon",&typeid(Polymorphic*)},
    {   "veps" ,  &typeid(double*) },
    { "comm", &typeid(pcommworld)} 
};


class Init { public:
    Init();
};

Init init;
Init::Init()
{ 

    Global.Add("MPILinearCG","(",new MPILinearCG<R>()); // old form  with rhs (must be zer
    //Global.Add("MPILinearGMRES","(",new LinearGMRES<R>()); // old form  with rhs (must be zer
    // Global.Add("LinearGMRES","(",new LinearGMRES<R>(1)); // old form  with rhs (must be zer
    Global.Add("MPILinearCG","(",new MPILinearCG<R>(1)); //  without right handsize
    Global.Add("MPINLCG","(",new MPILinearCG<R>(-1)); //  without right handsize
    
}

