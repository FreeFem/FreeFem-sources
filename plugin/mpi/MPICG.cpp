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


#include "mpi.h"
#include  <iostream>
using namespace std;

#include "ff++.hpp"
#include "CGNL.hpp"
//#include "gmres.hpp"


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
	  if (  kprint < nbitermax )
	      cout << "CGNL:" <<iter <<  ",  ro = " << ro << " ||g||^2 = " << g2 << endl; 
	  if (g2 < reps2) { 
	      if (kprint< nbitermax )
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
template<> struct MPI_TYPE<long>      {static MPI_Datatype  TYPE(){return MPI_LONG;}};
template<> struct MPI_TYPE<int>      {static MPI_Datatype TYPE(){return MPI_INT;}};
template<> struct MPI_TYPE<double>    {static MPI_Datatype TYPE(){return MPI_DOUBLE;}};
template<> struct MPI_TYPE<char>    {static MPI_Datatype TYPE(){return MPI_BYTE;}};

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
   // if (verbosity>99) kprint=1;
    R ro=1;
    Rn g(n),h(n),Ah(n), & Cg(Ah);  // on utilise Ah pour stocke Cg  
    g = A*x;
    g -= b;  
    Cg = C*g; // gradient preconditionne 
    h =-Cg; 
    R g2 = ReduceSum1((Cg,g),commworld);
    if (g2 < 1e-30) 
      { if(kprint<=nbitermax)
	  cout << "GC  g^2 =" << g2 << " < 1.e-30  Nothing to do " << endl;
	  return 2;  }
    if (kprint<5 ) 
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
	  if (std::norm(hAh)<1e-100) ExecError("CG2: Matrix is not defined (/0), sorry ");
	  ro =  -gh*rop/hAh ; // ro optimal (produit scalaire usuel)
	  x += (ro-rop) *h;
	  g += (ro/rop) *Ah; // plus besoin de Ah, on utilise avec Cg optimisation
	  Cg = C*g;
	  R g2p=g2; 
	  g2 = ReduceSum1((Cg,g),commworld);
	  if ( ( (iter%kprint) == kprint-1)  )
	      cout << "CG:" <<iter <<  "  ro = " << ro << " ||g||^2 = " << g2 << endl; 
	  if (g2 < reps2) { 
	      if (kprint <= nbitermax)
		  cout << "CG converges " << iter <<  "  ro = " << ro << " ||g||^2 = " << g2 << endl; 
	      return 1;// ok 
          }
	  R gamma = g2/g2p;       
	  h *= gamma;
	  h -= Cg;  //  h = -Cg * gamma* h       
      }
    //    if (itermax <= nbitermax  )
	cout << "CG does'nt converge: " << nbitermax <<   " ||g||^2 = " << g2 << " reps2= " << reps2 << endl; 
    return 0; 
}


template < class Operator, class Vector, class Preconditioner,
class Matrix, class Real >
int 
GMRES_MPI(const Operator &A, Vector &x, const Vector &b,
      const Preconditioner &M, Matrix &H, int &m, int &max_iter,
      Real &tol,MPI_Comm * commworld,long verbosity)
{
    Real resid;
    int i, j = 1, k;
    Vector s(m+1), cs(m+1), sn(m+1), w,r,Ax;
    r=M*b;
    Real normb = sqrt(ReduceSum1((r,r),commworld));
    
    Ax=A * x;
    Ax=b-Ax;
    r = M*(Ax);
    Real beta = sqrt(ReduceSum1((r,r),commworld));
    
    if ( abs(normb) < 1.e-30)
	normb = 1;
    
    if ((resid = beta / normb) <= tol) {
	tol = resid;
	max_iter = 0;
	return 0;
    }
    
    Vector *v = new Vector[m+1];
    
    while (j <= max_iter) {
	v[0] = r / beta;    
	s = 0.0;
	s(0) = beta;
	
	for (i = 0; i < m && j <= max_iter; i++, j++) {
	    w = M*(Ax=A * v[i]);
	    for (k = 0; k <= i; k++) {
		H(k, i) = ReduceSum1((w, v[k]),commworld);
		w -= H(k, i) * v[k];
	    }
	    H(i+1, i) = sqrt(ReduceSum1((w,w),commworld));
	    v[i+1] = w  / H(i+1, i) ; 
	    
	    for (k = 0; k < i; k++)
		ApplyPlaneRotation(H(k,i), H(k+1,i), cs(k), sn(k));
	    
	    GeneratePlaneRotation(H(i,i), H(i+1,i), cs(i), sn(i));
	    ApplyPlaneRotation(H(i,i), H(i+1,i), cs(i), sn(i));
	    ApplyPlaneRotation(s(i), s(i+1), cs(i), sn(i));
	    if(verbosity>5 || (verbosity>2 && j%100==0) )
		cout << "GMRES: " << j << " " << abs(s(i+1)) << " " <<  normb << " " 
		<<  abs(s(i+1)) / normb << " < " << tol << endl;
	    
	    if ((resid = abs(s(i+1)) / normb) < tol) {
		if(verbosity)
		    cout << "GMRES converges: " << j << " " << abs(s(i+1)) << " " <<  normb << " " 
		    <<  abs(s(i+1)) / normb << " < " << tol << endl;
		
		Update(x, i, H, s, v);
		tol = resid;
		max_iter = j;
		delete [] v;
		return 0;
	    }
	}
	if(!(j <= max_iter)) break;
	Update(x, i-1 , H, s, v);
	Ax = A*x;
	Ax = b-Ax;
	
	r = M*(Ax);
	beta = sqrt(ReduceSum1((r,r),commworld));
	if(verbosity>4)
	    cout << "GMRES: restart" << j << " " << beta << " " <<  normb << " " 
	    <<  beta / normb << " < " << tol << endl;
	if ((resid = beta / normb) < tol) {
	    tol = resid;
	    max_iter = j;
	    delete [] v;
	    return 0;
	}
    }
    
    if(verbosity)
	cout << "WARNING: GMRES do not converges: " << j <<"/" << max_iter << ",  resid = " << resid 
	<< ", tol=  " << tol << ", normb "<< normb << endl;
    tol = resid;
    delete [] v;
    
    return 1;
}


template<class R>
class MPILinearCG : public OneOperator 
{ 
public:
    typedef KN<R> Kn;
    typedef KN_<R> Kn_;
    const int cas,CG;
    
    class MatF_O: RNM_VirtualMatrix<R> { public:
	Stack stack;
	mutable  Kn x;       
	C_F0 c_x;
	Kn *b;
	
	Expression  mat1,mat;
	typedef  typename RNM_VirtualMatrix<R>::plusAx plusAx;
	MatF_O(int n,Stack stk,const OneOperator * op,Kn *bb=0) 
	: RNM_VirtualMatrix<R>(n),stack(stk),
	x(n),c_x(CPValue(x)),b(bb),
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
        void addMatMul(const KN_<R> &  xx, KN_<R> & Ax) const
	//void addMatMul(const  Kn_  & xx, Kn_ & Ax) const
        {
	    ffassert(xx.N()==Ax.N());
	    x =xx;
	    Ax  += GetAny<Kn_>((*mat)(stack));
	    if(b && &Ax!=b) Ax += *b; // Ax -b => add b (not in cas of init. b c.a.d  &Ax == b 
	    WhereStackOfPtr2Free(stack)->clean();
	} 
	plusAx operator*(const Kn &  x) const {return plusAx(this,x);} 
        bool ChecknbLine(int n) const { return true;}
        bool ChecknbColumn(int m) const { return true;}
	
    };  
    
    
    class E_LCG: public E_F0mps { public:
	const int cas;// <0 => Nolinear
	const int CG; 
	static const int n_name_param=7;
	
	static basicAC_F0::name_and_type name_param[] ;
	
	
	Expression nargs[n_name_param];
	
	const OneOperator *A, *C; 
	Expression X,B;
	
	E_LCG(const basicAC_F0 & args,int cc,int gc) :cas(cc),CG(gc) 
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
		double eps = 1.0e-6;
		int nbitermax=  100;
		long verb = verbosity;
		
		pcommworld vcommworld=0;
		long dKrylov=50; 
		if (nargs[0]) eps= GetAny<double>((*nargs[0])(stack));
		if (nargs[1]) nbitermax = GetAny<long>((*nargs[1])(stack));
		if (nargs[3]) eps= *GetAny<double*>((*nargs[3])(stack));
		if (nargs[4]) vcommworld = GetAny<pcommworld>((*nargs[4])(stack));
		if (nargs[5])  dKrylov= GetAny<long>((*nargs[5])(stack));
                if (nargs[6]) verb=Abs(GetAny<long>((*nargs[6])(stack)));
		long gcverb=51L-Min(Abs(verb),50L);
		if(verb==0) gcverb = 1000000000;// no print 

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
		KN<R> * bbgmres =0;
		if ( !B && !CG) bbgmres=bb; // none zero if gmres without B 		
		MatF_O AA(n,stack,A,bbgmres);
		if(bbgmres ){
                     AA.addMatMul(*bbgmres,*bbgmres); // *bbgmres= AA* *bbgmres; // Ok Ax == b -> not translation of b .
		    *bbgmres = - *bbgmres;
		    if(verbosity>1) cout << "  ** GMRES set b =  -A(0);  : max=" << bbgmres->max() << " " << bbgmres->min()<<endl;
		}
		
		if(CG)
		  {
		    
		  
		 if (cas<0) {
		    if (C) 
		      { MatF_O CC(n,stack,C);
			ret = NLCG(AA,CC,x,nbitermax,eps, gcverb ,commworld );}
		    else 
		      ret = NLCG(AA,MatriceIdentite<R>(n),x,nbitermax,eps, gcverb ,commworld);
		}
		 else 
		    if (C) 
		      { MatF_O CC(n,stack,C);
			ret = ConjuguedGradient2(AA,CC,x,*bb,nbitermax,eps, gcverb ,commworld);}
		    else 
		      ret = ConjuguedGradient2(AA,MatriceIdentite<R>(n),x,*bb,nbitermax,eps, gcverb ,commworld);
		  }
		else {// GMRES 
		    
		    KNM<R> H(dKrylov+1,dKrylov+1);
		    int k=dKrylov;//,nn=n;
		    if (cas<0) {
			ErrorExec("NL GMRES:  to do! sorry ",1);
			/*       if (C) 
			 { MatF_O CC(n,stack,C);
			 ret = NLGMRES(AA,CC,x,nbitermax,eps, 51L-Min(Abs(verbosity),50L) );}
			 else 
			 ret = NLGMRES(AA,MatriceIdentite<R>(n),x,nbitermax,eps, 51L-Min(Abs(verbosity),50L));
			 ConjuguedGradient  */
		    }
		    else 
		      {
			if (C)
			  { MatF_O CC(n,stack,C); 
			      ret=GMRES_MPI(AA,(KN<R> &)x, *bb,CC,H,k,nbitermax,eps,commworld,verb);}
			else
			    ret=GMRES_MPI(AA,(KN<R> &)x, *bb,MatriceIdentite<R>(n),H,k,nbitermax,eps,commworld,verb);       
		      }
		    
		}
		
		

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
	return new E_LCG(args,cas,CG);}

    MPILinearCG() :   OneOperator(atype<long>(),
			       atype<Polymorphic*>(),
 			       atype<KN<R> *>(),atype<KN<R> *>()),cas(2),CG(1){}
    
    MPILinearCG(int cc,int CGG) :   OneOperator(atype<long>(),
			       atype<Polymorphic*>(),
			       atype<KN<R> *>(),atype<KN<R> *>()),cas(cc),CG(CGG){}

    MPILinearCG(int cc,int CGG,int ) :   OneOperator(atype<long>(),
						atype<Polymorphic*>(),
						atype<KN<R> *>()),cas(cc),CG(CGG){}
    
    MPILinearCG(int cc) :   OneOperator(atype<long>(),
				     atype<Polymorphic*>(),
				     atype<KN<R> *>()),cas(cc),CG(1){}
    
};




template<class R>
basicAC_F0::name_and_type  MPILinearCG<R>::E_LCG::name_param[]= {
    {   "eps", &typeid(double)  },
    {   "nbiter",&typeid(long) },
    {   "precon",&typeid(Polymorphic*)},
    {   "veps" ,  &typeid(double*) },
    { "comm", &typeid(pcommworld)} ,
    {   "dimKrylov", &typeid(long) },
    {   "verbosity", &typeid(long) }
};





/* --FH:   class Init { public:
    Init();
};

LOADINIT(Init);
*/ 
static void Load_Init()
{ 

    Global.Add("MPILinearCG","(",new MPILinearCG<R>()); // old form  with rhs (must be zer
    Global.Add("MPIAffineCG","(",new MPILinearCG<R>(1)); //  without right handsize
    Global.Add("MPILinearGMRES","(",new MPILinearCG<R>(0,0)); //  with  right handsize
    Global.Add("MPIAffineGMRES","(",new MPILinearCG<R>(0,0,0)); //  with  right handsize
    Global.Add("MPINLCG","(",new MPILinearCG<R>(-1)); //  without right handsize
    
}

 LOADFUNC(Load_Init)
