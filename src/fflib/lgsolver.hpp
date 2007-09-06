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
#include "gmres.hpp"
namespace  Fem2D {

// hack ---  F. Hecht -----

//  une fonction pour change un tableau de 

// complex en tableau de real 

//  Idee faire

// Une class qui tranforme une matrice complex en matric real

// et faire de transformateur de vecteur

inline KN_<double> C2R(KN_<complex<double> > & vc)

{ 

  assert(vc.step==1); 

  complex<double> * pc=vc; // pointeur du tableau

  double *pr = static_cast<double*>(static_cast<void*>(pc));

  return KN_<double>(pr,vc.N()*2);

}


inline const KN_<double> C2R(const KN_<complex<double> > & vc)

{ 

  assert(vc.step==1); 

  complex<double> * pc=vc; // pointeur du tableau

  double *pr = static_cast<double*>(static_cast<void*>(pc));

  return KN_<double>(pr,vc.N()*2);

}




inline KN_<complex<double> > R2C(KN_<double>  & vr)

{

  assert(vr.step==1 && vr.N() %2 ==0);

  double *pr =vr; // pointeur du tableau

   complex<double> * pc  = static_cast<complex<double>* >(static_cast<void*>(pr));

  return KN_<complex<double> >(pc,vr.N()/2);

}



inline const KN_<complex<double> > R2C(const KN_<double>  & vr)

{

  assert(vr.step==1 && vr.N() %2 ==0);

  double *pr =vr; // pointeur du tableau

   complex<double> * pc  = static_cast<complex<double>* >(static_cast<void*>(pr));

  return KN_<complex<double> >(pc,vr.N()/2);

}



//  une classe pour transforme une Matrice complex en Matrice  real 

// -----------------------------------------------------------------

template <class MM>
class MatC2R : public VirtualMatrice<double> { public:

  typedef typename VirtualMatrice<double>::plusAx plusAx;

  //  typedef  VirtualMatrice<complex<double> > M;

  const MM &m;

  MatC2R(const MM &mm):
      VirtualMatrice<double>(mm.N*2,mm.M*2),m(mm) {}

  void addMatMul(const  KN_<double>  & x, KN_<double> & Ax) const {

    R2C(Ax) += m*R2C(x);

  }

  plusAx operator*(const KN<double> &  x) const {return plusAx(this,x);}
  virtual bool ChecknbLine(int n) const { return !N ||n==N;}  
  virtual bool ChecknbColumn(int m) const { return !M ||m==M;} 


};



template<class R>
class SolveGCPrecon :   public MatriceMorse<R>::VirtualSolver , public VirtualMatrice<R>{
  int n;
  int nbitermax;
  double eps;
  mutable double  epsr;
  const E_F0 * precon;
  KN<R> D1;   // pour le CL bloque (tgv)
  mutable KN<R> xx; 
  Expression xx_del, code_del; 
  Stack stack;
  typedef typename VirtualMatrice<R>::plusAx plusAx;
 
  public:
  SolveGCPrecon(const MatriceMorse<R> &A,const OneOperator * C,Stack stk,int itmax, double epsilon=1e-6) : 
    VirtualMatrice<R>(A.n),
    n(A.n),nbitermax(itmax?itmax: Max(100,n)),eps(epsilon),epsr(0),precon(0),
    D1(n),xx(n),stack(stk)
{
      assert(C); 

      WhereStackOfPtr2Free(stack)=new StackOfPtr2Free(stack);// FH mars 2005   
      throwassert(A.sym());
      Type_Expr te_xx(CPValue(xx));
      xx_del=te_xx.second;
      C_F0 e_xx(te_xx); // 1 undelete pointer
      code_del= C->code(basicAC_F0_wa(e_xx));
      precon =  to<KN_<R> >(C_F0(code_del,*C));// 2 undelete pointer
      
      throwassert(precon);
      R aii;
      A.getdiag(D1);
     for (int i=0;i<n;i++)
       D1[i] = (norm(aii=D1[i]) < 1e-20 ? R(1.) : R(1.)/aii);
      
}
   void Solver(const MatriceMorse<R> &a,KN_<R> &x,const KN_<R> &b) const  {
     epsr = (eps < 0) ? (epsr >0 ? -epsr : -eps ) : eps ;
    // cout << " epsr = " << epsr << endl;
     ConjuguedGradient<R,MatriceMorse<R>,SolveGCPrecon<R> >(a,*this,b,x,nbitermax,epsr);
   }
plusAx operator*(const KN_<R> &  x) const {return plusAx(this,x);} 


 void addMatMul(const KN_<R> & x, KN_<R> & Ax) const 
  { 
    assert(x.N()==Ax.N());
      
    xx=x;
   // cout << x[0] << "  ";
    xx=GetAny<KN_<R> >((*precon)(stack));
    WhereStackOfPtr2Free(stack)->clean(); 
  //  cout << (xx)[0] << "  " << endl;
    R dii;
    for (int i=0;i<n;i++) 
       Ax[i] += ((dii=D1[i])==1.0) ? (xx)[i] : x[i]*dii;
  }

    
  ~SolveGCPrecon(){

  WhereStackOfPtr2Free(stack)->clean(); // FH mars 2005    
  delete  xx_del;
  delete  code_del;
// cout << "~SolveGCPrecon " << endl;
 }
  virtual bool ChecknbLine(int n) const { return true;}  
  virtual bool ChecknbColumn(int m) const { return true;} 

};     

template<class R>
class SolveGMRESPrecon :   public MatriceMorse<R>::VirtualSolver , public VirtualMatrice<R>{
  int n;
  int nbitermax;
  double eps;
  mutable double  epsr;
  const E_F0 * precon;
  KN<R> D1;  
  mutable KN<R> xx; 
  Expression xx_del, code_del;    
  Stack stack;
  int dKrylov; 
  typedef typename VirtualMatrice<R>::plusAx plusAx;
  public:
  SolveGMRESPrecon(const MatriceMorse<R> &A,const OneOperator * C,Stack stk,int dk=50,int itmax=0,double epsilon=1e-6) : 
    VirtualMatrice<R>(A.n),
    n(A.n),nbitermax(itmax?itmax: Max(100,n)),eps(epsilon),epsr(0),precon(0),
    D1(n),xx(n),
    stack(stk),dKrylov(dk)

{
      assert(C); 
      WhereStackOfPtr2Free(stack)=new StackOfPtr2Free(stack);// FH mars 2005   
      // C_F0 e_xx(CPValue(xx));
      //precon =  to<KN_<R> >(C_F0(C->code(basicAC_F0_wa(e_xx)),*C));
      Type_Expr te_xx(CPValue(xx));
      xx_del=te_xx.second;
      C_F0 e_xx(te_xx); // 1 undelete pointer
      code_del= C->code(basicAC_F0_wa(e_xx));
      precon =  to<KN_<R> >(C_F0(code_del,*C));// 2 undelete pointer
      
      throwassert(precon);
      R aii;
      A.getdiag(D1);
      for (int i=0;i<n;i++)
        D1[i] = (norm(aii=D1[i]) < 1e-20 ? R(1.) : R(1.)/aii);
      
}
   void Solver(const MatriceMorse<R> &a,KN_<R> &x,const KN_<R> &b) const  {
     epsr = (eps < 0) ? (epsr >0 ? -epsr : -eps ) : eps ;
    // cout << " epsr = " << epsr << endl;
  //   ConjuguedGradient<R,MatriceMorse<R>,SolveGCPrecon<R> >(a,*this,b,x,nbitermax,epsr);
      KNM<R> H(dKrylov+1,dKrylov+1);
      int k=dKrylov,nn=nbitermax;

      //int res=
      GMRES(a,(KN<R> &)x, (const KN<R> &)b,*this,H,k,nn,epsr);

   }
plusAx operator*(const KN_<R> &  x) const {return plusAx(this,x);} 


 void addMatMul(const KN_<R> & x, KN_<R> & Ax) const 
  { 
    assert(x.N()==Ax.N());
      
    xx=x;
   // cout << x[0] << "  ";
    xx=GetAny<KN_<R> >((*precon)(stack));
    WhereStackOfPtr2Free(stack)->clean(); 

//    cout << (xx)[0] << "  " << endl;
    R dii;
    for (int i=0;i<n;i++) 
       Ax[i] += ((dii=D1[i])==1.0) ? (xx)[i] : x[i]*dii;
  }

    
  ~SolveGMRESPrecon(){
  delete  xx_del;
  delete  code_del;
   WhereStackOfPtr2Free(stack)->clean(); // FH mars 2005    
  // cout << "~SolveGMRESPrecon; " << endl;
 }
  virtual bool ChecknbLine(int n) const { return true;}  
  virtual bool ChecknbColumn(int m) const { return true;} 
 
};     

template<class R>
class SolveGMRESDiag :   public MatriceMorse<R>::VirtualSolver , public VirtualMatrice<R>{
  int n;
  int nbitermax;
  double eps;
  mutable double  epsr;
  int dKrilov;
  KN<R> D1;
  public:
  typedef typename VirtualMatrice<R>::plusAx plusAx;
  SolveGMRESDiag(const MatriceMorse<R> &A,int nbk=50,int itmax=0,double epsilon=1e-6) : 
     VirtualMatrice<R>(A.n),
    n(A.n),nbitermax(itmax?itmax: Max(100,n)),eps(epsilon),epsr(0),
    dKrilov(nbk) ,D1(n)
  { 
    R aii=0;
    A.getdiag(D1);
    for (int i=0;i<n;i++)
      D1[i] = (norm(aii=D1[i]) < 1e-20 ? R(1.) : R(1.)/aii);}

   void Solver(const MatriceMorse<R> &a,KN_<R> &x,const KN_<R> &b) const  {
      epsr = (eps < 0) ? (epsr >0 ? -epsr : -eps ) : eps ;
    // cout << " epsr = " << epsr << endl;
   //  ConjuguedGradient<R,MatriceMorse<R>,SolveGCDiag<R> >(a,*this,b,x,nbitermax,epsr);
         KNM<R> H(dKrilov+1,dKrilov+1);
      int k=dKrilov,nn=nbitermax;
      //int res=
      GMRES(a,(KN<R> &)x,(const KN<R> &)b,*this,H,k,nn,epsr);

 }

plusAx  operator*(const KN_<R> &  x) const {return plusAx(this,x);} 


 void addMatMul(const KN_<R> & x, KN_<R> & Ax) const 
  { 
     assert(x.N()==Ax.N());
   for (int i=0;i<n;i++) 
     Ax[i]+= D1[i]*x[i];}
     
  virtual bool ChecknbLine(int n) const { return true;}  
  virtual bool ChecknbColumn(int m) const { return true;} 

};   

template<>
class SolveGMRESDiag<Complex> :   public MatriceMorse<Complex>::VirtualSolver , public VirtualMatrice<Complex>{
   int n;
  int nbitermax;
  KN<Complex> D1;
  double eps;
  mutable double  epsr;
  int dKrilov;
  public:
  typedef  VirtualMatrice<Complex>::plusAx plusAx;
  SolveGMRESDiag(const MatriceMorse<Complex> &A,int nbk=50,int itmax=0,double epsilon=1e-6) : 
    VirtualMatrice<Complex>(A.n),
    n(A.n),nbitermax(itmax?itmax: Max(100,n)),D1(n),eps(epsilon),epsr(0),
    dKrilov(nbk) 
  { 
    Complex aii=0;
    A.getdiag(D1);
    for (int i=0;i<n;i++)
      D1[i] = (norm(aii=D1[i]) < 1e-20 ? Complex(1.) : Complex(1.)/aii);}

   void Solver(const MatriceMorse<Complex> &a,KN_<Complex> &x,const KN_<Complex> &b) const  {
      epsr = (eps < 0) ? (epsr >0 ? -epsr : -eps ) : eps ;
    // cout << " epsr = " << epsr << endl;
   //  ConjuguedGradient<Complex,MatriceMorse<Complex>,SolveGCDiag<Complex> >(a,*this,b,x,nbitermax,epsr);
         KNM<double> H(dKrilov+1,dKrilov+1);
      int k=dKrilov,nn=nbitermax;
      KN_<double> rx=C2R(x);
      const  KN_<double> rb=C2R(b);
      typedef MatC2R<MatriceMorse<Complex> > VA;
      typedef MatC2R<SolveGMRESDiag<Complex> > VC;
      VA AR(a);
      VC CR(*this);
      //int res=
      GMRES(AR,(KN<double> &)rx,(const KN<double> &)rb,CR,H,k,nn,epsr);

 }

plusAx  operator*(const KN_<Complex> &  x) const {return plusAx(this,x);} 


 void addMatMul(const KN_<Complex> & x, KN_<Complex> & Ax) const 
  { 
     assert(x.N()==Ax.N());
   for (int i=0;i<n;i++) 
     Ax[i]+= D1[i]*x[i];}
     
  virtual bool ChecknbLine(int n) const { return true;}  
  virtual bool ChecknbColumn(int m) const { return true;} 

};   

template<>
class SolveGMRESPrecon<Complex> :   public MatriceMorse<Complex>::VirtualSolver , public VirtualMatrice<Complex>{
  public:
  int n;
  int nbitermax;
  Expression xx_del, code_del; 
  const E_F0 * precon;
  Stack stack;
  double eps;
  mutable double  epsr;
  int dKrylov; 
  KN<Complex> D1;  
  mutable KN<Complex> xx;  
  typedef  VirtualMatrice<Complex>::plusAx plusAx;
  public:
  SolveGMRESPrecon(const MatriceMorse<Complex> &A,const OneOperator * C,Stack stk,int dk=50,int itmax=0,double epsilon=1e-6) : 
    VirtualMatrice<Complex>(A.n),   
    n(A.n),nbitermax(itmax?itmax: Max(100,n)),
    xx_del(0),code_del(0),
    precon(0),stack(stk),eps(epsilon),epsr(0),dKrylov(dk),
    D1(n),xx(n)
{
      assert(C); 
      WhereStackOfPtr2Free(stack)=new StackOfPtr2Free(stack);// FH mars 2005   
      Type_Expr te_xx(CPValue(xx));
      xx_del=te_xx.second;
      C_F0 e_xx(te_xx); // 1 undelete pointer
      code_del= C->code(basicAC_F0_wa(e_xx));
      precon =  to<KN_<R> >(C_F0(code_del,*C));// 2 undelete pointer
      //C_F0 e_xx(CPValue(xx));
      //precon =  to<KN_<Complex> >(C_F0(C->code(basicAC_F0_wa(e_xx)),*C));
      
      throwassert(precon);
      Complex aii;
      A.getdiag(D1);
      for (int i=0;i<n;i++)
       D1[i] = norm(aii=D1[i])<1e-20 ? Complex(1.0) : Complex(1.)/aii;
      
}
   void Solver(const MatriceMorse<Complex> &a,KN_<Complex> &x,const KN_<Complex> &b) const  {
     epsr = (eps < 0) ? (epsr >0 ? -epsr : -eps ) : eps ;
    // cout << " epsr = " << epsr << endl;
  //   ConjuguedGradient<Complex,MatriceMorse<Complex>,SolveGCPrecon<Complex> >(a,*this,b,x,nbitermax,epsr);
      KNM<double> H(dKrylov+1,dKrylov+1);
      int k=dKrylov,nn=nbitermax;

      KN_<double> rx=C2R(x);
      const  KN_<double> rb=C2R(b);
      typedef MatC2R<MatriceMorse<Complex> > VA;
      typedef MatC2R<SolveGMRESPrecon<Complex> > VC;
      VA AR(a);
      VC CR(*this);
      //int res=
	GMRES(AR,(KN<double> &)rx,(const KN<double> &)rb,CR,H,k,nn,epsr);
      
      // assert(0); // a faire 
      //int res=GMRES(a,(KN<double> &)x, (const KN<double> &)b,*this,H,k,nn,epsr);

   }
plusAx operator*(const KN_<Complex> &  x) const {return plusAx(this,x);} 


 void addMatMul(const KN_<Complex> & x, KN_<Complex> & Ax) const 
  { 
    assert(x.N()==Ax.N());
      
    xx=x;
   // cout << x[0] << "  ";
    xx=GetAny<KN_<Complex> >((*precon)(stack));
       WhereStackOfPtr2Free(stack)->clean(); 

//    cout << (xx)[0] << "  " << endl;
    Complex dii;
    for (int i=0;i<n;i++) 
       Ax[i] += ((dii=D1[i])==1.0) ? (xx)[i] : x[i]*dii;
  }

    
  ~SolveGMRESPrecon(){
    WhereStackOfPtr2Free(stack)->clean(); // FH mars 2005    
    delete  xx_del;
    delete  code_del;
  
  // cout << "~SolveGMRESPrecon; " << endl;
 }
   virtual bool ChecknbLine(int n) const { return true;}  
  virtual bool ChecknbColumn(int m) const { return true;} 

 
};     


template<class R>
typename MatriceMorse<R>::VirtualSolver *
BuildSolverGMRES(const MatriceMorse<R> *A,int strategy,
		 double tgv, double eps, double tol_pivot,double tol_pivot_sym,
		 int NbSpace,int itmax , const void * precon, void * stack )
{
    typename MatriceMorse<R>::VirtualSolver * ret=0;
    if (precon)
	ret=new SolveGMRESPrecon<R>(*A,(const OneOperator *)precon,stack,NbSpace,itmax,eps);
    else 
	ret=new SolveGMRESDiag<R>(*A,NbSpace,itmax,eps);
    
    return ret;
}

template<class R>
typename MatriceMorse<R>::VirtualSolver *
BuildSolverCG(const MatriceMorse<R> *A,int strategy,
		 double tgv, double eps, double tol_pivot,double tol_pivot_sym,
		 int NbSpace,int itmax , const void * precon, void * stack )
{
    typename MatriceMorse<R>::VirtualSolver * ret=0;
    if (precon)
	ret=new SolveGCPrecon<R>(*A,(const OneOperator *)precon,stack,itmax,eps);
    else 
	ret=new SolveGCDiag<R>(*A,itmax,eps);
    
    return ret;
}


}

