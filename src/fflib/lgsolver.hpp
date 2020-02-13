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
#ifndef lgsolver_hpp_
#define lgsolver_hpp_
#include "gmres.hpp"

typedef void *    pcommworld;

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

    template<class R>
    struct  RNM_VirtualMatrixTrans:public RNM_VirtualMatrix<R> { public:
        const RNM_VirtualMatrix<R> * A;
        int N,M;
        RNM_VirtualMatrixTrans(const RNM_VirtualMatrix<R> & AA): RNM_VirtualMatrix<R>(AA.M,AA.N),A(&AA)
        { N=AA.M; M=AA.N;
          //  cout <<"Assertion fail " <<  N << " " << M << " "<< AA.N << " " << AA.M << endl;
        }
        void addMatMul(const KN_<R> &  x, KN_<R> & y) const {A->addMatTransMul(x, y);}
        void addMatTransMul(const KN_<R> & x , KN_<R> & y) const  {A->addMatMul(x, y);}
        
        virtual bool ChecknbLine  (int n) const {return N==n;}
        virtual bool ChecknbColumn  (int m) const {return M==m;}
  /*      struct  plusAx { const RNM_VirtualMatrix * A; const KN_<R>   x;
            plusAx( const RNM_VirtualMatrix * B,const KN_<R> &  y) :A(B),x(y)
            { if(B) { ffassert(B->ChecknbColumn(y.N())); } }
        };
        
        plusAx operator*(const KN_<R> &  x) const {return plusAx(this,x);}
        
        struct  plusAtx { const RNM_VirtualMatrix * A; const KN_<R>   x;
            plusAtx( const RNM_VirtualMatrix * B,const KN_<R> &  y) :A(B),x(y)
            { if(B) { ffassert(B->ChecknbLine(y.N())); } } };
        
        struct  solveAxeqb { const RNM_VirtualMatrix * A; const KN_<R>   b;
            solveAxeqb( const RNM_VirtualMatrix * B,const KN_<R> &  y) :A(B),b(y)
            { if(B) { ffassert(B->ChecknbColumn(y.N())); } } };
        
        struct  solveAtxeqb { const RNM_VirtualMatrix * A; const KN_<R>   b;
            solveAtxeqb( const RNM_VirtualMatrix * B,const KN_<R> &  y) :A(B),b(y)
            { if(B) { ffassert(B->ChecknbColumn(y.N())); } } };
        */
         RNM_VirtualMatrixTrans(){}
    };
    

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
class MatC2R : public VirtualMatrix<double> { public:

  typedef typename VirtualMatrix<double>::plusAx plusAx;

  //  typedef  RNM_VirtualMatrix<complex<double> > M;

  const MM &m;

  MatC2R(const MM &mm):
      VirtualMatrix<double>(mm.N*2,mm.M*2),m(mm) {}

  void addMatMul(const  KN_<double>  & x, KN_<double> & Ax) const {

    R2C(Ax) += m*R2C(x);

  }

  plusAx operator*(const KN<double> &  x) const {return plusAx(this,x);}
  virtual bool ChecknbLine(int nn) const { return !n ||n==nn;}
  virtual bool ChecknbColumn(int mm) const { return !n ||m==mm;}


};

#ifdef REOMVE_CODE_V4

template<class R>
class SolveGCPrecon :   public MatriceCreuseOld<R>::VirtualSolver , public RNM_VirtualMatrix<R>{
  int n;
  int nbitermax;
  double eps;
  mutable double  epsr;
  const E_F0 * precon;
  KN<R> D1;   // pour le CL bloque (tgv)
  mutable KN<R> xx; 
  Expression xx_del, code_del; 
  Stack stack;
  typedef typename RNM_VirtualMatrix<R>::plusAx plusAx;
 
  public:
  SolveGCPrecon(const MatriceCreuseOld<R> &A,const OneOperator * C,Stack stk,int itmax, double epsilon=1e-6) :
    RNM_VirtualMatrix<R>(A.n),
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
    double tgv=std::norm(D1.linfty() )  ;
    if( std::norm(tgv) < 1e10) tgv=1e100;
    double tgv1 = 1./tgv; // Corrige fH 11 mai 2015 ...
     for (int i=0;i<n;i++)
       D1[i] = std::norm(D1[i]) < tgv  ? R(1.) : tgv1;
      
}
   void Solver(const MatriceCreuseOld<R> &a,KN_<R> &x,const KN_<R> &b) const  {
     epsr = (eps < 0) ? (epsr >0 ? -epsr : -eps ) : eps ;
    // cout << " epsr = " << epsr << endl;
     ConjuguedGradient<R,MatriceCreuseOld<R>,SolveGCPrecon<R>,StopGC<R> >(a,*this,b,x,nbitermax,epsr);
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
class SolveGMRESPrecon :   public MatriceCreuseOld<R>::VirtualSolver , public RNM_VirtualMatrix<R>{
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
  typedef typename RNM_VirtualMatrix<R>::plusAx plusAx;
  public:
  SolveGMRESPrecon(const MatriceCreuseOld<R> &A,const OneOperator * C,Stack stk,int dk=50,int itmax=0,double epsilon=1e-6) :
    RNM_VirtualMatrix<R>(A.n),
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
      double tgv = D1.linfty(); 
      if( tgv < 1e5) tgv = 1e100; // if no tgv remove .. 
     // if(verbosity>10 ) cout << "        in Precon GMRES, find tgv = " << tgv << endl;
      for (int i=0;i<n;i++)
        D1[i] = (std::abs(D1[i]) <  tgv ) ? R(1.)  : R() ; // remove the tgv ...
      
}
    void Solver(const MatriceCreuseOld<R> &a,KN_<R> &x,const KN_<R> &b) const  {
     epsr = (eps < 0) ? (epsr >0 ? -epsr : -eps ) : eps ;
      KNM<R> H(dKrylov+1,dKrylov+1);
      int k=dKrylov,nn=nbitermax;
      GMRES(a,(KN<R> &)x, (const KN<R> &)b,*this,H,k,nn,epsr,verbosity);

   }
    void SolverT(const MatriceCreuseOld<R> &a,KN_<R> &x,const KN_<R> &b) const  {
        RNM_VirtualMatrixTrans<R> at(a);
        epsr = (eps < 0) ? (epsr >0 ? -epsr : -eps ) : eps ;
        KNM<R> H(dKrylov+1,dKrylov+1);
        int k=dKrylov,nn=nbitermax;
        GMRES(at ,(KN<R> &)x, (const KN<R> &)b,*this,H,k,nn,epsr,verbosity);
    }
plusAx operator*(const KN_<R> &  x) const {return plusAx(this,x);} 


 void addMatMul(const KN_<R> & x, KN_<R> & Ax) const 
  { 
    assert(x.N()==Ax.N());
      
    xx=x;
   // cout << x[0] << "  ";
    xx=GetAny<KN_<R> >((*precon)(stack)); // xx value of the preoco 
    WhereStackOfPtr2Free(stack)->clean(); 

//    cout << (xx)[0] << "  " << endl;
 //   R dii;
  //   for (int i=0;i<n;i++) // take on value 
  //     Ax[i] += (D1[i]) ? xx[i] : x[i];//  remove dii ... mars 2011
    Ax += xx; 
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
class SolveGMRESDiag :   public MatriceCreuseOld<R>::VirtualSolver , public RNM_VirtualMatrix<R>{
  int n;
  int nbitermax;
  double eps;
  mutable double  epsr;
  int dKrilov;
  KN<R> D1;
  public:
  typedef typename RNM_VirtualMatrix<R>::plusAx plusAx;
  SolveGMRESDiag(const MatriceCreuseOld<R> &A,int nbk=50,int itmax=0,double epsilon=1e-6) :
     RNM_VirtualMatrix<R>(A.n),
    n(A.n),nbitermax(itmax?itmax: Max(100,n)),eps(epsilon),epsr(0),
    dKrilov(nbk) ,D1(n)
  { 
    R aii=0;
    A.getdiag(D1);
    for (int i=0;i<n;i++)
      D1[i] = (norm(aii=D1[i]) < 1e-20 ? R(1.) : R(1.)/aii);}

   void Solver(const MatriceCreuseOld<R> &a,KN_<R> &x,const KN_<R> &b) const  {
      epsr = (eps < 0) ? (epsr >0 ? -epsr : -eps ) : eps ;
      KNM<R> H(dKrilov+1,dKrilov+1);
      int k=dKrilov,nn=nbitermax;
      GMRES(a,(KN<R> &)x,(const KN<R> &)b,*this,H,k,nn,epsr,verbosity);

 }
    void SolverT(const MatriceCreuseOld<R> &a,KN_<R> &x,const KN_<R> &b) const  {
        RNM_VirtualMatrixTrans<R> at(a);
        epsr = (eps < 0) ? (epsr >0 ? -epsr : -eps ) : eps ;
        KNM<R> H(dKrilov+1,dKrilov+1);
        int k=dKrilov,nn=nbitermax;
        GMRES(at,(KN<R> &)x,(const KN<R> &)b,*this,H,k,nn,epsr,verbosity);
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
class SolveGMRESDiag<Complex> :   public MatriceMorseOld<Complex>::VirtualSolver , public RNM_VirtualMatrix<Complex>
{
   int n;
  int nbitermax;
  KN<Complex> D1;
  double eps;
  mutable double  epsr;
  int dKrilov;
  public:
  typedef  RNM_VirtualMatrix<Complex>::plusAx plusAx;
  SolveGMRESDiag(const MatriceMorseOld<Complex> &A,int nbk=50,int itmax=0,double epsilon=1e-6) :
    RNM_VirtualMatrix<Complex>(A.n),
    n(A.n),nbitermax(itmax?itmax: Max(100,n)),D1(n),eps(epsilon),epsr(0),
    dKrilov(nbk) 
  { 
    Complex aii=0;
    A.getdiag(D1);
    for (int i=0;i<n;i++)
      D1[i] = (std::norm(aii=D1[i]) < 1e-20 ? Complex(1.) : Complex(1.)/aii);}

   void Solver(const MatriceMorseOld<Complex> &a,KN_<Complex> &x,const KN_<Complex> &b) const  {
      epsr = (eps < 0) ? (epsr >0 ? -epsr : -eps ) : eps ;
          KNM<double> H(dKrilov+1,dKrilov+1);
      int k=dKrilov,nn=nbitermax;
      KN_<double> rx=C2R(x);
      const  KN_<double> rb=C2R(b);
      typedef MatC2R<MatriceMorseOld<Complex> > VA;
      typedef MatC2R<SolveGMRESDiag<Complex> > VC;
      VA AR(a);
      VC CR(*this);
      //int res=
      GMRES(AR,(KN<double> &)rx,(const KN<double> &)rb,CR,H,k,nn,epsr,verbosity);

 }

plusAx  operator*(const KN_<Complex> &  x) const {return plusAx(this,x);} 

    void SolverT(const MatriceCreuseOld<Complex> &a,KN_<Complex> &x,const KN_<Complex> &b) const  {
        epsr = (eps < 0) ? (epsr >0 ? -epsr : -eps ) : eps ;
        KNM<double> H(dKrilov+1,dKrilov+1);
        int k=dKrilov,nn=nbitermax;
        KN_<double> rx=C2R(x);
        const  KN_<double> rb=C2R(b);
        typedef MatC2R<MatriceCreuseOld<Complex> > VA;
        typedef MatC2R<SolveGMRESDiag<Complex> > VC;
        VA AR(a);
        VC CR(*this);
       RNM_VirtualMatrixTrans<double> ARt(AR);
         GMRES(ARt,(KN<double> &)rx,(const KN<double> &)rb,CR,H,k,nn,epsr,verbosity);
    }

 void addMatMul(const KN_<Complex> & x, KN_<Complex> & Ax) const 
  { 
     assert(x.N()==Ax.N());
   for (int i=0;i<n;i++) 
     Ax[i]+= D1[i]*x[i];}
     
  virtual bool ChecknbLine(int n) const { return true;}  
  virtual bool ChecknbColumn(int m) const { return true;} 

};


template<>
class SolveGMRESPrecon<Complex> :   public   MatriceCreuseOld<Complex>::VirtualSolver , public RNM_VirtualMatrix<Complex>
{
  public:
    typedef Complex R ;
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
  typedef  RNM_VirtualMatrix<Complex>::plusAx plusAx;
  public:
  SolveGMRESPrecon(const MatriceCreuseOld<Complex> &A,const OneOperator * C,Stack stk,int dk=50,int itmax=0,double epsilon=1e-6) :
    RNM_VirtualMatrix<Complex>(A.n),   
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
	D1[i] = std::norm(aii=D1[i])<1e-20 ? Complex(1.0) : Complex(1.)/aii;
      
}
   void Solver(const MatriceCreuseOld<Complex> &a,KN_<Complex> &x,const KN_<Complex> &b) const  {
     epsr = (eps < 0) ? (epsr >0 ? -epsr : -eps ) : eps ;
       KNM<double> H(dKrylov+1,dKrylov+1);
      int k=dKrylov,nn=nbitermax;
      KN_<double> rx=C2R(x);
      const  KN_<double> rb=C2R(b);
      typedef MatC2R<MatriceCreuseOld<Complex> > VA;
      typedef MatC2R<SolveGMRESPrecon<Complex> > VC;
      VA AR(a);
      VC CR(*this);
	GMRES(AR,(KN<double> &)rx,(const KN<double> &)rb,CR,H,k,nn,epsr,verbosity);
      
   }
    
    void SolverT(const MatriceCreuseOld<R> &a,KN_<R> &x,const KN_<R> &b) const  {
        epsr = (eps < 0) ? (epsr >0 ? -epsr : -eps ) : eps ;
        KNM<double> H(dKrylov+1,dKrylov+1);
        int k=dKrylov,nn=nbitermax;
        KN_<double> rx=C2R(x);
        const  KN_<double> rb=C2R(b);
        typedef MatC2R<MatriceCreuseOld<Complex> > VA;
        typedef MatC2R<SolveGMRESPrecon<Complex> > VC;
        VA AR(a);
        VC CR(*this);
        RNM_VirtualMatrixTrans<double> ARt(AR);
       GMRES(ARt,(KN<double> &)rx,(const KN<double> &)rb,CR,H,k,nn,epsr,verbosity);

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
   virtual bool ChecknbLine(int ) const { return true;}  
  virtual bool ChecknbColumn(int ) const { return true;} 

 
};     
#endif


#define LIST_NAME_PARM_MAT \
    {  "init", &typeid(bool)}, \
    {  "solver", &typeid(string*)}, \
    {  "eps", &typeid(double)  }, \
    {  "precon",&typeid(Polymorphic*)}, \
    {  "dimKrylov",&typeid(long)}, \
    {  "tgv",&typeid(double )}, \
    {  "factorize",&typeid(long)}, \
    {  "strategy",&typeid(long )}, \
    {  "tolpivot",&typeid(double )}, \
    {  "tolpivotsym",&typeid(double )}, \
    {  "nbiter", &typeid(long)}, \
    { "datafilename", &typeid(string*)} , \
    { "lparams",&typeid(KN_<long>)} , \
    { "dparams", &typeid(KN_<double>)},  \
    { "smap", &typeid(map<string,string>*)}, \
    { "permr", &typeid(KN_<long>)}, \
    { "permc", &typeid(KN_<long>)}, \
    { "scaler", &typeid(KN_<double>)}, \
    { "scalec", &typeid(KN_<double>)}, \
    { "sparams", &typeid(string*)}, \
    { "commworld", &typeid(pcommworld)}, \
    { "master", &typeid(long)}, \
    { "rinfo", &typeid(KN<double>*)}, \
    { "info", &typeid(KN<long>*)}, \
    { "kerneln", &typeid(  KNM<double> *)}, \
    { "kernelt", &typeid(  KNM<double> *)}, \
    { "kerneldim", &typeid(long*)}, \
    { "verb", &typeid(long)}, \
    { "x0", &typeid(bool)}, \
    { "veps", &typeid(double*)  }, \
    { "rightprecon", &typeid(bool)  }, \
    { "sym", &typeid(bool)  }, \
    { "positive", &typeid(bool)  }




const int NB_NAME_PARM_MAT =  24 +6+3  ;
    
    
template<class R>
inline void SetEnd_Data_Sparse_Solver(Stack stack,Data_Sparse_Solver & ds,Expression const *nargs ,int n_name_param,int syma=-1)
    {
        bool unset_eps=true;
        ds.initmat=true;
        ds.factorize=0;
	int kk = n_name_param-NB_NAME_PARM_MAT-1;
        if (nargs[++kk]) ds.initmat= ! GetAny<bool>((*nargs[kk])(stack));
	if (nargs[++kk]) ds.solver= * GetAny<string*>((*nargs[kk])(stack));
        ds.Init_sym_positive_var<R>(syma);//  set def value of sym and posi
	if (nargs[++kk]) ds.epsilon= GetAny<double>((*nargs[kk])(stack)),unset_eps=false;
	if (nargs[++kk])
	{// modif FH fev 2010 ...
	  const  Polymorphic * op=  dynamic_cast<const  Polymorphic *>(nargs[kk]);
            if(op)
            {
	      ds.precon = op->Find("(",ArrayOfaType(atype<KN<R>* >(),false)); // strange bug in g++ is R become a double
                ffassert(ds.precon);
            } // add miss 
	}
      
	if (nargs[++kk]) ds.NbSpace= GetAny<long>((*nargs[kk])(stack));
	if (nargs[++kk]) ds.tgv= GetAny<double>((*nargs[kk])(stack));
	if (nargs[++kk]) ds.factorize= GetAny<long>((*nargs[kk])(stack));	
	if (nargs[++kk]) ds.strategy = GetAny<long>((*nargs[kk])(stack)); 
	if (nargs[++kk]) ds.tol_pivot = GetAny<double>((*nargs[kk])(stack)); 
	if (nargs[++kk]) ds.tol_pivot_sym = GetAny<double>((*nargs[kk])(stack)); 
	if (nargs[++kk]) ds.itmax = GetAny<long>((*nargs[kk])(stack)); //  frev 2007 OK
	if (nargs[++kk]) ds.data_filename = *GetAny<string*>((*nargs[kk])(stack));
	if (nargs[++kk]) ds.lparams = GetAny<KN_<long> >((*nargs[kk])(stack));
	if (nargs[++kk]) ds.dparams = GetAny<KN_<double> >((*nargs[kk])(stack));
	if (nargs[++kk]) ds.smap = GetAny<MyMap<String,String> *>((*nargs[kk])(stack));
	if (nargs[++kk]) ds.perm_r = GetAny<KN_<long> >((*nargs[kk])(stack));
	if (nargs[++kk]) ds.perm_c = GetAny<KN_<long> >((*nargs[kk])(stack));
	if (nargs[++kk]) ds.scale_r = GetAny<KN_<double> >((*nargs[kk])(stack));
	if (nargs[++kk]) ds.scale_c = GetAny<KN_<double> >((*nargs[kk])(stack));
	if (nargs[++kk]) ds.sparams = *GetAny<string*>((*nargs[kk])(stack));
	if (nargs[++kk]) ds.commworld = GetAny<pcommworld>((*nargs[kk])(stack));
#ifdef VDATASPARSESOLVER
        if (nargs[++kk]) ds.master = GetAny<long>((*nargs[kk])(stack));
#else
        ++kk;
#endif
        // add FH nov 2015 ..
        if (nargs[++kk]) ds.rinfo = GetAny<KN<double>* >((*nargs[kk])(stack));
        if (nargs[++kk]) ds.info = GetAny<KN<long>* >((*nargs[kk])(stack));
         // add FH juin 2018   ..
        if (nargs[++kk]) ds.kerneln = GetAny< KNM<double>* >((*nargs[kk])(stack));
        if (nargs[++kk]) ds.kernelt = GetAny< KNM<double>* >((*nargs[kk])(stack));
        if (nargs[++kk]) ds.kerneldim = GetAny<long * >((*nargs[kk])(stack));
        if (nargs[++kk]) ds.verb = GetAny<long  >((*nargs[kk])(stack));
        if (nargs[++kk]) ds.x0 = GetAny<bool>((*nargs[kk])(stack));
        if (nargs[++kk]) ds.veps= GetAny<double*>((*nargs[kk])(stack));
        if( unset_eps && ds.veps) ds.epsilon = *ds.veps;//  if veps  and no def value  => veps def value of epsilon.
        if (nargs[++kk]) ds.rightprecon= GetAny<bool>((*nargs[kk])(stack));
        if (nargs[++kk]) ds.sym= GetAny<bool>((*nargs[kk])(stack));
        if (nargs[++kk]) ds.positive= GetAny<bool>((*nargs[kk])(stack));
        if(ds.solver == "")
        { // SET DEFAULT SOLVER TO HRE ... 
            if( ds.sym && ds.positive ) ds.solver=*def_solver_sym_dp;
            else if( ds.sym ) ds.solver=*def_solver_sym;
            else  ds.solver=*def_solver;
            if(verbosity>4) cout << "  **Warning: set default solver to " << ds.solver << endl;
        }

        ffassert(++kk == n_name_param);
    }
} // end of namespace Fem2D
#endif
