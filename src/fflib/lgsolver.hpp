#include "gmres.hpp"
namespace  Fem2D {
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
 
  public:
  SolveGCPrecon(const MatriceMorse<R> &A,const OneOperator * C,Stack stk,double epsilon=1e-6) : 
    n(A.n),nbitermax(Max(n,100)),precon(0),stack(stk),eps(epsilon),epsr(0),
    D1(n),xx(n)
{
      assert(C); 
      throwassert(A.sym());
      Type_Expr te_xx(CPValue(xx));
      xx_del=te_xx.second;
      C_F0 e_xx(te_xx); // 1 undelete pointer
      code_del= C->code(basicAC_F0_wa(e_xx));
      precon =  to<KN<R> *>(C_F0(code_del,*C));// 2 undelete pointer
      
      throwassert(precon);
      R aii;
      for (int i=0;i<n;i++)
       D1[i] = Abs(aii=A(i,i))<1e10 ? 1.0 : 1/aii;
      
}
   void Solver(const MatriceMorse<R> &a,KN_<R> &x,const KN_<R> &b) const  {
     epsr = (eps < 0) ? (epsr >0 ? -epsr : -eps ) : eps ;
    // cout << " epsr = " << epsr << endl;
     ConjuguedGradient<R,MatriceMorse<R>,SolveGCPrecon<R> >(a,*this,b,x,nbitermax,epsr);
   }
typename VirtualMatrice<R>::plusAx operator*(const KN_<R> &  x) const {return plusAx(this,x);} 


 void addMatMul(const KN_<R> & x, KN_<R> & Ax) const 
  { 
    assert(x.N()==Ax.N());
      
    xx=x;
   // cout << x[0] << "  ";
    xx=*GetAny<KN<R> *>((*precon)(stack));
  //  cout << (xx)[0] << "  " << endl;
    R dii;
    for (int i=0;i<n;i++) 
       Ax[i] += ((dii=D1[i])==1.0) ? (xx)[i] : x[i]*dii;
  }

    
  ~SolveGCPrecon(){
  delete  xx_del;
  delete  code_del;
// cout << "~SolveGCPrecon " << endl;
 }
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
  Stack stack;
  int dKrylov; 
  public:
  SolveGMRESPrecon(const MatriceMorse<R> &A,const OneOperator * C,Stack stk,int dk=50,int itmax=0,double epsilon=1e-6) : 
    n(A.n),nbitermax(itmax?itmax: Max(100,n)),precon(0),stack(stk),eps(epsilon),epsr(0),dKrylov(dk),
    D1(n),xx(n)
{
      assert(C); 
      C_F0 e_xx(CPValue(xx));
      precon =  to<KN<R> *>(C_F0(C->code(basicAC_F0_wa(e_xx)),*C));
      
      throwassert(precon);
      R aii;
      for (int i=0;i<n;i++)
       D1[i] = Abs(aii=A(i,i))<1e10 ? 1.0 : 1/aii;
      
}
   void Solver(const MatriceMorse<R> &a,KN_<R> &x,const KN_<R> &b) const  {
     epsr = (eps < 0) ? (epsr >0 ? -epsr : -eps ) : eps ;
    // cout << " epsr = " << epsr << endl;
  //   ConjuguedGradient<R,MatriceMorse<R>,SolveGCPrecon<R> >(a,*this,b,x,nbitermax,epsr);
      KNM<R> H(dKrylov+1,dKrylov+1);
      int k=dKrylov,nn=nbitermax;

      int res=GMRES(a,(KN<R> &)x, (const KN<R> &)b,*this,H,k,nn,epsr);

   }
typename VirtualMatrice<R>::plusAx operator*(const KN_<R> &  x) const {return plusAx(this,x);} 


 void addMatMul(const KN_<R> & x, KN_<R> & Ax) const 
  { 
    assert(x.N()==Ax.N());
      
    xx=x;
   // cout << x[0] << "  ";
    xx=*GetAny<KN<R> *>((*precon)(stack));
//    cout << (xx)[0] << "  " << endl;
    R dii;
    for (int i=0;i<n;i++) 
       Ax[i] += ((dii=D1[i])==1.0) ? (xx)[i] : x[i]*dii;
  }

    
  ~SolveGMRESPrecon(){
  // cout << "~SolveGMRESPrecon; " << endl;
 }
};     

template<class R>
class SolveGMRESDiag :   public MatriceMorse<R>::VirtualSolver , public VirtualMatrice<R>{
  int n;
  int nbitermax;
  double eps;
  mutable R  epsr;
  int dKrilov;
  KN<R> D1;
  public:
  SolveGMRESDiag(const MatriceMorse<R> &A,int nbk=50,int itmax=0,double epsilon=1e-6) : 
    n(A.n),nbitermax(itmax?itmax: Max(100,n)),D1(n),eps(epsilon),epsr(0),
    dKrilov(nbk) { 
    R aii=0;
    for (int i=0;i<n;i++)
      D1[i] = (fabs(aii=A(i,i)) < 1e-10 ? 1. : 1./aii);}
   void Solver(const MatriceMorse<R> &a,KN_<R> &x,const KN_<R> &b) const  {
      epsr = (eps < 0) ? (epsr >0 ? -epsr : -eps ) : eps ;
    // cout << " epsr = " << epsr << endl;
   //  ConjuguedGradient<R,MatriceMorse<R>,SolveGCDiag<R> >(a,*this,b,x,nbitermax,epsr);
         KNM<R> H(dKrilov+1,dKrilov+1);
      int k=dKrilov,nn=nbitermax;
      int res=GMRES(a,(KN<R> &)x,(const KN<R> &)b,*this,H,k,nn,epsr);

   }
typename VirtualMatrice<R>::plusAx operator*(const KN_<R> &  x) const {return plusAx(this,x);} 


 void addMatMul(const KN_<R> & x, KN_<R> & Ax) const 
  { 
     assert(x.N()==Ax.N());
   for (int i=0;i<n;i++) 
     Ax[i]+= D1[i]*x[i];}
     

};   
}

