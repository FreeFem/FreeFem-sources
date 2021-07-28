#ifndef MatriceCreuse_h_
#define MatriceCreuse_h_

#define VERSION_MATRICE_CREUSE 4

// add Sep 2007 for generic Space solver
#define DCL_ARG_SPARSE_SOLVER(T,A)  Stack stack, HashMatrix<int,T> *A, Data_Sparse_Solver & ds
#define ARG_SPARSE_SOLVER(A) stack,A, ds

#ifdef HAVE_LIBUMFPACK_XXXXXXXXXXXXX
extern "C" {
#ifdef HAVE_UMFPACK_H
#include <umfpack.h>
#else
#ifdef HAVE_UMFPACK_UMFPACK_H
#include <umfpack/umfpack.h>
#else
#ifdef HAVE_BIG_UMFPACK_UMFPACK_H
#include <UMFPACK/umfpack.h>
#else
#ifdef HAVE_UFSPARSE_UMFPACK_H
#include <ufsparse/umfpack.h>
#else
#ifdef HAVE_SUITESPARSE_UMFPACK_H
#include <suitesparse/umfpack.h>
#else

  // Defaults to a local version of the UMFPACK headers
#include "../../3rdparty/include/umfpack.h"

#endif // HAVE_SUITESPARSE_UMFPACK_H
#endif // HAVE_UFSPARSE_UMFPACK_H
#endif // HAVE_BIG_UMFPACK_UMFPACK_H
#endif // HAVE_UMFPACK_UMFPACK_H
#endif // HAVE_UMFPACK_H
}
#endif

#include "RNM.hpp"
#include "fem.hpp"
#include "FESpace.hpp"
#include "DOperator.hpp"
#include "QuadratureFormular.hpp"
#include "HashMatrix.hpp"

extern double ff_tgv;
using  Fem2D::Mesh;
using  Fem2D::FESpace;
using Fem2D::FElement;
using Fem2D::baseFElement;
using Fem2D::FMortar;
using Fem2D::TypeOfMortar;
using Fem2D::QuadratureFormular;
using Fem2D::QuadratureFormular1d;
using Fem2D::QuadratureFormular_T_5;
using Fem2D::QF_GaussLegendre3;
const double  EPSILON=1e-20;
using Fem2D::onWhatIsEdge;

//#define APROGRAMMER(a) {cerr << "A PROGRAMMER "  #a << endl; exit (1) ;}
#define ERREUR(a,b) {cerr << "ERROR " #a<< b <<endl; throw(ErrorExec("FATAL ERROR in " __FILE__  "\n" #a  " line: ",__LINE__)) ;}

template<class TypeIndex=int,class TypeScalar=double> class VirtualMatrix;
template<class TypeIndex,class TypeScalaire> class HashMatrix ;
template<class K>  using MatriceCreuse=VirtualMatrix<int, K>     ;
template<class K>  using MatriceMorse=HashMatrix<int,K>;

//template <class R> class MatriceCreuse;
template <class R> class MatriceElementaire;
template <class R,class FESpace> class MatriceElementaireSymetrique;
template <class R,class FESpace> class MatriceElementairePleine;
//template <class R> class MatriceMorseOld;
template <class R> class MatriceProdTensoriel;

//template <class R> R Square(R x){ return x*x;}

template <class T> T* docpyornot(bool nocpy,T* p,int n)
{
  T * r=p;
   if( !nocpy) { // do copy
      r= new T[n]; ffassert(r);
      if(p)
       for(int i=0;i<n;i++)
        r[i]=p[i];
      }
    else if( r==0) // always do allocation  July 2009 for mpi matrice
      { r= new T[n]; ffassert(r);}
   return r;
 }
 template <class T,class TT> T* docpy(TT* p,int n)
{
   T * r=0;
   if(p && n) { // do copy
      r= new T[n]; ffassert(r);
      for(int i=0;i<n;i++)
        r[i]=T(p[i]); // pour cadna ???? FH
      }
   return r;
 }




template <typename Z,typename R> class  VirtualMatrix;
template <typename Z,typename R>  class HashMatrix;



template<class R> class StopGC { public: virtual  bool Stop(int iter, R *, R * ){cout << " Stop !!!!!\n"; return false;} };
template<class R,class M,class P,class S >// S=StopGC<Real>
int ConjuguedGradient(const M & A,const P & C,const KN_<R> &b,KN_<R> &x,const int nbitermax, double &eps,long kprint=1000000000,S *Stop=0)
{

//  ConjuguedGradient lineare   A*x est appele avec des conditions au limites
//  non-homogene puis homogene  pour calculer le gradient
   if (verbosity>50)
     kprint=2;
   if (verbosity>99)  cout << A << endl;
//   throwassert(&x && &b && &A && &C);
   typedef KN<R> Rn;
   int n=b.N();
   throwassert(n==x.N());
   Rn g(n), h(n), Ah(n), & Cg(Ah);  // on utilise Ah pour stocke Cg
   g = A*x;
   double xx= std::real((x,conj(x)));
   double epsold=eps;
   g -= b;// g = Ax-b
   Cg = C*g; // gradient preconditionne
   h =-Cg;
   double g2 = std::real((Cg,conj(g)));
   if (g2 < 1e-30)
    { if(verbosity>1 || (kprint<100000))
       cout << "GC  g^2 =" << g2 << " < 1.e-30  Nothing to do " << endl;
     return 2;  }
   double reps2 =eps >0 ?  eps*eps*g2 : -eps; // epsilon relatif
   eps = reps2;
   for (int iter=0;iter<=nbitermax;iter++)
     {
       Ah = A*h;
       double hAh =std::real((h,conj(Ah)));
      // if (Abs(hAh)<1e-30) ExecError("CG2: Matrix non defined, sorry ");
       R ro =  - std::real((g,conj(h)))/ hAh; // ro optimal (produit scalaire usuel)
       x += ro *h;
       g += ro *Ah; // plus besoin de Ah, on utilise avec Cg optimisation
       Cg = C*g;
       double g2p=g2;
       g2 = std::real((Cg,conj(g)));
       bool stop = Stop && Stop->Stop(iter,x,g);

       if ( !(iter%kprint) && iter && (verbosity>3) )
         cout << "CG:" <<iter <<  "  ro = " << ro << " ||g||^2 = " << g2 << " " << stop << endl;

       if (g2 < reps2 || stop) {
         if ( !iter && !xx && g2 && epsold >0 ) {
             // change fo eps converge to fast due to the
             // penalization of boundary condition.
             eps =  epsold*epsold*g2;
             if (verbosity>3 || (kprint<3))
             cout << "CG converge to fast (pb of BC)  restart: " << iter <<  "  ro = "
                  << ro << " ||g||^2 = " << g2 << " <= " << reps2 << "  new eps2 =" <<  eps <<endl;
              reps2=eps;
           }
         else
          {
           if (verbosity>1 || (kprint<100000) )
            cout << "CG converge: " << iter <<  "  ro = " << ro << " ||g||^2 = " << g2 << endl;
           return 1;// ok
          }
          }
       double gamma = g2/g2p;
       h *= gamma;
       h -= Cg;  //  h = -Cg * gamma* h
     }
   if(verbosity)
   cout << " GC: method doesn't converge in " << nbitermax
        <<  " iteration , xx= "  << xx<< endl;
   return 0;
}

template<class R,class M,class P,class S>// S=StopGC<Real>
int ConjuguedGradient2(const M & A,const P & C,KN_<R> &x,const KN_<R> &b,const int nbitermax, double &eps,long kprint=1000000000,S *Stop=0)
{
//  ConjuguedGradient2 affine A*x = 0 est toujours appele avec les condition aux limites
//  -------------

  // throwassert(&x  && &A && &C);
   typedef KN<R> Rn;
   int n=x.N();
  // if (verbosity>99) kprint=1;
   R ro=1;
   Rn g(n),h(n),Ah(n), & Cg(Ah);  // on utilise Ah pour stocke Cg
   g = A*x;
   g -= b;
   Cg = C*g; // gradient preconditionne
   h =-Cg;
   R g2 = (Cg,g);
   if (g2 < 1e-30)
    { if(verbosity>1 || kprint< 1000000)
       cout << "GC  g^2 =" << g2 << " < 1.e-30  Nothing to do " << endl;
     return 2;  }
   if (verbosity>5 || (kprint<2))
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
       R hAh =(h,Ah);
         if (std::norm(hAh)<1e-60) ExecError("CG2: Matrix is not defined (/0), sorry ");
       ro =  - (g,h)*rop/hAh ; // ro optimal (produit scalaire usuel)
       if ( ro != ro) ExecError("CG2: Bug : ro is NaN ???  ");
       x += (ro-rop) *h;
       g += (ro/rop) *Ah; // plus besoin de Ah, on utilise avec Cg optimisation
       Cg = C*g;
       R g2p=g2;
       g2 = (Cg,g);
         bool stop = Stop && Stop->Stop(iter,x,g);
       if ( ( (iter%kprint) == kprint-1)  /*&&  verbosity >1*/ )
         cout << "CG:" <<iter <<  "  ro = " << ro << " ||g||^2 = " << g2 << " " << stop <<  endl;
       if (stop || g2 < reps2 ) {
         if (kprint <= nbitermax )
            cout << "CG converges " << iter <<  "  ro = " << ro << " ||g||^2 = " << g2
                 << " stop=" << stop << " /Stop= " << Stop <<  endl;
          return 1;// ok
          }
       R gamma = g2/g2p;
       h *= gamma;
       h -= Cg;  //  h = -Cg * gamma* h
     }
     //if (nbitermax <= nbitermax )
      cout << "CG doesn't converge: " << nbitermax <<   " ||g||^2 = " << g2 << " reps2= " << reps2 << endl;
   return 0;
}

template <class R>
class MatriceIdentite:public  RNM_VirtualMatrix<R> { public:
 typedef typename RNM_VirtualMatrix<R>::plusAx plusAx;
    MatriceIdentite(int n) :RNM_VirtualMatrix<R>(n) {};
 void addMatMul(const  KN_<R>  & x, KN_<R> & Ax) const {
     ffassert(x.N()==Ax.N());
   Ax+=x; }
    void addMatTransMul(const  KN_<R>  & x, KN_<R> & Ax) const {
    ffassert(x.N()==Ax.N());
    Ax+=x; }
     void Solve( KN_<R> & y ,const KN_<R> & x) const { y=x; }
      bool WithSolver() const {return true;}
   typename  RNM_VirtualMatrix<R>::plusAx operator*(const KN<R> &  x) const {return typename RNM_VirtualMatrix<R>::plusAx(this,x);}
  bool ChecknbLine(int n) const { return true;}
  bool ChecknbColumn(int m) const { return true;}

};

template<class R>
class SolveGCDiag :   public MatriceMorse<R>::VirtualSolver , public RNM_VirtualMatrix<R>{
  int n;
  int nbitermax;
  double eps;
  mutable double  epsr;
  KN<R> D1;
  public:
  typedef typename RNM_VirtualMatrix<R>::plusAx plusAx;
  SolveGCDiag(const MatriceMorse<R> &A,int itmax,double epsilon=1e-6) :
    RNM_VirtualMatrix<R>(A.n),
    n(A.n),nbitermax(itmax?itmax: Max(100,n)),eps(epsilon),epsr(0),D1(n)
  { //throwassert(A.sym());
    for (int i=0;i<n;i++)
      D1[i] = 1./A(i,i);}
   void Solver(const MatriceMorse<R> &a,KN_<R> &x,const KN_<R> &b) const  {
     epsr = (eps < 0) ? (epsr >0 ? -epsr : -eps ) : eps ;
    // cout << " epsr = " << epsr << endl;
     ConjuguedGradient<R,MatriceMorse<R>,SolveGCDiag<R>,StopGC<R> >(a,*this,b,x,nbitermax,epsr  );
   }
plusAx operator*(const KN_<R> &  x) const {return plusAx(this,x);}


 void addMatMul(const KN_<R> & x, KN_<R> & Ax) const
  {  ffassert(x.N()==Ax.N());
     for (int i=0;i<n;i++)
     Ax[i]+= D1[i]*x[i];}

   bool ChecknbLine(int nn) const { return n==nn;}
   bool ChecknbColumn(int mm) const { return n==mm;}

};




// alloc solver and solevr def
//template<class R,int sympos>   typename DefSparseSolverNew<R,sympos>::SparseMatSolver DefSparseSolverNew<R,sympos>::solver=0;
//template<class R,int sympos>   typename DefSparseSolverNew<R,sympos>::SparseMatSolver DefSparseSolverNew<R,sympos>::solverdef=0;

// End Sep 2007 for generic Space solver



inline void C2RR(int n,Complex *c,double *cr,double *ci)
{
 for (int i=0;i<n;i++)
  {
   cr[i]=real(c[i]);
   ci[i]=imag(c[i]);
  }
}

inline void RR2C(int n,double *cr,double *ci,Complex *c)
{
 for (int i=0;i<n;i++)
  {
    c[i]=Complex(cr[i],ci[i]);
  }
}



#endif
