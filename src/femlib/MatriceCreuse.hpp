#ifndef MatriceCreuse_h_
#define MatriceCreuse_h_

template<class T>
 T Square(const T & r){return r*r;}

 
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
#include "../../download/include/umfpack.h"

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
template <class R> class MatriceMorseOld;
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




template <class R>
class MatriceCreuseOld : public RefCounter,public RNM_VirtualMatrix<R> {
public:
    MatriceCreuseOld(int NbOfDF,int mm,int ddummy)
    : RNM_VirtualMatrix<R>(NbOfDF,mm),n(NbOfDF),m(mm),dummy(ddummy){}
    MatriceCreuseOld(int NbOfDF)
    : RNM_VirtualMatrix<R>(NbOfDF),n(NbOfDF),m(NbOfDF),dummy(1){}
    int n,m,dummy;
    virtual int size() const =0;
    
    virtual MatriceCreuseOld & operator +=(MatriceElementaire<R> & )=0;
    virtual void operator=(const R & v) =0; // Mise a zero
    KN_<R> & MatMul(KN_<R> &ax,const KN_<R> &x) const {
        ax= R();
        addMatMul(x,ax);
        return ax;}
    virtual ostream& dump (ostream&)  const =0;
    virtual void Solve(KN_<R> & x,const KN_<R> & b) const =0;
    virtual void SolveT(KN_<R> & x,const KN_<R> & b) const {  ERREUR(SolveT, "No Solver of A^t in this matrix type ???"); }
    virtual ~MatriceCreuseOld(){}
    virtual R & diag(int i)=0;
    virtual void SetBC(int i,double tgv)=0;
    virtual void SetBC(char *wbc,double tgv) { for (int i=0; i<n; ++i)  if(wbc[i]) SetBC(i,tgv);}
    virtual R & operator()(int i,int j)=0;
    virtual R * pij(int i,int j) const =0; // Add FH
    virtual  void  resize(int n,int m)  {AFAIRE("MatriceCreuseOld::resize");}  // a faire dans les classe derive ... // add march 2009  FH
    virtual MatriceMorseOld<R> *toMatriceMorse(bool transpose=false,bool copy=false) const {return 0;} // not
    virtual bool addMatTo(R coef,std::map< pair<int,int>, R> &mij,bool trans=false,int ii00=0,int jj00=0,bool cnj=false,double threshold=0.,const bool keepSym=false)=0;
    // Add FH april 2005
    virtual R pscal(const KN_<R> & x,const KN_<R> & y) =0 ; // produit scalaire
    virtual double psor(KN_<R> & x,const  KN_<R> & gmin,const  KN_<R> & gmax , double omega) =0;
    virtual void setdiag(const KN_<R> & x)=0 ;
    virtual void getdiag( KN_<R> & x) const =0 ;
    // end add
    virtual int NbCoef() const {return 0;};
    virtual void setcoef(const KN_<R> & x)=0 ;
    virtual void getcoef( KN_<R> & x) const =0 ;
    // Add FH oct 2005
    bool ChecknbLine(int nn) const { return n==nn;}
    bool ChecknbColumn(int mm) const { return m==mm;}
    virtual R trace() const {ffassert(n==m);  R t=R(), *p;  for(int i=0; i<n; ++i)  { p=pij(i,i);  if(p) t+= *p;} return t; }
    // end ADD
    virtual bool sym() const {return false;}
    
};


template <typename Z,typename R> class  VirtualMatrix;
template <typename Z,typename R>  class HashMatrix;
#include "SparseLinearSolver.hpp"

#ifdef REMOVE_CODE_OBSO
template <class R>  class MatriceProfile;
template <class R> 
<<<<<<< ours
inline ostream& operator <<(ostream& f,const MatriceCreuseOld<R> & m)
=======
class MatriceCreuse : public RefCounter,public VirtualMatrice<R> {
public:
  MatriceCreuse(int NbOfDF,int mm,int ddummy)
         : VirtualMatrice<R>(NbOfDF,mm),n(NbOfDF),m(mm),dummy(ddummy){}
  MatriceCreuse(int NbOfDF)
         : VirtualMatrice<R>(NbOfDF),n(NbOfDF),m(NbOfDF),dummy(1){}
  int n,m,dummy;
  virtual int size() const =0;

  virtual MatriceCreuse & operator +=(MatriceElementaire<R> & )=0;
  virtual void operator=(const R & v) =0; // Mise a zero 
  KN_<R> & MatMul(KN_<R> &ax,const KN_<R> &x) const { 
    ax= R();
    addMatMul(x,ax);
    return ax;}
  virtual ostream& dump (ostream&)  const =0;
  virtual void Solve(KN_<R> & x,const KN_<R> & b) const =0;
  virtual void SolveT(KN_<R> & x,const KN_<R> & b) const {  ERREUR(SolveT, "No Solver of A^t in this matrix type ???"); }
  virtual ~MatriceCreuse(){}
  virtual R & diag(int i)=0;
  virtual void SetBC(int i,double tgv)=0;
  virtual void SetBC(char *wbc,double tgv) { for (int i=0; i<n; ++i)  if(wbc[i]) SetBC(i,tgv);}
  virtual R & operator()(int i,int j)=0;
  virtual R * pij(int i,int j) const =0; // Add FH 
  virtual  void  resize(int n,int m)  {AFAIRE("MatriceCreuse::resize");}  // a faire dans les classe derive ... // add march 2009  FH 
  virtual MatriceMorse<R> *toMatriceMorse(bool transpose=false,bool nocopy=false) const {return 0;} // not 
  virtual bool addMatTo(R coef,std::map< pair<int,int>, R> &mij,bool trans=false,int ii00=0,int jj00=0,bool cnj=false,double threshold=0.,const bool keepSym=false)=0;
  // Add FH april 2005
  virtual R pscal(const KN_<R> & x,const KN_<R> & y) =0 ; // produit scalaire  
  virtual double psor(KN_<R> & x,const  KN_<R> & gmin,const  KN_<R> & gmax , double omega) =0;
  virtual void setdiag(const KN_<R> & x)=0 ;
  virtual void getdiag( KN_<R> & x) const =0 ;
  // end add
  virtual int NbCoef() const {return 0;};
  virtual void setcoef(const KN_<R> & x)=0 ;
  virtual void getcoef( KN_<R> & x) const =0 ;
  // Add FH oct 2005
   bool ChecknbLine(int nn) const { return n==nn;}  
   bool ChecknbColumn(int mm) const { return m==mm;}
   virtual R trace() const {ffassert(n==m);  R t=R(), *p;  for(int i=0; i<n; ++i)  { p=pij(i,i);  if(p) t+= *p;} return t; }
  // end ADD
  virtual bool sym() const {return false;}

};

template <class R> 
inline ostream& operator <<(ostream& f,const MatriceCreuse<R> & m) 
>>>>>>> theirs
    {return m.dump(f);}

template <class R> 
KN_<R> & operator/=(KN_<R> & x ,const MatriceProfile<R> & a) ;

/*
enum FactorizationType {
    FactorizationNO=0,
    FactorizationCholeski=1,
    FactorizationCrout=2,
    FactorizationLU=3};
*/
template <class R> 
class MatriceProfile:public MatriceCreuseOld<R> {
public:
  mutable R *L;  // lower 
  mutable R *U;  // upper
  mutable R *D;  // diagonal 
  int *pL; // profile L 
  int *pU; // profile U 
  mutable  FactorizationType typefac;
  FactorizationType typesolver; 
  ostream& dump (ostream&) const ;
  MatriceProfile(const int  n,const R *a);
  
  template<class FESpace>
  MatriceProfile(const FESpace &,bool VF=false);
  MatriceProfile(int NbOfDF,R* d,
		 R* u, int * pu,
		 R* l, int * pl,
		 FactorizationType tf=FactorizationNO) 
    : MatriceCreuseOld<R>(NbOfDF),L(l),U(u),D(d),pL(pl),pU(pu),
      typefac(tf),typesolver(FactorizationNO){}

  const MatriceProfile t() const   
     {return MatriceProfile(this->n,D,L,pL,U,pU,typefac);}
  const MatriceProfile lt() const  
     {return MatriceProfile(this->n,0,L,pL,0,0);}
  const MatriceProfile l() const  
     {return MatriceProfile(this->n,0,0,0,L,pL);}
  const MatriceProfile d() const   
     {return MatriceProfile(this->n,D,0,0,0,0);}
  const MatriceProfile ld() const  
     {return MatriceProfile(this->n,D,0,0,L,pL);}
  const MatriceProfile ldt() const 
     {return MatriceProfile(this->n,D,L,pL,0,0);}
  const MatriceProfile dut() const
    {return MatriceProfile(this->n,D,0,0,U,pU);}
  const MatriceProfile du() const
     {return MatriceProfile(this->n,D,U,pU,0,0);}
  const MatriceProfile u() const   
     {return MatriceProfile(this->n,0,U,pU,0,0);}
  const MatriceProfile ut() const
    {return MatriceProfile(this->n,0,0,0,U,pU);}
   R trace() const {ffassert(this->n==this->m);  R t=R();  for(int i=0; i<this->n; ++i) t+= D[i]; return t; }
  void Solve(KN_<R> &x,const KN_<R> &b) const {
    /*if (typefac==0)  code faux // FH   nov 2006
      switch(typefac) 
	{
	FactorizationCholeski: cholesky() ; break;
	FactorizationCrout:   crout(); break;
	FactorizationLU:      LU(); break; 
	}*/ 
    if (&x != &b) x=b;x/=*this;} 
    void SolveT(KN_<R> &x,const KN_<R> &b) const {
        if( L==U) Solve(x,b);// symetrique
        else {
           if( typefac == FactorizationNO && (U && L))
             {cerr << "APROGRAMMER (solveT /MatriceProfile)";throw(ErrorExec("exit",2));}
           if (&x != &b) x=b;
            if( U && L )
            {//  no trans:  x  /= l(); x  /= du();=> trans => x  /= dut(); x  /= lt();
               // cout << "\n x  /= l(); x  /= du();=> trans => x  /= dut(); x  /= lt();" << endl;
              x  /= dut();
              x  /= lt();
            }
            else if (U && D == 0)
                 x  /= ut();
            else if (L && D )
                 x  /= ldt();
            else {ERREUR(SolveT, "SolveT profile unknow case ");}

        }
    }

  int size() const ;
  void  resize(int n,int m)  { AFAIRE("MatriceProfile::resize");}  // a faire ...  add march 2009  FH 
  ~MatriceProfile();
  //  KN_<R>         operator* (const KN_<R> & ) const ;
  void addMatMul(const KN_<R> &x,KN_<R> &ax) const;
  void addMatTransMul(const KN_<R> &x,KN_<R> &ax) const ;
  //  { this->t().addMatMul(x,ax);}
  MatriceCreuseOld<R> & operator +=(MatriceElementaire<R> &);
  void operator=(const R & v); // Mise a zero 
  void cholesky(double = EPSILON/8.) const ; //
  void crout(double = EPSILON/8.) const ; //
  void LU(double = EPSILON/8.) const ; //
  R & diag(int i) { return D[i];}
  void SetBC (int i,double tgv) {
      if( tgv>=0) D[i]=tgv;
      else  { ffassert(tgv<0); }  // to hard ..
  }
    
  R & operator()(int i,int j) { if(i!=j) ffassert(0); return D[i];} // a faire 
  R * pij(int i,int j) const { if(i!=j) ffassert(0); return &D[i];} // a faire  Modif FH 31102005
<<<<<<< ours
  MatriceMorseOld<R> *toMatriceMorse(bool transpose=false,bool copy=false) const ;
=======
  MatriceMorse<R> *toMatriceMorse(bool transpose=false,bool nocopy=false) const ;
>>>>>>> theirs
  
  template<class F> void map(const  F & f)
  {
    for(int i=0;i<this->n;++i)
      D[i]=f(D[i]);
    if (L)
    for(int i=0;i<pL[this->n];++i)
      L[i]=f(L[i]);
    if (L && (L != U) )
    for(int i=0;i<pL[this->m];++i)
      U[i]=f(U[i]);
  }
  
  template<class RR> MatriceProfile(const MatriceProfile<RR> & A)
    : MatriceCreuseOld<R>(A.n,A.m,0)
  {
    
    typefac=A.typefac;
    pL=  docpy<int,int>(A.pL,this->n+1);
    D = docpy<R,RR>(A.D,this->n);
    if ( A.pL == A.pU ) pU=pL;
    else pU=  docpy<int,int>(A.pU,this->m+1);
    
      L= docpy<R,RR>(A.L,pL[this->n]);
      
    if ( A.L == A.U ) U=L;
    else  U= docpy<R,RR>(A.U,pU[this->m]);
    
  
  }

  
  bool addMatTo(R coef,std::map< pair<int,int>, R> &mij,bool trans=false,int ii00=0,int jj00=0,bool cnj=false,double threshold=0.,const bool keepSym=false);

  // Add FH april 2005
  R pscal(const KN_<R> & x,const KN_<R> & y); // produit scalaire  
  double psor(KN_<R> & x,const  KN_<R> & gmin,const  KN_<R> & gmax , double omega);
  void setdiag(const KN_<R> & x) ;
  void getdiag( KN_<R> & x) const ;
    // end add
  // Add FH oct 2005
   int NbCoef() const ;
   void setcoef(const KN_<R> & x);
   void getcoef( KN_<R> & x) const ;
  // end add
  
  /*----------------------------------------------------------------
    D[i] = A[ii]
    L[k] = A[ij]  j < i avec:   pL[i]<= k < pL[i+1] et j = pL[i+1]-k
    U[k] = A[ij]  i < j avec:   pU[j]<= k < pU[j+1] et i = pU[i+1]-k 
    remarque pL = pU generalement 
    si L = U => la matrice est symetrique 
    -------------------------------------------------------------------
  */ 
  private:
   void operator=(const MatriceProfile & A);
};



template <class R> 
class MatriceMorseOld:public MatriceCreuseOld<R> {
//  numebering  is no-symmetric
//  the all line  i :  
//     k=   lg[i] .. lg[i+1]+1
//        j = cl[k]
//        aij=a[k]
// otherwise  symmetric  case
// same but just the  LOWER part is store     (j <= i) 
// and aii exist always in symmetric case
//  -----------------------------------------

public:
  int nbcoef;
  bool symetrique;  
  R * a;
  int * lg; 
 
  int * cl;  
public:

    
 class VirtualSolver :public RefCounter { 
   friend class MatriceMorseOld;
   virtual void Solver(const MatriceMorseOld<R> &a,KN_<R> &x,const KN_<R> &b) const  =0;
     virtual void SolverT(const MatriceMorseOld<R> &a,KN_<R> &x,const KN_<R> &b) const  {ExecError("No Solver of Transpose matrix");};
};

  MatriceMorseOld():MatriceCreuseOld<R>(0),nbcoef(0),symetrique(true),a(0),lg(0),cl(0),solver(0) {};
  MatriceMorseOld(KNM_<R> & A, double tol) ;
  MatriceMorseOld(const int  n,const R *a);
//    :MatriceCreuseOld<R>(n),solver(0) {}
  MatriceMorseOld(istream & f);

  template<class FESpace>   explicit 
  MatriceMorseOld(const FESpace & Uh,bool sym,bool VF=false)
    :MatriceCreuseOld<R>(Uh.NbOfDF),solver(0) {Build(Uh,Uh,sym,VF);}

  template<class FESpace>   explicit 
  MatriceMorseOld(const FESpace & Uh,const FESpace & Vh,bool VF=false)
    :MatriceCreuseOld<R>(Uh.NbOfDF,Vh.NbOfDF,0),solver(0)
  {Build(Uh,Vh,false,VF);}

  template<class FESpace>   explicit 
  MatriceMorseOld(const FESpace & Uh,const FESpace & Vh,
               void (*build)(MatriceMorseOld *,const FESpace & Uh,const FESpace & Vh,void *data),void *data=0
	       )
    :MatriceCreuseOld<R>(Uh.NbOfDF,Vh.NbOfDF,0),solver(0)
  {build(this,Uh,Vh,data);           
  }
 
MatriceMorseOld(int nn,int mm,int nbc,bool sym,R *aa=0,int *ll=0,int *cc=0,bool dd=false, const VirtualSolver * s=0,bool transpose=false )
    :MatriceCreuseOld<R>(nn,mm,dd && !transpose),
     nbcoef(nbc),
     symetrique(sym), // transpose = true => dummy false (new matrix)
     a(docpyornot(this->dummy,aa,nbc)),
     lg(docpyornot(this->dummy,ll,nn+1)),
     cl(docpyornot(this->dummy,cc,nbc)),
     solver(s)
  { if(transpose) dotransposition(); };
  void Solve(KN_<R> &x,const KN_<R> &b) const;
  void SolveT(KN_<R> &x,const KN_<R> &b) const;
  int size() const ;
  void addMatMul(const KN_<R> &x,KN_<R> &ax) const;
  void addMatTransMul(const KN_<R> &x,KN_<R> &ax) const;
  MatriceMorseOld & operator +=(MatriceElementaire<R> &);
  void operator=(const R & v) 
    { for (int i=0;i< nbcoef;i++) a[i]=v;}
  virtual ~MatriceMorseOld(){ if (!this->dummy) { delete [] a; delete [] cl;delete [] lg;}}
  ostream&  dump(ostream & f) const ;
  R * pij(int i,int j) const ;
  R  operator()(int ii,int jj) const {R * p= pij(ii,jj) ;throwassert(p); return *p;}
  R & operator()(int ii,int jj)  {R * p= pij(ii,jj) ;throwassert(p); return *p;}
  R & diag(int ii)  {R * p= pij(ii,ii) ;throwassert(p); return *p;}
    R trace() const {ffassert(this->n==this->m);  R t=R(),*p;  for(int i=0; i<this->n; ++i) {p=pij(i,i) ;if(p) t+= *p; } return t; }
  void SetBC (int i,double tgv) {
	R * p= pij(i,i) ;
	ffassert(p);
	if( tgv>=0) *p=tgv;
	else  {
            ffassert(!symetrique);// 
	    for (int k=lg[i];k<lg[i+1]; ++k) a[k]=0;// put the line to Zero.
	    *p = 1. ; // and the diag coef to 1.
	}  
    }
    
void SetBC (char * wbc,double tgv)
    {
        for(int i=0; i< this->n; ++i)
            if(tgv<0)
            {
                
                if( wbc[i] )
                {
                    for (int k=lg[i];k<lg[i+1]; ++k)
                        if( cl[k]==i)
                            a[k] = 1.;
                        else
                            a[k]=0;// put the line to Zero.
                }
                else if(tgv < -1.999)
                    for (int k=lg[i];k<lg[i+1]; ++k)
                    {
                        int j = cl[k];
                        if( wbc[j] ) a[k]=0;//
                    }
            }
            else  if( wbc[i] ) { // tgv >= 0
                R * p= pij(i,i) ;
                ffassert(p);
                *p=tgv;
            }
    }
  
  void SetSolver(const VirtualSolver & s){solver=&s;}
  template<class T>
  void GetSolver(const T*& s) { if(solver) s = dynamic_cast<const T*>(&*solver); else s = 0; }
  void SetSolverMaster(const VirtualSolver * s){solver.master(s);}
  bool sym() const {return symetrique;}
  // Add FH april 2005
  R pscal(const KN_<R> & x,const KN_<R> & y); // produit scalaire  
  double psor(KN_<R> & x,const  KN_<R> & gmin,const  KN_<R> & gmax , double omega);
  void setdiag(const KN_<R> & x) ;
  void getdiag( KN_<R> & x) const ;
  // end add
  
  // Add FH oct 2005
   int NbCoef() const ;
   void setcoef(const KN_<R> & x);
   void getcoef( KN_<R> & x) const ;
  // end add
void  resize(int n,int m) ; // add march 2009 ...
template<class K>
  MatriceMorseOld(int nn,int mm, std::map< pair<int,int>, K> & m, bool sym);
  
 template<class RB,class RAB>
 void  prod(const MatriceMorseOld<RB> & B, MatriceMorseOld<RAB> & AB);
 
<<<<<<< ours
 MatriceMorseOld<R> *toMatriceMorse(bool transpose=false,bool copy=false) const {
     return new MatriceMorseOld(this->n,this->m,nbcoef,symetrique,a,lg,cl,copy, solver,transpose);}
=======
 MatriceMorse<R> *toMatriceMorse(bool transpose=false,bool nocopy=false) const {
     return new MatriceMorse(this->n,this->m,nbcoef,symetrique,a,lg,cl,nocopy, solver,transpose);}
>>>>>>> theirs
  bool  addMatTo(R coef,std::map< pair<int,int>, R> &mij,bool trans=false,int ii00=0,int jj00=0,bool cnj=false,double threshold=0.,const bool keepSym=false);
  
  template<typename RR,typename K> static  RR CastTo(K  b){return b;}
                                                  
  template<class K>
    MatriceMorseOld(const MatriceMorseOld<K> & , R (*f)(K) );
  template<class K>
    MatriceMorseOld(const MatriceMorseOld<K> & );

  private:
  void dotransposition ()  ;  // do the transposition 
  CountPointer<const VirtualSolver> solver;
  
  void operator=(const MatriceMorseOld & );

  template<class FESpace>
  void  Build(const FESpace & Uh,const FESpace & Vh,bool sym,bool VF=false);
  
};
#endif

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
   double xx= RNM::real((x,conj(x)));
   double epsold=eps;
   g -= b;// g = Ax-b
   Cg = C*g; // gradient preconditionne 
   h =-Cg; 
   double g2 = RNM::real((Cg,conj(g)));
   if (g2 < 1e-30) 
    { if(verbosity>1 || (kprint<100000))
       cout << "GC  g^2 =" << g2 << " < 1.e-30  Nothing to do " << endl;
     return 2;  }
   double reps2 =eps >0 ?  eps*eps*g2 : -eps; // epsilon relatif 
   eps = reps2;
   for (int iter=0;iter<=nbitermax;iter++)
     {      
       Ah = A*h;
       double hAh =RNM::real((h,conj(Ah)));
      // if (Abs(hAh)<1e-30) ExecError("CG2: Matrix non defined, sorry ");
       R ro =  - RNM::real((g,conj(h)))/ hAh; // ro optimal (produit scalaire usuel)
       x += ro *h;
       g += ro *Ah; // plus besoin de Ah, on utilise avec Cg optimisation
       Cg = C*g;
       double g2p=g2; 
       g2 = RNM::real((Cg,conj(g)));
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
                  << ro << " ||g||^2 = " << g2 << " <= " << reps2  <<endl;
              reps2=eps;
           } 
         else 
          { 
           if (verbosity>1 || (kprint<100000) )
            cout << "CG converge: " << iter <<  "  ro = " << ro << " ||g||^2 = " << g2 << " eps2 " << reps2  <<endl;
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
         if (RNM::norm2(hAh)<1e-60) ExecError("CG2: Matrix is not defined (/0), sorry ");
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
      cout << "CG does'nt converge: " << nbitermax <<   " ||g||^2 = " << g2 << " reps2= " << reps2 << endl; 
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
    bool x0;
    double *veps;
  public:
<<<<<<< ours
  typedef typename RNM_VirtualMatrix<R>::plusAx plusAx;
  SolveGCDiag(const MatriceMorse<R> &A,int itmax,double epsilon=1e-6) : 
    RNM_VirtualMatrix<R>(A.n),
    n(A.n),nbitermax(itmax?itmax: Max(100,n)),eps(epsilon),epsr(0),D1(n)
=======
  typedef typename VirtualMatrice<R>::plusAx plusAx;
  SolveGCDiag(const MatriceMorse<R> &A,int itmax,double epsilon=1e-6,bool x00=true,double *vveps=0) :
    VirtualMatrice<R>(A.n),
    n(A.n),nbitermax(itmax?itmax: Max(100,n)),eps(epsilon),epsr(0),D1(n),x0(x00),veps(vveps)
>>>>>>> theirs
  { //throwassert(A.sym());
    for (int i=0;i<n;i++)
      D1[i] = 1./A(i,i);}
   void Solver(const MatriceMorse<R> &a,KN_<R> &x,const KN_<R> &b) const  {
       epsr = veps ? *veps : eps;
   //    cout << " GC0 " << epsr<< " " << *veps << endl;
   //  epsr = (eps < 0) ? (epsr >0 ? -epsr : -eps ) : eps ;
    // cout << " epsr = " << epsr << endl;
       if(x0) x=R();
       double xl2 = x.l2() ;
     ConjuguedGradient<R,MatriceMorse<R>,SolveGCDiag<R>,StopGC<R> >(a,*this,b,x,nbitermax,epsr  );
    //  if(veps)  cout <<"  CG:" << epsr << " " << *veps << " " << xl2 <<endl;
       if(veps) *veps=epsr; 
   }
plusAx operator*(const KN_<R> &  x) const {return plusAx(this,x);} 


 void addMatMul(const KN_<R> & x, KN_<R> & Ax) const 
  {  ffassert(x.N()==Ax.N());
     for (int i=0;i<n;i++) 
     Ax[i]+= D1[i]*x[i];}
     
   bool ChecknbLine(int nn) const { return n==nn;}  
   bool ChecknbColumn(int mm) const { return n==mm;} 

};

#ifdef REMOVE_CODE
#define VDATASPARSESOLVER  1
struct TypeSolveMat {
    enum TSolveMat { NONESQUARE=0, LU=1, CROUT=2, CHOLESKY=3, GC = 4 , GMRES = 5, SparseSolver=6, SparseSolverSym=7 };
    TSolveMat t;
    bool sym;
    bool profile;
    TypeSolveMat(TSolveMat tt=LU) :t(tt),
    sym(t == CROUT || t ==CHOLESKY  ||  t==GC || t==SparseSolverSym ),
    profile(t == CROUT || t ==CHOLESKY  || t ==LU ) {}
    bool operator==(const TypeSolveMat & a) const { return t == a.t;}                               
    bool operator!=(const TypeSolveMat & a) const { return t != a.t;}
    static TSolveMat defaultvalue;
};

// add FH , JM  avril 2009 
template<class K,class V> class MyMap;
class String; 
typedef void *    pcommworld; // to get the pointeur to the comm word ... in mpi
//  to build 

int Data_Sparse_Solver_version() ; //{ return VDATASPARSESOLVER;}


struct Data_Sparse_Solver {
  bool initmat;
  TypeSolveMat* typemat;
  double epsilon;
  const void * precon;
  int NbSpace;
  int strategy;
  double tgv;
  bool factorize;
  double tol_pivot;
  double tol_pivot_sym;
  int itmax ;
  string data_filename;  
  KN<long> lparams;  //  copy arry more secure ...
  KN<double> dparams;   
  
  MyMap<String,String> * smap;   
  
  KN<long> perm_r; 
  KN<long> perm_c;     
  KN<double> scale_r; 
  KN<double> scale_c; 
  string sparams;  
  pcommworld commworld;  // pointeur sur le commworld
    int master; //  master rank in comm add FH 02/2013 for MUMPS ... => VDATASPARSESOLVER exist
    // array for return information for mumps ...
    KN<double> * rinfo;
    KN<long> * info;
    
    KNM<double>* kerneln;
    KNM<double> * kernelt;
    long *kerneldim;
    bool x0; //  init by 0 the inital data the solution
    double * veps; //    to get and set value of eps
 /*   
  int *param_int;
  double *param_double;
  string *param_char;
  int *perm_r;
  int *perm_c;
  string *file_param_int;
  string *file_param_double;
  string *file_param_char;
  string *file_param_perm_r;
  string *file_param_perm_c;
  */
    
    Data_Sparse_Solver()
	:
    initmat(1),
    typemat(0),
    strategy(0),
    tgv(ff_tgv),
    factorize(0),
    epsilon(1e-6),
    precon(0),
    tol_pivot(-1),
    tol_pivot_sym(-1),
    NbSpace(50),
    itmax(0),
 //   lparams(0,0),
 //   dparams(0,0),
    smap(0) ,
//    perm_r(0,0),
//    perm_c(0,0),
//    scale_r(0,0),
 //   scale_c(0,0)
    
    /*
    param_int(0),
    param_double(0),
    param_char(0),
    perm_r(0),
    perm_c(0),
    file_param_int(0),
    file_param_double(0),
    file_param_perm_r(0),
    file_param_perm_c(0),
     */
    //sparams, 
    commworld(0),
    master(0),
    rinfo(0),
    info(0),
    kerneln(0), kernelt(0), kerneldim(0),x0(true),veps(0) 
    {}
    
private:
    Data_Sparse_Solver(const Data_Sparse_Solver& ); // pas de copie 
};
#endif

// add Sep 2007 for generic Space solver
#define DCL_ARG_SPARSE_SOLVER(T,A)  Stack stack, MatriceMorse<T> *A, Data_Sparse_Solver & ds
#define ARG_SPARSE_SOLVER(A) stack,A, ds					

template<class K>  using VirtualSolverN=typename VirtualMatrix<int,K>::VSolver;
typedef VirtualSolverN<double>  * (*SparseRMatSolve)(DCL_ARG_SPARSE_SOLVER(double,A) );
typedef VirtualSolverN<Complex> * (*SparseCMatSolve)(DCL_ARG_SPARSE_SOLVER(Complex,A) );


template<class R,int sympos> struct DefSparseSolverNew {
  typedef VirtualSolverN<R>  *
  (*SparseMatSolver)(DCL_ARG_SPARSE_SOLVER(R,A) );
  static SparseMatSolver solver;
  static SparseMatSolver solverdef;
  static   VirtualSolverN<R> * Build( DCL_ARG_SPARSE_SOLVER(R,A) )
  {
     VirtualSolverN<R> *ret=0;
    if(solver)
      ret =(solver)(ARG_SPARSE_SOLVER(A));
    return ret;	
  }
};

// alloc solver and solevr def 
template<class R,int sympos>   typename DefSparseSolverNew<R,sympos>::SparseMatSolver DefSparseSolverNew<R,sympos>::solver=0;
template<class R,int sympos>   typename DefSparseSolverNew<R,sympos>::SparseMatSolver DefSparseSolverNew<R,sympos>::solverdef=0;

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
