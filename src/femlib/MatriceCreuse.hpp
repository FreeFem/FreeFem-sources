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
#define ERREUR(a,b) {cerr << "ERREUR " #a<< b <<endl; throw(ErrorExec("FATAL ERREUR dans " __FILE__  "\n" #a  " line: ",__LINE__)) ;}

template <class R> class MatriceCreuse; 
template <class R> class MatriceElementaire; 
template <class R,class FESpace> class MatriceElementaireSymetrique;
template <class R,class FESpace> class MatriceElementairePleine;
template <class R> class MatriceMorse;
template <class R> class MatriceProdTensoriel;

//template <class R> R Square(R x){ return x*x;}

template <class T> T* docpyornot(bool nocpy,T* p,int n)
{ 
  T * r=p;
   if( !nocpy) { // do copy 
      r= new T[n]; ffassert(r);
      for(int i=0;i<n;i++) 
        r[i]=p[i];
      }
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
class MatriceElementaire {

public:
  enum TypeOfMatriceElementaire {Full=1,Symmetric=2};


       
  int lga;  // size of array a    
  R* a;          // array  coef --
  int *ni,*nj;   //  list of df   
  // to build matrice on face or edge -----

  
  int n,m;       // n,m number of df
  const TypeOfMatriceElementaire mtype;
  KN<double> data; // to store value of basic function 

  const bool onFace ; //  true if do int on face or edge with jump (VF or GD : Galerkin Discontinus)
  // in with case  add ... 
  const int lnk; // size of the 4 next array
  int *nik,*nikk;  //  number of df in element   k,kk for VF and GD methode 
  int *njk,*njkk;  //  number of df in element   k,kk  for VF and GD methode


  MatriceElementaire(int datasize,int llga
                     ,int *nnj,int * nni,TypeOfMatriceElementaire t=Full)
    
    :   lga(llga),a(new R[lga]),
        ni(nni),nj(nnj),n(0),m(0),mtype(t),data(datasize),
        onFace(false),lnk(0),nik(0),nikk(0),njk(0),njkk(0) 
        {}
       

 //  for discontinous Galerkine method
  MatriceElementaire(int datasize,int llga,int *nni,
                     int lk,
                     TypeOfMatriceElementaire t=Symmetric
                     ) 
    :
    lga(llga),a(new R[lga]),
    ni(nni),nj(nni),n(0),m(0),mtype(t),data(datasize*(lk?2:1)) ,
       onFace(lk!=0),
       lnk(lk),
       nik(lk? new int[lk*2]:0),
       nikk(nik+lk),
       njk(nik),
       njkk(nik+lk)
       { ffassert(lk>=0);}
 
  virtual ~MatriceElementaire() {
    if(ni != nj) 
      delete [] nj;
         
    delete [] ni;
    delete [] a;
    if ( nik) delete[] nik;
     }
      
  virtual R & operator() (int i,int j) =0;
  virtual void call(int ,int ie,int label,void * data) =0;  // 
  const LinearComb<pair<MGauche,MDroit>,C_F0> * bilinearform;
  
  MatriceElementaire & operator()(int k,int ie,int label,void * s=0) {
    call(k,ie,label,s);
    return *this;}
};

template <class FES> 
class MatDataFES { 
public:
  typedef FES FESpace;
  typedef typename FESpace::FElement FElement;

  typedef typename  FESpace::QFElement QFElement; 
  typedef typename  FESpace::QFBorderElement QFBorderElement; 
  CountPointer<const FESpace> cUh,cVh;
  const FESpace &Uh;
  const FESpace &Vh;
  const QFElement & FIT;
  const QFBorderElement & FIE;
  MatDataFES(const FESpace & UUh,const QFElement & fit, const QFBorderElement & fie)
    :Uh(UUh),Vh(UUh),FIT(fit),FIE(fie) {}
  MatDataFES(const FESpace & UUh,const FESpace & VVh,const QFElement & fit, const QFBorderElement & fie)
    :Uh(UUh),Vh(VVh),FIT(fit),FIE(fie) {}
    

};

template <class R,class FES> 
class MatriceElementaireFES :   public MatDataFES<FES> ,   public MatriceElementaire<R> 
{  
public:
  typedef MatriceElementaire<R> MElm ;
  using MElm::Full;
  using MElm::Symmetric;

  typedef typename MElm::TypeOfMatriceElementaire TypeOfMatriceElementaire;
  typedef FES FESpace;

  typedef typename  FESpace::FElement FElement; 
  typedef typename  FESpace::QFElement QFElement; 
  typedef typename  FESpace::QFBorderElement QFBorderElement; 

  MatriceElementaireFES(const FESpace & UUh,const FESpace & VVh,int llga
			,int *nnj,int * nni,TypeOfMatriceElementaire t=Full,
			const QFElement & fit=*QFElement::Default,
			const QFBorderElement & fie =*QFBorderElement::Default) 
                     
    :
    MatDataFES<FES>(UUh,VVh,fit,fie),
    MatriceElementaire<R>(UUh.esize()+VVh.esize(),llga,nnj,nni,t)
  {}
       
  MatriceElementaireFES(const FESpace & UUh,int llga,int *nni,
			TypeOfMatriceElementaire t=Symmetric,
			const QFElement & fit=*QFElement::Default,
			const QFBorderElement & fie =*QFBorderElement::Default)
    :
    MatDataFES<FES>(UUh,UUh,fit,fie),
    MatriceElementaire<R>(UUh.esize(),llga,nni,nni,t)
  {}

  //  for discontinous Galerkine method
  MatriceElementaireFES(const FESpace & UUh,int llga,int *nni,
			int lk,
			TypeOfMatriceElementaire t=Symmetric,
			const QFElement & fit=*QFElement::Default,
			const QFBorderElement & fie =*QFBorderElement::Default) 
    :
    MatDataFES<FES>(UUh,UUh,fit,fie),
    MatriceElementaire<R>(UUh.esize(),llga,nni,lk,t)
  {}
  ~MatriceElementaireFES() {}
  const LinearComb<pair<MGauche,MDroit>,C_F0> * bilinearform;
  
  MatriceElementaireFES & operator()(int k,int ie,int label,void * s=0) {
    this->call(k,ie,label,s);
    return *this;}
};

template <class R,class FES=FESpace> 
class MatriceElementairePleine:public MatriceElementaireFES<R,FES> {

  /* --- stockage --
     //  n = 4 m = 5
     //  0  1  2  3  4
     //  5  6  7  8  9
     // 10 11 12 13 14
     // 15 16 17 18 19
     ------------------*/

public:
  typedef FES FESpace;
  typedef typename  FESpace::Mesh Mesh;
  typedef typename  FESpace::QFElement QFElement;
  typedef typename  FESpace::QFBorderElement QFBorderElement;
  typedef typename  FESpace::FElement FElement; 

  R & operator() (int i,int j) {return this->a[i*this->m+j];}
  // MatPleineElementFunc element;
  void  (* element)(MatriceElementairePleine &,const FElement &,const FElement &, double*,int ie,int label,void *) ; 
   void  (* faceelement)(MatriceElementairePleine &,const FElement &,const FElement &,const FElement &,const FElement &, double*,int ie,int iee, int label,void *) ; 
 void call(int k,int ie,int label,void *);
  
  MatriceElementairePleine & operator()(int k,int ie,int label,void * stack=0)
  {call(k,ie,label,stack);return *this;}
  MatriceElementairePleine(const FESpace & VVh,
                           const QFElement & fit=*QFElement::Default,
                           const QFBorderElement & fie =*QFBorderElement::Default)  
    :MatriceElementaireFES<R,FES>(VVh,
			Square(VVh.MaximalNbOfDF()),
			new int[VVh.MaximalNbOfDF()],this->Full,fit,fie),
    element(0),faceelement(0) {}
 
   //  matrice for VF or Galerkin Discontinus
   MatriceElementairePleine(const FESpace & VVh,bool VF,
                           const QFElement & fit=*QFElement::Default,
                           const QFBorderElement & fie =*QFBorderElement::Default)  
     :MatriceElementaireFES<R,FES>(VVh,
			Square(VVh.MaximalNbOfDF()*2),
			new int[VVh.MaximalNbOfDF()*2],
			VF?VVh.MaximalNbOfDF()*2:0,
			this->Full,fit,fie),			
    element(0),faceelement(0) {}

  MatriceElementairePleine(const FESpace & UUh,const FESpace & VVh,
                               const QFElement & fit=*QFElement::Default,
                               const QFBorderElement & fie =*QFBorderElement::Default) 
    :MatriceElementaireFES<R,FES>(UUh,VVh,
				  UUh.MaximalNbOfDF()*VVh.MaximalNbOfDF(),
				  new int[UUh.MaximalNbOfDF()],
				  new int[VVh.MaximalNbOfDF()],this->Full,fit,fie),
     element(0),faceelement(0) {}

}; 

template <class R,class FES=FESpace> 
class MatriceElementaireSymetrique:public MatriceElementaireFES<R,FES> {



  // --- stockage --
  //   0
  //   1 2
  //   3 4 5
  //   6 7 8 9
  //  10 . . . .
  //

public:
  typedef FES FESpace;
  typedef typename  FESpace::Mesh Mesh;
  typedef typename  FESpace::QFElement QFElement;
  typedef typename  FESpace::QFBorderElement QFBorderElement;
  typedef typename  FESpace::FElement FElement; 

  R & operator()(int i,int j) 
  {return j < i ? this->a[(i*(i+1))/2 + j] : this->a[(j*(j+1))/2 + i] ;}
  void (* element)(MatriceElementaireSymetrique &,const FElement &, double*,int ie,int label,void *) ; 
  void (* mortar)(MatriceElementaireSymetrique &,const FMortar &,void *) ; 
  void call(int k,int ie,int label,void * stack);
  MatriceElementaireSymetrique(const FESpace & VVh,
                               const QFElement & fit=*QFElement::Default,
                               const QFBorderElement & fie =*QFBorderElement::Default) 
    :MatriceElementaireFES<R,FES>(
           VVh,
	   int(VVh.MaximalNbOfDF()*(VVh.MaximalNbOfDF()+1)/2),
	   new int[VVh.MaximalNbOfDF()],this->Symmetric,
       fit,fie),
       element(0),mortar(0) {}
  MatriceElementaireSymetrique & operator()(int k,int ie,int label,void * stack=0) 
  {this->call(k,ie,label,stack);return *this;};
};


template <class R>  class MatriceProfile;

//  classe modele pour matrice creuse
//  ---------------------------------
template <class R> 
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
  virtual ~MatriceCreuse(){}
  virtual R & diag(int i)=0;
  virtual R & operator()(int i,int j)=0;
  virtual R * pij(int i,int j) const =0; // Add FH   
  virtual MatriceMorse<R> *toMatriceMorse(bool transpose=false,bool copy=false) const {return 0;} // not 
  virtual bool addMatTo(R coef,std::map< pair<int,int>, R> &mij,bool trans=false,int ii00=0,int jj00=0,bool cnj=false)=0;
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

  // end ADD

};

template <class R> 
inline ostream& operator <<(ostream& f,const MatriceCreuse<R> & m) 
    {return m.dump(f);}

template <class R> 
KN_<R> & operator/=(KN_<R> & x ,const MatriceProfile<R> & a) ;


enum FactorizationType {
    FactorizationNO=0,
    FactorizationCholeski=1,
    FactorizationCrout=2,
    FactorizationLU=3};

template <class R> 
class MatriceProfile:public MatriceCreuse<R> {
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
    : MatriceCreuse<R>(NbOfDF),L(l),U(u),D(d),pL(pl),pU(pu),
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
  const MatriceProfile du() const  
     {return MatriceProfile(this->n,D,U,pU,0,0);}
  const MatriceProfile u() const   
     {return MatriceProfile(this->n,0,U,pU,0,0);}
  const MatriceProfile ut() const  
    
    
    {return MatriceProfile(this->n,0,0,0,U,pU);}

  void Solve(KN_<R> &x,const KN_<R> &b) const {
    /*if (typefac==0)  code faux // FH   nov 2006
      switch(typefac) 
	{
	FactorizationCholeski: cholesky() ; break;
	FactorizationCrout:   crout(); break;
	FactorizationLU:      LU(); break; 
	}*/ 
    if (&x != &b) x=b;x/=*this;} 
  
  int size() const ;
  ~MatriceProfile();
  //  KN_<R>         operator* (const KN_<R> & ) const ;
  void addMatMul(const KN_<R> &x,KN_<R> &ax) const;
  void addMatTransMul(const KN_<R> &x,KN_<R> &ax) const 
    { this->t().addMatMul(x,ax);}
  MatriceCreuse<R> & operator +=(MatriceElementaire<R> &);
  void operator=(const R & v); // Mise a zero 
  void cholesky(double = EPSILON/8.) const ; //
  void crout(double = EPSILON/8.) const ; //
  void LU(double = EPSILON/8.) const ; //
  R & diag(int i) { return D[i];}
  R & operator()(int i,int j) { if(i!=j) ffassert(0); return D[i];} // a faire 
  R * pij(int i,int j) const { if(i!=j) ffassert(0); return &D[i];} // a faire  Modif FH 31102005
  MatriceMorse<R> *toMatriceMorse(bool transpose=false,bool copy=false) const ;
  
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
    : MatriceCreuse<R>(A.n,A.m,0)
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

  
  bool addMatTo(R coef,std::map< pair<int,int>, R> &mij,bool trans=false,int ii00=0,int jj00=0,bool cnj=false);

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
class MatriceMorse:public MatriceCreuse<R> {
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
   friend class MatriceMorse;
   virtual void Solver(const MatriceMorse<R> &a,KN_<R> &x,const KN_<R> &b) const  =0;
};

  MatriceMorse():MatriceCreuse<R>(0),nbcoef(0),symetrique(true),a(0),lg(0),cl(0),solver(0) {};
  MatriceMorse(KNM_<R> & A, double tol) ;
  MatriceMorse(const int  n,const R *a);
//    :MatriceCreuse<R>(n),solver(0) {}
  MatriceMorse(istream & f);

  template<class FESpace>   explicit 
  MatriceMorse(const FESpace & Uh,bool sym,bool VF=false)
    :MatriceCreuse<R>(Uh.NbOfDF),solver(0) {Build(Uh,Uh,sym,VF);}

  template<class FESpace>   explicit 
  MatriceMorse(const FESpace & Uh,const FESpace & Vh,bool VF=false)
    :MatriceCreuse<R>(Uh.NbOfDF,Vh.NbOfDF,0),solver(0) 
  {Build(Uh,Vh,false,VF);}

  template<class FESpace>   explicit 
  MatriceMorse(const FESpace & Uh,const FESpace & Vh,
               void (*build)(MatriceMorse *,const FESpace & Uh,const FESpace & Vh,void *data),void *data=0
	       )
    :MatriceCreuse<R>(Uh.NbOfDF,Vh.NbOfDF,0),solver(0) 
  {build(this,Uh,Vh,data);           
  }
  MatriceMorse(int nn,int mm,int nbc,bool sym,R *aa,int *ll,int *cc,bool dd, const VirtualSolver * s=0,bool transpose=false )
    :MatriceCreuse<R>(nn,mm,dd && !transpose),
     nbcoef(nbc),
     symetrique(sym), // transpose = true => dummy false (new matrix)
     a(docpyornot(this->dummy,aa,nbc)),
     lg(docpyornot(this->dummy,ll,nn+1)),
     cl(docpyornot(this->dummy,cc,nbc)),
     solver(s)
  { if(transpose) dotransposition(); };
  void Solve(KN_<R> &x,const KN_<R> &b) const;
  int size() const ;
  void addMatMul(const KN_<R> &x,KN_<R> &ax) const;
  void addMatTransMul(const KN_<R> &x,KN_<R> &ax) const;
  MatriceMorse & operator +=(MatriceElementaire<R> &);
  void operator=(const R & v) 
    { for (int i=0;i< nbcoef;i++) a[i]=v;}
  virtual ~MatriceMorse(){ if (!this->dummy) { delete [] a; delete [] cl;delete [] lg;}}
  ostream&  dump(ostream & f) const ;
  R * pij(int i,int j) const ;
  R  operator()(int i,int j) const {R * p= pij(i,j) ;throwassert(p); return *p;}
  R & operator()(int i,int j)  {R * p= pij(i,j) ;throwassert(p); return *p;}
  R & diag(int i)  {R * p= pij(i,i) ;throwassert(p); return *p;}
  
  void SetSolver(const VirtualSolver & s){solver=&s;}
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
  
template<class K>
  MatriceMorse(int nn,int mm, std::map< pair<int,int>, K> & m, bool sym);
  
 template<class RB,class RAB>
 void  prod(const MatriceMorse<RB> & B, MatriceMorse<RAB> & AB);
 
 MatriceMorse<R> *toMatriceMorse(bool transpose=false,bool copy=false) const {
     return new MatriceMorse(this->n,this->m,nbcoef,symetrique,a,lg,cl,copy, solver,transpose);}
  bool  addMatTo(R coef,std::map< pair<int,int>, R> &mij,bool trans=false,int ii00=0,int jj00=0,bool cnj=false);
  

  template<class K>
    MatriceMorse(const MatriceMorse<K> & );

  private:
  void dotransposition ()  ;  // do the transposition 
  CountPointer<const VirtualSolver> solver;
  
  void operator=(const MatriceMorse & );

  template<class FESpace>
  void  Build(const FESpace & Uh,const FESpace & Vh,bool sym,bool VF=false);
  
};




template<class R,class M,class P> 
int ConjuguedGradient(const M & A,const P & C,const KN_<R> &b,KN_<R> &x,const int nbitermax, double &eps,long kprint=1000000000)
{
//  ConjuguedGradient lineare   A*x est appele avec des conditions au limites 
//  non-homogene puis homogene  pour calculer le gradient  
   if (verbosity>50) 
     kprint=2;
   if (verbosity>99)  cout << A << endl;
   throwassert(&x && &b && &A && &C);
   typedef KN<R> Rn;
   int n=b.N();   
   throwassert(n==x.N());
   Rn g(n), h(n), Ah(n), & Cg(Ah);  // on utilise Ah pour stocke Cg  
   g = A*x;  
   double xx= real((x,conj(x)));
   double epsold=eps;
   g -= b;// g = Ax-b
   Cg = C*g; // gradient preconditionne 
   h =-Cg; 
   double g2 = real((Cg,conj(g)));
   if (g2 < 1e-30) 
    { if(verbosity>1)
       cout << "GC  g^2 =" << g2 << " < 1.e-30  Nothing to do " << endl;
     return 2;  }
   double reps2 =eps >0 ?  eps*eps*g2 : -eps; // epsilon relatif 
   eps = reps2;
   for (int iter=0;iter<=nbitermax;iter++)
     {      
       Ah = A*h;
       double hAh =real((h,conj(Ah)));
      // if (Abs(hAh)<1e-30) ExecError("CG2: Matrix non defined, sorry ");
       R ro =  - real((g,conj(h)))/ hAh; // ro optimal (produit scalaire usuel)
       x += ro *h;
       g += ro *Ah; // plus besoin de Ah, on utilise avec Cg optimisation
       Cg = C*g;
       double g2p=g2; 
       g2 = real((Cg,conj(g)));
       if ( !(iter%kprint) && iter && (verbosity>3) )
         cout << "CG:" <<iter <<  "  ro = " << ro << " ||g||^2 = " << g2 << endl; 
       if (g2 < reps2) { 
         if ( !iter && !xx && g2 && epsold >0 ) {  
             // change fo eps converge to fast due to the 
             // penalization of boundary condition.
             eps =  epsold*epsold*g2; 
             if (verbosity>3 )
             cout << "CG converge to fast (pb of BC)  restart: " << iter <<  "  ro = " 
                  << ro << " ||g||^2 = " << g2 << " <= " << reps2 << "  new eps2 =" <<  eps <<endl; 
              reps2=eps;
           } 
         else 
          { 
           if (verbosity>1 )
            cout << "CG converge: " << iter <<  "  ro = " << ro << " ||g||^2 = " << g2 << endl; 
           return 1;// ok 
          }
          }
       double gamma = g2/g2p;       
       h *= gamma;
       h -= Cg;  //  h = -Cg * gamma* h       
     }
   cout << " GC: method doesn't converge in " << nbitermax 
        <<  " iteration , xx= "  << xx<< endl;
   return 0; 
}

template<class R,class M,class P> 
int ConjuguedGradient2(const M & A,const P & C,KN_<R> &x,const KN_<R> &b,const int nbitermax, double &eps,long kprint=1000000000)
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
   R g2 = (Cg,g);
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
       R hAh =(h,Ah);
       if (norm(hAh)<1e-60) ExecError("CG2: Matrix is not defined (/0), sorry ");
       ro =  - (g,h)*rop/hAh ; // ro optimal (produit scalaire usuel)
       x += (ro-rop) *h;
       g += (ro/rop) *Ah; // plus besoin de Ah, on utilise avec Cg optimisation
       Cg = C*g;
       R g2p=g2; 
       g2 = (Cg,g);
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

template <class R> 
class MatriceIdentite: VirtualMatrice<R> { public:
 typedef typename VirtualMatrice<R>::plusAx plusAx;
    MatriceIdentite() :VirtualMatrice<R>(0) {}; 
 void addMatMul(const  KN_<R>  & x, KN_<R> & Ax) const { 
     ffassert(x.N()==Ax.N());
   Ax+=x; } 
 plusAx operator*(const KN<R> &  x) const {return plusAx(this,x);} 
  bool ChecknbLine(int n) const { return true;}  
  bool ChecknbColumn(int m) const { return true;} 

};  

template<class R>
class SolveGCDiag :   public MatriceMorse<R>::VirtualSolver , public VirtualMatrice<R>{
  int n;
  int nbitermax;
  double eps;
  mutable double  epsr;
  KN<R> D1;
  public:
  typedef typename VirtualMatrice<R>::plusAx plusAx;
  SolveGCDiag(const MatriceMorse<R> &A,int itmax,double epsilon=1e-6) : 
    VirtualMatrice<R>(A.n),
    n(A.n),nbitermax(itmax?itmax: Max(100,n)),eps(epsilon),epsr(0),D1(n)
  { throwassert(A.sym());
    for (int i=0;i<n;i++)
      D1[i] = 1./A(i,i);}
   void Solver(const MatriceMorse<R> &a,KN_<R> &x,const KN_<R> &b) const  {
     epsr = (eps < 0) ? (epsr >0 ? -epsr : -eps ) : eps ;
    // cout << " epsr = " << epsr << endl;
     ConjuguedGradient<R,MatriceMorse<R>,SolveGCDiag<R> >(a,*this,b,x,nbitermax,epsr);
   }
plusAx operator*(const KN_<R> &  x) const {return plusAx(this,x);} 


 void addMatMul(const KN_<R> & x, KN_<R> & Ax) const 
  {  ffassert(x.N()==Ax.N());
     for (int i=0;i<n;i++) 
     Ax[i]+= D1[i]*x[i];}
     
   bool ChecknbLine(int nn) const { return n==nn;}  
   bool ChecknbColumn(int mm) const { return n==mm;} 

};

// add Sep 2007 for generic Space solver
typedef MatriceMorse<double>::VirtualSolver *
   (*SparseRMatSolve)(const MatriceMorse<double> *A,int strategy,
		   double ttgv, double epsilon, double pivot,double pivot_sym,
		   int NbSpace,int itmax , const void * precon, void * stack);
typedef MatriceMorse<Complex>::VirtualSolver *
   (*SparseCMatSolve)(const MatriceMorse<Complex> *A,int strategy,
		   double ttgv, double epsilon, double pivot,double pivot_sym,
		   int NbSpace,int itmax , const void * precon, void * stack);


template<class R> struct DefSparseSolver {
    typedef typename MatriceMorse<R>::VirtualSolver * 
            (*SparseMatSolver)(const MatriceMorse<R> *A,int strategy,
			       double ttgv, double epsilon, double pivot,double pivot_sym ,
			       int NbSpace,int itmax , const void * precon, void * stack);
    static SparseMatSolver solver;
    
  static  typename MatriceMorse<R>::VirtualSolver * 
      Build(const MatriceMorse<R> *A,int strategy,double tgv, double eps, double tol_pivot,double tol_pivot_sym,
	    int NbSpace,int itmax ,const  void * precon, void * stack)

    {
      typename MatriceMorse<R>::VirtualSolver *ret=0;
	if(solver)
	    ret =(solver)(A,strategy,tgv,eps,tol_pivot,tol_pivot_sym,NbSpace,itmax,(const void *) precon,stack);
	return ret;	
    }
};

// End Sep 2007 for generic Space solver

struct TypeSolveMat {
    enum TSolveMat { NONESQUARE=0, LU=1, CROUT=2, CHOLESKY=3, GC = 4 , GMRES = 5, SparseSolver=6 };
    TSolveMat t;
    bool sym;
    bool profile;
    TypeSolveMat(TSolveMat tt=LU) :t(tt),
	sym(t == CROUT || t ==CHOLESKY  ||  t==GC ),
	profile(t != GC && t != GMRES && t != NONESQUARE && t != SparseSolver ) {}
    bool operator==(const TypeSolveMat & a) const { return t == a.t;}                               
    bool operator!=(const TypeSolveMat & a) const { return t != a.t;}
    static TSolveMat defaultvalue;
};


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


#ifdef HAVE_LIBUMFPACK__XXXXXXXXXXXX
template<class R>
class SolveUMFPack :   public MatriceMorse<R>::VirtualSolver  {
  double eps;
  mutable double  epsr;
  double tgv;
  void *Symbolic, *Numeric ;
  int umfpackstrategy;
  double tol_pivot_sym,tol_pivot; //Add 31 oct 2005
public:
  SolveUMFPack(const MatriceMorse<R> &A,int strategy,double ttgv, double epsilon=1e-6,
	       double pivot=-1.,double pivot_sym=-1.  ) : 
    eps(epsilon),epsr(0),
    tgv(ttgv),
    Symbolic(0),Numeric(0)  ,
    umfpackstrategy(strategy),
    tol_pivot_sym(pivot_sym),tol_pivot(pivot)
  { 
    
    int status;
    throwassert( !A.sym() && Numeric == 0 && Symbolic==0 );
    int n=A.n;
    double Control[UMFPACK_CONTROL];
    double Info[UMFPACK_INFO];
    
    for(int i=0;i<UMFPACK_CONTROL;i++) Control[i]=0;
    for(int i=0;i<UMFPACK_INFO;i++) Info[i]=0;
    
    umfpack_di_defaults (Control) ;
    Control[UMFPACK_PRL]=1;
   // Control[UMFPACK_PIVOT_TOLERANCE]=1E-10;
    
    if(verbosity>4) Control[UMFPACK_PRL]=2;
    if(tol_pivot_sym>0) Control[UMFPACK_SYM_PIVOT_TOLERANCE]=pivot_sym;
    if(tol_pivot>0) Control[UMFPACK_PIVOT_TOLERANCE]=pivot;
    if(umfpackstrategy>=0)   Control[UMFPACK_STRATEGY]=umfpackstrategy;
    if(verbosity>3) { 
      cout << "  UMFpack real  Solver Control :" ;
      cout << "\n\t SYM_PIVOT_TOLERANCE "<< Control[UMFPACK_SYM_PIVOT_TOLERANCE];
      cout << "\n\t PIVOT_TOLERANCE     "<< Control[UMFPACK_PIVOT_TOLERANCE];
      cout << "\n\t PRL                 "<< Control[UMFPACK_PRL];
      cout << "\n";      
    }
    
    status = umfpack_di_symbolic (n, n, A.lg, A.cl, A.a, &Symbolic,Control,Info) ;
    if (status !=  0)
    {
      (void) umfpack_di_report_matrix (n, n, A.lg, A.cl, A.a, 1, Control) ;

	umfpack_di_report_info (Control, Info) ;
	umfpack_di_report_status (Control, status) ;
	cerr << "umfpack_di_symbolic failed" << endl;
	ExecError("umfpack_di_symbolic failed");
	//ffassert(0);
    }

    status = umfpack_di_numeric (A.lg, A.cl, A.a, Symbolic, &Numeric,Control,Info) ;
    if (status !=  0)
    {
	umfpack_di_report_info (Control, Info) ;
	umfpack_di_report_status (Control, status) ;
	cerr << "umfpack_di_numeric failed" << endl;
	ExecError("umfpack_di_numeric failed");
	ffassert(0);
    }

    if (Symbolic) umfpack_di_free_symbolic (&Symbolic),Symbolic=0; 
    if(verbosity>3)
    cout << "  -- umfpack_di_build LU " << n <<  endl;
    if(verbosity>5)     (void)  umfpack_di_report_info(Control,Info);

  }
  void Solver(const MatriceMorse<R> &A,KN_<R> &x,const KN_<R> &b) const  {
    ffassert ( &x[0] != &b[0]);
    epsr = (eps < 0) ? (epsr >0 ? -epsr : -eps ) : eps ;
    // cout << " epsr = " << epsr << endl;
    double Control[UMFPACK_CONTROL];
    double Info[UMFPACK_INFO];
    for(int i=0;i<UMFPACK_CONTROL;i++) Control[i]=0;
    for(int i=0;i<UMFPACK_INFO;i++) Info[i]=0;
    int n= b.N(); 
     ffassert(A.ChecknbLine( n) && n == x.N() && A.ChecknbColumn(n) );
    
     umfpack_di_defaults (Control) ;
     // change UMFPACK_At to UMFPACK_Aat in complex 
    int status = umfpack_di_solve (UMFPACK_Aat, A.lg, A.cl, A.a, x, b, Numeric,Control,Info) ;
    if (status != 0)
    {
	umfpack_di_report_info (Control, Info) ;
	umfpack_di_report_status (Control, status) ;
	cerr << "umfpack_di_solve failed" << endl;
	ExecError("umfpack_di_solve failed");
	
	ffassert(0);
    }
     if(verbosity>2)
    cout << " -- umfpack_di_solve " << endl;
    if(verbosity>3)
    cout << "   b min max " << b.min() << " " <<b.max() << endl;
    if(verbosity>3)     (void)  umfpack_di_report_info(Control,Info);
     if(verbosity>1) cout << "   x min max " << x.min() << " " <<x.max() << endl;
  }

  ~SolveUMFPack() { 
   if(verbosity>3)
    cout << "~SolveUMFPack S:" << Symbolic << " N:" << Numeric <<endl;
    if (Symbolic)   umfpack_di_free_symbolic  (&Symbolic),Symbolic=0; 
    if (Numeric)    umfpack_di_free_numeric (&Numeric),Numeric=0;
  }
  void addMatMul(const KN_<R> & x, KN_<R> & Ax) const 
  {  
    ffassert(x.N()==Ax.N());
    Ax +=  (const MatriceMorse<R> &) (*this) * x; 
  }
     
}; 

template<>
class SolveUMFPack<Complex> :   public MatriceMorse<Complex>::VirtualSolver  {
  double eps;
  mutable double  epsr;
  int umfpackstrategy;
  double tgv;
  void *Symbolic, *Numeric ;
  double *ar,*ai;


    double tol_pivot_sym,tol_pivot; //Add 31 oct 2005

public:
  SolveUMFPack(const MatriceMorse<Complex> &A,int strategy,double ttgv, double epsilon=1e-6,
     double pivot=-1.,double pivot_sym=-1.
) : 
    eps(epsilon),epsr(0),umfpackstrategy(strategy),tgv(ttgv),
    Symbolic(0),Numeric(0),
    ar(0),ai(0),
    tol_pivot_sym(pivot_sym), 
    tol_pivot(pivot)
   { 
    int status;
    throwassert( !A.sym());
    int n=A.n;
    //  copy the coef of the matrice ---
     ar= new double[A.nbcoef];
     ai= new double[A.nbcoef];
     ffassert(ar && ai);
     C2RR(A.nbcoef,A.a,ar,ai);
        
    double Control[UMFPACK_CONTROL];
    double Info[UMFPACK_INFO];
    umfpack_zi_defaults (Control) ;
    Control[UMFPACK_PRL]=1;
    if(verbosity>4) Control[UMFPACK_PRL]=2;
   //    Control[UMFPACK_SYM_PIVOT_TOLERANCE]=1E-10;
  //  Control[UMFPACK_PIVOT_TOLERANCE]=1E-10;
    if(tol_pivot_sym>0) Control[UMFPACK_SYM_PIVOT_TOLERANCE]=pivot_sym;
    if(tol_pivot>0) Control[UMFPACK_PIVOT_TOLERANCE]=pivot;
    if(umfpackstrategy>=0) Control[UMFPACK_STRATEGY]=umfpackstrategy;
    if(verbosity>3) { 
      cout << "  UMFpack complex Solver Control :" ;
      cout << "\n\t SYM_PIVOT_TOLERANCE "<< Control[UMFPACK_SYM_PIVOT_TOLERANCE];
      cout << "\n\t PIVOT_TOLERANCE     "<< Control[UMFPACK_PIVOT_TOLERANCE];
      cout << "\n\t PRL                 "<< Control[UMFPACK_PRL];
      cout << "\n";      
    }
    status = umfpack_zi_symbolic (n, n, A.lg, A.cl, ar,ai, &Symbolic,Control,Info) ;
    if (status < 0)
    {
      (void) umfpack_zi_report_matrix (n, n, A.lg, A.cl, ar,ai, 1, Control) ;

	umfpack_zi_report_info (Control, Info) ;
	umfpack_zi_report_status (Control, status) ;
	cerr << "umfpack_zi_symbolic failed" << endl;
	ExecError("umfpack_zi_symbolic failed");
	ffassert(0);
	exit(2);
    }

    status = umfpack_zi_numeric (A.lg, A.cl, ar,ai, Symbolic, &Numeric,Control,Info) ;
    if (status < 0)
    {
	umfpack_zi_report_info (Control, Info) ;
	umfpack_zi_report_status (Control, status) ;
	cerr << "umfpack_zi_numeric failed" << endl;
	ExecError("umfpack_zi_numeric failed");
	ffassert(0);
	exit(2);
    }

    if (Symbolic) umfpack_zi_free_symbolic (&Symbolic),Symbolic=0; 
    if(verbosity>3)
    cout << "umfpack_zi_build LU " << n <<  endl;
    if(verbosity>5)     (void)  umfpack_zi_report_info(Control,Info);

  }
  void Solver(const MatriceMorse<Complex> &A,KN_<Complex> &x,const KN_<Complex> &b) const  {
        ffassert ( &x[0] != &b[0]);
    epsr = (eps < 0) ? (epsr >0 ? -epsr : -eps ) : eps ;
    // cout << " epsr = " << epsr << endl;
    double Control[UMFPACK_CONTROL];
    double Info[UMFPACK_INFO];
     umfpack_zi_defaults (Control) ;
     int n = b.N();
     ffassert(A.ChecknbLine( n) && n == x.N() && A.ChecknbColumn(n) );
     KN<double> xr(n),xi(n),br(n),bi(n);
     C2RR(n,b,br,bi);
     // change UMFPACK_At to UMFPACK_Aat in complex  oct 2005 
    int status = umfpack_zi_solve (UMFPACK_Aat, A.lg, A.cl, ar,ai, xr, xi, br,bi, Numeric,Control,Info) ;
    if (status < 0)
    {
	umfpack_zi_report_info (Control, Info) ;
	umfpack_zi_report_status (Control, status) ;
	cerr << "umfpack_zi_solve failed" << endl;
	ExecError("umfpack_zi_numeric failed");
	ffassert(0);
	exit(2);
    }
    RR2C(n,xr,xi,x);
    if(verbosity>1)
    {
     cout << "  -- umfpack_zi_solve " << endl;
     if(verbosity>3)     (void)  umfpack_zi_report_info(Control,Info);
    
      cout << "   b min max " << b.min() << " " <<b.max() << endl;
      cout << "   x min max " << x.min() << " " <<x.max() << endl;
    }
  }

  ~SolveUMFPack() { 
    if(verbosity>5)
    cout << "~SolveUMFPack " << endl;
    if (Symbolic)   umfpack_zi_free_symbolic  (&Symbolic),Symbolic=0; 
    if (Numeric)    umfpack_zi_free_numeric (&Numeric),Numeric=0;
    delete [] ar;
    delete [] ai;   
  }
  void addMatMul(const KN_<Complex> & x, KN_<Complex> & Ax) const 
  {  
    ffassert(x.N()==Ax.N());
    Ax +=  (const MatriceMorse<Complex> &) (*this) * x; 
  }
     

}; 


inline MatriceMorse<double>::VirtualSolver *
BuildSolverUMFPack(const MatriceMorse<double> *A,int strategy,double tgv, double eps, double tol_pivot,double tol_pivot_sym,
		   int NbSpace,int itmax , const void * precon, void * stack )
{
    return new SolveUMFPack<double>(*A,strategy,tgv,eps,tol_pivot,tol_pivot_sym);
}

inline MatriceMorse<Complex>::VirtualSolver *
BuildSolverUMFPack(const MatriceMorse<Complex> *A,int strategy,double tgv, double eps, double tol_pivot,double tol_pivot_sym,
		   int NbSpace,int itmax , const void * precon, void * stack )
{
    return new SolveUMFPack<Complex>(*A,strategy,tgv,eps,tol_pivot,tol_pivot_sym);
}




#endif  
//  endif ---- UMFPACK ----
#endif
