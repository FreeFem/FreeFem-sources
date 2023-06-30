#include "rgraph.hpp"
#include "RNM.hpp"
#include "fem.hpp"
#include "FESpacen.hpp"
#include "FESpace.hpp"
extern long verbosity ;

template<class I,class R> class LinearComb;
class MGauche ;
class MDroit;
class C_F0 ;
using FESpace=EF23::FESpace;


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
  const int lnki,lnkj; // size of the 4 next array
  int *nik,*nikk;  //  number of df in element   k,kk for VF and GD methode 
  int *njk,*njkk;  //  number of df in element   k,kk  for VF and GD methode
  const int optim;

  MatriceElementaire(int datasize,int llga
                     ,int *nnj,int * nni,TypeOfMatriceElementaire t=Full,int ooptim=1)
    
    :   lga(llga),a(new R[lga]),
        ni(nni),nj(nnj),n(0),m(0),mtype(t),data(datasize),
        onFace(false),lnki(0),lnkj(0),nik(0),nikk(0),njk(0),njkk(0),optim(ooptim)
        {}
       

 //  for discontinuous Galerkine method
  MatriceElementaire(int datasize,int llga,int *nni,
                     int lk,
                     TypeOfMatriceElementaire t=Symmetric, int ooptim=1
                     ) 
    :
    lga(llga),a(new R[lga]),
    ni(nni),nj(nni),n(0),m(0),mtype(t),data(datasize*(lk?2:1)) ,
       onFace(lk!=0),
       lnki(lk),lnkj(lk),
       nik(lk? new int[lk*2]:0),
       nikk(nik+lk),
       njk(nik),
       njkk(nik+lk),
        optim(ooptim)
       { ffassert(lk>=0);}

    //  for discontinuous Galerkine method
    MatriceElementaire(int datasize,int llga,int *nni,int lki,int *nnj,int lkj,
                       TypeOfMatriceElementaire t=Full,int ooptim=1
                       )
    :
    lga(llga),a(new R[lga]),
    ni(nni),nj(nnj),n(0),m(0),mtype(t),data(datasize*(lki+lkj?2:1)) ,
    onFace(lki+lkj),
    lnki(lki),lnkj(lkj),
    nik(lki? new int[lki*2]:0),
    nikk(nik+lki),
    njk(lkj? new int[lkj*2]:0),
    njkk(njk+lkj),
    optim(ooptim)
    {  ffassert(lki>=0);}// non teste ??? .... F. hecht ...

  virtual ~MatriceElementaire() {
    if(ni != nj) 
      delete [] nj;
         
    delete [] ni;
    delete [] a;
    if ( nik) delete[] nik;
     }
      
  virtual R & operator() (int i,int j) =0;
    virtual void call(int ,int ie,int label,void * data,void *Q=0) =0;  //
  const LinearComb<pair<MGauche,MDroit>,C_F0> * bilinearform;
  
  MatriceElementaire & operator()(int k,int ie,int label,void * s=0,void *B=0) {
    call(k,ie,label,s,B);
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
			const QFBorderElement & fie =*QFBorderElement::Default,int optim=1)
                     
    :
    MatDataFES<FES>(UUh,VVh,fit,fie),
    MatriceElementaire<R>(UUh.esize()+VVh.esize(),llga,nnj,nni,t,optim)
  {}
       
  MatriceElementaireFES(const FESpace & UUh,int llga,int *nni,
			TypeOfMatriceElementaire t=Symmetric,
			const QFElement & fit=*QFElement::Default,
			const QFBorderElement & fie =*QFBorderElement::Default,int optim=1)
    :
    MatDataFES<FES>(UUh,UUh,fit,fie),
    MatriceElementaire<R>(UUh.esize(),llga,nni,nni,t,optim)
  {}

  //  for discontinuous Galerkine method
  MatriceElementaireFES(const FESpace & UUh,int llga,int *nni,
			int lk,
			TypeOfMatriceElementaire t=Symmetric,
			const QFElement & fit=*QFElement::Default,
			const QFBorderElement & fie =*QFBorderElement::Default,int optim=1)
    :
    MatDataFES<FES>(UUh,UUh,fit,fie),
    MatriceElementaire<R>(UUh.esize(),llga,nni,lk,t,optim)
  {}
    
    MatriceElementaireFES(const FESpace & UUh,const FESpace & VVh,int llga
                          ,int *nnj,int lkj,int * nni,int lki,TypeOfMatriceElementaire t=Full,
                          const QFElement & fit=*QFElement::Default,
                          const QFBorderElement & fie =*QFBorderElement::Default,int optim=1)
    
    :
    MatDataFES<FES>(UUh,VVh,fit,fie),
    MatriceElementaire<R>(UUh.esize()+VVh.esize(),llga,nni,lki,nnj,lkj,t,optim)// correct swap i,j may 2023 
    {}
    
  ~MatriceElementaireFES() {}
  const LinearComb<pair<MGauche,MDroit>,C_F0> * bilinearform;
  
  MatriceElementaireFES & operator()(int k,int ie,int label,void * s=0,void *Q=0) {
    this->call(k,ie,label,s,Q);
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
  typedef typename  FESpace::Mesh::Rd Rd;
    

  R & operator() (int i,int j) {return this->a[i*this->m+j];}
  // MatPleineElementFunc element;
  void  (* element)(MatriceElementairePleine &,const FElement &,const FElement &, double*,int ie,int label,void *,Rd *) ;
  void  (* faceelement)(MatriceElementairePleine &,const FElement &,const FElement &,const FElement &,const FElement &, double*,int ie,int iee, int label,void *,Rd *) ;
    void call(int k,int ie,int label,void *,void *B);
  
  MatriceElementairePleine & operator()(int k,int ie,int label,void * stack=0,Rd *Q=0)
  {call(k,ie,label,stack,Q);return *this;}
  MatriceElementairePleine(const FESpace & VVh,
                           const QFElement & fit=*QFElement::Default,
                           const QFBorderElement & fie =*QFBorderElement::Default,int optim=1)
    :MatriceElementaireFES<R,FES>(VVh,
			Square(VVh.MaximalNbOfDF()),
			new int[VVh.MaximalNbOfDF()],this->Full,fit,fie,optim),
    element(0),faceelement(0) {}
 
   //  matrice for VF or Galerkin Discontinus
   MatriceElementairePleine(const FESpace & VVh,bool VF,
                           const QFElement & fit=*QFElement::Default,
                           const QFBorderElement & fie =*QFBorderElement::Default,int optim=1)
     :MatriceElementaireFES<R,FES>(VVh,
                                   RNM::Square(VVh.MaximalNbOfDF()*2),
			new int[VVh.MaximalNbOfDF()*2],
			VF?VVh.MaximalNbOfDF()*2:0,
                                   this->Full,fit,fie,optim),
    element(0),faceelement(0) {}

  MatriceElementairePleine(const FESpace & UUh,const FESpace & VVh,
                               const QFElement & fit=*QFElement::Default,
                               const QFBorderElement & fie =*QFBorderElement::Default,int optim=1)
    :MatriceElementaireFES<R,FES>(UUh,VVh,
				  UUh.MaximalNbOfDF()*VVh.MaximalNbOfDF(),
				  new int[UUh.MaximalNbOfDF()],
				  new int[VVh.MaximalNbOfDF()],this->Full,fit,fie,optim),
     element(0),faceelement(0) {}

    MatriceElementairePleine(const FESpace & UUh,const FESpace & VVh,bool VF,
                             const QFElement & fit=*QFElement::Default,
                             const QFBorderElement & fie =*QFBorderElement::Default,int optim=1)
    :MatriceElementaireFES<R,FES>(UUh,VVh,
                                  UUh.MaximalNbOfDF()*VVh.MaximalNbOfDF()*4,
                                  new int[UUh.MaximalNbOfDF()*2],VF?UUh.MaximalNbOfDF()*2:0,
                                  new int[VVh.MaximalNbOfDF()*2],VF?VVh.MaximalNbOfDF()*2:0,this->Full,fit,fie,optim),
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
    typedef void * FMortar ;
  typedef FES FESpace;
  typedef typename  FESpace::Mesh Mesh;
  typedef typename  FESpace::QFElement QFElement;
  typedef typename  FESpace::QFBorderElement QFBorderElement;
  typedef typename  FESpace::FElement FElement; 
  typedef typename  FESpace::Mesh::Rd Rd;
  R & operator()(int i,int j) 
  {return j < i ? this->a[(i*(i+1))/2 + j] : this->a[(j*(j+1))/2 + i] ;}
  void (* element)(MatriceElementaireSymetrique &,const FElement &, double*,int ie,int label,void *,Rd *) ;
  void (* mortar)(MatriceElementaireSymetrique &,const FMortar &,void *) ;
  void call(int k,int ie,int label,void * stack,void *B);
  MatriceElementaireSymetrique(const FESpace & VVh,
                               const QFElement & fit=*QFElement::Default,
                               const QFBorderElement & fie =*QFBorderElement::Default,int optim=1)
    :MatriceElementaireFES<R,FES>(
           VVh,
	   int(VVh.MaximalNbOfDF()*(VVh.MaximalNbOfDF()+1)/2),
	   new int[VVh.MaximalNbOfDF()],this->Symmetric,
       fit,fie,optim),
       element(0),mortar(0) {}
  MatriceElementaireSymetrique & operator()(int k,int ie,int label,void * stack=0,Rd *B=0)
  {this->call(k,ie,label,stack,B);return *this;};
};




//  classe modele pour matrice creuse
//  ---------------------------------
template <class R>  class MatriceElementaire;

