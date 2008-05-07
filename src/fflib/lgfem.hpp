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
// ATTENTION pfes est la classe qui modelise un pointeur
// sur un FESpace (donc un espace d'element fini 
//  mais si la maillage change alors 
//  l'espace change 
//   cf la fonction set qui reconstruit FESpace


void init_lgmat(); // initialisation for sparse mat functionnallity

class v_fes; 
typedef v_fes *pfes;

namespace  Fem2D {
  class Mesh3;

}
using Fem2D::Mesh;
using Fem2D::Mesh3;

typedef Mesh * pmesh;
typedef Mesh3 * pmesh3;

using  Fem2D::FESpace;
using  Fem2D::TypeOfFE;
using  Fem2D::R;
//  typedef double R;
namespace {
  using namespace Fem2D;
  using  Fem2D::Vertex;
  
  class lgVertex {
  public:
    typedef double R;
    CountPointer<Mesh> pTh;
    Vertex *v;
    void Check() const {   if (!v || !pTh) { ExecError("Too bad! Unset Vertex!"); } }
    void init() { v=0;pTh.init();}
    lgVertex(Mesh * Th,long kk): pTh(Th),v( &(*pTh)(kk)) {}
    lgVertex(Mesh * Th,Vertex * kk): pTh(Th),v(kk) {}
    operator int() const { Check(); return (* pTh)(v);} 
    operator R2*(){ Check(); return v;} 
    R x() const {Check() ; return v->x;}
    R y() const {Check() ; return v->y;}
    //  R z() const {Check() ; return v->z;}
    long lab() const {Check() ; return v->lab;}
    void destroy()  {pTh.destroy();}
  };
  
  class lgElement { public:
      CountPointer<Mesh> pTh;
    Triangle *k;
    
  lgElement():  k(0) {}
  void  Check() const  {   if (!k || !pTh) { ExecError("Unset Triangle,Sorry!"); } }
  void init() { k=0;pTh.init();}
  void destroy() {pTh.destroy();}
  lgElement(Mesh * Th,long kk): pTh(Th),k( &(*pTh)[kk]) {}
  lgElement(Mesh * Th,Triangle * kk): pTh(Th),k(kk) {}
  operator int() const { Check(); return (* pTh)(k);} 
  lgVertex operator [](const long & i) const { Check(); return lgVertex(pTh,&(*k)[i]);}   
  long lab() const {Check() ; return k ? k->lab : 0;}
  double area() const {Check() ; return k->area ;}
  long n() const { return k ? 3: 0 ;}

};

 } // end namespace blanc

void GetPeriodic(Expression perio,    int & nbcperiodic ,    Expression * &periodic);

bool BuildPeriodic( 
  int nbcperiodic,
  Expression *periodic,
  const Mesh &Th,Stack stack,
  int & nbdfv, KN<int> & ndfv,int & nbdfe, KN<int> & ndfe);
  
class v_fes : public RefCounter { public:
    const int d;
  const int N;
  const pmesh* ppTh; // adr du maillage
  CountPointer<FESpace>  pVh;
  Stack stack; // the stack is use whith periodique expression
  
  int nbcperiodic;
  Expression *periodic;

  
  operator FESpace * ()  { 
     throwassert(this);
    if  ( !pVh || *ppTh !=  &pVh->Th )
       pVh=CountPointer<FESpace>(update(),true);
    return  pVh   ;} 
  
  FESpace * update() ;
  
  
  v_fes(int NN,const pmesh* t,Stack s, int n,Expression *p)
    : d(2),N(NN), ppTh(t),pVh(0),stack(s), nbcperiodic(n),periodic(p) {}
  v_fes(int NN,const v_fes *f,Stack s,int n,Expression *p) 
    :  d(2),N(NN),ppTh(f->ppTh),pVh(0),stack(s), nbcperiodic(n),periodic(p)
    {}
  void destroy(){ ppTh=0;pVh=0; delete this;}
  virtual ~v_fes() {}
  bool buildperiodic(Stack stack,int & nbdfv, KN<int> & ndfv,int & nbdfe, KN<int> & ndfe)  ;
  virtual  FESpace * buildupdate(int & nbdfv, KN<int> & ndfv,int & nbdfe, KN<int> & ndfe) =0;
  virtual  FESpace * buildupdate() =0;
 
};  

  
class pfes_tef : public v_fes { public:

  const TypeOfFE * tef ;  
  pfes_tef(const pmesh* t,const TypeOfFE * tt,Stack s=NullStack, int n=0,Expression *p=0 ) 
    : v_fes(tt->N,t,s,n,p),tef(tt) { operator FESpace * ();} 
   FESpace * buildupdate(int & nbdfv, KN<int> & ndfv,int & nbdfe, KN<int> & ndfe) 
    { return  new FESpace(**ppTh,*tef,nbdfv,(int *) ndfv,nbdfe,(int*)ndfe);}
   FESpace * buildupdate()   {  return  new FESpace(**ppTh,*tef);}
  
 };
 
class pfes_tefk : public v_fes { public:

  const TypeOfFE ** tef ;
  const int k;  
  pfes_tefk(const pmesh* t,const TypeOfFE ** tt,int kk,Stack s=NullStack,int n=0,Expression *p=0 ) 
    : v_fes(sum(tt,&Fem2D::TypeOfFE::N,kk),t,s,n,p),tef(tt),k(kk)  { 
    // cout << "pfes_tefk const" << tef << " " << this << endl; 
     operator FESpace * ();} 
  FESpace * buildupdate() { 
  // cout << "pfes_tefk upd:" << tef << " " << this <<  endl; 
   assert(tef);
   return  new FESpace(**ppTh,tef,k);}
   virtual ~pfes_tefk() { delete [] tef;}
   FESpace * buildupdate(int & nbdfv, KN<int> & ndfv,int & nbdfe, KN<int> & ndfe) 
   {
   assert(tef);
   return  new FESpace(**ppTh,tef,k,nbdfv,ndfv,nbdfe,ndfe);}


 }; 
 
class pfes_fes : public v_fes { public:

  pfes * Vh;
  int n;
  pfes_fes( pfes * Vhh, int nn,Stack s=NullStack,int n=0,Expression *p=0) :v_fes((**Vhh).N*nn,static_cast<const v_fes *>(*Vhh),s,n,p),Vh(Vhh),n(nn)  
    { operator FESpace * ();}; 
  FESpace * buildupdate() {  return  new FESpace(***Vh,n);  }
  FESpace * buildupdate(int & nbdfv, KN<int> & ndfv,int & nbdfe, KN<int> & ndfe) 
  {
   InternalError(" No way to define a periodic BC in this case: tensorisation of FEspace ");
  //  return  new FESpace(***Vh,n,nbdfv,ndfv,nbdfe,ndfe);
  }
 };
 
template<class K> class FEbase;
template<class K> class FEcomp { public:
 friend class FEbase<K>;
 FEbase<K> * base;
 int comp;
 FEcomp(FEbase<K> * b,int c) :base(b),comp(c) {};
  private: // rule of programming 
     FEcomp(const FEcomp &);
     void operator= (const FEcomp &); 
};

template<class K>
class FEbase { public:
  v_fes *const*pVh; // pointeur sur la variable stockant FESpace;
  KN<K> * xx; // value
  CountPointer<FESpace> Vh; // espace courant 
  
  KN<K> *x() {return xx;}
  
  FEbase(const pfes  *ppVh) :pVh(ppVh), xx(0),Vh(0) {}
  
  ~FEbase() { delete xx;}  
  void destroy() { // cout << "~FEbase  destroy " << this << endl; 
    delete this;}
  
  void operator=( KN<K> *y) { Vh=**pVh; 
  throwassert((bool) Vh);
  if (xx) delete xx;xx=y;
  throwassert( y->N() == Vh->NbOfDF);}
  FESpace * newVh() { 
    throwassert(pVh  );
    const pfes pp= *pVh;
    // cout << pVh << " " << *pVh << endl;
    return *pp;}  
  
  operator  FESpace &() { throwassert(Vh); return *Vh;}

private: // rule of programming 
  FEbase(const FEbase &);
  void operator= (const FEbase &); 
};



template<class K>
class FEbaseArray { public:
  int N;
  FEbase<K>  **xx;
  FEbaseArray(const pfes  *ppVh,int NN) :N(NN),xx(new FEbase<K> * [NN])
    {
     for (int i=0;i<N;i++)
      xx[i]=new FEbase<K>(ppVh);
    }
 ~FEbaseArray() { 
     //  cout << " ~FEbaseArray " << endl;
      for (int i=0;i<N;i++)
       xx[i]->destroy();
      delete [] xx;} 
  void destroy() { //cout << " destroy ~FEbaseArray " << endl; 
      delete this;}         
 FEbase<K>** operator[](int i)  {
   if(xx==0 || i <0 || i>=N) 
     ExecError("Out of bound in FEbaseArray");
 return xx+i;}  
  private: // rule of programming 
     FEbaseArray(const FEbaseArray &);
     void operator= (const FEbaseArray &); 
};

void GetPeriodic(Expression perio,    int & nbcperiodic ,    Expression * &periodic);
int GetPeriodic(Expression  bb, Expression & b,Expression & f);

/*
template<class K>
class FE  { public:

  const pfes  *pVh; // pointeur sur la variable stockant FESpace;
  CountPointer<FESpace> Vh; // espace courant 
  virtual KN<K> *x() {return xx;}
  KN<K> * xx; // value
  int componante; 
  FE(const pfes  *ppVh,int comp=0) :pVh(ppVh), xx(0),Vh(0),componante(comp) {}
  void destroy() { delete this;}
  virtual ~FE() { delete xx;}
  virtual  void operator=( KN<K> *y) { Vh=**pVh; 
       throwassert((bool) Vh);
       if (xx) delete xx;xx=y;
       throwassert( y->N() == Vh->NbOfDF);}
 virtual  FESpace * newVh() { throwassert(*pVh);return **pVh;}
 virtual FE<K> * set()  {return this;};
 operator  FESpace &() {set(); throwassert(Vh); return *Vh;}
   
   private:
     FE(const FE &);
     void operator= (const FE &); 
};
// bof bof pour test 

template<class K,int comp>
class FE_ : public FE<K>{public:

  FE<K>  ** pere;
  FE_(FE<K> ** p) : FE<K>(0,comp),pere(p) {}
  operator FE<K>  * () { return new FE<K>(*pere,comp);}
  void destroy() { delete this;}    
  FE<K> * set() {throwassert(pere && *pere); xx=(**pere).xx;pVh=(**pere).pVh; Vh =(**pere).Vh;return this; }
  FESpace * newVh() { return **(**pere).pVh;}
  virtual ~FE_() { xx=0; }// pour ne pas  detruire les valeurs qui seront detruit quant le pere serat detruit
  virtual KN<K> *x() {return (**pere).x() ;}
  
  void operator=( KN<K> *y) { 
       **pere=y;
       set();}
  
};

*/
C_F0 NewFEarray(const char * id,Block *currentblock,C_F0 & fespacetype,CC_F0 init,bool cplx);
C_F0 NewFEarray(ListOfId * ids,Block *currentblock,C_F0 & fespacetype,CC_F0 init,bool cplx);
C_F0 NewFEvariable(const char * id,Block *currentblock,C_F0 & fespacetype,CC_F0 init,bool cplx);
C_F0 NewFEvariable(ListOfId * ids,Block *currentblock,C_F0 & fespacetype,CC_F0 init,bool cplx);
inline C_F0 NewFEvariable(const char * id,Block *currentblock,C_F0 & fespacetype,bool cplx)
  { CC_F0 init;init=0;return NewFEvariable( id,currentblock,fespacetype, init, cplx);}
    
inline C_F0 NewFEvariable(ListOfId * ids,Block *currentblock,C_F0 & fespacetype,bool cplx)
  { CC_F0 init;init=0;return NewFEvariable( ids,currentblock,fespacetype, init,cplx);}
 
 size_t dimFESpaceImage(const basicAC_F0 &args) ;

template<class K,class FE=FEbase<K> >
class E_FEcomp : public E_F0mps { public:
  typedef pair< FE * ,int> Result;
  Expression a0;
  const int comp, N;

  AnyType operator()(Stack s)  const {
    return SetAny<Result>( Result( *GetAny<FE **>((*a0)(s)),comp) );}  
  E_FEcomp(const C_F0 & x,const int cc,int NN) : a0(x.LeftValue()),comp(cc),N(NN)
      {if(x.left()!=atype<FE **>() ) 
        cout << "E_FEcomp: Bug " << *x.left() << " != " << *atype<FE **>() << "  case " <<typeid(K).name() << endl;
        //CompileError("E_FEcomp: Bug ?");
       throwassert(x.left()==atype<FE **>() &&a0);} 
    operator aType () const { return atype<Result>();}         
         
};

//typedef double R;
typedef  pair< FEbase<double> * ,int> aFEvarR;   
typedef  pair< FEbaseArray<double> * ,int> aFEArrayR;   
typedef  pair< FEbase<Complex> * ,int> aFEvarC;   

template<class K>
class interpolate_f_X_1 : public OneOperator  { public:
 //  to interpolate f o X^{-1} 
typedef FEbase<K> * pferbase ;
typedef FEbaseArray<K> * pferbasearray ;
typedef pair<pferbase,int> pfer ;
typedef pair<pferbasearray,int> pferarray ;
  class type {};
  class CODE : public E_F0mps { public:
    Expression f;
    Expression x,y;

    CODE(const basicAC_F0 & args) : f(args[0].LeftValue())  
      { const E_Array * X(dynamic_cast<const E_Array*>(args[1].LeftValue()));
        if ( ! X || X->size() != 2) 
           {
             CompileError("array of 2 double x,y   f([xx,yy]) = ... ");
           }
        x=to<double>((*X)[0]);
        y=to<double>((*X)[1]);
        assert(f && x && y ); 
      }
    AnyType operator()(Stack)  const  { ExecError("No evaluation"); return Nothing;} 
          operator aType () const { return atype<void>();}         

   }; //  end class CODE 
   
   E_F0 * code(const basicAC_F0 & args) const  
     { return  new CODE(args);}
     
   interpolate_f_X_1(): 
      OneOperator(map_type[typeid(type).name()],map_type[typeid(pfer).name()],map_type[typeid(E_Array).name()])
      {}
};

inline FESpace * v_fes::update() {     
    if (nbcperiodic ) {
       assert(periodic);
       const Mesh &Th(**ppTh);
       KN<int> ndfv(Th.nv);
       KN<int> ndfe(Th.neb);
       int nbdfv,nbdfe;    
       buildperiodic(stack,nbdfv,ndfv,nbdfe,ndfe);
       return   buildupdate(nbdfv,ndfv,nbdfe,ndfe);
      }
     else 
       return  buildupdate();
}

template<class A,class B>  A Build(B b) {  return A(b);}
