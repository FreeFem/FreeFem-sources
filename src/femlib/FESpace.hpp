#ifndef FESpace_h_
#define FESpace_h_
#include <cstdarg>

namespace  Fem2D {

class TypeOfFE;

// numbering of derivative 
enum operatortype { op_id=0,
   op_dx=1,op_dy=2,
   op_dxx=3,op_dyy=4,
   op_dyx=5,op_dxy=5,   
   op_dz=6,
   op_dzz=7,     
   op_dzx=8,op_dxz=8, 
   op_dzy=9,op_dyz=9   
   }; 
    
const int last_operatortype=10;


inline void initwhatd(bool *whatd,int k)
{
  for (int i=0;i<last_operatortype;i++)
    whatd[i]=false;    
  whatd[k]=true;
}

typedef KN_<R> RN_;
typedef KN<R>  RN;
typedef KNM_<R> RNM_;
typedef KNM<R>  RNM;
typedef KNMK_<R> RNMK_;
typedef KNMK<R>  RNMK;

inline int FindIn(int v,int *a,int n)
{
  for (int i=0;i<n;i++)
    if( v==a[i]) return i;
  return -1;
};


//  Methode boeuf 
//    on suppose que pour chaque elment fini
// on sait calculer 
//   un tableau de point x,y 
// les valeur de la fonction et des derives
// FB(RN_<R2> 
// i  =  df
// j  = [0 N[  (ou N est la dim de l'espace d'arrive N=3 
// k = 0,1,2   f,fx,fy 


class FElement;
class baseFElement;
class FMortar;
class TypeOfMortar;

//  for FE
//typedef void  (* basicFunction)(const FElement & FE,const R2 &P, RNMK_ & val);
typedef void  (* InterpolFunction)(R* v, const R2 & P,const baseFElement &,int where,const R2 & Phat, void * arg);
//typedef void (*InterpolantFunction)(const FElement & FE,RN_ & val, InterpolFunction f, R* v);
//  for FM   ( Finite Mortar  Mortar + interpolation)
typedef void  (* basicFunctionMortar)(const FMortar & FE,const R2 &P, RNMK_ & val);
typedef void  (* InterpolFunctionMortar)(R* v, const R2 & P,const Mortar &,int where,const R2 & Phat);
typedef void  (*InterpolantFunctionMortar)(const FMortar & FE,RN_ & val, InterpolFunctionMortar f, R* v);

//void P1Functions(const FElement &FE,const R2 & P,RNMK_ & val);
//void P2Functions(const FElement &FE,const R2 & P,RNMK_ & val);


class VofScalarFE { 
  R v[3];
  public:
  VofScalarFE() {v[0]=v[1]=v[2]=0;}
  VofScalarFE(R f,R fx,R fy) {v[0]=f;v[1]=fx;v[2]=fy;}

  R & operator[](int i){ return v[i];}
  const R & operator[](int i) const { return v[i];}  
  R f()  { return v[0];}
  R fx() { return v[1];}
  R fy() { return v[2];}
 VofScalarFE &operator+=(const VofScalarFE & a) { v[0]+=a.v[0]; v[1]+=a.v[1]; v[2]+=a.v[2];return *this;}  
};

class ConstructDataFElement {
  friend class FESpace;
  int thecounter;
  int * counter;
  int MaxNbNodePerElement;
  int MaxNbDFPerElement;
  int *NodesOfElement;
  int *FirstNodeOfElement;
  int *FirstDfOfNode;
  int NbOfElements;
  int NbOfDF;
  int NbOfNode;
  int Nproduit;
  ConstructDataFElement(const Mesh &Th,/*int NbDfOnSommet,int NbDfOnEdge,int NbDfOnElement,*/
  const  KN<const TypeOfFE *> & TFEs,const TypeOfMortar *tm=0,
  int nbdfv=0,const int *ndfv=0,int nbdfe=0,const int *ndfe=0);
  ConstructDataFElement(const ConstructDataFElement *,int k);
  ConstructDataFElement(const FESpace ** l,int k,const  KN<const TypeOfFE *> & TFEs) ;
  void renum(const long *r,int l)  ; 
  ~ConstructDataFElement();
 void Make(const Mesh &Th,/*int NbDfOnSommet,int NbDfOnEdge,int NbDfOnElement*/const  KN<const TypeOfFE *> & TFEs,const TypeOfMortar *tm=0,
   int nbdfv=0,const int *ndfv=0,int nbdfe=0,const int *ndfe=0);
  
};

template<class T>
inline int sum(const T ** l,int const T::*p,int n)
  {
    int r=0;
    for (int i=0;i<n;i++)
     r += l[i]->*p;
     return r;
  }
  
template<class T>
inline int max(const T ** l,int const T::*p,int n)
  {
    int r=0;
    for (int i=0;i<n;i++)
     r =Max( l[i]->*p,r);
     return r;
  }
 
 
class TypeOfFE { public:
//  The FEM is in R^N
//  The FEM is compose from nb_sub_fem
//  dim_which_sub_fem[N] give
  friend class FESpace;
  friend class FElement;
  friend class FEProduitConstruct;
  const int NbDoF;
  const int NbNodeOnVertex, NbNodeOnEdge, NbNodeOnElement;
  const int NbDfOnVertex, NbDfOnEdge, NbDfOnElement, N,nb_sub_fem;
  const  int NbNode;
  //  remark 
 // virtual void FB(const Mesh & Th,const Triangle & K,const R2 &P, RNMK_ & val) const =0;
  virtual void FB(const bool *,const Mesh & Th,const Triangle & K,const R2 &P, RNMK_ & val) const =0;
//  virtual void D2_FB(const Mesh & Th,const Triangle & K,const R2 &P, RNMK_ & val) const=0;
 // virtual void Pi_h(const baseFElement & K,RN_ & val, InterpolFunction f, R* v,int , void *) const=0;
  
  // soit 
  //  $(U_pj)_{{j=0,N-1}; {p=0,nbpoint_Pi_h-1}}$ les valeurs de U au point Phat[i];
  //   p est le numero du point d'integration
  //   j est la composante 
  //  l'interpole est  defini par
  //  Pi_h u = \sum_k u_k w^k , et u_i = \sum_pj alpha_ipj U_pj
  //  la matrice de alpha_ipj est tres creuse
  struct IPJ { 
    int i,p,j; // i is DoF, p is Point, j is componante 
    IPJ(int ii=0,int pp=0,int jj=0):i(ii),p(pp),j(jj) {}
  };
  
  friend KN<IPJ> Makepij_alpha(const TypeOfFE **,int );
  friend KN<R2> MakeP_Pi_h(const TypeOfFE **,int );
 
   const KN<IPJ > & Ph_ijalpha() const {return pij_alpha;} // la struct de la matrice
   const KN<R2> & Pi_h_R2() const { return P_Pi_h;}   // les points
   virtual void Pi_h_alpha(const baseFElement & K,KN_<double> & v) const
     {  assert(coef_Pi_h_alpha);
        v=KN_<double>(coef_Pi_h_alpha,pij_alpha.N());}

   //  ----  
  
  const int nbsubdivision; // nb of subdivision for plot
   int const * const DFOnWhat; //  0,1,2 vertex 3,4,5 edge 6  triangle
   int const * const DFOfNode; // n\u00c9 du df on Node
   int const * const NodeOfDF; // 
   int const * const fromFE;   //  the df  come from df of FE
   int const * const fromDF;   //  the df  come from with FE
   int const * const dim_which_sub_fem;
   KN<IPJ > pij_alpha ;
   KN<R2 > P_Pi_h ;
  double *coef_Pi_h_alpha;
   
  // if 0  no plot 
  public: 
  
    TypeOfFE(const TypeOfFE & t,int k,const int * data) 
    : 
      NbDoF(t.NbDoF*k),
      NbNodeOnVertex(t.NbNodeOnVertex),NbNodeOnEdge(t.NbNodeOnEdge),NbNodeOnElement(t.NbNodeOnElement),
      NbDfOnVertex(t.NbDfOnVertex*k),NbDfOnEdge(t.NbDfOnEdge*k),NbDfOnElement(t.NbDfOnElement*k),
      N(t.N*k),nb_sub_fem(t.nb_sub_fem*k),
      NbNode(t.NbNode), 
     nbsubdivision(t.nbsubdivision),
    DFOnWhat(data),
    DFOfNode(data+NbDoF),
    NodeOfDF(data+2*NbDoF),
    fromFE(data+3*NbDoF),
    fromDF(data+4*NbDoF),
    dim_which_sub_fem(data+5*NbDoF),
    pij_alpha(t.pij_alpha.N()*k),P_Pi_h(t.P_Pi_h),
    coef_Pi_h_alpha(0)
    
    { throwassert(dim_which_sub_fem[N-1]>=0 && dim_which_sub_fem[N-1]< nb_sub_fem);
     // Warning the componant is moving first 
        for (int j=0,l=0;j<t.pij_alpha.N();j++) // for all sub DF
          for(int i=0,i0=0; i<k; i++,l++) // for componate
          {
          pij_alpha[l].i=t.pij_alpha[j].i*k+i;   //  DoF number
          pij_alpha[l].p=t.pij_alpha[j].p;       //  point of interpolation
          pij_alpha[l].j=t.pij_alpha[j].j+i*t.N; //  componante of interpolation
          }                         
    } 
       
    TypeOfFE(const TypeOfFE ** t,int k,const int * data) 
    : 
      NbDoF(sum(t,&TypeOfFE::NbDoF,k)),
      NbNodeOnVertex(NbNodebyWhat(data,NbDoF,0)),
      NbNodeOnEdge(NbNodebyWhat(data,NbDoF,3)),
      NbNodeOnElement(NbNodebyWhat(data,NbDoF,6)),
      NbDfOnVertex(sum(t,&TypeOfFE::NbDfOnVertex,k)),
      NbDfOnEdge(sum(t,&TypeOfFE::NbDfOnEdge,k)),
      NbDfOnElement(sum(t,&TypeOfFE::NbDfOnElement,k)),
      N(sum(t,&TypeOfFE::N,k)),nb_sub_fem(sum(t,&TypeOfFE::nb_sub_fem,k)),
      NbNode( (NbDfOnVertex ? 3 :0) + (NbDfOnEdge ? 3 :0 ) +(NbDfOnElement? 1 :0) ),      
      nbsubdivision(max(t,&TypeOfFE::nbsubdivision,k)),
      
    DFOnWhat(data),
    DFOfNode(data+NbDoF),
    NodeOfDF(data+2*NbDoF),
    fromFE(data+3*NbDoF),
    fromDF(data+4*NbDoF),
    dim_which_sub_fem(data+5*NbDoF),
    pij_alpha(Makepij_alpha(t,k)),
    P_Pi_h(MakeP_Pi_h(t,k)),
    coef_Pi_h_alpha(0)

    
   { throwassert(dim_which_sub_fem[N-1]>=0 && dim_which_sub_fem[N-1]< nb_sub_fem);} 

  TypeOfFE(const int i,const int j,const int k,const int NN,const  int  *   data,int nsub,int nbsubfem,
    int kPi,int npPi,double * coef_Pi_h_a=0) 
   : NbDoF(3*(i+j)+k),
     NbNodeOnVertex(NbNodebyWhat(data,NbDoF,0)),
     NbNodeOnEdge(NbNodebyWhat(data,NbDoF,3)),
     NbNodeOnElement(NbNodebyWhat(data,NbDoF,6)),     
/*     NbDfOnVertex(Count(data,NbDoF,0)),
     NbDfOnEdge(Count(data,NbDoF,3)),
     NbDfOnElement(Count(data,NbDoF,6)),
 */  
     NbDfOnVertex(i),NbDfOnEdge(j),NbDfOnElement(k),N(NN),nb_sub_fem(nbsubfem),
     NbNode( (NbDfOnVertex ? 3 :0) + (NbDfOnEdge ? 3 :0 ) +(NbDfOnElement? 1 :0) ),
    nbsubdivision(nsub),
    DFOnWhat(data),
    DFOfNode(data+NbDoF),
    NodeOfDF(data+2*NbDoF),
    fromFE(data+3*NbDoF),
    fromDF(data+4*NbDoF),
    dim_which_sub_fem(data+5*NbDoF),
    pij_alpha(kPi),P_Pi_h(npPi),
    coef_Pi_h_alpha(coef_Pi_h_a)
    
     { 
     // cout << "TypeOfFE " <<NbDoF << " : " << NbDfOnVertex << " " << NbDfOnEdge << " " << NbDfOnElement << 
     // " : " << NbNodeOnVertex << " " << NbNodeOnEdge << " " << NbNodeOnElement << endl;
      assert(NbDfOnVertex==Count(data,NbDoF,0));
      assert(NbDfOnVertex==Count(data,NbDoF,1));
      assert(NbDfOnVertex==Count(data,NbDoF,2));
      assert(NbDfOnEdge==Count(data,NbDoF,3));
      assert(NbDfOnEdge==Count(data,NbDoF,4));
      assert(NbDfOnEdge==Count(data,NbDoF,5));
      assert(NbDfOnElement==Count(data,NbDoF,6));
      
      throwassert(dim_which_sub_fem[N-1]>=0 && dim_which_sub_fem[N-1]< nb_sub_fem);} 
      
  TypeOfFE(const int nbdf,const int NN,const  int  *   data,int nsub,int nbsubfem,
    int kPi,int npPi,double * coef_Pi_h_a=0) 
   : 
     NbDoF(nbdf),
     NbNodeOnVertex(NbNodebyWhat(data,NbDoF,0)),
     NbNodeOnEdge(NbNodebyWhat(data,NbDoF,3)),
     NbNodeOnElement(NbNodebyWhat(data,NbDoF,6)),     
     NbDfOnVertex(Count(data,NbDoF,0)),
     NbDfOnEdge(Count(data,NbDoF,3)),
     NbDfOnElement(Count(data,NbDoF,6)),
     N(NN),
     nb_sub_fem(nbsubfem),
     NbNode( (NbDfOnVertex ? 3 :0) + (NbDfOnEdge ? 3 :0 ) +(NbDfOnElement? 1 :0) ),
    nbsubdivision(nsub),
    DFOnWhat(data),
    DFOfNode(data+NbDoF),
    NodeOfDF(data+2*NbDoF),
    fromFE(data+3*NbDoF),
    fromDF(data+4*NbDoF),
    dim_which_sub_fem(data+5*NbDoF),
    pij_alpha(kPi),P_Pi_h(npPi),
    coef_Pi_h_alpha(coef_Pi_h_a)
    
     { 
      assert(NbDfOnVertex==Count(data,NbDoF,0));
      assert(NbDfOnVertex==Count(data,NbDoF,1));
      assert(NbDfOnVertex==Count(data,NbDoF,2));
      assert(NbDfOnEdge==Count(data,NbDoF,3));
      assert(NbDfOnEdge==Count(data,NbDoF,4));
      assert(NbDfOnEdge==Count(data,NbDoF,5));
      assert(NbDfOnElement==Count(data,NbDoF,6));
      
      throwassert(dim_which_sub_fem[N-1]>=0 && dim_which_sub_fem[N-1]< nb_sub_fem);} 

  virtual ~TypeOfFE() { }
  virtual R operator()(const FElement & K,const  R2 & PHat,const KN_<R> & u,int componante,int op)
 const ;
 
 private:
 static int Count(const int *data,int n,int which) 
   {
      int kk=0;
      for (int i=0;i<n;++i)
        if (which == data[i]) ++kk;
      return kk;}
      
 static int LastNode(const int *data,int n) 
   {
      int kk=0,i0=2*n;
      for(int i=0;i<n;i++)
        kk=Max(kk,data[i+i0]);
      return kk;}
      
static int NbNodebyWhat(const int *data,int n,int on) 
  { 
      const int nmax=100;
      int t[nmax];
      for (int i=0;i<n;i++)
      	t[i]=0;
       int k=0,i0=2*n;
       for(int i=0;i<n;i++)
          if ( data[i]==on)
           { int node= data[i+i0];
           //  cout << " node " << node << endl;
             assert(node < nmax);
            if( ! t[node] )
               t[node] = 1,++k;
            } 
           
//     cout << " on " << on << " k = " << k << endl;
     return k;
  }      
      
} ; 



class aSubFMortar;

class TypeOfMortar { 
  friend class FESpace;
  friend class FMortar;
  friend class ConstructDataFElement;
  friend class TheSubFMortars;
  virtual int NbLagrangeMult(const Mesh &Th,const Mortar &M) const =0;
  virtual int NbDoF(const Mesh &Th,const Mortar &M,int i) const =0;
  virtual int NbOfNodes(const Mesh &Th,const Mortar &M) const =0;
  virtual int NbDoF(const Mesh &Th,const Mortar &M) const =0;
  virtual int NodeOfDF(const FESpace &,const Mortar &M,int i) const =0;
  virtual int DFOfNode(const FESpace &,const Mortar &M,int i) const =0;
  virtual void ConstructionOfNode(const Mesh &Th,int im,int * NodesOfElement,int *FirstNodeOfElement,int &lastnodenumber) const=0; 
  virtual void ConsTheSubMortar(FMortar &) const =0; 
  protected:
  const int NbDfOnVertex, NbDfOnEdge , N;
  public: 
    TypeOfMortar (int i,int j): NbDfOnVertex(i),NbDfOnEdge(j),N(1){};
     
} ;

class FElement; 

class baseFElement { public:
  const FESpace &Vh;  
  const Triangle &T;
  const TypeOfFE * tfe; 
  const int N;
  const int number;  
  baseFElement(  const FESpace &aVh, int k) ;
  baseFElement(const   baseFElement & K,  const TypeOfFE & atfe) ;
  R EdgeOrientation(int i) const {return T.EdgeOrientation(i);}
  //  : Vh(aVh),T(Vh.Th[k]),tfe(aVh.TFE[k]),N(aVh.N),number(k){}

};

class FElement : public baseFElement { public:
  typedef const KN<TypeOfFE::IPJ> &  aIPJ;
  typedef  TypeOfFE::IPJ   IPJ;
  typedef const KN<R2> &  aR2;
  
  
   friend class FESpace;
  const int *p;
  const int nb;
  FElement(const FESpace * VVh,int k) ;
  public:
  int NbOfNodes()const {return nb;}
  int  operator[](int i) const ;//{ return  p ? p+i :  ((&T[i])-Vh.Th.vertices);}  N\u00c9 du noeud 
  int  NbDoF(int i) const ;  // number of DF 
  int  operator()(int i,int df) const ;// { N\u00c9 du DoF du noeud i de df local df 
  int  operator()(int df) const { return operator()(NodeOfDF(df),DFOfNode(df));}
  void Draw(const KN_<R> & U, const  KN_<R> & VIso,int j=0) const ;  
  void Drawfill(const KN_<R> & U, const  KN_<R> & VIso,int j=0,double rapz=1) const ;  
  void Draw(const RN_& U,const RN_& V, const  KN_<R> & Viso,R coef,int i0,int i1) const;
  R2   MinMax(const RN_& U,const RN_& V,int i0,int i1) const  ;
  R2   MinMax(const RN_& U,int i0) const  ;
  void BF(const R2 & P,RNMK_ & val) const;// { tfe->FB(Vh.Th,T,P,val);}
  void BF(const bool * whatd, const R2 & P,RNMK_ & val) const;// { tfe->FB(Vh.Th,T,P,val);}
 // void D2_BF(const R2 & P,RNMK_ & val) const ; //{ tfe->D2_FB(Vh.Th,T,P,val);}
  void Pi_h(RN_ val,InterpolFunction f,R *v,  void * arg) const;//  {tfe->Pi_h(Vh.Th,T,val,f,v);}
  aIPJ Pi_h_ipj() const { return tfe->Ph_ijalpha();}
  aR2 Pi_h_R2() const { return tfe->Pi_h_R2();}
  void Pi_h(KN_<R> & v) const { return tfe->Pi_h_alpha(*this,v);}
  
  int NbDoF() const { return tfe->NbDoF;}
  int DFOnWhat(int i) const { return tfe->DFOnWhat[i];}
  int FromDF(int i) const { return tfe->fromDF[i];}
  int FromFE(int i) const { return tfe->fromFE[i];}
 // df is the df in element 
  int NodeOfDF(int df) const { return tfe->NodeOfDF[df];} // a node
  int DFOfNode(int df) const { return tfe->DFOfNode[df];} // the df number on the node 
 
  R operator()(const R2 & PHat,const KN_<R> & u,int i,int op)  const ;
 // FElementGlobalToLocal operator()(const KN_<R> & u ) const { return FElementGlobalToLocal(*this,u);}
  private:
  int nbsubdivision() const { return tfe->nbsubdivision;} // for draw 
  };
  
  
class aSubFMortar { public:
  R a,b; 
  R2 A,B; 
  int left,right; //
  int TLeft() const { return left/3;}
  int ELeft() const{ return left%3;}
  
  int TRight() const{ return right/3;}
  int ERight()const{ return right%3;}
  
  int  Nbmul; 
  int * DfNumberOFmul;
 // int *whatnode[2];
  R (**f)(const FESpace *,const aSubFMortar * ,R); // array of function of mul on aSubMortar 
  R lg1(){throwassert(a<b);return b-a;}  
  aSubFMortar():f(0),DfNumberOFmul(0),Nbmul(0){}
 
  R2 operator()(const Mesh & Th,R x,int side) const
   {
      int kk = side==0 ? left : right;
      int k=kk/3;
      int e=kk%3;
      const Triangle & K(Th[k]);
      int c = VerticesOfTriangularEdge[e][0];
      int d = VerticesOfTriangularEdge[e][1];
      
      const R2 &C=K[c];
      const R2 &D=K[d];
      // $ [AB] \include [CD] $
      R2 CD(C,D);
      R  CD2=(CD,CD);
      R la= (CD,A-C)/CD2;
      R lb= (CD,B-C)/CD2;
      //  A = C*(1-la)+ la*D
      //  B = C*(1-lb)+ lb*D
      //  X = (1-x)*A+x*B;
      //  X = (1-x)*(C*(1-la)+ la*D) + x*(C*(1-lb)+ lb*D)
      //  X = C*((1-x)*(1-la)+x*(1-lb))
      R xx=((1-x)*(1-la)+x*(1-lb));
      throwassert(xx>=0 && xx <= 1);
      return TriangleHat[c]*(xx)+TriangleHat[d]*(1-xx);
   }
  
};  



class FMortar { public:
  friend class FESpace;
  const FESpace &Vh;  
  const Mortar &M;
  const int *p;
  const int nbn; //  nb of node
//  const int nbnl,nbnr; // nb of node left, right
  
  //   we supposse 
  //  the node numbering in a mortar is
  //  1   the lagrange mul.
  //  2   the node of left side
  //  3   the node of right side
  
  const int N;  
   TypeOfMortar const *  const tom; 
  FMortar(const FESpace * VVh,int k) ;
  public:
  int  operator[](int i) const ;
  int  operator()(int i,int df) const ;// { N\u00c9 du DoF du noeud i de df local df 
  int  operator()(int df) const { return operator()(NodeOfDF(df),DFOfNode(df));}
  int NbDoF(int i) const ;//{return tom->NbDoF(Vh.Th,M,i);};  // number of DF 
  int NbOfNodes()const {return nbn;}
  int NbDoF() const ;
  int NodeOfDF(int i) const ;
  int DFOfNode(int i) const ;
  
  int nbsm; // nb of submortar
  aSubFMortar * sm;
  ~FMortar() { 
    delete [] datai;
    delete [] dataf;}
  
  int *datai;
   R (**dataf)(const FESpace *,const aSubFMortar *,R);
  private:
   FMortar(const FMortar &);
   void operator=(const FMortar &);

};

extern TypeOfFE & P1Lagrange;
extern TypeOfFE & P2Lagrange;
extern TypeOfFE & RTLagrange;
extern TypeOfFE & RTLagrangeOrtho;
extern TypeOfFE & P0Lagrange;
extern TypeOfFE & P1ncLagrange;

class FESpace : public RefCounter {
  public:
  TypeOfFE const * const ptrTFE;
  KN<const TypeOfFE *>  TFE; 
  private:
  ConstructDataFElement * cdef; //  juste pour les constantes 
  public:
  CountPointer<const Mesh> cmesh;
  const int N; // dim espace d'arrive
  const int Nproduit; // dim de l'espace produit generalement 1
  const int NbOfDF;
  const int NbOfElements;
  const int NbOfNodes;
  const Mesh &Th;
  const int nb_sub_fem; // nb de sous elements finis tensorise (independe au niveau des composantes)
  int const* const dim_which_sub_fem;// donne les dependant des composantes liee a un meme sous element fini
   //   exemple si N=5,  
   // dim_which_sub_fem[0]=0;
   // dim_which_sub_fem[1] =1;
   // dim_which_sub_fem[2]= 2
   // dim_which_sub_fem[3] =2 
   // dim_which_sub_fem[4] =3
   //  =>
   //  le  sous  elements fini 0 est lie a la composante 0 
   //  le  sous  elements fini 1 est lie a la composante 1 
   //  le  sous  elements fini 2 est lie aux composantes 2,3 
   //  le  sous  elements fini 3 est lie a la composante 4 
   //  donc pour les CL. les composante 2 et 3 sont lie car elle sont utiliser pour definir un
   //  meme degre de libert\u00e9.
  
  int const*const NodesOfElement;
  int const*const FirstNodeOfElement;
  int const*const FirstDfOfNodeData;
  const TypeOfMortar * tom; 
  const int MaxNbNodePerElement;
  const int MaxNbDFPerElement;
  int FirstDFOfNode(int i) const {return FirstDfOfNodeData ? FirstDfOfNodeData[i] : i*Nproduit;}
  int LastDFOfNode(int i)  const {return FirstDfOfNodeData ? FirstDfOfNodeData[i+1] : (i+1)*Nproduit;}
  int NbDFOfNode(int i)  const {return FirstDfOfNodeData ? FirstDfOfNodeData[i+1]-FirstDfOfNodeData[i] : Nproduit;}
  int MaximalNbOfNodes() const {return MaxNbNodePerElement;};
  int MaximalNbOfDF() const {return MaxNbDFPerElement;};
  const int * PtrFirstNodeOfElement(int k) const {
      return NodesOfElement 
               ? NodesOfElement + (FirstNodeOfElement ? FirstNodeOfElement[k] : k*MaxNbNodePerElement)                 
               : 0;}   
               
  int SizeToStoreAllNodeofElement() const {  
       return  FirstNodeOfElement 
                ?  FirstNodeOfElement[NbOfElements] 
                : MaxNbNodePerElement*NbOfElements;}   
                      
  int NbOfNodesInElement(int k)   const {             
      return FirstNodeOfElement 
               ?  FirstNodeOfElement[k+1] - FirstNodeOfElement[k] 
               : MaxNbNodePerElement ;}
  int  esize() const { return  MaxNbDFPerElement*N*last_operatortype;}   // size to store all value of B. function        
//  par defaut P1                
  FESpace(const Mesh & TTh) 
   :
      ptrTFE(0),
      TFE(1,0,&P1Lagrange),tom(0),
      cmesh(TTh),cdef(0),N(1),Nproduit(1),Th(TTh),NbOfDF(TTh.nv),NbOfElements(TTh.nt),NbOfNodes(TTh.nv),
      FirstDfOfNodeData(0),
      MaxNbNodePerElement(3),MaxNbDFPerElement(3*Nproduit),
      NodesOfElement(0),FirstNodeOfElement(0),
      nb_sub_fem(TFE[0]->nb_sub_fem),
      dim_which_sub_fem(TFE[0]->dim_which_sub_fem)
 {}
      
    FESpace(const FESpace &,int k );
    FESpace(const FESpace **,int k ); 
    FESpace(const Mesh & TTh,const TypeOfFE **,int k,int nbdfv=0,const int *ndfv=0,int nbdfe=0,const int *ndfe=0 );//int NbDfOnSommet,int NbDfOnEdge,int NbDfOnElement,int NN=1); 
 
    FESpace(const Mesh & TTh,const TypeOfFE & ,
    int nbdfv=0,const int *ndfv=0,int nbdfe=0,const int *ndfe=0);//int NbDfOnSommet,int NbDfOnEdge,int NbDfOnElement,int NN=1); 
    
    FESpace(const Mesh & TTh,const TypeOfFE &,const TypeOfMortar & );//int NbDfOnSommet,int NbDfOnEdge,int NbDfOnElement,int NN=1); 
  ~FESpace();    
 // FESpace(Mesh & TTh,int NbDfOnSommet,int NbDfOnEdge,int NbDfOnElement,int NN=1); 
  int  renum();
      
  FElement operator[](int k) const { return FElement(this,k);} 
  FElement operator[](const Triangle & K) const { return FElement(this,Th.number(K));} 
  int  operator()(int k)const {return NbOfNodesInElement(k);}
  int  operator()(int k,int i) const { //  the node i of element k
     return NodesOfElement ?  *(PtrFirstNodeOfElement(k) + i)  : Th(k,i)  ;}
  
  void Draw(const KN_<R>& U,const KN_<R>& Viso,int j=0) const ; // Draw iso line
  void Drawfill(const KN_<R>& U,const KN_<R>& Viso,int j=0,double rapz=1) const ; // Draw iso line
 
  void Draw(const KN_<R>& U,const KN_<R> & Viso, R coef,int j0=0,int j1=1) const  ; // Arrow
  void Draw(const KN_<R>& U,const KN_<R>& V,const KN_<R> & Viso, R coef,int iu=0,int iv=0) const  ; // Arrow
  R2 MinMax(const KN_<R>& U,const KN_<R>& V,int j0,int j1,bool bb=true) const ;
  R2 MinMax(const KN_<R>& U,int j0, bool bb=true) const ;
 // void destroy() {RefCounter::destroy();}
  bool isFEMesh() const { return !cdef && ( N==1) ;} // to make optim
 private: // for gibbs  
  int gibbsv (long* ptvoi,long* vois,long* lvois,long* w,long* v);
};

inline baseFElement::baseFElement(  const FESpace &aVh, int k) 
    : Vh(aVh),T(Vh.Th[k]),tfe(aVh.TFE[k]),N(aVh.N),number(k){}
    
inline baseFElement::baseFElement(const   baseFElement & K,  const TypeOfFE & atfe) 
    : Vh(K.Vh),T(K.T),tfe(&atfe),N(Vh.N),number(K.number){}


inline FElement::FElement(const FESpace * VVh,int k) 
  : baseFElement(*VVh,k) ,
    p(Vh.PtrFirstNodeOfElement(k)),
    nb(Vh.NbOfNodesInElement(k))
    
     {} 
  
inline   int  FElement::operator[](int i) const {
   return  p ? p[i] :  ((&T[i])-Vh.Th.vertices);}  
   
inline   int  FElement::operator()(int i,int df) const {
   return  Vh.FirstDFOfNode(p ? p[i] :  ((&T[i])-Vh.Th.vertices)) + df;}  
   
inline   int  FMortar::operator()(int i,int df) const {throwassert(p);
   return  Vh.FirstDFOfNode(p[i]) + df;} 
    
inline   int  FMortar::operator[](int i) const {throwassert(p);
   return  p[i];}  
   
inline   int  FElement::NbDoF(int i) const {
   int node =p ? p[i] :  ((&T[i])-Vh.Th.vertices);
   return  Vh.LastDFOfNode(node)-Vh.FirstDFOfNode(node);}  

void  SetDefaultIsoValue(const KN_<R>& U,KN_<R> & Viso);
void  SetDefaultIsoValue(const KN_<R>& u,const KN_<R>& v,KN_<R> & Viso);
void MoveTo(R2 P); 
void LineTo(R2 P) ;

/*
void operator=(  KN_<R> & u,const FElementGlobalToLocal & x) 
{
  int n=u.N();
  throwassert(n==x.S.NbDoF());
  for (int i=0;i<n;i++) // get the local value
     v[i] = x.U[x.S(i)];
}
*/
inline  int FMortar::NbDoF(int i) const {
   int node = p[i];
   return  Vh.LastDFOfNode(node)-Vh.FirstDFOfNode(node);
};  // number of DF 
inline  int FMortar::NbDoF() const { return tom->NbDoF(Vh.Th,M);}
//inline  int FMortar::NbOfNodes()const {return }
inline  int FMortar::NodeOfDF(int i) const { return tom->NodeOfDF(Vh,M,i);}
inline  int FMortar::DFOfNode(int i) const { return tom->DFOfNode(Vh,M,i);}

inline ostream & operator << (ostream & f,const FElement & FE)
   {
     f << FE.number << "," <<FE.nb << ":" ;
     for (int i=0;i<FE.nb;i++) f << "\t"<<FE.p[i];
   return f;
     
   }
inline ostream & operator << (ostream & f,const FESpace & Vh)
   {
     cout << " list of nodes per elements :" << endl;
     for (int k=0;k< Vh.NbOfElements;k++)
      { f <<setw(3) <<  k << ":";
      for (int j=0;j<Vh(k);j++)
        f << " " << setw(3) << Vh(k,j);
        cout << endl;
      }
      
      f << endl << " FirstDFOfNode :" ;       
      for (int i=0;i<=Vh.NbOfNodes;i++)
        {if (i%10==0)  cout << "\n" << setw(3) << i << " : ";
        cout << setw(3) << Vh.FirstDFOfNode(i) << " ";}
    
   return f;
     
   }
   
inline void FElement::BF(const R2 & P,RNMK_ & val) const {
 static bool whatdold[last_operatortype]={true,true,true,false,false,false,false,false,false,false};
 tfe->FB(whatdold,Vh.Th,T,P,val);}
inline void FElement::BF(const bool * whatd,const R2 & P,RNMK_ & val) const { tfe->FB(whatd,Vh.Th,T,P,val);}
//inline  void FElement::D2_BF(const R2 & P,RNMK_ & val) const { tfe->D2_FB(Vh.Th,T,P,val);}


//  -------
 extern const TypeOfMortar & TheMortarCas1P2; 
 
void PlotValue(const RN_ & Viso,int  k0,const char * cmm);
 
// to store all the type of TFE
// the problem is the TFE can be define on lot of file.cpp
struct ListOfTFE { 
  const char * name;
  TypeOfFE * tfe;
  ListOfTFE * next;

  static ListOfTFE * all ; // list of all object of this type 
  ListOfTFE (const char * n,TypeOfFE *t);
};
// to get a unique list of TypeOfFE 
// local variable of TypeOfFE
ListOfTFE & GetListOfTFE() ;


inline   R FElement::operator()(const R2 & PHat,
                                const KN_<R> & u,int i,int op)  const
{
 return (*tfe)(*this,PHat,u,i,op);
}


}


#endif
