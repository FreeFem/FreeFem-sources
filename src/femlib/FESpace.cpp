#include <cmath>
#include <cstdlib>
#include "error.hpp"
#include <iostream>
#include <fstream>
#include <map>
#include "rgraph.hpp"
using namespace std;  

#include "RNM.hpp"
#include "fem.hpp"
#include "FESpace.hpp"

namespace  Fem2D {

 int  Make(const TypeOfFE ** t,int k,KN<R2> & P,KN<int> & I)
 {
   typedef  TypeOfFE::IPJ IPJ;

   int n=0,nn=0;

   for (int i=0;i<k;i++)
    {
     const KN<R2> p(t[i]->P_Pi_h);
     for (int j=0;j<p.N();j++,nn++)
     {
       P[n]=p[j]; // par defaut un nouveau => ajout
       I[nn]=n++;       
       for (int l=0;l<n-1;l++)
        {
          R2 QP(p[j],P[l]);
          if ( (QP,QP) < 1.0e-12 ) { 
            I[nn]=l;
            n--; // on a retrouver =>  detruit
            break;}
        }
       
     }
    }
   return n; // nombre de point commun
 }
 
 KN<TypeOfFE::IPJ > Makepij_alpha(const TypeOfFE ** t,int k)
 {
 // Attention les df est numerote de facon croissant 
 // en faisant une boucle sur les TypeOfFE
 // comme dans la class TypeOfFESum
  typedef TypeOfFE::IPJ IPJ;
  int n=0,m=0;
  for (int i=0;i< k;i++) {
    n += t[i]->pij_alpha.N();
    m += t[i]->P_Pi_h.N();
    }
  KN<TypeOfFE::IPJ > ij(n);
  KN<int> I(m);
  KN<R2> P(m);
  Make(t,k,P,I);
  int p0=0,i0=0,N0=0,nn=0;
  for (int i=0;i< k;i++) {
    const KN<IPJ > p(t[i]->pij_alpha);
     for (int j=0;j<p.N();j++,nn++)  
      {
        ij[nn].i= p[j].i+i0; //  comme dans TypeOfFESum
        ij[nn].j= N0+ p[j].j;
        ij[nn].p= I[p[j].p+p0];
       
      //   cout << nn << "Makepij_alpha: " << ij[nn].i << " " << ij[nn].p << " " << ij[nn].j << endl;
        
      }
    i0+=t[i]->NbDoF;
    p0+=t[i]->P_Pi_h.N();
    N0+=t[i]->N;}
  return ij; 
 }
 KN<R2 > MakeP_Pi_h(const TypeOfFE **t,int k)
 {
  int np=0;
  for (int i=0;i< k;i++) 
    np += t[i]->P_Pi_h.N();
    
  KN< R2 >  yy(np);
  KN<int> zz(np);
  int kk=Make(t,k,yy,zz);
 // cout << " MakeP_Pi_h: " << kk << " from " << np << endl; 
  return yy(SubArray(kk));
 
 }

ListOfTFE * ListOfTFE::all ; // list of all object of this type 

ListOfTFE::ListOfTFE (const char * n,TypeOfFE *t) : name(n),tfe(t) 
{
  if(!t)
  assert(t);
  static int count=0;
  if (count++==0) 
    all=0; // init of all in dependant of the ordre of the objet file   
  next=all;
  all=this;
}

const TypeOfFE ** Make(const FESpace **l,int k) {
  const TypeOfFE** p=new const TypeOfFE*[k];
  for (int i=0;i<k;i++)
    p[i]=l[i]->TFE[0];
  return p;
}
const TypeOfFE ** Make(const TypeOfFE **l,int k) {
  const TypeOfFE** p=new const TypeOfFE*[k];
  for (int i=0;i<k;i++)
    p[i]=l[i];
  return p;
}

  
bool Same(const FESpace **l,int k)
{
   for (int i=1;i<k;i++)
    if (l[0] != l[i] ) return false;
   return true;
}   

class FESumConstruct { protected:
   const TypeOfFE ** teb;
   const int k;
   int nbn; // nb of node 
   int  *  data;
   int  * const NN; //  NN[ i:i+1[ dimension de l'element i
   int  * const DF; // DF[i:i+1[  df associe a l'element i
   int  * const comp; //  
   FESumConstruct(int kk,const TypeOfFE **t);   
   virtual ~FESumConstruct(){
     delete [] DF;
     delete [] NN;
     delete [] comp;
     delete [] data;}   
};

  void FElement::Pi_h(RN_ val,InterpolFunction f,R *v, void * arg=0) const {
   // routine: a  tester FH.
    FElement::aIPJ ipj(Pi_h_ipj()); 
    FElement::aR2  PtHat(Pi_h_R2()); 
    KN<R>   Aipj(ipj.N());
    KNM<R>  Vp(N,PtHat.N());
        
     Pi_h(Aipj);
     for (int p=0;p<PtHat.N();p++)
          { 
            f(v,T(PtHat[p]),*this,T.lab,PtHat[p],arg);
            KN_<double> Vpp(Vp('.',p));
            for (int j=0;j<N;j++)          
               Vpp[j]=v[j];              
           }
           
         for (int i=0;i<Aipj.N();i++)
          { 
           const FElement::IPJ &ipj_i(ipj[i]);
           val[ipj_i.i] += Aipj[i]*Vp(ipj_i.j,ipj_i.p);           
          }
 }

class TypeOfFESum: public FESumConstruct, public  TypeOfFE { public:
   TypeOfFESum(const FESpace **t,int kk): 
     FESumConstruct(kk,Make(t,kk)),TypeOfFE(teb,kk,data) {}
       TypeOfFESum(const TypeOfFE **t,int kk): 
     FESumConstruct(kk,Make(t,kk)),TypeOfFE(teb,kk,data) {}

  // void FB(const Mesh & Th,const Triangle & K,const R2 &P, RNMK_ & val) const;
   void FB(const bool * whatd,const Mesh & Th,const Triangle & K,const R2 &P, RNMK_ & val) const;
//   void D2_FB(const Mesh & Th,const Triangle & K,const R2 &P, RNMK_ & val) const;
//  void Pi_h(const baseFElement & K,RN_ & val, InterpolFunction f, R* v,int, void * arg ) const; 
   virtual void Pi_h_alpha(const baseFElement & K,KN_<double> & v) const
    { 
      for (int i=0,k0=0;i<k;i++) {
          int n=teb[i]->NbDoF;
          KN_<R> sv(v(SubArray(n,k0)));
          teb[i]->Pi_h_alpha(K,sv);
          k0+= n;}
    }
   ~TypeOfFESum(){  delete []  teb;}
} ;

class FEProduitConstruct { protected:
   const TypeOfFE & teb;
   int k;
   int * data;
   FEProduitConstruct(int kk,const TypeOfFE &t)  ;   
   ~FEProduitConstruct(){delete [] data;}   
};

class TypeOfFEProduit: protected FEProduitConstruct, public  TypeOfFE { public:  
   TypeOfFEProduit(int kk,const TypeOfFE &t): 
     FEProduitConstruct(kk,t),TypeOfFE(t,kk,data)  {}
     
  // void FB(const Mesh & Th,const Triangle & K,const R2 &P, RNMK_ & val) const;
   void FB(const bool * whatd,const Mesh & Th,const Triangle & K,const R2 &P, RNMK_ & val) const;
//   void D2_FB(const Mesh & Th,const Triangle & K,const R2 &P, RNMK_ & val) const;
//   void Pi_h(const baseFElement & K,RN_ & val, InterpolFunction f, R* v,int, void * arg ) const;  
   virtual void Pi_h_alpha(const baseFElement & K,KN_<double> & v) const
    { int nbof=teb.NbDoF;
      for (int i=0,k0=0;i<k;i++,k0+=nbof)
        {          
          KN_<R> sv(v(SubArray(nbof,k0)));
          teb.Pi_h_alpha(K,sv);
        }
    }
 
   ~TypeOfFEProduit(){}
} ;

FEProduitConstruct::FEProduitConstruct(int kk,const TypeOfFE &t)
 :k(kk),teb(t) 
{  
  int m= teb.NbDoF;
  KN<int> nn(teb.NbNode);
  nn=0; // nb de dl par noeud 
  for (int i=0;i<m;i++)
    nn[teb.NodeOfDF[i]]++;   
    
  int n= m*kk;
  int N= teb.N*kk;
  data = new int [n*5+N];
  int c=0;
  
  for (int i=0;i<m;i++)
   for (int j=0;j<kk;j++)
    data[c++] = teb.DFOnWhat[i];
  
  for (int i=0;i<m;i++)
   for (int j=0;j<kk;j++) // num of df on node  for the df = j 
      data[c++] = teb.DFOfNode[i]+j*nn[teb.NodeOfDF[i]];
    
  for (int i=0;i<m;i++)
   for (int j=0;j<kk;j++) 
    data[c++] = teb.NodeOfDF[i]; //  node of df

  for (int i=0;i<m;i++)
   for (int j=0;j<kk;j++) 
    data[c++] = j; //  node from of FE
    
  for (int i=0;i<m;i++)
   for (int j=0;j<kk;j++) 
    data[c++] = i; //  node from of df in FE
    
   for (int j=0;j<kk;j++) 
    for (int i=0;i<teb.N;i++) 
      data[c++]= teb.dim_which_sub_fem[i] + teb.nb_sub_fem*j ;
}

FESumConstruct::FESumConstruct(int kk,const TypeOfFE **t)
 :k(kk),teb(t),NN(new int[kk+1]),DF(new int[kk+1]) , comp(new int[kk])
{  
   map<const TypeOfFE *,int> m;
   int i=k,j;    
   while(i--) // on va a l'envert pour avoir comp[i] <=i 
      m[teb[i]]=i;
    // l'ordre comp est important comp est croissant  mais pas de pb. 
   i=k;    
   while(i--) 
     comp[i]=m[teb[i]]; //  comp[i] <=i
        
  // reservatition des intervalles en espaces
  int n=0,N=0;
   for ( j=0;j<kk;j++)
     {NN[j]=N;N+=teb[j]->N;}
   NN[kk] = N;
 //  reservation des interval en df   
   n=0;
   for ( j=0;j<kk;j++)
    { DF[j]=n;n+=teb[j]->NbDoF;}
   DF[kk] = n;
//  n = nb de DF total   
//  N the fem is in R^N 
   
  data = new int [n*5 + N];
  int c=0;
  int ki= 0; 
// recherche des noeuds
   KN<int> w(7),nn(7);
   w=0;
   nn=0; 
   
         
   for ( j=0;j<kk;j++)
     for ( i=0;i<teb[j]->NbDoF;i++)
         nn[teb[j]->DFOnWhat[i]]++;
   nbn=0;      
   for( j=0;j<7;j++)
     if (nn[j]) nn[j]=nbn++;
     else nn[j]=-1;
   KN<int> dln(7);
   dln=0;
  // nn donne numero de noeud sur what            
   for ( j=0;j<kk;j++)
     for ( i=0;i<teb[j]->NbDoF;i++)
       data[c++] = teb[j]->DFOnWhat[i];
    
   for ( j=0;j<kk;j++)
    {
     int  cc=c;
     for ( i=0;i<teb[j]->NbDoF;i++)
       data[c++] = teb[j]->DFOfNode[i]+dln[teb[j]->DFOnWhat[i]];
     for ( i=0;i<teb[j]->NbDoF;i++)
       dln[teb[j]->DFOnWhat[i]]=Max(dln[teb[j]->DFOnWhat[i]],data[cc++]+1);      
    }
        
       
   for ( j=0;j<kk;j++)
    { 
     //  w renumerotation des noeuds 
     //  Ok si un noeud par what 
     for ( i=0;i<teb[j]->NbDoF;i++)
       data[c++] = nn[teb[j]->DFOnWhat[i]];
    }
    
   for ( j=0;j<kk;j++)
     for ( i=0;i<teb[j]->NbDoF;i++)
       data[c++] = j; //  node from of FE

    
   for ( j=0;j<kk;j++)
     for ( i=0;i<teb[j]->NbDoF;i++)
       data[c++] = i; //  node from of df in FE
       
   int xx=0;
   for (j=0;j<kk;j++)
     { 
      int xxx=xx;
      for (i=0;i<teb[j]->N;i++)
       { 
         data[c] = teb[j]->dim_which_sub_fem[i]+xx;
         xxx=Max(xxx,data[c]+1);
         c++;
       }
       xx=xxx;
     }
  
  throwassert(c== 5*n+N);      
/*  int cc=0;
   cout << " Data : " << endl;
  for ( i=0;i<5;i++)    {
    for (j=0;j<n;j++)
      cout << " " << data[cc++];
     cout << endl;}
 cout << " which " ;
 for (i=0;i<N;i++)
   cout << " " << data[cc++];
  cout << endl;*/
}

class TypeOfFE_P1Lagrange : public  TypeOfFE { public:  
  static int Data[];
  static double Pi_h_coef[];
   TypeOfFE_P1Lagrange(): TypeOfFE(1,0,0,1,Data,1,1,3,3,Pi_h_coef)
    { const R2 Pt[] = { R2(0,0), R2(1,0), R2(0,1) }; 
      for (int i=0;i<NbDoF;i++) {
       pij_alpha[i]= IPJ(i,i,0);
       P_Pi_h[i]=Pt[i]; }
     }
  // void FB(const Mesh & Th,const Triangle & K,const R2 &P, RNMK_ & val) const;
   void FB(const bool * whatd,const Mesh & Th,const Triangle & K,const R2 &P, RNMK_ & val) const;
   
//   void D2_FB(const Mesh & Th,const Triangle & K,const R2 &P, RNMK_ & val) const;
  // void Pi_h(const baseFElement & K,RN_ & val, InterpolFunction f, R* v,int, void *) const;
virtual R operator()(const FElement & K,const  R2 & PHat,const KN_<R> & u,int componante,int op) const ;
   
} ;

///////////////////////////////////////////////////////////////////////////////
////////////////////////////////// NEW ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

class TypeOfFE_P1Bubble : public  TypeOfFE { public:  
  static int Data[];
  static double Pi_h_coef[];
   TypeOfFE_P1Bubble(): TypeOfFE(1,0,1,1,Data,1,1,4,4,Pi_h_coef)
    { const R2 Pt[] = { R2(0,0), R2(1,0), R2(0,1), R2(1./3.,1./3.) }; 
      for (int i=0;i<NbDoF;i++) {
       pij_alpha[i]= IPJ(i,i,0);
       P_Pi_h[i]=Pt[i]; }
     }
  // void FB(const Mesh & Th,const Triangle & K,const R2 &P, RNMK_ & val) const;
   void FB(const bool * whatd,const Mesh & Th,const Triangle & K,const R2 &P, RNMK_ & val) const;
   
//   void D2_FB(const Mesh & Th,const Triangle & K,const R2 &P, RNMK_ & val) const;
 //  void Pi_h(const baseFElement & K,RN_ & val, InterpolFunction f, R* v,int, void *) const;
  //virtual R operator()(const FElement & K,const  R2 & PHat,const KN_<R> & u,int componante,int op) const ;
   
} ;
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

class TypeOfFE_P2Lagrange : public  TypeOfFE { public:  
  static int Data[];
  static double Pi_h_coef[];
  
   TypeOfFE_P2Lagrange(): TypeOfFE(1,1,0,1,Data,3,1,6,6,Pi_h_coef)
    { const R2 Pt[] = { R2(0,0), R2(1,0), R2(0,1),R2(0.5,0.5),R2(0,0.5),R2(0.5,0) };
      for (int i=0;i<NbDoF;i++) {
       pij_alpha[i]= IPJ(i,i,0);
       P_Pi_h[i]=Pt[i]; }
     }
   
  // void FB(const Mesh & Th,const Triangle & K,const R2 &P, RNMK_ & val) const;
   void FB(const bool * whatd,const Mesh & Th,const Triangle & K,const R2 &P, RNMK_ & val) const;
 //  void D2_FB(const Mesh & Th,const Triangle & K,const R2 &P, RNMK_ & val) const;
  // void Pi_h(const baseFElement & K,RN_ & val, InterpolFunction f, R* v,int, void *) const;
} ;

int TypeOfFE_P1Lagrange::Data[]={0,1,2,       0,0,0,       0,1,2,       0,0,0,        0,1,2,       0};
int TypeOfFE_P1Bubble::Data[]={0,1,2,6,     0,0,0,0,     0,1,2,3,     0,0,0,0,        0,1,2,3,     0};
int TypeOfFE_P2Lagrange::Data[]={0,1,2,3,4,5, 0,0,0,0,0,0, 0,1,2,3,4,5, 0,0,0,0,0,0,  0,1,2,3,4,5, 0};
double TypeOfFE_P1Lagrange::Pi_h_coef[]={1.,1.,1.};
double TypeOfFE_P1Bubble::Pi_h_coef[]={1.,1.,1.,1.};
double TypeOfFE_P2Lagrange::Pi_h_coef[]={1.,1.,1.,1.,1.,1.};

inline void dump(char *m,int n,int * p)
{
  cout << m ;
  for (int i=0;i<n;i++) cout << " " << p[i] ;
  cout << endl;
}



ConstructDataFElement::~ConstructDataFElement()
{
  if(*counter==0) 
   {
    delete [] NodesOfElement;
    delete []  FirstNodeOfElement;
    delete [] FirstDfOfNode;
  }
 else counter--;
}

 ConstructDataFElement::ConstructDataFElement(const ConstructDataFElement * t,int k)
  :thecounter(0), 
   counter(t->counter),
   MaxNbNodePerElement(t->MaxNbNodePerElement),
   MaxNbDFPerElement(t->MaxNbDFPerElement*k),
   NodesOfElement(t->NodesOfElement),
   FirstNodeOfElement(t->FirstNodeOfElement),
   FirstDfOfNode(0),
   NbOfElements(t->NbOfElements),
   NbOfDF(t->NbOfDF*k),
   NbOfNode(t->NbOfNode),
   Nproduit(t->Nproduit*k)
 {
   throwassert(t==0 || t->FirstDfOfNode==0);
   *counter++;
 }

ConstructDataFElement::ConstructDataFElement (const Mesh &Th,/*int NbDfOnSommet,int NbDfOnEdge,int NbDfOnElement*/
const  KN<const TypeOfFE *> & TFEs,const TypeOfMortar *tm,
int nbdfv,const int *ndfv,int nbdfe,const int *ndfe)
: counter(&thecounter),thecounter(0) 
{ 
 Make(Th,TFEs,/*NbDfOnSommet,NbDfOnEdge,NbDfOnElement,*/ tm,nbdfv,ndfv,nbdfe,ndfe);
}

ConstructDataFElement::ConstructDataFElement(const FESpace ** l,int k,const KN<const TypeOfFE *>  & TFEs) 
: thecounter(0),counter(&thecounter) 
{
 int NbDfOnSommet=0;
 int NbDfOnEdge=0;
 int NbDfOnElement=0;
 const Mesh & Th(l[0]->Th);
 for  (int i=0;i<k;i++)
   {
     NbDfOnSommet += l[i]->TFE[0]->NbDfOnVertex;
     NbDfOnEdge += l[i]->TFE[0]->NbDfOnEdge;
     NbDfOnElement += l[i]->TFE[0]->NbDfOnElement;
     throwassert( &Th== &l[i]->Th); 
     throwassert( l[i]->TFE.constant());
   }
   
 Make(Th,TFEs);//NbDfOnSommet,NbDfOnEdge,NbDfOnElement,0);   
}
 

void ConstructDataFElement::Make(const Mesh &Th,
const  KN<const TypeOfFE *> & TFEs,

/*int NbDfOnSommet,int NbDfOnEdge,int NbDfOnElement,*/const TypeOfMortar *tm,
int nb_dfv,const int *ndfv,int nb_dfe,const int *ndfe) 
/* pour le condition de periodicité ..
  nb_dfv :  nombre de sommets recolle par periodicité
  ndfv:  numerotation des sommets pour CL de periodicite
  ndfv[i] = numero du sommet i (i<Th.nv)  recolle par periodicité
  nb_dfe:  nombre de arete frontiere  recolle par periodicité
  ndfe[i]:  numero de l'arete frontiere (i < Th.neb) recollé. 
  F. H. 
*/  
{
  *counter=0;
  assert(TFEs.constant());
  const TypeOfFE & TFE(*TFEs[0]);
  int nbdfe=TFE.NbDoF; 
  int NbDfOnSommet=TFE.NbDfOnVertex;
  int NbDfOnEdge=TFE.NbDfOnEdge;
  int NbDfOnElement=TFE.NbDfOnElement;
  int NbNodeonVertex=0;
  int NbNodeonEdge=0;
  int NbNodeonElement=0;
  int NbNodes=TFE.NbNode;


   assert( nbdfe == 3*NbDfOnSommet+3*NbDfOnEdge+NbDfOnElement);
/*  Vieux code 
  { // construct du tableau NodeOn et calcul de NbNodeon.. 
   KN<int>NodeOn(NbNodes);
   NodeOn=-1;
   int nb[7];
   for (int i=0;i<7;i++)
     nb[i]=0;
   int kkk=0;
   for (int df=0;df<TFE.NbDoF;df++)
    {
     int node=TFE.NodeOfDF[df];
     int w=TFE.DFOnWhat[df];
     if ( NodeOn[node] >=0)
      assert( NodeOn[node] ==w); 
     else {
       NodeOn[node]=w;
  	     ++kkk;
  	     ++nb[w]; }// on vertex 0  	    
  	 }
   assert(nb[0]==nb[1] && nb[1] == nb[2]);
   assert(nb[3]==nb[4] && nb[4] == nb[5]);
   NbNodeonVertex=nb[0];
   NbNodeonEdge=nb[3];
   NbNodeonElement=nb[6]; 
   assert(kkk==NbNodes);
  }*/
  
   KN<int> NbDFonNode(NbNodes), NodeIsOn(NbNodes);
   NbDFonNode=0;
   for (int df=0;df<nbdfe;df++)
    {
     int node=TFE.NodeOfDF[df];
     int what=TFE.DFOnWhat[df];
     int ndfonnode = TFE.DFOfNode[df];
     NbDFonNode[node]=Max(NbDFonNode[node],ndfonnode+1);
     NodeIsOn[node]=what;
    }
   assert(NbDFonNode.sum() == nbdfe);
   
   NbNodeonVertex=TFE.NbNodeOnVertex;
   NbNodeonEdge=TFE.NbNodeOnEdge;
   NbNodeonElement=TFE.NbNodeOnElement; 
  
  assert(NbNodeonVertex<=1);
  assert(NbNodeonEdge<=1);
  assert(NbNodeonElement<=1);
  Nproduit =1;
  const int ndf=NbDFonNode[0];
  int samendf=1;
  
  for (int i=1;i<NbNodes;i++)
    if ( ndf != NbDFonNode[i]) 
     { samendf = 0;
       break;}
  
  int nbne = NbNodes;
  int nn=0;  
  int firstmul=0;
  throwassert( tm || Th.NbMortars==0);
  NbOfElements = Th.nt; //  by default
  // if mortars 
  //
  int NbOfNodeL=0;
  NbOfElements += Th.NbMortars;  
  FirstDfOfNode =0;
  FirstNodeOfElement=0;
  MaxNbDFPerElement=nbdfe; 
  assert(3*NbDfOnSommet+3*NbDfOnEdge+NbDfOnElement==MaxNbDFPerElement);

  int ks=TFE.NbNodeOnVertex>0,
      ke=TFE.NbNodeOnEdge>0,
      kt=TFE.NbNodeOnElement>0;
 /* Vieux code 
   if(NbDfOnSommet) { nbne+=3;
     ks=1;
     ndf=NbDfOnSommet;}

  if(NbDfOnEdge) {  nbne+=3;
     ke=1;
     samendf &= !ndf || ndf == NbDfOnEdge;
     ndf=NbDfOnEdge;}
     
  if(NbDfOnElement) {  nbne+=1;
     kt=1;
     samendf &= !ndf || ndf == NbDfOnElement;
     ndf=NbDfOnElement;}

  int NbDFonNode[7],NodeIsOn[7];
   {
     int j=0,k=0;
     if(ks) { NbDFonNode[j++]=NbDfOnSommet; NbDFonNode[j++]=NbDfOnSommet; NbDFonNode[j++]=NbDfOnSommet;}
     if(ke) { NbDFonNode[j++]=NbDfOnEdge; NbDFonNode[j++]=NbDfOnEdge; NbDFonNode[j++]=NbDfOnEdge;}
     if(kt) { NbDFonNode[j++]=NbDfOnElement;}

     if (ks) {NodeIsOn[k++]=0;NodeIsOn[k++]=1;NodeIsOn[k++]=2;}
     if (ke) {NodeIsOn[k++]=3;NodeIsOn[k++]=4;NodeIsOn[k++]=5;}
     if (kt) {NodeIsOn[k++]=6;}
     
     throwassert(j == nbne);
  }
*/     
    MaxNbNodePerElement=nbne;

//  
   if ( ks && (!ke && ! kt) && (ndfv==0 && ndfe==0))
     {nn=Th.nv;
      NodesOfElement=0;
      }
   else {
//    constuction du tableau  NodesOfElement  bofbof 

//     computation of the length lne of array NodesOfElement        
       int lne=  Th.nt*nbne;
       if (Th.NbMortars)
        {
          samendf= false;
          NbOfNodeL=Th.NbMortars;
          throwassert(tm);
          FirstNodeOfElement = new int[NbOfElements+1];
          int k=0,kk=0;          
          for (k=0;k<Th.nt;k++,kk+=nbne)
            FirstNodeOfElement[k] = kk;
           
          for (int im=0;im<Th.NbMortars;im++) // (www) code 
            {
             FirstNodeOfElement[k++]=lne;           
             lne += tm->NbOfNodes(Th,Th.mortars[im]);
            }
          FirstNodeOfElement[k++]=lne;           
        }       
        
       NodesOfElement = new int[lne];
       
       for (int i=0;i<lne;i++) 
         NodesOfElement[i]=-1;
       int i=0;
       int oe=0; 
       if(ks) { oe=3*NbNodeonVertex;nn= (ndfv ? nb_dfv : Th.nv)*NbNodeonVertex;}
       
       if (ke && ndfe) { // il y a des aretes avec cl de periodicite. 
        for (int be=0;be<Th.neb;be++)
         {
          int j,k=Th.BoundaryTriangle(be,j);
          int jj=j;
          int kk=Th.TriangleAdj(k,jj);
          for (int kkk=0;kkk<NbNodeonEdge;kkk++)
           {
             int kj=kkk,kjj=-1;
            // assert(kk != k);
             if ( kk >= 0 && jj >=0  &&  !(( kk == k ) && ( jj=j ) ) ) 
              if (k < kk ) kjj = NbNodeonEdge-kj-1; // 
              else   kj = NbNodeonEdge-kj-1, kjj = NbNodeonEdge-kj-1;
          
          if (kjj >=0)   
            NodesOfElement[kk*nbne+oe+jj] = nn + ndfe[be]*NbNodeonEdge+ kjj   ; // adj
          NodesOfElement[k*nbne+oe+j]   = nn + ndfe[be]*NbNodeonEdge+ kj   ; // new  
          }       
         }
         nn += nb_dfe;
        }
       for (int k=0;k<Th.nt;k++)
        {
          int iold=i;
          if(ks) {
            for (int j=0;j<3;j++)
             for (int jj=0;jj<NbNodeonVertex;jj++)
              NodesOfElement[i++]= (ndfv ?  ndfv[Th(k,j)] : Th(k,j))*NbNodeonVertex+jj ;
          }
          if(ke) {
            for (int j=0;j<3;j++)
             for (int ll=0;ll<NbNodeonEdge;ll++)
             if (NodesOfElement[i]<0) {
               int jj=j;
               int kk=Th.TriangleAdj(k,jj);
               int jjj=jj*NbNodeonEdge+NbNodeonEdge-ll-1; // autre sens  
               assert(kk>=k && jjj == jj); 
               NodesOfElement[kk*nbne+oe+jjj] = nn   ; // adj   
               NodesOfElement[i++]           = nn++ ; // new
             }  else i++;
          }
          for (int jj=0;jj<NbNodeonElement;jj++) 
             NodesOfElement[i++]= nn++;
         // cout << k ;
         // dump(" ",i-iold, NodesOfElement+iold);
        }
       // cout << i << " " << Th.nt*nbne << endl;
        firstmul=nn;
        if (Th.NbMortars)
          {  
            //  construction of the mortars element 
           int * color= new int [firstmul]; 
           //  
           int thecolor=0;
           for (int j=0;j<firstmul;j++) 
             color[j]=thecolor;
             
           for (int im=0;im<Th.NbMortars;im++)
            {   
             int iold=i;         
             thecolor++; // get a new color
             // tm->ConstructionOfNode(Th,im,NodesOfElement,FirstNodeOfElement,nn);  
             Mortar & M(Th.mortars[im]);
              NodesOfElement[i++] = nn++; // the first node is the lag. mul.
             int n=M.NbT();
             for (int kk=0;kk<n;kk++)
             {
              int K,e;
               K=M.T_e(kk,e); 
               int kb=FirstNodeOfElement[K];
               int ke=FirstNodeOfElement[K+1];             
               for (int j=kb,jj=0;j<ke;j++,jj++)
                if (onWhatIsEdge[e][NodeIsOn[jj]])
                 { // the node jj is on edge e
                    int node=NodesOfElement[j];
                   // cout << "." << jj << " K=" << K <<" "<<e << " n=" << node ;
                    throwassert(node<firstmul);
                    if (color[node] != thecolor) //  new node => adding
                     { 
                    //   cout << "a, ";
                       color[node] = thecolor;
                       NodesOfElement[i++] = node;
                     }
                   // else cout << ", ";
                 }
               }
               //cout << endl;
               
               //cout << im ;
               //dump(": ",i-iold, NodesOfElement+iold);
               throwassert(i==FirstNodeOfElement[im+Th.nt+1]);
              }
             delete [] color;
          } 
        else
          throwassert(i==Th.nt*nbne && i );
        NbOfNode=nn;
                
        
   }
  NbOfNode=nn;  
  int NbOfDFL=0;  
   if (! samendf) 
       {
         
         throwassert(NodesOfElement);
         FirstDfOfNode= new int [nn+1];
         for (int i=0;i<=nn;i++) FirstDfOfNode[i]=-1; 
         int i=0;
         //  the classical part (FEM)
         for (int k=0;k<Th.nt;k++)
           for (int j=0;j<nbne;j++) // thanks to student
             FirstDfOfNode[ NodesOfElement[i++]+1]=NbDFonNode[j];
         //  the mortars parts juste the mulplicator 
         
         for (int km=0,k=Th.nt;km<Th.NbMortars;km++,k++) 
            {  //  the lag. mult. is the first node of the mortar --
              throwassert(FirstNodeOfElement);
              //  hack 
              int fk=FirstNodeOfElement[k];
              int lk=FirstNodeOfElement[k+1];
              int ndlmul = tm->NbLagrangeMult(Th,Th.mortars[km]);  //  On node par 
              int nodemul = NodesOfElement[fk]; // the first node is the lagr. mul.
              throwassert(FirstDfOfNode[nodemul+1]==-1);      
              FirstDfOfNode[nodemul+1]= ndlmul;
              NbOfDFL += ndlmul;
              int nbdle=0;
              int nbnm=lk-fk;
              for (int j=fk;j<lk;j++)
               nbdle+=FirstDfOfNode[NodesOfElement[j]+1];
              MaxNbDFPerElement = Max(MaxNbDFPerElement,nbdle);
            }
             
         FirstDfOfNode[0]=0;
         for (int i=0;i<=nn;i++) throwassert(FirstDfOfNode[i]!=-1); 

         for (int i=0;i<nn;i++) 
              FirstDfOfNode[i+1] += FirstDfOfNode[i] ;
           NbOfDF=  FirstDfOfNode[nn];       
        }
    else
       {
         NbOfDF = nn*ndf; 
         Nproduit = ndf;        
       }
   MaxNbNodePerElement=nbne;
   
   cout << " Nb Of Nodes = " << nn << endl;   
   if(NbOfNodeL)        
     cout << " Nb of Lagrange Mul Node = " << NbOfNodeL << endl;        
   cout << " Nb of DF = " << NbOfDF << endl;   
   if(NbOfDFL) {  
     cout << " Nb of Lagrange Mul DF = "   << NbOfDFL << endl;  
     cout << " MaxNbDFPerElement     =   " << MaxNbDFPerElement << endl;
   };

}

FESpace::FESpace(const FESpace & Vh,int k )
 :
     ptrTFE(new TypeOfFEProduit(k,*Vh.TFE[0])),
     TFE(1,0,ptrTFE),
     cmesh(Vh.Th),
     cdef(Vh.cdef?new ConstructDataFElement(Vh.cdef,k):0),
     N(Vh.N*k),
     Nproduit(Vh.Nproduit*k),
     nb_sub_fem(TFE[0]->nb_sub_fem),
     dim_which_sub_fem(TFE[0]->dim_which_sub_fem),
     
     Th(Vh.Th),
     NbOfDF(Vh.NbOfDF*k),
     NbOfElements(Vh.NbOfElements),
     NbOfNodes(Vh.NbOfNodes),
     MaxNbNodePerElement(Vh.MaxNbNodePerElement),
     MaxNbDFPerElement(Vh.MaxNbDFPerElement*k),
     NodesOfElement(Vh.NodesOfElement),
     FirstDfOfNodeData(cdef?cdef->FirstDfOfNode:0),
     FirstNodeOfElement(Vh.FirstNodeOfElement),
     tom(0) {
     if(cdef) renum();}
 
FESpace::FESpace(const FESpace ** Vh,int k )
 :
     ptrTFE(new TypeOfFESum(Vh,k)),
     TFE(1,0,ptrTFE),
     cmesh((**Vh).Th),
     cdef(new ConstructDataFElement(Vh,k,TFE)),
     N(sum(Vh,&FESpace::N,k)),
     Nproduit(cdef->Nproduit),     
     nb_sub_fem(TFE[0]->nb_sub_fem),
     dim_which_sub_fem(TFE[0]->dim_which_sub_fem),

     Th((**Vh).Th),
     NbOfDF(cdef->NbOfDF),
     NbOfElements(cdef->NbOfElements),
     NbOfNodes(cdef->NbOfNode),
     MaxNbNodePerElement(cdef->MaxNbNodePerElement),
     MaxNbDFPerElement(cdef->MaxNbDFPerElement),
     NodesOfElement(cdef->NodesOfElement),
     FirstDfOfNodeData(cdef->FirstDfOfNode),
     FirstNodeOfElement(cdef->FirstNodeOfElement),
     tom(0) {
     if(cdef) renum(); }
     
FESpace::FESpace(const Mesh & TTh,const TypeOfFE ** tef,int k,int nbdfv,const int *ndfv,int nbdfe,const int *ndfe )
 :
     ptrTFE(new TypeOfFESum(tef,k)),
     TFE(1,0,ptrTFE),
     cmesh(TTh),
     cdef(new ConstructDataFElement(TTh,TFE,//sum(tef,&TypeOfFE::NbDfOnVertex,k),
                                       // sum(tef,&TypeOfFE::NbDfOnEdge,k),
                                        //sum(tef,&TypeOfFE::NbDfOnElement,k),
                                        0,nbdfv,ndfv,nbdfe,ndfe)),
     N(sum(tef,&TypeOfFE::N,k)),
     Nproduit(cdef->Nproduit),     
     nb_sub_fem(TFE[0]->nb_sub_fem),
     dim_which_sub_fem(TFE[0]->dim_which_sub_fem),

     Th(TTh),
     NbOfDF(cdef->NbOfDF),
     NbOfElements(cdef->NbOfElements),
     NbOfNodes(cdef->NbOfNode),
     MaxNbNodePerElement(cdef->MaxNbNodePerElement),
     MaxNbDFPerElement(cdef->MaxNbDFPerElement),
     NodesOfElement(cdef->NodesOfElement),
     FirstDfOfNodeData(cdef->FirstDfOfNode),
     FirstNodeOfElement(cdef->FirstNodeOfElement),
     tom(0) {
     if(cdef) renum(); }


 FESpace::FESpace(const Mesh & TTh,const TypeOfFE & tef,int nbdfv,const int *ndfv,int nbdfe,const int *ndfe)
   :  
     ptrTFE(0),
     TFE(1,0,&tef),
     cmesh(TTh),
     cdef(new ConstructDataFElement(TTh,TFE,0,nbdfv,ndfv,nbdfe,ndfe)),
     //tef.NbDfOnVertex,tef.NbDfOnEdge,tef.NbDfOnElement,0,nbdfv,ndfv,nbdfe,ndfe)),
     N(tef.N),
     Nproduit(cdef->Nproduit),
     nb_sub_fem(TFE[0]->nb_sub_fem),
     dim_which_sub_fem(TFE[0]->dim_which_sub_fem),
     Th(TTh),
     NbOfDF(cdef->NbOfDF),
     NbOfElements(cdef->NbOfElements),
     NbOfNodes(cdef->NbOfNode),
     MaxNbNodePerElement(cdef->MaxNbNodePerElement),
     MaxNbDFPerElement(cdef->MaxNbDFPerElement),
     NodesOfElement(cdef->NodesOfElement),
     FirstDfOfNodeData(cdef->FirstDfOfNode),
     FirstNodeOfElement(cdef->FirstNodeOfElement),
     tom(0) {
     if(tef.NbDfOnVertex || tef.NbDfOnEdge) renum();
     }
     
 FESpace::~FESpace()
   {
     SHOWVERB(cout << " FESpace::~FESpace() " << endl);
      delete  cdef;
      if(ptrTFE) 
        delete  ptrTFE;
   }

 FESpace::FESpace(const Mesh & TTh,const TypeOfFE & tef,const TypeOfMortar & tm)
   :  
     ptrTFE(0),
     TFE(1,0,&tef),
     cmesh(TTh),
     cdef(new ConstructDataFElement(TTh,TFE,&tm)),//tef.NbDfOnVertex,tef.NbDfOnEdge,tef.NbDfOnElement,&tm)),
     N(tef.N),
     Nproduit(1),
     nb_sub_fem(TFE[0]->nb_sub_fem),
     dim_which_sub_fem(TFE[0]->dim_which_sub_fem),     
     Th(TTh),
     NbOfDF(cdef->NbOfDF),
     NbOfElements(cdef->NbOfElements),
     NbOfNodes(cdef->NbOfNode),
     MaxNbNodePerElement(cdef->MaxNbNodePerElement),
     MaxNbDFPerElement(cdef->MaxNbDFPerElement),
     NodesOfElement(cdef->NodesOfElement),
     FirstDfOfNodeData(cdef->FirstDfOfNode),
     FirstNodeOfElement(cdef->FirstNodeOfElement),
     tom(&tm) { 
     // cout << "avant renum ="<< *this <<endl;
       renum();
     // cout << "apres renum ="<< *this <<endl;
     }
     
void ConstructDataFElement::renum(const long *r,int l)   
 { 
/*   cout << "renu=" << l << ":" << endl;
   for (int i=0;i<NbOfNode;i++)
      if (i%10) cout << r[i] << "\t";
      else cout << "\n " << i << ":\t" << r[i] << "\t";
      cout << endl; 
*/
   throwassert(this);
   if (NodesOfElement) 
     for (int i=0;i< l ; i++)
       NodesOfElement[i]=r[NodesOfElement[i]];
   if(FirstDfOfNode)  
    { int k,i,*n=new int[NbOfNode];
      for ( i=0;i<NbOfNode;i++)        
         n[r[i]]=FirstDfOfNode[i+1]-FirstDfOfNode[i]; 
      FirstDfOfNode[0]=k=0;
      for(i=0;i< NbOfNode;)
        {k+=n[i];
         FirstDfOfNode[++i]=k;}
       delete [] n;
    }
 }     
 
 
/*
 void TypeOfFEProduit::D2_FB(const Mesh & Th,const Triangle & K,const R2 & P,RNMK_ & val) const
 {
   int n=teb.NbDoF;
   int m=teb.N;   
   val=0.0;
   SubArray t(3);
   RNMK_ v(val(SubArray(n,0,k),SubArray(m),t));
   teb.D2_FB(Th,K,P,v);
   for (int i=1;i<k;i++)
     val(SubArray(n,i,k),SubArray(m,m*i),t)=v; 
 } 
*/ 
/*
 void TypeOfFEProduit::FB(const Mesh & Th,const Triangle & K,const R2 & P,RNMK_ & val) const
 {
   int n=teb.NbDoF;
   int m=teb.N;   
   val=0.0;
   SubArray t(3);
   RNMK_ v(val(SubArray(n,0,k),SubArray(m),t));
   teb.FB(Th,K,P,v);
   for (int i=1;i<k;i++)
     val(SubArray(n,i,k),SubArray(m,m*i),t)=v; 
 }
 */
 void TypeOfFEProduit::FB(const bool * whatd,const Mesh & Th,const Triangle & K,const R2 & P,RNMK_ & val) const
 {
   int n=teb.NbDoF;
   int m=teb.N;   
   val=0.0;
   SubArray t(val.K());
   RNMK_ v(val(SubArray(n,0,k),SubArray(m),t));
   teb.FB(whatd,Th,K,P,v);
   for (int i=1;i<k;i++)
     val(SubArray(n,i,k),SubArray(m,m*i),t)=v; 
 }
 
/* 
 void TypeOfFEProduit::Pi_h(const baseFElement & K,RN_ & val, InterpolFunction f, R* v,int j, void * arg) const
 {
   const baseFElement  KK(K,teb);
   int m=teb.N;   
    for(int i=0;i<k;i++)
     { RN_  vv(val(SubArray(m,m*i)));
     teb.Pi_h(KK,vv,f,v,j+i*m,arg);}
 }
 
/* 
 void TypeOfFESum::D2_FB(const Mesh & Th,const Triangle & K,const R2 & P,RNMK_ & val) const
 {
   val=0.0;
   SubArray t(3);
   for (int i=0;i<k;i++)
    {
     int j=comp[i];
     int ni=NN[i];
     int di=DF[i];  
     int i1=i+1; 
     int nii=NN[i1];
     int dii=DF[i1];
     throwassert(ni<nii && di < dii);
     RNMK_ v(val(SubArray(dii-di,di),SubArray(nii-ni,ni),t));     
     if (j<=i)
       teb[i]->D2_FB(Th,K,P,v);       
     else
       v=val(SubArray(DF[j+1]-DF[j],DF[j]),SubArray(NN[j+1]-NN[j],NN[j]),t);     
    } }
*/    
/*
 void TypeOfFESum::FB(const Mesh & Th,const Triangle & K,const R2 & P,RNMK_ & val) const
 {
   val=0.0;
   SubArray t(3);
   for (int i=0;i<k;i++)
    {
     int j=comp[i];
     int ni=NN[i];
     int di=DF[i];  
     int i1=i+1; 
     int nii=NN[i1];
     int dii=DF[i1];
     throwassert(ni<nii && di < dii);
     RNMK_ v(val(SubArray(dii-di,di),SubArray(nii-ni,ni),t));     
     if (j<=i)
       teb[i]->FB(Th,K,P,v);       
     else
       v=val(SubArray(DF[j+1]-DF[j],DF[j]),SubArray(NN[j+1]-NN[j],NN[j]),t);     
    }
 }
*/ 
  void TypeOfFESum::FB(const bool * whatd,const Mesh & Th,const Triangle & K,const R2 & P,RNMK_ & val) const
 {
   val=0.0;
   SubArray t(val.K());
   for (int i=0;i<k;i++)
    {
     int j=comp[i];
     int ni=NN[i];
     int di=DF[i];  
     int i1=i+1; 
     int nii=NN[i1];
     int dii=DF[i1];
     throwassert(ni<nii && di < dii);
     RNMK_ v(val(SubArray(dii-di,di),SubArray(nii-ni,ni),t));     
     if (j<=i)
       teb[i]->FB(whatd,Th,K,P,v);       
     else
       v=val(SubArray(DF[j+1]-DF[j],DF[j]),SubArray(NN[j+1]-NN[j],NN[j]),t);     
    }
 }
/* 
 void TypeOfFESum::Pi_h(const baseFElement & K,RN_ & val, InterpolFunction f, R* v,int jjj, void * arg) const
 {
    for(int i=0;i<k;i++)
     { 
       const baseFElement  KK(K,*teb[i]);
       int dfii=DF[i+1],dfi=DF[i];
       RN_  vv(val(SubArray(dfii-dfi,dfi)));
       teb[i]->Pi_h(KK,vv,f,v,jjj+NN[i],arg);
     }
    // cout << val(SubArray(NbDoF)) << endl;
     
 }
 /*
 void TypeOfFE_P1Lagrange::D2_FB(const Mesh & ,const Triangle & ,const R2 & ,RNMK_ & val) const
{ //  
  val=0;
}
*/
/*
 void TypeOfFE_P2Lagrange::D2_FB(const Mesh & ,const Triangle & K,const R2 & P,RNMK_ & val) const
{ // 2 times derivatives  for error indicator
//  const Triangle & K(FE.T);
  R2 A(K[0]), B(K[1]),C(K[2]);
  R l0=1-P.x-P.y,l1=P.x,l2=P.y; 
  R2 Dl0(K.H(0)), Dl1(K.H(1)), Dl2(K.H(2));
  R l4_0=(4*l0-1),l4_1=(4*l1-1),l4_2=(4*l2-1); 
  
  throwassert(val.N() >=6);
  throwassert(val.M()==1 );
  throwassert(val.K()==3 );
  
  val=0; 
  RN_ fxx(val('.',0,0)); 
  RN_ fxy(val('.',0,1)); 
  RN_ fyy(val('.',0,2)); 
  
  fxx[0] = 4*Dl0.x*Dl0.x;
  fxx[1] = 4*Dl1.x*Dl1.x;
  fxx[2] = 4*Dl2.x*Dl2.x;
  fxx[3] =  8*Dl1.x*Dl2.x;
  fxx[4] =  8*Dl0.x*Dl2.x;
  fxx[5] =  8*Dl0.x*Dl1.x;

  fyy[0] = 4*Dl0.y*Dl0.y;
  fyy[1] = 4*Dl1.y*Dl1.y;
  fyy[2] = 4*Dl2.y*Dl2.y;
  fyy[3] =  8*Dl1.y*Dl2.y;
  fyy[4] =  8*Dl0.y*Dl2.y;
  fyy[5] =  8*Dl0.y*Dl1.y;

  fxy[0] = 4*Dl0.y*Dl0.y;
  fxy[1] = 4*Dl1.y*Dl1.y;
  fxy[2] = 4*Dl2.y*Dl2.y;
  fxy[3] =  4*(Dl1.x*Dl2.y + Dl1.y*Dl2.x);
  fxy[4] =  4*(Dl0.x*Dl2.y + Dl0.y*Dl2.x);
  fxy[5] =  4*(Dl0.x*Dl1.y + Dl0.y*Dl1.x);

}
*/
 R TypeOfFE_P1Lagrange::operator()(const FElement & K,const  R2 & PHat,const KN_<R> & u,int componante,int op) const 
{ 
   R u0(u(K(0))), u1(u(K(1))), u2(u(K(2)));
   R r=0;
   if (op==0)
    {
      R l0=1-PHat.x-PHat.y,l1=PHat.x,l2=PHat.y; 
      r = u0*l0+u1*l1+l2*u2;
    }
   else
    { 
       const Triangle & T=K.T;
       R2 D0 = T.H(0) , D1 = T.H(1)  , D2 = T.H(2) ;
       if (op==1)
         r =  D0.x*u0 + D1.x*u1 + D2.x*u2 ;
        else 
         r =  D0.y*u0 + D1.y*u1 + D2.y*u2 ;
    }
 //  cout << r << "\t";
   return r;
}


void TypeOfFE_P1Lagrange::FB(const bool *whatd,const Mesh & ,const Triangle & K,const R2 & P,RNMK_ & val) const
{
//  const Triangle & K(FE.T);
  R2 A(K[0]), B(K[1]),C(K[2]);
  R l0=1-P.x-P.y,l1=P.x,l2=P.y; 
  
  if (val.N() <3) 
   throwassert(val.N() >=3);
  throwassert(val.M()==1 );
//  throwassert(val.K()==3 );
  
  val=0; 
  RN_ f0(val('.',0,op_id)); 
  
  if (whatd[op_id]) 
   {
    f0[0] = l0;
    f0[1] = l1;
    f0[2] = l2;}
 if (whatd[op_dx] || whatd[op_dy])
  {
  R2 Dl0(K.H(0)), Dl1(K.H(1)), Dl2(K.H(2));
  
  if (whatd[op_dx]) 
   {
    RN_ f0x(val('.',0,op_dx)); 
   f0x[0] = Dl0.x;
   f0x[1] = Dl1.x;
   f0x[2] = Dl2.x;
  }
  
  if (whatd[op_dy]) {
    RN_ f0y(val('.',0,op_dy)); 
   f0y[0] = Dl0.y;
   f0y[1] = Dl1.y;
   f0y[2] = Dl2.y;
  }
  }
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////// NEW /////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/*
void TypeOfFE_P1Bubble::FB(const Mesh & Th,const Triangle & K,const R2 &P, RNMK_ & val) const
{
  assert(0);
}


void TypeOfFE_P1Bubble::Pi_h(const baseFElement & K,RN_ & val, InterpolFunction f, R* v,int, void *) const
{
  assert(0);
}

/*
R TypeOfFE_P1Bubble::operator()(const FElement & K,const  R2 & PHat,const KN_<R> & u,int componante,int op) const
{
  assert(0);
}
*/

void TypeOfFE_P1Bubble::FB(const bool *whatd,const Mesh & ,const Triangle & K,const R2 & P,RNMK_ & val) const
{
//  const Triangle & K(FE.T);
  R2 A(K[0]), B(K[1]),C(K[2]);
  R l0=1-P.x-P.y, l1=P.x, l2=P.y, lb=l0*l1*l2*9.; 
  
  if (val.N() <4) 
   throwassert(val.N() >=4);
  throwassert(val.M()==1 );
//  throwassert(val.K()==3 );
  
  val=0; 
  RN_ f0(val('.',0,op_id)); 
  
  if (whatd[op_id]) 
   {
    f0[0] = l0-lb;
    f0[1] = l1-lb;
    f0[2] = l2-lb;
    f0[3] = 3.*lb;
   }
  if(  whatd[op_dx] || whatd[op_dy] || whatd[op_dxx] || whatd[op_dyy] ||  whatd[op_dxy])
 {
  R2 Dl0(K.H(0)), Dl1(K.H(1)), Dl2(K.H(2)), 
     Dlb((Dl0*l1*l2+Dl1*l0*l2+Dl2*l0*l1)*9.);
  
  if (whatd[op_dx]) 
   {
    RN_ f0x(val('.',0,op_dx)); 
   f0x[0] = Dl0.x-Dlb.x;
   f0x[1] = Dl1.x-Dlb.x;
   f0x[2] = Dl2.x-Dlb.x;
   f0x[3] = 3.*Dlb.x;
  }
  
  if (whatd[op_dy]) {
    RN_ f0y(val('.',0,op_dy)); 
   f0y[0] = Dl0.y-Dlb.y;
   f0y[1] = Dl1.y-Dlb.y;
   f0y[2] = Dl2.y-Dlb.y;
   f0y[3] = 3.*Dlb.y;
  }
 if (whatd[op_dxx])
  {  
    RN_ fxx(val('.',0,op_dxx)); 
    R lbdxx= 2*((Dl0.x*Dl1.x)*l2+(Dl1.x*Dl2.x)*l0+(Dl2.x*Dl0.x)*l1);
    fxx[0] = -lbdxx;
    fxx[1] = -lbdxx;
    fxx[2] = -lbdxx;
    fxx[3] = 3*lbdxx;
  }

 if (whatd[op_dyy])
  {  
    RN_ fyy(val('.',0,op_dyy));
    R lbdyy= 2*((Dl0.y*Dl1.y)*l2+(Dl1.y*Dl2.y)*l0+(Dl2.y*Dl0.y)*l1);
     
    fyy[0] =  -lbdyy;
    fyy[1] =  -lbdyy;
    fyy[2] =  -lbdyy;
    fyy[3] =  3*lbdyy;
  }
 if (whatd[op_dxy])
  {  
    assert(val.K()>op_dxy);
    RN_ fxy(val('.',0,op_dxy)); 
    R lbdxy= (Dl0.x*Dl1.y+ Dl0.y*Dl1.x)*l2+(Dl1.x*Dl2.y+Dl1.y*Dl2.x)*l0+(Dl2.x*Dl0.y+Dl2.y*Dl0.x)*l1;  
    fxy[0] = 4*Dl0.x*Dl0.y-9.*(l0-l1-l2);
    fxy[1] = 4*Dl1.x*Dl1.y-9.*(l0-l1-l2);
    fxy[2] = 4*Dl2.x*Dl2.y-9.*(l0-l1-l2);
    fxy[3] = 27.*(l0-l1-l2);
  }
  }

}
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/* void TypeOfFE_P1Lagrange::FB(const Mesh & ,const Triangle & K,const R2 & P,RNMK_ & val) const
{
//  const Triangle & K(FE.T);
  R2 A(K[0]), B(K[1]),C(K[2]);
  R l0=1-P.x-P.y,l1=P.x,l2=P.y; 
  R2 Dl0(K.H(0)), Dl1(K.H(1)), Dl2(K.H(2));
  
  if (val.N() <3) 
   throwassert(val.N() >=3);
  throwassert(val.M()==1 );
  throwassert(val.K()==3 );
  
  val=0; 
  RN_ f0(val('.',0,0)); 
  RN_ f0x(val('.',0,1)); 
  RN_ f0y(val('.',0,2)); 
  
  f0[0] = l0;
  f0[1] = l1;
  f0[2] = l2;
  
  f0x[0] = Dl0.x;
  f0x[1] = Dl1.x;
  f0x[2] = Dl2.x;
  
  f0y[0] = Dl0.y;
  f0y[1] = Dl1.y;
  f0y[2] = Dl2.y;
}

 void TypeOfFE_P2Lagrange::FB(const Mesh & ,const Triangle & K,const R2 & P,RNMK_ & val) const
{
//  const Triangle & K(FE.T);
  R2 A(K[0]), B(K[1]),C(K[2]);
  R l0=1-P.x-P.y,l1=P.x,l2=P.y; 
  R2 Dl0(K.H(0)), Dl1(K.H(1)), Dl2(K.H(2));
  R l4_0=(4*l0-1),l4_1=(4*l1-1),l4_2=(4*l2-1); 
  
//  throwassert(FE.N == 1);  
  throwassert( val.N()>=6);
  throwassert(val.M()==1);
  throwassert(val.K()==3 );
  
  val=0; 
  RN_ f0(val('.',0,0)); 
  RN_ f0x(val('.',0,1)); 
  RN_ f0y(val('.',0,2)); 
// --     
  f0[0] = l0*(2*l0-1);
  f0[1] = l1*(2*l1-1);
  f0[2] = l2*(2*l2-1);
  f0[3] = 4*l1*l2; // oppose au sommet 0
  f0[4] = 4*l0*l2; // oppose au sommet 1
  f0[5] = 4*l1*l0; // oppose au sommet 3
  
  
  f0x[0] = Dl0.x*l4_0;
  f0x[1] = Dl1.x*l4_1;
  f0x[2] = Dl2.x*l4_2;
  f0x[3] = 4*(Dl1.x*l2 + Dl2.x*l1) ;
  f0x[4] = 4*(Dl2.x*l0 + Dl0.x*l2) ;
  f0x[5] = 4*(Dl0.x*l1 + Dl1.x*l0) ;
  
  
  f0y[0] = Dl0.y*l4_0;
  f0y[1] = Dl1.y*l4_1;
  f0y[2] = Dl2.y*l4_2;
  f0y[3] = 4*(Dl1.y*l2 + Dl2.y*l1) ;
  f0y[4] = 4*(Dl2.y*l0 + Dl0.y*l2) ;
  f0y[5] = 4*(Dl0.y*l1 + Dl1.y*l0) ;
  
}
*/
void TypeOfFE_P2Lagrange::FB(const bool *whatd,const Mesh & ,const Triangle & K,const R2 & P,RNMK_ & val) const
{
//  const Triangle & K(FE.T);
  R2 A(K[0]), B(K[1]),C(K[2]);
  R l0=1-P.x-P.y,l1=P.x,l2=P.y; 
  R l4_0=(4*l0-1),l4_1=(4*l1-1),l4_2=(4*l2-1); 
  
//  throwassert(FE.N == 1);  
  throwassert( val.N()>=6);
  throwassert(val.M()==1);
//  throwassert(val.K()==3 );
  
  val=0; 
// --     
 if (whatd[op_id])
  {
   RN_ f0(val('.',0,op_id)); 
  f0[0] = l0*(2*l0-1);
  f0[1] = l1*(2*l1-1);
  f0[2] = l2*(2*l2-1);
  f0[3] = 4*l1*l2; // oppose au sommet 0
  f0[4] = 4*l0*l2; // oppose au sommet 1
  f0[5] = 4*l1*l0; // oppose au sommet 3
  }
 if(  whatd[op_dx] || whatd[op_dy] || whatd[op_dxx] || whatd[op_dyy] ||  whatd[op_dxy])
 {
   R2 Dl0(K.H(0)), Dl1(K.H(1)), Dl2(K.H(2));
  if (whatd[op_dx])
  {
    RN_ f0x(val('.',0,op_dx)); 
  f0x[0] = Dl0.x*l4_0;
  f0x[1] = Dl1.x*l4_1;
  f0x[2] = Dl2.x*l4_2;
  f0x[3] = 4*(Dl1.x*l2 + Dl2.x*l1) ;
  f0x[4] = 4*(Dl2.x*l0 + Dl0.x*l2) ;
  f0x[5] = 4*(Dl0.x*l1 + Dl1.x*l0) ;
  }

 if (whatd[op_dy])
  {  
    RN_ f0y(val('.',0,op_dy)); 
  f0y[0] = Dl0.y*l4_0;
  f0y[1] = Dl1.y*l4_1;
  f0y[2] = Dl2.y*l4_2;
  f0y[3] = 4*(Dl1.y*l2 + Dl2.y*l1) ;
  f0y[4] = 4*(Dl2.y*l0 + Dl0.y*l2) ;
  f0y[5] = 4*(Dl0.y*l1 + Dl1.y*l0) ;
  }
 
 if (whatd[op_dxx])
  {  
    RN_ fxx(val('.',0,op_dxx)); 

    fxx[0] = 4*Dl0.x*Dl0.x;
    fxx[1] = 4*Dl1.x*Dl1.x;
    fxx[2] = 4*Dl2.x*Dl2.x;
    fxx[3] =  8*Dl1.x*Dl2.x;
    fxx[4] =  8*Dl0.x*Dl2.x;
    fxx[5] =  8*Dl0.x*Dl1.x;
  }

 if (whatd[op_dyy])
  {  
    RN_ fyy(val('.',0,op_dyy)); 
    fyy[0] = 4*Dl0.y*Dl0.y;
    fyy[1] = 4*Dl1.y*Dl1.y;
    fyy[2] = 4*Dl2.y*Dl2.y;
    fyy[3] =  8*Dl1.y*Dl2.y;
    fyy[4] =  8*Dl0.y*Dl2.y;
    fyy[5] =  8*Dl0.y*Dl1.y;
  }
 if (whatd[op_dxy])
  {  
    assert(val.K()>op_dxy);
    RN_ fxy(val('.',0,op_dxy)); 
  
    fxy[0] = 4*Dl0.x*Dl0.y;
    fxy[1] = 4*Dl1.x*Dl1.y;
    fxy[2] = 4*Dl2.x*Dl2.y;
    fxy[3] =  4*(Dl1.x*Dl2.y + Dl1.y*Dl2.x);
    fxy[4] =  4*(Dl0.x*Dl2.y + Dl0.y*Dl2.x);
    fxy[5] =  4*(Dl0.x*Dl1.y + Dl0.y*Dl1.x);
  }
 
 }
 
}
/*
 void TypeOfFE_P1Lagrange::Pi_h(const baseFElement & K,RN_ & val, InterpolFunction f, R* v,int j,  void * arg) const
{
  const R2 Pt[] = { R2(0,0), R2(1,0), R2(0,1) };
   for (int i=0;i<3;i++)
     {  
     f(v,K.T(Pt[i]),K,i,Pt[i],arg),val[i]=*(v+j);}
 
}
 void TypeOfFE_P2Lagrange::Pi_h(const baseFElement & K,RN_ & val, InterpolFunction f, R* v,int j, void * arg) const
{
  const R2 Pt[] = { R2(0,0), R2(1,0), R2(0,1),R2(0.5,0.5),R2(0,0.5),R2(0.5,0) };
   for (int i=0;i<6;i++)
     {
     f(v,K.T(Pt[i]),K,i,Pt[i],arg),val[i]=*(v+j);}
 
}
*/
 
 
//TypeOfFE  P1Lagrange(1,0,0,P1Functions,D2_P1Functions,P1Interpolant,DataP1Lagrange);
//TypeOfFE  P2Lagrange(1,1,0,P2Functions,D2_P2Functions,P2Interpolant,DataP2Lagrange,3);

//  case of   fine mesh   
class TypeOfMortarCas1: public TypeOfMortar { 
  friend class FESpace;
  friend class FMortar;
  friend class ConstructDataFElement;
  protected:
  int NbLagrangeMult(const Mesh &,const Mortar &M) const ;
 
   int NbDoF(const Mesh &,const Mortar &M,int i) const 
     { int l(M.NbLeft()),r(M.NbRight());
       int n =Max(l,r);
       int mn=Min(l,r);
       return (l+r)*(NbDfOnVertex + NbDfOnEdge) + (n+1)*NbDfOnVertex + n*NbDfOnEdge -mn-1; 
      }
  int NbOfNodes(const Mesh &,const Mortar &M) const // call one time  
     {int l(M.NbLeft()),r(M.NbRight()); return (l+r)*(vertex_is_node+edge_is_node)+1;}
  int NbDoF(const Mesh &,const Mortar &M) const 
     { int l(M.NbLeft()),r(M.NbRight());
       int n =Max(l,r);
       int mn=Min(l,r);
       return (l+r)*(NbDfOnVertex + NbDfOnEdge) + (n+1)*NbDfOnVertex + n*NbDfOnEdge -mn-1; 
      }
  
   int NodeOfDF(const FESpace &Vh,const Mortar &M,int i) const 
     {throwassert(0);return 0;}
   int DFOfNode(const FESpace &Vh,const Mortar &M,int i) const 
     {throwassert(0);return 0;}
   TheSubFMortars  * Constructor(const FESpace &Vh,const Mortar &M);  
   void ConstructionOfNode(const Mesh &Th,int im,int * NodesOfElement,int *FirstNodeOfElement,int &lastnodenumber) const;
   void ConsTheSubMortar(FMortar & ) const; 
     
   const int vertex_is_node,edge_is_node;
  public: 
    TypeOfMortarCas1 (int i,int j): TypeOfMortar(i,j),
      vertex_is_node(i?1:0),edge_is_node(j?1:0) {};
     
}  MortarCas1P2(1,1) ;

 const TypeOfMortar & TheMortarCas1P2(MortarCas1P2); 
 

void TypeOfMortarCas1::ConstructionOfNode(const Mesh &Th,int im,int * NodesOfElement,int *FirstNodeOfElement,int &lastnodenumber) const
{  
  // im   mortar number 
 // trop complique on change 
             const Mortar &M(Th.mortars[im]);
             int k = Th.nt+im;
             int  kk=FirstNodeOfElement[k]; //  begin   
             // lagrange  multiplicator one new node 
              NodesOfElement[kk++] = lastnodenumber++;
/*                               
             int il = M.NbLeft();
             int ir = M.NbRight();
             int ir1 = ir-1;
             //  left
             
             for( int j=0;j<il;j++)  //  numbering vertex0 edge vertex1
              { 
                int K = M.TLeft(j);  // triangle
                int e = M.ELeft(j);  //  edge
                int nbneK = FirstNodeOfElement[K];
                int oe = vertex_is_node ? 3 : 0;
                int i0 = VerticesOfTriangularEdge[e][0];
                int i1 = VerticesOfTriangularEdge[e][1];
                if (vertex_is_node && !j)   //  just the first time 
                   NodesOfElement[kk++]=NodesOfElement[nbneK +i0];
                if (edge_is_node)
                   NodesOfElement[kk++]=NodesOfElement[nbneK+oe+e];
                if (vertex_is_node )  
                   NodesOfElement[kk++]=NodesOfElement[nbneK +i1];
              }

             //  right 
             for( int j=0;j<ir;j++)  //  numbering vertex0 edge vertex1
              { 
                int K = M.TRight(j);  // triangle
                int e = M.ERight(j);  //  edge
                int nbneK = FirstNodeOfElement[K];
                int oe = vertex_is_node ? 3 : 0;
                
                int i0 = VerticesOfTriangularEdge[e][1]; //  change the sens because right side
                int i1 = VerticesOfTriangularEdge[e][0];
              //  if (vertex_is_node &&  !j)   // never 
               //    NodesOfElement[kk++]=NodesOfElement[nbneK +i0];
                if (edge_is_node) 
                   NodesOfElement[kk++]=NodesOfElement[nbneK+oe+e];
                if (vertex_is_node  && (j != ir1) )  //  skip the last 
                   NodesOfElement[kk++]=NodesOfElement[nbneK +i1];
              } */
              
              throwassert(FirstNodeOfElement[k+1]==kk);
}

 R  d1P1F0(const FESpace *,const aSubFMortar *,R x) {return 1-x;}// 1 on 0
 R  d1P1F1 (const FESpace *,const aSubFMortar *,R x) {return x;}//  1 on 1
 
 R  d1P2F0 (const FESpace *,const aSubFMortar *,R x) {return (1-x)*(1-2*x);}// 1 on x=0
 R  d1P2F1(const FESpace *,const aSubFMortar *,R x) {return (1-x)*x*4;} // 1 on x=1/2
 R  d1P2F2(const FESpace *,const aSubFMortar *,R x) {return x*(2*x-1);} // 1 on x=1
 
 void  TypeOfMortarCas1::ConsTheSubMortar(FMortar & sm) const
   { //  constuction of 
   /* 
     int nbsm; // nb of submortar
  aSubFMortar * sm;
  ~FMortar() { delete [] dataDfNumberOFmul; delete [] dataf;}
  private:
  
  int *dataDfNumberOFmul;
   R (**dataf)(const FESpace *,const aSubFMortar *,R);

   */
   //  typedef
     const Mesh &Th(sm.Vh.Th);
     const Mortar & M(sm.M);
     int nl=M.NbLeft();
     int nr=M.NbRight();
     int nbsm= nl+nr-1;
     sm.nbsm = nbsm;
     int ldata = 6*nbsm;// 3 gauche+ 3 droite 
     sm.sm = new aSubFMortar[nbsm];
     sm.datai = new int [ldata];
     sm.dataf = new (R (*[ldata])(const FESpace *,const aSubFMortar *,R))  ;
     int *dataDfNumberOFmul=sm.datai;
     
     R (**dataf)(const FESpace *,const aSubFMortar *,R) ;
     dataf=sm.dataf;
     for (int i=0;i<ldata;i++) sm.dataf[i]=0;
   //  int * data0=sm.data+ldata;
   //  int * data1=data0+ldata;
     
     //  now the construction 
     int l=0,g=0;
     R2 A(M.VLeft(0));
     R2 B(M.VLeft(nl));
     throwassert(&M.VLeft(0) == &M.VRight(0));
     throwassert(&M.VLeft(nl) == &M.VRight(nr));
    
     R2 AB(A,B);
     R lg=Norme2(AB);
    // cout << " Mortar from " << A << " to " << B << " = " <<lg << endl;
     R2 AB1(AB/lg);
     int il=0,ir=0;
     int k=0;
     R la=0.0;
     R lb=0.0;
     R2 AA(A),BB(A);
   //  cout << "lg : " <<lg ;
     do {
       sm.sm[k].left  = M.left[il];
       sm.sm[k].right =  M.right[ir];
       R2 Bl(M.VLeft(il+1));
       R2 Br(M.VRight(ir+1));
       R ll=(AB1,Bl-A), lr=(AB1,Br-A);
     //  throwassert ( ll >=0 && lr >= 0);
     //  throwassert ( ll <=lg  && lr <= lg);
       
   //    cout << "AA , BB = " << AA << "," << BB << endl;
    //   cout << " " << ll << " " << lr << " ll=" << sm.sm[k].left << ", ";
       if (ll<lr) {BB=Bl,lb=ll,il++;} else {BB=Br, lb=lr,ir++;}
  //     cout << k << " " << k << " " << la/lg << " " << lb/lg << endl;
       sm.sm[k].a = la/lg;
       sm.sm[k].b = lb/lg;
       sm.sm[k].A=AA;
       sm.sm[k].B=BB;       
       la=lb;
       AA=BB;
       k++;
       throwassert(k<=nbsm);
     } while (il<nl && ir < nr);
     
  //   cout << "k=" << k <<endl;
     throwassert(nbsm==Max(nl,nr)); 
     //throwassert(nbsm<=k);
    nbsm=k;
    sm.nbsm=k;
//   construction of interpolation 
//  1) on a P1 on P2 
//   P2  si les longueurs des aSubMortar precedende et suivant  sont les meme 
//   sinon P1 
   //  calcul de leps 
     R leps=1.0/(1048576.0) ; // 1/2^20 
     R lgp=0;
     R lgc=0;
     R lgs=sm.sm[0].lg1();
    // cout << lgp << " " << lgc << " " << lgs << endl;
     
     int nmul=0;
     for (int k=0;k<nbsm;k++) 
       {
         
         lgp=lgc;
         lgc=lgs;
         lgs=  k+1 == nbsm  ? 0 : sm.sm[k+1].lg1();
         sm.sm[k].DfNumberOFmul= dataDfNumberOFmul;
         sm.sm[k].f=dataf;
        // cout << lgp << " " << lgc << " " << lgs << " ";
         if ( Abs(lgp-lgc) < leps && Abs(lgs-lgc) < leps )
          { // P2
           sm.sm[k].Nbmul=3;
           *dataDfNumberOFmul++=nmul++;
           *dataDfNumberOFmul++=nmul++;
           *dataDfNumberOFmul++=nmul;
           *dataf++ = d1P2F0;
           *dataf++ = d1P2F1;
           *dataf++ = d1P2F2;
          // cout << "P2 " << nmul << " " ;
           
          }
         else 
          { // P1
                   sm.sm[k].Nbmul=2;
           *dataDfNumberOFmul++=nmul++;
           *dataDfNumberOFmul++=nmul;
           *dataf++ = d1P1F0;
           *dataf++ = d1P1F1;
          // cout << "P1 " << nmul << " " ;
           

          }
       }
      nmul++;
     // cout << " " << nmul << " " <<  sm.NbDoF(0) << endl;
      throwassert(nmul==sm.NbDoF(0));
     
   }
   
  int TypeOfMortarCas1::NbLagrangeMult(const Mesh &,const Mortar &M) const 
     { 
       int nl = M.NbLeft();
       int nr = M.NbRight();
       R2 A(M.VLeft(0));
       R2 B(M.VLeft(nl));
       R2 AB(A,B);
       R lg=Norme2(AB);
       R leps = lg/1048576.0;
       throwassert(nl==1 || nr==1);
       R lgp=0,lgc=0,lgs=0;
       int nbmul=3; 
       if (nr==1) 
        {
        R2 AA(M.VLeft(0)),BB(M.VLeft(1));
        lgp= Norme2(BB-AA); // preced
        AA=BB;
        BB=M.VLeft(2); 
        lgc= Norme2(BB-AA); // courant 
        
        for (int i=1;i<nl-1;i++)
         { 
            AA=BB;
            BB=M.VLeft(i+2); 
            lgs=Norme2(AA-BB); // le suivant 
            if ( Abs(lgp-lgc) < leps && Abs(lgs-lgc) < leps )
              nbmul+=2; // P2
            else 
              nbmul+=1;// P1;
            lgp=lgc;
            lgc=lgs;
            
         }
        }
        else
        {
        R2 AA(M.VRight(0)),BB(M.VRight(1));
        lgp= Norme2(BB-AA); // preced
        AA=BB;
        BB=M.VRight(2); 
        lgc= Norme2(BB-AA); // courant 
        
        for (int i=1;i<nr-1;i++)
         { 
            AA=BB;
            BB=M.VRight(i+2); 
            lgs=Norme2(AA-BB); // le suivant 
            if ( Abs(lgp-lgc) < leps && Abs(lgs-lgc) < leps )
              nbmul+=2; // P2
            else 
              nbmul+=1;// P1;
            lgp=lgc;
            lgc=lgs;
            
         }
        }
       throwassert(nbmul>2);
       return nbmul;
      }  

 
// --- 
 FMortar::FMortar(const FESpace * VVh,int k)
  : 
    Vh(*VVh),
    M(Vh.Th.mortars[k-Vh.Th.nt]),
    N(VVh->N),
    p(Vh.PtrFirstNodeOfElement(k)),
    nbn(Vh.NbOfNodesInElement(k)),
    tom(Vh.tom)
    
 { throwassert(k>=Vh.Th.nt && k <Vh.Th.nt + Vh.Th.NbMortars);
   VVh->tom->ConsTheSubMortar(*this);}
   
 R TypeOfFE::operator()(const FElement & K,const  R2 & PHat,const KN_<R> & u,int componante,int op) const 
  {
   R v[1000],vf[100];
   assert(N*3*NbDoF<=1000 && NbDoF <100 );
   KNMK_<R> fb(v,NbDoF,N,op+1); //  the value for basic fonction
   KN_<R> fk(vf,NbDoF);
   for (int i=0;i<NbDoF;i++) // get the local value
    fk[i] = u[K(i)];
    //  get value of basic function
   bool whatd[last_operatortype];
   for (int i=0;i<last_operatortype;i++) 
     whatd[i]=false;
   whatd[op]=true;
   FB(whatd,K.Vh.Th,K.T,PHat,fb);  
   R r = (fb('.',componante,op),fk);  
   return r;
  }
 
static TypeOfFE_P1Lagrange P1LagrangeP1;
static TypeOfFE_P1Bubble P1BubbleP1;
static TypeOfFE_P2Lagrange P2LagrangeP2;

TypeOfFE  & P2Lagrange(P2LagrangeP2);
TypeOfFE  & P1Bubble(P1BubbleP1);
TypeOfFE  & P1Lagrange(P1LagrangeP1);

static ListOfTFE typefemP1("P1", &P1LagrangeP1);
static ListOfTFE typefemP1b("P1b", &P1BubbleP1);
static ListOfTFE typefemP2("P2", &P2LagrangeP2);
static  ListOfTFE typefemRT("RT0", &RTLagrange);
static  ListOfTFE typefemRTOrtho("RT0Ortho", &RTLagrangeOrtho);
 
 extern  TypeOfFE & RTmodifLagrange, & P1ttdc, & P2ttdc;
 static  ListOfTFE typefemRTmodif("RTmodif", &RTmodifLagrange);
 static ListOfTFE typefemP0("P0", &P0Lagrange);
 static ListOfTFE typefemP1nc("P1nc", &P1ncLagrange);
 static ListOfTFE typefemP1ttdc("P1dc", &P1ttdc);
 static ListOfTFE typefemP2ttdc("P2dc", &P2ttdc);
 
} // fin de namespace Fem2D 
