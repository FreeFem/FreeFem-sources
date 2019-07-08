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
#include "FESpacen.hpp"
#include "FESpace.hpp"
extern long verbosity ;
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
  return yy(SubArray(kk));

 }

ListOfTFE * ListOfTFE::all ; // list of all object of this type

void init_static_FE(); //   to correct so probleme with static Library FH aout 2004
//  the list of other FE file to force the link

ListOfTFE::ListOfTFE (const char * n,TypeOfFE *t) : name(n),tfe(t)
{
  if(!t)
  assert(t);
  static int count=0;
  if (count++==0)
    all=0; // init of all in dependant of the ordre of the objet file
  next=all;
  all=this;
 //  to correct so probleme with static Library FH aout 2004
 init_static_FE();
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
   const int k;
   const TypeOfFE ** teb;
   int nbn; // nb of node
   int  *  data;
   int  * data1;
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
     FESumConstruct(kk,Make(t,kk)),TypeOfFE(teb,kk,data,data1) {}
       TypeOfFESum(const TypeOfFE **t,int kk):
     FESumConstruct(kk,Make(t,kk)),TypeOfFE(teb,kk,data,data1) {}

   void FB(const bool * whatd,const Mesh & Th,const Triangle & K,const RdHat &PHat, RNMK_ & val) const;
   virtual void Pi_h_alpha(const baseFElement & K,KN_<double> & v) const
    {
      int k0=0;
      for (int i=0;i<k;i++) {
          int n=teb[i]->pij_alpha.N(); // ici BUG
          KN_<R> sv(v(SubArray(n,k0)));
          teb[i]->Pi_h_alpha(K,sv);
          k0+= n;}
      assert(pij_alpha.N()==k0);
    }
   ~TypeOfFESum(){  delete []  teb;}
} ;

class FEProduitConstruct { protected:
   int k;
   const TypeOfFE & teb;
   int * data;
   int * data1;
   FEProduitConstruct(int kk,const TypeOfFE &t)  ;
   ~FEProduitConstruct(){delete [] data;}
};

class TypeOfFEProduit: protected FEProduitConstruct, public  TypeOfFE { public:
  TypeOfFEProduit(int kk,const TypeOfFE &t):
    FEProduitConstruct(kk,t),TypeOfFE(t,kk,data,data1)  {}


  void FB(const bool * whatd,const Mesh & Th,const Triangle & K,const RdHat &PHat, RNMK_ & val) const;
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
  int No=teb.N;
  int N= No*kk;
  data = new int [n*(5+2)+3*N];
  data1 = data + n*(5)+N; // april 2006  add 2 array ????
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

   int ci=n;
   int cj=0;

  //  ou dans la partie miminal element finite atomic
  for (int i=0;i<m;i++)
   for (int j=0;j<kk;j++)
    {
      int il= teb.fromASubDF[i];
      int jl= teb.fromASubFE[i];
      data1[ci++]=il;
      data1[cj++]=j*teb.nb_sub_fem+jl;
    }
  //  warning the numbering of
  for(int j=0;j<kk;++j)
    for(int i=0;i<No;++i)
      data1[ci++]=0;//j*m+teb.begin_dfcomp[i];
  for(int j=0;j<kk;++j)
    for(int i=0;i<No;++i)
      data1[ci++]=m*kk;//j*m+teb.end_dfcomp[i];
  if (verbosity) cout << " kk "<< kk << " " << m << " : ";
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

  data = new int [n*(5+2) + 3*N];
  data1 = data + n*5+N; // april 2006  add 2 array ????

  int c=0;
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
  // error -- here
  //in case of [P2,P2],P1
   // we expect 0,0,1   and we get 0 1 2
   // => wrong BC ????
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


  //  ou dans la partie miminal element finite atomic

   int ci=n;
   int cf=2*n;
   int cl=cf+N;;
   int cj=0;
   int ccc=0;
   for ( j=0;j<kk;ccc+=teb[j++]->nb_sub_fem)
     for ( i=0;i<teb[j]->NbDoF;i++)
       {
	 int il= teb[j]->fromASubDF[i];
	 int jl= teb[j]->fromASubFE[i];
	 data1[ci++]=il;
	 data1[cj++]=ccc+jl;
       }

   for (int  j=0,ccn=0 ; j<kk ; ccn += teb[j++]->NbDoF)
     for(int k=0;k<teb[j]->N;++k)
       {
	 data1[cf++] = ccn + teb[j]->begin_dfcomp[k];
	 data1[cl++] = ccn + teb[j]->end_dfcomp[k];
       }
   ffassert(cl==2*n+2*N);

  ffassert(c== 5*n+N);
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
   void FB(const bool * whatd,const Mesh & Th,const Triangle & K,const RdHat &PHat, RNMK_ & val) const;

virtual R operator()(const FElement & K,const  R2 & PHat,const KN_<R> & u,int componante,int op) const ;

} ;
///////////////////////////////////////////////////////////////////////////////
// FH pour tester des idee de schema ----  Juin 2005 ---
///////////////////////////////////////////////////////////////////////////////

//  un VF cell centre
class TypeOfFE_P0VF : public  TypeOfFE { public:
  static int Data[];
  static double Pi_h_coef[];
   TypeOfFE_P0VF(): TypeOfFE(1,0,0,1,Data,1,1,3,3,Pi_h_coef)
    { const R2 Pt[] = { R2(0,0), R2(1,0), R2(0,1) };
      for (int i=0;i<NbDoF;i++) {
       pij_alpha[i]= IPJ(i,i,0);
       P_Pi_h[i]=Pt[i]; }
     }
   void FB(const bool * whatd,const Mesh & Th,const Triangle & K,const RdHat &PHat, RNMK_ & val) const;
   virtual R operator()(const FElement & K,const  RdHat & PHat,const KN_<R> & u,int componante,int op) const ;

} ;
int TypeOfFE_P0VF::Data[]={0,1,2,       0,0,0,       0,1,2,       0,0,0,        0,1,2,       0, 0,3};
double TypeOfFE_P0VF::Pi_h_coef[]={1.,1.,1.}; //  bofbof a verifier ...

 R TypeOfFE_P0VF::operator()(const FElement & K,const  R2 & PHat,const KN_<R> & u,int componante,int op) const
{
   R u0(u(K(0))), u1(u(K(1))), u2(u(K(2)));
   R r=0;
   if (op==0)
    {
      R l0=0,l1=PHat.x,l2=PHat.y;
      l1 = l1 * 3. < 1;
      l2 = l2 * 3. < 1;
      l0 = 1 - l0 -l2;

      r = u0*l0+u1*l1+l2*u2;
    }
   else
    {
      r =0;
    }
   return r;
}


void TypeOfFE_P0VF::FB(const bool *whatd,const Mesh & ,const Triangle & K,const RdHat & PHat,RNMK_ & val) const
{
//  const Triangle & K(FE.T);
  if (whatd[op_id])
   {
  R2 A(K[0]), B(K[1]),C(K[2]);
  R l0=1-PHat.x-PHat.y,l1=PHat.x,l2=PHat.y;
  l1 = l1 * 3. < 1;
  l2 = l2 * 3. < 1;
  l0 = 1 - l0 -l2;

  if (val.N() <3)
   throwassert(val.N() >=3);
  throwassert(val.M()==1 );

  val=0;
  RN_ f0(val('.',0,op_id));

    f0[0] = l0;
    f0[1] = l1;
    f0[2] = l2;}

}


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
   void FB(const bool * whatd,const Mesh & Th,const Triangle & K,const RdHat &PHat, RNMK_ & val) const;

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
   void FB(const bool * whatd,const Mesh & Th,const Triangle & K,const RdHat &PHat, RNMK_ & val) const;
} ;


class TypeOfFE_P2bLagrange : public  TypeOfFE { public:
  static int Data[];
  static double Pi_h_coef[];

   TypeOfFE_P2bLagrange(): TypeOfFE(1,1,1,1,Data,3,1,7,7,Pi_h_coef)
    { const R2 Pt[] = { R2(0,0), R2(1,0), R2(0,1),R2(0.5,0.5),R2(0,0.5),R2(0.5,0), R2(1./3.,1./3.) };
      for (int i=0;i<NbDoF;i++) {
       pij_alpha[i]= IPJ(i,i,0);
       P_Pi_h[i]=Pt[i]; }
     }

   void FB(const bool * whatd,const Mesh & Th,const Triangle & K,const RdHat &PHat, RNMK_ & val) const;
} ;

int TypeOfFE_P1Lagrange::Data[]={0,1,2,       0,0,0,       0,1,2,       0,0,0,        0,1,2,       0, 0,3};
int TypeOfFE_P1Bubble::Data[]={0,1,2,6,     0,0,0,0,     0,1,2,3,     0,0,0,0,        0,1,2,3,     0, 0,4};
int TypeOfFE_P2Lagrange::Data[]={0,1,2,3,4,5, 0,0,0,0,0,0, 0,1,2,3,4,5, 0,0,0,0,0,0,  0,1,2,3,4,5, 0 ,0,6};
int TypeOfFE_P2bLagrange::Data[]={0,1,2,3,4,5,6, 0,0,0,0,0,0,0, 0,1,2,3,4,5,6, 0,0,0,0,0,0,0,  0,1,2,3,4,5,6, 0,0,7};
double TypeOfFE_P1Lagrange::Pi_h_coef[]={1.,1.,1.};
double TypeOfFE_P1Bubble::Pi_h_coef[]={1.,1.,1.,1.};
double TypeOfFE_P2Lagrange::Pi_h_coef[]={1.,1.,1.,1.,1.,1.};
double TypeOfFE_P2bLagrange::Pi_h_coef[]={1.,1.,1.,1.,1.,1.,1.};

inline void dump(char *m,int n,int * p)
{
  cout << m ;
  for (int i=0;i<n;i++) cout << " " << p[i] ;
  cout << endl;
}



ConstructDataFElement::~ConstructDataFElement()
{
  if((*counter)--==0)
   {
    delete [] NodesOfElement;
    delete []  FirstNodeOfElement;
    delete [] FirstDfOfNode;
    (*counter)--; // correct bug oct 2008
    delete counter;
  }
}

 ConstructDataFElement::ConstructDataFElement(const ConstructDataFElement * t,int k)
  ://thecounter(0),
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
   (*counter)++;      // correction mai 2006 bug in counter incrementation

 }

ConstructDataFElement::ConstructDataFElement (const Mesh &Th,/*int NbDfOnSommet,int NbDfOnEdge,int NbDfOnElement*/
const  KN<const TypeOfFE *> & TFEs,const TypeOfMortar *tm,
int nbdfv,const int *ndfv,int nbdfe,const int *ndfe)
: counter(NewCounter())
{
 Make(Th,TFEs,/*NbDfOnSommet,NbDfOnEdge,NbDfOnElement,*/ tm,nbdfv,ndfv,nbdfe,ndfe);
}

ConstructDataFElement::ConstructDataFElement(const FESpace ** l,int k,const KN<const TypeOfFE *>  & TFEs)
: counter(NewCounter())
{
 const Mesh & Th(l[0]->Th);
 for  (int i=0;i<k;i++)
   {
     ffassert( &Th== &l[i]->Th);
     ffassert( l[i]->TFE.constant());
   }

 Make(Th,TFEs);//NbDfOnSommet,NbDfOnEdge,NbDfOnElement,0);
}


void ConstructDataFElement::Make(const Mesh &Th,
const  KN<const TypeOfFE *> & TFEs,

/*int NbDfOnSommet,int NbDfOnEdge,int NbDfOnElement,*/const TypeOfMortar *tm,
int nb_dfv,const int *ndfv,int nb_dfe,const int *ndfe)
/* pour le condition de periodicit� ..
  nb_dfv :  nombre de sommets recolle par periodicit�
  ndfv:  numerotation des sommets pour CL de periodicite
  ndfv[i] = numero du sommet i (i<Th.nv)  recolle par periodicit�
  nb_dfe:  nombre de arete frontiere  recolle par periodicit�
  ndfe[i]:  numero de l'arete frontiere (i < Th.neb) recoll�.
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
  ffassert( tm || Th.NbMortars==0);
  NbOfElements = Th.nt; //  by default
  // if mortars
  int NbOfNodeL=0;
  NbOfElements += Th.NbMortars;
  FirstDfOfNode =0;
  FirstNodeOfElement=0;
  MaxNbDFPerElement=nbdfe;
  assert(3*NbDfOnSommet+3*NbDfOnEdge+NbDfOnElement==MaxNbDFPerElement);

  int ks=TFE.NbNodeOnVertex>0,
      ke=TFE.NbNodeOnEdge>0,
      kt=TFE.NbNodeOnElement>0;
    MaxNbNodePerElement=nbne;

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
          ffassert(tm);
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
          int j,k=Th.BoundaryElement(be,j);
          int jj=j;
          int kk=Th.ElementAdj(k,jj);
          for (int kkk=0;kkk<NbNodeonEdge;kkk++)
           {
             int kj=kkk,kjj=-1;
            // assert(kk != k);
             if ( kk >= 0 && jj >=0  &&  !(( kk == k ) && ( jj=j ) ) )
             {
              if (k < kk ) kjj = NbNodeonEdge-kj-1; //
              else   kj = NbNodeonEdge-kj-1, kjj = NbNodeonEdge-kj-1;
             }
          if (kjj >=0)
            NodesOfElement[kk*nbne+oe+jj] = nn + ndfe[be]*NbNodeonEdge+ kjj   ; // adj
          NodesOfElement[k*nbne+oe+j]   = nn + ndfe[be]*NbNodeonEdge+ kj   ; // new
          }
         }
         nn += nb_dfe;
        }
       for (int k=0;k<Th.nt;k++)
        {
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
               int kk=Th.ElementAdj(k,jj);
               int jjj=jj*NbNodeonEdge+NbNodeonEdge-ll-1; // autre sens
               assert(kk>=k && jjj == jj);
               NodesOfElement[kk*nbne+oe+jjj] = nn   ; // adj
               NodesOfElement[i++]           = nn++ ; // new
             }  else i++;
          }
          for (int jj=0;jj<NbNodeonElement;jj++)
             NodesOfElement[i++]= nn++;
        }
        firstmul=nn;
        if (Th.NbMortars)
          {
            //  construction of the mortars element
           int * color= new int [firstmul];
           int thecolor=0;
           for (int j=0;j<firstmul;j++)
             color[j]=thecolor;

           for (int im=0;im<Th.NbMortars;im++)
            {
             thecolor++; // get a new color
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
                    throwassert(node<firstmul);
                    if (color[node] != thecolor) //  new node => adding
                     {
                       color[node] = thecolor;
                       NodesOfElement[i++] = node;
                     }
                 }
               }
               ffassert(i==FirstNodeOfElement[im+Th.nt+1]);
              }
             delete [] color;
          }
        else
          ffassert(i==Th.nt*nbne && i );
        NbOfNode=nn;


   }
  NbOfNode=nn;
  int NbOfDFL=0;
   if (! samendf)
       {

         ffassert(NodesOfElement);
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
   if(verbosity>2)
   {
     cout << "  FESpace: Nb Of Nodes = " << nn ;
     if(NbOfNodeL)
       cout << " Nb of Lagrange Mul Node = " << NbOfNodeL  ;
     cout << " Nb of DF = " << NbOfDF << endl;
     if(NbOfDFL) {
       cout << " Nb of Lagrange Mul DF = "   << NbOfDFL ;
       cout << " MaxNbDFPerElement     =   " << MaxNbDFPerElement ;
     };
       cout << endl;
   }
}

FESpace::FESpace(const FESpace & Vh,int k )
 :
     Th(Vh.Th),
     ptrTFE(new TypeOfFEProduit(k,*Vh.TFE[0])),
     TFE(1,0,ptrTFE),
     cdef(Vh.cdef?new ConstructDataFElement(Vh.cdef,k):0),
     cmesh(Vh.Th),
     N(Vh.N*k),
     Nproduit(Vh.Nproduit*k),
     NbOfDF(Vh.NbOfDF*k),
     NbOfElements(Vh.NbOfElements),
     NbOfNodes(Vh.NbOfNodes),
     nb_sub_fem(TFE[0]->nb_sub_fem),
     dim_which_sub_fem(TFE[0]->dim_which_sub_fem),
     NodesOfElement(Vh.NodesOfElement),
     FirstNodeOfElement(Vh.FirstNodeOfElement),
     FirstDfOfNodeData(cdef?cdef->FirstDfOfNode:0),
     tom(0),
     MaxNbNodePerElement(Vh.MaxNbNodePerElement),
     MaxNbDFPerElement(Vh.MaxNbDFPerElement*k)
{
      // correction mai 2006 no renumbering of existing cdef
     if(cdef && (Vh.cdef && Vh.cdef->counter != cdef->counter)) {
       renum(); // correction mai 2006 no renumbering of existing cdef
       }
    Show();
     }

FESpace::FESpace(const FESpace ** Vh,int k )
 :
     Th((**Vh).Th),
     ptrTFE(new TypeOfFESum(Vh,k)),
     TFE(1,0,ptrTFE),
     cdef(new ConstructDataFElement(Vh,k,TFE)),
     cmesh((**Vh).Th),
     N(sum(Vh,&FESpace::N,k)),
     Nproduit(cdef->Nproduit),
     NbOfDF(cdef->NbOfDF),
     NbOfElements(cdef->NbOfElements),
     NbOfNodes(cdef->NbOfNode),
     nb_sub_fem(TFE[0]->nb_sub_fem),
     dim_which_sub_fem(TFE[0]->dim_which_sub_fem),
     NodesOfElement(cdef->NodesOfElement),
     FirstNodeOfElement(cdef->FirstNodeOfElement),
     FirstDfOfNodeData(cdef->FirstDfOfNode),
     tom(0) ,
     MaxNbNodePerElement(cdef->MaxNbNodePerElement),
     MaxNbDFPerElement(cdef->MaxNbDFPerElement)
{
     if(cdef) renum();
    Show();
    // verification
    long snbdf=0;
    for(int i=0;i<k;++i)
        snbdf += Vh[i]->NbOfDF;
    if( snbdf !=NbOfDF)
        cerr << " Problem build of FEspace (2d) (may be : due to periodic Boundary condition missing ) FH " << endl
             << " The number of DF must be " << snbdf << "  and it is " << NbOfDF <<endl;
    ffassert(snbdf == NbOfDF );
}

FESpace::FESpace(const Mesh & TTh,const TypeOfFE ** tef,int k,int nbdfv,const int *ndfv,int nbdfe,const int *ndfe )
 :
     Th(TTh),
     ptrTFE(new TypeOfFESum(tef,k)),
     TFE(1,0,ptrTFE),
     cdef(new ConstructDataFElement(TTh,TFE,//sum(tef,&TypeOfFE::NbDfOnVertex,k),
                                       // sum(tef,&TypeOfFE::NbDfOnEdge,k),
                                        //sum(tef,&TypeOfFE::NbDfOnElement,k),
                                        0,nbdfv,ndfv,nbdfe,ndfe)),
     cmesh(TTh),
     N(sum(tef,&TypeOfFE::N,k)),
     Nproduit(cdef->Nproduit),

     NbOfDF(cdef->NbOfDF),
     NbOfElements(cdef->NbOfElements),
     NbOfNodes(cdef->NbOfNode),
     nb_sub_fem(TFE[0]->nb_sub_fem),
     dim_which_sub_fem(TFE[0]->dim_which_sub_fem),
     NodesOfElement(cdef->NodesOfElement),
     FirstNodeOfElement(cdef->FirstNodeOfElement),
     FirstDfOfNodeData(cdef->FirstDfOfNode),
     tom(0) ,
     MaxNbNodePerElement(cdef->MaxNbNodePerElement),
     MaxNbDFPerElement(cdef->MaxNbDFPerElement)
{
  if(cdef) renum();
  Show();
}


 FESpace::FESpace(const Mesh & TTh,const TypeOfFE & tef,int nbdfv,const int *ndfv,int nbdfe,const int *ndfe)
   :
     Th(TTh),
     ptrTFE(0),
     TFE(1,0,&tef),
     cdef(new ConstructDataFElement(TTh,TFE,0,nbdfv,ndfv,nbdfe,ndfe)),
     cmesh(TTh),
     N(tef.N),
     Nproduit(cdef->Nproduit),
     NbOfDF(cdef->NbOfDF),
     NbOfElements(cdef->NbOfElements),
     NbOfNodes(cdef->NbOfNode),
     nb_sub_fem(TFE[0]->nb_sub_fem),
     dim_which_sub_fem(TFE[0]->dim_which_sub_fem),
     NodesOfElement(cdef->NodesOfElement),
     FirstNodeOfElement(cdef->FirstNodeOfElement),
     FirstDfOfNodeData(cdef->FirstDfOfNode),
     tom(0),
     MaxNbNodePerElement(cdef->MaxNbNodePerElement),
     MaxNbDFPerElement(cdef->MaxNbDFPerElement)
{
  if(tef.NbDfOnVertex || tef.NbDfOnEdge) renum();
  Show();
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
     Th(TTh),
     ptrTFE(0),
     TFE(1,0,&tef),
     cdef(new ConstructDataFElement(TTh,TFE,&tm)),//tef.NbDfOnVertex,tef.NbDfOnEdge,tef.NbDfOnElement,&tm)),
     cmesh(TTh),
     N(tef.N),
     Nproduit(1),
     NbOfDF(cdef->NbOfDF),
     NbOfElements(cdef->NbOfElements),
     NbOfNodes(cdef->NbOfNode),
     nb_sub_fem(TFE[0]->nb_sub_fem),
     dim_which_sub_fem(TFE[0]->dim_which_sub_fem),
     NodesOfElement(cdef->NodesOfElement),
     FirstNodeOfElement(cdef->FirstNodeOfElement),
     FirstDfOfNodeData(cdef->FirstDfOfNode),
     tom(&tm),
     MaxNbNodePerElement(cdef->MaxNbNodePerElement),
     MaxNbDFPerElement(cdef->MaxNbDFPerElement)
{
       renum();
    Show();
     }

void ConstructDataFElement::renum(const long *r,int l)
 {
   ffassert(this);
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

 void TypeOfFEProduit::FB(const bool * whatd,const Mesh & Th,const Triangle & K,const RdHat & PHat,RNMK_ & val) const
 {
   int n=teb.NbDoF;
   int m=teb.N;
   val=0.0;
   SubArray t(val.K());
   RNMK_ v(val(SubArray(n,0,k),SubArray(m),t));
   teb.FB(whatd,Th,K,PHat,v);
   for (int i=1;i<k;i++)
     val(SubArray(n,i,k),SubArray(m,m*i),t)=v;
 }

  void TypeOfFESum::FB(const bool * whatd,const Mesh & Th,const Triangle & K,const RdHat & PHat,RNMK_ & val) const
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
       teb[i]->FB(whatd,Th,K,PHat,v);
     else
       v=val(SubArray(DF[j+1]-DF[j],DF[j]),SubArray(NN[j+1]-NN[j],NN[j]),t);
    }
 }
 R TypeOfFE_P1Lagrange::operator()(const FElement & K,const  R2 & PHat,const KN_<R> & u,int componante,int op) const
{
   R u0(u(K(0))), u1(u(K(1))), u2(u(K(2)));
   R r=0;
   if (op==op_id)
    {
      R l0=1-PHat.x-PHat.y,l1=PHat.x,l2=PHat.y;
      r = u0*l0+u1*l1+l2*u2;
    }
   else
    {
       const Triangle & T=K.T;
       R2 D0 = T.H(0) , D1 = T.H(1)  , D2 = T.H(2) ;
       if (op==op_dx)
         r =  D0.x*u0 + D1.x*u1 + D2.x*u2 ;
        else if(op==op_dy)
         r =  D0.y*u0 + D1.y*u1 + D2.y*u2 ;
    }
   return r;
}


void TypeOfFE_P1Lagrange::FB(const bool *whatd,const Mesh & ,const Triangle & K,const RdHat & PHat,RNMK_ & val) const
{
  R2 A(K[0]), B(K[1]),C(K[2]);
  R l0=1-PHat.x-PHat.y,l1=PHat.x,l2=PHat.y;

  if (val.N() <3)
   throwassert(val.N() >=3);
  throwassert(val.M()==1 );

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


void TypeOfFE_P1Bubble::FB(const bool *whatd,const Mesh & ,const Triangle & K,const RdHat & PHat,RNMK_ & val) const
{
  R2 A(K[0]), B(K[1]),C(K[2]);
  R l0=1-PHat.x-PHat.y, l1=PHat.x, l2=PHat.y, lb=l0*l1*l2*9.;

  if (val.N() <4)
   throwassert(val.N() >=4);
  throwassert(val.M()==1 );

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
    R lbdxx= 18*((Dl0.x*Dl1.x)*l2+(Dl1.x*Dl2.x)*l0+(Dl2.x*Dl0.x)*l1);
    fxx[0] = -lbdxx;
    fxx[1] = -lbdxx;
    fxx[2] = -lbdxx;
    fxx[3] = 3*lbdxx;
  }

 if (whatd[op_dyy])
  {
    RN_ fyy(val('.',0,op_dyy));
    R lbdyy= 18*((Dl0.y*Dl1.y)*l2+(Dl1.y*Dl2.y)*l0+(Dl2.y*Dl0.y)*l1);

    fyy[0] =  -lbdyy;
    fyy[1] =  -lbdyy;
    fyy[2] =  -lbdyy;
    fyy[3] =  3*lbdyy;
  }
 if (whatd[op_dxy])
  {
    assert(val.K()>op_dxy);
    RN_ fxy(val('.',0,op_dxy));
    R lbdxy= 9*(Dl0.x*Dl1.y+ Dl0.y*Dl1.x)*l2+(Dl1.x*Dl2.y+Dl1.y*Dl2.x)*l0+(Dl2.x*Dl0.y+Dl2.y*Dl0.x)*l1;
    fxy[0] = 4*Dl0.x*Dl0.y-lbdxy;
    fxy[1] = 4*Dl1.x*Dl1.y-lbdxy;
    fxy[2] = 4*Dl2.x*Dl2.y-lbdxy;
    fxy[3] = +3*lbdxy;
  }
  }

}
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

void TypeOfFE_P2Lagrange::FB(const bool *whatd,const Mesh & ,const Triangle & K,const RdHat & PHat,RNMK_ & val) const
{
  R2 A(K[0]), B(K[1]),C(K[2]);
  R l0=1-PHat.x-PHat.y,l1=PHat.x,l2=PHat.y;
  R l4_0=(4*l0-1),l4_1=(4*l1-1),l4_2=(4*l2-1);

  throwassert( val.N()>=6);
  throwassert(val.M()==1);

  val=0;
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


void TypeOfFE_P2bLagrange::FB(const bool *whatd,const Mesh & ,const Triangle & K,const RdHat & PHat,RNMK_ & val) const
{
  R2 A(K[0]), B(K[1]),C(K[2]);
  R l0=1-PHat.x-PHat.y,l1=PHat.x,l2=PHat.y,lb=l0*l1*l2*3.;
  R l4_0=(4*l0-1),l4_1=(4*l1-1),l4_2=(4*l2-1);

  throwassert( val.N()>=7);
  throwassert(val.M()==1);

  val=0;
 if (whatd[op_id])
  {
   R lb4=lb*4;
   RN_ f0(val('.',0,op_id));
  f0[0] = l0*(2*l0-1)+lb;
  f0[1] = l1*(2*l1-1)+lb;
  f0[2] = l2*(2*l2-1)+lb;
  f0[3] = 4*l1*l2-lb4; // oppose au sommet 0
  f0[4] = 4*l0*l2-lb4; // oppose au sommet 1
  f0[5] = 4*l1*l0-lb4; // oppose au sommet 3
  f0[6] = 9*lb;

  }
 if(  whatd[op_dx] || whatd[op_dy] || whatd[op_dxx] || whatd[op_dyy] ||  whatd[op_dxy])
 {
   R2 Dl0(K.H(0)), Dl1(K.H(1)), Dl2(K.H(2)),     Dlb((Dl0*l1*l2+Dl1*l0*l2+Dl2*l0*l1)*3.),
    Dlb4(Dlb*4.);

  if (whatd[op_dx])
  {
    RN_ f0x(val('.',0,op_dx));
  f0x[0] = Dl0.x*l4_0 +Dlb.x;
  f0x[1] = Dl1.x*l4_1 +Dlb.x;
  f0x[2] = Dl2.x*l4_2 +Dlb.x;
  f0x[3] = 4*(Dl1.x*l2 + Dl2.x*l1) -Dlb4.x;
  f0x[4] = 4*(Dl2.x*l0 + Dl0.x*l2) -Dlb4.x;
  f0x[5] = 4*(Dl0.x*l1 + Dl1.x*l0) -Dlb4.x;
  f0x[6] = 9.*Dlb.x;
  }

 if (whatd[op_dy])
  {
   RN_ f0y(val('.',0,op_dy));
  f0y[0] = Dl0.y*l4_0  +Dlb.y;
  f0y[1] = Dl1.y*l4_1  +Dlb.y;
  f0y[2] = Dl2.y*l4_2  +Dlb.y;
  f0y[3] = 4*(Dl1.y*l2 + Dl2.y*l1)  -Dlb4.y;
  f0y[4] = 4*(Dl2.y*l0 + Dl0.y*l2)  -Dlb4.y;
  f0y[5] = 4*(Dl0.y*l1 + Dl1.y*l0)  -Dlb4.y;
  f0y[6] = 9*Dlb.y;

  }

 if (whatd[op_dxx])
  {
    RN_ fxx(val('.',0,op_dxx));
    R lbdxx= 6*((Dl0.x*Dl1.x)*l2+(Dl1.x*Dl2.x)*l0+(Dl2.x*Dl0.x)*l1);
    R lbd4xx= 4*lbdxx;

    fxx[0] = 4*Dl0.x*Dl0.x  + lbdxx;
    fxx[1] = 4*Dl1.x*Dl1.x  + lbdxx;;
    fxx[2] = 4*Dl2.x*Dl2.x  + lbdxx;;
    fxx[3] = 8*Dl1.x*Dl2.x  - lbd4xx;
    fxx[4] = 8*Dl0.x*Dl2.x  - lbd4xx;
    fxx[5] = 8*Dl0.x*Dl1.x  - lbd4xx;
    fxx[6] = 9*lbdxx;
      }

 if (whatd[op_dyy])
  {
    RN_ fyy(val('.',0,op_dyy));
    R lbdyy= 6*((Dl0.y*Dl1.y)*l2+(Dl1.y*Dl2.y)*l0+(Dl2.y*Dl0.y)*l1);
    R lbd4yy= lbdyy*4;

    fyy[0] = 4*Dl0.y*Dl0.y + lbdyy;
    fyy[1] = 4*Dl1.y*Dl1.y + lbdyy;
    fyy[2] = 4*Dl2.y*Dl2.y + lbdyy;
    fyy[3] = 8*Dl1.y*Dl2.y - lbd4yy;
    fyy[4] = 8*Dl0.y*Dl2.y - lbd4yy;
    fyy[5] = 8*Dl0.y*Dl1.y - lbd4yy;
    fyy[6] = 9*lbdyy;

  }
 if (whatd[op_dxy])
  {
    assert(val.K()>op_dxy);
    RN_ fxy(val('.',0,op_dxy));
    R lbdxy= 3*(Dl0.x*Dl1.y+ Dl0.y*Dl1.x)*l2+(Dl1.x*Dl2.y+Dl1.y*Dl2.x)*l0+(Dl2.x*Dl0.y+Dl2.y*Dl0.x)*l1;
    R lbd4xy= lbdxy*4;

    fxy[0] = 4*Dl0.x*Dl0.y + lbdxy;
    fxy[1] = 4*Dl1.x*Dl1.y + lbdxy;
    fxy[2] = 4*Dl2.x*Dl2.y + lbdxy;
    fxy[3] =  4*(Dl1.x*Dl2.y + Dl1.y*Dl2.x) - lbd4xy;
    fxy[4] =  4*(Dl0.x*Dl2.y + Dl0.y*Dl2.x) - lbd4xy;
    fxy[5] =  4*(Dl0.x*Dl1.y + Dl0.y*Dl1.x) - lbd4xy;
    fxy[6] =  9.*lbdxy;
  }

 }

}

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
     {ffassert(0);return 0;}
   int DFOfNode(const FESpace &Vh,const Mortar &M,int i) const
     {ffassert(0);return 0;}
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
  //  const Mortar &M(Th.mortars[im]);
             int k = Th.nt+im;
             int  kk=FirstNodeOfElement[k]; //  begin
             // lagrange  multiplicator one new node
              NodesOfElement[kk++] = lastnodenumber++;
              throwassert(FirstNodeOfElement[k+1]==kk);
}

 R  d1P1F0(const FESpace *,const aSubFMortar *,R x) {return 1-x;}// 1 on 0
 R  d1P1F1 (const FESpace *,const aSubFMortar *,R x) {return x;}//  1 on 1

 R  d1P2F0 (const FESpace *,const aSubFMortar *,R x) {return (1-x)*(1-2*x);}// 1 on x=0
 R  d1P2F1(const FESpace *,const aSubFMortar *,R x) {return (1-x)*x*4;} // 1 on x=1/2
 R  d1P2F2(const FESpace *,const aSubFMortar *,R x) {return x*(2*x-1);} // 1 on x=1

 void  TypeOfMortarCas1::ConsTheSubMortar(FMortar & sm) const
   { //  constuction of
     const Mortar & M(sm.M);
     int nl=M.NbLeft();
     int nr=M.NbRight();
     int nbsm= nl+nr-1;
     sm.nbsm = nbsm;
     int ldata = 6*nbsm;// 3 gauche+ 3 droite
     sm.sm = new aSubFMortar[nbsm];
     sm.datai = new int [ldata];
     typedef  R (* Fdataf)(const FESpace *,const aSubFMortar *,R);
     sm.dataf =new Fdataf[ldata]; //  new (R (*[ldata])(const FESpace *,const aSubFMortar *,R))  ;
     ffassert( sm.dataf ); //  remove previous line FH, PB comp
     int *dataDfNumberOFmul=sm.datai;

     R (**dataf)(const FESpace *,const aSubFMortar *,R) ;
     dataf=sm.dataf;
     for (int i=0;i<ldata;i++) sm.dataf[i]=0;
     R2 A(M.VLeft(0));
     R2 B(M.VLeft(nl));
     ffassert(&M.VLeft(0) == &M.VRight(0));
     ffassert(&M.VLeft(nl) == &M.VRight(nr));

     R2 AB(A,B);
     R lg=Norme2(AB);
     R2 AB1(AB/lg);
     int il=0,ir=0;
     int k=0;
     R la=0.0;
     R lb=0.0;
     R2 AA(A),BB(A);
     do {
       sm.sm[k].left  = M.left[il];
       sm.sm[k].right =  M.right[ir];
       R2 Bl(M.VLeft(il+1));
       R2 Br(M.VRight(ir+1));
       R ll=(AB1,Bl-A), lr=(AB1,Br-A);
       if (ll<lr) {BB=Bl,lb=ll,il++;} else {BB=Br, lb=lr,ir++;}
       sm.sm[k].a = la/lg;
       sm.sm[k].b = lb/lg;
       sm.sm[k].A=AA;
       sm.sm[k].B=BB;
       la=lb;
       AA=BB;
       k++;
       throwassert(k<=nbsm);
     } while (il<nl && ir < nr);

     ffassert(nbsm==Max(nl,nr));
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

     int nmul=0;
     for (int k=0;k<nbsm;k++)
       {

         lgp=lgc;
         lgc=lgs;
         lgs=  k+1 == nbsm  ? 0 : sm.sm[k+1].lg1();
         sm.sm[k].DfNumberOFmul= dataDfNumberOFmul;
         sm.sm[k].f=dataf;
         if ( Abs(lgp-lgc) < leps && Abs(lgs-lgc) < leps )
          { // P2
           sm.sm[k].Nbmul=3;
           *dataDfNumberOFmul++=nmul++;
           *dataDfNumberOFmul++=nmul++;
           *dataDfNumberOFmul++=nmul;
           *dataf++ = d1P2F0;
           *dataf++ = d1P2F1;
           *dataf++ = d1P2F2;

          }
         else
          { // P1
                   sm.sm[k].Nbmul=2;
           *dataDfNumberOFmul++=nmul++;
           *dataDfNumberOFmul++=nmul;
           *dataf++ = d1P1F0;
           *dataf++ = d1P1F1;


          }
       }
      nmul++;
      ffassert(nmul==sm.NbDoF(0));

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
       ffassert(nl==1 || nr==1);
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
       ffassert(nbmul>2);
       return nbmul;
      }


// ---
 FMortar::FMortar(const FESpace * VVh,int k)
  :
    Vh(*VVh),
    M(Vh.Th.mortars[k-Vh.Th.nt]),
    p(Vh.PtrFirstNodeOfElement(k)),
    nbn(Vh.NbOfNodesInElement(k)),
    N(VVh->N),
    tom(Vh.tom)

 { ffassert(k>=Vh.Th.nt && k <Vh.Th.nt + Vh.Th.NbMortars);
   VVh->tom->ConsTheSubMortar(*this);}

 R TypeOfFE::operator()(const FElement & K,const  R2 & PHat,const KN_<R> & u,int componante,int op) const
  {
   R v[10000],vf[1000];
   assert(N*3*NbDoF<=10000 && NbDoF <1000 );
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
static TypeOfFE_P0VF VFP0VF;
static TypeOfFE_P1Bubble P1BubbleP1;
static TypeOfFE_P2Lagrange P2LagrangeP2;
static TypeOfFE_P2bLagrange P2bLagrangeP2;

TypeOfFE  & P2Lagrange(P2LagrangeP2);
TypeOfFE  & P2bLagrange(P2bLagrangeP2);
TypeOfFE  & P1Bubble(P1BubbleP1);
TypeOfFE  & P1Lagrange(P1LagrangeP1);
TypeOfFE  & P0VF(VFP0VF);

static ListOfTFE typefemP1("P1", &P1LagrangeP1);
static ListOfTFE typefemP0VF("P0VF", &P0VF);  //
static ListOfTFE typefemP1b("P1b", &P1BubbleP1);
static ListOfTFE typefemP2("P2", &P2LagrangeP2);
static  ListOfTFE typefemRT("RT0", &RTLagrange);
static  ListOfTFE typefemRTOrtho("RT0Ortho", &RTLagrangeOrtho);

 extern  TypeOfFE & RTmodifLagrange, & P1ttdc, & P2ttdc,   & P0edge;
 static  ListOfTFE typefemRTmodif("RTmodif", &RTmodifLagrange);
 static ListOfTFE typefemP0("P0", &P0Lagrange);
 static ListOfTFE typefemP1nc("P1nc", &P1ncLagrange);
 static ListOfTFE typefemP1ttdc("P1dc", &P1ttdc);
 static ListOfTFE typefemP2ttdc("P2dc", &P2ttdc);
 static ListOfTFE typefemP2b("P2b", &P2bLagrangeP2);

 static ListOfTFE typefemP0edge("P0edge", &P0edge);

// correct Probleme of static library link with new make file
void init_static_FE()
{ //  list of other FE file.o
   extern void init_FE_P2h() ;
  init_FE_P2h() ;
}

} // fin de namespace Fem2D
