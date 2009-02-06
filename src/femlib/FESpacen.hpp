// ********** DO NOT REMOVE THIS BANNER **********
// ORIG-DATE:     Jan 2008
// -*- Mode : c++ -*-
//
// SUMMARY  : Generic Fiinite Element header  1d, 2d, 3d  
// USAGE    : LGPL      
// ORG      : LJLL Universite Pierre et Marie Curi, Paris,  FRANCE 
// AUTHOR   : Frederic Hecht
// E-MAIL   : frederic.hecht@ann.jussieu.fr
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

 Thank to the ARN ()  FF2A3 grant
 ref:ANR-07-CIS7-002-01 
 */

#ifndef FESpacen_HPP_
#define FESpacen_HPP_
/*
 *  FEspacen.hpp
 *  EF23n
 *
 *  Created by Frédéric Hecht on 04/12/07.
 *  Copyright 2007 Universite Pierre et marie Curie  All rights reserved.
 *
 */
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <complex>
#include <map>
using namespace std;
#include "error.hpp"
#include "ufunction.hpp"
#include "Mesh3dn.hpp"
#include "Mesh2dn.hpp"
#include "Mesh1dn.hpp"

#include "RNM.hpp"


#include "QuadratureFormular.hpp"

namespace Fem2D {

template<class Mesh> class GFESpace;
template<class Mesh> class GFElement;
template<class Mesh> class GbaseFElement;
template<class Mesh> class GTypeOfFE;
 
// numbering of derivative 
enum operatortype { 
    op_id=0,
    op_dx=1,op_dy=2,
    op_dxx=3,op_dyy=4,
    op_dyx=5,op_dxy=5,   
    op_dz=6,
    op_dzz=7,     
    op_dzx=8,op_dxz=8, 
    op_dzy=9,op_dyz=9   
}; 

typedef unsigned int What_d;
 
const unsigned int Fop_id= 1<< op_id;

const unsigned int Fop_dx= 1<< op_dx;
const unsigned int Fop_dy= 1<< op_dy;
const unsigned int Fop_dz= 1<< op_dz;

const unsigned int Fop_dxx= 1<< op_dxx;
const unsigned int Fop_dxy= 1<< op_dxy;
const unsigned int Fop_dxz= 1<< op_dxz;

const unsigned int Fop_dyx= 1<< op_dyx;
const unsigned int Fop_dyy= 1<< op_dyy;
const unsigned int Fop_dyz= 1<< op_dyz;

const unsigned int Fop_dzx= 1<< op_dzx;
const unsigned int Fop_dzy= 1<< op_dzy;
const unsigned int Fop_dzz= 1<< op_dzz;

const unsigned int Fop_D0 = Fop_id;
const unsigned int Fop_D1 = Fop_dx | Fop_dy | Fop_dz;
const unsigned int Fop_D2 = Fop_dxx | Fop_dyy | Fop_dzz | Fop_dxy | Fop_dxz | Fop_dyz;
const unsigned int Fop_Dall = Fop_D0| Fop_D1| Fop_D2;

inline What_d  Fwhatd(const operatortype op) { return 1<< op;}


const int last_operatortype=10;
const bool operatortypeValue[last_operatortype]= {true,false,false,false,false,false,false,false,false,false} ; 


inline void initwhatd(bool *whatd,int k)
{
    for (int i=0;i<last_operatortype;i++)
	whatd[i]=false;    
    whatd[k]=true;
}

typedef double R;
typedef KN_<R> RN_;
typedef KN<R>  RN;
typedef KNM_<R> RNM_;
typedef KNMK_<R> RNMK_;
typedef KNMK<R>  RNMK;

template<class Mesh> class GFElement; 
template<class Mesh> class GFESpace; 
template<class Mesh>class GTypeOfFE ;

class dataTypeOfFE {
    

private:
  const int * data; 
  const int * dataalloc; 
    
public:
    const int ndfonVertex; 
    const int ndfonEdge;
    const int ndfonFace;
    const int ndfonVolume;
     const int NbDoF;
    const int NbNode;
    int  N,nb_sub_fem;  
    const int nbsubdivision; // nb of subdivision for plot    
    const bool discontinue;

    int const * const DFOnWhat; 
    int const * const DFOfNode; // nu  du df on Node
    int const * const NodeOfDF; // nu du node du df
    int const * const fromFE;   //  the df  come from df of FE
    int const * const fromDF;   //  the df  come from with FE
    int const * const fromASubFE; //  avril 2006 for CL
    int const * const fromASubDF; //  avril 2006 for CL
    int const * const dim_which_sub_fem; // from atomic sub FE for CL
    const  int * ndfOn() const { return & ndfonVertex;}
    
  dataTypeOfFE(const int *nnitemdim,const int dfon[4],int NN,int nbsubdivisionn,int nb_sub_femm=1,bool discon=true);
    // pour evite un template 
    // nitemdim : nbitem : si d==2   3,3,1,0 , si d=3: 4,6,4,1 , si d==1 = 2,1,0,0 
    // dfon : nombre de df par item 
    // NN 
    dataTypeOfFE(const int nitemdim[4],const KN< dataTypeOfFE const *>  & tef);
    
    virtual ~dataTypeOfFE(){ if(dataalloc) delete [] dataalloc;}
 };

template<class RdHat>
class InterpolationMatrix {
public:
  const int N,np,ncoef; 
  bool invariant;
  int k; 
  KN<RdHat> P; 
  KN<R> coef;
  KN<int> comp;
  KN<int> p;
  KN<int> dofe; 

  template<class Mesh>
  InterpolationMatrix(const GFESpace<Mesh> &Vh);

  template<class Mesh>
  InterpolationMatrix(const GTypeOfFE<Mesh> & tef);

  void set(int k);


private: // copie interdit ...
  InterpolationMatrix(const InterpolationMatrix &);
  void operator=(const InterpolationMatrix &);
};

template <class RdHat>
ostream & operator<<(ostream& f,const InterpolationMatrix<RdHat> &M)
{ f<< M.N << " "<< M.k << " "<< M.np << " "<< M.ncoef << endl;
  f<< " = " << M.P ;
  f << "coef=" <<M.coef ;
  f << "comp " << M.comp;
  f << "p = " << M.p;
  f << "dofe = " << M.dofe;
  return f;
}

/*
struct IPJ { 
    int i,p,j; // i is DoF, p is Point, j is componante 
    IPJ(int ii=0,int pp=0,int jj=0):i(ii),p(pp),j(jj) {}
};
inline ostream & operator<<( ostream & f,const IPJ & ipj)
{ f << ipj.i << ' '<< ipj.p << ' '<< ipj.j << " ,  ";
  return f;}
*/
/*
template<class RdHat>
class GFEInterpole 
{
public:
  KN<IPJ> pij_alpha ;
  KN<RdHat> P_Pi_h ;
  double *coef_Pi_h_alpha;
  bool clean;
  
  GFEInterpole(int kPi,int npPi,double * coef_Pi_h_a,bool cc)
    : pij_alpha(kPi),
      P_Pi_h(npPi),
      coef_Pi_h_alpha(coef_Pi_h_a),
      clean(cc)
  {}
  
  GFEInterpole(const KN< GFEInterpole< RdHat> const  *>  & tef, const KN<int> &tefN)
    : pij_alpha(),P_Pi_h(),coef_Pi_h_alpha(0),clean(true)
  {
    
    map<RdHat,int,lessRd> mpt;

    int np=0,nc=0;
    for(int i=0;i<tef.N();++i)
      {
	np += tef[i]->P_Pi_h.N();
	nc += tef[i]->pij_alpha.N();
      }
    
    KN<int>  kpt(np);
    int npp=0,kkk=0;
    for(int i=0,k=0;i<tef.N();++i)
      for(int j=0;j<tef[i]->P_Pi_h.N();++j,++k)
	{	
	  RdHat P(tef[i]->P_Pi_h[j]);
	  if( mpt.find(P) == mpt.end())
	    mpt[P]=npp++; 
	  kpt[kkk++]=mpt[P]; 
	}
    
    pij_alpha.init(nc);
    P_Pi_h.init(npp);
    for(int i=0,k=0;i<tef.N();++i)
      for(int j=0;j<tef[i]->P_Pi_h.N();++j,++k)
	P_Pi_h[kpt[k]]=tef[i]->P_Pi_h[j];
    int op=0,oj=0,oi=0;
    for(int i=0,k=0;i<tef.N();++i)
      {
	for(int j=0;j<tef[i]->pij_alpha.N();++j,++k)
	  {
	    IPJ a =tef[i]->pij_alpha[j];
	    pij_alpha[k]=IPJ(a.i+oi,kpt[a.p+op],a.j+oj);
	  }
	op += tef[i]->P_Pi_h.N();
	oi += tef[i]->pij_alpha.N();
	oj += tefN[i];
      }
    
    cout << "GFEInterpole pij_alpha : n Pt " << npp << " " << pij_alpha.N() <<  endl;
    for(int i=0;i<nc;++i)
      cout << pij_alpha[i] << endl;
    coef_Pi_h_alpha = new double[nc]; 
    assert(op==np);
    assert(oi==nc);
  }
      
  const KN<IPJ > & Ph_ijalpha() const {return pij_alpha;} // la struct de la matrice
  const KN<RdHat> & Pi_h_Rd() const { return P_Pi_h;}   // les points
  
  
  virtual ~GFEInterpole() 
  {
    if(clean) delete[] coef_Pi_h_alpha;
  }
  
  };
  
*/


template<class Mesh>
class GTypeOfFE : public  dataTypeOfFE
{ 
public:
  typedef typename Mesh::Element Element;
  typedef  Element E;
  typedef typename Element::RdHat RdHat;
  typedef typename Element::Rd Rd;
  typedef GFElement<Mesh> FElement;

  // generic data build interpolation
  bool  invariantinterpolationMatrix;
  int NbPtforInterpolation;
  int NbcoefforInterpolation;
  KN<RdHat> PtInterpolation; 
  KN<int> pInterpolation,cInterpolation,dofInterpolation; 
  KN<R>  coefInterpolation;
  
    // soit 
    //  $(U_pj)_{{j=0,N-1}; {p=0,nbpoint_Pi_h-1}}$ les valeurs de U au point Phat[i];
    //   p est le numero du point d'integration
    //   j est la composante 
    //  l'interpole est  defini par
    //  Pi_h u = \sum_k u_k w^k , et u_i = \sum_pj alpha_ipj U_pj
    //  la matrice de alpha_ipj est tres creuse
 


  KN<GTypeOfFE<Mesh> *> Sub_ToFE; 
  KN<int> begin_dfcomp, end_dfcomp;
  KN<int> first_comp,last_comp; // for each SubFE 
    
  virtual void init(InterpolationMatrix<RdHat> & M,FElement * pK,int odf,int ocomp,int *pp) const
  {
    assert(pK==0);
    assert(M.np==NbPtforInterpolation);
    assert(M.ncoef==NbcoefforInterpolation);
    M.P=PtInterpolation;
    M.coef=coefInterpolation;
    M.comp=cInterpolation;
    M.p=pInterpolation;
    M.dofe=dofInterpolation;
  }

  
  // the full constructor ..
  GTypeOfFE(const int dfon[4],const int NN,int  nsub,int kPi,int npPi,bool invar,bool discon) 
    : 
    dataTypeOfFE(Element::nitemdim,dfon, NN, nsub,1,discon),
    invariantinterpolationMatrix(invar),
    NbPtforInterpolation(npPi),
    NbcoefforInterpolation(kPi),
    PtInterpolation(NbPtforInterpolation), 
    pInterpolation(NbcoefforInterpolation),
    cInterpolation(NbcoefforInterpolation),
    dofInterpolation(NbcoefforInterpolation),
    coefInterpolation(NbcoefforInterpolation),

    Sub_ToFE(this->nb_sub_fem),
    begin_dfcomp(N,0),
    end_dfcomp(N,this->NbDoF),
    first_comp(this->nb_sub_fem,0),
    last_comp(this->nb_sub_fem,N)

    { 
      Sub_ToFE=this;
      assert(this->dim_which_sub_fem[this->N-1]>=0 && this->dim_which_sub_fem[this->N-1]< this->nb_sub_fem);
    } 

    //  simple constructeur of lagrange type Finite element   1 point par Node for interpolation
  GTypeOfFE(const int dfon[4],const int NN,int  nsub,bool invar,bool discon)
    : 
    dataTypeOfFE(Element::nitemdim,dfon, NN, nsub,1,discon),
    
    invariantinterpolationMatrix(invar),
    NbPtforInterpolation(this->NbDoF),
    NbcoefforInterpolation(this->NbDoF),
    PtInterpolation(NbPtforInterpolation),
    pInterpolation(NbcoefforInterpolation),
    cInterpolation(NbcoefforInterpolation),
    dofInterpolation(NbcoefforInterpolation),
    coefInterpolation(NbcoefforInterpolation),
    
    Sub_ToFE(this->nb_sub_fem),
    begin_dfcomp(N,0),
    end_dfcomp(N,this->NbDoF),
    first_comp(this->nb_sub_fem,0),
    last_comp(this->nb_sub_fem,N)    
  { 
    Sub_ToFE=this;
    assert(this->dim_which_sub_fem[this->N-1]>=0 && this->dim_which_sub_fem[this->N-1]< this->nb_sub_fem);
  } 

  //  simple constructeur pour sum direct d'espace
  GTypeOfFE(const KN<GTypeOfFE const  *>  & tef)
    : 
    dataTypeOfFE(Element::nitemdim,tef),

    invariantinterpolationMatrix(0),
    NbPtforInterpolation(0),
    NbcoefforInterpolation(ncoeftef(tef)),
    PtInterpolation(NbPtforInterpolation),
    pInterpolation(NbcoefforInterpolation),
    cInterpolation(NbcoefforInterpolation),
    dofInterpolation(NbcoefforInterpolation),
    coefInterpolation(NbcoefforInterpolation),


    Sub_ToFE(this->nb_sub_fem),
    begin_dfcomp(this->N,0),
    end_dfcomp(this->N,this->NbDoF),
    first_comp(this->nb_sub_fem,0),
    last_comp(this->nb_sub_fem,N)        
  { 
    Sub_ToFE=this;
    assert(this->dim_which_sub_fem[this->N-1]>=0 && this->dim_which_sub_fem[this->N-1]< this->nb_sub_fem);
  } 
  
  
  virtual R operator()(const GFElement<Mesh>  & K,const  RdHat & PHat,const KN_<R> & u,int componante,int op) const  ;
  virtual void FB(const What_d whatd,const Mesh & Th,const Element & K,const Rd &P, KNMK_<R> & val) const =0;
  
  static KN<int> tefN(const KN<GTypeOfFE const  *>  & tef) {
    int n=tef.N();
    KN<int> ntef(n);
    for(int i=0;i<n;++i) ntef[i]=tef[i]->N;
    return ntef;}

  static int ncoeftef(const KN<GTypeOfFE const  *>  & tef) {
    int k=0,n=tef.N();
    for(int i=0;i<n;++i)
      k+= tef[i]->NbcoefforInterpolation;
    return k;}


private:
  static int Count(const int *data,int n,int which) 
  {
    int kk=0;
    for (int i=0;i<n;++i)
      if (which == data[i]) ++kk;
    return kk;
  }
  
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

  template<class mesh>
  struct DataFE
  {
    static GTypeOfFE<mesh> & P0; 
    static GTypeOfFE<mesh> & P1; 
    static GTypeOfFE<mesh> & P2; 
  };
  

 
  template<class MMesh>
  class GbaseFElement
  {
  public:
    typedef MMesh  Mesh;
    typedef GFESpace<Mesh>  FESpace;
    typedef typename Mesh::Element Element;
    typedef  Element E;
    typedef typename E::Rd Rd;
    typedef typename E::RdHat RdHat;
    typedef Fem2D::GQuadratureFormular<RdHat>  QF;
    
  const GFESpace<Mesh> &Vh;  
  const Element &T;
  const GTypeOfFE<Mesh> * tfe; 
  const int N;
  const int number;  
  GbaseFElement(const GFESpace<Mesh> &aVh, int k) ;
  GbaseFElement(const   GbaseFElement & K,  const GTypeOfFE<Mesh> & atfe) ;
  R EdgeOrientation(int i) const {return T.EdgeOrientation(i);}
  
};

template<class Mesh>
class GFElement : public GbaseFElement<Mesh> 
{
public:
    
  typedef typename Mesh::Element Element;
  typedef  Element E;
  typedef typename E::Rd Rd;
  typedef typename E::RdHat RdHat;
  typedef Fem2D::GQuadratureFormular<RdHat>  QF;
  
  
  friend class GFESpace<Mesh>;
  const int *p;
  const int nb;
  
  GFElement(const GFESpace<Mesh> * VVh,int k) ;
  
  int NbOfNodes()const {return nb;}
  int  operator[](int i) const ;//{ return  p ? p+i :  ((&T[i])-Vh.Th.vertices);}  N\u00c9 du noeud 
  int  NbDoF(int i) const ;  // number of DF 
  int  operator()(int i,int df) const ;// { N\u00c9 du DoF du noeud i de df local df 
  int  operator()(int df) const { return operator()(NodeOfDF(df),DFOfNode(df));}
  // void Draw(const KN_<R> & U, const  KN_<R> & VIso,int j=0) const ;  
  //void Drawfill(const KN_<R> & U, const  KN_<R> & VIso,int j=0,double rapz=1) const ;  
  //void Draw(const RN_& U,const RN_& V, const  KN_<R> & Viso,R coef,int i0,int i1) const;
  //Rd   MinMax(const RN_& U,const RN_& V,int i0,int i1) const  ;
  //Rd   MinMax(const RN_& U,int i0) const  ;
  void BF(const Rd & P,RNMK_ & val) const;// { tfe->FB(Vh.Th,T,P,val);}
  void BF(const What_d whatd, const Rd & P,RNMK_ & val) const;// { tfe->FB(Vh.Th,T,P,val);}
  
  // add april 08   begin end number for df of the componante ic 
  int dfcbegin(int ic) const { return this->tfe->begin_dfcomp[ic];}
  int dfcend(int ic) const { return this->tfe->end_dfcomp[ic];}
  // the fist and last composant of a sub finite element
  //  int firstcomp(int isfe) const {return this->tfe->first_comp[isfe];}
  // int lastcomp(int isfe)  const {return this->tfe->last_comp[isfe];}
  int subFE(int df)       const {return this->tfe->fromASubFE[df] ;}

  template<class RR>
  KN_<RR> & Pi_h(KNM_<RR> vpt,KN_<RR> & vdf,InterpolationMatrix<RdHat> &M)    const 
  { 
    // compute  the interpolation  
    // in : vpt  value of componant j at point p : vpt(p,j) 
    // out: vdf  value du the degre of freedom
    vdf=RR();
    
    if(M.k != this->number) 
      M.set( this->number); 
    
    for (int k=0;k<M.ncoef;++k)
      vdf[M.dofe[k]] += M.coef[k]*vpt(M.p[k],M.comp[k]);            	
    
    return  vdf;     
  }
  
  
  int NbDoF() const { return this->tfe->NbDoF;}
  int DFOnWhat(int i) const { return this->tfe->DFOnWhat[i];}
  int FromDF(int i) const { return this->tfe->fromDF[i];}
  int FromFE(int i) const { return this->tfe->fromFE[i];}
  
  // df is the df in element 
  int NodeOfDF(int df) const { return this->tfe->NodeOfDF[df];} // a node
  int FromASubFE(int i) const { return this->tfe->fromASubFE[i];}
  int FromASubDF(int i) const { return this->tfe->fromASubDF[i];}
  int DFOfNode(int df) const { return this->tfe->DFOfNode[df];} // the df number on the node 
  
  R operator()(const Rd & PHat,const KN_<R> & u,int i,int op)  const ;
  complex<R> operator()(const RdHat & PHat,const KN_<complex<R> > & u,int i,int op)  const ;
  
  // GFElementGlobalToLocal operator()(const KN_<R> & u ) const { return GFElementGlobalToLocal(*this,u);}
private:
  int nbsubdivision() const { return this->tfe->nbsubdivision;} // for draw 
};


template<class MMesh> 
    class BuildTFE { protected:
	GTypeOfFE<MMesh> const * const tfe;
    };

template<class MMesh>
struct GFESpacePtrTFE {
    GTypeOfFE<MMesh> const * const ptrTFE;
    GFESpacePtrTFE(GTypeOfFE<MMesh> const * const pptrTFE=0) : ptrTFE(pptrTFE) {}
    virtual  ~GFESpacePtrTFE() { if(ptrTFE) delete ptrTFE;}
};
    
template<class MMesh> 
class GFESpace :  public RefCounter,protected GFESpacePtrTFE<MMesh>,  public DataFENodeDF, public UniqueffId  {
public:
  typedef MMesh Mesh;
  typedef GFElement<Mesh> FElement;
  typedef typename Mesh::Element  Element;
  typedef typename Mesh::BorderElement  BorderElement;
  typedef typename Mesh::Rd  Rd;
  typedef GTypeOfFE<Mesh> TypeOfFE;
  typedef GQuadratureFormular<typename Element::RdHat> QFElement;
  typedef GQuadratureFormular<typename BorderElement::RdHat>  QFBorderElement;

  const Mesh &Th;

  KN<const GTypeOfFE<Mesh> *>  TFE; 
private:
  
public:
  CountPointer<const Mesh> cmesh;
  const int N; // dim espace d'arrive
  const int Nproduit; // dim de l'espace produit generalement 1
  const int nb_sub_fem; // nb de sous elements finis tensorise (independe au niveau des composantes)
  int const* const dim_which_sub_fem;// donne les dependant des composantes liee a un meme sous element fini
  const int   maxNbPtforInterpolation;  
  const int   maxNbcoefforInterpolation;  

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
  
  //  par defaut P1                
  
    GFESpace(const Mesh & TTh,const GTypeOfFE<Mesh> & tfe=DataFE<Mesh>::P1)
    :
    GFESpacePtrTFE<MMesh>(0),
    DataFENodeDF(TTh.BuildDFNumbering(tfe.ndfonVertex,tfe.ndfonEdge,tfe.ndfonFace,tfe.ndfonVolume)),
    Th(TTh),
    TFE(1,0,&tfe), 
    cmesh(TTh),
    N(TFE[0]->N),
    Nproduit(TFE[0]->N),
    nb_sub_fem(TFE[0]->nb_sub_fem),
    dim_which_sub_fem(TFE[0]->dim_which_sub_fem),
    maxNbPtforInterpolation(TFE[0]->NbPtforInterpolation),
    maxNbcoefforInterpolation(TFE[0]->NbcoefforInterpolation)
    
    {
	if(verbosity) cout << "  -- FESpace: Nb of Nodes " << NbOfNodes 
	    << " Nb of DoF " << NbOfDF << endl;
    }
    
  GFESpace(const GFESpace & Vh,int kk);
  GFESpace(const GFESpace ** Vh,int kk);
    
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

  
  ~GFESpace()
  {
  }
  // GFESpace(Mesh & TTh,int NbDfOnSommet,int NbDfOnEdge,int NbDfOnElement,int NN=1); 
  int  renum();
  
  FElement operator[](int k) const { return FElement(this,k);} 
  FElement operator[](const Element & K) const { return FElement(this,Th(K));} 
  int  operator()(int k)const {return NbOfNodesInElement(k);}
  int  operator()(int k,int i) const { //  the node i of element k
      return NodesOfElement ?  *(PtrFirstNodeOfElement(k) + i)  : Th(k,i)  ;}
  
    /*
      void Draw(const KN_<R>& U,const KN_<R>& Viso,int j=0,float *colors=0,int nbcolors=0,bool hsv=true,bool drawborder=true) const ; // Draw iso line
    void Drawfill(const KN_<R>& U,const KN_<R>& Viso,int j=0,double rapz=1,float *colors=0,int nbcolors=0,bool hsv=true,bool drawborder=true) const ; // Draw iso line
    
    void Draw(const KN_<R>& U,const KN_<R> & Viso, R coef,int j0=0,int j1=1,float *colors=0,int nbcolors=0,bool hsv=true,bool drawborder=true) const  ; // Arrow
    void Draw(const KN_<R>& U,const KN_<R>& V,const KN_<R> & Viso, R coef,int iu=0,int iv=0,float *colors=0,int nbcolors=0,bool hsv=true,bool drawborder=true) const  ; // Arrow
    Rd MinMax(const KN_<R>& U,const KN_<R>& V,int j0,int j1,bool bb=true) const ;
    Rd MinMax(const KN_<R>& U,int j0, bool bb=true) const ;
    // void destroy() {RefCounter::destroy();}
    */
    bool isFEMesh() const { return ! NodesOfElement  && ( N==1) ;} // to make optim
  KN<R>  newSaveDraw(const KN_<R> & U,int composante,int & lg,KN<Rd> &Psub,KN<int> &Ksub,int op_U=0) const  ; 


private: // for gibbs  
  int gibbsv (long* ptvoi,long* vois,long* lvois,long* w,long* v);
};

template<class Mesh>
inline GbaseFElement<Mesh>::GbaseFElement(  const GFESpace<Mesh> &aVh, int k) 
  : Vh(aVh),T(Vh.Th[k]),tfe(aVh.TFE[k]),N(aVh.N),number(k){}

template<class Mesh>    
inline GbaseFElement<Mesh>::GbaseFElement(const   GbaseFElement & K,  const GTypeOfFE<Mesh> & atfe) 
  : Vh(K.Vh),T(K.T),tfe(&atfe),N(Vh.N),number(K.number){}

template<class Mesh>
GFElement<Mesh>::GFElement(const GFESpace<Mesh> * VVh,int k) 
  : GbaseFElement<Mesh>(*VVh,k) ,
    p(this->Vh.PtrFirstNodeOfElement(k)),
    nb(this->Vh.NbOfNodesInElement(k))    
{} 

template<class Mesh>
inline   int  GFElement<Mesh>::operator[](int i) const {
  return  p ? p[i] :  ((&this->T[i])-this->Vh.Th.vertices);}  

template<class Mesh>
inline   int  GFElement<Mesh>::operator()(int i,int df) const {
  return  this->Vh.FirstDFOfNode(p ? p[i] :  ((&this->T[i])-this->Vh.Th.vertices)) + df;}  

template<class Mesh>
inline   int  GFElement<Mesh>::NbDoF(int i) const {
  int node =p ? p[i] :  ((&this->T[i])-this->Vh.Th.vertices);
  return  this->Vh.LastDFOfNode(node)-this->Vh.FirstDFOfNode(node);}  

template<class Mesh>
inline void GFElement<Mesh>::BF(const Rd & P,RNMK_ & val) const {
  this->tfe->FB(Fop_D0|Fop_D1,this->Vh.Th,this->T,P,val);}

template<class Mesh>
inline void GFElement<Mesh>::BF(const What_d whatd,const Rd & P,RNMK_ & val) const { this->tfe->FB(whatd,this->Vh.Th,this->T,P,val);}

template<class Mesh>
inline   R GFElement<Mesh>::operator()(const Rd & PHat,
                                const KN_<R> & u,int i,int op)  const
{
    return (*this->tfe)(*this,PHat,u,i,op);
}

template<class Mesh>
inline  complex<R> GFElement<Mesh>::operator()(const RdHat & PHat,const KN_<complex<R> > & u,int i,int op)  const 
{
    complex<double> * pu=u; // pointeur du tableau
    double *pr = static_cast<double*>(static_cast<void*>(pu));
    
    const KN_<R>  ur(pr,u.n,u.step*2);
    const KN_<R>  ui(pr+1,u.n,u.step*2);
    
    return complex<R>((*this->tfe)(*this,PHat,ur,i,op),(*this->tfe)(*this,PHat,ui,i,op));
}



template<class Mesh>
R GTypeOfFE<Mesh>::operator()(const GFElement<Mesh> & K,const  RdHat & PHat,const KN_<R> & u,int componante,int op) const 
{
    R v[5000],vf[100];
    assert(N*last_operatortype*NbDoF<=5000 && NbDoF <100 );
    KNMK_<R> fb(v,NbDoF,N,last_operatortype); //  the value for basic fonction
    KN_<R> fk(vf,NbDoF);
    for (int i=0;i<NbDoF;i++) // get the local value
	fk[i] = u[K(i)];
    //  get value of basic function
    FB( 1 << op ,K.Vh.Th,K.T,PHat,fb);  
    R r = (fb('.',componante,op),fk);  
    return r;
}

int nbdf_d(const int ndfitem[4],const  int nd[4]);
int nbnode_d(const int ndfitem[4],const  int nd[4]);
int *builddata_d(const int ndfitem[4],const int nd[4]);


    
    
template<class Mesh>
class GTypeOfFESum:  public  GTypeOfFE<Mesh> {
    
public:
  typedef GFElement<Mesh> FElement;
  typedef typename Mesh::Element  Element;
  typedef typename Element::RdHat  RdHat;
  typedef typename Mesh::Rd  Rd;
  const int k; 
  KN<const  GTypeOfFE<Mesh> *> teb;
  KN<int> NN,DF,comp,numPtInterpolation;
 
  GTypeOfFESum(const KN< GTypeOfFE<Mesh> const *> & t);
  GTypeOfFESum(const GFESpace<Mesh> **,int kk);
  GTypeOfFESum(const GFESpace<Mesh> &,int kk);
    
  void Build();  // the true constructor 
    
  void init(InterpolationMatrix<RdHat> & M,FElement * pK=0,int odf=0,int ocomp=0,int *pp=0) const;
   
  void FB(const What_d whatd,const Mesh & Th,const Element & K,const Rd &P, KNMK_<R> & val) const ;
  ~GTypeOfFESum(){}
} ;


template<class Mesh>
void GTypeOfFESum<Mesh>::FB(const What_d whatd,const Mesh & Th,const Element & K,const Rd &P, KNMK_<R> & val) const  
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
	  assert(ni<nii && di < dii);
	  RNMK_ v(val(SubArray(dii-di,di),SubArray(nii-ni,ni),t));     
	  if (j<=i)
	      teb[i]->FB(whatd,Th,K,P,v);       
	  else
	      v=val(SubArray(DF[j+1]-DF[j],DF[j]),SubArray(NN[j+1]-NN[j],NN[j]),t);     
      }
} 

/* il faut reflechir   mais a faire
 // une class qui mettre une chapeau sur les type d'element finis. 
template<class Mesh>
struct TypeOfFE {
    
    GTypeOfFE<Mesh> *tef; 
    bool clean;
    TypeOfFE(GTypeOfFE<Mesh> & t):tef(&t),clean(false){}
    TypeOfFE(GTypeOfFE<Mesh> * t):tef(t),clean(false){}
    TypeOfFE(GTypeOfFE<Mesh> & t):tef(t),clean(false){}
    
    ~TypeOfFE() {if(clean) delete tef;}
}    
*/

template<class RdHat>
template<class Mesh>
InterpolationMatrix<RdHat>::InterpolationMatrix(const GFESpace<Mesh> &Vh)
  :
  N(Vh.N),np(Vh.maxNbPtforInterpolation),ncoef(Vh.maxNbcoefforInterpolation),
  invariant(Vh.TFE.constant() ? Vh.TFE[0]->invariantinterpolationMatrix: false),
  k(-1),
  P(np),
  comp(ncoef),
  p(ncoef),
  dofe(ncoef)
{ 
  Vh.TFE[0]->GTypeOfFE<Mesh>::init(*this,0,0,0,0);    
}

template<class RdHat>
template<class Mesh>
InterpolationMatrix<RdHat>::InterpolationMatrix(const GTypeOfFE<Mesh> & tef)
  :
  N(tef.N),np(tef.NbPtforInterpolation),ncoef(tef.NbcoefforInterpolation),
  invariant(tef.invariantinterpolationMatrix),
  k(-1),
  P(np),
  comp(ncoef),
  p(ncoef),
  dofe(ncoef)  
{
  //  virtual void init(InterpolationMatrix<RdHat> & M,FElement * pK=0,int odf=0,int ocomp=0,int *pp=0) const
  tef.GTypeOfFE<Mesh>::init(*this,0,0,0,0);
}

template<class RdHat>
void InterpolationMatrix<RdHat>::set(int kk)
{
  if(k==kk) return;
  k=kk;
  if(invariant) return;
  assert(invariant);
  // a faire ...
}

typedef  GTypeOfFE<Mesh3> TypeOfFE3;
typedef  GTypeOfFE<Mesh3> TypeOfFE3;
typedef  GFESpace<Mesh3> FESpace3;
typedef  GFESpace<Mesh2> FESpace2;
typedef  GFElement<Mesh3> FElement3;
typedef  GFElement<Mesh2> FElement2;
typedef  GFElement<Mesh3> FElement3;
typedef  GbaseFElement<Mesh2> baseFElement2;
typedef  GbaseFElement<Mesh3> baseFElement3;

}


#endif
