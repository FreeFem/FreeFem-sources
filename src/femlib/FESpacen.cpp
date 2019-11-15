// ********** DO NOT REMOVE THIS BANNER **********
// ORIG-DATE:     Jan 2008
// -*- Mode : c++ -*-
//
// SUMMARY  : Generic Fiinite Element   1d, 2d, 3d  
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


#include <map>

#include "ufunction.hpp"

#include "error.hpp"
#include "RNM.hpp"

#include "Mesh3dn.hpp"
#include "Mesh2dn.hpp"

#include "FESpacen.hpp"

#include "splitsimplex.hpp"
 int UniqueffId::count=0;
 namespace Fem2D {

//template<class Element>
int nbdf_d(const int ndfitem[4],const  int nd[4])
{
  const int ndf = ndfitem[0]*nd[0] + ndfitem[1]*nd[1]+  ndfitem[2]*nd[2]  + ndfitem[3]*nd[3];	
  return ndf;
}

//template<class Element>
int nbnode_d(const int ndfitem[4],const  int nd[4])
{
  //  const  int nd[]= {Element::nv, Element::ne,Element::nf,Element::nt};
    const int ndf = nd[0]*(ndfitem[0]!=0) + nd[1]*(ndfitem[1]!=0)+  nd[2]*(ndfitem[2]!=0)  + nd[3]*(ndfitem[3]!=0);	
    return ndf;
}

//template<class Element>
int *builddata_d(const int ndfitem[4],const int nd[4],int N)
{
  //    const int d=Element::Rd::d;
  //    const int nwhat=Element::nitem;
  //    const  int nd[]= {Element::nv, Element::ne,Element::nf,Element::nt};
  //    const int nitem=nd[0]+nd[1]+nd[2]+nd[3];
  //    cout << " nitem="<< nitem<< endl;
  const int ndf = nbdf_d(ndfitem,nd);
  const int nnode=nbnode_d(ndfitem,nd);
  int lgdata= ndf*5+N;
  int * data = new int[lgdata];
  int p=0;
  for(int i=0,nw=0;i<=3;++i)
    for(int j=0;j<nd[i];++j,nw++)// pour les what (support de node)
      for(int k=0;k<ndfitem[i];++k)// 	    
	data[p++] = nw;
 // cout << p << " " <<ndfitem[0]<< ndfitem[1]<<ndfitem[2]<<ndfitem[3]<< " " << ndf  << endl;
  assert(p==ndf);
  for(int i=0,nw=0;i<=3;++i)
    for(int j=0;j<nd[i];++j,nw++)// pour les what (support de node)
      for(int k=0;k<ndfitem[i];++k)// 	    
	data[p++] = k;
  // cout << p << " " << 2*ndf << " " << nitem << endl;
  int nn=0;
  for(int i=0;i<=3;++i)
    {
      int in=ndfitem[i]?1:0;
      for(int j=0;j<nd[i];++j,nn+=in )// pour les what (support de node)
	for(int k=0;k<ndfitem[i];++k)// 	    
	  data[p++] = nn;
    }
  // cout << p << " " << 3*ndf << " " << nitem << endl;
  for(int i=0;i<ndf*2+N;++i)	    
    data[p++] = 0;
  
 // data[p++] = 0;
 // data[p++] = 0;
  //cout << p << " == " << lgdata << endl;
  assert(p== lgdata);
  //cout << nn << " " << nnode << endl;
  p =0;
  
  /*  
      for(int j=0;j<4;j++)
      {
      for (int i=0;i<ndf;i++)
      cout << data[p++] << " ";
      cout << endl;
      } 
  */
  assert(nn==nnode); 
  return data;  
}

   dataTypeOfFE::dataTypeOfFE(const int nitemdim[4],const int dfon[4],int NN,int nbsubdivisionn,int nb_sub_femm,bool discon)
     :
     data(builddata_d(dfon,nitemdim,NN)),  
     dataalloc(data),
     ndfonVertex(dfon[0]),
     ndfonEdge(dfon[1]),
     ndfonFace(dfon[2]),
     ndfonVolume(dfon[3]),
     NbDoF(nbdf_d(dfon,nitemdim)),
     NbNode(nbnode_d(dfon,nitemdim)),
     N(NN),
     nb_sub_fem(nb_sub_femm),
     nbsubdivision(nbsubdivisionn),
     discontinue(discon),
     DFOnWhat(data+0*NbDoF),
     DFOfNode(data+1*NbDoF),
     NodeOfDF(data+2*NbDoF),
     fromFE(data+3*NbDoF),
     fromDF(data+4*NbDoF),
     fromASubFE(data+3*NbDoF),
     fromASubDF(data+4*NbDoF) ,
     dim_which_sub_fem(data+5*NbDoF)
   {}


int *builddata_d(const int nitemdim[4],const KN< dataTypeOfFE const  *> &teb)
{
    const int k = teb.N(); 
    KN<int> NN(k+1), DF(k+1) , comp(k+1);
    map< dataTypeOfFE const  *,int> m;
    int i=k,j;    
    while(i--) // on va a l'envert pour avoir comp[i] <=i 
	m[teb[i]]=i;
    // l'ordre comp est important comp est croissant  mais pas de pb. 
    i=k;    
    while(i--) 
	comp[i]=m[teb[i]]; //  comp[i] <=i
    int n=0,N=0;
    for ( j=0;j<k;j++)
      {NN[j]=N;N+=teb[j]->N;}
    NN[k] = N;
    //  reservation des interval en df   
    n=0;
    for ( j=0;j<k;j++)
      { DF[j]=n;n+=teb[j]->NbDoF;}
    DF[k] = n;
    
    int NbDoF=0;
    int dfon[4]={0,0,0,0};
    int nbsubdivision=0;
    int discon=0; 
    for (int i=0;i<k;++i)
      {
	NbDoF += teb[i]->NbDoF;
	  dfon[0] += teb[i]->ndfonVertex;
	dfon[1] += teb[i]->ndfonEdge;
	dfon[2] += teb[i]->ndfonFace;
	dfon[3] += teb[i]->ndfonVolume;
	nbsubdivision = max(nbsubdivision,teb[i]->nbsubdivision);
	discon = discon || teb[i]->discontinue; // bof bof 1 FE discontinue => discontinue
      }
    int nwhat=15; // 15 = 4+6+1+1 (nb of  support item  (what) : vertex, edges, fqces, tet)
    
    int ostart=nwhat;
    int * data0=new int[ostart+7*NbDoF+N];
    int * data=data0+ostart;
    int * data1=data+5*NbDoF;
   
      int c=0;
    KN<int> w(nwhat),nn(nwhat); 
    
    w=0;
    nn=0; 

    
    for ( j=0;j<k;j++)
	for ( i=0;i<teb[j]->NbDoF;i++)
	    nn[teb[j]->DFOnWhat[i]]++;
    int nbn=0;      
    for( j=0;j<nwhat;j++)
	if (nn[j]) nn[j]=nbn++;
	else nn[j]=-1;
    KN<int> dln(nwhat);
    dln=0;
    // nn donne numero de noeud sur what            
    for ( j=0;j<k;j++)
	for ( i=0;i<teb[j]->NbDoF;i++)
	    data[c++] = teb[j]->DFOnWhat[i];
    
    for ( j=0;j<k;j++)
      {
	  int  cc=c;
	  for ( i=0;i<teb[j]->NbDoF;i++)
	      data[c++] = teb[j]->DFOfNode[i]+dln[teb[j]->DFOnWhat[i]];
	  for ( i=0;i<teb[j]->NbDoF;i++)
	      dln[teb[j]->DFOnWhat[i]]=Max(dln[teb[j]->DFOnWhat[i]],data[cc++]+1);      
      }
    
    
    for ( j=0;j<k;j++)
      { 
	  //  w renumerotation des noeuds 
	  //  Ok si un noeud par what 
	  for ( i=0;i<teb[j]->NbDoF;i++)
	      data[c++] = nn[teb[j]->DFOnWhat[i]];
      }
    
    for ( j=0;j<k;j++)
	for ( i=0;i<teb[j]->NbDoF;i++)
	    data[c++] = j; //  node from of FE
    
    
    for ( j=0;j<k;j++)
	for ( i=0;i<teb[j]->NbDoF;i++)
	    data[c++] = i; //  node from of df in FE
    // error -- here 
    //in case of [P2,P2],P1  
    // we expect 0,0,1   and we get 0 1 2 
    // => wrong BC ???? 
    c+=2*n; // on saute le deux tableau en plus (cf data1.)
    
    
    int xx=0;
    for (j=0;j<k;j++)
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
    int cj=0;
    int ccc=0;
    for ( j=0;j<k;ccc+=teb[j++]->nb_sub_fem)
	
	for ( i=0;i<teb[j]->NbDoF;i++)
	  {
	      int il= teb[j]->fromASubDF[i];
	      int jl= teb[j]->fromASubFE[i];
	      data1[ci++]=il;
	      data1[cj++]=ccc+jl;      
	  }
    
    int nb_sub_fem=ccc;
    
    ffassert(c== 7*n+N);      
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
    
    
    for(int i=0;i<4;++i)
	data0[i]=dfon[i];
    
    data0[4]=NbDoF;
    data0[5]=nbn;// NbNode
    data0[6]=N;
    data0[7]=nb_sub_fem;
    data0[8]=nbsubdivision;
    data0[9]=discon;
    
    return data0;
}

  
dataTypeOfFE::dataTypeOfFE(const int nitemdim[4],const KN< dataTypeOfFE const *>  &  tef)
: 
data(builddata_d(nitemdim,tef)),  
dataalloc(data),
ndfonVertex(data[0]),
ndfonEdge(data[1]),
ndfonFace(data[2]),
ndfonVolume(data[3]),
NbDoF(data[4]),
NbNode(data[5]),
N(data[6]),
nb_sub_fem(data[7]),
nbsubdivision(data[8]),
discontinue(data[9]),
DFOnWhat(data+15+0*NbDoF),
DFOfNode(data+15+1*NbDoF),
NodeOfDF(data+15+2*NbDoF),
fromFE(data+15+3*NbDoF),
fromDF(data+15+4*NbDoF),
fromASubFE(data+15+5*NbDoF),
fromASubDF(data+15+6*NbDoF) ,
dim_which_sub_fem(data+15+7*NbDoF)
{}

template<class Mesh>
void GTypeOfFESum<Mesh>::init(InterpolationMatrix<RdHat> & M,FElement * pK,int odf,int ocomp,int *pp) const
{
  // a faire ..... cas matrix invariante 
  assert(0);
}

template<class Mesh>
     GTypeOfFESum<Mesh>::GTypeOfFESum(const KN< GTypeOfFE<Mesh> const *> & t)
     : 
     GTypeOfFE<Mesh>(t),
     k(t.N()),
     teb(t),
     NN(k+1),
     DF(k+1) ,
     comp(k+1) {Build();}
     
template<class Mesh> 
static  KN< GTypeOfFE<Mesh> const *> kn(const GFESpace<Mesh> ** tt,int kk)
     {
	 KN< GTypeOfFE<Mesh> const *> r(kk);
	 for(int i=0;i<kk;++i)
	   { r[i]=tt[i]->TFE[0];ffassert(tt[i]->TFE.constant());}
	 return r;
     }
template<class Mesh> 
static     KN< GTypeOfFE<Mesh> const *> kn(const GFESpace<Mesh> & tt,int kk)
     {
	 return  KN< GTypeOfFE<Mesh> const *> (kk,tt.TFE[0]);
     }
     
template<class Mesh>
     GTypeOfFESum<Mesh>::GTypeOfFESum(const GFESpace<Mesh> ** tt,int kk)
     :	
     GTypeOfFE<Mesh>(kn(tt,kk)),
     k(kk),
     teb(kn(tt,kk)),
     NN(k+1),
     DF(k+1) ,
     comp(k+1) {Build();}
     
template<class Mesh>
     GTypeOfFESum<Mesh>::GTypeOfFESum(const GFESpace<Mesh> & tt,int kk)
     :	
     GTypeOfFE<Mesh>(kn(tt,kk)),
     k(kk),
     teb(kn(tt,kk)),
     NN(k+1),
     DF(k+1) ,
     comp(k+1) {Build();}
     
template<class Mesh>
void GTypeOfFESum<Mesh>::Build()
{
    bool debug=verbosity>5;;
  {
    const KN< GTypeOfFE<Mesh> const *> & t=teb;
    map<const GTypeOfFE<Mesh> *,int> m;
    int i=k,j;    
    while(i--) // on va a l'envert pour avoir comp[i] <=i 
      m[teb[i]]=i;
    // l'ordre comp est important comp est croissant  mais pas de pb. 
    i=k;    
    while(i--) 
      comp[i]=m[teb[i]]; //  comp[i] <=i
    
    // reservatition des intervalles en espaces
    int n=0,N=0;
    for ( j=0;j<k;j++)
      {NN[j]=N;N+=teb[j]->N;}
    NN[k] = N;
    //  reservation des interval en df   
    n=0;
    for ( j=0;j<k;j++)
      { DF[j]=n;n+=teb[j]->NbDoF;}
    DF[k] = n;
  }
  int ii=0;
  for (int i=0;i<k;++i)
    {
      for (int j=0;j<teb[i]->nb_sub_fem;++j)
	this->Sub_ToFE[ii++]=teb[i]->Sub_ToFE[j];
    }
  assert(ii==this->nb_sub_fem );
  
  int c=0,c0=0, fcom=0;
  for (int i=0;i<this->nb_sub_fem;i++) 
    { 
      int N=this->Sub_ToFE[i]->N;
      int ndofi=this->Sub_ToFE[i]->NbDoF;
      this->first_comp[i]= fcom;
      this->last_comp[i]= fcom+N;
      fcom += N;

      for(int j=0;j<N;++j)
	{
	  this->begin_dfcomp[c] = c0 + this->Sub_ToFE[i]->begin_dfcomp[j] ; 
	  this->end_dfcomp[c]   = c0 + this->Sub_ToFE[i]->end_dfcomp[j] ;
	  c++;
	}
      c0+=ndofi;
      
    }
  if(debug)
    {
      cout <<" NbDoF : " << this->NbDoF <<endl;
      for(int i=0;i<this->N;++i)
	cout << "      comp " << i << " ["<<this->begin_dfcomp[i]<<", "<< this->end_dfcomp[i]<< "[\n";
    }
  
  // construction de l'interpolation .
  
  int npi=0;
  int nci=0;
  bool var=true;
  for (int i=0;i<this->nb_sub_fem;i++)
    {
      npi +=this->Sub_ToFE[i]->NbPtforInterpolation;
      nci +=this->Sub_ToFE[i]->NbcoefforInterpolation;
      var = var && this->Sub_ToFE[i]->invariantinterpolationMatrix;
    }
  assert(this->NbcoefforInterpolation== nci);
  this->invariantinterpolationMatrix=var;
  // this->pInterpolation.init(nci);
  // this->cInterpolation.init(nci);
  // this->dofInterpolation.iniy(nci);
  KN<int> opi(this->nb_sub_fem);// offset numumber point intgartion
  {
    map<RdHat,int,lessRd> mpt;
    numPtInterpolation.init(npi);
    int npp=0,kkk=0;
    KN<RdHat> Ptt(npi);
    for (int i=0;i<this->nb_sub_fem;i++)
      {
        opi[i]=kkk;
	const GTypeOfFE<Mesh> &ti=*this->Sub_ToFE[i];
	
	for(int p=0;p<ti.NbPtforInterpolation;++p,++kkk)
	  {
	    Ptt[kkk]=ti.PtInterpolation[p];
	    if( mpt.find(Ptt[kkk]) == mpt.end())
	      mpt[Ptt[kkk]]=npp++;
	    numPtInterpolation[kkk]=mpt[Ptt[kkk]];
              if(verbosity>100)
                  cout << "    p= "<< p << " [ " << Ptt[kkk]<< "] ,  "<< kkk<< " "<< npp<< " " << numPtInterpolation[kkk]<< endl;;

	  }
      }
    assert(this->NbPtforInterpolation==0);
    if(verbosity>5)
    cout << npp;
    this->NbPtforInterpolation=npp;
    this->PtInterpolation.init(npp);
    for(int i=0;i<npp;++i)
      this->PtInterpolation[numPtInterpolation[i]]=Ptt[i];
  }
  
  int oc=0,odof=0;
  for (int i=0,k=0;i<this->nb_sub_fem;i++)
    {
      const GTypeOfFE<Mesh> &ti=*this->Sub_ToFE[i];
      for(int j=0;j<ti.NbcoefforInterpolation; ++j,++k)
	{
	  this->pInterpolation[k]   = numPtInterpolation[opi[i]+ti.pInterpolation[j]];
	  this->cInterpolation[k]   = ti.cInterpolation[j]+oc;
	  this->dofInterpolation[k] = ti.dofInterpolation[j]+odof;
	  this->coefInterpolation[k]=ti.coefInterpolation[j];
	}
      oc += ti.N;
      odof += ti.NbDoF; 
    }
    if(verbosity>100)
    cout << " **GTypeOfFESum<Mesh>::Build() " <<this->pInterpolation <<endl;
  assert(c==this->N);
}

     
template<class Mesh> void GTypeOfFESum<Mesh>::set(const Mesh & Th,const Element & K,InterpolationMatrix<RdHat> & M,int oocoef,int oodf,int *nnump ) const
     {
	 int op=0,oc=0,odof=oodf,ocoef=oocoef;
	 assert(nnump==0);
	 for (int i=0,k=0;i<this->nb_sub_fem;i++)
	   {
	       const GTypeOfFE<Mesh> &ti=*this->Sub_ToFE[i];
	       if(!ti.invariantinterpolationMatrix)
		   ti.set(Th,K,M,ocoef,odof,&numPtInterpolation[op]);
	       oc += ti.N;
	       odof += ti.NbDoF; 
	       ocoef += ti.NbcoefforInterpolation;
	       op += ti.NbPtforInterpolation;
	       
	   }
         if( verbosity > 100) cout << " GTypeOfFESum set "<< this->coefInterpolation << endl;
     }
     
template<class MMesh> 
     GFESpace<MMesh>::GFESpace(const GFESpace & Vh,int kk,int nbequibe,int *equibe)
     :
     GFESpacePtrTFE<MMesh>(new GTypeOfFESum<MMesh>(Vh,kk)),
     DataFENodeDF(Vh.Th.BuildDFNumbering(this->ptrTFE->ndfonVertex,this->ptrTFE->ndfonEdge,this->ptrTFE->ndfonFace,this->ptrTFE->ndfonVolume,nbequibe,equibe)),
     Th(Vh.Th),
     TFE(1,0,this->ptrTFE), 
     cmesh(Th),
     N(TFE[0]->N),
     Nproduit(kk),
     nb_sub_fem(TFE[0]->nb_sub_fem),
     dim_which_sub_fem(TFE[0]->dim_which_sub_fem),
     maxNbPtforInterpolation(TFE[0]->NbPtforInterpolation),
     maxNbcoefforInterpolation(TFE[0]->NbcoefforInterpolation)
     
     {
     }
    
template<class MMesh> 
     GFESpace<MMesh>::GFESpace(const GFESpace ** pVh,int kk,int nbequibe,int *equibe)
     :
     GFESpacePtrTFE<MMesh>(new GTypeOfFESum<MMesh>(pVh,kk)),
     DataFENodeDF((**pVh).Th.BuildDFNumbering(this->ptrTFE->ndfonVertex,this->ptrTFE->ndfonEdge,this->ptrTFE->ndfonFace,this->ptrTFE->ndfonVolume,nbequibe,equibe)),
     Th((**pVh).Th),
     TFE(1,0,this->ptrTFE), 
     cmesh(Th),
     N(TFE[0]->N),
     Nproduit(FirstDfOfNodeData ? 1 :MaxNbDFPerNode),
     nb_sub_fem(TFE[0]->nb_sub_fem),
     dim_which_sub_fem(TFE[0]->dim_which_sub_fem),
     maxNbPtforInterpolation(TFE[0]->NbPtforInterpolation),
     maxNbcoefforInterpolation(TFE[0]->NbcoefforInterpolation)
     
     {
         long snbdf=0;
         for(int i=0;i<kk;++i)
             snbdf += pVh[i]->NbOfDF;
         if( snbdf !=NbOfDF)
             cerr << " Problem build of GFESpace (3d) (may be : due to periodic Boundary condition missing ) FH " << endl
             << " The number of DF must be " << snbdf << "  and it is " << NbOfDF <<endl; 
         ffassert(snbdf == NbOfDF );

         
	     for(int i=0;i<kk;++i)
	     ffassert(&Th==&pVh[i]->Th);
     }
   
     template<class MMesh>
     template<class R>
     KN<R>  GFESpace<MMesh>::newSaveDraw(const KN_<R> & U,int componante,int & lg,KN<typename MMesh::RdHat> &Psub,KN<int> &Ksub,int op_U) const
     {
         const int d =  Rd::d;
         Rd *Ps=0;
         int *Ks=0;
         int nsb = TFE[0]->nbsubdivision;
         int nvsub,nksub;
         SplitSimplex<Rd>(nsb, nvsub,  Ps,  nksub ,  Ks);
         ffassert( Psub.unset());
         ffassert( Ksub.unset());
         Psub.set(Ps,nvsub);
         Ksub.set(Ks,nksub*(d+1));
         lg= nvsub*Th.nt;
         KN<R> v(lg);
         for (int k=0,i=0;k<Th.nt;k++)
         {
             FElement K=(*this)[k];
             for(int l=0;l<nvsub;l++)
                 v[i++] =   K(Psub[l], U, componante, op_U)  ;
             
         }
         return KN<R>(true,v);// to remove the copy.
     }
    
     
     template< >
     template<class R>
     KN<R>  GFESpace<MeshS>::newSaveDraw(const KN_<R> & U,int componante,int & lg,KN<typename MeshS::RdHat> &Psub,KN<int> &Ksub,int op_U) const
     {
         typedef typename MeshS::RdHat RdHat;
         const int dHat =  RdHat::d;
         MeshS::RdHat  *Ps=0;
         int *Ks=0;
         int nsb = TFE[0]->nbsubdivision;
         int nvsub,nksub;
         SplitSimplex<RdHat>(nsb, nvsub,  Ps,  nksub ,  Ks); 
         ffassert( Psub.unset());
         ffassert( Ksub.unset());
         Psub.set(Ps,nvsub);
         Ksub.set(Ks,nksub*(dHat+1));
         lg= nvsub*Th.nt;
         KN<R> v(lg);
         for (int k=0,i=0;k<Th.nt;k++)
         {
             FElementS K=(*this)[k];
            for(int l=0;l<nvsub;l++)
                 v[i++] =   K(Psub[l], U, componante, op_U)  ;
             
         }
         return KN<R>(true,v);// to remove the copy.
     }
   
     
     template< >
     template<class R>
     KN<R>  GFESpace<MeshL>::newSaveDraw(const KN_<R> & U,int componante,int & lg,KN<typename MeshL::RdHat> &Psub,KN<int> &Ksub,int op_U) const
     {
         typedef typename MeshL::RdHat RdHat;
         const int dHat =  RdHat::d;
         MeshL::RdHat  *Ps=0;
         int *Ks=0;
         int nsb = TFE[0]->nbsubdivision;
         int nvsub,nksub;
         SplitSimplex<RdHat>(nsb, nvsub,  Ps,  nksub ,  Ks);
         ffassert( Psub.unset());
         ffassert( Ksub.unset());
         Psub.set(Ps,nvsub);
         Ksub.set(Ks,nksub*(dHat+1));
         lg= nvsub*Th.nt;
         KN<R> v(lg);
         for (int k=0,i=0;k<Th.nt;k++)
         {
             FElementL K=(*this)[k];
             for(int l=0;l<nvsub;l++)
                 v[i++] =   K(Psub[l], U, componante, op_U)  ;
             
         }
         return KN<R>(true,v);// to remove the copy.
     }
     
     
     
     /*
      template<class MMesh>
      KN<double>  GFESpace<MMesh>::newSaveDraw(const KN_<R> & U,int composante,int & lg,int & nsb) const
      {
      nsb = TFE[0]->nbsubdivision;
      int nsbv = NbOfSubInternalVertices(nsb,d);
      lg = nsbv*Th.nt;
      cout << "newSaveDraw what: nt " << Th.nt << " " << nsbv << " " << lg << endl;
      KN<double> v(lg);
      ffassert(v);
      for (int k=0,i=0;k<Th.nt;k++)
      {
      (*this)[k].SaveDraw( U,composante,&v[i]);
      i+=nsbv;
      }
      return KN<double>(true,v);// to remove the copy.
      }
      */
     // explicite instance..
     template class GTypeOfFESum<Mesh2>;
     template class GTypeOfFESum<Mesh3>;
     template class GTypeOfFESum<MeshS>;
     template class GTypeOfFESum<MeshL>;
     template class GFESpace<Mesh1>;
     template class GFESpace<Mesh2>;
     template class GFESpace<Mesh3>;
     template class GFESpace<MeshS>;
     template class GFESpace<MeshL>;
     
     template  KN<double>  GFESpace<MeshL>::newSaveDraw<double>(const KN_<double> & U,int componante,int & lg,KN<typename MeshL::RdHat> &Psub,KN<int> &Ksub,int op_U) const ;
     template  KN<double>  GFESpace<MeshS>::newSaveDraw<double>(const KN_<double> & U,int componante,int & lg,KN<typename MeshS::RdHat> &Psub,KN<int> &Ksub,int op_U) const ;
     template  KN<double>  GFESpace<Mesh3>::newSaveDraw<double>(const KN_<double> & U,int componante,int & lg,KN<typename Mesh3::RdHat> &Psub,KN<int> &Ksub,int op_U) const ;
     template  KN<double>  GFESpace<Mesh2>::newSaveDraw<double>(const KN_<double> & U,int componante,int & lg,KN<typename Mesh2::RdHat> &Psub,KN<int> &Ksub,int op_U) const ;
     template  KN<double>  GFESpace<Mesh1>::newSaveDraw<double>(const KN_<double> & U,int componante,int & lg,KN<typename Mesh1::RdHat> &Psub,KN<int> &Ksub,int op_U) const  ;
     
     typedef std::complex<double> Complex;
     template  KN<Complex>  GFESpace<MeshL>::newSaveDraw<Complex>(const KN_<Complex> & U,int componante,int & lg,KN<typename MeshL::RdHat> &Psub,KN<int> &Ksub,int op_U) const ;
     template  KN<Complex>  GFESpace<MeshS>::newSaveDraw<Complex>(const KN_<Complex> & U,int componante,int & lg,KN<typename MeshS::RdHat> &Psub,KN<int> &Ksub,int op_U) const ;
     template  KN<Complex>  GFESpace<Mesh3>::newSaveDraw<Complex>(const KN_<Complex> & U,int componante,int & lg,KN<typename Mesh3::RdHat> &Psub,KN<int> &Ksub,int op_U) const ;
     template  KN<Complex>  GFESpace<Mesh2>::newSaveDraw<Complex>(const KN_<Complex> & U,int componante,int & lg,KN<typename Mesh2::RdHat> &Psub,KN<int> &Ksub,int op_U) const ;
     template  KN<Complex>  GFESpace<Mesh1>::newSaveDraw<Complex>(const KN_<Complex> & U,int componante,int & lg,KN<typename Mesh1::RdHat> &Psub,KN<int> &Ksub,int op_U) const  ;
     
     
 }
