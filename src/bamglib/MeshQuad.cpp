// -*- Mode : c++ -*-
//
// SUMMARY  :      
// USAGE    :        
// ORG      : 
// AUTHOR   : Frederic Hecht
// E-MAIL   : hecht@ann.jussieu.fr
//
// ********** DO NOT REMOVE THIS BANNER **********
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

//
// SUMMARY: Bamg: Bidimensional Anisotrope Mesh Generator
// RELEASE: 0 
// AUTHOR:   F. Hecht,    
// ORG    :  UMPC
// E-MAIL :   Frederic.Hecht@Inria.fr   
//
// ORIG-DATE:     frev 98
//---------------------------------------------------------
//  to make quad 
// -------------------

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "Meshio.h"
#include "Mesh2.h"
#include "QuadTree.h"
#include "SetOfE4.h"

namespace bamg {

static const  Direction NoDirOfSearch=Direction();

#ifdef __MWERKS__
#pragma global_optimizer on
#pragma optimization_level 1
//#pragma inline_depth 0
#endif

class DoubleAndInt4 {
public: double q; Int4 i3j;
  int operator<(DoubleAndInt4 a){return q > a.q;}
};// to sort by decroissant



template<class T>  inline void  HeapSort(T *c,long n)
{
  long l,j,r,i;
  T crit;
  c--; // on decale de 1 pour que le tableau commence a 1
  if( n <= 1) return;
  l = n/2 + 1;
  r = n;
  while (1) { // label 2
    if(l <= 1 ) { // label 20
      crit = c[r];
      c[r--] = c[1];
    if ( r == 1 ) { c[1]=crit; return;}
    } else  crit = c[--l]; 
    j=l;
    while (1) {// label 4
      i=j;
      j=2*j;
      if  (j>r) {c[i]=crit;break;} // L8 -> G2
      if ((j<r) && (c[j] < c[j+1])) j++; // L5
      if (crit < c[j]) c[i]=c[j]; // L6+1 G4
      else {c[i]=crit;break;} //L8 -> G2
    }
  }
}

Triangle * swapTest(Triangle *t1,Int2 a);
// 


double QuadQuality(const Vertex & a,const Vertex &b,const Vertex &c,const Vertex &d)
{
  // calcul de 4 angles --
  R2 A((R2)a),B((R2)b),C((R2)c),D((R2)d);
  R2 AB(B-A),BC(C-B),CD(D-C),DA(A-D);
  //  Move(A),Line(B),Line(C),Line(D),Line(A);
  const Metric & Ma  = a;
  const Metric & Mb  = b;
  const Metric & Mc  = c;
  const Metric & Md  = d;
    
  double lAB=Norme2(AB);
  double lBC=Norme2(BC);
  double lCD=Norme2(CD);
  double lDA=Norme2(DA);
  AB /= lAB;  BC /= lBC;  CD /= lCD;  DA /= lDA;
  // version aniso 
  double cosDAB= Ma(DA,AB)/(Ma(DA)*Ma(AB)),sinDAB= Det(DA,AB);
  double cosABC= Mb(AB,BC)/(Mb(AB)*Mb(BC)),sinABC= Det(AB,BC);
  double cosBCD= Mc(BC,CD)/(Mc(BC)*Mc(CD)),sinBCD= Det(BC,CD);
  double cosCDA= Md(CD,DA)/(Md(CD)*Md(DA)),sinCDA= Det(CD,DA);
  double sinmin=Min(Min(sinDAB,sinABC),Min(sinBCD,sinCDA));
  // cout << A << B << C << D ;
  // cout << " sinmin " << sinmin << " " << cosDAB << " " << cosABC << " " << cosBCD << " " << cosCDA << endl;
  // rattente(1);
  if (sinmin<=0) return sinmin;
  return 1.0-Max(Max(Abs(cosDAB),Abs(cosABC)),Max(Abs(cosBCD),Abs(cosCDA)));
}

 GeometricalEdge *   Triangles::ProjectOnCurve( Edge & BhAB, Vertex &  vA, Vertex & vB,
						Real8 theta,
						Vertex & R,VertexOnEdge &  BR,VertexOnGeom & GR) 
                      
{
   void *pA=0,*pB=0;
   Real8 tA=0,tB=0;
   R2 A=vA,B=vB;
    Vertex * pvA=&vA, * pvB=&vB;
  if (vA.vint == IsVertexOnVertex)
    {
      //  cout << " Debut vertex = " << BTh.Number(vA.onbv) ;
      pA=vA.onbv;
    }
  else if (vA.vint == IsVertexOnEdge)
    {
      pA=vA.onbe->be;
      tA=vA.onbe->abcisse;
      // cout << " Debut edge = " << BTh.Number(vA.onbv) << " " << tA ;

    }
  else
    {cerr << "ProjectOnCurve On Vertex " << BTh.Number(vA) << " " << endl;
     cerr << " forget call to SetVertexFieldOnBTh" << endl;
     MeshError(-1);
    } 
    
  if (vB.vint == IsVertexOnVertex)
    {
      // cout << " Fin vertex = " << BTh.Number(vB.onbv) << endl;
      pB=vB.onbv;
    }
  else if(vB.vint == IsVertexOnEdge)
    {
      pB=vB.onbe->be;
      tB=vB.onbe->abcisse;
      // cout << " Fin edge = " << BTh.Number(vB.onbe->be) << " " <<  tB ;

    }
  else
    {cerr << "ProjectOnCurve On Vertex " << BTh.Number(vB) << " " << endl;
     cerr << " forget call to SetVertexFieldOnBTh" << endl;
     MeshError(-1);
    } 
  Edge * e = &BhAB;
 assert( pA && pB && e);
 // be carefull the back ground edge e is on same geom edge 
 // of the initiale edge def by the 2 vertex A B;
 assert(e>=BTh.edges && e<BTh.edges+BTh.nbe);// Is a background Mesh;   
// walk on BTh edge 
 //assert(0 /* not finish ProjectOnCurve with BackGround Mesh*/);
// 1 first find a back ground edge contening the vertex A
// 2 walk n back gound boundary to find the final vertex B


#ifdef DEBUG
// we suppose if the vertex A is on a background edge and
// the abcisse is 0 or 1 when the related point is not
// a end of curve <=>  !IsRequiredVertex
    if (vA.vint == IsVertexOnEdge)
      if (tA<=0)
	assert(! (*vA.onbe->be)[0].on->IsRequiredVertex());
      else if (tA>=1) 
      	assert(!(*vA.onbe->be)[1].on->IsRequiredVertex());
#endif

     if( vA.vint == IsVertexOnEdge) 
       { // find the start edge 
	     e = vA.onbe->be;	 

       } 
     else if (vB.vint == IsVertexOnEdge) 
       {
         theta = 1-theta;
	 Exchange(tA,tB);
	 Exchange(pA,pB);
	 Exchange(pvA,pvB);
	 Exchange(A,B);
	 e =  vB.onbe->be;

	 // cout << " EXCHANGE  A et B) " << endl;
       } 
     else
       { // do the search by walking 
	 assert(0 /* A FAIRE */);
       }

     // find the direction of walking with sens of edge and pA,PB;
   R2 AB=B-A;

   Real8 cosE01AB = (( (R2) (*e)[1] - (R2) (*e)[0] ) , AB);
   int kkk=0;
   int sens = (cosE01AB>0) ? 1 : 0;
 
   //   Real8 l=0; // length of the edge AB
   Real8 abscisse = -1;
 
   for (int cas=0;cas<2;cas++)
     {// 2 times algo:
       //    1 for computing the length l
       //    2 for find the vertex 
       int  iii;
        Vertex  *v0=pvA,*v1; 
       Edge *neee,*eee;
       Real8 lg =0; // length of the curve 
       Real8 te0;
       // we suppose take the curve's abcisse 
       // cout << kkk << " e = " << BTh.Number(e) << "  v0=  " 
       //    << BTh.Number(v0) << " v1 = " << BTh.Number((*e)[sens]) << endl;
       for ( eee=e,iii=sens,te0=tA;
             eee && ((( void*) eee) != pB) && (( void*) (v1=&((*eee)[iii]))) != pB ;
             neee = eee->adj[iii],iii = 1-neee->Intersection(*eee),eee = neee,v0=v1,te0=1-iii )   
	      { 
		//	cout << kkk << " eee = " << BTh.Number(eee) << "  v0=  " 
		//     << BTh.Number(v0) << " v1 = " << BTh.Number(v1) << endl;
		
		assert(kkk++<100);
		assert(eee);
		Real8 lg0 = lg;
		Real8 dp = LengthInterpole(v0->m,v1->m,(R2) *v1 - (R2) *v0);
	     	lg += dp;
		if (cas && abscisse <= lg)
		  { // ok we find the geom edge 
		    Real8 sss  =   (abscisse-lg0)/dp;
		    Real8 thetab = te0*(1-sss)+ sss*iii;
		    assert(thetab>=0 && thetab<=1);
		    BR = VertexOnEdge(&R,eee,thetab);

		    // cout << Number(R) << " = " <<  thetab << " on  " <<  BTh.Number(eee)
		    //	 << " = " << R << endl;

		    return  Gh.ProjectOnCurve(*eee,thetab,R,GR);

		  }
	      }
       // we find the end 
       if (v1 != pvB) 
	 {
	   if (( void*) v1 == pB)
	      tB = iii;
	     
	   Real8 lg0 = lg;
	   assert(eee);
	   v1 = pvB;
	   Real8 dp = LengthInterpole(v0->m,v1->m,(R2) *v1 - (R2) *v0);
	   lg += dp;	
	   abscisse = lg*theta;
	   if (abscisse <= lg && abscisse >= lg0 ) // small optimisation we know the lenght because end
	     { // ok we find the geom edge 
	       Real8 sss  =   (abscisse-lg0)/dp;
	       Real8 thetab = te0*(1-sss)+ sss*tB;
	       assert(thetab>=0 && thetab<=1);
	       BR = VertexOnEdge(&R,eee,thetab);
	      
	       //	cout << kkk << " eee = " << BTh.Number(eee) << "  v0=  " 
	       //     << BTh.Number(v0) << " " << te0
	       //      << " v1 = " << BTh.Number(v1) <<  " " << tB  << endl;

	       //out << Number(R) << " Opt  = " <<  thetab << " on  " <<  BTh.Number(eee) 
	       //	    << " = " << R << endl;

	       return  Gh.ProjectOnCurve(*eee,thetab,R,GR);
	     }
	 }
       abscisse = lg*theta;
       
     }
   cerr << " Big Bug" << endl;
   MeshError(678);
   return 0; // just for the compiler 
      
}                  


void Triangles::MakeQuadrangles(double costheta)
{
  if (verbosity>2) 
    cout << "  -- MakeQuadrangles costheta = " << costheta << endl;
  if (verbosity>5)  
    cout << "    (in)  Nb of Quadrilaterals = " << NbOfQuad 
	 << " Nb Of Triangles = " << nbt-NbOutT- NbOfQuad*2 
         << " Nb of outside triangles = " << NbOutT << endl;

  if (costheta >1) {
    if (verbosity>5)
      cout << "     do nothing costheta >1 "<< endl;
    return;}

  
  Int4 nbqq = (nbt*3)/2;
  DoubleAndInt4  *qq = new DoubleAndInt4[nbqq];

  Int4 i,ij;
  int j;
  Int4 k=0;
  for (i=0;i<nbt;i++)
    for (j=0;j<3;j++)
      if ((qq[k].q=triangles[i].QualityQuad(j))>=costheta)
 	   qq[k++].i3j=i*3+j;
//  sort  qq
  HeapSort(qq,k);
  
  Int4 kk=0;
  for (ij=0;ij<k;ij++)
    { 
      //      cout << qq[ij].q << " " << endl;
      i=qq[ij].i3j/3;
      j=(int) (qq[ij].i3j%3);
      // optisamition no float computation  
      if (triangles[i].QualityQuad(j,0) >=costheta) 
        triangles[i].SetHidden(j),kk++;
    }
  NbOfQuad = kk;
  if (verbosity>2) 
    {
    cout << "    (out)  Nb of Quadrilaterals = " << NbOfQuad 
	 << " Nb Of Triangles = " << nbt-NbOutT- NbOfQuad*2 
         << " Nb of outside triangles = " << NbOutT << endl;
    }
  delete [] qq;
#ifdef DRAWING2 
  Draw();
  inquire();
#endif

}
/*
Triangles::BThBoundary(Edge e,Real 8) const 
{
  // pointeur of the background must be on 
  // 
  Edge be = e.on;
}
*/
int  Triangles::SplitElement(int choice)
{
  Direction NoDirOfSearch;
  const  int withBackground = &BTh != this && &BTh;
 if (verbosity>2) 
   cout << "  -- SplitElement " << (choice? " Q->4Q and T->4T " : " Q->4Q or T->3Q " ) << endl;;
 if (verbosity>5)
   cout << endl << "    (in)  Nb of Quadrilaterals = " << NbOfQuad 
	<< " Nb Of Triangles = " << nbt-NbOutT- NbOfQuad*2 
	<< " Nb of outside triangles = " << NbOutT << endl;
 
 ReNumberingTheTriangleBySubDomain();
#ifdef DRAWING2 
  Draw();
 inquire();
#endif
 //int nswap =0;
  const Int4 nfortria( choice ? 4 : 6);
    if(withBackground) 
    {
      BTh.SetVertexFieldOn();
      SetVertexFieldOnBTh();
    }
   else
     BTh.SetVertexFieldOn();
   
  Int4 newnbt=0,newnbv=0;
  Int4 * kedge = 0;
 Int4 newNbOfQuad=NbOfQuad;
  Int4 * ksplit= 0, * ksplitarray=0;
  Int4 kkk=0;
  int ret =0;
  if (nbvx<nbv+nbe) return 1;//   
  Triangles *  OCurrentTh= CurrentTh;
  CurrentTh = this;
  // 1) create  the new points by spliting the internal edges 
  // set the 
  Int4 nbvold = nbv;
  Int4 nbtold = nbt;
  Int4 NbOutTold  = NbOutT;
  Int4  NbEdgeOnGeom=0;
  Int4 i;
  
  nbt = nbt - NbOutT; // remove all the  the ouside triangles 
  Int4 nbtsave = nbt;
  Triangle * lastT = triangles + nbt;
  for (i=0;i<nbe;i++)
    if(edges[i].on) NbEdgeOnGeom++;
  Int4 newnbe=nbe+nbe;
  //  Int4 newNbVerticesOnGeomVertex=NbVerticesOnGeomVertex;
  Int4 newNbVerticesOnGeomEdge=NbVerticesOnGeomEdge+NbEdgeOnGeom;
  // Int4 newNbVertexOnBThVertex=NbVertexOnBThVertex;
  Int4 newNbVertexOnBThEdge=withBackground ? NbVertexOnBThEdge+NbEdgeOnGeom:0;

  // do allocation for pointeur to the geometry and background
  VertexOnGeom * newVerticesOnGeomEdge = new VertexOnGeom[newNbVerticesOnGeomEdge];
  VertexOnEdge *newVertexOnBThEdge = newNbVertexOnBThEdge ?  new VertexOnEdge[newNbVertexOnBThEdge]:0;
  if (NbVerticesOnGeomEdge)
    memcpy(newVerticesOnGeomEdge,VerticesOnGeomEdge,sizeof(VertexOnGeom)*NbVerticesOnGeomEdge);
  if (NbVertexOnBThEdge)
    memcpy(newVertexOnBThEdge,VertexOnBThEdge,sizeof(VertexOnEdge)*NbVertexOnBThEdge);
  Edge *newedges = new Edge [newnbe];
  //  memcpy(newedges,edges,sizeof(Edge)*nbe);
  SetOfEdges4 * edge4= new SetOfEdges4(nbe,nbv);
#ifdef DEBUG
  for (i=0;i<nbt;i++)
    triangles[i].check();
#endif
#ifdef DRAWING2  
  reffecran();
#endif  
  Int4 k=nbv;
  Int4 kk=0;
  Int4 kvb = NbVertexOnBThEdge;
  Int4 kvg = NbVerticesOnGeomEdge;
  Int4 ie =0;
  Edge ** edgesGtoB=0;
  if (withBackground)
   edgesGtoB= BTh.MakeGeometricalEdgeToEdge();
  Int4 ferr=0;
  for (i=0;i<nbe;i++)
     newedges[ie].on=0;

  for (i=0;i<nbe;i++)
    {
      GeometricalEdge *ong =  edges[i].on;

      newedges[ie]=edges[i];
      newedges[ie].adj[0]=newedges+(edges[i].adj[0]-edges) ;
      newedges[ie].adj[1]=newedges + ie +1;
      R2 A = edges[i][0],B = edges[i][1];
      // cout << " ie = " << ie <<"  v0 = " <<  Number(newedges[ie][0]) << endl;


      kk += (i == edge4->addtrie(Number(edges[i][0]),Number(edges[i][1])));
      if (ong) // a geometrical edges 
	{ 
	  if (withBackground)
	    {
	      // walk on back ground mesh 
	      //  newVertexOnBThEdge[ibe++] = VertexOnEdge(vertices[k],bedge,absicsseonBedge); 
	      // a faire -- difficile 
	      // the first PB is to now a background edge between the 2 vertices
	      assert(edgesGtoB); 
	      // cout << " ie = " << ie <<"  v0 = " <<  Number(newedges[ie][0]) << endl;
             ong= ProjectOnCurve(*edgesGtoB[Gh.Number(edges[i].on)],
			     edges[i][0],edges[i][1],0.5,vertices[k],
                             newVertexOnBThEdge[kvb],
                             newVerticesOnGeomEdge[kvg++]);
	     vertices[k].ReferenceNumber= edges[i].ref;
	     vertices[k].DirOfSearch =   NoDirOfSearch;        
;
	     // get the Info on background mesh 
	     Real8 s =        newVertexOnBThEdge[kvb];
	     Vertex &  bv0  = newVertexOnBThEdge[kvb][0];
	     Vertex &  bv1  = newVertexOnBThEdge[kvb][1];
	     // compute the metrix of the new points 
	     vertices[k].m =  Metric(1-s,bv0,s,bv1); 
	     kvb++;
	     // cout << " ie = " << ie <<"  v0 = " <<  Number(newedges[ie][0]) << endl;
	    }
	  else 
	    {
	      ong=Gh.ProjectOnCurve(edges[i],
				    0.5,vertices[k],newVerticesOnGeomEdge[kvg++]);
	      // vertices[k].i = toI2( vertices[k].r);
	      vertices[k].ReferenceNumber = edges[i].ref;
	      vertices[k].DirOfSearch = NoDirOfSearch;
	      vertices[k].m =  Metric(0.5,edges[i][0],0.5,edges[i][1]);	      
	    }  
	}
      else // straigth line edge ---
	{ 
	  vertices[k].r = ((R2) edges[i][0] + (R2)  edges[i][1] )*0.5;
	  vertices[k].m =  Metric(0.5,edges[i][0],0.5,edges[i][1]);
	  vertices[k].on = 0;
	}
      //vertices[k].i = toI2( vertices[k].r);
      R2 AB =  vertices[k].r;
      R2 AA = (A+AB)*0.5;
      R2 BB = (AB+B)*0.5;
      vertices[k].ReferenceNumber = edges[i].ref;
	    vertices[k].DirOfSearch = NoDirOfSearch;
	    
      newedges[ie].on = Gh.Contening(AA,ong);
      newedges[ie++].v[1]=vertices+k;

      newedges[ie]=edges[i];
      newedges[ie].adj[0]=newedges + ie -1;
      newedges[ie].adj[1]=newedges+(edges[i].adj[1]-edges) ;
      newedges[ie].on =  Gh.Contening(BB,ong);
      newedges[ie++].v[0]=vertices+k;
      // cout << " ie = " << ie-2 << " vm " << k << " v0 = " <<  Number(newedges[ie-2][0])
      //	   << " v1 = " << Number(newedges[ie-1][1])  
      //	   << " ong =" << ong-Gh.edges 
      //	   << " on 0 =" <<  newedges[ie-2].on -Gh.edges << AA
      //	   << " on 1 =" <<  newedges[ie-1].on -Gh.edges << BB 
      //	   << endl;
      k++;
    }
#ifdef DEBUG
  assert(kvb ==  newNbVertexOnBThEdge);
  // verif edge 
  { Vertex *v0 = vertices, *v1 = vertices+ k;
    for (Int4  i=0;i<ie;i++)
     {
       assert( &newedges[i][0] >= v0 &&  &newedges[i][0] < v1);
       assert( &newedges[i][1] >= v0 &&  &newedges[i][1] < v1);
     }
  }
#endif
  if (edgesGtoB) delete [] edgesGtoB;
  edgesGtoB=0;

  newnbv=k;
  newNbVerticesOnGeomEdge=kvg;
  if (newnbv> nbvx) goto Error;// bug 
    
  nbv = k;


  kedge = new Int4[3*nbt+1];
  ksplitarray = new Int4[nbt+1];
  ksplit = ksplitarray +1; // because ksplit[-1] == ksplitarray[0]
 
  for (i=0;i<3*nbt;i++)
    kedge[i]=-1;

  //  

 for (i=0;i<nbt;i++)
   {  

     Triangle & t = triangles[i];
     assert(t.link);
     for(int j=0;j<3;j++)
       {
	 const TriangleAdjacent ta = t.Adj(j);
	 const Triangle & tt = ta;
	 if (&tt >= lastT)
	   t.SetAdj2(j,0,0);// unset adj
	 const Vertex & v0 = t[VerticesOfTriangularEdge[j][0]];
	 const Vertex & v1 = t[VerticesOfTriangularEdge[j][1]];
	 Int4  ke =edge4->findtrie(Number(v0),Number(v1));
	 if (ke>0) 
	   {
	     Int4 ii = Number(tt);
	     int  jj = ta;
	     Int4 ks = ke + nbvold;
	     kedge[3*i+j] = ks;
	     if (ii<nbt) // good triangle
	       kedge[3*ii+jj] = ks;
	     Vertex &A=vertices[ks];
	     Real8 aa,bb,cc,dd;
	     if ((dd=Area2(v0.r,v1.r,A.r)) >=0)
	       { // warning PB roundoff error 
		 if (t.link && ( (aa=Area2( A.r    , t[1].r , t[2].r )) < 0.0 
				||   (bb=Area2( t[0].r , A.r    , t[2].r )) < 0.0  
				||   (cc=Area2( t[0].r , t[1].r , A.r    )) < 0.0))
		   ferr++, cerr << " Error : " <<  ke + nbvold << " not in triangle " 
				<< i << " In=" << !!t.link
				<< " " <<  aa  << " " << bb << " " << cc << " " << dd << endl;
		 
	       }
	     
	     else
	       {
		 if (tt.link && ( (aa=Area2( A.r     , tt[1].r , tt[2].r )) < 0 
				 ||   (bb=Area2( tt[0].r , A.r     , tt[2].r )) < 0 
				 ||   (cc=Area2( tt[0].r , tt[1].r , A.r     )) < 0)) 
		   ferr++, cerr << " Warning : " <<  ke + nbvold << " not in triangle " << ii 
				<< " In=" << !!tt.link 
				<< " " <<  aa  << " " << bb << " " << cc << " " << dd << endl;
				
	       } 
	     
	   }
       }
   }
  if(ferr)
    {
      cerr << " Number of triangles with P2 interpolation Probleme " << ferr << endl;;
      MeshError(9);
    }

  for (i=0;i<nbt;i++)
    {
      ksplit[i]=1; // no split by default
      const Triangle & t = triangles[ i];
      // cout << " Triangle " << i << " " << t  << !!t.link << ":: " ;
      int nbsplitedge =0;
      int nbinvisible =0;
      int invisibleedge=0;
      int kkk[3];      
      for (int j=0;j<3;j++)
	{
	  if (t.Hidden(j)) invisibleedge=j,nbinvisible++;
	  
	  const TriangleAdjacent ta = t.Adj(j);
	  const Triangle & tt = ta;

	  
	  const Vertex & v0 = t[VerticesOfTriangularEdge[j][0]];
	  const Vertex & v1 = t[VerticesOfTriangularEdge[j][1]];
	 //  cout << " ke = " << kedge[3*i+j]  << " " << Number(v0) << " " << Number(v1) << "/ ";
	  if ( kedge[3*i+j] < 0) 
	    {
	      Int4  ke =edge4->findtrie(Number(v0),Number(v1));
	      //  cout << ":" << ke << "," << !!t.link << " " <<  &tt ;
	      if (ke<0) // new 
		{
		  if (&tt) // internal triangles all the boundary 
		      { // new internal edges 
			Int4 ii = Number(tt);
			int  jj = ta;
	      
			kedge[3*i+j]=k;// save the vertex number 
			kedge[3*ii+jj]=k;
			if (k<nbvx) 
			  {
			    vertices[k].r = ((R2) v0+(R2) v1 )/2;
			    //vertices[k].i = toI2( vertices[k].r);
			    vertices[k].ReferenceNumber=0;
	        vertices[k].DirOfSearch =NoDirOfSearch;
			    vertices[k].m =  Metric(0.5,v0,0.5,v1);
			  }
			k++;
			kkk[nbsplitedge++]=j;		      
		    } // tt 
		  else
		    cerr <<endl <<  " Bug " <<i<< " " << j << " t=" << t << endl;
		  
		} // ke<0	       
	      else
		{ // ke >=0
		  kedge[3*i+j]=nbvold+ke;
		  kkk[nbsplitedge++]=j;// previously splited
		}
	    }
	  else 
	    kkk[nbsplitedge++]=j;// previously splited
	  
	} 
      assert (nbinvisible<2);
     // cout << " " <<  nbinvisible << " " <<  nbsplitedge << endl;
      switch (nbsplitedge) {
      case 0: ksplit[i]=10; newnbt++; break;   // nosplit
      case 1: ksplit[i]=20+kkk[0];newnbt += 2; break; // split in 2 
      case 2: ksplit[i]=30+3-kkk[0]-kkk[1];newnbt += 3; break; // split in 3 
      case 3:
	if (nbinvisible) ksplit[i]=40+invisibleedge,newnbt += 4;
	else   ksplit[i]=10*nfortria,newnbt+=nfortria;
	break;
      } 
    assert(ksplit[i]>=40);
    }
  //  now do the element split
  newNbOfQuad = 4*NbOfQuad;
  nbv = k;
#ifdef DRAWING2  
  inquire();
#endif  
//  cout << " Nbv = " << nbv << endl;
  kkk = nbt;
  ksplit[-1] = nbt;
  // look on  old true  triangles 

  for (i=0;i<nbtsave;i++)
    {
      //     cout << "Triangle " << i << " " << ksplit[i] << ":" << triangles[i]
      //	   << "  ----------------------------------------------- " <<endl;
      // Triangle * tc=0;
      int  nbmkadj=0;
      Int4 mkadj [100];
      mkadj[0]=i;
      Int4 kk=ksplit[i]/10;
      int  ke=(int) (ksplit[i]%10);
      assert(kk<7 && kk >0);
      
      // def the numbering   k (edge) i vertex 
      int k0 = ke;
      int k1 = NextEdge[k0];
      int k2 = PreviousEdge[k0];
      int i0 = OppositeVertex[k0];
      int i1 = OppositeVertex[k1];
      int i2 = OppositeVertex[k2];
       
       Triangle &t0=triangles[i];
       Vertex * v0=t0(i0);           
       Vertex * v1=t0(i1);           
       Vertex * v2=t0(i2);

       // cout << "nbmkadj " << nbmkadj << " it=" << i <<endl;
       assert(nbmkadj< 10);
       // --------------------------
       TriangleAdjacent ta0(t0.Adj(i0)),ta1(t0.Adj(i1)),ta2(t0.Adj(i2));
       // save the flag Hidden
       int hid[]={t0.Hidden(0),t0.Hidden(1),t0.Hidden(2)};
       // un set all adj -- save Hidden flag --
       t0.SetAdj2(0,0,hid[0]);
       t0.SetAdj2(1,0,hid[1]);
       t0.SetAdj2(2,0,hid[2]);
       // --  remake 
       switch  (kk) {
       case 1: break;// nothing 
       case 2: // 
	 {
	   Triangle &t1=triangles[kkk++];
	   t1=t0;
	   assert (kedge[3*i+i0]>=0);
	   Vertex * v3 = vertices + kedge[3*i+k0];
	   
	   t0(i2) = v3;
	   t1(i1) = v3;
	   t0.SetAllFlag(k2,0);
	   t1.SetAllFlag(k1,0);
	 } 
	 break; 
       case 3: //
	 {
	   Triangle &t1=triangles[kkk++];
            Triangle &t2=triangles[kkk++];
            t2=t1=t0;
            assert (kedge[3*i+k1]>=0);
            assert (kedge[3*i+k2]>=0);
            
            Vertex * v01 = vertices + kedge[3*i+k2];
            Vertex * v02 = vertices + kedge[3*i+k1]; 
            t0(i1) = v01; 
            t0(i2) = v02; 
            t1(i2) = v02;
            t1(i0) = v01; 
            t2(i0) = v02; 
	    t0.SetAllFlag(k0,0);
	    t1.SetAllFlag(k1,0);
	    t1.SetAllFlag(k0,0);
	    t2.SetAllFlag(k2,0);
	 } 
	 break;
       case 4: // 
       case 6: // split in 4 
	 {
	   Triangle &t1=triangles[kkk++];
	   Triangle &t2=triangles[kkk++];
	   Triangle &t3=triangles[kkk++];
	   t3=t2=t1=t0;
	   assert(kedge[3*i+k0] >=0 && kedge[3*i+k1] >=0 && kedge[3*i+k2] >=0);
	   Vertex * v12 = vertices + kedge[3*i+k0];
	   Vertex * v02 = vertices + kedge[3*i+k1]; 
	   Vertex * v01 = vertices + kedge[3*i+k2];
	   // cout << Number(t0(i0))  << " " << Number(t0(i1)) 
	   //     << " " <<  Number(t0(i2)) 
	   //     << " " <<  kedge[3*i+k0] 
	   //     << " " <<  kedge[3*i+k1] 
	   //     << " " <<  kedge[3*i+k2] << endl;
	   t0(i1) = v01;
	   t0(i2) = v02;
	   t0.SetAllFlag(k0,hid[k0]);

	   t1(i0) = v01;
	   t1(i2) = v12;
	   t0.SetAllFlag(k1,hid[k1]);
	   
	   t2(i0) = v02;
	   t2(i1) = v12;
	   t2.SetAllFlag(k2,hid[k2]);
	   
	   t3(i0) = v12;
	   t3(i1) = v02;
	   t3(i2) = v01;
	   	   
	   t3.SetAllFlag(0,hid[0]);	   
	   t3.SetAllFlag(1,hid[1]);	   
	   t3.SetAllFlag(2,hid[2]);

	   if ( kk == 6)
	     {
	       
	       Triangle &t4=triangles[kkk++];
	       Triangle &t5=triangles[kkk++];
	       
	       t4 = t3;
	       t5 = t3;

	       t0.SetHidden(k0);
	       t1.SetHidden(k1);
	       t2.SetHidden(k2);
	       t3.SetHidden(0);
	       t4.SetHidden(1);
	       t5.SetHidden(2);
	       
		if (nbv < nbvx ) 
		  {
		    vertices[nbv].r = ((R2) *v01 + (R2) *v12  + (R2) *v02 ) / 3.0;
		    vertices[nbv].ReferenceNumber =0;
	      vertices[nbv].DirOfSearch =NoDirOfSearch;
		    //vertices[nbv].i = toI2(vertices[nbv].r);
		    Real8 a3[]={1./3.,1./3.,1./3.};
		    vertices[nbv].m = Metric(a3,v0->m,v1->m,v2->m);
		    Vertex * vc =  vertices +nbv++;
		    t3(i0) = vc;
		    t4(i1) = vc;
		    t5(i2) = vc;

		  }
		else
		  goto Error; 
	     }
		    
	 } 
	 break;         
       }
 
       // cout << "  -- " << i << " " << nbmkadj << " " << kkk << " " << tc << endl;
       //  t0.SetDetf();
       // save all the new triangles
       mkadj[nbmkadj++]=i;
       Int4 jj;
       if (t0.link) 
	 for (jj=nbt;jj<kkk;jj++)
	   {
	     triangles[jj].link=t0.link;
	     t0.link= triangles+jj;
	     mkadj[nbmkadj++]=jj;
	     // triangles[jj].SetDet();
	   }
       // cout << "  -- " << i << " " << nbmkadj << endl;
       assert(nbmkadj<=13);// 13 = 6 + 4 + 3
  
       if (kk==6)  newNbOfQuad+=3;
	 //	 triangles[i].Draw();       

       for (jj=ksplit[i-1];jj<kkk;jj++)
	 // triangles[jj].SetDet();
       //	   triangles[jj].Draw();
	 


       nbt = kkk;
       ksplit[i]= nbt; // save last adresse of the new triangles
       kkk = nbt;
       
    }

//  cout << " nv = " << nbv << " nbt = " << nbt << endl;
  for (i=0;i<nbv;i++)
    vertices[i].m = vertices[i].m*2.;
  //
  if(withBackground)
    for (i=0;i<BTh.nbv;i++)
      BTh.vertices[i].m =  BTh.vertices[i].m*2.;
#ifdef DRAWING2
  Draw();
  inquire();
#endif
  
  
  ret = 2;
  if (nbt>= nbtx) goto Error; // bug 
  if (nbv>= nbvx) goto Error; // bug 
  // generation of the new triangles 

  SetIntCoor("In SplitElement"); 

  ReMakeTriangleContainingTheVertex();
  if(withBackground)
    BTh.ReMakeTriangleContainingTheVertex();

  delete [] edges;
  edges = newedges;
  nbe = newnbe;
  NbOfQuad = newNbOfQuad;

  for (i=0;i<NbSubDomains;i++)
    { 
      Int4 k = subdomains[i].edge- edges;
      subdomains[i].edge =  edges+2*k; // spilt all edge in 2 
    }
    
  if (ksplitarray) delete [] ksplitarray;
  if (kedge) delete [] kedge;
  if (edge4) delete edge4;
  if (VerticesOnGeomEdge) delete [] VerticesOnGeomEdge;
  VerticesOnGeomEdge= newVerticesOnGeomEdge;
  if(VertexOnBThEdge) delete []  VertexOnBThEdge;
  VertexOnBThEdge = newVertexOnBThEdge;
  NbVerticesOnGeomEdge = newNbVerticesOnGeomEdge;
  NbVertexOnBThEdge=newNbVertexOnBThEdge;
  //  ReMakeTriangleContainingTheVertex();

  FillHoleInMesh();

#ifdef DEBUG
  for (i=0;i<nbt;i++)
    triangles[i].check();
#endif
  
 if (verbosity>2)
   cout << "    (out) Nb of Quadrilaterals = " << NbOfQuad 
	<< " Nb Of Triangles = " << nbt-NbOutT- NbOfQuad*2 
	<< " Nb of outside triangles = " << NbOutT << endl;

 CurrentTh=OCurrentTh;
 return 0; //ok

 Error:
  nbv = nbvold;
  nbt = nbtold;
  NbOutT = NbOutTold;
  // cleaning memory ---
  delete newedges;
  if (ksplitarray) delete [] ksplitarray;
  if (kedge) delete [] kedge;
  if (newVerticesOnGeomEdge) delete [] newVerticesOnGeomEdge;
  if (edge4) delete edge4;
  if(newVertexOnBThEdge) delete []  newVertexOnBThEdge;
  

  CurrentTh= OCurrentTh;
  return ret; // ok 
}

} //  end of namespcae bamg 
