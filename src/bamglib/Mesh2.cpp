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
#ifdef __MWERKS__
#ifdef __INTEL__
//#pragma global_optimizer off
//#pragma inline_depth(0)
//#pragma optimization_level 2
#endif
//#pragma inline_depth 0
#endif
extern bool withrgraphique;
#include <stdio.h>
#include <string.h>
#include <math.h> 
#include <time.h>
#include <iostream>
using namespace std; 

#include "Mesh2.h"
#include "QuadTree.h"
#include "SetOfE4.h"

namespace bamg {


#ifdef DEBUG1
extern int SHOW ; // for debugging 
int SHOW = 0; // for debugging 

#endif

int  Triangles::counter = 0;

Triangles * CurrentTh =0;

int hinterpole=1;


long NbUnSwap =0;
int  ForDebugging = 0;
const Direction NoDirOfSearch = Direction();
#ifndef NDEBUG 
inline void MyAssert(int i,char*ex,char * file,long line) 
{
  if( i) {
    cerr << "Error Assert:" << ex << " in " << file << " line: " << line << endl;
#ifdef  NOTFREEFEM
    exit(1); 
#else
    throw(ErrorExec("exit",1000));
#endif
  }
}
#endif

Int4 AGoodNumberPrimeWith(Int4 n)
{
  const Int4 BigPrimeNumber[] ={ 567890359L,
				 567890431L,  567890437L,  567890461L,  567890471L,
				 567890483L,  567890489L,  567890497L,  567890507L,
				 567890591L,  567890599L,  567890621L,  567890629L , 0};
  
  Int4 o = 0;
  Int4 pi = BigPrimeNumber[1];
  for (int i=0; BigPrimeNumber[i]; i++) {
    Int4 r = BigPrimeNumber[i] % n;
    Int4 oo = Min(Min(r,n-r),Min(Abs(n-2*r),Abs(n-3*r)));
    if ( o < oo) 
      o=oo,pi=BigPrimeNumber[i];}
  //  cout << " AGoodNumberPrimeWith " << n << " " <<pi << " "<< o << endl;
  return pi; 
}

class Triangles;
void MeshError(int Err,Triangles *Th){ 
 cerr << " Fatal error in the meshgenerator " << Err << endl ;
#ifdef  NOTFREEFEM
    exit(1); 
#else
  throw(ErrorMesh("Bamg",Err,Th));
#endif
}

 ostream& operator <<(ostream& f, const  Triangle & ta)
  {
    if(CurrentTh)
      f << "[" << CurrentTh->Number(ta) << "::" 
     <<  CurrentTh->Number(ta.ns[0]) << "," 
     <<  CurrentTh->Number(ta.ns[1]) << "," 
     <<  CurrentTh->Number(ta.ns[2]) << "," 
     << "{" <<  CurrentTh->Number(ta.at[0]) << " " << ta.aa[0] << "} " 
     << "{" <<  CurrentTh->Number(ta.at[1]) << " " << ta.aa[1] << "} " 
     << "{" <<  CurrentTh->Number(ta.at[2]) << " " << ta.aa[2] << "} " 
     << "]" ;
     else
       f << "[" 
     << ta.ns[0] << "," 
     << ta.ns[1] << "," 
     << ta.ns[2] << "," 
     << "{" << ta.at[0] << " " << ta.aa[0] << "} " 
     << "{" << ta.at[1] << " " << ta.aa[1] << "} " 
     << "{" << ta.at[2] << " " << ta.aa[2] << "} " 
     << "]" ;
   return f;}

void  swap(Triangle *t1,Int1 a1,
                 Triangle *t2,Int1 a2,
                 Vertex *s1,Vertex *s2,Icoor2 det1,Icoor2 det2)
{ // swap 
  // --------------------------------------------------------------
  // Int1 a2=aa[a];// les 2 numero de l arete dans les 2 triangles
  //                               
  //               sb                     sb    
  //             / | \                   /   \                      !
  //         as1/  |  \                 /a2   \                     !
  //           /   |   \               /    t2 \                    !
  //       s1 /t1  | t2 \s2  -->   s1 /___as2___\s2                 !
  //          \  a1|a2  /             \   as1   /  
  //           \   |   /               \ t1    /   
  //            \  |  / as2             \   a1/    
  //             \ | /                   \   /     
  //              sa                       sa   
  //  -------------------------------------------------------------
  int as1 = NextEdge[a1];
  int as2 = NextEdge[a2];
  int ap1 = PreviousEdge[a1];
  int ap2 = PreviousEdge[a2];
#ifdef DRAWING1
  couleur(0);
  t1->Draw();
  t2->Draw();
#endif
#ifdef DEBUG1
  t1->check();
  t2->check();
#endif
  (*t1)(VerticesOfTriangularEdge[a1][1]) = s2 ; // avant sb
  (*t2)(VerticesOfTriangularEdge[a2][1]) = s1  ; // avant sa
  // mise a jour des 2 adjacences externes 
  TriangleAdjacent taas1 = t1->Adj(as1),
    taas2 = t2->Adj(as2),
    tas1(t1,as1), tas2(t2,as2),
    ta1(t1,a1),ta2(t2,a2);
#ifdef DEBUG
  assert( ! ta1.Locked());
  assert( ! ta2.Locked());
#endif
  // externe haut gauche
  taas1.SetAdj2(ta2, taas1.GetAllFlag_UnSwap());
   // externe bas droite
  taas2.SetAdj2(ta1, taas2.GetAllFlag_UnSwap());
  // remove the Mark  UnMarkSwap 
  t1->SetUnMarkUnSwap(ap1);
  t2->SetUnMarkUnSwap(ap2);
  // interne 
  tas1.SetAdj2(tas2);

  t1->det = det1;
  t2->det = det2;

  t1->SetTriangleContainingTheVertex();
  t2->SetTriangleContainingTheVertex();
#ifdef DEBUG1
  t1->check();
  t2->check();
#endif
#ifdef DRAWING1 
  couleur(1);
  t1->Draw();
  t2->Draw();
#endif 
#ifdef DRAWING1
  if(  CurrentTh)
    CurrentTh->inquire();
#endif

} // end swap 





Int4 FindTriangle(Triangles &Th, Real8 x, Real8 y, double* a,int & inside)
 {
   CurrentTh=&Th;
   assert(&Th);
   I2 I = Th.toI2(R2(Min(Max(Th.pmin.x,x),Th.pmax.x),Min(Max(Th.pmin.y,y),Th.pmax.y))); 
   Icoor2 dete[3];
   Triangle & tb = *Th.FindTriangleContening(I,dete);
   
   if  (tb.link) 
     { // internal point in a true triangles
       a[0]= (Real8) dete[0]/ tb.det;
       a[1]= (Real8) dete[1] / tb.det;
       a[2] = (Real8) dete[2] / tb.det;
	 inside = 1;	 
	 return Th.Number(tb);
     } 
   else 
     {
       inside = 0; 
       double aa,bb;
       TriangleAdjacent  ta=CloseBoundaryEdgeV2(I,&tb,aa,bb);	 
       int k = ta;
       Triangle * tc = ta;
       if (!tc->link) 
	 { ta = ta.Adj();
	 tc=ta;
	 k = ta;
	 Exchange(aa,bb);
	 assert(tc->link);
	 }
       a[VerticesOfTriangularEdge[k][0]] = aa;
       a[VerticesOfTriangularEdge[k][1]] = bb;
       a[OppositeVertex[k]] = 1- aa -bb;
       return Th.Number(tc);
     }
 }


TriangleAdjacent CloseBoundaryEdge(I2 A,Triangle *t, double &a,double &b) {
// 
  //  cout << " - ";   	 
  int k=(*t)(0) ?  ((  (*t)(1) ? ( (*t)(2) ? -1 : 2) : 1  )) : 0;
  int dir=0;
  assert(k>=0);
  int kkk=0;  
  Icoor2 IJ_IA,IJ_AJ;
  TriangleAdjacent edge(t,OppositeEdge[k]);          
  for (;;edge = dir >0 ? Next(Adj(Next(edge))) : Previous(Adj(Previous(edge)))) 
   {  
   
    assert(kkk++<1000);      
    Vertex  &vI =  *edge.EdgeVertex(0);
    Vertex  &vJ =  *edge.EdgeVertex(1);
    I2 I=vI, J=vJ, IJ= J-I;
    IJ_IA = (IJ ,(A-I));
    //   cout << A << vI.i << vJ.i << edge << " " <<  IJ_IA << " dir " << dir <<endl;
    if (IJ_IA<0) {
     if (dir>0) {a=1;b=0;return edge;}// change of signe => I
     else {dir=-1;
        continue;}};// go in direction i 
    IJ_AJ = (IJ ,(J-A));
    if (IJ_AJ<0) {
    if(dir<0)  {a=0;b=1;return edge;}            
    else {dir = 1;
       continue;}}// go in direction j
    double IJ2 = IJ_IA + IJ_AJ;
    assert(IJ2);
    a= IJ_AJ/IJ2;
    b= IJ_IA/IJ2;
    //    cout<< "CloseBoundaryEdge a = " << a << " b= " << b << endl;
    return edge;
  } 
}

TriangleAdjacent Triangle::FindBoundaryEdge(int i) const
{
  // turn around  the vertex ns[i] also call  s
#ifdef DEBUG
  register Vertex * s  =  ns[i];
#endif
  Triangle   *t = (Triangle *) this , *ttc;
  int k=0,j = EdgesVertexTriangle[i][0],jc;
  int exterieur  = !link  ;
  
  do 
    {
      int exterieurp = exterieur;
      k++; 
#ifdef DEBUG
      assert( s == & (*t)[VerticesOfTriangularEdge[j][1]] );
#endif
      ttc =  t->at[j];
      exterieur = !ttc->link;
      if (exterieur+exterieurp == 1) 
	return TriangleAdjacent(t,j);
      jc = NextEdge[t->aa[j]&3];
      t = ttc;
      j = NextEdge[jc];
      assert(k<2000);
    } while ( (this!= t)); 
  return TriangleAdjacent(0,0);
 
}


TriangleAdjacent CloseBoundaryEdgeV2(I2 C,Triangle *t, double &a,double &b) 
{ 
 // walk around the vertex 
 // version 2 for remove the probleme if we fill the hole
  //int bug=1;
  //  Triangle *torigine = t;
  // restart:
  //   int dir=0;
  assert(t->link == 0);
  // to have a starting edges 
  // try the 3 edge bourna-- in case of internal hole 
  // and choice  the best 
  // 
  // 
  // the probleme is in case of  the fine and long internal hole
  // for exemple neart the training edge of a wing
  // 
  Vertex * s=0,*s1=0, *s0=0;
  Icoor2 imax = MaxICoor22;
  Icoor2 l0 = imax,l1 = imax;
  double dd2 =  imax;// infinity
  TriangleAdjacent er; 
  int  cas=-2;
  for (int j=0;j<3;j++)
    { 
      TriangleAdjacent ta=t->FindBoundaryEdge(j);
      if  (! (Triangle *) ta) continue;
      s0 = ta.EdgeVertex(0);
      s1 = ta.EdgeVertex(1);
      I2 A = * s0;
      I2 B = *ta.EdgeVertex(1);
      I2 AB = B-A,AC=C-A,BC=B-C;
      Icoor2  ACAC = (AC,AC), BCBC = (BC,BC);
      Icoor2  AB2  =   Norme2_2(AB); //  ||AB||^2
      Icoor2  ABAC  =   (AB,AC);         //  AB.AC|
      
      double d2;
      if ( ABAC < 0 )   // DIST A
        {
           if ( (d2=(double) ACAC)  <  dd2) 
             {
	       //  cout << " A "  << d2  << " " <<  dd2;
               er = ta;
               l0 = ACAC;
               l1 = BCBC;
               cas = 0;
               s = s0;
             }
        }
      else if (ABAC > AB2)  // DIST B
        {
           if ( (d2=(double) BCBC)  <  dd2) 
             {
	       // cout << " B "  << d2  << " " <<  dd2;
               dd2 = d2;
               er = Adj(ta); // other direction
               l0 = BCBC;
               l1 = ACAC;
               cas = 1;
               s = s1;
             }
        }
      else  // DIST AB
        { 

          double det_2 =  (double) Det(AB,AC); 
          det_2 *= det_2; // square of area*2 of triangle ABC
          d2 = det_2/ (double) AB2; // hauteur^2 in C of of triangle ABC      
	  //	  cout << " AB " << d2 << " " << dd2 
	  //      << " " << CurrentTh->Number(ta.EdgeVertex(0)) 
	  //     << " " << CurrentTh->Number(ta.EdgeVertex(1)) << " " ;

          if (d2 < dd2) 
	       {
	         dd2 = d2;
	         er = ta;
	         l0 = (AC,AC);
	         l1 = (BC,BC);
	         s = 0;
                 cas = -1;
		 //	 cout << " ABAC " <<  ABAC << " ABAC " << ABAC
		 //	      << " AB2 " << AB2 << endl;
		 b = ((double) ABAC/(double) AB2);
		 a = 1 - b;
	       }
        }
     }
   assert(cas !=-2);
   // l1 = ||C s1||  , l0 = ||C s0||
   // where s0,s1 are the vertex of the edge er

   if ( s) 
     { 
       t=er;
       TriangleAdjacent edge(er); 
       
       int kkk=0;  
       int linkp = t->link == 0;
       
       Triangle * tt=t=edge=Adj(Previous(edge));
       //  cout << CurrentTh->Number(t) << " " << linkp << endl;
       do  {  // loop around vertex s
	 
	 assert(edge.EdgeVertex(0)==s && kkk++<10000);
	 
	 int link = tt->link == 0;
	 //	 cout << CurrentTh->Number(tt) << " " << link << " " << CurrentTh->Number(s) 
	 //	      << " " << CurrentTh->Number(er.EdgeVertex(0)) 
	 //	      << " " << CurrentTh->Number(er.EdgeVertex(1)) 
	 //	      << " " << CurrentTh->Number(edge.EdgeVertex(0)) 
	 //	      << " " << CurrentTh->Number(edge.EdgeVertex(1)) 
	 //	      <<  endl;
	 if ((link + linkp) == 1) 
	   { // a boundary edge 
	     Vertex * st = edge.EdgeVertex(1);
	     I2 I=*st;
	     Icoor2  ll = Norme2_2 (C-I);
	     if (ll < l1) {  // the other vertex is neart 
	       s1=st;
	       l1=ll;
	       er = edge;
	       if(ll<l0) { // change of direction --
		 s1=s;
		 l1=l0;
		 s=st;
		 l0=ll;
		 t=tt;
		 edge=Adj(edge);
		 link=linkp;
		 er = edge;
	       }
	     }
	   }
	 
	 linkp=link;
	 edge=Adj(Previous(edge));
	 tt = edge;
       } while (t!=tt);

       assert((Triangle *) er);
       I2 A((I2)*er.EdgeVertex(0));
       I2 B((I2)*er.EdgeVertex(1));
       I2 AB=B-A,AC=C-A,CB=B-C;
       double aa =  (double) (AB,AC);
       double bb =  (double) (AB,CB);
       //  cout << " " << aa << " " << bb 
       //    << " " << CurrentTh->Number(er.EdgeVertex(0)) 
       //	    << " " << CurrentTh->Number(er.EdgeVertex(1)) ;
       if (aa<0)       a=1,b=0;
       else if(bb<0)   a=0,b=1;
       else  
	 {
	   a  = bb/(aa+bb);
	   b  = aa/(aa+bb);
	 }
     }
   
   //   cout <<" return= " <<  CurrentTh->Number(er.EdgeVertex(0)) << " " 
   //	<<  CurrentTh->Number(er.EdgeVertex(1)) << " " << a 
   //	<< " " << b <<" " << l0 << " " <<l1 <<endl;
   return er;
} 



Metric Triangles::MetricAt  (const R2 & A) const
  { //if ((vertices <= &v) && (vertices < v+nbv)) return v.m;
    I2 a = toI2(A);
    Icoor2 deta[3];
    Triangle * t =FindTriangleContening(a,deta);
    if (t->det <0) { // outside
      double ba,bb;
      TriangleAdjacent edge= CloseBoundaryEdge(a,t,ba,bb) ;
      return Metric(ba,*edge.EdgeVertex(0),bb,*edge.EdgeVertex(1));}
     else { // inside
      Real8   aa[3];
      Real8 s = deta[0]+deta[1]+deta[2];
      aa[0]=deta[0]/s;
      aa[1]=deta[1]/s;
      aa[2]=deta[2]/s;
      return Metric(aa,(*t)[0],(*t)[1],(*t)[2]);
     }
  }


void ListofIntersectionTriangles::SplitEdge(const Triangles & Bh,
       const R2 &A,const R2  &B,int nbegin)
{ //  SplitEdge
  //  if(SHOW)  cout << " splitedge " << A << B << " " <<  nbegin << endl;
  Triangle *tbegin, *t;

  Icoor2 deta[3], deti,detj;
  Real8 ba[3];
  int nbt =0,ifirst=-1,ilast;
  int i0,i1,i2;
  int ocut,i,j,k=-1;
  //  int OnAVertices =0;
  Icoor2 dt[3];
  I2 a = Bh.toI2(A) ,b= Bh.toI2(B);// compute  the Icoor a,b
  I2 vi,vj;  
  int iedge =-1;// not a edge

  if(nbegin)  {// optimisation 
    // we suppose  knowing the starting  triangle
    t=tbegin=lIntTria[ilast=(Size-1)].t;
    if (tbegin->det>=0) 
    ifirst = ilast;}  
  else {// not optimisation 
    init();
    t=tbegin = Bh.FindTriangleContening(a,deta);
    //    if(SHOW) cout <<t << " " << Real8(deta[0])/t->det<< " " << Real8(deta[1])/t->det
    //		  << " " << Real8(deta[2])/t->det << endl;
    if( t->det>=0)
      ilast=NewItem(t,Real8(deta[0])/t->det,Real8(deta[1])/t->det,Real8(deta[2])/t->det);
    else 
     {// find the nearest boundary edge  of the vertex A
      // find a edge or such normal projection a the edge IJ is on the edge
      //   <=> IJ.IA >=0 && IJ.AJ >=0
      ilast=ifirst;
      double ba,bb;
      TriangleAdjacent edge=CloseBoundaryEdge(a,t,ba,bb);
      Vertex & v0 = *edge.EdgeVertex(0), & v1 = *edge.EdgeVertex(1);
      NewItem(A,Metric(ba,v0,bb,v1));
      t=edge;
      // test if the point b is in the same side
      if (det(v0.i,v1.i,b)>=0) {
	//cout << " All the edge " << A << B << endl;
	TriangleAdjacent edge=CloseBoundaryEdge(a,t,ba,bb);
	Vertex & v0 = *edge.EdgeVertex(0), & v1 = *edge.EdgeVertex(1);
	NewItem(A,Metric(ba,v0,bb,v1));
	return;
      }
     } // find the nearest boundary edge  of the vertex A
   } // end not optimisation 
  if (t->det<0) {  // outside departure
    while (t->det <0) { // intersection boundary edge and a,b,
      k=(*t)(0) ?  ((  (*t)(1) ? ( (*t)(2) ? -1 : 2) : 1  )) : 0;
      assert(k>=0);
      ocut = OppositeEdge[k];
      i=VerticesOfTriangularEdge[ocut][0];
      j=VerticesOfTriangularEdge[ocut][1];
      vi=(*t)[i];
      vj=(*t)[j];
      deti = bamg::det(a,b,vi);
      detj = bamg::det(a,b,vj);
     //  if(SHOW) {  penthickness(3);
// 	Move(vi);Line(vj);CurrentTh->inquire();penthickness(1);
//         cout << Bh.Number(tbegin) << " " << Bh.Number(t) << " i= " << i <<" j= " <<  j << " k=" << k 
//       	   << " deti= " << deti << " detj= " << detj 
// 	     << " v = " << Bh.Number((*t)[i]) << (*t)[i].r <<  " " << Bh.Number((*t)[j]) << (*t)[j].r  << endl;}
      if (deti>0) // go to  i direction on gamma
	ocut = PreviousEdge[ocut];      
      else if (detj<=0) // go to j direction on gamma
	ocut = NextEdge[ocut];         
      TriangleAdjacent tadj =t->Adj(ocut);
      t = tadj;
      iedge= tadj; 
      if (t == tbegin) { // 
	double ba,bb;
	if (verbosity>7) 
	  cout << "       SplitEdge: All the edge " << A << B << nbegin <<  det(vi,vj,b) 
	       << " deti= " << deti <<  " detj=" <<detj << endl;
	TriangleAdjacent edge=CloseBoundaryEdge(a,t,ba,bb);
	Vertex & v0 = *edge.EdgeVertex(0), & v1 = *edge.EdgeVertex(1);
	NewItem(A,Metric(ba,v0,bb,v1));
	return;
	/*
	cerr << nbegin <<  det(vi,vj,b) << " deti= " << deti <<  " detj=" <<detj << endl;
	cerr << "SplitEdge on boucle A" << A << " B = " << B << endl;

#ifdef DRAWING
	reffecran();
	Bh.Draw();
        penthickness(5);
	Move(A);
	Line(B);
        penthickness(1);

	Bh.inquire();
        penthickness(5);
	Move(A);
	Line(B);
        penthickness(1);
	Bh.inquire();
#endif	
	MeshError(997);*/
      }
    } //  end while (t->det <0)
    // theoriticaly we have: deti =<0 and detj>0
  
      // computation of barycentric coor
    // test if the point b is on size on t
    // we revert vi,vj because vi,vj is def in Adj triangle
    if ( det(vi,vj,b)>=0) {
      if (verbosity>7)
      cout << "       SplitEdge: all AB outside " << A << B << endl;
      t=tbegin;
      Real8 ba,bb;
      TriangleAdjacent edge=CloseBoundaryEdge(b,t,ba,bb);
      NewItem(B,Metric(ba,*edge.EdgeVertex(0),bb,*edge.EdgeVertex(1)));
      return;
    }
    else
      {
	k = OppositeVertex[iedge];
	i=VerticesOfTriangularEdge[iedge][0];
	j=VerticesOfTriangularEdge[iedge][1];
	Real8 dij = detj-deti;
	assert(i+j+k == 0 + 1 +2);
	ba[j] =  detj/dij;
	ba[i] = -deti/dij;
	ba[k] = 0;
// 	if(SHOW) cout << i << " " << j << " " << k << " " << ba[i] << " " << ba[j] << endl;
	ilast=NewItem(t,ba[0],ba[1],ba[2]); }
  }  //  outside departure

     
   
  // recherche the intersection of [a,b] with Bh Mesh.
  // we know  a triangle ta contening the vertex a
  // we have 2 case for intersection [a,b] with a edge [A,B] of Bh
  // 1) the intersection point is in ]A,B[
  // 2)                        is A or B
  // first version --- 
  for (;;) {
    //    t->Draw();
    if (iedge < 0) {
      i0 =0;i1=1;i2=2;
      dt[0] =bamg::det(a,b,(*t)[0]);
      dt[1] =bamg::det(a,b,(*t)[1]);
      dt[2] =bamg::det(a,b,(*t)[2]);}
    else {
      i2 = iedge;
      i0 = NextEdge[i2];
      i1 = NextEdge[i0]; 
      dt[VerticesOfTriangularEdge[iedge][0]] = detj;// we revert i,j because
      dt[VerticesOfTriangularEdge[iedge][1]] = deti;// we take the Triangle by the other side
      dt[iedge] = det(a,b,(*t)[OppositeVertex[iedge]]);}
    
    // so we have just to see the transition from - to + of the det0..2 on edge of t
    // because we are going from a to b
    if       ((dt[i=VerticesOfTriangularEdge[i0][0]] <  0) &&
              ( dt[j=VerticesOfTriangularEdge[i0][1]] > 0))
      ocut =i0;
    else  if ((dt[i=VerticesOfTriangularEdge[i1][0]] <  0) &&
              (dt[j=VerticesOfTriangularEdge[i1][1]] >  0))
      ocut =i1;
    else  if ((dt[i=VerticesOfTriangularEdge[i2][0]] <  0) && 
              (dt[j=VerticesOfTriangularEdge[i2][1]] >  0))
      ocut =i2;
    else if   ((dt[i=VerticesOfTriangularEdge[i0][0]] == 0) &&
              ( dt[j=VerticesOfTriangularEdge[i0][1]] >  0))
      ocut =i0;
    else  if ((dt[i=VerticesOfTriangularEdge[i1][0]] == 0) &&
              (dt[j=VerticesOfTriangularEdge[i1][1]] >  0))
      ocut =i1;
    else  if ((dt[i=VerticesOfTriangularEdge[i2][0]] == 0) && 
              (dt[j=VerticesOfTriangularEdge[i2][1]] >  0))
      ocut =i2;
    else if   ((dt[i=VerticesOfTriangularEdge[i0][0]] <  0) &&
              ( dt[j=VerticesOfTriangularEdge[i0][1]] == 0))
      ocut =i0;
    else  if ((dt[i=VerticesOfTriangularEdge[i1][0]] <  0) &&
              (dt[j=VerticesOfTriangularEdge[i1][1]] == 0))
      ocut =i1;
    else  if ((dt[i=VerticesOfTriangularEdge[i2][0]] <  0) && 
              (dt[j=VerticesOfTriangularEdge[i2][1]] == 0))
      ocut =i2;
    else { //  On a edge (2 zero)
      k =0;
      if (dt[0]) ocut=0,k++; 
      if (dt[1]) ocut=1,k++; 
      if (dt[2]) ocut=2,k++;
      if(k == 1) {
        if (dt[ocut] >0) // triangle upper AB
          ocut = NextEdge[ocut];
        i= VerticesOfTriangularEdge[ocut][0];
        j= VerticesOfTriangularEdge[ocut][1];
      }
      else  {
        cerr << " Bug Split Edge " << endl;
        cerr << " dt[0]= " << dt[0] 
             << " dt[1]= " << dt[1] 
             << " dt[2]= "<< dt[2] << endl;
        cerr << i0 << " " << i1 << " " << i2 << endl;
        cerr << " A = " << A << " B= " << B << endl;
        cerr << " Triangle t = " <<  *t << endl;
        cerr << (*t)[0] << (*t)[1] << (*t)[0] << endl;
        cerr << " nbt = " << nbt << endl;
        MeshError(100);}}
    
    k = OppositeVertex[ocut];

    Icoor2 detbij = bamg::det((*t)[i],(*t)[j],b);

    
    if (detbij >= 0) { //we find the triangle contening b
      dt[0]=bamg::det((*t)[1],(*t)[2],b);
      dt[1]=bamg::det((*t)[2],(*t)[0],b);
      dt[2]=bamg::det((*t)[0],(*t)[1],b);
#ifdef DEBUG 
      assert(dt[0] >= 0);
      assert(dt[1] >= 0);
      assert(dt[2] >= 0);
#endif
      Real8 dd = t->det;
      NewItem(t,dt[0]/dd,dt[1]/dd,dt[2]/dd);      
      return ;}
    else { // next triangle by  adjacent by edge ocut 
      deti = dt[i];
      detj = dt[j];
      Real4 dij = detj-deti;
      ba[i] =  detj/dij;
      ba[j] = -deti/dij;
      ba[3-i-j ] = 0;
      ilast=NewItem(t, ba[0],ba[1],ba[2]);      
      
      TriangleAdjacent ta =t->Adj(ocut);
      t = ta;
      iedge= ta; 
      if (t->det <= 0)  {
        double ba,bb;
        TriangleAdjacent edge=CloseBoundaryEdge(b,t,ba,bb);
        NewItem(B,Metric(ba,*edge.EdgeVertex(0),bb,*edge.EdgeVertex(1)));
	// 	cout << " return " << ba << " " << bb << endl;
	// ajoute le 03 frev 1997 par F. hecht
        return;
        }
     }// we  go outside of omega 
  } // for(;;)
 
   
} // routine SplitEdge


int  ListofIntersectionTriangles::NewItem(Triangle * tt,Real8 d0,Real8 d1,Real8 d2) { 
  register int n;
  R2 x(0,0);
  if ( d0) x =      (*tt)[0].r * d0;
  if ( d1) x = x +  (*tt)[1].r * d1;
  if ( d2) x = x +  (*tt)[2].r * d2;
  // newer add same point 
  if(!Size ||  Norme2_2(lIntTria[Size-1].x-x)) {
    if (Size==MaxSize) ReShape();
    lIntTria[Size].t=tt;
    lIntTria[Size].bary[0]=d0;
    lIntTria[Size].bary[1]=d1;
    lIntTria[Size].bary[2]=d2;
    lIntTria[Size].x = x;
    Metric m0,m1,m2;
    register Vertex * v;
    if ((v=(*tt)(0))) m0    = v->m;
    if ((v=(*tt)(1))) m1    = v->m;
    if ((v=(*tt)(2))) m2    = v->m;
    lIntTria[Size].m =  Metric(lIntTria[Size].bary,m0,m1,m2);
#ifdef DEBUG1
    if(SHOW) { cout << "SHOW       ++ NewItem =" << Size << x ;
    cout << " " << d0 << " " << d1 << " " << d2 <<endl;}
#endif
    n=Size++;}
  else n=Size-1;
  return n;
}
int ListofIntersectionTriangles::NewItem(R2 A,const Metric & mm) {  
  register int n;
  if(!Size ||  Norme2_2(lIntTria[Size-1].x-A)) {
    if (Size==MaxSize) ReShape();
    lIntTria[Size].t=0;
    lIntTria[Size].x=A;
    lIntTria[Size].m=mm;
#ifdef DEBUG1
    if (SHOW)  cout << "SHOW       ++ NewItem A" << Size << A << endl;
#endif
    n=Size++;
   }
  else  n=Size-1;
 return  n; 
}

Real8  ListofIntersectionTriangles::Length()
{
  //  cout << " n= " << Size << ":" ;
  assert(Size>0);
  // computation of the length      
  R2 C;
  Metric Mx,My;
  int ii,jj;
  R2 x,y,xy;
  
  SegInterpolation *SegI=lSegsI;
  SegI=lSegsI;
  lSegsI[NbSeg].last=Size;// improvement 
  
  int EndSeg=Size;
     
  y = lIntTria[0].x;
  Real8 sxy, s = 0;
  lIntTria[0].s =0;
  SegI->lBegin=s;

  for (jj=0,ii=1;ii<Size;jj=ii++) 
    {  
      // seg jj,ii
      x=y;
      y = lIntTria[ii].x;
      xy = y-x;
      Mx = lIntTria[ii].m;
      My = lIntTria[jj].m;
      //      Real8 &sx=  lIntTria[ii].sp; // previous seg
      //  Real8 &sy=  lIntTria[jj].sn; // next seg
      //      sx = Mx(xy);
      //      sy = My(xy);
      //   sxy = (Mx(xy)+ My(xy))/2.0;
      sxy =  LengthInterpole(Mx,My,xy);
      s += sxy;
      lIntTria[ii].s = s;
      if (ii == EndSeg) 
	SegI->lEnd=s,
	  SegI++,
	  EndSeg=SegI->last,
	  SegI->lBegin=s;
      
      //    cout << ii << " " << jj << x<< y <<xy << s <<  lIntTria[ii].m  ;
    }
  len = s;
  SegI->lEnd=s;

  //  cout << " len= " << s << endl;
  return s;
}

Int4 ListofIntersectionTriangles::NewPoints(Vertex * vertices,Int4 & nbv,Int4  nbvx)
{

  const Int4 nbvold = nbv;
  Real8 s = Length();
  if (s <  1.5 ) return 0;
  //////////////////////   
  int ii = 1 ;
  R2 y,x;
  Metric My,Mx ;
  Real8 sx =0,sy;
  int nbi = Max(2,(int) (s+0.5));
  Real8 sint = s/nbi;
  Real8 si = sint;

  int EndSeg=Size;
  SegInterpolation *SegI=0;
  if (NbSeg) 
    SegI=lSegsI,EndSeg=SegI->last;
  
  for (int k=1;k<nbi;k++)
    {
      while ((ii < Size) && ( lIntTria[ii].s <= si )) 
	if (ii++ == EndSeg) 
	  SegI++,EndSeg=SegI->last;

      int ii1=ii-1;
      x  =lIntTria[ii1].x;
      sx =lIntTria[ii1].s;
      Metric Mx=lIntTria[ii1].m;
#ifdef DEBUG    
      double lx = lIntTria[ii-1].sn;
#endif
      y  =lIntTria[ii].x;
      sy =lIntTria[ii].s;
      Metric My=lIntTria[ii].m;
#ifdef DEBUG    
      double ly =lIntTria[ii].sp;  
      assert( sx <= si);
      assert( si <= sy);
      assert( sy != sx);
#endif 

      Real8 lxy = sy-sx;
      Real8 cy = abscisseInterpole(Mx,My,y-x,(si-sx)/lxy);
      
      R2 C;
      Real8 cx = 1-cy;
      C = SegI ? SegI->F(si): x * cx + y *cy;

    si += sint;
    if ( nbv<nbvx) {
      vertices[nbv].r = C;
      vertices[nbv++].m = Metric(cx,lIntTria[ii-1].m,cy,lIntTria[ii].m);
      if((verbosity/100%10)==2)
      cout << "   -- Add point " << nbv-1 << " " << vertices[nbv-1] << " " << vertices[nbv-1].m << endl;

#ifdef DEBUG
      if(k>1) {
	R2 AB = vertices[nbv-2].r - vertices[nbv-1].r ;
	Real8 dp = LengthInterpole(vertices[nbv-2].m,vertices[nbv-1].m,AB);
	if (dp > 1.6) {
	  cerr << "PB calcul new Int.  points trop loin l=" << dp << " v=" << nbv-1 << " " << nbv-2 <<Mx<<My<<y-x << endl;
	}
	}
#endif
    }
    else return nbv-nbvold;
  }
  return nbv-nbvold;
}

int SwapForForcingEdge(Vertex   *  & pva ,Vertex  * &   pvb ,
      TriangleAdjacent & tt1,Icoor2 & dets1, Icoor2 & detsa,Icoor2 & detsb, int & NbSwap)
{ // l'arete ta coupe l'arete pva pvb
  // de cas apres le swap sa coupe toujours
  // on cherche l'arete suivante 
  // on suppose que detsa >0 et detsb <0
  // attention la routine echange pva et pvb 

   if(tt1.Locked()) return 0; // frontiere croise 

  TriangleAdjacent tt2 = Adj(tt1);
  Triangle *t1=tt1,*t2=tt2;// les 2 triangles adjacent
  Int1 a1=tt1,a2=tt2;// les 2 numero de l arete dans les 2 triangles
  assert ( a1 >= 0 && a1 < 3 );
   
  Vertex & sa= (* t1)[VerticesOfTriangularEdge[a1][0]];
  Vertex & s1= (*t1)[OppositeVertex[a1]];
  Vertex & s2= (*t2)[OppositeVertex[a2]];
  

  Icoor2 dets2 = det(*pva,*pvb,s2);

#ifdef DEBUG
  Vertex & sb= (*t1)[VerticesOfTriangularEdge[a1][1]];
  Icoor2 wdets1 = det(*pva,*pvb,s1);  
  Icoor2 wdetsa = det(*pva,*pvb,sa);
  Icoor2 wdetsb = det(*pva,*pvb,sb);
  assert(wdets1 == dets1);
  assert(wdetsa == detsa);
  assert(wdetsb == detsb);
#endif
  
  Icoor2 det1=t1->det , det2=t2->det ;
#ifdef DEBUG  
  assert(det1>0 && det2 >0);
  Icoor2 ddet1 = det((*t1)[0],(*t1)[1],(*t1)[2]);
  Icoor2 ddet2 = det((*t2)[0],(*t2)[1],(*t2)[2]);
   if ((det1 != ddet1) || (det2 != ddet2) )
    {
      assert(det1 == ddet1);
      assert(det2 == ddet2);
     }
   Icoor2 detvasasb = det(*pva,sa,sb);
   Icoor2 detvbsasb = det(*pvb,sa,sb);
   if (  CurrentTh && !  ( ( (detvasasb <= 0) && (detvbsasb >= 0)) || ( (detvasasb >= 0) && (detvbsasb <= 0)))) 
     {
       cout << " detvasasb =" <<  detvasasb << "detvbsasb = " <<  detvbsasb 
	    << " " << pva << " " << pvb << " "  <<CurrentTh <<endl;
#ifdef DRAWING1
       reffecran();
       CurrentTh->Draw();
       penthickness(10);
       pva->MoveTo();pvb->LineTo();
       penthickness(1);
       CurrentTh->inquire();
#endif       
   }
     assert( ( (detvasasb <= 0) && (detvbsasb >= 0)) || ( (detvasasb >= 0) && (detvbsasb <= 0)));
#endif

  Icoor2 detT = det1+det2;
  assert((det1>0 ) && (det2 > 0));
  assert ( (detsa < 0) && (detsb >0) ); // [a,b] cut infinite line va,bb
  Icoor2 ndet1 = bamg::det(s1,sa,s2);
  Icoor2 ndet2 = detT - ndet1;

  int ToSwap =0; //pas de swap
  if ((ndet1 >0) && (ndet2 >0)) 
    { // on peut swaper  
      if ((dets1 <=0 && dets2 <=0) || (dets2 >=0 && detsb >=0))
        ToSwap =1; 
      else // swap alleatoire 
        if (BinaryRand()) 
          ToSwap =2; 
    }
#ifdef DEBUG
  if (ForDebugging) {
    cerr << "swap = " << ToSwap << " ndet1 " << ndet1 << ", ndet2 " << ndet2 << "det1  " << det1 << " det2 " <<  det2  
	 << " if1 = " << ((ndet1 >0) && (ndet2 >0)) 
	 << " if2 = " << ((dets1 <=0 && dets2 <=0) || (dets2 >=0 && detsb >=0)) << endl;
#ifdef DRAWING
  couleur(0);
  t1->Draw();
  t2->Draw();
#endif
  }
#endif
  if (ToSwap) NbSwap++,
     bamg::swap(t1,a1,t2,a2,&s1,&s2,ndet1,ndet2);
  
#ifdef DEBUG
  if (ForDebugging) {
#ifdef DRAWING
  couleur(4);
  t1->Draw();
  t2->Draw();
  rattente(1);
#endif
  }
#endif
  int ret=1;

  if (dets2 < 0) {// haut
    dets1 = ToSwap ? dets1 : detsa ;
    detsa = dets2; 
    tt1 =  Previous(tt2) ;}
  else if (dets2 > 0){// bas 
    dets1 = ToSwap ? dets1 : detsb ;
    detsb = dets2;
    //xxxx tt1 = ToSwap ? tt1 : Next(tt2);
    if(!ToSwap) tt1 =  Next(tt2);
    }
  else { // changement de sens 
     if (ForDebugging)  cout << "changement de sens" << endl;
    ret = -1;
    Exchange(pva,pvb);
    Exchange(detsa,detsb);
    Exchange(dets1,dets2);
    Exchange(tt1,tt2);
    dets1=-dets1;
    dets2=-dets2;
    detsa=-detsa;
    detsb=-detsb;

    if (ToSwap) 
      if (dets2 < 0) {// haut
        dets1 = (ToSwap ? dets1 : detsa) ;
        detsa = dets2; 
        tt1 =  Previous(tt2) ;}
      else if (dets2 > 0){// bas 
        dets1 = (ToSwap ? dets1 : detsb) ;
        detsb =  dets2;
        if(!ToSwap) tt1 =  Next(tt2);
       }
      else {// on a fin ???
        tt1 = Next(tt2);
        ret =0;}

  }
  return ret;
}

int ForceEdge(Vertex &a, Vertex & b,TriangleAdjacent & taret)  
{ 
#ifdef DEBUG
 restart: // for debug 
#endif
  int NbSwap =0;
  assert(a.t && b.t); // the 2 vertex is in a mesh 
  int k=0;
  taret=TriangleAdjacent(0,0); // erreur 
  
  TriangleAdjacent tta(a.t,EdgesVertexTriangle[a.vint][0]);
  Vertex   *v1, *v2 = tta.EdgeVertex(0),*vbegin =v2;
  // we turn around a in the  direct sens  
   
  Icoor2 det2 = v2 ? det(*v2,a,b): -1 , det1;
  if(v2) // normal case 
    det2 = det(*v2,a,b);
  else { // no chance infini vertex try the next
    tta= Previous(Adj(tta));
    v2 = tta.EdgeVertex(0);
    vbegin =v2;
    assert(v2);
    det2 = det(*v2,a,b);
 //   cout << " No Change try the next" << endl;
  }

#ifdef DRAWING1
  a.MoveTo();b.LineTo();
#endif

  while (v2 != &b) {
    TriangleAdjacent tc = Previous(Adj(tta));    
    v1 = v2; 
    v2 = tc.EdgeVertex(0);
    det1 = det2;
#ifdef DEBUG
    assert( v1 ==  tta.EdgeVertex(0));
    assert( &a ==  tc.EdgeVertex(1) );
#endif
    det2 =  v2 ? det(*v2,a,b): det2; 
    
    if((det1 < 0) && (det2 >0)) { 
      // try to force the edge 
      Vertex * va = &a, *vb = &b;
      tc = Previous(tc);
      assert ( v1 && v2);
      Icoor2 detss = 0,l=0,ks;
      // cout << "Real ForcingEdge " << *va << *vb << detss << endl;
#ifdef DEBUG
      Icoor2 dettt1 =  det(*v1,a,b);
      Icoor2 dettt2 =  det(*v2,a,b);

      if (!(dettt1==det1 && dettt2==det2))
	{
	  assert(ForDebugging==0);
          ForDebugging=1;
	  goto restart;
	}
	
#endif 
      while ((ks=SwapForForcingEdge(  va,  vb, tc, detss, det1,det2,NbSwap)))
	if(l++ > 10000000) {
	  cerr << " Loop in forcing Egde AB" 
               <<"\n vertex A " << a 
               <<"\n vertex B " <<  b 
               <<"\n nb de swap " << NbSwap 
               <<"\n nb of try  swap too big = " <<  l << " gearter than " <<  1000000 << endl;
   
         if ( CurrentTh ) 
            cerr << " vertex number " << CurrentTh->Number(a) << " " <<  CurrentTh->Number(b) << endl;
#ifdef DEBUG
	  ForDebugging = 1;
#endif
#ifdef DRAWING1
          if (  CurrentTh ) {
	    reffecran();
            couleur(6);
	    CurrentTh->Draw();
            couleur(1);
	    penthickness(10);
	    a.MoveTo();b.LineTo();
	    penthickness(1);
	    CurrentTh->inquire();
            couleur(6);
	    l=0;
	    reffecran();
            while (ks=SwapForForcingEdge(  va,  vb, tc, detss, det1,det2,NbSwap) && (l++ < 1000))
             cerr << " " << CurrentTh->Number(tc.EdgeVertex(0))<<" " <<CurrentTh->Number(tc.EdgeVertex(1)) << " ";
	  }
#endif	  
	  MeshError(990);
	}
      Vertex *aa = tc.EdgeVertex(0), *bb = tc.EdgeVertex(1);
      if (( aa == &a ) && (bb == &b) ||  (bb ==  &a ) && (aa == &b)) {
	tc.SetLock();
	a.Optim(1,0);
	b.Optim(1,0);
	taret = tc;
	return NbSwap;
      }
      else 
      {
	  taret = tc;
	  return -2; // error  boundary is crossing
/*	  cerr << "Fatal Error  boundary is crossing  ";
	  if(CurrentTh)
	  {
	      cerr << " edge:  [" << CurrentTh->Number(a) << ", " << CurrentTh->Number(b) <<  " ] and [ ";
	      cerr    << CurrentTh->Number(aa) << " " << CurrentTh->Number(bb) << " ] " << endl;
	  }
	  MeshError(991);
*/	  
      }
    }
    tta = tc;
    assert(k++<2000);
    if ( vbegin == v2 ) return -1;// error 
  }

  tta.SetLock();
  taret=tta;
  a.Optim(1,0);
  b.Optim(1,0);
  return NbSwap; 
}


int Triangle::swap(Int2 a,int koption){
#ifdef DEBUG
  if(a &4 ) return 0;// arete lock 
  int munswap1 = a/4;
  a &=3;
#else
 if(a/4 !=0) return 0;// arete lock or MarkUnSwap
#endif

  register Triangle *t1=this,*t2=at[a];// les 2 triangles adjacent
  register Int1 a1=a,a2=aa[a];// les 2 numero de l arete dans les 2 triangles
#ifdef DEBUG
  if(a2 & 4) return 0; // arete lock
  int munswap2 = a2/4;
  a2 &= 3;
#else
  if(a2/4 !=0) return 0; // arete lock or MarkUnSwap
#endif
  
  register Vertex  *sa=t1->ns[VerticesOfTriangularEdge[a1][0]];
  register Vertex  *sb=t1->ns[VerticesOfTriangularEdge[a1][1]];
  register Vertex  *s1=t1->ns[OppositeVertex[a1]];
  register Vertex  *s2=t2->ns[OppositeVertex[a2]];

#ifdef DEBUG
  assert ( a >= 0 && a < 3 );  
#endif
  
   Icoor2 det1=t1->det , det2=t2->det ;
   Icoor2 detT = det1+det2;
   Icoor2 detA = Abs(det1) + Abs(det2);
   Icoor2 detMin = Min(det1,det2);

   int OnSwap = 0;       
   // si 2 triangle infini (bord) => detT = -2;
   if (sa == 0) {// les deux triangles sont frontieres
     det2=bamg::det(s2->i,sb->i,s1->i);
     OnSwap = det2 >0;}
   else if (sb == 0) { // les deux triangles sont frontieres
     det1=bamg::det(s1->i,sa->i,s2->i);
     OnSwap = det1 >0;}
   else if(( s1 != 0) && (s2 != 0) ) {
     det1 = bamg::det(s1->i,sa->i,s2->i);
     det2 = detT - det1;
     OnSwap = (Abs(det1) + Abs(det2)) < detA;
   
     Icoor2 detMinNew=Min(det1,det2);
     //     if (detMin<0 && (Abs(det1) + Abs(det2) == detA)) OnSwap=BinaryRand();// just for test   
     if (! OnSwap &&(detMinNew>0)) {
       OnSwap = detMin ==0;
       if (! OnSwap) {
	 int  kopt = koption;
        while (1)
	 if(kopt) {
	 // critere de Delaunay pure isotrope
	 register Icoor2 xb1 = sb->i.x - s1->i.x,
	   x21 = s2->i.x - s1->i.x,
	   yb1 = sb->i.y - s1->i.y,
	   y21 = s2->i.y - s1->i.y,
	   xba = sb->i.x - sa->i.x, 
	   x2a = s2->i.x - sa->i.x,
	   yba = sb->i.y - sa->i.y,
	   y2a = s2->i.y - sa->i.y;
	 register double
	   cosb12 =  double(xb1*x21 + yb1*y21),
	   cosba2 =  double(xba*x2a + yba*y2a) ,
	   sinb12 = double(det2),
	   sinba2 = double(t2->det);

	 
	 // angle b12 > angle ba2 => cotg(angle b12) < cotg(angle ba2)
	 OnSwap =  ((double) cosb12 * (double)  sinba2) <  ((double) cosba2 * (double) sinb12);
//  	 if(CurrentTh) 
//  	   cout << "swap " << CurrentTh->Number(sa) << " " << CurrentTh->Number(sb) << " " ;
//  	 cout <<  cosb12 << " " <<  sinba2 << " "  <<  cosba2 << " " << sinb12 
//  	      << " Onswap = " <<  OnSwap << endl;
	  break;
	 }
	 else 
	   {	
	     // critere de Delaunay anisotrope 
	     Real8 som;
	     I2 AB=(I2) *sb - (I2) *sa;
	     I2 MAB2=((I2) *sb + (I2) *sa);
	     R2 MAB(MAB2.x*0.5,MAB2.y*0.5);
	     I2 A1=(I2) *s1 - (I2) *sa;
	     I2 D = (I2) * s1 - (I2) * sb ;
	     R2 S2(s2->i.x,s2->i.y);
	     R2 S1(s1->i.x,s1->i.y);
	     {
	       Metric M=s1->m;
	       R2 ABo = M.Orthogonal(AB);
	       R2 A1o = M.Orthogonal(A1);
	       // (A+B)+ x ABo = (S1+B)/2+ y A1 
	       // ABo x - A1o y =  (S1+B)/2-(A+B)/2 = (S1-B)/2 = D/2
	       double dd = Abs(ABo.x*A1o.y)+Abs(ABo.y*A1o.x);
	       double d = (ABo.x*A1o.y - ABo.y*A1o.x)*2; // because D/2
	       if (Abs(d) > dd*1.e-3) {
		 R2 C(MAB+ABo*((D.x*A1o.y - D.y*A1o.x)/d));
		 som  = M(C - S2)/M(C - S1);
	       } else 
		{kopt=1;continue;}
		
	     }
	     {
	       Metric M=s2->m;
	       R2 ABo = M.Orthogonal(AB);
	       R2 A1o = M.Orthogonal(A1);
	       // (A+B)+ x ABo = (S1+B)/2+ y A1 
	       // ABo x - A1o y =  (S1+B)/2-(A+B)/2 = (S1-B)/2 = D/2 
	       double dd = Abs(ABo.x*A1o.y)+Abs(ABo.y*A1o.x);
	       double d = (ABo.x*A1o.y - ABo.y*A1o.x)*2; // because D/2
	       if(Abs(d) > dd*1.e-3) {
		 R2 C(MAB+ABo*((D.x*A1o.y - D.y*A1o.x)/d));
		 som  += M(C - S2)/M(C -  S1);
	       } else 
		{kopt=1;continue;}
	     }
	     OnSwap = som < 2;
	     break;
	 }
  
       } // OnSwap 
     } // (! OnSwap &&(det1 > 0) && (det2 > 0) )
   }
#ifdef  DEBUG1   
   if (OnSwap &&  ( munswap1  || munswap2)) {
     cout << " erreur Mark unswap T " << CurrentTh->Number(t1) << " " <<  CurrentTh->Number(t2) << endl
	  << *t1 << endl
          << *t2 << endl;
     return 0;
   }
#endif 
   if( OnSwap ) 
       bamg::swap(t1,a1,t2,a2,s1,s2,det1,det2);
   else {
     NbUnSwap ++;
      t1->SetMarkUnSwap(a1);     
   }
   return OnSwap;
}

Real8  Vertex::Smoothing(Triangles & Th,const Triangles & BTh,Triangle  * & tstart ,Real8 omega)
{
#ifdef DEBUG  
  register  Int4 NbSwap =0;
#endif
  register Vertex * s  = this;
  Vertex &vP = *s,vPsave=vP;
  //  if (vP.on) return 0;// Don't move boundary vertex  
  
  register Triangle * tbegin= t , *tria = t , *ttc;
 
  register int k=0,kk=0,j = EdgesVertexTriangle[vint][0],jc;
  R2 P(s->r),PNew(0,0);
  //  cout << BTh.quadtree << " " <<  BTh.quadtree->root << endl;
  // assert(BTh.quadtree && BTh.quadtree->root);
  do {
	k++; 
	
#ifdef DEBUG
    assert( s == & (*tria)[VerticesOfTriangularEdge[j][1]] );
    assert( tria->det >0);
#endif
    if (!tria->Hidden(j))
      {
	Vertex &vQ = (*tria)[VerticesOfTriangularEdge[j][0]]; 
	
	R2 Q = vQ,QP(P-Q);
	Real8 lQP = LengthInterpole(vP,vQ,QP);
	PNew += Q+QP/Max(lQP,1e-20);
	kk ++;
      }
     ttc =  tria->TriangleAdj(j);
     jc = NextEdge[tria->NuEdgeTriangleAdj(j)];
     tria = ttc;
     j = NextEdge[jc];
     assert(k<2000);
  } while ( tbegin != tria); 
  if (kk<4) return 0;
  PNew = PNew/(Real8)kk;
  R2 Xmove((PNew-P)*omega);
  PNew = P+Xmove;
  Real8 delta=Norme2_2(Xmove); 
  
  
  // 
  Icoor2 deta[3];
  I2 IBTh  = BTh.toI2(PNew);
  
  tstart=BTh.FindTriangleContening(IBTh,deta,tstart);  
  
  if (tstart->det <0) 
    { // outside
      double ba,bb;
      TriangleAdjacent edge= CloseBoundaryEdge(IBTh,tstart,ba,bb) ;
      tstart = edge;
      vP.m= Metric(ba,*edge.EdgeVertex(0),bb,*edge.EdgeVertex(1));
    }
  else 
    { // inside
      Real8   aa[3];
      Real8 s = deta[0]+deta[1]+deta[2];
      aa[0]=deta[0]/s;
      aa[1]=deta[1]/s;
      aa[2]=deta[2]/s;
      vP.m = Metric(aa,(*tstart)[0],(*tstart)[1],(*tstart)[2]);
    }
  
  // recompute the det of the triangle
  vP.r = PNew;
  
  vP.i = Th.toI2(PNew);

  Vertex vPnew = vP;
  
  int ok=1;
  int loop=1;
  k=0;
  while (ok) 
    {
      ok =0;
      do {
	k++; 
	double detold = tria->det;
	tria->det =  bamg::det( (*tria)[0],(*tria)[1]  ,(*tria)[2]);
	if (loop) 
	  {
	    Vertex *v0,*v1,*v2,*v3;
	    if (tria->det<0) ok =1;			       
	    else if (tria->Quadrangle(v0,v1,v2,v3))
	      {
		vP = vPsave;
		Real8 qold =QuadQuality(*v0,*v1,*v2,*v3);
		vP = vPnew;
		Real8 qnew = QuadQuality(*v0,*v1,*v2,*v3);
		if (qnew<qold) ok = 1;
	      }
	    else if ( (double)tria->det < detold/2 ) ok=1;
	    
	  }
        tria->SetUnMarkUnSwap(0);
        tria->SetUnMarkUnSwap(1);
        tria->SetUnMarkUnSwap(2);
	ttc =  tria->TriangleAdj(j);
	jc = NextEdge[tria->NuEdgeTriangleAdj(j)];
	tria = ttc;
	j = NextEdge[jc];
	assert(k<2000);
      } while ( tbegin != tria); 
      if (ok && loop) vP=vPsave; // no move 
      loop=0;
    }
  return delta;
}




void Triangles::Add( Vertex & s,Triangle * t, Icoor2 * det3) 
{
  // -------------------------------------------
  //             s2
  //                                            !
  //             /|\                            !
  //            / | \                           !
  //           /  |  \                          !
  //    tt1   /   |   \ tt0                     !
  //         /    |s   \                        !
  //        /     .     \                       !
  //       /  .      `   \                      !
  //      / .           ` \                     !
  //      ----------------                      !
  //   s0       tt2       s1
  //-------------------------------------------- 
  
  Triangle * tt[3]; // the 3 new Triangles
  Vertex &s0 = (* t)[0], &s1=(* t)[1], &s2=(* t)[2];
  Icoor2  det3local[3];
  int infv = &s0 ?  ((  &s1 ? ( &s2  ? -1 : 2) : 1  )) : 0;
  // infv = ordre of the infini vertex (null)
  register int nbd0 =0; // number of zero det3
  register int izerodet=-1,iedge; // izerodet = egde contening the vertex s
  Icoor2 detOld = t->det;
  
  if ( ( infv <0 ) && (detOld <0) ||  ( infv >=0  ) && (detOld >0) ) 
    {
      cerr << "  infv " << infv << " det = " << detOld << endl;
      cerr << Number(s) << " "<< Number(s0) << " "  
	   << Number(s1) << " "  << Number(s2) << endl;
      MeshError(3);
    }
  
  // if det3 do not exist then constuct det3
  if (!det3) { 
    det3 = det3local; // alloc 
    if ( infv<0 ) {
      det3[0]=bamg::det(s ,s1,s2);
      det3[1]=bamg::det(s0,s ,s2);
      det3[2]=bamg::det(s0,s1,s );}
    else { 
      // one of &s1  &s2  &s0 is NULL so (&si || &sj) <=> !&sk
      det3[0]=  &s0 ? -1  : bamg::det(s ,s1,s2) ;
      det3[1]=  &s1 ? -1 : bamg::det(s0,s ,s2) ;
      det3[2]=  &s2 ? -1 : bamg::det(s0,s1,s ) ;}}

  
  if (!det3[0]) izerodet=0,nbd0++;
  if (!det3[1]) izerodet=1,nbd0++;
  if (!det3[2]) izerodet=2,nbd0++;
  
  if  (nbd0 >0 ) // point s on a egde or on a vertex 
    if (nbd0 == 1) {
      iedge = OppositeEdge[izerodet];
      TriangleAdjacent ta = t->Adj(iedge);

#ifdef DEBUG1      
      cout << " the point " << Number(s) << " is the edge " <<  izerodet 
	   << " of " << Number(t)	   << " det3 = " 
	   << det3[0] << " " <<  det3[1] << " " <<  det3[2] << " " <<  endl;
      cout  << " ta = "  <<  ta << "ta->det =" << ((Triangle*) ta)->det  
	    << " "<< t->det<< endl;
#endif

      // the point is on the edge 
      // if the point is one the boundary 
      // add the point in outside part 
      if ( t->det >=0) { // inside triangle
	if ((( Triangle *) ta)->det < 0 ) {
	  // add in outside triangle 
	  Add(s,( Triangle *) ta);
	  return;}
      }}
    else {
      cerr << " bug  " << nbd0 <<endl;
      cerr << " Bug double points in " << endl ;
      cerr << " s = " << Number(s) << " " <<  s << endl;
      cerr << " s0 = "<< Number(s0) << " "  << s0 << endl;
      cerr << " s1 = "<< Number(s1) << " "  << s1 << endl;
      cerr << " s2 = "<< Number(s2) << " "  << s2 << endl;
      MeshError(5,this);}

  // remove de MarkUnSwap edge
  t->SetUnMarkUnSwap(0);     
  t->SetUnMarkUnSwap(1);     
  t->SetUnMarkUnSwap(2);

  tt[0]= t;
  tt[1]= &triangles[nbt++];
  tt[2]= &triangles[nbt++];
  
  if (nbt>nbtx) {
    cerr << " No enougth triangles " << endl;
    MeshError(999,this);
  }

  *tt[1]=   *tt[2]= *t;
// gestion of the link
   tt[0]->link=tt[1];
   tt[1]->link=tt[2]; 
  
  (* tt[0])(OppositeVertex[0])=&s;
  (* tt[1])(OppositeVertex[1])=&s;
  (* tt[2])(OppositeVertex[2])=&s;

  tt[0]->det=det3[0];
  tt[1]->det=det3[1];
  tt[2]->det=det3[2];         

  //  update adj des triangles externe 
  tt[0]->SetAdjAdj(0);
  tt[1]->SetAdjAdj(1);
  tt[2]->SetAdjAdj(2);
  //  update des adj des 3 triangle interne
  const int i0 = 0;
  const int i1= NextEdge[i0];
  const int i2 = PreviousEdge[i0];

  tt[i0]->SetAdj2(i2,tt[i2],i0);
  tt[i1]->SetAdj2(i0,tt[i0],i1);
  tt[i2]->SetAdj2(i1,tt[i1],i2);
  
  tt[0]->SetTriangleContainingTheVertex();
  tt[1]->SetTriangleContainingTheVertex();
  tt[2]->SetTriangleContainingTheVertex();
 
 
  // swap if the point s is on a edge
  if(izerodet>=0) {
    //  cout << " the point s is on a edge =>swap " << iedge << " "  << *tt[izerodet] << endl;
    int rswap =tt[izerodet]->swap(iedge);
    
    if (!rswap) 
     {
       cout << " Pb swap the point s is on a edge =>swap " << iedge << " "  << *tt[izerodet] << endl;
#ifdef DRAWING
       if(  CurrentTh &&  withrgraphique) 
        {
       reffecran();
       
       DrawMark(s.r);
          CurrentTh->inquire();
       DrawMark(s.r);
       rattente(1); 
        }      
#endif
     }
    assert(rswap);
  }
 
#ifdef DEBUG 
  tt[0]->check();
  tt[1]->check();
  tt[2]->check();
#endif
#ifdef DRAWING1 
  tt[0]->Draw();
  tt[1]->Draw();
  tt[2]->Draw();
#endif
  
}


Int4  Triangles::SplitInternalEdgeWithBorderVertices()
{
  Int4 NbSplitEdge=0;
  SetVertexFieldOn();  
  Int4 it;
  Int4 nbvold=nbv;
  for (it=0;it<nbt;it++)
    {
      Triangle &t=triangles[it];
      if (t.link)
	for (int j=0;j<3;j++)
	  if(!t.Locked(j) && !t.Hidden(j)){
	    Triangle &tt = *t.TriangleAdj(j);
	    if ( &tt && tt.link && it < Number(tt)) 
	      { // an internal edge 
		Vertex &v0 = t[VerticesOfTriangularEdge[j][0]];
		Vertex &v1 = t[VerticesOfTriangularEdge[j][1]];
		if (v0.on && v1.on)
		  {
		    R2 P= ((R2) v0 + (R2) v1)*0.5;
		    if ( nbv<nbvx) {
		      vertices[nbv].r = P;
		      vertices[nbv++].m = Metric(0.5,v0.m,0.5,v1.m);
		      vertices[nbv].ReferenceNumber=0;
		      vertices[nbv].DirOfSearch = NoDirOfSearch ;
		    }
		    NbSplitEdge++;
		    if (verbosity>7)
		      cout <<" Internal edge with two vertices on boundary" 
			   << Number(v0) << " " << Number(v1) << " by " <<  endl;
		  }
	      }
	  }
    }
  ReMakeTriangleContainingTheVertex();    
  if (nbvold!=nbv) 
    {
      Int4  iv = nbvold;
      Int4 NbSwap = 0;
      Icoor2 dete[3];  
      for (Int4 i=nbvold;i<nbv;i++) 
	{// for all the new point
	  Vertex & vi = vertices[i];
	  vi.i = toI2(vi.r);
      vi.r = toR2(vi.i);
      //      if (!quadtree->ToClose(vi,seuil,hi,hj)) {
      // a good new point 
      vi.ReferenceNumber=0; 
      vi.DirOfSearch =NoDirOfSearch;
      //	cout << " Add " << Number(vi) << " " << vi 
      // << "   " <<  Number(vi) << " <--> " << Number(vi) <<endl;
      Triangle *tcvi = FindTriangleContening(vi.i,dete);
      if (tcvi && !tcvi->link) {
	cout << i <<  " PB insert point " << Number(vi) << vi << Number(vi) 
	     << " tcvi = " << tcvi << " " << tcvi->link << endl;
	cout << (*tcvi)[1] <<  (*tcvi)[2] << endl;
	tcvi = FindTriangleContening(vi.i,dete);
	cout << (*tcvi)[1] <<  (*tcvi)[2] << endl;
#ifdef DRAWING1
	inquire();
	penthickness(5);
	DrawMark(vi.r);
	penthickness(1);
	inquire();
#endif
	
	MeshError(1001,this);
      }
      
      
      quadtree->Add(vi);
#ifdef DRAWING1
      DrawMark(vi.r);
#endif
      assert (tcvi && tcvi->det >= 0) ;// internal 
      Add(vi,tcvi,dete);
      NbSwap += vi.Optim(1);          
      iv++;
      //      }
	}
      if (verbosity>3) 
	{
	  cout << "    Nb Of New Point " << iv ;
	  cout << " Nb swap = " << NbSwap << " to  split internal edges with border vertices" ;}
      
      nbv = iv;
    }
  if (NbSplitEdge >  nbv-nbvold)
    cout << " Warning not enough vertices  to split all internal edges "  << endl
	 << "    we lost " << NbSplitEdge - ( nbv-nbvold) << " Edges Sorry " << endl;
  if (verbosity>2)
  cout << "SplitInternalEdgeWithBorderVertices: Number of splited edge " << NbSplitEdge << endl;
  return  NbSplitEdge;
}
Int4 Triangles::InsertNewPoints(Int4 nbvold,Int4 & NbTSwap)
 {
  Real8 seuil= 1.414/2 ;// for two close point 
  Int4 i;
 // insertion part --- 

  const Int4 nbvnew = nbv-nbvold;
  if (verbosity>5) 
    cout << "    Try to Insert the " <<nbvnew<< " new points " << endl;  
  Int4 NbSwap=0;
  Icoor2 dete[3];
  
  // construction d'un ordre aleatoire 
  if (! nbvnew) 
    return 0; 
  if (nbvnew) {
  const Int4 PrimeNumber= AGoodNumberPrimeWith(nbv)  ;
  Int4 k3 = rand()%nbvnew ; 
  for (Int4 is3=0; is3<nbvnew; is3++) {
    register Int4 j = nbvold +(k3 = (k3 + PrimeNumber)% nbvnew);
    register Int4 i = nbvold+is3; 
    ordre[i]= vertices + j;
    ordre[i]->ReferenceNumber=i;
  }
  // be carefull 
  Int4  iv = nbvold;
  for (i=nbvold;i<nbv;i++) 
    {// for all the new point
      Vertex & vi = *ordre[i];
      vi.i = toI2(vi.r);
      vi.r = toR2(vi.i);
      Real4 hx,hy;
      vi.m.Box(hx,hy);
      Icoor1 hi=(Icoor1) (hx*coefIcoor),hj=(Icoor1) (hy*coefIcoor);
      if (!quadtree->ToClose(vi,seuil,hi,hj)) 
        {
			// a good new point 
			Vertex & vj = vertices[iv];
			Int4 j = vj.ReferenceNumber; 
			assert( &vj== ordre[j]);
			if(i!=j)
			  { //  for valgring 
			    Exchange(vi,vj);
			    Exchange(ordre[j],ordre[i]);
			  }
		      vj.ReferenceNumber=0; 
			//	cout << " Add " << Number(vj) << " " << vj 
			// << "   " <<  Number(vi) << " <--> " << Number(vj) <<endl;
			Triangle *tcvj = FindTriangleContening(vj.i,dete);
			if (tcvj && !tcvj->link) 
			 {
			  cerr << i <<  " PB insert point " << Number(vj) << vj << Number(vi) 
			       << " tcvj = " << tcvj << " " << tcvj->link << endl;
			  cerr << (*tcvj)[1] <<  (*tcvj)[2] << endl;
			  tcvj = FindTriangleContening(vj.i,dete);
			  cout << (*tcvj)[1] <<  (*tcvj)[2] << endl;
#ifdef DRAWING1
		 	  inquire();
			  penthickness(5);
			  DrawMark(vj.r);
			  penthickness(1);
			  
			  inquire();
#endif
	  
	          MeshError(1001,this);
	         }
	
	
	        quadtree->Add(vj);
#ifdef DRAWING1
	        DrawMark(vj.r);
#endif
			assert (tcvj && tcvj->det >= 0) ;// internal 
			Add(vj,tcvj,dete);
			NbSwap += vj.Optim(1);          
			iv++;
         }
    } 
  if (verbosity>3) {
    cout << "    Nb Of New Point " << iv << " Nb Of To close Points " << nbv-iv ;
    cout << " Nb swap = " << NbSwap << " after " ;}

    nbv = iv;
  }

#ifdef DRAWING1
  inquire();
#endif

  for (i=nbvold;i<nbv;i++) 
    NbSwap += vertices[i].Optim(1);  
   if (verbosity>3) 
     cout << " NbSwap = " <<  NbSwap << endl;


  NbTSwap +=  NbSwap ;
#ifdef DEBUG
{  
  Int4 NbErr=0;
  Int4 i;
  for (i=0;i<nbt;i++)
    if (triangles[i].link) 
      {
      double dd =Det(triangles[i][1].r-triangles[i][0].r,triangles[i][2].r-triangles[i][0].r);
      if(dd <=0) 
	{
	  NbErr++;
	  cerr << " det triangle i " << i << " = " << triangles[i].det ;
	  cerr << " det triangle  " << dd ;
	    cerr << " Les trois sommets " ;
	    cerr << Number(triangles[i][0]) << " "  << Number(triangles[i][1]) << " " 
		 << Number(triangles[i][2]) << endl;
	    cerr << "I2     " <<triangles[i][0].r << triangles[i][1].r << triangles[i][2].r << endl;
	    cerr << "R2     " << triangles[i][0].i << triangles[i][1].i << triangles[i][2].i << endl;
	    cerr << "I2-R2 =" <<toR2(triangles[i][0].i)-triangles[i][0].r 
		 << toR2(triangles[i][1].i)-triangles[i][1].r
		 << toR2(triangles[i][2].i)-triangles[i][2].r << endl;
	}
      }
  if(NbErr) { 
#ifdef DRAWING
    Int4 kkk=0;
    //  UnMarkUnSwapTriangle();
    //  for (i=0;i<nbv;i++)
    //  kkk += vertices[i].Optim(0);
    if(verbosity>3)
    cout << "    Nb of swap louche " << kkk << endl;
    if(kkk) {
  for (i=0;i<nbt;i++)
    if (triangles[i].link) 
      {
      double dd =Det(triangles[i][1].r-triangles[i][0].r,triangles[i][2].r-triangles[i][0].r);
      if(dd <=0) 
	{
	  NbErr++;
	  cerr << " xxxdet triangle i " << i << " = " << triangles[i].det ;
	  cerr << " xxxdet triangle  " << dd ;
	    cerr << " xxxLes trois sommets " ;
	    cerr << Number(triangles[i][0]) << " "  << Number(triangles[i][1]) << " " 
		 << Number(triangles[i][2]) << endl;
	    cerr << triangles[i][0].r << triangles[i][1].r << triangles[i][2].r << endl;
	    cerr << triangles[i][0].i << triangles[i][1].i << triangles[i][2].i << endl;
	}
      } }
   inquire();
#endif
    //   MeshError(11);
  }
  
}
#endif
  return nbv-nbvold;
 
 }
void  Triangles::NewPoints(Triangles & Bh,int KeepBackVertex)
{ // Triangles::NewPoints
  Int4 nbtold(nbt),nbvold(nbv);
  if (verbosity>2) 
    cout << "  -- Triangles::NewPoints ";
  if (verbosity>3)cout <<  " nbv (in)  on Boundary  = " << nbv  <<endl;
  Int4 i,k;
  int j;
  Int4 *first_np_or_next_t = new Int4[nbtx];
  Int4 NbTSwap =0;
//    insert old point
    nbtold = nbt;
    
  if (KeepBackVertex && (&Bh != this) && (nbv+Bh.nbv< nbvx)) 
   {
  //   Bh.SetVertexFieldOn();
     for (i=0;i<Bh.nbv;i++)
      { 
        Vertex & bv = Bh[i];
        if (!bv.on) {
          vertices[nbv].r = bv.r;
          vertices[nbv++].m = bv.m;}
      }
     int nbv1=nbv;
     Bh.ReMakeTriangleContainingTheVertex();     
     InsertNewPoints(nbvold,NbTSwap)   ;            
     if (verbosity>2)
        cout <<  "      (Nb of Points from background mesh  = " 
             << nbv-nbvold  << " / " << nbv1-nbvold << ")" << endl;
   }  
  else 
    Bh.ReMakeTriangleContainingTheVertex();     

  Triangle *t;
  // generation of the list of next Triangle 
  // at 1 time we test all the triangles
  Int4 Headt =0,next_t;
  for(i=0;i<nbt;i++)
    first_np_or_next_t[i]=-(i+1);
  // end list i >= nbt 
  // the list of test triangle is 
  // the next traingle on i is  -first_np_or_next_t[i]
  int iter=0;
  // Big loop 
  do {
    iter++;
    nbtold = nbt;
    nbvold = nbv;
#ifdef DRAWING1  
  inquire();
#endif  

  // default size of  IntersectionTriangle

  i=Headt;
  next_t=-first_np_or_next_t[i];
  for(t=&triangles[i];i<nbt;t=&triangles[i=next_t],next_t=-first_np_or_next_t[i]) 
    { // for each triangle  t
      // we can change first_np_or_next_t[i]
      //      cout << " Do the triangle " << i << " Next_t=" << next_t << endl;
      assert(i>=0 && i < nbt );
      first_np_or_next_t[i] = iter; 
      for(j=0;j<3;j++)
	{ // for each edge 
	  TriangleAdjacent tj(t,j);
	  Vertex & vA = * tj.EdgeVertex(0);
	  Vertex & vB = * tj.EdgeVertex(1);
	  
	  if (!t->link) continue;// boundary
	  if (t->det <0) continue;
	  if (t->Locked(j)) continue;

	  TriangleAdjacent tadjj = t->Adj(j);	  
	  Triangle * ta= tadjj;

	  if (ta->det <0) continue;	  
	  
	  R2 A = vA;
	  R2 B = vB;
	  
	  k=Number(ta);
	  
	  if(first_np_or_next_t[k]==iter)  // this edge is done before 
	    continue; // next edge of the triangle 
	  
	  //const Int4 NbvOld = nbv;
	  lIntTria.SplitEdge(Bh,A,B);
	  lIntTria.NewPoints(vertices,nbv,nbvx);
	} // end loop for each edge 
      
    }// for triangle   
  
#ifdef DRAWING1
  cout << "  -------------------------------------------- " << endl;
  inquire();
  reffecran();
  Draw();
  penthickness(2);
#endif

   if (!InsertNewPoints(nbvold,NbTSwap)) 
     break;
     
   for (i=nbtold;i<nbt;i++)
     first_np_or_next_t[i]=iter;

  Headt = nbt; // empty list 
  for (i=nbvold;i<nbv;i++) 
    { // for all the triangle contening the vertex i
       Vertex * s  = vertices + i;
       TriangleAdjacent ta(s->t, EdgesVertexTriangle[s->vint][1]);
       Triangle * tbegin= (Triangle*) ta;
       Int4 kt;
       do { 
	 kt = Number((Triangle*) ta);
	 if (first_np_or_next_t[kt]>0) 
	   first_np_or_next_t[kt]=-Headt,Headt=kt;
	 assert( ta.EdgeVertex(0) == s);
	 ta = Next(Adj(ta));
       } while ( (tbegin != (Triangle*) ta)); 
    }   
   
  } while (nbv!=nbvold);
  
  delete []  first_np_or_next_t;

  Int4 NbSwapf =0,NbSwp;
  
  // bofbof 
  
  
  NbSwp = NbSwapf;
  for (i=0;i<nbv;i++)
    NbSwapf += vertices[i].Optim(0);
  /*
  for (i=0;i<nbv;i++)
    NbSwapf += vertices[i].Optim(0);
  for (i=0;i<nbv;i++)
    NbSwapf += vertices[i].Optim(0);
  for (i=0;i<nbv;i++)
    NbSwapf += vertices[i].Optim(0);
  for (i=0;i<nbv;i++)
    NbSwapf += vertices[i].Optim(0);
    */
  NbTSwap +=  NbSwapf ;
  if (verbosity>3) cout << "   " ;
  if (verbosity>2) 
     cout << " Nb Of Vertices =" << nbv << " Nb of triangles = " << nbt-NbOutT 
	  << " NbSwap final = " << NbSwapf << " Nb Total Of Swap = " << NbTSwap << endl;


}

void  Triangles::NewPointsOld(Triangles & Bh)
{ // Triangles::NewPointsOld
  Real8 seuil= 0.7 ;// for two neart point 
  if (verbosity>1) 
    cout << " begin :  Triangles::NewPointsOld " << endl;
  Int4 i,k;
  int j;
  Int4 BeginNewPoint[3];
  Int4 EndNewPoint[3];
#ifdef TRACETRIANGLE
  Int4 trace=0;
#endif
  int step[3];
  Int4 *first_np_or_next_t = new Int4[nbtx];
  Int4 ColorEdge[3];
  Int4 color=-1;
  Triangle *t;
  // generation of the list of next Triangle 
  // at 1 time we test all the triangles
  Int4 Headt =0,next_t;
  for(i=0;i<nbt;i++)
    first_np_or_next_t[i]=-(i+1);
  // end list i >= nbt 
  // the list of test triangle is 
  // the next Triangle on i is  -first_np_or_next_t[i]
  Int4 nbtold,nbvold;

  // Big loop 
  do {
    nbtold = nbt;
    nbvold = nbv;
#ifdef DRAWING1  
  inquire();
#endif  

  // default size of  IntersectionTriangle

  i=Headt;
  next_t=-first_np_or_next_t[i];
  for(t=&triangles[i];i<nbt;t=&triangles[i=next_t],next_t=-first_np_or_next_t[i]) 
    { // for each triangle  t
      // we can change first_np_or_next_t[i]
#ifdef TRACETRIANGLE
      trace =  TRACETRIANGLE <0 ? 1 : i == TRACETRIANGLE;
#endif
      //      cout << " Do the triangle " << i << " Next_t=" << next_t << endl;
      assert(i>=0 && i < nbt );
      first_np_or_next_t[i] = nbv; // to save the fist new point of triangle
      for(j=0;j<3;j++)
	{ // for each edge 
	  TriangleAdjacent tj(t,j);
	  // color++;// the color is 3i+j
          color = 3*i + j ;;
	  ColorEdge[j]=color;
	  BeginNewPoint[j]=nbv;
	  EndNewPoint[j]=nbv-1;
	  step[j]=1;// right sens 
	  Vertex & vA = * tj.EdgeVertex(0);
	  Vertex & vB = * tj.EdgeVertex(1);
	  
#ifdef TRACETRIANGLE
	  if(trace) {
	    cout << " i " << Number(vA) <<" j "<<  Number(vB) 
		 << " "  << t->Locked(j) ;
	  }
#endif
	  if (!t->link) continue;// boundary
	  if (t->det <0) continue;
	  if (t->Locked(j)) continue;

	  TriangleAdjacent tadjj = t->Adj(j);	  
	  Triangle * ta= tadjj;

	  if (ta->det <0) continue;	  
	  
	  R2 A = vA;
	  R2 B = vB;
	  
	  k=Number(ta);
	  // the 2 opposite vertices 
	  const Vertex & vC1  =  *tj.OppositeVertex();
	  const Vertex & vC2 = *tadjj.OppositeVertex();
	  
#ifdef TRACETRIANGLE
	  trace = trace || k == TRACETRIANGLE;
	  if(trace) {
	    cout << "Test Arete " << i << " AB = " << A << B 
		 << "i "  <<Number(vA)<< "j" <<Number(vB); 
	    cout << " link" <<(int)t->link << " ta=" << Number( ta) 
		 << " det " <<ta->det ;
	    cout << " hA = " <<vA.m.h << " hB = " <<vB.m.h ;
	    cout << " loc " << ta->Locked(j) << endl;
	  }
#endif
	  
	  if(first_np_or_next_t[k]>0) { // this edge is done before 
	    // find the color of the edge and begin , end of newpoint
	    register int kk = t->NuEdgeTriangleAdj(j);
	    assert ((*t)(VerticesOfTriangularEdge[j][0]) == (*ta)(VerticesOfTriangularEdge[kk][1]));
	    assert ((*t)(VerticesOfTriangularEdge[j][1]) == (*ta)(VerticesOfTriangularEdge[kk][0]));
	    register Int4 kolor =3*k + kk;
	    ColorEdge[j]=kolor;
	    register Int4 kkk= 1;
	    step[j]=-1;// other sens 
	    BeginNewPoint[j]=0;
	    EndNewPoint[j]=-1; // empty list          
	    for (Int4 iv=first_np_or_next_t[k];iv<nbv;iv++) 
	      if (vertices[iv].color > kolor) 
		break; // the color is passed            
	      else if (vertices[iv].color == kolor) {
		EndNewPoint[j]=iv; 
		if (kkk) // one time test
		  kkk=0,BeginNewPoint[j]=iv;}
	    continue; // next edge of the triangle 
	  } // end  if( k < i) 
       
	  
#ifdef DRAWING1
	  penthickness(2); Move(A);Line(B);   penthickness(1);
#endif
	  const Int4 NbvOld = nbv;
	  lIntTria.SplitEdge(Bh,A,B);
	  // Int4 NbvNp =
	  lIntTria.NewPoints(vertices,nbv,nbvx);
	  Int4 nbvNew=nbv;
	  nbv = NbvOld;
	  for (Int4 iv=NbvOld;iv<nbvNew;iv++) {
	    vertices[nbv].color = color;
	    vertices[nbv].ReferenceNumber=nbv;// circular link
	    R2 C =  vertices[iv].r;
	    vertices[nbv].r =  C;
	    vertices[nbv].m =  vertices[iv].m ;
	    // test if the new point is not to close to the 2 opposite vertex
	    R2 CC1 = C-vC1 , CC2 = C-vC2;
	    if (   (  (vC1.m(CC1) + vertices[nbv].m(CC1)) >  seuil)
		&& (  (vC2.m(CC2) + vertices[nbv].m(CC2)) >  seuil) )
	      nbv++;
	  }
	    
	  EndNewPoint[j] = nbv-1;
	} // end loop for each edge 

#ifdef TRACETRIANGLE
      if(trace) {
	// verification des point cree 
	cout << "\n ------------ " << t->link << " " << t->det 
	     << " b " << BeginNewPoint[0] << " " << BeginNewPoint[1]
	     << " " << BeginNewPoint[2] << " " 
	     << " e " << EndNewPoint[0] << " " << EndNewPoint[1] 
	     << " " << EndNewPoint[2] << " " 
	     << " s " << step[0] << " " << step[1] << " " << step[2] << " " 
	     <<  endl;
      }
#endif

      if (!t->link) continue;// boundary
      if (t->det<=0) continue;// outside 
      //      continue;
      for(int j0=0;j0<3;j0++)
	for (Int4 i0= BeginNewPoint[j0]; i0 <= EndNewPoint[j0];i0++)
	  {
	   // find the neart  point one the opposite edge 
           // to compute i1
	    Vertex & vi0 = vertices[i0];
	    int kstack = 0;
	    Int4 stack[10];
	 //   Int4 savRef[10];
	    int j1 = j0; 
	    while (j0 != (j1 = NextEdge[j1])) {//loop on the 2 other edge
	      // computation of the intersection of edge j1 and DOrto
	      // take the good sens 
	      
	      if (BeginNewPoint[j1]> EndNewPoint[j1])
	        continue; // 
              else if (EndNewPoint[j1] - BeginNewPoint[j1] <1) {
                for (Int4 ii1= BeginNewPoint[j1];ii1<=EndNewPoint[j1];ii1++)
                    stack[kstack++] = ii1;
                 continue;}
	      
                
	      int k0,k1;
	      if (step[j1]<0) k0=1,k1=0; // reverse
	      else k0=0,k1=1; 
	      R2 V10 = (R2)(*t)[VerticesOfTriangularEdge[j1][k0]];
	      R2 V11 = (R2)(*t)[VerticesOfTriangularEdge[j1][k1]];
	      R2 D = V11-V10;
	      Real8 c0 =  vi0.m(D,(R2) vi0);

	      Real8 c10 =  vi0.m(D,V10);
	      Real8 c11 =  vi0.m(D,V11);

	      Real8 s;
	      //cout << " --i0 = " << i0  << D  << V10 << V11 << endl ;
	      //cout << "   c10 " <<  c10 << " c0 " << c0 << " c11 " << c11 << endl;
	      if (( c10 < c0 ) && (c0 < c11)) 
		s = (c11-c0)/(c11-c10);
	      else if  (( c11 < c0 ) && (c0 < c10)) 
		s = (c11-c0) /(c11-c10);
	      else break;
	      R2 VP = V10*s + V11*(1-s);
	      int sss = (c11-c10) >0 ? 1 : -1;
#ifdef DRAWING1
	      penthickness(2);
	      Move((R2) vi0);
	      Line(VP);
	      penthickness(1);
#endif
	      // find the 2 point by dichotomie
	      //cout << "   t =" << Number(t) << " c0 " << c0  ;
	      Int4 ii0 =  BeginNewPoint[j1];
	      Int4 ii1 =  EndNewPoint[j1];	     
              Real8 ciii=-1,cii0=-1,cii1=-1  ;
 	      if ( sss * ((cii0=vi0.m(D,(R2) vertices[ii0]))- c0) >0 )  
		stack[kstack++] = ii0;//cout << " add+0 " << ii0;
 	      else if ( sss * ((cii1=  vi0.m(D ,(R2) vertices[ii1]))- c0) < 0 )  
		stack[kstack++] = ii1;//cout << " add+1 " << ii1;
              else {
	        while ((ii1-ii0)> 1) {
		  Int4 iii = (ii0+ii1)/2;
	          ciii = vi0.m( D ,(R2) vertices[iii]);
		  //cout << " (iii = " << iii << " " <<  ciii << ") ";
		  if ( sss * (ciii - c0) <0 )  ii0 = iii;
		  else ii1 = iii;}	        	      
	        stack[kstack++] = ii0;// cout << " add0 " << ii0;
	        if (ii1 != ii0)  stack[kstack++] = ii1;//cout << " add1 " << ii1;
              }
#ifdef DEBUG2
	      cout << "ii1 = " << ii1 
		   << " ii0 = " << ii0 << endl;
              cout << "   cccc = " << cii0 << " " << ciii 
		   << " " << cii1 << " sss=" << sss << endl;
#endif
	      if (kstack >5) // bug ?
	    	cout << "NewPoints: bug????? " << kstack << " stack  " << stack[kstack]<< endl;
	    }
	    
	    stack[kstack++] = -1; // to stop
	    Int4 i1;
	    kstack =0; 
	    while( (i1=stack[kstack++]) >= 0) 
	      { // the two parameter is i0 and i1 
	      assert(i1 < nbv && i1 >= 0);
	      assert(i0 < nbv && i0 >= 0);
	      assert(i1 != i0);
	      R2 v01 = (R2) vertices[i1]- (R2) vertices[i0];
	      Real8 d01 = (vertices[i0].m(v01) + vertices[i1].m(v01));
		

#ifdef DRAWING1
	      Move(vertices[i0].r);
	      Line(vertices[i1].r);
#endif
#ifdef TRACETRIANGLE
	      if(trace) {
		cout << "\n test j" << j <<" "  << i0 
		     << " " << i1 << " d01=" << d01 <<endl;}
#endif
	      assert (i0 >= nbvold);
	      assert (i1 >= nbvold);
	      assert(i0 != i1);
	      if (d01 == 0) 
		break; 
	      if ( d01 < seuil) 
		if (i1<nbvold) {
		  // remove all the points i0;
		  register Int4 ip,ipp;
		  for  (ip=i0;i0 != (ipp = vertices[ip].ReferenceNumber);ip=ipp)
		    vertices[ip].ReferenceNumber = -1;// mark remove
		  vertices[ip].ReferenceNumber = -1;// mark remove
		}
	      else {
		// remove on of two points
		register Int4 ip0, ip1, ipp0,ipp1;
		register int kk0=1,kk1=1;
		// count the number of common points to compute weight w0,w1
		for  (ip0=i0;i0 != (ipp0 = vertices[ip0].ReferenceNumber);ip0=ipp0) kk0++;
		for  (ip1=i1;i1 != (ipp1 = vertices[ip1].ReferenceNumber);ip1=ipp1) kk1++;
		
		register Real8 w0 = ((Real8) kk0)/(kk0+kk1);
		register Real8 w1 = ((Real8) kk1)/(kk0+kk1);

		// make a circular link
		Exchange(vertices[i0].ReferenceNumber,vertices[i1].ReferenceNumber);
		// the new coordinate 
		R2 C //= vertices[i0] ;
		  =  vertices[i0].r *w0 + vertices[i1].r* w1;

#ifdef TRACETRIANGLE
		if(trace) {
		  cout << "\n ref = "<< vertices[i0].ref << " " <<vertices[i1].ref <<endl;
		}
#endif    
#ifdef DRAWING1  
		Move(vertices[i0].r);
		Line(vertices[i1].r);
		DrawMark(C);
#endif		
		// update the new point points of the list 
		for  (ip0=i0;i0 != (ipp0 = vertices[ip0].ReferenceNumber);ip0=ipp0)
		  vertices[ip0].r = C;	      
		vertices[ip0].r = C;
	      }
	    }
	} // for (i0= ....
  }// for triangle   

#ifdef DRAWING1
  cout << " -------------------------------------------- " << endl;
  inquire();
  reffecran();
  Draw();
  penthickness(2);
#endif
  
  // remove of all the double points

  Int4 ip,ipp,kkk=nbvold;
  for (i=nbvold;i<nbv;i++) 
     if (vertices[i].ReferenceNumber>=0) {// good points
    //  cout <<" i = " << i ;
      for  (ip=i;i != (ipp = vertices[ip].ReferenceNumber);ip=ipp)
         vertices[ip].ReferenceNumber = -1;// mark remove
      vertices[ip].ReferenceNumber = -1;// mark remove
      // cout << i << " ---> " << kkk << endl;        
      vertices[kkk] = vertices[i];
      vertices[kkk].i = toI2(vertices[kkk].r);
#ifdef DRAWING1
      DrawMark(vertices[kkk]);
#endif
      vertices[kkk++].ReferenceNumber = 0;
      
     } 
#ifdef DRAWING1  
  penthickness(1);
#endif
    
 // insertion part --- 

  const Int4 nbvnew = kkk-nbvold;

  cout << "    Remove " << nbv - kkk  << " to close  vertex " ;
  cout << " and   Insert the " <<nbvnew<< " new points " << endl;  
  nbv = kkk;
  Int4 NbSwap=0;
  Icoor2 dete[3];
  
  // construction d'un ordre aleatoire 
  if (! nbvnew) 
    break; 
  if (nbvnew) {
  const Int4 PrimeNumber= AGoodNumberPrimeWith(nbv)  ;
  Int4 k3 = rand()%nbvnew ; 
  for (Int4 is3=0; is3<nbvnew; is3++) 
    ordre[nbvold+is3]= &vertices[nbvold +(k3 = (k3 + PrimeNumber)% nbvnew)];

  for (i=nbvold;i<nbv;i++) 
    { Vertex * vi = ordre[i];
      Triangle *tcvi = FindTriangleContening(vi->i,dete);
      //     Vertex * nv =  quadtree->NearestVertex(vi->i.x,vi->i.y);
      //      cout << " Neart Vertex of "  << Number(vi)<< vi->i << " is " 
      //   << Number(nv) << nv->i  << endl;
      // Int4  kt = Number(tcvi);
      // 
      
      quadtree->Add(*vi); //
#ifdef DRAWING1
      DrawMark(vi->r);
#endif
      assert (tcvi->det >= 0) ;// internal 
      Add(*vi,tcvi,dete),NbSwap += vi->Optim(1);          
    }  
  }
  cout << " Nb swap = " << NbSwap << " after " ;
#ifdef DRAWING1
  inquire();
#endif

  for (i=nbvold;i<nbv;i++) 
    NbSwap += vertices[i].Optim(1);  
  cout << NbSwap << endl;

  for (i=nbtold;i<nbt;i++)
     first_np_or_next_t[i]=1;

  Headt = nbt; // empty list 
  for (i=nbvold;i<nbv;i++) 
    { // for all the triangle contening the vertex i
       Vertex * s  = vertices + i;
       TriangleAdjacent ta(s->t, EdgesVertexTriangle[s->vint][1]);
       Triangle * tbegin= (Triangle*) ta;
       Int4 kt;
       do { 
	 kt = Number((Triangle*) ta);
	 if (first_np_or_next_t[kt]>0) 
	   first_np_or_next_t[kt]=-Headt,Headt=kt;
	 assert( ta.EdgeVertex(0) == s);
	 ta = Next(Adj(ta));
       } while ( (tbegin != (Triangle*) ta)); 
    }


  } while (nbv!=nbvold);
  delete []  first_np_or_next_t;
#ifdef  DEBUG
  int nberr=0;
  for (int it=0;it<nbt;it++)
    if(triangles[it].det==0) 
      if(nberr++<10) cerr << "Bug Triangle nul" << it << triangles[it] << endl;
  if(nberr) MeshError(992,this);
#endif
  cout << " end :  Triangles::NewPoints old  nbv=" << nbv << endl;

}

void Triangles::Insert() 
{
  if (verbosity>2) cout << "  -- Insert initial " << nbv << " vertices " << endl ;
  Triangles * OldCurrentTh =CurrentTh;

  CurrentTh=this;
  double time0=CPUtime(),time1,time2,time3;
  SetIntCoor();
  Int4 i;
  for (i=0;i<nbv;i++) 
    ordre[i]= &vertices[i] ;

  // construction d'un ordre aleatoire 
  const Int4 PrimeNumber= AGoodNumberPrimeWith(nbv) ;
  Int4 k3 = rand()%nbv ; 
  for (int is3=0; is3<nbv; is3++) 
    ordre[is3]= &vertices[k3 = (k3 + PrimeNumber)% nbv];




  for (i=2 ; det( ordre[0]->i, ordre[1]->i, ordre[i]->i ) == 0;) 
    if  ( ++i >= nbv) {
      cerr << " All the vertices are aline " << endl;
      MeshError(998,this); }
          
  // echange i et 2 dans ordre afin 
  // que les 3 premiers ne soit pas aligne
  Exchange( ordre[2], ordre[i]);
  
  // on ajoute un point a l'infini pour construire le maillage
  // afin d'avoir une definition simple des aretes frontieres
  nbt = 2;
  
  
  // on construit un maillage trivale forme
  // d'une arete et de 2 triangles
  // construit avec le 2 aretes orientes et 
  Vertex *  v0=ordre[0], *v1=ordre[1];
  
  triangles[0](0) = 0; // sommet pour infini 
  triangles[0](1) = v0;
  triangles[0](2) = v1;
  
  triangles[1](0) = 0;// sommet pour infini 
  triangles[1](2) = v0;
  triangles[1](1) = v1;
  const int e0 = OppositeEdge[0];
  const int e1 = NextEdge[e0];
  const int e2 = PreviousEdge[e0];
  triangles[0].SetAdj2(e0, &triangles[1] ,e0);
  triangles[0].SetAdj2(e1, &triangles[1] ,e2);
  triangles[0].SetAdj2(e2, &triangles[1] ,e1);
  
  triangles[0].det = -1;  // faux triangles
  triangles[1].det = -1;  // faux triangles
  
  triangles[0].SetTriangleContainingTheVertex();
  triangles[1].SetTriangleContainingTheVertex();
  
  triangles[0].link=&triangles[1];
  triangles[1].link=&triangles[0];
  
#ifdef DEBUG 
  triangles[0].check();
  triangles[1].check();
#endif  
  //  nbtf = 2;
   if (  !quadtree )  quadtree = new QuadTree(this,0);
   quadtree->Add(*v0);
   quadtree->Add(*v1);
   
  // on ajoute les sommets un Ò un 
  Int4 NbSwap=0;

  time1=CPUtime();

    if (verbosity>3) cout << "  -- Begin of insertion process " << endl;

  for (Int4 icount=2; icount<nbv; icount++) {
    Vertex *vi  = ordre[icount];
    //    cout << " Insert " << Number(vi) << endl;
    Icoor2 dete[3];
    Triangle *tcvi = FindTriangleContening(vi->i,dete);
    quadtree->Add(*vi); 
    Add(*vi,tcvi,dete);
    NbSwap += vi->Optim(1,0);
#ifdef DRAWING1
    inquire();
#endif
    
  }// fin de boucle en icount
  time2=CPUtime();
  if (verbosity>3) 
    cout << "    NbSwap of insertion " <<    NbSwap 
	 << " NbSwap/Nbv " <<  (float) NbSwap / (float) nbv 
	 << " NbUnSwap " << NbUnSwap << " Nb UnSwap/Nbv " 
	 << (float)NbUnSwap /(float) nbv 
	 <<endl;
  NbUnSwap = 0;
  // construction d'un ordre aleatoire 
  //  const int PrimeNumber= (nbv % 999983) ? 1000003: 999983 ;
#ifdef NBLOOPOPTIM

  k3 = rand()%nbv ; 
  for (int is4=0; is4<nbv; is4++) 
    ordre[is4]= &vertices[k3 = (k3 + PrimeNumber)% nbv];
  
  double timeloop = time2 ;
  for(int Nbloop=0;Nbloop<NBLOOPOPTIM;Nbloop++) 
    {
      double time000 = timeloop;
      Int4  NbSwap = 0;
      for (int is1=0; is1<nbv; is1++) 
	NbSwap += ordre[is1]->Optim(0,0);
      timeloop = CPUtime();
  if (verbosity>3) 
      cout << "    Optim Loop "<<Nbloop<<" NbSwap: " <<  NbSwap 
	   << " NbSwap/Nbv " 	   <<  (float) NbSwap / (float) nbv 
	   << " CPU=" << timeloop - time000 << "  s, " 
	   << " NbUnSwap/Nbv " << (float)NbUnSwap /(float) nbv  
	   <<  endl;
      NbUnSwap = 0;
      if(!NbSwap) break;
    }
  ReMakeTriangleContainingTheVertex(); 
  // because we break the TriangleContainingTheVertex
#endif
#ifdef  DEBUG
  int nberr=0;
  for (int it=0;it<nbt;it++)
    if(triangles[it].det==0) 
      if(nberr++<10) cerr << "Bug Triangle nul" << it << triangles[it] << endl;
  if(nberr) MeshError(991,this);
#endif
  time3=CPUtime();
  if (verbosity>4) 
  cout << "    init " << time1 - time0 << " initialisation,  " 
       << time2 - time1 << "s, insert point  " 
       << time3 -time2 << "s, optim " << endl
       << "     Init Total Cpu Time = " << time3 - time0 << "s " << endl;
  
#ifdef DRAWING1
  inquire();
#endif  
  CurrentTh=OldCurrentTh;
}


void Triangles::ForceBoundary()
{
  if (verbosity > 2)
    cout << "  -- ForceBoundary  nb of edge " << nbe << endl;
  int k=0;
  Int4  nbfe=0,nbswp=0,Nbswap=0;
  for (Int4 t = 0; t < nbt; t++)  
    if (!triangles[t].det)
      k++,cerr << " det T" << t << " = " << 0 << endl;
  if (k!=0) {
    cerr << " ther is  " << k << "  triangles of mes = 0 " << endl;
    MeshError(11,this);}
  
  TriangleAdjacent ta(0,0);
  for (Int4 i = 0; i < nbe; i++) 
    {
      nbswp =  ForceEdge(edges[i][0],edges[i][1],ta);
      
      if ( nbswp < 0) 	k++;
      else Nbswap += nbswp;
      if (nbswp) nbfe++;
      if ( nbswp < 0 && k < 5)
      {
         cerr << " Missing  Edge " << i << " v0 =  " << Number(edges[i][0]) << edges[i][0].r
      	   <<" v1= " << Number(edges[i][1]) << edges[i][1].r << " " << edges[i].on->Cracked() << "  " << (Triangle *) ta ;
	  if(ta.t)
	  {
	      Vertex *aa = ta.EdgeVertex(0), *bb = ta.EdgeVertex(1);  
	      cerr << " crossing with  [" << Number(aa) << ", " << Number(bb) << "]\n";
	  }
	  else cerr << endl;
	   
      }
      if ( nbswp >=0 && edges[i].on->Cracked())
	  ta.SetCracked();
    }
		  
  
  if (k!=0) {
    cerr << " they is " << k << " lost edges " << endl;
    cerr << " The boundary is crossing may be!" << endl;
    MeshError(10,this);
  }
  for (Int4 j=0;j<nbv;j++)
    Nbswap +=  vertices[j].Optim(1,0);
  if (verbosity > 3)
    cout << "     Nb of inforced edge = " << nbfe << " Nb of Swap " << Nbswap << endl;

}

void Triangles::FindSubDomain(int OutSide=0)
{
    //#define DRAWING1
    
    if (verbosity >2)
    {
	if (OutSide)
	    cout << "  -- Find all external sub-domain ";	
	else
	    cout << "  -- Find all internal sub-domain ";
      if(verbosity>99)
	{
	  
	  for(int i=0;i<nbt;++i)
	      cout << i<< " " << triangles[i] << endl;
	}
	
    }
    // if (verbosity > 4) cout << " OutSide=" << OutSide << endl;
    short * HeapArete = new short[nbt];
    Triangle  **  HeapTriangle = new Triangle*  [nbt];
    Triangle *t,*t1;
    Int4 k,it;
    
    for (Int4 itt=0;itt<nbt;itt++) 
	triangles[itt].link=0; // par defaut pas de couleur
#ifdef DRAWING1
    reffecran();
#endif
    
    Int4  NbSubDomTot =0;
    for ( it=0;it<nbt;it++)  { 
	if ( ! triangles[it].link  ) {
	    t = triangles + it;
	    NbSubDomTot++;; // new composante connexe
	    Int4 i = 0; // niveau de la pile 
	    t->link = t ; // sd forme d'un triangle cicular link
#ifdef DRAWING1
	    t->Draw(NbSubDomTot-1);
#endif
	    
	    
	    HeapTriangle[i] =t ; 
	    HeapArete[i] = 3;
	    
	    while (i >= 0) // boucle sur la pile
	    { while ( HeapArete[i]--) // boucle sur les 3 aretes 
	    { 
		int na =  HeapArete[i];
		Triangle * tc =  HeapTriangle[i]; // triangle courant
		if( ! tc->Locked(na)) // arete non frontiere
		{
		    Triangle * ta = tc->TriangleAdj(na) ; // næ triangle adjacent
		    if (ta->link == 0 ) // non deja chainer => on enpile
		    { 
			i++;
#ifdef DRAWING1
			ta->Draw(NbSubDomTot-1);
#endif
			ta->link = t->link ;  // on chaine les triangles
			t->link = ta ;  // d'un meme sous domaine          
			HeapArete[i] = 3; // pour les 3 triangles adjacents
			HeapTriangle[i] = ta;
		    }}
	    } // deplie fin de boucle sur les 3 adjacences
		i--;
	    }          
	}      
    }
    
    // supression de tous les sous domaine infini <=>  contient le sommet NULL
    it =0;
    NbOutT = 0;
    while (it<nbt) {
	if (triangles[it].link) 
	{ 
	    if (!( triangles[it](0) &&  triangles[it](1) &&  triangles[it](2) )) 
	    {
		// infini triangle 
		NbSubDomTot --;
		//  cout << " triangle infini " << it << triangles[it] << endl;
		t=&triangles[it];
		NbOutT--;  // on fait un coup de trop. 
		while  (t){ // cout << Number(t) << " " << endl;
		    NbOutT++;
		    t1=t;
		    t=t->link;
		    t1->link=0;}//while (t)
	    }
	}   
	it++;} // end while (it<nbt)
    if (nbt == NbOutT ||  !NbSubDomTot) 
    {
	cout << "\n error : " <<  NbOutT << " " << NbSubDomTot <<" " << nbt << endl;
	cerr << "Error: The boundary is not close => All triangles are outside " << endl;
	MeshError(888,this);
    }
    
    delete [] HeapArete;
    delete [] HeapTriangle;
    
    
    if (OutSide|| !Gh.subdomains || !Gh.NbSubDomains ) 
    { // No geom sub domain
	Int4 i;
	if (subdomains) delete [] subdomains;
	subdomains = new SubDomain[ NbSubDomTot];
	NbSubDomains=  NbSubDomTot;
	for ( i=0;i<NbSubDomains;i++) {
	    subdomains[i].head=NULL;
	    subdomains[i].ref=i+1;
	}
	Int4 * mark = new Int4[nbt];
	for (it=0;it<nbt;it++)
	    mark[it]=triangles[it].link ? -1 : -2;
	
	it =0;
	k = 0;
	while (it<nbt) {
	    if (mark[it] == -1) {
		t1 = & triangles[it];
		t = t1->link;
		mark[it]=k;
#ifdef DRAWING1  
		t1->Draw(k);
#endif
		subdomains[k].head = t1;
		// cout << " New -- " << Number(t1) << " " << it << endl;
		do {// cout << " k " << k << " " << Number(t) << endl;
		    mark[Number(t)]=k;
#ifdef DRAWING1  
		    t->Draw(k);
#endif
		    t=t->link;
		} while (t!=t1);
#ifdef DRAWING1  
		t1->Draw(k);
#endif
		mark[it]=k++;}
	    //    else if(mark[it] == -2 ) triangles[it].Draw(999);
	    it++;} // end white (it<nbt)
	assert(k== NbSubDomains);
	if(OutSide) 
	{
	    //  to remove all the sub domain by parity adjacents
	    //  because in this case we have only the true boundary edge
	    // so teh boundary is manifold
	    Int4 nbk = NbSubDomains;
	    while (nbk)
		for (it=0;it<nbt && nbk ;it++)
		    for (int na=0;na<3 && nbk ;na++)
		    {
			Triangle *ta = triangles[it].TriangleAdj(na);
			Int4 kl = ta ? mark[Number(ta)] : -2;
			Int4 kr = mark[it];
			if(kr !=kl) {
			    //cout << kl << " " << kr << " rl "  << subdomains[kl].ref
			    // << " rr " << subdomains[kr].ref ;
			    if (kl >=0 && subdomains[kl].ref <0 && kr >=0 && subdomains[kr].ref>=0)
				nbk--,subdomains[kr].ref=subdomains[kl].ref-1;
			    if (kr >=0 && subdomains[kr].ref <0 && kl >=0 && subdomains[kl].ref>=0)
				nbk--,subdomains[kl].ref=subdomains[kr].ref-1;
			    if(kr<0 && kl >=0 && subdomains[kl].ref>=0)
				nbk--,subdomains[kl].ref=-1;
			    if(kl<0 && kr >=0 && subdomains[kr].ref>=0)
				nbk--,subdomains[kr].ref=-1;
			    //   cout << " after \t "   
			    //	 << kl << subdomains[kl].ref << " rr " << kr 
			    // << subdomains[kr].ref << endl;
			}
		    }
			//  cout << subdomains[0].ref << subdomains[1].ref << endl;
			Int4  j=0;
	    for ( i=0;i<NbSubDomains;i++)
		if((-subdomains[i].ref) %2) { // good 
					      //cout << " sudom ok  = " << i << " " << subdomains[i].ref
					      // << " " << (-subdomains[i].ref) %2 << endl;
		    if(i != j) 
			Exchange(subdomains[i],subdomains[j]);
		    j++;}
		    else
		    { //cout << " remove sub domain " << i << endl;
			t= subdomains[i].head;
			while  (t){// cout << Number(t) << " " << endl;
			    NbOutT++;
			    t1=t;
			    t=t->link;
			    t1->link=0;}//while (t)
		    }
		    
		    if(verbosity>4)
			cout << " Number of remove sub domain (OutSideMesh) =" << NbSubDomains-j << endl;
	    NbSubDomains=j;
	}
	
	delete []  mark; 
	
    }
    else
    { // find the head for all sub domaine
	if (Gh.NbSubDomains != NbSubDomains && subdomains)
	    delete [] subdomains, subdomains=0;
	if (! subdomains  ) 
	    subdomains = new SubDomain[ Gh.NbSubDomains];
	NbSubDomains =Gh.NbSubDomains;
	if(verbosity>4)
	    cout << "     find the " << NbSubDomains << " sub domain " << endl;
	Int4 err=0;
	ReMakeTriangleContainingTheVertex();
	Int4 * mark = new Int4[nbt];
	Edge **GeometricalEdgetoEdge = MakeGeometricalEdgeToEdge();
        
	for (it=0;it<nbt;it++)
	    mark[it]=triangles[it].link ? -1 : -2;
	Int4 inew =0;
	for (Int4 i=0;i<NbSubDomains;i++) 
	{
	    GeometricalEdge &eg = *Gh.subdomains[i].edge;
	    subdomains[i].ref = Gh.subdomains[i].ref;
	    int ssdlab = subdomains[i].ref;
	    // by carefull is not easy to find a edge create from a GeometricalEdge 
	    // see routine MakeGeometricalEdgeToEdge
	    Edge &e = *GeometricalEdgetoEdge[Gh.Number(eg)];
	    assert(&e);
	    Vertex * v0 =  e(0),*v1 = e(1);
	    Triangle *t  = v0->t;
	    int sens = Gh.subdomains[i].sens;
	    // test if ge and e is in the same sens 
	    //	cout << " geom edge = " <<  Gh.Number(eg) <<" @" << &eg << " ref = " << subdomains[i].ref 
	    //     << " ref edge =" << eg.ref << " sens " << sens ;
	    if (((eg[0].r-eg[1].r),(e[0].r-e[1].r))<0)
		sens = -sens ;
	    subdomains[i].sens = sens;
	    subdomains[i].edge = &e;
	    //	cout << " sens " << sens << " in geom " << eg[0].r << eg[1].r << " in mesh  " << e[0].r << e[1].r << endl;
	    //	cout << "  v0 , v1 = " << Number(v0) << " "  << Number(v1) << endl;
	    assert(t && sens);
	    
	    TriangleAdjacent  ta(t,EdgesVertexTriangle[v0->vint][0]);// previous edges
		
	    while (1) 
	    {
		assert( v0 == ta.EdgeVertex(1) );
		//	 cout << " recherche " << Number( ta.EdgeVertex(0)) << endl;
		if (ta.EdgeVertex(0) == v1) { // ok we find the edge
		    if (sens>0)  
			subdomains[i].head=t=Adj(ta);
		    else 
			subdomains[i].head=t=ta;
		    //cout << "      triangle  =" << Number(t) << " = " << (*t)[0].r <<  (*t)[1].r <<  (*t)[2].r << endl;
		    if(t<triangles || t >= triangles+nbt || t->det < 0 
		        || t->link == 0) // Ajoute aout 200 
		     {
			cerr << " Error in the def of sub domain "<<i
			     << " form border " << NbSubDomains - i  << "/" << NbSubDomains
			     << ": Bad sens  " << Gh.Number(eg) <<" "<< sens <<  endl;  
			err++;
			break;}
		    Int4 it = Number(t);
		    if (mark[it] >=0) {
			if(verbosity>10)
			    cerr << "     Warning: the sub domain " << i << " ref = " << subdomains[i].ref 
				<< " is previouly defined with "  <<mark[it] << " ref = " << subdomains[mark[it]].ref
				<< " skip this def " << endl;
			break;}
		    if(i != inew) 
			Exchange(subdomains[i],subdomains[inew]);
		    inew++;
		    Triangle *tt=t;
		    Int4 kkk=0;
		    do 
		    {
			kkk++;
			assert(mark[Number(tt)]<0);
#ifdef DRAWING1
			tt->Draw(i);
#endif
			mark[Number(tt)]=i;
			tt=tt->link;
		    } while (tt!=t);
		    if(verbosity>7)
			cout << "     Nb de triangles dans le sous domaine " << i << " de ref " << subdomains[i].ref << " = " << kkk << endl;
		    break;}
		ta = Previous(Adj(ta));         
		if(t == (Triangle *) ta) {
		    err++;
		    cerr << " Error in the def of sub domain " << i 
			<< " edge=" << Gh.Number(eg) << " " << sens << endl;
		    break;}
		//         cout << " NB of remove subdomain " << NbSubDomTot-NbSubDomains<< endl;
		
	    }
		
	}
	if (err) MeshError(777,this);
	
	if (inew < NbSubDomains) {
	    if (verbosity>5) 
		cout << "     Warning: We remove " << NbSubDomains-inew << " SubDomains " << endl;
	    NbSubDomains=inew;}
	
	
	for (it=0;it<nbt;it++)
	    if ( mark[it] ==-1 ) 
		NbOutT++,triangles[it].link =0;
	delete [] GeometricalEdgetoEdge;
	delete [] mark;
	
    }
    
#ifdef DRAWING1
    inquire();
#endif
    NbOutT=0;
    for (it=0;it<nbt;it++) 
	if(!triangles[it].link)  NbOutT++;
    if (verbosity> 4)
	cout << "    " ;
    if (verbosity> 2)
	cout << " Nb of Sub borned Domain  = " <<  NbSubDomTot << " NbOutTriangles = " << NbOutT <<endl;
#ifdef DRAWING1
    inquire();
#endif
    
    //#undef DRAWING1
    
    
}
void Triangles::ReNumberingVertex(Int4 * renu)
{
  // warning be carfull because pointeur
  // from on mesh to over mesh 
  //  --  so do ReNumbering a the beginning
  Vertex * ve = vertices+nbv;
  Int4 it,ie,i;
  
  for ( it=0;it<nbt;it++) 
    triangles[it].ReNumbering(vertices,ve,renu);
  
  for ( ie=0;ie<nbe;ie++) 
    edges[ie].ReNumbering(vertices,ve,renu);
  
  for (i=0;i< NbVerticesOnGeomVertex;i++)
    {
      Vertex *v = VerticesOnGeomVertex[i].mv;
      if (v>=vertices && v < ve)
	VerticesOnGeomVertex[i].mv=vertices+renu[Number(v)];
    }
  
  for (i=0;i< NbVerticesOnGeomEdge;i++)
    {
      Vertex *v =VerticesOnGeomEdge[i].mv;
      if (v>=vertices && v < ve)
	 VerticesOnGeomEdge[i].mv=vertices+renu[Number(v)];
    }

  for (i=0;i< NbVertexOnBThVertex;i++)
    {
      Vertex *v=VertexOnBThVertex[i].v;
      if (v>=vertices && v < ve)
	VertexOnBThVertex[i].v=vertices+renu[Number(v)];
    }
  
  for (i=0;i< NbVertexOnBThEdge;i++)
    {
      Vertex *v=VertexOnBThEdge[i].v;
      if (v>=vertices && v < ve)
	VertexOnBThEdge[i].v=vertices+renu[Number(v)];
    }

  // move the Vertices without a copy of the array 
  // be carefull not trivial code 
  Int4 j;
  for ( it=0;it<nbv;it++) // for all sub cycles of the permutation renu
    if (renu[it] >= 0) // a new sub cycle
      { 
        i=it;
        Vertex ti=vertices[i],tj;
        while ( (j=renu[i]) >= 0) 
         { // i is old, and j is new 
           renu[i] = -1-renu[i]; // mark 
           tj = vertices[j]; // save new
           vertices[j]= ti; // new <- old
           i=j;     // next 
           ti = tj;
         }  
     }
  if (quadtree) 
  {  delete quadtree;
   quadtree = new QuadTree(this);
  }
 for ( it=0;it<nbv;it++)
   renu[i]= -renu[i]-1;
  
}
void Triangles::ReNumberingTheTriangleBySubDomain(bool justcompress)
 {
  Int4 *renu= new Int4[nbt];
  register Triangle *t0,*t,*te=triangles+nbt;
  register Int4 k=0,it,i,j;
      
  for ( it=0;it<nbt;it++) 
    renu[it]=-1; // outside triangle 
  for ( i=0;i<NbSubDomains;i++)
   { 
     t=t0=subdomains[i].head;
     assert(t0); // no empty sub domain
     do { 
       Int4 kt = Number(t);
       assert(kt>=0 && kt < nbt );
       assert(renu[kt]==-1);
       renu[kt]=k++;
       }
     while (t0 != (t=t->link));
   }
  if (verbosity>9)
    cout << " number of inside triangles " << k << " nbt = " << nbt << endl;
  // take is same numbering if possible    
  if(justcompress)
      for ( k=0,it=0;it<nbt;it++) 
     if(renu[it] >=0 ) 
       renu[it]=k++;
   
   // put the outside triangles at the end
   for ( it=0;it<nbt;it++) 
    if (renu[it]==-1) 
      renu[it]=k++;
      
    assert(k == nbt);
   // do the change on all the pointeur 
   for ( it=0;it<nbt;it++)
     triangles[it].ReNumbering(triangles,te,renu);
     
   for ( i=0;i<NbSubDomains;i++)
     subdomains[i].head=triangles+renu[Number(subdomains[i].head)];
   
  // move the Triangles  without a copy of the array 
  // be carefull not trivial code 
  for ( it=0;it<nbt;it++) // for all sub cycles of the permutation renu
    if (renu[it] >= 0) // a new sub cycle
      { 
        i=it;
        Triangle ti=triangles[i],tj;
        while ( (j=renu[i]) >= 0) 
         { // i is old, and j is new 
           renu[i] = -1; // mark 
           tj = triangles[j]; // save new
           triangles[j]= ti; // new <- old
           i=j;     // next 
           ti = tj;
         }  
     }
   delete [] renu;
   nt = nbt - NbOutT;

#ifdef DEBUG
// verif 
     for ( it=0;it<nbt;it++)
      triangles[it].check();
#endif   
 }
Int4  Triangles::ConsRefTriangle(Int4 *reft) const
{
  assert(reft);
  register Triangle *t0,*t;
  register Int4 k=0, num;   
  for (Int4 it=0;it<nbt;it++) 
    reft[it]=-1; // outside triangle 
  for (Int4 i=0;i<NbSubDomains;i++)
   { 
     t=t0=subdomains[i].head;
     assert(t0); // no empty sub domain
     // register Int4 color=i+1;// because the color 0 is outside triangle
     do { k++;
       num = Number(t);
       assert(num>=0 &&num < nbt);
       reft[num]=i;
       //  cout << Number(t0) << " " <<Number(t)<< " "  << i << endl;
       }
     while (t0 != (t=t->link));
   }
  //  NbOutT = nbt - k;
  if (verbosity>5) 
    cout << " Nb of Sub Domain =" << NbSubDomains  << " Nb of In Triangles " << k 
	 << " Nbt = " << nbt << " Out Triangles = " << nbt - k <<  endl;
   
  return k;   

}

/*
void Triangles::ConsLinkTriangle()
{
  for (Int4 i=0;i<NbSubDomains;i++)
    subdomains[i].head=0;
  register Triangle * t=triangles,*tend = triangles+nbt,*hst;
  for (;t<tend;t++)
   {  
       if (t->link) 
        {
          register Int4 color = t->color-1;
          assert(color<NbSubDomains && color>=0);
          if (hst=subdomains[color].head) {
            t->link=hst->link;
            hst->link=t; }
          else {
            subdomains[color].head = t;
            t->link=t;}// circular link         
        }
     }
 {
   for (Int4 i=0;i<NbSubDomains;i++)
     assert(subdomains[i].head);
 }
   
}
*/
/* void Triangles::RandomInit() 
{ 
 //  Meshbegin = vertices;
 //  Meshend  = vertices + nbvx;
   nbv = nbvx;
   for (int i = 0; i < nbv; i++)
     {
        vertices[i].r.x= rand();
        vertices[i].r.y= rand();
        vertices[i].ref = 0;
     }
}
void Triangles::CubeInit(int n,int m) 
{ 
  //  nbvx = n*m;
 //  Meshbegin = vertices;
  // Meshend  = vertices + nbvx;
   nbv = n*m;
   assert(nbv <= nbvx);
   int k =0;
   for (int i = 0; i < n; i++)
   for (int j = 0; j < m; j++)
     {
        vertices[k].r.x= i;
        vertices[k].r.y= j;
        vertices[k++].ref = 0;
     }
}
*/
Vertex * Triangles::NearestVertex(Icoor1 i,Icoor1 j) 
   { return  quadtree->NearestVertex(i,j); } 


void Triangles::PreInit(Int4 inbvx,char *fname)
{
  srand(19999999);
  OnDisk =0;
  NbRef=0;
  //  allocGeometry=0;
  identity=0;
  NbOfTriangleSearchFind =0;
  NbOfSwapTriangle =0;
  nbiv=0;
  nbv=0;
  nbvx=inbvx;
  nbt=0;
  NbOfQuad = 0;
  nbtx=2*inbvx-2;
  NbSubDomains=0;
  NbVertexOnBThVertex=0;
  NbVertexOnBThEdge=0;
  VertexOnBThVertex=0;
  VertexOnBThEdge=0;
  
  NbCrackedVertices=0;
  NbCrackedEdges =0;
  CrackedEdges  =0;  
  nbe = 0; 
  name = fname ;

  if (inbvx) {
    vertices=new Vertex[nbvx];
    assert(vertices);
    ordre=new (Vertex* [nbvx]);
    assert(ordre);
    triangles=new Triangle[nbtx];
    assert(triangles);}
  else {
    vertices=0;
    ordre=0;
    triangles=0;
    nbtx=0;
   }
  if ( name || inbvx)
    {
      time_t timer =time(0);
      char buf[70];     
      strftime(buf ,70,", Date: %y/%m/%d %H:%M %Ss",localtime(&timer));
      counter++; 
      char countbuf[30];   
      sprintf(countbuf,"%d",counter);
      int lg =0 ;
      if (&BTh != this && BTh.name)
	lg = strlen(BTh.name)+4;
      identity = new char[ lg + strlen(buf) + strlen(countbuf)+ 2  + 10 + ( Gh.name ? strlen(Gh.name) + 4 : 0)];
      identity[0]=0;
      if (lg)
	strcat(strcat(strcat(identity,"B="),BTh.name),", ");

      if (Gh.name)
	strcat(strcat(identity,"G="),Gh.name);
      strcat(strcat(identity,";"),countbuf);
      strcat(identity,buf);
      // cout << "New MAILLAGE "<< identity << endl;
    } 
  
  quadtree=0;
//  edgescomponante=0;
  edges=0;
  VerticesOnGeomVertex=0;
  VerticesOnGeomEdge=0;
  NbVerticesOnGeomVertex=0;
  NbVerticesOnGeomEdge=0;
//  nbMaxIntersectionTriangles=0;
//  lIntTria;
  subdomains=0;
  NbSubDomains=0;
//  Meshbegin = vertices;
//  Meshend  = vertices + nbvx;
  if (verbosity>98) 
    cout << "Triangles::PreInit() " << nbvx << " " << nbtx 
	 << " " << vertices 
	 << " " << ordre << " " <<  triangles <<endl;

}

void Triangles::GeomToTriangles1(Int4 inbvx,int KeepBackVertices) 
{ 
//#define DRAWING1
  Gh.NbRef++;// add a ref to Gh

  
/************************************************************************* 
// methode in 2 step
// 1 - compute the number of new edge to allocate
// 2 - construct the edge
 remark: 
  in this part we suppose to have a background mesh with the same
  geometry 
  
  To construct the discretisation of the new mesh we have to 
  rediscretize the boundary of background Mesh 
  because we have only the pointeur from the background mesh to the geometry.
  We need the abcisse of the background mesh vertices on geometry
  so a vertex is 
  0 on GeometricalVertex ;
  1 on GeometricalEdge + abcisse
  2 internal 
  
  *************************************************************************/
  assert(&BTh.Gh == &Gh);

 // if(verbosity==100) Gh.Write("/tmp/gg.gmsh");
  BTh.NbRef++; // add a ref to BackGround Mesh
  PreInit(inbvx);
  BTh.SetVertexFieldOn();
#ifdef DRAWING
  if (withrgraphique)
   { BTh.InitDraw();
  reffecran(); 
  CurrentTh = this;}
#endif
  int * bcurve = new int[Gh.NbOfCurves]; // 
  
  // we have 2 ways to make the loop 
  // 1) on the geometry 
  // 2) on the background mesh
  //  if you do the loop on geometry, we don't have the pointeur on background,
  //  and if you do the loop in background we have the pointeur on geometry
  // so do the walk on  background
  //  Int4 NbVerticesOnGeomVertex;
  //  VertexOnGeom * VerticesOnGeomVertex;  
  //  Int4 NbVerticesOnGeomEdge;
  //  VertexOnGeom * VerticesOnGeomEdge;
  
  
  NbVerticesOnGeomVertex=0;
  NbVerticesOnGeomEdge=0;
  //1 copy of the  Required vertex
  int i; 
  for ( i=0;i<Gh.nbv;i++)
    if (Gh[i].Required()) NbVerticesOnGeomVertex++;
  
  VerticesOnGeomVertex = new VertexOnGeom[NbVerticesOnGeomVertex];
  VertexOnBThVertex = new VertexOnVertex[NbVerticesOnGeomVertex];
  //
  if( NbVerticesOnGeomVertex >= nbvx) 
    {
      cerr << " Too much vertices on geometry " << NbVerticesOnGeomVertex << " >= " << nbvx << endl; 
      MeshError(1,this);
    }
  assert(vertices);
  for (i=0;i<Gh.nbv;i++)
    if (Gh[i].Required()) {//Gh  vertices Required
        vertices[nbv] =  Gh[i];
        vertices[nbv].i = I2(0,0);
        Gh[i].to = vertices + nbv;// save Geom -> Th
        VerticesOnGeomVertex[nbv]= VertexOnGeom(vertices[nbv],Gh[i]);
	// cout << "--------- "  <<Number(Gh[i].to) << " " << Gh[i].to << " " << i << endl;
	nbv++;}
    else Gh[i].to=0;
  // 
  for (i=0;i<BTh.NbVerticesOnGeomVertex;i++)
    { 
      VertexOnGeom & vog = BTh.VerticesOnGeomVertex[i];
      if (vog.IsRequiredVertex()) {
	 GeometricalVertex * gv = vog;
	 Vertex *bv = vog;
	 assert(gv->to);// use of Geom -> Th
	 VertexOnBThVertex[NbVertexOnBThVertex++] = VertexOnVertex(gv->to,bv);
	 gv->to->m = bv->m; // for taking the metrix of the background mesh
	 ;}
    }
  assert(NbVertexOnBThVertex == NbVerticesOnGeomVertex);
// new stuff FH with curve
//  find the begin of the curve in BTh
{
  Gh.UnMarkEdges();	
  int bfind=0;
/*
   cout << " nb curves = " << Gh.NbOfCurves << endl;
   for(int i=0;i<Gh.NbOfCurves ;i++)
     {
        cout << " Curve " << i << " begin e=" << Gh.Number(Gh.curves[i].be) << " k=" << Gh.curves[i].kb 
             << "  end e= " << Gh.Number(Gh.curves[i].ee) << " k=" << Gh.curves[i].ke << endl;
     }*/
  for (int i=0;i<Gh.NbOfCurves;i++)
   {
    bcurve[i]=-1; 
   }
  
  for (int iedge=0;iedge<BTh.nbe;iedge++) 
    {      
      Edge & ei = BTh.edges[iedge];
      for(int je=0;je<2;je++) // for the 2 extremites
	 if (!ei.on->Mark() && ei[je].on->IsRequiredVertex() )
	    {
	      // a begin of curve 
	      int nc = ei.on->CurveNumber;
	      
	      //cout << "curve " <<  nc << " v " << Gh.Number((GeometricalVertex *) *ei[je].on) << " "
	      //     << " e "  << " " << Gh.Number(ei.on) << " vc " << Gh.Number((*Gh.curves[nc].be)[Gh.curves[nc].kb]) << endl;
	      if(
	         ei.on==Gh.curves[nc].be    && 
	         (GeometricalVertex *) *ei[je].on == &(*Gh.curves[nc].be)[Gh.curves[nc].kb] //  same extremity 
	         )     
	        { 
	        // cout << " find " << endl;
	         bcurve[nc]=iedge*2+je;
	         bfind++;	
	        }      
            }
    } 
   assert( bfind==Gh.NbOfCurves);
}          
// method in 2 + 1 step 
//  0.0) compute the length and the number of vertex to do allocation
//  1.0)  recompute the length
//  1.1)   compute the  vertex 
  Int4 nbex=0,NbVerticesOnGeomEdgex=0;
  for (int step=0; step <2;step++)
    {
      Int4 NbOfNewPoints=0;
      Int4 NbOfNewEdge=0;
      Int4 iedge;
      Gh.UnMarkEdges();	
/*   add Curve loop  FH    
      // find a starting back groud edges to walk 
      for (iedge=0;iedge<BTh.nbe;iedge++) {      
	    Edge & ei = BTh.edges[iedge];
	    for(int jedge=0;jedge<2;jedge++) // for the 2 extremites
	     if (!ei.on->Mark() && ei[jedge].on->IsRequiredVertex() )
	    {
*/  
// new code FH 2004 
       Real8 L=0;
       for (int icurve=0;icurve<Gh.NbOfCurves;icurve++)
            { 
              iedge=bcurve[icurve]/2;
              int jedge=bcurve[icurve]%2;
              if( ! Gh.curves[icurve].master) continue; // we skip all equi curve
	      Edge & ei = BTh.edges[iedge];
	      // warning: ei.on->Mark() can be change in
	      // loop for(jedge=0;jedge<2;jedge++) 
	      // new curve  
	      // good the find a starting edge 
	      Real8 Lstep=0,Lcurve=0;// step between two points   (phase==1) 
	      Int4 NbCreatePointOnCurve=0;// Nb of new points on curve     (phase==1) 

	      
	      //    cout.precision(16);
	      for(int phase=0;phase<=step;phase++) 
		{

		  for(Curve * curve= Gh.curves+icurve;curve;curve= curve->next)
		     {

		    int icurveequi= Gh.Number(curve);
  
                   if( phase == 0 &&  icurveequi != icurve)  continue;

                  int k0=jedge,k1;
                  Edge * pe=  BTh.edges+iedge;
		  //GeometricalEdge *ong = ei.on;
                   int iedgeequi=bcurve[icurveequi]/2;
                   int jedgeequi=bcurve[icurveequi]%2;

                  int k0equi=jedgeequi,k1equi;		  
                  Edge * peequi= BTh.edges+iedgeequi;
		  GeometricalEdge *ongequi = peequi->on;
		  
                  Real8 sNew=Lstep;// abcisse of the new points (phase==1) 
                  L=0;// length of the curve
                  Int4 i=0;// index of new points on the curve
                  register GeometricalVertex * GA0 = *(*peequi)[k0equi].on;
                  Vertex *A0;
                  A0 = GA0->to;  // the vertex in new mesh
                  Vertex *A1;
		  VertexOnGeom *GA1;
                  Edge * PreviousNewEdge = 0;
		  //  cout << "  --------------New Curve phase " << phase 
		  //       << "---------- A0=" << *A0 << ei[k0]  <<endl;
                  assert (A0-vertices>=0 && A0-vertices <nbv);
                  if(ongequi->Required() ) 
                        {
                         GeometricalVertex *GA1 = *(*peequi)[1-k0equi].on;
                         A1 = GA1->to;  //
                        }       
                  else 
                  for(;;) 
		    {
		      //   assert(pe && BTh.Number(pe)>=0 && BTh.Number(pe)<=BTh.nbe);
		      Edge &ee=*pe; 
		      Edge &eeequi=*peequi; 
		      k1 = 1-k0; // next vertex of the edge 
		      k1equi= 1 - k0equi;
		      
		      assert(pe  && ee.on);
		      ee.on->SetMark();
		      Vertex & v0=ee[0], & v1=ee[1];
		      R2 AB= (R2) v1 - (R2) v0;
		      Real8 L0=L,LAB;
		      LAB =  LengthInterpole(v0.m,v1.m,AB);
		      L+= LAB;    
		      if (phase) {// computation of the new points
			while ((i!=NbCreatePointOnCurve) && sNew <= L) { 
			  //    cout  << " L0= " << L0 << " L " << L << " sN=" 
			  //         << sNew << " LAB=" << LAB << " NBPC =" <<NbCreatePointOnCurve<< " i " << i  << endl;
			  assert (sNew >= L0);
			  assert(LAB);
			  
			  
			  assert(vertices && nbv<nbvx);
			  assert(edges && nbe < nbex);
			  assert(VerticesOnGeomEdge && NbVerticesOnGeomEdge < NbVerticesOnGeomEdgex);
			  // new vertex on edge
			  A1=vertices+nbv++;
			  GA1=VerticesOnGeomEdge+NbVerticesOnGeomEdge;
			  Edge *e = edges + nbe++;
			  Real8 se= (sNew-L0)/LAB;
			  assert(se>=0 && se < 1.000000001);
#ifdef DEBUG
			  se =  abscisseInterpole(v0.m,v1.m,AB,se); // because code \ref(xxx)
#else
			  se =  abscisseInterpole(v0.m,v1.m,AB,se,1);
#endif
			  assert(se>=0 && se <= 1);
			  //((k1==1) != (k1==k1equi))
			  se = k1 ? se : 1. - se;
			  se = k1==k1equi ? se : 1. - se;
			  VertexOnBThEdge[NbVerticesOnGeomEdge++] = VertexOnEdge(A1,&eeequi,se); // save 
			  ongequi = Gh.ProjectOnCurve(eeequi,se,*A1,*GA1); 
			  A1->ReferenceNumber = eeequi.ref;
			  A1->DirOfSearch =NoDirOfSearch;
			  //cout << icurveequi << " " << i << " " <<  *A1 << endl;
			  e->on = ongequi;
			  e->v[0]=  A0;
			  e->v[1]=  A1;
			  if(verbosity>99)
				cout << i << "+ New P "<< nbv-1 << " "  <<sNew<< " L0=" << L0 
				<< " AB=" << LAB << " s=" << (sNew-L0)/LAB << " se= "  
				<< se <<" B edge " << BTh.Number(ee) << " signe = " << k1 <<" " << A1->r <<endl;
			    
#ifdef DEBUG
			  // code \label(xxx)
			    R2  A1A0 = A1->r - A0->r;
			    Real8 dp = LengthInterpole(A0->m,A1->m,A1A0);
			    if (dp > 1.4) {
			      cerr << " PB new Points "<< nbv-1 ;
			      cerr << " AB=" << LAB << " s=" << (sNew-L0)/LAB << " se= "  ;
			      cerr << se <<" B edge " << BTh.Number(ee) << " signe = " << k1 <<endl;	      
			      cerr << " PB calcul new on cuver points trop loin l=" << dp 
				   << " v=" << nbv-1 << " " << nbv-2 << " Lcurve = " << Lcurve << AB <<v0.m<< v1.m <<endl;
			    }

#endif
			  e->ref = eeequi.ref;
			  e->adj[0]=PreviousNewEdge;
			  
			  if (PreviousNewEdge)
			    PreviousNewEdge->adj[1] =  e;
			  PreviousNewEdge = e;
#ifdef DRAWING1
			  e->Draw();
			  //         A0->Draw();
			  A1->m.Draw(*A1);
			  A1->Draw(Number(A1));
#endif
			  A0=A1;
			  sNew += Lstep;
			  //   cout << " sNew = " << sNew << " L = " << L 
			  //        << "  ------" <<NbCreatePointOnCurve << " " << i <<  endl;
			  if (++i== NbCreatePointOnCurve) break;
			}
			
		      }               
		      assert(ee.on->CurveNumber==ei.on->CurveNumber);
		      if(verbosity>98) cout <<  BTh.Number(ee) << " " << " on=" << *ee[k1].on << " "<< ee[k1].on->IsRequiredVertex() <<  endl;
		      if ( ee[k1].on->IsRequiredVertex()) {
		         assert(eeequi[k1equi].on->IsRequiredVertex());
			register GeometricalVertex * GA1 = *eeequi[k1equi].on;
			A1=GA1->to;// the vertex in new mesh
			assert (A1-vertices>=0 && A1-vertices <nbv);
			break;}
		      if (!ee.adj[k1])
			{cerr << "Error adj edge " << BTh.Number(ee) << ", nbe = "  << nbe 
			      << " Gh.vertices " << Gh.vertices 
			      << " k1 = " << k1 << " on=" << *ee[k1].on << endl;
			cerr << ee[k1].on->gv-Gh.vertices << endl;
			}
		      pe = ee.adj[k1]; // next edge
		      k0 = pe->Intersection(ee); 
		      peequi= eeequi.adj[k1equi];  // next edge
		      k0equi=peequi->Intersection(eeequi);            
		    }// for(;;) end of the curve
		  

		  if (phase) // construction of the last edge
                    {
		      Edge *e = edges + nbe++;
		      if (verbosity>10) 
			cout << " Fin curve A1" << *A1 << " " << icurve << " <=> " << icurveequi <<"-----" <<
			  NbCreatePointOnCurve << " == " <<i<<endl;
                      e->on  = ongequi;
                      e->v[0]=  A0;
                      e->v[1]=	A1;
                      e->ref = peequi->ref;
                      e->adj[0]=PreviousNewEdge;
                      e->adj[1]=0;
                      if (PreviousNewEdge)
			PreviousNewEdge->adj[1] =  e;
                      PreviousNewEdge = e;
		      //		      cout << "Last new edge " << nbe << " " << " on " << Gh.Number(pe->on) 
		      //   << " of curve =" <<pe->on->CurveNumber <<endl;
	    

#ifdef DRAWING1 
                      e->Draw();
                      A1->Draw();
                      A0->Draw();
		      //                      inquire();
#endif
                      assert(i==NbCreatePointOnCurve);
  
                    }
		   } //  end loop on equi curve 
                     
		  if (!phase)  { // 
		    Int4 NbSegOnCurve = Max((Int4)(L+0.5),(Int4) 1);// nb of seg
		    Lstep = L/NbSegOnCurve; 
		    Lcurve = L;
		    NbCreatePointOnCurve = NbSegOnCurve-1;
		    
		    for(Curve * curve= Gh.curves+icurve;curve;curve= curve->next)
		     {
		       NbOfNewEdge += NbSegOnCurve;
		       NbOfNewPoints += NbCreatePointOnCurve;
		     }
		    if(verbosity>5)
		      cout << icurve << " NbSegOnCurve = " <<  NbSegOnCurve << " Lstep=" 
			   << Lstep <<" " << NbOfNewPoints<< " NBPC= " << NbCreatePointOnCurve <<endl;
		    // do'nt 
		    //  if(NbCreatePointOnCurve<1) break;
		  }
		}//for(phase;;)
/*  modif FH add Curve class  		  
	    }}//for (iedge=0;iedge<BTh.nbe;iedge++) 
*/
// new code Curve class  	
     } //  end of curve loop 
// end new code	    
       // do the allocation
    if(step==0) 
      {
	//if(!NbOfNewPoints) break;// nothing ????? bug 
	if(nbv+NbOfNewPoints > nbvx) 
	  {
	    cerr << " Too much vertices on geometry " << nbv+NbOfNewPoints  << " >= " << nbvx << endl;
	    MeshError(3,this);
	  }
	//cout << " NbOfNewEdge" << NbOfNewEdge << " NbOfNewPoints " << NbOfNewPoints << endl;
	edges = new Edge[NbOfNewEdge];
	nbex = NbOfNewEdge;
	if(NbOfNewPoints) { // 
	   VerticesOnGeomEdge = new VertexOnGeom[NbOfNewPoints];
	   NbVertexOnBThEdge =NbOfNewPoints;
	   VertexOnBThEdge = new  VertexOnEdge[NbOfNewPoints];
	   NbVerticesOnGeomEdgex = NbOfNewPoints; }
	NbOfNewPoints =0;
	NbOfNewEdge = 0;
      }
    } // for(step;;)
  assert(nbe);

 delete [] bcurve;
 
  
#ifdef DRAWING1
  reffecran();
  InitDraw();
  Draw();
  inquire();
#endif
  
  Insert();
  ForceBoundary();
  FindSubDomain();
  
#ifdef DRAWING1
  reffecran();
  Draw();
  inquire();
#endif
  // NewPointsOld(*this) ;
//  BTh.ReMakeTriangleContainingTheVertex(); //  FH change => put in NewPoints
  //  for (Int4 iv=0;iv<BTh.nbv;iv++)
  //    BTh[iv].i = toI2(BTh[iv].r);
  NewPoints(BTh,KeepBackVertices) ;
  CurrentTh = 0;
//#undef  DRAWING1 
}

void Triangles::GeomToTriangles0(Int4 inbvx) 
{
  Gh.NbRef++;// add a ref to GH


  Int4 i,NbOfCurves=0,NbNewPoints,NbEdgeCurve;
  Real8 lcurve, lstep,s;
#ifdef DRAWING
 if (withrgraphique)
   {
  Gh.InitDraw() ;
  CurrentTh = this; }
#endif
  
  R2 AB;
  GeometricalVertex *a,*b;
  Vertex *va,*vb;
  GeometricalEdge * e;
  PreInit(inbvx);
  int  background = &BTh != this;
  //  int  SameGeom = background && (&BTh.Gh == &Gh);
  nbv = 0;
  NbVerticesOnGeomVertex=0;
  NbVerticesOnGeomEdge=0;
  for (i=0;i<Gh.nbv;i++)
    if (Gh[i].Required() && Gh[i].IsThe() ) NbVerticesOnGeomVertex++;
  VerticesOnGeomVertex = new VertexOnGeom[NbVerticesOnGeomVertex];  
//
  if( NbVerticesOnGeomVertex >= nbvx) 
    {
      cerr << " Too much vertices on geometry " << NbVerticesOnGeomVertex << " >= " << nbvx << endl;
      MeshError(1,this);
    }
  for (i=0;i<Gh.nbv;i++)
    if (Gh[i].Required()&& Gh[i].IsThe()  ) {//Gh  vertices Required
      if (nbv < nbvx)
	vertices[nbv] = Gh[i];
      Gh[i].to = vertices + nbv;// save Geom -> Th
      VerticesOnGeomVertex[nbv]= VertexOnGeom(*Gh[i].to,Gh[i]);
    //  cout << "--------- "  <<Number(Gh[i].to) << " " << Gh[i].to << " " << i << endl;
      nbv++;
    }
//  assert( Gh.nbv < nbvx);
  
  // Method in 2 step:  0 and 1 
  // 1) compute de nb of edge 
  // 2) construct the edge    
  // generation of the curves
  assert(! edges);
#ifdef DRAWING1
  reffecran();
#endif
// 2 step 
// --step=0 to compute the number of edges + alloc at end
// --step=1 to construct the edges
  for (int step=0;step<2;step++) 
    {//  for (int step=0;step<2;step++) 
      Int4 nbex = 0;
      nbe = 0;
      Int4 NbVerticesOnGeomEdge0=NbVerticesOnGeomEdge;
    //  cout <<  "  -------------- step =" << step << endl;
      Gh.UnMarkEdges();	
      NbOfCurves = 0;
      for (i=0;i<Gh.nbe;i++) {
	GeometricalEdge & ei = Gh.edges[i];   
	if (!ei.Dup()) // a good curve (not dup )
	for(int j=0;j<2;j++) 
	  if (!ei.Mark() && ei[j].Required()) { 
	    // warning ei.Mark() can be change in loop for(j=0;j<2;j++) 
	    //  cout << " New curve = " << NbOfCurves << endl;
	    Int4 nbvend  =0;
	  
           Edge * PreviousNewEdge=0;

          lstep = -1;//to do not create points
	  if(ei.Required())
	    {
	      if (j==0)
		if(step==0)
		  nbe++;
		else
		  { 
		    e = & ei;
		    a=ei(0)->The();
		    b=ei(1)->The();
		    assert(edges);
		    edges[nbe].v[0]=a->to;
		    edges[nbe].v[1]=b->to;;
		    edges[nbe].ref = e->ref;
		    edges[nbe].on = e;
		    edges[nbe].adj[0] = 0;
		    edges[nbe].adj[1] = 0;
#ifdef DRAWING1
		    edges[nbe].Draw();
#endif
		    nbe++;}
	    }
          else 
	    { // on curve ------
	      for ( int kstep=0;kstep<= step;kstep++)
		{ // begin  for ( int kstep=0;kstep<= step;kstep++)
		  // if 2nd step where 2 step
		  // -- 1 compute le length of the curve
		  // -- create the points and edge
		  PreviousNewEdge=0;
		  NbNewPoints=0;
		  NbEdgeCurve=0;
		  assert(nbvend < nbvx); 
		  lcurve =0;
		  s = lstep;
		  int k=j;
		  e = & ei;
		  a=ei(k)->The();
		  va = a->to;
		  e->SetMark();
		  //  cout << " New curve " ;
		  
		  // if SameGeo  We have go in the background geometry 
		  // to find the discretisation of the curve
		  
		  for(;;) 
		    { 
		      k = 1-k;
		      b= (*e)(k)->The();
		      AB = b->r - a->r;
		      Metric MA = background ? BTh.MetricAt(a->r) :a->m ;
		      Metric MB =  background ? BTh.MetricAt(b->r) :b->m ;
		      Real8 ledge = (MA(AB) + MB(AB))/2;
		      // 
		      const int MaxSubEdge = 10;
		      int NbSubEdge = 1;
		      Real8 lSubEdge[MaxSubEdge];
		      R2 A,B;
		      if (ledge < 1.5) 
			lSubEdge[0] = ledge;
		      else {
			NbSubEdge = Min( MaxSubEdge, (int) (ledge +0.5));
			A= a->r;
			Metric MAs =MA,MBs;
			// cout << " lSubEdge old=" << ledge 
			//      << " new " << A << MA << endl;
			ledge = 0; 
			Real8 x =0, xstep= 1. /  NbSubEdge;
			for (int kk=0; kk < NbSubEdge; kk++,A=B,MAs=MBs ) {
			  x += xstep;
			  B =  e->F(k ? x : 1-x);
			  MBs= background ? BTh.MetricAt(B) :Metric(1-x, MA, x ,MB);
			  AB = A-B;
			  lSubEdge[kk]= (ledge += (MAs(AB)+MBs(AB))/2);
			  // cout << "     " << lSubEdge[kk] << " x " << x  
			  //      << " " << A << B << MA << MB<< endl ;
			}
			//  cout << endl;
		      }
		      
		      Real8 lcurveb = lcurve+ ledge ;
		      while (lcurve<=s && s <= lcurveb && nbv < nbvend)
			{
			  // New points
			  
			  // Real8 aa=(lcurveb-s)/ledge;
			  // Real8 bb=(s-lcurve)/ledge;
			  
			  Real8 ss = s-lcurve;
			  // 1) find the SubEdge containing ss by dichotomie
			  int kk0=-1,kk1=NbSubEdge-1,kkk;
			  Real8 ll0=0,ll1=ledge,llk;
			  while (kk1-kk0>1)
			    {
			      if (ss < (llk=lSubEdge[kkk=(kk0+kk1)/2]))
				kk1=kkk,ll1=llk;
			      else
				kk0=kkk,ll0=llk;}
			  assert(kk1 != kk0);
			  
			  Real8 sbb = (ss-ll0  )/(ll1-ll0);
			  Real8 bb = (kk1+sbb)/NbSubEdge, aa=1-bb;
			  
			  // new vertex on edge
			  vb = &vertices[nbv++];
			  vb->m = Metric(aa,a->m,bb,b->m);
			  vb->ReferenceNumber = e->ref;
			  vb->DirOfSearch =NoDirOfSearch;
			  Real8 abcisse = k ? bb : aa;
			  vb->r =  e->F( abcisse );
			  VerticesOnGeomEdge[NbVerticesOnGeomEdge++]= VertexOnGeom(*vb,*e,abcisse);        
			  
			  // to take in account the sens of the edge
			  
			  s += lstep;
			  edges[nbe].v[0]=va;
			  edges[nbe].v[1]=vb;
			  edges[nbe].ref = e->ref;
			  edges[nbe].on = e;
			  edges[nbe].adj[0] = PreviousNewEdge;
			  if(PreviousNewEdge)
			    PreviousNewEdge->adj[1] = &edges[nbe];
#ifdef DRAWING1
			  vb->Draw();
			  edges[nbe].Draw();
#endif
			  PreviousNewEdge = edges + nbe;
			  nbe++;
#ifdef DEBUG1                 
			  cout << " new points " << nbv-1 << " " << vb->r ;
			  cout << " new edge " << nbe-1 << " " ;
			  cout << va << vb <<  " kk0 = " << kk0 
			       << " " << kk1 << " ss=" << ss ;
			  cout << " " << sbb << endl;
			  cout << "      " << aa << va->r << bb << vb->r 
			       <<" length=" << Norme(va->r-vb->r) << endl;
			  cout << "      s " << s << " lstep= " << lstep 
			       << " ledge= " << ledge 
			       << " lcurve= " << lcurve << endl;
#endif
			  va = vb;
			}
		      lcurve = lcurveb;
		      e->SetMark();
		      // cout << e-Gh.edges << ", " << k << " " 
		      //      <<(*e)[k].r <<" " <<(*e)[1-k].r <<" " 
		      //      << lcurve<< ";; " ;                          
		      a=b;
		      if (b->Required() ) break;
		      int kprev=k;
		      k = e->SensAdj[kprev];// next vertices
		      e = e->Adj[kprev];
		      assert(e);
		    }// for(;;)
		  vb = b->to;
		  //            cout << endl;
		  NbEdgeCurve = Max((Int4) (lcurve +0.5), (Int4) 1);
		  NbNewPoints = NbEdgeCurve-1;
		  if(!kstep)
		    { NbVerticesOnGeomEdge0 += NbNewPoints;
		    NbOfCurves++;}
		  
		  nbvend=nbv+NbNewPoints; 
		  
		  lstep = lcurve / NbEdgeCurve;
		  //   cout <<"lstep " << lstep << " lcurve " 
		  //    << lcurve << " NbEdgeCurve " << NbEdgeCurve << " " <<NbVerticesOnGeomEdge0<<" " <<NbVerticesOnGeomEdge<<" step =" <<step<<  endl;
		} 
	      // end of curve --
	      if (edges) { // last edges of the curves 
		edges[nbe].v[0]=va;
		edges[nbe].v[1]=vb;
		edges[nbe].ref = e->ref;
		edges[nbe].on = e;
		edges[nbe].adj[0] = PreviousNewEdge;
		edges[nbe].adj[1] = 0;
		if(PreviousNewEdge)
		  PreviousNewEdge->adj[1] = & edges[nbe];
		
		
#ifdef DRAWING1
		edges[nbe].Draw();
#endif
		nbe++;}
	      else
		nbe += NbEdgeCurve;
	    } // end on  curve ---
	} // if (edges[i][j].Corner())  
    } // for (i=0;i<nbe;i++)
    if(!step) {
     // cout << "edges " << edges << " VerticesOnGeomEdge " <<VerticesOnGeomEdge << endl;
      assert(!edges);
      assert(!VerticesOnGeomEdge);
      edges = new Edge[nbex=nbe];
      if(NbVerticesOnGeomEdge0)
	VerticesOnGeomEdge = new VertexOnGeom[NbVerticesOnGeomEdge0];
      assert(edges);
      assert(VerticesOnGeomEdge || NbVerticesOnGeomEdge0 ==0);
        // do the vertex on a geometrical vertex
       NbVerticesOnGeomEdge0 = NbVerticesOnGeomEdge;       
     }
     else 
       assert(NbVerticesOnGeomEdge == NbVerticesOnGeomEdge0);
    //     cout << " Nb of Curves = " << NbOfCurves << "nbe = " << nbe 
    //	  << "== " << nbex << "  nbv = " << nbv <<  endl;
    assert(nbex=nbe);
   } // for (step=0;step<2;step++)

#ifdef DRAWING1
  reffecran();
  InitDraw();
  Draw();
  inquire();
#endif
 
  Insert();
  ForceBoundary();
  FindSubDomain();

#ifdef DRAWING1
  reffecran();
  Draw();
  inquire();
#endif
   // NewPointsOld(*this) ;
    NewPoints(*this,0) ;
  CurrentTh = 0;
}

Edge** Triangles::MakeGeometricalEdgeToEdge()
 {
  assert(Gh.nbe);
  Edge **e= new (Edge* [Gh.nbe]);
  
  Int4 i;
  for ( i=0;i<Gh.nbe ; i++)
    e[i]=NULL;
  for ( i=0;i<nbe ; i++) 
    { 
      Edge * ei = edges+i;
      GeometricalEdge *on = ei->on; 
      e[Gh.Number(on)] = ei;    
    }
  for ( i=0;i<nbe ; i++) 
    for (int ii=0;ii<2;ii++) { 
       Edge * ei = edges+i;
       GeometricalEdge *on = ei->on;
      int j= ii;
      while (!(*on)[j].Required()) { 
	//	cout << i << " " << ii << " j= " << j << " curve = " 
	//           <<  on->CurveNumber << "  " << Gh.Number(on) << " on " << j 
	//   << " s0 " << Gh.Number( (*on)[0]) << " s1  " << Gh.Number( (*on)[1]) 
	//   << " ->  " ;
       Adj(on,j); // next geom edge
       j=1-j;
       //       cout << Gh.Number(on) << "  " << j  << " e[ON] =  " <<  e[Gh.Number(on)] 
       //    << " s0 " << Gh.Number( (*on)[0]) << " s1  " << Gh.Number( (*on)[1]) << endl; 
       if (e[Gh.Number(on)])  break; // optimisation     
       e[Gh.Number(on)] = ei; 
       }
   }

  int kk=0;
     for ( i=0;i<Gh.nbe ; i++)
       if (!e[i]) 
	 if(kk++<10) {
	   cerr << " Bug -- the geometrical edge " << i << " is on no edge curve = " << Gh.edges[i].CurveNumber 
		<< " s0 " << Gh.Number( Gh.edges[i][0]) << " s1  " << Gh.Number( Gh.edges[i][1]) << endl; 
	 //	 assert( e[i]);
       }
  if(kk) MeshError(997,this);

  return e;
 }

Triangles::~Triangles() 
{
  assert(NbRef<=0);
  if (CurrentTh == this) CurrentTh=0;
  if(verbosity>10)
    cout << " ~Triangles "<< this  <<" "<< identity << endl;
  if(vertices)  delete [] vertices;
  if(edges)     delete [] edges;
  if(triangles) delete [] triangles;
  if(quadtree)  delete  quadtree;
  if(ordre)     delete [] ordre;
  if( subdomains) delete []  subdomains;
  if (VerticesOnGeomEdge) delete [] VerticesOnGeomEdge;
  if (VerticesOnGeomVertex) delete [] VerticesOnGeomVertex;
  if (name) delete [] name;
  if (identity) delete [] identity;
  if (VertexOnBThVertex) delete [] VertexOnBThVertex;
  if (VertexOnBThEdge) delete [] VertexOnBThEdge;
  
  if (&Gh) 
    {
      if (Gh.NbRef>0) Gh.NbRef--;
      else if (Gh.NbRef==0) delete &Gh;
    }
  if (&BTh && (&BTh != this))
    {
      if (BTh.NbRef>0) BTh.NbRef--;
      else if (BTh.NbRef==0) delete &BTh;
    }
  PreInit(0); // set all to zero 
  
}

void Triangles::SetIntCoor(const char * strfrom)
{
    pmin =  vertices[0].r;
    pmax =  vertices[0].r;

    // recherche des extrema des vertices pmin,pmax
    Int4 i;
    for (i=0;i<nbv;i++) {
      pmin.x = Min(pmin.x,vertices[i].r.x);
      pmin.y = Min(pmin.y,vertices[i].r.y);
      pmax.x = Max(pmax.x,vertices[i].r.x);
      pmax.y = Max(pmax.y,vertices[i].r.y);
    }
    R2 DD = (pmax-pmin)*0.05;
    pmin = pmin-DD;
    pmax = pmax+DD; 
    coefIcoor= (MaxICoor)/(Max(pmax.x-pmin.x,pmax.y-pmin.y));
    assert(coefIcoor >0);

    // generation of integer coord  
    for (i=0;i<nbv;i++) {
      vertices[i].i = toI2(vertices[i].r);    
    }
#ifdef DRAWING
    xGrafCoef = coefIcoor;
    yGrafCoef = coefIcoor;
    xGrafOffSet = pmin.x;
    yGrafOffSet = pmin.y;
#ifdef DRAWING1
    rattente(1);
#endif
#endif

    // computation of the det 
    int Nberr=0;
    for (i=0;i<nbt;i++)
      {
	Vertex & v0 = triangles[i][0];
	Vertex & v1 = triangles[i][1];
	Vertex & v2 = triangles[i][2];
      if ( &v0 && &v1 &&  &v2 ) // a good triangles;
      {
	triangles[i].det= det(v0,v1,v2);
	if (triangles[i].det <=0 && Nberr++ <10)
	  {
	    if(Nberr==1)
	      if (strfrom)
		cerr << "+++ Fatal Error " << strfrom << "(SetInCoor)  Error :  area of Triangle < 0 " << endl; 
	      else 
		cerr << "+++  Fatal Error Triangle (in SetInCoor) area of Triangle < 0" << endl;
	    cerr << " Triangle " << i << "  det  (I2) = " << triangles[i].det ;
	    cerr << " (R2) " << Det(v1.r-v0.r,v2.r-v0.r);
	    cerr << "; The 3  vertices " << endl;
	    cerr << Number(v0) << " "  << Number(v1) << " " 
		 << Number(v2) << " : " ;
	    cerr << v0.r << v1.r << v2.r << " ; ";
	    cerr << v0.i << v1.i << v2.i << endl;
	  }
      }
    else
      triangles[i].det= -1; // Null triangle; 
      }
    if (Nberr) MeshError(899,this);
	
}

void Triangles::FillHoleInMesh() 
{
  Triangles * OldCurrentTh =CurrentTh;
  CurrentTh=this;
  //  Int4 NbTold = nbt;
  // generation of the integer coor
  {
 
    //  double coef = coefIcoor;
    // recherche des extrema des vertices pmin,pmax
    Int4 i;
    if(verbosity>2)
      cout << "  -- FillHoleInMesh: Nb of vertices =" << nbv 
	   << " Pmin = "<< pmin << " Pmax = "<< pmax << endl;
    
    assert(ordre);
    for (i=0;i<nbv;i++) 
      ordre[i]= 0 ;
    

    NbSubDomains =0;
    
  // generation of the adjacence of the triangles
    SetOfEdges4 * edge4= new SetOfEdges4(nbt*3,nbv);
    Int4 * st = new Int4[nbt*3];
    for (i=0;i<nbt*3;i++)
      st[i]=-1;
    Int4 kk =0;
    for (i=0;i<nbe;i++)
      kk += (i == edge4->addtrie(Number(edges[i][0]),Number(edges[i][1])));
    if (kk != nbe)
      { 
	cerr << " Some Double edge in the mesh, the number is " << kk-nbe << endl;
	MeshError(1002,this);
      }
    for (i=0;i<nbt;i++)
      for (int j=0;j<3;j++)
	{
	  // Int4 i0,i1;
	  Int4 k =edge4->addtrie(Number(triangles[i][VerticesOfTriangularEdge[j][0]]),
				 Number(triangles[i][VerticesOfTriangularEdge[j][1]]));
	  Int4 invisible = triangles[i].Hidden(j);
	  if(st[k]==-1)
	    st[k]=3*i+j;
	  else if(st[k]>=0) {
	    assert( ! triangles[i].TriangleAdj(j) && !triangles[st[k] / 3].TriangleAdj((int) (st[k]%3)));
	    
	    triangles[i].SetAdj2(j,triangles + st[k] / 3,(int) (st[k]%3));
	    if (invisible)  triangles[i].SetHidden(j);
	    if (k<nbe) {
	      triangles[i].SetLocked(j);
	    }
	    st[k]=-2-st[k]; }
	  else {
	    cerr << " The edge (" 
	     << Number(triangles[i][VerticesOfTriangularEdge[j][0]])
		 << " , " 
		 << Number(triangles[i][VerticesOfTriangularEdge[j][1]])
		 << " ) is in more than 2 triangles " <<k <<endl;
	    cerr << " Edge " << j << " Of Triangle " << i << endl;
	    cerr << " Edge " << (-st[k]+2)%3 << " Of Triangle " << (-st[k]+2)/3  << endl;
	    cerr << " Edge " << triangles[(-st[k]+2)/3].NuEdgeTriangleAdj((int)((-st[k]+2)%3))
		 << " Of Triangle " <<  Number(triangles[(-st[k]+2)/3].TriangleAdj((int)((-st[k]+2)%3))) << endl;
	    MeshError(9999,this);}	


	}
    if(verbosity>5) {
      cout << "    On Mesh " << name << endl;
      cout << "    - The number of Vertices  = " << nbv << endl;
      cout << "    - The number of Triangles = " << nbt << endl;
      cout << "    - The number of given edge = " << nbe << endl;
      cout << "    - The number of all edges = " << edge4->nb() << endl;
      cout << "    - The Euler number = 1-Nb Of Hole = " << nbt-edge4->nb()+nbv << endl; }
    
    
    // check the consistant of edge[].adj and the geometrical required  vertex
     Int4 k=0;
     for (i=0;i<edge4->nb();i++)
       if (st[i] >=0) // edge alone 
	 if (i < nbe) 
	   {
	     Int4 i0=edge4->i(i);ordre[i0] = vertices+i0;
	     Int4 i1=edge4->j(i);ordre[i1] = vertices+i1;
	   }
     	 else {
	   k++;
	   if (verbosity>20 && k <20) 
	     {
	     Int4 i0=edge4->i(i);
	     Int4 i1=edge4->j(i);
	     cerr << " Lose boundary edges " << i << " : " << i0 << " " << i1 << endl;
	     }
	 }
     	 
      if(k != 0) {
	if (verbosity>20)
	  {
	    cout << " The given edge are " << endl;
	      for (int i=0;i< nbe;i++)
		cout <<  " Edge " << i << " : " <<  Number(edges[i][0]) << " " <<  Number(edges[i][1]) 
		     << " " << edges[i].ref << endl; 
	  }
	cerr << k << " boundary edges  are not defined as edges " << endl;
	MeshError(9998,this);
      }
      // generation of the mesh with boundary points   
      Int4 nbvb = 0;
      for (i=0;i<nbv;i++)
	{ 
	  vertices[i].t=0;
	  vertices[i].vint=0;
	if (ordre[i]) 
	  ordre[nbvb++] = ordre[i];
	}

      Triangle *savetriangles= triangles;
      Int4 savenbt=nbt;
      Int4 savenbtx=nbtx;
      SubDomain * savesubdomains = subdomains;
      subdomains = 0;

      Int4  Nbtriafillhole = 2*nbvb;
      Triangle * triafillhole =new Triangle[Nbtriafillhole];
      if (verbosity>9)
	cout << " Nbtriafillhole triafillhole*" << triafillhole << endl; 
      triangles =  triafillhole;

      nbt=2;
      nbtx= Nbtriafillhole;
      
	for (i=2 ; det( ordre[0]->i, ordre[1]->i, ordre[i]->i ) == 0;) 
	  if  ( ++i >= nbvb) {
	    cerr << "FillHoleInMesh: All the vertices are aline " << nbvb << endl;
	    MeshError(998,this); }
	Exchange( ordre[2], ordre[i]);

	Vertex *  v0=ordre[0], *v1=ordre[1];


	triangles[0](0) = 0; // sommet pour infini 
	triangles[0](1) = v0;
	triangles[0](2) = v1;
	
	triangles[1](0) = 0;// sommet pour infini 
	triangles[1](2) = v0;
	triangles[1](1) = v1;
	const int e0 = OppositeEdge[0];
	const int e1 = NextEdge[e0];
	const int e2 = PreviousEdge[e0];
	triangles[0].SetAdj2(e0, &triangles[1] ,e0);
	triangles[0].SetAdj2(e1, &triangles[1] ,e2);
	triangles[0].SetAdj2(e2, &triangles[1] ,e1);
	
	triangles[0].det = -1;  // faux triangles
	triangles[1].det = -1;  // faux triangles
	
	triangles[0].SetTriangleContainingTheVertex();
	triangles[1].SetTriangleContainingTheVertex();
	
	triangles[0].link=&triangles[1];
	triangles[1].link=&triangles[0];
	
#ifdef DEBUG 
	triangles[0].check();
	triangles[1].check();
#endif  
	//  nbtf = 2;
	if (  !quadtree ) 
	   delete  quadtree; // ->ReInitialise();
	
	  quadtree = new QuadTree(this,0);
	quadtree->Add(*v0);
	quadtree->Add(*v1);
	
	// on ajoute les sommets un a un 
	Int4 NbSwap=0;
	for (Int4 icount=2; icount<nbvb; icount++) {
	  
	  Vertex *vi  = ordre[icount];
	  //	  cout << " Add vertex " <<  Number(vi) << endl;
	  Icoor2 dete[3];
	  Triangle *tcvi = FindTriangleContening(vi->i,dete);
	  quadtree->Add(*vi); 
	  Add(*vi,tcvi,dete);
	  NbSwap += vi->Optim(1,1);
	  
#ifdef DRAWING2
	  cout << Number(vi) << " " <<  NbSwap <<  endl;
	reffecran();
	Draw();
	vi->Draw();
	inquire();
#endif
	}// end loop on  icount	
#ifdef DRAWING1
	inquire();
#endif
	
	//Int4 nbtfillhole = nbt;
	 // inforce the boundary 
	 TriangleAdjacent ta(0,0);
	 Int4 nbloss = 0,knbe=0;
	 for ( i = 0; i < nbe; i++) 
	  if (st[i] >=0)  // edge alone => on border ...  FH oct 2009
	   {
	     Vertex & a=edges[i][0], & b =    edges[i][1];
	     if (a.t && b.t) // le bug est la si maillage avec des bod non raffine 1.
	       {
		 knbe++;
		 if (ForceEdge(a,b,ta)<0)
		   nbloss++;
	       }
	   }
	 if(nbloss)
	   {
	     cerr << " we loss some  " << nbloss << " "  << " edges other " << knbe << endl;
	     MeshError(1100,this);
	   }
	 FindSubDomain(1);
         // remove all the hole 
	 // remove all the good sub domain
	 Int4 krm =0;
	 for (i=0;i<nbt;i++)
	   if (triangles[i].link) // remove triangles
	     {
	       krm++;
	       for (int j=0;j<3;j++)
		 {
		   TriangleAdjacent ta =  triangles[i].Adj(j);
		   Triangle & tta = * (Triangle *) ta;
		   if(! tta.link) // edge between remove and not remove 
		     { // change the link of ta;
		       int ja = ta;
		       Vertex *v0= ta.EdgeVertex(0);
		       Vertex *v1= ta.EdgeVertex(1);
		       Int4 k =edge4->addtrie(v0?Number(v0):nbv,v1? Number(v1):nbv);
		       assert(st[k] >=0); 
		       tta.SetAdj2(ja,savetriangles + st[k] / 3,(int) (st[k]%3));
		       ta.SetLock();
		       st[k]=-2-st[k]; 
		     }
		 }
	     }
	 Int4 NbTfillHoll =0;
	 for (i=0;i<nbt;i++)
	   if (triangles[i].link) {
	     triangles[i]=Triangle((Vertex *) NULL,(Vertex *) NULL,(Vertex *) NULL);
	     triangles[i].color=-1;
	   }
	   else
	     {
	      triangles[i].color= savenbt+ NbTfillHoll++;
#ifdef DEBUG 
	     triangles[i].check();
#endif
      }
	 // cout <<      savenbt+NbTfillHoll << " " <<  savenbtx  << endl;
     assert(savenbt+NbTfillHoll <= savenbtx );
     // copy of the outside triangles in saveTriangles 
     for (i=0;i<nbt;i++)
       if(triangles[i].color>=0) 
	 {
          savetriangles[savenbt]=triangles[i];
	  savetriangles[savenbt].link=0;
	  savenbt++;
	 }
     // gestion of the adj
      k =0;
      Triangle * tmax = triangles + nbt;
      for (i=0;i<savenbt;i++)  
	{ 
	  Triangle & ti = savetriangles[i];
	  for (int j=0;j<3;j++)
	    {
	      Triangle * ta = ti.TriangleAdj(j);
	      int aa = ti.NuEdgeTriangleAdj(j);
	      int lck = ti.Locked(j);
	      if (!ta) k++; // bug 
	      else if ( ta >= triangles && ta < tmax) 
		{
		  ta= savetriangles + ta->color;
		  ti.SetAdj2(j,ta,aa);
		  if(lck) ti.SetLocked(j);
		}
	    }
	}
     //	 OutSidesTriangles = triangles;
      //	Int4 NbOutSidesTriangles = nbt;
	 
	 // restore triangles;
	 nbt=savenbt;
	 nbtx=savenbtx;
	 delete [] triangles;
	 delete [] subdomains;
	 triangles = savetriangles;
	 subdomains = savesubdomains;
	 //	 cout <<  triangles << " <> " << OutSidesTriangles << endl; 
	     /*	 k=0;
	 for (i=0;i<nbt;i++)
	   for (int j=0;j<3;j++)
	     if (!triangles[i].TriangleAdj(j))
	     k++;
	     */
	 if (k) {
	   cerr << "Error Nb of triangles edge alone = " << k << endl;
	   MeshError(9997,this);
	 }
	 FindSubDomain();
	 // cout << " NbTOld = " << NbTold << " ==  " << nbt - NbOutT << " " << nbt << endl;
  
    // 

    delete edge4;
    delete [] st;
    for (i=0;i<nbv;i++)
      quadtree->Add(vertices[i]);
    
    SetVertexFieldOn();

    for (i=0;i<nbe;i++)
      if(edges[i].on) 
	for(int j=0;j<2;j++)
	  if (!edges[i].adj[j])
	    if(!edges[i][j].on->IsRequiredVertex()) {
	      cerr << " Erreur adj et sommet requis edges [" << i <<  "][ " << j << "]= "
		   <<  Number(edges[i][j]) << " : "  << " on = " << Gh.Number(edges[i].on) ;
	      if (edges[i][j].on->OnGeomVertex())
		cerr << " vertex " << Gh.Number(edges[i][j].on->gv);
	      else if (edges[i][j].on->OnGeomEdge())
		cerr << "Edges " << Gh.Number(edges[i][j].on->ge);
	      else
		cerr << " = " << edges[i][j].on ;
	      cerr << endl;
	  }
    
#ifdef DRAWING1
    InitDraw();
#endif
    
  }
  CurrentTh=OldCurrentTh;
}

Triangles::Triangles(Triangles & Th,Geometry * pGh,Triangles * pBth,Int4 nbvxx) // COPY OPERATOR
: Gh(*(pGh?pGh:&Th.Gh)), BTh(*(pBth?pBth:this))
{
  Gh.NbRef++;
  nbvxx = Max(nbvxx,Th.nbv); 
  Int4 i;
  // do all the allocation to be sure all the pointer existe
  
  char * cname = 0;
  if (Th.name) 
    {
      cname = new char[strlen(Th.name)+1];
      strcpy(cname,Th.name);
    }
  PreInit(nbvxx,cname);// to make the allocation 
  // copy of triangles
  nt=Th.nt;
  nbv = Th.nbv;
  nbt = Th.nbt;
  nbiv = Th.nbiv;
  nbe = Th.nbe;
  NbSubDomains = Th.NbSubDomains;
  NbOutT = Th.NbOutT;
  NbOfQuad =  Th.NbOfQuad ;
  NbOfSwapTriangle =0;
  NbVerticesOnGeomVertex = Th.NbVerticesOnGeomVertex;
  if(NbVerticesOnGeomVertex)
    VerticesOnGeomVertex = new VertexOnGeom[NbVerticesOnGeomVertex];
  NbVerticesOnGeomEdge = Th.NbVerticesOnGeomEdge;
  if (NbVerticesOnGeomEdge)
     VerticesOnGeomEdge = new VertexOnGeom[NbVerticesOnGeomEdge] ;
  if (& BTh == & Th.BTh) // same back ground 
    {
      BTh.NbRef++;
      NbVertexOnBThVertex = Th.NbVertexOnBThVertex;
      if(NbVertexOnBThVertex)
	VertexOnBThVertex = new VertexOnVertex[NbVertexOnBThVertex];
      NbVertexOnBThEdge = Th.NbVertexOnBThEdge;
      if(NbVertexOnBThEdge)
	VertexOnBThEdge = new VertexOnEdge[NbVertexOnBThEdge];
    }
   else 
     { // no add on back ground mesh 
       BTh.NbRef++;
       NbVertexOnBThVertex=0;
       VertexOnBThVertex=0;
       NbVertexOnBThEdge=0;
       VertexOnBThEdge=0;
       //       assert (& BTh == this); // --- a voir 
		     
    }


  if(nbe)
    edges = new Edge[nbe];
  if(NbSubDomains)
    subdomains = new SubDomain[NbSubDomains];
  pmin = Th.pmin;
  pmax = Th.pmax;
  coefIcoor = Th.coefIcoor;
  for(i=0;i<nbt;i++)
     triangles[i].Set(Th.triangles[i],Th,*this);
  for(i=0;i<nbe;i++)
     edges[i].Set(Th,i,*this);
  for(i=0;i<nbv;i++)
     vertices[i].Set(Th.vertices[i],Th,*this);
  for(i=0;i<NbSubDomains;i++)  
    subdomains[i].Set(Th,i,*this);
  for (i=0;i<NbVerticesOnGeomVertex;i++)
    VerticesOnGeomVertex[i].Set(Th.VerticesOnGeomVertex[i],Th,*this);
  for (i=0;i<NbVerticesOnGeomEdge;i++)
    VerticesOnGeomEdge[i].Set(Th.VerticesOnGeomEdge[i],Th,*this);
  quadtree=0;


  //  assert(!OutSidesTriangles);
}

/** -- old with a bug we loss some time last swap

Int4  Triangle::Optim(Int2 i,int koption)
{
  // turn in the positif sens around vertex s  
  register  Int4 NbSwap =0;
  register Vertex * s  = ns[i];
  register Triangle * tbegin=0 , *t = this , *ttc;
  register int k=0,j = EdgesVertexTriangle[i][0],jc;
  tbegin=t;
  do {
    k++; 
#ifdef DEBUG
    assert( s == & (*t)[VerticesOfTriangularEdge[j][1]] );
#endif
#ifdef DRAWING1 
    t->Draw();
    DrawMark( s->r);
#endif
     ttc =  t->at[j];
     jc = NextEdge[t->aa[j]&3];
     cout << *t << " " <<  VerticesOfTriangularEdge[j][1] << "\n\t try swap " << * ttc << " " << jc ;
    while ( ttc->swap(jc,koption)) {
      NbSwap++,assert(k++<20000);
      ttc =  t->at[j];
      jc = NextEdge[t->aa[j]&3];
      cout << "\n\t s  " <<  *ttc << " " << jc << endl;
    }
    cout << endl;
    t = ttc;
    j = NextEdge[jc];
    assert(k<20000);
  } while ( (tbegin != t)); 
  
  return NbSwap;
}
*/
Int4  Triangle::Optim(Int2 i,int koption)
{
  // turne around in positif sens
  Int4 NbSwap =0;
#ifdef DEBUG
   Vertex * s  = ns[i];
#endif
   Triangle  *t = this;
   int k=0,j =OppositeEdge[i];
   int jp = PreviousEdge[j];
   // initialise   tp, jp the previous triangle & edge
   Triangle *tp= at[jp];
   jp = aa[jp]&3;
#ifdef DEBUG
   assert(tp->at[jp] == this);
#endif
   do {
#ifdef DEBUG
    assert(k++<20000);
    assert( s == & (*t)[OppositeVertex[j]] );
#endif
    //    cout << *t << " " <<  j  << "\n\t try swap " ;
    while (t->swap(j,koption))
      {
	NbSwap++;
	assert(k++<20000);
	t=  tp->at[jp];      // set unchange t qnd j for previous triangles
	j=  NextEdge[tp->aa[jp]&3];
	//   cout << "\n\t s  " <<  *t << " " << j << endl;
#ifdef DEBUG
	assert( s == & (*t)[OppositeVertex[j]] );
#endif
      }
    // end on this  Triangle 
    tp = t;
    jp = NextEdge[j];

    t=  tp->at[jp];      // set unchange t qnd j for previous triangles
    j=  NextEdge[tp->aa[jp]&3];

   } while( t != this);
   return NbSwap;
}

 void Triangles::SmoothingVertex(int nbiter,Real8 omega )
  { 
  //  if quatree exist remove it end reconstruct
    if (quadtree) delete quadtree;
    quadtree=0;
    ReMakeTriangleContainingTheVertex();
    Triangle vide; // a triangle to mark the boundary vertex
    Triangle   ** tstart= new Triangle* [nbv];
    Int4 i,j,k;
    //   attention si Background == Triangle alors on ne peut pas utiliser la rechech rapide 
    if ( this == & BTh)
     for ( i=0;i<nbv;i++)
      tstart[i]=vertices[i].t;     
    else 
     for ( i=0;i<nbv;i++)
      tstart[i]=0;
    for ( j=0;j<NbVerticesOnGeomVertex;j++ ) 
      tstart[ Number(VerticesOnGeomVertex[j].mv)]=&vide;
    for ( k=0;k<NbVerticesOnGeomEdge;k++ ) 
      tstart[ Number(VerticesOnGeomEdge[k].mv)]=&vide;
    if(verbosity>2) 
      cout << "  -- SmoothingVertex: nb Iteration = " << nbiter << " Omega = " << omega << endl;
    for (k=0;k<nbiter;k++)
     {
      Int4 i,NbSwap =0;
      Real8 delta =0;
      for ( i=0;i<nbv;i++)
        if (tstart[i] != &vide) // not a boundary vertex 
	  delta=Max(delta,vertices[i].Smoothing(*this,BTh,tstart[i],omega));
      if (!NbOfQuad)
      for ( i=0;i<nbv;i++)
        if (tstart[i] != &vide) // not a boundary vertex 
	  NbSwap += vertices[i].Optim(1);
       if (verbosity>3)
	 cout << "    Move max = " <<  sqrt(delta) << " iteration = " 
	      << k << " Nb of Swap = " << NbSwap << endl;
       }

    delete [] tstart;
    if (quadtree) quadtree= new QuadTree(this);
  }
void Triangles::MakeQuadTree()
{  
  if(verbosity>8)
    cout << "      MakeQuadTree" << endl;
  if (  !quadtree )  quadtree = new QuadTree(this);

  
#ifdef DRAWING1
  quadtree->Draw();
  rattente(1);
  reffecran();
  quadtree->Draw();
  rattente(1);
#endif
  
}
void  Triangles::ShowRegulaty() const// Add FH avril 2007 
{
   const  Real8  sqrt32=sqrt(3.)*0.5; 
   const Real8  aireKh=sqrt32*0.5;
   D2  Beq(1,0),Heq(0.5,sqrt32);
   D2xD2 Br(D2xD2(Beq,Heq).t());
   D2xD2 B1r(Br.inv());
   /*   D2xD2 BB = Br.t()*Br;
    cout << " BB = " << BB << " " << Br*B1r <<  endl; 
    MetricAnIso MMM(BB.x.x,BB.x.y,BB.y.y);
    MatVVP2x2 VMM(MMM);
    cout << " " << VMM.lambda1 << " " << VMM.lambda2 <<  endl; 
   */
    double gammamn=1e100,hmin=1e100;
    double gammamx=0,hmax=0;
    double beta=1e100;
    double beta0=0;
    double  alpha2=0;
    double area=0,Marea=0;
   // Real8 cf= Real8(coefIcoor);
   // Real8 cf2= 6.*cf*cf;
    int nt=0;
    for (int it=0;it<nbt;it++)
      if ( triangles[it].link) 
	{
	  nt++;
	  Triangle &K=triangles[it];
	  Real8  area3= Area2((R2) K[0],(R2) K[1],(R2) K[2])/6.;
	  area+= area3;
	  D2xD2 B_Kt(K[0],K[1],K[2]);
	  D2xD2 B_K(B_Kt.t());
	  D2xD2 B1K = Br*B_K.inv();
	  D2xD2 BK =  B_K*B1r;
	  D2xD2 B1B1 = B1K.t()*B1K;
	  MetricAnIso MK(B1B1.x.x,B1B1.x.y,B1B1.y.y);
	  MatVVP2x2 VMK(MK);
	  alpha2 = Max(alpha2,Max(VMK.lambda1/VMK.lambda2,VMK.lambda2/VMK.lambda1));
	  // cout << B_K << " * " << B1r << " == " << BK << " " << B_K*B_K.inv() << endl;
	  Real8 betaK=0;
	  
	  for (int j=0;j<3;j++)
	    {
	      Real8 he= Norme2(R2(K[j],K[(j+1)%3]));
	      hmin=Min(hmin,he);
	      hmax=Max(hmax,he);
	      Vertex & v=K[j];
	      D2xD2 M((MetricAnIso)v);
	      betaK += sqrt(M.det());
	      
	      D2xD2 BMB = BK.t()*M*BK;
	      MetricAnIso M1(BMB.x.x,BMB.x.y,BMB.y.y);
	      MatVVP2x2 VM1(M1);
	      //cout << B_K <<" " <<  M << " " <<  he << " " << BMB << " " << VM1.lambda1 << " " << VM1.lambda2<<   endl; 
	      gammamn=Min3(gammamn,VM1.lambda1,VM1.lambda2);
	      gammamx=Max3(gammamx,VM1.lambda1,VM1.lambda2);		
	    }
	  betaK *= area3;//  1/2 (somme sqrt(det))* area(K)
	  Marea+= betaK;
	  // cout << betaK << " " << area3 << " " << beta << " " << beta0 << " " << area3*3*3*3 <<endl;
	  beta=min(beta,betaK);
	  beta0=max(beta0,betaK);
	}   
    area*=3; 
    gammamn=sqrt(gammamn);
    gammamx=sqrt(gammamx);    
    cout << "  -- adaptmesh Regulary:  Nb triangles " << nt <<  " , h  min " << hmin  << " , h max " << hmax << endl;  
    cout << "     area =  " << area << " , M area = " << Marea << " , M area/( |Khat| nt) " << Marea/(aireKh*nt) << endl; 
    cout << "     infiny-regulaty:  min " << gammamn << "  max " << gammamx << endl;  
    cout << "     anisomax  "<< sqrt(alpha2) << ", beta max = " << 1./sqrt(beta/aireKh) 
	 << " min  "<<  1./sqrt(beta0/aireKh)  << endl;
}
void  Triangles::ShowHistogram() const
 {

    const Int4 kmax=10;
    const Real8 llmin = 0.5,llmax=2;
     const Real8 lmin=log(llmin),lmax=log(llmax),delta= kmax/(lmax-lmin);
    Int4 histo[kmax+1];
    Int4 i,it,k, nbedges =0;
    for (i=0;i<=kmax;i++) histo[i]=0;
    for (it=0;it<nbt;it++)
	if ( triangles[it].link) 
	{
	     
	    for (int j=0;j<3;j++)
	    {
		Triangle *ta = triangles[it].TriangleAdj(j);
		if ( !ta || !ta->link || Number(ta) >= it) 
		{ 
		    Vertex & vP = triangles[it][VerticesOfTriangularEdge[j][0]];
		    Vertex & vQ = triangles[it][VerticesOfTriangularEdge[j][1]];
		    if ( !& vP || !&vQ) continue;
		    R2 PQ = vQ.r - vP.r;
		    Real8 l = log(LengthInterpole(vP,vQ,PQ));
#ifdef DRAWING2             
		    if (l>1.4)  {
			penthickness(3);
			vP.MoveTo(),vQ.LineTo();
			penthickness(1);
			cout << "   l = " << l << Number(vP) << " edge = " << Number(vQ) << endl;
		    }
#endif             
		    nbedges++;
		    k = (int) ((l - lmin)*delta);
		    k = Min(Max(k,0L),kmax);
		    histo[k]++;
		}
	    }
	}  
	    cout << "  -- Histogram of the unit mesh,  nb of edges" << nbedges << endl <<endl;
 
    cout << "        length of edge in   | % of edge  | Nb of edges " << endl;
    cout << "        ------------------- | ---------- | ----------- " << endl;
    for   (i=0;i<=kmax;i++)
      { 
       cout << "    " ;
        cout.width(10);
        if (i==0) cout  << " 0 " ;
        else cout  << exp(lmin+i/delta) ;
        cout.width(); cout << "," ;
        cout.width(10);
        if (i==kmax) cout << " +infty " ;
        else cout  << exp(lmin+(i+1)/delta) ;
        cout.width();cout << "   |   " ;
        
       cout.precision(4);
       cout.width(6);
       cout <<  ((long)  ((10000.0 * histo[i])/ nbedges))/100.0 ;
       cout.width();
       cout.precision();
       cout <<  "   |   " << histo[i] <<endl;
      }
    cout << "        ------------------- | ---------- | ----------- " << endl <<endl;
    
 }

int  Triangles::Crack()
  { 
    assert(NbCrackedEdges ==0 || NbCrackedVertices >0); 
    for (int i=0;i<NbCrackedEdges;i++)
       CrackedEdges[i].Crack();
    return NbCrackedEdges;
  }
  
int Triangles::UnCrack() 
{ 
  assert(NbCrackedEdges ==0 || NbCrackedVertices >0); 
  for (int i=0;i<NbCrackedEdges;i++)
    CrackedEdges[i].UnCrack();
  return NbCrackedEdges;
}

int Triangles::CrackMesh()
{
    Triangles *CurrentThOld = CurrentTh;
  //  computed the number of cracked edge
  int i,k;
  for (k=i=0;i<nbe;i++)
    if(edges[i].on->Cracked()) k++;
  if( k==0) return 0;
    CurrentTh = this;
  cout << " Nb of Cracked Edges = " << k << endl;
  NbCrackedEdges =k;
  CrackedEdges = new  CrackedEdge[k];
  //  new edge
  Edge * e = new Edge[ nbe + k];

  // copy
  for (i=0;i<nbe;i++) 
    e[i] = edges[i];
  delete edges;
  edges = e;

  const int  nbe0  = nbe;
  for (k=i=0;i<nbe0;i++) // on double les arete cracked 
    if(edges[i].on->Cracked())
      {
	e[nbe] = e[i];
	//  return the edge 
	e[nbe].v[0] =  e[i].v[1];
	e[nbe].v[1] =  e[i].v[0];
	e[nbe].on = e[i].on->link ; // fqux 
	CrackedEdges[k++]=CrackedEdge(edges,i,nbe);
	nbe++;
      }
  ReMakeTriangleContainingTheVertex() ; 
  //  
  int nbcrakev  =0;
  Vertex *vlast = vertices + nbv;
  Vertex *vend = vertices + nbvx; // end of array
  for (int iv=0;iv<nbv;iv++) // vertex 
    {
      Vertex & v= vertices[iv];
      Vertex * vv = & v;  
      int kk=0; // nb cracked
      int kc=0; 
      int kkk =0; // nb triangle  with same number 
      Triangle * tbegin = v.t;
      int i  = v.vint;       
      assert(tbegin && (i >= 0 ) && (i <3));
      // turn around the vertex v
      TriangleAdjacent ta(tbegin,EdgesVertexTriangle[i][0]);// previous edge
      int k=0;
      do {
	int kv = VerticesOfTriangularEdge[ta][1];
	k++; 
	Triangle * tt (ta);
	if ( ta.Cracked() ) 
	  {   
	    TriangleAdjacent tta=(ta.Adj());
	    assert(tta.Cracked());
	    if ( kk == 0) tbegin=ta,kkk=0;  //  begin by a cracked edge  => restart                
	    if (  kkk ) { kc =1;vv = vlast++;  kkk = 0; } // new vertex if use 
	    kk++;// number of cracked edge view                 
	  }
	if ( tt->link ) { // if good triangles store the value 
	  int it = Number(tt);
	  assert(it < nt);
	  (*tt)(kv)= vv; //   Change the vertex of triangle 
	  if(vv<vend) {*vv= v;vv->ReferenceNumber=iv;} // copy the vertex value + store the old vertex number in ref 
	  //	  tt->SetTriangleContainingTheVertex();
	  kkk++;
	} else if (kk) { // crack + boundary 
	  if (  kkk ) { kc =1;vv = vlast++;  kkk = 0; } // new vertex if use 
	}
	
	ta = Next(ta).Adj(); 
      } while ( (tbegin != ta)); 
      assert(k);
      if (kc)  nbcrakev++;
    }
  
  if ( nbcrakev ) 
      for (int iec =0;iec < NbCrackedEdges; iec ++)
	  CrackedEdges[iec].Set();
  
  //  set the ref 
  cout << " set the ref " <<  endl ;
  NbCrackedVertices =   nbcrakev;
  // int nbvo = nbv;
  nbv = vlast - vertices;
  int nbnewv =  nbv - nbv; // nb of new vrtices 
  if (nbcrakev && verbosity > 1 )
    cout << " Nb of craked vertices = " << nbcrakev << " Nb of created vertices " <<   nbnewv<< endl;
  // all the new vertices are on geometry 
  //  BOFBO--  A VOIR
  if (nbnewv)
    { // 
      Int4 n = nbnewv+NbVerticesOnGeomVertex;
      Int4 i,j,k;
      VertexOnGeom * vog = new VertexOnGeom[n];
      for ( i =0; i<NbVerticesOnGeomVertex;i++) 
	vog[i]=VerticesOnGeomVertex[i];
      delete [] VerticesOnGeomVertex;
      VerticesOnGeomVertex = vog;
      // loop on cracked edge 
      Vertex * LastOld = vertices + nbv - nbnewv;
      for (int iec =0;iec < NbCrackedEdges; iec ++)
	for (k=0;k<2;k++)
	  {
	    Edge & e = *( k ? CrackedEdges[iec].a.edge : CrackedEdges[iec].b.edge);
	    for (j=0;j<2;j++) 
	      {
		Vertex * v = e(j);
		if ( v >=  LastOld)
		  { // a new vertex 
		    Int4 old = v->ReferenceNumber ; // the old same vertex 
		    Int4 i  = ( v - LastOld);
		    //  if the old is on vertex => warning
		    // else the old is on edge => ok 
		    vog[i] = vog[old];
				//  		    vog[i].mv = v;
				//g[i].ge = ;
				//og[i].abcisse = ;
		  }
		
	      }
	  }

	NbVerticesOnGeomVertex = n;
  }
  SetVertexFieldOn();

 
  if (vlast >= vend)
    {  
      cerr << " Not enougth vertices to crack the mesh we need " << nbv << " vertices " << endl;
      MeshError(555,this);
    }
  cout << "  NbCrackedVertices " <<  NbCrackedVertices << endl;
  CurrentTh = CurrentThOld;
  return  NbCrackedVertices;
}
      
Triangles::Triangles(const Triangles & Tho,const int *flag ,const int *bb)
  : Gh(*(new Geometry())), BTh(*this)
{ // truncature
  // 
  
  char cname[] = "trunc";

  int i,k,itadj;
  int kt=0;
  int * kk    = new int [Tho.nbv];
  Int4 * reft = new Int4[Tho.nbt];
  Int4 nbInT =    Tho.ConsRefTriangle(reft);
  Int4 * refv = new Int4[Tho.nbv];

  for (i=0;i<Tho.nbv;i++)
    kk[i]=-1;
  for (i=0;i<Tho.nbv;i++)
    refv[i]=0;
  int nbNewBedge =0;
  //  int nbOldBedge =0;  
  for (i=0;i<Tho.nbt;i++)
    if(  reft[i] >=0 && flag[i]) 
      {
        const Triangle & t = Tho.triangles[i];
        kt++;
        kk[Tho.Number(t[0])]=1;
        kk[Tho.Number(t[1])]=1;
        kk[Tho.Number(t[2])]=1;
        itadj=Tho.Number(t.TriangleAdj(0));
        if (  reft[itadj] >=0 && !flag[itadj])
          { nbNewBedge++;
          refv[Tho.Number(t[VerticesOfTriangularEdge[0][0]])]=bb[i];
          refv[Tho.Number(t[VerticesOfTriangularEdge[0][1]])]=bb[i];
          }
        itadj=Tho.Number(t.TriangleAdj(1));
        if (  reft[itadj] >=0 && !flag[itadj])
          { nbNewBedge++;
          refv[Tho.Number(t[VerticesOfTriangularEdge[1][0]])]=bb[i];
          refv[Tho.Number(t[VerticesOfTriangularEdge[1][1]])]=bb[i];}
        itadj=Tho.Number(t.TriangleAdj(2));
        if (  reft[itadj] >=0 && !flag[itadj])
          { nbNewBedge++;
          refv[Tho.Number(t[VerticesOfTriangularEdge[2][0]])]=bb[i];
          refv[Tho.Number(t[VerticesOfTriangularEdge[2][1]])]=bb[i];}
      }
  k=0;
  for (i=0;i<Tho.nbv;i++)
    if (kk[i]>=0) 
      kk[i]=k++;
  cout << " number of vertices " << k << " remove = " << Tho.nbv - k << endl;
  cout << " number of triangles " << kt << " remove = " << nbInT-kt << endl;
  cout << " number of New boundary edge " << nbNewBedge << endl;
  Int4 inbvx =k;
  PreInit(inbvx,cname);
  for (i=0;i<Tho.nbv;i++)
    if (kk[i]>=0) 
      {
        vertices[nbv] = Tho.vertices[i];
        if (!vertices[nbv].ref())
          vertices[nbv].ReferenceNumber = refv[i];
        nbv++;
      }
  assert(inbvx == nbv);
  for (i=0;i<Tho.nbt;i++)
    if(  reft[i] >=0 && flag[i]) 
      {
        const Triangle & t = Tho.triangles[i];
        int i0 = Tho.Number(t[0]);
        int i1 = Tho.Number(t[1]);
        int i2 = Tho.Number(t[2]);
        assert(i0>=0 && i1 >= 0 && i2  >= 0);
        assert(i0<Tho.nbv && i1 <Tho.nbv && i2  <Tho.nbv);
        // cout <<i<< " F" <<  flag[i] << " T " << nbt << "   = " <<  kk[i0] << " " << kk[i1] << " " << kk[i2] ;
        // cout << " OT  " <<  i0 << " "  << i1 << " " << i2  << " " << reft[i] << endl;
        triangles[nbt] = Triangle(this,kk[i0],kk[i1],kk[i2]);
        triangles[nbt].color = Tho.subdomains[reft[i]].ref; 
        nbt++;           
      }
  assert(kt==nbt);
  if (nbt ==0 && nbv ==0) {
    cout << "Error all triangles was remove " << endl;
    MeshError(999,this);
  }
  delete [] kk;
  delete [] reft;
  delete [] refv;
  double cutoffradian = 10.0/180.0*Pi;
  ConsGeometry(cutoffradian);
  Gh.AfterRead(); 
  SetIntCoor();
  FillHoleInMesh();
   
  assert(NbSubDomains);
  assert(subdomains[0].head && subdomains[0].head->link);
             
}
  
Triangle * Triangles::FindTriangleContening(const I2 & B,Icoor2 dete[3], Triangle *tstart) const
{ // in: B 
  // out: t
  // out : dete[3]
  // t the triangle and s0,s1,s2 the 3 vertices of t
  // in dete[3] = { det(B,s1,s2) , det(s0,B,s2), det(s0,s1,B)}
  // with det(a,b,c ) = -1 if one of 3 vertices a,b,c is NULL 
  Triangle * t=0;	
  int j,jp,jn,jj;
  if (tstart) 
    t=tstart;
  else 
   {
   assert(quadtree);
   Vertex *a = quadtree->NearestVertex(B.x,B.y) ;
  
  if (! a || !a->t ) {
    if (a) 
      {cerr << " Attention PB TriangleConteningTheVertex  vertex number=" << Number(a) << endl;
       cerr  << "We forget a call to ReMakeTriangleContainingTheVertex" << endl;}
    cerr << " Pb with " << B << toR2(B) << endl;
    MeshError(7777);
  }
  assert(a>= vertices && a < vertices+nbv);
#ifdef DRAWING1 
  a->Draw();
#endif 
  //  int k=0;
   t = a->t;
  assert(t>= triangles && t < triangles+nbt);
   
   }
  Icoor2  detop ;
  int kkkk =0; // number of test triangle 

  while ( t->det < 0) 
    { // the initial triangles is outside  
      int k0=(*t)(0) ?  ((  (*t)(1) ? ( (*t)(2) ? -1 : 2) : 1  )) : 0;
      assert(k0>=0); // k0 the NULL  vertex 
      int k1=NextVertex[k0],k2=PreviousVertex[k0];
      dete[k0]=det(B,(*t)[k1],(*t)[k2]);
      dete[k1]=dete[k2]=-1;     
      if (dete[k0] > 0) // outside B 
        return t; 
      t = t->TriangleAdj(OppositeEdge[k0]);
      assert(kkkk++ < 2);
    }

  jj=0;
  detop = det(*(*t)(VerticesOfTriangularEdge[jj][0]),*(*t)(VerticesOfTriangularEdge[jj][1]),B);
 
  while(t->det  > 0 ) 
    { 
      assert( kkkk++ < 2000 ); 
      j= OppositeVertex[jj];
      
#ifdef DRAWING1
      t->Draw();
#endif 
      dete[j] = detop;  //det(*b,*s1,*s2);
      jn = NextVertex[j];
      jp = PreviousVertex[j];
      dete[jp]= det(*(*t)(j),*(*t)(jn),B);
      dete[jn] = t->det-dete[j] -dete[jp];
      
#ifdef DEBUG
      const Vertex * s0 = (*t)(0);
      const Vertex * s1 = (*t)(1);
      const Vertex * s2 = (*t)(2);
      assert(dete[0] == det(B ,*s1,*s2));
      assert(dete[1] == det(*s0,B ,*s2));
      assert(dete[2] == det(*s0,*s1,B ));
      assert(t->det== (dete[0] + dete[1] +dete[2]));
#endif
      // count the number k of  dete <0
      int k=0,ii[3];
      if (dete[0] < 0 ) ii[k++]=0; 
      if (dete[1] < 0 ) ii[k++]=1;
      if (dete[2] < 0 ) ii[k++]=2;
      // 0 => ok
      // 1 => go in way 1
      // 2 => two way go in way 1 or 2 randomly
      
      if (k==0) 
	break;
      if (k == 2 && BinaryRand())
	Exchange(ii[0],ii[1]);
      assert ( k  < 3);
      TriangleAdjacent t1 = t->Adj(jj=ii[0]);
      if ((t1.det() < 0 ) && (k == 2))
	t1 = t->Adj(jj=ii[1]);
      t=t1;
      j=t1;// for optimisation we now the -det[OppositeVertex[j]];
      detop = -dete[OppositeVertex[jj]];
      jj = j;
    }
  
  if (t->det<0) // outside triangle 
    dete[0]=dete[1]=dete[2]=-1,dete[OppositeVertex[jj]]=detop;
  //  NbOfTriangleSearchFind += kkkk;  
  return t;
}

}

