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
#ifdef DRAWING

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <iostream>
using namespace std;

#include "Mesh2.h"
#include "QuadTree.h"
#include "SetOfE4.h"
extern bool withrgraphique;



namespace bamg {
Real4  xGrafCoef,yGrafCoef,xGrafOffSet,yGrafOffSet;
R2 GrafPMin,GrafPMax;
Real8 Grafh=0;



void  IMoveTo(long i,long j)
{  if (!withrgraphique) return;
   rmoveto(float(i)/xGrafCoef+xGrafOffSet ,float(j)/yGrafCoef+yGrafOffSet );}

void  ILineTo(long i,long j)
{      if (!withrgraphique) return;
  rlineto(float(i)/xGrafCoef+xGrafOffSet ,float(j)/yGrafCoef+yGrafOffSet );}


void   Triangles::Draw(   ) const {
    if (!withrgraphique) return;
  showgraphic();
//  if (!init) InitDraw();
  Int4 i;
//  for (i=0;i<nbv;i++)
//    vertices[i].Draw(i);    
  for (i=0;i<nbe;i++)
    edges[i].Draw();
  for (i=0;i<nbt;i++) 
      if (triangles[i].link) 
	triangles[i].Draw(-2);
//  for (i=0;i<nbv;i++)
//   ((Metric) vertices[i]).Draw(vertices[i]);    
   
  //  rattente(1);
  
}


void    Edge::Draw(Int4  i) const 
{
      if (!withrgraphique) return;
  if (InRecScreen(Min(v[0]->r.x,v[1]->r.x),
		  Max(v[0]->r.y,v[1]->r.y),
		  Min(v[0]->r.x,v[1]->r.x),
		  Max(v[0]->r.y,v[1]->r.y))) {
    v[0]->MoveTo();
    v[1]->LineTo();
    R2 M= ((R2) *v[0] + (R2) * v[1])*0.5;
    Move(M);
    char VertexDraw_i10[20];
    if (i<0)
      sprintf(VertexDraw_i10,"%p",this);
    else 
      sprintf(VertexDraw_i10,"%ld",i);
    if (i>=0)
      plotstring(&VertexDraw_i10[0]); 
    
  }
}
void    Vertex::Draw(Int4 i) const 
{
      if (!withrgraphique) return;
  if (CurrentTh && i<0&& i != -2 )
    {
      if (CurrentTh->vertices <= this && this < CurrentTh->vertices + CurrentTh->nbv)
	i = CurrentTh->Number(this);
    }
  if (InPtScreen(r.x,r.y)) {
   char VertexDraw_i10[20];
   if (i<0)
     sprintf(VertexDraw_i10,"%p",this);
   else 
     sprintf(VertexDraw_i10,"%ld",i);
  
  showgraphic();
  //  float eps = (MaxICoor/yGrafCoef)/100;
  DrawMark(r);
  
  if (i>=0)
    plotstring(&VertexDraw_i10[0]); }
}

void  DrawMark(R2 r) 
{
      if (!withrgraphique) return;
  if (InPtScreen(r.x,r.y)) {
  float eps = Grafh/100;
  rmoveto(r.x+eps,r.y);
  rlineto(r.x,r.y+eps);
  rlineto(r.x-eps,r.y);
  rlineto(r.x,r.y-eps);
  rlineto(r.x+eps,r.y);}
 
}
void  Triangle::Draw(Int4 i ) const
{
      if (!withrgraphique) return;
  //  int cc=LaCouleur();
  if (CurrentTh && i<0 && i != -2)
    {
      if (CurrentTh->triangles <= this && this < CurrentTh->triangles + CurrentTh->nbt)
	i = CurrentTh->Number(this);
    }	
  char i10[20];
  if (i<0)   sprintf(i10,"%p",this);
  else  sprintf(i10,"%ld",i);
  showgraphic();

  if (ns[0] == 0) {
  if (InRecScreen(Min(ns[1]->r.x,ns[2]->r.x),
		  Min(ns[1]->r.y,ns[2]->r.y),
		  Max(ns[1]->r.x,ns[2]->r.x),
		  Max(ns[1]->r.y,ns[2]->r.y)))
    { rmoveto(ns[1]->r.x,ns[1]->r.y);
      rlineto(ns[2]->r.x,ns[2]->r.y);
      rmoveto( (ns[1]->r.x + ns[2]->r.x)/2.0 + (ns[1]->r.y - ns[2]->r.y)*0.1,
              (ns[1]->r.y + ns[2]->r.y)/2.0 - (ns[1]->r.x - ns[2]->r.x)*0.1 );
      if(i>=0)  plotstring(&i10[0]);    
    }}
  else
  if (InRecScreen(Min3(ns[0]->r.x,ns[1]->r.x,ns[2]->r.x),
		  Min3(ns[0]->r.y,ns[1]->r.y,ns[2]->r.y),
		  Max3(ns[0]->r.x,ns[1]->r.x,ns[2]->r.x),
		  Max3(ns[0]->r.y,ns[1]->r.y,ns[2]->r.y))) {
    { 
      const int i0=0                               , j01=EdgesVertexTriangle[i0][1];//     j01==2
      const int i1=VerticesOfTriangularEdge[j01][1], j12=EdgesVertexTriangle[i1][1];//i1=1,j12=0
      const int i2=VerticesOfTriangularEdge[j12][1], j20=EdgesVertexTriangle[i2][1];//i1=2,j20=1

      rmoveto(ns[i0]->r.x,ns[i0]->r.y);
      if(Hidden(j01)) rmoveto(ns[i1]->r.x,ns[i1]->r.y); else  rlineto(ns[i1]->r.x,ns[i1]->r.y);
      if(Hidden(j12)) rmoveto(ns[i2]->r.x,ns[i2]->r.y); else  rlineto(ns[i2]->r.x,ns[i2]->r.y);
      if(Hidden(j20)) rmoveto(ns[i0]->r.x,ns[i0]->r.y); else  rlineto(ns[i0]->r.x,ns[i0]->r.y);
      
    rmoveto( (ns[0]->r.x + ns[1]->r.x + ns[2]->r.x)/3.0 ,
             (ns[0]->r.y + ns[1]->r.y + ns[2]->r.y)/3.0);
    
    }
  
  if(i>=0)  plotstring(&i10[0]);    
  }
  //  LaCouleur(cc);
}

void  Triangles::InitDraw() const
{ 
   if (!withrgraphique) return; 
   couleur(1);
   GrafPMin =  vertices[0].r;
   GrafPMax =  vertices[0].r;
  // recherche des extrema des vertices GrafPMin,GrafPMax
  Int4 i;
  for (i=0;i<nbv;i++) {
    //    GrafPMax = Max2(GrafPMax,vertices[i].r);
    //    GrafPMin = Min2(GrafPMin,vertices[i].r);
    GrafPMin.x = Min(GrafPMin.x,vertices[i].r.x);
    GrafPMin.y = Min(GrafPMin.y,vertices[i].r.y);
    GrafPMax.x = Max(GrafPMax.x,vertices[i].r.x);
    GrafPMax.y = Max(GrafPMax.y,vertices[i].r.y);
  }
  float hx = (GrafPMax.x-GrafPMin.x);
  float hy = (GrafPMax.y-GrafPMin.y);
  Grafh = Max(hx,hy)*0.55;
  cadreortho((GrafPMin.x+GrafPMax.x)*0.5,(GrafPMin.y+GrafPMax.y)*0.5,Grafh);
  
}
void Show(const char * s,int k=1)
{
  if (!withrgraphique) return;
  if(k) {
  couleur(1);
  float xmin,xmax,ymin,ymax;
       getcadre(xmin,xmax,ymin,ymax);
       rmoveto(xmin+(xmax-xmin)/100,ymax-(k)*(ymax-ymin)/30);
       plotstring(s);
       //  couleur(1);	
  }
}
void Triangles::inquire() 
{    
  if (! withrgraphique) return;
  int PS=0;
  cout << flush;
   Triangles *OldCurrentTh=CurrentTh;
  CurrentTh = this;
  if (! Grafh ) InitDraw();
  unsigned  char c;
  int rd;
  float x,y;
  int setv =0;
  Int4 nbtrel=0;
  Int4 * reft = new Int4 [Max(nbt,1L)];
  if (NbSubDomains && subdomains[0].head)
    nbtrel = ConsRefTriangle(reft);
  else
    for (Int4 kkk=0;kkk<nbt;kkk++) reft[kkk]=0;
  // cout << "inquire **********************************************************************"<<endl ;
      if (verbosity>2)
	cout << "inquire: Nb de Triangle reel = " <<ConsRefTriangle(reft)<<endl;

  while (Show("Enter ? for help",PS==0), c=::Getxyc(x,y), ( c && c != 'F'  && c != 'f')  && c < 250 )
    { rd = 0;
      couleur(1);
    //  cout << " #"<< c << "# " << (int) c << " xy=" << x << " "<< y  <<endl ;
      //     penthickness(2);
      if (c=='?') { int i=3;
        reffecran();
	Show("enter a keyboard character in graphic window to do",i++);
	i++;
	Show("f or mouse click: to continue the process",i++);
	Show("p: openPS file to save plot, P: close PS file ",i++);
        Show("B: set backgound mesh has drawing mesh",i++);
	Show("H: show histogramme",i++);
	Show("i: initDraw",i++);
	Show("+ or - : zoom + or -",i++);
	Show("r: redraw  , =: reinit the viewport ",i++);
	Show("q: quit / abort / stop ",i++);
	Show("Q: show quadtree",i++);
	Show("g: show geometry",i++);
	Show("d: show triangle with det = 0",i++);
	Show("m: show all the metric",i++);
	Show("V: show all the vertices",i++);
	Show("T: show all the sub domain ref number",i++);
	Show("V: show all the vertices",i++);
	Show("D: Inforce the mesh to be Delaunay",i++);
	Show("b: show all the boundary edges ref number",i++);
	Show("k: find the triangle contening the point",i++);
	Show("v,o,s : print the nearest vertex  s=> draw metric o=> optim around (debug)",i++);
	Show("t: find triangle contening the point  with brute force",i++);	
	Show("e: find the nearest edge of triangle contening the point",i++);
	Show("C: construct the in-circle in anistrope way",i++);
	Show("n: show normal for find ",i++);
      }
      if (c=='p') openPS(0),PS=1;
      if (c=='P') closePS(),PS=0;
      if (c=='B') BTh.inquire();
      if ( c=='D') {
	for(int iter=0;iter < 50;iter++)
	  { int k = 0;
	  for (Int4 icount=0; icount<nbv; icount++) {
	    k += vertices[icount].Optim(1,0);
	  }
	  if (k !=0) cout << " Bizarre le maillage n'est  pas delaunay  nbswap = " << k << endl;
	  else break;
	  }
	rd =1;
      }
      
      if (c=='H') ShowHistogram();
      if (c=='i') {InitDraw();
      }
      if (c=='+') Grafh /= 2.,rd=1,cadreortho(x,y,Grafh);
      if (c=='-') Grafh *= 2.,rd=1,cadreortho(x,y,Grafh);
      if (c=='r') rd =1;
      if (c=='Q') CurrentTh=0,MeshError(1);
      if (c=='q' && quadtree) penthickness(2),
				couleur(6), 
				quadtree->Draw(),
				couleur(1),
				penthickness(1);
      if (c=='g') couleur(6),Gh.Draw();
      if (c=='n') {
        for (int i=0;i<nbv;i++)
         {
           vertices[i].MoveTo();
           vertices[1].DirOfSearch.Draw();;
         }
      }
      if (c=='d')
	{      
          couleur(4);
	  for(int i=0;i<nbt;i++)
	    if(triangles[i].det == 0) {
	      triangles[i].Draw();
	      cout << " Bizzare " << i << endl;
	    }
          couleur(1);
       
	}
      if (c=='m')  {Int4 i;
        couleur(2);
        for (i=0;i<nbv;i++)
	((Metric) vertices[i]).Draw(vertices[i]); 
        couleur(1);}
      if (c=='V')  {Int4 i;
        for (i=0;i<nbv;i++)
	 vertices[i].Draw(); }
      if (c=='T') {
	Int4 i;
	for( i=0;i<nbt;i++)
	  if(reft[i]>=0)
	    { couleur(2+(reft[i])%6);
	      triangles[i].Draw(reft[i]);
              couleur(1); 
	    }
      }
	
      if (c=='b') {
	Int4 i;
	reffecran();
	for (i=0;i<nbe;i++)
	  edges[i].Draw(edges[i].ref);
	
      }
      if (c=='k') 
	{ 
	  if(!setv) ReMakeTriangleContainingTheVertex();
	  R2 P(x,y);
	  I2 I = toI2(P);
	  Icoor2 dete[3];
	  Triangle * tb = FindTriangleContening(I,dete);
	  penthickness(1);
	  tb->Draw(Number(tb));
	  penthickness(1);
	  if(tb->link==0) {
	    double a,b;
	    TriangleAdjacent ta=CloseBoundaryEdgeV2(I,tb,a,b);
	    R2 A = *ta.EdgeVertex(0);
	    R2 B = *ta.EdgeVertex(1);
	    //Triangle * tt  = ta;
	    //    tt->Draw(Number(tt));
	    penthickness(5);
	    //   ta.EdgeVertex(0)->MoveTo();
	    //  ta.EdgeVertex(1)->LineTo();
	    DrawMark(A*a+B*b);
	    penthickness(1);
	      }
	}
      if (c=='v' || c=='o' || c=='s' )
	{ 
	  couleur(3);
	  if(!setv) ReMakeTriangleContainingTheVertex();
	  setv=1;
	  
	  R2 XY(x,y),P;
	  Real8 d;
	  Int4 i,j=0;
	  P  = XY - vertices[0].r;
	  d = (P,P);
	  for (i=0;i<nbv;i++) {
	    P  = XY - vertices[i].r;
	  Real8 dd = (P,P);
	  if(dd < d) {
	    j=i;
	    d=dd;}}
	  cout << " sommet " << j << "= " <<	vertices[j]  << ", d = " << d << " " << vertices[j].m << endl;
	  vertices[j].Draw(j);
	  DrawMark(vertices[j].r);
	  if( c=='s')  vertices[j].m.Draw(vertices[j]);
	  if (c=='o')   cout << " Nb Swap = " << vertices[j].Optim(1) << endl;
	  couleur(1);
      }
      if (c=='t' || c=='c'||c=='e') {
	    couleur(4);
            R2 XY(x,y);
	    Real8 a12,a02,a01;
	    for (Int4 i=0;i<nbt;i++) {
	      if(triangles[i].det<=0)  continue;
	      //	      cout << " T = " << i << " " << triangles[i].det << " ";
	      //	      cout << Number(triangles[i][0]) << " "
	      //		   << Number(triangles[i][1]) << " "
	      //		   << Number(triangles[i][2]) << " ";
	      //    cout << area2(triangles[i][0].r,triangles[i][1].r,triangles[i][2].r) << " " ;
	      //	      cout << area2(XY,triangles[i][1].r,triangles[i][2].r) << " " ;
	      //      cout << area2(triangles[i][0].r,XY,triangles[i][2].r) << " " ;
	      //     cout << area2(triangles[i][0].r,triangles[i][1].r,XY) << endl;
			   
	      if( (a12=Area2(XY,triangles[i][1].r,triangles[i][2].r)) < 0) continue;
	      if( (a02=Area2(triangles[i][0].r,XY,triangles[i][2].r)) < 0) continue;
	      if( (a01=Area2(triangles[i][0].r,triangles[i][1].r,XY)) < 0) continue;
	      if(c=='e'|| c=='E') {
		int ie =0;
		Real8 am = a12;
		if(a02 < am) am=a02,ie=1;
		if(a01 < am) am=a01,ie=2;
		TriangleAdjacent tta(triangles+i,ie);
		Vertex *v0 =tta.EdgeVertex(0);
		Vertex *v1 =tta.EdgeVertex(1);
		tta.EdgeVertex(0)->MoveTo();
		tta.EdgeVertex(1)->LineTo();
		cout << " Edge " << Number(tta.EdgeVertex(0)) << " " <<  Number(tta.EdgeVertex(1)) << endl;
		for (Int4 k=0;k<nbe;k++)
		  if ( ( edges[k](0) == v0 &&  edges[k](1) == v1) ||
		       ( edges[k](0) == v1 &&  edges[k](1) == v0)	)			
		    cout << " Edge " << k << "  on Geo = " << Gh.Number(edges[k].on)<< endl;
		
		  
		if(c=='e') {
		  triangles[i].SetUnMarkUnSwap(ie);
		  triangles[i].swapDRAW(ie);}
		else if(c=='E')
		  {
		    ;
		  }

		break;
	      }
	      cout << " In triangle " << i << triangles[i];
	      triangles[i].Draw(i);
	      if(c=='c') {
		Vertex *s1=&triangles[i][0];
		Vertex *sa=&triangles[i][1];
		Vertex *sb=&triangles[i][2];
		D2 S1(s1->r.x,s1->r.y);
		D2 SA(sa->r.x,sa->r.y);
		D2 SB(sb->r.x,sb->r.y);
		D2 AB= SB-SA;
		D2 MAB2=SB + SA;
		D2 M1A=(S1 + SA)*0.5;
		D2 MAB(MAB2.x*0.5,MAB2.y*0.5);
		D2 A1=S1-SA;
		D2 D = S1 - SB ;
		{
		  Metric M=s1->m;
		  D2 ABo = M.Orthogonal(AB);
		  D2 A1o = M.Orthogonal(A1);
		  penthickness(1);
		  Move(MAB);
		  Line(MAB+ABo);
		  Line(MAB-ABo);
		  Move(M1A);
		  Line(M1A+A1o);
		  Line(M1A-A1o);
		  penthickness(3);
		  
		  // (A+B)+ x ABo = (S1+B)/2+ y A1 
		  // ABo x - A1o y =  (S1+B)/2-(A+B)/2 = (S1-B)/2 = D/2
		  double dd = Abs(ABo.x*A1o.y)+Abs(ABo.y*A1o.x);
		  double d = (ABo.x*A1o.y - ABo.y*A1o.x)*2; // because D/2
		  cout << " d = " << d << " dd= " << dd << endl;
		  if (Abs(d) > dd*1.e-10) {
		    D2 C(MAB+ABo*((D.x*A1o.y - D.y*A1o.x)/d));
		    cout << C << s1->r <<sa->r  <<sb->r << endl;
		    DrawMark(C);
		    cout << M << " l = "<<  M(C-S1) << " lnew = 1 == " ;
		    //M.Draw(R2(C.x,C.y));
		    M=M/M(C-S1);	
		    cout << M(C-S1) << M << endl;	 
		    M.Draw(R2(C.x,C.y));
		  }
		  
		}
	      }
	      break;
	    }
	    cout << " fin recherche triangle " << endl;
	   couleur(1); 
      }            
      if (c=='='||c==249) {
	rd=1;
	float hx = (GrafPMax.x-GrafPMin.x);
	float hy = (GrafPMax.y-GrafPMin.y);
	Grafh = Max(hx,hy)*0.55;
	cadreortho((GrafPMin.x+GrafPMax.x)*0.5,(GrafPMin.y+GrafPMax.y)*0.55,Grafh);
      }
       penthickness(1);
      if (rd) 
        reffecran(),Draw();

    }
  // cout << endl;
  delete [] reft;
  CurrentTh = OldCurrentTh;
 }
void Draw(float x,float y)
{
 if (!withrgraphique) return;

  if (InPtScreen(x,y)) {
  float eps = Max(GrafPMax.x-GrafPMin.x,GrafPMax.y-GrafPMin.y)/100;
  rmoveto(x+eps,y);
  rlineto(x,y+eps);
  rlineto(x-eps,y);
  rlineto(x,y-eps);
  rlineto(x+eps,y);}
}
char Getxyc(long &i,long &j)
{ 
   if (!withrgraphique) return 0;
 
   float x,y;
   char c=::Getxyc(x,y);
   i = (long)( (x - xGrafOffSet)*xGrafCoef);
   j = (long)((y - yGrafOffSet)*yGrafCoef);
   return c;
}


void Draw(long i,long j)
{
  if (!withrgraphique) return;
  Draw(((float) i)/xGrafCoef+xGrafOffSet ,((float) j)/yGrafCoef+yGrafOffSet );
}

int Triangle::swapDRAW(Int2 a){
  int NbUnSwap=0;
  if(a/4 !=0)  {cout << a << "arete lock"<< a <<endl;return 0;}// arete lock or MarkUnSwap

  register Triangle *t1=this,*t2=at[a];// les 2 triangles adjacent
  register Int1 a1=a,a2=aa[a];// les 2 numero de l arete dans les 2 triangles
  
  if(a2/4 !=0) {cout << a2 << "arete lock adj"<< a2 << endl; return 0;} // arete lock
  
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
     cout << " detMin = " << detMin << "  detMinNew " <<  detMinNew << endl;
     if (! OnSwap &&(detMinNew>0)) {
       OnSwap = detMin ==0;
       if (! OnSwap) {
	 int  kopt = 1;
        while (1)
	 if(kopt) {
	 // critere de Delaunay pure 
	 register Real8 xb1 = sb->i.x - s1->i.x,
	   x21 = s2->i.x - s1->i.x,
	   yb1 = sb->i.y - s1->i.y,
	   y21 = s2->i.y - s1->i.y,
	   xba = sb->i.x - sa->i.x, 
	   x2a = s2->i.x - sa->i.x,
	   yba = sb->i.y - sa->i.y,
	   y2a = s2->i.y - sa->i.y,
	   cosb12 =  xb1*x21 + yb1*y21 ,
	   cosba2 =  xba*x2a + yba*y2a ,
	   sinb12 = det2,
	   sinba2 = t2->det,
	   //zsinba2 = t2->det;
	 
	 // angle b12 > angle ba2 => cotg(angle b12) < cotg(angle ba2)
	 OnSwap =  (cosb12 *  sinba2) <  (cosba2 *  sinb12);
	 if(CurrentTh) 
	   cout << "swap s1=" << CurrentTh->Number(sa) << " s2=" << CurrentTh->Number(sb) 
		<< " t1= " <<  CurrentTh->Number(t1) << " t2=" <<  CurrentTh->Number(t2) << " " ;


	 cout <<  cosb12 << " " <<  sinba2 << " "  <<  cosba2 << " " << sinb12 
	      << " Onswap = " <<  OnSwap << endl;

	  break;
	 }
	 else 
	   {	
	     // critere de Delaunay pure 
	     Real8 som;
	     I2 AB=(I2) *sb - (I2) *sa;
	     I2 MAB2=((I2) *sb + (I2) *sa);
	     D2 MAB(MAB2.x*0.5,MAB2.y*0.5);
	     I2 A1=(I2) *s1 - (I2) *sa;
	     I2 D = (I2) * s1 - (I2) * sb ;
	     D2 S2(s2->i.x,s2->i.y);
	     D2 S1(s1->i.x,s1->i.y);
	     D2 SA(sa->i.x,sa->i.y);
	     D2 SB(sb->i.x,sb->i.y);
	     DrawMark(s2->r);
	     DrawMark(s1->r);
	     DrawMark(sa->r);
	     DrawMark(sb->r);
	     {
	       Metric M=s1->m;
	       D2 ABo = M.Orthogonal(SB-SA);
	       D2 A1o = M.Orthogonal(S1-SA);
	       // (A+B)+ x ABo = (S1+B)/2+ y A1 
	       // ABo x - A1o y =  (S1+B)/2-(A+B)/2 = (S1-B)/2 = D/2
	       double dd = Abs(ABo.x*A1o.y)+Abs(ABo.y*A1o.x);
	       double d = (ABo.x*A1o.y - ABo.y*A1o.x)*2; // because D/2
	       if (Abs(d) > dd*1.e-3) {
		 D2 C(MAB+ABo*((D.x*A1o.y - D.y*A1o.x)/d));
		 cout << "M1 r2 =" <<  M(C - S2) << " r1 = " << M(C - S1) 
		      << "ra = " << M(C - SA) << " rb = " << M(C-SB) ;
		 som  = M(C - S2)/M(C - S1);
		 cout << " r1/r2 = " << som << endl;
	       } else 
		{kopt=1;continue;}
		
	     }
	     {
	       Metric M=s2->m;
	       D2 ABo = M.Orthogonal(SB-SA);
	       D2 A1o = M.Orthogonal(S1-SA);
	       // (A+B)+ x ABo = (S1+B)/2+ y A1 
	       // ABo x - A1o y =  (S1+B)/2-(A+B)/2 = (S1-B)/2 = D/2 
	       double dd = Abs(ABo.x*A1o.y)+Abs(ABo.y*A1o.x);
	       double d = (ABo.x*A1o.y - ABo.y*A1o.x)*2; // because D/2
	       cout << " d = " << Abs(d) << " dd " << dd << endl;
	       if(Abs(d) > dd*1.e-3) {
		 D2 C(MAB+ABo*((D.x*A1o.y - D.y*A1o.x)/d));
		 cout << "M2 r1 =" <<  M(C - S2) << " r2 = " << M(C - S1) 
		      << "ra = " << M(C - SA) << " rb = " << M(C-SB) 
		      << " r1/r2 =  " << M(C - S2)/M(C -  S1)  << endl;
		 som  += M(C - S2)/M(C -  S1);
	       } else 
		{kopt=1;continue;}
	     }
	     {
	       Metric M=sa->m;
	       D2 ABo = M.Orthogonal(SB-SA);
	       D2 A1o = M.Orthogonal(S1-SA);
	       // (A+B)+ x ABo = (S1+B)/2+ y A1 
	       // ABo x - A1o y =  (S1+B)/2-(A+B)/2 = (S1-B)/2 = D/2 
	       double dd = Abs(ABo.x*A1o.y)+Abs(ABo.y*A1o.x);
	       double d = (ABo.x*A1o.y - ABo.y*A1o.x)*2; // because D/2
	       cout << " d = " << Abs(d) << " dd " << dd << endl;
	       if(Abs(d) > dd*1.e-3) {
		 D2 C(MAB+ABo*((D.x*A1o.y - D.y*A1o.x)/d));
		 cout << "M2 r1 =" <<  M(C - S2) << " r2 = " << M(C - S1) 
		      << "ra = " << M(C - SA) << " rb = " << M(C-SB) 
		      << " r1/r2 =  " << M(C - S2)/M(C -  S1)  << endl;
		 som  += M(C - S2)/M(C -  S1);
	       } else 
		{kopt=1;continue;}
	     }
  	     {
	       Metric M=sb->m;
	       D2 ABo = M.Orthogonal(SB-SA);
	       D2 A1o = M.Orthogonal(S1-SA);
	       // (A+B)+ x ABo = (S1+B)/2+ y A1 
	       // ABo x - A1o y =  (S1+B)/2-(A+B)/2 = (S1-B)/2 = D/2 
	       double dd = Abs(ABo.x*A1o.y)+Abs(ABo.y*A1o.x);
	       double d = (ABo.x*A1o.y - ABo.y*A1o.x)*2; // because D/2
	       cout << " d = " << Abs(d) << " dd " << dd << endl;
	       if(Abs(d) > dd*1.e-3) {
		 D2 C(MAB+ABo*((D.x*A1o.y - D.y*A1o.x)/d));
		 cout << "M2 r1 =" <<  M(C - S2) << " r2 = " << M(C - S1) 
		      << "ra = " << M(C - SA) << " rb = " << M(C-SB) 
		      << " r1/r2 =  " << M(C - S2)/M(C -  S1)  << endl;
		 som  += M(C - S2)/M(C -  S1);
	       } else 
		{kopt=1;continue;}
	     }
	     cout <<  som  << endl;
	     OnSwap = som < 4;
	     break;
	 }

       } // OnSwap 
     } // (! OnSwap &&(det1 > 0) && (det2 > 0) )
   }
   cout << OnSwap << endl;
   if( OnSwap ) {
     couleur(0);
     t1->Draw();
     t2->Draw();    
     bamg::swap(t1,a1,t2,a2,s1,s2,det1,det2);		
     couleur(1);
     t1->Draw();
     t2->Draw();    
   }
   else {
     NbUnSwap ++;
     t1->SetMarkUnSwap(a1);     
   }
   return OnSwap;
}


	   

void GeometricalEdge::Draw(Int4  i)
{ 
    if (!withrgraphique) return;

  if (CurrentTh && i<0 && i != -2)
    {
      if (CurrentTh->Gh.edges <= this && this < CurrentTh->Gh.edges + CurrentTh->Gh.nbe)
	i = CurrentTh->Gh.Number((this));
    }	

    v[0]->MoveTo();
    R2 x,x50;
    int k=0,k50=0;
    for (int ii = 0;ii<100;ii++) {
      x= F( Real4(ii)/100.0);
      if(ii==50) x50=x,k50=1;
      if (InPtScreen(x.x,x.y) )
	if(k) rlineto(x.x,x.y);
	else k=1,rmoveto(x.x,x.y);
    }
      if (InPtScreen(v[1]->r.x,v[1]->r.y))      
	v[1]->LineTo();  

      char VertexDraw_i10[20];
      if( k50) {
	if (i<0)
	  sprintf(VertexDraw_i10,"Eg%p",this);
	else 
	  sprintf(VertexDraw_i10,"Eg%ld",i);
	rmoveto(x50.x,x50.y);
	plotstring(&VertexDraw_i10[0]); }  

}

void  Geometry::InitDraw() const
{ 
   GrafPMin =  vertices[0].r;
   GrafPMax =  vertices[0].r;
  // recherche des extrema des vertices GrafPMin,GrafPMax
  Int4 i;
  for (i=0;i<nbv;i++) {
    //    GrafPMax = Max2(GrafPMax,vertices[i].r);
    //    GrafPMin = Min2(GrafPMin,vertices[i].r);
    GrafPMin.x = Min(GrafPMin.x,vertices[i].r.x);
    GrafPMin.y = Min(GrafPMin.y,vertices[i].r.y);
    GrafPMax.x = Max(GrafPMax.x,vertices[i].r.x);
    GrafPMax.y = Max(GrafPMax.y,vertices[i].r.y);
  }
  float hx = (GrafPMax.x-GrafPMin.x);
  float hy = (GrafPMax.y-GrafPMin.y);
  Grafh = Max(hx,hy)*7/10;
  if (withrgraphique)
   cadreortho((GrafPMin.x+GrafPMax.x)*0.5,(GrafPMin.y+GrafPMax.y)*0.55,Grafh);

  
}

void   Geometry::Draw() const {
    if (!withrgraphique) return;

  showgraphic();
  //  InitDraw();
  Int4 i;
  for (i=0;i<nbv;i++)
    vertices[i].Draw();    
  for (i=0;i<nbe;i++)
    {
	if(edges[i].Cracked()) couleur(4);
	else couleur(6);
    edges[i].Draw(i);
    }
   couleur(6);
  for (i=0;i<nbv;i++)
    if (vertices[i].Required()) {
      char i10[40];
      sprintf(i10,"%ld:%d",i,vertices[i].Required());
      Move(vertices[i].r);
      if(vertices[i].Corner()) couleur(2);
      plotstring(i10);  
      if(vertices[i].Corner()) couleur(6);
      }
  // rattente(1);
 }

}
#endif
