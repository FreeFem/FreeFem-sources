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

#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <time.h>
#if  (defined(unix) || defined(__unix)) && !defined(__AIX)
#define SYSTIMES
#include <sys/times.h>
#include <unistd.h>
#endif
#ifdef  DRAWING
#include "rgraph.hpp"
#endif

extern long verbosity;
extern int SHOW;
#include "meshtype.h"

#include "error.hpp"


#include "R2.h"

namespace bamg {



const  double Pi =  3.14159265358979323846264338328;
const  float fPi =  3.14159265358979323846264338328;


class MeshIstream;
class OFortranUnFormattedFile;
class IFortranUnFormattedFile;

extern int hinterpole;


typedef P2<Icoor1,Icoor2> I2;

inline int BinaryRand(){
#ifdef RAND_MAX
 const long HalfRandMax = RAND_MAX/2;
 return rand() <HalfRandMax;
#else
 return rand() & 16384; // 2^14 (for sun because RAND_MAX is not def in stdlib.h)
#endif

} 
typedef P2<Real8,Real8> R2;
typedef P2xP2<Int2,Int4> I2xI2;
typedef P2<Real4,Real8> R2xR2;

}

#include "Metric.h"

namespace bamg {
inline float OppositeAngle(float a)
 {return a<0 ? fPi + a :a - fPi ;}
inline double OppositeAngle(double a)
 {return a<0 ? Pi + a :a - Pi ;}
 
#ifdef DRAWING
extern Real4  xGrafCoef,yGrafCoef,xGrafOffSet,yGrafOffSet;
extern R2 GrafPMin,GrafPMax;
extern Real8 Grafh;
#endif

Icoor2 inline det(const I2 &a,const I2 & b,const I2 &c)
{
  register  Icoor2 bax = b.x - a.x ,bay = b.y - a.y; 
  register  Icoor2 cax = c.x - a.x ,cay = c.y - a.y; 
  return  bax*cay - bay*cax;}



// def de numerotation dans un triangles 
static const Int2 VerticesOfTriangularEdge[3][2] = {{1,2},{2,0},{0,1}};
static const Int2 EdgesVertexTriangle[3][2] = {{1,2},{2,0},{0,1}};
static const Int2 OppositeVertex[3] = {0,1,2};
static const Int2 OppositeEdge[3] =  {0,1,2};
static const Int2 NextEdge[3] = {1,2,0};
static const Int2 PreviousEdge[3] = {2,0,1};
static const Int2 NextVertex[3] = {1,2,0};
static const Int2 PreviousVertex[3] = {2,0,1};

Int4 AGoodNumberPrimeWith(Int4 n);

// remark all the angle are in radian beetwen [-Pi,Pi]


class Geometry;
//static Geometry *NULLGeometry=0;
class Triangles;
class Triangle;
class QuadTree;
class GeometricalEdge;
class VertexOnGeom;
class VertexOnEdge;
/////////////////////////////////////////////////////////////////////////////////////
const int IsVertexOnGeom = 8;
const int IsVertexOnVertex = 16;
const int IsVertexOnEdge = 32;
/////////////////////////////////////////////////////////////////////////////////////
#ifndef  NOTFREEFEM
class ErrorMesh : public Error
{  
public:
    Triangles *Th;
    ErrorMesh(const char * Text,int l,Triangles * TTh=0, const char *t2="") :
	Error(MESH_ERROR,"Meshing error: ",Text,"\n number : ",l,", ",t2),Th(TTh)  {}
};
#endif

class Direction { //   
  private:
  Icoor1 dir;
  public:
  Direction(): dir(MaxICoor){}; //  no direction set
  Direction(Icoor1 i,Icoor1 j) { Icoor2  n2 = 2*(Abs(i)+Abs(j));  
                                 Icoor2 r = MaxICoor* (Icoor2) i;
                                 Icoor1 r1 = (Icoor1) (2*(r/ n2)); // odd number 
                                 dir = (j>0) ? r1 : r1+1; //  odd -> j>0 even -> j<0
                               }
  int sens(    Icoor1 i,Icoor1 j) { int r =1; 
                                   if (dir!= MaxICoor) {
                                     Icoor2 x(dir/2),y1(MaxICoor/2-Abs(x)),y(dir%2?-y1:y1);
                                     r = (x*i + y*j) >=0;}
                                   return r;}
#ifdef DRAWING
 void  Draw() {
                         if (dir!= MaxICoor) {
                                     Icoor2 x(dir/2),y1(MaxICoor/2-Abs(x)),y(dir%2?-y1:y1);
                                     R2 D(x,y);
                                     double eps = Grafh/Norme2(D)/20;
                                     D = D*eps;
                                     rmoveto(D.x,D.y);
                                    }
  }
#endif                                   
                                   
      
                              
  
};
/////////////////////////////////////////////////////////////////////////////////////
class Vertex {public:
  I2 i;  // allow to use i.x, and i.y in long int (beware of scale and centering)
  R2 r;  // allow to use r.x, and r.y in double
  Metric m;
  Int4 ReferenceNumber;
  Direction DirOfSearch;
  union {
    Triangle * t; // one triangle which contained  the vertex
    Int4 color;
    Vertex * to;// use in geometry Vertex to now the Mesh Vertex associed 
    VertexOnGeom * on;     // if vint 8; // set with Triangles::SetVertexFieldOn()
    Vertex * onbv; // if vint == 16 on Background vertex Triangles::SetVertexFieldOnBTh()
    VertexOnEdge * onbe;   // if vint == 32 on Background edge
  };
  Int1 vint;  // the vertex number in triangle; varies between 0 and 2 in t
  operator  I2   () const {return i;} // operator de cast 
  operator  const R2 & () const {return r;}// operator de cast 
//  operator  R2 & () {return r;}// operator de cast 
  Real8 operator()(R2 x) const { return m(x);}
  operator Metric () const {return m;}// operator de cast 
  Int4  Optim(int  = 1,int =0); 
  //  Vertex(){}
  //  ~Vertex(){}
  Real8  Smoothing(Triangles & ,const Triangles & ,Triangle  * & ,Real8 =1);
  int ref() const { return ReferenceNumber;}

  friend ostream& operator <<(ostream& f, const  Vertex & v)
  {f << "(" << v.i  << "," << v.r << MatVVP2x2(v.m) << ")" ;   return f;}
  inline void Set(const Vertex & rec,const Triangles &,Triangles &);
  
#ifdef DRAWING
  void  Draw(Int4 =-1) const ;
  void MoveTo() const  {    rmoveto(r.x,r.y);  }
  void LineTo() const {    rlineto(r.x,r.y);  }
#endif  
};

double QuadQuality(const Vertex &,const Vertex &,const Vertex &,const Vertex &);

// extern Vertex *Meshend , *Meshbegin;

/////////////////////////////////////////////////////////////////////////////////////
class TriangleAdjacent {
  friend ostream& operator <<(ostream& f, const  TriangleAdjacent & ta)
  {f << "{" << ta.t << "," << ((int) ta.a) << "}" ;
   return f;}

public:
  Triangle * t; // le triangle 
  int  a; // le numero de l arete
  
  TriangleAdjacent(Triangle  * tt,int  aa): t(tt),a(aa &3) {};
  TriangleAdjacent() {};
  
  operator Triangle * () const {return t;}
  operator Triangle & () const {return *t;}
    operator int() const {return a;}
  TriangleAdjacent & operator++() 
  {
    a= NextEdge[a];
    return *this;}
  TriangleAdjacent operator--()
  { 
    a= PreviousEdge[a];
    return *this;}
  inline  TriangleAdjacent  Adj() const ;
  int swap();
  inline void SetAdj2(const TriangleAdjacent& , int =0);
  inline Vertex *  EdgeVertex(const int &) const ;
  inline Vertex *  OppositeVertex() const ;
  inline Icoor2 & det() const;
  inline int Locked() const  ;
  inline int GetAllFlag_UnSwap() const ;
  inline void SetLock();
  inline int MarkUnSwap()  const;
  inline void SetMarkUnSwap();
  inline void SetCracked();
  inline int Cracked() const ;
};// end of Vertex class  


/////////////////////////////////////////////////////////////////////////////////////
class Edge { public:
   Vertex * v[2];
   Int4 ref;
    GeometricalEdge * on;
   Vertex & operator[](int i){return *v[i];};
   Vertex * operator()(int i){return v[i];};

  void ReNumbering(Vertex *vb,Vertex *ve, Int4 *renu) 
   {
    if (v[0] >=vb && v[0] <ve) v[0] = vb + renu[v[0]-vb];
    if (v[1] >=vb && v[1] <ve) v[1] = vb + renu[v[1]-vb];
   }

  const Vertex & operator[](int i) const { return *v[i];};
  R2 operator()(double t) const; // return the point 
  //                                on the curve edge a t in [0:1]
  Edge * adj[2]; // the 2 adj edges if on the same curve 
  int Intersection(const  Edge & e) const { 
    if (!(adj[0]==&e || adj[1]==&e)) 
      cerr << "Bug : Intersection " << (void*) &e <<  "  " 
	   << adj[0] << " " <<  adj[1] << endl;
    assert(adj[0]==&e || adj[1]==&e);
    return adj[0]==&e ? 0 : 1;}
  Real8 MetricLength() const ;  
  inline void Set(const Triangles &,Int4,Triangles &);
  
#ifdef DRAWING
  void  Draw(Int4 = -1) const ;
#endif
}; // end of Edge class 

/////////////////////////////////////////////////////////////////////////////////////
class GeometricalVertex :public Vertex {
  int cas;
  friend class Geometry;
  GeometricalVertex * link; //  link all the same GeometricalVertex circular (Crack) 
public:
  int Corner() const {return cas&4;}
  int Required()const {return cas&6;}// a corner is required
  void  SetCorner(){ cas |= 4;}
  void  SetRequired(){ cas |= 2;}
  void  Set(){cas=0;}
  GeometricalVertex() :cas(0), link(this) {};
  GeometricalVertex * The() { assert(link); return link;}// return a unique vertices
  int IsThe() const { return link == this;}  
  
inline void Set(const GeometricalVertex & rec,const Geometry & Gh ,const Geometry & GhNew);
  inline friend ostream& operator <<(ostream& f, const  GeometricalVertex & s) 
     { f << s.r << "," << s.cas << ".";return f; }
};

/////////////////////////////////////////////////////////////////////////////////////
class GeometricalEdge {
  public:
   GeometricalVertex * v[2];
   Int4 ref;
   Int4  CurveNumber;
  R2 tg[2]; // the 2 tangente 
  //   if tg[0] =0 => no continuite 
  GeometricalEdge * Adj [2]; 
  int SensAdj[2];
//  private:
  int flag ;
  public: 
  GeometricalEdge * link; // if   Cracked() or Equi()

// end of data 
  
  GeometricalVertex & operator[](int i){return *v[i];};
  const GeometricalVertex & operator[](int i) const { return *v[i];};
  GeometricalVertex * operator()(int i){return v[i];};  
 // inline void Set(const Geometry &,Int4,Geometry &);

  R2 F(Real8 theta) const ; // parametrization of the curve edge
  Real8 R1tg(Real8 theta,R2 &t) const ; // 1/radius of curvature + tangente
  int Cracked() const {return flag & 1;}
  int Dup() const { return flag & 32;}
  int Equi()const {return flag & 2;}
  int ReverseEqui()const {return flag & 128;}
  int TgA()const {return flag &4;}
  int TgB()const {return flag &8;}
  int Tg(int i) const{return i==0 ? TgA() : TgB();}
  int Mark()const {return flag &16;}
  int Required() { return flag & 64;}
  void SetCracked() { flag |= 1;}
  void SetDup()     { flag |= 32;} // not a real edge 
  void SetEqui()    { flag |= 2;}
  void SetTgA()     { flag|=4;}
  void SetTgB()     { flag|=8;}
  void SetMark()    { flag|=16;}
  void SetUnMark()  { flag &= 1007 /* 1023-16*/;}
  void SetRequired() { flag |= 64;}
  void SetReverseEqui() {flag |= 128;}
  
  inline void Set(const GeometricalEdge & rec,const Geometry & Th ,Geometry & ThNew);

#ifdef DRAWING 
  void Draw(Int4  =-1);
#endif
  
};
  
class Curve {public:
 GeometricalEdge * be,*ee; // begin et end edge
 int kb,ke;  //  begin vetex and end vertex
 Curve *next; // next curve equi to this
 bool master; // true => of equi curve point on this curve  
 inline void Set(const Curve & rec,const Geometry & Th ,Geometry & ThNew);
 Curve() : be(0),ee(0),kb(0),ke(0),next(0),master(true) {} 
 void Reverse() { Exchange(be,ee); Exchange(kb,ke);} //  revese the sens of the curse 
};
  
   

/////////////////////////////////////////////////////////////////////////////////////
class Triangle {
  friend class TriangleAdjacent;
  friend ostream& operator <<(ostream& f, const  Triangle & ta);


  private: // les arete sont opposes a un sommet
  Vertex * ns [3]; // 3 vertices if t is triangle, t[i] allowed by access function, (*t)[i] if pointer
  Triangle * at [3]; // nu triangle adjacent  
  Int1  aa[3];  // les nu des arete dans le triangles (mod 4)
  public: 
  Icoor2 det; // determinant du triangle (2 fois l aire des vertices entieres)
  union { 
    Triangle * link ;
    Int4 color;
  };
  void SetDet() {
    if(ns[0] && ns[1] && ns[2])    det = bamg::det(*ns[0],*ns[1],*ns[2]);
    else det = -1; }
  Triangle() {}
  Triangle(Triangles *Th,Int4 i,Int4 j,Int4 k);
  Triangle(Vertex *v0,Vertex *v1,Vertex *v2);
  inline void Set(const Triangle &,const Triangles &,Triangles &);
  inline int In(Vertex *v) const { return ns[0]==v || ns[1]==v || ns[2]==v ;}
  TriangleAdjacent FindBoundaryEdge(int ) const;

  void ReNumbering(Triangle *tb,Triangle *te, Int4 *renu) 
   {
    if (link  >=tb && link  <te) link  = tb + renu[link -tb];
    if (at[0] >=tb && at[0] <te) at[0] = tb + renu[at[0]-tb];
    if (at[1] >=tb && at[1] <te) at[1] = tb + renu[at[1]-tb];
    if (at[2] >=tb && at[2] <te) at[2] = tb + renu[at[2]-tb];    
   }
  void ReNumbering(Vertex *vb,Vertex *ve, Int4 *renu) 
   {
    if (ns[0] >=vb && ns[0] <ve) ns[0] = vb + renu[ns[0]-vb];
    if (ns[1] >=vb && ns[1] <ve) ns[1] = vb + renu[ns[1]-vb];
    if (ns[2] >=vb && ns[2] <ve) ns[2] = vb + renu[ns[2]-vb];    
   }


  const Vertex & operator[](int i) const {return *ns[i];};
  Vertex & operator[](int i)  {return *ns[i];};
  
  const Vertex  *  operator()(int i) const {return ns[i];};
  Vertex  * & operator()(int i)  {return ns[i];};
  
  TriangleAdjacent Adj(int  i) const  // triangle adjacent + arete 
  { return TriangleAdjacent(at[i],aa[i]&3);};

  Triangle * TriangleAdj(int  i) const 
  {return at[i&3];} // triangle adjacent + arete 
  Int1  NuEdgeTriangleAdj(int  i) const 
  {return aa[i&3]&3;} // Number of the  adjacent edge in adj tria  

  inline Real4 qualite() ;
  

  void SetAdjAdj(Int1 a) 
  { a &= 3;
   register  Triangle *tt=at[a];
   aa [a] &= 55; // remove MarkUnSwap
   register Int1 aatt = aa[a] & 3;
   if(tt){ 
     tt->at[aatt]=this;
     tt->aa[aatt]=a + (aa[a] & 60 ) ;}// Copy all the mark 
  }
  
  void SetAdj2(Int1 a,Triangle *t,Int1 aat)
    {  at[a]=t;aa[a]=aat;
    if(t) {t->at[aat]=this;t->aa[aat]=a;}
    }
    
  void SetTriangleContainingTheVertex()
   { 
     if (ns[0]) (ns[0]->t=this,ns[0]->vint=0);
     if (ns[1]) (ns[1]->t=this,ns[1]->vint=1);
     if (ns[2]) (ns[2]->t=this,ns[2]->vint=2);
   }
   
  int swap(Int2 a1,int=0);
  Int4  Optim(Int2 a,int =0);

  int  Locked(int a)const { return aa[a]&4;} 
  int  Hidden(int a)const { return aa[a]&16;} 
  int  Cracked(int a) const { return aa[a] & 32;}
  // for optimisation 
  int  GetAllflag(int a){return aa[a] & 1020;}
  void SetAllFlag(int a,int f){aa[a] = (aa[a] &3) + (1020 & f);}

  void SetHidden(int a){
   register Triangle * t = at[a];
    if(t) t->aa[aa[a] & 3] |=16;
    aa[a] |= 16;}
  void SetCracked(int a){
   register Triangle * t = at[a];
    if(t) t->aa[aa[a] & 3] |=32;
    aa[a] |= 32;}
  
  double   QualityQuad(int a,int option=1) const;
  Triangle * Quadrangle(Vertex * & v0,Vertex * & v1,Vertex * & v2,Vertex * & v3) const ;

  void SetLocked(int a){
    register Triangle * t = at[a];
    t->aa[aa[a] & 3] |=4;
    aa[a] |= 4;}

  void SetMarkUnSwap(int a){
    register Triangle * t = at[a];
    t->aa[aa[a] & 3] |=8;
    aa[a] |=8 ;}


  void SetUnMarkUnSwap(int a){ 
    register Triangle * t = at[a];
    t->aa[aa[a] & 3] &=55; // 23 + 32 
    aa[a] &=55 ;}
  

 
#ifdef DEBUG       
  void inline checka(Int1 a); 
  void inline check();
#endif

#ifdef DRAWING
  void  Draw(Int4 i=-1) const;
  int swapDRAW(Int2 a1);

#endif

};  // end of Triangle class 




class ListofIntersectionTriangles {
/////////////////////////////////////////////////////////////////////////////////////
class IntersectionTriangles {
public: 
  Triangle *t;
  Real8  bary[3];  // use if t != 0
  R2 x;
  Metric m;
  Real8 s;// abscisse curviline
  Real8 sp; // len of the previous seg in m
  Real8 sn;// len of the  next seg in m
};
/////////////////////////////////////////////////////////////////////////////////////
class SegInterpolation {
 public:
  GeometricalEdge * e;
  Real8 sBegin,sEnd; // abscisse of the seg on edge parameter
  Real8 lBegin,lEnd; // length abscisse  set in ListofIntersectionTriangles::Length
  int last;// last index  in ListofIntersectionTriangles for this Sub seg of edge
  R2 F(Real8 s){ 
    Real8 c01=lEnd-lBegin, c0=(lEnd-s)/c01, c1=(s-lBegin)/c01;
    assert(lBegin<= s && s <=lEnd);
    return e->F(sBegin*c0+sEnd*c1);}
};
 
 int MaxSize; // 
 int Size; //
 Real8 len; //
 int state;
 IntersectionTriangles * lIntTria;
 int NbSeg;
 int MaxNbSeg;
 SegInterpolation * lSegsI;
 public:
 IntersectionTriangles & operator[](int i) {return lIntTria[i];}
 operator int&() {return Size;}
 ListofIntersectionTriangles(int n=256,int m=16)
   :   MaxSize(n), Size(0), len(-1),state(-1),lIntTria(new IntersectionTriangles[n]) ,
    NbSeg(0), MaxNbSeg(m), lSegsI(new SegInterpolation[m])  
  { if (verbosity>9) 
    cout << "      construct ListofIntersectionTriangles"
			  << MaxSize << " " <<  MaxNbSeg<< endl;};
  ~ListofIntersectionTriangles(){
    if (lIntTria) delete [] lIntTria,lIntTria=0;
    if (lSegsI) delete [] lSegsI,lSegsI=0;} 
  void init(){state=0;len=0;Size=0;}
  
  int NewItem(Triangle * tt,Real8 d0,Real8 d1,Real8 d2);
  int NewItem(R2,const Metric & );
  void NewSubSeg(GeometricalEdge *e,Real8 s0,Real8 s1) 
    { 
      if (NbSeg>=MaxNbSeg) {
           int mneo= MaxNbSeg;
           MaxNbSeg *= 2;
          if (verbosity>3) 
	    cout <<" reshape lSegsI from " << mneo << " to " 
		 << MaxNbSeg <<endl;
           SegInterpolation * lEn =  new SegInterpolation[MaxNbSeg];
           assert(lSegsI && NbSeg < MaxNbSeg);
           for (int i=0;i< NbSeg;i++) 
              lEn[i] = lSegsI[MaxNbSeg]; // copy old to new            
           delete []  lSegsI; // remove old
           lSegsI = lEn;        
           }
       if (NbSeg) 
         lSegsI[NbSeg-1].last=Size;
       lSegsI[NbSeg].e=e;
       lSegsI[NbSeg].sBegin=s0;
       lSegsI[NbSeg].sEnd=s1;     
       NbSeg++;           
    }
    
//  void CopyMetric(int i,int j){ lIntTria[j].m=lIntTria[i].m;}
//  void CopyMetric(const Metric & mm,int j){ lIntTria[j].m=mm;}

  void ReShape() { 
    register int newsize = MaxSize*2;
    IntersectionTriangles * nw = new IntersectionTriangles[newsize];
    assert(nw);
    for (int i=0;i<MaxSize;i++) // recopy
      nw[i] = lIntTria[i];       
    if(verbosity>3)
      cout << " ListofIntersectionTriangles  ReShape MaxSize " 
	   << MaxSize << " -> " 
	 <<  newsize << endl;
    MaxSize = newsize; 
    delete [] lIntTria;// remove old
    lIntTria = nw; // copy pointer
  }
  
  void SplitEdge(const Triangles & ,const R2 &,const R2  &,int nbegin=0); 
  Real8 Length(); 
  Int4 NewPoints(Vertex *,Int4 & nbv,Int4 nbvx);
};


/////////////////////////////////////////////////////////////////////////////////////
class GeometricalSubDomain {
public:
  GeometricalEdge *edge;
  int sens; // -1 or 1
  Int4 ref;
  inline void Set(const GeometricalSubDomain &,const Geometry & ,const Geometry &);
  
};
/////////////////////////////////////////////////////////////////////////////////////
class SubDomain {
public:
  Triangle * head;
  Int4  ref;  
  int sens; // -1 or 1
  Edge * edge; // to  geometrical 	
  inline void Set(const Triangles &,Int4,Triangles &);
  	 
};
/////////////////////////////////////////////////////////////////////////////////////
class VertexOnGeom {  public:

  Vertex * mv;
  Real8 abscisse;  
  union{ 
    GeometricalVertex * gv; // if abscisse <0; 
    GeometricalEdge * ge;  // if abscisse in [0..1]
  };
  inline void Set(const VertexOnGeom&,const Triangles &,Triangles &);  
  int OnGeomVertex()const {return this? abscisse <0 :0;}
  int OnGeomEdge() const {return this? abscisse >=0 :0;}
  VertexOnGeom(): mv(0),abscisse(0){gv=0;} 
  VertexOnGeom(Vertex & m,GeometricalVertex &g) : mv(&m),abscisse(-1){gv=&g;}
   //  cout << "        mv = " <<mv << " gv = "  << gv << endl;} 
  VertexOnGeom(Vertex & m,GeometricalEdge &g,Real8 s) : mv(&m),abscisse(s){ge=&g;}
    //cout << &g << " "  << ge << endl;} 
  operator Vertex * () const  {return mv;}
  operator GeometricalVertex * () const  {return gv;}
  operator GeometricalEdge * () const  {return ge;}
//  operator Real8 & () {return abscisse;}
  operator const Real8 & () const {return abscisse;}
  int IsRequiredVertex(){ return this? (( abscisse<0 ? (gv?gv->Required():0):0 )) : 0;}
  void SetOn(){mv->on=this;mv->vint=IsVertexOnGeom;}
  friend ostream& operator <<(ostream& f, const  VertexOnGeom & vog){
    f << vog.abscisse << " " << vog.mv << " " << vog.gv << " ; ";
    if (vog.abscisse < 0) f << *vog.gv << " ;; " ;
    //    else f << *vog.ge << " ;; " ;
    return f;}
  inline void Set(const Triangles &,Int4,Triangles &);
    
};
/////////////////////////////////////////////////////////////////////////////////////
class VertexOnVertex {public:
  Vertex * v, *bv;
  VertexOnVertex(Vertex * w,Vertex *bw) :v(w),bv(bw){}
  VertexOnVertex() {};
  inline void Set(const Triangles &,Int4,Triangles &);
  void SetOnBTh(){v->onbv=bv;v->vint=IsVertexOnVertex;}
};
/////////////////////////////////////////////////////////////////////////////////////
class VertexOnEdge {public:
  Vertex * v;
  Edge * be;
  Real8 abcisse;
  VertexOnEdge( Vertex * w, Edge *bw,Real8 s) :v(w),be(bw),abcisse(s) {}
  VertexOnEdge(){}
  inline void Set(const Triangles &,Int4,Triangles &);  
  void SetOnBTh(){v->onbe=this;v->vint=IsVertexOnEdge;}  
  Vertex & operator[](int i) const { return (*be)[i];}
  operator Real8 () const { return abcisse;}
  operator  Vertex *  () const { return v;}  
  operator  Edge *  () const { return be;}  
};

 inline TriangleAdjacent FindTriangleAdjacent(Edge &E);
 inline Vertex * TheVertex(Vertex * a); // for remove crak in mesh 
/////////////////////////////////////////////////////////////////////////////////////
 
class CrackedEdge { // a small class to store on crack an uncrack information 
  friend class Triangles;
  friend ostream& operator <<(ostream& f, const   Triangles & Th) ;  
 class CrackedTriangle {
  friend class Triangles;
  friend class CrackedEdge;
  friend ostream& operator <<(ostream& f, const   Triangles & Th) ;  
  Triangle * t; // edge of triangle t
  int i; //  edge number of in triangle
  Edge *edge; // the  2 edge 
  Vertex *New[2]; // new vertex number 
  CrackedTriangle() : t(0),i(0),edge(0) { New[0]=New[1]=0;} 
  CrackedTriangle(Edge * a) : t(0),i(0),edge(a) { New[0]=New[1]=0;} 
  void Crack(){ 
    Triangle & T(*t); 
    int i0=VerticesOfTriangularEdge[i][0];
    int i1=VerticesOfTriangularEdge[i][0];
    assert(New[0] && New[1]);
    T(i0) = New[0];
    T(i1) = New[1];}    
  void UnCrack(){ 
    Triangle & T(*t); 
    int i0=VerticesOfTriangularEdge[i][0];
    int i1=VerticesOfTriangularEdge[i][0];
    assert(New[0] && New[1]);
    T(i0) = TheVertex(T(i0));
    T(i1) = TheVertex(T(i1));} 
  void Set() {
     TriangleAdjacent ta ( FindTriangleAdjacent(*edge));
     t = ta;
     i = ta;
     
     New[0] = ta.EdgeVertex(0);
     New[1] = ta.EdgeVertex(1);
     // warning the ref 
     
     }    
     
  }; // end of class CrackedTriangle
 public:  
  CrackedTriangle a,b; 
  CrackedEdge() :a(),b() {}
  CrackedEdge(Edge * start, Int4  i,Int4 j) : a(start+i),b(start+j) {};
  CrackedEdge(Edge * e0, Edge * e1 ) : a(e0),b(e1) {};

  void Crack() { a.Crack(); b.Crack();}
  void UnCrack() { a.UnCrack(); b.UnCrack();}
  void Set() { a.Set(), b.Set();}
};

/////////////////////////////////////////////////////////////////////////////////////

class Triangles { 
public:

  enum TypeFileMesh {
    AutoMesh=0,BDMesh=1,NOPOMesh=2,amMesh=3,am_fmtMesh=4,amdbaMesh=5,
    ftqMesh=6,mshMesh=7};

  int static counter; // to kown the number of mesh in memory 
  int OnDisk;       // true if on disk 
  Geometry & Gh;   // Geometry
  Triangles & BTh; // Background Mesh Bth==*this =>no  background 
  
  Int4 NbRef; // counter of ref on the this class if 0 we can delete
  Int4 nbvx,nbtx;  // nombre max  de sommets , de  triangles
  
  Int4 nt,nbv,nbt,nbiv,nbe; // nb of legal triangles, nb of vertex, of triangles, 
  // of initial vertices, of edges with reference,
  Int4 NbOfQuad; // nb of quadrangle 

  Int4 NbSubDomains; // 
  Int4 NbOutT; // Nb of oudeside triangle
  Int4 NbOfTriangleSearchFind;
  Int4 NbOfSwapTriangle;
  char * name, *identity;
  Vertex * vertices;   // data of vertices des sommets
  
  Int4 NbVerticesOnGeomVertex;
  VertexOnGeom * VerticesOnGeomVertex;
  
  Int4 NbVerticesOnGeomEdge;
  VertexOnGeom * VerticesOnGeomEdge;

  Int4 NbVertexOnBThVertex;
  VertexOnVertex *VertexOnBThVertex;

  Int4 NbVertexOnBThEdge;
  VertexOnEdge *VertexOnBThEdge;

  
  Int4 NbCrackedVertices;
  

  Int4 NbCrackedEdges;
  CrackedEdge *CrackedEdges;
  
  
  R2 pmin,pmax; // extrema
  Real8 coefIcoor;  // coef to integer Icoor1;

  Triangle * triangles;
  Edge * edges; 

  QuadTree *quadtree;
  Vertex ** ordre;
  SubDomain * subdomains;
  ListofIntersectionTriangles  lIntTria;
// end of variable
  
  Triangles(Int4 i);//:BTh(*this),Gh(*new Geometry()){PreInit(i);}
  
  ~Triangles(); 
  Triangles(const char * ,Real8=-1) ;
  
  Triangles(Int4 nbvx,Triangles & BT,int keepBackVertices=1)
         :Gh(BT.Gh),BTh(BT) {
	     try {GeomToTriangles1(nbvx,keepBackVertices);}
	      catch(...) { this->~Triangles(); throw; } }
  
  Triangles(Int4 nbvx,Geometry & G)
         :Gh(G),BTh(*this){
	     try { GeomToTriangles0(nbvx);}
	     catch(...) { this->~Triangles(); throw; } }
  Triangles(Triangles &,Geometry * pGh=0,Triangles* pBTh=0,Int4 nbvxx=0 ); // COPY OPERATEUR
  //  Triangles(Triangles &){ cerr << " BUG call copy opretor of Triangles" << endl;MeshError(111);}
  Triangles(const Triangles &,const int *flag,const int *bb); // truncature
  

  void SetIntCoor(const char * from =0);

 // void  RandomInit();
 // void  CubeInit(int ,int);
  
  Real8 MinimalHmin() {return 2.0/coefIcoor;}
  Real8 MaximalHmax() {return Max(pmax.x-pmin.x,pmax.y-pmin.y);}
  const Vertex & operator[]  (Int4 i) const { return vertices[i];};
  Vertex & operator[](Int4 i) { return vertices[i];};
  const Triangle & operator()  (Int4 i) const { return triangles[i];};
  Triangle & operator()(Int4 i) { return triangles[i];};
  I2 toI2(const R2 & P) const {
          return  I2( (Icoor1) (coefIcoor*(P.x-pmin.x))
	                 ,(Icoor1) (coefIcoor*(P.y-pmin.y)) );}
  R2 toR2(const I2 & P) const {
          return  R2( (double) P.x/coefIcoor+pmin.x, (double) P.y/coefIcoor+pmin.y);}
  void Add( Vertex & s,Triangle * t,Icoor2 *  =0) ;
  void Insert();
  //  void InsertOld();
  void ForceBoundary();
  void Heap();
  void FindSubDomain(int );
  Int4  ConsRefTriangle(Int4 *) const;
  void ShowHistogram() const;
  void  ShowRegulaty() const; // Add FH avril 2007 
//  void ConsLinkTriangle();

  void ReMakeTriangleContainingTheVertex();
  void UnMarkUnSwapTriangle();
  void SmoothMetric(Real8 raisonmax) ;
  void BoundAnisotropy(Real8 anisomax,double hminaniso= 1e-100) ;
  void MaxSubDivision(Real8 maxsubdiv);
  void WriteMetric(ostream &,int iso) ;
  Edge** MakeGeometricalEdgeToEdge();
  void  SetVertexFieldOn();  
  void  SetVertexFieldOnBTh();
  Int4 SplitInternalEdgeWithBorderVertices();
  void MakeQuadrangles(double costheta);
  int SplitElement(int choice);
  void MakeQuadTree();
  void NewPoints( Triangles &,int KeepBackVertex =1 );
  Int4 InsertNewPoints(Int4 nbvold,Int4 & NbTSwap) ; 
  void NewPointsOld( Triangles & );
  void NewPoints(int KeepBackVertex=1){ NewPoints(*this,KeepBackVertex);}
  void ReNumberingTheTriangleBySubDomain(bool justcompress=false);
  void ReNumberingVertex(Int4 * renu);
  void SmoothingVertex(int =3,Real8=0.3);
  Metric MetricAt (const R2 &) const;
  GeometricalEdge * ProjectOnCurve( Edge & AB, Vertex &  A, Vertex & B,Real8 theta,
                      Vertex & R,VertexOnEdge & BR,VertexOnGeom & GR);
   
  
  void WriteElements(ostream& f,Int4 * reft ,Int4 nbInT) const;

  
  Int4 Number(const Triangle & t) const  { return &t - triangles;}
  Int4 Number(const Triangle * t) const  { return t - triangles;}
  Int4 Number(const Vertex & t) const  { return &t - vertices;}
  Int4 Number(const Vertex * t) const  { return t - vertices;}
  Int4 Number(const Edge & t) const  { return &t - edges;}
  Int4 Number(const Edge * t) const  { return t - edges;}
  Int4 Number2(const Triangle * t) const  {
    //   if(t>= triangles && t < triangles + nbt )
      return t - triangles;
      //  else  return t - OutSidesTriangles;
  }
  
  Vertex * NearestVertex(Icoor1 i,Icoor1 j) ;
  Triangle * FindTriangleContening(const I2 & ,Icoor2 [3],Triangle *tstart=0) const;
  void Write(const char * filename,const TypeFileMesh type = AutoMesh);
  void Write_am_fmt(ostream &) const ;
  void Write_am(ostream &) const ;
  void Write_ftq(ostream &) const ;
  void Write_nopo(ostream &) const ;
  void Write_msh(ostream &) const ;
  void Write_amdba(ostream &) const ;

  void Read(MeshIstream &,int version,Real8 cutoffradian);
  void Read_am_fmt(MeshIstream &);
  void Read_amdba(MeshIstream &);
  void Read_am(MeshIstream &);
  void Read_nopo(MeshIstream &);
  void Read_ftq(MeshIstream &);
  void Read_msh(MeshIstream &);

  void ReadMetric(const char * fmetrix,const Real8 hmin,const Real8 hmax,const Real8 coef);
  void IntersectConsMetric(const double * s,const Int4 nbsol,const int * typsols,
			   const  Real8 hmin,const Real8 hmax, const Real8 coef,
			   const Real8  anisomax,const Real8 CutOff=1.e-4,const int NbJacobi=1,
			   const int DoNormalisation=1,
			   const double power=1.0,
			   const int choise=0);
  void IntersectGeomMetric(const Real8 err,const int iso);
  
  
  int  isCracked() const {return NbCrackedVertices != 0;}
  int  Crack();
  int UnCrack();
  
#ifdef DEBUG
  void inline Check(); 
#endif
#ifdef DRAWING
  void  Draw() const ;
  void  InitDraw() const ;
  void   inquire()  ;
#endif
 friend ostream& operator <<(ostream& f,  const  Triangles & Th); 
  void  Write(const char * filename);
  void ConsGeometry(Real8 =-1.0,int *equiedges=0); // construct a geometry if no geo 
  void FillHoleInMesh() ;
  int CrackMesh();
 private:
  void GeomToTriangles1(Int4 nbvx,int KeepBackVertices=1);// the  real constructor mesh adaption
  void GeomToTriangles0(Int4 nbvx);// the  real constructor mesh generator
  void PreInit(Int4,char * =0 );
  //
  void Write_nop5(OFortranUnFormattedFile * f,
			     Int4 &lnop5,Int4 &nef,Int4 &lgpdn,Int4 ndsr) const ;

  
}; // End Class Triangles
/////////////////////////////////////////////////////////////////////////////////////
class Geometry { 
public:
  int OnDisk; 
  Int4 NbRef; // counter of ref on the this class if 0 we can delete

  char *name;
  Int4 nbvx,nbtx; // nombre max  de sommets , de  Geometry
  Int4 nbv,nbt,nbiv,nbe; // nombre de sommets, de Geometry, de sommets initiaux,
  Int4 NbSubDomains; // 
  Int4 NbEquiEdges;
  Int4 NbCrackedEdges; 
//  Int4 nbtf;//  de triangle frontiere
  Int4 NbOfCurves;
  int empty(){return (nbv ==0) && (nbt==0) && (nbe==0) && (NbSubDomains==0); }
  GeometricalVertex * vertices;   // data of vertices des sommets 
  Triangle * triangles; 
  GeometricalEdge * edges;
  QuadTree *quadtree;
  GeometricalSubDomain *subdomains;
  Curve *curves;
  ~Geometry(); 
  Geometry(const Geometry & Gh); //Copy  Operator 
  Geometry(int nbg,const Geometry ** ag); // intersection operator 
  
  R2 pmin,pmax; // extrema
  Real8 coefIcoor;  // coef to integer Icoor1;
  Real8 MaximalAngleOfCorner;
  
//  end of data 

  
  I2 toI2(const R2 & P) const {
          return  I2( (Icoor1) (coefIcoor*(P.x-pmin.x))
	                 ,(Icoor1) (coefIcoor*(P.y-pmin.y)) );}
  
  Real8 MinimalHmin() {return 2.0/coefIcoor;}
  Real8 MaximalHmax() {return Max(pmax.x-pmin.x,pmax.y-pmin.y);}
  void ReadGeometry(const char * ) ;
  void ReadGeometry(MeshIstream & ,const char *)  ;

  void EmptyGeometry();
  Geometry() {EmptyGeometry();}// empty Geometry
  void AfterRead();
  Geometry(const char * filename) {EmptyGeometry();OnDisk=1;ReadGeometry(filename);AfterRead();}

  void ReadMetric(const char *,Real8 hmin,Real8 hmax,Real8 coef);
  const GeometricalVertex & operator[]  (Int4 i) const { return vertices[i];};
  GeometricalVertex & operator[](Int4 i) { return vertices[i];};
  const  GeometricalEdge & operator()  (Int4 i) const { return edges[i];};
  GeometricalEdge & operator()(Int4 i) { return edges[i];}; 
  Int4 Number(const GeometricalVertex & t) const  { return &t - vertices;}
  Int4 Number(const GeometricalVertex * t) const  { return t - vertices;}
  Int4 Number(const GeometricalEdge & t) const  { return &t - edges;}
  Int4 Number(const GeometricalEdge * t) const  { return t - edges;}
  Int4 Number(const Curve * c) const  { return c - curves;}
  
  void UnMarkEdges() {
    for (Int4 i=0;i<nbe;i++) edges[i].SetUnMark();}

 GeometricalEdge *  ProjectOnCurve(const Edge & ,Real8,Vertex &,VertexOnGeom &) const ;
  GeometricalEdge *  Contening(const R2 P,  GeometricalEdge * start) const;
 friend ostream& operator <<(ostream& f, const   Geometry & Gh); 
 void Write(const char * filename);
 
#ifdef DEBUG
  void inline Check();   
#endif
#ifdef DRAWING
  void  Draw() const ;
  void  InitDraw() const ;
#endif
  
}; // End Class Geometry

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
///////////////////           END CLASS          ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

inline Triangles::Triangles(Int4 i) :Gh(*new Geometry()),BTh(*this){PreInit(i);}

extern Triangles * CurrentTh;

TriangleAdjacent CloseBoundaryEdge(I2 ,Triangle *, double &,double &) ;
TriangleAdjacent CloseBoundaryEdgeV2(I2 A,Triangle *t, double &a,double &b);

Int4 FindTriangle(Triangles &Th, Real8 x, Real8 y, double* a,int & inside);



inline Triangle *    Triangle::Quadrangle(Vertex * & v0,Vertex * & v1,Vertex * & v2,Vertex * & v3) const
{
// return the other triangle of the quad if a quad or 0 if not a quat
  Triangle * t =0;
  if (link) {
    int a=-1;
    if (aa[0] & 16 ) a=0;
    if (aa[1] & 16 ) a=1;
    if (aa[2] & 16 ) a=2;
    if (a>=0) {
      t = at[a];
      //  if (t-this<0) return 0;
      v2 = ns[VerticesOfTriangularEdge[a][0]];
      v0 = ns[VerticesOfTriangularEdge[a][1]];
      v1 = ns[OppositeEdge[a]];
      v3 = t->ns[OppositeEdge[aa[a]&3]];
    }
  }
  return t;
}

inline   double   Triangle::QualityQuad(int a,int option) const
{ // first do the logique part 
  double q;
  if (!link || aa[a] &4)
    q=  -1;
  else {
    Triangle * t = at[a];
    if (t-this<0) q=  -1;// because we do 2 times 
    else if (!t->link ) q=  -1;
    else if (aa[0] & 16 || aa[1] & 16  || aa[2] & 16 || t->aa[0] & 16 || t->aa[1] & 16 || t->aa[2] & 16 )
      q= -1;
    else if(option) 
      { 
	const Vertex & v2 = *ns[VerticesOfTriangularEdge[a][0]];
	const Vertex & v0 = *ns[VerticesOfTriangularEdge[a][1]];
	const Vertex & v1 = *ns[OppositeEdge[a]];
	const Vertex & v3 = * t->ns[OppositeEdge[aa[a]&3]];
	q =  QuadQuality(v0,v1,v2,v3); // do the float part
      }
    else q= 1;
  }
  return  q;
}


inline void Vertex::Set(const Vertex & rec,const Triangles & ,Triangles & )
 { 
   *this  = rec;
 }
inline void GeometricalVertex::Set(const GeometricalVertex & rec,const Geometry & ,const Geometry & )
 { 
   *this  = rec;
 }
inline void Edge::Set(const Triangles & Th ,Int4 i,Triangles & ThNew)
 { 
   *this = Th.edges[i];
   v[0] = ThNew.vertices + Th.Number(v[0]);    
   v[1] = ThNew.vertices + Th.Number(v[1]);
   if (on) 
     on =  ThNew.Gh.edges+Th.Gh.Number(on);
   if (adj[0]) adj[0] =   ThNew.edges +   Th.Number(adj[0]);
   if (adj[1]) adj[1] =   ThNew.edges +   Th.Number(adj[1]);

 }
inline void GeometricalEdge::Set(const GeometricalEdge & rec,const Geometry & Gh ,Geometry & GhNew)
 { 
   *this = rec;
   v[0] = GhNew.vertices + Gh.Number(v[0]);    
   v[1] = GhNew.vertices + Gh.Number(v[1]); 
   if (Adj[0]) Adj[0] =  GhNew.edges + Gh.Number(Adj[0]);     
   if (Adj[1]) Adj[1] =  GhNew.edges + Gh.Number(Adj[1]);     
 }
 
inline void Curve::Set(const Curve & rec,const Geometry & Gh ,Geometry & GhNew)
{
  *this = rec;
   be = GhNew.edges + Gh.Number(be);    
   ee = GhNew.edges + Gh.Number(ee); 
   if(next) next= GhNew.curves + Gh.Number(next); 
}

inline void Triangle::Set(const Triangle & rec,const Triangles & Th ,Triangles & ThNew)
 { 
   *this = rec;
   if ( ns[0] ) ns[0] = ThNew.vertices +  Th.Number(ns[0]);
   if ( ns[1] ) ns[1] = ThNew.vertices +  Th.Number(ns[1]);
   if ( ns[2] ) ns[2] = ThNew.vertices +  Th.Number(ns[2]);
   if(at[0]) at[0] =  ThNew.triangles + Th.Number(at[0]);
   if(at[1]) at[1] =  ThNew.triangles + Th.Number(at[1]);
   if(at[2]) at[2] =  ThNew.triangles + Th.Number(at[2]);
   if (link  >= Th.triangles && link  < Th.triangles + Th.nbt)
     link = ThNew.triangles + Th.Number(link);
 }
inline void VertexOnVertex::Set(const Triangles & Th ,Int4 i,Triangles & ThNew)
{ 
  *this = Th.VertexOnBThVertex[i];  
  v = ThNew.vertices + Th.Number(v);

}
inline void SubDomain::Set(const Triangles & Th ,Int4 i,Triangles & ThNew)
{
  *this = Th.subdomains[i];
  assert( head - Th.triangles >=0 && head - Th.triangles < Th.nbt);
  head = ThNew.triangles + Th.Number(head) ; 
  assert(edge - Th.edges >=0 && edge - Th.edges < Th.nbe); 
  edge = ThNew.edges+ Th.Number(edge);
}
 inline void GeometricalSubDomain::Set(const GeometricalSubDomain & rec,const Geometry & Gh ,const Geometry & GhNew)
{
   *this = rec;
   edge = Gh.Number(edge) + GhNew.edges;
}
inline void VertexOnEdge::Set(const Triangles & Th ,Int4 i,Triangles & ThNew)
{
  *this = Th.VertexOnBThEdge[i];  
  v = ThNew.vertices + Th.Number(v);
}

inline void VertexOnGeom::Set(const VertexOnGeom & rec,const Triangles & Th ,Triangles & ThNew)
{
  *this = rec;  
  mv = ThNew.vertices + Th.Number(mv);
  if (gv)
    if (abscisse < 0 )
      gv = ThNew.Gh.vertices + Th.Gh.Number(gv);
    else
      ge = ThNew.Gh.edges + Th.Gh.Number(ge);
  
}
inline Real8 Edge::MetricLength() const
  { 
    return LengthInterpole(v[0]->m,v[1]->m,v[1]->r - v[0]->r) ;
  }

inline  void  Triangles::ReMakeTriangleContainingTheVertex()
 {
  register Int4 i;
  for ( i=0;i<nbv;i++) 
    {
	vertices[i].vint = 0;
	vertices[i].t=0;
    }
  for ( i=0;i<nbt;i++) 
    triangles[i].SetTriangleContainingTheVertex();
 }

inline  void  Triangles::UnMarkUnSwapTriangle()
 {
  register Int4 i;
  for ( i=0;i<nbt;i++) 
    for(int  j=0;j<3;j++)
      triangles[i].SetUnMarkUnSwap(j);
 }

inline  void   Triangles::SetVertexFieldOn()
  {
    for (Int4 i=0;i<nbv;i++) 
       vertices[i].on=0;
    for (Int4 j=0;j<NbVerticesOnGeomVertex;j++ ) 
       VerticesOnGeomVertex[j].SetOn();
    for (Int4 k=0;k<NbVerticesOnGeomEdge;k++ ) 
       VerticesOnGeomEdge[k].SetOn();
    }	       
inline  void   Triangles::SetVertexFieldOnBTh()
  {
    for (Int4 i=0;i<nbv;i++) 
       vertices[i].on=0;
    for (Int4 j=0;j<NbVertexOnBThVertex;j++ ) 
       VertexOnBThVertex[j].SetOnBTh();
    for (Int4 k=0;k<NbVertexOnBThEdge;k++ ) 
       VertexOnBThEdge[k].SetOnBTh();
       
    }	       

inline  void  TriangleAdjacent::SetAdj2(const TriangleAdjacent & ta, int l  )
{ // set du triangle adjacent 
  if(t) {
    t->at[a]=ta.t;
    t->aa[a]=ta.a|l;}
  if(ta.t) {
    ta.t->at[ta.a] = t ;
    ta.t->aa[ta.a] = a| l ;
  }
}


inline int  TriangleAdjacent::Locked() const
{ return t->aa[a] &4;}
inline int  TriangleAdjacent::Cracked() const
{ return t->aa[a] &32;}
inline int  TriangleAdjacent::GetAllFlag_UnSwap() const
{ return t->aa[a] & 1012;} // take all flag except MarkUnSwap

inline int  TriangleAdjacent::MarkUnSwap() const
{ return t->aa[a] &8;}

inline void  TriangleAdjacent::SetLock(){ t->SetLocked(a);}

inline void  TriangleAdjacent::SetCracked() { t->SetCracked(a);}

inline  TriangleAdjacent TriangleAdjacent::Adj() const
{ return  t->Adj(a);}

inline Vertex  * TriangleAdjacent::EdgeVertex(const int & i) const
 {return t->ns[VerticesOfTriangularEdge[a][i]]; }
inline Vertex  * TriangleAdjacent::OppositeVertex() const
{return t->ns[bamg::OppositeVertex[a]]; }
inline Icoor2 &  TriangleAdjacent::det() const
{ return t->det;}
inline  TriangleAdjacent Adj(const TriangleAdjacent & a)
{ return  a.Adj();}

inline TriangleAdjacent Next(const TriangleAdjacent & ta) 
{ return TriangleAdjacent(ta.t,NextEdge[ta.a]);}

inline TriangleAdjacent Previous(const TriangleAdjacent & ta) 
{ return TriangleAdjacent(ta.t,PreviousEdge[ta.a]);}
 
inline void Adj(GeometricalEdge * & on,int &i) 
  {int j=i;i=on->SensAdj[i];on=on->Adj[j];}
  
inline Real4 qualite(const Vertex &va,const Vertex &vb,const Vertex &vc)
{
  Real4 ret; 
  I2 ia=va,ib=vb,ic=vc;
  I2 ab=ib-ia,bc=ic-ib,ac=ic-ia;
  Icoor2 deta=Det(ab,ac);
  if (deta <=0) ret = -1;
   else {
     Real8 a = sqrt((Real8) (ac,ac)),
       b = sqrt((Real8) (bc,bc)),
       c = sqrt((Real8) (ab,ab)),
       p = a+b+c;
     Real8 h= Max(Max(a,b),c),ro=deta/p;
   ret = ro/h;}
  return ret;
}


inline  Triangle::Triangle(Triangles *Th,Int4 i,Int4 j,Int4 k) {
  Vertex *v=Th->vertices;
  Int4 nbv = Th->nbv;
  assert(i >=0 && j >=0 && k >=0);
  assert(i < nbv && j < nbv && k < nbv);
  ns[0]=v+i;
  ns[1]=v+j;
  ns[2]=v+k;
  at[0]=at[1]=at[2]=0;
  aa[0]=aa[1]=aa[2]=0;
  det=0;
}

inline  Triangle::Triangle(Vertex *v0,Vertex *v1,Vertex *v2){
  ns[0]=v0;
  ns[1]=v1;
  ns[2]=v2;
  at[0]=at[1]=at[2]=0;
  aa[0]=aa[1]=aa[2]=0;
  if (v0) det=0;
  else {
    det=-1;
    link=NULL;};  
}

inline    Real4 Triangle::qualite()
{
  return det < 0 ? -1 :  bamg::qualite(*ns[0],*ns[1],*ns[2]);
}

Int4 inline  Vertex::Optim(int i,int koption)
{ 
  Int4 ret=0;
  if ( t && (vint >= 0 ) && (vint <3) )
    {
      ret = t->Optim(vint,koption);
      if(!i) 
	{
	  t =0; // for no future optime 
	  vint= 0; }
    }
  return ret;
}

Icoor2 inline det(const Vertex & a,const Vertex & b,const Vertex & c)
{
  register  Icoor2 bax = b.i.x - a.i.x ,bay = b.i.y - a.i.y; 
  register  Icoor2 cax = c.i.x - a.i.x ,cay = c.i.y - a.i.y; 
  return  bax*cay - bay*cax;}


void  swap(Triangle *t1,Int1 a1,
	   Triangle *t2,Int1 a2,
	   Vertex *s1,Vertex *s2,Icoor2 det1,Icoor2 det2);



int inline TriangleAdjacent::swap()
{ return  t->swap(a);}



int SwapForForcingEdge(Vertex   *  & pva ,Vertex  * &   pvb ,
		       TriangleAdjacent & tt1,Icoor2 & dets1,
		       Icoor2 & detsa,Icoor2 & detsb, int & nbswap);

int ForceEdge(Vertex &a, Vertex & b,TriangleAdjacent & taret) ;

// inline bofbof   FH 
inline  TriangleAdjacent FindTriangleAdjacent(Edge &E)
  {
    Vertex * a = E.v[0];
    Vertex * b = E.v[1];
    
    Triangle * t = a->t;
    int i = a->vint;
    TriangleAdjacent ta(t,EdgesVertexTriangle[i][0]); // Previous edge
    assert(t && i>=0 && i < 3);
    assert( a == (*t)(i));
    int k=0;
    do { // turn around vertex in direct sens (trigo)
      k++;assert(k< 20000);
      //  in no crack => ta.EdgeVertex(1) == a otherwise ??? 
      if (ta.EdgeVertex(1) ==  a && ta.EdgeVertex(0) ==  b) return ta; // find 
      ta = ta.Adj();
      if (ta.EdgeVertex(0) ==  a && ta.EdgeVertex(1) ==  b) return ta; // find 
      --ta;
      } while (t != (Triangle *)ta);
    assert(0);
    return TriangleAdjacent(0,0);// error 
  }
  
inline Vertex * TheVertex(Vertex * a) // give a unique vertex with smallest number
{ // in case on crack in mesh 
    Vertex * r(a), *rr;
    Triangle * t = a->t;
    int i = a->vint;
    TriangleAdjacent ta(t,EdgesVertexTriangle[i][0]); // Previous edge
    assert(t && i>=0 && i < 3);
    assert( a == (*t)(i));
    int k=0;
    do { // turn around vertex in direct sens (trigo)
      k++;assert(k< 20000);
      //  in no crack => ta.EdgeVertex(1) == a
      if ((rr=ta.EdgeVertex(0)) < r) r = rr;
      ta = ta.Adj();
      if ((rr=ta.EdgeVertex(1)) < r) r =rr;
      --ta;
     } while (t != (Triangle*) ta);  
    return r;
}

inline double CPUtime(){
#ifdef SYSTIMES
  struct tms buf;
  if (times(&buf)!=-1)
    return ((double)buf.tms_utime+(double)buf.tms_stime)/(long) sysconf(_SC_CLK_TCK);
  else
#endif
    return ((double) clock())/CLOCKS_PER_SEC;
}

#ifdef DEBUG
void inline Triangle::checka(Int1 a) {
  // verif de la coherence des adjacences de l arete a
  a = a%4;
  assert(a < 3 && a >= 0 );
  Triangle *t1=this,*t2=at[a];
  Int2 a1=a,a2=aa[a]%4;
  
  assert(a2 < 3 && a2 >= 0 );
  if (t2 && ( ((*t1).ns[VerticesOfTriangularEdge[a1][0]] != (*t2).ns[VerticesOfTriangularEdge[a2][1]])
	      || ((*t1).ns[VerticesOfTriangularEdge[a1][1]] != (*t2).ns[VerticesOfTriangularEdge[a2][0]])))
    {
      if (CurrentTh) cerr << " In Triangles beetween Triangle " << CurrentTh->Number(t1) << " and " 
		     <<  CurrentTh->Number(t2) <<  endl;
      cerr << "---- t1="<< t1 << " " << a1 <<",  t2="<< t2 << " " << a2 << endl;
      cerr <<"t1="<< t1 << " " << a1 << " " << t1->ns[VerticesOfTriangularEdge[a1][0]] 
	   << " " << t1->ns[VerticesOfTriangularEdge[a1][1]] <<endl;
      if (CurrentTh)
	cerr <<"t1="<< t1 << " " << a1 << " " << CurrentTh->Number(t1->ns[VerticesOfTriangularEdge[a1][0]])
	   << " " << CurrentTh->Number(t1->ns[VerticesOfTriangularEdge[a1][1]]) <<endl;
      if (t2) cerr <<"t2="<< t2 << " " << a2 << " " 
		   <<  t2->ns[VerticesOfTriangularEdge[a2][0]] 
		   << " " << t2->ns[VerticesOfTriangularEdge[a2][1]] <<endl;
      if (t2 &&CurrentTh)
	cerr <<"t2="<< t2 << " " << a2 << " " 
	     <<  CurrentTh->Number(t2->ns[VerticesOfTriangularEdge[a2][0]])
	     << " " << CurrentTh->Number(t2->ns[VerticesOfTriangularEdge[a2][1]]) <<endl;
      assert(0); 
    } 
  if (t2)   assert(t1->aa[a1]/4 == t2->aa[a2]/4); // lock compatibite
}


void inline  Triangle::check() {
  Icoor2 det2=0;
  //  cout << " check " << this << endl;
  int  infv=ns[0] ?  ((  ns[1] ? ( ns[2] ? -1 : 2) : 1  )) : 0;
  if (det<0) {
      if (infv<0 )
	{  if (CurrentTh) cerr << " In Triangles " << CurrentTh->Number(this) << endl;
	cerr << " det = " <<  det << " and " << infv << endl;
	MeshError(5);
	}}
  else  if (infv>=0 )
 	{  if (CurrentTh) cerr << " In Triangles " << CurrentTh->Number(this) << endl;
	cerr << " det = " << det << " and " << infv << endl;
	MeshError(5);
	}  
  
  if (det >=0) 
    if( det != (det2=bamg::det(*ns[0],*ns[1],*ns[2])))
	{ // penthickness(4);Draw();
	  if (CurrentTh) cerr << " In Triangles" << CurrentTh->Number(this) 
			 << endl;
      cerr << *ns[0] << *ns[1] << " " << *ns[2]  << " " << endl;
	   cerr << " Bug in triangle " << this 
		<< ":" << det << " !=  " << det2 << endl;
	   MeshError(5);
	}
  checka(0);
  checka(1);
  checka(2);
//  if (ns[0]) assert( ns[0] - Meshbegin  >= 0 );
//  if (ns[0]) assert( Meshend  - ns[0] >= 0 );
//  if (ns[1]) assert( ns[1] - Meshbegin  >= 0 );
//  if (ns[1]) assert( Meshend  - ns[1] >= 0 );
//  if (ns[2]) assert( ns[2] - Meshbegin  >= 0 );
//  if (ns[2]) assert( Meshend  - ns[2] >= 0 );
  assert(ns[0] != ns[2]);
  assert(ns[1] != ns[2]);
  assert(ns[0] != ns[1]);
}


#endif




#ifdef DRAWING 
extern Real4 xGrafCoef,yGrafCoef,xGrafOffSet,yGrafOffSet; // R2 -> I2 transform
extern R2 Gpmin,Gpmax;
//extern Real8 Gh;
// cf routine ILineTo IMoveto

extern void  IMoveTo(long i,long j);
extern void  ILineTo(long i,long j);
extern char Getxyc(long &i,long &j);
extern void Draw(float ,float );
extern void Draw(long ,long );
extern void DrawMark(R2 r);
//inline void DrawMark(D2 r) {DrawMark(R2(r.x,r.y));}
inline void Move(I2 x) {IMoveTo(x.x,x.y);}
inline void Move(R2 x) {rmoveto(x.x,x.y);}
//inline void Move(D2 x) {rmoveto(x.x,x.y);}
inline void Line(I2 x){ILineTo(x.x,x.y);}
inline void Line(R2 x) {rlineto(x.x,x.y);}
//inline void Line(D2 x) {rlineto(x.x,x.y);}
#endif

}






