// ********** DO NOT REMOVE THIS BANNER **********
// ORIG-DATE:     Jan 2008
// -*- Mode : c++ -*-
//
// SUMMARY  : Generic Tree header  (binairy, Quad, Oct)   
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

namespace Fem2D {
#include "R3.hpp"
}
namespace EF23 {
  using namespace Fem2D;
  using Fem2D::R;

  static const int MaxDeep = 30;
  typedef  int  IntQuad;
  typedef  long long Int8;
  static const IntQuad MaxISize = ( 1L << MaxDeep);
  static const IntQuad MaxISize1 =   MaxISize-1; 
class Z0 {
public:
  //  static bool INTER_SEG1d(int a,int b,int x,int y) { return (((y) > (a)) && ((x) <(b)));}
    static bool INTER_SEG1dn(int a,int b,int x0,int d) {
        int x= x0-d, y = x0+d;
        if(verbosity>999 && y < x0  ) cout << "overflow ?? inter_1dn" << " ] "<< a << "," << b << "[ inter ] " << x << "," << y << "[\n";
        if(y<0) y = std::numeric_limits<int>::max();// overflow in ..

        return (((y) > (a)) && ((x) <(b)));}

};
  class Z1 : public Z0 { public:
    //  static bool INTER_SEG1d(int a,int b,int x,int y) { return (((y) > (a)) && ((x) <(b)));}
    int x;
    Z1():x(0){}
    Z1(R1 P) : x((int)P.x) {}
    Z1( int i) : x(i){}
    Z1(const Z1 &pp,int k,int l): x(pp.x+(( k&1) ? l : 0)) {}
    void Add(int k,int l) { x+= (( k&1) ? l : 0) ;}
    Z1(const Z1 &A,const Z1 &B) : x(B.x-A.x){}
	
    int Case(int l) const  { return  ( x & l)? 1 : 0 ;}
    int norm() const { return abs(x);}
    Int8 norm2() const { return (Int8) x* (Int8) x;}
    void Bound() {   x = max(min(x,MaxISize1),0);}
    
    bool less(Z1 h) const  { return abs(x) <h.x ;}
    bool interseg(Z1 pp,int hb,int h) const { 
      return INTER_SEG1dn(x,x+hb,pp.x,h) ;}
    bool interseg(const Z1 &pp,int hb,const Z1 &h) const {
      return INTER_SEG1dn(x,x+hb,pp.x,h.x) ;}
    operator R1 () const { return R1(x);}
  };
  
 
  
  class Z2 : public Z0 { public:
 
    int x,y;
    Z2():x(0),y(0) {}
    Z2(R2 P) : x((int)P.x),y((int)P.y) {}
    Z2( int i) : x(i),y(i){}
    //Z2( int i,int j) : x(i),y(j) {}
    Z2(const Z2 &pp,int k,int l): x(pp.x+(( k&1) ? l : 0)),y(pp.y+(( k&2) ? l : 0)) {}
    void Add(int k,int l) { x+= (( k&1) ? l : 0) ; y+= (( k&2) ? l : 0);}
    Z2(const Z2 &A,const Z2 &B) : x(B.x-A.x),y(B.y-A.y) {}
	
    int Case(int l) const  { return ( ( y & l) ? (( x & l) ? 3 : 2 ) : ( ( x & l)? 1 : 0 )  ) ;}
    int norm() const { return Max(abs(x),abs(y));}
    Int8 norm2() const { return (Int8) x*(Int8) x + (Int8) y*(Int8)y;}
    void Bound() {   x = max(min(x,MaxISize1),0);
      y = max(min(y,MaxISize1),0);}
    
    bool less(Z2 h) const  { return abs(x) <h.x && abs(y) <h.y;}
    bool interseg(Z2 pp,int hb,int h) const { 
      return INTER_SEG1dn(x,x+hb,pp.x,h) && INTER_SEG1dn(y,y+hb,pp.y,h);}
    bool interseg(const Z2 &pp,int hb,const Z2 &h) const {
      return INTER_SEG1dn(x,x+hb,pp.x,h.x)
          && INTER_SEG1dn(y,y+hb,pp.y,h.y);}
    operator R2 () const { return R2(x,y);}
  };
  
  class Z3 : public Z0{ public:
  
    int x,y,z;
    Z3():x(0),y(0),z(0) {}
    
    Z3(R3 P) : x((int)P.x),y((int)P.y),z((int) P.z) {}
    Z3( int i) : x(i),y(i),z(i){}
    
    Z3(const Z3 &pp,int k,int l): x(pp.x+(( k&1) ? l : 0)),y(pp.y+(( k&2) ? l : 0)),z(pp.z+(( k&4) ? l : 0)) {}
    void Add(int k,int l) { x+= (( k&1) ? l : 0) ; y+= (( k&2) ? l : 0); z+= (( k&4) ? l : 0);}
    Z3(const Z3 &A,const Z3 &B) : x(B.x-A.x),y(B.y-A.y),z(B.z-A.z) {}
    void Bound() {  x = max(min(x,MaxISize1),0);
      y = max(min(y,MaxISize1),0);
      z = max(min(z,MaxISize1),0);}
    
    int Case(int l) const  {// cout << " case = "<< int((x&l)!=0)+(int((y&l)!=0)<<1) + (int((z&l)!=0)<<2);
      return int( (x&l)!=0) + ( int((y&l)!=0)<<1 ) + ( int( (z&l)!=0) <<2 ) ;}
    int norm() const { return Max(abs(x),abs(y),abs(z));}
    Int8 norm2() const { return (Int8) x*(Int8) x + (Int8) y*(Int8)y + (Int8) z*(Int8)z; }
    bool less(Z3 h) const  { return abs(x) <h.x && abs(y) <h.y && abs(z) < h.z ;}
    bool interseg(Z3 pp,int hb,int h) const { 
      return INTER_SEG1dn(x,x+hb,pp.x,h) && INTER_SEG1dn(y,y+hb,pp.y,h) && INTER_SEG1dn(z,z+hb,pp.z,h) ;
      }
    bool interseg(const Z3 &pp,int hb,const Z3 &h) const { 
      return INTER_SEG1dn(x,x+hb,pp.x,h.x) 
          && INTER_SEG1dn(y,y+hb,pp.y,h.y)
          && INTER_SEG1dn(z,z+hb,pp.z,h.z);
      }
    operator R3 () const { return R3(x,y,z);} 
    
  }; 

  inline  ostream& operator <<(ostream& f, const Z3 & P )   { f << P.x << ' ' << P.y << ' ' << P.z   ; return f; }
  inline  ostream& operator <<(ostream& f, const Z2 & P )   { f << P.x << ' ' << P.y   ; return f; }
  inline  ostream& operator <<(ostream& f, const Z1 & P )   { f << P.x    ; return f; }
  
  template<class Rd>    struct Traits_Zd {  typedef void Zd;};
  template<>    struct Traits_Zd<R1> {  typedef Z1 Zd;};
  template<>    struct Traits_Zd<R2> {  typedef Z2 Zd;};
  template<>    struct Traits_Zd<R3> {  typedef Z3 Zd;};
  
  template<class Vertex>    
      class GTree {
    typedef typename Vertex::Rd Rd;
    typedef typename Traits_Zd<Rd>::Zd Zd;
    
    
  public:
  
    static  const int d =Rd::d;
    static const int N = 1 << d;  // N=2^(d-1)
    
    
    class QuadTreeBox { 
    public:
      
      int n; // if n < 4 => Vertex else =>  QuadTreeBox;
      union {
	QuadTreeBox *b[N];
	Vertex * v[N];
      };
      // void init() { for(int i=0;i<N;++i) b[i]=0;}
      
    }; // end class QuadTreeBox  /////////////////
    
    class StorageQuadTreeBox {
    public:
      QuadTreeBox *b,*bc,*be;
      int len;
      StorageQuadTreeBox *n; // next StorageQuadTreeBox
      StorageQuadTreeBox(int ,StorageQuadTreeBox * =0);
      ~StorageQuadTreeBox();
      int  SizeOf() const {
	return len*sizeof(QuadTreeBox)+sizeof(StorageQuadTreeBox)+ (n?n->SizeOf():0);
      }
    }; // end class  StorageQuadTreeBox 
    
    StorageQuadTreeBox * sb;
    
    
    int  lenStorageQuadTreeBox;
    
  public:
    QuadTreeBox * root;
    // Mesh *th;
    
    long NbQuadTreeBoxSearch,NbVerticesSearch,NbSearch;
    long NbQuadTreeBox,NbVertices;
    
    Rd cMin,cMax; //  box of QuadTree
    R coef; //	
    
    
    Zd  RdtoZd(const Rd &P)  const {return Zd((Minc(Maxc(P,cMin),cMax)-cMin)*coef);} 
    Zd  VtoZd(const Vertex * v) const {return RdtoZd( (const Rd&) *v);} 
    Zd  VtoZd(const Vertex & v) const {return RdtoZd( (const Rd&) v);} 
    
    Rd  ZdtoRd(const Zd &I) const { return ( (Rd) I )/coef+cMin;}
    
    Vertex * NearestVertex(const Rd & P,bool trueNearest=false) {
      return NearestVertex(RdtoZd(P),trueNearest);} //XtoI(P.x),YtoJ(P.y));}
    Vertex * NearestVertexWithNormal(const Rd & P);
    Vertex * NearestVertex(Zd i2,bool trueNearest=false);
    
    Vertex *  ToClose(const Rd & ,R ,Zd, bool nearest=false );
    Vertex *  ToClose(const Rd & P,R delta,bool nearest=false){
      int hx = (int) (coef*delta);
       hx = hx>0 ? hx:1; // bof bof ....
       ffassert(hx>0);// bug if too small dec. 2019 FH.
      //if(verbosity > 5 ) cout << "hx=" << hx << " coef=" << coef << endl;
      return ToClose(P,delta,Zd(hx),nearest);}
    int SizeOf() const {return sizeof(GTree)+sb->SizeOf();}
    
    void  Add( Vertex & w);
    
    QuadTreeBox* NewQuadTreeBox()
  {
    ///cout << "NewQuadTreeBox " << sb << " " << sb->bc << " " 
    //<< sb->be << " " <<lenStorageQuadTreeBox <<endl;
    if(! (sb->bc<sb->be)) 
      sb=new StorageQuadTreeBox(lenStorageQuadTreeBox,sb);
    
    assert(sb && (sb->bc->n == 0));
    NbQuadTreeBox++;
    return sb->bc++;
  }
    ~GTree();
    GTree(Vertex * v,Rd Pmin,Rd Pmax,int nbv);
    GTree();
    template<class V>     
    friend ostream& operator <<(ostream& f, const  GTree<V> & qt);
  // add FH mars 2020 , for new search .. FH...
    int  ListNearestVertex(Vertex **lnv,int nvn,double delta,Rd P)
          {     int hx = (int) (coef*delta);
              hx = hx>0 ? hx:1; // bof bof ....
              return ListNearestVertex(lnv,nvn,hx,RdtoZd(P));}
          Vertex * TrueNearestVertex(Zd i2) {return NearestVertex(i2,true);}
    int ListNearestVertex(Vertex **lnv,int nlvnx,int dh,Zd xyi);
          

    Vertex * TrueNearestVertex(const Rd & P) {
        return TrueNearestVertex(RdtoZd(P));}

  };
  
  template<class Mesh>
  const typename  Mesh::Element * Find(const Mesh & Th,
				       GTree<typename Mesh::Vertex> *quadtree,
				       typename Mesh::Rd P,
				       typename Mesh::RdHat & Phat,
				       bool & outside,
				       const typename  Mesh::Element * tstart);
  
  
} // name space 

