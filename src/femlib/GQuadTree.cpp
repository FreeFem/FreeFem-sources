// ORIG-DATE:     Jan 2008
// -*- Mode : c++ -*-
//
// SUMMARY  : Generic Tree (binairy, Quad, Oct)   
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


#include <cmath>
#include <cstdlib>
#include "error.hpp"
#include <iostream>
#include <climits>
#include  <cfloat>
#include <cstring>
#include "ufunction.hpp"
#include "HeapSort.hpp"
using namespace std;


#include "GenericMesh.hpp"
#include "Mesh1dn.hpp"
#include "Mesh2dn.hpp"
#include "Mesh3dn.hpp"

#include  <set>

extern   long npichon2d, npichon3d;
extern   long npichon2d1, npichon3d1;



namespace EF23 {
  
  //  new version ----------
  //  ----------------------
  template<class Rd>
  void OrthoProject(const Rd &P,const Rd &A,const Rd & B,R * l)
  {
    Rd AB(A,B),AP(A,P),BP(B,P);
    R pa=(AB,AP),pb=(AB,BP);
    //  (l0 PA + l1 PB , AB)=0
    //   l0 pa + l1 pb = 0
    //   l0 + l1 =1
    //   l0 (pa - pb) = - pb;
    l[0] = - pb/(pa-pb);
    l[1] = 1-l[0];
  }
  
  template<class Rd>
  void OrthoProject(const Rd &P,const Rd &A,const Rd &B,const Rd &C,R * l)
  {
  Rd AB(A,B),AC(A,C),AP(A,P),BP(B,P),CP(C,P);
  R2 p0( (AB,AP) , (AC,AP) );
  R2 p1( (AB,BP) , (AC,BP) );
  R2 p2( (AB,CP) , (AC,CP) );
  // sum  li pi = 0
  //   l0    + l1    + l2     = 1
  //  =>     
  R2 O;
  R d = det(p0,p1,p2);
  l[0] = det(O ,p1,p2)/ d;
  l[1] = det(p0,O ,p2)/ d;
  l[2] = 1. -l[0] -l[1];
  //  
  //     
  }
    
template<class Vertex>
    int  GTree<Vertex>::ListNearestVertex(Vertex **lnv,int nlvnx,int dh,Zd xyi)
    {
        // warning this function return the NearestVertex in the first
        // none empty box contening the point xyi.
        //   They do not return the true nearest point in classical norm.
        int nlnv =0;
        QuadTreeBox * pb[ MaxDeep ];
        int  pi[ MaxDeep  ];
        Zd pp[  MaxDeep ];
        int l=0; // level
        QuadTreeBox * b;
        Int8  h2=(Int8)dh*dh,h0;
        IntQuad h=dh,hb =  MaxISize;
        Zd   p0(0,0,0);
        Zd  plus(xyi) ; plus.Bound();
        
        // xi<MaxISize?(xi<0?0:xi):MaxISize-1,yj<MaxISize?(yj<0?0:yj):MaxISize-1);
        
        Vertex *vn=0;
        
        // init for optimisation ---
        b = root;
        long  n0=0;
        if (!root->n)
            return 0; // empty tree
      
        if(verbosity>2000)
            cout << "        general case : NearVertex" << endl;
        // general case -----
        pb[0]= b;
        pi[0]=b->n>0 ?(int)  b->n : N  ;
        pp[0]=p0;
        
        do {
            b= pb[l];
            while (pi[l]--)
            {
                int k = pi[l];
               
                if (b->n>0) // Vertex QuadTreeBox none empty
                {
                    NbVerticesSearch++;
                    Zd i2 =  VtoZd(b->v[k]);
                    h0 = Zd(i2,plus).norm2();//  NORM(iplus,i2.x,jplus,i2.y);
                    if (h0 <h2)
                    {// on stock ..
                        if( nlnv < nlvnx)
                          lnv[nlnv++] = b->v[k];
                    }
                }
                else // Pointer QuadTreeBox
                {
                    QuadTreeBox *b0=b;
                    if ((b=b->b[k]))
                    {
                        hb >>=1 ; // div by 2
                        Zd ppp(pp[l],k,hb);
                        
                        if  ( ppp.interseg(plus,hb,h) )//(INTER_SEG(iii,iii+hb,iplus-h,iplus+h) && INTER_SEG(jjj,jjj+hb,jplus-h,jplus+h))
                        {
                            NbQuadTreeBoxSearch++;// add box search 
                            pb[++l]=  b;
                            pi[l]= b->n>0 ?(int)  b->n : N  ;
                            pp[l]=ppp;
                        }
                        else
                            b=b0, hb <<=1 ;
                    }
                    else
                        b=b0;
                }
            }
            hb <<= 1; // mul by 2
        } while (l--);
        
        return nlnv;
    }
    

  template<class Vertex>    
  Vertex *  GTree<Vertex>::NearestVertex(Zd xyi,bool trueNearest)//long xi,long yj)
  {
    // warning this function return the NearestVertex in the first
    // none empty box contening the point xyi.
    //   They do not return the true nearest point in classical norm.
      
    QuadTreeBox * pb[ MaxDeep ];
    int  pi[ MaxDeep  ];
    Zd pp[  MaxDeep ];
    int l=0; // level
    QuadTreeBox * b;
    Int8  h2=(Int8) MaxISize*(Int8)MaxISize*3,h0;
    IntQuad h=MaxISize,hb =  MaxISize;
    Zd   p0(0,0,0);
    Zd  plus(xyi) ; plus.Bound();
      int c2infty = 1+trueNearest;// 2 if trueNearest
    // xi<MaxISize?(xi<0?0:xi):MaxISize-1,yj<MaxISize?(yj<0?0:yj):MaxISize-1);
    
    Vertex *vn=0;
    
    // init for optimisation ---
    b = root;
    long  n0=0;
    if (!root->n)
      return vn; // empty tree 
    if( ! trueNearest )
    while( (n0 = b->n) < 0) 
      {
		// search the non empty 
		// QuadTreeBox containing  the point (i,j)
		long hb2 = hb >> 1 ;
		int k = plus.Case(hb2);//(iplus,jplus,hb2);// QuadTreeBox number of size hb2 contening i;j

		QuadTreeBox * b0= b->b[k];
 
		if ( ( b0 == 0) || (b0->n == 0) ){
		  break; // null box or empty box   => break 
		}   
		NbQuadTreeBoxSearch++;
		b=b0;
		p0.Add(k,hb2);
		hb = hb2; 
		
      }
      // n0 number of boxes of in b ("b0")
    if(verbosity>2000)
    cout << "        n0=" << n0 << " " << trueNearest <<endl;
    
    if ( n0 > 0) 
    {  
      for( int k=0;k<n0;k++)
	{
	  Zd i2 =  VtoZd(b->v[k]);
	  h0 = Zd(i2,plus).norm2();  //warning ..  norm sup and  not norm 2 ...
	  if (h0 <h2) {
	    h2 = h0;
	    vn = b->v[k];
              ffassert(vn);
	  }
	  NbVerticesSearch++;
	}
        if(verbosity>2000)
            cout << "        find " << vn << " " << h0 << " " << h2 << " " << endl;

      return vn;
    }
    
if(verbosity>2000)
    cout << "        general case : NearVertex" << endl;
    // general case -----
  pb[0]= b;
  pi[0]=b->n>0 ?(int)  b->n : N  ;
  pp[0]=p0;
  h=hb;
  do {    
    b= pb[l];
    while (pi[l]--)
      { 	      
        int k = pi[l];
	
	if (b->n>0) // Vertex QuadTreeBox none empty
	  { 
	    NbVerticesSearch++;
	    Zd i2 =  VtoZd(b->v[k]);
	    h0 = Zd(i2,plus).norm2();// norme 2 ^2
	    if (h0 <h2)
	      {
		h2 = h0;
                h =Zd(i2,plus).norm()*c2infty;// norm infty -> norm 2 big enought
		vn = b->v[k];
                  if(verbosity>2000)
                      cout << "        find   " << vn << " " << h0 << " " << h2 << " " << h << endl;

	      }
	  }
	else // Pointer QuadTreeBox 
	  { 
	    QuadTreeBox *b0=b;
	    NbQuadTreeBoxSearch++;
	    if ((b=b->b[k]))
	      {
		hb >>=1 ; // div by 2
		Zd ppp(pp[l],k,hb);
		
		if  ( ppp.interseg(plus,hb,h) )//(INTER_SEG(iii,iii+hb,iplus-h,iplus+h) && INTER_SEG(jjj,jjj+hb,jplus-h,jplus+h)) 
		  {
		    pb[++l]=  b;
		    pi[l]= b->n>0 ?(int)  b->n : N  ;
		    pp[l]=ppp;
		  }
		else
		  b=b0, hb <<=1 ;
	      }
	    else
	      b=b0;
	  }
      }
    hb <<= 1; // mul by 2 
  } while (l--);
  
  return vn;
  }
  
  
  template<class Vertex>
  Vertex *  GTree<Vertex>::ToClose(const Rd & v,R seuil,Zd H,bool nearest)
  {
    const Rd X(v);
    const Zd p(RdtoZd(v));
    R seuil2 = seuil*seuil;
    // const Metric  Mx(v.m);
    Vertex * pvr=0; 
    QuadTreeBox * pb[ MaxDeep ];
    int  pi[ MaxDeep  ];
    Zd pp[  MaxDeep ];
    
    int l=0; // level
    QuadTreeBox * b;
    long h=MaxISize;
    long hb =  MaxISize;
    Zd p0( 0 );
    
    if (!root->n)
      return 0; // empty tree 
   
    // general case -----
    pb[0]= root;
    pi[0]= root->n>0 ?(int)  root->n : N ;
    pp[0]= p0;
    h=hb;
    do {    
      b= pb[l];
      while (pi[l]--)
	{ 	      
	  int k = pi[l];
	  //cout << "b" << b << ", k= " << k << endl;
	  //cout << " b->n " << b->n << endl;
	  if (b->n>0) // Vertex QuadTreeBox none empty
	    { 
	      NbVerticesSearch++;
	      Vertex & V(*b->v[k]);
	      Zd i2 =  VtoZd(V);
	   
	      if ( Zd(i2,p).less(H) )
		{
		  
		  Rd XY(X,V);
		  R dd;
		  if( (dd= (XY,XY) ) < seuil2 ) // LengthInterpole(Mx(XY), b->v[k]->m(XY)))  < seuil )
		    { 
                        if( nearest )  // modif FH  aug. 2019 For P. Jolivet.
                        {
                            seuil2=dd;
                            pvr = & V;
                        }
                        else
                            return &V;
		    }
		}
	    }
	  else // Pointer QuadTreeBox 
	    { 
	      QuadTreeBox *b0=b;
	      NbQuadTreeBoxSearch++;
	      if( (b=b->b[k]) ) 
		{
			
		  hb >>=1; // div by 2
		  Zd ppp(pp[l],k,hb);
		  	
		  if( ppp.interseg(p,hb,H) )
		    {
		      	
		      pb[++l]=  b;
		      pi[l]= b->n>0 ?(int)  b->n : N;  
		      pp[l]=ppp;
		    }
		  else
		    b=b0, hb <<=1 ;
		    	
		}
	      else
		b=b0;
	    }
	} // fin: while(pi[l]--)
      hb <<= 1; // mul by 2 
    } while (l--);

    return pvr;
    
  }
  
  template<class Vertex>
  void  GTree<Vertex>::Add( Vertex & w)
  {
    QuadTreeBox ** pb , *b;
    Zd p(VtoZd(w));
    long l=MaxISize;
    pb = &root;
   
    while( (b=*pb) && (b->n<0))
      { 
	b->n--;
	l >>= 1;
	pb = &b->b[p.Case(l)];
      }
    if  (b) {
      for(int i=N-1;i>=0;--i)
	if (b->n > i &&  b->v[i] == &w)
	  {
	    //if( abs(w.x+0.5)<1e-10 ) cout << "if (b->n > i &&  b->v[i] == &w)" <<endl;
	    return;
	  }
    }
    assert(l);
    while ((b= *pb) && (b->n == N)) // the QuadTreeBox is full
      { 
	//if(verbosity > 5) cout << " b = " << b << b->n << "  " << l << endl;
	Vertex *v4[N]; // copy of the QuadTreeBox vertices      
	for(int i=0;i<N;++i)
	  { v4[i]= b->v[i];
	  b->v[i]=0;}
	
	b->n = -b->n; // mark is pointer QuadTreeBox
	
	l >>= 1;    // div the size by 2
	ffassert(l);
	for (int k=0;k<N;k++) // for the 4 vertices find the sub QuadTreeBox ij
	  { 
	    int ij;
	    QuadTreeBox * bb =  b->b[ij=VtoZd(v4[k]).Case(l)];
	    //if(verbosity > 5)  cout << "ij= "<< ij <<  " " << VtoZd(v4[k])<< endl;
	    if (!bb) 
	      bb=b->b[ij]=NewQuadTreeBox(); // alloc the QuadTreeBox 
	    //if(verbosity > 5)  cout << bb << " " << k << " "  << ij <<  endl;
	    bb->v[bb->n++] = v4[k];
	  }
	pb = &b->b[p.Case(l)];
      }
    if (!(b = *pb))
      { //if(verbosity > 5)  cout << "Qbox \n";
	b=*pb= NewQuadTreeBox(); //  alloc the QuadTreeBox 
      }
    //if(verbosity > 5)  cout << b << " " << b->n << endl;
    b->v[b->n++]=&w; // we add the vertex 
    NbVertices++;    
    
  }

    
template<class Vertex>    
GTree<Vertex>::GTree(Vertex * v,Rd Pmin,Rd Pmax,int nbv) : 
 lenStorageQuadTreeBox(nbv/8+100),
 // th(t),
 NbQuadTreeBoxSearch(0),
 NbVerticesSearch(0),
 NbQuadTreeBox(0),
 NbVertices(0),
 cMin(Pmin-(Pmax-Pmin)/2),
 cMax(Pmax+(Pmax-Pmin)/2),
 coef( MaxISize/Norme_infty(cMax-cMin))
 
{ 
  if(verbosity>5){
    cout << "  GTree: box: "<<  Pmin << " " << Pmax << " " << cMin << " "<< cMax << " nbv : " << nbv <<endl;
    cout << "MaxISize " << MaxISize << endl; 
    //cout << "  RdtoZd(cMin)" << RdtoZd(cMin) << " RdtoZd(cMax)" << RdtoZd(cMax) << endl;
    //cout << "  RdtoZd(Pmin)" << RdtoZd(Pmin) << " RdtoZd(Pmax)" << RdtoZd(Pmax) << endl;
    }
  sb =new StorageQuadTreeBox(lenStorageQuadTreeBox);
  root=NewQuadTreeBox(); 
  //  throwassert( MaxISize > MaxICoor);
  if (v)
    for (long i=0;i<nbv;i++) 
     Add(v[i]);
}
  
  template<class Vertex>
  GTree<Vertex>::GTree() : 
  lenStorageQuadTreeBox(100),
  // th(0),
  NbQuadTreeBoxSearch(0),
  NbVerticesSearch(0),
  NbQuadTreeBox(0),
  NbVertices(0),
  cMin(),cMax(),coef(0)
  {
    sb =new StorageQuadTreeBox(lenStorageQuadTreeBox);
    root=NewQuadTreeBox();
  }
  
  template<class Vertex>    
  GTree<Vertex>::StorageQuadTreeBox::StorageQuadTreeBox(int ll,StorageQuadTreeBox *nn)
  {
    len = ll;
    n = nn;
    b = new QuadTreeBox[ll];
    for (int i = 0; i <ll;i++)
      {
	b[i].n =0;
	for(int j=0;j<N;++j)
	  b[i].b[j]=0;
      }
    bc =b;
    be = b +ll;
    assert(b);
  }
  
  template<class Vertex>    GTree<Vertex>::StorageQuadTreeBox::~StorageQuadTreeBox()
  { //cout <<  "~StorageQuadTreeBox " << this << " n " << n << " b " << b << endl;
    if(n) delete n;
    delete [] b;
  }
  
  template<class Vertex> GTree<Vertex>::~GTree()
  {
    delete sb; 
  }
  
template<class Vertex> ostream& operator <<(ostream& f, const  GTree<Vertex> & qt)
{ 
  f << " the tree "  << endl;
  f << " NbTreeBox = " << qt.NbQuadTreeBox 
    << " Nb Vertices = " <<  qt.NbVertices << endl;
  f << " NbTreeBoxSearch " << qt.NbQuadTreeBoxSearch  
    << " NbVerticesSearch " << qt.NbVerticesSearch << endl;
  f << " SizeOf QuadTree" << qt.SizeOf() << endl;
  //     return  dump(f,*qt.root);
  return  f;
}
  
  template<class Vertex> Vertex *  GTree<Vertex>::NearestVertexWithNormal(const Rd &P)//(long xi,long yj)
  {
    QuadTreeBox * pb[ MaxDeep ];
    int  pi[ MaxDeep  ];
    Zd pp[ MaxDeep];
    int l; // level
    QuadTreeBox * b;
    IntQuad  h=MaxISize,h0;
    IntQuad hb =  MaxISize;
    Zd   p0;
    Zd  plus(RdtoZd(P) );//xi<MaxISize?(xi<0?0:xi):MaxISize-1,yj<MaxISize?(yj<0?0:yj):MaxISize-1);
    
    Vertex *vn=0;
    // init for optimisation ---
    b = root;
    long  n0;
    if (!root->n)
      return vn; // empty tree 
    
    while( (n0 = b->n) < 0) 
      {
	// search the non empty 
	// QuadTreeBox containing  the point (i,j)
	long hb2 = hb >> 1 ;
	int k = plus.Case(hb2);//(iplus,jplus,hb2);// QuadTreeBox number of size hb2 contening i;j
	QuadTreeBox * b0= b->b[k];
	if ( ( b0 == 0) || (b0->n == 0) ) 
	  break; // null box or empty   => break 	    
	NbQuadTreeBoxSearch++;	
	p0.Add(k,hb2);	
	hb = hb2; 
      }
    
    
    if ( n0 > 0) 
    {  
      for(int k=0;k<n0;k++)
	{
	  Vertex * v=b->v[k];
	  if (v->ninside(P)) {
	    Zd i2 =  VtoZd(v);
	    //   try if is in the right sens -- 
	    h0 = Zd(i2,plus).norm();// h0 = NORM(iplus,i2.x,jplus,i2.y);
	    if (h0 <h) {
	      h = h0;
	      vn = v;}
	    NbVerticesSearch++;}
	}
      if (vn) return vn; 
    }
    // general case -----
    // INITIALISATION OF THE STACK 
    l =0; // level 
    pb[0]= b;
    pi[0]= b->n>0 ?(int)  b->n : N  ;
    pp[0]=p0;
    h=hb;
  L1: 
    do {   // walk on the tree  
      b= pb[l];
      while (pi[l]--) // loop on 4 element of the box
	{ 	      
	  int k = pi[l];
	  
	  if (b->n>0) // Vertex QuadTreeBox none empty
	    { 
	      Vertex * v=b->v[k];
	      if (v->ninside(P) ) {	    
		NbVerticesSearch++;
		Zd i2 =  VtoZd(v);
		// if good sens when try -- 
		h0 = Zd(i2,plus).norm();//  NORM(iplus,i2.x,jplus,i2.y);
		if (h0 <h) 
		{
		  h = h0;
		  vn =v;
		}}
	    }
	  else // Pointer QuadTreeBox 
	    { 
	      QuadTreeBox *b0=b;
	      NbQuadTreeBoxSearch++;
	      if ((b=b->b[k])) 
		{
		  hb >>=1 ; // div by 2
		  Zd ppp(pp[l],k,hb);
		
		  if  ( ppp.interseg(plus,hb,h) )//(INTER_SEG(iii,iii+hb,iplus-h,iplus+h) && INTER_SEG(jjj,jjj+hb,jplus-h,jplus+h)) 
		    {
		      pb[++l]=  b;
		      pi[l]= b->n>0 ?(int)  b->n : N  ;
		      pp[l]=ppp;		    
		    }
		  else
		    b=b0, hb <<=1 ;
		}
	      else
		b=b0;
	    }
	}
      hb <<= 1; // mul by 2 
    } while (l--);
    if (!vn && b != root )
   {// cas particulier on repart du sommet on avais rien trouver 
     b=root;
     hb =  MaxISize;
     p0=Zd();
     l=0;
     pb[0]= b;
     pi[0]= b->n>0 ?(int)  b->n : N  ;
     pp[0]=Zd();
     
     goto L1;
   }
    return vn;
}




  //  static int kfind=0;
  // static int kthrough=0;
  
  inline void CoorBary(const Triangle2 & K,const  R2  & P, R *l)
  {
    R detK = 2.*K.mesure() ;
    l[1]=det(K[0],P,K[2])/detK;
    l[2]=det(K[0],K[1],P)/detK;
    l[0]=1-l[1]-l[2];
  }
  
    

  inline void CoorBary(const Tet & K,const  R3  & P, R *l)
  {
    R detK = 6.*K.mesure() ;
    l[1]=det(K[0],P,K[2],K[3])/detK;
    l[2]=det(K[0],K[1],P,K[3])/detK;
    l[3]=det(K[0],K[1],K[2],P)/detK;
    l[0]=1-l[1]-l[2]-l[3];
  }
    
    
  inline int  CoorBaryPos(const Triangle2 & K,const  R2  & P, R *l)
    {
        CoorBary(K,P,l);
        int n=0,m=-1;
        int nl[Tet::nv+1];
        R eps = -1e-10;
        for(int i=0;i<Triangle2::nv;++i)// correction FH 2/04/20
            if( l[i] < eps){
                nl[n++]=i;
            }
            else m = i;
        if( m>=0)
        { // proj P on face .. m
            int i0=Triangle2::nvadj[m][0];
            int i1=Triangle2::nvadj[m][1];
            R2 A(K[i0]),B(K[i1]);
            R2 AB(A,B),AP(A,P);
            double l2= AB.norme2();
            // L = (AB,AP)/l2
           
            l[i0] = max(0.,min(1.,(AB,AP)/l2));
            l[i1] = 1-l[i0];
            l[m]=0;
        }
        return n;
        
        
        
    }
    inline int  CoorBaryPos(const Tet & K,const  R3  & P, R *l)
    {
        CoorBary(K,P,l);
        int n=0,m=-1;
        int nl[Tet::nv+1];
        R eps = -1e-10;
        for(int i=0;i<Tet::nv;++i)// correction FH 2/04/20 
            if( l[i] < eps){
                nl[n++]=i;
            }
            else m = i;
        if( m>=0)
        { // proj P on face .. m
            int i0=Tet::nvadj[m][0];
            int i1=Tet::nvadj[m][1];
            int i2=Tet::nvadj[m][2];
            R3 A(K[i0]),B(K[i1]),C(K[i2]);
            R3 AB(A,B), AC(A,C);
            R3 N = AB^AC ;
            double N2 = (N,N);
            l[i0] = max(0.,min(1.,det(P,B,C,N)/N2));
            l[i1] = max(0.,min(1.,det(A,P,C,N)/N2));
            l[i2] = 1-l[i0]-l[i1];
            l[m]=0;
        }
        return n;
        

        
    }
  
//  inline    int nRand(int n) {
//    return  rand()%n; //avant random()%n;
//  }
  
  inline int find5(int i,int *k,int l)
  {
    if(l<5)
      {
	for(int j=0;j<l;++j)
	  if(i==k[j]) return j;
      }
    else  if(i>=k[0] && i<=k[l-1])
      {
	int i0=0,i1=l-1;
	while(i0<=i1)
	  { int im=(i0+i1)/2;
	    if(i<k[im]) i1=im-1;
	    else 
	      if(i>k[im]) i0=im+1;
	      else
		if(i==k[im]){
		  return im;
		}
	  }
      }
    return -1;
  } 
  
  template<class Mesh>
  const typename  Mesh::Element * Find(const Mesh & Th,
				       GTree<typename Mesh::Vertex> *quadtree,
				       typename Mesh::Rd P,
				       typename Mesh::RdHat & Phat,
				       bool & outside,
				       const typename  Mesh::Element * tstart)
{
  typedef  typename Mesh::Element Element;
  typedef  typename Mesh::Vertex Vertex;
  typedef  typename Mesh::Rd Rd;
  const int nkv=Element::nv;
  const int d=Rd::d;
  R dP=DBL_MAX, nddd=0;
  Rd PPhat,Delta;
  R eps = -1e-10;// def test in tet ..

  int k=0;
    const int nReStartMax = 1 << Rd::d; //  Nb vertex of the d-cube/
    int nReStart = 0;
    int itstart[100],itout[100],kstart=0;
    Rd Pout[100];
  int it,j,it00,nbdeja=0,nbdejax=0;
  const int mxbord=1000;
  int kbord[mxbord+1];
  int nbord=0;
  map<int,int> deja;
  if(searchMethod>1) goto PICHON;
  if ( tstart )
    it00=it =  Th(tstart);
  else  if(quadtree)
    {  
      const Vertex * v=quadtree->NearestVertex(P);// Change Juil 2017 ..
/*
      if (!v) 
	{ 
	  v=quadtree->NearestVertex(P);
	  assert(v);
	}*/
        if( !v) {
            verbosity = 10000000;
            cout << "Bug search  " << P << endl;
            v=quadtree->NearestVertex(P);
        }
        ffassert(v); // bug !!!!
      it00=it=Th.Contening(v);
      nddd=  Norme2(P-*v);
      if(verbosity>200)
          cout <<  "   Find: Close : "<<  *v << " , " << Th(v) << " dist " << nddd  << endl;
      
    }
  else ffassert(0);
RESTART:
  ffassert(kstart<100);
  itstart[kstart++]=it;
  if(verbosity>199)
    cout << "    " << nReStart << " tstart= " << tstart << " , it=" << it << " P="<< P << " " << kstart << endl;
  outside=true; 
  Mesh::kfind++;
  while (1)
    { 
      //if(verbosity>199) cout << "it " << it <<endl;
      const Element & K(Th[it]);
      Mesh::kthrough++;
      ffassert(k++<2000);
      if( k> 500)
      { // boucle ????  change eps
          nbdeja=++deja[it];// we stocke le deja view in map
        
          nbdejax = max(nbdeja,nbdeja);
          if(nbdejax> 5) eps=1e-7;
          else if(nbdejax> 4) eps=1e-8;
          else if(nbdejax> 3) eps=1e-9;
          
          ffassert(nbdeja<100);
      }
    
      int kk,n=0,nl[nkv];
      R l[nkv];
      for(int iii=0; iii<nkv; iii++)
	l[iii]=0.;
      
      CoorBary(K,P,l);
      
      // CoorBary :: donner entre 0 et 1
      // Pb si K.mesure*1.e-10 precision machine ==> bug
      
      // avant:
      // R eps =  -K.mesure()*1e-10;
        if( nbdeja>3) eps=-1e-8;
        
        
      for(int i=0;i<nkv;++i)
	if( l[i] < eps){
	  nl[n++]=i;
	}
        if(verbosity>19 && nbdeja >1) {
           cout << " Bizarre loop in search "<< nbdeja << "       tet it=" << it ;
           cout << "  K.mesure=" << K.mesure() ;
           cout << " eps=" << eps << " : " ;
           for(int i=0;i<nkv;++i)
              cout<< "  l["<< i <<"]=" <<  l[i] ;
           cout << " n=" << n << endl;
        }
      else if(verbosity>200){
	cout << "       tet it=" << it ;
	cout << "  K.mesure=" << K.mesure() ;
	cout << " eps=" << eps << " : " ;
	for(int i=0;i<nkv;++i)
	  cout<< "  l["<< i <<"]=" <<  l[i] ;
	cout << " n=" << n << endl;
      }
      
      if (n==0)
	{  // interior => return
	  outside=false; 
	  Phat=Rd(l +1);
	  // cout << Phat <<endl;
#ifndef NDEBUG
	  Rd pp=K(Phat)-P;
	  if( pp.norme()>1e-5)
	    {
	      cout << "     Bug find P " << P << " ==" << K(Phat) ;
	      cout << "Phat== " << Phat << " diff= " << pp << endl;
	      assert(0);
	    }
	  
#endif	  
	  return &K;
	}
      
      
      
      kk=n==1 ? 0 : randwalk(n);
      j= nl[ kk ];
      int itt =  Th.ElementAdj(it,j);
      if(itt!=it && itt >=0)  
	{
	  dP=DBL_MAX;
	  it=itt;
	  continue;
	}  
      int  inkbord=find5(it,kbord,nbord);
      if(inkbord<0 && nbord < mxbord)
	{
	  kbord[nbord++]=it;
	  if(nbord>=5)  HeapSort(kbord,nbord);
	}
      if(verbosity>1001)
	{
	  cout << "       bord "<< it<< "   nbf < 0 : " <<n << " (inb) " << inkbord << " nfb" << nbord<<endl;
	  R ss=0; 
	  for(int i=0;i<nkv;++i)
	    {  ss += l[i];
	      cout << l[i] << " ";}
	  cout << " s=" << ss << endl;;
	  
	}
      
      if(verbosity>2000)
	cout << "      GQuadTree::value of n " << n << endl;
      
      if ( n!=1 )  // on est sur le bord, mais plusieurs face <0 => on test les autre
	{  // 1) existe t'il un adj interne
	  int nn=0,ii;
	  int nadj[d+1],ni[d+1];
	  for(int i=0;i<nkv;++i)
	    //avant :: if (l[i] < eps && (itt=Th.ElementAdj(it,ii=i)) != it && itt && find5(itt,kbord,nbord) < -1 ) 
	    if (l[i] < eps && (itt=Th.ElementAdj(it,ii=i)) != it && (itt>=0) && find5(itt,kbord,nbord) < 0 ) 
	      ni[nn++]=i,nadj[i]=itt;
	  if(verbosity>1000)
	    cout << " nn : "<< nn << endl;
	  if (nn>0)
	    {
	      //j=nadj[nRand(nn)];
	      j=ni[randwalk(nn)];
	      it=nadj[j];
	      dP=DBL_MAX;
	      //cout << "new it= " << it << endl;
	      continue;
	    }
	}
      // toutes les faces <0  sont sur le bord.
      //   pour l'instant on s'arête
      // le point est externe. (mais trop cher pour faire mieux)
      // on projet le points sur le bord via le coordonne 
      //  barycentrique 
      { // a ridge on border  (to hard to do the correct stuff) 
	//  or a corner    just do the projection on lambda  
	R s=0.;
	for(int i=0;i<nkv;++i)
	  s += (l[i]= max(l[i],0.));
	for(int i=0;i<nkv;++i)
	  l[i]/=s;
	Phat=Rd(l +1);
	if(verbosity>1000)
	  {
	    cout << "      "<< P << " " << n << " l: ";
	    R ss=0; 
	    for(int i=0;i<nkv;++i)
	      {  ss += l[i];
		cout << l[i] << " ";}
	    cout << "   s=" << ss <<" " << s <<" exit by bord " << it << " "  << Phat << endl;;
          }
	outside=true;
// New Method mars 2020
        if(1) // searchMethod==0)
          {
          GenericDataFindBoundary<Mesh> * gdfb=Th.Buildgdfb();
          if(gdfb )
          {
              double l[4],lK[4];
              int na[4],nna=0,mma=0;
              int loutside;// of border ??? 0 inside, 1 out close, 2, out fare, , -1 inside
              int itt =gdfb->Find(P,l,loutside);
             // warning l is projecton on boundary ????
              CoorBary(Th[itt],P,lK);
              const double eps = 1e-10;
              for (int i=0; i<= Rd::d;++i )
                  if(lK[i] < -eps) na[nna++]=i;
              for( int k=0; k<nna;++k)
              {
                  int ii=na[k],ittt=Th.ElementAdj(itt,ii);
                  if( ittt<0 || ittt == itt) mma++;
              }
              it = itt;
          
              if( nna) //
              {
               if( nReStart++ < 1) goto RESTART;
               else
                  if(searchMethod) goto PICHON;
               Phat=Rd(l+1);
              }
              else // in itt
                Phat=Rd(lK+1);// correction FH 1/0>4/2020
              outside=nna>0;
              
              const Element &K=Th[it];
              if( verbosity > 9)
                  cout << "   - Find "<< P << " -> " << K(Phat) << " " << loutside << " k= " << itt
                  << " dist =" << (P-K(Phat)).norme() << " :: " << Phat << " nddd= " << nddd <<" " << nReStart <<  " " << nna <<  endl;
              return  Th.elements + it; // outside
              
          }
          }
// end New Methode
        // on retest with a other stating point??????
          // Mod.  23/02/2016 F. H
       itout[kstart-1]= it;
       Pout[kstart-1]=Phat;
       while(nReStart++ < nReStartMax)
        {// Loop on the vertex of a d-hyper-cude
            {
                int k= nReStart-1;
                int i = (nReStart-2)/2;
                
                for(int i=0; i<Rd::d; ++i)
                    Delta[i]=  ((1<<i) & k) ?  nddd : -nddd;
            }
            if( verbosity>199) {
                Rd pp=Th[it](Phat)-P;
                cout << "    restart find P " << P << " !=" <<Th[it](Phat)  << "\t" ;
                cout << "Phat== " << Phat << " diff= " << pp << " it = " << it << endl;
            }
            
            Rd PP= P + Delta;
            Vertex* v=quadtree->NearestVertex(PP);
            if( nddd ==0)  nddd= Norme2(P-*v);
            it=Th.Contening(v);
            bool same=false;

            for(int j=0;j<kstart ; ++j)
                if( it == itstart[j]) {same=true; break;}
            if( verbosity>199)
                cout << "   loop Search "<<nReStart << P << " Delta" <<  Delta << " it " << it << " same "<< same << endl;
            if(same) continue;
            if(Rd::d==2)  npichon2d1++;
            if(Rd::d==3)  npichon3d1++;
            
            goto RESTART;
        }
          
	if(searchMethod) goto PICHON;
          // Search best out ....
        int ko=0,k=0;
        double dmin=Norme2_2(P-Th[itout[k]](Pout[k]));
        for(int k=1; k<kstart; ++k)
          {
              double dk=Norme2_2(P-Th[itout[k]](Pout[k]));
              if( dk< dmin) { dmin=dk; ko=k;}
          }
        Phat=Pout[ko];
        it=itout[ko];
        if( verbosity>149)
            cout << " -- Find(out) " << it << " Outside " << outside << " it "<< it << " err = "<< sqrt(dmin) << "/" << Phat << " " << nddd << endl;
        return &Th[it] ;
      }		    
    }
  
 PICHON:
  {
      if(Rd::d==2)  npichon2d++;
      if(Rd::d==3)  npichon3d++;

  /*==============================PICHON=================*/
  // Brute force ** */    
  R l[4], eps = -1e-6;  //pichon
  double dist_baryC = 0.0, min_dist_baryC = 1e31;      
  long closestTet = -1;
  l[0]=l[1]=l[2]=l[3]=1; //  for d < 3
  for(int tet=0;tet<Th.nt;tet++) // nkv=4
    { 
      const Element & K(Th[tet]);
      CoorBary(K,P,l);
      
      // measure dist by sum ( |lambda_i| )
      dist_baryC = 0.0;      
      for(int i=0; i<nkv; i++)
	dist_baryC += abs(l[i]);
      
      // define closest Tetrahedron !!! TO VERIFY THE HYPOTHESE !!! 
      if( dist_baryC < min_dist_baryC )   
	{
	  min_dist_baryC = dist_baryC; 
	  closestTet = tet; 
	}
      
      
      if  ( (l[0] >= eps) && (l[1] >= eps) && (l[2] >= eps) &&  (l[3] >= eps) )
	{
         if( verbosity>199)
         {
             cout << "   #### brute force " << tet << " "<< endl;
         }
          outside=false;
	  Phat=Rd(l+1);
	  return &K;	
	}
    }
  
  const Element & K(Th[closestTet]);
  outside=true;
  CoorBaryPos(K,P,l);
  
  Phat=Rd(l+1);
  if(verbosity>2)
    cout << "  --vertex:" << P << " NOT in DOMAIN. use closestTet. Phat:" << Phat << endl;
  return &K;
  }
  /*==============================PICHON=================*/
}


// Instantiation  manuel des templates

template class GTree<Vertex2>;
template class GTree<Vertex3>;
template class GTree<Vertex1>;
    
///typedef MeshS::GMesh GMeshS;
typedef Mesh3::GMesh GMesh3;
typedef Mesh2::GMesh GMesh2;
typedef Mesh1::GMesh GMesh1;
 
    
/*template
const   GMeshS::Element * Find<GMeshS>(const GMeshS & Th,
                       GTree< GMeshS::Vertex> *quadtree,
                       GMeshS::Rd P,
                       GMeshS::RdHat & Phat,
                       bool & outside,
                       const   GMeshS::Element * tstart);
  */
 
template
const   GMesh3::Element * Find<GMesh3>(const GMesh3 & Th,
				       GTree< GMesh3::Vertex> *quadtree,
				       GMesh3::Rd P,
				       GMesh3::RdHat & Phat,
				       bool & outside,
				       const   GMesh3::Element * tstart);
template
const   GMesh2::Element * Find<GMesh2>(const GMesh2 & Th,
				       GTree< GMesh2::Vertex> *quadtree,
				       GMesh2::Rd P,
				       GMesh2::RdHat & Phat,
				       bool & outside,
				       const   GMesh2::Element * tstart);

/*
  template
  const   GMesh1::Element * Find<GMesh1>(const GMesh1 & Th,
				       GTree< GMesh1::Vertex> *quadtree,
				       GMesh1::Rd P,
				       GMesh1::RdHat & Phat,
				       bool & outside,
				       const   GMesh1::Element * tstart);
*/
}

template<typename Mesh>
GenericDataFindBoundary<Mesh>::~GenericDataFindBoundary()
{
    if(verbosity>1) cout << "\n -- ~GenericDataFindBoundary: Nb find " << nbfind << " nb element  "
    << nbelement << " ratio " << (double) nbelement/std::max(nbfind,1L)
    << " mpirank " << mpirank << endl;
    delete tree;
}
template<typename Mesh>
void GenericDataFindBoundary<Mesh>::gnuplot(const string & fn)
{ // for debugging ..
    if(dHat < 3)
    {
        ofstream gp(fn.c_str());
        ffassert(gp);
        
        const Mesh &Th = *pTh;
        if( bborder )
            for(int be=0; be<Th.nbe; ++be)
            {
                const BorderElement &B=Th.be(be);
                int e,k = Th.BoundaryElement(be,e);
                {
                    int ee=e, kk=  Th.ElementAdj(k,ee);
                    if ( kk >=0 || k != kk)
                    {
                        for(int j=0; j< BorderElement::nv;++j)
                            gp  << (Rd) B[j] << endl;
                        gp  << "\n\n";
                    }
                }
            }
        else
            for(int k=0; k<Th.nt; ++k)
            {
                const Element &K=Th[k];
                for(int j=0; j<= Element::nv;++j)
                    gp  << (Rd) K[j%Element::nv] << endl;
                gp  << "\n\n";
            }
        if( dHat==2)
        {
            ofstream gp((fn+"c").c_str());
            for(int i=0; i<P.N(); ++i)
            {
                int N=100;
                const double pi=3.14159265358979323846264338327950288 ;
                double dt = pi*2./N, r = delta[i];
                for(int j=0;j<=N; ++j)
                {
                    Fem2D::R3 Q=P[i];
                    double x = Q.x+r*cos(dt*j);
                    double y=  Q.y+r*sin(dt*j);
                    double z=  Q.z;
                    gp << x << " " << y << " " << z << endl;
                }
                gp << "\n\n";
            }}
    }
}
inline double Plambla(int dhat,const Fem2D::R2 *Q,Fem2D::R2 &P,double *l)
{
    ffassert(0); //  to do
}
inline double Plambla(int dhat,const Fem2D::R1 *Q,Fem2D::R1 &P,double *l)
{
    ffassert(0); //  to do
}

inline double  Plambla(int dhat,const Fem2D::R3 *Q,Fem2D::R3 &P,double *l)
{
    using Fem2D::R3 ;
    double dd=0;
    if(dhat==1)
    {
        R3 AB(Q[0],Q[1]), AP(Q[0],P);
        double ab2=(AB,AB);
        l[1] = (AP,AB)/ab2;
        l[0] = 1- l[1];
        dd= ab2;
    }
    else if(dhat==2)
    {
        R3 AB(Q[0],Q[1]), AC(Q[0],Q[2]);
        R3 N = AB^AC ;
        double N2 = (N,N);
        dd = N2;
        l[0] = det(Q[1]-P   ,Q[2]-P   ,N)/N2;
        l[1] = det(P   -Q[0],Q[2]-Q[0],N)/N2;
        l[2] = 1-l[0]-l[1];
    }
    else if( dhat==3)
    {
        double d=det(Q[0],Q[1],Q[2],Q[3]);
        dd = d;
        l[0]=det(P,Q[1],Q[2],Q[3])/d;
        l[1]=det(Q[0],P,Q[2],Q[3])/d;
        l[2]=det(Q[0],Q[1],P,Q[3])/d;
        l[3]= 1- l[0]-l[1]-l[2];
    }
    else assert(0);
    return dd;
}
inline double dist2(int dhat,const Fem2D::R2 *Q,Fem2D::R2 &P,double *l,double *dl)
{
    ffassert(0);
    return  0;
}
inline double dist2(int dhat,const Fem2D::R1 *Q,Fem2D::R1 &P,double *l,double *dl)
{
    ffassert(0);
    return  0;

}
inline double dist2(int dhat,const Fem2D::R3 *Q,Fem2D::R3 &P,double *l,double *dl)
{
    using Fem2D::R3 ;
    double d2=0;
    if(dhat==1)
    {
        dl[0]=min(1.,max(l[0],0.));
        dl[1]=1.-dl[0];
        R3 Pj = dl[0]*Q[0]+dl[1]*Q[1];
        d2 = R3(Pj,P).norme2();
        
    }
    else if(dhat==2)
    {
        int i=-1;
        if(l[0]<0) i=0;
        else if(l[1]<0) i=1;
        else if(l[2]<0) i=2;
        
        if( i>=0)
        {
            // project on edge i
            int i1= (i+1)%3,i2=(i+2)%3;
            R3 AB(Q[i1],Q[i2]), AP(Q[i1],P);
            dl[i2] = (AP,AB)/(AB,AB);
            dl[i2]=min(1.,max(dl[i2],0.));// proj on AB ..
            dl[i1] = 1- dl[i2];
            dl[i]=0;
        }
        else
        {
            dl[0]=l[0];
            dl[1]=l[1];
            dl[2]=l[2];
        }
        R3 Pj = dl[0]*Q[0]+dl[1]*Q[1]+dl[2]*Q[2];
        d2 = R3(Pj,P).norme2();
        
    }
    else if( dhat==3)
    {
        ffassert(0); // to do FH... non ...
    }
    else assert(0);
    return d2;
}
template<typename Mesh>
int GenericDataFindBoundary<Mesh>::Find(typename Mesh::Rd PP,double *l,int & outside) const
{  // FH: outside : 0 inside, 1 out close, 2, out fare, , -1 inside
    // warning l
    nbfind++;
    typedef double R;
    int nu=-1,ne=-1;
    R dnu= 1e200,deps=0;
    R dl[dHat+1];
    outside = 0;
    Vertex *p0= &P[0];
    Vertex *p =tree->TrueNearestVertex(PP);
    int ip = p-P;
    Rd Q[dHat+1];
    R ll[dHat+1],lpj[dHat+1];
//#define DEBUGGING
#ifdef DEBUGGING
    int err=0;
    double dispp=Rd(PP,*p).norme();
    for(int i=0; i<P.N(); ++i)
    {
        double l=Rd(PP,P[i]).norme();
        if( l < dispp*(1-1e-14) ) err++,cout << " bug "<< l << " " << dispp << " i "<< i << " ip" << ip << " " <<dispp-l << endl;
    }
    if(err)
    {// relay for debug
        p =tree->TrueNearestVertex(PP);
       ffassert(err==0);
    }
#endif
    Vertex  **plp = &lp[0];
    long lvp=tree->ListNearestVertex(plp,lp.N(), delta[ip],P[ip]);
    HeapSort(plp,lvp );
//  verif ListNearestVertex
    
#ifdef DEBUGGING
    {
        int err=0;
        set<int> lst;
        int nnnm=0,nnnp=0;
        double unp =(1+1e-14),unm=(1-1e-14);
        for(int j=0; j<lvp; ++j)
            lst.insert(lp[j]-p0);
        for(int i=0; i<P.N(); ++i)
        {
            double l=Rd(P[ip],P[i]).norme();
            if (l < delta[ip]*unp) nnnp++;
            if (l < delta[ip]*unm) nnnm++;
               
            if ( (l < delta[ip]*unp)   && (lst.find(i) == lst.end() ))
                {
                    err++;
                    cout << " bug "<< i << " not in set " << endl;
                }
        }
        ffassert(lvp<=nnnp && lvp >= nnnm);
        ffassert(err==0);
    }
#endif
    if( verbosity>19)
    {
        cout << ip << " , " << delta[ip] << " | " << lvp << " : ";
        for(int j=0; j<lvp; ++j)
        {
            int i = lp[j]-p0;
            int k = lp[j]->lab/Element::ne;
            int e = lp[j]->lab%Element::ne;
            cout << " "<< k << " ";
           // ffassert(i == k || );
        }
        cout << endl;
    }
    for(int j=0; j<lvp; ++j)
    {
        nbelement++;
        int k = lp[j]->lab/Element::ne;
        int e = lp[j]->lab%Element::ne;
        if(debug) cout << "    -- k = "<< k << " " << e << " " << j << endl;
        
        const Element & K=(*pTh)[k];
        int I[4]={0,1,2,3};
        int nI =dHat+1;
        if( bborder)
        {  // take just of part of Element ..
            nI= Element::nva;
            ffassert(nI==BorderElement::RdHat::d+1);
            for(int i=0; i< nI;++i)
                I[i]=Element::nvadj[e][i];
        }
        
        
        for(int i=0; i< nI;++i)
            Q[i]=K[I[i]];
        
        double ddeps=Plambla(nI-1,Q,PP,ll);//  return une taille ^2, ^4, ^3 suivant la dim  nI-1
        R d2 = dist2(nI-1,Q,PP,ll,lpj);
         if( dnu > d2)
        {
            if( nI==2)
                deps =ddeps/10000.;
            else if (nI==3)
                deps = sqrt(ddeps)/10000.; // epaisseur de l'objet au carre 1.100 de
            else ffassert(0);
            nu = k;
            ne= e;
            dnu=d2;
           
           
            if(  bborder)
            {  //  
                for(int i=0;i<=dHat;++i)
                    dl[i]=0;
                for(int i=0;i<nI;++i)
                    dl[I[i]]=lpj[i];
            }
            else
            {
            for(int i=0;i<=dHat;++i)
                dl[i]=lpj[i];
            }
        }
        if(verbosity>99) cout << "    Find " << k << " " << e << " / " << dnu  << endl;
    }
    ffassert(nu>=0);
    outside = dnu < deps ? 0
    : ((dnu< delta[ip]*delta[ip])? 1: 2) ;// 0 in, 1 near, 2 fare ..BofBof ..
    for(int i=0; i<= dHat;++i)
      l[i]=dl[i];

    if(debug)   cout << "  -- out nu "<< nu << " "<< ne << " , " << dnu <<" d_i " << delta[ip] << " :  "
        << l[1] << " " << l[2] << " "<< outside<<  endl;
    return nu;
}
template<typename Mesh>
int  TrueBorder(const Mesh &Th,typename Mesh::Vertex *P,double *delta)
{
    typedef typename Mesh::Vertex Vertex;
    
    typedef typename Mesh::Element Element;
    typedef typename Mesh::BorderElement BorderElement;
    typedef  typename Mesh::Rd Rd;
    typedef  typename Mesh::BorderElement::RdHat RdHat;
    static const int d = Rd::d;
    static const int dHat = RdHat::d;
    
    int nv =0;
    RdHat GHat(RdHat::diag(1./(dHat+1)));
    
    for(int be=0; be<Th.nbe; ++be)
    {
        const BorderElement &E=Th.be(be);
        int e,k = Th.BoundaryElement(be,e);
        {
            int ee=e, kk=  Th.ElementAdj(k,ee);
            if ( kk >=0 || k != kk)
            {
                E(GHat);
                Rd G(E(GHat));
                double l = 0;// 1.5 to be sure .. FH
                for(int i=0; i< BorderElement::nv ;++i)
                    l = max(l, Rd(G,E[i]).norme2()) ;
                delta[nv]=l;
                P[nv].lab= Element::ne*k+e;//  element and edge
                (Rd &) P[nv++]=G;
                
                
            }
        }
    }
    return nv;
}

template<typename Mesh>
GenericDataFindBoundary<Mesh>::GenericDataFindBoundary(Mesh const * _pTh,int ddebug)
: pTh(_pTh),tree(0),nbfind(0), nbelement(0), P(bborder ? pTh->nbe: pTh->nt),delta(P.N()),lp(0),debug(ddebug)
{
    //cout << " enter in GenericDataFindBoundary"<< endl;
    const int nvE= Element::nv;
    const int nvB = BorderElement::nv;
    const int nvK = bborder ? nvB : nvE;
    const Mesh &Th = *pTh;
    // extract true Border if d != dHat
    // othesize keep all mesh
    RdHat GHat(RdHat::diag(1./(dHat+1)));
    long ncount=0, mcount=0;
    int nv =0;
    //  warning in case of meshL ,  bord is points  => code bborder stupide..
    if(bborder)
        nv =  TrueBorder(Th,(Vertex *)P,delta);
    else
    { //
        for(int k=0; k<Th.nt; ++k)
        {
            const Element& K= Th[k];
            {
                Rd G(K(GHat));
                double l = 0;// 1.5 to be sur .. FH
                for(int i=0; i< Element::nv ;++i)
                    l = max(l, Rd(G,K[i]).norme2()) ;
                delta[nv]=l;
                P[nv].lab= Element::ne*k;//  element and edge
                (Rd &) P[nv++]=G;
                
                
            }
        }
    }
    
    //P.resize(nv); no resize because no copy of vertices ...
    delta.resize(nv);
    lp.resize(nv);
    if(debug>7)  gnuplot("dfb0.gp");
    Vertex * P0= &P[0];
    KN<double> d0(delta);
    for(int i=0; i< nv;++i)
        d0[i]=sqrt(d0[i]);
    delta= 0.;
    Rd Pn, Px;
    Th.BoundingBox(Pn, Px);
    double col=0;
    tree=new EF23::GTree<Vertex> (&P[0], Pn, Px,nv);// build quadtre
    int lvpx=0;
   // cout << " next step in GenericDataFindBoundary"<< endl;
    long  NbQuadTreeBoxSearch=tree->NbQuadTreeBoxSearch,NbVerticesSearch=tree->NbVerticesSearch;
    for(int i=0;i<nv; ++i)
    {
        ncount++;
        if(debug>9)   cout << i << " " << d0[i] << endl;
        long  nbsg=tree->NbQuadTreeBoxSearch,nvsg=tree->NbVerticesSearch;
        int lvp=tree->ListNearestVertex(lp,nv, d0[i]*2.1,P[i]);
        if(debug>9) cout <<i << " " << tree->NbVerticesSearch-nvsg << " " << tree->NbQuadTreeBoxSearch-nbsg << " / " << lvp << " " << d0[i] << endl;
        lvpx=max(lvp,lvpx);
        for(int j=0,k; j<lvp; ++j)
        {
            mcount++;
            k= lp[j]-P0;
            // double dij = Rd(P[i],*lp[j]).norme();
            // warning must be Proof FH.
            delta[k]=max(delta[k],d0[i]+d0[k]);
            delta[i]=max(delta[i],d0[i]+d0[k]);
            
            if(debug>9) cout << k << " "<< delta[k] << ", ";
        }
    }
    if(debug>9)
        for(int i=0;i<nv; ++i)
            cout  << i << " " << d0[i] << " " <<delta[i] << endl;
    if(debug>5)
        gnuplot("dfb1.gp");
    if(verbosity>4)
    {
    cout << " ** BuildDataFindBoundary " << nv << " " << delta.max() << " " << delta.sum()/nv <<  " " << delta.min()  << endl;
    cout << "    count  "<< mcount << " / " <<ncount<< " ratio "
    << (double) mcount / max(ncount,1L) << " max lvp " << lvpx
    << " gtree search box " << tree->NbQuadTreeBoxSearch - NbQuadTreeBoxSearch
    << " v " << tree->NbVerticesSearch-NbVerticesSearch << endl;
    }
}
// Bof Bof pas sur du tout.
/*
 template<typename Mesh>
 void BuildDataFindBoundary<Mesh>() const
 {
 static int count =0;
 if( dfb ==0) {
 dfb=new DataFindBoundary(this);//,count++==0?9:0);
 dfb->debug=0;
 }
 
 }
 */
template class GenericDataFindBoundary<Fem2D::MeshS::GMesh>;
template class GenericDataFindBoundary<Fem2D::Mesh3::GMesh>;
template class GenericDataFindBoundary<Fem2D::Mesh2::GMesh>;
template class GenericDataFindBoundary<Fem2D::Mesh1::GMesh>;
template class GenericDataFindBoundary<Fem2D::MeshL::GMesh>;
