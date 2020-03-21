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
	    h0 = Zd(i2,plus).norm2();//  NORM(iplus,i2.x,jplus,i2.y);
	    if (h0 <h2)
	      {
		h2 = h0;
                h =Zd(i2,plus).norm();
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
  int k=0;
    const int nReStartMax = 1 << Rd::d; //  Nb vertex of the d-cube/
    int nReStart = 0;
    int itstart[100],itout[100],kstart=0;
    Rd Pout[100];
  int it,j,it00;
  const int mxbord=1000;
  int kbord[mxbord+1];
  int nbord=0;

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
    cout << "    " << nReStart << " tstart= " << tstart << " , it=" << it << " P="<< P << endl;
  outside=true; 
  Mesh::kfind++;
  while (1)
    { 
      //if(verbosity>199) cout << "it " << it <<endl;
      const Element & K(Th[it]);
      Mesh::kthrough++;
      assert(k++<1000);
      int kk,n=0,nl[nkv];
      R l[nkv];
      for(int iii=0; iii<nkv; iii++)
	l[iii]=0.;
      
      CoorBary(K,P,l);
      
      // CoorBary :: donner entre 0 et 1
      // Pb si K.mesure*1.e-10 precision machine ==> bug
      
      // avant:
      // R eps =  -K.mesure()*1e-10;
      R eps = -1e-10;
      for(int i=0;i<nkv;++i)
	if( l[i] < eps){
	  nl[n++]=i;
	}
      if(verbosity>200){
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
  CoorBary(K,P,l);
  
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
