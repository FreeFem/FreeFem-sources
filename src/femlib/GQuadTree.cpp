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
// E-MAIL :   Frederic.Hecht@upmc.fr   
//
// ORIG-DATE:    jan 2008
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
  Vertex *  GTree<Vertex>::NearestVertex(Zd xyi)//long xi,long yj)
  {
    QuadTreeBox * pb[ MaxDeep ];
    int  pi[ MaxDeep  ];
    Zd pp[  MaxDeep ];
    int l=0; // level
    QuadTreeBox * b;
    IntQuad  h=MaxISize,h0;
    IntQuad hb =  MaxISize;
    Zd   p0;
    Zd  plus(xyi) ; plus.Bound();
    // xi<MaxISize?(xi<0?0:xi):MaxISize-1,yj<MaxISize?(yj<0?0:yj):MaxISize-1);
    
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
      b=b0;
      p0.Add(k,hb2);	
      hb = hb2; 
      }
    
    
  if ( n0 > 0) 
    {  
      for( int k=0;k<n0;k++)
	{
	  Zd i2 =  VtoZd(b->v[k]);
	  h0 = Zd(i2,plus).norm();//NORM(iplus,i2.x,jplus,i2.y);
	  if (h0 <h) {
	    h = h0;
	    vn = b->v[k];}
	  NbVerticesSearch++;
	}
      return vn;
    }
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
	    h0 = Zd(i2,plus).norm();//  NORM(iplus,i2.x,jplus,i2.y);
	    if (h0 <h) 
	      {
		h = h0;
		vn = b->v[k];
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
  Vertex *  GTree<Vertex>::ToClose(const Rd & v,R seuil,Zd H)
  {
    const Rd X(v);
    const Zd p(RdtoZd(v));
    R seuil2 = seuil*seuil;
    // const Metric  Mx(v.m);
    
    QuadTreeBox * pb[ MaxDeep ];
    int  pi[ MaxDeep  ];
    Zd pp[  MaxDeep ];
    
    int l=0; // level
    QuadTreeBox * b;
    long h=MaxISize;
    long hb =  MaxISize;
    Zd p0;
    
    if (!root->n)
      return 0; // empty tree 
    
    // general case -----
    pb[0]= root;
    pi[0]=  root->n>0 ?(int)  root->n : N ;
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
	      Vertex & V(*b->v[k]);
	      Zd i2 =  VtoZd(V);
	      if ( Zd(i2,p).less(H) )
		{
		  Rd XY(X,V);
		  R dd;
		  if( (dd= (XY,XY) ) < seuil2 ) // LengthInterpole(Mx(XY), b->v[k]->m(XY)))  < seuil )
		    {// cout << dd << " " << XY << " ";
		      return &V; }
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
		  if (ppp.interseg(p,hb,H))
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
    
    return 0;
  }
  
  template<class Vertex>
  void  GTree<Vertex>::Add( Vertex & w)
  {
    QuadTreeBox ** pb , *b;
    Zd p(VtoZd(w));
    long l=MaxISize;
    pb = &root;
    //cout << "add : "<< w << " : "<< p  << ", " << pb << " " << &root << endl;
    while( (b=*pb) && (b->n<0))
      { 
	b->n--;
	l >>= 1;
	pb = &b->b[p.Case(l)];
      }
    if  (b) {
      for(int i=N-1;i>=0;--i)
	if (b->n > i &&  b->v[i] == &w) return;
    }
    assert(l);
    while ((b= *pb) && (b->n == N)) // the QuadTreeBox is full
      { 
	//cout << " b = " << b << b->n << "  " << l << endl;
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
	    //cout << "ij= "<< ij <<  " " << VtoZd(v4[k])<< endl;
	    if (!bb) 
	      bb=b->b[ij]=NewQuadTreeBox(); // alloc the QuadTreeBox 
	    //cout << bb << " " << k << " "  << ij <<  endl;
	     bb->v[bb->n++] = v4[k];
	  }
	pb = &b->b[p.Case(l)];
      }
    if (!(b = *pb))
      { //cout << "Qbox \n";
	b=*pb= NewQuadTreeBox(); //  alloc the QuadTreeBox 
      }
    //   cout << b << " " << b->n << endl;
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
  if(verbosity>5)
    cout << "  GTree: box: "<<  Pmin << " " << Pmax << " " << cMin << " "<< cMax << " nbv : " << nbv <<endl;
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
	b=b0;	
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
  
  
  inline    int nRand(int n) {
    return  (rand()*n)/(RAND_MAX+1);
  }
  
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
	    else if(i>k[im]) i0=im+1;
	    else return im;
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
    int it,j;
    const int mxbord=100;
    int kbord[mxbord];
    int nbord=0;
    if ( tstart )
      it =  Th(tstart);
    else  
      {  
	const Vertex * v=quadtree->NearestVertexWithNormal(P);
	if (!v) 
	  { 
	  v=quadtree->NearestVertex(P);
	  assert(v);
	  }
	it=Th.Contening(v);
      }
    
    //     int itdeb=it;     
    //     int count=0;
    //     L1: 
  outside=true; 
  //int its=it;
  //dPdP	int iib=-1;//,iit=-1;
  R dP=DBL_MAX;
  Rd PPhat;
  int k=0;    
  Mesh::kfind++;
  while (1)
    { 
      // cout << "it " << it <<endl;
      const Element & K(Th[it]);
      Mesh::kthrough++;
      assert(k++<1000);
      int kk,n=0,nl[nkv];
      R l[nkv+1];
      CoorBary(K,P,l);
      
      R eps =  -K.mesure() *1e-10;
      for(int i=0;i<nkv;++i)
	if (l[i] < eps) nl[n++]=i;
      
      if (n==0)
	{  // interior => return
	  outside=false; 
	  Phat=Rd(l +1);
	  // cout << Phat <<endl;
	  return &K;
	}
      
      kk=n==1 ? 0 : nRand(n);
      j= nl[ kk ];
      int itt =  Th.ElementAdj(it,j);
      if(itt!=it && itt >=0)  
	{
	  dP=DBL_MAX;
	  it=itt;
	  continue;
	}  
      int  inkbord=find5(it,kbord,nbord);
      if(inkbord<0 )  
	{
	  assert(nbord < mxbord);
	  kbord[nbord++]=it;
	  if(nbord>=5)  HeapSort(kbord,nbord);
	}
      if(verbosity>100)
	cout << " bord "<< it<< "   nbf < 0 : " <<n << " (inb) " << inkbord << " nfb" << nbord<<endl;
      if ( n!=1 )  // on est sur le bord, mais plusieurs face <0 => on test les autre
	{  // 1) existe t'il un adj interne
	  int nn=0,ii;
	  int nadj[d+1],ni[d+1];
	  for(int i=0;i<nkv;++i)
	    if (l[i] < eps && (itt=Th.ElementAdj(it,ii=i)) != it && itt && find5(itt,kbord,nbord) < -1 ) 
	      ni[nn++]=i,nadj[i]=itt;
	  if(verbosity>100)
	    cout << " nn : "<< nn << endl;
	  if (nn>0)
	    {
	      j=nadj[nRand(nn)];
	      it=nadj[j];
	      dP=DBL_MAX;
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
	if(verbosity>100)
	  cout << " sortie on bord "<<endl;
	R s=0;
	for(int i=0;i<=nkv;++i)
	  s += (l[i]<0.) ?  (l[i]=0) : l[i];
	for(int i=0;i<=nkv;++i)
	  l[i]/=s;
	Phat=Rd(l +1);
	outside=true;
	return &Th[it] ;
      }		    
    }
}
// Instantiation du manuel des templates

  template class GTree<Vertex2>;
  template class GTree<Vertex3>;
  template class GTree<Vertex1>;
  typedef Mesh3::GMesh GMesh3;
  typedef Mesh2::GMesh GMesh2;
  typedef Mesh1::GMesh GMesh1;
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
