// ********** DO NOT REMOVE THIS BANNER **********
//
// SUMMARY:  Bamg: Bidimensional Anisotrope Mesh Generator
// RELEASE: 0 
// USAGE  : You may copy freely these files and use it for    
//          teaching or research. These or part of these may   
//          not be sold or used for a commercial purpose with- 
//          out our consent : fax (33) 1 39 63 55 14       
//
// AUTHOR:   F. Hecht,    
// ORG    :  INRIA
// E-MAIL :   Frederic.Hecht@Inria.fr   
//
// ORIG-DATE:     Dec 97
#include <cmath>
#include <cstdlib>
#include "error.hpp"
#include <iostream>

#include <limits.h>
//#include <stdio.h>
#include <string.h>
#include <cstdlib>
using namespace std;

#include "RNM.hpp"
#include "rgraph.hpp"

#include "fem.hpp"
using namespace Fem2D;

#define INTER_SEG(a,b,x,y) (((y) > (a)) && ((x) <(b)))
#define ABS(i) ((i)<0 ?-(i) :(i))
#define MAX1(i,j) ((i)>(j) ?(i) :(j))
#define NORM(i1,j1,i2,j2) MAX1(ABS((i1)-(j1)),ABS((i2)-(j2)))

#define IJ(i,j,l) ( ( j & l) ? (( i & l) ? 3 : 2 ) :( ( i & l)? 1 : 0 ))
#define I_IJ(k,l)  (( k&1) ? l : 0)
#define J_IJ(k,l)  (( k&2) ? l : 0)


#ifdef DRAWING

void  FQuadTree::Draw()
{
  QuadTreeBox * pb[ MaxDeep ];
  int  pi[ MaxDeep  ];
  long ii[  MaxDeep ], jj [ MaxDeep];
  int l=0; // level
  QuadTreeBox * b;
  IntQuad hb =  MaxISize;
  if (!root) return ;
  long kkk =0;
  pb[0]=  root;
  pi[0]= root->n>0 ?(int)  root->n : 4  ;;
  ii[0]=jj[0]=0;
  do{    
    b= pb[l];

    while (pi[l]--)
      { 
	if (b->n>0) // Vertex QuadTreeBox none empty
	    { // 
	      for (int k=0;k<b->n;k++)
		{
		  DrawMark(*b->v[k],0.002);

		}
	      break;
	    }
	else // Pointer QuadTreeBox 
	    { 
	      int lll = pi[l];
	      QuadTreeBox *b0=b;
	      
	      if ((b=b->b[lll])) 
		{ 
		  hb >>=1 ; // div by 2
		  long iii = ii[l]+I_IJ(lll,hb);
		  long jjj = jj[l]+J_IJ(lll,hb);
		  
		  pb[++l]=  b;
		  pi[l]= 4;
		  ii[l]= iii;
		  jj[l]= jjj;
		  
		  IMoveTo(ii[l],jj[l]);
		  ILineTo(ii[l]+hb,jj[l]);
		  ILineTo(ii[l]+hb,jj[l]+hb);
		  ILineTo(ii[l]   ,jj[l]+hb);
		  ILineTo(ii[l]   ,jj[l]);
		  
		  
		}
	      else
		{
		  long iii = ii[l]+I_IJ(lll,hb/2);
		  long jjj = jj[l]+J_IJ(lll,hb/2);
		  b=b0;

		  IMoveTo(iii,     jjj);
		  ILineTo(iii+hb/2,jjj+hb/2);
		  IMoveTo(iii+hb/2,jjj);
		  ILineTo(iii     ,jjj+hb/2);

		}
	    }
      }
    hb <<= 1; // mul by 2 
  } while (l--);

}

#endif

Vertex *  FQuadTree::NearestVertex(long xi,long yj)
{
  QuadTreeBox * pb[ MaxDeep ];
  int  pi[ MaxDeep  ];
  long ii[  MaxDeep ], jj [ MaxDeep];
  int l=0; // level
  QuadTreeBox * b;
  IntQuad  h=MaxISize,h0;
  IntQuad hb =  MaxISize;
  long  i0=0,j0=0;
  long  iplus( xi<MaxISize?(xi<0?0:xi):MaxISize-1);
  long  jplus( yj<MaxISize?(yj<0?0:yj):MaxISize-1);
  
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
       int k = IJ(iplus,jplus,hb2);// QuadTreeBox number of size hb2 contening i;j
      QuadTreeBox * b0= b->b[k];
      if ( ( b0 == 0) || (b0->n == 0) ) 
	break; // null box or empty   => break 	    
      NbQuadTreeBoxSearch++;
      b=b0;	
      i0 += I_IJ(k,hb2); // i orign of QuadTreeBox
      j0 += J_IJ(k,hb2); // j orign of QuadTreeBox 
      hb = hb2; 
    }
  
  
  if ( n0 > 0) 
    {  
      for( int k=0;k<n0;k++)
	{
	  I2 i2 =  R2ToI2(b->v[k]);
	  h0 = NORM(iplus,i2.x,jplus,i2.y);
	  if (h0 <h) {
	    h = h0;
	    vn = b->v[k];}
	  NbVerticesSearch++;
	}
      return vn;
    }
  // general case -----
  pb[0]= b;
  pi[0]=b->n>0 ?(int)  b->n : 4  ;
  ii[0]=i0;
  jj[0]=j0;
  h=hb;
  do {    
    b= pb[l];
    while (pi[l]--)
      { 	      
        int k = pi[l];
	
	if (b->n>0) // Vertex QuadTreeBox none empty
	  { 
	    NbVerticesSearch++;
	    I2 i2 =  R2ToI2(b->v[k]);
	    h0 = NORM(iplus,i2.x,jplus,i2.y);
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
	        long iii = ii[l]+I_IJ(k,hb);
	        long jjj = jj[l]+J_IJ(k,hb);
		
		if  (INTER_SEG(iii,iii+hb,iplus-h,iplus+h) && INTER_SEG(jjj,jjj+hb,jplus-h,jplus+h)) 
		  {
		    pb[++l]=  b;
		    pi[l]= b->n>0 ?(int)  b->n : 4  ;
		    ii[l]= iii;
		    jj[l]= jjj;
		    
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



Vertex *  FQuadTree::ToClose(const R2 & v,R seuil,long hx,long hy)
{
  const long i=XtoI(v.x);
  const long j=YtoJ(v.y);
  const R2 X(v);
  R seuil2 = seuil*seuil;
 // const Metric  Mx(v.m);

  QuadTreeBox * pb[ MaxDeep ];
  int  pi[ MaxDeep  ];
  long ii[  MaxDeep ], jj [ MaxDeep];
  int l=0; // level
  QuadTreeBox * b;
  long h=MaxISize;
  long hb =  MaxISize;
  long i0=0,j0=0;
  
  //  Vertex *vn=0;
  
  if (!root->n)
    return 0; // empty tree 
  
  // general case -----
  pb[0]= root;
  pi[0]=  root->n>0 ?(int)  root->n : 4 ;
  ii[0]=i0;
  jj[0]=j0;
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
	    I2 i2 =  R2ToI2(V);
	    if ( ABS(i-i2.x) <hx && ABS(j-i2.y) <hy )
	      {
		R2 XY(X,V);
		R dd;
	      // old code	        if( Mx(XY) + b->v[k]->m(XY) < seuil )
	        if( (dd= (XY,XY) ) < seuil ) // LengthInterpole(Mx(XY), b->v[k]->m(XY)))  < seuil )
		  {
		    //  cout <<  CurrentTh->Number(v) << "is To Close " 
		    // << CurrentTh->Number( b->v[k]) << " l=" <<dd<<endl;
		    return &V; 
		  }
	      }
	  }
	else // Pointer QuadTreeBox 
	  { 
	    QuadTreeBox *b0=b;
	    NbQuadTreeBoxSearch++;
	    if ((b=b->b[k])) 
	      {
		hb >>=1 ; // div by 2
	        long iii = ii[l]+I_IJ(k,hb);
		long jjj = jj[l]+J_IJ(k,hb);
		
		if  (INTER_SEG(iii,iii+hb,i-hx,i+hx) && INTER_SEG(jjj,jjj+hb,j-hy,j+hy)) 
		  {
		    pb[++l]=  b;
		    pi[l]= b->n>0 ?(int)  b->n : 4  ;
		    ii[l]= iii;
		    jj[l]= jjj;
		    
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


void  FQuadTree::Add( Vertex & w)
{
  QuadTreeBox ** pb , *b;
  long i= XtoI(w.x), j=YtoJ(w.y),l=MaxISize;
  pb = &root;
  //    cout << pb << " " << &root << endl;
  while( (b=*pb) && (b->n<0))
    { 
      b->n--;
      l >>= 1;
      pb = &b->b[IJ(i,j,l)];
    }
  if  (b) {      
    if (b->n > 3 &&  b->v[3] == &w) return;
    if (b->n > 2 &&  b->v[2] == &w) return;
    if (b->n > 1 &&  b->v[1] == &w) return;
    if (b->n > 0 &&  b->v[0] == &w) return;
  }
  throwassert(l);
  while ((b= *pb) && (b->n == 4)) // the QuadTreeBox is full
    { 
      Vertex *v4[4]; // copy of the QuadTreeBox vertices
      
      v4[0]= b->v[0];
      v4[1]= b->v[1];
      v4[2]= b->v[2];
      v4[3]= b->v[3];
      b->n = -b->n; // mark is pointer QuadTreeBox
      b->b[0]=b->b[1]=b->b[2]=b->b[3]=0; // set empty QuadTreeBox ptr
      l >>= 1;    // div the size by 2
      for (int k=0;k<4;k++) // for the 4 vertices find the sub QuadTreeBox ij
	{ 
	  int ij;
	  QuadTreeBox * bb =  b->b[ij=IJ(XtoI(v4[k]->x),YtoJ(v4[k]->y),l)];
	  if (!bb) 
	    bb=b->b[ij]=NewQuadTreeBox(); // alloc the QuadTreeBox 
	  //    cout << bb << " " << k << " "  << ij <<  endl;
	  bb->v[bb->n++] = v4[k];
	}
      pb = &b->b[IJ(i,j,l)];
    }
  if (!(b = *pb))
    b=*pb= NewQuadTreeBox(); //  alloc the QuadTreeBox 
  //   cout << b << " " << b->n << endl;
  b->v[b->n++]=&w; // we add the vertex 
  NbVertices++;    
}

FQuadTree::FQuadTree(Mesh * t,R2 Pmin,R2 Pmax,long nbv) : 
  th(t),
  lenStorageQuadTreeBox(t->nv/8+100),
  NbQuadTreeBox(0),
  NbVertices(0),
  NbQuadTreeBoxSearch(0),
  NbVerticesSearch(0),
  cMin(Pmin-(Pmax-Pmin)/2),
  cMax(Pmax+(Pmax-Pmin)/2),
  coef( MaxISize/Norme_infty(cMax-cMin))
  
{ 
  if (nbv == -1) nbv = t->nv;
  sb =new StorageQuadTreeBox(lenStorageQuadTreeBox);
  root=NewQuadTreeBox();
//  throwassert( MaxISize > MaxICoor);
  if (t)
  for (long i=0;i<nbv;i++) 
    Add(t->vertices[i]);
#ifdef DRAWING1
  Draw();
#endif
}

FQuadTree::FQuadTree() : 
  th(0),
  lenStorageQuadTreeBox(100),
  NbQuadTreeBox(0),
  NbVertices(0),
  NbQuadTreeBoxSearch(0),
  NbVerticesSearch(0),
  coef(0),cMin(0,0),cMax(0,0)
{
  sb =new StorageQuadTreeBox(lenStorageQuadTreeBox);
  root=NewQuadTreeBox();
}
FQuadTree::StorageQuadTreeBox::StorageQuadTreeBox(long ll,StorageQuadTreeBox *nn)
{
  len = ll;
  n = nn;
  b = new QuadTreeBox[ll];
  for (int i = 0; i <ll;i++)
    b[i].n =0,b[i].b[0]=b[i].b[1]=b[i].b[2]=b[i].b[3]=0;
  bc =b;
  be = b +ll;
  throwassert(b);
}

FQuadTree::~FQuadTree()
{
  delete sb; 
}

ostream& operator <<(ostream& f, const  FQuadTree & qt)
{ 
  f << " the quadtree "  << endl;
  f << " NbQuadTreeBox = " << qt.NbQuadTreeBox 
    << " Nb Vertices = " <<  qt.NbVertices << endl;
  f << " NbQuadTreeBoxSearch " << qt.NbQuadTreeBoxSearch  
    << " NbVerticesSearch " << qt.NbVerticesSearch << endl;
  f << " SizeOf QuadTree" << qt.SizeOf() << endl;
  //     return  dump(f,*qt.root);
  return  f;
}

Vertex *  FQuadTree::NearestVertexWithNormal(const R2 &P)//(long xi,long yj)
{
  long xi(XtoI(P.x)),yj(YtoJ(P.y));
  QuadTreeBox * pb[ MaxDeep ];
  int  pi[ MaxDeep  ];
  long ii[  MaxDeep ], jj [ MaxDeep];
  int l; // level
  QuadTreeBox * b;
  IntQuad  h=MaxISize,h0;
  IntQuad hb =  MaxISize;
  long  i0=0,j0=0;
  long  iplus( xi<MaxISize?(xi<0?0:xi):MaxISize-1);
  long  jplus( yj<MaxISize?(yj<0?0:yj):MaxISize-1);
  
  Vertex *vn=0;
 // Vertex *vc;  // the current vertex
  
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
      int k = IJ(iplus,jplus,hb2);// QuadTreeBox number of size hb2 contening i;j
      QuadTreeBox * b0= b->b[k];
      if ( ( b0 == 0) || (b0->n == 0) ) 
	break; // null box or empty   => break 	    
      NbQuadTreeBoxSearch++;
      b=b0;	
      i0 += I_IJ(k,hb2); // i orign of QuadTreeBox
      j0 += J_IJ(k,hb2); // j orign of QuadTreeBox 
      hb = hb2; 
    }
  
  
  if ( n0 > 0) 
    {  
      for(int k=0;k<n0;k++)
	{
	  Vertex * v=b->v[k];
	  if (v->ninside(P)) {
	   I2 i2 =  R2ToI2(v);
	  //   try if is in the right sens -- 
	   h0 = NORM(iplus,i2.x,jplus,i2.y);
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
  pi[0]= b->n>0 ?(int)  b->n : 4  ;
  ii[0]=i0;
  jj[0]=j0;
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
	     I2 i2 =  R2ToI2(v);
	     // if good sens when try -- 
	     h0 = NORM(iplus,i2.x,jplus,i2.y);
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
		long iii = ii[l]+I_IJ(k,hb);
		long jjj = jj[l]+J_IJ(k,hb);
		
		if  (INTER_SEG(iii,iii+hb,iplus-h,iplus+h) && INTER_SEG(jjj,jjj+hb,jplus-h,jplus+h)) 
		  {
		    pb[++l]=  b;
		    pi[l]= b->n>0 ?(int)  b->n : 4  ;
		    ii[l]= iii;
		    jj[l]= jjj;
		    
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
   i0=0;
   j0=0;
   l=0;
   pb[0]= b;
   pi[0]= b->n>0 ?(int)  b->n : 4  ;
   ii[0]=i0;
   jj[0]=j0;
   
     goto L1;
   }
  return vn;
}

