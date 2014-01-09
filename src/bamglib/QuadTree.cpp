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

#include <limits.h>
//#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "Meshio.h"
#include "Mesh2.h"
#include "QuadTree.h"

namespace bamg {

#define INTER_SEG(a,b,x,y) (((y) > (a)) && ((x) <(b)))
#define ABS(i) ((i)<0 ?-(i) :(i))
#define MAX1(i,j) ((i)>(j) ?(i) :(j))
#define NORM(i1,j1,i2,j2) MAX1(ABS((i1)-(j1)),ABS((i2)-(j2)))

#define IJ(i,j,l) ( ( j & l) ? (( i & l) ? 3 : 2 ) :( ( i & l)? 1 : 0 ))
#define I_IJ(k,l)  (( k&1) ? l : 0)
#define J_IJ(k,l)  (( k&2) ? l : 0)


#ifdef DRAWING

void  QuadTree::Draw()
{
  QuadTreeBox * pb[ MaxDeep ];
  int  pi[ MaxDeep  ];
  Icoor1 ii[  MaxDeep ], jj [ MaxDeep];
  register int l=0; // level
  register QuadTreeBox * b;
  IntQuad hb =  MaxISize;
  if(!root) return;
  //  Int4 kkk =0;
  pb[0]=  root;
  pi[0]= root->n>0 ?(int)  root->n : 4  ;
  ii[0]=jj[0]=0;
  do{    
    b= pb[l];

    while (pi[l]--)
      { 
	if (b->n>0) // Vertex QuadTreeBox none empty
	    { // 
	      for (int k=0;k<b->n;k++)
		{
		  I2 i2 =  b->v[k]->i;
		  IMoveTo(i2.x,i2.y+50);
		  ILineTo(i2.x,i2.y-50);
		  IMoveTo(i2.x+50,i2.y);
		  ILineTo(i2.x-50,i2.y);

		  assert(ii[l] <= i2.x);
		  assert(jj[l] <= i2.y);
		  assert(ii[l] +hb > i2.x);
		  assert(jj[l] +hb > i2.y);

		}
	      break;
	    }
	else // Pointer QuadTreeBox 
	    { 
	      register int lll = pi[l];
	      register QuadTreeBox *b0=b;
	      
	      if ((b=b->b[lll])) 
		{ 
		  hb >>=1 ; // div by 2
		  register Icoor1 iii = ii[l]+I_IJ(lll,hb);
		  register Icoor1 jjj = jj[l]+J_IJ(lll,hb);
		  
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
		  register Icoor1 iii = ii[l]+I_IJ(lll,hb/2);
		  register Icoor1 jjj = jj[l]+J_IJ(lll,hb/2);
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

Vertex *  QuadTree::NearestVertex(Icoor1 i,Icoor1 j)
{
  QuadTreeBox * pb[ MaxDeep ];
  int  pi[ MaxDeep  ];
  Icoor1 ii[  MaxDeep ], jj [ MaxDeep];
  register int l=0; // level
  register QuadTreeBox * b;
  IntQuad  h=MaxISize,h0;
  IntQuad hb =  MaxISize;
  Icoor1  i0=0,j0=0;
  Icoor1  iplus( i<MaxISize?(i<0?0:i):MaxISize-1);
  Icoor1  jplus( j<MaxISize?(j<0?0:j):MaxISize-1);
  
  Vertex *vn=0;
  
  // init for optimisation ---
  b = root;
  register Int4  n0;
  if (!root->n)
    return vn; // empty tree 
  
  while( (n0 = b->n) < 0) 
    {
      // search the non empty 
      // QuadTreeBox containing  the point (i,j)
      register Icoor1 hb2 = hb >> 1 ;
      register  int k = IJ(iplus,jplus,hb2);// QuadTreeBox number of size hb2 contening i;j
      register QuadTreeBox * b0= b->b[k];
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
      for(register int k=0;k<n0;k++)
	{
	  I2 i2 =  b->v[k]->i;
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
	register int k = pi[l];
	
	if (b->n>0) // Vertex QuadTreeBox none empty
	  { 
	    NbVerticesSearch++;
	    I2 i2 =  b->v[k]->i;
	    h0 = NORM(iplus,i2.x,jplus,i2.y);
	    if (h0 <h) 
	      {
		h = h0;
		vn = b->v[k];
	      }
	  }
	else // Pointer QuadTreeBox 
	  { 
	    register QuadTreeBox *b0=b;
	    NbQuadTreeBoxSearch++;
	    if ((b=b->b[k])) 
	      {
		hb >>=1 ; // div by 2
		register Icoor1 iii = ii[l]+I_IJ(k,hb);
		register Icoor1 jjj = jj[l]+J_IJ(k,hb);
		
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



Vertex *   QuadTree::ToClose(Vertex & v,Real8 seuil,Icoor1 hx,Icoor1 hy)
{
  const Icoor1 i=v.i.x;
  const Icoor1 j=v.i.y;
  const R2 X(v.r);
  const Metric  Mx(v.m);

  QuadTreeBox * pb[ MaxDeep ];
  int  pi[ MaxDeep  ];
  Icoor1 ii[  MaxDeep ], jj [ MaxDeep];
  register int l=0; // level
  register QuadTreeBox * b;
  Icoor1 h=MaxISize;
  Icoor1 hb =  MaxISize;
  Icoor1 i0=0,j0=0;
  
  //  Vertex *vn=0;
  
  if (!root->n)
    return 0; // empty tree 
  
  // general case -----
  pb[0]=root;
  pi[0]=root->n>0 ?(int)  root->n : 4  ;
  ii[0]=i0;
  jj[0]=j0;
  h=hb;
  do {    
    b= pb[l];
    while (pi[l]--)
      { 	      
	register int k = pi[l];
	
	if (b->n>0) // Vertex QuadTreeBox none empty
	  { 
	    NbVerticesSearch++;
	    I2 i2 =  b->v[k]->i;
	    if ( ABS(i-i2.x) <hx && ABS(j-i2.y) <hy )
	      {
		R2 XY(X,b->v[k]->r);
		Real8 dd;
	      // old code	        if( Mx(XY) + b->v[k]->m(XY) < seuil )
	        if( (dd= LengthInterpole(Mx(XY), b->v[k]->m(XY)))  < seuil )
		  {
		    //  cout <<  CurrentTh->Number(v) << "is To Close " 
		    // << CurrentTh->Number( b->v[k]) << " l=" <<dd<<endl;
		    return b->v[k]; 
		  }
	      }
	  }
	else // Pointer QuadTreeBox 
	  { 
	    register QuadTreeBox *b0=b;
	    NbQuadTreeBoxSearch++;
	    if ((b=b->b[k]))
	      {
		hb >>=1 ; // div by 2
		register long iii = ii[l]+I_IJ(k,hb);
		register long jjj = jj[l]+J_IJ(k,hb);
		
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


void  QuadTree::Add( Vertex & w)
{
  QuadTreeBox ** pb , *b;
  register long i=w.i.x, j=w.i.y,l=MaxISize;
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
  assert(l);
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
      for (register int k=0;k<4;k++) // for the 4 vertices find the sub QuadTreeBox ij
	{ 
	  register int ij;
	  register QuadTreeBox * bb =  b->b[ij=IJ(v4[k]->i.x,v4[k]->i.y,l)];
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

QuadTree::QuadTree(Triangles * t,long nbv) : 
  lenStorageQuadTreeBox(t->nbvx/8+10),
  th(t),
  NbQuadTreeBox(0),
  NbVertices(0),
  NbQuadTreeBoxSearch(0),
  NbVerticesSearch(0)
{ 
  if (nbv == -1) nbv = t->nbv;
  sb =new StorageQuadTreeBox(lenStorageQuadTreeBox);
  root=NewQuadTreeBox();
  assert( MaxISize > MaxICoor);
  for (Int4 i=0;i<nbv;i++) 
    Add(t->vertices[i]);
#ifdef DRAWING1
  Draw();
#endif
}

QuadTree::QuadTree() : 
  lenStorageQuadTreeBox(100),
  th(0),
  NbQuadTreeBox(0),
  NbVertices(0),
  NbQuadTreeBoxSearch(0),
  NbVerticesSearch(0)
{
  sb =new StorageQuadTreeBox(lenStorageQuadTreeBox);
  root=NewQuadTreeBox();
}
QuadTree::StorageQuadTreeBox::StorageQuadTreeBox(long ll,StorageQuadTreeBox *nn)
{
  len = ll;
  n = nn;
  b = new QuadTreeBox[ll];
  for (int i = 0; i <ll;i++)
    b[i].n =0,b[i].b[0]=b[i].b[1]=b[i].b[2]=b[i].b[3]=0;
  bc =b;
  be = b +ll;
  assert(b);
}

QuadTree::~QuadTree()
{
  delete sb; 
  root=0;
}

ostream& operator <<(ostream& f, const  QuadTree & qt)
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

Vertex *  QuadTree::NearestVertexWithNormal(Icoor1 i,Icoor1 j)
{
  QuadTreeBox * pb[ MaxDeep ];
  int  pi[ MaxDeep  ];
  Icoor1 ii[  MaxDeep ], jj [ MaxDeep];
  int l; // level
  QuadTreeBox * b;
  IntQuad  h=MaxISize,h0;
  IntQuad hb =  MaxISize;
  Icoor1  i0=0,j0=0;
  Icoor1  iplus( i<MaxISize?(i<0?0:i):MaxISize-1);
  Icoor1  jplus( j<MaxISize?(j<0?0:j):MaxISize-1);
  
  Vertex *vn=0;
  
  // init for optimisation ---
  b = root;
  register Int4  n0;
  if (!root->n)
    return vn; // empty tree 
  
  while( (n0 = b->n) < 0) 
    {
      // search the non empty 
      // QuadTreeBox containing  the point (i,j)
      register Icoor1 hb2 = hb >> 1 ;
      register  int k = IJ(iplus,jplus,hb2);// QuadTreeBox number of size hb2 contening i;j
      register QuadTreeBox * b0= b->b[k];
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
      for(register int k=0;k<n0;k++)
	{
	  I2 i2 =  b->v[k]->i;
	  //   try if is in the right sens -- 
	  h0 = NORM(iplus,i2.x,jplus,i2.y);
	  if (h0 <h) {
	    h = h0;
	    vn = b->v[k];}
	  NbVerticesSearch++;
	}
	if (vn) return vn; 
    }
  // general case -----
  // INITIALISATION OF THE HEAP 
  l =0; // level 
  pb[0]= b;
  pi[0]=b->n>0 ?(int)  b->n : 4  ;
  ii[0]=i0;
  jj[0]=j0;
  h=hb;
  do {   // walk on the tree  
    b= pb[l];
    while (pi[l]--) // loop on 4 element of the box
      { 	      
       int k = pi[l];
	
	if (b->n>0) // Vertex QuadTreeBox none empty
	  { 
	    NbVerticesSearch++;
	    I2 i2 =  b->v[k]->i;
	    // if good sens when try -- 
	    
	    h0 = NORM(iplus,i2.x,jplus,i2.y);
	    if (h0 <h) 
	      {
		h = h0;
		vn = b->v[k];
	      }
	  }
	else // Pointer QuadTreeBox 
	  { 
	    register QuadTreeBox *b0=b;
	    NbQuadTreeBoxSearch++;
	    if ((b=b->b[k])) 
	      {
		hb >>=1 ; // div by 2
		register Icoor1 iii = ii[l]+I_IJ(k,hb);
		register Icoor1 jjj = jj[l]+J_IJ(k,hb);
		
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


}  // end of namespace bamg

