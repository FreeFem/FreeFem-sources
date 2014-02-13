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

namespace bamg {

const int MaxDeep = 30;
typedef  long  IntQuad;
const IntQuad MaxISize = ( 1L << MaxDeep);


class Triangles;
class Vertex;

class QuadTree {
 public:

  class QuadTreeBox { 
  public:

    long n; // if n < 4 => Vertex else =>  QuadTreeBox;
    union {
      QuadTreeBox *b[4];
      Vertex * v[4];
    };
    

  }; // end class QuadTreeBox  /////////////////

  class StorageQuadTreeBox {
  public:
    QuadTreeBox *b,*bc,*be;
    long len;
    StorageQuadTreeBox *n; // next StorageQuadTreeBox
    StorageQuadTreeBox(long ,StorageQuadTreeBox * =0);
    ~StorageQuadTreeBox()
    { //cout <<  "~StorageQuadTreeBox " << this << " n " << n << " b " << b << endl;
      if(n) delete n;
      delete [] b;
    }
    long  SizeOf() const {
      return len*sizeof(QuadTreeBox)+sizeof(StorageQuadTreeBox)+ (n?n->SizeOf():0);
    }
  }; // end class  StorageQuadTreeBox 
  
  StorageQuadTreeBox * sb;
  
  
  long  lenStorageQuadTreeBox;

public:
  QuadTreeBox * root;
  Triangles *th;
  long NbQuadTreeBox,NbVertices;
  long NbQuadTreeBoxSearch,NbVerticesSearch;
  Vertex * NearestVertex(Icoor1 i,Icoor1 j);
  Vertex *  NearestVertexWithNormal(Icoor1 i,Icoor1 j); // new version  
  Vertex * ToClose(Vertex & ,Real8 ,Icoor1,Icoor1);
  long SizeOf() const {return sizeof(QuadTree)+sb->SizeOf();}

#ifdef DRAWING
    void Draw();
#endif

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
  ~QuadTree();
  QuadTree(Triangles * t,long nbv=-1);
  QuadTree();
  friend ostream& operator <<(ostream& f, const  QuadTree & qt);


	
};
}
//#undef IJ
