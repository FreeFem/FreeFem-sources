// ********** DO NOT REMOVE THIS BANNER **********
//
// SUMMARY: Bamg: Bidimensional Anisotrope Mesh Generator
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

#ifndef _SetOfEdge4_h
#define _SetOfEdge4_h

namespace bamg {

class SetOfEdges4 ;
class Int4Edge{
friend class SetOfEdges4;
public:
  Int4 i,j;
  Int4 next; 
};

class SetOfEdges4 {
  Int4 nx,nbax,NbOfEdges;
  Int4 * tete; 
  Int4Edge * Edges;

public:
  SetOfEdges4(Int4 ,Int4);// nb Edges mx , nb de sommet 
  ~SetOfEdges4() {// cout << " delete SetofArete " << endl ;
  delete [] tete; delete [] Edges;}
   Int4 add (Int4 ii,Int4 jj);
  Int4 addtrie (Int4 ii,Int4 jj) {return ii <=jj ? add (ii,jj)  : add (jj,ii) ;}
  Int4  nb(){return NbOfEdges;}
  Int4 find (Int4 ii,Int4 jj);
  Int4 findtrie (Int4 ii,Int4 jj) {return ii <=jj ? find (ii,jj)  : find (jj,ii) ;}
  // inline void close();
  Int4 i(Int4 k){return Edges[k].i;}
  Int4 j(Int4 k){return Edges[k].j;}
  Int4 newarete(Int4 k){return NbOfEdges == k+1;}
  Int4Edge & operator[](Int4 k){return  Edges[k];}
};
}

#endif 
