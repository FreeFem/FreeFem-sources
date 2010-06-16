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
