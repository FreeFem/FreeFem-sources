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
#include <iostream> 
using namespace std;
#include "meshtype.h"
#include "SetOfE4.h"

namespace bamg {

 SetOfEdges4::SetOfEdges4(Int4 mmx,Int4 nnx)
   {nx=nnx;
    nbax=mmx;
    NbOfEdges = 0;
    tete= new Int4 [nx];
    Int4 i=nx;
    while (i--)
      tete[i] = -1;// vide 
    Edges =new Int4Edge[nbax];
   }
    
 Int4 SetOfEdges4::find(Int4 ii,Int4 jj)
{ 
  if (tete == 0 ) {
    cerr <<"SetOfEdges4::find \nplus de tete de liste\n";
    MeshError(888);}
  Int4 n = tete[ Abs( ii ) % nx ];
  
  while (n >= 0) 
    if (ii == Edges[n].i && jj == Edges[n].j)
      return n;
    else n = Edges[n].next;
  return -1; // n'existe pas
}

 Int4 SetOfEdges4::add(Int4 ii,Int4 jj)
{
  if (tete == 0 ) {
    cerr << "SetOfEdges4::add\n plus de tete de liste \n" << endl;
    MeshError(888);}
  
  Int4 h;
  Int4 n = tete[ h = Abs( ii ) % nx ];
  while (n >= 0) 
   if (ii == Edges[n].i && jj == Edges[n].j)
            return n;
   else n = Edges[n].next;
  if (nbax <=NbOfEdges ) {
    cerr << " SetOfEdges4::add\noverflow de la pile "  << nbax << " " << NbOfEdges << endl;
    MeshError(888);}
  
   Edges[NbOfEdges].i=ii;
   Edges[NbOfEdges].j=jj;
   Edges[NbOfEdges].next= tete[h];
   tete[h] = NbOfEdges;
   return NbOfEdges ++;
}


}  // end of namespace bamg 
