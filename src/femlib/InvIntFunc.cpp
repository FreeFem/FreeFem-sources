// used by splitsimplex.cpp 
//  to inverse numering ... 
//   F. Hecht 
// ORIG-DATE:     fev 2009
// -*- Mode : c++ -*-
//
// SUMMARY  :  Model  mesh 2d   
// USAGE    : LGPL      
// ORG      : LJLL Universite Pierre et Marie Curie, Paris,  FRANCE 
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
 ref:ANR-07-CIS7-002-01  */

inline int InvIntFunction(int l)
// calcul de inverse de la fonction F 
{
  // inverse la function F 
  int i=0,j,k=F_1(l);
  int Fi=F(i),Fj,Fk=F(k);
  assert(l<=Fk); 
  while (1)
    {
      j = (i+k)/2;
      if(j==i) break;
      Fj=F(j);      
      // cout << i<< j << k << " " << (l < Fj)  << " : "; 
      if( l < Fj ) { k=j; Fk=Fj;}
      else if ( l > Fj  ) { i=j; Fi=Fj;}
      else { i=j;}      
      // cout <<  "                ** " << l << " : " << i<< " "<< j << " "<< k << " : Fi  " << Fi << " " << Fj << " "<< Fk << endl;
      
    }


  if( Fk==l) i=k;
  //  cout << "           i =" << i << " l=  " << l << " in  [  " << F(i) << ", " << F(i+1) << "[ " << endl;
  assert( (F(i) <= l)  && (l < F(i+1)  )  ); 
  return i; 
}
