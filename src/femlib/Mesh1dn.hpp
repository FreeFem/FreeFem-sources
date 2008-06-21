// ORIG-DATE:     Dec 2007
// -*- Mode : c++ -*-
//
// SUMMARY  :  Model  mesh 1d   
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

#ifndef MESH1DN_HPP_
#define MESH1DN_HPP_

namespace Fem2D {
#include "R1.hpp"
}

#include "GenericMesh.hpp" 
#include <cstdlib>
using namespace std;

namespace Fem2D {
  
typedef GenericVertex<R1> Vertex1;

struct DataSeg1  {
  static const int NbOfVertices =2;
  static const int NbOfFaces =0;
  static const int NbOfEdges =1;
  static const int NT =0;
  static const int NbOfAdjElem =NbOfVertices;
  static const int NbOfVertexOnHyperFace =NbOfVertices-1;
  typedef Vertex1 V;
  typedef  V::Rd Rd ;
  static R mesure(  V *  pv[NbOfVertices]) {    
    return pv[1]->x-pv[0]->x;
  } 
  typedef R0 RdHatBord;
  typedef R1 RdHat;
  static RdHat PBord(const int * nvb,const RdHatBord & P)  { return R1(*nvb) ;}  
  //  static const int (* const nvface)[3];// = nvfaceTria  ;
  //static const int (* const nvedge)[2];// = nvedgeTrai;

};


struct DataPoint1  {
  static const int NbOfVertices =1;
  static const int NbOfEdges =0;
  static const int NbOfFaces =0;
  static const int NT =0;
  static const int NbOfAdjElem =1;
  static const int NbOfVertexOnHyperFace =1;
  typedef Vertex1 V;
  typedef  V::Rd Rd;
  static R mesure(  V * pv[NbOfVertices]  ) {    
    return 1.;
  }
  typedef R0 RdHatBord;
  typedef R0 RdHat;
  static RdHat PBord(const int * nvb,const RdHatBord & P)  { return R0() ;}  

};



class Seg1: public GenericElement<DataSeg1>  {
public: 
  Seg1() {}; // constructor empty for array


  R1 H(int i) const { ASSERTION(i>=0 && i <1);
    return (2-i)/mesure();} // heigth 
  
  void Gradlambda(R1 * GradL) const
  {
    GradL[1]= H(1);
    GradL[0]=-GradL[1];
  }
  
  

};

class BoundaryPoint1: public GenericElement<DataPoint1>  {
public: 
  BoundaryPoint1() {}; // constructor empty for array


};


class Mesh1 : public GenericMesh<Seg1,BoundaryPoint1,Vertex1> { 
public:
  Mesh1(const char *); // 
  const Element * Find( R1 P, R1 & Phat,bool & outside,const Element * tstart) const;
 private:
  int load(const string & filename); 
   Mesh1(const Mesh1 &); // pas de construction par copie
   void operator=(const Mesh1 &);// pas affectation par copy 
};
}
#endif
 
