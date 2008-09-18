// ORIG-DATE:     Dec 2007
// -*- Mode : c++ -*-
//
// SUMMARY  :  Model  mesh 2d   
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

#ifndef MESH2DN_HPP_
#define MESH2DN_HPP_

 namespace Fem2D {
#include "R2.hpp"
 }

#include "GenericMesh.hpp" 
#include <cstdlib>
using namespace std;

 namespace Fem2D {


  
typedef GenericVertex<R2> Vertex2;

struct DataTriangle2  {
  static const int NbOfVertices =3;
  static const int NbOfFaces =1;
  static const int NbOfEdges =3;
  static const int NT =0;
  static const int NbOfAdjElem =NbOfVertices;
  static const int NbOfVertexOnHyperFace =NbOfVertices-1;
  typedef Vertex2 V;
  typedef  V::Rd Rd ;
  static R mesure(  V *  pv[NbOfVertices]) {    
    return det(*pv[0],*pv[1],*pv[2])*0.5;
  } 
  typedef R2 RdHat;
  typedef R1 RdHatBord;
  static RdHat PBord(const int * nvb,const RdHatBord & P)  { return RdHat::KHat[nvb[0]]*(1-P.x)+R2::KHat[nvb[1]]*(P.x) ;}  

  //  static const int (* const nvface)[3];// = nvfaceTria  ;
  //static const int (* const nvedge)[2];// = nvedgeTrai;

};


struct DataSeg2  {
  static const int NbOfVertices =2;
  static const int NbOfEdges =1;
  static const int NbOfFaces =0;
  static const int NT =0;
  static const int NbOfAdjElem =NbOfVertices;
  static const int NbOfVertexOnHyperFace =NbOfVertices-1;
  typedef Vertex2 V;
  typedef  V::Rd Rd;
  static R mesure(  V *  pv[NbOfVertices]) {    
    return R2(*pv[0],*pv[1]).norme();
  }
  typedef R1 RdHat;
  typedef R0 RdHatBord;
  static RdHat PBord(const int * nvb,const RdHatBord &P)  { return RdHat(*nvb) ;}  

  //static const int (* const nvface)[3];// = nvfaceSeg ;
  //static const int (* const nvedge)[2];//  = nvedgeSeg;

};



class Triangle2: public GenericElement<DataTriangle2>  {
public: 
  Triangle2() {}; // constructor empty for array


  R2 H(int i) const { ASSERTION(i>=0 && i <3);
    R2 E=Edge(i);return E.perp()/(2.*this->mesure());} // heigth 
  
  void Gradlambda(R2 * GradL) const
  {
    GradL[1]= H(1);
    GradL[2]= H(2);
    GradL[0]=-GradL[1]-GradL[2];
  }
  
  

};

class BoundaryEdge2: public GenericElement<DataSeg2>  {
public: 
  BoundaryEdge2() {}; // constructor empty for array


};


class Mesh2 : public GenericMesh<Triangle2,BoundaryEdge2,Vertex2> { 
public:
  Mesh2(const char *); // 
  Mesh2(int nnv, int nnt, int nnbe, Vertex2 *vv, Triangle2 *tt, BoundaryEdge2 *bb); 
  const Element * Find( R2 P, R2 & Phat,bool & outside,const Element * tstart) const;
  int Save(const string & filename);
  //int Popen(const FILE *namestream);
 private:
  int load(const string & filename); 
   Mesh2(const Mesh2 &); // pas de construction par copie
   void operator=(const Mesh2 &);// pas affectation par copy 
};
 }

#endif
 
