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


#ifndef MESH3DN_HPP_
#define MESH3DN_HPP_



// definition R
#include <cstdlib>
namespace Fem2D {

#include "R3.hpp"
}
 
 using namespace ::std;
#include "GenericMesh.hpp" 

namespace Fem2D {
  
typedef GenericVertex<R3> Vertex3;

struct DataTriangle3  {
  static const int NbOfVertices =3;
  static const int NbOfEdges =3;
  static const int NbOfFaces =1;
  static const int NT =0;
  static const int NbOfAdjElem =NbOfVertices;
  static const int NbOfVertexOnHyperFace =NbOfVertices-1;
  typedef Vertex3 V;
  typedef  V::Rd Rd ;
  typedef R2 RdHat;
  typedef R1 RdHatBord;
  static RdHat PBord(const int * nvb,const RdHatBord &P)  { 
  return RdHat::KHat[nvb[0]]*(1-P.x)+R2::KHat[nvb[1]]*(P.x) ;}  

  static R mesure(  V *  pv[NbOfVertices]) {    
    return (R3(*pv[0],*pv[1])^R3(*pv[0],*pv[2])).norme()*0.5;
  }
};

struct DataTet  {
  static const int NbOfVertices =4;
  static const int NbOfEdges =6;
  static const int NbOfFaces =4;
  static const int NT =1;
  static const int NbOfAdjElem =NbOfVertices;
  static const int NbOfVertexOnHyperFace =NbOfVertices-1;
  typedef Vertex3 V;
  typedef  V::Rd Rd ;
  static R mesure(  V *  pv[NbOfVertices]) 
  {    
    R3 AB(*pv[0],*pv[1]);
    R3 AC(*pv[0],*pv[2]);
    R3 AD(*pv[0],*pv[3]);
    return det(AB,AC,AD)/6.;
  }
  static const int (* const nvface)[3];// = nvfaceTet;
  static const int (* const nvedge)[2];//  = nvedgeTet;
  typedef R3 RdHat;
  typedef R2 RdHatBord;
  static RdHat PBord(const int * nvb,const RdHatBord& P)  { 
 //     cout << "PBORD : " << nvb[0] << " " <<  nvb[1] <<  nvb[2] << " " << P<< " -> " <<  RdHat::KHat[nvb[0]]*(1-P.x-P.y)+RdHat::KHat[nvb[1]]*(P.x)+RdHat::KHat[nvb[2]]*(P.y) 
//	<< "," <<  RdHat::KHat[nvb[0]] << "," <<  RdHat::KHat[nvb[1]] << "," << RdHat::KHat[nvb[2]] <<endl;
  return RdHat::KHat[nvb[0]]*(1-P.x-P.y)+RdHat::KHat[nvb[1]]*(P.x)+RdHat::KHat[nvb[2]]*(P.y) ;}  

};

class Tet: public GenericElement<DataTet>  {
public: 
  Tet() {}; // constructor empty for array

  
  R3 H(int i) const 
  { ASSERTION(i>=0 && i <4);
    R3 AB(at(this->nvface[i][0]),at(this->nvface[i][1]));
    R3 AC(at(this->nvface[i][0]),at(this->nvface[i][2]));
    return AB^AC/(6.*this->mesure());} // heigth 
 
    R3 n(int i) const 
    { ASSERTION(i>=0 && i <4);
	R3 AB(at(this->nvface[i][0]),at(this->nvface[i][1]));
	R3 AC(at(this->nvface[i][0]),at(this->nvface[i][2]));
	R3 N=AB^AC;
    return N/N.norme();} //  exterior normal  
    
  void Gradlambda(R3 * GradL) const
  {
    R3 V1(at(0),at(1));
    R3 V2(at(0),at(2));
    R3 V3(at(0),at(3));
    R det1=1./(6.*mesure());
    GradL[1]= (V2^V3)*det1;
    GradL[2]= (V3^V1)*det1;
    GradL[3]= (V1^V2)*det1;
    GradL[0]=-GradL[1]-GradL[2]-GradL[3];
  }

};

class Triangle3: public GenericElement<DataTriangle3>  {
public: 
  Triangle3() {}; // constructor empty for array

  Rd Edge(int i) const {ASSERTION(i>=0 && i <3);
    return Rd(this->at((i+1)%3),this->at((i+2)%3));}// opposite edge vertex i
  /*
  Rd H(int i) const { ASSERTION(i>=0 && i <3);
    Rd E=Edge(i);return E.perp()/(2.*this->mesure());} // heigth 
  
  void Gradlambda(Rd * GradL) const
  {
    GradL[1]= H(1);
    GradL[2]= H(2);
    GradL[0]=-GradL[1]-GradL[2];
  }
  */ 

};


class Mesh3 : public GenericMesh<Tet,Triangle3,Vertex3> { 
public:
  Mesh3(){}
  Mesh3(const string); 
  Mesh3(FILE *f);     
  Mesh3(int nnv, int nnt, int nnbe, Vertex3 *vv, Tet *tt, Triangle3 *bb); 
  Mesh3(int nnv, int nnbe, Vertex3 *vv, Triangle3 *bb);  // surface mesh 

  void GSave(FILE * f) const ;
  void GRead(FILE * f);
  int Save(const string & filename) const ;  
  int SaveSurface(const string & filename) const ;  
  int SaveSurface(const string & filename1, const string & filename2) const ;  
  void flipSurfaceMesh3(int surface_orientation);
  void read(istream &);
  //int Popen(const FILE *popenstream);

private:
  int load(const string & filename); 
  Mesh3(const Mesh3 &); // pas de construction par copie
  void operator=(const Mesh3 &);// pas affectation par copy 
};
// for the caracteristic method.
    int  WalkInTet(const Mesh3 & Th,int it, R3 & Phat,const R3 & U, R & dt);

} 
#endif
