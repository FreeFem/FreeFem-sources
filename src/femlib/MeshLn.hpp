// ORIG-DATE:     Dec 2007
// -*- Mode : c++ -*-
//
// SUMMARY  :  Model  mesh curve 3d   
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


#ifndef MESHLN_HPP_
#define MESHLN_HPP_


#include <cstdio>

// definition R
#include <cstdlib>
namespace Fem2D {
    
#include "R3.hpp"
#include "R2.hpp"
#include "R1.hpp"

}

using namespace ::std;
#include "GenericMesh.hpp"
//#include "MeshSn.hpp"

namespace Fem2D {
    
  typedef GenericVertex<R3> Vertex3;

  inline void Add(int *p,int n,int o) {
    for(int i=0;i<n;++i)
      p[i] += o;
  }
    
  static R1 PointHat[2] = { R1(0.), R1(1.) } ;
 
  struct DataPoint3  {
    static const int NbOfVertices =1;
    static const int NbOfEdges =0;
    static const int NbOfFaces =0;
    static const int NT =0;
    static const int NbOfAdjElem =1;
    static const int NbOfVertexOnHyperFace =1;
    typedef Vertex3 V;
    typedef  V::Rd Rd;
    static R mesure(  V * pv[NbOfVertices]  ) {
      return 0.;
    }
    typedef R0 RdHat;
    typedef R0 RdHatBord;   //hack no defined
    static RdHat PBord(const int * nvb,const RdHatBord & P)  { return R0() ;}
        
        
  };
    
    
    struct DataSeg3  {
        static const int NbOfVertices =2;
        static const int NbOfEdges =1;
        static const int NbOfFaces =0;
        static const int NT =0;
        static const int NbOfAdjElem =NbOfVertices;
        static const int NbOfVertexOnHyperFace =NbOfVertices-1;
        typedef Vertex3 V;
        typedef  V::Rd Rd;
        static R mesure(  V *  pv[NbOfVertices]) {
            return R3(*pv[0],*pv[1]).norme();
        }
        typedef R1 RdHat;
        typedef R0 RdHatBord;
        static RdHat PBord(const int * nvb,const RdHatBord &P)  { return RdHat(*nvb) ;}
        
        
    };
    
 
  class BoundaryPointL: public GenericElement<DataPoint3>
  {
  public: 
    BoundaryPointL() {}; // constructor empty for array
      
      
  };

  
  
  class EdgeL: public GenericElement<DataSeg3> {
  public:
    EdgeL() {}; // constructor empty for array
        
    R1 H(int i) const { ASSERTION(i>=0 && i <1);
      return (2-i)/mesure();} // heigth
        
    void Gradlambda(R3 * GradL) const
    {
        R3 V(at(0),at(1));
        GradL[1]= V/(V.norme2());
        GradL[0]= -GradL[1];
    }

    R3 NormalSUnitaire() const {
        R3 V(at(0),at(1));
        R3 N = V^R3(0,0,1);
        return N/N.norme();
    }

  };
    
    
    
  template<typename Mesh> void GSave2(FILE * ff,const Mesh & Th) ;
    
  class MeshL : public GenericMesh<EdgeL,BoundaryPointL,Vertex3> {
  public:
    // mapping for surface/line vertices
    int *mapSurf2Curv;
    int *mapCurv2Surf;
    MeshL():mapSurf2Curv(0),mapCurv2Surf(0) {};
    MeshL(const string);
    MeshL(const string filename, bool cleanmesh, bool removeduplicate=false, bool rebuildboundary=false, int orientation=1, double precis_mesh=1e-7, double ridgeangledetection=8.*atan(1.)/9.);
      
    void read(istream &f);
    void readmsh(ifstream & f,int offset);
    MeshL(FILE *f,int offset=0);
    MeshL(int nnv, int nnt, int nnbe, Vertex3 *vv, EdgeL *tt, BoundaryPointL *bb, bool cleanmesh=false, bool removeduplicate=false, bool rebuildboundary=false, int orientation=1, double precis_mesh=1e-7, double ridgeangledetection=8.*atan(1.)/9.);
    MeshL(const Serialize&);

    int load(const string & filename);
    const Element * Find( Rd P, R1 & Phat,bool & outside,const Element * tstart=0) const;
    int Save(const string & filename) const;
    //void flipSurfaceMeshS(int surface_orientation);
    void GSave(FILE * f,int offset=0) const ;
    void GRead(FILE * f,int offset);
    double hmin() const;
    //int Save(const string & filename) const;
    //Serialize serialize_withBorderMesh() const;
    void BuildBorderPt(const double angle=8.*atan(1.)/9.);

    ~MeshL() {
      delete [] mapSurf2Curv ;
      delete [] mapCurv2Surf ;
            
      SHOWVERB(cout << " %%%% delete MeshL"<< this << endl) ; }
  private:
    MeshL(const MeshL &); // pas de construction par copie
    void operator=(const MeshL &);// pas affectation par copy
  };
    
}

#endif

