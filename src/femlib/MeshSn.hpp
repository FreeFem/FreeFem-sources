// ORIG-DATE:     Dec 2007
// -*- Mode : c++ -*-
//
// SUMMARY  :  Model  mesh surface 3d   
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


#ifndef MESHSN_HPP_
#define MESHSN_HPP_


#include <cstdio>

// definition R
#include <cstdlib>
namespace Fem2D {
    
#include "R3.hpp"
// #include "R2.hpp"

}

using namespace ::std;
#include "GenericMesh.hpp"
#include "MeshLn.hpp"

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
 

    class BoundaryEdgeS: public GenericElement<DataSeg3>  
  {
  public: 
      BoundaryEdgeS() {}; // constructor empty for array
      
      
  };

  
  
    class TriangleS: public GenericElement<DataTriangle3> { //public Triangle3  {
    public:
        TriangleS() {}; // constructor empty for array
        
        void Gradlambda(Rd * GradL) const
        {
            R3 Normal = Edge(2)^Edge(1);
            R N = Normal.norme2();
            for(int i=0 ; i<3 ; i++)
                GradL[i]= (Normal^Edge(i)) / N;
        }
        
        R3 NFrenet() const {
            ASSERTION(i>=0 && i <3);
            return R3( Edge(2)^Edge(1) );
        }
        R3 NFrenetUnitaire() const {// Add FH. Dec. 2019 
            R3 N(Edge(2)^Edge(1) );
            return N/N.norme();
        }
    };
    
    
    
    template<typename Mesh> void GSave2(FILE * ff,const Mesh & Th) ;
    
    class MeshS : public GenericMesh<TriangleS,BoundaryEdgeS,Vertex3> {
    public:
         // mapping for volume/surface vertices
        int *mapSurf2Vol;
        int *mapVol2Surf;
        MeshL *meshL;
        int nPoints;
        MeshS():mapSurf2Vol(0),mapVol2Surf(0),meshL(0) {};
        MeshS(const string, double ridgeangledetection=8.*atan(1.)/9.);
        //MeshS(const string, const long);
        MeshS(const string filename, bool cleanmesh, bool removeduplicate=false, bool rebuildboundary=false, int orientation=1, double precis_mesh=1e-7, double ridgeangledetection=8.*atan(1.)/9.);
        void read(istream &f);
        void readmsh(ifstream & f,int offset);
        MeshS(FILE *f,int offset=0);
        //MeshS(const string);
        MeshS(int nnv, int nnt, int nnbe, Vertex3 *vv, TriangleS *tt, BoundaryEdgeS *bb, bool cleanmesh=false, bool removeduplicate=false, bool rebuildboundary=false, int orientation=1, double precis_mesh=1e-6);
        MeshS(const Serialize&);
        MeshS(const Serialize &serialized, int withSurface);
        int load(const string & filename);
        const Element * Find( Rd P, R2 & Phat,bool & outside,const Element * tstart=0) const;
        int Save(const string & filename);
        void flipSurfaceMeshS(int surface_orientation);
        void GSave(FILE * f,int offset=0) const ;
        void GRead(FILE * f,int offset);
        double hmin() const;
        int Save(const string & filename) const;
        void BuildBdElem(double angle=8.*atan(1.)/9.);  // default angle = 40 deg);
        Serialize serialize_withBorderMesh() const;
        void BuildMeshL(double angle=8.*atan(1.)/9.);  // default angle = 40 deg);
        void BuildCurvBasis();
        ~MeshS() {
            if (verbosity>4) cout << "destroy meshS" << this << " " << this->mapSurf2Vol << " " << this->mapVol2Surf << " destroy meshL " << this->meshL << endl;
            delete [] mapSurf2Vol ;
            delete [] mapVol2Surf ;
            if (meshL)
                meshL->destroy();//  Add clean mesh if necessary ...FH and AF. april 2019
            
            
            SHOWVERB(cout << " %%%% delete MeshS"<< this << endl) ; }
    private:
        MeshS(const MeshS &); // pas de construction par copie
        void operator=(const MeshS &);// pas affectation par copy
    };
    
}

#endif

