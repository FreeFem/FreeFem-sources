#ifndef FREEFEMBAMG_H_
#define FREEFEMBAMG_H_

namespace bamg {
extern  void (*MeshIstreamErrorHandler)(ios & );
class Triangles;
}


  Fem2D::Mesh *  ReadMeshbamg( string * const & s) ;
  Fem2D::Mesh *  ReadTriangulate( string * const & s) ;
  Fem2D::Mesh *bamg2msh( bamg::Triangles* tTh,bool renumbering=false); 
  bamg::Triangles * msh2bamg(const Fem2D::Mesh & Th, double cutoffradian=-1.0 ) ;
  bamg::Triangles * msh2bamg(const Fem2D::Mesh & Th,double cutoffradian,  int  nbdfv, int * ndfv,int  nbdfe, int * ndfe);

  Fem2D::Mesh *  BuildMesh(Stack stack  , E_BorderN const * const & b) ;
  Fem2D::Mesh *  BuildMeshBorder(Stack stack  , E_BorderN const * const & b) ;
  Fem2D::Mesh * MoveTheMesh(const Fem2D::Mesh &Th,const KN_<double> & u,const KN_<double> &v);

#endif
