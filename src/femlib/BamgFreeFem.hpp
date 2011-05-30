// ORIG-DATE:     Dec 2007
// -*- Mode : c++ -*-
//
// SUMMARY  :  Bamglib interface header  
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

#ifndef FREEFEMBAMG_H_
#define FREEFEMBAMG_H_

namespace bamg {
extern  void (*MeshIstreamErrorHandler)(ios & );
class Triangles;
}


  Fem2D::Mesh *  ReadMeshbamg( string * const & s) ;
  Fem2D::Mesh *  ReadTriangulate( string * const & s) ;
  Fem2D::Mesh *  Triangulate( const  KN_<double> & xx,const  KN_<double> & yy) ;
  Fem2D::Mesh *bamg2msh( bamg::Triangles* tTh,bool renumbering=false); 
  bamg::Triangles * msh2bamg(const Fem2D::Mesh & Th, double cutoffradian=-1.0,long * reqedgeslab=0,int nreqedgeslab =0) ;
  bamg::Triangles * msh2bamg(const Fem2D::Mesh & Th,double cutoffradian,  int  nbdfv, int * ndfv,int  nbdfe, int * ndfe,
     long * reqedgeslab=0,int nreqedgeslab=0);

  Fem2D::Mesh *  BuildMesh(Stack stack, E_BorderN const * const & b,bool justboundary,int nbvmax=0,bool Requiredboundary=true);
  Fem2D::Mesh *  BuildMesh(Stack stack  , E_BorderN const * const & b,bool Requiredboundary) ;
  Fem2D::Mesh *  BuildMeshBorder(Stack stack  , E_BorderN const * const & b) ;
  Fem2D::Mesh * MoveTheMesh(const Fem2D::Mesh &Th,const KN_<double> & u,const KN_<double> &v);
  Fem2D::Mesh *  buildmeshbamg( string * const & s,int =0);
#endif
