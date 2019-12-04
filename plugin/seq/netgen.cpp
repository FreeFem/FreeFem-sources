/****************************************************************************/
/* This file is part of FreeFEM.                                            */
/*                                                                          */
/* FreeFEM is free software: you can redistribute it and/or modify          */
/* it under the terms of the GNU Lesser General Public License as           */
/* published by the Free Software Foundation, either version 3 of           */
/* the License, or (at your option) any later version.                      */
/*                                                                          */
/* FreeFEM is distributed in the hope that it will be useful,               */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of           */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            */
/* GNU Lesser General Public License for more details.                      */
/*                                                                          */
/* You should have received a copy of the GNU Lesser General Public License */
/* along with FreeFEM. If not, see <http://www.gnu.org/licenses/>.          */
/****************************************************************************/
// SUMMARY : ...
// LICENSE : LGPLv3
// ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
// AUTHORS : Jacques Morice
// E-MAIL  : jacques.morice@ann.jussieu.fr

/* clang-format off */
//ff-c++-LIBRARY-dep:
//ff-c++-cpp-dep:
/* clang-format on */

/*
 * Thank to the ARN ()  FF2A3 grant
 * ref:ANR-07-CIS7-002-01
 */

// FH   July 2009
// comment all
// Th3_t->BuildBound();
// Th3_t->BuildAdj();
// Th3_t->Buildbnormalv();
// Th3_t->BuildjElementConteningVertex();
// is now in the constructor of Mesh3 to be consistante.
//

#include <iostream>
#include <cfloat>
#include <cmath>
#include <complex>
using namespace std;
#include "error.hpp"
#include "AFunction.hpp"
using namespace std;
#include "rgraph.hpp"
#include "RNM.hpp"
#include "fem.hpp"

#include "FESpacen.hpp"
#include "FESpace.hpp"

#include "MatriceCreuse_tpl.hpp"
#include "MeshPoint.hpp"
#include "Operator.hpp"
#include "lex.hpp"

#include "lgfem.hpp"
#include "lgmesh3.hpp"
#include "lgsolver.hpp"
#include "problem.hpp"
#include "msh3.hpp"

#include <set>
#include <vector>
#include <fstream>

namespace nglib {
#include "./includenetgen/nglib.h"
}
using namespace nglib;

using namespace Fem2D;

class RemplissageNetgen_Op : public E_F0mps {
 public:
  Expression eTh;
  Expression xx, yy, zz;
  static const int n_name_param = 5;    //
  static basicAC_F0::name_and_type name_param[];
  Expression nargs[n_name_param];
  KN_< long > arg(int i, Stack stack, KN_< long > a) const {
    return nargs[i] ? GetAny< KN_< long > >((*nargs[i])(stack)) : a;
  }

  KN_< double > arg(int i, Stack stack, KN_< double > a) const {
    return nargs[i] ? GetAny< KN_< double > >((*nargs[i])(stack)) : a;
  }

  double arg(int i, Stack stack, double a) const {
    return nargs[i] ? GetAny< double >((*nargs[i])(stack)) : a;
  }

  int arg(int i, Stack stack, int a) const {
    return nargs[i] ? GetAny< long >((*nargs[i])(stack)) : a;
  }

  string *arg(int i, Stack stack, string *a) const {
    return nargs[i] ? GetAny< string * >((*nargs[i])(stack)) : a;
  }

 public:
  RemplissageNetgen_Op(const basicAC_F0 &args, Expression tth) : eTh(tth) {
    if (verbosity) {
      cout << "construction par RemplissageNetgen_Op" << endl;
    }

    args.SetNameParam(n_name_param, name_param, nargs);
  }

  AnyType operator( )(Stack stack) const;
};

basicAC_F0::name_and_type RemplissageNetgen_Op::name_param[] = {
  {"reftet", &typeid(long)},
  {"refface", &typeid(KN_< long >)},
  // Netgen Options
  {"maxh", &typeid(double)},
  {"secondorder", &typeid(long)},
  {"meshsizefilename", &typeid(string *)}
  /* // Parametres de netgen non d�finis encore
   * // other parameters of netgen : see libsrc/meshing/meshtype.hpp of netgen directory for more
   * information
   * // parameters defined by default is given in MeshingParameters :: MeshingParameters () in
   * libsrc/meshing/meshtype.cpp.
   *
   * {  "optimize3d", &typeid(string*)},
   * {  "optsteps3d", &typeid(long)},
   * {  "optimize2d", &typeid(string*)},
   * {  "optsteps2d", &typeid(long)},
   * {  "opterrpow", &typeid(double)},
   * {  "blockfill", &typeid(long)},
   * {  "filldist", &typeid(double)},
   * {  "safety", &typeid(double)},
   * {  "relinnersafety", &typeid(double)},
   * {  "uselocalh", &typeid(long)},
   * {  "grading", &typeid(double)},
   * {  "delaunay", &typeid(long)},
   *
   * {  "minh", &typeid(double)},
   * {  "startinsurface", &typeid(long)},
   * {  "checkoverlap", &typeid(long)}, // debug
   * {  "checkoverlappingboundary", &typeid(long)},
   * {  "checkchartboundary", &typeid(long)},
   * {  "curvaturesafety", &typeid(double)},
   * {  "segmentsperedge", &typeid(double)},
   * {  "parthread",&typeid(long)}, //use parallel threads
   * {  "elsizeweight", &typeid(double)}, // weight of element size w.r.t element shape
   * // from mp3:
   * {  "giveuptol", &typeid(long)},
   * {  "maxoutersteps", &typeid(long)},
   * {  "starshapeclass", &typeid(long)},
   * {  "baseelnp", &typeid(long)},
   * {  "sloppy", &typeid(long)},
   * {  "badellimit", &typeid(double)}, /// limit for max element angle (150-180)
   * {  "check_impossible", &typeid(bool)},
   * {  "elementorder", &typeid(long)},
   * {  "quad", &typeid(long)},
   * {  "inverttets", &typeid(long)},
   * {  "inverttrigs", &typeid(long)},
   * {  "autozrefine", &typeid(long)}
   */
};

class RemplissageNetgen : public OneOperator {
 public:
  RemplissageNetgen( ) : OneOperator(atype< pmesh3 >( ), atype< pmesh3 >( )) {}

  E_F0 *code(const basicAC_F0 &args) const {
    return new RemplissageNetgen_Op(args, t[0]->CastTo(args[0]));
  }
};

AnyType RemplissageNetgen_Op::operator( )(Stack stack) const {
  MeshPoint *mp(MeshPointStack(stack)), mps = *mp;
  Mesh3 *pTh3 = GetAny< Mesh3 * >((*eTh)(stack));

  ffassert(pTh3);
  Mesh3 &Th3 = *pTh3;
  Mesh3 *m = pTh3;    // question a quoi sert *m ??

  // lecture des arguments
  KN< long > zzempty;
  long nrtet(arg(0, stack, 1));
  KN< long > nrf(arg(1, stack, zzempty));

  ffassert(nrf.N( ) % 2 == 0);

  int nbv = Th3.nv;     // nombre de sommet
  int nbt = Th3.nt;     // nombre de tetrahedre
  int nbe = Th3.nbe;    // nombre de surfaces

  // check consistency of surface mesh
  if (verbosity > 1) {
    cout << "check :: orientation des surfaces" << endl;
  }

  Th3.BuildSurfaceAdj( );
  if (verbosity > 1) {
    cout << "fin check :: orientation des surfaces" << endl;
  }

  if (nbt != 0) {
    cerr << " The mesh must be a 3D surface mesh " << endl;
    exit(1);
  }

  if (verbosity > 1) {
    cout << " ======================= " << endl;
  }

  if (verbosity > 1) {
    cout << " == RemplissageNetgen == " << endl;
  }

  Ng_Mesh *netgen_mesh;
  Ng_Init( );
  // creates mesh structure
  netgen_mesh = Ng_NewMesh( );

  double point[3];
  int trig[3], tet[4];

  for (int ii = 0; ii < Th3.nv; ii++) {
    point[0] = Th3.vertices[ii].x;
    point[1] = Th3.vertices[ii].y;
    point[2] = Th3.vertices[ii].z;

    Ng_AddPoint(netgen_mesh, point);
  }

  for (int ii = 0; ii < Th3.nbe; ii++) {
    const Triangle3 &K(Th3.be(ii));
    int label = K.lab;

    for (int jj = 0; jj < 3; jj++) {
      trig[jj] = Th3.operator( )(K[jj]) + 1;
    }

    Ng_AddSurfaceElement(netgen_mesh, NG_TRIG, trig);    // , &label);
  }

  Ng_Meshing_Parameters netgen_mp;

  if (nargs[2]) {
    double netgen_maxh = GetAny< double >((*nargs[2])(stack));
    netgen_mp.maxh = netgen_maxh;
  }

  if (nargs[3]) {
    int netgen_secondorder = GetAny< long >((*nargs[3])(stack));
    netgen_mp.secondorder = netgen_secondorder;
  }

  if (nargs[4]) {
    string *netgen_meshsize_filename = GetAny< string * >((*(nargs[4]))(stack));
    size_t size_filename = netgen_meshsize_filename->size( ) + 1;
    char *netgen_meshsize_filename_char = new char[size_filename];
    strncpy(netgen_meshsize_filename_char, netgen_meshsize_filename->c_str( ), size_filename);

    netgen_mp.meshsize_filename = netgen_meshsize_filename_char;
  }

  // Essai des restrictions pour tetgen

  double pmin[3], pmax[3];
  pmin[0] = 0.9;
  pmin[1] = -0.1;
  pmin[2] = -0.1;

  pmax[0] = 1.5;
  pmax[1] = 3.;
  pmax[2] = 1.5;

  Ng_RestrictMeshSizeBox(netgen_mesh, pmin, pmax, 0.05);

  cout << "start remeshing" << endl;
  Ng_GenerateVolumeMesh(netgen_mesh, &netgen_mp);
  cout << "meshing done" << endl;

  /* Transformation netgen -> freefem++ */

  map< int, int > mapface;

  for (int i = 0; i < nrf.N( ); i += 2) {
    if (nrf[i] != nrf[i + 1]) {
      mapface[nrf[i]] = nrf[i + 1];
    }
  }

  // read information of netgen
  int netgen_nv = Ng_GetNP(netgen_mesh);
  int netgen_nt = Ng_GetNE(netgen_mesh);
  int netgen_nbe = Ng_GetNSE(netgen_mesh);
  Vertex3 *v = new Vertex3[netgen_nv];
  Tet *t = new Tet[netgen_nt];
  Triangle3 *b = new Triangle3[netgen_nbe];
  // generation des nouveaux sommets
  Vertex3 *vv = v;
  Tet *tt = t;
  Triangle3 *bb = b;

  cout << " donnee sortie netgen:  Vertex" << netgen_nv << " Tetrahedre " << netgen_nt << "  "
       << netgen_nbe << endl;

  for (int ii = 0; ii < netgen_nv; ii++) {
    Ng_GetPoint(netgen_mesh, ii + 1, point);
    vv->x = point[0];
    vv->y = point[1];
    vv->z = point[2];
    vv->lab = 1;
    vv++;
  }

  for (int ii = 0; ii < netgen_nt; ii++) {
    int iv[4];
    Ng_GetVolumeElement(netgen_mesh, ii + 1, iv);

    for (int jj = 0; jj < 4; jj++) {
      iv[jj] = iv[jj] - 1;
    }

    (*tt++).set(v, iv, nrtet);
  }

  for (int ii = 0; ii < netgen_nbe; ii++) {
    const Triangle3 &K(Th3.be(ii));
    int label = K.lab;
    // int label=0;
    int iv[3];
    Ng_GetSurfaceElement(netgen_mesh, ii + 1, iv);

    for (int jj = 0; jj < 3; jj++) {
      iv[jj] = iv[jj] - 1;
    }

    (*bb++).set(v, iv, label);
  }

  Ng_DeleteMesh(netgen_mesh);
  Ng_Exit( );

  Mesh3 *Th3_t = new Mesh3(netgen_nv, netgen_nt, netgen_nbe, v, t, b);

  Th3_t->BuildGTree( );
  Add2StackOfPtr2FreeRC(stack, Th3_t);

  *mp = mps;

  return Th3_t;
}

class Netgen_STL_Op : public E_F0mps {
 public:
  Expression filename;
  static const int n_name_param = 3;    //
  static basicAC_F0::name_and_type name_param[];
  Expression nargs[n_name_param];
  KN_< long > arg(int i, Stack stack, KN_< long > a) const {
    return nargs[i] ? GetAny< KN_< long > >((*nargs[i])(stack)) : a;
  }

  KN_< double > arg(int i, Stack stack, KN_< double > a) const {
    return nargs[i] ? GetAny< KN_< double > >((*nargs[i])(stack)) : a;
  }

  double arg(int i, Stack stack, double a) const {
    return nargs[i] ? GetAny< double >((*nargs[i])(stack)) : a;
  }

  int arg(int i, Stack stack, int a) const {
    return nargs[i] ? GetAny< long >((*nargs[i])(stack)) : a;
  }

  string *arg(int i, Stack stack, string *a) const {
    return nargs[i] ? GetAny< string * >((*nargs[i])(stack)) : a;
  }

 public:
  Netgen_STL_Op(const basicAC_F0 &args, Expression ffname) : filename(ffname) {
    if (verbosity) {
      cout << "construction par RemplissageNetgen_Op" << endl;
    }

    args.SetNameParam(n_name_param, name_param, nargs);
  }

  AnyType operator( )(Stack stack) const;
};

basicAC_F0::name_and_type Netgen_STL_Op::name_param[] = {
  // Netgen Options
  {"maxh", &typeid(double)},
  {"secondorder", &typeid(long)},
  {"meshsizefilename", &typeid(string *)}};

class Netgen_STL : public OneOperator {
 public:
  Netgen_STL( ) : OneOperator(atype< pmesh3 >( ), atype< string * >( )) {}

  E_F0 *code(const basicAC_F0 &args) const {
    return new Netgen_STL_Op(args, t[0]->CastTo(args[0]));
  }
};

AnyType Netgen_STL_Op::operator( )(Stack stack) const {
  MeshPoint *mp(MeshPointStack(stack)), mps = *mp;
  string *pffname = GetAny< string * >((*filename)(stack));
  size_t size_pffname = pffname->size( ) + 1;
  char *char_pffname = new char[size_pffname];

  strncpy(char_pffname, pffname->c_str( ), size_pffname);

  // lecture des arguments

  if (verbosity > 1) {
    cout << " ===================== " << endl;
  }

  if (verbosity > 1) {
    cout << " ==   Netgen_STL    == " << endl;
  }

  int i, j, rv;
  Ng_Mesh *netgen_mesh;
  Ng_STL_Geometry *netgen_geom;
  Ng_Meshing_Parameters netgen_mp;

  if (nargs[0]) {
    double netgen_maxh = GetAny< double >((*nargs[0])(stack));
    netgen_mp.maxh = netgen_maxh;
  }

  if (nargs[1]) {
    int netgen_secondorder = GetAny< long >((*nargs[1])(stack));
    netgen_mp.secondorder = netgen_secondorder;
  }

  if (nargs[2]) {
    string *netgen_meshsize_filename = GetAny< string * >((*(nargs[2]))(stack));
    size_t size_filename = netgen_meshsize_filename->size( ) + 1;
    char *netgen_meshsize_filename_char = new char[size_filename];
    strncpy(netgen_meshsize_filename_char, netgen_meshsize_filename->c_str( ), size_filename);

    netgen_mp.meshsize_filename = netgen_meshsize_filename_char;
  }

  Ng_Init( );

  netgen_geom = Ng_STL_LoadGeometry(char_pffname);
  if (!netgen_geom) {
    cerr << "Ng_STL_LoadGeometry return NULL" << endl;
    exit(1);
  }

  rv = Ng_STL_InitSTLGeometry(netgen_geom);
  cout << "InitSTLGeometry: NG_result=" << rv << endl;

  netgen_mesh = Ng_NewMesh( );

  rv = Ng_STL_MakeEdges(netgen_geom, netgen_mesh, &netgen_mp);
  cout << "Make Edges: Ng_result=" << rv << endl;

  rv = Ng_STL_GenerateSurfaceMesh(netgen_geom, netgen_mesh, &netgen_mp);
  cout << "Generate Surface Mesh: Ng_result=" << rv << endl;

  rv = Ng_GenerateVolumeMesh(netgen_mesh, &netgen_mp);
  cout << "Generate Volume Mesh: Ng_result=" << rv << endl;

  // read information of netgen
  int netgen_nv = Ng_GetNP(netgen_mesh);
  int netgen_nt = Ng_GetNE(netgen_mesh);
  int netgen_nbe = Ng_GetNSE(netgen_mesh);
  Vertex3 *v = new Vertex3[netgen_nv];
  Tet *t = new Tet[netgen_nt];
  Triangle3 *b = new Triangle3[netgen_nbe];
  // generation des nouveaux sommets
  Vertex3 *vv = v;
  Tet *tt = t;
  Triangle3 *bb = b;

  cout << " donnee sortie netgen:  Vertex" << netgen_nv << " Tetrahedre " << netgen_nt << "  "
       << netgen_nbe << endl;

  for (int ii = 0; ii < netgen_nv; ii++) {
    double point[3];
    Ng_GetPoint(netgen_mesh, ii + 1, point);
    vv->x = point[0];
    vv->y = point[1];
    vv->z = point[2];
    vv->lab = 1;
    vv++;
  }

  for (int ii = 0; ii < netgen_nt; ii++) {
    int nrtet = 1;
    int iv[4];
    Ng_GetVolumeElement(netgen_mesh, ii + 1, iv, &nrtet);

    for (int jj = 0; jj < 4; jj++) {
      iv[jj] = iv[jj] - 1;
    }

    (*tt++).set(v, iv, nrtet);
  }

  for (int ii = 0; ii < netgen_nbe; ii++) {
    int label;
    int iv[3];
    Ng_GetSurfaceElement(netgen_mesh, ii + 1, iv, &label);

    for (int jj = 0; jj < 3; jj++) {
      iv[jj] = iv[jj] - 1;
    }

    (*bb++).set(v, iv, label);
  }

  Ng_DeleteMesh(netgen_mesh);
  Ng_Exit( );
  Mesh3 *Th3_t = new Mesh3(netgen_nv, netgen_nt, netgen_nbe, v, t, b);

  Th3_t->BuildGTree( );
  Add2StackOfPtr2FreeRC(stack, Th3_t);

  *mp = mps;

  return Th3_t;
}

class Netgen_Face {
 public:
  int label;
  int domin;
  int domout;
  int bcd;
  Netgen_Face *pp;
  Netgen_Face( );
  Netgen_Face(int a, int b, int c, int d) : label(a), domin(b), domout(c), bcd(d), pp(0){};

  int Test_Face(int a, int b, int c, int d) {
    if (label == a && domin == b && domout == c && bcd == d) {
      return 1;
    } else {
      return 0;
    }
  }

  Netgen_Face *Next_Face( ) { return pp; }

  void Initialisation_Face(int a, int b, int c, int d) {
    label = a;
    domin = b;
    domout = c;
    bcd = d;
  }

  void Initialisation_Next_Face(Netgen_Face *e) { pp = e; }

 private:
  ~Netgen_Face( );
};

class Netgen_LoadMesh_Op : public E_F0mps {
 public:
  Expression filename;
  static const int n_name_param = 2;    //
  static basicAC_F0::name_and_type name_param[];
  Expression nargs[n_name_param];

 public:
  Netgen_LoadMesh_Op(const basicAC_F0 &args, Expression ffname) : filename(ffname) {
    if (verbosity) {
      cout << "Load mesh given by Netgen " << endl;
    }

    args.SetNameParam(n_name_param, name_param, nargs);
  }

  AnyType operator( )(Stack stack) const;
};

basicAC_F0::name_and_type Netgen_LoadMesh_Op::name_param[] = {{"reftet", &typeid(long)},
                                                              {"renum", &typeid(long)}};

class Netgen_LoadMesh : public OneOperator {
 public:
  Netgen_LoadMesh( ) : OneOperator(atype< pmesh3 >( ), atype< string * >( )) {}

  E_F0 *code(const basicAC_F0 &args) const {
    return new Netgen_LoadMesh_Op(args, t[0]->CastTo(args[0]));
  }
};

Mesh3 *NETGEN_Load(const string &filename, const int &flagsnewlabelsurf) {
  char str[100];
  int i, n;
  double scale = 1;    // globflags.GetNumFlag ("scale", 1);

  // freefempp
  int nv, nt, nbe;
  int dimension;
  Vertex3 *v;

  // reading first the filex
  ifstream infile(filename.c_str( ));

  if (!infile.good( )) {
    cerr << "probleme in reading file" << filename << endl;
  }

  bool endmesh = false;

  while (infile.good( ) && !endmesh) {
    infile >> str;

    if (strcmp(str, "dimension") == 0) {
      infile >> dimension;
      if (dimension != 3) {
        cerr << "dimension will be 3 for reading a netgen file" << endl;
      }
    }

    if (strcmp(str, "surfaceelements") == 0) {
      infile >> n;
      cout << "surface element" << n << endl;
      nbe = n;
    }

    if (strcmp(str, "surfaceelementsgi") == 0) {
      infile >> n;
      cout << "surface element" << n << endl;
      nbe = n;
    }

    if (strcmp(str, "volumeelements") == 0) {
      infile >> n;
      cout << "volume element" << n << endl;
      nt = n;
    }

    if (strcmp(str, "edgesegments") == 0) {
      cout << "this parameter is not taken in consideration in freefem++" << endl;
    }

    if (strcmp(str, "edgesegmentsgi") == 0) {
      cout << "this parameter is not taken in consideration in freefem++" << endl;
    }

    if (strcmp(str, "edgesegmentsgi2") == 0) {
      cout << "this parameter is not taken in consideration in freefem++" << endl;
    }

    if (strcmp(str, "points") == 0) {
      infile >> n;
      cout << " vertex " << n << endl;
      nv = n;
      // read vertex coordinate
      v = new Vertex3[nv];

      for (i = 0; i < n; i++) {
        infile >> v[i].x >> v[i].y >> v[i].z;
        v[i].lab = 1;
      }
    }

    if (strcmp(str, "identifications") == 0) {
      cout << "this parameter is not taken in consideration in freefem++" << endl;
    }

    if (strcmp(str, "identificationtypes") == 0) {
      cout << "this parameter is not taken in consideration in freefem++" << endl;
    }

    if (strcmp(str, "materials") == 0) {
      cout << "this parameter is not taken in consideration in freefem++" << endl;
    }

    if (strcmp(str, "bcnames") == 0) {
      cout << "this parameter is not taken in consideration in freefem++" << endl;
    }

    if (strcmp(str, "singular_points") == 0) {
      cout << "this parameter is not taken in consideration in freefem++" << endl;
    }

    if (strcmp(str, "singular_edge_left") == 0) {
      cout << "this parameter is not taken in consideration in freefem++" << endl;
    }

    if (strcmp(str, "singular_edge_right") == 0) {
      cout << "this parameter is not taken in consideration in freefem++" << endl;
    }

    if (strcmp(str, "singular_face_inside") == 0) {
      cout << "this parameter is not taken in consideration in freefem++" << endl;
    }

    if (strcmp(str, "singular_face_outside") == 0) {
      cout << "this parameter is not taken in consideration in freefem++" << endl;
    }

    if (strcmp(str, "endmesh") == 0) {
      endmesh = true;
    }

    strcpy(str, "");
  }

  infile.close( );

  // Allocate Array
  Tet *t;
  if (nt != 0) {
    t = new Tet[nt];
  }

  Tet *tt = t;
  Triangle3 *b = new Triangle3[nbe];
  Triangle3 *bb = b;

  // second lecture
  ifstream infile2(filename.c_str( ));
  if (!infile2.good( )) {
    cerr << "probleme in reading file" << filename << endl;
  }

  bool endmesh2 = false;

  while (infile2.good( ) && !endmesh2) {
    infile2 >> str;

    if (strcmp(str, "dimension") == 0) {
      infile2 >> dimension;
    }

    if (strcmp(str, "surfaceelements") == 0) {
      infile2 >> n;
      cout << "surface elements " << n << endl;

      int flagsnewlabelsurf = 1;
      if (flagsnewlabelsurf == 0) {
        for (i = 1; i <= n; i++) {
          int j;
          int surfnr, bcp, domin, domout, nep, faceind = 0;

          infile2 >> surfnr >> bcp >> domin >> domout;
          /*
           * // details des entrees
           * surfnr: label de la surface
           * bcp   : condition de bord
           * domin : label du domaine interieure
           * domout: label du domaine exterieure
           */

          infile2 >> nep;
          if (nep != 3) {
            cerr << "only triangle elements are considered in Freefem++" << endl;
          }

          nep = 3;

          int label = faceind;
          int iv[3];

          for (j = 0; j < nep; j++) {
            infile2 >> iv[j];    // tri.PNum(j);
          }

          for (j = 0; j < nep; j++) {
            iv[j]--;
          }

          (bb++)->set(v, iv, label);    // add boundary
        }
      }

      if (flagsnewlabelsurf == 1) {
        int maxlabelsurf = 10000;
        int *ngf = new int[maxlabelsurf * 4];
        int nbfacesnetg = 0;

        for (i = 1; i <= n; i++) {
          int j;
          int surfnr, bcp, domin, domout, nep, faceind = 0;

          infile2 >> surfnr >> bcp >> domin >> domout;
          /*
           * // details des entrees
           * surfnr: label de la surface
           * bcp   : condition de bord
           * domin : label du domaine interieure
           * domout: label du domaine exterieure
           */

          int label;

          label = -1;

          for (j = 0; j < nbfacesnetg; j++) {
            if (ngf[4 * j] == surfnr && ngf[4 * j + 1] == bcp && ngf[4 * j + 2] == domin &&
                ngf[4 * j + 3] == domout) {
              label = j + 1;
            }
          }

          if (label == -1) {
            ngf[4 * nbfacesnetg] = surfnr;
            ngf[4 * nbfacesnetg + 1] = bcp;
            ngf[4 * nbfacesnetg + 2] = domin;
            ngf[4 * nbfacesnetg + 3] = domout;
            label = nbfacesnetg + 1;
            nbfacesnetg++;
            if (nbfacesnetg == 10000) {
              cerr << "number of different surfaces is larger than 10000" << endl;
              exit(1);
            }
          }

          infile2 >> nep;
          if (nep != 3) {
            cerr << "only triangle elements are considered in Freefem++" << endl;
          }

          nep = 3;

          int iv[3];

          for (j = 0; j < nep; j++) {
            infile2 >> iv[j];    // tri.PNum(j);
          }

          for (j = 0; j < nep; j++) {
            iv[j]--;
          }

          (bb++)->set(v, iv, label);    // add boundary
        }

        delete[] ngf;
      }

      if (flagsnewlabelsurf == 2) {
        Netgen_Face *debutpngdf;
        Netgen_Face *iteratorpngdf;
        int nbfacesnetg = 0;

        for (i = 1; i <= n; i++) {
          int j;
          int surfnr, bcp, domin, domout, nep, faceind = 0;

          infile2 >> surfnr >> bcp >> domin >> domout;
          /*
           * // details des entrees
           * surfnr: label de la surface
           * bcp   : condition de bord
           * domin : label du domaine interieure
           * domout: label du domaine exterieure
           */

          int label;

          label = -1;
          iteratorpngdf = debutpngdf;

          for (j = 0; j < nbfacesnetg; j++) {
            if (iteratorpngdf->Test_Face(surfnr, domin, domout, bcp) == 1) {
              label = j + 1;
            }

            iteratorpngdf = iteratorpngdf->Next_Face( );
          }

          if (label == -1) {
            if (nbfacesnetg == 0) {
              debutpngdf = new Netgen_Face(surfnr, domin, domout, bcp);
            } else {
              // iteratorpngdf = debutpngdf;
              Netgen_Face *bb = debutpngdf;

              for (j = 0; j < nbfacesnetg; j++) {
                iteratorpngdf = bb;
                bb = iteratorpngdf->Next_Face( );
              }

              bb = new Netgen_Face(surfnr, domin, domout, bcp);
              iteratorpngdf->Initialisation_Next_Face(bb);
            }

            label = nbfacesnetg + 1;
            nbfacesnetg++;
          }

          infile2 >> nep;
          if (nep != 3) {
            cerr << "only triangle elements are considered in Freefem++" << endl;
          }

          nep = 3;

          int iv[3];

          for (j = 0; j < nep; j++) {
            infile2 >> iv[j];    // tri.PNum(j);
          }

          for (j = 0; j < nep; j++) {
            iv[j]--;
          }

          (bb++)->set(v, iv, label);    // add boundary
        }
      }
    }

    if (strcmp(str, "surfaceelementsgi") == 0) {
      infile2 >> n;
      cout << "surface element sgi" << n << endl;

      if (flagsnewlabelsurf == 0) {
        for (i = 1; i <= n; i++) {
          int j;
          int surfnr, bcp, domin, domout, nep, faceind = 0;
          infile2 >> surfnr >> bcp >> domin >> domout;

          infile2 >> nep;
          if (nep != 3) {
            cerr << "only triangle elements are considered in Freefem++" << endl;
          }

          nep = 3;

          int label = surfnr;
          int iv[3], infoP[3];

          for (j = 0; j < nep; j++) {
            infile2 >> iv[j];    // tri.PNum(j);
          }

          for (j = 0; j < nep; j++) {
            infile2 >> infoP[j];    // tri.PNum(j); // No Need By Freefem++
          }

          for (j = 0; j < nep; j++) {
            iv[j]--;
          }

          (bb++)->set(v, iv, label);    // add boundary
        }
      }

      if (flagsnewlabelsurf == 1) {
        int maxlabelsurf = 10000;
        int *ngf = new int[maxlabelsurf * 4];
        int nbfacesnetg = 0;

        for (i = 1; i <= n; i++) {
          int j;
          int surfnr, bcp, domin, domout, nep, faceind = 0;
          infile2 >> surfnr >> bcp >> domin >> domout;

          int label;

          label = -1;

          for (j = 0; j < nbfacesnetg; j++) {
            if (ngf[4 * j] == surfnr && ngf[4 * j + 1] == bcp && ngf[4 * j + 2] == domin &&
                ngf[4 * j + 3] == domout) {
              label = j + 1;
            }
          }

          if (label == -1) {
            ngf[4 * nbfacesnetg] = surfnr;
            ngf[4 * nbfacesnetg + 1] = bcp;
            ngf[4 * nbfacesnetg + 2] = domin;
            ngf[4 * nbfacesnetg + 3] = domout;
            label = nbfacesnetg + 1;
            nbfacesnetg++;
            if (nbfacesnetg == 10000) {
              cerr << "numb" << endl;
              exit(1);
            }
          }

          infile2 >> nep;
          if (nep != 3) {
            cerr << "only triangle elements are considered in Freefem++" << endl;
          }

          nep = 3;

          int iv[3], infoP[3];

          for (j = 0; j < nep; j++) {
            infile2 >> iv[j];    // tri.PNum(j);
          }

          for (j = 0; j < nep; j++) {
            infile2 >> infoP[j];    // tri.PNum(j); // No Need By Freefem++
          }

          for (j = 0; j < nep; j++) {
            iv[j]--;
          }

          (bb++)->set(v, iv, label);    // add boundary
        }

        delete[] ngf;
      }

      if (flagsnewlabelsurf == 2) {
        Netgen_Face *debutpngdf;
        Netgen_Face *iteratorpngdf;
        int nbfacesnetg = 0;

        for (i = 1; i <= n; i++) {
          int j;
          int surfnr, bcp, domin, domout, nep, faceind = 0;
          infile2 >> surfnr >> bcp >> domin >> domout;

          int label;

          label = -1;
          iteratorpngdf = debutpngdf;

          for (j = 0; j < nbfacesnetg; j++) {
            if (iteratorpngdf->Test_Face(surfnr, domin, domout, bcp) == 1) {
              label = j + 1;
            }

            iteratorpngdf = iteratorpngdf->Next_Face( );
          }

          if (label == -1) {
            if (nbfacesnetg == 0) {
              // debutpngdf->Initialisation_Face( surfnr, domin, domout, bcp );
              debutpngdf = new Netgen_Face(surfnr, domin, domout, bcp);
            } else {
              // iteratorpngdf = debutpngdf;
              Netgen_Face *bb = debutpngdf;

              for (j = 0; j < nbfacesnetg; j++) {
                iteratorpngdf = bb;
                bb = iteratorpngdf->Next_Face( );
              }

              bb = new Netgen_Face(surfnr, domin, domout, bcp);
              iteratorpngdf->Initialisation_Next_Face(bb);
            }

            label = nbfacesnetg + 1;
            nbfacesnetg++;
          }

          infile2 >> nep;
          if (nep != 3) {
            cerr << "only triangle elements are considered in Freefem++" << endl;
          }

          nep = 3;

          int iv[3], infoP[3];

          for (j = 0; j < nep; j++) {
            infile2 >> iv[j];    // tri.PNum(j);
          }

          for (j = 0; j < nep; j++) {
            infile2 >> infoP[j];    // tri.PNum(j); // No Need By Freefem++
          }

          for (j = 0; j < nep; j++) {
            iv[j]--;
          }

          (bb++)->set(v, iv, label);    // add boundary
        }
      }
    }

    if (strcmp(str, "volumeelements") == 0) {
      infile2 >> n;
      nt = n;
      if (nt != 0) {
        for (i = 0; i < n; i++) {
          // Element el;
          int label, nep;
          int iv[4];

          infile2 >> label;
          infile2 >> nep;
          if (nep != 4) {
            cerr << "freefem++ doesn't support second order element" << endl;
          }

          for (int j = 0; j < nep; j++) {
            infile2 >> iv[j];
          }

          for (int j = 0; j < nep; j++) {
            iv[j]--;
          }

          (tt++)->set(v, iv, label);    // add element
        }
      }
    }

    if (strcmp(str, "edgesegments") == 0) {
      cout << "this parameter is not taken in consideration in freefem++" << endl;
    }

    if (strcmp(str, "edgesegmentsgi") == 0) {
      cout << "this parameter is not taken in consideration in freefem++" << endl;
    }

    if (strcmp(str, "edgesegmentsgi2") == 0) {
      cout << "this parameter is not taken in consideration in freefem++" << endl;
    }

    if (strcmp(str, "points") == 0) {
      cout << "this parameter was taken in consideration before" << endl;
    }

    if (strcmp(str, "identifications") == 0) {
      cout << "this parameter is not taken in consideration in freefem++" << endl;
    }

    if (strcmp(str, "identificationtypes") == 0) {
      cout << "this parameter is not taken in consideration in freefem++" << endl;
    }

    if (strcmp(str, "materials") == 0) {
      cout << "this parameter is not taken in consideration in freefem++" << endl;
    }

    if (strcmp(str, "bcnames") == 0) {
      cout << "this parameter is not taken in consideration in freefem++" << endl;
    }

    if (strcmp(str, "singular_points") == 0) {
      cout << "this parameter is not taken in consideration in freefem++" << endl;
    }

    if (strcmp(str, "singular_edge_left") == 0) {
      cout << "this parameter is not taken in consideration in freefem++" << endl;
    }

    if (strcmp(str, "singular_edge_right") == 0) {
      cout << "this parameter is not taken in consideration in freefem++" << endl;
    }

    if (strcmp(str, "singular_face_inside") == 0) {
      cout << "this parameter is not taken in consideration in freefem++" << endl;
    }

    if (strcmp(str, "singular_face_outside") == 0) {
      cout << "this parameter is not taken in consideration in freefem++" << endl;
    }

    if (strcmp(str, "endmesh2") == 0) {
      endmesh2 = true;
    }

    strcpy(str, "");
  }

  infile2.close( );

  if (nt == 0) {
    Mesh3 *Th3 = new Mesh3(nv, nbe, v, b);
    return Th3;
  } else {
    Mesh3 *Th3 = new Mesh3(nv, nt, nbe, v, t, b);
    return Th3;
  }
}

AnyType Netgen_LoadMesh_Op::operator( )(Stack stack) const {
  string *pffname = GetAny< string * >((*filename)(stack));
  int renumsurf = 0;

  if (nargs[1]) {
    renumsurf = GetAny< long >((*nargs[1])(stack));
  }

  assert(renumsurf <= 1 && renumsurf >= 0);

  Mesh3 *Th3_t = NETGEN_Load(*pffname, renumsurf);

  if (Th3_t->nt == 0) {
    Th3_t->BuildSurfaceAdj( );
  }

  Th3_t->BuildGTree( );
  Add2StackOfPtr2FreeRC(stack, Th3_t);

  return Th3_t;
}

static void Load_Init( ) {    // le constructeur qui ajoute la fonction "splitmesh3"  a freefem++
  if (verbosity) {
    cout << " load: netgen  " << endl;
  }

  Global.Add("netg", "(", new RemplissageNetgen);
  Global.Add("netgstl", "(", new Netgen_STL);
  if (verbosity) {
    cout << " load: netgloadmesh  " << endl;
  }

  Global.Add("netgload", "(", new Netgen_LoadMesh);
  if (verbosity) {
    cout << " load: netgload  " << endl;
  }
}

LOADFUNC(Load_Init)
