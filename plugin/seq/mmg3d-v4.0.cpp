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
// Frederic Hecht
// E-MAIL  : jacques.morice@ann.jussieu.fr
// frederic.hecht@sorbonne-universite.fr

/* clang-format off */
//ff-c++-LIBRARY-dep: mmg3d-v4
//ff-c++-cpp-dep:
/* clang-format on */

/*
 * Thank to the ARN ()  FF2A3 grant
 * ref:ANR-07-CIS7-002-01
 */

// ./ff-c++ mmg3dv4.cpp  -I../3rdparty/include/mmg3d/  -lmmg3d4

#include "ff++.hpp"
#include "msh3.hpp"
// #define ADAPTLIBRARY"
#include "dataff.h"

using namespace Fem2D;
// using namespace  mmg3d;

inline void add2(int *k, int n, int a) {
  for (int i = 0; i < n; ++i) {
    k[i] += a;
  }
}

void set_mesh(void *dataff, int *data, int ldata) {
  DataFF *dff = (DataFF *)dataff;
  int nnv = data[ff_id_vertex];
  int nnt = data[ff_id_tet];
  int nnbe = data[ff_id_tria];
  Vertex3 *vv = new Vertex3[nnv];
  Tet *tt = new Tet[nnt];
  Triangle3 *bb = new Triangle3[nnbe];
  Mesh3 *pTh = new Mesh3( );

  pTh->nv = nnv;
  pTh->nt = nnt;
  pTh->nbe = nnbe;

  pTh->vertices = vv;
  pTh->elements = tt;
  pTh->borderelements = bb;

  pTh->mes = 0.;
  pTh->mesb = 0.;

  dff->mesh = pTh;
  if (verbosity > 5) {
    cout << " Set_mesh nv=" << nnv << " nTet " << nnt << " NTria " << nnbe << endl;
  }
}

void end_mesh(void *dataff) {
  DataFF *dff = (DataFF *)dataff;
  Mesh3 &Th = *(Mesh3 *)dff->mesh;

  Th.mes = 0.;
  Th.mesb = 0.;

  for (int i = 0; i < Th.nbe; i++) {
    Th.mesb += Th.be(i).mesure( );
  }

  // Add FH to be consitant we all constructor ...  July 09
  Th.BuildBound( );
  // Th.Save("TTTh.mesh");
  if (verbosity > 5) {
    cout << "end_mesh:  Th.mes = " << Th.mes << " Th.mesb = " << Th.mesb << endl;
  }

  // if (Th.nt > 0) {
  Th.BuildAdj( );
  Th.Buildbnormalv( );
  Th.BuildjElementConteningVertex( );
  //}

  // end add

  if (verbosity > 1) {
    cout << "  -- End of Construct  mesh3: mesure = " << Th.mes << " border mesure " << Th.mesb
         << endl;
  }

  ffassert(Th.mes >= 0);    // add F. Hecht sep 2009.
}

void set_v(void *dataff, int i, double *xyz, int lab) {
  i--;
  DataFF *dff = (DataFF *)dataff;
  ffassert(dff->mesh);
  Mesh3 &Th = *(Mesh3 *)dff->mesh;
  Th.vertices[i].x = xyz[0];
  Th.vertices[i].y = xyz[1];
  Th.vertices[i].z = xyz[2];
  Th.vertices[i].lab = lab;
  if (verbosity > 10) {
    cout << " set_v3 " << i << " " << xyz[0] << " " << xyz[1] << " " << xyz[02] << " " << lab
         << endl;
  }
}

void set_elmt(void *dataff, int id, int i, int *k, int lab) {
  i--;
  int n = 0;
  DataFF *dff = (DataFF *)dataff;
  ffassert(dff->mesh);
  Mesh3 &Th = *(Mesh3 *)dff->mesh;
  if (id == 2) {
    n = 3;
    Mesh3::BorderElement &K(Th.be(i));
    add2(k, 3, -1);
    K.set(Th.vertices, k, lab);
  } else if (id == 3) {
    n = 4;
    Mesh3::Element &K(Th.t(i));
    add2(k, 4, -1);
    K.set(Th.vertices, k, lab);
  } else {
    cout << " unknows id = " << id << " not 2 or 3 " << endl;
    ffassert(0);
  }

  if (verbosity > 10) {
    cout << " set_ele" << n << " " << i << " ";

    for (int j = 0; j < n; j++) {
      cout << k[j] << " ";
    }

    cout << lab << endl;
  }
}

void get_mesh(void *dataff, int *data, int ldata) {
  DataFF *dff = (DataFF *)dataff;

  assert(ldata > 5);

  for (int i = 0; i < ldata; ++i) {
    data[i] = 0;
  }

  ffassert(dff->mesh);
  Mesh3 &Th = *(Mesh3 *)dff->mesh;
  data[ff_id_vertex] = Th.nv;
  data[ff_id_tria] = Th.nbe;
  data[ff_id_tet] = Th.nt;
  if (verbosity > 9) {
    cout << " get_mesh " << Th.nv << " " << Th.nbe << " " << Th.nt << endl;
  }
}

void get_v3(void *dataff, int i, double *xyz, int *lab) {
  i--;
  DataFF *dff = (DataFF *)dataff;
  ffassert(dff->mesh);
  Mesh3 &Th = *(Mesh3 *)dff->mesh;
  xyz[0] = Th.vertices[i].x;
  xyz[1] = Th.vertices[i].y;
  xyz[2] = Th.vertices[i].z;
  *lab = Th.vertices[i].lab;
  if (verbosity > 10) {
    cout << " get_v3 " << i << " " << xyz[0] << " " << xyz[1] << " " << xyz[2] << " " << *lab
         << endl;
  }
}

void get_elmt(void *dataff, int id, int i, int *k, int *lab) {
  i--;
  DataFF *dff = (DataFF *)dataff;
  ffassert(dff->mesh);
  Mesh3 &Th = *(Mesh3 *)dff->mesh;
  int n = 0;
  if (id == 2) {
    n = 3;
    Mesh3::BorderElement &K(Th.be(i));

    for (int j = 0; j < n; ++j) {
      k[j] = Th(K[j]);
    }

    *lab = K.lab;
  } else if (id == 3) {
    n = 4;
    Mesh3::Element &K(Th.t(i));

    for (int j = 0; j < n; ++j) {
      k[j] = Th(K[j]);
    }

    *lab = K.lab;
  } else {
    cout << "  id != 2, 3 , id = = " << id << endl;
    ffassert(0);
  }

  add2(k, n, +1);

  if (verbosity > 10) {
    cout << " get_ele" << n << " " << i << " ";

    for (int j = 0; j < n; j++) {
      cout << k[j] << " ";
    }

    cout << *lab << endl;
  }
}

class mmg3d_Op : public E_F0mps {
 public:
  Expression eTh, xx, yy, zz;
  static const int n_name_param = 5;    //
  static basicAC_F0::name_and_type name_param[];
  Expression nargs[n_name_param];
  KN_< long > karg(int i, Stack stack) const {
    return nargs[i] ? GetAny< KN_< long > >((*nargs[i])(stack)) : KN_< long >((long *)0, 0L);
  }

  double arg(int i, Stack stack, double a) const {
    return nargs[i] ? GetAny< double >((*nargs[i])(stack)) : a;
  }

  string arg(int i, Stack stack, const char *a) const {
    return nargs[i] ? *GetAny< string * >((*nargs[i])(stack)) : a;
  }

  long arg(int i, Stack stack, long a) const {
    return nargs[i] ? GetAny< long >((*nargs[i])(stack)) : a;
  }

 public:
  mmg3d_Op(const basicAC_F0 &args, Expression tth) : eTh(tth), xx(0), yy(0), zz(0) {
    if (verbosity > 1) {
      cout << "mmg3d v4 " << endl;
    }

    args.SetNameParam(n_name_param, name_param, nargs);

    const E_Array *a1 = 0;
    if (nargs[1]) {
      a1 = dynamic_cast< const E_Array * >(nargs[1]);
    }

    if (a1) {
      if (a1->size( ) != 3) {
        CompileError("mmg3d(Th,displacement=[X,Y,Z],) ");
      }

      xx = to< double >((*a1)[0]);
      yy = to< double >((*a1)[1]);
      zz = to< double >((*a1)[2]);
    } else if (nargs[1]) {
      CompileError("mmg3d(Th,displacement=[X,Y,Z], .... ) ");
    }
  }

  AnyType operator( )(Stack stack) const;
};

basicAC_F0::name_and_type mmg3d_Op::name_param[] = {
  {"metric", &typeid(KN< double > *)},      // 0
  {"displacement", &typeid(E_Array)},       // 1
  {"displVect", &typeid(KN_< double >)},    // 2
  {"opt", &typeid(string *)},               // 3
  {"Mb", &typeid(long)}                     // 4
};

class mmg3d_ff : public OneOperator {
 public:
  mmg3d_ff( ) : OneOperator(atype< pmesh3 >( ), atype< pmesh3 >( )) {}

  E_F0 *code(const basicAC_F0 &args) const { return new mmg3d_Op(args, t[0]->CastTo(args[0])); }
};

AnyType mmg3d_Op::operator( )(Stack stack) const {
  // initialisation
  MeshPoint *mp(MeshPointStack(stack)), mps = *mp;
  Mesh3 *pTh = GetAny< Mesh3 * >((*eTh)(stack));

  ffassert(pTh);
  Mesh3 &Th3 = *pTh;
  string sarg = arg(3, stack, "");
  DataFF dff;
  dff.memory = arg(4, stack, 128L);    // 128 Mb .. ????
  ffassert(dff.memory < 2048);         // 2 GiGa bytes   limite of integer ..
  dff.typesol = 0;
  dff.np = pTh->nv;
  dff.mesh = pTh;
  dff.meshname = "Th";
  dff.imprim = verbosity;
  dff.sol = 0;
  dff.mov = 0;
  dff.set_mesh = set_mesh;
  dff.end_mesh = end_mesh;
  dff.set_v = set_v;
  dff.set_elmt = set_elmt;
  dff.get_mesh = get_mesh;
  dff.get_v3 = get_v3;
  dff.get_elmt = get_elmt;

  KN< double > *pmetric = 0;

  if (nargs[0]) {
    pmetric = GetAny< KN< double > * >((*nargs[0])(stack));
    ffassert(pmetric);
  }

  KN< double > cmetric;
  if (pmetric) {
    int m = pmetric->N( );
    if (m == Th3.nv * 6) {    //
      cmetric = (*pmetric);
      dff.typesol = 6;
      dff.np = Th3.nv;

      for (int i = 0; i < m; i += 6) {
        std::swap(cmetric[i + 2], cmetric[i + 3]);
      }

      dff.sol = &cmetric[0];
      dff.solname = "metrix-aniso";
    } else if (m == Th3.nv) {
      dff.typesol = 1;
      dff.np = Th3.nv;
      dff.sol = &(*pmetric)[0];
      dff.solname = "metrix-iso";
    } else {
      ExecError(" mmg3d v4: incompatibility  metric array  mesh ");
    }
  }

  bool BoolMoving = 0;
  KN< double > Moving(0);

  if (nargs[1] || nargs[2]) {
    BoolMoving = 1;
    if (nargs[2]) {
      Moving = GetAny< double >((*nargs[2])(stack));
      assert(Moving.N( ) == 3 * Th3.nv);
      dff.movename = "move";
      dff.mov = Moving;
      if (Moving.N( ) != 3 * Th3.nv) {
        cerr << " Displacement vector is of size 3*Th.nv" << endl;
        ExecError(" mmg3d v4");
      }
    } else {
      MeshPoint *mp3(MeshPointStack(stack));
      Moving.resize(3 * Th3.nv);

      for (int i = 0; i < Th3.nv; ++i) {
        mp3->set(Th3.vertices[i].x, Th3.vertices[i].y, Th3.vertices[i].z);
        if (xx) {
          Moving[3 * i] = GetAny< double >((*xx)(stack));
        }

        if (yy) {
          Moving[3 * i + 1] = GetAny< double >((*yy)(stack));
        }

        if (zz) {
          Moving[3 * i + 2] = GetAny< double >((*zz)(stack));
        }
      }

      dff.mov = Moving;
      dff.movename = "move";
    }

    if (verbosity > 2) {
      cout << "displacement vector is loading" << endl;
    }
  }

  int argc = 1;
  char *argv[1000];
  char ff[10] = "ff++";
  KN< char > args(sarg.size( ) + 1);
  argv[0] = ff;
  argv[1] = &args[0];
  strcpy(args, sarg.c_str( ));
  char cc = '\0';

  for (int i = 0; i < args.N( ); cc = args[i++]) {
    if (isspace(args[i]) && cc != '\\') {
      cc = args[i] = '\0';
    } else if (!cc) {
      argv[argc++] = &args[i];
    }

    ffassert(argc < 1000);
  }

  int res = mainmmg3d(argc, argv, &dff);
  Mesh3 *pTh3 = (Mesh3 *)dff.mesh;
  if (res > 0) {
    dff.mesh = 0;
    cout << " problem of remeshing with mmg3d :: error" << res << endl;
  }

  if (!pTh3) {
    cout << " problem of remeshing with mmg3d v 4 (no mesh)  :: error" << res << endl;
    ExecError(" Error mmg3d");
  } else {
    // end build of TH3...
    if (verbosity > 10) {
      cout << "buildGtree" << endl;
    }

    pTh3->BuildGTree( );
  }

  *mp = mps;
  Add2StackOfPtr2FreeRC(stack, pTh3);
  return pTh3;
}

static void Load_Init( ) {    // le constructeur qui ajoute la fonction "splitmesh3"  a freefem++
  // if (verbosity)
  if (verbosity) {
    cout << " load: mmg3d  " << endl;
  }

  Global.Add("mmg3d", "(", new mmg3d_ff);
}

#define WITH_NO_INIT
#include "msh3.hpp"
LOADFUNC(Load_Init)
