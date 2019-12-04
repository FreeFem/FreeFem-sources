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
//ff-c++-LIBRARY-dep: mmg3d
//ff-c++-cpp-dep:
/* clang-format on */

/*
 * Thank to the ARN ()  FF2A3 grant
 * ref:ANR-07-CIS7-002-01
 */

// ./ff-c++ mmg3d.cpp -I../src/libMesh/ -I../3rdparty/include/mmg3d/ -L../3rdparty/lib/mmg3d/
// -lmmg3d -L../src/libMesh/ -lMesh

#include "ff++.hpp"
#include "msh3.hpp"
#include "memory.h"
#include "libmmg3d.h"
#include "mmg3d_defines.h"

using namespace Fem2D;
using namespace mmg3d;
extern "C" {
int MMG_markBdry(MMG_pMesh);
}
Mesh3 *MMG_pMesh_to_msh3(MMG_pMesh meshMMG) {
  int i;

  // determination of exact :: vertices, triangles, tetrahedrons
  int nv = 0, nt = 0, nbe = 0;
  MMG_pPoint ppt;
  MMG_pTetra ptetra;
  MMG_pTria ptriangle;

  for (int k = 1; k <= meshMMG->np; k++) {
    ppt = &meshMMG->point[k];
    if (ppt->tag & M_UNUSED) {
      continue;
    }

    ppt->tmp = nv++;
  }

  nt = 0;

  for (int k = 1; k <= meshMMG->ne; k++) {
    ptetra = &meshMMG->tetra[k];
    if (!ptetra->v[0]) {
      continue;
    }

    nt++;
  }

  meshMMG->nt = 0;
  if (MMG_markBdry(meshMMG)) {
    nbe = meshMMG->nt;
  }

  if (verbosity > 1) {
    cout << "transformation maillage --> msh3 " << endl;
    cout << "meshMMG->np =" << meshMMG->np << " -> nv  = " << nv << endl;
    cout << "meshMMG->ne =" << meshMMG->ne << "  -> nt  = " << nt << endl;
    cout << "meshMMG->nt =" << meshMMG->nt << "  -> nbe = " << nbe << endl;
  }

  Vertex3 *v = new Vertex3[nv];
  Tet *t = new Tet[nt];
  Tet *tt = t;
  Triangle3 *b = new Triangle3[nbe];
  Triangle3 *bb = b;
  int k;
  Vertex3 *v3 = v;

  for (k = 1; k <= meshMMG->np; k++) {
    ppt = &meshMMG->point[k];
    if (ppt->tag & M_UNUSED) {
      continue;
    }

    v3->x = ppt->c[0];
    v3->y = ppt->c[1];
    v3->z = ppt->c[2];
    v3->lab = ppt->ref;
    v3++;
  }

  for (k = 1; k <= meshMMG->ne; k++) {
    int iv[4], lab;
    ptetra = &meshMMG->tetra[k];
    if (!ptetra->v[0]) {
      continue;
    }

    for (int i = 0; i < 4; ++i) {
      iv[i] = meshMMG->point[ptetra->v[i]].tmp;
    }

    lab = ptetra->ref;
    tt++->set(v, iv, lab);
  }

  for (k = 1; k <= meshMMG->nt; k++) {
    int iv[3], lab;
    ptriangle = &meshMMG->tria[k];

    for (int i = 0; i < 3; ++i) {
      iv[i] = meshMMG->point[ptriangle->v[i]].tmp;
    }

    lab = ptriangle->ref;
    bb++->set(v, iv, lab);
  }

  Mesh3 *T_TH3 = new Mesh3(nv, nt, meshMMG->nt, v, t, b);
  if (verbosity > 1) {
    cout << "transformation maillage --> msh3 " << endl;
    cout << "vertices =" << nv << endl;
    cout << "tetrahedrons =" << nt << endl;
    cout << "triangles =" << meshMMG->nt << endl;
    cout << "T_TH3" << T_TH3->nv << " " << T_TH3->nt << " " << T_TH3->nbe << endl;
  }

  return T_TH3;
}

MMG_pMesh mesh3_to_MMG_pMesh(const Mesh3 &Th3, const int &nvmax, const int &ntrimax,
                             const int &ntetmax, const bool &boolMoving, KN< double > &Moving) {
  MMG_pMesh meshMMG;

  meshMMG = (MMG_pMesh)calloc(1, sizeof(MMG_Mesh));

  meshMMG->np = Th3.nv;
  meshMMG->nt = Th3.nbe;
  meshMMG->ne = Th3.nt;

  meshMMG->npmax = nvmax;
  meshMMG->ntmax = ntrimax;
  meshMMG->nemax = ntetmax;

  meshMMG->point = (MMG_pPoint)calloc(meshMMG->npmax + 1, sizeof(MMG_Point));
  meshMMG->tetra = (MMG_pTetra)calloc(meshMMG->nemax + 1, sizeof(MMG_Tetra));
  meshMMG->tria = (MMG_pTria)calloc(meshMMG->ntmax + 1, sizeof(MMG_Tria));

  if (boolMoving) {
    MMG_pDispl ppd;
    meshMMG->disp = (MMG_pDispl)calloc(1, sizeof(MMG_Displ));
    ppd = meshMMG->disp;
    ppd->np = meshMMG->np;
    ppd->mv = (double *)calloc(3 * (meshMMG->np + 1), sizeof(double));
    assert(meshMMG->disp->mv);
    ppd->alpha = (short *)calloc(meshMMG->np + 1, sizeof(short));
    assert(meshMMG->disp->alpha);

    for (int ii = 0; ii < meshMMG->np; ii++) {
      ppd->mv[3 * ii + 1] = Moving[3 * ii];
      ppd->mv[3 * ii + 2] = Moving[3 * ii + 1];
      ppd->mv[3 * ii + 3] = Moving[3 * ii + 2];
    }
  } else {
    meshMMG->disp = NULL;
  }

  meshMMG->adja = (int *)calloc(4 * meshMMG->nemax + 5, sizeof(int));

  int k;
  MMG_pPoint ppt;

  for (k = 1; k <= meshMMG->np; k++) {
    ppt = &meshMMG->point[k];
    ppt->c[0] = Th3.vertices[k - 1].x;
    ppt->c[1] = Th3.vertices[k - 1].y;
    ppt->c[2] = Th3.vertices[k - 1].z;
    ppt->ref = Th3.vertices[k - 1].lab;
  }

  MMG_pTetra ptetra;

  for (k = 1; k <= meshMMG->ne; k++) {
    const Tet &K(Th3.elements[k - 1]);
    ptetra = &meshMMG->tetra[k];
    ptetra->v[0] = Th3.operator( )(K[0]) + 1;
    ptetra->v[1] = Th3.operator( )(K[1]) + 1;
    ptetra->v[2] = Th3.operator( )(K[2]) + 1;
    ptetra->v[3] = Th3.operator( )(K[3]) + 1;
    ptetra->ref = K.lab;
  }

  MMG_pTria ptriangle;

  for (k = 1; k <= meshMMG->nt; k++) {
    const Triangle3 &K(Th3.be(k - 1));
    ptriangle = &meshMMG->tria[k];
    ptriangle->v[0] = Th3.operator( )(K[0]) + 1;
    ptriangle->v[1] = Th3.operator( )(K[1]) + 1;
    ptriangle->v[2] = Th3.operator( )(K[2]) + 1;
    ptriangle->ref = K.lab;
  }
  return meshMMG;
}

MMG_pSol metric_mmg3d(const int &nv, const int &nvmax, const KN< double > *pmetric) {
  static const int wrapperMetric[6] = {0, 1, 3, 2, 4, 5};
  MMG_pSol sol;

  sol = (MMG_pSol)calloc(1, sizeof(MMG_Sol));

  if (pmetric && pmetric->N( ) > 0) {
    const KN< double > &metric = *pmetric;
    sol->np = nv;
    sol->npmax = nvmax;

    if (metric.N( ) == nv) {
      const int ic = 1;
      char newvalue[sizeof(int)];
      sprintf(newvalue, "%s", (char *)&ic);
      sol->offset = *newvalue;
    } else {
      const int ic = 6;
      char newvalue[sizeof(int)];
      sprintf(newvalue, "%s", (char *)&ic);
      sol->offset = *newvalue;
    }

    sol->met = (double *)calloc(sol->npmax + 1, (int)sol->offset * sizeof(double));
    // need to have anisotropic meshes
    sol->metold = (double *)calloc(sol->npmax + 1, (int)sol->offset * sizeof(double));

    int k, isol, i;
    double tmp;

    for (k = 1; k <= sol->np; k++) {
      isol = (k - 1) * sol->offset + 1;

      for (i = 0; i < sol->offset; i++) {
        sol->met[isol + i] = metric[(isol - 1) + i];
      }

      // On a besoin de swap car on utilise celui donner par mshmet
      if (sol->offset == 6) {
        // MMG_swap data
        tmp = sol->met[isol + 2];
        sol->met[isol + 2] = sol->met[isol + 3];
        sol->met[isol + 3] = tmp;
      }
    }
  } else {
    sol->met = 0;
    sol->metold = 0;
    sol->np = 0;
    sol->offset = 1;
  }

  cout << " sol->offset " << (int)sol->offset << endl;
  return sol;
}

void metric_mmg3d_to_ff_metric(MMG_pSol sol, KN< double > &metric) {
  static const int wrapperMetric[6] = {0, 1, 3, 2, 4, 5};
  int k, isol, jsol, i;

  for (k = 1; k <= sol->np; k++) {
    isol = (k - 1) * sol->offset + 1;
    jsol = (wrapperMetric[k] - 1) * sol->offset;

    for (i = 0; i < sol->offset; i++) {
      metric[jsol + i] = sol->met[isol + i];
    }
  }
}

class mmg3d_Op : public E_F0mps {
 public:
  Expression eTh, xx, yy, zz;
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

  long arg(int i, Stack stack, long a) const {
    return nargs[i] ? GetAny< long >((*nargs[i])(stack)) : a;
  }

 public:
  mmg3d_Op(const basicAC_F0 &args, Expression tth) : eTh(tth) {
    if (verbosity > 1) {
      cout << "mmg3d" << endl;
    }

    args.SetNameParam(n_name_param, name_param, nargs);

    const E_Array *a1 = 0;
    if (nargs[2]) {
      a1 = dynamic_cast< const E_Array * >(nargs[2]);
    }

    if (a1) {
      if (a1->size( ) != 3) {
        CompileError("mmg3d(Th,displacement=[X,Y,Z],) ");
      }

      xx = to< double >((*a1)[0]);
      yy = to< double >((*a1)[1]);
      zz = to< double >((*a1)[2]);
    }
  }

  AnyType operator( )(Stack stack) const;
};

basicAC_F0::name_and_type mmg3d_Op::name_param[] = {
  {"options", &typeid(KN_< long >)},        // 3 -> 0
  {"metric", &typeid(KN< double > *)},      // 4 -> 1
  {"displacement", &typeid(E_Array)},       // 5 -> 2
  {"displVect", &typeid(KN_< double >)},    // 6 ->3
  {"memory", &typeid(long)}                 // 7->4
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
  int nv = Th3.nv;
  int nt = Th3.nt;
  int nbe = Th3.nbe;

  // default value for max nv,ntri,ntet  max
  // division of tetrahedrons in the middle of edges
  int defaultnvmax = 500000;
  int defaultntrimax = 1000000;
  int defaultntetmax = 3000000;

  KN< long > defaultopt(9);
  /*
   * defaultopt(0)= 1;
   * defaultopt(1)= 0;
   * defaultopt(2)= 64;
   * defaultopt(3)= 0;
   * defaultopt(4)= 0;
   * defaultopt(5)= 3;
   */
  defaultopt[0] = 4;      // splitting
  defaultopt[1] = 0;      // debug
  defaultopt[2] = 64;     // par default 64
  defaultopt[3] = 0;      // noswap
  defaultopt[4] = 0;      // noinsert
  defaultopt[5] = 0;      // nomove
  defaultopt[6] = 5;      // imprim
  defaultopt[7] = 3;      // renum (not use scotch ..)
  defaultopt[8] = 500;    // renum (not use scotch ..)

  int nvmax;      // (arg(0,stack,-1L));
  int ntrimax;    // (arg(1,stack,-1L));
  int ntetmax;    // (arg(2,stack,-1L));
  KN< int > opt(arg(0, stack, defaultopt));

  int memory(arg(4, stack, -1L));

  cout << "memory =" << memory << "Mb " << endl;
  if (memory < 0) {
    nvmax = max((int)1.5 * nv, defaultnvmax);
    ntrimax = max((int)1.5 * nbe, defaultntrimax);
    ntetmax = max((int)1.5 * nt, defaultntetmax);
  } else {
    int million = 1048576L;
    int bytes = sizeof(MMG_Point) + 0.2 * sizeof(MMG_Tria) + 6 * sizeof(MMG_Tetra) +
                4 * sizeof(int) + sizeof(MMG_Sol) + sizeof(MMG_Displ) + sizeof(int) +
                5 * sizeof(int);
    int npask = (double)memory / bytes * million;
    nvmax = max((int)1.5 * nv, npask);
    ntetmax = max((int)1.5 * nt, 6 * npask);
    ntrimax = max((int)1.5 * nbe, (int)(0.3 * npask));

    if (verbosity > 10) {
      cout << " mmg3d : npask=" << npask << endl;
      cout << "      memory is given " << endl;
      cout << "        nvmax " << nvmax << endl;
      cout << "        ntrimax " << ntrimax << endl;
      cout << "        ntetmax " << ntetmax << endl;
    }
  }

  if (verbosity > 10) {
    cout << " mmg3d: nvmax " << nvmax << endl;
    cout << "        ntrimax " << ntrimax << endl;
    cout << "        ntetmax " << ntetmax << endl;
  }

  KN< double > *metric = 0;

  if (nargs[1]) {
    metric = GetAny< KN< double > * >((*nargs[1])(stack));
  }

  bool BoolMoving = 0;
  KN< double > Moving(0);

  if (nargs[2] || nargs[3]) {
    BoolMoving = 1;
    if (nargs[3]) {
      Moving = GetAny< double >((*nargs[3])(stack));
      assert(Moving.N( ) == 3 * Th3.nv);
      if (Moving.N( ) != 3 * Th3.nv) {
        cerr << " Displacement vector is of size 3*Th.nv" << endl;
        exit(1);
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
    }

    // if(verbosity > 2)
    if (verbosity > 2) {
      cout << "displacement vector is loading" << endl;
    }
  }

  MMG_pMesh MMG_Th3 = mesh3_to_MMG_pMesh(Th3, nvmax, ntrimax, ntetmax, BoolMoving, Moving);
  MMG_pSol sol = metric_mmg3d(nv, nvmax, metric);
  int res = MMG_mmg3dlib(opt, MMG_Th3, sol);

  if (res > 0) {
    cout << " problem of remeshing with mmg3d :: error" << res << endl;
    free(MMG_Th3->point);
    free(MMG_Th3->tria);
    free(MMG_Th3->tetra);
    /*la desallocation de ce pointeur plante dans certains cas...*/
    if (verbosity > 10) {
      cout << "mesh: adja" << endl;
    }

    free(MMG_Th3->adja);
    if (verbosity > 10) {
      cout << "mesh: disp" << endl;
    }

    if (BoolMoving) {
      free(MMG_Th3->disp->alpha);
      free(MMG_Th3->disp->mv);
    }

    free(MMG_Th3->disp);
    free(MMG_Th3);

    free(sol->met);
    if (abs(opt[0]) != 9) {
      free(sol->metold);    // 9 -> free in mmg3dlib.c
    }

    free(sol);

    ExecError("mmg3d ??? ");
  }

  Mesh3 *Th3_T = MMG_pMesh_to_msh3(MMG_Th3);
  // return new metric in the parameter nargs
  int nvv = Th3_T->nv;    // number of vertices
  if (metric) {
    metric->resize(sol->offset * nvv);
    if (verbosity > 1) {
      cout << " sol->met " << nvv * sol->offset << endl;
    }

    {
      int k, isol, i;
      MMG_pPoint ppt;

      for (k = 1; k <= MMG_Th3->np; k++) {
        ppt = &MMG_Th3->point[k];
        if (ppt->tag & M_UNUSED) {
          continue;
        }

        isol = (k - 1) * sol->offset + 1;

        for (i = 0; i < sol->offset; i++) {
          (*metric)[(isol - 1) + i] = sol->met[isol + i];
        }
      }
    }
  }

  if (verbosity > 10) {
    cout << "buildGtree" << endl;
  }

  Th3_T->BuildGTree( );

  /* free mem */
  if (verbosity > 10) {
    cout << "mesh" << endl;
  }

  free(MMG_Th3->point);
  free(MMG_Th3->tria);
  free(MMG_Th3->tetra);
  /*la desallocation de ce pointeur plante dans certains cas...*/
  if (verbosity > 10) {
    cout << "mesh: adja" << endl;
  }

  free(MMG_Th3->adja);
  if (verbosity > 10) {
    cout << "mesh: disp" << endl;
  }

  if (BoolMoving) {
    free(MMG_Th3->disp->alpha);
    free(MMG_Th3->disp->mv);
  }

  free(MMG_Th3->disp);
  free(MMG_Th3);

  free(sol->met);
  if (abs(opt[0]) != 9) {
    free(sol->metold);    // 9 -> free in mmg3dlib.c
  }

  free(sol);

  *mp = mps;
  Add2StackOfPtr2FreeRC(stack, Th3_T);
  return Th3_T;
}

static void Load_Init( ) {    // le constructeur qui ajoute la fonction "splitmesh3"  a freefem++
  if (verbosity) {
    cout << " load: mmg3d  " << endl;
  }

  Global.Add("mmg3d", "(", new mmg3d_ff);
}

#define WITH_NO_INIT
#include "msh3.hpp"
LOADFUNC(Load_Init)
