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
//ff-c++-LIBRARY-dep: mshmet libMesh
//ff-c++-cpp-dep:
/* clang-format on */

/*
 * Thank to the ARN ()  FF2A3 grant
 * ref:ANR-07-CIS7-002-01
 */

// ./ff-c++ mshmet.cpp -I/Users/morice/Desktop/adaptmesh3d/mshmet.2009.09.16/sources
// -L/Users/morice/Desktop/adaptmesh3d/mshmet.2009.09.16/objects/i386/ -lmshmet
// -L/Users/morice/work/postdoc/freefem++prod/src/libMesh/ -lMesh

// ./ff-c++ mshmet.cpp -I../3rdparty/include/mshmet/ -L../3rdparty/lib/mshmet/ -lmshmet
// -L/Users/morice/work/postdoc/freefem++prod/src/libMesh/ -lMesh

#include "ff++.hpp"
#include "msh3.hpp"
#include "mshmetlib.h"

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include <assert.h>

/* prototype (re)definitions */
void *M_malloc(size_t size, char *call);
void *M_calloc(size_t nelem, size_t elsize, const char *call);
void *M_realloc(void *ptr, size_t size, const char *call);
void M_free(void *ptr);

/* ptototypes : tools */
int M_memLeak( );
void M_memDump( );
size_t M_memSize( );

#ifdef __cplusplus
}
#endif

using namespace Fem2D;
using namespace mshmet;
typedef const Mesh *pmesh;
typedef const Mesh3 *pmesh3;
// 2d mesh function
// Add FH thank to I. Bajc.  (iztok.bajc@fmf.uni-lj.si) 03/14
//
static void myMSHMET_free(MSHMET_pMesh mesh, MSHMET_pSol sol) {
  /* free mem */
  M_free(mesh->point);
  if (mesh->nt) {
    M_free(mesh->tria);
  }

  if (mesh->ne) {
    M_free(mesh->tetra);
  }

  M_free(mesh->adja);
  M_free(mesh);
  M_free(sol->sol);
  M_free(sol->met);
  M_free(sol);
}

MSHMET_pMesh mesh_to_MSHMET_pMesh(const Mesh &Th) {
  MSHMET_pMesh meshMSHMET;

  meshMSHMET = (MSHMET_pMesh)M_calloc(1, sizeof(MSHMET_Mesh), "mesh2");

  meshMSHMET->dim = 2;
  meshMSHMET->np = Th.nv;
  meshMSHMET->nt = Th.nt;
  meshMSHMET->ne = 0;

  meshMSHMET->point = (MSHMET_pPoint)M_calloc(meshMSHMET->np + 1, sizeof(MSHMET_Point), "point");
  meshMSHMET->tria = (MSHMET_pTria)M_calloc(meshMSHMET->nt + 1, sizeof(MSHMET_Tria), "tria");
  meshMSHMET->adja = (int *)M_calloc(3 * meshMSHMET->nt + 5, sizeof(int), "adja");

  int k;
  MSHMET_pPoint ppt;

  for (k = 1; k <= meshMSHMET->np; k++) {
    ppt = &meshMSHMET->point[k];
    ppt->c[0] = Th.vertices[k - 1].x;
    ppt->c[1] = Th.vertices[k - 1].y;
    ppt->c[2] = 0.;
  }

  MSHMET_pTria ptriangle;
  MSHMET_pPoint p0, p1, p2;
  double ux, uy, h1, h2, h3, pe, rins;
  int i;

  for (k = 1; k <= meshMSHMET->nt; k++) {
    const Mesh::Triangle &K(Th.t(k - 1));
    ptriangle = &meshMSHMET->tria[k];
    ptriangle->v[0] = Th.operator( )(K[0]) + 1;
    ptriangle->v[1] = Th.operator( )(K[1]) + 1;
    ptriangle->v[2] = Th.operator( )(K[2]) + 1;

    for (i = 0; i < 3; i++) {
      ppt = &meshMSHMET->point[ptriangle->v[i]];
      if (!ppt->s) {
        ppt->s = k;
      }
    }

    p0 = &meshMSHMET->point[ptriangle->v[0]];
    p1 = &meshMSHMET->point[ptriangle->v[1]];
    p2 = &meshMSHMET->point[ptriangle->v[2]];

    ux = p1->c[0] - p0->c[0];
    uy = p1->c[1] - p0->c[1];
    h1 = sqrt(ux * ux + uy * uy);

    ux = p2->c[0] - p0->c[0];
    uy = p2->c[1] - p0->c[1];
    h2 = sqrt(ux * ux + uy * uy);

    ux = p2->c[0] - p1->c[0];
    uy = p2->c[1] - p1->c[1];
    h3 = sqrt(ux * ux + uy * uy);
  }

  return meshMSHMET;
}

// 3d mesh function

MSHMET_pMesh mesh3_to_MSHMET_pMesh(const Mesh3 &Th3) {
  MSHMET_pMesh meshMSHMET;

  meshMSHMET = (MSHMET_pMesh)M_calloc(1, sizeof(MSHMET_Mesh), "Mesh3");

  meshMSHMET->dim = 3;
  meshMSHMET->np = Th3.nv;
  meshMSHMET->nt = 0;
  meshMSHMET->ne = Th3.nt;

  meshMSHMET->point = (MSHMET_pPoint)M_calloc(meshMSHMET->np + 1, sizeof(MSHMET_Point), "point3");
  meshMSHMET->tetra = (MSHMET_pTetra)M_calloc(meshMSHMET->ne + 1, sizeof(MSHMET_Tetra), "tetra");
  meshMSHMET->adja = (int *)M_calloc(4 * meshMSHMET->ne + 5, sizeof(int), "adja3");

  int k;
  MSHMET_pPoint ppt;

  for (k = 1; k <= meshMSHMET->np; k++) {
    ppt = &meshMSHMET->point[k];
    ppt->c[0] = Th3.vertices[k - 1].x;
    ppt->c[1] = Th3.vertices[k - 1].y;
    ppt->c[2] = Th3.vertices[k - 1].z;
  }

  int i;
  MSHMET_pTetra ptetra;

  for (k = 1; k <= meshMSHMET->ne; k++) {
    const Tet &K(Th3.elements[k - 1]);
    ptetra = &meshMSHMET->tetra[k];
    ptetra->v[0] = Th3.operator( )(K[0]) + 1;
    ptetra->v[1] = Th3.operator( )(K[1]) + 1;
    ptetra->v[2] = Th3.operator( )(K[2]) + 1;
    ptetra->v[3] = Th3.operator( )(K[3]) + 1;

    for (i = 0; i < 4; i++) {
      ppt = &meshMSHMET->point[ptetra->v[i]];
      if (meshMSHMET->dim == 3 && !ppt->s) {
        ppt->s = k;
      }
    }
  }

  return meshMSHMET;
}

MSHMET_pSol sol_mshmet(const int &dim, const int &np, const int &type, const int &size, int *typtab,
                       const KN< double > &solutions) {
  static const int wrapperMetric[6] = {0, 1, 3, 2, 4, 5};
  MSHMET_pSol sol;
  int k, ia, i;

  sol = (MSHMET_pSol)M_calloc(1, sizeof(MSHMET_Sol), "sol");
  sol->ver = 0;
  sol->np = np;
  sol->dim = dim;
  sol->type = type;    // nombre de solutions differentes
  sol->size = size;

  for (i = 0; i < sol->type; i++) {
    sol->typtab[i] = typtab[i];    // types des differentes solutions
  }

  sol->sol = (double *)M_calloc(sol->np + 1, sol->size * sizeof(double), "sol->sol");
  assert(sol->sol);

  for (k = 1; k <= sol->np; k++) {
    ia = (k - 1) * sol->size + 1;

    for (i = 0; i < sol->size; i++) {
      sol->sol[ia + i] = solutions[(ia - 1) + i];
    }
  }

  return sol;
}

void metric_mshmet(MSHMET_pSol sol, MSHMET_Info *info, const KN< double > &metric) {
  static const int wrapperMetric[6] = {0, 1, 3, 2, 4, 5};
  int k, ia, i;

  cout << " info->iso " << info->iso << endl;
  if (info->iso == 1) {
    cout << " info->iso 11 " << info->iso << endl;
    sol->met = (double *)M_calloc(sol->np + 1, sizeof(double), "sol->met");
    assert(sol->met);

    // isotrope
    for (k = 1; k <= sol->np; k++) {
      sol->met[k] = metric[k - 1];
    }
  } else {
    // anisotropie :: Hessian
    sol->met = (double *)M_calloc(sol->np + 1, 6 * sizeof(double), "sol->met6");
    assert(sol->met);

    for (k = 1; k <= sol->np; k++) {
      ia = (k - 1) * 6 + 1;

      for (i = 0; i < 6; i++) {
        sol->met[ia + i] = metric[(ia - 1) + wrapperMetric[i]];
      }
    }
  }
}

void metric_mshmet_to_ff_metric(MSHMET_pSol sol, MSHMET_Info *info, KN< double > &metric) {
  static const int invwrapperMetric[6] = {0, 1, 3, 2, 4, 5};
  int k, ia, i;

  if (info->iso == 1) {
    cout << " info->iso "
         << " metric " << metric.N( ) << " " << sol->np << endl;

    // isotrope
    for (k = 1; k <= sol->np; k++) {
      metric[k - 1] = sol->met[k];
    }
  } else {
    for (k = 1; k <= sol->np; k++) {
      ia = (k - 1) * 6 + 1;

      for (i = 0; i < 6; i++) {
        metric[(ia - 1) + i] = sol->met[ia + invwrapperMetric[i]];
      }
    }
  }
}

class mshmet3d_Op : public E_F0mps {
 public:
  typedef KN_< double > Result;
  Expression eTh;
  int nbsol;
  int nbsolsize;
  int typesol[GmfMaxTyp];
  int dim;
  vector< Expression > sol;

  static const int n_name_param = 12;    //
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

  int arg(int i, Stack stack, bool a) const {
    return nargs[i] ? GetAny< bool >((*nargs[i])(stack)) : a;
  }

 public:
  mshmet3d_Op(const basicAC_F0 &args) : sol(args.size( ) - 1) {
    // cout << "mshmet3d"<< endl;
    args.SetNameParam(n_name_param, name_param, nargs);
    eTh = to< pmesh3 >(args[0]);
    dim = 3;
    nbsol = args.size( ) - 1;
    int ksol = 0;
    ffassert(nbsol < GmfMaxTyp);

    for (int i = 1; i < nbsol + 1; i++) {
      if (args[i].left( ) == atype< E_Array >( )) {
        const E_Array *a = dynamic_cast< const E_Array * >(args[i].LeftValue( ));
        ffassert(a);
        ksol += a->size( );
      } else {
        ksol++;
      }
    }

    sol.resize(ksol);

    // typesol :: 1 sca, 2 vector, 3 symtensor

    ksol = 0;
    nbsolsize = 0;

    for (int i = 1; i < nbsol + 1; i++) {
      if (args[i].left( ) == atype< E_Array >( )) {
        const E_Array *a = dynamic_cast< const E_Array * >(args[i].LeftValue( ));
        ffassert(a);
        int N = a->size( );
        nbsolsize = nbsolsize + N;

        switch (N) {
          case 3:
            typesol[i - 1] = 2;

            for (int j = 0; j < N; j++) {
              sol[ksol++] = to< double >((*a)[j]);
            }

            break;
          case 6:
            typesol[i - 1] = 3;

            for (int j = 0; j < N; j++) {
              sol[ksol++] = to< double >((*a)[j]);
            }

            break;
          default:
            CompileError(" 3D solution for mshmest is vector (3 comp) or symetric tensor (6 comp)");
            break;
        }
      } else {
        typesol[i - 1] = 1;
        nbsolsize = nbsolsize + 1;
        sol[ksol++] = to< double >(args[i]);
      }
    }
  }

  static ArrayOfaType typeargs( ) { return ArrayOfaType(atype< pmesh3 >( ), true); }    // all type

  static E_F0 *f(const basicAC_F0 &args) { return new mshmet3d_Op(args); }

  AnyType operator( )(Stack stack) const;
  operator aType( ) const { return atype< KN_< double > >( ); }
};

basicAC_F0::name_and_type mshmet3d_Op::name_param[] = {
  {"loptions", &typeid(KN_< long >)},    // 0
  {"doptions", &typeid(KN_< double >)},
  {"metric", &typeid(KN_< double >)},
  {"normalization", &typeid(bool)},
  {"aniso", &typeid(bool)},
  {"levelset", &typeid(bool)},    // 5
  {"verbosity", &typeid(long)},
  {"nbregul", &typeid(long)},
  {"hmin", &typeid(double)},
  {"hmax", &typeid(double)},    // 9
  {"err", &typeid(double)},     // 10
  {"width", &typeid(double)}    // 11
};

template< class T >
ostream &dumpp(const T *p, int n, ostream &f) {
  for (int i = 0; i < n; ++i) {
    f << p[i] << " ";
  }

  return f;
}

AnyType mshmet3d_Op::operator( )(Stack stack) const {
  // initialisation
  MeshPoint *mp(MeshPointStack(stack)), mps = *mp;
  Mesh3 *pTh = GetAny< Mesh3 * >((*eTh)(stack));

  ffassert(pTh);
  Mesh3 &Th3 = *pTh;
  int nv = Th3.nv;
  int nt = Th3.nt;
  int nbe = Th3.nbe;

  KN< double > defaultfopt(4);
  /*
   * info->hmin   = fopt[0]; // 0.01;
   * info->hmax   = fopt[1]; // 1.0;
   * info->eps    = fopt[2]; // 0.01;
   * info->width  = fopt[3]; // 0.05;
   */
  defaultfopt(0) = 0.01;
  defaultfopt(1) = 1.0;
  defaultfopt(2) = 0.01;
  defaultfopt(3) = 0.00;
  KN< long > defaultintopt(7);
  /*
   * info->nnu    = intopt[0]; //  0;
   * info->iso    = intopt[1]; //  1;
   * info->ls     = intopt[2]; //  0;
   * info->ddebug = intopt[3]; //  0;
   * info->imprim = intopt[4]; // 10;
   * info->nlis   = intopt[5]; //  0;
   * info->metric =  intopt[6]; // 0; // metric given besoin ???
   */
  defaultintopt(0) = 1;
  defaultintopt(1) = 1;
  defaultintopt(2) = 0;
  defaultintopt(3) = 0;
  defaultintopt(4) = verbosity;
  defaultintopt(5) = 0;
  defaultintopt(6) = 0;
  if (nargs[11]) {
    defaultintopt[2] = 1;
    defaultfopt[3] = 0;
  }    // level set ...

  KN< int > intopt(arg(0, stack, defaultintopt));
  KN< double > fopt(arg(1, stack, defaultfopt));
  KN< double > *pmetric = new KN< double >(Th3.nv);
  KN< double > &metric = *pmetric;

  intopt[0] = arg(3, stack, intopt[0] != 0);     // normaliz
  intopt[1] = !arg(4, stack, intopt[1] == 0);    // aniso
  intopt[2] = arg(5, stack, intopt[2] != 0);     // levelset
  intopt[4] = arg(6, stack, intopt[4]);          // verbo
  intopt[5] = arg(7, stack, intopt[5]);          // nbregul
  fopt[0] = arg(8, stack, fopt[0]);              // hmin
  fopt[1] = arg(9, stack, fopt[1]);              // hmax
  fopt[2] = arg(10, stack, fopt[2]);             // err
  fopt[3] = arg(11, stack, fopt[3]);             // width

  if (verbosity > 2) {
    cout << "    -- mshmet : lopt ";
    dumpp((int *)intopt, intopt.N( ), cout) << endl;
    cout << "              : dopt ";
    dumpp((double *)fopt, fopt.N( ), cout) << endl;
  }

  metric = 0.;
  if (intopt.N( ) != 7) {
    ExecError(" Size of array of loption are wrong != 7");
  }

  if (fopt.N( ) != 4) {
    ExecError(" Size of array of doption are wrong != 4");
  }

  if (intopt[1] == 0) {
    metric.resize(6 * Th3.nv);
    metric = 0.;
  }

  // definiton d'une metric par default
  if (nargs[2]) {
    metric = GetAny< KN_< double > >((*nargs[2])(stack));
    assert(metric.N( ) == Th3.nv || metric.N( ) == 6 * Th3.nv);
    intopt[6] = 1;
    if (metric.N( ) == Th3.nv) {
      intopt[1] = 1;
    }

    if (metric.N( ) == 6 * Th3.nv) {
      intopt[1] = 0;
    }
  }

  MSHMET_pMesh mshmetmesh = mesh3_to_MSHMET_pMesh(Th3);
  int TypTab[nbsol];

  for (int ii = 0; ii < nbsol; ii++) {
    TypTab[ii] = typesol[ii];
  }

  KN< double > tabsol(nbsolsize * nv);
  tabsol = 0.;
  {
    MeshPoint *mp3(MeshPointStack(stack));

    KN< bool > takemesh(nv);
    takemesh = false;

    for (int it = 0; it < nt; it++) {
      for (int iv = 0; iv < 4; iv++) {
        int i = Th3(it, iv);

        if (takemesh[i] == false) {
          mp3->setP(&Th3, it, iv);

          for (int ii = 0; ii < nbsolsize; ii++) {
            tabsol[i * nbsolsize + ii] = GetAny< double >((*sol[ii])(stack));
          }

          takemesh[i] = true;
        }
      }
    }
  }
  if (verbosity > 5) {
    cout << "    min/max tabsol:  " << tabsol.min( ) << " " << tabsol.max( ) << endl;
  }

  MSHMET_pSol mshmetsol = sol_mshmet(dim, nv, nbsol, nbsolsize, TypTab, tabsol);
  if (intopt[1] == 1) {
    mshmetmesh->info.iso = 1;
  } else {
    mshmetmesh->info.iso = 0;
  }

  if (nargs[2]) {
    metric_mshmet(mshmetsol, &mshmetmesh->info, metric);
  }

  int res = MSHMET_mshmet(intopt, fopt, mshmetmesh, mshmetsol);

  if (res > 0) {
    cout << " problem with mshmet :: error " << res << endl;
    exit(1);
  }

  metric_mshmet_to_ff_metric(mshmetsol, &mshmetmesh->info, metric);

  // faire les free

  myMSHMET_free(mshmetmesh, mshmetsol);
  if (verbosity > 1000) {
    M_memDump( );
  }

  Add2StackOfPtr2Free(stack, pmetric);
  *mp = mps;
  return SetAny< KN< double > >(metric);
}

// mshmet2d
class mshmet2d_Op : public E_F0mps {
 public:
  typedef KN_< double > Result;
  Expression eTh;
  int nbsol;
  int nbsolsize;
  int typesol[GmfMaxTyp];
  int dim;
  vector< Expression > sol;

  static const int n_name_param = 12;    //
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
  mshmet2d_Op(const basicAC_F0 &args) : sol(args.size( ) - 1) {
    cout << "mshmet2d" << endl;
    args.SetNameParam(n_name_param, name_param, nargs);
    eTh = to< pmesh >(args[0]);
    dim = 3;
    nbsol = args.size( ) - 1;
    int ksol = 0;
    ffassert(nbsol < GmfMaxTyp);

    for (int i = 1; i < nbsol + 1; i++) {
      if (args[i].left( ) == atype< E_Array >( )) {
        const E_Array *a = dynamic_cast< const E_Array * >(args[i].LeftValue( ));
        ffassert(a);
        ksol += a->size( );
      } else {
        ksol++;
      }
    }

    sol.resize(ksol);

    // typesol :: 1 sca, 2 vector, 3 symtensor

    ksol = 0;
    nbsolsize = 0;

    for (int i = 1; i < nbsol + 1; i++) {
      if (args[i].left( ) == atype< E_Array >( )) {
        const E_Array *a = dynamic_cast< const E_Array * >(args[i].LeftValue( ));
        ffassert(a);
        int N = a->size( );
        nbsolsize = nbsolsize + N;

        switch (N) {
          case 2:
            typesol[i - 1] = 2;

            for (int j = 0; j < N; j++) {
              sol[ksol++] = to< double >((*a)[j]);
            }

            break;
          case 3:
            typesol[i - 1] = 3;

            for (int j = 0; j < N; j++) {
              sol[ksol++] = to< double >((*a)[j]);
            }

            break;
          default:
            CompileError(" 2D solution for mshmest is vector (2 comp) or symetric tensor (3 comp)");
            break;
        }
      } else {
        typesol[i - 1] = 1;
        nbsolsize = nbsolsize + 1;
        sol[ksol++] = to< double >(args[i]);
      }
    }
  }

  static ArrayOfaType typeargs( ) { return ArrayOfaType(atype< pmesh >( ), true); }    // all type

  static E_F0 *f(const basicAC_F0 &args) { return new mshmet2d_Op(args); }

  AnyType operator( )(Stack stack) const;
  operator aType( ) const { return atype< KN_< double > >( ); }
};

basicAC_F0::name_and_type mshmet2d_Op::name_param[] = {{"loptions", &typeid(KN_< long >)},
                                                       {"doptions", &typeid(KN_< double >)},
                                                       {"metric", &typeid(KN_< double >)}};
AnyType mshmet2d_Op::operator( )(Stack stack) const {
  // initialisation
  MeshPoint *mp(MeshPointStack(stack)), mps = *mp;
  const Mesh *pTh = GetAny< const Mesh * >((*eTh)(stack));

  ffassert(pTh);
  const Mesh &Th = *pTh;
  int nv = Th.nv;
  int nt = Th.nt;
  int nbe = Th.neb;

  KN< double > defaultfopt(4);
  /*
   * info->hmin   = fopt[0]; // 0.01;
   * info->hmax   = fopt[1]; // 1.0;
   * info->eps    = fopt[2]; // 0.01;
   * info->width  = fopt[3]; // 0.05;
   */
  defaultfopt(0) = 0.01;
  defaultfopt(1) = 1.0;
  defaultfopt(2) = 0.01;
  defaultfopt(3) = 0.05;
  KN< long > defaultintopt(7);
  /*
   * info->nnu    = intopt[0]; //  0;
   * info->iso    = intopt[1]; //  1;
   * info->ls     = intopt[2]; //  0;
   * info->ddebug = intopt[3]; //  0;
   * info->imprim = intopt[4]; // 10;
   * info->nlis   = intopt[5]; //  0;
   * info->metric =  intopt[6]; // 0; // metric given besoin ???
   */
  defaultintopt(0) = 0;
  defaultintopt(1) = 1;
  defaultintopt(2) = 0;
  defaultintopt(3) = 1;
  defaultintopt(4) = 10;
  defaultintopt(5) = 0;
  defaultintopt(6) = 0;

  KN< int > intopt(arg(0, stack, defaultintopt));
  KN< double > fopt(arg(1, stack, defaultfopt));

  KN< double > *pmetric = new KN< double >(Th.nv);
  KN< double > &metric = *pmetric;

  if (intopt[1] == 1) {
    metric.resize(6 * Th.nv);
    metric = 0.;
  }

  // definiton d'une metric par default
  if (nargs[2]) {
    metric = GetAny< KN_< double > >((*nargs[2])(stack));
    assert(metric.N( ) == Th.nv || metric.N( ) == 6 * Th.nv);
    intopt[6] = 1;
    if (metric.N( ) == Th.nv) {
      intopt[1] = 1;
    }

    if (metric.N( ) == 6 * Th.nv) {
      intopt[1] = 0;
    }
  }

  MSHMET_pMesh mshmetmesh = mesh_to_MSHMET_pMesh(Th);
  int TypTab[nbsol];

  for (int ii = 0; ii < nbsol; ii++) {
    TypTab[ii] = typesol[ii];
  }

  KN< double > tabsol(nbsolsize * nv);
  tabsol = 0.;
  {
    MeshPoint *mp3(MeshPointStack(stack));

    KN< bool > takemesh(nv);
    takemesh = false;

    for (int it = 0; it < nt; it++) {
      for (int iv = 0; iv < 3; iv++) {
        int i = Th(it, iv);

        if (takemesh[i] == false) {
          mp3->setP(&Th, it, iv);

          for (int ii = 0; ii < nbsolsize; ii++) {
            tabsol[i * nbsolsize + ii] = GetAny< double >((*sol[ii])(stack));
          }

          takemesh[i] = true;
        }
      }
    }
  }
  MSHMET_pSol mshmetsol = sol_mshmet(dim, nv, nbsol, nbsolsize, TypTab, tabsol);
  if (intopt[1] == 1) {
    mshmetmesh->info.iso = 1;
  } else {
    mshmetmesh->info.iso = 0;
  }

  if (nargs[2]) {
    metric_mshmet(mshmetsol, &mshmetmesh->info, metric);
  }

  int res = MSHMET_mshmet(intopt, fopt, mshmetmesh, mshmetsol);

  if (res > 0) {
    cout << " problem with mshmet :: error " << res << endl;
    exit(1);
  }

  metric_mshmet_to_ff_metric(mshmetsol, &mshmetmesh->info, metric);

  // faire les free
  myMSHMET_free(mshmetmesh, mshmetsol);
  *mp = mps;
  if (verbosity > 1000) {
    M_memDump( );
  }

  Add2StackOfPtr2Free(stack, pmetric);
  return SetAny< KN< double > >(metric);
}

static void Load_Init( ) {    // le constructeur qui ajoute la fonction "splitmesh3"  a freefem++
  if (verbosity) {
    cout << " load: mshmet  " << endl;
  }

  Global.Add("mshmet", "(", new OneOperatorCode< mshmet2d_Op >);
  Global.Add("mshmet", "(", new OneOperatorCode< mshmet3d_Op >);
}

LOADFUNC(Load_Init)
