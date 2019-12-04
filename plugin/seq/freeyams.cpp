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
// SUMMARY : Freeyams - FreeFem++ link
// LICENSE : LGPLv3
// ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
// AUTHORS : Jacques Morice
// E-MAIL  : jacques.morice@ann.jussieu.fr

/* clang-format off */
//ff-c++-LIBRARY-dep: freeyams libMesh
//ff-c++-cpp-dep:
/* clang-format on */

/*
 * Thank to the ARN ()  FF2A3 grant
 * ref:ANR-07-CIS7-002-01
 */

/*
 * ff-c++ -auto freeyams.cpp
 */

// ./ff-c++ yams.cpp -I../src/libMesh/ -I../3rdparty/include/yams/ -L../3rdparty/lib/yams/ -lyams2
// -L/Users/morice/work/postdoc/freefem++prod/src/libMesh/ -lMesh
#define _USE_MATH_DEFINES
#include "ff++.hpp"
#include "msh3.hpp"
#include "memory.h"
#include "freeyamslib.h"
#include "eigenv.h"    // include dans libMesh

using namespace Fem2D;
using namespace yams;

// 3d mesh function
// mesh3 -> yams_pSurfMesh
void mesh3_to_yams_pSurfMesh(const Mesh3 &Th, int memory, int choix, yams_pSurfMesh meshyams) {
  int k;
  int npinit, neinit;

  meshyams->dim = 3;
  meshyams->npfixe = Th.nv;
  meshyams->nefixe = Th.nbe;
  meshyams->ntet = Th.nt;
  meshyams->nafixe = 0;    // Edges
  meshyams->nvfixe = 0;    // Normals
  meshyams->ntfixe = 0;    // Tangents
  npinit = meshyams->npfixe;
  neinit = meshyams->nefixe;
  // cette fonction change la taille des tableaux en fonctions des options : choix, memory, sm->type
  zaldy1(meshyams->nefixe, meshyams->npfixe, meshyams->nvfixe, memory, meshyams, choix);

  yams_pPoint ppt;

  for (k = 1; k <= npinit; k++) {
    ppt = &meshyams->point[k];
    ppt->c[0] = Th.vertices[k - 1].x;
    ppt->c[1] = Th.vertices[k - 1].y;
    ppt->c[2] = Th.vertices[k - 1].z;
    ppt->ref = Th.vertices[k - 1].lab & 0x7fff;

    ppt->tag = M_UNUSED;
    ppt->color = 0;
    ppt->size = -1.;
    ppt->tge = 0;
    ppt->geom = M_CURVE;
  }

  meshyams->npfixe = npinit;

  /* read mesh triangles */
  yams_pTriangle ptriangle;

  for (k = 1; k <= neinit; k++) {
    const Triangle3 &K(Th.be(k - 1));
    ptriangle = &meshyams->tria[k];
    ptriangle->v[0] = Th.operator( )(K[0]) + 1;
    ptriangle->v[1] = Th.operator( )(K[1]) + 1;
    ptriangle->v[2] = Th.operator( )(K[2]) + 1;
    ptriangle->ref = K.lab & 0x7fff;
  }

  /* tetrahedra */
  if (meshyams->ntet) {
    yams_pTetra ptetra;
    meshyams->tetra = (yams_Tetra *)calloc((meshyams->ntet + 1), sizeof(yams_Tetra));
    assert(meshyams->tetra);

    for (k = 1; k <= meshyams->ntet; k++) {
      const Tet &K(Th.elements[k - 1]);
      ptetra = &meshyams->tetra[k];
      ptetra->v[0] = Th.operator( )(K[0]) + 1;
      ptetra->v[1] = Th.operator( )(K[1]) + 1;
      ptetra->v[2] = Th.operator( )(K[2]) + 1;
      ptetra->v[3] = Th.operator( )(K[3]) + 1;
      ptetra->ref = K.lab & 0x7fff;
    }
  }

  meshyams->ne = meshyams->nefixe;
  meshyams->np = meshyams->npfixe;
}

// meshS -> yams_pSurfMesh
void meshS_to_yams_pSurfMesh(const MeshS &Th, int memory, int choix, yams_pSurfMesh meshyams) {

  int k;
  int npinit, neinit;

  meshyams->dim = 3;
  meshyams->npfixe = Th.nv;
  meshyams->nefixe = Th.nt;    // Th.nbe;
  meshyams->ntet = 0;          // Th.nt;
  meshyams->nafixe = 0;        // Edges
  meshyams->nvfixe = 0;        // Normals
  meshyams->ntfixe = 0;        // Tangents
  npinit = meshyams->npfixe;
  neinit = meshyams->nefixe;
  // cette fonction change la taille des tableaux en fonctions des options : choix, memory, sm->type
  zaldy1(meshyams->nefixe, meshyams->npfixe, meshyams->nvfixe, memory, meshyams, choix);

  yams_pPoint ppt;

  for (k = 1; k <= npinit; k++) {
    ppt = &meshyams->point[k];
    ppt->c[0] = Th.vertices[k - 1].x;
    ppt->c[1] = Th.vertices[k - 1].y;
    ppt->c[2] = Th.vertices[k - 1].z;
    ppt->ref = Th.vertices[k - 1].lab & 0x7fff;

    ppt->tag = M_UNUSED;
    ppt->color = 0;
    ppt->size = -1.;
    ppt->tge = 0;
    ppt->geom = M_CURVE;
  }

  meshyams->npfixe = npinit;

  /* read mesh triangles */
  yams_pTriangle ptriangle;

  for (k = 1; k <= neinit; k++) {
    const TriangleS &K(Th.elements[k - 1]);
    ptriangle = &meshyams->tria[k];
    ptriangle->v[0] = Th.operator( )(K[0]) + 1;
    ptriangle->v[1] = Th.operator( )(K[1]) + 1;
    ptriangle->v[2] = Th.operator( )(K[2]) + 1;
    ptriangle->ref = K.lab & 0x7fff;
  }

  meshyams->ne = meshyams->nefixe;
  meshyams->np = meshyams->npfixe;
}

MeshS *yams_pSurfMesh_to_meshS(yams_pSurfMesh sm, int infondang, int infocc, int choix) {
  /*
   * Mesh3  :: maillage initiale
   * memory :: memoire pour yams
   * choix  :: option du remaillage
   * ref    ::
   */
  // variable a enlever par la suite
  yams_pGeomSupp gs;
  yams_pGeomtge gt;
  yams_pPoint ppt;
  yams_pTriangle pt1;
  yams_pTetra ptt;
  yams_pEdge pte;
  int i, k, np, ne, nn, nt, nav, natv, tatv, nbl;
  int nedge, nridge, ndang, nrequis;
  int is1, is2, ncorner, prequis;

  // freefempp variable
  int ff_nv, ff_nt, ff_nbe;
  ff_nbe = 0;
  /* mark connected component */
  ne = 0;

  for (k = 1; k <= sm->npmax; k++) {
    ppt = &sm->point[k];
    ppt->tag |= M_UNUSED;
    ppt->flag = ppt->color = 0;
  }

  // a enlever pour l'instant
  if (sm->connex > 0) {
    for (k = 1; k <= sm->ne; k++) {
      pt1 = &sm->tria[k];
      if (pt1->v[0] > 0 && pt1->cc == sm->connex) {
        ne++;

        for (i = 0; i < 3; i++) {
          ppt = &sm->point[pt1->v[i]];
          ppt->tag &= ~M_UNUSED;
        }
      }
    }
  } else {
    /* mark used faces */
    for (k = 1; k <= sm->ne; k++) {
      pt1 = &sm->tria[k];
      if (!pt1->v[0]) {
        continue;
      }

      ++ne;

      for (i = 0; i < 3; i++) {
        ppt = &sm->point[pt1->v[i]];
        ppt->tag &= ~M_UNUSED;
      }
    }
  }

  cout << "sm->ntet=" << sm->ntet << endl;

  /* mark used vertices */
  np = nav = 0;
  ncorner = prequis = 0;

  for (k = 1; k <= sm->npmax; k++) {
    ppt = &sm->point[k];
    if (ppt->tag & M_UNUSED) {
      continue;
    }

    ppt->tmp = ++np;
    if (ppt->tag == M_NOTAG) {
      nav++;
    }
  }

  ff_nv = np;    // number of vertex
  //
  Vertex3 *ff_v = new Vertex3[ff_nv];
  int kk = 0;

  for (k = 1; k <= sm->npmax; k++) {
    ppt = &sm->point[k];
    if (ppt->tag & M_UNUSED) {
      continue;
    }

    ff_v[kk].x = ppt->c[0];
    ff_v[kk].y = ppt->c[1];
    ff_v[kk].z = ppt->c[2];
    ff_v[kk].lab = ppt->ref;
    kk++;
    if (ppt->tag & M_CORNER) {
      ncorner++;
    }

    if (ppt->tag & M_REQUIRED) {
      prequis++;
    }
  }

  assert(kk == ff_nv);
  // write triangle
  nedge = sm->dim == 3 ? infondang : 0;
  nridge = nrequis = nn = nt = natv = tatv = 0;

  for (k = 1; k <= sm->ne; k++) {
    pt1 = &sm->tria[k];
    if (!pt1->v[0]) {
      continue;
    } else if (sm->connex > 0 && pt1->cc != sm->connex) {
      continue;
    }

    nt++;
  }

  ff_nt = nt;
  TriangleS *ff_t = new TriangleS[ff_nt];
  TriangleS *ff_tt = ff_t;

  for (k = 1; k <= sm->ne; k++) {
    int iv[3], lab;
    pt1 = &sm->tria[k];
    // lab = pt1->ref;
    if (!pt1->v[0]) {
      continue;
    } else if (sm->connex > 0 && pt1->cc != sm->connex) {
      continue;
    }

    iv[0] = sm->point[pt1->v[0]].tmp - 1;
    iv[1] = sm->point[pt1->v[1]].tmp - 1;
    iv[2] = sm->point[pt1->v[2]].tmp - 1;
    lab = pt1->ref;    // change fh 02/2013
    (*ff_tt++).set(ff_v, iv, lab);

    for (i = 0; i < 3; i++) {
      ppt = &sm->point[pt1->v[i]];
      gs = &sm->geom[pt1->vn[i]];
      gt = &sm->tgte[ppt->tge];
      if (ppt->tag > M_NOTAG) {
        natv++;
        if (ppt->tag & M_CORNER) {
          tatv++;
        }
      }

      if (!gs->newnum) {
        gs->newnum = ++nn;
      }

      if (!gt->newnum) {
        gt->newnum = ++nt;
      }

      if (!pt1->edg[i] && pt1->tag[i] == M_NOTAG) {
        continue;
      } else if (pt1->adj[i] && (k > pt1->adj[i])) {
        continue;
      }

      nedge++;
      if (pt1->tag[i] & M_RIDGE_GEO) {
        nridge++;
      }

      if (pt1->tag[i] & M_REQUIRED) {
        nrequis++;
      }
    }
  }

  BoundaryEdgeS *ff_b = new BoundaryEdgeS[ff_nbe];
  BoundaryEdgeS *ff_bb = ff_b;

  // les autres avoir par la suite
  if (verbosity > 1) {
    cout << " nv " << ff_nv << " nt" << ff_nt << " nbe" << ff_nbe << endl;
  }

  MeshS *THS_T = new MeshS(ff_nv, ff_nt, ff_nbe, ff_v, ff_t, ff_b);

  return THS_T;
}

// TODO CHECK
void solyams_pSurfMesh(yams_pSurfMesh sm, const int &type, const KN< double > &tabsol, float hmin,
                       float hmax) {
  yams_pPoint ppt;
  yams_pMetric pm;
  int i, k;
  double sizeh, m[6], lambda[3], vp[2][2], vp3[3][3];

  hmin = FLT_MAX;
  hmax = -FLT_MAX;
  float vpmin = FLT_MAX, vpmax = -FLT_MAX, mmin = FLT_MAX, mmax = -FLT_MAX;

  if (type == 1) {
    for (k = 1; k <= sm->npfixe; k++) {
      ppt = &sm->point[k];
      ppt->size = (float)tabsol[k - 1];    // change FH nov 2010: k -> k-1
      hmin = min(ppt->size, hmin);
      hmax = max(ppt->size, hmax);
    }
  } else if (type == 3) {
    if (!sm->metric && !zaldy3(sm, 3)) {
      ExecError("Pb alloc metric in freeyam ??? ");
    }

    for (k = 1; k <= sm->npfixe; k++) {
      ppt = &sm->point[k];
      pm = &sm->metric[k];    // coorrection FH dec 2010..
      memset(pm->m, 6 * sizeof(float), 0.);

      for (i = 0; i < 6; i++) {
        m[i] = (float)tabsol[(k - 1) * 6 + i];
      }

      pm->m[0] = m[0];
      pm->m[1] = m[1];
      pm->m[2] = m[3];
      pm->m[3] = m[2];
      pm->m[4] = m[4];
      pm->m[5] = m[5];
      pm->k1 = pm->k2 = (float)FLT_MAX;

      for (i = 0; i < 6; i++) {
        m[i] = pm->m[i];
      }

      if (!eigenv(1, m, lambda, vp3)) {
        fprintf(stderr, "  ## ERR 9201, inbbf, Not a metric tensor. Discarded\n");
        free(sm->metric);
        sm->metric = 0;
        ExecError("freeyamerr: ## ERR 9201, inbbf, Not a metric tensor. Discarded");
      }

      float vmn = min(min(lambda[0], lambda[1]), lambda[2]);
      float vmx = max(max(lambda[0], lambda[1]), lambda[2]);

      vpmin = min(vpmin, vmn);
      vpmax = max(vpmax, vmx);
      sizeh = vpmax;
      ppt->size = max(1.0 / sqrt(sizeh), EPS);
      hmin = min(ppt->size, hmin);
      hmax = max(ppt->size, hmax);
    }
  }

  // if(verbosity>4)
  {
    cout << " freeyams (metric in) :  hmin " << hmin << " , hmax " << hmax << endl;
    if (type == 3) {
      cout << "             min max of eigen val  " << vpmin << " " << vpmax << endl;
    }
  }

  if (type == 3 && vpmin < 0) {
    cout << "   Error Freeyam :  metric    min max of eigen val  " << vpmin << " " << vpmax << endl;
    ExecError("Error in metric definition freeyams (negative eigen value");
  }
}

static const int wrapper_intopt[13] = {0, 3, 7, 8, 9, 11, 12, 13, 14, 15, 17, 18, 22};

/*
 * static const int wrapper_fopt[12] = {  0, 1, 3,  4,  6,
 *                                     7, 8, 9, 10, 11,
 *                                    12, 13};
 */
static const int wrapper_fopt[11] = {1, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13};

void yams_inival(int intopt[23], double fopt[14]) {
  /*
   *  intopt : 0  !! anisotropie
   *           1  !! ecp  // enl
   *           2  !! extended out put file // enl
   *           3  !! FE correction
   *           4  !! Formatted (ascii) output file // enl
   *           5  !! save metric file // enl
   *           6  !! msh2  // enl
   *           7  !! Split multiple connected points
   *           8  !! memory
   *           9  !! connected component
   *          10  !! vrml  //enl
   *          11  !! imprim
   *          12  !! nm : Create point on straight edge (no mapping)
   *          13  !! nc : No validity check during smoothing (opt. 9)
   *          14  !! np : Specify number of points desired
   *          15  !! nit : Nb Iter
   *          16  !! nq  : Output quads // enl
   *          17  !! nr  : No ridge detection
   *          18  !! ns  : No point smoothing
   *          19  !! no  : No output file  // enl
   *          20  !! ref : Ignore face references // enl
   *          // rajouter lors de l'ouverture du fichiers yams
   *          21  !! absolute : opts.ctrl &= ~REL; par default 1 // enl
   *          22  !! set optim option
   *
   *  fopt   : 0  !! iso
   *           1  !! eps
   *           pas de 2
   *           3  !! opts.lambda
   *           4  !! opts.mu
   *           pas de 5
   *           6  !! hgrad  :: opts.shock
   *           7  !! hmin   :: opts.hmin
   *           8  !! hmax   :: opts.hmax
   *           // rajouter lors de l'ouverture du fichiers yams
   *           9  !! tolerance :: opts.bande
   *           10 !! degrad :: opts.degrad
   *           11 !! declic :: opts.declic
   *           12 !! walton :: opts.walton = cos(dummy/180.0*M_PI);
   *           13 !! ridge  :: opts.ridge
   */

  /* Set default values for options */
  // fopt 5,
  fopt[7] = -2.0;
  fopt[8] = -2.0;
  fopt[6] = 1.3;  /* default mesh gradation     */
  fopt[1] = 0.01; /* geometric approximation    */
  fopt[0] = 0.0;
  fopt[11] = 1.0 / BETAC;
  fopt[3] = -1.0;
  fopt[4] = -1.0;
  fopt[13] = 45.;      // default RIDG = 45.
  fopt[12] = COS45DEG; /* Walton limitation          */
  fopt[9] = -2;        /* default = 1 unit           */
  fopt[10] = QUALCOE;  /* quality degradation        */
  // opts.ctrl   =   REL | ISO;  initialisation by default

  // intopt :: 3,7,13,14,15,20
  intopt[15] = -1;
  intopt[13] = 0;
  intopt[14] = -1;

  /* get decimation parameters */
  intopt[20] = 0;
  intopt[3] = 0;
  intopt[7] = 0;     // Split multiple connected points  (no manifold)
  intopt[22] = 1;    // set optim option

  // demander P. Frey
  intopt[0] = 0;    // anisotropie
  intopt[1] = 0;    //
  intopt[2] = 0;

  intopt[4] = 0;
  intopt[5] = 0;
  intopt[6] = 0;

  intopt[8] = -1;    // memory
  intopt[9] = -1;    // par default   connex connected component (tout)
  intopt[10] = 0;    // vrml
  intopt[11] = verbosity;
  intopt[12] = 0;    // nm

  intopt[16] = 0;    // quad
  intopt[17] = 0;    // noridge
  intopt[18] = 0;    // nosmooth
  intopt[19] = 1;    // 1
  intopt[21] = 1;
}

// version with meshS in arg

class yams_Op_meshS : public E_F0mps {
 public:
  typedef pmeshS Result;
  Expression eTh;
  int nbsol;
  int nbsolsize;
  int type;
  int dim;
  vector< Expression > sol;

  static const int n_name_param = 14;    //
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

  int arg(int i, Stack stack, int a) const {
    return nargs[i] ? GetAny< long >((*nargs[i])(stack)) : a;
  }

  bool arg(int i, Stack stack, bool a) const {
    return nargs[i] ? GetAny< bool >((*nargs[i])(stack)) : a;
  }

 public:
  yams_Op_meshS(const basicAC_F0 &args) : sol(args.size( ) - 1) {
    cout << "yams" << endl;
    args.SetNameParam(n_name_param, name_param, nargs);
    eTh = to< pmeshS >(args[0]);
    dim = 3;
    nbsol = args.size( ) - 1;
    if (nbsol > 1) {
      CompileError(" yams accept only one solution ");
    }

    int ksol = 0;

    if (nbsol == 1) {
      int i = 1;
      if (args[i].left( ) == atype< E_Array >( )) {
        const E_Array *a = dynamic_cast< const E_Array * >(args[i].LeftValue( ));
        ffassert(a);
        ksol += a->size( );
      } else {
        ksol++;
      }

      sol.resize(ksol);

      // type :: 1 sca, 2 vector, 3 symtensor

      ksol = 0;
      nbsolsize = 0;
      type = 0;

      if (args[i].left( ) == atype< E_Array >( )) {
        const E_Array *a = dynamic_cast< const E_Array * >(args[i].LeftValue( ));
        ffassert(a);
        int N = a->size( );
        nbsolsize = nbsolsize + N;

        switch (N) {
            /*
             * case 3 :
             *  type[i-1]=2;
             *  for (int j=0;j<N;j++)
             *  sol[ksol++]=to<double>((*a)[j]);
             *  break;
             */
          case 6:
            type = 3;

            for (int j = 0; j < N; j++) {
              sol[ksol++] = to< double >((*a)[j]);
            }

            break;
          default:
            CompileError(
              " 3D solution for yams is a scalar (1 comp) or a symetric tensor (6 comp)");
            break;
        }
      } else {
        type = 1;
        nbsolsize = nbsolsize + 1;
        sol[ksol++] = to< double >(args[i]);
      }

      if (nargs[2]) {
        CompileError(" we give two metric for yams ");
      }
    }
  }

  static ArrayOfaType typeargs( ) { return ArrayOfaType(atype< pmeshS >( ), true); }    // all type

  static E_F0 *f(const basicAC_F0 &args) { return new yams_Op_meshS(args); }

  AnyType operator( )(Stack stack) const;
  operator aType( ) const { return atype< pmeshS >( ); }
};

basicAC_F0::name_and_type yams_Op_meshS::name_param[] = {
  {"loptions", &typeid(KN_< long >)},    // 0
  {"doptions", &typeid(KN_< double >)},
  {"metric", &typeid(KN_< double >)},
  {"aniso", &typeid(bool)},    // 3
  {"mem", &typeid(long)},
  {"hmin", &typeid(double)},
  {"hmax", &typeid(double)},    // 6
  {"gradation", &typeid(double)},
  {"option", &typeid(long)},          // 8
  {"ridgeangle", &typeid(double)},    // 9
  {"absolute", &typeid(bool)},        // 10
  {"verbosity", &typeid(long)},       // 11

  {"nr", &typeid(long)},    // 12 no ridge
  {"ns", &typeid(long)}     // 13 no point smoothing
};

AnyType yams_Op_meshS::operator( )(Stack stack) const {
  // initialisation
  MeshPoint *mp(MeshPointStack(stack)), mps = *mp;
  MeshS *pTh = GetAny< MeshS * >((*eTh)(stack));

  ffassert(pTh);
  MeshS &ThS = *pTh;
  int nv = ThS.nv;
  int nt = ThS.nt;
  int nbe = ThS.nbe;

  KN< int > defaultintopt(23);
  KN< double > defaultfopt(14);
  defaultintopt = 0;
  defaultfopt = 0.;
  yams_inival(defaultintopt, defaultfopt);

  KN< int > intopt(23);

  for (int ii = 0; ii < 23; ii++) {
    intopt[ii] = defaultintopt[ii];
  }

  KN< double > fopt(14);

  for (int ii = 0; ii < 14; ii++) {
    fopt[ii] = defaultfopt[ii];
  }

  assert(fopt.N( ) == 14);

  if (nargs[0]) {
    KN< int > intopttmp = GetAny< KN_< long > >((*nargs[0])(stack));
    if (intopttmp.N( ) != 13) {
      cerr << "the size of vector loptions is 13 " << endl;
      exit(1);
    } else {
      for (int ii = 0; ii < 13; ii++) {
        intopt[wrapper_intopt[ii]] = intopttmp[ii];
      }
    }
  }

  if (nargs[1]) {
    KN< double > fopttmp = GetAny< KN_< double > >((*nargs[1])(stack));
    if (fopttmp.N( ) != 11) {
      cerr << "the size of vector loptions is 11 not " << fopttmp.N( ) << endl;
      ExecError("FreeYams");
    } else {
      for (int ii = 0; ii < 11; ii++) {
        fopt[wrapper_fopt[ii]] = fopttmp[ii];
      }
    }
  }

  intopt[0] = arg(3, stack, intopt[0] != 1);
  intopt[8] = arg(4, stack, intopt[8]);
  fopt[7] = arg(5, stack, fopt[7]);
  fopt[8] = arg(6, stack, fopt[7]);
  fopt[6] = arg(7, stack, fopt[6]);
  intopt[22] = arg(8, stack, intopt[22]);    // optim option
  if (nargs[9]) {
    intopt[17] = 1;
  }

  fopt[13] = arg(9, stack, fopt[13]);             // ridge angle
  intopt[21] = arg(10, stack, intopt[21]);        // absolue
  intopt[11] = arg(11, stack, (int)verbosity);    // verbosity
  intopt[17] = arg(12, stack, intopt[17]);        // no ridge
  intopt[18] = arg(13, stack, intopt[18]);        // nb smooth
  if (verbosity > 1) {
    cout << " fopt = [";

    for (int i = 0; i < 11; ++i) {
      cout << fopt[wrapper_fopt[i]] << (i < 10 ? "," : "];\n");
    }

    cout << " intopt = [";

    for (int i = 0; i < 13; ++i) {
      cout << intopt[wrapper_intopt[i]] << (i < 12 ? "," : "];\n");
    }
  }

  /*
   * KN<int> intopt(arg(0,stack,defaultintopt));
   * assert( intopt.N() == 23 );
   * KN<double> fopt(arg(1,stack,defaultfopt));
   * assert( fopt.N() == 14 );
   */
  KN< double > metric;

  int mtype = type;
  if (nargs[2]) {
    metric = GetAny< KN_< double > >((*nargs[2])(stack));
    if (metric.N( ) == ThS.nv) {
      mtype = 1;
      intopt[1] = 0;
    } else if (metric.N( ) == 6 * ThS.nv) {
      intopt[1] = 1;
      mtype = 3;
    } else {
      cerr << "sizeof vector metric is incorrect, size will be Th.nv or 6*Th.nv" << endl;
    }
  } else if (nbsol > 0) {
    if (type == 1) {
      intopt[1] = 0;
      metric.resize(ThS.nv);
      metric = 0.;
    } else if (type == 3) {
      intopt[1] = 1;
      metric.resize(6 * ThS.nv);
      metric = 0.;
    }
  } else {
    if (intopt[1] == 0) {
      metric.resize(ThS.nv);
      metric = 0.;
    } else if (intopt[1] == 1) {
      metric.resize(6 * ThS.nv);
      metric = 0.;
    }
  }

  // mesh for yams
  yams_pSurfMesh yamsmesh;
  yamsmesh = (yams_pSurfMesh)calloc(1, sizeof(yams_SurfMesh));
  if (!yamsmesh) {
    cerr << "allocation error for SurfMesh for yams" << endl;
  }

  yamsmesh->infile = NULL;
  yamsmesh->outfile = NULL;
  yamsmesh->type = M_SMOOTH | M_QUERY | M_DETECT | M_BINARY | M_OUTPUT;

  meshS_to_yams_pSurfMesh(ThS, intopt[8], intopt[22], yamsmesh);

  // solution for freeyams2
  if (nbsol) {
    MeshPoint *mp3(MeshPointStack(stack));

    KN< bool > takemesh(nv);
    takemesh = false;

    for (int it = 0; it < nt; it++) {
      for (int iv = 0; iv < 3; iv++) {
        int i = ThS(it, iv);

        if (takemesh[i] == false) {
          mp3->setP(&ThS, it, iv);

          for (int ii = 0; ii < nbsolsize; ii++) {
            metric[i * nbsolsize + ii] = GetAny< double >((*sol[ii])(stack));
          }

          takemesh[i] = true;
        }
      }
    }
  }

  if (verbosity > 10) {
    cout << "nbsol  " << nargs[2] << endl;
  }

  if (nargs[2] || (nbsol > 0)) {
    float hmin, hmax;
    solyams_pSurfMesh(yamsmesh, mtype, metric, hmin, hmax);
    yamsmesh->nmfixe = yamsmesh->npfixe;
    if (fopt[7] < 0.0) {
      fopt[7] = max(fopt[7], hmin);
    }

    if (fopt[8] < 0.0) {
      fopt[8] = max(fopt[8], hmax);
    }
  } else {
    yamsmesh->nmfixe = 0;
  }

  int infondang = 0, infocc = 0;
  int res = yams_main(yamsmesh, intopt, fopt, infondang, infocc);
  if (verbosity > 10) {
    cout << " yamsmesh->dim " << yamsmesh->dim << endl;
  }

  if (res > 0) {
    cout << " problem with yams :: error " << res << endl;
    ExecError("Freeyams error");
  }

  MeshS *ThS_T = yams_pSurfMesh_to_meshS(yamsmesh, infondang, infocc, intopt[22]);

  // recuperer la solution ????
  if (verbosity > 10) {
    cout << &yamsmesh->point << " " << &yamsmesh->tria << " " << &yamsmesh->geom << " "
         << &yamsmesh->tgte << endl;
    cout << &yamsmesh << endl;
  }

  free(yamsmesh->point);
  free(yamsmesh->tria);
  free(yamsmesh->geom);
  free(yamsmesh->tgte);
  if (yamsmesh->metric) {
    free(yamsmesh->metric);
  }

  if (yamsmesh->edge) {
    free(yamsmesh->edge);
  }

  if (yamsmesh->tetra) {
    free(yamsmesh->tetra);
  }

  free(yamsmesh);

  *mp = mps;
  Add2StackOfPtr2FreeRC(stack, ThS_T);
  return SetAny< pmeshS >(ThS_T);
}

// version mesh3 in arg

class yams_Op_mesh3 : public E_F0mps {
 public:
  typedef pmeshS Result;
  Expression eTh;
  int nbsol;
  int nbsolsize;
  int type;
  int dim;
  vector< Expression > sol;

  static const int n_name_param = 14;    //
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

  int arg(int i, Stack stack, int a) const {
    return nargs[i] ? GetAny< long >((*nargs[i])(stack)) : a;
  }

  bool arg(int i, Stack stack, bool a) const {
    return nargs[i] ? GetAny< bool >((*nargs[i])(stack)) : a;
  }

 public:
  yams_Op_mesh3(const basicAC_F0 &args) : sol(args.size( ) - 1) {
    cout << "yams" << endl;
    args.SetNameParam(n_name_param, name_param, nargs);
    eTh = to< pmesh3 >(args[0]);
    dim = 3;
    nbsol = args.size( ) - 1;
    if (nbsol > 1) {
      CompileError(" yams accept only one solution ");
    }

    int ksol = 0;

    if (nbsol == 1) {
      int i = 1;
      if (args[i].left( ) == atype< E_Array >( )) {
        const E_Array *a = dynamic_cast< const E_Array * >(args[i].LeftValue( ));
        ffassert(a);
        ksol += a->size( );
      } else {
        ksol++;
      }

      sol.resize(ksol);

      // type :: 1 sca, 2 vector, 3 symtensor

      ksol = 0;
      nbsolsize = 0;
      type = 0;

      if (args[i].left( ) == atype< E_Array >( )) {
        const E_Array *a = dynamic_cast< const E_Array * >(args[i].LeftValue( ));
        ffassert(a);
        int N = a->size( );
        nbsolsize = nbsolsize + N;

        switch (N) {
          /*
           * case 3 :
           *  type[i-1]=2;
           *  for (int j=0;j<N;j++)
           *  sol[ksol++]=to<double>((*a)[j]);
           *  break;
           */
          case 6:
            type = 3;

            for (int j = 0; j < N; j++) {
              sol[ksol++] = to< double >((*a)[j]);
            }

            break;
          default:
            CompileError(
              " 3D solution for yams is a scalar (1 comp) or a symetric tensor (6 comp)");
            break;
        }
      } else {
        type = 1;
        nbsolsize = nbsolsize + 1;
        sol[ksol++] = to< double >(args[i]);
      }

      if (nargs[2]) {
        CompileError(" we give two metric for yams ");
      }
    }
  }

  static ArrayOfaType typeargs( ) { return ArrayOfaType(atype< pmesh3 >( ), true); }    // all type

  static E_F0 *f(const basicAC_F0 &args) { return new yams_Op_mesh3(args); }

  AnyType operator( )(Stack stack) const;
  operator aType( ) const { return atype< pmeshS >( ); }
};

basicAC_F0::name_and_type yams_Op_mesh3::name_param[] = {
  {"loptions", &typeid(KN_< long >)},    // 0
  {"doptions", &typeid(KN_< double >)},
  {"metric", &typeid(KN_< double >)},
  {"aniso", &typeid(bool)},    // 3
  {"mem", &typeid(long)},
  {"hmin", &typeid(double)},
  {"hmax", &typeid(double)},    // 6
  {"gradation", &typeid(double)},
  {"option", &typeid(long)},          // 8
  {"ridgeangle", &typeid(double)},    // 9
  {"absolute", &typeid(bool)},        // 10
  {"verbosity", &typeid(long)},       // 11

  {"nr", &typeid(long)},    // 12 no ridge
  {"ns", &typeid(long)}     // 13 no point smoothing
};

AnyType yams_Op_mesh3::operator( )(Stack stack) const {
  // initialisation
  MeshPoint *mp(MeshPointStack(stack)), mps = *mp;
  Mesh3 *pTh = GetAny< Mesh3 * >((*eTh)(stack));

  ffassert(pTh);
  Mesh3 &Th = *pTh;
  int nv = Th.nv;
  int nt = Th.nt;
  int nbe = Th.nbe;

  KN< int > defaultintopt(23);
  KN< double > defaultfopt(14);
  defaultintopt = 0;
  defaultfopt = 0.;
  yams_inival(defaultintopt, defaultfopt);

  KN< int > intopt(23);

  for (int ii = 0; ii < 23; ii++) {
    intopt[ii] = defaultintopt[ii];
  }

  KN< double > fopt(14);

  for (int ii = 0; ii < 14; ii++) {
    fopt[ii] = defaultfopt[ii];
  }

  assert(fopt.N( ) == 14);

  if (nargs[0]) {
    KN< int > intopttmp = GetAny< KN_< long > >((*nargs[0])(stack));
    if (intopttmp.N( ) != 13) {
      cerr << "the size of vector loptions is 13 " << endl;
      exit(1);
    } else {
      for (int ii = 0; ii < 13; ii++) {
        intopt[wrapper_intopt[ii]] = intopttmp[ii];
      }
    }
  }

  if (nargs[1]) {
    KN< double > fopttmp = GetAny< KN_< double > >((*nargs[1])(stack));
    if (fopttmp.N( ) != 11) {
      cerr << "the size of vector loptions is 11 not " << fopttmp.N( ) << endl;
      ExecError("FreeYams");
    } else {
      for (int ii = 0; ii < 11; ii++) {
        fopt[wrapper_fopt[ii]] = fopttmp[ii];
      }
    }
  }

  intopt[0] = arg(3, stack, intopt[0] != 1);
  intopt[8] = arg(4, stack, intopt[8]);
  fopt[7] = arg(5, stack, fopt[7]);
  fopt[8] = arg(6, stack, fopt[7]);
  fopt[6] = arg(7, stack, fopt[6]);
  intopt[22] = arg(8, stack, intopt[22]);    // optim option
  if (nargs[9]) {
    intopt[17] = 1;
  }

  fopt[13] = arg(9, stack, fopt[13]);             // ridge angle
  intopt[21] = arg(10, stack, intopt[21]);        // absolue
  intopt[11] = arg(11, stack, (int)verbosity);    // verbosity
  intopt[17] = arg(12, stack, intopt[17]);        // no ridge
  intopt[18] = arg(13, stack, intopt[18]);        // nb smooth
  if (verbosity > 1) {
    cout << " fopt = [";

    for (int i = 0; i < 11; ++i) {
      cout << fopt[wrapper_fopt[i]] << (i < 10 ? "," : "];\n");
    }

    cout << " intopt = [";

    for (int i = 0; i < 13; ++i) {
      cout << intopt[wrapper_intopt[i]] << (i < 12 ? "," : "];\n");
    }
  }

  KN< double > metric;

  int mtype = type;
  if (nargs[2]) {
    metric = GetAny< KN_< double > >((*nargs[2])(stack));
    if (metric.N( ) == Th.nv) {
      mtype = 1;
      intopt[1] = 0;
    } else if (metric.N( ) == 6 * Th.nv) {
      intopt[1] = 1;
      mtype = 3;
    } else {
      cerr << "sizeof vector metric is incorrect, size will be Th.nv or 6*Th.nv" << endl;
    }
  } else if (nbsol > 0) {
    if (type == 1) {
      intopt[1] = 0;
      metric.resize(Th.nv);
      metric = 0.;
    } else if (type == 3) {
      intopt[1] = 1;
      metric.resize(6 * Th.nv);
      metric = 0.;
    }
  } else {
    if (intopt[1] == 0) {
      metric.resize(Th.nv);
      metric = 0.;
    } else if (intopt[1] == 1) {
      metric.resize(6 * Th.nv);
      metric = 0.;
    }
  }

  // mesh for yams
  yams_pSurfMesh yamsmesh;
  yamsmesh = (yams_pSurfMesh)calloc(1, sizeof(yams_SurfMesh));
  if (!yamsmesh) {
    cerr << "allocation error for SurfMesh for yams" << endl;
  }

  yamsmesh->infile = NULL;
  yamsmesh->outfile = NULL;
  yamsmesh->type = M_SMOOTH | M_QUERY | M_DETECT | M_BINARY | M_OUTPUT;

  mesh3_to_yams_pSurfMesh(Th, intopt[8], intopt[22], yamsmesh);

  // solution for freeyams2
  if (nbsol) {
    MeshPoint *mp3(MeshPointStack(stack));

    KN< bool > takemesh(nv);
    takemesh = false;

    for (int it = 0; it < nt; it++) {
      for (int iv = 0; iv < 3; iv++) {
        int i = Th(it, iv);

        if (takemesh[i] == false) {
          mp3->setP(&Th, it, iv);

          for (int ii = 0; ii < nbsolsize; ii++) {
            metric[i * nbsolsize + ii] = GetAny< double >((*sol[ii])(stack));
          }

          takemesh[i] = true;
        }
      }
    }
  }

  if (verbosity > 10) {
    cout << "nbsol  " << nargs[2] << endl;
  }

  if (nargs[2] || (nbsol > 0)) {
    float hmin, hmax;
    solyams_pSurfMesh(yamsmesh, mtype, metric, hmin, hmax);
    yamsmesh->nmfixe = yamsmesh->npfixe;
    if (fopt[7] < 0.0) {
      fopt[7] = max(fopt[7], hmin);
    }

    if (fopt[8] < 0.0) {
      fopt[8] = max(fopt[8], hmax);
    }
  } else {
    yamsmesh->nmfixe = 0;
  }

  int infondang = 0, infocc = 0;
  int res = yams_main(yamsmesh, intopt, fopt, infondang, infocc);
  if (verbosity > 10) {
    cout << " yamsmesh->dim " << yamsmesh->dim << endl;
  }

  if (res > 0) {
    cout << " problem with yams :: error " << res << endl;
    ExecError("Freeyams error");
  }

  MeshS *ThS_T = yams_pSurfMesh_to_meshS(yamsmesh, infondang, infocc, intopt[22]);
  // Th3_T->getTypeMesh3()=1;
  // recuperer la solution ????
  if (verbosity > 10) {
    cout << &yamsmesh->point << " " << &yamsmesh->tria << " " << &yamsmesh->geom << " "
         << &yamsmesh->tgte << endl;
    cout << &yamsmesh << endl;
  }

  free(yamsmesh->point);
  free(yamsmesh->tria);
  free(yamsmesh->geom);
  free(yamsmesh->tgte);
  if (yamsmesh->metric) {
    free(yamsmesh->metric);
  }

  if (yamsmesh->edge) {
    free(yamsmesh->edge);
  }

  if (yamsmesh->tetra) {
    free(yamsmesh->tetra);
  }

  free(yamsmesh);

  *mp = mps;
  Add2StackOfPtr2FreeRC(stack, ThS_T);
  return SetAny< pmeshS >(ThS_T);
}

static void Load_Init( ) {    // le constructeur qui ajoute la fonction "splitmesh3"  a freefem++
  if (verbosity) {
    cout << " load: freeyams  " << endl;
  }

  Global.Add("freeyams", "(", new OneOperatorCode< yams_Op_mesh3 >);    //
  Global.Add("freeyams", "(", new OneOperatorCode< yams_Op_meshS >);
}

#define WITH_NO_INIT
#include "msh3.hpp"

LOADFUNC(Load_Init)
