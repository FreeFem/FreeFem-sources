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
// SUMMARY : Medit FreeFem++ link: popen
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

#include "mode_open.hpp"
#include "ff++.hpp"
#define WrdSiz 4

#ifdef _WIN32
string stringffmedit = "ffmedit.exe";
#else
string stringffmedit = "ffmedit";
#endif

const char *medit_popen = "-popen";    // 1";  // depend de l endroit ou se trouve medit
const char *medit_bin = "-filebin";
const char *medit_addsol = "-addsol";
const char *medit_debug = "-d";
static bool TheWait = false;
extern bool NoGraphicWindow;
using namespace std;
using namespace Fem2D;

// *******************************
//
// read solution of .sol or .solb
//
// *******************************

class readsol_Op : public E_F0mps {
 public:
  typedef KN_< double > Result;
  Expression eTh;
  Expression filename;
  static const int n_name_param = 1;
  static basicAC_F0::name_and_type name_param[];
  Expression nargs[n_name_param];
  long arg(int i, Stack stack, long a) const {
    return nargs[i] ? GetAny< long >((*nargs[i])(stack)) : a;
  }

 public:
  readsol_Op(const basicAC_F0 &args) {
    if (verbosity > 2) {
      cout << "readsol" << endl;
    }

    args.SetNameParam(n_name_param, name_param, nargs);
    if (BCastTo< string * >(args[0])) {
      filename = CastTo< string * >(args[0]);
    } else {
      CompileError("no filename given");
    }
  }

  static ArrayOfaType typeargs( ) { return ArrayOfaType(atype< string * >( )); }    // all type

  static E_F0 *f(const basicAC_F0 &args) { return new readsol_Op(args); }

  AnyType operator( )(Stack stack) const;
  operator aType( ) const { return atype< KN_< double > >( ); }
};

basicAC_F0::name_and_type readsol_Op::name_param[] = {{"number", &typeid(long)}};
AnyType readsol_Op::operator( )(Stack stack) const {
  string *ffname = GetAny< string * >((*filename)(stack));
  int k, i, isol, type,  ver, dim, typtab[GmfMaxTyp], offset;
  int64_t inm;
  char *ptr, data[128];
  // rajout freefem++
  int nv = 0, ntet = 0, ntri = 0;
  int nsol = 0;
  int key = 0;
  int numsol(arg(0, stack, -1L));

  assert(abs(numsol) >= 1);

  char *charfile = new char[ffname->size( ) + 1];
  strncpy(charfile, ffname->c_str( ), ffname->size( ) + 1);

  strcpy(data, charfile);
  ptr = strstr(data, ".sol");
  if (ptr) {
    *ptr = '\0';
  }

  strcat(data, ".solb");
  if (!(inm = GmfOpenMesh(data, GmfRead, &ver, &dim))) {
    ptr = strstr(data, ".solb");
    *ptr = '\0';
    strcat(data, ".sol");
    if (!(inm = GmfOpenMesh(data, GmfRead, &ver, &dim))) {
      cerr << "  ** " << (char *)data << " NOT FOUND.\n" << endl;
      exit(1);
    }
  }

  if (verbosity > 2) {
    cout << "  %%%%" << (char *)data << " OPENED\n" << endl;
  }

  nv = GmfStatKwd(inm, GmfSolAtVertices, &type, &offset, &typtab);
  if (nv) {
    key = GmfSolAtVertices;
    nsol = nv;
  } else {
    ntri = GmfStatKwd(inm, GmfSolAtTriangles, &type, &offset, &typtab);
    if (ntri) {
      key = GmfSolAtTriangles;
      nsol = ntri;
    } else {
      ntet = GmfStatKwd(inm, GmfSolAtTetrahedra, &type, &offset, &typtab);
    }

    if (ntet) {
      key = GmfSolAtTetrahedra;
      nsol = ntet;
    }
  }

  if ((nv == 0) && (ntri == 0) && (ntet == 0)) {
    cerr << "  ** MISSING DATA" << endl;
    exit(1);
  }

  int nbsol = nsol * offset;
  int offsettab = 0;
  int firstelem = 0;
  if (numsol != -1) {
    if (typtab[numsol - 1] == 1) {
      nbsol = nsol;
      offsettab = 1;
    } else if (typtab[numsol - 1] == 2) {
      nbsol = nsol * dim;
      offsettab = dim;
    } else if (typtab[numsol - 1] == 3) {
      nbsol = nsol * dim * (dim + 1) / 2;
      offsettab = dim * (dim + 1) / 2;
    } else {
      cerr << "bug in the definition of type of solution: 1 scalar, 2 vector, 3 symmetric tensor"
           << endl;
      exit(1);
    }

    for (int ii = 0; ii < (numsol - 1); ii++) {
      if (typtab[ii] == 1) {
        firstelem = firstelem + 1;
      } else if (typtab[ii] == 2) {
        firstelem = firstelem + dim;
      } else if (typtab[ii] == 3) {
        firstelem = firstelem + dim * (dim + 1) / 2;
        ;
      } else {
        cerr << "bug in the definition of type of solution: 1 scalar, 2 vector, 3 symmetric tensor"
             << endl;
        exit(1);
      }
    }
  }

  if (verbosity > 5) {
    cout << "nbsol " << nbsol << " offset " << offset << "  " << nsol << " " << endl;
  }

  float *buf = new float[offset];
  double tmp;
  double *bufd = new double[offset];

  KN< double > *ptabsol = new KN< double >(nbsol);
  KN< double > &tabsol = *ptabsol;

  if (numsol == -1) {
    GmfGotoKwd(inm, key);
    if (ver == GmfFloat) {
      for (k = 1; k <= nsol; k++) {
        isol = (k - 1) * offset;
        GmfGetLin(inm, key, buf);

        for (i = 0; i < offset; i++) {
          tabsol[isol + i] = (double)buf[i];
        }
      }
    } else {
      for (k = 1; k <= nsol; k++) {
        isol = (k - 1) * offset;
        GmfGetLin(inm, key, bufd);

        for (i = 0; i < offset; i++) {
          tabsol[isol + i] = bufd[i];
        }
      }
    }
  } else {
    GmfGotoKwd(inm, key);
    if (ver == GmfFloat) {
      for (k = 1; k <= nsol; k++) {
        isol = (k - 1) * offsettab;
        GmfGetLin(inm, key, buf);

        for (i = 0; i < offsettab; i++) {
          tabsol[isol + i] = buf[i + firstelem];
        }
      }
    } else {
      for (k = 1; k <= nsol; k++) {
        isol = (k - 1) * offsettab;
        GmfGetLin(inm, key, bufd);

        for (i = 0; i < offset; i++) {
          tabsol[isol + i] = bufd[i + firstelem];
        }
      }
    }
  }

  GmfCloseMesh(inm);
  delete[] buf;
  delete[] bufd;

  Add2StackOfPtr2Free(stack, ptabsol);
  return SetAny< KN< double > >(tabsol);
}

// *************************
//
// creation point sol
//
// *************************
// datasolMesh2

class datasolMesh2_Op : public E_F0mps {
 public:
  typedef long Result;
  Expression eTh;
  Expression filename;
  struct Expression2 {
    long what;       // 1 scalar, 2 vector, 3 symtensor
    long nbfloat;    // 1 scalar, 2 vector (3D), 3 symtensor(3D)
    string *dataname;
    Expression e[3];
    Expression2( ) {
      e[0] = 0;
      e[1] = 0;
      e[2] = 0;
      what = 0;
      nbfloat = 0;
    };
    Expression &operator[](int i) { return e[i]; }

    double eval(int i, Stack stack) const {
      if (e[i]) {
        return GetAny< double >((*e[i])(stack));
      } else {
        return 0;
      }
    }
  };
  vector< Expression2 > l;
  static const int n_name_param = 1;
  static basicAC_F0::name_and_type name_param[];
  Expression nargs[n_name_param];
  long arg(int i, Stack stack, long a) const {
    return nargs[i] ? GetAny< long >((*nargs[i])(stack)) : a;
  }

 public:
  datasolMesh2_Op(const basicAC_F0 &args) : l(args.size( ) - 2) {
    int ddim = 2;
    int stsize = 3;

    args.SetNameParam(n_name_param, name_param, nargs);

    if (BCastTo< string * >(args[0])) {
      filename = CastTo< string * >(args[0]);
    }

    if (BCastTo< pmesh >(args[1])) {
      eTh = CastTo< pmesh >(args[1]);
    }

    for (size_t i = 2; i < args.size( ); i++) {
      size_t jj = i - 2;

      if (BCastTo< double >(args[i])) {
        l[jj].what = 1;
        l[jj].nbfloat = 1;
        l[jj][0] = to< double >(args[i]);
      } else if (args[i].left( ) == atype< E_Array >( )) {
        const E_Array *a0 = dynamic_cast< const E_Array * >(args[i].LeftValue( ));
        if (a0->size( ) != ddim && a0->size( ) != stsize) {
          CompileError(
            "savesol in 2D: vector solution is 2 composant, tensor solution is 3 composant");
        }

        if (a0->size( ) == ddim) {
          // vector solution
          l[jj].what = 2;
          l[jj].nbfloat = ddim;

          for (int j = 0; j < ddim; j++) {
            l[jj][j] = to< double >((*a0)[j]);
          }
        } else if (a0->size( ) == stsize) {
          // symmetric tensor solution
          l[jj].what = 3;
          l[jj].nbfloat = stsize;

          for (int j = 0; j < stsize; j++) {
            l[jj][j] = to< double >((*a0)[j]);
          }
        }
      } else {
        cout << " arg " << i << " " << args[i].left( ) << endl;
        CompileError("savesol in 2D: Sorry no way to save this kind of data");
      }
    }
  }

  static ArrayOfaType typeargs( ) {
    return ArrayOfaType(atype< string * >( ), atype< pmesh >( ), true);
  }    // all type

  static E_F0 *f(const basicAC_F0 &args) { return new datasolMesh2_Op(args); }

  AnyType operator( )(Stack stack) const;
};

basicAC_F0::name_and_type datasolMesh2_Op::name_param[] = {{"order", &typeid(long)}};
AnyType datasolMesh2_Op::operator( )(Stack stack) const {
  MeshPoint *mp(MeshPointStack(stack));
  Mesh *pTh = GetAny< Mesh * >((*eTh)(stack));
  string *ffname = GetAny< string * >((*filename)(stack));

  ffassert(pTh);
  Mesh &Th = *pTh;
  int nt = Th.nt;
  int nv = Th.nv;
  int nbtype = l.size( );
  int nbsol;
  int solnbfloat;
  int TypTab[l.size( )];
  int resultorder = arg(0, stack, 1L);
  long longdefault;
    int ver = GmfFloat;
  int64_t outm;
  // determination de TypTab
  solnbfloat = 0;

  for (size_t i = 0; i < l.size( ); i++) {
    TypTab[i] = l[i].what;
    solnbfloat = solnbfloat + l[i].nbfloat;
  }

  float *OutSolTab = new float[solnbfloat];

  // determination de OutSolTab
  char *ret = new char[ffname->size( ) + 1];
  strcpy(ret, ffname->c_str( ));

  if (!(outm = GmfOpenMesh(ret, GmfWrite, ver, 2))) {
    cerr << "  -- Mesh3::Save  UNABLE TO OPEN  :" << ret << endl;
    exit(1);
  }

  if (resultorder == 0) {
    // ordre 0
    nbsol = nt;

    KN< double > valsol(solnbfloat * nbsol);

    MeshPoint *mp3(MeshPointStack(stack));
    R2 Cdg_hat = R2(1. / 3., 1. / 3.);

    for (int it = 0; it < nt; it++) {
      int h = 0;
      const Mesh::Triangle &K(Th.t(it));
      mp3->set(Th, K(Cdg_hat), Cdg_hat, K, K.lab);

      for (size_t i = 0; i < l.size( ); i++) {
        for (size_t j = 0; j < l[i].nbfloat; j++) {
          valsol[it * solnbfloat + h] = l[i].eval(j, stack);
          h = h + 1;
        }
      }

      assert(solnbfloat == h);
    }

    GmfSetKwd(outm, GmfSolAtTriangles, nbsol, nbtype, TypTab);

    for (int k = 0; k < nbsol; k++) {
      for (int i = 0; i < solnbfloat; i++) {
        OutSolTab[i] = valsol(k * solnbfloat + i);
      }

      GmfSetLin(outm, GmfSolAtTriangles, OutSolTab);
    }
  }

  if (resultorder == 1) {
    // ordre 1
    nbsol = nv;

    KN< double > valsol(solnbfloat * nbsol);
    valsol = 0.;
    KN< int > takemesh(nbsol);
    MeshPoint *mp3(MeshPointStack(stack));
    takemesh = 0;

    for (int it = 0; it < nt; it++) {
      for (int iv = 0; iv < 3; iv++) {
        int i = Th(it, iv);

        mp3->setP(&Th, it, iv);
        int h = 0;

        for (size_t ii = 0; ii < l.size( ); ii++) {
          for (size_t j = 0; j < l[ii].nbfloat; j++) {
            valsol[i * solnbfloat + h] = valsol[i * solnbfloat + h] + l[ii].eval(j, stack);
            h = h + 1;
          }
        }

        assert(solnbfloat == h);
        takemesh[i] = takemesh[i] + 1;
      }
    }

    for (int i = 0; i < nv; i++) {
      for (int h = 0; h < solnbfloat; h++) {
        valsol[i * solnbfloat + h] = valsol[i * solnbfloat + h] / takemesh[i];
      }
    }

    GmfSetKwd(outm, GmfSolAtVertices, nbsol, nbtype, TypTab);

    for (int k = 0; k < nbsol; k++) {
      for (int i = 0; i < solnbfloat; i++) {
        OutSolTab[i] = valsol(k * solnbfloat + i);
      }

      GmfSetLin(outm, GmfSolAtVertices, OutSolTab);
    }
  }

  GmfCloseMesh(outm);
  delete[] ret;
  delete[] OutSolTab;
  return longdefault;
}

// datasolMeshS for 3D surface
template< class MMesh, class v_fes >
class datasolMeshT_Op : public E_F0mps {
 public:
  typedef long Result;
  typedef const MMesh *ppmesh;
  Expression eTh;
  Expression filename;
  struct Expression2 {
    long what;       // 1 scalar, 2 vector, 3 symtensor
    long nbfloat;    // 1 scalar, 3 vector (3D), 6 symtensor(3D)
    Expression e[6];
    Expression2( ) {
      e[0] = 0;
      e[1] = 0;
      e[2] = 0;
      e[3] = 0;
      what = 0;
      nbfloat = 0;
    };
    Expression &operator[](int i) { return e[i]; }

    double eval(int i, Stack stack) const {
      if (e[i]) {
        return GetAny< double >((*e[i])(stack));
      } else {
        return 0;
      }
    }
  };
  vector< Expression2 > l;
  static const int n_name_param = 1;
  static basicAC_F0::name_and_type name_param[];
  Expression nargs[n_name_param];
  long arg(int i, Stack stack, long a) const {
    return nargs[i] ? GetAny< long >((*nargs[i])(stack)) : a;
  }

 public:
  datasolMeshT_Op(const basicAC_F0 &args) : l(args.size( ) - 2) {
    int ddim = 3;
    int stsize = 3;

    args.SetNameParam(n_name_param, name_param, nargs);

    filename = CastTo< string * >(args[0]);
    eTh = CastTo< ppmesh >(args[1]);

    for (size_t i = 2; i < args.size( ); i++) {
      size_t jj = i - 2;

      if (BCastTo< double >(args[i])) {
        l[jj].what = 1;
        l[jj].nbfloat = 1;
        l[jj][0] = to< double >(args[i]);
      } else if (args[i].left( ) == atype< E_Array >( )) {
        const E_Array *a0 = dynamic_cast< const E_Array * >(args[i].LeftValue( ));
        if (a0->size( ) != ddim && a0->size( ) != stsize) {
          if (a0 == nullptr) printf("Dynamic cast failed\n");
          CompileError(
            "savesol in 3D surface: vector solution is 3 composant, vector solution is 6 "
            "composant");
        }

        if (a0->size( ) == ddim) {
          // vector solution
          l[jj].what = 2;
          l[jj].nbfloat = ddim;

          for (int j = 0; j < ddim; j++) {
            l[jj][j] = to< double >((*a0)[j]);
          }
        } else if (a0->size( ) == stsize) {
          // symmetric tensor solution
          l[jj].what = 3;
          l[jj].nbfloat = stsize;

          for (int j = 0; j < stsize; j++) {
            l[jj][j] = to< double >((*a0)[j]);
          }
        }
      } else {
        CompileError("savesol in 3D surface: Sorry no way to save this kind of data");
      }
    }
  }

  static ArrayOfaType typeargs( ) {
    return ArrayOfaType(atype< string * >( ), atype< ppmesh >( ), true);
  }    // all type

  static E_F0 *f(const basicAC_F0 &args) { return new datasolMeshT_Op< MMesh, v_fes >(args); }

  AnyType operator( )(Stack stack) const;
};

template<>
basicAC_F0::name_and_type datasolMeshT_Op< MeshS, v_fesS >::name_param[] = {
  {"order", &typeid(long)}};
template<>
basicAC_F0::name_and_type datasolMeshT_Op< MeshL, v_fesL >::name_param[] = {
  {"order", &typeid(long)}};

template< class MMesh, class v_fes >
AnyType datasolMeshT_Op< MMesh, v_fes >::operator( )(Stack stack) const {
  MeshPoint *mp(MeshPointStack(stack));
  MMesh *pTh = GetAny< MMesh * >((*eTh)(stack));
  string *ffname = GetAny< string * >((*filename)(stack));
  typedef typename MMesh::Element T;
  typedef typename MMesh::BorderElement B;
  typedef typename MMesh::Vertex V;
  typedef typename MMesh::Element::RdHat TRdHat;
  typedef typename MMesh::BorderElement::RdHat BRdHat;

  ffassert(pTh);
  MMesh &Th = *pTh;
  int nt = Th.nt;
  int nv = Th.nv;
  int nbtype = l.size( );
  int nbsol;
  int solnbfloat;
  int TypTab[l.size( )];
  int resultorder = arg(0, stack, 1);
  long longdefault;
    int ver = GmfFloat;
  int64_t outm;
  // determination de TypTab
  solnbfloat = 0;

  for (size_t i = 0; i < l.size( ); i++) {
    TypTab[i] = l[i].what;
    solnbfloat = solnbfloat + l[i].nbfloat;
  }

  float *OutSolTab = new float[solnbfloat];

  // determination de OutSolTab
  char *ret = new char[ffname->size( ) + 1];
  strcpy(ret, ffname->c_str( ));
  if (verbosity > 2) {
    cout << ret << endl;
  }

  if (!(outm = GmfOpenMesh(ret, GmfWrite, ver, 3))) {
    cerr << "  -- MeshS::Save  UNABLE TO OPEN  :" << filename << endl;
    exit(1);
  }
  double k = T::nv;
  TRdHat Cdg_hat = TRdHat::diag(1. / k);

  if (resultorder == 0) {
    // Element
    // ordre 0
    nbsol = nt;
    KN< double > valsol(solnbfloat * nbsol);
    MeshPoint *mp3(MeshPointStack(stack));

    for (int it = 0; it < nt; it++) {
      int h = 0;
      const T &K(Th.elements[it]);
      mp3->set(Th, K(Cdg_hat), Cdg_hat, K, K.lab );

      for (size_t i = 0; i < l.size( ); i++)
        for (size_t j = 0; j < l[i].nbfloat; j++) {
          valsol[it * solnbfloat + h] = l[i].eval(j, stack);
          h = h + 1;
        }

      assert(solnbfloat == h);
    }
    if (std::is_same< MMesh, MeshS >::value)
      GmfSetKwd(outm, GmfSolAtTriangles, nbsol, nbtype, TypTab);
    else if (std::is_same< MMesh, MeshL >::value)
      GmfSetKwd(outm, GmfSolAtEdges, nbsol, nbtype, TypTab);

    for (int k = 0; k < nbsol; k++) {
      for (int i = 0; i < solnbfloat; i++) {
        OutSolTab[i] = valsol(k * solnbfloat + i);
      }
      if (std::is_same< MMesh, MeshS >::value)
        GmfSetLin(outm, GmfSolAtTriangles, OutSolTab);
      else if (std::is_same< MMesh, MeshL >::value)
        GmfSetLin(outm, GmfSolAtEdges, OutSolTab);
    }
  }

  if (resultorder == 1) {
    // ordre 1
    nbsol = nv;

    KN< double > valsol(solnbfloat * nbsol);
    KN< int > takemesh(nbsol);
    MeshPoint *mp3(MeshPointStack(stack));
    takemesh = 0;

    for (int it = 0; it < nt; it++)
      for (int iv = 0; iv < T::nv; iv++) {
        int i = Th(it, iv);

        if (takemesh[i] == 0) {
          mp3->setP(&Th, it, iv);
          int h = 0;

          for (size_t ii = 0; ii < l.size( ); ii++) {
            for (size_t j = 0; j < l[ii].nbfloat; j++) {
              valsol[i * solnbfloat + h] = l[ii].eval(j, stack);
              h = h + 1;
            }
          }

          assert(solnbfloat == h);
          takemesh[i] = takemesh[i] + 1;
        }
      }

    GmfSetKwd(outm, GmfSolAtVertices, nbsol, nbtype, TypTab);

    for (int k = 0; k < nbsol; k++) {
      for (int i = 0; i < solnbfloat; i++) {
        OutSolTab[i] = valsol(k * solnbfloat + i);
      }

      GmfSetLin(outm, GmfSolAtVertices, OutSolTab);
    }
  }

  GmfCloseMesh(outm);
  delete[] ret;
  delete[] OutSolTab;
  return longdefault;
}

// datasolMesh3 for 3D volume
template< class v_fes >
class datasolMesh3_Op : public E_F0mps {
 public:
  typedef long Result;
  Expression eTh;
  Expression filename;
  struct Expression2 {
    long what;       // 1 scalar, 2 vector, 3 symtensor
    long nbfloat;    // 1 scalar, 3 vector (3D), 6 symtensor(3D)
    Expression e[6];
    Expression2( ) {
      e[0] = 0;
      e[1] = 0;
      e[2] = 0;
      e[3] = 0;
      e[4] = 0;
      e[5] = 0;
      what = 0;
      nbfloat = 0;
    };
    Expression &operator[](int i) { return e[i]; }

    double eval(int i, Stack stack) const {
      if (e[i]) {
        return GetAny< double >((*e[i])(stack));
      } else {
        return 0;
      }
    }
  };
  vector< Expression2 > l;
  static const int n_name_param = 1;
  static basicAC_F0::name_and_type name_param[];
  Expression nargs[n_name_param];
  long arg(int i, Stack stack, long a) const {
    return nargs[i] ? GetAny< long >((*nargs[i])(stack)) : a;
  }

 public:
  datasolMesh3_Op(const basicAC_F0 &args) : l(args.size( ) - 2) {
    int ddim = 3;
    int stsize = 6;

    args.SetNameParam(n_name_param, name_param, nargs);

    filename = CastTo< string * >(args[0]);
    eTh = CastTo< pmesh3 >(args[1]);

    for (size_t i = 2; i < args.size( ); i++) {
      size_t jj = i - 2;

      if (BCastTo< double >(args[i])) {
        l[jj].what = 1;
        l[jj].nbfloat = 1;
        l[jj][0] = to< double >(args[i]);
      } else if (args[i].left( ) == atype< E_Array >( )) {
        const E_Array *a0 = dynamic_cast< const E_Array * >(args[i].LeftValue( ));
        if (a0->size( ) != ddim && a0->size( ) != stsize) {
          if (a0 == nullptr) cout << "dynamic cast failed" << endl;
          CompileError(
            "savesol in 3D: vector solution is 3 composant, vector solution is 6 composant");
        }

        if (a0->size( ) == ddim) {
          // vector solution
          l[jj].what = 2;
          l[jj].nbfloat = ddim;

          for (int j = 0; j < ddim; j++) {
            l[jj][j] = to< double >((*a0)[j]);
          }
        } else if (a0->size( ) == stsize) {
          // symmetric tensor solution
          l[jj].what = 3;
          l[jj].nbfloat = stsize;

          for (int j = 0; j < stsize; j++) {
            l[jj][j] = to< double >((*a0)[j]);
          }
        }
      } else {
        CompileError("savesol in 3D: Sorry no way to save this kind of data");
      }
    }
  }

  static ArrayOfaType typeargs( ) {
    return ArrayOfaType(atype< string * >( ), atype< pmesh3 >( ), true);
  }    // all type

  static E_F0 *f(const basicAC_F0 &args) { return new datasolMesh3_Op(args); }

  AnyType operator( )(Stack stack) const;
};

template< class v_fes >
basicAC_F0::name_and_type datasolMesh3_Op< v_fes >::name_param[] = {{"order", &typeid(long)}};
template< class v_fes >
AnyType datasolMesh3_Op< v_fes >::operator( )(Stack stack) const {
  MeshPoint *mp(MeshPointStack(stack));
  Mesh3 *pTh = GetAny< Mesh3 * >((*eTh)(stack));
  string *ffname = GetAny< string * >((*filename)(stack));

  ffassert(pTh);
  Mesh3 &Th = *pTh;
  int nt = Th.nt;
  int nv = Th.nv;
  int nbtype = l.size( );
  int nbsol;
  int solnbfloat;
  int TypTab[l.size( )];
  int resultorder = arg(0, stack, 1);
  long longdefault;
    int ver = GmfFloat;
  int64_t outm;
  // determination de TypTab
  solnbfloat = 0;

  for (size_t i = 0; i < l.size( ); i++) {
    TypTab[i] = l[i].what;
    solnbfloat = solnbfloat + l[i].nbfloat;
  }

  float *OutSolTab = new float[solnbfloat];

  // determination de OutSolTab
  char *ret = new char[ffname->size( ) + 1];
  strcpy(ret, ffname->c_str( ));
  if (verbosity > 2) {
    cout << ret << endl;
  }

  if (!(outm = GmfOpenMesh(ret, GmfWrite, ver, 3))) {
    cerr << "  -- Mesh3::Save  UNABLE TO OPEN  :" << filename << endl;
    exit(1);
  }

  if (resultorder == 0) {
    // Tetrahedra
    // ordre 0
    nbsol = nt;
    KN< double > valsol(solnbfloat * nbsol);

    MeshPoint *mp3(MeshPointStack(stack));
    R3 Cdg_hat = R3(1. / 4., 1. / 4., 1. / 4.);

    for (int it = 0; it < nt; it++) {
      int h = 0;
      const Tet &K(Th.elements[it]);
      mp3->set(Th, K(Cdg_hat), Cdg_hat, K, K.lab);

      for (size_t i = 0; i < l.size( ); i++) {
        for (size_t j = 0; j < l[i].nbfloat; j++) {
          valsol[it * solnbfloat + h] = l[i].eval(j, stack);
          h = h + 1;
        }
      }

      assert(solnbfloat == h);
    }

    GmfSetKwd(outm, GmfSolAtTetrahedra, nbsol, nbtype, TypTab);

    for (int k = 0; k < nbsol; k++) {
      for (int i = 0; i < solnbfloat; i++) {
        OutSolTab[i] = valsol(k * solnbfloat + i);
      }

      GmfSetLin(outm, GmfSolAtTetrahedra, OutSolTab);
    }
  }

  if (resultorder == 1) {
    // ordre 1
    nbsol = nv;

    KN< double > valsol(solnbfloat * nbsol);
    KN< int > takemesh(nbsol);
    MeshPoint *mp3(MeshPointStack(stack));
    R3 Cdg_hat = R3(1. / 4., 1. / 4., 1. / 4.);
    takemesh = 0;

    for (int it = 0; it < nt; it++) {
      for (int iv = 0; iv < 4; iv++) {
        int i = Th(it, iv);

        if (takemesh[i] == 0) {
          mp3->setP(&Th, it, iv);
          int h = 0;

          for (size_t ii = 0; ii < l.size( ); ii++) {
            for (size_t j = 0; j < l[ii].nbfloat; j++) {
              valsol[i * solnbfloat + h] = l[ii].eval(j, stack);
              h = h + 1;
            }
          }

          assert(solnbfloat == h);
          takemesh[i] = takemesh[i] + 1;
        }
      }
    }

    GmfSetKwd(outm, GmfSolAtVertices, nbsol, nbtype, TypTab);

    for (int k = 0; k < nbsol; k++) {
      for (int i = 0; i < solnbfloat; i++) {
        OutSolTab[i] = valsol(k * solnbfloat + i);
      }

      GmfSetLin(outm, GmfSolAtVertices, OutSolTab);
    }
  }

  GmfCloseMesh(outm);
  delete[] ret;
  delete[] OutSolTab;
  return longdefault;
}

// *************************
//
// medit
//
// *************************

static char *meditcmd(long filebin, int nbsol, int smedit, const string &meditff,
                      const string &ffnn) {
  string meditcmm = meditff;
  int ddebug = 0;

  if (ddebug) {
    meditcmm += ' ';
    meditcmm += medit_debug;
  }

  meditcmm += ' ';
  meditcmm += medit_popen;
  if (filebin) {
    meditcmm += ' ';
    meditcmm += medit_bin;
  }

  if (nbsol) {
    meditcmm += ' ';
    meditcmm += medit_addsol;
  }

  char meditsol[5];
  sprintf(meditsol, " %i", smedit);
  meditcmm += meditsol;

  meditcmm += ' ';

  KN< char > ret1(ffnn.size( ) + 1);
  strcpy(ret1, ffnn.c_str( ));

  int nbstrings = 1;
  char *tictac;
  tictac = strtok(ret1, " \n");

  meditcmm += ' ';
  meditcmm += tictac;

  while (tictac != NULL && nbstrings < nbsol) {
    tictac = strtok(NULL, " \n");
    meditcmm += ' ';
    meditcmm += tictac;
    nbstrings++;
  }

  if (nbstrings != smedit) {
    cout
      << "The number of string defined in string parameter is different of the number of solution"
      << endl;
    if (nbstrings < smedit) {
      // Add strings
      while (nbstrings < smedit) {
        nbstrings++;
        char newsol[10];
        sprintf(newsol, " ffsol%i", nbstrings);
        meditcmm += newsol;
      }
    }
  }

  char *ret = new char[meditcmm.size( ) + 1];
  // char ret[meditcmm.size()+1];
  strcpy(ret, meditcmm.c_str( ));
  return ret;
}

void writetabsol(const int &tsize, const int &nbofsol, const KN< double > &v1, KNM< double > &vv) {
  for (int i = 0; i < tsize; i++) {
    vv(nbofsol, i) = v1(i);
  }
}

void writetabsol(const int &tsize, const int &nbofsol, const KN< double > &v1,
                 const KN< double > &v2, KNM< double > &vv) {
  for (int i = 0; i < tsize; i++) {
    vv(nbofsol, i) = v1(i);
    vv(nbofsol + 1, i) = v2(i);
  }
}

void writetabsol(const int &tsize, const int &nbofsol, const KN< double > &v1,
                 const KN< double > &v2, const KN< double > &v3, KNM< double > &vv) {
  for (int i = 0; i < tsize; i++) {
    vv(nbofsol, i) = v1(i);
    vv(nbofsol + 1, i) = v2(i);
    vv(nbofsol + 2, i) = v3(i);
  }
}

void writetabsol(const int &tsize, const int &nbofsol, const KN< double > &v1,
                 const KN< double > &v2, const KN< double > &v3, const KN< double > &v4,
                 const KN< double > &v5, const KN< double > &v6, KNM< double > &vv) {
  for (int i = 0; i < tsize; i++) {
    vv(nbofsol, i) = v1(i);
    vv(nbofsol + 1, i) = v2(i);
    vv(nbofsol + 2, i) = v3(i);
    vv(nbofsol + 3, i) = v4(i);
    vv(nbofsol + 4, i) = v5(i);
    vv(nbofsol + 5, i) = v6(i);
  }
}

class PopenMeditMesh_Op : public E_F0mps {
 public:
  typedef long Result;
  Expression eTh;
  Expression filename;
  long offset;
  long nbTh;
  struct Expression2 {
    long what;       // 0 mesh, 1 scalar, 2 vector, 3 symtensor
    long nbfloat;    // 1 scalar, 2 vector (2D), 3 symtensor(2D)
    Expression e[3];
    Expression2( ) {
      e[0] = 0;
      e[1] = 0;
      e[2] = 0;
      what = 0;
      nbfloat = 0;
    };
    Expression &operator[](int i) { return e[i]; }

    double eval(int i, Stack stack) const {
      if (e[i]) {
        return GetAny< double >((*e[i])(stack));
      } else {
        return 0;
      }
    }

    const Mesh &evalm(int i, Stack stack) const {
      throwassert(e[i]);
      return *GetAny< pmesh >((*e[i])(stack));
    }
  };
  vector< Expression2 > l;

  static const int n_name_param = 5;
  static basicAC_F0::name_and_type name_param[];
  Expression nargs[n_name_param];
  long arg(int i, Stack stack, long a) const {
    return nargs[i] ? GetAny< long >((*nargs[i])(stack)) : a;
  }

  string *arg(int i, Stack stack, string *a) const {
    return nargs[i] ? GetAny< string * >((*nargs[i])(stack)) : a;
  }

 public:
  PopenMeditMesh_Op(const basicAC_F0 &args) : l(args.size( ) - 1) {
    int nbofsol;
    int ddim = 2;
    int stsize = 3;
    char *tictac;

    args.SetNameParam(n_name_param, name_param, nargs);

    if (BCastTo< string * >(args[0])) {
      filename = CastTo< string * >(args[0]);
    }

    for (size_t i = 1; i < args.size( ); i++) {
      size_t jj = i - 1;

      if (BCastTo< double >(args[i])) {
        l[jj].what = 1;
        l[jj].nbfloat = 1;
        l[jj][0] = to< double >(args[i]);
      } else if (args[i].left( ) == atype< E_Array >( )) {
        const E_Array *a0 = dynamic_cast< const E_Array * >(args[i].LeftValue( ));
        if (a0 == nullptr) cout << "dynamic cast error" << endl;
        if (a0->size( ) != ddim && a0->size( ) != stsize) {
          CompileError(
            "medit in 2D: vector solution is 2 composant, tensor solution is 3 composant");
        }

        if (a0->size( ) == ddim) {
          // vector solution
          l[jj].what = 2;
          l[jj].nbfloat = ddim;

          for (int j = 0; j < ddim; j++) {
            l[jj][j] = to< double >((*a0)[j]);
          }
        } else if (a0->size( ) == stsize) {
          // symmetric tensor solution
          l[jj].what = 3;
          l[jj].nbfloat = stsize;

          for (int j = 0; j < stsize; j++) {
            l[jj][j] = to< double >((*a0)[j]);
          }
        }
      } else if (BCastTo< pmesh >(args[i])) {
        l[jj].what = 0;
        l[jj].nbfloat = 0;
        l[jj][0] = CastTo< pmesh >(args[i]);
      } else {
        CompileError("medit in 2D: Sorry no way to save this kind of data");
      }
    }

    offset = 0;
    nbTh = 1;

    // determination of the number of solutions
    // =============================
    // 0 2 3 2 ! 0 2 3 2 ! 0 2  3  2
    // 0 1 2 3 ! 4 5 6 7 ! 8 9 10 11
    // =============================
    for (size_t jj = 1; jj < l.size( ); jj++) {
      if (l[jj].what == 0 && offset == 0) {
        offset = jj;
      }

      if (l[jj].what == 0) {
        nbTh++;
        if (jj != (nbTh - 1) * offset) {
          CompileError("the number of solution by mesh is different");
        }
      }
    }

    /*
     * if( offset-1 != nbstrings ){
     * CompileError("The number of string defined in string parameter is different of the number of
     * solution");
     * }
     */

    if (nbTh == 1) {
      // case of one mesh
      offset = l.size( );
    } else {
      // analyse of the different solution
      for (size_t jj = offset; jj < l.size( ); jj++) {
        if (l[jj].what != l[jj % offset].what) {
          char StringError[256];
          snprintf(StringError, 256,
                   "compile error ::  The solution %ld of mesh 1 and mesh %ld is not the same type",
                   jj % offset, jj / offset + 1);
          CompileError(StringError);
        }
      }
    }

    // verification que la nature des solutions sont identiques pour les differents maillages.
  }

  static ArrayOfaType typeargs( ) {
    return ArrayOfaType(atype< string * >( ), atype< pmesh >( ), true);
  }    // all type

  static E_F0 *f(const basicAC_F0 &args) { return new PopenMeditMesh_Op(args); }

  AnyType operator( )(Stack stack) const;
};

basicAC_F0::name_and_type PopenMeditMesh_Op::name_param[] = {{"order", &typeid(long)},
                                                             {"meditff", &typeid(string *)},
                                                             {"save", &typeid(string *)},
                                                             {"wait", &typeid(bool)},
                                                             {"bin", &typeid(long)}};
AnyType PopenMeditMesh_Op::operator( )(Stack stack) const {
  if (NoGraphicWindow) {
    return Nothing;
  }

  MeshPoint *mp(MeshPointStack(stack));
  long order(arg(0, stack, 1));
  int ver = GmfFloat;
  int dimp = 2;
  float fx, fy;
  long valsortie = 0;
  int typsol = 0, nbsol;
  nbsol = offset - 1;

  int TypTab[l.size( ) - 1];

  for (size_t i = 0; i < l.size( ) - 1; i++) {
    TypTab[i] = l[i + 1].what;
  }

  string *ffname = GetAny< string * >((*filename)(stack));
  string *meditff(arg(1, stack, &stringffmedit));
  long filebin(arg(4, stack, 1));
  int smedit = max(1, nbsol);
  char *commandline = meditcmd(filebin, nbsol, smedit, *meditff, *ffname);
  printf("version de medit %s\n", commandline);

  // lecture des differents maillages
  int nv = 0, nt = 0, nbe = 0;    // sommet, triangles, arretes du maillage unifies

  for (size_t i = 0; i < l.size( ); i = i + offset) {
    if (l[i].what != 0) {
      cerr << "this element is not a mesh" << i << endl;
    }

    const Mesh &Thtmp = l[i].evalm(0, stack);
    nv += Thtmp.nv;
    nt += Thtmp.nt;
    nbe += Thtmp.neb;
  }

  Mesh::Vertex *v = new Mesh::Vertex[nv];
  Mesh::Triangle *t = new Mesh::Triangle[nt];
  Mesh::BorderElement *b = new Mesh::BorderElement[nbe];
  Mesh::Triangle *tt = t;
  Mesh::BorderElement *bb = b;
  int iv, it, ibe;
  iv = 0;
  it = 0;
  ibe = 0;

  int numTht[nt];    // numero of Th assoctiated with a triangles
  int jt = 0;

  for (size_t i = 0; i < l.size( ); i = i + offset) {
    int nvtmp = iv;
    const Mesh &Thtmp = l[i].evalm(0, stack);

    for (int ii = 0; ii < Thtmp.nv; ii++) {
      const Mesh::Vertex &vi(Thtmp(ii));
      v[iv] = vi;
      iv++;
    }

    for (int ii = 0; ii < Thtmp.nt; ii++) {
      const Mesh::Triangle &vi(Thtmp.t(ii));
      int i0 = nvtmp + Thtmp(ii, 0);
      int i1 = nvtmp + Thtmp(ii, 1);
      int i2 = nvtmp + Thtmp(ii, 2);
      (*tt++).set(v, i0, i1, i2, vi.lab);
      numTht[it] = jt;
      it++;
    }

    for (int ii = 0; ii < Thtmp.neb; ii++) {
      const Mesh::BorderElement &vi(Thtmp.be(ii));    // const BoundaryEdge &vi(Thtmp(ii));
      int i0 = nvtmp + Thtmp.operator( )(vi[0]);
      int i1 = nvtmp + Thtmp.operator( )(vi[1]);
      (*bb++).set(v, i0, i1, vi.lab);
      ibe++;
    }

    jt++;
  }

  assert(it == nt);
  assert(iv == nv);
  assert(ibe = nbe);
  if (verbosity > 2) {
    cout << "Popen medit : vertex " << nv << " triangle " << nt << " edge " << nbe << endl;
  }

  Mesh *pTh = new Mesh(nv, nt, nbe, v, t, b);
  Mesh &Th = *pTh;

  // determination of the number of elements to represent the solution
  int datasize;
  if (order == 0) {
    datasize = nt;
  }

  if (order == 1) {
    datasize = nv;
  }

  // cas de sauvegarde
  bool boolsave = false;
  int solnbfloat = 0;
  KNM< double > solsave(1, 1);
  string *saveff;
  KN< double > vxx = 0., vyx = 0., vyy = 0.;

  if (nbsol > 0) {
    vxx.init(datasize);
    vyx.init(datasize);
    vyy.init(datasize);

    if (nargs[2]) {
      boolsave = true;
      saveff = GetAny< string * >((*nargs[2])(stack));
      int ddim = 2;

      for (size_t i = 0; i < offset; i++) {
        solnbfloat = solnbfloat + l[i].nbfloat;
      }

      solsave.init(solnbfloat, datasize);
      solsave = 0.;
    }
  }

  int nboftmp = 0;
  FILE *popenstream = popen(commandline, MODE_WRITE_BINARY);
  if (!popenstream) {
    cerr << " Error popen : " << commandline << endl;
    exit(1);
  }

  // mesh
  int jojo1;

  for (int jojo = 0; jojo < smedit; jojo++) {
    if (filebin) {
      int cod = 1;
      int KwdCod;
      int NulPos = 0;

      // determination of number solutions associated with a mesh
      fwrite((unsigned char *)&cod, WrdSiz, 1, popenstream);
      fwrite((unsigned char *)&ver, WrdSiz, 1, popenstream);
      KwdCod = GmfDimension;
      fwrite((unsigned char *)&KwdCod, WrdSiz, 1, popenstream);
      fwrite((unsigned char *)&NulPos, WrdSiz, 1, popenstream);
      fwrite((unsigned char *)&dimp, WrdSiz, 1, popenstream);

      // vertex
      KwdCod = GmfVertices;
      fwrite((unsigned char *)&KwdCod, WrdSiz, 1, popenstream);
      fwrite((unsigned char *)&NulPos, WrdSiz, 1, popenstream);
      fwrite((unsigned char *)&nv, WrdSiz, 1, popenstream);

      for (int k = 0; k < nv; k++) {
        const Mesh::Vertex &P = Th.vertices[k];
        fx = P.x;
        fy = P.y;
        fwrite((unsigned char *)&fx, WrdSiz, 1, popenstream);
        fwrite((unsigned char *)&fy, WrdSiz, 1, popenstream);
        fwrite((unsigned char *)&(P.lab), WrdSiz, 1, popenstream);
      }

      // triangles
      KwdCod = GmfTriangles;
      fwrite((unsigned char *)&KwdCod, WrdSiz, 1, popenstream);
      fwrite((unsigned char *)&NulPos, WrdSiz, 1, popenstream);
      fwrite((unsigned char *)&nt, WrdSiz, 1, popenstream);

      for (int k = 0; k < nt; k++) {
        const Mesh::Triangle &K(Th.t(k));
        int i0 = Th.operator( )(K[0]) + 1;
        int i1 = Th.operator( )(K[1]) + 1;
        int i2 = Th.operator( )(K[2]) + 1;
        int lab = K.lab;
        fwrite((unsigned char *)&i0, WrdSiz, 1, popenstream);
        fwrite((unsigned char *)&i1, WrdSiz, 1, popenstream);
        fwrite((unsigned char *)&i2, WrdSiz, 1, popenstream);
        fwrite((unsigned char *)&lab, WrdSiz, 1, popenstream);
      }

      // Edges
      KwdCod = GmfEdges;
      fwrite((unsigned char *)&KwdCod, WrdSiz, 1, popenstream);
      fwrite((unsigned char *)&NulPos, WrdSiz, 1, popenstream);
      fwrite((unsigned char *)&nbe, WrdSiz, 1, popenstream);

      for (int k = 0; k < nbe; k++) {
        const Mesh::BorderElement &K(Th.be(k));
        int i0 = Th.operator( )(K[0]) + 1;
        int i1 = Th.operator( )(K[1]) + 1;
        int lab = K.lab;
        fwrite((unsigned char *)&i0, WrdSiz, 1, popenstream);
        fwrite((unsigned char *)&i1, WrdSiz, 1, popenstream);
        fwrite((unsigned char *)&lab, WrdSiz, 1, popenstream);
      }

      // End
      KwdCod = GmfEnd;
      fwrite((unsigned char *)&KwdCod, WrdSiz, 1, popenstream);
      fwrite((unsigned char *)&NulPos, WrdSiz, 1, popenstream);
    } else {
      // determination of number solutions associated with a mesh
      fprintf(popenstream, "MeshVersionFormatted\n");
      fprintf(popenstream, "%i\n", ver);
      fprintf(popenstream, "Dimension\n");
      fprintf(popenstream, "%i\n", dimp);

      fprintf(popenstream, "Vertices\n");
      fprintf(popenstream, "%i\n", nv);

      for (int k = 0; k < nv; k++) {
        const Mesh::Vertex &P = Th.vertices[k];
        fx = P.x;
        fy = P.y;
        fprintf(popenstream, "%f %f %i\n", fx, fy, P.lab);
      }

      fprintf(popenstream, "Triangles\n");
      fprintf(popenstream, "%i\n", nt);

      for (int k = 0; k < nt; k++) {
        const Mesh::Triangle &K(Th.t(k));
        int i0 = Th.operator( )(K[0]) + 1;
        int i1 = Th.operator( )(K[1]) + 1;
        int i2 = Th.operator( )(K[2]) + 1;
        int lab = K.lab;
        fprintf(popenstream, "%i %i %i %i\n", i0, i1, i2, lab);
      }

      fprintf(popenstream, "Edges\n");
      fprintf(popenstream, "%i\n", nbe);

      for (int k = 0; k < nbe; k++) {
        const Mesh::BorderElement &K(Th.be(k));
        int i0 = Th.operator( )(K[0]) + 1;
        int i1 = Th.operator( )(K[1]) + 1;
        int lab = K.lab;
        fprintf(popenstream, "%i %i %i\n", i0, i1, lab);
      }

      fprintf(popenstream, "End\n");
    }

    // solution with a mesh
    if (nbsol > 0) {
      if (filebin) {
        int cod = 1;
        int NulPos = 0;
        int KwdCod;
        int codtypjm = 1;
        // determination of number solutions associated with a mesh
        fwrite((unsigned char *)&cod, WrdSiz, 1, popenstream);
        fwrite((unsigned char *)&ver, WrdSiz, 1, popenstream);
        KwdCod = GmfDimension;
        fwrite((unsigned char *)&KwdCod, WrdSiz, 1, popenstream);
        fwrite((unsigned char *)&NulPos, WrdSiz, 1, popenstream);
        fwrite((unsigned char *)&dimp, WrdSiz, 1, popenstream);

        if (order == 0) {
          KwdCod = GmfSolAtTriangles;
          fwrite((unsigned char *)&KwdCod, WrdSiz, 1, popenstream);
          fwrite((unsigned char *)&NulPos, WrdSiz, 1, popenstream);
          fwrite((unsigned char *)&nt, WrdSiz, 1, popenstream);

          printf("SolAtTriangles nt=%i\n", nt);
        }

        if (order == 1) {
          KwdCod = GmfSolAtVertices;
          fwrite((unsigned char *)&KwdCod, WrdSiz, 1, popenstream);
          fwrite((unsigned char *)&NulPos, WrdSiz, 1, popenstream);
          fwrite((unsigned char *)&nv, WrdSiz, 1, popenstream);

          printf("SolAtVertices nv=%i\n", nv);
        }

        typsol = TypTab[jojo];

        fwrite((unsigned char *)&codtypjm, WrdSiz, 1, popenstream);
        fwrite((unsigned char *)&typsol, WrdSiz, 1, popenstream);
      } else {
        fprintf(popenstream, "MeshVersionFormatted %i\n", ver);
        fprintf(popenstream, "Dimension %i\n", dimp);
        if (order == 0) {
          fprintf(popenstream, "SolAtTriangles\n");
          fprintf(popenstream, "%i\n", nt);
        }

        if (order == 1) {
          fprintf(popenstream, "SolAtVertices\n");
          fprintf(popenstream, "%i\n", nv);
        }

        typsol = TypTab[jojo];

        fprintf(popenstream, "%i %i\n", 1, typsol);
      }

      if (typsol == 1) {
        if (order == 0) {
          vxx = 0.;
          MeshPoint *mp3(MeshPointStack(stack));
          R2 Cdg_hat = R2(1. / 3., 1. / 3.);

          for (int it = 0; it < Th.nt; ++it) {
            jojo1 = jojo + 1 + offset * numTht[it];
            const Mesh::Triangle &K(Th.t(it));
            mp3->set(Th, K(Cdg_hat), Cdg_hat, K, K.lab);
            vxx[it] = l[jojo1].eval(0, stack);    // GetAny< double >( (*nargs[1])(stack) );
          }

          if (filebin) {
            for (int k = 0; k < nt; k++) {
              fwrite((unsigned char *)&(vxx[k]), WrdSiz, 2, popenstream);
            }
          } else {
            for (int k = 0; k < nt; k++) {
              fprintf(popenstream, "%f\n", vxx[k]);
            }
          }
        } else if (order == 1) {
          vxx = 0.;
          KN< int > takemesh(nv);
          MeshPoint *mp3(MeshPointStack(stack));

          takemesh = 0;

          for (int it = 0; it < Th.nt; ++it) {
            jojo1 = jojo + 1 + offset * numTht[it];

            for (int iv = 0; iv < 3; ++iv) {
              int i = Th(it, iv);

              mp3->setP(&Th, it, iv);
              vxx[i] =
                vxx[i] + l[jojo1].eval(0, stack);    // GetAny< double >( (*nargs[1])(stack) );
              takemesh[i] = takemesh[i] + 1;
            }
          }

          if (filebin) {
            for (int k = 0; k < nv; k++) {
              vxx[k] = vxx[k] / takemesh[k];
              fwrite((unsigned char *)&(vxx[k]), WrdSiz, 2, popenstream);
            }
          } else {
            for (int k = 0; k < nv; k++) {
              vxx[k] = vxx[k] / takemesh[k];
              fprintf(popenstream, "%f\n", vxx[k]);
            }
          }
        }
      } else if (typsol == 2) {
        if (order == 0) {
          // KN<double>  vxx(nt),vyy(nt);
          MeshPoint *mp3(MeshPointStack(stack));
          R2 Cdg_hat = R2(1. / 3., 1. / 3.);

          vxx = 0.;
          vyy = 0.;

          for (int it = 0; it < Th.nt; ++it) {
            jojo1 = jojo + 1 + offset * numTht[it];
            const Mesh::Triangle &K(Th.t(it));
            mp3->set(Th, K(Cdg_hat), Cdg_hat, K, K.lab);
            vxx[it] = l[jojo1].eval(0, stack);    // GetAny< double >( (*xx)(stack) );
            vyy[it] = l[jojo1].eval(1, stack);    // GetAny< double >( (*yy)(stack) );
          }

          if (filebin) {
            for (int k = 0; k < nt; k++) {
              fwrite((unsigned char *)&(vxx[k]), WrdSiz, 2, popenstream);
              fwrite((unsigned char *)&(vyy[k]), WrdSiz, 2, popenstream);
            }
          } else {
            for (int k = 0; k < nt; k++) {
              fprintf(popenstream, "%f %f\n", vxx[k], vyy[k]);
            }
          }
        } else if (order == 1) {
          KN< int > takemesh(nv);
          MeshPoint *mp3(MeshPointStack(stack));

          takemesh = 0;
          vxx = 0.;
          vyy = 0.;

          for (int it = 0; it < Th.nt; ++it) {
            jojo1 = jojo + 1 + offset * numTht[it];

            for (int iv = 0; iv < 3; ++iv) {
              int i = Th(it, iv);

              mp3->setP(&Th, it, iv);
              vxx[i] = vxx[i] + l[jojo1].eval(0, stack);    // GetAny< double >( (*xx)(stack) );
              vyy[i] = vyy[i] + l[jojo1].eval(1, stack);    // GetAny< double >( (*yy)(stack) );

              takemesh[i] = takemesh[i] + 1;
            }
          }

          if (filebin) {
            for (int k = 0; k < nv; k++) {
              vxx[k] = vxx[k] / takemesh[k];
              vyy[k] = vyy[k] / takemesh[k];
              fwrite((unsigned char *)&(vxx[k]), WrdSiz, 2, popenstream);
              fwrite((unsigned char *)&(vyy[k]), WrdSiz, 2, popenstream);
            }
          } else {
            for (int k = 0; k < nv; k++) {
              vxx[k] = vxx[k] / takemesh[k];
              vyy[k] = vyy[k] / takemesh[k];
              fprintf(popenstream, "%f %f\n", vxx[k], vyy[k]);
            }
          }
        }
      } else if (typsol == 3) {
        if (order == 0) {
          vxx = 0.;
          vyx = 0.;
          vyy = 0.;
          MeshPoint *mp3(MeshPointStack(stack));
          R2 Cdg_hat = R2(1. / 3., 1. / 3.);

          for (int it = 0; it < Th.nt; ++it) {
            jojo1 = jojo + 1 + offset * numTht[it];
            const Mesh::Triangle &K(Th.t(it));
            mp3->set(Th, K(Cdg_hat), Cdg_hat, K, K.lab);
            vxx[it] = l[jojo1].eval(0, stack);    // GetAny< double >( (*tsxx)(stack) );
            vyx[it] = l[jojo1].eval(1, stack);    // GetAny< double >( (*tsyx)(stack) );
            vyy[it] = l[jojo1].eval(2, stack);    // GetAny< double >( (*tsyy)(stack) );
          }

          if (filebin) {
            for (int k = 0; k < nt; k++) {
              fwrite((unsigned char *)&(vxx[k]), WrdSiz, 2, popenstream);
              fwrite((unsigned char *)&(vyx[k]), WrdSiz, 2, popenstream);
              fwrite((unsigned char *)&(vyy[k]), WrdSiz, 2, popenstream);
            }
          } else {
            for (int k = 0; k < nt; k++) {
              fprintf(popenstream, "%f %f %f\n", vxx[k], vyx[k], vyy[k]);
            }
          }
        } else if (order == 1) {
          KN< int > takemesh(nv);
          MeshPoint *mp3(MeshPointStack(stack));

          takemesh = 0;
          vxx = 0.;
          vyx = 0.;
          vyy = 0.;

          for (int it = 0; it < Th.nt; ++it) {
            jojo1 = jojo + 1 + offset * numTht[it];

            for (int iv = 0; iv < 3; ++iv) {
              int i = Th(it, iv);

              mp3->setP(&Th, it, iv);
              vxx[i] = vxx[i] + l[jojo1].eval(0, stack);    // GetAny< double >( (*tsxx)(stack) );
              vyx[i] = vyx[i] + l[jojo1].eval(0, stack);    // GetAny< double >( (*tsyx)(stack) );
              vyy[i] = vyy[i] + l[jojo1].eval(0, stack);    // GetAny< double >( (*tsyy)(stack) );

              takemesh[i] = takemesh[i] + 1;
            }
          }

          if (filebin) {
            for (int k = 0; k < nv; k++) {
              vxx[k] = vxx[k] / takemesh[k];
              vyx[k] = vyx[k] / takemesh[k];
              vyy[k] = vyy[k] / takemesh[k];

              fwrite((unsigned char *)&(vxx[k]), WrdSiz, 2, popenstream);
              fwrite((unsigned char *)&(vyx[k]), WrdSiz, 2, popenstream);
              fwrite((unsigned char *)&(vyy[k]), WrdSiz, 2, popenstream);
            }
          } else {
            for (int k = 0; k < nv; k++) {
              vxx[k] = vxx[k] / takemesh[k];
              vyx[k] = vyx[k] / takemesh[k];
              vyy[k] = vyy[k] / takemesh[k];
              fprintf(popenstream, "%f %f %f\n", vxx[k], vyx[k], vyy[k]);
            }
          }
        }
      }

      if (filebin) {
        int NulPos = 0;
        int KwdCod = GmfEnd;
        fwrite((unsigned char *)&KwdCod, WrdSiz, 1, popenstream);
        fwrite((unsigned char *)&NulPos, WrdSiz, 1, popenstream);
      } else {
        fprintf(popenstream, "End\n");
      }

      if (boolsave) {
        if (verbosity > 2) {
          cout << "writing solution in file" << endl;
        }

        if (typsol == 1) {
          writetabsol(datasize, nboftmp, vxx, solsave);
          nboftmp = nboftmp + 1;
        } else if (typsol == 2) {
          writetabsol(datasize, nboftmp, vxx, vyy, solsave);
          nboftmp = nboftmp + 2;
        } else if (typsol == 3) {
          writetabsol(datasize, nboftmp, vxx, vyx, vyy, solsave);
          nboftmp = nboftmp + 3;
        }
      }
    }
  }

  // drawing part  ------------------------------

  // suppression du no wait car plante sur mac
  // il faut faire thread qui lance ca de maniere asyncrone ...
  // a faire .. FH ...
  // if (wait && !NoWait)
  { pclose(popenstream); }

  // rajout pout la sauvegarde de la solution
  if (boolsave) {
      int64_t outm;
    int nbtype = nbsol;
    float *OutSolTab = new float[solnbfloat];

    if (!(outm = GmfOpenMesh(saveff->c_str( ), GmfWrite, ver, 2))) {
      cerr << "  -- Mesh3::Save  UNABLE TO OPEN  :" << saveff << endl;
      exit(1);
    }

    if (order == 0) {
      GmfSetKwd(outm, GmfSolAtTriangles, datasize, nbtype, TypTab);

      for (int k = 0; k < datasize; k++) {
        for (int i = 0; i < solnbfloat; i++) {
          OutSolTab[i] = solsave(i, k);
        }

        GmfSetLin(outm, GmfSolAtTriangles, OutSolTab);
      }
    } else if (order == 1) {
      GmfSetKwd(outm, GmfSolAtVertices, datasize, nbtype, TypTab);

      for (int k = 0; k < datasize; k++) {
        for (int i = 0; i < solnbfloat; i++) {
          OutSolTab[i] = solsave(i, k);
        }

        GmfSetLin(outm, GmfSolAtVertices, OutSolTab);
      }
    }

    GmfCloseMesh(outm);

    delete[] OutSolTab;
  }

  delete[] commandline;
  delete pTh;

  return valsortie;
}

template< class v_fes >
class PopenMeditMesh3_Op : public E_F0mps {
 public:
  typedef long Result;
  Expression eTh;
  Expression filename;
  long offset;
  long nbTh;
  struct Expression2 {
    long what;       // 0 mesh, 1 scalar, 2 vector, 3 symtensor
    long nbfloat;    // 0 mesh(3D), 1 scalar, 2 vector (3D), 3 symtensor(3D)
    Expression e[6];
    Expression2( ) {
      e[0] = 0;
      e[1] = 0;
      e[2] = 0;
      e[3] = 0;
      e[4] = 0;
      e[5] = 0;
      what = 0;
      nbfloat = 0;
    };
    Expression &operator[](int i) { return e[i]; }

    double eval(int i, Stack stack) const {
      if (e[i]) {
        return GetAny< double >((*e[i])(stack));
      } else {
        return 0;
      }
    }

    const Mesh3 &evalm(int i, Stack stack) const {
      throwassert(e[i]);
      return *GetAny< pmesh3 >((*e[i])(stack));
    }
  };
  vector< Expression2 > l;

  static const int n_name_param = 5;
  static basicAC_F0::name_and_type name_param[];
  Expression nargs[n_name_param];
  long arg(int i, Stack stack, long a) const {
    return nargs[i] ? GetAny< long >((*nargs[i])(stack)) : a;
  }

  string *arg(int i, Stack stack, string *a) const {
    return nargs[i] ? GetAny< string * >((*nargs[i])(stack)) : a;
  }

 public:
  PopenMeditMesh3_Op(const basicAC_F0 &args) : l(args.size( ) - 1) {
    int nbofsol;
    int ddim = 3;
    int stsize = 6;

    args.SetNameParam(n_name_param, name_param, nargs);

    if (BCastTo< string * >(args[0])) {
      filename = CastTo< string * >(args[0]);
    }

    for (size_t i = 1; i < args.size( ); i++) {
      size_t jj = i - 1;
      if (BCastTo< double >(args[i])) {
        l[jj].what = 1;
        l[jj].nbfloat = 1;
        l[jj][0] = to< double >(args[i]);
      } else if (args[i].left( ) == atype< E_Array >( )) {
        const E_Array *a0 = dynamic_cast< const E_Array * >(args[i].LeftValue( ));
        if (a0 == nullptr) cout << "Dynamic cast error" << endl;
        if (a0->size( ) != ddim && a0->size( ) != stsize) {
          CompileError(
            "medit in 3D: vector solution is 3 composant, tensor solution is 6 composant");
        }

        if (a0->size( ) == ddim) {
          // vector solution
          l[jj].what = 2;
          l[jj].nbfloat = ddim;

          for (int j = 0; j < ddim; j++) {
            l[jj][j] = to< double >((*a0)[j]);
          }
        } else if (a0->size( ) == stsize) {
          // symmetric tensor solution
          l[jj].what = 3;
          l[jj].nbfloat = stsize;

          for (int j = 0; j < stsize; j++) {
            l[jj][j] = to< double >((*a0)[j]);
          }
        }
      } else if (BCastTo< pmesh3 >(args[i])) {
        l[jj].what = 0;
        l[jj].nbfloat = 0;
        l[jj][0] = CastTo< pmesh3 >(args[i]);
      } else {
        CompileError("medit 3d: Sorry no way to save this kind of data");
      }
    }

    // determination of the number of solutions.
    size_t lastTh = 0;
    long offset1;
    offset = 0;
    nbTh = 0;

    for (size_t jj = 0; jj < l.size( ); jj++) {
      if (l[jj].what == 0) {
        nbTh++;
        offset1 = jj - lastTh;
        if (offset == 0) {
          offset = offset1;
        } else if (offset != offset1) {
          CompileError("the number of solution by mesh is different");
        }
        lastTh = jj;
      }
    }

    if (offset == 0) {
      offset = l.size( );
    }
  }

  static ArrayOfaType typeargs( ) {
    return ArrayOfaType(atype< string * >( ), atype< pmesh3 >( ), true);
  }    // all type

  static E_F0 *f(const basicAC_F0 &args) { return new PopenMeditMesh3_Op(args); }

  AnyType operator( )(Stack stack) const;
};

template< class v_fes >
basicAC_F0::name_and_type PopenMeditMesh3_Op< v_fes >::name_param[] = {
  {"order", &typeid(long)},
  {"meditff", &typeid(string *)},
  {"save", &typeid(string *)},
  {"wait", &typeid(bool)},
  {"bin", &typeid(long)}};

template< class v_fes >
AnyType PopenMeditMesh3_Op< v_fes >::operator( )(Stack stack) const {
  if (NoGraphicWindow) {
    return Nothing;
  }

  MeshPoint *mp(MeshPointStack(stack));
  long order(arg(0, stack, 1L));
  //
  int ver = GmfFloat;
  int dimp = 3;
  float fx, fy, fz;
  //
  long valsortie = 0;
  int typsol, nbsol;
  nbsol = offset - 1;

  int TypTab[l.size( ) - 1];

  for (size_t i = 0; i < l.size( ) - 1; i++) {
    TypTab[i] = l[i + 1].what;
  }

  string *ffname = GetAny< string * >((*filename)(stack));
  string *meditff(arg(1, stack, &stringffmedit));
  long filebin(arg(4, stack, 1L));
  int smedit = max(1, nbsol);
  char *commandline = meditcmd(filebin, nbsol, smedit, *meditff, *ffname);

  printf("version de medit %s\n", commandline);
  if (verbosity > 2) {
    cout << "number of solution = " << offset - 1 << endl;
  }

  if (verbosity > 2) {
    cout << "number of mesh     = " << nbTh << endl;
  }

  // lecture des differents maillages
  int nv = 0, nt = 0, nbe = 0;    // sommet, triangles, arretes du maillage unifies

  for (size_t i = 0; i < l.size( ); i = i + offset) {
    if (l[i].what != 0) {
      cerr << "this element is not a mesh" << i << endl;
    }

    const Mesh3 &Thtmp = l[i].evalm(0, stack);
    nv += Thtmp.nv;
    nt += Thtmp.nt;
    nbe += Thtmp.nbe;
  }

  Vertex3 *v = new Vertex3[nv];
  Tet *t = 0;
  if (nt != 0) {
    t = new Tet[nt];
  }

  Triangle3 *b = new Triangle3[nbe];
  Tet *tt = t;
  Triangle3 *bb = b;
  int iv = 0, it = 0, ibe = 0;
  int *numTht = new int[nt];    // numero of Th assoctiated with a tetrahedra
  int jt = 0;

  for (size_t i = 0; i < l.size( ); i = i + offset) {
    int nvtmp = iv;
    const Mesh3 &Thtmp = l[i].evalm(0, stack);

    for (int ii = 0; ii < Thtmp.nv; ii++) {
      const Vertex3 &vi(Thtmp.vertices[ii]);
      v[iv].x = vi.x;
      v[iv].y = vi.y;
      v[iv].z = vi.z;
      v[iv].lab = vi.lab;
      iv++;
    }

    for (int ii = 0; ii < Thtmp.nt; ii++) {
      const Tet &vi(Thtmp.elements[ii]);
      int iv[4];
      iv[0] = nvtmp + Thtmp.operator( )(vi[0]);
      iv[1] = nvtmp + Thtmp.operator( )(vi[1]);
      iv[2] = nvtmp + Thtmp.operator( )(vi[2]);
      iv[3] = nvtmp + Thtmp.operator( )(vi[3]);
      (*tt++).set(v, iv, vi.lab);
      numTht[it] = jt;
      it++;
    }

    for (int ii = 0; ii < Thtmp.nbe; ii++) {
      const Triangle3 &vi(Thtmp.be(ii));
      int iv[3];
      iv[0] = nvtmp + Thtmp.operator( )(vi[0]);
      iv[1] = nvtmp + Thtmp.operator( )(vi[1]);
      iv[2] = nvtmp + Thtmp.operator( )(vi[2]);
      (*bb++).set(v, iv, vi.lab);
      ibe++;
    }

    jt++;
  }

  assert(it == nt);
  assert(iv == nv);
  assert(ibe = nbe);
  if (verbosity > 2) {
    cout << "meditff :: Value of elements: vertex " << nv << " Tet " << nt << " triangle " << nbe
         << endl;
  }

  Mesh3 *pTh = new Mesh3(nv, nt, nbe, v, t, b);
  Mesh3 &Th = *pTh;

  // determination of the number of elements to represent the solution
  int datasize;
  if (order == 0) {
    datasize = nt;
  }

  if (order == 1) {
    datasize = nv;
  }

  // cas de sauvegarde
  bool boolsave = false;
  int solnbfloat = 0;
  KNM< double > solsave(1, 1);
  string *saveff;
  KN< double > vxx = 0., vyx = 0., vyy = 0., vzx = 0., vzy = 0., vzz = 0.;

  if (nbsol > 0) {
    vxx.init(datasize);
    vyx.init(datasize);
    vyy.init(datasize);
    vzx.init(datasize);
    vzy.init(datasize);
    vzz.init(datasize);

    if (nargs[2]) {
      boolsave = true;
      saveff = GetAny< string * >((*nargs[2])(stack));

      for (size_t i = 0; i < offset; i++) {
        solnbfloat = solnbfloat + l[i].nbfloat;
      }

      solsave.init(solnbfloat, datasize);
      solsave = 0.;
    }
  }

  int nboftmp = 0;
  FILE *popenstream = popen(commandline, MODE_WRITE_BINARY);
  if (!popenstream) {
    cerr << " Error popen : " << commandline << endl;
    exit(1);
  }

  // mesh
  int jojo1;

  for (int jojo = 0; jojo < smedit; jojo++) {
    if (filebin) {
      int cod = 1;
      int KwdCod;
      int NulPos = 0;

      // determination of number solutions associated with a mesh
      fwrite((unsigned char *)&cod, WrdSiz, 1, popenstream);
      fwrite((unsigned char *)&ver, WrdSiz, 1, popenstream);
      KwdCod = GmfDimension;
      fwrite((unsigned char *)&KwdCod, WrdSiz, 1, popenstream);
      fwrite((unsigned char *)&NulPos, WrdSiz, 1, popenstream);
      fwrite((unsigned char *)&dimp, WrdSiz, 1, popenstream);

      // vertex
      KwdCod = GmfVertices;
      fwrite((unsigned char *)&KwdCod, WrdSiz, 1, popenstream);
      fwrite((unsigned char *)&NulPos, WrdSiz, 1, popenstream);
      fwrite((unsigned char *)&nv, WrdSiz, 1, popenstream);

      for (int k = 0; k < nv; k++) {
        const Vertex3 &P = Th.vertices[k];
        fx = P.x;
        fy = P.y;
        fz = P.z;
        fwrite((unsigned char *)&fx, WrdSiz, 1, popenstream);
        fwrite((unsigned char *)&fy, WrdSiz, 1, popenstream);
        fwrite((unsigned char *)&fz, WrdSiz, 1, popenstream);
        fwrite((unsigned char *)&(P.lab), WrdSiz, 1, popenstream);
      }

      // tetrahedra
      KwdCod = GmfTetrahedra;
      fwrite((unsigned char *)&KwdCod, WrdSiz, 1, popenstream);
      fwrite((unsigned char *)&NulPos, WrdSiz, 1, popenstream);
      fwrite((unsigned char *)&nt, WrdSiz, 1, popenstream);

      for (int k = 0; k < nt; k++) {
        const Tet &K(Th.elements[k]);
        int i0 = Th.operator( )(K[0]) + 1;
        int i1 = Th.operator( )(K[1]) + 1;
        int i2 = Th.operator( )(K[2]) + 1;
        int i3 = Th.operator( )(K[3]) + 1;
        int lab = K.lab;
        fwrite((unsigned char *)&i0, WrdSiz, 1, popenstream);
        fwrite((unsigned char *)&i1, WrdSiz, 1, popenstream);
        fwrite((unsigned char *)&i2, WrdSiz, 1, popenstream);
        fwrite((unsigned char *)&i3, WrdSiz, 1, popenstream);
        fwrite((unsigned char *)&lab, WrdSiz, 1, popenstream);
      }

      // triangles
      KwdCod = GmfTriangles;
      fwrite((unsigned char *)&KwdCod, WrdSiz, 1, popenstream);
      fwrite((unsigned char *)&NulPos, WrdSiz, 1, popenstream);
      fwrite((unsigned char *)&nbe, WrdSiz, 1, popenstream);

      for (int k = 0; k < nbe; k++) {
        const Triangle3 &K(Th.be(k));
        int i0 = Th.operator( )(K[0]) + 1;
        int i1 = Th.operator( )(K[1]) + 1;
        int i2 = Th.operator( )(K[2]) + 1;
        int lab = K.lab;
        fwrite((unsigned char *)&i0, WrdSiz, 1, popenstream);
        fwrite((unsigned char *)&i1, WrdSiz, 1, popenstream);
        fwrite((unsigned char *)&i2, WrdSiz, 1, popenstream);
        fwrite((unsigned char *)&lab, WrdSiz, 1, popenstream);
      }

      // End
      KwdCod = GmfEnd;
      fwrite((unsigned char *)&KwdCod, WrdSiz, 1, popenstream);
      fwrite((unsigned char *)&NulPos, WrdSiz, 1, popenstream);
    } else {
      fprintf(popenstream, "MeshVersionFormatted\n");
      fprintf(popenstream, "%i\n", ver);
      fprintf(popenstream, "Dimension\n");
      fprintf(popenstream, "%i\n", dimp);
      fprintf(popenstream, "Vertices\n");
      fprintf(popenstream, "%i\n", nv);

      for (int k = 0; k < nv; k++) {
        const Vertex3 &P = Th.vertices[k];
        fx = P.x;
        fy = P.y;
        fz = P.z;
        fprintf(popenstream, "%f %f %f %i\n", fx, fy, fz, P.lab);
      }

      fprintf(popenstream, "Tetrahedra\n");
      fprintf(popenstream, "%i\n", nt);

      for (int k = 0; k < nt; k++) {
        const Tet &K(Th.elements[k]);
        int i0 = Th.operator( )(K[0]) + 1;
        int i1 = Th.operator( )(K[1]) + 1;
        int i2 = Th.operator( )(K[2]) + 1;
        int i3 = Th.operator( )(K[3]) + 1;
        int lab = K.lab;
        fprintf(popenstream, "%i %i %i %i %i\n", i0, i1, i2, i3, lab);
      }

      fprintf(popenstream, "Triangles\n");
      fprintf(popenstream, "%i\n", nbe);

      for (int k = 0; k < nbe; k++) {
        const Triangle3 &K(Th.be(k));
        int i0 = Th.operator( )(K[0]) + 1;
        int i1 = Th.operator( )(K[1]) + 1;
        int i2 = Th.operator( )(K[2]) + 1;
        int lab = K.lab;
        fprintf(popenstream, "%i %i %i %i\n", i0, i1, i2, lab);
      }

      fprintf(popenstream, "End");
    }

    if (nbsol > 0) {
      if (filebin) {
        int cod = 1;
        int NulPos = 0;
        int KwdCod;
        fwrite((unsigned char *)&cod, WrdSiz, 1, popenstream);
        fwrite((unsigned char *)&ver, WrdSiz, 1, popenstream);
        KwdCod = GmfDimension;
        fwrite((unsigned char *)&KwdCod, WrdSiz, 1, popenstream);
        fwrite((unsigned char *)&NulPos, WrdSiz, 1, popenstream);
        fwrite((unsigned char *)&dimp, WrdSiz, 1, popenstream);
      } else {
        fprintf(popenstream, "MeshVersionFormatted %i\n", ver);
        fprintf(popenstream, "Dimension %i\n", dimp);

        // detemination of solution
        // default scalaire // faire tableau pour plusieurs
      }

      typsol = TypTab[jojo];

      if (order == 0) {
        if (filebin) {
          int NulPos = 0;
          int KwdCod = GmfSolAtTetrahedra;
          int codtypjm = 1;
          fwrite((unsigned char *)&KwdCod, WrdSiz, 1, popenstream);
          fwrite((unsigned char *)&NulPos, WrdSiz, 1, popenstream);
          fwrite((unsigned char *)&datasize, WrdSiz, 1, popenstream);
          fwrite((unsigned char *)&codtypjm, WrdSiz, 1, popenstream);
          fwrite((unsigned char *)&typsol, WrdSiz, 1, popenstream);
        } else {
          fprintf(popenstream, "SolAtTetrahedra\n");
          fprintf(popenstream, "%i\n", datasize);
          fprintf(popenstream, "%i %i\n", 1, typsol);
        }

        if (typsol == 1) {
          MeshPoint *mp3(MeshPointStack(stack));
          R3 Cdg_hat = R3(1. / 4., 1. / 4., 1. / 4.);

          vxx = 0.;

          for (int it = 0; it < nt; it++) {
            jojo1 = jojo + 1 + offset * numTht[it];
            const Tet &K(Th.elements[it]);
            mp3->set(Th, K(Cdg_hat), Cdg_hat, K, K.lab);

            vxx[it] = l[jojo1].eval(0, stack);    // GetAny< double >( (*nargs[1])(stack) );
          }

          if (filebin) {
            for (int k = 0; k < nt; k++) {
              fwrite((unsigned char *)&(vxx[k]), WrdSiz, 2, popenstream);
            }
          } else {
            for (int k = 0; k < nt; k++) {
              fprintf(popenstream, "%f\n", vxx[k]);
            }
          }
        } else if (typsol == 2) {
          vyy = 0.;
          vxx = 0.;
          vzz = 0.;
          MeshPoint *mp3(MeshPointStack(stack));
          R3 Cdg_hat = R3(1. / 4., 1. / 4., 1. / 4.);

          for (int it = 0; it < Th.nt; ++it) {
            jojo1 = jojo + 1 + offset * numTht[it];
            const Tet &K(Th.elements[it]);
            mp3->set(Th, K(Cdg_hat), Cdg_hat, K, K.lab);

            vxx[it] = l[jojo1].eval(0, stack);    // GetAny< double >( (*xx)(stack) );
            vyy[it] = l[jojo1].eval(1, stack);    // GetAny< double >( (*yy)(stack) );
            vzz[it] = l[jojo1].eval(2, stack);    // GetAny< double >( (*zz)(stack) );
          }

          if (filebin) {
            for (int k = 0; k < nt; k++) {
              fwrite((unsigned char *)&(vxx[k]), WrdSiz, 2, popenstream);
              fwrite((unsigned char *)&(vyy[k]), WrdSiz, 2, popenstream);
              fwrite((unsigned char *)&(vzz[k]), WrdSiz, 2, popenstream);
            }
          } else {
            for (int k = 0; k < nt; k++) {
              fprintf(popenstream, "%f %f %f\n", vxx[k], vyy[k], vzz[k]);
            }
          }
        } else if (typsol == 3) {
          vxx(nv) = 0.;
          vyx(nv) = 0.;
          vyy(nv) = 0.;
          vzx(nv) = 0.;
          vzy(nv) = 0.;
          vzz(nv) = 0.;
          MeshPoint *mp3(MeshPointStack(stack));
          R3 Cdg_hat = R3(1. / 4., 1. / 4., 1. / 4.);

          for (int it = 0; it < Th.nt; ++it) {
            jojo1 = jojo + 1 + offset * numTht[it];
            const Tet &K(Th.elements[it]);
            mp3->set(Th, K(Cdg_hat), Cdg_hat, K, K.lab);

            vxx[it] = l[jojo1].eval(0, stack);    // GetAny< double >( (*tsxx)(stack) );
            vyx[it] = l[jojo1].eval(1, stack);    // GetAny< double >( (*tsyx)(stack) );
            vyy[it] = l[jojo1].eval(2, stack);    // GetAny< double >( (*tsyy)(stack) );
            vzx[it] = l[jojo1].eval(3, stack);    // GetAny< double >( (*tszx)(stack) );
            vzy[it] = l[jojo1].eval(4, stack);    // GetAny< double >( (*tszy)(stack) );
            vzz[it] = l[jojo1].eval(5, stack);    // GetAny< double >( (*tszz)(stack) );
          }

          if (filebin) {
            for (int k = 0; k < nt; k++) {
              fwrite((unsigned char *)&(vxx[k]), WrdSiz, 2, popenstream);
              fwrite((unsigned char *)&(vyx[k]), WrdSiz, 2, popenstream);
              fwrite((unsigned char *)&(vyy[k]), WrdSiz, 2, popenstream);
              fwrite((unsigned char *)&(vzx[k]), WrdSiz, 2, popenstream);
              fwrite((unsigned char *)&(vzy[k]), WrdSiz, 2, popenstream);
              fwrite((unsigned char *)&(vzz[k]), WrdSiz, 2, popenstream);
            }
          } else {
            for (int k = 0; k < nt; k++) {
              fprintf(popenstream, "%f %f %f %f %f %f\n", vxx[k], vyx[k], vyy[k], vzx[k], vzy[k],
                      vzz[k]);
            }
          }
        }
      } else if (order == 1) {
        if (filebin) {
          int NulPos = 0;
          int KwdCod = GmfSolAtVertices;
          int codtypjm = 1;
          fwrite((unsigned char *)&KwdCod, WrdSiz, 1, popenstream);
          fwrite((unsigned char *)&NulPos, WrdSiz, 1, popenstream);
          fwrite((unsigned char *)&datasize, WrdSiz, 1, popenstream);
          fwrite((unsigned char *)&codtypjm, WrdSiz, 1, popenstream);
          fwrite((unsigned char *)&typsol, WrdSiz, 1, popenstream);
        } else {
          fprintf(popenstream, "SolAtVertices\n");
          fprintf(popenstream, "%i\n", datasize);
          fprintf(popenstream, "%i %i\n", 1, typsol);
        }

        if (typsol == 1) {
          KN< int > takemesh(nv);
          MeshPoint *mp3(MeshPointStack(stack));

          takemesh = 0;
          vxx = 0.;

          for (int it = 0; it < Th.nt; it++) {
            jojo1 = jojo + 1 + offset * numTht[it];

            for (int iv = 0; iv < 4; iv++) {
              int i = Th(it, iv);

              mp3->setP(&Th, it, iv);
              vxx[i] =
                vxx[i] + l[jojo1].eval(0, stack);    // GetAny< double >( (*nargs[1])(stack) );

              takemesh[i] = takemesh[i] + 1;
            }
          }

          if (filebin) {
            for (int k = 0; k < nv; k++) {
              vxx[k] = vxx[k] / takemesh[k];
              fwrite((unsigned char *)&(vxx[k]), WrdSiz, 2, popenstream);
            }
          } else {
            for (int k = 0; k < nv; k++) {
              vxx[k] = vxx[k] / takemesh[k];
              fprintf(popenstream, "%f\n", vxx[k]);
            }
          }
        } else if (typsol == 2) {
          KN< int > takemesh(nv);
          MeshPoint *mp3(MeshPointStack(stack));

          vxx = 0.;
          vyy = 0.;
          vzz = 0.;
          takemesh = 0;

          for (int it = 0; it < Th.nt; ++it) {
            jojo1 = jojo + 1 + offset * numTht[it];

            for (int iv = 0; iv < 4; ++iv) {
              int i = Th(it, iv);

              mp3->setP(&Th, it, iv);
              vxx[i] = vxx[i] + l[jojo1].eval(0, stack);    // GetAny< double >( (*xx)(stack) );
              vyy[i] = vyy[i] + l[jojo1].eval(1, stack);    // GetAny< double >( (*yy)(stack) );
              vzz[i] = vzz[i] + l[jojo1].eval(2, stack);    // GetAny< double >( (*zz)(stack) );

              takemesh[i] = takemesh[i] + 1;
            }
          }

          if (filebin) {
            for (int k = 0; k < nv; k++) {
              vxx[k] = vxx[k] / takemesh[k];
              vyy[k] = vyy[k] / takemesh[k];
              vzz[k] = vzz[k] / takemesh[k];
              fwrite((unsigned char *)&(vxx[k]), WrdSiz, 2, popenstream);
              fwrite((unsigned char *)&(vyy[k]), WrdSiz, 2, popenstream);
              fwrite((unsigned char *)&(vzz[k]), WrdSiz, 2, popenstream);
            }
          } else {
            for (int k = 0; k < nv; k++) {
              vxx[k] = vxx[k] / takemesh[k];
              vyy[k] = vyy[k] / takemesh[k];
              vzz[k] = vzz[k] / takemesh[k];
              fprintf(popenstream, "%f %f %f\n", vxx[k], vyy[k], vzz[k]);
            }
          }
        } else if (typsol == 3) {
          KN< int > takemesh(nv);
          MeshPoint *mp3(MeshPointStack(stack));

          vxx = 0.;
          vyx = 0.;
          vyy = 0.;
          vzx = 0.;
          vzy = 0.;
          vzz = 0.;
          takemesh = 0;

          for (int it = 0; it < Th.nt; ++it) {
            jojo1 = jojo + 1 + offset * numTht[it];

            for (int iv = 0; iv < 4; ++iv) {
              int i = Th(it, iv);

              mp3->setP(&Th, it, iv);
              vxx[i] = vxx[i] + l[jojo1].eval(0, stack);    // GetAny< double >( (*tsxx)(stack) );
              vyx[i] = vyx[i] + l[jojo1].eval(1, stack);    // GetAny< double >( (*tsyx)(stack) );
              vyy[i] = vyy[i] + l[jojo1].eval(2, stack);    // GetAny< double >( (*tsyy)(stack) );
              vzx[i] = vzx[i] + l[jojo1].eval(3, stack);    // GetAny< double >( (*tszx)(stack) );
              vzy[i] = vzy[i] + l[jojo1].eval(4, stack);    // GetAny< double >( (*tszy)(stack) );
              vzz[i] = vzz[i] + l[jojo1].eval(5, stack);    // GetAny< double >( (*tszz)(stack) );

              takemesh[i] = takemesh[i] + 1;
            }
          }

          if (filebin) {
            for (int k = 0; k < nv; k++) {
              vxx[k] = vxx[k] / takemesh[k];
              vyx[k] = vyx[k] / takemesh[k];
              vyy[k] = vyy[k] / takemesh[k];
              vzx[k] = vzx[k] / takemesh[k];
              vzy[k] = vzy[k] / takemesh[k];
              vzz[k] = vzz[k] / takemesh[k];

              fwrite((unsigned char *)&(vxx[k]), WrdSiz, 2, popenstream);
              fwrite((unsigned char *)&(vyx[k]), WrdSiz, 2, popenstream);
              fwrite((unsigned char *)&(vyy[k]), WrdSiz, 2, popenstream);
              fwrite((unsigned char *)&(vzx[k]), WrdSiz, 2, popenstream);
              fwrite((unsigned char *)&(vzy[k]), WrdSiz, 2, popenstream);
              fwrite((unsigned char *)&(vzz[k]), WrdSiz, 2, popenstream);
            }
          } else {
            for (int k = 0; k < nv; k++) {
              vxx[k] = vxx[k] / takemesh[k];
              vyx[k] = vyx[k] / takemesh[k];
              vyy[k] = vyy[k] / takemesh[k];
              vzx[k] = vzx[k] / takemesh[k];
              vzy[k] = vzy[k] / takemesh[k];
              vzz[k] = vzz[k] / takemesh[k];

              fprintf(popenstream, "%f %f %f %f %f %f\n", vxx[k], vyx[k], vyy[k], vzx[k], vzy[k],
                      vzz[k]);
            }
          }
        }
      }

      if (filebin) {
        int NulPos = 0;
        int KwdCod = GmfEnd;
        fwrite((unsigned char *)&KwdCod, WrdSiz, 1, popenstream);
        fwrite((unsigned char *)&NulPos, WrdSiz, 1, popenstream);
      } else {
        fprintf(popenstream, "End");
      }

      if (boolsave) {
        if (verbosity > 2) {
          cout << "writing solution in file" << endl;
        }

        if (typsol == 1) {
          writetabsol(datasize, nboftmp, vxx, solsave);
          nboftmp = nboftmp + 1;
        } else if (typsol == 2) {
          writetabsol(datasize, nboftmp, vxx, vyy, vzz, solsave);
          nboftmp = nboftmp + 3;
        } else if (typsol == 3) {
          writetabsol(datasize, nboftmp, vxx, vyx, vyy, vzx, vzy, vzz, solsave);
          nboftmp = nboftmp + 6;
        }
      }
    }
  }

  delete[] numTht;
  // fermeture du stream pour popen

  // drawing part  ------------------------------

  { pclose(popenstream); }

  if (boolsave) {
      int64_t outm;
    int nbtype = nbsol;
    float *OutSolTab = new float[solnbfloat];

    if (!(outm = GmfOpenMesh(saveff->c_str( ), GmfWrite, ver, 3))) {
      cerr << "  -- Mesh3::Save  UNABLE TO OPEN  :" << saveff << endl;
      exit(1);
    }

    if (order == 0) {
      GmfSetKwd(outm, GmfSolAtTetrahedra, datasize, nbtype, TypTab);

      for (int k = 0; k < datasize; k++) {
        for (int i = 0; i < solnbfloat; i++) {
          OutSolTab[i] = solsave(i, k);
        }

        GmfSetLin(outm, GmfSolAtTetrahedra, OutSolTab);
      }
    } else if (order == 1) {
      GmfSetKwd(outm, GmfSolAtVertices, datasize, nbtype, TypTab);

      for (int k = 0; k < datasize; k++) {
        for (int i = 0; i < solnbfloat; i++) {
          OutSolTab[i] = solsave(i, k);
        }

        GmfSetLin(outm, GmfSolAtVertices, OutSolTab);
      }
    }

    delete[] OutSolTab;
    GmfCloseMesh(outm);
  }

  delete[] commandline;
  delete pTh;

  return valsortie;
}

template< class MMesh, class v_fes >
class PopenMeditMeshT_Op : public E_F0mps {
 public:
  typedef long Result;
  typedef const MMesh *ppmesh;
  Expression eTh;
  Expression filename;
  long offset;
  long nbTh;
  struct Expression2 {
    long what;       // 0 mesh, 1 scalar, 2 vector, 3 symtensor
    long nbfloat;    // 0 mesh(3D), 1 scalar, 2 vector (3D), 3 symtensor(3D)
    Expression e[6];
    Expression2( ) { e[0]=0 ; e[1]=0 ; e[2]=0; e[3]=0 ; e[4]=0; e[5]=0 ; what=0; nbfloat=0;
    };
    Expression &operator[](int i) { return e[i]; }

    double eval(int i, Stack stack) const {
      if (e[i]) {
        return GetAny< double >((*e[i])(stack));
      } else {
        return 0;
      }
    }

    const MMesh &evalm(int i, Stack stack) const {
      throwassert(e[i]);
      return *GetAny< ppmesh >((*e[i])(stack));
    }
  };
  vector< Expression2 > l;

  static const int n_name_param = 5;
  static basicAC_F0::name_and_type name_param[];
  Expression nargs[n_name_param];
  long arg(int i, Stack stack, long a) const {
    return nargs[i] ? GetAny< long >((*nargs[i])(stack)) : a;
  }

  string *arg(int i, Stack stack, string *a) const {
    return nargs[i] ? GetAny< string * >((*nargs[i])(stack)) : a;
  }

 public:
  PopenMeditMeshT_Op(const basicAC_F0 &args) : l(args.size( ) - 1) {
    int nbofsol;
    int ddim = 3;
    int stsize = 3;

    args.SetNameParam(n_name_param, name_param, nargs);

    if (BCastTo< string * >(args[0])) {
      filename = CastTo< string * >(args[0]);
    }

    for (size_t i = 1; i < args.size( ); i++) {
      size_t jj = i - 1;
      if (BCastTo< double >(args[i])) {
        l[jj].what = 1;
        l[jj].nbfloat = 1;
        l[jj][0] = to< double >(args[i]);
      } else if (args[i].left( ) == atype< E_Array >( )) {
        const E_Array *a0 = dynamic_cast< const E_Array * >(args[i].LeftValue( ));
        if (a0 == nullptr) cout << "Dynamic cast error" << endl;
        if (a0->size( ) != ddim && a0->size( ) != stsize) {
          CompileError(
            "medit in 3D: vector solution is 3 composant, tensor solution is 6 composant");
        }

        if (a0->size( ) == ddim) {
          // vector solution
          l[jj].what = 2;
          l[jj].nbfloat = ddim;

          for (int j = 0; j < ddim; j++) {
            l[jj][j] = to< double >((*a0)[j]);
          }
        } else if (a0->size( ) == stsize) {
          // symmetric tensor solution
          l[jj].what = 3;
          l[jj].nbfloat = stsize;

          for (int j = 0; j < stsize; j++) {
            l[jj][j] = to< double >((*a0)[j]);
          }
        }
      } else if (BCastTo< ppmesh >(args[i])) {
        l[jj].what = 0;
        l[jj].nbfloat = 0;
        l[jj][0] = CastTo< ppmesh >(args[i]);
      } else {
        CompileError("medit 3d: Sorry no way to save this kind of data");
      }
    }

    // determination of the number of solutions.
    size_t lastTh = 0;
    long offset1;
    offset = 0;
    nbTh = 0;

    for (size_t jj = 0; jj < l.size( ); jj++) {
      if (l[jj].what == 0) {
        nbTh++;
        offset1 = jj - lastTh;
        if (offset == 0) {
          offset = offset1;
        } else if (offset != offset1) {
          CompileError("the number of solution by mesh is different");
        }
        lastTh = jj;
      }
    }

    if (offset == 0) {
      offset = l.size( );
    }
  }

  static ArrayOfaType typeargs( ) {
    return ArrayOfaType(atype< string * >( ), atype< ppmesh >( ), true);
  }    // all type

  static E_F0 *f(const basicAC_F0 &args) { return new PopenMeditMeshT_Op< MMesh, v_fes >(args); }

  AnyType operator( )(Stack stack) const;
};

template<>
basicAC_F0::name_and_type PopenMeditMeshT_Op< MeshS, v_fesS >::name_param[] = {
  {"order", &typeid(long)},
  {"meditff", &typeid(string *)},
  {"save", &typeid(string *)},
  {"wait", &typeid(bool)},
  {"bin", &typeid(long)}};
template<>
basicAC_F0::name_and_type PopenMeditMeshT_Op< MeshL, v_fesL >::name_param[] = {
  {"order", &typeid(long)},
  {"meditff", &typeid(string *)},
  {"save", &typeid(string *)},
  {"wait", &typeid(bool)},
  {"bin", &typeid(long)}};

template< class MMesh, class v_fes >
AnyType PopenMeditMeshT_Op< MMesh, v_fes >::operator( )(Stack stack) const {
  if (NoGraphicWindow) {
    return Nothing;
  }

  typedef typename MMesh::Element T;
  typedef typename MMesh::BorderElement B;
  typedef typename MMesh::Vertex V;
  typedef typename MMesh::Element::RdHat TRdHat;
  typedef typename MMesh::BorderElement::RdHat BRdHat;

  MeshPoint *mp(MeshPointStack(stack));
  long order(arg(0, stack, 1L));
  //
  int ver = GmfFloat;
  int dimp = 3;
  float fx, fy, fz;
  //
  long valsortie = 0;
  int typsol, nbsol;
  nbsol = offset - 1;

  int TypTab[l.size( ) - 1];

  for (size_t i = 0; i < l.size( ) - 1; i++) {
    TypTab[i] = l[i + 1].what;
  }
  // string stringemptymedit= string("medit");
  string *ffname = GetAny< string * >((*filename)(stack));
  string *meditff(arg(1, stack, &stringffmedit));
  long filebin(arg(4, stack, 1L));
  int smedit = max(1, nbsol);
  char *commandline = meditcmd(filebin, nbsol, smedit, *meditff, *ffname);
  printf("version de medit %s\n", commandline);
  if (verbosity > 2) cout << "number of solution = " << offset - 1 << endl;

  if (verbosity > 2) cout << "number of mesh     = " << nbTh << endl;

  // lecture des differents maillages
  int nv = 0, nt = 0, nbe = 0;    // sommet, triangles, arretes du maillage unifies

  for (size_t i = 0; i < l.size( ); i = i + offset) {
    if (l[i].what != 0) cerr << "this element is not a mesh" << i << endl;

    const MMesh &Thtmp = l[i].evalm(0, stack);
    nv += Thtmp.nv;
    nt += Thtmp.nt;
    nbe += Thtmp.nbe;
  }
  V *v = new V[nv];
  T *t;
  if (nt != 0) t = new T[nt];

  B *b = new B[nbe];
  T *tt = t;
  B *bb = b;
  int iv = 0, it = 0, ibe = 0;
  int *numTht = new int[nt];    // numero of Th assoctiated with a surface triangle
  int jt = 0;

  for (size_t i = 0; i < l.size( ); i = i + offset) {
    int nvtmp = iv;
    const MMesh &Thtmp = l[i].evalm(0, stack);

    for (int ii = 0; ii < Thtmp.nv; ii++) {
      const V &vi(Thtmp.vertices[ii]);
      v[iv].x = vi.x;
      v[iv].y = vi.y;
      v[iv].z = vi.z;
      v[iv++].lab = vi.lab;
    }

    for (int ii = 0; ii < Thtmp.nt; ++ii) {
      const T &vi(Thtmp.elements[ii]);
      int iv[T::nv];
      for (int itt = 0; itt < T::nv; itt++) iv[itt] = nvtmp + Thtmp.operator( )(vi[itt]);
      (*tt++).set(v, iv, vi.lab);
      numTht[it] = jt;
      it++;
    }
    for (int ii = 0; ii < Thtmp.nbe; ++ii) {
      const B &vi(Thtmp.be(ii));
      int iv[B::nv];
      for (int ib = 0; ib < B::nv; ++ib) iv[ib] = nvtmp + Thtmp.operator( )(vi[ib]);
      (*bb++).set(v, iv, vi.lab);
      ibe++;
    }
    jt++;
  }

  assert(it == nt);
  assert(iv == nv);
  assert(ibe == nbe);
  if (verbosity > 2) {
    cout << "meditff :: Value of elements: vertex " << nv << " surface Triangle " << nt
         << " boundary edge " << nbe << endl;
  }

  MMesh *pTh = new MMesh(nv, nt, nbe, v, t, b);
  MMesh &Th = *pTh;

  // determination of the number of elements to represent the solution
  int datasize;
  if (order == 0) datasize = nt;

  if (order == 1) datasize = nv;

  // cas de sauvegarde
  bool boolsave = false;
  int solnbfloat = 0;
  KNM< double > solsave(1, 1);
  string *saveff;
  KN< double > vxx = 0., vyx = 0., vyy = 0., vzx = 0., vzy = 0., vzz = 0.;

  if (nbsol > 0) {
    vxx.init(datasize);
    vyx.init(datasize);
    vyy.init(datasize);
    vzx.init(datasize);
    vzy.init(datasize);
    vzz.init(datasize);

    if (nargs[2]) {
      boolsave = true;
      saveff = GetAny< string * >((*nargs[2])(stack));

      for (size_t i = 0; i < offset; i++) solnbfloat = solnbfloat + l[i].nbfloat;

      solsave.init(solnbfloat, datasize);
      solsave = 0.;
    }
  }
  int nboftmp = 0;
  FILE *popenstream = popen(commandline, MODE_WRITE_BINARY);
  if (!popenstream) {
    cerr << " Error popen : " << commandline << endl;
    exit(1);
  }

  // mesh
  int jojo1;
  for (int jojo = 0; jojo < smedit; jojo++) {
    if (filebin) {
      int cod = 1;
      int KwdCod;
      int NulPos = 0;

      // determination of number solutions associated with a mesh
      fwrite((unsigned char *)&cod, WrdSiz, 1, popenstream);
      fwrite((unsigned char *)&ver, WrdSiz, 1, popenstream);
      KwdCod = GmfDimension;
      fwrite((unsigned char *)&KwdCod, WrdSiz, 1, popenstream);
      fwrite((unsigned char *)&NulPos, WrdSiz, 1, popenstream);
      fwrite((unsigned char *)&dimp, WrdSiz, 1, popenstream);

      // vertex
      KwdCod = GmfVertices;
      fwrite((unsigned char *)&KwdCod, WrdSiz, 1, popenstream);
      fwrite((unsigned char *)&NulPos, WrdSiz, 1, popenstream);
      fwrite((unsigned char *)&nv, WrdSiz, 1, popenstream);

      for (int k = 0; k < nv; k++) {
        const V &P = Th.vertices[k];
        fx = P.x;
        fy = P.y;
        fz = P.z;
        fwrite((unsigned char *)&fx, WrdSiz, 1, popenstream);
        fwrite((unsigned char *)&fy, WrdSiz, 1, popenstream);
        fwrite((unsigned char *)&fz, WrdSiz, 1, popenstream);
        fwrite((unsigned char *)&(P.lab), WrdSiz, 1, popenstream);
      }

      // element
      if (std::is_same< MMesh, MeshS >::value)
        KwdCod = GmfTriangles;
      else if (std::is_same< MMesh, MeshL >::value)
        KwdCod = GmfEdges;
      fwrite((unsigned char *)&KwdCod, WrdSiz, 1, popenstream);
      fwrite((unsigned char *)&NulPos, WrdSiz, 1, popenstream);
      fwrite((unsigned char *)&nt, WrdSiz, 1, popenstream);

      for (int k = 0; k < nt; k++) {
        const T &K(Th.elements[k]);
        int iv[T::nv];
        for (int it = 0; it < T::nv; it++) {
          iv[it] = Th.operator( )(K[it]) + 1;
          fwrite((unsigned char *)&iv[it], WrdSiz, 1, popenstream);
        }
        fwrite((unsigned char *)&(K.lab), WrdSiz, 1, popenstream);
      }

      if (std::is_same< MMesh, MeshS >::value) {
        KwdCod = GmfEdges;
        fwrite((unsigned char *)&KwdCod, WrdSiz, 1, popenstream);
        fwrite((unsigned char *)&NulPos, WrdSiz, 1, popenstream);
        fwrite((unsigned char *)&nbe, WrdSiz, 1, popenstream);

        for (int k = 0; k < nbe; k++) {
          const B &K(Th.be(k));
          int iv[B::nv];
          for (int ib = 0; ib < B::nv; ib++) {
            iv[ib] = Th.operator( )(K[ib]) + 1;
            fwrite((unsigned char *)&iv[ib], WrdSiz, 1, popenstream);
          }
          fwrite((unsigned char *)&(K.lab), WrdSiz, 1, popenstream);
        }
      }
      // End
      KwdCod = GmfEnd;
      fwrite((unsigned char *)&KwdCod, WrdSiz, 1, popenstream);
      fwrite((unsigned char *)&NulPos, WrdSiz, 1, popenstream);
    } else {
      fprintf(popenstream, "MeshVersionFormatted\n");
      fprintf(popenstream, "%i\n", ver);
      fprintf(popenstream, "Dimension\n");
      fprintf(popenstream, "%i\n", dimp);
      fprintf(popenstream, "Vertices\n");
      fprintf(popenstream, "%i\n", nv);

      for (int k = 0; k < nv; k++) {
        const V &P = Th.vertices[k];
        fx = P.x;
        fy = P.y;
        fz = P.z;
        fprintf(popenstream, "%f %f %f %i\n", fx, fy, fz, P.lab);
      }
      if (std::is_same< MMesh, MeshS >::value)
        fprintf(popenstream, "Triangles\n");
      else if (std::is_same< MMesh, MeshL >::value)
        fprintf(popenstream, "Edges\n");
      fprintf(popenstream, "%i\n", nt);

      for (int k = 0; k < nt; k++) {
        const T &K(Th.elements[k]);
        int iv[T::nv];
        for (int it = 0; it < T::nv; it++) iv[it] = Th.operator( )(K[it]) + 1;
        if (std::is_same< T, TriangleS >::value)
          fprintf(popenstream, "%i %i %i %i\n", iv[0], iv[1], iv[2], K.lab);
        else if (std::is_same< T, EdgeL >::value)
          fprintf(popenstream, "%i %i %i\n", iv[0], iv[1], K.lab);
      }

      if (std::is_same< MMesh, MeshS >::value) {
        fprintf(popenstream, "Edges\n");

        fprintf(popenstream, "%i\n", nbe);

        for (int k = 0; k < nbe; k++) {
          const B &K(Th.be(k));
          int iv[B::nv];
          for (int ib = 0; ib < B::nv; ib++) iv[ib] = Th.operator( )(K[ib]) + 1;
          fprintf(popenstream, "%i %i %i\n", iv[0], iv[1], K.lab);
        }
      }

      fprintf(popenstream, "End");
    }

    if (nbsol > 0) {
      if (filebin) {
        int cod = 1;
        int NulPos = 0;
        int KwdCod;
        fwrite((unsigned char *)&cod, WrdSiz, 1, popenstream);
        fwrite((unsigned char *)&ver, WrdSiz, 1, popenstream);
        KwdCod = GmfDimension;
        fwrite((unsigned char *)&KwdCod, WrdSiz, 1, popenstream);
        fwrite((unsigned char *)&NulPos, WrdSiz, 1, popenstream);
        fwrite((unsigned char *)&dimp, WrdSiz, 1, popenstream);
      } else {
        fprintf(popenstream, "MeshVersionFormatted %i\n", ver);
        fprintf(popenstream, "Dimension %i\n", dimp);

        // detemination of solution
        // default scalaire // faire tableau pour plusieurs
      }

      typsol = TypTab[jojo];

      if (order == 0) {
        if (filebin) {
          int NulPos = 0;
          int KwdCod = GmfSolAtTetrahedra;
          int codtypjm = 1;
          fwrite((unsigned char *)&KwdCod, WrdSiz, 1, popenstream);
          fwrite((unsigned char *)&NulPos, WrdSiz, 1, popenstream);
          fwrite((unsigned char *)&datasize, WrdSiz, 1, popenstream);
          fwrite((unsigned char *)&codtypjm, WrdSiz, 1, popenstream);
          fwrite((unsigned char *)&typsol, WrdSiz, 1, popenstream);
        } else {
          if (std::is_same< MMesh, MeshS >::value)
            fprintf(popenstream, "SolAtTriangles\n");
          else if (std::is_same< MMesh, MeshL >::value)
            fprintf(popenstream, "SolAtEdges\n");

          fprintf(popenstream, "%i\n", datasize);
          fprintf(popenstream, "%i %i\n", 1, typsol);
        }

        double k = T::nv;
        TRdHat Cdg_hat = TRdHat::diag(1. / k);

        if (typsol == 1) {
          MeshPoint *mp3(MeshPointStack(stack));
          vxx = 0.;
          for (int it = 0; it < nt; it++) {
            jojo1 = jojo + 1 + offset * numTht[it];
            const T &K(Th.elements[it]);
            mp3->set(Th, K(Cdg_hat), Cdg_hat, K, K.lab);

            vxx[it] = l[jojo1].eval(0, stack);    // GetAny< double >( (*nargs[1])(stack) );
          }

          if (filebin) {
            for (int k = 0; k < nt; k++) {
              fwrite((unsigned char *)&(vxx[k]), WrdSiz, 2, popenstream);
            }
          } else {
            for (int k = 0; k < nt; k++) {
              fprintf(popenstream, "%f\n", vxx[k]);
            }
          }
        } else if (typsol == 2) {
          vxx = 0.;
          vyy = 0.;
          vzz = 0.;
          MeshPoint *mp3(MeshPointStack(stack));
          for (int it = 0; it < Th.nt; ++it) {
            jojo1 = jojo + 1 + offset * numTht[it];
            const T &K(Th.elements[it]);
            mp3->set(Th, K(Cdg_hat), Cdg_hat, K, K.lab);

            vxx[it] = l[jojo1].eval(0, stack);    // GetAny< double >( (*xx)(stack) );
            vyy[it] = l[jojo1].eval(1, stack);    // GetAny< double >( (*yy)(stack) );
            vzz[it] = l[jojo1].eval(2, stack);    // GetAny< double >( (*zz)(stack) );
          }

          if (filebin) {
            for (int k = 0; k < nt; k++) {
              fwrite((unsigned char *)&(vxx[k]), WrdSiz, 2, popenstream);
              fwrite((unsigned char *)&(vyy[k]), WrdSiz, 2, popenstream);
              fwrite((unsigned char *)&(vzz[k]), WrdSiz, 2, popenstream);
            }
          } else
            for (int k = 0; k < nt; k++) fprintf(popenstream, "%f %f %f\n", vxx[k], vyy[k], vzz[k]);
        } else if (typsol == 3) {
          vxx(nv) = 0.;
          vyx(nv) = 0.;
          vyy(nv) = 0.;
          vzx(nv) = 0.;
          vzy(nv) = 0.;
          vzz(nv) = 0.;
          MeshPoint *mp3(MeshPointStack(stack));

          for (int it = 0; it < Th.nt; ++it) {
            jojo1 = jojo + 1 + offset * numTht[it];
            const T &K(Th.elements[it]);
            mp3->set(Th, K(Cdg_hat), Cdg_hat, K, K.lab);

            vxx[it] = l[jojo1].eval(0, stack);    // GetAny< double >( (*tsxx)(stack) );
            vyx[it] = l[jojo1].eval(1, stack);    // GetAny< double >( (*tsyx)(stack) );
            vyy[it] = l[jojo1].eval(2, stack);    // GetAny< double >( (*tsyy)(stack) );
            vzx[it] = l[jojo1].eval(3, stack);    // GetAny< double >( (*tszx)(stack) );
            vzy[it] = l[jojo1].eval(4, stack);    // GetAny< double >( (*tszy)(stack) );
            vzz[it] = l[jojo1].eval(5, stack);    // GetAny< double >( (*tszz)(stack) );
          }

          if (filebin) {
            for (int k = 0; k < nt; k++) {
              fwrite((unsigned char *)&(vxx[k]), WrdSiz, 2, popenstream);
              fwrite((unsigned char *)&(vyx[k]), WrdSiz, 2, popenstream);
              fwrite((unsigned char *)&(vyy[k]), WrdSiz, 2, popenstream);
              fwrite((unsigned char *)&(vzx[k]), WrdSiz, 2, popenstream);
              fwrite((unsigned char *)&(vzy[k]), WrdSiz, 2, popenstream);
              fwrite((unsigned char *)&(vzz[k]), WrdSiz, 2, popenstream);
            }
          } else {
            for (int k = 0; k < nt; k++) {
              fprintf(popenstream, "%f %f %f %f %f %f\n", vxx[k], vyx[k], vyy[k], vzx[k], vzy[k],
                      vzz[k]);
            }
          }
        }
      } else if (order == 1) {
        if (filebin) {
          int NulPos = 0;
          int KwdCod = GmfSolAtVertices;
          int codtypjm = 1;
          fwrite((unsigned char *)&KwdCod, WrdSiz, 1, popenstream);
          fwrite((unsigned char *)&NulPos, WrdSiz, 1, popenstream);
          fwrite((unsigned char *)&datasize, WrdSiz, 1, popenstream);
          fwrite((unsigned char *)&codtypjm, WrdSiz, 1, popenstream);
          fwrite((unsigned char *)&typsol, WrdSiz, 1, popenstream);
        } else {
          fprintf(popenstream, "SolAtVertices\n");
          fprintf(popenstream, "%i\n", datasize);
          fprintf(popenstream, "%i %i\n", 1, typsol);
        }

        if (typsol == 1) {
          KN< int > takemesh(nv);
          MeshPoint *mp3(MeshPointStack(stack));

          takemesh = 0;
          vxx = 0.;

          for (int it = 0; it < Th.nt; it++) {
            jojo1 = jojo + 1 + offset * numTht[it];

            for (int iv = 0; iv < T::nv; iv++) {
              int i = Th(it, iv);

              mp3->setP(&Th, it, iv);
              vxx[i] =
                vxx[i] + l[jojo1].eval(0, stack);    // GetAny< double >( (*nargs[1])(stack) );

              takemesh[i] = takemesh[i] + 1;
            }
          }

          if (filebin) {
            for (int k = 0; k < nv; k++) {
              vxx[k] = vxx[k] / takemesh[k];
              fwrite((unsigned char *)&(vxx[k]), WrdSiz, 2, popenstream);
            }
          } else {
            for (int k = 0; k < nv; k++) {
              vxx[k] = vxx[k] / takemesh[k];
              fprintf(popenstream, "%f\n", vxx[k]);
            }
          }
        } else if (typsol == 2) {
          KN< int > takemesh(nv);
          MeshPoint *mp3(MeshPointStack(stack));

          vxx = 0.;
          vyy = 0.;
          vzz = 0.;
          takemesh = 0;

          for (int it = 0; it < Th.nt; ++it) {
            jojo1 = jojo + 1 + offset * numTht[it];

            for (int iv = 0; iv < T::nv; ++iv) {
              int i = Th(it, iv);

              mp3->setP(&Th, it, iv);
              vxx[i] = vxx[i] + l[jojo1].eval(0, stack);    // GetAny< double >( (*xx)(stack) );
              vyy[i] = vyy[i] + l[jojo1].eval(1, stack);    // GetAny< double >( (*yy)(stack) );
              vzz[i] = vzz[i] + l[jojo1].eval(2, stack);    // GetAny< double >( (*zz)(stack) );

              takemesh[i] = takemesh[i] + 1;
            }
          }

          if (filebin) {
            for (int k = 0; k < nv; k++) {
              vxx[k] = vxx[k] / takemesh[k];
              vyy[k] = vyy[k] / takemesh[k];
              vzz[k] = vzz[k] / takemesh[k];
              fwrite((unsigned char *)&(vxx[k]), WrdSiz, 2, popenstream);
              fwrite((unsigned char *)&(vyy[k]), WrdSiz, 2, popenstream);
              fwrite((unsigned char *)&(vzz[k]), WrdSiz, 2, popenstream);
            }
          } else {
            for (int k = 0; k < nv; k++) {
              vxx[k] = vxx[k] / takemesh[k];
              vyy[k] = vyy[k] / takemesh[k];
              vzz[k] = vzz[k] / takemesh[k];
              fprintf(popenstream, "%f %f %f\n", vxx[k], vyy[k], vzz[k]);
            }
          }
        } else if (typsol == 3) {
          KN< int > takemesh(nv);
          MeshPoint *mp3(MeshPointStack(stack));

          vxx = 0.;
          vyx = 0.;
          vyy = 0.;
          vzx = 0.;
          vzy = 0.;
          vzz = 0.;
          takemesh = 0;

          for (int it = 0; it < Th.nt; ++it) {
            jojo1 = jojo + 1 + offset * numTht[it];

            for (int iv = 0; iv < T::nv; ++iv) {
              int i = Th(it, iv);

              mp3->setP(&Th, it, iv);
              vxx[i] = vxx[i] + l[jojo1].eval(0, stack);    // GetAny< double >( (*tsxx)(stack) );
              vyx[i] = vyx[i] + l[jojo1].eval(1, stack);    // GetAny< double >( (*tsyx)(stack) );
              vyy[i] = vyy[i] + l[jojo1].eval(2, stack);    // GetAny< double >( (*tsyy)(stack) );
              vzx[i] = vzx[i] + l[jojo1].eval(3, stack);    // GetAny< double >( (*tszx)(stack) );
              vzy[i] = vzy[i] + l[jojo1].eval(4, stack);    // GetAny< double >( (*tszy)(stack) );
              vzz[i] = vzz[i] + l[jojo1].eval(5, stack);    // GetAny< double >( (*tszz)(stack) );

              takemesh[i] = takemesh[i] + 1;
            }
          }

          if (filebin) {
            for (int k = 0; k < nv; k++) {
              vxx[k] = vxx[k] / takemesh[k];
              vyx[k] = vyx[k] / takemesh[k];
              vyy[k] = vyy[k] / takemesh[k];
              vzx[k] = vzx[k] / takemesh[k];
              vzy[k] = vzy[k] / takemesh[k];
              vzz[k] = vzz[k] / takemesh[k];

              fwrite((unsigned char *)&(vxx[k]), WrdSiz, 2, popenstream);
              fwrite((unsigned char *)&(vyx[k]), WrdSiz, 2, popenstream);
              fwrite((unsigned char *)&(vyy[k]), WrdSiz, 2, popenstream);
              fwrite((unsigned char *)&(vzx[k]), WrdSiz, 2, popenstream);
              fwrite((unsigned char *)&(vzy[k]), WrdSiz, 2, popenstream);
              fwrite((unsigned char *)&(vzz[k]), WrdSiz, 2, popenstream);
            }
          } else {
            for (int k = 0; k < nv; k++) {
              vxx[k] = vxx[k] / takemesh[k];
              vyx[k] = vyx[k] / takemesh[k];
              vyy[k] = vyy[k] / takemesh[k];
              vzx[k] = vzx[k] / takemesh[k];
              vzy[k] = vzy[k] / takemesh[k];
              vzz[k] = vzz[k] / takemesh[k];

              fprintf(popenstream, "%f %f %f %f %f %f\n", vxx[k], vyx[k], vyy[k], vzx[k], vzy[k],
                      vzz[k]);
            }
          }
        }
      }

      if (filebin) {
        int NulPos = 0;
        int KwdCod = GmfEnd;
        fwrite((unsigned char *)&KwdCod, WrdSiz, 1, popenstream);
        fwrite((unsigned char *)&NulPos, WrdSiz, 1, popenstream);
      } else {
        fprintf(popenstream, "End");
      }

      if (boolsave) {
        if (verbosity > 2) {
          cout << "writing solution in file" << endl;
        }

        if (typsol == 1) {
          writetabsol(datasize, nboftmp, vxx, solsave);
          nboftmp = nboftmp + 1;
        } else if (typsol == 2) {
          writetabsol(datasize, nboftmp, vxx, vyy, vzz, solsave);
          nboftmp = nboftmp + 3;
        } else if (typsol == 3) {
          writetabsol(datasize, nboftmp, vxx, vyx, vyy, vzx, vzy, vzz, solsave);
          nboftmp = nboftmp + 6;
        }
      }
    }
  }

  delete[] numTht;
  // fermeture du stream pour popen

  // drawing part  ------------------------------

  { pclose(popenstream); }

  if (boolsave) {
      int64_t outm;
    int nbtype = nbsol;
    float *OutSolTab = new float[solnbfloat];

    if (!(outm = GmfOpenMesh(saveff->c_str( ), GmfWrite, ver, 3))) {
      cerr << "  -- MeshS::Save  UNABLE TO OPEN  :" << saveff << endl;
      exit(1);
    }

    if (order == 0) {
      if (std::is_same< MMesh, MeshS >::value)
        GmfSetKwd(outm, GmfSolAtTriangles, datasize, nbtype, TypTab);
      else if (std::is_same< MMesh, MeshL >::value)
        GmfSetKwd(outm, GmfSolAtEdges, datasize, nbtype, TypTab);
      for (int k = 0; k < datasize; k++) {
        for (int i = 0; i < solnbfloat; i++) {
          OutSolTab[i] = solsave(i, k);
        }
        if (std::is_same< MMesh, MeshS >::value)
          GmfSetLin(outm, GmfSolAtTriangles, OutSolTab);
        else if (std::is_same< MMesh, MeshL >::value)
          GmfSetLin(outm, GmfSolAtEdges, OutSolTab);
      }
    } else if (order == 1) {
      GmfSetKwd(outm, GmfSolAtVertices, datasize, nbtype, TypTab);

      for (int k = 0; k < datasize; k++) {
        for (int i = 0; i < solnbfloat; i++) {
          OutSolTab[i] = solsave(i, k);
        }

        GmfSetLin(outm, GmfSolAtVertices, OutSolTab);
      }
    }

    delete[] OutSolTab;
    GmfCloseMesh(outm);
  }

  delete[] commandline;
  delete pTh;

  return valsortie;
}

// truc pour que la fonction
// static void Load_Init() soit appele a moment du chargement dynamique
// du fichier
//

// <<dynamic_loading>>

static void Load_Init( ) {    // le constructeur qui ajoute la fonction "splitmesh3"  a freefem++
  typedef Mesh *pmesh;
  typedef Mesh3 *pmesh3;
  // typedef MeshS *pmeshS;
  // typedef MeshL *pmeshL;

  if (verbosity > 2) {
    cout << " load:popen.cpp  " << endl;
  }

  // 2D
  Global.Add("medit", "(", new OneOperatorCode< PopenMeditMesh_Op >);
  Global.Add("savesol", "(", new OneOperatorCode< datasolMesh2_Op >);

  // 3D curve
  Global.Add("medit", "(", new OneOperatorCode< PopenMeditMeshT_Op< MeshL, v_fesL > >);
  Global.Add("savesol", "(", new OneOperatorCode< datasolMeshT_Op< MeshL, v_fesL > >);

  // 3D surface/curve
  Global.Add("medit", "(", new OneOperatorCode< PopenMeditMeshT_Op< MeshS, v_fesS > >);
  Global.Add("savesol", "(", new OneOperatorCode< datasolMeshT_Op< MeshS, v_fesS > >);

  // 3D volume
  Global.Add("medit", "(", new OneOperatorCode< PopenMeditMesh3_Op< v_fes3 > >);
  Global.Add("savesol", "(", new OneOperatorCode< datasolMesh3_Op< v_fes3 > >);
  Global.Add("readsol", "(", new OneOperatorCode< readsol_Op >);
}

// <<medit_Load_Init>> static loading: calling Load_Init() from a function which is accessible from
// [[file:~/ff/src/fflib/load.cpp::static_load_medit]]

void medit_Load_Init( ) { Load_Init( ); }

// [[file:include/InitFunct.hpp::LOADFUNC]]
LOADFUNC(Load_Init)
