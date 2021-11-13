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
// SUMMARY : Add interface with partionning library scotch
// LICENSE : LGPLv3
// ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
// AUTHORS : Mathieu Cloirec
// E-MAIL  : cloirec@cines.fr

/* clang-format off */
//ff-c++-LIBRARY-dep: hdf5
//ff-c++-cpp-dep:
/* clang-format on */

#include "ff++.hpp"
#include "write_xdmf.hpp"
#include "write_hdf5.hpp"

#ifndef H5_NO_NAMESPACE
#ifndef H5_NO_STD
using std::cout;
using std::endl;
#endif    // H5_NO_STD
#endif

#include "H5Cpp.h"

#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif
#ifdef DEBUG
int debug = 1;
#else
int debug = 0;
#endif

using namespace std;

class datasolHDF5Mesh2_Op : public E_F0mps {
 public:
  typedef long Result;
  Expression eTh;
  Expression filename;
  struct Expression2 {
    long what;       // 1 scalar, 2 vector, 3 symtensor
    long nbfloat;    // 1 scalar, 2 vector (3D), 3 symtensor(3D)
    Expression e[3];
    Expression lename;
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
  datasolHDF5Mesh2_Op(const basicAC_F0 &args) : l((args.size( ) - 2) / 2) {
    int nbofsol;
    int ddim = 2;
    int stsize = 3;

    cout << " " << endl;

    if (debug == 1) {
      cout << "construction data hdf5 solution avec datasolHDF5Mesh2_Op" << endl;
      cout << "taille de args " << args.size( ) << endl;
    }

    args.SetNameParam(n_name_param, name_param, nargs);

    if (BCastTo< string * >(args[0])) {
      filename = CastTo< string * >(args[0]);
    }

    if (BCastTo< pmesh >(args[1])) {
      eTh = CastTo< pmesh >(args[1]);
    }

    nbofsol = l.size( );
    if (debug == 1) {
      cout << "hdf5 solution 2d nb sol: " << nbofsol << endl;
    }

    size_t kk = 0;

    for (size_t i = 2; i < (unsigned int)args.size( ); i = i + 2) {
      size_t jj = i - 2 - kk;
      if (BCastTo< double >(args[i])) {
        l[jj].what = 1;
        l[jj].nbfloat = 1;
        l[jj][0] = to< double >(args[i]);
        if (debug == 1) {
          cout << "hdf5 solution 2d N° " << jj << " is scalar type " << endl;
        }
      } else if (args[i].left( ) == atype< E_Array >( )) {
        const E_Array *a0 = dynamic_cast< const E_Array * >(args[i].LeftValue( ));
        if (a0->size( ) != ddim && a0->size( ) != stsize) {
          CompileError(
            "savesol in 2D: vector solution is 2 composant, tensor solution is 3 composant");
        }

        if (a0->size( ) == ddim) {
          // vector solution
          if (debug == 1) {
            cout << "hdf5 solution 2d N° " << jj << " is vector type" << endl;
          }

          l[jj].what = 2;
          l[jj].nbfloat = ddim;

          for (int j = 0; j < ddim; j++) {
            l[jj][j] = to< double >((*a0)[j]);
          }
        } else if (a0->size( ) == stsize) {
          // symmetric tensor solution
          if (debug == 1) {
            cout << "hdf5 solution 2d N° " << jj << " is tensor type" << endl;
          }

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

      if (BCastTo< string * >(args[i + 1])) {
        l[jj].lename = CastTo< string * >(args[i + 1]);
      }

      kk++;
    }
  }

  static ArrayOfaType typeargs( ) {
    return ArrayOfaType(atype< string * >( ), atype< pmesh >( ), true);
  }

  static E_F0 *f(const basicAC_F0 &args) { return new datasolHDF5Mesh2_Op(args); }

  AnyType operator( )(Stack stack) const;
};

basicAC_F0::name_and_type datasolHDF5Mesh2_Op::name_param[] = {{"order", &typeid(long)}};
AnyType datasolHDF5Mesh2_Op::operator( )(Stack stack) const {
  // Hyp A
  // -----
  // A priori, paraview - hdf5 -xdmf impose qu un vecteur possede 3 composantes.
  // Donc, pour stocker la solution, on transforme le resultat produit par le
  // code, i.e. vecteur 2D, en vecteur 3D en initialisant le vecteur 3D a zero.
  //

  // Hyp B
  // -----
  // Etant donne que le tenseur est un tenseur 2d symetrique,
  // on fait le choix, ici, pour paraview, de le representer sous forme
  // d'un vecteur a 3 composantes Exx,Eyy,Exy=Eyx
  // un autre choix possible serait de stocker les valeurs sous
  // un tenseur 3x3 avec Ezz=Exz=Ezx=Eyz=Ezy=0

  // mp & mps set but not used
  // MeshPoint *mp(MeshPointStack(stack)), mps=*mp;

  const Mesh *pTh = GetAny< const Mesh * >((*eTh)(stack));
  string *ffname = GetAny< string * >((*filename)(stack));

  ffassert(pTh);
  const Mesh &Th = *pTh;
  int nt = Th.nt;
  int nv = Th.nv;
  int nbsol;
  int solnbfloat;
  int resultorder = arg(0, stack, 1L);
  string *datafieldname;
  long longdefault = 0;

  if (verbosity > 2) {
    cout << "filename data hdf5 solution () : " << ffname << endl;
    cout << "hdf5 solution () nb vertices: " << nv << endl;
    cout << "hdf5 solution () nb triangles: " << nt << endl;
    cout << "hdf5 solution () nb of fields: " << l.size( ) << endl;
  }

  // write xdmf sol file
  WriteXdmf *XdmfSolFile2D = new WriteXdmf(ffname->c_str( ), nt, nv);
  XdmfSolFile2D->WriteXdmfSolFile2DInit( );

  for (size_t i = 0; i < l.size( ); i++) {
    // Hyp A
    int trans = -1;
    if (l[i].nbfloat == 2) {
      trans = 3;
    } else {
      trans = l[i].nbfloat;
    }

    datafieldname = GetAny< string * >((*(l[i].lename))(stack));
    XdmfSolFile2D->WriteXdmfSolFile2DAddField(datafieldname, (l[i].what - 1), resultorder, trans);
  }

  XdmfSolFile2D->WriteXdmfSolFile2DFinalize( );
  delete XdmfSolFile2D;

  // write hdf5 sol file
  WriteHdf5 *Hdf5SolFile2D = new WriteHdf5(ffname->c_str( ), nt, nv);
  Hdf5SolFile2D->WriteHdf5SolFile2DInit( );

  solnbfloat = 0;

  for (size_t ii = 0; ii < l.size( ); ii++) {
    datafieldname = GetAny< string * >((*(l[ii].lename))(stack));

    if (resultorder == 0) {
      // ordre 0
      // a priori la solution est par triangle

      // Hyp A
      int trans = -1;
      if (l[ii].nbfloat == 2) {
        trans = 3;
      } else {
        trans = l[ii].nbfloat;
      }

      // initialisation a 0 du field tab
      float *tab_vfield;
      tab_vfield = new float[nt * trans];
      memset(tab_vfield, 0, sizeof(float) * nt * trans);

      solnbfloat = l[ii].nbfloat;
      nbsol = nt;
      KN< double > valsol(solnbfloat * nbsol);
      MeshPoint *mp3(MeshPointStack(stack));
      R2 Cdg_hat = R2(1. / 3., 1. / 3.);

      // boucle sur les triangles
      for (int it = 0; it < nt; it++) {
        int h = 0;
        const Mesh::Triangle &K(Th.t(it));
        mp3->set(Th, K(Cdg_hat), Cdg_hat, K, K.lab);

        // boucle sur chaque champ des triangles
        for (size_t i = 0; i < l.size( ); i++) {
          // boucle sur les valeurs de chaque champ des triangles
          for (size_t j = 0; j < (unsigned int)l[ii].nbfloat; j++) {
            valsol[it * solnbfloat + h] = l[ii].eval(j, stack);
            h = h + 1;
          }
        }

        assert(solnbfloat == h);
      }

      // creation du tableau de valeur des champs par element (idem valsol mais pour le cas vecteur
      // il faut 0 sur la troisieme composante : Hyp A
      for (int i = 0; i < nt; i++) {
        for (int h = 0; h < solnbfloat; h++) {
          tab_vfield[i * trans + h] = valsol[i * solnbfloat + h];
        }
      }

      Hdf5SolFile2D->WriteHdf5SolFile2DAddField(datafieldname, resultorder, trans, (l[ii].what - 1),
                                                tab_vfield);
      delete[] tab_vfield;
    } else {
      // Hyp A
      int trans = -1;
      if (l[ii].nbfloat == 2) {
        trans = 3;
      } else {
        trans = l[ii].nbfloat;
      }

      // initialisation a 0 du field tab
      float *tab_vfield;
      tab_vfield = new float[nv * trans];
      memset(tab_vfield, 0, sizeof(float) * nv * trans);

      solnbfloat = l[ii].nbfloat;
      nbsol = nv;
      KN< double > valsol(solnbfloat * nbsol);
      valsol = 0.;
      KN< int > takemesh(nbsol);
      MeshPoint *mp3(MeshPointStack(stack));
      takemesh = 0;

      // boucle sur les triangles
      for (int it = 0; it < nt; it++) {
        // boucle sur les 3 noeuds du triangle
        for (int iv = 0; iv < 3; iv++) {
          // i = numéro du noeud
          // ex : triangle 0 avec noeuds : 2 3 9
          int i = Th(it, iv);
          mp3->setP(&Th, it, iv);
          int h = 0;

          for (size_t j = 0; j < (unsigned int)l[ii].nbfloat; j++) {
            // calcul de la somme des valeurs du champ Ux (par exemple) sur un noeud
            // appartenant à plusieurs triangles
            // u_noeud_2_appartenant_au_triangle_0 + u_noeud_2_appartenant_au_triangle_23 + ...
            valsol[i * solnbfloat + h] = valsol[i * solnbfloat + h] + l[ii].eval(j, stack);
            h = h + 1;
          }

          assert(solnbfloat == h);
          takemesh[i] = takemesh[i] + 1;
        }
      }

      for (int i = 0; i < nv; i++) {
        for (int h = 0; h < solnbfloat; h++) {
          // calcul de la moyenne des valeurs du champ Ux (par exemple) sur un noeud
          // appartenant à plusieurs triangles
          valsol[i * solnbfloat + h] = valsol[i * solnbfloat + h] / takemesh[i];
        }
      }

      for (int i = 0; i < nv; i++) {
        for (int h = 0; h < solnbfloat; h++) {
          tab_vfield[i * trans + h] = valsol[i * solnbfloat + h];
        }
      }

      if (verbosity > 10) {
        for (int m = 0; m < (solnbfloat * nbsol); m++) {
          cout << "valsol[m] hdf5 solution : " << valsol[m] << " for lii : " << ii << endl;
        }
      }

      Hdf5SolFile2D->WriteHdf5SolFile2DAddField(datafieldname, resultorder, trans, (l[ii].what - 1),
                                                tab_vfield);
      delete[] tab_vfield;
    }
  }

  Hdf5SolFile2D->WriteHdf5SolFile2DFinalize( );
  delete Hdf5SolFile2D;
  return longdefault;
};

// datasolMesh3
template< class v_fes >
class datasolHDF5Mesh3_Op : public E_F0mps {
 public:
  typedef long Result;
  Expression eTh;
  Expression filename;
  struct Expression2 {
    long what;       // 1 scalar, 2 vector, 3 symtensor
    long nbfloat;    // 1 scalar, 3 vector (3D), 6 symtensor(3D)
    Expression e[6];
    Expression lename;
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
  datasolHDF5Mesh3_Op(const basicAC_F0 &args) : l((args.size( ) - 2) / 2) {
    int nbofsol;
    int ddim = 3;
    int stsize = 6;

    cout << " " << endl;
    if (debug == 1) {
      cout << "construction data hdf5 solution avec datasolHDF5Mesh3_Op" << endl;
      cout << "taille de args " << args.size( ) << endl;
    }

    args.SetNameParam(n_name_param, name_param, nargs);
    if (BCastTo< string * >(args[0])) {
      filename = CastTo< string * >(args[0]);
    }

    if (BCastTo< pmesh3 >(args[1])) {
      eTh = CastTo< pmesh3 >(args[1]);
    }

    nbofsol = l.size( );
    if (verbosity > 1) {
      cout << "hdf5 solution 3d nb sol: " << nbofsol << endl;
    }

    size_t kk = 0;

    for (size_t i = 2; i < (unsigned int)args.size( ); i = i + 2) {
      size_t jj = i - 2 - kk;
      if (BCastTo< double >(args[i])) {
        l[jj].what = 1;
        l[jj].nbfloat = 1;
        l[jj][0] = to< double >(args[i]);
        if (verbosity > 9) {
          cout << "hdf5 solution 3d N° " << jj << " is scalar type " << endl;
        }
      } else if (args[i].left( ) == atype< E_Array >( )) {
        const E_Array *a0 = dynamic_cast< const E_Array * >(args[i].LeftValue( ));
        if (a0->size( ) != ddim && a0->size( ) != stsize) {
          CompileError(
            "savesol in 3D: vector solution is 3 composant, tensor solution is 6 composant");
        }

        if (a0->size( ) == ddim) {
          // vector solution
          if (verbosity > 9) {
            cout << "hdf5 solution 3d N° " << jj << " is vector type" << endl;
          }

          l[jj].what = 2;
          l[jj].nbfloat = ddim;

          for (int j = 0; j < ddim; j++) {
            l[jj][j] = to< double >((*a0)[j]);
          }
        } else if (a0->size( ) == stsize) {
          // symmetric tensor solution
          if (verbosity > 9) {
            cout << "hdf5 solution 3d N° " << jj << " is tensor type" << endl;
          }

          l[jj].what = 3;
          l[jj].nbfloat = stsize;

          for (int j = 0; j < stsize; j++) {
            l[jj][j] = to< double >((*a0)[j]);
          }
        }
      } else {
        CompileError("savesol in 3D: Sorry no way to save this kind of data");
      }

      if (BCastTo< string * >(args[i + 1])) {
        l[jj].lename = CastTo< string * >(args[i + 1]);
      }

      kk++;
    }
  }

  static ArrayOfaType typeargs( ) {
    return ArrayOfaType(atype< string * >( ), atype< pmesh3 >( ), true);
  }

  static E_F0 *f(const basicAC_F0 &args) { return new datasolHDF5Mesh3_Op(args); }

  AnyType operator( )(Stack stack) const;
};

template< class v_fes >
basicAC_F0::name_and_type datasolHDF5Mesh3_Op< v_fes >::name_param[] = {{"order", &typeid(long)}};

template< class v_fes >
AnyType datasolHDF5Mesh3_Op< v_fes >::operator( )(Stack stack) const {
  const Mesh3 *pTh = GetAny< const Mesh3 * >((*eTh)(stack));
  string *ffname = GetAny< string * >((*filename)(stack));

  ffassert(pTh);
  const Mesh3 &Th = *pTh;
  int trans = -1;
  int nt = Th.nt;
  int nv = Th.nv;
  int nbsol;
  int solnbfloat;
  int resultorder = arg(0, stack, 1);
  long longdefault = 0;
  string *datafieldname;

  if (verbosity > 2) {
    cout << "filename data hdf5 solution () : " << ffname << endl;
    cout << "hdf5 solution () nb vertices: " << nv << endl;
    cout << "hdf5 solution () nb tetrahedrons: " << nt << endl;
    cout << "hdf5 solution () nb of fields: " << l.size( ) << endl;
  }

  // write xdmf sol file
  WriteXdmf *XdmfSolFile3D = new WriteXdmf(ffname->c_str( ), nt, nv);
  XdmfSolFile3D->WriteXdmfSolFile3DInit( );

  for (size_t i = 0; i < l.size( ); i++) {
    // Hyp A
    trans = l[i].nbfloat;
    datafieldname = GetAny< string * >((*(l[i].lename))(stack));
    XdmfSolFile3D->WriteXdmfSolFile3DAddField(datafieldname, (l[i].what - 1), resultorder, trans);
  }

  if (verbosity > 2) {
    cout << "save xdmf file solution : " << ffname << endl;
  }
  XdmfSolFile3D->WriteXdmfSolFile3DFinalize( );
  delete XdmfSolFile3D;

  // write hdf5 sol file
  WriteHdf5 *Hdf5SolFile3D = new WriteHdf5(ffname->c_str( ), nt, nv);
  Hdf5SolFile3D->WriteHdf5SolFile3DInit( );

  solnbfloat = 0;

  for (size_t ii = 0; ii < l.size( ); ii++) {
    datafieldname = GetAny< string * >((*(l[ii].lename))(stack));
    trans = l[ii].nbfloat;

    if (resultorder == 0) {
      // Tetrahedra
      // ordre 0
      float *tab_vfield;
      tab_vfield = new float[nt * trans];
      memset(tab_vfield, 0, sizeof(float) * nt * trans);
      solnbfloat = l[ii].nbfloat;

      nbsol = nt;
      KN< double > valsol(solnbfloat * nbsol);
      MeshPoint *mp3(MeshPointStack(stack));
      R3 Cdg_hat = R3(1. / 4., 1. / 4., 1. / 4.);

      // boucle sur les tetrahedres
      for (int it = 0; it < nt; it++) {
        int h = 0;
        const Tet &K(Th.elements[it]);
        mp3->set(Th, K(Cdg_hat), Cdg_hat, K, K.lab);

        // boucle sur chaque champ des tetrahedres
        for (size_t i = 0; i < l.size( ); i++) {
          // boucle sur les valeurs de chaque champ des tetrahedres
          for (size_t j = 0; j < (unsigned int)l[ii].nbfloat; j++) {
            valsol[it * solnbfloat + h] = l[ii].eval(j, stack);
            h = h + 1;
          }
        }

        assert(solnbfloat == h);
      }

      // creation du tableau de valeur des champs par element
      for (int i = 0; i < nt; i++) {
        for (int h = 0; h < solnbfloat; h++) {
          tab_vfield[i * trans + h] = valsol[i * solnbfloat + h];
        }
      }

      Hdf5SolFile3D->WriteHdf5SolFile3DAddField(datafieldname, resultorder, trans, (l[ii].what - 1),
                                                tab_vfield);
      delete[] tab_vfield;
    }

    if (resultorder == 1) {
      // ordre 1
      float *tab_vfield;
      tab_vfield = new float[nv * trans];
      memset(tab_vfield, 0, sizeof(float) * nv * trans);
      solnbfloat = l[ii].nbfloat;
      nbsol = nv;
      KN< double > valsol(solnbfloat * nbsol);
      KN< int > takemesh(nbsol);
      MeshPoint *mp3(MeshPointStack(stack));
      // R3 Cdg_hat = R3(1./4.,1./4.,1./4.);
      takemesh = 0;

      // boucle sur les tetrahedres
      for (int it = 0; it < nt; it++) {
        // boucle sur les 4 noeuds du tetrahedre
        for (int iv = 0; iv < 4; iv++) {
          int i = Th(it, iv);
          if (takemesh[i] == 0) {
            mp3->setP(&Th, it, iv);
            int h = 0;

            for (size_t j = 0; j < (unsigned int)l[ii].nbfloat; j++) {
              valsol[i * solnbfloat + h] = l[ii].eval(j, stack);
              h = h + 1;
            }

            assert(solnbfloat == h);
            takemesh[i] = takemesh[i] + 1;
          }
        }
      }

      for (int i = 0; i < nv; i++) {
        for (int h = 0; h < solnbfloat; h++) {
          tab_vfield[i * trans + h] = valsol[i * solnbfloat + h];
        }
      }

      if (verbosity > 9) {
        for (int m = 0; m < (solnbfloat * nbsol); m++) {
          cout << "valsol[m] hdf5 solution : " << valsol[m] << " for lii : " << ii << endl;
        }
      }

      Hdf5SolFile3D->WriteHdf5SolFile3DAddField(datafieldname, resultorder, trans, (l[ii].what - 1),
                                                tab_vfield);
      delete[] tab_vfield;
    }
  }

  Hdf5SolFile3D->WriteHdf5SolFile3DFinalize( );
  delete Hdf5SolFile3D;
  return longdefault;
}

static void Load_Init( ) {
  cout << " " << endl;
  cout << " ---------------------- " << endl;

  typedef const Mesh *pmesh;
  typedef const Mesh3 *pmesh3;

  if (verbosity > 2) {
    cout << " load:popen.cpp  " << endl;
  }

  // 2D
  Global.Add("savehdf5sol", "(", new OneOperatorCode< datasolHDF5Mesh2_Op >);

  // 3D
  Global.Add("savehdf5sol", "(", new OneOperatorCode< datasolHDF5Mesh3_Op< v_fes3 > >);
}

LOADFUNC(Load_Init)
