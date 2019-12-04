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
// AUTHORS : ...
// E-MAIL  : ...

/* clang-format off */
//ff-c++-LIBRARY-dep:
//ff-c++-cpp-dep:
/* clang-format on */

using namespace std;
#include <fstream>
#include <iostream>
#include <cstring>
#include <map>
#include "libmesh5.h"
#include "ufunction.hpp"
#include "error.hpp"
#include "RNM.hpp"
#include <stdlib.h>
#include "ufunction.hpp"
#include "Mesh2dn.hpp"
#include "Mesh2dn.hpp"
#include "Mesh3dn.hpp"
#include "rgraph.hpp"
#include "RNM.hpp"
#include "fem.hpp"

using namespace std;
using namespace Fem2D;

#include "TransfoMesh_v2.hpp"

Mesh3 *Transfo_Mesh3(const double &precis_mesh, const Mesh3 &Th3, const double *tab_XX,
                     const double *tab_YY, const double *tab_ZZ, int &border_only,
                     int &recollement_element, int &recollement_border, int &point_confondus_ok) {
  // cas besoin memoire important

  Mesh3 *T_Th3 = new Mesh3;
  int nv_t, nt_t, nbe_t;
  int *Numero_Som;
  int *ind_nv_t;
  int *ind_nt_t;
  int *ind_nbe_t;
  int *label_nt_t;
  int *label_nbe_t;
  int i_som, i_elem, i_border;

  Numero_Som = new int[Th3.nv];

  ind_nv_t = new int[Th3.nv];
  ind_nt_t = new int[Th3.nt];
  ind_nbe_t = new int[Th3.nbe];

  label_nt_t = new int[Th3.nt];
  label_nbe_t = new int[Th3.nbe];

  cout << "Vertex, Tetrahedra, Border : " << Th3.nv << ", " << Th3.nt << ", " << Th3.nbe << endl;

  for (int ii = 0; ii < Th3.nv; ii++) {
    Numero_Som[ii] = ii;
  }

  cout << " debut: SamePointElement " << endl;

  SamePointElement(precis_mesh, tab_XX, tab_YY, tab_ZZ, Th3, recollement_element,
                   recollement_border, point_confondus_ok, Numero_Som, ind_nv_t, ind_nt_t,
                   ind_nbe_t, label_nt_t, label_nbe_t, nv_t, nt_t, nbe_t);

  cout << " fin: SamePointElement " << endl;

  // set size of Mesh T_Th3
  T_Th3->set(nv_t, nt_t, nbe_t);

  if (verbosity > 1) {
    cout << "Transfo TH3 : Vertex, Tetrahedra, Border : "
         << "nv_t=" << nv_t << " nt_t=" << nt_t << " nbe_t=" << nbe_t << endl;
  }

  // determination of vertex
  i_som = 0;

  for (int i = 0; i < nv_t; i++) {
    int &ii = ind_nv_t[i];
    assert(Numero_Som[ii] == i_som);

    const Vertex3 &K(Th3.vertices[ii]);

    T_Th3->vertices[i_som].x = tab_XX[ii];
    T_Th3->vertices[i_som].y = tab_YY[ii];
    T_Th3->vertices[i_som].z = tab_ZZ[ii];
    T_Th3->vertices[i_som].lab = K.lab;

    i_som = i_som + 1;
  }

  cout << "i_som, nv_t=" << i_som << " " << nv_t << endl;
  assert(i_som == nv_t);

  cout << " Transfo volume elements " << endl;
  // determination of volume elements
  i_elem = 0;

  for (int i = 0; i < nt_t; i++) {
    int &ii = ind_nt_t[i];

    // creation of elements
    const Tet &K(Th3.elements[ii]);
    int iv[4];
    int lab;
    lab = label_nt_t[i];

    for (int jj = 0; jj < 4; jj++) {
      iv[jj] = Numero_Som[Th3.operator( )(K[jj])];
      assert(iv[jj] >= 0 && iv[jj] < nv_t);
    }

    T_Th3->elements[i_elem].set(T_Th3->vertices, iv, lab);
    i_elem = i_elem + 1;
  }

  assert(i_elem == nt_t);

  cout << " Transfo border elements " << endl;
  // determination of border elements
  i_border = 0;

  for (int i = 0; i < nbe_t; i++) {
    int &ii = ind_nbe_t[i];

    // creation of elements
    const Triangle3 &K(Th3.be(ii));
    int iv[3];
    int lab;

    lab = label_nbe_t[i];

    for (int jj = 0; jj < 3; jj++) {
      iv[jj] = Numero_Som[Th3.operator( )(K[jj])];
      assert(iv[jj] >= 0 && iv[jj] < nv_t);
    }

    T_Th3->be(i_border).set(T_Th3->vertices, iv, lab);
    i_border = i_border + 1;
  }

  assert(i_border == nbe_t);

  delete[] label_nt_t;
  delete[] label_nbe_t;

  delete[] ind_nv_t;
  delete[] ind_nt_t;
  delete[] ind_nbe_t;

  delete[] Numero_Som;

  return T_Th3;
}

void SamePointElement(const double &precis_mesh, const double *tab_XX, const double *tab_YY,
                      const double *tab_ZZ, const Mesh3 &Th3, int &recollement_element,
                      int &recollement_border, int &point_confondus_ok, int *Numero_Som,
                      int *ind_nv_t, int *ind_nt_t, int *ind_nbe_t, int *label_nt_t,
                      int *label_nbe_t, int &nv_t, int &nt_t, int &nbe_t) {
  double hmin;
  R3 bmin, bmax;

  // int recollement_element=1,recollement_border=1;

  cout << "  BuilBound " << endl;
  BuildBoundMinDist_th3(precis_mesh, tab_XX, tab_YY, tab_ZZ, Th3, bmin, bmax, hmin);
  cout << " =============================== " << endl;

  double bmin3[3], bmax3[3];
  bmin3[0] = bmin.x;
  bmin3[1] = bmin.y;
  bmin3[2] = bmin.z;

  bmax3[0] = bmax.x;
  bmax3[1] = bmax.y;
  bmax3[2] = bmax.z;
  cout << "  OrderVertexTransfo_hcode gtree " << endl;
  OrderVertexTransfo_hcode_nv_gtree(Th3.nv, bmin, bmax, hmin, tab_XX, tab_YY, tab_ZZ, Numero_Som,
                                    ind_nv_t, nv_t);
  cout << "fin order vertex gtree: nv_t=" << nv_t << endl;
  cout << " =============================== " << endl;

  /* determination de nt_t et de nbe_t*/
  int i_elem, i_border;

  i_elem = 0;

  for (int ii = 0; ii < Th3.nt; ii++) {
    const Tet &K(Th3.elements[ii]);
    int iv[4];
    int Elem_ok;

    Elem_ok = 1;

    for (int jj = 0; jj < 4; jj++) {
      iv[jj] = Numero_Som[Th3.operator( )(K[jj])];
    }

    for (int jj = 0; jj < 4; jj++) {
      for (int kk = jj + 1; kk < 4; kk++) {
        if (iv[jj] == iv[kk]) {
          Elem_ok = 0;
        }
      }
    }

    if (Elem_ok == 1) {
      ind_nt_t[i_elem] = ii;
      label_nt_t[i_elem] = K.lab;
      i_elem = i_elem + 1;
    }
  }

  nt_t = i_elem;

  if (recollement_element == 1) {
    cout << "debut recollement : nt_t= " << nt_t << endl;

    int np, dim = 3;
    int *ind_np = new int[nt_t];
    int *label_t = new int[nt_t];
    double hmin_elem;
    double **Cdg_t = new double *[nt_t];

    for (int i = 0; i < nt_t; i++) {
      Cdg_t[i] = new double[dim];
    }

    for (int i_elem = 0; i_elem < nt_t; i_elem++) {
      int &ii = ind_nt_t[i_elem];
      const Tet &K(Th3.elements[ii]);
      int iv[4];

      for (int jj = 0; jj < 4; jj++) {
        iv[jj] = Th3.operator( )(K[jj]);
      }

      Cdg_t[i_elem][0] = (tab_XX[iv[0]] + tab_XX[iv[1]] + tab_XX[iv[2]] + tab_XX[iv[3]]) / 4.;
      Cdg_t[i_elem][1] = (tab_YY[iv[0]] + tab_YY[iv[1]] + tab_YY[iv[2]] + tab_YY[iv[3]]) / 4.;
      Cdg_t[i_elem][2] = (tab_ZZ[iv[0]] + tab_ZZ[iv[1]] + tab_ZZ[iv[2]] + tab_ZZ[iv[3]]) / 4.;
      label_t[i_elem] = K.lab;
    }

    hmin_elem = hmin / 4;
    PointCommun_hcode_gtree(dim, nt_t, 0, Cdg_t, label_t, bmin, bmax, hmin_elem, ind_np, label_nt_t,
                            np);    // nv

    assert(np <= nt_t);

    int *ind_nt_t_tmp = new int[np];

    for (int i_elem = 0; i_elem < np; i_elem++) {
      assert(ind_np[i_elem] >= 0 && ind_np[i_elem] <= nt_t);
      ind_nt_t_tmp[i_elem] = ind_nt_t[ind_np[i_elem]];
    }

    for (int i_elem = 0; i_elem < np; i_elem++) {
      ind_nt_t[i_elem] = ind_nt_t_tmp[i_elem];
    }

    delete[] ind_np;
    delete[] label_t;
    for (int i = 0; i < nt_t; i++) {
      delete[] Cdg_t[i];
    }
    delete[] Cdg_t;
    delete[] ind_nt_t_tmp;

    nt_t = np;

    cout << "fin recollement : nt_t= " << nt_t << endl;
  }

  // determination of border elements
  i_border = 0;

  for (int ii = 0; ii < Th3.nbe; ii++) {
    int Border_ok = 1;

    const Triangle3 &K(Th3.be(ii));
    int iv[3];

    for (int jj = 0; jj < 3; jj++) {
      iv[jj] = Numero_Som[Th3.operator( )(K[jj])];
      assert(iv[jj] >= 0 && iv[jj] < nv_t);
    }

    for (int jj = 0; jj < 3; jj++) {
      for (int kk = jj + 1; kk < 3; kk++) {
        if (iv[jj] == iv[kk]) {
          Border_ok = 0;
        }
      }
    }

    if (Border_ok == 1) {
      ind_nbe_t[i_border] = ii;
      label_nbe_t[i_border] = K.lab;
      i_border = i_border + 1;
    }
  }

  nbe_t = i_border;

  if (recollement_border == 1) {
    cout << "debut recollement : nbe_t= " << nbe_t << endl;

    int np, dim = 3;
    int *ind_np = new int[nbe_t];
    double hmin_border;
    double **Cdg_be = new double *[nbe_t];
    int *label_be = new int[nbe_t];

    for (int i = 0; i < nbe_t; i++) {
      Cdg_be[i] = new double[dim];
    }

    for (int i_border = 0; i_border < nbe_t; i_border++) {
      int &ii = ind_nbe_t[i_border];
      const Triangle3 &K(Th3.be(ii));
      int iv[3];

      for (int jj = 0; jj < 3; jj++) {
        iv[jj] = Th3.operator( )(K[jj]);
      }

      Cdg_be[i_border][0] =
        (tab_XX[iv[0]] + tab_XX[iv[1]] + tab_XX[iv[2]]) /
        3.;    // ( Th3.vertices[iv[0]].x + Th3.vertices[iv[1]].x + Th3.vertices[iv[2]].x )/3.;
      Cdg_be[i_border][1] =
        (tab_YY[iv[0]] + tab_YY[iv[1]] + tab_YY[iv[2]]) /
        3.;    // ( Th3.vertices[iv[0]].y + Th3.vertices[iv[1]].y + Th3.vertices[iv[2]].y )/3.;
      Cdg_be[i_border][2] =
        (tab_ZZ[iv[0]] + tab_ZZ[iv[1]] + tab_ZZ[iv[2]]) /
        3.;    // ( Th3.vertices[iv[0]].z + Th3.vertices[iv[1]].z + Th3.vertices[iv[2]].z )/3.;

      label_be[i_border] = K.lab;
    }

    hmin_border = hmin / 3.;
    cout << "hmin_border=" << hmin_border << endl;

    cout << "appele de PointCommun_hcode := " << point_confondus_ok << endl;
    // PointCommun_hcode( dim, nbe_t, point_confondus_ok, Cdg_be, bmin3, bmax3, hmin_border, ind_np,
    // np);
    PointCommun_hcode_gtree(dim, nbe_t, point_confondus_ok, Cdg_be, label_be, bmin, bmax,
                            hmin_border, ind_np, label_nbe_t, np);
    cout << "fin appele de PointCommun_hcode" << endl;

    assert(np <= nbe_t);

    int *ind_nbe_t_tmp = new int[np];

    for (int i_border = 0; i_border < np; i_border++) {
      ind_nbe_t_tmp[i_border] = ind_nbe_t[ind_np[i_border]];
    }

    for (int i_border = 0; i_border < np; i_border++) {
      ind_nbe_t[i_border] = ind_nbe_t_tmp[i_border];
    }

    delete[] ind_np;
    for (int i = 0; i < nbe_t; i++) {
      delete[] Cdg_be[i];
    }
    delete[] Cdg_be;
    delete[] label_be;
    delete[] ind_nbe_t_tmp;

    nbe_t = np;

    cout << "fin recollement : nbe_t= " << nbe_t << endl;

    // Affectation de la nouvelle valeur du label
  }
}

// 3D surface

Mesh3 *Transfo_Mesh3_surf(const double &precis_mesh, const Mesh3 &Th3, const double *tab_XX,
                          const double *tab_YY, const double *tab_ZZ, int &recollement_border,
                          int &point_confondus_ok) {
  // cas besoin memoire important

  Mesh3 *T_Th3 = new Mesh3;
  int nv_t, nbe_t;
  int nt_t = 0;
  int *Numero_Som;
  int *ind_nv_t;
  int *ind_nbe_t;
  int *label_nbe_t;
  int i_som, i_elem, i_border;

  assert(Th3.nt == 0);

  Numero_Som = new int[Th3.nv];
  ind_nv_t = new int[Th3.nv];
  ind_nbe_t = new int[Th3.nbe];
  label_nbe_t = new int[Th3.nbe];

  cout << "Vertex, Tetrahedra, Border : " << Th3.nv << ", " << Th3.nt << ", " << Th3.nbe << endl;

  for (int ii = 0; ii < Th3.nv; ii++) {
    Numero_Som[ii] = ii;
  }

  cout << " debut: SamePointElement " << endl;

  SamePointElement_surf(precis_mesh, tab_XX, tab_YY, tab_ZZ, Th3, recollement_border,
                        point_confondus_ok, Numero_Som, ind_nv_t, ind_nbe_t, label_nbe_t, nv_t,
                        nbe_t);

  cout << " fin: SamePointElement " << endl;

  // set size of Mesh T_Th3
  T_Th3->set(nv_t, nt_t, nbe_t);

  if (verbosity > 1) {
    cout << "Transfo TH3 : Vertex, Tetrahedra, Border : "
         << "nv_t=" << nv_t << " nt_t=" << nt_t << " nbe_t=" << nbe_t << endl;
  }

  // determination of vertex
  i_som = 0;

  for (int i = 0; i < nv_t; i++) {
    int &ii = ind_nv_t[i];
    assert(Numero_Som[ii] == i_som);

    const Vertex3 &K(Th3.vertices[ii]);

    T_Th3->vertices[i_som].x = tab_XX[ii];
    T_Th3->vertices[i_som].y = tab_YY[ii];
    T_Th3->vertices[i_som].z = tab_ZZ[ii];
    T_Th3->vertices[i_som].lab = K.lab;

    i_som = i_som + 1;
  }

  cout << "i_som, nv_t=" << i_som << " " << nv_t << endl;
  assert(i_som == nv_t);

  cout << " Transfo border elements " << endl;
  // determination of border elements
  i_border = 0;

  for (int i = 0; i < nbe_t; i++) {
    int &ii = ind_nbe_t[i];

    // creation of elements
    const Triangle3 &K(Th3.be(ii));
    int iv[3];
    int lab;

    lab = label_nbe_t[i];

    for (int jj = 0; jj < 3; jj++) {
      iv[jj] = Numero_Som[Th3.operator( )(K[jj])];
      assert(iv[jj] >= 0 && iv[jj] <= nv_t);
    }

    T_Th3->be(i_border).set(T_Th3->vertices, iv, lab);
    i_border = i_border + 1;
  }

  assert(i_border == nbe_t);

  delete[] Numero_Som;
  delete[] ind_nv_t;
  delete[] ind_nbe_t;
  delete[] label_nbe_t;

  return T_Th3;
}

void SamePointElement_surf(const double &precis_mesh, const double *tab_XX, const double *tab_YY,
                           const double *tab_ZZ, const Mesh3 &Th3, int &recollement_border,
                           int &point_confondus_ok, int *Numero_Som, int *ind_nv_t, int *ind_nbe_t,
                           int *label_nbe_t, int &nv_t, int &nbe_t) {
  int Elem_ok;
  double hmin, hmin_elem;
  R3 bmin, bmax;

  cout << "  OrderVertexTransfo_hcode gtree " << endl;
  BuildBoundMinDist_th3(precis_mesh, tab_XX, tab_YY, tab_ZZ, Th3, bmin, bmax, hmin);
  cout << " =============================== " << endl;

  double bmin3[3], bmax3[3];
  bmin3[0] = bmin.x;
  bmin3[1] = bmin.y;
  bmin3[2] = bmin.z;

  bmax3[0] = bmax.x;
  bmax3[1] = bmax.y;
  bmax3[2] = bmax.z;

  cout << "  OrderVertexTransfo_hcode gtree " << endl;
  OrderVertexTransfo_hcode_nv_gtree(Th3.nv, bmin, bmax, hmin, tab_XX, tab_YY, tab_ZZ, Numero_Som,
                                    ind_nv_t, nv_t);
  cout << "fin order vertex gtree: nv_t=" << nv_t << endl;

  cout << " =============================== " << endl;

  /* determination de nt_t et de nbe_t*/
  int i_border;

  // determination of border elements
  i_border = 0;

  for (int ii = 0; ii < Th3.nbe; ii++) {
    int Border_ok = 1;

    const Triangle3 &K(Th3.be(ii));
    int iv[3];

    for (int jj = 0; jj < 3; jj++) {
      iv[jj] = Numero_Som[Th3.operator( )(K[jj])];
    }

    for (int jj = 0; jj < 3; jj++) {
      for (int kk = jj + 1; kk < 3; kk++) {
        if (iv[jj] == iv[kk]) {
          Border_ok = 0;
        }
      }
    }

    if (Border_ok == 1) {
      ind_nbe_t[i_border] = ii;
      label_nbe_t[i_border] = K.lab;
      i_border = i_border + 1;
    }
  }

  nbe_t = i_border;

  if (recollement_border == 1) {
    // int point_confondus_ok = 1;
    cout << "debut recollement : nbe_t= " << nbe_t << endl;

    int np, dim = 3;
    int *ind_np = new int[nbe_t];
    double hmin_border;
    double **Cdg_be = new double *[nbe_t];
    int *label_be = new int[nbe_t];

    for (int i = 0; i < nbe_t; i++) {
      Cdg_be[i] = new double[dim];
    }

    for (int i_border = 0; i_border < nbe_t; i_border++) {
      int &ii = ind_nbe_t[i_border];
      const Triangle3 &K(Th3.be(ii));
      int iv[3];

      for (int jj = 0; jj < 3; jj++) {
        iv[jj] = Th3.operator( )(K[jj]);
      }

      Cdg_be[i_border][0] =
        (tab_XX[iv[0]] + tab_XX[iv[1]] + tab_XX[iv[2]]) /
        3.;    // ( Th3.vertices[iv[0]].x + Th3.vertices[iv[1]].x + Th3.vertices[iv[2]].x )/3.;
      Cdg_be[i_border][1] =
        (tab_YY[iv[0]] + tab_YY[iv[1]] + tab_YY[iv[2]]) /
        3.;    // ( Th3.vertices[iv[0]].y + Th3.vertices[iv[1]].y + Th3.vertices[iv[2]].y )/3.;
      Cdg_be[i_border][2] =
        (tab_ZZ[iv[0]] + tab_ZZ[iv[1]] + tab_ZZ[iv[2]]) /
        3.;    // ( Th3.vertices[iv[0]].z + Th3.vertices[iv[1]].z + Th3.vertices[iv[2]].z )/3.;

      label_be[i_border] = K.lab;
    }

    hmin_border = hmin / 3.;
    cout << "hmin_border=" << hmin_border << endl;

    cout << "appele de PointCommun_hcode := " << point_confondus_ok << endl;
    // PointCommun_hcode( dim, nbe_t, point_confondus_ok, Cdg_be, bmin3, bmax3, hmin_border, ind_np,
    // np);
    PointCommun_hcode_gtree(dim, nbe_t, point_confondus_ok, Cdg_be, label_be, bmin, bmax,
                            hmin_border, ind_np, label_nbe_t, np);
    cout << "fin appele de PointCommun_hcode" << endl;

    assert(np <= nbe_t);

    int *ind_nbe_t_tmp = new int[np];

    for (int i_border = 0; i_border < np; i_border++) {
      ind_nbe_t_tmp[i_border] = ind_nbe_t[ind_np[i_border]];
    }

    for (int i_border = 0; i_border < np; i_border++) {
      ind_nbe_t[i_border] = ind_nbe_t_tmp[i_border];
    }

    delete[] ind_np;
    for (int i = 0; i < nbe_t; i++) {
      delete[] Cdg_be[i];
    }
    delete[] Cdg_be;
    delete[] label_be;

    delete[] ind_nbe_t_tmp;

    nbe_t = np;

    cout << "fin recollement : nbe_t= " << nbe_t << endl;

    // Affectation de la nouvelle valeur du label
  }
}

// CAS 2D :

/* Creation de maillage 2D */
/*
 * Mesh3 * Transfo_Mesh2(const Mesh2 & Th2, const double *tab_XX, const double *tab_YY, const double
 * *tab_ZZ, int &border_only){
 *      // cas besoin memoire important
 *      Mesh3 *T_Th3= new Mesh3;
 *      int nv_t,nt_t,nbe_t;
 *      int* Numero_Som;
 *
 *      int* ind_nv_t;
 *      int* ind_nt_t;
 *      int* ind_nbe_t;
 *
 *      int i_som, i_border;
 *
 *      int recollement_border, point_confondus_ok;
 *
 *
 *      Numero_Som = new int[Th2.nv];
 *
 *      ind_nv_t   = new int[Th2.nv];
 *  ind_nbe_t  = new int[Th2.nt];
 *
 *  cout << Th2.nv << " "<<Th2.nt<< " " << Th2.nbe<< endl;
 *
 *      for(int ii=0; ii<Th2.nv; ii++){
 *              Numero_Som[ii]=ii;
 *      }
 *
 *      cout <<" debut: SamePointElement " <<endl;
 *
 *      SamePointElement_Mesh2( tab_XX, tab_YY, tab_ZZ, Th2, recollement_border, point_confondus_ok,
 * Numero_Som, ind_nv_t, ind_nt_t, ind_nbe_t, nv_t, nt_t, nbe_t);
 *
 *      cout <<" fin: SamePointElement " <<endl;
 *
 *      // set size of Mesh T_Th3
 *      T_Th3->set(nv_t,nt_t,nbe_t);
 *
 *      cout << "nv_t="<< nv_t << " nt_t=" << nt_t << " nbe_t=" << nbe_t << endl;
 *      // determination of vertex
 *      i_som = 0;
 *      for(int i=0; i<nv_t; i++){
 *              int & ii = ind_nv_t[i];
 *              assert( Numero_Som[ii] == i_som );
 *              const Mesh::Vertex & K(Th2.vertices[ii]);
 *              T_Th3->vertices[i_som].x = tab_XX[ii];
 *              T_Th3->vertices[i_som].y = tab_YY[ii];
 *              T_Th3->vertices[i_som].z = tab_ZZ[ii];
 *              T_Th3->vertices[i_som].lab = K.lab;
 *
 *              i_som = i_som + 1;
 *      }
 *      cout << "i_som, nv_t=" <<i_som << " "<<nv_t << endl;
 *      assert( i_som == nv_t);
 *
 *      cout << " pas creation volume elements " << endl;
 *
 *      cout << " Transfo border elements " << endl;
 *      // determination of border elements
 *      i_border= 0;
 *      for( int i=0; i< nbe_t; i++){
 *              int & ii=ind_nbe_t[i];
 *
 *              // creation of elements
 *              const Triangle2 & K(Th2.elements[ii]);
 *              int iv[3];
 *              int lab=K.lab;
 *              for(int jj=0; jj <3; jj++){
 *                      iv[jj] = Numero_Som[ Th2.operator()(K[jj]) ];
 *              }
 *              T_Th3->be(i_border).set(T_Th3->vertices, iv, lab);
 *              i_border=i_border+1;
 *      }
 *      assert( i_border == nbe_t);
 *
 *      return T_Th3;
 * }
 */

void Transfo_Mesh2_map_face(const Mesh &Th2, map< int, int > &maptri) {
  int numero_label = 0;

  for (int ii = 0; ii < Th2.nt; ii++) {
    const Mesh::Triangle &K(Th2.t(ii));
    map< int, int >::const_iterator imap = maptri.find(K.lab);

    if (imap == maptri.end( )) {
      maptri[K.lab] = numero_label;
      numero_label = numero_label + 1;
    }
  }
}

Mesh3 *MoveMesh2_func(const double &precis_mesh, const Mesh &Th2, const double *tab_XX,
                      const double *tab_YY, const double *tab_ZZ, int &border_only,
                      int &recollement_border, int &point_confondus_ok) {
  Mesh3 *T_Th3;    //= new Mesh3;
  int nv_t, nt_t, nbe_t;
  int *Numero_Som;
  int *ind_nv_t;
  int *ind_nt_t;
  int *ind_nbe_t;
  int *label_nbe_t;

  Numero_Som = new int[Th2.nv];
  ind_nv_t = new int[Th2.nv];
  ind_nbe_t = new int[Th2.nt];
  label_nbe_t = new int[Th2.nt];

  cout << "2D: Mesh::Vertex  triangle2  border " << Th2.nv << " " << Th2.nt << " " << Th2.neb
       << endl;

  for (int ii = 0; ii < Th2.nv; ii++) {
    Numero_Som[ii] = ii;
  }

  cout << " debut: SamePointElement " << endl;

  SamePointElement_Mesh2(precis_mesh, tab_XX, tab_YY, tab_ZZ, Th2, recollement_border,
                         point_confondus_ok, Numero_Som, ind_nv_t, ind_nt_t, ind_nbe_t, label_nbe_t,
                         nv_t, nt_t, nbe_t);

  cout << " fin: SamePointElement " << endl;

  cout << "2D transfo: Mesh::Vertex  triangle2  border " << nv_t << " " << nt_t << " " << nbe_t
       << endl;

  T_Th3->set(nv_t, 0, nbe_t);

  for (int nnv = 0; nnv < nv_t; nnv++) {
    int ii = ind_nv_t[nnv];
    assert(Numero_Som[ii] == nnv);
    const Mesh::Vertex &K =
      Th2.vertices[ii];    // const Vertex2 & K(Th2.vertices[ii]); //Version Mesh2
    T_Th3->vertices[nnv].x = tab_XX[ii];
    T_Th3->vertices[nnv].y = tab_YY[ii];
    T_Th3->vertices[nnv].z = tab_ZZ[ii];
    T_Th3->vertices[nnv].lab = K.lab;
  }

  for (int ibe = 0; ibe < nbe_t; ibe++) {
    int lab;
    int iv[3];
    int ii = ind_nbe_t[ibe];
    // creation of elements
    const Mesh::Triangle &K(
      Th2.t(ii));    // const Triangle2 & K(Th2.elements[ii]); // Version Mesh2
    iv[0] = Numero_Som[Th2.operator( )(K[0])];
    iv[1] = Numero_Som[Th2.operator( )(K[1])];
    iv[2] = Numero_Som[Th2.operator( )(K[2])];

    T_Th3->be(ibe).set(T_Th3->vertices, iv, K.lab);
  }

  delete[] Numero_Som;
  delete[] ind_nv_t;
  delete[] ind_nbe_t;
  delete[] label_nbe_t;

  return T_Th3;
}

/*
 * Mesh3 * RemplissageSurf3D_tetgen(const Mesh3 & Th3, const int & label_tet){
 *
 *      Mesh3 *T_Th3= new Mesh3;
 *
 *      assert(Th3.nt == 0 );
 *      int nv_t = Th3.nv;
 *      int nt_t = Th3.nt;
 *      int nbe_t = Th3.nbe;
 *
 *      cout << "3D RemplissageSurf3D:: Vertex  triangle2  border "
 *      << nv_t << " "<< nt_t << " " << nbe_t<< endl;
 *      // Creation des tableau de tetgen
 *
 *      tetgenio in,out;
 *      //tetgenio::facet *f;
 *      //tetgenio::polygon *p;
 *
 *      cout << " tetgenio: vertex " << endl;
 *      int itet,jtet;
 * // All indices start from 1.
 *      in.firstnumber = 1;
 *      in.numberofpoints = nv_t;
 *      in.pointlist = new REAL[in.numberofpoints*3];
 *      in.pointmarkerlist = new int[in.numberofpoints];
 *      itet=0;
 *      jtet=0;
 *      for(int nnv=0; nnv < nv_t; nnv++)
 *      {
 *              in.pointlist[itet]   = Th3.vertices[nnv].x;
 *              in.pointlist[itet+1] = Th3.vertices[nnv].y;
 *              in.pointlist[itet+2] = Th3.vertices[nnv].z;
 *              in.pointmarkerlist[nnv] =  Th3.vertices[nnv].lab;
 *              itet=itet+3;
 *      }
 *      assert(itet==in.numberofpoints*3);
 *
 *      cout << " tetgenio: facet " << endl;
 *      // Version avec des facettes
 *      in.numberoffacets = nbe_t;
 *      in.facetlist = new tetgenio::facet[in.numberoffacets];
 *      in.facetmarkerlist = new int[in.numberoffacets];
 *
 *      for(int ibe=0; ibe < nbe_t; ibe++){
 *              tetgenio::facet *f;
 *              tetgenio::polygon *p;
 *              f = &in.facetlist[ibe];
 *              f->numberofpolygons = 1;
 *              f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
 *              f->numberofholes = 0;
 *              f->holelist = NULL;
 *
 *              p = &f->polygonlist[0];
 *              p->numberofvertices = 3;
 *              p->vertexlist = new int[3];
 *
 *              // creation of elements
 *              const Triangle3 & K(Th3.be(ibe)); // const Triangle2 & K(Th2.elements[ii]); //
 * Version Mesh2 p->vertexlist[0] = Th3.operator()(K[0])+1; p->vertexlist[1] =
 * Th3.operator()(K[1])+1; p->vertexlist[2] = Th3.operator()(K[2])+1;
 *
 *              for( int kkk=0; kkk<3; kkk++){
 *                      assert( p->vertexlist[kkk]<= in.numberofpoints && p->vertexlist[kkk]>0);
 *              }
 *
 *              in.facetmarkerlist[ibe] = K.lab;
 *
 *      }
 *      cout << "debut de tetrahedralize( , &in, &out);" << endl;
 *
 *      tetrahedralize("pqCVV", &in, &out);
 *
 *      cout << "fin de tetrahedralize( , &in, &out);" << endl;
 *      mesh3_tetgenio_out( out, label_tet, *T_Th3);
 *
 *      cout <<" Finish Mesh3 tetgen :: Vertex, Element, Border" << T_Th3->nv << " "<< T_Th3->nt <<
 * " " << T_Th3->nbe << endl;
 *
 *      return T_Th3;
 * }
 *
 *
 * Mesh3 * Transfo_Mesh2_tetgen(const Mesh & Th2, const double *tab_XX, const double *tab_YY, const
 * double *tab_ZZ, int &border_only, int &recollement_border, int &point_confondus_ok, const int
 * &label_tet, const map<int, int> &maptri ){ Mesh3 *T_Th3= new Mesh3; int nv_t,nt_t,nbe_t; int*
 * Numero_Som;
 *
 *      int* ind_nv_t;
 *      int* ind_nt_t;
 *      int* ind_nbe_t;
 *
 *      int* label_nbe_t;
 *
 *      //int i_som;
 *      Numero_Som = new int[Th2.nv];
 *      ind_nv_t   = new int[Th2.nv];
 *      ind_nbe_t  = new int[Th2.nt];
 *
 *      label_nbe_t = new int[Th2.nt];
 *
 *  cout << "2D: Mesh::Vertex  triangle2  border " << Th2.nv << " "<<Th2.nt<< " " << Th2.neb<< endl;
 *
 *      for(int ii=0; ii<Th2.nv; ii++){
 *              Numero_Som[ii]=ii;
 *      }
 *      cout <<" debut: SamePointElement " <<endl;
 *
 *      SamePointElement_Mesh2( tab_XX, tab_YY, tab_ZZ, Th2, recollement_border, point_confondus_ok,
 *              Numero_Som, ind_nv_t, ind_nt_t, ind_nbe_t, label_nbe_t, nv_t, nt_t, nbe_t);
 *
 *      cout <<" fin: SamePointElement " <<endl;
 *
 *      cout << "2D transfo: Mesh::Vertex  triangle2  border " << nv_t << " "<< nt_t << " " <<
 * nbe_t<< endl;
 *      // Creation des tableau de tetgen
 *
 *
 *      tetgenio in,out;
 *      //tetgenio::facet *f;
 *      //tetgenio::polygon *p;
 *
 *      cout << " tetgenio: vertex " << endl;
 *      int itet,jtet;
 * // All indices start from 1.
 *      in.firstnumber = 1;
 *      in.numberofpoints = nv_t;
 *      in.pointlist = new REAL[in.numberofpoints*3];
 *      in.pointmarkerlist = new int[in.numberofpoints];
 *      itet=0;
 *      jtet=0;
 *      for(int nnv=0; nnv < nv_t; nnv++)
 *      {
 *              int & ii = ind_nv_t[nnv];
 *              //cout << "nnv ,  ii  =" << nnv << "  " << ii << endl;
 *              //cout << "tab_XX[ii], tab_YY[ii], tab_ZZ[ii]=" <<  tab_XX[ii] << " "<< tab_YY[ii]
 * << " "<< tab_ZZ[ii] << endl; assert( Numero_Som[ii] == nnv ); const Mesh::Vertex & K =
 * Th2.vertices[ii];//const Mesh::Vertex & K(Th2.vertices[ii]); //Version Mesh2 in.pointlist[itet]
 * = tab_XX[ii]; in.pointlist[itet+1] = tab_YY[ii]; in.pointlist[itet+2] = tab_ZZ[ii];
 *              in.pointmarkerlist[nnv] =  K.lab;
 *              itet=itet+3;
 *      }
 *      assert(itet==in.numberofpoints*3);
 *
 *      cout << " tetgenio: facet " << endl;
 *      // Version avec des facettes
 *      in.numberoffacets = nbe_t;
 *      in.facetlist = new tetgenio::facet[in.numberoffacets];
 *      in.facetmarkerlist = new int[in.numberoffacets];
 *
 *      for(int ibe=0; ibe < nbe_t; ibe++){
 *              tetgenio::facet *f;
 *              tetgenio::polygon *p;
 *              f = &in.facetlist[ibe];
 *              f->numberofpolygons = 1;
 *              f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
 *              f->numberofholes = 0;
 *              f->holelist = NULL;
 *
 *              p = &f->polygonlist[0];
 *              p->numberofvertices = 3;
 *              p->vertexlist = new int[3];
 *
 *              int & ii=ind_nbe_t[ibe];
 *              // creation of elements
 *              const Mesh::Triangle & K(Th2.t(ii)); // const Triangle2 & K(Th2.elements[ii]); //
 * Version Mesh2 p->vertexlist[0] = Numero_Som[ Th2.operator()(K[0]) ]+1; p->vertexlist[1] =
 * Numero_Som[ Th2.operator()(K[1]) ]+1; p->vertexlist[2] = Numero_Som[ Th2.operator()(K[2]) ]+1;
 *
 *              for( int kkk=0; kkk<3; kkk++){
 *                      assert( p->vertexlist[kkk]<= in.numberofpoints && p->vertexlist[kkk]> 0);
 *              }
 *              map< int, int>:: const_iterator imap;
 *              imap = maptri.find(K.lab); // imap= maptri.find( label_nbe_t[ibe] );
 *              assert( imap != maptri.end());
 *              in.facetmarkerlist[ibe] = imap->second; // K.lab; // before
 *
 *      }
 *      cout << "debut de tetrahedralize( , &in, &out);" << endl;
 *
 *      tetrahedralize("pqCVV", &in, &out);
 *
 *      cout << "fin de tetrahedralize( , &in, &out);" << endl;
 *      mesh3_tetgenio_out( out, label_tet, *T_Th3);
 *
 *      cout <<" Finish Mesh3 :: Vertex, Element, Border" << T_Th3->nv << " "<< T_Th3->nt << " " <<
 * T_Th3->nbe << endl;
 *
 *      return T_Th3;
 * }
 */

void SamePointElement_Mesh2(const double &precis_mesh, const double *tab_XX, const double *tab_YY,
                            const double *tab_ZZ, const Mesh &Th2, int &recollement_border,
                            int &point_confondus_ok, int *Numero_Som, int *ind_nv_t, int *ind_nt_t,
                            int *ind_nbe_t, int *label_nbe_t, int &nv_t, int &nt_t, int &nbe_t) {
  R3 bmin, bmax;
  double hmin;

  cout << "calculus of bound and minimal distance" << endl;
  BuildBoundMinDist_th2(precis_mesh, tab_XX, tab_YY, tab_ZZ, Th2, bmin, bmax, hmin);
  // assertion pour la taille de l octree
  assert(hmin > Norme2(bmin - bmax) / 1e9);

  double bmin3[3], bmax3[3];
  bmin3[0] = bmin.x;
  bmin3[1] = bmin.y;
  bmin3[2] = bmin.z;

  bmax3[0] = bmax.x;
  bmax3[1] = bmax.y;
  bmax3[2] = bmax.z;
  cout << "debut: OrderVertexTransfo_hcode_gtree " << endl;
  OrderVertexTransfo_hcode_nv_gtree(Th2.nv, bmin, bmax, hmin, tab_XX, tab_YY, tab_ZZ, Numero_Som,
                                    ind_nv_t, nv_t);
  cout << "fin: OrderVertexTransfo_hcode_gtree " << endl;

  /* determination de nt_t et de nbe_t*/
  nt_t = 0;
  int i_border;

  // determination of border elements
  i_border = 0;

  for (int ii = 0; ii < Th2.nt; ii++) {
    int Border_ok = 1;
    const Mesh::Triangle &K(Th2.t(ii));    // const Triangle2 & K(Th2.elements[ii]); // avant Mesh2
    int iv[3];

    for (int jj = 0; jj < 3; jj++) {
      iv[jj] = Numero_Som[Th2.operator( )(K[jj])];
    }

    for (int jj = 0; jj < 3; jj++) {
      for (int kk = jj + 1; kk < 3; kk++) {
        if (iv[jj] == iv[kk]) {
          Border_ok = 0;
        }
      }
    }

    if (Border_ok == 1) {
      ind_nbe_t[i_border] = ii;
      label_nbe_t[i_border] = K.lab;
      i_border = i_border + 1;
    }
  }

  nbe_t = i_border;

  if (recollement_border == 1) {
    // int point_confondus_ok=1;
    cout << "debut recollement : nbe_t= " << nbe_t << endl;

    int np, dim = 3;
    int *ind_np = new int[nbe_t];
    int *label_be = new int[nbe_t];
    double hmin_border;
    double **Cdg_be = new double *[nbe_t];

    for (int i = 0; i < nbe_t; i++) {
      Cdg_be[i] = new double[dim];
    }

    for (int i_border = 0; i_border < nbe_t; i_border++) {
      int &ii = ind_nbe_t[i_border];
      const Mesh::Triangle &K(
        Th2.t(ii));    // const Triangle2 & K(Th2.elements[ii]);  // avant Mesh2
      int iv[3];

      for (int jj = 0; jj < 3; jj++) {
        iv[jj] = Th2.operator( )(K[jj]);
      }

      Cdg_be[i_border][0] = (tab_XX[iv[0]] + tab_XX[iv[1]] + tab_XX[iv[2]]) / 3.;
      Cdg_be[i_border][1] = (tab_YY[iv[0]] + tab_YY[iv[1]] + tab_YY[iv[2]]) / 3.;
      Cdg_be[i_border][2] = (tab_ZZ[iv[0]] + tab_ZZ[iv[1]] + tab_ZZ[iv[2]]) / 3.;

      label_be[i_border] = K.lab;
    }

    hmin_border = hmin / 3.;
    cout << "points commun " << endl;
    // PointCommun_hcode( dim, nbe_t, point_confondus_ok, Cdg_be, bmin3, bmax3, hmin_border, ind_np,
    // np); // ancien
    PointCommun_hcode_gtree(dim, nbe_t, point_confondus_ok, Cdg_be, label_be, bmin, bmax,
                            hmin_border, ind_np, label_nbe_t, np);    // new
    cout << "points commun finis " << endl;
    assert(np <= nbe_t);

    int *ind_nbe_t_tmp = new int[np];

    for (int i_border = 0; i_border < np; i_border++) {
      ind_nbe_t_tmp[i_border] = ind_nbe_t[ind_np[i_border]];
    }

    for (int i_border = 0; i_border < np; i_border++) {
      ind_nbe_t[i_border] = ind_nbe_t_tmp[i_border];
    }

    delete[] ind_np;
    delete[] label_be;
    for (int i = 0; i < nbe_t; i++) {
      delete[] Cdg_be[i];
    }
    delete[] Cdg_be;

    delete[] ind_nbe_t_tmp;

    nbe_t = np;
    cout << "fin recollement : nbe_t= " << nbe_t << endl;
  }
}

/*
 * // Fonction pour tetgen
 *
 * void mesh3_tetgenio_out(const tetgenio &out, const int & label_tet, Mesh3 & Th3)
 * {
 * int i;
 *
 * // All indices start from 1.
 * if(out.firstnumber != 1){
 *  cout << " probleme ???" << endl;
 *  exit(1);
 * }
 *
 * if(out.numberoffacets !=0){
 *  cout << "tetgen: faces non triangulaire" << endl;
 *  exit(1);
 * }
 *
 * if(out.numberofcorners !=4){
 *  cout << "tetgen: element subparametric of order 2" <<endl;
 *  exit(1);
 * }
 *
 * cout << "Th3 :: Vertex Element Border :: " << out.numberofpoints << " " <<out.numberoftetrahedra
 * << " " << out.numberoftrifaces << endl; Th3.set(out.numberofpoints, out.numberoftetrahedra,
 * out.numberoftrifaces);
 *
 * i=0;
 * for(int nnv=0; nnv < Th3.nv; nnv++){
 *  Th3.vertices[nnv].x=out.pointlist[i];
 *  Th3.vertices[nnv].y=out.pointlist[i+1];
 *  Th3.vertices[nnv].z=out.pointlist[i+2];
 *  Th3.vertices[nnv].lab=out.pointmarkerlist[nnv];
 *  i=i+3;
 * }
 *
 * i=0;
 * for(int nnt=0; nnt < Th3.nt; nnt++){
 *  int iv[4],lab;
 *  iv[0] = out.tetrahedronlist[i]-1;
 *  iv[1] = out.tetrahedronlist[i+1]-1;
 *  iv[2] = out.tetrahedronlist[i+2]-1;
 *  iv[3] = out.tetrahedronlist[i+3]-1;
 *  lab   = label_tet;
 *  //lab = out.tetrahedronmarkerlist[nnt];
 *  Th3.elements[nnt].set( Th3.vertices, iv, lab);
 *  i=i+4;
 * }
 *
 * for(int ibe=0; ibe < Th3.nbe; ibe++){
 *  int iv[3];
 *  iv[0] = out.trifacelist[3*ibe]-1;
 *  iv[1] = out.trifacelist[3*ibe+1]-1;
 *  iv[2] = out.trifacelist[3*ibe+2]-1;
 *  Th3.be(ibe).set( Th3.vertices, iv, out.trifacemarkerlist[ibe]);
 * }
 * }
 */

//= =====================
// Fin cas 2D
//= =====================
// version Mesh2
void BuildBoundMinDist_th2(const double &precis_mesh, const double *tab_XX, const double *tab_YY,
                           const double *tab_ZZ, const Mesh &Th2, R3 &bmin, R3 &bmax,
                           double &hmin) {
  // determination de la boite englobante
  double precispt;

  bmin.x = tab_XX[0];
  bmin.y = tab_YY[0];
  bmin.z = tab_ZZ[0];

  bmax.x = bmin.x;
  bmax.y = bmin.y;
  bmax.z = bmin.z;

  cout << " determination of bmin and bmax" << endl;

  for (int ii = 1; ii < Th2.nv; ii++) {
    bmin.x = min(bmin.x, tab_XX[ii]);
    bmin.y = min(bmin.y, tab_YY[ii]);
    bmin.z = min(bmin.z, tab_ZZ[ii]);

    bmax.x = max(bmax.x, tab_XX[ii]);
    bmax.y = max(bmax.y, tab_YY[ii]);
    bmax.z = max(bmax.z, tab_ZZ[ii]);
  }

  double longmini_box;    // = 1e10;

  longmini_box = pow(bmax.x - bmin.x, 2) + pow(bmax.y - bmin.y, 2) + pow(bmax.z - bmin.z, 2);
  longmini_box = sqrt(longmini_box);

  // determination de hmin
  if (precis_mesh < 0) {
    precispt = longmini_box * 1e-7;
  } else {
    precispt = precis_mesh;
  }

  hmin = 1e10;

  for (int ii = 0; ii < Th2.nt; ii++) {
    const Mesh::Triangle &K(Th2.t(ii));    // const Triangle2 & K(Th2.elements[ii]);
    double longedge;
    int iv[3];

    for (int jj = 0; jj < 3; jj++) {
      iv[jj] = Th2.operator( )(K[jj]);
    }

    for (int jj = 0; jj < 3; jj++) {
      for (int kk = jj + 1; kk < 3; kk++) {
        int &i1 = iv[jj];
        int &i2 = iv[kk];
        longedge = pow(tab_XX[i1] - tab_XX[i2], 2) + pow(tab_YY[i1] - tab_YY[i2], 2) +
                   pow(tab_ZZ[i1] - tab_ZZ[i2], 2);
        longedge = sqrt(longedge);
        // cout << "longedge=" << longedge << endl;
        if (longedge > precispt) {
          hmin = min(hmin, longedge);
        }
      }
    }
  }

  cout << "longmin_box=" << longmini_box << endl;
  cout << "hmin =" << hmin << endl;
  cout << "Norme2(bmin-bmax)=" << Norme2(bmin - bmax) << endl;
  assert(hmin < longmini_box);

  // assertion pour la taille de l octree
  assert(hmin > Norme2(bmin - bmax) / 1e9);
}

// version Mesh3

void BuildBoundMinDist_th3(const double &precis_mesh, const double *tab_XX, const double *tab_YY,
                           const double *tab_ZZ, const Mesh3 &Th3, R3 &bmin, R3 &bmax,
                           double &hmin) {
  // determination de la boite englobante
  double precispt;

  bmin.x = tab_XX[0];
  bmin.y = tab_YY[0];
  bmin.z = tab_ZZ[0];

  bmax.x = bmin.x;
  bmax.y = bmin.y;
  bmax.z = bmin.z;

  cout << " determination of bmin and bmax" << endl;

  for (int ii = 1; ii < Th3.nv; ii++) {
    bmin.x = min(bmin.x, tab_XX[ii]);
    bmin.y = min(bmin.y, tab_YY[ii]);
    bmin.z = min(bmin.z, tab_ZZ[ii]);

    bmax.x = max(bmax.x, tab_XX[ii]);
    bmax.y = max(bmax.y, tab_YY[ii]);
    bmax.z = max(bmax.z, tab_ZZ[ii]);
  }

  double longmini_box;

  longmini_box = pow(bmax.x - bmin.x, 2) + pow(bmax.y - bmin.y, 2) + pow(bmax.z - bmin.z, 2);
  longmini_box = sqrt(longmini_box);

  cout << " bmin := " << bmin.x << " " << bmin.y << " " << bmin.z << endl;
  cout << " bmax := " << bmax.x << " " << bmax.y << " " << bmax.z << endl;
  cout << " box volume :=" << longmini_box << endl;

  if (precis_mesh < 0) {
    precispt = longmini_box * 1e-7;
  } else {
    precispt = precis_mesh;
  }

  // determination de hmin

  hmin = 1e10;

  for (int ii = 0; ii < Th3.nt; ii++) {
    const Tet &K(Th3.elements[ii]);
    double longedge;
    int iv[4];

    for (int jj = 0; jj < 4; jj++) {
      iv[jj] = Th3.operator( )(K[jj]);
    }

    for (int jj = 0; jj < 4; jj++) {
      for (int kk = jj + 1; kk < 4; kk++) {
        int &i1 = iv[jj];
        int &i2 = iv[kk];
        longedge = pow(tab_XX[i1] - tab_XX[i2], 2) + pow(tab_YY[i1] - tab_YY[i2], 2) +
                   pow(tab_ZZ[i1] - tab_ZZ[i2], 2);
        longedge = sqrt(longedge);
        if (longedge > precispt) {
          hmin = min(hmin, longedge);
        }
      }
    }
  }

  cout << "longmini_box" << longmini_box << endl;
  cout << "hmin =" << hmin << endl;
  assert(hmin < longmini_box);
  cout << "longmini_box" << longmini_box << endl;
  cout << "hmin =" << hmin << endl;
  cout << "Norme2(bmin-bmax)=" << Norme2(bmin - bmax) << endl;
  // assertion pour la taille de l octree
  assert(hmin > Norme2(bmin - bmax) / 1e9);
}

//= =====================
//
//= =====================
void OrderVertexTransfo_hcode_nv(const int &tab_nv, const double *tab_XX, const double *tab_YY,
                                 const double *tab_ZZ, const double *bmin, const double *bmax,
                                 const double hmin, int *Numero_Som, int *ind_nv_t, int &nv_t) {
  size_t j[3];
  size_t k[3];
  size_t NbCode = 100000;
  int *tcode;    //= new int[NbCode];
  int *posv = new int[tab_nv];
  double epsilon = hmin / 10.;

  k[0] = int((bmax[0] - bmin[0]) / epsilon);
  k[1] = int((bmax[1] - bmin[1]) / epsilon);
  k[2] = int((bmax[2] - bmin[2]) / epsilon);

  int numberofpoints = 0;

  for (int ii = 0; ii < tab_nv; ii++) {
    int numberofpointsdiff;

    numberofpointsdiff = 0;

    for (int jj = ii + 1; jj < tab_nv; jj++) {
      double dist;    // = 0.;
      dist = pow(tab_XX[jj] - tab_XX[ii], 2) + pow(tab_YY[jj] - tab_YY[ii], 2) +
             pow(tab_ZZ[jj] - tab_ZZ[ii], 2);    // pow(Coord_Point[jj][kk]-Coord_Point[ii][kk],2);
      if (sqrt(dist) < epsilon) {
        numberofpointsdiff = 1;
      }
    }

    if (numberofpointsdiff == 0) {
      numberofpoints = numberofpoints + 1;
    }
  }

  cout << "numberofpoints " << numberofpoints << endl;
  cout << "taille boite englobante =" << endl;

  for (int ii = 0; ii < 3; ii++) {
    cout << "ii=" << ii << " " << bmin[ii] << " " << bmax[ii] << endl;
  }

  for (int ii = 0; ii < 3; ii++) {
    cout << "k[" << ii << "]= " << k[ii] << endl;
  }

  NbCode = min(4 * (k[0] + k[1] + k[2]), NbCode);
  tcode = new int[NbCode];

  /* initialisation des codes */
  for (int ii = 0; ii < NbCode; ii++) {
    tcode[ii] = -1;
  }

  for (int ii = 0; ii < tab_nv; ii++) {
    size_t i;
    // boucle dans l autre sens pour assurer l'ordre des elements pour la suite
    j[0] = int((tab_XX[ii] - bmin[0]) / epsilon);
    j[1] = int((tab_YY[ii] - bmin[1]) / epsilon);
    j[2] = int((tab_ZZ[ii] - bmin[2]) / epsilon);

    assert(j[0] <= k[0] && j[0] >= 0);
    assert(j[1] <= k[1] && j[1] >= 0);
    assert(j[2] <= k[2] && j[2] >= 0);
    i = (j[2] * (k[1] + 1) + j[1] * (k[0] + 1) + j[0]);
    i = i % NbCode;
    assert(i < NbCode);
    posv[ii] = tcode[i];
    tcode[i] = ii;
  }

  cout << " boucle numero de Sommet " << endl;

  for (int ii = 0; ii < tab_nv; ii++) {
    Numero_Som[ii] = -1;
  }

  cout << " determinations des points confondus et numerotation " << endl;

  nv_t = 0;

  for (int icode = 0; icode < NbCode; icode++) {
    double dist;

    for (int ii = tcode[icode]; ii != -1; ii = posv[ii]) {
      if (Numero_Som[ii] != -1) {
        continue;
      }

      Numero_Som[ii] = nv_t;

      for (int jj = posv[ii]; jj != -1; jj = posv[jj]) {
        if (Numero_Som[jj] != -1) {
          continue;
        }

        dist = pow(tab_XX[jj] - tab_XX[ii], 2) + pow(tab_YY[jj] - tab_YY[ii], 2) +
               pow(tab_ZZ[jj] - tab_ZZ[ii], 2);

        if (sqrt(dist) < epsilon) {
          // point semblable
          Numero_Som[jj] = Numero_Som[ii];
        }
      }

      ind_nv_t[nv_t] = ii;    // Remarque on donne a nv_t le plus grand
      nv_t++;                 // nv_t = nvt+1;
    }
  }

  cout << "nv_t = " << nv_t << " / "
       << "nv_t(anc)" << tab_nv << endl;
  assert(nv_t == numberofpoints);

  delete[] tcode;
  delete[] posv;
}

void PointCommun_hcode(const int &dim, const int &NbPoints, const int &point_confondus_ok,
                       double **Coord_Point, const double *bmin, const double *bmax,
                       const double hmin, int *ind_np, int &np) {
  size_t j[dim];
  size_t k[dim];
  size_t NbCode = 100000;
  int *tcode;    //= new int[NbCode];
  int *posv = new int[NbPoints];
  int *Numero_Som = new int[NbPoints];
  double epsilon = hmin / 10.;

  assert(dim > 1);

  for (int jj = 0; jj < dim; jj++) {
    k[jj] = int((bmax[jj] - bmin[jj]) / epsilon);
  }

  int numberofpoints = 0;

  for (int ii = 0; ii < NbPoints; ii++) {
    int numberofpointsdiff;
    numberofpointsdiff = 0;

    for (int jj = ii + 1; jj < NbPoints; jj++) {
      double dist = 0.;

      for (int kk = 0; kk < 3; kk++) {
        dist = dist + pow(Coord_Point[jj][kk] - Coord_Point[ii][kk], 2);
      }

      if (sqrt(dist) < 1e-10) {
        numberofpointsdiff = 1;
      }
    }

    if (numberofpointsdiff == 0) {
      numberofpoints = numberofpoints + 1;
    }
  }

  cout << "numberofpoints " << numberofpoints << endl;

  NbCode = min(4 * (k[0] + k[1] + k[2]), NbCode);
  cout << "NbCode=" << NbCode << endl;
  tcode = new int[NbCode];

  /* initialisation des codes */
  for (int ii = 0; ii < NbCode; ii++) {
    tcode[ii] = -1;
  }

  for (int ii = 0; ii < NbPoints; ii++) {
    size_t i;
    // boucle dans l autre sens pour assurer l'ordre des elements pour la suite

    for (int jj = 0; jj < dim; jj++) {
      j[jj] = int((Coord_Point[ii][jj] - bmin[jj]) / epsilon);
    }

    assert(j[0] <= k[0] && j[0] >= 0);
    assert(j[1] <= k[1] && j[1] >= 0);
    assert(j[2] <= k[2] && j[2] >= 0);

    i = j[0];

    for (int jj = 1; jj < dim; jj++) {
      i = i + j[jj] * (k[jj - 1] + 1);
    }

    i = i % NbCode;

    assert(i < NbCode);
    posv[ii] = tcode[i];
    tcode[i] = ii;
  }

  for (int ii = 0; ii < NbPoints; ii++) {
    ind_np[ii] = -1;
    Numero_Som[ii] = -1;
  }

  /* Resolution probleme dans le cas o� le maillage se colle */

  /* maintenant determinations des points confondus et numerotation*/

  switch (point_confondus_ok) {
    case 0:
      np = 0;

      for (int icode = 0; icode < NbCode; icode++) {
        double dist;

        for (int ii = tcode[icode]; ii != -1; ii = posv[ii]) {
          if (Numero_Som[ii] != -1) {
            continue;
          }

          Numero_Som[ii] = np;

          for (int jj = posv[ii]; jj != -1; jj = posv[jj]) {
            if (Numero_Som[jj] != -1) {
              continue;
            }

            dist = 0.;

            for (int kk = 0; kk < dim; kk++) {
              dist = dist + pow(Coord_Point[jj][kk] - Coord_Point[ii][kk], 2);
            }

            if (sqrt(dist) < epsilon) {
              // point semblable
              Numero_Som[jj] = Numero_Som[ii];
              // minimum_np = min( jj, minimum_np);
            }
          }

          ind_np[np] =
            ii;    // min(ii,minimum_np);	// Remarque on donne a np le plus petit element
          np++;    // nv_t = nvt+1;
        }
      }

      break;

    case 1:
      int point_multiple;
      np = 0;

      for (int icode = 0; icode < NbCode; icode++) {
        double dist;

        for (int ii = tcode[icode]; ii != -1; ii = posv[ii]) {
          if (Numero_Som[ii] != -1) {
            continue;
          }

          Numero_Som[ii] = np;
          point_multiple = 0;

          for (int jj = posv[ii]; jj != -1; jj = posv[jj]) {
            if (Numero_Som[jj] != -1) {
              continue;
            }

            dist = 0.;

            for (int kk = 0; kk < dim; kk++) {
              dist = dist + pow(Coord_Point[jj][kk] - Coord_Point[ii][kk], 2);
            }

            if (sqrt(dist) < epsilon) {
              // point semblable
              Numero_Som[jj] = Numero_Som[ii];
              point_multiple = 1;
            }
          }

          if (point_multiple == 0) {
            ind_np[np] =
              ii;    // min(ii,minimum_np);	// Remarque on donne a np le plus petit element
            np++;    // nv_t = nvt+1;
          }
        }
      }

      break;
    default:
      cout << " point_confondus_ok dans fonction PointCommun_hcode vaut 1 ou 0." << endl;
      exit(-1);
  }

  delete[] tcode;
  delete[] posv;
  delete[] Numero_Som;
}

// fonction avec Gtree pour FreeFem++

void OrderVertexTransfo_hcode_nv_gtree(const int &tab_nv, const R3 &bmin, const R3 &bmax,
                                       const double &hmin, const double *tab_XX,
                                       const double *tab_YY, const double *tab_ZZ, int *Numero_Som,
                                       int *ind_nv_t, int &nv_t) {
  size_t i;
  size_t j[3];
  size_t k[3];

  // parametre interne pour debugger le code
  int verifnumberofpoints;

  verifnumberofpoints = 0;

  // hmin a determiner plus haut
  assert(hmin > Norme2(bmin - bmax) / 1e9);
  double hseuil = hmin / 10.;

  Vertex3 *v = new Vertex3[tab_nv];
  EF23::GTree< Vertex3 > *gtree = new EF23::GTree< Vertex3 >(v, bmin, bmax, 0);

  cout << "taille de la boite " << endl;
  cout << bmin.x << " " << bmin.y << " " << bmin.z << endl;
  cout << bmax.x << " " << bmax.y << " " << bmax.z << endl;

  // creation of octree
  nv_t = 0;

  for (int ii = 0; ii < tab_nv; ii++) {
    const R3 r3vi(tab_XX[ii], tab_YY[ii], tab_ZZ[ii]);
    const Vertex3 &vi(r3vi);
    Vertex3 *pvi = gtree->ToClose(vi, hseuil);
    if (!pvi) {
      v[nv_t].x = vi.x;
      v[nv_t].y = vi.y;
      v[nv_t].z = vi.z;
      v[nv_t].lab = vi.lab;    // lab mis a zero par default
      ind_nv_t[nv_t] = ii;
      Numero_Som[ii] = nv_t;
      gtree->Add(v[nv_t]);
      nv_t = nv_t + 1;
    } else {
      Numero_Som[ii] = pvi - v;
    }
  }

  cout << "hseuil=" << hseuil << endl;
  cout << "nv_t = " << nv_t << " / "
       << "nv_t(anc)" << tab_nv << endl;

  if (verifnumberofpoints == 1) {
    int numberofpoints = 0;
    int numberofpointsdiff;

    for (int ii = 0; ii < tab_nv; ii++) {
      numberofpointsdiff = 0;

      for (int jj = ii + 1; jj < tab_nv; jj++) {
        double dist;    // = 0.;
        dist =
          pow(tab_XX[jj] - tab_XX[ii], 2) + pow(tab_YY[jj] - tab_YY[ii], 2) +
          pow(tab_ZZ[jj] - tab_ZZ[ii], 2);    // pow(Coord_Point[jj][kk]-Coord_Point[ii][kk],2);
        if (sqrt(dist) < hseuil) {
          numberofpointsdiff = 1;
        }
      }

      if (numberofpointsdiff == 0) {
        numberofpoints = numberofpoints + 1;
      }
    }

    cout << "numberofpoints " << numberofpoints << endl;
    cout << "taille boite englobante =" << endl;
    assert(nv_t == numberofpoints);
  }

  delete[] v;
  delete gtree;
}

void PointCommun_hcode_gtree(const int &dim, const int &NbPoints, const int &point_confondus_ok,
                             double **Coord_Point, const int *label_point, const R3 &bmin,
                             const R3 &bmax, const double &hmin, int *ind_np, int *ind_label,
                             int &np) {
  double hseuil = hmin / 10.;
  Vertex3 *v = new Vertex3[NbPoints];
  EF23::GTree< Vertex3 > *gtree = new EF23::GTree< Vertex3 >(v, bmin, bmax, 0);

  cout << "verif hmin vertex3 GTree switch" << point_confondus_ok << endl;

  int int_point_confondus_ok = point_confondus_ok;

  if (int_point_confondus_ok == 0) {
    // accepte les points double
    np = 0;

    for (int ii = 0; ii < NbPoints; ii++) {
      const R3 r3vi(Coord_Point[ii][0], Coord_Point[ii][1], Coord_Point[ii][2]);
      const Vertex3 &vi(r3vi);
      Vertex3 *pvi = gtree->ToClose(vi, hseuil);

      if (!pvi) {
        v[np].x = vi.x;
        v[np].y = vi.y;
        v[np].z = vi.z;
        v[np].lab = vi.lab;    // lab mis a zero par default
        ind_np[np] = ii;
        ind_label[np] = label_point[ii];
        gtree->Add(v[np++]);
      } else {
        ind_label[pvi - v] = min(ind_label[pvi - v], label_point[ii]);
      }
    }

    cout << "np=" << np << endl;
  }

  if (int_point_confondus_ok == 1) {
    // accepte les points double sont enleves
    np = 0;

    for (int ii = 0; ii < NbPoints; ii++) {
      const R3 r3vi(Coord_Point[ii][0], Coord_Point[ii][1], Coord_Point[ii][2]);
      // int label = label_point[ii];
      const Vertex3 &vi(r3vi);
      Vertex3 *pvi = gtree->ToClose(vi, hseuil);

      if (!pvi) {
        v[np].x = vi.x;
        v[np].y = vi.y;
        v[np].z = vi.z;
        v[np].lab = vi.lab;    // lab mis a zero par default
        ind_np[np] = ii;
        ind_label[np] = label_point[ii];
        gtree->Add(v[np++]);
      } else {
        ind_label[pvi - v] = min(ind_label[pvi - v], label_point[ii]);
      }
    }

    int ind_multiple[np];

    for (int ii = 0; ii < np; ii++) {
      ind_multiple[ii] = -1;
    }

    for (int ii = 0; ii < NbPoints; ii++) {
      const R3 r3vi(Coord_Point[ii][0], Coord_Point[ii][1], Coord_Point[ii][2]);
      // int label =  label_point[ii];
      const Vertex3 &vi(r3vi);
      Vertex3 *pvi = gtree->ToClose(vi, hseuil);
      ind_multiple[pvi - v] = ind_multiple[pvi - v] + 1;
    }

    int jnp;
    jnp = 0;

    for (int ii = 0; ii < np; ii++) {
      if (ind_multiple[ii] == 0) {
        assert(jnp <= ii);
        ind_np[jnp] = ind_np[ii];
        ind_label[jnp] = ind_label[ii];
        jnp++;
      }
    }

    np = jnp;
  }

  if (int_point_confondus_ok != 0 && int_point_confondus_ok != 1) {
    cout << " point_confondus_ok dans fonction PointCommun_hcode vaut 1 ou 0." << endl;
    exit(1);
  }

  delete[] v;
  delete gtree;
}
