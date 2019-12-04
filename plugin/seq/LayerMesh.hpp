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

#ifndef LAYERMESH_HPP_
#define LAYERMESH_HPP_

#include "map"
using namespace std;

void recollement_maillage_mesh(const int Nmax, const int N, int *tab_recollement);
void dpent1_mesh(int idl[3], int nu[12], int &nbe, int &option);

double zmin_func_mesh(const int choix, const double x, const double y);
double zmax_func_mesh(const int choix, const double x, const double y);
int Ni_func_mesh(const int choix, const double x, const double y);

void tab_zmin_zmax_Ni_mesh(const int choix, const Mesh &Th2, int &Nmax, double *tab_zmin,
                           double *tab_zmax, int *tab_Ni);

void Tet_mesh3_mes_neg(Mesh3 &Th3);

void build_layer_map_tetrahedra(const Mesh &Th2, map< int, int > &maptet);
void build_layer_map_triangle(const Mesh &Th2, map< int, int > &maptrimil,
                              map< int, int > &maptrizmax, map< int, int > &maptrizmin);
void build_layer_map_edge(const Mesh &Th2, map< int, int > &mapemil, map< int, int > &mapezmax,
                          map< int, int > &mapezmin);

void NbSom3D_NbElem3D_NbBord2D_mesh_product_mesh_tab(const int Nmax, const int *tab_Ni,
                                                     const Mesh &Th, int &MajSom, int &MajElem,
                                                     int &MajBord2D);

void Som3D_mesh_product_Version_Sommet_mesh_tab(
  const int Nmax, const int *tab_Ni, const double *tab_zmin, const double *tab_zmax,
  const Mesh &Th2, const map< int, int > &maptet, const map< int, int > &maptrimil,
  const map< int, int > &maptrizmax, const map< int, int > &maptrizmin,
  const map< int, int > &mapemil, const map< int, int > &mapezmax, const map< int, int > &mapezmin,
  Mesh3 &Th3);

void transformation_2D_3D_maillage_mesh_tab(const Mesh &Th2, const int Nmax, const int *tab_Ni,
                                            const double *tab_zmin, const double *tab_zmax,
                                            Mesh3 &Th3);

Mesh3 *build_layer(const Mesh &Th2, const int Nmax, const int *tab_Ni, const double *tab_zmin,
                   const double *tab_zmax, const map< int, int > &maptet,
                   const map< int, int > &maptrimil, const map< int, int > &maptrizmax,
                   const map< int, int > &maptrizmin, const map< int, int > &mapemil,
                   const map< int, int > &mapezmax, const map< int, int > &mapezmin);

#endif
