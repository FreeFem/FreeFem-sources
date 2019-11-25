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

#ifndef TRANSFOMESH_V2_HPP_
#define TRANSFOMESH_V2_HPP_

#include <iostream>
#include <cmath>
#include <cassert>
#include <map>
#include "Mesh2dn.hpp"
#include "Mesh3dn.hpp"

using namespace std;

void BuildBoundMinDist_th2(const double &precis_mesh, const double *tab_XX, const double *tab_YY,
                           const double *tab_ZZ, const Mesh &Th2, R3 &bmin, R3 &bmax, double &hmin);
void BuildBoundMinDist_th3(const double &precis_mesh, const double *tab_XX, const double *tab_YY,
                           const double *tab_ZZ, const Mesh3 &Th3, R3 &bmin, R3 &bmax,
                           double &hmin);

void PointCommun_hcode(const int &dim, const int &NbPoints, const int &point_confondus_ok,
                       double **Coord_Point, const double *bmin, const double *bmax,
                       const double hmin, int *ind_np, int &np);

void PointCommun_hcode_gtree(const int &dim, const int &NbPoints, const int &point_confondus_ok,
                             double **Coord_Point, const int *label_point, const R3 &bmin,
                             const R3 &bmax, const double &hmin, int *ind_np, int *ind_label,
                             int &np);

void OrderVertexTransfo_hcode_nv(const int &tab_nv, const double *tab_XX, const double *tab_YY,
                                 const double *tab_ZZ, const double *bmin, const double *bmax,
                                 const double hmin, int *Numero_Som, int *ind_nv_t, int &nv_t);

void OrderVertexTransfo_hcode_nv_gtree(const int &tab_nv, const R3 &bmin, const R3 &bmax,
                                       const double &hmin, const double *tab_XX,
                                       const double *tab_YY, const double *tab_ZZ, int *Numero_Som,
                                       int *ind_nv_t, int &nv_t);

void SamePointElement(const double &precis_mesh, const double *tab_XX, const double *tab_YY,
                      const double *tab_ZZ, const Mesh3 &Th3, int &recollement_element,
                      int &recollement_border, int &point_confondus_ok, int *Numero_Som,
                      int *ind_nv_t, int *ind_nt_t, int *ind_nbe_t, int *label_nt_t,
                      int *label_nbe_t, int &nv_t, int &nt_t, int &nbe_t);

Mesh3 *Transfo_Mesh3(const double &precis_mesh, const Mesh3 &Th3, const double *tab_XX,
                     const double *tab_YY, const double *tab_ZZ, int &border_only,
                     int &recollement_element, int &recollement_border, int &point_confondus_ok);

// CAS 3D surfacique

void SamePointElement_surf(const double &precis_mesh, const double *tab_XX, const double *tab_YY,
                           const double *tab_ZZ, const Mesh3 &Th3, int &recollement_border,
                           int &point_confondus_ok, int *Numero_Som, int *ind_nv_t, int *ind_nbe_t,
                           int *label_nbe_t, int &nv_t, int &nbe_t);

Mesh3 *Transfo_Mesh3_surf(const double &precis_mesh, const Mesh3 &Th3, const double *tab_XX,
                          const double *tab_YY, const double *tab_ZZ, int &recollement_border,
                          int &point_confondus_ok);

// fonction pour le cas 2D

void SamePointElement_Mesh2(const double &precis_mesh, const double *tab_XX, const double *tab_YY,
                            const double *tab_ZZ, const Mesh &Th2, int &recollement_border,
                            int &point_confondus_ok, int *Numero_Som, int *ind_nv_t, int *ind_nt_t,
                            int *ind_nbe_t, int *label_nbe_t, int &nv_t, int &nt_t, int &nbe_t);

void Transfo_Mesh2_map_face(const Mesh &Th2, map< int, int > &maptri);

Mesh3 *MoveMesh2_func(const double &precis_mesh, const Mesh &Th2, const double *tab_XX,
                      const double *tab_YY, const double *tab_ZZ, int &border_only,
                      int &recollement_border, int &point_confondus_ok);

#endif
