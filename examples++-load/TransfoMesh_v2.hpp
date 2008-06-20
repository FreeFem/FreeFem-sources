#ifndef TRANSFOMESH_V2_HPP_
#define TRANSFOMESH_V2_HPP_

#include<iostream>
#include<cmath>
#include<cassert>
#include <map>
//#include "R3.hpp"
#include "Mesh2dn.hpp"
#include "Mesh3dn.hpp"
#include "tetgen.h"

using namespace std;

void BuildBoundMinDist_th2(  const double *tab_XX, const double *tab_YY, const double *tab_ZZ, const Mesh  & Th2, R3 &bmin, R3 &bmax, double &hmin);
void BuildBoundMinDist_th3(  const double *tab_XX, const double *tab_YY, const double *tab_ZZ, const Mesh3 & Th3, R3 &bmin, R3 &bmax, double &hmin);

//void PointCommun_hcode( const int &dim, const int &NbPoints, const int &point_confondus_ok,double **Coord_Point, int * ind_np, int & np);
void PointCommun_hcode( const int &dim, const int &NbPoints, const int &point_confondus_ok,double **Coord_Point, 
const double *bmin, const double *bmax, const double hmin, int * ind_np, int & np);

void PointCommun_hcode_gtree( const int &dim, const int &NbPoints, const int &point_confondus_ok, 
	double **Coord_Point, const int * label_point,
	const R3 & bmin, const R3 & bmax, const double &hmin, int * ind_np, int * ind_label, int & np);
//void OrderVertexTransfo_base(const double *tab_XX, const double *tab_YY, const double *tab_ZZ, const Mesh3 & Th3, int *Numero_Som, int * ind_nv_t, int & nv_t );

void OrderVertexTransfo_hcode_nv(const int &tab_nv,const double *tab_XX, const double *tab_YY, const double *tab_ZZ, 
const double *bmin, const double *bmax, const double hmin, int *Numero_Som, int * ind_nv_t, int & nv_t );

void OrderVertexTransfo_hcode_nv_gtree( const int & tab_nv, const R3 &bmin, const R3 &bmax, const double & hmin, 
    const double *tab_XX, const double *tab_YY, const double *tab_ZZ, int *Numero_Som, int * ind_nv_t, int & nv_t);

void SamePointElement( const double *tab_XX, const double *tab_YY, const double *tab_ZZ, const Mesh3 & Th3, 
	int &recollement_element, int &recollement_border, int &point_confondus_ok,
	int *Numero_Som, int *ind_nv_t, int *ind_nt_t, int *ind_nbe_t, int *label_nt_t, int *label_nbe_t, int & nv_t, int & nt_t,int & nbe_t );

Mesh3 * Transfo_Mesh3(const Mesh3 & Th3, const double *tab_XX, const double *tab_YY, const double *tab_ZZ, int &border_only, 
	int &recollement_element, int &recollement_border, int &point_confondus_ok);

// fonction pour le cas 2D

void SamePointElement_Mesh2( const double *tab_XX, const double *tab_YY, const double *tab_ZZ, const Mesh & Th2, 
	int &recollement_border, int &point_confondus_ok, int *Numero_Som, int *ind_nv_t, int *ind_nt_t, int *ind_nbe_t, 
	int *label_nbe_t, int & nv_t, int & nt_t,int & nbe_t );

// Mesh3 * Transfo_Mesh2(const Mesh2 & Th2, const double *tab_XX, const double *tab_YY, const double *tab_ZZ, int &border_only);

void Transfo_Mesh2_map_face(const Mesh &Th2, map<int, int> &maptri );

Mesh3 * Transfo_Mesh2_tetgen(const Mesh & Th2, const double *tab_XX, const double *tab_YY, const double *tab_ZZ, 
		int &border_only, int &recollement_border, int &point_confondus_ok, 
		const int &label_tet,const map<int, int> &maptri );

void mesh3_tetgenio_out(const tetgenio &out, const int &label_tet, Mesh3 & Th3);

#endif


