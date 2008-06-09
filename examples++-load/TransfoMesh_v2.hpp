#ifndef TRANSFOMESH_V2_HPP_
#define TRANSFOMESH_V2_HPP_

#include<iostream>
#include<cmath>
#include<cassert>

#include "Mesh2dn.hpp"
#include "Mesh3dn.hpp"

using namespace std;

void PointCommun_hcode( const int dim, const int NbPoints, const int point_confondus_ok,double **Coord_Point, int * ind_np, int & np);

//void OrderVertexTransfo_base(const double *tab_XX, const double *tab_YY, const double *tab_ZZ, const Mesh3 & Th3, int *Numero_Som, int * ind_nv_t, int & nv_t );

void OrderVertexTransfo_hcode_nv(const int tab_nv,const double *tab_XX, const double *tab_YY, const double *tab_ZZ, int *Numero_Som, int * ind_nv_t, int & nv_t );

void SamePointElement( const double *tab_XX, const double *tab_YY, const double *tab_ZZ, const Mesh3 & Th3, 
	int *Numero_Som, int *ind_nv_t, int *ind_nt_t, int *ind_nbe_t, int & nv_t, int & nt_t,int & nbe_t );

Mesh3 * Transfo_Mesh3(const Mesh3 & Th3, const double *tab_XX, const double *tab_YY, const double *tab_ZZ, int &border_only);

void SamePointElement_Mesh2( const double *tab_XX, const double *tab_YY, const double *tab_ZZ, const Mesh2 & Th2, 
	int *Numero_Som, int *ind_nv_t, int *ind_nt_t, int *ind_nbe_t, int & nv_t, int & nt_t,int & nbe_t );

void Transfo_Mesh2(const Mesh2 & Th2, const double *tab_XX, const double *tab_YY, const double *tab_ZZ, int &border_only, Mesh3 & T_Th3 );


#endif


