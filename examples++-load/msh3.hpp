#ifndef MSH3_HPP_
#define MSH3_HPP_


#include<iostream>
#include<cmath>
#include<cassert>
#include <map>
#include "Mesh2dn.hpp"
#include "Mesh3dn.hpp"

using namespace std;

/* avant TransfoMesh_v2.cpp */
void BuildBoundMinDist_th2( const double &precis_mesh, const double *tab_XX, const double *tab_YY, const double *tab_ZZ, const Mesh  & Th2, R3 &bmin, R3 &bmax, double &hmin);
void BuildBoundMinDist_th3( const double &precis_mesh, const double *tab_XX, const double *tab_YY, const double *tab_ZZ, const Mesh3 & Th3, R3 &bmin, R3 &bmax, double &hmin);

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

void SamePointElement( const double &precis_mesh, const double *tab_XX, const double *tab_YY, const double *tab_ZZ, const Mesh3 & Th3, 
	int &recollement_element, int &recollement_border, int &point_confondus_ok,
	int *Numero_Som, int *ind_nv_t, int *ind_nt_t, int *ind_nbe_t, int *label_nt_t, int *label_nbe_t, int & nv_t, int & nt_t,int & nbe_t );

Mesh3 * Transfo_Mesh3(const double &precis_mesh, const Mesh3 & Th3, const double *tab_XX, const double *tab_YY, const double *tab_ZZ, int &border_only, 
	int &recollement_element, int &recollement_border, int &point_confondus_ok);

// CAS 3D surfacique

void SamePointElement_surf( const double &precis_mesh, const double *tab_XX, const double *tab_YY, const double *tab_ZZ, const Mesh3 & Th3, 
	int &recollement_border, int &point_confondus_ok, int *Numero_Som, 
	int *ind_nv_t, int *ind_nbe_t, int *label_nbe_t, int & nv_t,int & nbe_t);

Mesh3 * Transfo_Mesh3_surf(const double &precis_mesh, const Mesh3 & Th3, const double *tab_XX, const double *tab_YY, const double *tab_ZZ, 
	int &recollement_border, int &point_confondus_ok);

	
// fonction pour le cas 2D

void SamePointElement_Mesh2( const double &precis_mesh, const double *tab_XX, const double *tab_YY, const double *tab_ZZ, const Mesh & Th2, 
	int &recollement_border, int &point_confondus_ok, int *Numero_Som, int *ind_nv_t, int *ind_nt_t, int *ind_nbe_t, 
	int *label_nbe_t, int & nv_t, int & nt_t,int & nbe_t );

// Mesh3 * Transfo_Mesh2(const Mesh2 & Th2, const double *tab_XX, const double *tab_YY, const double *tab_ZZ, int &border_only);

void Transfo_Mesh2_map_face(const Mesh &Th2, map<int, int> &maptri );

Mesh3 * MoveMesh2_func( const double &precis_mesh, const Mesh & Th2, const double *tab_XX, const double *tab_YY, const double *tab_ZZ, 
	int &border_only, int &recollement_border, int &point_confondus_ok);


/* avant LayerMesh.cpp */
void recollement_maillage_mesh(const int Nmax, const int N, int *tab_recollement);
void dpent1_mesh(int idl[3],int nu[12],int &nbe,int &option);

double  zmin_func_mesh( const int choix, const double x, const double y  );
double  zmax_func_mesh( const int choix, const double x, const double y  );
int       Ni_func_mesh( const int choix, const double x, const double y  );

void tab_zmin_zmax_Ni_mesh(const int choix, const Mesh & Th2, int & Nmax,double *tab_zmin, double *tab_zmax,int *tab_Ni);

void Tet_mesh3_mes_neg(Mesh3 & Th3);

void build_layer_map_tetrahedra(const Mesh &Th2, map<int, int> &maptet );
void build_layer_map_triangle(const Mesh &Th2, map<int, int> &maptrimil, map<int, int> &maptrizmax, map<int, int> &maptrizmin );
void build_layer_map_edge(const Mesh &Th2, map<int, int> &mapemil, map<int, int> &mapezmax, map<int, int> &mapezmin );

void NbSom3D_NbElem3D_NbBord2D_mesh_product_mesh_tab(const int Nmax, const int *tab_Ni, 
						const Mesh &Th,  int &MajSom, int &MajElem, int &MajBord2D);

void Som3D_mesh_product_Version_Sommet_mesh_tab(const int Nmax, const int *tab_Ni, const double *tab_zmin, 
						const double *tab_zmax, const Mesh &Th2, 
						const map<int, int> &maptet, 
						const map<int, int> &maptrimil, const map<int, int> &maptrizmax, const map<int, int> &maptrizmin, 
						const map<int, int> &mapemil, const map<int, int> &mapezmax, const map<int, int> &mapezmin,  Mesh3 & Th3);

void transformation_2D_3D_maillage_mesh_tab(const Mesh & Th2, const int Nmax, const int *tab_Ni, const double *tab_zmin, const double *tab_zmax, Mesh3 & Th3);

Mesh3 * build_layer (const Mesh & Th2, const int Nmax, const int *tab_Ni, const double *tab_zmin, const double *tab_zmax,
			const map<int, int> &maptet, 
			const map<int, int> &maptrimil, const map<int, int> &maptrizmax, const map<int, int> &maptrizmin, 
			const map<int, int> &mapemil, const map<int, int> &mapezmax, const map<int, int> &mapezmin );


// add 04/01/09
void TestSameVertexMesh3( const Mesh3 & Th3, const double & hseuil, const R3 & Psup, const R3 &Pinf, int & nv_t, int *Numero_Som);
void TestSameTetrahedraMesh3( const Mesh3 & Th3, const double & hseuil, const R3 & Psup, const R3 &Pinf, int & nt_t );
void TestSameTetrahedraMesh3( const Mesh3 & Th3, const double & hseuil, const R3 & Psup, const R3 &Pinf, int *Elem_ok, int & nt_t );
void TestSameTriangleMesh3( const Mesh3 & Th3, const double & hseuil, const R3 & Psup, const R3 &Pinf, int & nbe_t );
void TestSameTriangleMesh3( const Mesh3 & Th3, const double & hseuil, const R3 & Psup, const R3 &Pinf, int *Border_ok ,int & nbe_t );
int TestElementMesh3( const Mesh3 & Th3 );
Mesh3 *TestElementMesh3_patch( const Mesh3 & Th3 );


#endif
