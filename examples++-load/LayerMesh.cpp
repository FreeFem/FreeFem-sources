#include <fstream>
#include <iostream>
#include <cstring>
#include "libmesh5.h"
#include "ufunction.hpp"
#include "error.hpp"
#include "RNM.hpp"
#include<stdlib.h>
#include "ufunction.hpp"
#include "Mesh2dn.hpp"
#include "Mesh2dn.hpp"
#include "Mesh3dn.hpp"
#include "rgraph.hpp"
#include "RNM.hpp"
#include "fem.hpp"
using namespace std;
using namespace Fem2D;
#include "LayerMesh.hpp"



// valeur absolue
// chercher inferieure ou egale + pi
// exemple : C0(20): 20 intervalles ou 20 points de discrétisations
// puissance


// remarque choix 2 est a encore a determiner
double  zmin_func_mesh( const int choix, const double x, const double y  )
{
  
  switch(choix){
  case 0: 
    return 0.;
    break;
  case 1:
    return 0.;   
    break;
  case 2:
    return sqrt(pow(x,2)+pow(y,2));
    break;
  default :
    cout << "zmin_func pas définis" << endl;
    return 0.;
  }
}

double  zmax_func_mesh( const int choix, const double x, const double y  ){
  
  switch(choix){
  case 0: 
    return 1.;
    break;
  case 1:
    return 1.;   
    break;
  case 2:
    return 3.+sqrt(pow(x,2)+pow(y,2));
    break;
  default :
    cout << "zmaxfunc pas définis" << endl;
    return 0.;
  }
}

int Ni_func_mesh( const int choix, const double x, const double y  ){
  const int multi=1;
  int res;
  switch(choix){
  case 0:
    if( x==0. && y==0.){
      res = 3;
    }
    if( x==1. && y==0.){
      res = 5;
    }
    if( x==0. && y==1.){
      res = 7;
    }
    if( x==0.5 && y==0.5){
      res = 6;
    }
    return res ;
    //return multi;
    break;
  case 1:
    return 2;   
    break;
  case 2:
    return int(multi*(3+sqrt(pow(x,2)+pow(y,2))));
    break;
  default :
    cout << "Ni_func pas définis" << endl;
    return 0;
  }
}

void discretisation_max_mesh(const int choix,  const Mesh & Th2, int & Nmax){  
  int Ni;

  Nmax = 0;  
  /*for(int ii=0; ii < A2D.NbSommet2D;ii++){
    Ni   = Ni_func( choix, A2D.CoorSommet2D[ii][0], A2D.CoorSommet2D[ii][1]); 
  	Nmax = max(Ni,Nmax);
  }
  Nmax=4;*/
  for(int ii=0; ii < Th2.nv; ii++){
     const  Vertex2 & P = Th2.vertices[ii];
     Ni   = Ni_func_mesh( choix, P.x, P.y); 
  	 Nmax = max(Ni,Nmax);
  }
}

void tab_zmin_zmax_Ni_mesh(const int choix, const Mesh & Th2, int & Nmax,double *tab_zmin, double *tab_zmax,int *tab_Ni){
   Nmax = 0;	
   for(int ii=0; ii < Th2.nv; ii++){
	 const  Vertex2 & P = Th2.vertices[ii];
     tab_Ni[ii] = Ni_func_mesh( choix, P.x, P.y);
	 tab_zmin[ii] = zmin_func_mesh( choix, P.x, P.y );
     tab_zmax[ii] = zmax_func_mesh( choix, P.x, P.y );    
     Nmax = max(tab_Ni[ii],Nmax);
   }
}
/* Fonction permettant de transformer maillage 2D en maillage 3D*/

void Tet_mesh3_mes_neg(Mesh3 & Th3){
    int iv[4];
    int lab;
    
	for(int ii=0; ii< Th3.nt; ii++){
	  const Tet & K(Th3.t(ii));
 	lab   = K.lab;
 
    iv[0] = Th3.operator()(K[0]);
    iv[2] = Th3.operator()(K[1]);
    iv[1] = Th3.operator()(K[2]);
    iv[3] = Th3.operator()(K[3]);
    R3 A(Th3.vertices[iv[0]]);
    R3 B(Th3.vertices[iv[1]]);
    R3 C(Th3.vertices[iv[2]]);
    R3 D(Th3.vertices[iv[3]]);
    double mes=det(A,B,C,D)/6.;
    Th3.t(ii).set(Th3.vertices, iv, lab,mes);	
	}	
}



Mesh3 * build_layer (const Mesh & Th2, const int Nmax, const int *tab_Ni,
		     const double *tab_zmin, const double *tab_zmax){
  int MajSom, MajElem, MajBord2D;     
  Mesh3 *Th3=new Mesh3;
  NbSom3D_NbElem3D_NbBord2D_mesh_product_mesh_tab( Nmax, tab_Ni, Th2, MajSom, MajElem, MajBord2D);   
  cout << "MajSom = " <<  MajSom  << "  "  << "MajElem = " <<  MajElem  << " " << "MajBord2D =" << MajBord2D << endl;
  
  cout << "debut :   Th3.set(MajSom, MajElem, MajBord2D);     "<< endl;
  Th3->set(MajSom,MajElem,MajBord2D);
  
  cout << "debut :   Som3D_mesh_product_Version_Sommet_mesh_tab( Nmax, tab_Ni, tab_zmin, tab_zmax, Th2, Th3);   "<< endl;
  Som3D_mesh_product_Version_Sommet_mesh_tab( Nmax, tab_Ni, tab_zmin, tab_zmax, Th2, *Th3);    
  return Th3;
}

void NbSom3D_NbElem3D_NbBord2D_mesh_product_mesh_tab(const int Nmax, const int *tab_Ni, const Mesh &Th2,  int &MajSom, int &MajElem, int &MajBord2D){  
  int i;
  
  MajSom = 0;  
  for(int ii=0; ii < Th2.nv;ii++){  
    MajSom = MajSom + (tab_Ni[ii]+1);   
      
  	assert(tab_Ni[ii]<=Nmax);
  }

  MajElem = 0;  
  for(int ii=0; ii < Th2.nt; ii++){   
    const Mesh::Triangle & K(Th2.t(ii));
    for(int jj=0; jj < 3; jj++){ 
      //i  = A2D.ElemPoint2D[ii][jj]; 
      i  = Th2.operator()(K[jj]); 
      MajElem = MajElem + tab_Ni[i]; 
    }
  }
  
  // determination of NbBord2D
  MajBord2D = 2*Th2.nt;  
 
  for(int ii=0; ii < Th2.neb;ii++)
    {
      const Mesh::BorderElement  & K(Th2.be(ii));
      for(int jj=0; jj < 2; jj++)
	{  
	  // i  = A2D.ElemBord1D[ii][jj];
	  i=Th2(K[jj]); 
	  MajBord2D = MajBord2D + tab_Ni[i];
	  assert( tab_Ni[i] <= Nmax);	
	}
    }
}

void Som3D_mesh_product_Version_Sommet_mesh_tab(const int Nmax, const int *tab_Ni, const double *tab_zmin, const double *tab_zmax, const Mesh &Th2, Mesh3 & Th3){
  // intent(in)  Nmax,Mesh &A2D
  // intent(out) Mesh3 &A3D
  
  double val_zmin,val_zmax,val_dz;
  int    Ni;
  int    NumSommet;
  int    NumElement;
  int *  tab_NumSommet= new int[Th2.nv+1];

  // variable tet
  int    SommetPrisme[6];
 
  // variable creer pour le bord
  int    i_ind1,Ni_ind1;
  int    i_ind2,Ni_ind2; 
  int    i_recoll_1pp,i_recoll_2pp; 
  int    i_recoll_1, i_recoll_2; 
  //int    pas_recoll_1, pas_recoll_2; 
  int    type_dec_border; 
  
  // avec data
  int i_recoll_jMax,i_recoll_jMaxpp;
  int cas_decoupage; //, cas_data; 
  int int_decoup[3] = {1,2,4};
  int Ni_elem[3]; 
  int DiagMax1,DiagMax2;

 
  NumSommet = 0;

  for( int ii=0; ii < Th2.nv; ii++){
    const  Vertex2 & P = Th2.vertices[ii];

    val_zmin = tab_zmin[ii];
    val_zmax = tab_zmax[ii];    
    Ni       = tab_Ni[ii];

   
    val_dz = (val_zmax - val_zmin)/Ni;
    if( Ni == 0) val_dz = 0.;

    tab_NumSommet[ii] = NumSommet; // Numero du premier sommet 3D associé au sommet 2D ii.
	//cout << "ii, tab_NumSommet[ii]= "<< ii <<" "<< tab_NumSommet[ii] << endl;
    for(int j=0; j <= Ni; j++){ //changer
      Th3.vertices[NumSommet].x = P.x; 
      Th3.vertices[NumSommet].y = P.y; 
      Th3.vertices[NumSommet].z = val_zmin + val_dz*j;

      Th3.vertices[NumSommet].lab = P.lab; 
      if(j==0)  Th3.vertices[NumSommet].lab = 5;
      if(j==Ni) Th3.vertices[NumSommet].lab = 6;      
    
      NumSommet = NumSommet+1;
    }
   	
  }
  tab_NumSommet[Th2.nv] = NumSommet;

  assert( NumSommet == Th3.nv );
  
  
  //=======================================================================
  // creation des bord du maillage 3D a partir du bord 1D et du maillage 2D
  //=======================================================================

  cout << "calcul element du bord " << endl;
 
  // A mettre plus haut
  int ElemBord;

  ElemBord = 0;

  // bord définies en zmax
  
  for(int ii=0; ii < Th2.nt; ii++){
    int ijj[3];
    const Mesh::Element & K(Th2.t(ii));
    int lab=6;
    
    for(int kk=0; kk < 3; kk++){
	  ijj[kk] = Th2.operator()(K[kk]);  
      ijj[kk] = tab_NumSommet[ijj[kk]+1]-1;
	}
        
    Th3.be(ElemBord).set(Th3.vertices,ijj,lab);
    
    ElemBord = ElemBord+1;
  }

  // bord definis en zmin

  for(int ii=0; ii < Th2.nt; ii++){
    int ijj[3];
    const Mesh::Element & K(Th2.t(ii));
    int lab=5;
    
    for(int kk=0; kk < 3; kk++){
	  ijj[2-kk] = Th2.operator()(K[kk]); 
      ijj[2-kk] = tab_NumSommet[ijj[2-kk]];
	}
	
    Th3.be(ElemBord).set(Th3.vertices,ijj,lab);

    ElemBord = ElemBord+1;
  }
  
  cout << "bord sur le cote" << endl;

  for(int ii=0; ii < Th2.neb; ii++){
    int ijj[3];
	
    const Mesh::BorderElement & K(Th2.be(ii));
	int lab = K.lab;
	
	/*int lab = 0;*/

    i_ind1  = Th2.operator()(K[0]);
    i_ind2  = Th2.operator()(K[1]);

    Ni_ind1 =  tab_Ni[i_ind1];  
    Ni_ind2 =  tab_Ni[i_ind2]; 
	  
	assert( Ni_ind1 <= Nmax);
	assert( Ni_ind2 <= Nmax);
	
    for(int jNmax=Nmax-1; jNmax >=0; jNmax--){
	    
      /*
      i_recoll_1pp = int((jNmax+1)*Ni_ind1/Nmax);
      i_recoll_2pp = int((jNmax+1)*Ni_ind2/Nmax); 
      
      i_recoll_1 = int(jNmax*Ni_ind1/Nmax);
      i_recoll_2 = int(jNmax*Ni_ind2/Nmax); 
      */
      
      i_recoll_1 = int((jNmax+1)*Ni_ind1/Nmax);
      i_recoll_2 = int((jNmax+1)*Ni_ind2/Nmax); 
      
      i_recoll_1pp = int(jNmax*Ni_ind1/Nmax);
      i_recoll_2pp = int(jNmax*Ni_ind2/Nmax);
      
      /*
	type_dec_border = 0  tous les points sont confondus
	type_dec_border = 1  les points 1pp et 1 sont differents
	type_dec_border = 2  les points 2pp et 2 sont differents
	type_dec_border = 3  les points 1pp et 1 et les points 2pp et 2 sont differents 

	rappel :  1pp(0) 2pp(1) 2(2) 1(3) 
	data_dec_border 1 :   {3 1 0}
	data_dec_border 2 :   {2 1 0}
	data_dec_border 3 : type1 : { {2 1 0}{0 3 2} }
	                  : type2 : { {3 1 0}{1 3 2} }
      */
      type_dec_border = 0; 
      
      if( i_recoll_1pp != i_recoll_1){ 
	type_dec_border = type_dec_border + 1;
      }
      
      if( i_recoll_2pp != i_recoll_2){ 
	type_dec_border = type_dec_border + 2;
      }

      //cout << "type decoupage bord= " <<  type_dec_border <<endl;
           
      switch( type_dec_border ){
      case 0:
	// rien n a faire
	break;
      case 1:
	
	ijj[2] = tab_NumSommet[i_ind1]+i_recoll_1pp;
	ijj[1] = tab_NumSommet[i_ind2]+i_recoll_2pp;
	ijj[0] = tab_NumSommet[i_ind1]+i_recoll_1; 
	
	Th3.be(ElemBord).set(Th3.vertices,ijj,lab);
	
	ElemBord = ElemBord+1;
	break;
      case 2:
     
	ijj[2] = tab_NumSommet[i_ind1]+i_recoll_1pp;
	ijj[1] = tab_NumSommet[i_ind2]+i_recoll_2pp;
	ijj[0] = tab_NumSommet[i_ind2]+i_recoll_2; 
	
	Th3.be(ElemBord).set(Th3.vertices,ijj,lab);
	
	ElemBord = ElemBord+1;
	break;
      case 3:
	int idl;

	// determination de la diagonale Max
	DiagMax1 = max( tab_NumSommet[i_ind1]+i_recoll_1pp, tab_NumSommet[i_ind2]+i_recoll_2 );
	DiagMax2 = max( tab_NumSommet[i_ind2]+i_recoll_2pp, tab_NumSommet[i_ind1]+i_recoll_1 );	

	
	if(DiagMax1 > DiagMax2){  
	  idl = 1; 
      
      ijj[2] = tab_NumSommet[i_ind1]+i_recoll_1pp;
	  ijj[1] = tab_NumSommet[i_ind2]+i_recoll_2pp;
	  ijj[0] = tab_NumSommet[i_ind2]+i_recoll_2; 
	
	  Th3.be(ElemBord).set(Th3.vertices,ijj,lab);
	  
	  ijj[2] = tab_NumSommet[i_ind2]+i_recoll_2;
	  ijj[1] = tab_NumSommet[i_ind1]+i_recoll_1;
	  ijj[0] = tab_NumSommet[i_ind1]+i_recoll_1pp; 
	
	  Th3.be(ElemBord+1).set(Th3.vertices,ijj,lab);
	}
	else{
	  idl = 2;
	  
      ijj[2] = tab_NumSommet[i_ind1]+i_recoll_1pp;
	  ijj[1] = tab_NumSommet[i_ind2]+i_recoll_2pp;
	  ijj[0] = tab_NumSommet[i_ind1]+i_recoll_1; 

	  Th3.be(ElemBord).set(Th3.vertices,ijj,lab);
	  
	  ijj[2] = tab_NumSommet[i_ind2]+i_recoll_2;
	  ijj[1] = tab_NumSommet[i_ind1]+i_recoll_1;
	  ijj[0] = tab_NumSommet[i_ind2]+i_recoll_2pp; 
	  
	  Th3.be(ElemBord+1).set(Th3.vertices,ijj,lab);
	}
	//cout << "idl=" << idl << endl; 
	ElemBord = ElemBord+2;
	break;
      default:
	break;
      }     
    }    
  }
  
  assert( ElemBord == Th3.nbe );
  //=========================================
  // Creation + determination tetraedre

  cout << "calcul element tetraedre " << endl;
  
  NumElement =  0;

  for(int ii=0; ii < Th2.nt; ii++){
      /*  
	  nouvelle numerotation : 
	  -----------------------
	  Valeur de cas_deoupage
	  -----------------------
	  1 : sommet 0 et 3 differents
	  2 : sommet 1 et 4 differents
	  4 : sommet 2 et 5 differents
	  ============================
	  3 : sommet 0 et 3 differents + sommet 1 et 4 differents
	  5 : sommet 0 et 3 differents + sommet 2 et 5 differents
	  6 : sommet 1 et 4 differents + sommet 2 et 5 differents
	  ============================
	  7 : aucun sommets confondus

	  data_tetraedre 
	  ==============
	  1: 0++,1++,2++,SomDiff :  {0 1 2 3}        :: data 1
	  2: 0++,1++,2++,SomDiff :  {0 1 2 4}        :: data 2
	  4: 0++,1++,2++,SomDiff :  {0 1 2 5}        :: data 3
	  ==============
	  = deux cas possible depend du sommet le plus grand : Sommet le plus grand est un ++
	  = 0++,1++,2++,SomDiffMin  ||  SomDiffMax++, j_SomDiff_py[j_SomEgal][0], j_SomDiff_py[j_SomEgal][1], Som_Egal 
	  3:a: SommetMax diag 04  {0,1,2,4} {5,4,3,0} :: data 4
	  3:b: SommetMax diag 13  {0,1,2,3} {5,4,3,1} :: data 5
	  =============================================
	  5:a: SommetMax diag 05  {0,1,2,5} {5,4,3,0} :: data 6
	  5:b: SommetMax diag 23  {0,1,2,3} {5,4,3,2} :: data 7
	  =============================================
	  6:a: SommetMax diag 15  {0,1,2,5} {5,4,3,1} :: data 8
	  6:b: SommetMax diag 24  {0,1,2,4} {5,4,3,2} :: data 9
	  =============================================
	  7: aller chercher dans la fonction          :: data 10 a data 
	  == voir hecht routine

      */
    const Mesh::Element & K(Th2.t(ii));
      int somv[4];
      int K_jj[3]; 
      int lab=0;
     
      // valeur de Nombre de points
      for(int jj=0; jj <3; jj++){
        K_jj[jj] = Th2.operator()(K[jj]);   
        Ni_elem[jj] = tab_Ni[ K_jj[jj] ];    
      }
     
    for(int jNmax=Nmax-1; jNmax >=0; jNmax--){
      // determination des sommets + cas decoupage
      cas_decoupage = 0;
      for(int jj=0; jj<3; jj++){
	
	i_recoll_jMax   = int( (jNmax)*Ni_elem[jj]/Nmax );
	i_recoll_jMaxpp = int( (jNmax+1)*Ni_elem[jj]/Nmax );

	SommetPrisme[jj+3]   = tab_NumSommet[ K_jj[jj] ] + i_recoll_jMaxpp;
	SommetPrisme[jj] = tab_NumSommet[ K_jj[jj] ] + i_recoll_jMax;
   
	assert( SommetPrisme[jj]   <= Th3.nv);	
   	assert( SommetPrisme[jj+3] <= Th3.nv);
   	
   	//cout << "Sommet du prisme jj=" << jj << " " << "jj  "<< SommetPrisme[jj] << "  jj+3  " << SommetPrisme[jj+3] << endl;
       
	if( i_recoll_jMax != i_recoll_jMaxpp) cas_decoupage = cas_decoupage + int_decoup[jj];
      }

      //cout << "cas du decoupage= " << cas_decoupage << endl;

      switch( cas_decoupage ){
	
      case 0 : 
	// les points sont tous confondus pas d ajout element : rien a faire
	break;
	/*
	  CAS CREATION D UN TETRAEDRE : cas decoupage 1 2 4
	  
	 */
      case 1 :
	// On a un tetraedre
	
	somv[0] = SommetPrisme[0];
	somv[1] = SommetPrisme[1]; 
	somv[2] = SommetPrisme[2]; 
	somv[3] = SommetPrisme[3];

	Th3.t(NumElement).set(Th3.vertices, somv, lab);	
	
	NumElement = NumElement+1;
	break;
      case 2 :
	// On a un tetraedre
	
	somv[0] = SommetPrisme[0];
	somv[1] = SommetPrisme[1]; 
	somv[2] = SommetPrisme[2]; 
	somv[3] = SommetPrisme[4];

	Th3.t(NumElement).set(Th3.vertices, somv, lab);	
	
	NumElement = NumElement+1;
	break;
      case 4 :
	// On a un tetraedre
	
	somv[0] = SommetPrisme[0];
	somv[1] = SommetPrisme[1]; 
	somv[2] = SommetPrisme[2]; 
	somv[3] = SommetPrisme[5];

	Th3.t(NumElement).set(Th3.vertices, somv, lab);	
	
	NumElement = NumElement+1;
	break;
	/*
	  On a une pyramide a base rectangle: decoupe deux tetraedres
	  cas decoupage 3 5 6
	*/
      case 3 :
	// determination de la diagonale dominante
	DiagMax1 = max( SommetPrisme[0], SommetPrisme[4] );
	DiagMax2 = max( SommetPrisme[1], SommetPrisme[3] );
	
	//cout << "DiagMax1=" << DiagMax1 << " "<< SommetPrisme[0]<<" " <<SommetPrisme[4] << endl;

	if( DiagMax1 > DiagMax2){
	  //------------------
	  // premier tetraedre 
	  somv[0] = SommetPrisme[0];
	  somv[1] = SommetPrisme[1];	
	  somv[2] = SommetPrisme[2];
	  somv[3] = SommetPrisme[4]; 
	  
	  Th3.t(NumElement).set(Th3.vertices, somv, lab);	
	  // deuxieme tetraedre
	  somv[0] = SommetPrisme[5];
	  somv[1] = SommetPrisme[4];	
	  somv[2] = SommetPrisme[3];
	  somv[3] = SommetPrisme[0]; 
	  
	  Th3.elements[NumElement+1].set(Th3.vertices, somv, lab);	
	}			  
	else{
	  //------------------
	  // premier tetraedre 
	  somv[0] = SommetPrisme[0];
	  somv[1] = SommetPrisme[1];	
	  somv[2] = SommetPrisme[2];
	  somv[3] = SommetPrisme[3];
	  
	  Th3.t(NumElement).set(Th3.vertices, somv, lab);	 
	  // deuxieme tetraedre
	  somv[0] = SommetPrisme[5];
	  somv[1] = SommetPrisme[4];	
	  somv[2] = SommetPrisme[3];
	  somv[3] = SommetPrisme[1]; 
	  
	  Th3.elements[NumElement+1].set(Th3.vertices, somv, lab);	
	}
      
	NumElement = NumElement+2;
	break;

      case 5 :
	// determination de la diagonale dominante
	DiagMax1 = max( SommetPrisme[0], SommetPrisme[5] );
	DiagMax2 = max( SommetPrisme[2], SommetPrisme[3] );
	
	//cout << "DiagMax1=" << DiagMax1 << " "<< SommetPrisme[0]<<" " <<SommetPrisme[5] << endl;

	if( DiagMax1 > DiagMax2){
	  //------------------
	  // premier tetraedre 
	  somv[0] = SommetPrisme[0];
	  somv[1] = SommetPrisme[1];	
	  somv[2] = SommetPrisme[2];
	  somv[3] = SommetPrisme[5];
	  
	  Th3.t(NumElement).set(Th3.vertices, somv, lab);	 
	  // deuxieme tetraedre
	  somv[0] = SommetPrisme[5];
	  somv[1] = SommetPrisme[4];	
	  somv[2] = SommetPrisme[3];
	  somv[3] = SommetPrisme[0]; 
	  
	  Th3.elements[NumElement+1].set(Th3.vertices, somv, lab);	
	}			  
	else{
	  //------------------
	  // premier tetraedre 
	  somv[0] = SommetPrisme[0];
	  somv[1] = SommetPrisme[1];	
	  somv[2] = SommetPrisme[2];
	  somv[3] = SommetPrisme[3];
	  
	  Th3.t(NumElement).set(Th3.vertices, somv, lab);	 
	  // deuxieme tetraedre
	  somv[0] = SommetPrisme[5];
	  somv[1] = SommetPrisme[4];	
	  somv[2] = SommetPrisme[3];
	  somv[3] = SommetPrisme[2]; 
	  
	  Th3.elements[NumElement+1].set(Th3.vertices, somv, lab);	
	}
      
	NumElement = NumElement+2;
	break;

      case 6 :
	// determination de la diagonale dominante
	DiagMax1 = max( SommetPrisme[1], SommetPrisme[5] );
	DiagMax2 = max( SommetPrisme[2], SommetPrisme[4] );

	//cout << "DiagMax1=" << DiagMax1 << " "<< SommetPrisme[1]<<" " <<SommetPrisme[5] << endl;
	
	if( DiagMax1 > DiagMax2){
	  //------------------
	  // premier tetraedre 
	  somv[0] = SommetPrisme[0];
	  somv[1] = SommetPrisme[1];	
	  somv[2] = SommetPrisme[2];
	  somv[3] = SommetPrisme[5];
	  
	  Th3.t(NumElement).set(Th3.vertices, somv, lab);	 
	  // deuxieme tetraedre
	  somv[0] = SommetPrisme[5];
	  somv[1] = SommetPrisme[4];	
	  somv[2] = SommetPrisme[3];
	  somv[3] = SommetPrisme[1]; 
	  
	  Th3.elements[NumElement+1].set(Th3.vertices, somv, lab);	
	}			  
	else{
	  //------------------
	  // premier tetraedre 
	  somv[0] = SommetPrisme[0];
	  somv[1] = SommetPrisme[1];	
	  somv[2] = SommetPrisme[2];
	  somv[3] = SommetPrisme[4];
	  
	  Th3.t(NumElement).set(Th3.vertices, somv, lab);	 
	  // deuxieme tetraedre
	  somv[0] = SommetPrisme[5];
	  somv[1] = SommetPrisme[4];	
	  somv[2] = SommetPrisme[3];
	  somv[3] = SommetPrisme[2]; 
	  
	  Th3.elements[NumElement+1].set(Th3.vertices, somv, lab);	
	}
      
	NumElement = NumElement+2;
	break;

      case 7 :
	// on a un prisme 
	int nbe;
	int option=1;
	int idl[3];
	int nu[12];
	
    DiagMax1 = max( SommetPrisme[0], SommetPrisme[5] );
	DiagMax2 = max( SommetPrisme[2], SommetPrisme[3] );	
	
	// determination de idl
	// idl[0] : choix sommet 0 ou 2 (dpent1 equivalent 1 ou 3)
	
	if(DiagMax1 > DiagMax2){  
	  idl[0]=1;
	}
	else{
	  idl[0]=2;
	}
    
	DiagMax1 = max( SommetPrisme[0], SommetPrisme[4] );
	DiagMax2 = max( SommetPrisme[1], SommetPrisme[3] );	
	
	// idl[1] : choix sommet 0 ou 1 (dpent1 equivalent 1 ou 2)	
    if(DiagMax1 > DiagMax2){
	  idl[1]=1;
	}
	else{
	  idl[1]=2;
	}
	
	DiagMax1 = max( SommetPrisme[1], SommetPrisme[5] );
	DiagMax2 = max( SommetPrisme[2], SommetPrisme[4] );	
	
	// idl[2] : choix sommet 1 ou 2 (dpent1 equivalent 2 ou 3)	
	if(DiagMax1 > DiagMax2){
	  idl[2]=1;
	}
	else{
	  idl[2]=2;
	}
	//cout << "idl[0] << << idl[1] << << idl[2]" << endl;
    //cout << idl[0] << " " << idl[1] << "  "<< idl[2] << endl;
	
	nbe = 0;

	dpent1_mesh( idl, nu, nbe, option);
      
	if(nbe!=3){cout << nbe << endl; cerr << "probleme dans dpent1_mesh" << endl;  };
	
	//------------------
	// premier tetraedre 
	somv[0] = SommetPrisme[nu[0]];
	somv[1] = SommetPrisme[nu[1]];	
	somv[2] = SommetPrisme[nu[2]];	
	somv[3] = SommetPrisme[nu[3]];
	
	Th3.t(NumElement).set(Th3.vertices, somv, lab);		
	// deuxieme tetraedre
	somv[0] = SommetPrisme[nu[4]]; 
	somv[1] = SommetPrisme[nu[5]]; 
	somv[2] = SommetPrisme[nu[6]]; 
	somv[3] = SommetPrisme[nu[7]];
	
	Th3.t(NumElement+1).set(Th3.vertices, somv, lab);	 	
	// troisieme tetraedre
	somv[0] = SommetPrisme[nu[8]]; 
	somv[1] = SommetPrisme[nu[9]]; 
	somv[2] = SommetPrisme[nu[10]]; 
	somv[3] = SommetPrisme[nu[11]]; 
	
	Th3.t(NumElement+2).set(Th3.vertices, somv, lab);	
	
	NumElement = NumElement+3;
	break;
      }
      
   }
  // Au final : les sommers des tetraedres et la conectivité des tetraedres finaux
  assert(NumElement <= Th3.nt);
  }
}

void dpent1_mesh(int idl[3],int nu[12],int &nbe,int &option){
  // intent(inout)  :: idl
  // intent(out)    :: nu,nbe,option
  // option ne sert à rien
  //* version simplifie pour le mailleur par couche 2D 3D
  //-----------------------------------------------------------------------
  //      subroutine dpent1 (idl,nu,nbe,option)
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //                s.p. dpent1
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //  but : decoupe un pentaedre en 3 tetreadres suivant la decoupe des 3
  // ---   faces frontieres a 4 cotes
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //  parametres en entre :
  //   idl : parametre de decoupe de face calculer comme ceci :
  //       si idl(i) = 0 alors la face n'est pas  decoupee
  //       idl(1)= 1 si la face 1463 est decoupe par l'arete 16 ,sinon 2
  //       idl(2)= 1 si la face 1254 est decoupe par l'arete 15 ,sinon 2
  //       idl(3)= 1 si la face 2365 est decoupe par l'arete 26 ,sinon 2
  //          id = i1 + i2 * 2 + i3 * 4
  //  parametres en sortie :
  //   nbe : nbe de tetraedre de la decoupe
  //         nbe = 0 => decoup impossible
  //         nbe = 3 => decoup possible le tableau nu est genere
  //   nu(1:4,1:nbe) : tableau de numero des sommet 3 tetraedres dans le
  //                   pentaedre
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //  programmation : f77 ->c++ subroutine de f. hecht upmc 

  int idp[8];
  int i1,i2,i3,i,nbdp,idf,idecou;
  const int pdd[8] = {1,0,2,3,4,5,0,6};

  int mu[6][12];
  const int mu0[12] = {1,6,2,3, 1,5,2,6, 1,6,4,5};
  const int mu1[12] = {1,6,2,3, 1,4,2,6, 2,6,4,5};
  const int mu2[12] = {1,4,2,3, 2,6,3,4, 2,6,4,5};
  const int mu3[12] = {1,5,2,3, 1,5,3,6, 1,6,4,5};
  const int mu4[12] = {1,5,2,3, 1,5,3,4, 3,6,4,5};
  const int mu5[12] = {1,4,2,3, 2,5,3,4, 3,6,4,5};

  for(int jj=0; jj<12; jj++){
    mu[0][jj] = mu0[jj];
    mu[1][jj] = mu1[jj];
    mu[2][jj] = mu2[jj];
    mu[3][jj] = mu3[jj];
    mu[4][jj] = mu4[jj];
    mu[5][jj] = mu5[jj];
  }

  // calcul des descoupes possible du pentaedre
  idf  = -1;
  nbdp =  0;

  for(i3=1; i3<=2; i3++){
    for(i2=1; i2<=2; i2++){
      for(i1=1; i1<=2; i1++){
	idf=idf+1;
	if( (pdd[idf] != 0) 
	    && (  idl[0]==0  || idl[0]==i1 )
	    && (  idl[1]==0  || idl[1]==i2 )
	    && (  idl[2]==0  || idl[2]==i3 ) ){
	  //nbdp=nbdp+1;
	  idp[nbdp]=idf;
	  nbdp=nbdp+1;
	}
      }
    }
  }
  
  if(nbdp == 0){
    nbe=0;
  }
  else{
    nbe=3;
    idf=idp[0];
    idecou=pdd[idf];
   /* i=idf;  
    j=i/4;
    i=i-4*j;
    idl[2]=j+1;
    j=i/2;
    idl[1]=j+1;
    idl[0]=i-2*j+1;
    //cout << "idecou= " << idecou << endl;*/
    for(i=0; i<12;i++){
      nu[i]=mu[idecou-1][i]-1;
      //cout << "i, nu[i] "<< i <<" " << nu[i] << endl;
    }
  }

  }
//----------------------------------------------------------------------- 


/*
void dpent1_mesh(int idl[3],int nu[12],int &nbe,int &option){
  // intent(inout)  :: idl
  // intent(out)    :: nu,nbe,option
  // option ne sert à rien

  //-----------------------------------------------------------------------
  //      subroutine dpent1 (idl,nu,nbe,option)
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //                s.p. dpent1
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //  but : decoupe un pentaedre en 3 tetreadres suivant la decoupe des 3
  // ---   faces frontieres a 4 cotes
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //  parametres en entre :
  //   idl : parametre de decoupe de face calculer comme ceci :
  //       si idl(i) = 0 alors la face n'est pas  decoupee
  //       idl(1)= 1 si la face 1463 est decoupe par l'arete 16 ,sinon 2
  //       idl(2)= 1 si la face 1254 est decoupe par l'arete 15 ,sinon 2
  //       idl(3)= 1 si la face 2365 est decoupe par l'arete 26 ,sinon 2
  //          id = i1 + i2 * 2 + i3 * 4
  //  parametres en sortie :
  //   nbe : nbe de tetraedre de la decoupe
  //         nbe = 0 => decoup impossible
  //         nbe = 3 => decoup possible le tableau nu est genere
  //   nu(1:4,1:nbe) : tableau de numero des sommet 3 tetraedres dans le
  //                   pentaedre
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //  programmation : f77 ->c++ subroutine de f. hecht upmc 

  int idp[8];
  int i1,i2,i3,i,j,nbdp,idf,idecou;
  const int pdd[8] = {1,0,2,3,4,5,0,6};

  int mu[6][12];
  const int mu0[12] = {1,6,2,3, 1,5,2,6, 1,6,4,5};
  const int mu1[12] = {1,6,2,3, 1,4,2,6, 2,6,4,5};
  const int mu2[12] = {1,4,2,3, 2,6,3,4, 2,6,4,5};
  const int mu3[12] = {1,5,2,3, 1,5,3,6, 1,6,4,5};
  const int mu4[12] = {1,5,2,3, 1,5,3,4, 3,6,4,5};
  const int mu5[12] = {1,4,2,3, 2,5,3,4, 3,6,4,5};

  for(int jj=0; jj<12; jj++){
    mu[0][jj] = mu0[jj];
    mu[1][jj] = mu1[jj];
    mu[2][jj] = mu2[jj];
    mu[3][jj] = mu3[jj];
    mu[4][jj] = mu4[jj];
    mu[5][jj] = mu5[jj];
  }

  // calcul des descoupes possible du pentaedre
  idf  = -1;
  nbdp =  0;

  for(i3=1; i3<=2; i3++){
    for(i2=1; i2<=2; i2++){
      for(i1=1; i1<=2; i1++){
	idf=idf+1;
	if( (pdd[idf] != 0) 
	    && (  idl[0]==0  || idl[0]==i1 )
	    && (  idl[1]==0  || idl[1]==i2 )
	    && (  idl[2]==0  || idl[2]==i3 ) ){
	  //nbdp=nbdp+1;
	  idp[nbdp]=idf;
	  nbdp=nbdp+1;
	}
      }
    }
  }
  
  if(nbdp == 0){
    nbe=0;
  }
  else{
    nbe=3;
    idf=idp[0];
    idecou=pdd[idf];
    i=idf;  
    j=i/4;
    i=i-4*j;
    idl[2]=j+1;
    j=i/2;
    idl[1]=j+1;
    idl[0]=i-2*j+1;
    //cout << "idecou= " << idecou << endl;
    for(i=0; i<12;i++){
      nu[i]=mu[idecou-1][i]-1;
      //cout << "i, nu[i] "<< i <<" " << nu[i] << endl;
    }
  }

  }
//----------------------------------------------------------------------- 
*/
