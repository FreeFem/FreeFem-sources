/*#include<iostream>
#include<fstream>
#include<cmath>
#include<cassert>
#include<stdlib.h>
#include "Mesh2dn.hpp"
#include "Mesh3dn.hpp"
#include "TransfoMesh_v2.hpp"
using namespace std;*/
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
#include "TransfoMesh_v2.hpp"



Mesh3 * Transfo_Mesh3(const Mesh3 & Th3, const double *tab_XX, const double *tab_YY, const double *tab_ZZ, int &border_only){
	// cas besoin memoire important
	Mesh3 *T_Th3=new Mesh3;
	int nv_t,nt_t,nbe_t;
	int* Numero_Som;
	
	int* ind_nv_t;
	int* ind_nt_t;
	int* ind_nbe_t;
	
	
	int i_som,i_elem, i_border;	
	
	
	
	Numero_Som = new int[Th3.nv];
	
	ind_nv_t   = new int[Th3.nv];
	ind_nt_t   = new int[Th3.nt];
    ind_nbe_t  = new int[Th3.nbe];
	
       
    cout << "Vertex, Tetrahedra, Border : "<<Th3.nv << ", "<<Th3.nt<< ", " << Th3.nbe<< endl;
    
	for(int ii=0; ii<Th3.nv; ii++){
		Numero_Som[ii]=ii;
	}
	/*
	cout << tab_XX[9] << " "<< tab_YY[9] << " "<< tab_ZZ[9] << endl;
		
	cout << "Th:" << Th3.vertices[9].x << " "<< Th3.vertices[9].y  << " "<< Th3.vertices[9].z << endl;
		
	cout << tab_XX[139] << " "<< tab_YY[139] << " "<< tab_ZZ[139] << endl;
	cout << "Th:" << Th3.vertices[139].x << " "<< Th3.vertices[139].y  << " "<< Th3.vertices[139].z << endl;
		
	cout << tab_XX[410] << " "<< tab_YY[410] << " "<< tab_ZZ[410] << endl;
	cout << "Th:" << Th3.vertices[410].x << " "<< Th3.vertices[410].y  << " "<< Th3.vertices[410].z << endl;
	exit(-1);
	*/	
	cout <<" debut: SamePointElement " <<endl;
		
	SamePointElement( tab_XX, tab_YY, tab_ZZ, Th3, Numero_Som, ind_nv_t, ind_nt_t, ind_nbe_t, nv_t, nt_t, nbe_t);
	
	cout <<" fin: SamePointElement " <<endl;
	
	
	// set size of Mesh T_Th3 
	T_Th3->set(nv_t,nt_t,nbe_t);

	cout << "Transfo TH3 : Vertex, Tetrahedra, Border : "<< "nv_t="<< nv_t << " nt_t=" << nt_t << " nbe_t=" << nbe_t << endl;
		
	// determination of vertex		
	i_som = 0;
	for(int i=0; i<nv_t; i++){
		
		int & ii = ind_nv_t[i];		
		assert( Numero_Som[ii] == i_som );
		
		const Vertex3 & K(Th3.vertices[ii]);
			 
		T_Th3->vertices[i_som].x = tab_XX[ii];
		T_Th3->vertices[i_som].y = tab_YY[ii];
		T_Th3->vertices[i_som].z = tab_ZZ[ii];
		T_Th3->vertices[i_som].lab = K.lab;
			
		i_som = i_som + 1;		
	}	
	cout << "i_som, nv_t=" <<i_som << " "<<nv_t << endl;
	assert( i_som == nv_t);
	
	
	cout << " Transfo volume elements " << endl;
	// determination of volume elements
	i_elem = 0;
	for( int i=0; i< nt_t; i++){
		int & ii=ind_nt_t[i];
	
		// creation of elements
		
		const Tet & K(Th3.elements[ii]);
		int iv[4];
		int lab=K.lab;
		for(int jj=0; jj <4; jj++){
			iv[jj] = Numero_Som[ Th3.operator()(K[jj]) ]; 
			assert( iv[jj] >= 0 && iv[jj] <= nv_t);  	
		}
		//if(i_elem == 65) cout << "border 5:" << iv[0] << " "<< iv[1] << " " <<iv[2] << endl;
					
		
		T_Th3->elements[i_elem].set(T_Th3->vertices, iv, lab);
		i_elem=i_elem+1;
	} 
	
	assert( i_elem == nt_t);
	
	cout << " Transfo border elements " << endl;
	// determination of border elements
	i_border= 0;
	for( int i=0; i< nbe_t; i++){
		int & ii=ind_nbe_t[i];
	
		// creation of elements
		const Triangle3 & K(Th3.be(ii));
		int iv[3];
		int lab=K.lab;
		for(int jj=0; jj <3; jj++){
			iv[jj] = Numero_Som[ Th3.operator()(K[jj]) ];
			assert( iv[jj] >= 0 && iv[jj] <= nv_t);
		}
		//if(i_border == 24) cout << "border 24:" << iv[0] << " "<< iv[1] << " " <<iv[2] << endl;
			
		T_Th3->be(i_border).set(T_Th3->vertices, iv, lab);
		i_border=i_border+1;
	} 
	assert( i_border == nbe_t);
	
	return T_Th3;
}

void SamePointElement( const double *tab_XX, const double *tab_YY, const double *tab_ZZ, const Mesh3 & Th3, 
	int *Numero_Som, int *ind_nv_t, int *ind_nt_t, int *ind_nbe_t, int & nv_t, int & nt_t,int & nbe_t ){
	
		int Elem_ok, Border_ok;
		int recollement_element=1,recollement_border=1;
		
		
		
		cout << tab_XX[9] << " "<< tab_YY[9] << " "<< tab_ZZ[9] << endl;
		
		cout << "Th:" << Th3.vertices[9].x << " "<< Th3.vertices[9].y  << " "<< Th3.vertices[9].z << endl;
		
		cout << tab_XX[139] << " "<< tab_YY[139] << " "<< tab_ZZ[139] << endl;
		cout << "Th:" << Th3.vertices[139].x << " "<< Th3.vertices[139].y  << " "<< Th3.vertices[139].z << endl;
		
		cout << tab_XX[410] << " "<< tab_YY[410] << " "<< tab_ZZ[410] << endl;
		cout << "Th:" << Th3.vertices[410].x << " "<< Th3.vertices[410].y  << " "<< Th3.vertices[410].z << endl;
		
		
		//cout << "  OrderVertexTransfo_base " <<endl;
		//OrderVertexTransfo_base( tab_XX, tab_YY, tab_ZZ,  Th3, Numero_Som, ind_nv_t, nv_t );
		cout << "  OrderVertexTransfo_hcode " <<endl;
		OrderVertexTransfo_hcode_nv( Th3.nv, tab_XX, tab_YY, tab_ZZ, Numero_Som, ind_nv_t, nv_t );
		
		cout << "fin order vertex: nv_t=" << nv_t << endl;
		
		
		
		/* determination de nt_t et de nbe_t*/ 
		int i_elem, i_border;
		
		i_elem = 0;
		for(int ii=0; ii< Th3.nt; ii++){
			const Tet & K(Th3.elements[ii]);
			int iv[4];			
	
			Elem_ok = 1;
				
			for(int jj=0; jj <4; jj++){
				iv[jj] = Numero_Som[ Th3.operator()(K[jj]) ];   
			}
			
			for(int jj=0; jj<4; jj++){
				for(int kk=jj+1; kk<4; kk++){				
					if( iv[jj]==iv[kk] ){
						//cout << "iv[jj]==iv[kk]" << "  ";
						Elem_ok = 0;
						//cout << "Elem_ok[ii]" << ii << " " << Elem_ok[ii] << endl;
						
					 }
				}
			}
			
			if(Elem_ok==1){
				ind_nt_t[i_elem] = ii;
				i_elem = i_elem + 1;
			}

		} 
		nt_t=i_elem;
		
		if(recollement_element ==1){
			int point_confondus_ok = 0;
			cout << "debut recollement : nt_t= "<< nt_t << endl; 
			
			int np,dim=3;
			int *ind_np = new int [nt_t];
			double **Cdg_t=new double *[nt_t];
			for(int i=0; i<nt_t; i++) Cdg_t[i] = new double[dim];
		
		
			for( int i_elem=0; i_elem< nt_t; i_elem++){
				int & ii=ind_nt_t[i_elem];
				const Tet & K(Th3.elements[ii]);
				int iv[4];
			
				for(int jj=0; jj <4; jj++){
					iv[jj] = Th3.operator()(K[jj]) ;
				}	
	
				Cdg_t[i_elem][0] = ( tab_XX[iv[0]] + tab_XX[iv[1]] + tab_XX[iv[2]] + tab_XX[iv[3]] )/4.;
				Cdg_t[i_elem][1] = ( tab_YY[iv[0]] + tab_YY[iv[1]] + tab_YY[iv[2]] + tab_YY[iv[3]] )/4.;
				Cdg_t[i_elem][2] = ( tab_ZZ[iv[0]] + tab_ZZ[iv[1]] + tab_ZZ[iv[2]] + tab_ZZ[iv[3]] )/4.;	
			}
				
			PointCommun_hcode( dim, nt_t, point_confondus_ok, Cdg_t, ind_np, np);
		
			assert( np <= nt_t );
			
			int *ind_nt_t_tmp= new int [np];
		
			for( int i_elem=0; i_elem< np; i_elem++){
				assert( ind_np[i_elem] >=0 && ind_np[i_elem] <= nt_t );
				ind_nt_t_tmp[ i_elem ] = ind_nt_t[ ind_np[i_elem] ]; 
			}	
		
			for( int i_elem=0; i_elem< np; i_elem++){
				ind_nt_t[ i_elem ] = ind_nt_t_tmp[ i_elem ]; 
			}	
		
			nt_t = np;
			
			//delete [] ind_nbe_t_tmp;
			//delete [] ind_np;
		
			cout << "fin recollement : nt_t= "<< nt_t << endl; 
		
		}	
		
		// determination of border elements
		i_border= 0;
		for( int ii=0; ii< Th3.nbe; ii++){
			Border_ok=1;
		
			const Triangle3 & K(Th3.be(ii));
			int iv[3];
			
			for(int jj=0; jj <3; jj++){
				iv[jj] = Numero_Som[ Th3.operator()(K[jj]) ];
			}
			
			for(int jj=0; jj<3; jj++){
				for(int kk=jj+1; kk<3; kk++){				
					if( iv[jj]==iv[kk] ) Border_ok=0;
				}
			}		
			
			if(Border_ok==1){
				ind_nbe_t[i_border] = ii;
				i_border=i_border+1;
			}
		}
		nbe_t = i_border;
		
		if( recollement_border == 1){
			int point_confondus_ok = 1;
			cout << "debut recollement : nbe_t= "<< nbe_t << endl; 
			
			int np,dim=3;
			int *ind_np = new int [nbe_t];
			double **Cdg_be=new double *[nbe_t];
			for(int i=0; i<nbe_t; i++) Cdg_be[i] = new double[dim];
		
		
			for( int i_border=0; i_border< nbe_t; i_border++){
				
				int & ii=ind_nbe_t[i_border];
				const Triangle3 & K(Th3.be(ii));
				int iv[3];
			
				for(int jj=0; jj <3; jj++){
					iv[jj] = Th3.operator()(K[jj]);
				}	
								
				Cdg_be[i_border][0] = ( tab_XX[iv[0]] + tab_XX[iv[1]] + tab_XX[iv[2]] )/3.; //( Th3.vertices[iv[0]].x + Th3.vertices[iv[1]].x + Th3.vertices[iv[2]].x )/3.;
				Cdg_be[i_border][1] = ( tab_YY[iv[0]] + tab_YY[iv[1]] + tab_YY[iv[2]] )/3.; //( Th3.vertices[iv[0]].y + Th3.vertices[iv[1]].y + Th3.vertices[iv[2]].y )/3.;
				Cdg_be[i_border][2] = ( tab_ZZ[iv[0]] + tab_ZZ[iv[1]] + tab_ZZ[iv[2]] )/3.; //( Th3.vertices[iv[0]].z + Th3.vertices[iv[1]].z + Th3.vertices[iv[2]].z )/3.;		
			}
	
			cout << "appele de PointCommun_hcode" << endl;
			PointCommun_hcode( dim, nbe_t, point_confondus_ok, Cdg_be, ind_np, np);
			cout << "fin appele de PointCommun_hcode" << endl;
			
			assert( np <= nbe_t );
			
			int *ind_nbe_t_tmp= new int [np];
		
			for( int i_border=0; i_border<np; i_border++){
				ind_nbe_t_tmp[ i_border ] = ind_nbe_t[ ind_np[i_border] ]; 
			}	
		
			for( int i_border=0; i_border< np; i_border++){
				ind_nbe_t[ i_border ] = ind_nbe_t_tmp[ i_border ]; 
			}	
		
			nbe_t = np;
			
			//delete [] ind_nbe_t_tmp;
			//delete [] ind_np;
		
			cout << "fin recollement : nbe_t= "<< nbe_t << endl; 
		}
}

/* Creation de maillage 2D */ 
void Transfo_Mesh2(const Mesh2 & Th2, const double *tab_XX, const double *tab_YY, const double *tab_ZZ, int &border_only, Mesh3 & T_Th3 ){
	// cas besoin memoire important
	int nv_t,nt_t,nbe_t;
	int* Numero_Som;
	
	int* ind_nv_t;
	int* ind_nt_t;
	int* ind_nbe_t;
		
	int i_som, i_border;	
	
	Numero_Som = new int[Th2.nv];
	
	ind_nv_t   = new int[Th2.nv];
    ind_nbe_t  = new int[Th2.nt];
	     
    cout << Th2.nv << " "<<Th2.nt<< " " << Th2.nbe<< endl;
  
	for(int ii=0; ii<Th2.nv; ii++){ 
		Numero_Som[ii]=ii;
	}
	
	cout <<" debut: SamePointElement " <<endl;
		
	SamePointElement_Mesh2( tab_XX, tab_YY, tab_ZZ, Th2, Numero_Som, ind_nv_t, ind_nt_t, ind_nbe_t, nv_t, nt_t, nbe_t);
	
	cout <<" fin: SamePointElement " <<endl;
	
	// set size of Mesh T_Th3 
	T_Th3.set(nv_t,nt_t,nbe_t);

	cout << "nv_t="<< nv_t << " nt_t=" << nt_t << " nbe_t=" << nbe_t << endl;
	// determination of vertex	
	i_som = 0;
	for(int i=0; i<nv_t; i++){		
		int & ii = ind_nv_t[i];	assert( Numero_Som[ii] == i_som );
		const Vertex2 & K(Th2.vertices[ii]);
			 
		T_Th3.vertices[i_som].x = tab_XX[ii];
		T_Th3.vertices[i_som].y = tab_YY[ii];
		T_Th3.vertices[i_som].z = tab_ZZ[ii];
		T_Th3.vertices[i_som].lab = K.lab;
			
		i_som = i_som + 1;		
	}	
	cout << "i_som, nv_t=" <<i_som << " "<<nv_t << endl;
	assert( i_som == nv_t);
		
	cout << " pas creation volume elements " << endl;
	
	cout << " Transfo border elements " << endl;
	// determination of border elements
	i_border= 0;
	for( int i=0; i< nbe_t; i++){
		int & ii=ind_nbe_t[i];
		
		// creation of elements
		const Triangle2 & K(Th2.elements[ii]);
		int iv[3];
		int lab=K.lab;
		for(int jj=0; jj <3; jj++){
			iv[jj] = Numero_Som[ Th2.operator()(K[jj]) ];
		}
		
		T_Th3.be(i_border).set(T_Th3.vertices, iv, lab);
		i_border=i_border+1;
	} 
	assert( i_border == nbe_t);
	
}

void SamePointElement_Mesh2( const double *tab_XX, const double *tab_YY, const double *tab_ZZ, const Mesh2 & Th2, 
	int *Numero_Som, int *ind_nv_t, int *ind_nt_t, int *ind_nbe_t, int & nv_t, int & nt_t,int & nbe_t ){
	
		int Border_ok;
		int recollement_border=1;
	
		cout << "debut: OrderVertexTransfo_hcode " <<endl;
		OrderVertexTransfo_hcode_nv( Th2.nv, tab_XX, tab_YY, tab_ZZ, Numero_Som, ind_nv_t, nv_t );
		cout << "fin order vertex: nv_t=" << nv_t << endl;
				
		/* determination de nt_t et de nbe_t*/ 
		nt_t = 0;
		int i_border;
		
		// determination of border elements
		i_border= 0;
		for( int ii=0; ii< Th2.nt; ii++){
			Border_ok=1;
		
			const Triangle2 & K(Th2.elements[ii]);
			int iv[3];
			
			for(int jj=0; jj <3; jj++){
				iv[jj] = Numero_Som[ Th2.operator()(K[jj]) ];
			}
			
			for(int jj=0; jj<3; jj++){
				for(int kk=jj+1; kk<3; kk++){				
					if( iv[jj]==iv[kk] ) Border_ok=0;
				}
			}					
			if(Border_ok==1){
				ind_nbe_t[i_border] = ii;
				i_border=i_border+1;
			}
		}
		nbe_t = i_border;
		
		if( recollement_border == 1){
			int point_confondus_ok=1;
			cout << "debut recollement : nbe_t= "<< nbe_t << endl; 
			
			int np,dim=3;
			int *ind_np = new int [nbe_t];
			double **Cdg_be=new double *[nbe_t];
			for(int i=0; i<nbe_t; i++) Cdg_be[i] = new double[dim];
		
			for( int i_border=0; i_border< nbe_t; i_border++){
				//if(i_border>161) cout << "i_border=" << i_border << " " << nbe_t << endl;
				int & ii=ind_nbe_t[i_border];
				//if(i_border>161) cout <<"ii=" << ii << endl;
				const Triangle2 & K(Th2.elements[ii]);
				int iv[3];
				
				for(int jj=0; jj <3; jj++){
					iv[jj] = Th2.operator()(K[jj]);
				}	
				
								
				Cdg_be[i_border][0] = ( tab_XX[iv[0]] + tab_XX[iv[1]] + tab_XX[iv[2]] )/3.; 
				Cdg_be[i_border][1] = ( tab_YY[iv[0]] + tab_YY[iv[1]] + tab_YY[iv[2]] )/3.; 
				Cdg_be[i_border][2] = ( tab_ZZ[iv[0]] + tab_ZZ[iv[1]] + tab_ZZ[iv[2]] )/3.; 
			}
	
			cout << "points commun " << endl;
			PointCommun_hcode( dim, nbe_t, point_confondus_ok, Cdg_be, ind_np, np);
			cout << "points commun finis " <<endl;
			assert( np <= nbe_t );
		
			int *ind_nbe_t_tmp= new int [np];
		
			for( int i_border=0; i_border<np; i_border++){
				ind_nbe_t_tmp[ i_border ] = ind_nbe_t[ ind_np[i_border] ]; 
			}	
					
			for( int i_border=0; i_border< np; i_border++){
				ind_nbe_t[ i_border ] = ind_nbe_t_tmp[ i_border ]; 
			}	
		
			nbe_t = np;

			cout << "fin recollement : nbe_t= "<< nbe_t << endl; 
		}
}

void OrderVertexTransfo_hcode_nv( const int tab_nv, const double *tab_XX, const double *tab_YY, const double *tab_ZZ,  
	int *Numero_Som, int * ind_nv_t, int & nv_t){

	size_t i;
	size_t j[3];
	size_t k[3];	
	int NbCode = 10000;
	int *tcode; //= new int[NbCode];	

	int *posv = new int[tab_nv];	
	double epsilon=0.000001;
			
	// determination de boite englobante
	double bmin[3],bmax[3];
	
	bmin[0] = tab_XX[0];
	bmin[1] = tab_YY[0];
	bmin[2] = tab_ZZ[0];
	
	bmax[0] = bmin[0];
	bmax[1] = bmin[1];
	bmax[2] = bmin[2];
	
	cout << " determination bmin et bmax" << endl;
	
	for(int ii=1; ii<tab_nv; ii++){
 		bmin[0] = min(bmin[0],tab_XX[ii]);
 		bmin[1] = min(bmin[1],tab_YY[ii]);
 		bmin[2] = min(bmin[2],tab_ZZ[ii]);
 		
 		bmax[0] = max(bmax[0],tab_XX[ii]);
 		bmax[1] = max(bmax[1],tab_YY[ii]);
 		bmax[2] = max(bmax[2],tab_ZZ[ii]);
	}
	
	k[0] = int( (bmax[0]-bmin[0])/epsilon ); 
	k[1] = int( (bmax[1]-bmin[1])/epsilon ); 
	k[2] = int( (bmax[2]-bmin[2])/epsilon ); 
	
	
	//cout << "taille boite englobante =" << endl;
	//for(int ii=0; ii<3; ii++){
 	//	cout << "ii=" << ii << " " << bmin[ii] << " " << bmax[ii] <<endl;
	//}
		
	//for(int ii=0; ii<3; ii++){
 	//	cout << "k[" << ii << "]= " << k[ii]<<endl;
	//}
	
	NbCode = 10*(k[0]);		
	tcode = new int[NbCode];	
	/* initialisation des codes */
	for(int ii=0; ii< NbCode; ii++){
		 tcode[ii] = -1;
	}
	
	
	for(int ii=0; ii < tab_nv; ii++){
		// boucle dans l autre sens pour assurer l'ordre des elements pour la suite
		//cout << "vertex ii " << ii << "  max : " << tab_nv;  
		j[0] = int( (tab_XX[ii]-bmin[0])/epsilon  );
		j[1] = int( (tab_YY[ii]-bmin[1])/epsilon  );
		j[2] = int( (tab_ZZ[ii]-bmin[2])/epsilon  );  
		
		//for(int jj=0; jj<3; jj++){
 		//cout << "j[" << jj << "]= " << j[jj] <<endl;
 		//if( j[jj] > k[jj]) cout << "ii= " << ii <<" " << "tab_XX[ii]= " << tab_XX[ii] <<endl;
		//}
		
		assert( j[0] <=k[0] && j[0]>=0);
		assert( j[1] <=k[1] && j[1]>=0);
		assert( j[2] <=k[2] && j[2]>=0);
		i = (j[2]*(k[1]+1)+j[1]*(k[0]+1)+j[0])%NbCode;	
		//cout << i << endl;	
		i = i%NbCode;    
		//if(i<0) cout << k[0] << " " << k[1] << endl;
		//if(i<0) cout << j[0] << " " << j[1] << " "<< j[2] << endl;
		//cout << "ii= " << ii << " " << "i= " << i << endl;
			
		assert( i < abs(NbCode) );
		posv[ii] = tcode[i];
		tcode[i] = ii;
	}	
	
	cout << " boucle numero de Sommet " << endl;
	for(int ii=0; ii<tab_nv; ii++){
		Numero_Som[ii]=-1;
	}		
	
	
	cout << " determinations des points confondus et numerotation " << endl;
	nv_t=0;
	// maintenant determinations des points confondus et numerotation
	
	for(int	icode =0; icode < NbCode; icode++){
		int ii,jj;		
		double dist;		
		ii = tcode[icode];	
		if(ii==-1) continue;
		assert(ii!=-1);
		
		do{			
			
			if( Numero_Som[ii] != -1 ){
				ii=posv[ii];
				continue;
			}
			else{
				Numero_Som[ii] = nv_t;		
			}
			jj = posv[ii];			
		
			if( jj == -1){
				ind_nv_t[nv_t] = ii;	
				nv_t = nv_t+1;
		
				ii=posv[ii];
				continue;
			}

			do{
				if( Numero_Som[jj] != -1) {
					jj=posv[jj];	
					continue;
				}
				dist=pow(tab_XX[jj]-tab_XX[ii],2)+pow(tab_YY[jj]-tab_YY[ii],2)+pow(tab_ZZ[jj]-tab_ZZ[ii],2);	
				
				if( dist < epsilon ){
					// point semblable
					if( Numero_Som[ii] == 0 ){
						cout << "dist= " << dist << endl;
						cout << "composante du point et du semblable" << ii << " "<< jj << endl;
						cout << tab_XX[ii] << " "<<  tab_YY[ii] << " " << tab_ZZ[ii] << endl;
						cout << tab_XX[jj] << " "<<  tab_YY[jj] << " " << tab_ZZ[jj] << endl;
					}
											 
					Numero_Som[jj] = Numero_Som[ii];
				}	
				
				jj=posv[jj];				
			}while(jj!=-1);	
			
			
			ind_nv_t[nv_t] = ii;	
			nv_t = nv_t+1;
			
			ii=posv[ii];		
				
		}while( ii!=-1 ); 	
	}
	/*
	nv_t=0;
	for(int	icode =0; icode < NbCode; icode++){
		//int ii,jj;		
		double dist;		
		
		for(int ii=tcode[icode]; ii!=-1; ii=posv[ii]){
			if(Numero_Som[ii] != -1) continue; 
			Numero_Som[ii] = nv_t;				
			for(int jj=posv[ii]; jj!=-1; jj=posv[jj]){
				if(Numero_Som[jj] != -1) continue; 
				dist=pow(tab_XX[jj]-tab_XX[ii],2)+pow(tab_YY[jj]-tab_YY[ii],2)+pow(tab_ZZ[jj]-tab_ZZ[ii],2);	
			
				if( dist < epsilon ){
					// point semblable
					Numero_Som[jj] = Numero_Som[ii];
					//cout << "point semblable" << endl;
					//exit(-1); 
				}				
				
			}
			ind_nv_t[nv_t] = ii;	// Remarque on donne a nv_t le plus grand
			nv_t++; //nv_t = nvt+1;
		}
	}
	*/
	cout << "nv_t = " << nv_t << " / " << "nv_t(anc)" << tab_nv <<endl;   
	//exit(-1);
}

void PointCommun_hcode( const int dim, const int NbPoints, const int point_confondus_ok, double **Coord_Point,  
	int * ind_np, int & np){

	size_t i;
	size_t j[dim];
	size_t k[dim];	
	int NbCode = 10000;
	int *tcode; //= new int[NbCode];	
	int *posv = new int[NbPoints];
	int *Numero_Som =new int[NbPoints];
	
	double epsilon=0.0001;
			
	// determination de boite englobante
	double bmin[dim],bmax[dim];
	
	assert( dim > 1);
	
	for(int jj=0; jj<dim; jj++){ 
		bmin[jj] = Coord_Point[0][jj];	
		bmax[jj] = bmin[jj];
	}
	for(int ii=1; ii<NbPoints; ii++){
		for(int jj=0; jj<dim; jj++){ 
 			bmin[jj] = min(bmin[jj],Coord_Point[ii][jj]);
 			bmax[jj] = max(bmax[jj],Coord_Point[ii][jj]);	
		}
	}
	
	for(int jj=0; jj<dim; jj++){ 
		k[jj] = int( (bmax[jj]-bmin[jj])/epsilon );
	} 
	
	NbCode = 10*(k[0]);		
	tcode = new int[NbCode];	
	/* initialisation des codes */
	for(int ii=0; ii< NbCode; ii++){
		 tcode[ii] = -1;
	}
	
	for(int ii=0; ii < NbPoints; ii++){
		// boucle dans l autre sens pour assurer l'ordre des elements pour la suite
		
		for( int jj=0; jj<dim; jj++){
			j[jj] = int( (Coord_Point[ii][jj]-bmin[jj])/epsilon  );
		}
	
		assert( j[0] <=k[0] && j[0]>=0);
		assert( j[1] <=k[1] && j[1]>=0);
		assert( j[2] <=k[2] && j[2]>=0);
		
		i = j[0];
		for(int jj=1; jj<dim; jj++){
			i=i+j[jj]*(k[jj-1]+1);
		}	
		i = i%NbCode;    
		
		assert( i < NbCode );
		posv[ii] = tcode[i];
		tcode[i] = ii;
	}	

	for(int ii=0; ii < NbPoints; ii++){
		ind_np[ii] = -1;
		Numero_Som[ii] = -1;
	}
	
	/* Resolution probleme dans le cas où le maillage se colle */
	int minimum_np;
		
	np=0;
	/* maintenant determinations des points confondus et numerotation*/
	
	/*
	for(int	icode =0; icode < NbCode; icode++){
		int ii,jj;		
		double dist;	
			
		ii = tcode[icode];
		assert( ii < NbPoints);	
		if(ii==-1) continue;
			
		do{				
			if( Numero_Som[ii]!= -1 ){
				ii=posv[ii];				
				continue;
			}
			else{				
				Numero_Som[ii] = np;
				minimum_np=ii;
				//ind_np[np] = ii;						
			}
			
			jj = posv[ii];			
		
			if( jj == -1){
				ind_np[np] = ii;	
				np = np+1;
		
				ii=posv[ii];
				continue;
			}

			do{
			
				if( Numero_Som[jj] != -1) {
					jj=posv[jj];	
					continue;
				}
				
				dist = 0.;
				for( int kk=0; kk<dim; kk++){
					dist = dist +  pow(Coord_Point[jj][kk]-Coord_Point[ii][kk],2);
				}
			
				if( dist < epsilon ){
					// point semblable
					Numero_Som[jj] = Numero_Som[ii];
					minimum_np = min(minimum_np,jj);
					cout << "(minimum_np,jj)= " << minimum_np << " " << jj << " info=" << ii << endl;
 					// endroit pour savoir combien de points sont double 					
				}	
				
				jj=posv[jj];				
			}while(jj!=-1);	
			
			ind_np[np] = min(ii,minimum_np);	
			np = np+1;
			
			ii=posv[ii];		
			//cout << "fin: ii = " << ii << endl;	
		}while( ii!=-1 ); 
			
	}*/	
	
	switch( point_confondus_ok ){
	
	case 0:
		np=0;
		for(int	icode =0; icode < NbCode; icode++){
			//int ii,jj;		
			double dist;		
		
			for(int ii=tcode[icode]; ii!=-1; ii=posv[ii]){
			
				if( Numero_Som[ii]!= -1 ) continue;										
				//minimum_np=ii;
				Numero_Som[ii] = np;
			
				for(int jj=posv[ii]; jj!=-1; jj=posv[jj]){
				
					if(Numero_Som[jj] != -1) continue; 				
					dist = 0.;
				
					for( int kk=0; kk<dim; kk++){
						dist = dist +  pow(Coord_Point[jj][kk]-Coord_Point[ii][kk],2);
					}
							
					if( dist < epsilon ){
						// point semblable
						Numero_Som[jj] = Numero_Som[ii];
						//minimum_np = min( jj, minimum_np);
					}				
				}
				ind_np[np] = ii; //min(ii,minimum_np);	// Remarque on donne a np le plus petit element
				np++; //nv_t = nvt+1;
			}
		}
		break;
		
	case 1:
		int point_multiple;
		np=0;
		for(int	icode =0; icode < NbCode; icode++){
			//int ii,jj;		
			double dist;		
		
			for(int ii=tcode[icode]; ii!=-1; ii=posv[ii]){
			
				if( Numero_Som[ii]!= -1 ) continue;										
				//minimum_np=ii;
				Numero_Som[ii] = np;
				point_multiple = 0; 
				
				
				for(int jj=posv[ii]; jj!=-1; jj=posv[jj]){
				
					if(Numero_Som[jj] != -1) continue; 				
					dist = 0.;
				
					for( int kk=0; kk<dim; kk++){
						dist = dist +  pow(Coord_Point[jj][kk]-Coord_Point[ii][kk],2);
					}
							
					if( dist < epsilon ){
						// point semblable
						Numero_Som[jj] = Numero_Som[ii];
						point_multiple = 1; 
						//minimum_np = min( jj, minimum_np);
					}				
				}
				if(point_multiple ==0){
					ind_np[np] = ii; //min(ii,minimum_np);	// Remarque on donne a np le plus petit element
					np++; //nv_t = nvt+1;
				}
			}
		}
		break;
	default:
	cout << " point_confondus_ok dans fonction PointCommun_hcode vaut 1 ou 0." << endl;
	exit(-1);
	
	}
}

