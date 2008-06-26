// $Id$

#include  <iostream>
#include  <cfloat>
#include <cmath>
#include <complex>
using namespace std;
#include "error.hpp"
#include "AFunction.hpp"
using namespace std;  
#include "rgraph.hpp"
#include "RNM.hpp"
#include "fem.hpp"


#include "FESpacen.hpp" 
#include "FESpace.hpp" 

#include "MatriceCreuse_tpl.hpp"
#include "MeshPoint.hpp"
#include "Operator.hpp" 
#include "lex.hpp"

#include "lgfem.hpp"
#include "lgmesh3.hpp"
#include "lgsolver.hpp"
#include "problem.hpp"
#include "LayerMesh.hpp"
#include "TransfoMesh_v2.hpp"
#include "tetgen.h"
//#include "GQuadTree.hpp"

#include <set>
#include <vector>
#include <fstream>


using namespace  Fem2D;

// subroutine use for tetegen call

void mesh3_tetgenio_out(const tetgenio &out, const int & label_tet, Mesh3 & Th3);

Mesh3 * RemplissageSurf3D_tetgen(const Mesh3 & Th3, const int & label_tet);
	
Mesh3 * Transfo_Mesh2_tetgen(const Mesh & Th2, const double *tab_XX, const double *tab_YY, const double *tab_ZZ, 
		int &border_only, int &recollement_border, int &point_confondus_ok, 
		const int &label_tet,const map<int, int> &maptri );



class Build2D3D_Op : public E_F0mps 
{
public:
  Expression eTh;
  Expression xx,yy,zz;
  static const int n_name_param =3; // 
  static basicAC_F0::name_and_type name_param[] ;
  Expression nargs[n_name_param];
  KN_<long>  arg(int i,Stack stack,KN_<long> a ) const
  { return nargs[i] ? GetAny<KN_<long> >( (*nargs[i])(stack) ): a;}
  double  arg(int i,Stack stack,double a ) const{ return nargs[i] ? GetAny< double >( (*nargs[i])(stack) ): a;}
  
public:
  Build2D3D_Op(const basicAC_F0 &  args,Expression tth) 
    : eTh(tth),xx(0),yy(0),zz(0)
  {
    cout << "construction par BuilLayeMesh_Op" << endl;
    args.SetNameParam(n_name_param,name_param,nargs);
    const E_Array * a1=0 ;
    if(nargs[0])  a1  = dynamic_cast<const E_Array *>(nargs[0]);
    int err =0;
   
    if(a1) {
      if(a1->size() !=3) 
      CompileError("Build2D3D (Th,transfo=[X,Y,Z],) ");
      xx=to<double>( (*a1)[0]);
      yy=to<double>( (*a1)[1]);
      zz=to<double>( (*a1)[2]);
    }    
  } 
  
  AnyType operator()(Stack stack)  const ;
};


basicAC_F0::name_and_type Build2D3D_Op::name_param[]= {
  {  "transfo", &typeid(E_Array)},
  {  "reftet", &typeid(KN_<long>)}, 
  {  "refface", &typeid(KN_<long> )}
};

class  Build2D3D : public OneOperator { public:  
    Build2D3D() : OneOperator(atype<pmesh3>(),atype<pmesh>()) {}
  
  E_F0 * code(const basicAC_F0 & args) const 
  { 
	return  new Build2D3D_Op( args,t[0]->CastTo(args[0]) ); 
  }
};

AnyType Build2D3D_Op::operator()(Stack stack)  const 
{
  MeshPoint *mp(MeshPointStack(stack)) , mps=*mp;
  Mesh * pTh= GetAny<Mesh *>((*eTh)(stack));
  ffassert( pTh );
  Mesh &Th=*pTh;
  Mesh *m= pTh;   // question a quoi sert *m ??
  int nbv=Th.nv; // nombre de sommet 
  int nbt=Th.nt; // nombre de triangles
  int neb=Th.neb; // nombre d'aretes fontiere
  cout << " Vertex Triangle Border " << nbv<< "  "<< nbt << " nbe "<< neb << endl; 
 
  KN<long> zzempty;
  KN<long> nrt (arg(1,stack,zzempty));  
  KN<long> nrf (arg(2,stack,zzempty));
  
  int label_tet;
  
  if( nrt.N() < 0 || nrt.N() > 1){
	  cout << "tetgen allow one label for tetrahedra " << endl;
	  assert(  nrt.N() < 0 || nrt.N() > 1 );
  } 
   
  if( nrt.N() == 1 ){
	  label_tet=nrt[0];
  }
  else{
	  label_tet=0;
  }  
  
  ffassert( nrf.N() %2 ==0);
  
  map<int,int> mapf;
  for(int i=0;i<nrf.N();i+=2)
  {
	  if(nrf[i] != nrf[i+1]){
		mapf[nrf[i]]=nrf[i+1];
		}
  }
  
  map<int, int> mapfme;  
  
  Transfo_Mesh2_map_face( Th, mapfme );  
  
  // Map utilisateur
  map< int, int > :: iterator imap;
  for( int ii=0; ii < nrf.N(); ii+=2){
	imap = mapfme.find(nrf[ii]);
	if( imap != mapfme.end()){
		imap -> second = nrf[ii+1];
	}
  }  
  
	KN<double> txx(nbv), tyy(nbv), tzz(nbv);
	KN<int> takemesh(nbv);
	MeshPoint *mp3(MeshPointStack(stack)); 
	
	takemesh=0;  
	Mesh &rTh = Th;
	for (int it=0; it<nbt; ++it){
		  for( int iv=0; iv<3; ++iv){
				int i=Th(it,iv);
				if(takemesh[i]==0){
					mp3->setP(&Th,it,iv);
					if(xx){ 
						txx[i]=GetAny<double>((*xx)(stack));
					}
					if(yy){ 
						tyy[i]=GetAny<double>((*yy)(stack));
					}
					if(zz){ 
						tzz[i]=GetAny<double>((*zz)(stack));
					}
					takemesh[i] = takemesh[i]+1;
				}
	 		}
		}
		int border_only = 0;
		int recollement_border=1, point_confondus_ok=1;
		Mesh3 *Th3=Transfo_Mesh2_tetgen( Th, txx, tyy, tzz, border_only, 
			recollement_border, point_confondus_ok, label_tet, mapfme);  
		
		Th3->BuildBound();
		Th3->BuildAdj();
		Th3->Buildbnormalv();  
		Th3->BuildjElementConteningVertex();
		Th3->BuildGTree();
		Th3->decrement();    
		
		*mp=mps;
		return Th3;
}

// Fonction pour tetgen

void mesh3_tetgenio_out(const tetgenio &out, const int & label_tet, Mesh3 & Th3)
{ 
  int i;

// All indices start from 1.
  if(out.firstnumber != 1){
    cout << " probleme ???" << endl;
    exit(1);
  }   
  
  if(out.numberoffacets !=0){
    cout << "tetgen: faces non triangulaire" << endl;
    exit(1);
  }
  
  if(out.numberofcorners !=4){
    cout << "tetgen: element subparametric of order 2" <<endl;
    exit(1);
  }
  
  cout << "Th3 :: Vertex Element Border :: " << out.numberofpoints << " " <<out.numberoftetrahedra  << " " << out.numberoftrifaces << endl;
  Th3.set(out.numberofpoints, out.numberoftetrahedra, out.numberoftrifaces);
   
  i=0;
  for(int nnv=0; nnv < Th3.nv; nnv++){
    Th3.vertices[nnv].x=out.pointlist[i];
    Th3.vertices[nnv].y=out.pointlist[i+1];
    Th3.vertices[nnv].z=out.pointlist[i+2];       
    Th3.vertices[nnv].lab=out.pointmarkerlist[nnv];
    i=i+3;    
  }
    
  i=0;
  for(int nnt=0; nnt < Th3.nt; nnt++){
    int iv[4],lab;
    iv[0] = out.tetrahedronlist[i]-1;
    iv[1] = out.tetrahedronlist[i+1]-1;
    iv[2] = out.tetrahedronlist[i+2]-1;
    iv[3] = out.tetrahedronlist[i+3]-1;
    lab   = label_tet;
    //lab = out.tetrahedronmarkerlist[nnt];
    Th3.elements[nnt].set( Th3.vertices, iv, lab);
    i=i+4;
  }

  for(int ibe=0; ibe < Th3.nbe; ibe++){
    int iv[3];
    iv[0] = out.trifacelist[3*ibe]-1;
    iv[1] = out.trifacelist[3*ibe+1]-1;
    iv[2] = out.trifacelist[3*ibe+2]-1;  
    Th3.be(ibe).set( Th3.vertices, iv, out.trifacemarkerlist[ibe]);
  }
}  

Mesh3 * RemplissageSurf3D_tetgen(const Mesh3 & Th3, const int & label_tet){
	
	Mesh3 *T_Th3= new Mesh3;
	
	assert(Th3.nt == 0 );
	int nv_t = Th3.nv;
	int nt_t = Th3.nt;
	int nbe_t = Th3.nbe;
	
	cout << "3D RemplissageSurf3D:: Vertex  triangle2  border " 
	<< nv_t << " "<< nt_t << " " << nbe_t<< endl;
	// Creation des tableau de tetgen
	
	tetgenio in,out;
	//tetgenio::facet *f;
	//tetgenio::polygon *p;
	
	cout << " tetgenio: vertex " << endl;
	int itet,jtet;
// All indices start from 1.
	in.firstnumber = 1;
	in.numberofpoints = nv_t;
	in.pointlist = new REAL[in.numberofpoints*3];
	in.pointmarkerlist = new int[in.numberofpoints];
	itet=0;
	jtet=0;
	for(int nnv=0; nnv < nv_t; nnv++)
	{
		in.pointlist[itet]   = Th3.vertices[nnv].x;
		in.pointlist[itet+1] = Th3.vertices[nnv].y;
		in.pointlist[itet+2] = Th3.vertices[nnv].z;       
		in.pointmarkerlist[nnv] =  Th3.vertices[nnv].lab;   
		itet=itet+3;
	}
	assert(itet==in.numberofpoints*3);
	
	cout << " tetgenio: facet " << endl;
	// Version avec des facettes
	in.numberoffacets = nbe_t;
	in.facetlist = new tetgenio::facet[in.numberoffacets];
	in.facetmarkerlist = new int[in.numberoffacets];
	  
	for(int ibe=0; ibe < nbe_t; ibe++){
		tetgenio::facet *f;
		tetgenio::polygon *p;
		f = &in.facetlist[ibe];
		f->numberofpolygons = 1;
		f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
		f->numberofholes = 0;
		f->holelist = NULL;
		
		p = &f->polygonlist[0];
		p->numberofvertices = 3;
		p->vertexlist = new int[3];
		
		// creation of elements
		const Triangle3 & K(Th3.be(ibe)); // const Triangle2 & K(Th2.elements[ii]); // Version Mesh2  
		p->vertexlist[0] = Th3.operator()(K[0])+1;
		p->vertexlist[1] = Th3.operator()(K[1])+1;
		p->vertexlist[2] = Th3.operator()(K[2])+1;
		
		for( int kkk=0; kkk<3; kkk++){ 
			assert( p->vertexlist[kkk]<= in.numberofpoints && p->vertexlist[kkk]>0);
		}
		
		in.facetmarkerlist[ibe] = K.lab; 
		
	}  
	cout << "debut de tetrahedralize( , &in, &out);" << endl;
	
	tetrahedralize("pqCVV", &in, &out);
	 
	cout << "fin de tetrahedralize( , &in, &out);" << endl;
	mesh3_tetgenio_out( out, label_tet, *T_Th3);
	
	cout <<" Finish Mesh3 tetgen :: Vertex, Element, Border" << T_Th3->nv << " "<< T_Th3->nt << " " << T_Th3->nbe << endl;
	
	return T_Th3;
}

Mesh3 * Transfo_Mesh2_tetgen(const Mesh & Th2, const double *tab_XX, const double *tab_YY, const double *tab_ZZ, 
	int &border_only, int &recollement_border, int &point_confondus_ok, 
	const int &label_tet, const map<int, int> &maptri ){
	Mesh3 *T_Th3= new Mesh3;
	int nv_t,nt_t,nbe_t;
	int* Numero_Som;
	
	int* ind_nv_t;
	int* ind_nt_t;
	int* ind_nbe_t;
	
	int* label_nbe_t;
	
	//int i_som;
	Numero_Som = new int[Th2.nv];
	ind_nv_t   = new int[Th2.nv];
	ind_nbe_t  = new int[Th2.nt];
	     
	label_nbe_t = new int[Th2.nt];
	
    cout << "2D: Mesh::Vertex  triangle2  border " << Th2.nv << " "<<Th2.nt<< " " << Th2.neb<< endl;
	
	for(int ii=0; ii<Th2.nv; ii++){ 
		Numero_Som[ii]=ii;
	}
	cout <<" debut: SamePointElement " <<endl;
	
	SamePointElement_Mesh2( tab_XX, tab_YY, tab_ZZ, Th2, recollement_border, point_confondus_ok, 
		Numero_Som, ind_nv_t, ind_nt_t, ind_nbe_t, label_nbe_t, nv_t, nt_t, nbe_t);
	
	cout <<" fin: SamePointElement " <<endl;
	
	cout << "2D transfo: Mesh::Vertex  triangle2  border " << nv_t << " "<< nt_t << " " << nbe_t<< endl;
	// Creation des tableau de tetgen
	
	
	tetgenio in,out;
	//tetgenio::facet *f;
	//tetgenio::polygon *p;
	
	cout << " tetgenio: vertex " << endl;
	int itet,jtet;
// All indices start from 1.
	in.firstnumber = 1;
	in.numberofpoints = nv_t;
	in.pointlist = new REAL[in.numberofpoints*3];
	in.pointmarkerlist = new int[in.numberofpoints];
	itet=0;
	jtet=0;
	for(int nnv=0; nnv < nv_t; nnv++)
	{
		int & ii = ind_nv_t[nnv];
		//cout << "nnv ,  ii  =" << nnv << "  " << ii << endl;
		//cout << "tab_XX[ii], tab_YY[ii], tab_ZZ[ii]=" <<  tab_XX[ii] << " "<< tab_YY[ii] << " "<< tab_ZZ[ii] << endl;
		assert( Numero_Som[ii] == nnv );
		const Mesh::Vertex & K = Th2.vertices[ii];//const Mesh::Vertex & K(Th2.vertices[ii]); //Version Mesh2   
		in.pointlist[itet]   = tab_XX[ii];
		in.pointlist[itet+1] = tab_YY[ii];
		in.pointlist[itet+2] = tab_ZZ[ii];       
		in.pointmarkerlist[nnv] =  K.lab;   
		itet=itet+3;
	}
	assert(itet==in.numberofpoints*3);
	
	cout << " tetgenio: facet " << endl;
	// Version avec des facettes
	in.numberoffacets = nbe_t;
	in.facetlist = new tetgenio::facet[in.numberoffacets];
	in.facetmarkerlist = new int[in.numberoffacets];
	  
	for(int ibe=0; ibe < nbe_t; ibe++){
		tetgenio::facet *f;
		tetgenio::polygon *p;
		f = &in.facetlist[ibe];
		f->numberofpolygons = 1;
		f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
		f->numberofholes = 0;
		f->holelist = NULL;
		
		p = &f->polygonlist[0];
		p->numberofvertices = 3;
		p->vertexlist = new int[3];
		
		int & ii=ind_nbe_t[ibe];
		// creation of elements
		const Mesh::Triangle & K(Th2.t(ii)); // const Triangle2 & K(Th2.elements[ii]); // Version Mesh2  
		p->vertexlist[0] = Numero_Som[ Th2.operator()(K[0]) ]+1;
		p->vertexlist[1] = Numero_Som[ Th2.operator()(K[1]) ]+1;
		p->vertexlist[2] = Numero_Som[ Th2.operator()(K[2]) ]+1;
		
		for( int kkk=0; kkk<3; kkk++){ 
			assert( p->vertexlist[kkk]<= in.numberofpoints && p->vertexlist[kkk]> 0);
		}
		map< int, int>:: const_iterator imap;
		imap = maptri.find(K.lab); // imap= maptri.find( label_nbe_t[ibe] );
		assert( imap != maptri.end());
		in.facetmarkerlist[ibe] = imap->second; // K.lab; // before 
		
	}  
	cout << "debut de tetrahedralize( , &in, &out);" << endl;
	
	tetrahedralize("pqCV", &in, &out);
	 
	cout << "fin de tetrahedralize( , &in, &out);" << endl;
	mesh3_tetgenio_out( out, label_tet, *T_Th3);
	
	cout <<" Finish Mesh3 :: Vertex, Element, Border" << T_Th3->nv << " "<< T_Th3->nt << " " << T_Th3->nbe << endl;
	
	return T_Th3;
}


class Remplissage_Op : public E_F0mps 
{
public:
  Expression eTh;
  static const int n_name_param =3; // 
  static basicAC_F0::name_and_type name_param[] ;
  Expression nargs[n_name_param];
  KN_<long>  arg(int i,Stack stack,KN_<long> a ) const
  { return nargs[i] ? GetAny<KN_<long> >( (*nargs[i])(stack) ): a;}
  double  arg(int i,Stack stack,double a ) const{ return nargs[i] ? GetAny< double >( (*nargs[i])(stack) ): a;}
  
public:
  Remplissage_Op(const basicAC_F0 &  args,Expression tth) 
    : eTh(tth)
  {
    cout << "Remplissage du bord" << endl;
    args.SetNameParam(n_name_param,name_param,nargs);
  } 
  
  AnyType operator()(Stack stack)  const ;
};


basicAC_F0::name_and_type Remplissage_Op::name_param[]= {
  {  "methode",&typeid(KN_<long>)},
  {  "reftet", &typeid(KN_<long>)}, 
  {  "refface", &typeid(KN_<long> )}
};

class  Remplissage : public OneOperator { public:  
    Remplissage() : OneOperator(atype<pmesh3>(),atype<pmesh3>()) {}
  
  E_F0 * code(const basicAC_F0 & args) const 
  { 
	return  new Remplissage_Op( args,t[0]->CastTo(args[0]) ); 
  }
};

AnyType Remplissage_Op::operator()(Stack stack)  const 
{
  MeshPoint *mp(MeshPointStack(stack)) , mps=*mp;
  Mesh3 * pTh= GetAny<Mesh3 *>((*eTh)(stack));
  ffassert( pTh );
  Mesh3 &Th=*pTh;
  Mesh3 *m= pTh;   // question a quoi sert *m ??
  int nbv=Th.nv; // nombre de sommet 
  int nbt=Th.nt; // nombre de triangles
  int nbe=Th.nbe; // nombre d'aretes fontiere
  cout << " Vertex Triangle Border " << nbv<< "  "<< nbt << " nbe "<< nbe << endl; 
 
  KN<long> zzempty;
  KN<long> num_methode(arg(0,stack,zzempty));
  KN<long> nrt (arg(1,stack,zzempty));  
  KN<long> nrf (arg(2,stack,zzempty));
  
  int label_tet;
  
  if( nrt.N() < 0 || nrt.N() > 1){
	  cout << "tetgen allow one label for tetrahedra " << endl;
	  assert(  nrt.N() < 0 || nrt.N() > 1 );
  } 
   
  if( nrt.N() == 1 ){
	  label_tet=nrt[0];
  }
  else{
	  label_tet=0;
  }  
  
  ffassert( nrf.N() %2 ==0);
  
  
  map<int,int> mapf;
  for(int i=0;i<nrf.N();i+=2)
  {
	  if(nrf[i] != nrf[i+1]){
		mapf[nrf[i]]=nrf[i+1];
		}
  }
	Mesh3 *Th3 =new Mesh3;
	if(num_methode[0]==1) 
		Th3 = RemplissageSurf3D_tetgen( Th, label_tet);
	
	// changement de label 
	if( nrf.N() > 0){
		for(int ii=0; ii< Th3->nbe; ii++){
			const Triangle3 & K(Th3->be(ii));
			int lab;
			int iv[3];
			
			iv[0] = Th3->operator()(K[0]);
			iv[1] = Th3->operator()(K[1]);
			iv[2] = Th3->operator()(K[2]);
			
			map< int, int> :: const_iterator imap;
			imap = mapf.find(K.lab);
			if( imap != mapf.end() ){
				lab = imap -> second;
			}
			else{
				lab = K.lab;
			}
			Th3->be(ii).set(Th3->vertices,iv,lab);
		}
	}
	
	Th3->BuildBound();
	Th3->BuildAdj();
	Th3->Buildbnormalv();  
	Th3->BuildjElementConteningVertex();
	Th3->BuildGTree();
	Th3->decrement();    
	
	*mp=mps;
	return Th3;
}

class Init { public:
  Init();
};

static Init init;  //  une variable globale qui serat construite  au chargement dynamique 

Init::Init(){  // le constructeur qui ajoute la fonction "splitmesh3"  a freefem++ 
  
  if (verbosity)
    cout << " load: tetgen  " << endl;
   
  Global.Add("tetgtransfo","(",new Build2D3D);
  Global.Add("tetg","(",new Remplissage);
}
