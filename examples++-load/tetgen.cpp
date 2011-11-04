// ORIG-DATE:     Juin 2008
// -*- Mode : c++ -*-
//
// SUMMARY  : liaison medit freefem++ : popen  
// USAGE    : LGPL      
// ORG      : LJLL Universite Pierre et Marie Curie, Paris,  FRANCE 
// AUTHOR   : Jacques Morice
// E-MAIL   : jacques.morice@ann.jussieu.fr
//
//   for automatic  compilation with ff-c++
//ff-c++-LIBRARY-dep:   tetgen 
//ff-c++-cpp-dep: 
//  

/* 
 This file is part of Freefem++
 
 Freefem++ is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation; either version 2.1 of the License, or
 (at your option) any later version.
 
 Freefem++  is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.
 
 You should have received a copy of the GNU Lesser General Public License
 along with Freefem++; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

 Thank to the ARN ()  FF2A3 grant
 ref:ANR-07-CIS7-002-01 
 */
// $Id$

//  FH   July 2009
//   comment all
//       Th3_t->BuildBound();
//       Th3_t->BuildAdj();
//       Th3_t->Buildbnormalv();  
//       Th3_t->BuildjElementConteningVertex();
//   is now in the constructor of Mesh3 to be consistante. 
//   


#include "ff++.hpp" 
#include "msh3.hpp"
#define TETLIBRARY
#include "tetgen.h"

using namespace  Fem2D;

/*
// function to return inside point in a volume mesh 
// A rajouter par la suite //
void insidepoint( const Mesh3 &Th3

*/
// subroutine use for tetegen call

void mesh3_tetgenio_out(const tetgenio &out, Mesh3 & Th3);

void mesh3_tetgenio_out(const tetgenio &out, const int & label_tet, Mesh3 & Th3);

void mesh3_tetgenio_out(const tetgenio &out, const int & label_tet, const int & label_face, Mesh3 & Th3);

Mesh3 * mesh3_tetgenio_out(const tetgenio &out);

Mesh3 *mesh3_tetgenio_out(const tetgenio &out, const int & label_tet);

Mesh3 * mesh3_tetgenio_out(const tetgenio &out, const int & label_tet, const int & label_face);


Mesh3 * Convexhull_3Dpoints( char* switch_tetgen, const int &nv_t, const double *Xcoord, const double *Ycoord, const double *Zcoord, const int &label_tet);

Mesh3 * RemplissageSurf3D_tetgen( char* switch_tetgen, const Mesh3 & Th3, const int & label_tet);

Mesh3 * RemplissageSurf3D_tetgen_new( char* switch_tetgen, const Mesh3 & Th3, const int & label_tet,
				  const int &nbhole, const double *tabhole, 
				  const int & nbregion, const double *tabregion, 
				  const int &nbfacecl, const double *tabfacecl);

Mesh3 * Transfo_Mesh2_tetgen( const double &precis_mesh, char* switch_tetgen, const Mesh & Th2, const double *tab_XX, const double *tab_YY, const double *tab_ZZ, 
		int &border_only, int &recollement_border, int &point_confondus_ok, 
		const int &label_tet,const map<int, int> &maptri );

Mesh3 *Transfo_Mesh2_tetgen_new(const double &precis_mesh, char *switch_tetgen,const Mesh & Th2, const double *tab_XX, const double *tab_YY, const double *tab_ZZ, 
				 int &border_only, int &recollement_border, int &point_confondus_ok, 
				 const int &label_tet, const map<int, int> &maptri, 
				 const int &nbhole, const double *tabhole, 
				 const int & nbregion, const double *tabregion, 
				 const int &nbfacecl, const double *tabfacecl);

Mesh3 *  ReconstructionRefine_tetgen(char *switch_tetgen,const Mesh3 & Th3, 
				     const int &nbhole, const double *tabhole, 
				     const int & nbregion, const double *tabregion, 
		                     const int &nbfacecl, const double *tabfacecl, const double *tsizevol);


class Build2D3D_Op : public E_F0mps 
{
public:
  Expression eTh;
  Expression xx,yy,zz;
  static const int n_name_param =13+2; // 
  static basicAC_F0::name_and_type name_param[] ;
  Expression nargs[n_name_param];
  KN_<long>  arg(int i,Stack stack,KN_<long> a ) const
  { return nargs[i] ? GetAny<KN_<long> >( (*nargs[i])(stack) ): a;}
  KN_<double>  arg(int i,Stack stack,KN_<double> a ) const
  { return nargs[i] ? GetAny<KN_<double> >( (*nargs[i])(stack) ): a;}
  double  arg(int i,Stack stack,double a) const{ return nargs[i] ? GetAny< double >( (*nargs[i])(stack) ): a;}
  long   arg(int i,Stack stack, long  a) const{ return nargs[i] ? GetAny< long >( (*nargs[i])(stack) ): a;}
  string*  arg(int i,Stack stack, string* a) const{ return nargs[i] ? GetAny< string* >( (*nargs[i])(stack) ): a;}
public:
  Build2D3D_Op(const basicAC_F0 &  args,Expression tth) 
    : eTh(tth),xx(0),yy(0),zz(0)
  {
    if(verbosity) cout << "construction par BuilLayeMesh_Op" << endl;
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
    if( nargs[2] && nargs[13] ) 
	CompileError("uncompatible movemesh3 (Th, region= , reftet=  ");
    if( nargs[3] && nargs[14] ) 
	CompileError("uncompatible movemesh3 (Th, label= , refface=  ");
    
  } 
  
  AnyType operator()(Stack stack)  const ;
};


basicAC_F0::name_and_type Build2D3D_Op::name_param[]= {
  {  "transfo", &typeid(E_Array)},//0
  {  "switch", &typeid(string*)},
  {  "reftet", &typeid(long)}, //2
  {  "refface", &typeid(KN_<long>)},//3
  {  "facemerge", &typeid(long)},
  {  "ptmerge", &typeid(double)},
  // nouvelle variable
  {  "nbofholes", &typeid(long)},//6
  {  "holelist", &typeid(KN_<double>)},
  {  "nbofregions", &typeid(long)},
  {  "regionlist", &typeid(KN_<double>)},
  {  "nboffacetcl", &typeid(long)},
  {  "facetcl", &typeid(KN_<double>)},//11
  // mesure mesh
  {  "mesuremesh", &typeid(long)},
  {  "region", &typeid(long)}, //13
  {  "label", &typeid(KN_<long>)}//14
    
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
  if(verbosity) cout << " Vertex Triangle Border " << nbv<< "  "<< nbt << " " << neb << endl; 

  if(verbosity >1) cout <<" ======================= " << endl;
  if(verbosity >1) cout <<" == Build2D_3D_Op==" << endl;
 
  KN<long> zzempty;
  string  stringempty = string("pqaAAYCQ");
  string* switch_tet= (arg(1,stack,&stringempty));  
  int label_tet(arg(2,stack,arg(13,stack,0L)));  
  KN<long> nrf (arg(3,stack,arg(14,stack,zzempty)));
  int point_confondus_ok(arg(4,stack,0L));
  double precis_mesh(arg(5,stack,-1.));

  // new parameters
  KN<double> zdzempty;
  int nbhole (arg(6,stack,0L));
  KN<double> tabhole (arg(7,stack,zdzempty));
  int nbregion (arg(8,stack,0L));
  KN<double> tabregion (arg(9,stack,zdzempty));
  int nbfacecl (arg(10,stack,0L));
  KN<double> tabfacecl (arg(11,stack,zdzempty));
  if(nbhole && nbhole*3 != tabhole.N())
    {ExecError(" nbhole and holes are incompatibale ");}
  if(!nbhole) nbhole=tabhole.N()/3; // modif FH dec 2010...
  // mesuremesh parameters
  int mesureM(arg(12,stack,1L));
  int surface_orientation=1; 
  if( mesureM <0 ){
    surface_orientation=-1;
  }
    
    
    if(nbregion==0) nbregion=tabregion.N()/5;
    if(nbhole==0) nbhole=tabhole.N()/3;
    if(nbfacecl==0) nbfacecl=tabfacecl.N()/2;
    
  // assertion au niveau de la taille
  ffassert( tabhole.N()   == 3*nbhole);
  ffassert( tabregion.N() == 5*nbregion);
  ffassert( tabfacecl.N() == 2*nbfacecl);

  //====================================
  //  How to change string* into char* 
  //====================================
  cout << "string" << switch_tet << endl; 

  size_t size_switch_tet = switch_tet->size()+1;
  char* switch_tetgen =new char[size_switch_tet];
  strncpy(switch_tetgen, switch_tet->c_str(), size_switch_tet); 

  cout << "switch_tetgen=" << switch_tetgen << endl;
  //exit(1);
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
  
  //KN<double> txx(nbv), tyy(nbv), tzz(nbv);
  //KN<int> takemesh(nbv);
  double *txx=new double[nbv];
  double *tyy=new double[nbv];
  double *tzz=new double[nbv];
  int *takemesh=new int[nbv];

  MeshPoint *mp3(MeshPointStack(stack)); 
  
  for(int ii=0; ii<nbv; ii++) 
      takemesh[ii]=0;  

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

  delete [] takemesh;
  int border_only = 0;
  int recollement_border=1;
  /*
    Mesh3 *Th3=Transfo_Mesh2_tetgen( precis_mesh, switch_tetgen, Th, txx, tyy, tzz, border_only, 
				   recollement_border, point_confondus_ok, label_tet, mapfme);  
  */

  Mesh3 *Th3_tmp = MoveMesh2_func( precis_mesh, Th, txx, tyy, tzz, border_only, recollement_border, point_confondus_ok);

  /* delete array */
  delete [] txx;
  delete [] tyy;
  delete [] tzz;

  /* check orientation of the mesh and flip if necessary*/ 
  Th3_tmp->flipSurfaceMesh3(surface_orientation);
  
  int addcheckorientation=0;
  if( addcheckorientation==1  ){
    cout << "check :: orientation des surfaces" << endl;
    Th3_tmp->BuildBoundaryElementAdj();
    cout << "fin check :: orientation des surfaces" << endl;
  }
  /* set label of surface Th3_tmp */
  for(int ii=0; ii< Th3_tmp->nbe; ii++)
    {
      const Triangle3 & K(Th3_tmp->be(ii)); 
      int iv[3];
      int lab;
     
      iv[0] = Th3_tmp->operator()(K[0]);
      iv[1] = Th3_tmp->operator()(K[1]);
      iv[2] = Th3_tmp->operator()(K[2]);
		
      map< int, int>:: const_iterator imap;
      imap = mapfme.find(K.lab); 
		
      if(imap!= mapfme.end()){
	lab=imap->second;	
      } 
      else{
	lab=K.lab;
      }
		
      Th3_tmp->be(ii).set( Th3_tmp->vertices, iv, lab ) ;
    }
  /* mesh domains with tetgen */ 
  Mesh3 *Th3 = RemplissageSurf3D_tetgen_new( switch_tetgen, *Th3_tmp, label_tet,
					     nbhole, tabhole, nbregion, tabregion, 
					     nbfacecl, tabfacecl);


  /*
  Mesh3 *Th3=Transfo_Mesh2_tetgen_new( precis_mesh, switch_tetgen, Th, txx, tyy, tzz, border_only, 
				   recollement_border, point_confondus_ok, label_tet, mapfme, 
				   nbhole, tabhole, nbregion, tabregion, nbfacecl,tabfacecl);

  */

  delete Th3_tmp;
  
//Th3->BuildBound();
//  Th3->BuildAdj();
//  Th3->Buildbnormalv();  
//  Th3->BuildjElementConteningVertex();
  Th3->BuildGTree();
  //Th3->decrement();    
  Add2StackOfPtr2FreeRC(stack,Th3);

  delete [] switch_tetgen;
		
  *mp=mps;
  cout << "FreeFem++: End check mesh given by tetgen" << endl;  
  return Th3;
}

// Fonction pour tetgen
// new parameter
void mesh3_tetgenio_out(const tetgenio &out, Mesh3 & Th3)
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

  // new parameter
  if( out.numberoftetrahedronattributes !=  1){
    cout << "out.numberoftetrahedronattributes" << out.numberoftetrahedronattributes << endl;
  }
   
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
    //lab   = label_tet;
    for(int jj=0; jj<4; jj++){
      assert( iv[jj] >=0 && iv[jj] <Th3.nv );
    }
    //cout << "nnt= " <<  nnt << " " << lab  << " " << out.tetrahedronattributelist[nnt] << endl;
    lab =  out.tetrahedronattributelist[nnt];
    //cout << "nnt= " <<  lab  << " " << out.tetrahedronattributelist[nnt] << endl;
      
    Th3.elements[nnt].set( Th3.vertices, iv, lab);
    i=i+4;
  }
 
  for(int ibe=0; ibe < Th3.nbe; ibe++){
    int iv[3];
    iv[0] = out.trifacelist[3*ibe]-1;
    iv[1] = out.trifacelist[3*ibe+1]-1;
    iv[2] = out.trifacelist[3*ibe+2]-1; 
    for(int jj=0; jj<3; jj++){
      if(iv[jj]>= Th3.nv || iv[jj]< 0 ) cout << "iv[jj]=" << iv[jj] << " triangle" << ibe << endl;
      assert( iv[jj] >=0 && iv[jj] <Th3.nv );
    }
    Th3.be(ibe).set( Th3.vertices, iv, out.trifacemarkerlist[ibe]);
  }

  /*
  if( out.numberoftetrahedronattributes != 1 ){
    cout << "out.numberoftetrahedronattributes" << out.numberoftetrahedronattributes  << endl;
    exit(1);
   }
  */
}


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

void mesh3_tetgenio_out(const tetgenio &out, const int & label_tet, const int & label_face, Mesh3 & Th3)
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

  if(verbosity) cout << &out.trifacemarkerlist << endl;
  
  for(int ibe=0; ibe < Th3.nbe; ibe++){
    int iv[3];
    iv[0] = out.trifacelist[3*ibe]-1;
    iv[1] = out.trifacelist[3*ibe+1]-1;
    iv[2] = out.trifacelist[3*ibe+2]-1; 
    Th3.be(ibe).set( Th3.vertices, iv, label_face);
  }

}
// verison Mesh3 *

Mesh3 * mesh3_tetgenio_out(const tetgenio &out)
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
  //Th3.set(out.numberofpoints, out.numberoftetrahedra, out.numberoftrifaces);

  // new parameter
  if( out.numberoftetrahedronattributes !=  1){
    cout << "out.numberoftetrahedronattributes" << out.numberoftetrahedronattributes << endl;
  }
   
  Vertex3 *v = new Vertex3[out.numberofpoints];
  Tet *t  = new Tet[out.numberoftetrahedra];
  Tet *tt = t;
  Triangle3 *b  = new Triangle3[out.numberoftrifaces];
  Triangle3 *bb = b;

  i=0;
  for(int nnv=0; nnv<out.numberofpoints; nnv++){
    v[nnv].x=out.pointlist[i];
    v[nnv].y=out.pointlist[i+1];
    v[nnv].z=out.pointlist[i+2];       
    v[nnv].lab=out.pointmarkerlist[nnv];
    i=i+3;    
  }

  // test pour la distance minimale entre les points 
//   {
//     double dist,dist1;
//     dist = 1000000000000.;
//     for(int nnv=0; nnv<out.numberofpoints; nnv++){
//       for(int nnv1=nnv+1; nnv1< out.numberofpoints; nnv1++){
// 	dist1=(v[nnv].x-v[nnv1].x)*(v[nnv].x-v[nnv1].x)+(v[nnv].y-v[nnv1].y)*(v[nnv].y-v[nnv1].y)
// 	  +(v[nnv].z-v[nnv1].z)*(v[nnv].z-v[nnv1].z);
// 	dist=min(dist,sqrt(dist1));
// 	if( sqrt(dist1) < 1e-12){ 
// 	  cout << "point confondus" << nnv << "<--->" <<nnv1 << endl;
// 	  if( sqrt( pow(v[nnv].x,2)+pow(v[nnv].y,2)+pow(v[nnv].z,2) ) > 1e-10   )  cout << v[nnv] << " " << v[nnv1] << endl;
// 	}
//       }
//     }
//     cout << "dist entre les points du maillage tetgen" << dist << endl;
//   }
    
  i=0;
  for(int nnt=0; nnt < out.numberoftetrahedra; nnt++){
    int iv[4],lab;
    iv[0] = out.tetrahedronlist[i]-1;
    iv[1] = out.tetrahedronlist[i+1]-1;
    iv[2] = out.tetrahedronlist[i+2]-1;
    iv[3] = out.tetrahedronlist[i+3]-1;
    //lab   = label_tet;
    for(int jj=0; jj<4; jj++){
      assert( iv[jj] >=0 && iv[jj] < out.numberofpoints );
    }
    
    //cout << "nnt= " <<  nnt << " " << lab  << " " << out.tetrahedronattributelist[nnt] << endl;
    lab =  out.tetrahedronattributelist[nnt];
    //cout << "nnt= " <<  lab  << " " << out.tetrahedronattributelist[nnt] << endl;
    //Th3.elements[nnt].set( Th3.vertices, iv, lab);
    (*tt++).set( v, iv, lab);   
    i=i+4;
  }
 
  for(int ibe=0; ibe < out.numberoftrifaces; ibe++){
    int iv[3];
    iv[0] = out.trifacelist[3*ibe]-1;
    iv[1] = out.trifacelist[3*ibe+1]-1;
    iv[2] = out.trifacelist[3*ibe+2]-1; 
    for(int jj=0; jj<3; jj++){
      if(iv[jj]>= out.numberofpoints || iv[jj]< 0 ) cout << "iv[jj]=" << iv[jj] << " triangle" << ibe << endl;
      assert( iv[jj] >=0 && iv[jj] < out.numberofpoints );
    }
    //Th3.be(ibe).set( Th3.vertices, iv, out.trifacemarkerlist[ibe]);
    (*bb++).set( v, iv, out.trifacemarkerlist[ibe]);   
  }

  Mesh3 *T_TH3 = new Mesh3(out.numberofpoints, out.numberoftetrahedra, out.numberoftrifaces, v, t, b);
  cout << "FreeFem++: Check mesh given by tetgen" << endl;
  //return T_TH3;

  if( TestElementMesh3(*T_TH3) != 1){  
    return T_TH3;
  }
  else{
    //Mesh3 *T2_TH3 = TestElementMesh3_patch( *T_TH3 );
    //return T2_TH3;
    exit(1);
  }

}


Mesh3* mesh3_tetgenio_out(const tetgenio &out, const int & label_tet)
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
  //Th3.set(out.numberofpoints, out.numberoftetrahedra, out.numberoftrifaces);
   
  Vertex3 *v = new Vertex3[out.numberofpoints];
  Tet *t  = new Tet[out.numberoftetrahedra];
  Tet *tt = t;
  Triangle3 *b  = new Triangle3[out.numberoftrifaces];
  Triangle3 *bb = b;
    
   
  i=0;
  for(int nnv=0; nnv < out.numberofpoints ; nnv++){
    v[nnv].x=out.pointlist[i];
    v[nnv].y=out.pointlist[i+1];
    v[nnv].z=out.pointlist[i+2];       
    v[nnv].lab=out.pointmarkerlist[nnv];
    i=i+3;    
  }
    
  i=0;
  for(int nnt=0; nnt < out.numberoftetrahedra; nnt++){
    int iv[4],lab;
    iv[0] = out.tetrahedronlist[i]-1;
    iv[1] = out.tetrahedronlist[i+1]-1;
    iv[2] = out.tetrahedronlist[i+2]-1;
    iv[3] = out.tetrahedronlist[i+3]-1;
    lab   = label_tet;
    //lab = out.tetrahedronmarkerlist[nnt];
    //Th3.elements[nnt].set( Th3.vertices, iv, lab);
    (*tt++).set( v, iv, lab);   
    i=i+4;
  }
 
  for(int ibe=0; ibe < out.numberoftrifaces ; ibe++){
    int iv[3];
    iv[0] = out.trifacelist[3*ibe]-1;
    iv[1] = out.trifacelist[3*ibe+1]-1;
    iv[2] = out.trifacelist[3*ibe+2]-1; 
    //Th3.be(ibe).set( Th3.vertices, iv, out.trifacemarkerlist[ibe]);
    (*bb++).set( v, iv, out.trifacemarkerlist[ibe]);  
  }
  Mesh3 *T_TH3 = new Mesh3(out.numberofpoints, out.numberoftetrahedra, out.numberoftrifaces, v, t, b);
  //return T_TH3;
  cout << "FreeFem++: Check mesh given by tetgen" << endl;
  if( TestElementMesh3(*T_TH3) != 1){  
    return T_TH3;
  }
  else{ 
    //cout << "patch pour tetgen " << endl;
    //Mesh3 *T2_TH3 = TestElementMesh3_patch( *T_TH3 );
    //return T2_TH3;
    exit(1);
  }
}

Mesh3* mesh3_tetgenio_out(const tetgenio &out, const int & label_tet, const int & label_face)
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
  //Th3.set(out.numberofpoints, out.numberoftetrahedra, out.numberoftrifaces);

  Vertex3 *v = new Vertex3[out.numberofpoints];
  Tet *t  = new Tet[out.numberoftetrahedra];
  Tet *tt = t;
  Triangle3 *b  = new Triangle3[out.numberoftrifaces];
  Triangle3 *bb = b;

  i=0;
  for(int nnv=0; nnv < out.numberofpoints; nnv++){
    v[nnv].x=out.pointlist[i];
    v[nnv].y=out.pointlist[i+1];
    v[nnv].z=out.pointlist[i+2];       
    v[nnv].lab=out.pointmarkerlist[nnv];
    i=i+3;    
  }
    
  i=0;
  for(int nnt=0; nnt < out.numberoftetrahedra; nnt++){
    int iv[4],lab;
    iv[0] = out.tetrahedronlist[i]-1;
    iv[1] = out.tetrahedronlist[i+1]-1;
    iv[2] = out.tetrahedronlist[i+2]-1;
    iv[3] = out.tetrahedronlist[i+3]-1;
    lab   = label_tet;
    //lab = out.tetrahedronmarkerlist[nnt];
    //Th3.elements[nnt].set( Th3.vertices, iv, lab);
    (*tt++).set( v, iv, lab); 
    i=i+4;
  }

  if(verbosity) cout << &out.trifacemarkerlist << endl;
  
  for(int ibe=0; ibe < out.numberoftrifaces ; ibe++){
    int iv[3];
    iv[0] = out.trifacelist[3*ibe]-1;
    iv[1] = out.trifacelist[3*ibe+1]-1;
    iv[2] = out.trifacelist[3*ibe+2]-1; 
    //Th3.be(ibe).set( Th3.vertices, iv, label_face);
    (*bb++).set( v, iv, label_face);  
  }

  Mesh3 *T_TH3 = new Mesh3(out.numberofpoints, out.numberoftetrahedra, out.numberoftrifaces, v, t, b);

  if( TestElementMesh3( *T_TH3 ) != 1){  
    return T_TH3;
  }
  else{
    exit(1);
    //Mesh3 *T2_TH3 = TestElementMesh3_patch( *T_TH3 );
    //return T2_TH3;
  }
}



Mesh3 * Convexhull_3Dpoints(char *switch_tetgen, const int &nv_t, const double *Xcoord, const double *Ycoord, 
const double *Zcoord, const int &label_tet, const int &label_face){
	
  //Mesh3 *T_Th3= new Mesh3;
		
  tetgenio in,out;
  //tetgenio::facet *f;
  //tetgenio::polygon *p;

  cout << " tetgenio: vertex " << endl;
  int itet,jtet;

  in.firstnumber = 1;
  in.numberofpoints = nv_t;
  in.pointlist = new REAL[in.numberofpoints*3];
  in.pointmarkerlist = new int[in.numberofpoints];
  itet=0;
  jtet=0;
  for(int nnv=0; nnv < nv_t; nnv++)
    {
      in.pointlist[itet]   = Xcoord[nnv];
      in.pointlist[itet+1] = Ycoord[nnv];
      in.pointlist[itet+2] = Zcoord[nnv];       
      in.pointmarkerlist[nnv] =  0;
      
      itet=itet+3;
    }
  assert(itet==in.numberofpoints*3);
	
  in.numberoffacets = 0;

  cout << "tetgen: before tetrahedralize( , &in, &out): switch=" << switch_tetgen<< endl;
  tetrahedralize(switch_tetgen, &in, &out);
	 
  cout << "tetgen: finish tetrahedralize( , &in, &out);" << endl;
  //mesh3_tetgenio_out( out, label_tet, label_face,*T_Th3);
  Mesh3* T_Th3=mesh3_tetgenio_out( out, label_tet, label_face);
  cout <<" Finish Mesh3 tetgen :: Vertex, Element, Border" << T_Th3->nv << " "<< T_Th3->nt << " " << T_Th3->nbe << endl;

  cout << "FreeFem++: End check mesh given by tetgen" << endl;  
  return T_Th3; 
}


Mesh3 * RemplissageSurf3D_tetgen(char *switch_tetgen,const Mesh3 & Th3, const int & label_tet){
	
  //Mesh3 *T_Th3= new Mesh3;
	
  assert(Th3.nt == 0 );
  int nv_t = Th3.nv;
  int nt_t = Th3.nt;
  int nbe_t = Th3.nbe;
	
  if(verbosity) cout << "3D RemplissageSurf3D:: Vertex  triangle2  border "    << nv_t << " "<< nt_t << " " << nbe_t<< endl;
  

  // Creation des tableau de tetgen
	
  tetgenio in,out;
  //tetgenio::facet *f;
  //tetgenio::polygon *p;
	
  if(verbosity) cout << " tetgenio: vertex " << endl;
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
	
  if(verbosity) cout << " tetgenio: facet " << endl;
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
  cout << "tetgen: before tetrahedralize( , &in, &out);" << endl;
	
  tetrahedralize(switch_tetgen, &in, &out);
	 
  cout << "tetgen: after tetrahedralize( , &in, &out);" << endl;
  //mesh3_tetgenio_out( out, label_tet, *T_Th3);
  Mesh3 *T_Th3 = mesh3_tetgenio_out( out, label_tet);
  cout <<" Finish Mesh3 tetgen :: Vertex, Element, Border" << T_Th3->nv << " "<< T_Th3->nt << " " << T_Th3->nbe << endl;
  cout << "FreeFem++: End check mesh given by tetgen" << endl;  
  return T_Th3;
}

Mesh3 * RemplissageSurf3D_tetgen_new(char *switch_tetgen,const Mesh3 & Th3, const int & label_tet,
				     const int &nbhole, const double *tabhole, 
				     const int & nbregion, const double *tabregion, 
				     const int &nbfacecl, const double *tabfacecl){
	
  //Mesh3 *T_Th3= new Mesh3;
	
  assert(Th3.nt == 0 );
  int nv_t = Th3.nv;
  int nt_t = Th3.nt;
  int nbe_t = Th3.nbe;

 
  if(verbosity) cout << "3D RemplissageSurf3D:: Vertex  triangle2  border " << nv_t << " "<< nt_t << " " << nbe_t<< endl;
  // Creation des tableau de tetgen
	
  tetgenio in,out;
  //tetgenio::facet *f;
  //tetgenio::polygon *p;
	
  if(verbosity) cout << " tetgenio: vertex " << endl;
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
	
  if(verbosity) cout << " tetgenio: facet " << endl;
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
  // mise a jour des nouvelles variables
 
  in.numberofholes = nbhole;
  in.holelist = new REAL[3*nbhole];
  
  for(int ii=0; ii<3*in.numberofholes; ii++){
    in.holelist[ii]   = tabhole[ii];
  }

  in.numberofregions = nbregion;
  in.regionlist = new REAL[5*nbregion];
  for(int ii=0; ii<5*in.numberofregions; ii++){
    in.regionlist[ii] = tabregion[ii];
  }

  in.numberoffacetconstraints = nbfacecl;
  in.facetconstraintlist = new REAL[2*in.numberoffacetconstraints];

  for(int ii=0; ii<2*in.numberoffacetconstraints; ii++){
    in.facetconstraintlist[ii+1] = tabfacecl[ii+1];
  }


  cout << "tetgen: before tetrahedralize( , &in, &out);" << endl;
  cout << "numberof regions "<<  in.numberofregions  << endl;
  cout << "numberof hole "<<  in.numberofholes  << endl;
	
  tetrahedralize(switch_tetgen, &in, &out);
	 
  cout << "tetgen: after tetrahedralize( , &in, &out);" << endl;
  //mesh3_tetgenio_out( out, *T_Th3);
  Mesh3* T_Th3=mesh3_tetgenio_out( out);
  cout <<" Finish Mesh3 tetgen :: Vertex, Element, Border" << T_Th3->nv << " "<< T_Th3->nt << " " << T_Th3->nbe << endl;
  cout << "FreeFem++: End check mesh given by tetgen" << endl;  
  return T_Th3;
}

Mesh3 * RemplissageSurf3D_tetgen_new(char *switch_tetgen,const Mesh3 & Th3, const int & label_tet,
				     const int &nbhole, const double *tabhole, 
				     const int & nbregion, const double *tabregion, 
				     const int &nbfacecl, const double *tabfacecl, 
				     const int &nbinside, const double *InsidePoint, 
				     const int &sizeofmetric, const double *metric){
	
  //Mesh3 *T_Th3= new Mesh3;
	
  assert(Th3.nt == 0 );
  int nv_t = Th3.nv;
  int nt_t = Th3.nt;
  int nbe_t = Th3.nbe;

 
  if(verbosity) cout << "3D RemplissageSurf3D:: Vertex  triangle2  border " << nv_t << " "<< nt_t << " " << nbe_t<< endl;
  // Creation des tableau de tetgen
	
  tetgenio in,out;
  tetgenio addin;

  if(verbosity) cout << " tetgenio: vertex " << endl;
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
      in.pointlist[itet]      = Th3.vertices[nnv].x;
      in.pointlist[itet+1]    = Th3.vertices[nnv].y;
      in.pointlist[itet+2]    = Th3.vertices[nnv].z;       
      in.pointmarkerlist[nnv] = Th3.vertices[nnv].lab;   
      itet=itet+3;
    }
  assert(itet==in.numberofpoints*3);
   
  // Add inside point
  if( nbinside ){
    cout << "nbinside=" << nbinside << endl;
    addin.firstnumber = 1;
    addin.numberofpoints = nbinside;
    addin.pointlist= new REAL[3*nbinside];
    addin.pointmarkerlist = new int[addin.numberofpoints];
    for(int nnv=0; nnv < 3*nbinside; nnv++)
      addin.pointlist[nnv] = InsidePoint[nnv];
    for(int nnv=0; nnv < nbinside; nnv++)
      addin.pointmarkerlist[nnv] = 111;
  }

  // Add metric
  if( sizeofmetric ){
    cout << "sizeofmetric=" << sizeofmetric << endl;
    in.numberofpointmtrs = sizeofmetric;
    in.pointmtrlist = new REAL[in.numberofpointmtrs*in.numberofpoints];
    for(int nnv=0; nnv < in.numberofpointmtrs*in.numberofpoints; nnv++)
      in.pointmtrlist[nnv]=metric[nnv]; 
  }
 
 if(verbosity) cout << " tetgenio: facet " << endl;
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
  // mise a jour des nouvelles variables
 
  in.numberofholes = nbhole;
  in.holelist = new REAL[3*nbhole];
  
  for(int ii=0; ii<3*in.numberofholes; ii++){
    in.holelist[ii]   = tabhole[ii];
  }

  in.numberofregions = nbregion;
  in.regionlist = new REAL[5*nbregion];
  for(int ii=0; ii<5*in.numberofregions; ii++){
    in.regionlist[ii] = tabregion[ii];
  }

  in.numberoffacetconstraints = nbfacecl;
  in.facetconstraintlist = new REAL[2*in.numberoffacetconstraints];

  for(int ii=0; ii<2*in.numberoffacetconstraints; ii++){
    in.facetconstraintlist[ii+1] = tabfacecl[ii+1];
  }


  cout << "tetgen: before tetrahedralize( , &in, &out);" << endl;
  cout << "numberof regions "<<  in.numberofregions  << endl;
  cout << "numberof hole "<<  in.numberofholes  << endl;
	
  tetrahedralize(switch_tetgen, &in, &out, &addin);
	 
  cout << "tetgen: after tetrahedralize( , &in, &out);" << endl;
  //mesh3_tetgenio_out( out, *T_Th3);
  Mesh3* T_Th3=mesh3_tetgenio_out( out);
  cout <<" Finish Mesh3 tetgen :: Vertex, Element, Border" << T_Th3->nv << " "<< T_Th3->nt << " " << T_Th3->nbe << endl;
  cout << "FreeFem++: End check mesh given by tetgen" << endl;  
  return T_Th3;
}


Mesh3 * Transfo_Mesh2_tetgen(const double &precis_mesh, char *switch_tetgen,const Mesh & Th2, const double *tab_XX, const double *tab_YY, const double *tab_ZZ, 
	int &border_only, int &recollement_border, int &point_confondus_ok, 
	const int &label_tet, const map<int, int> &maptri ){
  
  //Mesh3 *T_Th3= new Mesh3;
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
	
        if(verbosity) cout << "2D: Mesh::Vertex  triangle2  border " << Th2.nv << " "<<Th2.nt<< " " << Th2.neb<< endl;
	
	for(int ii=0; ii<Th2.nv; ii++){ 
		Numero_Som[ii]=ii;
	}
	if(verbosity) cout <<" debut: SamePointElement " <<endl;
	
	SamePointElement_Mesh2( precis_mesh, tab_XX, tab_YY, tab_ZZ, Th2, recollement_border, point_confondus_ok, 
		Numero_Som, ind_nv_t, ind_nt_t, ind_nbe_t, label_nbe_t, nv_t, nt_t, nbe_t);
	
	if(verbosity) cout <<" fin: SamePointElement " <<endl;
	
	if(verbosity) cout << "2D transfo: Mesh::Vertex  triangle2  border " << nv_t << " "<< nt_t << " " << nbe_t<< endl;
	// Creation des tableau de tetgen
	
	
	tetgenio in,out;
	//tetgenio::facet *f;
	//tetgenio::polygon *p;
	
	if(verbosity) cout << " tetgenio: vertex " << endl;
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
	
	if(verbosity) cout << " tetgenio: facet " << endl;
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

	cout << "tetgen: before tetrahedralize( , &in, &out);" << endl;
	
	tetrahedralize(switch_tetgen, &in, &out);
	 
	cout << "tetgen: after tetrahedralize( , &in, &out);" << endl;
	//mesh3_tetgenio_out( out, label_tet, *T_Th3);
	Mesh3 *T_Th3=mesh3_tetgenio_out( out, label_tet);
	cout <<" Finish Mesh3 :: Vertex, Element, Border" << T_Th3->nv << " "<< T_Th3->nt << " " << T_Th3->nbe << endl;
	
	

	delete [] Numero_Som;
	delete [] ind_nv_t;
	delete [] ind_nbe_t;
	delete [] label_nbe_t;

	cout << "FreeFem++: End check mesh given by tetgen" << endl;  
	return T_Th3;


}

Mesh3 * Transfo_Mesh2_tetgen_new(const double &precis_mesh, char *switch_tetgen,const Mesh & Th2, const double *tab_XX, const double *tab_YY, const double *tab_ZZ, 
				 int &border_only, int &recollement_border, int &point_confondus_ok, 
				 const int &label_tet, const map<int, int> &maptri, 
				 const int &nbhole, const double *tabhole, const int & nbregion, const double *tabregion, 
				 const int &nbfacecl, const double *tabfacecl){
  //Mesh3 *T_Th3= new Mesh3;
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
	
  if(verbosity) cout << "2D: Mesh::Vertex  triangle2  border " << Th2.nv << " "<<Th2.nt<< " " << Th2.neb<< endl;
	
  for(int ii=0; ii<Th2.nv; ii++){ 
    Numero_Som[ii]=ii;
  }
  if(verbosity) cout <<" debut: SamePointElement " <<endl;
	
  SamePointElement_Mesh2( precis_mesh, tab_XX, tab_YY, tab_ZZ, Th2, recollement_border, point_confondus_ok, 
			  Numero_Som, ind_nv_t, ind_nt_t, ind_nbe_t, label_nbe_t, nv_t, nt_t, nbe_t);
	
  if(verbosity) cout <<" fin: SamePointElement " <<endl;
	
  if(verbosity) cout << "2D transfo: Mesh::Vertex  triangle2  border " << nv_t << " "<< nt_t << " " << nbe_t<< endl;
  // Creation des tableau de tetgen
	
	
  tetgenio in,out;
  //tetgenio::facet *f;
  //tetgenio::polygon *p;
	
  if(verbosity) cout << " tetgenio: vertex " << endl;
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
	
  if(verbosity) cout << " tetgenio: facet " << endl;
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

  // mise a jour des nouvelles variables

  in.numberofholes = nbhole;
  in.holelist = new REAL[3*nbhole];
  
  for(int ii=0; ii<3*in.numberofholes; ii++){
    in.holelist[ii]   = tabhole[ii];
  }


  in.numberofregions = nbregion;
  in.regionlist = new REAL[5*nbregion];
  for(int ii=0; ii<5*in.numberofregions; ii++){
    in.regionlist[ii] = tabregion[ii];
  }

  in.numberoffacetconstraints = nbfacecl;
  in.facetconstraintlist = new REAL[2*in.numberoffacetconstraints];

  for(int ii=0; ii<2*in.numberoffacetconstraints; ii++){
    in.facetconstraintlist[ii+1] = tabfacecl[ii+1];
  }
  
  cout << "tetgen: before tetrahedralize( , &in, &out);" << endl;
  
  tetrahedralize(switch_tetgen, &in, &out);
  
  cout << "tetgen: after tetrahedralize( , &in, &out);" << endl;
  //mesh3_tetgenio_out( out, *T_Th3);
  Mesh3 *T_Th3=mesh3_tetgenio_out( out);
  cout <<" Finish Mesh3 :: Vertex, Element, Border" << T_Th3->nv << " "<< T_Th3->nt << " " << T_Th3->nbe << endl;

  delete [] Numero_Som;
  delete [] ind_nv_t;
  delete [] ind_nbe_t;
  delete [] label_nbe_t;
  
  cout << "FreeFem++: End check mesh given by tetgen" << endl;  
  return T_Th3;
}

// Fonction Refine avec tetgen

Mesh3 * ReconstructionRefine_tetgen(char *switch_tetgen,const Mesh3 & Th3, 
				     const int &nbhole, const double *tabhole, 
				     const int & nbregion, const double *tabregion, 
				     const int &nbfacecl, const double *tabfacecl, const double *tsizevol){

  // verif option refine
  int i;	
  assert(Th3.nt != 0 );
  {
    size_t testr, testp;
    int lenswitch;
    const char* test_tetgen = switch_tetgen;
    
    testr = strcspn(test_tetgen,"r");
    testp = strcspn(test_tetgen,"p");
    
    if( testr == strlen(test_tetgen) )
      { 
	cout << "The option 'r' of tetgen is not used" << endl;
	exit(1);
      }
    testp = strcspn(test_tetgen,"p");
    if( testp != strlen(test_tetgen) )
      { 
	cout << "With TetGen :: the option 'p' is not possible to use with option 'r' " << endl;
	exit(1);
      }
  }

  int nv_t = Th3.nv;
  int nt_t = Th3.nt;
  int nbe_t = Th3.nbe;
	
  if(verbosity) cout << "3D RemplissageSurf3D:: Vertex  triangle2  border "  << nv_t << " "<< nt_t << " " << nbe_t<< endl;
  // Creation des tableau de tetgen
	
  tetgenio in,out;
  //tetgenio::facet *f;
  //tetgenio::polygon *p;
	
  if(verbosity) cout << " tetgenio: vertex " << endl;
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

  // Tetrahedrons
  if(verbosity) cout << "tetrahedrons" << endl;
  in.numberofcorners = 4;
  in.numberoftetrahedra = Th3.nt;
  in.tetrahedronlist = new int[in.numberofcorners*in.numberoftetrahedra];
  in.numberoftetrahedronattributes = 1;
  in.tetrahedronattributelist = new REAL[in.numberoftetrahedronattributes*in.numberoftetrahedra];
  
  in.tetrahedronvolumelist = new REAL[in.numberoftetrahedra];

  i=0;
  for(int nnt=0; nnt < Th3.nt; nnt++){
    const Tet & K(Th3.elements[nnt]);

    in.tetrahedronlist[i]   = Th3.operator()(K[0])+1;
    in.tetrahedronlist[i+1] = Th3.operator()(K[1])+1;
    in.tetrahedronlist[i+2] = Th3.operator()(K[2])+1;
    in.tetrahedronlist[i+3] = Th3.operator()(K[3])+1;

    in.tetrahedronvolumelist[nnt] = tsizevol[nnt];
    in.tetrahedronattributelist[nnt] = K.lab;
    
    i=i+4;
  }

  
  if(verbosity) cout << "lecture des facettes" << endl;
  in.numberoftrifaces = Th3.nbe;
  in.trifacelist = new int[3*in.numberoftrifaces];
  in.trifacemarkerlist = new int[in.numberoftrifaces];

  for(int ibe=0; ibe < Th3.nbe; ibe++){
    const Triangle3 & K(Th3.be(ibe));
    
    in.trifacelist[3*ibe]   = Th3.operator()(K[0])+1;
    in.trifacelist[3*ibe+1] = Th3.operator()(K[1])+1;
    in.trifacelist[3*ibe+2] = Th3.operator()(K[2])+1; 					      
    in.trifacemarkerlist[ibe] = K.lab;

  }
  
  // mise a jour des nouvelles variables
 
  in.numberofholes = nbhole;
  in.holelist = new REAL[3*nbhole];
  
  for(int ii=0; ii<3*in.numberofholes; ii++){
    in.holelist[ii]   = tabhole[ii];
    if(verbosity) cout << "in.holelist[ii]=" << in.holelist[ii] << endl;
  }

  in.numberofregions = nbregion;
  in.regionlist = new REAL[5*nbregion];
  for(int ii=0; ii<5*in.numberofregions; ii++){
    in.regionlist[ii] = tabregion[ii];
    if(verbosity) cout << "in.regionlist[ii]=" << in.regionlist[ii] << endl;
  }

  in.numberoffacetconstraints = nbfacecl;
  in.facetconstraintlist = new REAL[2*in.numberoffacetconstraints];

  for(int ii=0; ii<2*in.numberoffacetconstraints; ii++){
    in.facetconstraintlist[ii+1] = tabfacecl[ii+1];
  }

  if(verbosity > 0){
    cout << "tetgen: before tetrahedralize( , &in, &out);" << endl;
    cout << "numberof regions "<<  in.numberofregions  << endl;
    cout << "numberof hole "<<  in.numberofholes  << endl;
  }
  tetrahedralize(switch_tetgen, &in, &out);

  if(verbosity > 0)
    cout << "tetgen: after tetrahedralize( , &in, &out);" << endl;
 
  Mesh3 *T_Th3=mesh3_tetgenio_out( out);
  if(verbosity > 0){
    cout <<" Finish Mesh3 tetgen :: Vertex, Element, Border" << T_Th3->nv << " "<< T_Th3->nt << " " << T_Th3->nbe << endl;
    cout << "FreeFem++: End check mesh given by tetgen" << endl;  
  }
  return T_Th3;
}


// Fonction Refine avec tetgen à l'aide d'une metrique

Mesh3 * ReconstructionRefine_tetgen(char *switch_tetgen,const Mesh3 & Th3, 
				     const int &nbhole, const double *tabhole, 
				     const int & nbregion, const double *tabregion, 
				     const int &nbfacecl, const double *tabfacecl, 
				     const double *tsizevol, const int &sizeofmetric, const double *metric){

  // verif option refine
  int i;	
  assert(Th3.nt != 0 );
  {
    size_t testr, testp;
    int lenswitch;
    const char* test_tetgen = switch_tetgen;
    
    testr = strcspn(test_tetgen,"r");
    testp = strcspn(test_tetgen,"p");
    
    if( testr == strlen(test_tetgen) )
      { 
	cout << "The option 'r' of tetgen is not used" << endl;
	exit(1);
      }
    testp = strcspn(test_tetgen,"p");
    if( testp != strlen(test_tetgen) )
      { 
	cout << "With TetGen :: the option 'p' is not possible to use with option 'r' " << endl;
	exit(1);
      }
  }

  int nv_t = Th3.nv;
  int nt_t = Th3.nt;
  int nbe_t = Th3.nbe;
	
  if(verbosity) cout << "3D RemplissageSurf3D:: Vertex  triangle2  border "  << nv_t << " "<< nt_t << " " << nbe_t<< endl;
  // Creation des tableau de tetgen
	
  tetgenio in,out;
  //tetgenio::facet *f;
  //tetgenio::polygon *p;
	
  if(verbosity) cout << " tetgenio: vertex " << endl;
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

  if( verbosity ) 
    cout << "sizeofmetric=" << sizeofmetric << endl;
  in.numberofpointmtrs = sizeofmetric;
  in.pointmtrlist = new REAL[in.numberofpointmtrs*in.numberofpoints];
  for(int nnv=0; nnv < in.numberofpointmtrs*in.numberofpoints; nnv++)
    in.pointmtrlist[nnv]=metric[nnv]; 

  // Tetrahedrons
  if(verbosity) cout << "tetrahedrons" << endl;
  in.numberofcorners = 4;
  in.numberoftetrahedra = Th3.nt;
  in.tetrahedronlist = new int[in.numberofcorners*in.numberoftetrahedra];
  in.numberoftetrahedronattributes = 1;
  in.tetrahedronattributelist = new REAL[in.numberoftetrahedronattributes*in.numberoftetrahedra];
  
  in.tetrahedronvolumelist = new REAL[in.numberoftetrahedra];

  i=0;
  for(int nnt=0; nnt < Th3.nt; nnt++){
    const Tet & K(Th3.elements[nnt]);

    in.tetrahedronlist[i]   = Th3.operator()(K[0])+1;
    in.tetrahedronlist[i+1] = Th3.operator()(K[1])+1;
    in.tetrahedronlist[i+2] = Th3.operator()(K[2])+1;
    in.tetrahedronlist[i+3] = Th3.operator()(K[3])+1;

    in.tetrahedronvolumelist[nnt] = tsizevol[nnt];
    in.tetrahedronattributelist[nnt] = K.lab;
    
    i=i+4;
  }

  
  if(verbosity) cout << "lecture des facettes" << endl;
  in.numberoftrifaces = Th3.nbe;
  in.trifacelist = new int[3*in.numberoftrifaces];
  in.trifacemarkerlist = new int[in.numberoftrifaces];

  for(int ibe=0; ibe < Th3.nbe; ibe++){
    const Triangle3 & K(Th3.be(ibe));
    
    in.trifacelist[3*ibe]   = Th3.operator()(K[0])+1;
    in.trifacelist[3*ibe+1] = Th3.operator()(K[1])+1;
    in.trifacelist[3*ibe+2] = Th3.operator()(K[2])+1; 					      
    in.trifacemarkerlist[ibe] = K.lab;

  }
  
  // mise a jour des nouvelles variables
 
  in.numberofholes = nbhole;
  in.holelist = new REAL[3*nbhole];
  
  for(int ii=0; ii<3*in.numberofholes; ii++){
    in.holelist[ii]   = tabhole[ii];
    if(verbosity) cout << "in.holelist[ii]=" << in.holelist[ii] << endl;
  }

  in.numberofregions = nbregion;
  in.regionlist = new REAL[5*nbregion];
  for(int ii=0; ii<5*in.numberofregions; ii++){
    in.regionlist[ii] = tabregion[ii];
    if(verbosity) cout << "in.regionlist[ii]=" << in.regionlist[ii] << endl;
  }

  in.numberoffacetconstraints = nbfacecl;
  in.facetconstraintlist = new REAL[2*in.numberoffacetconstraints];

  for(int ii=0; ii<2*in.numberoffacetconstraints; ii++){
    in.facetconstraintlist[ii+1] = tabfacecl[ii+1];
  }

  if(verbosity > 0){
    cout << "tetgen: before tetrahedralize( , &in, &out);" << endl;
    cout << "numberof regions "<<  in.numberofregions  << endl;
    cout << "numberof hole "<<  in.numberofholes  << endl;
  }
  tetrahedralize(switch_tetgen, &in, &out);

  if(verbosity > 0)
    cout << "tetgen: after tetrahedralize( , &in, &out);" << endl;
 
  Mesh3 *T_Th3=mesh3_tetgenio_out( out);
  if(verbosity > 0){
    cout <<" Finish Mesh3 tetgen :: Vertex, Element, Border" << T_Th3->nv << " "<< T_Th3->nt << " " << T_Th3->nbe << endl;
    cout << "FreeFem++: End check mesh given by tetgen" << endl;  
  }
  return T_Th3;
}

// declaration pour FreeFem++

class Remplissage_Op : public E_F0mps 
{
public:
  //typedef pmesh3 Result;
  Expression eTh;    // Surface mesh
  // ====================
  // This parameter allow to add inside points of this initial volume mesh
  Expression eVolTh;  
  bool bVol;
  // ====================
  static const int n_name_param =9+2+1+1; // 
  static basicAC_F0::name_and_type name_param[] ;
  Expression nargs[n_name_param];
  KN_<long>  arg(int i,Stack stack,KN_<long> a ) const
  { return nargs[i] ? GetAny<KN_<long> >( (*nargs[i])(stack) ): a;}
  KN_<double>  arg(int i,Stack stack,KN_<double> a ) const
  { return nargs[i] ? GetAny<KN_<double> >( (*nargs[i])(stack) ): a;}
  double  arg(int i,Stack stack,double a ) const{ return nargs[i] ? GetAny< double >( (*nargs[i])(stack) ): a;}
  long   arg(int i,Stack stack, long  a ) const{ return nargs[i] ? GetAny< long >( (*nargs[i])(stack) ): a;}
  string*  arg(int i,Stack stack, string* a ) const{ return nargs[i] ? GetAny< string* >( (*nargs[i])(stack) ): a;}
public:
  Remplissage_Op(const basicAC_F0 &  args,Expression tth) 
    : eTh(tth)
  {
    if(verbosity >1) cout << "Remplissage du bord" << endl;
    args.SetNameParam(n_name_param,name_param,nargs);
    if( nargs[2] && nargs[9] ) 
	CompileError("uncompatible movemesh3 (Th, region= , reftet=  ");
    if( nargs[3] && nargs[10] ) 
	CompileError("uncompatible movemesh3 (Th, label= , refface=  ");
    
    bVol=false;
    /*
    if( BCastTo<Mesh3 *>(args[1]) ){
      eVolTh = CastTo<Mesh3 *>(args[1]);
      bVol=true;
    }
    else{
      bVol=false;
    } 
    */
  }
  Remplissage_Op(const basicAC_F0 &  args,Expression tth, Expression vth) 
    : eTh(tth),eVolTh(vth)
  {
    if(verbosity >1) cout << "Remplissage du bord" << endl;
    args.SetNameParam(n_name_param,name_param,nargs);
    if( nargs[2] && nargs[9] ) 
	CompileError("uncompatible movemesh3 (Th, region= , reftet=  ");
    if( nargs[3] && nargs[10] ) 
	CompileError("uncompatible movemesh3 (Th, label= , refface=  ");
    
    bVol=true;
  }
  /*
    static ArrayOfaType  typeargs() { return  ArrayOfaType( atype<pmesh3>(),true ); }// all type
    static  E_F0 * f(const basicAC_F0 & args) { return new Remplissage_Op(args);} 
    operator aType () const { return atype<pmesh3>();} 
  */
  AnyType operator()(Stack stack)  const ;
};


basicAC_F0::name_and_type Remplissage_Op::name_param[]= {
  {  "switch", &typeid(string*)},
  {  "reftet", &typeid(long)}, //1
  {  "refface", &typeid(KN_<long> )},//2
  // new parmameters
  {  "nbofholes", &typeid(long)},
  {  "holelist", &typeid(KN_<double>)},
  {  "nbofregions", &typeid(long)},
  {  "regionlist", &typeid(KN_<double>)},
  {  "nboffacetcl", &typeid(long)},
  {  "facetcl", &typeid(KN_<double>)},
  {  "region", &typeid(long)}, //9
  {  "label", &typeid(KN_<long>)},//10
  {  "addpointlist", &typeid(KN_<long>)}, // 11
  {  "metric", &typeid(KN_<long>)}
};

class  Remplissage : public OneOperator { public:  
    Remplissage() : OneOperator(atype<pmesh3>(),atype<pmesh3>()) {}
  
  E_F0 * code(const basicAC_F0 & args) const 
  { 
	return  new Remplissage_Op( args,t[0]->CastTo(args[0]) ); 
  }
};


class  RemplissageAddPoint : public OneOperator { public:  
    RemplissageAddPoint() : OneOperator(atype<pmesh3>(),atype<pmesh3>(),atype<pmesh3>()) {}
  
  E_F0 * code(const basicAC_F0 & args) const 
  { 
    return  new Remplissage_Op( args,t[0]->CastTo(args[0]),t[1]->CastTo(args[1]) ); 
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
  cout << "Tetgen : Vertex Triangle Border " << nbv<< "  "<< nbt << " nbe "<< nbe << endl; 
 
  KN<long> zzempty;
  //int intempty=0;
  string stringempty= string("pqaAAYQC");
  string* switch_tet(arg(0,stack,&stringempty));
  int label_tet(arg(1,stack,arg(9,stack,0L)));  
  KN<long> nrf (arg(2,stack,arg(10,stack,zzempty)));

  // new parameters
  KN<double> zdzempty;
  int nbhole (arg(3,stack,0L));
  KN<double> tabhole (arg(4,stack,zdzempty));
  int nbregion (arg(5,stack,0L));
  KN<double> tabregion (arg(6,stack,zdzempty));
  int nbfacecl (arg(7,stack,0L));
  KN<double> tabfacecl (arg(8,stack,zdzempty));
  // parameter inside point
  // need to add "i" to the switch
  KN<double> InsidePoint(arg(11,stack,zdzempty)); // Add inside point in the volume mesh generated by tetgen
  // need to add "m" to the switch
  KN<double> metric(arg(12,stack,zdzempty)); // Add metric for tetgen

  //=========================
  // Add  a metric
  int sizeofmetric= metric.N()/Th.nv;
  if( nargs[12] ){
    cout << " size of the metric " << metric.N()/Th.nv << endl;
    assert( (metric.N()/Th.nv)*Th.nv == metric.N() );
  }
  // fin add a metric
  //==========================
  
  //==========================
  // Add inside points

  if( nargs[11] )
    assert( ((InsidePoint.N()/3)*3) == InsidePoint.N() );
  // case with a inside meshes

  if( bVol ){
    // Inside point is given by a mesh
    Mesh3 *pvolTh = GetAny<Mesh3 *>((*eVolTh)(stack));
    Mesh3 &volTh=*pvolTh;
    
    KN<int> takevertex(volTh.nv);
    takevertex = 1;
    // determination of vertices in the border
    for(int ibe=0; ibe<volTh.nbe; ibe++){
      const Triangle3 & K(volTh.be(ibe));
      takevertex[volTh.operator()(K[0])] = 0;
      takevertex[volTh.operator()(K[1])] = 0;
      takevertex[volTh.operator()(K[2])] = 0;
    }
    int nvInside=0;
    // number of vertices inside the volume mesh
    for(int iv=0; iv<volTh.nv; iv++){ 
      if( takevertex[iv] == 1)
	nvInside++;      
    }
    InsidePoint.resize(3*nvInside);
    int loopnv=0;
    for(int iv=0; iv<volTh.nv; iv++){ 
      if( takevertex[iv] == 1){
	InsidePoint[loopnv]   = volTh.vertices[iv].x;
	InsidePoint[loopnv+1] = volTh.vertices[iv].y;
	InsidePoint[loopnv+2] = volTh.vertices[iv].z;
	loopnv=loopnv+3;
      }      
    }
    assert( loopnv/3 == nvInside );
  }
  
  if( !bVol && !nargs[11] ) assert( 	InsidePoint.N() == 0 ); 

  //  fin add inisde point
  //=========================
    
    if(nbregion==0) nbregion=tabregion.N()/5;
    if(nbhole==0) nbhole=tabhole.N()/3;
    if(nbfacecl==0) nbfacecl=tabfacecl.N()/2;
    

  // assertion au niveau de la taille 
  ffassert( tabhole.N()   == 3*nbhole);
  ffassert( tabregion.N() == 5*nbregion);
  ffassert( tabfacecl.N() == 2*nbfacecl);

  

  //====================================
  //  How to change string* into char* 
  //====================================
  cout << "string" << *switch_tet << endl; 
  size_t size_switch_tet = switch_tet->size()+1;
  char* switch_tetgen =new char[size_switch_tet];
  strncpy(switch_tetgen, switch_tet->c_str(), size_switch_tet); 
   
  cout << "char" << switch_tetgen << endl; 
  
  ffassert( nrf.N() %2 ==0);

  map<int,int> mapf;
  for(int i=0;i<nrf.N();i+=2)
    {
      if(nrf[i] != nrf[i+1]){
	mapf[nrf[i]]=nrf[i+1];
      }
    }
  
  if(verbosity>1) cout << "tetgen:" << "nbhole="   << nbhole << "nbregion=" << nbregion << endl;

  /*
    int addcheckorientation=0;
    if( addcheckorientation==1  ){
    cout << "check :: orientation des surfaces" << endl;
    Th.BuildBoundaryElementAdj();
    cout << "fin check :: orientation des surfaces" << endl;
    }
  */

  int nbinside=InsidePoint.N()/3;
  Mesh3 *Th3 = 0;

  if( nargs[11] || nargs[12] || bVol){
    Th3 = RemplissageSurf3D_tetgen_new( switch_tetgen, Th, label_tet, nbhole, tabhole, nbregion, tabregion, nbfacecl, tabfacecl, nbinside, InsidePoint, sizeofmetric, metric);
    // delete multiple vertex
    Th3->TrueVertex();
  }
  else
    Th3 = RemplissageSurf3D_tetgen_new( switch_tetgen, Th, label_tet, nbhole, tabhole, nbregion, tabregion, nbfacecl, tabfacecl);
  

  cout << "finish tetgen " << endl;
  // changement de label 
  if( nrf.N() > 0){
    cout << "changement de label" << endl;
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
  cout << "action sur le maillage" << endl;

 // Th3->BuildBound();
 // Th3->BuildAdj();
 // Th3->Buildbnormalv();  
 // Th3->BuildjElementConteningVertex();
  Th3->BuildGTree();
  //Th3->decrement();    
  Add2StackOfPtr2FreeRC(stack,Th3);
  
  *mp=mps;
  delete [] switch_tetgen; 
  cout << "FreeFem++: End check mesh given by tetgen" << endl;  
  return Th3;
}

// Refine et Recontruction 

class ReconstructionRefine_Op : public E_F0mps 
{
public:
  Expression eTh;
  static const int n_name_param =10+2+1; // 
  static basicAC_F0::name_and_type name_param[] ;
  Expression nargs[n_name_param];
  KN_<long>  arg(int i,Stack stack,KN_<long> a ) const
  { return nargs[i] ? GetAny<KN_<long> >( (*nargs[i])(stack) ): a;}
  KN_<double>  arg(int i,Stack stack,KN_<double> a ) const
  { return nargs[i] ? GetAny<KN_<double> >( (*nargs[i])(stack) ): a;}
  double  arg(int i,Stack stack,double a ) const{ return nargs[i] ? GetAny< double >( (*nargs[i])(stack) ): a;}
  long   arg(int i,Stack stack, long  a ) const{ return nargs[i] ? GetAny< long >( (*nargs[i])(stack) ): a;}
  string*  arg(int i,Stack stack, string* a ) const{ return nargs[i] ? GetAny< string* >( (*nargs[i])(stack) ): a;}
public:
  ReconstructionRefine_Op(const basicAC_F0 &  args,Expression tth) 
    : eTh(tth)
  {
    if(verbosity >1) cout << "ReconstructionRefine du bord"<< endl;
    args.SetNameParam(n_name_param,name_param,nargs);
    if( nargs[2] && nargs[10] ) 
	CompileError("uncompatible ... (Th, region= , reftet=  ");
    if( nargs[3] && nargs[11] ) 
	CompileError("uncompatible ... (Th, label= , refface=  ");
    
  } 
  
  AnyType operator()(Stack stack)  const ;
};


basicAC_F0::name_and_type ReconstructionRefine_Op::name_param[]= {
  {  "switch", &typeid(string*)},
  {  "reftet", &typeid(KN_<long>)}, 
  {  "refface", &typeid(KN_<long> )},
  // new parmameters
  {  "nbofholes", &typeid(long)},
  {  "holelist", &typeid(KN_<double>)},
  {  "nbofregions", &typeid(long)},
  {  "regionlist", &typeid(KN_<double>)},
  {  "nboffacetcl", &typeid(long)},
  {  "facetcl", &typeid(KN_<double>)},
  {  "sizeofvolume", &typeid(double)},
  {  "region", &typeid(KN_<long>)}, //10
  {  "label", &typeid(KN_<long>)},//11
  {  "metric", &typeid(KN_<double>)}//12  // parameter for tetgen
    
};

class  ReconstructionRefine : public OneOperator { public:  

  ReconstructionRefine() : OneOperator(atype<pmesh3>(),atype<pmesh3>()) {}
  
  E_F0 * code(const basicAC_F0 & args) const 
  { 
    return  new ReconstructionRefine_Op( args, t[0]->CastTo(args[0]) ); 
  }
};

AnyType ReconstructionRefine_Op::operator()(Stack stack)  const 
{
  MeshPoint *mp(MeshPointStack(stack)) , mps=*mp;
  Mesh3 * pTh= GetAny<Mesh3 *>((*eTh)(stack));
  //double msvol= GetAny<double>((*maxvol)(stack));
 
  ffassert( pTh );
  Mesh3 &Th=*pTh;
  Mesh3 *m= pTh;   // question a quoi sert *m ??
  int nbv=Th.nv; // nombre de sommet 
  int nbt=Th.nt; // nombre de triangles
  int nbe=Th.nbe; // nombre d'aretes fontiere
  cout << "refine tetgen: Vertex Triangle Border " << nbv<< "  "<< nbt << " "<< nbe << endl; 
 
  KN<long> zzempty;
  //int intempty=0;
  string stringempty= string("rqaAAYQC");
  string* switch_tet(arg(0,stack,&stringempty));
  KN<long> nrtet(arg(1,stack,arg(10,stack,zzempty)));  
  KN<long> nrf (arg(2,stack,arg(11,stack,zzempty)));

  // new parameters
  KN<double> zdzempty;
  int nbhole (arg(3,stack,0L));
  KN<double> tabhole (arg(4,stack,zdzempty));
  int nbregion (arg(5,stack,0L));
  KN<double> tabregion (arg(6,stack,zdzempty));
  int nbfacecl (arg(7,stack,0L));
  KN<double> tabfacecl (arg(8,stack,zdzempty));

  KN<double> metric(arg(12,stack,zdzempty)); // Add metric for tetgen

  //=========================
  // Add  a metric
  int sizeofmetric = metric.N()/Th.nv;
  if( nargs[12] ){
    cout << " size of the metric " << metric.N()/Th.nv << endl;
    assert( (metric.N()/Th.nv)*Th.nv == metric.N() );
  }
  // fin add a metric
  //==========================
    if(nbregion==0) nbregion=tabregion.N()/5;
    if(nbhole==0) nbhole=tabhole.N()/3;
    if(nbfacecl==0) nbfacecl=tabfacecl.N()/2;
    
  // assertion au niveau de la taille
    
  ffassert( tabhole.N()   == 3*nbhole);
  ffassert( tabregion.N() == 5*nbregion);
  ffassert( tabfacecl.N() == 2*nbfacecl);
  
  //====================================
  //  How to change string* into char* 
  //====================================
  size_t size_switch_tet = switch_tet->size()+1;
  char* switch_tetgen =new char[size_switch_tet];
  strncpy(switch_tetgen, switch_tet->c_str(), size_switch_tet); 
   
  ffassert( nrf.N() %2 ==0);
  map<int,int> mapf;
  for(int i=0;i<nrf.N();i+=2)
    {
      if(nrf[i] != nrf[i+1]){
	mapf[nrf[i]]=nrf[i+1];
      }
    }

  ffassert( nrtet.N() %2 ==0);
  map<int,int> maptet;
  for(int i=0;i<nrtet.N();i+=2)
    {
      if(nrtet[i] != nrtet[i+1]){
	maptet[nrtet[i]]=nrtet[i+1];
      }
    }

  KN<double> tsizevol(nbt);
  MeshPoint *mp3(MeshPointStack(stack)); 
	
  R3 Cdg_hat = R3(1./4.,1./4.,1./4.);

  for (int it=0; it<nbt; it++){
    Tet & K(Th.elements[it]);

    mp3->set( Th, K(Cdg_hat), Cdg_hat, K, K.lab);
    if( nargs[9]){
      tsizevol[it]=  GetAny< double >( (*nargs[9])(stack) );
    }
    else if(tabregion.N()==0){
      for(int i=0; i< nbregion;i++){
	if( K.lab == tabregion[3+5*i] ){
	  tsizevol[it] = tabregion[4+5*i];
	}
      }
      
    }
    else{
      tsizevol[it] = K.mesure();
    }
  }

  cout << "Before reconstruction:" << " nbhole="   << nbhole << " nbregion=" << nbregion << endl;
 
  Mesh3 *Th3=0; 
  
  int RefineMethod=-1;
  if( nargs[9] ) 
    RefineMethod=1;
  if( nargs[12] ) 
    RefineMethod=0;
  
  // Add parameter "perhaps' with add a metric which defined sizeofvolume

  if( RefineMethod == 1) 
    Th3 = ReconstructionRefine_tetgen(switch_tetgen, Th, nbhole, tabhole, nbregion, tabregion, nbfacecl,tabfacecl,tsizevol);
  else if( RefineMethod == 0) 
    Th3 = ReconstructionRefine_tetgen(switch_tetgen, Th, nbhole, tabhole, nbregion, tabregion, nbfacecl,tabfacecl,tsizevol,sizeofmetric,metric);
  else{
    cerr << " We can't refine the initial mesh with tetgen. No sizeofvolume or metric is given " << endl;
    exit(1);
  }

  cout << "finish reconstruction " << endl;
  // changement de label 1
 
  if( nrtet.N() > 0){
    for(int ii=0; ii< Th3->nt; ii++){
      const Tet & K(Th3->elements[ii]);
      int lab;
      int iv[4];
			
      iv[0] = Th3->operator()(K[0]);
      iv[1] = Th3->operator()(K[1]);
      iv[2] = Th3->operator()(K[2]);
      iv[3] = Th3->operator()(K[3]);

      map< int, int> :: const_iterator imap;
      imap = maptet.find(K.lab);
      if( imap != maptet.end() ){
	lab = imap -> second;
      }
      else{
	lab = K.lab;
      }
      Th3->elements[ii].set(Th3->vertices,iv,lab);
    }
  }

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
  //cout << "fin du changement de label " << endl;
  
//  Th3->BuildBound();
//  Th3->BuildAdj();
//  Th3->Buildbnormalv();  
//  Th3->BuildjElementConteningVertex();
  Th3->BuildGTree();
  //Th3->decrement();    
  Add2StackOfPtr2FreeRC(stack,Th3);

  delete [] switch_tetgen;
  *mp=mps;
  cout << "FreeFem++: End check mesh given by tetgen" << endl;  
  return Th3;
}



// ConvexHull3D_tetg_Op

class ConvexHull3D_tetg_Op : public E_F0mps 
{
public:
  Expression numofpts;
  Expression xx,yy,zz;
  static const int n_name_param =3+1; // 
  static basicAC_F0::name_and_type name_param[] ;
  Expression nargs[n_name_param];
  
  KN_<long>  arg(int i,Stack stack,KN_<long> a ) const
  { return nargs[i] ? GetAny<KN_<long> >( (*nargs[i])(stack) ): a;}
  KN_<double>  arg(int i,Stack stack,KN_<double> a ) const
  { return nargs[i] ? GetAny<KN_<double> >( (*nargs[i])(stack) ): a;}
  double  arg(int i,Stack stack,double a ) const{ return nargs[i] ? GetAny< double >( (*nargs[i])(stack) ): a;}
  long   arg(int i,Stack stack, long  a ) const{ return nargs[i] ? GetAny< long  >( (*nargs[i])(stack) ): a;}
  string*  arg(int i,Stack stack, string* a ) const{ return nargs[i] ? GetAny< string* >( (*nargs[i])(stack) ): a;}
  
public:
  ConvexHull3D_tetg_Op(const basicAC_F0 &  args, Expression nop, 
  Expression ffxx, Expression ffyy, Expression ffzz ) 
    : numofpts(nop), xx(ffxx), yy(ffyy), zz(ffzz)
  {
    if(verbosity) cout << "Convex Hull with TetGen" << endl;
    args.SetNameParam(n_name_param,name_param,nargs);
    if( nargs[2] && nargs[3] ) 
	CompileError("uncompatible ... (Th, region= , reftet=  ");
    if( nargs[3] && nargs[4] ) 
	CompileError("uncompatible ... (Th, label= , refface=  ");
    
  } 
  
  AnyType operator()(Stack stack)  const ;
};


basicAC_F0::name_and_type  ConvexHull3D_tetg_Op::name_param[]= {
  {  "switch", &typeid(string*)},
  {  "reftet", &typeid(long)}, 
  {  "refface", &typeid(long)},
    {  "region", &typeid(long)}, 
    {  "label", &typeid(long)}
};

class ConvexHull3D_tetg : public OneOperator { public:  
     ConvexHull3D_tetg() : OneOperator( atype<pmesh3>(), atype<long>(), 
     atype< KN_<double> >(), atype< KN_<double> >(), atype< KN_<double> >() ) {}
  
  E_F0 * code(const basicAC_F0 & args) const 
  { 
	return  new  ConvexHull3D_tetg_Op( args,t[0]->CastTo(args[0]), 
	      t[1]->CastTo(args[1]), t[2]->CastTo(args[2]), t[3]->CastTo(args[3]) ); 
  }
};

AnyType  ConvexHull3D_tetg_Op::operator()(Stack stack)  const 
{
  int nbv = (int) GetAny<long>((*numofpts)(stack));
  KN<double> cxx(nbv),cyy(nbv),czz(nbv);
  
  cxx = GetAny< KN<double> > ((*xx)(stack));
  cyy = GetAny< KN<double> > ((*yy)(stack));
  czz = GetAny< KN<double> > ((*zz)(stack));
  

  assert( cxx.N() == nbv );
  assert( cyy.N() == nbv );
  assert( czz.N() == nbv );
     
  KN<long> zzempty;
  //int intempty=0;
  string stringempty= string("YqaAAQC");
  string* switch_tet(arg(0,stack,&stringempty));
  int label_tet(arg(1,stack,arg(3,stack,0L)));  
  int label_face(arg(2,stack,arg(4,stack,1L)));

  //====================================
  //  How to change string* into char* 
  //====================================
  size_t size_switch_tet = switch_tet->size()+1;
  char* switch_tetgen =new char[size_switch_tet];
  strncpy(switch_tetgen, switch_tet->c_str(), size_switch_tet); 
  //======================================
 
  Mesh3 *Th3 =  Convexhull_3Dpoints ( switch_tetgen, nbv, cxx, cyy, czz, label_tet, label_face );
	
//  Th3->BuildBound();
//  Th3->BuildAdj();
//  Th3->Buildbnormalv();  
//  Th3->BuildjElementConteningVertex();
  Th3->BuildGTree();
  //Th3->decrement();   
  Add2StackOfPtr2FreeRC(stack,Th3);
  
  delete [] switch_tetgen;
  cout << "FreeFem++: End check mesh given by tetgen" << endl;  
  return Th3;
}

// ConvexHull3D_tetg_file_Op

class ConvexHull3D_tetg_file_Op: public E_F0mps 
{
public:
  Expression filename;
  static const int n_name_param =3; // 
  static basicAC_F0::name_and_type name_param[] ;
  Expression nargs[n_name_param];
  
  KN_<long>  arg(int i,Stack stack,KN_<long> a ) const
  { return nargs[i] ? GetAny<KN_<long> >( (*nargs[i])(stack) ): a;}
  KN_<double>  arg(int i,Stack stack,KN_<double> a ) const
  { return nargs[i] ? GetAny<KN_<double> >( (*nargs[i])(stack) ): a;}
  double  arg(int i,Stack stack,double a ) const{ return nargs[i] ? GetAny< double >( (*nargs[i])(stack) ): a;}
  long  arg(int i,Stack stack, long a ) const{ return nargs[i] ? GetAny< long >( (*nargs[i])(stack) ): a;}
  string*  arg(int i,Stack stack, string* a ) const{ return nargs[i] ? GetAny< string* >( (*nargs[i])(stack) ): a;}
  
public:
  ConvexHull3D_tetg_file_Op(const basicAC_F0 &  args, Expression zfilename ) 
    : filename(zfilename)
  {
    if(verbosity) cout << "Convex Hull with TetGen" << endl;
    args.SetNameParam(n_name_param,name_param,nargs);
    if( nargs[2] && nargs[3] ) 
	CompileError("uncompatible ... (Th, region= , reftet=  ");
    if( nargs[3] && nargs[4] ) 
	CompileError("uncompatible ... (Th, label= , refface=  ");
    
  } 
  
  AnyType operator()(Stack stack)  const ;
};


basicAC_F0::name_and_type  ConvexHull3D_tetg_file_Op::name_param[]= {
  {  "switch", &typeid(string*)},
  {  "reftet", &typeid(long)}, 
  {  "refface", &typeid(long)},
    {  "region", &typeid(long)}, 
    {  "label", &typeid(long)}
    
};

class ConvexHull3D_tetg_file : public OneOperator { public:  
     ConvexHull3D_tetg_file() : OneOperator( atype<pmesh3>(), atype<string *>() ) {}
  
  E_F0 * code(const basicAC_F0 & args) const 
  { 
	return  new  ConvexHull3D_tetg_file_Op( args,t[0]->CastTo(args[0]) ); 
  }
};

AnyType  ConvexHull3D_tetg_file_Op::operator()(Stack stack)  const 
{
  string*  pointsfile = GetAny<string*>((*filename)(stack));

  // lecture du fichier contenant les points
  int nbv;
  //int lec;
 
  ifstream fp( pointsfile->c_str() );

  if(!fp) {
    cerr << "  -- tetgconvexhull : Erreur openning " << pointsfile <<endl; 
    exit(1);
  }
  if(verbosity>1)
    cout << "  -- tetgconvexhull:  Read On file \"" << pointsfile <<"\""<<  endl;

  fp >> nbv;
 
  cout << "  -- Nb of Points " << nbv << endl;
  KN<double> cxx(nbv), cyy(nbv), czz(nbv);

  
  for(int lec=0;lec<nbv;lec++)
    {
      fp >> cxx[lec] >> cyy[lec] >> czz[lec];
    }

//  if( lec !=nbv ) {
//     cerr << "  -- tetgconvexhull : Erreur Reading File " << pointsfile <<endl; 
//     cerr << " number of points " << nbv << " " <<"number of lecture" << lec << endl;  
//     exit(1);
//   }
//   assert(lec==nbv);

  fp.close();    

  KN<long> zzempty;
  //int intempty=0;
  string stringempty= string("YQV");
  string* switch_tet(arg(0,stack,&stringempty));
  int label_tet(arg(1,stack,arg(3,stack,0L)));  
  int label_face(arg(2,stack,arg(4,stack,1L)));

  //====================================
  //  How to change string* into char* 
  //====================================
  size_t size_switch_tet = switch_tet->size()+1;
  char* switch_tetgen =new char[size_switch_tet];
  strncpy(switch_tetgen, switch_tet->c_str(), size_switch_tet); 
  //======================================
 
  Mesh3 *Th3 =new Mesh3; 
  
  Th3 = Convexhull_3Dpoints( switch_tetgen, nbv, cxx, cyy, czz, label_tet, label_face );
	
//  Th3->BuildBound();
//  Th3->BuildAdj();
//  Th3->Buildbnormalv();  
//  Th3->BuildjElementConteningVertex();
  Th3->BuildGTree();
  //Th3->decrement();    
  Add2StackOfPtr2FreeRC(stack,Th3);

  delete [] switch_tetgen;
  cout << "FreeFem++: End check mesh given by tetgen" << endl;  
  return Th3;
}


class Init1 { public:
  Init1();
};

LOADINIT(Init1)  //  une variable globale qui serat construite  au chargement dynamique 

Init1::Init1(){  // le constructeur qui ajoute la fonction "splitmesh3"  a freefem++ 
  
  //if (verbosity)
  if(verbosity) cout << " load: tetgen  " << endl;
  
  Global.Add("tetgconvexhull","(",new ConvexHull3D_tetg_file);
  Global.Add("tetgconvexhull","(",new ConvexHull3D_tetg);
  Global.Add("tetgtransfo","(",new Build2D3D);
  Global.Add("tetg","(",new Remplissage);
  Global.Add("tetg","(",new RemplissageAddPoint);
  Global.Add("tetgreconstruction","(",new ReconstructionRefine);

}
// because i include this file in tetgen.cpp (very bad) FH
// a will correct this in next version ...
#define  WITH_NO_INIT
#include "msh3.cpp" 
