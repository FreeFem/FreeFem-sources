// ORIG-DATE:     Aout 2008
// -*- Mode : c++ -*-
//
// SUMMARY  : liaison medit freefem++ : popen  
// USAGE    : LGPL      
// ORG      : LJLL Universite Pierre et Marie Curi, Paris,  FRANCE 
// AUTHOR   : Jacques Morice
// E-MAIL   : jacques.morice@ann.jussieu.fr
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

const char *medit_char= "meditff";
const char *medit_popen="-popen";// 1";  // depend de l endroit ou se trouve medit
const char *medit_addsol="-addsol";

//  const char *medit_charsol= "medit-2.3-win.exeCC -popen -addsol 1";

#include <iostream>
#include <cfloat>
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
#include "libmesh5.h"
#include "lgsolver.hpp"
#include "problem.hpp"
#include "LayerMesh.hpp"
#include "TransfoMesh_v2.hpp"
//#include "GQuadTree.hpp"

#include <set>
#include <vector>
#include <fstream>

//#include "Operator.hpp" 
//#include "array_resize.hpp"
//#include "lex.hpp"

//using Fem2D::Mesh;
//using Fem2D::MeshPoint;

//extern bool NoWait; 

using namespace std;
using namespace Fem2D;


class PopenMeditMesh_Op : public E_F0mps 
{
public:
  Expression eTh;
  Expression filename;
  Expression xx, yy; //expression  vector 
  Expression tsxx, tsyx, tsyy;  //expression sym tensor  
  static const int n_name_param =5;  
  static basicAC_F0::name_and_type name_param[] ;
  Expression nargs[n_name_param];
  long arg(int i,Stack stack,long a) const{ return nargs[i] ? GetAny< long >( (*nargs[i])(stack) ): a;}
public:
  PopenMeditMesh_Op(const basicAC_F0 &  args,Expression ffname,Expression tth) 
    : eTh(tth),filename(ffname),xx(0),yy(0),tsxx(0), tsyx(0), tsyy(0)
  {
    cout << "construction " << endl;
    args.SetNameParam(n_name_param,name_param,nargs);   
    const E_Array * a1=0 ;
    if(nargs[0])  a1  = dynamic_cast<const E_Array *>(nargs[2]);
    const E_Array * a2=0 ;
    if(nargs[0])  a2  = dynamic_cast<const E_Array *>(nargs[3]);
   
    if(a1) {
      if(a1->size() !=2) 
      CompileError("meditmesh (Th,vector=[fx,fy],) ");
      xx=to<double>( (*a1)[0]);
      yy=to<double>( (*a1)[1]);
    }    
    
    if(a2) {
      if(a2->size() !=3) 
      CompileError("meditmesh (Th,symtensor=[Txx,Tyx,Tyy],) ");
      tsxx=to<double>( (*a2)[0]);
      tsyx=to<double>( (*a2)[1]);
      tsyy=to<double>( (*a2)[2]);
    }   

    cout << "fin construction " << endl;
  } 
  
  AnyType operator()(Stack stack)  const ;
};

basicAC_F0::name_and_type PopenMeditMesh_Op::name_param[]= {
  {  "solution", &typeid(long)},
  {  "scalar", &typeid(double)},
  {  "vector", &typeid(E_Array)},
  {  "symtensor", &typeid(E_Array)},
  {  "order", &typeid(long)}
};

static char * meditcmd(int nbsol,const string &  ffnn)
{
  string  meditcmm=medit_char;
  meditcmm += ' ';
  meditcmm += medit_popen;
  if(nbsol) 
    {
      meditcmm += ' ';
      meditcmm += medit_addsol;
    }
  meditcmm += " 1 ";
  char * ret= new char[meditcmm.size()+1];
  strcpy( ret, meditcmm.c_str()); 
  return ret;
}
AnyType PopenMeditMesh_Op::operator()(Stack stack)  const 
{

  cout << "any type" << endl;
  MeshPoint *mp(MeshPointStack(stack)) , mps=*mp;
  Mesh * pTh= GetAny<Mesh *>((*eTh)(stack));
  ffassert(pTh);
  Mesh &Th=*pTh;
 
  int nv=Th.nv; // nombre de sommet 
  int nt=Th.nt; // nombre de triangles
  int nbe=Th.neb; // nombre d'aretes fontiere

  long sol  (arg(0,stack,0));
  long order (arg(4,stack,0));

  int ver = GmfFloat;
  int dimp =2;
  float fx,fy;
 
  //  const char *medit_char= "medit-2.3-win.exeCC -popen 1";  // depend de l endroit ou se trouve medit
  //  const char *medit_charsol= "medit-2.3-win.exeCC -popen -addsol 1";

  string * ffname= GetAny<string *>( (*filename)(stack) );
  char * commandline = meditcmd(sol,*ffname);
  long valsortie=0;
  int typsol,nbsol=sol;

  printf("version de medit %s\n",commandline);

  FILE *popenstream= popen(commandline,"wb");

  fprintf(popenstream,"MeshVersionFormatted\n");
  fprintf(popenstream,"%i\n",ver);
  fprintf(popenstream,"Dimension\n");
  fprintf(popenstream,"%i\n",dimp);
  fprintf(popenstream,"Vertices\n");
  fprintf(popenstream,"%i\n",nv);

  for (int k=0; k<nv; k++) {
    const  Mesh::Vertex & P = Th.vertices[k];
    fx=P.x; fy=P.y;
    fprintf(popenstream,"%f %f %i\n",fx,fy,P.lab);
  }
  
  fprintf(popenstream,"Triangles\n");
  fprintf(popenstream,"%i\n",nt);
  for (int k=0; k<nt; k++) {
    const Mesh::Triangle & K(Th.t(k));
    int i0=Th.operator()(K[0])+1;
    int i1=Th.operator()(K[1])+1;
    int i2=Th.operator()(K[2])+1;
    int lab=K.lab;
    fprintf(popenstream,"%i %i %i %i\n",i0,i1,i2,lab);
  }

  fprintf(popenstream,"Edges\n");
  fprintf(popenstream,"%i\n",nbe);
  for (int k=0; k<nbe; k++) {
    const Mesh::BorderElement & K(Th.be(k));
    int i0=Th.operator()(K[0])+1;
    int i1=Th.operator()(K[1])+1;
    int lab=K.lab;
    fprintf(popenstream,"%i %i %i\n",i0,i1,lab);
  }
  fprintf(popenstream,"End");

  if( sol > 0){
    fprintf(popenstream,"MeshVersionFormatted %i\n",ver);
    fprintf(popenstream,"Dimension %i\n",dimp);
    if(order==0){
      fprintf(popenstream,"SolAtTriangles\n");
      fprintf(popenstream,"%i\n",nv);
    }
    if(order==1){
      fprintf(popenstream,"SolAtVertices\n");
      fprintf(popenstream,"%i\n",nv);
    }
    // detemination of solution
    // default scalaire // faire tableau pour plusieurs
    
    nbsol = 1;
    typsol = 1;
    if ( xx ){
      typsol  = 2;
    }
    if( tsxx ){
      typsol  = 3;
    }
    fprintf(popenstream,"%i %i\n",nbsol,typsol);
       
    if(typsol==1){
      if(order==0){
	KN<double> solsca(nt);
	solsca=0.;
	MeshPoint *mp3(MeshPointStack(stack)); 
	R2 Cdg_hat = R2(1./3.,1./3.);  

	for (int it=0;it<Th.nt;++it){
	  const Mesh::Triangle & K(Th.t(it));
	  mp3->set( Th, K(Cdg_hat), Cdg_hat, K, K.lab);
	  solsca[it] =  GetAny< double >( (*nargs[1])(stack) );
	}
	
	for(int k=0; k<nt; k++){
	  fprintf(popenstream,"%f\n",solsca[k]);
	}
      }

      if(order==1){
	KN<double> solsca(nv);
	solsca=0.;
	KN<int> takemesh(nv);
	MeshPoint *mp3(MeshPointStack(stack)); 
	
	takemesh=0;
	for (int it=0;it<Th.nt;++it){
	  for( int iv=0; iv<3; ++iv){
	    int i=Th(it,iv);  
	    
	    if(takemesh[i]==0){
	      mp3->setP(&Th,it,iv);
	      solsca[i] =  GetAny< double >( (*nargs[1])(stack) );
	      
	      takemesh[i] = takemesh[i]+1;
	    }
	  }
	}
	
	for(int k=0; k<nv; k++){
	  fprintf(popenstream,"%f\n",solsca[k]);
	}
      }
    }
    else if(typsol==2){
      if(order==0){	
	KN<double>  vxx(nt),vyy(nt);
	MeshPoint *mp3(MeshPointStack(stack)); 
	R2 Cdg_hat = R2(1./3.,1./3.);  

	for (int it=0;it<Th.nt;++it){
	  const Mesh::Triangle & K(Th.t(it));
	  mp3->set( Th, K(Cdg_hat), Cdg_hat, K, K.lab);
	  vxx[it] =  GetAny< double >( (*xx)(stack) );
	  vyy[it] =  GetAny< double >( (*yy)(stack) );
	}
	
	for(int k=0; k<nt; k++){
	  fprintf(popenstream,"%f %f\n",vxx[k],vyy[k]);
	}
      }

      if(order==1){
	KN<double> vxx(nv),vyy(nv);
	KN<int> takemesh(nv);
	MeshPoint *mp3(MeshPointStack(stack)); 
	
	takemesh=0;
	for (int it=0;it<Th.nt;++it){
	  for( int iv=0; iv<3; ++iv){
	    int i=Th(it,iv);  
	    
	    if(takemesh[i]==0){
	      mp3->setP(&Th,it,iv);
	      vxx[i] = GetAny< double >( (*xx)(stack) );
	      vyy[i] = GetAny< double >( (*yy)(stack) );
	      
	      takemesh[i] = takemesh[i]+1;
	    }
	  }
	}
	for(int k=0; k<nv; k++){
	  fprintf(popenstream,"%f %f %f\n",vxx[k],vyy[k]);
	}
      }
    }
    else if(typsol==3){
      if(order==0){
      	KN<double>  vxx(nt),vyx(nt),vyy(nt);
	MeshPoint *mp3(MeshPointStack(stack)); 
	R2 Cdg_hat = R2(1./3.,1./3.);  

	for (int it=0;it<Th.nt;++it){
	  const Mesh::Triangle & K(Th.t(it));
	  mp3->set( Th, K(Cdg_hat), Cdg_hat, K, K.lab);
	  vxx[it] =  GetAny< double >( (*tsxx)(stack) );
	  vyx[it] =  GetAny< double >( (*tsyx)(stack) );
	  vyy[it] =  GetAny< double >( (*tsyy)(stack) );
	}
	
	for(int k=0; k<nt; k++){
	  fprintf(popenstream,"%f %f\n",vxx[k],vyx[k],vyy[k]);
	}
      }
      if(order==1){
	KN<double> vxx(nv),vyx(nv),vyy(nv);
	KN<int> takemesh(nv);
	MeshPoint *mp3(MeshPointStack(stack)); 
	
	takemesh=0;
	for (int it=0;it<Th.nt;++it){
	  for( int iv=0; iv<3; ++iv){
	    int i=Th(it,iv);  
	    
	    if(takemesh[i]==0){
	      mp3->setP(&Th,it,iv);
	      vxx[i] = GetAny< double >( (*tsxx)(stack) );
	      vyx[i] = GetAny< double >( (*tsyx)(stack) );
	      vyy[i] = GetAny< double >( (*tsyy)(stack) );
	    
	      takemesh[i] = takemesh[i]+1;
	    }
	  }
	}
      
	for(int k=0; k<nv; k++){
	  fprintf(popenstream,"%f %f %f\n",vxx[k],vyx[k],vyy[k]);
	}
      }
    }
    fprintf(popenstream,"End");
  }
  // fermeture du stream pour popen
  fclose(popenstream);
  delete [] commandline; 
  return valsortie;
}

class  PopenMeditMesh : public OneOperator { public:  
     PopenMeditMesh() : OneOperator(atype<long>(), atype<string *>(), atype<pmesh>()) {}
  
  E_F0 * code(const basicAC_F0 & args) const 
  {
	cout << " PopenMeditMesh " << endl; 
	//cout << "args: " << args << endl;   
	return  new  PopenMeditMesh_Op( args,t[0]->CastTo(args[0]),t[1]->CastTo(args[1]) ); 
  }
};

// solution 

class  PopenMeditMesh2ALL : public OneOperator { public:  
  PopenMeditMesh2ALL() : OneOperator(atype<long>(), atype<string*>(), atype<pmesh>()) {}
  
  class PopenMeditMesh2ALL_Op : public E_F0mps 
  {
  public:
    Expression eTh;
    Expression filename;
    Expression xx, yy; //expression  vector 
    Expression tsxx, tsyx, tsyy;  //expression sym tensor  
    static const int n_name_param = 5;  
    static basicAC_F0::name_and_type name_param[] ;
    Expression nargs[n_name_param];
    long    arg(int i,Stack stack,long a ) const{ return nargs[i] ? GetAny< long >( (*nargs[i])(stack) ): a;}
    
  public:
    PopenMeditMesh2ALL_Op(const basicAC_F0 &  args,Expression ffname,Expression tth) 
      : eTh(tth),filename(ffname),xx(0),yy(0),tsxx(0), tsyx(0), tsyy(0) 
    {
      //cout << "construction meditmesh" << args << endl;
      args.SetNameParam(n_name_param,name_param,nargs); 
      const E_Array * a1=0 ;
      if(nargs[0])  a1  = dynamic_cast<const E_Array *>(nargs[2]);
      const E_Array * a2=0 ;
      if(nargs[0])  a2  = dynamic_cast<const E_Array *>(nargs[3]);
      
      if(a1) {
	if(a1->size() !=2) 
	  CompileError("meditmesh (Th,vector=[fx,fy,fz],) ");
	xx=to< pair< FEbase<double,v_fes> *,int > >( (*a1)[0]);
	yy=to< pair< FEbase<double,v_fes> *,int > >( (*a1)[1]);
      }    
      
      if(a2) {
	if(a2->size() !=3) 
	  CompileError("meditmesh (Th,symtensor=[Txx,Tyx,Tyy,Tzx,Tzy,Tzz],) ");
	tsxx=to< pair< FEbase<double,v_fes> *,int > >( (*a2)[0]);
	tsyx=to< pair< FEbase<double,v_fes> *,int > >( (*a2)[1]);
	tsyy=to< pair< FEbase<double,v_fes> *,int > >( (*a2)[2]);
      }   
      
      //cout << "fin construction meditmesh" << args << endl;
    } 
    
    AnyType operator()(Stack stack)  const ;
  };


  E_F0 * code(const basicAC_F0 & args) const 
  {
    cout << " PopenMeditMeshVh const" << endl; 
    //cout << "args: " << args << endl;   
    return  new  PopenMeditMesh2ALL_Op(args,t[0]->CastTo(args[0]),t[1]->CastTo(args[1])); 
  }
};


basicAC_F0::name_and_type PopenMeditMesh2ALL::PopenMeditMesh2ALL_Op::name_param[]= {
  {  "solution", &typeid(long)},
  {  "scalar", &typeid(pair< FEbase<double,v_fes> *,int > )}, //pair< FEbase<double,v_fes> *,int >)},
  {  "vector", &typeid(E_Array)},
  {  "symtensor", &typeid(E_Array)},
  {  "order",&typeid(long)}
};


AnyType PopenMeditMesh2ALL::PopenMeditMesh2ALL_Op::operator()(Stack stack)  const 
{
  //typedef typename  v_fes::pfes pfes;
  //typedef typename  v_fes::FESpace FESpace;
  //typedef typename  FESpace::Mesh Mesh;
  //typedef typename  FESpace::FElement FElement;
  //typedef typename  Mesh::Element Element;
  //typedef typename  Mesh::Vertex Vertex;  
  //typedef typename  Mesh::RdHat RdHat;  
  //typedef typename  Mesh::Rd Rd;
  //cout << "anytype" << endl;

  MeshPoint *mp(MeshPointStack(stack)) , mps=*mp;
  Mesh * pTh= GetAny<Mesh *>((*eTh)(stack));
  string * ffname= GetAny<string *>( (*filename)(stack) );
  ffassert(pTh);
  Mesh &Th=*pTh;

  int nt = Th.nt;
  int nv = Th.nv;
  int nbe = Th.neb;
  long sol (arg(0,stack,0));
  long order (arg(4,stack,0));
  int ver  = GmfFloat;
  int dimp = 2;
  float fx,fy;

  
  //  string * ffname= GetAny<string *>( (*filename)(stack) );
  char * commandline = meditcmd(sol,*ffname);
  long valsortie=0;
  int typsol,nbsol=sol;
    
  FILE *popenstream= popen(commandline,"w");

  fprintf(popenstream,"MeshVersionFormatted\n");
  fprintf(popenstream,"%i\n",ver);
  fprintf(popenstream,"Dimension\n");
  fprintf(popenstream,"%i\n",dimp);
  fprintf(popenstream,"Vertices\n");
  fprintf(popenstream,"%i\n",nv);

  for (int k=0; k<nv; k++) {
    const  Mesh::Vertex & P = Th.vertices[k];
    fx=P.x; fy=P.y; 
    fprintf(popenstream,"%f %f %i\n",fx,fy,P.lab);
  }
  fprintf(popenstream,"Triangles\n");
  fprintf(popenstream,"%i\n",nt);
  for (int k=0; k<nt; k++) {
    const Mesh::Triangle & K(Th.t(k));
    int i0=Th.operator()(K[0])+1;
    int i1=Th.operator()(K[1])+1;
    int i2=Th.operator()(K[2])+1;
    int lab=K.lab;
    fprintf(popenstream,"%i %i %i %i \n",i0,i1,i2,lab);
  }
  fprintf(popenstream,"Edges\n");
  fprintf(popenstream,"%i\n",nbe);
  for (int k=0; k<nbe; k++) {
    const Mesh::BorderElement & K(Th.be(k));
    int i0=Th.operator()(K[0])+1;
    int i1=Th.operator()(K[1])+1;
    int lab=K.lab;
    fprintf(popenstream,"%i %i %i %i\n",i0,i1,lab);
  }
  fprintf(popenstream,"End");
  
  if(sol > 0){
    fprintf(popenstream,"MeshVersionFormatted %i\n",ver);
    fprintf(popenstream,"Dimension %i\n",dimp);
    if(order==1){
      fprintf(popenstream,"SolAtVertices\n");
      fprintf(popenstream,"%i\n",nv);
    }
    if(order==0){
      fprintf(popenstream,"SolAtTriangles\n");
      fprintf(popenstream,"%i\n",nt);
    }
    // detemination of solution
    // default scalaire // faire tableau pour plusieurs
    typsol = 1;
    nbsol = 1;
    if ( xx ){
      typsol  = 2;
    }
    if( tsxx ){
      typsol  = 3;
    }
    fprintf(popenstream,"%i %i\n",nbsol,typsol);
       
    if(typsol==1){
      pair< FEbase<double,v_fes> *  ,int> pp1=GetAny<pair< FEbase<double,v_fes> *,int> >( (*nargs[1])(stack) );
      FEbase<double,v_fes> & fe1( *pp1.first);
      int componante1=pp1.second;
      const  FESpace & Vh1(*fe1.Vh);

      if(order == 1){
	KN<double> solsca(nv);
	solsca=0.;
	KN<int> takemesh(nv);
	MeshPoint *mp3(MeshPointStack(stack)); 
	takemesh=0;
	for (int it=0;it<Th.nt;++it){
	  const Mesh::Triangle & K(Th.t(it));
	  for( int iv=0; iv<3; ++iv){
	    int i=Th(it,iv);  
	    
	    //if(takemesh[i]==0){
	    mp3->setP(&Th,it,iv);
	    R2 PHat(mp3->PHat.x,mp3->PHat.y);
	    const FElement KK1(Vh1[Th(K)]);
	    solsca[i] = solsca[i] + KK1(PHat,*fe1.x(),componante1,0);
	    
	    takemesh[i] = takemesh[i]+1;
	    //}
	  }
	}

	for(int i=0; i<nv; i++){
	  solsca[i] = solsca[i]/takemesh[i];
	}

	for(int k=0; k<nv; k++){
	  fprintf(popenstream,"%f\n",solsca[k]);
	}
      }

      if(order == 0){
	KN<double> solsca(Th.nt);
	MeshPoint *mp3(MeshPointStack(stack)); 
	R2 Cdg_hat = R2(1./3.,1./3.);  

	for (int it=0;it<Th.nt;++it){
	  const Mesh::Triangle & K(Th.t(it));
	  mp3->set( Th, K(Cdg_hat), Cdg_hat, K, K.lab);

	  const FElement KK1(Vh1[Th(K)]);
	  solsca[it] = KK1(Cdg_hat,*fe1.x(),componante1,0);	   
	}
	for(int k=0; k<Th.nt; k++){
	  fprintf(popenstream,"%f\n",solsca[k]);
	}
      }
    }
    else if(typsol==2){
      pair< FEbase<double,v_fes> *  ,int> pp1=GetAny<pair< FEbase<double,v_fes> *,int> >( (*xx)(stack) );
      FEbase<double,v_fes> & fe1( *pp1.first);
      int componante1=pp1.second;
      const  FESpace & Vh1(*fe1.Vh);
      
      pair< FEbase<double,v_fes> *  ,int> pp2=GetAny<pair< FEbase<double,v_fes> *,int> >( (*yy)(stack) );
      FEbase<double,v_fes> & fe2( *pp2.first);
      int componante2=pp2.second;
      const  FESpace & Vh2(*fe2.Vh);
     
      if(order==1){
	KN<double> vxx(nv),vyy(nv);
	vxx=0.;
	vyy=0.;
	KN<int> takemesh(nv);
	MeshPoint *mp3(MeshPointStack(stack)); 
	
	takemesh=0;
	for (int it=0;it<Th.nt;++it){
	  const Mesh::Triangle & K(Th.t(it));
	  for( int iv=0; iv<3; ++iv){
	    int i=Th(it,iv);  
	    
	    //if(takemesh[i]==0){
	    mp3->setP(&Th,it,iv);
	    R2 PHat(mp3->PHat.x,mp3->PHat.y);

	    const FElement KK1(Vh1[Th(K)]);
	    const FElement KK2(Vh2[Th(K)]);
	    
	    vxx[i] = vxx[i] + KK1(PHat,*fe1.x(),componante1,0);
	    vyy[i] = vyy[i] + KK2(PHat,*fe2.x(),componante2,0);
	    takemesh[i] = takemesh[i]+1;
	    //}
	  }
	}

	for(int i=0;i<nv;i++){
	  vxx[i] = vxx[i]/takemesh[i];
	  vyy[i] = vyy[i]/takemesh[i];
	}

	for(int k=0; k<nv; k++){
	  fprintf(popenstream,"%f %f\n",vxx[k],vyy[k]);
	}
      }
      if(order==0){
	KN<double> vxx(nt),vyy(nt);
	MeshPoint *mp3(MeshPointStack(stack)); 
	
	R2 Cdg_hat = R2(1./3.,1./3.);  
	for (int it=0;it<Th.nt;++it){
	  const Mesh::Triangle & K(Th.t(it));
	  mp3->set( Th, K(Cdg_hat), Cdg_hat, K, K.lab);
	  
	  const FElement KK1(Vh1[Th(K)]);
	  const FElement KK2(Vh2[Th(K)]);
	   
	  vxx[it] = KK1(Cdg_hat,*fe1.x(),componante1,0);
	  vyy[it] = KK2(Cdg_hat,*fe2.x(),componante2,0);
	}
	for(int k=0; k<Th.nt; k++){
	  fprintf(popenstream,"%f %f\n",vxx[k],vyy[k]);
	}
      }      
    }
    else if(typsol==3){
      pair< FEbase<double,v_fes> *  ,int> pp1=GetAny<pair< FEbase<double,v_fes> *,int> >( (*tsxx)(stack) );
      FEbase<double,v_fes> & fe1( *pp1.first);
      int componante1=pp1.second;
      const  FESpace & Vh1(*fe1.Vh);
      
      pair< FEbase<double,v_fes> *  ,int> pp2=GetAny<pair< FEbase<double,v_fes> *,int> >( (*tsyx)(stack) );
      FEbase<double,v_fes> & fe2( *pp2.first);
      int componante2=pp2.second;
      const  FESpace & Vh2(*fe2.Vh);
      
      pair< FEbase<double,v_fes> *  ,int> pp3=GetAny<pair< FEbase<double,v_fes> *,int> >( (*tsyy)(stack) );
      FEbase<double,v_fes> & fe3( *pp3.first);
      int componante3=pp3.second;
      const  FESpace & Vh3(*fe3.Vh);

    
      if(order==1){
	KN<double> vxx(nv),vyx(nv),vyy(nv);
	vxx=0.;
	vyx=0.;
	vyy=0.;
	KN<int> takemesh(nv);
	MeshPoint *mp3(MeshPointStack(stack)); 
	
	takemesh=0;
	for (int it=0;it<Th.nt;++it){
	  const Mesh::Triangle & K(Th.t(it));
	  for( int iv=0; iv<3; ++iv){
	    int i=Th(it,iv);  
	    
	    //if(takemesh[i]==0){
	    mp3->setP(&Th,it,iv);
	    R2 PHat(mp3->PHat.x,mp3->PHat.y);
	    const FElement KK1(Vh1[Th(K)]);
	    const FElement KK2(Vh2[Th(K)]);
	    const FElement KK3(Vh3[Th(K)]);
	     
	    vxx[i] = vxx[i] + KK1(PHat,*fe1.x(),componante1,0);
	    vyx[i] = vyx[i] + KK2(PHat,*fe2.x(),componante2,0); 
	    vyy[i] = vyy[i] + KK3(PHat,*fe3.x(),componante3,0); 
	    
	    takemesh[i] = takemesh[i]+1;
	    //}
	  }
	}
	
	for(int i=0;i<nv;i++){
	  vxx[i] = vxx[i]/takemesh[i];
	  vyx[i] = vyx[i]/takemesh[i];
	  vyy[i] = vyy[i]/takemesh[i];
	}

	for(int k=0; k<nv; k++){
	  fprintf(popenstream,"%f %f %f\n",vxx[k],vyx[k],vyy[k]);
	}
      }
      if(order==0){
	KN<double> vxx(nt),vyx(nt),vyy(nt);
	MeshPoint *mp3(MeshPointStack(stack)); 

	R2 Cdg_hat = R2(1./3.,1./3.);  
	for (int it=0;it<Th.nt;it++){
	  const Mesh::Triangle & K(Th.t(it));
	  mp3->set( Th, K(Cdg_hat), Cdg_hat, K, K.lab);
	      
	  const FElement KK1(Vh1[Th(K)]);
	  const FElement KK2(Vh2[Th(K)]);
	  const FElement KK3(Vh3[Th(K)]);
	 
	  vxx[it] = KK1(Cdg_hat,*fe1.x(),componante1,0);
	  vyx[it] = KK2(Cdg_hat,*fe2.x(),componante2,0); 
	  vyy[it] = KK3(Cdg_hat,*fe3.x(),componante3,0); 
	  
	}
	for(int k=0; k<Th.nt; k++){
	  fprintf(popenstream,"%f %f %f\n",vxx[k],vyx[k],vyy[k]);
	}
      }
    }
    fprintf(popenstream,"End");
  }
  
  // fermeture du stream pour popen
  fclose(popenstream);
  delete [] commandline; 
  return valsortie;
}


// datasolMesh2

class datasolMesh2_Op : public E_F0mps 
{
public:
  Expression eTh;
  Expression filename;
  
  struct Expression2 {
    long what; // 1 scalar, 2 vector, 3 symtensor
    long nbfloat; // 1 scalar, 2 vector (3D), 3 symtensor(3D)
    Expression e[3];
    Expression2() {e[0]=0; e[1]=0; e[2]=0;  what=0; nbfloat=0;};
    Expression &operator[](int i){return e[i];}
    double eval(int i,Stack stack) const  { 
    if (e[i]) {
      return GetAny< double >( (*e[i])(stack) );
    }
    else 
      return 0;
    }
  };
  vector<Expression2> l; 
  static const int n_name_param =7;  
  static basicAC_F0::name_and_type name_param[];
  Expression nargs[n_name_param];
  long    arg(int i,Stack stack,long a ) const{ return nargs[i] ? GetAny< long >( (*nargs[i])(stack) ): a;}
 
public:
  datasolMesh2_Op(const basicAC_F0 &  args,  Expression fname, Expression tth) 
    : eTh(tth), filename(fname)
  {
    int nbofsol;
    int ddim=2;
    int stsize=3;
    // cout << "construction data medit solution" << args << endl;
    args.SetNameParam(n_name_param,name_param,nargs); 

    const E_Array * a0=0;
    const E_Array * a1=0;
    const E_Array * a2=0;

    if(nargs[0])  a0  = dynamic_cast<const E_Array *>(nargs[1]);
    if(nargs[2])  a1  = dynamic_cast<const E_Array *>(nargs[3]);
    if(nargs[4])  a2  = dynamic_cast<const E_Array *>(nargs[5]);
   
    nbofsol =0;
    cout << "l.size()= "<< nbofsol << endl;
    if(a0){
      nbofsol = nbofsol+a0->size();
    }
    if(a1){
      if(a1->size()%ddim !=0)
	CompileError("the number of element of parameters vector is not correct");
       nbofsol = nbofsol+a1->size()/ddim;
    }    

    if(a2){
      if(a2->size()%stsize !=0)
	CompileError("the number of element of parameters symtensor is not correct");
       nbofsol = nbofsol+a2->size()/stsize;
    }    
    cout << "l.size()= "<< nbofsol << endl;
    
    if(nbofsol == 0)
      CompileError("there is no solution");
    
    //vector<Expression2> l1[nbofsol]; 
    l= vector<Expression2> (nbofsol);
    cout << "l.size()= " << l.size() << endl;
    assert(nbofsol==l.size());
    size_t i=0;

    if(a0){
      for(size_t ii=0; ii <a0->size(); ii++){
	l[i].what=1;
	l[i].nbfloat=1;
	l[i][0]=to<double>( (*a0)[ii]);
	i++;
      }
    }
    if(a1){
      for(size_t ii=0; ii <a1->size()%ddim; ii++){
	l[i].what=2;
	l[i].nbfloat=ddim;
	for(int j=0; j<ddim; j++){
	  l[i][j] = to<double>( (*a1)[ii*ddim+j]);
	}
	i++;
      }
    }
    if(a2){
      for(size_t ii=0; ii <a2->size()%stsize; ii++){
	l[i].what=3;
	l[i].nbfloat=stsize;
	for(int j=0; j<stsize; j++){
	  l[i][j] = to<double>( (*a2)[ii*stsize+j]);
	}
	i++;
      }
    }
    
    // cout << "fin construction meditmesh" << args << endl;
  }
  
  AnyType operator()(Stack stack)  const ;
};

basicAC_F0::name_and_type datasolMesh2_Op::name_param[]= {
  {  "nbscalar",&typeid(long)},
  {  "scalar", &typeid(E_Array)},
  {  "nbvector",&typeid(long)},
  {  "vector", &typeid(E_Array)},
  {  "nbsymtensor",&typeid(long)},
  {  "symtensor", &typeid(E_Array)},
  {  "order",&typeid(long)}
};

AnyType datasolMesh2_Op::operator()(Stack stack)  const 
{ 
  MeshPoint *mp(MeshPointStack(stack)) , mps=*mp;
  Mesh * pTh= GetAny<Mesh *>((*eTh)(stack));
  string * ffname= GetAny<string *>( (*filename)(stack) );
  ffassert(pTh);
  Mesh &Th=*pTh;

  int nt = Th.nt;
  int nv = Th.nv;
  cout << "l.size()="<< l.size() << endl;
  int nbtype=l.size();
  int nbsol;
  int solnbfloat;
  int TypTab[l.size()];
 
  int resultorder= arg(6, stack, 1);
  long longdefault;

  int ver = GmfFloat, outm;
  // determination de TypTab
  solnbfloat=0;
  for (size_t i=0;i<l.size();i++){
    TypTab[i]=l[i].what;
    solnbfloat=solnbfloat+l[i].nbfloat;
  }
  float *OutSolTab = new float[solnbfloat];

  // determination de OutSolTab

  if ( !(outm = GmfOpenMesh(ffname->c_str(),GmfWrite,ver,2)) ) {
    cerr <<"  -- Mesh3::Save  UNABLE TO OPEN  :"<< filename << endl;
    exit(1);
  }
 
  if(resultorder==0){
    // ordre 0
    nbsol = nt;  
    
    KN<double> valsol(solnbfloat*nbsol);
 
    MeshPoint *mp3(MeshPointStack(stack)); 
    R2 Cdg_hat = R2(1./3.,1./3.);  

    for (int it=0;it<nt;it++){
      int h=0;
      const Mesh::Triangle  & K(Th.t(it));
      mp3->set( Th, K(Cdg_hat), Cdg_hat, K, K.lab);
	  
      for(size_t i=0;i<l.size();i++){
	for(size_t j=0;j<l[i].nbfloat;j++){
	  valsol[it*solnbfloat+h] = l[i].eval(j,stack);
	  h=h+1;
	}
      } 
      assert(solnbfloat==h);
    }


    GmfSetKwd(outm,GmfSolAtTriangles, nbsol, nbtype, TypTab);
    for (int k=0; k<nbsol; k++){
      for (int i=0; i<solnbfloat ;i++){
	OutSolTab[i] =  valsol(k*solnbfloat+i);
      }
      GmfSetLin(outm, GmfSolAtTriangles, OutSolTab); 
    }
  }
  if(resultorder==1){
    // ordre 1
    nbsol = nv;  

    KN<double> valsol(solnbfloat*nbsol);
    valsol=0.;
    KN<int> takemesh(nbsol);
    MeshPoint *mp3(MeshPointStack(stack)); 
    //R2 Cdg_hat = R2(1./3.,1./3.);
    takemesh=0;
    for (int it=0;it<nt;it++){
      for(int iv=0;iv<3;iv++){
	int i=Th(it,iv);

	//if(takemesh[i]==0){
	mp3->setP(&Th,it,iv);
	int h=0;
	
	for(size_t ii=0;ii<l.size();ii++){
	  for(size_t j=0;j<l[ii].nbfloat;j++){
	    //cout << "ii=" << ii << " j=" << j<< endl;
	    valsol[i*solnbfloat+h] = valsol[i*solnbfloat+h] + l[ii].eval(j,stack);
	    h=h+1;
	  }
	} 
	assert(solnbfloat==h);
	takemesh[i] = takemesh[i]+1;
	//}
      }
    }
    for(int i=0; i<nv; i++){
      for(int h=0; h<solnbfloat; h++){
	valsol[i*solnbfloat+h] = valsol[i*solnbfloat+h]/takemesh[i]; 
      }
    }

    GmfSetKwd(outm,GmfSolAtVertices, nbsol, nbtype, TypTab);
    for (int k=0; k<nbsol; k++){
      for (int i=0; i<solnbfloat ;i++){
	OutSolTab[i] =  valsol(k*solnbfloat+i);
      }
      GmfSetLin(outm, GmfSolAtVertices, OutSolTab);
    }
  }
  GmfCloseMesh(outm);
  return longdefault;
}
  
class  datasolMesh2 : public OneOperator { public:  
  datasolMesh2() : OneOperator(atype<long>(),atype<string *>(),atype<pmesh>()) {}
  
  E_F0 * code(const basicAC_F0 & args) const 
  {
	cout << " PopenMeditMesh3 je suis dans code(const basicAC_F0 & args) const" << endl; 
	//cout << "args: " << args << endl;   
	return  new  datasolMesh2_Op(args,t[0]->CastTo(args[0]),t[1]->CastTo(args[1]) ); 
  }
};

// datasolMesh2VhALL

//template<class v_fes>  
class  datasolMesh2VhALL : public OneOperator { public:  
  datasolMesh2VhALL() : OneOperator(atype<long>(),atype<string *>(),atype<pmesh>()) {}
  
   class datasolMesh2VhALL_Op : public E_F0mps
   {
   public:
     Expression eTh;
     Expression filename;
     struct Expression2 { 
       long what; // 1 scalar, 2 vector, 3 symtensor
       long nbfloat; // 1 scalar, 2 vector (2D), 3 symtensor(2D)
       Expression e[3];
       Expression2() {e[0]=0; e[1]=0; e[2]=0; what=0; nbfloat=0;};
       Expression &operator[](int i){return e[i];}
       double eval(int &i,const Mesh &Th, const Mesh::Triangle &K, const R2 &PHat,Stack stack) const;  
       
     };
     vector<Expression2> l; 
     static const int n_name_param =7;  
     static basicAC_F0::name_and_type name_param[];
     Expression nargs[n_name_param];
     long    arg(int i,Stack stack,long a ) const{ return nargs[i] ? GetAny< long >( (*nargs[i])(stack) ): a;}
     
   public:
     datasolMesh2VhALL_Op(const basicAC_F0 &  args,  Expression fname, Expression tth) 
       : eTh(tth), filename(fname)
     {
        int nbofsol;
	int ddim=2;
	int stsize=3;
	//cout << "construction data medit solution" << args << endl;
       args.SetNameParam(n_name_param,name_param,nargs);    
       const E_Array * a0=0;
       const E_Array * a1=0;
       const E_Array * a2=0;
       if(nargs[0])  a0  = dynamic_cast<const E_Array *>(nargs[1]);
       if(nargs[2])  a1  = dynamic_cast<const E_Array *>(nargs[3]);
       if(nargs[4])  a2  = dynamic_cast<const E_Array *>(nargs[5]);

        nbofsol =0;
	cout << "l.size()= "<< nbofsol << endl;
	if(a0){
	  nbofsol = nbofsol+a0->size();
	}
	if(a1){
	  if(a1->size()%ddim !=0)
	    CompileError("the number of element of parameters vector is not correct");
	  nbofsol = nbofsol+a1->size()/ddim;
	}    
	
	if(a2){
	  if(a2->size()%stsize !=0)
	    CompileError("the number of element of parameters symtensor is not correct");
	  nbofsol = nbofsol+a2->size()/stsize;
	}    
	
	if(nbofsol == 0)
	  CompileError("there is no solution");
	
	//vector<Expression2> l1[nbofsol]; 
	l= vector<Expression2> (nbofsol);
	assert(nbofsol==l.size());

	size_t i=0;
	
	if(a0){
	  for(size_t ii=0; ii <a0->size(); ii++){
	    l[i].what=1;
	    l[i].nbfloat=1;
	    l[i][0] = to< pair< FEbase<double,v_fes> *,int > >( (*a0)[ii]);
	    i++;
	  }
	}
	if(a1){
	  for(size_t ii=0; ii <(a1->size()/ddim); ii++){
	    l[i].what=2;
	    l[i].nbfloat=ddim;
	    for(int j=0; j<ddim; j++){
	      l[i][j] = to< pair< FEbase<double,v_fes> *,int > >( (*a1)[ii*ddim+j]);
	    }
	    i++;
	  }
	}
	if(a2){
	  for(size_t ii=0; ii <(a2->size()/stsize); ii++){
	    l[i].what=3;
	    l[i].nbfloat=stsize;
	    for(int j=0; j<stsize; j++){
	      l[i][j] = to< pair< FEbase<double,v_fes> *,int > >( (*a2)[ii*stsize+j]);
	    }
	    i++;
	  }
	}

	// cout << "fin construction meditmesh" << args << endl;
     }
     
     AnyType operator()(Stack stack)  const ;
   };
  
  E_F0 * code(const basicAC_F0 & args) const 
  {
	cout << " PopenMeditMesh2 " << endl; 
	//cout << "args: " << args << endl;   
	return  new  datasolMesh2VhALL_Op(args,t[0]->CastTo(args[0]),t[1]->CastTo(args[1]) ); 
  }

};

//template<class v_fes>
double datasolMesh2VhALL::datasolMesh2VhALL_Op::Expression2::eval(int &i,const Mesh &Th, const Mesh::Triangle &K, const R2 &PHat, Stack stack) const{
  //typedef typename  v_fes::pfes pfes;
  //typedef typename  v_fes::FESpace FESpace;
  //typedef typename  FESpace::Mesh Mesh;
  //typedef typename  FESpace::FElement FElement;
  //typedef typename  Mesh::Element Element;
  //typedef typename  Mesh::Vertex Vertex;  
  //ypedef typename  Mesh::RdHat RdHat;  
  //typedef typename  Mesh::Rd Rd;  
  if (e[i]) {
    pair< FEbase<double,v_fes> *  ,int> pp1=GetAny< pair< FEbase<double,v_fes> *,int> >( (*e[i])(stack) );
    FEbase<double,v_fes> & fe1( *pp1.first);
    int componante1=pp1.second;
    const  FESpace & Vh1(*fe1.Vh);  
    const FElement KK1(Vh1[Th(K)]);
    
    return KK1(PHat,*fe1.x(),componante1,0);
  }
  else 
    return 0;
}

//template<class v_fes>
basicAC_F0::name_and_type datasolMesh2VhALL::datasolMesh2VhALL_Op::name_param[]= {
  {  "nbscalar",&typeid(long)},
  {  "scalar", &typeid(E_Array)},
  {  "nbvector",&typeid(long)},
  {  "vector",&typeid(E_Array)},
  {  "nbsymtensor",&typeid(long)},
  {  "symtensor", &typeid(E_Array)},
  {  "order",&typeid(long)},
};

//template<class v_fes>
AnyType  datasolMesh2VhALL::datasolMesh2VhALL_Op::operator()(Stack stack)  const 
{ 
  //typedef typename  v_fes::pfes pfes;
  //typedef typename  v_fes::FESpace FESpace;
  //typedef typename  FESpace::Mesh Mesh;
  //typedef typename  FESpace::FElement FElement;
  //typedef typename  Mesh::Element Element;
  //typedef typename  Mesh::Vertex Vertex;  
  //typedef typename  Mesh::RdHat RdHat;  
  //typedef typename  Mesh::Rd Rd;  


  MeshPoint *mp(MeshPointStack(stack)) , mps=*mp;
  Mesh * pTh= GetAny<Mesh *>((*eTh)(stack));
  string * ffname= GetAny<string *>( (*filename)(stack) );
  ffassert(pTh);
  Mesh &Th=*pTh;

  int nt = Th.nt;
  int nv = Th.nv;
 
  int nbsol;
  int solnbfloat;
  int resultorder= arg(6, stack, 1);
  long longdefault;

  int ver = GmfFloat, outm;
  int nbtype = l.size();
  int TypTab[nbtype];

  solnbfloat=0;
  for (size_t i=0;i<l.size();i++){
    TypTab[i]=l[i].what;
    solnbfloat=solnbfloat+l[i].nbfloat;
    cout << "i="<< i << " TypTab[i]=" << TypTab[i] << " solnbfloat=" << solnbfloat << endl; 
  }
  float *OutSolTab = new float[solnbfloat];

 
  // determination de OutSolTab
  if ( !(outm = GmfOpenMesh(ffname->c_str(),GmfWrite,ver,2)) ) {
    cerr <<"  -- Mesh::Save  UNABLE TO OPEN  :"<< filename << endl;
    exit(1);
  }
 
  if(resultorder==0){
    // ordre 0
    nbsol = nt;  
  
    KN<double> valsol(solnbfloat*nbsol);
    MeshPoint & mp3 = *MeshPointStack(stack); 
    R2 Cdg_hat = R2(1./3.,1./3.);  
    
    for (int it=0;it<nt;it++){
      int h=0;
      Mesh::Triangle & K(Th.t(it));
      mp3.set( Th, K(Cdg_hat), Cdg_hat, K, K.lab);
      
      for(int i=0;i<l.size();i++){
	for(int j=0;j<l[i].nbfloat;j++){
	  valsol[it*solnbfloat+h] = l[i].eval(j,Th,K,Cdg_hat,stack);
	  h=h+1;
	}
      } 
      assert(solnbfloat==h);
    }
  
    GmfSetKwd(outm,GmfSolAtTriangles, nbsol, nbtype, TypTab);
    for (int k=0; k<nbsol; k++){
      for (int i=0; i<solnbfloat ;i++){
	OutSolTab[i] =  valsol(k*solnbfloat+i);
      }
      GmfSetLin(outm, GmfSolAtTriangles, OutSolTab); 
    }
  }
  if(resultorder==1){
    // ordre 1
    nbsol = nv;  

    KN<double> valsol(solnbfloat*nbsol);
    valsol=0.;
    KN<int> takemesh(nbsol);
    MeshPoint *mp3(MeshPointStack(stack)); 
    takemesh=0;
    for (int it=0;it<nt;it++){
      const Mesh::Triangle & K(Th.t(it));
      for(int iv=0;iv<3;iv++){
	int i=Th(it,iv);

	if(takemesh[i]==0){
	    mp3->setP(&Th,it,iv);
	    const R2 PHat(mp3->PHat.x,mp3->PHat.y);
	    int h=0;

	    for(int ii=0;ii<l.size();ii++){
	      for(int j=0;j<l[ii].nbfloat;j++){
		//cout << "ii=" << ii << " j=" << j<< endl;
		valsol[i*solnbfloat+h] = l[ii].eval(j,Th,K,PHat,stack);
		h=h+1;
	      }
	    } 
	    assert(solnbfloat==h);
	    takemesh[i] = takemesh[i]+1;
	}
	/*
	  mp3->setP(&Th,it,iv);
	  int h=0;
	  
	  for(size_t ii=0;ii<l.size();ii++){
	  for(size_t j=0;j<l[ii].nbfloat;j++){
	  //cout << "ii=" << ii << " j=" << j<< endl;
	  valsol[i*solnbfloat+h] =  valsol[i*solnbfloat+h]+l[ii].eval(j,Th,K,mp3->PHat,stack);
	  h=h+1;
	  }
	  } 
	  assert(solnbfloat==h);
	  takemesh[i] = takemesh[i]+1;
	*/

      }
    }

    /*
      for(int i=0; i<nv ; i++){
      for(int h=0; h<solnbfloat; h++){
      valsol[i*solnbfloat+h]  = valsol[i*solnbfloat+h]/takemesh[i];
      }
      }
    */

    GmfSetKwd(outm,GmfSolAtVertices, nbsol, nbtype, TypTab);
    for (int k=0; k<nbsol; k++){
      for (int i=0; i<solnbfloat ;i++){
	OutSolTab[i] =  valsol(k*solnbfloat+i);
      }
      GmfSetLin(outm, GmfSolAtVertices, OutSolTab);
    }
  }

  GmfCloseMesh(outm);
  return longdefault;
}


// maillage 3D

class PopenMeditMesh3_Op : public E_F0mps 
{
public:
  Expression eTh;
  Expression filename;
  Expression xx, yy, zz; //expression  vector 
  Expression tsxx, tsyx, tsyy, tszx, tszy, tszz;  //expression sym tensor  
  static const int n_name_param =4;  
  static basicAC_F0::name_and_type name_param[] ;
  Expression nargs[n_name_param];
  long    arg(int i,Stack stack,long a ) const{ return nargs[i] ? GetAny< long >( (*nargs[i])(stack) ): a;}
 
public:
  PopenMeditMesh3_Op(const basicAC_F0 &  args,Expression ffname, Expression tth) 
    : eTh(tth),filename(ffname),xx(0),yy(0),zz(0),tsxx(0), tsyx(0), tsyy(0), tszx(0), tszy(0), tszz(0)
  {
    // cout << "construction meditmesh" << args << endl;
    args.SetNameParam(n_name_param,name_param,nargs); 
    const E_Array * a1=0 ;
    if(nargs[0])  a1  = dynamic_cast<const E_Array *>(nargs[2]);
    const E_Array * a2=0 ;
    if(nargs[0])  a2  = dynamic_cast<const E_Array *>(nargs[3]);
   
    if(a1) {
      if(a1->size() !=3) 
      CompileError("meditmesh (Th,vector=[fx,fy,fz],) ");
      xx=to<double>( (*a1)[0]);
      yy=to<double>( (*a1)[1]);
      zz=to<double>( (*a1)[2]);
    }    
    
    if(a2) {
      if(a2->size() !=6) 
      CompileError("meditmesh (Th,symtensor=[Txx,Tyx,Tyy,Tzx,Tzy,Tzz],) ");
      tsxx=to<double>( (*a2)[0]);
      tsyx=to<double>( (*a2)[1]);
      tsyy=to<double>( (*a2)[2]);
      tszx=to<double>( (*a2)[3]);
      tszy=to<double>( (*a2)[4]);
      tszz=to<double>( (*a2)[5]);
    }   

    // cout << "fin construction meditmesh" << args << endl;
  } 
  
  AnyType operator()(Stack stack)  const ;
};

basicAC_F0::name_and_type PopenMeditMesh3_Op::name_param[]= {
  {  "solution", &typeid(long)},
  {  "scalar", &typeid(double)},
  {  "vector", &typeid(E_Array)},
  {  "symtensor", &typeid(E_Array)}
};

AnyType PopenMeditMesh3_Op::operator()(Stack stack)  const 
{
  cout << "anytype" << endl;
  MeshPoint *mp(MeshPointStack(stack)) , mps=*mp;
  Mesh3 * pTh= GetAny<Mesh3 *>((*eTh)(stack));
  ffassert(pTh);
  Mesh3 &Th=*pTh;
  
  int nt = Th.nt;
  int nv = Th.nv;
  int nbe = Th.nbe;
  long sol (arg(0,stack,0));

  int ver = GmfFloat;
  int dimp =3;
  float fx,fy,fz;

  // char *medit_char= "medit-2.3-win.exeCC -popen 1";  // depend de l endroit ou se trouve medit
  //char *medit_charsol= "medit-2.3-win.exeCC -popen -addsol 1";
  
 //  char commandline[500];
//   long valsortie=0;
//   int typsol,nbsol;

//   string * ffname= GetAny<string *>( (*filename)(stack) );
//   size_t size_ffname = ffname->size()+1;
//   char* medit_name =new char[size_ffname];
//   strncpy( medit_name, ffname->c_str(), size_ffname); 

//   if(sol==0)
//     snprintf(commandline,500,"%s %s",medit_char,medit_name);
//   else
//     snprintf(commandline,500,"%s %s",medit_charsol,medit_name);
    
//   printf("version de medit %s :::  Vertices%i\n",commandline,nv);
  
//   FILE *popenstream= popen(commandline,"w");

  string * ffname= GetAny<string *>( (*filename)(stack) );
  char * commandline = meditcmd(sol,*ffname);
  long valsortie=0;
  int typsol,nbsol=sol;
    
  FILE *popenstream= popen(commandline,"w");

  fprintf(popenstream,"MeshVersionFormatted\n");
  fprintf(popenstream,"%i\n",ver);
  fprintf(popenstream,"Dimension\n");
  fprintf(popenstream,"%i\n",dimp);
  fprintf(popenstream,"Vertices\n");
  fprintf(popenstream,"%i\n",nv);

  for (int k=0; k<nv; k++) {
    const  Vertex3 & P = Th.vertices[k];
    fx=P.x; fy=P.y; fz=P.z;
    fprintf(popenstream,"%f %f %f %i\n",fx,fy,fz,P.lab);
  }
  fprintf(popenstream,"Tetrahedra\n");
  fprintf(popenstream,"%i\n",nt);
  for (int k=0; k<nt; k++) {
    const Tet & K(Th.elements[k]);
    int i0=Th.operator()(K[0])+1;
    int i1=Th.operator()(K[1])+1;
    int i2=Th.operator()(K[2])+1;
    int i3=Th.operator()(K[3])+1;
    int lab=K.lab;
    fprintf(popenstream,"%i %i %i %i %i\n",i0,i1,i2,i3,lab);
  }
  fprintf(popenstream,"Triangles\n");
  fprintf(popenstream,"%i\n",nbe);
  for (int k=0; k<nbe; k++) {
    const Triangle3 & K(Th.be(k));
    int i0=Th.operator()(K[0])+1;
    int i1=Th.operator()(K[1])+1;
    int i2=Th.operator()(K[2])+1;
    int lab=K.lab;
    fprintf(popenstream,"%i %i %i %i\n",i0,i1,i2,lab);
  }
  fprintf(popenstream,"End");
  
  if(sol > 0){
    fprintf(popenstream,"MeshVersionFormatted %i\n",ver);
    fprintf(popenstream,"Dimension %i\n",dimp);
    fprintf(popenstream,"SolAtVertices\n");
    fprintf(popenstream,"%i\n",nv);
    // detemination of solution
    // default scalaire // faire tableau pour plusieurs
    typsol = 1;
    nbsol = 1;
    if ( xx ){
      typsol  = 2;
    }
    if( tsxx ){
      typsol  = 3;
    }
    fprintf(popenstream,"%i %i\n",nbsol,typsol);
       
    if(typsol==1){
      KN<double> solsca(nv);
      KN<int> takemesh(nv);
      MeshPoint *mp3(MeshPointStack(stack)); 
      
      takemesh=0;
      for (int it=0;it<Th.nt;++it){
	for( int iv=0; iv<4; ++iv){
	  int i=Th(it,iv);  
	  
	  if(takemesh[i]==0){
	    mp3->setP(&Th,it,iv);
	    solsca[i] =  GetAny< double >( (*nargs[1])(stack) );
  
	    takemesh[i] = takemesh[i]+1;
	  }
	}
      }
      for(int k=0; k<nv; k++){
	fprintf(popenstream,"%f\n",solsca[k]);
      }
    }
    else if(typsol==2){
      KN<double> vxx(nv),vyy(nv),vzz(nv);
      KN<int> takemesh(nv);
      MeshPoint *mp3(MeshPointStack(stack)); 
      
      takemesh=0;
      for (int it=0;it<Th.nt;++it){
	for( int iv=0; iv<4; ++iv){
	  int i=Th(it,iv);  
	  
	  if(takemesh[i]==0){
	    mp3->setP(&Th,it,iv);
	    vxx[i] = GetAny< double >( (*xx)(stack) );
	    vyy[i] = GetAny< double >( (*yy)(stack) );
	    vzz[i] = GetAny< double >( (*zz)(stack) );

	    takemesh[i] = takemesh[i]+1;
	  }
	}
      }
      for(int k=0; k<nv; k++){
	fprintf(popenstream,"%f %f %f\n",vxx[k],vyy[k],vzz[k]);
      }
    }
    else if(typsol==3){
      KN<double> vxx(nv),vyx(nv),vyy(nv),vzx(nv),vzy(nv),vzz(nv);
      KN<int> takemesh(nv);
      MeshPoint *mp3(MeshPointStack(stack)); 
      
      takemesh=0;
      for (int it=0;it<Th.nt;++it){
	for( int iv=0; iv<4; ++iv){
	  int i=Th(it,iv);  
	  
	  if(takemesh[i]==0){
	    mp3->setP(&Th,it,iv);
	    vxx[i] = GetAny< double >( (*tsxx)(stack) );
	    vyx[i] = GetAny< double >( (*tsyx)(stack) );
	    vyy[i] = GetAny< double >( (*tsyy)(stack) );
	    vzx[i] = GetAny< double >( (*tszx)(stack) );
	    vzy[i] = GetAny< double >( (*tszy)(stack) );
	    vzz[i] = GetAny< double >( (*tszz)(stack) );
	    
	    takemesh[i] = takemesh[i]+1;
	  }
	}
      }
      for(int k=0; k<nv; k++){
	fprintf(popenstream,"%f %f %f %f %f %f\n",vxx[k],vyx[k],vyy[k],vzx[k],vzy[k],vzz[k]);
      }
    }
    fprintf(popenstream,"End");
  }
  
  // fermeture du stream pour popen
  fclose(popenstream);
  delete [] commandline; 
  return valsortie;
}

class  PopenMeditMesh3 : public OneOperator { public:  
     PopenMeditMesh3() : OneOperator(atype<long>(),atype<string*>(), atype<pmesh3>()) {}
  
  E_F0 * code(const basicAC_F0 & args) const 
  {
	cout << " PopenMeditMesh3 je suis dans code(const basicAC_F0 & args) const" << endl; 
	//cout << "args: " << args << endl;   
	return  new  PopenMeditMesh3_Op(args, t[0]->CastTo(args[0]), t[1]->CastTo(args[1]) ); 
  }
};



template<class v_fes>
class  PopenMeditMesh3ALL : public OneOperator { public:  
  PopenMeditMesh3ALL() : OneOperator(atype<long>(), atype<string*>(), atype<pmesh3>()) {}
  
  class PopenMeditMesh3ALL_Op : public E_F0mps 
  {
  public:
    Expression eTh;
    Expression filename;
    Expression xx, yy, zz; //expression  vector 
    Expression tsxx, tsyx, tsyy, tszx, tszy, tszz;  //expression sym tensor  
    static const int n_name_param = 5;  
    static basicAC_F0::name_and_type name_param[] ;
    Expression nargs[n_name_param];
    long    arg(int i,Stack stack,long a ) const{ return nargs[i] ? GetAny< long >( (*nargs[i])(stack) ): a;}
    
  public:
    PopenMeditMesh3ALL_Op(const basicAC_F0 &  args,Expression ffname,Expression tth) 
      : eTh(tth),filename(ffname),xx(0),yy(0),zz(0),tsxx(0), tsyx(0), tsyy(0), tszx(0), tszy(0), tszz(0)
    {
      // cout << "construction meditmesh" << args << endl;
      args.SetNameParam(n_name_param,name_param,nargs); 
      const E_Array * a1=0 ;
      if(nargs[0])  a1  = dynamic_cast<const E_Array *>(nargs[2]);
      const E_Array * a2=0 ;
      if(nargs[0])  a2  = dynamic_cast<const E_Array *>(nargs[3]);
      
      if(a1) {
	if(a1->size() !=3) 
	  CompileError("meditmesh (Th,vector=[fx,fy,fz],) ");
	xx=to< pair< FEbase<double,v_fes> *,int > >( (*a1)[0]);
	yy=to< pair< FEbase<double,v_fes> *,int > >( (*a1)[1]);
	zz=to< pair< FEbase<double,v_fes> *,int > >( (*a1)[2]);
      }    
      
      if(a2) {
	if(a2->size() !=6) 
	  CompileError("meditmesh (Th,symtensor=[Txx,Tyx,Tyy,Tzx,Tzy,Tzz],) ");
	tsxx=to< pair< FEbase<double,v_fes> *,int > >( (*a2)[0]);
	tsyx=to< pair< FEbase<double,v_fes> *,int > >( (*a2)[1]);
	tsyy=to< pair< FEbase<double,v_fes> *,int > >( (*a2)[2]);
	tszx=to< pair< FEbase<double,v_fes> *,int > >( (*a2)[3]);
	tszy=to< pair< FEbase<double,v_fes> *,int > >( (*a2)[4]);
	tszz=to< pair< FEbase<double,v_fes> *,int > >( (*a2)[5]);
      }   
      
      // cout << "fin construction meditmesh" << args << endl;
    } 
    
    AnyType operator()(Stack stack)  const ;
  };


  E_F0 * code(const basicAC_F0 & args) const 
  {
    cout << " PopenMeditMesh3 je suis dans code(const basicAC_F0 & args) const" << endl; 
    //cout << "args: " << args << endl;   
    return  new  PopenMeditMesh3ALL_Op(args,t[0]->CastTo(args[0]),t[1]->CastTo(args[1])); 
  }
};

template<class v_fes>
basicAC_F0::name_and_type PopenMeditMesh3ALL<v_fes>::PopenMeditMesh3ALL_Op::name_param[]= {
  {  "solution", &typeid(long)},
  {  "scalar", &typeid(pair< FEbase<double,v_fes> *,int >)},
  {  "vector", &typeid(E_Array)},
  {  "symtensor", &typeid(E_Array)},
  {  "order",&typeid(long)}
};

template<class v_fes>
AnyType PopenMeditMesh3ALL<v_fes>::PopenMeditMesh3ALL_Op::operator()(Stack stack)  const 
{
  typedef typename  v_fes::pfes pfes;
  typedef typename  v_fes::FESpace FESpace;
  typedef typename  FESpace::Mesh Mesh;
  typedef typename  FESpace::FElement FElement;
  //typedef typename  Mesh::Element Element;
  typedef typename  Mesh::Vertex Vertex;  
  typedef typename  Mesh::RdHat RdHat;  
  typedef typename  Mesh::Rd Rd;
  
  cout << "anytype" << endl;
  MeshPoint *mp(MeshPointStack(stack)) , mps=*mp;
  Mesh3 * pTh= GetAny<Mesh3 *>((*eTh)(stack));
  string * ffname= GetAny<string *>( (*filename)(stack) );
  ffassert(pTh);
  Mesh3 &Th=*pTh;

  int nt = Th.nt;
  int nv = Th.nv;
  int nbe = Th.nbe;
  long sol (arg(0,stack,0));
  long order (arg(4,stack,1));
  int ver = GmfFloat;
  int dimp =3;
  float fx,fy,fz;

  // char *medit_char= "medit-2.3-win.exeCC -popen 1";  // depend de l endroit ou se trouve medit
  //  char *medit_charsol= "medit-2.3-win.exeCC -popen -addsol 1";
  
//   char commandline[500];
//   long valsortie=0;
  
//   int typsol,nbsol;

//   size_t size_ffname = ffname->size()+1;
//   char* medit_name =new char[size_ffname];
//   strncpy( medit_name, ffname->c_str(), size_ffname); 

//   if(sol==0)
//     snprintf(commandline,500,"%s %s",medit_char,medit_name);
//   else
//     snprintf(commandline,500,"%s %s",medit_charsol,medit_name);
    
//   FILE *popenstream= popen(commandline,"w");


//  string * ffname= GetAny<string *>( (*filename)(stack) );
  char * commandline = meditcmd(sol,*ffname);
  long valsortie=0;
  int typsol,nbsol=sol;
    
  FILE *popenstream= popen(commandline,"w");

  fprintf(popenstream,"%s\n",ffname->c_str());
  fprintf(popenstream,"MeshVersionFormatted\n");
  fprintf(popenstream,"%i\n",ver);
  fprintf(popenstream,"Dimension\n");
  fprintf(popenstream,"%i\n",dimp);
  fprintf(popenstream,"Vertices\n");
  fprintf(popenstream,"%i\n",nv);

  for (int k=0; k<nv; k++) {
    const  Vertex3 & P = Th.vertices[k];
    fx=P.x; fy=P.y; fz=P.z;
    fprintf(popenstream,"%f %f %f %i\n",fx,fy,fz,P.lab);
  }
  fprintf(popenstream,"Tetrahedra\n");
  fprintf(popenstream,"%i\n",nt);
  for (int k=0; k<nt; k++) {
    const Tet & K(Th.elements[k]);
    int i0=Th.operator()(K[0])+1;
    int i1=Th.operator()(K[1])+1;
    int i2=Th.operator()(K[2])+1;
    int i3=Th.operator()(K[3])+1;
    int lab=K.lab;
    fprintf(popenstream,"%i %i %i %i %i\n",i0,i1,i2,i3,lab);
  }
  fprintf(popenstream,"Triangles\n");
  fprintf(popenstream,"%i\n",nbe);
  for (int k=0; k<nbe; k++) {
    const Triangle3 & K(Th.be(k));
    int i0=Th.operator()(K[0])+1;
    int i1=Th.operator()(K[1])+1;
    int i2=Th.operator()(K[2])+1;
    int lab=K.lab;
    fprintf(popenstream,"%i %i %i %i\n",i0,i1,i2,lab);
  }
  fprintf(popenstream,"End");
  
  if(sol > 0){
    fprintf(popenstream,"MeshVersionFormatted %i\n",ver);
    fprintf(popenstream,"Dimension %i\n",dimp);
    if(order==1){
      fprintf(popenstream,"SolAtVertices\n");
      fprintf(popenstream,"%i\n",nv);
    }
    if(order==0){
      fprintf(popenstream,"SolAtTetrahedra\n");
      fprintf(popenstream,"%i\n",nt);
    }
    // detemination of solution
    // default scalaire // faire tableau pour plusieurs
    typsol = 1;
    nbsol = 1;
    if ( xx ){
      typsol  = 2;
    }
    if( tsxx ){
      typsol  = 3;
    }
    fprintf(popenstream,"%i %i\n",nbsol,typsol);
       
    if(typsol==1){
      pair< FEbase<double,v_fes> *  ,int> pp1=GetAny<pair< FEbase<double,v_fes> *,int> >( (*xx)(stack) );
      FEbase<double,v_fes> & fe1( *pp1.first);
      int componante1=pp1.second;
      const  FESpace & Vh1(*fe1.Vh);

      if(order == 1){
	KN<double> solsca(nv);
	KN<int> takemesh(nv);
	MeshPoint *mp3(MeshPointStack(stack)); 
	takemesh=0;
	for (int it=0;it<Th.nt;++it){
	  const Tet & K(Th.elements[it]);
	  for( int iv=0; iv<4; ++iv){
	    int i=Th(it,iv);  
	    
	    if(takemesh[i]==0){
	      mp3->setP(&Th,it,iv);
	      const FElement KK1(Vh1[Th(K)]);
	      solsca[i] = KK1(mp3->PHat,*fe1.x(),componante1,0);
	      
	      takemesh[i] = takemesh[i]+1;
	    }
	  }
	}
	for(int k=0; k<nv; k++){
	  fprintf(popenstream,"%f\n",solsca[k]);
	}
      }

      if(order == 0){
	KN<double> solsca(Th.nt);
	MeshPoint *mp3(MeshPointStack(stack)); 
	R3 Cdg_hat = R3(1./4.,1./4.,1./4.);  

	for (int it=0;it<Th.nt;++it){
	  const Tet & K(Th.elements[it]);
	  mp3->set( Th, K(Cdg_hat), Cdg_hat, K, K.lab);

	  const FElement KK1(Vh1[Th(K)]);
	  solsca[it] = KK1(Cdg_hat,*fe1.x(),componante1,0);	   
	}
	for(int k=0; k<Th.nt; k++){
	  fprintf(popenstream,"%f\n",solsca[k]);
	}
      }
    }
    else if(typsol==2){
      pair< FEbase<double,v_fes> *  ,int> pp1=GetAny<pair< FEbase<double,v_fes> *,int> >( (*xx)(stack) );
      FEbase<double,v_fes> & fe1( *pp1.first);
      int componante1=pp1.second;
      const  FESpace & Vh1(*fe1.Vh);
      
      pair< FEbase<double,v_fes> *  ,int> pp2=GetAny<pair< FEbase<double,v_fes> *,int> >( (*yy)(stack) );
      FEbase<double,v_fes> & fe2( *pp2.first);
      int componante2=pp2.second;
      const  FESpace & Vh2(*fe2.Vh);
      
      pair< FEbase<double,v_fes> *  ,int> pp3=GetAny<pair< FEbase<double,v_fes> *,int> >( (*zz)(stack) );
      FEbase<double,v_fes> & fe3( *pp3.first);
      int componante3=pp3.second;
      const  FESpace & Vh3(*fe3.Vh);
     
      if(order==1){
	KN<double> vxx(nv),vyy(nv),vzz(nv);
	KN<int> takemesh(nv);
	MeshPoint *mp3(MeshPointStack(stack)); 
	
	takemesh=0;
	for (int it=0;it<Th.nt;++it){
	  const Tet & K(Th.elements[it]);
	  for( int iv=0; iv<4; ++iv){
	    int i=Th(it,iv);  
	    
	    if(takemesh[i]==0){
	      mp3->setP(&Th,it,iv);
	      const FElement KK1(Vh1[Th(K)]);
	      const FElement KK2(Vh2[Th(K)]);
	      const FElement KK3(Vh3[Th(K)]);
	      
	      vxx[i] = KK1(mp3->PHat,*fe1.x(),componante1,0);
	      vyy[i] = KK2(mp3->PHat,*fe2.x(),componante2,0);
	      vzz[i] = KK3(mp3->PHat,*fe3.x(),componante3,0);
	      takemesh[i] = takemesh[i]+1;
	    }
	  }
	}
	for(int k=0; k<nv; k++){
	  fprintf(popenstream,"%f %f %f\n",vxx[k],vyy[k],vzz[k]);
	}
      }
       if(order==0){
	KN<double> vxx(nt),vyy(nt),vzz(nt);
	MeshPoint *mp3(MeshPointStack(stack)); 
	
	R3 Cdg_hat = R3(1./4.,1./4.,1./4.);  
	for (int it=0;it<Th.nt;++it){
	  const Tet & K(Th.elements[it]);
	  mp3->set( Th, K(Cdg_hat), Cdg_hat, K, K.lab);
	      
	  const FElement KK1(Vh1[Th(K)]);
	  const FElement KK2(Vh2[Th(K)]);
	  const FElement KK3(Vh3[Th(K)]);
	 
	  vxx[it] = KK1(Cdg_hat,*fe1.x(),componante1,0);
	  vyy[it] = KK2(Cdg_hat,*fe2.x(),componante2,0);
	  vzz[it] = KK3(Cdg_hat,*fe3.x(),componante3,0);
	    
	}
	for(int k=0; k<Th.nt; k++){
	  fprintf(popenstream,"%f %f %f\n",vxx[k],vyy[k],vzz[k]);
	}
      }

      
    }
    else if(typsol==3){
      pair< FEbase<double,v_fes> *  ,int> pp1=GetAny<pair< FEbase<double,v_fes> *,int> >( (*tsxx)(stack) );
      FEbase<double,v_fes> & fe1( *pp1.first);
      int componante1=pp1.second;
      const  FESpace & Vh1(*fe1.Vh);
      
      pair< FEbase<double,v_fes> *  ,int> pp2=GetAny<pair< FEbase<double,v_fes> *,int> >( (*tsyx)(stack) );
      FEbase<double,v_fes> & fe2( *pp2.first);
      int componante2=pp2.second;
      const  FESpace & Vh2(*fe2.Vh);
      
      pair< FEbase<double,v_fes> *  ,int> pp3=GetAny<pair< FEbase<double,v_fes> *,int> >( (*tsyy)(stack) );
      FEbase<double,v_fes> & fe3( *pp3.first);
      int componante3=pp3.second;
      const  FESpace & Vh3(*fe3.Vh);

      pair< FEbase<double,v_fes> *  ,int> pp4=GetAny<pair< FEbase<double,v_fes> *,int> >( (*tszx)(stack) );
      FEbase<double,v_fes> & fe4( *pp4.first);
      int componante4=pp4.second;
      const  FESpace & Vh4(*fe4.Vh);
      
      pair< FEbase<double,v_fes> *  ,int> pp5=GetAny<pair< FEbase<double,v_fes> *,int> >( (*tszy)(stack) );
      FEbase<double,v_fes> & fe5( *pp5.first);
      int componante5=pp5.second;
      const  FESpace & Vh5(*fe5.Vh);
      
      pair< FEbase<double,v_fes> *  ,int> pp6=GetAny<pair< FEbase<double,v_fes> *,int> >( (*tszz)(stack) );
      FEbase<double,v_fes> & fe6( *pp6.first);
      int componante6=pp6.second;
      const  FESpace & Vh6(*fe6.Vh);

      if(order==1){
      KN<double> vxx(nv),vyx(nv),vyy(nv),vzx(nv),vzy(nv),vzz(nv);
      KN<int> takemesh(nv);
      MeshPoint *mp3(MeshPointStack(stack)); 
      
      takemesh=0;
      for (int it=0;it<Th.nt;++it){
	const Tet & K(Th.elements[it]);
	for( int iv=0; iv<4; ++iv){
	  int i=Th(it,iv);  
	  
	  if(takemesh[i]==0){
	    mp3->setP(&Th,it,iv);
	    const FElement KK1(Vh1[Th(K)]);
	    const FElement KK2(Vh2[Th(K)]);
	    const FElement KK3(Vh3[Th(K)]);
	    const FElement KK4(Vh4[Th(K)]);
	    const FElement KK5(Vh5[Th(K)]);
	    const FElement KK6(Vh6[Th(K)]);

	    vxx[i] = KK1(mp3->PHat,*fe1.x(),componante1,0);
	    vyx[i] = KK2(mp3->PHat,*fe2.x(),componante2,0); 
	    vyy[i] = KK3(mp3->PHat,*fe3.x(),componante3,0); 
	    vzx[i] = KK4(mp3->PHat,*fe4.x(),componante4,0);
	    vzy[i] = KK5(mp3->PHat,*fe5.x(),componante5,0);
	    vzz[i] = KK6(mp3->PHat,*fe6.x(),componante6,0);
	    
	    takemesh[i] = takemesh[i]+1;
	  }
	}
      }
      for(int k=0; k<nv; k++){
	fprintf(popenstream,"%f %f %f %f %f %f\n",vxx[k],vyx[k],vyy[k],vzx[k],vzy[k],vzz[k]);
      }
      }
      if(order==0){
	KN<double> vxx(nt),vyx(nt),vyy(nt),vzx(nt),vzy(nt),vzz(nt);
	MeshPoint *mp3(MeshPointStack(stack)); 

	R3 Cdg_hat = R3(1./4.,1./4.,1./4.);  
	for (int it=0;it<Th.nt;++it){
	  const Tet & K(Th.elements[it]);
	  mp3->set( Th, K(Cdg_hat), Cdg_hat, K, K.lab);
	      
	  const FElement KK1(Vh1[Th(K)]);
	  const FElement KK2(Vh2[Th(K)]);
	  const FElement KK3(Vh3[Th(K)]);
	  const FElement KK4(Vh4[Th(K)]);
	  const FElement KK5(Vh5[Th(K)]);
	  const FElement KK6(Vh6[Th(K)]);
	
	   vxx[it] = KK1(Cdg_hat,*fe1.x(),componante1,0);
	   vyx[it] = KK2(Cdg_hat,*fe2.x(),componante2,0); 
	   vyy[it] = KK3(Cdg_hat,*fe3.x(),componante3,0); 
	   vzx[it] = KK4(Cdg_hat,*fe4.x(),componante4,0);
	   vzy[it] = KK5(Cdg_hat,*fe5.x(),componante5,0);
	   vzz[it] = KK6(Cdg_hat,*fe6.x(),componante6,0);
	    
	}
	for(int k=0; k<Th.nt; k++){
	  fprintf(popenstream,"%f %f %f %f %f %f\n",vxx[k],vyx[k],vyy[k],vzx[k],vzy[k],vzz[k]);
	}
      }
    }
    fprintf(popenstream,"End");
  }
  
  // fermeture du stream pour popen
  fclose(popenstream);
  delete [] commandline; 
  return valsortie;
}


//class  PopenMeditMesh3 : public OneOperator { public:  
//     PopenMeditMesh3() : OneOperator(atype<long>(),atype<pmesh3>()) {}
//  
//  E_F0 * code(const basicAC_F0 & args) const 
//  {
//	cout << " PopenMeditMesh3 je suis dans code(const basicAC_F0 & args) const" << endl; 
//	//cout << "args: " << args << endl;   
//	return  new  PopenMeditMesh3_Op(args,t[0]->CastTo(args[0])); 
//  }
//};




// datasolMesh3

class datasolMesh3_Op : public E_F0mps 
{
public:
  Expression eTh;
  Expression filename;
  
  struct Expression2 {
    long what; // 1 scalar, 2 vector, 3 symtensor
    long nbfloat; // 1 scalar, 3 vector (3D), 6 symtensor(3D)
    Expression e[6];
    Expression2() {e[0]=0; e[1]=0; e[2]=0; e[3]=0; e[4]=0; e[5]=0; what=0; nbfloat=0;};
    Expression &operator[](int i){return e[i];}
    double eval(int i,Stack stack) const  { 
    if (e[i]) {
      return GetAny< double >( (*e[i])(stack) );
    }
    else 
      return 0;
    }
  };
  vector<Expression2> l; 
  static const int n_name_param =7;  
  static basicAC_F0::name_and_type name_param[];
  Expression nargs[n_name_param];
  long    arg(int i,Stack stack,long a ) const{ return nargs[i] ? GetAny< long >( (*nargs[i])(stack) ): a;}
 
public:
  datasolMesh3_Op(const basicAC_F0 &  args,  Expression fname, Expression tth) 
    : eTh(tth), filename(fname)
  {
    int nbofsol;
    int ddim=3;
    int stsize=6;
    // cout << "construction data medit solution" << args << endl;
    args.SetNameParam(n_name_param,name_param,nargs); 

    const E_Array * a0=0;
    const E_Array * a1=0;
    const E_Array * a2=0;

    if(nargs[0])  a0  = dynamic_cast<const E_Array *>(nargs[1]);
    if(nargs[2])  a1  = dynamic_cast<const E_Array *>(nargs[3]);
    if(nargs[4])  a2  = dynamic_cast<const E_Array *>(nargs[5]);
   
    nbofsol =0;
    cout << "l.size()= "<< nbofsol << endl;
    if(a0){
      nbofsol = nbofsol+a0->size();
    }
    if(a1){
      if(a1->size()%ddim !=0)
	CompileError("the number of element of parameters vector is not correct");
       nbofsol = nbofsol+a1->size()/ddim;
    }    

    if(a2){
      if(a2->size()%stsize !=0)
	CompileError("the number of element of parameters symtensor is not correct");
       nbofsol = nbofsol+a2->size()/stsize;
    }    
    cout << "l.size()= "<< nbofsol << endl;
    
    if(nbofsol == 0)
      CompileError("there is no solution");
    
    //vector<Expression2> l1[nbofsol]; 
    l= vector<Expression2> (nbofsol);
    cout << "l.size()= " << l.size() << endl;
    assert(nbofsol==l.size());
    size_t i=0;

    if(a0){
      for(size_t ii=0; ii <a0->size(); ii++){
	l[i].what=1;
	l[i].nbfloat=1;
	l[i][0]=to<double>( (*a0)[ii]);
	i++;
      }
    }
    if(a1){
      for(size_t ii=0; ii <a1->size()%ddim; ii++){
	l[i].what=2;
	l[i].nbfloat=ddim;
	for(int j=0; j<ddim; j++){
	  l[i][j] = to<double>( (*a1)[ii*ddim+j]);
	}
	i++;
      }
    }
    if(a2){
      for(size_t ii=0; ii <a2->size()%stsize; ii++){
	l[i].what=3;
	l[i].nbfloat=stsize;
	for(int j=0; j<stsize; j++){
	  l[i][j] = to<double>( (*a2)[ii*stsize+j]);
	}
	i++;
      }
    }
    
    //  cout << "fin construction meditmesh" << args << endl;
  }
  
  AnyType operator()(Stack stack)  const ;
};

basicAC_F0::name_and_type datasolMesh3_Op::name_param[]= {
  {  "nbscalar",&typeid(long)},
  {  "scalar", &typeid(E_Array)},
  {  "nbvector",&typeid(long)},
  {  "vector", &typeid(E_Array)},
  {  "nbsymtensor",&typeid(long)},
  {  "symtensor", &typeid(E_Array)},
  {  "order",&typeid(long)}
};

AnyType datasolMesh3_Op::operator()(Stack stack)  const 
{ 
  MeshPoint *mp(MeshPointStack(stack)) , mps=*mp;
  Mesh3 * pTh= GetAny<Mesh3 *>((*eTh)(stack));
  string * ffname= GetAny<string *>( (*filename)(stack) );
  ffassert(pTh);
  Mesh3 &Th=*pTh;

  int nt = Th.nt;
  int nv = Th.nv;

  int nbtype=l.size();
  int nbsol;
  int solnbfloat;
  int TypTab[l.size()];
 
  int resultorder= arg(6, stack, 1);
  long longdefault;

  int ver = GmfFloat, outm;
  // determination de TypTab
  solnbfloat=0;
  for (size_t i=0;i<l.size();i++){
    TypTab[i]=l[i].what;
    solnbfloat=solnbfloat+l[i].nbfloat;
  }
  float *OutSolTab = new float[solnbfloat];

  // determination de OutSolTab

  if ( !(outm = GmfOpenMesh(ffname->c_str(),GmfWrite,ver,3)) ) {
    cerr <<"  -- Mesh3::Save  UNABLE TO OPEN  :"<< filename << endl;
    exit(1);
  }
 
  if(resultorder==0){
    // ordre 0
    nbsol = nt;  
    
    KN<double> valsol(solnbfloat*nbsol);
 
    MeshPoint *mp3(MeshPointStack(stack)); 
    R3 Cdg_hat = R3(1./4.,1./4.,1./4.);  

    for (int it=0;it<nt;it++){
      int h=0;
      const Tet & K(Th.elements[it]);
      mp3->set( Th, K(Cdg_hat), Cdg_hat, K, K.lab);
	  
      for(size_t i=0;i<l.size();i++){
	for(size_t j=0;j<l[i].nbfloat;j++){
	  valsol[it*solnbfloat+h] = l[i].eval(j,stack);
	  h=h+1;
	}
      } 
      assert(solnbfloat==h);
    }


    GmfSetKwd(outm,GmfSolAtTetrahedra, nbsol, nbtype, TypTab);
    for (int k=0; k<nbsol; k++){
      for (int i=0; i<solnbfloat ;i++){
	OutSolTab[i] =  valsol(k*solnbfloat+i);
      }
      GmfSetLin(outm, GmfSolAtTetrahedra, OutSolTab); 
    }
  }
  if(resultorder==1){
    // ordre 1
    nbsol = nv;  

    KN<double> valsol(solnbfloat*nbsol);
    KN<int> takemesh(nbsol);
    MeshPoint *mp3(MeshPointStack(stack)); 
    R3 Cdg_hat = R3(1./4.,1./4.,1./4.);
    takemesh=0;
    for (int it=0;it<nt;it++){
      for(int iv=0;iv<4;iv++){
	int i=Th(it,iv);

	if(takemesh[i]==0){
	    mp3->setP(&Th,it,iv);
	    int h=0;

	    for(size_t ii=0;ii<l.size();ii++){
	      for(size_t j=0;j<l[ii].nbfloat;j++){
		//cout << "ii=" << ii << " j=" << j<< endl;
		valsol[i*solnbfloat+h] = l[ii].eval(j,stack);
		h=h+1;
	      }
	    } 
	    assert(solnbfloat==h);
	    takemesh[i] = takemesh[i]+1;
	}
      }
    }


    GmfSetKwd(outm,GmfSolAtVertices, nbsol, nbtype, TypTab);
    for (int k=0; k<nbsol; k++){
      for (int i=0; i<solnbfloat ;i++){
	OutSolTab[i] =  valsol(k*solnbfloat+i);
      }
      GmfSetLin(outm, GmfSolAtVertices, OutSolTab);
    }
  }
  GmfCloseMesh(outm);
  return longdefault;
}
  
class  datasolMesh3 : public OneOperator { public:  
  datasolMesh3() : OneOperator(atype<long>(),atype<string *>(),atype<pmesh3>()) {}
  
  E_F0 * code(const basicAC_F0 & args) const 
  {
	cout << " PopenMeditMesh3 je suis dans code(const basicAC_F0 & args) const" << endl; 
	//cout << "args: " << args << endl;   
	return  new  datasolMesh3_Op(args,t[0]->CastTo(args[0]),t[1]->CastTo(args[1]) ); 
  }
};
// datasolMesh3VhALL
template<class v_fes>  
class  datasolMesh3VhALL : public OneOperator { public:  
  datasolMesh3VhALL() : OneOperator(atype<long>(),atype<string *>(),atype<pmesh3>()) {}
  
   class datasolMesh3VhALL_Op : public E_F0mps
   {
   public:
     Expression eTh;
     Expression filename;
     struct Expression2 { 
       long what; // 1 scalar, 2 vector, 3 symtensor
       long nbfloat; // 1 scalar, 3 vector (3D), 6 symtensor(3D)
       Expression e[6];
       Expression2() {e[0]=0; e[1]=0; e[2]=0; e[3]=0; e[4]=0; e[5]=0; what=0; nbfloat=0;};
       Expression &operator[](int i){return e[i];}
       double eval(int i,const Mesh3 &Th, const Tet &K, const R3 &PHat,Stack stack) const;  
       
     };
     vector<Expression2> l; 
     static const int n_name_param =7;  
     static basicAC_F0::name_and_type name_param[];
     Expression nargs[n_name_param];
     long    arg(int i,Stack stack,long a ) const{ return nargs[i] ? GetAny< long >( (*nargs[i])(stack) ): a;}
     
   public:
     datasolMesh3VhALL_Op(const basicAC_F0 &  args,  Expression fname, Expression tth) 
       : eTh(tth), filename(fname)
     {
        int nbofsol;
	int ddim=3;
	int stsize=6;
	// cout << "construction data medit solution" << args << endl;
       args.SetNameParam(n_name_param,name_param,nargs);    
       const E_Array * a0=0;
       const E_Array * a1=0;
       const E_Array * a2=0;
       if(nargs[0])  a0  = dynamic_cast<const E_Array *>(nargs[1]);
       if(nargs[2])  a1  = dynamic_cast<const E_Array *>(nargs[3]);
       if(nargs[4])  a2  = dynamic_cast<const E_Array *>(nargs[5]);

        nbofsol =0;
	cout << "l.size()= "<< nbofsol << endl;
	if(a0){
	  nbofsol = nbofsol+a0->size();
	}
	if(a1){
	  if(a1->size()%ddim !=0)
	    CompileError("the number of element of parameters vector is not correct");
	  nbofsol = nbofsol+a1->size()/ddim;
	}    
	
	if(a2){
	  if(a2->size()%stsize !=0)
	    CompileError("the number of element of parameters symtensor is not correct");
	  nbofsol = nbofsol+a2->size()/stsize;
	}    
	cout << "l.size()= "<< nbofsol << endl;
	
	if(nbofsol == 0)
	  CompileError("there is no solution");
	
	//vector<Expression2> l1[nbofsol]; 
	l= vector<Expression2> (nbofsol);
	cout << "l.size()= " << l.size() << endl;
	assert(nbofsol==l.size());
	size_t i=0;
	
	if(a0){
	  for(size_t ii=0; ii <a0->size(); ii++){
	    l[i].what=1;
	    l[i].nbfloat=1;
	    l[i][0] = to< pair< FEbase<double,v_fes> *,int > >( (*a0)[ii]);
	    i++;
	  }
	}
	if(a1){
	  for(size_t ii=0; ii <(a1->size()/ddim); ii++){
	    l[i].what=2;
	    l[i].nbfloat=ddim;
	    for(int j=0; j<ddim; j++){
	      l[i][j] = to< pair< FEbase<double,v_fes> *,int > >( (*a1)[ii*ddim+j]);
	    }
	    i++;
	  }
	}
	if(a2){
	  for(size_t ii=0; ii <(a2->size()/stsize); ii++){
	    l[i].what=3;
	    l[i].nbfloat=stsize;
	    for(int j=0; j<stsize; j++){
	      l[i][j] = to< pair< FEbase<double,v_fes> *,int > >( (*a2)[ii*stsize+j]);
	    }
	    i++;
	  }
	}

	//cout << "fin construction meditmesh" << args << endl;
     }
     
     AnyType operator()(Stack stack)  const ;
   };
  
  E_F0 * code(const basicAC_F0 & args) const 
  {
	cout << " PopenMeditMesh3 je suis dans code(const basicAC_F0 & args) const" << endl; 
	//cout << "args: " << args << endl;   
	return  new  datasolMesh3VhALL_Op(args,t[0]->CastTo(args[0]),t[1]->CastTo(args[1]) ); 
  }

};
template<class v_fes>
double datasolMesh3VhALL<v_fes>::datasolMesh3VhALL_Op::Expression2::eval(int i,const Mesh3 &Th, const Tet &K, const R3 &PHat, Stack stack) const{
 typedef typename  v_fes::pfes pfes;
 typedef typename  v_fes::FESpace FESpace;
 typedef typename  FESpace::Mesh Mesh;
 typedef typename  FESpace::FElement FElement;
 typedef typename  Mesh::Element Element;
 typedef typename  Mesh::Vertex Vertex;  
 typedef typename  Mesh::RdHat RdHat;  
 typedef typename  Mesh::Rd Rd;  
 if (e[i]) {
   pair< FEbase<double,v_fes> *  ,int> pp1=GetAny< pair< FEbase<double,v_fes> *,int> >( (*e[i])(stack) );
   FEbase<double,v_fes> & fe1( *pp1.first);
   int componante1=pp1.second;
   const  FESpace & Vh1(*fe1.Vh);
   
   const FElement KK1(Vh1[Th(K)]);
   
   return KK1(PHat,*fe1.x(),componante1,0);
 }
 else 
   return 0;
}

template<class v_fes>
basicAC_F0::name_and_type datasolMesh3VhALL<v_fes>::datasolMesh3VhALL_Op::name_param[]= {
  {  "nbscalar",&typeid(long)},
  {  "scalar", &typeid(E_Array)},
  {  "nbvector",&typeid(long)},
  {  "vector",&typeid(E_Array)},
  {  "nbsymtensor",&typeid(long)},
  {  "symtensor", &typeid(E_Array)},
  {  "order",&typeid(long)},
};

template<class v_fes>
AnyType  datasolMesh3VhALL<v_fes>::datasolMesh3VhALL_Op::operator()(Stack stack)  const 
{ 
  typedef typename  v_fes::pfes pfes;
  typedef typename  v_fes::FESpace FESpace;
  typedef typename  FESpace::Mesh Mesh;
  typedef typename  FESpace::FElement FElement;
  typedef typename  Mesh::Element Element;
  typedef typename  Mesh::Vertex Vertex;  
  typedef typename  Mesh::RdHat RdHat;  
  typedef typename  Mesh::Rd Rd;  


  MeshPoint *mp(MeshPointStack(stack)) , mps=*mp;
  Mesh3 * pTh= GetAny<Mesh3 *>((*eTh)(stack));
  string * ffname= GetAny<string *>( (*filename)(stack) );
  ffassert(pTh);
  Mesh3 &Th=*pTh;

  int nt = Th.nt;
  int nv = Th.nv;
 
  int nbsol;
  int solnbfloat;
  int resultorder= arg(6, stack, 1);

  long longdefault;

  int ver = GmfFloat, outm;
  int nbtype = l.size();
  int TypTab[nbtype];

  solnbfloat=0;
  for (size_t i=0;i<l.size();i++){
    TypTab[i]=l[i].what;
    solnbfloat=solnbfloat+l[i].nbfloat;
    cout << "i="<< i << " TypTab[i]=" << TypTab[i] << " solnbfloat=" << solnbfloat << endl; 
  }
  float *OutSolTab = new float[solnbfloat];

 
  // determination de OutSolTab
  if ( !(outm = GmfOpenMesh(ffname->c_str(),GmfWrite,ver,3)) ) {
    cerr <<"  -- Mesh3::Save  UNABLE TO OPEN  :"<< filename << endl;
    exit(1);
  }
 
  if(resultorder==0){
    // ordre 0
    nbsol = nt;  
  
    KN<double> valsol(solnbfloat*nbsol);
    
    MeshPoint & mp3 = *MeshPointStack(stack); 
    R3 Cdg_hat = R3(1./4.,1./4.,1./4.);  
    
    for (int it=0;it<nt;it++){
      int h=0;
      const Tet & K(Th.elements[it]);
      mp3.set( Th, K(Cdg_hat), Cdg_hat, K, K.lab);
      
      for(size_t i=0;i<l.size();i++){
	for(size_t j=0;j<l[i].nbfloat;j++){
	  valsol[it*solnbfloat+h] = l[i].eval(j,Th,K,Cdg_hat,stack);
	  h=h+1;
	}
      } 
      assert(solnbfloat==h);
    }
  
    GmfSetKwd(outm,GmfSolAtTetrahedra, nbsol, nbtype, TypTab);
    for (int k=0; k<nbsol; k++){
      for (int i=0; i<solnbfloat ;i++){
	OutSolTab[i] =  valsol(k*solnbfloat+i);
      }
      GmfSetLin(outm, GmfSolAtTetrahedra, OutSolTab); 
    }
  }
  if(resultorder==1){
    // ordre 1
    nbsol = nv;  

    KN<double> valsol(solnbfloat*nbsol);
    KN<int> takemesh(nbsol);
    MeshPoint *mp3(MeshPointStack(stack)); 
    takemesh=0;
    for (int it=0;it<nt;it++){
      const Tet & K(Th.elements[it]);
      for(int iv=0;iv<4;iv++){
	int i=Th(it,iv);

	if(takemesh[i]==0){
	    mp3->setP(&Th,it,iv);
	    int h=0;

	    for(size_t ii=0;ii<l.size();ii++){
	      for(size_t j=0;j<l[ii].nbfloat;j++){
		//cout << "ii=" << ii << " j=" << j<< endl;
		valsol[i*solnbfloat+h] = l[ii].eval(j,Th,K,mp3->PHat,stack);
		h=h+1;
	      }
	    } 
	    assert(solnbfloat==h);
	    takemesh[i] = takemesh[i]+1;
	}
      }
    }


    GmfSetKwd(outm,GmfSolAtVertices, nbsol, nbtype, TypTab);
    for (int k=0; k<nbsol; k++){
      for (int i=0; i<solnbfloat ;i++){
	OutSolTab[i] =  valsol(k*solnbfloat+i);
      }
      GmfSetLin(outm, GmfSolAtVertices, OutSolTab);
    }
  }

  GmfCloseMesh(outm);
  return longdefault;
}



//  truc pour que la fonction 
// Init::Init() soit appele a moment du chargement dynamique
// du fichier 
//  

class Init { public:
  Init();
};

static Init init;  //  une variable globale qui serat construite  au chargement dynamique 

Init::Init(){  // le constructeur qui ajoute la fonction "splitmesh3"  a freefem++ 
  
  typedef Mesh *pmesh;
  //typedef Mesh2 *pmesh2;
  typedef Mesh3 *pmesh3;
  
  if (verbosity)
    cout << " load:popen.cpp  " << endl;
  
  // 2D
  Global.Add("meditmeshsol","(",new PopenMeditMesh);
  Global.Add("savesol","(",new datasolMesh2);
  Global.Add("meditmeshsolVh","(",new PopenMeditMesh2ALL);
  Global.Add("savesolVh","(",new datasolMesh2VhALL);
  
  // 3D
  Global.Add("meditmeshsol","(",new PopenMeditMesh3);
  Global.Add("savesol","(",new datasolMesh3);
  Global.Add("meditmeshsolVh","(",new PopenMeditMesh3ALL<v_fes3>);
  Global.Add("savesolVh","(",new datasolMesh3VhALL<v_fes3>);
  
}
