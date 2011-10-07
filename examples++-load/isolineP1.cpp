// ORIG-DATE: Novembre 2008
// -*- Mode : c++ -*-
//
// SUMMARY  :  
// USAGE    : LGPL      
// ORG      : LJLL Universite Pierre et Marie Curie, Paris,  FRANCE 
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

/*
calcul demander par F. Hecht
*/

#ifndef WITH_NO_INIT
#include "ff++.hpp"
#endif

//  TransfoMesh_v2.cpp
using namespace std;
// LayerMesh.cpp
// buildlayer.cpp
// rajout global

#include <set>
#include <vector>
#include "msh3.hpp"

using namespace  Fem2D;


// fonction determinant les points d'intersection
class ISOLINE_P1_Op : public E_F0mps 
{
public:
  Expression eTh,filename,ff;
  static const int n_name_param = 1; // 
  static basicAC_F0::name_and_type name_param[] ;
  Expression nargs[n_name_param];

  double  arg(int i,Stack stack,double a) const{ return nargs[i] ? GetAny< double >( (*nargs[i])(stack) ): a;}
 
public:
  ISOLINE_P1_Op(const basicAC_F0 &  args,Expression tth, Expression fff) 
    : eTh(tth), filename(fff)
  {
    args.SetNameParam(n_name_param,name_param,nargs);
 
    if (  BCastTo<double>(args[2])  )
      {	 
	ff=to<double>( args[2] );  
      }
    else {
      CompileError("no function to isolines \n");
    }    
    if( !nargs[0]) 
      CompileError("no isolines selected \n");
  } 
    
  AnyType operator()(Stack stack)  const ;
};

basicAC_F0::name_and_type ISOLINE_P1_Op::name_param[]= {
  {  "iso", &typeid(double)}
};

AnyType ISOLINE_P1_Op::operator()(Stack stack)  const 
{
  MeshPoint *mp(MeshPointStack(stack)) , mps=*mp;
  Mesh * pTh= GetAny<Mesh *>((*eTh)(stack));
  string * ffname= GetAny<string *>( (*filename)(stack) );
  ffassert(pTh);
  Mesh &Th=*pTh;
  Mesh *m= pTh;   // question a quoi sert *m ??
  int nbv=Th.nv; // nombre de sommet 
  int nbt=Th.nt; // nombre de triangles
  int nbe=Th.neb; // nombre d'aretes fontiere
  long valsortie;
  // value of isoline
  double isovalue;
  isovalue = GetAny< double >( (*nargs[0])(stack) ); 

  // evaluation de la fonction aux sommets
  KN<int> takevertex(Th.nv);
  KN<double> tff(Th.nv);

  MeshPoint *mp3(MeshPointStack(stack)); 
  
  takevertex=0;
  // loop over triangle 
  for (int it=0;it<Th.nt;++it){
    for( int iv=0; iv<3; ++iv){
      int i=Th(it,iv);  
      
      if(takevertex[i]==0){
	mp3->setP(&Th,it,iv);
	tff[i]=GetAny<double>((*ff)(stack));
	takevertex[i] = takevertex[i]+1;
      }
      
    }
  }

  
  // calcul des isolines dans les triangles
  KN<double> EdgeIter(3*Th.nt);
  KN<int> taketriangle(2*Th.nt);
  
  for(int ii=0; ii< 3*Th.nt; ii++){
    EdgeIter[ii] = -1;
  }
  for(int ii=0; ii< 2*Th.nt; ii++){
    taketriangle[ii] = -1;
  }
  
  int VertexIso[Th.nv];
  int VertexIsoTri[2*Th.nv];

  for(int ii=0; ii< Th.nv; ii++)
    VertexIso[ii] = 0;
  
  for(int ii=0; ii< 2*Th.nv; ii++)
    VertexIsoTri[ii] = 0;

  int NbIsoVertex    = 0;
  int NbIsoNonVertex = 0;

  for(int it=0; it< Th.nt; it++){
    const  Mesh::Triangle & K(Th.t(it));

    int nkeq[3];
    int ivertex = 0;
    int mark[3];

    int im=0;
    for(int ii=0; ii<3; ii++){ 
      int j0,j1;
      nkeq[ii] = 0;
      mark[ii] = 0;
      Th.VerticesNumberOfEdge( K, ii, j0, j1); 
      //cout << "it=" << it << " ii= " << ii << " j0= " << j0 << "j1= " << j1 << endl;

      double fi = tff[j0];
      double fj = tff[j1];
      double xf = isovalue;

      if( Abs( tff[j0]- isovalue) < 1.e-11*(Abs(isovalue)+Abs(tff[j0])) ) {
	nkeq[ii] = 1;
	mark[ii] = 1;
	EdgeIter[3*it+ii] =0.;
	ivertex++; 
      }
      else{
	if( Abs( tff[j1]- isovalue) < 1.e-11*(Abs(isovalue)+Abs(tff[j1])) ) {
	  mark[(ii+1)%3] = 1;
	}
	else
	  if( ((fi <= xf)&&(fj>=xf)) || ((fi>=xf)&&(fj<=xf)) ){
	    int eo;
	    eo=ii;
	    int ito=Th.ElementAdj(it,eo);
	    
	    mark[ii] = 1;
	    im++;
	    
	    double xlam = (fi-xf)/(fi-fj);
	   
	    if(ito<0){
	      EdgeIter[3*it+ii] = xlam;	     
	      NbIsoNonVertex++;
	    }
	    else if( 3*it+ii <= 3*ito+eo || it == ito  ){
	      EdgeIter[3*it+ii] = xlam;	     
	      NbIsoNonVertex++;
	    }       
	    
	    if(it <=10 && verbosity>10) cout << "vertex (it="<< it << ", i="<< ii << ") :: " << j0 << " " << j1 <<" xlam= "<< xlam << endl;
	  }

      }
    }
    if( ivertex == 3){
      cerr << " A triangle is a isovalue " << endl;
      exit(1);
    }

    else if( ivertex == 2){
      //*  search positive triangle *//
      // deux possibilités
      for(int iii=0;iii<3;iii++){
	if( nkeq[iii] != 1 ){
	  int j0,j1;
	  Th.VerticesNumberOfEdge( K, iii, j0, j1); 
	  R fi = tff[ j0 ];

	  if(  fi > isovalue ){
	     for(int i=0;i<3;i++) // for the  3 edges 
	       {		 
		 if( mark[i] == 1 && mark[(i+1)%3] == 1 ){
		   int jj0,jj1;
		   Th.VerticesNumberOfEdge( K, i, jj0, jj1); 

		   taketriangle[2*it+1] = i;
		   taketriangle[2*it] = (i+1)%3;
		   
		   VertexIso[ jj0 ]++;
		   VertexIso[ jj1 ]++;
		   
		   if( VertexIso[ jj0 ] > 2 || VertexIso[ jj1 ] > 2 )
		     cerr << " erreur: Le sommet au passe l'iso valeur est connecter à plus de 2 iso points " << endl;
		   else{
		     VertexIsoTri[ 2*jj0+VertexIso[ jj0 ] -1 ] = it;
		     VertexIsoTri[ 2*jj1+VertexIso[ jj1 ] -1 ] = it;
		   }
		     
		 }
	       }
	     
	   }
	}
      }
    }
    else if( im == 1 && ivertex == 1){
      for(int i=0;i<3;i++) // for the  3 edges 
	{
	  if( mark[i] == 1 && mark[(i+1)%3] == 1 ){
	    int jj0,jj1;
	    Th.VerticesNumberOfEdge( K, i, jj0, jj1);
	    
	    // determination
	    if(1==nkeq[i]){
	      EdgeIter[3*it+i] = 0.;
	      VertexIso[ jj0 ]++;
	      if( VertexIso[ jj0 ] > 2 )
		cerr << " erreur: Le sommet au passe l'iso valeur est connecter à plus de 2 iso points " << endl;
	      else{
		VertexIsoTri[ 2*jj0+VertexIso[ jj0 ] -1 ] = it;
	      }
	    }
	    else{
	      EdgeIter[3*it+(i+1)%3] = 0.;
	      VertexIso[ jj1 ]++;	
	      if( VertexIso[ jj1 ] > 2 )
		cerr << " erreur: Le sommet au passe l'iso valeur est connecter à plus de 2 iso points " << endl;
		else{
		  VertexIsoTri[ 2*jj1+VertexIso[ jj1 ] -1 ] = it;
		}
	    }
	    int jj00,jj10;
	    Th.VerticesNumberOfEdge( K, (i+2)%3, jj00, jj10);
	    R fi2 = tff[ jj00 ];
	    
	    if( fi2 > isovalue ){
	      taketriangle[2*it]   = (i+1)%3;
	      taketriangle[2*it+1] = i;
	    }
	    else{
	      taketriangle[2*it]   = i;
	      taketriangle[2*it+1] = (i+1)%3;
	    }
	    
	  

	  }
	}
      
    }
    else if( im == 2){
      
      for(int i=0;i<3;i++) // for the  3 edges 
	{
	  
	  if( mark[i] == 1 && mark[(i+1)%3] == 1 ){
	    int jj0,jj1;
	    Th.VerticesNumberOfEdge( K, i, jj0, jj1); 
	    R fi1 = tff[ jj0 ];
	   
	    if( it < 10 )  cout << "vertex (it="<< it << ", i="<< i << ") :: " << jj0 << " " << jj1 << endl;
	    if( fi1 > isovalue ){
	      taketriangle[2*it]   = (i+1)%3;
	      taketriangle[2*it+1] = i;
	    }
	    else{
	      taketriangle[2*it]   = i;
	      taketriangle[2*it+1] = (i+1)%3;
	    }
	    
	    
	    Th.VerticesNumberOfEdge( K, (i+1)%3, jj0, jj1); 	    
	    if( it < 10 ) cout << "vertex (it="<< it << ", i="<< (i+1)%3 << ") :: " << jj0 << " " << jj1 << endl;

	  }
	}
    }
    
    
    
  }
 
  
  //#################################  
  int NbInterBord=0;
  KN<int> ElementLink(Th.nt+Th.neb);
  for(int it=0; it< Th.nt+Th.neb; it++)
    ElementLink[it]=-1;

  
  for(int it=0; it< Th.nt; it++){
    if( taketriangle[2*it] < 0 ) continue;

    const  Mesh::Triangle & K(Th.t(it));
    

    int ii  = 2*it;
    int eT1 = taketriangle[2*it];
    int eT2 = taketriangle[2*it+1];
    int eo1 = eT1;
    int eo2 = eT2;
    int adjeT1 = Th.ElementAdj(it,eo1);
    int adjeT2 = Th.ElementAdj(it,eo2);
    
    /*
      if( VertexIso[ numv ] == 2 ){     
      if( VertexIsoTri[ 2*numv ] == it  ){
      newit = VertexIsoTri[ 2*numv+1 ];
      }
      else 
      newit = VertexIsoTri[ 2*numv ];
      }
    */
    
    // link with previous element 
    {
      int jj00,jj10;
      Th.VerticesNumberOfEdge( K, eT1, jj00, jj10);
      
      int numv=jj00;
      //cout << "eT1, jj00, jj10 " << eT1 << " " << jj00 << " " << jj10 << " VertexIso[ numv ] "<< VertexIso[ numv ] << endl;

      if( VertexIso[ numv ] == 2 ){     
	if( VertexIsoTri[ 2*numv ] == it  ){
	  ElementLink[ VertexIsoTri[ 2*numv+1 ] ] = it;
	  //newit = VertexIsoTri[ 2*numv+1 ];
	}
	else 
	  ElementLink[ VertexIsoTri[ 2*numv] ] = it;
	//newit = VertexIsoTri[ 2*numv ];
      }
      else if( VertexIso[ numv ] == 1 ){	
	if(adjeT1 >= 0 && it!=adjeT1 ){
	  
	  // search the edge of it in the border
	  int eT3 = 3-eT1-eT2;
	  int eo3 = eT3;
	  int adjeT3 = Th.ElementAdj(it,eo3);
	  
	  if(adjeT3 < 0 || it==adjeT3 ){
	    int jj000,jj100;
	    Th.VerticesNumberOfEdge( K, eT3, jj000, jj100);  
	    cout << " bug une iso vertex definis une fois doit etre sur le bord ::  vertex " << numv << endl;
	    int ibeV = Th.NumberOfTheBoundaryEdge(jj000,jj100); // old
	 
	    ElementLink[ ibeV+Th.nt ] = it;
	    //NbInterBord++;
	  }
	  //exit(1);
	}
	else{
	  // border 
	  // int ibeV = Th.BorderElementAdj( jj00, jj10); //old
	  int ibeV = Th.NumberOfTheBoundaryEdge( jj00, jj10);
	  ElementLink[ibeV+Th.nt] = it;
	  //NbInterBord++;
	}
      
      }
      else
	if(adjeT1 >= 0 && it!=adjeT1 ){
	  ElementLink[adjeT1] = it;
	}
	else{
	  // border 
	  int ibeV = Th.NumberOfTheBoundaryEdge(jj00,jj10);
	  //ElementLink[Th.NumberOfTheBoundaryEdge(jj00,jj10)+Th.nt] = it; // old
	  ElementLink[ibeV+Th.nt] = it;
	  //NbInterBord++;
	}
    }
    // link with the next element
    {
      int jj00,jj10;
      Th.VerticesNumberOfEdge( K, eT2, jj00, jj10);     
      int numv=jj00;

      if( VertexIso[ numv ] == 2 ){     
	if( VertexIsoTri[ 2*numv ] == it  ){
	  ElementLink[ it ] = VertexIsoTri[ 2*numv+1 ];
	  //newit = VertexIsoTri[ 2*numv+1 ];
	}
	else 
	  ElementLink[ it ] = VertexIsoTri[ 2*numv];
	//newit = VertexIsoTri[ 2*numv ];
      }
      else if( VertexIso[ numv ] == 1 ){
	
	if(adjeT2 >= 0 && it!=adjeT2 ){

	  // search the edge of it in the border
	  int eT3 = 3-eT1-eT2;
	  int eo3 = eT3;
	  int adjeT3 = Th.ElementAdj(it,eo3);
	  
	  if(adjeT3 < 0 || it==adjeT3 ){
	    int jj000,jj100;
	    Th.VerticesNumberOfEdge( K, eT3, jj000, jj100);  
	    cout << " bug une iso vertex definis une fois doit etre sur le bord ::  vertex " << numv << endl;
	    int ibeV = Th.NumberOfTheBoundaryEdge(jj000,jj100); // old
	 
	    ElementLink[ it ] = ibeV+Th.nt;
	    NbInterBord++;
	  }
	  //exit(1);
	}
	else{
	  // border 
	  int ibeV = Th.NumberOfTheBoundaryEdge(jj00,jj10); // old
	  //int ibeV =  Th.BorderElementAdj( jj00, jj10);
	  ElementLink[ it ] = ibeV+Th.nt;
	  NbInterBord++;
	}
      
      }
      else{
	if(adjeT2 >= 0 && it!=adjeT2) 
	  ElementLink[it] = adjeT2;
	else{
	  // border
	  int ibeV = Th.NumberOfTheBoundaryEdge(jj00,jj10);
	  ElementLink[it] = ibeV+Th.nt;
	  // ElementLink[it] = Th.BorderElementAdj( jj00, jj10)+Th.nt; // old
	  NbInterBord++;
	}

      //for(int iijj=0; iijj<10; iijj++){
      //cout << "ElementLink["<< iijj <<"]=" << ElementLink[iijj] <<endl;
      //}
      }
    }
    /*
      if(adjeT2 >= 0 && it!=adjeT2) 
      ElementLink[it] = adjeT2;
      else{
      // border
      const  Mesh::Triangle & K(Th.t(it));
      int jj00,jj10;
      Th.VerticesNumberOfEdge( K, taketriangle[ii], jj00, jj10);
      ElementLink[it] = Th.TheBoundaryEdge(jj00,jj10)+Th.nt;
      //NbInterBord++;
      }
    */
  }
 

  int NbBordVertex=0;
  if(verbosity>10)cout << " NbInterBord = " << NbInterBord << endl;
  if(NbInterBord>0){
    //#################################
    // boucle sur le bord    
    // determination des points sur le bord  
     
    for(int ii=0; ii < Th.neb; ii++){
      // determination du sens du bord
      int edgebid;
      int ffbid   = Th.BoundaryElement( ii, edgebid );     // ii : number of edge => sortie :: ffbid = numero triangles, edgebid = numero edges
      int j0bid,j1bid;
      Th.VerticesNumberOfEdge( Th.t(ffbid), edgebid, j0bid, j1bid);
      // j0bid --> j1bid sens de parcours sur le triangle
      double fi = tff[j0bid];
      double fj = tff[j1bid];
      double xf = isovalue;

      // sens fi -> fj (same as local triangle) 
 
      //  fi++ et fj++  ==> lien ibe ++ next border
      //  fi==0 et fjj++ ==> lien triangle contenant cette vertex ++ vers next border
      //  fj==0 et fi++  ==> lien ibe vers triangle contenant cette vertex
      //  fi==0 et fj==0 ==> lien triangle contenant cette vertex vers - next bord  ( fjj >= 0 )
      //                                                               - triangle contenant cette vertex (fjj < 0)

      if( VertexIso[ j0bid ] > 0 ){	
	// cas fi == isovalue
	if( VertexIso[ j1bid ] > 0 ){     
	  if(verbosity>10)cout << "the edge is a isovalue :: link is previously computed "<< endl;
	  //assert( VertexIso[ j0bid ] == 2);
	  //assert( VertexIso[ j1bid ] == 2);
	}
	else if( fj >=  isovalue){
	  // cas fi iso  fj++
	  //assert( VertexIso[ j0bid ] == 2);
	  if( VertexIso[ j0bid ] == 2 ){
	    if( VertexIsoTri[2*j0bid] == ffbid ){
	      //ElementLink[ VertexIsoTri[2*j0bid+1] ] = Th.nt + Th.BorderElementAdj(j1bid,j0bid); // old orientation 
	      ElementLink[  Th.nt + Th.NumberOfTheBoundaryEdge(j0bid,j1bid) ] = VertexIsoTri[2*j0bid+1];
	    }
	    else{
	      //ElementLink[ VertexIsoTri[2*j0bid] ] = Th.nt + Th.BorderElementAdj(j1bid,j0bid);  // old orientation 
	      ElementLink[  Th.nt + Th.NumberOfTheBoundaryEdge(j0bid,j1bid) ] = VertexIsoTri[2*j0bid];
	    }
	    NbBordVertex++;
	  }
	  if( VertexIso[ j0bid ] == 1 ){
	    if(verbosity>10) cout << "j0bid, j1bid, ffbid " << j0bid << " "<< j1bid << " " << ffbid << endl; 
	    //assert( VertexIsoTri[2*j0bid] == ffbid );
	    //if( VertexIsoTri[2*j0bid] == ffbid ){
	    ElementLink[  Th.nt + Th.NumberOfTheBoundaryEdge(j0bid,j1bid) ] = VertexIsoTri[2*j0bid];
	    //}
	    NbBordVertex++;
	  }
	}
      }
      else if( fi >=  isovalue){
	// cas fi == isovalue++
	if( VertexIso[ j1bid ] > 0 ){
	  assert( VertexIso[ j1bid ] == 1);
	  // ElementLink[ ii+Th.nt ] = VertexIsoTri[2*j1bid]; // old version
	  ElementLink[ VertexIsoTri[2*j1bid] ] = ii+Th.nt;
	  ElementLink[ ii+Th.nt ] = Th.BorderElementAdj(j0bid,j1bid)+Th.nt;
	}
	else if( fj >=  isovalue){
	  //ElementLink[ii+Th.nt] = Th.nt + Th.BorderElementAdj(j1bid,j0bid); // old orientation
	  ElementLink[ii+Th.nt] = Th.nt + Th.BorderElementAdj(j0bid,j1bid);
	  NbBordVertex++;
	}
	//old version ::
	else{
	  ElementLink[ii+Th.nt] = Th.nt + Th.BorderElementAdj(j0bid,j1bid);
	  NbBordVertex++;
	}
      }

      /*
	if( fi >=  isovalue && Abs( fi - isovalue) > 1.e-11*(Abs(isovalue)+Abs(fi)) ){
	if( fj >= isovalue && Abs( fj - isovalue) > 1.e-11*(Abs(isovalue)+Abs(fj)) ){
	ElementLink[ii+Th.nt] = Th.nt + Th.BorderElementAdj(ii,j0bid,j1bid);
	}
	
	NbBordVertex++;
	}
      */
    }
    
    // rappel: le lien a été effectuer entre les bords
    //   -1 : pas de suivant
    //   positive et nulle : le suivant border element
    //   -2 : pas d'élement à prendre en compte. 
  }

  if(verbosity>10)
    {
      cout << " NbInterBord = " << NbInterBord << endl;
      for(int iijj=0; iijj<10; iijj++){
	cout << "ElementLink["<< iijj <<"]=" << ElementLink[iijj] <<endl;
      }
      
      cout << " NbInterBord = " << NbInterBord << endl;
      for(int iijj=Th.nt; iijj<Th.nt+Th.neb; iijj++){
	cout << "ElementLink["<< iijj <<"]=" << ElementLink[iijj] <<endl;
      }
    }

  //#################################

  int NbVertex=NbIsoNonVertex+NbBordVertex;
  for(int i=0; i< Th.nv; i++){
    if( VertexIso[i] > 0) 
      NbVertex++;
  }
  
  Vertex *VertexIsoP = new Vertex[NbVertex];
  
  KN<int> TriangleVu(Th.nt);
  for(int iii=0; iii< Th.nt; iii++)
    TriangleVu[iii]= -1;

  int inv = 0;
  int label = 0;
  for(int it1=0; it1< Th.nt; it1++){
    if( TriangleVu[it1] == 1 || taketriangle[2*it1]== -1) continue;
    if(verbosity>10) cout << "it1 = " << it1 << " Th.nt "<< Th.nt << endl;
    // First point is taken
    TriangleVu[it1] = 1;
    {
      const  Mesh::Triangle & K(Th.t(it1));
      int ii = 2*it1;
      int eT1 = taketriangle[ii];
      int eT2 = taketriangle[ii+1];   
      assert(eT1>=0);
      assert(eT2>=0);
      int jj00,jj10;
      Th.VerticesNumberOfEdge( K, eT1, jj00, jj10);
      
      if( EdgeIter[3*it1+eT1] > -0.1){ 
	double xlam = EdgeIter[3*it1+eT1] ;
	
	VertexIsoP[inv].x   = (1-xlam)*Th.vertices[ jj00 ].x + xlam*Th.vertices[ jj10 ].x; 
	VertexIsoP[inv].y   = (1-xlam)*Th.vertices[ jj00 ].y + xlam*Th.vertices[ jj10 ].y;
	VertexIsoP[inv].lab = label;
	
	inv++;
      }
      else{
	int eo1;
	eo1=eT1;
	int ito1=Th.ElementAdj(it1,eo1);
	
	double xlam = EdgeIter[3*ito1+eo1];
	VertexIsoP[inv].x   = xlam*Th.vertices[ jj00 ].x + (1-xlam)*Th.vertices[ jj10 ].x; 
	VertexIsoP[inv].y   = xlam*Th.vertices[ jj00 ].y + (1-xlam)*Th.vertices[ jj10 ].y;
	VertexIsoP[inv].lab = label;
	inv++;
      }
    }
    
    int it2=ElementLink[it1];
      
    for( int it=it2; it != it1; it=ElementLink[it]){   
      assert(it >=0);
     
      if( it< Th.nt){
	// sur un triangle
	if(verbosity>10) cout << "it = " << it << " <--->  it2=" << it2 << " Th.nt "<< Th.nt << endl;
	assert( TriangleVu[it] == -1);  
	TriangleVu[it] = 1;

	int ii = 2*it;
	int eT1 = taketriangle[ii];
	assert( eT1 >=0);
	//int eT2 = taketriangle[ii+1];
	const  Mesh::Triangle & K(Th.t(it));
	int jj00,jj10;
	Th.VerticesNumberOfEdge( K, eT1, jj00, jj10);
	
	if( EdgeIter[3*it+eT1] > -0.1){ 
	  double xlam = EdgeIter[3*it+eT1] ;
	  
	  VertexIsoP[inv].x   = (1-xlam)*Th.vertices[ jj00 ].x + xlam*Th.vertices[ jj10 ].x; 
	  VertexIsoP[inv].y   = (1-xlam)*Th.vertices[ jj00 ].y + xlam*Th.vertices[ jj10 ].y;
	  VertexIsoP[inv].lab = label;
	  
	  inv++;
	}
	else{
	  int eo1;
	  eo1=eT1;
	  assert( eo1 >=0);
	  int ito1=Th.ElementAdj(it,eo1);
	  
	  double xlam = EdgeIter[3*ito1+eo1];
	  VertexIsoP[inv].x   = xlam*Th.vertices[ jj00 ].x + (1-xlam)*Th.vertices[ jj10 ].x; 
	  VertexIsoP[inv].y   = xlam*Th.vertices[ jj00 ].y + (1-xlam)*Th.vertices[ jj10 ].y;
	  VertexIsoP[inv].lab = label;
	  inv++;
	}
	 
      }
      else{
	int ibe=it-Th.nt;
	if(verbosity>10) cout << "ibe = " << ibe << " <--->  it2=" << it2 << " Th.nt "<< Th.nt << endl;
	// sur le bord
	int edgebid;
	//int newit;
	int ffbid   = Th.BoundaryElement( ibe, edgebid );     // ii : number of edge => sortie :: ffbid = numero triangles, edgebid = numero edges
	int j0bid,j1bid;
	Th.VerticesNumberOfEdge( Th.t(ffbid), edgebid, j0bid, j1bid);
	if(verbosity>10) cout << "Edge Vertex Number "<< j0bid+1 << " " << j1bid+1 << " number of triangle " << ffbid << endl;
	if( taketriangle[2*ffbid+1] == edgebid && VertexIso[j1bid] == 0){
	  //if( taketriangle[2*ffbid+1] == edgebid && VertexIso[j0bid] == 0){ // old version
	  if( EdgeIter[3*ffbid+edgebid] > -0.1){ 
	    double xlam = EdgeIter[3*ffbid+edgebid] ;
	    
	    VertexIsoP[inv].x   = (1-xlam)*Th.vertices[ j0bid ].x + xlam*Th.vertices[ j1bid ].x; 
	    VertexIsoP[inv].y   = (1-xlam)*Th.vertices[ j0bid ].y + xlam*Th.vertices[ j1bid ].y;
	    VertexIsoP[inv].lab = label;
	    
	    inv++;
	  }
	  else{
	    int eo1;
	    eo1=edgebid;
	    int ito1=Th.ElementAdj(ffbid,eo1);
	    
	    double xlam = EdgeIter[3*ito1+eo1];
	    VertexIsoP[inv].x   = xlam*Th.vertices[ j0bid ].x + (1-xlam)*Th.vertices[ j1bid ].x; 
	    VertexIsoP[inv].y   = xlam*Th.vertices[ j0bid ].y + (1-xlam)*Th.vertices[ j1bid ].y;
	    VertexIsoP[inv].lab = label;
	    inv++;
	  }
	  
	}
	else{
	  // old  version j0bid à la place de j1bid
	  VertexIsoP[inv].x   = Th.vertices[ j1bid ].x;
	  VertexIsoP[inv].y   = Th.vertices[ j1bid ].y;
	  VertexIsoP[inv].lab = label;
	  inv++;
	}
	
      }
    }
    label++;
  }
  

//   int TriangleVu[Th.nt];
//   for(int iii=0; iii< Th.nt; iii++)
//     TriangleVu[iii]= -1;
//   int label = 0;
//   int inv  = 0;

//   Vertex *VertexIsoP = new Vertex[NbVertex];
  
//   for(int it1=0; it1< Th.nt; it1++){
//     if( taketriangle[2*it1] < 0  || TriangleVu[it1] == 1 ) continue;

//     int it =it1;
//     do {
      
//       const  Mesh::Triangle & K(Th.t(it));

//       TriangleVu[it] = 1;
     
//       int ii = 2*it;
//       int eT1 = taketriangle[ii];
//       int eT2 = taketriangle[ii+1];
      
//       int jj00,jj10;
//       Th.VerticesNumberOfEdge( K, eT1, jj00, jj10);
      
//       if( EdgeIter[3*it+eT1] > -0.1){ 
// 	double xlam = EdgeIter[3*it+eT1] ;
	
// 	VertexIsoP[inv].x   = (1-xlam)*Th.vertices[ jj00 ].x + xlam*Th.vertices[ jj10 ].x; 
// 	VertexIsoP[inv].y   = (1-xlam)*Th.vertices[ jj00 ].y + xlam*Th.vertices[ jj10 ].y;
// 	VertexIsoP[inv].lab = label;

// 	inv++;
//       }
//       else{
// 	cerr << "bug in definition of vertex" << endl; 
// 	int eo1;
// 	eo1=eT1;
// 	int ito1=Th.ElementAdj(it,eo1);
	
// 	double xlam = EdgeIter[3*ito1+eo1];
// 	VertexIsoP[inv].x   = xlam*Th.vertices[ jj00 ].x + (1-xlam)*Th.vertices[ jj10 ].x; 
// 	VertexIsoP[inv].y   = xlam*Th.vertices[ jj00 ].y + (1-xlam)*Th.vertices[ jj10 ].y;
// 	VertexIsoP[inv].lab = label;
// 	inv++;
//       }
      
//       // on passe a eT2
//       int newit;
//       Th.VerticesNumberOfEdge( K, eT2, jj00, jj10);
//       int numv  =  jj00;
      
//       if( VertexIso[ numv ] == 2 ){     
// 	if( VertexIsoTri[ 2*numv ] == it  ){
// 	  newit = VertexIsoTri[ 2*numv+1 ];
// 	}
// 	else 
// 	  newit = VertexIsoTri[ 2*numv ];
//       }
//       else if( VertexIso[ numv ] == 1 ){     
// 	newit = -1;
//       }
//       else{
// 	int eo;
// 	eo=eT2;
// 	int ito=Th.ElementAdj(it,eo);
	 	
// 	if( (ito != it && taketriangle[2*ito] >= 0 ) && it >= 0){ 
// 	  newit = ito; 
// 	  if( taketriangle[2*ito] != eo){	    
// 	    exit(1);
// 	  }
// 	}
// 	else
// 	  newit = -1;
//       }         
//       if( newit < 0 ){
// 	 Th.VerticesNumberOfEdge( K, eT2, jj00, jj10);
      
// 	 if( EdgeIter[3*it+eT2] > -0.1){ 
// 	   double xlam = EdgeIter[3*it+eT2] ;
// 	   //cout << "vertex (it="<< it << ", eT1="<< eT1 << ") :: " << jj00 << " " << jj10 << endl; 
// 	   VertexIsoP[inv].x   = (1-xlam)*Th.vertices[ jj00 ].x + xlam*Th.vertices[ jj10 ].x; 
// 	   VertexIsoP[inv].y   = (1-xlam)*Th.vertices[ jj00 ].y + xlam*Th.vertices[ jj10 ].y;
// 	   VertexIsoP[inv].lab = label;
	  
// 	   inv++;
// 	 }
// 	 else{
// 	   //cerr << "bug in definition of vertex" << endl; 
// 	   // determination du triangle adjacents contenant e1
	   
// 	   //int ito1 = Th.TheAdjacencesLink[3*it+eT1]/3;
// 	   //int eo1  = Th.TheAdjacencesLink[3*it+eT1]%3;
// 	   int eo1;
// 	   eo1=eT2;
// 	   int ito1=Th.ElementAdj(it,eo1);
	   
// 	   double xlam = EdgeIter[3*ito1+eo1];
// 	   VertexIsoP[inv].x   = xlam*Th.vertices[ jj00 ].x + (1-xlam)*Th.vertices[ jj10 ].x; 
// 	   VertexIsoP[inv].y   = xlam*Th.vertices[ jj00 ].y + (1-xlam)*Th.vertices[ jj10 ].y;
// 	   VertexIsoP[inv].lab = label;
	
// 	   inv++;
// 	 }
//       }
//       //if(newit>=0)TriangleVu[it] = newit;
//       it = newit;
      
//     } while( it >= 0 && it!=it1 && taketriangle[2*it] >= 0  && TriangleVu[it] == -1 && inv <= NbVertex);

//     label=label+1;
//   } 
  if(verbosity) cout << " IsolineP1 :inv= " << inv << endl;
  if(verbosity) cout << "            NbVertex= " << NbVertex << endl;
  if(verbosity) cout << "            label =" << label << endl;
  assert(inv == NbVertex);
  if(verbosity>2) cout << "     file point \"" << ffname->c_str() <<"\""<< endl;
  FILE *fpoints = fopen(ffname->c_str(),"w");
  int lab=VertexIsoP[0].lab;
  for (int k=0; k<NbVertex; k++) {   
    //fprintf(fpoints,"%f %f %d\n",VertexIsoP[k].x,VertexIsoP[k].y,VertexIsoP[k].lab);
    if(VertexIsoP[k].lab != lab) fprintf(fpoints,"\n\n");
    fprintf(fpoints,"%f %f \n",VertexIsoP[k].x,VertexIsoP[k].y);
    lab=VertexIsoP[k].lab;
  }
  fclose(fpoints);

  delete [] VertexIsoP;

  return valsortie;
}

class  ISOLINE_P1: public OneOperator { public:  
typedef Mesh *pmesh;
  ISOLINE_P1() : OneOperator(atype<long>(),atype<pmesh>(),atype<string *>(),atype<double>() ) {}
  
  E_F0 * code(const basicAC_F0 & args) const 
  { 
    return  new ISOLINE_P1_Op( args, t[0]->CastTo(args[0]), t[1]->CastTo(args[1]) ); 
  }
};

class Init { public:
  Init();
};

static Init init;  //  une variable globale qui serat construite  au chargement dynamique 

Init::Init(){  // le constructeur qui ajoute la fonction "splitmesh3"  a freefem++ 
 
  typedef Mesh *pmesh;
  
  Global.Add("isolineP1","(",new ISOLINE_P1);

}
   
