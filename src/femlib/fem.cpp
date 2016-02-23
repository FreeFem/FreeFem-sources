// -*- Mode : c++ -*-
//
// SUMMARY  :      
// USAGE    :        
// ORG      : 
// AUTHOR   : Frederic Hecht
// E-MAIL   : hecht@ann.jussieu.fr
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
 */

extern long verbosity ;
extern long searchMethod; //pichon

#include <cmath>
#include  <cfloat>
#include <cstdlib>
#include "error.hpp"
#include <iostream>
#include <fstream>
//#include <strstream.h>
//using namespace std;  //introduces namespace std
#include "RNM.hpp"
#include "rgraph.hpp"
#include "Serialize.hpp"
#include "fem.hpp"
#include <set>
extern   long npichon2d, npichon3d;
extern   long npichon2d1, npichon3d1;

namespace Fem2D {
    
    
class SubMortar { 
	friend class Mesh;
	friend ostream & operator<<(ostream & f,const SubMortar & m);
	R  alpha; // angle in radian 
	R2 from,to;
	int k; //  triangle number
	int i; //  edge on triangle
	int sens; //       
		  //  int head;
public:
	    SubMortar() :alpha(0),k(0),i(0),sens(0){}
	SubMortar(const Vertex & s,const Vertex & ss,int kk,int ii,int se) 
	    : alpha(Theta((R2) ss - s)),from(s),to(ss),k(kk),i(ii),sens(se) {}
	R len2() const { return Norme2_2(to-from);}
	R len2(const R2 & A) const{ return (to-from,(A-from));}
	
	
	bool operator<(const SubMortar &  b){ return alpha < b.alpha ;}  // to sort
									 //  bool side(const Triangle &K) {}
    };
    
    ostream & operator<<(ostream & f,const SubMortar & m)
    { f << " a=" << m.alpha << " " << m.k << " " << m.i << " " << m.sens << " l^2=" <<m.len2() <<  ";" ;
	return f;}
    
    void Mesh::BuildBoundaryAdjacences()
    {
	if(!BoundaryAdjacencesHead)
	{
	    BoundaryAdjacencesHead = new int[nv];
	    BoundaryAdjacencesLink = new int[neb+neb];
	    for (int i=0;i<nv;i++) 
		BoundaryAdjacencesHead[i]=-1;
	    int j2=0;
	    for (int j=0;j<neb;j++)
		for (int k=0;k<2;k++,j2++)
		{  
		    int v = number(bedges[j][k]);
		    assert(v >=0 && v < nv);
		    BoundaryAdjacencesLink[j2]=BoundaryAdjacencesHead[v];
		    BoundaryAdjacencesHead[v]=j2;
		}
		    
	}
    }  
	void Mesh::ConsAdjacence()
    {
	    //  warning in the paper a mortar is the whole edge of the coarse triangle
	    //  here a mortar is a connected componand of he whole edge of the coarse triangle 
	    //   minus  the extremite of mortar
	    //  -----------
	    int NbCollision=0,NbOfEdges=0,NbOfBEdges=0,NbOfMEdges=0;
	    const char MaskEdge[]={1,2,4};
	    const char AddMortar[]={8,16,32};
	    //    reffecran();
	    //    cadreortho(0.22,0.22,.1);
	    if (neb) BuildBoundaryAdjacences();
            area=0;  // FH add nov 2010
            lenbord=0; // FH add nov 2010
	    if(TheAdjacencesLink) return; //
	    TheAdjacencesLink = new int[3*nt];
	    const int NbCode = 2*nv;
	    char * TonBoundary = new char [nt]; // the edge is 1 2 4   AddMortar = 8
	    
	    { int * Head = new int[NbCode];
		
		//  make the list
		int i,j,k,n,j0,j1;
		for ( i=0; i<NbCode; i++ )
		    Head[i]=-1; // empty list
		n=0; // make all the link
		for (i=0;i<nt;i++)
		{ 
		    Triangle & T=triangles[i];
		    area += T.area; // add FH nov 2010 
		    for( j=0; j<3; j++,n++ )
		    { 
			VerticesNumberOfEdge(T,j,j0,j1);
			k = j0+j1;
			TheAdjacencesLink[n]=Head[k]; 
			Head[k]=n; //            
		    }
		    
		}
		//
		if (neb==0) { // build boundary
		    for(int step=0;step<2;step++)
		    {
			neb=0;
			for (i=0;i<nt;i++)
			{ 
			    Triangle & T=triangles[i];
			    for( j=0; j<3; j++,n++ )
			    { 
				VerticesNumberOfEdge(T,j,j0,j1);
				int kk = 0,im=Min(j0,j1);             
				for (int n=Head[j0+j1]; n>=0; n=TheAdjacencesLink[n])
				{ int jj=n%3,ii=n/3, jj0,jj1;
				    VerticesNumberOfEdge(triangles[ii],jj,jj0,jj1);
				    if(im==Min(jj0,jj1)) // same edge 
					kk++;
				}
				if (kk==1) { 
				    if(step) bedges[neb].set(vertices,j0,j1,1);
				    neb++;
				}
				
			    }
			}
			if (step==0) {
			    if (verbosity) cout << " we build " << neb << " boundary edges" << endl;
			    bedges = new BoundaryEdge[neb];
			}
		    }
		    BuildBoundaryAdjacences(); 
		}
		for (int k=0;k<nt;k++) TonBoundary[k]=0;
		
		BoundaryEdgeHeadLink = new int[neb];       
		for (i=0;i<neb;i++)
		{  
		    BoundaryEdge & be(bedges[i]);
		    lenbord +=   be.length() ; // add Now 2010 FH 
		    int n;
		    int i0=number(be.vertices[0]);
		    int i1=number(be.vertices[1]);
		    throwassert(i0 >=0 && i0 < nv);
		    throwassert(i1 >=0 && i1 < nv);
		    throwassert(i1 != i0) ;
		    int im=Min(i0,i1);
		    BoundaryEdgeHeadLink[i]=-1; 
		    for ( n=Head[i0+i1]; n>=0; n=TheAdjacencesLink[n])
		    {
			int jj=n%3,ii=n/3, jj0,jj1;
			VerticesNumberOfEdge(triangles[ii],jj,jj0,jj1);
			if(im==Min(jj0,jj1)) // same edge 
			{
			    TonBoundary[n/3] += MaskEdge[n%3];
			    BoundaryEdgeHeadLink[i]=n;                  
			    if(i0==jj0) break; // FH 01072005 bon cote de l'arete 
					       // sinon on regard si cela existe?
			}
		    } 
		    if ( BoundaryEdgeHeadLink[i] <0 && verbosity) 
			cout << "   Attention l'arete frontiere " << i 
			    << " n'est pas dans le maillage " <<i0 << " " << i1 <<  endl;
		}
		
		//  find adj
		// reffecran();    
		
		for (i=0;i<nt;i++)
		{ 
		    Triangle & T=triangles[i];
		    for( j=0; j<3; j++,n++ )
		    { 
			VerticesNumberOfEdge(T,j,j0,j1);
			k = j0+j1; // code of current edge 
			int jm = Min(j0,j1), NbAdj=0, He=-1;
			int *pm=Head+k;
			while (*pm>=0) // be carefull  
			{                                 
			    int m=*pm,jj=m%3,ii=m/3, jj0,jj1;
			    VerticesNumberOfEdge(triangles[ii],jj,jj0,jj1);
			    if(jm==Min(jj0,jj1)) // same edge 
			    {                  
				NbAdj++;   
				// remove from  the liste 
				*pm=TheAdjacencesLink[m];
				TheAdjacencesLink[m]=He;  // link to He
				He = m;                  
			    }
			    else 
			    {
				NbCollision++;
				pm=TheAdjacencesLink+*pm; // next 
			    }
			}
			//  make a circular link 
			if (NbAdj>0) 
			{ 
			    int m=He;
			    while(TheAdjacencesLink[m]>=0)
				m=TheAdjacencesLink[m]; // find end of list     
							// close the List of common edges               
			    TheAdjacencesLink[m] = He;               
			}
			if (NbAdj >2) 
			{
			    int m=He;
			    do { 
				m=TheAdjacencesLink[m];
			    } while(TheAdjacencesLink[m]!=He);
			}
			
			if (NbAdj) NbOfEdges++;
			if(NbAdj==1)
                        {
                            if (! (TonBoundary[i]& MaskEdge[j]))
			    { NbOfMEdges++;
				if(verbosity>99) 
				    cout << " Edge (" << j0 << " "<< j1 << ") : "  << j  << " of Triangle " << &T-triangles << " on mortar \n"
				    <<" --- > " << number(T[0]) << " " << number(T[1]) << " " << number(T[2]) << " /" << int(TonBoundary[i])<< "\n" ;
				TonBoundary[i]+= AddMortar[j];
			    }
				else { NbOfBEdges++; }
                        }
		    }
		}
		    
		    if (verbosity>1 ) {
			cout << "    Nb of Vertices " << nv <<  " ,  Nb of Triangles " 
			<< nt << endl ;
			cout << "    Nb of edge on user boundary  " << neb 
			<< " ,  Nb of edges on true boundary  " << NbOfBEdges << endl;
			if(NbOfMEdges) cout << "    Nb of edges on Mortars  = " << NbOfMEdges << endl;

		    }
		    delete [] Head; // cleanning memory
		    NbMortars =0;
		    NbMortarsPaper=0;
	    }
		{
		    //  construct the mortar 
		    int * linkg = new int [nv]; //  to link
		    int * linkd = new int [nv]; //  to link
		    int * NextT3= new int[3*nt];
		    int * headT3= new int[nv];
		    ffassert( linkg && linkd);
		    for (int i=0;i<nv;i++) 
			headT3[i]=linkg[i]=linkd[i]=-1; // empty
							//   create the 2 link 
							//  reffecran();
							//  Draw(0);
		    for (int k=0;k<nt;k++)
			for (int j=0;j<3;j++)
			    if (TonBoundary[k] & AddMortar[j]) 
			    { 
				//   triangles[k].Draw(j,0.8);
				int s0,s1;
				VerticesNumberOfEdge(triangles[k],j,s0,s1);
				linkg[s0] = (linkg[s0] == -1) ? s1 : -2 ;
				linkd[s1] = (linkd[s1] == -1) ? s0 : -2 ; 
				//              throwassert(linkg[s0] != -1 && linkg[s1] != -1 );
				//   cout << "On Mortar " << s0 << " " << s1 << " link " <<  linkg[s0]  << " " << linkd[s1] <<endl  ;                              
			    } 
				//  we remove the  boundary link       
				for (int k=0;k<nt;k++)
				    for (int j=0;j<3;j++)
					if (TonBoundary[k] & MaskEdge[j]) 
					{ 
					    int s0,s1;
					    VerticesNumberOfEdge(triangles[k],j,s0,s1);
					    // cout << s0 << " " << s1 << " ld " << linkd[s0]  << " " << linkd[s1] << " lg " << linkg[s0]  << " " << linkg[s1] << " apres " ;              
					    linkg[s0] = linkg[s0] != -1 ?  -2 : -1;
					    linkg[s1] = linkg[s1] != -1 ?  -2 : -1;
					    
					    linkd[s1] = linkd[s1] != -1 ?  -2 : -1;                                                 
					    linkd[s0] = linkd[s0] != -1 ?  -2 : -1; 
					    // cout  << " ld " << linkd[s0]  << " " << linkd[s1] << " lg" << linkg[s0]  << " " << linkg[s1] << endl;
					    
					} 
					    
					    //    remark if   linkd[i]  == -2  extremities of mortars (more than 2 mortars)
					    //    if  ((linkd[i] != -1) || (linkd[i] != -1))   =>   i is on mortars 
					    //    make a link for all triangles contening a  mortars
					    const int k100=100;
		    SubMortar  bmortars[k100];
		    int k3j=0;
		    for (int k=0;k<nt;k++) {
			const  Triangle & T=triangles[k];
			for (int j=0;j<3;j++,k3j++) 
			{ 
			    NextT3[k3j]=-2; 
			    
			    int is= number(T[j]);
			    // if (linkg[is]<-1 || linkd[is]<-1) // extremite of bmortars
			    int j0 =EdgesVertexTriangle[j][0];  //  the 2 edges contening j 
			    int j1 =EdgesVertexTriangle[j][1];
			    if  (TonBoundary[k] & (AddMortar[j0]|AddMortar[j1]))  // (linkg[is]!=-1) 
			    { 
				throwassert(linkd[is]!=-1 || linkg[is]!=-1);
				NextT3[k3j]=headT3[is];
				headT3[is] =  k3j; 
			    }}
		    }
		    //    loop on extremite of bmortars 
		    // calcule du nombre de joint 
		    /*      rattente(1);      
		    for (int is=0;is<nv;is++)
			if ( (linkg[is]<-1 || linkd[is]<-1) )
		    {MoveTo( vertices[is]);
			    ostringstream ss;
			    ss << is ;
			    plotstring(ss.str().c_str());
		    }  */  
		    datamortars=0;
		    int kdmg = 0,kdmd=0; 
		    int kdmgaa = 0,kdmdaa=0; 
		    int step=0;
		    mortars=0;
		    NbMortars = 0;
		    NbMortarsPaper=0;
		    do { // two step  
			 // one to compute the NbMortars
			 // one to store in mortars and kdm 
			int * datag=0,*datad=0;
			ffassert(step++<2);
			if (NbMortars)
			{ // do allocation 
			    int kdm=kdmgaa+kdmdaa;
			    if (verbosity>2)
				cout << "   sizeof datamortars" << kdm << " NbMortars=" << NbMortars << endl;
			    mortars     = new Mortar [NbMortars];
			    datamortars = new int [kdm];
			    for (int i=0;i<kdm;i++)
				datamortars[i]=0;
			    datag=datamortars;
			    datad=datamortars+kdmgaa;
			    ffassert(kdmg && kdmd);
			}
			int onbm=NbMortars;
			// cout << "begin " << NbMortars << " g " << kdmg << " d " <<kdmd << endl;
			int kdmgo=kdmgaa;
			int kdmdo=kdmdaa;
			int kdm=kdmdaa+kdmgaa;
			//  reset 
			NbMortars =0;
			kdmg =0;
			kdmd =0;
			
			
			for (int is=0;is<nv;is++)
			    if (linkg[is] == -2) 
			    { // for all extremity of mortars 
				if(linkd[is] != -2)
				{
				  cout <<" Bug in mortar constrution : close to vertex "<< is << endl; 
				  ffassert(linkd[is] == -2);
				}
				const Vertex & S = vertices[is];
				R2  A(S);
				int km=0;
				int p;
				for ( p=headT3[is] ;p>=0; p=NextT3[p])
				{  //  for all nonconformitie around sm            
				    int k=p/3;
				    int i=p%3;
				    const Triangle & T(triangles[k]);
				    //throwassert( vertices + sm == &T[i]);
				    // for the 2 egdes contening the vertex j of the triangle k
				    
				    for (int jj=0;jj<2;jj++)   // 2 edge j contening i              
				    { 
					int j = EdgesVertexTriangle[i][jj];                     
					int s0,s1;
					VerticesNumberOfEdge(T,j,s0,s1);
					throwassert (s0  == is || s1 == is);
					if ( TonBoundary[k] & AddMortar[j])
					{
					    int ss,sens;
					    if ( s0 == is)
					    { ss=s1;sens=1;}
					    else
					    { ss=s0;sens=-1;}
					    const Vertex & vss( vertices[ss]);
					    bmortars[km++] = SubMortar(S,vss,k,i,sens);
					    throwassert(km<k100);
					} 
				    }
				}
				ffassert(p!=-2);
				
				// throwassert(km % 2 == 0);
				HeapSort(bmortars,km); 
				// cout <<" Nb of mortars around vertex "<< is << " = " << km<< endl;
				//  searh the same mortar (left and right)
				//  with same angle
				int im=km-1; // the last
				
				double eps = 1e-6;
				double thetai=bmortars[im].alpha-2*Pi;
				
				for (int jm=0;jm<km;im=jm++) 
				{ //  for all break of vertex 
				    
				    double theta= (bmortars[jm].alpha-thetai);
				    thetai=bmortars[jm].alpha;
				    if (theta < eps  || theta+eps > 2*Pi)
				    {//  same angle mod 2 * Pi
				     //  beginning of a mortars 
					if(mortars)  { //  set the pointeur
					    ffassert(NbMortars< onbm);
					    ffassert(datag && datad);
					    ffassert(datag< datamortars+ kdm);
					    mortars[NbMortars].left  = datag;
					    mortars[NbMortars].right = datad;
					    mortars[NbMortars].Th = this;
					}
					// cout << "   Begin mortar " << is << " " << im << " " << jm << " theta =" << thetai << endl; 
					R2 AM(A,bmortars[im].to);
					AM = AM/Norme2(AM);
					int ig(im),id(jm);
					if ( bmortars[im].sens <0) ig=jm,id=im;
					//SubMortar & mg(bmortars[ig]);
					//SubMortar & md(bmortars[id]);
					//  loop sur droite gauche
					//  meme sommet
					//   0 = gauche 1=droite 
					int sgd[2]={0,0};
					int *link[2];                   
					int nm[2]={0,0}; 
					R  ll[2]={0.0,0.0}; 
					// 
					sgd[0]=sgd[1]=is; 
					link[0] = linkg;
					link[1] = linkd;
					int gd=0; //  gd = 0 => left side  an gd=1 => right side
					
					
					int kkkk=0;
					do { //   for all 
					    
					    
					    int sm = sgd[gd]; 
					    int  dg = 1-gd;
					    // cout << " sm = " << sm << " h=" << headT3[sm] << " gb=" << gd << " autre " << sgd[dg ] << " " << kkkk << endl;
					    
					    R lAV,avam;
					    Vertex *pV=0;
					    int p,k,i,j;
					    //  search the the first start ( sens = gd )
					    throwassert(headT3[sm]>=0);// on a mortar ?? 
						for ( p=headT3[sm] ;p>=0; p=NextT3[p])
						{                
						    k=p/3;
						    i=p%3;
						    const Triangle & T(triangles[k]);
						    // couleur(1);T.Draw(0.3);
						    throwassert( vertices + sm == &T[i]);
						    // for the 2 egdes contening the vertex j of the triangle k
						    j = EdgesVertexTriangle[i][dg];
						    Vertex &V = T[VerticesOfTriangularEdge[j][dg]];
						    throwassert( &T[VerticesOfTriangularEdge[j][gd]] == vertices + sm);                
						    if ( TonBoundary[k] & AddMortar[j])
						    {  // check the sens and the direction
							
							//  cout << number(T[VerticesOfTriangularEdge[j][gd]])  << " pv=" 
							//      << number(V)  << " sm=" << sm << " " << headT3[number(V)] << endl; 
							
							R2 AV(A,V);
							lAV = Norme2(AV);
							avam = (AV,AM);
							// cout << " ---  " << avam << " > " << ll[gd] <<  " " << Abs((AM.perp(),AV) )
							//     << " " << ( avam > ll[gd] && Abs((AM.perp(),AV)) < lAV * 1e-6 ) << endl;
							// go ahead in direction AM 
							if ( (avam > ll[gd])  && Abs((AM.perp(),AV)) < lAV * 1e-6 )  
							{pV = &V;break;} //  ok good                         
						    }
						}
						if ( ! (p>=0 && pV))
						    throwassert(p>=0 && pV); //  PB reach the end without founding 
						if ( ! ( Abs((AM.perp(),A-*pV)) < 1e-5) )
						    // cout << Abs((AM.perp(),A-*pV)) <<*pV << endl, 
						    throwassert( Abs((AM.perp(),A-*pV)) < 1e-5);
						
						throwassert(sm != number(pV));
						int kkgd= 3*k + j;   
						ll[gd] = avam;
						//   find the SubMortar m of vertex  sm
						//   with same sens and direction 
						//   
						
						for (int s= number(pV);s>=0;s=link[gd][s],nm[gd]++) 
						{
						    //  on est sur l'arete kkgd
						    
						    throwassert( s == number(triangles[kkgd/3][VerticesOfTriangularEdge[kkgd%3][dg]]));
						    sgd[gd]=s;// save the last
							if ( ! ( Abs((AM.perp(),A-vertices[s])) < 1e-5) )
							    // cout << Abs((AM.perp(),A-vertices[s])) << vertices[s] << endl, 
							    throwassert( Abs((AM.perp(),A-vertices[s])) < 1e-5);
							//cout << " s=" << s << " h=" << headT3[s] << " " << link[gd][s] << " " <<  link[dg][s] << endl; ;
							throwassert(kkgd>=0 && kkgd < 3*nt);
							if (datamortars) 
							{
							    throwassert(datag - datamortars == nm[0] + kdmg);
							    throwassert(datad - datamortars == nm[1] + kdmd + kdmgo );
							    
							    if (gd == 0)  *datag++ = kkgd; // store 
							    else          *datad++ = kkgd; //
							    
							}  
							
							//                            cout  << " ++++ "<<  ll[gd] << " > " << ll[dg]  << " " << headT3[sgd[dg]] << " " <<sgd[dg] << endl;    
							
							int kk=kkgd,kkk=0,kkkk=0;
							if ( link[gd][s] >=0) { 
							    for (int pp=headT3[s] ;pp>=0; pp=NextT3[pp],kkk++)
							    {  throwassert(number(triangles[pp/3][pp%3]) == s);
								
								if( (pp/3)!= (kk/3))
								{
								    kkkk++,kkgd=pp - (pp%3) + EdgesVertexTriangle[pp%3][dg];                                  
								}
							    }
							    throwassert( kkgd/3 != kk/3);
							    throwassert(kkk==2 && kkkk==1);}  
						} 
						
						if (ll[gd]>ll[dg] &&  headT3[sgd[dg]]>=0) //changement de cote 
						    gd = dg; 
						throwassert(kkkk++<100);
						//cout <<" kkkk=" << kkkk << " " << sgd[0] << " " << sgd[1] << endl;
					} while (sgd[0] != sgd[1]); 
					
					kdmgaa = Max(kdmgaa,kdmg  + nm[0]);
					kdmdaa = Max(kdmdaa,kdmd  + nm[1]);                                                                    
					
					if (is < sgd[0]  &&  headT3[sgd[0]] >=0) {
					    //cout << "    Mortars from (on saute) " << is << " to " << sgd[0] << " " << nm[0] << " " << nm[1]<< " " <<  kdmg <<  " " << kdmd << endl;
					    if( mortars ) { //  restore 
						datag -= nm[0];
						datad -= nm[1];   } 
					    
					} 
					else {
					    // cout << "    Mortars from " << is << " to " << sgd[0] << " " << nm[0] << " " << nm[1]<< " " <<  kdmg <<  " " << kdmd << endl;
					    
					    if(mortars ) {
						ffassert(NbMortars< onbm);
						mortars[NbMortars].nleft  = nm[0];
						mortars[NbMortars].nright = nm[1];
						
						//  check
						for (int i=0;i<mortars[NbMortars].nleft;++i)
						    if ( mortars[NbMortars].left[i] <0 ||  mortars[NbMortars].left[i] < 3*nt)
							ffassert(0);
						for (int i=0;i<mortars[NbMortars].nright;++i)
						    if ( mortars[NbMortars].right[i] <0 ||  mortars[NbMortars].right[i] < 3*nt)
							ffassert(0);                             
						ffassert(datag <= datamortars + kdmgo + kdmdo); 
						ffassert(datad <= datamortars + kdmgo + kdmdo); 
					    }
					    kdmg += nm[0];
					    kdmd += nm[1]; 
					    NbMortars++;
					}
				    } // same angle 
				}//  for all break of vertex                  
			    } // for all extremity of mortars 
				if (verbosity>1 && NbMortars) 
				    cout << "    Nb Mortars " << NbMortars << /*" " << kdmg << " "<<  kdmd <<*/ endl;
			if (mortars) 
			{
			    // cout << "end " << NbMortars << " g " << kdmg << " d " <<kdmd << endl;
			    // cout << kdmgo << " " << kdmg << " " << kdmdo << " " << kdmd << endl;
			    ffassert(kdmgo == kdmgaa && kdmdo == kdmdaa);}
		    }  while (NbMortars && !mortars) ;
		    
		    //   rattente(1);
		    //   reffecran();
		    //    compute the NbMortarsPaper
		    int t3,nt3 = nt*3;
		    NbMortarsPaper=0;
		    for (int i=0;i<nt3;i++) 
			NextT3[i]=0;
		    for (int i=0;i<NbMortars;i++)
		    {
			
			if (mortars[i].nleft==1)
			{ t3=mortars[i].left[0];}
			else  if (mortars[i].nright==1)
			{ t3=mortars[i].right[0];}
			else { cout << "   -- Bizarre " << mortars[i].nleft << " " << mortars[i].nright << endl;
			}
			if (NextT3[t3]==0) NbMortarsPaper++;
			NextT3[t3]++;
		    }
		    delete [] linkg;
		    delete [] linkd;
		    delete [] NextT3;
		    delete [] headT3;
		    
		}//  end of mortar construction   
		
		
		
		delete [] TonBoundary;
		
		//  construction of TriangleConteningVertex
		TriangleConteningVertex = new int[nv];
		for (int it=0;it<nt;it++)
		    for (int j=0;j<3;j++)
			TriangleConteningVertex[(*this)(it,j)]=it;
		
		Buildbnormalv();
		if (verbosity>4) 
		{ 
		    cout << "    Number of Edges                 = " << NbOfEdges << endl;
		    cout << "    Number of Boundary Edges        = " << NbOfBEdges << " neb = " << neb << endl;
		    cout << "    Number of Mortars  Edges        = " << NbOfMEdges << endl;
		    cout << "    Nb Of Mortars with Paper Def    = " <<  NbMortarsPaper << " Nb Of Mortars = " << NbMortars;             
		    cout << "    Euler Number nt- NbOfEdges + nv = " 
			<< nt + NbMortars - NbOfEdges + nv << "= Nb of Connected Componant - Nb Of Hole " 
			<<endl;}
		
    }
	    void  Mesh::BoundingBox(R2 &Pmin,R2 &Pmax) const 
	    {
		ffassert(nv);
		Pmin=Pmax=vertices[0];
		for (int i=0;i<nv;i++)
		{ 
		    const R2 & P=vertices[i];
		    Pmin.x = Min(Pmin.x,P.x);
		    Pmin.y = Min(Pmin.y,P.y);
		    Pmax.x = Max(Pmax.x,P.x);
		    Pmax.y = Max(Pmax.y,P.y);
		} 
		//  cout << " Bounding Box = " << Pmin <<" " <<  Pmax << endl;       
	    }
    void Mesh::read(const char * filename)
    {
	ifstream f(filename);
	if (!f) {
	    cerr << "Erreur ouverture du fichier " << filename << endl;
	throw(ErrorExec("exit",1));}
	// ffassert(f);
	if(verbosity)
	    cout << "  -- Mesh::read On file \"" <<filename<<"\""<<  endl;
	read(f);
    }
    void Mesh::read(ifstream & f , bool bin )
    { // read the mesh
	dim=2;
	ne=0;
	ntet=0;
	volume=0;
	TriangleConteningVertex =0;
	BoundaryAdjacencesHead=0;
	BoundaryAdjacencesLink=0;
	BoundaryEdgeHeadLink=0;
	quadtree =0;
	NbMortars=0;
	tet=0;
	edges=0;    
	mortars=0;
	TheAdjacencesLink =0;
	area=0;
	bnormalv=0;
	
	int i,i0,i1,i2,ir;
	
	
	f >> nv >> nt >> neb ;
	if(verbosity>10)
	    cout << "    Nb of Vertex " << nv << " " << " Nb of Triangles " 
	    << nt << " Nb of boundary edge " << neb <<  endl;
	ffassert(f.good() && nt && nv) ;
	triangles = new Triangle [nt];
	vertices  = new Vertex[nv];
	bedges    = new BoundaryEdge[neb];
	area=0;
	ffassert(triangles && vertices && bedges);
	
	for (i=0;i<nv;i++)    
	    f >> vertices[i],ffassert(f.good());
	
	for (i=0;i<nt;i++) { 
	    f >> i0 >> i1 >> i2 >> ir;
	    ffassert(f.good() && i0>0 && i0<=nv && i1>0 && i1<=nv && i2>0 && i2<=nv);
	    triangles[i].set(vertices,i0-1,i1-1,i2-1,ir); 
	area += triangles[i].area;}
	
	for (i=0;i<neb;i++) { 
	    f >> i0 >> i1 >> ir;
	bedges[i] = BoundaryEdge(vertices,i0-1,i1-1,ir); }
	
	if(verbosity)
	    cout << "   End of read: area on mesh = " << area <<endl;  
	ConsAdjacence();
	//   BoundingBox(cMin,cMax);//  Set of cMin,Cmax
    }
    
	    Mesh::Mesh(int nbv,int nbt,int nbeb,Vertex *v,Triangle *t,BoundaryEdge  *b)
	    {
		TriangleConteningVertex =0;
		BoundaryAdjacencesHead=0;
		BoundaryAdjacencesLink=0;
		BoundaryEdgeHeadLink=0;
		quadtree =0;
		NbMortars=0;
		mortars=0;
		dim=2;
		tet=0;
		volume=0;
		edges=0;
		ntet=0;
		ne=0;    
		TheAdjacencesLink =0;
		nv=nbv;
		nt =nbt;
		neb=nbeb;
		triangles=t;
		vertices=v;
		bedges=b;
		area=0;
		bnormalv=0;
		if (t && nt >0) 
		{
		    for (int i=0;i<nt;i++)  
			area += triangles[i].area;
		    ConsAdjacence();  
		}  
		else
		{
		    bool removeouside=nbt>=0;
		    nt=0;
		    BuilTriangles(true,removeouside);
		    /*    if( ! removeouside)
		    {
			neb=0; // remove the edge
			delete [] bedges;
			bedges=0;        
			
		    } */
		    ConsAdjacence();
		}
		// cout << " build: Mesh : " << this << endl;
	    }  
	    
	    inline int BinaryRand(){
#ifdef RAND_MAX
		const long HalfRandMax = RAND_MAX/2;
		return rand() <HalfRandMax;
#else
		return rand() & 16384; // 2^14 (
#endif
		
} 
int  WalkInTriangle(const Mesh & Th,int it, double *lambda,
                    const  KN_<R> & U,const  KN_<R> & V, R & dt)
{
    const Triangle & T(Th[it]);
    const R2 Q[3]={(const R2) T[0],(const R2) T[1],(const R2) T[2]};
    int i0=Th.number(T[0]);
    int i1=Th.number(T[1]);
    int i2=Th.number(T[2]);
    
    R u   = lambda[0]*U[i0] + lambda[1]*U[i1] + lambda[2]*U[i2];
    R v   = lambda[0]*V[i0] + lambda[1]*V[i1] + lambda[2]*V[i2];
    R2 P  = lambda[0]*Q[0]  + lambda[1]*Q[1]  + lambda[2]*Q[2];
    
    //  cout << " " << u << " " << v ;
    R2 PF = P + R2(u,v)*dt;
    
    //  couleur(15);MoveTo( P); LineTo( PF);
    R l[3];
    l[0] = Area2(PF  ,Q[1],Q[2]);
    l[1] = Area2(Q[0],PF  ,Q[2]);
    l[2] = Area2(Q[0],Q[1],PF  );
    R Det = l[0]+l[1]+l[2];
    l[0] /= Det;
    l[1] /= Det;
    l[2] /= Det;
    const R eps = 1e-5;
    int neg[3],k=0;
    int kk=-1;
    if (l[0]>-eps && l[1]>-eps && l[2]>-eps) 
    {
	dt =0;
	lambda[0] = l[0];
	lambda[1] = l[1];
	lambda[2] = l[2];
    }
    else 
    {
	
	if (l[0]<eps && lambda[0] != l[0]) neg[k++]=0;
	if (l[1]<eps && lambda[1] != l[1]) neg[k++]=1;
	if (l[2]<eps && lambda[2] != l[2]) neg[k++]=2;
	R eps1 = T.area     * 1.e-5;
	
	if (k==2) // 2  
	{
	    // let j be the vertex beetween the 2 edges 
	    int j = 3-neg[0]-neg[1];
	    // 
	    R S = Area2(P,PF,Q[j]);
	    
	    if (S>eps1)
		kk = (j+1)%3;
	    else if (S<-eps1)
		kk = (j+2)%3;
	    else if (BinaryRand())
		kk = (j+1)%3;
	    else
		kk = (j+2)%3;
	    
	} 
	else if (k==1)
	    kk = neg[0];
	if(kk>=0)
	{
	    R d=lambda[kk]-l[kk];
	    
	    throwassert(d);
	    R coef =  lambda[kk]/d;
	    R coef1 = 1-coef;
	    dt        = dt*coef1;
	    lambda[0] = lambda[0]*coef1 + coef *l[0];
	    lambda[1] = lambda[1]*coef1 + coef *l[1];
	    lambda[2] = lambda[2]*coef1 + coef *l[2];
	    lambda[kk] =0;
	}
    }
    int jj=0;
    R lmx=lambda[0];
    if (lmx<lambda[1])  jj=1, lmx=lambda[1];
    if (lmx<lambda[2])  jj=2, lmx=lambda[2];
    if(lambda[0]<0) lambda[jj] += lambda[0],lambda[0]=0;
    if(lambda[1]<0) lambda[jj] += lambda[1],lambda[1]=0;
    if(lambda[2]<0) lambda[jj] += lambda[2],lambda[2]=0;
    return kk;
}
R2  SubInternalVertex(int N,int k)
    {
	if(N<0)
	  {
	      R eps = 1e-08;
	      R2 p[3][3]= { 
		  { R2(1-eps,+eps),R2(eps,1-eps),R2(1./3.+eps,1./3.+eps) },
		  { R2(0,1-eps)  ,R2(0,+eps),R2(1./3.-eps,1./3.) },
		  { R2(eps,0),R2(1-eps,0),R2(1./3.,1./3.-eps) } };
	      
	      int j=k%3;
	      R2 P=SubInternalVertex(-N,k/3);
	      R l0=1.-P.x-P.y,l1=P.x,l2=P.y;
	      return p[j][0]*l0+ p[j][1]*l1+ p[j][2]*l2;
	  }
	int i,j;
	num1SubTVertex(N,k,i,j);
	return R2( (double) i/ (double)N,(double) j/(double)N);
    }
    
R2 SubTriangle(const int N,const int n,const int l)
{
    // compute the subdivision of a triangle in N*N
    // N number of sub division
    // n number of the sub triangle
    // l vertex of the sub triangle
    if(N<0)
    {
	R eps = 1e-08;
	R2 p[3][3]= { 
	    { R2(1-eps,+eps),R2(eps,1-eps),R2(1./3.+eps,1./3.+eps) },
	    { R2(0,1-eps)  ,R2(0,+eps),R2(1./3.-eps,1./3.) },
	    { R2(eps,0),R2(1-eps,0),R2(1./3.,1./3.-eps) } };
	    
	int j=n%3;
	R2 P=SubTriangle(-N,n/3,l);
	R l0=1.-P.x-P.y,l1=P.x,l2=P.y;
	return p[j][0]*l0+ p[j][1]*l1+ p[j][2]*l2;
    }
    throwassert(n < N*N);
    int i = n % N;
    int j = n / N;
    int k = N - i - j;
    if(k<=0)
      {
	if(l==1) i++;
	else if(l==2) j++;
      }
    else
      if(l==1) j++;
      else if(l==2) i++;
    // if(l==1) i++;
    // if(l==2) j++;
    // if ( k <= 0 )cout << " - " << endl;
    return k >0 
	? R2( (double) i/ (double)N,(double) j/(double)N)
	: R2( (double) (N-j)/ (double)N,(double) (N-i)/(double)N);
    
} 

int numSubTriangle(const int N,const int n,const int l)
    {
	// compute the subdivision of a triangle in N*N
	// N number of sub division
	// n number of the sub triangle
	// l vertex of the sub triangle
	if(N<0)
	  {
	    int j=n%3;
	    return numSubTriangle(-N,n/3,l)*3+j;
	  }
	throwassert(n < N*N);
	int i = n % N;
	int j = n / N;
	int k = N - i - j;
	if(k<=0)
	  {
	    if(l==1) i++;
	    else if(l==2) j++;
	  }
	else
	  if(l==1) j++;
	  else if(l==2) i++;

	// if ( k <= 0 )cout << " - " << endl;
	return k >0 
	? numSubTVertex(N,i, j)
	: numSubTVertex(N,N-j,N-i);
	
} 
    

int Walk(const Mesh & Th,int& it, R *l,
         const KN_<R>  & U,const KN_<R>  & V, R dt) 
{
    
    int k=0;
    int j; 
    while ( (j=WalkInTriangle(Th,it,l,U,V,dt))>=0) 
    { 
	//int jj  = j;
	throwassert( l[j] == 0);
	R a= l[(j+1)%3], b= l[(j+2)%3];
	int itt =  Th.ElementAdj(it,j);
	if(itt==it || itt <0)  return -1;
	it = itt;
	l[j]=0;
	l[(j+1)%3] = b;
	l[(j+2)%3] = a;
	ffassert(k++<1000);
    }
    return it;
}
const Triangle *  Mesh::Find( R2 P, R2 & Phat,bool & outside,const Triangle * tstart) const
{
    int CasePichon=0;
    int it,j,it00;
    const Triangle *  rett=0;
    if ( tstart )
	it =  (*this)(tstart);
    else  
    {  
	const Vertex * v=quadtree->NearestVertexWithNormal(P);
	if (!v) 
	{ 
	    v=quadtree->NearestVertex(P);
	    assert(v);
	}
	/*
	 if (verbosity>100) 
	 cout << endl << (*this)(v) << *v << " " << Norme2(P-*v) << endl; 
	 */
	it00=it=Contening(v);
    }
RESTART:
    //     int itdeb=it;
    //     int count=0;
    //     L1: 
    outside=true; 
    //int its=it;
    int iib=-1;//,iit=-1;
	R dP=DBL_MAX;
	R2 PPhat;
	const Triangle * tt;
	int k=0,kout=0;    
	kfind++;
	while (1)
	{ 
	    const Triangle & K(triangles[it]);
	    kthrough++;
	    if (k++>=10000) 
	    {
		/* cout << P << endl;
		reffecran();
		Draw(0);
		triangles[its].Fill(2);
		DrawMark(P,0.01);
		rattente(1);  */
		ffassert(k++<10000);
	    }
	    int kk,n=0,nl[3];
	    
	    
	    
	    R2 & A(K[0]), & B(K[1]), & C(K[2]);
	    R l[3]={0,0,0};
	    R area2= K.area*2;
	    R eps =  -area2*1e-6;
	    l[0] = Area2(P,B,C);
	    l[1] = Area2(A,P,C);
	    l[2] = area2-l[0]-l[1];
	    if (l[0] < eps) nl[n++]=0;
	    if (l[1] < eps) nl[n++]=1;
	    if (l[2] < eps) nl[n++]=2;
	    
	    if (n==0)
	    {  // interior => return
		outside=false; 
		Phat=R2(l[1]/area2,l[2]/area2);
		return &K;
	    }
	    else if (n==1) 
		kk=0;
	    else  
		kk=BinaryRand() ? 1 : 0;
	    
	    j= nl[ kk ];
	    
	    int itt =  ElementAdj(it,j);
	    if(itt!=it && itt >=0)  
	    {
		dP=DBL_MAX;
		it=itt;
		continue;
	    }  
	    
	    //  edge j on border
	    l[j]=0;
	    
	    if ( n==2 ) 
	    {
		kk = 1-kk;
		j= nl[ kk ];
		itt =  ElementAdj(it,j);                  
		if (itt >=0 && itt != it) // correction FH oct 2009
		{
		    dP=DBL_MAX;
		    it=itt;
		    continue;
		}
		//  on a corner of the mesh 
		l[0]=l[1]=l[2]=0.;
		l[3-nl[0]-nl[1]]=1.; // correction april 2015 FH
		Phat=R2(l[1],l[2]);
	        rett=triangles +it;
	        if(searchMethod && outside) goto PICHON;
	      return rett;
	    }
	    //   on the border 
	    //   projection Ortho
	    
	    kout++;
	    
	    int j0=(j+1)%3,i0= &K[j0]-vertices;
	    int j1=(j+2)%3,i1= &K[j1]-vertices;
	    int ii,jj,iii;
	    bool ret=false;
	    
	    R2 AB=R2(K[j0],K[j1]),  AP(K[j0],P), BP(K[j1],P);
	    R la=  (AB,AP);
	    R lb= -(AB,BP);
	    if(la<0)     
		ii= i0, jj = j0,iii=i1;
	    else if ( lb <0)
		ii= i1, jj = j1,iii=i0;
	    else // PROJECTION between A,B
		ret = true;
	    if( ! ret)
	    { //  VERIF THE DISTANCE**2 Dicrease
		R2 Pjj(P,K[jj]);
		R dd = (Pjj,Pjj);
		if (dd >= dP ) {
		    Phat=PPhat;
		    // if(kout>1) cout << "        @ " << tt-triangles  << " " << Phat << " " << outside << endl;
		    rett=tt;
		    if(searchMethod && outside) goto PICHON;
		    return tt;
		}
		else
		{ 
		    l[0]=l[1]=l[2]=0;
		    l[jj]=1;
		    PPhat.x=l[1];
		    PPhat.y=l[2];                
		    dP=dd;
		    tt = triangles + it ;
		}
	    }  
	    /*
	     if(kout>1)
	     cout << "Find --- "<< P << " :  k=" << k << "  la lb " << la << " " << lb 
	     << " ii =" << ii << " iii= " << iii << " " << it << " dP= " 
	     << dP << " ret = " << ret <<endl;
	     */
	    if (ret || ii == iib) 
	    {  
		l[j]=0;
		
		l[j0]= Max(0.,Min(+lb/(la+lb),1.));
		l[j1]= 1-l[j0];
		Phat=R2(l[1],l[2]);
		//if(kout>1) cout << "        # " << it << " " << Phat << " " << outside << endl; 
		rett=triangles +it;
		if(searchMethod && outside) goto PICHON;
		return rett;
	    }
	    bool ok=false;
	    // next edge on true boundaryS
	    for (int p=BoundaryAdjacencesHead[ii];p>=0;p=BoundaryAdjacencesLink[p])
	    { int e=p/2, ie=p%2;// je=2-ie;
				// cout << number(bedges[e][0]) << " " << number(bedges[e][1]) << endl;
		if (!bedges[e].in( vertices+iii)) //  edge not equal  to i0 i1 
		{  
		    ok=true;
		    iib = ii;
		    it= BoundaryElement(e,ie);  //  next triangle                      
						 // cout << "  ------ " << it << " " << Phatt <<  endl;
		    break;
		}
	    }
	    ffassert(ok); 
	}
PICHON:	// Add dec 2010 ...
    CasePichon++;
    if(CasePichon==1) {// hack feb 2016  ..... ???????
        // change starting triangle  ????
        R2 PF=(*rett)(Phat);
        R2 PP=P + (P-PF);
        const Vertex * v=quadtree->NearestVertexWithNormal(PP);
        if (!v)
        {
            v=quadtree->NearestVertex(P);
            assert(v);
        }
        /*
         if (verbosity>100)
         cout << endl << (*this)(v) << *v << " " << Norme2(P-*v) << endl;
         */
        it=Contening(v);
        if( it != it00) goto RESTART;
    }
    
    npichon2d++;
	// Brute force .... bof bof ...
    double ddp=1e100;
    int pk=-1;
  
    for(int k=0;k<nt;++k)
      {
	int n=0,nl[3];
	Triangle & K=triangles[k];
	R2 & A(K[0]), & B(K[1]), & C(K[2]), G((A+B+C)/3.);
	R l[3]={0,0,0};
	R area2= K.area*2;
	R eps =  -area2*1e-6;
	l[0] = Area2(P,B,C);
	l[1] = Area2(A,P,C);
	l[2] = area2-l[0]-l[1];
	if (l[0] < eps) nl[n++]=0;
	if (l[1] < eps) nl[n++]=1;
	if (l[2] < eps) nl[n++]=2;
	if (n==0)
	  {  // interior => return
	      outside=false; 
	      Phat=R2(l[1]/area2,l[2]/area2);
	      return &K;
	  }
	R2 GP(G,P);
	double lgp2=(GP,GP);
	if(ddp > lgp2) {
	    ddp=lgp2;
	    pk=k;
	    
	}
      }
    
    return rett;
}


int  WalkInTriangle(const Mesh & Th,int it, double *lambda,
		    R u, R v, R & dt)
{
    const Triangle & T(Th[it]);
    const R2 Q[3]={(const R2) T[0],(const R2) T[1],(const R2) T[2]};
    // int i0=Th.number(T[0]);
    // int i1=Th.number(T[1]);
    // int i2=Th.number(T[2]);
    
    R2 P  = lambda[0]*Q[0]  + lambda[1]*Q[1]  + lambda[2]*Q[2];
    
    //  cout << " " << u << " " << v ;
    R2 PF = P + R2(u,v)*dt;
    
    //  couleur(15);MoveTo( P); LineTo( PF);
    R l[3];
    l[0] = Area2(PF  ,Q[1],Q[2]);
    l[1] = Area2(Q[0],PF  ,Q[2]);
    l[2] = Area2(Q[0],Q[1],PF  );
    R Det = l[0]+l[1]+l[2];
    l[0] /= Det;
    l[1] /= Det;
    l[2] /= Det;
    const R eps = 1e-5;
    int neg[3],k=0;
    int kk=-1;
    if (l[0]>-eps && l[1]>-eps && l[2]>-eps) 
    {
	dt =0;
	lambda[0] = l[0];
	lambda[1] = l[1];
	lambda[2] = l[2];
    }
    else 
    {
	
	if (l[0]<eps && lambda[0] != l[0]) neg[k++]=0;
	if (l[1]<eps && lambda[1] != l[1]) neg[k++]=1;
	if (l[2]<eps && lambda[2] != l[2]) neg[k++]=2;
	R eps1 = T.area     * 1.e-5;
	
	if (k==2) // 2  
	{
	    // let j be the vertex beetween the 2 edges 
	    int j = 3-neg[0]-neg[1];
	    // 
	    R S = Area2(P,PF,Q[j]);
	    
	    if (S>eps1)
		kk = (j+1)%3;
	    else if (S<-eps1)
		kk = (j+2)%3;
	    else if (BinaryRand())
		kk = (j+1)%3;
	    else
		kk = (j+2)%3;
	    
	} 
	else if (k==1)
	    kk = neg[0];
	if(kk>=0)
	{
	    R d=lambda[kk]-l[kk];
	    
	    throwassert(d);
	    R coef =  lambda[kk]/d;
	    R coef1 = 1-coef;
	    dt        = dt*coef1;
	    lambda[0] = lambda[0]*coef1 + coef *l[0];
	    lambda[1] = lambda[1]*coef1 + coef *l[1];
	    lambda[2] = lambda[2]*coef1 + coef *l[2];
	    lambda[kk] =0;
	}
    }
    int jj=0;
    R lmx=lambda[0];
    if (lmx<lambda[1])  jj=1, lmx=lambda[1];
    if (lmx<lambda[2])  jj=2, lmx=lambda[2];
    if(lambda[0]<0) lambda[jj] += lambda[0],lambda[0]=0;
    if(lambda[1]<0) lambda[jj] += lambda[1],lambda[1]=0;
    if(lambda[2]<0) lambda[jj] += lambda[2],lambda[2]=0;
    return kk;
}        
Mesh::~Mesh()
{
    //(cout << "   -- delete mesh " << this << endl);
    delete  quadtree;
    delete [] triangles;
    delete [] vertices;
    delete [] bedges;
    delete [] mortars;
    delete [] datamortars;
    // delete [] adj;
    delete [] TheAdjacencesLink;
    delete [] BoundaryEdgeHeadLink;
    delete [] BoundaryAdjacencesHead;
    delete [] BoundaryAdjacencesLink;
    delete []  TriangleConteningVertex;
    delete [] bnormalv;
    delete [] tet;
    delete [] edges;
}
//  for the  mortar elements
 int NbOfSubTriangle(int k)
{  
    if(k>0) return  k*k;
    else if(k<0) return 3*(k*k);
    ffassert(0);
    return 0;
}
    

 int NbOfSubInternalVertices(int kk)
{ 
    assert(kk);
    int k=Abs(kk);
    int  r= (k+2)*(k+1)/2;
    assert(r>=0);
    
    return kk<0 ? 3*r : r;
}



Mesh::Mesh(int nbv,R2 * P)
{
    
    TheAdjacencesLink=0;
    BoundaryEdgeHeadLink=0;
    BoundaryAdjacencesHead=0;
    BoundaryAdjacencesLink=0; 
    TriangleConteningVertex=0;              
    TriangleConteningVertex=0;
    dim=2;
    tet=0;
    edges=0;
    ntet=0;
    ne=0;
    volume=0;
    
    quadtree =0;
    NbMortars=0;
    NbMortarsPaper=0;
    nt=0;
    mortars=0;
    TheAdjacencesLink =0;
    datamortars=0;
    quadtree =0;
    bedges=0;
    neb =0;
    bnormalv=0;
    nv=nbv;
    vertices  = new Vertex[nv];
    triangles=0;
    area=0;
    for (int i=0;i<nv;i++)
	vertices[i]=P[i];
    bedges    =  0; 
    area=0;
    
    
    BuilTriangles(false) ;  
    
    ConsAdjacence();    
    
}
void Mesh::BuilTriangles(bool empty,bool removeouside)
{
    long nba = neb;
    // 
    long nbsd = 0; // bofbof 
    if(!removeouside) nbsd=1;
    // faux just pour un test 
    long  *sd;
    sd=new long[2];
    sd[0]=-1;
    sd[1]=0;        
    nbsd=removeouside?0:-1;
    
    long nbs=nv;
    long nbsmax=nv;
    long            err = 0;//, nbsold = nbs;       
	long           *c = 0;
	long           *tri = 0;
	long           *nu = 0;
	long           *reft = 0;
	typedef double Rmesh;
	Rmesh          *cr = 0;
	Rmesh          *h = 0;
	long nbtmax = 2 * nbsmax;
	long * arete  = nba ? new long[2*nba] : 0; 
	nu = new long[6*nbtmax];
	c = new long[2*nbsmax];
	tri = new long[(4 * nbsmax + 2 * nbsd)];
	reft = new long[nbtmax];
	cr = new Rmesh[(2 * nbsmax + 2)];
	h = new Rmesh[nbsmax];
	for(int i=0,k=0; i<nv; i++)
	{ 
	    cr[k++]  =vertices[i].x;
	    cr[k++]=vertices[i].y;
	    h[i]=1; 
	}
	for (int i=0,k=0 ;i<neb;i++)
	{
	    arete[k++] =number(bedges[i][0])+1;
	    arete[k++] =number(bedges[i][1])+1;
	}
	
	
	extern int 
	    mshptg8_ (Rmesh *cr, Rmesh *h, long *c, long *nu, long *nbs, long nbsmx, long *tri,
		      long *arete, long nba, long *sd,
		      long nbsd, long *reft, long *nbt, Rmesh coef, Rmesh puis, long *err);
	
	long nbt=0;
	mshptg8_ (cr, h, c, nu, &nbs, nbs, tri, arete, nba, (long *) sd, nbsd, reft, &nbt, .25, .75, &err);
	if(err) {
	    cerr << " Sorry Error build delaunay triangle   error = " << err << endl;
	    delete [] arete;
	    delete [] nu;
	    delete [] c;
	    delete [] tri;
	    delete [] reft;
	    delete [] cr;
	    delete [] h;
	    delete [] sd;	   
	    throw(ErrorExec("Error mshptg8_",1));
	   	}
	assert(err==0 && nbt !=0);
	delete [] triangles;
	nt = nbt;
	if(verbosity>1)
	    
	    cout << " Nb Triangles = " << nbt << endl;
	triangles = new Triangle[nt];
	for(int i=0,k=0;i<nt;i++,k+=3)
        {
	  if( verbosity>1000)
	    {
	      cout << vertices[nu[k]-1]<< "  " << endl;
	      cout << vertices[nu[k+1]-1]<< "  " << endl;
	      cout << vertices[nu[k+2]-1]<< "  " << endl;
	      cout << vertices[nu[k]-1]<< "  \n\n" << endl;
	    }
	    triangles[i].set(vertices,nu[k]-1,nu[k+1]-1,nu[k+2]-1,reft[i]);
	    area += triangles[i].area;
        }  
	
	delete [] arete;
	delete [] nu;
	delete [] c;
	delete [] tri;
	delete [] reft;
	delete [] cr;
	delete [] h;
	delete [] sd;
}
Mesh::Mesh(const Mesh & Th,int * split,bool WithMortar,int label)
{ //  routine complique 
  //  count the number of elements
    area=0; //Th.area;
    lenbord=0;
    volume=0;
    BoundaryAdjacencesHead=0;
    BoundaryAdjacencesLink=0;
    BoundaryEdgeHeadLink=0;
    quadtree =0;
    NbMortars=0;
    ntet=0;
    ne=0;
    dim=2;
    tet=0;
    edges=0;
    mortars=0;
    TheAdjacencesLink =0;
    quadtree =0;
    bedges=0;
    neb =0;
    bnormalv=0;
    R2 Pmin,Pmax;
    Th.BoundingBox(Pmin,Pmax);
    nt=0;
    int nebi=0; // nb arete interne 
    int nbsdd =0;
    int splitmin=100,splitmax=0;
    for (int i=0;i<Th.nt;i++)
	if(split[i]) 
	{
	    splitmin=Min(splitmin, split[i]);
	    splitmax=Max(splitmax, split[i]);	    
	    nt += NbOfSubTriangle(split[i]);
	    area+= Th[i].area; // error Nov 2010 FH ..
	}
    
    bool constsplit=splitmin==splitmax;
    bool noregenereration = constsplit || WithMortar;

    if(verbosity>2) 
	cout << "  -  Mesh construct : from " << &Th << " split min " << splitmin 
	    << "  max " << splitmax << ", recreate " << !noregenereration 
	    << " label =" << label << endl;
    
    triangles = new Triangle[nt];
    assert(triangles);
    //  computation of thee numbers of vertices
    //  on decoupe tous betement
    // et on recolle le points ensuite 
    // avec le quadtree qui sont sur les aretes ou les sommets 
    //  calcul du magorant du nombre de sommets 
    int nvmax= 0;// Th.nv;
    { //int v;
	KN<bool> setofv(Th.nv,false);
	for (int i=0;i<Th.nt;i++)
	    if ( split[i] ) 
	    {   nbsdd++;
		for (int j=0;j<3;j++)
		{
		    int jt=j,it=Th.ElementAdj(i,jt);
		    if(it==i || it <0) neb += split[i];  //on est sur la frontiere
		    else if  (!split[it]) neb += split[i];//le voisin ne doit pas etre decoupe
		    else  //on est dans le domaine et le voisin doit etre decoupe
		    {
			int ie0,ie1;
			Th.VerticesNumberOfEdge(Th[i],j,ie0,ie1);
			BoundaryEdge * pbe = Th.TheBoundaryEdge(ie0,ie1);
			if(pbe && &(*pbe)[0] == &Th(ie0))   
			    neb += max(split[i],split[it]); // aretes frontiere (FH juillet 2005)
			if (!pbe && (ie0 < ie1))
			{
			    nebi += max(split[i],split[it]); // arete interne a force ...  (FH jan 2007) 
			   // cout << nebi << " zzzz  t" << i << " a=" << j << " taj=" << it << " sa =" <<    ie0 << " " << ie1 << " + = " << max(split[i],split[it]) << endl;
			}
		    }
		}
		for (int j=0;j<3;j++)
		    //   if ( setofv.insert(Th(i,j) ).second ) 
		    if ( !setofv[Th(i,j)]) 
		    { 
			setofv[Th(i,j)]=true;
			//  cout << Th(i,j)  << " " ;
			nvmax++; 
		    }
	    }
		if(verbosity>4)
		    cout << "  - nv old " << nvmax << endl;
    }
		
	int nebmax=neb;
	int nebimax=nebi;   
	int nbsddmax=nbsdd;
	for (int i=0;i<Th.nt;i++)
	    if(split[i])
		nvmax += NbOfSubInternalVertices(split[i]) -3;
	
	//  compute the minimal Hsize of the new mesh  
	R hm = Th[0].h_min();
	//  cout << " hm " << hm << endl;
	// change h() in h_min() bug correct  july 2005 FH  ----
	for (int it=0;it<Th.nt;it++)
	{ 
	    assert(split[it]>=0 && split[it]<=64);
	    //      cout << " it = " <<it << " h " <<  Th[it].h() << " " << split[it] << " hm = " << hm <<  endl;
	    if (split[it]) 
		hm=Min(hm,Th[it].h_min()/(R) split[it]);
	}
	R seuil=hm/400.0;
	if(verbosity>5)   
	    cout << " seuil = " <<  seuil << " hmin = " << hm <<  endl; 
	assert(seuil>1e-15);
	vertices = new Vertex[nvmax];
	assert( vertices );
	
	nv =0;
	quadtree = new FQuadTree(this,Pmin,Pmax,nv); // build empty the quadtree
        {  // to keep the order of old vertices to have no problem in
             // interpolation on sub grid of RT finite element for example (Feb. 2016 F. Hecht)
        KN<bool> setofv(Th.nv,false);
        for (int i=0;i<Th.nt;++i)
          if (split[i])
              for (int j=0;j<3;++j)
              setofv[Th(i,j)]=true;
        for (int i=0;i<Th.nv;++i)
            if(setofv[i])
            {
                Vertex * pV=quadtree->ToClose(Th(i),seuil);
                if (pV ==0)
                {
                 vertices[nv] = Th(i);
                    if(verbosity>99) cout << " old to new: " << i << " -> " << nv << " / " << Th(i) <<endl;
                 vertices[nv].normal=0;
                 quadtree->Add(vertices[nv]);
                  nv++;
                }

            }
            if(verbosity>3)
            {
                cout << "  --- number of old vertices use: " << nv << endl;
                cout << "  --- number of  neb : " << nebmax << endl;
            }
        }
    /*
	for (int i=0;i<Th.nt;i++)
	    if (split[i]) 
		for (int j=0;j<3;j++)
		{
		    
		    Vertex * pV=quadtree->ToClose(Th[i][j],seuil);
		    if (pV ==0) { 
			// cout << Th(i,j) << " " ;
			vertices[nv] = Th[i][j];
			vertices[nv].normal=0;
			quadtree->Add(vertices[nv]);
			nv++;}
		}
		    // nv = Th.nv;
		    if(verbosity>3) 
		    {
			cout << "  -- number of old vertices use: " << nv << endl;      
			cout << "  -- number of  neb : " << nebmax << endl;
		    }
     */
	bedges = new BoundaryEdge[nebmax];
	BoundaryEdge * bedgesi= new BoundaryEdge[nebimax];  // jan 2007 FH
	int * sdd= new int[Th.nt]; 
	for (int i=0;i<Th.nt;++i)
		sdd[i]= 0; 
	/*   
	    for (int i=0;i<Th.neb;i++)
	{
		cout << i << " " << " " << Th(Th.bedges[i][0]) << " " << Th(Th.bedges[i][1]) << endl;
	} */
	assert(bedges && bedgesi);
	//  generation of the boundary edges
	neb =0;
	nebi=0;   //   generaztion des arete interne pour les force ...   
	nbsdd=0;
	//   for (int ieb=0,jj;ieb<Th.neb;ieb++)
	for (int it=0;it<Th.nt;it++)
	    if (  split[it] ) 
	    {
		// sdd[it]=-1;
		for (int jt=0;jt<3;jt++)
		{
		    int jtt=jt,itt=Th.ElementAdj(it,jtt);
		    //  cout << it <<  " " << jt << " " << jt << " " << itt << !split[itt] << endl;
		    int ie0,ie1;
		    Label  re(label); 
		    Th.VerticesNumberOfEdge(Th[it],jt,ie0,ie1);
		    //  cout << "++ ";
		    BoundaryEdge * pbe = Th.TheBoundaryEdge(ie0,ie1);
		    bool bbe= ( itt == it || itt <0 || !split[itt] || (pbe && &(*pbe)[0] == &Th(ie0))) ;
		    
		    BoundaryEdge * bbedges = bbe ? bedges : bedgesi;
		    int &  nneb = bbe ? neb : nebi;
		    int nnebmax = bbe ? nebmax : nebimax;
		    int offset= bbe ? 0 : nebmax;
		    if (bbe ||  (!pbe && (ie0 < ie1) ) ) // arete interne ou frontiere 
		    {
			int sens = 1; // par defaul le bon sens 
			int kold = it;   //Th.BoundaryElement(ieb,jj);
			int n=split[kold];
			if( itt>=0) n = max(n,split[itt]); //  pour les aretes internes (FH juillet 2005)
			if (!n) continue; 
			if (pbe ) {
			    re = *pbe;
			    if( & (*pbe)[0] == &Th(ie1) ) sens = -sens; // pour les aretes non decoupe avril 2007
			   // cout << " ### " << ie0 << " " << ie1 << " " <<  sens << endl;
			}
			// cout << " lab = " <<  re.lab << endl;
			Vertex *pva= quadtree->NearestVertex(Th(ie0));
			Vertex *pvb= quadtree->NearestVertex(Th(ie1));
			Vertex *pv0=pva;
			R2 A(*pva),B(*pvb);
			R la = 1,lb=0, delta=1.0/n;
			
			for (int j=1;j<n;j++) 
			{ 
			    sens = 1; //  arete decoupe => le sens change avril 2007
			    la-=delta;
			    lb+=delta;
			    assert(nv<nvmax);
			    Vertex *pv1= vertices + nv;
			    
			    (R2 &) *pv1 = A*la + B*lb;
			    (Label &) *pv1 = re ; //= (Label) be;
			    quadtree->Add(*pv1);
			    nv++;
			    
			    assert(nneb<nnebmax);
			    bbedges[nneb].vertices[0]=pv0;
			    bbedges[nneb].vertices[1]=pv1;
			    (Label &)  bbedges[nneb] = re;
			    nneb++;
			    
			    pv0=pv1;
			}
			//if (bbe)
			 //   cout << nneb+1 << " yyyy  t" << it << " a=" << jt << " taj=" << itt << " sa =" <<    ie0 << " " << ie1 << " + = " << n << endl;

			assert(nneb<nnebmax);
			bbedges[nneb].vertices[0]=pv0;
			bbedges[nneb].vertices[1]=pvb;
			(Label &) bbedges[nneb]= re ;		         			    
			sdd[it]= (1+ (nneb++ + offset))*sens; // numero de la derniere arete attention au sens si pas decoupe avril 2007
			//cout << " ### " << pv0 - vertices << " " << pvb - vertices << " " <<  sens << " it = " << it << " " <<  sdd[it] << " itt " << itt  << " -- " << sdd[406] << endl;

			if(  ( itt !=it || itt <0)  ) // interne 
			   sdd[itt]= -sdd[it]; 
		    }
		} 
		nbsdd++; 
	    }

		//   cout << "          (debug)  nb vertices on egdes " << nv << endl;  
		    //  cout << " " <<  nebmax << " " << neb << endl;
	assert(neb==nebmax);
	assert(nebi==nebimax);
	assert(nbsdd==nbsddmax);
	
	//   create the new vertices and the new triangle 
	int   kt=0;
	for (int K=0;K<Th.nt;K++)
	{ // cout << K << endl;
	    Triangle &T(Th[K]);
	    R2 A(T[0]);
	    R2 B(T[1]);
	    R2 C(T[2]);
	    long N=split[K];
	    if (!N) continue;
	    long N2=N*N;
	    int vt[3];
	    for (int n=0;n<N2;n++,kt++) //  loop on all sub triangle
	    {
		//(Label&) triangles[kt] = (Label&) T; 
		for(int j=0;j<3;j++) //  Loop on the 3 vertices
		{ 
		    R2 PTj=SubTriangle(N,n,j);
		    R la=1-PTj.x-PTj.y,lb=PTj.x,lc=PTj.y;
		    R lmin = Min(la,lb,lc);
		    R2 Pj= A*la+B*lb+C*lc; 
		    Vertex *pV;
		    pV=quadtree->ToClose(Pj,seuil);
		    // if !noregenereration => point du bord du triangle deja genere 
		    bool addv = !pV;
		    if(!noregenereration && pV==0) addv = lmin > 1e-5;
		    if ( addv )			
		    { // new vertex
		      // cout << "    -- " << nv << "  New Vertices " << Pj << " n=" << n << " j=" << j << " " << PTj << endl;
			assert(nv<nvmax);
			vertices[nv]=Pj;
			(Label&) vertices[nv]=0; //  Internal vertices 
			pV = vertices + nv;
			quadtree->Add(*pV);
			nv++;                  
		    }  //  end of new vertex
		    if(noregenereration)
		    {
		      assert(pV);
		    //cout << j <<  "  " << *pV << " " << Pj << " " << number(pV) << " " << Norme2(Pj-*pV) << " " << nv << endl;
		      vt[j]=number(pV);
		    }
		}
		   
		    //  triangles[kt].SetVertex(j,pV);
		 // for(int j=0;j<3;j++)
		  // cout << kt << " " << n << endl;
		if(noregenereration)
		{
		    R2 A=vertices[vt[0]];
		    R2 B=vertices[vt[1]];
		    R2 C=vertices[vt[2]]; 
		    R a = (( B-A)^(C-A))*0.5;
		    
		    if (a>0) 
			triangles[kt].set(vertices,vt[0],vt[1],vt[2],T.lab);
		    else 
			triangles[kt].set(vertices,vt[0],vt[2],vt[1],T.lab);
		}
	    }   // end loop on all sub triangle
	    
	} //  end 
	if (verbosity>3 ) 
	{ 
	    cout << "  - regeneration = " << ! noregenereration <<endl; 
	    cout << "  - Nb of vertices       " << nv << endl; 
	    cout << "  - Nb of triangle       " << nt << endl; 
	    cout << "  - Nb of boundary edges " << neb << endl; 
	
	}
	//  
	if (!noregenereration )   // REGENRATION DU MAILLAGE 
	{ // 
	   // nebi=0;
	    long nba = neb+nebi;
	    // 
	    long nbsd = nbsdd; // bofbof 
			   //ok,  with correction of mshptg 
	    long  *sd;
	    sd=new long[2*nbsd+2];
	    
	    sd[0]=-1;
	    sd[1]=Th[0].lab;         

	    
	    long nbs=nv;
	    long nbsmax=nv;
	    long            err = 0;//, nbsold = nbs;       
		long           *c = 0;
		long           *tri = 0;
		long           *nu = 0;
		long           *reft = 0;
		typedef double Rmesh ;
		Rmesh          *cr = 0;
		Rmesh          *h = 0;
		long nbtmax = 2 * nbsmax;
		long * arete  = new long[2*nba]; 
		nu = new long[6*nbtmax];
		c = new long[2*nbsmax];
		tri = new long[Max((4 * nbsmax + 2 * nbsd),nba)];
		reft = new long[nbtmax];
		cr = new Rmesh[(2 * nbsmax + 2)];
		h = new Rmesh[nbsmax];
		for(int i=0,k=0; i<nv; i++)
		{ 
		    cr[k++]  =vertices[i].x;
		    cr[k++]=vertices[i].y;
		    h[i]=1; 
		}
		{
		    int k=0;
		for (int i=0 ;i<neb;i++)
		{
		    arete[k++] =number(bedges[i][0])+1;
		    arete[k++] =number(bedges[i][1])+1;
		}
		//  ajoute des aretes interne pour les forces ..  FH 
		for (int i=0 ;i<nebi;i++)
		{
		    arete[k++] =number(bedgesi[i][0])+1;
		    arete[k++] =number(bedgesi[i][1])+1;
		}
		}
		// construction du tableau sd
		
		for (int it=0,j=0;it<Th.nt;it++)
		    if (  split[it]) // un ssd par vieux triangle 
		    {
			assert(sdd[it]);
			sd[j++]=sdd[it];
			sd[j++] = Th[it].lab; 
		    }
			
		extern int 
		    mshptg8_ (Rmesh *cr, Rmesh *h, long *c, long *nu, long *nbs, long nbsmx, long *tri,
			      long *arete, long nba, long *sd,
			      long nbsd, long *reft, long *nbt, Rmesh coef, Rmesh puis, long *err);
		
		long nbt=0;
		if(verbosity>10)
		{
		    cout << " mshptg8_ " << endl;
		    cout << "    nbs =" << nbs << endl;
		    cout << "    nba =" << nba << endl;
		    cout << "    nbt =" << nbt << endl;
		    cout << "    nbsd =" << nbsd << endl;
		    cout << " sommets : " << endl;
		    for(int i=0;i<nbs; ++i)
			cout << " " << i+1 << " -- " << cr[2*i] << " " << cr[2*i+1] << endl;
		    cout << " aretes : " << endl;
		    for(int i=0;i<nba; ++i)
			cout << " " << i+1 << " -- " << arete[2*i] << " " << arete[2*i+1] << endl;
		    cout << " Sd : " << endl;
		    for(int i=0;i<nbsd; ++i)
			cout << " " << i+1 << " -- " << sd[2*i] << "  lab " << sd[2*i+1] << endl;
		    //nbsd=0;
		}
		mshptg8_ (cr, h, c, nu, &nbs, nbs, tri, arete, nba, (long *) sd, nbsd, reft, &nbt, .25, .75, &err);
		if( !(err==0 && nbt !=0))
		{
		    cerr << " Error mesh generation in mshptg : err = " << err << " nb triangles = " << nbt << endl;
		    ffassert(err==0 && nbt !=0);
		}
		// 
		
		delete [] triangles;
		// Correction FH bug  trunc  mesh with hole 25032005 +july 2005 
		int kt=0;
		R dmin=1e100;
		for(int i=0,k=0;i<nbt;i++,k+=3)
		{
		    //int ir = reft[i]; 
		    reft[i]=-1;
		    R2 A=vertices[nu[k]-1],B=vertices[nu[k+1]-1],C=vertices[nu[k+2]-1];
		    R2 G=(A+B+C)/3.,PHat;
		    double d=Area2(A,B,C);
		    dmin=min(d,dmin);
		    //cout << A << "\n" << B << "\n" << C << "\n" << A << "\n\n";
		    if(d<=1e-5 && 0)
		    {
			cout<< " T = "<<  i << " det= " << d << "  ::  " << A << " " << B << " " << C  << endl;
			
		    }
		    bool outside;
		    const Triangle * t=Th.Find(G,PHat,outside,0);
		    if(!outside ) { 
			int k=Th(t);
			if( split[k] )
			{
			    reft[i] = k;
			    kt++;
			}
		    }
		}
		nt=kt;
		// cout << " " << dmin << endl;
		if(verbosity>3)      
		    cout << "  - Nb Triangles = " << nt <<  " remove triangle in hole :" <<  nbt - nt 
			<<endl;
		triangles = new Triangle[nt];
		kt=0;
		for(int i=0,k=0;i<nbt;i++,k+=3)
		    if(reft[i]>=0) 
			triangles[kt++].set(vertices,nu[k]-1,nu[k+1]-1,nu[k+2]-1,Th[reft[i]].lab);
		assert(kt==nt);
		// END  Correction FH bug  trunc  mesh with hole 25032005 + july 2005 
		
		delete [] arete;
		delete [] nu;
		delete [] c;
		delete [] tri;
		delete [] reft;
		delete [] cr;
		delete [] h;
		delete [] sd;
		
	}
	delete [] bedgesi;   
	delete [] sdd;   
	ConsAdjacence();
}


const char Mesh::magicmesh[8]="Mesh 2D";
int  Mesh::kthrough =0;
int  Mesh::kfind=0;

Mesh::Mesh(const  Serialize &serialized) 
{
    TriangleConteningVertex =0;
    BoundaryAdjacencesHead=0;
    BoundaryAdjacencesLink=0;
    BoundaryEdgeHeadLink=0;
    quadtree =0;
    volume=0;
    NbMortars=0;
    dim=0;
    tet=0;
    edges=0;
    ntet=0;
    ne=0;
    
    mortars=0;
    TheAdjacencesLink =0;
    nv=0;
    nt =0;
    neb=0;
    triangles=0;
    vertices=0;
    bedges=0;
    area=0;
    lenbord=0;
    bnormalv=0;
    //  ---  assert(serialized.samewhat(magicmesh));
    size_t  pp=0;;
    long long l;
    serialized.get(pp,l);
    serialized.get( pp,nt);
    serialized.get( pp,nv);
    serialized.get( pp,neb);
    if (verbosity>2) 
	cout << " mesh serialized : " << l << " " << nt << " " << nv << " " << neb << endl;
    assert ( nt > 0 && nv >0 && neb >=0);
    triangles = new Triangle [nt];
    vertices  = new Vertex[nv];
    bedges    = new BoundaryEdge[neb];
    area=0;
    ffassert(triangles && vertices && bedges);
    
    for (int i=0;i<nv;i++)    
    {
        serialized.get(pp,vertices[i].x);
        serialized.get(pp,vertices[i].y);
        serialized.get(pp,vertices[i].lab);
    }
    for (int i=0;i<nt;i++) { 
        int i0,i1,i2,ir;
        serialized.get(pp,i0);
        serialized.get(pp,i1);
        serialized.get(pp,i2);
        serialized.get(pp,ir);
	
        triangles[i].set(vertices,i0,i1,i2,ir); 
        area += triangles[i].area;}
    
    for (int i=0;i<neb;i++) { 
        int i0,i1,ir;
        serialized.get(pp,i0);
        serialized.get(pp,i1);
        serialized.get(pp,ir);
        bedges[i] = BoundaryEdge(vertices,i0,i1,ir);
	lenbord += bedges[i].length(); }
    assert( pp ==  serialized.size() );   
    if(verbosity>2) 
	cout << "   End of un serialize: area on mesh = " << area <<endl;  
    ConsAdjacence();	
    
}

Serialize  Mesh::serialize() const
{
    
    long long  l=0;
    l += sizeof(long long);
    l += 3*sizeof(int);
    l += nt*(sizeof(int)*4);
    l += nv*( sizeof(int) + sizeof(double)*2);
    l += neb*(sizeof(int)*3);
    
    // cout << l << magicmesh << endl;
    Serialize  serialized(l,magicmesh);
    // cout << l << magicmesh << endl;
    size_t pp=0;
    serialized.put(pp, l); 
    serialized.put( pp,nt);
    serialized.put( pp,nv);
    serialized.put( pp,neb);
    if (verbosity>2) 
      cout << " mesh Serialized : " << l << " "  << nt << " " << nv << " " << neb << endl;
    for (int i=0;i<nv;i++)
    {
	serialized.put(pp, vertices[i].x);
	serialized.put(pp, vertices[i].y);
	serialized.put(pp, vertices[i].lab);
    }
    //cout << " ---  " << endl;
    for (int i=0;i<nt;i++)
    {
	const Triangle & K(triangles[i]);
	int i0= &K[0]-vertices;
	int i1= &K[1]-vertices;
	int i2= &K[2]-vertices;
	serialized.put(pp, i0);
	serialized.put(pp, i1);
	serialized.put(pp, i2);
	serialized.put(pp, K.lab);
    }
    // cout << " ---  " << endl;
    for (int i=0;i<neb;i++)
    {
	const BoundaryEdge & K(bedges[i]);
	int i0= &K[0]-vertices;
	int i1= &K[1]-vertices;
	serialized.put(pp, i0);
	serialized.put(pp, i1);
	serialized.put(pp, K.lab);
    }
    //  cout << " ---  " << endl;
    // cout << pp << " == " << serialized.size() << endl;
    assert(pp==serialized.size());
    return serialized; 
}

void Mesh::Buildbnormalv()
{
    if (bnormalv) 
    {assert(0);return;}
    int nb=0;
    for (int k=0;k<nt;k++)
	for (int i=0;i<3;i++)
	{  
	    int ii(i),kk;
	    kk=ElementAdj(k,ii);
	    if (kk<0 || kk==k) nb++;
	}
	    if(verbosity>2)
		cout << " number of real boundary edges " << nb << endl;
    bnormalv= new R2[nb];
    R2 *n=bnormalv;
    for (int k=0;k<nt;k++)
	for (int i=0;i<3;i++)
	{  
	    int ii(i),kk;
	    kk=ElementAdj(k,ii);
	    if (kk<0 || kk==k) {
		Triangle & K(triangles[k]);
		R2 N=K.n(i);
		Vertex & v0= K.Edge(0,i);
		Vertex & v1= K.Edge(1,i);
		v0.SetNormal(n,N);
		v1.SetNormal(n,N);         
	    }
	}
	    // cout << nb << " == " << n-bnormalv << endl;
	    assert(n - bnormalv <= nb );
}
}
//  static const R2 Triangle::Hat[3]= {R2(0,0),R2(1,0),R2(0,1)} ; 
