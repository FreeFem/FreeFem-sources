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
#include <cmath>
#include <cstdlib>
#include "error.hpp"
#include <iostream>
#include <fstream>
#include <map>
#include "rgraph.hpp"
using namespace std;  

#include "RNM.hpp"
#include "fem.hpp"
#include "FESpace.hpp"

namespace  Fem2D {

class TypeOfFE_RT : public  TypeOfFE { public:  
  static int Data[];
   TypeOfFE_RT(): TypeOfFE(0,1,0,2,Data,1,1,6,3) 
     {const R2 Pt[] = { R2(0.5,0.5), R2(0.0,0.5), R2(0.5,0.0) };
      for (int p=0,kk=0;p<3;p++)
       { P_Pi_h[p]=Pt[p];   
        for (int j=0;j<2;j++) 
        pij_alpha[kk++]= IPJ(p,p,j);
       }}
  // void FB(const Mesh & Th,const Triangle & K,const R2 &P, RNMK_ & val) const;
   void FB(const bool * watdd, const Mesh & Th,const Triangle & K,const R2 &P, RNMK_ & val) const;
//   void D2_FB(const Mesh & Th,const Triangle & K,const R2 &P, RNMK_ & val) const;
  // void Pi_h(const baseFElement & K,RN_ & val, InterpolFunction f, R* v,int, void *) const;
   void Pi_h_alpha(const baseFElement & K,KN_<double> & v) const ;
} ; 
//                     on what     nu df on node node of df    
  int TypeOfFE_RT::Data[]={3,4,5,       0,0,0,       0,1,2,       0,0,0,        0,1,2,   0,0, 0,0,3,3};

/* void TypeOfFE_RT::D2_FB(const Mesh & ,const Triangle & ,const R2 & ,RNMK_ & val) const
{ //  
  val=0;
}*/
/*
 void TypeOfFE_RT::FB(const Mesh & Th,const Triangle & K,const R2 & PHat,RNMK_ & val) const
{ //  
//  const Triangle & K(FE.T);
  R2 P(K(PHat));
  R2 A(K[0]), B(K[1]),C(K[2]);
  R l0=1-P.x-P.y,l1=P.x,l2=P.y; 
 // R2 Dl0(K.H(0)), Dl1(K.H(1)), Dl2(K.H(2));
  if (val.N() <3) 
   throwassert(val.N() >=3);
  throwassert(val.M()==2 );
  throwassert(val.K()==3 );
  RN_ f0(val('.',0,0)); 
  RN_ f1(val('.',1,0)); 
  val=0;  
//  RN_ df0(val(0,'.',0)); 
//  RN_ fy(val('.','.',2)); 
  //     a_i ([x,y]-c_i) , ou  c_i = A,B , C si i= 0,1,2 
  //   int_T a_i div([x,y]-c_i) = 1
  //    div div([x,y]-c_i) = 2 
  //   donc a_i = 1/(2 area T)
  
  R a=1./(2*K.area);
  R a0=   K.EdgeOrientation(0) * a ;
  R a1=   K.EdgeOrientation(1) * a  ;
  R a2=   K.EdgeOrientation(2) * a ;
 // if (Th(K)< 2) cout << Th(K) << " " <<  A << " "  << B << " " << C << "; " <<  a0 << " " << a1 << " "<< a2 << endl;;

  //  ------------
  f0[0] = (P.x-A.x)*a0;
  f1[0] = (P.y-A.y)*a0;
  
  f0[1] = (P.x-B.x)*a1;
  f1[1] = (P.y-B.y)*a1;
  
  f0[2] = (P.x-C.x)*a2;
  f1[2] = (P.y-C.y)*a2;
  // ----------------
  // ----------------
  // BUG dans RT correct FH le 17 sept 2002 
  //  dx [x,y] = [1,0] et non [1,1]
  //  dy [x,y] = [0,1] et non [1,1] 
  // -------------------------------------
  
  val(0,0,1) =  a0;  
  val(1,0,1) =  a1;  
  val(2,0,1) =  a2;  
  val(0,1,2) =  a0;  
  val(1,1,2) =  a1;  
  val(2,1,2) =  a2;  
  
}
*/
 void TypeOfFE_RT::FB(const bool *whatd,const Mesh & Th,const Triangle & K,const R2 & PHat,RNMK_ & val) const
{ //  
//  const Triangle & K(FE.T);
  R2 P(K(PHat));
  R2 A(K[0]), B(K[1]),C(K[2]);
  // R l0=1-P.x-P.y,l1=P.x,l2=P.y; 
 // R2 Dl0(K.H(0)), Dl1(K.H(1)), Dl2(K.H(2));
  if (val.N() <3) 
   throwassert(val.N() >=3);
  throwassert(val.M()==2 );
//  throwassert(val.K()==3 );
  val=0;     
  R a=1./(2*K.area);
  R a0=   K.EdgeOrientation(0) * a ;
  R a1=   K.EdgeOrientation(1) * a  ;
  R a2=   K.EdgeOrientation(2) * a ;
 // if (Th(K)< 2) cout << Th(K) << " " <<  A << " "  << B << " " << C << "; " <<  a0 << " " << a1 << " "<< a2 << endl;;

  //  ------------
  if (whatd[op_id])
   {
   assert(val.K()>op_id);
  RN_ f0(val('.',0,0)); 
  RN_ f1(val('.',1,0)); 
  f0[0] = (P.x-A.x)*a0;
  f1[0] = (P.y-A.y)*a0;
  
  f0[1] = (P.x-B.x)*a1;
  f1[1] = (P.y-B.y)*a1;
  
  f0[2] = (P.x-C.x)*a2;
  f1[2] = (P.y-C.y)*a2;
  }
  // ----------------
  // BUG dans RT correct FH le 17 sept 2002 
  //  dx [x,y] = [1,0] et non [1,1]
  //  dy [x,y] = [0,1] et non [1,1] 
  // -------------------------------------
    if (whatd[op_dx])
   {
   assert(val.K()>op_dx);
   val(0,0,op_dx) =  a0;  
   val(1,0,op_dx) =  a1;  
   val(2,0,op_dx) =  a2; 
  } 
    if (whatd[op_dy])
   {
    assert(val.K()>op_dy);
    val(0,1,op_dy) =  a0;  
    val(1,1,op_dy) =  a1;  
    val(2,1,op_dy) =  a2;  
  }
}
/*
 void TypeOfFE_RT::Pi_h(const baseFElement & K,RN_ & val, InterpolFunction f, R* v,int j,  void * arg) const
{
  const R2 Pt[] = { R2(0.5,0.5), R2(0.0,0.5), R2(0.5,0.0) };
  const Triangle & T(K.T);
  //    if (K.number<2) cout << K.number << ": " ;
   for (int i=0;i<3;i++)
     {  
       f(v,T(Pt[i]),K,i,Pt[i],arg);
       R2 E(T.Edge(i));
       R signe = T.EdgeOrientation(i) ;
       val[i]= signe*(v[j]*E.y-v[j+1]*E.x); 
    //   if (K.number<2) cout <<  val[i] << " " ;
       }
   //   if (K.number<2) cout << endl;
}
*/
void TypeOfFE_RT::Pi_h_alpha(const baseFElement & K,KN_<double> & v) const 
{
  const Triangle & T(K.T);

   for (int i=0,k=0;i<3;i++)
     {  
        R2 E(T.Edge(i));
        R signe = T.EdgeOrientation(i) ;
        v[k++]= signe*E.y;
        v[k++]=-signe*E.x;
     }   
}
// -------------------



class TypeOfFE_RTmodif : public  TypeOfFE { public:  
  static int Data[];
   TypeOfFE_RTmodif(): TypeOfFE(0,1,0,2,Data,1,1,6,3) 
     {const R2 Pt[] = { R2(0.5,0.5), R2(0.0,0.5), R2(0.5,0.0) };
      for (int p=0,kk=0;p<3;p++)
       { P_Pi_h[p]=Pt[p];   
        for (int j=0;j<2;j++) 
        pij_alpha[kk++]= IPJ(p,p,j);
       }}
  // void FB(const Mesh & Th,const Triangle & K,const R2 &P, RNMK_ & val) const;
   void FB(const bool * whatd, const Mesh & Th,const Triangle & K,const R2 &P, RNMK_ & val) const;
 //  void D2_FB(const Mesh & Th,const Triangle & K,const R2 &P, RNMK_ & val) const;
  // void Pi_h(const baseFElement & K,RN_ & val, InterpolFunction f, R* v,int, void *) const;
   void Pi_h_alpha(const baseFElement & K,KN_<double> & v) const ;
} ; 
//                     on what     nu df on node node of df    
  int TypeOfFE_RTmodif::Data[]={3,4,5,       0,0,0,       0,1,2,       0,0,0,        0,1,2,   0,0,  0,0, 3,3};


 void TypeOfFE_RTmodif::FB(const bool * whatd,const Mesh & Th,const Triangle & K,const R2 & PHat,RNMK_ & val) const
{ //  
//  const Triangle & K(FE.T);
  R2 P(K(PHat));
  R2 A(K[0]), B(K[1]),C(K[2]);
  R la=1-PHat.x-PHat.y,lb=PHat.x,lc=PHat.y; 
  R2 Dla(K.H(0)), Dlb(K.H(1)), Dlc(K.H(2));
  if (val.N() <3) 
   throwassert(val.N() >=3);
  throwassert(val.M()==2 );

  R2 AB(A,B),AC(A,C),BA(B,A),BC(B,C),CA(C,A),CB(C,B);

  R aa0= 1./(((AB,Dlb) + (AC,Dlc))*K.area);
  R aa1= 1./(((BA,Dla) + (BC,Dlc))*K.area);
  R aa2= 1./(((CA,Dla) + (CB,Dlb))*K.area);
  int i=0;
  R a0=   &K[ (i+1)%3] < &K[ (i+2)%3] ? aa0 : -aa0 ;
  i=1;
  R a1=   &K[ (i+1)%3] < &K[ (i+2)%3] ? aa1 : -aa1 ;
  i=2;
  R a2=   &K[ (i+1)%3] < &K[ (i+2)%3] ? aa2 : -aa2 ;
 // if (Th(K)< 2) cout << Th(K) << " " <<  A << " "  << B << " " << C << "; " <<  a0 << " " << a1 << " "<< a2 << endl;;
 
  R2 Va= AB*(lb*a0) + AC*(lc*a0);
  R2 Vb= BA*(la*a1) + BC*(lc*a1);
  R2 Vc= CA*(la*a2) + CB*(lb*a2);
  R2 Va_x= AB*(Dlb.x*a0) + AC*(Dlc.x*a0);
  R2 Vb_x= BA*(Dla.x*a1) + BC*(Dlc.x*a1);
  R2 Vc_x= CA*(Dla.x*a2) + CB*(Dlb.x*a2);
  R2 Va_y= AB*(Dlb.y*a0) + AC*(Dlc.y*a0);
  R2 Vb_y= BA*(Dla.y*a1) + BC*(Dlc.y*a1);
  R2 Vc_y= CA*(Dla.y*a2) + CB*(Dlb.y*a2);
 
 if( whatd[op_id])
  {
    RN_ f0(val('.',0,0)); 
    RN_ f1(val('.',1,0)); 
    
    f0[0] = Va.x;
    f1[0] = Va.y;
    
    f0[1] = Vb.x;
    f1[1] = Vb.y;
    
    f0[2] = Vc.x;
    f1[2] = Vc.y;
  }
 // ----------------
 if( whatd[op_dx])
   {
     val(0,0,1) =  Va_x.x;  
     val(0,1,1) =  Va_x.y;  
     
     val(1,0,1) =  Vb_x.x;  
     val(1,1,1) =  Vb_x.y;  
     
     val(2,0,1) =  Vc_x.x;  
     val(2,1,1) =  Vc_x.y;
   }
 
 if( whatd[op_dy])
   {
     val(0,0,2) =  Va_y.x;  
     val(0,1,2) =  Va_y.y;  

     val(1,0,2) =  Vb_y.x;  
     val(1,1,2) =  Vb_y.y;  
     
  val(2,0,2) =  Vc_y.x;  
  val(2,1,2) =  Vc_y.y;  
   }
 
}

void TypeOfFE_RTmodif::Pi_h_alpha(const baseFElement & K,KN_<double> & v) const 
{
  const Triangle & T(K.T);

   for (int i=0,k=0;i<3;i++)
     {  
        R2 E(T.Edge(i));
        R signe = &T[ (i+1)%3] < &T[ (i+2)%3] ? 1.0 : -1.0 ;
        v[k++]= signe*E.y;
        v[k++]=-signe*E.x;
     }   
}


// ---------------------
class TypeOfFE_P0 : public  TypeOfFE { public:  
  static int Data[];
  static double Pi_h_coef[];

   TypeOfFE_P0(): TypeOfFE(0,0,1,1,Data,1,1,1,1,Pi_h_coef){
     pij_alpha[0]=IPJ(0,0,0);
     P_Pi_h[0]=R2(1./3.,1./3.);
     }
   void FB(const bool * watdd, const Mesh & Th,const Triangle & K,const R2 &P, RNMK_ & val) const;

};

//                     on what     nu df on node node of df    
  int TypeOfFE_P0::Data[]={6, 0, 0, 0 , 0 ,0, 0, 1 };
double TypeOfFE_P0::Pi_h_coef[]={1.0};


void TypeOfFE_P0::FB(const bool* whatd,const Mesh & ,const Triangle & K,const R2 & PHat,RNMK_ & val) const
{ //  
  //  const Triangle & K(FE.T);
  R2 P(K(PHat));
  R2 A(K[0]), B(K[1]),C(K[2]);
  //R l0=1-P.x-P.y,l1=P.x,l2=P.y; 
  R2 Dl0(K.H(0)), Dl1(K.H(1)), Dl2(K.H(2));
  throwassert(val.N() >=1);
  throwassert(val.M()==1 );
  // throwassert(val.K()==3 );
  val=0;
  if ( whatd[op_id])
    val(0,0,0) =1;
}

// ------ P1 non conforme --------
class TypeOfFE_P1ncLagrange : public  TypeOfFE { public:  
  static int Data[];
  static double Pi_h_coef[];
  TypeOfFE_P1ncLagrange(): TypeOfFE(0,1,0,1,Data,1,1,3,3,Pi_h_coef)
  {
    const R2 Pt[] = { R2(0.5,0.5), R2(0.0,0.5), R2(0.5,0.0) };
    for (int i=0;i<NbDoF;i++) {
      pij_alpha[i]= IPJ(i,i,0);
      P_Pi_h[i]=Pt[i]; }
  }

   void FB(const bool * whatd, const Mesh & Th,const Triangle & K,const R2 &P, RNMK_ & val) const;

} ;
    
//                     on what     nu df on node node of df    
  int TypeOfFE_P1ncLagrange::Data[]={3,4,5,       0,0,0,       0,1,2,       0,0,0,        0,1,2,       0, 0,3};
double TypeOfFE_P1ncLagrange::Pi_h_coef[]={1.,1.,1.};


class TypeOfFE_ConsEdge : public  TypeOfFE { public:  
  static int Data[];
    static double Pi_h_coef[];

   TypeOfFE_ConsEdge(): TypeOfFE(0,1,0,1,Data,3,1,3,3,Pi_h_coef)
    { const R2 Pt[] = { R2(0.5,0.5), R2(0.0,0.5), R2(0.5,0.0) };
      for (int i=0;i<NbDoF;i++) {
       pij_alpha[i]= IPJ(i,i,0);
       P_Pi_h[i]=Pt[i]; }
     }

   void FB(const bool * whatd, const Mesh & Th,const Triangle & K,const R2 &P, RNMK_ & val) const;
  
} ;
//                     on what     nu df on node node of df    
  int TypeOfFE_ConsEdge::Data[]={3,4,5,       0,0,0,       0,1,2,       0,0,0,        0,1,2,       0, 0,3};
double TypeOfFE_ConsEdge::Pi_h_coef[]={1.,1.,1.};
void TypeOfFE_ConsEdge::FB(const bool * whatd,const Mesh & ,const Triangle & K,const R2 & P,RNMK_ & val) const
{
  //  const Triangle & K(FE.T);
  R2 A(K[0]), B(K[1]),C(K[2]);
  R l0=1-P.x-P.y,l1=P.x,l2=P.y; 
  
  if (val.N() <3) 
   throwassert(val.N() >=3);
  throwassert(val.M()==1 );
  
  val=0; 
  if (whatd[op_id])
   {
     
     RN_ f0(val('.',0,0)); 
     //      
    f0[0] =  double(l0 <= min(l1,l2) ); // arete  
    f0[1] =  double(l1 <= min(l0,l2) );
    f0[2] =  double(l2 <= min(l0,l1) );
  }

}
/*    
class TypeOfFE_P1Edge : public  TypeOfFE { public:  
  static int Data[];
  static double Pi_h_coef[];
	
  TypeOfFE_P1Edge(): TypeOfFE(0,2,0,1,Data,3,1,12,6,Pi_h_coef)
    {  R2 Pt[6] ;
	
	int kk=0;
	for(int i=0;i<3;++i)
	  for(int j=0;i<QF_GaussLegendre2.n;++j)
	  { R2 A(TriangleHat[VerticesOfTriangularEdge[i][0]]);
	    R2 B(TriangleHat[VerticesOfTriangularEdge[i][1]]);
	     Pt[k++]=A*(QF_GaussLegendre2[j].x)+ B*(1-QF_GaussLegendre2[j].x)
	  }
	      
        int other[6]= { 1,0, 3,2,5,4 };
	k=0;
	  for (int i=0;i<NbDoF;i++) {
	    pij_alpha[i]= IPJ(i,i,0);
	    if(other[i]>=0)
		pij_alpha[kk++]= IPJ(i,other[i],0);
	       P_Pi_h[i]=Pt[i]; }
	}
	
	void FB(const bool * whatd, const Mesh & Th,const Triangle & K,const R2 &P, RNMK_ & val) const;
	
    } ;
    //                     on what     nu df on node node of df    
    int TypeOfFE_ConsEdge::Data[]={3,4,5,       0,0,0,       0,1,2,       0,0,0,        0,1,2,       0, 0,3};
    double TypeOfFE_ConsEdge::Pi_h_coef[]={1.,1.,1.};
    void TypeOfFE_ConsEdge::FB(const bool * whatd,const Mesh & ,const Triangle & K,const R2 & P,RNMK_ & val) const
    {
	//  const Triangle & K(FE.T);
	R2 A(K[0]), B(K[1]),C(K[2]);
	R l0=1-P.x-P.y,l1=P.x,l2=P.y; 
	
	if (val.N() <3) 
	    throwassert(val.N() >=3);
	throwassert(val.M()==1 );
	
	val=0; 
	if (whatd[op_id])
	{
	    
	    RN_ f0(val('.',0,0)); 
	    //      
	    f0[0] =  double(l0 <= min(l1,l2) ); // arete  
	    f0[1] =  double(l1 <= min(l0,l2) );
	    f0[2] =  double(l2 <= min(l0,l1) );
	}
	
    }
    
*/
    
void TypeOfFE_P1ncLagrange::FB(const bool * whatd,const Mesh & ,const Triangle & K,const R2 & P,RNMK_ & val) const
{
  //  const Triangle & K(FE.T);
  R2 A(K[0]), B(K[1]),C(K[2]);
  //  l1(  cshrink1*(cshrink*((1,0)-G)+G)-G)+G  = 1
  R l0=1-P.x-P.y,l1=P.x,l2=P.y; 
  
  if (val.N() <3) 
    throwassert(val.N() >=3);
  throwassert(val.M()==1 );
  // throwassert(val.K()==3 );
  
  val=0; 
  if (whatd[op_id])
   {
  RN_ f0(val('.',0,0)); 
  f0[0] = 1-l0*2;  
  f0[1] = 1-l1*2;
  f0[2] = 1-l2*2;
  }
  if (whatd[op_dx] || whatd[op_dy] )
   {
    R2 Dl0(K.H(0)), Dl1(K.H(1)), Dl2(K.H(2));
    if (whatd[op_dx]) 
      {
      RN_ f0x(val('.',0,op_dx)); 
      f0x[0] = -Dl0.x*2;
      f0x[1] = -Dl1.x*2;
      f0x[2] = -Dl2.x*2;
      }
    if (whatd[op_dy]) 
     {  
     RN_ f0y(val('.',0,op_dy)); 
     f0y[0] = -Dl0.y*2;
     f0y[1] = -Dl1.y*2;
     f0y[2] = -Dl2.y*2;
     }
   }
}

// The RT orthogonal FE

class TypeOfFE_RTortho : public  TypeOfFE { public:  
  static int Data[];
   TypeOfFE_RTortho(): TypeOfFE(0,1,0,2,Data,1,1,6,3) 
     {const R2 Pt[] = { R2(0.5,0.5), R2(0.0,0.5), R2(0.5,0.0) };
      for (int p=0,kk=0;p<3;p++)
       { P_Pi_h[p]=Pt[p];   
        for (int j=0;j<2;j++) 
        pij_alpha[kk++]= IPJ(p,p,j);
       }}

   void FB(const bool * watdd, const Mesh & Th,const Triangle & K,const R2 &P, RNMK_ & val) const;
   void Pi_h_alpha(const baseFElement & K,KN_<double> & v) const ;
} ; 
//                     on what     nu df on node node of df    
  int TypeOfFE_RTortho::Data[]={3,4,5,       0,0,0,       0,1,2,       0,0,0,        0,1,2,   0,0, 0,0, 3,3};



 void TypeOfFE_RTortho::FB(const bool *whatd,const Mesh & Th,const Triangle & K,const R2 & PHat,RNMK_ & val) const
{ //  
//  const Triangle & K(FE.T);
  R2 P(K(PHat));
  R2 A(K[0]), B(K[1]),C(K[2]);
  //R l0=1-P.x-P.y,l1=P.x,l2=P.y; 
 // R2 Dl0(K.H(0)), Dl1(K.H(1)), Dl2(K.H(2));
  if (val.N() <3) 
   throwassert(val.N() >=3);
  throwassert(val.M()==2 );
//  throwassert(val.K()==3 );
  val=0;  
  R a=1./(2*K.area);
  R a0=   K.EdgeOrientation(0) * a ;
  R a1=   K.EdgeOrientation(1) * a  ;
  R a2=   K.EdgeOrientation(2) * a ;
 // if (Th(K)< 2) cout << Th(K) << " " <<  A << " "  << B << " " << C << "; " <<  a0 << " " << a1 << " "<< a2 << endl;;

  //  ------------
  if (whatd[op_id])
   {
   assert(val.K()>op_id);
  RN_ f0(val('.',0,0)); 
  RN_ f1(val('.',1,0)); 
  f1[0] =  (P.x-A.x)*a0;
  f0[0] = -(P.y-A.y)*a0;
  
  f1[1] =  (P.x-B.x)*a1;
  f0[1] = -(P.y-B.y)*a1;
  
  f1[2] =  (P.x-C.x)*a2;
  f0[2] = -(P.y-C.y)*a2;
  }
  // ----------------
    if (whatd[op_dx])
   {
   assert(val.K()>op_dx);
   val(0,1,op_dx) =  a0;  
   val(1,1,op_dx) =  a1;  
   val(2,1,op_dx) =  a2; 
  } 
    if (whatd[op_dy])
   {
   assert(val.K()>op_dy);
    val(0,0,op_dy) =  -a0;  
    val(1,0,op_dy) =  -a1;  
    val(2,0,op_dy) =  -a2;  
  }
}

void TypeOfFE_RTortho::Pi_h_alpha(const baseFElement & K,KN_<double> & v) const 
{
  const Triangle & T(K.T);

   for (int i=0,k=0;i<3;i++)
     {  
        R2 E(T.Edge(i));
        R signe = &T[ (i+1)%3] < &T[ (i+2)%3] ? 1.0 : -1.0 ;
        v[k++]= signe*E.x;
        v[k++]= signe*E.y;
     }   
}
// -------------------
// ttdc finite element fully discontinue. 
// -------------------
class TypeOfFE_P1ttdc : public  TypeOfFE { public:  
  static int Data[];
  static double Pi_h_coef[];
    static const R2 G;
    static const R cshrink;
    static const R cshrink1;
    //  (1 -1/3)*
    
    static R2 Shrink(const R2& P){ return (P-G)*cshrink+G;}
    static R2 Shrink1(const R2& P){ return (P-G)*cshrink1+G;}
    
   TypeOfFE_P1ttdc(): TypeOfFE(0,0,3,1,Data,1,1,3,3,Pi_h_coef)
    { const R2 Pt[] = { Shrink(R2(0,0)), Shrink(R2(1,0)), Shrink(R2(0,1)) }; 
      for (int i=0;i<NbDoF;i++) {
       pij_alpha[i]= IPJ(i,i,0);
       P_Pi_h[i]=Pt[i];
       // cout << Pt[i] << " " ;
      }
      //	cout <<" cshrink: " << cshrink << " cshrink1 : "<< cshrink1 <<endl;
     }

   void FB(const bool * whatd,const Mesh & Th,const Triangle & K,const R2 &P, RNMK_ & val) const;
   
  
  virtual R operator()(const FElement & K,const  R2 & PHat,const KN_<R> & u,int componante,int op) const ;
   
} ;
    const R2 TypeOfFE_P1ttdc::G(1./3.,1./3.);   
    const R TypeOfFE_P1ttdc::cshrink=1-1e-2;   
    const R TypeOfFE_P1ttdc::cshrink1=1./TypeOfFE_P1ttdc::cshrink;   
    
class TypeOfFE_P2ttdc : public  TypeOfFE { public:  
  static int Data[];
  static double Pi_h_coef[];
    static const R2 G;
    static const R cshrink;
    static const R cshrink1;
    
    static R2 Shrink(const R2& P){ return (P-G)*cshrink+G;}
    static R2 Shrink1(const R2& P){ return (P-G)*cshrink1+G;}
    
   TypeOfFE_P2ttdc(): TypeOfFE(0,0,6,1,Data,3,1,6,6,Pi_h_coef)
    { const R2 Pt[] = { Shrink(R2(0,0)), Shrink(R2(1,0)), Shrink(R2(0,1)),
	                Shrink(R2(0.5,0.5)),Shrink(R2(0,0.5)),Shrink(R2(0.5,0)) };
      for (int i=0;i<NbDoF;i++) {
       pij_alpha[i]= IPJ(i,i,0);
       P_Pi_h[i]=Pt[i]; }
     }
   
  
   void FB(const bool * whatd,const Mesh & Th,const Triangle & K,const R2 &P, RNMK_ & val) const;
 
 
} ;
//                          on what   nu df on node  node of df    
  int TypeOfFE_P1ttdc::Data[]={6,6,6,       0,1,2,       0,0,0,       0,0,0,         0,1,2,       0, 0,3};
  int TypeOfFE_P2ttdc::Data[]={6,6,6,6,6,6, 0,1,2,3,4,5, 0,0,0,0,0,0,  0,0,0,0,0,0,  0,1,2,3,4,5, 0, 0,6};
double TypeOfFE_P1ttdc::Pi_h_coef[]={1.,1.,1.};
double TypeOfFE_P2ttdc::Pi_h_coef[]={1.,1.,1.,1.,1.,1.};
    
const R2 TypeOfFE_P2ttdc::G(1./3.,1./3.);   
const R TypeOfFE_P2ttdc::cshrink=1-1e-2;   
const R TypeOfFE_P2ttdc::cshrink1=1./TypeOfFE_P2ttdc::cshrink;   
    
 
R TypeOfFE_P1ttdc::operator()(const FElement & K,const  R2 & P1Hat,const KN_<R> & u,int componante,int op) const 
{ 

  R2 PHat=Shrink1(P1Hat);  
  R u0(u(K(0))), u1(u(K(1))), u2(u(K(2)));
  R r=0;
  if (op==0)
    {
      R l0=1-PHat.x-PHat.y,l1=PHat.x,l2=PHat.y; 
      r = u0*l0+u1*l1+l2*u2;
    }
   else
    { 
       const Triangle & T=K.T;
       R2 D0 = T.H(0)*cshrink1 , D1 = T.H(1)*cshrink1  , D2 = T.H(2)*cshrink1 ;
       if (op==1)
         r =  D0.x*u0 + D1.x*u1 + D2.x*u2 ;
        else 
         r =  D0.y*u0 + D1.y*u1 + D2.y*u2 ;
    }
 //  cout << r << "\t";
   return r;
}

void TypeOfFE_P1ttdc::FB(const bool *whatd,const Mesh & ,const Triangle & K,const R2 & P1,RNMK_ & val) const
{
  R2 P=Shrink1(P1);  
   
  //  const Triangle & K(FE.T);
  R2 A(K[0]), B(K[1]),C(K[2]);
  R l0=1-P.x-P.y,l1=P.x,l2=P.y; 
  
  if (val.N() <3) 
    throwassert(val.N() >=3);
  throwassert(val.M()==1 );
  //  throwassert(val.K()==3 );
  
  val=0; 
  RN_ f0(val('.',0,op_id)); 
  
  if (whatd[op_id]) 
   {
     f0[0] = l0;
    f0[1] = l1;
    f0[2] = l2;}
 if (whatd[op_dx] || whatd[op_dy])
  {
  R2 Dl0(K.H(0)*cshrink1), Dl1(K.H(1)*cshrink1), Dl2(K.H(2)*cshrink1);
  
  if (whatd[op_dx]) 
   {
    RN_ f0x(val('.',0,op_dx)); 
   f0x[0] = Dl0.x;
   f0x[1] = Dl1.x;
   f0x[2] = Dl2.x;
  }
  
  if (whatd[op_dy]) {
    RN_ f0y(val('.',0,op_dy)); 
   f0y[0] = Dl0.y;
   f0y[1] = Dl1.y;
   f0y[2] = Dl2.y;
  }
  }
}



 void TypeOfFE_P2ttdc::FB(const bool *whatd,const Mesh & ,const Triangle & K,const R2 & P1,RNMK_ & val) const
{
    R2 P=Shrink1(P1);  
   
//  const Triangle & K(FE.T);
  R2 A(K[0]), B(K[1]),C(K[2]);
  R l0=1-P.x-P.y,l1=P.x,l2=P.y; 
  R l4_0=(4*l0-1),l4_1=(4*l1-1),l4_2=(4*l2-1); 
  
//  throwassert(FE.N == 1);  
  throwassert( val.N()>=6);
  throwassert(val.M()==1);
//  throwassert(val.K()==3 );
  
  val=0; 
// --     
 if (whatd[op_id])
  {
   RN_ f0(val('.',0,op_id)); 
  f0[0] = l0*(2*l0-1);
  f0[1] = l1*(2*l1-1);
  f0[2] = l2*(2*l2-1);
  f0[3] = 4*l1*l2; // oppose au sommet 0
  f0[4] = 4*l0*l2; // oppose au sommet 1
  f0[5] = 4*l1*l0; // oppose au sommet 3
  }
 if(  whatd[op_dx] || whatd[op_dy] || whatd[op_dxx] || whatd[op_dyy] ||  whatd[op_dxy])
 {
   R2 Dl0(K.H(0)*cshrink1), Dl1(K.H(1)*cshrink1), Dl2(K.H(2)*cshrink1);
  if (whatd[op_dx])
  {
    RN_ f0x(val('.',0,op_dx)); 
  f0x[0] = Dl0.x*l4_0;
  f0x[1] = Dl1.x*l4_1;
  f0x[2] = Dl2.x*l4_2;
  f0x[3] = 4*(Dl1.x*l2 + Dl2.x*l1) ;
  f0x[4] = 4*(Dl2.x*l0 + Dl0.x*l2) ;
  f0x[5] = 4*(Dl0.x*l1 + Dl1.x*l0) ;
  }

 if (whatd[op_dy])
  {  
    RN_ f0y(val('.',0,op_dy)); 
  f0y[0] = Dl0.y*l4_0;
  f0y[1] = Dl1.y*l4_1;
  f0y[2] = Dl2.y*l4_2;
  f0y[3] = 4*(Dl1.y*l2 + Dl2.y*l1) ;
  f0y[4] = 4*(Dl2.y*l0 + Dl0.y*l2) ;
  f0y[5] = 4*(Dl0.y*l1 + Dl1.y*l0) ;
  }
 
 if (whatd[op_dxx])
  {  
    RN_ fxx(val('.',0,op_dxx)); 

    fxx[0] = 4*Dl0.x*Dl0.x;
    fxx[1] = 4*Dl1.x*Dl1.x;
    fxx[2] = 4*Dl2.x*Dl2.x;
    fxx[3] =  8*Dl1.x*Dl2.x;
    fxx[4] =  8*Dl0.x*Dl2.x;
    fxx[5] =  8*Dl0.x*Dl1.x;
  }

 if (whatd[op_dyy])
  {  
    RN_ fyy(val('.',0,op_dyy)); 
    fyy[0] = 4*Dl0.y*Dl0.y;
    fyy[1] = 4*Dl1.y*Dl1.y;
    fyy[2] = 4*Dl2.y*Dl2.y;
    fyy[3] =  8*Dl1.y*Dl2.y;
    fyy[4] =  8*Dl0.y*Dl2.y;
    fyy[5] =  8*Dl0.y*Dl1.y;
  }
 if (whatd[op_dxy])
  {  
    assert(val.K()>op_dxy);
    RN_ fxy(val('.',0,op_dxy)); 
  
    fxy[0] = 4*Dl0.x*Dl0.y;
    fxy[1] = 4*Dl1.x*Dl1.y;
    fxy[2] = 4*Dl2.x*Dl2.y;
    fxy[3] =  4*(Dl1.x*Dl2.y + Dl1.y*Dl2.x);
    fxy[4] =  4*(Dl0.x*Dl2.y + Dl0.y*Dl2.x);
    fxy[5] =  4*(Dl0.x*Dl1.y + Dl0.y*Dl1.x);
  }
 
 }
 
}

 
//
// end ttdc
// ------------------

static TypeOfFE_RTortho The_TypeOfFE_RTortho;
static TypeOfFE_RT The_TypeOfFE_RT;
static TypeOfFE_P0 The_TypeOfFE_P0;
static TypeOfFE_P1ttdc The_TypeOfFE_P1ttdc;
static TypeOfFE_P2ttdc The_TypeOfFE_P2ttdc;
static TypeOfFE_RTmodif The_TypeOfFE_RTmodif;
static TypeOfFE_P1ncLagrange The_TypeOfFE_P1nc;
static TypeOfFE_ConsEdge The_TypeOfFE_ConsEdge; // add FH 
TypeOfFE  & RTLagrangeOrtho(The_TypeOfFE_RTortho);
TypeOfFE  & RTLagrange(The_TypeOfFE_RT);
TypeOfFE  & RTmodifLagrange(The_TypeOfFE_RTmodif);
TypeOfFE  & P0Lagrange(The_TypeOfFE_P0);
TypeOfFE  & P1ncLagrange(The_TypeOfFE_P1nc);
TypeOfFE  & P1ttdc(The_TypeOfFE_P1ttdc);
TypeOfFE  & P2ttdc(The_TypeOfFE_P2ttdc);
TypeOfFE  & P0edge(The_TypeOfFE_ConsEdge);

}
