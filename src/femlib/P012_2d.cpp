// ********** DO NOT REMOVE THIS BANNER **********
// ORIG-DATE:     Jan 2008
// -*- Mode : c++ -*-
//
// SUMMARY  : P0,P1,P2 lagrange 2D 
// USAGE    : LGPL      
// ORG      : LJLL Universite Pierre et Marie Curi, Paris,  FRANCE 
// AUTHOR   : Frederic Hecht
// E-MAIL   : frederic.hecht@ann.jussieu.fr
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

 Thank to the ARN   FF2A3 grant
 ref:ANR-07-CIS7-002-01 
 */


#include "PkLagrange.hpp"

// P0 

 namespace Fem2D {
class TypeOfFE_P0Lagrange2d : public  TypeOfFE_Lagrange<Mesh2> { 
public:  
  TypeOfFE_P0Lagrange2d(): TypeOfFE_Lagrange<Mesh2>(0) {  }
  void FB(const What_d whatd,const Mesh & Th,const Element & K,const Rd &P, RNMK_ & val) const;
  virtual R operator()(const FElement & K,const  RdHat & PHat,const KN_<R> & u,int componante,int op) const ;
  
} ;



R TypeOfFE_P0Lagrange2d::operator()(const FElement & K,const  R2 & PHat,const KN_<R> & u,int componante,int op) const 
{ 
  R u0(u(K(0)));
  R r=0;
  if (op==0)
    r = u0;
  else  r=0;
  return r;
}

void TypeOfFE_P0Lagrange2d::FB(const What_d whatd,const Mesh & ,const Element & K,const R2 & P,RNMK_ & val) const
{
  assert(val.N() >=Element::nv);
  assert(val.M()==1 );
  
  val=0; 
  RN_ f0(val('.',0,op_id)); 
  
  if (whatd & Fop_D0) 
    {
      f0[0] = 1;
    }
}


// P1 
class TypeOfFE_P1Lagrange2d : public  TypeOfFE_Lagrange<Mesh2> { 
public:  
  TypeOfFE_P1Lagrange2d(): TypeOfFE_Lagrange<Mesh2>(1) {  }
  void FB(const What_d whatd,const Mesh & Th,const Element & K,const Rd &P, RNMK_ & val) const;
  virtual R operator()(const FElement & K,const  RdHat & PHat,const KN_<R> & u,int componante,int op) const ;
  
} ;



R TypeOfFE_P1Lagrange2d::operator()(const FElement & K,const  R2 & PHat,const KN_<R> & u,int componante,int op) const 
{ 
  R u0(u(K(0))), u1(u(K(1))), u2(u(K(2)));
  R r=0;
  if (op==0)
    {
      R l0=1-PHat.x-PHat.y,l1=PHat.x,l2=PHat.y; 
      r = u0*l0+u1*l1+l2*u2;
    }
  else if(op==op_dx || op==op_dx )
    { 
      const Element & T=K.T;
      R2 D[3];
      T.Gradlambda(D);
      if (op==op_dx)
	r =  D[0].x*u0 + D[1].x*u1 + D[2].x*u2;
      else 
	r =  D[0].y*u0 + D[1].y*u1 + D[2].y*u2 ;
    }
  //  cout << r << "\t";
  return r;
}

void TypeOfFE_P1Lagrange2d::FB(const What_d whatd,const Mesh & ,const Element & K,const R2 & P,RNMK_ & val) const
{
  //  const Triangle & K(FE.T);
  
  R l[]={1.-P.sum(),P.x,P.y}; 
  
  assert(val.N() >=Element::nv);
  assert(val.M()==1 );
  
  val=0; 
  RN_ f0(val('.',0,op_id)); 
  
  if (whatd & Fop_D0) 
    {
      f0[0] = l[0];
      f0[1] = l[1];
      f0[2] = l[2];
    }
  if (whatd & Fop_D1)
    {
      R2 Dl[3];
      K.Gradlambda(Dl);
      if (whatd & Fop_dx) 
	{
	  RN_ f0x(val('.',0,op_dx)); 
	  f0x[0] = Dl[0].x;
	  f0x[1] = Dl[1].x;
	  f0x[2] = Dl[2].x;
	  
	}
      
      if (whatd & Fop_dy) {
	RN_ f0y(val('.',0,op_dy)); 
	f0y[0] = Dl[0].y;
	f0y[1] = Dl[1].y;
	f0y[2] = Dl[2].y;
      }

    }
}




class TypeOfFE_P2Lagrange2d : public  TypeOfFE_Lagrange<Mesh2>  { 
public:  
  typedef Mesh2 Mesh;
  typedef GFElement<Mesh2> FElement;
  TypeOfFE_P2Lagrange2d(): TypeOfFE_Lagrange<Mesh2>(2) {  }
  void FB(const What_d whatd,const Mesh & Th,const Mesh2::Element & K,const Rd &P, RNMK_ & val) const;
  
} ;



void TypeOfFE_P2Lagrange2d::FB(const What_d whatd,const Mesh & ,const Element & K,const R2 & P,RNMK_ & val) const
{
  //  const Triangle & K(FE.T);
  
  R l[]={1.-P.sum(),P.x,P.y}; 
  
  assert(val.N() >=E::nv+E::ne);
  assert(val.M()==1 );
  
  val=0; 
  RN_ f0(val('.',0,op_id)); 
  //  
  if (whatd & Fop_D0) 
    {
	int k=0;
	for(int i=0;i<E::nv;++i)
          f0[k++] = l[i]*(2*l[i]-1.);
 	for(int i=0;i<E::ne;++i)
	  f0[k++] = 4.*l[E::nvedge[i][0]]*l[E::nvedge[i][1]];  
    }
  
  if (whatd & (Fop_D1|Fop_D2))
  {
      R2 Dl[3];
      R l4[3]={ (4*l[0]-1),(4*l[1]-1),(4*l[2]-1)}; 
    
      K.Gradlambda(Dl);
      RN_ f0x(val('.',0,op_dx));
      RN_ f0y(val('.',0,op_dy)); 
      RN_ f0z(val('.',0,op_dz)); 
      int k=0;
      for(int i=0;i<E::nv;++i,++k)
      {
	  f0x[k] = Dl[i].x*l4[i];
	  f0y[k] = Dl[i].y*l4[i];
      }
      for(int i=0;i<E::ne;++i,++k)
      {
	  int i0=E::nvedge[i][0],i1=E::nvedge[i][1];
	  f0x[k] = 4*(Dl[i1].x*l[i0] + Dl[i0].x*l[i1]) ;
	  f0y[k] = 4*(Dl[i1].y*l[i0] + Dl[i0].y*l[i1]) ;
      }
      assert(k==6);
      
      //cout << " D2 " << whatd <<  endl;
      if (whatd & Fop_D2)
      {
	  //cout << " D2 " << endl;
	  RN_ f0xx(val('.',0,op_dxx));
	  RN_ f0yy(val('.',0,op_dyy)); 
	  RN_ f0xy(val('.',0,op_dxy));
	  
	   k=0;
	  for(int i=0;i<E::nv;++i,++k)
	  {
	      f0xx[k] = 4.*Dl[i].x*Dl[i].x;
	      f0yy[k] = 4.*Dl[i].y*Dl[i].y;
	      f0xy[k] = 4.*Dl[i].x*Dl[i].y;
	  }
	  for(int i=0;i<E::ne;++i,++k)
	  {
	      int i0=E::nvedge[i][0],i1=E::nvedge[i][1];
	      f0xx[k] = 8.*Dl[i0].x*Dl[i1].x;
	      f0yy[k] = 8.*Dl[i0].y*Dl[i1].y;
	      f0xy[k] = 4.*(Dl[i0].x*Dl[i1].y+ Dl[i1].x*Dl[i0].y);
	  } 
	  assert(k==6);
      }
      
  } 
    
}



static TypeOfFE_P0Lagrange2d  P0_2d;
GTypeOfFE<Mesh2> & P0Lagrange2d(P0_2d);

static TypeOfFE_P1Lagrange2d  P1_2d;
GTypeOfFE<Mesh2> & P1Lagrange2d(P1_2d);

static TypeOfFE_P2Lagrange2d  P2_2d;
GTypeOfFE<Mesh2> & P2Lagrange2d(P2_2d);




template<> GTypeOfFE<Mesh2> & DataFE<Mesh2>::P0=P0_2d; 
template<> GTypeOfFE<Mesh2> & DataFE<Mesh2>::P1=P1_2d; 
template<> GTypeOfFE<Mesh2> & DataFE<Mesh2>::P2=P2_2d; 


 }
