// ********** DO NOT REMOVE THIS BANNER **********
// ORIG-DATE:     Jan 2008
// -*- Mode : c++ -*-
//
// SUMMARY  : P0,P1,P2 lagrange 3D 
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

 namespace Fem2D {
// P0 


class TypeOfFE_P0Lagrange3d : public  TypeOfFE_Lagrange<Mesh3> { 
public:  
  TypeOfFE_P0Lagrange3d(): TypeOfFE_Lagrange<Mesh3>(0) {  }
  void FB(const What_d whatd,const Mesh & Th,const Element & K,const Rd &P, RNMK_ & val) const;
  virtual R operator()(const FElement & K,const  RdHat & PHat,const KN_<R> & u,int componante,int op) const ;
  
} ;

     


R TypeOfFE_P0Lagrange3d::operator()(const FElement & K,const  R3 & PHat,const KN_<R> & u,int componante,int op) const 
{ 
  R u0(u(K(0)));
  R r=0;
  if (op==0)
    r = u0;
  else  r=0;
  return r;
}

void TypeOfFE_P0Lagrange3d::FB(const What_d whatd,const Mesh & ,const Element & K,const R3 & P,RNMK_ & val) const
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
class TypeOfFE_P1Lagrange3d : public  TypeOfFE_Lagrange<Mesh3> { 
public:  
  TypeOfFE_P1Lagrange3d(): TypeOfFE_Lagrange<Mesh3>(1) {  }
  void FB(const What_d whatd,const Mesh & Th,const Element & K,const Rd &P, RNMK_ & val) const;
  virtual R operator()(const FElement & K,const  RdHat & PHat,const KN_<R> & u,int componante,int op) const ;
  
} ;



R TypeOfFE_P1Lagrange3d::operator()(const FElement & K,const  R3 & PHat,const KN_<R> & u,int componante,int op) const 
{ 
  R u0(u(K(0))), u1(u(K(1))), u2(u(K(2))),u3(u(K(3)));
  R r=0;
  if (op==0)
    {
      R l0=1-PHat.x-PHat.y-PHat.z,l1=PHat.x,l2=PHat.y,l3=PHat.z; 
      r = u0*l0+u1*l1+l2*u2+l3*u3;
    }
  else if(op==op_dx || op==op_dy || op==op_dz) // dx => dy thank to Pichon 27/01/2008 (FH)
    { 
      const Element & T=K.T;
      R3 D[4];
      T.Gradlambda(D);
      if (op==op_dx)
	r =  D[0].x*u0 + D[1].x*u1 + D[2].x*u2+ D[3].x*u3 ;
      else if (op==op_dy) 
	r =  D[0].y*u0 + D[1].y*u1 + D[2].y*u2+ D[3].y*u3 ;
      else 
	r =  D[0].z*u0 + D[1].z*u1 + D[2].z*u2+ D[3].z*u3 ;
    }
  //  cout << r << "\t";
  return r;
}

void TypeOfFE_P1Lagrange3d::FB(const What_d whatd,const Mesh & ,const Element & K,const R3 & P,RNMK_ & val) const
{
  //  const Triangle & K(FE.T);
  
  R l[]={1.-P.sum(),P.x,P.y,P.z}; 
  
  assert(val.N() >=Element::nv);
  assert(val.M()==1 );
  
  val=0; 
  RN_ f0(val('.',0,op_id)); 
  
  if (whatd & Fop_D0) 
    {
      f0[0] = l[0];
      f0[1] = l[1];
      f0[2] = l[2];
      f0[3] = l[3];
    }
  if (whatd & Fop_D1)
    {
      R3 Dl[4];
      K.Gradlambda(Dl);
      //for(int i=0;i<4;++i)
	//      cout << Dl[i] << endl;
      if (whatd & Fop_dx) 
	{
	  RN_ f0x(val('.',0,op_dx)); 
	  f0x[0] = Dl[0].x;
	  f0x[1] = Dl[1].x;
	  f0x[2] = Dl[2].x;
	  f0x[3] = Dl[3].x;
	  
	}
      
      if (whatd & Fop_dy) {
	RN_ f0y(val('.',0,op_dy)); 
	f0y[0] = Dl[0].y;
	f0y[1] = Dl[1].y;
	f0y[2] = Dl[2].y;
	f0y[3] = Dl[3].y;
      }

      if (whatd & Fop_dz) {
	RN_ f0z(val('.',0,op_dz)); 
	f0z[0] = Dl[0].z;
	f0z[1] = Dl[1].z;
	f0z[2] = Dl[2].z;
	f0z[3] = Dl[3].z;
      }
    }
  //  cout << val << endl;
}




class TypeOfFE_P2Lagrange3d : public TypeOfFE_Lagrange<Mesh3>  { 
public:  
  typedef Mesh3 Mesh;
  typedef GFElement<Mesh3> FElement;
  TypeOfFE_P2Lagrange3d(): TypeOfFE_Lagrange<Mesh3>(2) {  }
  void FB(const What_d whatd,const Mesh & Th,const Mesh3::Element & K,const Rd &P, RNMK_ & val) const;
} ;

class TypeOfFE_P1bLagrange3d : public TypeOfFE_Lagrange<Mesh3>  { 
    public:  
	 typedef Mesh3 Mesh;
	 typedef GFElement<Mesh3> FElement;
	 TypeOfFE_P1bLagrange3d(): TypeOfFE_Lagrange<Mesh3>(-1) {  }
	 void FB(const What_d whatd,const Mesh & Th,const Mesh3::Element & K,const Rd &P, RNMK_ & val) const;
} ;
     
     
      


void TypeOfFE_P2Lagrange3d::FB(const What_d whatd,const Mesh & ,const Element & K,const R3 & P,RNMK_ & val) const
{
  //  const Triangle & K(FE.T);
  
  R l[]={1.-P.sum(),P.x,P.y,P.z}; 
  
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
      R3 Dl[4];
      R l4[4]={ (4*l[0]-1),(4*l[1]-1),(4*l[2]-1),(4*l[3]-1)}; 
      
      K.Gradlambda(Dl);
      RN_ f0x(val('.',0,op_dx));
      RN_ f0y(val('.',0,op_dy)); 
      RN_ f0z(val('.',0,op_dz)); 
      int k=0;
      for(int i=0;i<E::nv;++i,++k)
      {
	  f0x[k] = Dl[i].x*l4[i];
	  f0y[k] = Dl[i].y*l4[i];
	  f0z[k] = Dl[i].z*l4[i];
      }
      for(int i=0;i<E::ne;++i,++k)
      {
	  int i0=E::nvedge[i][0],i1=E::nvedge[i][1];
	  f0x[k] = 4*(Dl[i1].x*l[i0] + Dl[i0].x*l[i1]) ;
	  f0y[k] = 4*(Dl[i1].y*l[i0] + Dl[i0].y*l[i1]) ;
	  f0z[k] = 4*(Dl[i1].z*l[i0] + Dl[i0].z*l[i1]) ;
      }
      assert(k==10);
      
      //cout << " D2 " << whatd <<  endl;
      if (whatd & Fop_D2)
      {
	  //cout << " D2 " << endl;
	  RN_ f0xx(val('.',0,op_dxx));
	  RN_ f0yy(val('.',0,op_dyy)); 
	  RN_ f0zz(val('.',0,op_dzz)); 
	  RN_ f0xy(val('.',0,op_dxy));
	  RN_ f0xz(val('.',0,op_dxz));
	  RN_ f0yz(val('.',0,op_dyz));
	  
	   k=0;
	  for(int i=0;i<E::nv;++i,++k)
	  {
	      f0xx[k] = 4.*Dl[i].x*Dl[i].x;
	      f0yy[k] = 4.*Dl[i].y*Dl[i].y;
	      f0zz[k] = 4.*Dl[i].z*Dl[i].z;
	      f0xy[k] = 4.*Dl[i].x*Dl[i].y;
	      f0xz[k] = 4.*Dl[i].x*Dl[i].z;
	      f0yz[k] = 4.*Dl[i].y*Dl[i].z;
	  }
	  for(int i=0;i<E::ne;++i,++k)
	  {
	      int i0=E::nvedge[i][0],i1=E::nvedge[i][1];
	      f0xx[k] = 8.*Dl[i0].x*Dl[i1].x;
	      f0yy[k] = 8.*Dl[i0].y*Dl[i1].y;
	      f0zz[k] = 8.*Dl[i0].z*Dl[i1].z;
	      f0xy[k] = 4.*(Dl[i0].x*Dl[i1].y+ Dl[i1].x*Dl[i0].y);
	      f0xz[k] = 4.*(Dl[i0].x*Dl[i1].z+ Dl[i1].x*Dl[i0].z);
	      f0yz[k] = 4.*(Dl[i0].y*Dl[i1].z+ Dl[i1].y*Dl[i0].z);
	  } 
	  assert(k==10);
      }
      
  } 
/* 
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
 */
    
}

     /*
 R TypeOfFE_P1bLagrange3d::operator()(const FElement & K,const  R3 & PHat,const KN_<R> & u,int componante,int op) const 
     { 
	 R u0(u(K(0))), u1(u(K(1))), u2(u(K(2))),u3(u(K(3))),u4(u(K(4)));
	 R r=0;
	 if (op==0)
	   {
	       R l0=1-PHat.x-PHat.y-PHat.z,l1=PHat.x,l2=PHat.y,l3=PHat.z; 
	       R l0123=
	       r = u0*l0+u1*l1+l2*u2+l3*u3;
	   }
	 else if(op==op_dx || op==op_dy || op==op_dz) // dx => dy thank to Pichon 27/01/2008 (FH)
	   { 
	       const Element & T=K.T;
	       R3 D[4];
	       T.Gradlambda(D);
	       if (op==op_dx)
		   r =  D[0].x*u0 + D[1].x*u1 + D[2].x*u2+ D[3].x*u3 ;
	       else if (op==op_dy) 
		   r =  D[0].y*u0 + D[1].y*u1 + D[2].y*u2+ D[3].y*u3 ;
	       else 
		   r =  D[0].z*u0 + D[1].z*u1 + D[2].z*u2+ D[3].z*u3 ;
	   }
	 //  cout << r << "\t";
	 return r;
     }
 */    
     void TypeOfFE_P1bLagrange3d::FB(const What_d whatd,const Mesh & ,const Element & K,const R3 & P,RNMK_ & val) const
     {
	 //  const Triangle & K(FE.T);
	 const R d1=d+1.;
	 const R d13=d1*d1*d1;
	 const R d14=d13*d1;
	 R ll[]={1.-P.sum(),P.x,P.y,P.z}; 
	 R lb4= (ll[0]*ll[1]*ll[2]*ll[3])*d13; // d1^-4 d1^3 = 1/d1 in G
	 R lb=lb4*d1; // 1  in G     
	 R l[5];
	 for(int i=0;i<4;i++)
	    l[i]=ll[i]-lb4;  //  1/d1 in G - 1/d1 G =0 
	 l[5]=lb;
	 
	 assert(val.N() >=Element::nv);
	 assert(val.M()==1 );
	 
	 val=0; 
	 RN_ f0(val('.',0,op_id)); 
	 
	 if (whatd & Fop_D0) 
	   {
	       f0[0] = l[0];
	       f0[1] = l[1];
	       f0[2] = l[2];
	       f0[3] = l[3];
	       f0[4] = l[4];
	   }
	 if (whatd & Fop_D1)
	   {
	       R3 Dl[4];
	       K.Gradlambda(Dl);
	       R3 Dlb4 = (
			 + Dl[0]*(ll[1]*ll[2]*ll[3])
			 + Dl[1]*(ll[0]*ll[2]*ll[3])
			 + Dl[2]*(ll[0]*ll[1]*ll[3])
			 + Dl[3]*(ll[0]*ll[1]*ll[2]) )*d13;
	       
	       //for(int i=0;i<4;++i)
	       //      cout << Dl[i] << endl;
	      
	       if (whatd & Fop_dx) 
		 {
		     RN_ f0x(val('.',0,op_dx)); 
		     f0x[0] = Dl[0].x-Dlb4.x;
		     f0x[1] = Dl[1].x-Dlb4.x;
		     f0x[2] = Dl[2].x-Dlb4.x;
		     f0x[3] = Dl[3].x-Dlb4.x;
		     f0x[4] = Dlb4.x*d1;
		     
		 }
	       
	       if (whatd & Fop_dy) {
		   RN_ f0y(val('.',0,op_dy)); 
		   f0y[0] = Dl[0].y-Dlb4.y;
		   f0y[1] = Dl[1].y-Dlb4.y;
		   f0y[2] = Dl[2].y-Dlb4.y;
		   f0y[3] = Dl[3].y-Dlb4.y;
		   f0y[4] = Dlb4.y*d1;
	       }
	       
	       if (whatd & Fop_dz) {
		   RN_ f0z(val('.',0,op_dz)); 
		   f0z[0] = Dl[0].z-Dlb4.z;
		   f0z[1] = Dl[1].z-Dlb4.z;
		   f0z[2] = Dl[2].z-Dlb4.z;
		   f0z[3] = Dl[3].z-Dlb4.z;
		   f0z[4] = Dlb4.z*d1;
		   
	       }
	   }
	 else if (whatd & Fop_D2)
	     ffassert(0); // a faire ...
	 //  cout << val << endl;
     }
     
     

static TypeOfFE_P0Lagrange3d  P0_3d;
GTypeOfFE<Mesh3> & P0Lagrange3d(P0_3d);

static TypeOfFE_P1Lagrange3d  P1_3d;
GTypeOfFE<Mesh3> & P1Lagrange3d(P1_3d);

static TypeOfFE_P2Lagrange3d  P2_3d;
GTypeOfFE<Mesh3> & P2Lagrange3d(P2_3d);

static TypeOfFE_P1bLagrange3d  P1b_3d;
GTypeOfFE<Mesh3> & P1bLagrange3d(P1b_3d);




template<> GTypeOfFE<Mesh3> & DataFE<Mesh3>::P0=P0_3d; 
template<> GTypeOfFE<Mesh3> & DataFE<Mesh3>::P1=P1_3d; 
template<> GTypeOfFE<Mesh3> & DataFE<Mesh3>::P2=P2_3d; 


 }

