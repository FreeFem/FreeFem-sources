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
  void FB(const What_d whatd,const Mesh & Th,const Element & K,const RdHat &PHat, RNMK_ & val) const;
  virtual R operator()(const FElement & K,const  RdHat & PHat,const KN_<R> & u,int componante,int op) const ;
  
} ;

     


R TypeOfFE_P0Lagrange3d::operator()(const FElement & K,const  RdHat & PHat,const KN_<R> & u,int componante,int op) const
{ 
  R u0(u(K(0)));
  R r=0;
  if (op==0)
    r = u0;
  else  r=0;
  return r;
}

void TypeOfFE_P0Lagrange3d::FB(const What_d whatd,const Mesh & ,const Element & K,const RdHat & PHat,RNMK_ & val) const
{
  assert(val.N() >=1);
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
  void FB(const What_d whatd,const Mesh & Th,const Element & K,const RdHat &PHat, RNMK_ & val) const;
  virtual R operator()(const FElement & K,const  RdHat & PHat,const KN_<R> & u,int componante,int op) const ;
  
} ;

 class TypeOfFE_PkdcLagrange3d : public  TypeOfFE_LagrangeDC<Mesh3>  {
 public:
     TypeOfFE_PkdcLagrange3d(int k,double skrink=0.0001): TypeOfFE_LagrangeDC<Mesh3>(k,skrink) {  }

 };
 
 
R TypeOfFE_P1Lagrange3d::operator()(const FElement & K,const  RdHat & PHat,const KN_<R> & u,int componante,int op) const
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

void TypeOfFE_P1Lagrange3d::FB(const What_d whatd,const Mesh & ,const Element & K,const RdHat & PHat,RNMK_ & val) const
{
  //  const Triangle & K(FE.T);
  
  R l[]={1.-PHat.sum(),PHat.x,PHat.y,PHat.z};
  
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
  void FB(const What_d whatd,const Mesh & Th,const Mesh3::Element & K,const RdHat &PHat, RNMK_ & val) const;
} ;

class TypeOfFE_P1bLagrange3d : public TypeOfFE_Lagrange<Mesh3>  { 
    public:  
	 typedef Mesh3 Mesh;
	 typedef GFElement<Mesh3> FElement;
	 TypeOfFE_P1bLagrange3d(): TypeOfFE_Lagrange<Mesh3>(-1) {  }
	 void FB(const What_d whatd,const Mesh & Th,const Mesh3::Element & K,const RdHat &PHat, RNMK_ & val) const;
} ;
     

     

void TypeOfFE_P2Lagrange3d::FB(const What_d whatd,const Mesh & ,const Element & K,const RdHat & PHat,RNMK_ & val) const
{
  //  const Triangle & K(FE.T);
  
  R l[]={1.-PHat.sum(),PHat.x,PHat.y,PHat.z};
  
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

    if( whatd & Fop_dx)
      {
	RN_ f0x(val('.',0,op_dx));
 
	int k=0;
	for(int i=0;i<E::nv;++i,++k)
	  {
	    f0x[k] = Dl[i].x*l4[i];
	  }
	for(int i=0;i<E::ne;++i,++k)
	  {
	    int i0=E::nvedge[i][0],i1=E::nvedge[i][1];
	    f0x[k] = 4*(Dl[i1].x*l[i0] + Dl[i0].x*l[i1]) ;
	  }
	assert(k==10);
      }

    if( whatd & Fop_dy)
      {
	RN_ f0y(val('.',0,op_dy));
 
	int k=0;
	for(int i=0;i<E::nv;++i,++k)
	  {
	    f0y[k] = Dl[i].y*l4[i];
	  }
	for(int i=0;i<E::ne;++i,++k)
	  {
	    int i0=E::nvedge[i][0],i1=E::nvedge[i][1];
	    f0y[k] = 4*(Dl[i1].y*l[i0] + Dl[i0].y*l[i1]) ;
	  }
	assert(k==10);
      }

    if( whatd & Fop_dz)
      {
	RN_ f0z(val('.',0,op_dz));
 
	int k=0;
	for(int i=0;i<E::nv;++i,++k)
	  {
	    f0z[k] = Dl[i].z*l4[i];
	  }
	for(int i=0;i<E::ne;++i,++k)
	  {
	    int i0=E::nvedge[i][0],i1=E::nvedge[i][1];
	    f0z[k] = 4*(Dl[i1].z*l[i0] + Dl[i0].z*l[i1]) ;
	  }
	assert(k==10);
      }

    /* avant
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
    */
    //cout << " D2 " << whatd <<  endl;
    if (whatd & Fop_D2)
      {
	
	//cout << " D2 " << endl;
	if (whatd & Fop_dxx){
	  RN_ f0xx(val('.',0,op_dxx));
	
	  int k=0;
	  for(int i=0;i<E::nv;++i,++k)
	    {
	      f0xx[k] = 4.*Dl[i].x*Dl[i].x;
	    }
	  for(int i=0;i<E::ne;++i,++k)
	    {
	      int i0=E::nvedge[i][0],i1=E::nvedge[i][1];
	      f0xx[k] = 8.*Dl[i0].x*Dl[i1].x;
	    } 
	  assert(k==10);
	}
	//cout << " D2 " << endl;
	if (whatd & Fop_dyy){
	  RN_ f0yy(val('.',0,op_dyy)); 
	
	  int k=0;
	  for(int i=0;i<E::nv;++i,++k)
	    {	
	      f0yy[k] = 4.*Dl[i].y*Dl[i].y;
	    }
	  for(int i=0;i<E::ne;++i,++k)
	    {
	      int i0=E::nvedge[i][0],i1=E::nvedge[i][1];
	      f0yy[k] = 8.*Dl[i0].y*Dl[i1].y;
	    }
	  assert(k==10);
	}
	//cout << " D2 " << endl;
	if (whatd & Fop_dzz){	 
	  RN_ f0zz(val('.',0,op_dzz)); 
      
	  int k=0;
	  for(int i=0;i<E::nv;++i,++k)
	    {
	      f0zz[k] = 4.*Dl[i].z*Dl[i].z;
	    }
	  for(int i=0;i<E::ne;++i,++k)
	    {
	      int i0=E::nvedge[i][0],i1=E::nvedge[i][1];
	      f0zz[k] = 8.*Dl[i0].z*Dl[i1].z;
	      
	    }
	  assert(k==10);
	}
	//cout << " D2 " << endl;
	if (whatd & Fop_dxy){		  
	  RN_ f0xy(val('.',0,op_dxy));
	 
	  int k=0;
	  for(int i=0;i<E::nv;++i,++k)
	    {
	      f0xy[k] = 4.*Dl[i].x*Dl[i].y;
	    }
	  for(int i=0;i<E::ne;++i,++k)
	    {
	      int i0=E::nvedge[i][0],i1=E::nvedge[i][1];
	      f0xy[k] = 4.*(Dl[i0].x*Dl[i1].y+ Dl[i1].x*Dl[i0].y);
	    } 
	  assert(k==10);
	}

	//cout << " D2 " << endl;
	if (whatd & Fop_dxz){
	  RN_ f0xz(val('.',0,op_dxz));
	 
	  int k=0;
	  for(int i=0;i<E::nv;++i,++k)
	    {
	      f0xz[k] = 4.*Dl[i].x*Dl[i].z;
	    }
	  for(int i=0;i<E::ne;++i,++k)
	    {
	      int i0=E::nvedge[i][0],i1=E::nvedge[i][1];
	      f0xz[k] = 4.*(Dl[i0].x*Dl[i1].z+ Dl[i1].x*Dl[i0].z);
	    } 
	  assert(k==10);
	}
	//cout << " D2 " << endl;
	if (whatd & Fop_dyz){

	  RN_ f0yz(val('.',0,op_dyz));
	  
	  int k=0;
	  for(int i=0;i<E::nv;++i,++k)
	    {
	      f0yz[k] = 4.*Dl[i].y*Dl[i].z;
	    }
	  for(int i=0;i<E::ne;++i,++k)
	    {
	      int i0=E::nvedge[i][0],i1=E::nvedge[i][1];
	      f0yz[k] = 4.*(Dl[i0].y*Dl[i1].z+ Dl[i1].y*Dl[i0].z);
	    } 
	  assert(k==10);
	}
	/*
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
      
	*/
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
     void TypeOfFE_P1bLagrange3d::FB(const What_d whatd,const Mesh & ,const Element & K,const RdHat & PHat,RNMK_ & val) const
     {
	 //  const Triangle & K(FE.T);
	 const R d1=dHat+1.;
	 const R d13=d1*d1*d1;
	 const R d14=d13*d1;
	 R ll[]={1.-PHat.sum(),PHat.x,PHat.y,PHat.z};
	 R lb4= (ll[0]*ll[1]*ll[2]*ll[3])*d13; // d1^-4 d1^3 = 1/d1 in G
	 R lb=lb4*d1; // 1  in G     
	 R l[5];
	 for(int i=0;i<4;i++)
	    l[i]=ll[i]-lb4;  //  1/d1 in G - 1/d1 G =0 
	 l[4]=lb;
	 
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

     
     
     class TypeOfFE_RT0_3d : public    GTypeOfFE<Mesh3>  { 
     public:  
	 typedef Mesh3 Mesh;
	 typedef  Mesh3::Element  Element;
	 
	 typedef GFElement<Mesh3> FElement;
	 static int dfon[];
	 static const int d=Mesh::Rd::d;
	 TypeOfFE_RT0_3d();
	 int edgeface[4][3] ;
	void FB(const What_d whatd,const Mesh & Th,const Mesh3::Element & K,const RdHat &PHat, RNMK_ & val) const;
	void  set(const Mesh & Th,const Element & K,InterpolationMatrix<RdHat> & M,int ocoef,int odf,int *nump  ) const;
     } ;
     int TypeOfFE_RT0_3d::dfon[]={0,0,1,0}; 
     

     TypeOfFE_RT0_3d::TypeOfFE_RT0_3d(): GTypeOfFE<Mesh3>(dfon,d,1,3*4,4,false,true)
     { 	
       //  integration on middle of faces  (light ) on  each face ..
	R3 Pt[]= {R3(0.,0.,0.), R3(1.,0.,0.),R3(0.,1.,0.),R3(0.,0.,1.)};
         /*
	for (int i=0;i<Element::ne;++i)
          this->PtInterpolation[i]=(Pt[Element::nvedge[i][0]]+Pt[Element::nvedge[i][1]])*0.5;
	 */
         for (int i=0;i<Element::nf;++i)
             this->PtInterpolation[i]=(Pt[Element::nvface[i][0]]+Pt[Element::nvface[i][1]]+Pt[Element::nvface[i][2]])/3.;

//	 static const int  nvfaceTet[4][3]  ={{3,2,1}, {0,2,3},{ 3,1,0},{ 0,1,2}}  ;//{ {2,1,3},{0,2,3},{1,0,3},{0,1,2} };
//	    { {0,1},{0,2},{0,3},{1,2},{1,3},{2,3} };
	 //     0     1     2     3     4     5
     /*  {
	   int i=0;
	   for (int f=0;f<4;f++)
	       for (int e=0,i=0;e<6;e++)
		   if ((Element::nvedge[e][0] !=f)  && (Element::nvedge[e][1]!=f))
		       edgeface[f][i++]=e; 
      }*/
       {
	 int i=0;
	 for (int f=0;f<4;f++) 
	   {
	     //  cout << " face : " << f << endl;
	   //  for (int p=0;p<3;p++)
	     {
                 int e= f; //edgeface[f][p] ;
		// cout << "  , "  << this->PtInterpolation[e];
		 for (int c=0;c<3;c++,i++) 

	           {
	             this->pInterpolation[i]=e;
	             this->cInterpolation[i]=c;
		     this->dofInterpolation[i]=f;
	             this->coefInterpolation[i]=0.;	       
	           }
	     }
	   //cout <<  endl;
	   }
       }
       
     }
     void  TypeOfFE_RT0_3d::set(const Mesh & Th,const Element & K,InterpolationMatrix<RdHat> & M ,int ocoef,int odf,int *nump) const
     {
	 //   compute de coef d'interpolation
	// M.coef 
	 int i=ocoef;
	 for (int f=0;f<4;f++) 
	   {
	       R3 N=K.N(f);//  exterior and  ||N|| = 2* area f
	       N *= K.faceOrient(f)/2.;
	     //  for (int p=0;p<3;p++)
		 {
                     int e= f; //dgeface[f][p] ;
		     for (int c=0;c<3;c++,i++) 
			 
		       {
			   //this->pInterpolation[i]=e;
			   //this->cInterpolation[i]=c;
			   //this->dofInterpolation[i]=f;
			   M.coef[i]=N[c];	       
		       }
		 }}
		// cout << " M.coef :" << M.coef << endl;
	 //ffassert(i==M.ncoef && M.np == 6 );
 	 
     }
     void  TypeOfFE_RT0_3d::FB(const What_d whatd,const Mesh & Th,const Mesh3::Element & K,const RdHat &PHat, RNMK_ & val) const
     {
	 assert(val.N() >=4);
	 assert(val.M()==3 );
	 // wi = signe * (x - qi)/ (volume*d)   
	 val=0; 
        // cout << " TypeOfFE_RT0_3d "<< Th(K) << " " << Th.nt << " / "<<K.faceOrient(0)<< " " << K.faceOrient(1) << " " << K.faceOrient(2) << " " << K.faceOrient(3) <<endl;
	 R cc =1./(d*K.mesure());
	 R ci[4]={ cc*K.faceOrient(0),cc*K.faceOrient(1),cc*K.faceOrient(2),cc*K.faceOrient(3)};
	 
	 if (whatd & Fop_D0) 
	   {
	       R3 X=K(PHat);
	       int k=0;
	       for(int i=0;i<4;++i)
		 { R3 wi=(X-K[i])*ci[i];
		     val(i,0,op_id) = wi.x ;
		     val(i,1,op_id) = wi.y ;
		     val(i,2,op_id) = wi.z ;
		     //cout  << "RT 3d  "<<i << " "<< X << " " <<wi << " fo: " << K.faceOrient(i) <<endl;
		 }
	   }
	 
	 if (whatd & Fop_D1)
	   {
	       RN_ Ci(ci,4);
	       if (whatd & Fop_dx) 
		   val('.',0,op_dx) = Ci;
	       if (whatd & Fop_dy) 
		   val('.',1,op_dy) = Ci;	       
	       if (whatd & Fop_dz) 
		   val('.',2,op_dz) = Ci;
	        //cout  << "RT 3d  "<< val('.','.',op_dz) << endl;
	   }	 
	 
     }     
     

     class TypeOfFE_Edge0_3d  :public   GTypeOfFE<Mesh3>  { 
     public:  
	 typedef Mesh3 Mesh;
	 typedef  Mesh3::Element  Element;
	 
	 typedef GFElement<Mesh3> FElement;
	 static int dfon[];
	 static const int d=Mesh::Rd::d;
	 static const GQuadratureFormular<R1> QFe;
	 int edgeface[4][3] ;
	 TypeOfFE_Edge0_3d();
	 void FB(const What_d whatd,const Mesh & Th,const Mesh3::Element & K,const RdHat &PHat, RNMK_ & val) const;
	 void  set(const Mesh & Th,const Element & K,InterpolationMatrix<RdHat> & M ,int ocoef,int odf,int *nump  ) const;
     } ;
     int TypeOfFE_Edge0_3d::dfon[]={0,1,0,0}; 
     
     const  GQuadratureFormular<R1> TypeOfFE_Edge0_3d::QFe(-1+2*2,2,GaussLegendre(2),true);
     
     TypeOfFE_Edge0_3d::TypeOfFE_Edge0_3d(): GTypeOfFE<Mesh3>(TypeOfFE_Edge0_3d::dfon,d,1,Element::ne*3*QFe.n,Element::ne*QFe.n,false,true)
     { 	
	 assert(QFe.n);
	 //  integration on edge use QFe
	 R3 Pt[]= {R3(0.,0.,0.), R3(1.,0.,0.),R3(0.,1.,0.),R3(0.,0.,1.)};
	 for (int e=0,i=0;e<Element::ne;++e)
	     for(int q=0;q<QFe.n;++q,++i)
	       {
		   double x=QFe[q].x;
		   this->PtInterpolation[i]=Pt[Element::nvedge[e][0]]*x+Pt[Element::nvedge[e][1]]*(1-x);
	       }
       {
	   int i=0,p=0;
	   for (int e=0;e<Element::ne;e++) 
	     {
		 for(int q=0;q<QFe.n;++q,++p) 
		     for (int c=0;c<3;c++,i++) 
		       {
			   this->pInterpolation[i]=p;
			   this->cInterpolation[i]=c;
			   this->dofInterpolation[i]=e;
			   this->coefInterpolation[i]=0.;	       
		       }
		 
		 
	     }
       }
	// cout <<  " ++ TypeOfFE_Edge0_3d():"<< this->PtInterpolation << endl;
     }
     void  TypeOfFE_Edge0_3d::set(const Mesh & Th,const Element & K,InterpolationMatrix<RdHat> & M ,int ocoef,int odf,int *nump ) const
     {
	 // compute de coef d'interpolation
	 // M.coef 
	 int i=ocoef,p=0;
	 for (int e=0;e<Element::ne;++e) 
	   {
	       R3 E=K.EdgeOrientation(e)*K.Edge(e);//  exterior and  ||N|| = 2* area f
	       for(int q=0;q<QFe.n;++q,++p)  
		    for (int c=0;c<3;c++,i++) 
			     
			   {
			      // this->pInterpolation[i]=p;
			       //this->cInterpolation[i]=c;
			       //this->dofInterpolation[i]=e;
			       M.coef[i]=E[c]*QFe[q].a;	       
			   }
		 }
	 
	// ffassert(i==M.ncoef && M.np == p );
 	 
     }
     void  TypeOfFE_Edge0_3d::FB(const What_d whatd,const Mesh & Th,const Mesh3::Element & K,const RdHat &PHat, RNMK_ & val) const
     {
	 assert(val.N() >=6);
	 assert(val.M()==3 );
         R l[]={1.-PHat.sum(),PHat.x,PHat.y,PHat.z}; 
	 R3 D[4];
	 K.Gradlambda(D);
	 
	 // wi = signe * (x - qi)/ (volume*d)   
	 val=0; 
	 //  i,j : l1 grad lj - lj grad lj
	 // int_i^j  grad lj . t_ij = 1
	 
	 
	 int  se[]={ K.EdgeOrientation(0), K.EdgeOrientation(1), K.EdgeOrientation(2),
	 K.EdgeOrientation(3), K.EdgeOrientation(4), K.EdgeOrientation(5)};
	 
	 if (whatd & Fop_D0) 
	   {
	       R3 X=K(PHat);
	       int k=0;
	       for(int i=0;i<6;++i)
		 { 
		     int i0=Element::nvedge[i][0],i1=Element::nvedge[i][1];
		     if( se[i]<0) Exchange(i0,i1);
		     R3 wi = l[i0]*D[i1]-l[i1]*D[i0];
		     val(i,0,op_id) = wi.x ;
		     val(i,1,op_id) = wi.y ;
		     val(i,2,op_id) = wi.z ;
		    // cout  << "Edge0 3d  "<<i << " "<< X << " " <<wi << " fo: " <<se[i] <<endl;
		 }
	   }
	 
	 if (whatd & Fop_D1)
	     
	     for(int i=0;i<6;++i)
	       { 
		   int i0=Element::nvedge[i][0],i1=Element::nvedge[i][1];
		   if( se[i]<0) Exchange(i0,i1);
		   if (whatd & Fop_dx) 
		     {
			 R3 wi = D[i0].x*D[i1]-D[i1].x*D[i0];		    
			 val(i,0,op_dx) = wi.x ;
			 val(i,1,op_dx) = wi.y ;
			 val(i,2,op_dx) = wi.z ;
		     }
		   if (whatd & Fop_dy) 
		     {
			 R3 wi = D[i0].y*D[i1]-D[i1].y*D[i0];		    
			 val(i,0,op_dy) = wi.x ;
			 val(i,1,op_dy) = wi.y ;
			 val(i,2,op_dy) = wi.z ;
		     }
		   if (whatd & Fop_dz) 
		     {
			 R3 wi = D[i0].z*D[i1]-D[i1].z*D[i0];		    
			 val(i,0,op_dz) = wi.x ;
			 val(i,1,op_dz) = wi.y ;
			 val(i,2,op_dz) = wi.z ;
		     }
		   
	       }	 
	 
     }     
 class TypeOfFE_P0Face : public TypeOfFE_Lagrange<Mesh3>  {
 public:
   typedef Mesh3 Mesh;
   typedef GFElement<Mesh3> FElement;
     TypeOfFE_P0Face(): TypeOfFE_Lagrange<Mesh3>(2) {  }
   void FB(const What_d whatd,const Mesh & Th,const Mesh3::Element & K,const RdHat &PHat, RNMK_ & val) const;
 } ;

static TypeOfFE_P0Lagrange3d  P0_3d;
GTypeOfFE<Mesh3> & P0Lagrange3d(P0_3d);

static TypeOfFE_P1Lagrange3d  P1_3d;
GTypeOfFE<Mesh3> & P1Lagrange3d(P1_3d);

static TypeOfFE_P2Lagrange3d  P2_3d;
GTypeOfFE<Mesh3> & P2Lagrange3d(P2_3d);

static TypeOfFE_P1bLagrange3d  P1b_3d;
GTypeOfFE<Mesh3> & P1bLagrange3d(P1b_3d);

static TypeOfFE_RT0_3d  RT0_3d;
GTypeOfFE<Mesh3> & RT03d(RT0_3d);

static TypeOfFE_Edge0_3d  Edge0_3d;
GTypeOfFE<Mesh3> & Edge03d(Edge0_3d);

 // add 22 march 2021 FH ..
 static TypeOfFE_PkdcLagrange3d  P1dc_3d(1);
 static TypeOfFE_PkdcLagrange3d  P2dc_3d(2);
 static TypeOfFE_PkdcLagrange3d  P3dc_3d(3);
 static TypeOfFE_PkdcLagrange3d  P4dc_3d(4);

 
 
 static TypeOfFE_ConstDC<Mesh3>  P0Edge_3d(1,2);
 static TypeOfFE_ConstDC<Mesh3>  P0Edgedc_3d(1,1);
 static TypeOfFE_ConstDC<Mesh3>  P0Face_3d(2,2);
 static TypeOfFE_ConstDC<Mesh3>  P0Facedc_3d(2,1);
 static TypeOfFE_ConstDC<Mesh3>  P0VF_3d(3,2);
 static TypeOfFE_ConstDC<Mesh3>  P0VFdc_3d(3,1);
 
 GTypeOfFE<Mesh3> & G_P1dc_3d(P1dc_3d);
 GTypeOfFE<Mesh3> & G_P2dc_3d(P2dc_3d);
 GTypeOfFE<Mesh3> & G_P3dc_3d(P3dc_3d);
 GTypeOfFE<Mesh3> & G_P4dc_3d(P4dc_3d);
 GTypeOfFE<Mesh3> & G_P0Edge_3d (P0Edge_3d);
 GTypeOfFE<Mesh3> & G_P0Edgedc_3d (P0Edgedc_3d);
 GTypeOfFE<Mesh3> & G_P0Face_3d (P0Face_3d);
 GTypeOfFE<Mesh3> & G_P0Facedc_3d (P0Facedc_3d);
 GTypeOfFE<Mesh3> & G_P0VF_3d (P0VF_3d);
 GTypeOfFE<Mesh3> & G_P0VFdc_3d (P0VFdc_3d);

 
template<> GTypeOfFE<Mesh3> & DataFE<Mesh3>::P0=P0_3d;
template<> GTypeOfFE<Mesh3> & DataFE<Mesh3>::P1=P1_3d; 
template<> GTypeOfFE<Mesh3> & DataFE<Mesh3>::P2=P2_3d; 

 
     
     

 }

