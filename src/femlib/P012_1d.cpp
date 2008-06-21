// ********** DO NOT REMOVE THIS BANNER **********
// ORIG-DATE:     Jan 2008
// -*- Mode : c++ -*-
//
// SUMMARY  : P0,P1,P2 lagrange 1D 
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


class TypeOfFE_P0Lagrange1d : public  TypeOfFE_Lagrange<Mesh1> { 
public:  
  TypeOfFE_P0Lagrange1d(): TypeOfFE_Lagrange<Mesh1>(0) {  }
  void FB(const What_d whatd,const Mesh & Th,const Element & K,const Rd &P, RNMK_ & val) const;
  virtual R operator()(const FElement & K,const  RdHat & PHat,const KN_<R> & u,int componante,int op) const ;
  
} ;



R TypeOfFE_P0Lagrange1d::operator()(const FElement & K,const  R1 & PHat,const KN_<R> & u,int componante,int op) const 
{ 
  R u0(u(K(0)));
  R r=0;
  if (op==0)
    r = u0;
  else  r=0;
  return r;
}

void TypeOfFE_P0Lagrange1d::FB(const What_d whatd,const Mesh & ,const Element & K,const R1 & P,RNMK_ & val) const
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
class TypeOfFE_P1Lagrange1d : public  TypeOfFE_Lagrange<Mesh1> { 
public:  
  TypeOfFE_P1Lagrange1d(): TypeOfFE_Lagrange<Mesh1>(1) {  }
  void FB(const What_d whatd,const Mesh & Th,const Element & K,const Rd &P, RNMK_ & val) const;
  virtual R operator()(const FElement & K,const  RdHat & PHat,const KN_<R> & u,int componante,int op) const ;
  
} ;



R TypeOfFE_P1Lagrange1d::operator()(const FElement & K,const  R1 & PHat,const KN_<R> & u,int componante,int op) const 
{ 
  R u0(u(K(0))), u1(u(K(1)));
  R r=0;
  if (op==0)
    {
      R l0=1-PHat.x,l1=PHat.x; 
      r = u0*l0+u1*l1;
    }
  else if(op==op_dx  )
    { 
      const Element & T=K.T;
      R1 D[2];
      T.Gradlambda(D);
	r =  D[0].x*u0 + D[1].x*u1 ;
    }
  //  cout << r << "\t";
  return r;
}

void TypeOfFE_P1Lagrange1d::FB(const What_d whatd,const Mesh & ,const Element & K,const R1 & P,RNMK_ & val) const
{
  //  const Triangle & K(FE.T);
  
  R l[]={1.-P.sum(),P.x}; 
  
  assert(val.N() >=Element::nv);
  assert(val.M()==1 );
  
  val=0; 
  RN_ f0(val('.',0,op_id)); 
  
  if (whatd & Fop_D0) 
    {
      f0[0] = l[0];
      f0[1] = l[1];
    }
  if (whatd & Fop_D1)
    {
      R1 Dl[3];
      K.Gradlambda(Dl);
      if (whatd & Fop_dx) 
	{
	  RN_ f0x(val('.',0,op_dx)); 
	  f0x[0] = Dl[0].x;
	  f0x[1] = Dl[1].x;
	  
	}
      

    }
}




class TypeOfFE_P2Lagrange1d : public  TypeOfFE_Lagrange<Mesh1>  { 
public:  
  typedef Mesh1 Mesh;
  typedef GFElement<Mesh1> FElement;
  TypeOfFE_P2Lagrange1d(): TypeOfFE_Lagrange<Mesh1>(2) {  }
  void FB(const What_d whatd,const Mesh & Th,const Mesh1::Element & K,const Rd &P, RNMK_ & val) const;
  
} ;



void TypeOfFE_P2Lagrange1d::FB(const What_d whatd,const Mesh & ,const Element & K,const R1 & P,RNMK_ & val) const
{
  //  const Triangle & K(FE.T);
  
  R l[]={1.-P.sum(),P.x}; 
  
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
      R1 Dl[2];
      R l4[2]={ (4*l[0]-1),(4*l[1]-1)}; 
      
      K.Gradlambda(Dl);
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
      assert(k==1);
      
      //cout << " D2 " << whatd <<  endl;
      if (whatd & Fop_D2)
      {
	  //cout << " D2 " << endl;
	  RN_ f0xx(val('.',0,op_dxx));
	  
	   k=0;
	  for(int i=0;i<E::nv;++i,++k)
	  {
	      f0xx[k] = 4.*Dl[i].x*Dl[i].x;
	  }
	  for(int i=0;i<E::ne;++i,++k)
	  {
	      int i0=E::nvedge[i][0],i1=E::nvedge[i][1];
	      f0xx[k] = 8.*Dl[i0].x*Dl[i1].x;
	  } 
      }
      
  } 
    
}



static TypeOfFE_P0Lagrange1d  P0_1d;
GTypeOfFE<Mesh1> & P0Lagrange1d(P0_1d);

static TypeOfFE_P1Lagrange1d  P1_1d;
GTypeOfFE<Mesh1> & P1Lagrange1d(P1_1d);

static TypeOfFE_P2Lagrange1d  P2_1d;
GTypeOfFE<Mesh1> & P2Lagrange1d(P2_1d);




template<> GTypeOfFE<Mesh1> & DataFE<Mesh1>::P0=P0_1d; 
template<> GTypeOfFE<Mesh1> & DataFE<Mesh1>::P1=P1_1d; 
template<> GTypeOfFE<Mesh1> & DataFE<Mesh1>::P2=P2_1d; 


 }

