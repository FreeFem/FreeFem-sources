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
     
   // surface P0 
     
   class TypeOfFE_P0Lagrange_curve : public  TypeOfFE_Lagrange<MeshL> {
   public:
     TypeOfFE_P0Lagrange_curve(): TypeOfFE_Lagrange<MeshL>(0) {  }
     void FB(const What_d whatd,const Mesh & Th,const Element & K,const RdHat &PHat, RNMK_ & val) const;
     virtual R operator()(const FElement & K,const  RdHat & PHat,const KN_<R> & u,int componante,int op) const ;
   } ;
     
   R TypeOfFE_P0Lagrange_curve::operator()(const FElement & K,const  RdHat & PHat,const KN_<R> & u,int componante,int op) const
   {
     R u0(u(K(0)));
     R r=0;
     if (op==0)
       r = u0;
     else  r=0;
     return r;
   }
     
   void TypeOfFE_P0Lagrange_curve::FB(const What_d whatd,const Mesh & ,const Element & K,const RdHat & PHat,RNMK_ & val) const
   {
     assert(val.N() >=1);
     assert(val.M()==1 );
         
     val=0;
     RN_ f0(val('.',0,op_id));
         
     if (whatd & Fop_D0)
       f0[0] = 1;
         
   }
     
     
   // curve P1
   class TypeOfFE_P1Lagrange_curve : public  TypeOfFE_Lagrange<MeshL> {
   public:
     TypeOfFE_P1Lagrange_curve(): TypeOfFE_Lagrange<MeshL>(1) {  }
     void FB(const What_d whatd,const Mesh & Th,const Element & K,const RdHat &PHat, RNMK_ & val) const;
     virtual R operator()(const FElement & K,const  RdHat & PHat,const KN_<R> & u,int componante,int op) const ;
         
   } ;
     
     
     
   R TypeOfFE_P1Lagrange_curve::operator()(const FElement & K,const RdHat & PHat,const KN_<R> & u,int componante,int op) const
   {
     R u0(u(K(0))), u1(u(K(1)));
         
     R r=0;
     if (op==0) {
         R l0=1-PHat.x,l1=PHat.x;
         r = u0*l0+u1*l1;
     }
     else if(op==op_dx || op==op_dy || op==op_dz) {
         const Element & T=K.T;
         R3 D[2];
         T.Gradlambda(D);
         if (op==op_dx)
             r =  D[0].x*u0 + D[1].x*u1;
         else if (op==op_dy)
             r =  D[0].y*u0 + D[1].y*u1;
         else
             r =  D[0].z*u0 + D[1].z*u1;
     }
    return r;
   }
     
     
     
   void TypeOfFE_P1Lagrange_curve::FB(const What_d whatd,const Mesh & Th,const Element & K,const RdHat & PHat,RNMK_ & val) const
   {
         
     assert(val.N() >=Element::nv);
     assert(val.M()==1 );
       
     R l[]={1.-PHat.x,PHat.x};
         
     val=0;
     RN_ f0(val('.',0,op_id));
         
     if (whatd & Fop_D0) {
         f0[0] = l[0];
         f0[1] = l[1];
     }
         
     if (whatd & Fop_D1) {
         R3 Dl[2];
         K.Gradlambda(Dl);
             
         if (whatd & Fop_dx) {
             RN_ f0x(val('.',0,op_dx));
             f0x[0] = Dl[0].x;
             f0x[1] = Dl[1].x;
         }
             
         if (whatd & Fop_dy) {
             RN_ f0y(val('.',0,op_dy));
             f0y[0] = Dl[0].y;
             f0y[1] = Dl[1].y;
         }
             
         if (whatd & Fop_dz) {
             RN_ f0z(val('.',0,op_dz));
             f0z[0] = Dl[0].z;
             f0z[1] = Dl[1].z;
         }
     }
         
   }
     
     
     
   // curve P2
   class TypeOfFE_P2Lagrange_curve : public TypeOfFE_Lagrange<MeshL>  {
   public:
     typedef MeshL Mesh;
     typedef GFElement<MeshL> FElement;
     TypeOfFE_P2Lagrange_curve(): TypeOfFE_Lagrange<MeshL>(2) {  }
     void FB(const What_d whatd,const Mesh & Th,const MeshL::Element & K,const RdHat &PHat, RNMK_ & val) const;
   } ;
     
     
     
     
   void TypeOfFE_P2Lagrange_curve::FB(const What_d whatd,const Mesh & ,const Element & K,const RdHat & PHat,RNMK_ & val) const
   {
       
       
       R l[]={1.-PHat.x,PHat.x};
       
       assert(val.N() >=E::nv);
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
           assert(k==3);
       }
       
       if (whatd & (Fop_D1|Fop_D2))
       {
           
           R3 Dl[2];
           R l3[2]={ (4*l[0]-1),(4*l[1]-1)};
           
           K.Gradlambda(Dl);
           
           if( whatd & Fop_dx)
           {
               RN_ f0x(val('.',0,op_dx));
               
               int k=0;
               for(int i=0;i<E::nv;++i,++k)
               {
                   f0x[k] = Dl[i].x*l3[i];
               }
               for(int i=0;i<E::ne;++i,++k)
               {
                   int i0=E::nvedge[i][0],i1=E::nvedge[i][1];
                   f0x[k] = 4*(Dl[i1].x*l[i0] + Dl[i0].x*l[i1]) ;
               }
               assert(k==3);
           }
           
           if( whatd & Fop_dy)
           {
               RN_ f0y(val('.',0,op_dy));
               
               int k=0;
               for(int i=0;i<E::nv;++i,++k)
               {
                   f0y[k] = Dl[i].y*l3[i];
               }
               for(int i=0;i<E::ne;++i,++k)
               {
                   int i0=E::nvedge[i][0],i1=E::nvedge[i][1];
                   f0y[k] = 4*(Dl[i1].y*l[i0] + Dl[i0].y*l[i1]) ;
               }
               assert(k==3);
           }
           
           if( whatd & Fop_dz)
           {
               RN_ f0z(val('.',0,op_dz));
               
               int k=0;
               for(int i=0;i<E::nv;++i,++k)
               {
                   f0z[k] = Dl[i].z*l3[i];
               }
               for(int i=0;i<E::ne;++i,++k)
               {
                   int i0=E::nvedge[i][0],i1=E::nvedge[i][1];
                   f0z[k] = 4*(Dl[i1].z*l[i0] + Dl[i0].z*l[i1]) ;
               }
               assert(k==3);
           }
           if( whatd & Fop_D2)
               ffassert(0); // to DO FH. sorry
       }
       
   }
     
 
 
   static TypeOfFE_P0Lagrange_curve P0_curve;
   GTypeOfFE<MeshL> & P0Lagrange_curve(P0_curve);
     
   static TypeOfFE_P1Lagrange_curve P1_curve;
   GTypeOfFE<MeshL> & P1Lagrange_curve(P1_curve);
     
   static TypeOfFE_P2Lagrange_curve P2_curve;
   GTypeOfFE<MeshL> & P2Lagrange_curve(P2_curve);
      
   template<> GTypeOfFE<MeshL> & DataFE<MeshL>::P0=P0_curve;
   template<> GTypeOfFE<MeshL> & DataFE<MeshL>::P1=P1_curve;
   template<> GTypeOfFE<MeshL> & DataFE<MeshL>::P2=P2_curve;
     
  static TypeOfFE_LagrangeDC<MeshL>  P1dc_L(1,0.001);
  static TypeOfFE_LagrangeDC<MeshL>   P2dc_L(2,0.001);
  static TypeOfFE_LagrangeDC<MeshL>   P3dc_L(3,0.001);
  static TypeOfFE_LagrangeDC<MeshL>   P4dc_L(4,0.001);

  
  static TypeOfFE_ConstDC<MeshL>  P0VF_L(3,2);
  static TypeOfFE_ConstDC<MeshL>  P0VFdc_L(3,1);
  
  GTypeOfFE<MeshL> & G_P1dc_L(P1dc_L);
  GTypeOfFE<MeshL> & G_P2dc_L(P2dc_L);
  GTypeOfFE<MeshL> & G_P3dc_L(P3dc_L);
  GTypeOfFE<MeshL> & G_P4dc_L(P4dc_L);
  GTypeOfFE<MeshL> & G_P0VF_L (P0VF_L);
  GTypeOfFE<MeshL> & G_P0VFdc_L (P0VFdc_L);

     
     
     
     
     
     
     
     
     
     
     
     
 }

