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
     
   class TypeOfFE_P0Lagrange_surf : public  TypeOfFE_Lagrange<MeshS> {
   public:
     TypeOfFE_P0Lagrange_surf(): TypeOfFE_Lagrange<MeshS>(0) {  }
     void FB(const What_d whatd,const Mesh & Th,const Element & K,const RdHat &PHat, RNMK_ & val) const;
     virtual R operator()(const FElement & K,const  RdHat & PHat,const KN_<R> & u,int componante,int op) const ;
   } ;
     
   R TypeOfFE_P0Lagrange_surf::operator()(const FElement & K,const  RdHat & PHat,const KN_<R> & u,int componante,int op) const
   {
     R u0(u(K(0)));
     R r=0;
     if (op==0)
       r = u0;
     else  r=0;
     return r;
   }
     
   void TypeOfFE_P0Lagrange_surf::FB(const What_d whatd,const Mesh & ,const Element & K,const RdHat & PHat,RNMK_ & val) const
   {
     assert(val.N() >=1);
     assert(val.M()==1 );
         
     val=0;
     RN_ f0(val('.',0,op_id));
         
     if (whatd & Fop_D0)
       f0[0] = 1;
         
   }
     
     
   // surface P1 
   class TypeOfFE_P1Lagrange_surf : public  TypeOfFE_Lagrange<MeshS> {
   public:
     TypeOfFE_P1Lagrange_surf(): TypeOfFE_Lagrange<MeshS>(1) {  }
     void FB(const What_d whatd,const Mesh & Th,const Element & K,const RdHat &PHat, RNMK_ & val) const;
     virtual R operator()(const FElement & K,const  RdHat & PHat,const KN_<R> & u,int componante,int op) const ;
         
   } ;
     
     
     
   R TypeOfFE_P1Lagrange_surf::operator()(const FElement & K,const RdHat & PHat,const KN_<R> & u,int componante,int op) const
   {
     R u0(u(K(0))), u1(u(K(1))), u2(u(K(2)));
         
     R r=0;
     if (op==0)
       {
	 R l0=1-PHat.x-PHat.y,l1=PHat.x,l2=PHat.y;
	 r = u0*l0+u1*l1+l2*u2;
       }
     else if(op==op_dx || op==op_dy || op==op_dz) // dx => dy thank to Pichon 27/01/2008 (FH)
       {
	 const Element & T=K.T;
	 R3 D[3];
	 T.Gradlambda(D);
	 if (op==op_dx)
	   r =  D[0].x*u0 + D[1].x*u1 + D[2].x*u2;
	 else if (op==op_dy)
	   r =  D[0].y*u0 + D[1].y*u1 + D[2].y*u2;
	 else
	   r =  D[0].z*u0 + D[1].z*u1 + D[2].z*u2;
       }
     //  cout << r << "\t";
     return r;
   }
     
     
     
   void TypeOfFE_P1Lagrange_surf::FB(const What_d whatd,const Mesh & Th,const Element & K,const RdHat & PHat,RNMK_ & val) const
   {
         
     assert(val.N() >=Element::nv);
     assert(val.M()==1 );
     R l[]={1.-PHat.sum(),PHat.x,PHat.y};
         
     val=0;
     RN_ f0(val('.',0,op_id));
         
     if (whatd & Fop_D0) {
	   f0[0] = l[0];
	   f0[1] = l[1];
	   f0[2] = l[2];
     }
         
     if (whatd & Fop_D1) {
	   R3 Dl[3];
	   K.Gradlambda(Dl);
             
	 if (whatd & Fop_dx) {
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
	 if (whatd & Fop_dz) {
	   RN_ f0z(val('.',0,op_dz));
	   f0z[0] = Dl[0].z;
	   f0z[1] = Dl[1].z;
	   f0z[2] = Dl[2].z;
	 }
   }
   else if (whatd & Fop_D2)
     ffassert(0);
         
   }
     
     
     
   // surface P2
   class TypeOfFE_P2Lagrange_surf : public TypeOfFE_Lagrange<MeshS>  {
   public:
     typedef MeshS Mesh;
     typedef GFElement<MeshS> FElement;
     TypeOfFE_P2Lagrange_surf(): TypeOfFE_Lagrange<MeshS>(2) {  }
     void FB(const What_d whatd,const Mesh & Th,const MeshS::Element & K,const RdHat &PHat, RNMK_ & val) const;
   } ;
     
     
     
     
   void TypeOfFE_P2Lagrange_surf::FB(const What_d whatd,const Mesh & ,const Element & K,const RdHat & PHat,RNMK_ & val) const
   {
         
     R l[]={1.-PHat.sum(),PHat.x,PHat.y};
         
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
	 assert(k==6);
       }
         
     if (whatd & Fop_D1) {
             
	 R3 Dl[3];
	 R l4[3]={ (4*l[0]-1),(4*l[1]-1),(4*l[2]-1)};
             
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
	     assert(k==6);
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
	     assert(k==6);
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
	     assert(k==6);
	   }
     }
    else if (whatd & Fop_D2)
       ffassert(0);
   }
     
// RT0 surface

   class TypeOfFE_RT0_surf : public GTypeOfFE<MeshS>  {
     public:
       typedef MeshS Mesh;
          typedef MeshS::Element  Element;
          typedef GFElement<MeshS> FElement;
          static int dfon[];
          static const int d=Mesh::Rd::d;
 
          TypeOfFE_RT0_surf();
          void FB(const What_d whatd,const Mesh & Th,const Element & K,const RdHat &PHat, RNMK_ & val) const;
          void set(const Mesh & Th,const Element & K,InterpolationMatrix<RdHat> & M,int ocoef,int odf,int *nump) const;
    } ;
      
    int TypeOfFE_RT0_surf::dfon[]={0,1,0,0};   // dofs per vertice, edge, face, volume
      
      
    TypeOfFE_RT0_surf::TypeOfFE_RT0_surf(): GTypeOfFE<MeshS>(dfon,d,1,3*3,3,false,true) {
        
        R2 Pt[]={ R2(0.5,0.5), R2(0.0,0.5), R2(0.5,0.0) };
         
          for (int i=0;i<3;++i)
              this->PtInterpolation[i]=Pt[i];
         int i=0;
          for (int e=0;e<3;e++) //loop on edge
              for (int c=0;c<3;c++,i++) {
                  this->pInterpolation[i]=e;
                  this->cInterpolation[i]=c;
                  this->dofInterpolation[i]=e;
                  this->coefInterpolation[i]=0.;
              }
             
      }
      
      
      
      
   void TypeOfFE_RT0_surf::set(const Mesh & Th,const Element & K,InterpolationMatrix<RdHat> & M ,int ocoef,int odf,int *nump) const
   {
         
     int i=ocoef;
     for (int e=0;e<3;e++) {
       R signe = K.EdgeOrientation(e) ;
       R a=signe/(2.*K.mesure());
       R3 N=K.N(e)*a;
       for (int c=0;c<3;c++,i++)
         M.coef[i]=N[c];
     }
   }
  
      
   void  TypeOfFE_RT0_surf::FB(const What_d whatd,const Mesh & Th,const MeshS::Element & K,const RdHat &PHat, RNMK_ & val) const
   {
     assert(val.N()>=Element::ne);
     assert(val.M()==3 );
     val=0;

     R3 A(K[0]),B(K[1]),C(K[2]);
     R la=1-PHat.sum(),lb=PHat.x,lc=PHat.y;
     R3 D[3];
     K.Gradlambda(D);
     //loc basis
     R3 AB(A,B),AC(A,C),BA(B,A),BC(B,C),CA(C,A),CB(C,B);
     R c[3];
     c[0] = 1./(((AB,D[1]) + (AC,D[2])) *K.mesure()) ; c[0]*= K.EdgeOrientation(0) ;
     c[1] = 1./(((BA,D[0]) + (BC,D[2])) *K.mesure()) ; c[1]*= K.EdgeOrientation(1);
     c[2] = 1./(((CA,D[0]) + (CB,D[1])) *K.mesure()) ; c[2]*= K.EdgeOrientation(2);

     R3 f[3];
     f[0] = AB*(lb*c[0]) + AC*(lc*c[0]);
     f[1] = BA*(la*c[1]) + BC*(lc*c[1]);
     f[2] = CA*(la*c[2]) + CB*(lb*c[2]);
 
       
     if (whatd & Fop_D0)
       for(int i=0;i<3;i++){
         val(i,0,op_id) = f[i].x;
         val(i,1,op_id) = f[i].y;
         val(i,2,op_id) = f[i].z;
       }

     if (whatd & Fop_D1) {
       if (whatd & Fop_dx){
         R3 fx[3];
         fx[0] = AB*(D[1].x*c[0]) + AC*(D[2].x*c[0]);
         fx[1] = BA*(D[0].x*c[1]) + BC*(D[2].x*c[1]);
         fx[2] = CA*(D[0].x*c[2]) + CB*(D[1].x*c[2]);
         for(int i=0;i<3;i++){
           val(i,0,op_dx) = fx[i].x;
           val(i,1,op_dx) = fx[i].y;
           val(i,2,op_dx) = fx[i].z;
         }
       }
       if (whatd & Fop_dy) {
         R3 fy[3];
         fy[0] = AB*(D[1].y*c[0]) + AC*(D[2].y*c[0]);
         fy[1] = BA*(D[0].y*c[1]) + BC*(D[2].y*c[1]);
         fy[2] = CA*(D[0].y*c[2]) + CB*(D[1].y*c[2]);
         for(int i=0;i<3;i++){
           val(i,0,op_dy) = fy[i].x;
           val(i,1,op_dy) = fy[i].y;
           val(i,2,op_dy) = fy[i].z;
         }
       }
         if (whatd & Fop_dz) {
           R3 fz[3];
           fz[0] = AB*(D[1].z*c[0]) + AC*(D[2].z*c[0]);
           fz[1] = BA*(D[0].z*c[1]) + BC*(D[2].z*c[1]);
           fz[2] = CA*(D[0].z*c[2]) + CB*(D[1].z*c[2]);
           for(int i=0;i<3;i++){
             val(i,0,op_dz) = fz[i].x;
             val(i,1,op_dz) = fz[i].y;
             val(i,2,op_dz) = fz[i].z;
           }
         }
     }
     else if (whatd & Fop_D2)
       ffassert(0);
          
   }

  class TypeOfFE_P1bLagrange_surf : public TypeOfFE_Lagrange<MeshS>  {
	      public:  
	  	 typedef MeshS Mesh;
	  	 typedef GFElement<MeshS> FElement;
	  	 TypeOfFE_P1bLagrange_surf(): TypeOfFE_Lagrange<MeshS>(-1) {  }
	  	 void FB(const What_d whatd,const Mesh & Th,const MeshS::Element & K,const RdHat &PHat, RNMK_ & val) const;
	  } ;
     
 
 
 
      void TypeOfFE_P1bLagrange_surf::FB(const What_d whatd,const Mesh & ,const Element & K,const RdHat & PHat,RNMK_ & val) const
      {
 	 //  const Triangle & K(FE.T);
 
	 
	 
 	 R l[]={1.-PHat.sum(),PHat.x,PHat.y};
 	 R lb=l[0]*l[1]*l[2]*9.;
 	 
 	 assert(val.N() >=3);
 	 assert(val.M()==1 );
	 
 	 val=0; 
 	 RN_ f0(val('.',0,op_id)); 
	 
 	 if (whatd & Fop_D0) {
		 f0[0] = l[0]-lb;
		 f0[1] = l[1]-lb;
		 f0[2] = l[2]-lb;
		 f0[3] = 3.*lb;
 	 }
 	 if (whatd & Fop_D1) {
 	       R3 Dl[3];
 	       K.Gradlambda(Dl);
 	       R3 Dlb = (Dl[0]*l[1]*l[2] + Dl[1]*l[0]*l[2] + Dl[2]*l[0]*l[1])*9.;
	       
 	       //for(int i=0;i<4;++i)
 	       //      cout << Dl[i] << endl;
	      
 	       if (whatd & Fop_dx) 
 		 {
 		     RN_ f0x(val('.',0,op_dx)); 
 		     f0x[0] = Dl[0].x-Dlb.x;
 		     f0x[1] = Dl[1].x-Dlb.x;
 		     f0x[2] = Dl[2].x-Dlb.x;
 		     f0x[3] = Dlb.x*3.;
		     
 		 }
	       
 	       if (whatd & Fop_dy) {
 		   RN_ f0y(val('.',0,op_dy)); 
 		   f0y[0] = Dl[0].y-Dlb.y;
 		   f0y[1] = Dl[1].y-Dlb.y;
 		   f0y[2] = Dl[2].y-Dlb.y;
 		   f0y[3] = Dlb.y*3.;
 	       }
	       
 	       if (whatd & Fop_dz) {
 		   RN_ f0z(val('.',0,op_dz)); 
 		   f0z[0] = Dl[0].z-Dlb.z;
 		   f0z[1] = Dl[1].z-Dlb.z;
 		   f0z[2] = Dl[2].z-Dlb.z;
 		   f0z[3] = Dlb.z*3.;
		   
 	       }
 	   }
 	 else if (whatd & Fop_D2)
 	     ffassert(0);
      }

  
  class TypeOfFE_P2bLagrange_surf : public TypeOfFE_Lagrange<MeshS>  {
          public:
         typedef MeshS Mesh;
         typedef GFElement<MeshS> FElement;
         TypeOfFE_P2bLagrange_surf(): TypeOfFE_Lagrange<MeshS>(-2) {  }
         void FB(const What_d whatd,const Mesh & Th,const MeshS::Element & K,const RdHat &PHat, RNMK_ & val) const;
      } ;



  void TypeOfFE_P2bLagrange_surf::FB(const What_d whatd,const Mesh & Th,const Element & K,const RdHat & PHat,RNMK_ & val) const
  {
        
    assert(val.N() >=Element::nv+Element::ne+1);
    assert(val.M()==1 );
       
    val=0;
    
    R l[]={1.-PHat.sum(),PHat.x,PHat.y};
    R lb=l[0]*l[1]*l[2]*3.;
    R l4_0=(4*l[0]-1),l4_1=(4*l[1]-1),l4_2=(4*l[2]-1);

     

    if (whatd & Fop_D0) {
      R lb4=lb*4;
      val(0,0,op_id) = l[0]*(2.*l[0]-1.)+lb;
      val(1,0,op_id) = l[1]*(2.*l[1]-1.)+lb;
      val(2,0,op_id) = l[2]*(2.*l[2]-1.)+lb;
      val(3,0,op_id) = 4.*l[1]*l[2]-lb4;
      val(4,0,op_id) = 4.*l[0]*l[2]-lb4;
      val(5,0,op_id) = 4.*l[1]*l[0]-lb4;
      val(6,0,op_id) = 9.*lb;
    }
      
    if (whatd & Fop_D1) {
      R3 Dl[3];
      K.Gradlambda(Dl);
      R3 Dlb((Dl[0]*l[1]*l[2]+Dl[1]*l[0]*l[2]+Dl[2]*l[0]*l[1])*3.), Dlb4(Dlb*4.);
         
      if (whatd & Fop_dx) {
        val(0,0,op_dx) = Dl[0].x*l4_0 +Dlb.x;
        val(1,0,op_dx) = Dl[1].x*l4_1 +Dlb.x;
        val(2,0,op_dx) = Dl[2].x*l4_2 +Dlb.x;
        val(3,0,op_dx) = 4.*(Dl[1].x*l[2] + Dl[2].x*l[1]) -Dlb4.x;
        val(4,0,op_dx) = 4.*(Dl[2].x*l[0] + Dl[0].x*l[2]) -Dlb4.x;
        val(5,0,op_dx) = 4.*(Dl[0].x*l[1] + Dl[1].x*l[0]) -Dlb4.x;
        val(6,0,op_dx) = 9.*Dlb.x;
      }
      if (whatd & Fop_dy) {
        val(0,0,op_dy) = Dl[0].y*l4_0 +Dlb.y;
        val(1,0,op_dy) = Dl[1].y*l4_1 +Dlb.y;
        val(2,0,op_dy) = Dl[2].y*l4_2 +Dlb.y;
        val(3,0,op_dy) = 4.*(Dl[1].y*l[2] + Dl[2].y*l[1]) -Dlb4.y;
        val(4,0,op_dy) = 4.*(Dl[2].y*l[0] + Dl[0].y*l[2]) -Dlb4.y;
        val(5,0,op_dy) = 4.*(Dl[0].y*l[1] + Dl[1].y*l[0]) -Dlb4.y;
        val(6,0,op_dy) = 9.*Dlb.y;
      }
      if (whatd & Fop_dz) {
        val(0,0,op_dz) = Dl[0].z*l4_0 +Dlb.z;
        val(1,0,op_dz) = Dl[1].z*l4_1 +Dlb.z;
        val(2,0,op_dz) = Dl[2].z*l4_2 +Dlb.z;
        val(3,0,op_dz) = 4.*(Dl[1].z*l[2] + Dl[2].z*l[1]) -Dlb4.z;
        val(4,0,op_dz) = 4.*(Dl[2].z*l[0] + Dl[0].z*l[2]) -Dlb4.z;
        val(5,0,op_dz) = 4.*(Dl[0].z*l[1] + Dl[1].z*l[0]) -Dlb4.z;
        val(6,0,op_dz) = 9.*Dlb.z;
      }

    }
    else if (whatd & Fop_D2)
      ffassert(0);
        
  }
  
   static TypeOfFE_P0Lagrange_surf P0_surf;
   GTypeOfFE<MeshS> & P0Lagrange_surf(P0_surf);
   static TypeOfFE_P1Lagrange_surf P1_surf;
   GTypeOfFE<MeshS> & P1Lagrange_surf(P1_surf);
   static TypeOfFE_P2Lagrange_surf P2_surf;
   GTypeOfFE<MeshS> & P2Lagrange_surf(P2_surf);
   static TypeOfFE_RT0_surf  RT0_surf;
   GTypeOfFE<MeshS> & RT0surf(RT0_surf);
   static TypeOfFE_P1bLagrange_surf P1b_surf;
   GTypeOfFE<MeshS> & P1bLagrange_surf(P1b_surf);
   static TypeOfFE_P2bLagrange_surf P2b_surf;
   GTypeOfFE<MeshS> & P2bLagrange_surf(P2b_surf);
  
   template<> GTypeOfFE<MeshS> & DataFE<MeshS>::P0=P0_surf;
   template<> GTypeOfFE<MeshS> & DataFE<MeshS>::P1=P1_surf;
   template<> GTypeOfFE<MeshS> & DataFE<MeshS>::P2=P2_surf;

 }

