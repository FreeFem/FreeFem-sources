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
         
     if (whatd & Fop_D0)
       {
	 f0[0] = l[0];
	 f0[1] = l[1];
	 f0[2] = l[2];
       }
         
     if (whatd & Fop_D1)
       {
	 R3 Dl[3];
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
             
	 if (whatd & Fop_dz) {
	   RN_ f0z(val('.',0,op_dz));
	   f0z[0] = Dl[0].z;
	   f0z[1] = Dl[1].z;
	   f0z[2] = Dl[2].z;
	 }
       }
         
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
         
     if (whatd & (Fop_D1|Fop_D2))
       {
             
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
         
   }
     
// RT0 surface

      
      
      
      class TypeOfFE_RT0_surf : public GTypeOfFE<MeshS>  {
      public:
          typedef MeshS Mesh;
          typedef MeshS::Element  Element;
          
          typedef GFElement<MeshS> FElement;
          static int dfon[];
          static const int d=Mesh::Rd::d;
          //static const int d=Mesh::RdHat::d;
          TypeOfFE_RT0_surf();
          void FB(const What_d whatd,const Mesh & Th,const Element & K,const RdHat &PHat, RNMK_ & val) const;
          void set(const Mesh & Th,const Element & K,InterpolationMatrix<RdHat> & M,int ocoef,int odf,int *nump) const;
      } ;
      
      
      
      int TypeOfFE_RT0_surf::dfon[]={0,1,0,0};   // dofs per vertice, edge, face, volume
      
      
      TypeOfFE_RT0_surf::TypeOfFE_RT0_surf(): GTypeOfFE<MeshS>(dfon,d,1,3*3,3,false,true)
      {
          //  integration on middle of faces  (light ) on  each face ..
          R2 Pt[]={ R2(0.5,0.5), R2(0.0,0.5), R2(0.5,0.0) };
         
          for (int i=0;i<3;++i)
              this->PtInterpolation[i]=Pt[i];
       
          {
              int i=0;
              for (int e=0;e<3;e++) //loop on edge
                      for (int c=0;c<3;c++,i++) {
                        this->pInterpolation[i]=e;
                        this->cInterpolation[i]=c;
                        this->dofInterpolation[i]=e;
                        this->coefInterpolation[i]=0.;
                      }
             }
      }
      
      
      
      
      void TypeOfFE_RT0_surf::set(const Mesh & Th,const Element & K,InterpolationMatrix<RdHat> & M ,int ocoef,int odf,int *nump) const
      {
         
         cout << "RT0 surface ongoing" << endl;
          ffassert(0);
          //   compute de coef d'interpolation
          int i=ocoef;
          R3 Nk= K.Edge(0)^K.Edge(1);
          for (int e=0;e<3;e++)
          {
              
              R3 N= (K.EdgeOrientation(e)*Nk)^(K.Edge(e)); //^(K.NormalS(e)); //   K.N(e);//  exterior and  ||N|| = 2* area f
              N /=N.norme();
              N *= K.EdgeOrientation(e);//  exterior and  ||N|| = 2* area f
             
              for (int c=0;c<3;c++,i++)
                M.coef[i]=N[c];
              
              }
          
      }
  
      
      void  TypeOfFE_RT0_surf::FB(const What_d whatd,const MeshS & Th,const MeshS::Element & K,const RdHat &PHat, RNMK_ & val) const
      {
          assert(val.N()>=Element::ne); //cout << "RT0 surface ongoing" << endl;
          ffassert(0);
          assert(val.M()==3 );
          // wi = signe * (x - qi)/ (volume*d)
          val=0;
          R3 P=K(PHat);
          R3 A[3]={(K[0]), (K[1]), (K[2])};
          R cc =1./(2.*K.mesure());
          R ci[3]={cc*K.EdgeOrientation(0),cc*K.EdgeOrientation(1),cc*K.EdgeOrientation(2)};     /// orientation????????
          
          if (whatd & Fop_D0)
          {
              
              for(int i=0;i<3;++i)
              {
                  val(i,0,op_id) = (P.x-A[i].x)*ci[i] ; //cout << "test " << (P.x-A[i].x) << " " << (P.y-A[i].y) << " " << (P.z-A[i].z) << " " <<ci[i] << endl;
                  val(i,1,op_id) = (P.y-A[i].y)*ci[i] ;
                  val(i,2,op_id) = (P.z-A[i].z)*ci[i] ;
                  //cout  << "RT 3d  "<<i << " "<< X << " " <<wi << " fo: " << K.faceOrient(i) <<endl;
              }
          }
          
          if (whatd & Fop_D1)
          {
              RN_ Ci(ci,3);
              if (whatd & Fop_dx)
                  val('.',0,op_dx) = Ci;  // cout << " deriv 1 " << Ci << endl;
              if (whatd & Fop_dy)
                  val('.',1,op_dy) = Ci;//cout << " deriv 2" << Ci << endl;
              if (whatd & Fop_dz)
                  val('.',2,op_dz) = Ci;//cout << " deriv 3" << Ci << endl;
          }
          
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
 	 
 	 assert(val.N() >=4);
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
 	     ffassert(0); // a faire ...
 	 //  cout << val << endl;
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
   
   	  
   template<> GTypeOfFE<MeshS> & DataFE<MeshS>::P0=P0_surf;
   template<> GTypeOfFE<MeshS> & DataFE<MeshS>::P1=P1_surf;
   template<> GTypeOfFE<MeshS> & DataFE<MeshS>::P2=P2_surf;
   template<> GTypeOfFE<MeshS> & DataFE<MeshS>::RT0=RT0_surf;
     
     
     
     
     
     
     
     
     
     
     
     
     
     
 }

