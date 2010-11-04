// The  P2BR finite element : the Bernadi Raugel Finite Element
// F. Hecht, decembre 2005
// -------------
// See Bernardi, C., Raugel, G.: Analysis of some finite elements for the Stokes problem. Math. Comp. 44, 71-79 (1985).
//  It is  a 2d coupled FE 
// the Polynomial space is $ P1^2$ + 3 normals bubbles edges function $(P_2)$
// the degre of freedom is 6 values at of the 2 componantes at the  3 vertices
// and the 3 flux on the 3 edges  
//   So 9 degrees of freedom and  N= 2. 

// -----------------------  related files: 
//  to check  and validate  :  testFE.edp 
//  to get a real example   :  NSP2BRP0.edp
// ------------------------------------------------------------

// -----------------------

#include <ff++.hpp>
#include "AddNewFE.h"

namespace  Fem2D {
  
  class TypeOfFE_P2BRLagrange : public  TypeOfFE { public:  
    static int Data[];
    // double Pi_h_coef[];
    
    TypeOfFE_P2BRLagrange(): TypeOfFE(6+3+0,
				      2,
				      Data,
				      4,
				      1,
				      6+3*(2+2), // nb coef to build interpolation
				      9, // np point to build interpolation
				      0)
    {  
      const double gauss1=(1.-sqrt(1./3.))/2;
      const double gauss2=1.-gauss1;
      // faux 
      const R2 Pt[] = { R2(0,0),R2(1,0),R2(0,1)}; 
      // for the 3 vertices 6 coef 
      int kk=0;
      for (int p=0;p<3;p++)
	{ 
	  P_Pi_h[p]=Pt[p];
	  pij_alpha[kk]= IPJ(kk,p,0);
	  ++kk;
	  pij_alpha[kk]= IPJ(kk,p,1); 
	  ++kk;
	}
      // for 
      int p=3;
      for (int e=0;e<3;++e)
	{ // point d'integration sur l'arete e 
	  R2 A=Pt[VerticesOfTriangularEdge[e][0]];
	  R2 B=Pt[VerticesOfTriangularEdge[e][1]];
	  P_Pi_h[p]= A*gauss1+B*gauss2;
	  //	  cout <<"\n" <<  p << " --  " << P_Pi_h[p] << " ::  " << A << " " << B << endl;
	  pij_alpha[kk++]= IPJ(6+e,p,0); // coef = 0.5* l_e *ne_x * sge	  
	  pij_alpha[kk++]= IPJ(6+e,p,1); // coef = 0.5* l_e *ne_y * sge	  	  
	  p++;
	  P_Pi_h[p]= A*gauss2+B*gauss1;
	  // cout << p << " ++ " << P_Pi_h[p] << endl;
	  pij_alpha[kk++]= IPJ(6+e,p,0); // coef = 0.5* l_e *ne_x * sge	  	  
	  pij_alpha[kk++]= IPJ(6+e,p,1); // coef = 0.5* l_e *ne_y * sge	  	  
	  p++;
	}
      assert(P_Pi_h.N()==p);
      assert(pij_alpha.N()==kk);

     }
    void FB(const bool * whatd, const Mesh & Th,const Triangle & K,const R2 &P, RNMK_ & val) const;
    void Pi_h_alpha(const baseFElement & K,KN_<double> & v) const;
  } ;
  //                     on what     nu df on node node of df    
 int TypeOfFE_P2BRLagrange::Data[]={
   0,0, 1,1, 2,2,  3,4,5,
   0,1, 0,1, 0,1,  0,0,0,
   0,0, 1,1, 2,2,  3,4,5,
   0,0, 0,0, 0,0,  0,0,0,
   0,1, 2,3, 4,5,  6,7,8, 
   0,0, 
   0,0,
   0,9 
};

void TypeOfFE_P2BRLagrange::Pi_h_alpha(const baseFElement & K,KN_<double> & v) const
  {
    const Triangle & T(K.T);
    int k=0;
    // coef pour les 3 sommets  fois le 2 composantes 
    for (int i=0;i<6;i++)
      v[k++]=1; 
    //   integration sur les aretes 
    for (int i=0;i<3;i++)
      {
	
        R2 N(T.Edge(i).perp());
	N  *= T.EdgeOrientation(i)*0.5 ;
        v[k++]= N.x; 
        v[k++]= N.y;
        v[k++]= N.x;
        v[k++]= N.y;
      }
  }
  
  void TypeOfFE_P2BRLagrange::FB(const bool * whatd,const Mesh & ,const Triangle & K,const R2 & P,RNMK_ & val) const
  {
    R2 A(K[0]), B(K[1]),C(K[2]);
    R l0=1-P.x-P.y,l1=P.x,l2=P.y; 
    R l4_0=(4*l0-1),l4_1=(4*l1-1),l4_2=(4*l2-1); 
    //  $int_e_1 l0*l0 = |e_1| /3 $ et  $int_e_1 l0*l1 = |e_1| /6 $
    //  pour avoir  flux = 1   
    //  
    R2 E[3]={ K.Edge(0),K.Edge(1),K.Edge(2)};
    double l2E[3]={  (E[0],E[0]),  (E[1],E[1]),  (E[2],E[2]) };  
    double lE[3]={  sqrt(l2E[0]), sqrt(l2E[1]), sqrt(l2E[2]) };
    double sgE[3]={ K.EdgeOrientation(0), K.EdgeOrientation(1), K.EdgeOrientation(2)}; 
    R2 cN[3]= { 
      E[0].perp() *(6.*sgE[0]/l2E[0]), 
      E[1].perp() *(6.*sgE[1]/l2E[1]),
      E[2].perp() *(6.*sgE[2]/l2E[2]) 
    };

    val=0; 
    
    throwassert( val.N()>=9);
    throwassert(val.M()==2);
    
    
    val=0; 


  if (whatd[op_id])
    {
      RN_ f0(val('.',0,op_id)); 
      RN_ f1(val('.',1,op_id)); 
      
      f1[1]=f0[0] = l0;
      f1[3]=f0[2] = l1;
      f1[5]=f0[4] = l2;
      
      f0[6] = cN[0].x*l1*l2; // oppose au sommet 0
      f0[7] = cN[1].x*l0*l2; // oppose au sommet 1
      f0[8] = cN[2].x*l1*l0; // oppose au sommet 3
      
      f1[6] = cN[0].y*l1*l2; // oppose au sommet 0
      f1[7] = cN[1].y*l0*l2; // oppose au sommet 1
      f1[8] = cN[2].y*l1*l0; // oppose au sommet 3
      
    }
  
  if(  whatd[op_dx] || whatd[op_dy] || whatd[op_dxx] || whatd[op_dyy] ||  whatd[op_dxy])
    {
      R2 Dl0(K.H(0)), Dl1(K.H(1)), Dl2(K.H(2));
      if (whatd[op_dx])
	{
	  RN_ f0x(val('.',0,op_dx)); 
	  RN_ f1x(val('.',1,op_dx)); 
	  
	  f1x[1]=f0x[0] = Dl0.x;
	  f1x[3]=f0x[2] = Dl1.x;
	  f1x[5]=f0x[4] = Dl2.x;
	  
	  
	  f0x[6] = cN[0].x*(Dl1.x*l2 + Dl2.x*l1) ;
	  f0x[7] = cN[1].x*(Dl2.x*l0 + Dl0.x*l2) ;
	  f0x[8] = cN[2].x*(Dl0.x*l1 + Dl1.x*l0) ;
	  
	  f1x[6] = cN[0].y*(Dl1.x*l2 + Dl2.x*l1) ;
	  f1x[7] = cN[1].y*(Dl2.x*l0 + Dl0.x*l2) ;
	  f1x[8] = cN[2].y*(Dl0.x*l1 + Dl1.x*l0) ;
	  
	}
      
      if (whatd[op_dy])
	{  
	  RN_ f0y(val('.',0,op_dy)); 
	  RN_ f1y(val('.',1,op_dy)); 
	  
	  f1y[1]=f0y[0] = Dl0.y;
	  f1y[3]=f0y[2] = Dl1.y;
	  f1y[5]=f0y[4] = Dl2.y;
	  
	  f0y[6] = cN[0].x*(Dl1.y*l2 + Dl2.y*l1) ;
	  f0y[7] = cN[1].x*(Dl2.y*l0 + Dl0.y*l2) ;
	  f0y[8] = cN[2].x*(Dl0.y*l1 + Dl1.y*l0) ;
	  
	  f1y[6] = cN[0].y*(Dl1.y*l2 + Dl2.y*l1) ;
	  f1y[7] = cN[1].y*(Dl2.y*l0 + Dl0.y*l2) ;
	  f1y[8] = cN[2].y*(Dl0.y*l1 + Dl1.y*l0) ;
	}
      
      if (whatd[op_dxx])
	{  
	  
	  
	  RN_ f0xx(val('.',0,op_dxx)); 
	  RN_ f1xx(val('.',1,op_dxx)); 
	  
	  
	  f0xx[6] =  2*cN[0].x*Dl1.x*Dl2.x;
	  f0xx[7] =  2*cN[1].x*Dl0.x*Dl2.x;
	  f0xx[8] =  2*cN[2].x*Dl0.x*Dl1.x;
	  f1xx[6] =  2*cN[0].y*Dl1.x*Dl2.x;
	  f1xx[7] =  2*cN[1].y*Dl0.x*Dl2.x;
	  f1xx[8] =  2*cN[2].y*Dl0.x*Dl1.x;
	}
      
      if (whatd[op_dyy])
	{  
	  RN_ f0yy(val('.',0,op_dyy)); 
	  RN_ f1yy(val('.',1,op_dyy)); 
	  
	  f0yy[6] =  2*cN[0].x*Dl1.y*Dl2.y;
	  f0yy[7] =  2*cN[1].x*Dl0.y*Dl2.y;
	  f0yy[8] =  2*cN[2].x*Dl0.y*Dl1.y;
	  f1yy[6] =  2*cN[0].y*Dl1.y*Dl2.y;
	  f1yy[7] =  2*cN[1].y*Dl0.y*Dl2.y;
	  f1yy[8] =  2*cN[2].y*Dl0.y*Dl1.y;
	}
      if (whatd[op_dxy])
	{  
	  assert(val.K()>op_dxy);
	  RN_ f0xy(val('.',0,op_dxy)); 
	  RN_ f1xy(val('.',1,op_dxy)); 
	  
	  
	  f0xy[6] =  cN[0].x*(Dl1.x*Dl2.y + Dl1.y*Dl2.x);
	  f0xy[7] =  cN[1].x*(Dl0.x*Dl2.y + Dl0.y*Dl2.x);
	  f0xy[8] =  cN[2].x*(Dl0.x*Dl1.y + Dl0.y*Dl1.x);
	  f1xy[6] =  cN[0].y*(Dl1.x*Dl2.y + Dl1.y*Dl2.x);
	  f1xy[7] =  cN[1].y*(Dl0.x*Dl2.y + Dl0.y*Dl2.x);
	  f1xy[8] =  cN[2].y*(Dl0.x*Dl1.y + Dl0.y*Dl1.x);
	}
      
    }
  //  now remove the flux part on 6 first DL 
  //  w_i = w_i - a_i w_{k_i} - b_i w_{l_i} ;
  {
    int    k[6]={ 6+1 , 6+1,  6+2 , 6+2, 6+0, 6+0 };
    int    l[6]={ 6+2 , 6+2,  6+0 , 6+0, 6+1, 6+1 };

    R2 eN[3]= { 
      E[0].perp() *(0.5*sgE[0]),
      E[1].perp() *(0.5*sgE[1]),
      E[2].perp() *(0.5*sgE[2])
    };

    double a[6]={ eN[1].x , eN[1].y,   eN[2].x , eN[2].y, eN[0].x , eN[0].y};
    double b[6]={ eN[2].x , eN[2].y,   eN[0].x , eN[0].y, eN[1].x , eN[1].y};
    int nop=0;

    int vop[last_operatortype];
    for (int j=0;j<last_operatortype;j++)
      if (whatd[j])
	vop[nop++] = j;

    for(int i=0;i<6;++i)
      for(int jj=0;jj<nop;++jj)
      {
        int j=vop[jj];
	val(i,0,j) -= a[i]*val(k[i],0,j)  +    b[i]*val(l[i],0,j)  ;
	val(i,1,j) -= a[i]*val(k[i],1,j)  +    b[i]*val(l[i],1,j)  ;
      }
  }
  
}
//  ----   cooking to add the finite elemet to freefem table --------  
// a static variable to def the finite element 
  static TypeOfFE_P2BRLagrange P2LagrangeP2BR;
  //  now adding   FE in FreeFEm++  table
  static AddNewFE P2BR("P2BR",&P2LagrangeP2BR); 
// --- end cooking   
} // end FEM2d namespace 
  
