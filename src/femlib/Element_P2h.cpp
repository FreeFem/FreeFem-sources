#include "error.hpp"
#include "rgraph.hpp"
using namespace std;  
#include "RNM.hpp"
#include "fem.hpp"
#include "FESpace.hpp"

namespace  Fem2D {

 // ------ P2h  Hierarchical (just remove P1 node of the P2 finite element)  --------
 class TypeOfFE_P2hLagrange : public  TypeOfFE { public:  
  static int Data[];
    static double Pi_h_coef[];

   TypeOfFE_P2hLagrange(): TypeOfFE(0,1,0,1,Data,4,1,3,3,Pi_h_coef)
    {  
    
       const R2 Pt[] = { R2(0.5,0.5), R2(0.0,0.5), R2(0.5,0.0) };
      for (int i=0;i<NbDoF;i++) {
       pij_alpha[i]= IPJ(i,i,0);
       P_Pi_h[i]=Pt[i]; }
     }
   void FB(const bool * whatd, const Mesh & Th,const Triangle & K,const R2 &P, RNMK_ & val) const;
} ;
//                     on what     nu df on node node of df    
 int TypeOfFE_P2hLagrange::Data[]={3,4,5,       0,0,0,       0,1,2,       0,0,0,        0,1,2,       0};
  double TypeOfFE_P2hLagrange::Pi_h_coef[]={1.,1.,1.};

 void TypeOfFE_P2hLagrange::FB(const bool * whatd,const Mesh & ,const Triangle & K,const R2 & P,RNMK_ & val) const
{
  R2 A(K[0]), B(K[1]),C(K[2]);
  R l0=1-P.x-P.y,l1=P.x,l2=P.y; 
  R l4_0=(4*l0-1),l4_1=(4*l1-1),l4_2=(4*l2-1); 
  
  throwassert( val.N()>=3);
  throwassert(val.M()==1);
  
  val=0; 
// --     
 if (whatd[op_id])
  {
   RN_ f0(val('.',0,op_id)); 
  f0[0] = 4*l1*l2; // oppose au sommet 0
  f0[1] = 4*l0*l2; // oppose au sommet 1
  f0[2] = 4*l1*l0; // oppose au sommet 3
  }
 if(  whatd[op_dx] || whatd[op_dy] || whatd[op_dxx] || whatd[op_dyy] ||  whatd[op_dxy])
 {
   R2 Dl0(K.H(0)), Dl1(K.H(1)), Dl2(K.H(2));
  if (whatd[op_dx])
  {
    RN_ f0x(val('.',0,op_dx)); 
  f0x[0] = 4*(Dl1.x*l2 + Dl2.x*l1) ;
  f0x[1] = 4*(Dl2.x*l0 + Dl0.x*l2) ;
  f0x[2] = 4*(Dl0.x*l1 + Dl1.x*l0) ;
  }

 if (whatd[op_dy])
  {  
    RN_ f0y(val('.',0,op_dy)); 
  f0y[0] = 4*(Dl1.y*l2 + Dl2.y*l1) ;
  f0y[1] = 4*(Dl2.y*l0 + Dl0.y*l2) ;
  f0y[2] = 4*(Dl0.y*l1 + Dl1.y*l0) ;
  }
 
 if (whatd[op_dxx])
  {  
    RN_ fxx(val('.',0,op_dxx)); 

    fxx[0] =  8*Dl1.x*Dl2.x;
    fxx[1] =  8*Dl0.x*Dl2.x;
    fxx[2] =  8*Dl0.x*Dl1.x;
  }

 if (whatd[op_dyy])
  {  
    RN_ fyy(val('.',0,op_dyy)); 
    fyy[0] =  8*Dl1.y*Dl2.y;
    fyy[1] =  8*Dl0.y*Dl2.y;
    fyy[2] =  8*Dl0.y*Dl1.y;
  }
 if (whatd[op_dxy])
  {  
    assert(val.K()>op_dxy);
    RN_ fxy(val('.',0,op_dxy)); 
  
    fxy[0] =  4*(Dl1.x*Dl2.y + Dl1.y*Dl2.x);
    fxy[1] =  4*(Dl0.x*Dl2.y + Dl0.y*Dl2.x);
    fxy[2] =  4*(Dl0.x*Dl1.y + Dl0.y*Dl1.x);
  }
 
 }
 
}
// link with FreeFem++ 
extern  ListOfTFE typefem_P2h;

static TypeOfFE_P2hLagrange P2LagrangeP2h;
// given the name of the finite element in FreeFem++
 ListOfTFE typefem_P2h("P2h", &P2LagrangeP2h);


// --- fin -- 
} // FEM2d namespace 
