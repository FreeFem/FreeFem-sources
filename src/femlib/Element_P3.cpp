#include "error.hpp"
#include "rgraph.hpp"
using namespace std;  
#include "RNM.hpp"
#include "fem.hpp"
#include "FESpace.hpp"

namespace  Fem2D {

 // ------ P3  Hierarchical (just remove P1 node of the P2 finite element)  --------
 class TypeOfFE_P3Lagrange : public  TypeOfFE { public:  
   const int k=3; 
   const ndf = k*(k+1)/2;
   static int Data[];
   static double Pi_h_coef[];
   static int nn[ndf][k];    
   static int aa[ndf][k];    
   static int ff[ndf];    
   static int il[ndf];    
   static int jl[ndf];    
   static int kl[ndf];    

   TypeOfFE_P3Lagrange(): TypeOfFE(3+2*3+1,1,Data,4,1,3,3,Pi_h_coef)
   {  

     for (int i=0;i<NbDoF;i++) {
       pij_alpha[i]= IPJ(i,i,0);
       P_Pi_h[i]=Pt[i]; }
   }
   void FB(const bool * whatd, const Mesh & Th,const Triangle & K,const R2 &P, RNMK_ & val) const;
 } ;
  //                     on what     nu df on node node of df    
  int TypeOfFE_P3Lagrange::Data[]={
    0,1,2,3,3,4,4,5,5,6,     // the support number  of the node of the df 
    0,0,0,0,0,0,0,0,0,0,     // the number of the df on  the node  
    0,1,2,4,5,6,7,8,9,10,    // the node of the df 
    0,0,0,0,0,0,0,0,0,0,     //  the df come from which FE (generaly 0)
    0,1,2,4,5,6,7,8,9,10,    //  which are de df on sub FE
    0,0                      // for each compontant $j=0,N-1$ it give the sub FE associated 
    
    0,1,2,       0};
  double TypeOfFE_P3Lagrange::Pi_h_coef[]={1.,1.,1.,1.,1.,1.,1.,1.,1.,1.};
  
  void TypeOfFE_P3Lagrange::FB(const bool * whatd,const Mesh & ,const Triangle & K,const R2 & P,RNMK_ & val) const
  {
 
    R2 A(K[0]), B(K[1]),C(K[2]);
    R l0=1-P.x-P.y,l1=P.x,l2=P.y; 
    R L[3]={l0*k,l1*k,l2*k};
    throwassert( val.N()>=10);
    throwassert(val.M()==1);
    
    val=0; 
    /*
    //  le fonction de base du Pk Lagrange sont 
    //  
    //
    */
// -- 
    
    
    if (whatd[op_id])
      {
	RN_ f0(val('.',0,op_id)); 
	for (int df=0;df<ndf;df++)
	  f0[df] = (L[nn[df][0]]-aa[df][0])* (L[nn[df][1]]-aa[df][1])* (L[nn[df][2]]-aa[df][2])/ff[df];
      }
    
    
    if(  whatd[op_dx] || whatd[op_dy] || whatd[op_dxx] || whatd[op_dyy] ||  whatd[op_dxy])
      {
	R2 D[]={K.H(0)*k, K.H(1)*k,K.H(2)*k };
	if (whatd[op_dx])
	  {
	    RN_ f0x(val('.',0,op_dx)); 
	    for (int df=0;df<ndf;df++)
	      f0x[df] = ( + ( D[nn[df][0]].x        )* (L[nn[df][1]]-aa[df][1])*(L[nn[df][2]]-aa[df][2])
			 + ( L[nn[df][0]]-aa[df][0])* (D[nn[df][1]].x        )*(L[nn[df][2]]-aa[df][2])
			 + ( L[nn[df][0]]-aa[df][0])* (L[nn[df][1]]-aa[df][1])*(D[nn[df][2]].x        )
			 )/ff[df];
	  }

	if (whatd[op_dy])
	  {  
	    RN_ f0y(val('.',0,op_dy)); 
	    for (int df=0;df<ndf;df++)
	      f0y[df] = ( + ( D[nn[df][0]].y        )* (L[nn[df][1]]-aa[df][1])*(L[nn[df][2]]-aa[df][2])
			 + ( L[nn[df][0]]-aa[df][0])* (D[nn[df][1]].y        )*(L[nn[df][2]]-aa[df][2])
			 + ( L[nn[df][0]]-aa[df][0])* (L[nn[df][1]]-aa[df][1])*(D[nn[df][2]].y        )
			 )/ff[df];
	    
	  }
	
	if (whatd[op_dxx])
	  {  
	    RN_ fxx(val('.',0,op_dxx)); 
	    for (int df=0;df<ndf;df++)
	      f0y[df] = ( + ( D[nn[df][0]].y        )* (L[nn[df][1]]-aa[df][1])*(L[nn[df][2]]-aa[df][2])
			  + ( L[nn[df][0]]-aa[df][0])* (D[nn[df][1]].y        )*(L[nn[df][2]]-aa[df][2])
			  + ( L[nn[df][0]]-aa[df][0])* (L[nn[df][1]]-aa[df][1])*(D[nn[df][2]].y        )
			  + ( D[nn[df][0]].y        )* (L[nn[df][1]]-aa[df][1])*(L[nn[df][2]]-aa[df][2])
			  + ( L[nn[df][0]]-aa[df][0])* (D[nn[df][1]].y        )*(L[nn[df][2]]-aa[df][2])
			  + ( L[nn[df][0]]-aa[df][0])* (L[nn[df][1]]-aa[df][1])*(D[nn[df][2]].y        )
			 )/ff[df];

	  }
	
	if (whatd[op_dyy])
	  {  
	    RN_ fyy(val('.',0,op_dyy)); 
	    assert(0);
	  }
	if (whatd[op_dxy])
	  {  
	    assert(val.K()>op_dxy);
	    RN_ fxy(val('.',0,op_dxy)); 
	    assert(0);
	  }
	
      }
  }
#include "Element_P3.hpp"
  


// link with FreeFem++ 
extern  ListOfTFE typefem_P3;
//#pragma export list typefem_P3
//#pragma lib_export list typefem_P3

static TypeOfFE_P3Lagrange P3LagrangeP3;
// given the name of the finite element in FreeFem++
 ListOfTFE typefem_P3("P3", &P3LagrangeP3);


// --- fin -- 
} // Fem2D 

