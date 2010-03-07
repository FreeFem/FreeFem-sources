#include "ff++.hpp"
#include "AddNewFE.h"

// Attention probleme de numerotation des inconnues
// -------------------------------------------------
// dans freefem, il y a un noeud par objets  sommet, arete, element.
// et donc la numerotation des dl dans l'element depend 
// de l'orientation des aretes
// 
/// ---------------------------------------------------------------
namespace  Fem2D {

    
    // -------------------
    // ttdc1_ finite element fully discontinue. 
    // -------------------
    class TypeOfFE_P1ttdc1_ : public  TypeOfFE { public:  
	static int Data[];
	static double Pi_h_coef[];
	static const R2 G;
	static const R cshrink;
	static const R cshrink1;
	//  (1 -1/3)*
	
	static R2 Shrink(const R2& P){ return (P-G)*cshrink+G;}
	static R2 Shrink1(const R2& P){ return (P-G)*cshrink1+G;}
	
	TypeOfFE_P1ttdc1_(): TypeOfFE(0,0,3,1,Data,1,1,3,3,Pi_h_coef)
	{ const R2 Pt[] = { Shrink(R2(0,0)), Shrink(R2(1,0)), Shrink(R2(0,1)) }; 
	    for (int i=0;i<NbDoF;i++) {
		pij_alpha[i]= IPJ(i,i,0);
		P_Pi_h[i]=Pt[i];
		// cout << Pt[i] << " " ;
	    }
	    //	cout <<" cshrink: " << cshrink << " cshrink1 : "<< cshrink1 <<endl;
	}
	
	void FB(const bool * whatd,const Mesh & Th,const Triangle & K,const R2 &P, RNMK_ & val) const;
	
	
	virtual R operator()(const FElement & K,const  R2 & PHat,const KN_<R> & u,int componante,int op) const ;
	
    } ;
    const R2 TypeOfFE_P1ttdc1_::G(1./3.,1./3.);   
    const R TypeOfFE_P1ttdc1_::cshrink=1;   
    const R TypeOfFE_P1ttdc1_::cshrink1=1./TypeOfFE_P1ttdc1_::cshrink;   
    
    class TypeOfFE_P2ttdc1_ : public  TypeOfFE { public:  
	static int Data[];
	static double Pi_h_coef[];
	static const R2 G;
	static const R cshrink;
	static const R cshrink1;
	
	static R2 Shrink(const R2& P){ return (P-G)*cshrink+G;}
	static R2 Shrink1(const R2& P){ return (P-G)*cshrink1+G;}
	
	TypeOfFE_P2ttdc1_(): TypeOfFE(0,0,6,1,Data,3,1,6,6,Pi_h_coef)
	{ const R2 Pt[] = { Shrink(R2(0,0)), Shrink(R2(1,0)), Shrink(R2(0,1)),
	    Shrink(R2(0.5,0.5)),Shrink(R2(0,0.5)),Shrink(R2(0.5,0)) };
	    for (int i=0;i<NbDoF;i++) {
		pij_alpha[i]= IPJ(i,i,0);
		P_Pi_h[i]=Pt[i]; }
	}
	
	
	void FB(const bool * whatd,const Mesh & Th,const Triangle & K,const R2 &P, RNMK_ & val) const;
	
	
    } ;
    //                          on what   nu df on node  node of df    
    int TypeOfFE_P1ttdc1_::Data[]={6,6,6,       0,1,2,       0,0,0,       0,0,0,         0,1,2,       0, 0,3};
    int TypeOfFE_P2ttdc1_::Data[]={6,6,6,6,6,6, 0,1,2,3,4,5, 0,0,0,0,0,0,  0,0,0,0,0,0,  0,1,2,3,4,5, 0, 0,6};
    double TypeOfFE_P1ttdc1_::Pi_h_coef[]={1.,1.,1.};
    double TypeOfFE_P2ttdc1_::Pi_h_coef[]={1.,1.,1.,1.,1.,1.};
    
    const R2 TypeOfFE_P2ttdc1_::G(1./3.,1./3.);   
    const R TypeOfFE_P2ttdc1_::cshrink=1;   
    const R TypeOfFE_P2ttdc1_::cshrink1=1./TypeOfFE_P2ttdc1_::cshrink;   
    
    
    R TypeOfFE_P1ttdc1_::operator()(const FElement & K,const  R2 & P1Hat,const KN_<R> & u,int componante,int op) const 
    { 
	
	R2 PHat=Shrink1(P1Hat);  
	R u0(u(K(0))), u1(u(K(1))), u2(u(K(2)));
	R r=0;
	if (op==0)
	  {
	    R l0=1-PHat.x-PHat.y,l1=PHat.x,l2=PHat.y; 
	    r = u0*l0+u1*l1+l2*u2;
	  }
	else
	  { 
	      const Triangle & T=K.T;
	      R2 D0 = T.H(0)*cshrink1 , D1 = T.H(1)*cshrink1  , D2 = T.H(2)*cshrink1 ;
	      if (op==1)
		  r =  D0.x*u0 + D1.x*u1 + D2.x*u2 ;
	      else 
		  r =  D0.y*u0 + D1.y*u1 + D2.y*u2 ;
	  }
	//  cout << r << "\t";
	return r;
    }
    
    void TypeOfFE_P1ttdc1_::FB(const bool *whatd,const Mesh & ,const Triangle & K,const R2 & P1,RNMK_ & val) const
    {
      R2 P=Shrink1(P1);  
      
      //  const Triangle & K(FE.T);
      R2 A(K[0]), B(K[1]),C(K[2]);
      R l0=1-P.x-P.y,l1=P.x,l2=P.y; 
      
      if (val.N() <3) 
	  throwassert(val.N() >=3);
      throwassert(val.M()==1 );
      //  throwassert(val.K()==3 );
      
      val=0; 
      RN_ f0(val('.',0,op_id)); 
      
      if (whatd[op_id]) 
	{
	  f0[0] = l0;
	  f0[1] = l1;
	  f0[2] = l2;}
      if (whatd[op_dx] || whatd[op_dy])
	{
	  R2 Dl0(K.H(0)*cshrink1), Dl1(K.H(1)*cshrink1), Dl2(K.H(2)*cshrink1);
	  
	  if (whatd[op_dx]) 
	    {
	      RN_ f0x(val('.',0,op_dx)); 
	      f0x[0] = Dl0.x;
	      f0x[1] = Dl1.x;
	      f0x[2] = Dl2.x;
	    }
	  
	  if (whatd[op_dy]) {
	      RN_ f0y(val('.',0,op_dy)); 
	      f0y[0] = Dl0.y;
	      f0y[1] = Dl1.y;
	      f0y[2] = Dl2.y;
	  }
	}
    }
    
    
    
    void TypeOfFE_P2ttdc1_::FB(const bool *whatd,const Mesh & ,const Triangle & K,const R2 & P1,RNMK_ & val) const
    {
      R2 P=Shrink1(P1);  
      
      //  const Triangle & K(FE.T);
      R2 A(K[0]), B(K[1]),C(K[2]);
      R l0=1-P.x-P.y,l1=P.x,l2=P.y; 
      R l4_0=(4*l0-1),l4_1=(4*l1-1),l4_2=(4*l2-1); 
      
      //  throwassert(FE.N == 1);  
      throwassert( val.N()>=6);
      throwassert(val.M()==1);
      //  throwassert(val.K()==3 );
      
      val=0; 
      // --     
      if (whatd[op_id])
	{
	  RN_ f0(val('.',0,op_id)); 
	  f0[0] = l0*(2*l0-1);
	  f0[1] = l1*(2*l1-1);
	  f0[2] = l2*(2*l2-1);
	  f0[3] = 4*l1*l2; // oppose au sommet 0
	  f0[4] = 4*l0*l2; // oppose au sommet 1
	  f0[5] = 4*l1*l0; // oppose au sommet 3
	}
      if(  whatd[op_dx] || whatd[op_dy] || whatd[op_dxx] || whatd[op_dyy] ||  whatd[op_dxy])
	{
	  R2 Dl0(K.H(0)*cshrink1), Dl1(K.H(1)*cshrink1), Dl2(K.H(2)*cshrink1);
	  if (whatd[op_dx])
	    {
	      RN_ f0x(val('.',0,op_dx)); 
	      f0x[0] = Dl0.x*l4_0;
	      f0x[1] = Dl1.x*l4_1;
	      f0x[2] = Dl2.x*l4_2;
	      f0x[3] = 4*(Dl1.x*l2 + Dl2.x*l1) ;
	      f0x[4] = 4*(Dl2.x*l0 + Dl0.x*l2) ;
	      f0x[5] = 4*(Dl0.x*l1 + Dl1.x*l0) ;
	    }
	  
	  if (whatd[op_dy])
	    {  
		RN_ f0y(val('.',0,op_dy)); 
		f0y[0] = Dl0.y*l4_0;
		f0y[1] = Dl1.y*l4_1;
		f0y[2] = Dl2.y*l4_2;
		f0y[3] = 4*(Dl1.y*l2 + Dl2.y*l1) ;
		f0y[4] = 4*(Dl2.y*l0 + Dl0.y*l2) ;
		f0y[5] = 4*(Dl0.y*l1 + Dl1.y*l0) ;
	    }
	  
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
	  
	}
      
    }
    
    
    //
    // end ttdc1_
    // ------------------
    


// link with FreeFem++ 
static TypeOfFE_P1ttdc1_ P1dc1LagrangeP1dc1;
static TypeOfFE_P2ttdc1_ P2dc1LagrangeP2dc1;
// a static variable to add the finite element to freefem++
static AddNewFE  P1dcLagrange("P1dc1",&P1dc1LagrangeP1dc1); 
static AddNewFE  P2dcLagrange("P2dc1",&P2dc1LagrangeP2dc1); 
} // FEM2d namespace 


// --- fin -- 

