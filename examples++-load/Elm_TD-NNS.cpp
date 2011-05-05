//  MORLEY FINITE ELEMENT 
// F. Hecht  december 2005
// ----------------------------
//  the Polynomial space is $P¡2$, an the degree of freedom are
// the 3 values a the 3 vertices and the three 
// normal derivative at the middle at the tree the edges
// remark; 
//   to compute the interpolante, we need
//   the value , plus the value of the normal derivative
//   so I use the following hack, I say the is a tree dim vectorial 
//   finite element with give  the value, x derivative ,and the y derivative
//  Ref: chapter VII section 50  fig 50.2  of 
//   Ciarlet,   HandBook of Numerical Analysis, Volume II Finite elemet methodes (parts 1),
//             NORTH-HOLLAND 
// -----------------------  related files: 
//  to check  and validate  : testFETD_NNS.edp 
//  to get a real example   :  bilapMorley.edp
// ------------------------------------------------------------

#include "ff++.hpp"
#include "AddNewFE.h"

// #include "problem.hpp"
namespace  Fem2D {
  
  // ------ P2 TD_NNS  
  class TypeOfFE_TD_NNS : public  TypeOfFE { public:  
    static int Data[];
    // double Pi_h_coef[];
    
    TypeOfFE_TD_NNS(): TypeOfFE(3+3,
				  3,
				      Data,
				      2,
				      1,
				      3+6, // nb coef to build interpolation
				      6, // np point to build interpolation
				      0)
    {  
      const double gauss1=(1.-sqrt(1./3.))/2;
      const double gauss2=1.-gauss1;

      const R2 Pt[] = { R2(0,0),R2(1,0),R2(0,1),R2(0.5,0.5), R2(0,0.5),R2(0.5,0)}; 
      // for the 3 vertices 6 coef 
      int kk=0;
      for (int p=0;p<3;p++)
	{ 
	  P_Pi_h[p]=Pt[p];
	  pij_alpha[kk]= IPJ(kk,p,0);
	  kk++;
	}
      // for 
      int p=3;
      for (int e=0;e<3;++e)
	{ // point d'integration sur l'arete e 
	  P_Pi_h[p]= Pt[p];
	  //	  cout <<"\n" <<  p << " --  " << P_Pi_h[p] << " ::  " << A << " " << B << endl;
	  pij_alpha[kk++]= IPJ(3+e,p,1); // coef = 0.5* l_e *ne_x * sge	  
	  pij_alpha[kk++]= IPJ(3+e,p,2); // coef = 0.5* l_e *ne_y * sge	  	  
	  p++;
	}
      assert(P_Pi_h.N()==p);
      assert(pij_alpha.N()==kk);

     }
    void FB(const bool * whatd, const Mesh & Th,const Triangle & K,const R2 &P, RNMK_ & val) const;
    void Pi_h_alpha(const baseFElement & K,KN_<double> & v) const;
  } ;
  //                     on what     nu df on node node of df    
 int TypeOfFE_TD_NNS::Data[]={
   3,4,5,  6,6,6 ,//  support on what 
   0,0,0,  0,1,3,// df on node 
   0,1,2,  4,4,4,// th node of df 
   0,0,0,  0,0,0,//  df previou FE
   0,1,2,  3,4,5, //  which df on prev 
   0,0,0,
   0,0,0, 
   6,6,6 
};

void TypeOfFE_TD_NNS::Pi_h_alpha(const baseFElement & K,KN_<double> & v) const
  {
    const Triangle & T(K.T);
    int k=0;
    // coef pour les 3 sommets  fois le 2 composantes 
    for (int i=0;i<3;i++)
      v[k++]=1; 
    //   integration sur les aretes 
    for (int i=0;i<3;i++)
      {
	
        R2 N(T.Edge(i).perp());
	N  *= T.EdgeOrientation(i);
        v[k++]= N.x; 
        v[k++]= N.y;
      }
  }
  
  void TypeOfFE_TD_NNS::FB(const bool * whatd,const Mesh & ,const Triangle & K,const R2 & P,RNMK_ & val) const
  {
    typedef double R;
    R2 A(K[0]), B(K[1]),C(K[2]);
    R l0=1-P.x-P.y,l1=P.x,l2=P.y;
    cont R c3= 1./3.;
    R ll3[3]={l0-c3,l1-c3,l3-c3};
    R ll[3]={l0,l1,l3};
    R2 Dl[3]=  {K.H(0), K.H(1), K.H(2) };
    R2 Rl[3]= { K.E(0),   K.E(1),   K.E(2)};   //  
    /* bulle:
      $ B_i = ((Rot L_i+1 ) (Rot L_(i+2)' ))^s  L_i$
         s => symetrise ..
    */
    R S[3][3],S1[3][3];
    for(int i=0;i<3;++i)
      {
	int i1=(i+1)%3;
	int i2=(i+2)%3;
	S[0][i]= Rl[i1].x*Rl[i2].x;
	S[1][i]= (Rl[i1].x*Rl[i2].y+Rl[i1].y*Rl[i2].x)*0.5;
	S[2][i]= Rl[i1].y*Rl[i2].y;
      }
    // s[.] [i] = 
    { //     //  compute the inv of S 
    double det=S[0][0]*(S[1][1]*S[2][2]-S[2][1]*S[1][2])-S[0][1]*(S[1][0]*S[2][2]-S[1][2]*S[2][0])+S[0][2]*(S[1][0]*S[2][1]-S[1][1]*S[2][0]);//adjoin

    S1[0][0]=(S[1][1]*S[2][2]-S[2][1]*S[1][2])/det;
    S1[0][1]=-(S[1][0]*S[2][2]-S[1][2]*S[2][0])/det;
    S1[0][2]=(S[1][0]*S[2][1]-S[2][0]*S[1][1])/det;
    S1[1][0]=-(S[0][1]*S[2][2]-S[0][2]*S[2][1])/det;
    S1[1][1]=(S[0][0]*S[2][2]-S[0][2]*S[2][0])/det;
    S1[1][2]=-(S[0][0]*S[2][1]-S[2][0]*S[0][1])/det;
    S1[2][0]=(S[0][1]*S[1][2]-S[0][2]*S[1][1])/det;
    S1[2][1]=-(S[0][0]*S[1][2]-S[1][0]*S[0][2])/det;
    S1[2][2]=(S[0][0]*S[1][1]-S[1][0]*S[0][1])/det;
    }
    R B[3][3], BB[3][3]={ 0.,0.,0., 0.,0.,0., 0.,0.,0. };
    R cc = 3.;
    for(int j=0;j<3;++j)
      for(int i=0;i<3;++i)
	B[i][j]= S[i][j]*ll[i];
    
      for(int i=0;i<3;++i)
	for(int k=0;k<3;++k)
	  for(int j=0;j<3;++j)
	    BB[i][k] += cc*B[i]*[j]*S1[j]*[k];

    
    // the basic function are
    //  the space S_i * ( a_i+b_1lambda_i) 
    //  the base :
    //  tree egde function: 
    //   coefe*  S_i*( lambda_i - 1/3)   :  zero a barycenter 
    //   coefk*BB_i ,withh     B_i =   (S_i * lambda_i),   BB = B * S1 , ok because lambda_i = 1/3 a bary
    //  so BB_ij = coefk/3  delta_ij  at G the barycenter.
    // 
    KN<bool> wd(KN_<const bool>(whatd,last_operatortype));
    val=0; 
    
    throwassert( val.N()>=6);
    throwassert(val.M()==3);
    
    
    val=0; 
    
    
    if (wd[op_id])
      {
	
	RNM_ f(val('.','.',,op_id)); 
	for(int c=0;c<3;++c)
	  for(int i=0;i<3;++i){
	    f(i,c)    = S[c][i]*ll3[i];
	    f(i+3,c)  = BB[c][i];
	      
	      }
  
	if (wd[op_dx])
	  {
	    RN_ fx(val('.',wd_j[op_dx],wd_op[op_dx])); 
	    
	    fx[0] = Dl[0].x;
	    fx[1] = Dl[1].x;
	    fx[2] = Dl[2].x;
	    
	    fx[3] = dl3*Dl[0].x; 
	    fx[4] = dl4*Dl[1].x; 
	    fx[5] = dl5*Dl[2].x;
	    
	  }
	
	if (wd[op_dy])
	  {  
	    //      RN_ fy(val('.',0,op_dy)); 
	    RN_ fy(val('.',wd_j[op_dy],wd_op[op_dy]));       
	    fy[0] = Dl[0].y;
	    fy[1] = Dl[1].y;
	    fy[2] = Dl[2].y;
	    
	    fy[3] = dl3*Dl[0].y; 
	    fy[4] = dl4*Dl[1].y; 
	    fy[5] = dl5*Dl[2].y;
      
    }
  
  }


// a static variable to add the finite element to freefem++
static TypeOfFE_TD_NNS Elm_TD_NNS;
static AddNewFE FE__TD_NNS("TDNNS",&Elm_TD_NNS); 
} // FEM2d namespace 
  // --- fin -- 

  
