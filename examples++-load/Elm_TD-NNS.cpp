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
				      3+9, // nb coef to build interpolation
				      4, // np point to build interpolation
				      0)
    {  
      const double gauss1=(1.-sqrt(1./3.))/2;
      const double gauss2=1.-gauss1;

      const R2 Pt[] = { R2(1./3,1./3.),R2(0.5,0.5), R2(0,0.5),R2(0.5,0)}; 
      // for the 3 vertices 6 coef 
      P_Pi_h[0]=Pt[0];
	int kk=0;   
      for (int c=0;c<3;c++)
	{ 	  
	  pij_alpha[kk++]= IPJ(3+c,0,c);
	}
      // for 
      int p=1;
      for (int e=0;e<3;++e)
	{ // point d'integration sur l'arete e 
	  P_Pi_h[p]= Pt[p];
	  //	  cout <<"\n" <<  p << " --  " << P_Pi_h[p] << " ::  " << A << " " << B << endl;
	  pij_alpha[kk++]= IPJ(e,1+e,0); // coef = 0.5* l_e *ne_x * sge	  
	  pij_alpha[kk++]= IPJ(e,1+e,1); // coef = 0.5* l_e *ne_x * sge	  
	  pij_alpha[kk++]= IPJ(e,1+e,2); // coef = 0.5* l_e *ne_x * sge	  
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
	
        R2 N(T.Edge(i).perp()/ T.lenEdge(i));
        v[k++]= N.x*N.x; 
        v[k++]= N.y*N.x;
	v[k++]= N.y*N.y;
      }
  }
  
  void TypeOfFE_TD_NNS::FB(const bool * whatd,const Mesh & ,const Triangle & K,const R2 & P,RNMK_ & val) const
  {
    typedef double R;
    //R2 A(K[0]), B(K[1]),C(K[2]);
    R l0=1-P.x-P.y,l1=P.x,l2=P.y;
    const R c3= 1./3.;
    R ll3[3]={l0-c3,l1-c3,l2-c3};
    R ll[3]={l0,l1,l2};
    R2 Dl[3]=  {K.H(0), K.H(1), K.H(2) };
    R2 Rl[3]= { K.Edge(0)/ K.lenEdge(0),   K.Edge(1)/ K.lenEdge(1),   K.Edge(2)/ K.lenEdge(2)};   //  
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
    R B[3][3], BB[3][3];
    R cc = 3.;
    for(int j=0;j<3;++j)
      for(int i=0;i<3;++i)
	B[i][j]= S[i][j]*ll[j];
    
    for(int i=0;i<3;++i)
	for(int k=0;k<3;++k)
	  {  BB[i][k]=0.;	      
	      for(int j=0;j<3;++j)
		  BB[i][k] += cc*S[i][j]*ll[j]*S1[j][k];
	  }

    
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
	
	RNM_ f(val('.','.',op_id)); 
	for(int c=0;c<3;++c)
	    for(int i=0;i<3;++i){
		f(i,c)    = S[c][i]*ll3[i];
		f(i+3,c)  = BB[c][i];	      
	    }
      }
    if (wd[op_dx])
      {
	RNM_ f_x(val('.','.',op_dx));
	for(int i=0;i<3;++i)
	    for(int k=0;k<3;++k)
	      {  BB[i][k]=0.;	      
		  for(int j=0;j<3;++j)
		      BB[i][k] += cc*S[i][j]*Dl[j].x*S1[j][k];
	      }
	
	for(int c=0;c<3;++c)
	    for(int i=0;i<3;++i)	    
		f_x(i+3,c)  = BB[c][i];	      
	
	
      }
    
    if (wd[op_dy])
      {  
	  
	  RNM_ f_y(val('.','.',op_dy));
	  for(int i=0;i<3;++i)
	      for(int k=0;k<3;++k)
		{  BB[i][k]=0.	;      
		    for(int j=0;j<3;++j)
			BB[i][k] += cc*S[i][j]*Dl[j].y*S1[j][k];
		}
	  
	  for(int c=0;c<3;++c)
	      for(int i=0;i<3;++i)		    
		  f_y(i+3,c)  = BB[c][i];	      
	  
	  
	  
	  
      }
    
  }


// a static variable to add the finite element to freefem++
static TypeOfFE_TD_NNS Elm_TD_NNS;
static AddNewFE FE__TD_NNS("TDNNS",&Elm_TD_NNS); 
} // FEM2d namespace 
  // --- fin -- 

  
