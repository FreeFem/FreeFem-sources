#include "error.hpp"
#include "AFunction.hpp"
#include "rgraph.hpp"
using namespace std;  
#include "RNM.hpp"
#include "fem.hpp"
#include "FESpace.hpp"
#include "AddNewFE.h"
// Attention probleme de numerotation des inconnues
// -------------------------------------------------
// dans freefem, il y a un noeud par objets  sommet, arete, element.
// et donc la numerotation des dl dans l'element depend 
// de l'orientation des aretes
// 
/// ---------------------------------------------------------------
namespace  Fem2D {

 // ------ P3dc  
 class TypeOfFE_P3dcLagrange : public  TypeOfFE { public:  
   static const int k=3; 
   static const int ndf = (k+2)*(k+1)/2;
   static int Data[];
   static double Pi_h_coef[];
   static const int nn[10][3] ;
   static const int aa[10][3] ;
   static const int ff[10];
   static const int il[10];
   static const int jl[10];
   static const int kl[10];
    
   static const R2 G;
   static const R cshrink;
   static const R cshrink1;
   //  (1 -1/3)*
    
   static R2 Shrink(const R2& P){ return (P-G)*cshrink+G;}
   static R2 Shrink1(const R2& P){ return (P-G)*cshrink1+G;}
   
   TypeOfFE_P3dcLagrange(): TypeOfFE(3+2*3+1,1,Data,4,1,10,10,Pi_h_coef)
   {  
     static const  R2 Pt[10] = { 
       R2( 0/3. , 0/3. ) , 
       R2( 3/3. , 0/3. ) , 
       R2( 0/3. , 3/3. ) , 
       R2( 2/3. , 1/3. ) ,  
       R2( 1/3. , 2/3. ) , 
       R2( 0/3. , 2/3. ) , 
       R2( 0/3. , 1/3. ) , 
       R2( 1/3. , 0/3. ) , 
       R2( 2/3. , 0/3. ) , 
       R2( 1/3. , 1/3. ) } ;
     //   3,4,5,6,7,8
     for (int i=0;i<NbDoF;i++) {
       pij_alpha[i]= IPJ(i,i,0);
       P_Pi_h[i]=Shrink(Pt[i]); }

   }

  
   void FB(const bool * whatd, const Mesh & Th,const Triangle & K,const R2 &P, RNMK_ & val) const;
   /*   void Pi_h_alpha(const baseFElement & K,KN_<double> & v) const
   {
     for (int i=0;i<10;++i)
       v[i]=1;
     int e0=K.EdgeOrientation(0);
     int e1=K.EdgeOrientation(1);
     int e2=K.EdgeOrientation(2);
     int ooo[6]={e0,e0,e1,e1,e2,e2};
     int iii[6]; 
     int jjj[6];
     for(int i=0;i<6;++i)
       { 
	 iii[i]= 3+2*i; // si  orient = 1
         jjj[i]= 4+2*i; // si orient = -1
       }
     for(int i=0;i<6;++i)
       if(ooo[i]==1) v[jjj[i]]=0;
       else v[iii[i]]=0;
          
       }*/
 } ;


  const R2 TypeOfFE_P3dcLagrange::G(1./3.,1./3.);   
  const R TypeOfFE_P3dcLagrange::cshrink=1-1e-2;   
  const R TypeOfFE_P3dcLagrange::cshrink1=1./TypeOfFE_P3dcLagrange::cshrink;   

  //                     on what     nu df on node node of df    
  int TypeOfFE_P3dcLagrange::Data[]={
    6,6,6, 6,6, 6,6, 6,6, 6,    // the support number  of the node of the df  
    0,1,2,3,4,5,6,7,8,9,     // the number of the df on  the node  
    0,0,0,0,0,0,0,0,0,0,     //  the node of the df 
    0,0,0,0,0,0,0,0,0,0,     //  the df come from which FE (generaly 0)
    0,1,2,3,4,5,6,7,8,9,    //  which are de df on sub FE
    0,                      // for each compontant $j=0,N-1$ it give the sub FE associated 
    0, 10 };
  double TypeOfFE_P3dcLagrange::Pi_h_coef[]={1.,1.,1.,1.,1.,1.,1.,1.,1.,1.};
  
  void TypeOfFE_P3dcLagrange::FB(const bool * whatd,const Mesh & ,const Triangle & K,const R2 & P1,RNMK_ & val) const
  {
    R2 P=Shrink1(P1);
    R2 A(K[0]), B(K[1]),C(K[2]);
    R l0=1-P.x-P.y,l1=P.x,l2=P.y; 
    R L[3]={l0*k,l1*k,l2*k};
    throwassert( val.N()>=10);
    throwassert(val.M()==1);
    // Attention il faut renumeroter les fonction de bases
    //   car dans freefem++, il y a un node par sommet, arete or element
    //   et la numerotation naturelle  mais 2 noud pas arete
    //   donc p est la perumation
    //   echange de numerotation si les arete sont dans le mauvais sens 
    int p[10];
    for(int i=0;i<10;++i)
      p[i]=i;
    /*
    if(K.EdgeOrientation(0) <0) Exchange(p[3],p[4]);// 3,4
    if(K.EdgeOrientation(1) <0) Exchange(p[5],p[6]);// 5,6 
    if(K.EdgeOrientation(2) <0) Exchange(p[7],p[8]);// 7,8
    */
    //cout << KN_<int>(p,10) <<endl;
    val=0; 
    /*
    //  les fonction de base du Pk Lagrange sont 
    //  
    //
    */
// -- 
    
    
    if (whatd[op_id])
      {
	RN_ f0(val('.',0,op_id)); 
	for (int df=0;df<ndf;df++)
	  {
	    int pdf=p[df];
	    R f=1./ff[df];
	    
	    for( int i=0;i<k;++i)
	      {
		f *= L[nn[df][i]]-aa[df][i];
		//cout <<  L[nn[df][i]]-aa[df][i]<< " ";
	      }
	    f0[pdf] = f;
	    //cout << pdf<< " " << df << " f " <<f <<endl;
	  }
	//cout <<" L " << L[0] << " " << L[1] << " " << L[2] << endl;
	//cout << ndf << " nbf = "<< f0 <<endl;
      }
    
    
    if(  whatd[op_dx] || whatd[op_dy] || whatd[op_dxx] || whatd[op_dyy] ||  whatd[op_dxy])
      {
	R ks=k*cshrink1;
	R2 D[]={K.H(0)*ks, K.H(1)*ks,K.H(2)*ks };
	if (whatd[op_dx] || whatd[op_dy] )
	  {
	    for (int df=0;df<ndf;df++)
	      {
		int pdf=p[df];
		R fx=0.,fy=0.,f=1./ff[df];
		for( int i=0;i<k;++i) 
		  {
		    int n= nn[df][i];
		    R Ln=L[n]-aa[df][i];
		    fx= fx*Ln+f*D[n].x;
		    fy= fy*Ln+f*D[n].y;
		    f = f*Ln;
		  } 
		if(whatd[op_dx]) val(pdf,0,op_dx)=fx;
		if(whatd[op_dy]) val(pdf,0,op_dy)=fy;		
	      } 
	  }
	
	if (whatd[op_dyy] ||whatd[op_dxy]|| whatd[op_dxx] )
	  {  
	    for (int df=0;df<ndf;df++)
	      {
		int pdf=p[df];		
		R fx=0.,fy=0.,f=1./ff[df];
		R fxx=0.,fyy=0.,fxy=0.;
		for( int i=0;i<k;++i) 
		  {
		    int n= nn[df][i];
		    R Ln=L[n]-aa[df][i];
		    fxx=fxx*Ln+2.*fx*D[n].x;
		    fyy=fyy*Ln+2.*fy*D[n].y;
		    fxy=fxy*Ln+fx*D[n].y+fy*D[n].x;
		    fx= fx*Ln+f*D[n].x;
		    fy= fy*Ln+f*D[n].y;
		    f = f*Ln;
		  } 
		if(whatd[op_dxx]) val(pdf,0,op_dxx)=fxx;
		if(whatd[op_dyy]) val(pdf,0,op_dyy)=fyy;
		if(whatd[op_dxy]) val(pdf,0,op_dxy)=fxy;
	      }
	  }
	
      }
  }
#include "Element_P3dc.hpp"
  


// link with FreeFem++ 
static TypeOfFE_P3dcLagrange P3dcLagrangeP3dc;
// a static variable to add the finite element to freefem++
static AddNewFE  P3dcLagrange("P3dc",&P3dcLagrangeP3dc); 
} // FEM2d namespace 


// --- fin -- 

