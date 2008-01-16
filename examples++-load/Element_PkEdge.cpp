#include "error.hpp"
#include "AFunction.hpp"
#include "rgraph.hpp"
using namespace std;  
#include "RNM.hpp"
#include "fem.hpp"
#include "FESpace.hpp"
#include "QuadratureFormular.hpp"
#include "AddNewFE.h"
// Attention probleme de numerotation des inconnues
// -------------------------------------------------
// dans freefem, il y a un noeud par objets  sommet, arete, element.
// et donc la numerotation des dl dans l'element depend 
// de l'orientation des aretes
// 
/// ---------------------------------------------------------------
namespace  Fem2D {
  struct  InitTypeOfFE_PkEdge
  {
    int k;//  order poly on edge
    int npe; // nb point on edge
    int ndf; // nb dof
    KN<R> X;  //  point on edge
    //    KN<R> Pi_h_coef; // 1
    KN<int> Data; // data of TypeOfFE
    InitTypeOfFE_PkEdge(int KK) 
      :k(KK),npe(k+1),ndf(3*npe),X(npe),Data( 5*ndf+2)
    {
      //Pi_h_coef=1.;
      const QuadratureFormular1d QF(k+1);
      for (int i=0;i<k;++i)
	X[i]=QF[i].x;
      int j=0;
      int o[6];;
      o[0]=0;
      for(int i=1;i<6;++i)
	o[i]=o[i-1]+ndf;
      for(int df=0;df<ndf;++df)
	{
	  int e= df/3;
	  int n=df%3;
          Data[o[0]+df]=3+e;
          Data[o[1]+df]=n;
          Data[o[2]+df]=e;
          Data[o[3]+df]=0;
          Data[o[4]+df]=df;
	}
      Data[o[5]] =0;
      Data[o[5]=1] =0;
    }
  };

  class TypeOfFE_PkEdge :public InitTypeOfFE_PkEdge, public  TypeOfFE { 
  public:  
    static double Pi_h_coef[];
    
    
    TypeOfFE_PkEdge(int KK)
      :  InitTypeOfFE_PkEdge(KK),
	 TypeOfFE(ndf,1,Data,2+k,1,ndf*2,ndf,0)
    {  
      ffassert(k<2);
      int kkk=0;
      for (int i=0;i<NbDoF;i++) 
	{
	  int e= i/npe;
	  int j= i%npe;
	  int ii= e*npe+npe-j-1;
	  R2 A(TriangleHat[VerticesOfTriangularEdge[i][0]]);
	  R2 B(TriangleHat[VerticesOfTriangularEdge[i][1]]);	 
	  pij_alpha[kkk++]= IPJ(i,i,0);
	  pij_alpha[kkk++]= IPJ(i,ii,0);
	  P_Pi_h[i]= A*(1.-X[i])+ B*(X[i]);// X=0 => A  X=1 => B;       
	}
    }

   void Pi_h_alpha(const baseFElement & K,KN_<double> & v) const
   {
     int kkk=0;
     for (int e=0;e<3;++e)
       {
	 int i0=1;
	 if( !K.EdgeOrientation(e))
	   i0=1-i0;
	 int i1=1-i0;
	 for(int p=0;p<npe;++p)
	   {
	     v[kkk+i0]=0;
	     v[kkk+i1]=1;
	     kkk+=2;
	   }
       }          
   }
    
    void FB(const bool * whatd, const Mesh & Th,const Triangle & K,const R2 &P, RNMK_ & val) const;
    } ;
  // ENDOFCLASS TypeOfFE_PkEdge    
    
    void TypeOfFE_PkEdge::FB(const bool * whatd,const Mesh & ,const Triangle & K,const R2 & P,RNMK_ & val) const
    {
      
      R2 A(K[0]), B(K[1]),C(K[2]);
      R l0=1-P.x-P.y,l1=P.x,l2=P.y; 
      R L[3]={l0*k,l1*k,l2*k};
      throwassert( val.N()>=10);
      throwassert(val.M()==1);
      int ee=0;
      if (L[0] <= min(L[1],L[2]) ) ee=0; // arete  
      else if  (L[1] <= min(L[0],L[2]) ) ee=1;
      else ee=2;
      int e3=ee*npe;
      int s=1-L[ee];
      R xe = L[VerticesOfTriangularEdge[ee][0]]/s;//  go from 0 to 1 on edge 
      if(K.EdgeOrientation(ee) <0) 
	xe = 1-xe;
      val=0; 
      if (whatd[op_id])
	{
	  RN_ f0(val('.',0,op_id)); 
	  for (int l=0;l<npe;l++)
	    {
	      int df= e3+l;
	      R f=1.;
	      for (int i=0;i<npe;++i)
		if(i != l) 
		  f *= (xe-X[i])/(X[l]-X[i]);
	      f0[df] = f;
	    }
	}
      
      
      if(  whatd[op_dx] || whatd[op_dy] || whatd[op_dxx] || whatd[op_dyy] ||  whatd[op_dxy])
	{
	  cerr << " TO DO ???  FH " << endl;
	  ffassert(0);
	}
    }

  
  // link with FreeFem++ 
  static TypeOfFE_PkEdge PkEdgeP1(1);
  static TypeOfFE_PkEdge PkEdgeP2(2);
  static TypeOfFE_PkEdge PkEdgeP3(3);
  static TypeOfFE_PkEdge PkEdgeP4(4);
  static TypeOfFE_PkEdge PkEdgeP5(5);
  // a static variable to add the finite element to freefem++
  static AddNewFE  P1Edge("P1edge",&PkEdgeP1); 
  static AddNewFE  P2Edge("P2edge",&PkEdgeP2); 
  static AddNewFE  P3Edge("P2edge",&PkEdgeP3); 
  static AddNewFE  P4Edge("P2edge",&PkEdgeP4); 
  static AddNewFE  P5Edge("P2edge",&PkEdgeP5); 
} // FEM2d namespace 


// --- fin -- 

