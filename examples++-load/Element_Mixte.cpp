//  add some mixte finite element  
//  RT1 , RT1Ortho  ( Neledec type I order 1)
//  BDM1 (Brezzi, Douglas, Marini) ,  , BDM1rtho  ( Neledec type Iï order 1)
//   TDNNS  
// DM
// F. Hecht  May 2011
// ----------------------------

//             NORTH-HOLLAND 
// -----------------------  related files: 
//  to check  and validate  : test_MIXTE.edp 
//  to get a real example   :  bilapMorley.edp
// ------------------------------------------------------------
//ff-c++-LIBRARY-dep:   lapack
//ff-c++-LIBRARY-dep:   blas

#include "ff++.hpp"
#include "AddNewFE.h"

#ifdef __LP64__
typedef int intblas;
typedef int integer;
#else
typedef long intblas;
typedef long integer;
#endif


typedef integer  logical;
typedef float   LAPACK_real;
typedef double   doublereal;
typedef logical  (* L_fp)();
typedef integer      ftnlen;

typedef complex<float> LAPACK_complex;
typedef complex<double> doublecomplex;
typedef void VOID; 
#define complex LAPACK_complex 
#define real LAPACK_real 

#include "clapack.h"
#undef real
#undef complex 
// #include "problem.hpp"
namespace  Fem2D {
    
    // ------ P2 TD_NNS  
    class TypeOfFE_TD_NNS0 : public  TypeOfFE { public:  
	static int Data[];
	// double Pi_h_coef[];
	
	TypeOfFE_TD_NNS0(): TypeOfFE(3,
				     3,
				     Data,
				     1,
				     1,
				     9, // nb coef to build interpolation
				     3, // np point to build interpolation
				     0)
	{  
	    
	    
	    const R c3=1./3.;
	    const R2 Pt[] = {R2(0.5,0.5), R2(0,0.5),R2(0.5,0)}; 
	    // for the 3 vertices 6 coef 
	    // P_Pi_h[0]=Pt[0];
	    int kk=0;   
	    
	    for (int e=0;e<3;++e)
	      { // point d'integration sur l'arete e 
		  P_Pi_h[e]= Pt[e];
		  pij_alpha[kk++]= IPJ(e,e,0);   
		  pij_alpha[kk++]= IPJ(e,e,1); 	  
		  pij_alpha[kk++]= IPJ(e,e,2);   
		  
	      }
	    assert(P_Pi_h.N()==3);
	    assert(pij_alpha.N()==kk);
	    
	}
	void FB(const bool * whatd, const Mesh & Th,const Triangle & K,const R2 &P, RNMK_ & val) const;
	void Pi_h_alpha(const baseFElement & K,KN_<double> & v) const;
    } ;
    //                     on what     nu df on node node of df    
    int TypeOfFE_TD_NNS0::Data[]={
	3,4,5,  //  support on what 
	0,0,0,  // df on node 
	0,1,2,  // th node of df 
	0,0,0,  //  df previou FE
	0,1,2,  //  which df on prev 
	0,0,0,
	0,0,0, 
	3,3,3 
    };
    
    void TypeOfFE_TD_NNS0::Pi_h_alpha(const baseFElement & K,KN_<double> & v) const
    {
      const Triangle & T(K.T);
      int k=0;
      //   integration sur les aretes 
      for (int i=0;i<3;i++)
	{
	  
	  R2 N(T.Edge(i).perp());
	  v[k++]= N.x*N.x; 
	  v[k++]= 2*N.y*N.x;
	  v[k++]= N.y*N.y;
	}
      assert(k==9); 
    }
    
    void TypeOfFE_TD_NNS0::FB(const bool * whatd,const Mesh & ,const Triangle & K,const R2 & P,RNMK_ & val) const
    {
      typedef double R;
      //R2 A(K[0]), B(K[1]),C(K[2]);
      R l0=1-P.x-P.y,l1=P.x,l2=P.y;
      const R c3= 1./3.;
      R ll3[3]={l0-c3,l1-c3,l2-c3};
      R ll[3]={l0,l1,l2};
      R2 Dl[3]=  {K.H(0), K.H(1), K.H(2) };
      /* if T_i=Edge(i) ,N=T.perp :
       N_i' T_j T_k' N_i =0  if i=j or i=k
       N_i' T_j = det(T_i,T_j) = aire(K) 
       */
      R cK= 2* K.area; 
      R2 Rl[3]= { K.Edge(0)/cK,   K.Edge(1)/ cK,   K.Edge(2)/ cK};   //  
      /* bulle:
       $ B_i = ((Rot L_i+1 ) (Rot L_(i+2)' ))^s  L_i$
       s => symetrise ..
       */
      R S[3][3],S1[3][3];
      for(int i=0;i<3;++i)
	{
	  int i1=(i+1)%3;
	  int i2=(i+2)%3;
	  S[0][i]= -Rl[i1].x*Rl[i2].x;
	  S[1][i]= -(Rl[i1].x*Rl[i2].y+Rl[i1].y*Rl[i2].x)*0.5;
	  S[2][i]= -Rl[i1].y*Rl[i2].y;
	  
	}
      val=0; 
      KN<bool> wd(KN_<const bool>(whatd,last_operatortype));
      
      if (wd[op_id])
	{
	  for(int c=0;c<3;++c)
	      for(int i=0;i<3;++i){
		  val(i,c,op_id)    = S[c][i]; //  (c3-ll[i])/c3 
	      }
	}
      
      /*
       // s[.] [i] = 
       { //     //  compute the inv of S with lapack 
       for(int j=0;j<3;++j)
       for(int i=0;i<3;++i)
       S1[i][j]=S[i][j];
       
       int N=3,LWORK = 9;
       double WORK[9] ;
       int INFO,IPIV[4];
       
       dgetrf_(&N,&N,&(S1[0][0]),&N,IPIV,&INFO);
       ffassert(INFO==0);
       dgetri_(&N,&(S1[0][0]),&N,IPIV,WORK,&LWORK,&INFO);
       ffassert(INFO==0);
       
       }
       R B[3][3], BB[3][3];
       R cc = 3./K.area; 
       for(int j=0;j<3;++j)
       for(int i=0;i<3;++i)
       B[i][j]= S[i][j]*ll[j];
       
       for(int i=0;i<3;++i)
       for(int k=0;k<3;++k)
       {  BB[i][k]=0.;	      
       for(int j=0;j<3;++j)
       BB[i][k] += cc*S[i][j]*ll[j]*S1[j][k];
       }
       if(verbosity>1000)
       {
       
       cout << endl;
       cout <<  Rl[0] << " "<< Rl[1]  << ",  " <<Rl[2] << endl;	
       for(int i=0;i<3;++i)
       cout << " *****    " << BB[i][0] << " " << BB[i][1] << " " << BB[i][2] << "\t l " << ll[i] << endl;
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
       for(int c=0;c<3;++c)
       for(int i=0;i<3;++i){
       val(i,c,op_id)    = S[c][i]*(c3-ll[i])/c3; //  (c3-ll[i])/c3 
       val(i+3,c,op_id)  = BB[c][i];	      
       }
       }
       if (wd[op_dx])
       {
       for(int i=0;i<3;++i)
       for(int k=0;k<3;++k)
       {  BB[i][k]=0.;	      
       for(int j=0;j<3;++j)
       BB[i][k] += cc*S[i][j]*Dl[j].x*S1[j][k];
       }
       
       for(int c=0;c<3;++c)
       for(int i=0;i<3;++i)
       {
       val(i  ,c,op_dx)    = -S[c][i]*Dl[i].x/c3;
       val(i+3,c,op_dx)  = BB[c][i];	      
       
       }
       
       
       }
       
       if (wd[op_dy])
       {  
       
       for(int i=0;i<3;++i)
       for(int k=0;k<3;++k)
       {  BB[i][k]=0.	;      
       for(int j=0;j<3;++j)
       BB[i][k] += cc*S[i][j]*Dl[j].y*S1[j][k];
       }
       
       for(int c=0;c<3;++c)
       for(int i=0;i<3;++i)
       {
       val(i  ,c,op_dy)    = -S[c][i]*Dl[i].y/c3;
       val(i+3,c,op_dy)  = BB[c][i];	      
       
       }
       
       
       }	  
       
       
       */
      
    }
    
    
    // ------ P2 TD_NNS  
    class TypeOfFE_TD_NNS1 : public  TypeOfFE { public:  
	static int Data[];
	// double Pi_h_coef[];
	const QuadratureFormular1d & QFE;
	const  GQuadratureFormular<R2> & QFK;
	
	TypeOfFE_TD_NNS1(): TypeOfFE(3*2+3,
				     3,
				     Data,
				     2,
				     1,
				     3+6*3*QF_GaussLegendre2.n, // nb coef to build interpolation
				     QuadratureFormular_T_1.n+3*QF_GaussLegendre2.n, // np point to build interpolation
				     0),
	QFE(QF_GaussLegendre2), QFK(QuadratureFormular_T_1)	
	{  
	    
	    
	    int kk=0,kp=0;
	    for(int p=0;p<QFK.n;++p)
	      {
		P_Pi_h[kp++]=QFK[p];
	        for (int c=0;c<3;c++)
		    pij_alpha[kk++]= IPJ(3*2+c,p,c);
	      }
	    
	    for (int e=0;e<3;++e) 
	      {
		for(int p=0;p<QFE.n;++p)
		  {
		    R2 A(TriangleHat[VerticesOfTriangularEdge[e][0]]);
		    R2 B(TriangleHat[VerticesOfTriangularEdge[e][1]]);
		    P_Pi_h[kp++]= B*(QFE[p].x)+ A*(1.-QFE[p].x);// X=0 => A  X=1 => B;       
		    
		  }
	      }
	    
	    
	    for (int e=0;e<3;++e)
		for(int p=0;p<QFE.n;++p)	
		  { 
		      int pp=QFK.n+ e*QFE.n+p;
		      pij_alpha[kk++]= IPJ(2*e+0,pp,0);  
		      pij_alpha[kk++]= IPJ(2*e+1,pp,0);  
		      pij_alpha[kk++]= IPJ(2*e+0,pp,1);  
		      pij_alpha[kk++]= IPJ(2*e+1,pp,1);  
		      pij_alpha[kk++]= IPJ(2*e+0,pp,2);  
		      pij_alpha[kk++]= IPJ(2*e+1,pp,2);  
		      
		      
		  }
	    ffassert(P_Pi_h.N()==kp);
	    ffassert(pij_alpha.N()==kk);
	    
	}
	void FB(const bool * whatd, const Mesh & Th,const Triangle & K,const R2 &P, RNMK_ & val) const;
	void Pi_h_alpha(const baseFElement & K,KN_<double> & v) const;
    } ;
    //                     on what     nu df on node node of df    
    int TypeOfFE_TD_NNS1::Data[]={
	3,3, 4,4, 5,5, 6,6,6,//  support on what 
	0,1, 0,1, 0,1, 0,1,2, // df on node 
	0,0, 1,1, 2,2, 3,3,3,  // th node of df 
	0,0, 0,0, 0,0, 0,0,0, //  df previou FE
	0,1, 2,3, 4,5, 6,7,8, //  which df on prev 
	0,0,0,
	0,0,0, 
	9,9,9 
    };
    
    void TypeOfFE_TD_NNS1::Pi_h_alpha(const baseFElement & K,KN_<double> & v) const
    {
      
      const Triangle & T(K.T);
      int k=0;
      // coef pour les 3 sommets  fois le 2 composantes 
      // for (int i=0;i<3;i++)
      for (int p=0;p<QFK.n;++p) 	
	{ // wrong ... 
	    // 3 -1 -1    
	    // -1 3 -1    / 4  
	    // -1 -1 3 
	    /*R l[3]; QFK[p].toBary(l);
	     R c0 = QFK[p].a * T.area/4*(3*l[0]-l[1]-l[2]);
	     R c1 = QFK[p].a * T.area/4*(l[0]+3*l[1]-l[2]);
	     R c2 = QFK[p].a * T.area/4*(l[0]-l[1]+3*l[2]);*/
	    R cc=QFK[p].a * T.area;
	    v[k++]=cc; // * T.area; 
	    v[k++]=cc;// * T.area; 
	    v[k++]=cc;// * T.area; 
	}
      //   integration sur les aretes 
      for (int i=0;i<3;i++)
	{
	  R s = T.EdgeOrientation(i) ;
	  for (int p=0;p<QFE.n;++p) 
	    {
	      R l0 = QFE[p].x, l1 = 1-QFE[p].x;
	      R p0= (2*l0-l1)*2;// poly othogonaux to \lambda_1
	      R p1= (2*l1-l0)*2;// poly othogonaux to \lambda_0
	      R cc1 = p0*QFE[p].a; // 
	      R cc0 = p1*QFE[p].a; //
	      if(s<0) Exchange(cc1,cc0); // exch lambda0,lambda1
	      
	      
	      R2 N(T.Edge(i).perp());
	      v[k++]= cc0*N.x*N.x;
	      v[k++]= cc1*N.x*N.x;    
	      v[k++]= cc0*2*N.y*N.x;
	      v[k++]= cc1*2*N.y*N.x;    
	      v[k++]= cc0*N.y*N.y;
	      v[k++]= cc1*N.y*N.y;
	    }
	}
      ffassert(pij_alpha.N()==k);    
    } 
    void TypeOfFE_TD_NNS1::FB(const bool * whatd,const Mesh & ,const Triangle & K,const R2 & P,RNMK_ & val) const
    {
      typedef double R;
      R l0=1-P.x-P.y,l1=P.x,l2=P.y;
      const R c3= 1./3.;
      R ll3[3]={l0-c3,l1-c3,l2-c3};
      R ll[3]={l0,l1,l2};
      R2 Dl[3]=  {K.H(0), K.H(1), K.H(2) };
      /* if T_i=Edge(i) ,N=T.perp :
       N_i' T_j T_k' N_i =0  if i=j or i=k
       N_i' T_j = det(T_i,T_j) = aire(K) 
       */
      R cK= 2* K.area; 
      R2 Rl[3]= { K.Edge(0)/cK,   K.Edge(1)/ cK,   K.Edge(2)/ cK};   //  
      /* bulle:
       $ B_i = ((Rot L_i+1 ) (Rot L_(i+2)' ))^s  L_i$
       s => symetrise ..
       */
      R S[3][3],S1[3][3];
      int ei0[3]={1,2,0};
      int ei1[3]={ 2,0,1};
      
      for(int i=0;i<3;++i)
	{
	  if(K.EdgeOrientation(i) < 0) Exchange(ei0[i],ei1[i]);
	  int i1=(i+1)%3;
	  int i2=(i+2)%3;
	  S[0][i]= -Rl[i1].x*Rl[i2].x;
	  S[1][i]= -(Rl[i1].x*Rl[i2].y+Rl[i1].y*Rl[i2].x)*0.5;
	  S[2][i]= -Rl[i1].y*Rl[i2].y;
	  
	}
      // s[.] [i] = 
      { //     //  compute the inv of S with lapack 
	  for(int j=0;j<3;++j)
	      for(int i=0;i<3;++i)
		  S1[i][j]=S[i][j];
	  
	  int N=3,LWORK = 9;
	  double WORK[9] ;
	  int INFO,IPIV[4];
	  
	  dgetrf_(&N,&N,&(S1[0][0]),&N,IPIV,&INFO);
	  ffassert(INFO==0);
	  dgetri_(&N,&(S1[0][0]),&N,IPIV,WORK,&LWORK,&INFO);
	  ffassert(INFO==0);
	  
      }
      R B[3][3], BB[3][3];
      R cc = 3./K.area; 
      for(int j=0;j<3;++j)
	  for(int i=0;i<3;++i)
	      B[i][j]= S[i][j]*ll[j];
      
      for(int i=0;i<3;++i)
	  for(int k=0;k<3;++k)
	    {  BB[i][k]=0.;	      
		for(int j=0;j<3;++j)
		    BB[i][k] += cc*S[i][j]*ll[j]*S1[j][k];
	    }
      if(verbosity>1000)
	{
	  
	  cout << endl;
	  cout <<  Rl[0] << " "<< Rl[1]  << ",  " <<Rl[2] << endl;	
	  for(int i=0;i<3;++i)
	      cout << " *****    " << BB[i][0] << " " << BB[i][1] << " " << BB[i][2] << "\t l " << ll[i] << endl;
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
	  for(int c=0;c<3;++c)
	      for(int i=0;i<3;++i){
		  
		  val(2*i,c,op_id)      = S[c][i]*(ll[ei0[i]]-ll[i]);  
		  val(2*i+1,c,op_id)    = S[c][i]*(ll[ei1[i]]-ll[i]);  
		  val(i+6,c,op_id)  = BB[c][i];	      
	      }
	}
      if (wd[op_dx])
	{
	  for(int i=0;i<3;++i)
	      for(int k=0;k<3;++k)
		{  BB[i][k]=0.;	      
		    for(int j=0;j<3;++j)
			BB[i][k] += cc*S[i][j]*Dl[j].x*S1[j][k];
		}
	  
	  for(int c=0;c<3;++c)
	      for(int i=0;i<3;++i)
		{
		  val(2*i,c,op_dx)      = S[c][i]*(Dl[ei0[i]].x-Dl[i].x); 
		  val(2*i+1,c,op_dx)    = S[c][i]*(Dl[ei1[i]].x-Dl[i].x); 
		  val(i+6,c,op_dx)  = BB[c][i];	      
		  
		}
	  
	  
	}
      
      if (wd[op_dy])
	{  
	    
	    for(int i=0;i<3;++i)
		for(int k=0;k<3;++k)
		  {  BB[i][k]=0.	;      
		      for(int j=0;j<3;++j)
			  BB[i][k] += cc*S[i][j]*Dl[j].y*S1[j][k];
		  }
	    
	    for(int c=0;c<3;++c)
		for(int i=0;i<3;++i)
		  {
		    val(2*i,c,op_dy)      = S[c][i]*(Dl[ei0[i]].y-Dl[i].y); 
		    val(2*i+1,c,op_dy)    = S[c][i]*(Dl[ei1[i]].y-Dl[i].y); 
		    val(i+6,c,op_dy)  = BB[c][i];	      
		    
		  }
	    
	    
	}	  
      
      
      
      
    }
    
    
    struct  InitTypeOfRTk_2d
  {
    
    int k;//  order poly on edge
    int ndfi;// nb of internal dof 
    int npe; // nb point on edge
    int ndf; // nb dof
    
    KN<R> X;  //  point on edge
    //    KN<R> Pi_h_coef; // 1
    KN<int> Data; // data of TypeOfFE
    const QuadratureFormular1d QFE;
    const  GQuadratureFormular<R2> & QFK;
    InitTypeOfRTk_2d(int KK) 
    :k(KK),ndfi((k+1)*(k)), npe(k+1),ndf(3*npe+ndfi),Data( 5*ndf+6),
    QFE(-1+2*npe,npe,GaussLegendre(npe),true), QFK(QuadratureFormular_T_5)
    {
      int j=0;
      int ndfe=ndf-ndfi; // 
      int o[6];
      o[0]=0;
      for(int i=1;i<6;++i)
	  o[i]=o[i-1]+ndf;
      for(int df=0;df<ndf;++df)
	{
	  if( df < ndfe)
	    {
	      int e= df/npe;
	      int n= df%npe;
	      Data[o[0]+df]=3+e;
	      Data[o[1]+df]=n;
	      Data[o[2]+df]=e;
	      Data[o[3]+df]=0;
	      Data[o[4]+df]=df;
	    }
	  else {
	      int n= df-ndfe;
	      Data[o[0]+df]=6;
	      Data[o[1]+df]=n;
	      Data[o[2]+df]=3;
	      Data[o[3]+df]=0;
	      Data[o[4]+df]=df;	      
	  }
	  
	}
      Data[o[5]+0] =0;
      Data[o[5]+1] =0 ;
      Data[o[5]+2] =0;
      Data[o[5]+3] =0;
      Data[o[5]+4] =ndf ;
      Data[o[5]+5] =ndf;// end_dfcomp 
      
    }
  };
    
    class TypeOfFE_RT1_2d :public InitTypeOfRTk_2d, public  TypeOfFE { 
    public:  
	static double Pi_h_coef[];
	bool Ortho;
	
	TypeOfFE_RT1_2d(bool ortho)
	:  InitTypeOfRTk_2d(1),
	TypeOfFE(ndf,2,Data,2,1,
		 2*2*3*QFE.n+QFK.n*2,// nb coef mat interpole
		 3*QFE.n+QFK.n, // nb P interpolation 
		 0),
	Ortho(ortho)
	{  
	    //      cout << " Pk = " << k << endl;
	    int kkk=0,i=0;
	    for (int e=0;e<3;++e) 
		for (int p=0;p<QFE.n;++p) 
		  {
		    R2 A(TriangleHat[VerticesOfTriangularEdge[e][0]]);
		    R2 B(TriangleHat[VerticesOfTriangularEdge[e][1]]);
		    
		    pij_alpha[kkk++]= IPJ(2*e,i,0);
		    pij_alpha[kkk++]= IPJ(2*e,i,1);
		    pij_alpha[kkk++]= IPJ(2*e+1,i,0);
		    pij_alpha[kkk++]= IPJ(2*e+1,i,1);
		    
		    
		    P_Pi_h[i++]= B*(QFE[p].x)+ A*(1.-QFE[p].x);// X=0 => A  X=1 => B;       
		  }
	    int i6=6,i7=7;
	    if(Ortho) Exchange(i6,i7); // x,y -> -y, x 
	    for (int p=0;p<QFK.n;++p) 
	      {
		pij_alpha[kkk++]= IPJ(i6,i,0);
		pij_alpha[kkk++]= IPJ(i7,i,1);
		P_Pi_h[i++]= QFK[p]; 
	      }
	    //cout << kkk << " kkk == " << this->pij_alpha.N() << endl;
	    //cout << i << "  ii == " << this->P_Pi_h.N() << endl;
	    ffassert(kkk==this->pij_alpha.N());
	    ffassert(i==this->P_Pi_h.N() );
	}
	
	void Pi_h_alpha(const baseFElement & K,KN_<double> & v) const 
	{ // compute the coef of interpolation ...
	    const Triangle & T(K.T);
	    int k=0;
	    for (int i=0;i<3;i++)
	      {  
		  R2 E(Ortho? T.Edge(i) : -T.Edge(i).perp());
		  
		  R s = T.EdgeOrientation(i) ;
		  for (int p=0;p<QFE.n;++p) 
		    {
		      R l0 = QFE[p].x, l1 = 1-QFE[p].x;
		      R p0= (2*l0-l1)*2;// poly othogonaux to \lambda_1
		      R p1= (2*l1-l0)*2;// poly othogonaux to \lambda_0
		      R cc1 = s*p0*QFE[p].a; // 
		      R cc0 = s*p1*QFE[p].a; //
		      if(s<0) Exchange(cc1,cc0); // exch lambda0,lambda1
		      v[k++]= cc0*E.x;
		      v[k++]= cc0*E.y; 
		      v[k++]= cc1*E.x;
		      v[k++]= cc1*E.y; 
		    }
	      }
	    R sy= Ortho ? -1 : 1; 
	    
	    for (int p=0;p<QFK.n;++p) 	
	      {
		v[k++]=sy*QFK[p].a * T.area; 
		v[k++]=QFK[p].a * T.area; 
	      }
	    // cout << " k= " << k << " == " << this->pij_alpha.N() << endl;
	    assert(k==this->pij_alpha.N());
	}
	void FB(const bool * whatd, const Mesh & Th,const Triangle & K,const R2 &P, RNMK_ & val) const;
    } ;
    // ENDOFCLASS TypeOfFE_PkEdge    
    
    void TypeOfFE_RT1_2d::FB(const bool * whatd,const Mesh & ,const Triangle & K,const R2 & Phat,RNMK_ & val) const
    {
      R2 X=K(Phat);
      R2 Q[]={ R2(K[0]), R2(K[1]),R2(K[2])};
      R2 A[]={ R2(Q[1],Q[2]), R2(Q[2],Q[0]), R2(Q[0],Q[1])}; 
      R2 B[]={ R2(A[1],A[2]), R2(A[2],A[0]), R2(A[0],A[1])}; 
      R l0=1-Phat.x-Phat.y,l1=Phat.x,l2=Phat.y; 
      R L[3]={l0,l1,l2};
      
      static int count=10;
      
      if( count < 0)
	{
	  cout << "TypeOfFE_RT1_2d "<< " " << A[0]+A[1]+A[2] << " " <<  B[0]+B[1]+B[2] << endl;
	  cout << det(Q[0],Q[1],Q[2]) << " X = " << X << " Phat ="  << Phat << endl;
	  cout<< "Q="  << Q[0]<< "," << Q[1] << " , " << Q[2] <<endl; 
	  cout<< "A="  << A[0]<< "," << A[1] << " , " << A[2] << endl; 
	  cout<< "B="  << B[0]<< "," << B[1] << " , " << B[2] <<endl; 
	}
      /*
       THE 2 DOF k=0,1  are: on edge e   f -> \int_e f \lambda_{e+k} . n_e 
       THE 2 internal dof are : f -> \int_K f e_i  where e_i is the canonical basis of R^2 
       
       
       so the basis function are 
       
       let call \phi_i the basic fonction of RT0 (without orientation) so the normal is exterior.
       \phi_i (X) = ( X- Q_i ) / (2 |K|) =  \lambda_{i+1} Curl( \lambda_{i+2}) - \lambda_{i+2} Curl( \lambda_{i+1})
       
       edge function j=0,1
       i1= i+j+1, i2= i+2-j  remark : {i,i1,i2} <=> {i,i+1,i+2}  
       \phi_i ( \lambda_{i1} - 4/3 \lambda_i) + 1/3 \phi_{i1}\lambda_{i1}
       
       internal function are 
       \sum   bx_i \phi_{i}\lambda_{i}
       \sum   by_i \phi_{i}\lambda_{i}
       \sum bx_i = 1/c0
       \sum by_i = 1/c0 
       we have 
       \phi_{i} = A_{i+2}  \lambda_{i+1} - A_{i+1}  \lambda_{i+2}
       with
       A_i = Th.edge(i)/ ( 2 |K])    
       B_i = A_{i+2} - A_{i+1}  
       det( B_i ) = 9 *2 area 
       to be homogene
       cc0= |K]  sqrt(area)*sqrt(18)
       ccK= 9 *2 area *c0; 
       bx_0 = det(R2(cc0,0),B1,B2)/ ( cck)  
       
       
       so all basic d function are the sum of 3 function 
       
       sum_{k=0}^2  c_k  phi_{p_k} lambda_{l_k} 
       
       */
      
      
      assert( val.N()>=ndf);
      assert(val.M()==2);
      int ee=0;
      
      val=0; 
      
      R2 phi[3] = { X-Q[0], X-Q[1], X-Q[2] };// phi * area *2 
      
      int pI[8][3];// store p_k 
      int lI[8][3];// store l_k 
      R   cI[8][3];// store c_k 
      
      int df=0;
      R CKK = 2* K.area; 
      for(int e=0;e<3;++e)
	{
	  int i=e;
	  int ii[2]={(e+1)%3,(e+2)%3};
	  int i2=(e+2)%3;
	  R s = K.EdgeOrientation(e)/CKK;
	  if(s<0) Exchange(ii[0],ii[1]); // 
	  for(int k=0;k<2;++k,df++)
	    {
	      pI[df][0]= i;
	      lI[df][0]= ii[k];
	      cI[df][0]= s;
	      
	      pI[df][1]= i;
	      lI[df][1]= i;
	      cI[df][1]= -s*4./3.;
	      
	      pI[df][2]= ii[k];
	      lI[df][2]= ii[k];
	      cI[df][2]= s/3.;
	      
	      
	    }
	}
      /*     
       if(count<0)
       {
       // verif. 
       R2 PP[] ={ R2(0.5,0),R2(0.5,0.5),R2(0,0.5)};
       int err=0;
       for(int df = 0;df < 6;++df)
       {  
       cout << " df = " << df << " : \t";
       for(int k=0;k<3;++k)
       cout  <<"+ " << cI[df][k] << " *l" << lI[df][k]  << " *phi" << pI[df][k] << "  \t  ";
       cout << endl;    
       
       R2 fd;
       for(int p=0;p<3;++p)
       {
       R L[3]; PP[p].toBary(L);
       R2 X=K(PP[p]);
       cout << X << " ,\t " << L[0] << " " << L[1] << " " <<L[2] ;
       R2 phi[3] = { X-Q[0], X-Q[1], X-Q[2] };// phi * area *2 
       
       R2 ff = (cI[df][0] * L[lI[df][0]]) * phi[pI[df][0]]
       + (cI[df][1] * L[lI[df][1]]) * phi[pI[df][1]]
       + (cI[df][2] * L[lI[df][2]]) * phi[pI[df][2]]
       ;
       fd += ff;
       cout << " :::: " << 3*ff 
       << " :  " <<  3*(cI[df][0] * L[lI[df][0]]) * phi[pI[df][0]]
       << " ; " <<  3*(cI[df][1] * L[lI[df][1]]) * phi[pI[df][1]]
       << " ; " <<  3*(cI[df][2] * L[lI[df][2]]) * phi[pI[df][2]]
       << " :  " <<  3*(cI[df][0] * L[lI[df][0]])  <<"(" <<  phi[pI[df][0]] <<")"
       << " ; " <<  3*(cI[df][1] * L[lI[df][1]]) <<"(" << phi[pI[df][1]]<<")"
       << " ; " <<  3*(cI[df][2] * L[lI[df][2]]) <<"(" << phi[pI[df][2]]<<")"		    
       <<endl;
       
       }
       if( fd.norme() > 1e-5) err++;
       cout << " Verif " << df << " [ " << 3*fd << " ] " << fd.norme()  <<endl;		
       }
       ffassert(err==0);
       
       }
       */
      R cK = 18.* K.area;
      R c0 = sqrt(cK);
      R cb = 12/c0;
      R ccK = K.area*cK/c0; 
      for(int k=0;k<2;++k,df++)
	{
	  
	  R2 PB(0,0);
	  PB[k]=cb;
	  // if( count <5) cout << " PB = " << PB << " df = " << df << " " << cK << " ==" << det(B[0] ,B[1],B[2]) <<endl;
	  R b0 = det(PB   ,B[1],B[2])/ ccK;
	  R b1 = det(B[0],PB   ,B[2])/ ccK;
	  R b2 = det(B[0],B[1],PB   )/ ccK;
	  
	  // if( count <5) cout << " S= "<< b0*B[0]+b1*B[1]+b2*B[2] << " s= " << (b0+b1+b2) << " b=" << b0 << " " << b1 << " " << b2 <<  endl;
	  pI[df][0]= 0;
	  lI[df][0]= 0;
	  cI[df][0]= b0;
	  
	  pI[df][1]= 1;
	  lI[df][1]= 1;
	  cI[df][1]= b1;
	  
	  pI[df][2]= 2;
	  lI[df][2]= 2;
	  cI[df][2]= b2;
	  
	}
      /*
       if( count< 5)
       {
       cout << Phat << " " << X << endl; 
       for( int e=0;e<3;++e)
       {  int e1= (e+1)%3, e2=(e+2)%3 ;
       cout << " phi e " << phi[e] << " == " << L[e1]*A[e2] -  L[e2]*A[e1] << endl;		
       }
       for(int df=0;df< 8;++df)
       {
       cout << " df = " << df << " : \t";
       for(int k=0;k<3;++k)
       cout  <<"+ " << cI[df][k] << " *l" << lI[df][k]  << " *phi" << pI[df][k] << "  \t  ";
       cout << endl;    
       }
       
       }
       */      
      int ortho0=0,ortho1=1; R s1ortho=1;
      if(Ortho) { ortho0=1; ortho1=0; s1ortho=-1;}
      
      if (whatd[op_id])
	{
          
          for(int df=0;df< 8;++df)
	    {
	      R2 fd(0.,0.) ;
	      for(int k=0;k<3;++k) {
		  fd += (cI[df][k] * L[lI[df][k]]) * phi[pI[df][k]] ;
	      }
	      
	      val(df,ortho0,op_id)= fd.x;
	      val(df,ortho1,op_id)= s1ortho*fd.y;
	    }
	}
      
      
      if(  whatd[op_dx] || whatd[op_dy] || whatd[op_dxx] || whatd[op_dyy] ||  whatd[op_dxy])
	{
          R2 DL[3]={K.H(0),K.H(1),K.H(2)};
	  R2 Dphix(1,0);
	  R2 Dphiy(0,1);
	  
	  if (whatd[op_dx])
	    {
	      
	      for(int df=0;df< 8;++df)
		{
		  R2 fd(0.,0.);
		  for(int k=0;k<3;++k)
		      fd += cI[df][k] * (DL[lI[df][k]].x * phi[pI[df][k]] + L[lI[df][k]]* Dphix);
		  val(df,ortho0,op_dx)= fd.x;
		  val(df,ortho1,op_dx)= s1ortho*fd.y;	      
		}
	      
	    }
	  if (whatd[op_dy])
	    {
	      
	      for(int df=0;df< 8;++df)
		{
		  R2 fd(0.,0.);
		  for(int k=0;k<3;++k)
		      fd += cI[df][k] * (DL[lI[df][k]].y * phi[pI[df][k]] + L[lI[df][k]]* Dphiy);
		  val(df,ortho0,op_dy)= fd.x;
		  val(df,ortho1,op_dy)= s1ortho*fd.y;	      
		}
	      
	    }
	  if(whatd[op_dxx] || whatd[op_dyy] ||  whatd[op_dxy])
	    {
	      cout << " to do FH RT1 dxx, dyy dxy " << endl; 
	      ffassert(0); 
	    }
	  
	}
      count++;
    }
    
    
    class TypeOfFE_BDM1_2d : public  TypeOfFE { 
    public:  
	static int Data[];
	static double Pi_h_coef[];
	bool Ortho;
	const QuadratureFormular1d & QFE;
	TypeOfFE_BDM1_2d(bool ortho)
	: 
	TypeOfFE(6,2,Data,1,1,
		 2*2*3*2,// nb coef mat interpole
		 3*2, // nb P interpolation 
		 0),
	QFE(QF_GaussLegendre2),
	Ortho(ortho)
	{  
	    //      cout << " Pk = " << k << endl;
	    int kkk=0,i=0;
	    for (int e=0;e<3;++e) 
		for (int p=0;p<QFE.n;++p) 
		  {
		    R2 A(TriangleHat[VerticesOfTriangularEdge[e][0]]);
		    R2 B(TriangleHat[VerticesOfTriangularEdge[e][1]]);
		    
		    pij_alpha[kkk++]= IPJ(2*e,i,0);
		    pij_alpha[kkk++]= IPJ(2*e,i,1);
		    pij_alpha[kkk++]= IPJ(2*e+1,i,0);
		    pij_alpha[kkk++]= IPJ(2*e+1,i,1);
		    
		    
		    P_Pi_h[i++]= B*(QFE[p].x)+ A*(1.-QFE[p].x);// X=0 => A  X=1 => B;       
		  }
	    //cout << kkk << " kkk == " << this->pij_alpha.N() << endl;
	    //cout << i << "  ii == " << this->P_Pi_h.N() << endl;
	    ffassert(kkk==this->pij_alpha.N());
	    ffassert(i==this->P_Pi_h.N() );
	}
	
	void Pi_h_alpha(const baseFElement & K,KN_<double> & v) const 
	{ // compute the coef of interpolation ...
	    const Triangle & T(K.T);
	    int k=0;
	    for (int i=0;i<3;i++)
	      {  
		  R2 E(Ortho? T.Edge(i) : -T.Edge(i).perp());
		  
		  R s = T.EdgeOrientation(i) ;
		  for (int p=0;p<QFE.n;++p) 
		    {
		      R l0 = QFE[p].x, l1 = 1-QFE[p].x;
		      R p0= s; // poly othogonaux to \lambda_1
		      R p1= -3*(l0-l1);// poly othogonaux to \lambda_0
		      R cc0 = p0*QFE[p].a; // 
		      R cc1 = p1*QFE[p].a; //
		      //if(s<0) Exchange(cc1,cc0); // exch lambda0,lambda1
		      v[k++]= cc0*E.x;
		      v[k++]= cc0*E.y; 
		      v[k++]= cc1*E.x;
		      v[k++]= cc1*E.y; 
		    }
	      }
	    // cout << " k= " << k << " == " << this->pij_alpha.N() << endl;
	    assert(k==this->pij_alpha.N());
	}
	void FB(const bool * whatd, const Mesh & Th,const Triangle & K,const R2 &P, RNMK_ & val) const;
    } ;
    // ENDOFCLASS TypeOfFE_PkEdge
    
  int  TypeOfFE_BDM1_2d::Data[]={
	3,3, 4,4, 5,5,  //  support on what 
	0,1, 0,1, 1,0,  // df on node 
	0,0, 1,1, 2,2,  // th node of df 
	0,0,0, 0,0,0, //  df previou FE
	0,1,2,3,4,5,  //  which df on prev 
	0,0,
	0,0,
	6,6 
    };
    
    
    void TypeOfFE_BDM1_2d::FB(const bool * whatd,const Mesh & ,const Triangle & K,const R2 & Phat,RNMK_ & val) const
    {
      R2 X=K(Phat);
      R2 Q[]={ R2(K[0]), R2(K[1]),R2(K[2])};
      R l0=1-Phat.x-Phat.y,l1=Phat.x,l2=Phat.y; 
      R L[3]={l0,l1,l2};
      R2 Dl[3]={K.H(0),K.H(1),K.H(2)};
      
      static int count=10;
      
      
      
      assert( val.N()>=6);
      assert(val.M()==2);
      int ee=0;
      
      val=0; 
      R cK = 2* K.area;
      
      int ortho0=0,ortho1=1; R s1ortho=1;
      if(Ortho) { ortho0=1; ortho1=0; s1ortho=-1;}
      
      if (whatd[op_id])
	{
          
          for(int df=0,e =0;e< 3;++e)
	    {
	      int e1=(e+1)%3, e2=(e+2)%3;
	      R s = K.EdgeOrientation(e) ;
	      R2 f1= (X-Q[e]) * s / cK ;
	      R2 f2= -(Dl[e1]*L[e2]+Dl[e2]*L[e1]).perp() ;
	      
	      
	      val(df,ortho0,op_id)= f1.x;
	      val(df++,ortho1,op_id)= s1ortho*f1.y;
	      
	      val(df,ortho0,op_id)= f2.x;
	      val(df++,ortho1,op_id)= s1ortho*f2.y;
	    }
	}
      
      
      if(  whatd[op_dx] || whatd[op_dy] || whatd[op_dxx] || whatd[op_dyy] ||  whatd[op_dxy])
	{
 	  R2 Dphix(1,0);
	  R2 Dphiy(0,1);
	  
	  if (whatd[op_dx])
	    {
	      
	      for(int df=0,e =0;e< 3;++e)
		{
		  int e1=(e+1)%3, e2=(e+2)%3;
		  R s = K.EdgeOrientation(e) ;
		  R2 f1=  R2( s / cK,0.) ;
		  R2 f2= -(Dl[e1]*Dl[e2].x+Dl[e2]*Dl[e1].x).perp() ;
		  
		  
		  val(df,ortho0,op_dx)= f1.x;
		  val(df++,ortho1,op_dx)= s1ortho*f1.y;
		  
		  val(df,ortho0,op_dx)= f2.x;
		  val(df++,ortho1,op_dx)= s1ortho*f2.y;
		}
	      
	    }
	  if (whatd[op_dy])
	    {
	      
	      for(int df=0,e =0;e< 3;++e)
		{
		  int e1=(e+1)%3, e2=(e+2)%3;
		  R s = K.EdgeOrientation(e) ;
		  R2 f1=  R2(0.,  s / cK) ;
		  R2 f2= -(Dl[e1]*Dl[e2].y+Dl[e2]*Dl[e1].y).perp() ;
		  
		  
		  val(df,ortho0,op_dy)= f1.x;
		  val(df++,ortho1,op_dy)= s1ortho*f1.y;
		  
		  val(df,ortho0,op_dy)= f2.x;
		  val(df++,ortho1,op_dy)= s1ortho*f2.y;
		}
	      
	    }
	  
	}
      count++;
    }
    
    
    
    
    
    // a static variable to add the finite element to freefem++
    static TypeOfFE_RT1_2d Elm_TypeOfFE_RT1_2d(false);// RT1    
    static TypeOfFE_RT1_2d Elm_TypeOfFE_RT1_2dOrtho(true);// RT1ortho  
    static TypeOfFE_BDM1_2d Elm_TypeOfFE_BDM1_2d(false);// BDM1    
    static TypeOfFE_BDM1_2d Elm_TypeOfFE_BDM1_2dOrtho(true);// BDM1ortho  
    static TypeOfFE_TD_NNS0 Elm_TD_NNS;
    static TypeOfFE_TD_NNS1 Elm_TD_NNS1;
    static AddNewFE FE__TD_NNS("TDNNS0",&Elm_TD_NNS); 
    static AddNewFE FE__TD_NNS1("TDNNS1",&Elm_TD_NNS1);     
    static AddNewFE Elm__TypeOfFE_RT1_2d("RT1",&Elm_TypeOfFE_RT1_2d); 
    static AddNewFE Elm__TypeOfFE_RT1_2dOrtho("RT1Ortho",&Elm_TypeOfFE_RT1_2dOrtho); 
    static AddNewFE Elm__TypeOfFE_BDM1_2d("BDM1",&Elm_TypeOfFE_BDM1_2d); 
    static AddNewFE Elm__TypeOfFE_BDM1_2dOrtho("BDM1Ortho",&Elm_TypeOfFE_BDM1_2dOrtho); 
    
} // FEM2d namespace 
// --- fin -- 


