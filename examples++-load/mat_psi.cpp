// Example C++ function "myfunction", dynamically loaded into "load.edp"
// ---------------------------------------------------------------------
// $Id$

#include  <iostream>
#include  <cfloat>
#include  <cmath>
using namespace std;
#include "error.hpp"
#include "AFunction.hpp"
using namespace std;  
#include "rgraph.hpp"
#include "RNM.hpp"
#undef  HAVE_LIBUMFPACK
#include "MatriceCreuse_tpl.hpp"
#include "MeshPoint.hpp"
#include "lgfem.hpp"
#include "lgmesh3.hpp"
#include "lgsolver.hpp"
#include "problem.hpp"

class MatrixUpWind0 :  public E_F0 { public: 
  typedef Matrice_Creuse<R> * Result;
  Expression emat,expTh,expc,expu1,expu2; 
  MatrixUpWind0(const basicAC_F0 & args) 
  {   
    
    args.SetNameParam();
    emat =args[0];
    expTh= to<pmesh>(args[1]);
    expc = CastTo<double>(args[2]);
    const E_Array * a= dynamic_cast<const E_Array*>((Expression) args[3]);
    if (a->size() != 2) CompileError("syntax:  MatrixUpWind0(Th,rhi,[u1,u2])");
    int err =0;
    expu1= CastTo<double>((*a)[0]);
    expu2= CastTo<double>((*a)[1]);
    
  }
  
  ~MatrixUpWind0() 
  {  
  }     
  
  static ArrayOfaType  typeargs() { return  ArrayOfaType(atype<Matrice_Creuse<R>*>(),atype<pmesh>(),atype<double>(),atype<E_Array>());}
  static  E_F0 * f(const basicAC_F0 & args){ return new MatrixUpWind0(args);} 
  AnyType operator()(Stack s) const ;
  
};

int   gladys(double q[3][2], double u[2],double c[3], double a[3][3] )  //PSI Deconninck
{				// computes matrix a on a triangle for the Chacon-Reina Petrof-Galerkin upwind
  
  // working arrays
  double dw[3][2]; // basis function gradients times  area
  double ua[2], kk[3], beta[3]; // to define a[][]
  double udc=0;  // u.grad(w)*area
  bool oneaval=false;
  int i1=-1;
  
  for(int i=0;i<3;i++)
    {
      int ip=(i+1)%3, ipp=(ip+1)%3;
      for(int j=0;j<2;j++)
	dw[i][1-j]= (2*j-1)*(q[ipp][j]-q[ip][j])/2;
    }
  
  for(int i=0;i<3;i++){
    kk[i] = u[0]*dw[i][0]+u[1]*dw[i][1] ;
    udc += kk[i]*c[i];
  }
  
  for(int i=0;i<3;i++)
    {
      ua[0]=u[0]; ua[1]=u[1];
      int ip=(i+1)%3, ipp=(ip+1)%3;
      if(kk[i]>0 && kk[ip]<=0 && kk[ipp]<=0)
	{
	  beta[i]=1; beta[ip]=0; beta[ipp]=0; oneaval=true;
	}
      else if(kk[i]<=0 && kk[ip]>0 && kk[ipp]>0) i1=i;
    }
  
  if(!oneaval)
    {
      if(i1<0)cout<<"bug\n";
      int i=i1, ip=(i+1)%3, ipp=(i+2)%3;
      double lambda = (c[ip]-c[i])*(c[ipp]-c[i]);
      if (fabs(lambda) < -1e-20)
        {
	  return 0;
	}
      if(lambda < 0)
	{
	  if (udc>0)
	    {
	      beta[i]=0; beta[ip]=0; beta[ipp]=1;
	      ua[0] = udc*(q[ipp][0]-q[i][0])/(c[ipp]-c[i]); 
	      ua[1] = udc*(q[ipp][1]-q[i][1])/(c[ipp]-c[i]);
	    }
	  else
	    {
	      beta[i]=0; beta[ipp]=0; beta[ip]=1;
	      ua[0] = udc*(q[ip][0]-q[i][0])/(c[ip]-c[i]); 
	      ua[1] = udc*(q[ip][1]-q[i][1])/(c[ip]-c[i]);
	    }
	}
      else 
	{
	  beta[i]=0; 
	  beta[ip]=kk[ip]*(c[ip]-c[i])/udc; 
	  beta[ipp]=kk[ipp]*(c[ipp]-c[i])/udc;
	}
    }
  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)
      a[i][j]= beta[i]*(ua[0]*dw[j][0]+ua[1]*dw[j][1]);
  return 1;
}		


AnyType MatrixUpWind0::operator()(Stack stack) const 
{
  Matrice_Creuse<R> * sparce_mat =GetAny<Matrice_Creuse<R>* >((*emat)(stack)); 
  MatriceMorse<R> * amorse =0; 
  MeshPoint *mp(MeshPointStack(stack)) , mps=*mp;
  Mesh * pTh = GetAny<pmesh>((*expTh)(stack));
  ffassert(pTh);
  Mesh & Th (*pTh);
  {
    map< pair<int,int>, R> Aij;
    KN<double> cc(Th.nv);
    double infini=DBL_MAX;   
    cc=infini;
    for (int it=0;it<Th.nt;it++)
      for (int iv=0;iv<3;iv++)
	{
	  int i=Th(it,iv);
	  if ( cc[i]==infini) { // if nuset the set 
	    mp->setP(&Th,it,iv);
	    cc[i]=GetAny<double>((*expc)(stack));
	  }
	}
    
    for (int k=0;k<Th.nt;k++)
      {   
	const Triangle & K(Th[k]); 
	const Vertex & A(K[0]), &B(K[1]),&C(K[2]);
	R2 Pt(1./3.,1./3.);
	R u[2];
	MeshPointStack(stack)->set(Th,K(Pt),Pt,K,K.lab);
	u[0] = GetAny< R>( (*expu1)(stack) ) ;
	u[1] = GetAny< R>( (*expu2)(stack) ) ;
	
	int ii[3] ={  Th(A), Th(B),Th(C)};
	double q[3][2]= { { A.x,A.y} ,{B.x,B.y},{C.x,C.y} } ;  // coordinates of 3 vertices (input)
	double c[3]={cc[ii[0]],cc[ii[1]],cc[ii[2]]};
	double a[3][3];
	if (gladys(q,u,c,a) )
	  {
	    for (int i=0;i<3;i++)
	      for (int j=0;j<3;j++)
		if (fabs(a[i][j]) >= 1e-30)
		  Aij[make_pair(ii[i],ii[j])]+=a[i][j];	   
	  }
      }
    amorse=  new MatriceMorse<R>(Th.nv,Th.nv,Aij,false); 
  }
  sparce_mat->Uh=UniqueffId();
  sparce_mat->Vh=UniqueffId();
  sparce_mat->A.master(amorse);
  sparce_mat->typemat=(amorse->n == amorse->m) ? TypeSolveMat(TypeSolveMat::GMRES) : TypeSolveMat(TypeSolveMat::NONESQUARE); //  none square matrice (morse)
  *mp=mps;
  
  if(verbosity>3) { cout << "  End Build MatrixUpWind : " << endl;}
  
  return sparce_mat;  
}





class Init { public:
  Init();
};
 Init init;
 Init::Init()
   {
     cout << " lood: init Mat Chacon " << endl;
     Global.Add("MatUpWind0","(", new OneOperatorCode<MatrixUpWind0 >( ));
   }
