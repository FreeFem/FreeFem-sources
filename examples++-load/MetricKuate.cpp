//  Implementation of P1-P0 FVM-FEM
// ---------------------------------------------------------------------
// $Id$
// compile and link with ./load.link  mat\_dervieux.cpp 
#include  <iostream>
#include  <cfloat>
#include  <cmath>
using namespace std;
/*
#include "error.hpp"
#include "AFunction.hpp"
#include "rgraph.hpp"
#include "RNM.hpp"
// remove problem of include 
#undef  HAVE_LIBUMFPACK
#undef HAVE_CADNA
#include "MatriceCreuse_tpl.hpp"
#include "MeshPoint.hpp"
#include "lgfem.hpp"
#include "lgmesh3.hpp"
#include "problem.hpp"
//#include "ellipsemax.hpp"
*/
#include "ff++.hpp" 

using namespace std;


R Min (const R a,const R b);
R Max (const R a,const R b);
inline void Exchange (R a,R b) {R c=a;a=b;b=c;};
//const R precision = 1e-10;

int  LireTaille( const char * NomDuFichier,int &nbnoeuds);
int  Lire( const char * NomDuFichier,int n ,R2 noeuds[] );
template<typename T>
bool from_string( const string & Str, T & Dest )
{
  // créer un flux à partir de la chaine donnée
  istringstream iss( Str );
  // tenter la conversion vers Dest
  return iss >> Dest != 0;
};


void metrique(int nbpoints, R2 * Point,R &A, R &B,R &C,R epsilon)
{


  C=0.;
   
  R  epsilon0=1e-5,precision=1e-15,delta=1e-10;
  R inf =1e50;
 
  R Rmin=inf,Rmax=0.;
  
  int indiceX0=0;
  
 
  R2 PPoint[nbpoints];
  
  for(int i=0;i<nbpoints;i++)
    {
     
      Rmax=Max(Rmax,Point[i].norme());
      
      //---déplacement des points situées sur les axes--------------
      if(abs(Point[i].x)<=precision)
	{
	  if(Point[i].y<0)
	    {
	      Point[i].x=-delta;
	      Point[i].y=-sqrt(pow(Point[i].y,2)-pow(Point[i].x,2));
	    }
	  if(Point[i].y>0)
	    {
	      Point[i].x=delta;
	      Point[i].y=sqrt(pow(Point[i].y,2)-pow(Point[i].x,2));
	    }
	}
      
      if(abs(Point[i].y)<=precision)
	{
	  if(Point[i].x<0)
	    {
	      Point[i].y=-delta;
	      Point[i].x=-sqrt(pow(Point[i].x,2)-pow(Point[i].y,2));
	    }
	  if(Point[i].x>0)
	    {
	      Point[i].y=delta;
	      Point[i].x=sqrt(pow(Point[i].x,2)-pow(Point[i].y,2));
	    }
	}
      //-----------------------------------------------------------
      assert(abs(Point[i].x*Point[i].y)>=pow(precision,2));

      if(Rmin>Point[i].norme())
	{
	  indiceX0=i;
	  Rmin=Point[i].norme();
	}
    
      //cout<<Point[i]<<endl;	
    }
  
  
  
  //-------permutation des indices de la liste des points : 
  //ranger la liste en commençant par le point X0-------
  for(int k=0;k<nbpoints-indiceX0;k++)
    {
      PPoint[k]=Point[k+indiceX0];
    }
  
  for(int k=nbpoints-indiceX0;k<nbpoints;k++)
    {
      PPoint[k]=Point[k-nbpoints+indiceX0];
    }
  
  
  for(int i=0;i<nbpoints;i++) Point[i]=PPoint[i];	
 
  
  
  
  
  //----------------------------------------------------------------
  
 
  int test=-1;
  
  
  
 
  
  R X0;
  R Y0;
  

  R bmin= 0.,bmax=inf,b1,b2,aik=0.,bik=0.,cik=0.;
 R Xk=0.,Yk=0.,Ck=0.,Ak =0.,Bk=0.,Xi=0.,Yi=0.,ri,detXY=0.,Ri,R0,r0;

 
 X0=Point[0].x;
 Y0=Point[0].y;
 r0=Point[0].norme();
 assert(r0 == Rmin);

 //cout<<" Rmin = "<<Rmin<<" Rmax =  "<<Rmax<<endl;
 

 R EPS=0.;// pour recuperer la valeur de epsilon0 optimale
  R epsilonmax=r0*(1.-r0/Rmax)/20.;
 
  R * Tabepsilon;

  int neps =4;
 //--------- discertisation de epsilon0----------------------------------
 if(epsilonmax>1e-2)
   {
     neps=10;
     Tabepsilon=new R[neps];
     Tabepsilon[0]=1e-5;
     Tabepsilon[1]=1e-4;
     Tabepsilon[2]=1e-3;
     for(int i=3;i<neps;i++)
       {
	 Tabepsilon[i]=(i-3)*(epsilonmax-1e-2)/(neps-4.) +1e-2;
  
       }
   }


 else
   {
    
     Tabepsilon=new R[neps];
     Tabepsilon[0]=1e-5;
     Tabepsilon[1]=1e-4;
     Tabepsilon[2]=1e-3;
     Tabepsilon[3]=1e-2;
   }
 //------------------------------------------------------------------------

 int condition = -1;

 if(r0<=epsilon0)epsilon0=r0*epsilon0;
 B=A=1./((r0-epsilon0)*(r0-epsilon0));
  R epsilon0min=epsilon0;

   if(abs(Rmin-Rmax)>1e-5)
    {
     //cout<<" Rmax - Rmin  "<<Rmax-Rmin<<endl;
     for(int ee=0;ee<neps-1;ee++) //boucle sur epsilon0---------------
       {
	 epsilon0= Tabepsilon[ee];
	 if(r0<=epsilon0)epsilon0=r0*epsilon0;
	 assert(r0>epsilon0);
	 R0=r0/(r0-epsilon0);
	 
	 for(int i=1;i<nbpoints;i++)	 //boucle sur chaque noeud
	   {
	     
	     Xi=Point[i].x;
	     Yi=Point[i].y;
	     ri=Point[i].norme();
	     if(ri<=epsilon)epsilon=ri*epsilon;
	     assert(ri>epsilon);
	     Ri = ri/(ri-epsilon);
	     detXY=Xi*Y0-Yi*X0;
	     
	     
	     
	     
	     
	     
	     //------deplacement des points alignés avec l'origine et X0-----------
	     if(abs(detXY)<=precision)
	       {
		 
		 Xi += delta;
		 
		 if(Yi<0)Yi= -sqrt(pow(ri,2)- pow(Xi,2));
		 else Yi=sqrt(pow(ri,2)- pow(Xi,2));
		 Point[i].x=Xi;
		 Point[i].y=Yi;
		 ri=Point[i].norme();
		 if(ri<=epsilon)epsilon=ri*epsilon;
		 assert(ri>epsilon);
		 Ri = ri/(ri-epsilon);
		 
	       }
	     
	     detXY=Xi*Y0-Yi*X0;
	     
	     assert(abs(detXY)>=precision);   
	     
	     
	     //-----------------------------------------------------------------
	     
	     //-----racines du polynome en b à minimiser----------------------------
	    
	     R bb1=(1./pow(detXY,2))*(pow(X0*Ri,2)+pow(Xi*R0,2)-2.*abs(Xi*X0)*sqrt(pow(R0*Ri,2)-pow(detXY/(Rmax*(r0-epsilon0)),2)));
	    
	     R bb2=(1./pow(detXY,2))*(pow(X0*Ri,2)+pow(Xi*R0,2)+2.*abs(Xi*X0)*sqrt(pow(R0*Ri,2)-pow(detXY/(Rmax*(r0-epsilon0)),2))); 
	   
	    
	     //--fin----racines du polynome en b à minimiser--------------------
	     
	     
	     bmax=Min(bb2,pow(Rmax/pow((r0),2),2));
	     
	     
	     bmin=Max(1./(Rmax*Rmax),bb1);//minoration de b
	       R Cte = Max(1e-9,(bmax-bmin)*1e-9);
	       bmin=bmin*(1.+Cte);
	       bmax=bmax*(1.-Cte);
	     
	   
	     //bornes de b-----------------------------------------------------------
	     
	     //cas:  majoration de c --------------------------------------------
	     R Li=X0*Xi*(pow(Rmax/pow(r0-epsilon0min,2),2)-1./pow(Rmax,2))+(pow(Ri*X0,2)-pow(R0*Xi,2))/detXY;
	     R LiXY=Xi*Y0+Yi*X0;
	     
	     if(abs(LiXY)>=precision)
	       {
		 condition=1;
		 
		 if(Xi*X0>0)
		   {
		     if(LiXY>0) bmin=Max(bmin,-Li/LiXY);
		     else bmax=Min(bmax,-Li/LiXY);
		   }
		 else
		   {
		     
		     if(LiXY<0) bmin=Max(bmin,-Li/LiXY);
		     else bmax=Min(bmax,-Li/LiXY); 
		   }
		 
	       }
	     
	     else
	       {
		 if(Li<0) condition =0;
		 else condition =1;
	       }
	     
	     
	     
	     //cas  minoration de c --------------------------------------------
	     Li=X0*Xi*(-pow(Rmax/pow(r0-epsilon0min,2),2)+1./pow(Rmax,2))+(pow(Ri*X0,2)-pow(R0*Xi,2))/detXY;
	     LiXY=Xi*Y0+Yi*X0;
	     
	     if(abs(LiXY)>=precision)
	       {
		 condition=1;
		 
		 if(Xi*X0>0)
		   {
		     if(LiXY<0) bmin=Max(bmin,-Li/LiXY);
		     else bmax=Min(bmax,-Li/LiXY);
		   }
		 else
		   {
		     
		     if(LiXY>0) bmin=Max(bmin,-Li/LiXY);
		     else bmax=Min(bmax,-Li/LiXY); 
		   }
		 
	       }
	     
	     else
	       {
		 if(Li>0) condition =0;
		 else condition =1;
	       }
	     
	     if (condition==1)
	       {
		 
		 
		 //--cas : minoration de a-----------------------------------------------
		 
		 R Gi=((Xi*Yi*R0*R0-X0*Y0*Ri*Ri)/detXY +Xi*X0/(Rmax*Rmax))/(Yi*Y0);
		 
		 if(Xi*X0>0)
		   {
		     if(Yi*Y0>0) bmin=Max(bmin,Gi);
		     else bmax=Min(bmax,Gi);
		   }
		 else
		   {
		     if(Yi*Y0<0) bmin=Max(bmin,Gi);
		     else bmax=Min(bmax,Gi);
		   }
		 
		 //cas :majoration de a------------------------------------------------
		 
		 
		 R Hi=(Xi*X0*Rmax*Rmax/pow((r0-epsilon0min),4)+(Xi*Yi*R0*R0-X0*Y0*Ri*Ri)/detXY)/(Yi*Y0);
		 if(Xi*X0>0)
		   {
		     if(Yi*Y0>0) bmax=Min(bmax,Hi);
		     else bmin=Max(bmin,Hi);
		   }
		 else
		   {
		     if(Yi*Y0<0) bmax=Min(bmax,Hi);
		     else bmin=Max(bmin,Hi);
		     
		   }
		 
		 
		 
		 
		 
		 //------fin bornes de b------------------------------------------------
		 b2=bmax;
		 b1=bmin;
		 
		 for(int k=1;k<nbpoints ;k++)//on balaye les contraintes
		   {
		   
		     Xk=Point[k].x;
		     Yk=Point[k].y;
		     Bk=(Yk*Yk*Xi*X0 +Xk*(Xk*Yi*Y0-Yk*(Yi*X0+Xi*Y0)))/(Xi*X0);
		     Ck=(X0*Xi*detXY-Xk*(Xi*R0*R0*(Yk*Xi-Yi*Xk) +X0*Ri*Ri*(-Yk*X0+Y0*Xk)))/(Xi*X0*detXY);
		     
		     assert(abs(Xi*X0*Y0*Yi*Xk*Yk)>=pow(precision,5));
		     if(abs(Bk)>precision)//non nul
		       {
			 if(Bk<=0) bmax=Min(bmax,Ck/Bk);
			 
			 else  bmin=Max(bmin,Ck/Bk);
			 
			 if((bmax<b1)||(bmin>b2)||(bmin>bmax))
			   {
			     //cout<<" i = "<<i<<"  k = "<<k<<endl;	
			     test=0;
			     break;  
			   }
			 
			 else  test=1;		       
			 
		       }
		     else
		       {
			 
			 if (Ck>precision) 
			   { 
			     test=0;
			     break;  	
			   }
			 else //Ck<=0
			   {
			     test=-1;// 1 peut etre
			     
			   }
		       }
		     
		     
		   }
		 if(test==1)
		   { 
		   
		     R a0=-pow((detXY/(Xi*X0)),2);
		     R a1=2.*(pow(Ri/Xi,2)+pow(R0/X0,2));
		     if(((a0*bmax+a1)*bmax) <((a0*bmin+a1)*bmin)) bik=bmax;
		     else bik=bmin;
		     aik=(Ri*Ri*Y0*X0 -R0*R0*Yi*Xi+bik*Yi*Y0*detXY)/(detXY*Xi*X0);
		     //(Ri*Ri*Y0/Xi - R0*R0*Yi/X0)/detXY+bik*Yi*Y0/(Xi*X0);
		     cik=( -Ri*Ri*X0*X0 + R0*R0*Xi*Xi-bik*(Yi*X0+Y0*Xi)*detXY)/(detXY*Xi*X0);
		    
		     
		     
		   
		     assert((4.*aik*bik-cik*cik)>=0.);// aire positive	 
		     assert(abs((4.*aik*bik-cik*cik)-pow(2./(Rmax*(r0-epsilon0)),2))>0);// aire positive
		     if((4.*aik*bik-cik*cik)<=(4.*A*B-C*C)) 
		       {
			 A=aik;
			 B=bik;
			 C=cik;
			 EPS=epsilon0;
		       }
		     
		     
		     
		   }
		 
		 
		 
		 //-----------------------------------------------------------------  
		 
	       }  
	     
	     
	   }
       }
     
     }
 else
   {
     A=B=1./(Rmin*Rmin);
     C=0.;
     }
 
}














int  LireTaille( const char * NomDuFichier,int & nbnoeuds)
{ //Lire le maillage  sur le fichier de nom NomDuFichier
  
  //Ouverture du fichier  a partir de son nom
  ifstream f( NomDuFichier );
  //char  buffer[BUFSIZ];
  string buffer;
  
  nbnoeuds =0;
  if( !f )
    { cerr << "Erreur a l'ouverture du fichier " << NomDuFichier << endl;
      return 1; }
  
  while ( getline( f,buffer,'\n' ) )
    {
      if((buffer[0]!='#')&&(buffer!="")){ nbnoeuds +=1;
	// cout<<buffer<<endl;
      }
    }
  
  
  return 0;
}


int  Lire( const char * NomDuFichier,int n ,R2 noeuds[] ){
  
  ifstream f( NomDuFichier );
  
  string buffer;
  int i=0 ;
  
  
  while(i<n)
    {   
      f >> buffer;			
      if(buffer[0]=='#')
	{  getline( f,buffer );}
      else{
	
	
	
	// Lecture X Y Z de chacun des noeuds
	
	from_string( buffer,noeuds[i++].x ); 
	f>>noeuds[i-1].y>>buffer;
	
      }
      
      
      
    }
  
  
  
  return 0;
}





R Min (const R a,const R b){return a < b ? a : b;}
R Max (const R a,const R b){return a > b ? a : b;}







//  metrixkuate(Th,np,o,err,[m11,m12,m22]);
class MetricKuate :  public E_F0mps 
{
public:
  typedef bool  Result;
  Expression expTh;
  Expression expnp;
  Expression exphmin;
  Expression exphmax;
  Expression experr;
  Expression m11,m12,m22;
  Expression px,py;
  
  MetricKuate(const basicAC_F0 & args)
  {

    args.SetNameParam();
    expTh= to<pmesh>(args[0]);  // a the expression to get the mesh
    expnp= to<long>(args[1]);  // a the expression to get the mesh
    exphmin= to<double>(args[2]);  // a the expression to get the mesh
    exphmax= to<double>(args[3]);  // a the expression to get the mesh
    experr= to<double>(args[4]);  // a the expression to get the mesh
    //  a array expression [ a, b]
    const E_Array * ma= dynamic_cast<const E_Array*>((Expression) args[5]);
    const E_Array * mp= dynamic_cast<const E_Array*>((Expression) args[6]);
    if (ma->size() != 3) CompileError("syntax: MetricKuate(Th,np,o,err,[m11,m12,m22],[xx,yy])");
    if (mp->size() != 2) CompileError("syntax: MetricKuate(Th,np,o,err,[m11,m12,m22],[xx,yy])");
    int err =0;
    m11= CastTo<KN<double> * >((*ma)[0]); // fist exp of the array (must be a  double)
    m12= CastTo<KN<double> * >((*ma)[1]); // second exp of the array (must be a  double)
    m22= CastTo<KN<double> * >((*ma)[2]); // second exp of the array (must be a  double)
    px= CastTo<double * >((*mp)[0]); // fist exp of the array (must be a  double)
    py= CastTo<double * >((*mp)[1]); // second exp of the array (must be a  double)
  }

  ~MetricKuate()
  {
  }

  static ArrayOfaType  typeargs()
  { return  ArrayOfaType(
			 atype<pmesh>(),
			 atype<long>(),
			 atype<double>(),
			 atype<double>(),
			 atype<double>(),
			 atype<E_Array>(),
			 atype<E_Array>());
  }
  static  E_F0 * f(const basicAC_F0 & args){ return new MetricKuate(args);}
  AnyType operator()(Stack s) const ;

};

// the evaluation routine
AnyType MetricKuate::operator()(Stack stack) const
{
  MeshPoint *mp(MeshPointStack(stack)) , mps=*mp;
  Mesh * pTh = GetAny<pmesh>((*expTh)(stack));
  long np=  GetAny<long>((*expnp)(stack));
  double hmin=  GetAny<double>((*exphmin)(stack));
  double hmax=  GetAny<double>((*exphmax)(stack));
  KN<double> *pm11,*pm12,*pm22;
  double *pxx,*pyy;
  pm11=  GetAny<KN<double>*>((*m11)(stack));
  pm22=  GetAny<KN<double>*>((*m22)(stack));
  pm12=  GetAny<KN<double>*>((*m12)(stack));
  pxx =  GetAny<double*>((*px)(stack));
  pyy =  GetAny<double*>((*py)(stack));
  ffassert(pTh);
  KN<R2> Pt(np);
  Mesh & Th (*pTh);
  cout << " MetricKuate " << np << " hmin = "<< hmin << " hmax = " << hmax  << " nv = " << Th.nv <<  endl;
  R hmx2=1./(hmax*hmax);
  R hmn2=1./(hmin*hmin);

  ffassert(pm11->N()==Th.nv);
  ffassert(pm12->N()==Th.nv);
  ffassert(pm22->N()==Th.nv);
  {
    for (int iv=0;iv<Th.nv;iv++)
      {
	R2 P=Th(iv);
	MeshPointStack(stack)->set(P.x,P.y);
	double m11=1,m12=0,m22=1;
	for(int i=0;i<np;i++)
	  {
	    double t=(M_PI*2.*i+0.5)/np;
	    *pxx=cos(t);
	    *pyy=sin(t);
	    double ee =  fabs(GetAny<double>((*experr)(stack)));
	    *pxx*=M_E;
	    *pyy*=M_E;
	    double eee =  fabs(GetAny<double>((*experr)(stack)));
	    ee = max(ee,1e-30);
	    eee = max(eee,1e-30);
	    //  e^p  = eee/ee  
	    double p = Min(Max(log(eee)-log(ee),0.1),10);
	    // c^p ee = 1
	    // c = (1/ee)^1/p
	    double  c=pow(1./ee,1./p);
	    c = min(max(c,hmin),hmax);
	    Pt[i].x = *pxx*c/M_E;
	    Pt[i].y = *pyy*c/M_E;	   
	    if(iv==0) {
	      cout << Pt[i] << "  ++++ " << i <<" " << t <<  " " << p << " c = " << R2(*pxx*c/M_E,*pyy*c/M_E) << "e=  " << ee << " " << eee << " " << c << endl;
	      }
	     
	  }
	double epsilon=1e-5;
	metrique(np,Pt,m11,m22,m12,epsilon);
	if(iv==0)  cout << "  ---- 11,12,22 : " << m11 << " "<< m12/2.  << " "<< m22 << endl; 
	(*pm11)[iv]=m11;
	(*pm12)[iv]=m12/2.;;
	(*pm22)[iv]=m22;

      }
  }
  *mp = mps;
  return true;
}


class Init { public:
  Init();
};
LOADINIT(Init);
Init::Init()
{
  cout << "\n  -- lood: init MetricKuate\n";
  Global.Add("MetricKuate","(", new OneOperatorCode<MetricKuate >( ));
}
