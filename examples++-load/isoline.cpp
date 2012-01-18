// ORIG-DATE: Novembre 2008
// -*- Mode : c++ -*-
//
// SUMMARY  :  
// USAGE    : LGPL      
// ORG      : LJLL Universite Pierre et Marie Curie, Paris,  FRANCE 
// AUTHOR   : Jacques Morice, Modif par F. Hecht Dec. 2011
// E-MAIL   : jacques.morice@ann.jussieu.fr
//

/* 
 This file is part of Freefem++
 
 Freefem++ is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation; either version 2.1 of the License, or
 (at your option) any later version.
 
 Freefem++  is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.
 
 You should have received a copy of the GNU Lesser General Public License
 along with Freefem++; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

 Thank to the ARN ()  FF2A3 grant
 ref:ANR-07-CIS7-002-01 
 */

/*
calcul demander par F. Hecht
*/

#ifndef WITH_NO_INIT
#include "ff++.hpp"
#include "AFunction_ext.hpp"
#endif


using namespace std;

#include <set>
#include <vector>
#include <map> 
#include <algorithm>
//#include "msh3.hpp"

using namespace  Fem2D;


// fonction determinant les points d'intersection
  static int debug =0;
class ISOLINE_P1_Op : public E_F0mps 
{
public:
   
  Expression eTh,eff,emat,exx,eyy,exy,iso;
  static const int n_name_param = 7; // 
  static basicAC_F0::name_and_type name_param[] ;
  Expression nargs[n_name_param];

  double  arg(int i,Stack stack,double a) const{ return nargs[i] ? GetAny< double >( (*nargs[i])(stack) ): a;}
  long  arg(int i,Stack stack,long  a) const{ return nargs[i] ? GetAny< long  >( (*nargs[i])(stack) ): a;}
  
  KN<long> * arg(int i,Stack stack,KN<long> *  a) const{ return nargs[i] ? GetAny< KN<long> *  >( (*nargs[i])(stack) ): a;}
  string * arg(int i,Stack stack,string *  a) const{ return nargs[i] ? GetAny< string  *  >( (*nargs[i])(stack) ): a;}

 
 
public:
  ISOLINE_P1_Op(const basicAC_F0 &  args,Expression tth, Expression fff,Expression xxx,Expression yyy) 
    : eTh(tth),emat(0),eff(fff),exx(xxx),eyy(yyy),exy(0) 
  {
    args.SetNameParam(n_name_param,name_param,nargs);
  } 
  ISOLINE_P1_Op(const basicAC_F0 &  args,Expression tth, Expression fff,Expression xxyy) 
    : eTh(tth),emat(0),eff(fff),exx(0),eyy(0),exy(xxyy) 
  {
    args.SetNameParam(n_name_param,name_param,nargs);
  } 
    
  AnyType operator()(Stack stack)  const ;
};


basicAC_F0::name_and_type ISOLINE_P1_Op::name_param[]= {
  {  "iso", &typeid(double)},
  {  "close", &typeid(long)},
  {  "smoothing", &typeid(double)},
  {  "ratio", &typeid(double)},
  {  "eps", &typeid(double)},
  {  "beginend", &typeid(KN<long> *  )},
  {  "file", &typeid(string * )}
};

int IsoLineK(R2 *P,double *f,R2 *Q,int *i0,int *i1,double eps)
{
  
  int kv=0,ke=0,e=3;
  int tv[3],te[3],vk[3];
  for(int i=0;i<3;++i)
    {
      if( abs(f[i]) <= eps) {
	e -= tv[kv++]=i;
	vk[i]=1;
      }
      else 
	vk[i]=0;
    }
  if(kv>1) //  on vertex on the isoline ....
    {
      if(kv==2 && f[e] > 0.)
	{  // pb d'unicity, need to see the adj triangle ... 
	  return 10+e ; // edge number + 10
	}
      else if (kv == 1)
	{
	  int i = tv[0]; 
	  int i1= (i +1)%3;
	  int i2= (i +2)%3;
	  if( f[i1]*f[i2] < 0.)  // intersect 
	    {
	      ffassert(0); 
	    }
	  else return 0; 
	}
    }
  else // see internal edge .. 
    for(int e=0;e<3;++e)
      {
	int j0=(e+1)%3;
	int j1=(e+2)%3;
	if( vk[j0]) //  the intial  point on iso line
	  {
	    if(0. < f[j1])	    
	      te[ke]=e,i0[ke]=j0,i1[ke]=j0,++ke;     
	    else 
	      te[ke]=e+3,i0[ke]=j0,i1[ke]=j0,++ke;
	  }
	else if (vk[j1]); // skip the final point on iso line
	else if( f[j0] < 0. && 0. < f[j1])  // good  sens
	  te[ke]=e,i0[ke]=j0,i1[ke]=j1,++ke;     
	else if ( f[j0] > 0. && 0. > f[j1]) // inverse  sens 
	  te[ke]=e+3,i0[ke]=j1,i1[ke]=j0,++ke;
      }
  if( ke==2) 
    {
      // the  K[i1[0]] , Q[0], Q[1] must be direct ...  
      // the  K[i0[1]] , Q[0], Q[1] must be direct ...  
      // Warning   no trivail case ..  make a plot to see 
      //  with is good
      // the first edge must be
 
      if(te[0]<3)  // oriente the line 
	{ 
	  assert(te[1] >=3);
	  std::swap(te[0],te[1]);
	  std::swap(i0[0],i0[1]);
	  std::swap(i1[0],i1[1]);
	  if(debug) cout << " swap " << endl;
	}
      for(int i=0;i<2;++i)
	{
	  int j0=i0[i],j1=i1[i];
	  if( j0== j1)
	    Q[i] = P[j0];
	  else
	    Q[i] = (P[j0]*(f[j1]) -  P[j1]*(f[j0]) ) /(f[j1]-f[j0]);
	  if(debug) cout << i << " " << j0 << " " << j1 << " : " 
	  
	   << Q[i] << "***" << endl;
	}
	if(debug)
	{
	  cout << "i0 " << i0[0] << " " << i0[1] << " " << det(P[i1[0]],Q[0],Q[1]) <<endl;
	  cout << "i1 " << i1[0] << " " << i1[1] << " " << det(P[i0[1]],Q[1],Q[0]) <<endl;
	  cout << "f " << f[0] << " " << f[1] << " " << f[2] << endl;
	  cout << "P " << P[0] << ", " << P[1] << ", " << P[2] << endl;
	  cout << "Q " << Q[0] << ", " << Q[1]  << endl;
	}
	if(!vk[i1[0]])
	assert( det(P[i1[0]],Q[0],Q[1]) > 0);
	if(!vk[i0[1]])
	assert( det(P[i0[1]],Q[1],Q[0]) > 0);
	return 2;
    }
  // remark, the left of the line is upper . 
  return 0;
}
int LineBorder(R2 *P,double *f,long close,R2 *Q,int *i1,int *i2,double eps)
{
  int np=0;
  if(close)
    {
      if(f[0]>-eps)
	{ 
	  Q[np]=P[0];
	  i1[np]=i2[np]=0,np++;
	  
	}
	     if(f[1]>-eps) 
	{
	  Q[np]=P[1];
	  i1[np]=i2[np]=1,np++;
	}

      if (f[0]*f[1] <= - eps*eps)
	{
	  Q[np]= (P[0]*(f[1]) -  P[1]*(f[0]) ) /(f[1]-f[0]);
	  i1[np]=0,i2[np]=1,np++;
	}
     }
  else
    {
    }
  return np;
}
struct R2_I2 {
  R2 P;
  int nx;
  R2_I2(R2 A,int nxx=-1) : P(A),nx(nxx) {}
  bool add(int k0,int k1,multimap<int,int> & L) {
    if (nx=-1)  nx=k1;
    else {
      if(nx >0) {//  more than 2 seg ... put data in  the multi map .. 
	L.insert(make_pair(k0,nx));
	L.insert(make_pair(k0,k1));
	nx=-2;
      }
      else 
	L.insert(make_pair(k0,k1));
      return false;
    }
    return true;}
  
int next(int k0,multimap<int,int> & L,int rm=0) 
  {
    int nxx=-1;
    if(nx>=0) {
      nxx=nx;
      if(rm) nx=-2;}
    else {
      typedef multimap<int,int>::iterator IT;
      IT  f= L.find(k0);
      if(f == L.end()) 
	nxx=-1; // 
      else
	{
	  nxx = f->second;
	  if(rm) 
	    L.erase(f); 
	}
      
    }
    return  nxx; 
  }
  int count(int k0,const multimap<int,int> & L) const 
  {
    if(nx>=0) return 1;
    else return L.count(k0);
  }
};

AnyType ISOLINE_P1_Op::operator()(Stack stack)  const 
{
  MeshPoint *mp(MeshPointStack(stack)) , mps=*mp;
  KNM<double> * pxy=0;
  KN<double> * pxx=0; 
  KN<double> * pyy=0;
  if(exy)
    pxy = GetAny<KNM<double>*>((*exy)(stack));
  if(exx)
    pxx = GetAny<KN<double>*>((*exx)(stack));
  if(eyy) 
    pyy = GetAny<KN<double>*>((*eyy)(stack));
  //  cout << pxx << " " << pyy << " " << pxy << endl; 
  ffassert( (pxx || pyy) ==  !pxy ) ;  
  Mesh * pTh= GetAny<Mesh *>((*eTh)(stack));
  ffassert(pTh);
  Mesh &Th=*pTh;
  int nbv=Th.nv; // nombre de sommet 
  int nbt=Th.nt; // nombre de triangles
  int nbe=Th.neb; // nombre d'aretes fontiere
  long nbc;
  // value of isoline
  int ka=0;
  double isovalue = arg(0,stack,0.);
  long close      = arg(1,stack,1L);
  double smoothing   = arg(2,stack,0.);
  double  ratio  = arg(3,stack,1.);
  double  epsr    = arg(4,stack,1e-10);
  KN<long> *pbeginend = arg(5,stack, (KN<long> *) 0); 
  string * file = arg(6,stack,(string*) 0);
  vector< R2_I2 >  P;
  multimap<int,int> L; 

  map<pair<int,int>, int> FP;
  const   double unset = -1e-100;
  KN<double> tff(nbv, unset);

  // loop over triangle 
  for (int it=0;it<Th.nt;++it)
    {
      for( int iv=0; iv<3; ++iv)
	{
	  int i=Th(it,iv);  
	  if(tff[i]==unset){
	    mp->setP(&Th,it,iv);
	    tff[i]=GetAny<double>((*eff)(stack))-isovalue;
	  }	  
	}
    }
  if(close<0) 
    {
      tff=-tff;
      close=-close; 
    }

         
    
  double tffmax=tff.max(), tffmin=tff.min(); 
  if(verbosity)
    cout << " -- iosline close="<< close << " iso= "<< isovalue << " " << epsr << endl
	 << " bound - iso " << tffmin << " " << tffmax << endl;
  double eps = (tffmax-tffmin)*epsr;
  if( (tffmax <0.) || (tffmin > 0.)) return 0L; 
  if(epsr>0) eps = epsr;
  for (int k=0;k<Th.nt;++k)
    {
      int iK[3]={Th(k,0),Th(k,1),Th(k,2)};
      R2 Pk[3]={Th(iK[0]),Th(iK[1]),Th(iK[2]) };
      R  fk[3]={tff[iK[0]],tff[iK[1]],tff[iK[2]] };
      R2 Qk[6];
      int i1[6],i2[6];
      int np=IsoLineK(Pk,fk,Qk,i1,i2,eps);
      if(np==2)
	{
	  for(int i=0;i<np;++i) // sort i1,i2 .. 
	  {
	  	i1[i]=iK[i1[i]];
	  	i2[i]=iK[i2[i]];
	  	if(i2[i]<i1[i]) 
	  	  std::swap(i1[i],i2[i]);
	  	  	  	
	  }	
	  int p[2]; // point number
	  for(int i=0;i<2;++i)
	    {
	      pair< map<pair<int,int>, int>::iterator , bool > ii;
	      pair<int,int> e(i1[i],i2[i]);
	      ii=FP.insert(make_pair(e,P.size()));
	      if(ii.second)
		P.push_back(R2_I2(Qk[i])); 
	      if(debug) cout << i1[i] << " ---  " << i2[i] <<" ;     " << ii.second<<  endl;
 
	      p[i] = ii.first->second ;
	    }
	  // add line k[0], k[1]
	  P[p[0]].add(p[0],p[1],L);
	  if(debug)
	  cout <<  " +++ " << Qk[0] << " ->  " << Qk[1] << " :: " << p[0] << " -> " << p[1] << endl;

	}	  	
    }

  if(close)
    {
      if(debug) cout<< " Close path " << endl; 
      for (int k=0;k<Th.nt;++k) 
	{
	  Triangle &K=Th[k];
	  for(int e=0;e<3;++e)
	    {
	      int ee,kk=Th.ElementAdj(k,ee=e);
	      if( kk==k || kk<0) 
		{ // true border element edge
		  int iK[2]={Th(k,(e+1)%3),Th(k,(e+2)%3)};
		  R2 Pk[2]={Th(iK[0]),Th(iK[1])};
		  R  fk[2]={tff[iK[0]],tff[iK[1]]};
		  R2 Qk[2];
		  int i1[2],i2[2];
		  if(debug) cout << " LB : " << Pk[0] << ", " << fk[0] << " ->  "<< Pk[1] << ", " << fk[1] 
				 << " : " <<iK[0] << " " << iK[1]  << endl; 
		  int np=LineBorder(Pk,fk,close,Qk,i1,i2,eps);
		  //		  cout << np << endl; 
		  if(np>=10) 
		    { // full edge 
		      int ke =  0; 

		    }
		  else if(np==2)
		    {
		      for(int i=0;i<2;++i)
			{
			  i1[i]=iK[i1[i]];
			  i2[i]=iK[i2[i]];
			  if(i2[i]<i1[i]) 
			    std::swap(i1[i],i2[i]);
			}
		      if(debug) cout << " add  : " <<  Qk[0] << ", " << i1[0] <<',' << i2[0] 
				     << " ->   " << Qk[1] << ", " << i1[1] <<',' << i2[1]<< endl; 
		      int p[2]; // point number                                                                                                        
		      for(int i=0;i<2;++i)
			{
			  pair< map<pair<int,int>, int>::iterator , bool > ii;
			  pair<int,int> ee(i1[i],i2[i]);
			  ii=FP.insert(make_pair(ee,P.size()));
			  if(ii.second)
			    P.push_back(R2_I2(Qk[i]));
			  if(debug) cout << i1[i] << " ---  " << i2[i] <<" ;     " << ii.second<<  endl;
			  p[i] = ii.first->second ;
			}
		      // add line k[0], k[1]                                                                                                           
		      P[p[0]].add(p[0],p[1],L);
		      if(debug)
			cout <<  " +++ " << Qk[0] << " ->  " << Qk[1] << " :: " << p[0] << " -> " << p[1] << endl;
		    }

		}
	      
	    }
	}
      if(debug) cout<< " End Close path " << endl; 

    }
  if(verbosity>99)
    { // dump the data base
      cout << " IsolineP1 " << endl;
      for(int i=0;i<P.size();++i)
	cout << "\t" <<i << " :  " <<  P[i].P << " ->  " << P[i].nx  << endl;
      cout<< " multmap for execption " << endl;
      for( multimap<int,int>::const_iterator i= L.begin(); i!= L.end();++i)
	     cout <<  "\t" << i->first << " -> " << i->second << endl;
      cout<< " End multmap for execption " << endl;
    }
    
  vector<int> iQ, QQ; // QQ the list of curve
  int np= P.size();
  KN<int> start(np);
  start=-1;
  int kkk=0;
  while(1)
    {
      assert(kkk++ < 10);
      int kk=0,k=0;
      for(int i=0;i<np;++i)
	{ 
	  kk+= P[i].count(i,L);
	  if( P[i].nx>=0) 
	    start[P[i].nx]=i;
	  else if(start[i] == -1)
	    k++;
	}
	if(verbosity>19)
	  cout << "re starting:  k = " << k << " " << kk << endl;	
	if(kk==0) 
	  break;
      for(int i=0;i<np;++i)
	if( (start[i]==-1) || ( (k==0) && (start[i]>=0)) )
	  {
	   if(verbosity>9)
	    cout << "  isolineP1: start curve  = " << i << " -> " << P[i].next(i,L,0) ;
	   if(verbosity>99)
	    cout << " (" << P[i].nx <<endl;
	    iQ.push_back(QQ.size());
	    QQ.push_back(i);
	    start[i]=-2;
	    int i0=i,i1=0,ie=i;
	    while(1)
	      {
		i1= P[i0].next(i0,L,1);
		if(i1<0) break; 
		if(start[i1]<0) break; 
		
		QQ.push_back(i1);
		if(verbosity>99)
		cout << " " << i1;
		start[i1]=-2;  
		i0=i1;
	      }
	    if(i1==ie) // close the path .. 
	      {
	    	QQ.push_back(i1);
		if(verbosity>99)
		  cout << " " << i1;
	      }
	    
	    if(verbosity>99)
	    cout << ") "<< endl;
	    else if(verbosity>9) cout << endl; 
	    iQ.push_back(QQ.size());

	  }
    }
  // sort iQ
  if(iQ.size()>2)
    {
      vector< pair<int , pair<int,int> > > sQ(iQ.size()/2);
      for(int i=0,j=0; i<iQ.size();i+=2,++j)
	{
	  int i0=iQ[i];
	  int i1=iQ[i+1];
	  sQ[j] = make_pair(i0-i1,make_pair(i0,i1));
	}
      std::sort(sQ.begin(), sQ.end());
      for(int i=0,j=0; i<iQ.size();i+=2,++j)
	{
	  iQ[i]   =sQ[j].second.first;
	  iQ[i+1] =sQ[j].second.second;
	}
    }
if(smoothing>0)
  {
   KN<R2> P1(QQ.size()),P2(QQ.size());
   for(int i=0;i<QQ.size();++i)
	P1[i]=  P[QQ[i]].P ;  	
   
   
  	//  Smoothing the curve  
    double c1=1, c0=4, ct=2*c1+c0;
    c1/=ct;
    c0/=ct;

   
    for(int i=0; i<iQ.size();)
      {  
	int i0=iQ[i++];
	int i1=iQ[i++]-1; 
	int nbsmoothing = pow((i1-i0),ratio)*smoothing;
	if(verbosity>2) 
	  cout << "     curve " << i << " size = " << i1-i0 << " nbsmoothing = " << nbsmoothing << " " << i0 << " " << i1 << endl;
	P2=P1;
	for( int step=0;step<nbsmoothing;++step)
	  {
	    for(int j=i0+1;j<i1;++j)
	      P2[j] = c0*P1[j] +c1*P1[j-1] +c1*P1[j+1];
	    
	    if(QQ[i0]==QQ[i1]) // close curve 
	      {
		int j0= i0+1; // prec 
		int j1=  i1-1;  // next 
		P2[i0]=P2[i1] = c0*P1[i0] +c1*P1[j0] +c1*P1[j1];
	      }
	    P1=P2;
	  }
      }
    
    for(int i=0;i<QQ.size();++i)
      P[QQ[i]].P=P1[i]   ;  	
    
  }
 
  if(pbeginend)
  {
  	 pbeginend->resize(iQ.size());
  	 for(int i=0; i<iQ.size();++i)
  	  (*pbeginend)[i]=iQ[i];
  }
  
  if(pxx && pyy )
    {
      pxx->resize(QQ.size());      
      pyy->resize(QQ.size());
      for(int i=0;i<QQ.size();++i)
	{
	  int j=QQ[i];
	  (*pxx)[i]=  P[j].P.x ;
	  (*pyy)[i]=  P[j].P.y ;
	  
	}
    }
  else if(pxy) 
    {
      pxy->resize(3,QQ.size()); 

      
      for(int k=0; k<iQ.size();k+=2)
	{
	  int i0=iQ[k],i1=iQ[k+1];
	  double lg=0;
	  R2 Po=P[QQ[i0]].P;
	  for(int i=i0;i<i1;++i)
	    {
	      int j=QQ[i];
	      (*pxy)(0,i)=  P[j].P.x ;
	      (*pxy)(1,i)=  P[j].P.y ;
	      lg += R2(P[j].P,Po).norme() ; 
	      (*pxy)(2,i)=  lg; 
	      Po = P[j].P;
	    }
	}
    }
  else ffassert(0); 
  

  nbc = iQ.size()/2; 
  if(file) 
  {
    
    ofstream fqq(file->c_str());
    int i=0,i0,i1,n2= iQ.size(),k=0;
    while(i<n2)
      { 
	k++;
	i0=iQ[i++];
	i1=iQ[i++];
	cout<< i0 << " " << i1 << endl;
	for(int l=i0;l<i1;++l)
	{
	  int j=QQ[l];	
	  fqq << P[j].P.x << " " << P[j].P.y << " " << k << " "<< j <<  endl;
	}
	fqq << endl;
      }
  }
  /*
  int err=0;
  int pe[10];
  for(int i=0;i< P.size();++i)
    {
      int pb = P[i].count(i,L);
      if(pb) 
	if(err<10) 
	  pe[err++]=i;
	    else
	      err++;
	}
      
      if(err>0)
	{
	  for(int i=0;i<10;i++)
	    cout << " PB point = " << pe[i] << " " << P[pe[i]].P << " odd count = " << P[pe[i]].count(pe[i],L)  << endl;
	  ffassert(0);
	}
  */  
      // construction des courble 
  
  return nbc;
}

class  ISOLINE_P1: public OneOperator { public:  
typedef Mesh *pmesh;
  int cas;

  ISOLINE_P1() : OneOperator(atype<long>(),atype<pmesh>(),atype<double>(), atype<KN<double>*>(),atype<KN<double>* >() ) ,cas(4){}
  ISOLINE_P1(int ) : OneOperator(atype<long>(),atype<pmesh>(),atype<double>(), atype<KNM<double>*>() ),cas(3) {}
  
  E_F0 * code(const basicAC_F0 & args) const 
  { 
    if(cas==4)
      return  new ISOLINE_P1_Op( args,
				 t[0]->CastTo(args[0]),
				 t[1]->CastTo(args[1]),
				 t[2]->CastTo(args[2]),
				 t[3]->CastTo(args[3]) ); 
    else if(cas==3) 
      return  new ISOLINE_P1_Op( args,
				 t[0]->CastTo(args[0]),
				 t[1]->CastTo(args[1]),
				 t[2]->CastTo(args[2]) ); 
    else ffassert(0); // bug 
  }


  


};



R3  * Curve(Stack stack,const KNM_<double> &b,const  long &li0,const  long & li1,const double & ss)
{
  assert(b.N() >=3);
  int i0=li0,i1=li1,im;
  if(i0<0) i0=0;
  if(i1<0) i1=b.M()-1;
  double lg=b(2,i1);
  R3 Q;
  ffassert(lg>0 && b(2,0)==0.); 
  double s = ss*lg;
  int k=0,k1=i1;
  while(i0 < i1-1)
    {
      ffassert(k++ < k1);
      im = (i0+i1)/2;
      if(s <b(2,im)  )
	{ i1=im;
	}
      else if(s>b(2,im)  )
	{ i0= im;
	}
      else {  Q=R3(b(0,im),b(1,im),0);  i0=i1=im;break;}
    }
  if(i0<i1)
    {
      ffassert(b(2,i0) <= s );
      ffassert(b(2,i1) >= s );
      R2 A(b(0,i0),b(1,i0));
      R2 B(b(0,i1),b(1,i1));
      double l1=(b(2,i1)-s);
      double l0=s-b(2,i0);
      Q= (l1*A + l0*B)/(l1+l0);
    }
  R3 *pQ = Add2StackOfPtr2Free(stack,new R3(Q));
  // MeshPoint &mp= *MeshPointStack(stack); // the struct to get x,y, normal , value 
  //mp.P.x=Q.x; // get the current x value
  //mp.P.y=Q.y; // get the current y value
  return pQ; 
}

R3   * Curve(Stack stack,const KNM_<double> &b,const double & ss)
{
  return Curve(stack,b,-1,-1,ss);
}



void finit()
{  
 
  typedef Mesh *pmesh;
  
  Global.Add("isoline","(",new ISOLINE_P1);
  Global.Add("isoline","(",new ISOLINE_P1(1));
  Global.Add("Curve","(",new OneOperator2s_<R3*,KNM_<double>,double>(Curve));
  Global.Add("Curve","(",new OneOperator4s_<R3*,KNM_<double>,long,long,double>(Curve));
}

LOADFUNC(finit);  //  une variable globale qui serat construite  au chargement dynamique    
