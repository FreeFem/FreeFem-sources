// -*- Mode : c++ -*-
//
// SUMMARY  :      
// USAGE    :        
// ORG      : 
// AUTHOR   : Frederic Hecht
// E-MAIL   : hecht@ann.jussieu.fr
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
 */
#include <fstream>
#include <iostream>
#include <cassert>
#include <algorithm>
#include <complex>


using namespace std;
#include "error.hpp"
#include "arpackff.hpp"
#include "AFunction.hpp"
#include "rgraph.hpp"
#include "RNM.hpp"
#include "MatriceCreuse_tpl.hpp"
#include "Mesh3dn.hpp"
#include "MeshPoint.hpp"
#include "lgfem.hpp"
#include "lgmesh3.hpp"
#include "lgsolver.hpp"
#include "problem.hpp"
extern Block *currentblock;

typedef double R;
const bool dddd=false;

template<class K>
void  cB(long mode,int n,K *x,K *y,const MatriceCreuse<K>  & M)
{
  KN_<K> X(x,n);
  KN_<K> Y(y,n);
  switch(mode)
    {
    case 1:  
      Y= X; //
      break;
    case 2:  
    case 3:  
    case 4:  
    case 5:  
      Y = M*X;// OP = inv(M) A
      break;
    default:
      ffassert(0); // a faire 
    }
  
}
template<class K>
void  cOP(long mode,int n,K *x,K *y,KN<K> W,const MatriceCreuse<K> & OP,const MatriceCreuse<K>  & M)
{ 
  //  y entre , y = OP * x  sortie, w working arrayy
  /*
    c  Mode 1:  ===> OP = A  and  B = I.
    c
    c  Mode 2: 
    c           ===> OP = inv[M]*A  and  B = M.
    c
    c  Mode 3: 
    c           ===> OP = (inv[K - sigma*M])*M  and  B = M.  // DSAUPD
    c           ===> OP = Real_Part{ inv[A - sigma*M]*M }  and  B = M. //DNAUPD 
    c           ===> OP =  inv[A - sigma*M]*M   and  B = M.  // ZNAUPD
    c  
    c  Mode 4: 
    c           ===> OP = (inv[K - sigma*KG])*K  and  B = K. // DSAUPD  
    c           ===> OP = Imaginary_Part{ inv[A - sigma*M]*M }  and  B = M. //DNAUPD 
    c  Mode 5:
    c           ===> OP = inv[A - sigma*M]*[A + sigma*M]  and  B = M. // DSAUPD      
  */

  KN_<K> X(x,n);
  KN_<K> Y(y,n);
  switch(mode)
    {
    case 1:  
      Y= OP*X; //
      break;
    case 2:  
      W = OP*X;// OP = inv(M) A
      M.Solve(Y,W);
      break;
    case 3: // OP inv( A-sigma  B)    
      OP.Solve(Y,X); // 
      break;
    case 4: // OP image inv( A-sigma  B)    
      OP.Solve(Y,X); // 
      break;
    default:
      ffassert(0); // a faire 
    }
}

class EigenValue : public OneOperator
{ public:
  typedef R K;
  typedef KN<K> Kn;
  typedef KN_<K> Kn_;
  const int cas;
  class E_EV: public E_F0mps { public:
      const int cas;
    
    static basicAC_F0::name_and_type name_param[] ;
    static const int n_name_param =12;
    Expression nargs[n_name_param];
    Expression expOP1,expB;
    template<class T>
    T arg(int i,Stack stack,const T & a) const{ return nargs[i] ? GetAny<T>( (*nargs[i])(stack) ): a;}
    E_EV(const basicAC_F0 & args,int cc) :
      cas(cc)
    {
      // OP1 = (A-sigma*B)        
      //                int nbj= args.size()-1;
      args.SetNameParam(n_name_param,name_param,nargs);
      expOP1=to< Matrice_Creuse<K> *>(args[0]);
      expB=to< Matrice_Creuse<K> *>(args[1]);
      
    }
    
    AnyType operator()(Stack stack)  const;
    operator aType () const { return atype<long>();}         
    
  };
  
  E_F0 * code(const basicAC_F0 & args) const {
    return new E_EV(args,cas);}
  
  EigenValue(int c) :   
    OneOperator(atype<long>(),
		atype<Matrice_Creuse<K> *>(),
		atype<Matrice_Creuse<K> *>()),
    cas(c){}
    
};

class EigenValueC : public OneOperator
{ public:
  typedef Complex K;
  typedef double R;
  typedef KN<K> Kn;
  typedef KN_<K> Kn_;
  const int cas;
  class E_EV: public E_F0mps { public:
      const int cas;
    
    static basicAC_F0::name_and_type name_param[] ;
    static const int n_name_param =10;
    Expression nargs[n_name_param];
    Expression expOP1,expB;
    template<class T>
    T arg(int i,Stack stack,const T & a) const{ return nargs[i] ? GetAny<T>( (*nargs[i])(stack) ): a;}
    E_EV(const basicAC_F0 & args,int cc) :
      cas(cc)
    {
      // OP1 = (A-sigma*B)        
      //                int nbj= args.size()-1;
      args.SetNameParam(n_name_param,name_param,nargs);
      expOP1=to< Matrice_Creuse<K> *>(args[0]);
      expB=to< Matrice_Creuse<K> *>(args[1]);
      
    }
    
    AnyType operator()(Stack stack)  const;
    operator aType () const { return atype<long>();}         
    
  };
  
  E_F0 * code(const basicAC_F0 & args) const {
    return new E_EV(args,cas);}
  
  EigenValueC(int c) :   
    OneOperator(atype<long>(),
		atype<Matrice_Creuse<K> *>(),
		atype<Matrice_Creuse<K> *>()),
    cas(c){}
  
};

basicAC_F0::name_and_type  EigenValue::E_EV::name_param[]= {
  {   "tol", &typeid(double)  },
  {   "nev",&typeid(long) },
  {   "sym",&typeid(bool)},
  {   "sigma",&typeid(double)},
  {   "value",&typeid(KN<double> *)},
  {   "vector",&typeid(pferarray) }, 
  {   "ncv",&typeid(long) }, // the number of Arnoldi vectors generated 
  {   "maxit",&typeid(long)}, // the maximum number of Arnoldi iterations 
  {   "ivalue",&typeid(KN<double> *)},
  {   "rawvector",&typeid(KNM<double> *) },
  {   "resid",&typeid(KN<double> *) },
  {   "mode",&typeid(long) } // 12 ieme
  
  
};

basicAC_F0::name_and_type  EigenValueC::E_EV::name_param[]= {
  {  "tol", &typeid(double)  },
  {  "nev",&typeid(long) },
  {  "sigma",&typeid(K)},
  {  "value",&typeid(KN<Complex> *)},
  {  "vector",&typeid(pfecarray) }, 
  {  "ncv",&typeid(long) }, // the number of Arnoldi vectors generated 
  {  "maxit",&typeid(long)}, // the maximum number of Arnoldi iterations 
  {  "rawvector",&typeid(KNM<Complex> *) }, 
  {  "resid",&typeid(KN<Complex> *)},
  {   "mode",&typeid(long) } // 10 ieme   
  
};


AnyType EigenValue::E_EV::operator()(Stack stack)  const {
  double tol=1e-6;
  long nconv=0; 
  long nbev=1;
  bool sym=false;
  long ncv =0;  // the number of Arnoldi vectors generated 
  long maxit=0;  // the maximum number of Arnoldi iterations 
  double sigma=0;
  int mode = 3; //  
  KN<double> * evalue=0;
  KN<double> * evaluei=0;
  KN<double> * resid=0;
  KNM<double> * rawvector=0;
  double ws,vs;           // for debugging FH ++++++++++
  pferarray  evector2;
  pferbasearray   evector=0;
  tol=arg<double>(0,stack,0);
  nbev=arg<long>(1,stack,10);
  sym=arg<bool>(2,stack,false);
  sigma=arg<double>(3,stack,0.0);
  evalue=arg<KN<double> *>(4,stack,0);
  evector2 =arg<pferarray>(5,stack,make_pair<pferbasearray,int>(0,0)); 
  ncv= arg<long>(6,stack,0);
  maxit= arg<long>(7,stack,0);
  evaluei=arg<KN<double> *>(8,stack,0);
  rawvector=arg<KNM<double> *>(9,stack,0);
  resid=arg<KN<double> *>(10,stack,0);
  mode = arg<long>(11,stack,3);
    
  evector=evector2.first;
  Matrice_Creuse<K> *pOP1 =  GetAny<Matrice_Creuse<K> *>((*expOP1)(stack));
  Matrice_Creuse<K> *pB =  GetAny<Matrice_Creuse<K> *>((*expB)(stack));
  double * residptr=resid? (double*) *resid : 0;
  //cout << " residptr = " << residptr <<endl;
  
  if(evalue) nbev=Max( (long)evalue->N(),nbev);
  
  const MatriceCreuse<K> & OP1 = pOP1->A;
  const MatriceCreuse<K> & B = pB->A;
  
  long  n=OP1.n;
  if(sym)
    {
      nbev=min(n-1,nbev);
      if(!ncv) 	ncv = min(nbev*2+1,n);
    }
  else
    {
      nbev=min(nbev,n-2);
      if(!ncv)     ncv = nbev*2+1;        
    }
  ncv = max(nbev+2,ncv);
  ncv = min(ncv,n);

  if(!maxit)  maxit = 100*nbev;

    
    /* daupp

     c  Mode 1:  A*x = lambda*x, A symmetric 
     c           ===> OP = A  and  B = I.
     c
     c  Mode 2:  A*x = lambda*M*x, A symmetric, M symmetric positive definite
     c           ===> OP = inv[M]*A  and  B = M.
     c           ===> (If M can be factored see remark 3 below)
     c
     c  Mode 3:  K*x = lambda*M*x, K symmetric, M symmetric positive semi-definite
     c           ===> OP = (inv[K - sigma*M])*M  and  B = M. 
     c           ===> Shift-and-Invert mode
     c
     c  Mode 4:  K*x = lambda*KG*x, K symmetric positive semi-definite, 
     c           KG symmetric indefinite
     c           ===> OP = (inv[K - sigma*KG])*K  and  B = K.
     c           ===> Buckling mode
     c
     c  Mode 5:  A*x = lambda*M*x, A symmetric, M symmetric positive semi-definite
     c           ===> OP = inv[A - sigma*M]*[A + sigma*M]  and  B = M.
     c           ===> Cayley transformed mode
     
     c  WHICH   Character*2.  (INPUT)
     c          Specify which of the Ritz values of OP to compute.
     c
     c          'LA' - compute the NEV largest (algebraic) eigenvalues.
     c          'SA' - compute the NEV smallest (algebraic) eigenvalues.
     c          'LM' - compute the NEV largest (in magnitude) eigenvalues.
     c          'SM' - compute the NEV smallest (in magnitude) eigenvalues. 
     c          'BE' - compute NEV eigenvalues, half from each end of the
     c                 spectrum.  When NEV is odd, compute one more from the
     c                 high end than from the low end.

     c  NEV     Integer.  (INPUT)
     c          Number of eigenvalues of OP to be computed. 0 < NEV < N.

     c  NCV     Integer.  (INPUT)
     c          Number of columns of the matrix V (less than or equal to N).

     */
    
    /* Naupp
     Mode 1:  A*x = lambda*x.
     c           ===> OP = A  and  B = I.
     c
     c  Mode 2:  A*x = lambda*M*x, M symmetric positive definite
     c           ===> OP = inv[M]*A  and  B = M.
     c
     c  Mode 3:  A*x = lambda*M*x, M symmetric semi-definite
     c           ===> OP = Real_Part{ inv[A - sigma*M]*M }  and  B = M. 
     c           ===> shift-and-invert mode (in real arithmetic)

     c  Mode 4:  A*x = lambda*M*x, M symmetric semi-definite
     c           ===> OP = Imaginary_Part{ inv[A - sigma*M]*M }  and  B = M. 
     c           ===> shift-and-invert mode (in real arithmetic)



     c          IDO =  0: first call to the reverse communication interface
     c          IDO = -1: compute  Y = OP * X  where
     c                    IPNTR(1) is the pointer into WORKD for X,
     c                    IPNTR(2) is the pointer into WORKD for Y.
     c                    This is for the initialization phase to force the
     c                    starting vector into the range of OP.
     c          IDO =  1: compute  Y = OP * Z  and Z = B * X where
     c                    IPNTR(1) is the pointer into WORKD for X,
     c                    IPNTR(2) is the pointer into WORKD for Y,
     c                    IPNTR(3) is the pointer into WORKD for Z.
     c          IDO =  2: compute  Y = B * X  where
     c                    IPNTR(1) is the pointer into WORKD for X,
     c                    IPNTR(2) is the pointer into WORKD for Y.
     c          IDO =  3: compute the IPARAM(8) real and imaginary parts 
     c                    of the shifts where INPTR(14) is the pointer
     c                    into WORKL for placing the shifts. See Remark
     c                    5 below.
     c          IDO =  4: compute Z = OP * X
     c          IDO = 99: done
     c 
     
     c  WHICH   Character*2.  (INPUT)
     c          'LM' -> want the NEV eigenvalues of largest magnitude.
     c          'SM' -> want the NEV eigenvalues of smallest magnitude.
     c          'LR' -> want the NEV eigenvalues of largest real part.
     c          'SR' -> want the NEV eigenvalues of smallest real part.
     c          'LI' -> want the NEV eigenvalues of largest imaginary part.
     c          'SI' -> want the NEV eigenvalues of smallest imaginary part.
     c
     
     c  NCV     Integer.  (INPUT)
     c          Number of columns of the matrix V. NCV must satisfy the two
     c          inequalities 2 <= NCV-NEV and NCV <= N.
     
     c  NEV     Integer.  (INPUT)
     c          Number of eigenvalues of OP to be computed. 0 < NEV < N-1.
     
     
     */
    
    
    const char *serr[10];
    int err=0;
    if( ! (nbev < n) )
	serr[err++]="  Number of eigenvalues of OP to be computed nev <= n ";
    
    if( (mode < 1 || mode > 5) && sym) 
        serr[err++]="  the mode = 1 ,2 ,3, 4, 5  ";
    if( (mode < 1 || mode > 4) && !sym) 
        serr[err++]="  the mode = 1 ,2 ,3, 4  ";
    // 2 <= NCV-NEV and NCV <= N
    
    if( ! (   ncv <= n) && sym ) 
	serr[err++]="   ( ncv <= n) (symetric  case) ";
    if( ! ( (   ncv <= n) && 2 <= (ncv-nbev ) ) && !sym )
	serr[err++]="   ( ncv <= n)  2 <= (ncv-nev ) ( no-symetric  case) ";
    
    if (n != OP1.m) 
	serr[err++]=" the first matrix in EigneValue is not square.";
    if (n != B.n ) 
	serr[err++]="Sorry the row's number of the secand matrix in EigneValue is wrong.";
    if (n != B.m ) 
	serr[err++]="Sorry the colum's number of the secand matrix in EigneValue is wrong.";
	
   if(verbosity)
    if(sym)
      cout << "Real symmetric eigenvalue problem: A*x - B*x*lambda" << endl;
    else
      cout << "Real non symmetric eigenvalue problem: A*x - B*x*lambda" << endl;
  
        
    if(verbosity>9 || err)
	cout << "    n " << n << ", nev "<< nbev << ", tol =" << tol << ", maxit =" << maxit 
	<< ", ncv = " <<ncv << ", mode = " << mode << endl;
    if(err) 
      {
	  cerr << " list of the error " << endl;
	  for (int i=0;i<err;++i)
	      cerr << "\n-  " <<  serr[i] << endl;
	  ExecError("Fatal Error in EigenValue (real) ");
      }
    
  KN<K> work(n);

  if(sym)
    {
      int ido=0;
      char bmat= mode == 1 ? 'I' : 'G';
      char which[3]= "LM";	// larger value
      if(mode >2) which[0] ='S'; // smaller value
      int ishift=1; // Auto Shift true by default
      int iparam[12]= {0,ishift,0,maxit,1,nconv,0,mode,0,0,0,0};
      int ipntr[12]={ 0,0,0, 0,0,0,  0,0,0, 0,0,0 };
      KN<double> workd(3*n+1);
      int lworkl = ncv*(ncv+9);
      KN<double> workl(lworkl+1);
      KN<double> vp(ncv*n+1);

      int info= (residptr !=0);
      KN<double> vresid(residptr? 1: n);
      if(!residptr) residptr=&vresid[0];

      while (1) {
	
	saupp(ido,bmat,n,which,nbev,tol,  residptr,  ncv,  vp, n,
	      iparam,  ipntr,  workd,   workl,  lworkl, info);
/*
c          IDO =  0: first call to the reverse communication interface
c          IDO = -1: compute  Y = OP * X  where
c                    IPNTR(1) is the pointer into WORKD for X,
c                    IPNTR(2) is the pointer into WORKD for Y.
c                    This is for the initialization phase to force the
c                    starting vector into the range of OP.
c          IDO =  1: compute  Y = OP * Z  and Z = B * X where
c                    IPNTR(1) is the pointer into WORKD for X,
c                    IPNTR(2) is the pointer into WORKD for Y,
c                    IPNTR(3) is the pointer into WORKD for Z.
c          IDO =  2: compute  Y = B * X  where
c                    IPNTR(1) is the pointer into WORKD for X,
c                    IPNTR(2) is the pointer into WORKD for Y.
c          IDO =  3: compute the IPARAM(8) shifts where
c                    IPNTR(11) is the pointer into WORKL for
c                    placing the shifts. See remark 6 below.
*/
	if(verbosity>99) 
	  cout << "    saupp ido: " << ido << " info : " << info << endl;
	if(info<0) {cerr << " -- err arpack info = " << info << endl;}  
	sauppError(info);

	if(ido==99) break;
	
	K *xx=&workd[ipntr[1]];
	K *yy=&workd[ipntr[2]];
	K *zz=&workd[ipntr[3]];
	
	  switch (ido) {
	  case -1: 
	    cOP(mode,n,xx,yy,work,OP1,B);
	    break;
	  case  1:
	    cOP(mode,n,zz,yy,work,OP1,B);
	    cB(mode,n,xx,zz,B);
	    break;
	  case  2: 
	    cB(mode,n,xx,yy,B);	    
	    break;
	  default :
	    ffassert(0);
	  }
	  
      }
      nconv = iparam[5];
      if(nconv==0) cerr << " -- Strange: no eigens values ??? " << endl;
      // Finding eigenvalues and eigenvectors.
      if(nconv)
	{
	  KN<double> evr(nbev);
	  KNM<double> Z(n,nbev);
	  int ldz=n;
	  char HowMny ='A';
	  int rvec=1;
	  
	  seupp( rvec, HowMny, evr,  Z, ldz, sigma, bmat, n,
		 which,  nbev,  tol,  residptr, ncv, vp,  n, iparam,ipntr,  workd, workl,lworkl,info);
	  
	  if(verbosity>5)
	    {
	      
	      cout << "Dimension of the system            : " << n              << endl;
	      cout << "Number of 'requested' eigenvalues  : " << nbev  << endl;
	      cout << "Number of 'converged' eigenvalues  : " << nconv          << endl;
	      cout << "Number of Arnoldi vectors generated: " << ncv  << endl;
	      cout << "Number of iterations taken         : " << iparam[3]  << endl;
	      cout << endl;
	      
	      //if (prob.EigenvaluesFound()) {
	      cout << "Eigenvalues:" << endl;
	      for (int i=0; i<nconv; i++) {
		cout << "  lambda[" << (i+1) << "]: " << evr(i) << endl;
		  KN_<K> vi(Z(':',i)) ;
		  if(verbosity>99)
		    cout <<" Eigenvector: :" << vi << endl;
	      }
	      cout << endl;
	    }
	  
	  if (evalue)
	    {
	      KN<double> & ev(*evalue);
	      int m = Min(nconv,ev.N());
	      for(int i=0;i<m;i++)
		ev[i]=evr(i);
	      }
	  if(rawvector)
	    {
	      int m = Min(nconv,rawvector->M());
	      ffassert(rawvector->N()==n);
	      for(int i=0;i<m;i++)
		{
		  KN_<K> vi(Z(':',i)) ;
		  
		  //   cout << " ------ EV--raw " << vi.min() << " " << vi.max() << endl;
		  (*rawvector)(':',i)=vi;
		}
	      
	    }
	  
	  if (evector)
	    {
	      FEbaseArray<K,v_fes> & ev(*evector);
	      int m = Min(nconv,(long) ev.N);
	      for(int i=0;i<m;i++)
		{
		  FEbase<K,v_fes> & xx= **(ev[i]);
		    //if(xx.pVh->NbDoF != n)
		    //ExecError("Wrong Type size of FEspace to store the eigen vector ");
		    // if (xx.pVh != pOP1->pUh) 
		    //    ExecError("Wrong Type of FEspace to store the eigen vector ");
		    //xx.Vh = pOP1->Uh;
		  KN_<K> vi(Z(':',i)) ;
		  xx= new KN<K>(vi);
		  
		}
	    }
	}
    }
  else 
      { 
	// cas non symetrique ,
        // Finding an Arnoldi basis.                                                                                                               
        //int mode=3; //  Shift invert                                                                                                               
        int ido=0;
        char bmat='G';
        char which[]="LM";
        int ishift=1; // Auto Shift true by default                                                                                                
        int iparam[12]= {0,ishift,0,maxit,1,nconv,0,mode,0,0,0,0};
        int ipntr[15]={ 0,0,0, 0,0,0,  0,0,0, 0,0,0 ,0,0,0};
        KN<double> workd(3*n+1);
        int lworkl = 3*ncv*(ncv+2);
	KN<double> workl(lworkl+1);
        KN<double> vp(ncv*n+1);
        int info= (residptr !=0);
	KN<double> vresid(residptr? 1: n);
	if(!residptr) residptr=&vresid[0];
	if(verbosity>9)
	  cout << " n " << n << " nbev "<< nbev << " tol =" << tol << " maxit =" << maxit << " ncv = " <<ncv << endl;	
	
	while (1)
	{
	  
	  /*
	    
	    c          IDO =  0: first call to the reverse communication interface
	    c          IDO = -1: compute  Y = OP * X  where
	    c                    IPNTR(1) is the pointer into WORKD for X,
	    c                    IPNTR(2) is the pointer into WORKD for Y.
	    c                    This is for the initialization phase to force the
	    c                    starting vector into the range of OP.
	    c          IDO =  1: compute  Y = OP * Z  and Z = B * X where
	    c                    IPNTR(1) is the pointer into WORKD for X,
	    c                    IPNTR(2) is the pointer into WORKD for Y,
	    c                    IPNTR(3) is the pointer into WORKD for Z.
	    c          IDO =  2: compute  Y = B * X  where
	    c                    IPNTR(1) is the pointer into WORKD for X,
	    c                    IPNTR(2) is the pointer into WORKD for Y.
	    c          IDO =  3: compute the IPARAM(8) real and imaginary parts 
	    c                    of the shifts where INPTR(14) is the pointer
	    c                    into WORKL for placing the shifts. See Remark
	    c                    5 below.
	    c          IDO =  4: compute Z = OP * X
	    
	  */
	  naupp(ido,bmat,n,which,nbev,tol,  residptr,  ncv,  vp, n,
		iparam,  ipntr,  workd,   workl,  lworkl, info);
	  if(ido==99) break;

	  K *xx=&workd[ipntr[1]];
	  K *yy=&workd[ipntr[2]];
	  K *zz=&workd[ipntr[3]];

	  switch (ido) {
	  case -1: 
	    cOP(mode,n,xx,yy,work,OP1,B);
	    break;
	  case  1:
	    cOP(mode,n,zz,yy,work,OP1,B);
	    cB(mode,n,xx,zz,B);
	    break;
	  case  2: 
	    cB(mode,n,xx,yy,B);	    
	    break;
	    // case 5  je ne sais pas quoi faire... 
	  case 4:
	    cOP(mode,n,xx,zz,work,OP1,B);	    
	    break;
	  default :
	    cout << " bug ido " << ido << endl;
	    ffassert(0);
	  }
	sauppError(info);
	}
	nconv = iparam[5];
	if(nconv)
	  {
	    KN<double> evr(nbev+1), evi(nbev+1);
	    KNM<double> Z(n,nbev+1);
	    KN<double> workev(3*ncv);
	    int ldz=n;
	    char HowMny ='A';
	    int rvec=1;
	    double sigmai=0;
	    neupp( rvec, HowMny, evr,evi, Z , ldz, sigma,sigmai, workev,
		   bmat, n,   which,  nbev,  tol,  residptr, ncv,
		     vp, n, iparam , ipntr, workd, workl,lworkl,info);
	    
	    
	    if (verbosity)
	      {
		cout << "Real non-symmetric eigenvalue problem: A*x - B*x*lambda" << endl;
		cout << "mode " << mode << " sigma=" << sigma <<  endl << endl;
		
		
		cout << "Dimension of the system            : " << n              << endl;
		cout << "Number of 'requested' eigenvalues  : " << nbev  << endl;
		cout << "Number of 'converged' eigenvalues  : " << nconv          << endl;
		cout << "Number of Arnoldi vectors generated: " << ncv  << endl;
		cout << "Number of iterations taken         : " << iparam[3] << endl;
		cout << endl;
		  
		cout << "Eigenvalues:" << endl;
		for (int i=0; i<nconv; i++) {
		  cout << "  lambda[" << (i+1) << "]: " ;
		  ios::fmtflags oldflag= cout.flags();
		  cout.setf(ios::showpos);
		  cout << evr(i) 
		       <<  evi(i) << "i"<<  endl;
		  cout.flags(oldflag); // restore flags     
		  KN_<double> vi(Z(':',i)) ;
		  if(verbosity>99)
		    {
		      cout <<" Eigenvector: :" << vi << endl;
		      cout << endl;	
		    }    
		  }
	      }
	    
	    if (evalue)
	      {
		KN<double> & ev(*evalue);
		int m = Min(nconv,ev.N());
		for(int i=0;i<m;i++)
		  ev[i]=evr(i);
	      }
	    
	    if (evaluei)
	      {
		KN<double> & ev(*evaluei);
		int m = Min(nconv,ev.N());
		for(int i=0;i<m;i++)
		  ev[i]=evi(i);
	      }
	    if(rawvector)
	      {
		//K*  rawev(prob.RawEigenvectors());
		int m = Min(nconv,rawvector->M());
	    ffassert(rawvector->N()==n);
	    for(int i=0;i<m;i++)
	      {
		int k=i;
		KN_<K> vi(Z(':',i)) ;
		cout << " ------ EV--raw " << vi.min() << " " << vi.max() << endl;
		(*rawvector)(':',i)=vi;
	      }
	    
	      }
	
	    if (evector)
	      {
		// K*  rawev(prob.RawEigenvectors());
		// rawev + n*k is
		//  iev = prob.EigenvalueImag(k)
		//  iev==0 => the eigen vector 
		//  iev> 0 => real 
		//      start real :  rawev + n*k 
		//      start imag :  ramev +n*(k+1)
		//  iev < 0 =>  complex  
		//      start real :  rawev + n*(k-1) 
		//      -start imag :  ramev +n*(k)
		FEbaseArray<K,v_fes> & ev(*evector);
		int m = Min(nconv,(long) ev.N);
		for(int i=0;i<m;i++)
		  {
		    // K ev_i=
		    //prob.EigenvalueImag(i);
		    FEbase<K,v_fes> & xx= **(ev[i]);
		    // if (xx.pVh != pOP1->pUh) 
		    //    ExecError("Wrong Type of FEspace to store the eigen vector ");
		    // xx.Vh = pOP1->Uh;
		    // int  k=(ev_i < 0) ? i-1 : i;
		    //int k=i;
		    KN_<K> vi(Z(':',i));//rawev+n*k,n) ;
		    xx= new KN<K>(vi);
		    
		  }
	      }
	    
	  }
	} 
  
	return (long) nconv;
}
  

AnyType EigenValueC::E_EV::operator()(Stack stack)  const {
  
  double tol=1e-6;
  long nconv=0; 
  long nbev=1;
  long ncv =0;  // the number of Arnoldi vectors generated 
  long maxit=0;  // the maximum number of Arnoldi iterations 
  long mode=3;
  K sigma=0;
  KN<K> * evalue=0;
  KN<K> * resid=0;
  KNM<K> * rawvector=0;
  
  pfecarray  evector2;
  pfecbasearray   evector=0;
  tol=arg<double>(0,stack,0);
  nbev=arg<long>(1,stack,0);
  sigma=arg<K>(2,stack,0.0);
  evalue=arg<KN<K> *>(3,stack,0);
  evector2 =arg<pfecarray>(4,stack,make_pair<pfecbasearray,int>(0,0)); 
  ncv= arg<long>(5,stack,0);
  maxit= arg<long>(6,stack,0);
  rawvector=arg<KNM<K> *>(7,stack,0);
  resid=arg<KN<K> *>(8,stack,0);
  mode = arg<long>(9,stack,3);
  K * residptr= resid ? (K*) *resid : 0;
  evector=evector2.first;
  ffassert(mode>0 && mode <4) ; 
  Matrice_Creuse<K> *pOP1 =  GetAny<Matrice_Creuse<K> *>((*expOP1)(stack));
  Matrice_Creuse<K> *pB =  GetAny<Matrice_Creuse<K> *>((*expB)(stack));
  
  if(evalue) nbev=Max( (long)evalue->N(),nbev);
  if(!maxit)  maxit = 100*nbev;    
  
  const MatriceCreuse<K> & OP1 = pOP1->A;
  const MatriceCreuse<K> & B = pB->A;
  
  long  n=OP1.n;
  nbev=min(nbev,n-2);
  if(!ncv)     ncv = nbev*2+1;  
  ncv = max(nbev+2,ncv);
  ncv = min(ncv,n);
  const char *serr[10];
  int err=0;
  if(nbev>= n-1)
     serr[err++]="  Number of eigenvalues of OP to be computed <= n-2 ";
   if( mode < 1 || mode > 3) 
        serr[err++]="  the mode = 1 ,2 ,3  ";
  // 2 <= NCV-NEV and NCV <= N
   if( ! ( 2 <= nbev && ncv <= n)) 
       serr[err++]="   ( 2 <= nbve && nvc <= n) ";
  if (n != OP1.m) 
   serr[err++]=" the first matrix in EigneValue is not Hermitien.";
  if (n != B.n ) 
    serr[err++]="Sorry the row's number of the secand matrix in EigneValue is wrong.";
  if (n != B.m ) 
    serr[err++]="Sorry the colum's number of the secand matrix in EigneValue is wrong.";
    if(verbosity)
	cout << "Real complex eigenvalue problem: A*x - B*x*lambda" << endl;
  if(verbosity>9 ||  err)
	cout << "   n " << n << " nev "<< nbev << " tol =" << tol << " maxit =" << maxit << " ncv = " <<ncv << " mode = " << mode << endl;	
  if(err) 
    {
	cerr << " list of the error " << endl;
	for (int i=0;i<err;++i)
	    cerr << "\n-  " <<  serr[i] << endl;
	ExecError("Fatal Error in EigenValue (complex) ");
    }

      
  
 
  KN<K> work(n);
  //  ARrcSymGenEig is a class that requires the user to provide a
  //   way to perform the matrix-vector products w = OP*Bv =
  //   inv(A-sigma*B)*B*v and w = B*v, where sigma is the adopted shift.
  // OP1 = (A-sigma*B)
  // OP = inv(OP) 
  /*
  ffassert(0);
 
    
    ARrcCompGenEig<R> prob( n, nbev, sigma,"LM",ncv,tol,maxit,residptr);
    
    
    // OP = inv[A - sigma*I]
    
    
    // Finding an Arnoldi basis.
    
    while (!prob.ArnoldiBasisFound()) {
      
      // Calling ARPACK FORTRAN code. Almost all work needed to
      
      prob.TakeStep();
  */

  // cas non symetrique ,                                                                                                                    
  // Finding an Arnoldi basis.                                                                                                              \
                                                                                                                                                   
 // int mode=3; //  Shift invert			\
  
  int ido=0;
  char bmat='G';
  char which[]="LM";
  int ishift=1; // Auto Shift true by default                                                                                               \
                                                                                                                                                   
  int iparam[12]= {0,ishift,0,maxit,1,nconv,0,mode,0,0,0,0};
  int ipntr[15]={ 0,0,0, 0,0,0,  0,0,0, 0,0,0 ,0,0,0};
  KN<K> workd(3*n+1);
  int lworkl = 3*ncv*ncv+5*ncv;
  KN<K> workl(lworkl+1);
  KN<K> vp(ncv*n+1);
  int info= (residptr !=0);
  KN<double> rwork(ncv+1);
  KN<K> vresid(residptr? 1: n);
  if(!residptr) residptr=&vresid[0];

  while(1)
    {
      // GetVector supplies  a pointer to the input vector, v, 
      //  and PutVector a pointer  to the output vector, w.
      caupp(ido,bmat,n,which,nbev,
	    tol,  residptr,  ncv,
	    vp, n,   iparam,  ipntr,workd, workl,
	    lworkl,rwork,  info);
      
      if(ido==99) break;

      
      switch (ido) {
      case -1: {
	KN_<K> v(&workd[ipntr[1]],n);
	KN_<K> w(&workd[ipntr[2]],n);
	
	// Performing w <- OP*B*v for the first time.
	// This product must be performed only if GetIdo is equal to
	// -1. GetVector supplies a pointer to the input vector, v,
	// and PutVector a pointer to the output vector, w.
	work = B*v;
	if(mode==1)
	  w = OP1*v;
	else 	  
	  OP1.Solve(w,work);
	//    P.B.MultMv(&workd[ipntr[1]], temp);
	//    P.MultOPv(temp, &workd[ipntr[2]]);
	break;
      }
      case  1: {
	KN_<K> v(&workd[ipntr[3]],n);
	KN_<K> w(&workd[ipntr[2]],n);
	if(mode==1)
	  w = OP1*v;
	else 	  
	  OP1.Solve(w,v);
	break; }
	
      case  2: {
	KN_<K> v(&workd[ipntr[1]],n);
	KN_<K> w(&workd[ipntr[2]],n);
	w = B*v; 
	break;
      }
      case 3: {
	cout << " Bizzare " << 3 << endl;
      } 
      }
      sauppError(info);
    }
  nconv = iparam[5];
  if(nconv)
    {
      KN<K> evc(nbev);
      KNM<K> Z(n,nbev);
      KN<K> workev(3*ncv);
      int ldz=n;
      char HowMny ='A';
      int rvec=1;
      ceupp( rvec, HowMny, evc, Z , n, sigma, workev,	     bmat, n,   which, 
	     nbev,  tol,  residptr, ncv,
	     vp, n, iparam , ipntr, workd, workl,lworkl,rwork,info);
      
      if (verbosity)
	{
	  cout << "Complex eigenvalue problem: A*x - B*x*lambda" << endl;
	  cout << "mode =" << mode << "  sigma=" << sigma <<  endl << endl;
	  
	  cout << "Dimension of the system            : " << n              << endl;
	  cout << "Number of 'requested' eigenvalues  : " << nbev  << endl;
	  cout << "Number of 'converged' eigenvalues  : " << nconv          << endl;
	  cout << "Number of Arnoldi vectors generated: " << ncv << endl;
	  cout << "Number of iterations taken         : " << iparam[3]  << endl;
	  cout << endl;
	  
	  cout << "Eigenvalues:" << endl;
	  for (int i=0; i<nconv; i++) {
	    cout << "  lambda[" << (i+1) << "]: " << evc(i) << endl;
	    KN_<K> vi(Z(':',i)) ;
	    if(verbosity>99)
	      cout <<" Eigenvector: :" << vi << endl;
	  }
		       cout << endl;
	}
      
      
      if (evalue)
	{
	  KN<K> & ev(*evalue);
	  int m = Min(nconv,ev.N());
	  for(int i=0;i<m;i++)
	    ev[i]=evc(i);
	}
      if (evector)
	{
	  FEbaseArray<K,v_fes> & ev(*evector);
	  int m = Min(nconv,(long) ev.N);
	  for(int i=0;i<m;i++)
	    {
	      FEbase<K,v_fes> & xx= **(ev[i]);
	      KN_<K> vi(Z(':',i)) ;
	      xx= new KN<K>(vi);
	      
	    }
	  if(rawvector)
	    {
	      int m = Min(nconv,rawvector->M());
	      ffassert(rawvector->N()==n);
	      for(int i=0;i<m;i++)
		{
		  KN_<K> vi(Z(':',i)) ;
		  (*rawvector)(':',i)=vi;
		}
	      
	    }
	  
	}
    }
  
return (long) nconv;
}

#ifndef DYNM_LOAD
void init_eigenvalue()
{
    if(verbosity) cout << "eigenvalue ";
    Global.Add("EigenValue","(",new EigenValue(1));  //  j + dJ
    Global.Add("EigenValue","(",new EigenValueC(1));  //  j + dJ
    
}
#else
class Init {
public:
  Init()
  {
    if(verbosity) cout << "eigenvalue ";
    Global.Add("EigenValue2","(",new EigenValue(1));  //  j + dJ
    Global.Add("EigenValue2","(",new EigenValueC(1));  //  j + dJ
  }
};
Init init; 
#endif
