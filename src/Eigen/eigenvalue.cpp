#include <fstream>
#include <iostream>
#include <cassert>
#include <algorithm>

//#include <ngenprbc.h>
#include <lapackc.h>
#include <arrgsym.h>
#include <arrgnsym.h>
using namespace std;
#include "error.hpp"
#include "AFunction.hpp"
#include "rgraph.hpp"
#include "RNM.hpp"
#include "MatriceCreuse_tpl.hpp"
#include "MeshPoint.hpp"
#include "lgfem.hpp"
#include "lgsolver.hpp"
#include "problem.hpp"
extern Block *currentblock;

typedef double R;


class EigenValue : public OneOperator
{ public:
    typedef R K;
    typedef KN<K> Kn;
    typedef KN_<K> Kn_;
     const int cas;
    class E_EV: public E_F0mps { public:
        const int cas;
        
        static basicAC_F0::name_and_type name_param[] ;
        static const int n_name_param =9;
        Expression nargs[n_name_param];
        Expression expOP1,expB;
        template<class T>
        T arg(int i,Stack stack,const T & a) const{ return nargs[i] ? GetAny<T>( (*nargs[i])(stack) ): a;}
        E_EV(const basicAC_F0 & args,int cc) :
             cas(cc)
        {
              // OP1 = (A-sigma*B)        
                int nbj= args.size()-1;
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

basicAC_F0::name_and_type  EigenValue::E_EV::name_param[]= {
    "tol", &typeid(double)  ,
    "nbev",&typeid(long) ,
    "sym",&typeid(bool),
    "sigma",&typeid(double),
    "value",&typeid(KN<double> *),
    "vector",&typeid(pferarray) , 
    "ncv",&typeid(long) , // the number of Arnoldi vectors generated 
    "maxit",&typeid(long), // the maximum number of Arnoldi iterations 
    "ivalue",&typeid(KN<double> *),
};

        
 AnyType EigenValue::E_EV::operator()(Stack stack)  const {
           double tol=1e-6;
           long nconv=0; 
           long nbev=0;
           bool sym=false;
           long ncv =0;  // the number of Arnoldi vectors generated 
           long maxit=0;  // the maximum number of Arnoldi iterations 
           double sigma=0;
           KN<double> * evalue=0;
           KN<double> * evaluei=0;
           
           pferarray  evector2,eivector2;
           pferbasearray   evector=0,eivector=0;;
           tol=arg<double>(0,stack,0);
           nbev=arg<long>(1,stack,0);
           sym=arg<bool>(2,stack,false);
           sigma=arg<double>(3,stack,0.0);
           evalue=arg<KN<double> *>(4,stack,0);
           evector2 =arg<pferarray>(5,stack,make_pair<pferbasearray,int>(0,0)); 
           ncv= arg<long>(6,stack,0);
           maxit= arg<long>(7,stack,0);
           evaluei=arg<KN<double> *>(8,stack,0);
           eivector2 =arg<pferarray>(9,stack,make_pair<pferbasearray,int>(0,0)); 

           evector=evector2.first;
           eivector=eivector2.first;

           Matrice_Creuse<K> *pOP1 =  GetAny<Matrice_Creuse<K> *>((*expOP1)(stack));
           Matrice_Creuse<K> *pB =  GetAny<Matrice_Creuse<K> *>((*expB)(stack));
             
           if(evalue) nbev=Max( (long)evalue->N(),nbev);
           
           
           const MatriceCreuse<K> & OP1 = pOP1->A;
           const MatriceCreuse<K> & B = pB->A;
           
           int n=OP1.n;
           if (n != OP1.m) 
              ExecError("Sorry the first matrix in EigneValue is not symetric.");
           if (n != B.n ) 
              ExecError("Sorry the row's number of the secand matrix in EigneValue is wrong.");
           if (n != B.m ) 
              ExecError("Sorry the colum's number of the secand matrix in EigneValue is wrong.");
           if(verbosity)
             if(sym)
	       cout << "Real symmetric eigenvalue problem: A*x - B*x*lambda" << endl;
	     else
	       cout << "Real non symmetric eigenvalue problem: A*x - B*x*lambda" << endl;


	   KN<K> work(n);
	   //  ARrcSymGenEig is a class that requires the user to provide a
	   //   way to perform the matrix-vector products w = OP*Bv =
	   //   inv(A-sigma*B)*B*v and w = B*v, where sigma is the adopted shift.
	   // OP1 = (A-sigma*B)
	   // OP = inv(OP) 
	   // cas symetrique
	   try {
	     if(sym)
	       {
		 ARrcSymGenEig<K> prob('S', n, nbev, sigma,"LM",ncv,tol,maxit);
		 
		 // ARrcNonSymGenEig<K> prob(n, nbev, sigma);
		 
		 // OP = inv[A - sigma*I]
		 
		 
		 // Finding an Arnoldi basis.
		 
		 while (!prob.ArnoldiBasisFound()) {
		   
		   // Calling ARPACK FORTRAN code. Almost all work needed to
		   // find an Arnoldi basis is performed by TakeStep.
		   
		   prob.TakeStep();
		   
		   // GetVector supplies  a pointer to the input vector, v, 
		   //  and PutVector a pointer  to the output vector, w.
		   int kkk;
		   switch (kkk=prob.GetIdo()) {
		   case -1: {
		     KN_<K> v(prob.GetVector(),n);
		     KN_<K> w(prob.PutVector(),n);
		     
		     // Performing w <- OP*B*v for the first time.
		     // This product must be performed only if GetIdo is equal to
		     // -1. GetVector supplies a pointer to the input vector, v,
		     // and PutVector a pointer to the output vector, w.
		     work = B*v;
		     OP1.Solve(w,work);
		     //    P.B.MultMv(prob.GetVector(), temp);
		     //    P.MultOPv(temp, prob.PutVector());
		     break;
		   }
		   case  1: {
		     KN_<K> v(prob.GetProd(),n);
		     KN_<K> w(prob.PutVector(),n);
		     
		     // Performing w <- OP*B*v when Bv is available.
		     // This product must be performed whenever GetIdo is equal to
		     // 1. GetProd supplies a pointer to the previously calculated
		     // product Bv and PutVector a pointer to the output vector w.
		     OP1.Solve(w,v);
		     //P.MultOPv(prob.GetProd(), prob.PutVector());
		     break; }
		     
		   case  2: {
		     KN_<K> v(prob.GetVector(),n);
		     KN_<K> w(prob.PutVector(),n);
		     
		     // Performing w <- B*v.
		     //P.B.MultMv(prob.GetVector(), prob.PutVector());
		     w = B*v; 
		   }
		   }
		   // cout<< " GetIdo = " << kkk << endl;
		 }
		 
		 // Finding eigenvalues and eigenvectors.
		 
		 prob.FindEigenvectors();
		 
		 // Printing solution.
		 int    mode;
		 
		 nconv = prob.ConvergedEigenvalues();
		 mode  = prob.GetMode();
		 
		 if (verbosity)
		   {
		     cout << endl << endl << "Thanks to ARPACK++ class ARrcSymGenEig" << endl;
		     cout << "Real symmetric eigenvalue problem: A*x - B*x*lambda" << endl;
		     switch (mode) {
		     case 2:
		       cout << "Regular mode" << endl << endl;
		       break;
		     case 3:
		       cout << "Shift and invert mode  sigma=" << sigma <<  endl << endl;
		       break;
		     case 4:
		       cout << "Buckling mode" << endl << endl;
		       break;
		     case 5:
		       cout << "Cayley mode" << endl << endl;
		     }
		     
		     cout << "Dimension of the system            : " << n              << endl;
		     cout << "Number of 'requested' eigenvalues  : " << prob.GetNev()  << endl;
		     cout << "Number of 'converged' eigenvalues  : " << nconv          << endl;
		     cout << "Number of Arnoldi vectors generated: " << prob.GetNcv()  << endl;
		     cout << "Number of iterations taken         : " << prob.GetIter() << endl;
		     cout << endl;
		     
		     if (prob.EigenvaluesFound()) {
		       cout << "Eigenvalues:" << endl;
		       for (int i=0; i<nconv; i++) {
			 cout << "  lambda[" << (i+1) << "]: " << prob.Eigenvalue(i) << endl;
			 KN_<K> vi(prob.RawEigenvector(i),n) ;
			 if(verbosity>99)
			   cout <<" Eigenvector: :" << vi << endl;
		       }
		       cout << endl;
		     }
		   } 
		 
		 if (! prob.EigenvaluesFound()) 
		   nconv=0;  //  bizarre 
		 if (evalue)
		   {
		     KN<double> & ev(*evalue);
		     int m = Min(nconv,ev.N());
		     for(int i=0;i<m;i++)
		       ev[i]=prob.Eigenvalue(i);
		   }
		 if (evector)
		   {
		     FEbaseArray<K> & ev(*evector);
		     int m = Min(nconv,(long) ev.N);
		     for(int i=0;i<m;i++)
		       {
			 FEbase<K> & xx= **(ev[i]);
			 // if (xx.pVh != pOP1->pUh) 
			 //    ExecError("Wrong Type of FEspace to store the eigen vector ");
			 xx.Vh = pOP1->Uh;
			 KN_<K> vi(prob.RawEigenvector(i),n) ;
			 xx= new KN<K>(vi);
			 
		       }
		   }
		 
	       }
	     else 
	       {  // cas non symetric ,
		 
		 ARrcNonSymGenEig<K> prob( n, nbev, sigma,"LM",ncv,tol,maxit);
		 
		 
		 // Finding an Arnoldi basis.
		 
		 while (!prob.ArnoldiBasisFound()) {
		   
		   // Calling ARPACK FORTRAN code. Almost all work needed to
		   // find an Arnoldi basis is performed by TakeStep.
		   
		   prob.TakeStep();
		   
		   // GetVector supplies  a pointer to the input vector, v, 
		   //  and PutVector a pointer  to the output vector, w.
		   int kkk;
		   switch (kkk=prob.GetIdo()) {
		   case -1: {
		     KN_<K> v(prob.GetVector(),n);
		     KN_<K> w(prob.PutVector(),n);
		     
		     // Performing w <- OP*B*v for the first time.
		     // This product must be performed only if GetIdo is equal to
		     // -1. GetVector supplies a pointer to the input vector, v,
		     // and PutVector a pointer to the output vector, w.
		     work = B*v;
		     OP1.Solve(w,work);
		     //    P.B.MultMv(prob.GetVector(), temp);
		     //    P.MultOPv(temp, prob.PutVector());
		     break;
		   }
		   case  1: {
		     KN_<K> v(prob.GetProd(),n);
		     KN_<K> w(prob.PutVector(),n);
		     
		     // Performing w <- OP*B*v when Bv is available.
		     // This product must be performed whenever GetIdo is equal to
		     // 1. GetProd supplies a pointer to the previously calculated
		     // product Bv and PutVector a pointer to the output vector w.
		     OP1.Solve(w,v);
		     //P.MultOPv(prob.GetProd(), prob.PutVector());
		     break; }
		     
		   case  2: {
		     KN_<K> v(prob.GetVector(),n);
		     KN_<K> w(prob.PutVector(),n);
		     
		     // Performing w <- B*v.
		     //P.B.MultMv(prob.GetVector(), prob.PutVector());
		     w = B*v; 
		   }
		   }
		   // cout<< " GetIdo = " << kkk << endl;
		 }
		 
		 // Finding eigenvalues and eigenvectors.
		 
		 nconv =prob.FindEigenvectors();
		 
		 //	 RecoverEigenvalues(nconv,n, OP1,sigma,B, prob.RawEigenvectors(),
		 //		    prob.RawEigenvalues(), prob.RawEigenvaluesImag());
		 
		 
		 // Printing solution.
		 int    mode;
		 
		 nconv = prob.ConvergedEigenvalues();
		 mode  = prob.GetMode();
		 
		 if (verbosity)
		   {
		     cout << endl << endl << "Thanks to ARPACK++ class ARrcSymGenEig" << endl;
		     cout << "Real non-symmetric eigenvalue problem: A*x - B*x*lambda" << endl;
		     switch (mode) {
		     case 2:
		       cout << "Regular mode" << endl << endl;
		       break;
		     case 3:
		       cout << "Shift and invert mode  sigma=" << sigma <<  endl << endl;
		       break;
		     case 4:
		       cout << "Buckling mode" << endl << endl;
		       break;
		     case 5:
		       cout << "Cayley mode" << endl << endl;
		     }
		     
		     cout << "Dimension of the system            : " << n              << endl;
		     cout << "Number of 'requested' eigenvalues  : " << prob.GetNev()  << endl;
		     cout << "Number of 'converged' eigenvalues  : " << nconv          << endl;
		     cout << "Number of Arnoldi vectors generated: " << prob.GetNcv()  << endl;
		     cout << "Number of iterations taken         : " << prob.GetIter() << endl;
		     cout << endl;
		     
		     if (prob.EigenvaluesFound()) {
		       cout << "Eigenvalues:" << endl;
		       for (int i=0; i<nconv; i++) {
			 cout << "  lambda[" << (i+1) << "]: " ;
			 ios::fmtflags oldflag= cout.flags();
			 cout.setf(ios::showpos);
			 cout << prob.EigenvalueReal(i) 
			      <<  prob.EigenvalueImag(i) << "i"<<  endl;
			 cout.flags(oldflag); // restore flags     
			 KN_<K> vi(prob.RawEigenvector(i),n) ;
			 if(verbosity>99)
			   cout <<" Eigenvector: :" << vi << endl;
		       }
		       cout << endl;
		     }
		   } 
		 
		 if (! prob.EigenvaluesFound()) 
		   nconv=0;  //  bizarre 
		 if (evalue)
		   {
		     KN<double> & ev(*evalue);
		     int m = Min(nconv,ev.N());
		     for(int i=0;i<m;i++)
		       ev[i]=prob.EigenvalueReal(i);
		   }
		 
		 if (evaluei)
		   {
		     KN<double> & ev(*evaluei);
		     int m = Min(nconv,ev.N());
		     for(int i=0;i<m;i++)
		       ev[i]=prob.EigenvalueImag(i);
		   }
		 
		 if (evector)
		   {
              K*  rawev(prob.RawEigenvectors());
			       // rawev + n*k is
			       //  iev = prob.EigenvalueImag(k)
			       //  iev==0 => the eigen vector 
			       //  iev> 0 => real 
			       //      start real :  rawev + n*k 
			       //      start imag :  ramev +n*(k+1)
			       //  iev < 0 =>  complex  
			       //      start real :  rawev + n*(k-1) 
			       //      -start imag :  ramev +n*(k)
		     FEbaseArray<K> & ev(*evector);
		     int m = Min(nconv,(long) ev.N);
		     for(int i=0;i<m;i++)
		       {
			 K ev_i=prob.EigenvalueImag(i);
			 FEbase<K> & xx= **(ev[i]);
			 // if (xx.pVh != pOP1->pUh) 
			 //    ExecError("Wrong Type of FEspace to store the eigen vector ");
			 xx.Vh = pOP1->Uh;
            // int  k=(ev_i < 0) ? i-1 : i;
            int k=i;
			 KN_<K> vi(rawev+n*k,n) ;
			 xx= new KN<K>(vi);
			 
		       }
		   }
/*		   
		 if (eivector)
		   {
             K*  rawev(prob.RawEigenvectors());
		     FEbaseArray<K> & ev(*eivector);
		     int m = Min(nconv,(long) ev.N);
		     for(int i=0;i<m;i++)
		       {
			 K ev_i=prob.EigenvalueImag(i);
			 FEbase<K> & xx= **(ev[i]);
			 // if (xx.pVh != pOP1->pUh) 
			 //    ExecError("Wrong Type of FEspace to store the eigen vector ");
			 xx.Vh = pOP1->Uh;
			 int k= ev_i <=0 ? i : i+1;
			 KN_<K> vi(rawev+n*k,n) ;
			 xx= new KN<K>(vi);
			 if (ev_i<0) *xx.x() *= -1.;
			 if (!ev_i)  *xx.x() =0; 
		       }
		   }
*/		 
		 
		 
		 
	       }
	   } catch (ArpackError) 
	     {
	       ExecError("ArpackError");
	     }
	   
	   return (long) nconv;
 }


void init_eigenvalue()
{
  cout << "eigenvalue ";
   Global.Add("EigenValue","(",new EigenValue(1));  //  j + dJ
    
}
