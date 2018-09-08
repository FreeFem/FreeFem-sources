#include <iostream>
#include <cmath>
#include "HashMatrix.hpp"
#include <vector>
#include "VirtualSolver.hpp"

#include <complex>

template<class Z=int,class K=double>
class VirtualSolverSkyLine: public VirtualSolver<Z,K> {
public:
    typedef HashMatrix<Z,K>  HMat;
    HMat *A;
    
    Z *Ai,*Ap;
    Z *numi,*numj; //
    K *Ax;
    
    mutable int status;

    VirtualSolverSkyLine(HMat  *AA):A(AA) {}
    void dosolver(K *x,K*b,int N,int trans) {assert(0);}
    void fac_symbolic(){assert(0);}
    void fac_numeric(){assert(0);}
    ~VirtualSolverSkyLine(){}
    void SetState(){}
};

// specilisation
template<>
class  VirtualSolverSkyLine<int,double> : public VirtualSolver<int,double> {
public:
    typedef double K;
    typedef int Z;
    typedef HashMatrix<Z,K>  HMat;
    HMat *A;
    void *Symbolic, *Numeric ;
    Z *Ai,*Ap;
    K *Ax;
    int  cs,cn;
    VirtualSolverSkyLine(HMat  &AA,double tol_pivot)
    :A(&AA),Symbolic(0),Numeric(0),Ai(0),Ap(0),Ax(0),cs(0),cn(0)
    {
 
    }
    void dosolver(double *x,double*b,int N,int trans) {
        for(int k=0,oo=0; k<N;++k, oo+= A->n)
        {

        }

    }
    
    void SetState(){
        if( A->re_do_numerics ) cn++;
        if( A->re_do_symbolic) cs++;
        CheckState(A->n,cs,cn);
    }
                           
    void fac_symbolic(){
        A->CSC(Ap,Ai,Ax);
        
      }
    void fac_numeric(){
        if(Numeric)   umfpack_di_free_numeric (&Numeric) ;
       

    }
    ~VirtualSolverSkyLine()
    {
        if(Symbolic)  ;
         if(Numeric)   ;
    }
};

// specilisation
template<>
class  VirtualSolverSkyLine<int,std::complex<double> > : public VirtualSolver<int,std::complex<double> > {
public:
    typedef std::complex<double> K;
    typedef int Z;
    typedef HashMatrix<Z,K>  HMat;
    HMat *A;
    void *Symbolic, *Numeric ;
    Z *Ai,*Ap;
    K *Ac;
    double *Ax,*Az;
    int  cs,cn;
    mutable int status;
    
    VirtualSolverSkyLine(HMat  &AA,
                         double tol_pivot=-1.)
    :A(&AA)
    {
  
    }
    void dosolver(K *x,K*b,int N,int trans) {
        for(int k=0,oo=0; k<N;++k, oo+= A->n)
        {
            double * xx = (double *) (void*) x+oo,  *bb = (double *) (void*) b+oo, *zx=0;;
        }
        
    }
    
    void SetState(){
        if( A->re_do_numerics ) cn++;
        if( A->re_do_symbolic) cs++;
        CheckState(A->n,cs,cn);
    }
    
    void fac_symbolic(){
        A->CSC(Ap,Ai,Ac);
        Ax= (double *) (void *) Ac;
        Az=0;
        
        
     
        if(status) cout << " Error umpfack umfpack_zi_symbolic  status  " << status << endl;
    }
    void fac_numeric(){
        
        if(status) cout << " Error umpfack umfpack_zi_numeric  status  " << status << endl;
        
    }
    ~VirtualSolverSkyLine()
    {

    }
};


template<class Z=int,class K=double>
class VirtualSolverCHOLMOD: public VirtualSolver<Z,K> {
public:
    typedef HashMatrix<Z,K>  HMat;
    HMat *HA;
    cholmod_common Common, *cm ;
    cholmod_factor *L ;
    cholmod_sparse *A ;
    cholmod_dense *X = NULL, *B, *W, *R ;

    mutable int status;
    
    VirtualSolverCHOLMOD(HMat  *AA):A(AA) {}
    void dosolver(K *x,K*b,int N,int trans) {assert(0);}
    void fac_symbolic(){assert(0);}
    void fac_numeric(){assert(0);}
    ~VirtualSolverCHOLMOD(){}
    void SetState(){}
};

// specilisation
template<>
class  VirtualSolverCHOLMOD<int,double> : public VirtualSolver<int,double> {
public:
    typedef double K;
    typedef int Z;
    typedef HashMatrix<Z,K>  HMat;
    HMat *HA;
    int n;
    const int xtype=CHOLMOD_REAL;
    cholmod_common c ;
    cholmod_factor *L ;
    cholmod_sparse AA,*A ;
    Z *Ai,*Ap;
    K *Ax;
    cholmod_dense *Ywork, *Ework  ;
    int  cs,cn;
    
 
   void set_cholmod_dense(cholmod_dense & X,K *p,int m)
    {
        X.nrow=n;
        X.ncol=m;
        X.nzmax=n*m;
        X.d=n;
        X.x = p;
        X.z=0;
        X.xtype=CHOLMOD_REAL;
        X.dtype=CHOLMOD_DOUBLE;
    }
    
    mutable int status;
    VirtualSolverCHOLMOD(HMat  &HAA)
      :HA(&HAA),n(HAA.n),L(0),A(&AA),Ai(0),Ap(0),Ax(0),Ywork(0),Ework(0),cs(0),cn(0)
    {
       
        cholmod_start (&c) ;
        //CHOLMOD_FUNCTION_DEFAULTS (&c) ;
        AA.nrow=n;
        AA.ncol=n;
        AA.nzmax=HA->nnz;
        AA.p=0;
        AA.i=0;
        AA.nz=0;
        AA.x =0;
        AA.z=0;
        AA.stype=1;// U
        AA.itype=CHOLMOD_INT;
        AA.xtype=CHOLMOD_REAL;
        AA.dtype=CHOLMOD_DOUBLE;
        AA.sorted=1;
        AA.packed=1;
        
        //c.error_handler = my_handler ;

    }
    
    void dosolver(K *x,K*b,int N,int trans)
    {
        cholmod_dense XX,*X=&XX,B;
        set_cholmod_dense(*X,x,N);
        set_cholmod_dense(B,b,N);
        
        cout << " dosolver CHOLMoD double "<< endl;
       cholmod_solve2 (CHOLMOD_A, L, &B, NULL, &X, NULL,
                        &Ywork, &Ework, &c) ;
        if( X !=  &XX) cholmod_free_dense (&X, &c) ;
        
    }
    
    void SetState(){
        if( HA->re_do_numerics ) cn++;
        if( HA->re_do_symbolic) cs++;
        CheckState(HA->n,cs,cn);
    }
    
    void fac_symbolic()
    {
        AA.nzmax=HA->CSC_U(Ap,Ai,Ax);
        cout << "fac_symbolic cholmod R: nnz U=" << AA.nzmax << " nnz= "  << HA->nnz << endl;

        AA.p=Ap;
        AA.i=Ai;
        AA.x=Ax;
        if(L) cholmod_free_factor (&L, &c) ;		    /* free matrices */
        L = cholmod_analyze (A, &c) ;
     }
    void fac_numeric(){
        cout << " fac_numeric CHOLMoD double "<< endl;

        cholmod_factorize (A, L, &c) ;
      
    }
    ~VirtualSolverCHOLMOD()
    {
        if(L) cholmod_free_factor (&L, &c) ;		    /* free matrices */
 //w       if(A) cholmod_free_sparse (&A, &c) ;
        cholmod_finish (&c) ;

    }
};

// specialisation
template<>
class  VirtualSolverCHOLMOD<int,std::complex<double> > : public VirtualSolver<int,std::complex<double> >
{
public:
    typedef std::complex<double>  K;
    typedef int Z;
    typedef HashMatrix<Z,K>  HMat;
    HMat *HA;
    int n;
    static const int xtype=CHOLMOD_COMPLEX ;// a complex matrix (ANSI C99 compatible)
    static const int dtype =CHOLMOD_DOUBLE;
    cholmod_common c ;
    cholmod_factor *L ;
    cholmod_sparse AA,*A ;
    Z *Ai,*Ap;
    K *Ax;
    cholmod_dense *Ywork, *Ework  ;
    int  cs,cn;
    
    
    void set_cholmod_dense(cholmod_dense & X,K *p,int m)
    {
        X.nrow=n;
        X.ncol=m;
        X.nzmax=n*m;
        X.d=n;
        X.x = p;
        X.z=0;
        X.xtype=xtype;// a complex matrix (ANSI C99 compatible)
        X.dtype=dtype;
    }
    
    VirtualSolverCHOLMOD(HMat  &HAA)
    :HA(&HAA),n(HAA.n),L(0),A(&AA),Ai(0),Ap(0),Ax(0),Ywork(0),Ework(0),cs(0),cn(0)
    {
        
        cholmod_start (&c) ;
       // CHOLMOD_FUNCTION_DEFAULTS (&c) ;
        AA.nrow=n;
        AA.ncol=n;
        AA.nzmax=HA->nnz;
        AA.p=0;
        AA.i=0;
        AA.nz=0;
        AA.x =0;
        AA.z=0;
        AA.stype=1;// U
        AA.itype=CHOLMOD_INT;
        AA.xtype=xtype;
        AA.dtype=dtype;
        AA.sorted=1;
        AA.packed=1;
        
        //c.error_handler = my_handler ;
        
    }
    
    void dosolver(K *x,K*b,int N,int trans)
    {
        cout << " dosolver CHOLMoD Complex "<< endl;
        cholmod_dense XX,*X=&XX,B;
        set_cholmod_dense(*X,x,N);
        set_cholmod_dense(B,b,N);
        
        cholmod_solve2 (CHOLMOD_A, L, &B, NULL, &X, NULL,
                        &Ywork, &Ework, &c) ;
        if( X !=  &XX) cholmod_free_dense (&X, &c) ;
        
    }
    
    void SetState(){
        if( HA->re_do_numerics ) cn++;
        if( HA->re_do_symbolic) cs++;
        CheckState(HA->n,cs,cn);
    }
    
    void fac_symbolic()
    {
        AA.nzmax=HA->CSC_U(Ap,Ai,Ax);
        cout << "fac_symbolic cholmod C: nnz U=" << AA.nzmax << " nnz= "  << HA->nnz << endl;
 //       HA->CSC(Ap,Ai,Ax);
        AA.p=Ap;
        AA.i=Ai;
        AA.x=Ax;
        if(L) cholmod_free_factor (&L, &c) ;		    /* free matrices */
        L = cholmod_analyze (A, &c) ;
    }
    void fac_numeric(){
        cout << " fac_numeric CHOLMoD complex "<< endl;
        cholmod_factorize (A, L, &c) ;
        
    }
    ~VirtualSolverCHOLMOD()
    {
        if(L) cholmod_free_factor (&L, &c) ;		    /* free matrices */
        //w       if(A) cholmod_free_sparse (&A, &c) ;
        cholmod_finish (&c) ;
        
    }
};


