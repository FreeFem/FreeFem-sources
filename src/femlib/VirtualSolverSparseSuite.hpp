#ifndef VIRTUALSOLVERSPARSESUITE_HPP_
#define VIRTUALSOLVERSPARSESUITE_HPP_

#include <iostream>
#include <cmath>
#include "HashMatrix.hpp"
#ifdef HAVE_LIBUMFPACK
extern "C" {
#ifdef HAVE_UMFPACK_H
#include <umfpack.h>
#include <cholmod.h>
#else
#ifdef HAVE_UMFPACK_UMFPACK_H
#include <umfpack/umfpack.h>
#else
#ifdef HAVE_BIG_UMFPACK_UMFPACK_H
#include <UMFPACK/umfpack.h>
#else
#ifdef HAVE_UFSPARSE_UMFPACK_H
#include <ufsparse/umfpack.h>
#else
#ifdef HAVE_SUITESPARSE_UMFPACK_H
#include <suitesparse/umfpack.h>
#include <suitesparse/cholmod.h>
#else

    // Defaults to a local version of the UMFPACK headers
#include "umfpack.h"
#include "cholmod.h"

#endif // HAVE_SUITESPARSE_UMFPACK_H
#endif // HAVE_UFSPARSE_UMFPACK_H
#endif // HAVE_BIG_UMFPACK_UMFPACK_H
#endif // HAVE_UMFPACK_UMFPACK_H
#endif // HAVE_UMFPACK_H
}

#include <vector>
#include "VirtualSolver.hpp"
#include <complex>


template<class Z=int,class K=double>
class VirtualSolverUMFPACK: public VirtualSolver<Z,K> {
public:
    //  1 unsym , 2 sym, 4 pos , 8 nopos, 16  seq, 32  ompi, 64 mpi ,
    static const int orTypeSol = 1|8|16;
    typedef HashMatrix<Z,K>  HMat;
    HMat *A;
    void *Symbolic, *Numeric ;
    Z *Ai,*Ap;
    K *Ax;
    long verb;
    mutable int status;

    VirtualSolverUMFPACK(HMat  *AA):A(AA) {}
    void dosolver(K *x,K*b,int N,int trans) {assert(0);}
    void fac_symbolic(){assert(0);}
    void fac_numeric(){assert(0);}
    ~VirtualSolverUMFPACK(){}
    void UpdateState(){}
};

// specilisation
template<>
class  VirtualSolverUMFPACK<int,double> : public VirtualSolver<int,double> {
public:
    //  1 unsym , 2 sym, 4 pos , 8 nopos, 16  seq, 32  ompi, 64 mpi ,
    static const int orTypeSol = 1|8|16;

    typedef double K;
    typedef int Z;
    typedef HashMatrix<Z,K>  HMat;
    HMat *A;
    void *Symbolic, *Numeric ;
    Z *Ai,*Ap;
    K *Ax;
    int  cs,cn;
    long verb;
    mutable int status;
    double Control[UMFPACK_CONTROL];
    double Info[UMFPACK_INFO];
    VirtualSolverUMFPACK(HMat  &AA, const Data_Sparse_Solver & ds,Stack stack )
                        // int strategy=-1,
                        // double tol_pivot=-1.,double tol_pivot_sym=-1.,long vverb=verbosity)
    :A(&AA),Symbolic(0),Numeric(0),Ai(0),Ap(0),Ax(0),cs(0),cn(0), verb(ds.verb)
    {
         if(verb>2 || verbosity> 9) cout << " build solver UMFPACK double/int " << endl;
        for(int i=0;i<UMFPACK_CONTROL;i++) Control[i]=0;
        for(int i=0;i<UMFPACK_INFO;i++) Info[i]=0;

        umfpack_di_defaults (Control) ;

        if(ds.verb>4) Control[UMFPACK_PRL]=2;
        if(ds.tol_pivot_sym>0) Control[UMFPACK_SYM_PIVOT_TOLERANCE]=ds.tol_pivot_sym;
        if(ds.tol_pivot>0) Control[UMFPACK_PIVOT_TOLERANCE]=ds.tol_pivot;
        if(ds.strategy>=0)   Control[UMFPACK_STRATEGY]=ds.strategy;

    }
    void dosolver(double *x,double*b,int N,int trans) {
        if(verb>2 || verbosity> 9) cout << "dosolver UMFPACK double/int  "<< N << " " << trans << endl;
        int TS=trans ? UMFPACK_At: UMFPACK_A;
        for(int k=0,oo=0; k<N;++k, oo+= A->n)
        {

        status= umfpack_di_solve (TS, Ap, Ai, Ax, x+oo, b+oo, Numeric,Control,Info) ;
        if(status) cout << " Error umfpack_di_solve  status  " << status << endl;
        if(verbosity>3)     (void)  umfpack_di_report_info(Control,Info);
        }

    }

    void UpdateState(){
         if(verb>2 || verbosity> 9) std::cout << " UpdateState "<< A-> re_do_numerics << " " << A-> re_do_symbolic <<std::endl;
        if( A->GetReDoNumerics() ) cn++;
        if( A->GetReDoSymbolic()) cs++;
        ChangeCodeState(A->n,cs,cn);
    }

    void fac_symbolic(){
         if(verb>2 || verbosity> 9) cout << "fac_symbolic UMFPACK R: nnz U " << " nnz= "  << A->nnz << endl;

        A->CSC(Ap,Ai,Ax);
        if(Symbolic)  umfpack_di_free_symbolic (&Symbolic) ;
        status = umfpack_di_symbolic (A->n, A->m, Ap, Ai, Ax, &Symbolic,Control,Info) ;
        if(status) cout << " Error umpfack umfpack_di_symbolic  status  " << status << endl;
      }
    void fac_numeric(){
        if(Numeric)   umfpack_di_free_numeric (&Numeric) ;
         if(verb>2 || verbosity> 9) cout << "fac_numeric UMFPACK R: nnz U " << " nnz= "  << A->nnz << endl;

        status = umfpack_di_numeric (Ap, Ai, Ax, Symbolic, &Numeric, Control,Info) ;
        if(status) cout << " Error umpfack umfpack_di_numeric  status  " << status << endl;

    }
    ~VirtualSolverUMFPACK()
    {
     if(Symbolic)  umfpack_di_free_symbolic (&Symbolic) ;
     if(Numeric)   umfpack_di_free_numeric (&Numeric) ;
    }
};

// specilisation
template<>
class  VirtualSolverUMFPACK<int,std::complex<double> > : public VirtualSolver<int,std::complex<double> > {
public:
    //  1 unsym , 2 sym, 4 pos , 8 nopos, 16  seq, 32  ompi, 64 mpi ,
    static const int orTypeSol = 1|8|16;

    typedef std::complex<double> K;
    typedef int Z;
    typedef HashMatrix<Z,K>  HMat;
    HMat *A;
    void *Symbolic, *Numeric ;
    Z *Ai,*Ap;
    K *Ac;
    double *Ax,*Az;
    int  cs,cn;
    long verb;
    mutable int status;
    double Control[UMFPACK_CONTROL];
    double Info[UMFPACK_INFO];

    VirtualSolverUMFPACK(HMat  &AA,  const Data_Sparse_Solver & ds,Stack stack )//int strategy=-1,
                       //  double tol_pivot=-1.,double tol_pivot_sym=-1., long vverb=verbosity)
    :A(&AA),Symbolic(0),Numeric(0),Ai(0),Ap(0),Ax(0),cs(0),cn(0),verb(ds.verb)
    {
        if(verb>2 || verbosity> 9) cout << " build solver UMFPACK complex/int " << endl;

        for(int i=0;i<UMFPACK_CONTROL;i++) Control[i]=0;
        for(int i=0;i<UMFPACK_INFO;i++) Info[i]=0;

        umfpack_zi_defaults (Control) ;

        if(ds.verb>4) Control[UMFPACK_PRL]=2;
        if(ds.tol_pivot_sym>0) Control[UMFPACK_SYM_PIVOT_TOLERANCE]=ds.tol_pivot_sym;
        if(ds.tol_pivot>0) Control[UMFPACK_PIVOT_TOLERANCE]=ds.tol_pivot;
        if(ds.strategy>=0)   Control[UMFPACK_STRATEGY]=ds.strategy;

    }
    void dosolver(K *x,K*b,int N,int trans) {
        int ts = UMFPACK_A ;
        if(trans) ts =UMFPACK_At;
         if(verb>2 || verbosity> 9) cout << "dosolver UMFPACK complex/int "<< N << " " << trans << endl;
        for(int k=0,oo=0; k<N;++k, oo+= A->n)
        {
            double * xx = (double *) (void*) x+oo,  *bb = (double *) (void*) b+oo, *zx=0;;
            status= umfpack_zi_solve (ts, Ap, Ai, Ax,Az, xx,zx, bb, zx , Numeric, 0, 0) ;
            if(status) cout << " Error umfpack_di_solve  status  " << status << endl;
        }

    }

    void UpdateState(){
        if( A->GetReDoNumerics() ) cn++;
        if( A->GetReDoSymbolic()) cs++;
        ChangeCodeState(A->n,cs,cn);

    }

    void fac_symbolic(){
        A->CSC(Ap,Ai,Ac);
        Ax= (double *) (void *) Ac;
        Az=0;
         if(verb>2 || verbosity> 9) cout << "fac_symbolic UMFPACK C:  nnz= "  << A->nnz << endl;
        if(Symbolic)  umfpack_zi_free_symbolic (&Symbolic) ;
        status = umfpack_zi_symbolic (A->n, A->m, Ap, Ai, Ax,Az, &Symbolic, 0, 0) ;
        if(status) cout << " Error umpfack umfpack_zi_symbolic  status  " << status << endl;
    }
    void fac_numeric(){
        if(Numeric)   umfpack_zi_free_numeric (&Numeric) ;
         if(verb>2 || verbosity> 9) cout << "fac_numeric UMFPACK C:  nnz= "  << A->nnz << endl;
        status = umfpack_zi_numeric (Ap, Ai, Ax,Az, Symbolic, &Numeric, 0, 0) ;
        if(status) cout << " Error umpfack umfpack_zi_numeric  status  " << status << endl;

    }
    ~VirtualSolverUMFPACK()
    {
        if(Symbolic)  umfpack_zi_free_symbolic (&Symbolic) ;
        if(Numeric)   umfpack_zi_free_numeric (&Numeric) ;
    }
};


template<class Z=int,class K=double>
class VirtualSolverCHOLMOD: public VirtualSolver<Z,K> {
public:
    //  1 unsym , 2 sym, 4 pos , 8 nopos, 16  seq, 32  ompi, 64 mpi ,
    static const int orTypeSol = 2|4|8|16;

    typedef HashMatrix<Z,K>  HMat;
    HMat *HA;
    cholmod_common Common, *cm ;
    cholmod_factor *L ;
    cholmod_sparse *A ;
    cholmod_dense *X = NULL, *B, *W, *R ;

    mutable int status;

    VirtualSolverCHOLMOD(HMat  *AA, const Data_Sparse_Solver & ds,Stack stack ):A(AA) {}
    void dosolver(K *x,K*b,int N,int trans) {assert(0);}
    void fac_symbolic(){assert(0);}
    void fac_numeric(){assert(0);}
    ~VirtualSolverCHOLMOD(){}
    void UpdateState(){}
};

// specilisation
template<>
class  VirtualSolverCHOLMOD<int,double> : public VirtualSolver<int,double> {
public:
    //  1 unsym , 2 sym, 4 pos , 8 nopos, 16  seq, 32  ompi, 64 mpi ,
    static const int orTypeSol = 2|4|8|16;

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
    long verb;

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
    VirtualSolverCHOLMOD(HMat  &HAA,  const Data_Sparse_Solver & ds,Stack stack )
      :HA(&HAA),n(HAA.n),L(0),A(&AA),Ai(0),Ap(0),Ax(0),Ywork(0),Ework(0),cs(0),cn(0),verb(ds.verb)
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
    }

    void dosolver(K *x,K*b,int N,int trans)
    {
        cholmod_dense XX,*X=&XX,B;
        set_cholmod_dense(*X,x,N);
        set_cholmod_dense(B,b,N);

        cout << " dosolver CHOLMOD double "<< N << " " << trans <<  endl;
       cholmod_solve2 (CHOLMOD_A, L, &B, NULL, &X, NULL,
                        &Ywork, &Ework, &c) ;
        if( X !=  &XX) cholmod_free_dense (&X, &c) ;

    }

    void UpdateState(){
        if( HA->GetReDoNumerics() ) cn++;
        if( HA->GetReDoSymbolic()) cs++;
        ChangeCodeState(HA->n,cs,cn);

    }

    void fac_symbolic()
    {
        AA.nzmax=HA->nnz;
        if( HA->half)
         HA->CSR(Ap,Ai,Ax);
        else
         AA.nzmax=HA->CSC_U(Ap,Ai,Ax);

        if(verb>2 || verbosity> 9) cout << "fac_symbolic cholmod R: nnz U=" << AA.nzmax << " nnz= "  << HA->nnz << " " <<  HA->half << endl;

        AA.p=Ap;
        AA.i=Ai;
        AA.x=Ax;
        if(L) cholmod_free_factor (&L, &c) ;		    /* free matrices */
        L = cholmod_analyze (A, &c) ;
     }
    void fac_numeric(){
         if(verb>2 || verbosity> 9) cout << " fac_numeric CHOLMoD double "<< endl;

        cholmod_factorize (A, L, &c) ;

    }
    ~VirtualSolverCHOLMOD()
    {
        if(L) cholmod_free_factor (&L, &c) ;		    /* free matrices */
        cholmod_finish (&c) ;
    }
};

// specialisation
template<>
class  VirtualSolverCHOLMOD<int,std::complex<double> > : public VirtualSolver<int,std::complex<double> >
{
public:
    //  1 unsym , 2 sym, 4 pos , 8 nopos, 16  seq, 32  ompi, 64 mpi ,
    static const int orTypeSol = 2|4|8|16;
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
    long verb;

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

    VirtualSolverCHOLMOD(HMat  &HAA, const Data_Sparse_Solver & ds,Stack stack )
    :HA(&HAA),n(HAA.n),L(0),A(&AA),Ai(0),Ap(0),Ax(0),Ywork(0),Ework(0),cs(0),cn(0),verb(ds.verb)
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

    }

    void dosolver(K *x,K*b,int N,int trans)
    {
         if(verb>2 || verbosity> 9) cout << " dosolver CHOLMoD Complex "<< endl;
        cholmod_dense XX,*X=&XX,B;
        set_cholmod_dense(*X,x,N);
        set_cholmod_dense(B,b,N);
        cholmod_solve2 (CHOLMOD_A, L, &B, NULL, &X, NULL,
                        &Ywork, &Ework, &c) ;
        if( X !=  &XX) cholmod_free_dense (&X, &c) ;

    }

    void UpdateState(){
        if( HA->GetReDoNumerics() ) cn++;
        if( HA->GetReDoSymbolic()) cs++;
        ChangeCodeState(HA->n,cs,cn);

    }

    void fac_symbolic()
    {
        AA.nzmax=HA->CSC_U(Ap,Ai,Ax);
         if(verb>2 || verbosity> 9) cout << "fac_symbolic cholmod C: nnz U=" << AA.nzmax << " nnz= "  << HA->nnz << endl;
        AA.p=Ap;
        AA.i=Ai;
        AA.x=Ax;
        if(L) cholmod_free_factor (&L, &c) ;		    /* free matrices */
        L = cholmod_analyze (A, &c) ;
    }
    void fac_numeric(){
         if(verb>2 || verbosity> 9) cout << " fac_numeric CHOLMoD complex "<< endl;
        cholmod_factorize (A, L, &c) ;

    }
    ~VirtualSolverCHOLMOD()
    {
        if(L) cholmod_free_factor (&L, &c) ;		    /* free matrices */
        cholmod_finish (&c) ;
    }
};

// specilisation
template<>
class  VirtualSolverUMFPACK< SuiteSparse_long,double> : public VirtualSolver< SuiteSparse_long,double> {
public:
    //  1 unsym , 2 sym, 4 pos , 8 nopos, 16  seq, 32  ompi, 64 mpi ,
    static const int orTypeSol = 1|8|16;
    typedef double K;
    typedef  SuiteSparse_long Z;
    typedef HashMatrix<Z,K>  HMat;
    HMat *A;
    void *Symbolic, *Numeric ;
    Z *Ai,*Ap;
    K *Ax;
    int  cs,cn;
    long verb;
    mutable long status;
    double Control[UMFPACK_CONTROL];
    double Info[UMFPACK_INFO];
    VirtualSolverUMFPACK(HMat  &AA, const Data_Sparse_Solver & ds,Stack stack )
    :A(&AA),Symbolic(0),Numeric(0),Ai(0),Ap(0),Ax(0),cs(0),cn(0),verb(ds.verb)
    {
        if(verb>2 || verbosity> 9) cout << " build solver UMFPACK double/long " << endl;
       for(int i=0;i<UMFPACK_CONTROL;i++) Control[i]=0;
        for(int i=0;i<UMFPACK_INFO;i++) Info[i]=0;

        umfpack_di_defaults (Control) ;

        if(ds.verb>4) Control[UMFPACK_PRL]=2;
        if(ds.tol_pivot_sym>0) Control[UMFPACK_SYM_PIVOT_TOLERANCE]=ds.tol_pivot_sym;
        if(ds.tol_pivot>0) Control[UMFPACK_PIVOT_TOLERANCE]=ds.tol_pivot;
        if(ds.strategy>=0)   Control[UMFPACK_STRATEGY]=ds.strategy;

    }
    void dosolver(double *x,double*b,int N,int trans) {
        int ts = UMFPACK_A ;
        if(trans) ts =UMFPACK_At;
        if(verb>2 || verbosity> 9) cout << " dosolver UMFPACK double/long " << N << " " << trans <<endl;

        for(int k=0,oo=0; k<N;++k, oo+= A->n)
        {
            status= umfpack_dl_solve (ts, Ap, Ai, Ax, x+oo, b+oo, Numeric,Control,Info) ;
            if(status) cout << " Error umfpack_di_solve  status  " << status << endl;
            if(verbosity>3)     (void)  umfpack_di_report_info(Control,Info);
        }

    }

    void UpdateState(){
        if( A->GetReDoNumerics() ) cn++;
        if( A->GetReDoSymbolic()) cs++;
        ChangeCodeState(A->n,cs,cn);

    }

    void fac_symbolic(){
        A->CSC(Ap,Ai,Ax);
        if(verb>2 || verbosity> 9) cout << " fac_symbolic UMFPACK double/long " << endl;

        if(Symbolic)  umfpack_di_free_symbolic (&Symbolic) ;
        status = umfpack_dl_symbolic (A->n, A->m, Ap, Ai, Ax, &Symbolic,Control,Info) ;
        if(status) cout << " Error umpfack umfpack_di_symbolic  status  " << status << endl;
    }
    void fac_numeric(){
        if(Numeric)   umfpack_dl_free_numeric (&Numeric) ;
        if(verb>2 || verbosity> 9) cout << " fac_numeric UMFPACK double/long " << endl;

        status = umfpack_dl_numeric (Ap, Ai, Ax, Symbolic, &Numeric, Control,Info) ;
        if(status) cout << " Error umpfack umfpack_di_numeric  status  " << status << endl;
    }
    ~VirtualSolverUMFPACK()
    {
        if(Symbolic)  umfpack_dl_free_symbolic (&Symbolic) ;
        if(Numeric)   umfpack_dl_free_numeric (&Numeric) ;
    }
};

// specilisation
template<>
class  VirtualSolverUMFPACK<SuiteSparse_long,std::complex<double> > : public VirtualSolver< SuiteSparse_long,std::complex<double> > {
public:
    //  1 unsym , 2 sym, 4 pos , 8 nopos, 16  seq, 32  ompi, 64 mpi ,
    static const int orTypeSol = 1|8|16;

    typedef std::complex<double> K;
    typedef  SuiteSparse_long Z;
    typedef HashMatrix<Z,K>  HMat;
    HMat *A;
    void *Symbolic, *Numeric ;
    Z *Ai,*Ap;
    K *Ac;
    double *Ax,*Az;
    int  cs,cn;
    long verb;
    mutable Z status;
    double Control[UMFPACK_CONTROL];
    double Info[UMFPACK_INFO];

    VirtualSolverUMFPACK(HMat  &AA,  const Data_Sparse_Solver & ds,Stack stack )
    :A(&AA),Symbolic(0),Numeric(0),Ai(0),Ap(0),Ax(0),cs(0),cn(0),verb(ds.verb)
    {
        if(verb>2 || verbosity> 9) cout << " build solver UMFPACK complex/long " << endl;
        for(int i=0;i<UMFPACK_CONTROL;i++) Control[i]=0;
        for(int i=0;i<UMFPACK_INFO;i++) Info[i]=0;

        umfpack_zl_defaults (Control) ;

        if(ds.verb>4) Control[UMFPACK_PRL]=2;
        if(ds.tol_pivot_sym>0) Control[UMFPACK_SYM_PIVOT_TOLERANCE]=ds.tol_pivot_sym;
        if(ds.tol_pivot>0) Control[UMFPACK_PIVOT_TOLERANCE]=ds.tol_pivot;
        if(ds.strategy>=0)   Control[UMFPACK_STRATEGY]=ds.strategy;

    }
    void dosolver(K *x,K*b,int N,int trans) {
        int ts = UMFPACK_A ;
        if(trans) ts =UMFPACK_At;
        if(verb>2 || verbosity> 9) cout << " dosolver UMFPACK C/long " << endl;

        for(int k=0,oo=0; k<N;++k, oo+= A->n)
        {
            double * xx = (double *) (void*) x+oo,  *bb = (double *) (void*) b+oo, *zx=0;;
            status= umfpack_zl_solve (UMFPACK_A, Ap, Ai, Ax,Az, xx,zx, bb, zx , Numeric, 0, 0) ;
            if(status) cout << " Error umfpack_di_solve  status  " << status << endl;
        }

    }

    void UpdateState(){
        if( A->GetReDoNumerics() ) cn++;
        if( A->GetReDoSymbolic()) cs++;
        ChangeCodeState(A->n,cs,cn);

    }

    void fac_symbolic(){
        A->CSC(Ap,Ai,Ac);
        Ax= (double *) (void *) Ac;
        Az=0;
        if(verb>2 || verbosity> 9) cout << " fac_symbolic UMFPACK C/long " << endl;

        if(Symbolic)  umfpack_zl_free_symbolic (&Symbolic) ;
        status = umfpack_zl_symbolic (A->n, A->m, Ap, Ai, Ax,Az, &Symbolic, 0, 0) ;
        if(status) cout << " Error umpfack umfpack_zl_symbolic  status  " << status << endl;
    }
    void fac_numeric(){
        if(Numeric)   umfpack_zl_free_numeric (&Numeric) ;
        if(verb>2 || verbosity> 9) cout << " fac_numeric UMFPACK C/long " << endl;

        status = umfpack_zl_numeric (Ap, Ai, Ax,Az, Symbolic, &Numeric, 0, 0) ;
        if(status) cout << " Error umpfack umfpack_zl_numeric  status  " << status << endl;

    }
    ~VirtualSolverUMFPACK()
    {
        if(Symbolic)  umfpack_zl_free_symbolic (&Symbolic) ;
        if(Numeric)   umfpack_zl_free_numeric (&Numeric) ;
    }
};

// specilisation
template<>
class  VirtualSolverCHOLMOD< SuiteSparse_long,double> : public VirtualSolver< SuiteSparse_long,double> {
public:
    //  1 unsym , 2 sym, 4 pos , 8 nopos, 16  seq, 32  ompi, 64 mpi ,
    static const int orTypeSol = 2|4|8|16;
    typedef double K;
    typedef  SuiteSparse_long Z;
    typedef HashMatrix<Z,K>  HMat;
    HMat *HA;
    Z n;
    const int xtype=CHOLMOD_REAL;
    cholmod_common c ;
    cholmod_factor *L ;
    cholmod_sparse AA,*A ;
    Z *Ai,*Ap;
    K *Ax;
    cholmod_dense *Ywork, *Ework  ;
    int  cs,cn;
    long verb;


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
    VirtualSolverCHOLMOD(HMat  &HAA,  const Data_Sparse_Solver & ds,Stack stack )
    :HA(&HAA),n(HAA.n),L(0),A(&AA),Ai(0),Ap(0),Ax(0),Ywork(0),Ework(0),cs(0),cn(0),verb(ds.verb)
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
        AA.itype=CHOLMOD_LONG;
        AA.xtype=CHOLMOD_REAL;
        AA.dtype=CHOLMOD_DOUBLE;
        AA.sorted=1;
        AA.packed=1;

    }

    void dosolver(K *x,K*b,int N,int trans)
    {
        cholmod_dense XX,*X=&XX,B;
        set_cholmod_dense(*X,x,N);
        set_cholmod_dense(B,b,N);

         if(verb>2 || verbosity> 9)  cout << " dosolver CHOLMoD double "<< endl;
        cholmod_l_solve2 (CHOLMOD_A, L, &B, NULL, &X, NULL,
                        &Ywork, &Ework, &c) ;
        if( X !=  &XX) cholmod_free_dense (&X, &c) ;

    }

    void UpdateState(){
        if( HA->GetReDoNumerics() ) cn++;
        if( HA->GetReDoSymbolic()) cs++;
        ChangeCodeState(HA->n,cs,cn);

    }

    void fac_symbolic()
    {
        AA.nzmax=HA->CSC_U(Ap,Ai,Ax);
         if(verb>2 || verbosity> 9)  cout << "fac_symbolic cholmod R: nnz U=" << AA.nzmax << " nnz= "  << HA->nnz << endl;

        AA.p=Ap;
        AA.i=Ai;
        AA.x=Ax;
        if(L) cholmod_l_free_factor (&L, &c) ;            /* free matrices */
        L = cholmod_l_analyze (A, &c) ;
    }
    void fac_numeric(){
         if(verb>2 || verbosity> 9)  cout << " fac_numeric CHOLMoD double "<< endl;

        cholmod_l_factorize (A, L, &c) ;

    }
    ~VirtualSolverCHOLMOD()
    {
        if(L) cholmod_l_free_factor (&L, &c) ;            /* free matrices */
        cholmod_l_finish (&c) ;
    }
};

// specialisation
template<>
class  VirtualSolverCHOLMOD< SuiteSparse_long,std::complex<double> > : public VirtualSolver<int,std::complex<double> >
{
public:
    //  1 unsym , 2 sym, 4 pos , 8 nopos, 16  seq, 32  ompi, 64 mpi ,
    static const int orTypeSol = 2|4|8|16;
    typedef std::complex<double>  K;
    typedef  SuiteSparse_long Z;
    typedef HashMatrix<Z,K>  HMat;
    HMat *HA;
    Z n;
    static const int xtype=CHOLMOD_COMPLEX ;// a complex matrix (ANSI C99 compatible)
    static const int dtype =CHOLMOD_DOUBLE;
    cholmod_common c ;
    cholmod_factor *L ;
    cholmod_sparse AA,*A ;
    Z *Ai,*Ap;
    K *Ax;
    cholmod_dense *Ywork, *Ework  ;
    int  cs,cn;
    long verb;


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

    VirtualSolverCHOLMOD(HMat  &HAA, const Data_Sparse_Solver & ds,Stack stack )
    :HA(&HAA),n(HAA.n),L(0),A(&AA),Ai(0),Ap(0),Ax(0),Ywork(0),Ework(0),cs(0),cn(0),verb(ds.verb)
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
        AA.itype=CHOLMOD_LONG;
        AA.xtype=xtype;
        AA.dtype=dtype;
        AA.sorted=1;
        AA.packed=1;

    }
    void fac_init(){

    }
    void dosolver(K *x,K*b,int N,int trans)
    {
         if(verb>2 || verbosity> 9)  cout << " dosolver CHOLMoD Complex "<< endl;
        cholmod_dense XX,*X=&XX,B;
        set_cholmod_dense(*X,x,N);
        set_cholmod_dense(B,b,N);

        cholmod_l_solve2 (CHOLMOD_A, L, &B, NULL, &X, NULL,
                        &Ywork, &Ework, &c);
        if( X !=  &XX) cholmod_free_dense (&X, &c);
    }

    void UpdateState(){
        if( HA->GetReDoNumerics() ) cn++;
        if( HA->GetReDoSymbolic()) cs++;
        ChangeCodeState(HA->n,cs,cn);

    }

    void fac_symbolic()
    {
        AA.nzmax=HA->CSC_U(Ap,Ai,Ax);
         if(verb>2 || verbosity> 9)  cout << "fac_symbolic cholmod C: nnz U=" << AA.nzmax << " nnz= "  << HA->nnz << endl;
        AA.p=Ap;
        AA.i=Ai;
        AA.x=Ax;
        if(L) cholmod_l_free_factor (&L, &c) ;            /* free matrices */
        L = cholmod_l_analyze (A, &c) ;
    }
    void fac_numeric(){
         if(verb>2 || verbosity> 9)  cout << " fac_numeric CHOLMoD complex "<< endl;
        cholmod_l_factorize (A, L, &c) ;

    }
    ~VirtualSolverCHOLMOD()
    {
        if(L) cholmod_l_free_factor (&L, &c) ;            /* free matrices */
        cholmod_l_finish(&c) ;

    }
};
#endif

void init_UMFPack_solver();
#endif
