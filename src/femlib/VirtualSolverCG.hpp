#ifndef __VirtualSolverCG_HPP__
#define __VirtualSolverCG_HPP__

#include <fstream>
#include "CG.hpp"
#include <cmath>
#include "HashMatrix.hpp"
#include <vector>
#include <complex>
#include "VirtualSolver.hpp"

template<class I=int,class K=double>
class SolverCG: public VirtualSolver<I,K> {
public:
    
    typedef HashMatrix<I,K>  HMat;
    HMat *A;
    int verb,itermax,erronerr;
    double eps;
    SolverCG(HMat  &AA,double eeps=1e-6,int eoe=1,int v=1,int itx=0)
    :A(&AA),verb(v),itermax(itx>0?itx:A->n/2),erronerr(eoe),eps(eeps)
    {
        cout << " SolverCG  " << A->n << " "<<  A->m <<" eps " << eps << " eoe " << eoe << " v " << verb << " " << itermax <<endl;
        assert(A->n == A->m);
    }
 
    SolverCG(HMat  &AA,const Data_Sparse_Solver & ds)
    :A(&AA),verb(ds.verb),itermax(ds.itmax>0 ?ds.itmax:A->n),erronerr(1),eps(ds.epsilon)
    {
        std::cout << " SolverCG  " << A->n << "x"<<  A->m <<" eps " << eps << " eoe " << erronerr
                  << " v " << verb << "itsmx " << itermax <<endl;
        assert(A->n == A->m);
    }
    
    void SetState(){}
    
    struct HMatVirt: CGMatVirt<I,K> {
        HMat *A;
        int t;
        HMatVirt(HMat *AA,int trans) :CGMatVirt<I,K>(AA->n),A(AA),t(trans) {}
        K * addmatmul(K *x,K *Ax) const { return A->addMatMul(x,Ax,t);}
    };
    struct HMatVirtPreconDiag: CGMatVirt<I,K> {
        HMat *A;
        HMatVirtPreconDiag(HMat *AA) :CGMatVirt<I,K>(AA->n),A(AA){}
        K * addmatmul(K *x,K *Ax) const { for(int i=0; i<A->n; ++i) Ax[i] += std::abs((*A)(i,i)) ? x[i]/(*A)(i,i): x[i]; return Ax;}
    };
    void dosolver(K *x,K*b,int N,int trans)
    {
        std::cout <<" SolverCG::dosolver" << N<< " "<< eps << " "<< itermax << " "<< verb << std::endl;
        HMatVirt AA(A,trans);
        HMatVirtPreconDiag CC(A);
        int err=0;
        for(int k=0,oo=0; k<N; ++k, oo+= A->n )
        {
            int res=ConjugueGradient(AA,CC,b+oo,x+oo,itermax,eps,verb);
            if ( res==0 ) err++;
        }
        if(err && erronerr) {  std::cerr << "Error: ConjugueGradient do not converge nb end ="<< err << std::endl; assert(0); }
    }

};


template<class I=int,class K=double>
class SolverGMRES: public VirtualSolver<I,K> {
public:
    
    typedef HashMatrix<I,K>  HMat;
    HMat *A;
    long verb,itermax,restart,erronerr;
    double eps;
    SolverGMRES(HMat  &AA,double eeps=1e-6,int eoe=1,int v=1,int rrestart=50,int itx=0)
    :A(&AA),verb(v),itermax(itx>0?itx:A->n/2),restart(rrestart),
     erronerr(eoe),eps(eeps)
    {assert(A->n == A->m);}

    SolverGMRES(HMat  &AA,const Data_Sparse_Solver & ds)
    :A(&AA),verb(ds.verb),itermax(ds.itmax>0 ?ds.itmax:A->n),restart(ds.NbSpace),erronerr(1),eps(ds.epsilon)
    {
        std::cout << " SolverGMRES  " << A->n << "x"<<  A->m <<" eps " << eps << " eoe " << erronerr
        << " v " << verb << "itsmx " << itermax <<endl;
        assert(A->n == A->m);
    }

    
    void SetState(){}
    
    struct HMatVirt: CGMatVirt<I,K> {
        HMat *A;
        int t;
        HMatVirt(HMat *AA,int trans=0) :CGMatVirt<I,K>(AA->n),A(AA),t(trans){}
        K * addmatmul(K *x,K *Ax) const { return  A->addMatMul(x,Ax,t);}
    };
    struct HMatVirtPreconDiag: CGMatVirt<I,K> {
        HMat *A;
        HMatVirtPreconDiag(HMat *AA) :CGMatVirt<I,K>(AA->n),A(AA){}
        K * addmatmul(K *x,K *Ax) const { for(int i=0; i<A->n; ++i) Ax[i] += std::abs((*A)(i,i)) ? x[i]/(*A)(i,i): x[i]; return Ax;}
    };
    void dosolver(K *x,K*b,int N=0,int trans=1)
    {
        std::cout <<" SolverCG::dosolver" << N<< " "<< eps << " "<< itermax << " "<< verb << std::endl;
        HMatVirt AA(A,trans);
        HMatVirtPreconDiag CC(A);
        int err=0;
        for(int k=0,oo=0; k<N; ++k, oo+= A->n )
        {
            bool res=fgmres(AA,CC,b+oo,x+oo,eps,itermax,restart);
        
            if ( ! res ) err++;
        }
        if(err && erronerr) {  std::cerr << "Error: fgmres do not converge nb end ="<< err << std::endl; assert(0); }
    }
    
};

#endif
