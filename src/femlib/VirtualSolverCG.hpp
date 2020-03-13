#ifndef __VirtualSolverCG_HPP__
#define __VirtualSolverCG_HPP__

#include <fstream>
#include "CG.hpp"
#include <cmath>
#include "HashMatrix.hpp"
#include <vector>
#include <complex>
#include "VirtualSolver.hpp"
#include "AFunction.hpp"

template<class I=int,class K=double>
struct HMatVirtPrecon: CGMatVirt<I,K> {
    typedef HashMatrix<I,K>  HMat;
    HMat *A;
    //bool diag;
    //  Preco FF++
    Expression xx_del, code_del;
    const E_F0 * precon;
    Stack stack;
    KN<K> *xx;
    K *diag1;
    KN<int> *wcl;
    double tgv;
    int ntgv;
    HMatVirtPrecon(HMat *AA,const Data_Sparse_Solver * ds,Stack stk=0) :CGMatVirt<I,K>(AA->n),A(AA),//diag(!ds || !ds->precon|| !stk),
    xx_del(0),code_del(0),precon(0),stack(stk),wcl(0),xx(0),diag1(0),tgv(1e30),ntgv(0)
    {
        I n = A->n;
        if(ds) {
            tgv = ds->tgv;
            int ntgv1;
            double tgvm = A->gettgv(&ntgv1);
            if(ntgv1)
            {
                wcl = new KN<int>(n);
                double tgve =  tgvm;// ds->tgv;
                if( tgve <=0) tgve = 1e200;// no tgv
                ntgv =0;;
                for (int i=0;i<n;i++)
                ntgv += (*wcl)[i] = real((*A)(i,i))==tgve;
            }
            //cout << " HMatVirtPrecon: ntgv = " << ntgv << endl;
        }
        if(stack &&  ds->precon)
        {  cout << " with Preco " << endl;
            const OneOperator * C = static_cast<const OneOperator *>(ds->precon);
           
            xx = new KN<K>(n);
            
            WhereStackOfPtr2Free(stack)=new StackOfPtr2Free(stack);// FH mars 2005
            Type_Expr te_xx(CPValue(*xx));
            xx_del=te_xx.second;
            C_F0 e_xx(te_xx); // 1 undelete pointer
            code_del= C->code(basicAC_F0_wa(e_xx));
            precon =  to<KN_<K> >(C_F0(code_del,*C));// 2 undelete pointer
            throwassert(precon);
            if (verbosity>4 ) cout << " ## Precon  GC/GMRES : nb tgv in mat = "<< ntgv << " " << tgv << " " << this << endl;
        }
        else {stack=0;
            diag1 = new K[n];
            for(int i=0; i<n;++i)
                diag1[i]=(*A)(i,i);
              for(int i=0; i<n;++i)
                  if( std::norm(diag1[i]) < 1e-60) diag1[i]=1;
                  else diag1[i]=1./diag1[i];
            if (verbosity>4) cout << " ## Precon Diag GC/GMRES : nb tgv in mat = "<< ntgv << " " << tgv << this << endl;
        } // no freefem++ precon
    }
    K * addmatmul(K *x,K *Ax) const
    {
        int n = A->n;
        if(diag1)
         for(int i=0; i<n; ++i)
            Ax[i] += diag1[i]*x[i]; //std::norm((*A)(i,i))>1e-60 ? x[i]/(*A)(i,i): x[i];
        else {// Call Precon ff++
            KN<K> &ffx=*xx;
            KN_<K> ax(Ax,n);
            ffx=x;
            // cout << x[0] << "  ";
            ffx=GetAny<KN_<K> >((*precon)(stack));
            WhereStackOfPtr2Free(stack)->clean();
            //    cout << (xx)[0] << "  " << endl;
            K dii;
            if(wcl)
            for (int i=0;i<A->n;i++)
                if((*wcl)[i]) ffx[i] = x[i]/tgv ;
            ax += ffx ;
            
        }
        return Ax;}
    void  SetInitWithBC(K*rhs,K *x) const
    {
        if(wcl)
        for (int i=0;i<A->n;i++)
        if( (*wcl)[i])
        x[i] = rhs[i]/tgv;
    }
    int * pwcl() const {return wcl ?  (int*) *wcl  : 0; ;}
    ~HMatVirtPrecon()
    {
        if(verbosity>99) cout << " ## ~HMatVirtPrecon "<< this << endl;
        if(xx) delete xx;
        if( diag1) delete [] diag1;
        if(wcl) delete wcl;
        if(stack) WhereStackOfPtr2Free(stack)->clean(); // FH mars 2005
        if(xx_del) delete  xx_del;
        if(code_del) delete  code_del;
    }
};
template<class I=int,class K=double>
class SolverCG: public VirtualSolver<I,K> {
public:
    // 1 unsym , 2 herm, 4 sym, 8 pos , 16 nopos, 32  seq, 64  ompi, 128 mpi
    static const int orTypeSol = 1|2|4|8|32;

    typedef HashMatrix<I,K>  HMat;
    HMat *A;
    CGMatVirt<I,K> *pC;
    int verb,itermax,erronerr;
    double eps;
    SolverCG(HMat  &AA,double eeps=1e-6,int eoe=1,int v=1,int itx=0)
    :A(&AA),pC(0),verb(v),itermax(itx>0?itx:A->n/2),erronerr(eoe),eps(eeps)
    {
        if(verb>4)
            cout << " ## SolverCG  " << A->n << " "<<  A->m <<" eps " << eps << " eoe " << eoe << " v " << verb << " " << itermax <<endl;
        pC = new HMatVirtPreconDiag(A);
        assert(A->n == A->m);
    }
    
    SolverCG(HMat  &AA,const Data_Sparse_Solver & ds,Stack stack)
    :A(&AA),pC(0),verb(ds.verb),itermax(ds.itmax>0 ?ds.itmax:A->n),erronerr(1),eps(ds.epsilon)
    {
        if(verb>4)
            std::cout << " ## SolverCG  " << A->n << "x"<<  A->m <<" eps " << eps << " eoe " << erronerr
        << " v " << verb << " itmax " << itermax <<endl;
        assert(A->n == A->m);
        
        pC = new HMatVirtPrecon<I,K>(A,&ds,stack);
    }
    SolverCG() {delete pC;}
    void UpdateState(){}
    
    struct HMatVirt: CGMatVirt<I,K> {
        HMat *A;
        int t;
        HMatVirt(HMat *AA,int trans) :CGMatVirt<I,K>(AA->n),A(AA),t(trans) {}
        K * addmatmul(K *x,K *Ax) const { return A->addMatMul(x,Ax,t);}
    };
    struct HMatVirtPreconDiag: CGMatVirt<I,K> {
        HMat *A;
        K *d;
        HMatVirtPreconDiag(HMat *AA) :CGMatVirt<I,K>(AA->n),A(AA),d(A->n){
            for(int i=0; i<A->n; ++i) d[i] += std::abs((*A)(i,i)) ? 1./(*A)(i,i): 1.;
        }
        K * addmatmul(K *x,K *Ax) const {
            for(int i=0; i<A->n; ++i)
                Ax[i] += d[i]* x[i];
            return Ax;}
        ~HMatVirtPreconDiag() { delete d;}
    private:// not copy ..
        HMatVirtPreconDiag(const HMatVirtPreconDiag &AA);
        void operator=(const HMatVirtPreconDiag &AA);
    };
    void dosolver(K *x,K*b,int N,int trans)
    {
        if(verb>2|| verbosity>9)
        std::cout <<"   SolverCG::dosolver" << N<< " "<< eps << " "<< itermax << " "<< verb << std::endl;
        HMatVirt AA(A,trans);
        //HMatVirtPreconDiag CC(A);
        int err=0;
        for(int k=0,oo=0; k<N; ++k, oo+= A->n )
        {
            pC->SetInitWithBC(b+oo,x+oo);
            int res=ConjugueGradient(AA,*pC,b+oo,x+oo,itermax,eps,verb);
            if ( res==0 ) err++;
        }
        if(err && erronerr) {  std::cerr << "Error: ConjugueGradient do not converge nb end ="<< err << std::endl; assert(0); }
    }
    ~SolverCG() {delete pC;}
};


template<class I=int,class K=double>
class SolverGMRES: public VirtualSolver<I,K> {
public:
    // 1 unsym , 2 herm, 4 sym, 8 pos , 16 nopos, 32  seq, 64  ompi, 128 mpi
    static const int orTypeSol = 1|2|4|8|16|32;

    typedef HashMatrix<I,K>  HMat;
    HMat *A;
    CGMatVirt<I,K> *pC;
    long verb,itermax,restart,erronerr;
    double eps;
    SolverGMRES(HMat  &AA,double eeps=1e-6,int eoe=1,int v=1,int rrestart=50,int itx=0)
    :A(&AA),pC(0), verb(v),itermax(itx>0?itx:A->n/2),restart(rrestart),
    erronerr(eoe),eps(eeps)
    {assert(A->n == A->m);
        pC = new HMatVirtPrecon<I,K>(A);
    }
    
    SolverGMRES(HMat  &AA,const Data_Sparse_Solver & ds,Stack stack)
    :A(&AA),pC(0),verb(ds.verb),itermax(ds.itmax>0 ?ds.itmax:A->n),restart(ds.NbSpace),erronerr(1),eps(ds.epsilon)
    {
        if(verb>4)
            std::cout << " ## SolverGMRES  " << A->n << "x"<<  A->m <<" eps " << eps << " eoe " << erronerr
        << " v " << verb << "  itsmx " << itermax <<endl;
        assert(A->n == A->m);
        pC = new HMatVirtPrecon<I,K>(A,&ds,stack);
    }
    
    ~SolverGMRES() {delete pC;pC=0;}
    void UpdateState(){}
    
    struct HMatVirt: CGMatVirt<I,K> {
        HMat *A;
        int t;
        HMatVirt(HMat *AA,int trans=0) :CGMatVirt<I,K>(AA->n),A(AA),t(trans){}
        K * addmatmul(K *x,K *Ax) const { return  A->addMatMul(x,Ax,t);}
    };
    
    void dosolver(K *x,K*b,int N=0,int trans=1)
    {
        if(verbosity>9 || verb> 2)
            std::cout <<" ##  SolverGMRES::dosolver" << N<< " "<< eps << " "<< itermax << " "<< verb << std::endl;
        HMatVirt AA(A,trans);
        //HMatVirtPreconDiag CC(A);
        int err=0;
        for(int k=0,oo=0; k<N; ++k, oo+= A->n )
        {
            pC->SetInitWithBC(b+oo,x+oo);
            
            bool res=fgmres(AA,*pC,1,b+oo,x+oo,eps,itermax,restart,verb,pC->pwcl());
            
            if ( ! res ) err++;
        }
        if(err && erronerr) {  std::cerr << "Error: fgmres do not converge nb end ="<< err << std::endl; assert(0); }
    }
    
};

#endif
