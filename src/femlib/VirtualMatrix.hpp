#ifndef VIRTUALMATRIX_HPP__
#define VIRTUALMATRIX_HPP__
#include "MatriceElementaire.hpp"

inline void MATERROR(int i,const char *cmm)
{
    std::cerr << " MATERROR " << i << " : " << cmm << std::endl;
    throw(ErrorExec("MATERROR",1));
    std::abort();
}
template<class Z,class R> class HashMatrix ;
template<class TypeIndex,class TypeScalar>
class VirtualMatrix: public RefCounter, public RNM_VirtualMatrix<TypeScalar,TypeIndex> {
public:
    typedef TypeIndex I;
    typedef TypeScalar R;
    // 1 unsym , 2 herm, 4 sym, 8 pos , 16 nopos, 32  seq, 64  ompi, 128 mpi
    static const int  TS_unsym=1, TS_herm=2, TS_sym=4, TS_def_positif=8,  TS_not_def_positif=16, TS_sequental = 32, TS_mpi = 64;// for verification
    static const int TS_SYM=1,TS_DEF_POS=2,TS_PARA=4;
    typedef VirtualMatrix<I,R> VMat ;
    typedef void (*ERRORFunc)(int i,const char *cmm);
    ERRORFunc ERRORHandle;
    class VSolver {public:
        virtual R* solve(R *x,R*b,int N=0,int trans=0) =0;
        virtual void factorize(int st){ if(verbosity>4)  cout << " ** warning no factorisation" << st << endl;};
        virtual ~VSolver(){}
        int count;
        VSolver(): count(0) {}
        VSolver * copyptr(){ count++; return this;}
        void destroy() { if(count==0) delete this; else count --;}
       // virtual VSolver * clone() =0;  to hard
    };

    void ERROR(int i,const char *cmm) const
    {
        if(ERRORHandle) (*ERRORHandle)(i,cmm);
        else MATERROR(i,cmm);
    }
    I n,m; // size of matrix
    bool symetric,positive_definite;   // for cholesky or CG
    
    VSolver *vsolver;
    bool delvsolver;
    
   // long state,codeini,codesym,codenum;
    VirtualMatrix(I NN,I MM=-1,bool sym=false,bool dp=false) : RNM_VirtualMatrix<R,I> (NN,MM),
    n(NN),m(this->M),symetric(sym),positive_definite(dp),vsolver(0),delvsolver(false)
    {}
   virtual void setsdp(int sym,bool dp) { symetric=(sym>0); positive_definite=dp;}
   static R *Set2Const(I n,R *x,R c=R()) { std::fill(x,x+n,c); return x;}
    
   R* solve(R *x,R*b,int N=1,int transpo=0) const
    {  if(vsolver) return vsolver->solve(x,b,N,transpo);
        else
        ERROR(1,"VirtualMatrix:: no solver ?????");
        return 0;
    }
  virtual R* addMatMul(R *x,R*Ax,bool Transpose,I sx=1,I sAx=1) const  {  ERROR(1,"VirtualMatrix:: no AddMatMul ?????"); return 0;}
  R* addMatMul(R *x,R*Ax) const {  return addMatMul(x,Ax,0); }
  R* addMatTransMul(R *x,R*Atx) const {  return addMatMul(x,Atx,1); }
  R* MatMul(R *x,R*Ax) const { return addMatMul(x, Set2Const(n,Ax)); }
  R* MatTransMul(R *x,R*Atx) const { return addMatMul(x, Set2Const(m,Atx));}
   bool WithSolver() const {return vsolver;} // by default no solver
  //  void CloneSolver() { return vsolver ? vsolver->clone(): 0;}// for copy matrix  to hard ...
    void SetSolver(const VirtualMatrix<I,R> & A)
    {
        
        if(A.vsolver) SetSolver(vsolver->copyptr(),A.delvsolver);
        else SetSolver(0);
    }
  void SetSolver(VSolver *f=0, bool del = true)  {
      if(verbosity>4)  cout<< " ## SetSolver "<< this << " " << vsolver <<" "<<  f << endl;
      if(vsolver && delvsolver)
           vsolver->destroy();
      vsolver=f;
      delvsolver=del;
      // cout << "\n *** type SetSolver = " << typeid(f).name() << endl; 
  }

    KN_<R> & MatMul(KN_<R> &ax,const KN_<R> &x) const {
        MatMul((R*)x,(R*)ax);
        return ax;}
    void addMatMul(const KN_<R> &  x, KN_<R> & y) const { ffassert(Checknm(y.N(),x.N())); addMatMul(x,y,0,x.step,y.step);}
    void addMatTransMul(const KN_<R> & x , KN_<R> & y ) const {ffassert(Checknm(x.N(),y.N())); addMatMul(x,y,1,x.step,y.step);}
    void Solve(KN_<R> & x,const KN_<R> & b) const  {solve(x,b);}
    void SolveT(KN_<R> & x,const KN_<R> & b) const  {solve(x,b,1,1);}
                                                    
    bool ChecknbLine(int nn) const { return this->n==nn;}
    bool ChecknbColumn(int mm) const { return this->m==mm;}
   bool Checknm(int nn,int mm) const { return nn==this->n && mm==this->m;}
    virtual ~VirtualMatrix(){
        if(verbosity>99999) cout << " **  ~VirtualMatrix " << this << endl;
        if(vsolver && delvsolver)  vsolver->destroy();//   Warning the solver is del afer the Matrix
    } // clean solver
    
    virtual size_t size() const {return 0; };
    virtual VirtualMatrix  & operator +=(MatriceElementaire<R> & ){AFAIRE("VirtualMatrix::+=");}
    virtual void operator=(const R & v){AFAIRE("VirtualMatrix::=v");} // Mise a zero
    virtual ostream& dump (ostream&)  const {cout << mpirank << " BUG virtualmatrix " << this << endl; AFAIRE("VirtualMatrix::dump");}
    virtual R & diag(I i){AFAIRE("VirtualMatrix::diab");}
    virtual void SetBC(I i,double tgv){AFAIRE("VirtualMatrix::setbc");}
    virtual R & operator()(I i,I j){AFAIRE("VirtualMatrix::()(i,j)");}
    virtual R * pij(I i,I j) const {AFAIRE("VirtualMatrix::pij");} // Add FH

    virtual HashMatrix<I, R> *toMatriceMorse(bool transpose=false,bool copy=false) const {return 0;} // not
    virtual bool addMatTo(R coef,HashMatrix<I,R> &mij,bool trans=false,int ii00=0,int jj00=0,bool cnj=false,double threshold=0.,const bool keepSym=false)
    {AFAIRE("VirtualMatrix::addMatTo");};
    virtual R pscal(const KN_<R> & x,const KN_<R> & y) {AFAIRE("VirtualMatrix::pscal");} ; // produit scalaire
    virtual double psor(KN_<R> & x,const  KN_<R> & gmin,const  KN_<R> & gmax , double omega) {AFAIRE("VirtualMatrix::psor");}
    virtual void setdiag(const KN_<R> & x){AFAIRE("VirtualMatrix::setdiag");} ;
    virtual void getdiag( KN_<R> & x) const {AFAIRE("VirtualMatrix::getdiag");}
    virtual I NbCoef() const {return 0;};
    virtual void setcoef(const KN_<R> & x){AFAIRE("VirtualMatrix::setcoef");}
    virtual void getcoef( KN_<R> & x) const {AFAIRE("VirtualMatrix::getcoef");}
    virtual bool sym() const {return false;}
 
    virtual  void  resize(I n,I m)  {AFAIRE("VirtualMatrix::resize");}
    virtual void clear()
    { I NN=this->N, MM=this->M;
        resize(0,0);
        resize(NN,MM);
    }
   virtual R trace() const {ffassert(this->n==this->m);  R t=R(), *p;  for(int i=0; i<this->n; ++i)  { p=pij(i,i);  if(p) t+= *p;} return t; }
   virtual void SetBC(char *wbc,double tgv) { for (int i=0; i<this->n; ++i)  if(wbc[i]) SetBC(i,tgv);}
  //  void init(int nn=0,int mm=0) { VMat *p=new VMat(nn,mm);  }
  
    struct  plusAx { const VirtualMatrix * A; const KN_<R>   x;
        plusAx( const VirtualMatrix * B,const KN_<R> &  y) :A(B),x(y)
        { if(B) { ffassert(B->ChecknbColumn(y.N())); } }
        void call(KN_<R> &ax,int init=0) const { if(init) ax=R();  A->addMatMul(x,ax,0,x.step,ax.step); }
    };
    
    
    struct  plusAtx { const VirtualMatrix * A; const KN_<R>   x;
        plusAtx( const VirtualMatrix * B,const KN_<R> &  y) :A(B),x(y)
        { if(B) { ffassert(B->ChecknbLine(y.N())); } }
        void call(KN_<R> &ax,int init=0) const  {if(init) ax=R(); A->addMatMul(x,ax,1,x.step,ax.step); }
    };
    
    struct  solveAxeqb { const VirtualMatrix * A; const KN_<R>   b;
        solveAxeqb( const VirtualMatrix * B,const KN_<R> &  y) :A(B),b(y)
        { if(B) { ffassert(B->ChecknbColumn(y.N())); } }
            void call(KN_<R> &ax) const {ffassert(ax.contiguous() &&b.contiguous());   A->Solve(ax,b); }
    };
    
    struct  solveAtxeqb { const VirtualMatrix * A; const KN_<R>   b;
        solveAtxeqb( const VirtualMatrix * B,const KN_<R> &  y) :A(B),b(y)
        { if(B) { ffassert(B->ChecknbColumn(y.N())); } }
            void call(KN_<R> &ax) const {ffassert(ax.contiguous() &&b.contiguous());   A->SolveT(ax,b); }
    };

    VirtualMatrix(const VirtualMatrix<I,R> &A)
     : RNM_VirtualMatrix<TypeScalar,TypeIndex>(A)  { operator=(A); }
    void operator=(const VirtualMatrix<I,R> &A)
    {
        n=A.n;
        m=A.m; // size of matrix
        symetric=A.symetric;
        positive_definite=A.positive_definite;   // for cholesky or CG
        delvsolver = A.delvsolver;
        vsolver = A.vsolver && delvsolver ? A.vsolver->copyptr() : 0;
    }

};


template<class TypeIndex=int,class TypeScalar=double>
inline double * ProduitMatVec(const VirtualMatrix<TypeIndex,TypeScalar> *A,TypeScalar *x, TypeScalar *Ax) { return A->MatMul(x,Ax);}
template<class TypeIndex=int,class TypeScalar=double>
inline double * ProduitMatVec(const VirtualMatrix<TypeIndex,TypeScalar> &A,TypeScalar *x, TypeScalar *Ax) { return A.MatMul(x,Ax);}


#endif
