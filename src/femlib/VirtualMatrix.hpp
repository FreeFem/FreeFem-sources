
#include "MatriceElementaire.hpp"

inline void MATERROR(int i,const char *cmm)
{
    std::cerr << " MATERROR " << i << " : " << cmm << std::endl;
    std::abort();
}
template<class Z,class R> class HashMatrix ;
template<class TypeIndex,class TypeScalar>
class VirtualMatrix: public RefCounter {
public:
     typedef TypeIndex I;
    typedef TypeScalar R;
    typedef void (*ERRORFunc)(int i,const char *cmm);
    ERRORFunc ERRORHandle;
    class VSolver {public:
        virtual R* solve(R *x,R*b,int N=0,int trans=0) =0;
        virtual ~VSolver(){}
    };

    void ERROR(int i,const char *cmm) const
    {
        if(ERRORHandle) (*ERRORHandle)(i,cmm);
        else MATERROR(i,cmm);
    }
    I n,m; // size of matrix
    VSolver *vsolver;
    bool delvsolver;
   // long state,codeini,codesym,codenum;
    VirtualMatrix(I NN,I MM=-1) : n(NN),m(MM<0 ? n : MM),vsolver(0),delvsolver(false)
    //  state(0),codeini(0),codesym(0),codenum(0)
    {}
   static R *Set2Const(I n,R *x,R c=R()) { std::fill(x,x+n,c); return x;}
    
   R* solve(R *x,R*b,int N=1,int transpo=0) const
    {  if(vsolver) return vsolver->solve(x,b,N,transpo);
        else
        ERROR(1,"VirtualMatrix:: no solver ?????");
        return 0;
    }

  virtual R* addMatMul(R *x,R*Ax) const {  ERROR(1,"VirtualMatrix:: no AddMatMul ?????"); return 0;}
  virtual R* addMatTransMul(R *x,R*Atx) const {  ERROR(1,"VirtualMatrix:: no addMatTransMul ?????");return 0;}
  R* MatMul(R *x,R*Ax) const { return addMatMul(x, Set2Const(n,Ax)); }
  R* MatTransMul(R *x,R*Atx) const { return addMatMul(x, Set2Const(m,Atx));}
   bool WithSolver() const {return vsolver;} // by default no solver

  void SetSolver(VSolver *f=0, bool del = true)  {if(vsolver && delvsolver) delete vsolver; vsolver=f; delvsolver=del; }
/*
  void CheckState(long ci=0,long cs=0, long cn=0)
    {
      if(ci &&  ci != codeini) { codeini=ci; state=0;}// n, nzz fixe
      else if(cs &&  cs != codesym) { codesym=cs; state=1;}//  i,j fixe
      else if(cn &&  cn != codenum) { codenum=cn; state=2;}// a fixe
    };
    
  R* solve(R *x,R *b,int N=1,int trans=0)
  {
    SetState();
    if( state==0) {fac_init(); state=1;}
    if( state==1) {fac_symbolic(); state=2;}
    if( state==2) {fac_numeric();state=3;}
    dosolver(x,b,N,trans);
    return x;
  }
 */
    KN_<R> & MatMul(KN_<R> &ax,const KN_<R> &x) const {
        MatMul((R*)x,(R*)ax);
        return ax;}
     void Solve(KN_<R> & x,const KN_<R> & b) const  {solve(x,b);}
     void SolveT(KN_<R> & x,const KN_<R> & b) const  {solve(x,b,1,1);}
    bool ChecknbLine(int nn) const { return n==nn;}
    bool ChecknbColumn(int mm) const { return m==mm;}

    virtual ~VirtualMatrix(){ if(vsolver && delvsolver) delete vsolver; } // clean solver
    
    virtual size_t size() const =0;
    virtual VirtualMatrix  & operator +=(MatriceElementaire<R> & )=0;
    virtual void operator=(const R & v) =0; // Mise a zero
    virtual ostream& dump (ostream&)  const =0;
    virtual R & diag(I i)=0;
    virtual void SetBC(I i,double tgv)=0;
    virtual R & operator()(I i,I j)=0;
    virtual R * pij(I i,I j) const =0; // Add FH

    virtual HashMatrix<I, R> *toMatriceMorse(bool transpose=false,bool copy=false) const {return 0;} // not
    virtual bool addMatTo(R coef,std::map< pair<int,int>, R> &mij,bool trans=false,int ii00=0,int jj00=0,bool cnj=false,double threshold=0.,const bool keepSym=false)=0;
    virtual R pscal(const KN_<R> & x,const KN_<R> & y) =0 ; // produit scalaire
    virtual double psor(KN_<R> & x,const  KN_<R> & gmin,const  KN_<R> & gmax , double omega) =0;
    virtual void setdiag(const KN_<R> & x)=0 ;
    virtual void getdiag( KN_<R> & x) const =0 ;
    virtual I NbCoef() const {return 0;};
    virtual void setcoef(const KN_<R> & x)=0 ;
    virtual void getcoef( KN_<R> & x) const =0 ;
    virtual bool sym() const {return false;}
 
    virtual  void  resize(I n,I m)  {AFAIRE("~VirtualMatrix::resize");}
   virtual R trace() const {ffassert(n==m);  R t=R(), *p;  for(int i=0; i<n; ++i)  { p=pij(i,i);  if(p) t+= *p;} return t; }
   virtual void SetBC(char *wbc,double tgv) { for (int i=0; i<n; ++i)  if(wbc[i]) SetBC(i,tgv);}

};


template<class TypeIndex=int,class TypeScalar=double>
inline double * ProduitMatVec(const VirtualMatrix<TypeIndex,TypeScalar> *A,TypeScalar *x, TypeScalar *Ax) { return A->MatMul(x,Ax);}
template<class TypeIndex=int,class TypeScalar=double>
inline double * ProduitMatVec(const VirtualMatrix<TypeIndex,TypeScalar> &A,TypeScalar *x, TypeScalar *Ax) { return A.MatMul(x,Ax);}



