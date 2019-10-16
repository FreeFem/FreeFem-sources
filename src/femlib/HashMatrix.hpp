#ifndef __HashMatrix_HPP__
#define __HashMatrix_HPP__

#include <cstddef>
#include <algorithm>
#include <map>
#include <list>

#include <tuple>
using std::tuple;
#include <iostream>
#include <cassert>
#include <complex>
#include <cstdint>
#include "RNM.hpp"
#include "RefCounter.hpp"
#include "VirtualMatrix.hpp"

using std::max;
using std::min;
using std::cout;
using std::endl;
using std::pair;
using std::swap;
using std::complex;
using std::list;
extern  long  verbosity  ;

void init_HashMatrix ();
template<class Z,class ZZ>
inline uint64_t  roll64(Z y,ZZ r){uint64_t x=y; r %= 64; return (x<<r) | (x << (64-r)) ;}



int WhichMatrix(istream & f);
template<class TypeIndex,class TypeScalaire> class HashMatrix;

template<class I,class R,class K>
void Addto(HashMatrix<I,R> *P0,const HashMatrix<I,K> *PA,R (*f)(K) ,bool trans=false, I ii00=0,I jj00=0);



template<class TypeIndex,class TypeScalaire>
class HashMatrix : public VirtualMatrix<TypeIndex,TypeScalaire>
{
    static void  HeapSort(TypeIndex *ii,TypeIndex *jj, TypeScalaire *aij,long n);


public:

    template<class R> static void conj(R *x,TypeIndex nnz){}
    static void conj(complex<double> *x,TypeIndex nnz){ for(int k=0; k<nnz; ++k) x[k]=std::conj(x[k]) ; }
    static void conj(complex<float> *x,TypeIndex nnz){for(int k=0; k<nnz; ++k) x[k]=std::conj(x[k]) ; }
    static double conj(double x){return x;}
    static float conj(float  x){return x;}
    static complex<double> conj(const complex<double> &x){return std::conj(x);}
    static complex<float> conj(const complex<float> &x){return std::conj(x);}

    typedef TypeIndex I;
    typedef TypeScalaire R;
    typedef uint64_t uniquecodeInt;
    static const int type_isdeleted=-1,type_HM=0,type_COO=1, type_CSR=2,type_CSC=3;
    static const int  unsorted=0, sorted_ij=type_CSR,sorted_ji=type_CSC;
    typedef size_t iterator;
    typedef size_t Hash;
    typedef pair<size_t,size_t> Pair;
    size_t nnz,nnzmax,nhash;
    mutable size_t nbcollision,nbfind;
    mutable double matmulcpu;
    I * i,*j;
    I *p;
    R * aij;
    size_t * head;
    size_t * next;
    bool half;
    int state,type_state;
    size_t nbsort;
    I sizep;
    int lock;
    int fortran; // index start a one ..
    mutable int re_do_numerics,re_do_symbolic;
    static const  size_t empty= (I) -1;
     // for  Dirichlet BC
    double tgv;
    I ntgv;

    I NbCoef() const {return  (I) nnz;}
    void setcoef(const KN_<R> & x){KN_<R>c(this->aij,nnz); ffassert(x.SameShape(c));
        if( x.constant())  c=x[0];   else c = x;}
    void getcoef( KN_<R> & x) const {ffassert(x.N()==(I) nnz);x =KN_<R>(this->aij,nnz);}

    void setdiag(const KN_<R> & d);
    void getdiag( KN_<R> & d) const;
    R pscal(R *x,R *y,I sx=1,I sy=1);
    R pscal(const KN_<R> & x,const KN_<R> & y) { return pscal(x,y,(I) x.step,(I) y.step);}
    void SetMorse();
    void UnSetMorse();
    uniquecodeInt CodeIJ() const ;
    void init(I nn,I mm=-1,size_t nnnz=0,bool halff=false);
    HashMatrix(I nn,I mm=-1,I nnnz=0,bool halff=false);
    HashMatrix(istream & f,int cas=-1);
    HashMatrix(KNM_<R> F,double  threshold=1e-30);
    HashMatrix(bool Half,I nn);
    HashMatrix(const HashMatrix& A);
    HashMatrix(I nn,const R *diag);
    void RenumberingInv(KN_<I> II,KN_<I> JJ);
    void Renumbering(I nn,I mm,KN_<I> II,KN_<I> JJ);
    void RemoveDoubleij(int kk);// remove
    template<class R,class K> static R cast_funct(K x) { return (R) x;}

    template<class J,class K> HashMatrix(const HashMatrix<J,K> &A , R (*ff)(K) );
    template<class J> HashMatrix(const HashMatrix<J,R> &A );

    template<class II,class RR> HashMatrix & operator=(const HashMatrix<II,RR>& A )
    {
        if( (const void*)  this == (const void*) & A) return *this;

        set(A.n,A.m,A.half,A.nnz,A.i,A.j,A.aij,A.fortran,cast_funct<R,RR>);
        return *this;
    }

    template<class II,class RR> HashMatrix & operator+=(const HashMatrix<II,RR>& A );


    int IsTrianglulare() const ;

    void CheckUnLock(const char * cmm)
    {
        if( lock)
        {
            std::cerr << " Sorry Forbidder operation ( "
            << cmm<< " ) matix is lock " << endl;
            assert(0);
        }
    }
    void setp(I sp);
    void resize(I nn, I mm=0)  {resize(nn,mm,nnz); }

    void resize(I nn, I mm,size_t nnnz, double tol = -1., bool sym=false );
    void SymmetrizePattern(); // To do for Suzuki , Paradiso
    void clear();
    Hash hash(size_t ii,size_t jj) const{ return ( (ii-fortran)+ (jj-fortran)*this->n )%nhash; }

    void setfortran(int yes);

    size_t find(size_t ii,size_t jj) const { return find(ii,jj, hash(ii,jj)); }

    size_t  find(I ii,I jj,Hash h) const
    {
        nbfind++;
        for (size_t k=head[h];k!=empty;k=next[k])
        {
            ++nbcollision;
            if( ii== i[k] && jj==j[k] ) return k;
        }
        return empty;
    }

    size_t insert(I ii, I jj,const R & aa)
    {
        state=unsorted;
        Hash h= hash(ii,jj);
        size_t k=find(ii,jj,h);
        if(k==empty)
            k =simpleinsert(ii,jj,h);
        aij[k] += aa;
        return k;
    }

    template<typename T>                 static void HMresize(T *&t,size_t no,size_t nn);
    template<typename T,typename TT>     static void HMcopy( T *dst,const TT *from, size_t nn);
    template<typename T,typename TT>     static void HMcopy( T *dst,const TT *from, size_t nn, T (*ff)(TT) );

    bool do2Triangular(bool lower) ; //  put half tp lower or upper
    void dotranspose();
    void Increaze(size_t nnznew=0,size_t newnnz=0);// newnnz<0 => newnnz is set to nnz (change value of nnz)
    void ReHash();
    size_t size() const { return nnz;}
    size_t simpleinsert(I ii, I jj,Hash &h);
    R   *pij(I ii,I jj)  const;
    R   *npij(I ii,I jj);   // with add if no term ii,jj

    R & diag(I ii)  { return operator()(ii,ii);}
    R   diag(I ii) const  { return operator()(ii,ii);}
    R   operator()(I ii,I jj)  const
    {
        Hash h = hash(ii,jj);
        size_t k = find(ii,jj,h);
        if(k==empty) return R();
        else return aij[k];
    }

    R  & operator()(I ii,I jj)
    {
        return *npij(ii,jj);
    }

    R    operator()(pair<I,I> ij)  const {return operator()(ij.first,ij.second);}
    R &  operator()(pair<I,I> ij)  {return operator()(ij.first,ij.second);}
    R    operator[](pair<I,I> ij)  const {return operator()(ij.first,ij.second);}
    R &  operator[](pair<I,I> ij)  {return operator()(ij.first,ij.second);}
    ~HashMatrix();


    void Sortij();
    void Sortji();

    void set(I  nn,I  mm,bool hhalf,size_t nnnz, I  *ii, I *jj, R  *aa,int f77=0,int tcsr=0);
    template<class II>          void set(II nn,II mm,bool hhalf,size_t nnnz, II *ii, II*jj, R  *aa,int f77);
    template<class II,class RR> void set(II nn,II mm,bool hhalf,size_t nnnz, II *ii, II*jj, RR *aa,int f77,R (*ff)(RR));

    void Add(const HashMatrix<I,R> *PA,R coef=R(1),bool trans=false, I ii00=0,I jj00=0);

    HashMatrix &operator=(const HashMatrix &A) ;
    HashMatrix &operator+=(const HashMatrix &A) ;// {Add(&A); return *this;};
    HashMatrix &operator-=(const HashMatrix &A) ; //{Add(&A,R(-1.)); return *this;}



    void operator*=(R v);
    void operator=(const R & v);
    void HM();// un sorted ... Default  type ..
    void COO();
    void COO(I *& IA, I *& IJ, R *& A);
    void CSR(I *& IA, I *& JA, R *& A);
    void CSR();
    void CSC();
    void CSC(I *& JA, I *& IA, R *& A);
    static int addstateLU(int U) { return U>0 ? 4 : 5; };

    size_t SortLU(int U);
    size_t CSC_U(I *& JA, I *& IA, R *& A);
    size_t CSR_L(I *& IA, I *& JA, R *& A);
    void Buildp(I nn,I * IA,int type_m,size_t nnzz=0);
    Pair Row(I ii) { return RoworCol(ii, true /*!trans*/);}
    Pair Col(I jj) { return RoworCol(jj, false /*trans*/);}
    Pair RoworCol(I ii,bool row);
    void checksize(size_t nn,size_t mm=0) const
    { mm= mm ? mm : nn; assert( (nn ==this->n) && (mm==this->m));}
    R* addMatMul(R *x,R*Ax,bool Transpose,I sx=1,I sAx=1) const;
    R* addMatTransMul(R *x,R*Ax) const { return addMatMul(x,Ax,true); }
    R* addMatMul(R *x,R*Ax) const { return addMatMul(x,Ax,false);}
    R trace () const;
    double FrobeniusNorm() const;
    double norm1() const;
    double norminfty() const;
    bool sym() const {return half;}

    int typemat() const { return int(half)*VirtualMatrix<int,R>::TS_SYM ;}
    void SetBC(char *wbc,double ttgv);


    void addMap(R coef,std::map< pair<I,I>, R> &mij,bool trans=false,I ii00=0,I jj00=0,bool cnj=false,double threshold=0.);
    bool addMatTo(R coef,HashMatrix<I,R> & mij,bool trans=false,I ii00=0,I jj00=0,bool cnj=false,double threshold=0.,const bool keepSym=false) ;


    VirtualMatrix<I,R>  & operator +=(MatriceElementaire<R> & me) ;
     ostream& dump (ostream&f)  const { return f<<*this;}
    void SetBC(I ii,double ttgv) { diag(ii)=ttgv;};

    double gettgv(I * pntgv=0,double ratio=1e6) const ;
    bool GetReDoNumerics() const { bool b=re_do_numerics; re_do_numerics=0;return b;}
    bool GetReDoSymbolic() const { bool b=re_do_symbolic; re_do_symbolic=0;return  b;}


    HashMatrix<I, R> *toMatriceMorse(bool transpose=false,bool copy=false) const {ffassert(0); return 0;}
    double psor(KN_<R> & x,const  KN_<R> & gmin,const  KN_<R> & gmax , double omega) {ffassert(0); };

    void UnHalf();
    void Half() {resize(this->n,this->m,nnz,-1,true);}
    void RemoveHalf(int cas,double tol=-1) ;

    void setsdp(bool sym,bool dp); // set of unset to sym / defpos or not

    virtual bool ChecknbLine  (I n) const {return this->n==n;}
    virtual bool ChecknbColumn  (I m) const {return this->m==m;}
    static double CPUsecond() {
        return (double)clock()/CLOCKS_PER_SEC;
    }
};
// 0 good , -1 delete, ...
template<class I,class R> int GoodPtrHashMatrix(const HashMatrix<I,R> *p ) {
    if( p==0) return 1;
    if( p->N != p->n) return -2;
    if( p->M != p->m) return -3;
    if (p->nnz ==-1234567802) return  -4;
    if( p->i && p->j && p->aij ) return 0;
    return -5;
}
template<class I,class R> void CheckPtrHashMatrix(const HashMatrix<I,R> *p,const char * where )
{
    int gm=GoodPtrHashMatrix(p);
    if( gm !=0)
    {
        if(gm <0)
            cout << " n = " << p->n << " == " << p->N
                  << " , m= " << p->m << " "<< p->M
            << " nzz "<< p->nnz << endl;
        cerr << " Fatal Error " << where << "  invalide HashMatrix Ptr "<< gm << " "<< p << endl;
        ffassert(0);
    }
}

// END OF CLASS HashMatrix

template<class I,class R>
inline  size_t HashMatrix<I,R>::simpleinsert(I ii, I jj,Hash &h)
{
    state=unsorted;
    re_do_numerics=1;
    re_do_symbolic=1;

    type_state=type_HM;
    if(nnz==nnzmax) {
        Increaze();
        h = hash(ii,jj);
    }
    i[nnz] = ii;
    j[nnz] = jj;
    aij[nnz] = R();
    next[nnz]=head[h];
    head[h]=nnz;
    return nnz++;
}
template<class I,class R>
inline     R   * HashMatrix<I,R>::pij(I ii,I jj)  const
{
    re_do_numerics=1;
    Hash h = hash(ii,jj);
    size_t k = find(ii,jj,h);
    return k==empty ? 0 : aij+k;
}

template<class I,class R>
inline    R   *HashMatrix<I,R>::npij(I ii,I jj)   // with add if no term ii,jj
{
    re_do_numerics=1;
    Hash h = hash(ii,jj);
    size_t k = find(ii,jj,h);
    if(k==empty)
    {
        k =simpleinsert(ii,jj,h);
        aij[k]=0;
    }
    return aij+k;
}

//


template<class I,class RA,class RB=RA,class RAB=RA>
void AddMul(HashMatrix<I,RAB> &AB,HashMatrix<I,RA> &A, HashMatrix<I,RB> &B,bool ta=false,bool tb=false,R c=R(1))
{
    int An= A.n, Am =A.m;
     int Bn= B.n, Bm =B.m;
    if(ta) swap(An,Am);
    if(tb) swap(Bn,Bm);
    bool tcb = (std::is_same<RB,complex<double> >::value|| std::is_same<RB,complex<float> >::value ) && tb;
    bool tca = (std::is_same<RA,complex<double> >::value|| std::is_same<RA,complex<float> >::value ) && ta;
    AB.checksize(An,Bm);
    ffassert(Am == Bn);
    //  need A col sort  , b row sort
    if( tb)
      B.CSC(); // sort by COL nd build p.
    else
      B.CSR(); // sort by row... and build p.
    int * Bj = tb ? B.i : B.j;
    int * Bi = tb ? B.j : B.i;
    for(size_t l=0; l< A.nnz;++l)
    {
        I i=A.i[l],j=A.j[l];
        RA aij=A.aij[l];
        if(ta)  swap(i,j);
        if(tca) aij=HashMatrix<I,RA>::conj(aij);

        for(size_t ll=B.p[j]; ll<  B.p[j+1] ;++ll)
        {
            I k = Bj[ll];
            if(verbosity>1000000000) cout << " *** " << i<< " " << " " << k << " : " << j << " : "
                  << ll << " " << B.i[ll] <<" " << B.j[ll]<< " ::  " << A.aij[l]*B.aij[ll] <<endl;
            assert(j == Bi[ll]);
            RB bjk = tcb ? HashMatrix<I,RB>::conj(B.aij[ll]) : B.aij[ll];

            AB(i,k) += c* aij*bjk;
        }
    }
}

template<class I,class R>
std::ostream & operator<<(std::ostream & f,  const HashMatrix<I,R> &A)
{
    int p20=20;
    long pold= f.precision();
    if( pold > 20) p20= (int) pold;
    if(A.type_state==HashMatrix<I,R>::type_CSR)
    {
    using namespace std;
    f << "# Sparse Matrix (Morse)  " << &A << endl;
    f << "# first line: n m (is symmetic) nnz \n";
    f << "# after for each nonzero coefficient:   i j a_ij where (i,j) \\in  {1,...,n}x{1,...,m} \n";

    f << A.n << " " << A.m << " " << A.half << "  " << A.nnz <<endl;
    I k=A.p[0];

    for (I i=0;i<A.n;i++)
    {
        I ke=A.p[i+1];
        for (;k<ke;k++)
            f << setw(9) << i+1 << ' ' << setw(9) << A.j[k]+1 << ' ' << setprecision( p20) << A.aij[k]<< '\n' ;

    }
    }
    else
    {
        f << "#  HashMatrix Matrix (COO) "<< &A  << endl;
        f << "#    n       m        nnz     half     fortran   state  \n";
        f << A.n << " " << A.m << " " << A.nnz << " "<< A.half << " " << A.fortran
          << " " <<  A.state<< " " << A.type_state<< " " << endl;
        for(size_t k=0; k < A.nnz; ++k)
            f <<  setw(10) <<  A.i[k] << setw(10)  << A.j[k] << ' '<<  setprecision( p20)  << A.aij[k] << endl;
    }
    f.precision(pold);
    return f;
}

template<class R>
tuple<int,int,bool> BuildCombMat(HashMatrix<int,R> & mij,const list<tuple<R,VirtualMatrix<int,R>*,bool> >  &lM,bool trans,int ii00,int jj00,bool cnj)
{

    typedef typename list<tuple<R,VirtualMatrix<int,R> *,bool> >::const_iterator lconst_iterator;

    lconst_iterator begin=lM.begin();
    lconst_iterator end=lM.end();
    lconst_iterator i;


    int n=0,m=0;
    bool sym=true;
    for(i=begin;i!=end&&sym;i++)
    {
        if(get<1>(*i))// M == 0 => zero matrix
        {
            VirtualMatrix<int,R>& M=*get<1>(*i);
            if(!M.sym())
                sym = false;
        }
    }
    int iter=0;
    for(i=begin;i!=end;i++)
    {
        if(get<1>(*i)) // M == 0 => zero matrix
        {
            VirtualMatrix<int,R> & M=*get<1>(*i);
            bool transpose = get<2>(*i) !=  trans;
            bool conjuge = get<2>(*i) !=  cnj;

            ffassert( &M);
            R coef= get<0>(*i);
            if(verbosity>99)
                cout << "                "<< iter++<< " BuildCombMat + " << coef << "*" << &M << " " << sym << "  t = " << transpose << " " <<  get<2>(*i) << endl;
           { if(transpose) {m=max(m,M.n); n=max(n,M.m);} else{n=max(M.n,n); m=max(M.m,m);}}

            M.addMatTo(coef,mij,transpose,ii00,jj00,conjuge,0.0,sym);
        }
    }

    //V4 return new   MatriceMorseOld<R>(n,m,mij,sym);
    return make_tuple(n,m,sym);
}

template<class R>
tuple<int,int,bool> nmCombMat(const list<tuple<R,VirtualMatrix<int,R>*,bool> >  &lM,bool trans,int ii00,int jj00,bool cnj=false)
{

    typedef typename list<tuple<R,VirtualMatrix<int,R> *,bool> >::const_iterator lconst_iterator;

    lconst_iterator begin=lM.begin();
    lconst_iterator end=lM.end();
    lconst_iterator i;


    int n=0,m=0;
    bool sym=true;
    for(i=begin;i!=end&&sym;i++++)
    {
        if(std::get<1>(*i)) // M == 0 => zero matrix
        {
            VirtualMatrix<int,R>& M=*std::get<1>(*i);
            if(!M.sym())
                sym = false;
        }
    }

    for(i=begin;i!=end;i++++)
    {
        if(std::get<1>(*i)) // M == 0 => zero matrix
        {
            VirtualMatrix<int,R>& M=*std::get<1>(*i);
            bool transpose = std::get<2>(*i) !=  trans;
            ffassert( &M);
            R coef=std::get<0>(*i);
            if(verbosity>99)
                cout << "                BuildCombMat + " << coef << "*" << &M << " " << sym << "  t = " << transpose << " " << std::get<2>(*i) << endl;
            { if(transpose) {m=max(m,M.n); n=max(n,M.m);} else{n=max(M.n,n); m=max(M.m,m);}}

        }
    }

    return make_tuple(n,m,sym);
}

template<class R>
HashMatrix<int,R>* BuildCombMat(const list<tuple<R,VirtualMatrix<int,R>*,bool> >  &lM,bool trans=false,int ii00=0,int jj00=0)
{

    auto nmsym=nmCombMat(lM,trans,ii00,jj00);
    int n = std::get<0>(nmsym), m =std::get<1>(nmsym);
    bool half= std::get<2>(nmsym);
    HashMatrix<int,R> *  mij= new HashMatrix<int,R>(n,m,0,half);
    nmsym=BuildCombMat(*mij,lM,trans,ii00,jj00,trans);// remember trans => conj

    return mij; // V4 mij;

}

template<class I>
static   inline pair<I,I> ij_mat(bool trans,I ii00,I jj00,I i,I j) {
    // warning trans sub  matrix and not the block.
    return trans ? make_pair<I,I>(j+ii00,i+jj00)
    :  make_pair<I,I>(i+ii00,j+jj00) ; }



template<class I,class R>    template<typename T,typename TT>
void HashMatrix<I,R>::HMcopy( T *dst,const TT *from, size_t nn, T(*ff)(TT))
{
    for(size_t i=0; i< nn; ++i)
        dst[i]= ff(from[i]);
}
template<class I,class R>    template<typename T,typename TT>
void HashMatrix<I,R>::HMcopy( T *dst,const TT *from, size_t nn)
{
    for(size_t i=0; i< nn; ++i)
        dst[i]= (T) from[i];
}
template<class I,class R> template<class II,class RR> HashMatrix<I,R> &
HashMatrix<I,R>::operator+=(const HashMatrix<II,RR>& A )
{
    Addto(this,&A,cast_funct);
}


template<class I,class R> template<class II,class RR>
void HashMatrix<I,R>::set(II nn,II mm,bool hhalf,size_t nnnz, II *ii, II*jj, RR *aa,int f77,R(*ff)(RR))
{
    clear();
    this->n=nn;
    this->m=mm;
    this->N=nn;
    this->M=mm;
    fortran=f77;
    half=hhalf;
    Increaze(nnnz);
    nnz=nnnz;

    HMcopy(i,ii,nnnz);
    HMcopy(j,jj,nnnz);
    HMcopy(aij,aa,nnnz,ff);
    ReHash();
}
template<class I,class R> template<class II>
void HashMatrix<I,R>::set(II nn,II mm,bool hhalf,size_t nnnz, II *ii, II*jj, R *aa,int f77)
{
    clear();
    this->n=nn;
    this->m=mm;
    this->N=nn;
    this->M=mm;
    fortran=f77;
    half=hhalf;
    Increaze(nnnz);
    nnz=nnnz;

    HMcopy(i,ii,nnnz);
    HMcopy(j,jj,nnnz);
    HMcopy(aij,aa,nnnz);
    ReHash();
}
#endif
