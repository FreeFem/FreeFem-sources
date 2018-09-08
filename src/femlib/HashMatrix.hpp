#ifndef __HashMatrix_HPP__
#define __HashMatrix_HPP__

#include <cstddef>
#include <algorithm>
#include <map>
#include <iostream>
#include <cassert>
#include <complex>
#include <cstdint>
#include "RNM.hpp"
#include "RefCounter.hpp"
using std::max;
using std::min;
using std::cout;
using std::endl;
using std::pair;
using std::swap;
using std::complex;

extern  long  verbosity  ;

template<class Z>
inline uint64_t  roll64(Z y,int r){uint64_t x=y; r %= 64; return (x<<r) || (x << (64-r)) ;}

#include "VirtualMatrix.hpp"

template<class TypeIndex,class TypeScalaire>
class HashMatrix : public VirtualMatrix<TypeIndex,TypeScalaire>
{
    
    static double conj(double x){return x;}
    static float conj(float  x){return x;}
    static complex<double> conj(const complex<double> &x){return std::conj(x);}
    static complex<float> conj(const complex<float> &x){return std::conj(x);}

    static void  HeapSort(TypeIndex *ii,TypeIndex *jj, TypeScalaire *aij,long n)
    {
        long l,j,r,i;
        I criti,critj;
        R crita;
#define HM__criteq(i) criti=ii[i],critj=jj[i], crita=aij[i]
#define HM__eqcrit(i) ii[i]=criti,jj[i]=critj, aij[i]=crita
#define HM__eqij(i,j) ii[i]=ii[j],jj[i]=jj[j], aij[i]=aij[j]
#define HM__cmpij(i,j) ( ii[i] != ii[j] ? ii[i] < ii[j] : jj[i] < jj[j])
#define HM__cmpcrit(j) ( criti != ii[j] ? criti < ii[j] : critj < jj[j])
        
        if( n <= 1) return;
        l = n/2 + 1;
        r = n;
        while (1) { // label 2
            if(l <= 1 ) { // label 20
                HM__criteq(r-1);//crit = c[r];
                HM__eqij(r-1,0),--r; //c[r--] = c[1];
                if ( r == 1 ) { HM__eqcrit(0); /*c[1]=crit;*/ return;}
            } else  --l,  HM__criteq(l-1);// crit = c[--l];
            j=l;
            while (1) {// label 4
                i=j;
                j=2*j;
                if  (j>r) {HM__eqcrit(i-1);/*c[i]=crit;*/break;} // L8 -> G2
                if ((j<r) && (HM__cmpij(j-1,j) )) j++; // L5
                if (HM__cmpcrit(j-1)) HM__eqij(i-1,j-1);//c[i]=c[j]; // L6+1 G4
                else {HM__eqcrit(i-1); /*c[i]=crit;*/break;} //L8 -> G2
            }
        }
#undef HM__criteq
#undef HM__eqcrit
#undef HM__eqij
#undef HM__cmpij
#undef HM__cmpcrit
    }
    
public:
    typedef TypeIndex I;
    typedef TypeScalaire R;
    typedef uint64_t uniquecodeInt;
    static const int type_HM=0,type_COO=1, type_CSR=2,type_CSC=3;
    static const int  unsorted=0, sorted_ij=type_CSR,sorted_ji=type_CSC;
    typedef size_t iterator;
    typedef size_t Hash;
    typedef pair<size_t,size_t> Pair;
    size_t nnz,nnzmax,nhash,nbcollision,nbfind;
    I * i,*j;
    I *p;
    R * aij;
    size_t * head;
    size_t * next;
    // bool trans ; // bad idea at this level  F.H
    bool half; //
    int state,type_state;
    size_t nbsort;
    I sizep;
    int lock;
    int fortran; // index start a one ..
    int re_do_numerics,re_do_symbolic;
    static const  I empty= (I) -1;
    /*
     FOR compatibilite with MatriceMorse
     */
    I  nbcoef;
    bool  symetrique;
    R * a;
    I * lg;
    I * cl;

    int NbCoef() const {return  nnz;}
    void setcoef(const KN_<R> & x){ffassert(x.N()==nnz);KN_<R>(this->aij,nnz) = x;}
    void getcoef( KN_<R> & x){ffassert(x.N()==nnz);x =KN_<R>(this->aij,nnz);}
    
    void setdiag(const KN_<R> & d)
    {
        for(int ii=0; ii<this->n; ++ii)
            diag(ii)=d[ii];
    }
    void getdiag( KN_<R> & d) const
    {
        for(int ii=0; ii<this->n; ++ii)
            d[ii]=diag(ii);
    }

    R pscal(R *x,R *y,I sx=1,I sy=1)
    {
        R s=0;
        if( fortran)
        {
            x -= fortran;
            y -= fortran;
        }
        if(half)
          for (I k=0; k<nnz;++i)
          {
            if(i[k] != j[k] )
                s += conj(x[i[k]*sx])*aij[k]*y[j[k]*sy] + conj(x[j[k]*sx]*aij[k]) *y[i[k]*sy] ;
            else
                s += conj(x[i[k]*sx])*aij[k]*y[j[k]*sy];
          }
        else
         for (I k=0; k<nnz;++i)
            s += conj(x[i[k]*sx])*aij[k]*y[j[k]*sy];
        return s;
    }
    R pscal(const KN_<R> & x,const KN_<R> & y) { return pscal(x,y,x.step,y.step);}
    // produit scalaire
    
    // end add
   // void  resize(int n,int m) ; // add march 2009 ...

    /*
    class VirtualSolver :public RefCounter {
        friend class MatriceMorse;
        virtual void Solver(const MatriceMorse<R> &a,KN_<R> &x,const KN_<R> &b) const  =0;
    };
*/
    
    bool ismorse;
    void SetMorse(){
        CSR();
        nbcoef=nnz;
        symetrique= half;
        lg =p;
        aij= a;
        cl = j;
        
    }
    void UnSetMorse(){
        a=0;
        cl=0;
        lg=0;
        
    }
 
    uniquecodeInt CodeIJ() const  {
        uniquecodeInt code=this->n;
        code ^=roll64(this->n-this->m,24);
        code ^=roll64(nnz,48);
        for(I k=0,kk=0; k< nnz;++k)
        {
           code^=roll64(i[k],++kk);
           code^=roll64(j[k],++kk);
        }
    }
    
    HashMatrix(I nn,I mm=-1,size_t nnnz=0)
    :  VirtualMatrix<I,R>(nn,mm),  nnz(0),nnzmax(0),nhash(0),nbcollision(0),nbfind(0),i(0),j(0),p(0),aij(0),
    head(0), next(0),
    // trans(false),
    half(false), state(unsorted),type_state(type_HM),
    nbsort(0),sizep(0),lock(0), fortran(0) ,
    re_do_numerics(0),re_do_symbolic(0)
    {
        Increaze(nnnz);
    }
    
    HashMatrix(const HashMatrix& A)
    :  VirtualMatrix<I,R> (A), nnz(0),nnzmax(0),nhash(0),nbcollision(0),nbfind(0),
    i(0),j(0),aij(0)  ,
    head(0), next(0),
    //trans(false),
    half(false),  state(unsorted),type_state(type_HM),nbsort(0),
    sizep(0),lock(0), fortran(0),re_do_numerics(0),re_do_symbolic(0)
    {
        Increaze(A.nnz);
        operator=(A);
    }
    
    void CheckUnLock(const char * cmm)
    {
        if( lock)
        {
            std::cerr << " Sorry Forbidder operation ( "
            << cmm<< " ) matix is lock " << endl;
            assert(0);
        }
    }
    void setp(I sp)
    {
        if(sp == 0)
        {
            delete [] p;
            p=0;
        }
        
        else if( sp != sizep)
        {
            delete [] p;
            p = new I[sp];
            
        }
        if(p )
            for(size_t ii=0; ii<sp; ++ii)
                p[ii]=empty;
        sizep=sp;
    }
    void resize(int nn,int mm=0)
    {
        mm= mm ? mm : nn;
        if( mm>this->m ) this->m=mm;
        if (nn>this->n ) this->n = nn;
        if( nn == this->n && mm == this->m) return ;
        int kk=0;
        for(size_t k=0; k <nnz ;++k)
            if( i[k] < this->n && j[k] < this->m)
            {
                i[kk] = i[k];
                j[kk] = j[k];
                aij[kk] = aij[k];
                ++kk;
            }
        bool resize = kk*1.4 < nnz;
        nnz=kk;
        if(resize) Increaze();
        else ReHash();
        state= unsorted;
        type_state=type_COO;
    }
    
    void clear()
    {
        re_do_numerics=1;
        re_do_symbolic=1;
        half=false;
        lock=0;
        fortran = 0;
        state=unsorted;
        type_state = type_HM;
        if(nnz)
        {
            nnz=0;
            nbcollision=0;
            for (size_t h=0;h<nhash;++h)
                head[h]=empty;
        }
    }
    
    Hash hash(size_t ii,size_t jj)
    {
       // if(trans) swap(ii,jj);
        return ( (ii-fortran)+ (jj-fortran)*this->n )%nhash;
    }
    
    void setfortran(int yes)
    {
        
        int shift = fortran - yes;
        if( shift == 0) return ;
        CheckUnLock("setfortran");
        for( int k = 0; k<nnz ; ++ k)
        {
            i[k] += shift;
            j[k] += shift;
        }
        if (  type_CSC == type_state)
            for (int k=0; k<= this->m; ++k)
                p[k] += shift;
        else  if (  type_CSR == type_state)
            for (int k=0; k<= this->n; ++k)
                p[k] += shift;
        fortran = yes ;
        // do the shift
        return ;
    }
    
    size_t find(size_t ii,size_t jj)
    {
        return find(ii,jj, hash(ii,jj));
    }
    
    size_t  find(I ii,I jj,Hash h)
    {
       // if(trans) swap(ii,jj);
        nbfind++;
        for (size_t k=head[h];k!=empty;k=next[k])
        {
            ++nbcollision;
            if( ii== i[k] && jj==j[k] ) return k;
        }
        return empty;
    }
    
    iterator end()  { return iterator(this,nnz);}
    iterator begin(){ return iterator(this,0);}
    
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
    
    template<typename T>
    static void resize(T *&t,size_t no,size_t nn)
    {
        T * tt= new T[nn];
        for(size_t i=0; i< no; ++i)
            tt[i]= t[i];
        if(t) delete [] t;
        t=tt;
    }
    
    void Increaze(size_t nnznew=0)
    {
        if(!nnznew) nnznew = max(this->n,this->m)*10;
        size_t nzzx = max(size_t(nnzmax*1.2),nnznew);
        size_t nh = max(size_t(nhash*1.1),(size_t)max(this->n,this->m));
        resize(i,nnz,nzzx);
        resize(j,nnz,nzzx);
        resize(aij,nnz,nzzx);
        resize(next,0,nzzx);
        resize(head,0,nh);
        nnzmax=nzzx;
        nhash = nh;
        ReHash();
        state=unsorted;
    }
    
    void ReHash()
    {
        for(int h=0;h<nhash;++h)
            head[h]= empty;
        for(int k=0;k<nnz;++k)
        {
            Hash h= hash(i[k],j[k]);
            next[k]=head[h];
            head[h]=k;
        }
    }
    size_t size() const { return nnz;}
    size_t simpleinsert(I ii, I jj,Hash &h)
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
    R   *pij(I ii,I jj)  const
    {
        Hash h = hash(ii,jj);
        size_t k = find(ii,jj,h);
        return k==empty ? 0 : aij+k;
    }
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
        re_do_numerics=1;
        Hash h = hash(ii,jj);
        size_t k = find(ii,jj,h);
        if(k==empty) 
		{   
			k =simpleinsert(ii,jj,h);
			return aij[k]=R();
		}
        else return aij[k];
    }
    
    R    operator()(pair<I,I> ij)  const {return operator()(ij.first,ij.second);}
    R &  operator()(pair<I,I> ij)  {return operator()(ij.first,ij.second);}
    ~HashMatrix()
    {
        if(nbfind && verbosity>4)
            cout << "    ~HashMatrix:   Mean collision in hash: " << (double) nbcollision/ nbfind << endl;
        delete [] i;
        delete [] j;
        delete [] aij;
        delete [] next;
        delete [] head;
        if(p) delete [] p;
    }
    
    void Sortij()
    {
        
        if( state != sorted_ij)
        {
            CheckUnLock("Sortij");
            HeapSort(i,j,aij,nnz);
            ReHash();
            nbsort++;
            state = sorted_ij;
        }
    }
    
    void Sortji()
    {
        if( state != sorted_ji)
        {
            CheckUnLock("Sortji");
            HeapSort(j,i,aij,nnz);
            ReHash();
            state = sorted_ji;
            nbsort++;
        }
    }
    
    void operator=(const HashMatrix&A)
    {
        if( this == & A) return;
        clear();
        this->n=A.n;
        this->m=A.m;
        //nnz=0;
        nhash=A.nhash;// ????????
        //  nbcollision=0;
        nbfind=0;
        //trans=A.trans;
        fortran=A.fortran;
        //state=unsorted;
        //type_state=type_HM;
        nbsort=0;
        sizep=0;
        lock=0;
        re_do_numerics=1;
        re_do_symbolic=1;
        
        operator+=(A);
    }
    void operator+=(const HashMatrix&A)
    {
        this->n = max(this->n,A.n);
        this->m = max(this->m,A.m);
        
        if( this == & A)
            for(int k=0; k < nnz; ++k)
                aij[k] += aij[k];
        else
        {   if ( this->half!= this->half ) MATERROR(1,"+= on diff type of mat of type HashMatrix ");
            size_t nnzm= min(nnz,A.size()) ;
            for(size_t k=0;k<nnzm;++k)
                if( (i[k]==A.i[k]) && (j[k]==A.j[k]))
                    aij[k] += A.aij[k]; // Optim...
                else
                    operator()(A.i[k],A.j[k]) += A.aij[k];
            for(size_t  k=nnzm+1;k<A.size();++k)
                operator()(A.i[k],A.j[k]) += A.aij[k];
        }
    }
    
    void operator+=( pair<R,const HashMatrix*> pp)
    {
        const R  aa = pp.first;
        const HashMatrix &A= *pp.second;
        this->n = max(this->n,A.n);
        this->m = max(this->m,A.m);
        
        if( this == & A)
            for(int k=0; k < nnz; ++k)
                aij[k] += aa*aij[k];
        else
        {   if ( this->half!= this->half ) MATERROR(1,"+= on diff type of mat of type HashMatrix ");
            size_t nnzm= min(nnz,A.size()) ;
            for(size_t k=0;k<nnzm;++k)
                if( (i[k]==A.i[k]) && (j[k]==A.j[k]))
                    aij[k] += aa*A.aij[k]; // Optim...
                else
                    operator()(A.i[k],A.j[k]) += A.aij[k];
            for(size_t  k=nnzm+1;k<A.size();++k)
                operator()(A.i[k],A.j[k]) += aa*A.aij[k];
        }
    }
    
    void operator*=(R v)
    {
        re_do_numerics=1;
        for(int k=0; k < nnz; ++k)
            aij[k] *= v;
        
    }
    void operator=(R v)
    {
        re_do_numerics=1;
        for(int k=0; k < nnz; ++k)
            aij[k] = v;
        
    }
    
    void COO()
    {
        Sortij();
        setp(0);
        type_state=type_COO;
    }
    
    void COO(I *& IA, I *& IJ, R *& A)
    {
        Sortij();
        IA=i;
        IJ=j;
        A=aij;
        setp(0);
        type_state=type_COO;
    }
    
    void CSR(I *& IA, I *& JA, R *& A)
    {
        Sortij();
        Buildp(this->n,i,type_CSR);
        IA=p;
        JA=j;
        A=aij;
        type_state=type_CSR;
    }
    void CSR()
    {
        Sortij();
        Buildp(this->n,i,type_CSR);
        type_state=type_CSR;
    }
    
    void CSC()
    {
        Sortji();
        Buildp(this->m,j,type_CSC);
        type_state=type_CSC;
    }
    
    void CSC(I *& JA, I *& IA, R *& A)
    {
        
        Sortji();
        Buildp(this->m,j,type_CSC);
        JA=p;
        IA=i;
        A=aij;
        type_state=type_CSC;
    }
    static int addstateLU(int U) { return U>0 ? 4 : 5; }
    
    size_t SortLU(int U)
    {   // put U or L at the top and return the nnz term associated
        size_t nnzu =0,nnzd=0;
        for(int k=0; k<nnz; ++ k)
            if( i[k]<j[k]) ++nnzu; // 1 2 in U => i <= j
            else if (i[k]==j[k]) ++nnzd;
        size_t nnzl=nnz-nnzu-nnzd;
        cout << "SortLU " << U << ":s " << nnzl << " "<< nnzd << " " << nnzu << endl;
        size_t kR=0;
        if(U>0)
        {
            size_t kU=0L,k;
            
            for( k=0; k<nnz; ++ k)
                if( i[k]<=j[k])
                {
                    std::swap(aij[k],aij[kU]);
                    std::swap(i[k],i[kU]);
                    std::swap(j[k],j[kU++]);
                } // 1 2 in U => i <= j
            // verif
            cout << " SortLU "<< kU << " " << nnzu+nnzd << " " << k << endl;

            for( k=0; k<kU ; ++k)
                assert(i[k]<=j[k]);
            for( k=kU; k<nnz ; ++k)
                assert(i[k]>j[k]);
            
            kR = kU;
            assert(kU == nnzu+nnzd);
        }
        else
        {
            size_t kL=0L,k;
            
            for( k=0; k<nnz; ++ k)
                if( i[k]>=j[k])
                {
                    std::swap(aij[k],aij[kL]);
                    std::swap(i[k],i[kL]);
                    std::swap(j[k],j[kL++]);
                }
            
            kR = kL;
            

        }
        ReHash();
        nbsort++;
        state += addstateLU(U) ;
        return kR;
        
    }
    size_t CSC_U(I *& JA, I *& IA, R *& A)
    {
     
        if(type_CSC!=type_state)
            HeapSort(j,i,aij,nnz);
        state=type_CSC;
        size_t nnzu=SortLU(1);
        Buildp(this->m,j,type_CSC+addstateLU(1),nnzu);
     
        JA=p;
        IA=i;
        A=aij;
        return nnzu;
    }
    size_t CSR_L(I *& IA, I *& JA, R *& A)
    {
        if(type_CSR!=type_state)
           HeapSort(i,j,aij,nnz);
        state=type_CSR;
        size_t nnzl=SortLU(-1);
        Buildp(this->m,i,type_CSR+addstateLU(-1),nnzl);
        
        IA=p;
        JA=j;
        A=aij;
        return nnzl;
    }
    
    void Buildp(I nn,I * IA,int type_m,size_t nnzz=0)
    {
        if(nnzz==0) nnzz=nnz;
        if(type_m != type_state)
        {
            assert( state==type_m);
            
            setp(nn+1);
            I k=I(nnzz);
            assert( (nnzz-k) ==0);
            //int shift =  fortran;
            do
            {
                --k;
                p[IA[k]-fortran] = k+fortran;
            }
            while( k!=0);
            
            p[0]=0;
            p[nn] = (I) nnzz;
            // empty line
            for(size_t ii=0;ii<nn;++ii)
                if(p[ii+1]==empty)
                    p[ii+1]=p[ii];
            
            assert(p[nn]-nnzz ==0 );
            
            type_state=type_m;
            if(verbosity>10)
            {
            cout << "Buildp  p=" << nn<< " type_m =" << type_m << endl;
            for(int ii=0;ii<=nn;++ii)
                cout << ii << " " << p[ii] << endl;
            cout << " end Buildp\n";
            }
        }
        assert(p);
    }
    
    Pair Row(I ii) { return RoworCol(ii, true /*!trans*/);}
    Pair Col(I jj) { return RoworCol(jj, false /*trans*/);}
    
    Pair RoworCol(I ii,bool row)
    {
        size_t k0=0,k1=nnz-1,k=nnz-1;
        int * pp=0;
        if(row)
        {
            if(state != sorted_ij) Sortij();
            pp= i;
        }
        else
        {
            if(state != sorted_ji) Sortji();
            pp = j;
        }
        assert(nbsort < this->n+100);
        
        while(pp[k0] < ii)
        {
            if( k-k0 <=1) break;
            size_t kk = (k0+k)/2;
            //cout << kk << " " << k0 << endl;
            if(pp[kk]<ii)
                k0=kk;
            else
            {
                if(pp[kk]>ii) k1=kk;
                k=kk;
            }
        }
        if(pp[k0]<ii)
            ++k0;
        k=k0;
        while(ii<pp[k1])
        {
            if( ( k1-k ) <= 1) break;
            size_t kk = (k1+k)/2;
            //cout << kk << " " << k1 << endl;
            if(pp[kk]>ii)
                k1=kk;
            else
                k=kk;
        }
        if(pp[k1]>ii)
            --k1;
        //cout << ii << "k0 =" << k0 << " k1 " << k1 << endl;
        Pair  ret(k0,k1);
        return ret;
    }

    void checksize(size_t nn,size_t mm=0) const
    { mm= mm ? mm : nn; assert( (nn ==this->n) && (mm==this->m));}
    
    R* addMatMul(R *x,R*Ax,bool Transpose) const {
        I *ii=ii,*jj=j;
        R *aa=aij;
     //   if(Transpose != trans) {std::swap(ii,jj);}
        if(Transpose ) {std::swap(ii,jj);}
        if(fortran) {aa++;}
        if( half)
          for(int k=0; k<nnz;++k)
            {
                Ax[ii[k]] += aa[k]*x[jj[k]];
                if( ii[k] != jj[k]) Ax[jj[k]] += conj(aa[k])*x[ii[k]];
            }
        else
            for(int k=0; k<nnz;++k)
                Ax[ii[k]] += aa[k]*x[jj[k]];
        
         return Ax;}
    
    R* addMatTransMul(R *x,R*Ax) const {
        return addMatMul(x,Ax,true);
    }
 
    R* addMatMul(R *x,R*Ax) const {
        return addMatMul(x,Ax,false);
    }
    R trace () const {
        R t=R();
        for(int ii=0; ii<this->n; ++ii)
          t+= diag(ii);
        return t;
    }
  R FrobeniusNorm() const
    {
        R s=0;
        for(I k=0; k<nnz;++k)
            s += norm2(aij[k]);
        s = sqrt(s);
        return  s; 
    }
  R norm1() const
    {
        R s=0;
        for(I k=0; k<nnz;++k)
            s += abs(aij[k]);
        return  s;
    }
    R norminfty() const
    {
        R s=abs(aij[0]);
        for(I k=0; k<nnz;++k)
            s = std::max(abs(aij[k]),s);
        return  s;
    }
    bool sym() const {return half;}
 
    
    void SetBC(char *wbc,double tgv)
    {
        if ( this->n != this->m) MATERROR(1,"SetBC on none square matrix  ?????");
        for(I ii=0; ii< this->n; ++ii)
            if(tgv<0)
            {
                
                if( wbc[ii] )
                {
                    for (I k=lg[ii];k<lg[ii+1]; ++k)
                        if( cl[k]==ii)
                            a[k] = 1.;
                        else
                            a[k]=0;// put the line to Zero.
                }
                else if(tgv < -1.999)
                    for (I k=lg[ii];k<lg[ii+1]; ++k)
                    {
                        I jj = cl[k];
                        if( wbc[jj] ) a[k]=0;//
                    }
            }
            else  if( wbc[ii] ) { // tgv >= 0
                operator()(i,i)=tgv;
            }

    }

    bool  addMatTo(R coef,map<I,R> & A,bool trans=false,int ii00=0,int jj00=0,bool cnj=false,double threshold=0.,const bool keepSym=false)
    {
        
        ffassert(0); // to do ..
    }
    
    VirtualMatrix<I,R>  & operator +=(MatriceElementaire<R> & ){ ffassert(0); };
     void operator=(const R & v) =0; // Mise a zero
     ostream& dump (ostream&f)  const { return f<<*this;}
    void SetBC(int i,double tgv) { diag(i)=tgv;};
     HashMatrix<I, R> *toMatriceMorse(bool transpose=false,bool copy=false) const {return 0;} // not
     bool addMatTo(R coef,std::map< pair<int,int>, R> &mij,bool trans=false,int ii00=0,int jj00=0,bool cnj=false,double threshold=0.,const bool keepSym=false)=0;
    double psor(KN_<R> & x,const  KN_<R> & gmin,const  KN_<R> & gmax , double omega) {ffassert(0); };
    
    
};


// END OF CLASS HashMatrix
template<class I,class R>
void AddMul(HashMatrix<I,R> &A, HashMatrix<I,R> &B, HashMatrix<I,R> &AB,R c=R(1))
{
    AB.checksize(A.n,B.m);
    assert(A.m == B.n);
    //  need A col sort  , b row sort
    B.CSR(); // sort by row... and build p.
    
    cout << " B = " <<  B << " \n; \n\n;";
    for(size_t l=0; l< A.nnz;++l)
    {
        I i=A.i[l],j=A.j[l];
        //typename HashMatrix<I,R>::Pair rj(B.Row(j));
        //      for(size_t ll=rj.first; ll< rj.second ;++ll)
        
        for(size_t ll=B.p[j]; ll<  B.p[j+1] ;++ll)
        {
            cout << j << " : " << ll << " " << B.i[ll] <<" " << B.j[ll] <<endl;
            assert(j == B.i[ll]);
            I k = B.j[ll];
            AB(i,k) += c* A.aij[l]*B.aij[ll];
        }
    }
}

template<class I,class R>
std::ostream & operator<<(std::ostream & coutt,  const HashMatrix<I,R> &A)
{
    using namespace std;
    coutt << A.n << " x " << A.m << " " <<  " nnz =" << A.nnz << " states: "<<  A.state<< " " << A.type_state<< " " << endl;
    for(int k=0; k < A.nnz; ++k)
        coutt << k<< " : "<<  A.i[k] << " " << A.j[k] << " "<< A.aij[k] << endl;
    if(A.p)
        for(int k=0; k < A.sizep; ++k)
            coutt << k << " "<< A.p[k] << endl;
    return coutt;
}
#endif
