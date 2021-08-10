
#include "HashMatrix.hpp"
//   Warning the instation will be don a the add Of the File
// F. hecht
// version v0 ...

int WhichMatrix(istream & f)
{
    string line;
    while ( isspace(f.peek()))
        f.get();
    if  ( f.peek() =='#' )
    {
        line="";
        while ( f.good()  )
        {
            char c=f.get();
            if(c=='\n' || c=='\r') { break;}
            line += c;
        }
        if( line.find("(Morse)") != std::string::npos)
            return 2; // morse
        
        else  if( line.find("(COO)") != std::string::npos)
            return 3; // morse
        
        return 0;
    }
    return 0;
}


template<class I,class R>
void  HashMatrix<I,R>::HeapSort(I *ii,I *jj, R *aij,long n)
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
            HM__eqij(r-1,0);
            --r; //c[r--] = c[1];
            if ( r == 1 ) { HM__eqcrit(0); /*c[1]=crit;*/ return;}
        } else  {--l;  HM__criteq(l-1);}// crit = c[--l];
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

template<class I,class R>
void HashMatrix<I,R>::setdiag(const KN_<R> & d)
{
    for(int ii=0; ii<this->n; ++ii)
        diag(ii)=d[ii];
}

template<class I,class R>
void HashMatrix<I,R>::getdiag( KN_<R> & d) const
{
    for(int ii=0; ii<this->n; ++ii)
        d[ii]=diag(ii);
}

template<class I,class R>
R HashMatrix<I,R>::pscal(R *x,R *y,I sx,I sy)
{
    R s=0;
    if( fortran)
    {
        x -= fortran;
        y -= fortran;
    }
    if(half > 0)
        for (size_t k=0; k<nnz;++i)
        {
            if(i[k] != j[k] )
                s += conj(x[i[k]*sx])*aij[k]*y[j[k]*sy] + conj(x[j[k]*sx])* (half == 1?conj(aij[k]):aij[k]) *y[i[k]*sy] ;
            else
                s += conj(x[i[k]*sx])*aij[k]*y[j[k]*sy];
        }
    else
        for (size_t k=0; k<nnz;++i)
            s += conj(x[i[k]*sx])*aij[k]*y[j[k]*sy];
    return s;
}


template<class I,class R>
typename   HashMatrix<I,R>::uniquecodeInt HashMatrix<I,R>::CodeIJ() const  {
    uniquecodeInt code=this->n;
    code ^=roll64(this->n,0);
    code ^=roll64(this->m,24);
    code ^=roll64(nnz, 48);
    for(size_t k=0,kk=0; k< nnz;++k)
    {
        code^=roll64(i[k],++kk);
        code^=roll64(j[k],++kk);
    }
    return code;
}

template<class I,class R>
HashMatrix<I,R>::HashMatrix(I nn,I mm,I nnnz,int halff)
:  VirtualMatrix<I,R>(nn,mm),  nnz(0),nnzmax(0),nhash(0),nbcollision(0),nbfind(0),matmulcpu(0.),i(0),j(0),p(0),aij(0),
head(0), next(0),
// trans(false),
half(halff), state(unsorted),type_state(type_HM),
nbsort(0),sizep(0),lock(0), fortran(0) ,
re_do_numerics(0),re_do_symbolic(0)
{
    Increaze(nnnz);
}

template<class I,class R>
HashMatrix<I,R>::HashMatrix(istream & f,int cas):
VirtualMatrix<I,R>(0,0),  nnz(0),nnzmax(0),nhash(0),nbcollision(0),nbfind(0),matmulcpu(0.),i(0),j(0),p(0),aij(0),
head(0), next(0),
// trans(false),
half(0), state(unsorted),type_state(type_HM),
nbsort(0),sizep(0),lock(0), fortran(0) ,
re_do_numerics(0),re_do_symbolic(0)
{
    if(cas ==-1) cas=  WhichMatrix(f)  ;
    // eat lines with #
    string line;
    int k=0;
    while ( isspace(f.peek()))
        f.get();
    while ( f.peek() =='#' )
    {
        line="";
        while ( f.good()  )
        {
            char c=f.get();
            if(c=='\n' || c=='\r') { break;}
            line += c;
        }
        if( f.peek()=='\n' || f.peek()=='\r') f.get();
        if(verbosity>9)
            cout << "   --  Matrx: "<< k << " :"   << line << " cas " << cas << endl;
        k++;
    }
    if(cas== 2)
    {
        I rn,rm;
        int rsymetrique;
        size_t rnbcoef;
        ffassert(f.good() );
        f >> rn >> rm >> rsymetrique >>rnbcoef;
        
        if(verbosity>3)
            cout << "     -- Read Mat: " <<  this->n << " x " <<  this->m << " sym : " << rsymetrique << " nnz=" << rnbcoef <<endl;
        ffassert(f.good() && rn>=0 && rm>=0 && rnbcoef>=0 );
        resize(rn,rm,rnbcoef);
        half =rsymetrique;
        I ii,jj;
        
        R aaij;
        
        for (size_t kk =0;kk<rnbcoef; ++kk)
        {
            f >> ii >> jj >> aaij;
            ffassert(f.good() );
            operator()(ii-1,jj-1) = aaij;
        }
    }
    else if(cas== 3)
    {
        I nn,mm,f77, sstate, tstate;
        int hhalf;
        size_t nnnz;
        f >> nn >> mm >> nnnz >> hhalf >> f77 >> sstate >> tstate ;
        ffassert(f.good() && nn>0 && mm>0 && nnnz>0  );
        resize(nn,mm,nnnz);
        half =hhalf;
        for(size_t kk=0; kk < nnnz; ++kk)
        {
            I ii,jj;
            R aaij;
            f >>  ii >> jj >> aaij ;
            ffassert(f.good() );
            operator()(ii-f77,jj-f77) = aaij;
        }
        
    }
    else {
        cerr << " Unknown matrix type" << endl << endl;
        ffassert(0);
    }
    
    
}

template<class I,class R>
HashMatrix<I,R>::HashMatrix(KNM_<R> F,double  threshold)
: VirtualMatrix<I,R>(F.N(),F.M()),  nnz(0),nnzmax(0),nhash(0),nbcollision(0),nbfind(0),matmulcpu(0.),i(0),j(0),p(0),aij(0),
head(0), next(0),
half(0), state(unsorted),type_state(type_HM),
nbsort(0),sizep(0),lock(0), fortran(0) ,
re_do_numerics(0),re_do_symbolic(0),tgv(0), ntgv(0)

{
    Increaze();
    for(int ii=0; ii <F.N(); ++ii)
        for(int jj=0; jj <F.M(); ++jj)
        {
            R Fiijj = F(ii,jj);
            if(abs(Fiijj) > threshold)
                operator()(ii,jj)=Fiijj;
        }
}
//  PB ....
template<class I,class R>
HashMatrix<I,R>::HashMatrix(I nn,I mm,bool Half)
: VirtualMatrix<I,R>(nn,mm),  nnz(0),nnzmax(0),nhash(0),nbcollision(0),nbfind(0),matmulcpu(0.),i(0),j(0),p(0),aij(0),
head(0), next(0),
half(Half), state(unsorted),type_state(type_HM),
nbsort(0),sizep(0),lock(0), fortran(0) ,
re_do_numerics(0),re_do_symbolic(0)

{
    Increaze(max(nn,mm)*4);
}

template<class I,class R>
HashMatrix<I,R>::HashMatrix(const HashMatrix& A)
:  VirtualMatrix<I,R> (A),  nnz(0),nnzmax(0),nhash(0),nbcollision(0),nbfind(0),matmulcpu(0.),i(0),j(0),p(0),aij(0),
head(0), next(0),
half(A.half), state(unsorted),type_state(type_HM),
nbsort(0),sizep(0),lock(0), fortran(0) ,
re_do_numerics(0),re_do_symbolic(0)
{
    Increaze(A.nnz);
    operator=(A);
}

template<class I,class R> template<class II>
HashMatrix<I,R>::HashMatrix(const HashMatrix<II,R>& A)
:  VirtualMatrix<I,R> (A.n,A.m),  nnz(0),nnzmax(0),nhash(0),nbcollision(0),nbfind(0),matmulcpu(0.),i(0),j(0),p(0),aij(0),
head(0), next(0),
half(A.half), state(unsorted),type_state(type_HM),
nbsort(0),sizep(0),lock(0), fortran(0) ,
re_do_numerics(0),re_do_symbolic(0)
{
    HashMatrix<I,R>::set<II>(A.n,A.m,A.half,A.nnz,A.i,A.j,A.aij,fortran);
}

template<class I,class R> template<class II,class RR>
HashMatrix<I,R>::HashMatrix(const HashMatrix<II,RR>& A, R (*ff)(RR) )
:  VirtualMatrix<I,R> (A.n,A.m),  nnz(0),nnzmax(0),nhash(0),nbcollision(0),nbfind(0),matmulcpu(0.),i(0),j(0),p(0),aij(0),
head(0), next(0),
half(A.half), state(unsorted),type_state(type_HM),
nbsort(0),sizep(0),lock(0), fortran(0) ,
re_do_numerics(0),re_do_symbolic(0)
{
     HashMatrix<I,R>::set<II,RR>(A.n,A.m,A.half,A.nnz,A.i,A.j,A.aij,fortran,ff);
}

template<class I,class R>
HashMatrix<I,R>::HashMatrix(I nn,const R *diag)
:  VirtualMatrix<I,R> (nn,nn),  nnz(0),nnzmax(0),nhash(0),nbcollision(0),nbfind(0),matmulcpu(0.),i(0),j(0),p(0),aij(0),
head(0), next(0),
half(0), state(unsorted),type_state(type_HM),
nbsort(0),sizep(0),lock(0), fortran(0) ,
re_do_numerics(0),re_do_symbolic(0)
{
    Increaze(nn);
    for(int k=0;k<nn;++k)
        operator()(k,k) = diag[k];
    
}


template<class I,class R>
void HashMatrix<I,R>::setp(I sp)
{
 
    if(sp == 0)
    {
        delete [] p;
        p=0;
    }
    
    else if( sp != sizep)
    {
        if(p) delete [] p;
        p = new I[sp];
        
    }
    if(p )
        for(I ii=0; ii<sp; ++ii)
            p[ii]=-1;
     if(verbosity>999) cout << "  HashMatrix:: setp "<< this << " sp= " << sp << " p  "<< p << " / " << (sp != sizep) << endl;
    ffassert( (sp==0) ==  (p==0)  );
    sizep=sp;
    
}
template<class I,class R>
void  HashMatrix<I,R>::RemoveHalf(int cas,double tol)
{
    
    size_t kk=0;
    if( cas <=0) // L cas
    {
    for(size_t k=0; k <nnz ;++k)
        if(  abs(aij[k])> tol && ( i[k] >= j[k] ) )
        {
            i[kk] = i[k];
            j[kk] = j[k];
            aij[kk] = aij[k];
            ++kk;
        }
    }
    else
    {
     for(size_t k=0; k <nnz ;++k)
        if(  abs(aij[k])> tol && ( i[k] <= j[k] ) )
        {
            i[kk] = i[k];
            j[kk] = j[k];
            aij[kk] = aij[k];
            ++kk;
        }
    }
    half= cas==0;
    size_t nrt = nnz - kk;
    if (verbosity>2)
        cout << " RemoveHalf " << cas << " tol "  << tol << " , remove term : " << nrt << " half =" << half << endl;
    nnz = kk;
    if(nrt)  Increaze(kk);
    else ReHash();
    state= unsorted;
    type_state=type_COO;

}

template<class I,class R>
void  HashMatrix<I,R>::resize(I nn, I mm,size_t nnnz, double tol , int sym )
{
    
    mm= mm ? mm : nn;
    if( nn == this->n && mm == this->m && nnz == nnnz && sym == half && tol <0) return ;
    this->m=mm;
    this->n=nn;
    this->N=nn;
    this->M=mm;
    R mxt =0;
    size_t kk=0;
    if( (nn>0) && (mm >0))
    for(size_t k=0; k <nnz ;++k)
    {
        double t =abs(aij[k]);
        if( i[k] < this->n && j[k] < this->m && t > tol && (!sym || i[k] <= j[k] ) )
        {
            i[kk] = i[k];
            j[kk] = j[k];
            aij[kk] = aij[k];
            ++kk;
        }
    }
    if(verbosity>9 ) cout << "HashMatrix<I,R>::resize: new nnz  = " << kk << " " << nnz << " " << " sym " << sym << " "<< tol << " " << nn << " " << mm << endl;
    half = (half ? half : sym);
    nnnz = max(nnnz,kk);
    bool rresize = (nnnz < 0.8*nnz) || ((nnnz > 1.2*nnz)) ;
    nnz = kk; // forget after value ... 
    if(rresize) Increaze(nnnz);
    else ReHash();
    state= unsorted;
    type_state=type_COO;
}
template<class I,class R>
void HashMatrix<I,R>::SymmetrizePattern()
{
    ffassert( this->n==this->m);
    for(size_t k=0; k <nnz ;++k)
        npij( j[k], i[k]);
    
}

template<class I,class R>
void HashMatrix<I,R>::clear()
{
    nbsort=0;
    re_do_numerics=1;
    re_do_symbolic=1;
    half=0;
    lock=0;
    fortran = 0;
    state=unsorted;
    type_state = type_HM;
    setp(0); // unset
    if(nnz)
    {
        nnz=0;
        nbcollision=0;
        for (size_t h=0;h<nhash;++h)
            head[h]=empty;
    }
    
}

template<class I,class R>
void HashMatrix<I,R>::setfortran(int yes)
{
// trop casse geule
    int shift =  yes-fortran;
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




template<class I,class R>  template<typename T>
void HashMatrix<I,R>::HMresize(T *&t,size_t no,size_t nn)
{
    if( no != nn || t==0)
    {
    T * tt= new T[nn];
    if(t)
        for(size_t i=0; i< no; ++i)
           tt[i]= t[i];
    if(t) delete [] t;
    t=tt;
    }
}


template<class I,class R>   bool  HashMatrix<I,R>::do2Triangular(bool lower)
{
    ffassert(half); //
    size_t nc=0;
    if( lower )
       for(size_t k=0; k< nnz; ++k)
       {
           if(  (i[k] < j[k]) ) // Upper
               ++nc,swap(i[k],j[k]);
       }
    else
        for(size_t k=0; k< nnz; ++k)
        {
            if(  (i[k] > j[k]) ) // Lower
                ++nc,swap(i[k],j[k]);
        }
    if( nc)
    {
    ReHash();
    state=unsorted;
    }
    return  nc;
}
template<class I,class R>   int HashMatrix<I,R>::IsTrianglulare() const
{
    size_t nU, nL=0,nD=0;
    for(size_t k=0; k< nnz; ++k)
    {
        if( i[k] < j[k] ) ++nU;// 0,10
        else if( i[k] == j[k] ) ++nD;// 10,10
        else ++nL;// 10,0 
    }
    // 3 => no sym
    // 1 nU==0 =>   Trian Lower (no upper term)
    // 2 nL==0 =>   Trian Upper (no upper term)
    return 2*!nL +  !nU ;
}
template<class I,class R>
void HashMatrix<I,R>::dotranspose()
{
    swap(i,j);
    swap(this->n,this->m);
    swap(this->N,this->M);
    conj(aij,this->nnz);
    ReHash();
    state=unsorted;
}
template<class I,class R>
void HashMatrix<I,R>::Renumbering(I nn,I mm,KN_<I> ii,KN_<I> jj)
{
    //Do  B_ii(i),jj(j) : A_i,j
    size_t kk=0;
    for(size_t k=0; k<nnz; ++k)
    {
        I i1=ii[i[k]],j1=jj[i[k]];
        if( i1 >=0 && j1 >=0 && i1 <nn && j1 < mm) // coef existe
        {
            i[kk] =i1;
            j[kk] =j1;
            aij[kk++]=aij[k];
        }
    }
    nnz = kk;
    this->n=nn;
    this->m=mm;
    this->N=nn;
    this->M=mm;
    state=unsorted;
    
    RemoveDoubleij(kk); // remove double term 
  

    
}
template<class I,class R>
void HashMatrix<I,R>::RemoveDoubleij(int kkk)
{
    nnz=kkk;
    COO();
    I ip=-1,jp=-1;
    long kk=-1;
    for(size_t k=0; k<nnz; ++k)
    {
        if( ip != i[k] && jp != j[k] )
        {  // new term
            ++kk;
            i[kk] =ip=i[k];
            j[kk] =jp=j[k];
            aij[kk]=aij[k];
        }
        else
            aij[kk]=+aij[k];
    }
    if(verbosity>4)
        cout << " HashMatrix::RemoveDoubleij  remove "<< nnz - kk << " coef "<< endl;
    Increaze(kk,kk);
}
template<class I,class R>
void HashMatrix<I,R>::RenumberingInv(KN_<I> II,KN_<I> JJ)
{
       //Do  B_(i),(j) : A_ii(i),jj(j)
    I n = this->n, m = this->m;
    I nn= II.N(), mm= JJ.N();
    const I minus1 =-1;
    KN<I> ii(n,minus1), jj(m,minus1);
    // build inversion
    int notinjection=0;
    // build invertion ..
    for( I l=0; l<nn; ++l)
    {
        I IIl = II[l];
        if( IIl >= 0 && IIl < n)
        {
            if(ii[IIl]>=0 ) notinjection =1;
          ii[IIl]=l;
        }
    }
    for( int l=0; l<mm; ++l)
    {
        I IIl = JJ[l];
        if( IIl >= 0 && IIl < m)
        {
            if(jj[IIl]>=0 ) notinjection |=2 ;
            jj[IIl]=l;
        }
    }
    if( !notinjection)
    {
        cerr << " HashMatrix<I,R>::Renumbering not injection " <<notinjection <<endl;
        ffassert(0); // to do 
    }
    Renumbering(nn,mm,ii,jj);
}

template<class I,class R>
void HashMatrix<I,R>::Increaze(size_t nnzxnew,size_t newnnz)
{
    size_t mnx =(size_t)max(this->n,this->m);
    if(newnnz==0) newnnz = nnz;
    else nnz=newnnz; 
    if( nnzxnew==0) {
        nnzxnew = max(max(mnx*2,newnnz),(size_t) (nnzmax*1.2) );
    }
    
    size_t nzzx =  nnzxnew;
    double nnzl = min(max(1.,double(nzzx)/mnx),50.);// bof bof !!!! FH
    size_t nh = mnx*nnzl;
    if(verbosity>999) cout << "HMresize "<< nzzx << " " <<nnzxnew << " mpirank" << mpirank << " " << this << endl; 
    HMresize(i,nnz,nzzx);
    HMresize(j,nnz,nzzx);
    HMresize(aij,nnz,nzzx);
    HMresize(next,0,nzzx);
    HMresize(head,0,nh);
    nnzmax=nzzx;
    nhash = nh;
    ReHash();
    state=unsorted;
}

template<class I,class R>
void HashMatrix<I,R>::ReHash()
{
    
    for(size_t h=0;h<nhash;++h)
        head[h]= empty;
    for(size_t k=0;k<nnz;++k)
    {
        Hash h= hash(i[k],j[k]);
        next[k]=head[h];
        head[h]=k;
    }
}


template<class I,class R>
HashMatrix<I,R>::~HashMatrix()
{
    if(nbfind && verbosity>4)
        cout << "    ~HashMatrix:   Mean collision in hash: " << (double) nbcollision/ nbfind << " "
        << this << " rank: "<< mpirank << " matmul " << matmulcpu << "s" << endl;
    
    delete [] i;
    delete [] j;
    delete [] aij;
    delete [] next;
    delete [] head;
    if(p) delete [] p;
    type_state=type_isdeleted;// Mark Matrix is type_isdeleted
    i =0;
    j=0;
    aij=0;
    next=0;
    head =0;
    this->n=-1234567890;//  stupide number
    this->m=-1234567801;//  stupide number
    this->nnz=-1234567802;
    // because the solver is some time deleted after ...
}

template<class I,class R>
void HashMatrix<I,R>::Sortij()
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

template<class I,class R>
void HashMatrix<I,R>::Sortji()
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

template<class I,class R>
void HashMatrix<I,R>::set(I nn,I mm,int hhalf,size_t nnnz, I *ii, I*jj, R *aa,int f77,int tcsr)
{
//    tcsr >0 => CSR ii pointer size nn+1, jj col size nnnz
//    tcrs <0 = CSC  jj pointer size mm+1, ii row size nnnz
//    tcsr == 0 => COO

    clear();
    this->n=nn;
    this->m=mm;
    this->N=nn;
    this->M=mm;
    fortran=f77;
    half=hhalf;
    Increaze(nnnz);
    nnz=nnnz;
    if(tcsr >=0) //  input CSR
      HMcopy(j,jj,nnnz);// copy of col
    if(tcsr <=0) //  input CSC
      HMcopy(i,ii,nnnz);// copy of  row
    if( tcsr>0)//  input CSR
        for(I ip=0; ip< nn; ++ip)
            for(I k=ii[ip]; k<ii[ip+1]; ++k)
                i[k-f77]=ip+f77;
    else if( tcsr<0) //  input CSC
        for(I jp=0; jp< mm; ++jp)
            for(I k=jj[jp]; k<jj[jp+1]; ++k)
                j[k-f77]=jp+f77;

    HMcopy(aij,aa,nnnz);
    ReHash();
}


template<class I,class R>
HashMatrix<I,R> & HashMatrix<I,R>::operator=(const HashMatrix &A)
{
    if( (const void*)  this == (const void*) & A) return *this;
    set(A.n,A.m,A.half,A.nnz,A.i,A.j,A.aij,A.fortran);
    return *this;
}

template<class I,class R>
HashMatrix<I,R> & HashMatrix<I,R>::operator+=(const HashMatrix &A)
{
    Add(&A);
    return *this;
}

template<class I,class R>
HashMatrix<I,R> & HashMatrix<I,R>::operator-=(const HashMatrix &A)
{
    Add(&A,R(-1.));
    return *this;
}

template<class I,class R>
void HashMatrix<I,R>::Add(const HashMatrix<I,R> *PA,R coef,bool trans, I ii00,I jj00)
{
    const HashMatrix<I,R> &A=*PA;
    I nn=A.n,mm=A.m;
    I *ii=A.i;
    I *jj=A.j;
    size_t Annz= A.nnz;  ;
    if( trans ) swap(nn,mm),swap(ii,jj);
    nn = max(this->n,nn);
    mm = max(this->m,mm);
    resize(nn,mm);
    if( (const void*) this == (const void*) & A && ii00==0 && jj00==0 && !trans )
        for(I k=0; k < nnz; ++k)
            aij[k] += coef*aij[k];
    else
    {   ffassert((const void*) this != (const void*) & A); // not the same matrix
        if ( this->half!= A.half ) MATERROR(1,"+= on diff type of mat of type HashMatrix ");
        
        for(size_t k=0;k<Annz;++k)
            operator()(ii[k]+ii00,jj[k]+jj00) += coef*A.aij[k];
    }
}


template<class I,class R,class K>
void Addto(HashMatrix<I,R> *P0, const HashMatrix<I,K> *PA,R (*f)(K) ,bool trans, I ii00,I jj00)
{
    const HashMatrix<I,K> &A=*PA;
    I nn=A.n,mm=A.m;
    I *ii=A.i;
    I *jj=A.j;
    size_t Annz= A.nnz;  ;
    if( trans ) swap(nn,mm),swap(ii,jj);
    nn = max(P0->n,nn);
    mm = max(P0->m,mm);
    P0->resize(nn,mm);
    
    {
        if( (const void*) P0 == (const void*) & A && ii00==0 && jj00==0 && !trans )
            for(I k=0; k < P0->nnz; ++k)
                P0->aij[k] += f(A.aij[k]);
        else
        {   ffassert((const void*) P0 != (const void*) & A); // not the same matrix
            if ( P0->half!= A.half ) MATERROR(1,"+= on diff type of mat of type HashMatrix ");
            
            for(size_t k=0;k<Annz;++k)
                P0->operator()(ii[k]+ii00,jj[k]+jj00) += f(A.aij[k]);
        }
    }
}




template<class I,class R>
void HashMatrix<I,R>::operator*=(R v)
{
    re_do_numerics=1;
    for(int k=0; k < nnz; ++k)
        aij[k] *= v;
    
}

template<class I,class R>
void HashMatrix<I,R>::operator=(const R & v)
{
    re_do_numerics=1;
    for(int k=0; k < nnz; ++k)
        aij[k] = v;
    
}

template<class I,class R>
void HashMatrix<I,R>::HM()
{
    re_do_numerics++;
    re_do_symbolic++;
    setp(0);
    type_state=type_HM;
    
}

template<class I,class R>
void HashMatrix<I,R>::COO()
{
     if(verbosity>999) cout << " HashMatrix:: COO "<< this->n << " x " << this->m << " "<< p << " / "<< sizep << endl;
    Sortij();
    setp(0);
    type_state=type_COO;
}


template<class I,class R>
void HashMatrix<I,R>::COO(I *& IA, I *& IJ, R *& A)
{
    Sortij();
    IA=i;
    IJ=j;
    A=aij;
    setp(0);
    type_state=type_COO;
}


template<class I,class R>
void HashMatrix<I,R>::CSR(I *& IA, I *& JA, R *& A)
{
    Sortij();
    Buildp(this->n,i,type_CSR);
    IA=p;
    JA=j;
    A=aij;
    type_state=type_CSR;
}

template<class I,class R>
void HashMatrix<I,R>::CSR()
{
     if(verbosity>999) cout << " HashMatrix:: csr()  " << state << " " << sorted_ij<< endl;
    Sortij();
    Buildp(this->n,i,type_CSR);
    type_state=type_CSR;
}


template<class I,class R>
void HashMatrix<I,R>::CSC()
{
    Sortji();
    Buildp(this->m,j,type_CSC);
    type_state=type_CSC;
}


template<class I,class R>
void HashMatrix<I,R>::CSC(I *& JA, I *& IA, R *& A)
{
    
    Sortji();
    Buildp(this->m,j,type_CSC);
    JA=p;
    IA=i;
    A=aij;
    type_state=type_CSC;
}



template<class I,class R>
size_t HashMatrix<I,R>::SortLU(int U)
{   // put U or L at the top and return the nnz term associated
    if( half) {
        
    }
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

template<class I,class R>
size_t HashMatrix<I,R>::CSC_U(I *& JA, I *& IA, R *& A)
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

template<class I,class R>
size_t HashMatrix<I,R>::CSR_L(I *& IA, I *& JA, R *& A)
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


template<class I,class R>
void HashMatrix<I,R>::Buildp(I nn,I * IA,int type_m,size_t nnzz)
{
   if(verbosity>999) cout << " HashMatrix:: Buildp"<< this->n<< " x " << this->m << " " << nn << " " << IA<< " " << type_m << " " << nnzz << " / " << nnz << " / p=" << p <<endl;
    if(nnzz==0) nnzz=nnz;
    if(type_m != type_state)
    {
        
        assert( state==type_m);
        
        setp(nn+1);
         //int shift =  fortran;
        std::fill(p,p+nn+1,-1);
        p[nn] = fortran+ (I) nnzz;
        for( I k=I(nnzz)-1; k>=0 ; --k )
        {
            p[IA[k]-fortran] = k+fortran;
            ffassert( (IA[k]-fortran>=0 ) && (IA[k]-fortran<=nn));
        }
        p[nn]=nnzz+fortran;
        // remove empty row
        for(I ii=nn-1;ii>=0;--ii)
            if(p[ii]<0)//  empty row
                p[ii]=p[ii+1];
        
     
        
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



template<class I,class R>
typename HashMatrix<I,R>::Pair HashMatrix<I,R>::RoworCol(I ii,bool row)
{
    size_t k0=0,k1=nnz-1,k=nnz-1;
    I * pp=0;
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
    assert(nbsort < (size_t) this->n+100 );
    
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



template<class I,class R>
R* HashMatrix<I,R>::addMatMul(R *x,R*Ax,bool Transpose,I sx,I sAx) const {
    I *ii=i,*jj=j;
    R *aa=aij;
    double t0= CPUsecond();
    //   if(Transpose != trans) {std::swap(ii,jj);}
    const bool do_conj  = ((std::is_same<R,complex<double> >::value || std::is_same<R,complex<float> >::value )) && Transpose;
    if(Transpose ) {std::swap(ii,jj);}
    if(fortran) {aa++;}
    if(do_conj)// complex and transpose => hermitian
    {
    if(sx==1 && sAx==1)
    {
        if( half)
            for(size_t k=0; k<nnz;++k)
            {
                R aak=aa[k],caak=conj(aa[k]);
                Ax[ii[k]] += caak*x[jj[k]];
                if( ii[k] != jj[k])
                    Ax[jj[k]] += caak*x[ii[k]];
            }
        else
            for(size_t k=0; k<nnz;++k)
                Ax[ii[k]] += conj(aa[k])*x[jj[k]];
    }
    else
    if( half)
        for(size_t k=0; k<nnz;++k)
        {
            R aak=aa[k],caak=conj(aa[k]);
            Ax[ii[k]*sAx] += caak*x[jj[k]*sx];
            if( ii[k] != jj[k])
                Ax[jj[k]*sAx] += aak*x[ii[k]*sx];
        }
    else
        for(size_t k=0; k<nnz;++k)
            Ax[ii[k]*sAx] += conj(aa[k])*x[jj[k]*sx];
    }
    else // no transp or no complex =>  (no conj)
    {
    if(sx==1 && sAx==1)
    {
        if( half)
            for(size_t k=0; k<nnz;++k)
            {    
                Ax[ii[k]] += aa[k]*x[jj[k]];
                if( ii[k] != jj[k])
                    Ax[jj[k]] += aa[k] *x[ii[k]];
            }
        else
            for(size_t k=0; k<nnz;++k)
                Ax[ii[k]] += aa[k]*x[jj[k]];
    }
    else
    if( half)
        for(size_t k=0; k<nnz;++k)
        {
            Ax[ii[k]*sAx] += aa[k]*x[jj[k]*sx];
            if( ii[k] != jj[k])
                Ax[jj[k]*sAx] += aa[k]*x[ii[k]*sx];
        }
    else
        for(size_t k=0; k<nnz;++k)
            Ax[ii[k]*sAx] += aa[k]*x[jj[k]*sx];
    }
    matmulcpu+=  CPUsecond()-t0;
    return Ax;}


template<class I,class R>
R HashMatrix<I,R>::trace () const {
    R t=R();
    for(int ii=0; ii<this->n; ++ii)
        t+= diag(ii);
    return t;
}

template<class I,class R>
double HashMatrix<I,R>::FrobeniusNorm() const
{
    double s=0;
    for(size_t k=0; k<nnz;++k)
        s += norm(aij[k]);
    s = sqrt(s);
    return  s;
}

template<class I,class R>
double HashMatrix<I,R>::norm1() const
{
    double  s=0;
    for(size_t k=0; k<nnz;++k)
        s += abs(aij[k]);
    return  s;
}

template<class I,class R>
double HashMatrix<I,R>::norminfty() const
{
    double s=abs(aij[0]);
    for(size_t k=0; k<nnz;++k)
        s = std::max(abs(aij[k]),s);
    return  s;
}


template<class I,class R>
void HashMatrix<I,R>::SetBC(char *wbc,double ttgv)
{
    tgv = ttgv;
    ntgv =0;
    if(ttgv<0)
        CSR();
    if ( this->n != this->m) MATERROR(1,"SetBC on none square matrix  ?????");
    for(I ii=0; ii< this->n; ++ii)
        if( wbc[ii] )
        {
            ntgv++;
            if(ttgv<0)
            {
                
                if( wbc[ii] )
                {
                    for (I k=p[ii];k<p[ii+1]; ++k)
                        if( j[k]==ii )
                            aij[k] = (std::abs(ttgv+10.0) < 1.0e-10 ? 0.0 : 1.0);
                        else
                            aij[k]=0;// put the line to Zero.
                }
                
            }
            else
                operator()(ii,ii)=ttgv;
        }
    if(std::abs(ttgv+2.0) < 1.0e-10 || std::abs(ttgv+20.0) < 1.0e-10) //  remove also columm tgv == -2 .....
    {
        CSC();
        for(I jj=0; jj< this->n; ++jj)
            if( wbc[jj] ) {
                for (I k=p[jj];k<p[jj+1]; ++k)
                    if( i[k]!=jj || std::abs(ttgv+20.0) < 1.0e-10)
                        aij[k]=0;//
            }
    }
    
    
}



template<class I,class R>
void HashMatrix<I,R>::addMap(R coef,std::map< pair<I,I>, R> &mij,bool trans,I ii00,I jj00,bool cnj,double threshold)
{
    for (auto pm=mij.begin();  pm != mij.end(); ++pm)
    {
        I ii = pm->first.first+ii00, jj= pm->first.second+jj00;
        R cmij = coef* pm->second;
        
        if(trans) swap(ii,jj);
        if( abs(cmij)  > threshold )
        {
          if(cnj) cmij = conj(cmij);
          operator()(ii,jj) += cmij;
        }
    }
}

template<class I,class R>
bool HashMatrix<I,R>::addMatTo(R coef,HashMatrix<I,R> & mij,bool trans,I ii00,I jj00,bool cnj,double threshold,const bool keepSym)
{
    //  add a mij + = coef * [(this)^trans^cnj , 
    double eps0=max(numeric_limits<double>::min(),threshold);
    
    if (half)
    {
        for( size_t kk= 0; kk<nnz ;++kk)
        {
            I ii=i[kk], jj=j[kk];
            R cij =  coef* ( cnj ? RNM::conj(aij[kk]) : aij[kk]);
            if(threshold==0 || std::norm(cij)>eps0) // remove for IPOPT april 2019 FH only if threshold >0
            {
                mij[ij_mat(trans,ii00,jj00,ii,jj)] += cij ;
                if (ii!=jj&&!keepSym)
                    mij[ij_mat(trans,ii00,jj00,jj,ii)] += cij;
            }
        }
    }
    else
    {
        for(size_t  kk= 0; kk<nnz ;++kk)
        {
            I ii=i[kk], jj=j[kk];
            R cij =  coef* ( cnj ? RNM::conj(aij[kk]) : aij[kk]);
            if(threshold==0 || std::norm(cij)>eps0) // / remove for IPOPT april 2019 FH
                mij[ij_mat(trans,ii00,jj00,ii,jj)] += cij;
        }
    }
    
    return keepSym;
}


template<class I,class R>
VirtualMatrix<I,R>  & HashMatrix<I,R>::operator +=(MatriceElementaire<R> & me) {
    //  R zero=R();
    int il,jl,i,j;
    int * mi=me.ni, *mj=me.nj;
    if ((this->n==0) && (this->m==0))
    {
        
        
        cout << "  -- Bug: HashMat  is empty let's build it" << endl;
        ffassert(0);
        
    }
    R * al = me.a;
    R * aij;
    switch (me.mtype) {
        case MatriceElementaire<R>::Full : ffassert(!half);
            for (il=0; il<me.n; ++il)  { i=mi[il];
                for ( jl=0; jl< me.m ; ++jl,++al)  {j=mj[jl];
                    aij = npij(i,j);
                    {
                        throwassert(aij);
                        *aij += *al;}}}
            break;
            
        case MatriceElementaire<R>::Symmetric : ffassert(half);
            for (il=0; il<me.n; ++il) {  i=mi[il] ;
                for (jl=0;jl< il+1 ; ++jl) { j=mj[jl];
                    aij =    (j<i) ? npij(i,j) : npij(j,i);
                    throwassert(aij);
                    *aij += *al++;}}
            break;
        default:
            cerr << "Big bug type MatriceElementaire unknown" << (int) me.mtype << endl;
            ErrorExec("Bi-ug in  Hashmat += MatriceElementaire",1);
            break;
    }
    return *this;
}


template<class I,class R>
void SetBC(I ii,double ttgv) { diag(ii)=ttgv;};



template<class I,class R>
double HashMatrix<I,R>::gettgv(I * pntgv,double ratio) const
{
    if( this->n != this->m) return 0; // no ttgv
    double ttgv =0, max1=0;
    I ntgv=0;
    for (I ii=0; ii<this->n;++ii)
    {
        R * p=pij(ii,ii);
        if (p)
        {
            
            double a=real(*p);
            if( a> ttgv )
            {
                max1=ttgv;
                ttgv=a;
                ntgv =1;
            }
            else if (a== ttgv)
                ++ntgv;
            else if( a >max1)
                max1 = a;
        }
        
    }
    if( max1*ratio> ttgv) // ttgv to small => no ttgv ....
        ttgv=0,ntgv=0; // no ttgv
    if(pntgv) *pntgv = ntgv;
    return ttgv;
    
}

template<class I,class R>
void HashMatrix<I,R>::setsdp(int sym,bool dp) // sym dpos para
{
    this->symetric=(sym > 0);
    this->positive_definite=dp;
    if( (half>0) != sym)
    {
        if( sym)
            Half(sym);//  remove half part
        else
            UnHalf();
    }
    else {
      half = sym;
    }
}

template<class I,class R>
void HashMatrix<I,R>::UnHalf()
{
    
    if (half==0) return;
    HM();
    size_t nnz0=nnz,err=0;
    for(int k=0; k<nnz0; ++k)
      if( i[k] > j[k] )
       {
        size_t ki=insert(j[k],i[k],(half == 1?conj(aij[k]):aij[k]));
        err += ki < nnz0;
        }
    else if ( i[k] < j[k] )
        err++;
    half = 0;
    if( err )
        cerr << " Try of unsymmetrize no half matrix (Bug) ... ????" <<endl;
    ffassert(err==0);
    
}

typedef double R;
typedef complex<R> C;
//  because for UMFPACK 64  because long are only 32 bits under windows 
#ifdef _WIN32
typedef long long  int64;
#else
typedef long int64;
#endif

template class HashMatrix<int,R>;
template class HashMatrix<int,C >;
template class HashMatrix<int64,R>;
template class HashMatrix<int64,C >;

template HashMatrix<int64,R>::HashMatrix(const HashMatrix<int,R> & );
template HashMatrix<int64,C>::HashMatrix(const HashMatrix<int,C> & );
template HashMatrix<int,C>::HashMatrix(const HashMatrix<int64,C> & );
template HashMatrix<int,R>::HashMatrix(const HashMatrix<int,C> & , R(*ff)(C));
//template HashMatrix<int,C>::HashMatrix(const HashMatrix<int,R> & );

//template HashMatrix<int,R> & HashMatrix<int,R>::operator=(const HashMatrix<int,C> & );
//template HashMatrix<int,C> & HashMatrix<int,C>::operator=(const HashMatrix<int,R> & );
//template HashMatrix<int,R> & HashMatrix<int,R>::operator+=(const HashMatrix<int,C> & );
//template HashMatrix<int,C> & HashMatrix<int,C>::operator+=(const HashMatrix<int,R> & );

template  void Addto<int,R,C>(HashMatrix<int,R> *P0, const HashMatrix<int,C> *PA,R (*f)(C) ,bool trans, int ii00,int jj00);
template  void Addto<int,C,R>(HashMatrix<int,C> *P0, const HashMatrix<int,R> *PA,C (*f)(R) ,bool trans, int ii00,int jj00);

template void HashMatrix<int,R>::set<int64>(int64 nn,int64 mm,int hhalf,size_t nnnz, int64 *ii, int64*jj, R *aa,int f77);
//template void HashMatrix<int,R>::set<int,R>(int nn,int mm,bool hhalf,size_t nnnz, int *ii, int*jj, R *aa,int f77,R(*ff)(R));
//template void HashMatrix<int,R>::set<long,R>(long nn,long mm,bool hhalf,size_t nnnz, long *ii, long *jj, R *aa,int f77,R(*ff)(R));
template void HashMatrix<int,C>::set<int,C>(int nn,int mm,int hhalf,size_t nnnz, int *ii, int*jj, C *aa,int f77,C(*ff)(C));
template void HashMatrix<int,R>::set<int,C>(int nn,int mm,int hhalf,size_t nnnz, int *ii, int*jj, C *aa,int f77,R(*ff)(C));

//template void HashMatrix<int,C>::set<int,R>(int nn,int mm,bool hhalf,size_t nnnz, int *ii, int*jj, R *aa,int f77,C(*ff)(R));
//template void HashMatrix<int,C>::set<long,C>(long nn,long mm,bool hhalf,size_t nnnz, long *ii, long *jj, C *aa,int f77,C(*ff)(C));


 void init_HashMatrix ()
{
    
}
// just to test
/*
static void tttt() {
    HashMatrix<int,R> AiR(10);
    HashMatrix<long,R> AlR(AiR);//,HashMatrix<long,R>::cast_funct);
 
    HashMatrix<long,C> AlC(10);
    HashMatrix<int,C> AiC(10);

}

*/
