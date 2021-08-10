#include <iostream>
#include <cmath>
#include "HashMatrix.hpp"

//#include "cholmod_function.h"
#include <complex>
enum FactorizationType {
    FactorizationNO=0,
    FactorizationCholesky=1,
    FactorizationCrout=2,
    FactorizationLU=3};

template <class Z,class R>
class SkyLineMatrix {
public:
    typedef HashMatrix<Z,R>  HMat;
    Z n,m;
    FactorizationType whichfac;
    /*----------------------------------------------------------------
     D[i] = A[ii]
     L[k] = A[ij]  j < i avec:   pL[i]<= k < pL[i+1] et j = pL[i+1]-k
     U[k] = A[ij]  i < j avec:   pU[j]<= k < pU[j+1] et i = pU[i+1]-k
     remarque pL = pU generalement
     si L = U => la matrice est symetrique
     -------------------------------------------------------------------
     */

    mutable R *L,*oL;  // lower  if oL == L => no delete
    mutable R *U,*oU;  // upper   ..
    mutable R *D,*oD;  // diagonal   ..
    Z *pL,*opL; // profile L
    Z *pU,*opU; // profile U
    mutable  FactorizationType typefac;
    int verb;
    FactorizationType typesolver;
    std::ostream& dump (std::ostream&) const ;
    SkyLineMatrix(HMat *A,Z *p,int tf,int verbb);
    
    SkyLineMatrix(const SkyLineMatrix& A,bool copy=false,int vrb=verbosity)
    : n(A.n),m(n),whichfac(A.whichfac),
    L(A.L),oL(L),      U(A.U),oU(U),     D(A.D),oD(A.D),
    pL(A.pL),opL(pL),  pU(A.pU),opU(pU),
    typefac(A.typefac),verb(vrb)
    { if(copy) docopy(); }
    
    SkyLineMatrix(Z NbOfDF,R* d,
                   R* u, Z * pu,
                   R* l, Z * pl,
                   bool copy,int vrb=verbosity,
                   FactorizationType tf=FactorizationNO)
    : n(NbOfDF),m(n),
      L(l),oL(l),       U(u),oU(u),    D(d),oD(d),
      pL(pl),opL(pl),   pU(pu),opU(pu),
      typefac(tf),verb(vrb)
    { if(copy) docopy(); }
    
 
    ~SkyLineMatrix()
    {
        if(D &&  D != oD) delete [] D;
        if(L &&  L != oL) delete [] L;
        if(U &&  U != oU && U !=L ) delete [] U;
        if(pL &&  pL != opL  ) delete [] pL;
        if(pU &&  pU != opU && pU != pL  ) delete [] pU;
    }
    
    const SkyLineMatrix t(bool cpy=false) const  {return SkyLineMatrix(this->n,D,L,pL,U,pU,cpy,verb,typefac);}
    const SkyLineMatrix lt(bool cpy=false) const {return SkyLineMatrix(this->n,0,L,pL,0,0,cpy,verb);}
    const SkyLineMatrix l(bool cpy=false) const  {return SkyLineMatrix(this->n,0,0,0,L,pL,cpy,verb);}
    const SkyLineMatrix d(bool cpy=false) const  {return SkyLineMatrix(this->n,D,0,0,0,0,cpy,verb);}
    const SkyLineMatrix ld(bool cpy=false) const {return SkyLineMatrix(this->n,D,0,0,L,pL,cpy,verb);}
    const SkyLineMatrix ldt(bool cpy=false) const {return SkyLineMatrix(this->n,D,L,pL,0,0,cpy,verb);}
    const SkyLineMatrix du(bool cpy=false) const {return SkyLineMatrix(this->n,D,U,pU,0,0,cpy,verb);}
    const SkyLineMatrix u(bool cpy=false) const  {return SkyLineMatrix(this->n,0,U,pU,0,0,cpy,verb);}
    const SkyLineMatrix ut(bool cpy=false) const {return SkyLineMatrix(this->n,0,0,0,U,pU,cpy,verb);}
    const SkyLineMatrix dut(bool cpy=false) const {return SkyLineMatrix(this->n,D,0,0,U,pU,cpy,verb);}

    Z size() const ;
    
    void cholesky(double = 1e-15) const ; //
    void crout(double = 1e-15) const ; //
    void LU(double = 1e-15)  const ; //
    void factorize(double eps= 1e-15) const
    {
        if(whichfac==FactorizationLU) LU(eps);
        else if (whichfac==FactorizationCholesky) cholesky(eps);
        else if (whichfac==FactorizationCrout) crout(eps);
        else MATERROR(1,"UNKNOWN FACTORIZATION type SkyLineMatrix");
    }

  
    
   static void SumpT(Z *pT,Z n)
    {
        Z i=0, s=0, pTi=0;
        for( i=0; i<n;++i)
        {
            pTi=pT[i];
            pT[i]=s;
            s+= pTi;
        }
        pT[n]=s;
        if(verbosity > 1)
            cout << " SumpT="<< s << endl;
    }
    R* solve(R*x,int trans=1) const;
    template <class RR> static RR* newcpy(RR *p,size_t n)
    {
        RR *c = new RR[n];
        for(size_t i=0; i<n; ++i)
            c[i]=p[i];
        return c;
    }

private:
    void docopy()
    {
        bool copy=true;
        if(copy&& D ) D= newcpy(oD,n);
        if(copy&& L ) L= newcpy(oL,pL[n]);
        if(copy&& U ) L= newcpy(oU,pU[n]);
        if(copy&& pL) pL= newcpy(opL,n+1);
        if(opL==opU)  pU=pL;
        else if(copy&& pU) pU= newcpy(opU,n+1);
    }
    // no copy
    void operator=(const SkyLineMatrix & A);
};


template <class Z,class R>
 SkyLineMatrix<Z,R>::SkyLineMatrix(HashMatrix<Z,R> *A,Z *p,int typfact,int verbb)
: n(A->n),m(n),whichfac((FactorizationType)typfact),verb(verbb),
L(0),oL(0),U(0),oU(0),D(0),oD(0),
pL(0),opL(0),pU(0),opU(0),
typefac(FactorizationNO)
{
    bool sym=whichfac==FactorizationCrout || whichfac==FactorizationCholesky;
   //  calcul de pL
    //    A_p(i),p(j) = PAP'_ij
    // p1[p[i]=i
    //    i,j,aij   =>> p1(i),p1(j), ajj
    std::vector<Z> p1(n);
    for(Z i=0; i<n; ++i)
    {
         p[i]=i;//ZZZZZ
        p1[p[i]]=i;// perm inver
    }
    pU=pL= new Z[n+1];
    std::fill(pL,pL+n,0);
    if(!sym) {
        pU= new Z[n+1];
        std::fill(pU,pU+n,0);
    }
   
    for(int k=0; k<A->nnz;++k)
    {
        Z i=p1[A->i[k]],j=p1[A->j[k]];//    i,j,aij   =>> p1(i),p1(j), ajj
        if(i>j) pL[i] = max(i-j,pL[i]);// i>j => L
        if(i<j) pU[j] = max(j-i,pU[j]);// i < j => U
    }
   
    if(pL !=  pU)
    {
         Z diff=0;
        for(int i=0; i<n;++i)
            diff += abs(pL[i]-pU[i]);
        if(verb > 2 || verbosity >9)
            cout << "  - SkyLineMatrix : diff pL/pU: "<< n<< " " << diff<< endl;
        if(2*diff <= n)
        {
            // Skyline symetric near symetric
            for(int i=0; i<n;++i)
                pL[i]= max(pL[i],pU[i]);
            delete [] pU;
            pU=pL;
        }
    }
    SumpT(pL,n);
    if(pU && pU != pL) SumpT(pU,n);
    if(sym)
    {
        L = new R[pL[n]];
        D = new R[n];
        fill(L,L+pL[n],R());
        fill(D,D+n,R());

        U = L;
        for(int k=0; k<A->nnz;++k)
        {
            Z i=p[A->i[k]],j=p[A->j[k]];
            if( i ==j)
                D[i]= A->aij[k];
            else if( i > j )
                L[pL[i+1]+j-i]=A->aij[k];
        }
    }
    else
    {
        L = new R[pL[n]];
        D = new R[n];
        U = new R[pU[n]];
        fill(L,L+pL[n],R());
        fill(D,D+n,R());
        fill(U,U+pU[n],R());

        for(int k=0; k<A->nnz;++k)
        {
            Z i=p1[A->i[k]],j=p1[A->j[k]];
            if( i ==j) D[i]= A->aij[k];
            else if( i > j ) L[pL[i+1]+j-i]=A->aij[k];
            else U[pU[j+1]-j+i]=A->aij[k];
        }
    }
    if(verbosity>4 || verb)
        cout << "  SkyLineMatrix: size pL/pU: "<< n<< " " << pL[n] << " " << pU[n] << " moy=" << pU[n]/(double) n << endl;
}

template <class Z,class R>
void SkyLineMatrix<Z,R>::cholesky(double eps) const {
    double eps2=eps*eps;
    R  *ij , *ii  , *ik , *jk , xii;
    Z i,j,k;
    if (L != U) MATERROR(1,"cholesky Skyline matrix non symmetric");
    U = 0; //
    typefac = FactorizationCholesky;
    if(verbosity>3 || verb >1 )
        cout << "  -- SkyLineMatrix Factorize/Cholesky   " << endl;
    if ( std::norm(D[0]) <= 1.0e-60)
        MATERROR(2,"cholesky SkyLine pivot ");
        
        D[0] = sqrt(D[0]);
    ij = L ; // pointeur sur le terme ij de la matrice avec j<i
    for (i=1;i<this->n;i++) // boucle sur les lignes
    { ii = L+pL[i+1]; // pointeur sur le terme fin de la ligne +1 =>  ij < ii;
        xii = D[i] ;
        for ( ; ij < ii ; ij++) // pour les j la ligne i
        { j = i -(ii - ij);
            k = max( j - (pL[j+1]-pL[j]) ,  i-(pL[i+1]-pL[i]) );
            ik =  ii - (i - k);
            jk =  L + pL[j+1] -(j - k);
            k = j - k ;
            R s= -*ij;
#ifdef WITHBLAS
            s += blas_sdot(k,ik,1,jk,1);
#else
            while(k--) s += *ik++ * *jk++;
#endif
            *ij =  -s/D[j] ;
            xii -= *ij * *ij ;
        }
        // cout << std::norm(xii) << " " << Max(eps2*std::norm(D[i]),1.0e-60) << " " << sqrt(xii) <<endl;
        if ( std::norm(xii) <= max(eps2*std::norm(D[i]),1.0e-60))
            MATERROR(3,"cholesky SkyLine pivot ");
            D[i] = sqrt(xii);
    }
}
template <class Z,class R>
void SkyLineMatrix<Z,R>::crout(double eps) const  {
    R  *ij , *ii  , *ik , *jk , xii, *dkk;
    Z i,j,k;
    double eps2=eps*eps;
    if (L != U) MATERROR(1,"Skyline matrix  non symmetric");
    if(verbosity>3 || verb >1 )
        cout << "  -- SkyLineMatrix Factorize/Crout   " << endl;
    U = 0; //
    typefac = FactorizationCrout;
    
    ij = L ; // pointeur sur le terme ij de la matrice avec j<i
    for (i=1;i<this->n;i++) // boucle sur les lignes
    { ii = L+pL[i+1]; // pointeur sur le terme fin de la ligne +1 =>  ij < ii;
        xii = D[i] ;
        for ( ; ij < ii ; ij++) // pour les j la ligne i
        { j = i -(ii - ij);
            k = max( j - (pL[j+1]-pL[j]) ,  i-(pL[i+1]-pL[i]) );
            ik =  ii - (i - k);
            jk =  L + pL[j+1] -(j - k);
            dkk = D + k;
            k = j - k ;
            R s=-*ij;
            while ( k-- ) s += *ik++ * *jk++ * *dkk++;
            *ij = -s/ *dkk ; // k = j ici
            
            xii -= *ij * *ij * *dkk;
        }
        if (std::norm(xii) <= max(eps2*std::norm(D[i]),1.0e-60))
        {
            cout << " Crout: zero pivot (" << i << " )= " << abs(xii)<< " <= " << eps*abs(D[i])
            << " eps = " << eps <<endl;
            MATERROR(3,"Crout SkyLine pivot ");
        }
            D[i] = xii;
    }
}
template <class Z,class R>
void SkyLineMatrix<Z,R>::LU(double eps) const  {
    R s,uii;
    double eps2=eps*eps;
    Z i,j;
    if (L == U && ( pL[this->n]  || pU[this->n] ) ) MATERROR(3,"matrix LU  symmetric");
    if(verbosity>3 || verb >1 )
        cout << "  -- SkyLineMatrix Factorize/LU SkyLine  " << endl;
    typefac=FactorizationLU;
    
    for (i=1;i<this->n;i++) // boucle sur les sous matrice de rang i
    {
        // for L(i,j)  j=j0,i-1
        Z j0 = i-(pL[i+1]-pL[i]);
        for ( j = j0; j<i;j++)
        {
            Z k0 = max(j0,j-(pU[j+1]-pU[j]));
            R *Lik = L + pL[i+1]-i+k0; // lower
            R *Ukj = U + pU[j+1]-j+k0; // upper
            s =0;
#ifdef WITHBLAS
            s = blas_sdot(j-k0,Lik,1,Ukj,1);
            Lik += j-k0;
#else
            for (int k=k0;k<j;k++) // k < j < i ;
                s += *Lik++ * *Ukj++ ;     // a(i,k)*a(k,j);
#endif
            *Lik -= s;
            *Lik /= D[j]; //  k == j here
        }
        // for U(j,i) j=0,i-1
        j0=i-pU[i+1]+pU[i];
        for (j=j0;j<i;j++)
        {
            s = 0;
            Z k0 = max(j0,j-pL[j+1]+pL[j]);
            R *Ljk = L + pL[j+1]-j+k0;
            R *Uki = U + pU[i+1]-i+k0;
#ifdef WITHBLAS
            s = blas_sdot(j-k0,Ljk,1,Uki,1);
            Uki += j-k0;
#else
            for (int k=k0  ;k<j;k++)    //
                s +=  *Ljk++ * *Uki++ ;
#endif
            *Uki -= s;  // k = j here
        }
        // for D (i,i) in last because we need L(i,k) and U(k,i) for k<j
        Z k0 = i-min(pL[i+1]-pL[i],pU[i+1]-pU[i]);
        R *Lik = L + pL[i+1]-i+k0; // lower
        R *Uki = U + pU[i+1]-i+k0; // upper
        s =0;
#ifdef WITHBLAS
        s = blas_sdot(i-k0,Lik,1,Uki,1);
#else
        for (Z k=k0;k<i;k++) // k < i < i ;
            s += *Lik++ * *Uki++ ;     // a(i,k)*a(k,i);
#endif
        // cout << " k0 " << k0 << " i = " << i << " " <<  s << endl;
        uii = D[i] -s;
        
        if (std::norm(uii) <= max(eps2*std::norm(D[i]),1.0e-30))
            MATERROR(3,"LU SkyLine pivot ");
        
        D[i] = uii;
        
    }
}

// solve in place ...
template <class Z,class R>
R* SkyLineMatrix<Z,R>::solve(R*x,int trans) const
{
    
    // --------------------------------------------------------------------
    //   si La diagonal D n'existe pas alors on suppose 1 dessus (cf crout)
    // --------------------------------------------------------------------
    R * v = x;
    
   
    const R *ij ,*ii, *ik, *ki;
    R *xk,*xi;
    int i;
    switch (this->typefac) {
        case FactorizationNO:
            if (this->U && this->L) {std::cerr << "APROGRAMMER (KN_<R><R>::operator/SkyLineMatrix)";
                MATERROR(10,"SkyLineMatrix solve no factorized ");}
            
            if ( this->U && !this->L )
            { // matrice triangulaire superieure
               if( verbosity> 5 || verb>2) cout << "  -- backward substitution  " << (this->D ? "DU" : "U") << endl;
                ki = this->U + this->pU[n];
                i = n;
                while ( i-- )
                { ii = this->U + this->pU[i];
                    xi= xk  = v +  i ;
                    if (this->D) *xi /= this->D[i];// pour crout ou LU
                    while ( ki > ii)
                        *--xk  -=  *--ki *  *xi ;
                }
            }
            else if  ( !this->U && this->L )
            { // matrice triangulaire inferieure
                if( verbosity> 5 || verb>2) cout << "  -- forward substitution "  <<( this->D ? "LD" : "L" ) <<endl;
                ii = this->L;
                for (i=0; i<n; i++)
                { ij = ik = (this->L + this->pL[i+1]) ;  // ii =debut,ij=fin+1 de la ligne
                    xk = v + i;
                    R ss = v[i];
                    while ( ik > ii)
                        ss -= *--ik * *--xk ;
                    if ( this->D) ss /= this->D[i];// pour crout ou LU
                    v[i] = ss ;
                    ii = ij;
                }
            }
            else if (this->D)
            { // matrice diagonale
                if( verbosity> 5 || verb>2) cout << "  -- diagonal substitution D" <<endl;
                for (i=0;i<n;i++)
                    v[i]=v[i]/this->D[i];
            }
            break;
        case FactorizationCholesky:
            //if(verbosity>9|| verb>5 ) cout << " Factorization Choslesky" << endl;
           ld().solve(x);
           ldt().solve(x);
            break;
        case FactorizationCrout:
          //     if(verbosity>4 || verb>1) cout << " Factorization Crout" << endl;
           this->l().solve(x);
           this->d().solve(x);
           this->lt().solve(x);
           break;
        case FactorizationLU:
         //  if(verbosity>4|| verb>1) cout << "  LU" << endl;
          if(trans)
          {
              this->dut().solve(x);
              this->lt().solve(x);
          }
          else
          {
          this->l().solve(x);
          this->du().solve(x);
          }
            break;
        default:
            std::cerr << "  (operator /=(SkyLineMatrix, Error unknown type of Factorization  =" << typefac<< std::endl;
            MATERROR(1,"unknown type of Factorization");
    }
    return x;
}

template <class Z,class R>
Z SkyLineMatrix<Z,R>::size() const {
    Z s = sizeof(SkyLineMatrix<Z,R>);
    if (D) s += this->n*sizeof(R);
    if (pL) s += this->n*sizeof(int);
    if (pU && (pU != pL)) s += this->n*sizeof(Z);
    if (L) s += pL[this->n]*sizeof(int);
    if (U && (U != L)) s += pU[this->n]*sizeof(Z);
    return s;
}

template <class Z,class R>
std::ostream&  SkyLineMatrix<Z,R>::dump (std::ostream& f) const
{f<< " matrix skyline " << this->n << '\t' << this->m << '\t' ;
    f <<  "  this " << endl;
    f << " pL = " << pL << " L ="  << L << endl
    << " pU = " << pU << " U ="  << U << endl
    << " D = " << D << endl;
    if ( (pL == pU) &&  (U == L) )
        if (pL && L)
        {f << " skyline symmetric " <<endl;
            Z i,j,k;
            for (i = 0;i<this->n;i++)
            { f << i << " {" << pL[i+1]-pL[i] << "}" <<'\t' ;
                for (k=pL[i];k<pL[i+1];k++)
                { j=i-(pL[i+1]-k);
                    f << j << " " << L[k] << "; ";
                }
                f <<  i  << ":" << D[i] << endl  ;
            }
        }
        else f << "Skyline: pointeur null " <<endl;
        else
        {
            f << " Skyline  non symmetric " << endl;
            Z i,k;
            for (i = 0;i<this->n;i++)
            {
                f << i ;
                if (pL && L)
                {
                    f << " jO=" << i-pL[i+1]+pL[i] << " L= " <<'\t' ;
                    for (k=pL[i];k<pL[i+1];k++)
                    {
                        f <<  " " << L[k] ;
                    }
                }
                if (D)
                    f  << " D= " << D[i]  << '\t' ;
                else
                    f  << " D=0 => 1 ; ";
                if (pU && U)
                {
                    f << " i0=" << i-pU[i+1]+pU[i] << " U= " <<'\t' ;
                    for (k=pU[i];k<pU[i+1];k++)
                        f << " " << U[k] ;
                    
                }
                f << endl;
            }
            
        }
    return f;
}
template <class Z,class R>
inline std::ostream& operator <<(std::ostream& f,const SkyLineMatrix<Z,R> & m) {return m.dump(f);}


