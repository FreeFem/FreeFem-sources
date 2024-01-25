/*
 GradientConjugue / version C++
 */
#include <cassert>
#include <cmath>
#include <iostream>
#include "CG.hpp"
#include <complex>

//int verbosity =2  ;
typedef std::complex<double> Complex;


template<class I,class K>
K * myscopy(I n,const K *x,K *y)
{
    
    for(I i=0;i<n;++i)
        y[i]=x[i];
    return y;
}
template<class I,class K>
K mysdot(I n,const K *x,const K *y)
{
    K s=0;
    for(I i=0;i<n;++i)
        s+= x[i]*y[i];
    return  s;
}
template<class I>
Complex mysdot(I n,const Complex *x,const Complex *y)
{
    Complex s=0;
    for(I i=0;i<n;++i)
        s+= std::conj(x[i])*y[i];
    return  s;
}
template<class I,class K>
K * mysaxpy(I n,K a,const K *x,K *y)
{
    for(I i=0; i< n; ++i)
        y[i] += a* x[i];
    return  y;
}
template<class I,class K>
K * mysax2y(I n,K a,const K *x,K *y)
{
    for(I i=0; i< n; ++i)
        y[i] = a* x[i];
    return  y;
}

template<class I,class K>
K * myscal(I n,K a,K *x)
{
    for(I i=0; i<n;++i)
        x[i] *= a;
    return  x;
}

template<class I,class K>
double mysnrm2(I n,const K *x)
{
    double s=0.;
    for(I i=0; i<n;++i)
        s += std::norm(x[i]);
    return  std::sqrt(s);
}

#ifdef WITHBLAS

#ifdef __cplusplus
extern "C" {
#endif
    // BLAS function
    double ddot_(int *n,const double *x, int *incx,const double *y, int *incy);
    void  dcopy_(int *n,const double *x, int *incx,const double *y,int *incy);
    void  dscal_(int *n,const double *alpha, double *x, int *incx);
    void  daxpy_(int *n,const double *alpha,const double *x, int *incx,
                 double *y, int *incy);
    double dnrm2_(int *n,const double *x, int *incx);
    void   idmin_(int *n,const double *sx, int *incx);
    complex zdot_(int *n,const Complex *x, int *incx,const Complex *y, int *incy);
    void  zcopy_(int *n,const Complex *x, int *incx, Complex *y,int *incy);
    void  zscal_(int *n,const Complex *alpha, Complex *x, int *incx);
    void  zaxpy_(int *n,const Complex *alpha,const Complex *x, int *incx,
                 Complex *y, int *incy);
    double znrm2_(int *n,const Complex *x, int *incx);
    
#ifdef __cplusplus
}
#endif

template<int,double>
double * mysnrm2(int n,double *x)
{
    int un=1;
    return dnrm2_(&n,x,un);
}


template<int,double>
double * myscopy(int n,const double *x,double *y)
{
    int un=1;
    dcopy_(&n,x,&un,y,&un);
    return y;
}
template<int,double>
double mysdot(int n,const double *x,const double *y)
{
    int un=1;
    return  ddot_(&n,x,&un,y,&un);
}
template<int,double>
double * mysaxpy(int n,double a,const double *x,double *y)
{
    int un=1;
    daxpy_(&n,&a,x,&un,y,&un);
    return  y;
}
template<int,double>
double * mysax2y(int n,double a,const double *x,double *y)
{
    int un=1;
    std::fill(y,y+1,0.);
    daxpy_(&n,&a,x,&un,y,&un);
    return  y;
}

template<int,double>
double * myscal(int n,double a,double *x)
{
    int un=1;
    dscal_(&n,&a,x,&un);
    return  x;
}

template<int,Complex>
double * mysnrm2(int n,const Complex *x)
{
    int un=1;
    return znrm2_(&n,x,un);
}

template<int,Complex>
Complex * myscopy(int n,const Complex *x,Complex *y)
{
    int un=1;
    zcopy_(&n,x,&un,y,&un);
    return y;
}
template<int,Complex>
Complex mysdot(int n,const Complex *x,const Complex *y)
{
    int un=1;
    return  zdotc_(&n,x,&un,y,&un);
}
template<int,Complex>
Complex * mysaxpy(int n,Complex a,const Complex *x,Complex *y)
{
    int un=1;
    zaxpy_(&n,&a,x,&un,y,&un);
    return  y;
}
template<int,Complex>
Complex * mysax2y(int n,Complex a,const Complex *x,Complex *y)
{
    int un=1;
    std::fill(x,x+n,Complex(0.));
    zaxpy_(&n,&a,x,&un,y,&un);
    return  y;
}

template<int,Complex>
Complex * myscal(int n,Complex a,Complex *x)
{
    int un=1;
    zscal_(&n,&a,x,&un);
    return  x;
}

#else

#endif


template<class I,class K>
K * SetArrayGC(I *x,long n,K c)
{
    std::fill(x,x+n,c);
    return x;
}

template<class TypeIndex,class TypeScalar>
int ConjugueGradient(CGMatVirt<TypeIndex,TypeScalar> &A, // fonction et pointeur data pour A
                     CGMatVirt<TypeIndex,TypeScalar>  &C, // fonction et pointeur data pour C
                     TypeScalar * b, // second membre
                     TypeScalar * x, // solution qui contient une initialisation
                     int &nbitermax,// change FH mars 2020 add &
                     double &eps,// change FH mars 2020 add &
                     int niveauimpression)

{
    using namespace std;
    typedef TypeIndex Z;
    typedef TypeScalar K;
    
    // niveauimpression: 0: pas impression .. 10 a chaque iteration
    Z n = A.n, nret=0;
    K NaN = nan("");
    K * G = SetArrayGC(new K[n] ,n,NaN);
    K * CG= SetArrayGC(new K[n],n,NaN);
    K *AH=CG;// meme tableau que CG Attention ??????
    K *H=SetArrayGC(new K[n],n,NaN);
    K rho, gamma;
    K minus1=-1.;
    double gCgp, gCg, eps2=eps*eps;
    assert( A.m==n && C.m==n && C.n==n);
    niveauimpression=std::min(niveauimpression,10);
    int nprint =std::max((int)(max(nbitermax+1,1000)*(10.-niveauimpression)/10.),1);
    
    mysaxpy(n,minus1,b,A.matmul(x,G));// G = Ax -b
    gCg = real(mysdot(n,G,C.matmul(G,CG))) ;
    myscal(n,minus1,myscopy(n,CG,H)); // H =- CG;
    if( eps >0 ) { eps2 *=gCg ; eps = sqrt(eps2); }
    
    assert( !(gCg != gCg) ); //  verif if NaN => bad Matrix..
    if(gCg < 1e-30)
    {   if(niveauimpression)
        std::cout << " GC: converge after 0 iteration ||g||_C^2" << gCg  << std::endl;
        nret = 2;
        nbitermax=0;
    }
    else
        for(int iter=1; iter <= nbitermax; ++iter)
        {
            gCgp = gCg;
            rho = real(-mysdot(n,G,H)/mysdot(n,H,A.matmul(H,AH)));
            
            mysaxpy(n,rho,H,x);
            mysaxpy(n,rho,AH,G);
            gCg = real(mysdot(n,G,C.matmul(G,CG)));
            gamma = gCg/ gCgp;
            mysaxpy(n,minus1,CG,myscal(n,gamma,H));
            
            if(gCg < eps2) // We have converged ...
            {
                if(niveauimpression)
                    std::cout << " GC:  converge after " <<iter
                    << " g=" << gCg << " rho= " << rho << " gamma= " <<gamma<<std::endl;;
                nret= 1;
                nbitermax= iter;
                break;
            }
            else
                if ( ((iter+1) % nprint) == 0 )
                    std::cout <<"  GC:iteration "<< iter << " rho "<< rho << " gamma "
                    <<gamma<< " ||g||_C^2:" << gCg << " / " << eps2 <<std::endl;
        }
    delete[] H;
    delete[] CG;
    delete[] G;
    return nret;
}

template
int ConjugueGradient<int,double> (CGMatVirt<int,double> &A, // fonction et pointeur data pour A
                                  CGMatVirt<int,double>  &C, // fonction et pointeur data pour C
                                  double * b, // second membre
                                  double * x, // solution qui contient une initialisation
                                  int &nbitermax,
                                  double &eps,
                                  int niveauimpression)
;

template
int ConjugueGradient<int,Complex> (CGMatVirt<int,Complex> &A, // fonction et pointeur data pour A
                                   CGMatVirt<int,Complex>  &C, // fonction et pointeur data pour C
                                   Complex * b, // second membre
                                   Complex * x, // solution qui contient une initialisation
                                   int &nbitermax,
                                   double &eps,
                                   int niveauimpression)
;

template
int ConjugueGradient<long,double> (CGMatVirt<long,double> &A, // fonction et pointeur data pour A
                                  CGMatVirt<long,double>  &C, // fonction et pointeur data pour C
                                  double * b, // second membre
                                  double * x, // solution qui contient une initialisation
                                  int &nbitermax,
                                  double &eps,
                                  int niveauimpression)
;

template
int ConjugueGradient<long,Complex> (CGMatVirt<long,Complex> &A, // fonction et pointeur data pour A
                                   CGMatVirt<long,Complex>  &C, // fonction et pointeur data pour C
                                   Complex * b, // second membre
                                   Complex * x, // solution qui contient une initialisation
                                   int &nbitermax,
                                   double &eps,
                                   int niveauimpression)
;


typedef CGMatVirt<int,double> MatVirt;
#define DBL_EPSILON 2.2204460492503131e-16 // double epsilon

#define TINY 1.0e-20
//typedef double FLOAT;
//typedef double REAL;
//#defined WITH_MPI
template<typename K>
void UnDoPermute(int n,int *p,K *x,K *tmp)
{
    if(p)
    {
        for(int i=0; i< n; ++i)
            tmp[p[i]]=x[i];
        for(int i=0; i< n; ++i)
            x[i]= tmp[i];
    }
}
template<typename K>

void DoPermute(int n,int *p,K *x,K *tmp)
{
    if(p)
    {
        for(int i=0; i< n; ++i)
            tmp[i]=x[p[i]];
        for(int i=0; i< n; ++i)
            x[i]= tmp[i];
    }
}




#include "RNM.hpp"

    double ffconj(double x) {return x;}
    Complex ffconj(Complex x) {return std::conj(x);}

template<typename K,typename Z=int>
bool fgmres(CGMatVirt<Z,K> &A, // fonction et pointeur data pour A
            CGMatVirt<Z,K> &CC,int leftC,
            K *prhs,
            K *px,
            double &eps,
            int &nbitermx,
            int nbkrylov,
            int verbo,
            int *wbc)
{
    //  This code is base of  fgmres of F. Nataf, P. Jolivet and P-H Tournier.
    // witten in freefem++
    
    verbo=std::min(verbo,10);
    int nprint =std::max((int)(max(nbkrylov+1,100)*(10.-verbo)/10.),1);
    leftC=1;
    typedef double R;
    R relerr=1e100 , relres=1e100,normb=0.;
    Z n = A.n;
    K NaN = nan("");
    K * rot0 = SetArrayGC(new K[nbkrylov+2] ,nbkrylov+2,NaN);
    K * rot1 = SetArrayGC(new K[nbkrylov+2] ,nbkrylov+2,NaN);
    K * g = SetArrayGC(new K[nbkrylov+1] ,nbkrylov+1,NaN);
    K * g1 = SetArrayGC(new K[nbkrylov+1] ,nbkrylov+1,NaN);
    leftC = 0;
    CGMatVirtId<Z,K> MatId(n);
    CGMatVirt<Z,K> & C=CC;
    CGMatVirt<Z,K> *pCl = leftC ? &C  :  &MatId;
    CGMatVirt<Z,K> *pCr = !leftC ? &C :  &MatId;
    Z nrestart=0;
    
    
    KN<K> rwi(n),uni(n),x0(n),ri(n);
    KNM<K> Hn(nbkrylov+2,nbkrylov+1);
    KN_<K> x(px,n),rhs(prhs,n);
    KN< KN<K> > Vi(nbkrylov+1),Vpi(nbkrylov+1);
    KN<K> zi(n),vi(n),wi(n);
    Hn = K();
    wi = rhs; // remove tgv ????
//    ffassert(wbc) ;
    int kk=0;
    if( wbc)
     for(int i=0; i< n; ++i)
      if( wbc[i]) // remove TGV in RHS ...
        kk++,wi[i]=K();
   // cout << " nbcl " << kk<< endl;
    pCl->matmul(wi,vi);
    // remove BC part
    
    normb = mysnrm2(n,(K*) vi);
    //cout << " normb " << normb << " "<< endl; 
    uni = x;
    K minus1 = K(-1.);
    bool noconv = true;
    int iter =0;//
    while (noconv)
    {
        x0=uni;
        // ri[] = matA(uni);  ri[] -= rhs;  ri[] *= -1.0;
        myscal(n,minus1,mysaxpy(n,minus1,prhs,A.matmul(uni,ri)));
        pCl->matmul(ri,zi);
        // g[0] = sqrt(real(pr#scalprod(zi[],zi[])));
        g[0]=mysnrm2(n,(K*)zi);
        if(verbo>3)
            cout << "  ** fgmres: " << iter << " residus 0 " << abs(g[0]) << " Cl: " << normb << endl ;
        if (normb < 1.e-20 || eps < 0) normb = 1.;
        Vi[0]=(1./g[0])*zi;
        int it; // need for reconstruction
        for( it=0; it<nbkrylov; it++,iter ++)
        {
            pCr->matmul(Vi[it],vi);
            
            if (!leftC) {
                C.matmul(Vi[it],vi);// preCON(Vi[it][]);
                Vpi[it]=vi;
                //  vi[]=Vpi[it][];
                A.matmul(vi,wi);// wi[]=matA(vi[]);
                
            }
            else {
                A.matmul(Vi[it],vi);// vi[]=matA(Vi[it][]);
                C.matmul(vi,wi);// wi[] = preCON(vi[]);
            }
            // mofif Gram Schmidt
            for(int i=0; i<it+1; i++)
            {
                Hn(i,it) = mysdot(n, (K*)wi,(K*)Vi[i]);
                mysaxpy(n,-ffconj(Hn(i,it)),(K*)Vi[i],(K*)wi);//wi = ffconj(Hn(i,it))* Vi[i];
            }
            K aux = Hn(it+1,it) = mysnrm2(n,(K*)(K*)wi);
            
            Vi[it+1] = (1./aux)*wi;
            /* QR decomposition of Hn*/
            for(int i=0; i<it; i++)
            {      /* QR decomposition of Hn*/
                K aa = ffconj(rot0[i])*Hn(i,it )+ffconj(rot1[i])*Hn(i+1,it);
                K bb = -rot1[i]*Hn(i,it)+rot0[i]*Hn(i+1,it);
                Hn(i,it) = aa;
                Hn(i+1,it) = bb;
               // cout << i << " " << aa << " " << bb << endl;
            }
            K sq = sqrt( ffconj(Hn(it,it))*Hn(it,it) + Hn(it+1,it)*Hn(it+1,it) );
            rot0[it] = Hn(it,it)/sq;
            rot1[it] = Hn(it+1,it)/sq;
           // cout << " sq " << sq << " " << rot0[it] << " " << rot1[it] <<  endl;
            
            Hn(it,it) = ffconj(rot0[it])*Hn(it,it)+ffconj(rot1[it])*Hn(it+1,it);
            Hn(it+1,it) =  0.;
            g[it+1] = -rot1[it]*g[it];
            g[it] = ffconj(rot0[it])*g[it];
            
            
        
            relres = abs(g[it+1]);//  residu gmres ||Ax -b ||_2
            if ((iter+1) % nprint ==0)
                cout << "     fgmres "<< iter << " Res:  = " << relres << " Rel res = " << relres/normb <<   endl;
            
            if(relres/normb < abs(eps)) {
                noconv= false;
                if (verbo ) {
                    cout << "  **  fgmres has converged in " << (iter) << " iterations "
                    << "The relative residual is " <<  relres/normb << " Cl: " << normb << endl;
                }
                break;
            }
           if( it > nbitermx) break; // no converge
        }
        it = min(it,nbkrylov-1);
        /* Reconstruct the solution */
        // use g0, g1 , Hn , Vp, Vpi (fgmres)
        
        KN<K> y(it+1);
        for(int i=it; i>=0; i--)
        {
            g1[i] = g[i];
            for(int j=i+1; j<it+1; j++)
            g1[i] = g1[i]-Hn(i,j)*y[j];
            y[i]=g1[i]/Hn(i,i);
        }
        
        wi = K();
        for(int i=0;i<it+1;i++){
            if (!leftC)
            wi +=  ffconj(y[i])*Vpi[i];// ICI Fgmres ... pas de tableau Vpi ...
            else
            wi +=  ffconj(y[i])*Vi[i];
        }
        uni = wi; //
       // pCr->matmul(wi,uni);// uni[]= Precon(uni[]) if (!leffC) pas flexi 
        uni += x0;//
        x = uni;
        // Fin reconstruction de la solution

        if(!noconv) break; 
        if( iter > nbitermx) break; // no converge
        if( (nrestart++< verbo) && (verbo> 2) )
            cout << "  ** restart fgmres iter " <<iter << " res:  " << relres/normb << endl;
    }
    if(noconv && verbo  )
    {
        cout << " !!!!!!!! fgmres has not  converged in " << iter << " iterations "
        << "The relative residual is " <<  relres/normb << " Cl: " << normb << endl;
    }
    nbitermx=iter; //  to
    delete [] g1;
    delete [] g;
    delete [] rot1;
    delete [] rot0;
    return !noconv;
}

template
bool fgmres(CGMatVirt<int,double> &A, // fonction et pointeur data pour A
            CGMatVirt<int,double> &C,int leftC,
            double *y,
            double *x,
            double &tol,
            int &maxits,
            int restart,
            int verb,
            int *wbc);
template
bool fgmres(CGMatVirt<int,std::complex<double> > &A, // fonction et pointeur data pour A
            CGMatVirt<int,std::complex<double> > &C,int leftC,
            std::complex<double>  *y,
            std::complex<double>  *x,
            double &tol,
            int &maxits,
            int restart,
            int verb,
            int *wbc);
template
bool fgmres(CGMatVirt<long,double> &A, // fonction et pointeur data pour A
            CGMatVirt<long,double> &C,int leftC,
            double *y,
            double *x,
            double &tol,
            int &maxits,
            int restart,
            int verb,
            int *wbc );
template
bool fgmres(CGMatVirt<long,std::complex<double> > &A, // fonction et pointeur data pour A
            CGMatVirt<long,std::complex<double> > &C,int leftC,
            std::complex<double>  *y,
            std::complex<double>  *x,
            double &tol,
            int &maxits,
            int restart,
            int verb,
            int *wbc );


template double * myscopy<int,double>(int n,const double *x,double *y);
template double * myscal<int,double>(int n,double  a,double *x);
template double * myscopy<unsigned long,double>(unsigned long n,const double *x,double *y);
template double * myscal<unsigned long,double>(unsigned long n,double  a,double *x);
template double * mysaxpy(int n,double a,const double *x,double *y);

template Complex * myscopy<int,Complex>(int n,const Complex *x,Complex *y);
template Complex * myscal<int,Complex>(int n,Complex  a,Complex *x);
template Complex * myscopy<unsigned long,Complex>(unsigned long n,const Complex *x,Complex *y);
template Complex * myscal<unsigned long,Complex>(unsigned long n,Complex  a,Complex *x);
template Complex * mysaxpy(int n,Complex a,const Complex *x,Complex *y);
