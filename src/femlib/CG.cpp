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
Complex * myscal(int n,Complex a,Complex *x)
{
    int un=1;
    zscal_(&n,&a,x,&un);
    return  x;
}

#else

#endif


template<class I,class K>
K * SetArray(I *x,long n,K c)
{
    std::fill(x,x+n,c);
    return x;
}

template<class TypeIndex,class TypeScalar>
int ConjugueGradient(CGMatVirt<TypeIndex,TypeScalar> &A, // fonction et pointeur data pour A
                     CGMatVirt<TypeIndex,TypeScalar>  &C, // fonction et pointeur data pour C
                     TypeScalar * b, // second membre
                     TypeScalar * x, // solution qui contient une initialisation
                     int nbitermax,
                     double eps,
                     int niveauimpression)

{
    using namespace std;
    typedef TypeIndex Z;
    typedef TypeScalar K;
    
    // niveauimpression: 0: pas impression .. 10 a chaque iteration
    Z n = A.n, nret=0;
    K NaN = nan("");
    K * G = SetArray(new K[n] ,n,NaN);
    K * CG= SetArray(new K[n],n,NaN);
    K *AH=CG;// meme tableau que CG Attention ??????
    K *H=SetArray(new K[n],n,NaN);
    K rho, gamma;
    K minus1=-1.;
    double gCgp, gCg, eps2=eps*eps;
    assert( A.m==n && C.m==n && C.n==n);
    niveauimpression=std::min(niveauimpression,10);
    int nprint =std::max((int)(max(nbitermax+1,1000)*(10.-niveauimpression)/10.),1);
    
    mysaxpy(n,minus1,b,A.matmul(x,G));// G = Ax -b
    gCg = real(mysdot(n,G,C.matmul(G,CG))) ;
    myscal(n,minus1,myscopy(n,CG,H)); // H =- CG;
    if( eps >0 ) eps2 *=gCg ;
    
    assert( !(gCg != gCg) ); //  verif if NaN => bad Matrix..
    if(gCg < 1e-30)
    {   if(niveauimpression)
        std::cout << " GC: on a converge on 0 iteration ||g||_C^2" << gCg  << std::endl;
        nret = 2;}
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
                    std::cout << " GC: on a converge en iteration " <<iter
                    << " g=" << gCg << " rho=" << rho << " gamma=" <<gamma<<std::endl;;
                nret= 1;
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
                                  int nbitermax,
                                  double eps,
                                  int niveauimpression)
;

template
int ConjugueGradient<int,Complex> (CGMatVirt<int,Complex> &A, // fonction et pointeur data pour A
                                   CGMatVirt<int,Complex>  &C, // fonction et pointeur data pour C
                                   Complex * b, // second membre
                                   Complex * x, // solution qui contient une initialisation
                                   int nbitermax,
                                   double eps,
                                   int niveauimpression)
;

template
int ConjugueGradient<long,double> (CGMatVirt<long,double> &A, // fonction et pointeur data pour A
                                  CGMatVirt<long,double>  &C, // fonction et pointeur data pour C
                                  double * b, // second membre
                                  double * x, // solution qui contient une initialisation
                                  int nbitermax,
                                  double eps,
                                  int niveauimpression)
;

template
int ConjugueGradient<long,Complex> (CGMatVirt<long,Complex> &A, // fonction et pointeur data pour A
                                   CGMatVirt<long,Complex>  &C, // fonction et pointeur data pour C
                                   Complex * b, // second membre
                                   Complex * x, // solution qui contient une initialisation
                                   int nbitermax,
                                   double eps,
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


/*------------- classical Gram - Schmidt -------------------*/
// real ..
template<typename K>
void GramSchmidt(int nloc, int i1,int ptih, K *vv,K *hh)
{
    int  pti1 = i1*nloc;
    K alpha;
    for(int j=0; j<i1; j++)
    {
        hh[ptih+j] = mysdot(nloc, &vv[j*nloc], &vv[pti1]);
        alpha = -hh[ptih+j];
        mysaxpy(nloc, alpha, &vv[j*nloc], &vv[pti1]);
    }
}

/* Givens Rotation */

double sgn(double x, double y){
    /* Sign transfer function */
    if (y >= 0.0) return fabs(x);
    else return -fabs(x);
}
double abssq(Complex  x){
    return (pow(real(x), 2) + pow(imag(x),2));
    //  return cpow(cabs(x),2);
}

void zclartg(Complex  f, Complex g, double *cs, Complex  *sn, Complex *rot)
{
    const double  epsmac = 1.0e-16;
    double one=1.0, zero=0.0;
    Complex czero(0.0,0.0);
    double D, F2, G2, FG2;
    if (abs(g) <= epsmac)   {
        *cs = sgn(one, real(f));
        *sn = czero;
        *rot = f*(*cs);
    }
    else if (std::abs(f) <= epsmac){
        *cs = zero;
        *sn = std::conj(g)/std::abs(g);
        *rot = std::abs( g );
    }
    else{
        F2 = abssq(f);
        G2 = abssq( g );
        FG2 = F2 + G2;
        if(fabs(FG2) <= epsmac) FG2 = epsmac;
        D = 1/sqrt(F2*FG2);
        *cs = F2 * D;
        FG2 = FG2*D;
        *rot = f*FG2;
        *sn = f*D;
        *sn = conj(g)*(*sn);
    }
}

double PlaneRot(int nloc,int i,int ptih,Complex *hh,double *c,Complex *s,Complex *rs)
{
    Complex t1;
    int  i1=i+1;
    double ro=0;
    Complex rot;
    if (i != 0) {
        for(int k = 1; k <= i; k++) {
            int k1 = k-1;
            t1 = hh[ptih+k1];
            
            hh[ptih+k1] = c[k1]*t1 + s[k1]*hh[ptih+k];
            hh[ptih+k] = -conj(s[k1])*t1 + c[k1]*hh[ptih+k];
        }
    }
    /*-----------get next plane rotation------------ */
    zclartg(hh[ptih+i], hh[ptih+i1], &c[i], &s[i], &rot);
    rs[i1] = -conj(s[i])*rs[i];
    rs[i] = c[i]*rs[i];
    hh[ptih+i] = rot;
    ro = abs(rs[i1]);
    return ro;
}

double PlaneRot(int nloc,int i,int ptih,double *hh,double *c,double *s,double *rs)
{
    
    double gam,ro;
    double t1;
    int i1 = i + 1;
    if (i != 0) {
        for(int k=1; k<=i; k++){
            int k1 = k-1;
            t1 = hh[ptih+k1];
            hh[ptih+k1] = c[k1]*t1 + s[k1]*hh[ptih+k];
            hh[ptih+k] = -s[k1]*t1 + c[k1]*hh[ptih+k];
        }
    }
    /*-----------get next plane rotation------------ */
    gam = sqrt(hh[ptih+i]*hh[ptih+i] + hh[ptih+i1]*hh[ptih+i1]);
    /*
     if gamma is zero then any small value will do ...
     will affect only residual estimate
     */
    if (fabs(gam) <= TINY) gam = TINY;
    
    /* determine-next-plane-rotation */
    c[i] = hh[ptih+i]/gam;
    s[i] = hh[ptih+i1]/gam;
    
    rs[i1] = -s[i]*rs[i];
    rs[i] = c[i]*rs[i];
    /* determine res. norm and test for convergence */
    hh[ptih+i] = c[i]*hh[ptih+i] + s[i]*hh[ptih+i1];
    ro = fabs(rs[i1]);
    return ro;
}

template<typename K,typename Z=int>
bool fgmres(CGMatVirt<Z,K> &A, // fonction et pointeur data pour A
            CGMatVirt<Z,K> &C,
            K *y,
            K *x,
            double tol,
            int maxits,
            int restart,
            int verbo,
            int *perm)
{
   // typedef K FLOAT ;
    typedef double R ;
    
    bool outflag, intflag;
    Z i, i1, pti, pti1, ptih, j, its;
    Z ii, jj, k, k1, size;
    K oneK = 1.;
    K  *vv=0, *z=0, *hh=0, *s=0, *rs=0, t1=0.;
    R  eps1,  ro=1e100, t, *c;
    int n = A.n;
    int nloc= n;
    verbo=std::max(verbo,10);
    int nprint =std::max((int)((maxits+1)*(10.-verbo)/10.),1);
    
 #if defined(WITH_MPI)
    int rank;
    MPI_Comm comm;
    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
    eps1 = tol;
    
    /*----- allocate memory for working local arrays -------*/
    
    /* ---------allocate size for nx(m+1) matrix vv */
    size = nloc*(restart+1);
    vv=new K[size];
    /*-------- allocate size for nxm matrix z -----*/
    size = nloc*restart;
    z=new K[size];
    
    /*----- allocate memory for Hessenberg matrix (nonzeros only)
     *----- and rotation vectors s and rs -------------------*/
    size = (((restart+1)*(restart+2)/2) - 1) + 2*(restart+1);
    hh=new K[size];
    
    s = hh + ((restart+1)*(restart+2)/2) - 1;
    rs = s + restart + 1;
    /*--------- allocate memory for rotation matrix c -------
     * This is done separately since for complex valued systems
     * c is still a real-valued vector, and hence cannot be
     * allocated as sizeof(complex double) as is the case for
     * s and rs and hh above --------------------------------*/
    c=new R[ restart+1];
    
    DoPermute(n,perm,x,vv);
    DoPermute(n,perm,y,vv);
    
    /* outer loop starts here */
    its = 0;
    
    outflag = true;
    while(outflag) {
        /* compute vv[0] = A * x */
        A.matmul(x,vv);
        mysaxpy(n,-oneK,y,vv);// vv = A x - y
        /* compute the norm of the residual */
        ro=mysnrm2(n,vv);
      //  std::cout << " ro 0 = " << ro << std::endl;
        if(fabs(ro) <= DBL_EPSILON)  {
            outflag = false;
            break;
        }
        t1 = 1.0 / ro;
        myscal(n,t1,vv);
        if(its == 0)
            eps1 = tol*ro;
        
        /* ----------initialize 1-st term of rhs of hessenberg system ------------*/
        
        rs[0] = ro;
        
        i = -1;
        pti = 0;
        pti1 = 0;
        ptih = 0;  //  
        intflag = true;
        while (intflag) {
            i++;
            its++;
            // std::cout << i << " " << ptih << std::endl;
             i1 = i + 1;
            pti = i*nloc;
            pti1 = i1*nloc;
            
            /*------------- preconditioning operation z = K^{-1}vv ---------------*/
            C.matmul(&vv[pti], &z[pti]);
            /*------------- compute A*z -----------------*/
            A.matmul(&z[pti], &vv[pti1]);
            
            /*------------- classical Gram - Schmidt -------------------*/
            GramSchmidt( nloc,  i1, ptih, vv, hh);
 
            t=mysnrm2(n,&vv[pti1]);
            
            hh[ptih+i1] = t;
            
            if (fabs(t) > TINY) {
                t1 = 1.0 / t;
                myscal(n,t1,&vv[pti1]);
            }
            
            /* done with classical Gram-Schmidt and Arnoldi step. now update
             * factorization of hh */
            ro=  PlaneRot( nloc, i, ptih,  hh, c, s,rs);
            
            /*------------ Check for convergence ---------*/
            if ((i+1 >= restart) || (ro <= eps1) || its >= maxits)
                intflag = false;
            else
            /*------------ update hh pointer ptih ---------*/
                ptih += i+2;
            if ( (its % nprint) == 0 )
            std::cout <<"  fgmres:iteration "<< its <<"/"<<maxits << " ro "<< ro << " / " << eps1 << " flags: " << intflag<< "/"<< outflag <<std::endl;

        }
        
        /* now compute solution first solve upper triangular system */
        rs[i] = rs[i]/hh[ptih+i];
        for (ii = 1; ii <= i; ii++) {
            k = i-ii;
            k1 = k+1;
            t1 = rs[k];
            for (j = k1; j <= i; j++) {
                jj = ((j+1)*(j+2)/2) - 1;
                t1 = t1 - hh[jj+k]*rs[j];
            }
            jj = ((k+1)*(k+2)/2)-1;
            rs[k] = t1/hh[jj+k];
        }
        /* done with back substitution. now form linear combination to
         * get solution */
        for (j = 0; j <= i; j++) {
            t1 = rs[j];
            //parms_VecAXPY(x, &z[j*nloc], t1, is);
            mysaxpy(nloc,t1,&z[j*nloc],x);
        }
        /* test for return */
        if ((ro <= eps1) || (its >= maxits)) {
            outflag = false;
            if(verbo)
                std::cout << "fgmres  stop  " << ro << " " << eps1 << " its "<< its << " " << maxits << std::endl;
        }
        
        
    }
    
    
    
    /* reset isvecperm and do inverse permutation*/
    /* permutes x and y */
    UnDoPermute(n,perm,x,vv);
    UnDoPermute(n,perm,y,vv);
    
    delete [] vv;
    delete [] z;
    delete [] hh;
    delete [] c;
    
    return (ro <= eps1);
}

template
bool fgmres(CGMatVirt<int,double> &A, // fonction et pointeur data pour A
            CGMatVirt<int,double> &C,
            double *y,
            double *x,
            double tol,
            int maxits,
            int restart,
            int verb,
            int *perm );
template
bool fgmres(CGMatVirt<int,std::complex<double> > &A, // fonction et pointeur data pour A
            CGMatVirt<int,std::complex<double> > &C,
            std::complex<double>  *y,
            std::complex<double>  *x,
            double tol,
            int maxits,
            int restart,
            int verb,
            int *perm );
template
bool fgmres(CGMatVirt<long,double> &A, // fonction et pointeur data pour A
            CGMatVirt<long,double> &C,
            double *y,
            double *x,
            double tol,
            int maxits,
            int restart,
            int verb,
            int *perm );
template
bool fgmres(CGMatVirt<long,std::complex<double> > &A, // fonction et pointeur data pour A
            CGMatVirt<long,std::complex<double> > &C,
            std::complex<double>  *y,
            std::complex<double>  *x,
            double tol,
            int maxits,
            int restart,
            int verb,
            int *perm );


template double * myscopy<int,double>(int n,const double *x,double *y);
template double * myscal<int,double>(int n,double  a,double *x);
template double * myscopy<unsigned long,double>(unsigned long n,const double *x,double *y);
template double * myscal<unsigned long,double>(unsigned long n,double  a,double *x);

