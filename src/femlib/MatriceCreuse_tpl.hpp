#ifndef  MatriceCreuse_tpl_
#define MatriceCreuse_tpl_
#ifndef  MatriceCreuse_h_
#include "MatriceCreuse.hpp"
#include <limits>
#include <set>
#include <list>
#include <map>
#endif

#ifndef __MWERKS__
// test blas
//  on MacOS9 under MWERKS
//  cblas_ddot macos-9 is not
#ifdef HAVE_CBLAS_H_BUG
extern "C" {
#define FF_VERSION VERSION
#undef VERSION
#include <cblas.h>
//#undef VERSION
//#define VERSION VERSION
}
#define WITHBLAS 1
#elif HAVE_VECLIB_CBLAS_BUG
#include <vecLib/cblas.h>
#define WITHBLAS 1
#endif
#endif
#ifdef WITHBLAS
template<class R> inline R blas_sdot(const int n,const R *sx,const int incx,const R *sy,const int  incy)
{
  R s=R();

  if(incx == 1 && incy == 1)
   for (int k = 0; k< n; k++)
    s += *sx++ * * sy++;
  else
   for (int k = 0; k< n; k++, sx += incx, sy += incy)
    s += *sx * *sy;
   return  s;
}

template<class R> inline R blas_sdot( int n, R *sx, int incx, R *sy, int  incy)
{
    R s=R();

    if(incx == 1 && incy == 1)
        for (int k = 0; k< n; k++)
            s += *sx++ * * sy++;
    else
        for (int k = 0; k< n; k++, sx += incx, sy += incy)
            s += *sx * *sy;
    return  s;
}
// OpenBlas PB with constant  remove const ....
template<> inline float blas_sdot(const int n,  float *sx, int incx, float *sy, int  incy)
{
    return cblas_sdot(n,sx,incx,sy,incy);
}
template<> inline double blas_sdot( int n,  double *sx, int incx, double *sy, int  incy)
{
    return cblas_ddot(n,sx,incx,sy,incy);
}

#ifdef OPENBLAS_CONFIG_H
typedef  openblas_complex_double *BLAS_ptr_complex16;
#else
typedef void *BLAS_ptr_complex16;

#endif
template<> inline  complex<double> blas_sdot( int n,  complex<double> *sx, int incx, complex<double> *sy, int  incy)
{
    complex<double> s;
    cblas_zdotu_sub(n,( double *)sx,incx,( double *)sy,incy,(BLAS_ptr_complex16)&s);
    return s;
}
//template<> inline  complex<float> blas_sdot( int n,  complex<float> *sx, int incx, complex<float> *sy, int  incy)
//{
//    complex<float> s;
//    cblas_cdotu_sub(n,( void *)sx,incx,( void *)sy,incy,(BLAS_ptr_complex8)&s);
//    return s;
//}

#endif
//  end modif FH



using std::numeric_limits;

//  -----------
template<class FElement>
inline int  BuildMEK_KK(const int l,int *p,int *pk,int *pkk,const FElement * pKE,const FElement*pKKE)
{
  // routine which find common  dof on to adjacent element pKE and pKKE and make the link ..
  // if pKKE== 0   then no  adj element
  // the idea is find common dof, but this work only if all dot a different
  // in on elemnt, so we can have a bug
  //  in case of periodic boundary condition ..
  // not correct ... F.Hecht ...

  // -----
  // routine  build  les array p, pk,pkk
  // which return number of df int 2 element pKE an pKKE
  // max l size of array p, pk, pkk
  // p[i] is the global number of freedom
  // pk[i] is is the local number in pKE ( -1 if not in pKE element)
  // pkk[i] is is the local number in pKKE ( -1 if not in pKKE element)
     const FElement (*pK[2])={pKE,pKKE};

     int ndf=0; // number of dl
     int * qk=pk, *qkk=pkk;
     for (int k=0;k<2;k++)
      if(pK[k])
       {
         if(k) Exchange(qk,qkk);
         const FElement& FEK=*pK[k];
         int nbdf =FEK.NbDoF();

         for (int ii=0;ii<nbdf;ii++)
          {
           p[ndf] = 2*FEK(ii)+k; // copy the numbering
           qk[ndf] = ii;
           qkk[ndf++] = -1;
          } // end for ii
         }
      ffassert(ndf <=l);
    int bug=0;
   // compression suppression des doublons
    // attention un df peu aparaitre 2 fois (CL period) dans un element ..
       Fem2D::HeapSort(p,pk,pkk,ndf);
       int k=0;
        for(int ii=1;ii<ndf;++ii)
          if (p[k]/2==p[ii]/2) // doublons k,kk
            {
              if (pkk[ii]>=0) pkk[k]=pkk[ii];
              if (pk[ii]>=0) pk[k]=pk[ii];
              assert(pk[k] >=0 && pkk[k]>=0);
            }
           else { // copy
              if(p[k]==p[ii]) bug++;
              p[++k] =p[ii];
              pk[k]=pk[ii];
              pkk[k]=pkk[ii];
             }
        ndf=k+1;
  for(int ii=0;ii<ndf;++ii)
      p[ii]= p[ii]/2;// clean pp to revome bug(CL period)
    if( bug && pKKE) {
        static int count =0;
        if( count++ < 2 && verbosity )
        {
        cerr << "  May be a Bug in BuildMEK_KK , the code is not safe , periodic boundary condition  on 1 element . " << bug <<  endl;
        cerr << "  They is a problem   in this case (I am not sure) F.H.  ????" << endl;
            cerr << "   exempt if the associed  matrix voefficient is 0. "<< endl;
        }
        //ffassert(0); // bof bof ... remove of case of jump in internal edge ...
    }
   return ndf;
} //  BuildMEK_KK

template<class R,class FES>
void MatriceElementairePleine<R,FES>::call(int k,int ie,int label,void * stack,void *B) {
 //   cout << " BUG ?? " << k << " " << *((long*) (void *)  (this->data)-1) <<endl;
  for (int i=0;i<this->lga;i++)
     this->a[i]=0;
  if(this->onFace)
    {
     throwassert(faceelement);
     const Mesh &Th(this->Vh.Th);

     int iie=ie,kk=Th.ElementAdj(k,iie);
     if(kk==k|| kk<0) kk=-1;
     if ( &this->Vh == &this->Uh)
      {
       FElement Kv(this->Vh[k]);
       if(kk<0)
        { // return ; // on saute ????  bof bof
	  this->n=this->m=BuildMEK_KK<FElement>(this->lnki,this->ni,this->nik,this->nikk,&Kv,0);
         int n2 =this->m*this->n;
         for (int i=0;i<n2;i++) this->a[i]=0;
         faceelement(*this,Kv,Kv,Kv,Kv,this->data,ie,iie,label,stack,reinterpret_cast<Rd*>(B));
        }
        else
        {
         FElement KKv(this->Vh[kk]);
         this->n=this->m=BuildMEK_KK<FElement>(this->lnki,this->ni,this->nik,this->nikk,&Kv,&KKv);


         faceelement(*this,Kv,KKv,Kv,KKv,this->data,ie,iie,label,stack,reinterpret_cast<Rd*>(B));

        }
      }
     else
      {
          throwassert(faceelement);
          const Mesh &Th(this->Vh.Th);

          int iie=ie,kk=Th.ElementAdj(k,iie);
          if(kk==k|| kk<0) kk=-1;
          if ( &this->Vh.Th  == &this->Uh.Th )
          {
              FElement Kv(this->Vh[k]);
              FElement Ku(this->Uh[k]);
              if(kk<0)
              { // return ; // on saute ????  bof bof
                  this->n=BuildMEK_KK<FElement>(this->lnki,this->ni,this->nik,this->nikk,&Kv,0);
                  this->m=BuildMEK_KK<FElement>(this->lnkj,this->nj,this->njk,this->njkk,&Ku,0);
                  int n2 =this->m*this->n;
                  for (int i=0;i<n2;i++) this->a[i]=0;
                  faceelement(*this,Ku,Ku,Kv,Kv,this->data,ie,iie,label,stack,reinterpret_cast<Rd*>(B));
              }
              else
              {
                  FElement KKv(this->Vh[kk]);
                  FElement KKu(this->Uh[kk]);
                  this->n=BuildMEK_KK<FElement>(this->lnki,this->ni,this->nik,this->nikk,&Kv,&KKv);
                  this->m=BuildMEK_KK<FElement>(this->lnkj,this->nj,this->njk,this->njkk,&Ku,&KKu);

                  faceelement(*this,Ku,KKu,Kv,KKv,this->data,ie,iie,label,stack,reinterpret_cast<Rd*>(B));//  correct may 2023 FH
                  

              }
          }
          else
          {
              ERREUR(" No DG on diff meshes : ???? /Impossible  TO DO  (see F. hecht) ", 0);
              ffassert(0); // a faire F. Hecht desole
          }

      }
   }
  else {
  throwassert(element);
  const FElement&Kv(this->Vh[k]);
  int nbdf =Kv.NbDoF();
  for (int i=0;i<nbdf;i++)
     this->ni[i] = Kv(i); // copy the numbering
  this->m=this->n=nbdf;

  if(this->ni != this->nj) { //
    const FElement&Ku(this->Uh[k]);
    int nbdf =Ku.NbDoF();
    for (int i=0;i<nbdf;i++)
      this->nj[i] = Ku(i); // copy the numbering
     this->m=nbdf;
     int n2 =this->m*this->n;
     for (int i=0;i<n2;i++) this->a[i]=0;
     element(*this,Ku,Kv,this->data,ie,label,stack,reinterpret_cast<Rd*>(B));
  }
  else
    {
     int n2 =this->m*this->n;
     for (int i=0;i<n2;i++) this->a[i]=0;
     element(*this,Kv,Kv,this->data,ie,label,stack,reinterpret_cast<Rd*>(B));
   // call the elementary mat
    }
  }
}

template<class R,class FES>
void MatriceElementaireSymetrique<R,FES>::call(int k,int ie,int label,void * stack,void  *B) {
  // mise a zero de la matrice elementaire, plus sur
  for (int i=0;i<this->lga;i++)
    this->a[i]=0;
  if(this->onFace)
    {
      ffassert(0); // a faire
    }
  else {

    if (k< this->Uh.Th.nt)
      {
	throwassert(element);
	const FElement K(this->Uh[k]);
	int nbdf =K.NbDoF();
	for (int i=0;i<nbdf;i++)
	  this->ni[i] = K(i); // copy the numbering
	this->m=this->n = nbdf;

	element(*this,K,this->data,ie,label,stack,static_cast<Rd*>(B));
      }// call the elementary mat
    else
      {
	ffassert(0); // remove code for the 3d
      }
  }
}



#endif
