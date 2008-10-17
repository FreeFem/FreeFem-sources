#ifndef  MatriceCreuse_tpl_
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
#ifdef HAVE_CBLAS_H
extern "C" {
#include <cblas.h> 
}
#define WITHBLAS 1
#elif HAVE_VECLIB_CBLAS_H
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
template<> inline float blas_sdot(const int n,const  float *sx,const int incx,const float *sy,const int  incy)
{
  return cblas_sdot(n,sx,incx,sy,incy);
}
template<> inline double blas_sdot(const int n,const  double *sx,const int incx,const double *sy,const int  incy)
{
  return cblas_ddot(n,sx,incx,sy,incy);
}
template<> inline  complex<double> blas_sdot(const int n,const  complex<double> *sx,const int incx,const complex<double> *sy,const int  incy)
{
  complex<double> s;
   cblas_zdotu_sub(n,(const void *)sx,incx,(const void *)sy,incy,( void *)&s);
  return s;
}
template<> inline  complex<float> blas_sdot(const int n,const  complex<float> *sx,const int incx,const complex<float> *sy,const int  incy)
{
  complex<float> s;
   cblas_zdotu_sub(n,(const void *)sx,incx,(const void *)sy,incy,( void *)&s);
  return s;
}
#endif
//  end modif FH
using Fem2D::HeapSort;
using std::numeric_limits;

//  -----------
template<class FElement>
inline int  BuildMEK_KK(const int l,int *p,int *pk,int *pkk,const FElement * pKE,const FElement*pKKE)
{
  // routine  build  les array p, pk,pkk 
  // which return number of df int 2 element pKE an pKKE
  // max l size of array p, pk, pkk
  // p[i] is the global number of freedom
  // pk[i] is is the local number in pKE ( -1 if not in pKE element)
  // pkk[i] is is the local number in pKKE ( -1 if not in pKKE element)
  //  remark, if pKKE = 0 => 
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
           p[ndf] = FEK(ii); // copy the numbering 
           qk[ndf] = ii;
           qkk[ndf++] = -1;
          } // end for ii
         } 
      ffassert(ndf <=l);
   // compression suppression des doublons
       Fem2D::HeapSort(p,pk,pkk,ndf);
       int k=0;  
        for(int ii=1;ii<ndf;++ii)
          if (p[k]==p[ii]) // doublons 
            { 
              if (pkk[ii]>=0) pkk[k]=pkk[ii];
              if (pk[ii]>=0) pk[k]=pk[ii];
              assert(pk[k] >=0 && pkk[k]>=0);
            }
           else { // copy 
              p[++k] =p[ii];
              pk[k]=pk[ii];
              pkk[k]=pkk[ii];
             }
        ndf=k+1; 
         
   return ndf;
} //  BuildMEK_KK

template<class R,class FES>
void MatriceElementairePleine<R,FES>::call(int k,int ie,int label,void * stack) {
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
	  this->n=this->m=BuildMEK_KK<FElement>(this->lnk,this->ni,this->nik,this->nikk,&Kv,0);
         int n2 =this->m*this->n; 
         for (int i=0;i<n2;i++) this->a[i]=0;
         faceelement(*this,Kv,Kv,Kv,Kv,this->data,ie,iie,label,stack);
        }
        else
        {
         FElement KKv(this->Vh[kk]);
         this->n=this->m=BuildMEK_KK<FElement>(this->lnk,this->ni,this->nik,this->nikk,&Kv,&KKv);
        
         
         faceelement(*this,Kv,KKv,Kv,KKv,this->data,ie,iie,label,stack);

        }
      }
     else 
      {
        ERREUR("A FAIRE/ TO DO  (see F. hecht) ", 0); 
        ffassert(0); // a faire F. Hecht desole 
       
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
     element(*this,Ku,Kv,this->data,ie,label,stack);
  }
  else 
    {
     int n2 =this->m*this->n;
     for (int i=0;i<n2;i++) this->a[i]=0;
     element(*this,Kv,Kv,this->data,ie,label,stack);
   // call the elementary mat 
    }  
  }  
}

template<class R,class FES>
void MatriceElementaireSymetrique<R,FES>::call(int k,int ie,int label,void * stack) {
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
	
	element(*this,K,this->data,ie,label,stack); 
      }// call the elementary mat 
    else
      {
	ffassert(0); // remove code for the 3d 
	/*
	throwassert(mortar);
	{
	  const FMortar K(&(this->Uh),k);
	  int nbdf = K.NbDoF();
	  for (int i=0;i<nbdf;i++)
	    this->ni[i] = K(i); // copy the numbering 
	  this->m=this->n = nbdf; 
	  // mise a zero de la matrice elementaire, plus sur
	  
	  mortar(*this,K,stack);}
	*/
      }
  }
}
  
template<class R>
MatriceProfile<R>::~MatriceProfile() {
  if(!this->dummy) 
    { //cout << " del mat profile " << endl ;
    if (U && (U !=L))  delete [] U;
    if (D)  delete [] D;
    if (L)  delete [] L;
    if (pU && (pU != pL)) delete [] pU;
    if (pL) delete [] pL;
    //cout << " dl de MatriceProfile " << this << endl;
    }
}
template<class R>
int MatriceProfile<R>::size() const {
  int s = sizeof(MatriceProfile<R>);
  if (D) s += this->n*sizeof(R);
  if (pL) s += this->n*sizeof(int);
  if (pU && (pU != pL)) s += this->n*sizeof(int);
  if (L) s += pL[this->n]*sizeof(int);
  if (U && (U != L)) s += pU[this->n]*sizeof(int);
  return s;
}
/*
template<class R>
int MatriceProfile<R>::MatriceProfile(const MatriceProfile<RR> & A )
  : MatriceCreuse<R>(A.n,A.m,0)
  {
    
    typefac=A.typefac;
    pL=  docpy<int,int>(A.pL,n+1);
    D = docpy<R,RR>(A.D,n);
    if ( A.pL == A.pU ) pU=pL;
    else pU=  docpy<int,int>(A.pU,m+1);
    
      L= docpy<R,RR>(A.L,pL[n]);
      
    if ( A.L == A.U ) U=L;
    else  U= docpy<R,RR>(A.U,pU[m]);
    
  
  }*/
template<class R>
  MatriceMorse<R> *MatriceProfile<R>::toMatriceMorse(bool transpose,bool copy) const 
  {
  // A FAIRE;
    ffassert(0); // TODO
   return 0;
  }
  
 inline pair<int,int> ij_mat(bool trans,int ii00,int jj00,int i,int j) {
  // warning trans sub  matrix and not the block. 
  return trans ? make_pair<int,int>(j+ii00,i+jj00)
                :  make_pair<int,int>(i+ii00,j+jj00) ; }
 
template<class R>
bool MatriceProfile<R>::addMatTo(R coef,std::map< pair<int,int>, R> &mij,bool trans,int ii00,int jj00,bool cnj)
{
   double eps0=numeric_limits<double>::min();
 if( norm(coef)<eps0) return  L == U ;
 int i,j,kf,k;
  if(D)
   for( i=0;i<this->n;i++)
    if( norm(D[i])>eps0)
     mij[ij_mat(trans,ii00,jj00,i,i)] += coef*(cnj? conj(D[i]) : D[i]);
   else
   for(int i=0;i<this->n;i++) // no dia => identity dai
     mij[ij_mat(trans,ii00,jj00,i,i)] += coef;
     
 if (L && pL )    
   for (kf=pL[0],i=0;  i<this->n;   i++  )  
     for ( k=kf,kf=pL[i+1], j=i-kf+k;   k<kf; j++,  k++  )
        if(norm(L[k])>eps0)
        mij[ij_mat(trans,ii00,jj00,i,j)]= coef*(cnj? conj(L[k]) : L[k]);
 if (U && pU)     
   for (kf=pU[0],j=0;  j<this->m;  j++)  
     for (k=kf,kf=pU[j+1], i=j-kf+k;   k<kf; i++,  k++  )
      if(norm(U[k])>eps0)
        mij[ij_mat(trans,ii00,jj00,i,j)]= coef*(cnj? conj(U[k]) : U[k]);
 return L == U ; // symetrique               
}
template<class R>
MatriceProfile<R>::MatriceProfile(const int nn,const R *a)
  :MatriceCreuse<R>(nn,nn,0),typefac(FactorizationNO)
{
  int *pf = new int [this->n+1];
  int i,j,k;
  k=0;
  for (i=0;i<=this->n;k+=i++)
    {
      pf[i]=k;
      //  cout << " pf " << i<< " = " << k  << endl;
    }
  ffassert( pf[this->n]*2 == this->n*(this->n-1));
  pU = pf; // pointeur profile U
  pL = pf; // pointeur profile L
  U = new R[pf[this->n]];
  L = new R[pf[this->n]];
  D = new R[this->n];  
  const R *aij=a;
  for (i=0;i<this->n;i++)
    for (j=0;j<this->n;j++)
      if      (j<i)   L[pL[i+1]-i+j] = *aij++;
      else if (j>i)   U[pU[j+1]-j+i] = *aij++;
      else            D[i] = *aij++;
}

template<class R>
template<class FESpace>
MatriceProfile<R>::MatriceProfile(const FESpace & Vh,bool VF) 
  :MatriceCreuse<R>(Vh.NbOfDF,Vh.NbOfDF,0),typefac(FactorizationNO)
{
   // for galerkine discontinue ....
   // VF : true=> Finite Volume matrices (change the stencil) 
   // VF = false => Finite element 
   // F. Hecht nov 2003
   // -----
  this->dummy=0;
  this->n = this->m = Vh.NbOfDF;
  int i,j,k,ke,ie,mn,jl,iVhk;
  int itab,tabk[5]; 
  int *pf = new int [this->n+1];
  for (i=0;i<this->n;i++)  pf[i]=0;
  for (ke=0;ke<Vh.NbOfElements;ke++)
    { 
      itab=0;
      tabk[itab++]=ke;
      if(VF) itab += Vh.Th.GetAllElementAdj(ke,tabk+itab);
      tabk[itab]=-1;    
      mn = this->n;
      for( k=tabk[ie=0]; ie <itab; k=tabk[++ie])
	{ iVhk=(int) Vh(k);
        for (jl=0;jl<iVhk;jl++) // modif Oct 2008 valgrind
	  { 
	    j=Vh(k,jl) ;
	    mn = Min ( mn , Vh.FirstDFOfNode(j) ) ;}
        }
       //for( k=tabk[ie=0]; ie <itab; k=tabk[++ie])
        { k=ke; // bof bof a verifier finement .... FH
	  iVhk=(int) Vh(k);  
        //for (j=Vh(k,jl=0);jl<(int) Vh(k);j=Vh(k,++jl)) 
	for (jl=0;jl<iVhk;jl++) // modif Oct 2008 valgrind
	     {
		 j=Vh(k,jl);
	      int df1 = Vh.LastDFOfNode(j);
	      for (int df= Vh.FirstDFOfNode(j);  df < df1; df++  )
	       pf[df] = Max(pf[df],df-mn);
	     }
	     }
    }
  int l =0;
  for (i=0;i<this->n;i++)  {int tmp=l;l += pf[i]; pf[i]=tmp;}
  pf[this->n] = l;
  if(verbosity >3) 
    cout << "  -- SizeOfSkyline =" <<l << endl;

  pU = pf; // pointeur profile U
  pL = pf; // pointeur profile L
  D = 0; // diagonal
  U = 0; // upper part
  L = 0; // lower part 
}

template<class R>
void MatriceProfile<R>::addMatMul(const KN_<R> &x,KN_<R> &ax) const 
{if (x.n!= this->n ) ERREUR(MatriceProfile MatMut(xa,x) ," longueur incompatible x (in) ") ;
 if (ax.n!= this->n ) ERREUR(MatriceProfile MatMut(xa,x) ," longueur incompatible ax (out)") ;
 int i,j,k,kf;
 ffassert(this->n == this->m);
 if (D) 
   for (i=0;i<this->n;i++) 
     ax[i] += D[i]*x[i];
 else
   for (i=0;i<this->n;i++) // no dia => identyty dai
     ax[i] +=x[i];
      
 if (L && pL )    
   for (kf=pL[0],i=0;  i<this->n;   i++  )  
     for ( k=kf,kf=pL[i+1], j=i-kf+k;   k<kf; j++,  k++  )
       ax[i] += L[k]*x[j],throwassert(i>=0 && i <this->n && j >=0 && j < this->m && k>=0 && k < pL[this->n]);
       
 if (U && pU)     
   for (kf=pU[0],j=0;  j<this->m;  j++)  
     for (k=kf,kf=pU[j+1], i=j-kf+k;   k<kf; i++,  k++  )
       ax[i] += U[k]*x[j],throwassert(i>=0 && i <this->n && j >=0 && j < this->m &&  k>=0 && k < pU[this->n]);
 

}


template<class R>
void MatriceProfile<R>::operator=(const R & v) {
  if(v!=R())
    { cerr << " Mise a zero d'une matrice MatriceProfile<R>::operator=(R v) uniquement v=" << v << endl;
    throw(ErrorExec("exit",1));
    }
  typefac = FactorizationNO;
  delete [] U;
  delete [] L;
  delete [] D;
  U=L=D=0;
}
template<class R>
MatriceCreuse<R>  & MatriceProfile<R>::operator +=(MatriceElementaire<R> & me) {
  int il,jl,i,j,k;
  int * mi=me.ni, *mj=me.nj;
  if (!D)  // matrice vide 
    { D  = new R[this->n];
    L  = pL[this->n] ? new R[pL[this->n]] :0 ;
    for (i =0;i<this->n;i++) D[i] =0;
    for (k =0;k<pL[this->n];k++) L[k] =0;
    switch (me.mtype) {
    case MatriceElementaire<R>::Full :     
      U  = pU[this->n] ? new R[pU[this->n]] : 0;
      for (k =0;k<pU[this->n];k++) U[k] =0;
      break;
    case MatriceElementaire<R>::Symmetric :     
      U = L; 
      break;
    default:
      cerr << "Big bug type MatriceElementaire unknown" << (int) me.mtype << endl;
      throw(ErrorExec("exit",1));
      break; 
    }
    }
  R * al = me.a; 
  switch (me.mtype) {
  case MatriceElementaire<R>::Full : //throwassert(L !=U);
    for (il=0; il<me.n; ++il) // modif overflow FH win32  oct 2005
     { i=mi[il];
      for ( jl=0; jl< me.m ; ++jl,++al)  // modif overflow FH
        { j=mj[jl] ;
	  if      (j<i)  L[ pL[i+1] - (i-j) ] += *al;
	  else if (j>i)  U[ pU[j+1] - (j-i) ] += *al;
	  else           D[i] += *al;}}
    break;
     
  case MatriceElementaire<R>::Symmetric : //throwassert(L ==U);   
    for (il=0; il<me.n; ++il) // modif overflow FH win32
      { i=mi[il];
      for (jl=0;jl<= il;++jl)
       { j=mj[jl]  ;
	 if      (j<i)  L[ pL[i+1] - (i-j) ] += *al++;
	 else if (j>i)  U[ pU[j+1] - (j-i) ] += *al++;
	 else           D[i] += *al++;}}
    break;
  default:
    cerr << "Big bug type MatriceElementaire unknown" << (int) me.mtype << endl;
    exit(1);
    break; 
  }      
  return *this;
} 

template<class R>
ostream& MatriceProfile<R>::dump (ostream& f) const 
{f<< " matrix skyline " << this->n << '\t' << this->m << '\t' ;
 f <<  "  this " << endl;
 f << " pL = " << pL << " L ="  << L << endl
   << " pU = " << pU << " U ="  << U << endl
   << " D = " << D << endl;
 if ( (pL == pU) &&  (U == L) )
   if (pL && L) 
     {f << " skyline symmetric " <<endl;
     int i,j,k;
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
     int i,k;
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
template<class R>
void MatriceProfile<R>::cholesky(double eps) const {
  double eps2=eps*eps;
  R  *ij , *ii  , *ik , *jk , xii;
  int i,j,k;
  if (L != U) ERREUR(factorise,"Skyline matrix non symmetric");
  U = 0; // 
  typefac = FactorizationCholeski;
  if ( norm(D[0]) <= 1.0e-60)
      ERREUR(cholesky,"pivot (" << 0 << ")= " << D[0] )
  
  D[0] = sqrt(D[0]); 
  ij = L ; // pointeur sur le terme ij de la matrice avec j<i 
  for (i=1;i<this->n;i++) // boucle sur les lignes 
    { ii = L+pL[i+1]; // pointeur sur le terme fin de la ligne +1 =>  ij < ii;
    xii = D[i] ; 
    for ( ; ij < ii ; ij++) // pour les j la ligne i
      { j = i -(ii - ij); 
      k = Max( j - (pL[j+1]-pL[j]) ,  i-(pL[i+1]-pL[i]) ); 
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
    // cout << norm(xii) << " " << Max(eps2*norm(D[i]),1.0e-60) << " " << sqrt(xii) <<endl;
    if ( norm(xii) <= Max(eps2*norm(D[i]),1.0e-60)) 
      ERREUR(cholesky,"pivot (" << i << ")= " << xii << " < " << eps*abs(D[i]))
    D[i] = sqrt(xii);
    }
}
template<class R>
void MatriceProfile<R>::crout(double eps) const  {
  R  *ij , *ii  , *ik , *jk , xii, *dkk;
  int i,j,k;
  double eps2=eps*eps;
  if (L != U) ERREUR(factorise,"Skyline matrix  non symmetric");
  U = 0; // 
  typefac = FactorizationCrout;
   
  ij = L ; // pointeur sur le terme ij de la matrice avec j<i 
  for (i=1;i<this->n;i++) // boucle sur les lignes 
    { ii = L+pL[i+1]; // pointeur sur le terme fin de la ligne +1 =>  ij < ii;
    xii = D[i] ; 
    for ( ; ij < ii ; ij++) // pour les j la ligne i
      { j = i -(ii - ij); 
      k = Max( j - (pL[j+1]-pL[j]) ,  i-(pL[i+1]-pL[i]) ); 
      ik =  ii - (i - k); 
      jk =  L + pL[j+1] -(j - k); 
      dkk = D + k;
      k = j - k ; 
      R s=-*ij;
      while ( k-- ) s += *ik++ * *jk++ * *dkk++;  
      *ij = -s/ *dkk ; // k = j ici 

      xii -= *ij * *ij * *dkk;
      }
    if (norm(xii) <= Max(eps2*norm(D[i]),1.0e-60))
      ERREUR(crout,"pivot (" << i << " )= " << abs(xii)<< " <= " << eps*abs(D[i]) << " eps = " << eps)
	D[i] = xii;
    }
}
template<class R>
void MatriceProfile<R>::LU(double eps) const  {
  R s,uii;
  double eps2=eps*eps;
  int i,j;
  if (L == U && ( pL[this->n]  || pU[this->n] ) ) ERREUR(LU,"matrix LU  symmetric");
  if(verbosity>3)
  cout << " -- LU " << endl;
  typefac=FactorizationLU;

  for (i=1;i<this->n;i++) // boucle sur les sous matrice de rang i 
    { 
      // for L(i,j)  j=j0,i-1
      int j0 = i-(pL[i+1]-pL[i]);
      for ( j = j0; j<i;j++)
        {           
          int k0 = Max(j0,j-(pU[j+1]-pU[j]));
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
	  int k0 = Max(j0,j-pL[j+1]+pL[j]);
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
      int k0 = i-Min(pL[i+1]-pL[i],pU[i+1]-pU[i]);
      R *Lik = L + pL[i+1]-i+k0; // lower
      R *Uki = U + pU[i+1]-i+k0; // upper
      s =0;
#ifdef WITHBLAS
       s = blas_sdot(i-k0,Lik,1,Uki,1);
#else
      for (int k=k0;k<i;k++) // k < i < i ;
	s += *Lik++ * *Uki++ ;     // a(i,k)*a(k,i);
#endif
      // cout << " k0 " << k0 << " i = " << i << " " <<  s << endl;
      uii = D[i] -s;
      
      if (norm(uii) <= Max(eps2*norm(D[i]),1.0e-30))
	ERREUR(LU,"pivot (" << i << " )= " << abs(uii) << " <= " << eps*abs(D[i]) << " eps = " << eps);     
      
      D[i] = uii;
      
    }
}


template<class R>
KN_<R> & operator/=(KN_<R> & x ,const MatriceProfile<R> & a) 
{
  // --------------------------------------------------------------------
  //   si La diagonal D n'existe pas alors on suppose 1 dessus (cf crout)
  // --------------------------------------------------------------------
  R * v = &x[0];
  int n = a.n;  
  if (x.n != n ) 
    ERREUR (KN_<R> operator/(MatriceProfile<R>),"  matrice et KN_<R> incompatible");
  const R *ij ,*ii, *ik, *ki;
  R *xk,*xi;
  int i;
  switch (a.typefac) {
  case FactorizationNO:
    if (a.U && a.L) {cerr << "APROGRAMMER (KN_<R><R>::operator/MatriceProfile)";throw(ErrorExec("exit",2));}
   
    if ( a.U && !a.L ) 
      { // matrice triangulaire superieure
	// cout << " remonter " << (a.D ? "DU" : "U") << endl;
	ki = a.U + a.pU[n]; 
	i = n;
	while ( i-- )
	  { ii = a.U + a.pU[i];
          xi= xk  = v +  i ;
          if (a.D) *xi /= a.D[i];// pour crout ou LU
          while ( ki > ii) 
	    *--xk  -=  *--ki *  *xi ; 
	  }
      }
    else if  ( !a.U && a.L ) 
      { // matrice triangulaire inferieure
	// cout << " descente "  <<( a.D ? "LD" : "L" ) <<endl;
	ii = a.L;
	for (i=0; i<n; i++)
	  { ij = ik = (a.L + a.pL[i+1]) ;  // ii =debut,ij=fin+1 de la ligne 
          xk = v + i;
          R ss = v[i]; 
          while ( ik > ii) 
	    ss -= *--ik * *--xk ; 
          if ( a.D) ss /= a.D[i];// pour crout ou LU
          v[i] = ss ;
          ii = ij;
	  }
      }
    else if (a.D) 
      { // matrice diagonale
	// cout << " diagonal D" <<endl;
	for (i=0;i<n;i++) 
	  v[i]=v[i]/a.D[i];
      }
    break;
  case FactorizationCholeski:  
    //     cout << " FactorizationChosleski" << endl;
    x /= a.ld();
    x /= a.ldt();   
    break;
  case FactorizationCrout:
    //   cout << " FactorizationCrout" << endl;
    x /= a.l();
    x /= a.d();
    x /= a.lt();
    break;
  case FactorizationLU:
    //  cout << " FactorizationLU" << endl;
    x  /= a.l();
    x  /= a.du();
    break;
    /*   default:
	 ERREUR  (operator /=(MatriceProfile," Error unkown type of Factorization  =" << typefac);
    */
  }
  return x;
}

template <class R> 
 MatriceMorse<R>::MatriceMorse(KNM_<R> & A,double tol)
    :MatriceCreuse<R>(A.N(),A.M(),false),solver(0) 
      {
  double tol2=tol*tol;    
  symetrique = false;
  this->dummy=false;
  int nbcoeff=0;
  for(int i=0;i<this->n;i++)
    for(int j=0;j<this->m;j++)
      if(norm(A(i,j))>tol2) nbcoeff++;

  nbcoef=nbcoeff;
  nbcoeff=Max(nbcoeff,1); // pour toujours alloue quelque chose FH Bug dans CheckPtr
  a=new R[nbcoeff] ;
  lg=new int [this->n+1];
  cl=new int [nbcoeff];
  nbcoeff=0;
  R aij;
  for(int i=0;i<this->n;i++)
   { 
    lg[i]=nbcoeff;
    for(int j=0;j<this->m;j++)
     
      if(norm(aij=A(i,j))>tol2)
       {
         cl[nbcoeff]=j;
         a[nbcoeff]=aij;
         nbcoeff++;
       }
    }
   lg[this->n]=nbcoeff;

  
}
template <class R> 
 MatriceMorse<R>::MatriceMorse(const int  nn,const R *aa)
    :MatriceCreuse<R>(nn),solver(0) 
      {
  symetrique = true;
  this->dummy=false;
  this->n=nn;
  nbcoef=this->n;
  a=new R[this->n] ;
  lg=new int [this->n+1];
  cl=new int [this->n];
  for(int i=0;i<this->n;i++)
   {
    lg[i]=i;
    cl[i]=i;
    a[i]=aa[i];      
      }
lg[this->n]=this->n;
}

template<class R>
template<class K>
 MatriceMorse<R>::MatriceMorse(const MatriceMorse<K> & A)
   : MatriceCreuse<R>(A.n,A.m,A.dummy),nbcoef(A.nbcoef),      
     symetrique(A.symetrique),       
     a(new R[nbcoef]),
     lg(new int [this->n+1]),
     cl(new int[nbcoef]),
     solver(0)
{
  ffassert(a && lg &&  cl);
  for (int i=0;i<=this->n;i++)
    lg[i]=A.lg[i];
  for (int k=0;k<nbcoef;k++)
    {
      cl[k]=A.cl[k];
      a[k]=A.a[k];
    }
  
}



template <class R> 
int MatriceMorse<R>::size() const 
{
  return nbcoef*(sizeof(int)+sizeof(R))+ sizeof(int)*(this->n+1);
}

inline int WhichMatrix(istream & f)
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
	if( line.find("(Morse)")) 
	    return 2; // morse 
	else 
	    return 0; 
    }
  return 0;   
}
template <class R>
  MatriceMorse<R>::MatriceMorse(istream & f)
:  MatriceCreuse<R>(0,0,0),nbcoef(0),
a(0),
lg(0),
cl(0),

solver(0)
{
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
	 cout << "Read matrice: "<< k << " :"   << line << endl;
	k++;    
      }
      
      f >> this->n >> this->m >> symetrique >>nbcoef;
      if(verbosity>3)
      cout << " read mat: " <<  this->n << " " <<  this->m << " " << symetrique << " " << nbcoef <<endl;
      lg= new int [this->n+1];
      cl= new int[nbcoef];
      a= new R[nbcoef];
      ffassert(f.good() && lg && a && cl );
      int i,j,i0,j0;
      i0=-1;j0=2000000000;
      R aij;
      int imx=-2000000000, jmx=-2000000000;
      int imn= 2000000000, jmn= 2000000000;
      
      for (int k =0;k<nbcoef; ++k)
      {
	  f >> i >> j >> aij;
	  ffassert(f.good() );
	  i--;j--;
	  imx=max(imx,i);
	  imx=max(jmx,j);
	  imn=min(imn,i);
	  imn=min(jmn,j);
	  //cout << i << " " << j << " " << aij << endl;
	  if(i0!=i) {j0=-1;lg[i]=k;}
	  ffassert(i0<=i && j0<j);
	  lg[i+1]=k+1;
	  cl[k]=j;
	  a[k]=aij;
	  j0=j;i0=i;
      }
      ffassert( imx < this->n && jmx < this->m );
      ffassert( imn >=0 && jmn >=0);
      
}

template <class R> 
ostream& MatriceMorse<R>::dump(ostream & f) const 
{
  f << "# Sparse Matrix (Morse)  " << endl;
  f << "# first line: n m (is symmetic) nbcoef \n";
  f << "# after for each nonzero coefficient:   i j a_ij where (i,j) \\in  {1,...,n}x{1,...,m} \n";  
 
  f << this->n << " " << this->m << " " << symetrique << "  " << nbcoef <<endl;
  int k=lg[0];
  int pold= f.precision();
  for (int i=0;i<this->n;i++)
   { 
    
//    f << i << " : " << lg[i] <<","<< lg[i+1]-1 << " : " ;
    int ke=lg[i+1];
    for (;k<ke;k++)
      f << setw(9) << i+1 << ' ' << setw(9) << cl[k]+1 << ' ' << setprecision( 20) << a[k]<< '\n' ;
     // if (norm(a[k])) f  << cl[k] << " " << a[k]<< ", ";
     // else f  << cl[k] << " 0., " ;
   // f << endl;    
   }
   f.precision(pold);
  return f;
}
template <class R> 
inline R*  MatriceMorse<R>::pij(int i,int j) const 
 {
   if (! (i<this->n && j< this->m)) 
   throwassert(i<this->n && j< this->m);
   int i0=lg[i];
   int i1=lg[i+1]-1;
   while (i0<=i1) // dichotomie
    { 
      int im=(i0+i1)/2;
      if (j<cl[im]) i1=im-1;
      else if (j>cl[im]) i0=im+1;
      else return a+im;      
    }
   return 0;     
 }
template <class R>
template <class FESpace> 
void MatriceMorse<R>::Build(const FESpace & Uh,const FESpace & Vh,bool sym,bool VF)
{
  typedef typename FESpace::Mesh Mesh;
  
  // for galerkine discontinue ....
  // VF : true=> Finite Volume matrices (change the stencil) 
  // VF = false => Finite element 
  // F. Hecht nov 2003
  // -----
  symetrique = sym;
  this->dummy=false;
  a=0;
  lg=0;
  cl=0;
  //  bool same  = &Uh == & Vh;
  ffassert( &Uh.Th == &Vh.Th);  // same Mesh
  const Mesh & Th(Uh.Th);
  //int nbt = Th.nt;
  //int nbv = Th.nv;
  //int nbm = Th.NbMortars;
  int nbe = Uh.NbOfElements;
  int nbn_u = Uh.NbOfNodes;
  int nbn_v = Vh.NbOfNodes;
  
  KN<int> mark(nbn_v);
  KN<int> pe_u(nbn_u+1+Uh.SizeToStoreAllNodeofElement());
  //  les element du node i 
  // sont dans pe_u[k] pour k \in [ pe_u[i] , pe_u[i+1] [
  pe_u=0;
  for (int k=0;k<nbe;k++)
    { 
      int nbne=Uh(k);
      for (int in=0;in<nbne;in++)
        pe_u[(Uh(k,in)+1)]++;
   }
  int kk= nbn_u+1,kkk=kk;
  pe_u[0]=kk;
  for (int in1=1;in1<=nbn_u;in1++)
    { // in1 = in + 1
      kk += pe_u[in1];
      pe_u[in1] = kkk; // store the last of in 
      kkk=kk;
    } 
  if(verbosity>4) 
    cout <<" -- MatriceMorse<R>::Build " << kk << " " << nbn_u << " " << Uh.SizeToStoreAllNodeofElement() 
	 << " " <<  nbn_u+1+Uh.SizeToStoreAllNodeofElement() << endl;
  ffassert(kk== nbn_u+1+Uh.SizeToStoreAllNodeofElement());
  for (int k=0;k<nbe;k++)
    { 
      int nbne=Uh(k);
      for (int in=0;in<nbne;in++)
        pe_u[pe_u[(Uh(k,in)+1)]++] = k;
    }
  
  
  int color=0;
  mark=color++;
  lg = new int [this->n+1];
  ffassert(lg);
  for (int step=0;step<2;step++) 
    { 
      int ilg=0;
      lg[0]=ilg;
      int kij=0;
    for (int in=0;in<nbn_u;in++)
      {
	int nbj=0; // number of j
	int kijs=kij;
	// for all triangle contening node in
	for (int kk= pe_u[in];kk<pe_u[in+1];kk++)
	  {
	    int ke=pe_u[kk];// element of 
	    int tabk[10];
	    int ltab=0;
	    tabk[ltab++]=ke;
	    if( VF) // if Finite volume then add Triangle adj in stencil ...
	      ltab+= Th.GetAllElementAdj(ke,tabk+ltab);
	    tabk[ltab]=-1;
	    for(int ik=0,k=tabk[ik];ik<ltab;k=tabk[++ik])
	      {
		throwassert(k>=0 && k < nbe);
		int njloc = Vh(k);
		for (int jloc=0;jloc<njloc;jloc++)
		  { 
		    int  jn = Vh(k,jloc);
		    if (mark[jn] != color && (!sym ||  jn < in) ) 
		      {
			mark[jn] = color;
			int fdf=Vh.FirstDFOfNode(jn);
			int ldf=Vh.LastDFOfNode(jn);
			if (step)
			  for (int j=fdf;j<ldf;j++)
			    cl[kij++] = j;
			nbj += ldf-fdf;
		      }            
		  }} 
	  }
	int fdf=Uh.FirstDFOfNode(in);
	int ldf=Uh.LastDFOfNode(in);
	int kijl=kij;
	if (step)
	  {
	    HeapSort(cl+kijs,kij-kijs);
	    for (int i=fdf;i<ldf;i++)
	      { 
		if (i!=fdf) //  copy the ligne if not the first 
		  for (int k=kijs;k<kijl;k++)
		    cl[kij++]=cl[k]; 
		if (sym) // add block diag
		  for(int j=fdf;j<=i;j++)
		    cl[kij++]=j;            
		throwassert(kij==lg[i+1]);// verif           
	      }
	  }
	else
	  for (int i=fdf;i<ldf;i++)
	    { 
	      if (sym) ilg += ++nbj; // for the diag block
	      else ilg += nbj;             
	      lg[i+1]=ilg;
	    }
	color++; // change the color
      }
    if (step==0) { // do allocation 
      nbcoef=ilg;
      if (verbosity >3)
        cout << "  -- MorseMatrix: Nb coef !=0 " << nbcoef << endl;
      a = new R[nbcoef];
      cl = new int [nbcoef];}
    ffassert( a && cl);
    for (int i=0;i<nbcoef;i++) 
      a[i]=0;
    
   }
  
}
template<class R> inline void ConjArray( R  *v, int n) 
{
  for (int i=0;i<n;++i)
    v[i] = conj(v[i]);
}
template<> inline void ConjArray<double>(double *v, int n) {}
template<> inline void ConjArray<float>(float *v, int n) {}

template<class R>
 void  MatriceMorse<R>::dotransposition()
 {
   if(symetrique) return; 
   
   ffassert(this->dummy==false);  
   int *llg= new int[nbcoef];
   int *clg= new int[this->m+1];
   
   for (int i=0;i<this->n;i++)
     for (int k=lg[i];k<lg[i+1];k++)
        llg[k]=i;
 
  HeapSort(cl,llg,a,nbcoef);
  for(int k=0;k<this->m;k++)
    clg[k]=-1;

  // build new line end (old column)
  for(int k=0;k<nbcoef;k++)
    clg[cl[k]+1]=k+1;
      
   for(int kk=0, k=0;k<=this->m;k++)
   if (clg[k]==-1)
      clg[k]=kk;
    else kk=clg[k];
    
  clg[this->m]=nbcoef;
  // sort the new column (old line)
  for(int i=0;i<this->m;i++)  
    HeapSort(llg+clg[i],cl+clg[i],a+clg[i],clg[i+1]-clg[i]); 
  delete[] cl;
  delete[] lg;
  Exchange(this->n,this->m);       
  cl=llg;
  lg=clg;
  ConjArray(a,nbcoef);    
 }

template<class R>
 triplet<int,int,bool> BuildCombMat(std::map< pair<int,int>, R> & mij,const list<triplet<R,MatriceCreuse<R> *,bool> >  &lM,bool trans,int ii00,int jj00,bool cnj=false)
  {
    typedef typename list<triplet<R,MatriceCreuse<R> *,bool> >::const_iterator lconst_iterator;
    
    lconst_iterator begin=lM.begin();
    lconst_iterator end=lM.end();
    lconst_iterator i;
    
   // std::map< pair<int,int>, R> mij;
    
    int n=0,m=0;
    bool sym=true;
    for(i=begin;i!=end;i++++)
     {
	if(i->second) // M == 0 => zero matrix 
	{
	    MatriceCreuse<R> & M=*i->second;
	    bool transpose = i->third !=  trans;
	    ffassert( &M);
	    R coef=i->first;
	    if(verbosity>3)
		cout << "                BuildCombMat + " << coef << "*" << &M << " " << sym << "  t = " << transpose << " " <<  i->third << endl;
	    //  change to max FH dec 2007 to hard to satisfy
	   /* if (n==0)*/ { if(transpose) {m=max(m,M.n); n=max(n,M.m);} else{n=max(M.n,n); m=max(M.m,m);}}// Modif mars 2007 FH
	   /* else { if(transpose)  ffassert(n== M.m && m==M.n); else ffassert(n== M.n && m==M.m);}*/
	    sym = M.addMatTo(coef,mij,transpose,ii00,jj00,cnj) && sym;  
	}
     } 
    int nbcoef=mij.size();
    if(sym) nbcoef = (nbcoef+n)/2;

  // return new   MatriceMorse<R>(n,m,mij,sym);   
    return make_triplet(n,m,sym);
  }
  
template<class R>
  MatriceMorse<R> * BuildCombMat(const list<triplet<R,MatriceCreuse<R> *,bool> >  &lM,bool trans,int ii00,int jj00)
  {
   
    std::map< pair<int,int>, R> mij;
    triplet<int,int,bool> nmsym=BuildCombMat(mij,lM,trans,ii00,jj00);

   return new   MatriceMorse<R>(nmsym.first,nmsym.second,mij,nmsym.third);   
     
  }
template<class R>
bool MatriceMorse<R>::addMatTo(R coef,std::map< pair<int,int>, R> &mij,bool trans,int ii00,int jj00,bool cnj)
{
  double eps0=numeric_limits<double>::min();
  int i,j,k;
  if (symetrique)
   {
     for ( i=0;i<this->n;i++)
       for ( k=lg[i];k<lg[i+1];k++)
         {
           j=cl[k];
           R cij =  coef* ( cnj ? conj(a[k]) : a[k]);
           if(norm(cij)>eps0)
           {
            mij[ij_mat(trans,ii00,jj00,i,j)] += cij ;
           if (i!=j)
             mij[ij_mat(trans,ii00,jj00,j,i)] += cij;
           }
         }
           
   }
  else
   {
     for ( i=0;i<this->n;i++)
       for ( k=lg[i];k<lg[i+1];k++)
         {
           j=cl[k];
           R cij =  coef* ( cnj ? conj(a[k]) : a[k]);

           if(norm(cij)>eps0)
           mij[ij_mat(trans,ii00,jj00,i,j)] += cij;
         }
   }

return symetrique;
}

 
template<class R> 
template<class K>
MatriceMorse<R>::MatriceMorse(int nn,int mm, std::map< pair<int,int>, K> & m, bool sym):
  MatriceCreuse<R>(nn,mm,0),
  nbcoef(m.size()),symetrique(sym),
  a(new R[nbcoef]),
  lg(new int[nn+1]),
  cl(new int[nbcoef]),     
  solver(0)
{
     int k=0;
     bool nosym=!sym;
     typename std::map< pair<int,int>, R>::iterator iter=m.begin(), mend=m.end();
     //  remarque lg est croissant Bug trouver par 
     for(int i=0;i<=nn;i++) lg[i]=0; 
     while(iter!=mend)
      { 
        int i=iter->first.first;
        int j=iter->first.second;
        K & aij=iter->second;
        assert( i < nn && j < mm);
        if(j<=i || nosym)
        {
         cl[k]=j;
         a[k]=aij;
         lg[i+1]=++k;
        }
        ++iter;
       }
    // lg est croissant  on bouche les trou   
   for(int i=1;i<=nn;i++) lg[i]=Max(lg[i-1],lg[i]); 
      
   ffassert(nbcoef==k);  
  }

template<class RA>
 template<class RB,class RAB>
 void  MatriceMorse<RA>::prod(const MatriceMorse<RB> & B, MatriceMorse<RAB> & AB)
 {
   //  compute the s
  bool sym=this == & B &&symetrique;
  int *blg=B.lg;
  int *bcl=B.cl;
  ffassert(this->m==B.n); 
  bool delbl= B.symetrique;
  if (delbl)
    {
     int nn=B.n;
      blg = new int[nn+1];
     for (int i=0;i<B.n;i++)
         blg[i]=B.lg[i+1]-B.lg[i];
      blg[nn]=0;   
      
      for (int i=0;i<nn;i++)
        for (int k= B.lg[i];k<B.lg[i+1];k++)
          {  int j=B.cl[k];
              assert(j <= i);
             if (j!=i)  
               blg[j]++;
             }
             
      for (int i=1;i<=nn;i++)
       blg[i]+=blg[i-1];
      int nbnz = blg[nn];
      bcl= new int[nbnz];
      
      for (int i=0;i<B.n;i++)
        for (int k= B.lg[i];k<B.lg[i+1];k++)
          {  int j=B.cl[k];
             assert(j <= i);
             bcl[--blg[i] ]=j;
             if(i !=j)
               bcl[--blg[j]]=i;
          }
    }
   
   set<pair<int,int> > sij;
   double eps0=numeric_limits<double>::min();

     for (int i=0;i<this->n;i++)
       for (int k=lg[i];k<lg[i+1];k++)
         {    
           int j=cl[k];
           if(norm(a[k])<eps0) continue;
           int ii[2],jj[2];
           ii[0]=i;ii[1]=j;
           jj[0]=j;jj[1]=i;
           int kk=1;
           if(symetrique && i != j) kk=2;
           for (int ll=0;ll<kk;ll++)
            {
                int i=ii[ll];
                int j=jj[ll];
                for (int kkb=blg[j];kkb<blg[j+1];kkb++)
                  { 
                   int kz= bcl[kkb];
                   RB bjk;
                   if (B.symetrique && kz > j)
                     bjk=B(kz,j);
                   else
                      bjk=B(j,kz);
                   if( norm(bjk)>eps0 && (!sym || kz<=i))
                     sij.insert(make_pair(i,kz));
                  }
            }
           
         }
    int nn=this->n;
    int mm=B.m;
    int * llg=new int[nn+1];
    int * lcl=new int[sij.size()];  
    RAB * aa = new RAB[sij.size()];
    for(int i=0;i<=nn;i++)
        llg[i]=0;
        
    for (set<pair<int,int> >::iterator iter=sij.begin();iter!=sij.end();++iter)
      { 
        int i=iter->first;
	// int j=iter->second;
        llg[i]++;
       }
     for (int i=1;i<=nn;i++)
       llg[i]+=llg[i-1];
     ffassert(llg[this->n]==(long) sij.size());
     for (set<pair<int,int> >::iterator iter=sij.begin();iter!=sij.end();++iter)
      { 
        int i=iter->first;
        int j=iter->second;
       // cout << i << " , " << j << endl;
        lcl[--llg[i]]=j;
       }
     for(int i=0;i<nn;i++)  
       HeapSort(lcl+llg[i],llg[i+1]-llg[i]); 
       
     AB.n=nn;
     AB.m=mm;
     AB.N=nn;  // add missing jan 2008 FH
     AB.M=mm;  // add missing jan 2008 FH

     AB.lg=llg;
     AB.cl=lcl;
     AB.a=aa;        
     AB.nbcoef=sij.size();
     AB.symetrique=sym;
     AB.dummy=false;
     AB = RAB();
     for (int i=0;i<this->n;i++)
       for (int k=lg[i];k<lg[i+1];k++)
         {    
           int j=cl[k];
           RAB aij = a[k];
           if(norm(aij) <eps0 ) continue;
           int ii[2],jj[2];
           ii[0]=i;ii[1]=j;
           jj[0]=j;jj[1]=i;
           int kk=1;
           if(symetrique && i != j) kk=2;
           for (int ll=0;ll<kk;ll++)
            {
                int i=ii[ll];
                int j=jj[ll];
                for (int kb=blg[j];kb<blg[j+1];kb++)
                  { 
                   int k= bcl[kb];
                   RB bjk;
                   if (B.symetrique && k > j)
                     bjk=B(k,j);
                   else
                      bjk=B(j,k);
                //   cout << i << "," << "," << j << "," << k << " " << aij << " " << bjk << endl;
                   if( norm( bjk)> eps0  && (!sym || k<=i))
                       AB(i,k) += aij*bjk;
                  }
            }
           
         }

    if (delbl) {
      delete [] blg;
      delete [] bcl;
    }
     
     
 }

template<class R>
  void  MatriceMorse<R>::addMatMul(const KN_<R> &  x, KN_<R> & ax) const   
{
  int i,j,k;
  if( ! (this->n==ax.N() && this->m==x.N()))
    {cerr << " Err MatriceMorse<R>:  ax += A x" <<endl;
      cerr << " A.n " << this->n<< " !=  "<< ax.N() << " ax.n \n";
      cerr << " A.m " << this->m<< " != " <<x.N() << " x.n \n" ;
      ffassert(0); 
      abort();
    }
  if (symetrique)
   {
     for (i=0;i<this->n;i++)
       for (k=lg[i];k<lg[i+1];k++)
         {
           j=cl[k];
           ax[i] += a[k]*x[j];
           if (i!=j)
             ax[j] += a[k]*x[i];
         }
           
   }
  else
   {
     for (i=0;i<this->n;i++)
       for (k=lg[i];k<lg[i+1];k++)
         {
           j=cl[k];
           ax[i] += a[k]*x[j];
         }
   }
}

template<class R>
  void  MatriceMorse<R>::addMatTransMul(const KN_<R> &  x, KN_<R> & ax) const   
{
  int i,j,k;
  ffassert(this->m==ax.N());
  ffassert(this->n==x.N());  
  if (symetrique)
   {
     for (i=0;i<this->n;i++)
       for (k=lg[i];k<lg[i+1];k++)
         {
           j=cl[k];
           ax[j] += a[k]*x[i];
           if (i!=j)
             ax[i] += a[k]*x[j];
         }
           
   }
  else
   {
     for (i=0;i<this->n;i++)
       for (k=lg[i];k<lg[i+1];k++)
         {
           j=cl[k];
           ax[j] += a[k]*x[i];
         }
   }
}


template<class R>
MatriceMorse<R>  & MatriceMorse<R>::operator +=(MatriceElementaire<R> & me) {
  int il,jl,i,j;
  int * mi=me.ni, *mj=me.nj;
  if ((this->n==0) && (this->m==0))
   {
   
    //    if(verbosity>3)
    cout << " -- Morse Matrice is empt: let's build it" << endl;
    ffassert(0); 
    /*
    this->n=me.Uh.NbOfDF;
    this->m=me.Vh.NbOfDF;
    switch (me.mtype) {
    case MatriceElementaire<R>::Full : 
      Build(me.Uh,me.Vh,false);    
      break;
    case MatriceElementaire<R>::Symmetric :     
      Build(me.Uh,me.Vh,true);    
      break;
     default:
      cerr << "Big bug type MatriceElementaire is unknown" << (int) me.mtype << endl;
      throw(ErrorExec("exit",1));
      break; }     
    */
   }
  R * al = me.a; 
  R * aij;
  switch (me.mtype) { // modif FH overfloat in array mi and mj => trap on win32
  case MatriceElementaire<R>::Full : ffassert(!symetrique);
    for (il=0; il<me.n; ++il)  { i=mi[il]; 
      for ( jl=0; jl< me.m ; ++jl,++al)  {j=mj[jl];
        aij = pij(i,j);
        throwassert(aij);
	*aij += *al;}}
    break;
     
  case MatriceElementaire<R>::Symmetric : ffassert(symetrique);   
    for (il=0; il<me.n; ++il) {  i=mi[il] ;
      for (jl=0;jl< il+1 ; ++jl) { j=mj[jl];
	 aij =    (j<i) ? pij(i,j) : pij(j,i);
         throwassert(aij);
         *aij += *al++;}}
    break;
  default:
    cerr << "Big bug type MatriceElementaire unknown" << (int) me.mtype << endl;
    exit(1);
    break; 
  }      
  return *this;
} 

template<class R>
  void MatriceMorse<R>::Solve(KN_<R> &x,const KN_<R> &b) const{
    if (solver)    
      solver->Solver(*this,x,b);
    else
  {  cerr << "No Solver defined  for this Morse matrix " << endl;
    throw(ErrorExec("exit",1));}
  }


template<class R>
double MatriceMorse<R>::psor(KN_<R> & x,const  KN_<R> & gmin,const  KN_<R> & gmax , double omega) 
{
  double err=0;
  int n=this->n;
  ffassert(n==this->m);
  ffassert(n==x.N());
  ffassert(n==gmin.N());
  ffassert(n==gmax.N());
  if (symetrique)
   {
     ErrorExec("Error:sorry psor just for no symmetric Morse matrices",1);
   }
  else
   {
     for (int i=0;i<this->n;i++)
      {
       R xnew =x[i];
       R aii=R();
       for (int k=lg[i];k<lg[i+1];k++)
         {
           int j=cl[k];
           if(j!= i) 
             xnew -= a[k]*x[j];
            else aii=a[k];
         }
        if(aii != R())
           xnew /= aii;
         else ErrorExec("Error: psor diagonal coef = 0 ",1);
        R dx  = (xnew - x[i])*omega ;
        R xi = RNM::Min(RNM::Max(x[i]+dx,gmin[i]),gmax[i]);
        dx = x[i]- xi;
        err = Max(err, norm(dx));
        x[i] = xi;
        }
   }  return sqrt(err);
  
}

template<class R>
double MatriceProfile<R>::psor(KN_<R> & x,const  KN_<R> & gmin,const  KN_<R> & gmax , double omega) 
{
  double rr=0;
  ErrorExec("Error:sorry psor just for no symmetric Morse matrices (will do in future FH??? )",2);
  return rr;
  
}

template<class R>
void MatriceProfile<R>::setdiag(const KN_<R> & x) 
{
  ffassert(D);
 ffassert( this->n == x.N());
  KN_<R> d(D,this->n) ;
  d=x;
}
template<class R>
void MatriceProfile<R>::getdiag(KN_<R> & x) const 
{
  ffassert(D);
  ffassert( this->n == x.N());
  KN_<R> d(D,this->n) ;
  x=d;  
}
template<class R>
void MatriceMorse<R>::setdiag(const KN_<R> & x) 
{
 ffassert( this->n == this->m&& this->n == x.N());
 for (int i=0;i<this->n;++i)
    {
      R * p= pij(i,i);
      if(p)     *p = x[i];
      else ffassert( norm(x[i]) < 1e-30);}
}
template<class R>
void MatriceMorse<R>::getdiag(KN_<R> & x) const 
{
 ffassert( this->n == this->m && this->n == x.N());
 for (int i=0;i<this->n;++i)
    {
      R * p= pij(i,i);
      x[i]=  p ?  *p : R() ;
    }
  
}
template<class R>
R MatriceMorse<R>::pscal(const KN_<R> & x,const KN_<R> & y)
{ // (x, Ay)
  R sum=R();
  int i,j,k;
  ffassert(this->n==x.N());
  ffassert(this->m==y.N());  
  if (symetrique)
   {
     for (i=0;i<this->n;i++)
       for (k=lg[i];k<lg[i+1];k++)
         {
           j=cl[k];
           sum += a[k]*x[i]*y[j];
           if (i!=j)
             sum += a[k]*x[j]*y[i];
         }
           
   }
  else
   {
     for (i=0;i<this->n;i++)
       for (k=lg[i];k<lg[i+1];k++)
         {
           j=cl[k];
           sum += a[k]*x[i]*y[j];
         }
   }
  return sum;
}
template<class R>
R MatriceProfile<R>::pscal(const KN_<R> & x,const KN_<R> & y)
{
 if (y.n != this->n || x.n != this->n ) ERREUR(MatriceProfile pscal(xa,x) ," longueur incompatible c (out)") ;
 int i,j,k,kf;
 R sum = R();
 ffassert(this->n == this->m);
 if (D) 
   for (i=0;i<this->n;i++) 
     sum += D[i]*x[i]*y[i];
 else
   for (i=0;i<this->n;i++) // no dia => identyty dai
     sum +=x[i]*y[i];
      
 if (L && pL )    
   for (kf=pL[0],i=0;  i<this->n;   i++  )  
     for ( k=kf,kf=pL[i+1], j=i-kf+k;   k<kf; j++,  k++  )
       sum += L[k]*x[i]*y[j],throwassert(i>=0 && i <this->n && j >=0 && j < this->m && k>=0 && k < pL[this->n]);
       
 if (U && pU)     
   for (kf=pU[0],j=0;  j<this->m;  j++)  
     for (k=kf,kf=pU[j+1], i=j-kf+k;   k<kf; i++,  k++  )
       sum += U[k]*x[i]*y[j],throwassert(i>=0 && i <this->n && j >=0 && j < this->m &&  k>=0 && k < pU[this->n]);
 
 return sum;
}

template<class R>
void MatriceMorse<R>::getcoef(KN_<R> & x) const 
{
 ffassert(x.N()==this->nbcoef);
 x = KN_<R>(this->a,nbcoef);  
}
template<class R>
void MatriceMorse<R>::setcoef(const KN_<R> & x)  
{
 ffassert(x.N()==nbcoef);
  KN_<R>(this->a,nbcoef) = x;
}
template<class R>
int MatriceMorse<R>::NbCoef() const  
{
  return this->nbcoef;
}

template<class R>
void MatriceProfile<R>::getcoef(KN_<R> & x) const 
{
 ffassert(x.N()==this->NbCoef());
 int k=0,kk;
 if (D)
  {  kk=this->n;
     x(SubArray(kk,k))  = KN_<R>(D,kk);
     k += kk; }
 if (L)
  {  kk= pL[this->n];
     x(SubArray(kk,k))  = KN_<R>(L,kk);
     k += kk; }
  if (U && (U != L)) 
  {  kk=  pU[this->n];
     x(SubArray(kk,k))  = KN_<R>(U,kk);
     k += kk; }
   
}
template<class R>
void MatriceProfile<R>::setcoef(const KN_<R> & x)  
{
 ffassert(x.N()==this->NbCoef());
   int k=0,kk;
 if (D)
  {  kk=this->n;
     KN_<R>(D,kk)=x(SubArray(kk,k))   ;
     k += kk; }
 if (L)
  {  kk= pL[this->n];
     KN_<R>(L,kk)=x(SubArray(kk,k))   ;
     k += kk; }
  if (U && (U != L)) 
  {  kk=  pU[this->n];
     KN_<R>(U,kk)=x(SubArray(kk,k)) ;
     k += kk; }

}
template<class R>
int MatriceProfile<R>::NbCoef() const  
{
  int s=0;
  if (D) s += this->n;
  if (L) s += pL[this->n];
  if (U && (U != L)) s += pU[this->n];
  return s;
}
#endif

