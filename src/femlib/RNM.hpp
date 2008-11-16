// ********** DO NOT REMOVE THIS BANNER **********
// ORIG-DATE:    29 fev 2000  
// -*- Mode : c++ -*-
//
// SUMMARY  : array modelisation 
// USAGE    : LGPL      
// ORG      : LJLL Universite Pierre et Marie Curie, Paris,  FRANCE 
// AUTHOR   : Frederic Hecht
// E-MAIL   : frederic.hecht@ann.jussieu.fr
//

/*
 

 
 Freefem++ is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation; either version 2.1 of the License, or
 (at your option) any later version.
 
 Freefem++  is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.
 
 You should have received a copy of the GNU Lesser General Public License
 along with Freefem++; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 

 */
#ifndef KNM_H_
#define KNM_H_
// version sept 2008 FH.
// ----------------------
// une tentative qui ne marche pas 
// de tableau constant 
#include <complex>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cassert>


using namespace std;
#define const_R   R

#include <cstdlib>
inline void Check_Kn(const char * str,const char * file,int line)
{
 cerr << "CHECK_KN: " << str << " in file: " << file << ", line " << line <<endl;
 assert(0);
 abort();
}

#define K_bigassert(i)  if (!(i)) Check_Kn(#i,__FILE__,__LINE__);

#ifdef CHECK_KN

#define K_throwassert(i)  if (!(i)) Check_Kn(#i,__FILE__,__LINE__);

#else
#define K_throwassert(i) 
#endif
// version du 29 fev 2000  
//  correction   for (... lj++,ui++  qui apelle le produit scalaire
//  petite correction  throwassert 
// ajoute de operateur /= et *= sur des vecteurs
//   suppression de constructeur qui pose de probleme   
//   correction oper +=  ...  dans RNM_op.h ligne 56  change = en oper
// version de 25 nov 99  sans const R
// ajoute de '.' pour extraire une colonne, ou ligne  , ...
//  version du 22 nov 1999   cast to KN_<const R> 
//   version du 21 nov 1999  correction delete 
//  version du 13 nov 1999
//  version du 18 mars 99
//  F. Hecht 
// attention les indexations les indexations peuvent changer
//  puisque que l'on peut prendre la transposer d'une matrice
// tableau 
// mais ils partent de 0 
// version corrigee du 15/11/98
// version avec sous tableau  ---  mars 99 
// -------
//  remarque du 8 mars 99 FH
// class pour prendre des sous-tableau
// attention aux PB de continute dans les tableaux
// on a supposer que les tableaux multi indices pouvait est vue comme 
// un tableau continue ce qui est generalement faux quand l'on en 
// prend un sous tableau
//   exemple: un tableau 3,5 est numerote comme:
//    0  3  6  9 12
//    1  4  7 10 13
//    2  5  8 11 14
//             step
//   indexi  n 1
//   indexj  m n
//   est le sous tableau  3,3  n'est pas numeroter consecutivement 
//  
//    Donc la fonction  IsVector1() nous dit si un tableau 
//    a un 2 ou 3 indices est ou non consecutif en memoire
//    
//  --  ajoute d'une classe VirtualMatrice
// pour modeliser le produit matrice vecteur 
//  x = A*v; via des fonctions virtuelle 
//  ---------------------------------- 
//   version du 6 mars 2001 FH
//   ---  initialisation --
//   --------------------------------
//   version du 9 oct 2001 FH
//   ajoute de constructeur par defaut d'une vecteur
//   +  set , pour definir le vecteur 
//   ou l'affectation (bof bof) 
// ---------------------
//  version sep 2002
//  ajoute  operateur >> pour  KN<R> et KN_<R>
//  --------------------  
//  version  april 2003
//  ajoute un gestion auto de 
//  la fonction InternalError pour les matriceVirtuel
//  --------------------  
//   version jan 2004
//   correction pour go ++ 
//   des operateur  #=  pour les matrices et tenseurs
//  ----------------------
//   version  feb 2004 
//   v(i(:)) =  w   //  i(1:10) 
//   w=u(i(:))  //  
//   version mars 2004 make small correction
//   in  ITAB operator problem if non type R a defi 
//   -------------------
//   Modif pour version avec les Complex   mai 2004                                                                                                 
//   (u,v)  donne le produit complex utiliser dans le produit matrice vecteur
//   (u,conj(v))  donne le produit hermitiene pour le gradient conjugue
//
//   -- de fonction dans le cas real                                                                                                                 
// modif for g++ 4.0 and xlc++   mai 2005
//  adding some this-> 
//   mars 2007
// correction in operator operation:b -1*c 
// aout 2007, 
//  correct y = A*x ; when y is unset 
//  correct y += A*x ; when y is unset 
//  re-correct += sep 2007
//  add size of the matrix in VirtualMatrix class.
// ----------------
inline double  conj(const double & x){return x;}
inline float  conj(const float &x){return x;}
inline long  conj(const long &x){return x;}
inline double  real(const double &x){return x;}
inline float  real(const float &x){return x;}

namespace RNM {
template<class T> inline T Min (const T &a,const T &b){return a < b ? a : b;}
template<class T> inline T Max (const T &a,const T & b){return a > b ? a : b;}
template<class T> inline T Abs (const T &a){return a <0 ? -a : a;}

template<class T> inline void Exchange (T& a,T& b) {T c=a;a=b;b=c;}
template<class T> inline T Max (const T &a,const T & b,const T & c){return Max(Max(a,b),c);}
template<class T> inline T Min (const T &a,const T & b,const T & c){return Min(Min(a,b),c);}
// specialisation cas complex ---
template<class T> 
inline complex<T> Min(const complex<T> &a,complex<T> &b)
{ return complex<T>(min(a.real(),b.real()),min(a.imag(),b.imag()));}
template<class T> 
inline complex<T> Max(const complex<T> &a,const complex<T> &b)
{ return complex<T>(max(a.real(),b.real()),max(a.imag(),b.imag()));}

/*inline complex<double> Min(const complex<double> &a,complex<double> &b)
{ return complex<double>(Min(real(a),real(b)),Min(imag(a),imag(b)));}
inline complex<double> Max(const complex<double> &a,const complex<double> &b)
{ return complex<double>(Max(real(a),real(b)),Max(imag(a),imag(b)));}
*/ }
//  ----                                                                                                                                             

template<class R> class KNMK_ ;
template<class R> class KNM_ ;
template<class R> class KN_ ;
template<class R> class TKN_ ; // KN_ transpose 
template<class R> class notKN_ ; // KN_ not 
template<class R> class notnotKN_ ; // KN_ not not 

template<class R> class KNMK ;
template<class R> class KNM ;
template<class R> class KN ;

template<class R> class conj_KN_ ;
template<class R> class Add_KN_;
template<class R> class Sub_KN_;
template<class R> class Mulc_KN_;
template<class R> class Add_Mulc_KN_;
template<class R> class Mul_KNM_KN_; 
template<class R> class DotStar_KN_;
template<class R> class DotSlash_KN_;

template<class R> class outProduct_KN_;
template<class R> class if_KN_;
template<class R> class if_arth_KN_;
template<class R> class ifnot_KN_;
template<class R,class I> class KN_ITAB; 

template<class R,typename A,typename B> class F_KN_;


#ifndef ffassert
#define ffassert assert
#endif

// gestion des erreur interne --
#ifndef InternalError
typedef void (* TypeofInternalErrorRoutine)(const char *) ;
static TypeofInternalErrorRoutine &InternalErrorRoutinePtr()
{ 
  static TypeofInternalErrorRoutine routine=0;
  return routine;
}

static void InternalError(const char * str) {
  if (InternalErrorRoutinePtr() ) (*InternalErrorRoutinePtr())(str);
  cerr << str; 
  exit(1);
}
inline void  SetInternalErrorRoutine(TypeofInternalErrorRoutine f) 
{
  InternalErrorRoutinePtr()=f;
}
#endif
//  -- 
template<class P,class Q> 
   struct PplusQ { const P & p;const Q & q;
    PplusQ(const P & pp,const Q & qq) : p(pp),q(qq){}
    };

template<class R> 
struct  VirtualMatrice { public:
    int N,M;
    VirtualMatrice(int nn,int mm): N(nn),M(mm) {}
    VirtualMatrice(int nn): N(nn),M(nn) {}
  //  y += A x
  virtual void addMatMul(const KN_<R> &  x, KN_<R> & y) const =0; 
  virtual void addMatTransMul(const KN_<R> &  x, KN_<R> & y) const 
    { InternalError("VirtualMatrice::addMatTransMul not implemented "); }
  virtual bool WithSolver() const {return false;} // by default no solver          
  virtual void Solve( KN_<R> &  x,const KN_<R> & b) const 
    { InternalError("VirtualMatrice::solve not implemented "); } 

#ifdef VersionFreeFempp
  virtual bool ChecknbLine  (int n) const= 0; 
  virtual bool ChecknbColumn  (int m) const =0; 
#else
  virtual bool ChecknbLine  (int n) const {return true;} 
  virtual bool ChecknbColumn  (int m) const {return true;}
#endif
  struct  plusAx { const VirtualMatrice * A; const KN_<R>   x;
   plusAx( const VirtualMatrice * B,const KN_<R> &  y) :A(B),x(y) 
      { ffassert(B->ChecknbColumn(y.N())); }
    };
    
   plusAx operator*(const KN_<R> &  x) const {return plusAx(this,x);}
   
  struct  plusAtx { const VirtualMatrice * A; const KN_<R>   x;
   plusAtx( const VirtualMatrice * B,const KN_<R> &  y) :A(B),x(y) 
    {ffassert(B->ChecknbLine(y.N()));} };
    
  struct  solveAxeqb { const VirtualMatrice * A; const KN_<R>   b;
   solveAxeqb( const VirtualMatrice * B,const KN_<R> &  y) :A(B),b(y) 
    {ffassert(B->ChecknbColumn(y.N()));} };
  
  virtual ~VirtualMatrice(){} 
};

    

//template <class R> class MatriceCreuseMulKN_;
//template <class R> class MatriceCreuseDivKN_;

class ShapeOfArray;

class FromTo{ public:
  long from,to;
  FromTo(long i,long j):from(i),to(j) {K_throwassert(i<j);} 
 };
 
class SubArray{ public:
  const long n,step,start;
//  SubArray(char  nn): n(-1),step(1),start(0) {}  
 explicit SubArray(long nn,long sta=0,long s=1): n(nn),step(s),start(sta) {}
  SubArray(const FromTo& ft) : n(ft.to-ft.from+1),step(1),start(ft.from) {}
  SubArray(const ShapeOfArray & ); // all 
  long end()  const  { return start+ step*n;}
  long last() const  { return start+ step*(n-1);}
  long len1() const  { return step*(n-1);}  
};


class ShapeOfArray{ protected:
  public:

    long n;     //   n  nb of item
    long step;  //   step  nb of between 2 item
    long next;  //  the   next array of same type in matrix for subarray  
              // by default  no next 
  ShapeOfArray(const ShapeOfArray & s,long nn): n(s.n),step(s.n),next(nn) {}              
  ShapeOfArray(long nn): n(nn),step(1),next(-1) {}
  
  ShapeOfArray(long nn,long s): n(nn),step(s),next(-1) {}
  
  ShapeOfArray(long nn,long s,long nextt): n(nn),step(s),next(nextt) {}
  
  ShapeOfArray(const ShapeOfArray &old,const SubArray &sub) 
         : n(sub.n),step(old.step*sub.step),next(old.next)
          { K_throwassert((sub.last())*old.step <= old.last());} // a constructor 
          
  ShapeOfArray(const ShapeOfArray &old,long stepo,long start) 
         : n(old.n-start),step(old.step*stepo),next(old.next) 
          { K_throwassert(n>=0);}        
          
  long end()  const      { return n*step;}
  long last()      const      { return (n-1)*step;}
  long constant()  const { return step==0;}
  long index(long k) const { K_throwassert( (k>=0) && ( (k <n) || !step) );
           return step*k;}
  ShapeOfArray operator*(long stepp) const {return  ShapeOfArray(n,step*stepp,next);}  
  bool SameShape(const ShapeOfArray & a) const 
          { return  !step || !a.step || a.n == n ;} 
  long  N(const ShapeOfArray & a) { return step ? n : a.n;} // size of 2 shape 

  
// protected:
  long operator[](long k) const {  K_throwassert( (k>=0) && ( (k <n) || !step) );
           return step*k;}     
  void init(long nn,long s=1,long nextt=-1) { n=nn; step=s; next=nextt;}         
};

ostream & operator<<(ostream & f,const ShapeOfArray & s);

inline bool  SameShape(const ShapeOfArray & a,const ShapeOfArray & b) 
           { return  !a.step || !b.step || a.n == b.n ;} 
   
inline long N(const ShapeOfArray & a,const ShapeOfArray & b) 
           { K_throwassert(SameShape(a,b)); return  a.step ? a.n :  b.n ;}

inline SubArray::SubArray(const ShapeOfArray & s) 
           :n(s.n),step(s.step),start(0) {}
           
   

template<class R> 
ostream & operator<<(ostream & f,const KN_<const_R> & v) ;

template<class R> istream & operator>>(istream & f, KN_<R> & v);
template<class R> istream & operator>>(istream & f, KN<R> & v);


template<class R>
class KN_: public  ShapeOfArray {
protected:
  R *v;
public:
  typedef R K; // type of data 
  long N() const {return n;}
  bool unset() const { return !v;}
  void set(R * vv,int nn,int st=1,int nx=-1) {v=vv;n=nn;step=st;next=nx;}
  long size() const{return step?n*step:n;}
  operator R *() const {return v;}
  KN_(const KN_<R> & u) :ShapeOfArray(u),v(u.v){} 
  KN_(const KN_<R> & U,const SubArray & sa)  : ShapeOfArray(U,sa),v(U.v + U.index(sa.start)) {}

  KN_ operator()(const SubArray & sa) const { return KN_(*this,sa);} // sub array 
  
  R & operator[](long i) const {return v[index(i)];}
  R & operator()(long i) const {return v[index(i)];}
  R & operator[](int i) const {return v[index(i)];}
  R & operator()(int i) const {return v[index(i)];}

  R operator,(const KN_<const_R> & v) const; // dot  product 

   KN_& operator  =(const KN_<const_R> & u)  ;
   KN_& operator +=(const KN_<const_R> & u)  ;
   KN_& operator -=(const KN_<const_R> & u)  ;
  
   KN_& operator *=(const KN_<const_R> & u)  ;
   KN_& operator /=(const KN_<const_R> & u)  ;
  
  
   KN_& operator = (const_R  a) ;  
   KN_& operator +=(const_R  a) ;
   KN_& operator -=(const_R  a) ;  
   KN_& operator /=(const_R  a) ;
   KN_& operator *=(const_R  a) ;
  
   KN_& operator  = (R*  a) { return operator =(KN_<R>(a,n));}
   KN_& operator += (R*  a) { return operator+=(KN_<R>(a,n));}  
   KN_& operator -= (R*  a) { return operator-=(KN_<R>(a,n));}  
   KN_& operator *= (R*  a) { return operator*=(KN_<R>(a,n));}  
   KN_& operator /= (R*  a) { return operator/=(KN_<R>(a,n));}  


  
  R min() const ;
  R max() const ;
  R sum() const ;
  double norm() const ;
  double l2() const ;
  double l1() const ;
  double linfty() const ;
  double lp(double p) const ;
  
  template<class T> long last(const T &) const;
  template<class T> long first(const T &) const;
  
  void map(R (*f)(R )); // apply the f fonction a all element of the array
  void map(R (*f)(const  R& )); // apply the f fonction a all element of the array

 template<class T>
   void set(R (*f)(const  T& ),KN_<T> & u); // apply the f fonction a all element of the array u
  
   KN_& operator =(const DotStar_KN_<R> & u) ;
   KN_& operator+=(const DotStar_KN_<R> & u) ;
   KN_& operator-=(const DotStar_KN_<R> & u) ;
   KN_& operator*=(const DotStar_KN_<R> & u) ;
   KN_& operator/=(const DotStar_KN_<R> & u) ;

   KN_& operator =(const DotSlash_KN_<R> & u) ;
   KN_& operator+=(const DotSlash_KN_<R> & u) ;
   KN_& operator-=(const DotSlash_KN_<R> & u) ;
   KN_& operator*=(const DotSlash_KN_<R> & u) ;
   KN_& operator/=(const DotSlash_KN_<R> & u) ;

   KN_& operator =(const if_KN_<R> & u) ;
   KN_& operator+=(const if_KN_<R> & u) ;
   KN_& operator-=(const if_KN_<R> & u) ;
   KN_& operator*=(const if_KN_<R> & u) ;
   KN_& operator/=(const if_KN_<R> & u) ;

   KN_& operator =(const ifnot_KN_<R> & u) ;
   KN_& operator+=(const ifnot_KN_<R> & u) ;
   KN_& operator-=(const ifnot_KN_<R> & u) ;
   KN_& operator*=(const ifnot_KN_<R> & u) ;
   KN_& operator/=(const ifnot_KN_<R> & u) ;

   KN_& operator =(const Add_KN_<R> & u) ;
   KN_& operator+=(const Add_KN_<R> & u) ;
   KN_& operator-=(const Add_KN_<R> & u) ;
   KN_& operator*=(const Add_KN_<R> & u) ;
   KN_& operator/=(const Add_KN_<R> & u) ;
  
   template<class I,class T> KN_& operator =  (const KN_ITAB<T,I> & u);
   template<class I,class T> KN_& operator +=  (const KN_ITAB<T,I> & u);
   template<class I,class T> KN_& operator -=  (const KN_ITAB<T,I> & u);
   template<class I,class T> KN_& operator *=  (const KN_ITAB<T,I> & u);
   template<class I,class T> KN_& operator /=  (const KN_ITAB<T,I> & u);


   KN_ITAB< KN_<R>,const KN_<int> >  operator()(const KN_<int> &itab) ;
   KN_ITAB< KN_<R>,const KN_<long> >  operator()(const KN_<long> &itab) ;
   KN_ITAB<const KN_<R>,const KN_<int> > operator()(const KN_<int> &itab) const ;
   KN_ITAB<const KN_<R>,const KN_<long> >  operator()(const KN_<long> &itab) const ;


    

  
   KN_& operator =(const Sub_KN_<R> & u) ;
   KN_& operator-=(const Sub_KN_<R> & u) ;
   KN_& operator+=(const Sub_KN_<R> & u) ;
   KN_& operator*=(const Sub_KN_<R> & u) ;
   KN_& operator/=(const Sub_KN_<R> & u) ;
  
   KN_& operator =(const Mulc_KN_<R> & u) ;
   KN_& operator+=(const Mulc_KN_<R> & u) ;
   KN_& operator-=(const Mulc_KN_<R> & u) ;
   KN_& operator*=(const Mulc_KN_<R> & u) ;
   KN_& operator/=(const Mulc_KN_<R> & u) ;
  
   KN_& operator =(const Add_Mulc_KN_<R> & u) ;
   KN_& operator+=(const Add_Mulc_KN_<R> & u) ;
   KN_& operator-=(const Add_Mulc_KN_<R> & u) ;
   KN_& operator*=(const Add_Mulc_KN_<R> & u) ;
   KN_& operator/=(const Add_Mulc_KN_<R> & u) ;

   KN_& operator =(const if_arth_KN_<R> & u) ;
   KN_& operator+=(const if_arth_KN_<R> & u) ;
   KN_& operator-=(const if_arth_KN_<R> & u) ;
   KN_& operator*=(const if_arth_KN_<R> & u) ;
   KN_& operator/=(const if_arth_KN_<R> & u) ;

  
   KN_& operator =(const Mul_KNM_KN_<R> & u) ; 
   KN_& operator+=(const Mul_KNM_KN_<R> & u) ; 
   KN_& operator-=(const Mul_KNM_KN_<R> & u) ; 
   KN_& operator*=(const Mul_KNM_KN_<R> & u) ; 
   KN_& operator/=(const Mul_KNM_KN_<R> & u) ; 
  
 //  KN_& operator =(const MatriceCreuseMulKN_<R> & ) ;
 //  KN_& operator +=(const MatriceCreuseMulKN_<R> & ) ;
   KN_& operator =(const typename VirtualMatrice<R>::plusAx & Ax)  
    {*this=R(); Ax.A->addMatMul(Ax.x,*this);return *this;}
   KN_& operator =(const typename VirtualMatrice<R>::plusAtx & Ax)  
    {*this=R(); Ax.A->addMatTransMul(Ax.x,*this);return *this;}
   KN_& operator +=(const typename VirtualMatrice<R>::plusAx & Ax)  
    {  Ax.A->addMatMul(Ax.x,*this);return *this;}
   KN_& operator +=(const typename VirtualMatrice<R>::plusAtx & Ax)  
    {  Ax.A->addMatTransMul(Ax.x,*this);return *this;}
   KN_& operator =(const typename VirtualMatrice<R>::solveAxeqb & Ab)  
    {*this=R(); Ab.A->Solve(*this,Ab.b);return *this;}
    
  template<class  A,class B,class C> KN_&  operator =  (const F_KN_<A,B,C>  & u) ;
  template<class  A,class B,class C> KN_&  operator +=  (const F_KN_<A,B,C>  & u) ;
  template<class  A,class B,class C> KN_&  operator -=  (const F_KN_<A,B,C>  & u) ;
  template<class  A,class B,class C> KN_&  operator /=  (const F_KN_<A,B,C>  & u) ;
  template<class  A,class B,class C> KN_&  operator *=  (const F_KN_<A,B,C>  & u) ;
    
   
//   KN_& operator =(const MatriceCreuseDivKN_<R> &)  ;

 friend   ostream & operator<< <R>(ostream & f,const KN_<const_R> & v)  ;
 
  KN_(R *u,const ShapeOfArray & s):ShapeOfArray(s),v(u){}
  KN_(R *u,long nn,long s):ShapeOfArray(nn,s),v(u){}
  KN_(R *u,long nn,long s,long nextt):ShapeOfArray(nn,s,nextt),v(u){}
  KN_(R *u,long nn):ShapeOfArray(nn),v(u){}


  TKN_<R>  t() ; // transpose
  const TKN_<R>  t() const ; // transpose
  notKN_<R>  operator!()  ; //  not
  const  notKN_<R>  operator!() const ; // not

  //  operator KN<R> &();
  // operator const KN<R> &() const;

 private:
 
  KN_&  operator++(){K_throwassert(next>=0);v += next;return *this;} //    ++U
  KN_&  operator--(){K_throwassert(next>=0);v -= next;return *this;} //    --U
  KN_   operator++(int ){K_throwassert(next>=0); KN_ old=*this;v = v +next;return old;} // U++
  KN_   operator--(int ){K_throwassert(next>=0); KN_ old=*this;v = v -next;return old;} // U++
 
  KN_(const KN_<R> & u,long offset) :ShapeOfArray(u),v(&u[offset]){} 
  KN_(const KN_<R> & u,const ShapeOfArray &sh,long startv=0)
         :ShapeOfArray(sh*u.step),v(&u[startv]){}
  KN_(const KN_<R> & u,long nnext,const ShapeOfArray &sh,long startv=0)
         :ShapeOfArray(sh.n,sh.step*u.step,nnext),v(&u[startv]){ }

  
// friend class KN_<R>;   
 friend class KNM_<R>;   
 friend class KNMK_<R>;   
 friend class KN<R>;   
 friend class KNM<R>;   
 friend class KNMK<R>;   
  
};

template<class R>
class KNM_: public KN_<R> {
  public:
  ShapeOfArray shapei;
  ShapeOfArray shapej;
  public:
  long IsVector1() const  {  return (shapei.n*shapej.n) ==  this->n ;} 
  long N() const {return shapei.n;}
  long M() const {return shapej.n;}  
  long size() const { return shapei.n*shapej.n;}
  
  KNM_(R* u,const ShapeOfArray & s,
            const ShapeOfArray & si,
            const ShapeOfArray & sj)
             : KN_<R>(u,s),shapei(si),shapej(sj){} 
  KNM_(R* u,long n,long m)
             : KN_<R>(u,ShapeOfArray(n*m)),shapei(n,1,n),shapej(m,n,1){}
  KNM_(R* u,long n,long m,long s)
             : KN_<R>(u,ShapeOfArray(n*m,s)),shapei(n,1,n),shapej(m,n,1){}                     
  KNM_(KN_<R> u,long n,long m) 
             : KN_<R>(u,ShapeOfArray(m*n)),shapei(n,1,n),shapej(m,n,1){ }
             
  KNM_(const KN_<R> &u,const ShapeOfArray & si,const ShapeOfArray & sj,long offset=0) 
             : KN_<R>(&u[offset],si.last()+sj.last()+1,u.step),shapei(si),shapej(sj) 
             {K_throwassert( offset>=0 && this->n+ (this->v-(R*)u) <= u.n);}
  KNM_(const KN_<R> &u,const ShapeOfArray & si,const ShapeOfArray & sj,long offset,long nnext) 
             : KN_<R>(&u[offset],si.last()+sj.last()+1,u.step,nnext),shapei(si),shapej(sj) 
             {K_throwassert( offset>=0 && this->n+ (this->v-(R*)u) <= u.n);}
  
  KNM_(KNM_<R> U,const SubArray & si,const SubArray & sj)  
             :KN_<R>(U,SubArray(U.ij(si.len1(),sj.len1())+1,U.ij(si.start,sj.start))),
                     shapei(U.shapei,si),shapej(U.shapej,sj){} 

  KNM_(KNM_<R> U,const SubArray & sa,const SubArray & si,const SubArray & sj)  
             :KN_<R>(U,SubArray(sa)),shapei(U.shapei,si),shapej(U.shapej,sj){} 

  KNM_(const KNM_<R> & u) 
             :KN_<R>(u),shapei(u.shapei),shapej(u.shapej) {}

  KNM_ operator()(const SubArray & sa,const SubArray & sb) const 
            { return KNM_(*this,sa,sb);} // sub array 
  
  long ij(long i,long j) const   
            { return shapei.index(i)+shapej.index(j);}
  long indexij(long i,long j)        const   
            { return this->index(shapei.index(i)+shapej.index(j));}
  R & operator()(long i,long j)     const   
            { return this->v[indexij(i,j)];}
  R & operator()(int i,int j)     const   
            { return this->v[indexij(i,j)];}
            
            
  KN_<R> operator()(const char,long j    )  const   // une colonne j  ('.',j)
            { return KN_<R>(&this->v[this->index(shapej.index(j))],shapei*this->step);} 
  KN_<R> operator()(long i    ,const char)  const   // une ligne i  (i,'.')
            { return KN_<R>(&this->v[this->index(shapei.index(i))],shapej*this->step);}  
  KN_<R> operator()(const char,int j    )  const   // une colonne j  ('.',j)
            { return KN_<R>(&this->v[this->index(shapej.index(j))],shapei*this->step);} 
  KN_<R> operator()(int i    ,const char)  const   // une ligne i  (i,'.')
            { return KN_<R>(&this->v[this->index(shapei.index(i))],shapej*this->step);}  
  KN_<R> operator()(const char,const char)  const   // tous 
            { return *this;}
  KNM_<R> t() const
    { return KNM_<R>(this->v,*this,shapej,shapei);} //  before  { return KNM_<R>(*this,shapej,shapei,v);}
  
   KNM_& operator =(const KNM_<const_R> & u) ;
   KNM_& operator =(const_R a)               ;
   KNM_& operator+=(const_R a)               ;
   KNM_& operator-=(const_R a)               ; 
   KNM_& operator/=(const_R a)               ;
   KNM_& operator*=(const_R a)               ; 
   KNM_& operator+=(const KNM_<const_R> & u) ;
   KNM_& operator-=(const KNM_<const_R> & u) ;
   KNM_& operator*=(const KNM_<const_R> & u) ;
   KNM_& operator/=(const KNM_<const_R> & u) ;
   
   KNM_ &operator =(const outProduct_KN_<R> &);
   KNM_ &operator +=(const outProduct_KN_<R> &);
   KNM_ &operator -=(const outProduct_KN_<R> &);
   KNM_ &operator /=(const outProduct_KN_<R> &); // bofbof
   KNM_ &operator *=(const outProduct_KN_<R> &); // bofbof

private:  
  KNM_& operator++() {this->v += this->next;return *this;} // ++U
  KNM_& operator--() {this->v -= this->next;return *this;} // ++U
  KNM_  operator++(int ){KNM_<R> old=*this;this->v = this->v +this->next;return old;} // U++
  KNM_  operator--(int ){KNM_<R> old=*this;this->v = this->v -this->next;return old;} // U--


 friend class KN_<R>;   
// friend class KNM_<R>;   
 friend class KNMK_<R>;   
 friend class KN<R>;   
 friend class KNM<R>;   
 friend class KNMK<R>;   
};

template<class T,class I> 
struct KN_ITAB
{
  KN_ITAB(const T &vv,const I &iindex) : v(vv),index(iindex) {}
  T  v;
  I  index;
  KN_ITAB & operator=(const T & t);
  KN_ITAB & operator+=(const T & t);
  KN_ITAB & operator-=(const T & t);
  KN_ITAB & operator*=(const T & t);
  KN_ITAB & operator/=(const T & t);
  typename T::R & operator[](long i){ return v[index[i]];}
  const typename T::R &  operator[](long i) const { return v[index[i]];}
  long N() const { return index.N();}    
};

template<class R>  KN_ITAB<const KN_<R>,const KN_<int> >  KN_<R>::operator()(const KN_<int>  &itab) const { return KN_ITAB<const KN_<R>,const KN_<int> > (*this,itab);}
template<class R>  KN_ITAB<const KN_<R>,const KN_<long> >  KN_<R>::operator()(const KN_<long>  &itab) const { return KN_ITAB<const KN_<R>,const KN_<long> > (*this,itab);}
template<class R>  KN_ITAB< KN_<R>,const KN_<int> >  KN_<R>::operator()(const KN_<int>  &itab) { return KN_ITAB<KN_<R>,const KN_<int> > (*this,itab);}
template<class R>  KN_ITAB< KN_<R>,const KN_<long> >  KN_<R>::operator()(const KN_<long>  &itab) { return KN_ITAB<KN_<R>,const KN_<long> > (*this,itab);}


template<class R>
struct TKN_:public KN_<R> {
    TKN_(const KN_<R> &x) : KN_<R>(x) {}
};

template<class R>
struct notKN_:public KN_<R> {
    notKN_(const KN_<R> &x) : KN_<R>(x) {}
    notnotKN_<R>  operator!()  ; //  not
    const  notnotKN_<R>  operator!() const ; // not
};

template<class R>
struct notnotKN_:public KN_<R> {
    notnotKN_(const notKN_<R> &x) : KN_<R>(x) {}
    notKN_<R>  operator!()  ; //  notnot
    const  notKN_<R>  operator!() const ; // notnot
};

template<class R>
TKN_<R>  KN_<R>::t() { return *this;} // transpose

template<class R>
const TKN_<R>  KN_<R>::t() const { return *this;} // transpose

template<class R>
notKN_<R>  KN_<R>::operator!() { return *this;} // not

template<class R>
const notKN_<R>  KN_<R>::operator!() const { return *this;} // not

template<class R>
notnotKN_<R>  notKN_<R>::operator!() { return *this;} // not

template<class R>
const notnotKN_<R>  notKN_<R>::operator!() const { return *this;} // not


template<class R>
struct outProduct_KN_ {
    const KN_<R>  a,b;
    R c;
    outProduct_KN_(const KN_<R> & aa, const KN_<R> &bb,R cc=(R)1) : a(aa),b(bb),c(cc) {}
    outProduct_KN_(const KN_<R> * aa, const KN_<R> &bb,R cc=(R)1) : a(*aa),b(bb),c(cc) {}
    outProduct_KN_(const KN_<R> * aa, const KN_<R> *bb,R cc=(R)1) : a(*aa),b(*bb),c(cc) {}
    outProduct_KN_(const Mulc_KN_<R> & aa,const KN_<R> & bb) : a(aa.a),b(bb),c(aa.b) {}    
    outProduct_KN_ operator * (R cc) { return outProduct_KN_(a,b,c*cc);}    
};

template<class R>
struct if_KN_ {
    const KN_<R> & a,&b;
    R c;
    if_KN_(const KN_<R> & aa, const KN_<R> &bb,R cc=1.) : a(aa),b(bb),c(cc) {}
    if_KN_ operator * (R cc) { return if_KN_(a,b,c*cc);}    
};

template<class R>
struct ifnot_KN_ {
    const KN_<R> & a,&b;
    R c;
    ifnot_KN_(const KN_<R> & aa, const KN_<R> &bb,R cc=1.) : a(aa),b(bb),c(cc) {}
    ifnot_KN_ operator * (R cc) { return ifnot_KN_(a,b,c*cc);}    
};


template<class R> 
outProduct_KN_<R> operator*(const KN_<R> &a,const TKN_<R> &b) 
{ return outProduct_KN_<R>(a,b);}

template<class R> 
ifnot_KN_<R> operator*(const KN_<R> &a,const notKN_<R> &b) 
{ return ifnot_KN_<R>(b,a);}

template<class R> 
ifnot_KN_<R> operator*(const KN_<R> &a,const notnotKN_<R> &b) 
{ return if_KN_<R>(b,a);}

template<class R> 
ifnot_KN_<R> operator*(const notKN_<R> &b,const KN_<R> &a) 
{ return ifnot_KN_<R>(b,a);}

template<class R> 
ifnot_KN_<R> operator*(const notnotKN_<R> &b,const KN_<R> & a) 
{ return if_KN_<R>(b,a);}


template<class R> 
R operator*(const TKN_<R> &a,const KN_<R> &b) 
{ return (a,b);}

template<class R>
class KNMK_: public KN_<R> {
  friend class KNMK<R>;
  public:
  ShapeOfArray shapei;
  ShapeOfArray shapej;
  ShapeOfArray shapek;
  public:
  long IsVector1() const {  return (shapei.n*shapej.n*shapek.n) == this->n ;} 
  long N() const {return shapei.n;}
  long M() const {return shapej.n;}
  long K() const {return shapek.n;}
  long size() const { return shapei.n*shapej.n*shapek.n;}
  KNMK_(const ShapeOfArray & s,
        const ShapeOfArray & si,
        const ShapeOfArray & sj,
        const ShapeOfArray & sk,
	    R * u)
    : KN_<R>(u,s),shapei(si),shapej(sj),shapek(sk){} 
    
  KNMK_(R* u,long n,long m,long k)
    : KN_<R>(u, ShapeOfArray(n*m*k)),shapei(n,1,n),shapej(m,n,1),shapek(k,n*m,n*m){};
    
//  KNMK_(const KN_<R> & u,long n,long m,long k)
//   : KN_<R>(ShapeOfArray(n*m*k)),shapei(n,1,n),shapekj(m,n,1),u),
//     shapek(k,n*m,n*m){};

  KNMK_(const KNMK_<R> &U,const SubArray & si,const SubArray & sj,const SubArray & sk)  :
    KN_<R>(U,SubArray(U.ijk(si.len1(),sj.len1(),sk.len1())+1,
                       U.ijk(si.start,sj.start,sk.start))),
                       shapei(U.shapei,si),
                       shapej(U.shapej,sj),
                       shapek(U.shapek,sk){} 

  KNMK_(const KNMK_<R> & u) :KN_<R>(u),shapei(u.shapei),shapej(u.shapej),shapek(u.shapek) {}

    
  long ijk(long i,long j,long k) const 
              { return shapei.index(i)+shapej.index(j)+shapek.index(k);}
  long indexijk(long i,long j,long k) const 
              {return this->index(shapei.index(i)+shapej.index(j)+shapek.index(k));} 
                           
  R & operator()(long i,long j,long k)   const   {return this->v[indexijk(i,j,k)];}
  R & operator()(int i,int j,int k)   const   {return this->v[indexijk(i,j,k)];}
  
//  pas de tableau suivant
 KN_<R>  operator()(const char ,long j,long k)  const  { // le tableau (.,j,k) 
        return KN_<R>(*this,-1,shapei,shapej[j]+shapek[k]);}
 KN_<R>  operator()(long i,const char ,long k)  const  { // le tableau (i,.,k) 
        return KN_<R>(*this,-1,shapej,shapei[i]+shapek[k]);}
 KN_<R>  operator()(long i,long j,const char )  const  { // le tableau (i,j,.) 
        return KN_<R>(*this,-1,shapek,shapei[i]+shapej[j]);}

 KN_<R>  operator()(const char ,int j,int k)  const  { // le tableau (.,j,k) 
        return KN_<R>(*this,-1,shapei,shapej[j]+shapek[k]);}
 KN_<R>  operator()(int i,const char ,int k)  const  { // le tableau (i,.,k) 
        return KN_<R>(*this,-1,shapej,shapei[i]+shapek[k]);}
 KN_<R>  operator()(int i,int j,const char )  const  { // le tableau (i,j,.) 
        return KN_<R>(*this,-1,shapek,shapei[i]+shapej[j]);}
//                                              
 KNM_<R>  operator()(const char ,const char ,long k)  const  { // le tableau (.,.,k) 
        return KNM_<R>(*this,shapei,shapej,shapek[k],shapek.next);} // step = n*m
 //attention les suivants ne marche pas
 KNM_<R>  operator()(const char ,long j,const char )  const  { // le tableau (.,j,.) 
        return KNM_<R>(*this,shapei,shapek,shapej[j],-1/*shapej.next*/);} // step = n
        
 KNM_<R>  operator()(long i,const char ,const char )  const  { // le tableau (i,.,.) 
        return KNM_<R>(*this,shapej,shapek,shapei[i],-1/*shapei.next*/);}  // step = 1

 KNM_<R>  operator()(const char ,const char ,int k)  const  { // le tableau (.,.,k) 
        return KNM_<R>(*this,shapei,shapej,shapek[k],shapek.next);} // step = n*m
 //attention les suivants ne marche pas
 KNM_<R>  operator()(const char ,int j,const char )  const  { // le tableau (.,j,.) 
        return KNM_<R>(*this,shapei,shapek,shapej[j],-1/*shapej.next*/);} // step = n
        
 KNM_<R>  operator()(int i,const char ,const char )  const  { // le tableau (i,.,.) 
        return KNM_<R>(*this,shapej,shapek,shapei[i],-1/*shapei.next*/);}  // step = 1

   KNMK_& operator =(const KNMK_<const_R> & u) ;
   KNMK_& operator+=(const KNMK_<const_R> & u)  ;
   KNMK_& operator-=(const KNMK_<const_R> & u)  ;
   KNMK_& operator/=(const KNMK_<const_R> & u)  ;
   KNMK_& operator*=(const KNMK_<const_R> & u)  ;
   KNMK_& operator =(const_R a)  ; 
   KNMK_& operator+=(const_R a)  ;
   KNMK_& operator-=(const_R a)  ;
   KNMK_& operator/=(const_R a)  ;
   KNMK_& operator*=(const_R a)  ;

  KNMK_  operator()(SubArray si,SubArray sj,SubArray sk) const 
        {return KNMK_(*this,si,sj,sk);}

  private:
//  KNMK_&  operator++(){v += next;return *this;} // ++U
//  KNMK_&  operator--(){v -= next;return *this;} // --U
//  KNMK_  operator++(long ){KNMK_ old=*this;v = v +next;return old;} // U++ 
//  KNMK_  operator--(long ){KNMK_ old=*this;v = v -next;return old;} // U--
 
        
friend class KNM_<R>;   
friend class KN_<R>;   

};



template<class R>
class KN :public KN_<R> { public:

  typedef R K;

 // explicit  KN(const R & u):KN_<R>(new R(uu),1,0) {}
  KN() : KN_<R>(0,0) {}
  KN(long nn) : KN_<R>(new R[nn],nn)         {} 
  KN(long nn, R * p) : KN_<R>(new R[nn],nn)  
    { KN_<R>::operator=(KN_<R>(p,nn));}
  KN(long nn,R (*f)(long i) ) : KN_<R>(new R[nn],nn) 
        {for(long i=0;i<this->n;i++) this->v[i]=f(i);}  
  KN(long nn,const  R & a) : KN_<R>(new R[nn],nn) 
        { KN_<R>::operator=(a);} 
  KN(long nn,long s,const  R  a) : KN_<R>(new R[nn],nn,s) 
        { KN_<R>::operator=(a);} 
  template<class S>   KN(const KN_<S> & s):KN_<R>(new R[s.n],s.n) 
        {for (long i=0;i<this->n;i++) this->v[i] = s[i];}
  template<class S>  KN(const KN_<S> & s,R (*f)(S )):KN_<R>(new R[s.n],s.n) 
        {for (long i=0;i<this->n;i++) this->v[i] = f(s[i]);}
  KN(const KN<R> & u):KN_<R>(new R[u.n],u.n)
        { KN_<R>::operator=(u);}
  KN(bool ,KN<R> & u):KN_<R>(u) {u.v=0;u.n=0;}// remove copy for return of local KN. 
    
  //  explicit KN(const KN_<R> & u):KN_<R>(new R[u.n],u.n) 
  //      { KN_<R>::operator=(u);}
        
  ~KN(){delete [] this->v;}
  void CheckSet() { if(!(this->n)) {cerr << "Error RNM set array\n";K_throwassert(0); exit(1);}}
   KN& operator  = (R*  a) { CheckSet(); return operator =(KN_<R>(a,this->n));}
   KN& operator += (R*  a) { CheckSet(); return operator+=(KN_<R>(a,this->n));}  
   KN& operator -= (R*  a) { CheckSet(); return operator-=(KN_<R>(a,this->n));}  
   KN& operator *= (R*  a) { CheckSet(); return operator*=(KN_<R>(a,this->n));}  
   KN& operator /= (R*  a) { CheckSet(); return operator/=(KN_<R>(a,this->n));}  
  
   KN& operator =(const_R a)  
        { if(this->unset()) set(new R[1],1,0,0); KN_<R>::operator= (a);return *this;}
   KN& operator =(const KN_<R>& a)  
        { if(this->unset()) set(new R[a.N()],a.N()); KN_<R>::operator= (a);return *this;}                
   KN& operator =(const KN<R>& a)  
        { if(this->unset()) set(new R[a.N()],a.N()); KN_<R>::operator= (a);return *this;}                
   KN& operator =(const Add_KN_<R> & u)  
        { if(this->unset()) set(new R[u.a.N()],u.a.N());KN_<R>::operator=(u);return *this;}
   KN& operator =(const DotStar_KN_<R> & u)  
        { if(this->unset()) set(new R[u.a.N()],u.a.N());KN_<R>::operator=(u);return *this;}
   KN& operator =(const if_KN_<R> & u)  
        { if(this->unset()) set(new R[u.a.N()],u.a.N());KN_<R>::operator=(u);return *this;}
   KN& operator =(const ifnot_KN_<R> & u)  
        { if(this->unset()) set(new R[u.a.N()],u.a.N());KN_<R>::operator=(u);return *this;}
   KN& operator =(const DotSlash_KN_<R> & u)  
        { if(this->unset()) set(new R[u.a.N()],u.a.N());KN_<R>::operator=(u);return *this;}
   KN& operator =(const Sub_KN_<R> & u)  
        { if(this->unset()) set(new R[u.a.N()],u.a.N());KN_<R>::operator=(u);return *this;}
   KN& operator =(const Mulc_KN_<R> & u)  
        { if(this->unset()) set(new R[u.a.N()],u.a.N());KN_<R>::operator=(u);return *this;}
   KN& operator =(const Add_Mulc_KN_<R> & u)  
        { if(this->unset()) set(new R[u.a.N()],u.a.N());KN_<R>::operator=(u);return *this;}
   KN& operator =(const if_arth_KN_<R> & u)  
        { if(this->unset()) set(new R[u.a.N()],u.a.N());KN_<R>::operator=(u);return *this;}
        
        
   KN& operator =(const Mul_KNM_KN_<R> & u) 
        { if(this->unset()) set(new R[u.b.N()],u.b.N());KN_<R>::operator=(u);return *this;}
//   KN& operator =(const MatriceCreuseMulKN_<R> & Ax) 
//       {if(this->unset()) set(new R[Ax.v.N()],Ax.v.N()); KN_<R>::operator=(Ax);return *this;}
//   KN& operator +=(const MatriceCreuseMulKN_<R> & Ax) 
//       {if(this->unset()) set(new R[Ax.v.N()],Ax.v.N()); KN_<R>::operator+=(Ax);return *this;}
//   KN& operator =(const MatriceCreuseDivKN_<R> & A1x)  
//       { if(this->unset()) set(new R[A1x.v.N()],A1x.v.N());KN_<R>::operator=(A1x);return *this;}
  // correcton aout 2007 FH  add N,M flied in VirtualMatrice
   KN& operator =(const typename VirtualMatrice<R>::plusAx & Ax)  
        { if(this->unset() && Ax.A->N ) set(new R[Ax.A->N],Ax.A->N);KN_<R>::operator=(Ax);return *this;}
   KN& operator =(const typename VirtualMatrice<R>::solveAxeqb & Ab)  
        { if(this->unset()) set(new R[Ab.b.N()],Ab.b.N());KN_<R>::operator=(Ab);return *this;}
   KN& operator +=(const typename  VirtualMatrice<R>::plusAx & Ax)  
  { if(this->unset()  && Ax.A->N) {
        set(new R[Ax.A->N],Ax.A->N);
        KN_<R>::operator=(R());}
    KN_<R>::operator+=(Ax);
    return *this;}
   KN& operator =(const typename VirtualMatrice<R>::plusAtx & Ax)  
        { if(this->unset()&&Ax.A->M) set(new R[Ax.A->M],Ax.A->M);KN_<R>::operator=(Ax);return *this;}
   KN& operator +=(const typename VirtualMatrice<R>::plusAtx & Ax)  
  { if(this->unset()&&Ax.A->M) {
       set(new R[Ax.A->M],Ax.A->M);
      KN_<R>::operator=(R());}
      KN_<R>::operator+=(Ax);
     return *this;}
// end correcton FH
   template<class P,class Q> 
     KN& operator =(const  PplusQ<P,Q> & PQ)  
      { *this=PQ.p; *this+=PQ.q;return *this; } 
   template<class P,class Q> 
     KN& operator +=(const  PplusQ<P,Q> & PQ)  
      { *this+=PQ.p; *this+=PQ.q;return *this; } 
           
   KN& operator -=(const_R a)  
        { KN_<R>::operator-=(a);return *this;}
   KN& operator -=(const KN_<R>& a)  
        { KN_<R>::operator-= (a);return *this;}
   KN& operator -=(const Add_KN_<R> & u)  
        { KN_<R>::operator-=(u);return *this;}
   KN& operator -=(const DotStar_KN_<R> & u)  
        { KN_<R>::operator-=(u);return *this;}
   KN& operator -=(const DotSlash_KN_<R> & u)  
        { KN_<R>::operator-=(u);return *this;}
   KN& operator -=(const Sub_KN_<R> & u)  
        { KN_<R>::operator-=(u);return *this;}
   KN& operator -=(const Mulc_KN_<R> & u)  
        { KN_<R>::operator-=(u);return *this;}
   KN& operator -=(const Add_Mulc_KN_<R> & u)  
        { KN_<R>::operator-=(u);return *this;}
   KN& operator -=(const if_arth_KN_<R> & u)  
        { KN_<R>::operator-=(u);return *this;}
   KN& operator -=(const Mul_KNM_KN_<R> & u) 
        { KN_<R>::operator-=(u);return *this;}
 
   KN& operator +=(const_R a)  
        { KN_<R>::operator += (a);return *this;}
   KN& operator += (const KN_<R>& a)  
        { KN_<R>::operator+= (a);return *this;}
   KN& operator +=(const Add_KN_<R> & u)  
        { KN_<R>::operator+=(u);return *this;}
   KN& operator +=(const DotStar_KN_<R> & u)  
        { KN_<R>::operator+=(u);return *this;}
   KN& operator +=(const DotSlash_KN_<R> & u)  
        { KN_<R>::operator+=(u);return *this;}
   KN& operator +=(const Sub_KN_<R> & u)  
        { KN_<R>::operator+=(u);return *this;}
   KN& operator +=(const Mulc_KN_<R> & u)  
        { KN_<R>::operator+=(u);return *this;}
   KN& operator +=(const Add_Mulc_KN_<R> & u)  
        { KN_<R>::operator+=(u);return *this;}
   KN& operator +=(const if_arth_KN_<R> & u)  
        { KN_<R>::operator+=(u);return *this;}
   KN& operator +=(const Mul_KNM_KN_<R> & u) 
        { KN_<R>::operator+=(u);return *this;}
        

   KN& operator/=(const_R a)  
        { KN_<R>::operator/=(a);return *this;}
   KN& operator /= (const KN_<R>& a)  
        { KN_<R>::operator/= (a);return *this;}
   KN& operator /=(const Add_KN_<R> & u)  
        { KN_<R>::operator/=(u);return *this;}
   KN& operator /=(const Sub_KN_<R> & u)  
        { KN_<R>::operator/=(u);return *this;}
   KN& operator /=(const Mulc_KN_<R> & u)  
        { KN_<R>::operator/=(u);return *this;}
   KN& operator /=(const Add_Mulc_KN_<R> & u)  
        { KN_<R>::operator/=(u);return *this;}
   KN& operator /=(const if_arth_KN_<R> & u)  
        { KN_<R>::operator/=(u);return *this;}
        
   KN& operator /=(const Mul_KNM_KN_<R> & u) 
        { KN_<R>::operator/=(u);return *this;}
        
   KN& operator*=(const_R a)  
        { KN_<R>::operator*=(a);return *this;}
   KN& operator*=(const KN_<const_R>& a)  
        { KN_<R>::operator*= (a);return *this;}
   KN& operator *=(const Add_KN_<R> & u)  
        { KN_<R>::operator*=(u);return *this;}
   KN& operator *=(const Sub_KN_<R> & u)  
        { KN_<R>::operator*=(u);return *this;}
   KN& operator *=(const Mulc_KN_<R> & u)  
        { KN_<R>::operator*=(u);return *this;}
   KN& operator *=(const Add_Mulc_KN_<R> & u)  
        { KN_<R>::operator*=(u);return *this;}
   KN& operator *=(const if_arth_KN_<R> & u)  
        { KN_<R>::operator*=(u);return *this;}
   KN& operator *=(const Mul_KNM_KN_<R> & u) 
        { KN_<R>::operator*=(u);return *this;}
        
  
  template<class I,class T> KN& operator =   (const KN_ITAB<T ,I> & ui)
     {  KN_<R>::operator =(ui);  return *this;}
  template<class I,class T> KN& operator +=   (const KN_ITAB<T ,I> & ui)
     {  KN_<R>::operator +=(ui);  return *this;}
  template<class I,class T> KN& operator -=   (const KN_ITAB<T ,I> & ui)
     {  KN_<R>::operator -=(ui);  return *this;}
  template<class I,class T> KN& operator *=   (const KN_ITAB<T ,I> & ui)
     {  KN_<R>::operator *=(ui);  return *this;}
  template<class I,class T> KN& operator /=   (const KN_ITAB<T ,I> & ui)
     {  KN_<R>::operator /=(ui);  return *this;}
        
        
  //  two opertor to cast to an array of constant      
//    operator KN_<const_R> & ()  
//          { return *  (KN_<const_R>*) this;}
//    operator KN_<const_R> const & ()  const 
//          { return *(const KN_<const_R>*) this;}
//    operator KN<const_R> & () 
//          { return   (KN<const_R> &) *this;}
//    operator KN<const_R> const & ()  const 
//          { return (const KN<const_R>& ) *this;}
  void init(long nn) {this->n=nn;this->step=1;this->next=-1;this->v=new R[nn];}
  void init() {this->n=0;this->step=1;this->next=-1;this->v=0;}
  void init(const KN_<R> & a){init(a.N()); operator=(a);}
  void resize(long nn) {
    if ( nn != this->n) 
     {
       R *vo=this->v;
       long no=std::min(this->n,nn), so=this->step;
       ShapeOfArray::init(nn);
       this->v=new R[this->n];
       // copy
       if(this->v && vo) 
         for(long i=0,j=0;j<no;i++,j+=so) 
           this->v[i]=vo[j]; 
        delete [] vo;} }
  void destroy(){delete [] this->v; this->v=0;this->n=0;}
};

//  Array with 2 indices
//  ---------------------

template<class R>
class KNM: public KNM_<R>{ public:

  KNM(long n,long m) 
        :KNM_<R>(new R[n*m],n,m){}
   KNM(const KNM<R> & u)  // PB si stepi ou stepj nulle
        :KNM_<R>(new R[u.size()],u.N(),u.M()) 
       { KN_<R>::operator=(u);}
  explicit KNM(const KNM_<R> & u)
        :KNM_<R>(new R[u.size()],u.N(),u.M()) 
        { KNM_<R>::operator=(u);}
        
  ~KNM(){delete [] this->v;}
  
   KNM& operator=(const KNM_<const_R> & u)   
        { KNM_<R>::operator=(u);return *this;}
   KNM& operator=(const_R a)                 
        { KNM_<R>::operator=(a);return *this;}
   KNM& operator+=(const_R  a)               
        { KNM_<R>::operator+=(a);return *this;}
   KNM& operator-=(const_R  a)               
        { KNM_<R>::operator-=(a);return *this;}
   KNM& operator/=(const_R  a)               
        { KNM_<R>::operator/=(a);return *this;}
   KNM& operator*=(const_R  a)               
        { KNM_<R>::operator*=(a);return *this;}
   KNM& operator+=(const KNM_<const_R> & u)  
        { KNM_<R>::operator+=(u);return *this;}
   KNM& operator-=(const KNM_<const_R> & u)  
        { KNM_<R>::operator-=(u);return *this;}

   KNM& operator/=(const KNM_<const_R> & u)  
        { KNM_<R>::operator/=(u);return *this;}
   KNM& operator*=(const KNM_<const_R> & u)  
        { KNM_<R>::operator*=(u);return *this;}


   KNM &operator  =(const outProduct_KN_<R> & u)
        { KNM_<R>::operator =(u);return *this;}
   KNM &operator +=(const outProduct_KN_<R> & u)
        { KNM_<R>::operator+=(u);return *this;}
   KNM &operator -=(const outProduct_KN_<R> & u)
        { KNM_<R>::operator-=(u);return *this;}
   KNM &operator /=(const outProduct_KN_<R> & u) 
        { KNM_<R>::operator/=(u);return *this;}
   KNM &operator *=(const outProduct_KN_<R> & u)
        { KNM_<R>::operator*=(u);return *this;}
  
        
  //  two opertors to cast to un array of constant        
//    operator KNM_<const_R> & ()  
//          { return *  (KNM_<const_R>*) this;}
//    operator KNM_<const_R> const & ()  const 
//          { return *(const KNM_<const_R>*) this;}

//    operator KNM<const_R> & ()  
//          { return *  (KNM<const_R>*) this;}
//    operator KNM<const_R> const & ()  const 
//          { return *(const KNM<const_R>*) this;}
  void init(long nn,long mm) {
    ShapeOfArray::init(nn*mm);
    this->shapei.init(nn,1,nn);
    this->shapej.init(mm,nn,1),
    this->v=new R[nn*mm];}
    
  void resize(long nn,long mm) {     
    long kk=nn*mm;
    long n = this->shapei.n;
    long m = this->shapej.n;
    
    if( n !=nn && m != mm) 
     {
    KNM_ <R> old(*this); 
    long no=std::min(n,nn);
    long mo=std::min(m,mm);
    R *vo=this->v;
    
    // new mat 
    ShapeOfArray::init(kk);
    this->v=new R[this->n];
    this->shapei.init(nn,1,nn);
    this->shapej.init(mm,nn,1);
    
    if(this->v && vo)  // copy
        (*this)(SubArray(no),SubArray(mo)) = old(SubArray(no),SubArray(mo));
      
    delete []vo;
    }
        
   }
    
  void destroy(){delete [] this->v;this->n=0 ;}

};

//  Array with 3 indices
//  ---------------------
template<class R>
class KNMK: public KNMK_<R>{ public:

  KNMK(long n,long m,long k) 
     :KNMK_<R>(new R[n*m*k],n,m,k){}
  explicit KNMK(const KNMK_<R> & u)
     :KNMK_<R>(new R[u.size()],u.N(),u.M(),u.K()) 
     { KNMK_<R>::operator=(u);}
   KNMK(const KNMK<R> & u)
     :KNMK_<R>(new R[u.size()],u.N(),u.M(),u.K()) 
     { KNMK_<R>::operator=(u);}
     
  ~KNMK(){delete [] this->v;}
  
  KNMK& operator=(const KNMK_<const_R> & u)   
     { KNMK_<R>::operator=(u);return *this;}
  KNMK& operator=(const_R a)                  
     { KNMK_<R>::operator=(a);return *this;}
  KNMK& operator+=(const_R  a)                
     { KNMK_<R>::operator+=(a);return *this;}
  KNMK& operator-=(const_R  a)                
     { KNMK_<R>::operator-=(a);return *this;}
  KNMK& operator/=(const_R  a)                
     { KNMK_<R>::operator/=(a);return *this;}
  KNMK& operator*=(const_R  a)                
     { KNMK_<R>::operator*=(a);return *this;}
  KNMK& operator+=(const KNMK_<const_R> & u)  
     { KNMK_<R>::operator+=(u);return *this;}
  // ici jd 
  KNMK& operator-=(const KNMK_<const_R> & u)  
     { KNMK_<R>::operator-=(u);return *this;}
  KNMK& operator*=(const KNMK_<const_R> & u)   
     { KNMK_<R>::operator*=(u);return *this;}
  KNMK& operator/=(const KNMK_<const_R> & u)   
     { KNMK_<R>::operator/=(u);return *this;}
     
//  two opertor to cast to un array of constant          
//    operator KNMK_<const_R> & ()  
//       { return *  (KNMK_<const_R>*) this;}
//    operator KNMK_<const_R> const & ()  const 
//       { return *(const KNMK_<const_R>*) this;}  

//    operator KNMK<const_R> & ()  
//       { return *  (KNMK<const_R>*) this;}
//    operator KNMK<const_R> const & ()  const 
//       { return *(const KNMK<const_R>*) this;}  
};

//  -------------  optimization ---------------------
template<class R> 
class conj_KN_{public:
  const KN_<const_R> & a;
  conj_KN_(const KN_<const_R> & aa) : a(aa){}
};


inline const KN_<long> conj(const KN_<long> &a){ return a;}
inline const KN_<double> conj(const KN_<double> &a){ return a;}
inline const KN_<float> conj(const KN_<float> &a){ return a;}

//template<class R> conj_KN_<R> conj(const KN<R> &a){ return a;}
template<class R> conj_KN_<R> conj(const KN_<R> &a){ return a;}

template<class R> 
class DotStar_KN_{public: 
  const KN_<const_R>  a; const KN_<const_R>  b;
  DotStar_KN_(const KN_<const_R> & aa,const KN_<const_R> & bb) : a(aa),b(bb)  {}
 }; 

 
template<class R> 
class DotSlash_KN_{public: 
  const KN_<const_R>  a; const KN_<const_R>  b;
  DotSlash_KN_(const KN_<const_R> & aa,const KN_<const_R> & bb) : a(aa),b(bb)  {}
 }; 

template<class R> 
class Add_KN_{public: 
  const KN_<const_R>  a; const KN_<const_R>  b;
  Add_KN_(const KN_<const_R> & aa,const KN_<const_R> & bb) 
     : a(aa),b(bb)  { K_throwassert(SameShape(a,b));}
 };  
 
template<class R> 
class Sub_KN_{public: 
  const KN_<const_R>  a; const KN_<const_R>  b;
  Sub_KN_(const KN_<const_R> & aa,const KN_<const_R> & bb) 
    : a(aa),b(bb) { K_throwassert(SameShape(a,b));}
 };
 
template<class R> 
class Mulc_KN_ { public: 
  const KN_<const_R>  a;  const_R  b;
  Mulc_KN_(const KN_<const_R> & aa,const_R  bb) : a(aa),b(bb) {}
  Mulc_KN_(const Mulc_KN_<R> & aa,const_R  bb) : a(aa.a),b(aa.b*bb) {}
  Mulc_KN_ operator-() const {return Mulc_KN_(a,-b);}
  outProduct_KN_<R> operator*(const  TKN_<double> & bb)
{  return outProduct_KN_<R>(a,bb,b);} 

 };  

template<class R> 
class Add_Mulc_KN_ { public:
  const KN_<const_R> a,b;
  const R ca,cb; 
  Add_Mulc_KN_(const Mulc_KN_<R> & aa,const Mulc_KN_<R> & bb)  
        : a(aa.a),b(bb.a),ca(aa.b),cb(bb.b) { K_throwassert(SameShape(a,b));}
  Add_Mulc_KN_(const Mulc_KN_<R> & aa,const KN_<const_R> & bb,const R cbb) 
        : a(aa.a),b(bb),ca(aa.b),cb(cbb)  { K_throwassert(SameShape(a,b));}
  Add_Mulc_KN_(const KN_<const_R> & aa,const R caa,const KN_<const_R> & bb,const R cbb) 
        : a(aa),b(bb),ca(caa),cb(cbb) { K_throwassert(SameShape(a,b));}
 };  

template<class R> 
class if_arth_KN_ { public:
  const KN_<const_R> a,b,c;
  if_arth_KN_(const KN_<R> & aa,const KN_<R> & bb,const KN_<R> & cc)  
        : a(aa),b(bb),c(cc){ K_throwassert(SameShape(a,b)&&SameShape(a,c));}
 };  



template<class R> 
class Mul_KNM_KN_ { public:
  const KNM_<const_R> &A;
  const KN_<const_R> &b;
  Mul_KNM_KN_(const  KNM_<const_R>  &aa,const KN_<const_R>  &bb)  
        : A(aa),b(bb) {K_throwassert(SameShape(A.shapej,b));} 
};


ostream & operator<<(ostream & f,const ShapeOfArray & s);

template<class R> ostream & operator<<(ostream & f,const KN_<const_R>   & v);
template<class R> ostream & operator<<(ostream & f,const KNM_<const_R>  & v);
template<class R> ostream & operator<<(ostream & f,const KNMK_<const_R> & v);
template<class R> inline ostream & operator<<(ostream & f,const KN<const_R>   & v) 
    { return f << (const KN_<const_R> &) v;}
template<class R> inline ostream & operator<<(ostream & f,const KNM<const_R>  & v) 
    { return f << (const KNM_<const_R> &) v;}
template<class R> inline ostream & operator<<(ostream & f,const KNMK<const_R> & v) 
    { return f << (const KNMK_<const_R> &) v;}


template<class R> inline Add_KN_<R> operator+(const KN_<const_R> &a,const KN_<const_R> &b) 
    { return Add_KN_<R>(a,b);}
template<class R> inline Sub_KN_<R> operator-(const KN_<const_R> &a,const KN_<const_R> &b) 
    { return Sub_KN_<R>(a,b);}
template<class R> inline Mulc_KN_<R> operator*(const KN_<const_R> &a,const R &b) 
    { return Mulc_KN_<R>(a,b);}
template<class R> inline Mulc_KN_<R> operator/(const KN_<const_R> &a,const R &b) 
    { return Mulc_KN_<R>(a,1/b);}
template<class R> inline Mulc_KN_<R> operator*(const R &b,const KN_<const_R> &a) 
    { return Mulc_KN_<R>(a,b);}
template<class R> inline Mulc_KN_<R> operator-(const KN_<const_R> &a) 
    { return Mulc_KN_<R>(a,-1);}



template<class R> inline Add_Mulc_KN_<R> operator+(const  Mulc_KN_<R>& a,const Mulc_KN_<R> &b) 
    { return Add_Mulc_KN_<R>(a,b);}
template<class R> inline Add_Mulc_KN_<R> operator-(const  Mulc_KN_<R>& a,const Mulc_KN_<R> &b) 
    { return Add_Mulc_KN_<R>(a,b.a,-b.b);}

template<class R> inline Add_Mulc_KN_<R> operator+(const  Mulc_KN_<R>& a,const KN_<const_R> &b) 
   { return Add_Mulc_KN_<R>(a,b,R(1));}
template<class R> inline Add_Mulc_KN_<R> operator-(const  Mulc_KN_<R>& a,const KN_<const_R> &b) 
   { return Add_Mulc_KN_<R>(a,b,R(-1));}

template<class R> inline Add_Mulc_KN_<R> operator+(const KN_<const_R> & b,const  Mulc_KN_<R>& a) 
   { return Add_Mulc_KN_<R>(a,b,R(1));}

// modif FH mars 2007 
template<class R> inline Add_Mulc_KN_<R> operator-(const KN_<const_R> & a,const  Mulc_KN_<R>& b) 
   { return Add_Mulc_KN_<R>(a,R(1),b.a,-b.b);}// modif FH mars 2007  

template<class R> inline Mul_KNM_KN_<R> operator*(const  KNM_<const_R> & A,const  KN_<const_R> & b) 
    { return Mul_KNM_KN_<R>(A,b);}


template<class R> inline bool  SameShape(const ShapeOfArray & a,const Add_Mulc_KN_<R> & b) 
           { return SameShape(a,b.a) ;} 
template<class R> inline bool  SameShape(const ShapeOfArray & a,const if_arth_KN_<R> & b) 
           { return SameShape(a,b.a) ;} 
template<class R> inline bool  SameShape(const ShapeOfArray & a,const Add_KN_<R> & b) 
           { return SameShape(a,b.a) ;} 
template<class R> inline bool  SameShape(const ShapeOfArray & a,const Sub_KN_<R> & b) 
           { return SameShape(a,b.a) ;} 
template<class R> inline bool  SameShape(const ShapeOfArray & a,const Mulc_KN_<R> & b) 
           { return SameShape(a,b.a) ;} 
template<class R> inline bool  SameShape(const ShapeOfArray & a,const DotStar_KN_<R> & b) 
           { return SameShape(a,b.a) ;} 
template<class R> inline bool  SameShape(const ShapeOfArray & a,const DotSlash_KN_<R> & b) 
           { return SameShape(a,b.a) ;} 
template<class R> inline bool  SameShape(const ShapeOfArray & a,const Mul_KNM_KN_<R> & b) 
           { return a.n==b.A.N() ;} 
 inline bool  SameShape(const ShapeOfArray & a,const VirtualMatrice<double>::plusAx & b) 
           { return true ;} //  pas de test car la matrice peut etre rectangulaire
 inline bool  SameShape(const ShapeOfArray & a,const VirtualMatrice<double>::plusAtx & b) 
           { return true ;} //  pas de test car la matrice peut etre rectangulaire
 inline bool  SameShape(const ShapeOfArray & a,const VirtualMatrice<complex<double> >::plusAx & b) 
           { return true ;} //  pas de test car la matrice peut etre rectangulaire
 inline bool  SameShape(const ShapeOfArray & a,const VirtualMatrice<complex<double> >::plusAtx & b) 
           { return true ;} //  pas de test car la matrice peut etre rectangulaire

 inline bool  SameShape(const ShapeOfArray & a,const double) 
           { return true;} 
 inline bool  SameShape(const ShapeOfArray & a,const complex<double>) 
           { return true;} 
 inline bool  SameShape(const ShapeOfArray & a,const complex<float>) 
           { return true;}            

template<class R>
 inline bool SameShape(KNM<R>& m, const outProduct_KN_<R>& p)
 { return p.a.N()>=m.N() && m.M()>=p.b.N(); } 

template<class R> inline long SameAdress(const KN_<R> &a, const KN_<R> &b) { return &a[0]==&b[0];}
// bof -bof 
//template<class R> inline
//  KN_<R>::operator KN<R> &() { return *(KN<R> *) (void *) this;}
//template<class R> inline
//  KN_<R>::operator const KN<R> &() const { return *(const KN<R> *) ( const void *) this;}

//  operateur y=Ax-b ou y=Ax + b pour le GC
template<class R> 
   PplusQ< typename VirtualMatrice<R>::plusAx, Mulc_KN_<R> >  operator-(const typename VirtualMatrice<R>::plusAx & A,const KN_<R> & B)  
    { return PplusQ< typename VirtualMatrice<R>::plusAx, Mulc_KN_<R> >(A,Mulc_KN_<R>(B,R(-1.)));}

template<class R> 
   PplusQ< typename VirtualMatrice<R>::plusAx, KN_<R> >  operator+(const typename VirtualMatrice<R>::plusAx & A,const KN_<R> & B)  
    { return PplusQ< typename VirtualMatrice<R>::plusAx, KN_<R> >(A,B);}

template<class R> 
   PplusQ< typename VirtualMatrice<R>::plusAx, Mulc_KN_<R> >  operator-(const typename VirtualMatrice<R>::plusAx & A,const KN<R> & B)  
    { return PplusQ< typename VirtualMatrice<R>::plusAx, Mulc_KN_<R> >(A,Mulc_KN_<R>(B,R(-1.)));}

template<class R> 
   PplusQ< typename VirtualMatrice<R>::plusAx, KN_<R> >  operator+(const typename VirtualMatrice<R>::plusAx & A,const KN<R> & B)  
    { return PplusQ< typename VirtualMatrice<R>::plusAx, KN_<R> >(A,B);}
    

template<class R>
KN_<R> diagonal(const KNM<R> & A) { 
  K_throwassert(A.N() == A.M()); 
  return KN_<R>(A,SubArray(A.N(),0,A.N()+1));}

// to def  inv permutation FH mars 2006 
class Inv_KN_long{ public:
  KN_<long>  t;
  Inv_KN_long(const KN_<long> & v)
   : t(v) {}
  Inv_KN_long( KN_<long> const  * & v)
   : t(*v) {}
  operator const KN_<long> & () const {return t;}
};


template<class R,typename A,typename B=R> class  F_KN_ 
{ 
  public: 
  A (*f)(B);
  KN_<R> a;
  long N() const {return a.N();}
  F_KN_( A (*ff)(B),const KN_<R> & aa): f(ff),a(aa) {}
  A operator[](long i) const { return f(a[i]);}
  bool check(long n)  const { return  n <= a.N() || a.constant(); }
  bool constant() const {return a.constant();}
}; 

template<class R,typename A,typename B>
inline bool  SameShape(const ShapeOfArray & a,const F_KN_<R,A,B>  & b) 
           { return  !a.step || b.constant()  || a.n == b.N() ;} 
           
#include "RNM_tpl.hpp"
#ifdef K_throwassert
#undef K_throwassert
#endif
#endif
