#ifndef KNM_H_
#define KNM_H_

// une tentative qui ne marche pas 
// de tableau constant 
#include <complex>
#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;
#define const_R   R
#ifdef CHECK_KN
#include <cstdlib>
inline void Check_Kn(const char * str,const char * file,int line)
{
 cerr << "CHECK_KN: " << str << " in file: " << file << ", line " << line <<endl;
 assert(0);
 abort();
}
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
//    \u00f5 un 2 ou 3 indices est ou non consecutif en memoire
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
//  
// 
template<class R> class KNMK_ ;
template<class R> class KNM_ ;
template<class R> class KN_ ;
template<class R> class TKN_ ; // KN_ transpose 

template<class R> class KNMK ;
template<class R> class KNM ;
template<class R> class KN ;

template<class R> class Add_KN_;
template<class R> class Sub_KN_;
template<class R> class Mulc_KN_;
template<class R> class Add_Mulc_KN_;
template<class R> class Mul_KNM_KN_; 
template<class R> class DotStar_KN_;
template<class R> class DotSlash_KN_;
template<class R> class outProduct_KN_;
template<class R> class VirtualMatrice;


// gestion des erreur interne --
#ifndef InternalError
typedef void (* TypeofInternalErrorRoutine)(const char *) ;
inline TypeofInternalErrorRoutine &InternalErrorRoutinePtr()
{ 
  static TypeofInternalErrorRoutine routine=0;
  return routine;
}

inline void InternalError(const char * str) {
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

template<class R> 
struct  VirtualMatrice { public:

  //  y += A x 
  virtual void addMatMul(const KN_<R> &  x, KN_<R> & y) const =0; 
  virtual void addMatTransMul(const KN_<R> &  x, KN_<R> & y) const 
    { InternalError("VirtualMatrice::addMatTransMul not implemented "); }
  virtual bool WithSolver() const {return false;} // by default no solver          
  virtual void Solve( KN_<R> &  x,const KN_<R> & b) const 
    { InternalError("VirtualMatrice::solve not implemented "); } 
  struct  plusAx { const VirtualMatrice * A; const KN_<R> &  x;
   plusAx( const VirtualMatrice * B,const KN_<R> &  y) :A(B),x(y) {} };
   plusAx operator*(const KN_<R> &  x) const {return plusAx(this,x);}
  struct  plusAtx { const VirtualMatrice * A; const KN_<R> &  x;
   plusAtx( const VirtualMatrice * B,const KN_<R> &  y) :A(B),x(y) {} };
  struct  solveAxeqb { const VirtualMatrice * A; const KN_<R> &  b;
   solveAxeqb( const VirtualMatrice * B,const KN_<R> &  y) :A(B),b(y) {} };
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
  void init(long nn,long s=1,long nextt=-1) { n=nn; step=s; next=next;}         
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
  
  R min() const ;
  R max() const ;
  R sum() const ;
  template<class T> long last(const T &) const;
  template<class T> long first(const T &) const;
  
  KN_ & map(R (*f)(R ));
  
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

   KN_& operator =(const Add_KN_<R> & u) ;
   KN_& operator+=(const Add_KN_<R> & u) ;
   KN_& operator-=(const Add_KN_<R> & u) ;
   KN_& operator*=(const Add_KN_<R> & u) ;
   KN_& operator/=(const Add_KN_<R> & u) ;
  
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
  
   KN_& operator =(const Mul_KNM_KN_<R> & u) ; 
   KN_& operator+=(const Mul_KNM_KN_<R> & u) ; 
   KN_& operator-=(const Mul_KNM_KN_<R> & u) ; 
   KN_& operator*=(const Mul_KNM_KN_<R> & u) ; 
   KN_& operator/=(const Mul_KNM_KN_<R> & u) ; 
  
 //  KN_& operator =(const MatriceCreuseMulKN_<R> & ) ;
 //  KN_& operator +=(const MatriceCreuseMulKN_<R> & ) ;
   KN_& operator =(const typename VirtualMatrice<R>::plusAx & Ax)  
    {*this=0; Ax.A->addMatMul(Ax.x,*this);return *this;}
   KN_& operator =(const typename VirtualMatrice<R>::plusAtx & Ax)  
    {*this=0; Ax.A->addMatTransMul(Ax.x,*this);return *this;}
   KN_& operator +=(const typename VirtualMatrice<R>::plusAx & Ax)  
    {  Ax.A->addMatMul(Ax.x,*this);return *this;}
   KN_& operator +=(const typename VirtualMatrice<R>::plusAtx & Ax)  
    {  Ax.A->addMatTransMul(Ax.x,*this);return *this;}
   KN_& operator =(const typename VirtualMatrice<R>::solveAxeqb & Ab)  
    {*this=0; Ab.A->Solve(*this,Ab.b);return *this;}
  
//   KN_& operator =(const MatriceCreuseDivKN_<R> &)  ;

 friend   ostream & operator<< <R>(ostream & f,const KN_<const_R> & v)  ;
 
  KN_(R *u,const ShapeOfArray & s):ShapeOfArray(s),v(u){}
  KN_(R *u,long nn,long s):ShapeOfArray(nn,s),v(u){}
  KN_(R *u,long nn,long s,long nextt):ShapeOfArray(nn,s,nextt),v(u){}
  KN_(R *u,long nn):ShapeOfArray(nn),v(u){}


  TKN_<R>  t() ; // transpose
  const  TKN_<R>  t() const ; // transpose

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
  long IsVector1() const  {  return (shapei.n*shapej.n) ==  n ;} 
  long N() const {return shapei.n;}
  long M() const {return shapej.n;}  
  long size() const { return shapei.n*shapej.n;}
  
  KNM_(R* u,const ShapeOfArray & s,
            const ShapeOfArray & si,
            const ShapeOfArray & sj)
             : KN_<R>(u,s),shapei(si),shapej(sj){} 
  KNM_(R* u,long n,long m)
             : KN_<R>(u,ShapeOfArray(n*m)),shapei(n,1,n),shapej(m,n,1){}
  KNM_(R* u,long n,long m, long s)
             : KN_<R>(u,ShapeOfArray(n*m,s)),shapei(n,1,n),shapej(m,n,1){}                     
  KNM_(KN_<R> u,long n,long m) 
             : KN_<R>(u,ShapeOfArray(m*n)),shapei(n,1,n),shapej(m,n,1){ }
             
  KNM_(const KN_<R> &u,const ShapeOfArray & si,const ShapeOfArray & sj,long offset=0) 
             : KN_<R>(&u[offset],si.last()+sj.last()+1,u.step),shapei(si),shapej(sj) 
             {K_throwassert( offset>=0 && n+ (v-(R*)u) <= u.n);}
  KNM_(const KN_<R> &u,const ShapeOfArray & si,const ShapeOfArray & sj,long offset,long nnext) 
             : KN_<R>(&u[offset],si.last()+sj.last()+1,u.step,nnext),shapei(si),shapej(sj) 
             {K_throwassert( offset>=0 && n+ (v-(R*)u) <= u.n);}
  
  KNM_(KNM_<R> U,const SubArray & si,const SubArray & sj)  
             :KN_<R>(U,SubArray(U.ij(si.len1(),sj.len1())+1,U.ij(si.start,sj.start))),
                     shapei(U.shapei,si),shapej(U.shapej,sj){} 

  KNM_(KNM_<R> U,const SubArray & sa,const SubArray & si,const SubArray & sj)  
             :KN_<R>(U,SubArray(sa)),shapei(U.shapei,si),shapej(U.shapej,sj){} 

  KNM_(const KNM_<R> & u) 
             :KN_<R>(u),shapei(u.shapei),shapej(u.shapej) {}

  KNM_ operator()(const SubArray & sa,const SubArray & sb) const 
            { return KNM_(*this,sa,sb);} // sub array 
  
  long ij(const long i,const long j) const   
            { return shapei.index(i)+shapej.index(j);}
  long indexij(long i,long j)        const   
            { return index(shapei.index(i)+shapej.index(j));}
  R & operator()(long i,long j)     const   
            { return v[indexij(i,j)];}
  R & operator()(int i,int j)     const   
            { return v[indexij(i,j)];}
            
            
  KN_<R> operator()(const char,long j    )  const   // une colonne j  ('.',j)
            { return KN_<R>(&v[index(shapej.index(j))],shapei*step);} 
  KN_<R> operator()(long i    ,const char)  const   // une ligne i  (i,'.')
            { return KN_<R>(&v[index(shapei.index(i))],shapej*step);}  
  KN_<R> operator()(const char,int j    )  const   // une colonne j  ('.',j)
            { return KN_<R>(&v[index(shapej.index(j))],shapei*step);} 
  KN_<R> operator()(int i    ,const char)  const   // une ligne i  (i,'.')
            { return KN_<R>(&v[index(shapei.index(i))],shapej*step);}  
  KN_<R> operator()(const char,const char)  const   // tous 
            { return *this;}
  KNM_<R> t() const
    { return KNM_<R>(v,*this,shapej,shapei);} //  before  { return KNM_<R>(*this,shapej,shapei,v);}
  
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
  KNM_& operator++() {v += next;return *this;} // ++U
  KNM_& operator--() {v -= next;return *this;} // ++U
  KNM_  operator++(int ){KNM_<R> old=*this;v = v +next;return old;} // U++
  KNM_  operator--(int ){KNM_<R> old=*this;v = v -next;return old;} // U--


 friend class KN_<R>;   
// friend class KNM_<R>;   
 friend class KNMK_<R>;   
 friend class KN<R>;   
 friend class KNM<R>;   
 friend class KNMK<R>;   
};

template<class R>
struct TKN_:public KN_<R> {
    TKN_(const KN_<R> &x) : KN_<R>(x) {}
};

template<class R>
TKN_<R>  KN_<R>::t() { return *this;} // transpose

template<class R>
const TKN_<R>  KN_<R>::t() const { return *this;} // transpose

template<class R>
struct outProduct_KN_ {
    const KN_<R> & a,&b;
    R c;
    outProduct_KN_(const KN_<R> & aa, const KN_<R> &bb,R cc=1.) : a(aa),b(bb),c(cc) {}
    outProduct_KN_ operator * (R cc) { return outProduct_KN_(a,b,c*cc);}    
};

template<class R> 
outProduct_KN_<R> operator*(const KN_<R> &a,const TKN_<R> &b) 
{ return outProduct_KN_<R>(a,b);}

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
  long IsVector1() const {  return (shapei.n*shapej.n*shapek.n) == n ;} 
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
    
  KNMK_(R* u,const long n,const long m,const long k)
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

    
  long ijk(const long i,const long j,const long k) const 
              { return shapei.index(i)+shapej.index(j)+shapek.index(k);}
  long indexijk(long i,long j,long k) const 
              {return index(shapei.index(i)+shapej.index(j)+shapek.index(k));} 
                           
  R & operator()(long i,long j,long k)   const   {return v[indexijk(i,j,k)];}
  R & operator()(int i,int j,int k)   const   {return v[indexijk(i,j,k)];}
  
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
  KN(const long nn) : KN_<R>(new R[nn],nn)         {} 
  KN(const long nn, R * p) : KN_<R>(new R[nn],nn)  
    { KN_<R>::operator=(KN_<R>(p,nn));}
  KN(const long nn,R (*f)(long i) ) : KN_<R>(new R[nn],nn) 
        {for(long i=0;i<n;i++) v[i]=f(i);}  
  KN(const long nn,const  R & a) : KN_<R>(new R[nn],nn) 
        { KN_<R>::operator=(a);} 
  KN(const long nn,long s,const  R  a) : KN_<R>(new R[nn],nn,s) 
        { KN_<R>::operator=(a);} 
  template<class S>   KN(const KN_<S> & s):KN_<R>(new R[s.n],s.n) 
        {for (long i=0;i<n;i++) v[i] = s[i];}
  template<class S>  KN(const KN_<S> & s,R (*f)(S )):KN_<R>(new R[s.n],s.n) 
        {for (long i=0;i<n;i++) v[i] = f(s[i]);}
  KN(const KN<R> & u):KN_<R>(new R[u.n],u.n)
        { KN_<R>::operator=(u);}
  //  explicit KN(const KN_<R> & u):KN_<R>(new R[u.n],u.n) 
  //      { KN_<R>::operator=(u);}
        
  ~KN(){delete [] v;}
  
   KN& operator =(const_R a)  
        { if(unset()) set(new R[1],1,0,0); KN_<R>::operator= (a);return *this;}
   KN& operator =(const KN_<R>& a)  
        { if(unset()) set(new R[a.N()],a.N()); KN_<R>::operator= (a);return *this;}                
   KN& operator =(const KN<R>& a)  
        { if(unset()) set(new R[a.N()],a.N()); KN_<R>::operator= (a);return *this;}                
   KN& operator =(const Add_KN_<R> & u)  
        { if(unset()) set(new R[u.a.N()],u.a.N());KN_<R>::operator=(u);return *this;}
   KN& operator =(const DotStar_KN_<R> & u)  
        { if(unset()) set(new R[u.a.N()],u.a.N());KN_<R>::operator=(u);return *this;}
   KN& operator =(const DotSlash_KN_<R> & u)  
        { if(unset()) set(new R[u.a.N()],u.a.N());KN_<R>::operator=(u);return *this;}
   KN& operator =(const Sub_KN_<R> & u)  
        { if(unset()) set(new R[u.a.N()],u.a.N());KN_<R>::operator=(u);return *this;}
   KN& operator =(const Mulc_KN_<R> & u)  
        { if(unset()) set(new R[u.a.N()],u.a.N());KN_<R>::operator=(u);return *this;}
   KN& operator =(const Add_Mulc_KN_<R> & u)  
        { if(unset()) set(new R[u.a.N()],u.a.N());KN_<R>::operator=(u);return *this;}
   KN& operator =(const Mul_KNM_KN_<R> & u) 
        { if(unset()) set(new R[u.b.N()],u.b.N());KN_<R>::operator=(u);return *this;}
//   KN& operator =(const MatriceCreuseMulKN_<R> & Ax) 
//       {if(unset()) set(new R[Ax.v.N()],Ax.v.N()); KN_<R>::operator=(Ax);return *this;}
//   KN& operator +=(const MatriceCreuseMulKN_<R> & Ax) 
//       {if(unset()) set(new R[Ax.v.N()],Ax.v.N()); KN_<R>::operator+=(Ax);return *this;}
//   KN& operator =(const MatriceCreuseDivKN_<R> & A1x)  
//       { if(unset()) set(new R[A1x.v.N()],A1x.v.N());KN_<R>::operator=(A1x);return *this;}
   KN& operator =(const typename VirtualMatrice<R>::plusAx & Ax)  
        { if(unset()) set(new R[Ax.x.N()],Ax.x.N());KN_<R>::operator=(Ax);return *this;}
   KN& operator =(const typename VirtualMatrice<R>::solveAxeqb & Ab)  
        { if(unset()) set(new R[Ab.b.N()],Ab.b.N());KN_<R>::operator=(Ab);return *this;}
   KN& operator +=(const typename  VirtualMatrice<R>::plusAx & Ax)  
        { if(unset()) set(new R[Ax.x.N()],Ax.x.N());KN_<R>::operator+=(Ax);return *this;}
   KN& operator =(const typename VirtualMatrice<R>::plusAtx & Ax)  
        { if(unset()) set(new R[Ax.x.N()],Ax.x.N());KN_<R>::operator=(Ax);return *this;}
   KN& operator +=(const typename VirtualMatrice<R>::plusAtx & Ax)  
        { if(unset()) set(new R[Ax.x.N()],Ax.x.N());KN_<R>::operator+=(Ax);return *this;}
        
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
   KN& operator *=(const Mul_KNM_KN_<R> & u) 
        { KN_<R>::operator*=(u);return *this;}
        
        
  //  two opertor to cast to an array of constant      
//    operator KN_<const_R> & ()  
//          { return *  (KN_<const_R>*) this;}
//    operator KN_<const_R> const & ()  const 
//          { return *(const KN_<const_R>*) this;}
//    operator KN<const_R> & () 
//          { return   (KN<const_R> &) *this;}
//    operator KN<const_R> const & ()  const 
//          { return (const KN<const_R>& ) *this;}
  void init(long nn) {n=nn;step=1;next=-1;v=new R[nn];}
  void destroy(){delete [] v; v=0;n=0;}
};

//  Array with 2 indices
//  ---------------------

template<class R>
class KNM: public KNM_<R>{ public:

  KNM(const long n,const long m) 
        :KNM_<R>(new R[n*m],n,m){}
   KNM(const KNM<R> & u)  // PB si stepi ou stepj nulle
        :KNM_<R>(new R[u.size()],u.N(),u.M()) 
       { KN_<R>::operator=(u);}
  explicit KNM(const KNM_<R> & u)
        :KNM_<R>(new R[u.size()],u.N(),u.M()) 
        { KNM_<R>::operator=(u);}
        
  ~KNM(){delete [] v;}
  
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
    shapei.init(nn,1,nn);
    shapej.init(mm,nn,1),
    v=new R[nn*mm];}
    
  void destroy(){delete [] v;n=0 ;}

};

//  Array with 3 indices
//  ---------------------
template<class R>
class KNMK: public KNMK_<R>{ public:

  KNMK(const long n,const long m,const long k) 
     :KNMK_<R>(new R[n*m*k],n,m,k){}
  explicit KNMK(const KNMK_<R> & u)
     :KNMK_<R>(new R[u.size()],u.N(),u.M(),u.K()) 
     { KNMK_<R>::operator=(u);}
   KNMK(const KNMK<R> & u)
     :KNMK_<R>(new R[u.size()],u.N(),u.M(),u.K()) 
     { KNMK_<R>::operator=(u);}
     
  ~KNMK(){delete [] v;}
  
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
class DotStar_KN_{public: 
  const KN_<const_R> & a; const KN_<const_R> & b;
  DotStar_KN_( KN_<const_R> & aa, KN_<const_R> & bb) : a(aa),b(bb)  {}
 }; 
 
template<class R> 
class DotSlash_KN_{public: 
  const KN_<const_R> & a; const KN_<const_R> & b;
  DotSlash_KN_( KN_<const_R> & aa, KN_<const_R> & bb) : a(aa),b(bb)  {}
 }; 

template<class R> 
class Add_KN_{public: 
  const KN_<const_R> & a; const KN_<const_R> & b;
  Add_KN_(const KN_<const_R> & aa,const KN_<const_R> & bb) 
     : a(aa),b(bb)  { K_throwassert(SameShape(a,b));}
 };  
 
template<class R> 
class Sub_KN_{public: 
  const KN_<const_R> & a; const KN_<const_R> & b;
  Sub_KN_(const KN_<const_R> & aa,const KN_<const_R> & bb) 
    : a(aa),b(bb) { K_throwassert(SameShape(a,b));}
 };
 
template<class R> 
class Mulc_KN_ { public: 
  const KN_<const_R> & a;  const_R  b;
  Mulc_KN_(const KN_<const_R> & aa,const_R  bb) : a(aa),b(bb) {}
  Mulc_KN_(const Mulc_KN_<R> & aa,const_R  bb) : a(aa.a),b(aa.b*bb) {}
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
class Mul_KNM_KN_ { public:
  const KNM_<const_R> A;
  const KN_<const_R> b;
  Mul_KNM_KN_(const  KNM_<const_R> & aa,const KN_<const_R> & bb)  
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
   { return Add_Mulc_KN_<R>(a,b,1.0);}
template<class R> inline Add_Mulc_KN_<R> operator-(const  Mulc_KN_<R>& a,const KN_<const_R> &b) 
   { return Add_Mulc_KN_<R>(a,b,-1.0);}

template<class R> inline Add_Mulc_KN_<R> operator+(const KN_<const_R> & b,const  Mulc_KN_<R>& a) 
   { return Add_Mulc_KN_<R>(a,b,1.0);}
template<class R> inline Add_Mulc_KN_<R> operator-(const KN_<const_R> & b,const  Mulc_KN_<R>& a) 
   { return Add_Mulc_KN_<R>(a,b,-1.0);}
template<class R> inline Mul_KNM_KN_<R> operator*(const  KNM_<const_R> A,const  KN_<const_R> b) 
    { return Mul_KNM_KN_<R>(A,b);}


template<class R> inline bool  SameShape(const ShapeOfArray & a,const Add_Mulc_KN_<R> & b) 
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
 inline bool  SameShape(const ShapeOfArray & a,const VirtualMatrice<double>::plusAx & b) 
           { return true ;} //  pas de test car la matrice peut etre rectangulaire
 inline bool  SameShape(const ShapeOfArray & a,const VirtualMatrice<double>::plusAtx & b) 
           { return true ;} //  pas de test car la matrice peut etre rectangulaire
 inline bool  SameShape(const ShapeOfArray & a,const double) 
           { return true;} 
 inline bool  SameShape(const ShapeOfArray & a,const complex<double>) 
           { return true;} 
inline bool  SameShape(const ShapeOfArray & a,const complex<float>) 
           { return true;}            

template<class R> inline long SameAdress(const KN_<R> &a, const KN_<R> &b) { return &a[0]==&b[0];}
// bof -bof 
//template<class R> inline
//  KN_<R>::operator KN<R> &() { return *(KN<R> *) (void *) this;}
//template<class R> inline
//  KN_<R>::operator const KN<R> &() const { return *(const KN<R> *) ( const void *) this;}


template<class R>
KN_<R> diagonal(const KNM<R> & A) { 
  K_throwassert(A.N() == A.M()); 
  return KN_<R>(A,SubArray(A.N(),0,A.N()+1));}

#include "RNM_tpl.hpp"
#ifdef K_throwassert
#undef K_throwassert
#endif
#endif
