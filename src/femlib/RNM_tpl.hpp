#ifndef  RNM_tpl_
#define  RNM_tpl_

#include "RNM.hpp"

//   version du 22 nov 1999
//   Voila une debut de class vecteur + matrice 
//   les class avec termine avec un _ ne travail que sur 
//   un pointeur existant => pas de new et delete de pointeur
//   la class correspondant sans de _ genere les pointeurs
//  
//   avec ses classes on peut prendre une ligne 
//   ou une colonne d'une matrice
// -----------------------


template<class R>
void MatMul(KNM_<R> & ab, KNM_<R> &  a, KNM_<R> & b){
  // attention ne marche que si les adresses ne sont pas les memes
  long N= a.shapei.n;
  long M= a.shapej.n;
  K_throwassert(a.shapej.n == a.shapei.n);
  K_throwassert(a.shapei.n == ab.shapei.n);
  K_throwassert(b.shapej.n == ab.shapej.n);
  K_throwassert(b.v != ab.v);
  K_throwassert(a.v != ab.v);
  KN_<R> ai =a(0);
  for(long i=1;i<N;i++,++ai){
    KN_<R> bj=b[0];
    for(long j=1;j<M;j++,++bj)
      ab(i,j) = (ai , bj)  ;}
}



inline ostream & operator<<(ostream & f,const ShapeOfArray & s)
  { f << s.n ; 
     const int i10=10; 
     int prec=f.precision();
     if(prec<i10) f.precision(i10);
    if(s.step != 1)
      f << ":" << s.step ;
    if (s.step*s.n  !=  s.next ) 
      f << " n: " << setw(3) << s.next ;
    f << ",";
     if(prec<i10) f.precision(prec);
     return f;
   };


template<class R> ostream & operator<<(ostream & f,const KN_<const_R> & v)
  { //f <<  " KN_ : " << (ShapeOfArray) v << " "   <<  (const_R *) v << " :\n\t"  ;
    f << v.N() << "\t\n\t" ; f.precision(10);
    for (long i=0;i<v.N();i++)
      f   << setw(3) << v[i] << ((i % 5) == 4 ? "\n\t" : "\t");
     return f;
   };

template<class R> istream & operator>>(istream & f, KN_<R> & v)
 {
     int n;char c;
     f >> n;
     assert(fdat.good());
     assert(n==v.N());
     while (f.get(c) &&  (c!='\n' && c!='\r' ) ) 0; // eat until control (new line

     for (int i=0;i<n;i++)
       f >> v[i] ;
     assert(fdat.good());
     return f;
}

template<class R> istream & operator>>(istream & f, KN<R> & v)
 {
     int n;char c;
     f >> n;
     if (v.unset()) v.init(n);
     cout << n << " == " << v.N() << endl;
     assert(n==v.N());
     while (f.get(c) &&  (c!='\n' && c!='\r' ) ) 0; // eat until control (new line

     for (int i=0;i<n;i++)
       f >> v[i] ;
     return f;
}


template<class R> ostream & operator<<(ostream & f,const KNM_<const_R> & v)
  {  //f << " KNM_ "<<v.N()<<"x"<<v.M()<< ": " << (ShapeOfArray) v  
     //<< " i "  << v.shapei
     // << " j "  << v.shapej
     // << " " << &v(0,0) << " :\n\t";
      f.precision(10);
     f << v.N()<<' '<<v.M() /*<< "  n" << v.next<<" :"<< v.shapei.next << "," << v.shapej.next */<< "\t\n\t" ;
    for (long i=0;i<v.N();i++) {
      for (long j=0;j<v.M();j++) 
        f << " " << setw(3) << v(i,j);
       f << "\n\t";}
  return f;
    
   };

template<class R> ostream & operator<<(ostream & f,const KNMK_<const_R> & v)
  { //f << " KNM_" <<v.N()<<"x"<<v.M()<<"x"<<v.K()<< " : " << (ShapeOfArray) v  
    // << " i "  << v.shapei
    // << " j "  << v.shapej
    // << " k "  << v.shapek << endl;
    // << " " << (void *) & v(0,0,0) << "\n\t" ;
   f << v.N()<< 'x' <<v.M()<< 'x'<<v.K() << "\t:\n\t" ;    
  for (long i=0;i<v.shapei.n;i++){
    for (long j=0;j<v.shapej.n;j++){
      for (long k=0;k<v.shapek.n;k++)
	f << " " << setw(3) << v(i,j,k);
      f << "\n\t";}
    f << "\n\t";}
  return f;
    
   };

template<class R>
 R  KN_<R>::operator,(const KN_<const_R> & u) const {
    K_throwassert(u.n == n);
    R  s=0; 
    R * l(v);
    R  *r(u.v);    
    for (long i=0;i<n;i++,l += step, r += u.step) s += *l * *r;
    return s;
  }


template<class R>
R  KN_<R>::min() const {
    R minv = v[index(0)];
    for (long i=1;i<n;i++)
      minv = minv < v[index(i)] ? minv : v[index(i)];
    return minv;
  }
template<class R>
R  KN_<R>::max() const {
    R maxv = v[index(0)];
    for (long i=1;i<n;i++)
      maxv = maxv > v[index(i)] ? maxv : v[index(i)];
    return maxv;
  }
  
template<class R>
R  KN_<R>::sum() const {
    R s = v[index(0)];
    for (long i=1;i<n;i++)
      s +=  v[index(i)];
    return s;
  }

template<class R> template<class T>
long  KN_<R>::last(const T & a) const {
    for (long i=n;i-- >0;)
      if  (a(v[index(i)])) 
        return i;
    return -1;
 }
 
template<class R> template<class T>
long  KN_<R>::first(const T & a) const {
    for (long i=0;i<n;i++)
      if  (a(v[index(i)])) return i;
    return n;
 }


template<class R>
 KN_<R>&  KN_<R>::map(R (*f)(R )) {
    for (long i=0;i<n;i++)
      {  R & x(v[index(i)]);
          x =  f(x);}
   
  }


///////////////// definition des operateurs d'affectation /////////////////////////
#define oper =
#include "RNM_op.hpp"
#include "RNM_opc.hpp"
#define oper +=
#include "RNM_op.hpp"
#include "RNM_opc.hpp"
#define oper -=
#include "RNM_op.hpp"
#include "RNM_opc.hpp"
#define oper *=
#include "RNM_op.hpp"
#include "RNM_opc.hpp"
#define oper /=
#include "RNM_op.hpp"
#include "RNM_opc.hpp"

#endif
