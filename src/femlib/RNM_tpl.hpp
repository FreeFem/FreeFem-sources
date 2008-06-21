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

namespace RNM  {

template <class T> inline double norm(const T & x){return std::norm(x);}
inline double norm(double x){return x*x;} 
inline double norm(float x){return x*x;}
inline long norm(long x){return x*x;}
inline int norm(int x){return x*x;}

}

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
  for(long i=0;i<N;i++,++ai){
    KN_<R> bj=b[0];
    for(long j=0;j<M;j++,++bj)
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
    f << v.N() << "\t\n\t" ; 
    const int i10=10; 
    int prec=f.precision();
    if(prec<i10) f.precision(i10);    
    for (long i=0;i<v.N();i++)
      f   << setw(3) << v[i] << ((i % 5) == 4 ? "\n\t" : "\t");
    if(prec<i10) f.precision(prec); 
    return f;
  };

template<class R> istream & operator>>(istream & f, KN_<R> & v)
 {
     int n;char c;
     f >> n;
     ffassert(f.good());
     ffassert(n==v.N());
     while (f.get(c) &&  (c!='\n' && c!='\r' ) ) ; // eat until control (new line

     for (int i=0;i<n;i++)
      {  f >> v[i] ;
       ffassert(f.good());} // modif FH  main 2006
     return f;
}

template<class R> istream & operator>>(istream & f, KN<R> & v)
 {
     int n;char c;
     f >> n;
     if (v.unset()) v.init(n);
     cout << n << " == " << v.N() << endl;
     ffassert(n==v.N());
     while (f.get(c) &&  (c!='\n' && c!='\r' ) ) ; // eat until control (new line

     for (int i=0;i<n;i++)
       {
       f >> v[i] ;
       ffassert(f.good());}// modif FH  main 2006
     return f;
}


template<class R> ostream & operator<<(ostream & f,const KNM_<const_R> & v)
  {  //f << " KNM_ "<<v.N()<<"x"<<v.M()<< ": " << (ShapeOfArray) v  
     //<< " i "  << v.shapei
     // << " j "  << v.shapej
     // << " " << &v(0,0) << " :\n\t";
    const int i10=10; 
    int prec=f.precision();
    if(prec<i10) f.precision(i10);    
     f << v.N()<<' '<<v.M() /*<< "  n" << v.next<<" :"<< v.shapei.next << "," << v.shapej.next */<< "\t\n\t" ;
    for (long i=0;i<v.N();i++) {
      for (long j=0;j<v.M();j++) 
        f << " " << setw(3) << v(i,j);
       f << "\n\t";}
    if(prec<i10) f.precision(prec);
  return f;
    
   };

template<class R> ostream & operator<<(ostream & f,const KNMK_<const_R> & v)
  { //f << " KNM_" <<v.N()<<"x"<<v.M()<<"x"<<v.K()<< " : " << (ShapeOfArray) v  
    // << " i "  << v.shapei
    // << " j "  << v.shapej
    // << " k "  << v.shapek << endl;
    // << " " << (void *) & v(0,0,0) << "\n\t" ;
    f << v.N()<< 'x' <<v.M()<< 'x'<<v.K() << "\t:\n\t" ;    
    const int i10=10; 
    int prec=f.precision();
    if(prec<i10) f.precision(i10);    
    for (long i=0;i<v.shapei.n;i++){
      for (long j=0;j<v.shapej.n;j++){
	for (long k=0;k<v.shapek.n;k++)
	  f << " " << setw(3) << v(i,j,k);
	f << "\n\t";}
      f << "\n\t";}
    if(prec<i10) f.precision(prec);
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
 R  operator,(const KN_<const_R> & u,const conj_KN_<const_R> & vc) {
  int n=u.n;
    K_throwassert(n == vc.a.n);
    R  s=0; 
    R * l(u);
    R  *r(vc.a);
    int stepl= u.step, stepr=vc.a.step;    
    for (long i=0;i<n;i++,l += stepl, r += stepr) s += *l * conj(*r);
    return s;
  }

template<class R>
 R  operator,(const conj_KN_<const_R> & u,const KN_<const_R> & vc) {
  int n=u.a.n;
    K_throwassert(n == vc.n);
    R  s=0; 
    R * l(u.a);
    R  *r(vc);
    int stepl= u.a.step, stepr=vc.step;    
    for (long i=0;i<n;i++,l += stepl, r += stepr) s += conj(*l) * (*r);
    return s;
  }

template<class R>
 R  operator,(const KN<const_R> & u,const conj_KN_<const_R> & vc) {  return ( (KN_<R>) u,vc);}
template<class R>
 R  operator,(const conj_KN_<const_R> & u,const KN<const_R> & vc) {  return (  u, (KN_<R>) vc);}


template<class R>
R  KN_<R>::min() const {
    R minv = v[index(0)];
    for (long i=1;i<n;i++)
      minv = RNM::Min(minv, v[index(i)]) ;
    return minv;
  }
template<class R>
R  KN_<R>::max() const {
    R maxv = v[index(0)];
    for (long i=1;i<n;i++)
      maxv = RNM::Max(maxv ,v[index(i)]);
    return maxv;
  }


  
template<class R>
R  KN_<R>::sum() const {
    R s = v[index(0)];
    for (long i=1;i<n;i++)
      s +=  v[index(i)];
    //    cout << " sum = " << s << endl;
    return s;
  }

template<class R>
double  KN_<R>::norm() const {
  double s = 0.;
    for (long i=0;i<n;i++)
      s +=  RNM::norm(v[index(i)]);
    return s;
  }

template<class R>
double  KN_<R>::l2() const {
  double s = 0.;
    for (long i=0;i<n;i++)
      s +=  RNM::norm(v[index(i)]);
    return sqrt(s);
  }
template<class R>
double  KN_<R>::l1() const {
  double s = 0.;
    for (long i=0;i<n;i++)
      s +=  std::abs(v[index(i)]);
    return (s);
  }
template<class R>
double  KN_<R>::linfty() const {
  double s = 0.;
    for (long i=0;i<n;i++)
      s = std::max( (double) std::abs(v[index(i)]),s);
    return (s);
  }
template<class R>
double  KN_<R>::lp(double p) const {
  if( p==1.) return l1();
  else if (p==2.) return l2();
  else if(p>1.e10) return linfty();
  else
  {
  double s = 0.;
    for (long i=0;i<n;i++)
      s = pow(std::max( (double) std::abs(v[index(i)]),s),p);
    return pow(s,1./p);
   }
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
 void  KN_<R>::map(R (*f)(R )) {
    for (long i=0;i<n;i++)
      {  R & x(v[index(i)]);
          x =  f(x);}
   
  }

template<class R>
 void  KN_<R>::map(R (*f)(const R& )) {
    for (long i=0;i<n;i++)
      {  R & x(v[index(i)]);
          x =  f(x);}
   
  }
  
template<class R>
template<class T>
 void  KN_<R>::set(R (*f)(const  T& ),KN_<T> & u)
 { 
   K_throwassert(N() == u.N());
   for (long i=0;i<n;i++)
      {  R & x(v[index(i)]);
          v[index(i)]= f(u[i]);}
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
