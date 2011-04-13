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

template<class R>  
KNM_<R> & KNM_<R>::operator oper (const outProduct_KN_<R> & u)  
{
  //   *this  oper  A* t B 
    K_throwassert (shapei.SameShape(u.a) && shapej.SameShape(u.b) );
    long n= N(), m= M();
    
    R * ai(u.a),cc, c= u.c;
    long stepi=u.a.step;
    R * bj, *bb(u.b);
    long stepj=u.b.step;
    KN_<const_R>  li((*this)(0,'.')); //  first line
    int stepij= li.step;
    for (long i=0;i<n;i++,ai += stepi,++li)
      {
        cc= c * *ai;
        R * mij = li;
        bj = bb;
        for (long j=0;   j<m; j++, bj += stepj, mij += stepij )         
          *mij oper cc * conj(*bj) ; 
       }
    return *this;
 }


template<class R>
template<class  A,class B,class C>
 KN_<R>& KN_<R>::operator oper (const F_KN_<A,B,C> & u)  {
    K_throwassert ( u.check(this->N()) );
    R * l(v);  //  first line   
    for (long i=0;i<n;i++,l += step)  
      *l oper  u[i]; 
    return *this;}

template<class R>
KN_<R>& KN_<R>::operator oper (const SetArray<R> & u)  {
    R * l(v);  //  first line   
    for (long i=0;i<n;i++,l += step)  
	*l oper  u[i]; 
return *this;}



template<class R>
 KN_<R>& KN_<R>::operator oper (const Mul_KNM_KN_<R> & u)  {
    K_throwassert (SameShape(u.A.shapei) && !constant());
    R * l(v); KN_<const_R>  li(u.A(0,'.')); //  first line   
    for (long i=0;i<n;i++,l += step,++li)  
      *l oper (li,u.b); 
    return *this;}


template<class R>
 KN_<R>&  KN_<R>::operator oper (const DotStar_KN_<R> & u) {
    K_throwassert(u.a.N() == N()  );
    long stepa(u.a.step),stepb(u.b.step);
    R * l(v); const_R  *aa(u.a), *bb(u.b);    
    for (long i=0;i<n;i++,l += step, aa +=stepa, bb += stepb)
      *l oper *aa * *bb;
    return *this;
  }
template<class R>
 KN_<R>&  KN_<R>::operator oper (const DotSlash_KN_<R> & u) {
    K_throwassert(u.a.N() == N()  );
    long stepa(u.a.step),stepb(u.b.step);
    R * l(v); const_R  *aa(u.a), *bb(u.b);    
    for (long i=0;i<n;i++,l += step, aa +=stepa, bb += stepb)
      *l oper *aa / *bb;
    return *this;
  }

  
template<class R>
 KN_<R>&  KN_<R>::operator oper (const Add_KN_<R> & u) {
    K_throwassert(u.a.N() == N()  );
    long stepa(u.a.step),stepb(u.b.step);
    R * l(v); const_R  *aa(u.a), *bb(u.b);    
    for (long i=0;i<n;i++,l += step, aa +=stepa, bb += stepb)
      *l oper *aa+*bb;
    return *this;
  }

template<class R>
 KN_<R>&  KN_<R>::operator oper (const Sub_KN_<R> & u) {
    K_throwassert(u.a.N() == N()  );
    long stepa(u.a.step),stepb(u.b.step);
    R * l(v); const_R  *aa(u.a), *bb(u.b);    
    for (long i=0;i<n;i++,l += step, aa +=stepa, bb += stepb)
      *l oper  *aa-*bb;
    return *this;
  }
  
template<class R>
 KN_<R>&  KN_<R>::operator oper (const Mulc_KN_<R> & u) {
    K_throwassert(u.a.N() == N()  );
    long stepa(u.a.step);
    R * l(v); const_R  *aa(u.a),bb(u.b)  ;
    for (long i=0;i<n;i++,l += step, aa +=stepa)
      *l oper *aa * bb;
    return *this;
  }
  
template<class R>
 KN_<R>&  KN_<R>::operator oper (const Add_Mulc_KN_<R> & u) {
    K_throwassert(u.a.N() == N()  );
    const long stepa(u.a.step),stepb(u.b.step);
    const R ca(u.ca),cb(u.cb);    
    R * l(v);
    const R *aa(u.a),*bb(u.b);    
    for (long i=0;i<n;i++,l += step, aa +=stepa, bb += stepb)
      *l oper *aa*ca + *bb*cb;
    return *this;
  }

template<class R>
 KN_<R>&  KN_<R>::operator oper (const if_arth_KN_<R> & u) {
    K_throwassert(u.a.N() == N()  );
    R zero=R();
    const long stepa(u.a.step),stepb(u.b.step),stepc(u.c.step);
    R * l(v);
    const R *aa(u.a),*bb(u.b),*cc(u.c);    
    for (long i=0;i<n;i++,l += step, aa +=stepa, bb += stepb ,  cc += stepc)
      *l oper ( (*aa != zero) ?  *bb : *cc);
    return *this;
  }
