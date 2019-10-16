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

// #ifndef RNM_OPC_HPP_	//NO: multiple call
// #define RNM_OPC_HPP_

template<class R>
inline KN_<R>&  KN_<R>::operator oper (const_R a)  {
    R * l(v);
    for (long i=0;i<n;i++,l += step)
      *l oper a;
    return *this;
  }

template<class R>
inline    KNM_<R> & KNM_<R>::operator oper (const_R a)
{
  if(IsVector1() )
        KN_<R>::operator oper (a);
  else {
          KN_<R>  lj(operator()('.',0)); //  (.,.,O)
          for (long  j=0;j<M();++j,++lj)
             lj oper a;}
  return *this;
}

template<class R>
inline    KNMK_<R> & KNMK_<R>::operator oper (const_R a)
{
  if(IsVector1() )
        KN_<R>::operator oper (a);
  else {
          KNM_<R>  lj(operator()('.','.',0)); //  (.,.,O)
          long j=K();
           while(j--)
            {lj oper a;++lj;}
       }
  return *this;
}

template<class R>
inline KN_<R>&  KN_<R>::operator oper (const KN_<const_R> & u)   {
    K_throwassert(!u.step || u.n == n);// Add constant vector
    R * l(v);
    const R *r(u);
    for (long i=0;i<n;i++,l += step, r += u.step) *l oper *r;
    return *this;
  }

template<class R>
inline    KNM_<R> & KNM_<R>::operator oper (const KNM_<const_R> & u)
{
  K_throwassert( N() == u.N() && M() == u.M() );
  if(IsVector1() && u.IsVector1() && shapei.step == u.shapei.step ) // modif 2011 (thank to Oka)
        KN_<R>::operator oper(u); // modif FH jan 2004
  else {
      long n=N(),m=M();
      for(int j=0; j <m;++j)
          for(int i=0; i <n;++i)
            this->operator()(i,j) oper u(i,j);
        }
  return *this;
}


template<class R>
inline  KNMK_<R> & KNMK_<R>::operator oper (const KNMK_<const_R> & u)
{
  K_throwassert( N() == u.N() && M() == u.M() &&   K() == u.K() );

  if(IsVector1() && u.IsVector1() && u.N() == N() &&  shapei.step == u.shapei.step) // modif 2011 (thank to Oka)
        KN_<R>::operator oper(u); // modif FH 2004
  else {
          K_throwassert( K() == u.K());
          KNM_<R>  lj(operator()('.','.',0)); //  (.,O)
          KNM_<const_R>  uj(u('.','.',0));
          long j=K();
          while (j--)
           { lj oper uj;++lj;++uj;}
       }
  return *this;
}

#undef oper

// #endif //RNM_OPC_HPP_
