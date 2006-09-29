// -*- Mode : c++ -*-
//
// SUMMARY  :      
// USAGE    :        
// ORG      : 
// AUTHOR   : Frederic Hecht
// E-MAIL   : hecht@ann.jussieu.fr
//

/*
 
 This file is part of Freefem++
 
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

#include "array_tlp.hpp"
#include "array_init.hpp"


KN_<double> Get_C2R(KN_<complex<double> >  vc)
{
  ffassert(vc.step==1);
  complex<double> * pc=vc; // pointeur du tableau                                                                            
  double *pr = static_cast<double*>(static_cast<void*>(pc));
  return KN_<double>(pr,vc.N()*2);
}

template<int offset>
KN_<double> Get_C2R_(KN_<complex<double> >  vc)
{
  complex<double> * pc=vc; // pointeur du tableau                                                                            
  double *pr = static_cast<double*>(static_cast<void*>(pc));
  return KN_<double>(pr+offset,vc.N(),vc.step*2);
}


void initArrayDCLComplex()
{
     ArrayDCL<Complex>();
}

Complex square(const Complex & x){return x*x;}

void initArrayOperatorComplex()
{
    typedef Complex K;
    typedef const K & KK;
    ArrayOperator<Complex,long>();
    ArrayOperatorF<Complex,const Complex&>();
     // take the real or imag part of complex array
    Add<KN<Complex> *>("im",".",new OneOperator1<KN_<double>,KN_<Complex> >(Get_C2R_<1>,atype<KN<Complex> * >()));
    Add<KN_<Complex> >("im",".",new OneOperator1<KN_<double>,KN_<Complex> >(Get_C2R_<1>,atype<KN_<Complex>  >()));
    Add<KN<Complex> *>("re",".",new OneOperator1<KN_<double>,KN_<Complex> >(Get_C2R_<0>,atype<KN<Complex> * >()));
    Add<KN_<Complex> >("re",".",new OneOperator1<KN_<double>,KN_<Complex> >(Get_C2R_<0>,atype<KN_<Complex>  >()));

    Global.Add("square","(",new OneOperator1F_KN_<F_KN_<K,K,KK>,K,KK,KN_<K> >(square));
    Global.Add("conj","(",new OneOperator1F_KN_<F_KN_<K,K,KK>,K,KK,KN_<K> >(conj));


}
