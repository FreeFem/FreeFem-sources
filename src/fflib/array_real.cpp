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


double square(double x){return x*x;}

void initArrayDCLdouble()
{
//     ArrayOperator<long>();
     ArrayDCL<double>();
}
void initArrayOperatordouble()
{

    ArrayOperator<double,long>();
    ArrayOperatorF<double,double>();
    typedef double K;
    typedef double KK;
    
     Global.Add("abs","(",new OneOperator1F_KN_<F_KN_<K,K,KK>,K,KK,KN_<K> >(fabs));
     Global.Add("acos","(",new OneOperator1F_KN_<F_KN_<K,K,KK>,K,KK,KN_<K> >(acos));
     Global.Add("asin","(",new OneOperator1F_KN_<F_KN_<K,K,KK>,K,KK,KN_<K> >(asin));
     Global.Add("atan","(",new OneOperator1F_KN_<F_KN_<K,K,KK>,K,KK,KN_<K> >(atan));
     Global.Add("floor","(",new OneOperator1F_KN_<F_KN_<K,K,KK>,K,KK,KN_<K> >(floor));
     Global.Add("ceil","(",new OneOperator1F_KN_<F_KN_<K,K,KK>,K,KK,KN_<K> >(ceil));

    Global.Add("square","(",new OneOperator1F_KN_<F_KN_<K,K,KK>,K,KK,KN_<K> >(square));
    
//     ArrayDCL<long>();
}
