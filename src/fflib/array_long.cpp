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
const basicForEachType *aatypeknlongp;
/*
void initArrayOperators()
{
     ArrayOperator<double>();
     ArrayOperator<Complex>();
     ArrayOperator<long>();

}
void  initArrayDCL()
{
    ArrayDCL<double>();
    ArrayDCL<Complex>();
    ArrayDCL<long>();
}
*/

 aType aaaa_knlp;
void initArrayDCLlong()
{
//     ArrayOperator<long>();
     Dcl_Type<Inv_KN_long>(); // Add FH mars 2005
     ArrayDCL<long>();
     aaaa_knlp = atype<KN<long>*>();

}


class OneBinaryOperatorInv_KN_long : public OneOperator { public:  
  OneBinaryOperatorInv_KN_long(basicForEachType * ti) : OneOperator(atype<Inv_KN_long >(), ti ,atype<long>()) {}
    E_F0 * code(const basicAC_F0 & args) const 
     { Expression p=args[1];
       if ( ! p->EvaluableWithOutStack() ) 
        { 
          bool bb=p->EvaluableWithOutStack();
          cout << bb << " " <<  * p <<  endl;
          CompileError("Inverse:  int[int] I,  array, with  I^p, The p must be a constant == -1, sorry");}
       long pv = GetAny<long>((*p)(0));
        if (pv !=-1)   
         { char buf[100];
           sprintf(buf,"Inverse:  int[int] I,  array, I^%ld, The pow must be  == -1, sorry",pv);
           CompileError(buf);}     
       return  new E_F_F0<Inv_KN_long,KN_<long> >(Build<Inv_KN_long,KN_<long> >,to< KN_<long> >(args[0])); 
    }
};


void initArrayOperatorlong()
{
    typedef long K;
     ArrayOperator<long,long>();
     // to def inverse permutation // Add FH mars 2005
     TheOperators->Add("^", new OneBinaryOperatorInv_KN_long(atype<KN_<long> >() ));
    //- TheOperators->Add("^", new OneBinaryOperatorInv_KN_long(atype<KN<long> *>() )) ;
     aatypeknlongp= atype<KN<long>*>(); // for  compilation error with g++ 3.2.2

     Add<KN_<long> >("sort",".",new OneOperator1_<KN_<K>,KN_<K> >(SortKn<K, KN_<K> >));
    // Add<KN<long> >("sort",".",new OneOperator1_<KN<K>,KN<K> >(SortKn<K, KN<K> >));
   //  Add<KN<long> *>("sort",".",new OneOperator1_<KN<K>*,KN<K>* >(SortpKn<K>));
     
     
//     ArrayDCL<long>();
}

 void xxxx() {
}
