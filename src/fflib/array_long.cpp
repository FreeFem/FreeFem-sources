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
// end Hack 
// Add mars 2010 
template<class R>  R * set_init_init( R* const & a,const long & n){ 
    SHOWVERB( cout << " set_init " << typeid(R).name() << " " << n << endl);
    a->init(n);
    for (int i=0;i<n;i++)
	(*a)[i].init();
    return a;}
inline   string ** get_elements( KN<String> *  const  &  a,long  const   & b)
{   ffassert( a && b >=0 && b < a->size());
    String & Sret =  (*a)[b]; // correction FH feb 2004
    // delete b; la chaine est detruire automatiquement en fin d'instruction  FH jan 2010
    return Sret.getap();}

template<class A> inline AnyType Destroy_KN(Stack,const AnyType &x){
    KN<A> * a=GetAny<KN<A>*>(x);
    for (int i=0;i<a->N(); i++)
	(a)[i].destroy();
    a->destroy(); 
    return  Nothing;
}
// fin add


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
    
    Dcl_TypeandPtr_<KN_<String>,KN<String> *>(0,0,0,::Destroy<KN<String> >, ::ClearReturnKK_<K,KN<String>,KN_<String> >,::ClearReturnpKK<String,KN<String> >);
    atype<KN<String>* >()->Add("[","",new OneOperator2_<string**,KN<String>*,long >(get_elements));
    TheOperators->Add("<-", 
		      new OneOperator2_<KN<String> *,KN<String> *,long>(&set_init_init));
    map_type_of_map[make_pair(atype<long>(),atype<string*>())]=atype<KN<String>*>(); // vector
    Add<KN<String> *>("n",".",new OneOperator1<long,KN<String> *>(get_n));  
    extern   KN<String> *pkarg;
    Global.New("ARGV",CPValue<KN<String> >(*pkarg));// add FH mars 2010
    
}

 void xxxx() {
}
