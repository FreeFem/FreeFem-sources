#include "array_tlp.hpp"
#include "array_init.hpp"

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

void initArrayDCLlong()
{
//     ArrayOperator<long>();
     Dcl_Type<Inv_KN_long>(); // Add FH mars 2005
     ArrayDCL<long>();

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
           sprintf(buf,"Inverse:  int[int] I,  array, I^%d, The pow must be  == -1, sorry",pv);
           CompileError(buf);}     
       return  new E_F_F0<Inv_KN_long,KN_<long> >(Build<Inv_KN_long,KN_<long> >,to< KN_<long> >(args[0])); 
    }
};


void initArrayOperatorlong()
{
     
     ArrayOperator<long>();
     // to def inverse permutation // Add FH mars 2005
     TheOperators->Add("^", new OneBinaryOperatorInv_KN_long(atype<KN_<long> >() ));
     TheOperators->Add("^", new OneBinaryOperatorInv_KN_long(atype<KN<long> *>() )) ;

     
     
//     ArrayDCL<long>();
}
