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
     ArrayOperator<double>();
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
