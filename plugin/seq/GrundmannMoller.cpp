/* Authors: A. Chabib, P-H. Tournier & F. Hecht
This is an implementation of  Grundmann & Moller formula,
Article: Invariant Integration Formulas for the n-Simplex 
by Combinatorial Methods,SIAM Journal on Numerical Analysis,
V. 15, Iss 2 1978. 
This file generates three operators that are called from the DSL and 
take an order (integer) 's' as input. The exact quadrature formulas 
allow the exact integration of polynomials of degree up to (2s+1)
/!\ Warning: The method exhibits stability issues as the order increases
    the one can have some FreeFem warnings like : 
    Erreur dans la formule d'integration d=1 exact = 40 Nb Point = 861
    which can also affect the quality of the results

*/
#include <vector>
#include "AFunction.hpp"
#include "error.hpp"
#include "RNM.hpp"
#include "rgraph.hpp"
#include "fem.hpp"
#include "QuadratureFormular.hpp"

using namespace Fem2D;
using namespace std;

typedef long double REAL;
//using REAL = boost::multiprecision::cpp_dec_float_100;

// Compute the factorial
REAL factorial(int nbr) {
    if (nbr <= 1)
        return REAL(1.);
    else {
        REAL prod = REAL(1.);
        for (int i = 1; i <= nbr; i++) {
            prod *= i;
        }
        return prod;
    }
}

// Generates all compositions of the integer sum into nparts non-negative parts.
vector<vector<REAL>> generate_partitions(int sum, int nparts) {
    vector<vector<REAL>> partitions;

    if (nparts == 1) {
        partitions.push_back({(REAL)(sum)});
        return partitions;
    }

    for (int first = 0; first <= sum; first++) {
        vector<vector<REAL>> rest_partitions = generate_partitions(sum - first, nparts - 1);
        for (auto &rest : rest_partitions) {
            vector<REAL> comp;
            comp.push_back(first);
            comp.insert(comp.end(), rest.begin(), rest.end());
            partitions.push_back(comp);
        }
    }
    return partitions;
}

// Implementation of Grundmann & Moller formula,
// Article: Invariant Integration Formulas for the n-Simplex by Combinatorial Methods
//SIAM Journal on Numerical Analysis, V. 15, Iss 2 1978
vector<pair<REAL, vector<REAL>>> integration_weightspoints(int n, int s) {

    int d = 2 * s + 1;
    vector<pair<REAL, vector<REAL>>> data;
    REAL signe =REAL(1.);
    for (int i = 0; i < s + 1; i++) {
        REAL dn_2i=d+n-2*i;
        REAL coef = signe;
        int j_hi = max(n, max(d, d + n - i));
        // Loop to avoid the overflow for a better understanding check the formula 
        for(int j = 1; j <= j_hi; j++)
        {
        REAL rj= REAL(j);
        if(j <= n) { coef *= rj; }
        if(j <= d) { coef *= dn_2i; }
        if(j <= 2 * s) { coef /= REAL(2.0); }
        if(j <= i) { coef /= rj; }
        if(j <= d + n - i) { coef /= rj; }
        }
        signe = - signe;

        vector<vector<REAL>> betas = generate_partitions(s - i, n + 1);
        REAL denom = d + n - 2 * i;

        for (const auto &beta_orig : betas) {
            vector<REAL> beta;
            for (auto x : beta_orig) {
                beta.push_back((2 * x + 1) / denom);
            }
            data.push_back({coef, beta}); // weight, points for each beta
        }
    }

    return data;
}

// Static Quadrature formulas tables
static GQuadratureFormular<R1> **TabQuadrFormula1 = new GQuadratureFormular<R1> *[1000];
static GQuadratureFormular<R2> **TabQuadrFormula2 = new GQuadratureFormular<R2> *[1000];
static GQuadratureFormular<R3> **TabQuadrFormula3 = new GQuadratureFormular<R3> *[1000];

// Wrapper DSL

/*template<typename Ri>
const GQuadratureFormular<Ri> *GenerateQuadratureFormularForOperatorRi(Stack stack, const long &s) {
    // One Wrapper but requires  anew version > c++17, 
    // To solve the issue we generate 3 wrappers ( Dim 1,2& 3)  
    if  constexpr ( (std::is_same<Ri, R1>::value)){
        if( TabQuadrFormula1[s] != nullptr) 
        return TabQuadrFormula1[s];

    } 
    else if constexpr( (std::is_same<Ri, R2>::value) ) {
        if(TabQuadrFormula2[s] != nullptr)
           return TabQuadrFormula2[s];
    }
    else if constexpr( (std::is_same<Ri, R3>::value)){if  ( TabQuadrFormula3[s] != nullptr )  return TabQuadrFormula3[s];} 

    int n = std::is_same<Ri, R1>::value ? 1 : (std::is_same<Ri, R2>::value ? 2 : 3);

    auto Vec = integration_weightspoints(n, s);
    int NbrP = Vec.size();
    GQuadraturePoint<Ri> *Tab = new GQuadraturePoint<Ri>[NbrP];
    for (int i = 0; i < NbrP; i++) {
        if  constexpr (std::is_same<Ri, R1>::value) Tab[i] = GQuadraturePoint<R1>(R1(static_cast<long double>(Vec[i].second[0])), static_cast<long double>(Vec[i].first));
        else if constexpr (std::is_same<Ri, R2>::value) Tab[i] = GQuadraturePoint<R2>(R2(static_cast<long double>(Vec[i].second[0]), static_cast<long double>(Vec[i].second[1])), static_cast<long double>(Vec[i].first) );
        else if constexpr (std::is_same<Ri, R3>::value) Tab[i] = GQuadraturePoint<R3>(R3(static_cast<long double>(Vec[i].second[0]), static_cast<long double>(Vec[i].second[1]), static_cast<long double>(Vec[i].second[2])), 
        static_cast<long double>(Vec[i].first) );

    }
        
    if constexpr (std::is_same<Ri, R1>::value){
        TabQuadrFormula1[s] = new GQuadratureFormular<R1>(s, NbrP, Tab);
        return TabQuadrFormula1[s];
    }
    else if constexpr(std::is_same<Ri, R2>::value){
        TabQuadrFormula2[s] = new GQuadratureFormular<R2>(s, NbrP, Tab);
        return TabQuadrFormula2[s];
    } 
    else if constexpr(std::is_same<Ri, R3>::value){
        TabQuadrFormula3[s] = new GQuadratureFormular<R3>(s, NbrP, Tab);
        return TabQuadrFormula3[s];

    }
}*/

const GQuadratureFormular<R1> *GenerateQuadratureFormularForOperatorR1(Stack stack, const long &s) {
    if( TabQuadrFormula1[s] != nullptr) return TabQuadrFormula1[s];
    int n =1;
    auto Vec = integration_weightspoints(n, s);
    int NbrP = Vec.size();
    GQuadraturePoint<R1> *Tab = new GQuadraturePoint<R1>[NbrP];
    for (int i = 0; i < NbrP; i++) {
        Tab[i] = GQuadraturePoint<R1>(R1(static_cast<long double>(Vec[i].second[0])), static_cast<long double>(Vec[i].first));
    }
    TabQuadrFormula1[s] = new GQuadratureFormular<R1>(s, NbrP, Tab);
    return TabQuadrFormula1[s];
}

const GQuadratureFormular<R2> *GenerateQuadratureFormularForOperatorR2(Stack stack, const long &s) {
    if(TabQuadrFormula2[s] != nullptr) return TabQuadrFormula2[s];
    int n =2;
    auto Vec = integration_weightspoints(n, s);
    int NbrP = Vec.size();
    GQuadraturePoint<R2> *Tab = new GQuadraturePoint<R2>[NbrP];
    for (int i = 0; i < NbrP; i++) {
        Tab[i] = GQuadraturePoint<R2>(R2(static_cast<long double>(Vec[i].second[0]), static_cast<long double>(Vec[i].second[1])), static_cast<long double>(Vec[i].first) );
    }
        
    TabQuadrFormula2[s] = new GQuadratureFormular<R2>(s, NbrP, Tab);
    return TabQuadrFormula2[s];
}


const GQuadratureFormular<R3> *GenerateQuadratureFormularForOperatorR3(Stack stack, const long &s) {
    if  ( TabQuadrFormula3[s] != nullptr )  return TabQuadrFormula3[s];
    int n =3;
    auto Vec = integration_weightspoints(n, s);
    int NbrP = Vec.size();
    GQuadraturePoint<R3> *Tab = new GQuadraturePoint<R3>[NbrP];
    for (int i = 0; i < NbrP; i++) {
        Tab[i] = GQuadraturePoint<R3>(R3(static_cast<long double>(Vec[i].second[0]), static_cast<long double>(Vec[i].second[1]), static_cast<long double>(Vec[i].second[2])), 
        static_cast<long double>(Vec[i].first) );

    }
    TabQuadrFormula3[s] = new GQuadratureFormular<R3>(s, NbrP, Tab);
        return TabQuadrFormula3[s];

}

#include "lex.hpp"
extern mylex *zzzfff;
static void Load_Init() {
    // Initialization 
    for (int i = 0; i < 1000; i++){
        TabQuadrFormula1[i] = nullptr;
        TabQuadrFormula2[i] = nullptr;
        TabQuadrFormula3[i] = nullptr;

    }

    // Declare operators to be used in the DSL 
    Global.Add("GMQuadrature1D", "(", new OneOperator1s_<const GQuadratureFormular<R1> *, long>(GenerateQuadratureFormularForOperatorR1));
    Global.Add("GMQuadrature2D", "(", new OneOperator1s_<const GQuadratureFormular<R2> *, long>(GenerateQuadratureFormularForOperatorR2));
    Global.Add("GMQuadrature3D", "(", new OneOperator1s_<const GQuadratureFormular<R3> *, long>(GenerateQuadratureFormularForOperatorR3));

}

LOADFUNC(Load_Init)