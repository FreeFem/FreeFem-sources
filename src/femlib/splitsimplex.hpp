
#ifndef SPLITSIMPLEX_HPP_
#define SPLITSIMPLEX_HPP_
#include "throwassert.hpp"
template<class Rd>  void  SplitSimplex(int N,int & nv, Rd *& P, int & nk , int *& K);
// Add J. Morice for function trunc.
void SplitSurfaceSimplex(int N,int &ntri2, int *&tri);
void SplitEdgeSimplex(int N,int &ntri2, int *&tri);
// add for generation of Pkdc Finiet element mars 2021 FH.
void  invNumSimplex2(int n,int &i,int &j);
void  invNumSimplex3(int n,int &i,int &j,int &k);

#endif //SPLITSIMPLEX_HPP_
