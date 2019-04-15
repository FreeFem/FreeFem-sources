
#ifndef SPLITSIMPLEX_HPP_
#define SPLITSIMPLEX_HPP_

template<class Rd>  void  SplitSimplex(int N,int & nv, Rd *& P, int & nk , int *& K);
// Add J. Morice for function trunc.
void SplitSurfaceSimplex(int N,int &ntri2, int *&tri);
void SplitEdgeSimplex(int N,int &ntri2, int *&tri);

#endif //SPLITSIMPLEX_HPP_
