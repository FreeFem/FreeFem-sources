// used by splitsimplex.cpp 
// ----   build a the simplex decomposition 
/*
template void  SplitSimplex<R1>(int N,int & nv, R1 *& P, int & nk , int *& K);
template void  SplitSimplex<R2>(int N,int & nv, R2 *& P, int & nk , int *& K);
template void  SplitSimplex<R3>(int N,int & nv, R3 *& P, int & nk , int *& K);
$:
*/
//
//   
// ORIG-DATE:     fev 2009
// -*- Mode : c++ -*-
//
// SUMMARY  :  Model  mesh 2d   
// USAGE    : LGPL      
// ORG      : LJLL Universite Pierre et Marie Curie, Paris,  FRANCE 
// AUTHOR   : Frederic Hecht
// E-MAIL   : frederic.hecht@ann.jussieu.fr
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

 Thank to the ARN ()  FF2A3 grant
 ref:ANR-07-CIS7-002-01  */

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <cassert>
#include <fstream>
#include "ufunction.hpp"
using namespace std;
namespace Fem2D  {
#include "R3.hpp"
};
using Fem2D::R1;
using Fem2D::R2;
using Fem2D::R3;
extern  long verbosity ;

#include "splitsimplex.hpp"

/*
  construction an array of sub simplex for plot ... 
  see last template function 
     SplitSimplex(N,nv,P,nk,K);
  N > 1 -> classical split un N^d sub simplex
  nv : number of point on    
*/
//  d = 1   trivial
void SplitSimplex(int N,R1 *P,int *K,int op=0,R1 *AB=0)
{
  assert(N>0);
  double h = 1./ N;
  for(int i=0;i<=N;++i)
    if(AB)
      P[i+op] = R1(i*h).Bary(AB);
    else
      P[i+op] = R1(i*h);
    int l=0;
  for(int i=0;i<N;++i) {
    K[l++]=i+op;
    K[l++]=i+1+op;
   if(verbosity>200)
     cout << "l="<< l/2 <<  " "<< K[l-2] <<" "<<  K[l-1] <<" "<<  endl;
    }
}
/*
0124 1234 1345 3456 2346 3567 1 solutions de taille 6  trouvees.
0124 1245 1235 2356 2456 3567 Solution courte numero 7:
0124 1246 1236 1356 1456 3567 Solution courte numero 12:
*/



typedef int int4 [4] ;

// InvIntFunc

inline int NumSimplex2(int i) {return ((i)*(i+1))/2;}
inline int NumSimplex2_1(int l) { return sqrt(l*2)+3;}
inline int NumSimplex3(int i) {return ((i)*(i+1)*(i+2))/6;}
inline int NumSimplex3_1(int l) { return pow(l*6,1./3)+4;}

#define InvIntFunction  invNumSimplex2
#define F(i) NumSimplex2(i)
#define F_1(i) NumSimplex2_1(i)
#include "InvIntFunc.cpp"
#undef F
#define F_1
#undef InvIntFunction

#define InvIntFunction invNumSimplex3
#define F(i) NumSimplex3(i)
#define F_1(i) NumSimplex3_1(i)
#include "InvIntFunc.cpp"
#undef F
#define F_1
#undef InvIntFunction

inline int NumSimplex1(int i) { return i;}
inline int NumSimplex2(int i,int j) { return j+NumSimplex2(i+j);}
inline int NumSimplex3(int i,int j,int k) { return NumSimplex3(i+j+k)+NumSimplex2(j+k)+k;}


inline void  invNumSimplex2(int n,int &i,int &j) 
{
  int k= invNumSimplex2(n); //( i+j) 
  j=n-NumSimplex2(k);
  i= k-j;
  //  cout << n << " " << k << " -> " << i << " " << j << endl;
  assert( n == NumSimplex2(i,j));
}
inline void  invNumSimplex3(int n,int &i,int &j,int &k) 
{
  int l= invNumSimplex3(n); //=( i+j+k) 
  invNumSimplex2(n-NumSimplex3(l),j,k);
  assert(j>=0 && k>=0);
  i=l-k-j; 
  // cout << n << "   " << l << "-> " << i << " " << j << " " << k <<endl;
  assert( n == NumSimplex3(i,j,k)) ;
}


// d = 2
void SplitSimplex(int N,R2 *P,int *K,int op=0,R2 *ABC=0)
{
    // warning PB of oriention if (i+j<N)  , the transo is positif => no swap on 1, 2 ve
  assert(N>0);
  int nv = (N+1)*(N+2)/2;
  double h=1./N;
  //   loop sur les diag   i+j = k 
  //   num  ( i+j,j) lexico croissant  
  for(int l=0;l<nv;l++)
    {
      int i,j;
      invNumSimplex2(l,i,j);
      if(ABC)
	P[l+op]= R2(i*h,j*h).Bary(ABC);
      else 
	P[l+op]= R2(i*h,j*h);
      assert(l<nv);
    }
  //    generation des trianges
  // --------
  int l=0;
  for (int i=0;i<N;++i)
    for (int j=0;j<N;++j)
      if(i+j<N) 
	{
	  K[l++]= op+NumSimplex2(i,j);
	  K[l++]= op+NumSimplex2(i+1,j);
	  K[l++]= op+NumSimplex2(i,j+1);
	}
      else
	{
	  K[l++]= op+NumSimplex2(N-i,N-j);
      K[l++]= op+NumSimplex2(N-i-1,N-j);//  swap vertices  1 and 2 May 2019 FH (Thanks to AF)
	  K[l++]= op+NumSimplex2(N-i,N-j-1);
	}
}

/*
// d = 3 Surfacic
void SplitSimplex(int N,R3 *P,int *K,int op=0,R3 *ABC=0)
{
    assert(N>0);
    int nv = (N+1)*(N+2)/2;
    double h=1./N;
    //   loop sur les diag   i+j = k
    //   num  ( i+j,j) lexico croissant
    for(int l=0;l<nv;l++)
    {
        int i,j,k;
        invNumSimplex3(l,i,j,k);
        if(ABC)
            P[l+op]= R3(i*h,j*h,k*h).Bary(ABC);
        else
            P[l+op]= R3(i*h,j*h,k*h);
        assert(l<nv);
    }
    //    generation des trianges
    // --------
    int l=0;
    for (int i=0;i<N;++i)
        for (int j=0;j<N;++j)
            if(i+j<N)
            {
                K[l++]= op+NumSimplex2(i,j);
                K[l++]= op+NumSimplex2(i+1,j);
                K[l++]= op+NumSimplex2(i,j+1);
            }
            else
            {
                K[l++]= op+NumSimplex2(N-i,N-j);
                K[l++]= op+NumSimplex2(N-i,N-j-1);
                K[l++]= op+NumSimplex2(N-i-1,N-j);
            }
} */


void SplitSimplex(int N,R3 *P,int *tet,int op=0,R3* Khat=0)
{
    const int n=N;
  const int n2=n*n;
  const int n3=n2*n;
  const int ntc=6;
  int nv = (N+1)*(N+2)*(N+3)/6;
  
  int d1[6][4] = { {0,1,2,4} , {1,2,4,3},{1,3,4,5},{3,4,5,6},{ 2,3,6,4}, {3,5,7,6} };
  int ntt=0;
  int n8[8];
  int ptet=0;
  double h=1./N;
  for(int l=0;l<nv;l++)
    {
      int i,j,k;
      invNumSimplex3(l,i,j,k);
     assert( i>=0&& j >=0  && k >= 0);
     assert( i<=N && j <= N && k <= N);
      if(Khat)
	P[l+op]= R3(i*h,j*h,k*h).Bary(Khat);
      else 
	P[l+op]= R3(i*h,j*h,k*h);
      assert(l<nv);
    }
  // comment n3  i=m[n]
  //  m = i+j*n+k*n2 
  for (int m=0;m<n3;m++)
    {
      int i = m % n;
      int j = (m / n) % n;
      int k = m /( n2);
      for (int l=0;l<8;++l)
	{
	  int ii = ( i + ( (l & 1) != 0) );
	  int jj = ( j + ( (l & 2) != 0) );
	  int kk = ( k + ( (l & 4) != 0) );
	  int ll= NumSimplex3(ii,jj,kk);
	  n8[l]=ll;
	  if(ii+jj+kk>n) n8[l]=-1;
	}
      for (int l=0;l<ntc;++l)
	if(ntt<n3)
	  {
	    int out =0;
	    for(int m=0;m<4;++m)
	      if( (tet[ptet++]= n8[d1[l][m]]+op) < op) out++;
	    if(out == 0 )  ntt++;
	    else ptet -= 4; // remove tet
	  }
      
    }
 /* if(verbosity>199)
    {
      
      cout <<   "  SplitSimplex   " << endl;
      for (int i=0,l=0;i<n3;i++)
       for(int m=0;m<4;++m)
         cout << tet[l++] << (m==3 ? '\n' : ' ' );
       cout << ptet << "   " << tet << endl;
    }*/
  assert(ntt==n3);
}

// Add J. Morice (trunc functions)

void SplitSurfaceSimplex(int N,int &ntri2,int *&tri)
{
  const int n=N;
  const int n2=n*n;
  
  int ntri=3*ntri2;
  int op=0;
  
  tri = new int[ntri];
  //    generation des trianges
  // --------
  
  // face i=0
  int l=0;
  if(verbosity>200)
    cout << "face i=0" << endl;
    for (int i=0;i<N;++i)
        for (int j=0;j<N;++j){
            if(i+j<N)
	{
	  tri[l++]= op+NumSimplex3(0,i,j);
	  tri[l++]= op+NumSimplex3(0,i+1,j);
	  tri[l++]= op+NumSimplex3(0,i,j+1); 
	}
      else
	{ // correction orinatation ... 
	  tri[l++]= op+NumSimplex3(0,N-i,N-j);
	  tri[l++]= op+NumSimplex3(0,N-i-1,N-j);
          tri[l++]= op+NumSimplex3(0,N-i,N-j-1);
	}
      //cout << "i,j " << i << "," << j << endl;
      if(verbosity>200)
      cout << "l="<< l/3 <<  " "<< tri[l-3] <<" "<< tri[l-2] <<" "<<  tri[l-1] <<" "<<  endl;
    }
  // face j=0
  if(verbosity>200)
  cout << "face j=0" << endl;
  for (int i=0;i<N;++i)
    for (int j=0;j<N;++j){
      if(i+j<N) 
	{
	  tri[l++]= op+NumSimplex3(i,0,j);
	  tri[l++]= op+NumSimplex3(i,0,j+1); // inverser les deux lignes
	  tri[l++]= op+NumSimplex3(i+1,0,j);
	}
      else
	{ // CHANGE ORIENTATION FH MAY
	  tri[l++]= op+NumSimplex3(N-i,0,N-j);
	  tri[l++]= op+NumSimplex3(N-i,0,N-j-1);
          tri[l++]= op+NumSimplex3(N-i-1,0,N-j); // inverser les deux lignes
	}
      //cout << "i,j " << i << "," << j << endl;
      if(verbosity>200)
	cout << "l="<< l/3 <<  " "<< tri[l-3] <<" "<< tri[l-2] <<" "<<  tri[l-1] <<" "<<  endl;
    }
  // face k=0
  if(verbosity>200)
    cout << "face k=0" << endl;
  for (int i=0;i<N;++i)
    for (int j=0;j<N;++j){
      if(i+j<N) 
	{
	  tri[l++]= op+NumSimplex3(i,j,0);
	  tri[l++]= op+NumSimplex3(i+1,j,0);
	  tri[l++]= op+NumSimplex3(i,j+1,0);
	  //tri[l++]= op+NumSimplex3(i+1,j,0);
	}
      else
	{ // CHANGE ORIENTATION FH MAY
	  tri[l++]= op+NumSimplex3(N-i,N-j,0);
	  tri[l++]= op+NumSimplex3(N-i-1,N-j,0);
            tri[l++]= op+NumSimplex3(N-i,N-j-1,0);
	  //tri[l++]= op+NumSimplex3(N-i,N-j-1,0);
	}
      //cout << "i,j " << i << "," << j << endl;
      if(verbosity>200)
	cout << "l="<< l/3 <<  " "<< tri[l-3] <<" "<< tri[l-2] <<" "<<  tri[l-1] <<" "<<  endl;
    }
  // face i+j+k=1
  if(verbosity>200)
    cout << "dernier face " << endl;
  for (int k=0;k<N;++k)
    for (int j=0;j<N;++j){
      if(k+j<N) 
	{
	  int i=N-j-k;
	  tri[l++]= op+NumSimplex3(   i,   j,   k);
	  tri[l++]= op+NumSimplex3( i-1,   j, k+1); 
	  tri[l++]= op+NumSimplex3( i-1, j+1,   k);
	}
      else
	{
	  int i=N-(N-j-1)-(N-k);
	  tri[l++]= op+NumSimplex3(   i, N-j-1, N-k);
	  tri[l++]= op+NumSimplex3( i-1, N-j,   N-k); 
	  tri[l++]= op+NumSimplex3(   i, N-j, N-k-1);
	}
      if(verbosity>200)
	cout << "l="<< l/3 <<  " "<< tri[l-3] <<" "<< tri[l-2] <<" "<<  tri[l-1] <<" "<<  endl;
    }
  if(verbosity>200)
    cout << "l= " << l << " ntri=" << ntri << endl;
  assert( l == ntri);
}




void SplitEdgeSimplex(int N,int &nedge2,int *&edge)
{
  const int n=N;
  const int n2=n*n;
    
  int nedge=2*nedge2;
  int op=0;
    
  edge = new int[nedge];
  //    generation des edges
  // --------
  // (i,j,k) barycentric coordinates

  int l=0;

  for (int i=0;i<N;++i)
    for (int j=0;j<N;++j) {
      if(i+j<N) {
        edge[l++]= op+NumSimplex2(i+1,j);
        edge[l++]= op+NumSimplex2(i,j+1);
      }
      if(verbosity>200)
        cout << "l="<< l/2 <<" "<< edge[l-2] <<" "<<  edge[l-1] <<" "<<  endl;
    }
 
    for (int i=0;i<N;++i)
     for (int j=0;j<N;++j) {
       if(i+j<N) {
         edge[l++]= op+NumSimplex2(i,j);
         edge[l++]= op+NumSimplex2(i,j+1);
       }
       if(verbosity>200)
         cout << "l="<< l/2 <<" "<< edge[l-2] <<" "<<  edge[l-1] <<" "<<  endl;
     }
    
    for (int i=0;i<N;++i)
      for (int j=0;j<N;++j) {
        if(i+j<N) {
          edge[l++]= op+NumSimplex2(i,j);
          edge[l++]= op+NumSimplex2(i+1,j);
        }
        if(verbosity>200)
          cout << i+j << "l="<< l/2 <<" "<< edge[l-2] <<" "<<  edge[l-1] <<" "<<  endl;
      }

    if(verbosity>200)
      cout << "l= " << l << " nedge=" << nedge << endl;
    assert( l == nedge);
}



/*
void  SplitSimplex(int N,int & nv, R1 *& P, int & nk , int *& K)
{
  typedef R1 Rd;
  const int d = Rd::d;
  int cas = (N>0) ? 1 : d+1;
  N=abs(N);
  assert(N);
  int nv1=(N+1);
  int nk1= N;
  nv = cas*nv1;
  nk = cas*nk1;
  
  P = new Rd[nv];
  K = new int [nk*(d+1)];
  if( cas ==1) 
    SplitSimplex( N, P,K) ;
    else 
      {
	Rd AB1[2]= { Rd(0),Rd(0.5)};
	SplitSimplex( N, P,K,0,AB1)      ;
	Rd AB2[2]= { Rd(0.5),Rd(1)};
	SplitSimplex( N, P,K+nk1,nv1,AB2);      
      }
}


void  SplitSimplex(int N,int & nv, R2 *& P, int & nk , int *& K)
{
  typedef R2 Rd;
  const int d = Rd::d;
  int cas = (N>0) ? 1 : d+1;
  assert(N);
  N=abs(N);
  int nv1=N*(N+1)/2;
  int nk1= N*N;
  nv = cas*nv1;
  nk = cas*nk1;
  P = new Rd[nv];
  K = new int [nk*(d+1)];
  if( cas ==1) 
    SplitSimplex( N, P,K);
  else 
    { 
      Rd G=Rd::diag(1./(d+1));
      R2 Khat[d+1];
      for (int i=0;i<d;++i)
	Khat[i+1][i]=1; //  modit  25/2/2009
      for(int i=0;i<=d;++i)
	{
	  Rd S=Khat[i];
	  Khat[i]=G;
	  SplitSimplex( N, P,K+nk1*i,nv1*i,Khat);
	  Khat[i]=S;
	}     
    }
}

*/

template<class Rd>
void  SplitSimplex(int N,int & nv, Rd *& P, int & nk , int *& K)
{
  const int d = Rd::d,d1 =d+1;
  int cas = (N>0) ? 1 : d+1;
  assert(N);
  N=abs(N);
  int nv1=N+1; // nb simplexe
  int nk1=N; // nb vertices
  for(int i=2;i<=d;++i)
    {
      nk1 *=N;
      nv1 = (nv1)*(N+i)/i;
    }
  nv = cas*nv1;
  nk = cas*nk1;
  P = new Rd[nv]; // no bug correct jan 2024 FH (thank to OP)
  K = new int[nk*d1]; // no  bug correct jan 2024 FH (thank to OP)
  if( cas ==1) 
    SplitSimplex( N, P,K);
  else 
    { 
      Rd G=Rd::diag(1./(d1));
      for(int i=0;i<=d;++i)
	{
	  Rd Khat[d+1];
	  for (int j=1;j<=d;++j)
	    Khat[j][j-1]=1;// bug correct jan 2024 FH (thank to OP) 
	  Khat[i]=G;
	  SplitSimplex( N, P,K+(nk1*d1)*i,nv1*i,Khat); // FH  no recursion here ...
	}     
    }
  if(verbosity>99)
     {
         cout << "SplitSimplex : nv ="  << nv << " nk :" << nk << " " << N <<  endl ;
         for(int i=0; i< nv ; ++i)
             cout << i << " / " << P[i] <<  endl;
         
         for(int k=0,kk=0; k < nk; ++k)
         {
            cout << k << " " << kk << " : ";
             for(int m=0;m<d+1;++m)
                 cout << K[kk++] << " " ;
             cout << endl;
         }
         
     }
}

template void  SplitSimplex<R1>(int N,int & nv, R1 *& P, int & nk , int *& K);
template void  SplitSimplex<R2>(int N,int & nv, R2 *& P, int & nk , int *& K);
template void  SplitSimplex<R3>(int N,int & nv, R3 *& P, int & nk , int *& K);
/*
int main(int argc,const char ** argv)
{
  R3 *P;
  int *K,nv,nk;
  int N=2;
  if(argc>1) N=atoi(argv[1]);
  SplitSimplex(N,nv,P,nk,K);
  cout << P << " " << K << endl;
  cout << N << " nv " << nv << " nk =" << nk << endl;
  for(int i=0;i<nv;++i) cout << P[i] << endl;
  for(int i=0,l=0;i<nk;i++) 
    {
      for(int j=0;j<4;j++) 
	cout << K[l++] << " ";
      cout << endl;
    }

  
}
*/
