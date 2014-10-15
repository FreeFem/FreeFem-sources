// ********** DO NOT REMOVE THIS BANNER **********
// ORIG-DATE:     Jan 2008
// -*- Mode : c++ -*-
//
// SUMMARY  : Generic  Quadrature Formular 
// USAGE    : LGPL      
// ORG      : LJLL Universite Pierre et Marie Curi, Paris,  FRANCE 
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

 Thank to the ARN   FF2A3 grant
 ref:ANR-07-CIS7-002-01 
 */


#ifndef _QuadratureFormular_h
#define _QuadratureFormular_h

#include <iostream>

namespace Fem2D {
  using namespace std;
#include "R3.hpp"



  //class QuadratureFormular;
struct  QuadratureWeight {
  R a;
  QuadratureWeight(R aa): a(aa){}
};

template<class Rd>
class GQuadraturePoint: public QuadratureWeight,public Rd {
 public:
  typedef GQuadraturePoint QP ;
  GQuadraturePoint(): QuadratureWeight(0),Rd() {}
  GQuadraturePoint(R aa, const Rd &xx): QuadratureWeight(aa),Rd(xx) {}
  GQuadraturePoint(const Rd & xx,R aa): QuadratureWeight(aa),Rd(xx) {}
  operator R() const {return a;}
  GQuadraturePoint(R aa,R xx):QuadratureWeight(aa),Rd(xx) {}
  GQuadraturePoint(R aa,R x,R y):QuadratureWeight(aa),Rd(x,y) {}
  GQuadraturePoint(R aa,R x,R y,R z):QuadratureWeight(aa),Rd(x,y,z) {}
  template<class RD >
  GQuadraturePoint<RD>
    Bary(RD * K, double mes) { return GQuadraturePoint<RD>(this->Rd::Bary(K),a*mes);}
};

template<class Rdd>
class GQuadratureFormular {

public:
  typedef Rdd Rd;
  typedef  GQuadraturePoint<Rd> QuadraturePoint;
  typedef  GQuadraturePoint<Rd> QP;
  int exact;            // exact
  int n;                // nombre de point d'integration
  const int size;             //  size of the array
private:
  QP *p;  // les point d'integration 
  const bool clean;
public:
// -- les fonctions ------------------
  void Verification(); // for verification 
  GQuadratureFormular (int e,int NbOfNodes,QuadraturePoint *pp,bool c=false)
 :exact(e), n(NbOfNodes),size(n),p(pp),clean(c)  {Verification();}
  GQuadratureFormular (int e,int NbOfNodes,const QuadraturePoint *pp,bool c=false)
 :exact(e),n(NbOfNodes),p(pp),clean(c)   {Verification();}
  
  GQuadratureFormular(int ex,QP p0,QP p1,QP p2,QP p3,QP p4) 
    : exact(ex),n(5),size(n),p(new QP[5]),clean(true) { p[0]=p0;p[1]=p1;p[2]=p2;p[3]=p3;p[4]=p4;Verification();}
  GQuadratureFormular(int ex,QP p0,QP p1,QP p2,QP p3) 
    : exact(ex),n(4),size(n),p(new QP[4]) ,clean(true){ p[0]=p0,p[1]=p1,p[2]=p2;p[3]=p3;Verification();}
  GQuadratureFormular(int ex,QP p0,QP p1,QP p2) 
    : exact(ex),n(3),size(n),p(new QP[3]),clean(true) { p[0]=p0,p[1]=p1,p[2]=p2;Verification();}
  GQuadratureFormular(int ex,QP p0,QP p1) 
    : exact(ex),n(2),size(n),p(new QP[2]),clean(true) { p[0]=p0,p[1]=p1;Verification();}
  GQuadratureFormular(int ex,QP p0) 
    : exact(ex),n(1),size(n),p(new QP[1]),clean(true) { p[0]=p0;Verification();}
  // bluid a empty GQuadratureFormular
  GQuadratureFormular(int ssize):exact(0),n(0),size(ssize),p(new QP[size]),clean(true) {}

  const QP & operator [](int i) const  {return p[i];} 
  const QP  & operator ()(int i) const {return p[i];}
  ~GQuadratureFormular() {if(clean) delete [] p;}
    
    GQuadratureFormular(const GQuadratureFormular & QF, int mul=1)
    :exact(QF.exact),n(QF.n),size(QF.size*mul),p(new QP[size]),clean(true){ operator=(QF);}
  void operator=(const GQuadratureFormular &QF)
    {
      assert(size>=QF.n);
      n = QF.n;
      for(int i=0;i<n;++i)
           p[i]=QF.p[i];
    }
  void operator*=( double c)
    {
      for(int i=0;i<n;++i)
          p[i].a *=c;
    }
  //  Add new GQuadratureFormular on element K to the current Quadarture formular ..
    // FH   april 2014 ..
    // to compute int under levelset ..
    void reset() { n =0;exact=0;}
  template<class RD >
  void AddQFK(const GQuadratureFormular<RD> &QF ,Rd *K,double mes,int n0=0)
    {
        
        assert( size >=  n0  + QF.n );
        n0 += n;
        n = n0 + QF.n;
        for(int i=0;i<QF.n;++i)
          p[i+n0]=QF.p[i].Bary(K,mes);
        
    }
private:
 /* GQuadratureFormular(const GQuadratureFormular &)
    :exact(0),n(0),p(0){assert(0);}
  void operator=(const GQuadratureFormular &)
  {assert(0);}
  GQuadratureFormular()
    :exact(0),n(0),p(0){assert(0);}*/
  static const GQuadratureFormular * Default;
};



template<class Rd>
ostream& operator <<(ostream& , const GQuadratureFormular<Rd> & ) ;
template<class Rd>
ostream& operator <<(ostream& , GQuadraturePoint<Rd> & );

typedef GQuadratureFormular<R1> QuadratureFormular1d;


  extern const QuadratureFormular1d QF_GaussLegendre1; 
  extern const QuadratureFormular1d QF_GaussLegendre2; 
  extern const QuadratureFormular1d QF_GaussLegendre3; 
  extern const QuadratureFormular1d QF_GaussLegendre4; 
  extern const QuadratureFormular1d QF_GaussLegendre5; 
  extern const QuadratureFormular1d QF_LumpP1_1D; 
  
  
  extern const GQuadratureFormular<R2> QuadratureFormular_T_1;
  extern const GQuadratureFormular<R2> QuadratureFormular_T_1lump;
  extern const GQuadratureFormular<R2> QuadratureFormular_T_2;
  extern const GQuadratureFormular<R2> QuadratureFormular_T_5;
  extern const GQuadratureFormular<R2> QuadratureFormular_T_2_4P1;
  extern const GQuadratureFormular<R2> QuadratureFormular_T_7;
  extern const GQuadratureFormular<R2>  QuadratureFormular_T_9;

  extern const GQuadratureFormular<R3> QuadratureFormular_Tet_1;
  extern const GQuadratureFormular<R3> QuadratureFormular_Tet_1lump;
  extern const GQuadratureFormular<R3> QuadratureFormular_Tet_2;
  extern const GQuadratureFormular<R3> QuadratureFormular_Tet_5;
  


template<class Rd>
GQuadratureFormular<Rd> * QF_Simplex(int exact);
//  { return  QF_exact<GQuadratureFormular<Rd>,Rd::d+1>(exact);}

template<class Rd>
void  setQF( GQuadratureFormular<Rd> &FI,
             const GQuadratureFormular<Rd> & FI1 ,
             const GQuadratureFormular<Rd> & FI0 ,
             Rd Q[Rd::d][Rd::d+1],
             double *cmes,
             int n)
    {
        FI.reset();
        for(int i=0; i< n; ++i)
            if(cmes[i]==1)
                FI=FI1;
            else if(cmes[i] > 1e-4)
              FI.AddQFK(FI1,Q[i],cmes[i]);
            else if ( cmes[i]> 1e-8 )
              FI.AddQFK(FI0,Q[i],cmes[i]);

    }


    
}

namespace Fem2D {
  typedef  GQuadratureFormular<R2> QuadratureFormular;
  typedef  GQuadraturePoint<R2>    QuadraturePoint;
  typedef  GQuadratureFormular<R1> QuadratureFormular1d;
  typedef  GQuadraturePoint<R1>   QuadratureFormular1dPoint;
  GQuadraturePoint<R1>  *  GaussLegendre(int nn);

}
#endif
