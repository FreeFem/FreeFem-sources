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
};

template<class Rdd>
class GQuadratureFormular {

public:
  typedef Rdd Rd;
  typedef  GQuadraturePoint<Rd> QuadraturePoint;
  typedef  GQuadraturePoint<Rd> QP;
  const int exact;            // exact
  const int n;                // nombre de point d'integration
private:
  QP *p;  // les point d'integration 
  const bool clean;
public:
// -- les fonctions ------------------
  void Verification(); // for verification 
  GQuadratureFormular (int e,int NbOfNodes,QuadraturePoint *pp,bool c=false)
 :exact(e), n(NbOfNodes),p(pp),clean(c)  {Verification();}
  GQuadratureFormular (int e,int NbOfNodes,const QuadraturePoint *pp,bool c=false)
 :exact(e),n(NbOfNodes),p(pp),clean(c)   {Verification();}
  
  GQuadratureFormular(int ex,QP p0,QP p1,QP p2,QP p3,QP p4) 
    : exact(ex),n(5),p(new QP[5]),clean(true) { p[0]=p0;p[1]=p1;p[2]=p2;p[3]=p3;p[4]=p4;Verification();}
  GQuadratureFormular(int ex,QP p0,QP p1,QP p2,QP p3) 
    : exact(ex),n(4),p(new QP[4]) ,clean(true){ p[0]=p0,p[1]=p1,p[2]=p2;p[3]=p3;Verification();}
  GQuadratureFormular(int ex,QP p0,QP p1,QP p2) 
    : exact(ex),n(3),p(new QP[3]),clean(true) { p[0]=p0,p[1]=p1,p[2]=p2;Verification();}
  GQuadratureFormular(int ex,QP p0,QP p1) 
    : exact(ex),n(2),p(new QP[2]),clean(true) { p[0]=p0,p[1]=p1;Verification();}
  GQuadratureFormular(int ex,QP p0) 
    : exact(ex),n(1),p(new QP[1]),clean(true) { p[0]=p0;Verification();}


  const QP & operator [](int i) const  {return p[i];} 
  const QP  & operator ()(int i) const {return p[i];}
  ~GQuadratureFormular() {if(clean) delete [] p;}
private:
  GQuadratureFormular(const GQuadratureFormular &)
    :exact(0),n(0),p(0){assert(0);}
  void operator=(const GQuadratureFormular &)
  {assert(0);}
  GQuadratureFormular()
    :exact(0),n(0),p(0){assert(0);}
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

  extern const GQuadratureFormular<R3> QuadratureFormular_Tet_2;
  extern const GQuadratureFormular<R3> QuadratureFormular_Tet_5;
  


template<class Rd>
GQuadratureFormular<Rd> * QF_Simplex(int exact);
//  { return  QF_exact<GQuadratureFormular<Rd>,Rd::d+1>(exact);}



      
}

namespace Fem2D {
  typedef  GQuadratureFormular<R2> QuadratureFormular;
  typedef  GQuadraturePoint<R2>    QuadraturePoint;
  typedef  GQuadratureFormular<R1> QuadratureFormular1d;
  typedef  GQuadraturePoint<R1>   QuadratureFormular1dPoint;
  GQuadraturePoint<R1>  *  GaussLegendre(int nn);

}
#endif
