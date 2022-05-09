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

 Thank to the ARN   FF2A3 grant
 ref:ANR-07-CIS7-002-01 
 
For 2d hight order formula: 
 Thanks to http://xyz.lanl.gov/format/math.NA/0501496
  Mathematics, abstract
  math.NA/0501496
  From: Mark Taylor [view email]
  Date: Thu, 27 Jan 2005 19:17:37 GMT   (27kb)

    Several new quadrature formulas for polynomial integration in the  triangle
    Authors:  Mark A. Taylor,  Beth A. Wingate,  Len P. Bos
    Comments: 14 pages, 14 figures, 5 pages of tabulated quadrature points
    Report-no: SAND2005-0034J

For 3d  order formula:

  from: http://www.cs.kuleuven.be/~nines/research/ecf/mtables.html

---------------------------------------------------------------------------
| Name   : Ronald  Cools               |                                  |
|          Prof. Dr. ir.               |        Address:                  |
| Email  : Ronald.Cools@cs.kuleuven.be | Katholieke Universiteit Leuven   |
|                                      | Department of Computer Science   |
| http://www.cs.kuleuven.be/~ronald/   | Celestijnenlaan 200 A            |
| Fax    : +(32) 16 32 79 96           | B-3001 HEVERLEE                  |
| Phone  : +(32) 16 32 75 62           | BELGIUM                          |
---------------------------------------------------------------------------
 */



#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
using namespace std;
#include "error.hpp"

#include "ufunction.hpp"

#include "QuadratureFormular.hpp"

namespace Fem2D {
  
template<class QF,int ON> 
  QF  * QF_exact(int exactpm,QF * p=0)
  {
    int exact = abs(exactpm);
    exact=max(0,exact);
    const int N=100;
    assert(exact<N&& exact>=0);
    static QF ** a=0;
      static int M=0;
    if(a==0)
      { //    
	a = new  QF*[N];
	assert(a);
	for(int i=0;i<N;++i)
	  a[i]=0;
      }
    assert(a && exact >=0 && exact < N);
    if ( p  )
      {
	//cout << endl << " QF " << exact << " " << p->exact << " " << p->n << endl;;
        M = max(M,exact);
	for( int i=0;i<=exact;++i)
	  {
	    if( a[i]== 0 || a[i]->n > p->n)
	      a[i]= p;
	    //  cout << " QF: on " << ON << " exact P_" << i << " : "<< a[i]->n << endl;
	  }
      }
    else if(exactpm<0)
      p=a[min(exact,M)];
    else
      p=a[exact];
    return p;
  }



template<class Rd>
GQuadratureFormular<Rd> * QF_Simplex(int exact)
 {
   return  QF_exact<GQuadratureFormular<Rd>,Rd::d+1>(exact);
 }

// explict instantiation 
template  GQuadratureFormular<R1> * QF_Simplex<R1>(int exact);
template  GQuadratureFormular<R2> * QF_Simplex<R2>(int exact);
template  GQuadratureFormular<R3> * QF_Simplex<R3>(int exact);




template<class Rd>    
ostream& operator <<(ostream& f,const  GQuadraturePoint<Rd> & p) 
      { f << '{' << (const R) p << '\t' << (const Rd &) p << '}' ; 
        return f;}

template<class Rd>    
ostream& operator <<(ostream& f, const GQuadratureFormular<Rd>  & fi) 
      { f << "nb de point integration " << fi.n << ", adr = " << &f << endl;
        for (int i=0;i<fi.n;i++) f << '\t' << fi[i] << endl; 
        return f;}

// ----------------------------------------------------------------------


static  GQuadraturePoint<R2>  P_QuadratureFormular_T_1[1] = {
  GQuadraturePoint<R2>(1.,R2(1./3.,1./3.)) };

const GQuadratureFormular<R2> QuadratureFormular_T_1(1,1,P_QuadratureFormular_T_1);
// ----------------------------------------------------------------------
static  GQuadraturePoint<R2> P_QuadratureFormular_T_1lump[3] = {
  GQuadraturePoint<R2>(1./3.,R2(0.,0.)) ,
  GQuadraturePoint<R2>(1./3.,R2(1.,0.)) ,
  GQuadraturePoint<R2>(1./3.,R2(0.,1.)) };

GQuadratureFormular<R2> const QuadratureFormular_T_1lump(1,3,P_QuadratureFormular_T_1lump);
// ----------------------------------------------------------------------

static GQuadraturePoint<R2> P_QuadratureFormular_T_2[3] = {
  GQuadraturePoint<R2>(1./3.,R2(0.5,0.5)) ,
  GQuadraturePoint<R2>(1./3.,R2(0.0,0.5)) ,
  GQuadraturePoint<R2>(1./3.,R2(0.5,0.0)) };

GQuadratureFormular<R2> const QuadratureFormular_T_2(2,3,P_QuadratureFormular_T_2);
// ----------------------------------------------------------------------
static  GQuadraturePoint<R2> P_QuadratureFormular_T_2_4P1[9] = {
  GQuadraturePoint<R2>(1./12.,R2(0.25,0.75)) ,
      GQuadraturePoint<R2>(1./12.,R2(0.75,0.25)) ,
      GQuadraturePoint<R2>(1./12.,R2(0.0,0.25)) ,
      GQuadraturePoint<R2>(1./12.,R2(0.0,0.75)) ,
      GQuadraturePoint<R2>(1./12.,R2(0.25,0.0)) ,
      GQuadraturePoint<R2>(1./12.,R2(0.75,0.0)) ,
      GQuadraturePoint<R2>(1./6.,R2(0.25,0.25)) ,
      GQuadraturePoint<R2>(1./6.,R2(0.25,0.50)) ,
      GQuadraturePoint<R2>(1./6.,R2(0.50,0.25))       
      };

GQuadratureFormular<R2> const QuadratureFormular_T_2_4P1(2,9,P_QuadratureFormular_T_2_4P1);
// ----------------------------------------------------------------------
// STROUD page  314 
// -----------------------------
const R sqrt15 = 3.87298334620741688517926539978;
const R t_T5 =1.E0/3.E0        ,                           A_T5 = 0.225E0;
const R r_T5 = (6-sqrt15)/21   ,  s_T5 = (9+2*sqrt15)/21 , B_T5 = (155-sqrt15)/1200;
const R u_T5 = (6+sqrt15)/21   ,  v_T5 = (9-2*sqrt15)/21 , C_T5 = (155+sqrt15)/1200;
// OK cette  formule  est OK 
static  GQuadraturePoint<R2> P_QuadratureFormular_T_5[] = {
  GQuadraturePoint<R2>(A_T5,R2(t_T5,t_T5)),
  GQuadraturePoint<R2>(B_T5,R2(r_T5,r_T5)),
  GQuadraturePoint<R2>(B_T5,R2(r_T5,s_T5)),
  GQuadraturePoint<R2>(B_T5,R2(s_T5,r_T5)),
  GQuadraturePoint<R2>(C_T5,R2(u_T5,u_T5)),
  GQuadraturePoint<R2>(C_T5,R2(u_T5,v_T5)),
  GQuadraturePoint<R2>(C_T5,R2(v_T5,u_T5))
};
const GQuadratureFormular<R2> QuadratureFormular_T_5(5,7,P_QuadratureFormular_T_5);
// ------------------
//----

 
// Thanks to http://xyz.lanl.gov/format/math.NA/0501496
/*
Mathematics, abstract
math.NA/0501496
From: Mark Taylor [view email]
Date: Thu, 27 Jan 2005 19:17:37 GMT   (27kb)

Several new quadrature formulas for polynomial integration in the  triangle
Authors:  Mark A. Taylor,  Beth A. Wingate,  Len P. Bos
Comments: 14 pages, 14 figures, 5 pages of tabulated quadrature points
Report-no: SAND2005-0034J

*/
// ----------------------------------------------------------------------
// awk '/15:/,/21:/ {print "GQuadraturePoint<R2>(" $3 "/2," $1"," $2"),"}' coords.txt   
static  GQuadraturePoint<R2> P_QuadratureFormular_T_7[] = {  
GQuadraturePoint<R2>(0.0102558174092/2,R2(1.0000000000000,0.0000000000000)),
GQuadraturePoint<R2>(0.0102558174092/2,R2(0.0000000000000,0.0000000000000)),
GQuadraturePoint<R2>(0.0102558174092/2,R2(0.0000000000000,1.0000000000000)),
GQuadraturePoint<R2>(0.1116047046647/2,R2(0.7839656651012,0.0421382841642)),
GQuadraturePoint<R2>(0.1116047046647/2,R2(0.1738960507345,0.7839656651012)),
GQuadraturePoint<R2>(0.1116047046647/2,R2(0.1738960507345,0.0421382841642)),
GQuadraturePoint<R2>(0.1116047046647/2,R2(0.0421382841642,0.1738960507345)),
GQuadraturePoint<R2>(0.1116047046647/2,R2(0.7839656651012,0.1738960507345)),
GQuadraturePoint<R2>(0.1116047046647/2,R2(0.0421382841642,0.7839656651012)),
GQuadraturePoint<R2>(0.1679775595335/2,R2(0.4743880861752,0.4743880861752)),
GQuadraturePoint<R2>(0.1679775595335/2,R2(0.4743880861752,0.0512238276497)),
GQuadraturePoint<R2>(0.1679775595335/2,R2(0.0512238276497,0.4743880861752)),
GQuadraturePoint<R2>(0.2652238803946/2,R2(0.2385615300181,0.5228769399639)),
GQuadraturePoint<R2>(0.2652238803946/2,R2(0.5228769399639,0.2385615300181)),
GQuadraturePoint<R2>(0.2652238803946/2,R2(0.2385615300181,0.2385615300181))
};
const GQuadratureFormular<R2> QuadratureFormular_T_7(7,15,P_QuadratureFormular_T_7);

// awk '/21:/,/28:/ {print "GQuadraturePoint<R2>(" $3 "/2," $1"," $2"),"}' coords.txt
static GQuadraturePoint<R2> P_QuadratureFormular_T_9[] = {  
GQuadraturePoint<R2>(0.0519871420646/2,R2(0.0451890097844,0.0451890097844)),
GQuadraturePoint<R2>(0.0519871420646/2,R2(0.0451890097844,0.9096219804312)),
GQuadraturePoint<R2>(0.0519871420646/2,R2(0.9096219804312,0.0451890097844)),
GQuadraturePoint<R2>(0.0707034101784/2,R2(0.7475124727339,0.0304243617288)),
GQuadraturePoint<R2>(0.0707034101784/2,R2(0.2220631655373,0.0304243617288)),
GQuadraturePoint<R2>(0.0707034101784/2,R2(0.7475124727339,0.2220631655373)),
GQuadraturePoint<R2>(0.0707034101784/2,R2(0.2220631655373,0.7475124727339)),
GQuadraturePoint<R2>(0.0707034101784/2,R2(0.0304243617288,0.7475124727339)),
GQuadraturePoint<R2>(0.0707034101784/2,R2(0.0304243617288,0.2220631655373)),
GQuadraturePoint<R2>(0.0909390760952/2,R2(0.1369912012649,0.2182900709714)),
GQuadraturePoint<R2>(0.0909390760952/2,R2(0.6447187277637,0.2182900709714)),
GQuadraturePoint<R2>(0.0909390760952/2,R2(0.1369912012649,0.6447187277637)),
GQuadraturePoint<R2>(0.0909390760952/2,R2(0.2182900709714,0.6447187277637)),
GQuadraturePoint<R2>(0.0909390760952/2,R2(0.2182900709714,0.1369912012649)),
GQuadraturePoint<R2>(0.0909390760952/2,R2(0.6447187277637,0.1369912012649)),
GQuadraturePoint<R2>(0.1032344051380/2,R2(0.0369603304334,0.4815198347833)),
GQuadraturePoint<R2>(0.1032344051380/2,R2(0.4815198347833,0.0369603304334)),
GQuadraturePoint<R2>(0.1032344051380/2,R2(0.4815198347833,0.4815198347833)),
GQuadraturePoint<R2>(0.1881601469167/2,R2(0.4036039798179,0.1927920403641)),
GQuadraturePoint<R2>(0.1881601469167/2,R2(0.4036039798179,0.4036039798179)),
GQuadraturePoint<R2>(0.1881601469167/2,R2(0.1927920403641,0.4036039798179))
};
const GQuadratureFormular<R2> QuadratureFormular_T_9(9,21,P_QuadratureFormular_T_9);





  // Computation of nodes and weights for a Gauss-Legendre quadrature formula
  // on  [0,1] seg.
typedef GQuadraturePoint<R1> QP1;
QP1 *  GaussLegendre(int nn)
// : exact(nn*2-1),n(nn),p(new Point[nn])
  {
    QP1 *p=new QP1[nn];
    int n=nn;
    double r, r1, p1, p2, p3, dp3;
    double eps=1e-16;
      const double pi=3.141592653589793238463;//M_PI;
    
    for(int i = 0,ii=n; i <= (n+1)/2-1; i++)
      { 
	ii-=1;
	r = cos(pi*(4*i+3)/(4*n+2));
        do
	  {
            p2 = 0;
            p3 = 1;
            for(int j = 0; j <= n-1; j++)
	      {
                p1 = p2;
                p2 = p3;
                p3 = ((2*j+1)*r*p2-j*p1)/(j+1);
	      }
	    
            dp3 = n*(r*p3-p2)/(r*r-1);
            r1 = r;
            r = r-p3/dp3;
	  }
        while(fabs(r-r1)>=eps*(1+fabs(r))*100);
        p[i].x = 0.5 +r/2;
        p[ii].x =0.5  -r/2;
        p[i].a=p[ii].a= 1./((1.-r*r)*dp3*dp3);
	
      }
    
    //    Check();
    return p;
  }
/*
void QuadratureFormular1d::Check()
{
    QF_exact<QuadratureFormular1d,2>(exact,this);
    int err=0;
    for(int m=0;m<=exact;++m)
    {
	R ve = 1./ ( m+1), v=0.;
	for (int i=0;i<n;++i)
	    v += p[i].a*pow(p[i].x,m);
	if (abs( ve-v)/ve > 1.e-8)
	{
	    cout << " erreur QuadratureFormular1d  n= " << n  << " exact = " << exact << endl;
	    cout << " int x^" <<m << " == " << ve << " ! = " << v  << endl; 
	    err++;
	}
    }
    assert(err==0);	
}
*/
static const R gauss_n2_1=  (1-sqrt(1./3.))/2;
static const R gauss_n2_2=  1 - gauss_n2_1;

const R gauss_n3_0=  0.5 ;
const R gauss_n3_1=  (1-sqrt(3./5.)) /2  ;
const R gauss_n3_2 =  1 - gauss_n3_1 ;

const R pgauss_n3_0=  8./18.;
const R pgauss_n3_1=  5./18.;
const R pgauss_n3_2=  5./18.;

const R pgauss1_n4_a= (18.+sqrt(30.))/36.;
const R pgauss1_n4_b= (18.-sqrt(30.))/36.;
const R gauss1_n4_a=  sqrt( 525. - 70.* sqrt(30.))/35.;
const R gauss1_n4_b=  sqrt( 525. + 70.* sqrt(30.))/35.;

const R pgauss1_n5_0 = 128./225.;
const R pgauss1_n5_a= (322+13.*sqrt(70.))/900.;
const R pgauss1_n5_b= (322-13.*sqrt(70.))/900.;
const R gauss1_n5_a=  sqrt(245.-14.*sqrt(70.))/21.;
const R gauss1_n5_b=  sqrt(245.+14.*sqrt(70.))/21.;

const QuadratureFormular1d QF_GaussLegendre3(5,
                  QuadratureFormular1d::QP(pgauss_n3_0,gauss_n3_0),
                  QuadratureFormular1d::QP(pgauss_n3_1,gauss_n3_1),
                  QuadratureFormular1d::QP(pgauss_n3_2,gauss_n3_2)); 
                    
const QuadratureFormular1d QF_GaussLegendre2(3,
                  QuadratureFormular1d::QP(0.5,gauss_n2_1),
                  QuadratureFormular1d::QP(0.5,gauss_n2_2)); 

const QuadratureFormular1d QF_GaussLegendre1(1,QuadratureFormular1d::QP(1,0.5)); 

const QuadratureFormular1d QF_LumpP1_1D(1,
                  QuadratureFormular1d::QP(0.5,0.),
                  QuadratureFormular1d::QP(0.5,1.)); 


const QuadratureFormular1d QF_GaussLegendre4(7,
					     QuadratureFormular1d::QP(pgauss1_n4_a/2.,(1.+gauss1_n4_a)/2.),
					     QuadratureFormular1d::QP(pgauss1_n4_a/2.,(1.-gauss1_n4_a)/2.),
					     QuadratureFormular1d::QP(pgauss1_n4_b/2.,(1.+gauss1_n4_b)/2.),
					     QuadratureFormular1d::QP(pgauss1_n4_b/2.,(1.-gauss1_n4_b)/2.)
					     ); 

const QuadratureFormular1d QF_GaussLegendre5(-1+2*5,5,GaussLegendre(5),true);
const QuadratureFormular1d QF_GaussLegendre6(-1+2*6,6,GaussLegendre(6),true);
const QuadratureFormular1d QF_GaussLegendre7(-1+2*7,7,GaussLegendre(7),true);
const QuadratureFormular1d QF_GaussLegendre8(-1+2*8,8,GaussLegendre(8),true);
const QuadratureFormular1d QF_GaussLegendre9(-1+2*9,9,GaussLegendre(9),true);
const QuadratureFormular1d QF_GaussLegendre10(-1+2*10,10,GaussLegendre(10),true);



/*
Region: Simplex 
Dimension: 3 
Degree: 3 
Points: 5 
Structure: Fully symmetric 
Rule struct: 0 0 0 0 0 0 1 1 0 0 0 
Generator: [ Fully symmetric ] 
( 0.25, 0.25, 0.25, ) 
Corresponding weight: 
-0.133333333333333333333333333333333,
Generator: [ Fully symmetric ] 
( 0.166666666666666666666666666666666, 
0.166666666666666666666666666666666, 
0.166666666666666666666666666666666, 
) 
Corresponding weight: 
0.075,

Region: Simplex 
Dimension: 3 
Degree: 4 
Points: 11 
Structure: Fully symmetric 
Rule struct: 0 0 0 0 0 0 1 1 1 0 0 
Generator: [ Fully symmetric ] 
( 0.25, 0.25, 0.25, ) 
Corresponding weight: 
-0.0131555555555555555555555555555555,
Generator: [ Fully symmetric ] 
( 0.0714285714285714285714285714285714, 
0.0714285714285714285714285714285714, 
0.0714285714285714285714285714285714, 
) 
Corresponding weight: 
10 ^ -3 x 7.62222222222222222222222222222222,

Generator: [ Fully symmetric ] 
( 0.399403576166799204996102147461640, 
0.399403576166799204996102147461640, 
0.100596423833200795003897852538359, 
) 
Corresponding weight: 
0.0248888888888888888888888888888888,


Region: Simplex 
Dimension: 3 
Degree: 6 
Points: 24 
Structure: Fully symmetric 
Rule struct: 0 0 0 0 0 0 0 3 0 1 0 
Generator: [ Fully symmetric ] 
( 0.214602871259152029288839219386284, 
0.214602871259152029288839219386284, 
0.214602871259152029288839219386284, 
) 
Corresponding weight: 
10 ^ -3 x 6.65379170969458201661510459291332,
Generator: [ Fully symmetric ] 
( 0.0406739585346113531155794489564100, 
0.0406739585346113531155794489564100, 
0.0406739585346113531155794489564100, 
) 
Corresponding weight: 
10 ^ -3 x 1.67953517588677382466887290765614,

Generator: [ Fully symmetric ] 
( 0.322337890142275510343994470762492, 
0.322337890142275510343994470762492, 
0.322337890142275510343994470762492, 
) 
Corresponding weight: 
10 ^ -3 x 9.22619692394245368252554630895433,

Generator: [ Fully symmetric ] 
( 0.0636610018750175252992355276057269, 
0.0636610018750175252992355276057269, 
0.269672331458315808034097805727606, 
) 
Corresponding weight: 
10 ^ -3 x 8.03571428571428571428571428571428,



Region: Simplex 
Dimension: 3 
Degree: 7 
Points: 31 
Structure: Fully symmetric 
Rule struct: 0 1 0 0 0 0 1 3 0 1 0 
Generator: [ Fully symmetric ] 
( 0.5, 0.5, 0., ) 
Corresponding weight: 
10 ^ -4 x 9.70017636684303350970017636684303,
Generator: [ Fully symmetric ] 
( 0.25, 0.25, 0.25, ) 
Corresponding weight: 
0.0182642234661088202912015685649462,

Generator: [ Fully symmetric ] 
( 0.0782131923303180643739942508375545, 
0.0782131923303180643739942508375545, 
0.0782131923303180643739942508375545, 
) 
Corresponding weight: 
0.0105999415244136869164138748545257,

Generator: [ Fully symmetric ] 
( 0.121843216663905174652156372684818, 
0.121843216663905174652156372684818, 
0.121843216663905174652156372684818, 
) 
Corresponding weight: 
-0.0625177401143318516914703474927900,

Generator: [ Fully symmetric ] 
( 0.332539164446420624152923823157707, 
0.332539164446420624152923823157707, 
0.332539164446420624152923823157707, 
) 
Corresponding weight: 
10 ^ -3 x 4.89142526307349938479576303671027,

Generator: [ Fully symmetric ] 
( 0.1, 0.1, 0.2, ) 
Corresponding weight: 
0.0275573192239858906525573192239858,



Region: Simplex 
Dimension: 3 
Degree: 8 
Points: 43 
Structure: Fully symmetric 
Rule struct: 0 0 0 0 0 0 1 3 1 2 0 
Generator: [ Fully symmetric ] 
( 0.25, 0.25, 0.25, ) 
Corresponding weight: 
-0.0205001886586399158405865177642941,
Generator: [ Fully symmetric ] 
( 0.206829931610673204083980900024961, 
0.206829931610673204083980900024961, 
0.206829931610673204083980900024961, 
) 
Corresponding weight: 
0.0142503058228669012484397415358704,

Generator: [ Fully symmetric ] 
( 0.0821035883105467230906058078714215, 
0.0821035883105467230906058078714215, 
0.0821035883105467230906058078714215, 
) 
Corresponding weight: 
10 ^ -3 x 1.96703331313390098756280342445466,

Generator: [ Fully symmetric ] 
( 10 ^ -3 x 5.78195050519799725317663886414270, 
10 ^ -3 x 5.78195050519799725317663886414270, 
10 ^ -3 x 5.78195050519799725317663886414270, 
) 
Corresponding weight: 
10 ^ -4 x 1.69834109092887379837744566704016,

Generator: [ Fully symmetric ] 
( 0.0505327400188942244256245285579071, 
0.0505327400188942244256245285579071, 
0.449467259981105775574375471442092, 
) 
Corresponding weight: 
10 ^ -3 x 4.57968382446728180074351446297276,

Generator: [ Fully symmetric ] 
( 0.229066536116811139600408854554753, 
0.229066536116811139600408854554753, 
0.0356395827885340437169173969506114, 
) 
Corresponding weight: 
10 ^ -3 x 5.70448580868191850680255862783040,

Generator: [ Fully symmetric ] 
( 0.0366077495531974236787738546327104, 
0.0366077495531974236787738546327104, 
0.190486041934633455699433285315099, 
) 
Corresponding weight: 
10 ^ -3 x 2.14051914116209259648335300092023,

Region: Simplex 
Dimension: 3 
Degree: 9 
Points: 53 
Structure: Fully symmetric 
Rule struct: 0 0 0 0 0 0 1 4 0 3 0 
Generator: [ Fully symmetric ] 
( 0.25, 0.25, 0.25, ) 
Corresponding weight: 
-0.137799038326108641245781827743502,
Generator: [ Fully symmetric ] 
( 0.0483510385497367408791638973335652, 
0.0483510385497367408791638973335652, 
0.0483510385497367408791638973335652, 
) 
Corresponding weight: 
10 ^ -3 x 1.86533656908528954751732956352791,

Generator: [ Fully symmetric ] 
( 0.324579280117882365858796772348839, 
0.324579280117882365858796772348839, 
0.324579280117882365858796772348839, 
) 
Corresponding weight: 
10 ^ -3 x 4.30942396949340069685481460143672,

Generator: [ Fully symmetric ] 
( 0.114616540223995219683790300317451, 
0.114616540223995219683790300317451, 
0.114616540223995219683790300317451, 
) 
Corresponding weight: 
-0.0901847664812015251273931145424214,

Generator: [ Fully symmetric ] 
( 0.225489951911513918847703826550363, 
0.225489951911513918847703826550363, 
0.225489951911513918847703826550363, 
) 
Corresponding weight: 
0.0446725762025114446937850368088403,

Generator: [ Fully symmetric ] 
( 0.131627809246869809838976179197822, 
0.131627809246869809838976179197822, 
0.0836647016171849679665385541714202, 
) 
Corresponding weight: 
0.0347004058845507618227002576570271,

Generator: [ Fully symmetric ] 
( 0.433951461411406772411426685800524, 
0.433951461411406772411426685800524, 
0.107769859549428611825218473500043, 
) 
Corresponding weight: 
10 ^ -3 x 3.35258390266064697010720483063291,

Generator: [ Fully symmetric ] 
( 10 ^ -3 x -1.37627731813820071002030321419028, 
10 ^ -3 x -1.37627731813820071002030321419028, 
0.276553472636807342120192082186261, 
) 
Corresponding weight: 
10 ^ -4 x 4.31628875556996929641889902726182,


 */
template<class Rd>
void GQuadratureFormular<Rd>::Verification()
{
  QF_exact<GQuadratureFormular<Rd>,Rd::d+1>(exact,this);
  const int d=Rd::d;
  const double tol=1.e-12;
  R err=0;
  R a[d+1],h;
  
  for (int k=0;k<=exact;  k++)
    {

      //  formule magic: int_K  \lamda^k =   d! k! / ( d+k)! 
      
      R sa[d+1];
      for(int l=0;l<=d;++l)	      
	sa[l]=0.;
     
      for (int j=0;j<n;j++)
	{
	  h = p[j];
	  Rd P = p[j];
	  
	  for(int l=0;l<d;++l)	      
	    a[l]=p[j][l];
	  a[d] = 1.-p[j].sum();
	  
	  for(int l=0;l<=d;++l)	      	  
	    sa[l]+=h*pow(a[l],k);  
	}
      
      R se(1),see(1);
      for (int i=1;i<=k;i++) 
	se *= (R) i / (R) (i+d);
      see=se;
      
      for(int l=0;l<=d;++l)
	err = Max(err,Abs(se-sa[l]));
      
      if (err>tol)
	{
	  cerr << " d= " << d << "T Ordre= " << k << " d!k!/(d+k)!= " << se << " " ;
	  for(int l=0;l<=d;++l)
	    cerr << sa[l] << " ";
	  cerr << " err= " << err << endl;
	}
    }
  
  if(err>tol)
    {
      cerr << "Erreur dans la formule d'integration d=" <<d  << " exact = " << exact 
	   << " Nb Point = " << n << endl;
      assert(0);
    }
    
}
/*
  from:
http://www.cs.kuleuven.be/~nines/research/ecf/mtables.html
 */
typedef GQuadraturePoint<R3>  PQP3;
typedef GQuadratureFormular<R3>  PQF3;
// 2 4 (formule 2)
// A.H. Stroud, Approximate calculation of multiple integrals, Prentice-Hall, Englewood Cliffs, N.J., 1971.
// JavaScript:formule('t3-2-4b')
    
    PQP3  QF_TET_1[]=
    { 
    PQP3(R3(0.25,0.25,0.25), 1)
    };
    PQF3 const QuadratureFormular_Tet_1(1,1,QF_TET_1);
 
    
    PQP3  QF_TET_1lump[]=
    { 
    PQP3(R3(0,0,0),0.25),
    PQP3(R3(1.,0,0),0.25),
    PQP3(R3(0,1.,0),0.25),
    PQP3(R3(0,0,1.),0.25),
    };
    PQF3 const QuadratureFormular_Tet_1lump(1,4,QF_TET_1lump);
    
    
PQP3  QF_TET_2[]=
  { 
     PQP3(R3(0.58541019662496845446137605030968, 0.138196601125010515179541316563436, 0.138196601125010515179541316563436), 0.25),
     PQP3(R3(0.138196601125010515179541316563436,  0.58541019662496845446137605030968, 0.138196601125010515179541316563436), 0.25),
     PQP3(R3(0.138196601125010515179541316563436, 0.138196601125010515179541316563436,  0.58541019662496845446137605030968),0.25),
     PQP3(R3(0.138196601125010515179541316563436, 0.138196601125010515179541316563436, 0.138196601125010515179541316563436),0.25)
  };
  PQF3 const QuadratureFormular_Tet_2(2,4,QF_TET_2);
// 5  14  (formule 1)
/*
GM78
A. Grundmann and H.M. Möller, Invariant integration formulas for the n-simplex by combinatorial methods, SIAM J. Numer. Anal. 15 (1978), 282--290.
 */
PQP3  QF_TET_5[]=
   { 
     PQP3(R3(0.7217942490673263207930282587889082 , 0.0927352503108912264023239137370306 , 0.0927352503108912264023239137370306) , 0.0122488405193936582572850342477212*6.),
     PQP3(R3(0.0927352503108912264023239137370306 , 0.7217942490673263207930282587889082 , 0.0927352503108912264023239137370306) , 0.0122488405193936582572850342477212*6.),
     PQP3(R3(0.0927352503108912264023239137370306 , 0.0927352503108912264023239137370306 , 0.7217942490673263207930282587889082) , 0.0122488405193936582572850342477212*6.),
     PQP3(R3(0.0927352503108912264023239137370306 , 0.0927352503108912264023239137370306 , 0.0927352503108912264023239137370306) , 0.0122488405193936582572850342477212*6),
     PQP3(R3(0.067342242210098170607962798709629 , 0.310885919263300609797345733763457 , 0.310885919263300609797345733763457) , 0.0187813209530026417998642753888810*6.),
     PQP3(R3(0.310885919263300609797345733763457 , 0.067342242210098170607962798709629 , 0.310885919263300609797345733763457) , 0.0187813209530026417998642753888810*6.),
     PQP3(R3(0.310885919263300609797345733763457 , 0.310885919263300609797345733763457 , 0.067342242210098170607962798709629) , 0.0187813209530026417998642753888810*6.),
     PQP3(R3(0.310885919263300609797345733763457 , 0.310885919263300609797345733763457 , 0.310885919263300609797345733763457) , 0.0187813209530026417998642753888810*6.),
     PQP3(R3(0.454496295874350350508119473720660,  0.454496295874350350508119473720660, 0.045503704125649649491880526279339),  7.09100346284691107301157135337624E-3*6.),
     PQP3(R3(00.454496295874350350508119473720660,  0.045503704125649649491880526279339, 0.454496295874350350508119473720660), 7.09100346284691107301157135337624E-3*6.),
     PQP3(R3(00.045503704125649649491880526279339,  0.454496295874350350508119473720660, 0.454496295874350350508119473720660), 7.09100346284691107301157135337624E-3*6.),
     PQP3(R3(00.045503704125649649491880526279339,  0.045503704125649649491880526279339, 0.454496295874350350508119473720660), 7.09100346284691107301157135337624E-3*6.),
     PQP3(R3(00.045503704125649649491880526279339,  0.454496295874350350508119473720660, 0.045503704125649649491880526279339), 7.09100346284691107301157135337624E-3*6.),
     PQP3(R3(00.454496295874350350508119473720660,  0.045503704125649649491880526279339, 0.045503704125649649491880526279339), 7.09100346284691107301157135337624E-3*6.)
   };

  PQF3 const QuadratureFormular_Tet_5(5,14,QF_TET_5);

  template<> const GQuadratureFormular<R1> * GQuadratureFormular<R1>::Default= & QF_GaussLegendre3;
  template<> const GQuadratureFormular<R2> * GQuadratureFormular<R2>::Default= & QuadratureFormular_T_5;
  template<> const GQuadratureFormular<R3> * GQuadratureFormular<R3>::Default= & QuadratureFormular_Tet_5;

    // Instanciation explicite of icpc   compiler  Version 11.1   (FH  oct 2010) 
    template void  GQuadratureFormular<R1>::Verification();
    template void  GQuadratureFormular<R2>::Verification();
    template void  GQuadratureFormular<R3>::Verification();
    
   
    
} // end namespace Fem2D
