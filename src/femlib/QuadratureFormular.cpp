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

#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include "error.hpp"
using namespace std;


#include "RNM.hpp"
#include "rgraph.hpp"
#include "fem.hpp"
#include "QuadratureFormular.hpp"
namespace Fem2D {

ostream& operator <<(ostream& f,const  QuadraturePoint & p) 
      { f << '{' << (const R) p << '\t' << (const R2) p << '}' ; 
        return f;}

ostream& operator <<(ostream& f, const QuadratureFormular & fi) 
      { f << "nb de point integration " << fi.n << ", adr = " << &f << endl;
        for (int i=0;i<fi.n;i++) f << '\t' << fi[i] << endl; 
        return f;}

// ----------------------------------------------------------------------

static const QuadraturePoint P_QuadratureFormular_T_1[1] = {
      QuadraturePoint(1.,1./3.,1./3.) };

const QuadratureFormular QuadratureFormular_T_1(3,1,1,P_QuadratureFormular_T_1);
// ----------------------------------------------------------------------
static const QuadraturePoint P_QuadratureFormular_T_1lump[3] = {
      QuadraturePoint(1./3.,0.,0.) ,
      QuadraturePoint(1./3.,1.,0.) ,
      QuadraturePoint(1./3.,0.,1.) };

QuadratureFormular const QuadratureFormular_T_1lump(3,1,3,P_QuadratureFormular_T_1lump);
// ----------------------------------------------------------------------

static const QuadraturePoint P_QuadratureFormular_T_2[3] = {
      QuadraturePoint(1./3.,0.5,0.5) ,
      QuadraturePoint(1./3.,0.0,0.5) ,
      QuadraturePoint(1./3.,0.5,0.0) };

QuadratureFormular const QuadratureFormular_T_2(3,2,3,P_QuadratureFormular_T_2);
// ----------------------------------------------------------------------
static const QuadraturePoint P_QuadratureFormular_T_2_4P1[9] = {
      QuadraturePoint(1./12.,0.25,0.75) ,
      QuadraturePoint(1./12.,0.75,0.25) ,
      QuadraturePoint(1./12.,0.0,0.25) ,
      QuadraturePoint(1./12.,0.0,0.75) ,
      QuadraturePoint(1./12.,0.25,0.0) ,
      QuadraturePoint(1./12.,0.75,0.0) ,
      QuadraturePoint(1./6.,0.25,0.25) ,
      QuadraturePoint(1./6.,0.25,0.50) ,
      QuadraturePoint(1./6.,0.50,0.25)       
      };

QuadratureFormular const QuadratureFormular_T_2_4P1(3,2,9,P_QuadratureFormular_T_2_4P1);
// ----------------------------------------------------------------------
// STROUD page  314 
// -----------------------------
const R sqrt15 = 3.87298334620741688517926539978;
const R t_T5 =1.E0/3.E0        ,                           A_T5 = 0.225E0;
const R r_T5 = (6-sqrt15)/21   ,  s_T5 = (9+2*sqrt15)/21 , B_T5 = (155-sqrt15)/1200;
const R u_T5 = (6+sqrt15)/21   ,  v_T5 = (9-2*sqrt15)/21 , C_T5 = (155+sqrt15)/1200;
// OK cette  formule  est OK 
static const QuadraturePoint P_QuadratureFormular_T_5[] = {
      QuadraturePoint(A_T5,t_T5,t_T5),
      QuadraturePoint(B_T5,r_T5,r_T5),
      QuadraturePoint(B_T5,r_T5,s_T5),
      QuadraturePoint(B_T5,s_T5,r_T5),
      QuadraturePoint(C_T5,u_T5,u_T5),
      QuadraturePoint(C_T5,u_T5,v_T5),
      QuadraturePoint(C_T5,v_T5,u_T5)
};
const QuadratureFormular QuadratureFormular_T_5(3,5,7,P_QuadratureFormular_T_5);
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
// awk '/15:/,/21:/ {print "QuadraturePoint(" $3 "/2," $1"," $2"),"}' coords.txt   
static const QuadraturePoint P_QuadratureFormular_T_7[] = {  
QuadraturePoint(0.0102558174092/2,1.0000000000000,0.0000000000000),
QuadraturePoint(0.0102558174092/2,0.0000000000000,0.0000000000000),
QuadraturePoint(0.0102558174092/2,0.0000000000000,1.0000000000000),
QuadraturePoint(0.1116047046647/2,0.7839656651012,0.0421382841642),
QuadraturePoint(0.1116047046647/2,0.1738960507345,0.7839656651012),
QuadraturePoint(0.1116047046647/2,0.1738960507345,0.0421382841642),
QuadraturePoint(0.1116047046647/2,0.0421382841642,0.1738960507345),
QuadraturePoint(0.1116047046647/2,0.7839656651012,0.1738960507345),
QuadraturePoint(0.1116047046647/2,0.0421382841642,0.7839656651012),
QuadraturePoint(0.1679775595335/2,0.4743880861752,0.4743880861752),
QuadraturePoint(0.1679775595335/2,0.4743880861752,0.0512238276497),
QuadraturePoint(0.1679775595335/2,0.0512238276497,0.4743880861752),
QuadraturePoint(0.2652238803946/2,0.2385615300181,0.5228769399639),
QuadraturePoint(0.2652238803946/2,0.5228769399639,0.2385615300181),
QuadraturePoint(0.2652238803946/2,0.2385615300181,0.2385615300181)
};
const QuadratureFormular QuadratureFormular_T_7(3,7,15,P_QuadratureFormular_T_7);

// awk '/21:/,/28:/ {print "QuadraturePoint(" $3 "/2," $1"," $2"),"}' coords.txt
static const QuadraturePoint P_QuadratureFormular_T_9[] = {  
QuadraturePoint(0.0519871420646/2,0.0451890097844,0.0451890097844),
QuadraturePoint(0.0519871420646/2,0.0451890097844,0.9096219804312),
QuadraturePoint(0.0519871420646/2,0.9096219804312,0.0451890097844),
QuadraturePoint(0.0707034101784/2,0.7475124727339,0.0304243617288),
QuadraturePoint(0.0707034101784/2,0.2220631655373,0.0304243617288),
QuadraturePoint(0.0707034101784/2,0.7475124727339,0.2220631655373),
QuadraturePoint(0.0707034101784/2,0.2220631655373,0.7475124727339),
QuadraturePoint(0.0707034101784/2,0.0304243617288,0.7475124727339),
QuadraturePoint(0.0707034101784/2,0.0304243617288,0.2220631655373),
QuadraturePoint(0.0909390760952/2,0.1369912012649,0.2182900709714),
QuadraturePoint(0.0909390760952/2,0.6447187277637,0.2182900709714),
QuadraturePoint(0.0909390760952/2,0.1369912012649,0.6447187277637),
QuadraturePoint(0.0909390760952/2,0.2182900709714,0.6447187277637),
QuadraturePoint(0.0909390760952/2,0.2182900709714,0.1369912012649),
QuadraturePoint(0.0909390760952/2,0.6447187277637,0.1369912012649),
QuadraturePoint(0.1032344051380/2,0.0369603304334,0.4815198347833),
QuadraturePoint(0.1032344051380/2,0.4815198347833,0.0369603304334),
QuadraturePoint(0.1032344051380/2,0.4815198347833,0.4815198347833),
QuadraturePoint(0.1881601469167/2,0.4036039798179,0.1927920403641),
QuadraturePoint(0.1881601469167/2,0.4036039798179,0.4036039798179),
QuadraturePoint(0.1881601469167/2,0.1927920403641,0.4036039798179)
};
const QuadratureFormular QuadratureFormular_T_9(3,9,21,P_QuadratureFormular_T_9);


// ----------------------------------------------------------------------

static const QuadraturePoint P_QuadratureFormular_Q_1[1] = {
       QuadraturePoint(1.,0.5,0.5) };

const QuadratureFormular QuadratureFormular_Q_1(4,1,1,P_QuadratureFormular_Q_1);

// const gauss_n2_1=  0.7113248654051872 ;
// const gauss_n2_2=  0.2886751345948128 ;
static const R gauss_n2_1=  (1-sqrt(1./3.))/2;
static const R gauss_n2_2=  1 - gauss_n2_1;

static const QuadraturePoint P_QuadratureFormular_Q_3[4] = {
      QuadraturePoint(.25, gauss_n2_1, gauss_n2_1 ),
      QuadraturePoint(.25, gauss_n2_1, gauss_n2_2 ), 
      QuadraturePoint(.25, gauss_n2_2, gauss_n2_2 ),  
      QuadraturePoint(.25, gauss_n2_2, gauss_n2_1 )   }; 
// ----------------------------------------------------------------------

const QuadratureFormular QuadratureFormular_Q_3(4,3,4,P_QuadratureFormular_Q_3);
// from:  http://mathworld.wolfram.com/Legendre-GaussQuadrature.html
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



/* 
~ script degeneration en csh 

foreach i  (0 1 2)
foreach j  (0 1 2)
 echo "QuadraturePoint( pgauss_n3_$i * pgauss_n3_$j,  gauss_n3_$i, gauss_n3_$j ),"
end
end
*/ 

static const QuadraturePoint P_QuadratureFormular_Q_5[9] = {
       QuadraturePoint( pgauss_n3_0 * pgauss_n3_0,  gauss_n3_0, gauss_n3_0 ),
       QuadraturePoint( pgauss_n3_0 * pgauss_n3_1,  gauss_n3_0, gauss_n3_1 ),
       QuadraturePoint( pgauss_n3_0 * pgauss_n3_2,  gauss_n3_0, gauss_n3_2 ),
       QuadraturePoint( pgauss_n3_1 * pgauss_n3_0,  gauss_n3_1, gauss_n3_0 ),
       QuadraturePoint( pgauss_n3_1 * pgauss_n3_1,  gauss_n3_1, gauss_n3_1 ),
       QuadraturePoint( pgauss_n3_1 * pgauss_n3_2,  gauss_n3_1, gauss_n3_2 ),
       QuadraturePoint( pgauss_n3_2 * pgauss_n3_0,  gauss_n3_2, gauss_n3_0 ),
       QuadraturePoint( pgauss_n3_2 * pgauss_n3_1,  gauss_n3_2, gauss_n3_1 ),
       QuadraturePoint( pgauss_n3_2 * pgauss_n3_2,  gauss_n3_2, gauss_n3_2 ) };

const QuadratureFormular QuadratureFormular_Q_5(4,5,9,P_QuadratureFormular_Q_5);
// ----------------------------------------------------------------------

/* 
~ script degeneration en csh 

foreach i  (0 1 2 3 4)
foreach j  (0 1 2 3 4 )
 echo "QuadraturePoint( pgauss_n5_$i * pgauss_n5_$j,  gauss_n5_$i, gauss_n5_$j ),"
end
end
*/ 
// coef zienkiewicz -- la methode de element finis :  afnor technique page 180
const R g5_2  = 0.906179845938664;
const R g5_1  = 0.538469310105683;
const R pg5_2 = 0.236926885056189;
const R pg5_1 = 0.478628670499366;
const R pg5_0 = 512.0/900.0;

const R gauss_n5_0 = 0.5 - g5_2/2;  
const R gauss_n5_1 = 0.5 - g5_1/2;  
const R gauss_n5_2 = 0.5;
const R gauss_n5_3 = 0.5 + g5_1/2;
const R gauss_n5_4 = 0.5 + g5_2/2;


const R pgauss_n5_0=  pg5_2/2;
const R pgauss_n5_1=  pg5_1/2;
const R pgauss_n5_2=  pg5_0/2;
const R pgauss_n5_3=  pg5_1/2;
const R pgauss_n5_4=  pg5_2/2;


static const QuadraturePoint P_QuadratureFormular_Q_7[25] = {
QuadraturePoint( pgauss_n5_0 * pgauss_n5_0,  gauss_n5_0, gauss_n5_0 ),
QuadraturePoint( pgauss_n5_0 * pgauss_n5_1,  gauss_n5_0, gauss_n5_1 ),
QuadraturePoint( pgauss_n5_0 * pgauss_n5_2,  gauss_n5_0, gauss_n5_2 ),
QuadraturePoint( pgauss_n5_0 * pgauss_n5_3,  gauss_n5_0, gauss_n5_3 ),
QuadraturePoint( pgauss_n5_0 * pgauss_n5_4,  gauss_n5_0, gauss_n5_4 ),
QuadraturePoint( pgauss_n5_1 * pgauss_n5_0,  gauss_n5_1, gauss_n5_0 ),
QuadraturePoint( pgauss_n5_1 * pgauss_n5_1,  gauss_n5_1, gauss_n5_1 ),
QuadraturePoint( pgauss_n5_1 * pgauss_n5_2,  gauss_n5_1, gauss_n5_2 ),
QuadraturePoint( pgauss_n5_1 * pgauss_n5_3,  gauss_n5_1, gauss_n5_3 ),
QuadraturePoint( pgauss_n5_1 * pgauss_n5_4,  gauss_n5_1, gauss_n5_4 ),
QuadraturePoint( pgauss_n5_2 * pgauss_n5_0,  gauss_n5_2, gauss_n5_0 ),
QuadraturePoint( pgauss_n5_2 * pgauss_n5_1,  gauss_n5_2, gauss_n5_1 ),
QuadraturePoint( pgauss_n5_2 * pgauss_n5_2,  gauss_n5_2, gauss_n5_2 ),
QuadraturePoint( pgauss_n5_2 * pgauss_n5_3,  gauss_n5_2, gauss_n5_3 ),
QuadraturePoint( pgauss_n5_2 * pgauss_n5_4,  gauss_n5_2, gauss_n5_4 ),
QuadraturePoint( pgauss_n5_3 * pgauss_n5_0,  gauss_n5_3, gauss_n5_0 ),
QuadraturePoint( pgauss_n5_3 * pgauss_n5_1,  gauss_n5_3, gauss_n5_1 ),
QuadraturePoint( pgauss_n5_3 * pgauss_n5_2,  gauss_n5_3, gauss_n5_2 ),
QuadraturePoint( pgauss_n5_3 * pgauss_n5_3,  gauss_n5_3, gauss_n5_3 ),
QuadraturePoint( pgauss_n5_3 * pgauss_n5_4,  gauss_n5_3, gauss_n5_4 ),
QuadraturePoint( pgauss_n5_4 * pgauss_n5_0,  gauss_n5_4, gauss_n5_0 ),
QuadraturePoint( pgauss_n5_4 * pgauss_n5_1,  gauss_n5_4, gauss_n5_1 ),
QuadraturePoint( pgauss_n5_4 * pgauss_n5_2,  gauss_n5_4, gauss_n5_2 ),
QuadraturePoint( pgauss_n5_4 * pgauss_n5_3,  gauss_n5_4, gauss_n5_3 ),
QuadraturePoint( pgauss_n5_4 * pgauss_n5_4,  gauss_n5_4, gauss_n5_4 )};

const QuadratureFormular QuadratureFormular_Q_7(4,7,25,P_QuadratureFormular_Q_7);

void QuadratureFormular::Verification()
{
  const double tol=1.e-12;
  R err=0;
  // on triangle or Quad 
  R a,b,c,h;
  for (int k=0;k<=exact;  k++)
    {
      R sa(0),sb(0),sc(0);
      for (int j=0;j<n;j++)
	{
	  h = p[j];
	  R2 P = p[j];
	  a = P.x;
	  b = P.y;
	  c = 1.-a-b;
	  R hak(h),hbk(h),hck(h);
	  for (int i=0;i<k;i++) 
	    hak *= a,hbk *= b,hck *= c;
	  sa += hak;
	  sb += hbk;
	  sc += hck;
	  
	}
      
      R se(1),see(1);
      if (on == 3) {
	for (int i=1;i<=k;i++) 
	  se *= (R) i / (R) (i+2);
	see=se;
      }
      else if (on == 4)
	{
	  se = 1.0 / (R) ( k+1);
	  see = k % 2 ? 0 :  2 / (R) ((k+1 )*(k+2)) ;
	}
      err = Max(err,Abs(se-sa));
      err = Max(err,Abs(se-sb));
      err = Max(err,Abs(see-sc));
      if (err>tol)
	if (on == 3) 
	  cerr << "T Ordre= " << k << " 2!k!/(2+k)!= " << se << " " << sa << " " << sb << " " << sc 
	       << " err= " << err << endl;
	else
	  cerr << "Q Ordre= " << k << " 1/(k+1) =" << se << " " << sa << " " << sb << " " 
	       << sc << " == " << see
	       << " err= " << err << endl;

    }
  
  if(err>tol)
    {
      cerr << "Erreur dans la formule d'integration on=" << on << " exact = " << exact 
	   << " Nb Point = " << n << endl;
      throw(ErrorExec("exit",1));
    }
    
}

void QuadratureFormular1d::Check()
{
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
const QuadratureFormular1d QF_GaussLegendre3(5,
                  QuadratureFormular1d::Point(pgauss_n3_0,gauss_n3_0),
                  QuadratureFormular1d::Point(pgauss_n3_1,gauss_n3_1),
                  QuadratureFormular1d::Point(pgauss_n3_2,gauss_n3_2)); 
                    
const QuadratureFormular1d QF_GaussLegendre2(3,
                  QuadratureFormular1d::Point(0.5,gauss_n2_1),
                  QuadratureFormular1d::Point(0.5,gauss_n2_2)); 

const QuadratureFormular1d QF_GaussLegendre1(1,QuadratureFormular1d::Point(1,0.5)); 

const QuadratureFormular1d QF_LumpP1_1D(1,
                  QuadratureFormular1d::Point(0.5,0.),
                  QuadratureFormular1d::Point(0.5,1.)); 


const QuadratureFormular1d QF_GaussLegendre4(7,
					     QuadratureFormular1d::Point(pgauss1_n4_a/2.,(1.+gauss1_n4_a)/2.),
					     QuadratureFormular1d::Point(pgauss1_n4_a/2.,(1.-gauss1_n4_a)/2.),
					     QuadratureFormular1d::Point(pgauss1_n4_b/2.,(1.+gauss1_n4_b)/2.),
					     QuadratureFormular1d::Point(pgauss1_n4_b/2.,(1.-gauss1_n4_b)/2.)
					     ); 

const QuadratureFormular1d QF_GaussLegendre5(9,
					     QuadratureFormular1d::Point(pgauss1_n5_0/2.,0.5),
					     QuadratureFormular1d::Point(pgauss1_n5_a/2.,(1.+gauss1_n5_a)/2),
					     QuadratureFormular1d::Point(pgauss1_n5_a/2.,(1.-gauss1_n5_a)/2.),
					     QuadratureFormular1d::Point(pgauss1_n5_b/2.,(1.+gauss1_n5_b)/2.),
					     QuadratureFormular1d::Point(pgauss1_n5_b/2.,(1.-gauss1_n5_b)/2.)
					     ); 

}
