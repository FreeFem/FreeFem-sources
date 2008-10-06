// ORIG-DATE:     Dec 2007
// -*- Mode : c++ -*-
//
// SUMMARY  :  Model of $\mathbb{R}^1$    
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

 Thank to the ARN ()  FF2A3 grant
 ref:ANR-07-CIS7-002-01 
 */


#ifndef R1_HPP
#define  R1_HPP
#include <cmath>
#include <cstdlib>
#include <iostream>
typedef double R;

//  R0 R1, R2 , R3 to be uniforme. 
// The class R1
class R0 {
public:  
  typedef double R;
  static const int d=0;
  R0(){}
};

class R1 {
public:  
  typedef double R;
  static const int d=1;

  R x;  // declaration de membre 
  // les 3 constructeurs ---
  R1 () :x(0.){} // rappel : x(0)  sont initialiser via le constructeur de double 
  R1 (R a):x(a)  {}
  R1 (const R1 & a,const R1 & b):x(b.x-a.x)  {}
  static  R1 diag(R a){ return R1(a);}

  // rappel les operator definis dans une class on un parametre
  // cache qui est la classe elle meme (*this)

  // les operateurs affectation
  //  operateur affection (*this) = P est inutil par defaut il fait le travail correctement
  // les autre operoteur affectations
  R1 &  operator+=(const R1 & P)  {x += P.x;return *this;}
  R1 &  operator-=(const R1 & P) {x -= P.x;return *this;}
  R1 &  operator*=(R a) {x *= a;return *this;}
  R1 &  operator/=(R a) {x /= a;return *this;}
  // operateur binaire + - * , ^ /
  R1   operator+(const R1 & P)const   {return R1(x+P.x);}
  R1   operator-(const R1 & P)const   {return R1(x-P.x);}
  R    operator,(const R1 & P)const  {return  x*P.x;} // produit scalaire
  R1   operator*(R c)const {return R1(x*c);}
  R1   operator/(R c)const {return R1(x/c);}
  // operateur unaire 
  R1   operator-()const  {return R1(-x);} 
  R1   operator+()const  {return *this;}
  // un methode
  R sum() const { return x;}
  R * toBary(R * b) const { b[0]=1-x;b[1]=x;return b;}
  // les operators  tableaux
  // version qui peut modifie la class  via l'adresse de x
  R  &  operator[](int i){ return (&x)[i];}
  const R  &  operator[](int i) const { return (&x)[i];}


  R norme() const { return std::sqrt(x*x);}
  R norme2() const { return (x*x);}
  R1 Bary(R1 P[d+1]) const { return (1-x)*P[0]+x*P[1];}  // add FH 
  R1 Bary(const R1 *const *const P ) const { return (1-x)*(*P[0])+x*(*P[1]);}  // add FH 
friend  R1 operator*(R c,const R1 & P) {return P*c;} 

friend  std::ostream& operator <<(std::ostream& f, const R1 & P )
       { f << P.x   ; return f; }
friend  std::istream& operator >>(std::istream& f,  R1 & P)
       { f >>  P.x  ; return f; }
  static const R1 KHat[d+1];
};

inline R1 Minc(const R1 & A,const R1 &B){ return R1(min(A.x,B.x));}
inline R1 Maxc(const R1 & A,const R1 &B){ return R1(max(A.x,B.x));}
inline double Norme_infty(const R1 & A){return std::fabs(A.x);}
inline double Norme2_2(const R1 & A){ return A.x*A.x;}
inline double Norme2(const R1 & A){ return std::fabs(A.x);}


#endif
