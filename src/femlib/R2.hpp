#ifndef R2_HPP
#define  R2_HPP
#include <cmath>
#include <cstdlib>
#include <iostream>
// Definition de la class R2 
//  sans compilation separe toute les fonctions 
// sous defini dans ce R2.hpp avec des inline 
//  
// definition R (les nombres reals)
// remarque la fonction abort est defini dans 
// #include <cstdlib> 

// The class R2
class R2 {
public:  
  typedef double R;
  static const int d=2;

  R x,y;  // declaration de membre 
  // les 3 constructeurs ---
  R2 () :x(0.),y(0.) {} // rappel : x(0), y(0)  sont initialiser via le constructeur de double 
  R2 (R a,R b):x(a),y(b)  {}
  R2 (const R * a):x(a[0]),y(a[1])  {}
  R2 ( R * a):x(a[0]),y(a[1])  {}

  R2 (const R2 & a,const R2 & b):x(b.x-a.x),y(b.y-a.y)  {}
  static  R2 diag(R a){ return R2(a,a);}
  // le constucteur par defaut est inutile
  //R2 (const R2 & a) :x(a.x),y(a.y) {} 

  // rappel les operator definis dans une class on un parametre
  // cache qui est la classe elle meme (*this)

  // les operateurs affectation
  //  operateur affection (*this) = P est inutil par defaut il fait le travail correctement
  //R2 &  operator=(const R2 & P)  {x = P.x;y = P.y;return *this;}
  // les autre operoteur affectations
  R2 &  operator+=(const R2 & P)  {x += P.x;y += P.y;return *this;}
  R2 &  operator-=(const R2 & P) {x -= P.x;y -= P.y;return *this;}
  R2 &  operator*=(R a) {x *= a;y *= a;return *this;}
  R2 &  operator/=(R a) {x /= a;y /= a;return *this;}

  // operateur binaire + - * , ^ /
  R2   operator+(const R2 & P)const   {return R2(x+P.x,y+P.y);}
  R2   operator-(const R2 & P)const   {return R2(x-P.x,y-P.y);}
  R    operator,(const R2 & P)const  {return  x*P.x+y*P.y;} // produit scalaire
  R    operator^(const R2 & P)const {return  x*P.y-y*P.x;} // produit mixte
  R2   operator*(R c)const {return R2(x*c,y*c);}
  R2   operator/(R c)const {return R2(x/c,y/c);}
  // operateur unaire 
  R2   operator-()const  {return R2(-x,-y);} 
  R2   operator+()const  {return *this;}
  // un methode
  R2   perp() const {return R2(-y,x);} // la perpendiculaire
  R sum() const { return x+y;}

  // les operators  tableau
  // version qui peut modifie la class  via l'adresse de x ou y 
  R  &  operator[](int i){ return (&x)[i];}
  const R  &  operator[](int i) const { return (&x)[i];}

  R  X() const {return x;}
  R  Y() const {return y;}
  R  Z() const {return 0.;}

  R norme() const { return std::sqrt(x*x+y*y);}
  R norme2() const { return (x*x+y*y);}
  R2 Bary(R2 P[d+1]) const { return (1-x-y)*P[0]+x*P[1]+y*P[2];}  // add FH 
  R2 Bary(const R2 *const *const P ) const { return (1-x-y)*(*P[0])+x*(*P[1])+y*(*P[2]);}  // add FH 
friend  R2 operator*(R c,const R2 & P) {return P*c;} 
friend  R2 perp(const R2 & P) { return R2(-P.y,P.x) ; }
//inline R2 Perp(const R2 & P) { return P.perp(); }  // autre ecriture  de la fonction perp
friend R  det(const R2 & A,const R2 & B,const R2 &C) { return R2(A,B)^R2(A,C);}

friend  std::ostream& operator <<(std::ostream& f, const R2 & P )
       { f << P.x << ' ' << P.y   ; return f; }
friend  std::istream& operator >>(std::istream& f,  R2 & P)
       { f >>  P.x >>  P.y  ; return f; }
  static const R2 KHat[d+1];
};

inline R2 Minc(const R2 & A,const R2 &B){ return R2(min(A.x,B.x),min(A.y,B.y));}
inline R2 Maxc(const R2 & A,const R2 &B){ return R2(max(A.x,B.x),max(A.y,B.y));}
inline double Norme_infty(const R2 & A){return max(std::fabs(A.x),std::fabs(A.y));}
inline double Norme2_2(const R2 & A){ return (A,A);}
inline double  Norme2(const R2 & A){ return sqrt((A,A));}

#endif

