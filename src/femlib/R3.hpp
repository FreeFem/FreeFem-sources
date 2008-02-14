#ifndef R3_HPP
#define  R3_HPP

#include "R1.hpp" 
#include "R2.hpp" 

class R3  {

public:  
  typedef double R;

  static const int d=3;
  
  R x,y,z;
 
  R3 () :x(0),y(0),z(0) {};
  R3 (R a,R b,R c):x(a),y(b),z(c)  {}
  R3 (R2 P2):x(P2.x),y(P2.y),z(0)  {}
  R3 (R2 P2,R zz):x(P2.x),y(P2.y),z(zz)  {}
  R3 (const R3 & A,const R3 & B) :x(B.x-A.x),y(B.y-A.y),z(B.z-A.z)  {}
  R3 & operator=(R2 P2) {x=P2.x;y=P2.y;z=0;return *this;}
  R3   operator+(R3 P)const   {return R3(x+P.x,y+P.y,z+P.z);}
  R3   operator-(R3 P)const   {return R3(x-P.x,y-P.y,z-P.z);}

  R3 & operator+=(R3 P)  {x += P.x;y += P.y;z += P.z;return *this;}
  R3 & operator-=(R3 P) {x -= P.x;y -= P.y;z -= P.z;return *this;}
  R3 & operator*=(R c) {x *= c;y *= c;z *= c;return *this;}
  R3 & operator/=(R c) {x /= c;y /= c;z /= c;return *this;}

  R3   operator-()const  {return R3(-x,-y,-z);}
  R3   operator+()const  {return *this;}

  R    operator,(R3 P)const  {return  x*P.x+y*P.y+z*P.z;} // produit scalaire
  R3   operator^(R3 P)const {return R3(y*P.z-z*P.y ,P.x*z-x*P.z, x*P.y-y*P.x);} // produit vectoreil

  R3   operator*(R c)const {return R3(x*c,y*c,z*c);}
  R3   operator/(R c)const {return R3(x/c,y/c,z/c);}

  R  &  operator[](int i){ return (&x)[i];}
  const R  &  operator[](int i) const { return (&x)[i];}

  friend R3 operator*(R c,R3 P) {return P*c;}
  friend R3 operator/(R c,R3 P) {return P/c;}

  R norme() const { return std::sqrt(x*x+y*y+z*z);}
  R norme2() const { return (x*x+y*y+z*z);}
  R sum() const { return x+y+z;}
  R3 Bary(const R3 P[d+1]) const { return (1-x-y-z)*P[0]+x*P[1]+y*P[2]+z*P[3];}  // add FH 
  R3 Bary(const R3 **P ) const { return (1-x-y-z)*(*P[0])+x*(*P[1])+y*(*P[2])+z*(*P[3]);}  // add FH 
  friend ostream& operator <<(ostream& f, const R3 & P )
  { f << P.x << ' ' << P.y << ' ' << P.z   ; return f; }
  friend istream& operator >>(istream& f,  R3 & P)
  { f >>  P.x >>  P.y >>  P.z  ; return f; }

  friend R det(R3 A,R3 B, R3 C) 
  {
    R  s=1.;
    if(abs(A.x)<abs(B.x)) Exchange(A,B),s = -s;
    if(abs(A.x)<abs(C.x)) Exchange(A,C),s = -s;
    if(abs(A.x)>1e-50)
      {
	s *= A.x;
	A.y /= A.x; A.z /= A.x;
	B.y  -=  A.y*B.x  ;   B.z  -=  A.z*B.x ;
	C.y  -=  A.y*C.x  ;   C.z  -=  A.z*C.x ;
	return s* ( B.y*C.z - B.z*C.y) ;
      }
    else return 0.   ;
  }
  R2 p2() const { return R2(x,y);}
  R1 p1() const { return R1(x);}

  friend R  det(R3 A,R3 B, R3 C, R3 D) { return det(R3(A,B),R3(A,C),R3(A,D));}
  static const R3 KHat[d+1];
  

};

// inline R3 Min(const R3 & A,const R3 &B){ return R3(min(A.x,B.x),min(A.y,B.y),min(A.z,B.z));}
// inline R3 Max(const R3 & A,const R3 &B){ return R3(max(A.x,B.x),max(A.y,B.y),max(A.z,B.z));}

#endif

