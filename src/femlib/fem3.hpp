#ifndef FEM_3_H_
#define FEM_3_H_

// The class R3
class R3: public R2 {
  friend ostream& operator <<(ostream& f, const R3 & P )
  { f << P.x << ' ' << P.y << ' ' << P.z   ; return f; }
  friend istream& operator >>(istream& f,  R3 & P)
  { f >>  P.x >>  P.y >> P.z  ; return f; }

public:  
  R z;
 
  R3 () :z(0) {};
  R3 (R a,R b,R c):R2(a,b),z(c)  {}
  R3 (R3 A,R3 B) :R2(B.x-A.x,B.y-A.y),z(B.z-A.z)  {}
  R3   operator+(R3 P)const   {return R3(x+P.x,y+P.y,z+P.z);}
  R3   operator+=(R3 P)  {x += P.x;y += P.y;z += P.z;return *this;}
  R3   operator-(R3 P)const   {return R3(x-P.x,y-P.y,z -P.z);}
  R3   operator-=(R3 P) {x -= P.x;y -= P.y;z -= P.z;return *this;}
  R3   operator-()const  {return R3(-x,-y,-z);}
  R3   operator+()const  {return *this;}
  R   operator,(R3 P)const  {return  x*P.x+y*P.y+z*P.z;} // produit scalaire
  R3   operator^(R3 P)const {return R3(y*P.z-z*P.y ,P.x*z-x*P.z, x*P.y-y*P.x);} // produit mixte
  R3   operator*(R c)const {return R3(x*c,y*c,z*c);}
  R3   operator/(R c)const {return R3(x/c,y/c,z/c);}
  R  &  operator[](int i){ return (&x)[i];}
//  R2   perp() {return R2(-y,x);} // the perpendicular
  friend R3 operator*(R c,R3 P) {return P*c;}
  friend R3 operator/(R c,R3 P) {return P/c;}
};

//inline void MoveTo(R2 P) { rmoveto((float) P.x,(float)P.y);}
//inline void LineTo(R2 P) { rlineto((float)P.x,(float)P.y);}

 
inline R Norme2_2(const R3 & A){ return (A,A);}
inline R Norme2(const R3 & A){ return sqrt((A,A));}
inline R Norme_infty(const R3 & A){return Max(Abs(A.x),Abs(A.y),Abs(A.z));}
inline R Area2(const R3 A,const R3 B,const R3 C){return Norme2((B-A)^(C-A));}

#endif
