#ifndef _QuadratureFormular_h
#define _QuadratureFormular_h

namespace Fem2D {
class QuadratureFormular;

class QuadraturePoint : public R2{
  friend ostream& operator <<(ostream& , QuadraturePoint & );
public:
const R a;
//const R2 c;
  QuadraturePoint(R aa,R2 xx): a(aa),R2(xx) {}
  QuadraturePoint(R aa,R xx,R yy): a(aa),R2(xx,yy) {}
  operator R() const {return a;}
//  operator R2()   const {return c;}
};


class QuadratureFormular {
 friend  ostream& operator <<(ostream& , const QuadratureFormular & ) ;

public:
  const int n;                // nombre de point d'integration
  const int exact;            // exact
  const int on;               // on = 3 => triangle, on=4 => quad
  const QuadraturePoint *p;  // les point d'integration 
// -- les fonctions ------------------
  void Verification(); // for verification 
  QuadratureFormular (int o,int e,int NbOfNodes,QuadraturePoint *pp):on(o),exact(e), n(NbOfNodes),p(pp)
    {Verification();}
  QuadratureFormular (int o,int e,int NbOfNodes,const QuadraturePoint *pp):on(o),exact(e),n(NbOfNodes),p(pp)
    {Verification();}

const QuadraturePoint & operator [](int i) const  {return p[i];} 
const QuadraturePoint & operator ()(int i) const {return p[i];}

private:
 QuadratureFormular(const QuadratureFormular &)
   :n(0),exact(0),on(0),p(0){throwassert(0);}
 void operator=(const QuadratureFormular &)
   {throwassert(0);}
 QuadratureFormular()
   :n(0),exact(0),on(0),p(0){throwassert(0);}
 

};

ostream& operator <<(ostream& , const QuadratureFormular & ) ;
ostream& operator <<(ostream& , QuadraturePoint & );


class QuadratureFormular1d {public:
  const int n;
  class Point { public: R a,x ; Point(R aa=0,R xx=0): a(aa),x(xx) {} } *p;
  QuadratureFormular1d(Point p0,Point p1,Point p2) : n(3),p(new Point[3]) { p[0]=p0,p[1]=p1,p[2]=p2;}
  QuadratureFormular1d(Point p0,Point p1) : n(2),p(new Point[2]) { p[0]=p0,p[1]=p1;}
  QuadratureFormular1d(Point p0) : n(1),p(new Point[1]) { p[0]=p0;}
  ~QuadratureFormular1d(){ delete [] p;}
 Point operator[](int i) const { return p[i];}  
 private:
 // pas operator de copie pour les formules d'intergation
 QuadratureFormular1d()
  : n(0){throwassert(0);};
 QuadratureFormular1d(const QuadratureFormular1d &)
  :n(0) {throwassert(0);}
 void operator=(const QuadratureFormular1d &){throwassert(0);}
} ;

extern const QuadratureFormular1d QF_GaussLegendre1; 
extern const QuadratureFormular1d QF_GaussLegendre2; 
extern const QuadratureFormular1d QF_GaussLegendre3; 

extern const QuadratureFormular QuadratureFormular_T_1;
extern const QuadratureFormular QuadratureFormular_T_1lump;
extern const QuadratureFormular QuadratureFormular_T_2;
extern const QuadratureFormular QuadratureFormular_T_5;
extern const QuadratureFormular QuadratureFormular_T_2_4P1;
//extern const QuadratureFormular QuadratureFormular_T_7;
extern const QuadratureFormular QuadratureFormular_Q_1;
extern const QuadratureFormular QuadratureFormular_Q_3;
extern const QuadratureFormular QuadratureFormular_Q_5;
extern const QuadratureFormular QuadratureFormular_Q_7;

ostream& operator <<(ostream& f,const  QuadraturePoint & p) ;
ostream& operator <<(ostream& f, const QuadratureFormular & fi) ;
      
}

#endif
