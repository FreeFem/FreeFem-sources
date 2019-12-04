//
// RZ.h
// MeshGenQA
//
// Created by Jean-Marie Mirebeau on 18/06/11.
// Copyright 2011 UPMC. All rights reserved.
//
#ifndef RZ_H
#define RZ_H

#include <iostream>
#include <vector>
#include "math.h"
#include "float.h"
#include "BasicMath.h"
using namespace std;

/*
 * Strongly inspired by R2.hpp from F. Hecht
 *
 * The class RZ contains a double and an integer. Elements of this class are ordered
 * lexicographically.
 *
 * The template classes Bidim and TriDim describe points with coordinates in a ring (usually integer
 * or real). The classes sym2 and sym3 describe symmetric matrices with real entries.
 *
 * Lattice basis reduction is implemented under the following names :
 * bool sym2::reduced_basis(Z2 Basis[2]) const;
 * bool sym3::reduced_basis(Z3 Basis[3]) const;
 *
 */

class RZ;

template< class ring >
class BiDim;
typedef BiDim< double > R2;
typedef BiDim< int > Z2;

template< class ring >
class TriDim;
typedef TriDim< double > R3;
typedef TriDim< int > Z3;

R2 ZdtoRd(const Z2 &P);
R3 ZdtoRd(const Z3 &P);

class sym2;
class sym3;

class Metric2;
class FctMetric2;
class Metric3;
class FctMetric3;

// print conventions at the end of the file
template< class ring >
inline ostream &operator<<(ostream &f, const BiDim< ring > &P);
template< class ring >
inline ostream &operator<<(ostream &f, const TriDim< ring > &P);
inline ostream &operator<<(ostream &f, const sym2 &S);
inline ostream &operator<<(ostream &f, const sym3 &S);

/**************** operators *************/
// boost implementation is presumably better, but you need to install boost first.
#ifdef USE_BOOST_OPERATORS
#include <boost/operators.hpp>
using boost::equality_comparable;
using boost::totally_ordered;
template< class V, class K >
class vector_space : boost::additive< V >, boost::multiplicative< V, K > {};

#else
// Operators for a vector space V on a field K. (or a module on a ring, but be careful with
// divisions.) K elements are passed by value, since we assume that it is a machine type (e.g.
// double, int)
template< class V, class K >
class vector_space {
 public:
  friend V operator+(V u, const V &v) { return u += v; }

  friend V operator-(V u, const V &v) { return u -= v; }

  friend V operator*(V u, K k) { return u *= k; }

  friend V operator*(K k, V u) { return u *= k; }

  friend V operator/(V u, K k) { return u /= k; }
};

template< class A >
class totally_ordered {
  friend bool operator>(const A &a, const A &b) { return b < a; }

  friend bool operator<=(const A &a, const A &b) { return !(b < a); }

  friend bool operator>=(const A &a, const A &b) { return !(a < b); }

  friend bool operator!=(const A &a, const A &b) { return !(a == b); }
};

template< class A >
class equality_comparable {
  friend bool operator!=(const A &a, const A &b) { return !(a == b); }
};

#endif

/**************** classe RZ (distance, number - index) *************/

class RZ : totally_ordered< RZ > {    // could have used std::pair but too late
 public:
  double distance;
  int number;

  RZ( ) : distance(-DBL_MAX), number(INT_MIN) {}

  RZ(double dist, int num) : distance(dist), number(num) {}

  RZ(const RZ &L) : distance(L.distance), number(L.number) {}

  RZ &operator=(const RZ &L) {
    distance = L.distance;
    number = L.number;
    return *this;
  }

  bool operator==(const RZ &L) const { return distance == L.distance && number == L.number; }

  bool operator<(const RZ &L) const {
    return distance < L.distance || (distance == L.distance && number < L.number);
  }

  RZ MIN( ) { return RZ(-DBL_MAX, INT_MIN); }

  RZ MAX( ) { return RZ(DBL_MAX, INT_MAX); }
};

inline ostream &operator<<(ostream &f, const RZ &P) {
  f << " " << P.distance << " " << P.number;
  return f;
}

inline ostream_math operator<<(ostream_math f, const RZ &P) {
  f << "{" << P.distance << "," << P.number << "}";
  return f;
}

/***************** classes 2d ****************/

template< class ring >
class BiDim : vector_space< BiDim< ring >, ring >, totally_ordered< BiDim< ring > > {
 public:
  ring x, y;
  BiDim( ) : x(0), y(0){};
  BiDim(ring X, ring Y) : x(X), y(Y){};
  BiDim &operator*=(ring c) {
    x *= c;
    y *= c;
    return *this;
  }

  BiDim &operator/=(ring c) {
    const ring d = 1 / c;
    assert_msg(c != ring(0) && d != ring(0), name << " error : division by " << c);
    x *= d;
    y *= d;
    return *this;
  }

  BiDim &operator+=(const BiDim &P) {
    x += P.x;
    y += P.y;
    return *this;
  }

  BiDim &operator-=(const BiDim &P) {
    x -= P.x;
    y -= P.y;
    return *this;
  }

  bool operator==(const BiDim &P) const { return x == P.x && y == P.y; }

  bool operator<(const BiDim &P) const { return x < P.x || (x == P.x && y < P.y); }

  const BiDim operator-( ) const { return BiDim(-x, -y); }

  const BiDim lin_sum(const BiDim &P, const BiDim &Q) const {
    return BiDim(x * P.x + y * Q.x, x * P.y + y * Q.y);
  };    // produit [P,Q] *this
  const TriDim< ring > lin_sum(const TriDim< ring > &P, const TriDim< ring > &Q) const {
    return TriDim< ring >(x * P.x + y * Q.x, x * P.y + y * Q.y, x * P.z + y * Q.z);
  }

  const BiDim lin_solve(const BiDim &P, const BiDim &Q)
    const;    // coefficients de *this dans la base (P,Q) : [P,Q]^{-1} *this.

  ring prod( ) const { return x * y; }

  ring min( ) const { return x < y ? x : y; }

  ring norm2( ) const { return x * x + y * y; }

  double norm( ) const { return sqrt(norm2( )); }

  double norm1( ) const { return fabs(x) + fabs(y); }

  ring scal(const BiDim &Q) const { return x * Q.x + y * Q.y; }

  bool isPositive( ) const { return x >= ring(0) && y >= ring(0); }

  // bool isStrictlyPositive() const {return x>ring(0) && y>ring(0);}
  bool isZero( ) const { return x == ring(0) && y == ring(0); }

  ring operator[](int i) const {
    assert_msg(0 <= i && i <= 1, name << "::[] error : out of bounds");
    if (i == 0) {
      return x;
    } else {
      return y;
    }
  }

  ring &operator[](int i) {
    assert_msg(0 <= i && i <= 1, name << "::[] error : out of bounds");
    if (i == 0) {
      return x;
    } else {
      return y;
    }
  }

  static const int dim = 2;
  static const string name;
  static const BiDim NABiDim;
  BiDim &sort( ) {
    if (x > y) {
      const ring t = x;
      x = y;
      y = t;
    }

    return *this;
  }
};

template< class ring >
ring det(const BiDim< ring > &P, const BiDim< ring > &Q) {
  return P.x * Q.y - P.y * Q.x;
}

inline R2 ZdtoRd(const Z2 &P) { return R2(P.x, P.y); }

template< class ring >
const BiDim< ring > BiDim< ring >::lin_solve(const BiDim< ring > &P, const BiDim< ring > &Q)
  const {    // coefficients de *this dans la base (P,Q) : [P,Q]^{-1} *this.
  ring Det = P.x * Q.y - P.y * Q.x;

  if (Det == 0) {
    cout << name << "::lin_solve error : vectors are collinear " << P << "; " << Q << endl;
    return NABiDim;
  }

  ring invDet = 1 / Det;
  if (invDet == 0) {
    cout << name << "::lin_solve error : determinant is not invertible " << Det << "; " << P << "; "
         << Q << endl;
    return NABiDim;
  }

  return BiDim((x * Q.y - y * Q.x) * invDet, (-x * P.y + y * P.x) * invDet);
};

// matrices symétriques
class sym2 : vector_space< sym2, double > {
 public:
  double xx, xy, yy;
  sym2( ) : xx(1), xy(0), yy(1) {}

  sym2(double XX, double XY, double YY) : xx(XX), xy(XY), yy(YY) {}

  sym2(double lambda, double mu, const R2 &P);
  sym2(const sym2 &S) : xx(S.xx), xy(S.xy), yy(S.yy) {}

  explicit sym2(double s) : xx(s), xy(0), yy(s){};

  template< class ring >
  double norm2(const BiDim< ring > &P) const {
    return xx * (P.x * P.x) + 2 * xy * (P.x * P.y) + yy * (P.y * P.y);
  }

  template< class ring >
  double norm(const BiDim< ring > &P) const {
    return sqrt(norm2(P));
  }

  template< class ring >
  double scal(const BiDim< ring > &P, const BiDim< ring > &Q) const {
    return xx * (P.x * Q.x) + xy * (P.x * Q.y + P.y * Q.x) + yy * (P.y * Q.y);
  }

  template< class ring >
  double cos2(const BiDim< ring > &P, const BiDim< ring > &Q) const {
    return square(scal(P, Q)) / (norm2(P) * norm2(Q));
  }

  template< class ring >
  double cos(const BiDim< ring > &P, const BiDim< ring > &Q) const {
    return scal(P, Q) / (norm(P) * norm(Q));
  }

  double det( ) const { return xx * yy - xy * xy; }

  double trace( ) const { return xx + yy; }

  bool isDiagonal( ) const { return xy == 0; }

  bool isPositiveDefinite( ) const { return det( ) > 0 && trace( ) > 0; }

  void eigen(double E[2]) const {    // petite puis grande vap
    double hDiff = sqrt(square(xx - yy) / 4 + square(xy));
    double hSum = (xx + yy) / 2;

    E[0] = hSum - hDiff;
    E[1] = hSum + hDiff;
  }

  R2 eigensys(double E[2]) const;
  R2 eigensys( ) const {
    double E[2];
    return eigensys(E);
  }

  double norm( ) const {
    double E[2];
    eigen(E);
    return max(-E[0], E[1]);
  }

  double invNorm( ) const { return inverse( ).norm( ); }

  // sym2 &  operator+=(double mu)  {xx += mu;yy += mu;return *this;}

  sym2 &operator+=(sym2 S) {
    xx += S.xx;
    xy += S.xy;
    yy += S.yy;
    return *this;
  }

  sym2 &operator-=(sym2 S) {
    xx -= S.xx;
    xy -= S.xy;
    yy -= S.yy;
    return *this;
  }

  sym2 &operator*=(double lambda) {
    xx *= lambda;
    xy *= lambda;
    yy *= lambda;
    return *this;
  }

  sym2 &operator/=(double lambda) {
    assert_msg(lambda != 0, "sym2 error : division by 0");
    return operator*=(1 / lambda);
  }

  R2 solve(const R2 &P) const;
  sym2 comatrix( ) const { return sym2(yy, -xy, xx); }

  sym2 inverse( ) const;
  sym2 sqrtSym( ) const;
  sym2 conjugate(const R2 &P, const R2 &Q) const { return sym2(norm2(P), scal(P, Q), norm2(Q)); }

  sym2 exaggerate( ) const;
  sym2 tame( ) const;

  void reduced_basis(Z2 Basis[2]) const;
  R2 operator*(const R2 &P) const { return R2(xx * P.x + xy * P.y, xy * P.x + yy * P.y); }
};

inline sym2::sym2(double lambda, double mu, const R2 &P) {
  const double normP = P.norm( );    // sqrt(P.x*P.x+P.y*P.y); //norme euclidienne...

  if (normP == 0.) {
    xx = yy = sqrt(fabs(lambda * mu));
    xy = 0;
    return;
  }

  ;
  R2 Q = P / normP;
  const double nu = lambda - mu;
  xx = nu * Q.x * Q.x + mu;
  xy = nu * Q.x * Q.y;
  yy = nu * Q.y * Q.y + mu;
}

inline R2 sym2::eigensys(double E[2])
  const {    // petite puis grande vap. vecteur propre associé à la petite valeur propre
  eigen(E);
  const double div = E[0] * E[0] - E[1] * E[1];
  if (div == 0.) {
    return R2(0, 0);    // matrice isotrope : on renvoit 0;
  }

  double c = (E[0] * xx - E[1] * yy) / div;
  c = c >= 0 ? sqrt(c) : 0;    // normalement toujours positif
  double s = (E[0] * yy - E[1] * xx) / div;
  s = s >= 0 ? sqrt(s) : 0;
  s = (E[0] - E[1]) * xy > 0 ? s : -s;
  return R2(c, s);
}

inline R2 sym2::solve(const R2 &P) const {
  const double Det = det( );

  assert_msg(Det != 0, "sym2::solve error : determinant is zero. " << *this);
  return R2(yy * P.x - xy * P.y, -xy * P.x + xx * P.y) / Det;
}

inline sym2 sym2::inverse( ) const {
  double Det = det( );

  assert_msg(Det != 0, "sym2::inverse error : determinant is zero. " << *this);

  return sym2(yy / Det, -xy / Det, xx / Det);
}

inline sym2 sym2::sqrtSym( ) const {
  double E[2];

  eigen(E);
  sym2 Sq(*this);

  assert_msg(E[0] >= 0 && E[1] > 0,
             "sym2::sqrtSym error : non positive definite matrix. " << *this);
  double a = 1 / (sqrt(E[0]) + sqrt(E[1]));
  double b = a * sqrt(E[0] * E[1]);
  Sq *= a;
  Sq += sym2(b);
  return Sq;
}

inline sym2 sym2::exaggerate( ) const {    // returns a symmetric matrix with the same eigenvectors,
                                           // but a slightly increased (up to 4x) largest eigenvalue
  double E[2];

  eigen(E);
  sym2 Ex(*this);

  assert(E[0] >= 0 && E[1] > 0);
  if (E[0] == E[1]) {
    return Ex;
  }

  const double diff = E[1] - E[0], E1n = E[1] * square(2 / (1 + E[0] / E[1]));
  const double a = E[0] * (E[1] - E1n) / diff, b = (E1n - E[0]) / diff;
  return sym2(a) + b * Ex;
}

inline sym2 sym2::tame( )
  const {    // returns a symmetric matrix with the same eigenvectors, but with the largest
             // eigenvalue divided by 4 or equal to the smallest
  double E[2];

  eigen(E);
  sym2 Tm(*this);

  assert(E[0] >= 0 && E[1] > 0);
  if (E[0] == E[1]) {
    return Tm;
  }

  const double diff = E[1] - E[0], E1n = max(E[1] / square(2 / (1 + E[0] / E[1])), E[0]);
  const double a = E[0] * (E[1] - E1n) / diff, b = (E1n - E[0]) / diff;
  return sym2(a) + b * Tm;
}

inline void sym2::reduced_basis(Z2 Basis[2]) const {
  Z2 &b0 = Basis[0], &b1 = Basis[1];

  b0 = Z2(1, 0);
  b1 = Z2(0, 1);
  double n0 = norm2(b0), n1 = norm2(b1);    // Squared norms

  if (n1 < n0) {
    swap(b0, b1);
    swap(n0, n1);
  }    // Sorting the canonical basis by norm

  if (isDiagonal( )) {
    return;
  }

  while (true) {
    b1 -= round(scal(b0, b1) / n0) * b0;
    n1 = norm2(b1);
    if (n0 <= n1) {
      break;
    }

    swap(b0, b1);
    swap(n0, n1);
  }
}

/***************** classes 3d ****************/

template< class ring >
class TriDim : vector_space< TriDim< ring >, ring >, totally_ordered< TriDim< ring > > {
 public:
  ring x, y, z;
  TriDim( ) : x( ), y( ), z( ){};
  TriDim(ring X, ring Y, ring Z) : x(X), y(Y), z(Z){};
  explicit TriDim(int i) : x( ), y( ), z( ) {
    assert_msg(0 <= i && i <= 2, "TriDim::Tridim(int i) error : out of bounds " << i);
    if (i == 0) {
      x = 1;
    } else if (i == 1) {
      y = 1;
    } else {
      z = 1;
    }
  }

  TriDim &operator+=(const TriDim &P) {
    x += P.x;
    y += P.y;
    z += P.z;
    return *this;
  }

  TriDim &operator-=(const TriDim &P) {
    x -= P.x;
    y -= P.y;
    z -= P.z;
    return *this;
  }

  TriDim &operator*=(ring c) {
    x *= c;
    y *= c;
    z *= c;
    return *this;
  }

  TriDim &operator/=(ring c) {
    ring d = 1 / c;
    assert_msg(c != ring(0) && d != ring(0), name << " error : division by " << c);
    return operator*=(d);
  }

  TriDim operator-( ) const { return TriDim(-x, -y, -z); }

  bool operator==(const Z3 &P) const { return x == P.x && y == P.y && z == P.z; }

  TriDim operator^(const TriDim &P) const {
    return TriDim(y * P.z - z * P.y, z * P.x - x * P.z, x * P.y - y * P.x);
  }

  bool operator<(const TriDim &P) const {
    return x < P.x || (x == P.x && (y < P.y || (y == P.y && z < P.z)));
  }

  TriDim lin_sum(const TriDim &P, const TriDim &Q, const TriDim &R) const    // [P,Q,R]*this
  {
    return TriDim(x * P.x + y * Q.x + z * R.x, x * P.y + y * Q.y + z * R.y,
                  x * P.z + y * Q.z + z * R.z);
  }

  TriDim lin_solve(const TriDim &P, const TriDim &Q, const TriDim &R) const;    // [P,Q,R]^{-1}*this

  int prod( ) const { return x * y * z; }

  ring min( ) const { return x < y ? (x < z ? x : z) : (y < z ? y : z); }

  ring max( ) const { return x > y ? (x > z ? x : z) : (y > z ? y : z); }

  ring operator[](int i) const {
    assert_msg(0 <= i && i < 3, name << "::[] error : out of bounds");
    if (i == 0) {
      return x;
    } else if (i == 1) {
      return y;
    } else {
      return z;
    }
  }

  ring &operator[](int i) {
    assert_msg(0 <= i && i < 3, name << "::[] error : out of bounds");
    if (i == 0) {
      return x;
    } else if (i == 1) {
      return y;
    } else {
      return z;
    }
  }

  bool isZero( ) const { return x == ring(0) && y == ring(0) && z == ring(0); }

  bool isPositive( ) const { return x >= ring(0) && y >= ring(0) && z >= ring(0); }

  ring norm2( ) const { return x * x + y * y + z * z; }

  double norm( ) const { return sqrt(norm2( )); }

  double norm1( ) const { return fabs(x) + fabs(y) + fabs(z); }

  ring scal(const TriDim &Q) const { return x * Q.x + y * Q.y + z * Q.z; }

  static const int dim = 3;
  static const string name;
};

inline R3 ZdtoRd(const Z3 &P) { return R3(P.x, P.y, P.z); }

template< class ring >
inline double determinant(TriDim< ring > &P, TriDim< ring > &Q, TriDim< ring > &R) {
  return P.x * Q.y * R.z + P.y * Q.z * R.x + P.z * Q.x * R.y - P.z * Q.y * R.x - P.y * Q.x * R.z -
         P.x * Q.z * R.y;
}

template< class ring >
TriDim< ring > TriDim< ring >::lin_solve(const TriDim< ring > &P, const TriDim< ring > &Q,
                                         const TriDim< ring > &R) const {
  const ring Det = P.x * Q.y * R.z + P.y * Q.z * R.x + P.z * Q.x * R.y - P.z * Q.y * R.x -
                   P.y * Q.x * R.z - P.x * Q.z * R.y;

  assert_msg(Det != 0, name << "::lin_solve error : matrix is not invertible" << P << "; " << Q
                            << "; "
                            << "; " << R);
  const ring invDet = 1 / Det;
  assert_msg(invDet != 0, name << "::lin_solve error : determinant is not invertible" << Det << "; "
                               << P << "; " << Q << "; "
                               << "; " << R);
  return R3((Q.y * R.z - Q.z * R.y) * x + (Q.z * R.x - Q.x * R.z) * y + (Q.x * R.y - Q.y * R.x) * z,
            (R.y * P.z - R.z * P.y) * x + (R.z * P.x - R.x * P.z) * y + (R.x * P.y - R.y * P.x) * z,
            (P.y * Q.z - P.z * Q.y) * x + (P.z * Q.x - P.x * Q.z) * y +
              (P.x * Q.y - P.y * Q.x) * z) *
         invDet;
}

// matrices symétriques
class sym3 : vector_space< sym3, double >, equality_comparable< sym3 > {
 public:
  double xx, yy, zz, xy, yz, zx;
  sym3( ) : xx(1), yy(1), zz(1), xy(0), yz(0), zx(0) {}

  explicit sym3(double lambda) : xx(lambda), yy(lambda), zz(lambda), xy(0), yz(0), zx(0) {}

  sym3(double XX, double YY, double ZZ, double XY, double YZ, double ZX)
    : xx(XX), yy(YY), zz(ZZ), xy(XY), yz(YZ), zx(ZX) {}

  sym3(double lambda, double mu, const R3 &P) {
    const R3 Q = P / P.norm( );
    const double nu = lambda - mu;

    xx = nu * Q.x * Q.x + mu;
    yy = nu * Q.y * Q.y + mu;
    zz = nu * Q.z * Q.z + mu;
    xy = nu * Q.x * Q.y;
    yz = nu * Q.y * Q.z;
    zx = nu * Q.z * Q.x;
  }

  template< class ring >
  double norm2(const TriDim< ring > &P) const {
    return xx * (P.x * P.x) + yy * (P.y * P.y) + zz * (P.z * P.z) + 2 * xy * (P.x * P.y) +
           2 * yz * (P.y * P.z) + 2 * zx * (P.z * P.x);
  }

  template< class ring >
  double norm(const TriDim< ring > &P) const {
    return sqrt(norm2(P));
  }

  template< class ring >
  double scal(const TriDim< ring > &P, const TriDim< ring > &Q) const {
    return xx * (P.x * Q.x) + yy * (P.y * Q.y) + zz * (P.z * Q.z) + xy * (P.x * Q.y + P.y * Q.x) +
           yz * (P.y * Q.z + P.z * Q.y) + zx * (P.z * Q.x + P.x * Q.z);
  }

  double det( ) const {
    return xx * yy * zz + 2 * xy * yz * zx - xx * yz * yz - yy * zx * zx - zz * xy * xy;
  }

  double trace( ) const { return xx + yy + zz; }

  double invariant( ) const { return xx * xx + yy * yy + zz * zz - xy * xy - yz * yz - zx * zx; }

  bool isDiagonal( ) const { return xy == 0 && yz == 0 && zx == 0; }

  bool isPositiveDefinite( ) const { return det( ) > 0 && trace( ) > 0 && invariant( ) > 0; }

  void reduced_basis(Z3 Basis[3]) const;
  bool operator==(const sym3 &S) const {
    return xx == S.xx && yy == S.yy && zz == S.zz && xy == S.xy && yz == S.yz && zx == S.zx;
  }

  sym3 &operator+=(double mu) {
    xx += mu;
    yy += mu;
    zz += mu;
    return *this;
  }

  sym3 &operator*=(double lambda) {
    xx *= lambda;
    yy *= lambda;
    zz *= lambda;
    xy *= lambda;
    yz *= lambda;
    zx *= lambda;
    return *this;
  }

  sym3 &operator/=(double lambda) {
    assert_msg(lambda != 0, "sym3 error : division by 0");
    return operator*=(1 / lambda);
  }

  sym3 &operator+=(const sym3 &m) {
    xx += m.xx;
    yy += m.yy;
    zz += m.zz;
    xy += m.xy;
    yz += m.yz;
    zx += m.zx;
    return *this;
  }

  sym3 &operator-=(const sym3 &m) {
    xx -= m.xx;
    yy -= m.yy;
    zz -= m.zz;
    xy -= m.xy;
    yz -= m.yz;
    zx -= m.zx;
    return *this;
  }

  const R3 operator*(const R3 &P) const {
    return R3(xx * P.x + xy * P.y + zx * P.z, xy * P.x + yy * P.y + yz * P.z,
              zx * P.x + yz * P.y + zz * P.z);
  }

  const sym3 comatrix( ) const;
  const R3 KernelRep( ) const;    // matrix needs to be rank 2, otherwise a null vector is returned.
  const sym3 inverse( ) const;
  template< class ring >
  sym3 conjugate(const TriDim< ring > &u, const TriDim< ring > &v, const TriDim< ring > &w) const {
    return sym3(norm2(u), norm2(v), norm2(w), scal(u, v), scal(v, w), scal(w, u));
  }
};

inline const sym3 sym3::comatrix( ) const {
  return sym3(yy * zz - square(yz), xx * zz - square(zx), xx * yy - square(xy), yz * zx - xy * zz,
              zx * xy - xx * yz, xy * yz - yy * zx);
}

inline const sym3 sym3::inverse( ) const {
  const double Det = det( );

  assert_msg(Det != 0, "sym3::inverse error : Matrix is not invertible " << *this);

  return comatrix( ) / Det;
}

inline void sym3::reduced_basis(Z3 Basis[3]) const {
  assert_msg(isPositiveDefinite( ), "sym3::reduced_basis error : matrix is not positive definite.");

  Z3 &b0 = Basis[0], &b1 = Basis[1], &b2 = Basis[2];
  b0 = Z3(1, 0, 0);
  b1 = Z3(0, 1, 0);
  b2 = Z3(0, 0, 1);
  double n0 = norm2(b0), n1 = norm2(b1), n2 = norm2(b2);    // Squared norms

  if (n2 < n1) {
    swap(b1, b2);
    swap(n1, n2);
  }    // Sorting the canonical basis by norm

  if (n1 < n0) {
    swap(b0, b1);
    swap(n0, n1);
  }

  if (n2 < n1) {
    swap(b1, b2);
    swap(n1, n2);
  }

  if (isDiagonal( )) {
    return;
  }

  while (true) {
    while (true) {
      b1 -= round(scal(b0, b1) / n0) * b0;
      n1 = norm2(b1);
      if (n0 <= n1) {
        break;
      }

      swap(b0, b1);
      swap(n0, n1);
    }

    const sym2 Gram(n0, scal(b0, b1), n1);
    R2 P = Gram.solve(R2(scal(b0, b2), scal(b1, b2)));
    const Z2 Q(round(P.x), round(P.y));
    P -= ZdtoRd(Q);
    b2 -= Q.lin_sum(b0, b1);
    n2 = norm2(b2);
    const Z3 b20 = b2 - sign(P.x) * b0;
    const double n20 = norm2(b20);
    const Z3 b21 = b2 - sign(P.y) * b1;
    const double n21 = norm2(b21);
    const double n2p = min(n2, min(n20, n21));
    if (n20 == n2p) {
      b2 = b20;
    }

    if (n21 == n2p) {
      b2 = b21;
    }

    n2 = n2p;
    if (n1 <= n2) {
      break;
    }

    swap(b1, b2);
    swap(n1, n2);
    if (n0 <= n1) {
      continue;
    }

    swap(b0, b1);
    swap(n0, n1);
  }
}

inline const R3 sym3::KernelRep( ) const {
  const sym3 co = comatrix( );
  const double norms2[3] = {square(co.xx) + square(co.xy) + square(co.zx),
                            square(co.xy) + square(co.yy) + square(co.yz),
                            square(co.zx) + square(co.yz) + square(co.zz)};

  return norms2[1] > norms2[2]
           ? (norms2[0] > norms2[1] ? R3(co.xx, co.xy, co.zx) : R3(co.xy, co.yy, co.yz))
           : (norms2[0] > norms2[2] ? R3(co.xx, co.xy, co.zx) : R3(co.zx, co.yz, co.zz));
}

/************************** affichage de variables et tableaux **********************/

template< class ring >
inline ostream &operator<<(ostream &f, const BiDim< ring > &P) {
  f << P.x << " " << P.y;
  return f;
}

template< class ring >
inline ostream_math operator<<(ostream_math f, const BiDim< ring > &P) {
  if (f.format == Mathematica) {
    f << "{" << P.x << "," << P.y << "}";
  } else {
    f.os << P;
  }

  return f;
}

template< class ring >
inline ostream &operator<<(ostream &f, const TriDim< ring > &P) {
  f << P.x << " " << P.y << " " << P.z;
  return f;
}

template< class ring >
inline ostream_math operator<<(ostream_math f, const TriDim< ring > &P) {
  if (f.format == Mathematica) {
    f << "{" << P.x << "," << P.y << "," << P.z << "}";
  } else {
    f.os << P;
  }

  return f;
}

inline ostream &operator<<(ostream &f, const sym2 &S) {
  f << "xx : " << S.xx << "; xy : " << S.xy << "; yy : " << S.yy;
  return f;
}

inline ostream_math operator<<(ostream_math f, const sym2 &S) {
  if (f.format == Mathematica) {
    f << "{{" << S.xx << "," << S.xy << "},{" << S.xy << "," << S.yy << "}}";
  } else {
    f.os << S;
  }

  return f;
}

inline ostream &operator<<(ostream &f, const sym3 &S) {
  f << "xx : " << S.xx << "; yy : " << S.yy << "; zz : " << S.zz;
  f << "; xy : " << S.xy << "; yz : " << S.yz << "; zx : " << S.zx;
  return f;
}

inline ostream_math operator<<(ostream_math f, const sym3 &S) {
  if (f.format == Mathematica) {
    f << "{{" << S.xx << "," << S.xy << "," << S.zx << "},{" << S.xy << "," << S.yy << "," << S.yz
      << "},{" << S.zx << "," << S.yz << "," << S.zz << "}}";
  } else {
    f.os << S;
  }

  return f;
}

// affichage de tableaux
template< class ForwardIterator, class Zd >
void print_array(ostream_math f, ForwardIterator first, ForwardIterator last, Zd N,
                 bool one_per_line = false) {
  int prod[N.dim];

  prod[0] = N[0];

  for (int k = 1; k < N.dim; ++k) {
    prod[k] = N[k] * prod[k - 1];
  }

  for (int k = 0; k < N.dim; ++k) {
    f << "{";
  }

  if (first != last) {
    f << *first++;
  }

  int i = 1;

  for (; first != last; ++first, ++i) {
    for (int k = 0; k < N.dim; ++k) {
      if (i % (prod[k]) == 0) {
        f << "}";
      }
    }

    f << ",";
    if (one_per_line) {
      f << "\n";
    }

    for (int k = 0; k < N.dim; ++k) {
      if (i % (prod[k]) == 0) {
        f << "{";
      }
    }

    f << *first;
  }

  for (int k = 0; k < N.dim; ++k) {
    f << "}";
  }

  assert(i == prod[N.dim - 1]);    // size and format must match
}

// *********************** Metric classes **********************

// 2d
class Metric2 {
 public:
  double lip;
  Metric2( ) : lip(0){};    // ?? value -1 is reserved to identify the euclidean metric ?? Useful ??
  virtual const sym2 operator( )(const R2 &P) const { return sym2(1, 0, 1); }

  virtual ~Metric2( ) {}
};

class FctMetric2 : public Metric2 {
 public:
  const sym2 (*metric_)(const R2 &);
  FctMetric2(const sym2 (*metric)(const R2 &), double Lip = 5) : metric_(metric) { lip = Lip; }

  const sym2 operator( )(const R2 &P) const { return metric_(P); }
};

// 3d
class Metric3 {
 public:
  double lip;
  Metric3( ) : lip(0){};
  virtual const sym3 operator( )(const R3 &P) const { return sym3(1, 1, 1, 0, 0, 0); }

  virtual ~Metric3( ) {}
};

class FctMetric3 : public Metric2 {
 public:
  const sym3 (*metric_)(const R3 &);
  FctMetric3(const sym3 (*metric)(const R3 &), double Lip = 5) : metric_(metric) { lip = Lip; }

  const sym3 operator( )(const R3 &P) const { return metric_(P); }
};

#endif
