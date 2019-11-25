/****************************************************************************/
/* This file is part of FreeFEM.                                            */
/*                                                                          */
/* FreeFEM is free software: you can redistribute it and/or modify          */
/* it under the terms of the GNU Lesser General Public License as           */
/* published by the Free Software Foundation, either version 3 of           */
/* the License, or (at your option) any later version.                      */
/*                                                                          */
/* FreeFEM is distributed in the hope that it will be useful,               */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of           */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            */
/* GNU Lesser General Public License for more details.                      */
/*                                                                          */
/* You should have received a copy of the GNU Lesser General Public License */
/* along with FreeFEM. If not, see <http://www.gnu.org/licenses/>.          */
/****************************************************************************/
// SUMMARY : ...
// LICENSE : LGPLv3
// ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
// AUTHORS : Jean-Marie Mirebeau
// E-MAIL  : jean-marie.mirebeau@math.u-psud.fr

// *INDENT-OFF* //
// ff-c++-LIBRARY-dep:
// ff-c++-cpp-dep:
// *INDENT-ON* //

// This file contains two main routines :
// - one for the construction of a metric
// - one for the computing the derivatives of a polynomial sampled at lagrange points on a triangle

#ifndef _TensorK_h
#define _TensorK_h

#include <vector>
#include <cmath>
#include <algorithm>
#include "assert.h"
using std::vector;

// ***************** Quelques fonctions ************

// quelques fonctions

inline double square(double u) { return u * u; };
inline double max(double a, double b) { return a > b ? a : b; }

inline double max(double a, double b, double c) { return a > c ? max(a, b) : max(b, c); }

inline double max(double a, double b, double c, double d) {
  return a > d ? max(a, b, c) : max(b, c, d);
}

inline double max(double a, double b, double c, double d, double e) {
  return a > e ? max(a, b, c, d) : max(b, c, d, e);
}

// ***************** Prototypes ********************

class TensorK {
  vector< double > factorials;
  std::vector< double > exponents;
  double factorial(int n) const {
    assert(0 <= n <= factorials.size( ));
    return factorials[n];
  }

  double binomial(int n, int k) const {
    assert(0 <= n && n <= factorials.size( ) && 0 <= k && k <= n);
    return factorials[n] / (factorials[k] * factorials[n - k]);
  }

  const int t_deg;    // degree of the polynomials involved

 public:
  const int m_deg;
  const int r_deg;
  enum triangulation_type {
    Graded = 0,
    Quasi_Acute = 1,
    Quasi_Acute_Unrefined = 2,
    Quasi_Acute_Proved = 3
  };
  const triangulation_type ttype;
  const double p_exp;
  enum which_matrix { M0_alone = 0, M1_alone = 1, M0_M1_weighted_sum = 2 };
  const which_matrix wmat;
  static const int d_dim = 2;
  const double gamma_exp;
  const double homog_exp;
  const bool is_valid;

  TensorK(int m_deg_, int r_deg_, triangulation_type ttype_, which_matrix wmat, double p_exp_);

  void getM(const double *D, double M[3])
    const;    // L^Infinity metric        // size required : D[m_deg+1], m_deg+1 == m_dim
  void equilibrate(const double M[3],
                   double Me[3]) const;    // Me = (det M)^(-1/((m-r)p+d)) M = (det M)^gamma_exp M.

 private:
  void getMc(const double *D, double Mc[3]) const;    // size required : D[t_deg+1]
  void getM0(const double E[2], double c, double s,
             double M0[3]) const;    // input : eigenvalues and eigenvector of Mc
  void getM1(const double *D, double c, double s,
             double M1[3]) const;                     // size required : D[t_deg+1]
  void getMs(const double *D, double Ms[3]) const;    // size required : D[t_deg+1]

  void rotate(const double *D, double *Dr, double c,
              double s) const;                          // size required : D[t_deg+1]
  void getSquare(const double *D, double *Ds) const;    // get |d^r pi|^2. D[m_deg+1], Ds[t_deg+1]

 public:    // derivative estimation
  void getDerivatives(const std::vector< double > &DOFt, const R2 invHauteur[3], double *f) const;

 private:
  template< int m >
  void Derivatives(const std::vector< double > &DOFt, const R2 invHauteur[3], double f[m]) const;

 public:    // symmetric matrix utilities
  static void EigenSym(const double S[3], double E[2]);
  static void EigenSysSym(const double S[3], double Eigen[2], double &c, double &s);
  static void MakeEigenSym(double S[3], double vap[2], double c, double s);
  static void AffSym(double S[3], double a, double b);
  static void MaxSym(double S[3], double lambda);
  static double det(const double S[3]) { return S[0] * S[2] - S[1] * S[1]; }

  static void PowSym(double S[3], double p);

 private:    // debugging
  friend int main(int argc, const char *argv[]);
};

// **************** Constructor ********************

TensorK::TensorK(int m_deg_, int r_deg_, triangulation_type ttype_ = Graded,
                 which_matrix wmat_ = M1_alone, double p_exp_ = 2)
  : m_deg(m_deg_), r_deg(r_deg_), ttype(ttype_), wmat(wmat_), p_exp(p_exp_),
    t_deg((ttype_ == Quasi_Acute_Proved) ? 2 * (m_deg_ - r_deg_) : m_deg_),
    gamma_exp(-1. / ((m_deg_ - r_deg_) * p_exp_ + d_dim)),
    homog_exp(1. / ((m_deg_ - r_deg_) * (ttype_ == Quasi_Acute_Proved ? 2. : 1.))),
    is_valid(m_deg_ >= 2 && m_deg_ <= 5 && r_deg_ >= 0 && r_deg_ <= m_deg_ - 1 && ttype_ >= 0 &&
             ttype_ <= 3 && wmat_ >= 0 && wmat_ <= 2 && p_exp_ >= 0) {
  factorials.resize(t_deg + 1);
  factorials[0] = 1;

  for (int i = 1; i <= t_deg; ++i) {
    factorials[i] = i * factorials[i - 1];
  }

  exponents.resize(t_deg + 1);

  for (int k = 1; k <= t_deg; ++k) {
    switch (ttype) {
      case Graded:
        exponents[k] = 1. / k;
        break;

      case Quasi_Acute:
        exponents[k] = 1. / std::min(k, m_deg - r_deg);
        break;

      case Quasi_Acute_Unrefined:
        exponents[k] = (k <= m_deg - r_deg) ? 1. / k : 1. / (k - 1. / p_exp);
        break;

      case Quasi_Acute_Proved:
        exponents[k] = 1. / k;
        break;
    }
  }
}

// void Derivatives(const double DOFt[3],  const R2 invHauteur[3], double &fx,      double &fy);
// void Derivatives(const double DOFt[6],  const R2 invHauteur[3], double &fxx,     double &fxy,
// double &fyy); void Derivatives(const double DOFt[10], const R2 invHauteur[3], double &fxxx,
// double &fxxy,   double &fxyy,   double &fyyy); void Derivatives(const double DOFt[15], const R2
// invHauteur[3], double &fxxxx,   double &fxxxy,  double &fxxyy,  double &fxyyy, double &fyyyy);

// ********************* Matrices **********************

void TensorK::getMc(const double *D, double Mc[3]) const {
  Mc[0] = 0;
  Mc[1] = 0;
  Mc[2] = 0;

  for (int i = 0; i < t_deg; ++i) {
    Mc[0] += binomial(t_deg - 1, i) * D[i] * D[i];
    Mc[1] += binomial(t_deg - 1, i) * D[i] * D[i + 1];
    Mc[2] += binomial(t_deg - 1, i) * D[i + 1] * D[i + 1];
  }
}

void TensorK::getM0(const double E[2], double c, double s,
                    double M0[3]) const {    // size required : D[t_deg+1]
  double EPow[2] = {pow(2 * E[0], homog_exp), pow(2 * E[1], homog_exp)};

  MakeEigenSym(M0, EPow, c, s);
}

void TensorK::getM1(const double *D, double c, double s, double M1[3]) const {
  double DChg[t_deg + 1];

  rotate(D, DChg, c, -s);
  double DChgMax = 0;

  for (int r = 0; r <= t_deg; ++r) {
    DChg[r] = fabs(DChg[r]);
    DChgMax = max(DChgMax, DChg[r]);
  }

  if (DChgMax == 0) {
    M1[0] = 0;
    M1[1] = 0;
    M1[2] = 0;
    return;
  }

  double vap[2] = {0, 0};    // The eigenvalues of the matrix returned.

  for (int r = 0; r < t_deg; ++r) {
    vap[0] = max(vap[0], pow(DChg[r] / DChgMax, exponents[t_deg - r]));
  }

  for (int r = 1; r <= t_deg; ++r) {
    vap[1] = max(vap[1], pow(DChg[r] / DChgMax, exponents[r]));
  }

  const double scal = pow(DChgMax, homog_exp);
  const double scal2 = pow(2., m_deg * homog_exp);

  for (int i = 0; i < 2; ++i) {
    vap[i] *= scal;
    vap[i] = vap[i] * vap[i];
    vap[i] *= scal2;
  }

  MakeEigenSym(M1, vap, c, s);
}

void TensorK::getMs(const double *D, double Ms[3]) const {
  double Mc[3];
  double E[2];
  double c, s;

  getMc(D, Mc);
  EigenSysSym(Mc, E, c, s);

  switch (wmat) {
    case M0_M1_weighted_sum: {
      double M1[3];
      getM1(D, c, s, M1);
      double M0[3];
      getM0(E, c, s, M0);
      const double r = max(2 - E[1] / E[0], 0);

      for (int i = 0; i < 3; ++i) {
        Ms[i] = r * M0[i] + (1 - r) * M1[i];
      }

      return;
    }

    case M1_alone:
      getM1(D, c, s, Ms);
      return;

    case M0_alone:
      getM0(E, c, s, Ms);
      return;

    default:
      std::cout << "TensorK::getMs error ! Unsupplied case." << std::endl;
      break;
  }
}

void TensorK::getM(const double *D, double *M) const {
  if (ttype == Quasi_Acute_Proved) {
    double Ds[t_deg + 1];
    getSquare(D, Ds);
    getMs(Ds, M);
  } else {
    getMs(D, M);
  }
}

void TensorK::getSquare(const double *D, double *Ds) const {
  assert(ttype == Quasi_Acute_Proved);

  const int diff_deg = m_deg - r_deg;
  assert(2 * diff_deg == t_deg);

  for (int i = 0; i <= t_deg; ++i) {
    Ds[i] = 0;
  }

  for (int r = 0; r <= r_deg; ++r) {
    const double *const Dr = D + r;

    for (int p = 0; p <= diff_deg; ++p) {
      for (int q = 0; q <= diff_deg; ++q) {
        Ds[p + q] +=
          binomial(diff_deg, p) * binomial(diff_deg, q) / binomial(t_deg, p + q) * Dr[p] * Dr[q];
      }
    }
  }
}

void TensorK::rotate(const double *D, double *Dr, double c, double s) const {
  for (int i = 0; i <= t_deg; ++i) {
    Dr[i] = 0;
  }

  double cpow[t_deg + 1];
  cpow[0] = 1;

  for (int i = 1; i <= t_deg; ++i) {
    cpow[i] = c * cpow[i - 1];
  }

  double spow[t_deg + 1];
  spow[0] = 1;

  for (int i = 1; i <= t_deg; ++i) {
    spow[i] = s * spow[i - 1];
  }

  double parity[t_deg + 1];
  parity[0] = 1;

  for (int i = 1; i <= t_deg; ++i) {
    parity[i] = -parity[i - 1];
  }

  for (int i = 0; i <= t_deg; ++i) {
    const int j = t_deg - i;

    for (int p = 0; p <= i; ++p) {
      const int q = i - p;

      for (int u = 0; u <= j; ++u) {
        const int v = j - u;
        Dr[p + u] += D[i] * binomial(t_deg, i) * binomial(i, p) * binomial(j, u) /
                     binomial(t_deg, p + u) * cpow[p + v] * spow[q + u] * parity[q];
      }
    }
  }
}

void TensorK::equilibrate(const double M[3], double Me[3]) const {
  const double scal = pow(det(M), gamma_exp);

  for (int i = 0; i < 3; ++i) {
    Me[i] = scal * M[i];
  }
};    // Me = (det M)^(-1/((m-r)p+d)) M

// ******************** Matrix Utilities ********************
void TensorK::EigenSym(const double S[3], double E[2]) {
  double hDiff = sqrt(square(S[0] - S[2]) / 4 + square(S[1]));
  double hSum = (S[0] + S[2]) / 2;

  E[0] = hSum - hDiff;
  E[1] = hSum + hDiff;
  return;
}

void TensorK::EigenSysSym(const double S[3], double Eigen[2], double &c, double &s) {
  EigenSym(S, Eigen);    // calcul des valeurs propres de la matrice Mc, dans l'ordre croissant
  const double div = Eigen[0] * Eigen[0] - Eigen[1] * Eigen[1];
  if (div == 0.) {
    c = 1;
    s = 0;
  }    // si les vap sont égales (la matrice est positive)
  else {
    c = (Eigen[0] * S[0] - Eigen[1] * S[2]) / div;
    c = c >= 0 ? sqrt(c) : 0;    // normalement toujours positif
    s = (Eigen[0] * S[2] - Eigen[1] * S[0]) / div;
    s = s >= 0 ? sqrt(s) : 0;    // idem
    s = (Eigen[0] - Eigen[1]) * S[1] > 0 ? s : -s;
  }
}

void TensorK::MakeEigenSym(double S[3], double vap[2], double c, double s) {
  S[0] = vap[0] * c * c + vap[1] * s * s;
  S[1] = (vap[0] - vap[1]) * c * s;
  S[2] = vap[0] * s * s + vap[1] * c * c;
}

void TensorK::AffSym(double S[3], double a, double b) {
  S[0] = a * S[0] + b;
  S[1] = a * S[1];
  S[2] = a * S[2] + b;
}

void TensorK::MaxSym(double S[3], double lambda) {
  double E[2];

  EigenSym(S, E);
  if (lambda <= E[0]) {
    return;
  }

  if (E[1] <= lambda) {
    S[0] = lambda;
    S[1] = 0;
    S[2] = lambda;
    return;
  }

  ;
  AffSym(S, (E[1] - lambda) / (E[1] - E[0]), E[1] * (lambda - E[0]) / (E[1] - E[0]));
}

void TensorK::PowSym(double S[3], double p) {
  double E[2];    // old eigenvalues

  EigenSym(S, E);

  double Ep[2];    // new eigenvalues
  if (p == -2) {
    Ep[0] = 1 / (E[0] * E[0]);
    Ep[1] = 1 / (E[1] * E[1]);
  } else if (p == -0.5) {
    Ep[0] = 1 / sqrt(E[0]);
    Ep[1] = 1 / sqrt(E[1]);
  } else {
    Ep[0] = pow(E[0], p);
    Ep[1] = pow(E[1], p);
  }

  const double diff = E[1] - E[0];
  if (diff == 0) {
    S[0] = Ep[0];
    S[1] = 0;
    S[2] = Ep[0];
    return;
  }    // scalar case

  AffSym(S, (Ep[1] - Ep[0]) / diff, (Ep[0] * E[1] - Ep[1] * E[0]) / diff);
}

// ********************* Derivatives estimation ********************

template<>
void TensorK::Derivatives< 2 >(const std::vector< double > &DOFt, const R2 invHauteur[3],
                               double f[2]) const {
  f[0] = -DOFt[0] * invHauteur[0].x - DOFt[1] * invHauteur[1].x - DOFt[2] * invHauteur[2].x;

  f[1] = -DOFt[0] * invHauteur[0].y - DOFt[1] * invHauteur[1].y - DOFt[2] * invHauteur[2].y;
}

template<>
void TensorK::Derivatives< 3 >(const std::vector< double > &DOFt, const R2 invHauteur[3],
                               double f[3]) const {
  f[0] = 4 * DOFt[0] * invHauteur[0].x * invHauteur[0].x +
         4 * DOFt[1] * invHauteur[1].x * invHauteur[1].x +
         4 * DOFt[2] * invHauteur[2].x * invHauteur[2].x +
         8 * DOFt[3] * invHauteur[1].x * invHauteur[2].x +
         8 * DOFt[4] * invHauteur[2].x * invHauteur[0].x +
         8 * DOFt[5] * invHauteur[0].x *
           invHauteur[1].x;    // dérivée seconde en x de la fonction P2 sur le triangle d'intérêt.

  f[1] = 4 * DOFt[0] * invHauteur[0].x * invHauteur[0].y +
         4 * DOFt[1] * invHauteur[1].x * invHauteur[1].y +
         4 * DOFt[2] * invHauteur[2].x * invHauteur[2].y +
         4 * DOFt[3] * (invHauteur[1].x * invHauteur[2].y + invHauteur[1].y * invHauteur[2].x) +
         4 * DOFt[4] * (invHauteur[2].x * invHauteur[0].y + invHauteur[2].y * invHauteur[0].x) +
         4 * DOFt[5] * (invHauteur[0].x * invHauteur[1].y + invHauteur[0].y * invHauteur[1].x);

  f[2] = 4 * DOFt[0] * invHauteur[0].y * invHauteur[0].y +
         4 * DOFt[1] * invHauteur[1].y * invHauteur[1].y +
         4 * DOFt[2] * invHauteur[2].y * invHauteur[2].y +
         8 * DOFt[3] * invHauteur[1].y * invHauteur[2].y +
         8 * DOFt[4] * invHauteur[2].y * invHauteur[0].y +
         8 * DOFt[5] * invHauteur[0].y * invHauteur[1].y;
}

template<>
void TensorK::Derivatives< 4 >(const std::vector< double > &DOFt, const R2 invHauteur[3],
                               double f[4]) const {
  f[0] = -6. * ((9. / 2.) * DOFt[0] * invHauteur[0].x * invHauteur[0].x * invHauteur[0].x +
                (9. / 2.) * DOFt[1] * invHauteur[1].x * invHauteur[1].x * invHauteur[1].x +
                (9. / 2.) * DOFt[2] * invHauteur[2].x * invHauteur[2].x * invHauteur[2].x +
                (27. / 2.) * DOFt[3] * invHauteur[1].x * invHauteur[1].x * invHauteur[2].x +
                (27. / 2.) * DOFt[4] * invHauteur[1].x * invHauteur[2].x * invHauteur[2].x +
                (27. / 2.) * DOFt[5] * invHauteur[2].x * invHauteur[2].x * invHauteur[0].x +
                (27. / 2.) * DOFt[6] * invHauteur[2].x * invHauteur[0].x * invHauteur[0].x +
                (27. / 2.) * DOFt[7] * invHauteur[0].x * invHauteur[0].x * invHauteur[1].x +
                (27. / 2.) * DOFt[8] * invHauteur[0].x * invHauteur[1].x * invHauteur[1].x +
                27. * DOFt[9] * invHauteur[0].x * invHauteur[1].x * invHauteur[2].x);

  f[1] = -6. * ((9. / 2.) * DOFt[0] * invHauteur[0].x * invHauteur[0].x * invHauteur[0].y +
                (9. / 2.) * DOFt[1] * invHauteur[1].x * invHauteur[1].x * invHauteur[1].y +
                (9. / 2.) * DOFt[2] * invHauteur[2].x * invHauteur[2].x * invHauteur[2].y +
                (27. / 2.) * DOFt[3] *
                  (invHauteur[1].x * invHauteur[1].x * invHauteur[2].y * (1. / 3.) +
                   invHauteur[1].y * invHauteur[1].x * invHauteur[2].x * (2. / 3.)) +
                (27. / 2.) * DOFt[4] *
                  (invHauteur[1].x * invHauteur[2].x * invHauteur[2].y * (2. / 3.) +
                   invHauteur[1].y * invHauteur[2].x * invHauteur[2].x * (1. / 3.)) +
                (27. / 2.) * DOFt[5] *
                  (invHauteur[2].x * invHauteur[2].x * invHauteur[0].y * (1. / 3.) +
                   invHauteur[2].y * invHauteur[2].x * invHauteur[0].x * (2. / 3.)) +
                (27. / 2.) * DOFt[6] *
                  (invHauteur[2].x * invHauteur[0].x * invHauteur[0].y * (2. / 3.) +
                   invHauteur[2].y * invHauteur[0].x * invHauteur[0].x * (1. / 3.)) +
                (27. / 2.) * DOFt[7] *
                  (invHauteur[0].x * invHauteur[0].x * invHauteur[1].y * (1. / 3.) +
                   invHauteur[0].y * invHauteur[0].x * invHauteur[1].x * (2. / 3.)) +
                (27. / 2.) * DOFt[8] *
                  (invHauteur[0].x * invHauteur[1].x * invHauteur[1].y * (2. / 3.) +
                   invHauteur[0].y * invHauteur[1].x * invHauteur[1].x * (1. / 3.)) +
                27. * DOFt[9] *
                  (invHauteur[0].x * invHauteur[1].x * invHauteur[2].y / 3. +
                   invHauteur[0].x * invHauteur[1].y * invHauteur[2].x / 3. +
                   invHauteur[0].y * invHauteur[1].x * invHauteur[2].x / 3.));
  f[2] = -6. * ((9. / 2.) * DOFt[0] * invHauteur[0].x * invHauteur[0].y * invHauteur[0].y +
                (9. / 2.) * DOFt[1] * invHauteur[1].x * invHauteur[1].y * invHauteur[1].y +
                (9. / 2.) * DOFt[2] * invHauteur[2].x * invHauteur[2].y * invHauteur[2].y +
                (27. / 2.) * DOFt[3] *
                  (invHauteur[1].y * invHauteur[1].y * invHauteur[2].x * (1. / 3.) +
                   invHauteur[1].x * invHauteur[1].y * invHauteur[2].y * (2. / 3.)) +
                (27. / 2.) * DOFt[4] *
                  (invHauteur[1].y * invHauteur[2].y * invHauteur[2].x * (2. / 3.) +
                   invHauteur[1].x * invHauteur[2].y * invHauteur[2].y * (1. / 3.)) +
                (27. / 2.) * DOFt[5] *
                  (invHauteur[2].y * invHauteur[2].y * invHauteur[0].x * (1. / 3.) +
                   invHauteur[2].x * invHauteur[2].y * invHauteur[0].y * (2. / 3.)) +
                (27. / 2.) * DOFt[6] *
                  (invHauteur[2].y * invHauteur[0].y * invHauteur[0].x * (2. / 3.) +
                   invHauteur[2].x * invHauteur[0].y * invHauteur[0].y * (1. / 3.)) +
                (27. / 2.) * DOFt[7] *
                  (invHauteur[0].y * invHauteur[0].y * invHauteur[1].x * (1. / 3.) +
                   invHauteur[0].x * invHauteur[0].y * invHauteur[1].y * (2. / 3.)) +
                (27. / 2.) * DOFt[8] *
                  (invHauteur[0].y * invHauteur[1].y * invHauteur[1].x * (2. / 3.) +
                   invHauteur[0].x * invHauteur[1].y * invHauteur[1].y * (1. / 3.)) +
                27. * DOFt[9] *
                  (invHauteur[0].y * invHauteur[1].y * invHauteur[2].x / 3. +
                   invHauteur[0].y * invHauteur[1].x * invHauteur[2].y / 3. +
                   invHauteur[0].x * invHauteur[1].y * invHauteur[2].y / 3.));

  f[3] = -6. * ((9. / 2.) * DOFt[0] * invHauteur[0].y * invHauteur[0].y * invHauteur[0].y +
                (9. / 2.) * DOFt[1] * invHauteur[1].y * invHauteur[1].y * invHauteur[1].y +
                (9. / 2.) * DOFt[2] * invHauteur[2].y * invHauteur[2].y * invHauteur[2].y +
                (27. / 2.) * DOFt[3] * invHauteur[1].y * invHauteur[1].y * invHauteur[2].y +
                (27. / 2.) * DOFt[4] * invHauteur[1].y * invHauteur[2].y * invHauteur[2].y +
                (27. / 2.) * DOFt[5] * invHauteur[2].y * invHauteur[2].y * invHauteur[0].y +
                (27. / 2.) * DOFt[6] * invHauteur[2].y * invHauteur[0].y * invHauteur[0].y +
                (27. / 2.) * DOFt[7] * invHauteur[0].y * invHauteur[0].y * invHauteur[1].y +
                (27. / 2.) * DOFt[8] * invHauteur[0].y * invHauteur[1].y * invHauteur[1].y +
                27. * DOFt[9] * invHauteur[0].y * invHauteur[1].y * invHauteur[2].y);
}

template<>
void TensorK::Derivatives< 5 >(const std::vector< double > &DOFt, const R2 invHauteur[3],
                               double f[5]) const {
  f[0] =
    (32. / 3.) * DOFt[0] * 24 * invHauteur[0].x * invHauteur[0].x * invHauteur[0].x *
      invHauteur[0].x +
    (32. / 3.) * DOFt[1] * 24 * invHauteur[1].x * invHauteur[1].x * invHauteur[1].x *
      invHauteur[1].x +
    (32. / 3.) * DOFt[2] * 24 * invHauteur[2].x * invHauteur[2].x * invHauteur[2].x *
      invHauteur[2].x +
    (128. / 3.) * DOFt[3] * 24 * invHauteur[1].x * invHauteur[1].x * invHauteur[1].x *
      invHauteur[2].x +
    64 * DOFt[4] * 24 * invHauteur[1].x * invHauteur[1].x * invHauteur[2].x * invHauteur[2].x +
    (128. / 3.) * DOFt[5] * 24 * invHauteur[1].x * invHauteur[2].x * invHauteur[2].x *
      invHauteur[2].x +
    (128. / 3.) * DOFt[6] * 24 * invHauteur[2].x * invHauteur[2].x * invHauteur[2].x *
      invHauteur[0].x +
    64 * DOFt[7] * 24 * invHauteur[2].x * invHauteur[2].x * invHauteur[0].x * invHauteur[0].x +
    (128. / 3.) * DOFt[8] * 24 * invHauteur[2].x * invHauteur[0].x * invHauteur[0].x *
      invHauteur[0].x +
    (128. / 3.) * DOFt[9] * 24 * invHauteur[0].x * invHauteur[0].x * invHauteur[0].x *
      invHauteur[1].x +
    64 * DOFt[10] * 24 * invHauteur[0].x * invHauteur[0].x * invHauteur[1].x * invHauteur[1].x +
    (128. / 3.) * DOFt[11] * 24 * invHauteur[0].x * invHauteur[1].x * invHauteur[1].x *
      invHauteur[1].x +
    128. * DOFt[12] * 24 * invHauteur[0].x * invHauteur[1].x * invHauteur[2].x * invHauteur[2].x +
    128. * DOFt[13] * 24 * invHauteur[2].x * invHauteur[0].x * invHauteur[1].x * invHauteur[1].x +
    128. * DOFt[14] * 24 * invHauteur[1].x * invHauteur[2].x * invHauteur[0].x * invHauteur[0].x;

  f[1] = (32. / 3.) * DOFt[0] * 24 * invHauteur[0].x * invHauteur[0].x * invHauteur[0].x *
           invHauteur[0].y +
         (32. / 3.) * DOFt[1] * 24 * invHauteur[1].x * invHauteur[1].x * invHauteur[1].x *
           invHauteur[1].y +
         (32. / 3.) * DOFt[2] * 24 * invHauteur[2].x * invHauteur[2].x * invHauteur[2].x *
           invHauteur[2].y +

         (128. / 3.) * DOFt[3] *
           (6 * invHauteur[1].x * invHauteur[1].x * invHauteur[1].x * invHauteur[2].y +
            18 * invHauteur[1].x * invHauteur[1].x * invHauteur[1].y * invHauteur[2].x) +
         64 * DOFt[4] *
           (12 * invHauteur[1].x * invHauteur[1].x * invHauteur[2].x * invHauteur[2].y +
            12 * invHauteur[1].x * invHauteur[1].y * invHauteur[2].x * invHauteur[2].x) +
         (128. / 3.) * DOFt[5] *
           (18 * invHauteur[1].x * invHauteur[2].x * invHauteur[2].x * invHauteur[2].y +
            6 * invHauteur[1].y * invHauteur[2].x * invHauteur[2].x * invHauteur[2].x) +

         (128. / 3.) * DOFt[6] *
           (6 * invHauteur[2].x * invHauteur[2].x * invHauteur[2].x * invHauteur[0].y +
            18 * invHauteur[2].x * invHauteur[2].x * invHauteur[2].y * invHauteur[0].x) +
         64 * DOFt[7] *
           (12 * invHauteur[2].x * invHauteur[2].x * invHauteur[0].x * invHauteur[0].y +
            12 * invHauteur[2].x * invHauteur[2].y * invHauteur[0].x * invHauteur[0].x) +
         (128. / 3.) * DOFt[8] *
           (18 * invHauteur[2].x * invHauteur[0].x * invHauteur[0].x * invHauteur[0].y +
            6 * invHauteur[2].y * invHauteur[0].x * invHauteur[0].x * invHauteur[0].x) +

         (128. / 3.) * DOFt[9] *
           (6 * invHauteur[0].x * invHauteur[0].x * invHauteur[0].x * invHauteur[1].y +
            18 * invHauteur[0].x * invHauteur[0].x * invHauteur[0].y * invHauteur[1].x) +
         64 * DOFt[10] *
           (12 * invHauteur[0].x * invHauteur[0].x * invHauteur[1].x * invHauteur[1].y +
            12 * invHauteur[0].x * invHauteur[0].y * invHauteur[1].x * invHauteur[1].x) +
         (128. / 3.) * DOFt[11] *
           (18 * invHauteur[0].x * invHauteur[1].x * invHauteur[1].x * invHauteur[1].y +
            6 * invHauteur[0].y * invHauteur[1].x * invHauteur[1].x * invHauteur[1].x) +

         128. * DOFt[12] *
           (12 * invHauteur[0].x * invHauteur[1].x * invHauteur[2].x * invHauteur[2].y +
            6 * invHauteur[0].x * invHauteur[1].y * invHauteur[2].x * invHauteur[2].x +
            6 * invHauteur[0].y * invHauteur[1].x * invHauteur[2].x * invHauteur[2].x) +
         128. * DOFt[13] *
           (12 * invHauteur[2].x * invHauteur[0].x * invHauteur[1].x * invHauteur[1].y +
            6 * invHauteur[2].x * invHauteur[0].y * invHauteur[1].x * invHauteur[1].x +
            6 * invHauteur[2].y * invHauteur[0].x * invHauteur[1].x * invHauteur[1].x) +
         128. * DOFt[14] *
           (12 * invHauteur[1].x * invHauteur[2].x * invHauteur[0].x * invHauteur[0].y +
            6 * invHauteur[1].x * invHauteur[2].y * invHauteur[0].x * invHauteur[0].x +
            6 * invHauteur[1].y * invHauteur[2].x * invHauteur[0].x * invHauteur[0].x);

  f[2] = (32. / 3.) * DOFt[0] * 24 * invHauteur[0].x * invHauteur[0].x * invHauteur[0].y *
           invHauteur[0].y +
         (32. / 3.) * DOFt[1] * 24 * invHauteur[1].x * invHauteur[1].x * invHauteur[1].y *
           invHauteur[1].y +
         (32. / 3.) * DOFt[2] * 24 * invHauteur[2].x * invHauteur[2].x * invHauteur[2].y *
           invHauteur[2].y +

         (128. / 3.) * DOFt[3] *
           (12 * invHauteur[1].x * invHauteur[1].x * invHauteur[1].y * invHauteur[2].y +
            12 * invHauteur[1].x * invHauteur[1].y * invHauteur[1].y * invHauteur[2].x) +
         64 * DOFt[4] *
           (4 * invHauteur[1].x * invHauteur[1].x * invHauteur[2].y * invHauteur[2].y +
            16 * invHauteur[1].x * invHauteur[1].y * invHauteur[2].x * invHauteur[2].y +
            4 * invHauteur[1].y * invHauteur[1].y * invHauteur[2].x * invHauteur[2].x) +
         (128. / 3.) * DOFt[5] *
           (12 * invHauteur[1].x * invHauteur[2].x * invHauteur[2].y * invHauteur[2].y +
            12 * invHauteur[1].y * invHauteur[2].x * invHauteur[2].x * invHauteur[2].y) +

         (128. / 3.) * DOFt[6] *
           (12 * invHauteur[2].x * invHauteur[2].x * invHauteur[2].y * invHauteur[0].y +
            12 * invHauteur[2].x * invHauteur[2].y * invHauteur[2].y * invHauteur[0].x) +
         64 * DOFt[7] *
           (4 * invHauteur[2].x * invHauteur[2].x * invHauteur[0].y * invHauteur[0].y +
            16 * invHauteur[2].x * invHauteur[2].y * invHauteur[0].x * invHauteur[0].y +
            4 * invHauteur[2].y * invHauteur[2].y * invHauteur[0].x * invHauteur[0].x) +
         (128. / 3.) * DOFt[8] *
           (12 * invHauteur[2].x * invHauteur[0].x * invHauteur[0].y * invHauteur[0].y +
            12 * invHauteur[2].y * invHauteur[0].x * invHauteur[0].x * invHauteur[0].y) +

         (128. / 3.) * DOFt[9] *
           (12 * invHauteur[0].x * invHauteur[0].x * invHauteur[0].y * invHauteur[1].y +
            12 * invHauteur[0].x * invHauteur[0].y * invHauteur[0].y * invHauteur[1].x) +
         64 * DOFt[10] *
           (4 * invHauteur[0].x * invHauteur[0].x * invHauteur[1].y * invHauteur[1].y +
            16 * invHauteur[0].x * invHauteur[0].y * invHauteur[1].x * invHauteur[1].y +
            4 * invHauteur[0].y * invHauteur[0].y * invHauteur[1].x * invHauteur[1].x) +
         (128. / 3.) * DOFt[11] *
           (12 * invHauteur[0].x * invHauteur[1].x * invHauteur[1].y * invHauteur[1].y +
            12 * invHauteur[0].y * invHauteur[1].x * invHauteur[1].x * invHauteur[1].y) +

         128. * DOFt[12] *
           (4 * invHauteur[0].x * invHauteur[1].x * invHauteur[2].y * invHauteur[2].y +
            4 * invHauteur[0].y * invHauteur[1].y * invHauteur[2].x * invHauteur[2].x +
            8 * invHauteur[0].x * invHauteur[1].y * invHauteur[2].x * invHauteur[2].y +
            8 * invHauteur[0].y * invHauteur[1].x * invHauteur[2].x * invHauteur[2].y) +
         128. * DOFt[13] *
           (4 * invHauteur[2].x * invHauteur[0].x * invHauteur[1].y * invHauteur[1].y +
            4 * invHauteur[2].y * invHauteur[0].y * invHauteur[1].x * invHauteur[1].x +
            8 * invHauteur[2].x * invHauteur[0].y * invHauteur[1].x * invHauteur[1].y +
            8 * invHauteur[2].y * invHauteur[0].x * invHauteur[1].x * invHauteur[1].y) +
         128. * DOFt[14] *
           (4 * invHauteur[1].x * invHauteur[2].x * invHauteur[0].y * invHauteur[0].y +
            4 * invHauteur[1].y * invHauteur[2].y * invHauteur[0].x * invHauteur[0].x +
            8 * invHauteur[1].x * invHauteur[2].y * invHauteur[0].x * invHauteur[0].y +
            8 * invHauteur[1].y * invHauteur[2].x * invHauteur[0].x * invHauteur[0].y);

  f[3] = (32. / 3.) * DOFt[0] * 24 * invHauteur[0].y * invHauteur[0].y * invHauteur[0].y *
           invHauteur[0].x +
         (32. / 3.) * DOFt[1] * 24 * invHauteur[1].y * invHauteur[1].y * invHauteur[1].y *
           invHauteur[1].x +
         (32. / 3.) * DOFt[2] * 24 * invHauteur[2].y * invHauteur[2].y * invHauteur[2].y *
           invHauteur[2].x +

         (128. / 3.) * DOFt[3] *
           (6 * invHauteur[1].y * invHauteur[1].y * invHauteur[1].y * invHauteur[2].x +
            18 * invHauteur[1].y * invHauteur[1].y * invHauteur[1].x * invHauteur[2].y) +
         64 * DOFt[4] *
           (12 * invHauteur[1].y * invHauteur[1].y * invHauteur[2].y * invHauteur[2].x +
            12 * invHauteur[1].y * invHauteur[1].x * invHauteur[2].y * invHauteur[2].y) +
         (128. / 3.) * DOFt[5] *
           (18 * invHauteur[1].y * invHauteur[2].y * invHauteur[2].y * invHauteur[2].x +
            6 * invHauteur[1].x * invHauteur[2].y * invHauteur[2].y * invHauteur[2].y) +

         (128. / 3.) * DOFt[6] *
           (6 * invHauteur[2].y * invHauteur[2].y * invHauteur[2].y * invHauteur[0].x +
            18 * invHauteur[2].y * invHauteur[2].y * invHauteur[2].x * invHauteur[0].y) +
         64 * DOFt[7] *
           (12 * invHauteur[2].y * invHauteur[2].y * invHauteur[0].y * invHauteur[0].x +
            12 * invHauteur[2].y * invHauteur[2].x * invHauteur[0].y * invHauteur[0].y) +
         (128. / 3.) * DOFt[8] *
           (18 * invHauteur[2].y * invHauteur[0].y * invHauteur[0].y * invHauteur[0].x +
            6 * invHauteur[2].x * invHauteur[0].y * invHauteur[0].y * invHauteur[0].y) +

         (128. / 3.) * DOFt[9] *
           (6 * invHauteur[0].y * invHauteur[0].y * invHauteur[0].y * invHauteur[1].x +
            18 * invHauteur[0].y * invHauteur[0].y * invHauteur[0].x * invHauteur[1].y) +
         64 * DOFt[10] *
           (12 * invHauteur[0].y * invHauteur[0].y * invHauteur[1].y * invHauteur[1].x +
            12 * invHauteur[0].y * invHauteur[0].x * invHauteur[1].y * invHauteur[1].y) +
         (128. / 3.) * DOFt[11] *
           (18 * invHauteur[0].y * invHauteur[1].y * invHauteur[1].y * invHauteur[1].x +
            6 * invHauteur[0].x * invHauteur[1].y * invHauteur[1].y * invHauteur[1].y) +

         128. * DOFt[12] *
           (12 * invHauteur[0].y * invHauteur[1].y * invHauteur[2].y * invHauteur[2].x +
            6 * invHauteur[0].y * invHauteur[1].x * invHauteur[2].y * invHauteur[2].y +
            6 * invHauteur[0].x * invHauteur[1].y * invHauteur[2].y * invHauteur[2].y) +
         128. * DOFt[13] *
           (12 * invHauteur[2].y * invHauteur[0].y * invHauteur[1].y * invHauteur[1].x +
            6 * invHauteur[2].y * invHauteur[0].x * invHauteur[1].y * invHauteur[1].y +
            6 * invHauteur[2].x * invHauteur[0].y * invHauteur[1].y * invHauteur[1].y) +
         128. * DOFt[14] *
           (12 * invHauteur[1].y * invHauteur[2].y * invHauteur[0].y * invHauteur[0].x +
            6 * invHauteur[1].y * invHauteur[2].x * invHauteur[0].y * invHauteur[0].y +
            6 * invHauteur[1].x * invHauteur[2].y * invHauteur[0].y * invHauteur[0].y);
  f[4] =
    (32. / 3.) * DOFt[0] * 24 * invHauteur[0].y * invHauteur[0].y * invHauteur[0].y *
      invHauteur[0].y +
    (32. / 3.) * DOFt[1] * 24 * invHauteur[1].y * invHauteur[1].y * invHauteur[1].y *
      invHauteur[1].y +
    (32. / 3.) * DOFt[2] * 24 * invHauteur[2].y * invHauteur[2].y * invHauteur[2].y *
      invHauteur[2].y +
    (128. / 3.) * DOFt[3] * 24 * invHauteur[1].y * invHauteur[1].y * invHauteur[1].y *
      invHauteur[2].y +
    64 * DOFt[4] * 24 * invHauteur[1].y * invHauteur[1].y * invHauteur[2].y * invHauteur[2].y +
    (128. / 3.) * DOFt[5] * 24 * invHauteur[1].y * invHauteur[2].y * invHauteur[2].y *
      invHauteur[2].y +
    (128. / 3.) * DOFt[6] * 24 * invHauteur[2].y * invHauteur[2].y * invHauteur[2].y *
      invHauteur[0].y +
    64 * DOFt[7] * 24 * invHauteur[2].y * invHauteur[2].y * invHauteur[0].y * invHauteur[0].y +
    (128. / 3.) * DOFt[8] * 24 * invHauteur[2].y * invHauteur[0].y * invHauteur[0].y *
      invHauteur[0].y +
    (128. / 3.) * DOFt[9] * 24 * invHauteur[0].y * invHauteur[0].y * invHauteur[0].y *
      invHauteur[1].y +
    64 * DOFt[10] * 24 * invHauteur[0].y * invHauteur[0].y * invHauteur[1].y * invHauteur[1].y +
    (128. / 3.) * DOFt[11] * 24 * invHauteur[0].y * invHauteur[1].y * invHauteur[1].y *
      invHauteur[1].y +
    128. * DOFt[12] * 24 * invHauteur[0].y * invHauteur[1].y * invHauteur[2].y * invHauteur[2].y +
    128. * DOFt[13] * 24 * invHauteur[2].y * invHauteur[0].y * invHauteur[1].y * invHauteur[1].y +
    128. * DOFt[14] * 24 * invHauteur[1].y * invHauteur[2].y * invHauteur[0].y * invHauteur[0].y;
}

void TensorK::getDerivatives(const std::vector< double > &DOFt, const R2 invHauteur[3],
                             double *f) const {
  assert(DOFt.size( ) == ((m_deg + 1) * m_deg) / 2);

  // [((m+1)*(m+2))/2]
  switch (m_deg) {
    case 2:
      return Derivatives< 2 >(DOFt, invHauteur, f);
    case 3:
      return Derivatives< 3 >(DOFt, invHauteur, f);
    case 4:
      return Derivatives< 4 >(DOFt, invHauteur, f);
    case 5:
      return Derivatives< 5 >(DOFt, invHauteur, f);
  }
}

#endif
