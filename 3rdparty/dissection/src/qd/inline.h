/*
 * include/inline.h
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2000-2001
 *
 * This file contains the basic functions used both by double-double
 * and quad-double package.  These are declared as inline functions as
 * they are the smallest building blocks of the double-double and 
 * quad-double arithmetic.
 */
#ifndef _QD_INLINE_H
#define _QD_INLINE_H

#define _QD_SPLITTER 134217729.0               // = 2^27 + 1
#define _QD_SPLIT_THRESH 6.69692879491417e+299 // = 2^996

#ifdef QD_VACPP_BUILTINS_H
/* For VisualAge C++ __fmadd */
#include <builtins.h>
#endif

#include <cmath>
#include <limits>

namespace qd {

static const double _d_nan = std::numeric_limits<double>::quiet_NaN();
static const double _d_inf = std::numeric_limits<double>::infinity();

/*********** Basic Functions ************/
/* Computes fl(a+b) and err(a+b).  Assumes |a| >= |b|. */
inline double quick_two_sum(double a, double b, volatile double &err) {
  volatile double s = a + b;
  volatile double ss = s - a;
  err = b - ss;
  //  err = b - (s - a);
  return s;
}

/* Computes fl(a-b) and err(a-b).  Assumes |a| >= |b| */
inline double quick_two_diff(double a, double b, volatile double &err) {
  volatile double s = a - b;
  volatile double ss = a - s;
  err = ss - b;
  //  err = (a - s) - b;
  return s;
}

/* Computes fl(a+b) and err(a+b).  */
inline double two_sum(double a, double b, volatile double &err) {
  volatile double s = a + b;
  volatile double bb = s - a;
  volatile double sbb = s - bb;
  volatile double bbb = b - bb;
  err = (a - sbb) + bbb;
//  err = (a - (s - bb)) + (b - bb);
  return s;
}

/* Computes fl(a-b) and err(a-b).  */
inline double two_diff(double a, double b, volatile double &err) {
  volatile double s = a - b;
  volatile double bb = s - a;
  volatile double sbb = s - bb;
  volatile double bbb = b + bb;
  err = (a - sbb) - bbb;
  return s;
}

#ifndef QD_FMS
/* Computes high word and lo word of a */
inline void split(double a, double &hi, double &lo) {
  volatile double temp;
  if (a > _QD_SPLIT_THRESH || a < -_QD_SPLIT_THRESH) {
    a *= 3.7252902984619140625e-09;  // 2^-28
    temp = _QD_SPLITTER * a;
    volatile double tempa = temp - a;
    hi = temp - tempa;
    lo = a - hi;
    hi *= 268435456.0;          // 2^28
    lo *= 268435456.0;          // 2^28
  } else {
    temp = _QD_SPLITTER * a;
    volatile double tempa = temp - a;
    hi = temp - tempa;
    lo = a - hi;
  }
}
#endif

/* Computes fl(a*b) and err(a*b). */
inline double two_prod(double a, double b, volatile double &err) {
#ifdef QD_FMS
  volatile double p = a * b;
  err = QD_FMS(a, b, p);
  return p;
#else
  double a_hi, a_lo, b_hi, b_lo;
  volatile double p = a * b;
  split(a, a_hi, a_lo);
  split(b, b_hi, b_lo);
  err = ((a_hi * b_hi - p) + a_hi * b_lo + a_lo * b_hi) + a_lo * b_lo;
  return p;
#endif
}

/* Computes fl(a*a) and err(a*a).  Faster than the above method. */
inline double two_sqr(double a, volatile double &err) {
#ifdef QD_FMS
  volatile double p = a * a;
  err = QD_FMS(a, a, p);
  return p;
#else
  double hi, lo;
  volatile double q = a * a;
  split(a, hi, lo);
  err = ((hi * hi - q) + 2.0 * hi * lo) + lo * lo;
  return q;
#endif
}

/* Computes the nearest integer to d. */
inline double nint(double d) {
  if (d == std::floor(d))
    return d;
  return std::floor(d + 0.5);
}

/* Computes the truncated integer. */
inline double aint(double d) {
  return (d >= 0.0) ? std::floor(d) : std::ceil(d);
}

/* These are provided to give consistent 
   interface for double with double-double and quad-double. */
inline void sincosh(double t, double &sinh_t, double &cosh_t) {
  sinh_t = std::sinh(t);
  cosh_t = std::cosh(t);
}

inline double sqr(double t) {
  return t * t;
}

inline double to_double(double a) { return a; }
inline int    to_int(double a) { return static_cast<int>(a); }

}

#endif /* _QD_INLINE_H */
