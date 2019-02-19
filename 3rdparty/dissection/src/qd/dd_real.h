/*
 * include/dd_real.h
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2000-2007
 *
 * Double-double precision (>= 106-bit significand) floating point
 * arithmetic package based on David Bailey's Fortran-90 double-double
 * package, with some changes. See  
 *
 *   http://www.nersc.gov/~dhbailey/mpdist/mpdist.html
 *   
 * for the original Fortran-90 version.
 *
 * Overall structure is similar to that of Keith Brigg's C++ double-double
 * package.  See  
 *
 *   http://www-epidem.plansci.cam.ac.uk/~kbriggs/doubledouble.html
 *
 * for more details.  In particular, the fix for x86 computers is borrowed
 * from his code.
 *
 * Yozo Hida
 */
// operator int(), dd_real copysign(), fmax(), lobg(), scalbn() are added
// for complex class of LLVM Clang++ : 23 Jul.2015 Atsushi Suzuki

#ifndef _QD_DD_REAL_H
#define _QD_DD_REAL_H

#include <cmath>
#include <iostream>
#include <string>
#include <limits>
#include <qd/qd_config.h>
#include <qd/fpu.h>

// Some compilers define isnan, isfinite, and isinf as macros, even for
// C++ codes, which cause havoc when overloading these functions.  We undef
// them here.
#ifdef isnan
#undef isnan
#endif

#ifdef isfinite
#undef isfinite
#endif

#ifdef isinf
#undef isinf
#endif

#ifdef max
#undef max
#endif

#ifdef min
#undef min
#endif

struct QD_API dd_real {
  double x[2];

  dd_real(double hi, double lo) { x[0] = hi; x[1] = lo; }
  dd_real() {x[0] = 0.0; x[1] = 0.0; }
  dd_real(double h) { x[0] = h; x[1] = 0.0; }
  dd_real(int h) {
    x[0] = (static_cast<double>(h));
    x[1] = 0.0;
  }

  dd_real (const char *s);
  explicit dd_real (const double *d) {
    x[0] = d[0]; x[1] = d[1];
  }

  static void error(const char *msg);

  double _hi() const { return x[0]; }
  double _lo() const { return x[1]; }

  static const dd_real _2pi;
  static const dd_real _pi;
  static const dd_real _3pi4;
  static const dd_real _pi2;
  static const dd_real _pi4;
  static const dd_real _e;
  static const dd_real _log2;
  static const dd_real _log10;
  static const dd_real _nan;
  static const dd_real _inf;

  static const double _eps;
  static const double _min_normalized;
  static const dd_real _max;
  static const dd_real _safe_max;
  static const int _ndigits;

  bool isnan() const { return QD_ISNAN(x[0]) || QD_ISNAN(x[1]); }
  bool isfinite() const { return QD_ISFINITE(x[0]); }
  bool isinf() const { return QD_ISINF(x[0]); }

  static dd_real add(double a, double b);
  static dd_real ieee_add(const dd_real &a, const dd_real &b);
  static dd_real sloppy_add(const dd_real &a, const dd_real &b);

  dd_real &operator+=(double a);
  dd_real &operator+=(const dd_real &a);

  static dd_real sub(double a, double b);

  dd_real &operator-=(double a);
  dd_real &operator-=(const dd_real &a);

  dd_real operator-() const;
  dd_real operator+() const;

  static dd_real mul(double a, double b);

  dd_real &operator*=(double a);
  dd_real &operator*=(const dd_real &a);

  static dd_real div(double a, double b);
  static dd_real sloppy_div(const dd_real &a, const dd_real &b);
  static dd_real accurate_div(const dd_real &a, const dd_real &b);
  
  dd_real &operator/=(double a);
  dd_real &operator/=(const dd_real &a);

  dd_real &operator=(double a);
  dd_real &operator=(const char *s);

  dd_real operator^(int n);
  static dd_real sqr(double d);

  static dd_real sqrt(double a);
  
  bool is_zero() const;
  bool is_one() const;
  bool is_positive() const;
  bool is_negative() const;

  static dd_real rand(void);

  void to_digits(char *s, int &expn, int precision = _ndigits) const;
  void write(char *s, int len, int precision = _ndigits, 
      bool showpos = false, bool uppercase = false) const;
  std::string to_string(int precision = _ndigits, int width = 0, 
      std::ios_base::fmtflags fmt = static_cast<std::ios_base::fmtflags>(0), 
      bool showpos = false, bool uppercase = false, char fill = ' ') const;
  int read(const char *s, dd_real &a);

  /* Debugging Methods */
  void dump(const std::string &name = "", std::ostream &os = std::cerr) const;
  void dump_bits(const std::string &name = "", 
                 std::ostream &os = std::cerr) const;

  static dd_real debug_rand();

  // added for complex class of LLVM Clang++: 23 Jul.2015 Atsushi Suzuki
  operator int() { return int(x[0]); }
  operator int() const { return int(x[0]); }
};


namespace std {
  template <>
  class numeric_limits<dd_real> : public numeric_limits<double> {
  public:
    inline static double epsilon() { return dd_real::_eps; }
    inline static dd_real max() { return dd_real::_max; }
    inline static dd_real safe_max() { return dd_real::_safe_max; }
    inline static double min() { return dd_real::_min_normalized; }
    static const int digits = 104;
    static const int digits10 = 31;
  };
}

QD_API dd_real ddrand(void);
QD_API dd_real sqrt(const dd_real &a);

QD_API dd_real polyeval(const dd_real *c, int n, const dd_real &x);
QD_API dd_real polyroot(const dd_real *c, int n, 
    const dd_real &x0, int max_iter = 32, double thresh = 0.0);

QD_API inline bool isnan(const dd_real &a) { return a.isnan(); }
QD_API inline bool isfinite(const dd_real &a) { return a.isfinite(); }
QD_API inline bool isinf(const dd_real &a) { return a.isinf(); }

/* Computes  dd * d  where d is known to be a power of 2. */
QD_API dd_real mul_pwr2(const dd_real &dd, double d);

QD_API dd_real operator+(const dd_real &a, double b);
QD_API dd_real operator+(double a, const dd_real &b);
QD_API dd_real operator+(const dd_real &a, const dd_real &b);

QD_API dd_real operator-(const dd_real &a, double b);
QD_API dd_real operator-(double a, const dd_real &b);
QD_API dd_real operator-(const dd_real &a, const dd_real &b);

QD_API dd_real operator*(const dd_real &a, double b);
QD_API dd_real operator*(double a, const dd_real &b);
QD_API dd_real operator*(const dd_real &a, const dd_real &b);

QD_API dd_real operator/(const dd_real &a, double b);
QD_API dd_real operator/(double a, const dd_real &b);
QD_API dd_real operator/(const dd_real &a, const dd_real &b);

QD_API dd_real inv(const dd_real &a);

QD_API dd_real rem(const dd_real &a, const dd_real &b);
QD_API dd_real drem(const dd_real &a, const dd_real &b);
QD_API dd_real divrem(const dd_real &a, const dd_real &b, dd_real &r);

QD_API dd_real pow(const dd_real &a, int n);
QD_API dd_real pow(const dd_real &a, const dd_real &b);
QD_API dd_real npwr(const dd_real &a, int n);
QD_API dd_real sqr(const dd_real &a);

QD_API dd_real sqrt(const dd_real &a);
QD_API dd_real nroot(const dd_real &a, int n);

QD_API bool operator==(const dd_real &a, int b);
QD_API bool operator==(int a, const dd_real &b);
QD_API bool operator==(const dd_real &a, double b);
QD_API bool operator==(double a, const dd_real &b);
QD_API bool operator==(const dd_real &a, const dd_real &b);

QD_API bool operator<=(const dd_real &a, int b);
QD_API bool operator<=(int a, const dd_real &b);
QD_API bool operator<=(const dd_real &a, double b);
QD_API bool operator<=(double a, const dd_real &b);
QD_API bool operator<=(const dd_real &a, const dd_real &b);

QD_API bool operator>=(const dd_real &a, int b);
QD_API bool operator>=(int a, const dd_real &b);
QD_API bool operator>=(const dd_real &a, double b);
QD_API bool operator>=(double a, const dd_real &b);
QD_API bool operator>=(const dd_real &a, const dd_real &b);

QD_API bool operator<(const dd_real &a, int b);
QD_API bool operator<(int a, const dd_real &b);
QD_API bool operator<(const dd_real &a, double b);
QD_API bool operator<(double a, const dd_real &b);
QD_API bool operator<(const dd_real &a, const dd_real &b);

QD_API bool operator>(const dd_real &a, int b);
QD_API bool operator>(int a, const dd_real &b);
QD_API bool operator>(const dd_real &a, double b);
QD_API bool operator>(double a, const dd_real &b);
QD_API bool operator>(const dd_real &a, const dd_real &b);

QD_API bool operator!=(const dd_real &a, int b);
QD_API bool operator!=(int a, const dd_real &b);
QD_API bool operator!=(const dd_real &a, double b);
QD_API bool operator!=(double a, const dd_real &b);
QD_API bool operator!=(const dd_real &a, const dd_real &b);

QD_API dd_real nint(const dd_real &a);
QD_API dd_real floor(const dd_real &a);
QD_API dd_real ceil(const dd_real &a);
QD_API dd_real aint(const dd_real &a);

QD_API dd_real ddrand(void);

double to_double(const dd_real &a);
int    to_int(const dd_real &a);

// QD_API dd_real exp(const dd_real &a);
QD_API dd_real ldexp(const dd_real &a, int exp);
QD_API dd_real log(const dd_real &a);
QD_API dd_real log10(const dd_real &a);

QD_API dd_real sin(const dd_real &a);
QD_API dd_real cos(const dd_real &a);
QD_API dd_real tan(const dd_real &a);
QD_API void sincos(const dd_real &a, dd_real &sin_a, dd_real &cos_a);

QD_API dd_real asin(const dd_real &a);
QD_API dd_real acos(const dd_real &a);
QD_API dd_real atan(const dd_real &a);
QD_API dd_real atan2(const dd_real &y, const dd_real &x);

QD_API dd_real sinh(const dd_real &a);
QD_API dd_real cosh(const dd_real &a);
QD_API dd_real tanh(const dd_real &a);
QD_API void sincosh(const dd_real &a, 
                      dd_real &sinh_a, dd_real &cosh_a);

QD_API dd_real asinh(const dd_real &a);
QD_API dd_real acosh(const dd_real &a);
QD_API dd_real atanh(const dd_real &a);

QD_API dd_real fabs(const dd_real &a);
QD_API dd_real abs(const dd_real &a);   /* same as fabs */

QD_API dd_real fmod(const dd_real &a, const dd_real &b);

QD_API std::ostream& operator<<(std::ostream &s, const dd_real &a);
QD_API std::istream& operator>>(std::istream &s, dd_real &a);
#ifdef QD_INLINE
#include <qd/dd_inline.h>
#endif

// added for complex class of LLVM Clang++ : 23 Jul.2015 Atsushi Suzuki
inline dd_real copysign(const dd_real &x, const dd_real &y) {
  return (y.x[0] < 0.0) ? ((x.x[0] < 0.0) ? x : (-x)) : ((x.x[0] < 0.0) ? (-x) : x); 
}

inline dd_real fmax(const dd_real &x, const dd_real &y) {
  return x.x[0] < y.x[0] ? y : x;
}

inline dd_real logb(const dd_real &y) {
  return dd_real(logb(y.x[0]));
}
#if 0
inline dd_real scalbn(const dd_real &x, int n) {
  return dd_real(scalb(x.x[0], n));
}
#endif
#endif /* _QD_DD_REAL_H */

