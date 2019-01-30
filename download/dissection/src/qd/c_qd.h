/*
 * include/c_qd.h
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2000-2001
 *
 * Contains C wrapper function prototypes for quad-double precision 
 * arithmetic.  This can also be used from fortran code.
 */
#ifndef _QD_C_QD_H
#define _QD_C_QD_H

#include <qd/c_dd.h>
#include <qd/qd_config.h>

#ifdef __cplusplus
extern "C" {
#endif

/* add */
void c_qd_add(const double *a, const double *b, double *c);
void c_qd_add_dd_qd(const double *a, const double *b, double *c);
void c_qd_add_qd_dd(const double *a, const double *b, double *c);
void c_qd_add_d_qd(double a, const double *b, double *c);
void c_qd_add_qd_d(const double *a, double b, double *c);
void c_qd_selfadd(const double *a, double *b);
void c_qd_selfadd_dd(const double *a, double *b);
void c_qd_selfadd_d(double a, double *b);

/* sub */
void c_qd_sub(const double *a, const double *b, double *c);
void c_qd_sub_dd_qd(const double *a, const double *b, double *c);
void c_qd_sub_qd_dd(const double *a, const double *b, double *c);
void c_qd_sub_d_qd(double a, const double *b, double *c);
void c_qd_sub_qd_d(const double *a, double b, double *c);
void c_qd_selfsub(const double *a, double *b);
void c_qd_selfsub_dd(const double *a, double *b);
void c_qd_selfsub_d(double a, double *b);

/* mul */
void c_qd_mul(const double *a, const double *b, double *c);
void c_qd_mul_dd_qd(const double *a, const double *b, double *c);
void c_qd_mul_qd_dd(const double *a, const double *b, double *c);
void c_qd_mul_d_qd(double a, const double *b, double *c);
void c_qd_mul_qd_d(const double *a, double b, double *c);
void c_qd_selfmul(const double *a, double *b);
void c_qd_selfmul_dd(const double *a, double *b);
void c_qd_selfmul_d(double a, double *b);

/* div */
void c_qd_div(const double *a, const double *b, double *c);
void c_qd_div_dd_qd(const double *a, const double *b, double *c);
void c_qd_div_qd_dd(const double *a, const double *b, double *c);
void c_qd_div_d_qd(double a, const double *b, double *c);
void c_qd_div_qd_d(const double *a, double b, double *c);
void c_qd_selfdiv(const double *a, double *b);
void c_qd_selfdiv_dd(const double *a, double *b);
void c_qd_selfdiv_d(double a, double *b);

/* copy */
void c_qd_copy(const double *a, double *b);
void c_qd_copy_dd(const double *a, double *b);
void c_qd_copy_d(double a, double *b);

void c_qd_sqrt(const double *a, double *b);
void c_qd_sqr(const double *a, double *b);

void c_qd_abs(const double *a, double *b);

void c_qd_npwr(const double *a, int b, double *c);
void c_qd_nroot(const double *a, int b, double *c);

void c_qd_nint(const double *a, double *b);
void c_qd_aint(const double *a, double *b);
void c_qd_floor(const double *a, double *b);
void c_qd_ceil(const double *a, double *b);

void c_qd_exp(const double *a, double *b);
void c_qd_log(const double *a, double *b);
void c_qd_log10(const double *a, double *b);

void c_qd_sin(const double *a, double *b);
void c_qd_cos(const double *a, double *b);
void c_qd_tan(const double *a, double *b);

void c_qd_asin(const double *a, double *b);
void c_qd_acos(const double *a, double *b);
void c_qd_atan(const double *a, double *b);
void c_qd_atan2(const double *a, const double *b, double *c);

void c_qd_sinh(const double *a, double *b);
void c_qd_cosh(const double *a, double *b);
void c_qd_tanh(const double *a, double *b);

void c_qd_asinh(const double *a, double *b);
void c_qd_acosh(const double *a, double *b);
void c_qd_atanh(const double *a, double *b);

void c_qd_sincos(const double *a, double *s, double *c);
void c_qd_sincosh(const double *a, double *s, double *c);

void c_qd_read(const char *s, double *a);
void c_qd_swrite(const double *a, int precision, char *s, int len);
void c_qd_write(const double *a);
void c_qd_neg(const double *a, double *b);
void c_qd_rand(double *a);
void c_qd_comp(const double *a, const double *b, int *result);
void c_qd_comp_qd_d(const double *a, double b, int *result);
void c_qd_comp_d_qd(double a, const double *b, int *result);
void c_qd_pi(double *a);

#ifdef __cplusplus
}
#endif

#endif  /* _QD_C_QD_H */
