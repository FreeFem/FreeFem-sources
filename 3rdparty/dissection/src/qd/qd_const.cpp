/*
 * src/qd_const.cc
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2000-2001
 *
 * Defines constants used in quad-double package.
 */
#include "config.h"
#include <qd/qd_real.h>

/* Some useful constants. */
const qd_real qd_real::_2pi = qd_real(6.283185307179586232e+00,
                                      2.449293598294706414e-16,
                                      -5.989539619436679332e-33,
                                      2.224908441726730563e-49);
const qd_real qd_real::_pi = qd_real(3.141592653589793116e+00,
                                     1.224646799147353207e-16,
                                     -2.994769809718339666e-33,
                                     1.112454220863365282e-49);
const qd_real qd_real::_pi2 = qd_real(1.570796326794896558e+00,
                                      6.123233995736766036e-17,
                                      -1.497384904859169833e-33,
                                      5.562271104316826408e-50);
const qd_real qd_real::_pi4 = qd_real(7.853981633974482790e-01,
                                      3.061616997868383018e-17,
                                      -7.486924524295849165e-34,
                                      2.781135552158413204e-50);
const qd_real qd_real::_3pi4 = qd_real(2.356194490192344837e+00,
                                       9.1848509936051484375e-17,
                                       3.9168984647504003225e-33,
                                      -2.5867981632704860386e-49);
const qd_real qd_real::_e = qd_real(2.718281828459045091e+00,
                                    1.445646891729250158e-16,
                                    -2.127717108038176765e-33,
                                    1.515630159841218954e-49);
const qd_real qd_real::_log2 = qd_real(6.931471805599452862e-01,
                                       2.319046813846299558e-17,
                                       5.707708438416212066e-34,
                                       -3.582432210601811423e-50);
const qd_real qd_real::_log10 = qd_real(2.302585092994045901e+00,
                                        -2.170756223382249351e-16,
                                        -9.984262454465776570e-33,
                                        -4.023357454450206379e-49);
const qd_real qd_real::_nan = qd_real(qd::_d_nan, qd::_d_nan, 
                                      qd::_d_nan, qd::_d_nan);
const qd_real qd_real::_inf = qd_real(qd::_d_inf, qd::_d_inf, 
                                      qd::_d_inf, qd::_d_inf);

const double qd_real::_eps = 1.21543267145725e-63; // = 2^-209
const double qd_real::_min_normalized = 1.6259745436952323e-260; // = 2^(-1022 + 3*53)
const qd_real qd_real::_max = qd_real(
    1.79769313486231570815e+308, 9.97920154767359795037e+291, 
    5.53956966280111259858e+275, 3.07507889307840487279e+259);
const qd_real qd_real::_safe_max = qd_real(
    1.7976931080746007281e+308,  9.97920154767359795037e+291, 
    5.53956966280111259858e+275, 3.07507889307840487279e+259);
const int qd_real::_ndigits = 62;

