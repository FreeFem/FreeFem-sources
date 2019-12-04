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
// AUTHORS : ...
// E-MAIL  : ...

/*
 * //  -- missing type "const double[]"
 * missing: 5   gsl_ran_dirichlet_pdf  -> /usr/local/include/gsl/gsl_randist.h:double
 * gsl_ran_dirichlet_pdf (const size_t K, const double alpha[], const double theta[]); missing: 5
 * gsl_ran_dirichlet_lnpdf  -> /usr/local/include/gsl/gsl_randist.h:double gsl_ran_dirichlet_lnpdf
 * (const size_t K, const double alpha[], const double theta[]);
 * //  -- missing type "size_t"
 * missing: 4   gsl_ran_discrete_pdf  -> /usr/local/include/gsl/gsl_randist.h:double
 * gsl_ran_discrete_pdf (size_t k, const gsl_ran_discrete_t *g);
 * //  -- missing type "const gsl_mode_t"
 * missing: 5   gsl_sf_airy_Ai_e  -> /usr/local/include/gsl/gsl_sf_airy.h:int gsl_sf_airy_Ai_e(const
 * double x, const gsl_mode_t mode, gsl_sf_result * result);
 * //  -- missing type "gsl_sf_result *"
 * missing: 5   gsl_sf_airy_Bi_e  -> /usr/local/include/gsl/gsl_sf_airy.h:int gsl_sf_airy_Bi_e(const
 * double x, gsl_mode_t mode, gsl_sf_result * result); missing: 5   gsl_sf_airy_Ai_scaled_e  ->
 * /usr/local/include/gsl/gsl_sf_airy.h:int gsl_sf_airy_Ai_scaled_e(const double x, gsl_mode_t mode,
 * gsl_sf_result * result); missing: 5   gsl_sf_airy_Bi_scaled_e  ->
 * /usr/local/include/gsl/gsl_sf_airy.h:int gsl_sf_airy_Bi_scaled_e(const double x, gsl_mode_t mode,
 * gsl_sf_result * result); missing: 5   gsl_sf_airy_Ai_deriv_e  ->
 * /usr/local/include/gsl/gsl_sf_airy.h:int gsl_sf_airy_Ai_deriv_e(const double x, gsl_mode_t mode,
 * gsl_sf_result * result); missing: 5   gsl_sf_airy_Bi_deriv_e  ->
 * /usr/local/include/gsl/gsl_sf_airy.h:int gsl_sf_airy_Bi_deriv_e(const double x, gsl_mode_t mode,
 * gsl_sf_result * result); missing: 5   gsl_sf_airy_Ai_deriv_scaled_e  ->
 * /usr/local/include/gsl/gsl_sf_airy.h:int gsl_sf_airy_Ai_deriv_scaled_e(const double x, gsl_mode_t
 * mode, gsl_sf_result * result); missing: 5   gsl_sf_airy_Bi_deriv_scaled_e  ->
 * /usr/local/include/gsl/gsl_sf_airy.h:int gsl_sf_airy_Bi_deriv_scaled_e(const double x, gsl_mode_t
 * mode, gsl_sf_result * result); missing: 4   gsl_sf_airy_zero_Ai_e  ->
 * /usr/local/include/gsl/gsl_sf_airy.h:int gsl_sf_airy_zero_Ai_e(unsigned int s, gsl_sf_result *
 * result); missing: 4   gsl_sf_airy_zero_Bi_e  -> /usr/local/include/gsl/gsl_sf_airy.h:int
 * gsl_sf_airy_zero_Bi_e(unsigned int s, gsl_sf_result * result); missing: 4
 * gsl_sf_airy_zero_Ai_deriv_e  -> /usr/local/include/gsl/gsl_sf_airy.h:int
 * gsl_sf_airy_zero_Ai_deriv_e(unsigned int s, gsl_sf_result * result); missing: 4
 * gsl_sf_airy_zero_Bi_deriv_e  -> /usr/local/include/gsl/gsl_sf_airy.h:int
 * gsl_sf_airy_zero_Bi_deriv_e(unsigned int s, gsl_sf_result * result); missing: 4
 * gsl_sf_bessel_J0_e  -> /usr/local/include/gsl/gsl_sf_bessel.h:int gsl_sf_bessel_J0_e(const double
 * x,  gsl_sf_result * result); missing: 4   gsl_sf_bessel_J1_e  ->
 * /usr/local/include/gsl/gsl_sf_bessel.h:int gsl_sf_bessel_J1_e(const double x, gsl_sf_result *
 * result); missing: 5   gsl_sf_bessel_Jn_e  -> /usr/local/include/gsl/gsl_sf_bessel.h:int
 * gsl_sf_bessel_Jn_e(int n, double x, gsl_sf_result * result); missing: 4   gsl_sf_bessel_Y0_e  ->
 * /usr/local/include/gsl/gsl_sf_bessel.h:int gsl_sf_bessel_Y0_e(const double x, gsl_sf_result *
 * result); missing: 4   gsl_sf_bessel_Y1_e  -> /usr/local/include/gsl/gsl_sf_bessel.h:int
 * gsl_sf_bessel_Y1_e(const double x, gsl_sf_result * result); missing: 5   gsl_sf_bessel_Yn_e  ->
 * /usr/local/include/gsl/gsl_sf_bessel.h:int gsl_sf_bessel_Yn_e(int n,const double x, gsl_sf_result
 * * result); missing: 4   gsl_sf_bessel_I0_e  -> /usr/local/include/gsl/gsl_sf_bessel.h:int
 * gsl_sf_bessel_I0_e(const double x, gsl_sf_result * result); missing: 4   gsl_sf_bessel_I1_e  ->
 * /usr/local/include/gsl/gsl_sf_bessel.h:int gsl_sf_bessel_I1_e(const double x, gsl_sf_result *
 * result); missing: 5   gsl_sf_bessel_In_e  -> /usr/local/include/gsl/gsl_sf_bessel.h:int
 * gsl_sf_bessel_In_e(const int n, const double x, gsl_sf_result * result); missing: 4
 * gsl_sf_bessel_I0_scaled_e  -> /usr/local/include/gsl/gsl_sf_bessel.h:int
 * gsl_sf_bessel_I0_scaled_e(const double x, gsl_sf_result * result); missing: 4
 * gsl_sf_bessel_I1_scaled_e  -> /usr/local/include/gsl/gsl_sf_bessel.h:int
 * gsl_sf_bessel_I1_scaled_e(const double x, gsl_sf_result * result); missing: 5
 * gsl_sf_bessel_In_scaled_e  -> /usr/local/include/gsl/gsl_sf_bessel.h:int
 * gsl_sf_bessel_In_scaled_e(int n, const double x, gsl_sf_result * result); missing: 4
 * gsl_sf_bessel_K0_e  -> /usr/local/include/gsl/gsl_sf_bessel.h:int gsl_sf_bessel_K0_e(const double
 * x, gsl_sf_result * result); missing: 4   gsl_sf_bessel_K1_e  ->
 * /usr/local/include/gsl/gsl_sf_bessel.h:int gsl_sf_bessel_K1_e(const double x, gsl_sf_result *
 * result); missing: 5   gsl_sf_bessel_Kn_e  -> /usr/local/include/gsl/gsl_sf_bessel.h:int
 * gsl_sf_bessel_Kn_e(const int n, const double x, gsl_sf_result * result); missing: 4
 * gsl_sf_bessel_K0_scaled_e  -> /usr/local/include/gsl/gsl_sf_bessel.h:int
 * gsl_sf_bessel_K0_scaled_e(const double x, gsl_sf_result * result); missing: 4
 * gsl_sf_bessel_K1_scaled_e  -> /usr/local/include/gsl/gsl_sf_bessel.h:int
 * gsl_sf_bessel_K1_scaled_e(const double x, gsl_sf_result * result); missing: 5
 * gsl_sf_bessel_Kn_scaled_e  -> /usr/local/include/gsl/gsl_sf_bessel.h:int
 * gsl_sf_bessel_Kn_scaled_e(int n, const double x, gsl_sf_result * result); missing: 4
 * gsl_sf_bessel_j0_e  -> /usr/local/include/gsl/gsl_sf_bessel.h:int gsl_sf_bessel_j0_e(const double
 * x, gsl_sf_result * result); missing: 4   gsl_sf_bessel_j1_e  ->
 * /usr/local/include/gsl/gsl_sf_bessel.h:int gsl_sf_bessel_j1_e(const double x, gsl_sf_result *
 * result); missing: 4   gsl_sf_bessel_j2_e  -> /usr/local/include/gsl/gsl_sf_bessel.h:int
 * gsl_sf_bessel_j2_e(const double x, gsl_sf_result * result); missing: 5   gsl_sf_bessel_jl_e  ->
 * /usr/local/include/gsl/gsl_sf_bessel.h:int gsl_sf_bessel_jl_e(const int l, const double x,
 * gsl_sf_result * result);
 * //  -- missing type "double *"
 * missing: 5   gsl_sf_bessel_jl_array  -> /usr/local/include/gsl/gsl_sf_bessel.h:int
 * gsl_sf_bessel_jl_array(const int lmax, const double x, double * result_array); missing: 5
 * gsl_sf_bessel_jl_steed_array  -> /usr/local/include/gsl/gsl_sf_bessel.h:int
 * gsl_sf_bessel_jl_steed_array(const int lmax, const double x, double * jl_x_array); missing: 4
 * gsl_sf_bessel_y0_e  -> /usr/local/include/gsl/gsl_sf_bessel.h:int gsl_sf_bessel_y0_e(const double
 * x, gsl_sf_result * result); missing: 4   gsl_sf_bessel_y1_e  ->
 * /usr/local/include/gsl/gsl_sf_bessel.h:int gsl_sf_bessel_y1_e(const double x, gsl_sf_result *
 * result); missing: 4   gsl_sf_bessel_y2_e  -> /usr/local/include/gsl/gsl_sf_bessel.h:int
 * gsl_sf_bessel_y2_e(const double x, gsl_sf_result * result); missing: 5   gsl_sf_bessel_yl_e  ->
 * /usr/local/include/gsl/gsl_sf_bessel.h:int gsl_sf_bessel_yl_e(int l, const double x,
 * gsl_sf_result * result); missing: 5   gsl_sf_bessel_yl_array  ->
 * /usr/local/include/gsl/gsl_sf_bessel.h:int gsl_sf_bessel_yl_array(const int lmax, const double x,
 * double * result_array); missing: 4   gsl_sf_bessel_i0_scaled_e  ->
 * /usr/local/include/gsl/gsl_sf_bessel.h:int gsl_sf_bessel_i0_scaled_e(const double x,
 * gsl_sf_result * result); missing: 4   gsl_sf_bessel_i1_scaled_e  ->
 * /usr/local/include/gsl/gsl_sf_bessel.h:int gsl_sf_bessel_i1_scaled_e(const double x,
 * gsl_sf_result * result); missing: 4   gsl_sf_bessel_i2_scaled_e  ->
 * /usr/local/include/gsl/gsl_sf_bessel.h:int gsl_sf_bessel_i2_scaled_e(const double x,
 * gsl_sf_result * result); missing: 5   gsl_sf_bessel_il_scaled_e  ->
 * /usr/local/include/gsl/gsl_sf_bessel.h:int gsl_sf_bessel_il_scaled_e(const int l, double x,
 * gsl_sf_result * result); missing: 5   gsl_sf_bessel_il_scaled_array  ->
 * /usr/local/include/gsl/gsl_sf_bessel.h:int gsl_sf_bessel_il_scaled_array(const int lmax, const
 * double x, double * result_array); missing: 4   gsl_sf_bessel_k0_scaled_e  ->
 * /usr/local/include/gsl/gsl_sf_bessel.h:int gsl_sf_bessel_k0_scaled_e(const double x,
 * gsl_sf_result * result); missing: 4   gsl_sf_bessel_k1_scaled_e  ->
 * /usr/local/include/gsl/gsl_sf_bessel.h:int gsl_sf_bessel_k1_scaled_e(const double x,
 * gsl_sf_result * result); missing: 4   gsl_sf_bessel_k2_scaled_e  ->
 * /usr/local/include/gsl/gsl_sf_bessel.h:int gsl_sf_bessel_k2_scaled_e(const double x,
 * gsl_sf_result * result); missing: 5   gsl_sf_bessel_kl_scaled_e  ->
 * /usr/local/include/gsl/gsl_sf_bessel.h:int gsl_sf_bessel_kl_scaled_e(int l, const double x,
 * gsl_sf_result * result); missing: 5   gsl_sf_bessel_kl_scaled_array  ->
 * /usr/local/include/gsl/gsl_sf_bessel.h:int gsl_sf_bessel_kl_scaled_array(const int lmax, const
 * double x, double * result_array); missing: 5   gsl_sf_bessel_Jnu_e  ->
 * /usr/local/include/gsl/gsl_sf_bessel.h:int gsl_sf_bessel_Jnu_e(const double nu, const double x,
 * gsl_sf_result * result); missing: 5   gsl_sf_bessel_Ynu_e  ->
 * /usr/local/include/gsl/gsl_sf_bessel.h:int gsl_sf_bessel_Ynu_e(double nu, double x, gsl_sf_result
 * * result); missing: 5   gsl_sf_bessel_Inu_scaled_e  -> /usr/local/include/gsl/gsl_sf_bessel.h:int
 * gsl_sf_bessel_Inu_scaled_e(double nu, double x, gsl_sf_result * result); missing: 5
 * gsl_sf_bessel_Inu_e  -> /usr/local/include/gsl/gsl_sf_bessel.h:int gsl_sf_bessel_Inu_e(double nu,
 * double x, gsl_sf_result * result); missing: 5   gsl_sf_bessel_Knu_scaled_e  ->
 * /usr/local/include/gsl/gsl_sf_bessel.h:int gsl_sf_bessel_Knu_scaled_e(const double nu, const
 * double x, gsl_sf_result * result);
 * //  -- missing type "gsl_sf_result_e10 *"
 * missing: 5   gsl_sf_bessel_Knu_scaled_e10_e  -> /usr/local/include/gsl/gsl_sf_bessel.h:int
 * gsl_sf_bessel_Knu_scaled_e10_e(const double nu, const double x, gsl_sf_result_e10 * result);
 * missing: 5   gsl_sf_bessel_Knu_e  -> /usr/local/include/gsl/gsl_sf_bessel.h:int
 * gsl_sf_bessel_Knu_e(const double nu, const double x, gsl_sf_result * result); missing: 5
 * gsl_sf_bessel_lnKnu_e  -> /usr/local/include/gsl/gsl_sf_bessel.h:int gsl_sf_bessel_lnKnu_e(const
 * double nu, const double x, gsl_sf_result * result); missing: 4   gsl_sf_bessel_zero_J0_e  ->
 * /usr/local/include/gsl/gsl_sf_bessel.h:int gsl_sf_bessel_zero_J0_e(unsigned int s, gsl_sf_result
 * * result); missing: 4   gsl_sf_bessel_zero_J1_e  -> /usr/local/include/gsl/gsl_sf_bessel.h:int
 * gsl_sf_bessel_zero_J1_e(unsigned int s, gsl_sf_result * result); missing: 5
 * gsl_sf_bessel_zero_Jnu_e  -> /usr/local/include/gsl/gsl_sf_bessel.h:int
 * gsl_sf_bessel_zero_Jnu_e(double nu, unsigned int s, gsl_sf_result * result); missing: 4
 * gsl_sf_clausen_e  -> /usr/local/include/gsl/gsl_sf_clausen.h:int gsl_sf_clausen_e(double x,
 * gsl_sf_result * result); missing: 5   gsl_sf_hydrogenicR_1_e  ->
 * /usr/local/include/gsl/gsl_sf_coulomb.h:int gsl_sf_hydrogenicR_1_e(const double Z, const double
 * r, gsl_sf_result * result); missing: 5   gsl_sf_coulomb_CL_e  ->
 * /usr/local/include/gsl/gsl_sf_coulomb.h:int gsl_sf_coulomb_CL_e(double L, double eta,
 * gsl_sf_result * result); missing: 6   gsl_sf_coulomb_CL_e  ->
 * /usr/local/include/gsl/gsl_sf_coulomb.h:int gsl_sf_coulomb_CL_array(double Lmin, int kmax, double
 * eta, double * cl); missing: 4   gsl_sf_dawson_e  -> /usr/local/include/gsl/gsl_sf_dawson.h:int
 * gsl_sf_dawson_e(double x, gsl_sf_result * result); missing: 4   gsl_sf_debye_1_e  ->
 * /usr/local/include/gsl/gsl_sf_debye.h:int     gsl_sf_debye_1_e(const double x, gsl_sf_result *
 * result); missing: 4   gsl_sf_debye_2_e  -> /usr/local/include/gsl/gsl_sf_debye.h:int
 * gsl_sf_debye_2_e(const double x, gsl_sf_result * result); missing: 4   gsl_sf_debye_3_e  ->
 * /usr/local/include/gsl/gsl_sf_debye.h:int     gsl_sf_debye_3_e(const double x, gsl_sf_result *
 * result); missing: 4   gsl_sf_debye_4_e  -> /usr/local/include/gsl/gsl_sf_debye.h:int
 * gsl_sf_debye_4_e(const double x, gsl_sf_result * result); missing: 4   gsl_sf_debye_5_e  ->
 * /usr/local/include/gsl/gsl_sf_debye.h:int     gsl_sf_debye_5_e(const double x, gsl_sf_result *
 * result); missing: 4   gsl_sf_debye_6_e  -> /usr/local/include/gsl/gsl_sf_debye.h:int
 * gsl_sf_debye_6_e(const double x, gsl_sf_result * result); missing: 4   gsl_sf_dilog_e  ->
 * /usr/local/include/gsl/gsl_sf_dilog.h:int     gsl_sf_dilog_e(const double x, gsl_sf_result *
 * result); missing: 5   gsl_sf_multiply_e  -> /usr/local/include/gsl/gsl_sf_elementary.h:int
 * gsl_sf_multiply_e(const double x, const double y, gsl_sf_result * result); missing: 5
 * gsl_sf_ellint_Kcomp_e  -> /usr/local/include/gsl/gsl_sf_ellint.h:int gsl_sf_ellint_Kcomp_e(double
 * k, gsl_mode_t mode, gsl_sf_result * result); missing: 5   gsl_sf_ellint_Ecomp_e  ->
 * /usr/local/include/gsl/gsl_sf_ellint.h:int gsl_sf_ellint_Ecomp_e(double k, gsl_mode_t mode,
 * gsl_sf_result * result); missing: 5   gsl_sf_ellint_Dcomp_e  ->
 * /usr/local/include/gsl/gsl_sf_ellint.h:int gsl_sf_ellint_Dcomp_e(double k, gsl_mode_t mode,
 * gsl_sf_result * result); missing: 4   gsl_sf_erfc_e  -> /usr/local/include/gsl/gsl_sf_erf.h:int
 * gsl_sf_erfc_e(double x, gsl_sf_result * result); missing: 4   gsl_sf_log_erfc_e  ->
 * /usr/local/include/gsl/gsl_sf_erf.h:int gsl_sf_log_erfc_e(double x, gsl_sf_result * result);
 * missing: 4   gsl_sf_erf_e  -> /usr/local/include/gsl/gsl_sf_erf.h:int gsl_sf_erf_e(double x,
 * gsl_sf_result * result); missing: 4   gsl_sf_erf_Z_e  -> /usr/local/include/gsl/gsl_sf_erf.h:int
 * gsl_sf_erf_Z_e(double x, gsl_sf_result * result); missing: 4   gsl_sf_erf_Q_e  ->
 * /usr/local/include/gsl/gsl_sf_erf.h:int gsl_sf_erf_Q_e(double x, gsl_sf_result * result);
 * missing: 4   gsl_sf_hazard_e  -> /usr/local/include/gsl/gsl_sf_erf.h:int gsl_sf_hazard_e(double
 * x, gsl_sf_result * result); missing: 4   gsl_sf_exp_e  -> /usr/local/include/gsl/gsl_sf_exp.h:int
 * gsl_sf_exp_e(const double x, gsl_sf_result * result); missing: 4   gsl_sf_exp_e10_e  ->
 * /usr/local/include/gsl/gsl_sf_exp.h:int gsl_sf_exp_e10_e(const double x, gsl_sf_result_e10 *
 * result); missing: 5   gsl_sf_exp_mult_e  -> /usr/local/include/gsl/gsl_sf_exp.h:int
 * gsl_sf_exp_mult_e(const double x, const double y, gsl_sf_result * result); missing: 5
 * gsl_sf_exp_mult_e10_e  -> /usr/local/include/gsl/gsl_sf_exp.h:int gsl_sf_exp_mult_e10_e(const
 * double x, const double y, gsl_sf_result_e10 * result); missing: 4   gsl_sf_expm1_e  ->
 * /usr/local/include/gsl/gsl_sf_exp.h:int gsl_sf_expm1_e(const double x, gsl_sf_result * result);
 * missing: 4   gsl_sf_exprel_e  -> /usr/local/include/gsl/gsl_sf_exp.h:int gsl_sf_exprel_e(const
 * double x, gsl_sf_result * result); missing: 4   gsl_sf_exprel_2_e  ->
 * /usr/local/include/gsl/gsl_sf_exp.h:int gsl_sf_exprel_2_e(double x, gsl_sf_result * result);
 * missing: 5   gsl_sf_exprel_n_e  -> /usr/local/include/gsl/gsl_sf_exp.h:int
 * gsl_sf_exprel_n_e(const int n, const double x, gsl_sf_result * result); missing: 5
 * gsl_sf_exprel_n_CF_e  -> /usr/local/include/gsl/gsl_sf_exp.h:int gsl_sf_exprel_n_CF_e(const
 * double n, const double x, gsl_sf_result * result); missing: 5   gsl_sf_exp_err_e  ->
 * /usr/local/include/gsl/gsl_sf_exp.h:int gsl_sf_exp_err_e(const double x, const double dx,
 * gsl_sf_result * result); missing: 5   gsl_sf_exp_err_e10_e  ->
 * /usr/local/include/gsl/gsl_sf_exp.h:int gsl_sf_exp_err_e10_e(const double x, const double dx,
 * gsl_sf_result_e10 * result); missing: 7   gsl_sf_exp_err_e10_e  ->
 * /usr/local/include/gsl/gsl_sf_exp.h:int gsl_sf_exp_mult_err_e(const double x, const double dx,
 * const double y, const double dy, gsl_sf_result * result); missing: 7   gsl_sf_exp_err_e10_e  ->
 * /usr/local/include/gsl/gsl_sf_exp.h:int gsl_sf_exp_mult_err_e10_e(const double x, const double
 * dx, const double y, const double dy, gsl_sf_result_e10 * result); missing: 4   gsl_sf_expint_E1_e
 * -> /usr/local/include/gsl/gsl_sf_expint.h:int     gsl_sf_expint_E1_e(const double x,
 * gsl_sf_result * result); missing: 4   gsl_sf_expint_E2_e  ->
 * /usr/local/include/gsl/gsl_sf_expint.h:int     gsl_sf_expint_E2_e(const double x, gsl_sf_result *
 * result); missing: 5   gsl_sf_expint_En_e  -> /usr/local/include/gsl/gsl_sf_expint.h:int
 * gsl_sf_expint_En_e(const int n, const double x, gsl_sf_result * result); missing: 4
 * gsl_sf_expint_E1_scaled_e  -> /usr/local/include/gsl/gsl_sf_expint.h:int
 * gsl_sf_expint_E1_scaled_e(const double x, gsl_sf_result * result); missing: 4
 * gsl_sf_expint_E2_scaled_e  -> /usr/local/include/gsl/gsl_sf_expint.h:int
 * gsl_sf_expint_E2_scaled_e(const double x, gsl_sf_result * result); missing: 5
 * gsl_sf_expint_En_scaled_e  -> /usr/local/include/gsl/gsl_sf_expint.h:int
 * gsl_sf_expint_En_scaled_e(const int n, const double x, gsl_sf_result * result); missing: 4
 * gsl_sf_expint_Ei_e  -> /usr/local/include/gsl/gsl_sf_expint.h:int     gsl_sf_expint_Ei_e(const
 * double x, gsl_sf_result * result); missing: 4   gsl_sf_expint_Ei_scaled_e  ->
 * /usr/local/include/gsl/gsl_sf_expint.h:int     gsl_sf_expint_Ei_scaled_e(const double x,
 * gsl_sf_result * result); missing: 4   gsl_sf_Shi_e  -> /usr/local/include/gsl/gsl_sf_expint.h:int
 * gsl_sf_Shi_e(const double x, gsl_sf_result * result); missing: 4   gsl_sf_Chi_e  ->
 * /usr/local/include/gsl/gsl_sf_expint.h:int     gsl_sf_Chi_e(const double x, gsl_sf_result *
 * result); missing: 4   gsl_sf_expint_3_e  -> /usr/local/include/gsl/gsl_sf_expint.h:int
 * gsl_sf_expint_3_e(const double x, gsl_sf_result * result); missing: 4   gsl_sf_Si_e  ->
 * /usr/local/include/gsl/gsl_sf_expint.h:int     gsl_sf_Si_e(const double x, gsl_sf_result *
 * result); missing: 4   gsl_sf_Ci_e  -> /usr/local/include/gsl/gsl_sf_expint.h:int
 * gsl_sf_Ci_e(const double x, gsl_sf_result * result); missing: 4   gsl_sf_atanint_e  ->
 * /usr/local/include/gsl/gsl_sf_expint.h:int     gsl_sf_atanint_e(const double x, gsl_sf_result *
 * result); missing: 4   gsl_sf_fermi_dirac_m1_e  -> /usr/local/include/gsl/gsl_sf_fermi_dirac.h:int
 * gsl_sf_fermi_dirac_m1_e(const double x, gsl_sf_result * result); missing: 4
 * gsl_sf_fermi_dirac_0_e  -> /usr/local/include/gsl/gsl_sf_fermi_dirac.h:int
 * gsl_sf_fermi_dirac_0_e(const double x, gsl_sf_result * result); missing: 4 gsl_sf_fermi_dirac_1_e
 * -> /usr/local/include/gsl/gsl_sf_fermi_dirac.h:int     gsl_sf_fermi_dirac_1_e(const double x,
 * gsl_sf_result * result); missing: 4   gsl_sf_fermi_dirac_2_e  ->
 * /usr/local/include/gsl/gsl_sf_fermi_dirac.h:int     gsl_sf_fermi_dirac_2_e(const double x,
 * gsl_sf_result * result); missing: 5   gsl_sf_fermi_dirac_int_e  ->
 * /usr/local/include/gsl/gsl_sf_fermi_dirac.h:int     gsl_sf_fermi_dirac_int_e(const int j, const
 * double x, gsl_sf_result * result); missing: 4   gsl_sf_fermi_dirac_mhalf_e  ->
 * /usr/local/include/gsl/gsl_sf_fermi_dirac.h:int     gsl_sf_fermi_dirac_mhalf_e(const double x,
 * gsl_sf_result * result); missing: 4   gsl_sf_fermi_dirac_half_e  ->
 * /usr/local/include/gsl/gsl_sf_fermi_dirac.h:int     gsl_sf_fermi_dirac_half_e(const double x,
 * gsl_sf_result * result); missing: 4   gsl_sf_fermi_dirac_3half_e  ->
 * /usr/local/include/gsl/gsl_sf_fermi_dirac.h:int     gsl_sf_fermi_dirac_3half_e(const double x,
 * gsl_sf_result * result); missing: 5   gsl_sf_fermi_dirac_inc_0_e  ->
 * /usr/local/include/gsl/gsl_sf_fermi_dirac.h:int     gsl_sf_fermi_dirac_inc_0_e(const double x,
 * const double b, gsl_sf_result * result); missing: 4   gsl_sf_lngamma_e  ->
 * /usr/local/include/gsl/gsl_sf_gamma.h:int gsl_sf_lngamma_e(double x, gsl_sf_result * result);
 * missing: 5   gsl_sf_lngamma_sgn_e  -> /usr/local/include/gsl/gsl_sf_gamma.h:int
 * gsl_sf_lngamma_sgn_e(double x, gsl_sf_result * result_lg, double *sgn); missing: 4 gsl_sf_gamma_e
 * -> /usr/local/include/gsl/gsl_sf_gamma.h:int gsl_sf_gamma_e(const double x, gsl_sf_result *
 * result); missing: 4   gsl_sf_gammastar_e  -> /usr/local/include/gsl/gsl_sf_gamma.h:int
 * gsl_sf_gammastar_e(const double x, gsl_sf_result * result); missing: 4   gsl_sf_gammainv_e  ->
 * /usr/local/include/gsl/gsl_sf_gamma.h:int gsl_sf_gammainv_e(const double x, gsl_sf_result *
 * result); missing: 5   gsl_sf_taylorcoeff_e  -> /usr/local/include/gsl/gsl_sf_gamma.h:int
 * gsl_sf_taylorcoeff_e(const int n, const double x, gsl_sf_result * result); missing: 4
 * gsl_sf_fact_e  -> /usr/local/include/gsl/gsl_sf_gamma.h:int gsl_sf_fact_e(const unsigned int n,
 * gsl_sf_result * result); missing: 4   gsl_sf_doublefact_e  ->
 * /usr/local/include/gsl/gsl_sf_gamma.h:int gsl_sf_doublefact_e(const unsigned int n, gsl_sf_result
 * * result); missing: 4   gsl_sf_lnfact_e  -> /usr/local/include/gsl/gsl_sf_gamma.h:int
 * gsl_sf_lnfact_e(const unsigned int n, gsl_sf_result * result); missing: 4   gsl_sf_lndoublefact_e
 * -> /usr/local/include/gsl/gsl_sf_gamma.h:int gsl_sf_lndoublefact_e(const unsigned int n,
 * gsl_sf_result * result); missing: 5   gsl_sf_lnchoose_e  ->
 * /usr/local/include/gsl/gsl_sf_gamma.h:int gsl_sf_lnchoose_e(unsigned int n, unsigned int m,
 * gsl_sf_result * result); missing: 5   gsl_sf_choose_e  ->
 * /usr/local/include/gsl/gsl_sf_gamma.h:int gsl_sf_choose_e(unsigned int n, unsigned int m,
 * gsl_sf_result * result); missing: 5   gsl_sf_lnpoch_e  ->
 * /usr/local/include/gsl/gsl_sf_gamma.h:int gsl_sf_lnpoch_e(const double a, const double x,
 * gsl_sf_result * result); missing: 5   gsl_sf_poch_e  -> /usr/local/include/gsl/gsl_sf_gamma.h:int
 * gsl_sf_poch_e(const double a, const double x, gsl_sf_result * result); missing: 5
 * gsl_sf_pochrel_e  -> /usr/local/include/gsl/gsl_sf_gamma.h:int gsl_sf_pochrel_e(const double a,
 * const double x, gsl_sf_result * result); missing: 5   gsl_sf_gamma_inc_Q_e  ->
 * /usr/local/include/gsl/gsl_sf_gamma.h:int gsl_sf_gamma_inc_Q_e(const double a, const double x,
 * gsl_sf_result * result); missing: 5   gsl_sf_gamma_inc_P_e  ->
 * /usr/local/include/gsl/gsl_sf_gamma.h:int gsl_sf_gamma_inc_P_e(const double a, const double x,
 * gsl_sf_result * result); missing: 5   gsl_sf_gamma_inc_e  ->
 * /usr/local/include/gsl/gsl_sf_gamma.h:int gsl_sf_gamma_inc_e(const double a, const double x,
 * gsl_sf_result * result); missing: 5   gsl_sf_lnbeta_e  ->
 * /usr/local/include/gsl/gsl_sf_gamma.h:int gsl_sf_lnbeta_e(const double a, const double b,
 * gsl_sf_result * result); missing: 5   gsl_sf_beta_e  -> /usr/local/include/gsl/gsl_sf_gamma.h:int
 * gsl_sf_beta_e(const double a, const double b, gsl_sf_result * result); missing: 5
 * gsl_sf_gegenpoly_1_e  -> /usr/local/include/gsl/gsl_sf_gegenbauer.h:int
 * gsl_sf_gegenpoly_1_e(double lambda, double x, gsl_sf_result * result); missing: 5
 * gsl_sf_gegenpoly_2_e  -> /usr/local/include/gsl/gsl_sf_gegenbauer.h:int
 * gsl_sf_gegenpoly_2_e(double lambda, double x, gsl_sf_result * result); missing: 5
 * gsl_sf_gegenpoly_3_e  -> /usr/local/include/gsl/gsl_sf_gegenbauer.h:int
 * gsl_sf_gegenpoly_3_e(double lambda, double x, gsl_sf_result * result); missing: 5
 * gsl_sf_hyperg_0F1_e  -> /usr/local/include/gsl/gsl_sf_hyperg.h:int gsl_sf_hyperg_0F1_e(double c,
 * double x, gsl_sf_result * result); missing: 5   gsl_sf_laguerre_1_e  ->
 * /usr/local/include/gsl/gsl_sf_laguerre.h:int gsl_sf_laguerre_1_e(const double a, const double x,
 * gsl_sf_result * result); missing: 5   gsl_sf_laguerre_2_e  ->
 * /usr/local/include/gsl/gsl_sf_laguerre.h:int gsl_sf_laguerre_2_e(const double a, const double x,
 * gsl_sf_result * result); missing: 5   gsl_sf_laguerre_3_e  ->
 * /usr/local/include/gsl/gsl_sf_laguerre.h:int gsl_sf_laguerre_3_e(const double a, const double x,
 * gsl_sf_result * result); missing: 4   gsl_sf_lambert_W0_e  ->
 * /usr/local/include/gsl/gsl_sf_lambert.h:int     gsl_sf_lambert_W0_e(double x, gsl_sf_result *
 * result); missing: 4   gsl_sf_lambert_Wm1_e  -> /usr/local/include/gsl/gsl_sf_lambert.h:int
 * gsl_sf_lambert_Wm1_e(double x, gsl_sf_result * result); missing: 5   gsl_sf_legendre_Pl_e  ->
 * /usr/local/include/gsl/gsl_sf_legendre.h:int     gsl_sf_legendre_Pl_e(const int l, const double
 * x, gsl_sf_result * result); missing: 4   gsl_sf_legendre_P1_e  ->
 * /usr/local/include/gsl/gsl_sf_legendre.h:int gsl_sf_legendre_P1_e(double x, gsl_sf_result *
 * result); missing: 4   gsl_sf_legendre_P2_e  -> /usr/local/include/gsl/gsl_sf_legendre.h:int
 * gsl_sf_legendre_P2_e(double x, gsl_sf_result * result); missing: 4   gsl_sf_legendre_P3_e  ->
 * /usr/local/include/gsl/gsl_sf_legendre.h:int gsl_sf_legendre_P3_e(double x, gsl_sf_result *
 * result); missing: 4   gsl_sf_legendre_Q0_e  -> /usr/local/include/gsl/gsl_sf_legendre.h:int
 * gsl_sf_legendre_Q0_e(const double x, gsl_sf_result * result); missing: 4   gsl_sf_legendre_Q1_e
 * -> /usr/local/include/gsl/gsl_sf_legendre.h:int gsl_sf_legendre_Q1_e(const double x,
 * gsl_sf_result * result); missing: 5   gsl_sf_legendre_Ql_e  ->
 * /usr/local/include/gsl/gsl_sf_legendre.h:int gsl_sf_legendre_Ql_e(const int l, const double x,
 * gsl_sf_result * result); missing: 5   gsl_sf_conicalP_half_e  ->
 * /usr/local/include/gsl/gsl_sf_legendre.h:int gsl_sf_conicalP_half_e(const double lambda, const
 * double x, gsl_sf_result * result); missing: 5   gsl_sf_conicalP_mhalf_e  ->
 * /usr/local/include/gsl/gsl_sf_legendre.h:int gsl_sf_conicalP_mhalf_e(const double lambda, const
 * double x, gsl_sf_result * result); missing: 5   gsl_sf_conicalP_0_e  ->
 * /usr/local/include/gsl/gsl_sf_legendre.h:int gsl_sf_conicalP_0_e(const double lambda, const
 * double x, gsl_sf_result * result); missing: 5   gsl_sf_conicalP_1_e  ->
 * /usr/local/include/gsl/gsl_sf_legendre.h:int gsl_sf_conicalP_1_e(const double lambda, const
 * double x, gsl_sf_result * result); missing: 5   gsl_sf_legendre_H3d_0_e  ->
 * /usr/local/include/gsl/gsl_sf_legendre.h:int gsl_sf_legendre_H3d_0_e(const double lambda, const
 * double eta, gsl_sf_result * result); missing: 5   gsl_sf_legendre_H3d_1_e  ->
 * /usr/local/include/gsl/gsl_sf_legendre.h:int gsl_sf_legendre_H3d_1_e(const double lambda, const
 * double eta, gsl_sf_result * result); missing: 4   gsl_sf_log_e  ->
 * /usr/local/include/gsl/gsl_sf_log.h:int gsl_sf_log_e(const double x, gsl_sf_result * result);
 * missing: 4   gsl_sf_log_abs_e  -> /usr/local/include/gsl/gsl_sf_log.h:int gsl_sf_log_abs_e(const
 * double x, gsl_sf_result * result); missing: 4   gsl_sf_log_1plusx_e  ->
 * /usr/local/include/gsl/gsl_sf_log.h:int gsl_sf_log_1plusx_e(const double x, gsl_sf_result *
 * result); missing: 4   gsl_sf_log_1plusx_mx_e  -> /usr/local/include/gsl/gsl_sf_log.h:int
 * gsl_sf_log_1plusx_mx_e(const double x, gsl_sf_result * result);
 * //  -- missing type "gsl_sf_result"
 * missing: 5   gsl_sf_mathieu_a  -> /usr/local/include/gsl/gsl_sf_mathieu.h:int
 * gsl_sf_mathieu_a(int order, double qq, gsl_sf_result *result); missing: 5   gsl_sf_mathieu_b  ->
 * /usr/local/include/gsl/gsl_sf_mathieu.h:int gsl_sf_mathieu_b(int order, double qq, gsl_sf_result
 * *result); missing: 6   gsl_sf_mathieu_b  -> /usr/local/include/gsl/gsl_sf_mathieu.h:int
 * gsl_sf_mathieu_a_coeff(int order, double qq, double aa, double coeff[]); missing: 6
 * gsl_sf_mathieu_b  -> /usr/local/include/gsl/gsl_sf_mathieu.h:int gsl_sf_mathieu_b_coeff(int
 * order, double qq, double aa, double coeff[]); missing: 6   gsl_sf_mathieu_b  ->
 * /usr/local/include/gsl/gsl_sf_mathieu.h:int gsl_sf_mathieu_ce(int order, double qq, double zz,
 * gsl_sf_result *result); missing: 6   gsl_sf_mathieu_b  ->
 * /usr/local/include/gsl/gsl_sf_mathieu.h:int gsl_sf_mathieu_se(int order, double qq, double zz,
 * gsl_sf_result *result); missing: 5   gsl_sf_pow_int_e  ->
 * /usr/local/include/gsl/gsl_sf_pow_int.h:int     gsl_sf_pow_int_e(double x, int n, gsl_sf_result *
 * result); missing: 4   gsl_sf_psi_int_e  -> /usr/local/include/gsl/gsl_sf_psi.h:int
 * gsl_sf_psi_int_e(const int n, gsl_sf_result * result); missing: 4   gsl_sf_psi_e  ->
 * /usr/local/include/gsl/gsl_sf_psi.h:int     gsl_sf_psi_e(const double x, gsl_sf_result * result);
 * missing: 4   gsl_sf_psi_1piy_e  -> /usr/local/include/gsl/gsl_sf_psi.h:int
 * gsl_sf_psi_1piy_e(const double y, gsl_sf_result * result); missing: 4   gsl_sf_psi_1_int_e  ->
 * /usr/local/include/gsl/gsl_sf_psi.h:int     gsl_sf_psi_1_int_e(const int n, gsl_sf_result *
 * result); missing: 4   gsl_sf_psi_1_e  -> /usr/local/include/gsl/gsl_sf_psi.h:int
 * gsl_sf_psi_1_e(const double x, gsl_sf_result * result); missing: 5   gsl_sf_psi_n_e  ->
 * /usr/local/include/gsl/gsl_sf_psi.h:int     gsl_sf_psi_n_e(const int n, const double x,
 * gsl_sf_result * result);
 * //  -- missing type "const gsl_sf_result_e10 *"
 * missing: 4   gsl_sf_result_smash_e  -> /usr/local/include/gsl/gsl_sf_result.h:int
 * gsl_sf_result_smash_e(const gsl_sf_result_e10 * re, gsl_sf_result * r); missing: 4
 * gsl_sf_synchrotron_1_e  -> /usr/local/include/gsl/gsl_sf_synchrotron.h:int
 * gsl_sf_synchrotron_1_e(const double x, gsl_sf_result * result); missing: 4 gsl_sf_synchrotron_2_e
 * -> /usr/local/include/gsl/gsl_sf_synchrotron.h:int     gsl_sf_synchrotron_2_e(const double x,
 * gsl_sf_result * result); missing: 4   gsl_sf_transport_2_e  ->
 * /usr/local/include/gsl/gsl_sf_transport.h:int     gsl_sf_transport_2_e(const double x,
 * gsl_sf_result * result); missing: 4   gsl_sf_transport_3_e  ->
 * /usr/local/include/gsl/gsl_sf_transport.h:int     gsl_sf_transport_3_e(const double x,
 * gsl_sf_result * result); missing: 4   gsl_sf_transport_4_e  ->
 * /usr/local/include/gsl/gsl_sf_transport.h:int     gsl_sf_transport_4_e(const double x,
 * gsl_sf_result * result); missing: 4   gsl_sf_transport_5_e  ->
 * /usr/local/include/gsl/gsl_sf_transport.h:int     gsl_sf_transport_5_e(const double x,
 * gsl_sf_result * result); missing: 4   gsl_sf_sin_e  -> /usr/local/include/gsl/gsl_sf_trig.h:int
 * gsl_sf_sin_e(double x, gsl_sf_result * result); missing: 4   gsl_sf_cos_e  ->
 * /usr/local/include/gsl/gsl_sf_trig.h:int gsl_sf_cos_e(double x, gsl_sf_result * result); missing:
 * 5   gsl_sf_hypot_e  -> /usr/local/include/gsl/gsl_sf_trig.h:int gsl_sf_hypot_e(const double x,
 * const double y, gsl_sf_result * result); missing: 4   gsl_sf_sinc_e  ->
 * /usr/local/include/gsl/gsl_sf_trig.h:int gsl_sf_sinc_e(double x, gsl_sf_result * result);
 * missing: 4   gsl_sf_lnsinh_e  -> /usr/local/include/gsl/gsl_sf_trig.h:int gsl_sf_lnsinh_e(const
 * double x, gsl_sf_result * result); missing: 4   gsl_sf_lncosh_e  ->
 * /usr/local/include/gsl/gsl_sf_trig.h:int gsl_sf_lncosh_e(const double x, gsl_sf_result * result);
 * missing: 5   gsl_sf_sin_err_e  -> /usr/local/include/gsl/gsl_sf_trig.h:int gsl_sf_sin_err_e(const
 * double x, const double dx, gsl_sf_result * result); missing: 5   gsl_sf_cos_err_e  ->
 * /usr/local/include/gsl/gsl_sf_trig.h:int gsl_sf_cos_err_e(const double x, const double dx,
 * gsl_sf_result * result); missing: 3   gsl_sf_angle_restrict_symm_e  ->
 * /usr/local/include/gsl/gsl_sf_trig.h:int gsl_sf_angle_restrict_symm_e(double * theta); missing: 3
 * gsl_sf_angle_restrict_pos_e  -> /usr/local/include/gsl/gsl_sf_trig.h:int
 * gsl_sf_angle_restrict_pos_e(double * theta); missing: 4   gsl_sf_angle_restrict_symm_err_e  ->
 * /usr/local/include/gsl/gsl_sf_trig.h:int gsl_sf_angle_restrict_symm_err_e(const double theta,
 * gsl_sf_result * result); missing: 4   gsl_sf_angle_restrict_pos_err_e  ->
 * /usr/local/include/gsl/gsl_sf_trig.h:int gsl_sf_angle_restrict_pos_err_e(const double theta,
 * gsl_sf_result * result); missing: 4   gsl_sf_zeta_int_e  ->
 * /usr/local/include/gsl/gsl_sf_zeta.h:int gsl_sf_zeta_int_e(const int n, gsl_sf_result * result);
 * missing: 4   gsl_sf_zeta_e  -> /usr/local/include/gsl/gsl_sf_zeta.h:int gsl_sf_zeta_e(const
 * double s, gsl_sf_result * result); missing: 4   gsl_sf_zetam1_e  ->
 * /usr/local/include/gsl/gsl_sf_zeta.h:int gsl_sf_zetam1_e(const double s, gsl_sf_result * result);
 * missing: 4   gsl_sf_zetam1_int_e  -> /usr/local/include/gsl/gsl_sf_zeta.h:int
 * gsl_sf_zetam1_int_e(const int s, gsl_sf_result * result); missing: 5   gsl_sf_hzeta_e  ->
 * /usr/local/include/gsl/gsl_sf_zeta.h:int gsl_sf_hzeta_e(const double s, const double q,
 * gsl_sf_result * result); missing: 4   gsl_sf_eta_int_e  ->
 * /usr/local/include/gsl/gsl_sf_zeta.h:int gsl_sf_eta_int_e(int n, gsl_sf_result * result);
 * missing: 4   gsl_sf_eta_e  -> /usr/local/include/gsl/gsl_sf_zeta.h:int gsl_sf_eta_e(const double
 * s, gsl_sf_result * result);
 */
/*****************/
/*****************/
double gsl_cdf_ugaussian_P__(double const &x) { return gsl_cdf_ugaussian_P((const double)x); }

double gsl_cdf_ugaussian_Q__(double const &x) { return gsl_cdf_ugaussian_Q((const double)x); }

double gsl_cdf_ugaussian_Pinv__(double const &x) { return gsl_cdf_ugaussian_Pinv((const double)x); }

double gsl_cdf_ugaussian_Qinv__(double const &x) { return gsl_cdf_ugaussian_Qinv((const double)x); }

double gsl_cdf_gaussian_P__(double const &x, double const &y) {
  return gsl_cdf_gaussian_P((const double)x, (const double)y);
}

double gsl_cdf_gaussian_Q__(double const &x, double const &y) {
  return gsl_cdf_gaussian_Q((const double)x, (const double)y);
}

double gsl_cdf_gaussian_Pinv__(double const &x, double const &y) {
  return gsl_cdf_gaussian_Pinv((const double)x, (const double)y);
}

double gsl_cdf_gaussian_Qinv__(double const &x, double const &y) {
  return gsl_cdf_gaussian_Qinv((const double)x, (const double)y);
}

double gsl_cdf_gamma_P__(double const &x, double const &y, double const &z) {
  return gsl_cdf_gamma_P((const double)x, (const double)y, (const double)z);
}

double gsl_cdf_gamma_Q__(double const &x, double const &y, double const &z) {
  return gsl_cdf_gamma_Q((const double)x, (const double)y, (const double)z);
}

double gsl_cdf_gamma_Pinv__(double const &x, double const &y, double const &z) {
  return gsl_cdf_gamma_Pinv((const double)x, (const double)y, (const double)z);
}

double gsl_cdf_gamma_Qinv__(double const &x, double const &y, double const &z) {
  return gsl_cdf_gamma_Qinv((const double)x, (const double)y, (const double)z);
}

double gsl_cdf_cauchy_P__(double const &x, double const &y) {
  return gsl_cdf_cauchy_P((const double)x, (const double)y);
}

double gsl_cdf_cauchy_Q__(double const &x, double const &y) {
  return gsl_cdf_cauchy_Q((const double)x, (const double)y);
}

double gsl_cdf_cauchy_Pinv__(double const &x, double const &y) {
  return gsl_cdf_cauchy_Pinv((const double)x, (const double)y);
}

double gsl_cdf_cauchy_Qinv__(double const &x, double const &y) {
  return gsl_cdf_cauchy_Qinv((const double)x, (const double)y);
}

double gsl_cdf_laplace_P__(double const &x, double const &y) {
  return gsl_cdf_laplace_P((const double)x, (const double)y);
}

double gsl_cdf_laplace_Q__(double const &x, double const &y) {
  return gsl_cdf_laplace_Q((const double)x, (const double)y);
}

double gsl_cdf_laplace_Pinv__(double const &x, double const &y) {
  return gsl_cdf_laplace_Pinv((const double)x, (const double)y);
}

double gsl_cdf_laplace_Qinv__(double const &x, double const &y) {
  return gsl_cdf_laplace_Qinv((const double)x, (const double)y);
}

double gsl_cdf_rayleigh_P__(double const &x, double const &y) {
  return gsl_cdf_rayleigh_P((const double)x, (const double)y);
}

double gsl_cdf_rayleigh_Q__(double const &x, double const &y) {
  return gsl_cdf_rayleigh_Q((const double)x, (const double)y);
}

double gsl_cdf_rayleigh_Pinv__(double const &x, double const &y) {
  return gsl_cdf_rayleigh_Pinv((const double)x, (const double)y);
}

double gsl_cdf_rayleigh_Qinv__(double const &x, double const &y) {
  return gsl_cdf_rayleigh_Qinv((const double)x, (const double)y);
}

double gsl_cdf_chisq_P__(double const &x, double const &y) {
  return gsl_cdf_chisq_P((const double)x, (const double)y);
}

double gsl_cdf_chisq_Q__(double const &x, double const &y) {
  return gsl_cdf_chisq_Q((const double)x, (const double)y);
}

double gsl_cdf_chisq_Pinv__(double const &x, double const &y) {
  return gsl_cdf_chisq_Pinv((const double)x, (const double)y);
}

double gsl_cdf_chisq_Qinv__(double const &x, double const &y) {
  return gsl_cdf_chisq_Qinv((const double)x, (const double)y);
}

double gsl_cdf_exponential_P__(double const &x, double const &y) {
  return gsl_cdf_exponential_P((const double)x, (const double)y);
}

double gsl_cdf_exponential_Q__(double const &x, double const &y) {
  return gsl_cdf_exponential_Q((const double)x, (const double)y);
}

double gsl_cdf_exponential_Pinv__(double const &x, double const &y) {
  return gsl_cdf_exponential_Pinv((const double)x, (const double)y);
}

double gsl_cdf_exponential_Qinv__(double const &x, double const &y) {
  return gsl_cdf_exponential_Qinv((const double)x, (const double)y);
}

double gsl_cdf_exppow_P__(double const &x, double const &y, double const &z) {
  return gsl_cdf_exppow_P((const double)x, (const double)y, (const double)z);
}

double gsl_cdf_exppow_Q__(double const &x, double const &y, double const &z) {
  return gsl_cdf_exppow_Q((const double)x, (const double)y, (const double)z);
}

double gsl_cdf_tdist_P__(double const &x, double const &y) {
  return gsl_cdf_tdist_P((const double)x, (const double)y);
}

double gsl_cdf_tdist_Q__(double const &x, double const &y) {
  return gsl_cdf_tdist_Q((const double)x, (const double)y);
}

double gsl_cdf_tdist_Pinv__(double const &x, double const &y) {
  return gsl_cdf_tdist_Pinv((const double)x, (const double)y);
}

double gsl_cdf_tdist_Qinv__(double const &x, double const &y) {
  return gsl_cdf_tdist_Qinv((const double)x, (const double)y);
}

double gsl_cdf_fdist_P__(double const &x, double const &y, double const &z) {
  return gsl_cdf_fdist_P((const double)x, (const double)y, (const double)z);
}

double gsl_cdf_fdist_Q__(double const &x, double const &y, double const &z) {
  return gsl_cdf_fdist_Q((const double)x, (const double)y, (const double)z);
}

double gsl_cdf_fdist_Pinv__(double const &x, double const &y, double const &z) {
  return gsl_cdf_fdist_Pinv((const double)x, (const double)y, (const double)z);
}

double gsl_cdf_fdist_Qinv__(double const &x, double const &y, double const &z) {
  return gsl_cdf_fdist_Qinv((const double)x, (const double)y, (const double)z);
}

double gsl_cdf_beta_P__(double const &x, double const &y, double const &z) {
  return gsl_cdf_beta_P((const double)x, (const double)y, (const double)z);
}

double gsl_cdf_beta_Q__(double const &x, double const &y, double const &z) {
  return gsl_cdf_beta_Q((const double)x, (const double)y, (const double)z);
}

double gsl_cdf_beta_Pinv__(double const &x, double const &y, double const &z) {
  return gsl_cdf_beta_Pinv((const double)x, (const double)y, (const double)z);
}

double gsl_cdf_beta_Qinv__(double const &x, double const &y, double const &z) {
  return gsl_cdf_beta_Qinv((const double)x, (const double)y, (const double)z);
}

double gsl_cdf_flat_P__(double const &x, double const &y, double const &z) {
  return gsl_cdf_flat_P((const double)x, (const double)y, (const double)z);
}

double gsl_cdf_flat_Q__(double const &x, double const &y, double const &z) {
  return gsl_cdf_flat_Q((const double)x, (const double)y, (const double)z);
}

double gsl_cdf_flat_Pinv__(double const &x, double const &y, double const &z) {
  return gsl_cdf_flat_Pinv((const double)x, (const double)y, (const double)z);
}

double gsl_cdf_flat_Qinv__(double const &x, double const &y, double const &z) {
  return gsl_cdf_flat_Qinv((const double)x, (const double)y, (const double)z);
}

double gsl_cdf_lognormal_P__(double const &x, double const &y, double const &z) {
  return gsl_cdf_lognormal_P((const double)x, (const double)y, (const double)z);
}

double gsl_cdf_lognormal_Q__(double const &x, double const &y, double const &z) {
  return gsl_cdf_lognormal_Q((const double)x, (const double)y, (const double)z);
}

double gsl_cdf_lognormal_Pinv__(double const &x, double const &y, double const &z) {
  return gsl_cdf_lognormal_Pinv((const double)x, (const double)y, (const double)z);
}

double gsl_cdf_lognormal_Qinv__(double const &x, double const &y, double const &z) {
  return gsl_cdf_lognormal_Qinv((const double)x, (const double)y, (const double)z);
}

double gsl_cdf_gumbel1_P__(double const &x, double const &y, double const &z) {
  return gsl_cdf_gumbel1_P((const double)x, (const double)y, (const double)z);
}

double gsl_cdf_gumbel1_Q__(double const &x, double const &y, double const &z) {
  return gsl_cdf_gumbel1_Q((const double)x, (const double)y, (const double)z);
}

double gsl_cdf_gumbel1_Pinv__(double const &x, double const &y, double const &z) {
  return gsl_cdf_gumbel1_Pinv((const double)x, (const double)y, (const double)z);
}

double gsl_cdf_gumbel1_Qinv__(double const &x, double const &y, double const &z) {
  return gsl_cdf_gumbel1_Qinv((const double)x, (const double)y, (const double)z);
}

double gsl_cdf_gumbel2_P__(double const &x, double const &y, double const &z) {
  return gsl_cdf_gumbel2_P((const double)x, (const double)y, (const double)z);
}

double gsl_cdf_gumbel2_Q__(double const &x, double const &y, double const &z) {
  return gsl_cdf_gumbel2_Q((const double)x, (const double)y, (const double)z);
}

double gsl_cdf_gumbel2_Pinv__(double const &x, double const &y, double const &z) {
  return gsl_cdf_gumbel2_Pinv((const double)x, (const double)y, (const double)z);
}

double gsl_cdf_gumbel2_Qinv__(double const &x, double const &y, double const &z) {
  return gsl_cdf_gumbel2_Qinv((const double)x, (const double)y, (const double)z);
}

double gsl_cdf_weibull_P__(double const &x, double const &y, double const &z) {
  return gsl_cdf_weibull_P((const double)x, (const double)y, (const double)z);
}

double gsl_cdf_weibull_Q__(double const &x, double const &y, double const &z) {
  return gsl_cdf_weibull_Q((const double)x, (const double)y, (const double)z);
}

double gsl_cdf_weibull_Pinv__(double const &x, double const &y, double const &z) {
  return gsl_cdf_weibull_Pinv((const double)x, (const double)y, (const double)z);
}

double gsl_cdf_weibull_Qinv__(double const &x, double const &y, double const &z) {
  return gsl_cdf_weibull_Qinv((const double)x, (const double)y, (const double)z);
}

double gsl_cdf_pareto_P__(double const &x, double const &y, double const &z) {
  return gsl_cdf_pareto_P((const double)x, (const double)y, (const double)z);
}

double gsl_cdf_pareto_Q__(double const &x, double const &y, double const &z) {
  return gsl_cdf_pareto_Q((const double)x, (const double)y, (const double)z);
}

double gsl_cdf_pareto_Pinv__(double const &x, double const &y, double const &z) {
  return gsl_cdf_pareto_Pinv((const double)x, (const double)y, (const double)z);
}

double gsl_cdf_pareto_Qinv__(double const &x, double const &y, double const &z) {
  return gsl_cdf_pareto_Qinv((const double)x, (const double)y, (const double)z);
}

double gsl_cdf_logistic_P__(double const &x, double const &y) {
  return gsl_cdf_logistic_P((const double)x, (const double)y);
}

double gsl_cdf_logistic_Q__(double const &x, double const &y) {
  return gsl_cdf_logistic_Q((const double)x, (const double)y);
}

double gsl_cdf_logistic_Pinv__(double const &x, double const &y) {
  return gsl_cdf_logistic_Pinv((const double)x, (const double)y);
}

double gsl_cdf_logistic_Qinv__(double const &x, double const &y) {
  return gsl_cdf_logistic_Qinv((const double)x, (const double)y);
}

double gsl_cdf_binomial_P__(long const &x, double const &y, long const &z) {
  return gsl_cdf_binomial_P((const unsigned int)x, (const double)y, (const unsigned int)z);
}

double gsl_cdf_binomial_Q__(long const &x, double const &y, long const &z) {
  return gsl_cdf_binomial_Q((const unsigned int)x, (const double)y, (const unsigned int)z);
}

double gsl_cdf_poisson_P__(long const &x, double const &y) {
  return gsl_cdf_poisson_P((const unsigned int)x, (const double)y);
}

double gsl_cdf_poisson_Q__(long const &x, double const &y) {
  return gsl_cdf_poisson_Q((const unsigned int)x, (const double)y);
}

double gsl_cdf_geometric_P__(long const &x, double const &y) {
  return gsl_cdf_geometric_P((const unsigned int)x, (const double)y);
}

double gsl_cdf_geometric_Q__(long const &x, double const &y) {
  return gsl_cdf_geometric_Q((const unsigned int)x, (const double)y);
}

double gsl_cdf_negative_binomial_P__(long const &x, double const &y, double const &z) {
  return gsl_cdf_negative_binomial_P((const unsigned int)x, (const double)y, (const double)z);
}

double gsl_cdf_negative_binomial_Q__(long const &x, double const &y, double const &z) {
  return gsl_cdf_negative_binomial_Q((const unsigned int)x, (const double)y, (const double)z);
}

double gsl_cdf_pascal_P__(long const &x, double const &y, long const &z) {
  return gsl_cdf_pascal_P((const unsigned int)x, (const double)y, (const unsigned int)z);
}

double gsl_cdf_pascal_Q__(long const &x, double const &y, long const &z) {
  return gsl_cdf_pascal_Q((const unsigned int)x, (const double)y, (const unsigned int)z);
}

double gsl_ran_bernoulli_pdf__(long const &x, double const &y) {
  return gsl_ran_bernoulli_pdf((const unsigned int)x, (double)y);
}

double gsl_ran_beta__(gsl_rng **const &x, double const &y, double const &z) {
  return gsl_ran_beta((const gsl_rng *)*x, (const double)y, (const double)z);
}

double gsl_ran_beta_pdf__(double const &x, double const &y, double const &z) {
  return gsl_ran_beta_pdf((const double)x, (const double)y, (const double)z);
}

double gsl_ran_binomial_pdf__(long const &x, double const &y, long const &z) {
  return gsl_ran_binomial_pdf((const unsigned int)x, (const double)y, (const unsigned int)z);
}

double gsl_ran_exponential__(gsl_rng **const &x, double const &y) {
  return gsl_ran_exponential((const gsl_rng *)*x, (const double)y);
}

double gsl_ran_exponential_pdf__(double const &x, double const &y) {
  return gsl_ran_exponential_pdf((const double)x, (const double)y);
}

double gsl_ran_exppow__(gsl_rng **const &x, double const &y, double const &z) {
  return gsl_ran_exppow((const gsl_rng *)*x, (const double)y, (const double)z);
}

double gsl_ran_exppow_pdf__(double const &x, double const &y, double const &z) {
  return gsl_ran_exppow_pdf((const double)x, (const double)y, (const double)z);
}

double gsl_ran_cauchy__(gsl_rng **const &x, double const &y) {
  return gsl_ran_cauchy((const gsl_rng *)*x, (const double)y);
}

double gsl_ran_cauchy_pdf__(double const &x, double const &y) {
  return gsl_ran_cauchy_pdf((const double)x, (const double)y);
}

double gsl_ran_chisq__(gsl_rng **const &x, double const &y) {
  return gsl_ran_chisq((const gsl_rng *)*x, (const double)y);
}

double gsl_ran_chisq_pdf__(double const &x, double const &y) {
  return gsl_ran_chisq_pdf((const double)x, (const double)y);
}

double gsl_ran_erlang__(gsl_rng **const &x, double const &y, double const &z) {
  return gsl_ran_erlang((const gsl_rng *)*x, (const double)y, (const double)z);
}

double gsl_ran_erlang_pdf__(double const &x, double const &y, double const &z) {
  return gsl_ran_erlang_pdf((const double)x, (const double)y, (const double)z);
}

double gsl_ran_fdist__(gsl_rng **const &x, double const &y, double const &z) {
  return gsl_ran_fdist((const gsl_rng *)*x, (const double)y, (const double)z);
}

double gsl_ran_fdist_pdf__(double const &x, double const &y, double const &z) {
  return gsl_ran_fdist_pdf((const double)x, (const double)y, (const double)z);
}

double gsl_ran_flat__(gsl_rng **const &x, double const &y, double const &z) {
  return gsl_ran_flat((const gsl_rng *)*x, (const double)y, (const double)z);
}

double gsl_ran_flat_pdf__(double const &x, double const &y, double const &z) {
  return gsl_ran_flat_pdf((double)x, (const double)y, (const double)z);
}

double gsl_ran_gamma__(gsl_rng **const &x, double const &y, double const &z) {
  return gsl_ran_gamma((const gsl_rng *)*x, (const double)y, (const double)z);
}

double gsl_ran_gamma_int__(gsl_rng **const &x, long const &y) {
  return gsl_ran_gamma_int((const gsl_rng *)*x, (const unsigned int)y);
}

double gsl_ran_gamma_pdf__(double const &x, double const &y, double const &z) {
  return gsl_ran_gamma_pdf((const double)x, (const double)y, (const double)z);
}

double gsl_ran_gamma_mt__(gsl_rng **const &x, double const &y, double const &z) {
  return gsl_ran_gamma_mt((const gsl_rng *)*x, (const double)y, (const double)z);
}

double gsl_ran_gamma_knuth__(gsl_rng **const &x, double const &y, double const &z) {
  return gsl_ran_gamma_knuth((const gsl_rng *)*x, (const double)y, (const double)z);
}

double gsl_ran_gaussian__(gsl_rng **const &x, double const &y) {
  return gsl_ran_gaussian((const gsl_rng *)*x, (const double)y);
}

double gsl_ran_gaussian_ratio_method__(gsl_rng **const &x, double const &y) {
  return gsl_ran_gaussian_ratio_method((const gsl_rng *)*x, (const double)y);
}

double gsl_ran_gaussian_ziggurat__(gsl_rng **const &x, double const &y) {
  return gsl_ran_gaussian_ziggurat((const gsl_rng *)*x, (const double)y);
}

double gsl_ran_gaussian_pdf__(double const &x, double const &y) {
  return gsl_ran_gaussian_pdf((const double)x, (const double)y);
}

double gsl_ran_ugaussian__(gsl_rng **const &x) { return gsl_ran_ugaussian((const gsl_rng *)*x); }

double gsl_ran_ugaussian_ratio_method__(gsl_rng **const &x) {
  return gsl_ran_ugaussian_ratio_method((const gsl_rng *)*x);
}

double gsl_ran_ugaussian_pdf__(double const &x) { return gsl_ran_ugaussian_pdf((const double)x); }

double gsl_ran_gaussian_tail__(gsl_rng **const &x, double const &y, double const &z) {
  return gsl_ran_gaussian_tail((const gsl_rng *)*x, (const double)y, (const double)z);
}

double gsl_ran_gaussian_tail_pdf__(double const &x, double const &y, double const &z) {
  return gsl_ran_gaussian_tail_pdf((const double)x, (const double)y, (const double)z);
}

double gsl_ran_ugaussian_tail__(gsl_rng **const &x, double const &y) {
  return gsl_ran_ugaussian_tail((const gsl_rng *)*x, (const double)y);
}

double gsl_ran_ugaussian_tail_pdf__(double const &x, double const &y) {
  return gsl_ran_ugaussian_tail_pdf((const double)x, (const double)y);
}

double gsl_ran_landau__(gsl_rng **const &x) { return gsl_ran_landau((const gsl_rng *)*x); }

double gsl_ran_landau_pdf__(double const &x) { return gsl_ran_landau_pdf((const double)x); }

double gsl_ran_geometric_pdf__(long const &x, double const &y) {
  return gsl_ran_geometric_pdf((const unsigned int)x, (const double)y);
}

double gsl_ran_gumbel1__(gsl_rng **const &x, double const &y, double const &z) {
  return gsl_ran_gumbel1((const gsl_rng *)*x, (const double)y, (const double)z);
}

double gsl_ran_gumbel1_pdf__(double const &x, double const &y, double const &z) {
  return gsl_ran_gumbel1_pdf((const double)x, (const double)y, (const double)z);
}

double gsl_ran_gumbel2__(gsl_rng **const &x, double const &y, double const &z) {
  return gsl_ran_gumbel2((const gsl_rng *)*x, (const double)y, (const double)z);
}

double gsl_ran_gumbel2_pdf__(double const &x, double const &y, double const &z) {
  return gsl_ran_gumbel2_pdf((const double)x, (const double)y, (const double)z);
}

double gsl_ran_logistic__(gsl_rng **const &x, double const &y) {
  return gsl_ran_logistic((const gsl_rng *)*x, (const double)y);
}

double gsl_ran_logistic_pdf__(double const &x, double const &y) {
  return gsl_ran_logistic_pdf((const double)x, (const double)y);
}

double gsl_ran_lognormal__(gsl_rng **const &x, double const &y, double const &z) {
  return gsl_ran_lognormal((const gsl_rng *)*x, (const double)y, (const double)z);
}

double gsl_ran_lognormal_pdf__(double const &x, double const &y, double const &z) {
  return gsl_ran_lognormal_pdf((const double)x, (const double)y, (const double)z);
}

double gsl_ran_logarithmic_pdf__(long const &x, double const &y) {
  return gsl_ran_logarithmic_pdf((const unsigned int)x, (const double)y);
}

double gsl_ran_negative_binomial_pdf__(long const &x, double const &y, double const &z) {
  return gsl_ran_negative_binomial_pdf((const unsigned int)x, (const double)y, (double)z);
}

double gsl_ran_pascal_pdf__(long const &x, double const &y, long const &z) {
  return gsl_ran_pascal_pdf((const unsigned int)x, (const double)y, (unsigned int)z);
}

double gsl_ran_pareto__(gsl_rng **const &x, double const &y, double const &z) {
  return gsl_ran_pareto((const gsl_rng *)*x, (double)y, (const double)z);
}

double gsl_ran_pareto_pdf__(double const &x, double const &y, double const &z) {
  return gsl_ran_pareto_pdf((const double)x, (const double)y, (const double)z);
}

double gsl_ran_poisson_pdf__(long const &x, double const &y) {
  return gsl_ran_poisson_pdf((const unsigned int)x, (const double)y);
}

double gsl_ran_rayleigh__(gsl_rng **const &x, double const &y) {
  return gsl_ran_rayleigh((const gsl_rng *)*x, (const double)y);
}

double gsl_ran_rayleigh_pdf__(double const &x, double const &y) {
  return gsl_ran_rayleigh_pdf((const double)x, (const double)y);
}

double gsl_ran_rayleigh_tail__(gsl_rng **const &x, double const &y, double const &z) {
  return gsl_ran_rayleigh_tail((const gsl_rng *)*x, (const double)y, (const double)z);
}

double gsl_ran_rayleigh_tail_pdf__(double const &x, double const &y, double const &z) {
  return gsl_ran_rayleigh_tail_pdf((const double)x, (const double)y, (const double)z);
}

double gsl_ran_tdist__(gsl_rng **const &x, double const &y) {
  return gsl_ran_tdist((const gsl_rng *)*x, (const double)y);
}

double gsl_ran_tdist_pdf__(double const &x, double const &y) {
  return gsl_ran_tdist_pdf((const double)x, (const double)y);
}

double gsl_ran_laplace__(gsl_rng **const &x, double const &y) {
  return gsl_ran_laplace((const gsl_rng *)*x, (const double)y);
}

double gsl_ran_laplace_pdf__(double const &x, double const &y) {
  return gsl_ran_laplace_pdf((const double)x, (const double)y);
}

double gsl_ran_levy__(gsl_rng **const &x, double const &y, double const &z) {
  return gsl_ran_levy((const gsl_rng *)*x, (const double)y, (const double)z);
}

double gsl_ran_weibull__(gsl_rng **const &x, double const &y, double const &z) {
  return gsl_ran_weibull((const gsl_rng *)*x, (const double)y, (const double)z);
}

double gsl_ran_weibull_pdf__(double const &x, double const &y, double const &z) {
  return gsl_ran_weibull_pdf((const double)x, (const double)y, (const double)z);
}

double gsl_sf_airy_Ai__(double const &x, long const &y) {
  return gsl_sf_airy_Ai((const double)x, (gsl_mode_t)y);
}

double gsl_sf_airy_Bi__(double const &x, long const &y) {
  return gsl_sf_airy_Bi((const double)x, (gsl_mode_t)y);
}

double gsl_sf_airy_Ai_scaled__(double const &x, long const &y) {
  return gsl_sf_airy_Ai_scaled((const double)x, (gsl_mode_t)y);
}

double gsl_sf_airy_Bi_scaled__(double const &x, long const &y) {
  return gsl_sf_airy_Bi_scaled((const double)x, (gsl_mode_t)y);
}

double gsl_sf_airy_Ai_deriv__(double const &x, long const &y) {
  return gsl_sf_airy_Ai_deriv((const double)x, (gsl_mode_t)y);
}

double gsl_sf_airy_Bi_deriv__(double const &x, long const &y) {
  return gsl_sf_airy_Bi_deriv((const double)x, (gsl_mode_t)y);
}

double gsl_sf_airy_Ai_deriv_scaled__(double const &x, long const &y) {
  return gsl_sf_airy_Ai_deriv_scaled((const double)x, (gsl_mode_t)y);
}

double gsl_sf_airy_Bi_deriv_scaled__(double const &x, long const &y) {
  return gsl_sf_airy_Bi_deriv_scaled((const double)x, (gsl_mode_t)y);
}

double gsl_sf_airy_zero_Ai__(long const &x) { return gsl_sf_airy_zero_Ai((unsigned int)x); }

double gsl_sf_airy_zero_Bi__(long const &x) { return gsl_sf_airy_zero_Bi((unsigned int)x); }

double gsl_sf_airy_zero_Ai_deriv__(long const &x) {
  return gsl_sf_airy_zero_Ai_deriv((unsigned int)x);
}

double gsl_sf_airy_zero_Bi_deriv__(long const &x) {
  return gsl_sf_airy_zero_Bi_deriv((unsigned int)x);
}

double gsl_sf_bessel_J0__(double const &x) { return gsl_sf_bessel_J0((const double)x); }

double gsl_sf_bessel_J1__(double const &x) { return gsl_sf_bessel_J1((const double)x); }

double gsl_sf_bessel_Jn__(long const &x, double const &y) {
  return gsl_sf_bessel_Jn((const int)x, (const double)y);
}

double gsl_sf_bessel_Y0__(double const &x) { return gsl_sf_bessel_Y0((const double)x); }

double gsl_sf_bessel_Y1__(double const &x) { return gsl_sf_bessel_Y1((const double)x); }

double gsl_sf_bessel_Yn__(long const &x, double const &y) {
  return gsl_sf_bessel_Yn((const int)x, (const double)y);
}

double gsl_sf_bessel_I0__(double const &x) { return gsl_sf_bessel_I0((const double)x); }

double gsl_sf_bessel_I1__(double const &x) { return gsl_sf_bessel_I1((const double)x); }

double gsl_sf_bessel_In__(long const &x, double const &y) {
  return gsl_sf_bessel_In((const int)x, (const double)y);
}

double gsl_sf_bessel_I0_scaled__(double const &x) {
  return gsl_sf_bessel_I0_scaled((const double)x);
}

double gsl_sf_bessel_I1_scaled__(double const &x) {
  return gsl_sf_bessel_I1_scaled((const double)x);
}

double gsl_sf_bessel_In_scaled__(long const &x, double const &y) {
  return gsl_sf_bessel_In_scaled((const int)x, (const double)y);
}

double gsl_sf_bessel_K0__(double const &x) { return gsl_sf_bessel_K0((const double)x); }

double gsl_sf_bessel_K1__(double const &x) { return gsl_sf_bessel_K1((const double)x); }

double gsl_sf_bessel_Kn__(long const &x, double const &y) {
  return gsl_sf_bessel_Kn((const int)x, (const double)y);
}

double gsl_sf_bessel_K0_scaled__(double const &x) {
  return gsl_sf_bessel_K0_scaled((const double)x);
}

double gsl_sf_bessel_K1_scaled__(double const &x) {
  return gsl_sf_bessel_K1_scaled((const double)x);
}

double gsl_sf_bessel_Kn_scaled__(long const &x, double const &y) {
  return gsl_sf_bessel_Kn_scaled((const int)x, (const double)y);
}

double gsl_sf_bessel_j0__(double const &x) { return gsl_sf_bessel_j0((const double)x); }

double gsl_sf_bessel_j1__(double const &x) { return gsl_sf_bessel_j1((const double)x); }

double gsl_sf_bessel_j2__(double const &x) { return gsl_sf_bessel_j2((const double)x); }

double gsl_sf_bessel_jl__(long const &x, double const &y) {
  return gsl_sf_bessel_jl((const int)x, (const double)y);
}

double gsl_sf_bessel_y0__(double const &x) { return gsl_sf_bessel_y0((const double)x); }

double gsl_sf_bessel_y1__(double const &x) { return gsl_sf_bessel_y1((const double)x); }

double gsl_sf_bessel_y2__(double const &x) { return gsl_sf_bessel_y2((const double)x); }

double gsl_sf_bessel_yl__(long const &x, double const &y) {
  return gsl_sf_bessel_yl((const int)x, (const double)y);
}

double gsl_sf_bessel_i0_scaled__(double const &x) {
  return gsl_sf_bessel_i0_scaled((const double)x);
}

double gsl_sf_bessel_i1_scaled__(double const &x) {
  return gsl_sf_bessel_i1_scaled((const double)x);
}

double gsl_sf_bessel_i2_scaled__(double const &x) {
  return gsl_sf_bessel_i2_scaled((const double)x);
}

double gsl_sf_bessel_il_scaled__(long const &x, double const &y) {
  return gsl_sf_bessel_il_scaled((const int)x, (const double)y);
}

double gsl_sf_bessel_k0_scaled__(double const &x) {
  return gsl_sf_bessel_k0_scaled((const double)x);
}

double gsl_sf_bessel_k1_scaled__(double const &x) {
  return gsl_sf_bessel_k1_scaled((const double)x);
}

double gsl_sf_bessel_k2_scaled__(double const &x) {
  return gsl_sf_bessel_k2_scaled((const double)x);
}

double gsl_sf_bessel_kl_scaled__(long const &x, double const &y) {
  return gsl_sf_bessel_kl_scaled((const int)x, (const double)y);
}

double gsl_sf_bessel_Jnu__(double const &x, double const &y) {
  return gsl_sf_bessel_Jnu((const double)x, (const double)y);
}

double gsl_sf_bessel_Ynu__(double const &x, double const &y) {
  return gsl_sf_bessel_Ynu((const double)x, (const double)y);
}

double gsl_sf_bessel_Inu_scaled__(double const &x, double const &y) {
  return gsl_sf_bessel_Inu_scaled((double)x, (double)y);
}

double gsl_sf_bessel_Inu__(double const &x, double const &y) {
  return gsl_sf_bessel_Inu((double)x, (double)y);
}

double gsl_sf_bessel_Knu_scaled__(double const &x, double const &y) {
  return gsl_sf_bessel_Knu_scaled((const double)x, (const double)y);
}

double gsl_sf_bessel_Knu__(double const &x, double const &y) {
  return gsl_sf_bessel_Knu((const double)x, (const double)y);
}

double gsl_sf_bessel_lnKnu__(double const &x, double const &y) {
  return gsl_sf_bessel_lnKnu((const double)x, (const double)y);
}

double gsl_sf_bessel_zero_J0__(long const &x) { return gsl_sf_bessel_zero_J0((unsigned int)x); }

double gsl_sf_bessel_zero_J1__(long const &x) { return gsl_sf_bessel_zero_J1((unsigned int)x); }

double gsl_sf_bessel_zero_Jnu__(double const &x, long const &y) {
  return gsl_sf_bessel_zero_Jnu((double)x, (unsigned int)y);
}

double gsl_sf_clausen__(double const &x) { return gsl_sf_clausen((const double)x); }

double gsl_sf_hydrogenicR_1__(double const &x, double const &y) {
  return gsl_sf_hydrogenicR_1((const double)x, (const double)y);
}

double gsl_sf_dawson__(double const &x) { return gsl_sf_dawson((double)x); }

double gsl_sf_debye_1__(double const &x) { return gsl_sf_debye_1((const double)x); }

double gsl_sf_debye_2__(double const &x) { return gsl_sf_debye_2((const double)x); }

double gsl_sf_debye_3__(double const &x) { return gsl_sf_debye_3((const double)x); }

double gsl_sf_debye_4__(double const &x) { return gsl_sf_debye_4((const double)x); }

double gsl_sf_debye_5__(double const &x) { return gsl_sf_debye_5((const double)x); }

double gsl_sf_debye_6__(double const &x) { return gsl_sf_debye_6((const double)x); }

double gsl_sf_dilog__(double const &x) { return gsl_sf_dilog((const double)x); }

double gsl_sf_multiply__(double const &x, double const &y) {
  return gsl_sf_multiply((const double)x, (const double)y);
}

double gsl_sf_ellint_Kcomp__(double const &x, long const &y) {
  return gsl_sf_ellint_Kcomp((double)x, (gsl_mode_t)y);
}

double gsl_sf_ellint_Ecomp__(double const &x, long const &y) {
  return gsl_sf_ellint_Ecomp((double)x, (gsl_mode_t)y);
}

double gsl_sf_ellint_Pcomp__(double const &x, double const &y, long const &z) {
  return gsl_sf_ellint_Pcomp((double)x, (double)y, (gsl_mode_t)z);
}

double gsl_sf_ellint_Dcomp__(double const &x, long const &y) {
  return gsl_sf_ellint_Dcomp((double)x, (gsl_mode_t)y);
}

double gsl_sf_ellint_F__(double const &x, double const &y, long const &z) {
  return gsl_sf_ellint_F((double)x, (double)y, (gsl_mode_t)z);
}

double gsl_sf_ellint_E__(double const &x, double const &y, long const &z) {
  return gsl_sf_ellint_E((double)x, (double)y, (gsl_mode_t)z);
}

double gsl_sf_ellint_RC__(double const &x, double const &y, long const &z) {
  return gsl_sf_ellint_RC((double)x, (double)y, (gsl_mode_t)z);
}

double gsl_sf_erfc__(double const &x) { return gsl_sf_erfc((double)x); }

double gsl_sf_log_erfc__(double const &x) { return gsl_sf_log_erfc((double)x); }

double gsl_sf_erf__(double const &x) { return gsl_sf_erf((double)x); }

double gsl_sf_erf_Z__(double const &x) { return gsl_sf_erf_Z((double)x); }

double gsl_sf_erf_Q__(double const &x) { return gsl_sf_erf_Q((double)x); }

double gsl_sf_hazard__(double const &x) { return gsl_sf_hazard((double)x); }

double gsl_sf_exp__(double const &x) { return gsl_sf_exp((const double)x); }

double gsl_sf_exp_mult__(double const &x, double const &y) {
  return gsl_sf_exp_mult((const double)x, (const double)y);
}

double gsl_sf_expm1__(double const &x) { return gsl_sf_expm1((const double)x); }

double gsl_sf_exprel__(double const &x) { return gsl_sf_exprel((const double)x); }

double gsl_sf_exprel_2__(double const &x) { return gsl_sf_exprel_2((const double)x); }

double gsl_sf_exprel_n__(long const &x, double const &y) {
  return gsl_sf_exprel_n((const int)x, (const double)y);
}

double gsl_sf_expint_E1__(double const &x) { return gsl_sf_expint_E1((const double)x); }

double gsl_sf_expint_E2__(double const &x) { return gsl_sf_expint_E2((const double)x); }

double gsl_sf_expint_En__(long const &x, double const &y) {
  return gsl_sf_expint_En((const int)x, (const double)y);
}

double gsl_sf_expint_E1_scaled__(double const &x) {
  return gsl_sf_expint_E1_scaled((const double)x);
}

double gsl_sf_expint_E2_scaled__(double const &x) {
  return gsl_sf_expint_E2_scaled((const double)x);
}

double gsl_sf_expint_En_scaled__(long const &x, double const &y) {
  return gsl_sf_expint_En_scaled((const int)x, (const double)y);
}

double gsl_sf_expint_Ei__(double const &x) { return gsl_sf_expint_Ei((const double)x); }

double gsl_sf_expint_Ei_scaled__(double const &x) {
  return gsl_sf_expint_Ei_scaled((const double)x);
}

double gsl_sf_Shi__(double const &x) { return gsl_sf_Shi((const double)x); }

double gsl_sf_Chi__(double const &x) { return gsl_sf_Chi((const double)x); }

double gsl_sf_expint_3__(double const &x) { return gsl_sf_expint_3((double)x); }

double gsl_sf_Si__(double const &x) { return gsl_sf_Si((const double)x); }

double gsl_sf_Ci__(double const &x) { return gsl_sf_Ci((const double)x); }

double gsl_sf_atanint__(double const &x) { return gsl_sf_atanint((const double)x); }

double gsl_sf_fermi_dirac_m1__(double const &x) { return gsl_sf_fermi_dirac_m1((const double)x); }

double gsl_sf_fermi_dirac_0__(double const &x) { return gsl_sf_fermi_dirac_0((const double)x); }

double gsl_sf_fermi_dirac_1__(double const &x) { return gsl_sf_fermi_dirac_1((const double)x); }

double gsl_sf_fermi_dirac_2__(double const &x) { return gsl_sf_fermi_dirac_2((const double)x); }

double gsl_sf_fermi_dirac_int__(long const &x, double const &y) {
  return gsl_sf_fermi_dirac_int((const int)x, (const double)y);
}

double gsl_sf_fermi_dirac_mhalf__(double const &x) {
  return gsl_sf_fermi_dirac_mhalf((const double)x);
}

double gsl_sf_fermi_dirac_half__(double const &x) {
  return gsl_sf_fermi_dirac_half((const double)x);
}

double gsl_sf_fermi_dirac_3half__(double const &x) {
  return gsl_sf_fermi_dirac_3half((const double)x);
}

double gsl_sf_fermi_dirac_inc_0__(double const &x, double const &y) {
  return gsl_sf_fermi_dirac_inc_0((const double)x, (const double)y);
}

double gsl_sf_lngamma__(double const &x) { return gsl_sf_lngamma((const double)x); }

double gsl_sf_gamma__(double const &x) { return gsl_sf_gamma((const double)x); }

double gsl_sf_gammastar__(double const &x) { return gsl_sf_gammastar((const double)x); }

double gsl_sf_gammainv__(double const &x) { return gsl_sf_gammainv((const double)x); }

double gsl_sf_taylorcoeff__(long const &x, double const &y) {
  return gsl_sf_taylorcoeff((const int)x, (const double)y);
}

double gsl_sf_fact__(long const &x) { return gsl_sf_fact((const unsigned int)x); }

double gsl_sf_doublefact__(long const &x) { return gsl_sf_doublefact((const unsigned int)x); }

double gsl_sf_lnfact__(long const &x) { return gsl_sf_lnfact((const unsigned int)x); }

double gsl_sf_lndoublefact__(long const &x) { return gsl_sf_lndoublefact((const unsigned int)x); }

double gsl_sf_lnchoose__(long const &x, long const &y) {
  return gsl_sf_lnchoose((unsigned int)x, (unsigned int)y);
}

double gsl_sf_choose__(long const &x, long const &y) {
  return gsl_sf_choose((unsigned int)x, (unsigned int)y);
}

double gsl_sf_lnpoch__(double const &x, double const &y) {
  return gsl_sf_lnpoch((const double)x, (const double)y);
}

double gsl_sf_poch__(double const &x, double const &y) {
  return gsl_sf_poch((const double)x, (const double)y);
}

double gsl_sf_pochrel__(double const &x, double const &y) {
  return gsl_sf_pochrel((const double)x, (const double)y);
}

double gsl_sf_gamma_inc_Q__(double const &x, double const &y) {
  return gsl_sf_gamma_inc_Q((const double)x, (const double)y);
}

double gsl_sf_gamma_inc_P__(double const &x, double const &y) {
  return gsl_sf_gamma_inc_P((const double)x, (const double)y);
}

double gsl_sf_gamma_inc__(double const &x, double const &y) {
  return gsl_sf_gamma_inc((const double)x, (const double)y);
}

double gsl_sf_lnbeta__(double const &x, double const &y) {
  return gsl_sf_lnbeta((const double)x, (const double)y);
}

double gsl_sf_beta__(double const &x, double const &y) {
  return gsl_sf_beta((const double)x, (const double)y);
}

double gsl_sf_beta_inc__(double const &x, double const &y, double const &z) {
  return gsl_sf_beta_inc((const double)x, (const double)y, (const double)z);
}

double gsl_sf_gegenpoly_1__(double const &x, double const &y) {
  return gsl_sf_gegenpoly_1((double)x, (double)y);
}

double gsl_sf_gegenpoly_2__(double const &x, double const &y) {
  return gsl_sf_gegenpoly_2((double)x, (double)y);
}

double gsl_sf_gegenpoly_3__(double const &x, double const &y) {
  return gsl_sf_gegenpoly_3((double)x, (double)y);
}

double gsl_sf_gegenpoly_n__(long const &x, double const &y, double const &z) {
  return gsl_sf_gegenpoly_n((int)x, (double)y, (double)z);
}

double gsl_sf_hyperg_0F1__(double const &x, double const &y) {
  return gsl_sf_hyperg_0F1((const double)x, (const double)y);
}

double gsl_sf_hyperg_1F1_int__(long const &x, long const &y, double const &z) {
  return gsl_sf_hyperg_1F1_int((const int)x, (const int)y, (double)z);
}

double gsl_sf_hyperg_1F1__(double const &x, double const &y, double const &z) {
  return gsl_sf_hyperg_1F1((double)x, (double)y, (double)z);
}

double gsl_sf_hyperg_U_int__(long const &x, long const &y, double const &z) {
  return gsl_sf_hyperg_U_int((const int)x, (const int)y, (const double)z);
}

double gsl_sf_hyperg_U__(double const &x, double const &y, double const &z) {
  return gsl_sf_hyperg_U((const double)x, (const double)y, (const double)z);
}

double gsl_sf_hyperg_2F0__(double const &x, double const &y, double const &z) {
  return gsl_sf_hyperg_2F0((const double)x, (const double)y, (const double)z);
}

double gsl_sf_laguerre_1__(double const &x, double const &y) {
  return gsl_sf_laguerre_1((double)x, (double)y);
}

double gsl_sf_laguerre_2__(double const &x, double const &y) {
  return gsl_sf_laguerre_2((double)x, (double)y);
}

double gsl_sf_laguerre_3__(double const &x, double const &y) {
  return gsl_sf_laguerre_3((double)x, (double)y);
}

double gsl_sf_laguerre_n__(long const &x, double const &y, double const &z) {
  return gsl_sf_laguerre_n((int)x, (double)y, (double)z);
}

double gsl_sf_lambert_W0__(double const &x) { return gsl_sf_lambert_W0((double)x); }

double gsl_sf_lambert_Wm1__(double const &x) { return gsl_sf_lambert_Wm1((double)x); }

double gsl_sf_legendre_Pl__(long const &x, double const &y) {
  return gsl_sf_legendre_Pl((const int)x, (const double)y);
}

double gsl_sf_legendre_P1__(double const &x) { return gsl_sf_legendre_P1((const double)x); }

double gsl_sf_legendre_P2__(double const &x) { return gsl_sf_legendre_P2((const double)x); }

double gsl_sf_legendre_P3__(double const &x) { return gsl_sf_legendre_P3((const double)x); }

double gsl_sf_legendre_Q0__(double const &x) { return gsl_sf_legendre_Q0((const double)x); }

double gsl_sf_legendre_Q1__(double const &x) { return gsl_sf_legendre_Q1((const double)x); }

double gsl_sf_legendre_Ql__(long const &x, double const &y) {
  return gsl_sf_legendre_Ql((const int)x, (const double)y);
}

double gsl_sf_legendre_Plm__(long const &x, long const &y, double const &z) {
  return gsl_sf_legendre_Plm((const int)x, (const int)y, (const double)z);
}

double gsl_sf_legendre_sphPlm__(long const &x, long const &y, double const &z) {
  return gsl_sf_legendre_sphPlm((const int)x, (const int)y, (const double)z);
}

// long gsl_sf_legendre_array_size__(long const & x , long const & y ){ return
// gsl_sf_legendre_array_size( (const int) x , (const int) y );}
double gsl_sf_conicalP_half__(double const &x, double const &y) {
  return gsl_sf_conicalP_half((const double)x, (const double)y);
}

double gsl_sf_conicalP_mhalf__(double const &x, double const &y) {
  return gsl_sf_conicalP_mhalf((const double)x, (const double)y);
}

double gsl_sf_conicalP_0__(double const &x, double const &y) {
  return gsl_sf_conicalP_0((const double)x, (const double)y);
}

double gsl_sf_conicalP_1__(double const &x, double const &y) {
  return gsl_sf_conicalP_1((const double)x, (const double)y);
}

double gsl_sf_conicalP_sph_reg__(long const &x, double const &y, double const &z) {
  return gsl_sf_conicalP_sph_reg((const int)x, (const double)y, (const double)z);
}

double gsl_sf_conicalP_cyl_reg__(long const &x, double const &y, double const &z) {
  return gsl_sf_conicalP_cyl_reg((const int)x, (const double)y, (const double)z);
}

double gsl_sf_legendre_H3d_0__(double const &x, double const &y) {
  return gsl_sf_legendre_H3d_0((const double)x, (const double)y);
}

double gsl_sf_legendre_H3d_1__(double const &x, double const &y) {
  return gsl_sf_legendre_H3d_1((const double)x, (const double)y);
}

double gsl_sf_legendre_H3d__(long const &x, double const &y, double const &z) {
  return gsl_sf_legendre_H3d((const int)x, (const double)y, (const double)z);
}

double gsl_sf_log__(double const &x) { return gsl_sf_log((const double)x); }

double gsl_sf_log_abs__(double const &x) { return gsl_sf_log_abs((const double)x); }

double gsl_sf_log_1plusx__(double const &x) { return gsl_sf_log_1plusx((const double)x); }

double gsl_sf_log_1plusx_mx__(double const &x) { return gsl_sf_log_1plusx_mx((const double)x); }

double gsl_sf_pow_int__(double const &x, long const &y) {
  return gsl_sf_pow_int((const double)x, (const int)y);
}

double gsl_sf_psi_int__(long const &x) { return gsl_sf_psi_int((const int)x); }

double gsl_sf_psi__(double const &x) { return gsl_sf_psi((const double)x); }

double gsl_sf_psi_1piy__(double const &x) { return gsl_sf_psi_1piy((const double)x); }

double gsl_sf_psi_1_int__(long const &x) { return gsl_sf_psi_1_int((const int)x); }

double gsl_sf_psi_1__(double const &x) { return gsl_sf_psi_1((const double)x); }

double gsl_sf_psi_n__(long const &x, double const &y) {
  return gsl_sf_psi_n((const int)x, (const double)y);
}

double gsl_sf_synchrotron_1__(double const &x) { return gsl_sf_synchrotron_1((const double)x); }

double gsl_sf_synchrotron_2__(double const &x) { return gsl_sf_synchrotron_2((const double)x); }

double gsl_sf_transport_2__(double const &x) { return gsl_sf_transport_2((const double)x); }

double gsl_sf_transport_3__(double const &x) { return gsl_sf_transport_3((const double)x); }

double gsl_sf_transport_4__(double const &x) { return gsl_sf_transport_4((const double)x); }

double gsl_sf_transport_5__(double const &x) { return gsl_sf_transport_5((const double)x); }

double gsl_sf_sin__(double const &x) { return gsl_sf_sin((const double)x); }

double gsl_sf_cos__(double const &x) { return gsl_sf_cos((const double)x); }

double gsl_sf_hypot__(double const &x, double const &y) {
  return gsl_sf_hypot((const double)x, (const double)y);
}

double gsl_sf_sinc__(double const &x) { return gsl_sf_sinc((const double)x); }

double gsl_sf_lnsinh__(double const &x) { return gsl_sf_lnsinh((const double)x); }

double gsl_sf_lncosh__(double const &x) { return gsl_sf_lncosh((const double)x); }

double gsl_sf_angle_restrict_symm__(double const &x) {
  return gsl_sf_angle_restrict_symm((const double)x);
}

double gsl_sf_angle_restrict_pos__(double const &x) {
  return gsl_sf_angle_restrict_pos((const double)x);
}

double gsl_sf_zeta_int__(long const &x) { return gsl_sf_zeta_int((const int)x); }

double gsl_sf_zeta__(double const &x) { return gsl_sf_zeta((const double)x); }

double gsl_sf_zetam1__(double const &x) { return gsl_sf_zetam1((const double)x); }

double gsl_sf_zetam1_int__(long const &x) { return gsl_sf_zetam1_int((const int)x); }

double gsl_sf_hzeta__(double const &x, double const &y) {
  return gsl_sf_hzeta((const double)x, (const double)y);
}

double gsl_sf_eta_int__(long const &x) { return gsl_sf_eta_int((const int)x); }

double gsl_sf_eta__(double const &x) { return gsl_sf_eta((const double)x); }

/*****************/
/*****************/
void init_gsl_sf( ) {
  Global.Add("gslcdfugaussianP", "(", new OneOperator1_< double, double >(gsl_cdf_ugaussian_P__));
  Global.Add("gslcdfugaussianQ", "(", new OneOperator1_< double, double >(gsl_cdf_ugaussian_Q__));
  Global.Add("gslcdfugaussianPinv", "(",
             new OneOperator1_< double, double >(gsl_cdf_ugaussian_Pinv__));
  Global.Add("gslcdfugaussianQinv", "(",
             new OneOperator1_< double, double >(gsl_cdf_ugaussian_Qinv__));
  Global.Add("gslcdfgaussianP", "(",
             new OneOperator2_< double, double, double >(gsl_cdf_gaussian_P__));
  Global.Add("gslcdfgaussianQ", "(",
             new OneOperator2_< double, double, double >(gsl_cdf_gaussian_Q__));
  Global.Add("gslcdfgaussianPinv", "(",
             new OneOperator2_< double, double, double >(gsl_cdf_gaussian_Pinv__));
  Global.Add("gslcdfgaussianQinv", "(",
             new OneOperator2_< double, double, double >(gsl_cdf_gaussian_Qinv__));
  Global.Add("gslcdfgammaP", "(",
             new OneOperator3_< double, double, double, double >(gsl_cdf_gamma_P__));
  Global.Add("gslcdfgammaQ", "(",
             new OneOperator3_< double, double, double, double >(gsl_cdf_gamma_Q__));
  Global.Add("gslcdfgammaPinv", "(",
             new OneOperator3_< double, double, double, double >(gsl_cdf_gamma_Pinv__));
  Global.Add("gslcdfgammaQinv", "(",
             new OneOperator3_< double, double, double, double >(gsl_cdf_gamma_Qinv__));
  Global.Add("gslcdfcauchyP", "(", new OneOperator2_< double, double, double >(gsl_cdf_cauchy_P__));
  Global.Add("gslcdfcauchyQ", "(", new OneOperator2_< double, double, double >(gsl_cdf_cauchy_Q__));
  Global.Add("gslcdfcauchyPinv", "(",
             new OneOperator2_< double, double, double >(gsl_cdf_cauchy_Pinv__));
  Global.Add("gslcdfcauchyQinv", "(",
             new OneOperator2_< double, double, double >(gsl_cdf_cauchy_Qinv__));
  Global.Add("gslcdflaplaceP", "(",
             new OneOperator2_< double, double, double >(gsl_cdf_laplace_P__));
  Global.Add("gslcdflaplaceQ", "(",
             new OneOperator2_< double, double, double >(gsl_cdf_laplace_Q__));
  Global.Add("gslcdflaplacePinv", "(",
             new OneOperator2_< double, double, double >(gsl_cdf_laplace_Pinv__));
  Global.Add("gslcdflaplaceQinv", "(",
             new OneOperator2_< double, double, double >(gsl_cdf_laplace_Qinv__));
  Global.Add("gslcdfrayleighP", "(",
             new OneOperator2_< double, double, double >(gsl_cdf_rayleigh_P__));
  Global.Add("gslcdfrayleighQ", "(",
             new OneOperator2_< double, double, double >(gsl_cdf_rayleigh_Q__));
  Global.Add("gslcdfrayleighPinv", "(",
             new OneOperator2_< double, double, double >(gsl_cdf_rayleigh_Pinv__));
  Global.Add("gslcdfrayleighQinv", "(",
             new OneOperator2_< double, double, double >(gsl_cdf_rayleigh_Qinv__));
  Global.Add("gslcdfchisqP", "(", new OneOperator2_< double, double, double >(gsl_cdf_chisq_P__));
  Global.Add("gslcdfchisqQ", "(", new OneOperator2_< double, double, double >(gsl_cdf_chisq_Q__));
  Global.Add("gslcdfchisqPinv", "(",
             new OneOperator2_< double, double, double >(gsl_cdf_chisq_Pinv__));
  Global.Add("gslcdfchisqQinv", "(",
             new OneOperator2_< double, double, double >(gsl_cdf_chisq_Qinv__));
  Global.Add("gslcdfexponentialP", "(",
             new OneOperator2_< double, double, double >(gsl_cdf_exponential_P__));
  Global.Add("gslcdfexponentialQ", "(",
             new OneOperator2_< double, double, double >(gsl_cdf_exponential_Q__));
  Global.Add("gslcdfexponentialPinv", "(",
             new OneOperator2_< double, double, double >(gsl_cdf_exponential_Pinv__));
  Global.Add("gslcdfexponentialQinv", "(",
             new OneOperator2_< double, double, double >(gsl_cdf_exponential_Qinv__));
  Global.Add("gslcdfexppowP", "(",
             new OneOperator3_< double, double, double, double >(gsl_cdf_exppow_P__));
  Global.Add("gslcdfexppowQ", "(",
             new OneOperator3_< double, double, double, double >(gsl_cdf_exppow_Q__));
  Global.Add("gslcdftdistP", "(", new OneOperator2_< double, double, double >(gsl_cdf_tdist_P__));
  Global.Add("gslcdftdistQ", "(", new OneOperator2_< double, double, double >(gsl_cdf_tdist_Q__));
  Global.Add("gslcdftdistPinv", "(",
             new OneOperator2_< double, double, double >(gsl_cdf_tdist_Pinv__));
  Global.Add("gslcdftdistQinv", "(",
             new OneOperator2_< double, double, double >(gsl_cdf_tdist_Qinv__));
  Global.Add("gslcdffdistP", "(",
             new OneOperator3_< double, double, double, double >(gsl_cdf_fdist_P__));
  Global.Add("gslcdffdistQ", "(",
             new OneOperator3_< double, double, double, double >(gsl_cdf_fdist_Q__));
  Global.Add("gslcdffdistPinv", "(",
             new OneOperator3_< double, double, double, double >(gsl_cdf_fdist_Pinv__));
  Global.Add("gslcdffdistQinv", "(",
             new OneOperator3_< double, double, double, double >(gsl_cdf_fdist_Qinv__));
  Global.Add("gslcdfbetaP", "(",
             new OneOperator3_< double, double, double, double >(gsl_cdf_beta_P__));
  Global.Add("gslcdfbetaQ", "(",
             new OneOperator3_< double, double, double, double >(gsl_cdf_beta_Q__));
  Global.Add("gslcdfbetaPinv", "(",
             new OneOperator3_< double, double, double, double >(gsl_cdf_beta_Pinv__));
  Global.Add("gslcdfbetaQinv", "(",
             new OneOperator3_< double, double, double, double >(gsl_cdf_beta_Qinv__));
  Global.Add("gslcdfflatP", "(",
             new OneOperator3_< double, double, double, double >(gsl_cdf_flat_P__));
  Global.Add("gslcdfflatQ", "(",
             new OneOperator3_< double, double, double, double >(gsl_cdf_flat_Q__));
  Global.Add("gslcdfflatPinv", "(",
             new OneOperator3_< double, double, double, double >(gsl_cdf_flat_Pinv__));
  Global.Add("gslcdfflatQinv", "(",
             new OneOperator3_< double, double, double, double >(gsl_cdf_flat_Qinv__));
  Global.Add("gslcdflognormalP", "(",
             new OneOperator3_< double, double, double, double >(gsl_cdf_lognormal_P__));
  Global.Add("gslcdflognormalQ", "(",
             new OneOperator3_< double, double, double, double >(gsl_cdf_lognormal_Q__));
  Global.Add("gslcdflognormalPinv", "(",
             new OneOperator3_< double, double, double, double >(gsl_cdf_lognormal_Pinv__));
  Global.Add("gslcdflognormalQinv", "(",
             new OneOperator3_< double, double, double, double >(gsl_cdf_lognormal_Qinv__));
  Global.Add("gslcdfgumbel1P", "(",
             new OneOperator3_< double, double, double, double >(gsl_cdf_gumbel1_P__));
  Global.Add("gslcdfgumbel1Q", "(",
             new OneOperator3_< double, double, double, double >(gsl_cdf_gumbel1_Q__));
  Global.Add("gslcdfgumbel1Pinv", "(",
             new OneOperator3_< double, double, double, double >(gsl_cdf_gumbel1_Pinv__));
  Global.Add("gslcdfgumbel1Qinv", "(",
             new OneOperator3_< double, double, double, double >(gsl_cdf_gumbel1_Qinv__));
  Global.Add("gslcdfgumbel2P", "(",
             new OneOperator3_< double, double, double, double >(gsl_cdf_gumbel2_P__));
  Global.Add("gslcdfgumbel2Q", "(",
             new OneOperator3_< double, double, double, double >(gsl_cdf_gumbel2_Q__));
  Global.Add("gslcdfgumbel2Pinv", "(",
             new OneOperator3_< double, double, double, double >(gsl_cdf_gumbel2_Pinv__));
  Global.Add("gslcdfgumbel2Qinv", "(",
             new OneOperator3_< double, double, double, double >(gsl_cdf_gumbel2_Qinv__));
  Global.Add("gslcdfweibullP", "(",
             new OneOperator3_< double, double, double, double >(gsl_cdf_weibull_P__));
  Global.Add("gslcdfweibullQ", "(",
             new OneOperator3_< double, double, double, double >(gsl_cdf_weibull_Q__));
  Global.Add("gslcdfweibullPinv", "(",
             new OneOperator3_< double, double, double, double >(gsl_cdf_weibull_Pinv__));
  Global.Add("gslcdfweibullQinv", "(",
             new OneOperator3_< double, double, double, double >(gsl_cdf_weibull_Qinv__));
  Global.Add("gslcdfparetoP", "(",
             new OneOperator3_< double, double, double, double >(gsl_cdf_pareto_P__));
  Global.Add("gslcdfparetoQ", "(",
             new OneOperator3_< double, double, double, double >(gsl_cdf_pareto_Q__));
  Global.Add("gslcdfparetoPinv", "(",
             new OneOperator3_< double, double, double, double >(gsl_cdf_pareto_Pinv__));
  Global.Add("gslcdfparetoQinv", "(",
             new OneOperator3_< double, double, double, double >(gsl_cdf_pareto_Qinv__));
  Global.Add("gslcdflogisticP", "(",
             new OneOperator2_< double, double, double >(gsl_cdf_logistic_P__));
  Global.Add("gslcdflogisticQ", "(",
             new OneOperator2_< double, double, double >(gsl_cdf_logistic_Q__));
  Global.Add("gslcdflogisticPinv", "(",
             new OneOperator2_< double, double, double >(gsl_cdf_logistic_Pinv__));
  Global.Add("gslcdflogisticQinv", "(",
             new OneOperator2_< double, double, double >(gsl_cdf_logistic_Qinv__));
  Global.Add("gslcdfbinomialP", "(",
             new OneOperator3_< double, long, double, long >(gsl_cdf_binomial_P__));
  Global.Add("gslcdfbinomialQ", "(",
             new OneOperator3_< double, long, double, long >(gsl_cdf_binomial_Q__));
  Global.Add("gslcdfpoissonP", "(", new OneOperator2_< double, long, double >(gsl_cdf_poisson_P__));
  Global.Add("gslcdfpoissonQ", "(", new OneOperator2_< double, long, double >(gsl_cdf_poisson_Q__));
  Global.Add("gslcdfgeometricP", "(",
             new OneOperator2_< double, long, double >(gsl_cdf_geometric_P__));
  Global.Add("gslcdfgeometricQ", "(",
             new OneOperator2_< double, long, double >(gsl_cdf_geometric_Q__));
  Global.Add("gslcdfnegativebinomialP", "(",
             new OneOperator3_< double, long, double, double >(gsl_cdf_negative_binomial_P__));
  Global.Add("gslcdfnegativebinomialQ", "(",
             new OneOperator3_< double, long, double, double >(gsl_cdf_negative_binomial_Q__));
  Global.Add("gslcdfpascalP", "(",
             new OneOperator3_< double, long, double, long >(gsl_cdf_pascal_P__));
  Global.Add("gslcdfpascalQ", "(",
             new OneOperator3_< double, long, double, long >(gsl_cdf_pascal_Q__));
  Global.Add("gslranbernoullipdf", "(",
             new OneOperator2_< double, long, double >(gsl_ran_bernoulli_pdf__));
  Global.Add("gslranbeta", "(",
             new OneOperator3_< double, gsl_rng **, double, double >(gsl_ran_beta__));
  Global.Add("gslranbetapdf", "(",
             new OneOperator3_< double, double, double, double >(gsl_ran_beta_pdf__));
  Global.Add("gslranbinomialpdf", "(",
             new OneOperator3_< double, long, double, long >(gsl_ran_binomial_pdf__));
  Global.Add("gslranexponential", "(",
             new OneOperator2_< double, gsl_rng **, double >(gsl_ran_exponential__));
  Global.Add("gslranexponentialpdf", "(",
             new OneOperator2_< double, double, double >(gsl_ran_exponential_pdf__));
  Global.Add("gslranexppow", "(",
             new OneOperator3_< double, gsl_rng **, double, double >(gsl_ran_exppow__));
  Global.Add("gslranexppowpdf", "(",
             new OneOperator3_< double, double, double, double >(gsl_ran_exppow_pdf__));
  Global.Add("gslrancauchy", "(",
             new OneOperator2_< double, gsl_rng **, double >(gsl_ran_cauchy__));
  Global.Add("gslrancauchypdf", "(",
             new OneOperator2_< double, double, double >(gsl_ran_cauchy_pdf__));
  Global.Add("gslranchisq", "(", new OneOperator2_< double, gsl_rng **, double >(gsl_ran_chisq__));
  Global.Add("gslranchisqpdf", "(",
             new OneOperator2_< double, double, double >(gsl_ran_chisq_pdf__));
  Global.Add("gslranerlang", "(",
             new OneOperator3_< double, gsl_rng **, double, double >(gsl_ran_erlang__));
  Global.Add("gslranerlangpdf", "(",
             new OneOperator3_< double, double, double, double >(gsl_ran_erlang_pdf__));
  Global.Add("gslranfdist", "(",
             new OneOperator3_< double, gsl_rng **, double, double >(gsl_ran_fdist__));
  Global.Add("gslranfdistpdf", "(",
             new OneOperator3_< double, double, double, double >(gsl_ran_fdist_pdf__));
  Global.Add("gslranflat", "(",
             new OneOperator3_< double, gsl_rng **, double, double >(gsl_ran_flat__));
  Global.Add("gslranflatpdf", "(",
             new OneOperator3_< double, double, double, double >(gsl_ran_flat_pdf__));
  Global.Add("gslrangamma", "(",
             new OneOperator3_< double, gsl_rng **, double, double >(gsl_ran_gamma__));
  Global.Add("gslrangammaint", "(",
             new OneOperator2_< double, gsl_rng **, long >(gsl_ran_gamma_int__));
  Global.Add("gslrangammapdf", "(",
             new OneOperator3_< double, double, double, double >(gsl_ran_gamma_pdf__));
  Global.Add("gslrangammamt", "(",
             new OneOperator3_< double, gsl_rng **, double, double >(gsl_ran_gamma_mt__));
  Global.Add("gslrangammaknuth", "(",
             new OneOperator3_< double, gsl_rng **, double, double >(gsl_ran_gamma_knuth__));
  Global.Add("gslrangaussian", "(",
             new OneOperator2_< double, gsl_rng **, double >(gsl_ran_gaussian__));
  Global.Add("gslrangaussianratiomethod", "(",
             new OneOperator2_< double, gsl_rng **, double >(gsl_ran_gaussian_ratio_method__));
  Global.Add("gslrangaussianziggurat", "(",
             new OneOperator2_< double, gsl_rng **, double >(gsl_ran_gaussian_ziggurat__));
  Global.Add("gslrangaussianpdf", "(",
             new OneOperator2_< double, double, double >(gsl_ran_gaussian_pdf__));
  Global.Add("gslranugaussian", "(", new OneOperator1_< double, gsl_rng ** >(gsl_ran_ugaussian__));
  Global.Add("gslranugaussianratiomethod", "(",
             new OneOperator1_< double, gsl_rng ** >(gsl_ran_ugaussian_ratio_method__));
  Global.Add("gslranugaussianpdf", "(",
             new OneOperator1_< double, double >(gsl_ran_ugaussian_pdf__));
  Global.Add("gslrangaussiantail", "(",
             new OneOperator3_< double, gsl_rng **, double, double >(gsl_ran_gaussian_tail__));
  Global.Add("gslrangaussiantailpdf", "(",
             new OneOperator3_< double, double, double, double >(gsl_ran_gaussian_tail_pdf__));
  Global.Add("gslranugaussiantail", "(",
             new OneOperator2_< double, gsl_rng **, double >(gsl_ran_ugaussian_tail__));
  Global.Add("gslranugaussiantailpdf", "(",
             new OneOperator2_< double, double, double >(gsl_ran_ugaussian_tail_pdf__));
  Global.Add("gslranlandau", "(", new OneOperator1_< double, gsl_rng ** >(gsl_ran_landau__));
  Global.Add("gslranlandaupdf", "(", new OneOperator1_< double, double >(gsl_ran_landau_pdf__));
  Global.Add("gslrangeometricpdf", "(",
             new OneOperator2_< double, long, double >(gsl_ran_geometric_pdf__));
  Global.Add("gslrangumbel1", "(",
             new OneOperator3_< double, gsl_rng **, double, double >(gsl_ran_gumbel1__));
  Global.Add("gslrangumbel1pdf", "(",
             new OneOperator3_< double, double, double, double >(gsl_ran_gumbel1_pdf__));
  Global.Add("gslrangumbel2", "(",
             new OneOperator3_< double, gsl_rng **, double, double >(gsl_ran_gumbel2__));
  Global.Add("gslrangumbel2pdf", "(",
             new OneOperator3_< double, double, double, double >(gsl_ran_gumbel2_pdf__));
  Global.Add("gslranlogistic", "(",
             new OneOperator2_< double, gsl_rng **, double >(gsl_ran_logistic__));
  Global.Add("gslranlogisticpdf", "(",
             new OneOperator2_< double, double, double >(gsl_ran_logistic_pdf__));
  Global.Add("gslranlognormal", "(",
             new OneOperator3_< double, gsl_rng **, double, double >(gsl_ran_lognormal__));
  Global.Add("gslranlognormalpdf", "(",
             new OneOperator3_< double, double, double, double >(gsl_ran_lognormal_pdf__));
  Global.Add("gslranlogarithmicpdf", "(",
             new OneOperator2_< double, long, double >(gsl_ran_logarithmic_pdf__));
  Global.Add("gslrannegativebinomialpdf", "(",
             new OneOperator3_< double, long, double, double >(gsl_ran_negative_binomial_pdf__));
  Global.Add("gslranpascalpdf", "(",
             new OneOperator3_< double, long, double, long >(gsl_ran_pascal_pdf__));
  Global.Add("gslranpareto", "(",
             new OneOperator3_< double, gsl_rng **, double, double >(gsl_ran_pareto__));
  Global.Add("gslranparetopdf", "(",
             new OneOperator3_< double, double, double, double >(gsl_ran_pareto_pdf__));
  Global.Add("gslranpoissonpdf", "(",
             new OneOperator2_< double, long, double >(gsl_ran_poisson_pdf__));
  Global.Add("gslranrayleigh", "(",
             new OneOperator2_< double, gsl_rng **, double >(gsl_ran_rayleigh__));
  Global.Add("gslranrayleighpdf", "(",
             new OneOperator2_< double, double, double >(gsl_ran_rayleigh_pdf__));
  Global.Add("gslranrayleightail", "(",
             new OneOperator3_< double, gsl_rng **, double, double >(gsl_ran_rayleigh_tail__));
  Global.Add("gslranrayleightailpdf", "(",
             new OneOperator3_< double, double, double, double >(gsl_ran_rayleigh_tail_pdf__));
  Global.Add("gslrantdist", "(", new OneOperator2_< double, gsl_rng **, double >(gsl_ran_tdist__));
  Global.Add("gslrantdistpdf", "(",
             new OneOperator2_< double, double, double >(gsl_ran_tdist_pdf__));
  Global.Add("gslranlaplace", "(",
             new OneOperator2_< double, gsl_rng **, double >(gsl_ran_laplace__));
  Global.Add("gslranlaplacepdf", "(",
             new OneOperator2_< double, double, double >(gsl_ran_laplace_pdf__));
  Global.Add("gslranlevy", "(",
             new OneOperator3_< double, gsl_rng **, double, double >(gsl_ran_levy__));
  Global.Add("gslranweibull", "(",
             new OneOperator3_< double, gsl_rng **, double, double >(gsl_ran_weibull__));
  Global.Add("gslranweibullpdf", "(",
             new OneOperator3_< double, double, double, double >(gsl_ran_weibull_pdf__));
  Global.Add("gslsfairyAi", "(", new OneOperator2_< double, double, long >(gsl_sf_airy_Ai__));
  Global.Add("gslsfairyBi", "(", new OneOperator2_< double, double, long >(gsl_sf_airy_Bi__));
  Global.Add("gslsfairyAiscaled", "(",
             new OneOperator2_< double, double, long >(gsl_sf_airy_Ai_scaled__));
  Global.Add("gslsfairyBiscaled", "(",
             new OneOperator2_< double, double, long >(gsl_sf_airy_Bi_scaled__));
  Global.Add("gslsfairyAideriv", "(",
             new OneOperator2_< double, double, long >(gsl_sf_airy_Ai_deriv__));
  Global.Add("gslsfairyBideriv", "(",
             new OneOperator2_< double, double, long >(gsl_sf_airy_Bi_deriv__));
  Global.Add("gslsfairyAiderivscaled", "(",
             new OneOperator2_< double, double, long >(gsl_sf_airy_Ai_deriv_scaled__));
  Global.Add("gslsfairyBiderivscaled", "(",
             new OneOperator2_< double, double, long >(gsl_sf_airy_Bi_deriv_scaled__));
  Global.Add("gslsfairyzeroAi", "(", new OneOperator1_< double, long >(gsl_sf_airy_zero_Ai__));
  Global.Add("gslsfairyzeroBi", "(", new OneOperator1_< double, long >(gsl_sf_airy_zero_Bi__));
  Global.Add("gslsfairyzeroAideriv", "(",
             new OneOperator1_< double, long >(gsl_sf_airy_zero_Ai_deriv__));
  Global.Add("gslsfairyzeroBideriv", "(",
             new OneOperator1_< double, long >(gsl_sf_airy_zero_Bi_deriv__));
  Global.Add("gslsfbesselJ0", "(", new OneOperator1_< double, double >(gsl_sf_bessel_J0__));
  Global.Add("gslsfbesselJ1", "(", new OneOperator1_< double, double >(gsl_sf_bessel_J1__));
  Global.Add("gslsfbesselJn", "(", new OneOperator2_< double, long, double >(gsl_sf_bessel_Jn__));
  Global.Add("gslsfbesselY0", "(", new OneOperator1_< double, double >(gsl_sf_bessel_Y0__));
  Global.Add("gslsfbesselY1", "(", new OneOperator1_< double, double >(gsl_sf_bessel_Y1__));
  Global.Add("gslsfbesselYn", "(", new OneOperator2_< double, long, double >(gsl_sf_bessel_Yn__));
  Global.Add("gslsfbesselI0", "(", new OneOperator1_< double, double >(gsl_sf_bessel_I0__));
  Global.Add("gslsfbesselI1", "(", new OneOperator1_< double, double >(gsl_sf_bessel_I1__));
  Global.Add("gslsfbesselIn", "(", new OneOperator2_< double, long, double >(gsl_sf_bessel_In__));
  Global.Add("gslsfbesselI0scaled", "(",
             new OneOperator1_< double, double >(gsl_sf_bessel_I0_scaled__));
  Global.Add("gslsfbesselI1scaled", "(",
             new OneOperator1_< double, double >(gsl_sf_bessel_I1_scaled__));
  Global.Add("gslsfbesselInscaled", "(",
             new OneOperator2_< double, long, double >(gsl_sf_bessel_In_scaled__));
  Global.Add("gslsfbesselK0", "(", new OneOperator1_< double, double >(gsl_sf_bessel_K0__));
  Global.Add("gslsfbesselK1", "(", new OneOperator1_< double, double >(gsl_sf_bessel_K1__));
  Global.Add("gslsfbesselKn", "(", new OneOperator2_< double, long, double >(gsl_sf_bessel_Kn__));
  Global.Add("gslsfbesselK0scaled", "(",
             new OneOperator1_< double, double >(gsl_sf_bessel_K0_scaled__));
  Global.Add("gslsfbesselK1scaled", "(",
             new OneOperator1_< double, double >(gsl_sf_bessel_K1_scaled__));
  Global.Add("gslsfbesselKnscaled", "(",
             new OneOperator2_< double, long, double >(gsl_sf_bessel_Kn_scaled__));
  Global.Add("gslsfbesselj0", "(", new OneOperator1_< double, double >(gsl_sf_bessel_j0__));
  Global.Add("gslsfbesselj1", "(", new OneOperator1_< double, double >(gsl_sf_bessel_j1__));
  Global.Add("gslsfbesselj2", "(", new OneOperator1_< double, double >(gsl_sf_bessel_j2__));
  Global.Add("gslsfbesseljl", "(", new OneOperator2_< double, long, double >(gsl_sf_bessel_jl__));
  Global.Add("gslsfbessely0", "(", new OneOperator1_< double, double >(gsl_sf_bessel_y0__));
  Global.Add("gslsfbessely1", "(", new OneOperator1_< double, double >(gsl_sf_bessel_y1__));
  Global.Add("gslsfbessely2", "(", new OneOperator1_< double, double >(gsl_sf_bessel_y2__));
  Global.Add("gslsfbesselyl", "(", new OneOperator2_< double, long, double >(gsl_sf_bessel_yl__));
  Global.Add("gslsfbesseli0scaled", "(",
             new OneOperator1_< double, double >(gsl_sf_bessel_i0_scaled__));
  Global.Add("gslsfbesseli1scaled", "(",
             new OneOperator1_< double, double >(gsl_sf_bessel_i1_scaled__));
  Global.Add("gslsfbesseli2scaled", "(",
             new OneOperator1_< double, double >(gsl_sf_bessel_i2_scaled__));
  Global.Add("gslsfbesselilscaled", "(",
             new OneOperator2_< double, long, double >(gsl_sf_bessel_il_scaled__));
  Global.Add("gslsfbesselk0scaled", "(",
             new OneOperator1_< double, double >(gsl_sf_bessel_k0_scaled__));
  Global.Add("gslsfbesselk1scaled", "(",
             new OneOperator1_< double, double >(gsl_sf_bessel_k1_scaled__));
  Global.Add("gslsfbesselk2scaled", "(",
             new OneOperator1_< double, double >(gsl_sf_bessel_k2_scaled__));
  Global.Add("gslsfbesselklscaled", "(",
             new OneOperator2_< double, long, double >(gsl_sf_bessel_kl_scaled__));
  Global.Add("gslsfbesselJnu", "(",
             new OneOperator2_< double, double, double >(gsl_sf_bessel_Jnu__));
  Global.Add("gslsfbesselYnu", "(",
             new OneOperator2_< double, double, double >(gsl_sf_bessel_Ynu__));
  Global.Add("gslsfbesselInuscaled", "(",
             new OneOperator2_< double, double, double >(gsl_sf_bessel_Inu_scaled__));
  Global.Add("gslsfbesselInu", "(",
             new OneOperator2_< double, double, double >(gsl_sf_bessel_Inu__));
  Global.Add("gslsfbesselKnuscaled", "(",
             new OneOperator2_< double, double, double >(gsl_sf_bessel_Knu_scaled__));
  Global.Add("gslsfbesselKnu", "(",
             new OneOperator2_< double, double, double >(gsl_sf_bessel_Knu__));
  Global.Add("gslsfbessellnKnu", "(",
             new OneOperator2_< double, double, double >(gsl_sf_bessel_lnKnu__));
  Global.Add("gslsfbesselzeroJ0", "(", new OneOperator1_< double, long >(gsl_sf_bessel_zero_J0__));
  Global.Add("gslsfbesselzeroJ1", "(", new OneOperator1_< double, long >(gsl_sf_bessel_zero_J1__));
  Global.Add("gslsfbesselzeroJnu", "(",
             new OneOperator2_< double, double, long >(gsl_sf_bessel_zero_Jnu__));
  Global.Add("gslsfclausen", "(", new OneOperator1_< double, double >(gsl_sf_clausen__));
  Global.Add("gslsfhydrogenicR1", "(",
             new OneOperator2_< double, double, double >(gsl_sf_hydrogenicR_1__));
  Global.Add("gslsfdawson", "(", new OneOperator1_< double, double >(gsl_sf_dawson__));
  Global.Add("gslsfdebye1", "(", new OneOperator1_< double, double >(gsl_sf_debye_1__));
  Global.Add("gslsfdebye2", "(", new OneOperator1_< double, double >(gsl_sf_debye_2__));
  Global.Add("gslsfdebye3", "(", new OneOperator1_< double, double >(gsl_sf_debye_3__));
  Global.Add("gslsfdebye4", "(", new OneOperator1_< double, double >(gsl_sf_debye_4__));
  Global.Add("gslsfdebye5", "(", new OneOperator1_< double, double >(gsl_sf_debye_5__));
  Global.Add("gslsfdebye6", "(", new OneOperator1_< double, double >(gsl_sf_debye_6__));
  Global.Add("gslsfdilog", "(", new OneOperator1_< double, double >(gsl_sf_dilog__));
  Global.Add("gslsfmultiply", "(", new OneOperator2_< double, double, double >(gsl_sf_multiply__));
  Global.Add("gslsfellintKcomp", "(",
             new OneOperator2_< double, double, long >(gsl_sf_ellint_Kcomp__));
  Global.Add("gslsfellintEcomp", "(",
             new OneOperator2_< double, double, long >(gsl_sf_ellint_Ecomp__));
  Global.Add("gslsfellintPcomp", "(",
             new OneOperator3_< double, double, double, long >(gsl_sf_ellint_Pcomp__));
  Global.Add("gslsfellintDcomp", "(",
             new OneOperator2_< double, double, long >(gsl_sf_ellint_Dcomp__));
  Global.Add("gslsfellintF", "(",
             new OneOperator3_< double, double, double, long >(gsl_sf_ellint_F__));
  Global.Add("gslsfellintE", "(",
             new OneOperator3_< double, double, double, long >(gsl_sf_ellint_E__));
  Global.Add("gslsfellintRC", "(",
             new OneOperator3_< double, double, double, long >(gsl_sf_ellint_RC__));
  Global.Add("gslsferfc", "(", new OneOperator1_< double, double >(gsl_sf_erfc__));
  Global.Add("gslsflogerfc", "(", new OneOperator1_< double, double >(gsl_sf_log_erfc__));
  Global.Add("gslsferf", "(", new OneOperator1_< double, double >(gsl_sf_erf__));
  Global.Add("gslsferfZ", "(", new OneOperator1_< double, double >(gsl_sf_erf_Z__));
  Global.Add("gslsferfQ", "(", new OneOperator1_< double, double >(gsl_sf_erf_Q__));
  Global.Add("gslsfhazard", "(", new OneOperator1_< double, double >(gsl_sf_hazard__));
  Global.Add("gslsfexp", "(", new OneOperator1_< double, double >(gsl_sf_exp__));
  Global.Add("gslsfexpmult", "(", new OneOperator2_< double, double, double >(gsl_sf_exp_mult__));
  Global.Add("gslsfexpm1", "(", new OneOperator1_< double, double >(gsl_sf_expm1__));
  Global.Add("gslsfexprel", "(", new OneOperator1_< double, double >(gsl_sf_exprel__));
  Global.Add("gslsfexprel2", "(", new OneOperator1_< double, double >(gsl_sf_exprel_2__));
  Global.Add("gslsfexpreln", "(", new OneOperator2_< double, long, double >(gsl_sf_exprel_n__));
  Global.Add("gslsfexpintE1", "(", new OneOperator1_< double, double >(gsl_sf_expint_E1__));
  Global.Add("gslsfexpintE2", "(", new OneOperator1_< double, double >(gsl_sf_expint_E2__));
  Global.Add("gslsfexpintEn", "(", new OneOperator2_< double, long, double >(gsl_sf_expint_En__));
  Global.Add("gslsfexpintE1scaled", "(",
             new OneOperator1_< double, double >(gsl_sf_expint_E1_scaled__));
  Global.Add("gslsfexpintE2scaled", "(",
             new OneOperator1_< double, double >(gsl_sf_expint_E2_scaled__));
  Global.Add("gslsfexpintEnscaled", "(",
             new OneOperator2_< double, long, double >(gsl_sf_expint_En_scaled__));
  Global.Add("gslsfexpintEi", "(", new OneOperator1_< double, double >(gsl_sf_expint_Ei__));
  Global.Add("gslsfexpintEiscaled", "(",
             new OneOperator1_< double, double >(gsl_sf_expint_Ei_scaled__));
  Global.Add("gslsfShi", "(", new OneOperator1_< double, double >(gsl_sf_Shi__));
  Global.Add("gslsfChi", "(", new OneOperator1_< double, double >(gsl_sf_Chi__));
  Global.Add("gslsfexpint3", "(", new OneOperator1_< double, double >(gsl_sf_expint_3__));
  Global.Add("gslsfSi", "(", new OneOperator1_< double, double >(gsl_sf_Si__));
  Global.Add("gslsfCi", "(", new OneOperator1_< double, double >(gsl_sf_Ci__));
  Global.Add("gslsfatanint", "(", new OneOperator1_< double, double >(gsl_sf_atanint__));
  Global.Add("gslsffermidiracm1", "(",
             new OneOperator1_< double, double >(gsl_sf_fermi_dirac_m1__));
  Global.Add("gslsffermidirac0", "(", new OneOperator1_< double, double >(gsl_sf_fermi_dirac_0__));
  Global.Add("gslsffermidirac1", "(", new OneOperator1_< double, double >(gsl_sf_fermi_dirac_1__));
  Global.Add("gslsffermidirac2", "(", new OneOperator1_< double, double >(gsl_sf_fermi_dirac_2__));
  Global.Add("gslsffermidiracint", "(",
             new OneOperator2_< double, long, double >(gsl_sf_fermi_dirac_int__));
  Global.Add("gslsffermidiracmhalf", "(",
             new OneOperator1_< double, double >(gsl_sf_fermi_dirac_mhalf__));
  Global.Add("gslsffermidirachalf", "(",
             new OneOperator1_< double, double >(gsl_sf_fermi_dirac_half__));
  Global.Add("gslsffermidirac3half", "(",
             new OneOperator1_< double, double >(gsl_sf_fermi_dirac_3half__));
  Global.Add("gslsffermidiracinc0", "(",
             new OneOperator2_< double, double, double >(gsl_sf_fermi_dirac_inc_0__));
  Global.Add("gslsflngamma", "(", new OneOperator1_< double, double >(gsl_sf_lngamma__));
  Global.Add("gslsfgamma", "(", new OneOperator1_< double, double >(gsl_sf_gamma__));
  Global.Add("gslsfgammastar", "(", new OneOperator1_< double, double >(gsl_sf_gammastar__));
  Global.Add("gslsfgammainv", "(", new OneOperator1_< double, double >(gsl_sf_gammainv__));
  Global.Add("gslsftaylorcoeff", "(",
             new OneOperator2_< double, long, double >(gsl_sf_taylorcoeff__));
  Global.Add("gslsffact", "(", new OneOperator1_< double, long >(gsl_sf_fact__));
  Global.Add("gslsfdoublefact", "(", new OneOperator1_< double, long >(gsl_sf_doublefact__));
  Global.Add("gslsflnfact", "(", new OneOperator1_< double, long >(gsl_sf_lnfact__));
  Global.Add("gslsflndoublefact", "(", new OneOperator1_< double, long >(gsl_sf_lndoublefact__));
  Global.Add("gslsflnchoose", "(", new OneOperator2_< double, long, long >(gsl_sf_lnchoose__));
  Global.Add("gslsfchoose", "(", new OneOperator2_< double, long, long >(gsl_sf_choose__));
  Global.Add("gslsflnpoch", "(", new OneOperator2_< double, double, double >(gsl_sf_lnpoch__));
  Global.Add("gslsfpoch", "(", new OneOperator2_< double, double, double >(gsl_sf_poch__));
  Global.Add("gslsfpochrel", "(", new OneOperator2_< double, double, double >(gsl_sf_pochrel__));
  Global.Add("gslsfgammaincQ", "(",
             new OneOperator2_< double, double, double >(gsl_sf_gamma_inc_Q__));
  Global.Add("gslsfgammaincP", "(",
             new OneOperator2_< double, double, double >(gsl_sf_gamma_inc_P__));
  Global.Add("gslsfgammainc", "(", new OneOperator2_< double, double, double >(gsl_sf_gamma_inc__));
  Global.Add("gslsflnbeta", "(", new OneOperator2_< double, double, double >(gsl_sf_lnbeta__));
  Global.Add("gslsfbeta", "(", new OneOperator2_< double, double, double >(gsl_sf_beta__));
  Global.Add("gslsfbetainc", "(",
             new OneOperator3_< double, double, double, double >(gsl_sf_beta_inc__));
  Global.Add("gslsfgegenpoly1", "(",
             new OneOperator2_< double, double, double >(gsl_sf_gegenpoly_1__));
  Global.Add("gslsfgegenpoly2", "(",
             new OneOperator2_< double, double, double >(gsl_sf_gegenpoly_2__));
  Global.Add("gslsfgegenpoly3", "(",
             new OneOperator2_< double, double, double >(gsl_sf_gegenpoly_3__));
  Global.Add("gslsfgegenpolyn", "(",
             new OneOperator3_< double, long, double, double >(gsl_sf_gegenpoly_n__));
  Global.Add("gslsfhyperg0F1", "(",
             new OneOperator2_< double, double, double >(gsl_sf_hyperg_0F1__));
  Global.Add("gslsfhyperg1F1int", "(",
             new OneOperator3_< double, long, long, double >(gsl_sf_hyperg_1F1_int__));
  Global.Add("gslsfhyperg1F1", "(",
             new OneOperator3_< double, double, double, double >(gsl_sf_hyperg_1F1__));
  Global.Add("gslsfhypergUint", "(",
             new OneOperator3_< double, long, long, double >(gsl_sf_hyperg_U_int__));
  Global.Add("gslsfhypergU", "(",
             new OneOperator3_< double, double, double, double >(gsl_sf_hyperg_U__));
  Global.Add("gslsfhyperg2F0", "(",
             new OneOperator3_< double, double, double, double >(gsl_sf_hyperg_2F0__));
  Global.Add("gslsflaguerre1", "(",
             new OneOperator2_< double, double, double >(gsl_sf_laguerre_1__));
  Global.Add("gslsflaguerre2", "(",
             new OneOperator2_< double, double, double >(gsl_sf_laguerre_2__));
  Global.Add("gslsflaguerre3", "(",
             new OneOperator2_< double, double, double >(gsl_sf_laguerre_3__));
  Global.Add("gslsflaguerren", "(",
             new OneOperator3_< double, long, double, double >(gsl_sf_laguerre_n__));
  Global.Add("gslsflambertW0", "(", new OneOperator1_< double, double >(gsl_sf_lambert_W0__));
  Global.Add("gslsflambertWm1", "(", new OneOperator1_< double, double >(gsl_sf_lambert_Wm1__));
  Global.Add("gslsflegendrePl", "(",
             new OneOperator2_< double, long, double >(gsl_sf_legendre_Pl__));
  Global.Add("gslsflegendreP1", "(", new OneOperator1_< double, double >(gsl_sf_legendre_P1__));
  Global.Add("gslsflegendreP2", "(", new OneOperator1_< double, double >(gsl_sf_legendre_P2__));
  Global.Add("gslsflegendreP3", "(", new OneOperator1_< double, double >(gsl_sf_legendre_P3__));
  Global.Add("gslsflegendreQ0", "(", new OneOperator1_< double, double >(gsl_sf_legendre_Q0__));
  Global.Add("gslsflegendreQ1", "(", new OneOperator1_< double, double >(gsl_sf_legendre_Q1__));
  Global.Add("gslsflegendreQl", "(",
             new OneOperator2_< double, long, double >(gsl_sf_legendre_Ql__));
  Global.Add("gslsflegendrePlm", "(",
             new OneOperator3_< double, long, long, double >(gsl_sf_legendre_Plm__));
  Global.Add("gslsflegendresphPlm", "(",
             new OneOperator3_< double, long, long, double >(gsl_sf_legendre_sphPlm__));
  // Global.Add("gslsflegendrearraysize","(",new OneOperator2_<long,long,long>(
  // gsl_sf_legendre_array_size__));
  Global.Add("gslsfconicalPhalf", "(",
             new OneOperator2_< double, double, double >(gsl_sf_conicalP_half__));
  Global.Add("gslsfconicalPmhalf", "(",
             new OneOperator2_< double, double, double >(gsl_sf_conicalP_mhalf__));
  Global.Add("gslsfconicalP0", "(",
             new OneOperator2_< double, double, double >(gsl_sf_conicalP_0__));
  Global.Add("gslsfconicalP1", "(",
             new OneOperator2_< double, double, double >(gsl_sf_conicalP_1__));
  Global.Add("gslsfconicalPsphreg", "(",
             new OneOperator3_< double, long, double, double >(gsl_sf_conicalP_sph_reg__));
  Global.Add("gslsfconicalPcylreg", "(",
             new OneOperator3_< double, long, double, double >(gsl_sf_conicalP_cyl_reg__));
  Global.Add("gslsflegendreH3d0", "(",
             new OneOperator2_< double, double, double >(gsl_sf_legendre_H3d_0__));
  Global.Add("gslsflegendreH3d1", "(",
             new OneOperator2_< double, double, double >(gsl_sf_legendre_H3d_1__));
  Global.Add("gslsflegendreH3d", "(",
             new OneOperator3_< double, long, double, double >(gsl_sf_legendre_H3d__));
  Global.Add("gslsflog", "(", new OneOperator1_< double, double >(gsl_sf_log__));
  Global.Add("gslsflogabs", "(", new OneOperator1_< double, double >(gsl_sf_log_abs__));
  Global.Add("gslsflog1plusx", "(", new OneOperator1_< double, double >(gsl_sf_log_1plusx__));
  Global.Add("gslsflog1plusxmx", "(", new OneOperator1_< double, double >(gsl_sf_log_1plusx_mx__));
  Global.Add("gslsfpowint", "(", new OneOperator2_< double, double, long >(gsl_sf_pow_int__));
  Global.Add("gslsfpsiint", "(", new OneOperator1_< double, long >(gsl_sf_psi_int__));
  Global.Add("gslsfpsi", "(", new OneOperator1_< double, double >(gsl_sf_psi__));
  Global.Add("gslsfpsi1piy", "(", new OneOperator1_< double, double >(gsl_sf_psi_1piy__));
  Global.Add("gslsfpsi1int", "(", new OneOperator1_< double, long >(gsl_sf_psi_1_int__));
  Global.Add("gslsfpsi1", "(", new OneOperator1_< double, double >(gsl_sf_psi_1__));
  Global.Add("gslsfpsin", "(", new OneOperator2_< double, long, double >(gsl_sf_psi_n__));
  Global.Add("gslsfsynchrotron1", "(", new OneOperator1_< double, double >(gsl_sf_synchrotron_1__));
  Global.Add("gslsfsynchrotron2", "(", new OneOperator1_< double, double >(gsl_sf_synchrotron_2__));
  Global.Add("gslsftransport2", "(", new OneOperator1_< double, double >(gsl_sf_transport_2__));
  Global.Add("gslsftransport3", "(", new OneOperator1_< double, double >(gsl_sf_transport_3__));
  Global.Add("gslsftransport4", "(", new OneOperator1_< double, double >(gsl_sf_transport_4__));
  Global.Add("gslsftransport5", "(", new OneOperator1_< double, double >(gsl_sf_transport_5__));
  Global.Add("gslsfsin", "(", new OneOperator1_< double, double >(gsl_sf_sin__));
  Global.Add("gslsfcos", "(", new OneOperator1_< double, double >(gsl_sf_cos__));
  Global.Add("gslsfhypot", "(", new OneOperator2_< double, double, double >(gsl_sf_hypot__));
  Global.Add("gslsfsinc", "(", new OneOperator1_< double, double >(gsl_sf_sinc__));
  Global.Add("gslsflnsinh", "(", new OneOperator1_< double, double >(gsl_sf_lnsinh__));
  Global.Add("gslsflncosh", "(", new OneOperator1_< double, double >(gsl_sf_lncosh__));
  Global.Add("gslsfanglerestrictsymm", "(",
             new OneOperator1_< double, double >(gsl_sf_angle_restrict_symm__));
  Global.Add("gslsfanglerestrictpos", "(",
             new OneOperator1_< double, double >(gsl_sf_angle_restrict_pos__));
  Global.Add("gslsfzetaint", "(", new OneOperator1_< double, long >(gsl_sf_zeta_int__));
  Global.Add("gslsfzeta", "(", new OneOperator1_< double, double >(gsl_sf_zeta__));
  Global.Add("gslsfzetam1", "(", new OneOperator1_< double, double >(gsl_sf_zetam1__));
  Global.Add("gslsfzetam1int", "(", new OneOperator1_< double, long >(gsl_sf_zetam1_int__));
  Global.Add("gslsfhzeta", "(", new OneOperator2_< double, double, double >(gsl_sf_hzeta__));
  Global.Add("gslsfetaint", "(", new OneOperator1_< double, long >(gsl_sf_eta_int__));
  Global.Add("gslsfeta", "(", new OneOperator1_< double, double >(gsl_sf_eta__));
}

/*****************/
/*****************/
