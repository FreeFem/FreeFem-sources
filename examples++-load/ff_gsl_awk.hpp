 /*  
 minssing gsl_sf_coulomb.h:double gsl_sf_hydrogenicR(const int n, const int l, const double Z, const double r);
 minssing gsl_sf_coupling.h:double gsl_sf_coupling_3j(int two_ja, int two_jb, int two_jc,
 minssing gsl_sf_coupling.h:double gsl_sf_coupling_6j(int two_ja, int two_jb, int two_jc,
 minssing gsl_sf_coupling.h:double gsl_sf_coupling_RacahW(int two_ja, int two_jb, int two_jc,
 minssing gsl_sf_coupling.h:double gsl_sf_coupling_9j(int two_ja, int two_jb, int two_jc,
 minssing gsl_sf_coupling.h:double gsl_sf_coupling_6j_INCORRECT(int two_ja, int two_jb, int two_jc,
 minssing gsl_sf_ellint.h:double gsl_sf_ellint_Pcomp(double k, double n, gsl_mode_t mode);
 minssing gsl_sf_ellint.h:double gsl_sf_ellint_F(double phi, double k, gsl_mode_t mode);
 minssing gsl_sf_ellint.h:double gsl_sf_ellint_E(double phi, double k, gsl_mode_t mode);
 minssing gsl_sf_ellint.h:double gsl_sf_ellint_P(double phi, double k, double n, gsl_mode_t mode);
 minssing gsl_sf_ellint.h:double gsl_sf_ellint_D(double phi, double k, double n, gsl_mode_t mode);
 minssing gsl_sf_ellint.h:double gsl_sf_ellint_RC(double x, double y, gsl_mode_t mode);
 minssing gsl_sf_ellint.h:double gsl_sf_ellint_RD(double x, double y, double z, gsl_mode_t mode);
 minssing gsl_sf_ellint.h:double gsl_sf_ellint_RF(double x, double y, double z, gsl_mode_t mode);
 minssing gsl_sf_ellint.h:double gsl_sf_ellint_RJ(double x, double y, double z, double p, gsl_mode_t mode);
 minssing gsl_sf_gamma.h:double gsl_sf_beta_inc(const double a, const double b, const double x);
 minssing gsl_sf_gegenbauer.h:double gsl_sf_gegenpoly_n(int n, double lambda, double x);
 minssing gsl_sf_hyperg.h:double gsl_sf_hyperg_1F1_int(const int m, const int n, double x);
 minssing gsl_sf_hyperg.h:double gsl_sf_hyperg_1F1(double a, double b, double x);
 minssing gsl_sf_hyperg.h:double gsl_sf_hyperg_U_int(const int m, const int n, const double x);
 minssing gsl_sf_hyperg.h:double gsl_sf_hyperg_U(const double a, const double b, const double x);
 minssing gsl_sf_hyperg.h:double gsl_sf_hyperg_2F1(double a, double b, double c, double x);
 minssing gsl_sf_hyperg.h:double gsl_sf_hyperg_2F1_conj(double aR, double aI, double c, double x);
 minssing gsl_sf_hyperg.h:double gsl_sf_hyperg_2F1_renorm(double a, double b, double c, double x);
 minssing gsl_sf_hyperg.h:double gsl_sf_hyperg_2F1_conj_renorm(double aR, double aI, double c, double x);
 minssing gsl_sf_hyperg.h:double     gsl_sf_hyperg_2F0(const double a, const double b, const double x);
 minssing gsl_sf_laguerre.h:double     gsl_sf_laguerre_n(int n, double a, double x);
 minssing gsl_sf_legendre.h:double  gsl_sf_legendre_Plm(const int l, const int m, const double x);
 minssing gsl_sf_legendre.h:double  gsl_sf_legendre_sphPlm(const int l, const int m, const double x);
 minssing gsl_sf_legendre.h:double gsl_sf_conicalP_sph_reg(const int l, const double lambda, const double x);
 minssing gsl_sf_legendre.h:double gsl_sf_conicalP_cyl_reg(const int m, const double lambda, const double x);
 minssing gsl_sf_legendre.h:double gsl_sf_legendre_H3d(const int l, const double lambda, const double eta);
 */ 
/*****************/
/*****************/
double gsl_sf_airy_Ai( double x, long y   ) { return gsl_sf_airy_Ai( ( const double) x, ( gsl_mode_t)  y  );}
double gsl_sf_airy_Bi( double x, long y   ) { return gsl_sf_airy_Bi( ( const double) x, ( gsl_mode_t)  y  );}
double gsl_sf_airy_Ai_scaled( double x, long y   ) { return gsl_sf_airy_Ai_scaled( ( const double) x, ( gsl_mode_t)  y  );}
double gsl_sf_airy_Bi_scaled( double x, long y   ) { return gsl_sf_airy_Bi_scaled( ( const double) x, ( gsl_mode_t)  y  );}
double gsl_sf_airy_Ai_deriv( double x, long y   ) { return gsl_sf_airy_Ai_deriv( ( const double) x, ( gsl_mode_t)  y  );}
double gsl_sf_airy_Bi_deriv( double x, long y   ) { return gsl_sf_airy_Bi_deriv( ( const double) x, ( gsl_mode_t)  y  );}
double gsl_sf_airy_Ai_deriv_scaled( double x, long y   ) { return gsl_sf_airy_Ai_deriv_scaled( ( const double) x, ( gsl_mode_t)  y  );}
double gsl_sf_airy_Bi_deriv_scaled( double x, long y   ) { return gsl_sf_airy_Bi_deriv_scaled( ( const double) x, ( gsl_mode_t)  y  );}
double gsl_sf_airy_zero_Ai( long x ) { return gsl_sf_airy_zero_Ai( ( unsigned int) x );}
double gsl_sf_airy_zero_Bi( long x ) { return gsl_sf_airy_zero_Bi( ( unsigned int) x );}
double gsl_sf_airy_zero_Ai_deriv( long x ) { return gsl_sf_airy_zero_Ai_deriv( ( unsigned int) x );}
double gsl_sf_airy_zero_Bi_deriv( long x ) { return gsl_sf_airy_zero_Bi_deriv( ( unsigned int) x );}
double gsl_sf_bessel_Jn( long x, double y   ) { return gsl_sf_bessel_Jn( ( const int) x, ( const double)  y  );}
double gsl_sf_bessel_Yn( long x, double y   ) { return gsl_sf_bessel_Yn( ( const int) x, ( const double)  y  );}
double gsl_sf_bessel_In( long x, double y   ) { return gsl_sf_bessel_In( ( const int) x, ( const double)  y  );}
double gsl_sf_bessel_In_scaled( long x, double y   ) { return gsl_sf_bessel_In_scaled( ( const int) x, ( const double)  y  );}
double gsl_sf_bessel_Kn( long x, double y   ) { return gsl_sf_bessel_Kn( ( const int) x, ( const double)  y  );}
double gsl_sf_bessel_Kn_scaled( long x, double y   ) { return gsl_sf_bessel_Kn_scaled( ( const int) x, ( const double)  y  );}
double gsl_sf_bessel_jl( long x, double y   ) { return gsl_sf_bessel_jl( ( const int) x, ( const double)  y  );}
double gsl_sf_bessel_yl( long x, double y   ) { return gsl_sf_bessel_yl( ( const int) x, ( const double)  y  );}
double gsl_sf_bessel_il_scaled( long x, double y   ) { return gsl_sf_bessel_il_scaled( ( const int) x, ( const double)  y  );}
double gsl_sf_bessel_kl_scaled( long x, double y   ) { return gsl_sf_bessel_kl_scaled( ( const int) x, ( const double)  y  );}
double gsl_sf_bessel_Jnu( double x, long y   ) { return gsl_sf_bessel_Jnu( ( const double) x, ( const double)  y  );}
double gsl_sf_bessel_Ynu( double x, long y   ) { return gsl_sf_bessel_Ynu( ( const double) x, ( const double)  y  );}
double gsl_sf_bessel_Inu_scaled( double x, long y   ) { return gsl_sf_bessel_Inu_scaled( ( double) x, ( double)  y  );}
double gsl_sf_bessel_Inu( double x, long y   ) { return gsl_sf_bessel_Inu( ( double) x, ( double)  y  );}
double gsl_sf_bessel_Knu_scaled( double x, long y   ) { return gsl_sf_bessel_Knu_scaled( ( const double) x, ( const double)  y  );}
double gsl_sf_bessel_Knu( double x, long y   ) { return gsl_sf_bessel_Knu( ( const double) x, ( const double)  y  );}
double gsl_sf_bessel_lnKnu( double x, long y   ) { return gsl_sf_bessel_lnKnu( ( const double) x, ( const double)  y  );}
double gsl_sf_bessel_zero_J0( long x ) { return gsl_sf_bessel_zero_J0( ( unsigned int) x );}
double gsl_sf_bessel_zero_J1( long x ) { return gsl_sf_bessel_zero_J1( ( unsigned int) x );}
double gsl_sf_bessel_zero_Jnu( double x, long y   ) { return gsl_sf_bessel_zero_Jnu( ( double) x, ( unsigned int)  y  );}
double gsl_sf_hydrogenicR_1( double x, long y   ) { return gsl_sf_hydrogenicR_1( ( const double) x, ( const double)  y  );}
double gsl_sf_multiply( double x, long y   ) { return gsl_sf_multiply( ( const double) x, ( const double)  y  );}
double gsl_sf_ellint_Kcomp( double x, long y   ) { return gsl_sf_ellint_Kcomp( ( double) x, ( gsl_mode_t)  y  );}
double gsl_sf_ellint_Ecomp( double x, long y   ) { return gsl_sf_ellint_Ecomp( ( double) x, ( gsl_mode_t)  y  );}
double gsl_sf_ellint_Dcomp( double x, long y   ) { return gsl_sf_ellint_Dcomp( ( double) x, ( gsl_mode_t)  y  );}
double gsl_sf_exp_mult( double x, long y   ) { return gsl_sf_exp_mult( ( const double) x, ( const double)  y  );}
double gsl_sf_exprel_n( long x, double y   ) { return gsl_sf_exprel_n( ( const int) x, ( const double)  y  );}
double gsl_sf_expint_En( long x, double y   ) { return gsl_sf_expint_En( ( const int) x, ( const double)  y  );}
double gsl_sf_expint_En_scaled( long x, double y   ) { return gsl_sf_expint_En_scaled( ( const int) x, ( const double)  y  );}
double gsl_sf_fermi_dirac_int( long x, double y   ) { return gsl_sf_fermi_dirac_int( ( const int) x, ( const double)  y  );}
double gsl_sf_fermi_dirac_inc_0( double x, long y   ) { return gsl_sf_fermi_dirac_inc_0( ( const double) x, ( const double)  y  );}
double gsl_sf_taylorcoeff( long x, double y   ) { return gsl_sf_taylorcoeff( ( const int) x, ( const double)  y  );}
double gsl_sf_fact( long x ) { return gsl_sf_fact( ( const unsigned int) x );}
double gsl_sf_doublefact( long x ) { return gsl_sf_doublefact( ( const unsigned int) x );}
double gsl_sf_lnfact( long x ) { return gsl_sf_lnfact( ( const unsigned int) x );}
double gsl_sf_lndoublefact( long x ) { return gsl_sf_lndoublefact( ( const unsigned int) x );}
double gsl_sf_lnchoose( long x, long y   ) { return gsl_sf_lnchoose( ( unsigned int) x, ( unsigned int)  y  );}
double gsl_sf_choose( long x, long y   ) { return gsl_sf_choose( ( unsigned int) x, ( unsigned int)  y  );}
double gsl_sf_lnpoch( double x, long y   ) { return gsl_sf_lnpoch( ( const double) x, ( const double)  y  );}
double gsl_sf_poch( double x, long y   ) { return gsl_sf_poch( ( const double) x, ( const double)  y  );}
double gsl_sf_pochrel( double x, long y   ) { return gsl_sf_pochrel( ( const double) x, ( const double)  y  );}
double gsl_sf_gamma_inc_Q( double x, long y   ) { return gsl_sf_gamma_inc_Q( ( const double) x, ( const double)  y  );}
double gsl_sf_gamma_inc_P( double x, long y   ) { return gsl_sf_gamma_inc_P( ( const double) x, ( const double)  y  );}
double gsl_sf_gamma_inc( double x, long y   ) { return gsl_sf_gamma_inc( ( const double) x, ( const double)  y  );}
double gsl_sf_lnbeta( double x, long y   ) { return gsl_sf_lnbeta( ( const double) x, ( const double)  y  );}
double gsl_sf_beta( double x, long y   ) { return gsl_sf_beta( ( const double) x, ( const double)  y  );}
double gsl_sf_gegenpoly_1( double x, long y   ) { return gsl_sf_gegenpoly_1( ( double) x, ( double)  y  );}
double gsl_sf_gegenpoly_2( double x, long y   ) { return gsl_sf_gegenpoly_2( ( double) x, ( double)  y  );}
double gsl_sf_gegenpoly_3( double x, long y   ) { return gsl_sf_gegenpoly_3( ( double) x, ( double)  y  );}
double gsl_sf_hyperg_0F1( double x, long y   ) { return gsl_sf_hyperg_0F1( ( const double) x, ( const double)  y  );}
double gsl_sf_laguerre_1( double x, long y   ) { return gsl_sf_laguerre_1( ( double) x, ( double)  y  );}
double gsl_sf_laguerre_2( double x, long y   ) { return gsl_sf_laguerre_2( ( double) x, ( double)  y  );}
double gsl_sf_laguerre_3( double x, long y   ) { return gsl_sf_laguerre_3( ( double) x, ( double)  y  );}
double gsl_sf_legendre_Pl( long x, double y   ) { return gsl_sf_legendre_Pl( ( const int) x, ( const double)  y  );}
double gsl_sf_legendre_Ql( long x, double y   ) { return gsl_sf_legendre_Ql( ( const int) x, ( const double)  y  );}
double gsl_sf_conicalP_half( double x, long y   ) { return gsl_sf_conicalP_half( ( const double) x, ( const double)  y  );}
double gsl_sf_conicalP_mhalf( double x, long y   ) { return gsl_sf_conicalP_mhalf( ( const double) x, ( const double)  y  );}
double gsl_sf_conicalP_0( double x, long y   ) { return gsl_sf_conicalP_0( ( const double) x, ( const double)  y  );}
double gsl_sf_conicalP_1( double x, long y   ) { return gsl_sf_conicalP_1( ( const double) x, ( const double)  y  );}
double gsl_sf_legendre_H3d_0( double x, long y   ) { return gsl_sf_legendre_H3d_0( ( const double) x, ( const double)  y  );}
double gsl_sf_legendre_H3d_1( double x, long y   ) { return gsl_sf_legendre_H3d_1( ( const double) x, ( const double)  y  );}
double gsl_sf_pow_int( double x, long y   ) { return gsl_sf_pow_int( ( const double) x, ( const int)  y  );}
double gsl_sf_psi_int( long x ) { return gsl_sf_psi_int( ( const int) x );}
double gsl_sf_psi_1_int( long x ) { return gsl_sf_psi_1_int( ( const int) x );}
double gsl_sf_psi_n( long x, double y   ) { return gsl_sf_psi_n( ( const int) x, ( const double)  y  );}
double gsl_sf_hypot( double x, long y   ) { return gsl_sf_hypot( ( const double) x, ( const double)  y  );}
double gsl_sf_zeta_int( long x ) { return gsl_sf_zeta_int( ( const int) x );}
double gsl_sf_zetam1_int( long x ) { return gsl_sf_zetam1_int( ( const int) x );}
double gsl_sf_hzeta( double x, long y   ) { return gsl_sf_hzeta( ( const double) x, ( const double)  y  );}
double gsl_sf_eta_int( long x ) { return gsl_sf_eta_int( ( const int) x );}

/*****************/
/*****************/
 void init_gsl_sf() { 


   Global.Add("gslsfairyAi","(",new OneOperator2<double,double,long >( gsl_sf_airy_Ai )); 
   Global.Add("gslsfairyBi","(",new OneOperator2<double,double,long >( gsl_sf_airy_Bi )); 
   Global.Add("gslsfairyAiscaled","(",new OneOperator2<double,double,long >( gsl_sf_airy_Ai_scaled )); 
   Global.Add("gslsfairyBiscaled","(",new OneOperator2<double,double,long >( gsl_sf_airy_Bi_scaled )); 
   Global.Add("gslsfairyAideriv","(",new OneOperator2<double,double,long >( gsl_sf_airy_Ai_deriv )); 
   Global.Add("gslsfairyBideriv","(",new OneOperator2<double,double,long >( gsl_sf_airy_Bi_deriv )); 
   Global.Add("gslsfairyAiderivscaled","(",new OneOperator2<double,double,long >( gsl_sf_airy_Ai_deriv_scaled )); 
   Global.Add("gslsfairyBiderivscaled","(",new OneOperator2<double,double,long >( gsl_sf_airy_Bi_deriv_scaled )); 
   Global.Add("gslsfairyzeroAi","(",new OneOperator1<double ,long>( gsl_sf_airy_zero_Ai )); 
   Global.Add("gslsfairyzeroBi","(",new OneOperator1<double ,long>( gsl_sf_airy_zero_Bi )); 
   Global.Add("gslsfairyzeroAideriv","(",new OneOperator1<double ,long>( gsl_sf_airy_zero_Ai_deriv )); 
   Global.Add("gslsfairyzeroBideriv","(",new OneOperator1<double ,long>( gsl_sf_airy_zero_Bi_deriv )); 
   Global.Add("gslsfbesselJ0","(",new OneOperator1<double ,double>( gsl_sf_bessel_J0 )); 
   Global.Add("gslsfbesselJ1","(",new OneOperator1<double ,double>( gsl_sf_bessel_J1 )); 
   Global.Add("gslsfbesselJn","(",new OneOperator2<double,long,double >( gsl_sf_bessel_Jn )); 
   Global.Add("gslsfbesselY0","(",new OneOperator1<double ,double>( gsl_sf_bessel_Y0 )); 
   Global.Add("gslsfbesselY1","(",new OneOperator1<double ,double>( gsl_sf_bessel_Y1 )); 
   Global.Add("gslsfbesselYn","(",new OneOperator2<double,long,double >( gsl_sf_bessel_Yn )); 
   Global.Add("gslsfbesselI0","(",new OneOperator1<double ,double>( gsl_sf_bessel_I0 )); 
   Global.Add("gslsfbesselI1","(",new OneOperator1<double ,double>( gsl_sf_bessel_I1 )); 
   Global.Add("gslsfbesselIn","(",new OneOperator2<double,long,double >( gsl_sf_bessel_In )); 
   Global.Add("gslsfbesselI0scaled","(",new OneOperator1<double ,double>( gsl_sf_bessel_I0_scaled )); 
   Global.Add("gslsfbesselI1scaled","(",new OneOperator1<double ,double>( gsl_sf_bessel_I1_scaled )); 
   Global.Add("gslsfbesselInscaled","(",new OneOperator2<double,long,double >( gsl_sf_bessel_In_scaled )); 
   Global.Add("gslsfbesselK0","(",new OneOperator1<double ,double>( gsl_sf_bessel_K0 )); 
   Global.Add("gslsfbesselK1","(",new OneOperator1<double ,double>( gsl_sf_bessel_K1 )); 
   Global.Add("gslsfbesselKn","(",new OneOperator2<double,long,double >( gsl_sf_bessel_Kn )); 
   Global.Add("gslsfbesselK0scaled","(",new OneOperator1<double ,double>( gsl_sf_bessel_K0_scaled )); 
   Global.Add("gslsfbesselK1scaled","(",new OneOperator1<double ,double>( gsl_sf_bessel_K1_scaled )); 
   Global.Add("gslsfbesselKnscaled","(",new OneOperator2<double,long,double >( gsl_sf_bessel_Kn_scaled )); 
   Global.Add("gslsfbesselj0","(",new OneOperator1<double ,double>( gsl_sf_bessel_j0 )); 
   Global.Add("gslsfbesselj1","(",new OneOperator1<double ,double>( gsl_sf_bessel_j1 )); 
   Global.Add("gslsfbesselj2","(",new OneOperator1<double ,double>( gsl_sf_bessel_j2 )); 
   Global.Add("gslsfbesseljl","(",new OneOperator2<double,long,double >( gsl_sf_bessel_jl )); 
   Global.Add("gslsfbessely0","(",new OneOperator1<double ,double>( gsl_sf_bessel_y0 )); 
   Global.Add("gslsfbessely1","(",new OneOperator1<double ,double>( gsl_sf_bessel_y1 )); 
   Global.Add("gslsfbessely2","(",new OneOperator1<double ,double>( gsl_sf_bessel_y2 )); 
   Global.Add("gslsfbesselyl","(",new OneOperator2<double,long,double >( gsl_sf_bessel_yl )); 
   Global.Add("gslsfbesseli0scaled","(",new OneOperator1<double ,double>( gsl_sf_bessel_i0_scaled )); 
   Global.Add("gslsfbesseli1scaled","(",new OneOperator1<double ,double>( gsl_sf_bessel_i1_scaled )); 
   Global.Add("gslsfbesseli2scaled","(",new OneOperator1<double ,double>( gsl_sf_bessel_i2_scaled )); 
   Global.Add("gslsfbesselilscaled","(",new OneOperator2<double,long,double >( gsl_sf_bessel_il_scaled )); 
   Global.Add("gslsfbesselk0scaled","(",new OneOperator1<double ,double>( gsl_sf_bessel_k0_scaled )); 
   Global.Add("gslsfbesselk1scaled","(",new OneOperator1<double ,double>( gsl_sf_bessel_k1_scaled )); 
   Global.Add("gslsfbesselk2scaled","(",new OneOperator1<double ,double>( gsl_sf_bessel_k2_scaled )); 
   Global.Add("gslsfbesselklscaled","(",new OneOperator2<double,long,double >( gsl_sf_bessel_kl_scaled )); 
   Global.Add("gslsfbesselJnu","(",new OneOperator2<double,double,long >( gsl_sf_bessel_Jnu )); 
   Global.Add("gslsfbesselYnu","(",new OneOperator2<double,double,long >( gsl_sf_bessel_Ynu )); 
   Global.Add("gslsfbesselInuscaled","(",new OneOperator2<double,double,long >( gsl_sf_bessel_Inu_scaled )); 
   Global.Add("gslsfbesselInu","(",new OneOperator2<double,double,long >( gsl_sf_bessel_Inu )); 
   Global.Add("gslsfbesselKnuscaled","(",new OneOperator2<double,double,long >( gsl_sf_bessel_Knu_scaled )); 
   Global.Add("gslsfbesselKnu","(",new OneOperator2<double,double,long >( gsl_sf_bessel_Knu )); 
   Global.Add("gslsfbessellnKnu","(",new OneOperator2<double,double,long >( gsl_sf_bessel_lnKnu )); 
   Global.Add("gslsfbesselzeroJ0","(",new OneOperator1<double ,long>( gsl_sf_bessel_zero_J0 )); 
   Global.Add("gslsfbesselzeroJ1","(",new OneOperator1<double ,long>( gsl_sf_bessel_zero_J1 )); 
   Global.Add("gslsfbesselzeroJnu","(",new OneOperator2<double,double,long >( gsl_sf_bessel_zero_Jnu )); 
   Global.Add("gslsfclausen","(",new OneOperator1<double ,double>( gsl_sf_clausen )); 
   Global.Add("gslsfhydrogenicR1","(",new OneOperator2<double,double,long >( gsl_sf_hydrogenicR_1 )); 
   Global.Add("gslsfdawson","(",new OneOperator1<double ,double>( gsl_sf_dawson )); 
   Global.Add("gslsfdebye1","(",new OneOperator1<double ,double>( gsl_sf_debye_1 )); 
   Global.Add("gslsfdebye2","(",new OneOperator1<double ,double>( gsl_sf_debye_2 )); 
   Global.Add("gslsfdebye3","(",new OneOperator1<double ,double>( gsl_sf_debye_3 )); 
   Global.Add("gslsfdebye4","(",new OneOperator1<double ,double>( gsl_sf_debye_4 )); 
   Global.Add("gslsfdebye5","(",new OneOperator1<double ,double>( gsl_sf_debye_5 )); 
   Global.Add("gslsfdebye6","(",new OneOperator1<double ,double>( gsl_sf_debye_6 )); 
   Global.Add("gslsfdilog","(",new OneOperator1<double ,double>( gsl_sf_dilog )); 
   Global.Add("gslsfmultiply","(",new OneOperator2<double,double,long >( gsl_sf_multiply )); 
   Global.Add("gslsfellintKcomp","(",new OneOperator2<double,double,long >( gsl_sf_ellint_Kcomp )); 
   Global.Add("gslsfellintEcomp","(",new OneOperator2<double,double,long >( gsl_sf_ellint_Ecomp )); 
   Global.Add("gslsfellintDcomp","(",new OneOperator2<double,double,long >( gsl_sf_ellint_Dcomp )); 
   Global.Add("gslsferfc","(",new OneOperator1<double ,double>( gsl_sf_erfc )); 
   Global.Add("gslsflogerfc","(",new OneOperator1<double ,double>( gsl_sf_log_erfc )); 
   Global.Add("gslsferf","(",new OneOperator1<double ,double>( gsl_sf_erf )); 
   Global.Add("gslsferfZ","(",new OneOperator1<double ,double>( gsl_sf_erf_Z )); 
   Global.Add("gslsferfQ","(",new OneOperator1<double ,double>( gsl_sf_erf_Q )); 
   Global.Add("gslsfhazard","(",new OneOperator1<double ,double>( gsl_sf_hazard )); 
   Global.Add("gslsfexp","(",new OneOperator1<double ,double>( gsl_sf_exp )); 
   Global.Add("gslsfexpmult","(",new OneOperator2<double,double,long >( gsl_sf_exp_mult )); 
   Global.Add("gslsfexpm1","(",new OneOperator1<double ,double>( gsl_sf_expm1 )); 
   Global.Add("gslsfexprel","(",new OneOperator1<double ,double>( gsl_sf_exprel )); 
   Global.Add("gslsfexprel2","(",new OneOperator1<double ,double>( gsl_sf_exprel_2 )); 
   Global.Add("gslsfexpreln","(",new OneOperator2<double,long,double >( gsl_sf_exprel_n )); 
   Global.Add("gslsfexpintE1","(",new OneOperator1<double ,double>( gsl_sf_expint_E1 )); 
   Global.Add("gslsfexpintE2","(",new OneOperator1<double ,double>( gsl_sf_expint_E2 )); 
   Global.Add("gslsfexpintEn","(",new OneOperator2<double,long,double >( gsl_sf_expint_En )); 
   Global.Add("gslsfexpintE1scaled","(",new OneOperator1<double ,double>( gsl_sf_expint_E1_scaled )); 
   Global.Add("gslsfexpintE2scaled","(",new OneOperator1<double ,double>( gsl_sf_expint_E2_scaled )); 
   Global.Add("gslsfexpintEnscaled","(",new OneOperator2<double,long,double >( gsl_sf_expint_En_scaled )); 
   Global.Add("gslsfexpintEi","(",new OneOperator1<double ,double>( gsl_sf_expint_Ei )); 
   Global.Add("gslsfexpintEiscaled","(",new OneOperator1<double ,double>( gsl_sf_expint_Ei_scaled )); 
   Global.Add("gslsfShi","(",new OneOperator1<double ,double>( gsl_sf_Shi )); 
   Global.Add("gslsfChi","(",new OneOperator1<double ,double>( gsl_sf_Chi )); 
   Global.Add("gslsfexpint3","(",new OneOperator1<double ,double>( gsl_sf_expint_3 )); 
   Global.Add("gslsfSi","(",new OneOperator1<double ,double>( gsl_sf_Si )); 
   Global.Add("gslsfCi","(",new OneOperator1<double ,double>( gsl_sf_Ci )); 
   Global.Add("gslsfatanint","(",new OneOperator1<double ,double>( gsl_sf_atanint )); 
   Global.Add("gslsffermidiracm1","(",new OneOperator1<double ,double>( gsl_sf_fermi_dirac_m1 )); 
   Global.Add("gslsffermidirac0","(",new OneOperator1<double ,double>( gsl_sf_fermi_dirac_0 )); 
   Global.Add("gslsffermidirac1","(",new OneOperator1<double ,double>( gsl_sf_fermi_dirac_1 )); 
   Global.Add("gslsffermidirac2","(",new OneOperator1<double ,double>( gsl_sf_fermi_dirac_2 )); 
   Global.Add("gslsffermidiracint","(",new OneOperator2<double,long,double >( gsl_sf_fermi_dirac_int )); 
   Global.Add("gslsffermidiracmhalf","(",new OneOperator1<double ,double>( gsl_sf_fermi_dirac_mhalf )); 
   Global.Add("gslsffermidirachalf","(",new OneOperator1<double ,double>( gsl_sf_fermi_dirac_half )); 
   Global.Add("gslsffermidirac3half","(",new OneOperator1<double ,double>( gsl_sf_fermi_dirac_3half )); 
   Global.Add("gslsffermidiracinc0","(",new OneOperator2<double,double,long >( gsl_sf_fermi_dirac_inc_0 )); 
   Global.Add("gslsflngamma","(",new OneOperator1<double ,double>( gsl_sf_lngamma )); 
   Global.Add("gslsfgamma","(",new OneOperator1<double ,double>( gsl_sf_gamma )); 
   Global.Add("gslsfgammastar","(",new OneOperator1<double ,double>( gsl_sf_gammastar )); 
   Global.Add("gslsfgammainv","(",new OneOperator1<double ,double>( gsl_sf_gammainv )); 
   Global.Add("gslsftaylorcoeff","(",new OneOperator2<double,long,double >( gsl_sf_taylorcoeff )); 
   Global.Add("gslsffact","(",new OneOperator1<double ,long>( gsl_sf_fact )); 
   Global.Add("gslsfdoublefact","(",new OneOperator1<double ,long>( gsl_sf_doublefact )); 
   Global.Add("gslsflnfact","(",new OneOperator1<double ,long>( gsl_sf_lnfact )); 
   Global.Add("gslsflndoublefact","(",new OneOperator1<double ,long>( gsl_sf_lndoublefact )); 
   Global.Add("gslsflnchoose","(",new OneOperator2<double,long,long >( gsl_sf_lnchoose )); 
   Global.Add("gslsfchoose","(",new OneOperator2<double,long,long >( gsl_sf_choose )); 
   Global.Add("gslsflnpoch","(",new OneOperator2<double,double,long >( gsl_sf_lnpoch )); 
   Global.Add("gslsfpoch","(",new OneOperator2<double,double,long >( gsl_sf_poch )); 
   Global.Add("gslsfpochrel","(",new OneOperator2<double,double,long >( gsl_sf_pochrel )); 
   Global.Add("gslsfgammaincQ","(",new OneOperator2<double,double,long >( gsl_sf_gamma_inc_Q )); 
   Global.Add("gslsfgammaincP","(",new OneOperator2<double,double,long >( gsl_sf_gamma_inc_P )); 
   Global.Add("gslsfgammainc","(",new OneOperator2<double,double,long >( gsl_sf_gamma_inc )); 
   Global.Add("gslsflnbeta","(",new OneOperator2<double,double,long >( gsl_sf_lnbeta )); 
   Global.Add("gslsfbeta","(",new OneOperator2<double,double,long >( gsl_sf_beta )); 
   Global.Add("gslsfgegenpoly1","(",new OneOperator2<double,double,long >( gsl_sf_gegenpoly_1 )); 
   Global.Add("gslsfgegenpoly2","(",new OneOperator2<double,double,long >( gsl_sf_gegenpoly_2 )); 
   Global.Add("gslsfgegenpoly3","(",new OneOperator2<double,double,long >( gsl_sf_gegenpoly_3 )); 
   Global.Add("gslsfhyperg0F1","(",new OneOperator2<double,double,long >( gsl_sf_hyperg_0F1 )); 
   Global.Add("gslsflaguerre1","(",new OneOperator2<double,double,long >( gsl_sf_laguerre_1 )); 
   Global.Add("gslsflaguerre2","(",new OneOperator2<double,double,long >( gsl_sf_laguerre_2 )); 
   Global.Add("gslsflaguerre3","(",new OneOperator2<double,double,long >( gsl_sf_laguerre_3 )); 
   Global.Add("gslsflambertW0","(",new OneOperator1<double ,double>( gsl_sf_lambert_W0 )); 
   Global.Add("gslsflambertWm1","(",new OneOperator1<double ,double>( gsl_sf_lambert_Wm1 )); 
   Global.Add("gslsflegendrePl","(",new OneOperator2<double,long,double >( gsl_sf_legendre_Pl )); 
   Global.Add("gslsflegendreP1","(",new OneOperator1<double ,double>( gsl_sf_legendre_P1 )); 
   Global.Add("gslsflegendreP2","(",new OneOperator1<double ,double>( gsl_sf_legendre_P2 )); 
   Global.Add("gslsflegendreP3","(",new OneOperator1<double ,double>( gsl_sf_legendre_P3 )); 
   Global.Add("gslsflegendreQ0","(",new OneOperator1<double ,double>( gsl_sf_legendre_Q0 )); 
   Global.Add("gslsflegendreQ1","(",new OneOperator1<double ,double>( gsl_sf_legendre_Q1 )); 
   Global.Add("gslsflegendreQl","(",new OneOperator2<double,long,double >( gsl_sf_legendre_Ql )); 
   Global.Add("gslsfconicalPhalf","(",new OneOperator2<double,double,long >( gsl_sf_conicalP_half )); 
   Global.Add("gslsfconicalPmhalf","(",new OneOperator2<double,double,long >( gsl_sf_conicalP_mhalf )); 
   Global.Add("gslsfconicalP0","(",new OneOperator2<double,double,long >( gsl_sf_conicalP_0 )); 
   Global.Add("gslsfconicalP1","(",new OneOperator2<double,double,long >( gsl_sf_conicalP_1 )); 
   Global.Add("gslsflegendreH3d0","(",new OneOperator2<double,double,long >( gsl_sf_legendre_H3d_0 )); 
   Global.Add("gslsflegendreH3d1","(",new OneOperator2<double,double,long >( gsl_sf_legendre_H3d_1 )); 
   Global.Add("gslsflog","(",new OneOperator1<double ,double>( gsl_sf_log )); 
   Global.Add("gslsflogabs","(",new OneOperator1<double ,double>( gsl_sf_log_abs )); 
   Global.Add("gslsflog1plusx","(",new OneOperator1<double ,double>( gsl_sf_log_1plusx )); 
   Global.Add("gslsflog1plusxmx","(",new OneOperator1<double ,double>( gsl_sf_log_1plusx_mx )); 
   Global.Add("gslsfpowint","(",new OneOperator2<double,double,long >( gsl_sf_pow_int )); 
   Global.Add("gslsfpsiint","(",new OneOperator1<double ,long>( gsl_sf_psi_int )); 
   Global.Add("gslsfpsi","(",new OneOperator1<double ,double>( gsl_sf_psi )); 
   Global.Add("gslsfpsi1piy","(",new OneOperator1<double ,double>( gsl_sf_psi_1piy )); 
   Global.Add("gslsfpsi1int","(",new OneOperator1<double ,long>( gsl_sf_psi_1_int )); 
   Global.Add("gslsfpsi1","(",new OneOperator1<double ,double>( gsl_sf_psi_1 )); 
   Global.Add("gslsfpsin","(",new OneOperator2<double,long,double >( gsl_sf_psi_n )); 
   Global.Add("gslsfsynchrotron1","(",new OneOperator1<double ,double>( gsl_sf_synchrotron_1 )); 
   Global.Add("gslsfsynchrotron2","(",new OneOperator1<double ,double>( gsl_sf_synchrotron_2 )); 
   Global.Add("gslsftransport2","(",new OneOperator1<double ,double>( gsl_sf_transport_2 )); 
   Global.Add("gslsftransport3","(",new OneOperator1<double ,double>( gsl_sf_transport_3 )); 
   Global.Add("gslsftransport4","(",new OneOperator1<double ,double>( gsl_sf_transport_4 )); 
   Global.Add("gslsftransport5","(",new OneOperator1<double ,double>( gsl_sf_transport_5 )); 
   Global.Add("gslsfsin","(",new OneOperator1<double ,double>( gsl_sf_sin )); 
   Global.Add("gslsfcos","(",new OneOperator1<double ,double>( gsl_sf_cos )); 
   Global.Add("gslsfhypot","(",new OneOperator2<double,double,long >( gsl_sf_hypot )); 
   Global.Add("gslsfsinc","(",new OneOperator1<double ,double>( gsl_sf_sinc )); 
   Global.Add("gslsflnsinh","(",new OneOperator1<double ,double>( gsl_sf_lnsinh )); 
   Global.Add("gslsflncosh","(",new OneOperator1<double ,double>( gsl_sf_lncosh )); 
   Global.Add("gslsfanglerestrictsymm","(",new OneOperator1<double ,double>( gsl_sf_angle_restrict_symm )); 
   Global.Add("gslsfanglerestrictpos","(",new OneOperator1<double ,double>( gsl_sf_angle_restrict_pos )); 
   Global.Add("gslsfzetaint","(",new OneOperator1<double ,long>( gsl_sf_zeta_int )); 
   Global.Add("gslsfzeta","(",new OneOperator1<double ,double>( gsl_sf_zeta )); 
   Global.Add("gslsfzetam1","(",new OneOperator1<double ,double>( gsl_sf_zetam1 )); 
   Global.Add("gslsfzetam1int","(",new OneOperator1<double ,long>( gsl_sf_zetam1_int )); 
   Global.Add("gslsfhzeta","(",new OneOperator2<double,double,long >( gsl_sf_hzeta )); 
   Global.Add("gslsfetaint","(",new OneOperator1<double ,long>( gsl_sf_eta_int )); 
   Global.Add("gslsfeta","(",new OneOperator1<double ,double>( gsl_sf_eta )); 
 } 
/*****************/
/*****************/
