// Regression tests
// ----------------

// $Id$

include(../regtests.m4)

// The values tested here may not have a physical or mathematical
// meaning. Their main property is to gather numerical values from the
// whole domain, to be checked for consistency with previous runs.

// (syntax of ONETEST macro defined in ../regtests.m4)

ONETEST(aadaptation,u[]'*u[],0.13,0.15);
ONETEST(aalapacien,u[]'*u[],122058,122060);
ONETEST(aalaplace-nc,u[]'*u[],1713.47,1713.49);
ONETEST(aamove,u[]'*u[],0.068,0.070847);
ONETEST(aaRT,u1[]'*u1[],280.924,280.926);
ONETEST(arrayoFVh,u[]'*u[],4560.73,4560.75);
ONETEST(bilap,xx(0:n)'*xx(0:n),33.834,33.836);
ONETEST(D2,w[]'*w[],-1e-30,1e-30);
ONETEST(demo1,u[]'*u[],1.3,1.4);
ONETEST(demo,u[]'*u[],1.32,1.35);
ONETEST(funct,myfunction(1.0,3.),3.999999,4.000001);
ONETEST(lapacienprecon,u[]'*u[],122058,122060);
ONETEST(lap_mat,u1[]'*u1[],85.3775,85.3777);
ONETEST(NSP1P1b,u1[]'*u1[],85,121);
ONETEST(NSP1P1,u1x[]'*u1x[],3.99,5.59);
ONETEST(NSP1P2,u1[]'*u1[],161.2743646,195);
ONETEST(parareal,pu'*pu,1383.26,1383.28);

dnl alh - 30/6/04 - unstable?
dnl ONETEST(Richard,h[]'*h[],1383.26,1383.28);

ONETEST(teste,P.x,0.999999,1.000001);
ONETEST(testFE);
ONETEST(wafer-heating-laser-axi,xx'*xx,284,286);
