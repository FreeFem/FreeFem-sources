// Regression tests
// ----------------

// $Id$

include(../regtests.m4)

// The values tested here may not have a physical or mathematical
// meaning. Their main property is to gather numerical values from the
// whole domain, to be checked for consistency with previous runs.

dnl (syntax of ONETEST macro defined in ../regtests.m4)

ONETEST(aadaptation,u[]'*u[],0.1);
ONETEST(aalapacien,u[]'*u[],0.1);
ONETEST(aalaplace-nc,u[]'*u[],0.1);
ONETEST(aamove,u[]'*u[],0.1);
ONETEST(aaRT,u1[]'*u1[],0.1);
ONETEST(arrayoFVh,u[]'*u[],0.1);
ONETEST(bilap,xx(0:n)'*xx(0:n),0.1);
ONETEST(D2,w[]'*w[],1e-20,1e-20);
ONETEST(demo1,u[]'*u[],0.1);
ONETEST(demo,u[]'*u[],0.1);
ONETEST(funct,myfunction(1.0,3.),0.1);
ONETEST(lapacienprecon,u[]'*u[],0.1);
ONETEST(lap_mat,u1[]'*u1[],0.1);
ONETEST(NSP1P1b,u1[]'*u1[],0.2);
ONETEST(NSP1P1,u1x[]'*u1x[],0.1);
ONETEST(NSP1P2,u1[]'*u1[],0.1);
ONETEST(parareal,pu'*pu,0.1);

dnl alh - 30/6/04 - unstable?
dnl ONETEST(Richard,h[]'*h[],0.1);

ONETEST(teste,P.x,0.1);
ONETEST(testFE);
ONETEST(wafer-heating-laser-axi,xx'*xx,0.1);
