// Regression tests
// ----------------

// $Id$

include(../regtests.m4)

// The values tested here may not have a physical or mathematical
// meaning. Their main property is to gather numerical values from the
// whole domain, to be checked for consistency with previous runs.

dnl (syntax of ONETEST macro defined in ../regtests.m4)

ONETEST(adapt,u[].max-u[].min,0.01)
ONETEST(adaptindicatorP1,u[].max-u[].min,0.01)
ONETEST(adaptindicatorP2,u[].max-u[].min,0.01)
ONETEST(algo,umax,0.01)
ONETEST(algowithmacro,umax,0.01)
ONETEST(array)
ONETEST(a_tutorial,1+max(err[].max,-err[].min),0.001)
ONETEST(beam,uu[]'*uu[],5e-2)
ONETEST(BlackSchole,normvL2,0.1)
ONETEST(calculus)
ONETEST(cavity,psi[]'*psi[],1e-2)
ONETEST(convect2,v[]'*v[],1e-1)
ONETEST(convect-apt,error,0,1e-2)
ONETEST(convect,v[]'*v[],1e-1)
ONETEST(dumptable)
ONETEST(ex-vf)
ONETEST(FE,wdc[]'*wdc[],1e-2)
ONETEST(fluidStructAdapt,uu[]'*uu[],2e-1)
ONETEST(fluidStruct,uu[]'*uu[],2e-1)
ONETEST(freeboundary,u[]'*u[],5e-2)
ONETEST(freeboundary-weak,p[]'*p[],5e-2)
ONETEST(LapDG2,u[]'*u[],1e-2)
ONETEST(Laplace,uh[]'*uh[],1e-2)
ONETEST(LaplaceP1bis,u[]'*u[],1e-2)
ONETEST(LaplaceP1,uh[]'*uh[],1e-2)
ONETEST(LaplaceP1P2h,u2h[]'*u2h[],1e-2)
ONETEST(LaplaceRT,u1[]'*u1[],1e-2)
ONETEST(mesh)
ONETEST(movemesh,u[]'*u[],1e-2)
ONETEST(nolinear-elas,un[]'*un[],1e-2)
ONETEST(NSUzawaCahouetChabart,u1[]'*u1[],1e-2)
ONETEST(onde,u[]'*u[],1e-2)

dnl The following two tests have suspicious results (1e20 and bigger)
dnl ONETEST(periodic4,uh[]'*uh[],1e-2)
dnl ONETEST(Periodic,uh[]'*uh[],1e-2)

ONETEST(plot,uh[]'*uh[],1e-2)
ONETEST(readmesh,u[]'*u[],2e-1)
ONETEST(region,u[]'*u[],1e-1)
ONETEST(saverestore)
ONETEST(schwarz-gc,u1[]'*u1[],5e-2)
ONETEST(schwarz-no-overlap,u[]'*u[],5e-2)
ONETEST(schwarz-overlap,u[]'*u[],5e-2)
ONETEST(sparse-matrix,xx[]'*xx[],5e-2)
ONETEST(sparse-cmatrix,real(xx[]'*xx[]),5e-2)
ONETEST(StokesUzawa,u1[]'*u1[],5e-2)
ONETEST(tablefunction,fxy[]'*fxy[],1e-2)
