// Regression tests
// ----------------

// The tests are checked against reference values by "make check"

// "ref.edp" contains all reference values and may be rebuilt with
// "make Ref"

// $Id$

// The values tested here may not have a physical or mathematical
// meaning. Their main property is to gather numerical values from the
// whole domain, to be checked for consistency with previous runs.

NoUseOfWait=true;
int verbosityy=verbosity;

dnl May write or read a reference file
changequote([[,]])

dnl UMFPACK give results that are different enough
dnl to create a specific reference file
ifelse(UMFPACKLIB,,
	[[define(REFFILE,"refnoumf.edp")]],
	[[define(REFFILE,"refumf.edp")]])
ifdef([[ASSERT]],
	include REFFILE;,
	ofstream ref(REFFILE);)

dnl $1=file name
dnl $2=reference value (if there is one)
dnl $3=precision of reference value (if there is one)
dnl	or minimum absolute variation if $4 is defined
dnl $4=maximum absolute variation (if defined)

define(ONETEST,
[[cout << "--------- file : $1.edp -----------------" << endl;
verbosity=verbosityy;
{
	dnl Place the dash first to avoid any confusion with things like "a-z"
	define([[TESTVAR]],TEST[[]]translit($1,-_,XX))
	define([[REFVAR]],REF[[]]translit($1,-_,XX))
	include "$1.edp";
	ifelse($2,,,
		[[real TESTVAR=$2;
		ifdef([[ASSERT]],
			cout<<"$1 reference value = "<<REFVAR
				<<" test value ="<<TESTVAR<<endl;
			ifelse($4,,
				assert(TESTVAR<REFVAR*(1+$3));
				assert(TESTVAR>REFVAR*(1-$3));,
				assert(TESTVAR<REFVAR+$4);
				assert(TESTVAR>REFVAR-$3);),
			ref<<"real REFVAR="<<TESTVAR<<";"<<endl;)]])
};
]])

ONETEST(adapt,u[]'*u[],5e-2)
ONETEST(adaptindicatorP1,u[]'*u[],5e-2)
ONETEST(adaptindicatorP2,u[]'*u[],1,1)
ONETEST(algo)
ONETEST(algowithmacro)
ONETEST(array)
ONETEST(a_tutorial,u[]'*u[],5e-2)
ONETEST(beam,uu[]'*uu[],1e-2)
ONETEST(blakschol,u[]'*u[],1)
ONETEST(calculus)
ONETEST(cavity,psi[]'*psi[],1e-2)
ONETEST(convect2,v[]'*v[],5e-2)
ONETEST(convect-apt,v[]'*v[],5e-2)
ONETEST(convect,v[]'*v[],5e-2)
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
ONETEST(readmesh,u[]'*u[],5e-2)
ONETEST(region,u[]'*u[],1e-2)
ONETEST(saverestore)
ONETEST(schwarz-gc,u1[]'*u1[],5e-2)
ONETEST(schwarz-no-overlap,u[]'*u[],5e-2)
ONETEST(schwarz-overlap,u[]'*u[],5e-2)
ONETEST(sparse-matrix,xx[]'*xx[],5e-2)
ONETEST(StokesUzawa,u1[]'*u1[],5e-2)
ONETEST(tablefunction,fxy[]'*fxy[],1e-2)
