// Regression tests
// ----------------

// The tests are checked against reference values by "make check" in
// each examples subdirectory

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

define(REFFILE,"ref.edp")
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
