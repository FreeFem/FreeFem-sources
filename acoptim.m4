# Choosing debugging and/or optimization flags for compilation
# ------------------------------------------------------------

AC_ARG_ENABLE(profiling,[  --enable-profiling	Turn on profiling])
if test "$enable_profiling" = yes
then
	CXXFLAGS="$CXXFLAGS -pg"
	LDFLAGS="$LDFLAGS -pg"
fi

# Debug mode (no optimisation)
# ----------------------------

AC_MSG_CHECKING(whether to generate debugging information)

AC_ARG_ENABLE(debug,[  --enable-debug	Turn on debug versions of FreeFem++])
AC_ARG_ENABLE(optim,[  --enable-optim	Turn on compiler optimization])

# Autoconf always chooses -O2. -O2 in gcc makes some functions
# disappear. This is not ideal for debugging. And when we optimize, we
# do not use -O2 anyway.

CFLAGS="`echo $CFLAGS | sed 's/-O2//g'`"
FFLAGS="`echo $FFLAGS | sed 's/-O2//g'`"
CXXFLAGS="`echo $CXXFLAGS | sed 's/-O2//g'`"

if test "$enable_debug" = yes;
then
	AC_MSG_RESULT(yes)
else
	AC_MSG_RESULT(no)

	# No debugging information in optimized code

	CFLAGS="`echo $CFLAGS | sed 's/-g//g'` -DNDEBUG"
	FFLAGS="`echo $FFLAGS | sed 's/-g//g'` -DNDEBUG"
	CXXFLAGS="`echo $CXXFLAGS | sed 's/-g//g'` -DNDEBUG"
fi

# Hardware-independant optimization
# ---------------------------------

if test "$enable_debug" != yes -a "$enable_optim" != no;
then
	CHECK_COMPILE_FLAG(C,-O3,CFLAGS)
	CHECK_COMPILE_FLAG(C++,-O3,CXXFLAGS)
	CHECK_COMPILE_FLAG(Fortran 77,-O3,FFLAGS)
fi

AC_ARG_ENABLE(generic,
[  --enable-generic	Turn off hardware-dependant optimization options])

# Generic code
if test "$enable_debug" != yes \
    -a "$enable_optim" != no \
    -a "$enable_generic" = yes
then
	CHECK_COMPILE_FLAG(C,-mcpu=common,CFLAGS)
	CHECK_COMPILE_FLAG(C++,-mcpu=common,CXXFLAGS)
	CHECK_COMPILE_FLAG(Fortran 77,-mcpu=common,FFLAGS)
fi

# Hardware-dependant optimization
# -------------------------------

if test "$enable_debug" != yes \
    -a "$enable_optim" != no \
    -a "$enable_generic" != yes
then

    # MacOS X
    if test -x /usr/bin/hostinfo
	then


	# If we are on MacOS X to choise the optimisaztion 
	AC_MSG_CHECKING(GCC version)

        ff_gcc4=`gcc  --version |awk  ' NR==1 {print $3}'|sed -e 's/\..*$//'` 
	AC_MSG_RESULT($ff_gcc4)

	# At the moment, we do not know how to produce correct
	# optimizated code on G5.
	AC_MSG_CHECKING(PowerPC architecture)

        # CPU detection
	ff_cpu=unkown
	ff_optim_type=
	if test `/usr/bin/hostinfo|grep ppc7450|wc -l` -gt 0
	    then
	    ff_cpu=G4
	    ff_optim_type=-G4
	elif test `/usr/bin/hostinfo|grep ppc970|wc -l` -gt 0
	    then
	    ff_cpu=G5
	    ff_optim_type=-G5
	fi
	if test $ff_cpu = unknown;
	    then
	    AC_MSG_ERROR(cannot determine PowerPC cpu type)
	fi
	AC_MSG_RESULT($ff_cpu)


	if test `/usr/bin/hostinfo|grep Darwin|wc -l` -gt 0 \
		-a $ff_cpu != G5
	    then

	    # Optimization flags: -fast option do not work because the
	    # -malign-natural flags create wrong IO code
            if test  "$ff_gcc4" -eq 4 
	    then
             ff_fast='-fast '
            else
	    ff_fast='-funroll-loops -fstrict-aliasing -fsched-interblock -falign-loops=16 -falign-jumps=16 -falign-functions=16 -falign-jumps-max-skip=15 -falign-loops-max-skip=15 -ffast-math -mdynamic-no-pic -mpowerpc-gpopt -force_cpusubtype_ALL -fstrict-aliasing  -mpowerpc64 '
	    fi
	    # Specific PowerPC G5 optimization
	    if test "$ff_cpu" = G5;
		then

	        # remove -fstrict-aliasing on G5 to much optim the
	        # code cash in GC

		ff_fast="`echo $ff_fast| sed 's/-fstrict-aliasing //g'`"

	    fi

	    CHECK_COMPILE_FLAG(C,$ff_fast,CFLAGS)
	    CHECK_COMPILE_FLAG(C++,$ff_fast,CXXFLAGS)
	    CHECK_COMPILE_FLAG(Fortran 77,$ff_fast,FFLAGS)
	fi

	# CPU reference: PowerPC G4
	if test "$ff_cpu" = G4;
	    then
	    CHECK_COMPILE_FLAG(C,-mcpu=7450,CFLAGS)
	    CHECK_COMPILE_FLAG(C++,-mcpu=7450,CXXFLAGS)
	    CHECK_COMPILE_FLAG(Fortran 77,-mcpu=7450,FFLAGS)

	# CPU reference: PowerPC G5
	elif test "$ff_cpu" = G5;
	    then

	    # but at least this way we can see
	    # that the automatic detection worked.

	    CHECK_COMPILE_FLAG(C,-mcpu=G5,CFLAGS)
	    CHECK_COMPILE_FLAG(C++,-mcpu=G5,CXXFLAGS)
	    CHECK_COMPILE_FLAG(Fortran 77,-mcpu=G5,FFLAGS)
	fi

    # Linux
    elif test -f /proc/cpuinfo
	then

	# Specific processors
	proc_type=unknown
	ff_optim_type=
	if test `grep 'Pentium III (Coppermine)' /proc/cpuinfo|wc -l` -gt 0
	    then
	    proc_type=pentium3
	    ff_optim_type=-P3
	elif test `grep 'Intel(R) Pentium(R) III ' /proc/cpuinfo|wc -l` -gt 0
	    then
	    proc_type=pentium3
	    ff_optim_type=-P3
	elif test `grep 'Intel(R) Pentium(R) 4 ' /proc/cpuinfo|wc -l` -gt 0
	    then
	    proc_type=pentium4
	    ff_optim_type=-P4
	elif test `grep 'Intel(R) Xeon(TM) CPU' /proc/cpuinfo|wc -l` -gt 0
	    then
	    proc_type=pentium4
	    ff_optim_type=-P4
	elif test `grep 'AMD Athlon(tm) Processor' /proc/cpuinfo|wc -l` -gt 0
	    then
	    proc_type=athlon
	    ff_optim_type=-Athlon
	elif test `grep 'AMD Athlon(tm) XP' /proc/cpuinfo|wc -l` -gt 0
	    then
	    proc_type=athlon-xp
	    ff_optim_type=-AthlonXP
	fi

	if test "$proc_type" != unknown
	    then
	    CHECK_COMPILE_FLAG(C,-march=$proc_type,CFLAGS)
	    CHECK_COMPILE_FLAG(C++,-march=$proc_type,CXXFLAGS)
	    CHECK_COMPILE_FLAG(Fortran 77,-march=$proc_type,FFLAGS)
	fi

	# If we did not find a processor type (this happens with
	# cygwin), try and select separate capabilities instead.

	if test "$proc_type" = unknown
	    then
	    if test `grep -e '^flags.*mmx' /proc/cpuinfo|wc -l` -gt 0
		then
		CHECK_COMPILE_FLAG(C,-mmmx,CFLAGS)
		CHECK_COMPILE_FLAG(C++,-mmmx,CXXFLAGS)
		CHECK_COMPILE_FLAG(Fortran 77,-mmmx,FFLAGS)
	    fi
	    if test `grep -e '^flags.*sse ' /proc/cpuinfo|wc -l` -gt 0
		then
		CHECK_COMPILE_FLAG(C,-msse,CFLAGS)
		CHECK_COMPILE_FLAG(C++,-msse,CXXFLAGS)
		CHECK_COMPILE_FLAG(Fortran 77,-msse,FFLAGS)
	    fi
	    if test `grep -e '^flags.*sse2' /proc/cpuinfo|wc -l` -gt 0
		then
		CHECK_COMPILE_FLAG(C,-msse2,CFLAGS)
		CHECK_COMPILE_FLAG(C++,-msse2,CXXFLAGS)
		CHECK_COMPILE_FLAG(Fortran 77,-msse2,FFLAGS)
	    fi
	    if test `grep -e '^flags.*3dnow' /proc/cpuinfo|wc -l` -gt 0
		then
		CHECK_COMPILE_FLAG(C,-m3dnow,CFLAGS)
		CHECK_COMPILE_FLAG(C++,-m3dnow,CXXFLAGS)
		CHECK_COMPILE_FLAG(Fortran 77,-m3dnow,FFLAGS)
	    fi
	fi
    fi
fi

# Defines a variable containing the optimization type, to be used in
# binary archive names. It may be empty if only generic optimization
# is used.

AC_SUBST(OPTIM_TYPE,$ff_optim_type)
