# Choosing debugging and/or optimization flags for compilation
# ------------------------------------------------------------

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
	if test `/usr/bin/hostinfo|grep ppc970|wc -l` -gt 0
	    then

	    # -fast used by gcc on PowerPC G5
	    CHECK_COMPILE_FLAG(C,-fast,CFLAGS)
	    CHECK_COMPILE_FLAG(C++,-fast,CXXFLAGS)
	    CHECK_COMPILE_FLAG(Fortran 77,-fast,FFLAGS)
	fi

    # Linux
    elif test -f /proc/cpuinfo
	then

	# Processors
	proc_type=unknown
	if test `grep 'Intel(R) Pentium(R) 4 ' /proc/cpuinfo|wc -l` -gt 0
	    then
	    proc_type=pentium4
	    CHECK_COMPILE_FLAG(C,-march=pentium4,CFLAGS)
	    CHECK_COMPILE_FLAG(C++,-march=pentium4,CXXFLAGS)
	    CHECK_COMPILE_FLAG(Fortran 77,-march=pentium4,FFLAGS)
	elif test `grep 'AMD Athlon(tm) XP' /proc/cpuinfo|wc -l` -gt 0
	    then
	    proc_type=athlon-xp
	    CHECK_COMPILE_FLAG(C,-march=athlon-xp,CFLAGS)
	    CHECK_COMPILE_FLAG(C++,-march=athlon-xp,CXXFLAGS)
	    CHECK_COMPILE_FLAG(Fortran 77,-march=athlon-xp,FFLAGS)
	fi

	# If we did not find a processor type (this happens with
	# cygwin), try and select separate capabilities instead.

	if test "$proc_type" = no
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
