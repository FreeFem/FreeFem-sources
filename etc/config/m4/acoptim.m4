# Choosing debugging and/or optimization flags for compilation
# ------------------------------------------------------------
#   get CPU Type
cputype=unknow
if test -x /usr/bin/machine ; then
  cputype=`/usr/bin/machine`
elif test -x /usr/bin/arch ; then
  cputype=`/usr/bin/arch`
fi
cpuintel=no;
case "$cputype" in 
i386|i486| x86_64*)  cpuintel=yes;;
*) cpuintel=no;;
esac
AC_MSG_NOTICE([    -----   CPU kind: $cputype , intel/amd: $cpuintel ])
AC_ARG_ENABLE(profiling,[  --enable-profiling	Turn on profiling])
if test "$enable_profiling" = yes
then
	CXXFLAGS="$CXXFLAGS -pg"
	LDFLAGS="$LDFLAGS -pg"
fi

if test "$enable_m64" = yes -a "$enable_m32" 
then
  	    AC_MSG_ERROR([ Choose  32 or 64 architecture not the both ],1);  
fi
AC_ARG_ENABLE(m64,[  --enable-m64	Turn on 64 bits architecture])
if test "$enable_m64" = yes
then
	ff_m64=-m64	
	ff_ok=no
        CHECK_COMPILE_FLAG(C,$ff_m64,CFLAGS,ff_ok)
        if test "$ff_ok" = yes ;then  CNOFLAGS="$CFLAGS $ff_m64";fi
	CHECK_COMPILE_FLAG(C++,$ff_m64,CXXFLAGS)
	CHECK_COMPILE_FLAG(Fortran 77,$ff_m64,FFLAGS)	
#  add -fPIC on on 64 architecture 
   if test "$ff_ok" = yes -a "$ff_fpic" != "no" ;then    
        CHECK_COMPILE_FLAG(C,-fPIC,CFLAGS,ff_ok)
	CHECK_COMPILE_FLAG(C++,-fPIC,CXXFLAGS)
	CHECK_COMPILE_FLAG(Fortran 77,-fPIC,FFLAGS)	
  fi

fi
AC_ARG_ENABLE(m32,[  --enable-m32	Turn on 32 bits architecture])
if test "$enable_m32" = yes 
then
	ff_m32=-m32	
	ff_ok=no
        CHECK_COMPILE_FLAG(C,$ff_m32,CFLAGS,ff_ok)
        if test "$ff_ok" = yes ;then  CNOFLAGS="$CFLAGS $ff_m32";fi
        CHECK_COMPILE_FLAG(C,$ff_m32,CNOFLAGS)
	CHECK_COMPILE_FLAG(C++,$ff_m32,CXXFLAGS)
	CHECK_COMPILE_FLAG(Fortran 77,$ff_m32,FFLAGS)	
#  add -fPIC on on 64 architecture 
#        CHECK_COMPILE_FLAG(C,-fPIC,CFLAGS)
#	CHECK_COMPILE_FLAG(C++,-fPIC,CXXFLAGS)
#	CHECK_COMPILE_FLAG(Fortran 77,-fPIC,FFLAGS)	

fi

# Debug mode (no optimisation)
# ----------------------------

AC_MSG_CHECKING(whether to generate debugging information)

AC_ARG_ENABLE(debug,[  --enable-debug	Turn on debug versions of FreeFem++])
AC_ARG_ENABLE(optim,[  --enable-optim	Turn on compiler optimization])

if test "$enable_debug" = yes;
then

	AC_MSG_RESULT(yes)
	CFLAGS="`echo $CFLAGS | sed 's/-O2//g'`"
	FFLAGS="`echo $FFLAGS | sed 's/-O2//g'`"
	CXXFLAGS="`echo $CXXFLAGS | sed 's/-O2//g'`"

else
	AC_MSG_RESULT(no)

	# No debugging information in optimized code

	CFLAGS="$CFLAGS -DNDEBUG"
	FFLAGS="$FFLAGS -DNDEBUG"
	CXXFLAGS="$CXXFLAGS -DNDEBUG"
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

# FFCS: remove "-mcpu=common" to allow other hardware-dependant values of cpu for PowerPC - thank you Fred (20/02/11)
# FH 
    # Generic code
if test  "$enable_generic" = yes 
    then
	CHECK_COMPILE_FLAG(C,-mtune=generic,CFLAGS)
	CHECK_COMPILE_FLAG(C++,-mtune=generic,CXXFLAGS)
	CHECK_COMPILE_FLAG(Fortran 77,-mtune=generic,FFLAGS)
fi


# Hardware-dependant optimization
# -------------------------------

if test "$enable_debug" != yes \
    -a "$enable_optim" != no \
    -a "$enable_generic" != yes
then

# Autoconf always chooses -O2. -O2 in gcc makes some functions
# disappear. This is not ideal for debugging. And when we optimize, we
# do not use -O2 anyway.

CFLAGS="`echo $CFLAGS | sed 's/-O2//g'`"
FFLAGS="`echo $FFLAGS | sed 's/-O2//g'`"
CXXFLAGS="`echo $CXXFLAGS | sed 's/-O2//g'`"

	    if test `grep -e '^flags.*mmx' /proc/cpuinfo|wc -l` -gt 0
		then
		CHECK_COMPILE_FLAG(C,-mmmx,CFLAGS)
		CHECK_COMPILE_FLAG(C++,-mmmx,CXXFLAGS)
		CHECK_COMPILE_FLAG(Fortran 77,-mmmx,FFLAGS)
	    fi
	    if test `grep -e '^flags.*avx2' /proc/cpuinfo|wc -l` -gt 0
		then
		CHECK_COMPILE_FLAG(C,-mavx2,CFLAGS)
		CHECK_COMPILE_FLAG(C++,-mavx2,CXXFLAGS)
		CHECK_COMPILE_FLAG(Fortran 77,-mavx2,FFLAGS)
        fi
	    if test `grep -e '^flags.*avx\>' /proc/cpuinfo|wc -l` -gt 0
		then
		CHECK_COMPILE_FLAG(C,-mavx,CFLAGS)
		CHECK_COMPILE_FLAG(C++,-mavx,CXXFLAGS)
		CHECK_COMPILE_FLAG(Fortran 77,-mavx,FFLAGS)
        fi
	    if test `grep -e '^flags.*sse4_2' /proc/cpuinfo|wc -l` -gt 0
		then
		CHECK_COMPILE_FLAG(C,-msse4.2,CFLAGS)
		CHECK_COMPILE_FLAG(C++,-msse4.2,CXXFLAGS)
		CHECK_COMPILE_FLAG(Fortran 77,-msse4.2,FFLAGS)
	    fi
	    if test `grep -e '^flags.*sse2' /proc/cpuinfo|wc -l` -gt 0
	    then
	    CHECK_COMPILE_FLAG(C,-msse2,CFLAGS)
	    CHECK_COMPILE_FLAG(C++,-msse2,CXXFLAGS)
	    CHECK_COMPILE_FLAG(Fortran 77,-msse2,FFLAGS)
		fi
	    if test `grep -e '^flags.*sse\>' /proc/cpuinfo|wc -l` -gt 0
		then
		CHECK_COMPILE_FLAG(C,-msse,CFLAGS)
		CHECK_COMPILE_FLAG(C++,-msse,CXXFLAGS)
		CHECK_COMPILE_FLAG(Fortran 77,-msse,FFLAGS)
	    fi
fi

# Defines a variable containing the optimization type, to be used in
# binary archive names. It may be empty if only generic optimization
# is used.

AC_SUBST(OPTIM_TYPE,$ff_optim_type)
AC_MSG_NOTICE([     CXXFLAGS =   $CXXFLAGS  ])
