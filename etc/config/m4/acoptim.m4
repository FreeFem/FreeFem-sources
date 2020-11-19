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
        CHECK_COMPILE_FLAG(C,-g,CFLAGS)
	CHECK_COMPILE_FLAG(C++,-g,CXXFLAGS)
	CHECK_COMPILE_FLAG(Fortran 77,-g,FFLAGS)	

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

    # MacOS X Darwin
    if test -x /usr/bin/hostinfo
	then
        

	# If we are on MacOS X to choise the optimisaztion 
	AC_MSG_CHECKING(GCC version)

        ff_gcc4=`$CC  --version |awk  ' NR==1 {print $3}'|sed -e 's/\..*$//'` 
	ff_clang=`$CC  --version |awk  '/clang/  {print $4}'`
	if test -n "$ff_clang" ; then ff_gcc4="llvm"; fi
	AC_MSG_RESULT($ff_gcc4)

	# At the moment, we do not know how to produce correct
	# optimizated code on G5.
	AC_MSG_CHECKING(PowerPC architecture)
	ff_machine=`(test -x /usr/bin/machine && /usr/bin/machine) || echo unknow`
        ff_fast="-O3"
	if test	-n "$ff_clang" ; then
          ff_fast='-O3'
	elif test `uname` = Darwin 
	    then
	    # Optimization flags: -fast option do not work because the
	    # -malign-natural flags create wrong IO code
            if test  "$ff_gcc4" -eq 4 
	    then
               ff_fast='-fast'
            else
	      ff_fast='-O3 -funroll-loops -fstrict-aliasing -fsched-interblock -falign-loops=16 -falign-jumps=16 -falign-functions=16 -falign-jumps-max-skip=15 -falign-loops-max-skip=15 -ffast-math  -mpowerpc-gpopt -force_cpusubtype_ALL -fstrict-aliasing  -mpowerpc64 '
	    fi
	fi        


        # CPU detection

	case $ff_machine  in
	  ppc7450) # G4
		ff_fast="$ff_fast -mtune=G4 -mcpu=G4";;
          ppc970) # G5 
	        # remove -fstrict-aliasing on G5 to much optim the
	        # code cash in GC
		ff_fast="`echo $ff_fast -mtune=G5 -mcpu=G5| sed 's/-fstrict-aliasing //g'`";;
          ppc*) # G3 ????
	       ff_fast="-O3";;
	  i486)
	    ff_fast="-O3 $ff_fast";;
	  x86_64*)
	  ff_fast="-O3 $ff_fast";;
	  arm64*)
	  ff_fast="-O3 $ff_fast";;
	  *)
	    AC_MSG_ERROR(cannot determine apple cpu type )
	    ff_fast="-O3";;
	 esac


	AC_MSG_RESULT($ff_fast)

        CHECK_COMPILE_FLAG(C,$ff_fast,CFLAGS)
	CHECK_COMPILE_FLAG(C++,$ff_fast,CXXFLAGS)
	CHECK_COMPILE_FLAG(Fortran 77,$ff_fast,FFLAGS)


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
	    if test `grep -e '^flags.*avx' /proc/cpuinfo|wc -l` -gt 0
		then
		CHECK_COMPILE_FLAG(C,-mavx,CFLAGS)
		CHECK_COMPILE_FLAG(C++,-mavx,CXXFLAGS)
		CHECK_COMPILE_FLAG(Fortran 77,-mavx,FFLAGS)
	    else
	        if test `grep -e '^flags.*sse4_2' /proc/cpuinfo|wc -l` -gt 0
		    then
		    CHECK_COMPILE_FLAG(C,-msse4.2,CFLAGS)
		    CHECK_COMPILE_FLAG(C++,-msse4.2,CXXFLAGS)
		    CHECK_COMPILE_FLAG(Fortran 77,-msse4.2,FFLAGS)
	        else
	            if test `grep -e '^flags.*sse2' /proc/cpuinfo|wc -l` -gt 0
	                then
	                CHECK_COMPILE_FLAG(C,-msse2,CFLAGS)
	                CHECK_COMPILE_FLAG(C++,-msse2,CXXFLAGS)
	                CHECK_COMPILE_FLAG(Fortran 77,-msse2,FFLAGS)
				else
	            	if test `grep -e '^flags.*sse ' /proc/cpuinfo|wc -l` -gt 0
		        	then
		        	CHECK_COMPILE_FLAG(C,-msse,CFLAGS)
		        	CHECK_COMPILE_FLAG(C++,-msse,CXXFLAGS)
		        	CHECK_COMPILE_FLAG(Fortran 77,-msse,FFLAGS)
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
	fi
    fi
fi

# Defines a variable containing the optimization type, to be used in
# binary archive names. It may be empty if only generic optimization
# is used.

AC_SUBST(OPTIM_TYPE,$ff_optim_type)
AC_MSG_NOTICE([     CXXFLAGS =   $CXXFLAGS  ])
