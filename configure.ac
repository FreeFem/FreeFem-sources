# ------------------------------------------------------------

# Antoine Le Hyaric - LJLL Paris 6 - lehyaric@ann.jussieu.fr - 13/5/04
# $Id: configure.ac,v 1.338 2010/05/10 21:35:07 hecht Exp $

# Version numbering: x.xx-pp where "pp" is the package version (when
# the same FreeFem++ version is packaged several times). For
# coherency, this should have the same value as the top-most package
# number in debian/changelog.


AC_INIT(FreeFem++,3.47,hecht@ann.jussieu.fr,freefem++)
dnl : ${CFLAGS=""}
dnl : ${CXXFLAGS=""}
dnl : ${FCFLAGS=""}
dnl : ${FFLAGS=""}

# progg test  expact under windows 
ff_TEST_FFPP="../src/nw/FreeFem++"
ff_TEST_FFPP_MPI="../src/mpi/FreeFem++-mpi"

AC_PREREQ(2.69) dnl for AC_LANG_CASE and check stuff 
dnl AC_CONFIG_SRCDIR(src/FreeFem++-CoCoa)

# Automake 1.11 is too old for check ...

# ALH - 10/9/13 - FFCS needs to skip this automake version requirement when compiling on older platforms (eg cygwin or
# last Ubuntu LTS). So we fix a low hardcoded requirement and test for higher versions when not in FFCS.

AM_INIT_AUTOMAKE(1.11 dist-zip)

if test "$enable_ffcs" != yes
then
    if test `echo $am__api_version|awk '{if($1>1.13)print "ok";}'` != ok
    then
	AC_MSG_ERROR([Automake version needs to be 1.13 or later to enable "make check"])
    fi
fi

AC_CONFIG_HEADERS(config.h)
AC_PROG_MAKE_SET
dnl AM_COLOR_TESTS=always
AC_PROG_RANLIB

ff_where_lib_conf=examples++-load/WHERE_LIBRARY-config
dnl search of associad software
m4_define([AC_FF_ADDWHERELIB],
  [	if test -z "$ff_where_lib_conf_$1" ; then
	    echo "$1 LD $2"  >>$ff_where_lib_conf
	    test -n "$3" && echo "$1 INCLUDE $3 "  >>$ff_where_lib_conf
            ff_where_lib_conf_$1=1
	    ff_where_lib_$1="$2"
	    ff_where_inc_$1="$3"
	    AC_MSG_NOTICE([    ++ add $1 : $2 $3 in  $ff_where_lib_conf "])
	else
	    AC_MSG_NOTICE([    -- do not add $1 : $2 $3 in  $ff_where_lib_conf "])
        fi
]) dnl
m4_define([AC_FF_WHERELIB],
    [  	if test -z "$ff_where_lib_conf_$1" ; then
        AC_MSG_CHECKING(check $1)
	ff_save_libs=$LIBS
	LIBS="$LIBS $2 $4"
	AC_LINK_IFELSE([AC_LANG_PROGRAM([m4_if($3,,,[#include <$3>])],[])], [ff_WHERE=yes],[ff_WHERE=no])
	if test "$ff_WHERE" = "yes" ; then
	    echo "$1 LD $2"  >>$ff_where_lib_conf
	    ff_WHERE_INC=`AS_DIRNAME(["$3"])`
	    case "$ff_WHERE_INC" in 
		/*) echo "$1 INCLUDE -I$ff_WHERE_INC"  >>$ff_where_lib_conf ;;
	    esac
            ff_where_lib_conf_$1=1
	fi
	AC_MSG_RESULT($ff_WHERE)
        ff_$1_ok=$ff_WHERE;
	LIBS=$ff_save_libs
	fi
    ])dnl end m4_define AC_FF_WHERELIB
define([AC_FF_PETSC_AWK],[awk -F'=' '"$1 "==[$]1 {print [$]2}' [$]ffconfpetsc])
m4_define([AC_FF_PETSC_WHERELIB],
	  [
	    ff_$2=`AC_FF_PETSC_AWK($2)`
	    ff_$3=`AC_FF_PETSC_AWK($3)`
	   if test -n "[$]ff_$2"; then 
			    ff_$1_ok=yes;
			    AC_FF_ADDWHERELIB($1,[$]ff_$2,[$]ff_$3)
			    AC_SUBST([TOOL_COMPILE_$1],"")
	  fi	    
	  ]
)dnl end m4_defineAC_FF_PETSC_WHERELIB

# ALH - this is required by FFCS that needs to deactivate some tools that do not work on all platforms. Some FF users
# may also find interesting to specify a local version of a tool instead of downloading it.

# m4 macro parameters: $1 = tool name, $2 = dynamic library name, $3 = download directory name if different from $1

m4_define([TOOL_PARAMETERS],
    [AC_ARG_WITH($1_include,AC_HELP_STRING([--with-$1-include=],[Include directives for $1 instead of automatic download]))
    AC_ARG_WITH($1_ldflags,AC_HELP_STRING([--with-$1-ldflags=],[Link-time directives for $1 instead of automatic download]))

    if test "$with_$1_include" != "" || test "$with_$1_ldflags" != ""
    then

	# some directives have been specified, use them instead of downloading

	AC_FF_ADDWHERELIB($1,$with_$1_ldflags,$with_$1_include)
        AC_SUBST([TOOL_COMPILE_$1],"")
	enable_$1_download=no
    fi
])

m4_define([TOOL_DISABLE],
    [AC_ARG_ENABLE($1,AC_HELP_STRING([--disable-$1],[Do not use $1]))
    if test "$enable_$1" = "no"
    then
        AC_SUBST([TOOL_COMPILE_$1],"")
        AC_SUBST([TOOL_DYLIB_$1],"")
    elif  test   "$ff_$1_ok" = "yes" ; then
        AC_SUBST([TOOL_DYLIB_$1],$2)
        enable_$1_download=no
    else  
        AC_SUBST([TOOL_COMPILE_$1],ifelse($3,,$1,$3))
        AC_SUBST([TOOL_DYLIB_$1],$2)
	    enable_$1_download=yes
    fi

    # Also allow to disable the download of one tool if it is already locally installed

    TOOL_PARAMETERS($1,$2,$3)
])

ff_HOSTARCH_pastix=
echo "#  Build with freefem++ with ./configure " `date`  >$ff_where_lib_conf
# To allow anonymous CVS version to contain a "./configure" and
# Makefiles

AM_MAINTAINER_MODE

echo >config_LIB_INFO

# Necessary compilers

AC_PROG_CC
AC_PROG_CXX
AC_LANG(C++)
if test `uname` != Darwin; then
 ff_ldeg="-Wl,--end-group"
 ff_ldbg="-Wl,--start-group"
fi

#  suffix of  dynamic lib .. 
# ---------------------------
ff_uname=`uname`
#  flag to build window 32   version   ff_mingw = yes
# bof bof F. Hecht 
case $ff_uname in
    CYGWIN*|MINGW*|MSYS_NT*)
	ff_suffix_dylib="dll"
	ff_suffix_dylib_a="dll.a";;
     Darwin)
	ff_suffix_dylib="dylib"
	ff_suffix_dylib_a="dylib";;
    *) 
       ff_suffix_dylib="so";
       ff_suffix_dylib_a="so";
esac 
# end suffix ...



AC_CHECK_PROG(ff_m4,m4,yes,no)
AC_CHECK_PROG(ff_bison,bison,yes,no)
AC_CHECK_PROG(ff_flex,flex,yes,no)
dnl AC_CHECK_PROG(ff_cmake,cmake,yes,no)
AC_CHECK_PROG(ff_patch,patch,yes,no)
case "$ff_m4 $ff_bison $ff_patch $ff_flex"  in
*no*)
    AC_MSG_NOTICE([ to install missing package  under debian or ubuntu, try ])
    test "$ff_m4" = no && AC_MSG_NOTICE([ sudo apt-get install m4])
    test "$ff_bison" = no && AC_MSG_NOTICE([ sudo apt-get install bison])
    test "$ff_flex" = no && AC_MSG_NOTICE([ sudo apt-get install flex])
    test "$ff_patch" = no && AC_MSG_NOTICE([ sudo apt-get install patch])
    AC_MSG_ERROR([ Sorry missing  m4,bison,flex,patch command !]);;
esac


if test -z "$CXX"; then
  AC_MSG_NOTICE( [ fatal error : c++ compiler ! ] );
  AC_MSG_ERROR([ Sorry no c++ compiler !])
fi

AC_COMPUTE_INT(ff_size_long,[sizeof(long)])
AC_COMPUTE_INT(ff_size_int,[sizeof(int)])
AC_COMPUTE_INT(ff_size_ptr,[sizeof(int*)])

AC_SUBST(SIZEOF_LONG,$ff_size_long)
AC_SUBST(SIZEOF_INT,$ff_size_int)
AC_SUBST(SIZEOF_PTR,$ff_size_ptr)
if test "0$ff_size_ptr" -eq 4  ;then
 AC_SUBST(SIZEOF_PTRINBIT,32)
elif test  "0$ff_size_ptr" -eq 8  ;then
 AC_SUBST(SIZEOF_PTRINBIT,64)
else
  AC_MSG_NOTICE( [ fatal error : sizeof pointer $ff_size_ptr !  or no c++ compiler: $CXX] );
  AC_MSG_ERROR([ Sorry  sizeof c++ pointer $ff_size_ptr are not  4 or 8 ])
fi

# FFCS - build the code for FreeFem++-cs
AC_ARG_ENABLE(ffcs,AC_HELP_STRING([--enable-ffcs],[build FreeFem++ for use by FreeFem++-cs]))
if test "$enable_ffcs" = yes
then
	AC_DEFINE_UNQUOTED(ENABLE_FFCS,$enable_ffcs,[build FreeFem++ for use by FreeFem++-cs])
else
	enable_ffcs=no
fi
AC_SUBST(ENABLE_FFCS,"$enable_ffcs")
AM_CONDITIONAL([ENABLE_FFCS],[test $enable_ffcs = yes])

# dur dur car sous MacOsX le fortran n'est pas standard.
ff_AR="ar"
ff_ARFLAGS="rv"
ff_RANLIB="ranlib"
AC_ARG_ENABLE(fortran,AC_HELP_STRING([--disable-fortran],[No Fortran compiler available ( ARPACK need it)]))

ff_g2c_lib="";
if test "$enable_fortran" != no
then

    # ALH-FFCS-2/3/10: add gfortran-mp-4.4 for MacPorts on MacOS 10.6
    # FH add iforr
    AC_PROG_FC(gfortran  f90  xlf90 g95 ifort gfortran-mp-4.4)
    AC_PROG_F77(gfortran f90 xlf xlf90 g95 f77 fort77 "$FC" ifort gfortran-mp-4.4)

#	if test -n "$F77"
#	then
    ff_flibs=""
        # modif FH
    AC_F77_LIBRARY_LDFLAGS
    dnl	AC_F77_WRAPPERS
	# correct pb of double def under macos
    case $F77 in
	 *fort77)
	    if test -z "$FLIBS" ; then
#             FLIBS=`fort77 -v  a.out  2>&1|awk '/a.out/ && /fort77/  { print $(NF-2),$(NF-1)}'`
#  FH to remove " "
		FLIBS=`fort77 -v  a.out  2>&1|awk '/a.out/ && /fort77/  { print "echo",$(NF-2),$(NF-1)}'|sh`

		AC_MSG_WARN([  fort77   FLIBS :  $FLIBS ])
	    fi
	    ;;
	*g77)
	    for i in $FLIBS; do
		case $i in
# save last directory of the list
		    -L*)
			d=`expr $i : '-L\(.*\)'`;
			echo " try $d "
			if test -e "$d/libg2c.so" ; then
			    ff_flibs="$d/libg2c.so"
			elif test -e "$d/libg2c.dylib" ; then
			    ff_flibs="$d/libg2c.dylib"
			elif test -e "$d/libg2c.a" ; then
			    ff_flibs="$d/libg2c.a"
			elif test  -e $d/libg2c.so -o	-e $d/libg2c.a ; then
			    ff_flibs="$i"
			fi;;
		esac
	    done
	    if test -e "$ff_flibs"  ;then
		FLIBS="$ff_flibs"
	    else
		FLIBS="$ff_flibs -lg2c"
	    fi
	    AC_MSG_WARN([  get dir of -lg2c  FLIBS :  $FLIBS ])
	    ;;
# add FH sep 2006 / modif  2009
        *gfortran)
	    ff_okkk=0
            for i in $FLIBS; do
		case $i in
		    -L*)
			d=`expr $i : '-L\(.*\)'`;
			if test -e "$d/libgfortran.$ff_suffix_dylib_a" -a  -e "$d/libquadmath.$ff_suffix_dylib_a"  ; then
			   ff_flibs="$d/libgfortran.$ff_suffix_dylib_a $d/libquadmath.$ff_suffix_dylib_a"
			    ff_okkk=1
			elif test -e "$d/libgfortran.a"  -a -e "$d/libquadmath.a" ; then
			   ff_flibs="$d/libgfortran.a $d/libquadmath.a"
			    ff_okkk=1
			elif test -e "$d/libgfortran.$ff_suffix_dylib_a"   ; then
			   ff_flibs="$d/libgfortran.$ff_suffix_dylib_a"
			    ff_okkk=2
			elif test -e "$d/libgfortran.a" ; then
			   ff_flibs="$d/libgfortran.a"
			    ff_okkk=2
			fi;;
		esac
            done
	    if test  "$ff_okkk" -ge  1  ;then
		FLIBS="$ff_flibs"
	    else
		FLIBS="$ff_flibs -lgfortran"
	    fi
            AC_MSG_WARN([  get dir of -lgfortran  FLIBS :  $FLIBS ])
            ;;
    esac
    AC_ARG_WITH(flib,[  --with-flib= the fortran library ])
#  correct FH sep 2006 -o -> -a
    if test "$with_flib" != no -a -n  "$with_flib"
    then
	ff_g2c_lib="$with_flib"
	FLIBS="$with_flib"
    fi
# add FH oct 2007 for download f2c
    if test -z "$F77"  ; then
	ff_f77=`pwd`/download/bin/fort77
	ff_flibs=`pwd`/download/lib/libf2c.a
	if test -x $ff_f77  -a -f $ff_flibs ; then
	    AC_MSG_WARN([ no fortran, but find download f2c/fort]);
	    F77=$ff_f77
	    FLIBS=$ff_flibs
	else
	    AC_MSG_NOTICE( [ fatal error : no fortran ] );
	    AC_MSG_NOTICE( [add --disable-fortran  ] );
	    AC_MSG_NOTICE( [or try to compile f2c in directory download/f2c  ] );
	    AC_MSG_NOTICE( [  just do:  make install ] );
	    AC_MSG_ERROR([ Fatal error No Fortran compiler . ],1);
	fi
    fi
#  check if the FLIBS is correct
    ff_libs="$LIBS"
    LIBS="$ff_libs $FLIBS"
    AC_TRY_LINK_FUNC(exit,ff_err=,ff_err=ok);
    if test "$ff_err" = "ok" ; then  AC_MSG_ERROR([ Fatal FLIBS: $FLIBS is incorrect. ],1); fi
    LIBS="$ff_libs"
    echo "F77 LD $ff_libs" >config_LIB_INFO

    AC_MSG_CHECKING([ Size of fortran 77 integer ])



##   AC_LANG(fortran);
    AC_LANG_PUSH([Fortran 77])
    AC_RUN_IFELSE([AC_LANG_SOURCE([[
		program test
		integer  p,i
		p=1024*1024
		i= p*p
		open(7,file="conftest.out")
		if (i>0) then
		    write(7,*) 8
		else
         write(7,*) 4
	 endif
	 close(7)
	 end
 ]])] ,ff_f77intsize=`cat conftest.out`,ff_f77intsize=4,ff_f77intsize=4)
## AC_LANG(C++)
    AC_LANG_POP([Fortran 77])
    AC_MSG_RESULT($ff_f77intsize)
    AC_F77_WRAPPERS

    ff_cfagsf77=""
    case $ac_cv_f77_mangling in
#  "lower case, no underscore, no extra underscore")
#  "lower case, no underscore, extra underscore")
# ;;
	"lower case, underscore, no extra underscore") 	ff_cfagsf77="-DAdd_";;
	"lower case, underscore, extra underscore")       	ff_cfagsf77="-DAdd__ -Df77IsF2C";;
	"upper case, no underscore, no extra underscore")     ff_cfagsf77="-DUPPER";;
	"upper case, no underscore, extra underscore")        ff_cfagsf77="-DUPPER -DAdd__";;
	"upper case, underscore, no extra underscore")      ff_cfagsf77="-DUPPER -DAdd_";;
	"upper case, underscore, extra underscore")           ff_cfagsf77="-DUPPER -DAdd__";;
	*) ;;
    esac
    AC_SUBST(CFLAGSF77,$ff_cfagsf77)
    if test $ff_f77intsize -ne  4
    then
	AC_MSG_ERROR([ Fatal Error: Your fortran 77  compiler have by not 4 bytes integer ( $ff_f77intsize bytes) ],1);
    fi
    test -f /mingw/bin/libpthread-2.dll && ff_pthread_dll=/mingw/bin/libpthread-2.dll ; \
    AC_FF_ADDWHERELIB(f77,$FLIBS,$ff_cfagsf77 $ff_pthread_dll)
    AC_FF_ADDWHERELIB(fc,$FLIBS,$ff_cfagsf77 $ff_pthread_dll)
fi

#  fin test FORTRAN ..........
# ----------------------------

AC_ARG_ENABLE(c,AC_HELP_STRING([--disable-c],[No C compiler available (C BLAS need it)]))
if test "$enable_c" != no
then
    AC_PROG_CC
    AM_PROG_CC_C_O
else

	# We need to define this otherwise ./configure crashes with
	# the message "configure: error: conditional "am__fastdepCC"
	# was never defined".

    am__fastdepCC_TRUE=
    am__fastdepCC_FALSE='#'
fi

AM_PROG_LEX
AC_PROG_YACC
AC_LANG(C++)


# Some useful libraries

AC_CHECK_LIB(pthread,pthread_create,ff_pthread="-lpthread",ff_pthread="")
AC_CHECK_LIB(iomp5,omp_get_wtime,ff_iomp5="-liomp5",ff_iomp5="")
AC_CHECK_LIB([rt],[clock_gettime])

if test -n "$ff_pthread"  ; then
    AC_FF_ADDWHERELIB(pthread,$ff_pthread,)
fi

# Necessary absolute pathname for local directory when some libraries
# are used from several different locations (for instance locally in
# configure.ac and in a subdir).
curdir=`pwd`

# Configure options
# -----------------

echo $ac_configure_args > configure.param
echo "$PATH" > config.path
eval for i in $ac_configure_args\; do echo \$i\; done | sed -e "s/$/'/" -e "s/^/'/"   > config.param

# by default the suffix of the  .so lib file is .so
#      it is .dll   under windows
#      it is .dylib under macos X

# suffix of dynamic lib 
# Checking wether we can produce a MacIntosh-specific version
# -----------------------------------------------------------
AC_MSG_CHECKING(wether we are on a MacIntosh)
ff_mac=no
if test `uname` = Darwin;
then
	ff_HOSTARCH_pastix=i686_mac
	ff_suffix_dylib="dylib"
	ff_suffix_dylib_a="dylib"
	ff_mac=yes
AC_DEFINE_UNQUOTED(FF_PREFIX_DIR_APPLE,"/Applications/FreeFem++.app/Contents/",FreeFem prefix  dir)
fi
AC_MSG_RESULT($ff_mac)

AC_MSG_CHECKING(wether we are on  SunOS)
ff_sunos=no
if test `uname -s` = SunOS;
then
	ff_sunos=yes
fi
AC_MSG_RESULT($ff_sunos)

# glut ------------
# ----------------
ff_glut=""
ff_glut_ok=""
AC_ARG_WITH(glut,[  --with-glut  glutlib and include  ])
if test "$with_glut" != no -a  -n "$with_glut" ; then
    ff_libs_save="$LIBS"
    LIBS="$LIBS"
dnl not works .... so no check .... FH 
      ff_glut="$with_glut";ff_glut_ok=yes
dnl     AC_SEARCH_LIBS(abort,[$with_glut],[ff_glut="$with_glut";ff_glut_ok=yes],[ff_glut="";ff_glut_ok="no"])
    LIBS="$ff_libs_save"
dnl     AC_MSG_ERROR([ Sorry  bad --with-glut :  $with_glut !])
elif test "$with_glut" = no ; then
    ff_glut_ok="no"
fi

# Checking wether we can produce a Microsoft Windows-specific version
# -------------------------------------------------------------------

AC_ARG_ENABLE(cygwindll,[  --enable-cygwindll	Forces the use of the Cygwin DLL (not recommended)])
AC_ARG_ENABLE(mingw64,[  --enable-mingw64	Uses mingw64 compilers on Cygwin])
AC_MSG_CHECKING(wether we are on Microsoft Windows)
ff_uname=`uname`
#  flag to build window 32   version   ff_mingw = yes
ff_mingw=no
ff_fpic=yes
case $ff_uname in
    CYGWIN*)
	ff_fpic=no
	#    ff_nocygwin=-mno-cygwin
	AC_SUBST(GCCNOCYGWIN,$ff_nocygwin);;
    MINGW*|MSYS_NT*)
        enable_cygwindll=no;;


esac
case $ff_uname in
CYGWIN*|MINGW*)
	ff_fpic=no
	ff_suffix_dylib="dll";
	ff_suffix_dylib_a="dll.a";
        ff_win32=yes;

	# FFCS - 8/3/12 - remove -D_MSC_VER under MinGW64 because it forces system calls to be compiled into any object
	# (which creates thousands of duplicate definitions for sytem calls like time()).


	# We need Mingw to avoid Cygwin's extra DLLs
	if test "$enable_cygwindll" != yes
	then
# 	        CHECK_COMPILE_FLAG(C++,-mwindows,CXXFLAGS)

		# FFCS: on Windows, FF crashes when compiling GL/glut.h and the option "--disable-opengl" is not
		# operational because ff_glut_ok is forced to yes here.
		if test $enable_ffcs = yes
		then
	            ff_glut_ok=no
		    ff_glut=""
		    enable_opengl=no
		fi

                ff_mingw=yes
                enable_cygwindll=no;
		ff_pthread="-mthreads"

		# FFCS does not use the Cygwin MinGW compilers any more
		if test $enable_ffcs = no
		then
		    CXXFLAGS="$CXXFLAGS $ff_nocygwin -I/usr/include/mingw"
		    FFLAGS="$FFLAGS $ff_nocygwin"
		    CFLAGS="$CFLAGS $ff_nocygwin -I/usr/include/mingw"
                    AC_COMPILE_IFELSE([AC_LANG_SOURCE([int a;])],[],
			[ff_nocygwin="";
			    AC_MSG_NOTICE([Sorry $ff_nocygwin optio is wrong try whitout , but try with gcc-3.3])
			    ])
		    CXXFLAGS="$CXXFLAGS $ff_nocygwin -I/usr/include/mingw"
		    FFLAGS="$FFLAGS $ff_nocygwin"
		    CFLAGS="$CFLAGS $ff_nocygwin -I/usr/include/mingw"
		    CNOFLAGS="$CNOFLAGS $ff_nocygwin -I/usr/include/mingw"
		fi

		LIBS="$LIBS $ff_nocygwin -mthreads -lws2_32 -lcomdlg32"
		LIBSNOCONSOLE="-mwindows"

		# FFCS uses a specific compiler, so we specify its libraries explicitely
		if test $enable_ffcs = no
		then
		    test -z "$MPIRUN" &&  MPIRUN=`which mpiexe.exe`
		    if test "$enable_fortran" != no  -o  "$with_flib" != no ;   then
			case "$F77" in
	 		    *gfortran) FLIBS="$ff_nocygwin -lgfortran -lquadmath";;
	 		    *g77) FLIBS="$ff_nocygwin -lg2c";;
			    *)   AC_MSG_ERROR([ Sorry no known FLIBS with this $F77  !]) ;;
			esac
		    fi
		fi

		if test -z "$ff_glut" -a  "$ff_glut_ok" != "no" ; then 
                    ff_glutname="glut32"
		#  check abort a existing function just to find in glut32.dll exist in the path
		#  because glutInit is not the real symbol on win32 dur dur FH !!!!!!!!!
		    AC_CHECK_LIB(glut32,abort,ff_glut="-l$ff_glutname -mthreads -lglu32 -lopengl32",ff_glut="")
		    if test -z "$ff_glut" ; then
			ff_glutname="freeglut"
			if test -x  /usr/bin/pkg-config.exe ; then
			ff_glut="-lglu32 "`/usr/bin/pkg-config.exe --libs freeglut`
			else
			ff_glut="-l$ff_glutname -mthreads -lglu32 -lopengl32"
			fi
			AC_CHECK_LIB(freeglut,abort,ff_glut_ok=yes,ff_glut="")
		    fi
		fi
         	# Resources for FreeFem++-cs in Microsoft Windows format
	        AC_SUBST(FFGLUTNAME,$ff_glutname)
	        AC_SUBST(WINDRESOBJ,windres.o)
	        AC_SUBST(LIBSNOCONSOLE,$LIBSNOCONSOLE)
            AC_SUBST(WIN32DLLTARGET,win32-dll-target)
		    ff_TEST_FFPP="../src/bin-win32/FreeFem++.exe"
		    ff_TEST_FFPP_MPI="../src/bin-win32/FreeFem++-mpi.exe"

            AC_DEFINE(PURE_WIN32,1,A pure windows applications no cygwin dll)
		ff_dynload=yes
	fi


	;;
*) ff_win32=no;;
esac
AC_MSG_RESULT($ff_win32)

if test "$ff_win32" = no
then
	enable_cygwindll=no
fi

# FreeFem++-specific version information
# --------------------------------------

# Version numbering, converted to a floating point value
ff_numver="`echo $VERSION|sed 's/-\(.*\)/+\1.0*0.000001/'`"
AC_DEFINE_UNQUOTED(VersionFreeFempp,$ff_numver,FreeFem++ version as a float)

# Since src/fflib/strversionnumber.cpp is recreated at each build, this
# date is only useful for config-version.h

test "$prefix" = NONE && prefix="$ac_default_prefix"
ff_prefix_dir="${prefix}/lib/ff++/$VERSION"
AC_MSG_CHECKING(prefix dir freefem++  )
AC_MSG_RESULT($ff_prefix_dir)

FF_DATE=`date`
AC_DEFINE_UNQUOTED(VersionFreeFemDate,"$FF_DATE",FreeFem++ build date)
AC_DEFINE_UNQUOTED(FF_PREFIX_DIR,"${ff_prefix_dir}",FreeFem prefix  dir)
AC_SUBST(ff_prefix_dir,$ff_prefix_dir)

# Separate configuration header file for version information, included
# in config-macos9.h

cat << EOF > config-version.h
/* FreeFem++ version information for MacOS 9 configuration */
/* File generated by configure.ac */

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "$PACKAGE_STRING"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "$PACKAGE_TARNAME"

/* Define to the version of this package. */
#define PACKAGE_VERSION "$PACKAGE_VERSION"

/* Version number of package */
#define VERSION "$VERSION"

/* FreeFem++ build date */
#define VersionFreeFemDate "$FF_DATE"

/* FreeFem++ version as a float */
#define VersionFreeFempp ($ff_numver)
EOF

# Checking wether to produce a debug version
# ------------------------------------------

# Generic or hardware-dependant optimization
m4_include(acmacros.m4)
m4_include(acoptim.m4)

# [[file:ax_lib_hdf5.m4]]
m4_include(ax_lib_hdf5.m4)

# [[file:ax_lib_gsl.m4]]
m4_include(ax_lib_gsl.m4)

# Checking c++11 for schwarz/hpddm  interface of P. Jolivet ...
# F.Hecht:
# do not compile all with c++11 version to much 
# probleme to day (frev 2015) with g++ 4.8.1
ff_cxx_save=$CXXFLAGS
CHECK_COMPILE_FLAG(C++,[-std=c++11],CXXFLAGS,ff_cxx11)
dnl CXXFLAGS="$CXX11FLAGS $CXXFLAGS"
dnl ff_cxx11="no"
dnl if test -n "$CXX11FLAGS" ; then ff_cxx11="yes" ; fi
if test "$ff_cxx11" = yes ; then
    CXX11FLAGS="-std=c++11";
    AC_FF_ADDWHERELIB(cxx11,$CXX11FLAGS,$CXX11FLAGS)
fi
CXXFLAGS="$ff_cxx_save"
AC_SUBST(WITH_CXX11,$ff_cxx11)
AC_SUBST(CXX11FLAGS,$CXX11FLAGS)
#end Checking  c++11
# ------------------------------
# compilation flags of FreeFem++
# ------------------------------
#  -DDRAWING  bamg active some drawing facility in bamg (inquire mesh)
#  -DBAMG_LONG_LONG active the use of  long long in bamg
#  -DDEBUG active the  checking in bamg
#  -DNCHECKPTR remove check pointer facility
#  -DCHECK_KN active subscripting of some array
#  -DWITHCHECK of the dynamic type ckecking of the langague (very slow)
#  -DEIGENVALUE  to compile the eigen value part
#  ------------------------------------------------

# FreeFrem++-specific optimizations -fpascal-strings


if test "$enable_debug" = yes;
then
	CXXFLAGS="$CXXFLAGS  -DBAMG_LONG_LONG -DCHECK_KN"
	if test "$ff_mac" = yes;
	then
		CXXFLAGS="$CXXFLAGS -fno-inline  -fexceptions"
	fi
else
	CXXFLAGS="$CXXFLAGS -DBAMG_LONG_LONG  -DNCHECKPTR"
	if test "$ff_mac" = yes;
	then
		CXXFLAGS="$CXXFLAGS "
	fi
fi

#whether or not to add a suffix to the package name
#--------------------------------------------------

AC_ARG_WITH(suffix,[  --with-suffix=[G3/G4/G4,P4,..] or --without-suffix (default no suffix to package name)])

AC_MSG_CHECKING(suffix to add to package name)
ff_suffix=
if test "$with_suffix" = yes -o "$with_suffix" = no -o "$with_suffix" = ""
then
  AC_MSG_RESULT(none)
else
  ff_suffix="-$with_suffix"
  AC_MSG_RESULT($ff_suffix)
fi
AC_SUBST(ADD_PACKAGE_NAME,$ff_suffix)




# Checking that OpenGL is accessible
# ----------------------------------

ff_gl_ok=no
ff_glx_ok=no

ff_libs="$LIBS"
AC_ARG_WITH(flib,[  --with-flib= the fortran library ])
AC_ARG_ENABLE([opengl],AS_HELP_STRING([--enable-opengl],
                                      [Enable/disable support for OpenGL]))

# echo "  .... LIBS = $LIBS .... "
if test x$enable_opengl != xno; then
# GL libraries seem to depend on some dlopen capabilities
if test "$enable_static" != yes -a "$ff_glut_ok" != yes
then
	# If we cannot find gl.h, no OpenGL version can be produced
    AC_CHECK_HEADERS(OpenGL/gl.h,
        ff_gl_ok=yes
        ,
        ff_gl_ok=no
        )


    if test "$ff_gl_ok" = no
    then
	AC_CHECK_HEADERS(GL/gl.h,
	    ff_gl_ok=yes
	    ,
	    ff_gl_ok=no
	    )
    fi
    AC_CHECK_HEADERS(GLUT/glut.h,
        ff_glut_ok=yes,
        ff_glut_ok=no)

    if test "$ff_glut_ok" = no
    then
	AC_CHECK_HEADERS(GL/glut.h,
	    ff_glut_ok=yes,
	    ff_glut_ok=no)
    fi

# GLUT .....
    if test "$ff_glut_ok" =	yes
    then
	for glut in \
	    "-framework GLUT -framework OpenGL -framework Cocoa" \
	    "-lglut -lGLU -lGL"
	do
	    ff_glut="$glut"
	    ff_libs1="$LIBS"
     	    LIBS="$LIBS $ff_glut"
	    AC_TRY_LINK_FUNC(glutInit,
		ff_glut_ok=yes
	        LIBS="$ff_libs1"
		break,
		ff_glut_ok=next);
	    LIBS="$ff_libs1"
	done
    fi


fi
fi
AC_SUBST(LIBSPTHREAD,$ff_pthread)

# ALH - 16/9/13 - at the moment the requirements for medit (glut) are not compatible with FFCS so medit is
# only compiled with FF.

if test $enable_ffcs = no
then
    if test "$ff_glut_ok" = yes
    then
	AC_SUBST(LIBSGLUT,$ff_glut)
	if test -n  "$ff_pthread" ; then
	    ff_ffglutprog="ffglut${EXEEXT}"
	    AC_DEFINE_UNQUOTED(PROG_FFGLUT,"$ff_ffglutprog", the ffglut application for the new graphics )
	fi
	ff_meditprog="ffmedit${EXEEXT}"
    fi
fi

AC_ARG_WITH(suffix,[  --with-suffix=[G3/G4/G4,P4,..] or --without-suffix (default no suffix package name)])
if test "$with_suffix" = yes
then
  with_suffix=""
elif test "$with_suffix" = no -o -z "$with_suffix"
then
  with_suffix=""
else
  with_suffix="-$with_suffix"
fi
  AC_MSG_CHECKING(add suffix $with_suffix )
  AC_SUBST(ADD_PACKAGE_NAME,$with_suffix)
  AC_MSG_RESULT(yes)


# Checking wether we can produce a parallel version
# -------------------------------------------------

if test $enable_ffcs = no
then
    # FF case
    m4_include(acmpi.m4)
else
    # FFCS - use the same MPI configuration choices as FFCS
    if test -z "$MPICXX"
    then
	ff_mpi=no
    else
	ff_mpi=yes
	AC_SUBST(MPICXX,$MPICXX)
	AC_SUBST(MPICC,$MPICC)
	AC_SUBST(MPIF77,$MPIF77)
	AC_SUBST(MPIFC,$MPIFC)
	AC_SUBST(MPIPROG,"FreeFem++-mpi${EXEEXT}")
	AC_SUBST(MPI_INCLUDE,"-I $MPI_INC_DIR")
	AC_SUBST(MPI_INC_DIR,$MPI_INC_DIR)
	AC_SUBST(MPI_LIB_DIRS,"")
	AC_SUBST(MPI_LIB,$MPI_LIB)
	AC_SUBST(MPI_LIBC,"")
	AC_SUBST(MPI_LIBFC,"")

        # Extra MPI-dependant configuration options that are set by FF
        # during MPI configuration. FFCS - 25/2/13 - Fred noticed that if
        # PASTIX_HOSTARCH stays blank, pastix compilation breaks. At least
        # i686_pc_linux and i686_mac are required by pastix on the
        # corresponding platforms.

	if test"$ff_HOSTARCH_pastix" = ""
	then
	    ff_HOSTARCH_pastix=i686_pc_linux
	fi
	AC_SUBST(PASTIX_HOSTARCH,$ff_HOSTARCH_pastix)

        # these values should not be empty otherwise
        # examples++-load/ff-get-dep will think that they are not
        # defined
	AC_FF_ADDWHERELIB(mpi,-DDUMMY,-I$MPI_INC_DIR)

        # mpifc and mpif77 libraries should always be specified
        # because FF never calls the Fortran MPI compiler. It always
        # uses mpicxx in [[file:examples++-load/ff-c++]]. The
        # resulting Fortran libraries (eg Mumps) would compile even
        # without the proper Fortran libs, but they would not load
        # properly.

        # under Win32, libmpi_f77.a is not the right name and FFCS
        # mingw/mpicxx adds the right libraries by itself

        # With mpich2 Debian wheezy 32 bits, libmpi_f77.a is not available (request from Cico, 14/10/13)

	if test -f /usr/lib/openmpi/lib/libmpi_f77.a
	then
            AC_FF_ADDWHERELIB(mpifc,-lmpi_f77,)
    	    AC_FF_ADDWHERELIB(mpif77,-lmpi_f77,)
	else
            AC_FF_ADDWHERELIB(mpifc,-DDUMMY,)
    	    AC_FF_ADDWHERELIB(mpif77,-DDUMMY,)
	fi

        # FFCS - MPI_DOUBLE_COMPLEX kept from original FF configure script
	AC_MSG_CHECKING( MPI_DOUBLE_COMPLEX)
	AC_COMPILE_IFELSE(
	    [AC_LANG_SOURCE([
#include <mpi.h>
			MPI_Datatype xxxx=MPI_DOUBLE_COMPLEX;
			])],
	    ff_mpi_double_complex=yes,
	    ff_mpi_double_complex=no)
	AC_MSG_RESULT($ff_mpi_double_complex)
	if test "$ff_mpi_double_complex" = yes  ; then
	    AC_DEFINE(HAVE_MPI_DOUBLE_COMPLEX,1, mpi_double_complex)
	fi
    fi
fi

# FFCS needs to change some of the FF makefiles to compile without MPI
AM_CONDITIONAL([FFCS_MPIOK],[test $ff_mpi = yes])

# Looking for useful configuration utilities
# ------------------------------------------

AC_ARG_ENABLE(download,AC_HELP_STRING([--enable-download],[Download missing libraries (BLAS,ARPACK,UMFPACK,...)]))
## PETSc
AC_ARG_WITH(petsc,[  --with-petsc=/usr/local/petsc/conf/petscvariables   --without-petsc  ])
# /usr/local/petsc/conf/petscvariables 
#/usr/local/ff++/petsc/lib/petsc/conf/petscvariables
ff_petsc_ok=no
## echo --$with_petsc--
if test "$with_petsc" != no -a  $ff_mpi = yes; then
    for d in "$with_petsc" /usr/local/ff++/petsc/lib  /usr /usr/local /opt/usr  /opt/local ;do
	for dd in "" "/petsc/conf/petscvariables" "/lib/petsc-conf/petscvariables" ; do
	    ffconfpetsc="$d$dd"		
	    if test -f "$ffconfpetsc" ; then 
		if test  "$ff_petsc_ok" = no ; then
		    PETSC_MPIRUN=`awk -F' *= *' '"MPIEXEC"==$1 {print $2}' "$ffconfpetsc"`
		    test "$PETSC_MPIRUN" = "$MPIRUN" ||  cmp -s  `which "$PETSC_MPIRUN"` `which "$MPIRUN"`
		    if test "$?"  -eq 0 ; then
		        AC_FF_PETSC_WHERELIB(petsc,PETSC_WITH_EXTERNAL_LIB,PETSC_CC_INCLUDES)
		        AC_FF_PETSC_WHERELIB(parmetis,PARMETIS_LIB,PARMETIS_INCLUDE)
		        AC_FF_PETSC_WHERELIB(metis,METIS_LIB,METIS_INCLUDE)
		        AC_FF_PETSC_WHERELIB(ptscotch,PTSCOTCH_LIB,PTSCOTCH_INCLUDE)
		        AC_FF_PETSC_WHERELIB(scotch,SCOTCH_LIB,SCOTCH_INCLUDE)
		        AC_FF_PETSC_WHERELIB(fftw3_mpi,FFTW_LIB,FFTW_INCLUDE)
		        AC_FF_PETSC_WHERELIB(mumps,MUMPS_LIB,MUMPS_INCLUDE)
		        AC_FF_PETSC_WHERELIB(scalapack,SCALAPACK_LIB,SCALAPACK_INCLUDE)
			AC_FF_PETSC_WHERELIB(suitesparse,SUITESPARSE_LIB,SUITESPARSE_INCLUDE)
			AC_FF_PETSC_WHERELIB(blaslapack,BLASLAPACK_LIB,BLASLAPACK_INCLUDE)
			
		#	PETSC_WITH_EXTERNAL_LIB=`awk -F'=' '"PETSC_WITH_EXTERNAL_LIB "==$1 {print $2}' $ffconfpetsc`
		#	PETSC_CC_INCLUDES=`awk -F'=' '"PETSC_CC_INCLUDES "==$1 {print $2}' $ffconfpetsc`
		#	if test -n "$PETSC_WITH_EXTERNAL_LIB"; then 
		#	    ff_petsc_ok=yes;
		#	    AC_FF_ADDWHERELIB(petsc,$PETSC_WITH_EXTERNAL_LIB,$PETSC_CC_INCLUDES)
		#	fi
		    else
			AC_MSG_NOTICE([" Warning PETSC MPI and FF++ MPI not the same: $PETSC_MPIRUN != $MPIRUN." ]);			
		    fi
		fi
	    else
		# test  "$ff_petsc_ok" = no && echo "  *** try $ffconfpetsc"
		ffconfpetsc="";

	    fi
	done
    done
fi


if test "$ff_petsc_ok" != no ;   then
    AC_MSG_NOTICE([ with petsc  $ff_petsc_ok])
else
    AC_MSG_NOTICE([ without petsc ***** ])
fi

if test "$ff_where_lib_conf_fftw3_mpi" = 1  ; then
  ff_where_fftw3=`echo $ff_where_lib_fftw3_mpi| sed -e 's/-lfftw3_mpi//g'` 
  AC_FF_ADDWHERELIB(fftw3,$ff_where_fftw3,$ff_where_inc_fftw3_mpi)
  AC_MSG_NOTICE([ fftw3 from  petsc ***** ]) 
fi


#replacing wget with another command
ff_wget_command="wget --no-check-certificate"
ff_wget=no

AC_ARG_WITH(wget,
	[  --with-wget=command	Replace "wget" with another command. Implies --enable-download],
	   enable_download=yes
	   ff_wget=yes
	   ff_wget_command="${withval}")

ff_mkl_libpath=yes
AC_ARG_WITH(mkl,
	[  --with-mkl= the  MKL LIBPATH : (ie. /opt/intel/mkl/RR.r.y.xxx/lib/<arch>)],
	   ff_mkl_libpath="${withval}")

# if enabling downloads find wget or curl to do download
# ------------------------------------------------------

if test  "$ff_wget" =  no
then
	ff_wget=""
	AC_CHECK_PROG(ff_wget,wget --no-check-certificate,yes,no)
fi
if test  "$ff_wget" = no
then
        ff_wget_command="curl --fail -O "
  	ff_curl=""
	AC_CHECK_PROG(ff_curl,curl -O,yes,no)
	ff_wget="$ff_curl"
fi

# for automatique compilation of
# lib in download
if test  "$ff_wget" != yes  -a "$enable_download"  = "yes"  ; then
  enable_download=no
  AC_MSG_ERROR([   enable-download and no wget or curl. ],1);
fi
if test "$enable_download"  = "yes" ; then
DOWNLOADCOMPILE="compile-pkg install-other"
else
DOWNLOADCOMPILE=
fi

AC_SUBST(WGET,$ff_wget_command)
AC_SUBST(DOWNLOADCOMPILE,$DOWNLOADCOMPILE)

#  modif FH -----
#  -- looking for cadna
#  the round-off error propagation
# the web site http://www-anp.lip6.fr/cadna/
#  --------------------
AC_ARG_WITH(cadna,[  --with-cadna= cadna library --without-cadna  ])
if test "$with_cadna" != no -o -n  "$with_cadna"
    then
    ff_cadna="$with_cadna"
fi
ff_cadna_ok=no
if test "$with_cadna" != no
    then
    ff_libs_old="$LIBS"
    ff_ldflags_old="$LDFLAGS"
    ff_cadna_dir="${curdir}/download/cadna"
    LDFLAGS="$LDFLAGS -L$ff_cadna_dir"
    AC_CHECK_LIB(cadnafree,arit_zero,
	ff_cadna_ok=yes)
    AC_CHECK_HEADERS(${ff_cadna_dir}/cadnafree.h,
	ff_cadna_h=yes, ff_cadna_h=no)

    LIBS="$ff_libs_old"
    LDFLAGS="$ff_ldflags_old"

    if test "$ff_cadna_ok" = yes  -a "$ff_cadna_h" = yes
	then
	AC_DEFINE(HAVE_CADNA,1, freecadna is use to evalute the round-off error propagation )
	CPPFLAGS="$CPPFLAGS -I$ff_cadna_dir"
	LIBS=" -L$ff_cadna_dir -lcadnafree $ff_libs_old"
    else
	AC_MSG_NOTICE([ without cadna ***** ])
    fi
else
    AC_MSG_NOTICE([ without cadna ***** ])
fi

# Looking for FFTW
# ----------------

# ALH - 18/9/13 - deactivates FFTW detection for testing purposes

AC_ARG_ENABLE(system_fftw,AC_HELP_STRING([--disable-system-fftw],[Disable the automatic detection of FFTW]))
if test "$enable_system_fftw" != no
then

    # ALH - 4/9/13 - avoid recompiling FFTW if a copy was found on the system (request from Helmut Jarausch, 1/8/13)

    AC_CHECK_LIB(fftw3,fftw_execute,ff_fftw_ok=yes)
    AC_CHECK_HEADERS(fftw3.h,ff_fftw_h=yes, ff_fftw_h=no)
fi

if test "$ff_fftw_ok" = yes -a "$ff_fftw_h" = yes
then
    AC_FF_ADDWHERELIB(fftw3,-lfftw3,)
    fftw_download=
else
    if test "$enable_download" = yes
    then
	ff_DOWNLOAD_FFTW=fftw
    fi
fi

# used in [[file:download/fftw/Makefile.am::DOWNLOAD_FFTW]]

AC_SUBST(DOWNLOAD_FFTW,$ff_DOWNLOAD_FFTW)

dnl ---------------------------
dnl   Looking for the tetgen
dnl --------------------------
dnl AC_CHECK_LIB(tet,tetrahedralize,ff_fftet_ok=yes)
dnl AC_CHECK_HEADERS(tetgen.h,ff_tet_h=yes, ff_tet_h=no)
dnl if test "$ff_tet_ok" = yes -a "$ff_tet_h" = yes
dnl then
dnl    AC_SUBST([TOOL_COMPILE_tetgen],"")
dnl AC_FF_ADDWHERELIB(tetgen,-ltet,)
dnl    AC_DEFINE(HAVE_TETGEN,1, tetgen is compute tetrahedralize volume of an enclosed surface)
dnl fi

# Looking for the BLAS
# --------------------

ff_blas_ok=no
ff_blas_inc=""

# ALH - 18/9/13 - give the option to deactivate system blas for testing purposes
AC_ARG_ENABLE(system-blas,AC_HELP_STRING([--disable-system-blas],[Disable the search for a system-wide BLAS library]))
if test "$ff_where_lib_conf_blaslapack" = 1  ; then
  echo " use BLAS/Lapack of petsc "
#  echo "  lib: $ff_where_lib_blaslapack inc: $ff_where_inc_blaslapack"
	    ff_blas_ok=yes
	    ff_lapack_ok=yes
	    ff_blas_libs="$ff_where_lib_blaslapack"
	    ff_lapack_lib=
	    ff_blas_inc="$ff_where_inc_blaslapack"
       AC_FF_ADDWHERELIB(lapack,$ff_mkl_lapack,$ff_blas_inc)
	   AC_FF_ADDWHERELIB(blas,$ff_blas_libs,$ff_blas_inc)
	   AC_MSG_NOTICE([ BLAS and LAPACK  from  petsc ***** ]) 
	   if  echo  "$ff_where_lib_blaslapack"| grep -q  " -lmkl"  ; then
	     	AC_MSG_NOTICE([ MKL  from  petsc ***** ]) 
	        AC_FF_ADDWHERELIB(mkl,$ff_blas_libs,$ff_blas_inc)
	   fi

fi

if test "$enable_system_blas" != no -a "$ff_blas_ok" = no
then

    # User-specified location
    # add MKL  seach dec 2010  FH ....

    ff_mkl_flags=""
    if test "$ff_mkl_libpath" != "not" ; then
	ff_CFLAGS="$CFLAGS"
	CHECK_COMPILE_FLAG(C,-mkl,CFLAGS)
	if test  "$ff_CFLAGS" != "$CFLAGS" ; then
	    ff_ff_mkl_flags="-mkl"
	fi
	CFLAGS="$ff_CFLAGS"
    fi

    AC_MSG_CHECKING(for MKL)
    ff_mkl_root=""
    if test "$ff_mkl_libpath" != "not" -a -z "$ff_mkl_flags" ; then
	
	if test "$ff_mkl_libpath"  = "yes"  ; then
	    ff_IFS="$IFS"
	    IFS=":"
	    for i in $LD_LIBRARY_PATH:$LIBRARY_PATH ; do
		case $i in
		    */mkl/*) ff_mkl_libpath=$i;;
		esac ;
	    done;
	    IFS="$ff_IFS"
#	else
#	    ff_mkl_libpath=no
	fi

	if test "$ff_mkl_libpath" != "no" -a -d "$ff_mkl_libpath" ; then
	    ff_mkl_root=`expr "//$ff_mkl_libpath" : '//\(.*\)/lib.*'`
	    ff_mkl_arch=`expr "//$ff_mkl_libpath" : '//.*/lib/\(.*\)'`
            ff_mkl_lp=_lp64
            case "$ff_mkl_arch" in
		*64*)  ff_mkl_lp=_lp64
	    esac
            case $F77 in
		*ifort*) ff_mkl_cc=intel;;
		*) ff_mkl_cc=gf;;
            esac
            # bof bof ....
	    case "$MPIRUN" in
		*/sgi/*)         ff_mkl_mpi=_sgimpt;;
		*/intel/*)         ff_mkl_mpi=_intelmpi;;
		*)        ff_mkl_mpi=_openmpi;;
	    esac
#            echo ................  ff_mkl_root = $ff_mkl_root  .. $ff_mkl_arch
	    if test ! -d "$ff_mkl_libpath" ; then
		ff_mkl_libpath="$ff_mkl_root/lib/$ff_mkl_arch"
	    fi
	    
	    if test  -f "$ff_mkl_libpath/libmkl_rt.$ff_suffix_dylib" ; then
		mkl_blas=rt
		mkl_lapack=rt
		ff_mkl_mlt="-L$ff_mkl_libpath -lmkl_rt -lmkl_intel_thread -lmkl_core  $ff_iomp5 $ff_pthread"
		ff_mkl_blas="-L$ff_mkl_libpath -lmkl_rt -lmkl_sequential -lmkl_core   $ff_iomp5 $ff_pthread"
		ff_mkl_lapack="$ff_mkl_blas"
	    elif test -f "$ff_mkl_libpath/libmkl_lapack.$ff_suffix_dylib" ; then
		ff_mkl_mlt="-L$ff_mkl_libpath -lmkl_${ff_mkl_cc}${ff_mkl_lp} -lmkl_lapack -lmkl_intel${ff_mkl_lp} -lmkl_intel_thread  -lmkl_core -lguide  -lm -lpthread"
		ff_mkl_blas="-L$ff_mkl_libpath -lmkl_${ff_mkl_cc}${ff_mkl_lp} -lmkl_lapack -lmkl_intel${ff_mkl_lp} -lmkl_sequential -lmkl_core -lguide -lm -lpthread"
		ff_mkl_lapack="-L$ff_mkl_root/lib/$ff_mkl_arch -lmkl_lapack"
	    else
		ff_mkl_mlt="$ff_ldbg  -L$ff_mkl_libpath -lmkl_${ff_mkl_cc}${ff_mkl_lp}  -lmkl_intel${ff_mkl_lp}  -lmkl_intel_thread  -lmkl_core -lmkl_intel_thread $ff_iomp5  $ff_ldeg  -lm $ff_pthread"
		ff_mkl_blas="$ff_ldbg  -L$ff_mkl_libpath -lmkl_${ff_mkl_cc}${ff_mkl_lp}  -lmkl_intel${ff_mkl_lp} -lmkl_sequential -lmkl_core  -lmkl_sequential $ff_ldeg   -lm $ff_pthread"
		ff_mkl_lapack="$ff_mkl_blas"
	    fi
	    if test -f "-L$ff_mkl_libpathmkl_scalapack${ff_mkl_lp}" ; then 
            ff_mkl_scalapack="-L$ff_mkl_libpath -lmkl_blacs${ff_mkl_mpi}${ff_mkl_lp} -lmkl_scalapack${ff_mkl_lp} "
            ff_mkl_blacs="-L$ff_mkl_libpath -lmkl_blacs${ff_mkl_mpi}${ff_mkl_lp} -lmkl_scalapack${ff_mkl_lp}"
	    fi
	    ff_blas_ok=yes
	    ff_lapack_ok=yes

	    # <<enable_mkl_mlt>> ALH - 6/11/13 - request from Atsushi Suzuki - the default MKL library must be able to
	    # handle threads when FF is connected to AS' solver. So we need an option to configure FF with the
	    # multithreaded MKL by default.

	    AC_ARG_ENABLE(mkl_mlt,AC_HELP_STRING([--enable-mkl-mlt],[Link with the multithreaded instead of the monothreaded version of the MKL]))

	    if test "$enable_mkl_mlt" = yes
	    then
		ff_blas_libs="$ff_mkl_mlt"
	    else
		ff_blas_libs="$ff_mkl_blas"
	    fi

	    ff_blas_inc="-I$ff_mkl_root/include"
	    dnl  scalapack, and blacs do not not work with mlk and sgi
	    if test "$ff_mkl_mpi" = "_intelmpi" ; then
		AC_FF_ADDWHERELIB(scalapack,$ff_mkl_scalapack,$ff_blas_inc)
		AC_FF_ADDWHERELIB(blacs,$ff_mkl_blacs,$ff_blas_inc)
	    else
		ff_warm=" (We do not use MKL scalapack and blacs with sgi MPI),  "
	    fi
	    AC_FF_ADDWHERELIB(lapack,$ff_mkl_lapack,$ff_blas_inc)
	    AC_FF_ADDWHERELIB(mkl,$ff_mkl_mlt,$ff_blas_inc)
            AC_DEFINE(HAVE_MKL,1, the MKL intel lib is present for BLAS and LAPACK  )
	else
	    ff_mkl_libpath=no

	# FH - pardiso is there as soon as mkl is
            enable_pardiso=no
	fi
    fi


    AC_MSG_RESULT( [ $ff_warm  root: $ff_mkl_root , arch: $ff_mkl_arch , $ff_mkl_lp ...  ])

    AC_ARG_WITH(blas,
	AC_HELP_STRING([--with-blas=library],[Use a specific version of the Blas]),
	ff_blas_ok=yes
	ff_blas_libs="${withval}")

    # Specific BLAS library location for FreeBSD
    ff_freebsd_blas="-lf2c -lf77blas -latlas -lgslcblas"
    ff_dll_blas=`which  libopenblas.dll`
    if [ -z "$ff_dll_blas" ] ; then  ff_dll_blas="no"; fi
        #  zmach is in lapack
    # Trying  blas library
    # echo "LIBS =  $LIBS blas --- $ff_blas_ok"
    for iblas in "$ff_dll_blas" \
	"-framework Accelerate" "$ff_mkl_blas" "-lblas"  "-L/usr/lib/atlas -lblas" \
	$ff_freebsd_blas
    do
	if test "$ff_blas_ok" = no -a "$iblas" != "no" 
	then
            AC_MSG_CHECKING(for daxpy_ in $iblas)
            ff_save_libs="$LIBS"
            LIBS="$LIBS $iblas"
            AC_LINK_IFELSE(
		[AC_LANG_CALL(,daxpy_)],
		ff_blas_ok=yes
		ff_blas_libs="$iblas",)
            LIBS="$ff_save_libs"
            AC_MSG_RESULT($ff_blas_ok)
	fi
    done

    #cblas_zdotu_sub
    #echo "LIBS = " $LIBS
    ff_cblas_libs=no
    if test "$ff_blas_ok" = yes
    then
	AC_MSG_CHECKING(for blas_zdotu_sub in $ff_blas_libs)
	ff_save_libs="$LIBS"
	LIBS="$LIBS $ff_blas_libs $FLIBS"
	AC_LINK_IFELSE(
            [AC_LANG_CALL(,cblas_zdotu_sub)],
            ff_cblas_libs=yes)
	AC_MSG_RESULT($ff_cblas_libs)
	LIBS="$ff_save_libs"
    fi
fi

# <<OpenBLAS>> ALH - 18/9/13 - option to compile the OpenBLAS moved from the FFCS tree to FF tree

AC_ARG_ENABLE(openblas,AC_HELP_STRING([--disable-openblas],[Disable the automatic download of OpenBLAS]))
if  test "$ff_blas_ok" = no && test "$enable_openblas" != no && test "$enable_download" = yes
then
    AC_CHECK_PROG(ff_git,git,yes,no)
    AC_MSG_CHECKING([Activating the OpenBLAS (deactivate with --disable-openblas)])
    if test $ff_git = no
    then
	AC_MSG_RESULT([git not found])
    else
	ff_download_blas=openblas

	# it is necessary to split the library path into -L and -l otherwise ff/upstream/examples+++-load/ff-c++ will place -lgfortran in
	# LIBS, .../libgoto___.lib in OTHER, and call the compiler with $LIBS before $OTHER, which will fail. But libopenblas.a is placed in
	# the standard directory download/lib so the -L option is not required anymore.

	ff_blas_libs="-lopenblas -lgfortran -lquadmath $ff_pthread"
	ff_blas_inc="-I${curdir}/download/blas/CBLAS/src"

        # skip generic blas compilation and activate OpenBLAS (see [[file:download/blas/Makefile.am::OpenBLAS]])
	AC_SUBST(DOWNLOADED_BLAS,"")
	AC_SUBST(COMPILE_OPENBLAS,openblas)

	ff_blas_ok=yes
	ff_cblas_h=yes
	ff_cblas_libs=yes

	AC_MSG_RESULT([ok])
    fi
fi

# If all else fails, download a generic version

if test "$ff_blas_ok" = no -a "$enable_download" = yes -a "$enable_fortran" != no -a "$enable_c" != no
then
    AC_CHECK_PROG(ff_unzip,unzip,yes,no)

    AC_MSG_CHECKING(for BLAS version to download)

    # Do not update $LIBS, but create an extra LIB variable, because this lib does not exist yet, and this could make
    # the following tests fail.

    # When compiling a generic version, we do not need an optimized version of the BLAS.

    ff_download_blas=generic
    ff_blas_libs="-L${curdir}/download/blas -lcblas -lf77blas"
    ff_blas_inc="-I${curdir}/download/blas/CBLAS/src"

    AC_SUBST(DOWNLOADED_BLAS,"libf77blas.a libcblas.a")
    AC_SUBST(DOWNLOADED_BLAS_BUILT_SOURCES,'BLAS CBLAS $(F77BLAS_SOURCES) $(CBLAS_SOURCES)')

    AC_MSG_RESULT($ff_download_blas)
    ff_blas_ok=yes
    ff_cblas_h=yes
    ff_cblas_libs=yes
fi

AC_SUBST(BLASLIBS,$ff_blas_libs)
AC_SUBST(BLASINC,$ff_blas_inc)

#looking for cblas.h  FH

if test "$ff_blas_ok" = yes -a "$ff_cblas_libs" = yes ;
then
    AC_CHECK_HEADERS(cblas.h,
	ff_cblas_h=yes,
	ff_cblas_h=no)
    AC_CHECK_HEADERS(Accelerate/cblas.h,
	ff_cblas_h=yes
	ff_cblas_h=no)
    AC_CHECK_HEADERS(atlas/cblas.h,
	ff_cblas_h=yes
	ff_cblas_h=no)
fi

# ALH - 18/9/13 - [[file:download/blas/Makefile.am]] does not update the WHERE mechanism so it always needs to be set up
# here even in the case of a downloaded blas library.

AC_FF_ADDWHERELIB(blas,$ff_blas_libs,)

# end of BLAS -------------------

# Looking for ARPACK
# ------------------

# We need the following g77 libraries to connect to the Fortran 77
# Arpack.

if test "$ff_blas_ok" = yes;
then
  ff_g2c_lib="$FLIBS"
       if test  -n "$ff_g2c_lib"
       then
           LIBS="$LIBS $ff_g2c_lib"
       else
	  ff_g2c_ok=no
	  AC_CHECK_LIB(g2c,G77_second_0,
		LIBS="$LIBS -lg2c"
		ff_g2c_lib="-lg2c"
		ff_g2c_ok=yes)

	  if test "$ff_g2c_ok" = no;
	     then
	        AC_MSG_CHECKING(for G77_second_0 in /sw/lib/libg2c.a)
		ff_save_libs="$LIBS"
		LIBS="$LIBS -L/sw/lib -lg2c"
		AC_LINK_IFELSE(
			[AC_LANG_CALL(,G77_second_0)],
			ff_g2c_lib="-L/sw/lib -lg2c"
			ff_g2c_ok=yes,
			LIBS="$ff_save_libs")
		AC_MSG_RESULT($ff_g2c_ok)
	  fi
	fi
fi

# Copy the result of g2c investigations into a separate variable
# because BLAS compilation will need it.

AC_SUBST(G2CLIB,$ff_g2c_lib)
#looking of  lapack if no compile the arpack lapack
# warning $$ because the make eat one
ff_lapack_ok=no;
ff_lapack_lib=
AC_MSG_CHECKING([for lapack in $LIBS, $ff_blas_libs and -llapack] )
    if test "$ff_blas_ok" = yes ; then
        AC_LINK_IFELSE( [AC_LANG_CALL(,dgeqr2_)], ff_lapack_ok=yes)
	if test "$ff_lapack_ok" = no ; then
	    ff_save_libs=$LIBS
	    LIBS="$ff_save_libs $ff_blas_libs"
            AC_LINK_IFELSE( [AC_LANG_CALL(,dgeqr2_)], [ff_lapack_ok=yes;ff_lapack_lib="$ff_blas_libs"])
	    if test "$ff_lapack_ok" = no ; then
	    LIBS="$ff_save_libs -llapack $ff_blas_libs"
            AC_LINK_IFELSE( [AC_LANG_CALL(,dgeqr2_)], [ff_lapack_ok=yes
	    ff_lapack_lib="-llapack"] )
	    fi
	    LIBS="$ff_save_libs"
	fi
    fi

AC_MSG_RESULT($ff_lapack_ok)
if test "$ff_lapack_ok" != no ; then
AC_SUBST(LAPACKLIBS,$ff_lapack_lib)
fi
if test "$ff_lapack_ok" = no; then
    AC_ARG_WITH(lapack,
	AC_HELP_STRING([--with-lapack=library],[Use a specific version of Lapack]),
	ff_lapack_ok=yes
	ff_lapack_lib="${withval}"
	LIBS="$ff_lapack_lib $LIBS")
fi

# Lapack configuration for Arpack

ff_lapackdir='$$(LAPACKdir)'
if test "$ff_lapack_ok" = yes
then
    # no compilation of lapack in arpack
    ff_lapackdir=
else
    ff_lapack_lib="-L${curdir}/download/lib -llapack"
fi

# Arpack itself

ff_arpack_ok=no
ff_save_libs="$LIBS"
ff_arpack_libs=
if test "$ff_blas_ok" = yes;
then

	# User-specified location
	AC_ARG_WITH(arpack,
		[  --with-arpack=library	Use a specific version of Arpack],
		ff_arpack_ok=yes
		ff_arpack_libs="${withval}")

	# Default locations
	if test "$ff_arpack_ok" = no;
	then
		AC_CHECK_LIB(arpack,dsaupd_,
			ff_arpack_libs="-larpack -llapack"
			ff_arpack_ok=yes,,
			-llapack)
	fi

   # Trying to "locate" Arpack
   if test "$ff_arpack_ok" = no -a "$enable_download" != yes ;
   then
        AC_MSG_CHECKING(for libarpack with locate)
        ff_lib_arpack=`locate libarpack|grep  'libarpack.*.a$'|head -1`
	LIBS="$ff_lib_arpack $LIBS"
        AC_LINK_IFELSE(
                [AC_LANG_CALL(,dsaupd_)],
                ff_arpack_ok=yes
                ff_arpack_libs="$ff_lib_arpack")
        AC_MSG_RESULT($ff_arpack_ok)
    fi
  if test "$ff_arpack_ok" = yes
  then
      AC_FF_ADDWHERELIB(arpack,$ff_arpack_libs,)
      echo  arpack LD "'$ff_arpack_libs'"  >>$ff_where_lib_conf
  fi

	# If all else fails, download!
	if test "$ff_arpack_ok" = no -a "$enable_download" = yes \
	   -a "$enable_fortran" != no
	then
                ff_arpack_download=yes
		AC_MSG_NOTICE(using downloaded Arpack)

		# ALH - 6/11/13 - this install goal is the standard goal for all downloaded packages in
		# [[file:download/common.mak::install]]

		AC_SUBST(DOWNLOAD_ARPACK,install)

		AC_SUBST(FF_LAPACKdir,$ff_lapackdir)
		AC_SUBST(ARPACKLIB,${curdir}/download/lib/libarpack.a)
		AC_SUBST(LAPACK_arpack_LIB,${curdir}/download/lib/liblapack.a)

		# Do not update $LIBS, but create an extra LIB variable, because this lib does not exist yet, and this
		# could make the following tests fail.

		# ALH - 30/9/13 - do not use the "-L ${curdir}/download/lib" directive because it would allow other
		# following -l directives (eg -lumfpack) to pick an old locally compiled library instead of the system
		# ones.

                ff_arpack_libs="${curdir}/download/lib/libarpack.a $ff_lapack_lib"
		ff_arpack_ok=yes
	fi
fi
if test  "$ff_arpack_ok" != yes
 then
	    AC_MSG_NOTICE([-- NO ARPACK --  enable_download : $enable_download , wget: $ff_wget ])
fi

# Do not insert ARPACK libs in $LIBS yet, because they may not exist
# yet, and this could make the following tests fail.

LIBS="$ff_save_libs"

if test "$ff_arpack_ok" = yes;
then
        AC_SUBST(ARPACKLIBS,$ff_arpack_libs)
	EIGENOBJ='eigenvalue.$(OBJEXT)'
	AC_DEFINE(HAVE_LIBARPACK,1,Arpack is used for eigenvalue computation)

	# Determines whether to run the eigenvalue tests
else
# all eigen test fails
	AC_SUBST([SKIP_TESTS_EIGEN],[yes])
fi
AC_SUBST([EIGENOBJ])

# ALH - 25/9/13 - <<WHERE_lapack>> always run the WHERE LD statement for lapack because some libraries in
# examples++-load require it (eg [[file:examples++-load/Element_Mixte.cpp::lapack]]). Request from Fred.

AC_FF_ADDWHERELIB(lapack,$ff_lapack_lib,)

# Looking for UMFPACK
# -------------------
ff_umfpack_incs=""
ff_amd_ok=no
ff_umfpack_ok=no
# echo @@@@@@@@@ ff_where_lib_conf_suitesparse == $ff_where_lib_conf_suitesparse  @@@@
if test "$ff_where_lib_conf_suitesparse" = "1"  ; then
  echo "amd/umfpack/suitesparse  of petsc"
  echo "    lib: $ff_where_lib_suitesparse inc: $ff_where_inc_suitesparse"
  ff_amd_ok=yes
   ff_umfpack_ok=yes
   enable_system_umfpack=no
   ff_umfpack_libs="$ff_where_lib_suitesparse"
   ff_umfpack_incs="$ff_where_inc_suitesparse"
   AC_DEFINE(HAVE_LIBUMFPACK,1,Umfpack is used for sparse matrices computations)
   AC_DEFINE(HAVE_LIBCHOLMOD,1,Cholmod  is used for sparse matrices computations)
   AC_SUBST([UMFPACK_CPPFLAGS],[$ff_where_inc_suitesparse])
   AC_DEFINE(HAVE_UMFPACK_H,1,If umfpack.h is located in UMFPACK subdir)
fi

ff_save_libs="$LIBS"
AC_ARG_ENABLE(system_umfpack,AC_HELP_STRING([--disable-system-umfpack],[Disable the automatic detection of umfpack, colmod, amd, ....]))
if test "$ff_umfpack_ok"  = no -a "$enable_system_umfpack" != no -a  "$ff_blas_ok" = yes;
then
    # User-specified location

    AC_ARG_WITH(amd,
	AC_HELP_STRING([--with-amd=library],[Use a specific version of AMD]),
        ff_amd_ok=yes
	ff_umfpack_libs="${withval}")

    AC_ARG_WITH(umfpack,
        AC_HELP_STRING([--with-umfpack=library],[Use a specific version of Umfpack]),
        ff_umfpack_ok=yes
        ff_umfpack_libs="${withval} $ff_umfpack_libs")

    if test "${with_umfpack_include+set}" = set; then 
      	CPPFLAGS="$CPPFLAGS $with_umfpack_include"
    fi

    AC_CHECK_HEADERS(umfpack.h umfpack/umfpack.h ufsparse/umfpack.h suitesparse/umfpack.h,
	ff_umfpack_header=yes
	ff_umfpack_dir=`dirname $ac_header`)

    # Somes systems like FreeBSD hide umfpack.h in a directory called UMFPACK (all capitals). This breaks the
    # standard #define produced by autoconf in config.h.in.

    LIBS="$ff_blas_libs $LIBS"
    if test "$ff_umfpack_header" != yes;
    then
	AC_CHECK_HEADER(UMFPACK/umfpack.h,
	    ff_umfpack_header=yes
	    ff_umfpack_dir=UMFPACK
	    AC_DEFINE(HAVE_BIG_UMFPACK_UMFPACK_H,1,If umfpack.h is located in UMFPACK subdir))
    fi

    if test "$ff_amd_ok" = no;
    then
	AC_CHECK_LIB(amd,amd_info,
	    ff_umfpack_libs="$ff_umfpack_libs -lamd"
	    ff_amd_ok=yes)
    fi

    if test "$ff_umfpack_ok" = no -a "$ff_amd_ok" = yes;
    then

        # ALH - 30/9/13 - other libraries required by Umfpack

	AC_CHECK_LIB(cholmod,cholmod_add,ff_umfpack_libs="$ff_umfpack_libs -lcholmod")
	AC_CHECK_LIB(colamd,colamd_set_defaults,ff_umfpack_libs="$ff_umfpack_libs -lcolamd")

	AC_CHECK_LIB(umfpack,umf_i_malloc,
	    ff_umfpack_libs="-lumfpack $ff_umfpack_libs"
	    ff_umfpack_ok=yes,,$ff_umfpack_libs)
    fi

    if test "$ff_umfpack_header" != yes -o "$ff_umfpack_ok" != yes; then
        ff_umfpack_ok=no
	AC_MSG_WARN([Sorry, we could not find the UMFPACK lib or the UMFPACK headers])
    fi

    if test "$ff_umfpack_ok" = yes -a "$ff_amd_ok" = yes; then
	AC_DEFINE(HAVE_LIBUMFPACK,1,Umfpack is used for sparse matrices computations)
    fi
fi

LIBS="$ff_save_libs"

# If all else fails, download!
#  attention no /usr/include in WHERE 
if test "$ff_umfpack_ok" = yes -a -n "$ff_umfpack_dir" -a -z "$ff_umfpack_incs"
then
	ff_umfpack_incs="-I/usr/include/$ff_umfpack_dir"	
fi
if test "$ff_umfpack_ok" = no -a "$enable_download" = yes
then
    AC_MSG_NOTICE(using downloaded UMFPACK)
    AC_SUBST(DOWNLOAD_UMFPACK,"umfpack")
    ff_umfpack_download=yes

    # Do not update $LIBS, but create an extra LIB variable, because this lib does not exist yet, and this could make
    # the following tests fail.
    
    ff_umfpack_libs="-L${curdir}/download/lib -lumfpack -lcholmod -lcolamd -lamd -lsuitesparseconfig"
    ff_umfpack_incs="-I${curdir}/download/include"
    AC_DEFINE(HAVE_LIBUMFPACK,1,UMFPACK)
    if test "$ff_win32" = yes; then
  	AC_SUBST(FF_UMFPACK_CONFIG,-DCBLAS)
    fi
    ff_umfpack_ok=yes
fi

# ALH - 17/9/13 - moved UMFPACK configuration settings in wherelib to _after_ configuring the download version because
# [[file:download/umfpack/Makefile.am]] does not set the WHERE mechanism. Also removed -I/usr/include/$ff_umfpack_dir
# from include options for FFCS because it breaks the MingW64 compilation process.

if test "$ff_umfpack_ok" = yes
then
    if test $enable_ffcs = no
    then
	AC_FF_ADDWHERELIB(amd,$ff_umfpack_libs,$ff_umfpack_incs)
	AC_FF_ADDWHERELIB(umfpack,$ff_umfpack_libs,$ff_umfpack_incs)
    else
	AC_FF_ADDWHERELIB(amd,$ff_umfpack_libs,)
	AC_FF_ADDWHERELIB(umfpack,$ff_umfpack_libs,)
    fi
else
    AC_MSG_NOTICE([ -- NO UMFPACK (ff_wget = $ff_wget)])
fi

AC_SUBST(UMFPACKLIBS,$ff_umfpack_libs)

# If times() and sysconf() are not here, UMFPACK should know

if test "$ff_umfpack_ok" = yes
then
    AC_CHECK_FUNCS(times sysconf,
	ff_umfpack_posix_ok=yes,
	ff_umfpack_posix_ok=no)

    if test "$ff_umfpack_posix_ok" = no
    then
	AC_SUBST(FF_UMFPACK_CONFIG,"-DCBLAS -DNPOSIX")
    fi
fi

# Checking for some functions that may not appear everywhere
# ----------------------------------------------------------
AC_CHECK_HEADERS(cstddef)
AC_CHECK_HEADERS(stddef.h)

AC_CHECK_HEADERS(regex.h,ff_regex_h=1, ff_regex_h=0)
AC_SUBST(FF_HAVE_REGEX_H,$ff_regex_h)

AC_CHECK_HEADERS([unistd.h])
# asinh acosh atanh are not in Mingw yet gettimeofday
ff_malloc_h=""
AC_HEADER_TIME
AC_CHECK_HEADERS(malloc.h,ff_malloc_h=1)
AC_SUBST(FF_MALLOC_H,$ff_malloc_h)

AC_CHECK_FUNCS(asinh acosh atanh getenv jn erfc tgamma gettimeofday mallinfo mstats)
AC_CHECK_FUNCS(srandomdev)


AC_CHECK_FUNCS(second_,ff_second="",ff_second=second.o)
AC_SUBST(FF_SECOND,"$ff_second")

# Enable static linking (no shared libraries)
# -------------------------------------------
AC_CHECK_PROG(ff_libtool,libtool,yes,no)
if test  "$ff_mac" = "yes"  -a "$ff_libtool" = yes ; then
  ff_AR="libtool"
  ff_ARFLAGS="-static -o"
  ff_RANLIB="echo"
fi

AC_ARG_ENABLE(static,
	[  --enable-static	Build binaries with no shared library dependencies])
if test "$enable_static" = yes
then
	AC_CHECK_PROG(ff_libtool,libtool,yes,no)
	if test "$ff_libtool" = yes
	then
		LDFLAGS="$LDFLAGS -all-static"
		AC_SUBST(STATICTOOL,libtool)
	else
		AC_MSG_ERROR(libtool not found)
	fi
fi
#  for compiation of  arpack  use libtool to bluid universal library on MacOs.
AC_SUBST(AR,$ff_AR)
AC_SUBST(ARFLAGS,$ff_ARFLAGS)
AC_SUBST(RANLIN,$ff_RANLIB)


# Dynamic loading of compiled functions
# -------------------------------------

# Not if we don't want shared libraries (non FH  modif FH juin 2005)
ff_dynload=no
if test "$enable_static" != yes
then

	# Availability of dlopen(). Use AC_COMPILE rather than
	# AC_CHECK_HEADERS because the latter has problems seeing it (in
	# Cygwin) when it does not compile (in Mingw).

	AC_MSG_CHECKING(for dlfcn.h)
	AC_COMPILE_IFELSE([AC_LANG_SOURCE([[#include <dlfcn.h>]])],
		ff_dynload=yes,
		ff_dynload=no)
	AC_MSG_RESULT($ff_dynload)
fi

# FFCS - -lm missing for ffmedit link stage on Debian Testing
AC_CHECK_LIB(m,sin)

# Checks that we also have the corresponding library
if test "$ff_dynload" = yes
then
	AC_CHECK_LIB(dl,dlinfo)

	# Checks that everythings works ok
	AC_MSG_CHECKING(whether dlopen links ok)
	AC_LINK_IFELSE(
[AC_LANG_SOURCE([[#include <dlfcn.h>
int main(int argc,char **argv){
  dlopen("",RTLD_LAZY);
  return 0;
}]])],
	ff_dynload=yes,
	ff_dynload=no)
	AC_MSG_RESULT($ff_dynload)
fi

#  the -rdynamic don't exist on macos and sunOS
if test "$ff_dynload" = yes
then
	AC_DEFINE(HAVE_DLFCN_H,1,Dynamic loading - not mandatory)
	# Activate dynamic loading tests (see examples++-load/Makefile.am)
	AC_SUBST(LOAD_TESTS,../regtests.sh)
	AC_SUBST(LOAD_COMPILE,load_compile)

	# gcc on MacOS does not produce an error with "-rdynamic" but
	# still complains about it.
	if test "$ff_mac" = "no" -a "$ff_win32"  = "no" -a "$ff_sunos"  = "no" ;
	then
		CHECK_COMPILE_FLAG(C++,-rdynamic,LDFLAGS)
dnl		CHECK_COMPILE_FLAG(C,-rdynamic,CNOFLAGS)
dnl		CHECK_COMPILE_FLAG(C,-rdynamic,CFLAGS)
	fi
	if test "$ff_fpic" != "no" ; then
	CHECK_COMPILE_FLAG(C++,-fPIC,CXXFLAGS)
	CHECK_COMPILE_FLAG(C,-fPIC,CFLAGS)
	CHECK_COMPILE_FLAG(C,-fPIC,CNOFLAGS)
	fi 
	if test "$enable_fortran" != no ;  then
	if test "$ff_fpic" != "no" ; then
	 CHECK_COMPILE_FLAG(Fortran,-fPIC,FFLAGS)
	 CHECK_COMPILE_FLAG(Fortran,-fPIC,FNOFLAGS)
	 CHECK_COMPILE_FLAG(Fortran,-fPIC,FCFLAGS)
	 CHECK_COMPILE_FLAG(Fortran,-fPIC,FCNOFLAGS)
	 fi
dnl	 CHECK_COMPILE_FLAG(Fortran,-rdynamic,FNOFLAGS)
dnl	 CHECK_COMPILE_FLAG(Fortran,-rdynamic,FFLAGS)

	fi
fi
AC_SUBST(DYLIB_SUFFIX,$ff_suffix_dylib)

# Checking wether we can generate some documentation
# --------------------------------------------------

AC_CHECK_PROG(ff_latex,latex,yes,no)
AC_CHECK_PROG(ff_makeindex,makeindex,yes,no)
AC_CHECK_PROG(ff_dvips,dvips,yes,no)
# to translate the figure
AC_CHECK_PROG(ff_pdf2ps,pdf2ps,yes,no)
AC_CHECK_PROGS(EPSTOPDF,[epstopdf pstopdf],[false])
AC_CHECK_PROG(ff_convert,convert,yes,no)
if test "$ff_latex" = yes -a "$ff_makeindex" = yes -a "$ff_dvips" = yes -a $ff_pdf2ps = yes -a $ff_convert = yes;
then
	AC_SUBST(DOCPS,"freefem++doc.ps")

	AC_CHECK_PROG(ff_gzip,gzip,yes,no)
	if test "$ff_gzip" = yes;
	then
		AC_SUBST(DOCPSGZ,"freefem++doc.ps.gz")
	fi
fi

# PDF documentation building sometimes poses problems because of pdfsync.sty. So we need to be able to disable it.
# FFCS - 9/4/14 - need to disable PDF also during FFCS build (which does not use the resulting PDF).

AC_ARG_ENABLE(pdf,[  --disable-pdf	Disable PDF documentation building])
if test "$enable_pdf" != no
then
   AC_CHECK_PROG(ff_pdflatex,pdflatex,yes,no)
   if test "$ff_pdflatex" = yes -a $EPSTOPDF != false -a $ff_convert = yes;
   then
	AC_SUBST(DOCPDF,"freefem++doc.pdf")
   fi
fi

# Choosing compilation options for the standard version (in src/std)
# ------------------------------------------------------------------

# The "standard" configured version can use win32 (mingw)

if test "$ff_mingw" = yes
then

    # FFCS does not use FreeFem++-std, and Pcrgraph.cpp does not compile under mingwin64
    if test $enable_ffcs = no
    then
	ff_stdprog="FreeFem++-std${EXEEXT}"
	ff_std_graph_obj=Pcrgraph.$OBJEXT
    fi

    # ALH - FFCS - 30/11/8 - I need to get the output from FF for FFCS regression tests
    if test $enable_ffcs = yes
    then
	ff_std_ldflags="-mconsole -mwindows"
    else
	ff_std_ldflags=-mwindows
    fi

	ff_std_libs=
fi
AC_SUBST(STD_GRAPH_OBJ,$ff_std_graph_obj)
AC_SUBST(STD_LDFLAGS,$ff_std_ldflags)

AC_SUBST(STD_LIBS,$ff_std_libs)

# Allow some downloaded tools not to be compiled
# ----------------------------------------------


##  try to see pakage is hon computer  if the FH ZZZZ 
## FH to find gsl ... 
AX_PATH_GSL(1.15, ff_with_gsl=yes, ff_with_gsl=no)
if test "$ff_with_gsl" = "yes"; then  
  AC_FF_ADDWHERELIB(gsl,$GSL_LIBS,$GSL_CFLAGS)
fi

##


m4_map([AC_FF_WHERELIB],[
	[[mumps],[-ldmumps -lzmumps  -lmumps_common  -lpord],[dmumps_c.h],[]],	
	[[mumps_seq],[-ldmumps_seq -lzmumps_seq  -lmumps_common_seq  -lpord_seq -lmpiseq_seq],[dmumps_c.h],[]],
	[[libseq],[-lmpiseq_seq],[/usr/include/mumps_seq/mpi.h]],
	[[mumps_ptscotch],[-lpord_ptscotch -lmumps_common_ptscotch -ldmumps_ptscotch -lzmumps_ptscotch -lpord_ptscotch],[dmumps_c.h]],
	[[mumps_scotch],[-lpord_scotch -lmumps_common_scotch -ldmumps_scotch -lzmumps_scotch -lpord_scotch],[dmumps_c.h]],
	[[hypre],[-lHYPRE]],
	[[fftw3],[-lfftw3],[/usr/include/fftw3.h],[]],
	[[superlu_dist],[-lsuperlu-dist],[/usr/include/superlu-dist/superlu_defs.h],[]],
	[[superlu],[-lsuperlu],[/usr/include/superlu/superlu_enum_consts.h],[]],
	[[superlu4],[-lsuperlu4],[/usr/include/superlu4/superlu_enum_consts.h],[]],
        [[blacs],[ -lblacsCinit$ff_with_mpi	-lblacsF77init$ff_with_mpi	-lblacs$ff_with_mpi],[]],
        [[scalapack],[-lscalapack$ff_with_mpi],[]],
        [[scotch],[-lscotch -lscotcherr],[scotch.h]],
        [[ptscotch],[-lptscotch -lptscotcherr],[ptscotch.h]],
	[[metis],[-lmetis],[/usr/include/metis/metis.h],[]],
	[[metis],[-lmetis],[metis.h],[]],
	[[parmetis],[-lparmetis -lmetis],[],[]],
	[[freeyams],[-lfreeyams],[freeyamslib.h],[]],
	[[mmg3d],[-lmmg3d],[libmmg3d.h],[]],
	[[mshmet],[-lmshmet],[],[]],
dnl	[[gsl],[-lgsl -lgslcblas -lm],[gsl/gsl_sf.h],[]],
	[[parms],[-lparms -litsol -llapack -lblas -lm],[],[]],
	[[tetgen],[-ltet],[tetgen.h],[]],
	[[ipopt],[ -lipopt],[/usr/include/coin/IpTNLP.hpp],[]],
	[[nlopt],[ -lnlopt],[nlopt.hpp],[]]
    ]
    )
## before try 

TOOL_DISABLE(tetgen,tetgen.$DYLIB_SUFFIX)
TOOL_DISABLE(hips,hips_FreeFem.$DYLIB_SUFFIX)
TOOL_DISABLE(ipopt,ff-Ipopt.$DYLIB_SUFFIX)
TOOL_DISABLE(lapack,"lapack.$DYLIB_SUFFIX fflapack.$DYLIB_SUFFIX")
TOOL_DISABLE(parmetis,metis.$DYLIB_SUFFIX)
TOOL_DISABLE(mmg3d,mmg3d-v4.0.$DYLIB_SUFFIX)
TOOL_DISABLE(mshmet,"mshmet.$DYLIB_SUFFIX aniso.$DYLIB_SUFFIX")
TOOL_DISABLE(gmm,ilut.$DYLIB_SUFFIX)
TOOL_DISABLE(mumps,"MUMPS_FreeFem.$DYLIB_SUFFIX MUMPS.$DYLIB_SUFFIX")
TOOL_DISABLE(mumps_seq,"MUMPS_seq.$DYLIB_SUFFIX MUMPS.$DYLIB_SUFFIX",mumps-seq)
TOOL_DISABLE(nlopt,ff-NLopt.$DYLIB_SUFFIX)
dnl TOOL_DISABLE(parmetis)
TOOL_DISABLE(parms,parms_FreeFem.$DYLIB_SUFFIX)
TOOL_DISABLE(pastix,"interfacepastix.$DYLIB_SUFFIX complex_pastix_FreeFem.$DYLIB_SUFFIX real_pastix_FreeFem.$DYLIB_SUFFIX")
TOOL_DISABLE(pipe,pipe.$DYLIB_SUFFIX)
TOOL_DISABLE(scotch,scotch.$DYLIB_SUFFIX)
TOOL_DISABLE(superlu4,SuperLu.$DYLIB_SUFFIX,superlu)
TOOL_DISABLE(superludist,"complex_SuperLU_DIST_FreeFem.$DYLIB_SUFFIX real_SuperLU_DIST_FreeFem.$DYLIB_SUFFIX dSuperLU_DIST.$DYLIB_SUFFIX")
TOOL_DISABLE(umfpack,UMFPACK64.$DYLIB_SUFFIX)
TOOL_DISABLE(yams,freeyams.$DYLIB_SUFFIX)
TOOL_DISABLE(pipe,pipe.$DYLIB_SUFFIX)

# FFCS - MUMPS_seq has a different Win32 compiler setup from FFCS, so we need to add some extra parameters

if test "$OS" = Windows_NT
then
    CFLAGS="$CFLAGS -DWITHOUT_PTHREAD -DAdd_"

    # we also need to satisfy ff-c++ that the pthread are not a blocking point
    if test -n "$ff_pthread"  ; then
       AC_FF_ADDWHERELIB(pthread,"",)
    fi
fi

# ALH - pARMS needs "-fno-range-check" on Windows, but this options fails on MacOS 10.8.  Add no-range-check for Windows
# for hexadecimal parameter constants like:
#
# [[file:c:/cygwin/home/alh/ffcs/rel/mingw/mpif.h::PARAMETER MPI_SHORT_INT z 8c000003]]
#
# Such constants are rejected without [[file:download/parms/makefile-parms.in::NO_RANGE_CHECK]]

if test "$OS" = Windows_NT
then
    AC_SUBST(NO_RANGE_CHECK,-fno-range-check)
fi

# ALH - 4/9/13 - request from Helmut Jarausch - allow to change Scotch include path
if test "$with_scotch_include" = "" 
then
    with_scotch_include=$ac_pwd/download/include/scotch
fi
AC_SUBST(SCOTCH_INCLUDE,$with_scotch_include)

# Find out kernel and libc versions
# ---------------------------------

if test "$ff_win32" != yes -a "$ff_mac" != yes
then
	AC_MSG_CHECKING(kernel version)
	ff_kernel_version=`cat /proc/version|perl -e '<STDIN>=~/(\d+\.\d+\.\d+)/;print $1;'`
	AC_MSG_RESULT($ff_kernel_version)
	AC_SUBST(KERNEL_VERSION,$ff_kernel_version)

	AC_MSG_CHECKING(libc version)
	ff_libc_version=`ldd /bin/sh | awk '/libc/{print $3}' | xargs readlink | sed -e 's/\.so$//'`
	AC_MSG_RESULT($ff_libc_version)
	AC_SUBST(LIBC_VERSION,$ff_libc_version)
fi
#  def variable pour les makefiles

# creating all makefiles
# ----------------------

ff_bamgprog="bamg${EXEEXT} cvmsh2${EXEEXT}"

AC_SUBST(CNOFLAGS,$CNOFLAGS) dnl for superludist CFLAGS without optim  ...
AC_SUBST(FNOFLAGS,$FNOFLAGS) dnl for blacs CFLAGS without optim  ...

# The final list of executable programs

AC_SUBST(MEDITPROG,$ff_meditprog)
AC_SUBST(FFGLUTPROG,$ff_ffglutprog)
AC_SUBST(BAMGPROG,$ff_bamgprog)
AC_SUBST(STDPROG,$ff_stdprog)
ff_progs="FreeFem++-nw $ff_bamgprog  $ff_mpiprog $ff_meditprog $ff_ffglutprog"

AC_SUBST(TEST_FFPP,$ff_TEST_FFPP)
AC_SUBST(TEST_FFPP_MPI,$ff_TEST_FFPP_MPI)
AC_SUBST(TEST_FFPPMPI,"../src/mpi/ff-mpirun")

ff_with_mpi=-openmpi;

ff_blacs="-lblacsCinit$ff_with_mpi     -lblacsF77init$ff_with_mpi      -lblacs$ff_with_mpi"
ff_scalapack=-lscalapack

# change MKL interface ...
test -n "$ff_mkl_blacs" && ff_blacs="$ff_mkl_blacs"
test -n "$ff_mkl_scalapack" && ff_scalapack="$ff_mkl_scalapack"
test -n "$ff_mkl_root" && ff_winc="$ff_mkl_root/include/mkl_blas.h"

# FFCS - 27/10/11 - Some extra conditionals for things that do not work on certain systems (eg MPI libraries under
# Windows)

##  search of HDF5 ....
AX_LIB_HDF5()
if test "$with_hdf5" = "yes"; then  
  AC_FF_ADDWHERELIB(hdf5,$HDF5_LDFLAGS  $HDF5_LIBS,$HDF5_CPPFLAGS)
  LIBS="$LIBS $HDF5_LDFLAGS  $HDF5_LIBS"
else
  enable_iohdf5=no  
fi
TOOL_DISABLE(iohdf5,iohdf5.$DYLIB_SUFFIX)


AM_CONDITIONAL([FFCS_WINDOWS],[test "$OS" = Windows_NT])

#  remove gsl if not find ... FH
# correction FH .. 18/12/2013. 
# ALH - 7/1/14 - not able to compile gsl or pardiso with FFCS on Windows
if test "$enable_ffcs" != yes
then
    if test  "$ff_where_lib_conf_gsl" = 1  -a  "$enable_gsl" != "no" ; then enable_gsl=yes; fi;
    if test  "$ff_where_lib_conf_mkl" = 1  -a  "$enable_mkl" != "no" ; then enable_pardiso=yes; fi;
fi
TOOL_DISABLE(gsl,"gsl.$DYLIB_SUFFIX")
TOOL_DISABLE(NewSolver,"NewSolver.$DYLIB_SUFFIX")
TOOL_DISABLE(pardiso,"PARDISO.$DYLIB_SUFFIX")

 echo "debug cxxx11: $ff_cxx11 mpi: $ff_mpi petsc: $ff_petsc_ok "
test "$enable_hpddm" != no && enable_hpddm=yes  
test "$ff_cxx11" != yes -o  "$ff_mpi" != yes  -o "$ff_umfpack_ok" != "yes"  &&  enable_hpddm=no
FF_PETCS_DYLIB=""
if test "$ff_petsc_ok" != no ;   then
FF_PETCS_DYLIB="PETSc.$DYLIB_SUFFIX"
fi
##echo @@@@@@@@@@@@ $ff_petsc_ok @@ "hpddm.$DYLIB_SUFFIX $FF_PETCS_DYLIB" 
TOOL_DISABLE(hpddm,["hpddm.$DYLIB_SUFFIX hpddm_substructuring.$DYLIB_SUFFIX $FF_PETCS_DYLIB"])


# All makefiles
AC_OUTPUT(Makefile
    download/Makefile
    download/blas/Makefile
    download/arpack/Makefile
    download/umfpack/Makefile
    download/fftw/Makefile
    src/Makefile
    src/bamglib/Makefile
    src/Graphics/Makefile
    src/femlib/Makefile
    src/Algo/Makefile
    src/lglib/Makefile
    src/fflib/Makefile
    src/nw/Makefile
    src/mpi/Makefile
    src/bamg/Makefile
    src/libMesh/Makefile
    src/medit/Makefile
    src/bin-win32/Makefile
    examples++-load/Makefile
    examples++-tutorial/Makefile
    examples++/Makefile
    examples++-mpi/Makefile
    examples++-hpddm/Makefile
    examples++-eigen/Makefile
    examples++-chapt3/Makefile
    examples++-bug/Makefile
    examples++-other/Makefile
    examples++-3d/Makefile
    DOC/Makefile)

AC_MSG_NOTICE([  freefem++ used  download : $enable_download ])
AC_MSG_NOTICE([	 	 --  Dynamic load facility: $ff_dynload ])
AC_MSG_NOTICE([		 --  ARPACK (eigen value): $ff_arpack_ok ])
AC_MSG_NOTICE([		 --  UMFPACK (sparse solver) $ff_umfpack_ok ])
AC_MSG_NOTICE([		 --  BLAS $ff_blas_ok ])
AC_MSG_NOTICE([		 --  with MPI             $ff_mpi])
AC_MSG_NOTICE([          --  with PETSC $ff_petsc_ok ]);
AC_MSG_NOTICE([          --  with hpddm  $enable_hpddm (need MPI & c++11 :  $ff_cxx11 ) ]);


AC_MSG_NOTICE([    progs: $ff_progs ])
if test  "$ff_umfpack_download" = yes -a  -n "$TOOL_COMPILE_umfpack" ;then
    AC_MSG_NOTICE([      use of download UMFPACK see  download/umfpack/SuiteSparse/UMFPACK/README.txt for the License])
fi
if test "$ff_arpack_download" = yes -a  -n "$TOOL_COMPILE_arpack"  ; then
    AC_MSG_NOTICE([      use of download ARPACK see  download/arpack/ARPACK/README  no License ])
fi

if test "$enable_download" = yes  ;then
    if test  -n "$TOOL_COMPILE_fftw3"   ; then 
    AC_MSG_NOTICE([      use of download fftw   see  download/fftw/fftw-3.2/COPYRIGHT ])
    fi 
    if test  -n "$TOOL_COMPILE_tetgen"   ; then 
    AC_MSG_NOTICE([      use of download tetgen  see download/tetgen/tetgen1.5.1-beta/LICENSE ])
    fi
    if test  -n "$TOOL_COMPILE_yams"   ; then 
    AC_MSG_NOTICE([      use of download freeyams  see download/yams/freeyams.2011.02.22 (suface mesh adaptation)  ])
    fi
    if test  -n "$TOOL_COMPILE_mmg3d"   ; then     
    AC_MSG_NOTICE([      use of download mmg3d (v4)   see download/mmg3d/mmg3d4/LICENCE.txt  ])
    fi
    if test  -n "$TOOL_COMPILE_parmets"   ; then     
    AC_MSG_NOTICE([      use of download parmetis])
     fi
    # ALH - some tools may be activated but not downloaded if a local version is specified (see
    # [[TOOL_PARAMETERS]])

    if test "$enable_superlu_download" = yes
    then
	AC_MSG_NOTICE([      use of download superlu])
    fi

    if test -n "$MPICC" ; then
   	AC_MSG_NOTICE([      try to download: blacs parmetis  scalapack scotch superludist pastix hypre hips library])
    fi

    if test "$ff_download_blas" = "generic"  ; then
        AC_MSG_NOTICE([      use of download generic blas and cblas freefem may be slow  ])
        AC_MSG_NOTICE([      you can try to use the  Kazushige Goto s BLAS at http://www.cs.utexas.edu/users/flame/goto/  ])
        AC_MSG_NOTICE([          or at  http://www.tacc.utexas.edu/~kgoto/ for the best BLAS .])
        AC_MSG_NOTICE([      or try to download/compile the altas blas .])
    fi

    # ALH - 30/10/13 - (request from FH) all downloads are now part of a separate script ([[file:download/getall]])
    # which needs to be run before make for the user to specify whether he is ok with each of the package licences.

    AC_MSG_NOTICE([      Please run "download/getall" to download all necessary packages before running make])
fi

# FFCS does not use glut, so remove this message because it could make the user think that something is broken

if test $enable_ffcs = no
then
    if test "$ff_glut_ok" != yes  ; then
	AC_MSG_NOTICE([     *********************************************** ])
	AC_MSG_NOTICE([     WARNING: you do not have the new grachics tools ])
	AC_MSG_NOTICE([         because the configure do not find OpenGL, GLUT or pthread  developer stuff ])
	AC_MSG_NOTICE([         read the README to find missing package  ])
	AC_MSG_NOTICE([         F. Hecht  ])
	AC_MSG_NOTICE([ to install missing package  under debian or ubuntu, try: sudo apt-get install freeglut3-dev ])
	AC_MSG_NOTICE([     *********************************************** ])
    fi
fi

# Local Variables:
# mode:shell-script
# ispell-local-dictionary:"british"
# coding:utf-8
# End:
