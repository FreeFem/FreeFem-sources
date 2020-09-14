# ------------------------------------------------------------

# Antoine Le Hyaric - LJLL Paris 6 - lehyaric@ann.jussieu.fr - 13/5/04
# $Id: configure.ac,v 1.338 2010/05/10 21:35:07 hecht Exp $

# Version numbering: x.xx-pp where "pp" is the package version (when
# the same FreeFEM version is packaged several times). For
# coherency, this should have the same value as the top-most package
# number in debian/changelog.


AC_INIT(FreeFEM,4.7,frederic.hecht@sorbonne-universite.fr,FreeFEM)
dnl : ${CFLAGS=""}
dnl : ${CXXFLAGS=""}
dnl : ${FCFLAGS=""}
dnl : ${FFLAGS=""}

AC_CONFIG_MACRO_DIR([etc/config/m4])

# progg test  exec for windows
ff_TEST_FFPP="../../src/nw/FreeFem++"
ff_TEST_FFPP_MPI="../../src/mpi/FreeFem++-mpi"
ff_uncompile_plugin=""
ff_unlib_plugin=""
ff_download_lib=""
AC_PREREQ(2.63) dnl for AC_LANG_CASE and check stuff
dnl AC_CONFIG_SRCDIR(src/FreeFem++-CoCoa)

# Automake 1.11 is too old for check ...

# ALH - 10/9/13 - FFCS needs to skip this automake version requirement when compiling on older platforms (eg cygwin or
# last Ubuntu LTS). So we fix a low hardcoded requirement and test for higher versions when not in FFCS.

AM_INIT_AUTOMAKE(1.11 dist-zip foreign subdir-objects)

if test "$enable_ffcs" != yes
then
    if test "$am__api_version" \< 1.13 ;   then
	AC_MSG_ERROR([Automake version $am__api_version needs to be 1.13 or later to enable "make check"])
    fi
fi

AC_CONFIG_HEADERS(config.h)
AC_PROG_MAKE_SET
dnl AM_COLOR_TESTS=always
AC_PROG_RANLIB
inc_usr_include=""
ff_where_lib_conf=plugin/seq/WHERE_LIBRARY-config
dnl search of associad software
m4_define([AC_FF_ADDWHERELIB],
  [
    if test -z "$ff_where_lib_conf_$1"  ; then
	    echo "$1 LD $2"  >>$ff_where_lib_conf
	    test -n "$3" && echo "$1 INCLUDE $3 "  >>$ff_where_lib_conf
            ff_where_lib_conf_$1=1
	    ff_where_lib_$1="$2"
	    ff_where_inc_$1="$3"
	    AC_MSG_NOTICE([    ++ add $1 : $2 ;  $3 in  $ff_where_lib_conf "])
	else
	    AC_MSG_NOTICE([    -- do not add $1 : $2 $3 in  $ff_where_lib_conf "])
        fi
]) dnl
m4_define([AC_FF_WHERELIB],
    [  	if test -z "$ff_where_lib_conf_$1" -a "$enable_download_$1"  != "yes" ; then
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
define([AC_FF_PETSC_AWK],[awk -F'=' '"$1 "==[$]1 {print [$]2}' $2])
m4_define([AC_FF_PETSC_WHERELIB],
	  [
	    ff_$2=`AC_FF_PETSC_AWK($2,$4)`
	    ff_$3=`AC_FF_PETSC_AWK($3,$4)`
	   if test -n "[$]ff_$2"; then
			    ff_$1_ok=yes;
			    m4_ifval($5,ff_$2="$ff_$2 $5";
			    echo " **** add $5 to ff_$2 ******", ) 
			    AC_FF_ADDWHERELIB($1,[$]ff_$2,[$]ff_$3)
			    AC_SUBST([TOOL_COMPILE_$1],"")
	  fi
	  ]
)dnl end m4_defineAC_FF_PETSC_WHERELIB
m4_define([AC_FF_PETSC_WHERELIB_BASIC],# PETSC_WITH_EXTERNAL_LIB,PETSC_CC_INCLUDES
	  [
	    ff_PETSC_EXTERNAL_LIB_BASIC=`AC_FF_PETSC_AWK(PETSC_EXTERNAL_LIB_BASIC,$2)`
	    ff_PETSC_LIB_BASIC=`AC_FF_PETSC_AWK(PETSC_LIB_BASIC,$2)`
	    ff_PETSC_PREFIXDIR=`awk -F' ' '"PREFIXDIR"==[$]1 {print [$]3}' $2`
	    ff_PETSC_WITH_EXTERNAL_LIB="-Wl,-rpath,[$]ff_PETSC_PREFIXDIR/lib -L[$]ff_PETSC_PREFIXDIR/lib [$]ff_PETSC_LIB_BASIC [$]ff_PETSC_EXTERNAL_LIB_BASIC"
	    ff_PETSC_CC_INCLUDES=`AC_FF_PETSC_AWK(PETSC_CC_INCLUDES,$2)`
	   if test -n "[$]ff_PETSC_WITH_EXTERNAL_LIB"; then
			    ff_$1_ok=yes;
			    AC_FF_ADDWHERELIB($1,[$]ff_PETSC_WITH_EXTERNAL_LIB,[$]ff_PETSC_CC_INCLUDES)
			    AC_SUBST([TOOL_COMPILE_$1],"")
	  fi
	  ]
)dnl end m4_defineAC_FF_PETSC_WHERELIB_BASIC

# ALH - this is required by FFCS that needs to deactivate some tools that do not work on all platforms. Some FF users
# may also find interesting to specify a local version of a tool instead of downloading it.

# m4 macro parameters: $1 = tool name, $2 = dynamic library name, $3 = download directory name if different from $1

m4_define([TOOL_PARAMETERS],
    [
	AC_ARG_WITH($1_include,AC_HELP_STRING([--with-$1-include=],[Include directives for $1 instead of automatic download]))
    AC_ARG_WITH($1_ldflags,AC_HELP_STRING([--with-$1-ldflags=],[Link-time directives for $1 instead of automatic download]))
	AC_ARG_ENABLE(download-$1,AC_HELP_STRING([--enable-download_$1],[force the download of $1]))

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
		ff_uncompile_plugin+=" "$2
		ff_unlib_plugin+=" $1"
    elif  test   "$ff_$1_ok" = "yes" ; then
        AC_SUBST([TOOL_DYLIB_$1],$2)
        enable_$1_download=no
    else
	   if test "$enable_download" != "yes" -a "$enable_download_$1" != "yes" ; then
	    AC_SUBST([TOOL_COMPILE_$1],"")
		ff_uncompile_plugin+=" "$2
		ff_unlib_plugin+=" $1"
        AC_SUBST([TOOL_DYLIB_$1],"")
	   else
        if test "$ff_wget" = yes -a "$enable_download_$1" = "yes" -a "$DOWNLOADCOMPILE" = "" ; then
            DOWNLOADCOMPILE="compile-pkg install-other"
        fi
        AC_SUBST([TOOL_COMPILE_$1],ifelse($3,,$1,$3))
        AC_SUBST([TOOL_DYLIB_$1],$2)
		ifelse($4,,[ff_download_lib="$1 $ff_download_lib"])
	    ifelse($4,,[enable_$1_download=yes])		
	   fi
     fi

    # Also allow to disable the download of one tool if it is already locally installed

    TOOL_PARAMETERS($1,$2,$3)
])
m4_define([TOOL_DISABLE_NO],
    [
	AC_ARG_ENABLE($1,AC_HELP_STRING([--disable-$1],[Do not use $1]))
    if test "$enable_$1" = "no" -o "$enable_$1" = ""
    then
	    enable_$1=no
        AC_SUBST([TOOL_COMPILE_$1],"")
        AC_SUBST([TOOL_DYLIB_$1],"")
    elif  test   "$ff_$1_ok" = "yes" ; then
        AC_SUBST([TOOL_DYLIB_$1],$2)
        enable_$1_download=no
    else
        AC_SUBST([TOOL_COMPILE_$1],ifelse($3,,$1,$3))
        AC_SUBST([TOOL_DYLIB_$1],$2)
		ifelse($4,,[ff_download_lib="$1 $ff_download_lib"])
	    ifelse($4,,[enable_$1_download=yes])
    fi

    # Also allow to disable the download of one tool if it is already locally installed

    TOOL_PARAMETERS($1,$2,$3)
])

echo "#  Build with FreeFEM with ./configure " `date`  >$ff_where_lib_conf
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



AC_CHECK_PROG(ff_unzip,unzip,yes,no)
AC_CHECK_PROG(ff_m4,m4,yes,no)
AC_CHECK_PROG(ff_bison,bison,yes,no)
AC_CHECK_PROG(ff_flex,flex,yes,no)
AC_CHECK_PROG(ff_cmake,cmake,yes,no)
AC_CHECK_PROG(ff_patch,patch,yes,no)
case "$ff_m4 $ff_bison $ff_patch $ff_flex $ff_unzip"  in
*no*)
    AC_MSG_NOTICE([ to install missing package  under debian or ubuntu, try ])
    test "$ff_unzip" = no && AC_MSG_NOTICE([ sudo apt-get install unzip])
    test "$ff_m4" = no && AC_MSG_NOTICE([ sudo apt-get install m4])
    test "$ff_bison" = no && AC_MSG_NOTICE([ sudo apt-get install bison])
    test "$ff_flex" = no && AC_MSG_NOTICE([ sudo apt-get install flex])
    test "$ff_patch" = no && AC_MSG_NOTICE([ sudo apt-get install patch])
    AC_MSG_ERROR([ Sorry missing  unzip,m4,bison,flex,patch command !]);;
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
 ff_ptrbit=32
elif test  "0$ff_size_ptr" -eq 8  ;then
 AC_SUBST(SIZEOF_PTRINBIT,64)
 ff_ptrbit=64
else
  AC_MSG_NOTICE( [ fatal error : sizeof pointer $ff_size_ptr !  or no c++ compiler: $CXX] );
  AC_MSG_ERROR([ Sorry  sizeof c++ pointer $ff_size_ptr are not  4 or 8 ])
fi

# FFCS - build the code for FreeFem++-cs
AC_ARG_ENABLE(ffcs,AC_HELP_STRING([--enable-ffcs],[build FreeFEM for use by FreeFem++-cs]))
if test "$enable_ffcs" = yes
then
	AC_DEFINE_UNQUOTED(ENABLE_FFCS,$enable_ffcs,[build FreeFEM for use by FreeFem++-cs])
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
	ff_f77=`pwd`/3rdparty/bin/fort77
	ff_flibs=`pwd`/3rdparty/lib/libf2c.a
	if test -x $ff_f77  -a -f $ff_flibs ; then
	    AC_MSG_WARN([ no fortran, but find download f2c/fort]);
	    F77=$ff_f77
	    FLIBS=$ff_flibs
	else
	    AC_MSG_NOTICE( [ fatal error : no fortran ] );
	    AC_MSG_NOTICE( [add --disable-fortran  ] );
	    AC_MSG_NOTICE( [or try to compile f2c in directory 3rdparty/f2c  ] );
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
ff_PURE_WIN32=0
case $ff_uname in
CYGWIN*|MINGW*|MSYS_NT*)
	ff_fpic=no
	ff_suffix_dylib="dll";
	ff_suffix_dylib_a="dll.a";
        ff_win32=yes;
    # echo "ff win 32 --------------------------------------------"
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
		    export PATH="$PATH:`cygpath $MSMPI_BIN`"
		    test -z "$MPIRUN" &&  MPIRUN=`which mpiexe.exe`
		    if test "$enable_fortran" != no -a "$with_flib" != no ;   then
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
		    ff_TEST_FFPP="../../src/bin-win32/FreeFem++.exe"
		    ff_TEST_FFPP_MPI="../../src/bin-win32/FreeFem++-mpi.exe"

            AC_DEFINE(PURE_WIN32,1,A pure windows applications no cygwin dll)
		ff_PURE_WIN32=1
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
AC_DEFINE_UNQUOTED(VersionFreeFem,$ff_numver,FreeFEM version as a float)

# Since src/fflib/strversionnumber.cpp is recreated at each build, this
# date is only useful for config-version.h

test "$prefix" = NONE && prefix="$ac_default_prefix"
ff_prefix_dir="${prefix}/lib/ff++/$VERSION"
AC_MSG_CHECKING(prefix dir freefem++  )
AC_MSG_RESULT($ff_prefix_dir)

FF_DATE=`date`
AC_DEFINE_UNQUOTED(VersionFreeFemDate,"$FF_DATE",FreeFEM build date)
AC_DEFINE_UNQUOTED(FF_PREFIX_DIR,"${ff_prefix_dir}",FreeFem prefix  dir)
AC_SUBST(ff_prefix_dir,$ff_prefix_dir)

ff_prefix_dir_lib="$ff_prefix_dir/lib"
ff_prefix_dir_data="$datadir/$PACKAGE_NAME"

ff_prefix_dir_include="$ff_prefix_dir/include"
ff_prefix_dir_etc="$ff_prefix_dir/etc"
ff_prefix_dir_lib_mpi="$ff_prefix_dir/lib/mpi"
ff_prefix_dir_example="$ff_prefix_dir_data/$PACKAGE_VERSION"
AC_SUBST(ff_prefix_dir_lib,$ff_prefix_dir_lib)
AC_SUBST(ff_prefix_dir_include,$ff_prefix_dir_include)
AC_SUBST(ff_prefix_dir_etc,$ff_prefix_dir_etc)
AC_SUBST(ff_prefix_dir_lib_mpi,$ff_prefix_dir_lib_mpi)
AC_SUBST(ff_prefix_dir_example,$ff_prefix_dir_example)

# Separate configuration header file for version information, included
# in config-macos9.h

cat << EOF > config-version.h
/* FreeFEM version information for MacOS 9 configuration */
/* File generated by configure.ac */

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "$PACKAGE_STRING"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "$PACKAGE_TARNAME"

/* Define to the version of this package. */
#define PACKAGE_VERSION "$PACKAGE_VERSION"

/* Version number of package */
#define VERSION "$VERSION"

/* FreeFEM build date */
#define VersionFreeFemDate "$FF_DATE"

/* FreeFEM version as a float */
#define VersionFreeFem ($ff_numver)
EOF

# Checking wether to produce a debug version
# ------------------------------------------

# Generic or hardware-dependant optimization
m4_include(etc/config/m4/acmacros.m4)
m4_include(etc/config/m4/acoptim.m4)

# [[file:ax_lib_hdf5.m4]]
m4_include(etc/config/m4/ax_lib_hdf5.m4)

# [[file:ax_lib_gsl.m4]]
m4_include(etc/config/m4/ax_lib_gsl.m4)

# Checking c++11 for HashMatix (Madatory now)
# F.Hecht:
# do not compile all with c++11 version to much
# probleme to day (frev 2015) with g++ 4.8.1
ff_cxx_save=$CXXFLAGS
CHECK_COMPILE_FLAG(C++,[-std=c++11],CXXFLAGS,ff_cxx11)
dnl CXXFLAGS="$CXX11FLAGS $CXXFLAGS"
dnlff_cxx11="no"
dnl if test -n "$CXX11FLAGS" ; then ff_cxx11="yes" ; fi
if test "$ff_cxx11" = yes ; then
    CXX11FLAGS="-std=c++11";
##	CXXFLAGS="$CXXFLAGS -std=c++11"
    AC_FF_ADDWHERELIB(cxx11,$CXX11FLAGS,$CXX11FLAGS)
else
  AC_MSG_NOTICE( [ fatal error : c++ compiler ! ] );
  AC_MSG_ERROR([ Sorry no c++ compiler !])
  AC_MSG_ERROR([ Fatal error c++11  ],1);
fi
##CXXFLAGS="$ff_cxx_save"
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
    m4_include(etc/config/m4/acmpi.m4)
else
    # FFCS - use the same MPI configuration choices as FFCS
    if test -z "$MPICXX"
    then
	ff_mpi=no
    else
	ff_mpi=yes
	AC_SUBST(MPICXX,$MPICXX)
	AC_SUBST(MPICC,$MPICC)
	AC_SUBST(MPIFC,$MPIFC)
	AC_SUBST(MPIPROG,"FreeFem++-mpi${EXEEXT}")
	AC_SUBST(MPI_INCLUDE,"-I $MPI_INC_DIR")
	AC_SUBST(MPI_INC_DIR,$MPI_INC_DIR)
	AC_SUBST(MPI_LIB_DIRS,"")
	AC_SUBST(MPI_LIB,$MPI_LIB)
	AC_SUBST(MPI_LIBC,"")
	AC_SUBST(MPI_LIBFC,"")

        # these values should not be empty otherwise
        # plugin/seq/ff-get-dep will think that they are not
        # defined
	AC_FF_ADDWHERELIB(mpi,-DDUMMY,-I$MPI_INC_DIR)

        # mpifc and mpif77 libraries should always be specified
        # because FF never calls the Fortran MPI compiler. It always
        # uses mpicxx in [[file:plugin/seq/ff-c++]]. The
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

# Looking for useful configuration utilities
# ------------------------------------------

AC_ARG_ENABLE(summary,AC_HELP_STRING([--enable-summary],[Display activated libraries list at the end of the configure process]))

AC_ARG_ENABLE(download,AC_HELP_STRING([--enable-download],[Download missing libraries (BLAS,ARPACK,UMFPACK,...)]))
## PETSc
AC_ARG_WITH(petsc,[ --with-petsc=/usr/local/petsc/conf/petscvariables --without-petsc])
AC_ARG_WITH(prefix_petsc,[ --prefix_petsc=directory where all PETSc/SLEPc will be installed if you download ])
if test -z "${prefix_petsc}"; then prefix_petsc="${prefix}/ff-petsc"; fi
AC_SUBST([FF_prefix_petsc],"$prefix_petsc")
AC_SUBST([FF_prefix_petsc_real],"$prefix_petsc_real")
AC_SUBST([FF_prefix_petsc_complex],"$prefix_petsc_complex")
AC_SUBST([FF_generic_petsc],"$enable_generic")
# /usr/local/petsc/conf/petscvariables
#/usr/local/ff++/petsc/lib/petsc/conf/petscvariables
ff_petsc_ok=no
ff_slepc_ok=no
if test "$with_petsc" != no; then
    for d in "$with_petsc" "${prefix_petsc}/r/lib" /usr/local/ff++/petsc/lib  /usr /usr/local /opt/usr  /opt/local ;do
	for dd in "" "/petsc/conf/petscvariables" "/lib/petsc/conf/petscvariables" ; do
	    ffconfpetsc="$d$dd"
		ffconfslepc="$d/slepc/conf/slepcvariables"

		 echo " petsc ... $ffconfpetsc"
	    if test -f "$ffconfpetsc" ; then
		if test  "$ff_petsc_ok" = no ; then
		    PETSC_PREFIXDIR=`awk -F' *= *' '"PREFIXDIR"==$1 {print $2}' "$ffconfpetsc"`
		    PETSC_SCALAR=`awk -F' *= *' '"PETSC_SCALAR"==$1 {print $2}' "$ffconfpetsc"`
		    PETSC_MPIRUN=`awk -F' *= *' '"MPIEXEC"==$1 {print $2}' "$ffconfpetsc"`
		    PETSC_MPIRUN=`echo "$PETSC_MPIRUN" | awk '{gsub(" --oversubscribe", ""); print}'`	#remove --oversubscribe
		    AC_MSG_NOTICE([ find real ( $PETSC_SCALAR )  petsc in $ffconfpetsc]);
            if test "$ff_mpi" != yes; then
                MPIRUN=`echo $PETSC_MPIRUN`
                MPICC=`awk -F' *= *' '"CC"==$1 {print $2}' "$ffconfpetsc"`
                MPICXX=`awk -F' *= *' '"CXX"==$1 {print $2}' "$ffconfpetsc"`
                MPIFC=`awk -F' *= *' '"FC"==$1 {print $2}' "$ffconfpetsc"`
                m4_include(etc/config/m4/acmpi.m4)
            fi
            if test "$ff_mpi" = yes; then
                    WPETSC_MPIRUN=`which "$PETSC_MPIRUN"`
                    WMPIRUN=`which "$MPIRUN"` 
		    test  "$PETSC_SCALAR" = "real" && test "$ff_win32" = yes -o  "$PETSC_MPIRUN" = "$MPIRUN" ||  cmp -s "$WPETSC_MPIRUN" "$WMPIRUN"
		    if test "$?"  -eq 0 ; then
	                AC_SUBST([FF_conf_petsc_real],"$ffconfpetsc")
		        AC_SUBST([FF_prefix_petsc_real],"$PETSC_PREFIXDIR")
		        AC_FF_PETSC_WHERELIB_BASIC(petsc,$ffconfpetsc)
		        AC_FF_PETSC_WHERELIB(parmetis,PARMETIS_LIB,PARMETIS_INCLUDE,$ffconfpetsc)
		        AC_FF_PETSC_WHERELIB(metis,METIS_LIB,METIS_INCLUDE,$ffconfpetsc)
		        AC_FF_PETSC_WHERELIB(ptscotch,PTSCOTCH_LIB,PTSCOTCH_INCLUDE,$ffconfpetsc)
                        ff_SCOTCH_INCLUDE=$ff_PTSCOTCH_INCLUDE
                        ff_SCOTCH_LIB=`echo $ff_PTSCOTCH_LIB | sed -e 's/@<:@^ @:>@*ptscotch@<:@^ @:>@*//g' | sed -e 's/@<:@^ @:>@*ptes@<:@^ @:>@*//g'`
                        if test -n "$ff_SCOTCH_LIB"; then
                          ff_scotch_ok=yes;
                          AC_FF_ADDWHERELIB(scotch,$ff_SCOTCH_LIB,$ff_SCOTCH_INCLUDE)
                          AC_SUBST([TOOL_COMPILE_scotch],"")
                        fi
		        AC_FF_PETSC_WHERELIB(mumps,MUMPS_LIB,MUMPS_INCLUDE,$ffconfpetsc)
                        AC_FF_PETSC_WHERELIB(tetgen,TETGEN_LIB,TETGEN_INCLUDE,$ffconfpetsc)
                        AC_FF_PETSC_WHERELIB(hpddm,HPDDM_LIB,HPDDM_INCLUDE,$ffconfpetsc)
                        AC_FF_PETSC_WHERELIB(superlu,SUPERLU_LIB,SUPERLU_INCLUDE,$ffconfpetsc)
		        AC_FF_PETSC_WHERELIB(scalapack,SCALAPACK_LIB,SCALAPACK_INCLUDE,$ffconfpetsc)
  		        AC_FF_PETSC_WHERELIB(suitesparse,SUITESPARSE_LIB,SUITESPARSE_INCLUDE,$ffconfpetsc,$ff_METIS_LIB)
 dnl               if test -n "$ff_SUITESPARSE_LIB" -a -n "$ff_METIS_LIB"; then
 dnl                   ff_where_lib_suitesparse="$ff_METIS_LIB $ff_SUITESPARSE_LIB"
 dnl               fi
			AC_FF_PETSC_WHERELIB(blaslapack,BLASLAPACK_LIB,BLASLAPACK_INCLUDE,$ffconfpetsc)
                        if test -f $d/libslepc.$ff_suffix_dylib -o -f $d/libslepc.a ; then
			   ff_slepc_ok=yes
                           ff_SLEPC_LIB="-Wl,-rpath,$d -L$d -lslepc"
			   AC_FF_ADDWHERELIB(slepc,$ff_SLEPC_LIB,$ff_where_inc_petsc)
		        fi
            if test "$ff_slepc_ok" != yes ;   then
                AC_ARG_WITH(slepc_include,AC_HELP_STRING([--with-slepc-include=],[Include directives for slepc instead of automatic download]))
                AC_ARG_WITH(slepc_ldflags,AC_HELP_STRING([--with-slepc-ldflags=],[Link-time directives for slepc instead of automatic download]))
                if test "$with_slepc_include" != "" || test "$with_slepc_ldflags" != ""
                then
                    # some directives have been specified, use them instead of downloading
                    AC_FF_ADDWHERELIB(slepc,$with_slepc_ldflags,$with_slepc_include)
                    ff_slepc_ok=yes
                fi
            fi
		    else
			 AC_MSG_NOTICE([" Warning PETSC MPI and FF++ MPI not the same: $PETSC_MPIRUN != $MPIRUN or $PETSC_SCALAR != real ." ]);
		    fi
        fi
		fi
	    else
		ffconfpetsc="";
	    fi
	done
    done
fi

if test "$ff_petsc_ok" != no ;   then
    AC_MSG_NOTICE([ with petsc  $ff_petsc_ok,  with slepc $ff_slepc_ok ])
else
    AC_MSG_NOTICE([ without petsc, slepc  ***** ])
fi

## PETSc complex
AC_ARG_WITH(petsc_complex,[  --with-petsc_complex=/usr/local/petsc/conf/petscvariables   --without-petsc_complex  ])
ff_petsccomplex_ok=no
ff_slepccomplex_ok=no
if test "$with_petsc_complex" != no; then
    for d in "$with_petsc_complex" "${prefix_petsc}/c/lib" /usr/local/ff++/petsc/c/lib ;do
	for dd in "" "/petsc/conf/petscvariables" "/lib/petsc/conf/petscvariables" ; do
	    ffconfpetscc="$d$dd"
	    if test -f "$ffconfpetscc" -a  "$ff_petsccomplex_ok" = no ; then
		    PETSC_COMPLEX_SCALAR=`awk -F' *= *' '"PETSC_SCALAR"==$1 {print $2}' "$ffconfpetscc"`
		    PETSC_COMPLEX_MPIRUN=`awk -F' *= *' '"MPIEXEC"==$1 {print $2}' "$ffconfpetscc"`
		    PETSC_COMPLEX_PREFIXDIR=`awk -F' *= *' '"PREFIXDIR"==$1 {print $2}' "$ffconfpetscc"`
		    
			PETSC_COMPLEX_MPIRUN=`echo "$PETSC_COMPLEX_MPIRUN" | awk '{print $1}'`	#remove --oversubscribe
			AC_MSG_NOTICE([ find complex ($PETSC_COMPLEX_SCALAR) petsc in $ffconfpetscc]);
            if test "$ff_mpi" != yes; then
                MPIRUN=`echo $PETSC_COMPLEX_MPIRUN`
                MPICC=`awk -F' *= *' '"CC"==$1 {print $2}' "$ffconfpetscc"`
                MPICXX=`awk -F' *= *' '"CXX"==$1 {print $2}' "$ffconfpetscc"`
                MPIFC=`awk -F' *= *' '"FC"==$1 {print $2}' "$ffconfpetscc"`
                m4_include(etc/config/m4/acmpi.m4)
            fi
            if test "$ff_mpi" = yes; then
		    test  "$PETSC_COMPLEX_SCALAR" = "complex" &&  test  "$ff_win32" = yes -o  "$PETSC_COMPLEX_MPIRUN" = "$MPIRUN" ||  cmp -s `which "$PETSC_COMPLEX_MPIRUN"` `which "$MPIRUN"`
		   	if test "$?"  -eq 0 ; then
            AC_SUBST([FF_prefix_petsc_complex],"$PETSC_COMPLEX_PREFIXDIR")
            AC_SUBST([FF_conf_petsc_complex],"$ffconfpetscc")
	    
            AC_FF_PETSC_WHERELIB_BASIC(petsccomplex,$ffconfpetscc)
            AC_FF_PETSC_WHERELIB(parmetis,PARMETIS_LIB,PARMETIS_INCLUDE,$ffconfpetscc)
            AC_FF_PETSC_WHERELIB(metis,METIS_LIB,METIS_INCLUDE,$ffconfpetscc)
            AC_FF_PETSC_WHERELIB(ptscotch,PTSCOTCH_LIB,PTSCOTCH_INCLUDE,$ffconfpetscc)
            ff_SCOTCH_INCLUDE=$ff_PTSCOTCH_INCLUDE
            ff_SCOTCH_LIB=`echo $ff_PTSCOTCH_LIB | sed -e 's/@<:@^ @:>@*ptscotch@<:@^ @:>@*//g' | sed -e 's/@<:@^ @:>@*ptes@<:@^ @:>@*//g'`
            if test -n "$ff_SCOTCH_LIB"; then
                ff_scotch_ok=yes;
                AC_FF_ADDWHERELIB(scotch,$ff_SCOTCH_LIB,$ff_SCOTCH_INCLUDE)
                AC_SUBST([TOOL_COMPILE_scotch],"")
            fi
            AC_FF_PETSC_WHERELIB(mumps,MUMPS_LIB,MUMPS_INCLUDE,$ffconfpetscc)
            AC_FF_PETSC_WHERELIB(tetgen,TETGEN_LIB,TETGEN_INCLUDE,$ffconfpetscc)
            AC_FF_PETSC_WHERELIB(hpddm,HPDDM_LIB,HPDDM_INCLUDE,$ffconfpetscc)
            AC_FF_PETSC_WHERELIB(superlu,SUPERLU_LIB,SUPERLU_INCLUDE,$ffconfpetscc)
            AC_FF_PETSC_WHERELIB(scalapack,SCALAPACK_LIB,SCALAPACK_INCLUDE,$ffconfpetscc)
            AC_FF_PETSC_WHERELIB(suitesparse,SUITESPARSE_LIB,SUITESPARSE_INCLUDE,$ffconfpetscc)
            AC_FF_PETSC_WHERELIB(blaslapack,BLASLAPACK_LIB,BLASLAPACK_INCLUDE,$ffconfpetscc)
            if test -f $d/libslepc.$ff_suffix_dylib -o -f $d/libslepc.a ; then
				   ff_slepccomplex_ok=yes
                   ff_SLEPCCOMPLEX_LIB="-Wl,-rpath,$d -L$d -lslepc"
			       AC_FF_ADDWHERELIB(slepccomplex,$ff_SLEPCCOMPLEX_LIB,$ff_where_inc_petsccomplex)
			fi
            if test "$ff_slepccomplex_ok" != yes ;   then
                AC_ARG_WITH(slepccomplex_include,AC_HELP_STRING([--with-slepccomplex-include=],[Include directives for slepccomplex instead of automatic download]))
                AC_ARG_WITH(slepccomplex_ldflags,AC_HELP_STRING([--with-slepccomplex-ldflags=],[Link-time directives for slepccomplex instead of automatic download]))
                if test "$with_slepccomplex_include" != "" || test "$with_slepccomplex_ldflags" != ""
                then
                    # some directives have been specified, use them instead of downloading
                    AC_FF_ADDWHERELIB(slepccomplex,$with_slepccomplex_ldflags,$with_slepccomplex_include)
                    ff_slepccomplex_ok=yes
                fi
            fi
			else
			 AC_MSG_NOTICE([" Warning PETSC complex MPI and FF++ MPI not the same: $PETSC_COMPLEX_MPIRUN != $MPIRUN or $PETSC_COMPLEX_SCALAR != complex ." ]);
		    fi
            fi
		fi
	done
    done
fi


if test "$ff_petsccomplex_ok" != no ;   then
    AC_MSG_NOTICE([ with petsc complex  $ff_petsccomplex_ok, with slepc complex $ff_slepccomplex_ok ])
else
    AC_MSG_NOTICE([ without petsc complex ***** ])
fi
##  slepc

### end petsc ....

# FFCS needs to change some of the FF makefiles to compile without MPI
AM_CONDITIONAL([FFCS_MPIOK],[test $ff_mpi = yes])

#replacing wget with another command
ff_wget_command="wget --no-check-certificate"
ff_wget=no

AC_ARG_WITH(wget,
	[  --with-wget=command	Replace "wget" with another command. Implies --enable-download],
	   enable_download=yes
	   ff_wget=yes
	   ff_wget_command="${withval}")

ff_MKL_libpath=yes
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
enable_download=no
fi
if test "$enable_download"  = "yes"  ; then
 AC_CHECK_PROG(ff_git,git,yes,no)
 case "$ff_git"  in
 *no*)
     AC_MSG_NOTICE([ to install missing package  under debian or ubuntu, try ])
     test "$ff_git" = no && AC_MSG_NOTICE([ sudo apt-get install git])
     AC_MSG_ERROR([ Sorry missing  git command !]);;
 esac
fi

AC_SUBST(WGET,$ff_wget_command)

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
    ff_cadna_dir="${curdir}/3rdparty/cadna"
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

# used in [[file:3rdparty/fftw/Makefile.am::DOWNLOAD_FFTW]]

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
            ff_mkl_ok=yes
			ff_mkl=yes
dnl            AC_DEFINE(HAVE_MKL,1, the MKL intel lib is present for BLAS and LAPACK  )
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
		mkl_l="-Wl,-rpath,$ff_mkl_libpath -L$ff_mkl_libpath"
		mkl_piomp5=$(realpath $ff_mkl_libpath/../../compiler/lib/)
		ff_mkl_inc=$(realpath $ff_mkl_libpath/../include)
		ff_blas_inc="-I$ff_mkl_inc"
		
		if test -d "$mkl_piomp5" ; then mkl_l="$mkl_l -Wl,-rpath,$mkl_piomp5 -L$mkl_piomp5" ; fi
#		echo " ################### $mkl_l $ff_blas_inc " 
	    if test  -f "$ff_mkl_libpath/libmkl_rt.$ff_suffix_dylib" ; then
		 mkl_blas=rt
		 mkl_lapack=rt
		 ff_mkl_mlt="$mkl_l -lmkl_rt -lmkl_intel_thread -lmkl_core  $ff_iomp5 $ff_pthread"
		 ff_mkl_blas="$mkl_l -lmkl_rt -lmkl_sequential -lmkl_core   $ff_iomp5 $ff_pthread"
		 ff_mkl_lapack="$ff_mkl_blas"
	    elif test -f "$ff_mkl_libpath/libmkl_lapack.$ff_suffix_dylib" ; then
		 ff_mkl_mlt="$mkl_l -lmkl_${ff_mkl_cc}${ff_mkl_lp} -lmkl_lapack -lmkl_intel${ff_mkl_lp} -lmkl_intel_thread  -lmkl_core -lguide  -lm -lpthread"
		 ff_mkl_blas="$mkl_l -lmkl_${ff_mkl_cc}${ff_mkl_lp} -lmkl_lapack -lmkl_intel${ff_mkl_lp} -lmkl_sequential -lmkl_core -lguide -lm -lpthread"
		 ff_mkl_lapack="$mkl_l -lmkl_lapack"
	    else
		 ff_mkl_mlt="$ff_ldbg $mkl_l  -lmkl_${ff_mkl_cc}${ff_mkl_lp}  -lmkl_intel${ff_mkl_lp}  -lmkl_intel_thread  -lmkl_core -lmkl_intel_thread $ff_iomp5  $ff_ldeg  -lm $ff_pthread"
		 ff_mkl_blas="$ff_ldbg $mkl_l -lmkl_${ff_mkl_cc}${ff_mkl_lp}  -lmkl_intel${ff_mkl_lp} -lmkl_sequential -lmkl_core  -lmkl_sequential $ff_ldeg   -lm $ff_pthread"
		 ff_mkl_lapack="$ff_mkl_blas"
	    fi
	    if test -f "-L$ff_mkl_libpathmkl_scalapack${ff_mkl_lp}" ; then
            ff_mkl_scalapack="$mkl_l -lmkl_blacs${ff_mkl_mpi}${ff_mkl_lp} -lmkl_scalapack${ff_mkl_lp} "
            ff_mkl_blacs="$mkl_l -lmkl_blacs${ff_mkl_mpi}${ff_mkl_lp} -lmkl_scalapack${ff_mkl_lp}"
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
		ff_mkl_ok=yes
		ff_mkl=yes
dnl            AC_DEFINE(HAVE_MKL,1, the MKL intel lib is present for BLAS and LAPACK  )

	else
	    ff_mkl_libpath=no
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
    if test -z "$ff_dll_blas" ; then  ff_dll_blas="no"; fi
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

	# it is necessary to split the library path into -L and -l otherwise ff/upstream/plugin/seq/ff-c++ will place -lgfortran in
	# LIBS, .../libgoto___.lib in OTHER, and call the compiler with $LIBS before $OTHER, which will fail. But libopenblas.a is placed in
	# the standard directory 3rdparty/lib so the -L option is not required anymore.

	ff_blas_libs="-L${curdir}/3rdparty/lib -lopenblas $LIBS $ff_pthread"
	ff_blas_inc="-I${curdir}/3rdparty/include"

        # skip generic blas compilation and activate OpenBLAS (see [[file:3rdparty/blas/Makefile.am::OpenBLAS]])
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
    ff_blas_libs="-L${curdir}/3rdparty/blas -lcblas -lf77blas"
    ff_blas_inc="-I${curdir}/3rdparty/blas/CBLAS/src"

    AC_SUBST(DOWNLOADED_BLAS,"libf77blas.a libcblas.a")
    AC_SUBST(DOWNLOADED_BLAS_BUILT_SOURCES,'BLAS CBLAS $(F77BLAS_SRC) $(CBLAS_SRC)')

    AC_MSG_RESULT($ff_download_blas)
    ff_blas_ok=yes
    ff_cblas_h=yes
    ff_cblas_libs=yes
fi

AC_SUBST(BLASLIBS,$ff_blas_libs)
AC_SUBST(BLASINC,$ff_blas_inc)

#looking for cblas.h  FH

if test "$ff_blas_ok" = yes -a -z "$ff_download_blas" -a "$ff_mkl_ok" != yes; then
  AC_CHECK_HEADERS(cblas.h,
    ff_cblas_h=yes,
    ff_cblas_h=no)
  AC_CHECK_HEADERS(Accelerate/cblas.h,
    ff_cblas_h=yes,
    ff_cblas_h=no)
  AC_CHECK_HEADERS(atlas/cblas.h,
    ff_cblas_h=yes,
    ff_cblas_h=no)
  AC_CHECK_LIB([cblas], [cblas_zdotu_sub], [LIBS+=-lcblas])
fi

# ALH - 18/9/13 - [[file:3rdparty/blas/Makefile.am]] does not update the WHERE mechanism so it always needs to be set up
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

if test "$ff_lapack_ok" = no; then
    AC_ARG_WITH(lapack,
	AC_HELP_STRING([--with-lapack=library],[Use a specific version of Lapack]),
	ff_lapack_ok=yes
	ff_lapack_lib="${withval}"
	LIBS="$ff_lapack_lib $LIBS")
fi
AC_MSG_RESULT($ff_lapack_ok)
if test "$ff_lapack_ok" != no ; then
AC_SUBST(LAPACKLIBS,$ff_lapack_lib)
fi

# Lapack configuration for Arpack

ff_lapackdir='$$(LAPACKdir)'
if test "$ff_lapack_ok" = yes
then
    # no compilation of lapack in arpack
    ff_lapackdir=
else
    ff_lapack_lib="-L${curdir}/3rdparty/lib -llapack"
fi

# Arpack itself
AC_ARG_ENABLE(download-arpack,AC_HELP_STRING([--enable-download_arpack],[force the use  download of arpack]))
#echo "****** enable_download_arpack=$enable_download_arpack"

ff_arpack_ok=no
ff_save_libs="$LIBS"
ff_arpack_libs=
if test "$ff_blas_ok" = yes ;
then
    if test "$enable_download_arpack" != yes ; then
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
	fi
    fi
    # If all else fails, download!
    if test "$enable_fortran" != no -o -z "$enable_fortran"
    then
    if test "$ff_arpack_ok" = no -a "$enable_download" = yes -o "$enable_download_arpack" = yes
    then
	ff_arpack_download=yes
	AC_MSG_NOTICE(using downloaded Arpack)

		# ALH - 6/11/13 - this install goal is the standard goal for all downloaded packages in
		# [[file:3rdparty/common.mak::install]]

	AC_SUBST(DOWNLOAD_ARPACK,install)

	AC_SUBST(FF_LAPACKdir,$ff_lapackdir)
	AC_SUBST(ARPACKLIB,${curdir}/3rdparty/lib/libarpack.a)
	AC_SUBST(LAPACK_arpack_LIB,${curdir}/3rdparty/lib/liblapack.a)

		# Do not update $LIBS, but create an extra LIB variable, because this lib does not exist yet, and this
		# could make the following tests fail.

		# ALH - 30/9/13 - do not use the "-L ${curdir}/3rdparty/lib" directive because it would allow other
		# following -l directives (eg -lumfpack) to pick an old locally compiled library instead of the system
		# ones.

	ff_arpack_libs="${curdir}/3rdparty/lib/libarpack.a $ff_lapack_lib"
	ff_arpack_ok=yes
    fi
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
AM_CONDITIONAL([ENABLE_EIGEN], [test "$ff_arpack_ok" = "yes"])

AC_SUBST([EIGENOBJ])

# ALH - 25/9/13 - <<WHERE_lapack>> always run the WHERE LD statement for lapack because some libraries in
# plugin/seq require it (eg [[file:plugin/seq/Element_Mixte.cpp::lapack]]). Request from Fred.

AC_FF_ADDWHERELIB(lapack,$ff_lapack_lib,)

# Looking for UMFPACK
# -------------------
UMFPACK_CPPFLAGS="";
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
    UMFPACK_CPPFLAGS="$ff_where_inc_suitesparse"
    dnl   AC_SUBST([UMFPACK_CPPFLAGS],[$ff_where_inc_suitesparse])
    AC_DEFINE(HAVE_UMFPACK_H,1,If umfpack.h is located in UMFPACK subdir)
fi
#echo "@@@@@@@"
ff_save_libs="$LIBS"
AC_ARG_ENABLE(system_umfpack,AC_HELP_STRING([--disable-system-umfpack],[Disable the automatic detection of umfpack, colmod, amd, ....]))
if test "$ff_umfpack_ok"  = no -a "$enable_system_umfpack" != no -a  "$ff_blas_ok" = yes;
then
#    echo "@@@@@@@qqqq"
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
	[ff_umfpack_header=yes
	    ff_umfpack_dir=`dirname $ac_header`
	    break])

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
	ff_umfpack_incs=""
	UMFPACK_CPPFLAGS=""
    fi

    if test "$ff_umfpack_ok" = yes -a "$ff_amd_ok" = yes; then
	 if test -n "$ff_umfpack_incs" ; then
	  UMFPACK_CPPFLAGS="$ff_umfpack_incs"
	 fi
dnl    AC_SUBST([UMFPACK_CPPFLAGS],[$ff_umfpack_incs])
	  AC_DEFINE(HAVE_LIBUMFPACK,1,Umfpack is used for sparse matrices computations)
    fi
fi

LIBS="$ff_save_libs"
#echo "@@@@22@@@"
# If all else fails, download!
#  attention no /usr/include in WHERE
if test "$ff_umfpack_ok" = yes -a -n "$ff_umfpack_dir" -a -z "$ff_umfpack_incs"
then
	ff_umfpack_incs="-I/usr/include/$ff_umfpack_dir"
	UMFPACK_CPPFLAGS="$ff_umfpack_incs"
fi
if test "$ff_umfpack_ok" = no -a "$enable_download" = yes
then
    AC_MSG_NOTICE(using downloaded UMFPACK)
    AC_SUBST(DOWNLOAD_UMFPACK,"umfpack")
    ff_umfpack_download=yes

    # Do not update $LIBS, but create an extra LIB variable, because this lib does not exist yet, and this could make
    # the following tests fail.

    ff_umfpack_libs="-L${curdir}/3rdparty/lib -lumfpack -lcholmod -lcolamd -lamd -lsuitesparseconfig"
    ff_umfpack_incs="-I${curdir}/3rdparty/include"
	UMFPACK_CPPFLAGS="$ff_umfpack_incs"
    AC_DEFINE(HAVE_LIBUMFPACK,1,UMFPACK)
    if test "$ff_win32" = yes; then
  	AC_SUBST(FF_UMFPACK_CONFIG,-DCBLAS)
    fi
    ff_umfpack_ok=yes
fi
AC_SUBST([UMFPACK_CPPFLAGS],[$UMFPACK_CPPFLAGS])

# ALH - 17/9/13 - moved UMFPACK configuration settings in wherelib to _after_ configuring the download version because
# [[file:3rdparty/umfpack/Makefile.am]] does not set the WHERE mechanism. Also removed -I/usr/include/$ff_umfpack_dir
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
#echo "@@@@33@@@@"
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
AC_CHECK_HEADERS(semaphore.h,ff_sem=1,ff_sem=0)
AC_CHECK_HEADERS(sys/mman.h,ff_mmap=1,ff_mmap=0)


if test "$ff_sem" -eq  1 -a  "$ff_mmap" -eq  1 ; then
BIN_ffmaster=ffmaster
else
enable_MMAP=no
BIN_ffmaster=
fi

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
if test "$ff_PURE_WIN32" -ne 1 ; then
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
if test "$ff_dynload" = yes
then
	AC_DEFINE(HAVE_DLFCN_H,1,Dynamic loading - not mandatory)
fi
fi  # endif ff_PURE_WIN32 != 1
#  the -rdynamic don't exist on macos and sunOS
if test "$ff_dynload" = yes
then
	# Activate dynamic loading tests (see plugin/seq/Makefile.am)
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
#  the doc is now in https://github.com/FreeFem/FreeFem-doc-pdf/blob/master/freefem%2B%2Bdoc.pdf
# Checking wether we can generate some documentation
# --------------------------------------------------
dnl
dnl AC_CHECK_PROG(ff_latex,latex,yes,no)
dnl AC_CHECK_PROG(ff_makeindex,makeindex,yes,no)
dnl AC_CHECK_PROG(ff_dvips,dvips,yes,no)
dnl # to translate the figure
dnl AC_CHECK_PROG(ff_pdf2ps,pdf2ps,yes,no)
dnl AC_CHECK_PROGS(EPSTOPDF,[epstopdf pstopdf],[false])
dnl AC_CHECK_PROG(ff_convert,convert,yes,no)
dnl if test "$ff_latex" = yes -a "$ff_makeindex" = yes -a "$ff_dvips" = yes -a $ff_pdf2ps = yes -a $ff_convert = yes;
dnl then
dnl 	AC_SUBST(DOCPS,"freefem++doc.ps")
dnl
dnl 	AC_CHECK_PROG(ff_gzip,gzip,yes,no)
dnl 	if test "$ff_gzip" = yes;
dnl 	then
dnl 		AC_SUBST(DOCPSGZ,"freefem++doc.ps.gz")
dnl 	fi
dnl fi

# PDF documentation building sometimes poses problems because of pdfsync.sty. So we need to be able to disable it.
# FFCS - 9/4/14 - need to disable PDF also during FFCS build (which does not use the resulting PDF).

dnl AC_ARG_ENABLE(pdf,[  --disable-pdf	Disable PDF documentation building])
dnl if test "$enable_pdf" != no
dnl then
dnl    AC_CHECK_PROG(ff_pdflatex,pdflatex,yes,no)
dnl    if test "$ff_pdflatex" = yes -a $EPSTOPDF != false -a $ff_convert = yes;
dnl    then
dnl 	AC_SUBST(DOCPDF,"freefem++doc.pdf")
dnl    fi
dnl fi

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

m4_map([AC_FF_WHERELIB],[
	[[mumps],[-ldmumps -lzmumps  -lmumps_common  -lpord],[dmumps_c.h],[]],
	[[mumps_seq],[-ldmumps_seq -lzmumps_seq  -lmumps_common_seq  -lpord_seq -lmpiseq_seq],[dmumps_c.h],[]],
	[[libseq],[-lmpiseq_seq],[${inc_usr_include}mumps_seq/mpi.h]],
	[[mumps_ptscotch],[-lpord_ptscotch -lmumps_common_ptscotch -ldmumps_ptscotch -lzmumps_ptscotch -lpord_ptscotch],[dmumps_c.h]],
	[[mumps_scotch],[-lpord_scotch -lmumps_common_scotch -ldmumps_scotch -lzmumps_scotch -lpord_scotch],[dmumps_c.h]],
	[[fftw3],[-lfftw3],[${inc_usr_include}fftw3.h],[]],
	[[superlu],[-lsuperlu],[${inc_usr_include}superlu/superlu_enum_consts.h],[]],
        [[blacs],[ -lblacsCinit$ff_with_mpi	-lblacsF77init$ff_with_mpi	-lblacs$ff_with_mpi],[]],
        [[scalapack],[-lscalapack$ff_with_mpi],[]],
        [[scotch],[-lscotch -lscotcherr],[scotch.h]],
        [[ptscotch],[-lptscotch -lptscotcherr],[ptscotch.h]],
	[[metis],[-lmetis],[${inc_usr_include}metis/metis.h],[]],
	[[metis],[-lmetis],[metis.h],[]],
	[[parmetis],[-lparmetis -lmetis],[],[]],
	[[freeyams],[-lfreeyams],[freeyamslib.h],[]],
	[[mmg3d],[-lmmg3d-v4],[libmmg3d.h],[]],
        [[mmg],[-lmmg],[libmmg.h],[]],
        [[parmmg],[-lparmmg],[libparmmg.h],[]],
	[[mshmet],[-lmshmet],[mshmetlib.h],[]],
dnl	[[gsl],[-lgsl -lgslcblas -lm],[gsl/gsl_sf.h],[]],
	[[tetgen],[-ltet],[tetgen.h],[]],
	[[ipopt],[ -lipopt],[${inc_usr_include}coin/IpTNLP.hpp],[]],
	[[nlopt],[ -lnlopt],[nlopt.hpp],[]]
    ]
    )
## before try

TOOL_DISABLE(tetgen,tetgen.$DYLIB_SUFFIX)
TOOL_DISABLE(lapack,lapack.$DYLIB_SUFFIX)
TOOL_DISABLE(metis,metis.$DYLIB_SUFFIX)
TOOL_DISABLE(parmetis,parmetis.$DYLIB_SUFFIX)
TOOL_DISABLE(mmg3d,mmg3d-v4.0.$DYLIB_SUFFIX)
if test "$ff_cmake" = yes; then
     TOOL_DISABLE(mmg,mmg.$DYLIB_SUFFIX)
     TOOL_DISABLE(parmmg,parmmg.$DYLIB_SUFFIX)
fi
TOOL_DISABLE(mshmet,"mshmet.$DYLIB_SUFFIX aniso.$DYLIB_SUFFIX")
TOOL_DISABLE(gmm,ilut.$DYLIB_SUFFIX,,enable_gmm=no)
if test "$enable_fortran" != no -o -z "$enable_fortran"; then
TOOL_DISABLE(ipopt,ff-Ipopt.$DYLIB_SUFFIX)
TOOL_DISABLE(scalapack,"MUMPS.$DYLIB_SUFFIX MUMPS_mpi.$DYLIB_SUFFIX")
TOOL_DISABLE(mumps,"MUMPS.$DYLIB_SUFFIX MUMPS_mpi.$DYLIB_SUFFIX")
TOOL_DISABLE(mumps_seq,"MUMPS_seq.$DYLIB_SUFFIX MUMPS.$DYLIB_SUFFIX",mumps-seq)
AC_SUBST([TOOL_DYLIB_ffnewuoa],ffnewuoa.$DYLIB_SUFFIX)
else
TOOL_DISABLE_NO(ipopt,ff-Ipopt.$DYLIB_SUFFIX)
TOOL_DISABLE_NO(scalapack,"MUMPS.$DYLIB_SUFFIX MUMPS_mpi.$DYLIB_SUFFIX")
TOOL_DISABLE_NO(mumps,"MUMPS.$DYLIB_SUFFIX MUMPS_mpi.$DYLIB_SUFFIX")
TOOL_DISABLE_NO(mumps_seq,"MUMPS_seq.$DYLIB_SUFFIX MUMPS.$DYLIB_SUFFIX",mumps-seq)	
AC_SUBST([TOOL_DYLIB_ffnewuoa],"")
ff_uncompile_plugin+=" ffnewuoa"
fi
TOOL_DISABLE(nlopt,ff-NLopt.$DYLIB_SUFFIX)

TOOL_DISABLE(scotch,scotch.$DYLIB_SUFFIX)
TOOL_DISABLE(superlu,SuperLu.$DYLIB_SUFFIX)
TOOL_DISABLE(umfpack,UMFPACK64.$DYLIB_SUFFIX)
TOOL_DISABLE(yams,freeyams.$DYLIB_SUFFIX)
TOOL_DISABLE(pipe,pipe.$DYLIB_SUFFIX,,nodownload)
AC_SUBST([BIN_ffmaster],"$BIN_ffmaster")

if test "$ff_win32" = "yes" ; then
TOOL_DISABLE(libpthread_google,ff-mmap-semaphore.$DYLIB_SUFFIX,libpthread-google)	
else
TOOL_DISABLE(MMAP,ff-mmap-semaphore.$DYLIB_SUFFIX,,nodownload)
fi

ff_pthread_google=$ff_pthread
if test "$enable_libpthread_google_download" = "yes" ; then
    ff_pthread_google="-I${curdir}/3rdparty/include -L${curdir}/3rdparty/lib -lpthread-google"
fi
AC_SUBST(LIB_PTHREAD_GOOGLE,$ff_pthread_google)

# FFCS - MUMPS_seq has a different Win32 compiler setup from FFCS, so we need to add some extra parameters

if test "$OS" = Windows_NT
then
    CFLAGS="$CFLAGS -DWITHOUT_PTHREAD -DAdd_"

    # we also need to satisfy ff-c++ that the pthread are not a blocking point
    if test -n "$ff_pthread"  ; then
       AC_FF_ADDWHERELIB(pthread,"",)
    fi
# ALH  MPI  "-fno-range-check" on Windows, but this options fails on MacOS 10.8.  Add no-range-check for Windows
# for hexadecimal parameter constants like:
# [[file:c:/cygwin/home/alh/ffcs/rel/mingw/mpif.h::PARAMETER MPI_SHORT_INT z 8c000003]]
#
    AC_SUBST(NO_RANGE_CHECK,-fno-range-check)
fi


# ALH - 4/9/13 - request from Helmut Jarausch - allow to change Scotch include path
if test "$with_scotch_include" = ""
then
    with_scotch_include=$ac_pwd/3rdparty/include/scotch
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
AC_SUBST(TEST_FFPPMPI,"../../src/mpi/ff-mpirun")

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
  ff_save_cppflags="$CPPFLAGS";
  CPPFLAGS="$CPPFLAGS $HDF5_CPPFLAGS"
  AC_CHECK_HEADER([H5Cpp.h],[],[with_hdf=no])
  CPPFLAGS="$ff_save_cppflags"
  #if test "$with_hdf5" = "yes"; then
  AC_FF_ADDWHERELIB(hdf5,$HDF5_LDFLAGS  $HDF5_LIBS,$HDF5_CPPFLAGS)
  LIBS="$LIBS $HDF5_LDFLAGS  $HDF5_LIBS"
  # hack, if enable_download=no
  AC_SUBST([TOOL_COMPILE_iohdf5],ifelse(,,iohdf5,))
  AC_SUBST([TOOL_DYLIB_iohdf5],iohdf5.$DYLIB_SUFFIX)
  ifelse(nodownload,,[ff_download_lib="iohdf5 $ff_download_lib"])
  ifelse(nodownload,,[enable_iohdf5_download=yes])
  #fi
else
  enable_iohdf5=no
fi
#TOOL_DISABLE(iohdf5,iohdf5.$DYLIB_SUFFIX,,nodownload)


AM_CONDITIONAL([FFCS_WINDOWS],[test "$OS" = Windows_NT])

#  remove gsl if not find ... FH
# correction FH .. 18/12/2013.
# ALH - 7/1/14 - not able to compile gsl with FFCS on Windows
if test "$enable_ffcs" != yes
then
    if test  "$ff_where_lib_conf_gsl" = 1  -a  "$enable_gsl" != "no" ; then enable_gsl=yes; fi;
fi

# hack, if enable_download=no
if test "$ff_with_gsl" = "yes"; then
  AC_SUBST([TOOL_COMPILE_gsl],ifelse(,,gsl,))
  AC_SUBST([TOOL_DYLIB_gsl],gsl.$DYLIB_SUFFIX)
  ifelse(nodownload,,[ff_download_lib="gsl $ff_download_lib"])
  ifelse(nodownload,,[enable_gsl_download=yes])  
else
  enable_gsl=no
fi

TOOL_DISABLE(NewSolver,"NewSolver.$DYLIB_SUFFIX",,nodownload)
TOOL_DISABLE(mkl,"PARDISO.$DYLIB_SUFFIX",,nodownload)

# echo "debug cxxx11: $ff_cxx11 mpi: $ff_mpi petsc: $ff_petsc_ok "
test "$enable_download" != "yes" -a "$enable_download_hpddm" != "yes" -a \( "$enable_hpddm" != "yes" -o "$ff_hpddm_ok" != "yes" \) && enable_hpddm=no
test "$enable_download" != "yes" -a "$enable_htool" != yes && enable_htool=no

test "$enable_download" != "yes" -a "$enable_bemtool" != yes && enable_bemtool=no
test "$enable_download" != "yes" -a "$enable_boost" != yes && enable_boost=no

test \( \( "$enable_download" = "yes" -o "$enable_download_hpddm" = "yes" \) -a "$enable_hpddm" != "no" \) -o "$ff_hpddm_ok" = "yes" && enable_hpddm=yes
test "$enable_download" = yes -a "$enable_htool" != no && enable_htool=yes

test "$enable_download" = yes  -a "$enable_bemtool" != no && enable_bemtool=yes
test "$enable_download" = yes -a "$enable_boost" != no && enable_boost=yes

test "$ff_cxx11" != yes -o  "$ff_mpi" != yes  -o \( "$ff_umfpack_ok" != "yes" -a "$ff_mumps_ok" != "yes" \) &&  enable_hpddm=no
test "$ff_cxx11" != yes -o  "$ff_mpi" != yes   &&  enable_htool=no
test "$ff_cxx11" != yes -o  "$ff_mpi" != yes  -o \( "$enable_bemtool" != "yes" -a "$enable_boost" != "yes" \) &&  enable_bemtool=no 

FF_PETSC_DYLIB=""
if test "$ff_petsc_ok" != no ;   then
FF_PETSC_DYLIB="PETSc.$DYLIB_SUFFIX function-PETSc.$DYLIB_SUFFIX"
fi
if test "$ff_petsccomplex_ok" != no ;   then
FF_PETSC_DYLIB="$FF_PETSC_DYLIB PETSc-complex.$DYLIB_SUFFIX"
fi
if test "$ff_slepccomplex_ok" != no ;   then
FF_PETSC_DYLIB="$FF_PETSC_DYLIB SLEPc-complex.$DYLIB_SUFFIX"
fi
if test "$ff_slepc_ok" != no ;   then
FF_PETSC_DYLIB="$FF_PETSC_DYLIB SLEPc.$DYLIB_SUFFIX"
fi

##echo @@@@@@@@@@@@ $ff_petsc_ok @@ "hpddm.$DYLIB_SUFFIX $FF_PETCS_DYLIB"
TOOL_DISABLE(hpddm,["hpddm.$DYLIB_SUFFIX hpddm_substructuring.$DYLIB_SUFFIX $FF_PETSC_DYLIB"])
TOOL_DISABLE(bem,["bem.$DYLIB_SUFFIX"])

##
if test "$ff_mkl_ok" = yes ; then
    AC_DEFINE(HAVE_MKL,1, the MKL intel lib is present for BLAS and LAPACK  )
fi

## Condtional for check
test "$ff_fftw_ok" = yes -o -n "$ff_DOWNLOAD_FFTW" && ff_fftw_=yes || ff_fftw_=no
test "$ff_ipopt_ok" = yes -o -n "$TOOL_COMPILE_ipopt" && ff_ipopt_=yes || ff_ipopt_=no
test "$ff_metis_ok" = yes -o -n "$TOOL_COMPILE_metis" && ff_metis_=yes || ff_metis_=no
test "$ff_mmg_ok" = yes -o -n "$TOOL_COMPILE_mmg" && ff_mmg_=yes || ff_mmg_=no
test "$ff_mshmet_ok" = yes -o -n "$TOOL_COMPILE_mshmet" && ff_mshmet_=yes || ff_mshmet_=no
test "$ff_mumps_ok" = yes -o -n "$TOOL_COMPILE_mumps" && ff_mumps_=yes || ff_mumps_=no
test "$ff_mumps_seq_ok" = yes -o -n "$TOOL_COMPILE_mumps_seq" && ff_mumps_seq_=yes || ff_mumps_seq_=no
test "$ff_parmetis_ok" = yes -o -n "$TOOL_COMPILE_parmetis" && ff_parmetis_=yes || ff_parmetis_=no
test "$ff_scotch_ok" = yes -o -n "$TOOL_COMPILE_scotch" && ff_scotch_=yes || ff_scotch_=no
test "$ff_tetgen_ok" = yes -o -n "$TOOL_COMPILE_tetgen" && ff_tetgen_=yes || ff_tetgen_=no
test "$ff_umfpack_ok" = yes -o -n "$TOOL_COMPILE_umfpack"  && ff_umfpack_=yes || ff_umfpack_=no
test "$ff_parmmg_ok" = yes -o -n "$TOOL_COMPILE_parmmg"  && ff_parmmg_=yes || ff_parmmg_=no

AC_SUBST(DOWNLOADCOMPILE,$DOWNLOADCOMPILE)
AM_CONDITIONAL([DOWNLOAD], [test "$enable_download" = yes])
AM_CONDITIONAL([ARPACK], [test "$ff_arpack_ok" = yes])
AM_CONDITIONAL([BLACS], [test "$ff_blacs_ok" = yes])
AM_CONDITIONAL([BEMTOOL], [test "$enable_bemtool" = yes])
AM_CONDITIONAL([BOOST], [test "$enable_boost" = yes])
AM_CONDITIONAL([FFTW3], [test "$ff_fftw_" = yes])
AM_CONDITIONAL([FORTRAN], [test "$enable_fortran" != no ])
AM_CONDITIONAL([FREEYAMS], [test "$ff_freeyams_ok" = yes -o -n "$TOOL_COMPILE_yams"])
AM_CONDITIONAL([GMM], [test "$enable_gmm" = yes -o -n "$TOOL_COMPILE_gmm"])
AM_CONDITIONAL([GSL], [test "$ff_with_gsl" = yes])
AM_CONDITIONAL([HDF5], [test "$with_hdf5" = yes])
AM_CONDITIONAL([HTOOL], [test "$enable_htool" = yes])
AM_CONDITIONAL([HPDDM], [test "$enable_hpddm" = yes])
AM_CONDITIONAL([IPOPT], [test "$ff_ipopt_" = yes])
AM_CONDITIONAL([LIBSEQ], [test "$ff_libseq_ok" = yes])
AM_CONDITIONAL([METIS], [test "$ff_metis_" = yes])
AM_CONDITIONAL([MKL], [test "$ff_mkl" = yes])
AM_CONDITIONAL([MMG3D], [test "$ff_mmg3d_ok" = yes -o -n "$TOOL_COMPILE_mmg3d"])
AM_CONDITIONAL([MMG], [test "$ff_mmg_" = yes])
AM_CONDITIONAL([PARMMG], [test "$ff_parmmg_" = yes])
AM_CONDITIONAL([MPI], [test "$ff_mpi" = yes])
AM_CONDITIONAL([MSHMET], [test "$ff_mshmet_" = yes])
AM_CONDITIONAL([MUMPS], [test "$ff_mumps_" = yes])
AM_CONDITIONAL([MUMPS_ptscotch], [test "$ff_mumps_ptscotch_ok" = yes])
AM_CONDITIONAL([MUMPS_scotch], [test "$ff_mumps_scotch_ok" = yes])
AM_CONDITIONAL([MUMPS_seq], [test "$ff_mumps_seq_" = yes])
AM_CONDITIONAL([NLOPT], [test "$ff_nlopt_ok" = yes -o -n "$TOOL_COMPILE_nlopt"])
AM_CONDITIONAL([PARMETIS], [test "$ff_parmetis_" = yes])
AM_CONDITIONAL([PETSC], [test "$ff_petsc_ok" = yes])
AM_CONDITIONAL([PETSCCOMPLEX], [test "$ff_petsccomplex_ok" = yes])
AM_CONDITIONAL([PIPE], [test -n "$TOOL_COMPILE_pipe"])
AM_CONDITIONAL([PTHREAD], [test -n "$ff_pthread"])
AM_CONDITIONAL([PTSCOTCH], [test "$ff_ptscotch_ok" = yes])
AM_CONDITIONAL([SCALAPACK], [test "$ff_scalapack_ok" = yes])
AM_CONDITIONAL([SCOTCH], [test "$ff_scotch_" = yes])
AM_CONDITIONAL([SLEPC], [test "$ff_slepc_ok" = yes])
AM_CONDITIONAL([SLEPCCOMPLEX], [test "$ff_slepccomplex_ok" = yes])
AM_CONDITIONAL([SUPERLU], [test "$ff_superlu_ok" = yes])
AM_CONDITIONAL([TETGEN], [test "$ff_tetgen_" = yes])
AM_CONDITIONAL([UMFPACK], [test "$ff_umfpack_" = yes])

##
# All makefiles
AC_OUTPUT(
    Makefile
    3rdparty/Makefile
    3rdparty/blas/Makefile
    3rdparty/arpack/Makefile
    3rdparty/umfpack/Makefile
    3rdparty/fftw/Makefile
    src/Makefile
    src/bamglib/Makefile
    src/Graphics/Makefile
    src/femlib/Makefile
    src/Algo/Makefile
    src/lglib/Makefile
    src/ffgraphics/Makefile
    src/ffgraphics/server/Makefile
    src/ffgraphics/client/Makefile
    src/fflib/Makefile
    src/nw/Makefile
    src/mpi/Makefile
    src/bamg/Makefile
    src/libMesh/Makefile
    src/medit/Makefile
    src/bin-win32/Makefile
    plugin/Makefile
    plugin/seq/Makefile
    plugin/mpi/Makefile
    examples/Makefile
    examples/plugin/Makefile
    examples/tutorial/Makefile
    examples/misc/Makefile
    examples/mpi/Makefile
    examples/hpddm/Makefile
    examples/ffddm/Makefile
    examples/eigen/Makefile
    examples/examples/Makefile
    examples/bug/Makefile
    examples/3d/Makefile
    examples/3dSurf/Makefile
    examples/3dCurve/Makefile
)

AC_MSG_NOTICE([  FreeFEM used download: $enable_download ])
AC_MSG_NOTICE([  --  Dynamic load facility: $ff_dynload ])
AC_MSG_NOTICE([  --  ARPACK (eigen value): $ff_arpack_ok ])
AC_MSG_NOTICE([  --  UMFPACK (sparse solver): $ff_umfpack_ok ])
AC_MSG_NOTICE([  --  BLAS: $ff_blas_ok ])
AC_MSG_NOTICE([  --  with MPI:     $ff_mpi])
AC_MSG_NOTICE([    --  with PETSc: $ff_petsc_ok / PETSc complex: $ff_petsccomplex_ok ]);
AC_MSG_NOTICE([    --  with SLEPc: $ff_slepc_ok / SLEPc complex: $ff_slepccomplex_ok ]);
AC_MSG_NOTICE([    --  with hpddm: $enable_hpddm ]);
AC_MSG_NOTICE([    --  with htool: $enable_htool ]);
AC_MSG_NOTICE([    --  with bemtool: $enable_bemtool (need boost: $enable_boost and htool: $enable_htool) ]);
AC_MSG_NOTICE([  --  without libs: $ff_unlib_plugin ]);
AC_MSG_NOTICE([  --  without plugins: $ff_uncompile_plugin ]);


AC_MSG_NOTICE([    progs: $ff_progs ])
if test  "$ff_umfpack_download" = yes -a  -n "$TOOL_COMPILE_umfpack" ;then
    AC_MSG_NOTICE([      use of download UMFPACK see  3rdparty/umfpack/SuiteSparse/UMFPACK/README.txt for the License])
fi
if test "$ff_arpack_download" = yes -a  -n "$TOOL_COMPILE_arpack"  ; then
    AC_MSG_NOTICE([      use of download ARPACK see  3rdparty/arpack/ARPACK/README  no License ])
fi

if test "$enable_download" = yes; then
  if test -n "$ff_DOWNLOAD_FFTW"; then
    AC_MSG_NOTICE([      use of download fftw   see  3rdparty/fftw/fftw-3.2/COPYRIGHT ])
  fi
  if test -n "$TOOL_COMPILE_tetgen"; then
    AC_MSG_NOTICE([      use of download tetgen  see 3rdparty/tetgen/tetgen1.5.1-beta/LICENSE ])
  fi
  if test -n "$TOOL_COMPILE_yams"; then
    AC_MSG_NOTICE([      use of download freeyams  see 3rdparty/yams/freeyams.2011.02.22 (suface mesh adaptation)  ])
  fi
  if test -n "$TOOL_COMPILE_mmg3d"; then
    AC_MSG_NOTICE([      use of download mmg3d (v4)   see 3rdparty/mmg3d/mmg3d4/LICENCE.txt  ])
  fi
  if test -n "$TOOL_COMPILE_mmg"; then
    AC_MSG_NOTICE([      use of download mmg  see 3rdparty/mmg/mmg-sources/LICENSE ])
  fi
  if test -n "$TOOL_COMPILE_parmmg"; then
    AC_MSG_NOTICE([      use of download parmmg  see 3rdparty/parmmg/parmmg-sources/LICENSE ])
  fi
  if test -n "$TOOL_COMPILE_parmetis"; then
    AC_MSG_NOTICE([      use of download parmetis])
  fi

  # ALH - some tools may be activated but not downloaded if a local version is specified (see
  # [[TOOL_PARAMETERS]])

  if test "$enable_superlu_download" = yes; then
    AC_MSG_NOTICE([      use of download superlu])
  fi

  if test -n "$MPICC"; then
    AC_MSG_NOTICE([      try to download: $ff_download_lib])
  fi

  if test "$enable_boost" = yes; then
    AC_MSG_NOTICE([      use of download boost        see 3rdparty/boost/LICENSE_1_0.txt  ])
  fi

  if test "$ff_download_blas" = "generic"; then
    AC_MSG_NOTICE([      use of download generic blas and cblas freefem may be slow  ])
    AC_MSG_NOTICE([      you can try to use the  Kazushige Goto s BLAS at http://www.cs.utexas.edu/users/flame/goto/  ])
    AC_MSG_NOTICE([          or at  http://www.tacc.utexas.edu/~kgoto/ for the best BLAS .])
    AC_MSG_NOTICE([      or try to 3rdparty/compile the altas blas .])
  fi

  # ALH - 30/10/13 - (request from FH) all downloads are now part of a separate script ([[file:3rdparty/getall]])
  # which needs to be run before make for the user to specify whether he is ok with each of the package licences.

  AC_MSG_NOTICE([      Please run "3rdparty/getall" to download all necessary packages before running make])
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

if test "$enable_summary" = yes; then
  if test "$enable_fortran" != yes; then
    enable_fortran=no
  fi
  if test "$enable_gmm" != yes; then
    enable_gmm=no
  fi
  if test "$ff_mkl" != yes; then
    ff_mkl=no
  fi
  if test -n "$ff_pthread"; then
    ff_pthread_ok=yes
  else
    ff_pthread_ok=no
  fi
 	

  AC_MSG_NOTICE([ Activated libraries ])
  AC_MSG_NOTICE([ - ARPACK: $ff_arpack_ok])
  AC_MSG_NOTICE([ - BLACS: $ff_blacs_ok])
  AC_MSG_NOTICE([ - BEMTOOL: $enable_bemtool])
  AC_MSG_NOTICE([ - BOOST: $enable_boost])
  AC_MSG_NOTICE([ - FFTW3: $ff_fftw_])
  AC_MSG_NOTICE([ - FORTRAN: $enable_fortran])
  AC_MSG_NOTICE([ - FREEYAMS: $ff_freeyams_ok])
  AC_MSG_NOTICE([ - GMM: $enable_gmm])
  AC_MSG_NOTICE([ - GSL: $ff_with_gsl])
  AC_MSG_NOTICE([ - HDF5: $with_hdf5])
  AC_MSG_NOTICE([ - HTOOL: $enable_htool])
  AC_MSG_NOTICE([ - HPDDM: $enable_hpddm])
  AC_MSG_NOTICE([ - IPOPT: $ff_ipopt_])
  AC_MSG_NOTICE([ - LIBSEQ: $ff_libseq_ok])
  AC_MSG_NOTICE([ - METIS: $ff_metis_])
  AC_MSG_NOTICE([ - MKL: $ff_mkl])
  AC_MSG_NOTICE([ - MMG3D: $ff_mmg3d_ok])
  AC_MSG_NOTICE([ - MMG: $ff_mmg_])
  AC_MSG_NOTICE([ - PARMMG: $ff_parmmg_])
  AC_MSG_NOTICE([ - MPI: $ff_mpi])
  AC_MSG_NOTICE([ - MSHMET: $ff_mshmet_])
  AC_MSG_NOTICE([ - MUMPS: $ff_mumps_])
  AC_MSG_NOTICE([ - MUMPS_ptscotch: $ff_mumps_ptscotch_ok])
  AC_MSG_NOTICE([ - MUMPS_scotch: $ff_mumps_scotch_ok])
  AC_MSG_NOTICE([ - MUMPS_seq: $ff_mumps_seq_])
  AC_MSG_NOTICE([ - NLOPT: $ff_nlopt_ok])
  AC_MSG_NOTICE([ - PARMETIS: $ff_parmetis_])
  AC_MSG_NOTICE([ - PETSC: $ff_petsc_ok])
  AC_MSG_NOTICE([ - PETSCCOMPLEX: $ff_petsccomplex_ok])
  AC_MSG_NOTICE([ - PTHREAD: $ff_pthread_ok])
  AC_MSG_NOTICE([ - PTSCOTCH: $ff_ptscotch_ok])
  AC_MSG_NOTICE([ - SCALAPACK: $ff_scalapack_ok])
  AC_MSG_NOTICE([ - SCOTCH: $ff_scotch_])
  AC_MSG_NOTICE([ - SLEPC: $ff_slepc_ok])
  AC_MSG_NOTICE([ - SLEPCCOMPLEX: $ff_slepccomplex_ok])
  AC_MSG_NOTICE([ - SUPERLU: $ff_superlu_ok])
  AC_MSG_NOTICE([ - TETGEN: $ff_tetgen_])
  AC_MSG_NOTICE([ - UMFPACK: $ff_umfpack_])
fi
if test "${prefix_petsc}" != no -a   ! \( -d "${prefix_petsc}/real/lib"  -o  -d "${prefix_petsc}/r/lib" \) ; then
  AC_MSG_NOTICE([ PETSc directories do not exist, to build do:])
  AC_MSG_NOTICE([cd 3rdparty/ff-petsc && make])
fi
