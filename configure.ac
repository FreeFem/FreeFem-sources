# Configuration script using Automake + Autoconf for FreeFem++
# ------------------------------------------------------------

# Antoine Le Hyaric - LJLL Paris 6 - lehyaric@ann.jussieu.fr - 13/5/04
# $Id$

AC_INIT(FreeFem++,1.41)
AC_CONFIG_SRCDIR(src/FreeFem++-CoCoa)

# Automake 1.4 is too old
AM_INIT_AUTOMAKE(1.6.3)

AM_CONFIG_HEADER(config.h)
AC_PROG_MAKE_SET

AC_PROG_RANLIB

# Necessary compilers

AC_PROG_CC
AC_PROG_CXX
AC_PROG_YACC
AM_PROG_LEX
# AC_F77_LIBRARY_LDFLAGS

# Some useful libraries

AC_CHECK_LIB(pthread,pthread_create)

# FreeFem++-specific version information
# --------------------------------------

AC_DEFINE_UNQUOTED(VersionFreeFempp,$VERSION,FreeFem++ version as a float)
FF_DATE=`date`
AC_DEFINE_UNQUOTED(VersionFreeFemDate,"$FF_DATE",FreeFem++ build date)

# Checking wether to produce a debug version
# ------------------------------------------

# ------------------------------
# compilation flags of FreeFem++ 
# ------------------------------
#  -DDRAWING  bamg active some drawing facility in bamg (inquire mesh)
#  -DLONG_LONG active the use of  long long in bamg
#  -DDEBUG active the  checking in bamg  
#  -DNCHECKPTR remove check pointer facility 
#  -DCHECK_KN active subscripting of some array  
#  -DWITHCHECK of the dynamic type ckecking of the langague (very slow)
#  -DEIGENVALUE  to compile the eigen value part
#  ------------------------------------------------

AC_ARG_ENABLE(debug,
     [  --enable-debug    Turn on debug versions of FreeFem++],
     [case "${enableval}" in
       yes) ff_debug=true ;;
       no)  ff_debug=false ;;
       *) AC_MSG_ERROR(bad value ${enableval} for --enable-debug) ;;
     esac],
     [ff_debug=false])

if test "$ff_debug" = true;
then
	CXXFLAGS="-DDRAWING -DLONG_LONG -DCHECK_KN -g"
	if test "ff_mac" = yes;
	then
		CXXFLAGS="$CXXFLAGS -fno-inline -fpascal-strings -fexceptions"
	fi
	AC_SUBST(FF_PROG_EXT,-g)
else
	CXXFLAGS="-DDRAWING -DLONG_LONG  -DNCHECKPTR -O"
	if test "ff_mac" = yes;
	then
		CXXFLAGS="$CXXFLAGS  -fpascal-strings"
	fi
fi

# Checking wether we can produce a X11 version
# --------------------------------------------

AC_MSG_CHECKING(wether we can compile for X11)
ff_x11=no
if test -d /usr/X11R6;
then
	ff_x11=yes
	AC_SUBST(STDPROG,FreeFem++${FF_PROG_EXT}${EXEEXT})
	AC_SUBST(X11PROG,FreeFem++-x11${FF_PROG_EXT}${EXEEXT})
	AC_SUBST(INCLUDEX11,/usr/X11R6/include)
	AC_SUBST(LIBX11,"-L/usr/X11R6/lib -lX11 -lm")
fi
AC_MSG_RESULT($ff_x11)

# Checking wether we can produce a MacIntosh-specific version
# -----------------------------------------------------------

AC_MSG_CHECKING(wether we are on a MacIntosh)
ff_mac=no
ff_aglprog=
if test `uname` = Darwin;
then
	ff_mac=yes
	ff_aglprog=FreeFem++-agl${FF_PROG_EXT}${EXEEXT}
fi
AC_MSG_RESULT($ff_mac)

# Checking that OpenGL is accessible
# ----------------------------------

# If we cannot find gl.h, no OpenGL version can be produced
AC_CHECK_HEADER(GL/gl.h,
	ff_glxprog=FreeFem++-glx${FF_PROG_EXT}${EXEEXT},
	ff_glxprog=
	ff_aglprog=
)

# Looking for the BLAS
# --------------------

# Trying default locations
ff_blas_ok=no
ff_save_libs="$LIBS"
AC_CHECK_LIB(blas,daxpy_,
	LIBS="$LIBS -lblas"
	ff_blas_ok=yes,
	LIBS="$ff_save_libs")

# Trying Atlas (location corresponding to Debian package)
if test "$ff_blas_ok" = no;
then
	AC_MSG_CHECKING(for daxpy_ in /usr/lib/atlas/libblas.a)
	ff_save_libs="$LIBS"
	LIBS="$LIBS -L/usr/lib/atlas -lblas"
	AC_LINK_IFELSE(
		[AC_LANG_CALL(,daxpy_)],
		ff_blas_ok=yes,
		LIBS="$ff_save_libs")
	AC_MSG_RESULT($ff_blas_ok)
fi

# Trying  freebsd framework (on darwin/ MacOS X)
if test "$ff_blas_ok" = no;
then
        AC_MSG_CHECKING(for daxpy_ in framework vecLib)
        ff_save_libs="$LIBS"
        LIBS="$LIBS -framework vecLib"
        AC_LINK_IFELSE(
                [AC_LANG_CALL(,daxpy_)],
                ff_blas_ok=yes,
                LIBS="$ff_save_libs")
        AC_MSG_RESULT($ff_blas_ok)
fi

# Looking for ARPACK
# ------------------

# Trying default locations
#
ff_arpack_ok=no   
ff_save_libs="$LIBS"
ff_arpack_libs=
if test "$ff_blas_ok" = yes;
   then	
   AC_CHECK_LIB(arpack,dsaupd_,
	ff_arpack_libs="-larpack -llapack -rdynamic"
        ff_arpack_ok=yes,,
	-llapack -rdynamic)

   # Trying to "locate" Arpack
   if test "$ff_arpack_ok" = no;
   then
        AC_MSG_CHECKING(for libarpack with locate)
        ff_lib_arpack=`locate libarpack|grep  'libarpack.*.a$'|head -1`" -L/sw/lib -lg2c"
	LIBS="$ff_lib_arpack $LIBS"
        AC_LINK_IFELSE(
                [AC_LANG_CALL(,dsaupd_)],
                ff_arpack_ok=yes
                ff_arpack_libs=$ff_lib_arpack)
        AC_MSG_RESULT($ff_arpack_ok)
    fi
fi
LIBS="$ff_save_libs $ff_arpack_libs"

if test "$ff_arpack_ok" = yes;
then
	EIGENOBJ='eigenvalue.$(OBJEXT)'
	AC_DEFINE(HAVE_LIBARPACK,1,
		Arpack is used for eigenvalue computation)
fi
AC_SUBST([EIGENOBJ])

# Looking for UMFPACK
# -------------------

# Trying default locations

if test "$ff_blas_ok" = yes;
then
	AC_CHECK_HEADERS(umfpack.h umfpack/umfpack.h)
	AC_CHECK_LIB(amd,amd_info)
	AC_CHECK_LIB(umfpack,umf_i_malloc)
fi

# Checking wether we can produce a parallel version
# -------------------------------------------------

AC_CHECK_PROG(ff_mpi,mpiCC,yes,no)
if test "$ff_mpi" = yes;
then

	# We do not AC_DEFINE any special flag for parallel
	# computation here, because it must only be set when the
	# parallel program is compiled (see src/mpi/Makfile.am)

	AC_SUBST(MPIPROG,FreeFem++-mpi${FF_PROG_EXT}${EXEEXT})
fi

# Checking wether we can generate some documentation
# --------------------------------------------------

AC_CHECK_PROG(ff_latex,latex,yes,no)
AC_CHECK_PROG(ff_makeindex,makeindex,yes,no)
AC_CHECK_PROG(ff_dvips,dvips,yes,no)
if test "$ff_latex" = yes -a "$ff_makeindex" = yes -a "$ff_dvips" = yes;
then
	AC_SUBST(DOCPS,manual.ps)

	AC_CHECK_PROG(ff_gzip,gzip,yes,no)
	if test "$ff_gzip" = yes;
	then
		AC_SUBST(DOCPSGZ,manual.ps.gz)
	fi
fi

AC_CHECK_PROG(ff_pdflatex,pdflatex,yes,no)
if test "$ff_pdflatex" = yes;
then
	AC_SUBST(DOCPDF,manual.pdf)
fi

# creating all makefiles
# ----------------------

AC_SUBST(GLXPROG,$ff_glxprog)
AC_SUBST(AGLPROG,$ff_aglprog)
AC_OUTPUT(Makefile
	src/Makefile
	src/bamglib/Makefile
	src/femlib/Makefile
	src/Algo/Makefile
	src/lglib/Makefile
	src/fflib/Makefile
	src/std/Makefile
	src/x11/Makefile
	src/nw/Makefile
	src/mpi/Makefile
	src/agl/Makefile
	src/glx/Makefile
	examples++/Makefile
	examples++-eigen/Makefile
	examples++-tutorial/Makefile
	DOC/Makefile,
        echo timestamp > stamp-h)
