# Configuration script using Automake + Autoconf for FreeFem++
# ------------------------------------------------------------

# Antoine Le Hyaric - LJLL Paris 6 - lehyaric@ann.jussieu.fr - 13/5/04
# $Id$

AC_INIT(FreeFem++,1.41)
AC_CONFIG_SRCDIR(src/FreeFem++-CoCoa)
AM_INIT_AUTOMAKE
AM_CONFIG_HEADER(config.h)
AC_PROG_MAKE_SET

# Necessary compilers

AC_PROG_CC
AC_PROG_CXX
AC_PROG_YACC
AM_PROG_LEX
# AC_F77_LIBRARY_LDFLAGS

# We do not need the default optimization options that autoconf
# chooses

CXXFLAGS=

# Some useful libraries

AC_CHECK_LIB(pthread,pthread_create)


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
AC_SUBST(arpacklib,[])
if test "$ff_blas_ok" = yes;
   then	
   AC_CHECK_LIB(arpack,dsaupd_,
	AC_SUBST(arpacklib,-larpack)
        ff_arpack_ok=yes)
	
   # Trying to "locate" Arpack
   if test "$ff_arpack_ok" = no;
   then
        AC_MSG_CHECKING(for libarpack with locate)
        ff_lib_arpack=`locate libarpack|grep  'libarpack.*.a$'|head -1`" -L/sw/lib -lg2c"
	LIBS="$ff_lib_arpack $LIBS"
        AC_LINK_IFELSE(
                [AC_LANG_CALL(,dsaupd_)],
                ff_arpack_ok=yes
                AC_SUBST(arpacklib,$ff_lib_arpack))
        AC_MSG_RESULT($ff_arpack_ok)
    fi
fi
LIBS=$ff_save_libs

# Checking wether we can produce a parallel version
# -------------------------------------------------

AC_CHECK_PROG(ff_mpi,mpiCC,yes,no)
AM_CONDITIONAL(WITHMPI,test "$ff_mpi" = yes)

# Checking wether we can produce an optimized version
# ---------------------------------------------------

AC_ARG_ENABLE(ff_optim,
     [  --enable-optim    Turn on optimized versions of FreeFem++],
     [case "${enableval}" in
       yes) ff_optim=true ;;
       no)  ff_optim=false ;;
       *) AC_MSG_ERROR(bad value ${enableval} for --enable-optim) ;;
     esac],
     [ff_optim=true])
AM_CONDITIONAL(OPTIM,test "$ff_optim" = true)

# Creating all makefiles
# ----------------------

AC_OUTPUT(Makefile
	src/Makefile
	DOC/Makefile,
        echo timestamp > stamp-h)
