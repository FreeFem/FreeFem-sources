# Configuration script using Automake + Autoconf for FreeFem++
# ------------------------------------------------------------

# Antoine Le Hyaric - LJLL Paris 6 - lehyaric@ann.jussieu.fr - 13/5/04
# $Id$

AC_INIT(FreeFem++,1.41)
AC_CONFIG_SRCDIR(src/FreeFem++-CoCoa)
AM_INIT_AUTOMAKE
AM_CONFIG_HEADER(config.h)
AC_PROG_MAKE_SET

# Find out the architecture. This will be used as a name for the
# directory containing objet files for each different architecture.

AC_CANONICAL_BUILD

# Necessary compilers

AC_PROG_CXX
AC_PROG_YACC
AM_PROG_LEX
# AC_F77_LIBRARY_LDFLAGS

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


# Nothing can be compiled without the BLAS
if test "$ff_blas_ok" = no;
then
	AC_MSG_ERROR(Cannot find BLAS library)
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
	
   # Trying "locate"
   if test "$ff_arpack_ok" = no;
   then
        AC_MSG_CHECKING(for libarpack with locate)
        lib_arpack=`locate libarpack|grep  'libarpack.*.a$'|head -1`" -L/sw/lib -lg2c"
	AC_MSG_CHECKING("----------------  $lib_arpack");
	LIBS="$lib_arpack $LIBS"
        AC_LINK_IFELSE(
                [AC_LANG_CALL(,dsaupd_)],
                ff_arpack_ok=yes
                AC_SUBST(arpacklib,$lib_arpack))
        AC_MSG_RESULT($ff_arpack_ok)
    fi
fi
LIBS=$ff_save_libs
# Looking for the ARPACK                                                                                                      

# Checking wether we can produce a parallel version
# -------------------------------------------------

AC_SUBST(mpiprogs,[])
AC_CHECK_PROG(ff_mpi,mpiCC,yes,no)
if test "$ff_mpi" = yes;
then
	AC_SUBST(mpiprogs,mpi mpi-g)
fi

# Creating all makefiles
# ----------------------

AC_OUTPUT(Makefile
	src/Makefile
	DOC/Makefile,
        echo timestamp > stamp-h)
