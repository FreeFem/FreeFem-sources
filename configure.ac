# Configure script using Automake + Autoconf for FreeFem++
# --------------------------------------------------------

# Antoine Le Hyaric - LJLL Paris 6 - lehyaric@ann.jussieu.fr - 13/5/04
# $Id$

AC_INIT(FreeFem++,1.41)
AC_CONFIG_SRCDIR(src/FreeFem++-CoCoa)
AM_INIT_AUTOMAKE
AM_CONFIG_HEADER(config.h)
AC_PROG_MAKE_SET

# To locate objet files for different architectures in separate
# directories

AC_CANONICAL_BUILD

# Compilers

AC_PROG_CXX
AC_PROG_YACC
AM_PROG_LEX

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
	ff_save_libs="$LIBS"
	LIBS="$LIBS -L/usr/lib/atlas"
	AC_CHECK_LIB(blas,daxpy_,
		LIBS="$LIBS -lblas"
		ff_blas_ok=yes,
		LIBS="$ff_save_libs")
fi

AC_OUTPUT(Makefile
	src/Makefile
	DOC/Makefile,
        echo timestamp > stamp-h)
