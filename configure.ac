# Configuration script using Automake + Autoconf for FreeFem++
# ------------------------------------------------------------

# Antoine Le Hyaric - LJLL Paris 6 - lehyaric@ann.jussieu.fr - 13/5/04
# $Id$

AC_INIT(FreeFem++,1.41,hecht@ann.jussieu.fr,freefem++)
AC_CONFIG_SRCDIR(src/FreeFem++-CoCoa)

# Automake 1.4 is too old
AM_INIT_AUTOMAKE(1.6.3 dist-zip)

AM_CONFIG_HEADER(config.h)
AC_PROG_MAKE_SET

AC_PROG_RANLIB

# Necessary compilers

AC_PROG_F77
AC_PROG_CC
AC_PROG_CXX
AC_PROG_YACC
AM_PROG_LEX
AC_LANG(C++)

# We do not use the default C++ flags
CXXFLAGS=

# Some useful libraries

AC_CHECK_LIB(pthread,pthread_create)

# Necessary absolute pathname for local directory when some libraries
# are used from several different locations (for instance locally in
# configure.ac and in a subdir).
curdir=`pwd`

# Checking wether we can produce a MacIntosh-specific version
# -----------------------------------------------------------

AC_MSG_CHECKING(wether we are on a MacIntosh)
ff_mac=no
ff_aglprog=
if test `uname` = Darwin;
then
	ff_mac=yes
	ff_aglprog=FreeFem++-agl${EXEEXT}
fi
AC_MSG_RESULT($ff_mac)

# Checking wether we can produce a Microsoft Windows-specific version
# -------------------------------------------------------------------

AC_MSG_CHECKING(wether we are on Microsoft Windows)
ff_libdl="-ldl"
ff_uname=`uname`
case $ff_uname in
CYGWIN*)
	ff_win32=yes

	# No need for libdl under Windows
	ff_libdl=

	# We need Mingw (to avoid Cygwin's extra DLLs)
	CXXFLAGS="$CXXFLAGS -mno-cygwin -I/usr/include/mingw"
	LIBS="$LIBS -mno-cygwin -L/lib/mingw"

	;;
*)ff_win32=no;;
esac
AC_MSG_RESULT($ff_win32)
AC_SUBST(LIBDL,$ff_libdl)

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
     --enable-debug	Turn on debug versions of FreeFem++,
     [case "${enableval}" in
       yes) ff_debug=true ;;
       no)  ff_debug=false ;;
       *) AC_MSG_ERROR(bad value ${enableval} for --enable-debug) ;;
     esac],
     [ff_debug=false])

if test "$ff_debug" = true;
then
	CXXFLAGS="$CXXFLAGS -DDRAWING -DLONG_LONG -DCHECK_KN -g"
	if test "ff_mac" = yes;
	then
		CXXFLAGS="$CXXFLAGS -fno-inline -fpascal-strings -fexceptions"
	fi
else
	CXXFLAGS="$CXXFLAGS -DDRAWING -DLONG_LONG  -DNCHECKPTR -O"
	if test "ff_mac" = yes;
	then
		CXXFLAGS="$CXXFLAGS  -fpascal-strings"
	fi
fi

# Checking wether we can produce a X11 version
# --------------------------------------------

AC_MSG_CHECKING(wether we can compile for X11)
ff_x11=no

# At the moment, we cannot compile for X11 on Mingw because /usr/X11R6
# is part of Cygwin, not Mingw.

if test -d /usr/X11R6 -a "$ff_win32" = no;
then
	ff_x11=yes
	ff_x11prog=FreeFem++-x11${EXEEXT}
	AC_SUBST(INCLUDEX11,-I/usr/X11R6/include)
	AC_SUBST(LIBX11,"-L/usr/X11R6/lib -lX11 -lm")
fi
AC_MSG_RESULT($ff_x11)

# Checking that OpenGL is accessible
# ----------------------------------

# If we cannot find gl.h, no OpenGL version can be produced
AC_CHECK_HEADER(GL/gl.h,
	ff_glxprog=FreeFem++-glx${EXEEXT},
	ff_glxprog=
	ff_aglprog=
)

# If we cannot find glx.h, "FreeFem++-glx" cannot be produced
AC_CHECK_HEADER(GL/glx.h,,ff_glxprog=)

# Checking wether we can produce a parallel version
# -------------------------------------------------

AC_MSG_CHECKING(for mpiCC)
ff_save_cxx="$CXX"
CXX=mpiCC
AC_LINK_IFELSE(
[#include <stdio.h>
#include <mpi.h>
int main(int argc,char **argv){
  char name[[BUFSIZ]];
  int length;
  
  MPI_Init(&argc, &argv);
  MPI_Get_processor_name(name, &length);
  printf("%s: hello world\n", name);
  MPI_Finalize();
  return 0;
}],
ff_mpi=yes,
ff_mpi=no)
CXX="$ff_save_cxx"
AC_MSG_RESULT($ff_mpi)

if test "$ff_mpi" = yes;
then

	# We do not AC_DEFINE any special flag for parallel
	# computation here, because it must only be set when the
	# parallel program is compiled (see src/mpi/Makfile.am)

	AC_SUBST(MPIPROG,FreeFem++-mpi${EXEEXT})
fi

# Looking for the BLAS
# --------------------

# Trying default locations
ff_blas_ok=no
ff_save_libs="$LIBS"
AC_CHECK_LIB(blas,daxpy_,
	LIBS="$LIBS -lblas"
	ff_blas_ok=yes
	ff_blas_lib="-lblas",
	LIBS="$ff_save_libs")

# Trying Atlas (location corresponding to Debian package)
if test "$ff_blas_ok" = no;
then
	AC_MSG_CHECKING(for daxpy_ in /usr/lib/atlas/libblas.a)
	ff_save_libs="$LIBS"
	LIBS="$LIBS -L/usr/lib/atlas -lblas"
	AC_LINK_IFELSE(
		[AC_LANG_CALL(,daxpy_)],
		ff_blas_ok=yes
		ff_blas_lib="-L/usr/lib/atlas -lblas",
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
                ff_blas_ok=yes
		ff_blas_lib="-framework vecLib",
                LIBS="$ff_save_libs")
        AC_MSG_RESULT($ff_blas_ok)
fi

# If all else fails, download!
if test "$ff_blas_ok" = no;
then
	AC_MSG_NOTICE(using downloaded Blas)
	AC_SUBST(DOWNLOAD_BLAS,blas)
	ff_blas_lib="${curdir}/download/ATLAS/lib/ff++/libatlas.a"
	LIBS="$ff_blas_lib $LIBS"
	ff_blas_ok=yes
fi
AC_SUBST(BLASLIB,$ff_blas_lib)

# Looking for ARPACK
# ------------------

# We need the following g77 libraries to connect to the Fortran 77
# Arpack.

if test "$ff_blas_ok" = yes;
then
	ff_g2c_ok=no
	AC_CHECK_LIB(g2c,G77_second_0,
		LIBS="$LIBS -lg2c"
		ff_g2c_ok=yes)

	if test "$ff_g2c_ok" = no;
	then
	        AC_MSG_CHECKING(for G77_second_0 in /sw/lib/libg2c.a)
		ff_save_libs="$LIBS"
		LIBS="$LIBS -L/sw/lib -lg2c"
		AC_LINK_IFELSE(
			[AC_LANG_CALL(,G77_second_0)],
			ff_g2c_ok=yes,
			LIBS="$ff_save_libs")
		AC_MSG_RESULT($ff_g2c_ok)
	fi
fi

# Arpack itself

ff_arpack_ok=no   
ff_save_libs="$LIBS"
ff_arpack_libs=
if test "$ff_blas_ok" = yes;
then

	# User-specified location
	AC_ARG_WITH(arpack,
		--with-arpack=archive	Use a specific version of Arpack,
		ff_arpack_ok=yes
		ff_arpack_libs="${withval}")

	# Default locations
	if test "$ff_arpack_ok" = no;
	then
		AC_CHECK_LIB(arpack,dsaupd_,
			ff_arpack_libs="-larpack -llapack -rdynamic"
			ff_arpack_ok=yes,,
			-llapack -rdynamic)
	fi

   # Trying to "locate" Arpack
   if test "$ff_arpack_ok" = no;
   then
        AC_MSG_CHECKING(for libarpack with locate)
        ff_lib_arpack=`locate libarpack|grep  'libarpack.*.a$'|head -1`
	LIBS="$ff_lib_arpack $LIBS"
        AC_LINK_IFELSE(
                [AC_LANG_CALL(,dsaupd_)],
                ff_arpack_ok=yes
                ff_arpack_libs=$ff_lib_arpack)
        AC_MSG_RESULT($ff_arpack_ok)
    fi

	# If all else fails, download!
	if test "$ff_arpack_ok" = no;
	then
		AC_MSG_NOTICE(using downloaded Arpack)
		AC_SUBST(DOWNLOAD_ARPACK,arpack)
		ff_arpack_libs="${curdir}/download/arpack/ARPACK/libarpack_ff++.a"
		ff_arpack_ok=yes
	fi
fi
LIBS="$ff_arpack_libs $ff_save_libs"

if test "$ff_arpack_ok" = yes;
then
	EIGENOBJ='eigenvalue.$(OBJEXT)'
	AC_DEFINE(HAVE_LIBARPACK,1,Arpack is used for eigenvalue computation)
fi
AC_SUBST([EIGENOBJ])

# Looking for UMFPACK
# -------------------

ff_amd_ok=no
ff_umfpack_ok=no
if test "$ff_blas_ok" = yes;
then
        # User-specified location
        AC_ARG_WITH(umfpack,
                --with-umfpack=archive	Use a specific version of Umfpack,
                ff_umfpack_ok=yes
                LIBS="${withval} $LIBS"
		AC_DEFINE(HAVE_LIBUMFPACK,1,
			Umfpack is used for sparse matrices computations))

        AC_ARG_WITH(amd,
                --with-amd=archive	Use a specific version of AMD,
                ff_amd_ok=yes
                LIBS="${withval} $LIBS")

	# Trying default locations

	AC_CHECK_HEADERS(umfpack.h umfpack/umfpack.h)

	if test "$ff_amd_ok" = no;
	then
		AC_CHECK_LIB(amd,amd_info,
			AC_DEFINE(HAVE_LIBAMD,1,AMD for UMFPACK)
			ff_amd_ok=yes)
	fi

	if test "$ff_umfpack_ok" = no -a "$ff_amd_ok" = yes;
	then
		AC_CHECK_LIB(umfpack,umf_i_malloc,
			AC_DEFINE(HAVE_LIBUMFPACK,1,UMFPACK)
			ff_umfpack_ok=yes)
	fi
fi

# If all else fails, download!
if test "$ff_umfpack_ok" = no;
then
	AC_MSG_NOTICE(using downloaded Umfpack)
	AC_SUBST(DOWNLOAD_UMFPACK,umfpack)
	LIBS="${curdir}/download/umfpack/UMFPACKv4.3/UMFPACK/Lib/libumfpack.a ${curdir}/download/umfpack/UMFPACKv4.3/AMD/Lib/libamd.a $LIBS"
	AC_DEFINE(HAVE_LIBUMFPACK,1,UMFPACK)
	ff_umfpack_ok=yes
fi

# Checking for some functions that may not appear everywhere
# ----------------------------------------------------------

# asinh acosh atanh are not in Mingw yet
AC_CHECK_FUNCS(asinh acosh atanh)

# Availability of dlopen(). Use AC_COMPILE rather than
# AC_CHECK_HEADERS the latter has problems seeing it (in Cygwin) when
# it does not compile (in Mingw).

AC_MSG_CHECKING(for dlfcn.h)
AC_COMPILE_IFELSE([#include <dlfcn.h>],
	ff_dlfcn=yes,
	ff_dlfcn=no)
AC_MSG_RESULT($ff_dlfcn)
if test "$ff_dlfcn" = yes;
then
	AC_DEFINE(HAVE_DLFCN_H,1,Dynamic loading - not mandatory)
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

# Choosing compilation options for the standard version (in src/std)
# ------------------------------------------------------------------

# The "standard" configured version can use win32 (mingw) or x11

if test "$ff_win32" = yes;
then
	ff_std_graph_obj=Pcrgraph.$OBJEXT
	ff_std_x11_lib=
	ff_std_x11_inc=
	ff_std_ldflags=-mwindows
else
	ff_std_graph_obj=Xrgraph.$OBJEXT
	ff_std_x11_lib=$LIBX11
	ff_std_x11_inc=$INCLUDEX11
	ff_std_ldflags=

	# If the standard version already uses x11, we do not need to
	# compile the src/x11 directory.
	ff_x11prog=
fi
AC_SUBST(STD_GRAPH_OBJ,$ff_std_graph_obj)
AC_SUBST(STD_X11_LIB,$ff_std_x11_lib)
AC_SUBST(STD_X11_INC,$ff_std_x11_inc)
AC_SUBST(STD_LDFLAGS,$ff_std_ldflags)

# creating all makefiles
# ----------------------

AC_SUBST(GLXPROG,$ff_glxprog)
AC_SUBST(AGLPROG,$ff_aglprog)
AC_SUBST(X11PROG,$ff_x11prog)
AC_OUTPUT(Makefile
	download/Makefile
	download/blas/Makefile
	download/arpack/Makefile
	download/umfpack/Makefile
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
	examples++-mpi/Makefile
	examples++-load/Makefile
	examples++-bug/Makefile
	DOC/Makefile,
        echo timestamp > stamp-h)
