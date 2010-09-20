g# Configuration script using Automake + Autoconf for FreeFem++
# ------------------------------------------------------------

# Antoine Le Hyaric - LJLL Paris 6 - lehyaric@ann.jussieu.fr - 13/5/04
# $Id: configure.ac,v 1.338 2010/05/10 21:35:07 hecht Exp $

# Version numbering: x.xx-pp where "pp" is the package version (when
# the same FreeFem++ version is packaged several times). For
# coherency, this should have the same value as the top-most package
# number in debian/changelog.
# add check of FLIBS and X11_LIBS  v2.24-2

AC_INIT(FreeFem++,3.9-3,hecht@ann.jussieu.fr,freefem++)
AC_PREREQ(2.50) dnl for AC_LANG_CASE
AC_CONFIG_SRCDIR(src/FreeFem++-CoCoa)
# Automake 1.4 is too old
AM_INIT_AUTOMAKE(1.6.3 dist-zip)

AM_CONFIG_HEADER(config.h)
AC_PROG_MAKE_SET

AC_PROG_RANLIB

ff_where_lib_conf=examples++-load/WHERE_LIBRARY-config
ff_HOSTARCH_pastix=
echo "#  Build with freefem++ with ./configure " `date`  >$ff_where_lib_conf
# To allow anonymous CVS version to contain a "./configure" and
# Makefiles

AM_MAINTAINER_MODE

echo >config_LIB_INFO

# Necessary compilers
AC_PROG_CC
AC_PROG_CXX




if test -z "$CXX"; then
  AC_MSG_ERROR([ Sorry no c++ compiler !])
  AC_MSG_NOTICE( [ fatal error : c++ compiler ! ] );
fi

AC_COMPUTE_INT(ff_size_long,[sizeof(long])
AC_COMPUTE_INT(ff_size_int,[sizeof(int)])
AC_COMPUTE_INT(ff_size_ptr,[sizeof(int*)])

AC_SUBST(SIZEOF_LONG,$ff_size_long)
AC_SUBST(SIZEOF_INT,$ff_size_int)
AC_SUBST(SIZEOF_PTR,$ff_size_ptr)
if test "$ff_size_ptr" -eq 4  ;then
 AC_SUBST(SIZEOF_PTRINBIT,32)
elif test  "$ff_size_ptr" -eq 8  ;then
 AC_SUBST(SIZEOF_PTRINBIT,64)
else
  AC_MSG_ERROR([ Sorry no sizeof c++ pointer are not  4 or 8])
  AC_MSG_NOTICE( [ fatal error : sizeof pointer ! ] ); 
fi


# dur dur car sous MacOsX le fortran n'est pas standard. 
ff_AR="ar"
ff_ARFLAGS="rv"
ff_RANLIB="ranlib" 

AC_ARG_ENABLE(fortran,[  --disable-fortran	No Fortran compiler available ( ARPACK need it)])
ff_g2c_lib="";	 
if test "$enable_fortran" != no
then
    AC_PROG_F77(g77 f77 gfortran g77 f90 xlf xlf90 fort77)
    AC_PROG_FC(gfortran  f90  xlf90 )
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
            for i in $FLIBS; do
		case $i in
		    -L*)  
			d=`expr $i : '-L\(.*\)'`;
			if test -e "$d/libgfortran.so" ; then
			    ff_flibs="$d/libgfortran.so"
			elif test -e "$d/libgfortran.dylib" ; then
			    ff_flibs="$d/libgfortran.dylib"
			elif test -e "$d/libgfortran.a" ; then
			    ff_flibs="$d/libgfortran.a"
			elif test -e $d/libgfortran.so -o -e $d/libgfortran.a -o -e $d/libgfortran.dylib ; then   
			    ff_flibs="$i" 
			fi;;
		esac
            done
	    if test -e "$ff_flibs"  ;then
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
	    AC_MSG_NOTICE( [add ---without-fortran  ] );
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
    AC_RUN_IFELSE([[
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
 ]],ff_f77intsize=`cat conftest.out`,ff_f77intsize=4,ff_f77intsize=4)
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
    if test $ff_f77intsize -eq  8
    then
	AC_MSG_ERROR([ Fatal Error: Your fortran 77  compiler have by 8 bytes integer ],1);  
    fi
    [echo f77 LD  $FLIBS  >>$ff_where_lib_conf]
    [echo f77 INCLUDE $ff_cfagsf77  >>$ff_where_lib_conf]
#  bof bof -- FH ... 
    [echo fc LD  $FLIBS  >>$ff_where_lib_conf]
    [echo fc INCLUDE $ff_cfagsf77  >>$ff_where_lib_conf]

# AC_DEFINE_UNQUOTED(FF_F77_INTEGER_SIZE,$ff_f77intsize,[ size of F77 integer ])
fi
#  fin test FORTRAN ..........
# -----------------------------
AC_ARG_ENABLE(c,[  --disable-c		No C compiler available (C BLAS need it)])
if test "$enable_c" != no
then
    AC_PROG_CC
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

# Necessary absolute pathname for local directory when some libraries
# are used from several different locations (for instance locally in
# configure.ac and in a subdir).
curdir=`pwd`

# Configure options
# -----------------

echo $ac_configure_args > configure.param
# by default the suffix of the  .so lib file is .so
#      it is .dll   under windows
#      it is .dylib under macos X
ff_suffix_dylib="so";	
# Checking wether we can produce a MacIntosh-specific version
# -----------------------------------------------------------

AC_MSG_CHECKING(wether we are on a MacIntosh)
ff_mac=no
ff_aglprog=
if test `uname` = Darwin;
then
	ff_HOSTARCH_pastix=i686_mac
	ff_suffix_dylib="dylib"
	ff_mac=yes
#	ff_aglprog=FreeFem++-agl${EXEEXT}
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

# Checking wether we can produce a Microsoft Windows-specific version
# -------------------------------------------------------------------

# Sometimes (under cygwin for instance) the default version of FLTK
# does not fit our purpose

AC_ARG_ENABLE(default-fltk,
	[  --disable-default-fltk	Does not use default FLTK])

AC_ARG_ENABLE(cygwindll,
	[  --enable-cygwindll	Forces the use of the Cygwin DLL (not recommended)])

# Default FLTK version under Cygwin does not use threads. So we cannot
# use it.

if test "$enable_cygwindll" = yes
then
	enable_default_fltk=no
fi

AC_MSG_CHECKING(wether we are on Microsoft Windows)
ff_uname=`uname`
#  flag to build window 32   version   ff_mingw = yes 
ff_mingw=no
case $ff_uname in
    CYGWIN*)	
	ff_nocygwin=-mno-cygwin;
	AC_SUBST(GCCNOCYGWIN,$ff_nocygwin);;
    MINGW*)
        enable_cygwindll=no;;
esac
case $ff_uname in
CYGWIN*|MINGW*)
	ff_suffix_dylib="dll";
        ff_win32=yes;
	AC_SUBST(FFMETIS_CFLAGS,"-D__VC__")
	# We need Mingw to avoid Cygwin's extra DLLs
	if test "$enable_cygwindll" != yes
	then
	        ff_glut_ok=yes
                ff_mingw=yes 
                enable_cygwindll=no;
		ff_pthread="-mthread"
		CXXFLAGS="$CXXFLAGS $ff_nocygwin -I/usr/include/mingw"		
		FFLAGS="$FFLAGS $ff_nocygwin"
		CFLAGS="$CFLAGS $ff_nocygwin -I/usr/include/mingw"
                AC_COMPILE_IFELSE([int a;],[],[AC_MSG_ERROR([Sorry $ff_nocygwin, but try with gcc-3.3])])
		LIBS="$LIBS $ff_nocygwin -mthreads -lws2_32 -lcomdlg32" 
		LIBSNOCONSOLE="-mwindows"
		FLIBS="$ff_nocygwin -lg2c"
                ff_glutname="glut32"
		#  check abort a existing function just to find in glut32.dll exist in the path 
		#  because glutInit is not the real symbol on win32 dur dur FH !!!!!!!!!
		AC_CHECK_LIB(glut32,abort,ff_glut="-l$ff_glutname -mthreads -lglu32 -lopengl32",ff_glut="")
		if test -z "$ff_glut" ; then
                ff_glutname="freeglut"
		AC_CHECK_LIB(freeglut,abort,ff_glut="-l$ff_glutname -mthreads -lglu32 -lopengl32",ff_glut_ok=no)
                fi
         	# Resources for FreeFem++-cs in Microsoft Windows format
	        AC_SUBST(FFGLUTNAME,$ff_glutname)
	        AC_SUBST(WINDRESOBJ,windres.o)
	        AC_SUBST(LIBSNOCONSOLE,$LIBSNOCONSOLE)
                AC_SUBST(WIN32DLLTARGET,win32-dll-target)
                AC_DEFINE(PURE_WIN32,1,A pure windows applications no cygwin dll)
		ff_dynload=yes
	fi


	;;
*)ff_win32=no;;
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
	CXXFLAGS="$CXXFLAGS -DDRAWING -DBAMG_LONG_LONG -DCHECK_KN"
	if test "$ff_mac" = yes;
	then
		CXXFLAGS="$CXXFLAGS -fno-inline  -fexceptions"
	fi
else
	CXXFLAGS="$CXXFLAGS -DDRAWING -DBAMG_LONG_LONG  -DNCHECKPTR"
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

# Checking wether we can produce a X11 version
# --------------------------------------------

AC_MSG_CHECKING(wether we can compile for X11)
ff_x11=no

# At the moment, we cannot compile for X11 on Mingw because /usr/X11R6
# is part of Cygwin, not Mingw.

if test "$ff_mingw" = no && test "$with_x" != "no";
then
	ff_x11=yes
fi

AC_MSG_RESULT($ff_x11)

if test "$ff_x11" = yes
then
   ff_save_libs="$LIBS"
   AC_PATH_X
   AC_PATH_XTRA
   if test "$no_x" = "yes" ;then
      ff_x11_libs="";
      ff_x11=no;
   else
	CPPFLAGS="$CPPFLAGS $X_CFLAGS"
#     echo "X11 --   $X_LIBS ext  --   $X_EXTRA_LIBS -- $X_PRE_LIBS "
	# When libXxf86vm is present, it is necessary to add it to define
	# "XF86VidModeQueryVersion" (called by OpenGL) on newer X11 systems
	LIBS="$LIBS $X_LIBS"
	AC_CHECK_LIB(Xxf86vm,
		XF86VidModeQueryVersion,
		ff_Xxf86vm="-lXxf86vm",
		ff_Xxf86vm="")
        AC_CHECK_LIB(Xpm,
                XpmReadFileToImage,
                ff_Xpm="-lXpm",
                ff_Xpm86vm="")
	LIBS="$ff_save_libs"

	ff_x11_libs="$X_LIBS $ff_Xxf86vm -lXext  -lX11 $ff_Xpm  -lm"
	AC_MSG_CHECKING(for X11_LIBS)
	AC_MSG_RESULT(["$ff_x11_libs"])	
#	ff_x11_libs="-L/usr/X11R6/lib $ff_Xxf86vm -lXext -lXpm -lX11 -lm"
	AC_SUBST(X11_LIBS,$ff_x11_libs)
	ff_x11prog=FreeFem++-x11${EXEEXT}

#    check if the X11_LIBS is correct
     ff_libs="$LIBS"
     LIBS="$ff_libs $X11_LIBS"
     AC_TRY_LINK_FUNC(XOpenDisplay,ff_err=,ff_err=ok);
# ?????
     if test "$ff_err" = "ok" ; then  AC_MSG_ERROR([ Fatal X11_LIBS: $X11_LIBS is incorrect (missing package libX11-dev,libxext-dev or libxt-dev) . ],1); fi
      LIBS="$ff_libs"
  
   fi

fi


# Checking that OpenGL is accessible
# ----------------------------------

ff_gl_ok=no
ff_glx_ok=no
ff_agl_ok=no
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
        ff_agl_ok=yes,
        ff_gl_ok=no
        ff_agl_ok=no)
    
    
    if test "$ff_gl_ok" = no 
    then
	AC_CHECK_HEADERS(GL/gl.h,
	    ff_gl_ok=yes
	    ff_agl_ok=yes,
	    ff_gl_ok=no
	    ff_agl_ok=no)
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
    
    if test "$ff_gl_ok" = yes
    then
	AC_CHECK_HEADERS(GL/glx.h,
	    ff_glx_ok=yes,
	    ff_glx_ok=no)
    fi
    if test "$ff_glx_ok" = yes -a   "$ff_x11" = yes 
    then
	AC_CHECK_LIB(GL,
	    glBegin,
	    ff_glx_ok=yes,
	    ff_glx_ok=no)
    fi
    if test "$ff_glx_ok" = yes
    then
	AC_CHECK_LIB(GLU,
	    gluSphere,
	    ff_glx_ok=yes,
	    ff_glx_ok=no,
	    -lGL)
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

if test "$ff_glut_ok" = yes
then
    AC_SUBST(LIBSGLUT,$ff_glut)
    if test -n  "$ff_pthread" ; then 
	ff_ffglutprog="ffglut${EXEEXT}"
	AC_DEFINE_UNQUOTED(PROG_FFGLUT,"$ff_ffglutprog", the ffglut application for the new graphics )
	AC_SUBST(LIBSPTHREAD,$ff_pthread)
    fi
    ff_meditprog="ffmedit${EXEEXT}"
fi
# echo "  .... LIBS = $LIBS .... " 

# If we cannot find GL, we cannot compile FreeFem++-agl
if test "$ff_agl_ok" = no
then

	ff_aglprog=
fi

# If we cannot find GLX, we cannot compile FreeFem++ with it
if test "$ff_glx_ok" = yes
then
	ff_glxprog=FreeFem++-glx${EXEEXT}
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
dnl m4_include(ax_mpi.m4)
ff_save_path="$PATH"
# We need to choose between mpich, openmpi  and lam for the Debian package
AC_ARG_WITH(mpi,[  --with-mpi=[yes|no|mpic++|lam|mpich|openmpi|/usr/local/bin/mpic++|... ]	or --without-mpi	Choose MPI implementation (default is mpic++)])

#if test "$with_mpi" != no ; then
#AX_MPI(with_mpi=yes, with_mpi=no)
#fi

# Default is mpic++ 
ff_mpi_suffix="";
if test "$with_mpi" = yes -o -z "$with_mpi" 
then
   ff_mpicxx=mpic++
else 
  case "$with_mpi" in
 lam|mpich|openmpi)   ff_mpi_suffix=.$with_mpi;ff_mpicxx=mpic++.$with_mpi;;
 *)  ff_mpicxx="$with_mpi" ;;
 esac
fi

dnl  correct ff_mpi_path august 2010 -- FH ...  
if test "$with_mpi" != no
then
	ff_mpi_path=`AS_DIRNAME(["$ff_mpicxx"])`
dnl	echo "ff_mpi_path '$ff_mpi_path' .............."
	case "$ff_mpi_path" in
	    .|"") ff_mpi_path="$PATH";ff_defmpicxx="$ff_mpicxx";;
	    *) ff_mpi_path="$ff_mpi_path";ff_defmpicxx=`expr "//$ff_mpicxx" : '.*/\(.*\)'`;; 
dnl if also  add $PATH they  could be missing some different mpi version... 
	esac	
        
	AC_ARG_VAR(MPICXX,[MPI C++ compiler command])
	AC_PATH_PROGS(MPICXX,$ff_defmpicxx mpic++$ff_mpi_suffix mpicxx$ff_mpi_suffix mpiCC$ff_mpi_suffix mpCC hcp mpxlC mpxlC_r cmpic++,no,$ff_mpi_path)
	AC_MSG_CHECKING(for $MPICXX)
	ff_mpicxx="$MPICXX"
	ff_save_cxx="$CXX"
	CXX=$ff_mpicxx
	AC_LINK_IFELSE(
[
#include <mpi.h>
#include <stdio.h>
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
	AC_MSG_RESULT($ff_mpi)

	# Also check that mpirun is there. If it isn't, then MPI is
	# not fully installed.

	if test "$ff_mpi" = yes;
	then
	echo "ff_mpi_path::$ff_mpi_path:"
		AC_PATH_PROGS(ff_mpirun,mpirun mpiexec,no,$ff_mpi_path)
		if test "$ff_mpirun" = no
		then
			ff_mpi=no
		fi
	fi

	if test "$ff_mpi" = yes;
	then
	  echo "MPI CC $ff_mpi" >config_LIB_INFO

		# We do not AC_DEFINE any special flag for parallel
		# computation here, because it must only be set when the
 		# parallel program is compiled (see src/mpi/Makfile.am)
		ff_mpiprog="FreeFem++-mpi${EXEEXT}"
   		  AC_SUBST(MPIPROG,$ff_mpiprog)
   		  AC_SUBST(MPIRUN,$ff_mpirun)
                  AC_SUBST(MPICXX,$ff_mpicxx)
	else
	        AC_SUBST(MPICXX,$ff_save_cxx)
	fi
	if test "$ff_mpi" = yes;
	then
          ff_mpishow=`$MPICXX -show`
#          case $ff_mpishow in
#	  *openmpi*) echo which mpi : openmpi;;
#	  *lam*) echo lam mpi;;
#	  *mpich2*) echo mpich2 ;;
#	  *) echo unknow mpi;;
#	  esac


	    [ff_MPI_INCLUDE=`$MPICXX -show|tr ' ' '\n'| grep -E '^-[^WLlOgp]|^-Wp,'|tr '\n' ' '`]
#	    ff_MPI_INCLUDE=
	    ff_MPI_LIB_DIRS=""
	    [ff_MPI_LIB=`$MPICXX -show|tr ' ' '\n'| grep -E '^-[Llp]|^-Wl,'|tr '\n' ' '`]
	    [ff_mpi_idir=`$MPICXX -show|tr ' ' '\n'| grep -E '^-I'|sed s/^-I//|tr '\n' ' '`' /usr/include']
            ff_MPI_INCLUDE_DIR=""
#	    echo ff_MPI_INCLUDE $ff_MPI_INCLUDE
#	    echo ff_MPI_LIB $ff_MPI_LIB
	    for i in $ff_mpi_idir; do
	    if test -f "$i/mpi.h" -a -z "$ff_MPI_INCLUDE_DIR"  ;then
	      ff_MPI_INCLUDE_DIR=$i
	    fi
	    done

	    AC_SUBST(MPI_INCLUDE,$ff_MPI_INCLUDE)
	    AC_SUBST(MPI_LIB_DIRS,$ff_MPI_LIB_DIRS)
	    AC_SUBST(MPI_LIB,$ff_MPI_LIB)

	    AC_SUBST(MPICXX,$MPICXX)		
	    AC_ARG_VAR(MPICC,[MPI C compiler command in $ff_mpi_path])
	    AC_PATH_PROGS(MPICC,mpicc$ff_mpi_suffix hcc mpcc mpcc_r mpxlc cmpicc, "",$ff_mpi_path)
	    AC_SUBST(MPICC,$MPICC)
	    AC_SUBST(PASTIX_HOSTARCH,$ff_HOSTARCH_pastix)
	    if test "$enable_fortran" != no
	    then	
		
		AC_ARG_VAR(MPIF77,[MPI Fortran 77 compiler command])
		AC_PATH_PROGS(MPIF77, mpif90$ff_mpi_suffix mpif77$ff_mpi_suffix hf77 mpxlf mpf77 mpif90 mpf90 mpxlf90 mpxlf95 mpxlf_r cmpifc cmpif90c, "",$ff_mpi_path)
		AC_SUBST(MPIF77)
		
		AC_ARG_VAR(MPIFC,[MPI Fortran 90  compiler command])
		AC_PATH_PROGS(MPIFC, mpif90$ff_mpi_suffix mpxlf95_r mpxlf90_r mpxlf95 mpxlf90 mpf90 cmpif90c, "",$ff_mpi_path)
		AC_SUBST(MPIFC)
	    fi
	    if test ! -f "$ff_MPI_INCLUDE_DIR/mpif.h"  ; then
	        AC_MSG_NOTICE([ MPI without fortran no file "$ff_MPI_INCLUDE_DIR/mpif.h"  ])
	    else
		if test -n "$MPIFC" ; then
		    [echo mpifc LD `$MPIFC -show|tr ' ' '\n'| grep -E '^-[Llp]|^-Wl,'|tr '\n' ' '`  >>$ff_where_lib_conf]
		    [echo mpifc INCLUDE `$MPIFC -show|tr ' ' '\n'| grep -E '^-[ID]|^-Wp,'|tr '\n' ' '`  >>$ff_where_lib_conf]
		fi
		if test -n "$MPIF77" ; then
		    [echo mpif77 LD `$MPIF77 -show|tr ' ' '\n'| grep -E '^-[Llp]|^-Wl,'|tr '\n' ' '`  >>$ff_where_lib_conf]
		    [echo mpif77 INCLUDE `$MPIF77 -show|tr ' ' '\n'| grep -E '^-[ID]|^-Wp,'|tr '\n' ' '`  >>$ff_where_lib_conf]
		fi
  	    fi
	    if test -n "$MPICXX" ; then 	    
#	        [echo " ********* MPICXX $MPICXX " mpi LD `$MPICXX -show|tr ' ' '\n'| grep -E '^-[Llp]|^-Wl,'|tr '\n' ' '` ]
		[echo mpi LD `$MPICXX -show|tr ' ' '\n'| grep -E '^-[Llp]|^-Wl,'|tr '\n' ' '`  >>$ff_where_lib_conf]
		[echo mpi INCLUDE `$MPICXX -show|tr ' ' '\n'| grep -E '^-[ID]|^-Wp,'|tr '\n' ' '`  >>$ff_where_lib_conf]
	    fi
	    AC_SUBST(MPI_INC_DIR,$ff_MPI_INCLUDE_DIR)      		
	fi
	CXX="$ff_save_cxx"
fi
#export PATH="$ff_save_path"
# Looking for useful configuration utilities
# ------------------------------------------

#getting utility software archives through the internet
AC_ARG_ENABLE(generic_blas,[  --enable-generic_blas	Download generic blas (implie --enable-download)])

AC_ARG_ENABLE(download,[  --enable-download	Download missing libraries (BLAS,ARPACK,UMFPACK,FLTK)])

#replacing wget with another command
ff_wget_command="wget --no-check-certificate"
ff_wget=no

AC_ARG_WITH(wget,
	[  --with-wget=command	Replace "wget" with another command. Implies --enable-download],
	   enable_download=yes
	   ff_wget=yes
	   ff_wget_command="${withval}")

if test "$enable_generic_blas" = yes
then
	enable_download=yes
fi
#  echo " ----------- $enable_download --"
# if enabling downloads find wget or curl to do download 
# ------------------------------------------------------

if test  "$ff_wget" =  no 
then
	ff_wget=""	
	AC_CHECK_PROG(ff_wget,wget --no-check-certificate,yes,no)
fi	
if test  "$ff_wget" = no 
then
        ff_wget_command="curl -O"
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
DOWNLOADCOMPILE=compile-pkg
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

dnl  : Looking for the FFTW3
dnl --------------------------
AC_CHECK_LIB(fftw3,fftw_execute,ff_fftw_ok=yes)
AC_CHECK_HEADERS(fftw3.h,ff_fftw_h=yes, ff_fftw_h=no)
if test "$ff_fftw_ok" = yes -a "$ff_fftw_h" = yes
then
    [echo fftw3 LD  -lffw3   >>$ff_where_lib_conf]
fi
dnl ---------------------------
dnl   Looking for the tetgen
dnl --------------------------
AC_CHECK_LIB(tet,fftw_execute,ff_fftet_ok=yes)
AC_CHECK_HEADERS(tetgen.h,ff_tet_h=yes, ff_tet_h=no)
if test "$ff_tet_ok" = yes -a "$ff_tet_h" = yes
then
    [echo tetgen LD  -ltet   >>$ff_where_lib_conf]
fi
dnl ---------------------------


# Modif: morice 2008
# Looking for the tetgen
# --------------------
ff_tetgen_ok=no
ff_tetgen_dir="${curdir}/download"
ff_libs_old="$LIBS"
ff_ldflags_old="$LDFLAGS"
ff_cppflags_old="$CPPFLAGS"
#for ff_tetgendir in "" "${curdir}/download";
#do
if test "$ff_tetgen_ok" = no
	then
	  if test -n "$ff_tetgen_dir"
	  then
           LDFLAGS="$LDFLAGS -L${ff_tetgen_dir}/lib"
           CPPFLAGS="$CPPFLAGS -I${ff_tetgen_dir}/include"
      fi

  AC_CHECK_LIB(tet,tetrahedralize,
                ff_tetgen_ok=yes)
  AC_CHECK_HEADERS(tetgen.h,
             ff_tetgen_h=yes, ff_tetgen_h=no)

   if test "$ff_tetgen_ok" = yes  -a "$ff_tetgen_h" = yes 
    then
      AC_DEFINE(HAVE_TETGEN,1, tetgen is compute tetrahedralize volume of an enclosed surface)
   else
       LIBS="$ff_libs_old"
       LDFLAGS="$ff_ldflags_old"
       CPPFLAGS="$ff_cppflags_old"    
   fi
fi
#done
# Looking for the BLAS
# --------------------

# User-specified location

ff_blas_ok=no
AC_ARG_WITH(blas,
	[  --with-blas=library	Use a specific version of the Blas],
	ff_blas_ok=yes
	ff_blas_libs="${withval}")

ff_blas_inc=""

# Specific BLAS library location for FreeBSD
ff_freebsd_blas="-lf2c -lf77blas -latlas -lgslcblas"
#  zmach is in lapack 
# Trying  blas library
#echo "LIBS = " $LIBS
for iblas in \
"-framework vecLib" "-lblas"  "-L/usr/lib/atlas -lblas" \
$ff_freebsd_blas
do
#   
#"-L${curdir}/download/blas -lcblas -lf77blas" \
#"-L${curdir}/download/blas/ATLAS/lib/ff++ -lcblas -lf77blas -latlas" \ 
#"-L${curdir}/download/blas/WinNT_PPRO256 -lcblas -lf77blas -latlas" \
#
if test "$ff_blas_ok" = no;
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
if test "$ff_blas_ok" = yes 
then
    echo blas LD  \'$ff_blas_libs\'  >>$ff_where_lib_conf
fi
#echo "LIBS = " $LIBS
# If all else fails, download!

if test "$ff_blas_ok" = no -a "$enable_download" = yes \
   -a "$enable_fortran" != no -a "$enable_c" != no
then
	AC_CHECK_PROG(ff_unzip,unzip,yes,no)

	AC_MSG_CHECKING(for BLAS version to download)

	# Do not update $LIBS, but create an extra LIB variable,
	# because this lib does not exist yet, and this could make the
	# following tests fail.

	# When compiling a generic version, we do not need an
	# optimized version of the BLAS.

#	if test "$enable_generic" = yes -o "$enable_generic_blas" = yes 
#   supress the automatic compilation of atlas, to expencive  (FH version 3.0-5)
	if true 
	then
		ff_download_blas=generic
		ff_blas_libs="-L${curdir}/download/blas -lcblas -lf77blas"
		ff_blas_inc="-I${curdir}/download/blas/CBLAS/src"

		AC_SUBST(GENERIC_BLAS,"libf77blas.a libcblas.a")
		AC_SUBST(GENERIC_BLAS_BUILT_SOURCES,'BLAS CBLAS  $(F77BLAS_SOURCES) $(CBLAS_SOURCES)')
	else
		# There are precompiled versions for Windows. We cannot
		# distinguish the different versions easily so we just
		# download the most generic one. A more specific version can
		# always be downloaded by hand.

		# BUG: at the moment, some necessary symbols are missing from
		# the precompiled version (e.g. _dgemm_).

###	if test "$ff_win32" = yes -a "$enable_generic" = yes -a "$ff_unzip" = yes
		if false
		then
			ff_download_blas=atlas-binary
			ff_blas_libs="-L${curdir}/download/blas/WinNT_PPRO256 -lcblas -llapack -lf77blas -latlas"
			ff_blas_inc="-I${curdir}/download/blas/WinNT_PPRO256"
			ff_lapack_ok=ok
		else
			ff_download_blas=atlas-source
			ff_blas_libs="-L${curdir}/download/blas/ATLAS/lib/ff++ -lcblas -llapack -lf77blas -latlas"
			ff_blas_inc="-I${curdir}/download/blas/ATLAS/include"
			ff_lapack_ok=ok
		fi
	fi

	AC_MSG_RESULT($ff_download_blas)
	AC_SUBST(DOWNLOAD_BLAS,$ff_download_blas)
	ff_blas_ok=yes
        ff_cblas_h=yes
	ff_cblas_libs=yes
fi
AC_MSG_RESULT($ff_blas_inc)
AC_SUBST(BLASLIBS,$ff_blas_libs)
AC_SUBST(BLASINC,$ff_blas_inc)
#looking for cblas.h  FH

if test "$ff_blas_ok" = yes -a "$ff_cblas_libs" = yes ;
then
  AC_CHECK_HEADERS(cblas.h,
       ff_cblas_h=yes,
       ff_cblas_h=no)
  AC_CHECK_HEADERS(vecLib/cblas.h,
       ff_cblas_h=yes
       ff_cblas_h=no)
  AC_CHECK_HEADERS(atlas/cblas.h,
       ff_cblas_h=yes
       ff_cblas_h=no)

fi



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
            AC_LINK_IFELSE( [AC_LANG_CALL(,dgeqr2_)], ff_lapack_ok=yes)	
	    if test "$ff_lapack_ok" = no ; then
	    LIBS="$ff_save_libs -llapack $ff_blas_libs"
            AC_LINK_IFELSE( [AC_LANG_CALL(,dgeqr2_)], [ff_lapack_ok=yes
	    ff_lapack_lib="-llapack"] )	
	    fi		    
	    LIBS="$ff_save_libs"	   
	fi
    fi

AC_MSG_RESULT($ff_lapack_ok)

if test "$ff_lapack_ok" = no; then
AC_ARG_WITH(lapack,
	[  --with-lapack=library	Use a specific version of the Lapack],
	ff_lapack_ok=yes
	ff_lapack_lib="${withval}"
	LIBS="$ff_lapack_lib $LIBS"
)
fi

ff_lapackdir='$$(LAPACKdir)'
if  test "$ff_lapack_ok" = yes ; then
  # no compilation of lapack in  arpack 
  ff_lapackdir=
  if test "$ff_lapack_ok" = yes   
  then
      echo lapack LD  \'$ff_lapack_lib\'  >>$ff_where_lib_conf
  fi   

else
  ff_lapack_lib=-llapack
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
      echo  arpack LD \'$ff_arpack_libs\'  >>$ff_where_lib_conf
  fi   
  
	# If all else fails, download!
	if test "$ff_arpack_ok" = no -a "$enable_download" = yes \
	   -a "$enable_fortran" != no
	then
                ff_arpack_download=yes
		AC_MSG_NOTICE(using downloaded Arpack)
		AC_SUBST(DOWNLOAD_ARPACK,arpack)
		AC_SUBST(FF_LAPACKdir,$ff_lapackdir)
		AC_SUBST(ARPACKLIB,${curdir}/download/lib/libarpack.a)
		AC_SUBST(LAPACKLIB,${curdir}/download/lib/liblapack.a)
		# Do not update $LIBS, but create an extra LIB variable,
		# because this lib does not exist yet, and this could make the
		# following tests fail.
                ff_arpack_libs="-L${curdir}/download/lib -larpack  $ff_lapack_lib"
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
	AC_SUBST([EIGENTEST],../regtests.sh)
fi
AC_SUBST([EIGENOBJ])

# Looking for UMFPACK
# -------------------

ff_amd_ok=no
ff_umfpack_ok=no
ff_save_libs="$LIBS"
if test "$ff_blas_ok" = yes;
then
        # User-specified location
        AC_ARG_WITH(amd,
                [  --with-amd=library	Use a specific version of AMD],
                ff_amd_ok=yes
		ff_umfpack_libs="${withval}")

        AC_ARG_WITH(umfpack,
                [  --with-umfpack=library	Use a specific version of Umfpack],
                ff_umfpack_ok=yes
                ff_umfpack_libs="${withval} $ff_umfpack_libs"
               )
#		AC_DEFINE(HAVE_LIBUMFPACK,1,
#			Umfpack is used for sparse matrices computations))

	# Trying default locations

	AC_CHECK_HEADERS(umfpack.h umfpack/umfpack.h ufsparse/umfpack.h suitesparse/umfpack.h,
	    ff_umfpack_header=yes
	    ff_umfpack_dir=`dirname $ac_header`)

	# Somes systems like FreeBSD hide umfpack.h in a directory
	# called UMFPACK (all capitals). This breaks the standard
	# #define produced by autoconf in config.h.in.

	if test "$ff_umfpack_header" != yes;
	then
		AC_CHECK_HEADER(UMFPACK/umfpack.h,
			ff_umfpack_header=yes
			ff_umfpack_dir=UMFPACK
			AC_DEFINE(HAVE_BIG_UMFPACK_UMFPACK_H,1,
				If umfpack.h is located in UMFPACK subdir))
	fi
	LIBS="$ff_blas_libs $LIBS"
	if test "$ff_amd_ok" = no;
	then
		AC_CHECK_LIB(amd,amd_info,
			ff_umfpack_libs="$ff_umfpack_libs -lamd"
			ff_amd_ok=yes)
	fi

	if test "$ff_umfpack_ok" = no -a "$ff_amd_ok" = yes;
	then        
		AC_CHECK_LIB(umfpack,umf_i_malloc,
			ff_umfpack_libs=" -lumfpack $ff_umfpack_libs"
			ff_umfpack_ok=yes,,$ff_umfpack_libs)
		if test "$ff_umfpack_ok" != yes -o "$ff_umfpack_header" != yes ; then
			AC_MSG_WARN([		Sorry we find UMFPACK lib but not the include umfpack.h file])
		fi
	fi	
        if  test	 "$ff_umfpack_header" != yes -a "$ff_umfpack_ok" != no; then
            ff_umfpack_ok=no
         fi
	 if test "$ff_umfpack_ok" = yes -a "$ff_amd_ok" = yes; then
		AC_DEFINE(HAVE_LIBUMFPACK,1, Umfpack is used for sparse matrices computations )			
	 fi
fi
  if test "$ff_umfpack_ok" = yes 
  then
      echo  amd LD \'$ff_umfpack_libs\'  >>$ff_where_lib_conf
      echo  amd INCLUDE  \'-I/usr/include/$ff_umfpack_dir\'  >>$ff_where_lib_conf
      echo  umfpack LD  '$ff_umfpack_libs\'   >>$ff_where_lib_conf     
      echo  umfpack  INCLUDE  \'-I/usr/include/$ff_umfpack_dir\'  >>$ff_where_lib_conf
  fi   

LIBS="$ff_save_libs"


# If all else fails, download!
if test "$ff_umfpack_ok" = no -a "$enable_download" = yes 
   then
	AC_MSG_NOTICE(using downloaded UMFPACK)
	AC_SUBST(DOWNLOAD_UMFPACK,umfpack)
	ff_umfpack_download=yes
	# Do not update $LIBS, but create an extra LIB variable,
	# because this lib does not exist yet, and this could make the
	# following tests fail.

	ff_umfpack_libs="${curdir}/download/lib/libumfpack.a ${curdir}/download/lib/libamd.a"
	AC_DEFINE(HAVE_LIBUMFPACK,1,UMFPACK)
	if test "$ff_win32" = yes; then
  	AC_SUBST(FF_UMFPACK_CONFIG,-DCBLAS)
	fi
	ff_umfpack_ok=yes
fi

if test "$ff_umfpack_ok" = no
then
	 AC_MSG_NOTICE( -- NOT  UMFPACK  ff_wget = $ff_wget)
fi

AC_SUBST(UMFPACKLIBS,$ff_umfpack_libs)

# Checking for some functions that may not appear everywhere
# ----------------------------------------------------------

# asinh acosh atanh are not in Mingw yet
ff_malloc_h=""
AC_CHECK_HEADERS(malloc.h,ff_malloc_h=1)
AC_SUBST(FF_MALLOC_H,$ff_malloc_h)
AC_CHECK_FUNCS(asinh acosh atanh getenv jn erfc tgamma)
AC_CHECK_FUNCS(srandomdev)

AC_CHECK_FUNCS(second_,ff_second="",ff_second=second.o)
AC_SUBST(FF_SECOND,"$ff_second")

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
	AC_COMPILE_IFELSE([#include <dlfcn.h>],
		ff_dynload=yes,
		ff_dynload=no)
	AC_MSG_RESULT($ff_dynload)
fi

# Checks that we also have the corresponding library
if test "$ff_dynload" = yes
then
	AC_CHECK_LIB(dl,dlinfo)

	# Checks that everythings works ok
	AC_MSG_CHECKING(whether dlopen links ok)
	AC_LINK_IFELSE(
[#include <dlfcn.h>
int main(int argc,char **argv){
  dlopen("",RTLD_LAZY);
  return 0;
}],
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
		CHECK_COMPILE_FLAG(C++,-rdynamic,CXXFLAGS)
		CHECK_COMPILE_FLAG(C,-rdynamic,CNOFLAGS)
		CHECK_COMPILE_FLAG(C,-rdynamic,CFLAGS)
	fi
	CHECK_COMPILE_FLAG(C++,-fPIC,CXXFLAGS)
	CHECK_COMPILE_FLAG(C,-fPIC,CFLAGS)
	CHECK_COMPILE_FLAG(C,-fPIC,CNOFLAGS)
	if test "$enable_fortran" != no 
	 then 
	 CHECK_COMPILE_FLAG(Fortran,-fPIC,FFLAGS)	
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

# PDF documentation building sometimes poses problems because of
# pdfsync.sty. So we need to be able to disable it.
enable_pdf=yes
AC_ARG_ENABLE(pdf,[  --disable-pdf	Disable PDF documentation building])
if test "$enable_pdf" = yes 
then
   AC_CHECK_PROG(ff_pdflatex,pdflatex,yes,no)
   if test "$ff_pdflatex" = yes -a $EPSTOPDF != false -a $ff_convert = yes;
   then
	AC_SUBST(DOCPDF,"freefem++doc.pdf")
   fi
fi

# HISTORY logging through CVS, but only if we are in a CVS working area
if test -d CVS
then
	AC_CHECK_PROG(ff_history,cvs2cl,yes,no)
	if test "$ff_history" = yes
	then
	   AC_SUBST(HISTORY,history)
	fi
fi

# Choosing compilation options for the standard version (in src/std)
# ------------------------------------------------------------------

# The "standard" configured version can use win32 (mingw) or x11

if test "$ff_mingw" = yes 
then
	ff_stdprog="FreeFem++-std${EXEEXT}"
	ff_std_graph_obj=Pcrgraph.$OBJEXT
	ff_std_ldflags=-mwindows
	ff_std_libs=
elif test "$ff_x11" = yes ; then
	ff_std_graph_obj=Xrgraph.$OBJEXT
	ff_std_ldflags=
	ff_std_libs=$ff_x11_libs
	# If the standard version already uses x11, we do not need to
	# make the contents of the src/x11 directory.
	ff_x11prog="FreeFem++-x11${EXEEXT}"
fi
AC_SUBST(STD_GRAPH_OBJ,$ff_std_graph_obj)
AC_SUBST(STD_LDFLAGS,$ff_std_ldflags)

AC_SUBST(STD_LIBS,$ff_std_libs)

# Can we produce an IDE?
# ----------------------

# For development purposes, we may need to use X libraries instead of
# native libraries (mainly to run the programs remotely)
#  AC_ARG_WITH(fltk,[  --without-fltk	(to suppress explicity    --enable-download imply --with-fltk ) ])
#  AC_ARG_WITH(xfltk,[  --with-xfltk	When recompiling FLTK, use X libraries instead of native Win32 or MacOS])
#  remove all fltk stuff
# to remove of FLTK
#AC_MSG_NOTICE(without ide/ fltk ..... Please wait the next version.... ) 
if test 0 != 0 
 then  
if test "$with_fltk" = no 
then 
     ff_fltk=no
	AC_MSG_NOTICE(without FLTK sorry in next no ide today ) 
else
# under mingw use download fltk version
if test "$ff_mingw" = yes 
then
	AC_MSG_NOTICE(under mingw not using default FLTK)
 enable_default_fltk=no 
fi

# Default FLTK (in the machine's default PATH)
if test "$enable_default_fltk" = no
then
	AC_MSG_NOTICE(not using default FLTK)
	ff_fltk=no
else

	# Checks that fltk-config exists
	ff_fltk_config=fltk-config
	AC_CHECK_PROG(ff_fltk,fltk-config,yes,no)

	# Checks that FLTK is compiled with threads enabled
	AC_MSG_CHECKING(that FLTK is compiled with threads)
	AC_COMPILE_IFELSE(
[#include <FL/Fl.H>
main(){
  Fl::lock();
}
],
		ff_fltk_threads=yes,
		ff_fltk_threads=no)
	AC_MSG_RESULT($ff_fltk_threads)
	if test "$ff_fltk_threads" = no
	then
		ff_fltk=no
		AC_MSG_NOTICE(not using default FLTK because it does not know threads)
	fi	
fi

# Downloading and compiling FLTK
if test "$ff_fltk" = no -a "$enable_download" = yes
then
	AC_MSG_NOTICE(using downloaded FLTK)
	AC_SUBST(DOWNLOAD_FLTK,fltk)
	ff_fltk_version=1.1.6
	AC_SUBST(FLTK_VERSION,$ff_fltk_version)
	ff_fltk_config="$curdir/download//bin/fltk-config"
	ff_fltk=yes
fi

if test "$ff_fltk" = yes
then
	# For the moment, FreeFem++-ide${EXEEXT} is not built
	# because it first needs FreeFEM++ to accept several
	# successive runs without restart.
        
	ff_ideprog="FreeFem++-ide${EXEEXT} FreeFem++-cs${EXEEXT} FreeFem++-client${EXEEXT} FreeFem++-server${EXEEXT}"
	ff_idetest="./testhighlight.sh"
	ff_testideprog="testhighlight${EXEEXT}"
	if test "$enable_cygwindll" = yes 
	 then
           AC_SUBST(FLTK_CONFIG_PARAM,"--enable-localzlib --enable-cygwin")
        elif test "$ff_mingw" = yes 
	 then
           AC_SUBST(FLTK_CONFIG_PARAM,"--enable-localzlib --disable-cygwin")
	else
           AC_SUBST(FLTK_CONFIG_PARAM,"--disable-cygwin")
	fi

	AC_SUBST(FLTK_CONFIG,$ff_fltk_config)

	# Only links X11 libraries into FreeFem++-cs if we are not on
	# Windows (Mingw) or Macintosh (these do not need X11)

	if test "$ff_mac" != yes -a "$ff_mingw" != yes 
	then
		AC_SUBST(IDE_LIBS,$ff_x11_libs)
	fi

	if test "$with_xfltk" = yes
	then
		AC_SUBST(FLTK_XONLY,"--with-x")
	fi
fi
fi
# end of remove fltk
fi 
# Flex compatibility with Bison
# -----------------------------

# If Flex does not know Bison options, it cannot deal with
# src/ide/hl_lex.l++. So we need to make src/ide/hl_lex.c++ up-to-date
# by hand (under Cygwin, "cvs update" does not make hl_lex.c++ newer
# than hl_lex.l++ even if it was when a commit was done).
AC_MSG_CHECKING(flex compatibility with bison)
if test `flex --help|fgrep bison|wc -l` = 0
then
	AC_MSG_RESULT(no)
	touch src/ide/hl_lex.c++
else
	AC_MSG_RESULT(yes)
fi

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
ff_bamgprog="bamg${EXEEXT}  cvmsh2${EXEEXT}"
if test "$ff_x11" = yes -o  "$ff_win32" = yes   ; then
ff_bamgprog="$ff_bamgprog drawbdmesh${EXEEXT}"
fi

AC_SUBST(CNOFLAGS,$CNOFLAGS) dnl for superludist CFLAGS without optim  ...
# The final list of executable programs
AC_SUBST(MEDITPROG,$ff_meditprog)
AC_SUBST(FFGLUTPROG,$ff_ffglutprog)
AC_SUBST(BAMGPROG,$ff_bamgprog)
AC_SUBST(GLXPROG,$ff_glxprog)
AC_SUBST(STDPROG,$ff_stdprog)
AC_SUBST(AGLPROG,$ff_aglprog)
AC_SUBST(X11PROG,$ff_x11prog)
AC_SUBST(IDEPROG,$ff_ideprog)
AC_SUBST(TESTIDEPROG,$ff_testideprog)
# The final test program
AC_SUBST(IDE_TESTS,$ff_idetest)
ff_progs="$ff_stdprog FreeFem++-nw $ff_bamgprog $ff_glxprog $ff_aglprog $ff_x11prog $ff_ideprog  $ff_mpiprog $ff_meditprog $ff_ffglutprog"


#echo xxxxxxxxxxx
# hips LD -L/Users/hecht/work/freefem++/download/lib/hips -lio -lhips
# hips INCLUDE -I/Users/hecht/work/freefem++/download/include/hips
# hipssequential LD -L/Users/hecht/work/freefem++/download/lib/hips -lio -lhipssequential
# hipssequential INCLUDE -I/Users/hecht/work/freefem++/download/include/hips
# hypre LD -L/Users/hecht/work/freefem++/download/lib/hypre  -lHYPRE
# hypre INCLUDE -I/Users/hecht/work/freefem++/download/include/hypre
#metis LD -L/Users/hecht/work/freefem++/download/lib -lmetis
#metis INCLUDE -I/Users/hecht/work/freefem++/download/include/metis
# mumps LD '-L/Users/hecht/work/freefem++/download//lib  -lpord -lmumps_common -ldmumps -lzmumps -lpord'
# mumps INCLUDE '-I/Users/hecht/work/freefem++/download/include'
# parmetis LD -L/Users/hecht/work/freefem++/download/lib -lparmetis -lmetis
# parmetis INCLUDE -I/Users/hecht/work/freefem++/download/include/metis
# scalapack LD '-L/Users/hecht/work/freefem++/download/lib  -lscalapack '
# scalapack INCLUDE '-I/Users/hecht/work/freefem++/download/include'
# ptscotch LD '-L/Users/hecht/work/freefem++/download/lib  -lptesmumps -lptscotch -lptscotcherr '
# ptscotch INCLUDE '-I/Users/hecht/work/freefem++/download/include/scotch'
# scotch LD '-L/Users/hecht/work/freefem++/download/lib   -lscotch -lscotcherr '
# scotch INCLUDE '-I/Users/hecht/work/freefem++/download/include/scotch'
# superlu LD -L/Users/hecht/work/freefem++/download/lib -lsuperlu
# superlu INCLUDE -I/Users/hecht/work/freefem++/download/include
# superlu_dist LD '-L/Users/hecht/work/freefem++/download/lib -lsuperlu_dist_2.3'
# superlu_dist INCLUDE '-I/Users/hecht/work/freefem++/download/include/superludist'
# tetgen LD -L/Users/hecht/work/freefem++/download/lib -ltet
# tetgen INCLUDE -I/Users/hecht/work/freefem++/download/include/
# arpack LD '-L/Users/hecht/work/freefem++/download/lib' '-larpack'
# arpack INCLUDE '-I/Users/hecht/work/freefem++/download/include'
# umfpack LD '-L/Users/hecht/work/freefem++/download/lib' '-lumfpack'
# umfpack INCLUDE '-I/Users/hecht/work/freefem++/download/include'
# amd LD '-L/Users/hecht/work/freefem++/download/lib' '-lamd'
# amd INCLUDE '-I/Users/hecht/work/freefem++/download/include'
# fftw3 LD '-L/Users/hecht/work/freefem++/download/lib' '-lfftw3'
# fftw3 INCLUDE '-I/Users/hecht/work/freefem++/download/include'
# mmg3d LD '-L/Users/hecht/work/freefem++/download/lib' '-lmmg3d'
# mmg3d INCLUDE '-I/Users/hecht/work/freefem++/download/include'


if test -f download/pkg/mmg3dlib.tar.gz  ;then
    AC_MSG_NOTICE([      Compile : download/pkg/mmg3dlib.tar.gz  (3d mesh adapaton)])
else
    ff_MACOS=`expr "xx $CXXFLAGS" : '.*-mmacosx-version-min=\([0-9.]*\).*'`
    
    case $ff_MACOS  in
	[0-9.]*) ff_MACOS="$ff_MACOS";;
	*) ff_MACOS=`sw_vers -productVersion| awk -F. '{print $1"."$2}'`
    esac
    ff_mmg3dlib=""
    if test "$ff_wget" != no   ; then  
	if test "$ff_mac" = yes ; then
	    MACOS=10.4
	    ff_mmg3dlib="libmmg3d-1001-OSX-$ff_MACOS.tar.gz"
	elif test "$ff_win32" = yes ; then
	    ff_mmg3dlib="libmmg3d-1001-WIN32.tar.gz"
	else 
	    ff_mmg3dlib="libmmg3d-1001-LINUX-$ff_size_ptr.tar.gz"
	fi
	$ff_wget_command  http://www.freefem.org/ff++/ftp/libmmg3d/$ff_mmg3dlib
	tar zxvf $ff_mmg3dlib -C download 

	AC_MSG_CHECKING(for MMG_mmg3dlib in mmg3dlib )
	    ff_save_libs="$LIBS"
	    LIBS="$LIBS -Ldownload/lib -lmmg3d"
	    AC_LINK_IFELSE(
		[AC_LANG_CALL(,MMG_mmg3dlib)],
		ff_mmg3d_ok=yes 
		mv $ff_mmg3dlib download/pkg/.
		,
		ff_mmg3d_ok=no
		rm download/lib/libmmg3d.a download/include/libmmg3d.a
                ff_mmg3d_ok=
		)
	    LIBS="$ff_save_libs"
	    AC_MSG_RESULT($ff_mmg3d_ok)
    fi
    AC_SUBST(MMG3DLIB_TAR,$ff_mmg3dlib) 
fi

ff_with_mpi=-openmpi; 

dnl search of associad software 
m4_define([AC_FF_WHERELIB],
    [  	AC_MSG_CHECKING(check $1)
	ff_save_libs=$LIBS
	LIBS="$LIBS $2"
	AC_LINK_IFELSE([AC_LANG_PROGRAM([m4_if($3,,,[#include <$3>])],[])], [ff_WHERE=yes],[ff_WHERE=no]) 
	if test "$ff_where" = "yes" ; then
	    echo "$1 LD '$2'"  >$ff_where_lib_conf
	    if test -n "$3" ; then
		echo "$1 INCLUDE  /usr/include/AS_DIRNAME[$3]"  >$ff_where_lib_conf
	    fi 
	fi
	AC_MSG_RESULT($ff_WHERE)
	LIBS=$ff_save_libs
	
    ])dnl end m4_define AC_FF_WHERELIB

m4_map([AC_FF_WHERELIB],[ 
	[[mumps],[-lpord -lmumps_common -ldmumps -lzmumps -lpord],[dmumps_c.h]],
	[[fftw3],[-lfftw3],[fftw3.h]],
	[[superlu],[-lsuperlu],[superlu_defs.h]],
	[[blacs],[ -lblacsCinit$ff_with_mpi	-lblacsF77init$ff_with_mpi	-lblacs$ff_with_mpi],[]],
	[[scalapack],[-lscalapack],[]],
	[[scotch],[-lscotch -lscotcherr],[]],	
	[[metis],[-lmetis],[metis/metis.h]],
	[[parmetis],[-lparmetis],[parmetis/parmetis.h]]
    ]
    )


# All makefiles
AC_OUTPUT(Makefile
    download/Makefile
    download/blas/Makefile
    download/arpack/Makefile
    download/umfpack/Makefile
    src/Makefile
    src/bamglib/Makefile
    src/Graphics/Makefile
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
    src/bamg/Makefile
    src/libMesh/Makefile
    src/medit/Makefile
    src/bin-win32/Makefile
    examples++/Makefile
    examples++-eigen/Makefile
    examples++-tutorial/Makefile
    examples++-mpi/Makefile
    examples++-load/Makefile
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
AC_MSG_NOTICE([		 --  FLTK (ide package) $ff_fltk ])
AC_MSG_NOTICE([		 --  with X11             $ff_x11])
AC_MSG_NOTICE([		 --  with MPI             $ff_mpi])
AC_MSG_NOTICE([    progs: $ff_progs ])
if test  "$ff_umfpack_download" = yes  ;then
    AC_MSG_NOTICE([      use of download UMFPACK see  download/umfpack/SuiteSparse/UMFPACK/README.txt for the License])
fi
if test "$ff_arpack_download" = yes  ; then
    AC_MSG_NOTICE([      use of download ARPACK see  download/arpack/ARPACK/README  no License ])
fi

if test "$enable_download" = yes  ;then
    AC_MSG_NOTICE([      use of download fftw   see  download/fftw/fftw-3.2/COPYRIGHT ])
    AC_MSG_NOTICE([      use of download tetgen  see download/tetgen/tetgen1.4.3/LICENSE ])
    AC_MSG_NOTICE([      use of download freeyams  see download/yams/freeyams.2010.02.22 (suface mesh adaptation)  ])
    AC_MSG_NOTICE([      use of download superlu, metis])
    if test -n "$MPICC" ; then
   	AC_MSG_NOTICE([      try to download: blacs parmetis  scalapack scotch superludist pastix hypre hips library])
    fi
    if test -f download/pkg/MUMPS_4.9.2.tar.gz  ;then
	AC_MSG_NOTICE([      try to compile : MUMPS parallele software ])
    else
	AC_MSG_NOTICE([     missing download/pkg/MUMPS_4.9.2.tar.gz, see http://mumps.enseeiht.fr/ to download by hand ])
    fi
    if test -f download/lib/libmmg3d.a ; then
	echo "mmg3d LD '-L`pwd`/download/lib' '-lmmg3d'" >>$ff_where_lib_conf
	echo "mmg3d INCLUDE '-I`pwd`/download/include'" >>$ff_where_lib_conf
    else 
	AC_MSG_NOTICE([   sorry no download/pkg/mmg3dlib.tar.gz 3d mesh adapaton see Dobrzynski Cecile <cecile.dobrzynski@math.u-bordeaux1.fr> ])
    fi
    
    if test "$ff_download_blas" = "generic"  ; then
        AC_MSG_NOTICE([      use of download generic blas and cblas freefem may be slow  ])
        AC_MSG_NOTICE([      you can try to use the  Kazushige Goto s BLAS at http://www.cs.utexas.edu/users/flame/goto/  ])
        AC_MSG_NOTICE([          or at  http://www.tacc.utexas.edu/~kgoto/ for the best BLAS .])
        AC_MSG_NOTICE([      or try to download/compile the altas blas .])
    fi

    if test -f download/pkg/mmg3dlib.tar.gz  ; then
	AC_MSG_NOTICE([      Compile : download/pkg/mmg3dlib.tar.gz  (3d mesh adapaton)])
    elif test 	  "$ff_mmg3d_ok" = yes ; then
	AC_MSG_NOTICE([      Use of precompile version : download/pkg/$ff_mmg3dlib  (3d mesh adapaton)])
    else
	AC_MSG_NOTICE([   sorry no download/pkg/mmg3dlib.tar.gz 3d mesh adapaton see Dobrzynski Cecile <cecile.dobrzynski@math.u-bordeaux1.fr> ])
    fi	
fi  # ---  test "$enable_download" = yes  ;then


if test "$ff_glut_ok" != yes  ; then
    AC_MSG_NOTICE([     *********************************************** ])
    AC_MSG_NOTICE([     WARNING: you do not have the new grachics tools ])
    AC_MSG_NOTICE([         because the configure do not find OpenGL, GLUT or pthread  developer stuff ])
    AC_MSG_NOTICE([         read the README to find missing package  ])
    AC_MSG_NOTICE([         F. Hecht  ])
    AC_MSG_NOTICE([     *********************************************** ])
fi

