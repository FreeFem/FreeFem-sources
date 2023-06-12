# Checking wether we can produce a parallel version
# -------------------------------------------------

dnl m4_include(ax_mpi.m4)
ff_save_path="$PATH"
# We need to choose between mpich, openmpi  and lam for the Debian package
AC_ARG_WITH(mpipath,[  --with-mpipath= the path of mpich under windows (no command  mpic++, ... )])
AC_ARG_WITH(mpilibs,[  --with-mpilibs= the libs to add to c++,fc, ... (to link with c++ - ex:   -L/usr/local/lib -lmpi_f90  -lmpi_cxx -lmpi -lopen-rte -lopen-pal -lutil) ])
AC_ARG_WITH(mpilibsc,[  --with-mpilibsc= the libs to add to c  ... (to link with cc (for pastix lib)  ex:   -L/usr/local/lib -lmpi -lopen-rte -lopen-pal -lutil) ])
AC_ARG_WITH(mpiinc,[  --with-mpiinc= the include directory directive and preprocess directive  (no mpicc++, just use the compiler)) ])
AC_ARG_WITH(mpi,[  --with-mpi=[yes|no|mpic++|lam|mpich|openmpi|/usr/local/bin/mpic++|... ]	or --without-mpi	Choose MPI implementation (default is mpic++)])
if test "$with_mpi" != no ; then
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

dnl AC_MSG_NOTICE([ xxxxxxxxxxxxxxxxxxxx --$with_mpilibs--]);
if test -n "$with_mpiinc"  -a "$with_mpiinc" != no ; then
  if test  "$with_mpi" = 'no' ; then with_mpi='yes'; fi
  ff_MPI_INCLUDE="$with_mpiinc"
fi
if test -n "$with_mpilibs" -a "$with_mpilibs" != no ; then
    ff_MPI_LIB="$with_mpilibs"
    ff_MPI_LIBC="$with_mpilibs"
    ff_MPI_LIBFC="$with_mpilibs"
    MPICXX="$CXX $ff_MPI_INCLUDE"
    MPIFC="$FC  $ff_MPI_INCLUDE"
    MPICC="$CC  $ff_MPI_INCLUDE"
    AC_MSG_NOTICE([   ---  set  all MPI compile to compiler:   $MPICC, $MPIFC, $MPICC ])
fi

if test -n "$with_mpilibsc" -a "$with_mpilibsc" != no ; then
 ff_MPI_LIBC="$with_mpilibsc"
fi

AC_ARG_VAR(MPIRUN,[MPI run command ])
AC_MSG_CHECKING(for MPIRUN)

if test -z "$MPIRUN" ; then
  if test -n "$MSMPI_BIN" -a -x "$MSMPI_BIN/mpiexec.exe" ;  then
     MPIRUN="$MSMPI_BIN"\mpiexec.exe
  else
    AC_PATH_PROGS(MPIRUN,mpiexec mpirun mpiexec.exe,no)   
  fi
  if test "$MPIRUN" = no
    then
		 ff_mpi=no
	fi
fi
AC_MSG_RESULT($MPIRUN)
if test "ff_mpi" != "no" ; then 
 AC_MSG_CHECKING(for MPIRUN option: )
 ff_mpi_option=""
 ff_okkk=`"$MPIRUN"  -np 2 --oversubscribe echo ff__okkk 2>/dev/null| grep ff__okkk |wc -l`
 if test "$ff_okkk" -eq 2 ; then ff_mpi_option="--oversubscribe" ; fi 
 AC_MSG_RESULT($ff_mpi_option)
fi

AC_MSG_CHECKING(for mpipath )

if test "$with_mpi" != no -a ! -d  "$with_mpipath" -a "$MPIRUN" != no ; then
#   if "$MPIRUN" != no ; tehn
    with_mpipath=`AS_DIRNAME(["$MPIRUN"])`
    with_mpipath=`AS_DIRNAME(["$with_mpipath"])`
#    echo " ***** with_mpipath $with_mpipath \n"

#    else
#    for i in '/c/Program Files (x86)/MPICH2' '/c/Program Files/MPICH2' 'c:\Program Files (x86)\MPICH2' 'c:\Program Files\MPICH2' ; do
#	test -d "$i" &&  with_mpipath="$i" && break
#    done
#    fi
fi
#echo "****  with_mpipath  '$with_mpipath' $MPIRUN *****"
dnl if test "$with_mpilibs" != "no" ; then
dnl fi
case "$MPIRUN" in
 */sgi/mpt/*)
	ff_MPI_INCLUDE_DIR=
	ff_MPI_LIB_DIR=
        test -f "$with_mpipath/include/mpif.h" &&  ff_MPI_INCLUDE_DIR="$with_mpipath/include"
        test -f "$with_mpipath/lib/libmpi.so" &&  ff_MPI_LIB_DIR="$with_mpipath/lib"
        if test -n "$ff_MPI_INCLUDE_DIR" -a -n "$ff_MPI_LIB_DIR" ; then
            ff_MPI_INCLUDE="-I'$ff_MPI_INCLUDE_DIR' "
            with_mpiinc="$ff_MPI_INCLUDE"
            ff_MPI_LIBC="-L'$ff_MPI_LIB_DIR' -lmpi"
            ff_MPI_LIB="-L'$ff_MPI_LIB_DIR' -lmpi++ -lmpi"
            ff_MPI_LIBFC="-L'$ff_MPI_LIB_DIR'  -lmpi"
	    ff_mpitype=sgi
            test -z "$MPICXX" && MPICXX="$CXX $ff_MPI_INCLUDE"
            test -z "$MPIFC"  && MPIFC="$FC  $ff_MPI_INCLUDE"
            test -z "$MPICC"  && MPICC="$CC  $ff_MPI_INCLUDE"
#	    echo " *** MPI sgi ..... "
        fi
	;;
esac


AC_MSG_RESULT($ff_mpi_path)




dnl  correct ff_mpi_path august 2010 -- FH ...


ff_save_cxx="$CXX"
ff_save_libs="$LIBS"


if test "$with_mpi" != no
then
	ff_mpi_path=`AS_DIRNAME(["$MPIRUN"])`
	dnl	echo "ff_mpi_path '$ff_mpi_path' .............."
	case "$ff_mpi_path" in
	    .|"") ff_mpi_path="$PATH";ff_defmpicxx="$ff_mpicxx";;
	    *) ff_mpi_path="$ff_mpi_path";ff_defmpicxx=`expr "//$ff_mpicxx" : '.*/\(.*\)'`;;
	    dnl if also  add $PATH they  could be missing some different mpi version...
	esac
	AC_ARG_VAR(MPICXX,[MPI C++ compiler command])
	if test -z "$MPICXX" ; then
	    AC_PATH_PROGS(MPICXX,$ff_defmpicxx mpic++$ff_mpi_suffix mpicxx$ff_mpi_suffix mpiCC$ff_mpi_suffix mpCC hcp mpxlC mpxlC_r cmpic++,no,$ff_mpi_path)
	    AC_MSG_CHECKING(for MPICXX)
	fi
	ff_mpicxx="eval $MPICXX"
	CXX=$ff_mpicxx
	LIBS="$LIBS $ff_MPI_LIB"
	test -z "$ff_mpi" && ff_mpi=yes
	AC_LINK_IFELSE(
	    [AC_LANG_SOURCE([
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
}])],ff_mpi=yes,ff_mpi=no)
	AC_MSG_RESULT($ff_mpi)

	# Also check that mpirun is there. If it isn't, then MPI is
	# not fully installed.


	if test "$ff_mpi" = yes;
	then

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


	  echo "MPI CC $ff_mpi" >config_LIB_INFO

		# We do not AC_DEFINE any special flag for parallel
		# computation here, because it must only be set when the
 		# parallel program is compiled (see src/mpi/Makfile.am)
		ff_mpiprog="FreeFem++-mpi${EXEEXT}"
   		  AC_SUBST(MPIPROG,"$ff_mpiprog")
   		  AC_SUBST(MPISCRIPT,"ff-mpirun")
   		  AC_SUBST(MPIRUN,"$MPIRUN")
                  AC_SUBST(MPICXX,$MPICXX)
	else
	        AC_SUBST(MPICXX,$ff_save_cxx)
	fi

	if test "$ff_mpi" = yes;
	then
	  if test "$enable_fortran" != no
	  then

	      AC_ARG_VAR(MPIFC,[MPI Fortran 90  compiler command])
	      if test -z "$MPIFC" ; then
		  AC_PATH_PROGS(MPIFC, mpif90$ff_mpi_suffix mpxlf95_r mpxlf90_r mpxlf95 mpxlf90 mpf90 cmpif90c, "",$ff_mpi_path)
	      fi
	      AC_SUBST(MPIFC)
	  fi

#	echo " ********************ffmpi= '$ff_mpi' *************   "
	ff_MPI_INCLUDE="$with_mpiinc"
	if test -z "$ff_mpitype" ; then
           test -n "$MPICXX" && ff_mpishow=`$MPICXX -show` 2>/dev/null
           test -n "$MPICC" && ff_mpicshow=`$MPICC -show` 2>/dev/null
           test -n "$MPIFC" && ff_mpifcshow=`$MPIFC -show` 2>/dev/null
	   if test "$with_mpiinc" = no -o -z "$with_mpiinc" ; then
	      [ff_MPI_INCLUDE=`echo $ff_mpishow|tr ' ' '\n' | sed '1 d'| grep -E '^[-/][^WLlOgpf]|^-Wp,'|tr '\n' ' '`]
	      [ff_mpi_idir=`echo $ff_mpishow|tr ' ' '\n'| grep -E '^-I'|sed s/^-I//|tr '\n' ' '`' /usr/include']
	    else
	      [ff_mpi_idir=`echo $ff_MPI_INCLUDE|tr ' ' '\n'| grep -E '^-I'|sed s/^-I//|tr '\n' ' '`' /usr/include']
	    fi
	    if test "$with_mpilibs" = no -o -z "$with_mpilibs" ; then
		ff_MPI_LIB_DIRS=""
		[ff_MPI_LIB=`echo $ff_mpishow|tr ' ' '\n'| grep -E '^-[Llp]|^-Wl,'|tr '\n' ' '`]
		[ff_MPI_LIBC=`echo $ff_mpicshow|tr ' ' '\n'| grep -E '^-[Llp]|^-Wl,'|tr '\n' ' '`]
		[ff_MPI_LIBFC=`echo $ff_mpifcshow|tr ' ' '\n'| grep -E '^-[Llp]|^-Wl,'|grep -v 'commons,use_dylibs' |tr '\n' ' '`]
            else
	      ff_MPI_LIB="$with_mpilibs"
	      ff_MPI_LIBC="$with_mpilibs"
	      ff_MPI_LIBFC="$with_mpilibs"
	    fi
	    [ff_mpi_ldir=`echo $ff_MPI_LIB|tr ' ' '\n'| grep -E '^-[Llp]|^-Wl,'|sed -e 's/^-[Llp]//' -e 's/^-Wl,]//'  |tr '\n' ' '`' /usr/lib']

	    if  test -z "$ff_MPI_INCLUDE_DIR" ; then
		for i in $ff_mpi_idir; do
		    if test -f "$i/mpi.h" -a -z "$ff_MPI_INCLUDE_DIR"  ;then
			ff_MPI_INCLUDE_DIR=$i
		    fi
		done
	    fi
	    for i in $ff_mpi_ldir; do
	      if test -d $i ; then
		ff_tmp=`ls $i/libmpi.*|head -1`
		if test  -f "$ff_tmp"  -a -z "$ff_MPI_LIB_DIRS"  ;then
		    ff_MPI_LIB_DIRS=$i
		fi
	      fi
	    done
	fi
	AC_SUBST(MPICXX,$MPICXX)
	AC_ARG_VAR(MPICC,[MPI C compiler command in $ff_mpi_path])
	if test -z "$MPICC" ; then
	    AC_PATH_PROGS(MPICC,mpicc$ff_mpi_suffix hcc mpcc mpcc_r mpxlc cmpicc, "",$ff_mpi_path)
	fi
	AC_SUBST(MPICC,$MPICC)

	if test ! -f "$ff_MPI_INCLUDE_DIR/mpif.h"  ; then
	    AC_MSG_NOTICE([ MPI without fortran no file "$ff_MPI_INCLUDE_DIR/mpif.h"  ])
	else
	    if test -n "$MPIFC" ; then
	        AC_FF_ADDWHERELIB(mpifc,$ff_MPI_LIBFC,$ff_MPI_INCLUDE)
	        AC_FF_ADDWHERELIB(mpif77,$ff_MPI_LIBFC,$ff_MPI_INCLUDE)
dnl		  [echo mpifc LD "'$ff_MPI_LIBFC'"   >>$ff_where_lib_conf ]
dnl		  [echo mpifc INCLUDE "'$ff_MPI_INCLUDE'" >>$ff_where_lib_conf ]
dnl		  [echo mpif77 LD "'$ff_MPI_LIBFC'"   >>$ff_where_lib_conf ]
dnl		  [echo mpif77 INCLUDE "'$ff_MPI_INCLUDE'" >>$ff_where_lib_conf ]
	    fi
  	fi
	if test -n "$MPICXX" ; then
            AC_FF_ADDWHERELIB(mpi,$ff_MPI_LIB,$ff_MPI_INCLUDE)
dnl              [echo mpi LD "'$ff_MPI_LIB'"    >>$ff_where_lib_conf ]
dnl              [echo mpi INCLUDE "'$ff_MPI_INCLUDE'" >>$ff_where_lib_conf ]
	fi
	AC_SUBST(MPI_INC_DIR,$ff_MPI_INCLUDE_DIR)
	AC_SUBST(MPI_INCLUDE,$ff_MPI_INCLUDE)
	AC_SUBST(MPI_LIB_DIRS,$ff_MPI_LIB_DIRS)
	AC_SUBST(MPI_LIB,$ff_MPI_LIB)
	AC_SUBST(MPI_LIBC,$ff_MPI_LIBC)
	AC_SUBST(MPI_LIBFC,$ff_MPI_LIBFC)
	AC_SUBST(MPI_RUN_OPTION,$ff_mpi_option)
        AC_SUBST(SKIP_TESTS_MPI,"no")
	fi
	CXX="$ff_save_cxx"
	LIBS="$ff_save_libs"
fi
fi
##  clean on MPI variable if not MPI ...
if test "$ff_mpi" != yes ; then

	  AC_SUBST(MPIRUN,"")
	  AC_SUBST(MPICC,"")
	  AC_SUBST(MPICXX,"")
	  AC_SUBST(MPIFC,"")
	  AC_SUBST(MPI_INCLUDE,"")
	  AC_SUBST(MPI_LIB_DIRS,"")
	  AC_SUBST(MPI_LIB,"")
	  AC_SUBST(MPI_LIBC,"")
	  AC_SUBST(MPI_LIBFC,"")
          AC_SUBST(SKIP_TESTS_MPI,"yes")
	  ff_mpi=no
dnl    AC_MSG_ERROR([ Sorry nompi  compiler !])
fi

# Local Variables:
# mode:shell-script
# ispell-local-dictionary:"british"
# coding:utf-8
# End:
