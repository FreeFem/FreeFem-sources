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
    MPIF77="$F77 $ff_MPI_INCLUDE"
    MPIFC="$FC  $ff_MPI_INCLUDE"
    MPICC="$CC  $ff_MPI_INCLUDE"
    AC_MSG_NOTICE([   ---  set  all MPI compile to compiler:   $MPICC , $MPIF77, $MPIFC, $MPICC  ])
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
    AC_PATH_PROGS(MPIRUN,mpirun mpiexec mpiexec.exe,no)   
  fi
  if test "$MPIRUN" = no
    then
		 ff_mpi=no
	fi
fi
AC_MSG_RESULT($MPIRUN)
if test "ff_mpi" != "no" ; then 
 AC_MSG_CHECKING(for MPIRUN option: )
 ff_okkk=`"$MPIRUN"  -np 2 --oversubscribe echo ff__okkk| grep ff__okkk |wc -l`
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
            test -z "$MPIF77" && MPIF77="$F77 $ff_MPI_INCLUDE"
            test -z "$MPIFC"  && MPIFC="$FC  $ff_MPI_INCLUDE"
            test -z "$MPICC"  && MPICC="$CC  $ff_MPI_INCLUDE"
#	    echo " *** MPI sgi ..... "
        fi
	;;
esac
#	echo " ####   --$MSMPI_INC--$MSMPI_BIN--$ff_win32"
 if test -n "$MSMPI_INC" -a -n "$MSMPI_BIN" -a  "$ff_win32" = yes   ; then
   echo " ####  check  MSMPI"
		   # MSMPI_LIB64 MSMPI_LIB32   $ff_ptrbit is  32 or 64
		   ffMSMPI_BIN=`cygpath $MSMPI_BIN` 
		   ffMSMPI_INC=`cygpath $MSMPI_INC` 
		   ffMSMPI_LIB32=`cygpath $MSMPI_LIB32` 
		   ffMSMPI_LIB64=`cygpath $MSMPI_LIB64` 
		   
		   mkdir -p 3rdparty/include/msmpi
		   mkdir -p 3rdparty/lib/msmpi
		   ## add to msmpi 10.0
		   echo " hack MSMPI V10.0 "
		   echo "void __guard_check_icall_fptr(unsigned long ptr) { }" > 3rdparty/lib/msmpi/cfg_stub.c
		   gcc -o 3rdparty/lib/msmpi/cfg_stub.o -c 3rdparty/lib/msmpi/cfg_stub.c
		   gcc -shared -o 3rdparty/lib/msmpi/cfg_stub.dll 3rdparty/lib/msmpi/cfg_stub.o
		   rm 3rdparty/lib/msmpi/cfg_stub.o
		   rm 3rdparty/lib/msmpi/cfg_stub.c
		   
		   cp "$MSMPI_INC"/*.h 3rdparty/include/msmpi
		   grep -v INT_PTR_KIND "$MSMPI_INC"/mpif.h >3rdparty/include/msmpi/mpif.h
		   test "$ff_ptrbit" -eq 64 && cp "$MSMPI_INC"/x64/*.h 3rdparty/include/msmpi
		   test "$ff_ptrbit" -eq 32 && cp "$MSMPI_INC"/x86/*.h 3rdparty/include/msmpi
		   ff_MPI_INCLUDE_DIR=`pwd`/3rdparty/include/msmpi
		   ff_msmpi_lib="$MSMPI_LIB64"
		   test "$ff_ptrbit" -eq 32 && ff_msmpi_lib="$MSMPI_LIB32"
		   cp "$ff_msmpi_lib/msmpifec.lib" "$ff_msmpi_lib/msmpi.lib" 3rdparty/lib/msmpi
		   ff_msmpi_lib=`pwd`/3rdparty/lib/msmpi
		   
		   
		   # to reinstall msmpi .. 
		   
#  MSMPI
    if test -x "`which msmpi.dll`"
    then
#  Remove for scotch and parmetis
	ff_MPI_INCLUDE="-I$ff_MPI_INCLUDE_DIR  -D__int64=long\ long"
	with_mpiinc="$ff_MPI_INCLUDE"
	test -z "$MPIRUN" -a -x "$ffMSMPI_BIN/mpiexe.exe" && MPIRUN="$MSMPI_BIN\mpiexe.exe"
	ff_MPI_LIBC="'$ff_msmpi_lib/msmpi.lib'"
	ff_MPI_LIB="'$ff_msmpi_lib/msmpi.lib'"
	ff_MPI_LIBFC="'$ff_msmpi_lib/msmpifec.lib' '$ff_msmpi_lib/msmpi.lib' '$ff_msmpi_lib/cfg_stub.dll' "
	ff_mpiexec_win="C:\Program Files\Microsoft MPI\Bin\mpiexec.exe"
	test -z "$ff_mpiexec_win" && MPIRUN="$ff_mpiexec_win"
	test -z "$MPICXX" && MPICXX="$CXX $ff_MPI_INCLUDE"
	test -z "$MPIF77" && MPIF77="$F77 $ff_MPI_INCLUDE"
	test -z "$MPIFC"  && MPIFC="$FC  $ff_MPI_INCLUDE"
	test -z "$MPICC"  && MPICC="$CC  $ff_MPI_INCLUDE"
	ff_mpitype=MSMPI
    else
	echo " #### no msmpi.dll  => no mpi under windows .... (FH) " >&AS_MESSAGE_LOG_FD
	echo " #### no msmpi.dll  => no mpi under windows .... (FH) " >&AS_MESSAGE_FD
	with_mpipath=no
	with_mpi=no
    fi
elif test  -d "$with_mpipath" -a "$ff_win32" = yes  ; then
#    sed -e "s?@MPIDIR@?$with_mpipath?" -e "s?@F77@?$F77?" -e "s?@CC@?$CC?" -e "s?@CXX@?$CXX?"   -e "s?@FC@?$FC?"  <mpic++.in >mpic++
 #   chmod a+rx mpic++
  #  for i in mpicc mpif90 mpifc mpif77 ; do cp mpic++ $i; done
#    ff_pwd=`pwd`
 #   with_mpi="$ff_pwd"/mpic++
 #   MPICXX="$ff_pwd/mpic++"
 #   MPIF77="$ff_pwd/mpif77"
 #   MPIFC="$ff_pwd/mpif90"
 #   MPICC="$ff_pwd/mpicc" zzzzzzzzzzz
    if  with_mpilibs=`which msmpi.dll`
    then
	case "$ff_size_ptr"  in
	    4) with_mpipathlib="$with_mpipath/Lib/i386";;
	    8) with_mpipathlib="$with_mpipath/Lib/amd64";;
	    *) with_mpipath=no;;
	esac


	test -d "$with_mpipath/Inc" &&  ff_MPI_INCLUDE_DIR="$with_mpipath/Inc"
	test -d "$with_mpipath/Include" &&  ff_MPI_INCLUDE_DIR="$with_mpipath/Include"
#  Remove for scotch and parmetis
#	ff_MPI_INCLUDE="-I'$ff_MPI_INCLUDE_DIR' '-D_MSC_VER' '-D__int64=long long'"
	ff_MPI_INCLUDE="-I'$ff_MPI_INCLUDE_DIR'  '-D__int64=long long'"
	with_mpiinc="$ff_MPI_INCLUDE"
	test -z "$MPIRUN" && MPIRUN="$with_mpipath/bin/mpiexe.exe"
	ff_MPI_LIBC="$with_mpilibs"
	ff_MPI_LIB="$with_mpilibs"
	ff_MPI_LIBFC="$with_mpilibs"
	test -z "$MPICXX" && MPICXX="$CXX $ff_MPI_INCLUDE"
	test -z "$MPIF77" && MPIF77="$F77 $ff_MPI_INCLUDE"
	test -z "$MPIFC"  && MPIFC="$FC  $ff_MPI_INCLUDE"
	test -z "$MPICC"  && MPICC="$CC  $ff_MPI_INCLUDE"
    else
	echo " #### no msmpi.dll  => no mpi under windows .... (FH) " >&AS_MESSAGE_LOG_FD
	echo " #### no msmpi.dll  => no mpi under windows .... (FH) " >&AS_MESSAGE_FD
	with_mpipath=no
	with_mpi=no
    fi
else
    with_mpipath=no
fi


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

	      AC_ARG_VAR(MPIF77,[MPI Fortran 77 compiler command])
	      if test -z "$MPIF77" ; then
		  AC_PATH_PROGS(MPIF77, mpif90$ff_mpi_suffix mpif77$ff_mpi_suffix hf77 mpxlf mpf77 mpif90 mpf90 mpxlf90 mpxlf95 mpxlf_r cmpifc cmpif90c, "",$ff_mpi_path)
	      fi
	      AC_SUBST(MPIF77)
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
	    if test "$with_mpilibs" = no -o -z "$with_mpilibs" ; then
		[ff_MPI_INCLUDE=`echo $ff_mpishow|tr ' ' '\n'| grep -E '^[-/][^WLlOgp]|^-Wp,'|tr '\n' ' '`]
		ff_MPI_LIB_DIRS=""
		[ff_MPI_LIB=`echo $ff_mpishow|tr ' ' '\n'| grep -E '^-[Llp]|^-Wl,'|tr '\n' ' '`]
		[ff_MPI_LIBC=`echo $ff_mpicshow|tr ' ' '\n'| grep -E '^-[Llp]|^-Wl,'|tr '\n' ' '`]
		[ff_MPI_LIBFC=`echo $ff_mpifcshow|tr ' ' '\n'| grep -E '^-[Llp]|^-Wl,'|grep -v 'commons,use_dylibs' |tr '\n' ' '`]
		[ff_mpi_idir=`echo $ff_mpishow|tr ' ' '\n'| grep -E '^-I'|sed s/^-I//|tr '\n' ' '`' /usr/include']
	    fi
	    [ff_mpi_idir=`echo $ff_MPI_INCLUDE|tr ' ' '\n'| grep -E '^-I'|sed s/^-I//|tr '\n' ' '`' /usr/include']
	    [ff_mpi_ldir=`echo $ff_MPI_LIB|tr ' ' '\n'| grep -E '^-[Llp]|^-Wl,'|sed -e 's/^-[Llp]//' -e 's/^-Wl,]//'  |tr '\n' ' '`' /usr/lib']

	    if  test -z "$ff_MPI_INCLUDE_DIR" ; then
		for i in $ff_mpi_idir; do
		    if test -f "$i/mpi.h" -a -z "$ff_MPI_INCLUDE_DIR"  ;then
			ff_MPI_INCLUDE_DIR=$i
		    fi
		done
	    fi
	    for i in $ff_mpi_ldir; do
		ff_tmp=`ls $i/libmpi.*|head -1`
		if test  -f "$ff_tmp"  -a -z "$ff_MPI_LIB_DIRS"  ;then
		    ff_MPI_LIB_DIRS=$i
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
	  AC_SUBST(MPIF77,"")
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
