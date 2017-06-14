#!/bin/sh
# Create a loadable object from a C++ function defined in a .cpp file
# $Id$
dir=..
do="yes"
DEBUG=""
uu=`uname -s` 
INC=""
LIBS=""
DLL=""
bin="."
while true; do
case "$1" in
 -[h?]*) echo usage $0 "[-n] [-g] [-win32] [-l libfile] [-I includedir]  fileprefixe"
 echo "    -n  :  do nothing just print"
 echo "    -g  :  compile with -g option"
 echo "    -win33: compile for win32 OS (Window XP, ...)"
 echo "    -l files  add files to the ld process (link)"
 echo "    -I dir  add dir in include seach dir for compilation"
 echo "    -dll file  add  dll and copie  in $bin "
 exit 0;
;;
 -*) ;;
 *) break;
esac
  
if [ "$1" = "-n" ]; then
  shift
  do="no"
elif [ "$1" = "-g" ]; then
  shift
  DEBUG="-g"
elif [ "$1" = "-win32" ]; then
  shift
  bin="$dir"
  uu="WIN32"
elif  [ "$1" = "-I" ]; then
  INC="$INC -I$2"; shift;shift;
elif  [ "$1" = "-l" ]; then
  LIBS="$LIBS $2"; shift;shift;
elif  [ "$1" = "-dll" ]; then
  DLL="$DLL $2"; shift;shift;
else
  break;
fi  
done
# Default compiler
if [ "$CXX" = "" ];
then
    CXX=g++
fi

INC="-Iinclude $INC"  
SUF=so
if [ -f "$1.cpp" ] ; then
    case "$uu" in
	Darwin) 
	export MACOSX_DEPLOYMENT_TARGET=10.3
	SUF=dylib
	SHARED="-bundle -undefined dynamic_lookup"  ;;
	CYGWIN*|FreeBSD)
	    SHARED="-shared " ;;

	# 64 bit Linux needs -fPIC (ALH)
	Linux)
	    FLAGS='-fPIC'
	    SHARED="-shared " ;;

        WIN32) 
        echo " window "
        b=$bin
	SHARED="-v -shared --unresolved-symbols=ignore-all"
        FLAGS='  -mno-cygwin '
        LIBS="$b/libff0.dll $b/libff1.dll $b/libff2.dll $LIBS $DLL"
        SUF=dll;;
	*)
	echo "sorry unknown achitecture "`uname`
	exit 1
    esac
   FLAGS="$FLAGS -g"
   echo $CXX -c $FLAGS $INC $PIC $1.cpp 
   test $do = yes &&$CXX -c $INC $FLAGS $PIC  $1.cpp 

   echo $CXX $SHARED $FLAGS $1.o -o $bin/$1.$SUF $LIBS $DLL
   test $do = yes &&$CXX $SHARED $FLAGS $1.o -o $bin/$1.$SUF $LIBS $DLL
   if [ -n "$DLL" ] ; then
     echo cp $DLL $bin
     test $do = yes &&  cp $DLL $bin
   fi
else
    echo "sorry file $1.cpp does not exist"
fi
