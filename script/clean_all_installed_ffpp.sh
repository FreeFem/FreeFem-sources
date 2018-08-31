#!/bin/sh 
# clean all freefem++ version ...
rm="rm"

cd / 
if [ -d Applications/FreeFem++.app ] ; then
  $rm -rf Applications/FreeFem++.app
fi
if [ -d usr/local/share/freefem++ ] ; then
  $rm -rf usr/local/share/freefem++
fi
if [ -d  usr/local/lib/ff++/ ]; then
 $rm -rf usr/local/lib/ff++
fi

ffexe="FreeFem++*  ff-c++ ff-pkg-download  bamg ff-get-dep ffglut cvmsh2 ff-mpirun ffmedit"
( cd usr/local/bin &&echo cd "usr/local/bin" &&  $rm $ffexe; )

if [ -f  etc/paths.d/FreeFem++  ]; then
	ffbin=`cat etc/paths.d/FreeFem++`
	( cd "$ffbin" &&echo cd "$ffbin" && $rm $ffexe; )
fi
if [ -d usr/local/ff++ ] ; then
	echo " Warning dir usr/local/ff++ exist"
	ls usr/local/ff++
	echo " Warning the directory  /usr/local/ff++  no remove ? "
fi
exit 0 
