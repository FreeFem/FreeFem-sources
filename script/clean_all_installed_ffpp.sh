#!/bin/sh 
# clean all freefem++ version ...
cd / 
if [ -d Applications/FreeFem++.app ] ; then
  rm -rf Applications/FreeFem++.app
fi
if [ -d usr/local/share/freefem++ ] ; then
  rm -rf usr/local/share/freefem++
fi
if [ -d  usr/local/lib/ff++/ ]; then
 rm -rf usr/local/lib/ff++
fi
rm usr/local/bin/FreeFem++		usr/local/bin/FreeFem++-nw	usr/local/bin/ff-c++		usr/local/bin/ff-pkg-download 
rm usr/local/bin/FreeFem++-CoCoa		usr/local/bin/bamg	usr/local/bin/ff-get-dep
rm  usr/local/bin/ffglut usr/local/bin/FreeFem++-mpi		usr/local/bin/cvmsh2	
rm usr/local/bin/ff-mpirun  usr/local/bin/ffmedit
exit 0 
