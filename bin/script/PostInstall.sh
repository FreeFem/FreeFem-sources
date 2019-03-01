#!/bin/sh
# "-DFF__FVER=$(PACKAGE_VERSION)" 
# "-DFF_BINDIR=$(bindir)" 
# "-DFF__DATADIR=$(pkgdatadir)
#  "FFBIN="@prefix@"/bin
ff_desktop="$HOME/Desktop/FreeFem++-""4.0"
mkdir -p -m 0755 /etc/paths.d
ln -sf "/usr/local/share/freefem++"/"freefem++doc.pdf" "$HOME/Desktop"
test -e "$ff_desktop" || ln -sf "/usr/local/share/freefem++"/"4.0" "$ff_desktop"
echo Install /etc/paths.d/FreeFem++ file:  "/usr/local/bin"

echo "/usr/local/bin" > /etc/paths.d/FreeFem++
chmod a+r /etc/paths.d/FreeFem++
echo " Try to Clean old file version "
if [ -d  /usr/local/bin ] ; then  
  cd /usr/local/bin
  for i in  FreeFem++ FreeFem++-CoCoa FreeFem++-mpi FreeFem++-nw bamg cvmsh2 ff-c++ ff-get-dep ff-mpirun ff-pkg-download ffglut ffmedit; 
  do 

      if [  -f  "$i" ] ; then 
	  echo " clean $i "
	  rm "$i";
      fi
  done


echo ln -s /usr/local/bin/FreeFem++-CoCoa  /usr/local/bin/ 
ln -s /usr/local/bin/FreeFem++-CoCoa  /usr/local/bin/ 

fi
# bluid new link to new 


