#!/bin/sh
# "-DFF__FVER=$(PACKAGE_VERSION)" 
# "-DFF_BINDIR=$(bindir)" 
# "-DFF__DATADIR=$(pkgdatadir)
#  "FFBIN="@prefix@"/bin

if [ "$(uname)" = "Darwin" ]; then
  ff_desktop="$HOME/Desktop/FreeFem++-""FF__FVER"
  mkdir -p -m 0755 /etc/paths.d
  ln -sf "FF__DATADIR"/"FreeFEM-documentation.pdf" "$HOME/Desktop"
  test -e "$ff_desktop" || ln -sf "FF__DATADIR"/"FF__FVER" "$ff_desktop"
  echo Install /etc/paths.d/FreeFem++ file:  "FF_BINDIR"

  echo "FF_BINDIR" > /etc/paths.d/FreeFem++
  chmod a+r /etc/paths.d/FreeFem++
fi

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

  if [ "$(uname)" = "Darwin" ]; then
      echo ln -s FF_BINDIR/FreeFem++-CoCoa  /usr/local/bin/ 
      ln -s FF_BINDIR/FreeFem++-CoCoa  /usr/local/bin/ 
  fi
fi
# bluid new link to new 


