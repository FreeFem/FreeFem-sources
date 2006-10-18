#!/bin/sh
cd `dirname $0`
echo "Installtion of Freefem++ " 
if [ -f OtherMacOsLib.tgz ]; then
for i in `tar ztf OtherMacOsLib.tgz`; do
 if [ ! -f "/$i" ]; then
  echo " the Libary '/$i' don't exist => install (need of admin password)"
   sudo  tar zxvf OtherMacOsLib.tgz -C / $i
  else
  echo " the Libary '/$i'  exist "
 fi
done 
#  verif .... 
for i in `tar ztf OtherMacOsLib.tgz`; do
  if [ ! -f "/$i" ]; then
    echo " the Libary '/$i' don't exist FreeFEM cannot run (call you adminisator sorry)"
    echo "Sorry"
    exit 1;
 fi
done
fi

if [ -f  FreeFem++-CoCoa ] ;then
echo "  install FreeFem++-CoCoa script in /usr/local/bin (need of admin password)"
sudo cp FreeFem++-CoCoa /usr/local/bin
fi
echo " copy FreeFem++.app in /Applications "
if [ -d FreeFem++.app ] ; then
  rsync -avHE --delete   FreeFem++.app/  /Applications/FreeFem++.app
fi
echo "++  FreeFem++ is correctly install in /Application directory."