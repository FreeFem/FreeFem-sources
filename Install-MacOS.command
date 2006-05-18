#!/bin/sh
if [ -f OtherMacOsLib.tgz ]; then
for i in `tar ztf OtherMacOsLib.tgz`; do
 if [ ! -f "/$i" ]; then
  echo " the Libary '/$i' don't exist => install (need of admin password)"
  (cd /;sudo echo tar zxvf OtherMacOsLib.tgz $i)
  else
  echo " the Libary '/$i'  exist "
 fi
done 
fi
if [ -f  FreeFem++-CoCoa ] ;then
echo "  install FreeFem++-CoCoa script in /usr/local/bin (need of admin password)"
sudo cp FreeFem++-CoCoa /usr/local/bin
fi
echo " copy FreeFem++.app in /Apllications "
if [ -d FreeFem++.app ] ; then
  rsync -av  FreeFem++.app/  /Applications/FreeFem++.app
fi