#!/bin/sh
appl=/Applications
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


echo " copy FreeFem++.app in "$appl" "
if [ -d FreeFem++.app ] ; then
  rsync -avHE --delete   FreeFem++.app/  "$appl"/FreeFem++.app
fi


lbin=`cd $appl"/FreeFem++.app/Contents/bin/; echo *`
echo "  install $lbin    commands in /usr/local/bin (need of admin password)"

sudo mkdir -p  /usr/local/bin
sudo ln -s  "$appl"/FreeFem++.app/Contents/bin/* /usr/local/bin
sudo rm usr/local/bin/ff-c++
sudo sed <"$appl"/FreeFem++.app/Contents/bin/ff-c++ >/usr/local/bin/ff-c++ \
	-e 's;FFAPPLI_INC;$app/FreeFem++.app/Contents/include;' 

echo "++  FreeFem++ is correctly install in $appl  directory."