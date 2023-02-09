#!/bin/bash 
# Apple
# hdiutil create -srcfolder /Users/hecht/Desktop/FF-dmg/FreeFEM-4.12-Apple-M1 -fs HFS+  FreeFEM-4.12-Apple-M1.dmg
ppwd=`dirname $0`
#sudo cp -rf FreeFem++.app /Applications
#sudo xattr -rc /Applications/FreeFem++.app

cd 
filerc=".profile"
case "$SHELL" in
	*zsh) filerc=".zprofile";;
	*bash) filerc=".bash_profile";;
	*) echo erreur SHELL $SHELL; exit 1;;
esac;
echo filerc : $filerc 
ver=4.12
arch=M1
dirff=/Applications/FreeFem++.app/Contents/ff-$ver/bin
ffexec=/Applications/FreeFem++.app/Contents/ff-$ver/bin/FreeFem++
ff=$(which FreeFem++)
dd=$(dirname "$ff")

if [ ! -d   /Applications/FreeFem++.app/Contents/ff-$ver/ ] ; then
	echo " No FreeFem++.app  $ff-ver so copy FreeFem++.app in /Application (with sudo do be sure)"
	echo =============================
	sudo cp -rf /Volumes/FreeFEM-$ver-Apple-$arch/FreeFem++.app /Applications
	sudo cp /Applications/FreeFem++.app/Contents/ff-4.12/bin/FreeFem++-CoCoa  /usr/local/bin
fi

if test ! -d "$dd"  ; then
	echo " bug no directory  $dd !"
	exit 1; 
fi

if [ -d  $dirff  ] ; then
	# verification of com.apple.quarantine
	quarantine=`xattr -l /Applications/FreeFem++.app/Contents/ff-$ver/bin/FreeFem++|grep com.apple.quarantine|wc -l`
	if [  $quarantine -ne 0 ] ; then
		echo " remove: com.apple.quarantine from FreeFem++ app (need sudo )"
		echo do sudo xattr -rc /Applications/FreeFem++.app/
		sudo xattr -rc /Applications/FreeFem++.app/
		echo " verif quarantine "
		xattr -l /Applications/FreeFem++.app/Contents/ff-$ver/bin/FreeFem++
	fi
else 
	echo " OK : No quarantine "
fi
# echo test FreeFem++ 
/Applications/FreeFem++.app/Contents/ff-$ver/bin/FreeFem++  /Applications/FreeFem++.app/Contents/ff-$ver/share/FreeFEM/$ver/examples/tutorial/Laplace.edp 

if [ $dd != $dirff ] ; then 
	echo update $filerc 
	echo =============================
	
	test -e $filerc || touch $filerc 
	echo add export PATH=$dirff:\$PATH " >>$filerc"
	grep -v FreeFem++.app/Contents $filerc >$filerc.new
	echo export PATH=$dirff:\"\$PATH\" >>$filerc.new
	diff -u $filerc.new $filerc
	if test $? -eq  1; then
		mv -f $filerc.new $filerc
		chmod a+x $filerc
	else
		rm $filerc.new
	fi
	
fi

cd ~/Desktop
echo Create Link of example and Doc 
echo =============================
test -L  'FreeFEM doc and Examples' && rm 'FreeFEM doc and Examples'
test -e 'FreeFEM doc and Examples' || ln -s  /Applications/FreeFem++.app/Contents/ff-$ver/share/FreeFEM/ 'FreeFEM doc and Examples'
cd -


vvv=`$dirff/FreeFem++ -nw| grep 'version'| wc -l`
if [ $vvv -eq 0 ]; then
  echo Error missing install missing lib ??
else
	echo FreeFem++ work
fi