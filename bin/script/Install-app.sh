#!/bin/bash 
# Apple
# hdiutil create -srcfolder /Users/hecht/Desktop/FF-dmg/FreeFEM-4.12-Apple-M1 -fs HFS+  FreeFEM-4.12-Apple-M1.dmg
scriptpath=`realpath $0`
ppwd=`dirname $scriptpath`
vffapp=$ppwd/FreeFem++.app
vffver=$ppwd/FreeFem++.app/Contents/ff-*
ffver=`basename $vffver`
vffexe=$vffver/bin/FreeFem++	
tyexe=`file -b $vffexe`
case $tyexe-`arch`	 in
 *arm64*-*arm64*) echo  same arch arm64 arm64  ok;;
 *x86*-*x86*)  echo  same arch x86 x85  ok;;
 *) echo bad arch ; exit 1;;
 esac
#sudo cp -rf FreeFem++.app /Applications
#sudo xattr -rc /Applications/FreeFem++.app
# echo verif arch 
cd 
filerc=".profile"
case "$SHELL" in
	*zsh) filerc=".zprofile";;
	*bash) filerc=".bash_profile";;
	*) echo erreur SHELL $SHELL; exit 1;;
esac;
echo filerc : $filerc 

dirff=/Applications/FreeFem++.app/Contents/$ffver/bin
ffexec=/Applications/FreeFem++.app/Contents/$ffver/bin/FreeFem++
ff=$(which FreeFem++)
dd=$(dirname "$ff")
exit 0
if [ ! -d   /Applications/FreeFem++.app/Contents/$ffver/ ] ; then
	echo " No FreeFem++.app  $ffver so copy FreeFem++.app in /Application (with sudo do be sure)"
	echo =============================
	sudo cp -rf $vffapp/FreeFem++.app /Applications
	sudo cp /Applications/FreeFem++.app/Contents/$ffver/bin/FreeFem++-CoCoa  /usr/local/bin
fi

if test ! -d "$dd"  ; then
	echo " bug no directory  $dd !"
	exit 1; 
fi

if [ -d  $dirff  ] ; then
	# verification of com.apple.quarantine
	quarantine=`xattr -l /Applications/FreeFem++.app/Contents/$ffver/bin/FreeFem++|grep com.apple.quarantine|wc -l`
	if [  $quarantine -ne 0 ] ; then
		echo " remove: com.apple.quarantine from FreeFem++ app (need sudo )"
		echo do sudo xattr -rc /Applications/FreeFem++.app/
		sudo xattr -rc /Applications/FreeFem++.app/
		echo " verif quarantine "
		xattr -l /Applications/FreeFem++.app/Contents/$ffver/bin/FreeFem++
	fi
else 
	echo " OK : No quarantine "
fi
# echo test FreeFem++ 
/Applications/FreeFem++.app/Contents/$ffver/bin/FreeFem++  /Applications/FreeFem++.app/Contents/$ffver/share/FreeFEM/$ver/examples/tutorial/Laplace.edp 

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
test -e 'FreeFEM doc and Examples' || ln -s  /Applications/FreeFem++.app/Contents/$ffver/share/FreeFEM/ 'FreeFEM doc and Examples'
cd -


vvv=`$dirff/FreeFem++ -nw| grep 'version'| wc -l`
if [ $vvv -eq 0 ]; then
  echo Error missing install missing lib or quarantine ?
else
	echo FreeFem++ work
fi