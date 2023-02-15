#!/bin/bash
# FF_DMGNAME
# FF_APP
include(`bin/init-shell.sh')
# Apple
# hdiutil create -srcfolder /Users/hecht/Desktop/FF-dmg/FreeFEM-4.12-Apple-M1 -fs HFS+  FreeFEM-4.12-Apple-M1.dmg

cd 
filerc=".profile"
case "$SHELL" in
	*zsh) filerc=".zprofile";;
	*bash) filerc=".bash_profile";;
	*) echo erreur SHELL $SHELL; exit 1;;
esac;
echo filerc : $filerc 
dirapp=$(realpath $ff_prefix/../..)
case $dirapp in 
   *.app) ;;
   *)  echo "$dirapp error is not .app directory !!!!"
       exit 1;;
esac;

#  verif 
ffexec=$bindir/FreeFem++
if [ -x $ffexec ]; then
fi
  

if [ ! -d   $ff_prefix  ] ; then
	echo " No $ff_prefix  so copy to $ff_prefix  form dmg (with sudo do be sure)"
	echo =============================
	sudp mkdir -p $ff_prefix
	sudo cp -rf /Volumes/FF_DMGNAME/FF_APP  $ff_prefix/../..
	sudo cp $bindir/FreeFem++-CoCoa  /usr/local/bin
fi

if test ! -d "$bindir"  ; then
	echo " bug no directory  $bindir !"
	exit 1; 
fi

if [ -d  $dirff  ] ; then
	# verification of com.apple.quarantine
	quarantine=`xattr -l $bindir/FreeFem++|grep com.apple.quarantine|wc -l`
	if [  $quarantine -ne 0 ] ; then
		echo " remove: com.apple.quarantine from FF_APP  (need sudo )"
		echo do sudo xattr -rc $ff_prefix/../..
		sudo xattr -rc $ff_prefix/../..
		echo " verif quarantine "
		xattr -l $bindir/FreeFem++
	fi
else 
	echo " OK : No quarantine "
fi
# echo test FreeFem++
echo try  $bindir/FreeFem++  $ff_data/examples/tutorial/Laplace.edp 
$bindir/FreeFem++  $ff_data/examples/tutorial/Laplace.edp 
if [ $? -ne 0 ]; then
   echo open "Preference System / Confidentiality and Security menu to unlock"
   echo $bindir/FreeFem++ 
fi
if [ $dd != $bindir ] ; then 
	echo update $filerc 
	echo =============================
	
	test -e $filerc || touch $filerc 
	echo add export PATH=$dirff:\$PATH " >>$filerc"
	grep -v FreeFem++ $filerc >$filerc.new
	echo "# add FreeFem++ $version in PATH " >>$filerc.new
	echo export PATH=$bindir:\"\$PATH\" >>$filerc.new
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
test -e 'FreeFEM doc and Examples' || ln -s  $ff_data/..  'FreeFEM doc and Examples'
cd -


vvv=`$bindir/FreeFem++ -nw| grep 'version'| wc -l`
if [ $vvv -eq 0 ]; then
  echo Error missing install missing lib ?? or locked application
  echo open "Preference System / Confidentiality and Security menu to unlock"
  echo $bindir/FreeFem++ 
  
else
	echo FreeFem++ work's
fi