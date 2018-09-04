#!/bin/bash

############################################################################
# This file is part of FreeFem++.                                          #
#                                                                          #
# FreeFem++ is free software: you can redistribute it and/or modify        #
# it under the terms of the GNU Lesser General Public License as           #
# published by the Free Software Foundation, either version 3 of           #
# the License, or (at your option) any later version.                      #
#                                                                          #
# FreeFem++ is distributed in the hope that it will be useful,             #
# but WITHOUT ANY WARRANTY; without even the implied warranty of           #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
# GNU Lesser General Public License for more details.                      #
#                                                                          #
# You should have received a copy of the GNU Lesser General Public License #
# along with FreeFem++. If not, see <http://www.gnu.org/licenses/>.        #
############################################################################
# SUMMARY : Deploy FreeFem++ binary on Ubuntu 16.04
# LICENSE : LGPLv3
# ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
# AUTHORS : Simon Garnotel
#           Pierre Jolivet
#           Frederic Hecht
# E-MAIL  : simon.garnotel@gmail.com

## TODO: hdf5 does not work properly

# Colors
bold="\e[1m"
redBold="\e[1m\e[31m"
greenBold="\e[1m\e[32m"
yellowBold="\e[1m\e[33m"
reset="\e[0m"

# Error functions
errCopy=0
errConfig=0
errBashrc=0
checkErrorCopy () {
	if [ $? -eq 0 ]
	then
		echo -e "$greenBold INFO:$reset Copy success"
	else
		echo -e "CRITICAL:$reset Copy failed"
		((errCopy++))
	fi
}

checkErrorConfig () {
	if [ $? -eq 0 ]
	then
		echo -e "$greenBold INFO:$reset Configuration file success"
	else
		echo -e "$redBold CRITICAL:$reset Configuration file failed"
		echo -e "\tSome functionality will not work properly"
		((errConfig++))
	fi
}

checkErrorBashrc () {
	if [ $? -eq 0 ]
	then
		echo -e "$greenBold INFO:$reset bashrc file success"
	else
		echo -e "$redBold CRITICAL:$reset bashrc file failed"
		echo -e "\tSome functionality will not work properly"
		((errBashrc++))
	fi
}

# Copy binaries to /usr/local/*
echo -e "$bold Deploy$reset"
cp freefem-source-master/bin/* /usr/local/bin/
checkErrorCopy
cp -R freefem-source-master/ff-petsc /usr/local/
checkErrorCopy
cp -R freefem-source-master/lib/* /usr/local/lib/
checkErrorCopy
cp -R freefem-source-master/share/* /usr/local/share/
checkErrorCopy

# Create /etc/freefem++.pref
echo -e "$bold Create configuration file$reset"
echo -e "verbosity = 5" > /etc/freefem++.pref
checkErrorConfig
echo -e "#LOAD" >> /etc/freefem++.pref
echo -e "loadpath += \"/usr/local/lib/ff++/3.61/lib\"" >> /etc/freefem++.pref
checkErrorConfig
echo -e "#INCLUDE" >> /etc/freefem++.pref
echo "includepath += \"/usr/local/lib/ff++/3.61/idp\"" >> /etc/freefem++.pref
checkErrorConfig

# Append to bashrc
echo -e "$bold Append library path to bashrc$reset"
if grep -qe 'export LD_LIBRARY_PATH="*:/usr/local/ff-petsc/real/lib:/usr/local/ff-petsc/complex/lib"' /etc/bash.bashrc 
then
	echo -e "\tbashrc already contains FreeFem++ libraries"
else
	echo -e "" >> /etc/bash.bashrc
	checkErrorBashrc
	echo -e "export LD_LIBRARY_PATH=\"$LD_LIBRARY_PATH:/usr/local/ff-petsc/real/lib:/usr/local/ff-petsc/complex/lib\"" >> /etc/bash.bashrc
	checkErrorBashrc
fi

# Summary
echo -e "$bold Summary$reset"
echo -e "\tCopy error: $errCopy"
echo -e "\tFile error: $errConfig"
echo -e "\tBashrc error: $errBashrc"

if [ $errCopy -gt 0 ]
then
	echo -e "$redBold Some copy errors occurs$reset"
fi
if [ $errConfig -gt 0 ]
then
	echo -e "$redBold Some configuration errors occurs$reset"
fi
if [ $errBashrc -gt 0 ]
then
	echo -e "$redBold Some bashrc errors occurs$reset"
fi

if [ $errCopy -eq 0 ] && [ $errConfig -eq 0 ] && [ $errBashrc -eq 0 ]
then
	echo -e "$greenBold You can now use the FreeFem++ command directly in your terminal$reset"
	echo -e " Complete integration needs a restart"
else
	echo -e "$yellowBold Try to use the FreeFem++ command in your terminal$reset"
	echo -e "$redBold If that does not work, please post an issue in https://github.com/FreeFem/FreeFem-sources/issues$reset"
fi



