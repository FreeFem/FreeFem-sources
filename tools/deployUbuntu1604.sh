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
##       metis does not work (is it compiled in the CI)

# Error functions
errCopy=0
errConfig=0
errBashrc=0
checkErrorCopy () {
	if [ $? -eq 0 ]
	then
		echo -e "INFO: Copy success"
	else
		echo -e "CRITICAL: Copy failed"
		((errCopy++))
	fi
}

checkErrorConfig () {
	if [ $? -eq 0 ]
	then
		echo -e "INFO: Configuration file success"
	else
		echo -e "CRITICAL: Configuration file failed"
		((errConfig++))
	fi
}

checkErrorBashrc () {
	if [ $? -eq 0 ]
	then
		echo -e "INFO: bashrc file success"
	else
		echo -e "CRITICAL: bashrc file failed"
		((errBashrc++))
	fi
}

# Copy binaries to /usr/local/*
echo -e "Deploy"
cp freefem-source-master/bin/* /usr/local/bin/
checkErrorCopy
cp -R freefem-source-master/ff-petsc /usr/local/
checkErrorCopy
cp -R freefem-source-master/lib/* /usr/local/lib/
checkErrorCopy
cp -R freefem-source-master/share/* /usr/local/share/
checkErrorCopy

# Create /etc/freefem++.pref
echo -e "Create configuration file"
echo -e "verbosity = 5" > /etc/freefem++.pref
checkErrorConfig
echo -e "#LOAD" >> /etc/freefem++.pref
echo -e "loadpath += \"/usr/local/lib/ff++/3.61/lib\"" >> /etc/freefem++.pref
checkErrorConfig
echo -e "#INCLUDE" >> /etc/freefem++.pref
echo "includepath += \"/usr/local/lib/ff++/3.61/idp\"" >> /etc/freefem++.pref
checkErrorConfig

echo -e "Summary"
echo -e "\tCopy error: $errCopy"
echo -e "\tFile error: $errConfig"
echo -e "\tBashrc error: $errBashrc"

# Append to bashrc
echo -e "" >> /etc/bash.bashrc
checkErrorBashrc
echo -e "LD_LIBRARY_PATH=\"$LD_LIBRARY_PATH:/usr/local/ff-petsc/real/lib:/usr/local/ff-petsc/complex/lib\"" >> /etc/bash.bashrc
checkErrorBashrc

if [ $errCopy -gt 0 ]
then
	echo -e "Some copy errors occurs"
fi
if [ $errConfig -gt 0 ]
then
	echo -e "Some configuration errors occurs"
fi
if [ $errBashrc -gt 0 ]
then
	echo -e "Some bashrc errors occurs"
fi

if [ $errCopy -eq 0 ] && [ $errConfig -eq 0 ] && [ $errBashrc -eq 0 ]
then
	echo -e "You can now use the FreeFem++ command directly in your terminal"
else
	echo -e "Try to use the FreeFem++ command in your terminal"
	echo -e "If that does not work, please post an issue in https://github.com/FreeFem/FreeFem-sources/issues"
fi



