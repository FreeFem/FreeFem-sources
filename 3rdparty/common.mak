############################################################################
# This file is part of FreeFEM.                                            #
#                                                                          #
# FreeFEM is free software: you can redistribute it and/or modify          #
# it under the terms of the GNU Lesser General Public License as           #
# published by the Free Software Foundation, either version 3 of           #
# the License, or (at your option) any later version.                      #
#                                                                          #
# FreeFEM is distributed in the hope that it will be useful,               #
# but WITHOUT ANY WARRANTY; without even the implied warranty of           #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
# GNU Lesser General Public License for more details.                      #
#                                                                          #
# You should have received a copy of the GNU Lesser General Public License #
# along with FreeFEM. If not, see <http://www.gnu.org/licenses/>.          #
############################################################################
# SUMMARY : Common make rules for all downloaded packages (request from FH)
# LICENSE : LGPLv3
# ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
# AUTHORS : ...
# E-MAIL  : ...

# Common goals for all packages:
# download compile install reinstall clean veryclean

# <<download>>

# PKGCOMMON_PACKTITLE corresponds to package names in [[file:getall]]
download::
	../getall -o $(PKGCOMMON_PACKTITLE) -a
$(PKGCOMMON_PACKAGES):download

## regle qui force le telecharmeent a l'install je vide F. H 
compilepkg::

# <<install>>

install::compilepkg

# <<reinstall>>

reinstall::compilepkg

clean-local::

veryclean::clean
	-rm $(PKGCOMMON_PACKAGES)

# Local Variables:
# mode:makefile
# ispell-local-dictionary:"british"
# coding:utf-8
# End:
