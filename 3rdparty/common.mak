# Common make rules for all downloaded packages (request from FH)
# ======================================================================
# Written by Antoine Le Hyaric
# http://www.ljll.math.upmc.fr/lehyaric
# Laboratoire Jacques-Louis Lions
# Université Pierre et Marie Curie-Paris6, UMR 7598, Paris, F-75005 France
# ======================================================================
# This file is part of Freefem++
# 
# Freefem++ is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2.1 of
# the License, or (at your option) any later version.
# 
# Freefem++ is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public
# License along with Freefem++; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
# ======================================================================
# headeralh brief="Common make rules for all downloaded packages (request from FH)" default=0 freefem make start=06/11/2013 upmc written

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
