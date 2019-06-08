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
# SUMMARY : Makefile for FreeFem++
# LICENSE : LGPLv3
# ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
# AUTHORS : ...
# E-MAIL  : ...

EXTRA_DIST=glrgraph.hpp  macglrgraf.cpp	\
macrgraf.cpp Pcrgraph.cpp  rgraph.hpp		\
sansrgraph.cpp xglrgraf.cpp Xrgraph.cpp DefColor.cpp   \
getprog-unix.hpp mode_open.hpp ffglut.hpp ff-win32.cpp

# FFCS: no compilation here, so the tags file must be built by hand
tags:TAGS
TAGS:
	etags *.?pp
clean-local::
	-rm TAGS
