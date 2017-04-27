#! \file   Makefile
#  \brief  source files in Driver/ directory and depending libraries
#  \author Atsushi Suzuki, Laboratoire Jacques-Louis Lions
#  \date   Jul. 12th 2015
#  \date   Feb. 29th 2016
#  \date   Nov. 30th 2016

# This file is part of Dissection
# 
# Dissection is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Linking Dissection statically or dynamically with other modules is making
# a combined work based on Disssection. Thus, the terms and conditions of 
# the GNU General Public License cover the whole combination.
#
# As a special exception, the copyright holders of Dissection give you 
# permission to combine Dissection program with free software programs or 
# libraries that are released under the GNU LGPL and with independent modules 
# that communicate with Dissection solely through the Dissection-fortran 
# interface. You may copy and distribute such a system following the terms of 
# the GNU GPL for Dissection and the licenses of the other code concerned, 
# provided that you include the source code of that other code when and as
# the GNU GPL requires distribution of source code and provided that you do 
# not modify the Dissection-fortran interface.
#
# Note that people who make modified versions of Dissection are not obligated 
# to grant this special exception for their modified versions; it is their
# choice whether to do so. The GNU General Public License gives permission to 
# release a modified version without this exception; this exception also makes
# it possible to release a modified version which carries forward this
# exception. If you modify the Dissection-fortran interface, this exception 
# does not apply to your modified version of Dissection, and you must remove 
# this exception when you distribute your modified version.
#
# This exception is an additional permission under section 7 of the GNU 
# General Public License, version 3 ("GPLv3")
#
# Dissection is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Dissection.  If not, see <http://www.gnu.org/licenses/>.

DIRECTORY  = Interface
SOURCESCPP = fortranAPI.cpp 

LIBDEP     = -L../../lib -lDriver -lUtilities $(SHARED_BLAS) $(SHARED_INTEL) 

include ../Defs.inc
include ./Make.inc
