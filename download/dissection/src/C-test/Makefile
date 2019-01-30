#! \file   Makefile
#  \brief  to create MM-Dissection
#  \author Atsushi Suzuki, Laboratoire Jacques-Louis Lions
#  \date   Jun. 20th 2014
#  \date   Jul. 12th 2015
#  \date   Feb. 29th 2016

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
#

#
include ../Defs.inc
DISSECT_LIB_DIR = $(PROJ_LIBDIR)
#
# modifications to predefined variables in Defs.inc
#
ifeq ($(UNAME), Darwin)
   DISSECT_LIB_PATH = -Xlinker -rpath -Xlinker $(DISSECT_LIB_DIR)
   SHARED_BLAS_PATH = -Xlinker -rpath -Xlinker $(LIB_DIR_MKL)	
   SHARED_INTEL_PATH = -Xlinker -rpath -Xlinker $(LIB_DIR_INTEL)	
   LIB_SYSTEM =
endif 
ifeq ($(UNAME), Linux)
   DISSECT_LIB_PATH  = -Xlinker -rpath=$(DISSECT_LIB_DIR)
   SHARED_BLAS_PATH = -Xlinker -rpath=$(LIB_DIR_MKL)
   SHARED_INTEL_PATH = -Xlinker -rpath=$(LIB_DIR_INTEL)	
   LIB_SYSTEM = -lrt
endif
#
ifeq ($(COMPILERBASE), INTEL)
   SHARED_BLAS_PATH = -Xlinker -rpath -Xlinker $(LIB_DIR_MKL) -L$(LIB_DIR_MKL)
#   SHARED_BLAS = -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
   SHARED_INTEL_PATH = -Xlinker -rpath -Xlinker $(LIB_DIR_INTEL) -L$(LIB_DIR_INTEL)
#   SHARED_INTEL = -lifport -lifcoremt -lpthread

   BLASLIB    = $(SHARED_BLAS_PATH) $(SHARED_BLAS) \
                $(SHARED_INTEL_PATH) $(SHARED_INTEL) 
else
ifeq ($(COMPILER), GNU)
   ifeq ($(BLAS), ATLAS)
   LIB_DIR_ATLAS = /usr/local/atlas/lib
   SHARED_BLAS_PATH = -Xlinker -rpath -Xlinker $(LIB_DIR_ATLAS) -L$(LIB_DIR_ATLAS)
   SHARED_BLAS = -lcblas -lsatlas 
   endif
   BLASLIB    = $(SHARED_BLAS_PATH) $(SHARED_BLAS) $(SHARED_GNU) #-lrt
else
ifeq ($(COMPILER), SUN)
else
   BLASLIB    = $(SHARED_BLAS) 
endif
endif
endif
#
# DEBUGFLAG and FFFLAGS_OMP are defined in Defs.inc
FFLAGS     = $(DEBUGFLAG) 
#
LDOPT      = $(DEBUGFLAG)
#
DISSECTLIB= $(DISSECT_LIB_PATH) \
	    -L$(DISSECT_LIB_DIR) -lDissection 

all: MM-Dissection

.cpp.o:
	@echo Compiling $<
	$(CXX) $(CCFLAGS) -I../ -c $< -o $@

# MM-Dissection-Cinterface: MM-Dissection-Cinterface.o
# 	@echo Create executable $(DIRECTORY)
# 	@echo -n "**** setting of variables : "
# 	@echo -n "COMPILER=$(COMPILER) "
# 	@echo -n "BLAS=$(BLAS) "
# 	@echo -n "PARALLEL=$(PARALLEL) "
# 	@echo    "****"
# 	$(LD) -o $@ $(LDOPT) MM-Dissection-Cinterface
#.o $(DISSECTLIB) $(BLASLIB) $(LIBFORTRAN) $(LIB_SYSTEM) $(SHARED_COMPILER)

MM-Dissection: MM-Dissection.o
	@echo Create executable $(DIRECTORY)
	@echo -n "**** setting of variables : "
	@echo -n "COMPILER=$(COMPILER) "
	@echo -n "BLAS=$(BLAS) "
	@echo -n "PARALLEL=$(PARALLEL) "
	@echo    "****"
	$(LD) -o $@ $(LDOPT) MM-Dissection.o $(DISSECTLIB) $(BLASLIB) $(LIBFORTRAN) $(LIB_SYSTEM) $(SHARED_COMPILER)

MM-Dissection-quad: MM-Dissection-quad.o
	@echo Create executable $(DIRECTORY)
	@echo -n "**** setting of variables : "
	@echo -n "COMPILER=$(COMPILER) "
	@echo -n "BLAS=$(BLAS) "
	@echo -n "PARALLEL=$(PARALLEL) "
	@echo    "****"
	$(LD) -o $@ $(LDOPT) MM-Dissection-quad.o $(DISSECTLIB) $(BLASLIB) $(LIBFORTRAN) $(LIB_SYSTEM) $(SHARED_COMPILER)

MM-Dissection-double-quad: MM-Dissection-double-quad.o
	@echo Create executable $(DIRECTORY)
	@echo -n "**** setting of variables : "
	@echo -n "COMPILER=$(COMPILER) "
	@echo -n "BLAS=$(BLAS) "
	@echo -n "PARALLEL=$(PARALLEL) "
	@echo    "****"
	$(LD) -o $@ $(LDOPT) MM-Dissection-double-quad.o $(DISSECTLIB) $(BLASLIB) $(LIBFORTRAN) $(LIB_SYSTEM) $(SHARED_COMPILER)

MM-Dissection-mRHS: MM-Dissection-mRHS.o
	@echo Create executable $(DIRECTORY)
	@echo -n "**** setting of variables : "
	@echo -n "COMPILER=$(COMPILER) "
	@echo -n "BLAS=$(BLAS) "
	@echo -n "PARALLEL=$(PARALLEL) "
	@echo    "****"
	$(LD) -o $@ $(LDOPT) MM-Dissection-mRHS.o $(DISSECTLIB) $(BLASLIB) $(LIBFORTRAN) -lrt

clean:
	@echo Cleaning the trash...
	rm -fr *~ *.o *.so *.mod core *.d *.dll $(DIRECTORY) \
        MM-Dissection MM-Dissection-quad MM-Dissection-double-quad
