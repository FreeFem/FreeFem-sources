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
# SUMMARY : Autotools Makefile for parallel version of FreeFEM
# LICENSE : LGPLv3
# ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
# AUTHORS : ...
# E-MAIL  : ...

# change SG june 2019: remove file copy, use PHONY to force rebuild target with different DEFINE
# change FH mars 2010 for sgi mpiu ....
# FFCS - we need to call MPICXX directly because Windows options are too different to be detailed to the FF configuration script.

CXX = @MPICXX@

# To create statically linked executables (see configure.ac)
CXXLD = $(STATICTOOL) $(CXX)

# MPI may not be installed. Checked in configure.ac
bin_PROGRAMS = $(MPIPROG)
bin_SCRIPTS = $(MPISCRIPT)
EXTRA_PROGRAMS = FreeFem++-mpi
EXTRA_SCRIPTS = ff-mpirun
EXTRA_DIST = ff-mpirun.in

# FFCS: visualization stream redirection
FreeFem___mpi_SOURCES = /ffapi.cpp \
	../Graphics/sansrgraph.cpp \
	../lglib/mymain.cpp \
	../lglib/lg.tab.cpp \
	parallelempi.cpp
FreeFem___mpi_DEPENDENCIES = ../fflib/libff.a
LDADD = ../fflib/libff.a @UMFPACKLIBS@ @ARPACKLIBS@ @BLASLIBS@ @MPI_LIB@

AM_CXXFLAGS = -DPARALLELE
AM_CPPFLAGS = -I$(srcdir)/../fflib -I$(srcdir)/../Graphics -I$(srcdir)/../femlib @MPI_INCLUDE@ $(HDF5_CPPFLAGS) $(UMFPACK_CPPFLAGS) -I$(srcdir)/../ffgraphics/asio/asio/include -I$(srcdir)/../ffgraphics/json/single_include -I$(srcdir)/../ffgraphics/server
ff-mpirun: ff-mpirun.in Makefile
	../../config.status  --file=ff-mpirun:ff-mpirun.in
	chmod a+x ff-mpirun

# # Force rebuild with current flags: PARALLELE
# .PHONY: ../fflib/ffapi.cpp ## That does not work for Windows executable
## TODO: same compilation process for *Nix and Windows ?
# copy this file because we use different FLAGS (-DPARALLELE)
ffapi.cpp: ../fflib/ffapi.cpp
	cp $^ .
