# Makefile using Automake + Autoconf
# ----------------------------------
# $Id$
# change FH mars 2010 for sgi mpiu .... 

# FFCS - we need to call MPICXX directly because Windows options are too different to be detailed to the FF configuration script.

CXX=@MPICXX@

# To create statically linked executables (see configure.ac)
CXXLD=$(STATICTOOL) $(CXX)

# MPI may not be installed. Checked in configure.ac
bin_PROGRAMS=$(MPIPROG)
bin_SCRIPTS=$(MPISCRIPT)
EXTRA_PROGRAMS=FreeFem++-mpi
EXTRA_SCRIPTS=ff-mpirun
EXTRA_DIST=ff-mpirun.in

# FFCS: visualization stream redirection
FreeFem___mpi_SOURCES=../Graphics/sansrgraph.cpp ../lglib/mymain.cpp parallelempi.cpp ../lglib/lg.tab.cpp	\
	../fflib/ffapi.cpp
FreeFem___mpi_DEPENDENCIES=../fflib/libff.a
LDADD=../fflib/libff.a @UMFPACKLIBS@ @ARPACKLIBS@ @BLASLIBS@ @MPI_LIB@

AM_CXXFLAGS=-DPARALLELE
AM_CPPFLAGS=-I$(srcdir)/../fflib -I$(srcdir)/../Graphics -I$(srcdir)/../femlib @MPI_INCLUDE@
ff-mpirun:ff-mpirun.in Makefile
	../../config.status  --file=ff-mpirun:ff-mpirun.in
	chmod a+x  ff-mpirun

