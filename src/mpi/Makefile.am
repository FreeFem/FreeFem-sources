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
FreeFem___mpi_SOURCES=sansrgraph.cpp mymain.cpp parallelempi.cpp lg.tab.cpp	ffapi.cpp
FreeFem___mpi_DEPENDENCIES=../fflib/libff.a
LDADD=../fflib/libff.a @UMFPACKLIBS@ @ARPACKLIBS@ @BLASLIBS@ @MPI_LIB@

AM_CXXFLAGS=-DPARALLELE
AM_CPPFLAGS=-I$(srcdir)/../fflib -I$(srcdir)/../Graphics -I$(srcdir)/../femlib @MPI_INCLUDE@ $(HDF5_CPPFLAGS) $(UMFPACK_CPPFLAGS)
ff-mpirun:ff-mpirun.in Makefile
	../../config.status  --file=ff-mpirun:ff-mpirun.in
	chmod a+x  ff-mpirun


## copie file because no compile with same FLAGS (with -DPARALLELE)
sansrgraph.cpp:../Graphics/sansrgraph.cpp
	cp $^ .
mymain.cpp:../lglib/mymain.cpp
	cp $^ .
lg.tab.cpp:../lglib/lg.tab.cpp
	cp $^ .
ffapi.cpp:../fflib/ffapi.cpp
	cp $^ .