# Makefile using Automake + Autoconf
# ----------------------------------
# $Id$
# change FH 
CXX=@MPICXX@

# To create statically linked executables (see configure.ac)
CXXLD=$(STATICTOOL) $(CXX)

# MPI may not be installed. Checked in configure.ac
bin_PROGRAMS=$(MPIPROG)
EXTRA_PROGRAMS=FreeFem++-mpi

FreeFem___mpi_SOURCES=../Graphics/sansrgraph.cpp ../lglib/mymain.cpp   parallelempi.cpp	\
../lglib/lg.tab.cpp
FreeFem___mpi_DEPENDENCIES=../fflib/libff.a
LDADD=../fflib/libff.a @UMFPACKLIB@ @ARPACKLIB@ @BLASLIB@

AM_CXXFLAGS=-DPARALLELE
INCLUDES=-I../fflib -I../Graphics -I../femlib
