# Makefile using Automake + Autoconf
# ----------------------------------
# $Id$

CXX=mpiCC

# MPI may not be installed. Checked in configure.ac
bin_PROGRAMS=$(MPIPROG)
EXTRA_PROGRAMS=FreeFem++-mpi

FreeFem___mpi_SOURCES=../Graphics/sansrgraph.cpp parallelempi.cpp	\
../lglib/lg.tab.cpp
FreeFem___mpi_DEPENDENCIES=../fflib/libff.a
LDADD=../fflib/libff.a $(LIBDL)

AM_CXXFLAGS=-DPARALLELE
INCLUDES=-I../fflib -I../Graphics -I../femlib
