# Makefile using Automake + Autoconf
# ----------------------------------
# $Id$

CXX=mpiCC

bin_PROGRAMS=FreeFem++-mpi
FreeFem___mpi_SOURCES=../Graphics/sansrgraph.cpp parallelempi.cpp	\
../lglib/lg.tab.cpp
FreeFem___mpi_DEPENDENCIES=../fflib/libff.a
FreeFem___mpi_LDADD=../fflib/libff.a

CXXFLAGS=-DPARALLELE
INCLUDES=-I../fflib -I../Graphics -I../femlib
LIBS=-ldl
