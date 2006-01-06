# Makefile using Automake + Autoconf
# ----------------------------------
# $Id$
CXXLD=$(STATICTOOL) $(CXX)

bin_PROGRAMS=bamg  cvmsh2
bamg_SOURCES=bamg.cpp 
cvmsh2_SOURCES=cvmsh2.cpp
bamg_DEPENDENCIES=../fflib/libff.a 
cvmsh2_DEPENDENCIES=../fflib/libff.a 
LDADD= ../fflib/libff.a
INCLUDES=-I../fflib -I../bamglib -I../Graphics


