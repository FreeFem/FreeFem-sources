# Makefile using Automake + Autoconf
# ----------------------------------
# $Id$

bin_PROGRAMS=$(PCPROG)
EXTRA_PROGRAMS=FreeFem++-pc
FreeFem___pc_SOURCES=../Graphics/Pcrgraph.cpp
FreeFem___pc_DEPENDENCIES=../fflib/libff.a ../lglib/liblg.a
AM_LDFLAGS= -mwindows
LDADD=../lglib/liblg.a ../fflib/libff.a
INCLUDES=-I../fflib
