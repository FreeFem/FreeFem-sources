# Makefile using Automake + Autoconf
# ----------------------------------
# $Id$

noinst_LIBRARIES=liblg.a
liblg_a_SOURCES=lg.tab.hpp lg.tab.cpp  mymain.cpp
INCLUDES=-I$(srcdir)/../fflib -I$(srcdir)/../bamglib -I$(srcdir)/../Graphics -I$(srcdir)/../femlib

BUILT_SOURCES=lg.tab.hpp lg.tab.cpp

lg.tab.hpp  lg.tab.cpp: lg.ypp
	-rm lg.tab.hpp  lg.tab.cpp lg.tab.cpp.h
	bison -dtv -p lg  $< -o lg.tab.cpp
	-mv lg.tab.cpp.h  lg.tab.hpp # pour un  pb entre des versions de bison 
libfreefem++.dylib: lg.tab.o ../fflib/libff.a
	g++ -bundle lg.tab.o  ../fflib/libff.a $(UMFPACKLIB) $(ARPACKLIB) $(BLASLIB) $(LIBS) -o libfreefem++.dylib
EXTRA_DIST=lg.ypp
