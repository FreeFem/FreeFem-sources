# Makefile using Automake + Autoconf
# ----------------------------------
# $Id$

noinst_LIBRARIES=liblg.a
liblg_a_SOURCES=lg.tab.hpp lg.tab.cpp
INCLUDES=-I../fflib -I../bamglib -I../Graphics -I../femlib

BUILT_SOURCES=lg.tab.hpp lg.tab.cpp

lg.tab.hpp  lg.tab.cpp: lg.y
	-rm lg.tab.hpp  lg.tab.cpp lg.tab.cpp.h
	bison -dtv -p lg  $< -o lg.tab.cpp
	-mv lg.tab.cpp.h  lg.tab.hpp # pour un  pb entre des versions de bison 

EXTRA_DIST=lg.y
