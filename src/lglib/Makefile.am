# Makefile using Automake + Autoconf
# ----------------------------------

noinst_LIBRARIES=liblg.a
liblg_a_SOURCES=lg.tab.hpp lg.tab.cpp  mymain.cpp
LDADD = -lpthread
AM_CPPFLAGS=-I$(srcdir)/../fflib -I$(srcdir)/../bamglib -I$(srcdir)/../Graphics -I$(srcdir)/../femlib -I$(srcdir)/../ffgraphics/asio/asio/include -I$(srcdir)/../ffgraphics/json/single_include -I$(srcdir)/../ffgraphics/server

BUILT_SOURCES=lg.tab.hpp lg.tab.cpp

# <<lg.tab.?pp>> ALH - 10/9/13 - skip dependency on lg.ypp to make sure that lg.tab.?pp are not rebuilt on every new
# machine (the danger is to generate unnecessary version conflicts) on these automatically generated files

lg.tab.hpp lg.tab.cpp:
	-rm lg.tab.hpp  lg.tab.cpp lg.tab.cpp.h
	bison -dtv -p lg  lg.ypp  -o lg.tab.cpp
	-mv lg.tab.cpp.h  lg.tab.hpp # pour un  pb entre des versions de bison

libfreefem++.dylib: lg.tab.o ../fflib/libff.a
	g++ -bundle lg.tab.o  ../fflib/libff.a $(UMFPACKLIB) $(ARPACKLIB) $(BLASLIB) $(LIBS) -o libfreefem++.dylib
EXTRA_DIST=lg.ypp
