# Makefile using Automake + Autoconf
# ----------------------------------
# $Id$

bin_PROGRAMS=$(AGLPROG)
EXTRA_PROGRAMS=FreeFem++-agl

FreeFem___agl_SOURCES=../Graphics/glrgraph.hpp ../Graphics/macglrgraf.cpp
FreeFem___agl_DEPENDENCIES=../fflib/libff.a ../lglib/liblg.a
LDADD=../lglib/liblg.a ../fflib/libff.a @UMFPACKLIB@ @ARPACKLIB@ @BLASLIB@
INCLUDES=-I../fflib -I../Graphics -I../femlib

# Automake does not want to see "-framework" in LDADD
AM_LDFLAGS=-framework Carbon -framework OpenGL -framework AGL

install-exec-local:: FreeFem++-agl
		echo " Install $(top_srcdir)/FreeFem++.app/Contents/MacOS/FreeFem++"
		cp "FreeFem++-agl" $(top_srcdir)/FreeFem++.app/Contents/MacOS/FreeFem++ 
		echo " Install FreeFem++.app/Resources/MacOS/FreeFem++-agl"
		cp "FreeFem++-agl" $(top_srcdir)/FreeFem++.app/Contents/Resources/FreeFem++ 
		echo Install /Applications/FreeFem++.app; 
		sed <$(top_srcdir)/Info-plist.am >$(top_srcdir)/FreeFem++.app/Contents/Info.plist \
	            -e "s/@Version@/$(VERSION)/g"   -e "s/@DATE@/`date`/g" ; 
		(cd $(top_srcdir); tar cf - FreeFem++.app| (cd /Applications;tar xf -)) ; 
		echo Install FreeFem++-CoCoa; 
		cp $(top_srcdir)/src/FreeFem++-CoCoa  $(bindir);
	@echo "End of Installation"
