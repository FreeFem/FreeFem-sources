# Makefile using Automake + Autoconf
# ----------------------------------
# $Id$

bin_PROGRAMS=$(AGLPROG)
EXTRA_PROGRAMS=FreeFem++-agl
EXTRA_DIST=FreeFem++.app.tgz Info-plist.am 
FreeFem___agl_SOURCES=../Graphics/glrgraph.hpp ../Graphics/macglrgraf.cpp
FreeFem___agl_DEPENDENCIES=../fflib/libff.a ../lglib/liblg.a
LDADD=../lglib/liblg.a ../fflib/libff.a @UMFPACKLIB@ @ARPACKLIB@ @BLASLIB@
INCLUDES=-I../fflib -I../Graphics -I../femlib

# Automake does not want to see "-framework" in LDADD
AM_LDFLAGS=-framework Carbon -framework OpenGL -framework AGL

install-exec-local:: INSTALL_APP_$(AGLPROG)

INSTALL_APP_FreeFem++-agl::
	echo Install /Applications/FreeFem++.app;
	rm -rf /Applications/FreeFem++.app;
	zcat FreeFem++.app.tgz|(cd /Applications;tar xf -) ;
	sed <Info-plist.am >/Applications/FreeFem++.app/Contents/Info.plist           \
		-e "s/@VVERSION@/$(VERSION)-$(FREEFEM_PACKAGING)/g" -e "s/@DATE@/`date`/g" ;
	cp $(AGLPROG) /Applications/FreeFem++.app/Contents/FreeFem++.app/FreeFem++
	mkdir -p /Applications/FreeFem++.app/Contents/FreeFem++.app/Contents
	cp /Applications/FreeFem++.app/Contents/Info.plist /Applications/FreeFem++.app/Contents/FreeFem++.app/Contents
	echo Install FreeFem++-CoCoa;
	cp $(top_srcdir)/src/FreeFem++-CoCoa $(bindir);
	echo "End of Installation";
INSTALL_APP_::

