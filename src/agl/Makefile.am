# Makefile using Automake + Autoconf
# ----------------------------------
# $Id$

bin_PROGRAMS=$(AGLPROG)
EXTRA_PROGRAMS=FreeFem++-agl

FreeFem___agl_SOURCES=../Graphics/glrgraph.hpp ../Graphics/macglrgraf.cpp
FreeFem___agl_DEPENDENCIES=../fflib/libff.a ../lglib/liblg.a
LDADD=../lglib/liblg.a ../fflib/libff.a -ldl
INCLUDES=-I../fflib -I../Graphics -I../femlib

# Automake does not want to see "-framework" in LDADD
AM_LDFLAGS=-framework Carbon -framework OpenGL -framework AGL

install-exec-local::
	if [ \( -f $(COMPILE_DIR)/FreeFem++-agl \) -a \( -n "$(BUILD_COCOA)"  \)  ] ; then \
		echo " Install $(BIN_DIR)/../FreeFem++.app/Contents/MacOS/FreeFem++"; \
		cp $(COMPILE_DIR)/FreeFem++-agl $(BIN_DIR)/../FreeFem++.app/Contents/MacOS/FreeFem++ ;\
		echo " Install FreeFem++.app/Resources/MacOS/FreeFem++-agl"; \
		cp $(COMPILE_DIR)/FreeFem++-agl ./FreeFem++.app/Contents/Resources/FreeFem++ ;\
		echo Install /Applications/FreeFem++.app; \
		(tar cf - FreeFem++.app| (cd /Applications;tar xf -)) ;\
		echo Install FreeFem++-CoCoa; \
		cp src/FreeFem++-CoCoa  $(BIN_DIR);\
	fi;
	@echo "End of Installation"
