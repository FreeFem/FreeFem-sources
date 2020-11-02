# and compiling extra libraries
# -----------------------------------------

include cxxflags
all-local: tetgen

# Downloading and compiling Tetgen
# ------------------------------

# Tetgen information
DIRPKG=../pkg
tetgen_VERSION=1.5.1-beta1
SRCDIR=tetgen$(tetgen_VERSION)
PACKAGE=$(DIRPKG)/tetgen$(tetgen_VERSION).tar.gz

# FFCS - 6/11/12 - curl is not able to follow redirections from http://tetgen.org/files
#SERVER=http://wias-berlin.de/software/tetgen/files

INSTALL=../..

tetgen: FAIRE 


FAIRE:install.done 
	touch FAIRE


# ALH - FFCS - 18/12/8 - need '-fPIC' on Linux64 because the .a will be used in a .so (by examples++-load/tetgen.cpp).

# FFCS - 30/11/10 - need ranlib on Win64. The PATH is setup so that mingw/ranlib is called

fait-tags.done:tar-tags.done cxxflags 
	cd $(SRCDIR);$(CXX) $(CXXFLAGS) -O3 -fPIC -DSELF_CHECK  -DNDEBUG -DTETLIBRARY -c tetgen.cxx
	cd $(SRCDIR);$(CXX) $(CXXFLAGS) -O0 -fPIC -DSELF_CHECK  -DNDEBUG -DTETLIBRARY -c predicates.cxx
	touch fait-tags.done
WHERE:WHERE.done
WHERE.done: fait-tags.done
	echo tetgen LD -L@DIR@/lib -ltet  >$(SRCDIR)/$(INSTALL)/lib/WHERE.tetgen
	echo tetgen  INCLUDE -I@DIR@/include >> $(SRCDIR)/$(INSTALL)/lib/WHERE.tetgen
	touch $@
clean-local::
	-rm WHERE.done

# FFCS - avoid remaking install every time
install.done:fait-tags.done
	cd $(SRCDIR);$(AR) $(ARFLAGS)  $(INSTALL)/lib/libtet.a tetgen.o predicates.o
	ranlib $(SRCDIR)/$(INSTALL)/lib/libtet.a
	cp $(SRCDIR)/tetgen.h $(SRCDIR)/$(INSTALL)/include 
	touch $@
clean-local::
	-rm install.done


tar-tags.done: $(PACKAGE)
	gunzip -c  $(PACKAGE)|tar zxf - 
#       FFCS: needs to patch tetgen because mingw64 has 4-byte longs
	touch tar-tags.done
#	cd tetgen1.4.2;patch -p1 <../tetgen1.4.2.patch 

$(PACKAGE):
	../getall -o TetGen -a

clean-local::
	-rm -rf tetgen1.?.*  *.done
clean:
	-rm $(SRCDIR)
	-rm *.done
	-rm -rf tetgen1.?.* cxxflags
cxxflags: ../Makefile  Makefile ../getall 
	-mv -f cxxflags cxxflags.old
	grep 'CXX *=' ../Makefile >cxxflags
	grep 'CC *=' ../Makefile >>cxxflags
#	FFCS - 10/5/12 - bug under Windows if -O3 is specified tetgen never returns. It could also be the case under
#	Ubuntu. All optimisation options are removed for safety.
	grep 'CXXFLAGS *=' ../Makefile | sed 's/ -O[0-9]* / /g'  >>cxxflags
	grep 'WGET *=' ../Makefile >>cxxflags
	grep 'AR *=' ../Makefile >>cxxflags
	grep 'ARFLAGS *=' ../Makefile >>cxxflags
	../../bin/ff-md5 $(PACKAGE) >> cxxflags	
	if diff -q cxxflags cxxflags.old ; then echo No Modif skip compile of tetgen ; mv -f cxxflags.old cxxflags; else touch -c tar-tags.done; fi
	
.PHONY:$(SRCDIR)/$(INSTALL)