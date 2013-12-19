# ======================================================================
# Laboratoire Jacques-Louis Lions
# Université Pierre et Marie Curie-Paris6, UMR 7598, Paris, F-75005 France
# http://www.ljll.math.upmc.fr/lehyaric
# ======================================================================
# This file is part of Freefem++
# 
# Freefem++ is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2.1 of
# the License, or (at your option) any later version.
# 
# Freefem++ is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public
# License along with Freefem++; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
# ======================================================================
# headeralh default=0 freefem make multipleauthors start=19/03/10 upmc

# Downloading and compiling extra libraries
# -----------------------------------------

include Makefile.inc

all-local: scotch

# Downloading and compiling scotch
# -------------------------------
# http://gforge.inria.fr/frs/download.php/23391/scotch_5.1.7_esmumps.tar.gz
# Scotch information
DIRPKG=../pkg
SRCDIR=scotch_$(scotch_VERSION)
PACKAGE=$(DIRPKG)/scotch_$(scotch_VERSION_LOCAL).tar.gz
SERVER=https://gforge.inria.fr/frs/download.php/23391
INSTALL=../..
scotch_VERSION=5.1_esmumps
scotch_VERSION_LOCAL=5.1.6_esmumps
#  trick to no in MPI on NOT .. 
TESTMPI= test -n '$(CCP)'
scotch: FAIRE

$(SRCDIR)/AFAIRE:
	$(MAKE) install 

FAIRE: $(SRCDIR)/FAIT
	$(MAKE)  WHERE
	touch FAIRE
Makefile.inc: ../../config.status Makefile-scotch.inc Makefile
	../../config.status  --file="Makefile.inc:Makefile-scotch.inc"
	case `uname` in *CYGWIN*)  cp Makefile.inc Makefile.inn; sed "s/COPTIONS =/COPTIONS = -D__VC__/" <Makefile.inn > Makefile.inc; rm Makefile.inn ;; esac 
	case `uname` in *CYGWIN*|*MINGW*)  cp Makefile.inc Makefile.inn; sed -e "s/-DCOMMON_TIMING_OLD//" -e "s/-DCOMMON_RANDOM_FIXED_SEED/-DCOMMON_RANDOM_RAND  -DCOMMON_RANDOM_FIXED_SEED -D'pipe(pfds)=_pipe(pfds,1024,0x8000)'/" <Makefile.inn > Makefile.inc; rm Makefile.inn ;; esac 

# FFCS - 16/1/13 - what was "$(CCP)" used for?
$(SRCDIR)/FAIT: $(SRCDIR) Makefile.inc
	cp Makefile.inc $(SRCDIR)/src
	cd $(SRCDIR)/src && $(MAKE) scotch
#
#	FFCS - 16/1/13 - this crashes when run in parallel
#       FH: compile pt only if mpi is here 
	if $(TESTMPI) ; then  cd $(SRCDIR)/src && $(MAKE) -j1 ptscotch ; fi
	touch $(SRCDIR)/FAIT

install:$(SRCDIR)/FAIT
	mkdir -p $(SRCDIR)/$(INSTALL)/include/scotch
	cp $(SRCDIR)/include/*.h $(SRCDIR)/$(INSTALL)/include/scotch/
	cp $(SRCDIR)/lib/*.a $(SRCDIR)/$(INSTALL)/lib/

# ALH - 16/1/13 - simplifying make rules
WHERE:$(SRCDIR)/FAIT
	$(MAKE) install
	-rm $(SRCDIR)/$(INSTALL)/lib/WHERE.scotch
	-$(TESTMPI) &&echo ptscotch LD -L@DIR@/lib -lptesmumps -lptscotch -lptscotcherr >>$(SRCDIR)/$(INSTALL)/lib/WHERE.scotch
	-$(TESTMPI) &&echo ptscotch INCLUDE -I@DIR@/include/scotch >> $(SRCDIR)/$(INSTALL)/lib/WHERE.scotch
	echo scotch LD -L@DIR@/lib -lscotch -lscotcherr >>$(SRCDIR)/$(INSTALL)/lib/WHERE.scotch
	echo scotch INCLUDE -I@DIR@/include/scotch >> $(SRCDIR)/$(INSTALL)/lib/WHERE.scotch

$(SRCDIR): $(PACKAGE)
	tar xvzf $(PACKAGE)
	cd $(SRCDIR)/src/libscotch/; patch -p1 < ../../../scotch_5.1_esmumps.patch
	cd $(SRCDIR)/src/esmumps/; patch -p1 < ../../../Makefile.patch
	touch $(SRCDIR)

# FFCS - 28/3/13 - more cleaning
clean-local:
	-cd $(SRCDIR)/src && $(MAKE)  realclean  -C $(SRCDIR)/src
	-rm config.log
# FH 	-rm $(PACKAGE) 

clean: clean-local
	-rm Makefile.inc
	-rm -rf $(SRCDIR)
	-rm FAIRE 
.PHONY:$(SRCDIR)/$(INSTALL)

# Local Variables:
# mode:makefile
# ispell-local-dictionary:"british"
# coding:utf-8
# End:
