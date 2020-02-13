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
# https://gforge.inria.fr/frs/download.php/file/34099/scotch_6.0.4.tar.gz
# http://gforge.inria.fr/frs/download.php/23391/scotch_5.1.7_esmumps.tar.gz
# Scotch information
DIRPKG=../pkg
SRCDIR=scotch_$(scotch_VERSION)
PACKAGE=$(DIRPKG)/scotch_$(scotch_VERSION_LOCAL).tar.gz
SERVER=https://gforge.inria.fr/frs/download.php/file/34099
#https://gforge.inria.fr/frs/download.php/file/34099/scotch_6.0.4.tar.gz
#https://gforge.inria.fr/frs/download.php/23391
INSTALL=../..
scotch_VERSION=6.0.4
scotch_VERSION_LOCAL=6.0.4
#  trick to no in MPI on NOT .. 
NOTESTMPI= test -z "$(CCP)"
scotch: FAIRE

$(SRCDIR)/AFAIRE:
	$(MAKE) install 
	$(NOTESTMPI) || $(MAKE) install-mpi

FAIRE: $(SRCDIR)/FAIT
	$(MAKE)  WHERE
	$(NOTESTMPI) || $(MAKE) WHERE-mpi
	touch FAIRE

Makefile.inc: ../../config.status Makefile-scotch.inc Makefile
	../../config.status  --file="Makefile.inc:Makefile-scotch.inc"

# DCOMMON_TIMING_OLD
# FFCS - 16/1/13 - what was "$(CCP)" used for?
$(SRCDIR)/FAIT: $(SRCDIR) Makefile.inc
	cp Makefile.inc $(SRCDIR)/src
	cd $(SRCDIR)/src && $(MAKE) scotch esmumps
	touch $(SRCDIR)/FAIT
$(SRCDIR)/FAIT-mpi: $(SRCDIR) Makefile.inc
	cp Makefile.inc $(SRCDIR)/src
	cd $(SRCDIR)/src && $(MAKE) -j1 ptscotch ptesmumps 
	touch $(SRCDIR)/FAIT-mpi

install:$(SRCDIR)/FAIT
	mkdir -p $(SRCDIR)/$(INSTALL)/include/scotch
	cp $(SRCDIR)/include/*.h $(SRCDIR)/$(INSTALL)/include/scotch/
	cp $(SRCDIR)/lib/*.a $(SRCDIR)/$(INSTALL)/lib/
install-mpi:$(SRCDIR)/FAIT-mpi
	mkdir -p $(SRCDIR)/$(INSTALL)/include/scotch
	cp $(SRCDIR)/include/*.h $(SRCDIR)/$(INSTALL)/include/scotch/
	cp $(SRCDIR)/lib/*.a $(SRCDIR)/$(INSTALL)/lib/

# ALH - 16/1/13 - simplifying make rules
WHERE:$(SRCDIR)/FAIT
	$(MAKE) install
	-rm $(SRCDIR)/$(INSTALL)/lib/WHERE.scotch
	echo scotch LD -L@DIR@/lib  -lesmumps -lscotch -lscotcherr >>$(SRCDIR)/$(INSTALL)/lib/WHERE.scotch
#	echo scotchmetis LD -L@DIR@/lib  -lscotchmetis >>$(SRCDIR)/$(INSTALL)/lib/WHERE.scotch
	echo scotch INCLUDE -I@DIR@/include/scotch >> $(SRCDIR)/$(INSTALL)/lib/WHERE.scotch
WHERE-mpi:$(SRCDIR)/FAIT-mpi
	$(MAKE) install-mpi
	-rm $(SRCDIR)/$(INSTALL)/lib/WHERE.scotch
	echo ptscotch LD -L@DIR@/lib  -lesmumps -lptscotch -lptscotcherr >>$(SRCDIR)/$(INSTALL)/lib/WHERE.scotch
#	echo ptscotchparmetis  LD -L@DIR@/lib  -lesmumps -lptscotchparmetis >>$(SRCDIR)/$(INSTALL)/lib/WHERE.scotch
	echo ptscotch INCLUDE -I@DIR@/include/scotch >> $(SRCDIR)/$(INSTALL)/lib/WHERE.scotch
	echo scotch LD -L@DIR@/lib  -lesmumps -lscotch -lscotcherr >>$(SRCDIR)/$(INSTALL)/lib/WHERE.scotch
#	echo scotchmetis LD -L@DIR@/lib  -lscotchmetis >>$(SRCDIR)/$(INSTALL)/lib/WHERE.scotch
	echo scotch INCLUDE -I@DIR@/include/scotch >> $(SRCDIR)/$(INSTALL)/lib/WHERE.scotch

$(SRCDIR): $(PACKAGE)
	tar xvzf $(PACKAGE)
	patch -p1 <scotch_6.0.4.patch
	#cd $(SRCDIR)/src/libscotch/; patch -p1 < ../../../scotch_5.1_esmumps.patch
	#cd $(SRCDIR)/src/esmumps/; patch -p1 < ../../../Makefile.patch
	touch $(SRCDIR)

$(PACKAGE):
	../getall -o Scotch -a

# FFCS - 28/3/13 - more cleaning
clean-local:
	-cd $(SRCDIR)/src && $(MAKE)  realclean  -C $(SRCDIR)/src
	-rm config.log
	-rm -rf ../include/scotch/
	-rm -rf ../lib/*scotch*
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
